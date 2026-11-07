#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* One reported motif: bases covered under the screen rule and as tracts. */
typedef struct {
    int32_t id, screen, tract;
} Hit;
/* Per-motif counters. An epoch stamp resets a motif on first use in each scan. */
typedef struct {
    int32_t n_motifs;      /* number of canonical motifs */
    int32_t epoch;         /* scan counter, compared with stamp[] */
    int32_t n_active;      /* motifs touched in the current scan */
    int32_t *map;          /* motif id of every word of length 1-8 */
    int32_t *offset;       /* offset[k]: index in map of the first length-k word */
    int32_t *stamp;        /* epoch of each motif's last use */
    int32_t *active;       /* ids of the motifs touched in the current scan */
    int32_t *screen;       /* screen-rule bases covered, per motif */
    int32_t *tract;        /* tract bases covered, per motif */
    int32_t *screen_end;   /* end of the last screen interval added, per motif */
    int32_t *tract_end;    /* end of the last tract interval added, per motif */
    uint64_t *screen_bits; /* which of the last 64 positions are screen-covered, per motif */
    uint64_t *tract_bits;  /* which of the last 64 positions are tract-covered, per motif */
} Scanner;

/* 2-bit code of a base, or -1 for anything other than ACGT. */
static int base(unsigned char c) {
    switch (c) {
    case 'A':
        return 0;
    case 'C':
        return 1;
    case 'G':
        return 2;
    case 'T':
        return 3;
    default:
        return -1;
    }
}

/* Number of bases that differ between two 2-bit packed words. */
static int distance(uint64_t a, uint64_t b) {
    uint64_t x = a ^ b;
    return __builtin_popcountll((x | (x >> 1)) & UINT64_C(0x5555555555555555));
}

/* map: motif id of every word (-1 if not primitive). offset[k]: first length-k word. */
void *create_scanner(int32_t n_motifs, const int32_t *map, const int32_t *offset, int32_t n_words) {
    Scanner *p = calloc(1, sizeof(*p));
    if (!p)
        return NULL;
    p->n_motifs = n_motifs;
    p->map = malloc(n_words * sizeof(int32_t));
    p->offset = malloc(9 * sizeof(int32_t));
    p->stamp = calloc(n_motifs, sizeof(int32_t));
    p->active = malloc(n_motifs * sizeof(int32_t));
    p->screen = malloc(n_motifs * sizeof(int32_t));
    p->tract = malloc(n_motifs * sizeof(int32_t));
    p->screen_end = malloc(n_motifs * sizeof(int32_t));
    p->tract_end = malloc(n_motifs * sizeof(int32_t));
    p->screen_bits = malloc(n_motifs * sizeof(uint64_t));
    p->tract_bits = malloc(n_motifs * sizeof(uint64_t));
    if (!p->map || !p->offset || !p->stamp || !p->active || !p->screen || !p->tract ||
        !p->screen_end || !p->tract_end || !p->screen_bits || !p->tract_bits)
        abort();
    memcpy(p->map, map, n_words * sizeof(int32_t));
    memcpy(p->offset, offset, 9 * sizeof(int32_t));
    return p;
}

void free_scanner(void *v) {
    Scanner *p = v;
    if (!p)
        return;
    free(p->map);
    free(p->offset);
    free(p->stamp);
    free(p->active);
    free(p->screen);
    free(p->tract);
    free(p->screen_end);
    free(p->tract_end);
    free(p->screen_bits);
    free(p->tract_bits);
    free(p);
}

/* Add [s, e) to a motif's screen or tract count. bits remembers the last 64 positions,
   so overlapping bases are counted once. */
static void add(Scanner *p, int id, int s, int e, int tract) {
    if (id < 0)
        return;
    if (p->stamp[id] != p->epoch) {
        p->stamp[id] = p->epoch;
        p->active[p->n_active++] = id;
        p->screen[id] = p->tract[id] = 0;
        p->screen_end[id] = p->tract_end[id] = -1;
        p->screen_bits[id] = p->tract_bits[id] = 0;
    }
    int *last = tract ? &p->tract_end[id] : &p->screen_end[id];
    int *count = tract ? &p->tract[id] : &p->screen[id];
    uint64_t *bits = tract ? &p->tract_bits[id] : &p->screen_bits[id];
    int delta = e - *last;
    *bits = delta >= 64 ? 0 : *bits << delta;
    uint64_t mask = (UINT64_C(1) << (e - s)) - 1;
    *count += __builtin_popcountll(mask & ~*bits);
    *bits |= mask;
    *last = e;
}

/* Screen and tract coverage of every motif in s. Reports motifs above minimum_percent. */
int scan(void *v, const char *s, int n, Hit *hits, int minimum_percent) {
    Scanner *p = v;
    if (++p->epoch == INT32_MAX) {
        memset(p->stamp, 0, p->n_motifs * sizeof(int32_t));
        p->epoch = 1;
    }
    p->n_active = 0;
    int run = 0, previous = -1;
    /* period 1: runs of 12 or more identical bases are tracts */
    for (int i = 0; i < n; i++) {
        int b = base(s[i]);
        run = (b >= 0 && b == previous) ? run + 1 : (b >= 0);
        previous = b;
        if (run >= 12) {
            int id = p->map[p->offset[1] + b];
            add(p, id, i - 11, i + 1, 1);
        }
    }
    /* periods 2-8: a sliding window of `copies` units, at least 3 units and 12 bp */
    for (int k = 2; k <= 8; k++) {
        int copies = (12 + k - 1) / k;
        if (copies < 3)
            copies = 3;
        int width = copies * k, valid = 0, off = p->offset[k];
        uint64_t w = 0, mask = (UINT64_C(1) << (2 * k)) - 1;
        for (int i = 0; i < n; i++) {
            int b = base(s[i]);
            if (b < 0) {
                valid = 0;
                w = 0;
                continue;
            }
            w = (w << 2) | (unsigned)b;
            if (valid < width)
                valid++;
            if (valid < k)
                continue;
            /* screen rule: exact units, and adjacent unit pairs within one substitution */
            uint64_t a = w & mask;
            int ia = p->map[off + a];
            add(p, ia, i - k + 1, i + 1, 0);
            if (valid < 2 * k)
                continue;
            uint64_t z = (w >> (2 * k)) & mask;
            int iz = p->map[off + z];
            int dz = distance(a, z);
            if (dz > 1)
                continue;
            add(p, ia, i - 2 * k + 1, i + 1, 0);
            add(p, iz, i - 2 * k + 1, i + 1, 0);
            /* tract: `copies` units within one substitution, no two mismatched in a row */
            if (valid < width)
                continue;
            uint64_t patterns[2] = {a, z};
            int ids[2] = {ia, iz};
            for (int t = 0; t < (a == z ? 1 : 2); t++) {
                if (ids[t] < 0)
                    continue;
                int last = 0, ok = 1;
                for (int j = 0; j < copies; j++) {
                    int d = distance((w >> (2 * k * j)) & mask, patterns[t]);
                    if (d > 1 || (j && d + last > 1)) {
                        ok = 0;
                        break;
                    }
                    last = d;
                }
                if (ok)
                    add(p, ids[t], i - width + 1, i + 1, 1);
            }
        }
    }

    /* report motifs whose screen or tract coverage exceeds the threshold */
    int nh = 0;
    for (int i = 0; i < p->n_active; i++) {
        int id = p->active[i];
        if ((int64_t)p->screen[id] * 100 > (int64_t)minimum_percent * n ||
            (int64_t)p->tract[id] * 100 > (int64_t)minimum_percent * n)
            hits[nh++] = (Hit){id, p->screen[id], p->tract[id]};
    }
    return nh;
}

/* One motif (np rotations of k bp) in s: out[0] screen-rule bp, out[1] exact-unit bp
   with gaps of up to 3 bp closed. np <= 16: k rotations on each strand. */
void target_coverage(const char *s, int n, const char *patterns, int k, int np, int *out) {
    uint64_t exact[16];
    for (int j = 0; j < np; j++) {
        exact[j] = 0;
        for (int t = 0; t < k; t++)
            exact[j] = (exact[j] << 2) | base(patterns[j * k + t]);
    }
    int valid = 0, end = -1, closed = 0, total = 0, last = -1;
    uint64_t covered = 0, w = 0, mask = (UINT64_C(1) << (2 * k)) - 1;
    /* period 1: runs of 12 or more of a listed base count for both outputs */
    if (k == 1) {
        int run = 0, prev = -1;
        for (int i = 0; i < n; i++) {
            int b = base(s[i]);
            run = (b >= 0 && b == prev) ? run + 1 : (b >= 0);
            prev = b;
            int match = 0;
            for (int t = 0; t < np; t++)
                match |= (b == (int)exact[t]);
            if (run >= 12 && match) {
                int start = i - 11;
                total += i + 1 - (start > last ? start : last);
                last = i + 1;
            }
        }
        out[0] = out[1] = total;
        return;
    }
    for (int i = 0; i < n; i++) {
        int b = base(s[i]);
        if (b < 0) {
            valid = 0;
            w = 0;
            continue;
        }
        w = (w << 2) | (unsigned)b;
        if (valid < 2 * k)
            valid++;
        if (valid < k)
            continue;
        /* exact units add to the closed span and the screen count */
        uint64_t a = w & mask;
        int hit = 0;
        for (int t = 0; t < np; t++)
            if (a == exact[t]) {
                hit = 1;
                break;
            }
        if (hit) {
            int start = i - k + 1;
            closed += i + 1 - (start > end ? start : end);
            if (end >= 0 && start > end && start - end <= 3)
                closed += start - end;
            end = i + 1;
            int delta = i + 1 - last;
            covered = delta >= 64 ? 0 : covered << delta;
            uint64_t mm = (UINT64_C(1) << k) - 1;
            total += __builtin_popcountll(mm & ~covered);
            covered |= mm;
            last = i + 1;
        }
        /* two adjacent units at most one substitution apart, one of them exact */
        if (valid < 2 * k)
            continue;
        uint64_t z = (w >> (2 * k)) & mask;
        if (distance(a, z) <= 1) {
            hit = 0;
            for (int t = 0; t < np; t++)
                if (a == exact[t] || z == exact[t]) {
                    hit = 1;
                    break;
                }
            if (hit) {
                int delta = i + 1 - last;
                covered = delta >= 64 ? 0 : covered << delta;
                uint64_t mm = (UINT64_C(1) << (2 * k)) - 1;
                total += __builtin_popcountll(mm & ~covered);
                covered |= mm;
                last = i + 1;
            }
        }
    }
    out[0] = total;
    out[1] = closed;
}

/* Tract bp of one motif: out[0] in all of s, out[1] in its first 12 bp, out[2] in its last 12.
   Period 1 needs exact units; longer periods use the scan() tract rule. */
void target_tract(const char *s, int n, const char *patterns, int k, int np, int *out) {
    uint64_t exact[16];
    for (int j = 0; j < np; j++) {
        exact[j] = 0;
        for (int t = 0; t < k; t++)
            exact[j] = (exact[j] << 2) | base(patterns[j * k + t]);
    }
    int copies = (12 + k - 1) / k;
    if (copies < 3)
        copies = 3;
    int width = copies * k, valid = 0;
    uint64_t w = 0, mask = (UINT64_C(1) << (2 * k)) - 1;
    unsigned char *covered = calloc(n ? n : 1, 1);
    if (!covered)
        abort();
    for (int i = 0; i < n; i++) {
        int b = base(s[i]);
        if (b < 0) {
            valid = 0;
            w = 0;
            continue;
        }
        w = (w << 2) | (unsigned)b;
        if (valid < width)
            valid++;
        if (valid < width)
            continue;
        int hit = 0;
        for (int t = 0; t < np && !hit; t++) {
            int last = 0, ok = 1;
            for (int j = 0; j < copies; j++) {
                int d = distance((w >> (2 * k * j)) & mask, exact[t]);
                if ((k == 1 && d) || (k > 1 && (d > 1 || (j && d + last > 1)))) {
                    ok = 0;
                    break;
                }
                last = d;
            }
            hit = ok;
        }
        if (hit)
            memset(covered + i - width + 1, 1, width);
    }
    out[0] = out[1] = out[2] = 0;
    for (int i = 0; i < n; i++) {
        out[0] += covered[i];
        if (i < 12)
            out[1] += covered[i];
        if (i >= n - 12)
            out[2] += covered[i];
    }
    free(covered);
}
