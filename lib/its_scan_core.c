/* Telomeric tract scanner: maximal runs of at least four whole motif units in one phase, with at
 * most one substitution in any two adjacent units and at least one exact unit.
 * Coordinates are zero-based half-open. */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define MAX_PATTERNS 32
#define MAX_UNIT 8

typedef struct {
    int64_t start, end;
    int32_t pattern, mismatches, exact_units;
} Tract;

typedef struct {
    Tract *tracts;
    int64_t count, capacity;
    int failed;
} Result;

/* the open run of one (pattern, phase) */
typedef struct {
    int64_t last_unit, start;
    int32_t units, mismatches, exact_units, last_mismatches;
} Run;

static int base_code(unsigned char base) {
    switch (base) {
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

/* report a run of at least four units with at least one exact unit */
static void finish(Result *result, Run *run, int unit, int pattern) {
    if (result->failed || run->units < 4 || run->exact_units == 0)
        return;

    /* grow the tract array by doubling */
    if (result->count == result->capacity) {
        int64_t capacity = result->capacity ? result->capacity * 2 : 256;
        Tract *grown = realloc(result->tracts, (size_t)capacity * sizeof(Tract));
        if (!grown) {
            result->failed = 1;
            return;
        }
        result->tracts = grown;
        result->capacity = capacity;
    }
    result->tracts[result->count++] = (Tract){
        run->start, run->last_unit + unit, pattern, run->mismatches, run->exact_units};
}

/* add one unit at position. A gap, or >1 substitution across the adjacent pair, restarts the run */
static void extend(Result *result, Run *run, int64_t position, int unit, int pattern,
                   int mismatches) {
    int adjacent = run->units && run->last_unit + unit == position;
    if (!adjacent || run->last_mismatches + mismatches > 1) {
        finish(result, run, unit, pattern);
        run->start = position;
        run->units = 0;
        run->mismatches = 0;
        run->exact_units = 0;
    }
    run->last_unit = position;
    run->last_mismatches = mismatches;
    run->units++;
    run->mismatches += mismatches;
    run->exact_units += mismatches == 0;
}

int64_t scan_its(const char *sequence, int64_t length, const char *patterns, int unit,
                 Tract **out) {
    *out = NULL;
    if (unit < 2 || unit > MAX_UNIT || length < 0)
        return -2;
    int n_patterns = 2 * unit;
    uint32_t n_words = 1u << (2 * unit);
    uint32_t word_mask = n_words - 1;

    /* per 2-bit encoded word: bit p set when the word equals pattern p (exact) or is within
     * one substitution of it (near) */
    uint32_t *exact = calloc(n_words, sizeof(uint32_t));
    uint32_t *near = calloc(n_words, sizeof(uint32_t));
    if (!exact || !near) {
        free(exact);
        free(near);
        return -1;
    }
    for (int p = 0; p < n_patterns; p++) {
        uint32_t word = 0;
        for (int j = 0; j < unit; j++) {
            int code = base_code(patterns[p * unit + j]);
            if (code < 0) {
                free(exact);
                free(near);
                return -2;
            }
            word = (word << 2) | code;
        }
        exact[word] |= 1u << p;
        near[word] |= 1u << p;
        for (int j = 0; j < unit; j++) {
            for (int code = 0; code < 4; code++) {
                uint32_t variant = (word & ~(3u << (2 * j))) | ((uint32_t)code << (2 * j));
                near[variant] |= 1u << p;
            }
        }
    }

    Run runs[MAX_PATTERNS][MAX_UNIT];
    memset(runs, 0, sizeof(runs));
    Result result = {0};
    uint32_t word = 0;
    int valid = 0;

    for (int64_t i = 0; i < length; i++) {
        /* a non-ACGT base breaks every word that contains it */
        int code = base_code(sequence[i]);
        if (code < 0) {
            valid = 0;
            word = 0;
            continue;
        }
        word = ((word << 2) | code) & word_mask;
        if (valid < unit)
            valid++;
        if (valid < unit)
            continue;

        /* extend the run of every pattern within one substitution, in this phase */
        int64_t position = i - unit + 1;
        int phase = (int)(position % unit);
        for (uint32_t bits = near[word]; bits; bits &= bits - 1) {
            int p = __builtin_ctz(bits);
            int mismatches = (exact[word] & (1u << p)) ? 0 : 1;
            extend(&result, &runs[p][phase], position, unit, p, mismatches);
        }
    }

    /* close every open run */
    for (int p = 0; p < n_patterns; p++)
        for (int phase = 0; phase < unit; phase++)
            finish(&result, &runs[p][phase], unit, p);

    free(exact);
    free(near);
    if (result.failed) {
        free(result.tracts);
        return -1;
    }
    *out = result.tracts;
    return result.count;
}

void free_its(Tract *tracts) { free(tracts); }
