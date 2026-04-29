#!/usr/bin/env python3
"""
Build per-example DNA-letter MSA figures for §9f linker-origin findings.

Wrapped layout: 80 columns per row, multiple wrap-blocks stacked vertically.
"""

import csv, glob, json, subprocess, tempfile
from pathlib import Path
from collections import defaultdict, Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path("/scratch1/alex/telotrons")
PEUK = ROOT / "pan_euk_telotrons"
OUT_DIR = PEUK / "real_telotrons" / "eimeria_linker_blast" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

ACC_TO_NAME = {
    "GCF_000499385.1": "E. necatrix",
    "GCF_000499425.1": "E. acervulina",
    "GCF_000499545.2": "E. tenella",
    "GCF_000499605.1": "E. maxima",
    "GCF_000499745.2": "E. mitis",
    "GCF_002999335.1": "Cyclospora cayetanensis",
    "GCF_000006565.2": "Toxoplasma gondii",
    "GCF_002563875.1": "Besnoitia besnoiti",
    "GCF_000208865.1": "Neospora caninum",
    "GCA_000727475.1": "Sarcocystis neurona",
    "GCF_000006425.1": "Cryptosporidium hominis",
    "GCF_900005855.1": "Plasmodium gallinaceum",
}

FWD = {"TTAGGG", "TTTAGGG"}
REV_TBL = str.maketrans("ACGT", "TGCA")
ALL_TELO = set()
for u in FWD:
    for i in range(len(u)):
        ALL_TELO.add(u[i:] + u[:i])
    rev = u[::-1].translate(REV_TBL)
    for i in range(len(rev)):
        ALL_TELO.add(rev[i:] + rev[:i])


def telo_mask_seq(seq):
    s = seq.upper()
    n = len(s)
    m = [0] * n
    i = 0
    while i < n:
        hit = False
        for k in ALL_TELO:
            if s[i:i+len(k)] == k:
                for j in range(i, min(i+len(k), n)):
                    m[j] = 1
                i += len(k); hit = True; break
        if not hit:
            i += 1
    return m


def telo_mask_aligned(aln_seq):
    ung = aln_seq.replace("-", "")
    m_ung = telo_mask_seq(ung)
    m_aln = []
    idx = 0
    for c in aln_seq:
        if c == "-":
            m_aln.append(0)
        else:
            m_aln.append(m_ung[idx])
            idx += 1
    return m_aln


_genome_cache = {}

def get_contig_seq(acc, contig):
    if (acc, contig) in _genome_cache:
        return _genome_cache[(acc, contig)]
    candidates = list((PEUK / "genomes" / acc).rglob("*.fna")) + \
                 list((PEUK / "genomes_extra" / acc).rglob("*.fna"))
    if not candidates:
        return None
    fna = candidates[0]
    cur = None
    parts = []
    target = None
    with open(fna) as f:
        for line in f:
            if line.startswith(">"):
                if cur == contig:
                    target = "".join(parts).upper()
                    break
                cur = line[1:].rstrip().split()[0]
                parts = []
            else:
                parts.append(line.rstrip())
        if target is None and cur == contig:
            target = "".join(parts).upper()
    _genome_cache[(acc, contig)] = target
    return target


def parse_q(qseqid):
    p = qseqid.split("__")
    acc = p[0]; contig = p[1]; coords = p[2]
    start, end = coords.split("-")
    return acc, contig, int(start), int(end)


def pad_region(acc, contig, start, end, pad=200):
    seq = get_contig_seq(acc, contig)
    if seq is None:
        return None, None, None
    lo = max(0, start - pad)
    hi = min(len(seq), end + pad)
    return seq[lo:hi], start - lo, hi - end


def revcomp(s):
    return s[::-1].translate(REV_TBL)


def mafft_align(seqs_dict):
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        in_fa = td / "in.fa"
        with open(in_fa, "w") as f:
            for name, seq in seqs_dict.items():
                clean = name.replace("|", "_").replace(":", "_").replace(" ", "_")
                f.write(f">{clean}\n{seq}\n")
        r = subprocess.run(
            ["mafft", "--quiet", "--auto", str(in_fa)],
            capture_output=True, text=True, timeout=120,
        )
        out = {}
        cur, parts = None, []
        for line in r.stdout.split("\n"):
            if line.startswith(">"):
                if cur:
                    out[cur] = "".join(parts).upper()
                cur = line[1:].strip()
                parts = []
            else:
                parts.append(line.strip())
        if cur:
            out[cur] = "".join(parts).upper()
    name_map = {n.replace("|", "_").replace(":", "_").replace(" ", "_"): n
                for n in seqs_dict}
    return {name_map.get(k, k): v for k, v in out.items()}


def map_ungapped_to_aligned(aln_seq, ung_start, ung_end):
    a_idx = []
    for i, c in enumerate(aln_seq):
        if c != "-":
            a_idx.append(i)
    if ung_start >= len(a_idx):
        ung_start = len(a_idx) - 1
    if ung_end > len(a_idx):
        ung_end = len(a_idx)
    return a_idx[ung_start], a_idx[min(ung_end - 1, len(a_idx) - 1)] + 1


def render_msa_wrapped(title, subtitle, alignment, telo_span_host, hit_region_host,
                       host_label, sister_label, out_path, annotation_text="",
                       cols_per_row=80, base_offset=0):
    """Wrapped MSA rendering. cols_per_row columns per wrap block."""
    seq_names = list(alignment.keys())
    aln_len = len(next(iter(alignment.values())))

    # Consensus
    cols_data = [[alignment[n][i] for n in seq_names] for i in range(aln_len)]
    consensus = []
    for col in cols_data:
        vals = [c for c in col if c != "-"]
        consensus.append(Counter(vals).most_common(1)[0][0] if vals else "-")

    # Telo masks per row
    telo_mask_row = {n: telo_mask_aligned(alignment[n]) for n in seq_names}

    # Colors
    base_telo_match     = "#bce4bc"
    base_telo_mismatch  = "#fdd7a8"
    base_link_match     = "#e8e8e8"
    base_link_mismatch  = "#ffc8c8"

    n_seqs = len(seq_names)
    n_blocks = (aln_len + cols_per_row - 1) // cols_per_row

    # Layout dims (inches)
    cell_w = 0.10
    cell_h = 0.20
    block_h = (n_seqs + 2) * cell_h + 0.10  # +2 for telotron-bar above + coord ruler below
    fig_w = cols_per_row * cell_w + 4.0     # 4 inches for left labels
    fig_h = 1.4 + 0.4 + n_blocks * block_h + 1.0  # title + spacing + blocks + legend

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, fig_h)
    ax.invert_yaxis()
    ax.axis("off")

    # Title + subtitle
    ax.text(fig_w / 2, 0.25, title, ha="center", va="top", fontsize=12, fontweight="bold")
    ax.text(fig_w / 2, 0.55, subtitle, ha="center", va="top", fontsize=8.5, color="#444")
    if annotation_text:
        ax.text(fig_w / 2, 0.80, annotation_text, ha="center", va="top",
                fontsize=8.5, color="#222", style="italic")

    label_x = 0.05
    aln_x0 = 3.5  # left margin for species labels (inches)

    y_cursor = 1.20
    for blk in range(n_blocks):
        col_start = blk * cols_per_row
        col_end = min((blk + 1) * cols_per_row, aln_len)

        # Telotron-span bar (red) atop block 0
        bar_y = y_cursor + 0.05
        ts_lo, ts_hi = telo_span_host
        # Show portion of telotron span inside this block
        ts_lo_blk = max(ts_lo, col_start)
        ts_hi_blk = min(ts_hi, col_end)
        if ts_lo_blk < ts_hi_blk:
            x_lo = aln_x0 + (ts_lo_blk - col_start) * cell_w
            x_hi = aln_x0 + (ts_hi_blk - col_start) * cell_w
            ax.add_patch(Rectangle((x_lo, bar_y), x_hi - x_lo, 0.10,
                                   facecolor="#d62728", alpha=0.85, edgecolor="none"))

        # Each species row
        for ri, name in enumerate(seq_names):
            seq = alignment[name]
            tm = telo_mask_row[name]
            row_y = y_cursor + 0.20 + ri * cell_h
            # Label
            ax.text(aln_x0 - 0.05, row_y + cell_h / 2, name,
                    ha="right", va="center", fontsize=8, family="monospace")
            # Cells
            for ci_local, ci in enumerate(range(col_start, col_end)):
                x = aln_x0 + ci_local * cell_w
                base = seq[ci]
                cons = consensus[ci]
                is_match = (base != "-" and cons != "-" and base == cons)
                is_gap = (base == "-")
                is_telo = bool(tm[ci])
                if is_gap:
                    facecolor = "white"; edgecolor = "#dddddd"; lw = 0.2
                elif is_telo and is_match:
                    facecolor = base_telo_match; edgecolor = "none"; lw = 0
                elif is_telo and not is_match:
                    facecolor = base_telo_mismatch; edgecolor = "none"; lw = 0
                elif is_match:
                    facecolor = base_link_match; edgecolor = "none"; lw = 0
                else:
                    facecolor = base_link_mismatch; edgecolor = "none"; lw = 0
                ax.add_patch(Rectangle((x, row_y), cell_w, cell_h,
                                       facecolor=facecolor, edgecolor=edgecolor, linewidth=lw))
                if base != "-":
                    color = "black" if is_match else "#aa0000"
                    weight = "normal" if is_match else "bold"
                    ax.text(x + cell_w / 2, row_y + cell_h / 2, base,
                            ha="center", va="center", fontsize=5.5,
                            family="monospace", color=color, weight=weight)

        # Coordinate ruler (alignment col)
        ruler_y = y_cursor + 0.20 + n_seqs * cell_h + 0.04
        # Tick every 10
        for ci in range(col_start, col_end, 10):
            x = aln_x0 + (ci - col_start) * cell_w
            ax.text(x + cell_w / 2, ruler_y, str(ci),
                    ha="center", va="top", fontsize=6.5, color="#666")

        # End-of-block coordinates per row (ungapped)
        for ri, name in enumerate(seq_names):
            seq = alignment[name]
            row_y = y_cursor + 0.20 + ri * cell_h
            # ungapped count up to col_end
            ung_count = sum(1 for c in seq[:col_end] if c != "-")
            ax.text(aln_x0 + (col_end - col_start) * cell_w + 0.05,
                    row_y + cell_h / 2, str(ung_count),
                    ha="left", va="center", fontsize=6.5, color="#888",
                    family="monospace")

        y_cursor += block_h

    # Legend
    legend_y = y_cursor + 0.10
    legend_items = [
        (base_link_match, "match (linker/non-telomeric)"),
        (base_link_mismatch, "mismatch (linker)"),
        (base_telo_match, "match (telomeric kmer)"),
        (base_telo_mismatch, "mismatch (telomeric)"),
    ]
    x = 0.5
    for color, label in legend_items:
        ax.add_patch(Rectangle((x, legend_y), 0.25, 0.20,
                               facecolor=color, edgecolor="#888888", linewidth=0.4))
        ax.text(x + 0.30, legend_y + 0.10, label, ha="left", va="center", fontsize=8)
        x += 0.30 + len(label) * 0.07 + 0.4
    ax.add_patch(Rectangle((x, legend_y), 0.25, 0.10,
                           facecolor="#d62728", alpha=0.85, edgecolor="none"))
    ax.text(x + 0.30, legend_y + 0.10, "telotron span (host)", ha="left", va="center", fontsize=8)

    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()
    print(f"  wrote {out_path}")


def determine_strand_via_blast(host_telo_seq, sister_seq):
    """Run blastn host telotron vs sister to figure out which strand sister hit is on."""
    if len(host_telo_seq) < 30 or len(sister_seq) < 30:
        return "+"
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td/"q.fa").write_text(f">q\n{host_telo_seq}\n")
        (td/"s.fa").write_text(f">s\n{sister_seq}\n")
        r = subprocess.run(
            ["blastn", "-task", "blastn", "-query", str(td/"q.fa"),
             "-subject", str(td/"s.fa"), "-evalue", "1",
             "-outfmt", "6 pident length sstart send"],
            capture_output=True, text=True, timeout=30,
        )
    best_strand = "+"
    best_score = 0
    for line in r.stdout.strip().split("\n"):
        if not line: continue
        p = line.split("\t")
        score = float(p[0]) * int(p[1])
        if score > best_score:
            best_score = score
            best_strand = "-" if int(p[2]) > int(p[3]) else "+"
    return best_strand


def build_example(ex, category, ex_idx):
    q = ex["q"]
    h = ex["best"]
    qacc, qcontig, qstart, qend = parse_q(q)
    PAD = 200

    host_seq, host_pad_left, host_pad_right = pad_region(qacc, qcontig, qstart, qend, PAD)
    if host_seq is None:
        print(f"  [{q}] host genome missing"); return False

    sacc = h["sacc"]; scontig = h["contig"]
    s_lo, s_hi = h["sstart"], h["send"]
    s_seq, s_pad_left, s_pad_right = pad_region(sacc, scontig, s_lo, s_hi, PAD)
    if s_seq is None:
        print(f"  [{q}] sister genome missing"); return False

    # Determine strand
    host_telo = host_seq[host_pad_left:host_pad_left + (qend - qstart)]
    s_strand = determine_strand_via_blast(host_telo, s_seq)
    if s_strand == "-":
        s_seq = revcomp(s_seq)

    host_label = f"{ACC_TO_NAME.get(qacc, qacc):<25s} (HOST)"
    sister_label = f"{ACC_TO_NAME.get(sacc, sacc):<25s} (sister)"

    seqs = {host_label: host_seq, sister_label: s_seq}
    aln = mafft_align(seqs)
    if host_label not in aln or sister_label not in aln:
        print(f"  [{q}] mafft failed"); return False

    host_aln = aln[host_label]
    teln_aln_start, teln_aln_end = map_ungapped_to_aligned(
        host_aln, host_pad_left, host_pad_left + (qend - qstart))

    title = f"Category {category}, Example {ex_idx+1}"
    subtitle = (
        f"Host: {qacc} {qcontig}:{qstart}–{qend} "
        f"(telotron, {qend-qstart} bp)    "
        f"Sister hit: {sacc} {scontig}:{s_lo}–{s_hi} ({s_strand} strand)    "
        f"BLAST: {h['pident']}% identity over {h['length']} bp"
    )
    annotation = {
        "A": "Sister-species region is plain genomic DNA (no telotron at this locus) — the linker preserves ancestral coccidian DNA into which telomeric arrays were inserted.",
        "B": "Sister-species region is itself a telotron — the linker is shared between telotrons across species (recurrent capture or shared origin).",
        "C": "Sarcocystidae outgroup (~150 My to Eimeria) — deeply conserved sequence shared with Toxoplasma/Neospora/Besnoitia/Sarcocystis.",
        "D": "Plasmodium / Cryptosporidium hit (>400 My). Likely conserved coding-flank context near the telotron rather than novel linker DNA.",
    }.get(category, "")

    aln_ord = {host_label: host_aln, sister_label: aln[sister_label]}
    out_path = OUT_DIR / f"cat{category}_example{ex_idx+1}.png"
    render_msa_wrapped(
        title, subtitle, aln_ord,
        telo_span_host=(teln_aln_start, teln_aln_end),
        hit_region_host=None,
        host_label=host_label, sister_label=sister_label,
        out_path=out_path,
        annotation_text=annotation,
        cols_per_row=100,
    )
    return True


def main():
    ex_data = json.load(open("/tmp/examples.json"))
    for cat_key, examples in ex_data.items():
        category = cat_key.split("_")[0]
        print(f"\n=== Category {category} ===")
        for i, ex in enumerate(examples[:3]):
            print(f"  example {i+1}: {ex['q']}")
            try:
                build_example(ex, category, i)
            except Exception as e:
                print(f"  ERROR: {e}")
                import traceback; traceback.print_exc()


if __name__ == "__main__":
    main()
