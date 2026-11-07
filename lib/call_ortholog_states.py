#!/usr/bin/env python3
"""Call comparator intron states (NO_INTRON, TELOTRON, MIXED, ORDINARY, AMBIGUOUS) past gates."""

import argparse
import collections
import csv
import os
import statistics

from _common import flank_gate, intron_observed, its_spans, locus_stub, read_aln, repeat_spans
from paired_repeat_tracts import tracts

WINDOW_AA = 10  # CIWOG simRegion
MIN_SIM = 0.10  # CIWOG minSim, to call an intron present
MIN_ABSENT_SIM = 0.30  # CIWOG minAbsentSim, higher bar to call one absent
GATE_SLOP = 2  # residues either side of the junction where an indel is tolerated
NONTELO_MAX = 0.30  # repeat purity below this is ORDINARY
MIN_ITS_UNITS = 4  # ITS: >=4 units, <=1 substitution per adjacent two-unit pair
REPEAT_POLICY = "exact_or_1sub_per_2units_v1"
RBH_PASS = ("1", "True", "true")
COLUMNS = [
    "telotron",
    "focal",
    "ortholog",
    "state",
    "aligner_call",
    "local_ident",
    "protein_ident",
    "rbh",
    "synteny_frac",
    "repeat_purity",
    "longest_array_bp",
    "intron_bp",
    "repeat_policy",
]


def measure_repeat_profile(seq, motifs):
    """(covered fraction, longest single-phase tract bp, covered bp) over the given motifs."""
    if not seq:
        return 0.0, 0, 0
    labels = repeat_spans(seq.upper(), motifs)
    covered = sum(label in 'FfRr' for label in labels)

    # longest tract over all motifs, in bp
    longest = max(
        (end - start for motif in motifs for start, end, _, _ in tracts(seq, motif, min_units=2)),
        default=0,
    )
    return covered / len(seq), longest, covered


def repeat_profile(seq, motifs):
    """Profile of the single best motif by coverage, then longest tract."""
    profiles = (measure_repeat_profile(seq, [motif]) for motif in motifs)
    return max(profiles, key=lambda profile: (profile[0], profile[1]))


def synteny_value(row):
    """Synteny fraction when computed, else '' (absent or unparsable)."""
    value = (row.get("synteny_frac") or "").strip()
    if value in ("", "-"):
        return ""
    try:
        return float(value)
    except ValueError:
        return ""


def intron_state(seq, motifs, telo_purity, telo_min_array):
    """(state, purity, longest tract) of an observed comparator intron."""
    purity, longest, _ = repeat_profile(seq, motifs)
    if purity > telo_purity and longest >= telo_min_array:
        state = "TELOTRON"
    elif any(its_spans(seq, motif, MIN_ITS_UNITS) for motif in motifs):
        # ITS inside an otherwise ordinary intron (published ITS criterion)
        state = "MIXED"
    elif purity < NONTELO_MAX:
        state = "ORDINARY"
    else:
        state = "AMBIGUOUS"
    return state, purity, longest


def call_row(row, window_gates, intron_seqs, focal, motifs, args):
    """(state record, None) for one comparator row, or (None, drop reason)."""
    if row.get("frameshift_at_junction") == "1":
        return None, "frameshift or stop at junction"
    ortholog = row["ortholog"]
    window = window_gates.get(ortholog)
    if window is None:
        return None, "not in protein MSA"

    # use ok, not n_indel: ok also fails an incomplete window
    window_ok, local_id, _ = window
    if not window_ok:
        return None, "indel or incomplete local window"
    if row.get("rbh_pass") not in RBH_PASS:
        return None, "failed RBH"

    # ABSENT becomes NO_INTRON, an observed intron is classed by repeat content
    aligned = intron_seqs.get(ortholog, "")
    seq = aligned.replace("-", "")
    if row["call"] == "ABSENT":
        if local_id < MIN_ABSENT_SIM:
            return None, "below minAbsentSim"
        state, purity, longest = "NO_INTRON", "", ""
    elif intron_observed(row, aligned):
        if local_id < MIN_SIM:
            return None, "below minSim"
        state, purity, longest = intron_state(seq, motifs, args.telo_purity, args.telo_min_array)
    else:
        return None, "no measurable intron (NO_DATA)"

    record = dict(
        telotron=row["telotron"],
        focal=focal,
        ortholog=ortholog,
        state=state,
        repeat_policy=REPEAT_POLICY,
        aligner_call=row["call"],
        local_ident=round(local_id, 4),
        protein_ident=row.get("pct_identity", ""),
        rbh=row.get("rbh_pass", ""),
        synteny_frac=synteny_value(row),
        repeat_purity="" if purity == "" else round(purity, 4),
        longest_array_bp="" if purity == "" else longest,
        intron_bp=len(seq) if seq else "",
    )
    return record, None


def report(records, dropped, identities):
    """Console report: gates, state counts, drop reasons, identity ranges."""
    loci = {record["telotron"] for record in records}
    print("ORTHOLOGY: RBH required=True, synteny floor=off (applied only where computed)")
    print(
        f"POSITION:  +/-{WINDOW_AA} aa, minSim {MIN_SIM}, "
        f"minAbsentSim {MIN_ABSENT_SIM}, no indel outside +/-{GATE_SLOP}"
    )
    print(f"\nstates: {len(records)} over {len(loci)} loci")
    for state, n in collections.Counter(record["state"] for record in records).most_common():
        print(f"   {state:<12}{n:>6}")
    print("\ndropped:")
    for reason, n in dropped.most_common():
        print(f"   {reason:<34}{n:>7}")

    with_synteny = sum(1 for record in records if record["synteny_frac"] != "")
    print(
        f"\nsynteny available on {with_synteny}/{len(records)} kept rows "
        f"({100 * with_synteny / max(1, len(records)):.1f}%)"
    )

    if identities:
        identities.sort()
        print(
            f"protein identity of KEPT rows (reported, not gated): "
            f"min {identities[0]:.2f}  q1 {identities[len(identities)//4]:.2f}  "
            f"median {statistics.median(identities):.2f}  max {identities[-1]:.2f}"
        )

    # zero kept rows is valid and must not fail after the output is written
    local = sorted(record["local_ident"] for record in records)
    if local:
        print(
            f"local window identity: min {local[0]:.2f}  q1 {local[len(local)//4]:.2f}  "
            f"median {statistics.median(local):.2f}"
        )
    else:
        print("local window identity: no rows kept")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ortho-dir", required=True)
    ap.add_argument("--focal", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--motifs", required=True, help="comma-separated focal-lineage motifs")
    ap.add_argument(
        "--raw-dna",
        action="store_true",
        help="read intron_dna_raw/<locus>.fasta instead of intron_dna/<locus>.aln.fasta",
    )
    ap.add_argument("--telo-purity", type=float, default=0.50)
    ap.add_argument(
        "--telo-min-array",
        type=int,
        default=55,
        help="longest contiguous array (bp) a TELOTRON must carry",
    )
    ap.add_argument(
        "--min-focal-purity",
        type=float,
        default=0.90,
        help="focal telotron purity required to use the locus",
    )
    args = ap.parse_args()
    motifs = args.motifs.split(",")
    directory = args.ortho_dir

    # ortholog_map rows grouped by telotron locus
    rows_by_locus = collections.defaultdict(list)
    with open(f"{directory}/ortholog_map.tsv") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows_by_locus[row["telotron"]].append(row)

    # per locus: protein MSA and focal-purity gates, then each comparator row
    records = []
    dropped = collections.Counter()
    identities = []
    for tid, rows in rows_by_locus.items():
        msa_path = f"{directory}/protein_msa/{args.focal}__{rows[0]['protein']}.aln.fasta"
        if not os.path.exists(msa_path):
            dropped["no protein MSA"] += len(rows)
            continue

        # intron DNA per genome (raw or aligned), keyed by genome ID
        if args.raw_dna:
            dna_path = f"{directory}/intron_dna_raw/{locus_stub(tid)}.fasta"
        else:
            dna_path = f"{directory}/intron_dna/{locus_stub(tid)}.aln.fasta"
        intron_seqs = {}
        if os.path.exists(dna_path):
            intron_seqs = {name.split("|")[0]: seq for name, seq in read_aln(dna_path)}
        focal_seq = intron_seqs.get(args.focal, "").replace("-", "")
        if repeat_profile(focal_seq, motifs)[0] < args.min_focal_purity:
            dropped["focal below purity floor"] += len(rows)
            continue

        # min_ident 0: local identity is only read here and gated in call_row
        residue = int(rows[0]["intron_residue"] or 0)
        msa = read_aln(msa_path)
        window_gates = flank_gate(msa, args.focal, residue, WINDOW_AA, 0.0, GATE_SLOP)

        for row in rows:
            record, reason = call_row(row, window_gates, intron_seqs, args.focal, motifs, args)
            if record is None:
                dropped[reason] += 1
                continue

            # identity is reported for kept rows, not gated
            try:
                identities.append(float(row["pct_identity"]))
            except (TypeError, ValueError):
                pass
            records.append(record)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)

    report(records, dropped, identities)
    print(f"\n  wrote {args.out}")


if __name__ == "__main__":
    main()
