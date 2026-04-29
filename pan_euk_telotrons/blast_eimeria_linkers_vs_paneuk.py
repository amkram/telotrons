#!/usr/bin/env python3
"""
BLAST 759 pure non-repeat Eimeria linkers against the entire pan-eukaryotic
genome collection on disk. Goal: find any genomic source for these linkers,
or confirm they are unique to host telotron context.

Pipeline:
  1. Build single concatenated FNA (208 GB) with accession-tagged seqids.
  2. makeblastdb (-parse_seqids).
  3. Extract pure non-repeat Eimeria linkers from strict_linkers_eimeria.json.
  4. blastn vs the combined db. Filter for ≥80% pident, ≥30 bp.
  5. Classify hits as self-Eimeria / other-Apicomplexa / other-eukaryote / Tara MAG.
"""

import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path("/scratch1/alex/telotrons")
PEUK = ROOT / "pan_euk_telotrons"
WORK = PEUK / "real_telotrons" / "eimeria_linker_blast"
WORK.mkdir(parents=True, exist_ok=True)

LINKERS_JSON = PEUK / "real_telotrons" / "strict_linkers_eimeria.json"
COMBINED_FA = WORK / "apicomplexa_genomes.fa"
DB_PREFIX = WORK / "apicomplexa_db"

# Restricted set: only close relatives of Eimeria (Apicomplexa) — 12 genomes.
# Excludes distant outgroups to avoid diluting signal with cross-kingdom hits.
APICOMPLEXA_ACCS = {
    # Eimeria (host species — used as self-hit reference, NOT excluded from db)
    "GCF_000499385.1": "Eimeria_necatrix",
    "GCF_000499425.1": "Eimeria_acervulina",
    "GCF_000499545.2": "Eimeria_tenella",
    "GCF_000499605.1": "Eimeria_maxima",
    "GCF_000499745.2": "Eimeria_mitis",
    # Eimeriidae sister genus
    "GCF_002999335.1": "Cyclospora_cayetanensis",
    # Sarcocystidae (next family within Coccidia)
    "GCF_000006565.2": "Toxoplasma_gondii",
    "GCF_002563875.1": "Besnoitia_besnoiti",
    "GCF_000208865.1": "Neospora_caninum",
    "GCA_000727475.1": "Sarcocystis_neurona",
    # Other Apicomplexa (deeper outgroups within phylum)
    "GCF_000006425.1": "Cryptosporidium_hominis",
    "GCF_900005855.1": "Plasmodium_gallinaceum",
}

GENOME_DIRS_CANDIDATES = [
    PEUK / "genomes",
    PEUK / "genomes_extra",
]

EIMERIA_ACCS = {acc for acc, sp in APICOMPLEXA_ACCS.items() if sp.startswith("Eimeria_")}


FWD = {"TTAGGG", "TTTAGGG"}
REV_TBL = str.maketrans("ACGT", "TGCA")
def rotations(unit):
    return {unit[i:] + unit[:i] for i in range(len(unit))}
ALL_TELO = set()
for u in FWD:
    ALL_TELO |= rotations(u)
    ALL_TELO |= rotations(u[::-1].translate(REV_TBL))


def build_combined_fasta():
    if COMBINED_FA.exists() and COMBINED_FA.stat().st_size > 100_000_000:
        print(f"[build_combined_fasta] reusing {COMBINED_FA} ({COMBINED_FA.stat().st_size/1e6:.1f} MB)")
        return

    print(f"[build_combined_fasta] gathering {len(APICOMPLEXA_ACCS)} Apicomplexa genomes...")
    n_contigs = 0
    n_genomes = 0
    with open(COMBINED_FA, "w") as out:
        for acc, sp_name in APICOMPLEXA_ACCS.items():
            fna_path = None
            for d in GENOME_DIRS_CANDIDATES:
                cand = d / acc
                if cand.is_dir():
                    fnas = list(cand.rglob("*.fna"))
                    if fnas:
                        fna_path = fnas[0]
                        break
            if fna_path is None:
                print(f"[build_combined_fasta] MISSING: {acc} ({sp_name})", file=sys.stderr)
                continue
            n_genomes += 1
            with open(fna_path) as f:
                for line in f:
                    if line.startswith(">"):
                        n_contigs += 1
                        orig_id = line[1:].rstrip().split()[0]
                        out.write(f">{acc}::{orig_id}\n")
                    else:
                        out.write(line)
            print(f"  + {acc} ({sp_name}) -> {fna_path.name}")
    sz = COMBINED_FA.stat().st_size
    print(f"[build_combined_fasta] wrote {COMBINED_FA} ({sz/1e6:.1f} MB, {n_genomes} genomes, {n_contigs} contigs)")


def build_blast_db():
    nsq = Path(str(DB_PREFIX) + ".nsq")
    nhr = Path(str(DB_PREFIX) + ".nhr")
    if nsq.exists() and nhr.exists():
        print(f"[build_blast_db] reusing {DB_PREFIX}")
        return
    print("[build_blast_db] running makeblastdb (will take a while)...")
    subprocess.run(
        ["makeblastdb", "-in", str(COMBINED_FA), "-dbtype", "nucl",
         "-parse_seqids", "-out", str(DB_PREFIX)],
        check=True,
    )
    print(f"[build_blast_db] done")


def extract_linkers():
    out_fa = WORK / "eimeria_pure_linkers.fa"
    if out_fa.exists() and out_fa.stat().st_size > 1000:
        print(f"[extract_linkers] reusing {out_fa}")
        return out_fa

    d = json.load(open(LINKERS_JSON))
    n_pure = 0
    with open(out_fa, "w") as f:
        for r in d:
            seq = (r.get("linker_seq") or "").upper()
            if len(seq) < 25:
                continue
            if any(k in seq for k in ALL_TELO):
                continue
            rid = f"{r['acc']}__{r['contig']}__{r['start']}-{r['end']}__len{len(seq)}"
            f.write(f">{rid}\n{seq}\n")
            n_pure += 1
    print(f"[extract_linkers] wrote {out_fa} ({n_pure} pure-non-repeat linkers ≥25bp)")
    return out_fa


def run_blast(query_fa, n_threads=64):
    out_tsv = WORK / "eimeria_linkers_vs_apicomplexa.blastn.tsv"
    if out_tsv.exists() and out_tsv.stat().st_size > 0:
        print(f"[run_blast] {out_tsv} exists ({out_tsv.stat().st_size} bytes); reusing")
        return out_tsv
    print(f"[run_blast] BLASTing {query_fa} vs {DB_PREFIX} on {n_threads} threads...")
    fmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    subprocess.run(
        ["blastn", "-task", "blastn",
         "-query", str(query_fa),
         "-db", str(DB_PREFIX),
         "-evalue", "1e-3",
         "-perc_identity", "70",
         "-word_size", "11",
         "-dust", "no",
         "-num_threads", str(n_threads),
         "-max_target_seqs", "20",
         "-outfmt", fmt,
         "-out", str(out_tsv)],
        check=True,
    )
    print(f"[run_blast] wrote {out_tsv}")
    return out_tsv


SARCOCYSTIDAE = {"GCF_000006565.2", "GCF_002563875.1", "GCF_000208865.1", "GCA_000727475.1"}
EIMERIIDAE_SISTER = {"GCF_002999335.1"}  # Cyclospora
DEEPER_APICOMPLEXA = {"GCF_000006425.1", "GCF_900005855.1"}

def classify_acc(acc):
    if acc in EIMERIA_ACCS:
        return ("eimeria", APICOMPLEXA_ACCS[acc])
    if acc in EIMERIIDAE_SISTER:
        return ("eimeriidae", APICOMPLEXA_ACCS[acc])
    if acc in SARCOCYSTIDAE:
        return ("sarcocystidae", APICOMPLEXA_ACCS[acc])
    if acc in DEEPER_APICOMPLEXA:
        return ("apicomplexa_deep", APICOMPLEXA_ACCS[acc])
    return ("unknown", acc)


def analyze_hits(blast_tsv):
    rows = []
    with open(blast_tsv) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12: continue
            qseqid = p[0]
            sseqid = p[1]
            pident = float(p[2])
            length = int(p[3])
            evalue = float(p[10])
            bitscore = float(p[11])
            # Strip subject acc
            sacc = sseqid.split("::", 1)[0] if "::" in sseqid else sseqid
            # Strip query acc
            qacc = qseqid.split("__", 1)[0]
            if pident < 80 or length < 30:
                continue
            # Skip self-hits (same accession)
            if sacc == qacc:
                continue
            rows.append({
                "qseqid": qseqid, "qacc": qacc,
                "sseqid": sseqid, "sacc": sacc,
                "pident": pident, "length": length,
                "evalue": evalue, "bitscore": bitscore,
            })

    print(f"\n=== Strong cross-genome hits (≥80% ident, ≥30 bp, non-self) ===")
    print(f"Total: {len(rows)}")

    # Per-linker
    per_q = defaultdict(list)
    for r in rows:
        per_q[r["qseqid"]].append(r)

    n_linkers_with_hit = len(per_q)
    print(f"Linkers with ≥1 cross-genome hit: {n_linkers_with_hit}")

    # Lineage tally
    lineage_count = Counter()
    for q, hits in per_q.items():
        lineages = set()
        for h in hits:
            lineage, name = classify_acc(h["sacc"])
            lineages.add(lineage)
        for l in lineages:
            lineage_count[l] += 1

    print(f"\nLinkers with hit, by source-lineage class:")
    for lin, n in lineage_count.most_common():
        print(f"  {lin:<15s}: {n}")

    # Top genome contributors
    sacc_count = Counter(h["sacc"] for h in rows)
    print(f"\nTop 20 source genomes (by total hits):")
    for sacc, n in sacc_count.most_common(20):
        lin, name = classify_acc(sacc)
        print(f"  {sacc:<25s}  ({lin:<12s} {name[:40]:<40s})  hits={n}")

    # Save full hit table (passing filters)
    out_hits = WORK / "eimeria_linker_strong_hits.tsv"
    if rows:
        cols = list(rows[0].keys())
        with open(out_hits, "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in rows:
                f.write("\t".join(str(r[c]) for c in cols) + "\n")
        print(f"\nWrote {out_hits}")

    summary = {
        "n_linkers_input": sum(1 for _ in open(WORK / "eimeria_pure_linkers.fa") if _.startswith(">")),
        "n_strong_cross_genome_hits": len(rows),
        "n_linkers_with_any_hit": n_linkers_with_hit,
        "lineage_count": dict(lineage_count),
        "top_source_genomes": sacc_count.most_common(20),
    }
    out_json = WORK / "eimeria_linker_blast_summary.json"
    json.dump(summary, open(out_json, "w"), indent=2, default=str)
    print(f"\nSummary: {out_json}")
    print(json.dumps(summary, indent=2, default=str))


def main():
    build_combined_fasta()
    build_blast_db()
    q = extract_linkers()
    out_tsv = run_blast(q)
    analyze_hits(out_tsv)


if __name__ == "__main__":
    main()
