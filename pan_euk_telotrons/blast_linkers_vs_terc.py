#!/usr/bin/env python3
"""
BLAST pure non-repeat telotron linkers against a curated TERC sequence database.

Open-question #1 from real_telotrons/ITS_MECHANISM_REPORT.md: are any of the
truly non-repeat medium-length linkers TERC retrotranscript fragments
(Nergadze 2007 model)?

Inputs:
  real_telotrons/strict_linkers_eimeria.json  (825 Eimeria linkers)
  real_telotrons/tara_linker_meta.json        (599 haptophyte linkers)

Steps:
  1. Filter linkers to "pure non-repeat" (no TTAGGG/CCCTAA kmers).
  2. Pull a curated TERC database via NCBI Entrez (Apicomplexa + diverse euks).
  3. BLAST linkers vs TERC db (blastn-short, e<1e-3, ident>=80%, len>=20).
  4. Report hits per linker, per source organism.
"""

import json
import subprocess
import sys
import time
from pathlib import Path
from collections import Counter

REAL_DIR = Path("/scratch1/alex/telotrons/pan_euk_telotrons/real_telotrons")
WORK_DIR = REAL_DIR / "terc_blast"
WORK_DIR.mkdir(parents=True, exist_ok=True)

FWD_TELO_KMERS = {'TTAGGG', 'TAGGGT', 'AGGGTT', 'GGGTTA', 'GGTTAG', 'GTTAGG'}
REV_TELO_KMERS = {'CCCTAA', 'CCTAAC', 'CTAACC', 'TAACCC', 'AACCCT', 'ACCCTA'}
# Also include TTTAGGG-family for plants/Apicomplexa (Eimeria uses TTTAGGG)
PLANT_FWD = {'TTTAGGG', 'TTAGGGT', 'TAGGGTT', 'AGGGTTT', 'GGGTTTA', 'GGTTTAG', 'GTTTAGG'}
PLANT_REV = {'CCCTAAA', 'CCTAAAC', 'CTAAACC', 'TAAACCC', 'AAACCCT', 'AACCCTA', 'ACCCTAA'}
ALL_TELO_KMERS = FWD_TELO_KMERS | REV_TELO_KMERS | PLANT_FWD | PLANT_REV


def has_telo(seq):
    return any(k in seq for k in ALL_TELO_KMERS)


def extract_pure_linkers():
    """Read both linker JSONs, return list of (id, seq) for pure non-repeat."""
    out = []

    # Eimeria
    eim = json.load(open(REAL_DIR / "strict_linkers_eimeria.json"))
    for r in eim:
        seq = (r.get('linker_seq') or '').upper()
        if len(seq) < 20:
            continue
        if has_telo(seq):
            continue
        rid = f"EIM__{r['acc']}__{r['contig']}__{r['start']}-{r['end']}__len{len(seq)}"
        out.append((rid, seq))

    # Tara haptophytes
    tara = json.load(open(REAL_DIR / "tara_linker_meta.json"))
    for mag, rows in tara.items():
        for r in rows:
            seq = (r.get('linker') or '').upper()
            if len(seq) < 20:
                continue
            if has_telo(seq):
                continue
            rid = f"TARA__{mag}__{r['contig']}__{r['start']}-{r['end']}__len{len(seq)}"
            out.append((rid, seq))

    return out


def fetch_terc_db():
    """Pull TERC nucleotide sequences from NCBI Entrez. Returns path to FASTA."""
    fa_path = WORK_DIR / "terc_db.fa"
    if fa_path.exists() and fa_path.stat().st_size > 50_000:
        print(f"[fetch_terc_db] reusing existing {fa_path} ({fa_path.stat().st_size} bytes)")
        return fa_path

    # Curated TERC search terms - one per major lineage. We deliberately favor
    # diverse representation over completeness; want a few hundred sequences.
    queries = [
        '("telomerase RNA"[Title] OR "telomerase RNA component"[Title] OR TERC[Title]) AND eukaryota[organism] AND 100:5000[Sequence Length]',
        '"telomerase template"[Title] AND eukaryota[organism] AND 100:5000[Sequence Length]',
        'TERC[gene] AND apicomplexa[organism] AND 100:5000[Sequence Length]',
        'TERC[gene] AND alveolata[organism] AND 100:5000[Sequence Length]',
        'TERC[gene] AND haptophyta[organism] AND 100:5000[Sequence Length]',
        'TERC[gene] AND stramenopiles[organism] AND 100:5000[Sequence Length]',
    ]

    all_ids = set()
    for q in queries:
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=nucleotide&term={q}&retmax=200&retmode=json"
        )
        try:
            r = subprocess.run(
                ["curl", "-s", "--max-time", "30",
                 "--data-urlencode", f"term={q}",
                 "--data", "db=nucleotide",
                 "--data", "retmax=200",
                 "--data", "retmode=json",
                 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"],
                check=True, capture_output=True, text=True, timeout=60,
            )
            d = json.loads(r.stdout)
            ids = d.get('esearchresult', {}).get('idlist', [])
            print(f"[fetch_terc_db] query «{q[:60]}…» -> {len(ids)} ids")
            all_ids.update(ids)
        except Exception as e:
            print(f"[fetch_terc_db] query failed: {e}", file=sys.stderr)
        time.sleep(0.5)  # NCBI rate limit

    if not all_ids:
        print("[fetch_terc_db] ERROR: no IDs found", file=sys.stderr)
        sys.exit(1)

    print(f"[fetch_terc_db] total unique TERC ids: {len(all_ids)}")

    # Fetch FASTA in batches
    ids_list = sorted(all_ids)
    BATCH = 100
    with open(fa_path, "w") as out:
        for i in range(0, len(ids_list), BATCH):
            batch = ids_list[i:i+BATCH]
            try:
                r = subprocess.run(
                    ["curl", "-s", "--max-time", "60",
                     "--data", f"db=nucleotide",
                     "--data", f"id={','.join(batch)}",
                     "--data", "rettype=fasta",
                     "--data", "retmode=text",
                     "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"],
                    check=True, capture_output=True, text=True, timeout=120,
                )
                out.write(r.stdout)
                print(f"[fetch_terc_db] fetched batch {i//BATCH + 1} ({len(batch)} ids)")
            except Exception as e:
                print(f"[fetch_terc_db] batch fail: {e}", file=sys.stderr)
            time.sleep(0.5)

    sz = fa_path.stat().st_size
    print(f"[fetch_terc_db] wrote {fa_path} ({sz} bytes)")
    if sz < 1000:
        print("[fetch_terc_db] WARNING: FASTA looks empty", file=sys.stderr)
    return fa_path


def run_blast(query_fa, subject_fa):
    """Run blastn-short of query vs subject. Return list of hit dicts."""
    db_prefix = WORK_DIR / "terc_db"
    if not (db_prefix.with_suffix(".nsq").exists() or
            (Path(str(db_prefix) + ".nsq")).exists()):
        subprocess.run(
            ["makeblastdb", "-in", str(subject_fa), "-dbtype", "nucl",
             "-out", str(db_prefix)],
            check=True, capture_output=True,
        )
        print(f"[run_blast] built db {db_prefix}")

    out_tsv = WORK_DIR / "linkers_vs_terc.blastn.tsv"
    fmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
    subprocess.run(
        ["blastn", "-task", "blastn-short",
         "-query", str(query_fa),
         "-db", str(db_prefix),
         "-evalue", "1e-3",
         "-outfmt", fmt,
         "-num_threads", "8",
         "-max_target_seqs", "5",
         "-out", str(out_tsv)],
        check=True,
    )
    print(f"[run_blast] wrote {out_tsv}")

    hits = []
    with open(out_tsv) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 13:
                continue
            hits.append({
                'qseqid': p[0], 'sseqid': p[1],
                'pident': float(p[2]), 'length': int(p[3]),
                'mismatch': int(p[4]), 'gapopen': int(p[5]),
                'qstart': int(p[6]), 'qend': int(p[7]),
                'sstart': int(p[8]), 'send': int(p[9]),
                'evalue': float(p[10]), 'bitscore': float(p[11]),
                'stitle': p[12],
            })
    return hits


def main():
    # 1. Pure non-repeat linkers
    pure = extract_pure_linkers()
    print(f"\n=== Pure non-repeat linkers ===")
    print(f"Total: {len(pure)}")
    by_source = Counter(rid.split("__")[0] for rid, _ in pure)
    for src, n in by_source.items():
        print(f"  {src}: {n}")

    # Length distribution
    lens = [len(s) for _, s in pure]
    print(f"  length: min={min(lens)}, median={sorted(lens)[len(lens)//2]}, max={max(lens)}")
    medium = sum(1 for x in lens if 30 <= x <= 120)
    print(f"  medium (30-120bp, 'TERC-frag zone'): {medium}")

    query_fa = WORK_DIR / "pure_nonrepeat_linkers.fa"
    with open(query_fa, "w") as f:
        for rid, seq in pure:
            f.write(f">{rid}\n{seq}\n")
    print(f"  -> {query_fa}")

    # 2. Build TERC database
    terc_fa = fetch_terc_db()
    n_terc = sum(1 for line in open(terc_fa) if line.startswith(">"))
    print(f"\n=== TERC database ===")
    print(f"Sequences: {n_terc}")

    # 3. BLAST
    print(f"\n=== BLAST linkers vs TERC ===")
    hits = run_blast(query_fa, terc_fa)
    print(f"Total raw hits (e<1e-3): {len(hits)}")

    # 4. Filter: ident >= 80, length >= 20
    strong = [h for h in hits if h['pident'] >= 80 and h['length'] >= 20]
    print(f"Strong hits (pident>=80, len>=20): {len(strong)}")

    if strong:
        print("\n--- Strong hits ---")
        # Per linker
        per_q = {}
        for h in strong:
            per_q.setdefault(h['qseqid'], []).append(h)
        for q, hh in sorted(per_q.items())[:30]:
            top = max(hh, key=lambda x: x['bitscore'])
            print(f"  {q}")
            print(f"    -> {top['stitle'][:90]}")
            print(f"       pident={top['pident']}, len={top['length']}, "
                  f"qry={top['qstart']}-{top['qend']}, sbj={top['sstart']}-{top['send']}, "
                  f"e={top['evalue']:.2e}")

    # 5. Save summary
    summary = {
        'n_linkers_total': len(pure),
        'n_linkers_eimeria': by_source.get('EIM', 0),
        'n_linkers_tara': by_source.get('TARA', 0),
        'n_terc_seqs': n_terc,
        'n_raw_hits': len(hits),
        'n_strong_hits': len(strong),
        'n_unique_linkers_with_strong_hit': len({h['qseqid'] for h in strong}),
        'top_subjects': Counter(h['sseqid'] for h in strong).most_common(10),
    }
    out_json = WORK_DIR / "terc_blast_summary.json"
    json.dump(summary, open(out_json, "w"), indent=2)
    print(f"\nSummary: {out_json}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
