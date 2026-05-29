#!/usr/bin/env python3
"""For each target genome, build a per-genome proteome from its closest RefSeq relatives.

For each genome_id in --manifest we:
  1. Look up its row in the (cached) full RefSeq assembly_summary_refseq.txt
  2. Find sibling RefSeq assemblies of the same genus (organism_name first token),
     excluding self and excluding assemblies with no annotation. Prefer
     'reference genome' / 'representative genome' / 'Complete Genome' first;
     fall back to any annotated assembly. Cap at --max-siblings.
  3. Download each sibling's `_protein.faa.gz` from its FTP path.
  4. Concatenate the downloaded fastas into `work/manifests/close_proteomes/{gid}.faa.gz`.

For TARA MAGs (no taxonomy here) we write an empty file so downstream tools
don't have to special-case missing entries — they'll just fall back to the
base proteome.

Cached files:
  work/manifests/assembly_summary_refseq.txt           (full RefSeq assembly summary)
  work/manifests/close_proteomes/_cache/{accession}_protein.faa.gz   (per-sibling)
  work/manifests/close_proteomes/{gid}.faa.gz                        (per-target union)

The per-sibling cache means multiple targets that share a sibling
(e.g. five Eimeria spp. all pulling Toxoplasma) only download it once.
"""
import argparse
import os
import sys
import subprocess as sp

import pandas as pd

REFSEQ_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"


def curl_to(path, url, timeout=600):
    tmp = path + ".tmp"
    rc = sp.call(["curl", "-sSL", "--fail", "--retry", "5", "--retry-delay", "3",
                  "--retry-all-errors", "--max-time", str(timeout),
                  "-o", tmp, url])
    if rc == 0 and os.path.exists(tmp) and os.path.getsize(tmp) > 0:
        os.rename(tmp, path)
        return True
    if os.path.exists(tmp):
        os.remove(tmp)
    return False


def ensure_summary(cache_path):
    if os.path.exists(cache_path) and os.path.getsize(cache_path) > 100_000_000:
        return cache_path
    print(f"  fetching {REFSEQ_SUMMARY_URL} → {cache_path}", flush=True)
    if not curl_to(cache_path, REFSEQ_SUMMARY_URL, timeout=600):
        raise RuntimeError(f"curl failed to download {REFSEQ_SUMMARY_URL}")
    return cache_path


def load_summary(path):
    """Return dict accession→row (selected columns). Header line starts with '#assembly_accession'."""
    cols = None
    rows = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                cols = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if not cols:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < len(cols):
                f += [""] * (len(cols) - len(f))
            r = dict(zip(cols, f))
            acc = r.get("assembly_accession", "")
            if acc:
                rows[acc] = r
    return rows


def _na(s):
    s = (s or "").strip()
    return "" if s.lower() == "na" else s


def rank_sibling(r):
    """Lower is better."""
    cat = (r.get("refseq_category") or "").lower()
    lvl = (r.get("assembly_level") or "").lower()
    excl = _na(r.get("excluded_from_refseq"))
    annot = _na(r.get("annotation_release")) or _na(r.get("annotation_date"))
    score = 0
    if excl:
        score += 1000
    if not annot:
        score += 500  # no annotation → no protein file
    if cat == "reference genome":
        score -= 50
    elif cat == "representative genome":
        score -= 25
    if lvl in ("complete genome", "chromosome"):
        score -= 10
    elif lvl == "scaffold":
        score += 5
    elif lvl == "contig":
        score += 10
    return score


def genus_of(name):
    if not name:
        return ""
    return name.split()[0]


def find_siblings(target_row, all_rows, max_n):
    """Pick up to max_n best annotated assemblies in the same genus, excluding self."""
    self_acc = target_row["assembly_accession"]
    g = genus_of(target_row.get("organism_name", ""))
    if not g:
        return []
    cands = []
    for acc, r in all_rows.items():
        if acc == self_acc:
            continue
        if genus_of(r.get("organism_name", "")) != g:
            continue
        if _na(r.get("excluded_from_refseq")):
            continue
        if not (r.get("ftp_path") or "").startswith("http"):
            continue
        cands.append(r)
    cands.sort(key=rank_sibling)
    return cands[:max_n]


def download_protein(ftp_path, acc, cache_dir):
    out = os.path.join(cache_dir, f"{acc}_protein.faa.gz")
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return out
    base = ftp_path.rstrip("/").rsplit("/", 1)[-1]
    url = f"{ftp_path.rstrip('/')}/{base}_protein.faa.gz"
    if curl_to(out, url, timeout=180):
        return out
    print(f"    miss {acc}: curl failed for {url}", flush=True)
    return None


def concat_gz(paths, out_path):
    """Concatenate multiple .gz files into one .gz (gz is concat-safe)."""
    if not paths:
        # write a tiny placeholder so consumers can still test for file existence
        open(out_path, "wb").close()
        return 0
    total = 0
    with open(out_path, "wb") as fout:
        for p in paths:
            with open(p, "rb") as fin:
                while True:
                    chunk = fin.read(1 << 20)
                    if not chunk:
                        break
                    fout.write(chunk)
            total += os.path.getsize(p)
    return total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, help="work/results/all_species_raw_summary.tsv or any TSV with genome_id")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--summary-cache",
                    default="work/manifests/assembly_summary_refseq.txt")
    ap.add_argument("--max-siblings", type=int, default=5)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    cache_dir = os.path.join(args.outdir, "_cache")
    os.makedirs(cache_dir, exist_ok=True)

    summary_path = ensure_summary(args.summary_cache)
    print("  loading assembly_summary…", flush=True)
    all_rows = load_summary(summary_path)
    print(f"  {len(all_rows)} assemblies indexed", flush=True)

    manifest = pd.read_csv(args.manifest, sep="\t")
    for _, m in manifest.iterrows():
        gid = m.genome_id
        out_path = os.path.join(args.outdir, f"{gid}.faa.gz")
        if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
            print(f"  {gid}: already built ({os.path.getsize(out_path):,} bytes)", flush=True)
            continue

        if gid.startswith("TARA") or gid not in all_rows:
            print(f"  {gid}: no assembly_summary entry — empty placeholder", flush=True)
            concat_gz([], out_path)
            continue

        self_row = all_rows[gid]
        siblings = find_siblings(self_row, all_rows, args.max_siblings)
        if not siblings:
            print(f"  {gid}: no same-genus siblings ({genus_of(self_row.get('organism_name',''))})",
                  flush=True)
            concat_gz([], out_path)
            continue

        print(f"  {gid} ({self_row.get('organism_name','?')}): "
              f"{len(siblings)} siblings →", flush=True)
        downloaded = []
        for r in siblings:
            acc = r["assembly_accession"]
            name = r.get("organism_name", "?")
            p = download_protein(r["ftp_path"], acc, cache_dir)
            if p:
                downloaded.append(p)
                print(f"    + {acc} {name}", flush=True)
        n = concat_gz(downloaded, out_path)
        print(f"    → {out_path} ({n:,} bytes uncompressed-input)", flush=True)


if __name__ == "__main__":
    sys.exit(main())
