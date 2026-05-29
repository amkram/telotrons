#!/usr/bin/env python3
"""Build telomerase reference FASTAs from telomerase.us curated sequence links.

Pages:
  TR     https://telomerase.us/sequences_tr.html    — telomerase RNA component
  TERT   https://telomerase.us/sequences_tert.html  — reverse transcriptase
  Other  https://telomerase.us/sequences_other.html — accessory proteins

The curated sequences are linked as NCBI sequence-viewer URLs carrying a GI
number and a db:

  .../viewer.fcgi?db=nucleotide&val=<GI>
  .../viewer.fcgi?db=nuccore&id=<GI>
  .../viewer.fcgi?db=protein&id=<GI>

We follow ONLY those links — earlier text-scraping of accession-like strings
pulled in citation/genome accessions, and `fasta_cds_aa` on a whole-genome
record exploded into thousands of unrelated CDS proteins (96% contamination).

Per page we emit:
  TR     → {tr}.nucl.fa     (nucleotide GIs; RNA component)
  TERT   → {tert}.nucl.fa   (nucleotide GIs; mRNA/gene records)
           {tert}.prot.fa   (protein GIs direct + CDS translations of nucl GIs)
  Other  → {other}.nucl.fa  (nucleotide GIs)
           {other}.prot.fa  (protein GIs direct + CDS translations of nucl GIs)

Two guards keep genome records out:
  * MAX_NT_LEN — skip any nucleotide record longer than this (telomerase RNA
    < 2 kb, TERT mRNA < 5 kb; the only oversized links are whole genomes, e.g.
    the 178 kb Gallid-herpesvirus record that carries the viral TR).
  * MAX_CDS_PER_RECORD — when translating a nucleotide record's CDS, drop the
    whole record if it yields more than this many proteins (a real telomerase
    mRNA yields 1; a genome yields thousands).
"""
import argparse
import os
import re
import sys
import time
import urllib.parse
import urllib.request

PAGES = {
    "tr":    "https://telomerase.us/sequences_tr.html",
    "tert":  "https://telomerase.us/sequences_tert.html",
    "other": "https://telomerase.us/sequences_other.html",
}

# What to emit per category.
PLAN = {
    "tr":    {"nucl": True,  "prot": False},
    "tert":  {"nucl": True,  "prot": True},
    "other": {"nucl": True,  "prot": True},
}

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH = 200
SLEEP = 0.2

MAX_NT_LEN = 20000           # skip nucleotide records longer than this (genomes)
MAX_CDS_PER_RECORD = 3       # drop a nucl record if its CDS translation exceeds this

# Sequence-viewer links carry the GI in val= or id=, with db in {nucleotide,
# nuccore, protein}. nucleotide/nuccore are the same molecule class for efetch.
VIEWER_RE = re.compile(
    r"viewer\.fcgi\?db=(nucleotide|nuccore|protein)"
    r"(?:&amp;|&)(?:val|id)=(\d+)",
    re.IGNORECASE,
)


def http_get(url, ua="telotron-pipeline/1.0", retries=3, sleep=SLEEP):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": ua})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception as e:
            last = e
            time.sleep(sleep * (2 ** i))
    raise RuntimeError(f"GET failed for {url}: {last}")


def cache_html(name, url, cache_dir):
    path = os.path.join(cache_dir, f"{name}.html")
    if not os.path.exists(path):
        data = http_get(url)
        with open(path, "wb") as fh:
            fh.write(data)
        time.sleep(SLEEP)
    return path


def extract_gis(html_text):
    """Return (nucl_gis, prot_gis): ordered-unique GI strings by molecule class."""
    nucl, prot = [], []
    seen = set()
    for m in VIEWER_RE.finditer(html_text):
        db, gi = m.group(1).lower(), m.group(2)
        if gi in seen:
            continue
        seen.add(gi)
        (prot if db == "protein" else nucl).append(gi)
    return nucl, prot


def efetch(ids, db, rettype):
    params = {"db": db, "id": ",".join(ids), "rettype": rettype, "retmode": "text"}
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        params["api_key"] = api_key
    url = EUTILS + "?" + urllib.parse.urlencode(params)
    return http_get(url).decode("utf-8", errors="replace")


def parse_fasta(txt):
    """Yield (header_without_'>', sequence) records from FASTA text."""
    header, buf = None, []
    for line in txt.splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(buf)
            header, buf = line[1:], []
        elif header is not None:
            buf.append(line.strip())
    if header is not None:
        yield header, "".join(buf)


def fetch_nucl(gis, out_path, log_prefix=""):
    """Plain nucleotide FASTA, skipping records longer than MAX_NT_LEN."""
    n_kept = n_skip = 0
    with open(out_path, "w") as out:
        for i in range(0, len(gis), BATCH):
            batch = gis[i:i + BATCH]
            try:
                txt = efetch(batch, "nuccore", "fasta")
            except Exception as e:
                print(f"{log_prefix}  batch {i} failed: {e}", flush=True)
                time.sleep(SLEEP * 5)
                continue
            for header, seq in parse_fasta(txt):
                if len(seq) > MAX_NT_LEN:
                    n_skip += 1
                    print(f"{log_prefix}  skip {len(seq)} bp (>{MAX_NT_LEN}): "
                          f"{header[:70]}", flush=True)
                    continue
                out.write(f">{header}\n{seq}\n")
                n_kept += 1
            time.sleep(SLEEP)
    return n_kept, n_skip


def fetch_prot_direct(gis, out_fh, log_prefix=""):
    """Protein FASTA fetched directly by protein GI."""
    n = 0
    for i in range(0, len(gis), BATCH):
        batch = gis[i:i + BATCH]
        try:
            txt = efetch(batch, "protein", "fasta")
        except Exception as e:
            print(f"{log_prefix}  batch {i} failed: {e}", flush=True)
            time.sleep(SLEEP * 5)
            continue
        for header, seq in parse_fasta(txt):
            if seq:
                out_fh.write(f">{header}\n{seq}\n")
                n += 1
        time.sleep(SLEEP)
    return n


def fetch_prot_from_nucl(gis, out_fh, log_prefix=""):
    """CDS translations of nucleotide records. Fetched one GI at a time so a
    genome that slips past the link filter can be rejected by CDS count without
    discarding the rest of the batch."""
    n_kept = n_skip_records = 0
    for gi in gis:
        try:
            txt = efetch([gi], "nuccore", "fasta_cds_aa")
        except Exception as e:
            print(f"{log_prefix}  GI {gi} failed: {e}", flush=True)
            time.sleep(SLEEP * 5)
            continue
        recs = [(h, s) for h, s in parse_fasta(txt) if s]
        if len(recs) > MAX_CDS_PER_RECORD:
            n_skip_records += 1
            print(f"{log_prefix}  skip GI {gi}: {len(recs)} CDS "
                  f"(>{MAX_CDS_PER_RECORD}) — looks like a genome", flush=True)
            time.sleep(SLEEP)
            continue
        for header, seq in recs:
            out_fh.write(f">{header}\n{seq}\n")
            n_kept += 1
        time.sleep(SLEEP)
    return n_kept, n_skip_records


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.cache_dir, exist_ok=True)
    os.makedirs(args.outdir, exist_ok=True)

    summary = []
    for name, url in PAGES.items():
        print(f"[{name}] {url}", flush=True)
        html_path = cache_html(name, url, args.cache_dir)
        html_text = open(html_path).read()
        nucl_gis, prot_gis = extract_gis(html_text)
        with open(os.path.join(args.outdir, f"{name}.gis.txt"), "w") as fh:
            fh.write("# nucl\n" + "\n".join(nucl_gis) + "\n")
            fh.write("# prot\n" + "\n".join(prot_gis) + "\n")
        print(f"  GIs: nucl={len(nucl_gis)} prot={len(prot_gis)}", flush=True)

        n_nucl = n_prot = 0
        if PLAN[name]["nucl"] and nucl_gis:
            out = os.path.join(args.outdir, f"{name}.nucl.fa")
            n_nucl, n_skip = fetch_nucl(nucl_gis, out, log_prefix=f"[{name}.nucl]")
            print(f"  → {out}  ({n_nucl} kept, {n_skip} oversized skipped)",
                  flush=True)
        if PLAN[name]["prot"]:
            out = os.path.join(args.outdir, f"{name}.prot.fa")
            with open(out, "w") as oh:
                if prot_gis:
                    n_prot += fetch_prot_direct(prot_gis, oh,
                                                log_prefix=f"[{name}.prot.direct]")
                if nucl_gis:
                    n, n_skip = fetch_prot_from_nucl(nucl_gis, oh,
                                                     log_prefix=f"[{name}.prot.cds]")
                    n_prot += n
                    print(f"  prot.cds: {n} kept, {n_skip} genome-like skipped",
                          flush=True)
            print(f"  → {out}  ({n_prot} records)", flush=True)
        summary.append((name, len(nucl_gis), len(prot_gis), n_nucl, n_prot))

    print("\nSummary  page  nucl_gis prot_gis  nucl_recs prot_recs")
    for name, nn_gi, np_gi, nn, np_ in summary:
        print(f"  {name:6s} {nn_gi:6d} {np_gi:6d}   {nn:6d}  {np_:6d}")


if __name__ == "__main__":
    sys.exit(main())
