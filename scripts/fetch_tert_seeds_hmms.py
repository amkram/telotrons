#!/usr/bin/env python3
"""Fetch TERT search references: apicomplexan TERT seed proteins (NCBI E-utils)
and the two diagnostic Pfam HMMs (Telomerase_RBD PF12009, RVT_1 PF00078).

Cached: if the outputs already exist and are non-empty, the network is not hit.
Used by the find_tert rule (deep-homology TERT identification).
"""
import argparse
import gzip
import os
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
INTERPRO = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam"
PFAMS = {"PF12009": "Telomerase_RBD", "PF00078": "RVT_1"}


def _get(url, binary=False, tries=4):
    last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                data = r.read()
            return data if binary else data.decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"failed GET {url}: {last}")


def fetch_seeds(out_faa):
    if os.path.exists(out_faa) and os.path.getsize(out_faa) > 0:
        print(f"seeds cached: {out_faa}")
        return
    term = ('(telomerase reverse transcriptase[Protein Name]) '
            'AND Apicomplexa[Organism] NOT partial[All Fields]')
    q = urllib.parse.quote(term)
    ids_xml = _get(f"{EUTILS}/esearch.fcgi?db=protein&term={q}&retmax=60")
    ids = []
    for chunk in ids_xml.split("<Id>")[1:]:
        ids.append(chunk.split("</Id>")[0])
    if not ids:
        sys.exit("no apicomplexan TERT seed IDs returned by esearch")
    print(f"esearch: {len(ids)} TERT protein IDs")
    fasta = _get(f"{EUTILS}/efetch.fcgi?db=protein&id={','.join(ids)}"
                 f"&rettype=fasta&retmode=text")
    with open(out_faa, "w") as fh:
        fh.write(fasta)
    n = fasta.count(">")
    print(f"wrote {out_faa} ({n} sequences)")


def fetch_hmm(pf, out_hmm):
    if os.path.exists(out_hmm) and os.path.getsize(out_hmm) > 0:
        print(f"{pf} HMM cached: {out_hmm}")
        return
    gz = _get(f"{INTERPRO}/{pf}?annotation=hmm", binary=True)
    data = gzip.decompress(gz)
    with open(out_hmm, "wb") as fh:
        fh.write(data)
    print(f"wrote {out_hmm} ({PFAMS[pf]})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    fetch_seeds(os.path.join(args.outdir, "tert_seeds.faa"))
    for pf in PFAMS:
        fetch_hmm(pf, os.path.join(args.outdir, f"{pf}.hmm"))


if __name__ == "__main__":
    main()
