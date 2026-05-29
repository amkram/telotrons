#!/usr/bin/env python3
"""Generalized deep-homology gene finder (intron-aware miniprot + Pfam HMM confirm).

Adapted from the validated find_tert_deep_homology.py. Detects a divergent,
possibly-unannotated gene across a genome set and makes a present/absent/
pseudogenized call per genome. Designed for the telomere-restriction-factor
survey (Pif1, CST, shelterin, MRN, ...), where the key question is whether a
factor present in telotron-poor sisters is specifically lost in Eimeria.

Confirmation: a candidate is "confirmed" if it carries the diagnostic Pfam
domain(s) at i-Evalue < --hmm-evalue. Iterative within-clade seeding (round>1
adds confirmed proteins back as seeds) rescues divergent genomes that a distant
seed only partially maps -- this is what overturned the TERT "absent" call.

Per-genome call:
    PRESENT     : confirmed locus (domain below threshold), low internal stops
    PSEUDOGENE? : locus maps + domain near threshold but many internal stops/frameshifts
    ABSENT      : no domain-bearing locus found (only meaningful if sister genomes are PRESENT)
"""
import argparse, glob, gzip, os, shutil, subprocess as sp, sys, tempfile
import pandas as pd

ENV_BIN = "/home/alex/.snakemake-envs/telotrons/7aa0f2645427a08bdf054e56b1233123_/bin"


def which(name, fb):
    return shutil.which(name) or next((f for f in fb if os.path.exists(f)), None)


def find_genome_fasta(gid, refseq_dir, tara_dir):
    if gid.startswith("TARA"):
        hits = glob.glob(f"{tara_dir}/Contigs/{gid}.fa*")
    else:
        hits = glob.glob(f"{refseq_dir}/{gid}/**/*_genomic.fna*", recursive=True)
    hits = [h for h in hits if not h.endswith((".fai", ".gzi", ".bai"))]
    plain = [h for h in hits if h.endswith((".fna", ".fa", ".fasta"))]
    pick = plain or [h for h in hits if h.endswith(".gz")] or hits
    if not pick:
        raise FileNotFoundError(gid)
    return pick[0]


def ensure_plain(src, workdir, gid):
    if src.endswith(".gz"):
        dst = os.path.join(workdir, f"{gid}.fna")
        if not os.path.exists(dst) or os.path.getsize(dst) == 0:
            with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
                shutil.copyfileobj(fi, fo)
        return dst
    return src


def run_miniprot_trans(miniprot, genome, seeds, threads, outc, outn):
    proc = sp.run([miniprot, "-t", str(threads), "--trans",
                   f"--outc={outc}", f"--outn={outn}", genome, seeds],
                  check=True, stdout=sp.PIPE, stderr=sp.DEVNULL, text=True)
    cands, paf = [], None
    for line in proc.stdout.splitlines():
        if line.startswith("##STA"):
            if paf is None:
                continue
            prot = line.split("\t", 1)[1].strip() if "\t" in line else ""
            f = paf
            cands.append({"seed": f[0], "qlen": int(f[1]), "qstart": int(f[2]),
                          "qend": int(f[3]), "strand": f[4], "contig": f[5],
                          "start": min(int(f[7]), int(f[8])), "end": max(int(f[7]), int(f[8])),
                          "protein": prot.replace("*", "X")})
            paf = None
        elif not line.startswith("##"):
            paf = line.split("\t")
    return cands


def hmm_best(hmmsearch, hmm, faa, td, tag):
    if os.path.getsize(faa) == 0:
        return {}
    dt = os.path.join(td, f"{tag}.domtbl")
    sp.run([hmmsearch, "--max", "-E", "10", "--domtblout", dt, hmm, faa],
           check=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
    best = {}
    with open(dt) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) >= 19:
                nm, e = f[0], float(f[12])
                if nm not in best or e < best[nm]:
                    best[nm] = e
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True, help="gene label, e.g. Pif1")
    ap.add_argument("--genome-ids", required=True)
    ap.add_argument("--seeds", required=True)
    ap.add_argument("--hmm", required=True, nargs="+", help="diagnostic Pfam HMM(s)")
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--hmm-evalue", type=float, default=1e-5)
    ap.add_argument("--iterate", type=int, default=3)
    ap.add_argument("--outc", type=float, default=0.05)
    ap.add_argument("--outn", type=int, default=8)
    ap.add_argument("--threads", type=int, default=16)
    args = ap.parse_args()

    miniprot = which("miniprot", [f"{ENV_BIN}/miniprot"])
    hmmsearch = which("hmmsearch", ["/usr/bin/hmmsearch"])
    gids = [g.strip() for g in args.genome_ids.split(",") if g.strip()]
    os.makedirs(args.outdir, exist_ok=True)
    prot_dir = os.path.join(args.outdir, "proteins")
    os.makedirs(prot_dir, exist_ok=True)
    work = os.path.join(args.outdir, "_work")
    os.makedirs(work, exist_ok=True)

    genome_fa = {}
    for gid in gids:
        try:
            genome_fa[gid] = ensure_plain(find_genome_fasta(gid, args.refseq_dir, args.tara_dir), work, gid)
        except FileNotFoundError:
            print(f"skip {gid}: genome not found", file=sys.stderr)

    seeds = os.path.join(work, "seeds.faa")
    shutil.copyfile(args.seeds, seeds)
    seen = set()
    best_per_genome = {}

    for rnd in range(1, args.iterate + 1):
        nseed = sum(1 for _ in open(seeds) if _.startswith(">"))
        print(f"\n=== {args.gene} ROUND {rnd} (seeds {nseed}) ===")
        round_conf = []
        for gid, fa in genome_fa.items():
            cands = run_miniprot_trans(miniprot, fa, seeds, args.threads, args.outc, args.outn)
            if not cands:
                print(f"  {gid}: 0 candidates"); continue
            with tempfile.TemporaryDirectory(dir=work) as td:
                faa = os.path.join(td, "c.faa")
                names = [f"{gid}__r{rnd}c{i}" for i in range(len(cands))]
                with open(faa, "w") as o:
                    for nm, c in zip(names, cands):
                        o.write(f">{nm}\n{c['protein']}\n")
                evals = {}
                for hp in args.hmm:
                    for nm, e in hmm_best(hmmsearch, hp, faa, td, os.path.basename(hp)).items():
                        evals[nm] = min(evals.get(nm, 9e9), e)
            nconf = 0
            for nm, c in zip(names, cands):
                e = evals.get(nm, float("nan"))
                stops = c["protein"].count("X")
                conf = pd.notna(e) and e < args.hmm_evalue
                prev = best_per_genome.get(gid)
                cur = {"genome_id": gid, "contig": c["contig"], "start": c["start"],
                       "end": c["end"], "strand": c["strand"], "prot_len": len(c["protein"]),
                       "internal_X": stops, "hmm_evalue": e, "seed": c["seed"],
                       "round": rnd, "confirmed": conf, "protein": c["protein"]}
                # keep best by (confirmed, lowest evalue, longest)
                key = lambda r: (r["confirmed"], -(r["hmm_evalue"] if pd.notna(r["hmm_evalue"]) else 9e9), r["prot_len"])
                if prev is None or key(cur) > key(prev):
                    best_per_genome[gid] = cur
                if conf:
                    nconf += 1
                    round_conf.append((nm, c["protein"]))
            be = best_per_genome[gid]["hmm_evalue"]
            print(f"  {gid}: {len(cands)} cand, {nconf} confirmed; best HMM E={be:.1e}" if pd.notna(be)
                  else f"  {gid}: {len(cands)} cand, 0 confirmed")
        if rnd < args.iterate and round_conf:
            with open(seeds, "a") as o:
                for nm, p in round_conf:
                    if p not in seen:
                        seen.add(p); o.write(f">{nm}\n{p}\n")

    rows = []
    for gid in gids:
        b = best_per_genome.get(gid)
        if b is None:
            rows.append({"genome_id": gid, "call": "ABSENT", "hmm_evalue": None,
                         "prot_len": 0, "internal_X": 0, "contig": "", "start": 0, "end": 0})
            continue
        if b["confirmed"] and b["internal_X"] <= 5:
            call = "PRESENT"
        elif b["confirmed"]:
            call = "PRESENT_frameshifted"
        elif pd.notna(b["hmm_evalue"]) and b["hmm_evalue"] < 1e-2:
            call = "PSEUDOGENE?"
        else:
            call = "ABSENT"
        rows.append({"genome_id": gid, "call": call, "hmm_evalue": b["hmm_evalue"],
                     "prot_len": b["prot_len"], "internal_X": b["internal_X"],
                     "contig": b["contig"], "start": b["start"], "end": b["end"]})
        open(os.path.join(prot_dir, f"{gid}.faa"), "w").write(
            f">{gid}|{b['contig']}:{b['start']}-{b['end']}|{b['strand']}|E={b['hmm_evalue']}\n{b['protein']}\n")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.outdir, f"{args.gene}_calls.tsv"), sep="\t", index=False)
    print(f"\n=== {args.gene} present/absent ===")
    for r in rows:
        e = f"{r['hmm_evalue']:.1e}" if r['hmm_evalue'] is not None and pd.notna(r['hmm_evalue']) else "NA"
        print(f"  {r['genome_id']:18s} {r['call']:20s} HMM_E={e:>9} len={r['prot_len']}")
    print(f"\nwrote {args.outdir}/{args.gene}_calls.tsv")


if __name__ == "__main__":
    main()
