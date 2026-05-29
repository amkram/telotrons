#!/usr/bin/env python3
"""Identify telomerase reverse transcriptase (TERT) by deep homology.

BLAST/tblastn fail on apicomplexan TERT: the gene is ~25% identical to known
TERTs, heavily multi-exon, and unannotated across all of Coccidia. This script
recovers it with the field-standard approach for the twilight zone:

    1. intron-aware protein->genome mapping  (miniprot --trans)
       using phylogenetically-local apicomplexan TERT seeds
    2. domain-architecture confirmation       (hmmsearch vs Pfam)
       a real TERT carries the telomerase RBD (PF12009) UPSTREAM of the
       reverse-transcriptase domain (PF00078). Retroelement RTs (abundant in
       Eimeria) carry the RT but NEVER the RBD -- that is the discriminator.
    3. optional iterative recovery: confirmed proteins are added to the seed
       set and the search re-run, so a clean within-clade ortholog rescues the
       fragmented/divergent genomes a distant seed only partially mapped.

Outputs (under --outdir):
    confirmed_tert.tsv         one row per confirmed locus (all genomes)
    tert_proteins/<gid>.faa    best confirmed protein per genome
    all_confirmed.faa          every confirmed protein (for tree building)
    gff/<gid>.miniprot.gff     full miniprot GFF per genome (provenance)
"""
import argparse
import gzip
import os
import shutil
import subprocess as sp
import sys
import tempfile

import pandas as pd

from _common import find_genome_fasta, ensure_tool_on_path

# Make required binaries reachable for bare `python scripts/...` runs outside
# `snakemake --use-conda` (no-op when already on PATH).
ensure_tool_on_path("miniprot", "hmmsearch")


# ── tool resolution ─────────────────────────────────────────────────────────
def _resolve(name, override, fallbacks):
    if override:
        return override
    p = shutil.which(name)
    if p:
        return p
    for fb in fallbacks:
        if os.path.exists(fb):
            return fb
    sys.exit(f"ERROR: cannot find '{name}' on PATH or fallbacks {fallbacks}")



def ensure_plain_fasta(src, workdir, gid):
    """Decompress .gz to workdir; return a plain-fasta path miniprot can read."""
    if src.endswith(".gz"):
        dst = os.path.join(workdir, f"{gid}.fna")
        if not os.path.exists(dst) or os.path.getsize(dst) == 0:
            with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
                shutil.copyfileobj(fin, fout)
        return dst
    return src


# ── miniprot ─────────────────────────────────────────────────────────────────
def run_miniprot_trans(miniprot, genome_fa, seeds_fa, threads, outc, outn):
    """Return list of dicts: {contig,start,end,strand,seed,qlen,qstart,qend,nres,protein}."""
    cmd = [miniprot, "-t", str(threads), "--trans",
           f"--outc={outc}", f"--outn={outn}", genome_fa, seeds_fa]
    proc = sp.run(cmd, check=True, stdout=sp.PIPE, stderr=sp.DEVNULL, text=True)
    cands, paf = [], None
    for line in proc.stdout.splitlines():
        if line.startswith("##STA"):
            if paf is None:
                continue
            protein = line.split("\t", 1)[1].strip() if "\t" in line else ""
            f = paf
            cands.append({
                "seed": f[0], "qlen": int(f[1]), "qstart": int(f[2]), "qend": int(f[3]),
                "strand": f[4], "contig": f[5],
                "start": min(int(f[7]), int(f[8])), "end": max(int(f[7]), int(f[8])),
                "nres": int(f[9]), "protein": protein.replace("*", ""),
            })
            paf = None
        elif not line.startswith("##"):
            paf = line.split("\t")
    return cands


def run_miniprot_gff(miniprot, genome_fa, seeds_fa, threads, out_gff, outc, outn):
    proc = sp.run([miniprot, "-t", str(threads), "--gff",
                   f"--outc={outc}", f"--outn={outn}", genome_fa, seeds_fa],
                  check=True, stdout=sp.PIPE, stderr=sp.DEVNULL, text=True)
    with open(out_gff, "w") as fh:
        fh.write(proc.stdout)


# ── hmmsearch ─────────────────────────────────────────────────────────────────
def hmmsearch_best(hmmsearch, hmm, faa, tmpdir, tag):
    """Return {target_name: (best_dom_evalue, ali_from, ali_to)} over all domains."""
    if os.path.getsize(faa) == 0:
        return {}
    domtbl = os.path.join(tmpdir, f"{tag}.domtbl")
    sp.run([hmmsearch, "--max", "-E", "10", "--domtblout", domtbl, hmm, faa],
           check=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
    best = {}
    with open(domtbl) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 21:
                continue
            name = f[0]
            dom_e = float(f[12])          # i-Evalue (independent)
            ali_from, ali_to = int(f[17]), int(f[18])   # env coords on target
            if name not in best or dom_e < best[name][0]:
                best[name] = (dom_e, ali_from, ali_to)
    return best


def write_faa(records, path):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome-ids", required=True,
                    help="comma-separated genome_ids to search")
    ap.add_argument("--seeds", required=True, help="apicomplexan TERT seed proteins (FASTA)")
    ap.add_argument("--trbd-hmm", required=True, help="Pfam PF12009 Telomerase_RBD HMM")
    ap.add_argument("--rt-hmm", required=True, help="Pfam PF00078 RVT_1 HMM")
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--iterate", type=int, default=2,
                    help="rounds; >1 augments seeds with confirmed proteins (within-clade rescue)")
    ap.add_argument("--trbd-evalue", type=float, default=1e-5)
    ap.add_argument("--rt-evalue", type=float, default=1e-2)
    ap.add_argument("--outc", type=float, default=0.05)
    ap.add_argument("--outn", type=int, default=5)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--miniprot", default=None)
    ap.add_argument("--hmmsearch", default=None)
    args = ap.parse_args()

    miniprot = _resolve("miniprot", args.miniprot, [])
    hmmsearch = _resolve("hmmsearch", args.hmmsearch, ["/usr/bin/hmmsearch"])

    gids = [g.strip() for g in args.genome_ids.split(",") if g.strip()]
    os.makedirs(args.outdir, exist_ok=True)
    gff_dir = os.path.join(args.outdir, "gff")
    prot_dir = os.path.join(args.outdir, "tert_proteins")
    os.makedirs(gff_dir, exist_ok=True)
    os.makedirs(prot_dir, exist_ok=True)

    workdir = os.path.join(args.outdir, "_work")
    os.makedirs(workdir, exist_ok=True)

    # resolve genome fastas up front
    genome_fa = {}
    for gid in gids:
        try:
            src = find_genome_fasta(gid, args.refseq_dir, args.tara_dir)
        except FileNotFoundError:
            print(f"skip {gid}: FASTA not found", file=sys.stderr)
            continue
        genome_fa[gid] = ensure_plain_fasta(src, workdir, gid)

    # seed set grows across rounds
    seeds_path = os.path.join(workdir, "seeds_round.faa")
    shutil.copyfile(args.seeds, seeds_path)
    seen_seed_seqs = set()

    all_candidate_rows = []
    for rnd in range(1, args.iterate + 1):
        with open(seeds_path) as _sf:
            n_seeds = sum(1 for _ in _sf if _.startswith(">"))
        print(f"\n===== ROUND {rnd} (seeds: {n_seeds}) =====")
        round_confirmed_proteins = []
        for gid, fa in genome_fa.items():
            cands = run_miniprot_trans(miniprot, fa, seeds_path, args.threads, args.outc, args.outn)
            if rnd == 1:
                run_miniprot_gff(miniprot, fa, args.seeds, args.threads,
                                 os.path.join(gff_dir, f"{gid}.miniprot.gff"),
                                 args.outc, args.outn)
            if not cands:
                print(f"  {gid}: 0 candidates")
                continue
            # write candidate proteins, hmmsearch both domains
            with tempfile.TemporaryDirectory(dir=workdir) as td:
                faa = os.path.join(td, "cand.faa")
                names = [f"{gid}__r{rnd}c{i}" for i in range(len(cands))]
                write_faa(list(zip(names, [c["protein"] for c in cands])), faa)
                trbd = hmmsearch_best(hmmsearch, args.trbd_hmm, faa, td, "trbd")
                rt = hmmsearch_best(hmmsearch, args.rt_hmm, faa, td, "rt")

            ngc = 0
            for name, c in zip(names, cands):
                t = trbd.get(name)
                r = rt.get(name)
                has_trbd = bool(t and t[0] < args.trbd_evalue)
                has_rt = bool(r and r[0] < args.rt_evalue)
                # architecture: TRBD must start N-terminal of RT (when both present)
                order_ok = bool(t and r and t[1] < r[1])
                strict = has_trbd and has_rt and order_ok
                all_candidate_rows.append({
                    "genome_id": gid, "round": rnd, "cand": name,
                    "contig": c["contig"], "start": c["start"], "end": c["end"],
                    "strand": c["strand"], "locus_len": c["end"] - c["start"],
                    "protein_len": len(c["protein"]),
                    "best_seed": c["seed"], "seed_qcov": round((c["qend"] - c["qstart"]) / c["qlen"], 3),
                    "trbd_evalue": t[0] if t else float("nan"), "trbd_ali": f"{t[1]}-{t[2]}" if t else "",
                    "rt_evalue": r[0] if r else float("nan"), "rt_ali": f"{r[1]}-{r[2]}" if r else "",
                    "has_trbd": has_trbd, "has_rt": has_rt, "strict": strict,
                    "protein": c["protein"],
                })
                if strict:
                    ngc += 1
                    round_confirmed_proteins.append((name, c["protein"]))
            print(f"  {gid}: {len(cands)} candidates, {ngc} confirmed (TRBD+RT, N->C)")

        # augment seeds for the next round with this round's confirmed proteins
        if rnd < args.iterate and round_confirmed_proteins:
            with open(seeds_path, "a") as fh:
                added = 0
                for name, seq in round_confirmed_proteins:
                    if seq in seen_seed_seqs:
                        continue
                    seen_seed_seqs.add(seq)
                    fh.write(f">{name}\n{seq}\n")
                    added += 1
            print(f"  augmented seeds with {added} confirmed proteins for round {rnd+1}")
        elif rnd < args.iterate:
            print("  no confirmed proteins to augment with; stopping early")
            break

    if not all_candidate_rows:
        print("\nNo TERT candidates at all.", file=sys.stderr)
        pd.DataFrame(columns=["genome_id", "contig", "start", "end", "strand"]).to_csv(
            os.path.join(args.outdir, "confirmed_tert.tsv"), sep="\t", index=False)
        open(os.path.join(args.outdir, "all_confirmed.faa"), "w").close()
        return

    df = pd.DataFrame(all_candidate_rows)
    # Group overlapping candidates into loci per genome+contig, then assign a tier:
    #   strict  : a single protein has TRBD + RT in N->C order
    #   partial : RT confirmed (PF00078) AND >=2 distinct seeds map to the locus,
    #             but no single protein clears the TRBD HMM (divergent/fragmented TRBD)
    # Representative = strict member (best TRBD E) if any, else the longest RT-bearing
    # member (best RT E). Distinct-seed count uses the un-augmented apicomplexan seeds
    # plus round-1 within-genus proteins as independent evidence.
    df = df.sort_values(["genome_id", "contig", "start"])
    locus_recs = []
    for (gid, contig), sub in df.groupby(["genome_id", "contig"]):
        sub = sub.sort_values("start")
        grp, last_end, groups = -1, -1, []
        for _, row in sub.iterrows():
            if row["start"] > last_end:
                grp += 1
            groups.append(grp)
            last_end = max(last_end, row["end"])
        sub = sub.assign(_g=groups)
        for _, g in sub.groupby("_g"):
            n_seeds = g["best_seed"].nunique()
            strict_members = g[g["strict"]]
            rt_members = g[g["has_rt"]]
            if len(strict_members):
                rep = strict_members.sort_values(["trbd_evalue", "protein_len"],
                                                 ascending=[True, False]).iloc[0]
                tier = "strict_TRBD+RT"
            elif len(rt_members) and n_seeds >= 2:
                rep = rt_members.sort_values(["rt_evalue", "protein_len"],
                                             ascending=[True, False]).iloc[0]
                tier = "partial_RT+synteny"
            else:
                continue  # weak/noise locus
            locus_recs.append({**rep.to_dict(), "tier": tier,
                               "n_seeds_at_locus": int(n_seeds),
                               "locus_start": int(g["start"].min()),
                               "locus_end": int(g["end"].max())})

    if not locus_recs:
        print("\nNo TERT loci passed strict or partial tiers.", file=sys.stderr)
        pd.DataFrame(columns=["genome_id", "contig", "start", "end", "strand"]).to_csv(
            os.path.join(args.outdir, "confirmed_tert.tsv"), sep="\t", index=False)
        open(os.path.join(args.outdir, "all_confirmed.faa"), "w").close()
        return

    best_df = pd.DataFrame(locus_recs)
    # chosen-per-genome: prefer strict over partial; within tier, best TRBD then RT then length
    tier_rank = {"strict_TRBD+RT": 0, "partial_RT+synteny": 1}
    best_df["_tier_rank"] = best_df["tier"].map(tier_rank)
    best_df["_trbd_sort"] = best_df["trbd_evalue"].fillna(1.0)
    best_df["_rt_sort"] = best_df["rt_evalue"].fillna(1.0)
    best_df["chosen"] = False
    for gid, g in best_df.groupby("genome_id"):
        idx = g.sort_values(["_tier_rank", "_trbd_sort", "_rt_sort", "protein_len"],
                            ascending=[True, True, True, False]).index[0]
        best_df.loc[idx, "chosen"] = True

    cols = ["genome_id", "contig", "locus_start", "locus_end", "strand", "tier",
            "n_seeds_at_locus", "protein_len", "best_seed", "seed_qcov",
            "trbd_evalue", "trbd_ali", "rt_evalue", "rt_ali", "round", "chosen"]
    best_df.sort_values(["chosen", "_tier_rank", "_trbd_sort"],
                        ascending=[False, True, True])[cols].to_csv(
        os.path.join(args.outdir, "confirmed_tert.tsv"), sep="\t", index=False)

    all_conf = []
    for _, row in best_df.iterrows():
        te = f"{row.trbd_evalue:.0e}" if pd.notna(row.trbd_evalue) else "NA"
        nm = f"{row.genome_id}|{row.contig}:{row.locus_start}-{row.locus_end}|{row.strand}|{row.tier}|trbdE={te}"
        all_conf.append((nm, row["protein"]))
        if row.chosen:
            write_faa([(nm, row["protein"])], os.path.join(prot_dir, f"{row.genome_id}.faa"))
    write_faa(all_conf, os.path.join(args.outdir, "all_confirmed.faa"))

    print("\n=== TERT loci (chosen per genome) ===")
    for _, row in best_df[best_df.chosen].sort_values(["_tier_rank", "genome_id"]).iterrows():
        te = f"{row.trbd_evalue:.1e}" if pd.notna(row.trbd_evalue) else "NA"
        re_ = f"{row.rt_evalue:.1e}" if pd.notna(row.rt_evalue) else "NA"
        print(f"  [{row.tier:18s}] {row.genome_id}: {row.contig}:{row.locus_start}-{row.locus_end} "
              f"({row.strand}) prot={row.protein_len}aa TRBD_E={te} RT_E={re_} nseeds={row.n_seeds_at_locus}")
    nstrict = int((best_df.chosen & (best_df.tier == "strict_TRBD+RT")).sum())
    npartial = int((best_df.chosen & (best_df.tier == "partial_RT+synteny")).sum())
    print(f"\nwrote {args.outdir}/confirmed_tert.tsv "
          f"({best_df.chosen.sum()} genomes: {nstrict} strict, {npartial} partial; {len(best_df)} total loci)")


if __name__ == "__main__":
    main()
