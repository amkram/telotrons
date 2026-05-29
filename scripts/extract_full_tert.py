#!/usr/bin/env python3
"""Recover the longest homology-supported TERT protein per genome by translating
the COMPLETE miniprot gene-model CDS (not the frameshift-split --trans segments).

Strategy
    1. Run miniprot --gff with an augmented seed set (apicomplexan TERTs + any
       already-recovered coccidian proteins). The closest/longest seed extends
       the gene model furthest toward the N-terminal TEN domain.
    2. For each target locus (from confirmed_tert.tsv) pick the gene model with
       the largest aligned span (then lowest seed qstart = most N-terminal).
    3. Extract that model's CDS exons, splice, reverse-complement if needed, and
       translate the full ORF.
    4. Optionally iterate: the longest proteins re-seed the next round so
       N-terminal coverage propagates across the clade.

The TEN domain is fast-evolving; cross-genus homology may not reach the true
5' end. Output reports genomic span, seed N-terminal coverage, and domain
content so completeness can be judged honestly.
"""
import argparse
import glob
import gzip
import os
import subprocess as sp
import tempfile

import pandas as pd

CODON = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
ENV_BIN = "/home/alex/.snakemake-envs/telotrons/7aa0f2645427a08bdf054e56b1233123_/bin"


def which(name, fallbacks):
    import shutil
    return shutil.which(name) or next((f for f in fallbacks if os.path.exists(f)), None)


def rc(s):
    return s.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]


def translate(nt):
    aa = []
    for i in range(0, len(nt) - 2, 3):
        aa.append(CODON.get(nt[i:i+3].upper(), 'X'))
    return "".join(aa)


def load_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    seqs, sid, buf = {}, None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if sid is not None:
                    seqs[sid] = "".join(buf).upper()
                sid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if sid is not None:
            seqs[sid] = "".join(buf).upper()
    return seqs


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


def parse_gff_models(gff_path):
    """Return list of models: {id, contig, strand, qstart, qend, cds:[(s,e),...]}."""
    models = {}
    order = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] == "mRNA":
                mid = attr.get("ID")
                tgt = attr.get("Target", "").split()
                qs = int(tgt[1]) if len(tgt) >= 3 else 0
                qe = int(tgt[2]) if len(tgt) >= 3 else 0
                models[mid] = {"id": mid, "contig": f[0], "strand": f[6],
                               "qstart": qs, "qend": qe, "cds": []}
                order.append(mid)
            elif f[2] == "CDS":
                mid = attr.get("Parent")
                if mid in models:
                    models[mid]["cds"].append((int(f[3]), int(f[4])))
    return [models[m] for m in order]


def translate_model(model, chrom):
    cds = sorted(model["cds"])
    nt = "".join(chrom[s - 1:e] for s, e in cds)
    if model["strand"] == "-":
        nt = rc(nt)
    prot = translate(nt)
    return prot


def hmm_evalue(hmmsearch, hmm, faa, tmpdir, tag):
    if os.path.getsize(faa) == 0:
        return None
    dt = os.path.join(tmpdir, f"{tag}.domtbl")
    sp.run([hmmsearch, "--max", "-E", "10", "--domtblout", dt, hmm, faa],
           check=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
    best = None
    with open(dt) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) >= 13:
                e = float(f[12])
                best = e if best is None else min(best, e)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--confirmed", required=True, help="confirmed_tert.tsv (chosen loci)")
    ap.add_argument("--seeds", required=True, help="augmented seed FASTA (apicomplexan + recovered cores)")
    ap.add_argument("--trbd-hmm", required=True)
    ap.add_argument("--rt-hmm", required=True)
    ap.add_argument("--refseq-dir", required=True)
    ap.add_argument("--tara-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--iterate", type=int, default=2)
    ap.add_argument("--threads", type=int, default=16)
    args = ap.parse_args()

    miniprot = which("miniprot", [f"{ENV_BIN}/miniprot"])
    hmmsearch = which("hmmsearch", ["/usr/bin/hmmsearch"])
    os.makedirs(args.outdir, exist_ok=True)
    prot_dir = os.path.join(args.outdir, "tert_proteins_full")
    os.makedirs(prot_dir, exist_ok=True)

    conf = pd.read_csv(args.confirmed, sep="\t")
    conf = conf[conf["chosen"]] if "chosen" in conf.columns else conf
    loci = {r.genome_id: (r.contig, int(r.locus_start), int(r.locus_end))
            for r in conf.itertuples(index=False)}

    seeds = args.seeds
    rows = []
    for rnd in range(1, args.iterate + 1):
        nseed = sum(1 for _ in open(seeds) if _.startswith(">"))
        print(f"\n===== FULL-LENGTH ROUND {rnd} (seeds: {nseed}) =====")
        round_prots = []
        rows = []
        for gid, (contig, ls, le) in loci.items():
            try:
                fa = find_genome_fasta(gid, args.refseq_dir, args.tara_dir)
            except FileNotFoundError:
                print(f"  {gid}: genome not found"); continue
            chrom = load_fasta(fa).get(contig, "")
            if not chrom:
                print(f"  {gid}: contig {contig} not found"); continue
            with tempfile.TemporaryDirectory(dir=args.outdir) as td:
                gff = os.path.join(td, "m.gff")
                proc = sp.run([miniprot, "-t", str(args.threads), "--gff",
                               "--outc=0.05", "--outn=20", fa, seeds],
                              check=True, stdout=sp.PIPE, stderr=sp.DEVNULL, text=True)
                open(gff, "w").write(proc.stdout)
                models = [m for m in parse_gff_models(gff)
                          if m["contig"] == contig and not (m["cds"] and
                              (max(e for _, e in m["cds"]) < ls - 5000 or
                               min(s for s, _ in m["cds"]) > le + 5000))]
                if not models:
                    print(f"  {gid}: no model at locus"); continue
                # score: aligned aa span (qend-qstart), then most N-terminal (low qstart)
                models.sort(key=lambda m: (-(m["qend"] - m["qstart"]), m["qstart"]))
                best = models[0]
                prot = translate_model(best, chrom).strip("X")
                # trim trailing stop, keep internal (frameshift) info
                clean = prot.rstrip("*")
                ninternal_stop = clean.count("*")
                # validate domains
                faa = os.path.join(td, "p.faa")
                open(faa, "w").write(f">{gid}\n{clean.replace('*','X')}\n")
                trbd = hmm_evalue(hmmsearch, args.trbd_hmm, faa, td, "trbd")
                rt = hmm_evalue(hmmsearch, args.rt_hmm, faa, td, "rt")
                gs = min(s for s, _ in best["cds"]); ge = max(e for _, e in best["cds"])
                rows.append({"genome_id": gid, "contig": contig,
                             "gene_start": gs, "gene_end": ge, "strand": best["strand"],
                             "prot_len": len(clean), "internal_stops": ninternal_stop,
                             "seed_qstart": best["qstart"], "seed_qend": best["qend"],
                             "trbd_evalue": trbd, "rt_evalue": rt,
                             "protein": clean.replace("*", "X")})
                round_prots.append((gid, clean.replace("*", "X")))
                print(f"  {gid}: prot={len(clean)}aa span={gs}-{ge} qstart={best['qstart']} "
                      f"TRBD={trbd:.0e} RT={rt} stops={ninternal_stop}"
                      if trbd else f"  {gid}: prot={len(clean)}aa (no TRBD)")
        # augment seeds with this round's full proteins for the next round
        if rnd < args.iterate and round_prots:
            aug = os.path.join(args.outdir, f"_seeds_round{rnd+1}.faa")
            with open(aug, "w") as o:
                for line in open(seeds):
                    o.write(line)
                for gid, p in round_prots:
                    o.write(f">FULL_{gid}\n{p}\n")
            seeds = aug

    df = pd.DataFrame(rows)
    df.drop(columns=["protein"]).to_csv(os.path.join(args.outdir, "full_tert.tsv"),
                                        sep="\t", index=False)
    for r in rows:
        open(os.path.join(prot_dir, f"{r['genome_id']}.faa"), "w").write(
            f">{r['genome_id']}|{r['contig']}:{r['gene_start']}-{r['gene_end']}|"
            f"{r['strand']}|fulllen={r['prot_len']}aa\n{r['protein']}\n")
    with open(os.path.join(args.outdir, "all_full.faa"), "w") as o:
        for r in rows:
            o.write(f">{r['genome_id']}\n{r['protein']}\n")
    print(f"\nwrote {args.outdir}/full_tert.tsv and tert_proteins_full/ ({len(rows)} proteins)")


if __name__ == "__main__":
    main()
