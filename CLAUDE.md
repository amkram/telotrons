# CLAUDE.md

Guidance for Claude Code (claude.ai/code) when working in this repository.

## What this repo is

Research code for **telotrons** — introns composed of tandem telomeric repeats
(TTAGGG, TTTAGGG, TTAGG, …), hypothesized to form when telomerase heals a
double-strand break inside coding sequence and the inserted repeat array is then
spliced out. The pipeline surveys eukaryotic genomes (NCBI RefSeq + Tara Oceans
SMAGs) for telomeric-repeat arrays inside annotated introns, then runs downstream
analyses (architecture, motif discovery, telomerase/TERT homology,
interstitial-array comparison, MSA) and emits figures.

Everything from the prior generation of the project lives in [old/](old/) — a
read-only reference archive (genomes, manuscript drafts, ED figures, audit
notes, `pan_euk_telotrons/`, `tara_oceans_euk_mags/`, `AGENT_HANDOFF.md`). Treat
`old/` as archaeology, not as a place to import from or add code.

## Pipeline (Snakefile + scripts/ + config.yaml)

One Snakemake workflow (~1500 lines, ~60 rules) configured by
[config.yaml](config.yaml) (`configfile:`). It is **not** linear: a scan/filter
core fans out into many independent analysis and plotting arms. `scripts/` holds
~36 standalone Python CLIs plus a few `.sh` helpers; each rule shells out to one.
There are no tests or lints — iterate by running a stage directly with the same
args the rule uses.

### Core survey path
1. `manifests` / `refseq_urls` / `tara_archives` / `download_refseq` — build
   `manifests/all_genomes.tsv` (`genome_id, organism, group, source`) and
   download FASTAs+GFFs (RefSeq via `xargs -P`).
2. `canonical_motifs` — emit the curated per-genome/per-group telomere-motif
   table from config (overrides contig-end motif inference where set).
3. `scan_all` → [scripts/scan_telotrons.py](scripts/scan_telotrons.py) — derive
   introns (`gt gff3 -addintrons`) and scan each for the configured motifs
   (rotations + reverse complements) via `seqkit locate`; emit three TSVs
   (per-intron, candidate loci, per-species summary). Python only composes the
   per-locus record (orientation, splice signals, terminal-motif inference with
   canonical override, flank fraction).
4. `filter_final` → [scripts/filter_final_set.py](scripts/filter_final_set.py) —
   two admission pathways (single-array `filter.min_repeat_frac`=0.85;
   bidirectional `bidir_min_repeat_frac`=0.40 + `bidir_min_hits`=3) plus
   `require_terminal_motif_match`, `collapse_unique_loci`, optional
   `require_canonical_splice`. Splits positives from zero-telotron negative controls.
5. `dedup_telotrons` / `classify_architecture` — collapse doubly-assembled /
   close-paralog loci (flank blastn + union-find); classify splice architecture
   (GT-F-AG etc.) and emit per-architecture boundary k-mers.
6. `analyze` → [scripts/analyze_telotrons.py](scripts/analyze_telotrons.py) —
   boundary k-mer enrichment vs control introns, distance-to-contig-end test,
   architecture summary.
7. `figures` → [scripts/plot_telotrons.py](scripts/plot_telotrons.py) — counts +
   per-species panels. `package` — zip final TSVs + figures.

### Downstream analysis arms (mutually independent)
- **Extraction / alignment** — `extract_telotron_fasta`, `build_non_telotron_controls`
  + `extract_non_telotron_fasta`, `msa_telotron_regions` (MAFFT), `blast_linkers`.
- **Interstitial arrays** — `make_unannotated_masks` → `find_interstitial_arrays`
  + splice-candidate filtering (telomeric arrays *outside* introns, as a contrast).
- **Motif discovery** (`envs/meme.yaml`) — STREME on telotrons / non-telo introns /
  linkers / branchpoints; FIMO branchpoint scan.
- **Telomerase / TERT homology** — `fetch_telomerase_db` (UniProt Swiss-Prot),
  BLAST telomerase vs genomes; TERT deep-homology search (`fetch_tert_seeds_hmms`
  → `find_tert`, miniprot + Pfam TRBD/RVT).
- **Figures** — ~20 plotting rules: boundary-kmer plots, splice/sequence logos
  (telotron, control, composite, by-architecture, by-5′/3′-category), array-length
  distributions, terminal-motif density, pipeline-stage diagram.

### Running it
```bash
# Managed deps (builds envs/telotrons.yaml + envs/meme.yaml).
snakemake --use-conda -j 16

# Or with tools on PATH: gt (GenomeTools), seqkit, bedtools, samtools, blast+,
# mafft, miniprot, hmmer, the MEME suite, python + pandas/scipy/matplotlib.
snakemake -j 16

snakemake -n                                     # dry run
snakemake --use-conda -j 8 --forcerun scan_all   # single-stage rerun
```
Config (threads, `refseq_groups`, `accessions` whitelist, `telomere_motifs`,
`canonical_telomere_motifs`, `scan`/`filter` cutoffs, data-source URLs) lives in
`config.yaml`; override with `--config k=v` or `--configfile`. Hand-run drivers
also exist in `scripts/` (`run_full_survey.sh` for the production survey,
`run_test_survey.sh` for a small test manifest).

## Non-obvious constraints (still apply)

- **Misannotation is the dominant false positive.** Gene predictors split long
  telomeric arrays into fake exon/intron structure. `scan_telotrons.py` rejects
  loci with `flank_telomeric_frac > scan.max_flank_repeat_frac` (default 0.50).
  Don't loosen without thinking.
- **`require_terminal_motif_match` is doing real work** — it forces the intronic
  motif to equal the motif dominating that genome's contig ends (the real
  telomere), killing cross-motif noise.
- **Repeat type is species-specific.** TTAGGG (vertebrates/fungi) vs TTTAGGG
  (plants, Apicomplexa) vs TTAGG (insects), etc. The curated list is in
  `config.yaml` (`telomere_motifs`) — searching only TTAGGG misses major lineages.
- **NCBI download concurrency.** `download_refseq` fans out with `xargs -P`; keep
  threads ≤ ~8 against NCBI to avoid rate limiting.
- **Human intronic telomeric insertions are passengers, not telotrons.** CCCTAA
  lacks GT/AG, so reverse-strand insertions can't create splice signals — they
  land inside pre-existing introns.

## External data sources
- **NCBI RefSeq** assembly summary and **Tara Oceans SMAGs v1** — URLs in
  `config.yaml` (`refseq_url`, `tara_base`).
- **UniProt Swiss-Prot** FASTA for the telomerase DB — hard-coded in the
  `fetch_telomerase_db` rule.

## old/
Reference archive only; don't import from it. Notable: `old/AGENT_HANDOFF.md`,
`old/MANUSCRIPT_AUDIT_REPORT.md`, `old/pan_euk_telotrons/` (v1/v2 surveys,
ULTRA-based scanner, validated `real_telotrons/`), `old/tara_oceans_euk_mags/`
(606-MAG discovery pipeline), many ED figures and draft `.docx`. Also
`old/_deslop_2026-05-29/` — items pulled out of the live tree during the
2026-05-29 cleanup (a source backup, redundant AF3 zips, two superseded one-off
scripts, a stray `pangraph` binary).
