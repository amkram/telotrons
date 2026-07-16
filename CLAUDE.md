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

Everything from the prior generation of the project lives in
[work/old/](work/old/) — a read-only reference archive (genomes, manuscript
drafts, ED figures, audit notes, `pan_euk_telotrons/`, `tara_oceans_euk_mags/`,
`AGENT_HANDOFF.md`). Treat `work/old/` as archaeology, not as a place to import
from or add code.

**Layout (post-2026-06-05 deslop).** Exactly four root subdirs:
- `scripts/` — pipeline source (~74 files). Every Snakefile rule shells out to
  one of these. `_common.py` and `telomere_mask.py` are shared helpers
  (the latter closes the telomere-rotation contamination trap — always import
  before any composition/motif/logo analysis on telomere-adjacent sequence).
- `envs/` — conda env YAMLs (`telotrons.yaml`, `meme.yaml`). Materialised envs
  live in `~/.snakemake-envs/` and are gitignored.
- `data/raw/` — all downloaded inputs (NCBI RefSeq + Tara SMAGs +
  telomerase.us + long-read FASTQs). Snakefile rules and CLI args reference
  `data/raw/refseq` / `data/raw/tara` / etc.
- `work/` — every pipeline output, log, manifest, and the historical archive:
  - `work/results/` — Snakemake rule outputs (TSVs + figures)
  - `work/manifests/` — generated manifests
  - `work/logs/` — runtime logs from snakemake invocations
  - `work/old/` — read-only archive: prior generations, manuscript drafts,
    audit snapshots, decoupled experiments (see "work/old/" section below).

`Snakefile`, `config.yaml`, `CLAUDE.md`, `paper.md`, `.gitignore` are the only
tracked files at root. No shell orchestrators — invoke via `snakemake`.

## Pipeline (Snakefile + scripts/ + config.yaml)

One Snakemake workflow (~1700 lines, 67 rules) configured by
[config.yaml](config.yaml) (`configfile:`). It is **not** linear: a scan/filter
core fans out into many independent analysis and plotting arms. `scripts/` holds
~70 standalone Python CLIs; each rule shells out to one. **The Snakefile is
canonical** — all prior shell orchestrators (`run_full_survey.sh`,
`run_test_survey.sh`, `_survey_env.sh`, `blast_by_arch.sh`,
`restriction_factor_sweep.sh`, `apply_good_orthologs.sh`) were retired in the
2026-06-05 deslop because the Snakefile covers every stage they wrapped.

There are no tests or lints — iterate by running a stage directly with the same
args the rule uses (or `snakemake -n --rulegraph` for DAG inspection).

### Core survey path
1. `manifests` / `tara_archives` / `download_assemblies` — build
   `work/manifests/all_genomes.tsv` (columns:
   `genome_id, organism, group, ftp_path, source`) from three streams:
   **RefSeq** (curated GCF_), **GenBank** (annotated eukaryotes without a
   paired GCF_, `annotation_provider != "na"`), and **Tara SMAGs v1**. RefSeq
   downloads land in `data/raw/refseq/`; GenBank-only in `data/raw/genbank/`.
   URL derivation is inlined into `download_assemblies`.
2. `scan_all` → [scripts/scan_telotrons.py](scripts/scan_telotrons.py) — derive
   introns (`gt gff3 -addintrons`) and scan each for the configured motifs
   (rotations + reverse complements) via `seqkit locate`. Builds the canonical
   per-genome motif TSV inline from config. Emits per-intron, candidate loci,
   per-species summary, and canonical_motifs.tsv.
3. `filter_final` → [scripts/filter_final_set.py](scripts/filter_final_set.py) —
   two admission pathways (single-array `filter.min_repeat_frac`=0.85;
   bidirectional `bidir_min_repeat_frac`=0.40 + `bidir_min_hits`=3) plus
   `require_terminal_motif_match`, `collapse_unique_loci`, optional
   `require_canonical_splice`. Splits positives from zero-telotron negative controls.
4. `dedup_telotrons` / `classify_architecture` — collapse doubly-assembled /
   close-paralog loci (flank blastn + union-find); classify splice architecture
   (GT-F-AG etc.) and emit per-architecture boundary k-mers.
5. `analyze` → [scripts/analyze_telotrons.py](scripts/analyze_telotrons.py) —
   boundary k-mer enrichment vs control introns, distance-to-contig-end test,
   architecture summary.
6. `confident_species` → [scripts/confident_species.py](scripts/confident_species.py)
   — emit `work/results/confident_species.tsv`, the paper's central bearer set.
   A species is admitted when it has **≥3 telotrons passing filter_final**
   (`confident_species.min_n`) OR **≥2 bidirectional architectures**
   (`confident_species.min_bidir`; GT-F-R-AG, a linker variant, or
   Multi-junction — a distinctive telomerase-mediated signature). Both
   thresholds live in `config.yaml` and are wired through the Snakefile;
   `min_bidir` is 2 rather than 1 because motif rotations alone can
   occasionally hit both strands at `bidir_min_repeat_frac`=0.40, and
   noise-driven singletons must not seed downstream figures at the
   RefSeq+GenBank ~50k-assembly scale. **Every downstream analysis keys off
   this file** — new bearer species flow through automatically without
   touching Python.
7. `package` — zip final TSVs + confident-species set + manifest into the
   deliverable at `work/results/telotron_pipeline_outputs.zip`.

### Downstream analysis arms (mutually independent — invoke by rule name)
- **Extraction / control** — `extract_fasta` (one wildcarded rule for both
  `telotron` and `non_telotron` sets) + `build_non_telotron_controls`.
- **Interstitial arrays** — `find_interstitial_arrays` (builds its 6-frame ORF
  mask inline; telomeric arrays *outside* introns as a contrast).
- **TERT homology** — `fetch_tert_seeds_hmms` → `find_tert` (miniprot + Pfam
  TRBD/RVT deep-homology search). Recovers Eimeria and sister coccidia TERT
  below BLAST detection.
- **Telotron-gene orthologs** — `telotron_orthologs`: miniprot-maps the
  orthologous locus in a panel of telotron-negative sisters + cross-Eimeria,
  aligns the gene in protein space (introns removed), DNA-aligns the telotron
  against the orthologous intron; emits the compiled per-locus PDF in the
  same rule. Config: `telotron_ortholog` (focal_ids, ortholog_ids, cutoffs).
- **Architecture / linker** — `linker_analysis` (segment telotrons into arrays
  + linkers, cluster linkers by 7-mer Jaccard for cross-locus recurrence),
  `mask_telotron_arrays` (G/A/L architecture cartoon + MSA).
- **Analysis arm** (`analysis_arm` aggregate) —
  `nucleosome_analysis` (insertion-site composition / 10-bp WW periodicity /
  CpG panel + within-gene sibling-intron control, one rule),
  `telotron_gene_class` (host vs disjoint non-host + per-intron dissection,
  iterates confident_species dynamically), `length_distribution_by_arch`
  (BH-corrected, single-MAG caveat; emits per-arm burst-length figure too),
  `telotron_expr_figures` (tenella + necatrix, size-controlled OLS +
  per-intron rate vs expression quintile) + `rnaseq_gene_coverage`
  (per-species SRA→`samtools bedcov`, config-driven).

### Running it
```bash
# Managed deps (builds envs/telotrons.yaml + envs/meme.yaml).
snakemake --use-conda -j 16

# Or with tools on PATH: gt (GenomeTools), seqkit, bedtools, samtools, blast+,
# mafft, miniprot, hmmer, the MEME suite, python + pandas/scipy/matplotlib.
snakemake -j 16

snakemake -n                                     # dry run
snakemake --use-conda -j 8 --forcerun scan_all   # single-stage rerun

# Test on the 22-species set: set `accessions: [...]` in config.yaml or:
snakemake -j 8 --config 'accessions=["GCF_000499545.2","GCF_000499605.1",...]'
```
Config (threads, `refseq_groups`, `accessions` whitelist, `telomere_motifs`,
`canonical_telomere_motifs`, `scan`/`filter` cutoffs, data-source URLs) lives in
`config.yaml`; override with `--config k=v` or `--configfile`.

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
- **NCBI download concurrency.** `download_assemblies` fans out with `xargs -P`; keep
  threads ≤ ~8 against NCBI to avoid rate limiting.
- **Human intronic telomeric insertions are passengers, not telotrons.** CCCTAA
  lacks GT/AG, so reverse-strand insertions can't create splice signals — they
  land inside pre-existing introns.

## External data sources
- **NCBI RefSeq + GenBank** assembly summaries and **Tara Oceans SMAGs v1** —
  URLs in `config.yaml` (`refseq_url`, `genbank_url`, `tara_base`).

## work/old/
Reference archive only; don't import from it. Aggressive-prune pass on
2026-07-14 removed the 5 largest historical dirs (`pan_euk_telotrons/` 408G,
`good_set/` 172G, `analysis_2026-06/` 81G, `tara_oceans_euk_mags/` 28G,
`comparative_genomics_2026-06/` 5G) — total ~694G reclaimed; findings live in
memory files and prior commit history. Notable surviving subtrees:
- `work/old/AGENT_HANDOFF.md`, `MANUSCRIPT_AUDIT_REPORT.md` (top-level audit notes)
- `work/old/eimeria_rnaseq/` (22G), `work/old/toxo_rnaseq/` (16G) — RNA-seq
  workspaces used for the expression arm's original hand-runs.
- `work/old/its_comparison/` (5.5G) — cross-clade ITS load comparison inputs.
- `work/old/figs/` (3.7G) — legacy figure workspaces (superseded by
  `work/results/figures/`).
- `work/old/_deslop_2026-05-29/` (244M) — 2026-05-29 cleanup (source backup,
  redundant AF3 zips, superseded scripts, stray `pangraph` binary).
- `work/old/_deslop_2026-06-04/scripts/` (4.9G total dir) — 67 one-off
  mechanism-deep-dive scripts (kill-tests `k1–k4`, `q2/q3`, `probe_*`,
  `subtelo_*`, one-off `plot_*`/`telotron_intron_*`).
- `work/old/audit_2026-06/` (310M) — 2026-06-05 deslop: dated audit snapshots
  (`register_lock`, `pipeline_review`, `pipeline_fixes`, `intron_overwrite`,
  `artifact_2026-05-29` formerly `.artifact/` at root).
- `work/old/paper_backups/` — pre-round manuscript backups (adversarial-review,
  MAG-highlight, round-4).
