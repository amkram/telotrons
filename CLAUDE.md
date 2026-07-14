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
1. `manifests` / `refseq_urls` / `tara_archives` / `download_refseq` — build
   `work/manifests/all_genomes.tsv` (`genome_id, organism, group, source`) and
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
  linkers. (Branchpoint STREME/FIMO arm retired 2026-07-14 — never produced a
  usable per-species PWM at coccidian divergence.)
- **TERT homology** — `fetch_tert_seeds_hmms` → `find_tert` (miniprot + Pfam
  TRBD/RVT deep-homology search). Recovers Eimeria and sister coccidia TERT
  below BLAST detection. (The old proteome-BLAST arm — `fetch_telomerase_db`,
  `blast_telomerase_vs_genomes`, `process_telomerase_blast`,
  `find_gene_deep_homology` — was retired 2026-07-14; its per-genome
  BLAST-based i-Evalues were superseded by the HMM-based `find_tert` results.)
- **Telotron-gene orthologs** (`orthologs` target) — `telotron_ortholog_align`: for each
  telotron-host gene, miniprot-maps the orthologous locus in a panel of telotron-negative
  sisters + cross-Eimeria, aligns the gene in protein space (introns removed), and DNA-aligns
  the telotron against the orthologous intron where present. Resolves fill-vs-create per locus;
  flags (does **not** correct) frameshifts/stops at the homologous junction as bad intron-
  boundary alignment. Config: `telotron_ortholog` (focal_ids, ortholog_ids, cutoffs).
  Two viewers consume its output: `plot_telotron_ortholog_loci` (compiled per-locus PDF:
  flanking-exon protein MSA + flank DNA + intron DNA, poorly-flank-aligned orthologs dropped)
  and `telotron_ortholog_textdump` (per-locus text files — unaligned/aligned × DNA/aa for
  flanks+intron — grouped into `locus_text/{intron_present_nontelo,telotron_present,intron_absent,uncertain}/`).
- **Figures** — ~20 plotting rules: boundary-kmer plots, splice/sequence logos
  (telotron, control, composite, by-architecture, by-5′/3′-category), array-length
  distributions, terminal-motif density, pipeline-stage diagram.
- **Architecture / linker / interstitial-ITS** (`architecture_analyses` target, added
  2026-06-04) — `interstitial_ortholog_textdump` (non-genic ITS DNA-flank orthology,
  gold-standard ≥4-unit/<1mm-per-unit filter; chains off `find_interstitial_arrays`),
  `linker_segmentation` + `cluster_linkers` (telotron array/linker decomposition),
  `mask_telotron_arrays` (G/A/L architecture cartoon + MSA), `ortholog_review_html`
  (`build_ortholog_review.py` interactive reviewer). The linker/mask/review rules read
  `work/results/telotron_orthologs_v2/locus_text` (a **manual** v2 ortholog run — reconcile
  to fully connect the DAG). `curate_locus_text.py` (in-place review-decision curation) and
  `wrap_locus_text.py` are human-in-the-loop tools, run by hand not as DAG rules.
- **Analysis arm** (`analysis_arm` target, wired 2026-06-16) — the previously hand-run
  downstream scripts whose inputs are regenerable from the core pipeline: `nucleosome_features`
  (telomere-MASKED composition/periodicity panel + BH-FDR), `nucleosome_withingene` (within-gene
  sibling-intron control), `telotron_gene_bias` (host vs **disjoint** non-host gene class),
  `telotron_per_intron`, `length_distribution_by_arch` + `length_per_arm_figure` (BH-corrected,
  single-MAG caveat), `mechanism_diagrams`, and the RNA-seq **expression arm**
  (`telotron_expr_figures` + `rnaseq_gene_coverage`: per-species SRA→`samtools bedcov` gene coverage
  driven by `config["rnaseq"]`, writing `work/results/rnaseq/{species}_gene_cov.tsv` — replaces the old
  hand-run `run_pipeline.sh` + ephemeral `/tmp/eten_gene_cov.tsv`; the optional locus-level splice panel
  still needs the manual `data/raw/rnaseq_splice_2026/per_locus_counts.tsv` and is skipped if absent).
  **Not wired** (need external/non-regenerable inputs — documented in a Snakefile note after the
  `analysis_arm` rule): the
  Hi-C arm (`hic_*` — ENA FASTQs + cooltools), the ONT arm (`ont_*` — ONT BAM + `work/old/` JSON),
  the age-ladder / cross-strain scripts (`age_ladder_*`, `expanded_*`, `crossstrain_telotron`
  — all read the long-read archive under `data/raw/longread/`), and `build_telogator2_ref`
  (hand-curated cap survey). `characterize_arrays`, `extend_telomeres_from_reads`,
  `detect_telomere_boundaries` are CLI-parameterised helpers (the last is driven by `telomere_boundaries`).

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
- **NCBI download concurrency.** `download_refseq` fans out with `xargs -P`; keep
  threads ≤ ~8 against NCBI to avoid rate limiting.
- **Human intronic telomeric insertions are passengers, not telotrons.** CCCTAA
  lacks GT/AG, so reverse-strand insertions can't create splice signals — they
  land inside pre-existing introns.

## External data sources
- **NCBI RefSeq** assembly summary and **Tara Oceans SMAGs v1** — URLs in
  `config.yaml` (`refseq_url`, `tara_base`).

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
