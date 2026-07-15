# Telotron pipeline — Snakefile rule reference

Surveys eukaryotic genomes (**all annotated RefSeq + all annotated GenBank + all Tara Oceans SMAGs**) for **telotrons** — introns composed of tandem telomeric repeats, hypothesized to form when telomerase heals a double-strand break inside coding sequence and the inserted repeat array is then spliced out. A scan/filter core (`scan_all` → `filter_final` → `dedup_telotrons` → `classify_architecture`) feeds `analyze` (per-locus stats) and `confident_species` (species-level bearer set). Downstream analyses (TERT, orthologs, linker/mask, nucleosome, gene class, expression, length) key off the confident-species set. One Snakefile, **25 rules**; every rule shells out to one script under `scripts/`.

## 1. Download

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| manifests | Union RefSeq (curated GCF_) + GenBank (annotated eukaryotes w/o paired GCF_) + Tara SMAGs into one manifest with a `source` column | RefSeq + GenBank assembly summaries + Tara tarball index | work/manifests/all_genomes.tsv + refseq_euk + genbank_euk + tara_mags | n/a |
| tara_archives | Stream-extract only the wanted MAGs from the two ~50 GB Tara tarballs | work/manifests/tara_mags.tsv | data/raw/tara/.fna.done + .gff.done | n/a |
| download_assemblies | Derive per-genome URLs inline + fan out curls (`xargs -P 8`); RefSeq → data/raw/refseq/, GenBank → data/raw/genbank/ | refseq_euk.tsv + genbank_euk.tsv | data/raw/refseq/.done + data/raw/genbank/.done | n/a |

## 2. Core survey → confident set

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| scan_all | Per-intron motif coverage, splice signals, orientation, distance-to-end. Builds canonical per-genome motif TSV inline from config | all_genomes.tsv + refseq/genbank/tara FASTAs+GFFs | all_telotron_loci.tsv + all_introns_scanned.tsv + canonical_motifs.tsv | scripts/scan_telotrons.py |
| filter_final | Two admission pathways (single-array `min_repeat_frac`=0.85; bidirectional `bidir_min_repeat_frac`=0.40 + `bidir_min_hits`=3) plus terminal-motif matching | all_telotron_loci.tsv + all_introns_scanned.tsv | final_telotron_set.tsv + final_negative_controls.tsv | scripts/filter_final_set.py |
| dedup_telotrons | Collapse doubly-assembled exons / close paralogs by within-species flank blastn union-find | final_telotron_set.tsv + refseq/genbank/tara FASTAs | final_telotron_set_dedup.tsv + dedup_log.tsv | scripts/dedup_telotrons.py |
| classify_architecture | Assign each locus to GT-F-AG / GT-R-AG / linker architecture and dump boundary 6-mers by class | final_telotron_set_dedup.tsv + refseq/genbank/tara FASTAs | final_telotron_set_architecture.tsv + boundary_kmers_by_architecture.tsv | scripts/classify_telotron_architecture.py |
| analyze | Boundary k-mer enrichment vs controls, length-matched distance-to-end test, architecture summary | final_telotron_set_dedup.tsv + final_species_summary.tsv | boundary_kmer_enrichment.tsv + distance_to_end.tsv | scripts/analyze_telotrons.py |
| confident_species | **Paper's central bearer set** — species admitted when ≥3 telotrons pass filter OR ≥1 bidirectional architecture. Every downstream analysis keys off this file | final_telotron_set_architecture.tsv + all_genomes.tsv | work/results/confident_species.tsv | scripts/confident_species.py |
| package | Bundle final TSVs + confident-species set + manifest into a single deliverable zip | PACKAGE_INPUTS list | work/results/telotron_pipeline_outputs.zip | n/a (zip) |

## 3. Extraction + control

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| extract_fasta | Per-species×architecture FASTAs + flanked-text extracts. Wildcarded on `{set}` (telotron | non_telotron); one rule for both extract arms | final_telotron_set_architecture.tsv OR non_telotron_controls.tsv + refseq/genbank/tara FASTAs | work/results/{set}_fasta/.done | scripts/extract_telotron_fasta.py |
| build_non_telotron_controls | Sample low-telomeric-fraction introns from the same positive species as intron control mirror | all_introns_scanned.tsv + final_telotron_set.tsv | work/results/non_telotron_controls.tsv | scripts/build_non_telotron_controls.py |

## 4. Interstitial contrast

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| find_interstitial_arrays | Locate non-terminal, non-genic, non-intronic telomeric arrays; 6-frame ORF mask built inline | all_species_raw_summary.tsv + refseq/genbank/tara FASTAs | work/results/interstitial_arrays.tsv | scripts/make_unannotated_mask.py + scripts/find_interstitial_arrays.py |

## 5. TERT deep homology

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| fetch_tert_seeds_hmms | Download apicomplexan TERT seeds + Pfam TRBD (PF12009) + RT (PF00078) HMMs | (network) | tert_seeds.faa + PF12009.hmm + PF00078.hmm | scripts/fetch_tert_seeds_hmms.py |
| find_tert | Miniprot + Pfam-architecture TERT search (RBD upstream of RT), iterative re-seeding | tert seeds + HMMs + refseq/genbank/tara FASTAs | tert_deep_homology/confirmed_tert.tsv + all_confirmed.faa | scripts/find_tert_deep_homology.py |

## 6. Telotron-gene orthologs

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| telotron_orthologs | Miniprot-map ortholog in each panel genome, protein-align, DNA-align telotron vs orthologous intron to resolve fill-vs-create; emit compiled per-locus PDF (was 2 rules) | final_telotron_set_architecture.tsv + refseq/genbank/tara FASTAs | work/results/telotron_orthologs/ (sentinel) + telotron_ortholog_loci.pdf | scripts/telotron_ortholog_align.py + scripts/plot_telotron_ortholog_loci.py |

## 7. Architecture / linker analyses

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| linker_analysis | Segment each telotron into arrays + linkers, then cluster linker sequences by 7-mer Jaccard to quantify cross-locus recurrence (was 2 rules) | v3 locus_text | linker_segmentation.tsv + linker_clusters.tsv + recurrence_summary.txt | scripts/linker_segmentation.py + scripts/cluster_linkers.py |
| mask_telotron_arrays | Mask telomeric arrays, build G/A/L architecture cartoon MSA over surviving linker/flank | v3 locus_text | telotron_masked/telotron_masked.fasta + .msa.html | scripts/mask_telotron_arrays.py |

## 8. Analysis arm

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| telotron_gene_class | Host-gene-class dissection (host vs disjoint non-host intron count + gene length; per-intron rate vs gene properties + 5' preference). Iterates confident_species dynamically | ARCH_TSV + confident_species.tsv + refseq/genbank GFFs | work/results/figures/telotron_gene_class.png | scripts/telotron_gene_class.py |
| length_distribution_by_arch | Length distribution by lineage × architecture (BH-FDR); also emits per-arm burst-length figure | ARCH_TSV | work/results/length_distribution/summary_by_lineage_arch.tsv + per_arm_burst_length.pdf | scripts/length_distribution_by_arch_lineage.py |
| nucleosome_analysis | Insertion-site composition / 10-bp WW periodicity / CpG panel (telomere-masked, BH-FDR) + within-gene sibling-intron control. One rule (was 2). NuPoP occupancy retired | ARCH_TSV + non_telotron_controls.tsv + refseq/genbank/tara FASTAs | work/results/nucleosome/nucleosome_feature_summary.png + withingene_control.txt | scripts/nucleosome_inputs.py + scripts/nucleosome_control_inputs.py + scripts/nucleosome_features.py + scripts/nucleosome_withingene.py |
| rnaseq_gene_coverage | Per-species SRA→samtools bedcov gene-coverage TSV (config-driven) | config `rnaseq` + genome/GFF | work/results/rnaseq/{species}_gene_cov.tsv | scripts/rnaseq_gene_coverage.py |
| telotron_expr_figures | Expression vs telotron presence for E. tenella + E. necatrix (size-controlled OLS + per-intron rate vs expression quintile) | rnaseq TSVs + ARCH_TSV | work/results/figures/telotron_expression.png | scripts/telotron_expression.py |
| analysis_arm | Aggregate of every wired downstream figure (invoke with `snakemake analysis_arm`) | rules above | (target only) | n/a |

## Running it

```bash
# Managed deps (builds envs/telotrons.yaml).
snakemake --use-conda -j 16

# Or with tools on PATH: gt (GenomeTools), seqkit, bedtools, samtools,
# blast+, mafft, miniprot, hmmer, python + pandas/scipy/matplotlib.
snakemake -j 16

snakemake -n                                     # dry run
snakemake --use-conda -j 8 --forcerun scan_all   # single-stage rerun

# Test on a subset: set `accessions: [...]` in config.yaml or:
snakemake -j 8 --config 'accessions=["GCF_000499545.2","GCF_000499605.1"]'
```
