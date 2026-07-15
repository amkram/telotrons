# Telotron pipeline — Snakefile rule reference

This pipeline surveys eukaryotic genomes (NCBI RefSeq + Tara Oceans SMAGs) for **telotrons** — introns composed of tandem telomeric repeats, hypothesized to form when telomerase heals a double-strand break inside coding sequence and the inserted repeat array is then spliced out. A scan/filter core (`scan_telotrons` → `filter_final` → `dedup_telotrons` → `classify_architecture`) fans out into independent analysis arms (architecture, motif discovery, TERT deep homology, interstitial-array comparison, MSA, ortholog alignment, nucleosome features, RNA-seq expression, long-read Telogator2). Everything is orchestrated by a single Snakefile with ~56 rules; each rule shells out to one script under `scripts/`.

## 1. Genome-manifest + download

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| manifests | Build unified genome_id manifest from RefSeq assembly summary + Tara SMAGs index, optionally subset to a whitelist | RefSeq assembly summary + Tara tarball index (streamed) | work/manifests/all_genomes.tsv + refseq_euk.tsv | n/a (inline shell) |
| refseq_urls | Derive per-genome FNA + GFF download URLs from RefSeq FTP paths | work/manifests/refseq_euk.tsv | work/manifests/refseq_urls.tsv | n/a (awk) |
| tara_archives | Stream-extract only the wanted MAGs from the two ~50 GB Tara tarballs so disk doesn't blow up | work/manifests/tara_mags.tsv | data/raw/tara/.fna.done + .gff.done | n/a (curl/tar) |
| download_refseq | Fan out per-genome curls with `xargs -P 8` capped to avoid NCBI rate limits | work/manifests/refseq_urls.tsv | data/raw/refseq/.done | n/a (xargs curl) |
| canonical_motifs | Emit curated per-genome/per-group telomere-motif TSV so scan uses literature motif rather than inferring from contig ends | work/manifests/all_genomes.tsv + config `canonical_telomere_motifs` | work/manifests/canonical_motifs.tsv | n/a (inline Python `run:`) |

## 2. Core survey path

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| scan_all | Broad sieve: per-intron motif coverage, splice signals, orientation, distance-to-end for every RefSeq + Tara genome | all_genomes.tsv + refseq/tara FASTAs+GFFs | all_telotron_loci.tsv + all_introns_scanned.tsv | scripts/scan_telotrons.py |
| filter_final | Tighten scan candidates to the final admitted set using single-array + bidirectional pathways and terminal-motif matching | all_telotron_loci.tsv + all_introns_scanned.tsv | final_telotron_set.tsv + final_negative_controls.tsv | scripts/filter_final_set.py |
| dedup_telotrons | Collapse doubly-assembled exons / close paralogs by within-species flank blastn union-find, keeping longest-intron representative | final_telotron_set.tsv + refseq/tara FASTAs | final_telotron_set_dedup.tsv + dedup_log.tsv | scripts/dedup_telotrons.py |
| classify_architecture | Assign each locus to one of GT-F-AG / GT-R-AG / linker architectures and dump boundary 6-mers by class | final_telotron_set_dedup.tsv + refseq/tara FASTAs | final_telotron_set_architecture.tsv + boundary_kmers_by_architecture.tsv | scripts/classify_telotron_architecture.py |
| analyze | Boundary k-mer enrichment vs controls, length-matched distance-to-end test, architecture summary | final_telotron_set_dedup.tsv + final_species_summary.tsv | boundary_kmer_enrichment.tsv + distance_to_end.tsv | scripts/analyze_telotrons.py |
| figures | Four headline figures — counts, orientation, boundary k-mers, distance scatter | final_species_summary.tsv + analyze outputs | work/results/figures/telotron_counts.png (+ siblings) | scripts/plot_telotrons.py |
| package | Bundle final TSVs + headline figures + manifest into a single zip for sharing | PACKAGE_INPUTS list | work/results/telotron_pipeline_outputs.zip | n/a (zip) |

## 3. Extraction + alignment arm

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| extract_telotron_fasta | Per-species×architecture FASTAs + flanked [LEFT100][INTRON][RIGHT100] text extracts for downstream logos/MSAs | final_telotron_set_architecture.tsv + refseq/tara FASTAs | work/results/telotron_fasta/ + telotron_flanked/ (sentinel) | scripts/extract_telotron_fasta.py |
| build_non_telotron_controls | Sample low-telomeric-fraction introns from the same positive-species set as an intron control mirror | all_introns_scanned.tsv + final_telotron_set.tsv | work/results/non_telotron_controls.tsv | scripts/build_non_telotron_controls.py |
| extract_non_telotron_fasta | Same FASTA + flanked extraction for the control set (reuses the telotron extractor) | non_telotron_controls.tsv + refseq/tara FASTAs | work/results/non_telotron_fasta/ + non_telotron_flanked/ (sentinel) | scripts/extract_telotron_fasta.py |
| msa_telotron_regions | Per-species×architecture MAFFT MSAs of upstream\|intron\|downstream with a combined fixed-width view | final_telotron_set_architecture.tsv + refseq/tara FASTAs | work/results/msa_regions/ (sentinel) | scripts/msa_telotron_regions.py |
| blast_linkers | BLAST linker-architecture linker sequences vs own genome and vs all-species DB to find cross-species recurrence | final_telotron_set_architecture.tsv + refseq/tara FASTAs | linker_blast_hits_own_genome.tsv + linker_blast_hits_all_genomes.tsv | scripts/blast_linkers.py |

## 4. Interstitial-arrays arm

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| make_unannotated_masks | Per-genome 6-frame ATG→stop ORF mask (≥450 nt) to add to the gene/intron exclusion mask | all_species_raw_summary.tsv + refseq/tara FASTAs | work/results/masks/ (sentinel) | scripts/make_unannotated_mask.py |
| find_interstitial_arrays | Locate telomeric arrays that are non-terminal, non-genic, non-intronic (as contrast to telotrons) | manifest + refseq/tara FASTAs + masks | work/results/interstitial_arrays.tsv | scripts/find_interstitial_arrays.py |

## 5. Motif discovery arm

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| build_streme_inputs | Assemble positive FASTAs (telotrons + linkers) for STREME de novo motif discovery | telotron_fasta sentinel + all_introns_scanned.tsv + linker FASTA | streme_inputs/telotrons.fa + linkers.fa | scripts/build_streme_inputs.py |
| streme | Run STREME (--dna, w 4-12, 10 motifs, p<0.05) on each input set with STREME's own shuffled control | work/results/streme_inputs/{set}.fa | work/results/streme/{set}/streme.xml | n/a (streme CLI) |

## 6. TERT deep-homology arm

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| fetch_tert_seeds_hmms | Download apicomplexan TERT seeds + Pfam TRBD (PF12009) and RT (PF00078) HMMs used to confirm hits | (network) | tert_seeds.faa + PF12009.hmm + PF00078.hmm | scripts/fetch_tert_seeds_hmms.py |
| find_tert | Miniprot + Pfam-architecture TERT search (RBD upstream of RT), iterative re-seeding — replaces BLAST which fails on apicomplexans | tert seeds + HMMs + refseq/tara FASTAs | tert_deep_homology/confirmed_tert.tsv + all_confirmed.faa | scripts/find_tert_deep_homology.py |
| telomerase_search | Aggregate target that pulls in confirmed_tert.tsv | TERT_DEEP_HOMOLOGY_TSV | (target only) | n/a |

## 7. Telotron-gene orthologs arm

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| telotron_ortholog_align | Per telotron-host gene, miniprot-map the ortholog in each panel genome, protein-align it, then DNA-align telotron vs orthologous intron to resolve fill-vs-create | final_telotron_set_architecture.tsv + refseq/tara FASTAs | work/results/telotron_orthologs/ (sentinel) | scripts/telotron_ortholog_align.py |
| plot_telotron_ortholog_loci | Compiled per-locus figure PDF: flanking-exon protein MSA + intron DNA alignment | telotron_orthologs sentinel | work/results/figures/telotron_ortholog_loci.pdf | scripts/plot_telotron_ortholog_loci.py |
| telotron_ortholog_textdump | Per-locus text dump (unaligned/aligned × DNA/aa) grouped by within-Eimeria and outgroup panels | telotron_orthologs sentinel | telotron_orthologs_v2/locus_text/{within_eimeria,outgroup}/ | scripts/telotron_ortholog_textdump.py |
| orthologs | Aggregate target for the ortholog arm | ortho sentinel + PDF + textdump sentinel | (target only) | n/a |

## 8. Architecture / linker / interstitial-ITS analyses

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| interstitial_ortholog_textdump | Gold-standard-filtered non-genic ITS DNA-flank orthology dump (Giulotto ≥4-unit / <1mm criterion) | interstitial_arrays.tsv + refseq | work/results/interstitial_orthologs/locus_text/ | scripts/interstitial_ortholog_textdump.py |
| linker_segmentation | Decompose each telotron into array/linker segments with per-locus architecture calls | v2 locus_text | mechanism_deepdive/linker_segmentation.tsv + architecture_per_locus.tsv | scripts/linker_segmentation.py |
| cluster_linkers | Cluster segmented linkers to quantify cross-locus recurrence | linker_segmentation.tsv | linker_clusters.tsv + linker_recurrence_summary.txt | scripts/cluster_linkers.py |
| mask_telotron_arrays | Mask telomeric arrays and build a G/A/L architecture cartoon MSA over the survivor linker/flank | v2 locus_text | telotron_masked/telotron_masked.fasta + .msa.html | scripts/mask_telotron_arrays.py |
| architecture_analyses | Aggregate target chaining the four mechanism-deep-dive rules above | rules above | (target only) | n/a |

## 9. Analysis arm

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| nucleosome_inputs | Build insertion-site + flank FASTA manifests (telotron + matched control) feeding the composition/periodicity + within-gene panels | final_telotron_set_architecture.tsv + non_telotron_controls.tsv | nucleosome/manifest.tsv + control_manifest.tsv | scripts/nucleosome_inputs.py + nucleosome_control_inputs.py |
| nucleosome_features | Insertion-site composition/periodicity feature panel with telomere-masked flanks and BH-FDR | nucleosome manifests | nucleosome/nucleosome_feature_summary.png | scripts/nucleosome_features.py |
| nucleosome_withingene | Within-gene control: telotron introns vs same-gene sibling introns vs random non-host to separate gene-class confound | ARCH_TSV + non_telotron_controls.tsv | nucleosome/withingene_control.txt | scripts/nucleosome_withingene.py |
| telotron_gene_bias | Host-gene-class characterisation: intron-rich/long, host vs disjoint non-host, per-genome contrast | ARCH_TSV + refseq | work/results/figures/telotron_gene_bias.png | scripts/telotron_gene_bias.py |
| telotron_per_intron | Per-intron telotron-rate logistic model with descriptive panels and cluster caveat | ARCH_TSV + refseq | work/results/figures/telotron_per_intron.png | scripts/telotron_per_intron.py |
| length_distribution_by_arch | Length distribution by lineage×architecture (BH-FDR), single-MAG caveat; also emits per-arm burst-length subfigure | ARCH_TSV | length_distribution/summary_by_lineage_arch.tsv + length_hist PDF + per_arm_burst_length.pdf | scripts/length_distribution_by_arch_lineage.py |
| mechanism_diagrams | Capstone mechanism cartoon regenerated deterministically (no data inputs) | (none) | work/results/figures/proven_mechanism.png | scripts/plot_proven_mechanism.py |
| rnaseq_gene_coverage | Per-species SRA→samtools bedcov gene-coverage TSV — reproducible replacement for the old hand-run pipeline | config `rnaseq` + reference genome/GFF | work/results/rnaseq/{species}_gene_cov.tsv | scripts/rnaseq_gene_coverage.py |
| telotron_expr_figures | Expression-vs-telotron-presence figure (E. tenella + E. necatrix); one consolidated script | rnaseq TSVs + ARCH_TSV | work/results/figures/telotron_expression.png | scripts/telotron_expression.py |
| analysis_arm | Aggregate target for every wired downstream-analysis rule above | rules above | (target only) | n/a |

## 10. Telogator2 (long-read telomere lengths)

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| longread_data | Download every long-read FASTQ listed in data/raw/longread/manifest.tsv | manifest.tsv | data/raw/longread/{platform}/*.fastq.gz | n/a |
| download_longread_fastq | Generic per-file resumable downloader for a single FASTQ | manifest URL lookup | data/raw/longread/{platform}/{filename}.fastq.gz | n/a (wget -c) |
| telogator2 | Aggregate target that runs Telogator2 on every nanopore run in the manifest | per-accession tlens outputs | (target only) | n/a |
| telogator2_one | Per-accession Telogator2 allele-specific telomere-length / TVR call against a per-species Eimeria subtelomere reference | nanopore FASTQ + per-species telogator2 ref | telogator2/{species}_{accession}/tlens_by_allele.tsv | n/a (telogator2.py from external repo) |

Note: `build_telogator2_ref` is documented as a hand-run helper — no Snakefile rule exists because the per-species `*.telogator2.fasta` are hand-built from a curated cap-survey TSV.

## 11. Figures / plots

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| terminal_motif_figures | Per-species rolling hit-density plots of the terminal telomere motif with telotron loci overlaid | all_species_raw_summary.tsv + all_telotron_loci.tsv | work/results/figures/terminal_motif/ (sentinel) | scripts/plot_terminal_motif_density.py |
| plot_boundary_kmers_by_arch | Per-species top boundary 6-mers split by architecture (donor / acceptor / linker edges) | boundary_kmers_by_architecture.tsv + final_species_summary.tsv | figures/boundary_kmers_by_arch/ (sentinel) | scripts/plot_boundary_kmers_by_arch.py |
| plot_interstitial_boundary_kmers | Per-species 5'/3' boundary 6- and 12-mers of interstitial arrays, optionally split by 5'/3' orientation category | interstitial_arrays.tsv | figures/interstitial_boundary_kmers{cat}/ (sentinel) | scripts/plot_interstitial_boundary_kmers.py |
| plot_array_length_distribution | Per-species histograms comparing ITS, telotron, and non-telotron intron lengths (unfiltered vs ≥40 bp) | interstitial + final_telotron_set + non_telotron_controls | figures/array_length_distribution{suffix}.png | scripts/plot_array_length_distribution.py |
| plot_pipeline_stages | Pipeline-stage summary figures: manifest, scan, filter funnel, distance, architecture | all_genomes.tsv + raw/final summaries + arch loci | figures/pipeline_stages/ (sentinel) | scripts/plot_pipeline_stages.py |
| plot_splice_signal_logos | ±10 bp GT donor / AG acceptor sequence logos for telotron (optionally by architecture) or non-telotron control | flanked-text sentinel (+ arch table when filtering) | figures/{set}_splice_logos{cat}/ (sentinel) | scripts/plot_splice_signal_logos.py |
| telomere_boundaries | Motif sliding-window telomere-boundary detection on assembled long-read contigs — per-arm boundary table | assembled contigs + config `telomere_boundaries` | telomere_boundaries/per_arm.tsv + telomeres.bed | scripts/detect_telomere_boundaries.py |
| pipeline_report | Aggregate the whole pipeline into a single self-contained HTML report with inline base64 figures and colored MSAs | PACKAGE_INPUTS + fasta/MSA/figure sentinels | work/results/pipeline_report.html | scripts/build_pipeline_report.py |

## 12. Aggregates

| Rule | Purpose | Key inputs | Key outputs | Script |
|---|---|---|---|---|
| all | Default target — core survey through `package` + `pipeline_report` | zip + report | (target only) | n/a |

## Running it

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
