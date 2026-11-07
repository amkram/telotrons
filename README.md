# Telotrons: analysis code

Code for *Telotrons: spliceosomal introns composed of telomeric repeats* (Kramer and Corbett-Detig). It runs all analyses and builds the tables and figures.

## Inputs

All under `config/`. Genomes, annotations and reads are not included. Put assemblies and annotations under `data/raw/refseq/<accession>/`, `data/raw/genbank/<accession>/` or `data/raw/tara/`. `lib/_common.py` finds the FASTA and GFF by accession.

| File | Contents | Stage |
|---|---|---|
| `genomes.tsv` | 25 bearers and 25 controls: accession, lineage, motif, pairing | all |
| `strand_corrections.tsv` | 216 loci whose strand is set to an earlier curated call | 1 |
| `rna_runs.tsv` | SRA/ENA runs per genome | 4 |
| `comparator_panels.tsv`, `comparator_motifs.tsv` | comparator genomes and motifs; genomes go under `data/raw` | 5 |
| `orthology_review.tsv` | 14 telotron-comparator pairs set by manual review (12 UNRESOLVED_MAPPING, 2 AMBIGUOUS) | 5 |
| `alternative_requests.json` | alternative comparator introns to measure (format in `scripts/06_exact_position/extract_alternatives.py`) | 6 |
| `its_donors.tsv` | genomes whose proteins are projected onto each of the 50 to mask genes; the 21 that are not bearers or controls go under `data/raw` | 7 |
| `single_motif_genomes.tsv` | the 50 genomes with FASTA and GFF paths | 10 |
| `figures/` | tree topologies (row order), dated backbone and crown ages, exon phases, Figure 2B examples, MAG BUSCO scores | 8–9 |
| `supplement/` | display names, sources and intron counts for the supplementary tables | 8 |
| `screened_genomes.tsv` | the 10,310 assemblies of the initial survey and their lineage motifs; a record, not read by any stage | – |

## Run

    conda env create -f environment.yml && conda activate telotrons
    export TMPDIR=/path/to/job/scratch   # never /tmp
    workflow/run_all.sh

`run_all.sh` writes the results package to `TELO_PACKAGE` (default `work/package`: `tables/`, `figures/`, `results/`, `supplementary/`, `qa/`, `inputs/`) and stage 5 work to `work/05_orthology`. Each stage reads `config/`, `data/` and earlier stages' output under `work/`.

| Stage | Scripts | Output |
|---|---|---|
| 1 | `01_scan_introns.py` per genome, then `--collect` | introns >50% one lineage motif, with architecture, deduplication and strand (`standard` >85%, `additional` >50% to 85%); `work/catalogue/catalogue.tsv.gz` |
| 2 | `02_control_screen.py` | highest intron coverage in each control |
| 3 | `03_tsd.py` | putative TSDs at primary loci |
| 4 | `04_rna_support.py`, `04b_ordinary_rna_baseline.py` per genome, then `04c_collect_rna.py` | junction reads (minimap2 `splice:sr`, MAPQ ≥20, ≥10 bases each side) per locus, and the ordinary-intron baseline |
| 5 | `05_orthology/` | miniprot mapping, reciprocal best hits (≥40% identity, ≥12 matches per 20-residue flank), absence correction, manual review, closest comparator per locus |
| 6 | `06_exact_position/` | position and phase of ordinary comparator introns |
| 7 | `07_its/` | ITS tracts (≥4 units, ≤1 substitution per adjacent pair) outside annotated and projected genes, densities, exact Wilcoxon test on seven lineage pairs |
| 8 | `08_tables/` | evidence joined onto the catalogue and the supplementary tables (`.docx` via `node build_supplement_docx.js` after `npm install docx`) |
| 9 | `09_figures/` | Figures 1–3, S1, S2 and an unnumbered architecture-by-tier figure (`prepare_architecture.py`, then `render_all.py`) |
| 10 | `10_single_motif/` | every intron of the 50 genomes scored against every repeat motif of period 1–8 (`scan.py`, `deduplicate.py`), then Figure S3 (`09_figures/figure_s3_single_motif.py`) |

Figure 1 and S2 use TimeTree ages for nodes between lineages (`config/figures/backbone_dated.nwk`) and draw each lineage as a polytomy at its published crown age (`config/figures/crown_ages.tsv`, with sources).
