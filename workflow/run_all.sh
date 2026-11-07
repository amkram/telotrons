#!/usr/bin/env bash
# Every stage in order, from config/ and data/ to the results package. Per-genome loops can run as
# cluster array tasks instead (one genome per task). Each stage needs only the outputs of earlier
# stages under work/.
set -euo pipefail
cd "$(dirname "$0")/.."

# job scratch from the caller, never /tmp. Plotting caches go there too.
export TMPDIR="${TMPDIR:?set TMPDIR to a job-specific scratch directory}"
export TMP="$TMPDIR" TEMP="$TMPDIR" MPLCONFIGDIR="$TMPDIR/matplotlib"
# the output package (tables/ figures/ results/ supplementary/ qa/ inputs/) and the stage 5
# work directory
export TELO_PACKAGE="${TELO_PACKAGE:-$PWD/work/package}"
export TELO_WORKDIR="$PWD/work/05_orthology"

# genome ids of config/genomes.tsv with a given role (bearer, control or all)
genomes() {
  awk -F'\t' -v role="$1" 'NR > 1 && (role == "all" || $3 == role) {print $1}' config/genomes.tsv
}
# genomes with RNA-seq runs
rna_genomes() { awk -F'\t' 'NR > 1 {print $1}' config/rna_runs.tsv; }

# 1. admission under one lineage motif, per genome, then one catalogue
for g in $(genomes all); do
  python scripts/01_scan_introns.py "$g" --threads 8
done
python scripts/01_scan_introns.py --collect

# 2. control screen
python scripts/02_control_screen.py $(genomes control)

# 3. putative TSDs at the primary loci
python scripts/03_tsd.py --catalogue work/catalogue/catalogue.tsv.gz

# 4. RNA junction reads at the catalogue loci and at ordinary introns, per genome with runs, then
#    work/04_rna/RNA_loci.tsv (04_rna_support.py --bam reuses an existing alignment)
for g in $(rna_genomes); do
  python scripts/04_rna_support.py --genome-id "$g"
  python scripts/04b_ordinary_rna_baseline.py --genome-id "$g"
done
python scripts/04c_collect_rna.py

# 5. orthology: inputs, then per tier and bearer (array index into inputs/genomes.json) host
#    mapping, alignment, recall and the annotation-supported absence correction, then the curated
#    pair review (config/orthology_review.tsv), the correction summary and the closest comparator of
#    every locus
python scripts/05_orthology/prepare_inputs.py
n_bearers=$(genomes bearer | wc -l)
for tier in primary additional; do
  for i in $(seq 0 $((n_bearers - 1))); do
    TELO_TASK_INDEX=$i python scripts/05_orthology/host.py "$tier"
    TELO_TASK_INDEX=$i python scripts/05_orthology/orthology.py "$tier"
    TELO_TASK_INDEX=$i python scripts/05_orthology/recall.py "$tier"
    TELO_TASK_INDEX=$i python scripts/05_orthology/absence_correction.py "$tier"
  done
done
python scripts/05_orthology/apply_review.py
python scripts/05_orthology/summarize_correction.py
python scripts/05_orthology/closest.py

# 6. exact position and phase of ordinary comparator introns. extract_alternatives.py reads the
#    hand-written config/alternative_requests.json
for step in extract_evidence extract_annotations extract_predictions extract_alternatives \
  analyze_models analyze_predictions measure_alternatives summarize combine; do
  python "scripts/06_exact_position/$step.py"
done

# 7. ITS tracts and projected gene masks (donors from config/its_donors.tsv) for every genome,
#    then densities, pairs and the exact test
for g in $(genomes all); do
  python scripts/07_its/scan.py --genome-id "$g"
  python scripts/07_its/projection.py --genome-id "$g"
done
python scripts/07_its/collect.py
python scripts/07_its/exact_stats.py --pairs work/07_its/pairs.tsv \
  --out work/07_its/exact_wilcoxon.json

# 8. evidence joined onto the catalogue, source tables, supplementary tables
python scripts/08_tables/aggregate.py
python scripts/08_tables/prepare_final.py
python scripts/08_tables/build_table_s2.py
python scripts/08_tables/build_tables_s1_s4.py
python scripts/08_tables/build_supplement.py
#    Word files (needs `npm install docx`): node scripts/08_tables/build_supplement_docx.js

# 9. figures
python scripts/09_figures/prepare_architecture.py
python scripts/09_figures/render_all.py

# 10. introns made of any one repeat motif in the 50 genomes, then Figure S3
mkdir -p work/10_single_motif
cc -O3 -std=c11 -Wall -Wextra -Werror -shared -fPIC scripts/10_single_motif/repeats.c \
  -o work/10_single_motif/librepeats.so
for g in $(genomes all); do
  python scripts/10_single_motif/scan.py --genome "$g" --manifest config/single_motif_genomes.tsv \
    --library work/10_single_motif/librepeats.so --out "work/10_single_motif/$g"
  python scripts/10_single_motif/deduplicate.py --directory "work/10_single_motif/$g"
done
python scripts/09_figures/figure_s3_single_motif.py
