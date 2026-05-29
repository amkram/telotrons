#!/bin/bash
# Full pipeline run on the 22-species test set.
# Reads manifests from manifests.test_20260522_0450/, writes to results_test/.
# Safe to run alongside the ongoing full survey (separate output dirs).

set -euo pipefail

cd /scratch1/alex/telotrons

ENV_BIN=/home/alex/.snakemake-envs/telotrons/7aa0f2645427a08bdf054e56b1233123_/bin
export PATH="$ENV_BIN:$PATH"
PYTHON="$ENV_BIN/python"
MOTIFS=TTAGG,TCAGG,TTGGG,TTAGGG,TTAGGC,TTGGGG,TTTAGGG,TTTTAGGG,TTTTGGGG,TATTAGGG,TTATTGGG,TTATTAGGG,TTACTTGGG,TTATTGGGG,TTTTTTAGGG,TTTATTAGGG,TTAGGTTGGGG,CTCGGTTATGGG
SCAN_THREADS=32
DOWNLOAD_PARALLEL=8

MANIFEST_DIR=manifests.test_20260522_0450
OUTDIR=results_test
LOGDIR=logs/$(date +%Y%m%d_%H%M%S)_test_survey
mkdir -p "$LOGDIR" "$OUTDIR" raw/refseq raw/tara

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOGDIR/run.log"; }

log "=== Test survey: 22-species set ==="
log "manifest: $MANIFEST_DIR/all_genomes.tsv"

# ─────────────────────────────────────────────────────────────────────────
log "=== 1. Tara (already extracted) ==="
if [ ! -f raw/tara/.fna.done ] || [ ! -f raw/tara/.gff.done ]; then
    log "Tara not extracted — extracting now"
    TARA_BASE="https://www.genoscope.cns.fr/tara/localdata/data/SMAGs-v1"
    (curl -L --fail -s "${TARA_BASE}/SMAGs_contigs_individual.fna.tar.gz" \
        | tar -xzf - -C raw/tara) &
    PID_F=$!
    (curl -L --fail -s "${TARA_BASE}/SMAGs_v1_individual.gff.tar.gz" \
        | tar -xzf - -C raw/tara) &
    PID_G=$!
    wait $PID_F $PID_G
    touch raw/tara/.fna.done raw/tara/.gff.done
    log "tara extraction done"
else
    log "tara already extracted, skipping"
fi

# ─────────────────────────────────────────────────────────────────────────
log "=== 2. Download RefSeq test genomes ==="
REFSEQ_IDS=$(awk -F'\t' 'NR>1 && $3!="tara" {print $1}' "$MANIFEST_DIR/all_genomes.tsv")
for gid in $REFSEQ_IDS; do
    if ls raw/refseq/"$gid"/*.fna.gz &>/dev/null 2>&1; then
        log "  $gid already present, skipping"
        continue
    fi
    log "  downloading $gid"
    row=$(grep "^$gid	" "$MANIFEST_DIR/refseq_urls.tsv" || true)
    if [ -z "$row" ]; then
        log "  WARNING: no URL found for $gid"
        continue
    fi
    fna_url=$(echo "$row" | awk -F'\t' '{print $2}')
    gff_url=$(echo "$row" | awk -F'\t' '{print $3}')
    mkdir -p "raw/refseq/$gid"
    curl -L --fail --retry 5 --retry-delay 3 --retry-all-errors -s \
        "$fna_url" -o "raw/refseq/$gid/$(basename "$fna_url")" &
    curl -L --fail --retry 5 --retry-delay 3 --retry-all-errors -s \
        "$gff_url" -o "raw/refseq/$gid/$(basename "$gff_url")" &
done
wait
log "RefSeq downloads done"

# ─────────────────────────────────────────────────────────────────────────
log "=== 3. Scan all test genomes ==="
if [ ! -s "$OUTDIR/all_telotron_loci.tsv" ]; then
    "$PYTHON" scripts/scan_telotrons.py \
        --manifest "$MANIFEST_DIR/all_genomes.tsv" \
        --canonical-motifs "$MANIFEST_DIR/canonical_motifs.tsv" \
        --refseq-dir raw/refseq --tara-dir raw/tara \
        --motifs "$MOTIFS" \
        --threads "$SCAN_THREADS" \
        --loci "$OUTDIR/all_telotron_loci.tsv" \
        --introns "$OUTDIR/all_introns_scanned.tsv" \
        --summary "$OUTDIR/all_species_raw_summary.tsv" \
        2>&1 | tee -a "$LOGDIR/scan.log"
    log "scan done"
else
    log "scan outputs exist, skipping"
fi

# ─────────────────────────────────────────────────────────────────────────
log "=== 4. Filter ==="
"$PYTHON" scripts/filter_final_set.py \
    --loci "$OUTDIR/all_telotron_loci.tsv" \
    --introns "$OUTDIR/all_introns_scanned.tsv" \
    --summary "$OUTDIR/all_species_raw_summary.tsv" \
    --min-repeat-frac 0.85 --bidir-min-repeat-frac 0.40 --bidir-min-hits 3 \
    --collapse-unique-loci --require-terminal-motif-match \
    --final "$OUTDIR/final_telotron_set.tsv" \
    --species "$OUTDIR/final_species_summary.tsv" \
    --negatives "$OUTDIR/final_negative_controls.tsv" \
    2>&1 | tee -a "$LOGDIR/filter.log"
log "filter done"

# ─────────────────────────────────────────────────────────────────────────
log "=== 4b. Deduplicate (flank-similarity collapse) ==="
"$PYTHON" scripts/dedup_telotrons.py \
    --final "$OUTDIR/final_telotron_set.tsv" \
    --refseq-dir raw/refseq --tara-dir raw/tara \
    --out-final "$OUTDIR/final_telotron_set_dedup.tsv" \
    --out-log "$OUTDIR/dedup_log.tsv" \
    --threads "$SCAN_THREADS" \
    2>&1 | tee -a "$LOGDIR/dedup.log"
log "dedup done"

# Replace the working file in-place so all downstream stages see the dedup set,
# but keep the raw filter output beside it for provenance.
mv "$OUTDIR/final_telotron_set.tsv" "$OUTDIR/final_telotron_set.prededup.tsv"
cp "$OUTDIR/final_telotron_set_dedup.tsv" "$OUTDIR/final_telotron_set.tsv"

# ─────────────────────────────────────────────────────────────────────────
log "=== 5. Analyze ==="
"$PYTHON" scripts/analyze_telotrons.py \
    --final "$OUTDIR/final_telotron_set.tsv" \
    --introns "$OUTDIR/all_introns_scanned.tsv" \
    --species "$OUTDIR/final_species_summary.tsv" \
    --boundary-kmers "$OUTDIR/boundary_kmer_enrichment.tsv" \
    --distance "$OUTDIR/distance_to_end.tsv" \
    --architecture "$OUTDIR/architecture_summary.tsv" \
    --threads "$SCAN_THREADS" \
    2>&1 | tee -a "$LOGDIR/analyze.log"
log "analyze done"

# ─────────────────────────────────────────────────────────────────────────
log "=== 6. Architecture classification ==="
"$PYTHON" scripts/classify_telotron_architecture.py \
    --final "$OUTDIR/final_telotron_set.tsv" \
    --refseq-dir raw/refseq --tara-dir raw/tara \
    --out-loci "$OUTDIR/final_telotron_set_architecture.tsv" \
    --out-kmers "$OUTDIR/boundary_kmers_by_architecture.tsv" \
    2>&1 | tee -a "$LOGDIR/architecture.log"
log "architecture done"

# ─────────────────────────────────────────────────────────────────────────
log "=== 7. FASTA extracts ==="
"$PYTHON" scripts/extract_telotron_fasta.py \
    --final "$OUTDIR/final_telotron_set_architecture.tsv" \
    --refseq-dir raw/refseq --tara-dir raw/tara \
    --fasta-dir "$OUTDIR/telotron_fasta" \
    --flanked-dir "$OUTDIR/telotron_flanked" \
    2>&1 | tee -a "$LOGDIR/extract_fasta.log"
log "fasta extracts done"

# ─────────────────────────────────────────────────────────────────────────
log "=== 7b. BLAST linkers (own genome + combined db) ==="
"$PYTHON" scripts/blast_linkers.py \
    --arch-tsv "$OUTDIR/final_telotron_set_architecture.tsv" \
    --refseq-dir raw/refseq --tara-dir raw/tara \
    --workdir "$OUTDIR/blast_linkers" \
    --out-own "$OUTDIR/linker_blast_hits_own_genome.tsv" \
    --out-all "$OUTDIR/linker_blast_hits_all_genomes.tsv" \
    --threads "$SCAN_THREADS" \
    2>&1 | tee -a "$LOGDIR/blast_linkers.log"
log "blast linkers done"

# ─────────────────────────────────────────────────────────────────────────
log "=== 7c. MAFFT MSAs per species × architecture × region ==="
"$PYTHON" scripts/msa_telotron_regions.py \
    --arch-tsv "$OUTDIR/final_telotron_set_architecture.tsv" \
    --refseq-dir raw/refseq --tara-dir raw/tara \
    --outdir "$OUTDIR/msa_regions" \
    --threads 4 \
    2>&1 | tee -a "$LOGDIR/msa_regions.log"
log "MSAs done"

# ─────────────────────────────────────────────────────────────────────────
log "=== 7d. Interstitial telomeric arrays (non-genic, non-intronic) ==="
"$PYTHON" scripts/find_interstitial_arrays.py \
    --manifest "$MANIFEST_DIR/all_genomes.tsv" \
    --refseq-dir raw/refseq --tara-dir raw/tara \
    --motifs "$MOTIFS" \
    --out "$OUTDIR/interstitial_arrays.tsv" \
    --threads "$SCAN_THREADS" \
    2>&1 | tee -a "$LOGDIR/interstitial.log"

"$PYTHON" scripts/plot_interstitial_boundary_kmers.py \
    --arrays "$OUTDIR/interstitial_arrays.tsv" \
    --outdir "$OUTDIR/figures/interstitial_boundary_kmers" \
    2>&1 | tee -a "$LOGDIR/plot_interstitial.log"
log "interstitial done"

# ─────────────────────────────────────────────────────────────────────────
log "=== 8. Figures ==="
mkdir -p "$OUTDIR/figures"

"$PYTHON" scripts/plot_telotrons.py \
    --species "$OUTDIR/final_species_summary.tsv" \
    --final "$OUTDIR/final_telotron_set.tsv" \
    --kmers "$OUTDIR/boundary_kmer_enrichment.tsv" \
    --distance "$OUTDIR/distance_to_end.tsv" \
    --architecture "$OUTDIR/architecture_summary.tsv" \
    --outdir "$OUTDIR/figures" 2>&1 | tee -a "$LOGDIR/plot_telotrons.log"

"$PYTHON" scripts/plot_boundary_kmers_by_arch.py \
    --kmers "$OUTDIR/boundary_kmers_by_architecture.tsv" \
    --species "$OUTDIR/final_species_summary.tsv" \
    --outdir "$OUTDIR/figures/boundary_kmers_by_arch" \
    2>&1 | tee -a "$LOGDIR/plot_boundary.log"

"$PYTHON" scripts/plot_pipeline_stages.py \
    --manifests "$MANIFEST_DIR/all_genomes.tsv" \
    --raw-summary "$OUTDIR/all_species_raw_summary.tsv" \
    --species-final "$OUTDIR/final_species_summary.tsv" \
    --distance "$OUTDIR/distance_to_end.tsv" \
    --architecture-loci "$OUTDIR/final_telotron_set_architecture.tsv" \
    --outdir "$OUTDIR/figures/pipeline_stages" \
    2>&1 | tee -a "$LOGDIR/plot_stages.log"

log "=== ALL DONE ==="
log "results_test/ size: $(du -sh "$OUTDIR" | awk '{print $1}')"
