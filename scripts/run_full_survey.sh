#!/bin/bash
# End-to-end pan-eukaryotic survey with per-chunk download + scan + cleanup.
# Designed for a multi-day background run with peak raw disk capped at ~100 GB.
#
# Stages:
#   1. (assumes) work/manifests/all_genomes.tsv + work/manifests/refseq_urls.tsv exist
#   2. Tara: extract once (~60 GB), KEEP throughout (small enough)
#   3. Tara scan (one pass)
#   4. RefSeq: stream in chunks of 32 — download, scan, rm -rf chunk
#   5. Concatenate per-chunk TSVs into the canonical work/results/all_*.tsv
#   6. Re-download FASTAs for positive species only, run downstream analysis
#
# Resumable: each chunk writes its own TSVs; existing chunk outputs are skipped.

set -euo pipefail

cd /scratch1/alex/telotrons

source scripts/_survey_env.sh   # PATH/PYTHON, MOTIFS, SCAN_THREADS, DOWNLOAD_PARALLEL, log()

CHUNK_SIZE=32
MIN_FREE_GB=150
LOGDIR=logs/$(date +%Y%m%d_%H%M%S)_full_survey
mkdir -p "$LOGDIR" work/results/_chunks raw/refseq raw/tara

free_gb() { df -BG /scratch1 | tail -1 | awk '{gsub("G","",$4); print $4+0}'; }

guard_disk() {
    local fg=$(free_gb)
    if [ "$fg" -lt "$MIN_FREE_GB" ]; then
        log "ABORT: only ${fg}G free, below ${MIN_FREE_GB}G threshold"
        return 1
    fi
}

# ─────────────────────────────────────────────────────────────────────────
log "=== 1. Tara extraction (one-time, ~60GB kept) ==="
if [ ! -f raw/tara/.fna.done ] || [ ! -f raw/tara/.gff.done ]; then
    guard_disk
    TARA_BASE="https://www.genoscope.cns.fr/tara/localdata/data/SMAGs-v1"
    log "extracting full Tara SMAGs v1 (all MAGs)"
    (curl -L --fail -s "${TARA_BASE}/SMAGs_contigs_individual.fna.tar.gz" \
        | tar -xzf - -C raw/tara) &
    PID_F=$!
    (curl -L --fail -s "${TARA_BASE}/SMAGs_v1_individual.gff.tar.gz" \
        | tar -xzf - -C raw/tara) &
    PID_G=$!
    wait $PID_F $PID_G
    touch raw/tara/.fna.done raw/tara/.gff.done
    log "tara extraction done, $(du -sh raw/tara | awk '{print $1}') on disk"
else
    log "tara already extracted, skipping"
fi

# ─────────────────────────────────────────────────────────────────────────
log "=== 2. Scan Tara MAGs ==="
if [ ! -s work/results/_chunks/tara_loci.tsv ]; then
    awk -F'\t' 'NR==1 || $3=="tara"' work/manifests/all_genomes.tsv > "$LOGDIR/tara_manifest.tsv"
    "$PYTHON" scripts/scan_telotrons.py \
        --manifest "$LOGDIR/tara_manifest.tsv" \
        --canonical-motifs work/manifests/canonical_motifs.tsv \
        --refseq-dir raw/refseq --tara-dir raw/tara \
        --motifs "$MOTIFS" \
        --threads "$SCAN_THREADS" \
        --loci work/results/_chunks/tara_loci.tsv \
        --introns work/results/_chunks/tara_introns.tsv \
        --summary work/results/_chunks/tara_summary.tsv \
        2>&1 | tee -a "$LOGDIR/tara_scan.log"
    log "tara scan done"
else
    log "tara scan already complete, skipping"
fi

# ─────────────────────────────────────────────────────────────────────────
log "=== 3. RefSeq streaming scan (chunks of $CHUNK_SIZE) ==="

# Build the refseq-only URL+manifest pair, sorted so chunks are reproducible
awk -F'\t' 'NR>1 && $3!="tara"' work/manifests/all_genomes.tsv | sort -k1,1 > "$LOGDIR/refseq_manifest_rows.tsv"
sort -k1,1 work/manifests/refseq_urls.tsv > "$LOGDIR/refseq_urls_sorted.tsv"

N_REFSEQ=$(wc -l < "$LOGDIR/refseq_urls_sorted.tsv")
N_CHUNKS=$(( (N_REFSEQ + CHUNK_SIZE - 1) / CHUNK_SIZE ))
log "RefSeq genomes: $N_REFSEQ, chunks: $N_CHUNKS"

split -l "$CHUNK_SIZE" -d -a 5 "$LOGDIR/refseq_urls_sorted.tsv" "$LOGDIR/chunk_url_"

for chunk_url_file in "$LOGDIR"/chunk_url_*; do
    CHUNK_ID=$(basename "$chunk_url_file" | sed 's/chunk_url_//')
    LOCI_OUT="work/results/_chunks/refseq_loci_${CHUNK_ID}.tsv"
    INTRONS_OUT="work/results/_chunks/refseq_introns_${CHUNK_ID}.tsv"
    SUMMARY_OUT="work/results/_chunks/refseq_summary_${CHUNK_ID}.tsv"

    if [ -s "$SUMMARY_OUT" ]; then
        # already processed
        continue
    fi

    guard_disk || { log "disk guard tripped at chunk $CHUNK_ID"; exit 1; }
    N=$(wc -l < "$chunk_url_file")
    log "chunk $CHUNK_ID: downloading $N genomes (free $(free_gb)G)"

    # Download
    awk -F'\t' '{print $1"\t"$2"\n"$1"\t"$3}' "$chunk_url_file" \
      | xargs -P "$DOWNLOAD_PARALLEL" -n2 sh -c \
          'mkdir -p raw/refseq/$0 && curl -L --fail --retry 5 --retry-delay 3 --retry-all-errors -s "$1" -o raw/refseq/$0/$(basename "$1") 2>/dev/null || true' \
        2>&1 | tee -a "$LOGDIR/downloads.log" >/dev/null

    # Subset manifest to this chunk
    head -1 work/manifests/all_genomes.tsv > "$LOGDIR/chunk_manifest.tsv"
    awk -F'\t' 'NR==FNR{ids[$1]=1; next} ids[$1]' "$chunk_url_file" "$LOGDIR/refseq_manifest_rows.tsv" >> "$LOGDIR/chunk_manifest.tsv"

    # Scan
    log "chunk $CHUNK_ID: scanning"
    "$PYTHON" scripts/scan_telotrons.py \
        --manifest "$LOGDIR/chunk_manifest.tsv" \
        --canonical-motifs work/manifests/canonical_motifs.tsv \
        --refseq-dir raw/refseq --tara-dir raw/tara \
        --motifs "$MOTIFS" \
        --threads "$SCAN_THREADS" \
        --loci "$LOCI_OUT" \
        --introns "$INTRONS_OUT" \
        --summary "$SUMMARY_OUT" \
        2>&1 | tee -a "$LOGDIR/scan.log" | tail -3

    # Cleanup
    log "chunk $CHUNK_ID: cleaning up raw/refseq"
    awk -F'\t' '{print $1}' "$chunk_url_file" | while read gid; do
        rm -rf "raw/refseq/$gid"
    done
    log "chunk $CHUNK_ID: done (free $(free_gb)G after cleanup)"
done

# ─────────────────────────────────────────────────────────────────────────
log "=== 4. Concatenate chunk outputs ==="
concat() {
    local out="$1"; shift
    local first="$1"; shift
    head -1 "$first" > "$out"
    for f in "$first" "$@"; do
        tail -n +2 "$f" >> "$out"
    done
    log "wrote $out ($(wc -l < "$out") lines)"
}

concat work/results/all_telotron_loci.tsv \
       work/results/_chunks/tara_loci.tsv work/results/_chunks/refseq_loci_*.tsv
concat work/results/all_introns_scanned.tsv \
       work/results/_chunks/tara_introns.tsv work/results/_chunks/refseq_introns_*.tsv
concat work/results/all_species_raw_summary.tsv \
       work/results/_chunks/tara_summary.tsv work/results/_chunks/refseq_summary_*.tsv

# ─────────────────────────────────────────────────────────────────────────
log "=== 5. Filter + analyze ==="
"$PYTHON" scripts/filter_final_set.py \
    --loci work/results/all_telotron_loci.tsv \
    --introns work/results/all_introns_scanned.tsv \
    --summary work/results/all_species_raw_summary.tsv \
    --min-repeat-frac 0.85 --bidir-min-repeat-frac 0.40 --bidir-min-hits 3 \
    --collapse-unique-loci --require-terminal-motif-match \
    --final work/results/final_telotron_set.tsv \
    --species work/results/final_species_summary.tsv \
    --negatives work/results/final_negative_controls.tsv \
    2>&1 | tee -a "$LOGDIR/filter.log"

"$PYTHON" scripts/dedup_telotrons.py \
    --final work/results/final_telotron_set.tsv \
    --refseq-dir raw/refseq --tara-dir raw/tara \
    --out-final work/results/final_telotron_set_dedup.tsv \
    --out-log work/results/dedup_log.tsv \
    --threads "$SCAN_THREADS" \
    2>&1 | tee -a "$LOGDIR/dedup.log"

# Replace working file with dedup set; preserve pre-dedup for provenance.
mv work/results/final_telotron_set.tsv work/results/final_telotron_set.prededup.tsv
cp work/results/final_telotron_set_dedup.tsv work/results/final_telotron_set.tsv

"$PYTHON" scripts/analyze_telotrons.py \
    --final work/results/final_telotron_set.tsv \
    --introns work/results/all_introns_scanned.tsv \
    --species work/results/final_species_summary.tsv \
    --boundary-kmers work/results/boundary_kmer_enrichment.tsv \
    --distance work/results/distance_to_end.tsv \
    --architecture work/results/architecture_summary.tsv \
    --threads "$SCAN_THREADS" \
    2>&1 | tee -a "$LOGDIR/analyze.log"

# ─────────────────────────────────────────────────────────────────────────
log "=== 6. Re-download positive species for FASTA extracts + per-species figures ==="

# Positive species = those with at least one telotron after filtering, OR being scanned
# for terminal-motif figures, OR being non-trivial. To keep disk bounded, take only
# genomes with >0 final telotrons.
awk -F'\t' 'NR>1 && $12+0 > 0 {print $1}' work/results/final_species_summary.tsv > "$LOGDIR/positive_ids.txt"
N_POS=$(wc -l < "$LOGDIR/positive_ids.txt")
log "$N_POS positive species — re-downloading for downstream extracts"

if [ "$N_POS" -gt 0 ]; then
    awk -F'\t' 'NR==FNR{ids[$1]=1; next} ids[$1]' "$LOGDIR/positive_ids.txt" "$LOGDIR/refseq_urls_sorted.tsv" \
        > "$LOGDIR/positive_urls.tsv"
    awk -F'\t' '{print $1"\t"$2"\n"$1"\t"$3}' "$LOGDIR/positive_urls.tsv" \
      | xargs -P "$DOWNLOAD_PARALLEL" -n2 sh -c \
          'mkdir -p raw/refseq/$0 && curl -L --fail --retry 5 --retry-delay 3 --retry-all-errors -s "$1" -o raw/refseq/$0/$(basename "$1") 2>/dev/null || true'
    log "positive RefSeq genomes re-downloaded ($(du -sh raw/refseq | awk '{print $1}'))"
fi

log "=== 7. Architecture classification, FASTA extracts, figures ==="
"$PYTHON" scripts/classify_telotron_architecture.py \
    --final work/results/final_telotron_set.tsv \
    --refseq-dir raw/refseq --tara-dir raw/tara \
    --out-loci work/results/final_telotron_set_architecture.tsv \
    --out-kmers work/results/boundary_kmers_by_architecture.tsv \
    2>&1 | tee -a "$LOGDIR/architecture.log"

"$PYTHON" scripts/extract_telotron_fasta.py \
    --final work/results/final_telotron_set_architecture.tsv \
    --refseq-dir raw/refseq --tara-dir raw/tara \
    --fasta-dir work/results/telotron_fasta \
    --flanked-dir work/results/telotron_flanked \
    2>&1 | tee -a "$LOGDIR/extract_fasta.log"

"$PYTHON" scripts/blast_linkers.py \
    --arch-tsv work/results/final_telotron_set_architecture.tsv \
    --refseq-dir raw/refseq --tara-dir raw/tara \
    --workdir work/results/blast_linkers \
    --out-own work/results/linker_blast_hits_own_genome.tsv \
    --out-all work/results/linker_blast_hits_all_genomes.tsv \
    --threads "$SCAN_THREADS" \
    2>&1 | tee -a "$LOGDIR/blast_linkers.log"

"$PYTHON" scripts/msa_telotron_regions.py \
    --arch-tsv work/results/final_telotron_set_architecture.tsv \
    --refseq-dir raw/refseq --tara-dir raw/tara \
    --outdir work/results/msa_regions \
    --threads 4 \
    2>&1 | tee -a "$LOGDIR/msa_regions.log"

"$PYTHON" scripts/find_interstitial_arrays.py \
    --manifest work/manifests/all_genomes.tsv \
    --refseq-dir raw/refseq --tara-dir raw/tara \
    --motifs "$MOTIFS" \
    --out work/results/interstitial_arrays.tsv \
    --threads "$SCAN_THREADS" \
    2>&1 | tee -a "$LOGDIR/interstitial.log"

"$PYTHON" scripts/plot_interstitial_boundary_kmers.py \
    --arrays work/results/interstitial_arrays.tsv \
    --outdir work/results/figures/interstitial_boundary_kmers \
    2>&1 | tee -a "$LOGDIR/plot_interstitial.log"

mkdir -p work/results/figures
"$PYTHON" scripts/plot_telotrons.py \
    --species work/results/final_species_summary.tsv \
    --final work/results/final_telotron_set.tsv \
    --kmers work/results/boundary_kmer_enrichment.tsv \
    --distance work/results/distance_to_end.tsv \
    --architecture work/results/architecture_summary.tsv \
    --outdir work/results/figures 2>&1 | tee -a "$LOGDIR/plot_telotrons.log"

"$PYTHON" scripts/plot_boundary_kmers_by_arch.py \
    --kmers work/results/boundary_kmers_by_architecture.tsv \
    --species work/results/final_species_summary.tsv \
    --outdir work/results/figures/boundary_kmers_by_arch \
    2>&1 | tee -a "$LOGDIR/plot_boundary.log"

"$PYTHON" scripts/plot_pipeline_stages.py \
    --manifests work/manifests/all_genomes.tsv \
    --raw-summary work/results/all_species_raw_summary.tsv \
    --species-final work/results/final_species_summary.tsv \
    --distance work/results/distance_to_end.tsv \
    --architecture-loci work/results/final_telotron_set_architecture.tsv \
    --outdir work/results/figures/pipeline_stages \
    2>&1 | tee -a "$LOGDIR/plot_stages.log"

# Terminal motif figures (requires raw FASTAs — already re-downloaded for positives)
"$PYTHON" scripts/plot_terminal_motif_density.py \
    --summary work/results/all_species_raw_summary.tsv \
    --loci work/results/all_telotron_loci.tsv \
    --refseq-dir raw/refseq --tara-dir raw/tara \
    --motifs "$MOTIFS" \
    --outdir work/results/figures/terminal_motif \
    --threads "$SCAN_THREADS" \
    2>&1 | tee -a "$LOGDIR/plot_terminal.log"

# Final cleanup: drop the re-downloaded positives too
log "=== 8. Final cleanup: removing positive-species raw FASTAs ==="
if [ "$N_POS" -gt 0 ]; then
    while read gid; do
        rm -rf "raw/refseq/$gid"
    done < "$LOGDIR/positive_ids.txt"
fi

log "=== ALL DONE ==="
log "work/results/ size: $(du -sh results | awk '{print $1}')"
log "raw/ size:     $(du -sh raw | awk '{print $1}')"
log "free disk:     $(free_gb)G"
