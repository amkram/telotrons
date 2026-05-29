#!/usr/bin/env bash
# Run blastn for each genome × architecture category, then filter self-matches.
# Results: $OUTDIR/{genome_dir}/{arch}.raw.tsv  and  {arch}.filtered.tsv
set -euo pipefail

FASTA_DIR="${1:-results/telotron_fasta}"
DB_DIR="${2:-results/blast_linkers/dbs}"
OUTDIR="${3:-results/blast_telotrons_by_arch}"
THREADS="${4:-16}"
RUNNER="scripts/blast_with_alignment.py"

mkdir -p "$OUTDIR"

for genome_dir in "$FASTA_DIR"/*/; do
    dirname=$(basename "$genome_dir")

    # Find db by looking for any db slug that is a prefix of the dir name.
    db=""
    for candidate in "$DB_DIR"/*.ndb "$DB_DIR"/*.nin; do
        [[ -f "$candidate" ]] || continue
        slug=$(basename "$candidate")
        slug="${slug%%.*}"   # strip extension(s)
        if [[ "$dirname" == "${slug}"* ]]; then
            db="$DB_DIR/$slug"
            break
        fi
    done

    if [[ -z "$db" ]]; then
        echo "skip $dirname: no matching db in $DB_DIR"
        continue
    fi

    mkdir -p "$OUTDIR/$dirname"

    for fa in "$genome_dir"*.fa; do
        [[ -f "$fa" ]] || continue
        arch=$(basename "$fa" .fa)
        n=$(grep -c "^>" "$fa" 2>/dev/null || echo 0)
        [[ "$n" -eq 0 ]] && continue

        out_aln="$OUTDIR/$dirname/${arch}.aligned.tsv"
        out_log="$OUTDIR/$dirname/${arch}.log"

        echo "  $dirname / $arch  ($n queries)"
        python "$RUNNER" \
            -db "$db" \
            -query "$fa" \
            -num_threads "$THREADS" \
            -word_size 7 \
            -evalue 1000 \
            -qcov_hsp_perc 95 \
            -perc_identity 98 \
            -out "$out_aln" \
            2>"$out_log"
        echo "    $(cat "$out_log")"
    done
done

echo "done — results in $OUTDIR"