#!/bin/bash
# Quick status view for a running (or completed) telotron slurm run.
#
# Usage:
#   ./slurm/monitor.sh              # squeue + DAG counts + last log tail
#   ./slurm/monitor.sh -f           # tail -f the newest driver log
#   ./slurm/monitor.sh <rule>       # show recent sacct rows for a rule name

set -u
cd "$(dirname "$0")/.."

USER_ID=${USER:-$(id -un)}
LOG_DIR="work/logs/slurm"
mkdir -p "$LOG_DIR"

if [[ "${1-}" == "-f" ]]; then
    newest=$(ls -t "$LOG_DIR"/*.out 2>/dev/null | head -1)
    if [[ -z "$newest" ]]; then
        echo "no logs under $LOG_DIR yet"
        exit 1
    fi
    echo "tailing $newest (Ctrl-C to stop)"
    exec tail -F "$newest"
fi

if [[ $# -gt 0 && "$1" != -* ]]; then
    rule="$1"
    echo "recent sacct rows matching name=snakemake-$rule*"
    sacct -u "$USER_ID" --format=JobID,JobName%40,State,Elapsed,MaxRSS,NodeList%20,End \
          --starttime="$(date -d '2 days ago' +%Y-%m-%d)" \
        | awk -v r="snakemake-$rule" 'NR<=2 || $2 ~ r'
    exit 0
fi

echo "=== squeue ($USER_ID) ==="
if command -v squeue >/dev/null; then
    squeue -u "$USER_ID" -o "%.10i %.30j %.10T %.10M %.6D %R" || true
    n_pend=$(squeue -u "$USER_ID" -h -t PENDING 2>/dev/null | wc -l)
    n_run=$( squeue -u "$USER_ID" -h -t RUNNING 2>/dev/null | wc -l)
    echo "  pending=$n_pend  running=$n_run"
else
    echo "  (squeue not available on this host)"
fi

echo
echo "=== per-rule log tallies (last 24h) ==="
find "$LOG_DIR" -maxdepth 1 -name '*.out' -mmin -1440 2>/dev/null \
    | sed -E 's|.*/([^.]+)\..*|\1|' \
    | sort | uniq -c | sort -rn | head -20 || echo "  (no logs in last 24h)"

echo
echo "=== DAG progress (snakemake --summary tail) ==="
if command -v snakemake >/dev/null; then
    snakemake --quiet -n --summary 2>/dev/null \
        | awk 'NR==1 || $4 !~ /^ok$/' \
        | head -30 \
        || echo "  (no snakemake state — has the driver started?)"
fi

echo
echo "=== newest driver log tail ==="
newest=$(ls -t "$LOG_DIR"/*.out 2>/dev/null | head -1)
if [[ -n "$newest" ]]; then
    echo "  $newest"
    tail -n 20 "$newest"
else
    echo "  (no logs under $LOG_DIR yet)"
fi
