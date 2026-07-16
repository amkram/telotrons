#!/bin/bash
# Submit the telotron pipeline (or a specific target) to slurm.
#
# Every rule runs as its own sbatch job; per-rule cpus/mem/time come from
# profiles/slurm/config.yaml. The driver snakemake stays alive on the login
# node (or a small persistent job) to coordinate. See slurm/persistent.sh
# for running the driver itself as a slurm job.
#
# Usage:
#   ./slurm/submit.sh                           # default target: `all`
#   ./slurm/submit.sh scan_all                  # single target
#   ./slurm/submit.sh scan_all -- --dry-run     # extra snakemake args after --
#   ./slurm/submit.sh -n                        # dry-run driver locally
#
# Env overrides (see slurm/site.sh for defaults):
#   SLURM_PARTITION=cpu SLURM_ACCOUNT=proj ./slurm/submit.sh

set -euo pipefail
cd "$(dirname "$0")/.."           # repo root

source slurm/site.sh              # partition, account, module loads

TARGET="all"
EXTRA=()
if [[ $# -gt 0 ]] && [[ "$1" != -* ]]; then
  TARGET="$1"; shift
fi
if [[ "${1-}" == "--" ]]; then shift; EXTRA=("$@"); fi

mkdir -p work/logs/slurm

# Rewrite $VAR placeholders in the profile from live env (partition, account).
# This lets a user swap clusters without editing the profile yaml.
tmp_profile=$(mktemp -d)
trap 'rm -rf "$tmp_profile"' EXIT
cp -r profiles/slurm/* "$tmp_profile/"
sed -i \
  -e "s|\\\$SLURM_PARTITION|$SLURM_PARTITION|g" \
  -e "s|\\\$SLURM_ACCOUNT|${SLURM_ACCOUNT:-DUMMY}|g" \
  "$tmp_profile/config.yaml"

# If no account set, strip the --account flag entirely (avoids "sbatch: error").
if [[ -z "$SLURM_ACCOUNT" ]]; then
  sed -i "/--account=/d" "$tmp_profile/config.yaml"
fi

echo "=== telotron slurm submit ==="
echo "  target:    $TARGET"
echo "  partition: $SLURM_PARTITION"
echo "  account:   ${SLURM_ACCOUNT:-(none)}"
echo "  extra:     ${EXTRA[*]:-(none)}"
echo "  logs:      work/logs/slurm/{rule}.{jobid}.{out,err}"
echo "============================="

exec snakemake \
  --profile "$tmp_profile" \
  "${EXTRA[@]}" \
  "$TARGET"
