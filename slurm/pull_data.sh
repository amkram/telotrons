#!/bin/bash
# Pull the already-downloaded genome data from a source host instead of
# re-fetching ~1.6 TB from NCBI (which took ~4 h of rate-limited curls).
#
# RUN THIS ON THE DESTINATION (the new/remote machine), from the repo root.
#
#   ./slurm/pull_data.sh                          # pull everything
#   ./slurm/pull_data.sh --dry-run                # show what would transfer
#   TELO_SRC_HOST=other.host ./slurm/pull_data.sh # different source
#
# Safe to re-run: rsync skips files already transferred, and --partial-dir
# resumes interrupted large files rather than restarting them. Run it under
# tmux/screen — a 1.6 TB pull will outlive your SSH session.
#
# WHY mtimes MATTER (do not drop -a): snakemake decides staleness by mtime.
# On the source, work/manifests/*.tsv are OLDER than data/raw/*/.done, which
# is what tells snakemake the downloads are current. rsync -a replays those
# timestamps exactly, so download_assemblies stays satisfied. Copy with
# anything that stamps fresh mtimes (scp, cp, tar without -p) and snakemake
# will re-download all 1.6 TB.
#
# WHY NO -z: 1,425 of the 1,600 GB are already .gz. Compressing compressed
# bytes burns CPU on both ends and transfers no faster.

set -euo pipefail
cd "$(dirname "$0")/.."           # repo root

SRC_HOST="${TELO_SRC_HOST:-silverbullet.ucsc.edu}"
SRC_USER="${TELO_SRC_USER:-alex}"
SRC_ROOT="${TELO_SRC_ROOT:-/scratch1/alex/telotrons}"
SRC="${SRC_USER}@${SRC_HOST}:${SRC_ROOT}"

DRY=()
if [[ "${1-}" == "--dry-run" || "${1-}" == "-n" ]]; then
  DRY=(--dry-run)
  echo "=== DRY RUN — nothing will be written ==="
fi

# -a  archive: recurse + preserve mtimes/perms  (mtimes are load-bearing, see above)
# -h  human-readable sizes
# --partial-dir  resume interrupted files instead of restarting them
# --info=progress2  one aggregate progress line rather than per-file spam
RSYNC_OPTS=(
  -ah
  --partial-dir=.rsync-partial
  --info=progress2
  "${DRY[@]}"
)
# Optional throttle on a shared uplink, e.g. TELO_BWLIMIT=50M
if [[ -n "${TELO_BWLIMIT:-}" ]]; then
  RSYNC_OPTS+=("--bwlimit=${TELO_BWLIMIT}")
fi

echo "=== telotron data pull ==="
echo "  source: ${SRC}"
echo "  dest:   $(pwd)"
echo "  bwlimit: ${TELO_BWLIMIT:-(none)}"
echo "=========================="

# 1. Manifests first (~3 MB). These are download_assemblies' INPUTS; without
#    them the manifests rule regenerates them with fresh mtimes, which makes
#    the .done sentinels look stale and triggers a full re-download.
echo
echo "--- [1/2] work/manifests (~3 MB) ---"
mkdir -p work/manifests
rsync "${RSYNC_OPTS[@]}" "${SRC}/work/manifests/" work/manifests/

# 2. The genome data (~1.6 TB: refseq 692G + genbank 736G + tara 30G).
#    The trailing slash on the source copies its CONTENTS into dest.
#    Hidden .done / .fna.done / .gff.done sentinels are included by -a and
#    are what keep snakemake from re-running the download rules.
echo
echo "--- [2/2] data/raw (~1.6 TB — expect hours; resumable) ---"
mkdir -p data/raw
rsync "${RSYNC_OPTS[@]}" "${SRC}/data/raw/" data/raw/

if [[ ${#DRY[@]} -gt 0 ]]; then
  echo
  echo "=== dry run complete ==="
  exit 0
fi

echo
echo "=== verifying ==="
for f in data/raw/refseq/.done data/raw/genbank/.done \
         data/raw/tara/.fna.done data/raw/tara/.gff.done \
         work/manifests/all_genomes.tsv; do
  if [[ -e "$f" ]]; then
    echo "  OK      $f"
  else
    echo "  MISSING $f   <-- snakemake will re-run its rule"
  fi
done
echo "  refseq dirs:  $(ls data/raw/refseq 2>/dev/null | wc -l)   (source had 2412)"
echo "  genbank dirs: $(ls data/raw/genbank 2>/dev/null | wc -l)   (source had 8108)"

cat <<'EOF'

=== next step ===
Confirm snakemake agrees the downloads are done — this must NOT list
download_assemblies / tara_archives / manifests:

    snakemake -n --quiet

If it wants to re-run them, mtimes were not preserved; re-pull with -a.
EOF
