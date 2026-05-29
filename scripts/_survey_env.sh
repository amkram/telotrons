# Shared environment for the hand-run survey drivers (run_full_survey.sh,
# run_test_survey.sh). Source this AFTER `cd`-ing to the repo root:
#     source scripts/_survey_env.sh
# Sets PATH/PYTHON to the telotrons conda env, the telomere MOTIFS string,
# default thread counts, and a log() helper. Set LOGDIR before the first log().

# Locate the telotrons conda env by the tool it ships (no hardcoded env hash,
# which rotates on every env rebuild). Falls back to ambient PATH.
for _d in ~/.snakemake-envs/telotrons/*/bin; do
    if [ -x "$_d/seqkit" ]; then
        export PATH="$_d:$PATH"
        PYTHON="$_d/python"
        break
    fi
done
: "${PYTHON:=python3}"

# Curated telomere motifs (G-rich strand); rotations + reverse complements are
# expanded downstream. Mirrors config.yaml: telomere_motifs.
MOTIFS=TTAGG,TCAGG,TTGGG,TTAGGG,TTAGGC,TTGGGG,TTTAGGG,TTTTAGGG,TTTTGGGG,TATTAGGG,TTATTGGG,TTATTAGGG,TTACTTGGG,TTATTGGGG,TTTTTTAGGG,TTTATTAGGG,TTAGGTTGGGG,CTCGGTTATGGG

SCAN_THREADS=${SCAN_THREADS:-32}
DOWNLOAD_PARALLEL=${DOWNLOAD_PARALLEL:-8}

log() { echo "[$(date '+%F %T')] $*" | tee -a "${LOGDIR:-logs}/run.log"; }
