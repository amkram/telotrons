# Site-specific slurm settings. Edit these for your cluster.
#
# Every telotron slurm invocation sources this file; nothing else in the repo
# hardcodes partition/account names.

# Slurm partition that per-rule jobs will be submitted to.
export SLURM_PARTITION="${SLURM_PARTITION:-cpu}"

# Slurm account / allocation to charge. Leave empty on clusters that don't
# require one; snakemake's cluster template drops the --account flag if unset.
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"

# Optional: notification email + events (BEGIN,END,FAIL,ALL,NONE).
# Set to empty to disable.
export SLURM_MAIL_USER="${SLURM_MAIL_USER:-}"
export SLURM_MAIL_TYPE="${SLURM_MAIL_TYPE:-END,FAIL}"

# Environment module or conda activation. Uncomment / edit for your site:
#   module load snakemake gt seqkit bedtools samtools hmmer miniprot blast+
# or:
#   conda activate telotrons
