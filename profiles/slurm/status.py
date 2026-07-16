#!/usr/bin/env python3
"""Cluster-status hook for `snakemake --profile profiles/slurm`.

Snakemake calls this once per running job with the parsable jobid it got back
from `sbatch --parsable`. The script must print exactly one of:
    running | success | failed
so snakemake can decide whether to move on, mark done, or retry.
"""
import shutil
import subprocess
import sys

RUN_STATES = {
    "PENDING", "RUNNING", "REQUEUED", "SUSPENDED",
    "CONFIGURING", "COMPLETING", "RESIZING",
}
OK_STATES = {"COMPLETED"}


def status(jobid: str) -> str:
    # First try `sacct` — authoritative even after the job leaves squeue.
    if shutil.which("sacct"):
        r = subprocess.run(
            ["sacct", "-j", jobid, "--format=State", "--noheader", "--parsable2"],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
        )
        line = (r.stdout or "").strip().splitlines()
        # The primary step is the first non-empty line; ignore .batch/.extern.
        for L in line:
            s = L.split("+", 1)[0].split(" ", 1)[0].strip()
            if not s:
                continue
            if s in RUN_STATES:
                return "running"
            if s in OK_STATES:
                return "success"
            return "failed"
    # Fallback: squeue (only sees live jobs).
    if shutil.which("squeue"):
        r = subprocess.run(
            ["squeue", "-j", jobid, "-h", "-o", "%T"],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
        )
        s = (r.stdout or "").strip()
        if s in RUN_STATES:
            return "running"
        if not s:
            # Not in squeue and sacct absent — assume completed.
            return "success"
    return "failed"


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: status.py <jobid>")
    print(status(sys.argv[1]))
