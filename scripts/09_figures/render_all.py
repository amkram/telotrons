"""Render every figure but S3 (stage 10), plus the captions. Run after prepare_architecture.py."""

from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
# each script runs in its own interpreter, and any failure stops the run
for name in [
    'trees.py',
    'architecture.py',
    'its.py',
    'supplements.py',
    'figure_s1_unclassified.py',
    'captions.py',
]:
    subprocess.run([sys.executable, str(HERE / name)], check=True)
