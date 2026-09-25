#!/usr/bin/env python3
"""Replay saved pinned-NCBI Stage D diagnostic streams in fresh temporary directories."""
from pathlib import Path
import hashlib
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
RUNNERS = (
    ("masking_options_20260925/run_ncbi_masking.py", "masking_options_20260925/run_20260925"),
    ("option_cross_20260925/run_ncbi_options.py", "option_cross_20260925/run_20260925"),
    ("alternate_matrix_20260925/run_ncbi_alternate.py", "alternate_matrix_20260925/run_20260925"),
    ("long_subject_20260925/run_ncbi_long.py", "long_subject_20260925/run_20260925"),
    ("remaining_local_20260925/run_ncbi_remaining.py", "remaining_local_20260925/run_20260925"),
    ("extended_chunks_20260925/run_ncbi_extended_chunks.py", "extended_chunks_20260925/run_20260925"),
)

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    lines = []
    with tempfile.TemporaryDirectory(prefix="tlosan-d-replay-") as temp:
        for index, (runner, saved) in enumerate(RUNNERS):
            fresh = Path(temp) / f"run{index}"
            subprocess.run([sys.executable, str(ROOT / runner), str(fresh)], check=True)
            expected = ROOT / saved
            names = sorted(path.name for path in expected.iterdir() if path.is_file())
            assert names == sorted(path.name for path in fresh.iterdir() if path.is_file()), saved
            for name in names:
                assert sha(expected / name) == sha(fresh / name), (saved, name)
            lines.append(f"{saved}: {len(names)} files byte-identical")
    print("\n".join(lines))

if __name__ == "__main__": main()
