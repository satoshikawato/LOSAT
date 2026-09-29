#!/usr/bin/env python3
"""Compare LOSAT with NCBI in outfmt 6 on every outfmt 0 fixture (comparison only).

A fixture can only be byte-identical in outfmt 0 when the search finds the same HSPs,
so this runs each manifest search with -outfmt 6 in both programs (from LOSAT/) and
reports the fixtures whose tabular output differs or whose run fails.

Usage: precheck_hits.py --bin-dir DIR --losat LOSAT
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from run_oracle import ENGINE, read_manifest, search_argv


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--losat", type=Path, required=True)
    args = parser.parse_args()
    differing = []
    print("fixture_id\tncbi_exit\tncbi_lines\tlosat_exit\tlosat_lines\tresult\tlosat_stderr")
    for row in read_manifest()[2]:
        argv = [*search_argv(row)[:-1], "6"]
        ncbi = subprocess.run([str(args.bin_dir / row["program"]), *argv], cwd=ENGINE, capture_output=True)
        losat = subprocess.run([str(args.losat), row["program"], *argv], cwd=ENGINE, capture_output=True)
        same = ncbi.returncode == losat.returncode == 0 and ncbi.stdout == losat.stdout
        if not same:
            differing.append(row["fixture_id"])
        stderr = losat.stderr.decode(errors="replace").strip().splitlines()
        print("\t".join([row["fixture_id"], str(ncbi.returncode), str(len(ncbi.stdout.splitlines())),
                         str(losat.returncode), str(len(losat.stdout.splitlines())),
                         "same" if same else "DIFF", stderr[0] if stderr else ""]))
    print(f"# fixtures={len(read_manifest()[2])} differing={differing}")
    return 1 if differing else 0


if __name__ == "__main__":
    sys.exit(main())
