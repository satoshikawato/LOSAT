#!/usr/bin/env python3
"""Compare LOSAT with NCBI under CTOOLKIT_COMPATIBLE for every outfmt 0/7 fixture (comparison only).

Runs every row of LOSAT/tests/outfmt0_manifest.tsv (all programs; the rows of the approved
db_gencode deviation are left out, as NCBI's local -subject search does not apply their
subject code) with CTOOLKIT_COMPATIBLE unset, set to 1 and set to the empty string, with
NCBI BLAST+ and with LOSAT, and compares stdout, stderr and the exit status. NCBI's
showdefline.cpp writes "(bits)" in the description table header when the variable is set,
even to the empty string; every other variable that changes NCBI's report is unset.

Usage: ctoolkit_compare.py --bin-dir DIR --losat LOSAT [--jobs N]
"""
from __future__ import annotations

import argparse
import concurrent.futures
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "losat_web_e2a"))
import run_oracle  # noqa: E402

REPORT_ENV = ("BL2SEQ_LEGACY", "CTOOLKIT_COMPATIBLE", "OLD_FSC", "BATCH_SIZE", "CHUNK_SIZE",
              "OVERLAP_CHUNK_SIZE", "ADAPTIVE_CBS", "PRE_FETCH_SEQS_LIMIT")


def environment(value: str | None) -> dict[str, str]:
    env = {key: val for key, val in os.environ.items()
           if key not in REPORT_ENV and not key.startswith(("LOSAT_", "RAYON_"))}
    if value is not None:
        env["CTOOLKIT_COMPATIBLE"] = value
    return env


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--losat", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    args = parser.parse_args()
    if (Path.home() / ".ncbirc").exists():
        raise SystemExit("remove ~/.ncbirc before comparing")
    rows = run_oracle.read_manifest()[2]
    compared = [row for row in rows if row.get("contract") != run_oracle.DEVIATION]
    runs = [(row, value) for row in compared for value in (None, "1", "")]

    def one(run):
        row, value = run
        argv = run_oracle.search_argv(row)
        ncbi = subprocess.run([str(args.bin_dir / row["program"]), *argv], cwd=run_oracle.ENGINE,
                              capture_output=True, env=environment(value))
        losat = subprocess.run([str(args.losat.resolve()), row["program"], *argv], cwd=run_oracle.ENGINE,
                               capture_output=True, env=environment(value))
        same = (ncbi.stdout, ncbi.stderr, ncbi.returncode) == (losat.stdout, losat.stderr, losat.returncode)
        label = "unset" if value is None else repr(value)
        return f"{row['fixture_id']}\t{label}\t{'same' if same else 'DIFF'}\t(bits)={b'(bits)' in ncbi.stdout}"

    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        lines = list(pool.map(one, runs))
    print("fixture_id\tCTOOLKIT_COMPATIBLE\tresult\tncbi_header")
    print("\n".join(lines))
    differing = [line for line in lines if "\tDIFF\t" in line]
    print(f"# rows={len(compared)} runs={len(lines)} differing={len(differing)} "
          f"left_out_deviation_rows={len(rows) - len(compared)}")
    return 1 if differing else 0


if __name__ == "__main__":
    sys.exit(main())
