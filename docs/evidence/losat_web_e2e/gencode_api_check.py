#!/usr/bin/env python3
"""Check TBLASTN's non-default `-db_gencode` rows of the S08+ sweep with NCBI's C++ API (comparison only).

AGENTS.md's approved TBLASTN decision (`PD-TLOSAN-LOCAL-GENCODE-32`) makes a local
`-subject` search translate the subjects with the selected `-db_gencode`, which NCBI's
CLI does only for a BLAST database. A database search is not an exact oracle for it: NCBI
`-db` and `-subject` differ even at code 1 for some options (the S08+ sweep found 39
gencode rows where LOSAT differs from NCBI's `-db` search). TLOSAN Stage E/G checked the
decision with a comparison-only C++ API oracle that runs NCBI's local-subject search
(`CLocalDbAdapter(..., dbscan_mode=true)`, tblastn_app.cpp's path) with the selected code
supplied to the sequence source (`docs/evidence/tlosan_stage_e/tblastn_stage_e_local_oracle.cpp`,
built by `gates/build_api_oracle.sh`), at the CLI's query batches
(`docs/evidence/tlosan_stage_g/run_batch_api_oracle.py`) and with Stage E's outfmt 0
calibration (`calibrate_pairwise`).

The check uses the sweep's TBLASTN inputs and every `-db_gencode` spelling of the sweep
that LOSAT runs, in outfmt 0, 6 and 7:

1. calibration: the API at code 1 equals NCBI's CLI `-subject` search (stdout);
2. for each code: LOSAT's stdout equals the API's at that code, and LOSAT's stderr equals
   NCBI's CLI stderr for the same arguments (code 1's for ID 32, which the CLI rejects).

Prints one row per case ("same" or "DIFF ..."); exits 1 when any row is not "same".

Usage: gencode_api_check.py --api ORACLE --bin-dir DIR --losat LOSAT
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "docs/evidence/losat_web_e2e"))
sys.path.insert(0, str(REPO / "docs/evidence/tlosan_stage_e"))
from option_sweep import GENCODES, INPUTS  # noqa: E402
from run_oracle import ENGINE  # noqa: E402
from run_stage_e_codes import calibrate_pairwise  # noqa: E402

BATCH_WRAPPER = REPO / "docs/evidence/tlosan_stage_g/run_batch_api_oracle.py"


def ncbi_integer(spelling: str) -> int | None:
    """The value NCBI's integer argument reads (decimal, or hex with 0x), or None."""
    text = spelling.lstrip("+")
    try:
        return int(text, 16) if text.lower().startswith("0x") else int(text, 10)
    except ValueError:
        return None


def api_output(api: Path, argv: list[str], code: int, fmt: str) -> bytes:
    """NCBI's local-subject API oracle's stdout for the query and subject of `argv` at
    `code`, at the CLI's query batches (calibrated for outfmt 0)."""
    query, subject = argv[argv.index("-query") + 1], argv[argv.index("-subject") + 1]
    result = subprocess.run([sys.executable, str(BATCH_WRAPPER), str(api.resolve()), query, subject, str(code), fmt],
                            cwd=ENGINE, capture_output=True, stdin=subprocess.DEVNULL)
    if result.returncode:
        raise RuntimeError(f"API oracle failed at code {code} outfmt {fmt}: {result.stderr.decode(errors='replace')}")
    return calibrate_pairwise(result.stdout) if fmt == "0" else result.stdout


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--api", type=Path, required=True)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--losat", type=Path, required=True)
    args = parser.parse_args()
    inputs = INPUTS["tblastn"]
    query, subject = inputs[inputs.index("-query") + 1], inputs[inputs.index("-subject") + 1]
    failures = 0

    def run(argv: list[str]) -> subprocess.CompletedProcess:
        return subprocess.run(argv, cwd=ENGINE, capture_output=True, stdin=subprocess.DEVNULL)

    def api(code: int, fmt: str) -> bytes:
        return api_output(args.api, ["-query", query, "-subject", subject], code, fmt)

    for fmt in ("0", "6", "7"):
        cli = run([str(args.bin_dir / "tblastn"), *inputs, "-outfmt", fmt])
        control = api(1, fmt)
        status = "same" if cli.returncode == 0 and control == cli.stdout else "DIFF calibration"
        failures += status != "same"
        print(f"calibration\t1\t{fmt}\t{status}", flush=True)
        for spelling in GENCODES:
            code = ncbi_integer(spelling)
            ours = run([str(args.losat.resolve()), "tblastn", *inputs, "-db_gencode", spelling, "-outfmt", fmt])
            if ours.returncode or code in (None, 1):
                continue
            ncbi = run([str(args.bin_dir / "tblastn"), *inputs, "-db_gencode", "1" if code == 32 else spelling,
                        "-outfmt", fmt])
            expected = api(code, fmt)
            if ours.stdout != expected:
                status = "DIFF stdout"
            elif ncbi.returncode or ours.stderr != ncbi.stderr:
                status = f"DIFF stderr (NCBI CLI exit {ncbi.returncode})"
            else:
                status = "same"
            failures += status != "same"
            print(f"code\t{spelling}\t{fmt}\t{status}", flush=True)
    print(f"# failures {failures}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
