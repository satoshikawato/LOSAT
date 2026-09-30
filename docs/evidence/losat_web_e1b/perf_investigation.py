#!/usr/bin/env python3
"""The S03 performance investigation: docs/evidence/losat_web_e1a/measure_perf.py with
one more pair of cases, the full AvCLPV proteins (120 queries) against the PsCLPV
genome (the Stage G full real code-1 cross, docs/evidence/tlosan_stage_g/), in
outfmt 6 and outfmt 0. The Stage G benchmark fixture of measure_perf.py takes about
0.05 s natively, so process start-up dominates it; this fixture takes about 2 s.

Usage: the same as measure_perf.py (run / check); the extra cases are
tblastn-full and tblastn-full-fmt0.
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("measure_perf", HERE.parent / "losat_web_e1a" / "measure_perf.py")
measure_perf = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(measure_perf)

FULL = ["tblastn", "-task", "tblastn", "-query", str(measure_perf.FASTA / "AvCLPV.faa"),
        "-subject", str(measure_perf.FASTA / "PsCLPV.fasta")]
measure_perf.FIXTURES["tblastn-full"] = FULL + ["-outfmt", "6"]
measure_perf.FIXTURES["tblastn-full-fmt0"] = FULL + ["-outfmt", "0"]

if __name__ == "__main__":
    sys.exit(measure_perf.main())
