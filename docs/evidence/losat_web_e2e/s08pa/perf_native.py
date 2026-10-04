#!/usr/bin/env python3
"""Native-only timing of the S08 performance cases (docs/evidence/losat_web_e2b/perf_cases.py)
for the S08+a changes, before the formal V-PERF of S08+b.

Same protocol as measure_perf.py (one untimed warmup per side, the two builds alternating on
every repetition, median and full range, output SHA-256 of both sides), but only the native
mode, so a single binary per side is enough.

Usage:
  perf_native.py run --before LOSAT --after LOSAT --out FILE.json [--cases blastp,...] [--repeat N]
  perf_native.py check FILE.json
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "perf_cases_e2b", HERE.parents[1] / "losat_web_e2b" / "perf_cases.py")
perf_cases = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(perf_cases)
measure_perf = perf_cases.measure_perf
measure_perf.MODES = (("native", 0, 1),)
# A longer TBLASTN search (300 subjects, about 4 s) than the 0.1-s `tblastn` case, so the
# gapped alignment shared by BLASTP, TBLASTN and BLASTX (TN-4) is a measurable share.
OUTFMT0 = measure_perf.FASTA / "outfmt0"
measure_perf.FIXTURES["tblastn-many-e2e"] = [
    "tblastn", "-query", str(OUTFMT0 / "e2e_protein_query.faa"),
    "-subject", str(OUTFMT0 / "e2e_many_subject.fna"), "-outfmt", "6"]

if __name__ == "__main__":
    argv = sys.argv[1:]
    for flag in ("--before", "--after"):
        if flag in argv:
            index = argv.index(flag) + 1
            argv[index] = ",".join([argv[index]] * 3)
    sys.argv = [sys.argv[0], *argv]
    sys.exit(measure_perf.main())
