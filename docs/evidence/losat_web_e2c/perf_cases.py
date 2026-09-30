#!/usr/bin/env python3
"""The S07+ performance cases: docs/evidence/losat_web_e1c/perf_cases.py (which adds the
Gate A EDL933 x Sakai megablast search to docs/evidence/losat_web_e1a/measure_perf.py)
with two more BLASTN cases. S07+ reads the FASTA bytes once more for the input checks and
computes a Karlin block per query context, so the added cases are the same large search
in outfmt 0 (`blastn-large-fmt0`) and a search of 260 queries (`blastn-many`, the
many_subject records of the outfmt 0 fixtures against LC738884 with -task blastn).

Usage: the same as measure_perf.py (run / check).
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("perf_cases", HERE.parent / "losat_web_e1c" / "perf_cases.py")
perf_cases = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(perf_cases)
measure_perf = perf_cases.measure_perf

measure_perf.FIXTURES["blastn-large-fmt0"] = perf_cases.LARGE + ["-outfmt", "0"]
measure_perf.FIXTURES["blastn-many"] = [
    "blastn", "-task", "blastn", "-query", str(measure_perf.FASTA / "outfmt0" / "many_subject.fasta"),
    "-subject", str(measure_perf.FASTA / "LC738884.fasta"), "-outfmt", "6"]

if __name__ == "__main__":
    sys.exit(measure_perf.main())
