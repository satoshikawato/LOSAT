#!/usr/bin/env python3
"""The S08 performance cases: docs/evidence/losat_web_e2c/perf_cases.py with two more
TBLASTX cases. S08 searches TBLASTX queries in NCBI's 10002-nt batches, orders the hit
list as Blast_HitListUpdate and resolves ambiguous subject bases for the preliminary
search, so the added cases are a search of 4 queries in 3 batches (`tblastx-multi`, the
multi_query records of the outfmt 0 fixtures against LC741431) and a search of 260
subjects (`tblastx-many`, the many records of the outfmt 0 fixtures). Both give the same
outfmt 6 bytes before and after S08.

Usage: the same as measure_perf.py (run / check).
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("perf_cases", HERE.parent / "losat_web_e2c" / "perf_cases.py")
perf_cases = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(perf_cases)
measure_perf = perf_cases.measure_perf
OUTFMT0 = measure_perf.FASTA / "outfmt0"

measure_perf.FIXTURES["tblastx-multi"] = [
    "tblastx", "-query", str(OUTFMT0 / "tblastx_multi_query.fasta"),
    "-subject", str(measure_perf.FASTA / "LC741431.fasta"), "-outfmt", "6"]
measure_perf.FIXTURES["tblastx-many"] = [
    "tblastx", "-query", str(OUTFMT0 / "tblastx_many_query.fasta"),
    "-subject", str(OUTFMT0 / "tblastx_many_subject.fasta"), "-outfmt", "6"]

if __name__ == "__main__":
    sys.exit(measure_perf.main())
