#!/usr/bin/env python3
"""The E2h performance cases: docs/evidence/losat_web_e2b/perf_cases.py (E2b -> E2c -> E1c -> measure_perf)
with the read-heavy cases of stage E2h, which the new FASTA reader (all four programs) must not slow down:

  blastn-q100k          100000 nucleotide queries of 300 nt against LC738884 (megablast, outfmt 6)
  blastn-genome-1line   AP027152 against a 5 Mb subject on ONE line
  blastn-genome-80col   the same subject in 80-column lines (same output bytes as the one-line case)
  blastp-many           20000 protein queries of 300 aa against PajaWSV.faa

Inputs come from gen_perf_inputs.py ($BUILD_ROOT/sfb-e2h/perf-inputs/ or $PERF_INPUTS). PERF_MODES
(comma list of native, serial-wasi, threaded-wasi) restricts the modes of a run.

Usage: the same as measure_perf.py (run / check).
"""
import importlib.util
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("perf_cases_e2b", HERE.parent.parent / "losat_web_e2b" / "perf_cases.py")
perf_cases = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(perf_cases)
measure_perf = perf_cases.measure_perf
FASTA = measure_perf.FASTA
INPUTS = Path(os.environ.get("PERF_INPUTS") or Path(os.environ.get("BUILD_ROOT", "/home/kawato/.cache/losat-work")) / "sfb-e2h" / "perf-inputs")

measure_perf.FIXTURES["blastn-q100k"] = [
    "blastn", "-task", "megablast", "-query", str(INPUTS / "q100k_300nt.fna"),
    "-subject", str(FASTA / "LC738884.fasta"), "-outfmt", "6"]
for name, file in (("blastn-genome-1line", "genome5mb_1line.fna"), ("blastn-genome-80col", "genome5mb_80col.fna")):
    measure_perf.FIXTURES[name] = [
        "blastn", "-task", "megablast", "-query", str(FASTA / "AP027152.fasta"),
        "-subject", str(INPUTS / file), "-outfmt", "6"]
measure_perf.FIXTURES["blastp-many"] = [
    "blastp", "-query", str(INPUTS / "protein_many.faa"), "-subject", str(FASTA / "PajaWSV.faa"), "-outfmt", "6"]

if os.environ.get("PERF_MODES"):
    wanted = set(os.environ["PERF_MODES"].split(","))
    measure_perf.MODES = tuple(m for m in measure_perf.MODES if m[0] in wanted)

if __name__ == "__main__":
    sys.exit(measure_perf.main())
