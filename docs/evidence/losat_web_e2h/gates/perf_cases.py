#!/usr/bin/env python3
"""The E2h performance cases: docs/evidence/losat_web_e2b/perf_cases.py (E2b -> E2c -> E1c -> measure_perf)
with the read-heavy cases of stage E2h, which the new FASTA reader (all four programs) must not slow down:

  blastn-q100k          100000 nucleotide queries of 300 nt against LC738884 (megablast, outfmt 6)
  blastn-genome-1line   AP027152 against a 50 Mb subject on ONE line
  blastn-genome-80col   the same subject in 80-column lines (same output bytes as the one-line case)
  blastp-many           2000 protein queries of 300 aa against PajaWSV.faa
  blastn-q100k-stdin-file, blastn-q100k-stdin-pipe
                        blastn-q100k with `-query -` (standard input redirected from the file, or a pipe): the
                        reader refills a stream without `in_avail` one byte at a time. Native only (the WASI
                        runners give the program no standard input): run them with PERF_MODES=native.

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
for name, file in (("blastn-genome-1line", "genome_1line.fna"), ("blastn-genome-80col", "genome_80col.fna")):
    measure_perf.FIXTURES[name] = [
        "blastn", "-task", "megablast", "-query", str(FASTA / "AP027152.fasta"),
        "-subject", str(INPUTS / file), "-outfmt", "6"]
measure_perf.FIXTURES["blastp-many"] = [
    "blastp", "-query", str(INPUTS / "protein_many.faa"), "-subject", str(FASTA / "PajaWSV.faa"), "-outfmt", "6"]

STDIN = {"blastn-q100k-stdin-file": 'f=$1; shift; exec "$@" < "$f"', "blastn-q100k-stdin-pipe": 'f=$1; shift; cat "$f" | "$@"'}
for name in STDIN:
    measure_perf.FIXTURES[name] = [
        "blastn", "-task", "megablast", "-query", "-", "-subject", str(FASTA / "LC738884.fasta"), "-outfmt", "6"]
_command = measure_perf.command


def command(mode, binary, argv, threads, out):
    """measure_perf.command, with the query of the standard-input cases fed through sh."""
    cmd = _command(mode, binary, argv, threads, out)
    for name, script in STDIN.items():
        if argv is measure_perf.FIXTURES[name]:
            if mode != "native":
                raise SystemExit(f"{name} runs in native mode only (PERF_MODES=native)")
            return ["sh", "-c", script, "sh", str(INPUTS / "q100k_300nt.fna")] + cmd
    return cmd


measure_perf.command = command

if os.environ.get("PERF_MODES"):
    wanted = set(os.environ["PERF_MODES"].split(","))
    measure_perf.MODES = tuple(m for m in measure_perf.MODES if m[0] in wanted)

if __name__ == "__main__":
    sys.exit(measure_perf.main())
