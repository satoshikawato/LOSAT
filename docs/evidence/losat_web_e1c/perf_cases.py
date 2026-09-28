#!/usr/bin/env python3
"""The S04 performance cases: docs/evidence/losat_web_e1a/measure_perf.py with one more
BLASTN fixture. The `blastn` case of measure_perf.py (AP027152 x LC738884 megablast)
finds no hit and takes about 0.05 s natively, so it measures start-up only. The added
cases are the Gate A EDL933 x Sakai megablast search (LOSAT/tests/blastn_parity_manifest.tsv,
about 1 s natively, 5.6 Mb query), in outfmt 6 and outfmt 7.

Usage: the same as measure_perf.py (run / check); the extra cases are
blastn-large and blastn-large-fmt7.
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("measure_perf", HERE.parent / "losat_web_e1a" / "measure_perf.py")
measure_perf = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(measure_perf)

LARGE = ["blastn", "-task", "megablast", "-query", str(measure_perf.FASTA / "EDL933.fna"),
         "-subject", str(measure_perf.FASTA / "Sakai.fna")]
measure_perf.FIXTURES["blastn-large"] = LARGE + ["-outfmt", "6"]
measure_perf.FIXTURES["blastn-large-fmt7"] = LARGE + ["-outfmt", "7"]

if __name__ == "__main__":
    sys.exit(measure_perf.main())
