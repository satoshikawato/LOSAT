#!/usr/bin/env python3
"""Compare BLASTN input handling with NCBI (Session S07+++b; comparison only).

This is docs/evidence/losat_web_e2f/check_inputs.py with the expectations of E2g. LOSAT
now reads BATCH_SIZE as NCBI (E2g T7), and reports a Karlin-Altschul table error of a later
query batch after the reports of the batches before (T8), so these cases that S07+ and
S07++ rejected are compared with NCBI ("same" or "same-error"). A few cases more check T7,
T8 and T10.

Usage: check_inputs.py --bin-dir DIR --losat LOSAT --work DIR
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("e2f_check_inputs", HERE.parent / "losat_web_e2f" / "check_inputs.py")
e2f = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(e2f)
e2c = e2f.e2c

E2G_EXPECT = {
    "audit2.batch_size_env": "same",
    "scoring_error.first_batch_invalid": "same-error",
}
e2f_cases = e2f.cases
F = "tests/fasta/outfmt0"
COMPACT = "tests/fasta/blastn_parity_compact.fasta"


def cases(work: Path) -> list[tuple]:
    rows = []
    for name, argv, expect, *extra in e2f_cases(work):
        rows.append((name, argv, E2G_EXPECT.get(name, expect), *extra))
    batch_invalid = ["-query", f"{F}/edge_batch_allN.fasta", "-subject", COMPACT, "-reward", "1", "-penalty", "-6"]
    rows += [
        ("e2g.t8.later_batch.fmt0", batch_invalid, "same-error"),
        ("e2g.t8.later_batch.fmt7", [*batch_invalid, "-outfmt", "7"], "same-error"),
        ("e2g.t7.batch_size_100", ["-query", f"{F}/multi_query.fasta", "-subject", f"{F}/multi_subject.fasta",
                                   "-outfmt", "6"], "same", {"BATCH_SIZE": "100"}),
        ("e2g.t7.chunk_size_500", ["-query", f"{F}/multi_query.fasta", "-subject", f"{F}/multi_subject.fasta",
                                   "-task", "blastn", "-outfmt", "6"], "same", {"CHUNK_SIZE": "500"}),
        ("e2g.t7.batch_size_text", ["-query", f"{F}/multi_query.fasta", "-subject", f"{F}/multi_subject.fasta",
                                    "-outfmt", "6"], "losat-rejects", {"BATCH_SIZE": "abc"}),
        ("e2g.t10.no_subject", ["-query", f"{F}/multi_query.fasta", "-outfmt", "6"], "same-error"),
    ]
    return rows


e2c.cases = cases

if __name__ == "__main__":
    sys.exit(e2c.main())
