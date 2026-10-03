#!/usr/bin/env python3
"""Compare BLASTN input handling and query batches for three tasks with NCBI (Session E2i; comparison only).

This runs the cases of docs/evidence/losat_web_e2g/check_inputs.py and, for each case that
gives no `-task`, two more: the same arguments with `-task dc-megablast` (case name suffixed
`.dc`) and with `-task blastn-short` (`.short`). The expectation of a case stays the one of
E2g (a case that expects `losat-rejects` still expects an explicit rejection by LOSAT), except
those that the tasks change (E2I_EXPECT, from the first run against 90c5f0181). The
output is that of E2c and E2g. The `.short` variants of the audit15.small_evalue cases are left out:
NCBI's blastn-short (word size 7, e-value 1000 by default) on the two full genomes of those cases
needs more than 14 GB of memory and exceeds the 60 s limit of the case runner.

Usage: check_inputs.py --bin-dir DIR --losat LOSAT --work DIR
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("e2g_check_inputs", HERE.parent / "losat_web_e2g" / "check_inputs.py")
e2g = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(e2g)
e2c = e2g.e2c
e2g_cases = e2g.cases
TASKS = (("dc", "dc-megablast"), ("short", "blastn-short"))
SKIP_SHORT = "audit15.small_evalue."
# The expectations that the tasks change (all are NCBI's own results, which LOSAT gives):
# - NCBI's errors where the base case searches: the scores of the case have no gap costs 5/2
#   in NCBI's tables (the DP tasks' defaults), dc-megablast's word sizes other than 11 and 12,
#   zero gap costs without greedy extension;
# - LOSAT's limit on greedy gap costs (`-task megablast` only) does not apply (the cases of
#   E2c that give the tasks themselves now expect `same` there);
# - with blastn-short's e-value 1000 the punctuation title has hits, where NCBI crashes
#   (approved exception 2 of PD-LOSAT-NCBI-DEFECTS).
E2I_EXPECT = {
    **{f"{case}.{suffix}": "same-error" for case in (
        "audit.reward_65538", "audit.large_divisible", "audit5.bare_0x.gaps", "audit7.bit_score_99.fmt0",
        "audit7.bit_score_99.fmt6", "audit7.bit_score_99.fmt7") for suffix in ("dc", "short")},
    "audit4.hex.word_size_16.dc": "same-error",
    "audit14.iupac_seed.word_size_24.dc": "same-error",
    **{f"{case}.{suffix}": "same" for case in ("audit.megablast_gap_max", "audit2.greedy_gap_limit")
       for suffix in ("dc", "short")},
    "audit12.crash_title_no_hit.fmt0.short": "exception-2",
}


def cases(work: Path) -> list[tuple]:
    rows = []
    for name, argv, expect, *extra in e2g_cases(work):
        rows.append((name, argv, E2I_EXPECT.get(name, expect), *extra))
        if "-task" not in argv:
            for suffix, task in TASKS:
                if suffix == "short" and name.startswith(SKIP_SHORT):
                    continue
                variant = f"{name}.{suffix}"
                rows.append((variant, [*argv, "-task", task], E2I_EXPECT.get(variant, expect), *extra))
    return rows


e2c.cases = cases

if __name__ == "__main__":
    sys.exit(e2c.main())
