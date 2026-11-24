#!/usr/bin/env python3
"""Compare BLASTN input handling and query batches for three tasks with NCBI (Session E2i; comparison only).

This runs the cases of docs/evidence/losat_web_e2g/check_inputs.py and, for each case that
gives no `-task`, two more: the same arguments with `-task dc-megablast` (case name suffixed
`.dc`) and with `-task blastn-short` (`.short`). The expectation of a case stays the one of
E2g (a case that expects `losat-rejects` still expects an explicit rejection by LOSAT), except
those that the tasks change (E2I_EXPECT, from the first run against 90c5f0181) and those that
the E2c runner's reading of NCBI in LOSAT's words and the ports since E2g changed (E2I_NOW, from
the run against 2cd5a1906 on 2026-11-24). For the hit list sizes of approved exception 19, NCBI
runs at the size whose output LOSAT gives (E2I_NEAR). The output is that of E2c and E2g. The
`.short` variants of the audit15.small_evalue cases are left out: NCBI's blastn-short (word size 7, e-value 1000 by default) on the two full genomes of those cases
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

# The expectations that the E2c runner's reading of NCBI in LOSAT's words and the later ports changed
# (NCBI's own results, measured against dev 2cd5a1906 plus that runner on 2026-11-24). A case of E2g
# that expected `losat-rejects` or `arg-error` now has LOSAT give NCBI's result or NCBI's error (the
# FASTA reader and the options that E2h's FIXED list and the later lines ported); the names are those
# of the E2g case, and each also stands for its `.dc` and `.short` variants unless listed in the
# second tuple of its class (only the variants that changed).
NOW_SAME = (
    "audit.blank_defline", "audit.empty_defline", "audit.empty_query.subject_x", "audit.empty_record",
    "audit.empty_subject_record", "audit.leading_blank_line", "audit.piped_empty_query",
    "audit10.ideographic_space.subject", "audit10.nbsp.query", "audit10.title_9.fmt0", "audit3.evalue_hex",
    "audit4.fifo.empty_query", "audit4.outfmt.delimiter", "audit4.stdin_both.file", "audit4.stdin_both.pipe",
    "audit4.tab_subject", "audit5.empty_record_subject", "audit5.outfmt.tabular_delim",
    "audit6.option.num_alignments", "audit9.gap_line.query", "audit9.gap_line.subject",
    "audit9.gap_line.underscore", "audit9.title_warning.trailing_space", "digits", "hyphen", "leading_space",
    "tab_after_id", "tab_in_title", "utf8_subject", "x_residue",
)
NOW_SAME_VARIANTS = (
    "audit6.task.rmblastn", "audit2.reward_16bit_0.penalty_0.dc", "audit2.reward_16bit_0.penalty_0.short",
    "audit6.option.no_greedy.zero_gaps.dc", "audit6.option.no_greedy.zero_gaps.short",
    "audit6.option.off_diagonal_range.onehit.dc",
)
NOW_SAME_ERROR = (
    "audit.byte_order_mark", "audit.empty_records_only", "audit12.bom_query", "audit3.evalue_inf",
    "audit3.evalue_nan", "audit3.evalue_space", "audit3.piped_white_space_query", "audit4.hex.no_digits",
    "audit4.leading_blank_subject.penalty_0", "audit4.option.limit_lookup", "audit4.option.max_db_word_count",
    "audit4.option.min_hit_length", "audit4.option.verbose", "audit5.bare_0x.penalty", "audit5.bare_0x.reward",
    "audit5.bare_0x.word_size", "audit5.long_name.out", "e2g.t7.batch_size_text",
)
NOW_SAME_ERROR_VARIANTS = (
    "audit6.task.capitals",
)
E2I_NOW = {
    **{f"{case}{suffix}": "same" for case in NOW_SAME for suffix in ("", ".dc", ".short")},
    **{name: "same" for name in NOW_SAME_VARIANTS},
    **{f"{case}{suffix}": "same-error" for case in NOW_SAME_ERROR for suffix in ("", ".dc", ".short")},
    **{name: "same-error" for name in NOW_SAME_ERROR_VARIANTS},
}
# Approved exception 19 of PD-LOSAT-NCBI-DEFECTS (hit list sizes past Int4): NCBI wraps or crashes at
# these `-max_target_seqs` values and LOSAT keeps every hit, which is NCBI's output at 1073741798. As
# the other harnesses do, NCBI's run takes that size; a `same` result is reported as `exception-19`.
NEAR_SIZE = "1073741798"
E2I_NEAR = {
    f"{case}{suffix}": (word, NEAR_SIZE, "exception-19")
    for case, word in (("audit5.hitlist.wrap", "1073741824"), ("audit5.hitlist.wrap.hex.fmt0", "0x40000000"),
                       ("audit5.hitlist.wrap.fmt7", "2147483597"), ("audit5.hitlist.negative", "2147483598"))
    for suffix in ("", ".dc", ".short")
}
E2I_NOW.update({name: near[2] for name, near in E2I_NEAR.items()})


def expectation(name: str, e2g_expect: str) -> str:
    return E2I_NOW.get(name, E2I_EXPECT.get(name, e2g_expect))


def cases(work: Path) -> list[tuple]:
    rows = []
    for name, argv, expect, *extra in e2g_cases(work):
        rows.append((name, argv, expectation(name, expect), *extra))
        if "-task" not in argv:
            for suffix, task in TASKS:
                if suffix == "short" and name.startswith(SKIP_SHORT):
                    continue
                variant = f"{name}.{suffix}"
                rows.append((variant, [*argv, "-task", task], expectation(variant, expect), *extra))
    return rows


e2c.cases = cases
e2c.NCBI_NEAR.update(E2I_NEAR)

if __name__ == "__main__":
    sys.exit(e2c.main())
