#!/usr/bin/env python3
"""Frozen NCBI BLAST+ 2.17.0 TBLASTX outputs for fast pull-request regressions.

The cases cover TBLASTX paths that the Gate A manifest (outfmt 6 of whole genomes) and
the outfmt 0/7 fixtures (LOSAT/tests/outfmt0_manifest.tsv) do not reach: query batches
set by BATCH_SIZE, batches without a valid query, the order of the warnings and the
report (`.merged` cases), hit lists that overflow with ties, ambiguous subjects at other
thread counts, CTOOLKIT_COMPATIBLE, empty inputs and failed writes.

The format is that of blastn_regression_fixtures.py (whose helpers it uses): a case may
set environment variables (`env`) for both programs, every other variable that changes
NCBI's batches or report is unset, `oracle_env` is NCBI's configuration for an approved
exception, and a case id ending in `.merged` runs with standard error merged into
standard output (`2>&1`).

- generate: writes the inputs to LOSAT/tests/fixtures/tblastx_regression/inputs/ (empty
  files and the code4 windows as RNA; the files are committed); the other inputs are those
  of the outfmt 0/7 fixtures.
- freeze --oracle TBLASTX: runs NCBI BLAST+ (comparison oracle only) from LOSAT/ and
  writes <case>.out, <case>.err and the hash columns of manifest.tsv.
- check --losat LOSAT: runs `LOSAT tblastx <argv> <losat_extra>` from LOSAT/ and compares
  stdout, stderr and the exit status with the frozen files.

Usage:
  tblastx_regression_fixtures.py generate
  tblastx_regression_fixtures.py freeze --oracle /path/to/tblastx
  tblastx_regression_fixtures.py check --losat BIN [--jobs N] [--out TSV]
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from blastn_regression_fixtures import FIELDS, REPORT_ENV, run_case, sha256  # noqa: E402

ENGINE = Path(__file__).resolve().parents[1]
FIXTURES = ENGINE / "tests/fixtures/tblastx_regression"
INPUTS = FIXTURES / "inputs"
MANIFEST = FIXTURES / "manifest.tsv"

I = "tests/fixtures/tblastx_regression/inputs"
O = "tests/fasta/outfmt0"
SUBJECTS = f"-subject {O}/tblastx_multi_subject.fasta"
M = f"-query {O}/tblastx_multi_query.fasta {SUBJECTS}"
B = f"-query {O}/tblastx_batch_query.fasta {SUBJECTS}"
V = f"-query {O}/tblastx_invalid_query.fasta {SUBJECTS}"
U = f"-query {O}/tblastx_unsearched_query.fasta {SUBJECTS}"
A = f"-query {O}/tblastx_ambig_query.fasta -subject {O}/tblastx_ambig_subject.fasta"
C = f"-query {O}/tblastx_code4_query.fasta -subject {O}/tblastx_code4_subject.fasta -query_gencode 4"
Y = f"-query {O}/tblastx_many_query.fasta -subject {O}/tblastx_many_subject.fasta"
# (case_id, NCBI argv after `tblastx`, extra LOSAT-only arguments[, environment[, NCBI environment]])
_CASES = [
    # Query batches (blast_input_aux.cpp GetQueryBatchSize; the linking cutoffs use each
    # batch's average query length and smallest Lambda).
    # Queries of 1500, 20000, 800 and 6000 nt: the default batches are (1500 20000) and
    # (800 6000); 700 searches each query alone, 22000 the first three together, 100000 all.
    ("env.batch700", f"{B} -outfmt 6", "", "BATCH_SIZE=700"),
    ("env.batch22000_fmt0", f"{B} -outfmt 0", "", "BATCH_SIZE=22000"),
    ("env.batch100000_fmt7", f"{B} -outfmt 7", "", "BATCH_SIZE=100000"),
    ("env.batch700_threads4", f"{B} -outfmt 6", "-num_threads 4", "BATCH_SIZE=700"),
    ("env.batch_negative", f"{B} -outfmt 6", "", "BATCH_SIZE=-1"),
    # NCBI checks for an empty query before it reads BATCH_SIZE (tblastx_app.cpp:132-137);
    # LOSAT rejects a BATCH_SIZE that is not an integer (NCBI's CStringException), after it.
    ("env.batch_text_empty_query", f"-query {I}/empty.fa {SUBJECTS} -outfmt 6", "", "BATCH_SIZE=abc"),
    # A batch size of 0: the empty first batch fails after the outfmt 0 prolog (exit 3).
    ("env.batch0_fmt0", f"{B} -outfmt 0", "", "BATCH_SIZE=0"),
    ("env.batch0_fmt7", f"{B} -outfmt 7", "", "BATCH_SIZE=0"),
    ("env.batch0_empty_query_fmt0", f"-query {I}/empty.fa {SUBJECTS} -outfmt 0", "", "BATCH_SIZE=0"),
    # A batch of only the all-N query between searched batches (no search, warnings), then a
    # searched batch with the 2-nt query (no warning).
    ("env.batch300_invalid_fmt0.merged", f"{V} -outfmt 0", "", "BATCH_SIZE=300"),
    ("env.batch300_invalid_fmt7.merged", f"{V} -outfmt 7", "", "BATCH_SIZE=300"),
    ("env.batch300_invalid_fmt6.merged", f"{V} -outfmt 6", "", "BATCH_SIZE=300"),
    # Warnings and the report in one stream: the unsearched batch's warnings before its
    # reports, the -max_target_seqs warning before the prolog.
    ("warnings.unsearched_fmt0.merged", f"{U} -outfmt 0", ""),
    ("warnings.unsearched_fmt7.merged", f"{U} -outfmt 7", ""),
    ("warnings.unsearched_fmt6.merged", f"{U} -outfmt 6", ""),
    ("warnings.mts2_fmt0.merged", f"{C} -max_target_seqs 2 -outfmt 0", ""),
    ("warnings.mts4_fmt7.merged", f"{C} -max_target_seqs 4 -outfmt 7", ""),
    # Empty inputs: "Query is Empty!" (tblastx_app.cpp) and an empty subject (exit 3).
    ("empty.query_fmt0", f"-query {I}/empty.fa {SUBJECTS} -outfmt 0", ""),
    ("empty.query_fmt7", f"-query {I}/empty.fa {SUBJECTS} -outfmt 7", ""),
    ("empty.blank_query_fmt6", f"-query {I}/blank.fa {SUBJECTS} -outfmt 6", ""),
    ("empty.subject_fmt0", f"-query {O}/tblastx_code4_query.fasta -subject {I}/empty.fa -outfmt 0", ""),
    ("empty.subject_fmt6", f"-query {O}/tblastx_code4_query.fasta -subject {I}/empty.fa -outfmt 6", ""),
    # Hit lists of 260 subjects: overflow (e-value order, ties by subject order) and threads.
    ("hitlist.max1", f"{Y} -max_target_seqs 1 -outfmt 6", ""),
    ("hitlist.max2_fmt7", f"{Y} -max_target_seqs 2 -outfmt 7", ""),
    ("hitlist.max5", f"{Y} -max_target_seqs 5 -outfmt 6", ""),
    ("hitlist.max250", f"{Y} -max_target_seqs 250 -outfmt 6", ""),
    ("hitlist.default_threads4", f"{Y} -outfmt 6", "-num_threads 4"),
    ("hitlist.max5_fmt0_threads2", f"{Y} -max_target_seqs 5 -outfmt 0", "-num_threads 2"),
    # Ambiguous subjects (random ncbi2na bases in the preliminary search) at other thread counts.
    ("ambig.fmt6", f"{A} -outfmt 6", ""),
    ("ambig.fmt6_threads4", f"{A} -outfmt 6", "-num_threads 4"),
    ("ambig.fmt0_threads2", f"{A} -outfmt 0", "-num_threads 2"),
    # Word and window options in the other formats.
    ("options.thrwin_fmt6", f"{M} -threshold 12 -window_size 10 -outfmt 6", ""),
    ("options.thrwin_fmt7", f"{M} -threshold 16 -window_size 60 -outfmt 7", ""),
    ("options.segno_fmt7", f"{C} -seg no -outfmt 7", ""),
    # showdefline.cpp kBits is "(bits)" when CTOOLKIT_COMPATIBLE is set (also empty); the
    # tabular formats do not read it.
    ("ctoolkit.fmt0", f"{C} -outfmt 0", "", "CTOOLKIT_COMPATIBLE=1"),
    ("ctoolkit.empty_fmt0", f"{C} -max_target_seqs 3 -outfmt 0", "", "CTOOLKIT_COMPATIBLE="),
    ("ctoolkit.fmt7", f"{C} -outfmt 7", "", "CTOOLKIT_COMPATIBLE=1"),
    # A .ncbirc (found through $HOME) with entries that change no output.
    ("ncbirc.harmless_fmt0", f"{C} -outfmt 0", "",
     "HOME=tests/fixtures/blastn_regression/ncbirc_home BLAST_USAGE_REPORT=0"),
    # U is read as T (CFastaReader keeps U; the search and the display translate it as T,
    # on both strands): the windows of the code4 fixtures with every T as U.
    ("input.rna_query", f"-query {I}/rna_query.fa -subject {O}/tblastx_code4_subject.fasta -query_gencode 4"
                        " -outfmt 6", ""),
    ("input.rna_subject_fmt0", f"-query {O}/tblastx_code4_query.fasta -subject {I}/rna_subject.fa"
                               " -query_gencode 4 -outfmt 0", ""),
    ("input.rna_both_fmt7", f"-query {I}/rna_query.fa -subject {I}/rna_subject.fa -query_gencode 4 -outfmt 7", ""),
    # Sum statistics link the HSPs per query and strand (link_hsps.c context/3): a 3- or
    # 4-nt query before another query in the batch.
    ("batch.short3_first", f"-query {I}/short3_query.fa {SUBJECTS} -outfmt 6", ""),
    ("batch.short3_first_fmt0", f"-query {I}/short3_query.fa {SUBJECTS} -max_target_seqs 1 -outfmt 0", ""),
    ("batch.short4_between_fmt7", f"-query {I}/short4_query.fa {SUBJECTS} -outfmt 7", ""),
    # Ties of HSPs with equal scores in different query frames, found in S08 (the inputs
    # come from the inventory's and the investigations' runs, not from `generate`): the
    # init hit list sorted by score_compare_match per subject chunk (aa_ungapped.c:234-235;
    # tie_frames: a cut of LC738874 20001-100000 against itself, tie_init), the HSP list
    # sorted by score after the first BLAST_LinkHsps (link_hsps.c:1802-1803; tie_nisland:
    # an N island in the subject, tie_trim: a two-hit tail trimmed by the re-evaluation),
    # and BLAST_LargeGapSumE in NCBI's order of evaluation (blast_stat.c:4560-4561, sume).
    ("tie.frames", f"-query {I}/tie_frames_query.fa -subject {I}/tie_frames_subject.fa -outfmt 6", ""),
    ("tie.init_seg_no", f"-query {I}/tie_init_query.fa -subject {I}/tie_init_subject.fa -seg no -outfmt 6", ""),
    ("tie.nisland_seg_no", f"-query {I}/tie_nisland_query.fa -subject {I}/tie_nisland_subject.fa -seg no"
                           " -outfmt 6", ""),
    ("tie.trim_seg_no_fmt0", f"-query {I}/tie_trim_query.fa -subject {I}/tie_trim_subject.fa -seg no -outfmt 0",
     ""),
    ("sume.large_gap_seg_no", f"-query {I}/sume_query.fa -subject {I}/sume_subject.fa -seg no -outfmt 6", ""),
    # SEG keeps only the head of the segments of a left recursion (blast_seg.c:2086-2101):
    # a low-complexity run in frame +1 (the inputs come from the S08 investigation).
    ("seg.left_recursion", f"-query {I}/seg_query.fa -subject {I}/seg_subject.fa -outfmt 6", ""),
    ("seg.left_recursion_fmt0", f"-query {I}/seg_query.fa -subject {I}/seg_subject.fa -outfmt 0", ""),
    # A failed outfmt 0 write ("BLAST failed to write output", exit 6; Linux /dev/full),
    # also with the warnings of an unsearched batch (the stream fails before the query is read).
    ("write.devfull_fmt0", f"{C} -outfmt 0 -out /dev/full", ""),
    ("write.devfull_unsearched_fmt0", f"{U} -outfmt 0 -out /dev/full", ""),
    # NCBI's check of the hit saving options (blast_options.c:1518-1523): an -evalue of 0 (or
    # one that reads as 0) fails after the formatting warning and before `Query is Empty!`.
    ("options.evalue0_fmt0", f"{C} -evalue 0 -outfmt 0", ""),
    ("options.evalue0_mts1_fmt6.merged", f"{C} -evalue 1e-400 -max_target_seqs 1 -outfmt 6", ""),
    ("options.evalue0_empty_query_fmt7", f"-query {I}/empty.fa {SUBJECTS} -evalue 0 -outfmt 7", ""),
    # NCBI reads a blank line before the first defline of the subjects without a message,
    # so an empty query still gives `Query is Empty!`.
    ("input.subject_blank_first_empty_query", f"-query {I}/empty.fa -subject {I}/blank_first_subject.fa"
                                              " -outfmt 0", ""),
    # Subject titles that NCBI's HtmlDecode leaves as they are (not a name of its table, a
    # final `;` trimmed before the decoding, `&xi;` read from the `i`), and titles that NCBI
    # decodes or reads past on subjects without hits (NCBI makes only the shown titles).
    ("title.kept_fmt0", f"-query {O}/tblastx_many_query.fasta -subject {I}/titles_kept_subject.fa"
                        " -outfmt 0", ""),
    ("title.hitless_fmt0", f"-query {O}/tblastx_many_query.fasta -subject {I}/titles_hitless_subject.fa"
                           " -outfmt 0", ""),
    # A DEL (0x7f) in the deflines: NCBI's title ends at the first byte below a space only
    # (fasta_reader_utils.cpp:215-225).
    ("input.del_deflines_fmt0", f"-query {I}/del_query.fa -subject {I}/del_subject.fa -query_gencode 4"
                                " -outfmt 0", ""),
    ("input.del_deflines_fmt7", f"-query {I}/del_query.fa -subject {I}/del_subject.fa -query_gencode 4"
                                " -outfmt 7", ""),
    # The query splitter of every batch reads CHUNK_SIZE (a size_t that must be divisible by
    # 3 for a translated query; -1 is 2^64 - 1) and OVERLAP_CHUNK_SIZE, but does not split
    # an ungapped search (split_query_cxx.cpp:55-61, local_blast.cpp:98-103).
    ("env.chunk10_fmt0", f"{C} -outfmt 0", "", "CHUNK_SIZE=10"),
    ("env.chunk_minus2_fmt6", f"{C} -outfmt 6", "", "CHUNK_SIZE=-2"),
    ("env.chunk9_fmt7", f"{C} -outfmt 7", "", "CHUNK_SIZE=9"),
    ("env.chunk_minus1_fmt6", f"{C} -outfmt 6", "", "CHUNK_SIZE=-1"),
    ("env.overlap_negative_fmt6", f"{C} -outfmt 6", "", "OVERLAP_CHUNK_SIZE=-5"),
    # A SEG window, locut or hicut that is not above 0 keeps NCBI's default (blast_filter.c
    # 1147-1154), before SEG's own check of the parameters.
    ("options.seg_locut0_fmt0", f"-query {I}/seg_query.fa -subject {I}/seg_subject.fa -seg '12 0 2.5'"
                                " -outfmt 0", ""),
    ("options.seg_nonpositive_fmt6", f"-query {I}/seg_query.fa -subject {I}/seg_subject.fa -seg '0 -1 0'"
                                     " -outfmt 6", ""),
    ("options.seg_window_negative_fmt7", f"-query {I}/seg_query.fa -subject {I}/seg_subject.fa"
                                         " -seg '-5 2.2 2.5' -outfmt 7", ""),
    # BLAST_Cutoffs returns at least the caller's 1 (blast_stat.c:4097, 4126-4129): a large
    # -evalue against the small search space of a 150-nt query (a copy of LC738875 at
    # 76964; audit (b) F-1) gives a cutoff of 1, not a negative one, to the linking.
    ("cutoff.floor_evalue_1e10", f"-query {I}/cutoff_q150.fa -subject tests/fasta/LC738875.fasta"
                                 " -evalue 1e10 -outfmt 6", ""),
]
CASES = [(*case, *[""] * (5 - len(case))) for case in _CASES]


def rna(source: Path, target: Path) -> None:
    """The records of `source` with every T read as U (u for t), deflines unchanged."""
    lines = source.read_text().splitlines()
    target.write_text("".join((line if line.startswith(">") else line.replace("T", "U").replace("t", "u")) + "\n"
                              for line in lines))


def command_generate(_args) -> int:
    INPUTS.mkdir(parents=True, exist_ok=True)
    (INPUTS / "empty.fa").write_bytes(b"")
    (INPUTS / "blank.fa").write_bytes(b"\n  \n\t\n")
    rna(ENGINE / O / "tblastx_code4_query.fasta", INPUTS / "rna_query.fa")
    rna(ENGINE / O / "tblastx_code4_subject.fasta", INPUTS / "rna_subject.fa")
    # Queries of 3 and 4 nt (2 and 4 frames with a residue; NCBI still has 6 contexts per
    # query) before and between the multi queries, in the same batch.
    multi = (ENGINE / O / "tblastx_multi_query.fasta").read_text()
    records = multi.split(">")[1:]
    (INPUTS / "short3_query.fa").write_text(">s3 three\nACG\n" + ">" + ">".join(records))
    (INPUTS / "short4_query.fa").write_text(">" + records[1] + ">s4 four\nACGT\n>" + records[2]
                                            + ">sN3 three N\nNNN\n>" + records[3])
    code4_query = (ENGINE / O / "tblastx_code4_query.fasta").read_text()
    code4_subject = (ENGINE / O / "tblastx_code4_subject.fasta").read_text()
    (INPUTS / "blank_first_subject.fa").write_text("\n" + code4_subject)
    # The first records of the `many` subjects (each with hits of the `many` query) with
    # new titles, and subjects of N only (no hits).
    many = [chunk.split("\n", 1)[1] for chunk in (ENGINE / O / "tblastx_many_subject.fasta").read_text().split(">")[1:6]]
    kept = ("R&D; x", "a&foo;b c", "a&amp;", "s &xi;t", "q &X41; r")
    (INPUTS / "titles_kept_subject.fa").write_text("".join(f">{title}\n{seq}" for title, seq in zip(kept, many)))
    unknown = "N" * 60 + "\n"
    (INPUTS / "titles_hitless_subject.fa").write_text(
        f">hit one\n{many[0]}>, ,\n{unknown * 15}>s &amp; t\n{unknown * 15}>x &#38; y\n{unknown * 15}")
    (INPUTS / "del_query.fa").write_text(code4_query.replace(">", ">q\x7fid del\x7f ", 1))
    (INPUTS / "del_subject.fa").write_text(code4_subject.replace(">", ">s\x7fid del\x7f ", 1))
    return 0


def read_manifest() -> list[dict[str, str]]:
    lines = [line for line in MANIFEST.read_text().splitlines() if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def command_freeze(args) -> int:
    if any(key in os.environ for key in REPORT_ENV) or (Path.home() / ".ncbirc").exists():
        raise SystemExit(f"unset {REPORT_ENV} and remove ~/.ncbirc before freezing")
    oracle = str(Path(args.oracle).resolve())
    version = subprocess.run([oracle, "-version"], capture_output=True, text=True).stdout.strip().replace("\n", "; ")
    rows = []
    for case_id, argv, extra, case_env, oracle_env in CASES:
        result = run_case(oracle, [], case_id, argv, oracle_env or case_env)
        (FIXTURES / f"{case_id}.out").write_bytes(result.stdout)
        err = FIXTURES / f"{case_id}.err"
        if result.stderr:
            err.write_bytes(result.stderr)
        elif err.exists():
            err.unlink()
        rows.append({"case_id": case_id, "argv": argv, "losat_extra": extra, "exit": str(result.returncode),
                     "stdout_sha256": sha256(result.stdout), "stdout_bytes": str(len(result.stdout)),
                     "stderr_sha256": sha256(result.stderr) if result.stderr else "", "env": case_env,
                     "oracle_env": oracle_env})
        print(f"{case_id}\texit {result.returncode}\t{result.stdout.count(b'\n')} lines\t{len(result.stderr)} stderr bytes", flush=True)
    with open(MANIFEST, "w", newline="") as handle:
        handle.write(f"# Frozen TBLASTX outputs of {version} (comparison oracle only; see tblastx_regression_fixtures.py).\n")
        handle.write("# Run from LOSAT/: env <env> tblastx <argv> (LOSAT: LOSAT tblastx <argv> <losat_extra>). Written by `freeze`; do not edit.\n")
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return 0


def check_one(losat: str, row: dict[str, str]) -> tuple[str, str]:
    result = run_case(losat, ["tblastx"], row["case_id"], row["argv"], row.get("env") or "", row["losat_extra"])
    expected_out = (FIXTURES / f"{row['case_id']}.out").read_bytes()
    err_path = FIXTURES / f"{row['case_id']}.err"
    expected_err = err_path.read_bytes() if err_path.exists() else b""
    problems = []
    if sha256(expected_out) != row["stdout_sha256"]:
        problems.append("frozen stdout file does not match the manifest")
    if str(result.returncode) != row["exit"]:
        problems.append(f"exit {result.returncode}, expected {row['exit']}")
    if result.stdout != expected_out:
        expected_lines, actual_lines = expected_out.split(b"\n"), result.stdout.split(b"\n")
        line = next((n for n, (a, b) in enumerate(zip(expected_lines, actual_lines), 1) if a != b),
                    min(len(expected_lines), len(actual_lines)))
        problems.append(f"stdout differs (expected {len(expected_lines) - 1} lines, got {len(actual_lines) - 1}; first difference at line {line})")
    if result.stderr != expected_err:
        problems.append(f"stderr differs: {result.stderr[:200]!r}")
    return row["case_id"], "; ".join(problems) or "same"


def command_check(args) -> int:
    losat = str(Path(args.losat).resolve())
    rows = read_manifest()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(lambda row: check_one(losat, row), rows))
    lines = [f"{case_id}\t{verdict}" for case_id, verdict in results]
    if args.out:
        Path(args.out).write_text("case_id\tresult\n" + "\n".join(lines) + "\n")
    print("\n".join(lines))
    differing = [case_id for case_id, verdict in results if verdict != "same"]
    print(f"{len(results)} cases, {len(differing)} differ")
    return 1 if differing else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("generate")
    freeze = sub.add_parser("freeze")
    freeze.add_argument("--oracle", required=True)
    check = sub.add_parser("check")
    check.add_argument("--losat", required=True)
    check.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    check.add_argument("--out")
    args = parser.parse_args()
    return {"generate": command_generate, "freeze": command_freeze, "check": command_check}[args.action](args)


if __name__ == "__main__":
    sys.exit(main())
