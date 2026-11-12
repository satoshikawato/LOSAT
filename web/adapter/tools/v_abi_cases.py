#!/usr/bin/env python3
"""The V-ABI cases (web/adapter/tests/v_abi.js) as JSON.

Every case is one search: a program, its argv without -outfmt, -out and -num_threads
(the adapter writes every format, the host sets the threads), the working directory,
and for each supported output format the frozen SHA-256 where one exists. A search with
`"outfmt"` adds that -outfmt (a format besides 0, 6 and 7, or 6 or 7 with fields: abi_v2.md
§7) to the argv; the output of it (stream 4) is compared with the CLI's for the same argv, and
its frozen SHA-256 is `frozen["custom"]`. For SAM (`"sam": true`) both sides are compared
without the CL tag of the @PG line, where each writes its own command line, and the frozen
SHA-256 is that of NCBI's fixture without it (`db_search_fixtures.strip_sam_cl`). Outfmt 13 and 14
(`"files": true`) write a ZIP archive on stream 4 (abi_v2.md §5): its entries are compared with the
files of the CLI's run with `-out results`, and the frozen SHA-256 are those of NCBI's master file
with the names of `-out results` (`frozen["custom"]`) and of its side files (`frozen["sides"]`,
joined with `,`).

- full: the regression cases of docs/evidence/losat_web_e1a/capture_outputs.py (Gate A,
  the Stage G matrix, the BLASTP manifest in outfmt 0/6/7), grouped into searches.
  Their inputs use the recorded path spellings, so this needs the Gate A lexical root
  that capture_outputs.py documents; the Stage G searches run as capture_outputs.py runs
  them, with LOSAT/tests/frozen_paths.py's run root in place of the recorded prefix, and
  carry `"frozen_prefix"` (`frozen_paths.case_table_entry`): v_abi.js puts the prefix back
  into an output before comparing it with a frozen hash. It adds the outfmt 0 and 7 fixtures of
  LOSAT/tests/outfmt0_manifest.tsv, whose frozen hashes are NCBI's.
  The rows with -num_descriptions / -num_alignments run as searches of one or two hit list
  sizes (abi_v2.md §7). The group `fasta_input` (session SFc, S10) adds SF fixture inputs of
  LOSAT/tests/fasta_input_fixtures.py that `register` accepts (deflines, sequence lines,
  line ends, text before the first defline, records without residues, ...), FASTA_INPUT_SEARCHES
  per program, with NCBI's frozen hashes, and per program one search whose first query batch
  fails after the outfmt 0 prolog (`"expect": "failure"`: queries without residues only).
- quick: a few searches per program with repository paths, compared with the native
  CLI only (for CI), the NCBI-frozen fixtures of the outfmt 0 manifest in 0, 6 and 7 (with
  fields too) and, of the other formats, QUICK_PER_FORMAT rows of each program and format and
  the archive rows with ranges and titles (the full suite has every row), the NCBI-frozen
  numdesc searches with two hit list sizes, the
  NCBI-frozen BLASTN `-ungapped` searches (four tasks, two input pairs) and three NCBI-frozen
  BLASTN X-drop / `-no_greedy` searches (DI-6).

Both suites add a TBLASTX search over several subjects, which the threaded reactor
reduces in the wasm-threads-only path, and the BLASTN `-window_masker_db` searches
(`window_masker`), whose counts files the host registers by name (abi_v2.md §4 `file_begin`). `--group fasta_input` writes that group of the full
suite alone.

Usage: v_abi_cases.py --suite quick|full [--group fasta_input] --out FILE.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shlex
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "docs/evidence/losat_web_e1a"))
import capture_outputs as capture  # noqa: E402  (the verified case builders)
import frozen_paths  # noqa: E402  (on the path through capture_outputs)

FORMATS = {"blastp": ["0", "6", "7"], "tblastn": ["0", "6", "7"], "blastn": ["0", "6", "7"], "tblastx": ["0", "6", "7"],
           "blastx": ["0", "6", "7"]}
ADAPTER_OWNED = {"-outfmt", "-out", "-num_threads"}
# NCBI's white space around an -outfmt value (`NStr::TruncateSpaces`).
NCBI_SPACE = " \t\n\x0b\x0c\r"


def format_number(outfmt: str) -> int | None:
    """The format number of an -outfmt value (its first word, as NCBI reads it)."""
    number = outfmt.strip(NCBI_SPACE).partition(" ")[0]
    return int(number) if re.fullmatch(r"\+?\d+", number) else None


def custom_outfmt(outfmt: str) -> bool:
    """Whether the -outfmt value adds an output of its own to a run (`run::parse_run`): a format
    besides 0, 6 and 7, or 6 or 7 with words after the number."""
    number = format_number(outfmt)
    has_words = " " in outfmt.strip(NCBI_SPACE)
    return number is not None and number != 0 and (has_words or number not in (6, 7))


# Of the other formats, the quick suite takes this many rows of each program and format, and the
# archive rows whose Bioseqs are the whole records or show the deflines as read.
QUICK_PER_FORMAT = 2
QUICK_ROWS = re.compile(r"^sa\.arch\.\w+\.(qloc|sloc|titles|stitles)\.11$")


def split(command: list[str]) -> tuple[str | None, list[str]]:
    """Returns the -outfmt value and the argv without the adapter-owned options."""
    outfmt, rest, i = None, [], 0
    while i < len(command):
        if command[i] in ADAPTER_OWNED:
            if command[i] == "-outfmt":
                outfmt = command[i + 1]
            i += 2
            continue
        rest.append(command[i])
        i += 1
    return outfmt, rest


def group(cases: list[capture.Case], frozen: bool) -> list[dict]:
    searches: dict[tuple, dict] = {}
    for case in cases:
        if case.program not in FORMATS:
            continue
        command = [frozen_paths.to_run(arg) for arg in case.command]
        outfmt, argv = split(command)
        outfmt = outfmt or "0"
        if outfmt not in FORMATS[case.program]:
            continue  # the other formats and the field lists: the NCBI-frozen rows (outfmt0_fixtures)
        key = (case.program, tuple(argv), str(case.cwd))
        search = searches.setdefault(key, {"program": case.program, "argv": argv, "cwd": str(case.cwd),
                                           "cases": [], "frozen": {}})
        if command != case.command:
            search["frozen_prefix"] = frozen_paths.case_table_entry()
        search["cases"].append(case.case_id)
        if frozen and case.expected_sha256:
            previous = search["frozen"].setdefault(outfmt, case.expected_sha256)
            if previous != case.expected_sha256:
                raise SystemExit(f"{case.case_id}: two frozen hashes for outfmt {outfmt}")
    return list(searches.values())


def outfmt0_fixtures(every_format: bool = True) -> list[dict]:
    """The NCBI-frozen fixtures of the outfmt 0 manifest (docs/evidence/losat_web_e2a/, e2b/), run
    from LOSAT/: one search per argv, with the frozen hash of each fixture's format (none for an
    approved genetic-code deviation, whose oracle is NCBI's database search). A row of another
    format than 0, 6 and 7 (or of 6 or 7 with fields) is a search of its own, which writes that
    output besides; without `every_format`, QUICK_PER_FORMAT rows of each program and format and
    the QUICK_ROWS."""
    sys.path.insert(0, str(REPO / "docs/evidence/losat_web_e2a"))
    import run_oracle  # noqa: E402  (the manifest reader and the oracle's arguments)
    searches: dict[tuple, dict] = {}
    taken: dict[tuple, int] = {}
    for row in run_oracle.read_manifest()[2]:
        if row["program"] not in FORMATS:
            continue
        outfmt = run_oracle.outfmt(row)
        # The rows of outfmt 13 and 14 (the `files` column: side files and a master file next to
        # `-out`) are searches whose output is the ZIP archive of those files.
        files = row.get("files") == "xinclude"
        if row.get("files") and not files:
            continue
        argv = run_oracle.search_argv(row)
        # A hit list of 0 (-num_alignments 0) fails the CLI's outfmt 6 and 7 while the Web run
        # succeeds with the error on stream 3 (abi_v2.md §5, §7), so there is no native output of
        # those formats to compare: the adapter's tests cover the row (run.rs).
        if any(argv[k:k + 2] == ["-num_alignments", "0"] for k in range(len(argv) - 1)):
            continue
        # -parse_deflines works in the CLI only: the adapter refuses it (`run::parse`).
        if "-parse_deflines" in argv:
            continue
        # `register` refuses a `>?` gap line (maintainer decision 3, as in fasta_input): a row
        # whose input has one is a CLI case.
        if any(b">?" in (run_oracle.ENGINE / row[side]).read_bytes() for side in ("query", "subject")):
            continue
        assert argv[-2:] == ["-outfmt", outfmt], row["fixture_id"]
        # A row of another format, or of 6 or 7 with fields, is a search of its own, which writes
        # that output besides.
        custom = custom_outfmt(outfmt)
        sam = run_oracle.is_sam_outfmt(outfmt)
        frozen_kept = True
        if custom and not every_format and " " not in outfmt.strip(NCBI_SPACE):
            kind = (row["program"], format_number(outfmt))
            if taken.get(kind, 0) >= QUICK_PER_FORMAT and not QUICK_ROWS.match(row["fixture_id"]):
                # The search of 0, 6 and 7 alone, as the row was before W11.
                custom = sam = frozen_kept = False
            else:
                taken[kind] = taken.get(kind, 0) + 1
        key = (row["program"], tuple(argv[:-2]), outfmt if custom else None)
        search = searches.setdefault(key, {"program": row["program"], "argv": [row["program"], *argv[:-2]],
                                           "cwd": str(run_oracle.ENGINE), "cases": [], "frozen": {},
                                           **({"outfmt": outfmt} if custom else {}),
                                           **({"sam": True} if sam else {}),
                                           **({"files": True} if files and custom else {})})
        search["cases"].append(row["fixture_id"])
        if frozen_kept and row.get("contract") != run_oracle.DEVIATION:
            frozen = row["stdout_sha256"]
            fixture = run_oracle.ENGINE / "tests/fixtures/outfmt0" / f"{row['fixture_id']}.out"
            if sam:
                frozen = hashlib.sha256(run_oracle.strip_sam_cl(fixture.read_bytes())).hexdigest()
            if files:
                # The master lists the side files by the base of `-out`: `<fixture_id>.out_<n>` in
                # NCBI's run, `results_<n>` in the CLI's and the browser's.
                master = fixture.read_bytes().replace(f"{row['fixture_id']}.out_".encode(), b"results_")
                frozen = hashlib.sha256(master).hexdigest()
                search["frozen"]["sides"] = row["side_sha256"]
            search["frozen"]["custom" if custom else outfmt] = frozen
    return list(searches.values())


FASTA_INPUT_SEARCHES = 32
# Queries without residues only: NCBI's first batch has no query to search, so the run fails
# after the outfmt 0 prolog (on a file sink: AUTHORITY.md of losat_web_e2h, S5 decision 5).
FIRST_BATCH_FAILURE = "@/rec_empty_all.q.{}"


def fasta_input() -> list[dict]:
    """SF fixture inputs through the ABI (gate NOTES open item 3): the rows of
    LOSAT/tests/fasta_input_fixtures.py that NCBI runs on files (no standard input, no
    environment, one thread) and completes in every format of the search, whose inputs
    `register` accepts (no `>?` gap line: maintainer decision 3). One search per program,
    task, inputs and options, with NCBI's frozen hash of each fixture's format; per program
    FASTA_INPUT_SEARCHES of them, taken in turn from each family of inputs (the variant
    name's first word, and the BLASTN task), and the first-batch failure."""
    sys.path.insert(0, str(REPO / "LOSAT/tests"))
    import fasta_input_fixtures as fixtures  # noqa: E402  (the row definitions and frozen hashes)
    grouped: dict[tuple, list] = {}
    for row in fixtures.read_manifest():
        if row.program not in FORMATS or row.kind != "ncbi" or row.stdin or row.env or row.threads != "1":
            continue
        if row.query in ("", "-") or row.subject in ("", "-"):
            continue
        grouped.setdefault((row.program, row.task, row.query, row.subject, row.args), []).append(row)
    failures, families = [], {}
    for (program, task, query, subject, args), rows in sorted(grouped.items()):
        paths = [fixtures.ENGINE / fixtures.rp(query), fixtures.ENGINE / fixtures.rp(subject)]
        if not all(path.is_file() for path in paths) or any(b">?" in path.read_bytes() for path in paths):
            continue
        argv = [program, "-query", fixtures.rp(query), "-subject", fixtures.rp(subject)]
        argv += ["-task", task] if task else []
        search = {"program": program, "argv": argv + shlex.split(args), "cwd": str(fixtures.ENGINE),
                  "cases": [row.row_id for row in rows], "frozen": {}, "group": "fasta_input"}
        rcs = {row.rc for row in rows}
        if query.startswith(FIRST_BATCH_FAILURE.format("")):
            if (rcs == {"3"} and any(row.outfmt == "0" for row in rows) and task in ("", "megablast")
                    and not args and subject == fixtures.base_pair(program)[1]):
                failures.append({**search, "expect": "failure"})
            continue
        if rcs != {"0"}:
            continue
        for row in rows:
            search["frozen"][row.outfmt] = row.stdout_sha256
        family = (program, task, rows[0].row_id.removeprefix(f"{program}.").removeprefix(f"{task}.").split("_")[0])
        families.setdefault(program, {}).setdefault(family, []).append(search)
    searches = []
    for program in FORMATS:
        queues = [list(queue) for _, queue in sorted(families.get(program, {}).items())]
        chosen: list[dict] = []
        while len(chosen) < FASTA_INPUT_SEARCHES and any(queues):
            for queue in queues:
                if queue and len(chosen) < FASTA_INPUT_SEARCHES:
                    chosen.append(queue.pop(0))
        searches += chosen
    return searches + failures


def direct_tblastx() -> dict:
    run = sorted(REPO.glob("docs/evidence/losat_web_e1c/run-*/tblastx-threaded-direct"))[-1]
    rel = run.relative_to(REPO)
    return {"program": "tblastx", "argv": ["tblastx", "-query", str(rel / "q.fna"), "-subject", str(rel / "s.fna")],
            "cwd": str(REPO), "cases": ["tblastx-threaded-direct"], "frozen": {}}


def window_masker() -> list[dict]:
    """BLASTN `-window_masker_db` (line NP WW1): the host registers the counts file that the argv
    names, under the argv's value, before the run (`file_begin`, abi_v2.md §4), and keeps it for the
    next searches of the same name; a value that is not a regular file is not registered, as a file
    that the CLI cannot open. The NCBI-frozen rows of LOSAT/tests/fixtures/blastn_regression/
    manifest.tsv (`wm.*`) without an environment, in their format (0, 6 and 7 compared with the CLI,
    another format as the search's `-outfmt`), and the rows of LOSAT/tests/fixtures/winmask_blastn/
    cases.tsv with `-window_masker_db` and without an environment whose counts are a file (the counts that cannot be read, from
    `make_counts`) or a name of no file, in outfmt 0, 6 and 7: those that fail, with the CLI's error
    and its place among the other errors."""
    sys.path.insert(0, str(REPO / "LOSAT/tests"))
    sys.path.insert(0, str(REPO / "LOSAT/tests/winmask"))
    import blastn_regression_fixtures as regression  # noqa: E402  (the manifest reader)
    import make_blastn_error_fixtures as errors  # noqa: E402  (the unreadable counts files)

    engine = REPO / "LOSAT"
    searches: dict[tuple, dict] = {}

    def add(case_id: str, argv: list[str], frozen: str | None, failure: bool, prolog: bool = True) -> None:
        outfmt, argv = split(argv)
        outfmt = outfmt or "0"
        custom = custom_outfmt(outfmt)
        if failure and outfmt not in FORMATS["blastn"]:
            return  # a failure is compared in the formats of every run
        key = (tuple(argv), outfmt if custom else None)
        search = searches.setdefault(key, {"program": "blastn", "argv": argv, "cwd": str(engine), "cases": [],
                                           "frozen": {}})
        if custom:
            search["outfmt"] = outfmt
        if failure:
            search["expect"] = "failure"
            if not prolog:
                search["prolog"] = False
        search["cases"].append(case_id)
        if frozen:
            search["frozen"]["custom" if custom else outfmt] = frozen

    for row in regression.read_manifest():
        if (row["case_id"].startswith("wm.") and not row["env"] and not row["losat_extra"]
                and "-window_masker_db" in shlex.split(row["argv"])):
            add(row["case_id"], ["blastn", *shlex.split(row["argv"])], row["stdout_sha256"], row["exit"] != "0")
    errors.make_counts()
    for case, env, argv in errors.read_cases():
        if "-window_masker_db" not in argv:
            continue  # -window_masker_taxid (line NP WW2)
        value = argv[argv.index("-window_masker_db") + 1]
        counts = engine / value
        if env or "-db" in argv or (value and counts.exists() and not counts.is_file()) or case == "unreadable":
            continue  # LOSAT Web gives a file's bytes only: no directory, no file it cannot read, no -db
        rc = int((errors.FIXTURES / f"{case}.rc").read_text())
        # Exit 255 and 1 are the errors of the options and their files, before the first query
        # batch: the CLI writes no output. Exit 3 is the first batch's error, after the outfmt 0 prolog.
        add(f"winmask_blastn.{case}", ["blastn", *argv], None, rc != 0, prolog=rc == 3)
        # The counts file's error before the options checked after it (the query range, the
        # expect value of each search's checks), in validate as in run.
        if case == "ascii":
            for option in (["-query_loc", "5-1"], ["-evalue", "0"]):
                at = argv.index("-outfmt")
                add(f"winmask_blastn.ascii{option[0]}", ["blastn", *argv[:at], *option, *argv[at:]], None, True,
                    prolog=False)
    return list(searches.values())


def quick() -> list[dict]:
    cases = []
    blastp = capture.audit_blastp_v010.load_manifest(capture.TESTS / "blastp_v010_parity_manifest.tsv", REPO)
    for base in blastp[:2]:
        _, losat = capture.audit_blastp_v010.build_commands(base, Path("unused"), Path("losat"), Path("{OUTDIR}"))
        cases.append(capture.Case("blastp", base.case_id, capture.strip_binary(losat), REPO))
    for row in capture.compare_blastn_parity.read_manifest(capture.TESTS / "blastn_parity_manifest.tsv"):
        if row["case_id"].startswith(("compact.", "MjPMNV.")):
            absolute = capture.certify._blastn_command_row(REPO, row)
            command = capture.compare_blastn_parity.build_losat_command(absolute, Path("losat"), "{OUT}")
            cases.append(capture.Case("blastn", row["case_id"], capture.strip_binary(command), REPO))
    tblastx = capture.audit_tblastx_v010.load_manifest(capture.TESTS / "tblastx_v010_parity_manifest.tsv", REPO)
    for case in tblastx[1:2]:
        _, losat = capture.audit_tblastx_v010.build_commands(case, Path("unused"), Path("losat"), Path("unused"), Path("{OUT}"))
        cases.append(capture.Case("tblastx", case.case_id, capture.strip_binary(losat), REPO))
    # BLASTX (session SX): the smaller pair of capture_outputs.py's BLASTX_CASES.
    cases += [case for case in capture.blastx_cases() if case.case_id.startswith("small_vs_paja.")]
    prefix = frozen_paths.FROZEN_PREFIX + "/"
    with open(REPO / "docs/evidence/tlosan_stage_g/matrix_162.jsonl") as handle:
        for line in handle:
            row = json.loads(line)
            if row["case"] in {"code1", "code4", "code11", "code32"} and row["threads"] == 1:
                command = [arg.removeprefix(prefix) for arg in row["losat_command"][1:]]
                cases.append(capture.Case("tblastn", row["case"], command, REPO))
    for stem in ("valid_then_skipped", "skipped_then_valid"):
        cases.append(capture.Case("tblastn", stem, [
            "tblastn", "-query", f"docs/evidence/tlosan_stage_g/batch_boundary/{stem}.faa",
            "-subject", "docs/evidence/tlosan_stage_d/all_codes_20260925/fixtures/code1.fna"], REPO))
    searches = group(cases, frozen=False)
    # BLASTP -dbsize (DB-3): an option of the search that LOSAT Web takes with subject records; the
    # effective database length replaces the total length in the statistics.
    first_blastp = next(search for search in searches if search["program"] == "blastp")
    searches.append({**first_blastp, "argv": first_blastp["argv"] + ["-dbsize", "1000000"],
                     "cases": [case + ".dbsize" for case in first_blastp["cases"]], "frozen": {}})
    # BLASTP in blocks (LOSAT_X_BPQBLOCK): a query of about 200,000 letters, four blocks of
    # LOSAT_X_BPQBLOCK's 50,000 letters, with the compressed lookup table of blastp-fast.
    searches.append({"program": "blastp", "argv": ["blastp", "-task", "blastp-fast",
                     "-query", "LOSAT/tests/fasta/NZ_CP006932.faa", "-subject", "LOSAT/tests/fasta/SicyWSV.faa"],
                     "cwd": str(REPO), "cases": ["blastp_fast_query_blocks"], "frozen": {}})
    # Every NCBI-frozen fixture of the outfmt 0 manifest (outfmt0_fixtures), each also compared with the
    # native CLI: among them BLASTN dc-megablast and blastn-short (Session SD), -num_descriptions /
    # -num_alignments (O-webopt), BLASTN -ungapped (O4a-1) and the X-drop and window options (DI-6). From
    # 2026-10-25 an operator precedence in the selection of the -num_descriptions fixtures took them all
    # by accident (`set & {...} | {...}`); the suite keeps them on purpose (line SX, DECISIONS 7).
    searches += outfmt0_fixtures(every_format=False)
    # The frozen outfmt 6 of BLASTN -ungapped on the t2 pair (O4a-1; the four tasks, of the
    # LOSAT/tests/fixtures/blastn_regression/manifest.tsv rows without stderr).
    ungapped = {f"ungapped.{task}.t2.0" for task in ("mb", "bn", "dc", "short")}
    sys.path.insert(0, str(REPO / "LOSAT/tests"))
    import blastn_regression_fixtures as regression  # noqa: E402  (the manifest reader)
    rows = {row["case_id"]: row for row in regression.read_manifest()}
    for search in searches:
        row = next((rows[case.removesuffix(".0")] for case in search["cases"] if case in ungapped), None)
        if row and row["exit"] == "0" and not row["stderr_sha256"]:
            search["frozen"]["6"] = row["stdout_sha256"]
    return searches


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--suite", choices=["quick", "full"], default="quick")
    parser.add_argument("--group", choices=["fasta_input"], help="one group of the full suite alone")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.group and args.suite != "full":
        parser.error("--group selects a group of --suite full")
    if args.group == "fasta_input":
        searches = fasta_input()
    elif args.suite == "full":
        searches = group(capture.all_cases(set(FORMATS)), frozen=True) + outfmt0_fixtures() + fasta_input()
    else:
        searches = quick()
    if not args.group:
        searches.append(direct_tblastx())
        searches += window_masker()
    args.out.write_text(json.dumps(searches, indent=2) + "\n")
    print(f"{len(searches)} searches -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
