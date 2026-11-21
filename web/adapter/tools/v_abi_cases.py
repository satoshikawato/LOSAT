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
reduces in the wasm-threads-only path, the BLASTN `-window_masker_db` and `-window_masker_taxid`
searches (`window_masker`), whose counts files the host registers by name (abi_v2.md §4 `file_begin`),
the runs whose outcome and standard error are those of the argv's own `-outfmt` (`own_search`), and the
inputs of line WQA's D13 and D16 (`scan_d13`: gap lines, joined lines, Seq-id first lines, and
`-conffile` with `[BLAST] DATA_LOADERS=none`, whose registry file the host registers as kind 2).
`--group fasta_input` writes that group of the full suite alone.

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
        # A PSSM query (tblastn -in_pssm, PSI-TBLASTN; the query cell `pssm:` or `pssm-gz:` names
        # the PSSM file) is not a Web input yet (line PS, item PSW1): the adapter refuses it.
        if row["query"].startswith(("pssm:", "pssm-gz:")):
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
        # -parse_deflines runs in LOSAT Web since line WQA's sweep, D7, 2026-11-19 (the host reads the
        # inputs together, `losat_web2_register_inputs`). The pd.self rows other than outfmt 0 are a
        # record that is both a query and a subject and has masks, whose outfmt 0 LOSAT refuses, that
        # format's alone (`believed::refuse_shared_masks_in`): the CLI's run of outfmt 0 fails and the
        # row's format runs, so each row is a search of its own whose argv names its format (6 and 7
        # too), with outfmt 0 failing alone (`"fails_alone"`, web/adapter/tests/v_abi.js).
        fails_alone = "-parse_deflines" in argv and row["fixture_id"].startswith("pd.self.") \
            and not row["fixture_id"].endswith(".0")
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
        own = custom or fails_alone
        key = (row["program"], tuple(argv[:-2]), outfmt if own else None)
        search = searches.setdefault(key, {"program": row["program"], "argv": [row["program"], *argv[:-2]],
                                           "cwd": str(run_oracle.ENGINE), "cases": [], "frozen": {},
                                           **({"outfmt": outfmt} if own else {}),
                                           **({"fails_alone": [0]} if fails_alone else {}),
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
    environment, one thread) and completes in every format of the search (a `>?` gap line too since
    line WQA's D13, 2026-11-20, which lifted maintainer decision 3's refusal). One search per program,
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
        if not all(path.is_file() for path in paths):
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
    """BLASTN `-window_masker_db` (line NP WW1) and `-window_masker_taxid` (WW2): the host registers
    the counts file that the argv names, under the argv's value, before the run (`file_begin`,
    abi_v2.md §4), and keeps it for the next searches of the same name; a value that is not a regular
    file is not registered, as a file that the CLI cannot open. A search with `"window_masker_path"`
    (the CLI's `WINDOW_MASKER_PATH`) has the host register every readable regular file
    `<path>/<n>/wmasker.obinary` as `windowmasker/<n>/wmasker.obinary`, LOSAT Web's directory of
    `-window_masker_taxid` (a directory under that name is not registered: the CLI's path that is no
    file). The NCBI-frozen rows of LOSAT/tests/fixtures/blastn_regression/manifest.tsv (`wm.*`)
    without an environment, in their format (0, 6 and 7 compared with the CLI, another format as the
    search's `-outfmt`), and the rows of LOSAT/tests/fixtures/winmask_blastn/cases.tsv in outfmt 0, 6
    and 7, with those that fail with the CLI's error and its place among the other errors: with
    `-window_masker_db` whose counts are a file (the counts that cannot be read, from `make_counts`) or
    a name of no file, and with `-window_masker_taxid` whose environment is `WINDOW_MASKER_PATH`
    alone, or none with the taxid 0. LOSAT Web has no other environment (`CHUNK_SIZE`), no registry
    (`.ncbirc`, `NCBI_CONFIG*`), no current directory for the taxid (its directory is always
    `windowmasker`), no file it cannot read (`taxid_unreadable`, `unreadable`), no directory as a
    counts file, and no `-db`. Two searches of its own (`winmask_web.taxid_unknown_*`): a taxid whose
    file has no known format; the first again with `-num_descriptions 5`, a run of two searches
    (abi_v2.md §7) that fails, whose stream 3 and error are the CLI's standard error for the argv's
    own outfmt 0, without the tabular search's `-num_descriptions` warning (line WQA's sweep, D3
    review, 2026-11-19)."""
    sys.path.insert(0, str(REPO / "LOSAT/tests"))
    sys.path.insert(0, str(REPO / "LOSAT/tests/winmask"))
    import blastn_regression_fixtures as regression  # noqa: E402  (the manifest reader)
    import make_blastn_error_fixtures as errors  # noqa: E402  (the unreadable counts files)

    engine = REPO / "LOSAT"
    searches: dict[tuple, dict] = {}

    def add(case_id: str, argv: list[str], frozen: str | None, failure: bool, prolog: bool = True,
            window_masker_path: str | None = None) -> None:
        outfmt, argv = split(argv)
        outfmt = outfmt or "0"
        custom = custom_outfmt(outfmt)
        if failure and outfmt not in FORMATS["blastn"]:
            return  # a failure is compared in the formats of every run
        key = (tuple(argv), outfmt if custom else None, window_masker_path)
        search = searches.setdefault(key, {"program": "blastn", "argv": argv, "cwd": str(engine), "cases": [],
                                           "frozen": {}})
        if custom:
            search["outfmt"] = outfmt
        if window_masker_path is not None:
            search["window_masker_path"] = window_masker_path
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
        if "-db" in argv or case in ("unreadable", "taxid_unreadable"):
            continue  # LOSAT Web has no -db and no file it cannot read
        variables = dict(item.split("=", 1) for item in env.split())
        directory = None
        if "-window_masker_taxid" in argv:
            taxid = argv[argv.index("-window_masker_taxid") + 1]
            if set(variables) - {"WINDOW_MASKER_PATH"}:
                continue  # no other environment, no registry
            directory = variables.get("WINDOW_MASKER_PATH")
            if not directory and taxid != "0":
                continue  # the CLI's current directory: LOSAT Web's directory is always windowmasker
        elif variables:
            continue
        if "-window_masker_db" in argv:
            value = argv[argv.index("-window_masker_db") + 1]
            counts = engine / value
            if value and counts.exists() and not counts.is_file():
                continue  # LOSAT Web gives a file's bytes only: no directory
        rc = int((errors.FIXTURES / f"{case}.rc").read_text())
        # Exit 255 and 1 are the errors of the options and their files, before the first query
        # batch: the CLI writes no output. Exit 3 is the first batch's error, after the outfmt 0 prolog.
        add(f"winmask_blastn.{case}", ["blastn", *argv], None, rc != 0, prolog=rc == 3,
            window_masker_path=directory)
        # The counts file's error before the options checked after it (the query range, the
        # expect value of each search's checks), in validate as in run.
        if case == "ascii":
            for option in (["-query_loc", "5-1"], ["-evalue", "0"]):
                at = argv.index("-outfmt")
                add(f"winmask_blastn.ascii{option[0]}", ["blastn", *argv[:at], *option, *argv[at:]], None, True,
                    prolog=False)
    # The file of a taxid is read without NCBI's check of the format (`WindowMaskerTaxidToDb`): a file
    # of no known format (empty, 3 bytes) is the reader's "unrecognized unit counts format" in each
    # query batch (exit 3), not the refusal of `-window_masker_db`'s check (line NP WW2 review F2).
    unknown = "target/winmask_blastn/taxid_unknown"
    for taxid, content in (("11", b""), ("12", b"abc")):
        (engine / unknown / taxid).mkdir(parents=True, exist_ok=True)
        (engine / unknown / taxid / "wmasker.obinary").write_bytes(content)
        add(f"winmask_web.taxid_unknown_{taxid}",
            ["blastn", "-query", "tests/fixtures/blastn_regression/inputs/wm_query.fa",
             "-subject", "tests/fixtures/blastn_regression/inputs/s60k.fa", "-task", "blastn",
             "-window_masker_taxid", taxid, "-outfmt", "6"], None, True, prolog=True, window_masker_path=unknown)
    add("winmask_web.taxid_unknown_11.num_descriptions",
        ["blastn", "-query", "tests/fixtures/blastn_regression/inputs/wm_query.fa",
         "-subject", "tests/fixtures/blastn_regression/inputs/s60k.fa", "-task", "blastn",
         "-window_masker_taxid", "11", "-num_descriptions", "5", "-outfmt", "6"], None, True, prolog=True,
        window_masker_path=unknown)
    return list(searches.values())


def own_search() -> list[dict]:
    """Runs whose outcome and standard error are those of the argv's own `-outfmt`, as the CLI's
    run of that argv (line WQA's sweep, D3 reviews, 2026-11-19).

    Runs of two searches (`-num_descriptions`, `-num_alignments`; abi_v2.md §7) that fail:
    `-num_alignments 0` with `-outfmt 6` in every program (NCBI's `No hits are being saved` after
    the formatting handler's warning, where the CLI's outfmt 0 completes), also with `-evalue 0`,
    whose error the CLI gives after the hit list's (`validate` too); and two searches that both
    fail, with and without `-outfmt 6`, whose error and stream 3 are those of the argv's own
    `-outfmt` (BLASTX's `-num_alignments 0 -evalue 0`). Runs of two searches that complete where
    the other search fails, whose formats are empty and whose warnings and error stream 3 adds
    after the argv's (review r2 L1): `-num_alignments 0` without `-outfmt` in every program, also
    with BLASTN's `-num_descriptions 2147483647`, whose format 0 keeps every hit (LOSAT computes the
    preliminary hit list size without NCBI's 32-bit overflow, PD-LOSAT-NCBI-DEFECTS); BLASTN's
    `-num_descriptions 2147483647` with `-outfmt 6` fails as the argv's `-num_alignments 0`, and
    without `-num_alignments 0` completes with every format, with and without `-outfmt 6`.

    The formatting handler's warnings of the formats that the argv does not ask for (review r2
    M1): runs that fail with `-line_length`, `-sorthits`, `-sorthsps` or equal `-num_descriptions`
    and `-num_alignments` and an `-outfmt` of 0 (or none) to 4, whose standard error has none of the
    tabular formats' lines (subjects without residues: NCBI's `The average subject length is too
    short`; queries without residues: the first batch's `Sequence contains no data`), and runs that
    complete, whose stream 3 adds those lines after the argv's standard error (after a subject's
    warning where the argv asks for outfmt 0 or 1, before it for outfmt 6).

    The application start's error without a line end (review r2 L2): `-conffile` of a file that is
    not registered (NCBI's `Critical: (106.23) Registry file "x.ini" not found, aborting.`).

    Runs whose own `-outfmt` is structured (5, 8, 11, 13; BLASTN's 17), with a query whose first
    record has no residues (review r3 M2): the CLI's run of the structured format does not write
    the query's `Sequence contains no data` warning, which those of outfmt 0, 6 and 7 write, so
    stream 3 adds it after the argv's standard error; once, also in a run of two searches
    (`-outfmt 5 -num_descriptions 5`), whose outfmt 0 writes it in the other search (review of
    line WQA-LW, r2 L3, 2026-11-20).

    TBLASTX runs whose messages go nowhere (`-logfile /dev/null`, `-logfile ""`; review r3 M1):
    stream 3 is empty, as the CLI's standard error for every format, also where the other formats
    would add a query's warning to it (`-outfmt 5`, the query with an empty first record)."""
    inputs = {
        "blastp": ["-query", "tests/fasta/outfmt0/e2e_many_query.faa",
                   "-subject", "tests/fasta/outfmt0/e2e_many_subject.faa"],
        "tblastn": ["-query", "tests/fasta/outfmt0/e2e_many_query.faa",
                    "-subject", "tests/fasta/outfmt0/e2e_many_subject.fna"],
        "blastn": ["-query", "tests/fasta/outfmt0/multi_query.fasta", "-subject", "tests/fasta/outfmt0/multi_subject.fasta",
                   "-task", "blastn", "-word_size", "7", "-evalue", "1000"],
        "tblastx": ["-query", "tests/fasta/outfmt0/tblastx_multi_query.fasta",
                    "-subject", "tests/fasta/outfmt0/tblastx_multi_subject.fasta"],
        "blastx": ["-query", "tests/db_search/bslice_q.fna", "-subject", "tests/db_search/bslice.faa"],
    }
    empty = "tests/fixtures/fasta_input/rec_empty_all"
    # Subjects without residues (BLASTN, TBLASTX) and queries without residues (the others).
    without_residues = {
        "blastn": ["-query", "tests/fasta/outfmt0/multi_query.fasta", "-subject", f"{empty}.s.nuc.fa"],
        "tblastx": ["-query", "tests/fasta/outfmt0/tblastx_multi_query.fasta", "-subject", f"{empty}.s.nuc.fa"],
        "blastp": ["-query", f"{empty}.q.prot.faa", "-subject", "tests/fasta/outfmt0/e2e_many_subject.faa"],
        "tblastn": ["-query", f"{empty}.q.prot.faa", "-subject", "tests/fasta/outfmt0/e2e_many_subject.fna"],
        "blastx": ["-query", f"{empty}.q.nuc.fa", "-subject", "tests/db_search/bslice.faa"],
    }
    # Subjects with a record without residues among others: NCBI's `Subject sequence contains no
    # data` warning in a run that completes.
    nucleotide, protein = "tests/fasta/outfmt11/arch_tbx_s_emptyrec.fna", "tests/fasta/outfmt11/arch_bp_emptysubj.faa"
    warned = {program: [*argv[:argv.index("-subject") + 1], protein if program in ("blastp", "blastx") else nucleotide,
                        *argv[argv.index("-subject") + 2:]]
              for program, argv in inputs.items()}
    searches = []

    def add(case_id: str, program: str, options: list[str], outfmt: str | None = None, prolog: bool = True,
            failure: bool = True, given: dict | None = None) -> None:
        searches.append({"program": program, "argv": [program, *(given or inputs)[program], *options],
                         "cwd": str(REPO / "LOSAT"), "cases": [case_id], "frozen": {},
                         **({"expect": "failure"} if failure else {}),
                         **({"outfmt": outfmt} if outfmt else {}), **({} if prolog else {"prolog": False})})

    for program in FORMATS:
        add(f"own_search.{program}.na0.outfmt6", program, ["-num_alignments", "0"], "6")
        add(f"own_search.{program}.na0", program, ["-num_alignments", "0"], failure=False)
    overflow = ["-num_descriptions", "2147483647"]
    add("own_search.blastn.nd_overflow", "blastn", overflow, failure=False)
    add("own_search.blastn.nd_overflow.outfmt6", "blastn", overflow, "6", failure=False)
    add("own_search.blastn.nd_overflow.na0.outfmt0", "blastn", [*overflow, "-num_alignments", "0"], failure=False)
    add("own_search.blastn.nd_overflow.na0.outfmt6", "blastn", [*overflow, "-num_alignments", "0"], "6")
    for outfmt, name in ((None, "outfmt0"), ("6", "outfmt6")):
        for program in ("blastx", "blastp", "tblastn", "tblastx"):
            add(f"own_search.{program}.na0.evalue0.{name}", program, ["-num_alignments", "0", "-evalue", "0"],
                outfmt, prolog=False)
    handler = {
        "line_length": ["-line_length", "80"],
        "sorthits": ["-sorthits", "1"],
        "sorthsps": ["-sorthsps", "1"],
        "nd5_na5": ["-num_descriptions", "5", "-num_alignments", "5"],
    }
    for program in FORMATS:
        for name, options in handler.items():
            add(f"handler.{program}.{name}.no_residues", program, options, given=without_residues)
        add(f"handler.{program}.sorthsps.no_residues.outfmt1", program, handler["sorthsps"], "1",
            given=without_residues)
        add(f"handler.{program}.line_length.sorthsps", program, [*handler["line_length"], *handler["sorthsps"]],
            failure=False)
        add(f"handler.{program}.line_length.sorthsps.outfmt1", program,
            [*handler["line_length"], *handler["sorthsps"]], "1", failure=False)
        for outfmt in (None, "1", "6"):
            add(f"handler.{program}.line_length.sorthsps.warned{'.outfmt' + outfmt if outfmt else ''}", program,
                [*handler["line_length"], *handler["sorthsps"]], outfmt, failure=False, given=warned)
    add("start.blastp.conffile_missing", "blastp", ["-conffile", "x.ini"], prolog=False)
    arch = "tests/fasta/outfmt11"
    empty_first = {
        "blastn": ["-query", f"{arch}/arch_tbx_emptyfirst.fna", "-subject", f"{arch}/arch_tbx_s.fna"],
        "tblastx": ["-query", f"{arch}/arch_tbx_emptyfirst.fna", "-subject", f"{arch}/arch_tbx_s.fna"],
        "blastp": ["-query", f"{arch}/arch_bp_emptyfirst.faa", "-subject", f"{arch}/arch_bp_s.faa"],
        "tblastn": ["-query", f"{arch}/arch_tbn_emptyfirst.faa", "-subject", f"{arch}/arch_tbn_s.fna"],
        "blastx": ["-query", f"{arch}/arch_bx_emptyfirst.fna", "-subject", f"{arch}/arch_bx_s.faa"],
    }
    for program in FORMATS:
        for outfmt in ("5", "8", "11", "13", *(("17",) if program == "blastn" else ())):
            add(f"query_warnings.{program}.outfmt{outfmt}", program, [], outfmt, failure=False, given=empty_first)
            searches[-1].update({"files": True} if outfmt == "13" else {"sam": True} if outfmt == "17" else {})
        # A run of two searches (review r2 of 2026-11-20, L3): outfmt 6 and 7 of the argv's search
        # add the warning, and outfmt 0's search, which writes it too, does not add it again.
        add(f"query_warnings.{program}.outfmt5.nd5", program, ["-num_descriptions", "5"], "5", failure=False,
            given=empty_first)
    for name, logfile in (("devnull", "/dev/null"), ("empty", "")):
        add(f"query_warnings.tblastx.outfmt5.logfile.{name}", "tblastx", ["-logfile", logfile], "5", failure=False,
            given=empty_first)
    for name, logfile in (("devnull", "/dev/null"), ("empty", "")):
        add(f"logfile.tblastx.na0.{name}", "tblastx", ["-num_alignments", "0", "-logfile", logfile], failure=False)
        add(f"logfile.tblastx.line_length.warned.{name}", "tblastx", [*handler["line_length"], "-logfile", logfile],
            failure=False, given=warned)
    return searches


# Line WQA's D13 and D16 (2026-11-20): rows of the SF fixtures that NCBI ran (its frozen hashes), by program, query
# and subject: `>?` gap lines (a), which LOSAT Web refused where the CLI reads them, and lines that NCBI's line reader
# joins (b), which LOSAT Web read alike; BLASTX has no such row and runs against the native CLI alone.
SCAN_FIXTURES = [
    ("blastn", "gap_100.q.nuc.fa", "base.s.nuc.fa"), ("blastn", "base.q.nuc.fa", "gap_mods_bad.s.nuc.fa"),
    ("blastn", "eol_mixed_lf_crlf_cr.q.nuc.fa", "base.s.nuc.fa"), ("blastn", "base.q.nuc.fa", "eol_mixed_cr_lf.s.nuc.fa"),
    ("blastp", "gap_100.q.prot.faa", "base.s.prot.faa"), ("blastp", "base.q.prot.faa", "gap_unk100.s.prot.faa"),
    ("blastp", "eol_lfcr.q.prot.faa", "base.s.prot.faa"), ("blastp", "base.q.prot.faa", "eol_mixed_lf_crlf_cr.s.prot.faa"),
    ("tblastx", "gap_unk100.q.nuc.fa", "base.s.nuc.fa"), ("tblastx", "base.q.nuc.fa", "eol_lfcr.s.nuc.fa"),
    ("tblastn", "base.q.prot.faa", "gap_100.s.nuc.fa"), ("tblastn", "eol_mixed_lf_crlf_cr.q.prot.faa", "base.s.nuc.fa"),
]
# NCBI BLAST+ 2.17.0's outfmt 6 of the D16 search (-conffile nol.ini), run by line WQA's survey of D12 on 2026-11-20
# (/home/kawato/losat-baselines/web-qa-20261119/survey-d12/probe2.sh, `sid.pdb`, output o/sid.pdb.n).
D16_NCBI_OUTFMT6 = "64a3b01939bf72e8d0b9a3dbc3f4b6ad6630a3bb898bcfad241848c414554df5"


def scan_d13() -> list[dict]:
    """The inputs of line WQA's D13 and D16 (2026-11-20), each compared with the native CLI: the SCAN_FIXTURES with
    NCBI's frozen hashes, a BLASTX gap line on both sides; a first line that may be a Seq-id, which the CLI refuses
    with the data loaders on (c), in the query and in the subject (`"expect": "failure"`: the CLI's error, no output);
    a record whose CR line end NCBI's line reader joins to a `!` line in a file of LF line ends (b; LOSAT Web refused
    it: web/adapter/tests/fixtures/d13/joined.nuc.fa, made from base.q.nuc.fa), as the query and as the subject;
    and the D16 search, whose `-conffile nol.ini` (`[BLAST] DATA_LOADERS=none`, web/adapter/tests/fixtures/d16)
    turns the data loaders off, so that the reader reads the query's first line `2HHB` as data, as the CLI and NCBI
    do, and the same search without it, which the CLI refuses."""
    sys.path.insert(0, str(REPO / "LOSAT/tests"))
    import fasta_input_fixtures as fixtures  # noqa: E402  (the row definitions and frozen hashes)
    rows: dict[tuple, list] = {}
    for row in fixtures.read_manifest():
        # BLASTN's rows of its default task (megablast), which the argv leaves out.
        if row.kind == "ncbi" and not row.stdin and not row.env and row.threads == "1" and not row.args \
                and row.task in ("", "megablast"):
            key = (row.program, row.query.removeprefix("@/"), row.subject.removeprefix("@/"))
            rows.setdefault(key, []).append(row)
    searches = []
    for program, query, subject in SCAN_FIXTURES:
        found = rows[(program, query, subject)]
        assert {row.rc for row in found} == {"0"}, (program, query, subject)
        searches.append({"program": program, "cwd": str(fixtures.ENGINE), "group": "d13_scan",
                         "argv": [program, "-query", fixtures.rp(f"@/{query}"), "-subject", fixtures.rp(f"@/{subject}")],
                         "cases": [row.row_id for row in found],
                         "frozen": {row.outfmt: row.stdout_sha256 for row in found if row.outfmt in ("0", "6", "7")}})
    fasta = "tests/fixtures/fasta_input"
    native = {"cwd": str(REPO / "LOSAT"), "frozen": {}, "group": "d13_scan"}
    searches.append({**native, "program": "blastx", "cases": ["d13.blastx.gap"],
                     "argv": ["blastx", "-query", f"{fasta}/gap_100.q.nuc.fa", "-subject", f"{fasta}/gap_unk100.s.prot.faa"]})
    searches.append({**native, "program": "blastn", "cases": ["d13.blastn.seqid_query"], "expect": "failure",
                     "argv": ["blastn", "-query", f"{fasta}/seqid0.q.nuc.fa", "-subject", f"{fasta}/base.s.nuc.fa"]})
    searches.append({**native, "program": "blastp", "cases": ["d13.blastp.seqid_subject"], "expect": "failure",
                     "argv": ["blastp", "-query", f"{fasta}/base.q.prot.faa", "-subject", f"{fasta}/seqid1.s.prot.faa"]})
    joined = str(REPO / "web/adapter/tests/fixtures/d13/joined.nuc.fa")
    searches.append({**native, "program": "blastn", "cases": ["d13.blastn.joined_query"],
                     "argv": ["blastn", "-query", joined, "-subject", f"{fasta}/base.s.nuc.fa"]})
    searches.append({**native, "program": "blastn", "cases": ["d13.blastn.joined_subject"],
                     "argv": ["blastn", "-query", f"{fasta}/base.q.nuc.fa", "-subject", joined]})
    d16 = {"program": "blastp", "cwd": str(REPO / "web/adapter/tests/fixtures/d16"), "group": "d13_scan"}
    searches.append({**d16, "cases": ["d16.blastp.conffile_no_loaders"], "frozen": {"6": D16_NCBI_OUTFMT6},
                     "argv": ["blastp", "-conffile", "nol.ini", "-query", "sid.q.faa", "-subject", "sid.s.faa"]})
    searches.append({**d16, "cases": ["d16.blastp.loaders_on"], "frozen": {}, "expect": "failure",
                     "argv": ["blastp", "-query", "sid.q.faa", "-subject", "sid.s.faa"]})
    return searches


def failed_outputs() -> list[dict]:
    """Runs that fail after the CLI wrote part of its output (line WQA's sweep, D15, 2026-11-20): the
    run sends what the output of the argv's own `-outfmt` held, so that it is the CLI's file byte for
    byte, and drops the other outputs (abi_v2.md §5). BLASTP with a query of two batches whose second
    has no residues (`tests/fixtures/failed_outputs`): the CLI writes the first batch's hits, then
    NCBI's `Sequence contains no data` ends the run, in outfmt 0, 6, 7, and 5 and `6 qseqid sseqid`
    (the output of the argv's own, stream 4); and BLASTN with `-parse_deflines` on a record that is
    both a query and a subject and has masks, whose outfmt 0 LOSAT refuses after its prolog, where the
    CLI's outfmt 6 and 7 complete (their streams are dropped with the failed run's)."""
    batch = "web/adapter/tests/fixtures/failed_outputs"
    searches = []
    # The archives of outfmt 13 and 14: the CLI leaves `results` empty (no side file, the first batch's report
    # is not written before the error), and the run's archive holds that master file as an empty entry.
    for outfmt, name in ((None, "outfmt0"), ("6", "outfmt6"), ("7", "outfmt7"), ("5", "outfmt5"),
                         ("6 qseqid sseqid", "outfmt6_fields"), ("13", "outfmt13"), ("14", "outfmt14")):
        searches.append({"program": "blastp",
                         "argv": ["blastp", "-query", f"{batch}/later_batch.q.faa", "-subject", f"{batch}/later_batch.s.faa"],
                         "cwd": str(REPO), "cases": [f"failed_outputs.blastp.later_batch.{name}"], "frozen": {},
                         "expect": "failure", **({"outfmt": outfmt} if outfmt else {}),
                         **({"files": True} if outfmt in ("13", "14") else {})})
    searches.append({"program": "blastn",
                     "argv": ["blastn", "-query", f"{batch}/later_batch.q.fna", "-subject", f"{batch}/later_batch.s.fna"],
                     "cwd": str(REPO), "cases": ["failed_outputs.blastn.later_batch.outfmt13"], "frozen": {},
                     "expect": "failure", "outfmt": "13", "files": True})
    searches.append({"program": "blastn",
                     "argv": ["blastn", "-query", "tests/fasta/outfmt0/pd_self_n.fna", "-subject", "tests/fasta/outfmt0/pd_self_n.fna",
                              "-parse_deflines", "-dust", "yes"],
                     "cwd": str(REPO / "LOSAT"), "cases": ["failed_outputs.blastn.pd_self.outfmt0"], "frozen": {},
                     "expect": "failure"})
    return searches


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
        searches += own_search()
        searches += scan_d13()
        searches += failed_outputs()
    args.out.write_text(json.dumps(searches, indent=2) + "\n")
    print(f"{len(searches)} searches -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
