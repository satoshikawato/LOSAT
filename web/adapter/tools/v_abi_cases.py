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
reduces in the wasm-threads-only path, the searches of a query whose reader fails in its second
query batch (`reader_error`), the BLASTN `-window_masker_db` and `-window_masker_taxid`
searches (`window_masker`), whose counts files the host registers by name (abi_v2.md §4 `file_begin`),
the runs whose outcome and standard error are those of the argv's own `-outfmt` (`own_search`), the
`-db` searches of the db_search rows that the wasm-threading check runs and of two filters (`database`),
whose files the engine opens in the run's working directory (abi_v2.md §3; line WQA D12a, 2026-11-20), with
a query read with the data loaders on and off (`-conffile` of `DATA_LOADERS=none`), and the
inputs of line WQA's D13 and D16 (`scan_d13`: gap lines, joined lines, Seq-id first lines, and
`-conffile` with `[BLAST] DATA_LOADERS=none`, whose registry file the host registers as kind 2), and BLASTN
`-filtering_db` (line RG item D4): with `-db` (`REPEAT_ROWS`, `REPEAT_FAILURE_ROWS`) and with `-subject`
(`repeat_subject`, `"run_directory": true`: the subject is registered as for every search, and the run goes to
the instance whose `/work` holds the repeat database).
`--group fasta_input` writes that group of the full suite alone.

Both suites add the tool runs of `losat_web2_tool` and `losat_web2_validate` for `windowmasker` (`tool`; line NP
WW3a, 2026-11-22): a search with `"tool"` is one run of the command line `losat windowmasker <argv>` (`"tool": "run"`)
or one check of its argv (`"tool": "validate"`), in a directory that holds the files of `"files"` (published path -> source
file). v_abi.js runs it on an instance whose `/work` is a folder of those files, and the native CLI in another folder of
the same files, and compares exit status, standard output, standard error and every file the run wrote.

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


# The db_search rows of LOSAT/tests/check_wasm_threading.py (one outfmt 6 row per program) and two filters of a
# BLASTP database: -seqidlist (a binary list) and -taxids (line WQA D12a, 2026-11-20).
DATABASE_ROWS = ["db.blastn.mega.6", "db.d1.6", "db.blastx.default.6", "db.tblastn.default.6", "db.tblastx.default.6",
                 "db.sl.bp.bsl.6", "db.tf.bp.nosql.6"]
# BLASTN `-filtering_db` with `-db` (line RG item D4, D4-web): a second database in the run's directory, whose
# files the engine opens by `-filtering_db`'s name. Rows of db.rpt.* in LOSAT/tests/db_search_manifest/blastn.rpt.tsv
# without an environment: four tasks and both strands in outfmt 6, outfmt 0, the archive (11) with the lower case
# joined to the repeats, a -query_loc edge, and the XML of the filter (5).
REPEAT_ROWS = ["db.rpt.mb.6", "db.rpt.bn.6", "db.rpt.short.6", "db.rpt.dc.6", "db.rpt.plus.6", "db.rpt.minus.6", "db.rpt.mb.0",
               "db.rpt.edge_lo.0", "db.rpt.lc.11", "db.rpt.edge_lo.11", "db.rpt.mb.5"]
# The repeat search's errors through the adapter (D4-web; review N4): a repeat database whose header counts are
# wrong and a query without letters, whose NCBI messages name no path (exit 3; the native CLI's run is compared).
REPEAT_FAILURE_ROWS = ["db.rpt.badcount.6", "db.rpt.noletters.6"]
# The -subject run of a repeat database (D4-web): rpt.* rows of LOSAT/tests/fixtures/blastn_regression/manifest.tsv.
REPEAT_SUBJECT_ROWS = ["rpt.megablast", "rpt.blastn_fmt0", "rpt.lc", "rpt.fmt11"]


def database() -> list[dict]:
    """The `-db` searches (`DATABASE_ROWS`), run from LOSAT/ as the manifest's paths are, with NCBI's frozen
    outfmt 6 (LOSAT/tests/db_search_manifest). The host gives the engine LOSAT/ as the run's working
    directory, and the run no subject."""
    sys.path.insert(0, str(REPO / "LOSAT/tests"))
    import db_search_fixtures  # noqa: E402  (the manifest reader)
    rows = {row["row_id"]: row for row in db_search_fixtures.read_manifest()[1]}
    searches = []
    for row_id in DATABASE_ROWS:
        row = rows[row_id]
        assert row["kind"] == "ncbi" and row["outfmt"] == "6", row_id
        searches.append({"program": row["program"],
                         "argv": [row["program"], "-db", row["db"], "-query", row["query"], *shlex.split(row["extra_args"])],
                         "cwd": str(REPO / "LOSAT"), "cases": [row_id], "frozen": {"6": row["stdout_sha256"]},
                         "database": True})
    for row_id in REPEAT_ROWS:
        row = rows[row_id]
        assert row["kind"] == "ncbi" and row["rc"] == "0" and row["program"] == "blastn", row_id
        argv = ["blastn", "-db", row["db"], "-query", row["query"], *shlex.split(row["extra_args"])]
        search = {"program": "blastn", "argv": argv, "cwd": str(REPO / "LOSAT"), "cases": [row_id], "frozen": {},
                  "database": True}
        if custom_outfmt(row["outfmt"]):
            search["outfmt"] = row["outfmt"]
            search["frozen"]["custom"] = row["stdout_sha256"]
        else:
            search["frozen"][row["outfmt"]] = row["stdout_sha256"]
        searches.append(search)
    for row_id in REPEAT_FAILURE_ROWS:
        row = rows[row_id]
        assert row["kind"] == "ncbi" and row["rc"] == "3" and row["program"] == "blastn", row_id
        searches.append({"program": "blastn", "argv": ["blastn", "-db", row["db"], "-query", row["query"], *shlex.split(row["extra_args"])],
                         "cwd": str(REPO / "LOSAT"), "cases": [row_id], "frozen": {}, "expect": "failure", "database": True})
    searches += repeat_subject()
    # With -parse_deflines the run registers its query with `register_inputs` (its subject handle 0, the
    # response's subject null), and the records' `s_idx` are the database's OIDs (abi_v2.md §4, §8); the
    # native CLI's run of the argv is the reference (review of D12a, low 3).
    parsed = dict(searches[DATABASE_ROWS.index("db.d1.6")])
    searches.append({**parsed, "argv": [*parsed["argv"], "-parse_deflines"], "cases": ["db.d1.6.parse_deflines"], "frozen": {}})
    # A run of -db reads its query with the data loaders of its argv, as every run's inputs are read (line WQA's D13
    # and D16 with D12a, 2026-11-21): with the `-conffile` of `[BLAST] DATA_LOADERS=none` (web/adapter/tests/fixtures/d16)
    # the query's first line `2HHB` is data, and without it the CLI refuses that line. Run from the repository's root,
    # which holds the database and the D16 files; the native CLI's run of each argv is the reference.
    d16 = "web/adapter/tests/fixtures/d16"
    argv = ["blastp", "-db", f"LOSAT/{rows['db.d1.6']['db']}", "-query", f"{d16}/sid.q.faa"]
    searches.append({"program": "blastp", "argv": [argv[0], "-conffile", f"{d16}/nol.ini", *argv[1:]], "cwd": str(REPO),
                     "cases": ["db.d16.blastp.conffile_no_loaders"], "frozen": {}, "database": True})
    searches.append({"program": "blastp", "argv": argv, "cwd": str(REPO), "cases": ["db.d16.blastp.loaders_on"],
                     "frozen": {}, "expect": "failure", "database": True})
    # A run of -db reads its query as every run reads it (line WQA's D10 with D12a, 2026-11-22): 30 copies of the
    # db.d1 query (358 letters each) fill more than BLASTP's first batch of 10000 letters, and a reader's error on
    # the record after them ends the run in its batch, after the database search of the batches before, as on the
    # CLI (the native CLI's run of the argv is the reference). Written under LOSAT/target, as reader_error's inputs.
    directory = Path("target/v_abi_reader_error")
    (REPO / "LOSAT" / directory).mkdir(parents=True, exist_ok=True)
    d1 = rows["db.d1.6"]
    query = (REPO / "LOSAT" / d1["query"]).read_bytes() * 30 + b">bad\nMKVL%%%%\n"
    (REPO / "LOSAT" / directory / "blastp_db.q.fa").write_bytes(query)
    searches.append({"program": "blastp", "argv": ["blastp", "-db", d1["db"], "-query", str(directory / "blastp_db.q.fa")],
                     "cwd": str(REPO / "LOSAT"), "cases": ["db.d1.reader_error"], "frozen": {}, "expect": "failure",
                     "database": True})
    return searches


def repeat_subject() -> list[dict]:
    """`-subject` searches with `-filtering_db` (line RG item D4, D4-web): the subject is registered as for every
    subject run, and the repeat database is a set of files in the run's working directory (LOSAT/, as for `-db`):
    `"run_directory"` has the harness run it on the instance whose `/work` is that directory. The rows are
    NCBI-frozen (`rpt.*` of the blastn regression manifest)."""
    sys.path.insert(0, str(REPO / "LOSAT/tests"))
    import blastn_regression_fixtures as regression  # noqa: E402  (the manifest reader)
    rows = {row["case_id"]: row for row in regression.read_manifest()}
    searches = []
    for case_id in REPEAT_SUBJECT_ROWS:
        row = rows[case_id]
        assert row["exit"] == "0" and not row["env"], case_id
        outfmt, argv = split(["blastn", *shlex.split(row["argv"])])
        outfmt = outfmt or "0"
        search = {"program": "blastn", "argv": argv, "cwd": str(REPO / "LOSAT"), "cases": [case_id], "frozen": {},
                  "run_directory": True}
        if custom_outfmt(outfmt):
            search["outfmt"] = outfmt
            search["frozen"]["custom"] = row["stdout_sha256"]
        else:
            search["frozen"][outfmt] = row["stdout_sha256"]
        searches.append(search)
    return searches


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


# The copies of each program's subject sequence in `reader_error`: the outfmt 0 report of the first
# batch is then over 1 MiB (the CLI's file before its error: 1.26 to 1.46 MB, measured 2026-11-20).
COPIES = {"blastn": 72, "blastp": 30, "tblastn": 16, "blastx": 72, "tblastx": 10}


def reader_error() -> list[dict]:
    """A query whose reader fails on a record of its second query batch (line WQA's sweep, D10,
    2026-11-19): the CLI searches the first batch, writes its report and ends with the reader's error
    (`CFastaReader: Near line N, ...`, exit 1) in the batch that reads the record; `register` takes
    the records before it with the error (abi_v2.md §4), and the run ends with the same error. The
    records fill more than each program's first batch (BLASTP 10000 letters, TBLASTN 20000, BLASTX
    and TBLASTX 10002, BLASTN's `CBatchSizeMixer` batch for a subject of 5000 letters), as in
    LOSAT/src/api/reader_end_tests.rs; the inputs are written under LOSAT/target, as
    window_masker's. A BLASTP search with -parse_deflines: two records of one ID in the first batch
    end the CLI there, before the record of the second batch whose Seq-id names a nucleotide
    sequence (`register_inputs` reads it, and the run ends with the first batch's error, after the
    CLI writes the outfmt 0 prolog, as NCBI writes it when it stops in its first batch: line WQA D10
    round 4, 2026-11-20; before, LOSAT refused before the prolog). The subject holds copies of its
    sequence, so that the outfmt 0 report of the first batch is over the adapter's chunk of 1 MiB
    (`emit::CHUNK`): a run sends whole chunks as they fill and nothing smaller before it fails, so
    the run's first chunk of outfmt 0, which the CLI's file starts with, shows the batch before the
    error (`sends`; review L3, 2026-11-20)."""
    import random

    engine = REPO / "LOSAT"
    directory = Path("target/v_abi_reader_error")
    (engine / directory).mkdir(parents=True, exist_ok=True)
    rng = random.Random(20261119)
    codons = dict(zip("ACDEFGHIKLMNPQRSTVWY", ["GCT", "TGT", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT", "AAA", "CTG",
                                                "ATG", "AAT", "CCG", "CAG", "CGT", "AGC", "ACC", "GTG", "TGG", "TAT"]))

    def letters(alphabet: str, length: int) -> str:
        return "".join(rng.choice(alphabet) for _ in range(length))

    def protein(length: int) -> str:
        return letters("ACDEFGHIKLMNPQRSTVWY", length)

    def coding(residues: str) -> str:
        return "".join(codons[residue] for residue in residues)

    def fasta(names: list[str], sequence: str) -> bytes:
        rows = [sequence[i:i + 60] for i in range(0, len(sequence), 60)]
        lines = [line for name in names for line in (f">{name}", *rows)]
        return ("\n".join(lines) + "\n").encode()

    bad = b">bad\nACGTAC%%%%\n"
    nucleotide, residues, long_residues = letters("ACGT", 5000), protein(1100), protein(2100)
    # The query records, the query and subject sequences, and the copies of the subject sequence.
    inputs = {
        "blastn": (40, nucleotide, nucleotide, COPIES["blastn"]),
        "blastp": (12, residues, residues, COPIES["blastp"]),
        "tblastn": (12, long_residues, coding(long_residues), COPIES["tblastn"]),
        "blastx": (5, coding(residues), residues, COPIES["blastx"]),
        "tblastx": (5, coding(residues), coding(residues), COPIES["tblastx"]),
    }
    searches = []
    for program, (count, query, subject, copies) in inputs.items():
        names = [f"q{i}" for i in range(1, count + 1)]
        (engine / directory / f"{program}.q.fa").write_bytes(fasta(names, query) + bad)
        subjects = b"".join(fasta([f"s{i}"], subject) for i in range(1, copies + 1))
        (engine / directory / f"{program}.s.fa").write_bytes(subjects)
        searches.append({"program": program, "argv": [program, "-query", str(directory / f"{program}.q.fa"),
                                                      "-subject", str(directory / f"{program}.s.fa")],
                         "cwd": str(engine), "cases": [f"reader_error.{program}"], "frozen": {}, "expect": "failure",
                         "sends": [0]})
    names = ["q1", "q1", *(f"q{i}" for i in range(3, 12)), "ref|NM_000001.1|"]
    (engine / directory / "blastp_pd.q.fa").write_bytes(fasta(names, residues))
    (engine / directory / "blastp_pd.s.fa").write_bytes(fasta(["s1"], residues[:300]))
    searches.append({"program": "blastp", "argv": ["blastp", "-query", str(directory / "blastp_pd.q.fa"),
                                                   "-subject", str(directory / "blastp_pd.s.fa"), "-parse_deflines"],
                     "cwd": str(engine), "cases": ["reader_error.blastp_pd_order"], "frozen": {}, "expect": "failure"})
    return searches


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


# The files that a tool run sees in its directory, published path -> source file under LOSAT/tests (`tool`):
# the windowmasker fixtures of line NP M5 and the counts of winmask/, small ones, and a set of the same names for every row,
# as every file of a folder is there for the command line.
TOOL_FILES = {
    "fa/dm2L_3k.fa": "winmask/fasta/dm2L_3k.fa", "fa/dm2L_60k.fa": "winmask/fasta/dm2L_60k.fa",
    "fa/dm_edge.fa": "winmask/fasta/dm_edge.fa", "fa/multi.fa": "windowmasker/in/multi.fa",
    "fa/cd.fa": "windowmasker/in/cd.fa", "fa/gap.fa": "windowmasker/in/gap.fa", "fa/short.fa": "windowmasker/in/short.fa",
    "fa/letters.fa": "windowmasker/in/letters.fa",
    "counts/dm2L_u9.binary": "winmask/counts/dm2L_u9.binary", "counts/mg1655_u9.obinary": "winmask/counts/mg1655_u9.obinary",
    "counts/mg1655_u9.ascii": "winmask/counts/mg1655_u9.ascii", "counts/mg1655_u9.oascii": "winmask/counts/mg1655_u9.oascii",
    "counts/dm_u2.obinary": "winmask/counts/dm_u2.obinary",
    "counts/u4.ascii": "windowmasker/counts/u4.ascii", "counts/u4.binary": "windowmasker/counts/u4.binary",
    "counts/u4.oascii": "windowmasker/counts/u4.oascii", "counts/u4.obinary": "windowmasker/counts/u4.obinary",
    "counts/garbage.counts": "windowmasker/in/garbage.counts", "counts/empty.counts": "windowmasker/in/empty.counts",
    "ids/a.txt": "windowmasker/in/ids_a.txt", "ids/b.txt": "windowmasker/in/ids_b.txt", "ids/db.txt": "windowmasker/in/ids_db.txt",
    **{f"db/{name}.{extension}": f"windowmasker/db/{name}.{extension}"
       for name, extensions in (("wm", "ndb nhr nin njs not nsq ntf nto"), ("wm_ps", "ndb nhr nin njs nog nos not nsq ntf nto"))
       for extension in extensions.split()},
}
# The bytes of a tool run's FASTA input of more than 1 MiB (the file it writes arrives in two frames at least).
TOOL_BIG_LETTERS = 1_500_000


def tool() -> list[dict]:
    """The tool runs of `windowmasker` (`"tool": "run"`) and the checks of its argv (`"tool": "validate"`), compared
    with the command line's run of the same words in a folder of the same files (line NP WW3a, 2026-11-22,
    D-WW3-6, D-WW3-7, D-WW3-8): the modes `-mk_counts` (each `-sformat`, to the standard output and to a file),
    `-convert` and `-ustat` (each `-outfmt`), their options, the help and version words, and the errors that the
    command line has before and after it opens a file. `"expect"` is `"failure"` for a run whose exit status is
    not 0 (a check of v_abi.js against the native run), `"late_failure"` for an argv that `validate` accepts though
    the command line's run fails after it opens a file (the checks stop there, D-WW3-10), and `"refusal"` for LOSAT's
    refusal of an output name that LOSAT Web does not write (the command line would write it; exit status 1, no
    file); `"frames"` is the least number of frames of bytes that the largest file of the run arrives in. Rows are written with the words of the command line (the mode's own
    options, the files by their published names)."""
    tests = REPO / "LOSAT/tests"
    directory = REPO / "LOSAT/target/v_abi_tool"
    directory.mkdir(parents=True, exist_ok=True)
    files = {name: str(tests / source) for name, source in TOOL_FILES.items()}
    # A list of FASTA files (`-fa_list T`), written for the published names; and a FASTA file of 1.5 million letters
    # (a fixed generator: the same bytes in every run), which a masked output of more than 1 MiB follows.
    (directory / "falist.txt").write_text("fa/dm2L_3k.fa\n\nfa/multi.fa\n")
    files["falist.txt"] = str(directory / "falist.txt")
    (directory / "falist_missing.txt").write_text("m_list.fa\n")
    files["falist_missing.txt"] = str(directory / "falist_missing.txt")
    # A registry file of the run's directory (`-conffile wm.ini`), which the run and its validation read.
    (directory / "wm.ini").write_text("[WINDOW_MASKER]\n")
    files["wm.ini"] = str(directory / "wm.ini")
    (directory / "keep").write_text("k\n")
    with_logs = {**files, "logs/keep": str(directory / "keep")}
    state, lines = 12345, [">big generated record\n"]
    letters = []
    for _ in range(TOOL_BIG_LETTERS):
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        letters.append("ACGT"[(state >> 16) & 3])
    text = "".join(letters)
    lines += [text[at:at + 60] + "\n" for at in range(0, len(text), 60)]
    (directory / "big.fa").write_text("".join(lines))
    big = {**files, "big.fa": str(directory / "big.fa")}
    rows: list[dict] = []

    def add(case_id: str, words: str, failure: bool | str = False, kind: str = "run", frames: int | None = None,
            published: dict | None = None) -> None:
        row = {"program": "windowmasker", "tool": kind, "argv": shlex.split(words), "cwd": str(tests),
               "files": published or files, "cases": [f"tool.{case_id}"], "frozen": {}}
        if failure:
            row["expect"] = failure if isinstance(failure, str) else "failure"
        if frames is not None:
            row["frames"] = frames
        rows.append(row)

    mask = "-ustat counts/dm2L_u9.binary"
    # -mk_counts: each -sformat to the standard output (the default, `-out -` and without it) and to a file.
    add("mk.default.stdout", "-mk_counts -in fa/dm2L_60k.fa -unit 5")
    for sformat in ("ascii", "binary", "oascii", "obinary"):
        base = f"-mk_counts -in fa/dm2L_60k.fa -unit 5 -sformat {sformat}"
        add(f"mk.{sformat}.stdout", base)
        add(f"mk.{sformat}.out_dash", f"{base} -out -")
        add(f"mk.{sformat}.file", f"{base} -out out.counts")
    add("mk.t_low_pct.file", "-mk_counts -in fa/dm2L_60k.fa -unit 5 -t_low_pct 80 -t_extend_pct 95 -t_thres_pct 97 "
        "-t_high_pct 99.9 -sformat binary -out out.counts")
    add("mk.t_low_pct.stdout", "-mk_counts -in fa/dm2L_60k.fa -unit 5 -t_low_pct 80 -sformat obinary")
    add("mk.checkdup", "-mk_counts -in fa/dm2L_60k.fa,fa/cd.fa -checkdup T -unit 3 -sformat binary")
    add("mk.checkdup.file", "-mk_counts -in fa/dm2L_60k.fa,fa/cd.fa -checkdup T -unit 3 -sformat obinary -out out.counts")
    add("mk.fa_list", "-mk_counts -in falist.txt -fa_list T -unit 5 -sformat binary")
    add("mk.fa_list.missing", "-mk_counts -in nonexistent.txt -fa_list T", failure=True)
    add("mk.comma", "-mk_counts -in fa/dm2L_3k.fa,fa/multi.fa -unit 5 -sformat binary")
    add("mk.blastdb", "-mk_counts -infmt blastdb -in db/wm -sformat binary")
    add("mk.blastdb.file", "-mk_counts -infmt blastdb -in db/wm,db/wm_ps -unit 5 -sformat obinary -out out.counts")
    add("mk.blastdb.missing", "-mk_counts -infmt blastdb -in db/none -unit 5", failure=True)
    add("mk.ids", "-mk_counts -in fa/multi.fa -ids ids/a.txt -unit 5 -sformat binary")
    add("mk.exclude_ids", "-mk_counts -in fa/multi.fa -exclude_ids ids/a.txt -unit 5 -sformat binary")
    add("mk.parse_seqids", "-mk_counts -in fa/dm2L_3k.fa -parse_seqids -sformat binary")
    add("mk.meta", "-mk_counts -in fa/dm2L_3k.fa -unit 5 -meta something -sformat ascii")
    add("mk.genome_size", "-mk_counts -in fa/dm2L_3k.fa -genome_size 1000 -unit 3 -sformat binary")
    add("mk.logfile", "-mk_counts -in fa/dm2L_3k.fa -unit 5 -sformat binary -logfile run.log")
    add("mk.logfile.file", "-mk_counts -in fa/dm2L_3k.fa -unit 5 -out out.counts -logfile run.log")
    # -convert: each target format, from each source format, to the standard output and to a file.
    for source in ("ascii", "binary", "oascii", "obinary"):
        for target in ("ascii", "binary", "oascii", "obinary"):
            if source == "ascii" or target == "ascii" or source[0] == target[0]:
                add(f"convert.{source}.{target}", f"-convert -in counts/u4.{source} -sformat {target}")
    for target in ("ascii", "binary", "oascii", "obinary"):
        add(f"convert.ascii.{target}.file", f"-convert -in counts/u4.ascii -sformat {target} -out out.counts")
    add("convert.default", "-convert -in counts/u4.obinary")
    add("convert.meta", "-convert -in counts/u4.obinary -sformat ascii -meta 'replaced note'")
    add("convert.logfile", "-convert -in counts/u4.obinary -sformat oascii -out out.counts -logfile run.log")
    add("convert.mg1655", "-convert -in counts/mg1655_u9.obinary -sformat ascii -out out.counts")
    # -ustat: every -outfmt on the standard output and in a file, and with -parse_seqids.
    outfmts = ("interval", "fasta", "seqloc_asn1_bin", "seqloc_asn1_text", "seqloc_xml", "maskinfo_asn1_bin",
               "maskinfo_asn1_text", "maskinfo_xml")
    for outfmt in outfmts:
        add(f"ustat.{outfmt}.stdout", f"{mask} -in fa/multi.fa -outfmt {outfmt}")
        add(f"ustat.{outfmt}.file", f"{mask} -in fa/multi.fa -outfmt {outfmt} -out out.res")
        add(f"ustat.{outfmt}.parse_seqids", f"{mask} -in fa/multi.fa -outfmt {outfmt} -parse_seqids")
    add("ustat.default", f"{mask} -in fa/dm2L_60k.fa")
    add("ustat.default.file", f"{mask} -in fa/dm2L_60k.fa -out out.res")
    add("ustat.dust", f"{mask} -in fa/dm2L_60k.fa -dust T")
    add("ustat.dust.fasta.file", f"{mask} -in fa/dm2L_60k.fa -dust T -dust_level 10 -outfmt fasta -out out.fa")
    add("ustat.ids", f"{mask} -in fa/multi.fa -ids ids/a.txt")
    add("ustat.exclude_ids", f"{mask} -in fa/multi.fa -exclude_ids ids/b.txt -outfmt fasta")
    add("ustat.ids_and_exclude_ids", f"{mask} -in fa/dm2L_3k.fa -ids ids/a.txt -exclude_ids ids/b.txt", failure=True)
    add("ustat.text_match", f"{mask} -in fa/multi.fa -ids ids/a.txt -text_match F", failure=True)
    add("ustat.seqids", f"{mask} -infmt seqids -in ids/a.txt", failure=True)
    add("ustat.parse_seqids", f"{mask} -in fa/multi.fa -parse_seqids -exclude_ids ids/a.txt")
    add("ustat.blastdb", f"{mask} -infmt blastdb -in db/wm -outfmt fasta")
    add("ustat.blastdb.file", f"{mask} -infmt blastdb -in db/wm_ps -parse_seqids -outfmt maskinfo_asn1_text -out out.res")
    add("ustat.blastdb.missing", f"{mask} -infmt blastdb -in db/none", failure=True)
    add("ustat.comma", f"{mask} -in fa/dm2L_3k.fa,fa/multi.fa", failure=True)
    add("ustat.thresholds", f"{mask} -in fa/dm2L_60k.fa -t_extend 60 -t_thres 90 -t_high 300 -t_low 50 -window 13")
    add("ustat.t_low_pct.file", f"{mask} -in fa/dm2L_60k.fa -t_low_pct 90 -t_extend_pct 99 -out out.res")
    add("ustat.window5", f"{mask} -in fa/dm2L_3k.fa -window 5 -outfmt maskinfo_xml", failure=True)
    add("ustat.meta", f"{mask} -in fa/dm2L_3k.fa -meta something")
    add("ustat.obinary", "-ustat counts/mg1655_u9.obinary -in fa/dm_edge.fa -outfmt fasta")
    add("ustat.ascii", "-ustat counts/mg1655_u9.ascii -in fa/dm_edge.fa -out out.res")
    add("ustat.oascii", "-ustat counts/mg1655_u9.oascii -in fa/dm_edge.fa -outfmt maskinfo_xml -out out.res")
    add("ustat.u2", "-ustat counts/dm_u2.obinary -in fa/dm2L_3k.fa")
    add("ustat.gap", f"{mask} -in fa/gap.fa -dust T -outfmt fasta")
    add("ustat.logfile", f"{mask} -in fa/dm2L_3k.fa -logfile run.log")
    add("ustat.logfile.file", f"{mask} -in fa/dm2L_3k.fa -outfmt fasta -out out.fa -logfile run.log")
    add("ustat.stdin", f"{mask} -in -")
    # The output file of more than 1 MiB arrives in more than one frame.
    add("ustat.big.file", "-ustat counts/mg1655_u9.obinary -in big.fa -outfmt fasta -out out.fa", frames=2, published=big)
    add("ustat.big.stdout", "-ustat counts/mg1655_u9.obinary -in big.fa -outfmt fasta", published=big)
    add("mk.big.file", "-mk_counts -in big.fa -unit 9 -sformat obinary -out out.counts", published=big)
    # Errors, before and after the run opens a file.
    add("error.no_mode", "-in fa/dm2L_3k.fa", failure=True)
    add("error.no_arguments", "", failure=True)
    add("error.unknown_option", f"{mask} -nosuchopt", failure=True)
    add("error.extra_positional", f"{mask} extra", failure=True)
    add("error.duplicate_option", f"{mask} -window 20 -window 30", failure=True)
    add("error.window_abc", f"{mask} -in fa/dm2L_3k.fa -window abc", failure=True)
    add("error.sformat_with_ustat", f"{mask} -sformat ascii", failure=True)
    add("error.mk_and_convert", "-mk_counts -convert", failure=True)
    add("error.mk_and_ustat", f"-mk_counts {mask}", failure=True)
    add("error.outfmt_bad", f"{mask} -in fa/dm2L_3k.fa -outfmt nosuch", failure=True)
    add("error.in_missing", f"{mask} -in nonexistent.fa", failure=True)
    add("error.in_missing.mk", "-mk_counts -in nonexistent.fa -unit 5", failure=True)
    add("error.counts_missing", "-ustat nonexistent.counts -in fa/dm2L_3k.fa", failure=True)
    add("error.counts_garbage", "-ustat counts/garbage.counts -in fa/dm2L_3k.fa", failure=True)
    add("error.counts_empty", "-ustat counts/empty.counts -in fa/dm2L_3k.fa", failure=True)
    add("error.counts_garbage.logfile", "-ustat counts/garbage.counts -in fa/dm2L_3k.fa -logfile run.log", failure=True)
    add("error.convert_garbage", "-convert -in counts/garbage.counts -sformat ascii", failure=True)
    add("error.out_missing_directory", f"{mask} -in fa/dm2L_3k.fa -out nodir/x.txt", failure=True)
    add("mk.out_missing_directory", "-mk_counts -in fa/dm2L_60k.fa -unit 5 -out nodir/c.counts")
    add("error.out_directory", f"{mask} -in fa/dm2L_3k.fa -out fa", failure=True)
    add("error.unit_0", "-mk_counts -in fa/dm2L_3k.fa -unit 0", failure=True)
    add("ustat.ids_missing", f"{mask} -in fa/multi.fa -ids nonexistent.txt")
    add("error.mk.ids_and_exclude_ids", "-mk_counts -in fa/multi.fa -ids ids/a.txt -exclude_ids ids/b.txt -unit 5",
        failure=True)
    add("error.mk.seqids", "-mk_counts -infmt seqids -in ids/a.txt -unit 5", failure=True)
    add("error.convert.stdin", "-convert -in - -sformat ascii", failure=True)
    add("error.convert.stdin.out_dash", "-convert -in - -out -", failure=True)
    add("error.conffile_missing", "-conffile none.ini -mk_counts -in fa/dm2L_3k.fa -unit 5", failure=True)
    # A folder named as an input opens, and its read fails (EISDIR), as on Linux.
    add("error.in_directory", f"{mask} -in fa", failure=True)
    add("error.in_directory.mk", "-mk_counts -in fa -unit 5", failure=True)
    add("error.counts_directory", "-ustat fa -in fa/dm2L_3k.fa", failure=True)
    # One file that the run opens twice (`-out X -logfile X`, two spellings of one path): one file, as the system
    # shares it (D-WW3-9).
    add("mk.shared.logfile", "-mk_counts -in fa/dm2L_3k.fa -unit 5 -out out.counts -logfile out.counts")
    add("mk.shared.dot", "-mk_counts -in fa/dm2L_3k.fa -unit 5 -out ./out.counts -logfile out.counts")
    add("mk.shared.ascii", "-mk_counts -in fa/dm2L_60k.fa -unit 5 -sformat ascii -out out.counts -logfile out.counts")
    add("convert.shared.logfile", "-convert -in counts/u4.obinary -sformat oascii -out out.counts -logfile out.counts")
    add("ustat.shared.logfile", f"{mask} -in fa/dm2L_3k.fa -outfmt fasta -out out.fa -logfile out.fa")
    add("error.counts_garbage.shared", "-ustat counts/garbage.counts -in fa/dm2L_3k.fa -out out.res -logfile out.res",
        failure=True)
    # The registry file of -conffile, and the output names of the system's devices.
    add("mk.conffile", "-conffile wm.ini -mk_counts -in fa/dm2L_3k.fa -unit 5 -sformat binary")
    add("ustat.dev_stdout", f"{mask} -in fa/dm2L_3k.fa -out /dev/stdout")
    add("mk.dev_stdout", "-mk_counts -in fa/dm2L_3k.fa -unit 5 -sformat binary -out /dev/stdout")
    add("ustat.dev_null", f"{mask} -in fa/dm2L_3k.fa -out /dev/null")
    add("refusal.out_existing", f"{mask} -in fa/dm2L_3k.fa -out fa/short.fa", failure="refusal")
    add("refusal.out_dev_stderr", f"{mask} -in fa/dm2L_3k.fa -out /dev/stderr", failure="refusal")
    # -logfile against the run's real files (D19 cases t0-t4).
    add("d19.t0_nofolder", "-ustat counts/garbage.counts -in fa/dm2L_3k.fa -logfile logs/run.log", failure=True)
    add("d19.t0_nofolder_clean", "-mk_counts -in fa/dm2L_3k.fa -unit 5 -out x.counts -logfile logs/run.log")
    add("d19.t0_folder", "-ustat counts/garbage.counts -in fa/dm2L_3k.fa -logfile logs/run.log",
        failure=True, published=with_logs)
    add("d19.t1_folder", "-ustat counts/garbage.counts -in fa/dm2L_3k.fa -logfile fa", failure=True)
    add("d19.t1_folder_slash", "-ustat counts/garbage.counts -in fa/dm2L_3k.fa -logfile fa/", failure=True)
    add("d19.t2_input", "-ustat counts/garbage.counts -in fa/dm2L_3k.fa -logfile fa/dm2L_3k.fa", failure="refusal")
    add("d19.t2_missing_input", "-mk_counts -in m.fa -logfile m.fa", failure="refusal")
    add("d19.t2_falist_missing", "-mk_counts -fa_list true -in falist_missing.txt -logfile m_list.fa", failure="refusal")
    add("d19.t2_blastdb_missing", "-mk_counts -infmt blastdb -in db/none -logfile db/none.nal", failure="refusal")
    add("d19.t3_stdout", "-ustat counts/garbage.counts -in fa/dm2L_3k.fa -logfile /dev/stdout", failure=True)
    add("d19.t3_stdout_counts", "-mk_counts -in fa/dm2L_3k.fa -unit 5 -logfile /dev/stdout")
    add("d19.t3_stdout_both", "-mk_counts -in fa/dm2L_3k.fa -unit 5 -out /dev/stdout -logfile /dev/stdout")
    add("d19.t3_stdout_binary", "-mk_counts -in fa/dm2L_3k.fa -unit 5 -sformat binary -out - -logfile /dev/stdout")
    add("d19.t4_space", "-mk_counts -in fa/dm2L_3k.fa -unit 5 -out x.counts -logfile ' run.log'")
    add("d19.t4_slash", "-mk_counts -in fa/dm2L_3k.fa -unit 5 -out x.counts -logfile logs/")
    add("d19.t4_out_dot", "-mk_counts -in fa/dm2L_3k.fa -unit 5 -out ./x.counts -logfile ./run.log")
    add("d19.t4_out_up", "-mk_counts -in fa/dm2L_3k.fa -unit 5 -out fa/../x.counts -logfile fa/../run.log")
    # The help and the version words.
    for word in ("-h", "-help", "-xmlhelp", "-version", "-version-full", "-version-full-xml", "-version-full-json"):
        add(f"print{word}", word)
    add("print.help_with_options", f"{mask} -h")
    # losat_web2_validate of a tool's argv: valid -> 0 and nothing emitted; invalid -> -1 and the command line's
    # standard error (then the file of -logfile) without the last line end.
    for case_id, words, failure in (
            ("ustat", f"{mask} -in fa/dm2L_3k.fa", False),
            ("ustat.options", f"{mask} -in fa/multi.fa -outfmt maskinfo_xml -dust T -ids ids/a.txt -out out.res", False),
            ("mk_counts", "-mk_counts -in fa/dm2L_3k.fa -unit 5 -sformat obinary -out out.counts", False),
            ("convert", "-convert -in counts/u4.ascii -sformat obinary", False),
            ("logfile", f"{mask} -in fa/dm2L_3k.fa -logfile run.log", False),
            ("help", "-h", False), ("xmlhelp", "-xmlhelp", False), ("version", "-version", False),
            ("no_mode", "-in fa/dm2L_3k.fa", True), ("no_arguments", "", True),
            ("unknown_option", f"{mask} -nosuchopt", True), ("extra_positional", f"{mask} extra", True),
            ("duplicate_option", f"{mask} -window 20 -window 30", True), ("window_abc", f"{mask} -window abc", True),
            ("unit_17", "-mk_counts -unit 17", True), ("mk_and_ustat", f"-mk_counts {mask}", True),
            ("no_mode.logfile", "-in fa/dm2L_3k.fa -logfile run.log", True),
            ("unknown_option.logfile", f"{mask} -nosuchopt -logfile run.log", True),
            # The registry file of the run's directory, which validation reads (D-WW3-10).
            ("conffile", "-conffile wm.ini -mk_counts -in fa/dm2L_3k.fa -out out.counts", False),
            ("conffile_missing", "-conffile none.ini -mk_counts -in fa/dm2L_3k.fa -out out.counts", True),
            # The command line's refusals before it opens a file, and one after (D-WW3-10).
            ("mk.ids_and_exclude_ids", "-mk_counts -in fa/multi.fa -ids ids/a.txt -exclude_ids ids/b.txt", True),
            ("mk.seqids", "-mk_counts -infmt seqids -in ids/a.txt", True),
            ("convert.stdin", "-convert -in - -sformat ascii", True),
            ("convert.stdin.out_dash", "-convert -in - -out -", True),
            ("ustat.ids_and_exclude_ids", f"{mask} -in fa/dm2L_3k.fa -ids ids/a.txt -exclude_ids ids/b.txt",
             "late_failure"),
            ("shared", "-mk_counts -in fa/dm2L_3k.fa -unit 5 -out ./out.counts -logfile out.counts", False),
            ("d19.t2_input", "-ustat counts/garbage.counts -in fa/dm2L_3k.fa -logfile fa/dm2L_3k.fa", "refusal"),
            ("d19.t2_missing_input", "-mk_counts -in m.fa -logfile m.fa", "refusal"),
            ("d19.t2_falist_missing", "-mk_counts -fa_list true -in falist_missing.txt -logfile m_list.fa", "refusal"),
            ("d19.t2_blastdb_missing", "-mk_counts -infmt blastdb -in db/none -logfile db/none.nal", "refusal"),
            ("d19.t3_stdout", "-in fa/dm2L_3k.fa -logfile /dev/stdout", True)):
        add(f"validate.{case_id}", words, failure, kind="validate")

    # blast_formatter (`losat_web2_tool` and `losat_web2_validate`; line RG-F3p7):
    # -subject archives, local -db archives (with in/db/s3.* in /work), multi-object archives,
    # -out / -logfile / shared files, help/version, errors and LOSAT Web refusals.
    bf_tests = REPO / "LOSAT/tests/blast_formatter"
    bf_files: dict[str, str] = {}
    for name in ("bn.asn", "db.asn", "batch.asn", "empty.asn", "cut2.asn", "garbage.asn", "warn.asn", "nodb.asn"):
        bf_files[f"in/{name}"] = str(bf_tests / "in" / name)
    for db_file in sorted((bf_tests / "in/db").iterdir()):
        if db_file.is_file():
            bf_files[f"in/db/{db_file.name}"] = str(db_file)

    def add_bf(case_id: str, words: str, failure: bool | str = False, kind: str = "run") -> None:
        row = {"program": "blast_formatter", "tool": kind, "argv": shlex.split(words), "cwd": str(bf_tests),
               "files": bf_files, "cases": [f"tool.bf.{case_id}"], "frozen": {}}
        if failure:
            row["expect"] = failure if isinstance(failure, str) else "failure"
        rows.append(row)

    for outfmt in ("8", "9", "12"):
        add_bf(f"subject.of{outfmt}.stdout", f"-archive in/bn.asn -outfmt {outfmt}")
        add_bf(f"subject.of{outfmt}.file", f"-archive in/bn.asn -outfmt {outfmt} -out out.dat")
        add_bf(f"db.of{outfmt}.stdout", f"-archive in/db.asn -outfmt {outfmt}")
        add_bf(f"db.of{outfmt}.file", f"-archive in/db.asn -outfmt {outfmt} -out out.dat")
    add_bf("batch.of8.stdout", "-archive in/batch.asn -outfmt 8")
    add_bf("db.mts2", "-archive in/db.asn -outfmt 8 -max_target_seqs 2")
    add_bf("db.na0", "-archive in/db.asn -outfmt 8 -num_alignments 0")
    add_bf("warn.mts2", "-archive in/warn.asn -outfmt 8 -max_target_seqs 2")
    add_bf("logfile.file", "-archive in/bn.asn -outfmt 8 -out out.asn -logfile run.log")
    add_bf("shared.logfile", "-archive in/bn.asn -outfmt 8 -out shared.asn -logfile shared.asn")
    for word in ("-h", "-help", "-xmlhelp", "-version", "-version-full"):
        add_bf(f"print{word}", word)
    # Errors and refusals.
    add_bf("error.no_args", "", failure=True)
    # An empty last word on the wire (`blast_formatter\0-rid\0\0`, abi_v2.md §7): the CLI's "Empty RID string specified".
    add_bf("error.rid_empty", "-rid ''", failure=True)
    add_bf("error.rid_refused", "-rid ABCDEFGH01R", failure=True)
    add_bf("error.rid_archive", "-rid ABCDEFGH01R -archive in/bn.asn", failure=True)
    add_bf("error.of5_refused", "-archive in/bn.asn -outfmt 5", failure=True)
    add_bf("error.missing_archive", "-archive in/no_such.asn -outfmt 8", failure=True)
    add_bf("error.empty_archive", "-archive in/empty.asn -outfmt 8", failure=True)
    add_bf("error.cut_archive", "-archive in/cut2.asn -outfmt 8", failure=True)
    add_bf("trailing_garbage", "-archive in/garbage.asn -outfmt 8")
    add_bf("error.missing_db", "-archive in/nodb.asn -outfmt 8", failure=True)
    add_bf("refusal.out_existing", "-archive in/bn.asn -outfmt 8 -out in/bn.asn", failure="refusal")
    add_bf("refusal.logfile_existing", "-archive in/bn.asn -outfmt 8 -logfile in/bn.asn", failure="refusal")
    add_bf("refusal.logfile_missing_archive", "-archive missing.asn -outfmt 8 -logfile missing.asn", failure="refusal")
    # validate
    for case_id, words, failure in (
            ("subject", "-archive in/bn.asn -outfmt 8", False),
            ("db", "-archive in/db.asn -outfmt 12 -out report.json", False),
            ("help", "-h", False),
            ("xmlhelp", "-xmlhelp", False),
            ("version", "-version", False),
            ("no_args", "", True),
            ("rid_empty", "-rid ''", True),
            ("rid_refused", "-rid ABCDEFGH01R", True),
            ("rid_archive", "-rid X -archive in/bn.asn", True),
            ("unknown_opt", "-archive in/bn.asn -nosuchopt", True),
            ("late_of5", "-archive in/bn.asn -outfmt 5", "late_failure"),
            ("refusal.out_existing", "-archive in/bn.asn -outfmt 8 -out in/bn.asn", "refusal"),
            ("refusal.logfile_existing", "-archive in/bn.asn -outfmt 8 -logfile in/bn.asn", "refusal"),
            ("refusal.logfile_missing_archive", "-archive missing.asn -outfmt 8 -logfile missing.asn", "refusal")):
        add_bf(f"validate.{case_id}", words, failure, kind="validate")
    return rows


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
        searches += reader_error()
        searches += own_search()
        searches += database()
        searches += scan_d13()
        searches += failed_outputs()
        searches += tool()
    args.out.write_text(json.dumps(searches, indent=2) + "\n")
    print(f"{len(searches)} searches -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
