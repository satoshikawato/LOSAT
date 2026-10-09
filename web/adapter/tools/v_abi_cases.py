#!/usr/bin/env python3
"""The V-ABI cases (web/adapter/tests/v_abi.js) as JSON.

Every case is one search: a program, its argv without -outfmt, -out and -num_threads
(the adapter writes every format, the host sets the threads), the working directory,
and for each supported output format the frozen SHA-256 where one exists.

- full: the regression cases of docs/evidence/losat_web_e1a/capture_outputs.py (Gate A,
  the Stage G matrix, the BLASTP manifest in outfmt 0/6/7), grouped into searches.
  Their inputs use the recorded path spellings, so this needs the Gate A lexical root
  that capture_outputs.py documents. It adds the outfmt 0 and 7 fixtures of
  LOSAT/tests/outfmt0_manifest.tsv, whose frozen hashes are NCBI's.
  The group `fasta_input` (session SFc, S10) adds SF fixture inputs of
  LOSAT/tests/fasta_input_fixtures.py that `register` accepts (deflines, sequence lines,
  line ends, text before the first defline, records without residues, ...), FASTA_INPUT_SEARCHES
  per program, with NCBI's frozen hashes, and per program one search whose first query batch
  fails after the outfmt 0 prolog (`"expect": "failure"`: queries without residues only).
- quick: a few searches per program with repository paths, compared with the native
  CLI only (for CI).

Both suites add a TBLASTX search over several subjects, which the threaded reactor
reduces in the wasm-threads-only path. `--group fasta_input` writes that group of the full
suite alone.

Usage: v_abi_cases.py --suite quick|full [--group fasta_input] --out FILE.json
"""
from __future__ import annotations

import argparse
import json
import shlex
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "docs/evidence/losat_web_e1a"))
import capture_outputs as capture  # noqa: E402  (the verified case builders)

FORMATS = {"blastp": ["0", "6", "7"], "tblastn": ["0", "6", "7"], "blastn": ["0", "6", "7"], "tblastx": ["0", "6", "7"]}
ADAPTER_OWNED = {"-outfmt", "-out", "-num_threads"}


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
        outfmt, argv = split(case.command)
        outfmt = outfmt or "0"
        if outfmt not in FORMATS[case.program]:
            continue  # custom field lists are CLI-only
        key = (case.program, tuple(argv), str(case.cwd))
        search = searches.setdefault(key, {"program": case.program, "argv": argv, "cwd": str(case.cwd),
                                           "cases": [], "frozen": {}})
        search["cases"].append(case.case_id)
        if frozen and case.expected_sha256:
            previous = search["frozen"].setdefault(outfmt, case.expected_sha256)
            if previous != case.expected_sha256:
                raise SystemExit(f"{case.case_id}: two frozen hashes for outfmt {outfmt}")
    return list(searches.values())


def outfmt0_fixtures() -> list[dict]:
    """The NCBI-frozen outfmt 0 and 7 fixtures (docs/evidence/losat_web_e2a/, e2b/), run from
    LOSAT/: one search per argv, with the frozen hash of each fixture's format (none for an
    approved genetic-code deviation, whose oracle is NCBI's database search)."""
    sys.path.insert(0, str(REPO / "docs/evidence/losat_web_e2a"))
    import run_oracle  # noqa: E402  (the manifest reader and the oracle's arguments)
    searches: dict[tuple, dict] = {}
    for row in run_oracle.read_manifest()[2]:
        argv = run_oracle.search_argv(row)
        outfmt = run_oracle.outfmt(row)
        assert argv[-2:] == ["-outfmt", outfmt], row["fixture_id"]
        key = (row["program"], tuple(argv[:-2]))
        search = searches.setdefault(key, {"program": row["program"], "argv": [row["program"], *argv[:-2]],
                                           "cwd": str(run_oracle.ENGINE), "cases": [], "frozen": {}})
        search["cases"].append(row["fixture_id"])
        if row.get("contract") != run_oracle.DEVIATION:
            search["frozen"][outfmt] = row["stdout_sha256"]
    return list(searches.values())


FASTA_INPUT_SEARCHES = 32
# Queries without residues only: NCBI's first batch has no query to search, so the run fails
# after the outfmt 0 prolog (on a file sink: AUTHORITY.md of losat_web_e2h, S5 decision 5).
FIRST_BATCH_FAILURE = "@/rec_empty_all.q.{}"


def fasta_input() -> list[dict]:
    """SF fixture inputs through the ABI (gate NOTES open item 3): the rows of
    LOSAT/tests/fasta_input_fixtures.py that NCBI runs on files (no standard input, no
    environment, one thread) and completes in every format of the search, whose inputs
    `register` accepts (no `>?` gap line: maintainer decision 3), without the
    `-lcase_masking` rows that LOSAT's TBLASTX and BLASTP reject. One search per program,
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
        if "-lcase_masking" in args and program in ("tblastx", "blastp"):
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


def quick() -> list[dict]:
    cases = []
    blastp = capture.audit_blastp_v010.load_manifest(capture.TESTS / "blastp_v010_parity_manifest.tsv", REPO)
    for base in blastp[:2]:
        _, losat = capture.audit_blastp_v010.build_commands(base, Path("unused"), Path("LOSAT"), Path("{OUTDIR}"))
        cases.append(capture.Case("blastp", base.case_id, capture.strip_binary(losat), REPO))
    for row in capture.compare_blastn_parity.read_manifest(capture.TESTS / "blastn_parity_manifest.tsv"):
        if row["case_id"].startswith(("compact.", "MjPMNV.")):
            absolute = capture.certify._blastn_command_row(REPO, row)
            command = capture.compare_blastn_parity.build_losat_command(absolute, Path("LOSAT"), "{OUT}")
            cases.append(capture.Case("blastn", row["case_id"], capture.strip_binary(command), REPO))
    tblastx = capture.audit_tblastx_v010.load_manifest(capture.TESTS / "tblastx_v010_parity_manifest.tsv", REPO)
    for case in tblastx[1:2]:
        _, losat = capture.audit_tblastx_v010.build_commands(case, Path("unused"), Path("LOSAT"), Path("unused"), Path("{OUT}"))
        cases.append(capture.Case("tblastx", case.case_id, capture.strip_binary(losat), REPO))
    prefix = "/mnt/c/Users/genom/GitHub/LOSAT/"
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
    # BLASTN dc-megablast and blastn-short (Session SD): their default fixtures, with the
    # NCBI-frozen hashes of outfmt 0 and 7.
    searches += [search for search in outfmt0_fixtures()
                 if set(search["cases"]) & {"dc.default", "short.default"}]
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
    args.out.write_text(json.dumps(searches, indent=2) + "\n")
    print(f"{len(searches)} searches -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
