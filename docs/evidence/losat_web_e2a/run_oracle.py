#!/usr/bin/env python3
"""Run the NCBI oracle over LOSAT/tests/outfmt0_manifest.tsv (comparison only).

Every fixture runs from LOSAT/ as `<program> -query Q -subject S [-task T] EXTRA -outfmt F`
(F is the `outfmt` column, 0 when it is empty) with the output on stdout. The run writes
DIR/<fixture_id>.out and DIR/<fixture_id>.err (when stderr is not empty), and prints the
binaries, then one line per fixture (wall time, sizes, hashes, command); keep that as the
run's log.

A row whose `contract` is `approved_db_gencode_deviation` (a TBLASTX search with a
non-default `-db_gencode`, AGENTS.md's approved exception) also runs NCBI with the subject
as a BLAST database (`makeblastdb -dbtype nucl` without `-parse_seqids`, then `-db`
instead of `-subject`), whose search applies `-db_gencode` to the subject as LOSAT's local
search does; it writes DIR/<fixture_id>.db.out and checks `db_stdout_sha256`, the hash
of that output with its `Posted date:` line (the time makeblastdb ran) replaced by a
fixed text (check_losat.py compares LOSAT with that output outside the database lines).

A row whose `contract` is `approved_punct_title:<stand-in subject>` (approved exception 2
of PD-LOSAT-NCBI-DEFECTS: subject titles made only of punctuation, on which NCBI's outfmt 0
crashes when the subject has hits) runs NCBI with the stand-in subject instead: the same
records, with a placeholder defline (a run of one letter) of the length of the title that
NCBI's `x_CleanAndCompress` gives when it stops at the end of the string, for each defline
that differs. Its expected output is NCBI's report of the stand-in with each placeholder
replaced by that title and the stand-in's path by the subject's (the two paths have the
same length), the check of docs/evidence/losat_web_e2g/title_sweep.py and
docs/evidence/losat_web_e2e/title_sweep.py.

- default: compare the outputs with the manifest's hash columns; exit 1 on a difference.
- --freeze: write the hash and size columns of the manifest instead (run with
  --out LOSAT/tests/fixtures/outfmt0 so the fixtures are the frozen bytes).
- --only PREFIX: run (and freeze) only the fixtures whose id starts with PREFIX; the other
  rows of the manifest keep their columns.

NCBI reads these environment variables and .ncbirc, and each one changes the report
(AUTHORITY.md §A.9), so the run refuses to start when one is set or present.

Usage: run_oracle.py --bin-dir DIR --out DIR [--freeze] [--only PREFIX]
"""
from __future__ import annotations

import argparse
import hashlib
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
ENGINE = REPO / "LOSAT"
MANIFEST = ENGINE / "tests/outfmt0_manifest.tsv"
HASHES = ("stdout_sha256", "stdout_bytes", "stderr_sha256")
DEVIATION = "approved_db_gencode_deviation"
PUNCT_TITLE = "approved_punct_title:"
# BATCH_SIZE and CHUNK_SIZE change the query batches (blast_input_aux.cpp:86-90,
# local_blast.cpp:59-62), and with them the report of an invalid query.
REPORT_ENV = ("BL2SEQ_LEGACY", "CTOOLKIT_COMPATIBLE", "OLD_FSC", "BATCH_SIZE", "CHUNK_SIZE")


def read_manifest() -> tuple[list[str], list[str], list[dict[str, str]]]:
    lines = MANIFEST.read_text().splitlines()
    comments = [line for line in lines if line.startswith("#")]
    header, *body = [line.split("\t") for line in lines if not line.startswith("#")]
    return comments, header, [dict(zip(header, row)) for row in body]


def outfmt(row: dict[str, str]) -> str:
    """The row's output format (0 when the column is empty)."""
    return row.get("outfmt") or "0"


def search_argv(row: dict[str, str]) -> list[str]:
    """The search arguments after the program name (the same for NCBI and LOSAT)."""
    task = ["-task", row["task"]] if row["task"] else []
    return ["-query", row["query"], "-subject", row["subject"], *task, *shlex.split(row["extra_args"]),
            "-outfmt", outfmt(row)]


def db_argv(row: dict[str, str], db: str) -> list[str]:
    """The arguments of the database oracle of an approved-deviation row."""
    argv = search_argv(row)
    index = argv.index("-subject")
    return [*argv[:index], "-db", db, *argv[index + 2:]]


def punct_title_run(row: dict[str, str], binary: Path) -> subprocess.CompletedProcess:
    """NCBI's run of an `approved_punct_title` row: the stand-in subject's report with the
    placeholders replaced by the titles (see the module)."""
    sys.path.insert(0, str(REPO / "docs/evidence/losat_web_e2g"))
    from title_sweep import title  # noqa: E402  (NCBI's title stopped at the end of the string)
    stand_in = row["contract"].removeprefix(PUNCT_TITLE)
    if len(stand_in) != len(row["subject"]):
        raise SystemExit(f"{row['fixture_id']}: the stand-in path must have the subject path's length")
    deflines = [line[1:] for line in (ENGINE / row["subject"]).read_text().splitlines() if line.startswith(">")]
    placeholders = [line[1:] for line in (ENGINE / stand_in).read_text().splitlines() if line.startswith(">")]
    argv = search_argv(row)
    argv[argv.index("-subject") + 1] = stand_in
    run = subprocess.run([str(binary), *argv], cwd=ENGINE, capture_output=True)
    titles = {}
    for defline, placeholder in zip(deflines, placeholders, strict=True):
        if defline != placeholder:
            expected = title(defline)
            if not expected or placeholder != placeholder[0] * len(expected) or not placeholder.isalpha():
                raise SystemExit(f"{row['fixture_id']}: placeholder {placeholder!r} for {defline!r}")
            titles[placeholder.encode()] = expected.encode()
    # A title starts a line of the description table and follows "> " in an alignment
    # heading; only there is a placeholder replaced.
    lines = run.stdout.replace(stand_in.encode(), row["subject"].encode()).split(b"\n")
    for index, line in enumerate(lines):
        for placeholder, expected in titles.items():
            for prefix in (b"", b"> "):
                rest = line[len(prefix) + len(placeholder):]
                if line.startswith(prefix + placeholder) and (not rest or rest.startswith(b" ")):
                    lines[index] = prefix + expected + rest
    return subprocess.CompletedProcess(run.args, run.returncode, b"\n".join(lines), run.stderr)


def without_posted_date(report: bytes) -> bytes:
    """The report with the time of the database (outfmt 0 `Posted date:`) as a fixed text."""
    return re.sub(rb"(?m)^(    Posted date:  ).*$", rb"\1(the time makeblastdb ran)", report)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--freeze", action="store_true")
    parser.add_argument("--only", default="")
    args = parser.parse_args()
    found = [name for name in REPORT_ENV if name in os.environ]
    found += [str(path) for path in (ENGINE / ".ncbirc", Path.home() / ".ncbirc") if path.exists()]
    if "NCBI" in os.environ:
        found.append("NCBI")
    if found:
        raise SystemExit(f"report-changing NCBI configuration present: {found}")
    args.out.mkdir(parents=True, exist_ok=True)
    comments, header, rows = read_manifest()
    for program in sorted({row["program"] for row in rows}):
        binary = args.bin_dir / program
        version = subprocess.run([binary, "-version"], capture_output=True, text=True, check=True).stdout
        print(f"# {binary} sha256={digest(binary.read_bytes())} {' / '.join(version.splitlines())}")
    if any(row.get("contract") == DEVIATION for row in rows):
        binary = args.bin_dir / "makeblastdb"
        version = subprocess.run([binary, "-version"], capture_output=True, text=True, check=True).stdout
        print(f"# {binary} sha256={digest(binary.read_bytes())} {' / '.join(version.splitlines())}")

    records, differing = [], []
    selected = [row for row in rows if row["fixture_id"].startswith(args.only)]
    for row in selected:
        command = [str(args.bin_dir / row["program"]), *search_argv(row)]
        start = time.monotonic()
        if row.get("contract", "").startswith(PUNCT_TITLE):
            run = punct_title_run(row, args.bin_dir / row["program"])
        else:
            run = subprocess.run(command, cwd=ENGINE, capture_output=True)
        wall = time.monotonic() - start
        if run.returncode != 0:
            raise SystemExit(f"{row['fixture_id']}: exit {run.returncode}\n{run.stderr.decode()}")
        (args.out / f"{row['fixture_id']}.out").write_bytes(run.stdout)
        err = args.out / f"{row['fixture_id']}.err"
        if run.stderr:
            err.write_bytes(run.stderr)
        elif err.exists():
            err.unlink()
        observed = {"stdout_sha256": digest(run.stdout), "stdout_bytes": str(len(run.stdout)),
                    "stderr_sha256": digest(run.stderr) if run.stderr else ""}
        hashes = list(HASHES)
        if row.get("contract") == DEVIATION:
            # outfmt 7 prints the -db argument, so the database has a fixed path.
            db_dir = Path(tempfile.gettempdir()) / "losat_outfmt0_db" / row["fixture_id"]
            shutil.rmtree(db_dir, ignore_errors=True)
            db_dir.mkdir(parents=True)
            db = str(db_dir / "subject")
            subprocess.run([str(args.bin_dir / "makeblastdb"), "-in", row["subject"], "-dbtype", "nucl",
                            "-out", db], cwd=ENGINE, capture_output=True, check=True)
            db_run = subprocess.run([str(args.bin_dir / row["program"]), *db_argv(row, db)], cwd=ENGINE,
                                    capture_output=True, check=True)
            shutil.rmtree(db_dir)
            (args.out / f"{row['fixture_id']}.db.out").write_bytes(db_run.stdout)
            observed["db_stdout_sha256"] = digest(without_posted_date(db_run.stdout))
            hashes.append("db_stdout_sha256")
        if args.freeze:
            row.update(observed)
        elif any(row.get(key, "") != observed[key] for key in hashes):
            differing.append(row["fixture_id"])
        shown = shlex.join([row["program"], *search_argv(row)])
        records.append([row["fixture_id"], f"{wall:.2f}", *(observed[key] for key in HASHES), shown])

    for record in [["fixture_id", "wall_s", *HASHES, "command"], *records]:
        print("\t".join(record))
    if args.freeze:
        MANIFEST.write_text("\n".join([*comments, "\t".join(header),
                                       *("\t".join(row.get(key, "") for key in header) for row in rows)]) + "\n")
    print(f"# fixtures={len(selected)} differing={differing}")
    return 1 if differing else 0


if __name__ == "__main__":
    sys.exit(main())
