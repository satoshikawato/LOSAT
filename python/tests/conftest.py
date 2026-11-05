"""Fixtures of the binding's tests: LOSAT's own FASTA files and, for the comparisons with the
command line, the `losat` executable named by `LOSAT_BIN`."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

FASTA = Path(__file__).resolve().parents[2] / "LOSAT" / "tests" / "fasta"

# (program, query, subject, options as `Search` takes them, the same options as words)
CASES = [
    (
        "blastn",
        "small_test.fasta",
        "small_test.fasta",
        {"task": "megablast"},
        ["-task", "megablast"],
    ),
    (
        "blastp",
        "outfmt0/e2e_many_query.faa",
        "outfmt0/e2e_many_subject.faa",
        {"evalue": 1e-5},
        ["-evalue", "1e-5"],
    ),
    (
        "blastx",
        "outfmt0/sml_bx_split_q.fna",
        "outfmt0/sml_bp_split_s.faa",
        {"num_threads": 2},
        ["-num_threads", "2"],
    ),
    (
        "tblastn",
        "outfmt0/e2e_many_query.faa",
        "outfmt0/tblastx_many_subject.fasta",
        {"evalue": 10.0, "db_gencode": 11},
        ["-evalue", "10", "-db_gencode", "11"],
    ),
    (
        "tblastx",
        "outfmt0/tblastx_many_query.fasta",
        "outfmt0/tblastx_many_subject.fasta",
        {"query_gencode": 11, "db_gencode": 11},
        ["-query_gencode", "11", "-db_gencode", "11"],
    ),
]


def fasta(name: str) -> str:
    return str(FASTA / name)


@pytest.fixture(scope="session")
def command_line():
    """Runs `losat <program> -query <query> -subject <subject> <words>`: stdout and stderr."""
    executable = os.environ.get("LOSAT_BIN")
    if not executable:
        pytest.skip("LOSAT_BIN names no losat executable")

    def run(program: str, query: str, subject: str, words: list[str]) -> tuple[bytes, bytes]:
        command = [executable, program, "-query", query, "-subject", subject, *words]
        done = subprocess.run(command, capture_output=True, check=True)
        return done.stdout, done.stderr

    return run


@pytest.fixture(scope="session")
def command_line_files():
    """Runs `losat <program> -query <query> -subject <subject> <words>`, which writes files."""
    executable = os.environ.get("LOSAT_BIN")
    if not executable:
        pytest.skip("LOSAT_BIN names no losat executable")

    def run(program: str, query: str, subject: str, words: list[str]) -> None:
        command = [executable, program, "-query", query, "-subject", subject, *words]
        subprocess.run(command, capture_output=True, check=True)

    return run


@pytest.fixture(scope="session")
def many_batches(tmp_path_factory) -> tuple[str, str]:
    """A BLASTP query file of several query batches, and its subject."""
    records = (FASTA / "outfmt0" / "e2e_many_query.faa").read_text()
    directory = tmp_path_factory.mktemp("batches")
    query = directory / "query.faa"
    query.write_text(records * max(4, 30_000 // len(records)))
    return str(query), fasta("outfmt0/e2e_many_subject.faa")
