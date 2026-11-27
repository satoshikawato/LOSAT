"""Byte-parity tests for `pylosat.blast_formatter` against the native `losat blast_formatter` CLI."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pylosat
import pytest
from conftest import FASTA

BF_DIR = FASTA.parent / "blast_formatter"


@pytest.fixture(scope="session")
def losat_bin() -> str:
    executable = os.environ.get("LOSAT_BIN")
    if not executable:
        pytest.skip("LOSAT_BIN names no losat executable")
    return executable


def run_cli_in(binary: str, cwd: Path, *words: str, stdin: bytes = b"") -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        [binary, "blast_formatter", *words],
        cwd=cwd,
        input=stdin,
        capture_output=True,
        check=False,
    )


@pytest.mark.parametrize("outfmt", ["8", "9", "12"])
def test_subject_archive_file_and_stdin_match_cli(outfmt, losat_bin):
    archive = BF_DIR / "in/bn.asn"
    cli = run_cli_in(losat_bin, BF_DIR, "-archive", str(archive), "-outfmt", outfmt)
    assert cli.returncode == 0

    from_file = pylosat.blast_formatter("-archive", str(archive), "-outfmt", outfmt)
    assert from_file.bytes == cli.stdout
    assert from_file.diagnostics_bytes == cli.stderr

    raw = archive.read_bytes()
    from_stdin_bytes = pylosat.blast_formatter("-archive", "-", "-outfmt", outfmt, input=raw)
    assert from_stdin_bytes.bytes == cli.stdout
    assert from_stdin_bytes.diagnostics_bytes == cli.stderr

    from_stdin_str = pylosat.blast_formatter("-archive", "-", "-outfmt", outfmt, input=raw.decode("utf-8"))
    assert from_stdin_str.bytes == cli.stdout
    assert from_stdin_str.diagnostics_bytes == cli.stderr


@pytest.mark.parametrize("outfmt", ["8", "9", "12"])
def test_database_archive_matches_cli(outfmt, losat_bin, monkeypatch):
    monkeypatch.chdir(BF_DIR)
    cli = run_cli_in(losat_bin, BF_DIR, "-archive", "in/db.asn", "-outfmt", outfmt)
    assert cli.returncode == 0

    output = pylosat.blast_formatter("-archive", "in/db.asn", "-outfmt", outfmt)
    assert output.bytes == cli.stdout
    assert output.diagnostics_bytes == cli.stderr


def test_batch_options_and_out_file_match_cli(tmp_path, losat_bin, monkeypatch):
    monkeypatch.chdir(BF_DIR)
    for words in (
        ("-archive", "in/batch.asn", "-outfmt", "8"),
        ("-archive", "in/db.asn", "-outfmt", "8", "-max_target_seqs", "2"),
        ("-archive", "in/db.asn", "-outfmt", "8", "-num_alignments", "0"),
        ("-archive", "in/db.asn", "-outfmt", "8", "-parse_deflines"),
        ("-archive", "in/warn.asn", "-outfmt", "8", "-max_target_seqs", "2"),
    ):
        cli = run_cli_in(losat_bin, BF_DIR, *words)
        assert cli.returncode == 0, cli.stderr
        output = pylosat.blast_formatter(*words)
        assert output.bytes == cli.stdout
        assert output.diagnostics_bytes == cli.stderr

    cli_out = tmp_path / "cli.out"
    py_out = tmp_path / "py.out"
    cli = run_cli_in(losat_bin, BF_DIR, "-archive", "in/bn.asn", "-outfmt", "12", "-out", str(cli_out))
    assert cli.returncode == 0
    output = pylosat.blast_formatter("-archive", "in/bn.asn", "-outfmt", "12", "-out", str(py_out))
    assert output.bytes == b""
    assert output.diagnostics_bytes == cli.stderr
    assert py_out.read_bytes() == cli_out.read_bytes()


def test_roundtrip_from_pylosat_search_outfmt_11(losat_bin):
    query = pylosat.Sequences.from_text(">q1\nACGTTGCAAGGCTTAACCGGTTAACCGGATCGATCGATCG\n")
    subject = pylosat.Sequences.from_text(">s1\nACGTTGCAAGGCTTAACCGGTTAACCGGATCGATCGATCG\n")
    archive = pylosat.Search("blastn").run(query, subject, outfmt="11").reports[0].bytes

    for outfmt in ("8", "9", "12"):
        cli = run_cli_in(losat_bin, BF_DIR, "-archive", "-", "-outfmt", outfmt, stdin=archive)
        assert cli.returncode == 0, cli.stderr
        output = pylosat.blast_formatter("-archive", "-", "-outfmt", outfmt, input=archive)
        assert output.bytes == cli.stdout
        assert output.diagnostics_bytes == cli.stderr


@pytest.mark.parametrize(
    "words",
    [
        (),
        ("-rid", ""),
        ("-rid", "ABCDEFGH01R"),
        ("-rid", "ABCDEFGH01R", "-archive", "in/bn.asn"),
        ("-archive", "in/bn.asn", "-outfmt", "6"),
        ("-archive", "in/empty.asn", "-outfmt", "8"),
        ("-archive", "in/cut2.asn", "-outfmt", "8"),
        ("-archive", "in/no_such_file.asn", "-outfmt", "8"),
    ],
)
def test_errors_and_refusals_match_cli(words, losat_bin, monkeypatch):
    monkeypatch.chdir(BF_DIR)
    cli = run_cli_in(losat_bin, BF_DIR, *words)
    assert cli.returncode != 0
    with pytest.raises(pylosat.LosatError) as raised:
        pylosat.blast_formatter(*words)
    assert raised.value.message_bytes == cli.stderr.rstrip(b"\n")
