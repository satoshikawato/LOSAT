"""Byte-parity tests for `pylosat.blast_formatter` against the native `losat blast_formatter` CLI."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pylosat
import pytest
from conftest import FASTA

BF_DIR = FASTA.parent / "blast_formatter"
OUTFMTS = ("8", "9", "12")


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


def assert_tool_matches_cli(
    binary: str,
    words: tuple[str, ...],
    *,
    stdin: bytes = b"",
    input_data: bytes | str | None = None,
) -> subprocess.CompletedProcess[bytes]:
    cli = run_cli_in(binary, BF_DIR, *words, stdin=stdin)
    assert cli.returncode == 0, cli.stderr
    output = pylosat.blast_formatter(*words, input=input_data) if input_data is not None else pylosat.blast_formatter(*words)
    assert output.bytes == cli.stdout
    assert output.diagnostics_bytes == cli.stderr
    return cli


@pytest.mark.parametrize("outfmt", OUTFMTS)
def test_subject_and_database_archives_match_cli(outfmt, losat_bin, monkeypatch):
    monkeypatch.chdir(BF_DIR)
    archive = BF_DIR / "in/bn.asn"
    cli = assert_tool_matches_cli(losat_bin, ("-archive", str(archive), "-outfmt", outfmt))

    raw = archive.read_bytes()
    for payload in (raw, raw.decode("utf-8")):
        streamed = pylosat.blast_formatter("-archive", "-", "-outfmt", outfmt, input=payload)
        assert streamed.bytes == cli.stdout
        assert streamed.diagnostics_bytes == cli.stderr

    assert_tool_matches_cli(losat_bin, ("-archive", "in/db.asn", "-outfmt", outfmt))


def test_batch_options_and_out_file_match_cli(tmp_path, losat_bin, monkeypatch):
    monkeypatch.chdir(BF_DIR)
    for words in (
        ("-archive", "in/batch.asn", "-outfmt", "8"),
        ("-archive", "in/db.asn", "-outfmt", "8", "-max_target_seqs", "2"),
        ("-archive", "in/db.asn", "-outfmt", "8", "-num_alignments", "0"),
        ("-archive", "in/db.asn", "-outfmt", "8", "-parse_deflines"),
        ("-archive", "in/warn.asn", "-outfmt", "8", "-max_target_seqs", "2"),
    ):
        assert_tool_matches_cli(losat_bin, words)

    cli_out = tmp_path / "cli.out"
    py_out = tmp_path / "py.out"
    cli = run_cli_in(losat_bin, BF_DIR, "-archive", "in/bn.asn", "-outfmt", "12", "-out", str(cli_out))
    assert cli.returncode == 0
    written = pylosat.blast_formatter("-archive", "in/bn.asn", "-outfmt", "12", "-out", str(py_out))
    assert written.bytes == b""
    assert written.diagnostics_bytes == cli.stderr
    assert py_out.read_bytes() == cli_out.read_bytes()


def test_roundtrip_from_pylosat_search_outfmt_11(losat_bin):
    query = pylosat.Sequences.from_text(">q1\nACGTTGCAAGGCTTAACCGGTTAACCGGATCGATCGATCG\n")
    subject = pylosat.Sequences.from_text(">s1\nACGTTGCAAGGCTTAACCGGTTAACCGGATCGATCGATCG\n")
    archive = pylosat.Search("blastn").run(query, subject, outfmt="11").reports[0].bytes

    for outfmt in OUTFMTS:
        assert_tool_matches_cli(losat_bin, ("-archive", "-", "-outfmt", outfmt), stdin=archive, input_data=archive)


@pytest.mark.parametrize(
    "words",
    [
        (),
        ("-rid", ""),
        ("-rid", "ABCDEFGH01R"),
        ("-rid", "ABCDEFGH01R", "-archive", "in/bn.asn"),
        ("-archive", "in/bn.asn", "-outfmt", "5"),
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
