from __future__ import annotations

import bz2
import gzip
import pickle
import shutil
from pathlib import Path

import pytest
from conftest import fasta
from pylosat import InputError, Search, Sequences

RECORDS = [
    ("q1 a nucleotide query", "ACGTTGCAACGGATTACAGGCTTAGCATCGATCGGATCCTAGCTAGGATCGATTACG"),
    ("q2", "TTGACCGATAGCATCAGGCTAGGATCCATCGATTTAGCGGATCGAACTGGCATCGAT"),
]
FASTA_TEXT = "".join(f">{title}\n{sequence}\n" for title, sequence in RECORDS)


def hsps(query: Sequences, subject: Sequences):
    return Search("blastn", task="blastn-short").run(query, subject).hsps


def test_the_three_constructors_read_one_input(tmp_path: Path):
    path = tmp_path / "records.fa"
    path.write_text(FASTA_TEXT)
    subject = Sequences.from_text(FASTA_TEXT)
    expected = hsps(Sequences.from_fasta(path), subject)
    assert len(expected) > 0
    assert hsps(Sequences.from_text(FASTA_TEXT), subject) == expected
    assert hsps(Sequences.from_text(FASTA_TEXT.encode()), subject) == expected
    assert hsps(Sequences.from_pairs(RECORDS), subject) == expected


@pytest.mark.parametrize(
    ("suffix", "pack"), [(".gz", gzip.compress), (".bz2", bz2.compress), ("", bytes)]
)
def test_a_compressed_file_reads_as_on_the_command_line(tmp_path: Path, suffix, pack):
    path = tmp_path / f"records.fa{suffix}"
    path.write_bytes(pack(FASTA_TEXT.encode()))
    subject = Sequences.from_text(FASTA_TEXT)
    assert hsps(Sequences.from_fasta(path), subject) == hsps(Sequences.from_text(FASTA_TEXT), subject)


def test_a_file_that_does_not_open_is_the_command_lines_error(tmp_path: Path):
    missing = tmp_path / "none.fa"
    subject = Sequences.from_text(FASTA_TEXT)
    with pytest.raises(InputError, match='Argument "query". File is not accessible') as caught:
        hsps(Sequences.from_fasta(missing), subject)
    assert str(missing) in str(caught.value)


def test_an_input_survives_pickle():
    for sequences in (Sequences.from_text(FASTA_TEXT, name="n"), Sequences.from_fasta("f.fa")):
        assert pickle.loads(pickle.dumps(sequences)) == sequences


def test_a_line_end_in_a_pair_is_refused():
    with pytest.raises(ValueError, match="line end"):
        Sequences.from_pairs([("title", "ACGT\n>another\nACGT")])
    with pytest.raises(ValueError, match="line end"):
        Sequences.from_pairs([("two\nlines", "ACGT")])


def test_a_path_is_not_an_input_until_it_is_named_one():
    text = Sequences.from_text(FASTA_TEXT)
    search = Search("blastn")
    with pytest.raises(TypeError, match="from_fasta"):
        search.run(fasta("small_test.fasta"), text)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="subject"):
        search.batches(text, ">s\nACGT\n")  # type: ignore[arg-type]


def test_a_fasta_error_is_the_readers():
    with pytest.raises(InputError, match="BLAST query error"):
        hsps(Sequences.from_text("xyz{}[]|~`^&*()=+\n"), Sequences.from_text(FASTA_TEXT))


def test_the_name_shows_where_the_command_line_shows_its_file_name(tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    path = Path("subject.fa")
    shutil.copy(fasta("small_test.fasta"), path)
    query = Sequences.from_fasta(path)

    def report_with(subject: Sequences) -> str:
        # The report breaks a long line.
        return " ".join(Search("blastn").run(query, subject, outfmt="0").reports[0].text.split())

    assert f"sequence set (Input: {path})" in report_with(Sequences.from_fasta(path))
    assert "sequence set (Input: other)" in report_with(Sequences.from_fasta(path, name="other"))
    assert "sequence set (Input: subject)" in report_with(Sequences.from_text(path.read_text()))
