"""A BLAST database as the subject, reports in files, and the two database programs."""

from __future__ import annotations

from pathlib import Path

import pylosat
import pytest

from conftest import CASES, FASTA, fasta

ENGINE_TESTS = FASTA.parent
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ"
FASTA_TEXT = f">one a protein\n{PROTEIN}\n>two the same protein\n{PROTEIN}A\n"

# (program, query, database, options) under LOSAT/tests.
DATABASES = [
    ("blastn", "seqdb/nucl.fna", "seqdb/db/nucl_v5", {"task": "blastn"}),
    ("blastp", "seqdb/prot.faa", "seqdb/db/prot_v5", {}),
    ("blastx", "db_search/bslice_q.fna", "db_search/bslice", {}),
    ("tblastn", "db_search/tslice_q.faa", "db_search/tslice", {}),
    ("tblastx", "db_search/xslice_q.fna", "db_search/xslice", {}),
    ("psiblast", "fasta/psiblast/tk.faa", "db_search/paja", {}),
]


def engine_file(name: str) -> str:
    return str(ENGINE_TESTS / name)


@pytest.mark.parametrize(("program", "query", "database", "options"), DATABASES)
def test_a_database_gives_the_rows_of_its_report(program, query, database, options):
    search = pylosat.Search(program, **options)
    query = pylosat.Sequences.from_fasta(engine_file(query))
    outcome = search.run(query, pylosat.Database(engine_file(database)), outfmt="6")
    rows = outcome.reports[0].text.splitlines()
    assert rows
    assert len(rows) == len(outcome.hsps)
    for row, hsp in zip(rows, outcome.hsps, strict=True):
        fields = row.split("\t")
        assert (fields[0], fields[1]) == (hsp.query_id, hsp.subject_id)
        assert int(fields[9]) == hsp.s_end


def test_batches_of_a_database_are_the_whole_search():
    search = pylosat.Search("blastp")
    query = pylosat.Sequences.from_fasta(engine_file("seqdb/prot.faa"))
    database = pylosat.Database(engine_file("seqdb/db/prot_v5"))
    whole = search.run(query, database)
    rows = [hsp.s_idx for batch in search.batches(query, database) for hsp in batch.hsps]
    assert rows == [hsp.s_idx for hsp in whole.hsps]


def test_a_missing_database_is_the_command_lines_error(tmp_path):
    search = pylosat.Search("blastp")
    query = pylosat.Sequences.from_text(f">q\n{PROTEIN}\n")
    with pytest.raises(pylosat.SearchError, match="No alias or index file found"):
        search.run(query, pylosat.Database(tmp_path / "none"))


def test_a_database_takes_a_path_and_only_a_database_takes_a_taxid_list(tmp_path):
    assert pylosat.Database(Path("a") / "b").name == str(Path("a") / "b")
    search = pylosat.Search("blastp", taxids="9606")
    sequences = pylosat.Sequences.from_text(f">q\n{PROTEIN}\n")
    with pytest.raises(pylosat.ArgumentError, match="-taxids"):
        search.run(sequences, sequences)


@pytest.mark.parametrize(("program", "query", "subject", "options", "words"), CASES)
def test_a_report_in_a_file_has_its_side_files(
    program, query, subject, options, words, tmp_path, command_line_files
):
    out = tmp_path / "report.out"
    search = pylosat.Search(program, **options)
    query_file, subject_file = fasta(query), fasta(subject)
    search.run_to_file(
        pylosat.Sequences.from_fasta(query_file),
        pylosat.Sequences.from_fasta(subject_file),
        outfmt="14",
        path=out,
    )
    written = {path.name: path.read_bytes() for path in tmp_path.iterdir()}
    assert len(written) > 1
    for path in tmp_path.iterdir():
        path.unlink()
    command_line_files(program, query_file, subject_file, [*words, "-outfmt", "14", "-out", str(out)])
    assert {path.name: path.read_bytes() for path in tmp_path.iterdir()} == written


def test_a_database_made_here_is_read_and_searched_here(tmp_path):
    name = str(tmp_path / "made")
    made = pylosat.makeblastdb("-dbtype", "prot", "-title", "made", "-out", name, input=FASTA_TEXT)
    assert "added 2 sequences" in made.text
    assert made.diagnostics == ""
    entries = pylosat.blastdbcmd("-db", name, "-entry", "all", "-outfmt", "%o %l")
    assert entries.text == f"0 {len(PROTEIN)}\n1 {len(PROTEIN) + 1}\n"
    assert isinstance(entries.bytes, bytes)
    outcome = pylosat.Search("blastp", max_hsps=1).run(
        pylosat.Sequences.from_text(FASTA_TEXT), pylosat.Database(name)
    )
    assert sorted((hsp.q_idx, hsp.s_idx) for hsp in outcome.hsps) == [
        (0, 0),
        (0, 1),
        (1, 0),
        (1, 1),
    ]


def test_a_tool_fails_with_the_command_lines_text(tmp_path):
    with pytest.raises(pylosat.SearchError, match="No alias or index file found"):
        pylosat.blastdbcmd("-db", str(tmp_path / "none"), "-info")
    with pytest.raises(pylosat.ArgumentError):
        pylosat.makeblastdb("-dbtype", "neither")


def test_psiblast_is_a_program():
    assert "psiblast" in pylosat.PROGRAMS
    flags = {parameter.flag for parameter in pylosat.describe("psiblast").parameters}
    assert "-inclusion_ethresh" in flags
    subject = pylosat.Sequences.from_text(f">s\n{PROTEIN}\n")
    outcome = pylosat.Search("psiblast").run(pylosat.Sequences.from_text(f">q\n{PROTEIN}\n"), subject)
    assert len(outcome.hsps) == 1
