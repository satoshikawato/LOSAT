"""A BLAST database as the subject, reports in files, and the two database programs."""

from __future__ import annotations

import gzip
from pathlib import Path

import pylosat
import pytest
from conftest import CASES, FASTA, fasta

ENGINE_TESTS = FASTA.parent
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ"
FASTA_TEXT = f">one a protein\n{PROTEIN}\n>two the same protein\n{PROTEIN}A\n"

# (program, query, database, options, the same options as words), the first two under LOSAT/tests.
DATABASES = [
    ("blastn", "seqdb/nucl.fna", "seqdb/db/nucl_v5", {"task": "blastn"}, ["-task", "blastn"]),
    ("blastp", "seqdb/prot.faa", "seqdb/db/prot_v5", {}, []),
    ("blastx", "db_search/bslice_q.fna", "db_search/bslice", {}, []),
    ("tblastn", "db_search/tslice_q.faa", "db_search/tslice", {}, []),
    ("tblastx", "db_search/xslice_q.fna", "db_search/xslice", {}, []),
    ("psiblast", "fasta/psiblast/tk.faa", "db_search/paja", {}, []),
]


def engine_file(name: str) -> str:
    return str(ENGINE_TESTS / name)


@pytest.mark.parametrize("case", DATABASES, ids=lambda case: case[0])
@pytest.mark.parametrize("outfmt", ["6", "7", "0"])
def test_a_database_gives_the_bytes_of_the_command_line(case, outfmt, losat):
    program, query, database, options, words = case
    query, database = engine_file(query), engine_file(database)
    stdout, stderr = losat(program, "-query", query, "-db", database, *words, "-outfmt", outfmt)
    outcome = pylosat.Search(program, **options).run(
        pylosat.Sequences.from_fasta(query), pylosat.Database(database), outfmt=outfmt
    )
    assert outcome.reports[0].bytes == stdout
    assert outcome.diagnostics == stderr.decode()


@pytest.mark.parametrize("case", DATABASES, ids=lambda case: case[0])
def test_a_database_gives_the_rows_of_its_report(case):
    program, query, database, options, _ = case
    outcome = pylosat.Search(program, **options).run(
        pylosat.Sequences.from_fasta(engine_file(query)),
        pylosat.Database(engine_file(database)),
        outfmt="6",
    )
    rows = outcome.reports[0].text.splitlines()
    assert rows
    assert len(rows) == len(outcome.hsps)
    for row, hsp in zip(rows, outcome.hsps, strict=True):
        fields = row.split("\t")
        assert (fields[0], fields[1]) == (hsp.query_id, hsp.subject_id)
        assert int(fields[9]) == hsp.s_end


# Series Q item Q6b: a PSSM query (`-in_pssm` of psiblast and tblastn), against a database and
# against records: (program, PSSM under LOSAT/tests, subject, whether it is a database).
PSSM_QUERIES = [
    ("psiblast", "fixtures/pssm/ps4_tk_paja.asn.gz", "db_search/paja", True),
    ("psiblast", "fixtures/pssm/ps4_tk_paja.asn.gz", "fasta/PajaWSV.faa", False),
    ("tblastn", "fixtures/pssm/ps7_ts1_tslice.asn.gz", "db_search/tslice", True),
    ("tblastn", "db_search/of1/tnq1.pssm", "db_search/of1/tns.fa", False),
]


def pssm_file(name: str, directory: Path) -> str:
    """The PSSM file under LOSAT/tests, a `.gz` one inflated into `directory`."""
    path = ENGINE_TESTS / name
    if path.suffix != ".gz":
        return str(path)
    inflated = directory / path.stem
    inflated.write_bytes(gzip.decompress(path.read_bytes()))
    return str(inflated)


@pytest.mark.parametrize("case", PSSM_QUERIES, ids=lambda case: f"{case[0]}-{case[2]}")
@pytest.mark.parametrize("outfmt", ["6", "0"])
def test_a_pssm_query_gives_the_bytes_of_the_command_line(case, outfmt, losat, tmp_path):
    program, pssm, subject, database = case
    pssm, subject = pssm_file(pssm, tmp_path), engine_file(subject)
    flag = "-db" if database else "-subject"
    stdout, stderr = losat(program, "-in_pssm", pssm, flag, subject, "-outfmt", outfmt)
    inputs = pylosat.Database(subject) if database else pylosat.Sequences.from_fasta(subject)
    outcome = pylosat.Search(program).run(pylosat.Pssm.from_file(pssm), inputs, outfmt=outfmt)
    assert outcome.reports[0].bytes == stdout
    assert outcome.diagnostics == stderr.decode()


@pytest.mark.parametrize("case", PSSM_QUERIES, ids=lambda case: f"{case[0]}-{case[2]}")
def test_a_pssm_query_is_query_0_of_the_rows(case, tmp_path):
    program, pssm, subject, database = case
    pssm, subject = pssm_file(pssm, tmp_path), engine_file(subject)
    inputs = pylosat.Database(subject) if database else pylosat.Sequences.from_fasta(subject)
    search = pylosat.Search(program)
    outcome = search.run(pylosat.Pssm.from_file(pssm), inputs, outfmt="6")
    rows = outcome.reports[0].text.splitlines()
    assert rows and len(rows) == len(outcome.hsps)
    for row, hsp in zip(rows, outcome.hsps, strict=True):
        assert row.split("\t")[0] == hsp.query_id
        assert hsp.q_idx == 0
    assert len(outcome.query_stats) == 1 and outcome.query_stats[0] is not None
    # The file of `run_to_file` holds the same report.
    out = tmp_path / "out.6"
    search.run_to_file(pylosat.Pssm.from_file(Path(pssm)), inputs, outfmt="6", path=out)
    assert out.read_bytes() == outcome.reports[0].bytes
    # One batch: NCBI reads no query batches for a PSSM query.
    batches = search.batches(pylosat.Pssm.from_file(pssm), inputs)
    assert batches.query_count == 1
    [batch] = list(batches)
    assert batch.query_indices == range(0, 1)
    assert len(batch.hsps) == len(outcome.hsps)


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


@pytest.mark.parametrize("case", CASES, ids=lambda case: case[0])
def test_a_report_in_a_file_has_its_side_files(case, tmp_path, command_line):
    program, query, subject, options, words = case
    out = tmp_path / "report.out"
    query, subject = fasta(query), fasta(subject)
    pylosat.Search(program, **options).run_to_file(
        pylosat.Sequences.from_fasta(query),
        pylosat.Sequences.from_fasta(subject),
        outfmt="14",
        path=out,
    )
    written = {path.name: path.read_bytes() for path in tmp_path.iterdir()}
    assert len(written) > 1
    for path in tmp_path.iterdir():
        path.unlink()
    command_line(program, query, subject, [*words, "-outfmt", "14", "-out", str(out)])
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


def test_a_tool_is_the_command_line(losat):
    name = engine_file("seqdb/db/prot_v5")
    for words in (
        ("-db", name, "-info"),
        ("-db", name, "-entry", "all"),
        ("-db", name, "-entry", "all", "-outfmt", "%o %l %t"),
    ):
        stdout, stderr = losat("blastdbcmd", *words)
        output = pylosat.blastdbcmd(*words)
        assert output.bytes == stdout
        assert output.diagnostics == stderr.decode()


def test_a_tool_reads_its_input_as_the_command_line_reads_standard_input(losat):
    words = ("-db", engine_file("blastdbcmd/db/p5"), "-entry_batch", "-", "-outfmt", "%a")
    batch = pylosat.blastdbcmd(*words, input="P99999.2\n")
    assert batch.text == "XP_000001.1\nP99999.2\nAAA00001.1\n"
    fasta = engine_file("dustmasker/in/basic.fa")
    with open(fasta, "rb") as file:
        data = file.read()
    for words in ((), ("-outfmt", "fasta"), ("-outfmt", "maskinfo_asn1_text")):
        stdout, stderr = losat("dustmasker", "-in", fasta, *words)
        output = pylosat.dustmasker(*words, input=data)
        assert output.bytes == stdout
        assert output.diagnostics == stderr.decode()
    fasta = engine_file("segmasker/in/basic.fa")
    with open(fasta, "rb") as file:
        data = file.read()
    for words in ((), ("-outfmt", "fasta"), ("-outfmt", "maskinfo_asn1_text", "-locut", "2e")):
        stdout, stderr = losat("segmasker", "-in", fasta, *words)
        output = pylosat.segmasker(*words, input=data)
        assert output.bytes == stdout
        assert output.diagnostics == stderr.decode()
    counts = engine_file("windowmasker/in/nopct.ascii")
    fasta = engine_file("windowmasker/in/multi.fa")
    with open(fasta, "rb") as file:
        data = file.read()
    for words in (("-ustat", counts), ("-ustat", counts, "-outfmt", "fasta", "-dust", "T")):
        stdout, stderr = losat("windowmasker", "-in", fasta, *words)
        output = pylosat.windowmasker(*words, input=data)
        assert output.bytes == stdout
        assert output.diagnostics == stderr.decode()
    words = ("-masking_algorithm", "dust", "-masking_options", "o", "-outfmt", "interval")
    output = pylosat.convert2blastmask(*words, input=">a\nACGTacgt\n")
    assert output.bytes == b"4 - 7\n"
    assert output.diagnostics == ""


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
