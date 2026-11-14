"""The binding against the command line: the reports and the records of the same search."""

from __future__ import annotations

import pytest
from conftest import CASES, fasta
from pylosat import Search, SearchError, Sequences

COLUMNS = (
    "qaccver saccver length mismatch gapopen qstart qend sstart send score nident positive "
    "gaps qlen slen sstrand pident qframe sframe"
)


def search_of(case) -> Search:
    program, _, _, options, _ = case
    return Search(program, **options)


def inputs(case) -> tuple[Sequences, Sequences]:
    return Sequences.from_fasta(fasta(case[1])), Sequences.from_fasta(fasta(case[2]))


@pytest.mark.parametrize("case", CASES, ids=[case[0] for case in CASES])
@pytest.mark.parametrize("outfmt", ["6", "7", "0", "10", "6 qseqid sseqid evalue pident"])
def test_a_report_is_the_bytes_of_the_command_line(case, outfmt, command_line):
    program, query, subject, _, words = case
    stdout, stderr = command_line(program, fasta(query), fasta(subject), [*words, "-outfmt", outfmt])
    outcome = search_of(case).run(*inputs(case), outfmt=outfmt)
    assert outcome.reports[0].outfmt == outfmt
    assert outcome.reports[0].bytes == stdout
    assert outcome.diagnostics.encode() == stderr


@pytest.mark.parametrize("case", CASES, ids=[case[0] for case in CASES])
def test_the_table_holds_the_columns_of_the_command_line(case, command_line):
    program, query, subject, _, words = case
    stdout, _ = command_line(
        program, fasta(query), fasta(subject), [*words, "-outfmt", f"6 {COLUMNS}"]
    )
    rows = [line.split("\t") for line in stdout.decode().splitlines()]
    outcome = search_of(case).run(*inputs(case))
    assert len(outcome.hsps) == len(rows)
    assert len(rows) > 0
    for hsp, row in zip(outcome.hsps, rows):
        strand = hsp.subject_strand or "N/A"
        frames = row[-2:]
        del row[-2:]
        # Printed for a side that the program does not translate as well (`README.md`).
        if hsp.query_frame is not None:
            assert str(hsp.query_frame) == frames[0]
        if hsp.subject_frame is not None:
            assert str(hsp.subject_frame) == frames[1]
        assert [
            hsp.query_id,
            hsp.subject_id,
            *map(
                str,
                (
                    hsp.length,
                    hsp.mismatch,
                    hsp.gap_open,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                    hsp.raw_score,
                    hsp.num_ident,
                    hsp.num_positives,
                    hsp.gaps,
                    hsp.query_length,
                    hsp.subject_length,
                ),
            ),
            strand,
            f"{hsp.percent_identity:.3f}",
        ] == row


def test_the_keywords_mean_what_the_words_mean(command_line):
    query, subject = fasta("outfmt0/e2e_many_query.faa"), fasta("outfmt0/e2e_many_subject.faa")
    stdout, _ = command_line(
        "blastp", query, subject, ["-evalue", "0.001", "-max_target_seqs", "3", "-outfmt", "6"]
    )
    search = Search("blastp", evalue=0.001, max_target_seqs=3)
    outcome = search.run(Sequences.from_fasta(query), Sequences.from_fasta(subject), outfmt="6")
    assert outcome.reports[0].bytes == stdout


# The NCBI C++ Toolkit's `-logfile`, `-conffile` and `-dryrun` run in each search as on the command
# line (line RG item E1w-2): `losat`'s files and streams are the binding's.

TOOLKIT_CASES = [
    # tblastx sends its run's messages to the `-logfile` file, the other programs to stderr.
    ("tblastx", "outfmt0/tblastx_many_query.fasta", "outfmt0/tblastx_many_subject.fasta"),
    ("blastp", "outfmt0/e2e_many_query.faa", "outfmt0/e2e_many_subject.faa"),
]


@pytest.mark.parametrize("case", TOOLKIT_CASES, ids=[case[0] for case in TOOLKIT_CASES])
def test_the_logfile_is_the_command_line_s(case, command_line, tmp_path):
    program, query, subject = case
    words = ["-max_target_seqs", "1", "-outfmt", "6"]
    cli_log, library_log = tmp_path / "cli.log", tmp_path / "library.log"
    for log in (cli_log, library_log):
        log.write_bytes(b"OLD LINE\n")
    stdout, stderr = command_line(
        program, fasta(query), fasta(subject), [*words, "-logfile", str(cli_log)]
    )
    search = Search.from_argv(program, ["-max_target_seqs", "1", "-logfile", str(library_log)])
    assert not search.dry_run
    outcome = search.run(
        Sequences.from_fasta(fasta(query)), Sequences.from_fasta(fasta(subject)), outfmt="6"
    )
    assert outcome.reports[0].bytes == stdout
    assert outcome.diagnostics.encode() == stderr
    assert library_log.read_bytes() == cli_log.read_bytes()
    warning = f"Warning: [{program}] Examining 5 or more matches is recommended\n".encode()
    if program == "tblastx":
        assert (stderr, cli_log.read_bytes()) == (b"", b"OLD LINE\n" + warning)
    else:
        assert (stderr, cli_log.read_bytes()) == (warning, b"OLD LINE\n")


def test_a_dry_run_searches_nothing_and_reads_no_input(command_line, tmp_path):
    missing = str(tmp_path / "missing.faa")
    log = tmp_path / "dry.log"
    words = ["-evalue", "1", "-dryrun", "-logfile", str(log)]
    stdout, stderr = command_line("blastp", missing, missing, words)
    assert (stdout, stderr, log.read_bytes()) == (b"", b"", b"")
    log.unlink()
    search = Search.from_argv("blastp", words)
    assert search.dry_run
    assert not Search("blastp", evalue=1.0).dry_run
    inputs = (Sequences.from_fasta(missing), Sequences.from_fasta(missing))
    outcome = search.run(*inputs, outfmt=["6", "0"])
    assert [report.outfmt for report in outcome.reports] == ["6", "0"]
    assert [report.bytes for report in outcome.reports] == [b"", b""]
    assert len(outcome.hsps) == 0
    assert outcome.diagnostics == ""
    assert log.read_bytes() == b""
    out = tmp_path / "out.txt"
    assert search.run_to_file(*inputs, outfmt="6", path=out) == ""
    assert not out.exists()
    batches = list(search.batches(*inputs))
    assert len(batches) == 1
    assert batches[0].query_indices == range(0, 0)
    assert len(batches[0].hsps) == 0
    assert batches[0].diagnostics == ""


def test_a_missing_conffile_is_ncbi_s_critical_post(tmp_path):
    missing = tmp_path / "missing.ini"
    search = Search.from_argv("blastp", ["-conffile", str(missing)])
    inputs = (
        Sequences.from_fasta(fasta("outfmt0/e2e_many_query.faa")),
        Sequences.from_fasta(fasta("outfmt0/e2e_many_subject.faa")),
    )
    critical = f'Critical: (106.23) Registry file "{missing}" not found, aborting.'
    with pytest.raises(SearchError) as error:
        search.run(*inputs, outfmt="6")
    assert str(error.value) == critical
    # A file that loads changes nothing.
    conffile = tmp_path / "c.ini"
    conffile.write_text("[BLAST]\n")
    loaded = Search.from_argv("blastp", ["-conffile", str(conffile)]).run(*inputs, outfmt="6")
    assert loaded.reports == Search("blastp").run(*inputs, outfmt="6").reports
