"""The binding against the command line: the reports and the records of the same search."""

from __future__ import annotations

import pytest
from conftest import CASES, fasta
from pylosat import Search, Sequences

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
