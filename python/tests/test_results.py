from __future__ import annotations

import multiprocessing
import pickle
from concurrent.futures import ProcessPoolExecutor

import pylosat
import pytest
from conftest import CASES, fasta
from pylosat import Hsp, HspTable, Outcome, Search, Sequences

DEFAULT_REPORT = (
    "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore"
)


def outcome_of(case, outfmt=()) -> Outcome:
    program, query, subject, options, _ = case
    return Search(program, **options).run(
        Sequences.from_fasta(fasta(query)), Sequences.from_fasta(fasta(subject)), outfmt=outfmt
    )


@pytest.fixture(scope="module")
def blastp() -> Outcome:
    return outcome_of(CASES[1], outfmt=("6", "0"))


def test_the_columns_start_with_the_default_report(blastp):
    columns = blastp.hsps.columns
    assert columns[:12] == DEFAULT_REPORT.split()
    assert "evalue" in columns and "e_value" not in columns


def test_a_table_is_a_sequence_of_hsps(blastp):
    hsps = blastp.hsps
    assert len(hsps) > 3
    assert isinstance(hsps[0], Hsp)
    assert hsps[-1] == hsps[len(hsps) - 1]
    assert [hsp.index for hsp in hsps] == list(range(len(hsps)))
    assert hsps[1:3][0] == hsps[1]
    assert len(hsps[1:3]) == 2
    assert len(hsps[5:5]) == 0
    with pytest.raises(IndexError):
        hsps[len(hsps)]
    with pytest.raises(ValueError, match="step 1"):
        hsps[::2]
    assert repr(hsps).startswith("<HspTable:")
    assert "query_id" in repr(hsps[0])


def test_the_views_hold_each_hsp_once_in_the_order_of_the_hit_list(blastp):
    seen = []
    for query in blastp.queries():
        assert query == blastp.query(query.q_idx)
        assert all(hsp.q_idx == query.q_idx for hsp in query.hsps)
        for subject in query.subjects():
            assert subject == query.subject(subject.s_idx)
            assert all(hsp.s_idx == subject.s_idx for hsp in subject.hsps)
            seen.extend(hsp.index for hsp in subject.hsps)
    assert seen == list(range(len(blastp.hsps)))
    assert blastp.queries()[0].subject(10**6) is None


def test_each_query_has_its_statistics(blastp):
    assert len(blastp.query_stats) == len(blastp.queries())
    stats = blastp.query_stats[0]
    assert stats.eff_space > 0
    assert stats.lambda_ > 0 and stats.kappa > 0 and stats.entropy > 0


def test_the_ranges_hold_the_letters_of_the_rows():
    case = CASES[0]
    records = Sequences.from_fasta(fasta(case[1]))
    outcome = Search("blastn", task="megablast").run(records, records)
    text = open(fasta(case[1])).read().split("\n", 1)[1].replace("\n", "")
    hsp = outcome.hsps[0]
    assert text[hsp.query_range].upper() == hsp.query_aligned.replace("-", "").upper()


def test_a_span_is_the_row_of_the_hsp(blastp):
    report = blastp.reports[0]
    rows = report.bytes.decode().splitlines()
    for hsp, row in zip(blastp.hsps, rows):
        assert report.bytes[report.spans[hsp.index]].decode().rstrip("\n") == row
        assert row.startswith(f"{hsp.query_id}\t{hsp.subject_id}\t")
    pairwise = blastp.reports[1]
    span = pairwise.spans[0]
    assert b"Score =" in pairwise.bytes[span]
    assert b"Length=" in pairwise.bytes[pairwise.subject_headings[0]]


def test_the_frames_exist_where_the_program_translates():
    blastp_row = outcome_of(CASES[1]).hsps[0]
    assert (blastp_row.query_frame, blastp_row.subject_frame) == (None, None)
    assert (blastp_row.query_strand, blastp_row.subject_strand) == (None, None)
    tblastx_row = outcome_of(CASES[4]).hsps[0]
    assert tblastx_row.query_frame is not None and tblastx_row.subject_frame is not None
    assert tblastx_row.sum_n is not None
    assert outcome_of(CASES[0]).hsps[0].subject_strand in ("plus", "minus")


@pytest.mark.parametrize("case", CASES, ids=[case[0] for case in CASES])
def test_results_survive_pickle(case):
    outcome = outcome_of(case, outfmt="0")
    assert pickle.loads(pickle.dumps(outcome)) == outcome
    assert pickle.loads(pickle.dumps(outcome.query(0))) == outcome.query(0)
    hsp = outcome.hsps[-1]
    assert pickle.loads(pickle.dumps(hsp)) == hsp
    part = outcome.hsps[1:3]
    assert pickle.loads(pickle.dumps(part)) == part


def test_a_pickled_row_does_not_carry_the_table(blastp):
    assert len(pickle.dumps(blastp.hsps[:1])) < len(pickle.dumps(blastp.hsps)) / 2


def test_a_table_of_no_hits():
    search = Search("blastp", evalue=1e-200)
    empty = search.run(Sequences.from_text(">q\nMKTAYIAKQR\n"), Sequences.from_text(">s\nWWWWWWWWWW\n"))
    assert len(empty.hsps) == 0
    assert list(empty.hsps) == []
    assert pickle.loads(pickle.dumps(empty)) == empty
    assert len(empty.queries()[0].subjects()) == 0


def test_bytes_that_are_no_table_are_refused():
    with pytest.raises(ValueError):
        HspTable(b"not a stream")


def search_in_another_process(fasta_text: str) -> Outcome:
    records = Sequences.from_text(fasta_text)
    return Search("blastn").run(records, records)


def test_a_process_pool_returns_results():
    text = open(fasta("small_test.fasta")).read()
    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(1, mp_context=context) as pool:
        outcome = pool.submit(search_in_another_process, text).result()
    assert outcome == search_in_another_process(text)


def test_arrow_consumers_read_the_table(blastp):
    pyarrow = pytest.importorskip("pyarrow")
    table = pyarrow.table(blastp.hsps)
    assert table.column_names == blastp.hsps.columns
    assert table.num_rows == len(blastp.hsps)
    assert table.column("evalue").to_pylist() == [hsp.e_value for hsp in blastp.hsps]
    assert table.column("bitscore").to_pylist() == [hsp.bit_score for hsp in blastp.hsps]
    assert table.column("qaccver").to_pylist() == [hsp.query_id for hsp in blastp.hsps]
    assert table.column("qframe").null_count == table.num_rows


def test_polars_reads_the_table(blastp):
    polars = pytest.importorskip("polars")
    frame = polars.DataFrame(blastp.hsps)
    assert frame.columns == blastp.hsps.columns
    assert frame["pident"].to_list() == [hsp.percent_identity for hsp in blastp.hsps]


def test_pandas_reads_the_table(blastp):
    pandas = pytest.importorskip("pandas", minversion="3.1")
    frame = pandas.DataFrame.from_arrow(blastp.hsps)
    assert list(frame.columns) == blastp.hsps.columns
    assert list(frame["length"]) == [hsp.length for hsp in blastp.hsps]


# The fields of `Hsp` and the columns that hold them.
FIELDS = {
    "query_id": "qaccver",
    "subject_id": "saccver",
    "percent_identity": "pident",
    "length": "length",
    "mismatch": "mismatch",
    "gap_open": "gapopen",
    "q_start": "qstart",
    "q_end": "qend",
    "s_start": "sstart",
    "s_end": "send",
    "e_value": "evalue",
    "bit_score": "bitscore",
    "raw_score": "score",
    "num_ident": "nident",
    "num_positives": "positive",
    "gaps": "gaps",
    "query_length": "qlen",
    "subject_length": "slen",
    "query_frame": "qframe",
    "subject_frame": "sframe",
    "subject_strand": "sstrand",
    "query_aligned": "qseq",
    "subject_aligned": "sseq",
    "query_strand": "query_strand",
    "sum_n": "sum_n",
    "comp_adjust_method": "comp_adjust_method",
    "index": "index",
    "q_idx": "q_idx",
    "s_idx": "s_idx",
    "rank": "rank",
}


@pytest.mark.parametrize("case", CASES, ids=[case[0] for case in CASES])
def test_each_field_of_an_hsp_reads_its_column(case):
    pyarrow = pytest.importorskip("pyarrow")
    hsps = outcome_of(case).hsps
    assert sorted(FIELDS.values()) == sorted(hsps.columns)
    rows = pyarrow.table(hsps).to_pylist()
    assert len(rows) > 0
    for hsp, row in zip(hsps, rows):
        for name, column in FIELDS.items():
            assert getattr(hsp, name) == row[column], name


def test_the_versions():
    assert pylosat.ENGINE_VERSION.count(".") == 2
