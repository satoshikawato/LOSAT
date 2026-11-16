from __future__ import annotations

import _thread
import itertools
import threading
import time

import pytest
from conftest import CASES, fasta
from pylosat import Batch, Search, Sequences


def inputs(pair: tuple[str, str]) -> tuple[Sequences, Sequences]:
    return Sequences.from_fasta(pair[0]), Sequences.from_fasta(pair[1])


@pytest.mark.parametrize("case", CASES, ids=[case[0] for case in CASES])
def test_the_batches_one_after_the_other_are_the_whole_search(case):
    program, query, subject, options, _ = case
    search = Search(program, **options)
    records = Sequences.from_fasta(fasta(query)), Sequences.from_fasta(fasta(subject))
    whole = search.run(*records, outfmt="6")
    batches = search.batches(*records)
    found = list(batches)
    hsps = [hsp for batch in found for hsp in batch.hsps]
    assert hsps == list(whole.hsps)
    assert [stats for batch in found for stats in batch.query_stats] == whole.query_stats
    assert "".join(batch.diagnostics for batch in found) == whole.diagnostics
    assert found[-1].query_indices.stop == batches.query_count == len(whole.query_stats)


def test_a_query_file_of_several_batches(many_batches):
    search = Search("blastp", evalue=1e-5)
    records = inputs(many_batches)
    batches = search.batches(*records)
    found = list(batches)
    assert len(found) > 1
    position = 0
    for batch in found:
        assert batch.query_indices.start == position
        assert len(batch.query_stats) == len(batch.query_indices)
        assert [query.q_idx for query in batch.queries()] == list(batch.query_indices)
        position = batch.query_indices.stop
    assert position == batches.query_count
    hsps = [hsp for batch in found for hsp in batch.hsps]
    assert hsps == list(search.run(*records).hsps)
    assert [hsp.index for hsp in hsps] == list(range(len(hsps)))


def test_a_caller_stops_by_leaving_the_loop(many_batches):
    batches = Search("blastp").batches(*inputs(many_batches))
    first = next(batches)
    assert isinstance(first, Batch)
    assert first.query_indices.stop < batches.query_count
    del batches


def test_a_batch_that_fails_ends_the_search():
    from pylosat import SearchError

    search = Search("blastp")
    batches = search.batches(Sequences.from_text(">q\nMKTAYIAKQR\n"), Sequences.from_text(""))
    with pytest.raises(SearchError):
        next(batches)
    with pytest.raises(StopIteration):
        next(batches)


def test_ctrl_c_stops_between_two_batches(many_batches):
    batches = Search("blastp").batches(*inputs(many_batches))
    seen = 0
    with pytest.raises(KeyboardInterrupt):
        for _ in batches:
            seen += 1
            _thread.interrupt_main()
            # The next line runs or not: the interrupt comes at the latest with the next batch.
    assert seen == 1


def test_ctrl_c_stops_a_search_without_reports(many_batches):
    search = Search("blastp")
    records = inputs(many_batches)
    started = time.monotonic()
    search.run(*records)
    whole = time.monotonic() - started
    timer = threading.Timer(whole / 10, _thread.interrupt_main)
    timer.start()
    started = time.monotonic()
    with pytest.raises(KeyboardInterrupt):
        search.run(*records)
    assert time.monotonic() - started < whole * 0.8
    timer.join()


def test_python_runs_during_a_search(many_batches):
    counter = itertools.count()
    done = threading.Event()

    def count() -> None:
        while not done.is_set():
            next(counter)
            time.sleep(0.001)

    thread = threading.Thread(target=count)
    thread.start()
    Search("blastp").run(*inputs(many_batches), outfmt="6")
    done.set()
    thread.join()
    assert next(counter) > 10


def test_two_threads_search_at_once(many_batches):
    search = Search("blastp", evalue=1e-5, num_threads=2)
    records = inputs(many_batches)
    results: list[object] = [None, None]

    def run(position: int) -> None:
        results[position] = search.run(*records, outfmt="6")

    threads = [threading.Thread(target=run, args=(position,)) for position in range(2)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert results[0] == results[1]
