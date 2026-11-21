from __future__ import annotations

import pickle
import subprocess
import sys
from pathlib import Path

import pylosat
import pytest
from conftest import fasta
from pylosat import Alignment, ArgumentError, Pssm, Search, Sequences


def test_keywords_are_the_words_of_the_command_line():
    search = Search("blastn", task="megablast", evalue=1e-5, max_target_seqs=10)
    assert search.program == "blastn"
    assert search.argv == ("-task", "megablast", "-evalue", "1e-05", "-max_target_seqs", "10")
    assert search == Search.from_argv("blastn", search.argv)


def test_a_real_keeps_its_value():
    value = 0.1 + 0.2
    words = Search("blastp", evalue=value).argv
    assert float(words[1]) == value


def test_each_kind_of_value():
    assert Search("blastp", ungapped=True).argv == ("-ungapped",)
    assert Search("blastp", ungapped=False).argv == ()
    assert Search("blastn", soft_masking=False).argv == ("-soft_masking", "false")
    assert Search("blastn", word_size=11).argv == ("-word_size", "11")
    assert Search("blastp", evalue=1).argv == ("-evalue", "1.0")
    assert Search("blastp", matrix="PAM30").argv == ("-matrix", "PAM30")


@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("ungapped", 1),
        ("soft_masking", "true"),
        ("word_size", True),
        ("evalue", True),
        ("word_size", 1.5),
        ("evalue", "1e-5"),
        ("task", 1),
    ],
)
def test_a_value_of_another_type_is_refused(option, value):
    with pytest.raises(TypeError):
        Search("blastn", **{option: value})


@pytest.mark.parametrize("option", ["query", "subject", "db", "out", "outfmt", "no_such_option"])
def test_what_is_not_an_option_here_is_refused(option):
    with pytest.raises(TypeError, match=option):
        Search("blastn", **{option: "x"})
    if option != "no_such_option":
        with pytest.raises(ArgumentError):
            Search.from_argv("blastn", [f"-{option}", "x"])


@pytest.mark.parametrize("option", ["out_pssm", "out_ascii_pssm", "save_pssm_after_last_round", "save_each_pssm"])
def test_psiblast_takes_its_pssm_files(option):
    takes_a_file = option in ("out_pssm", "out_ascii_pssm")
    words = (f"-{option}", "p.txt") if takes_a_file else (f"-{option}",)
    assert Search("psiblast", **{option: "p.txt" if takes_a_file else True}).argv == words
    assert Search.from_argv("psiblast", list(words)).argv == words


def test_a_pssm_query_is_an_input_not_an_option():
    # Series Q item Q6b: the PSSM query is `Pssm`, which the search gives as `-in_pssm`.
    for program in ("psiblast", "tblastn"):
        with pytest.raises(TypeError, match="in_pssm"):
            Search(program, in_pssm="q.asn")
        with pytest.raises(ArgumentError, match="a PSSM query is given as the inputs"):
            Search.from_argv(program, ["-in_pssm", "q.asn"])
    with pytest.raises(TypeError, match="or an Alignment"):
        Search("psiblast").run("q.asn", Sequences.from_text(">s\nMKT\n"))
    with pytest.raises(ArgumentError, match="psiblast and tblastn, not by blastp"):
        Search("blastp").run(Pssm.from_file("q.asn"), Sequences.from_text(">s\nMKT\n"))


def test_an_alignment_query_is_an_input_and_its_master_choice_an_option():
    # Line LIB item 1.6: the alignment query is `Alignment`, which the search gives as `-in_msa`;
    # `-msa_master_idx` and `-ignore_msa_master` are options of the search.
    with pytest.raises(TypeError, match="in_msa"):
        Search("psiblast", in_msa="m.aln")
    with pytest.raises(ArgumentError, match="a query alignment is given as the inputs"):
        Search.from_argv("psiblast", ["-in_msa", "m.aln"])
    assert Search("psiblast", msa_master_idx=2).argv == ("-msa_master_idx", "2")
    assert Search("psiblast", ignore_msa_master=True).argv == ("-ignore_msa_master",)
    assert Search.from_argv("psiblast", ["-msa_master_idx", "2"]).argv == ("-msa_master_idx", "2")
    flags = {parameter.flag for parameter in pylosat.describe("psiblast").parameters}
    assert {"-msa_master_idx", "-ignore_msa_master"} <= flags and "-in_msa" not in flags
    # The other programs have no such options.
    with pytest.raises(TypeError, match="msa_master_idx"):
        Search("blastp", msa_master_idx=2)
    with pytest.raises(TypeError, match="or an Alignment"):
        Search("psiblast").run("m.aln", Sequences.from_text(">s\nMKT\n"))
    assert Alignment.from_file(Path("m.aln")) == Alignment("m.aln")


ROUNDS_QUERY = "outfmt0/e2e_many_query.faa"
ROUNDS_SUBJECT = "outfmt0/e2e_many_subject.faa"


def rounds_inputs() -> tuple[Sequences, Sequences]:
    return Sequences.from_fasta(fasta(ROUNDS_QUERY)), Sequences.from_fasta(fasta(ROUNDS_SUBJECT))


def test_psiblast_runs_its_rounds_here():
    # Line RG item Q5: psiblast's rounds. NCBI 2.17.0 finds 300 subjects in round 1 and 299 in
    # round 2, then the rounds converge.
    search = Search("psiblast", evalue=1e-5, num_iterations=3)
    outcome = search.run(*rounds_inputs())
    one = Search("psiblast", evalue=1e-5).run(*rounds_inputs())
    [rounds] = outcome.rounds
    assert rounds.converged and not one.rounds[0].converged
    assert len(rounds.stats) == 2 and rounds.stats[0] == outcome.query_stats[0]
    assert rounds.stats[0] != rounds.stats[1]
    hsps = outcome.hsps
    assert [hsp.index for hsp in hsps] == list(range(len(hsps)))
    assert {hsp.round for hsp in hsps} == {1, 2}
    first = [hsp.e_value for hsp in hsps if hsp.round == 1]
    assert first == [hsp.e_value for hsp in one.hsps]
    # The views: each round with its HSPs and statistics, a subject once a round.
    [view] = outcome.queries()
    assert view.round is None and view.stats == rounds.stats[0]
    by_round = view.rounds()
    assert [part.round for part in by_round] == [1, 2]
    assert [part.stats for part in by_round] == rounds.stats
    assert sum(len(part.hsps) for part in by_round) == len(hsps)
    assert all(hsp.round == part.round for part in by_round for hsp in part.hsps)
    assert by_round[1].rounds() == [by_round[1]]
    pairs = [(subject.s_idx, subject.round) for subject in view.subjects()]
    assert len(pairs) == len(set(pairs)) and {round for _, round in pairs} == {1, 2}
    found = view.subject(pairs[-1][0])
    assert found is not None and found.round == 1
    # The batches give the same records and rounds.
    batches = list(search.batches(*rounds_inputs()))
    assert [batch.rounds for batch in batches] == [outcome.rounds]
    assert [hsp.round for hsp in batches[0].hsps] == [hsp.round for hsp in hsps]


def test_psiblast_rounds_are_the_command_lines(command_line):
    words = ["-evalue", "1e-5", "-num_iterations", "3", "-outfmt", "7"]
    expected, _ = command_line("psiblast", fasta(ROUNDS_QUERY), fasta(ROUNDS_SUBJECT), words)
    search = Search("psiblast", evalue=1e-5, num_iterations=3)
    outcome = search.run(*rounds_inputs(), outfmt="7")
    assert outcome.reports[0].bytes == expected
    text = outcome.reports[0].bytes.decode()
    for hsp in outcome.hsps:
        span = outcome.reports[0].spans[hsp.index]
        assert span is not None
        heading = text.rindex("# Iteration: ", 0, span.start)
        assert text[heading:].splitlines()[0] == f"# Iteration: {hsp.round}"


def test_the_parser_refuses_as_on_the_command_line():
    with pytest.raises(ArgumentError):
        Search.from_argv("blastp", ["-evalue"])
    with pytest.raises(ArgumentError, match="unknown program"):
        Search("nosuchprogram")  # type: ignore[call-overload]


def test_an_error_survives_pickle():
    with pytest.raises(ArgumentError) as caught:
        Search.from_argv("blastp", ["-evalue"])
    copy = pickle.loads(pickle.dumps(caught.value))
    assert type(copy) is ArgumentError
    assert isinstance(copy, pylosat.LosatError)
    assert str(copy) == str(caught.value)


def test_a_search_is_a_value():
    search = Search("tblastx", query_gencode=11)
    assert pickle.loads(pickle.dumps(search)) == search
    assert hash(Search.from_argv("tblastx", search.argv)) == hash(search)
    assert search != Search("tblastx")
    assert "tblastx" in repr(search)


def test_describe_is_the_engines_table():
    assert pylosat.PROGRAMS == ("blastn", "blastp", "blastx", "tblastn", "tblastx", "psiblast")
    description = pylosat.describe("blastp")
    parameters = {parameter.name: parameter for parameter in description.parameters}
    assert parameters["evalue"].flag == "-evalue"
    assert parameters["evalue"].value_kind == "real"
    assert parameters["ungapped"].value_kind == "flag"
    assert not parameters["ungapped"].takes_value
    assert "blastp-fast" in parameters["task"].choices
    assert "outfmt" not in parameters
    assert description.query_gencodes == ()
    assert 11 in pylosat.describe("tblastx").subject_gencodes


def test_the_option_stubs_are_those_of_the_table():
    tool = Path(__file__).resolve().parents[1] / "tools" / "generate_option_stubs.py"
    done = subprocess.run([sys.executable, str(tool), "--check"], capture_output=True)
    assert done.returncode == 0, "run tools/generate_option_stubs.py\n" + done.stdout.decode()
