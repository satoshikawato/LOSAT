from __future__ import annotations

import pickle
import subprocess
import sys
from pathlib import Path

import pylosat
import pytest
from pylosat import ArgumentError, Search


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
def test_psiblast_writes_no_pssm_file_here(option):
    takes_a_file = option in ("out_pssm", "out_ascii_pssm")
    with pytest.raises(TypeError, match=option):
        Search("psiblast", **{option: "p.txt" if takes_a_file else True})
    words = [f"-{option}", "p.txt"] if takes_a_file else [f"-{option}"]
    with pytest.raises(ArgumentError, match="writes no PSSM file"):
        Search.from_argv("psiblast", words)


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
    assert done.returncode == 0, "run tools/generate_option_stubs.py"
