"""The Python side of the binding: the shapes of `losat::api`, in Python's terms.

Nothing here computes a BLAST value or formats one
(docs/product_decisions/PD-LOSAT-PYTHON-BINDING-BOUNDARY.md): options become the command line's
words, and results are what the engine returned.
"""

from __future__ import annotations

import builtins
import functools
import numbers
import operator
import os
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal, overload

from . import _pylosat
from ._pylosat import HspTable

if TYPE_CHECKING:
    from typing_extensions import Self, Unpack

    from ._options import (
        BlastnOptions,
        BlastpOptions,
        BlastxOptions,
        PsiblastOptions,
        TblastnOptions,
        TblastxOptions,
    )

PROGRAMS = ("blastn", "blastp", "blastx", "tblastn", "tblastx", "psiblast")
"""The search programs, by their command-line names."""


@dataclass(frozen=True)
class Parameter:
    """One option of a program, as the command line defines it."""

    name: str
    """The keyword argument of `Search` (`evalue`)."""
    flag: str
    """The command line's flag (`-evalue`)."""
    help: str
    takes_value: bool
    value_kind: str
    """`flag`, `boolean`, `integer`, `real`, `path` or `text`."""
    default: str | None
    """The default of the parser; `None` for an option without one, whatever its help says."""
    choices: tuple[str, ...]
    """The values the parser accepts, when it lists them."""


@dataclass(frozen=True)
class Description:
    """The options of a program and the genetic codes its parser accepts."""

    program: str
    parameters: tuple[Parameter, ...]
    query_gencodes: tuple[int, ...]
    subject_gencodes: tuple[int, ...]


@functools.cache
def describe(program: str) -> Description:
    """The options of `program`, read from the engine's own table of arguments."""
    rows, query_gencodes, subject_gencodes = _pylosat.describe(program)
    parameters = tuple(
        Parameter(flag.lstrip("-"), flag, help_text, takes_value, kind, default, tuple(choices))
        for flag, help_text, takes_value, kind, default, choices in rows
    )
    return Description(program, parameters, tuple(query_gencodes), tuple(subject_gencodes))


def _words(parameter: Parameter, value: object) -> list[str]:
    """The command line's words for one option."""
    kind = parameter.value_kind
    if kind in ("flag", "boolean"):
        if not isinstance(value, bool):
            raise TypeError(f"{parameter.name} takes a bool, not {value!r}")
        if kind == "boolean":
            return [parameter.flag, "true" if value else "false"]
        return [parameter.flag] if value else []
    if kind in ("integer", "real"):
        if isinstance(value, bool) or not isinstance(value, numbers.Real):
            raise TypeError(f"{parameter.name} takes a number, not {value!r}")
        if kind == "integer":
            return [parameter.flag, str(operator.index(value))]  # type: ignore[arg-type]
        return [parameter.flag, repr(float(value))]
    if kind == "path":
        return [parameter.flag, os.fspath(value)]  # type: ignore[call-overload]
    if not isinstance(value, str):
        raise TypeError(f"{parameter.name} takes a str, not {value!r}")
    return [parameter.flag, value]


@dataclass(frozen=True)
class Sequences:
    """A FASTA input of a search.

    The search reads it with the command line's reader, by the rules of its program. Use one of
    the constructors: a path and a sequence are never the same argument.
    """

    _source: bytes | str = field(repr=False)
    """The FASTA bytes, or the name of a file."""
    name: str | None = None
    """What the reports show as the file's name (`-query`, `-subject`)."""

    @classmethod
    def from_fasta(cls, path: str | os.PathLike[str], *, name: str | None = None) -> Self:
        """The FASTA file `path`, read when a search starts.

        A file named `.gz`, `.bz2` or `.zst` is inflated, as on the command line. The reports show
        `name`, by default the path.
        """
        path = os.fspath(path)
        return cls(path, path if name is None else name)

    @classmethod
    def from_text(cls, text: str | bytes, *, name: str | None = None) -> Self:
        """The FASTA text itself."""
        return cls(text.encode() if isinstance(text, str) else bytes(text), name)

    @classmethod
    def from_pairs(cls, pairs: Iterable[tuple[str, str]], *, name: str | None = None) -> Self:
        """Records of `(title, sequence)`, as the FASTA file `>title` and `sequence` lines."""
        lines = []
        for title, sequence in pairs:
            if any(line_end in title + sequence for line_end in "\r\n"):
                raise ValueError(f"a line end in the record {title!r}")
            lines.append(f">{title}\n{sequence}\n")
        return cls("".join(lines).encode(), name)

    def _named(self, default_name: str) -> tuple[bytes | str, str]:
        return self._source, default_name if self.name is None else self.name


@dataclass(frozen=True)
class Database:
    """A BLAST database as the subject of a search: the value of `-db`.

    The search opens it by its name, one or more space-separated names, and finds it as the
    command line finds it (`BLASTDB`, `.ncbirc`).
    """

    name: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", os.fspath(self.name))


@dataclass(frozen=True)
class Pssm:
    """A PSSM file as the query of a search: the value of `-in_pssm` of psiblast and tblastn.

    The search reads the file as the command line reads it (psiblast's checkpoint file, as
    `out_pssm` writes it); its one query is the PSSM's sequence, query 0 of the results.
    """

    path: str

    @classmethod
    def from_file(cls, path: str | os.PathLike[str]) -> Self:
        """The PSSM file `path`, read when a search starts."""
        return cls(os.fspath(path))

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", os.fspath(self.path))


@dataclass(frozen=True)
class Alignment:
    """An alignment file as the query of a search: the value of `-in_msa` of psiblast.

    The search reads the file as the command line reads it and makes the query's PSSM from it, as
    NCBI BLAST+ does; its one query is the alignment's master, query 0 of the results. The master
    is the first row, or the row that the search's `msa_master_idx` names; `ignore_msa_master`
    builds the PSSM without counting it. The two options need an `Alignment` as the query.
    """

    path: str

    @classmethod
    def from_file(cls, path: str | os.PathLike[str]) -> Self:
        """The alignment file `path`, read when a search starts."""
        return cls(os.fspath(path))

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", os.fspath(self.path))


def _query(query: Sequences | Pssm | Alignment) -> tuple[bytes | str, str] | str | dict[str, str]:
    if isinstance(query, Pssm):
        return query.path
    if isinstance(query, Alignment):
        return {"alignment": query.path}
    if not isinstance(query, Sequences):
        raise TypeError(
            "the query is a Sequences (from_fasta, from_text or from_pairs), "
            f"a Pssm (from_file) or an Alignment (from_file), not {query!r}"
        )
    return _input(query, "query")


def _input(sequences: Sequences, role: str) -> tuple[bytes | str, str]:
    if not isinstance(sequences, Sequences):
        raise TypeError(
            f"the {role} is a Sequences (from_fasta, from_text or from_pairs), "
            f"not {sequences!r}"
        )
    return sequences._named(role)


def _subject(subject: Sequences | Database) -> tuple[bytes | str | None, str]:
    if isinstance(subject, Database):
        return None, subject.name
    return _input(subject, "subject")


@dataclass(frozen=True)
class QueryStats:
    """The statistics of one query (NCBI's `CBlastAncillaryData`): what fixed its E-values."""

    eff_space: int
    """The effective search space, in protein letters where a side is translated."""
    hsp_len: int
    """The length adjustment."""
    lambda_: float
    """Karlin-Altschul lambda."""
    kappa: float
    """Karlin-Altschul K."""
    entropy: float
    """Karlin-Altschul H."""


def _stats(rows: Iterable[tuple[int, int, float, float, float] | None]) -> list[QueryStats | None]:
    return [None if row is None else QueryStats(*row) for row in rows]


@dataclass(frozen=True)
class QueryRounds:
    """The rounds of one query: psiblast's rounds (`num_iterations` other than 1), one round for the
    other searches."""

    stats: list[QueryStats | None]
    """The statistics of each round, round 1 first: what fixed the E-values of the HSPs of that round
    (`Hsp.round`). The first is the query's `query_stats` entry."""
    converged: bool
    """Whether the rounds converged (NCBI's `Search has CONVERGED!`); false without rounds."""


_RoundsRow = tuple[list[tuple[int, int, float, float, float] | None], bool]


def _rounds(rows: Iterable[_RoundsRow]) -> list[QueryRounds]:
    return [QueryRounds(_stats(stats), converged) for stats, converged in rows]


@dataclass(frozen=True)
class SubjectHits:
    """The HSPs of one query against one subject in one round."""

    q_idx: int
    s_idx: int
    hsps: HspTable
    round: int = 1
    """The round of the HSPs (`Hsp.round`): 1 in a search without rounds."""


@dataclass(frozen=True)
class QueryHits:
    """The results of one query: of all its rounds, or of one (`rounds`)."""

    q_idx: int
    stats: QueryStats | None
    """The query's statistics (with psiblast's rounds round 1's); in the view of one round, that
    round's."""
    hsps: HspTable
    round: int | None = None
    """`None` for the HSPs of all the query's rounds, the round in the view of one round."""
    _round_stats: list[QueryStats | None] = field(default_factory=list, repr=False)
    """The statistics of the rounds of the view, from its first round (`rounds`)."""

    def subjects(self) -> list[SubjectHits]:
        """The subjects with HSPs for this query, in the order of the hit list, round after round
        (a subject found again in a later round comes again there)."""
        return [
            SubjectHits(self.q_idx, s_idx, hsps, round)
            for s_idx, round, hsps in self.hsps._subjects()
        ]

    def subject(self, s_idx: int) -> SubjectHits | None:
        """The HSPs of this query against subject `s_idx` in the first round where it has any."""
        return next((hits for hits in self.subjects() if hits.s_idx == s_idx), None)

    def rounds(self) -> list[QueryHits]:
        """The results of each round of this query, round 1 first, each with its HSPs and its
        statistics: one for a search without rounds. The view of one round gives itself."""
        first = 1 if self.round is None else self.round
        return [
            QueryHits(
                self.q_idx, stats, self.hsps._round(first + offset), first + offset, [stats]
            )
            for offset, stats in enumerate(self._round_stats)
        ]


def _query_hits(
    hsps: HspTable, q_idx: int, stats: QueryStats | None, rounds: QueryRounds | None
) -> QueryHits:
    return QueryHits(
        q_idx, stats, hsps._query(q_idx), None, [] if rounds is None else rounds.stats
    )


@dataclass(frozen=True)
class Report:
    """One report: the bytes that the command line writes to `-out` for this `-outfmt`."""

    outfmt: str
    bytes: bytes = field(repr=False)
    spans: list[slice | None] = field(repr=False)
    """By `Hsp.index`: the HSP's row or section in `bytes`, if this report shows it."""
    subject_headings: list[slice | None] = field(repr=False)
    """By `Hsp.index`: the heading of the HSP's subject in `bytes`, in the pairwise format."""

    @property
    def text(self) -> str:
        """The report as text. Raises `UnicodeDecodeError` for a title that is not UTF-8."""
        return self.bytes.decode()


@dataclass(frozen=True)
class Outcome:
    """What a search returned."""

    reports: list[Report]
    """One report per requested format, in the requested order."""
    hsps: HspTable
    """The final hit list, one row per HSP."""
    query_stats: list[QueryStats | None] = field(repr=False)
    """By query record; `None` for a query that has no statistics. With psiblast's rounds, round
    1's (`rounds`)."""
    rounds: list[QueryRounds] = field(repr=False)
    """By query record: the statistics of each of its rounds and whether they converged."""
    diagnostics_bytes: bytes
    """What the command line writes to its standard error for this search: its bytes."""

    @property
    def diagnostics(self) -> str:
        """`diagnostics_bytes` as text, with one U+FFFD per maximal ill-formed subpart (a byte
        that is not UTF-8, from a title read from a Latin-1 file, for example)."""
        return self.diagnostics_bytes.decode("utf-8", "replace")

    def queries(self) -> list[QueryHits]:
        """The results of each query record, also of those without hits."""
        return [
            _query_hits(self.hsps, q_idx, stats, rounds)
            for q_idx, (stats, rounds) in enumerate(zip(self.query_stats, self.rounds))
        ]

    def query(self, q_idx: int) -> QueryHits:
        """The results of query record `q_idx`."""
        return _query_hits(self.hsps, q_idx, self.query_stats[q_idx], self.rounds[q_idx])


@dataclass(frozen=True)
class FileOutcome:
    """What `Search.run_to_file` returned; the report is in its file."""

    diagnostics_bytes: bytes
    """What the command line writes to its standard error for this search: its bytes."""

    @property
    def diagnostics(self) -> str:
        """`diagnostics_bytes` as text, with one U+FFFD per maximal ill-formed subpart (a byte
        that is not UTF-8, from a title read from a Latin-1 file, for example)."""
        return self.diagnostics_bytes.decode("utf-8", "replace")


@dataclass(frozen=True)
class Batch:
    """The results of one of NCBI's query batches."""

    query_indices: range
    """The query records of the batch."""
    hsps: HspTable
    """The batch's HSPs. `q_idx` and `index` count as in the whole search."""
    query_stats: list[QueryStats | None] = field(repr=False)
    """One entry per query of `query_indices`."""
    rounds: list[QueryRounds] = field(repr=False)
    """One entry per query of `query_indices`, as in `Outcome.rounds`."""
    diagnostics_bytes: bytes
    """What the command line writes to its standard error for this batch: its bytes."""

    @property
    def diagnostics(self) -> str:
        """`diagnostics_bytes` as text, with one U+FFFD per maximal ill-formed subpart."""
        return self.diagnostics_bytes.decode("utf-8", "replace")

    def queries(self) -> list[QueryHits]:
        """The results of each query of the batch."""
        return [
            _query_hits(self.hsps, q_idx, stats, rounds)
            for q_idx, stats, rounds in zip(self.query_indices, self.query_stats, self.rounds)
        ]


class Batches(Iterator[Batch]):
    """A search of one query batch at a time; each step searches one batch.

    Leaving the loop stops the search, and Ctrl-C does so between two batches. The GIL is released
    while a batch runs.
    """

    def __init__(self, running: _pylosat.Batches) -> None:
        self._running = running

    @property
    def query_count(self) -> int:
        """How many query records the inputs hold: with `Batch.query_indices`, the progress."""
        return self._running.query_count

    def __iter__(self) -> Self:
        return self

    def __next__(self) -> Batch:
        start, end, hsps, stats, rounds, diagnostics = next(self._running)
        return Batch(range(start, end), hsps, _stats(stats), _rounds(rounds), diagnostics)


class Search:
    """A program and its options, checked by the command line's parser.

    Each keyword is an option of the command line without its dash: `evalue=1e-5` is
    `-evalue 1e-5`. `describe` lists them. `-query`, `-subject`, `-db`, `-in_pssm`, `-in_msa`, `-out` and
    `-outfmt` are not options here: the inputs are `Sequences`, `Database`, `Pssm` (the PSSM query of
    psiblast and tblastn) and `Alignment` (the alignment query of psiblast, whose master
    `msa_master_idx` and `ignore_msa_master` choose); psiblast's PSSM files (`out_pssm`, `out_ascii_pssm`) are written where the
    command line writes them. The search itself
    is `run`, `batches` or `run_to_file`. `-logfile`, `-conffile` and
    `-dryrun` work as on the command line (`from_argv`): each search opens the file of `-logfile`
    and loads the file of `-conffile`, and `dry_run` tells a search that searches nothing.
    """

    __slots__ = ("_argv", "_program")

    _argv: tuple[str, ...]
    _program: str

    @overload
    def __init__(self, program: Literal["blastn"], /, **options: Unpack[BlastnOptions]) -> None: ...
    @overload
    def __init__(self, program: Literal["blastp"], /, **options: Unpack[BlastpOptions]) -> None: ...
    @overload
    def __init__(self, program: Literal["blastx"], /, **options: Unpack[BlastxOptions]) -> None: ...
    @overload
    def __init__(
        self, program: Literal["tblastn"], /, **options: Unpack[TblastnOptions]
    ) -> None: ...
    @overload
    def __init__(
        self, program: Literal["tblastx"], /, **options: Unpack[TblastxOptions]
    ) -> None: ...
    @overload
    def __init__(
        self, program: Literal["psiblast"], /, **options: Unpack[PsiblastOptions]
    ) -> None: ...

    def __init__(self, program: str, /, **options: object) -> None:
        parameters = {parameter.name: parameter for parameter in describe(program).parameters}
        argv: list[str] = []
        for name, value in options.items():
            if name not in parameters:
                raise TypeError(f"{program} has no option {name}")
            argv.extend(_words(parameters[name], value))
        _pylosat.validate(program, argv)
        self._program = program
        self._argv = tuple(argv)

    @classmethod
    def from_argv(cls, program: str, argv: Iterable[str]) -> Self:
        """The search of the words after the program name on a command line."""
        argv = tuple(argv)
        _pylosat.validate(program, list(argv))
        search = cls.__new__(cls)
        search._program = program
        search._argv = argv
        return search

    @property
    def program(self) -> str:
        return self._program

    @property
    def argv(self) -> tuple[str, ...]:
        """The command line's words after the program name."""
        return self._argv

    @property
    def dry_run(self) -> bool:
        """Whether the words have `-dryrun`: a search then runs NCBI's application start (the file
        of `-logfile`, the file of `-conffile`) and checks the settings, and searches nothing.
        `run` gives an empty report per format, `run_to_file` makes no file and `batches` one
        batch of no query; the query and subject are not read."""
        return _pylosat.dry_run(self._program, list(self._argv))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Search):
            return NotImplemented
        return (self._program, self._argv) == (other._program, other._argv)

    def __hash__(self) -> int:
        return hash((self._program, self._argv))

    def __repr__(self) -> str:
        return f"Search.from_argv({self._program!r}, {list(self._argv)!r})"

    def __reduce__(self) -> tuple[object, tuple[str, tuple[str, ...]]]:
        return (Search.from_argv, (self._program, self._argv))

    def run(
        self,
        query: Sequences | Pssm | Alignment,
        subject: Sequences | Database,
        *,
        outfmt: str | Sequence[str] = (),
    ) -> Outcome:
        """Runs the search; `outfmt` names the reports, as the values of `-outfmt`.

        The GIL is released. Without reports, Ctrl-C stops the search between two query batches.
        With reports the engine runs in one call, and Ctrl-C takes effect when it returns.
        """
        formats = [outfmt] if isinstance(outfmt, str) else list(outfmt)
        reports, hsps, stats, rounds, diagnostics = _pylosat.run(
            self._program,
            list(self._argv),
            _query(query),
            _subject(subject),
            formats,
        )
        return Outcome(
            [Report(*report) for report in reports], hsps, _stats(stats), _rounds(rounds), diagnostics
        )

    def run_to_file(
        self,
        query: Sequences | Pssm | Alignment,
        subject: Sequences | Database,
        *,
        outfmt: str,
        path: str | os.PathLike[str],
    ) -> FileOutcome:
        """Runs the search and writes the report of `outfmt` to the file `path`, as `-out path`.

        Outfmt 13 and 14 come through here: their side files lie next to `path`. Returns what the
        command line writes to its standard error (`FileOutcome`). The GIL is released.
        """
        return FileOutcome(
            _pylosat.run_to_file(
                self._program,
                list(self._argv),
                _query(query),
                _subject(subject),
                outfmt,
                os.fspath(path),
            )
        )

    def batches(self, query: Sequences | Pssm | Alignment, subject: Sequences | Database) -> Batches:
        """The search as an iterator of `Batch`, one of NCBI's query batches at a time."""
        return Batches(
            _pylosat.Batches(
                self._program,
                list(self._argv),
                _query(query),
                _subject(subject),
            )
        )


@dataclass(frozen=True)
class ToolOutput:
    """What a database program writes: its standard output and its standard error."""

    bytes: bytes
    """The standard output's bytes."""
    # `builtins.bytes`: in the class body, the field above has the type's name.
    diagnostics_bytes: builtins.bytes
    """The standard error's bytes."""

    @property
    def text(self) -> str:
        """`bytes` as text, with one U+FFFD per maximal ill-formed subpart."""
        return self.bytes.decode("utf-8", "replace")

    @property
    def diagnostics(self) -> str:
        """`diagnostics_bytes` as text, with one U+FFFD per maximal ill-formed subpart."""
        return self.diagnostics_bytes.decode("utf-8", "replace")


def makeblastdb(*words: str, input: bytes | str = b"") -> ToolOutput:
    """Runs `makeblastdb` with the command line's words, which writes the database's files.

    `input` is what the command line reads as standard input (`-in -`, and without `-in`), so a
    database can be built from FASTA in memory. A failure raises `SearchError` with what the
    command line writes to its standard error. The GIL is released.
    """
    data = input.encode() if isinstance(input, str) else bytes(input)
    return ToolOutput(*_pylosat.makeblastdb(list(words), data))


def blastdbcmd(*words: str, input: bytes | str = b"") -> ToolOutput:
    """Runs `blastdbcmd` with the command line's words.

    `input` is what the command line reads as standard input (`-entry_batch -`, `-ipg_batch -`,
    `-taxidlist -`), so a batch of ids can be given from memory. A failure raises `SearchError`
    with what the command line writes to its standard error. The GIL is released.
    """
    data = input.encode() if isinstance(input, str) else bytes(input)
    return ToolOutput(*_pylosat.blastdbcmd(list(words), data))


def dustmasker(*words: str, input: bytes | str = b"") -> ToolOutput:
    """Runs `dustmasker` with the command line's words, which masks the low-complexity regions
    of nucleotide sequences.

    `input` is what the command line reads as standard input (`-in -`, and without `-in`), so
    sequences in memory can be masked. A failure raises `SearchError` with what the command line
    writes to its standard error. The GIL is released.
    """
    data = input.encode() if isinstance(input, str) else bytes(input)
    return ToolOutput(*_pylosat.dustmasker(list(words), data))


def segmasker(*words: str, input: bytes | str = b"") -> ToolOutput:
    """Runs `segmasker` with the command line's words, which masks the low-complexity regions
    of protein sequences.

    `input` is what the command line reads as standard input (`-in -`, and without `-in`), so
    sequences in memory can be masked. A failure raises `SearchError` with what the command line
    writes to its standard error. The GIL is released.
    """
    data = input.encode() if isinstance(input, str) else bytes(input)
    return ToolOutput(*_pylosat.segmasker(list(words), data))


def windowmasker(*words: str, input: bytes | str = b"") -> ToolOutput:
    """Runs `windowmasker` with the command line's words, which masks the repeats of nucleotide
    sequences from a file of unit counts (`-ustat`).

    `input` is what the command line reads as standard input (`-in -`, and without `-in`), so
    sequences in memory can be masked. A failure raises `SearchError` with what the command line
    writes to its standard error. The GIL is released.
    """
    data = input.encode() if isinstance(input, str) else bytes(input)
    return ToolOutput(*_pylosat.windowmasker(list(words), data))


def convert2blastmask(*words: str, input: bytes | str = b"") -> ToolOutput:
    """Runs `convert2blastmask` with the command line's words, which turns the lower case
    letters of FASTA input into the masking information that `makeblastdb -mask_data` reads.

    `input` is what the command line reads as standard input (`-in -`, and without `-in`), so
    sequences in memory can be converted. A failure raises `SearchError` with what the command
    line writes to its standard error. The GIL is released.
    """
    data = input.encode() if isinstance(input, str) else bytes(input)
    return ToolOutput(*_pylosat.convert2blastmask(list(words), data))


def makembindex(*words: str, input: bytes | str = b"") -> ToolOutput:
    """Runs `makembindex` with the command line's words, which writes the megablast index of a
    BLAST database (`-input`) or, with `-old_style_index true`, of FASTA input.

    `input` is what the command line reads as standard input (FASTA input without `-input`).
    The index files are written where the command line writes them; `diagnostics` holds its
    `creating ...` lines. A failure raises `SearchError` with what the command line writes to its
    standard error. The GIL is released.
    """
    data = input.encode() if isinstance(input, str) else bytes(input)
    return ToolOutput(*_pylosat.makembindex(list(words), data))


def blast_formatter(*words: str, input: bytes | str = b"") -> ToolOutput:
    """Runs `blast_formatter` with the command line's words, which writes the reports of a
    search again from its BLAST archive (`-outfmt 11`) without searching again.

    `input` is what the command line reads as standard input (`-archive -`), so an archive in
    memory can be formatted. A failure raises `SearchError` with what the command line writes to
    its standard error. The GIL is released.
    """
    data = input.encode() if isinstance(input, str) else bytes(input)
    return ToolOutput(*_pylosat.blast_formatter(list(words), data))
