"""The Python side of the binding: the shapes of `losat::api`, in Python's terms.

Nothing here computes a BLAST value or formats one
(docs/product_decisions/PD-LOSAT-PYTHON-BINDING-BOUNDARY.md): options become the command line's
words, and results are what the engine returned.
"""

from __future__ import annotations

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
class SubjectHits:
    """The HSPs of one query against one subject."""

    q_idx: int
    s_idx: int
    hsps: HspTable


@dataclass(frozen=True)
class QueryHits:
    """The results of one query."""

    q_idx: int
    stats: QueryStats | None
    hsps: HspTable

    def subjects(self) -> list[SubjectHits]:
        """The subjects with HSPs for this query, in the order of the hit list."""
        return [
            SubjectHits(self.q_idx, s_idx, hsps) for s_idx, hsps in self.hsps._subjects()
        ]

    def subject(self, s_idx: int) -> SubjectHits | None:
        """The HSPs of this query against subject `s_idx`, if it has any."""
        return next((hits for hits in self.subjects() if hits.s_idx == s_idx), None)


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
    """By query record; `None` for a query that has no statistics."""
    diagnostics: str
    """What the command line writes to its standard error for this search."""

    def queries(self) -> list[QueryHits]:
        """The results of each query record, also of those without hits."""
        return [
            QueryHits(q_idx, stats, self.hsps._query(q_idx))
            for q_idx, stats in enumerate(self.query_stats)
        ]

    def query(self, q_idx: int) -> QueryHits:
        """The results of query record `q_idx`."""
        return QueryHits(q_idx, self.query_stats[q_idx], self.hsps._query(q_idx))


@dataclass(frozen=True)
class Batch:
    """The results of one of NCBI's query batches."""

    query_indices: range
    """The query records of the batch."""
    hsps: HspTable
    """The batch's HSPs. `q_idx` and `index` count as in the whole search."""
    query_stats: list[QueryStats | None] = field(repr=False)
    """One entry per query of `query_indices`."""
    diagnostics: str
    """What the command line writes to its standard error for this batch."""

    def queries(self) -> list[QueryHits]:
        """The results of each query of the batch."""
        return [
            QueryHits(q_idx, stats, self.hsps._query(q_idx))
            for q_idx, stats in zip(self.query_indices, self.query_stats)
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
        start, end, hsps, stats, diagnostics = next(self._running)
        return Batch(range(start, end), hsps, _stats(stats), diagnostics)


class Search:
    """A program and its options, checked by the command line's parser.

    Each keyword is an option of the command line without its dash: `evalue=1e-5` is
    `-evalue 1e-5`. `describe` lists them. `-query`, `-subject`, `-db`, `-out` and `-outfmt` are
    not options here, nor are psiblast's PSSM files (the library writes none). The search itself
    is `run`, `batches` or `run_to_file`. The NCBI C++ Toolkit's `-logfile`, `-conffile` and
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
        query: Sequences,
        subject: Sequences | Database,
        *,
        outfmt: str | Sequence[str] = (),
    ) -> Outcome:
        """Runs the search; `outfmt` names the reports, as the values of `-outfmt`.

        The GIL is released. Without reports, Ctrl-C stops the search between two query batches.
        With reports the engine runs in one call, and Ctrl-C takes effect when it returns.
        """
        formats = [outfmt] if isinstance(outfmt, str) else list(outfmt)
        reports, hsps, stats, diagnostics = _pylosat.run(
            self._program,
            list(self._argv),
            _input(query, "query"),
            _subject(subject),
            formats,
        )
        return Outcome([Report(*report) for report in reports], hsps, _stats(stats), diagnostics)

    def run_to_file(
        self,
        query: Sequences,
        subject: Sequences | Database,
        *,
        outfmt: str,
        path: str | os.PathLike[str],
    ) -> str:
        """Runs the search and writes the report of `outfmt` to the file `path`, as `-out path`.

        Outfmt 13 and 14 come through here: their side files lie next to `path`. Returns what the
        command line writes to its standard error. The GIL is released.
        """
        return _pylosat.run_to_file(
            self._program,
            list(self._argv),
            _input(query, "query"),
            _subject(subject),
            outfmt,
            os.fspath(path),
        )

    def batches(self, query: Sequences, subject: Sequences | Database) -> Batches:
        """The search as an iterator of `Batch`, one of NCBI's query batches at a time."""
        return Batches(
            _pylosat.Batches(
                self._program,
                list(self._argv),
                _input(query, "query"),
                _subject(subject),
            )
        )


@dataclass(frozen=True)
class ToolOutput:
    """What a database program writes: its standard output and its standard error."""

    bytes: bytes
    diagnostics: str

    @property
    def text(self) -> str:
        return self.bytes.decode(errors="replace")


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


def convert2blastmask(*words: str, input: bytes | str = b"") -> ToolOutput:
    """Runs `convert2blastmask` with the command line's words, which turns the lower case
    letters of FASTA input into the masking information that `makeblastdb -mask_data` reads.

    `input` is what the command line reads as standard input (`-in -`, and without `-in`), so
    sequences in memory can be converted. A failure raises `SearchError` with what the command
    line writes to its standard error. The GIL is released.
    """
    data = input.encode() if isinstance(input, str) else bytes(input)
    return ToolOutput(*_pylosat.convert2blastmask(list(words), data))
