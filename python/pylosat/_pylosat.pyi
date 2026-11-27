"""The native module. `pylosat` is its public face; only `Hsp`, `HspTable` and the exceptions show."""

from collections.abc import Iterator
from typing import Any, overload

ENGINE_VERSION: str

class LosatError(Exception):
    """A failure that LOSAT reports. `str(error)` is the command line's message for the same
    failure as text, with one U+FFFD per maximal ill-formed subpart of its bytes, as
    `bytes.decode("utf-8", "replace")` gives; `message_bytes` holds the bytes. A search's
    message of NCBI's keeps its newlines; LOSAT's own refusals have no `Error: ` and no last
    newline; an input's error, a tool's and that of a start that ends the program are the
    command line's text without its last newlines. Every exception that pylosat raises has
    `message_bytes`."""

    message_bytes: bytes

class ArgumentError(LosatError): ...
class InputError(LosatError): ...
class SearchError(LosatError): ...

class Hsp:
    """One HSP: a row of its table, with the field names of LOSAT's `HspRecord`."""

    def __new__(cls, table: HspTable, row: int) -> Hsp: ...
    @property
    def index(self) -> int:
        """0-based position in the final, sorted hit list: the key of `Report.spans`."""
    @property
    def q_idx(self) -> int:
        """0-based index of the query record."""
    @property
    def s_idx(self) -> int:
        """0-based index of the subject record; of a database, the sequence's ordinal id."""
    @property
    def rank(self) -> int:
        """0-based position among the HSPs of its query in its round."""
    @property
    def round(self) -> int:
        """psiblast's round that found the HSP (1-based); 1 in a search without rounds."""
    @property
    def query_id(self) -> str:
        """`qaccver`, as text: its bytes with one U+FFFD per maximal ill-formed subpart, as
        `bytes.decode("utf-8", "replace")` gives; `report.bytes[report.spans[hsp.index]]` holds
        the row with the exact bytes."""
    @property
    def subject_id(self) -> str:
        """`saccver`, as text, as `query_id`."""
    @property
    def raw_score(self) -> int:
        """`score`."""
    @property
    def bit_score(self) -> float:
        """`bitscore`, not rounded."""
    @property
    def e_value(self) -> float:
        """`evalue`, not rounded. A value of the whole search (`QueryStats.eff_space`)."""
    @property
    def percent_identity(self) -> float:
        """`pident`, not rounded."""
    @property
    def q_start(self) -> int:
        """`qstart`, 1-based."""
    @property
    def q_end(self) -> int:
        """`qend`, 1-based: before `q_start` on the minus strand."""
    @property
    def s_start(self) -> int:
        """`sstart`, 1-based."""
    @property
    def s_end(self) -> int:
        """`send`, 1-based: before `s_start` on the minus strand."""
    @property
    def query_strand(self) -> str | None:
        """`plus` or `minus`; `None` for a protein query."""
    @property
    def subject_strand(self) -> str | None:
        """`sstrand`: `plus` or `minus`; `None` for a protein subject."""
    @property
    def query_frame(self) -> int | None:
        """`qframe` where the program translates the query, `None` otherwise."""
    @property
    def subject_frame(self) -> int | None:
        """`sframe` where the program translates the subject, `None` otherwise."""
    @property
    def query_length(self) -> int:
        """`qlen`."""
    @property
    def subject_length(self) -> int:
        """`slen`."""
    @property
    def length(self) -> int:
        """`length`: the alignment's columns, gaps included."""
    @property
    def num_ident(self) -> int:
        """`nident`."""
    @property
    def num_positives(self) -> int:
        """`positive`."""
    @property
    def mismatch(self) -> int:
        """`mismatch`."""
    @property
    def gap_open(self) -> int:
        """`gapopen`."""
    @property
    def gaps(self) -> int:
        """`gaps`."""
    @property
    def sum_n(self) -> int | None:
        """How many HSPs NCBI's sum statistics linked with this one, where the program links."""
    @property
    def comp_adjust_method(self) -> int | None:
        """The composition-based statistics method, where the program has one."""
    @property
    def query_aligned(self) -> str | None:
        """`qseq`: the aligned query row, `-` for gaps."""
    @property
    def subject_aligned(self) -> str | None:
        """`sseq`: the aligned subject row, `-` for gaps."""
    @property
    def query_range(self) -> slice:
        """The query interval as a 0-based slice on the record's plus strand."""
    @property
    def subject_range(self) -> slice:
        """The subject interval as a 0-based slice on the record's plus strand."""

class HspTable:
    """The HSPs of a search: one row per HSP, in the order of the final hit list.

    Any Arrow consumer reads a table without a copy: `pyarrow.table(hsps)`,
    `polars.DataFrame(hsps)`, `pandas.DataFrame.from_arrow(hsps)`. The columns have the names of
    `-outfmt`, the default report's first.

    The ids (`qaccver`, `saccver`) are text, with one U+FFFD per maximal ill-formed subpart of
    their bytes, as `bytes.decode("utf-8", "replace")` gives; the exact bytes are those of the
    HSP's row in a report, `report.bytes[report.spans[hsp.index]]`.
    """

    def __new__(cls, ipc_stream: bytes) -> HspTable: ...
    @property
    def columns(self) -> list[str]: ...
    def __len__(self) -> int: ...
    @overload
    def __getitem__(self, key: int) -> Hsp: ...
    @overload
    def __getitem__(self, key: slice) -> HspTable: ...
    def __iter__(self) -> Iterator[Hsp]: ...
    def __arrow_c_stream__(self, requested_schema: object | None = None) -> object: ...
    def _query(self, q_idx: int) -> HspTable: ...
    def _subjects(self) -> list[tuple[int, int, HspTable]]: ...
    def _round(self, round: int) -> HspTable: ...

class Batches:
    def __new__(
        cls,
        program: str,
        argv: list[str],
        query: tuple[bytes | str, str] | str | dict[str, str],
        subject: tuple[bytes | str | None, str],
    ) -> Batches: ...
    @property
    def query_count(self) -> int: ...
    def __iter__(self) -> Batches: ...
    def __next__(
        self,
    ) -> tuple[
        int,
        int,
        HspTable,
        list[tuple[int, int, float, float, float] | None],
        list[tuple[list[tuple[int, int, float, float, float] | None], bool]],
        bytes,
    ]: ...

def describe(
    program: str,
) -> tuple[list[tuple[str, str, bool, str, str | None, list[str]]], list[int], list[int]]: ...
def validate(program: str, argv: list[str]) -> None: ...
def dry_run(program: str, argv: list[str]) -> bool: ...
def run(
    program: str,
    argv: list[str],
    query: tuple[bytes | str, str] | str | dict[str, str],
    subject: tuple[bytes | str | None, str],
    formats: list[str],
) -> tuple[
    list[tuple[str, bytes, list[slice | None], list[slice | None]]],
    HspTable,
    list[tuple[int, int, float, float, float] | None],
    list[tuple[list[tuple[int, int, float, float, float] | None], bool]],
    bytes,
]: ...
def run_to_file(
    program: str,
    argv: list[str],
    query: tuple[bytes | str, str] | str | dict[str, str],
    subject: tuple[bytes | str | None, str],
    outfmt: str,
    path: str,
) -> bytes: ...
def makeblastdb(argv: list[str], input: bytes) -> tuple[bytes, bytes]: ...
def blastdbcmd(argv: list[str], input: bytes) -> tuple[bytes, bytes]: ...
def dustmasker(argv: list[str], input: bytes) -> tuple[bytes, bytes]: ...
def segmasker(argv: list[str], input: bytes) -> tuple[bytes, bytes]: ...
def windowmasker(argv: list[str], input: bytes) -> tuple[bytes, bytes]: ...
def convert2blastmask(argv: list[str], input: bytes) -> tuple[bytes, bytes]: ...
def makembindex(argv: list[str], input: bytes) -> tuple[bytes, bytes]: ...
def blast_formatter(argv: list[str], input: bytes) -> tuple[bytes, bytes]: ...
