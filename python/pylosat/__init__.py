"""LOSAT, a pure-Rust reimplementation of NCBI BLAST+, as a Python library.

>>> import pylosat
>>> search = pylosat.Search("blastn", task="megablast", evalue=1e-5)
>>> query = pylosat.Sequences.from_fasta("query.fna")
>>> subject = pylosat.Sequences.from_fasta("subject.fna")
>>> outcome = search.run(query, subject, outfmt="6")
>>> outcome.reports[0].text      # the bytes of `losat blastn ... -outfmt 6`
>>> polars.DataFrame(outcome.hsps)
"""

from importlib.metadata import PackageNotFoundError, version

from ._api import (
    PROGRAMS,
    Alignment,
    Batch,
    Batches,
    Database,
    Description,
    FileOutcome,
    Outcome,
    Parameter,
    Pssm,
    QueryHits,
    QueryRounds,
    QueryStats,
    Report,
    Search,
    Sequences,
    SubjectHits,
    ToolOutput,
    blast_formatter,
    blastdbcmd,
    convert2blastmask,
    describe,
    dustmasker,
    makeblastdb,
    makembindex,
    segmasker,
    windowmasker,
)
from ._pylosat import (
    ENGINE_VERSION,
    ArgumentError,
    Hsp,
    HspTable,
    InputError,
    LosatError,
    SearchError,
)

try:
    __version__ = version("pylosat")
except PackageNotFoundError:
    __version__ = "unknown"

__all__ = [
    "ENGINE_VERSION",
    "PROGRAMS",
    "Alignment",
    "ArgumentError",
    "Batch",
    "Batches",
    "Database",
    "Description",
    "FileOutcome",
    "Hsp",
    "HspTable",
    "InputError",
    "LosatError",
    "Outcome",
    "Parameter",
    "Pssm",
    "QueryHits",
    "QueryRounds",
    "QueryStats",
    "Report",
    "Search",
    "SearchError",
    "Sequences",
    "SubjectHits",
    "ToolOutput",
    "__version__",
    "blast_formatter",
    "blastdbcmd",
    "convert2blastmask",
    "describe",
    "dustmasker",
    "makeblastdb",
    "makembindex",
    "segmasker",
    "windowmasker",
]
