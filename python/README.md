# pylosat

The Python binding of [LOSAT](https://github.com/satoshikawato/LOSAT), a pure-Rust reimplementation
of NCBI BLAST+. A search runs in the Python process, with the command line's options and the
command line's output bytes, and no executable is started.

```python
import polars
import pylosat

search = pylosat.Search("blastn", task="megablast", evalue=1e-5, num_threads=4)
query = pylosat.Sequences.from_fasta("query.fna")
subject = pylosat.Sequences.from_fasta("subject.fna")

outcome = search.run(query, subject, outfmt=["6", "0"])
outcome.reports[0].text           # the bytes of `losat blastn ... -outfmt 6`, as text
outcome.reports[1].bytes          # the bytes of `-outfmt 0`
polars.DataFrame(outcome.hsps)    # one row per HSP, the numbers not rounded
outcome.query_stats[0].eff_space  # what fixed the E-values of the first query
outcome.diagnostics               # what the command line writes to its standard error
```

## Options

`Search` takes the program and the command line's options, each as a keyword without its dash:
`evalue=1e-5` is `-evalue 1e-5`, and `ungapped=True` is `-ungapped`. `pylosat.describe("blastn")`
lists them, from the engine's own table, and the type stubs list them by program.
`Search.from_argv("blastn", ["-evalue", "1e-5"])` takes the words themselves. The command line's
parser checks every option, so an option means, and fails, as it does there (`ArgumentError`).

`-query`, `-subject`, `-db`, `-out` and `-outfmt` are not options here: the inputs and the formats
are arguments of `run`. A database subject (`-db`) is not supported yet.

## Inputs

A `Sequences` is one FASTA input, and the search reads it with the command line's reader:

- `Sequences.from_fasta(path)`: a file, inflated when it is named `.gz`, `.bz2` or `.zst`.
- `Sequences.from_text(text)`: FASTA text.
- `Sequences.from_pairs([(title, sequence), ...])`: the FASTA file of those records.

## Results

`outcome.hsps` is an `HspTable`. Any Arrow consumer reads it without a copy:
`polars.DataFrame(hsps)`, `pyarrow.table(hsps)`, `pandas.DataFrame.from_arrow(hsps)`. No Arrow
library is required to use pylosat. Its columns have the names of `-outfmt`, the default report's
first: `qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore`.

`hsps[i]` is an `Hsp` with the field names of LOSAT's Rust API (`bit_score`, `e_value`, `q_start`).
`outcome.queries()` and `QueryHits.subjects()` group the rows by query and by subject. Everything
that a search returns can be pickled, so a process pool can return it.

The numbers are the engine's, not rounded. Text to show comes from a report: `report.spans[hsp.index]`
is the slice of `report.bytes` that holds the row or section of an HSP.

`qframe` and `sframe` are null where the program does not translate that side; the reports print
a number there.

## Long searches

```python
batches = search.batches(query, subject)
for batch in batches:
    print(f"{batch.query_indices.stop} of {batches.query_count} queries")
    save(polars.DataFrame(batch.hsps))
```

`batches` searches one of NCBI's query batches at a time, so a run holds the results of one batch.
Their rows and statistics, one batch after the other, are those of the whole search. Leaving the
loop stops the search.

The GIL is released while a search runs, so other Python threads go on. Python cannot stop a thread
inside its process, so Ctrl-C takes effect between two query batches, in `batches` and in a `run`
without reports. A `run` with reports is one call of the engine, and Ctrl-C takes effect when it
returns. One query longer than a batch is one batch.

## Errors

`ArgumentError` (the options), `InputError` (a FASTA input) and `SearchError` (the search) derive
from `LosatError`, each with the command line's message for the same failure.

## Development

See [`AGENTS.md`](AGENTS.md). The boundary of this package is
[`PD-LOSAT-PYTHON-BINDING-BOUNDARY`](../docs/product_decisions/PD-LOSAT-PYTHON-BINDING-BOUNDARY.md).
