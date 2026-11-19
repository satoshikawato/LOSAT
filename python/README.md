# pylosat

Python bindings for [LOSAT](https://github.com/satoshikawato/LOSAT), a pure-Rust reimplementation of NCBI BLAST+ (`2.17.0` compatible CLI and output).

`pylosat` runs sequence searches directly inside your Python process without spawning subprocesses, exposing the exact CLI option parser, byte-identical BLAST reports, Karlin-Altschul statistics, and zero-copy [Apache Arrow](https://arrow.apache.org/) HSP tables (`polars`, `pyarrow`, `pandas`).

- **Documentation:** [satoshikawato.github.io/LOSAT/libraries/python.html](https://satoshikawato.github.io/LOSAT/libraries/python.html)
- **Repository:** [github.com/satoshikawato/LOSAT](https://github.com/satoshikawato/LOSAT)
- **Issues:** [github.com/satoshikawato/LOSAT/issues](https://github.com/satoshikawato/LOSAT/issues)

## Quick Start

```python
import polars as pl
import pylosat

search = pylosat.Search("blastn", task="megablast", evalue=1e-5, num_threads=4)
query = pylosat.Sequences.from_fasta("query.fna")
subject = pylosat.Sequences.from_fasta("subject.fna")

outcome = search.run(query, subject, outfmt=["6", "0"])

# Exact CLI report bytes / text
print(outcome.reports[0].text)

# Zero-copy Arrow HSP table into Polars (or pyarrow.table / pandas.DataFrame.from_arrow)
df = pl.DataFrame(outcome.hsps)

# Search-space statistics per query
print(outcome.query_stats[0].eff_space)
```

## Installation

Pre-built wheels for Linux (`x86_64`, `aarch64`), macOS (Apple silicon and Intel), and Windows (`x64`) on Python 3.10+ are scheduled for PyPI immediately after v0.2.0:

```bash
pip install pylosat
```

You can also download wheels built for every release architecture from the [Python wheels GitHub Actions workflow](https://github.com/satoshikawato/LOSAT/actions/workflows/python-wheels.yml), or build locally from the `python/` directory with [maturin](https://www.maturin.rs/) and Rust 1.92+:

```bash
cd python
pip install maturin
maturin build --release
pip install target/wheels/pylosat-*.whl
```

## Supported Programs and Tools

`pylosat.PROGRAMS` lists the search programs available through `pylosat.Search`:

- **Search**: `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `psiblast`
- **Databases and Masking**: `makeblastdb`, `blastdbcmd`, `dustmasker`, `segmasker`, `convert2blastmask`

For `rpsblast`, `rpstblastn`, and `rpsbproc`, use the [`losat` CLI](https://satoshikawato.github.io/LOSAT/installation.html).

## Building Searches and Reading Options

`Search` accepts any CLI option as a keyword argument without the leading hyphen (`evalue=1e-5` corresponds to `-evalue 1e-5`; boolean flags take `True`/`False`, such as `ungapped=True`):

```python
# Inspect every parameter accepted by a program directly from the engine
desc = pylosat.describe("blastp")
for param in desc.parameters[:5]:
    print(param.flag, param.value_kind, param.default)

# Or construct directly from CLI tokens (excluding I/O routing flags)
search = pylosat.Search.from_argv("blastn", ["-task", "dc-megablast", "-evalue", "1e-5"])
```

Input/output routing options (`-query`, `-subject`, `-db`, `-in_pssm`, `-out`, `-outfmt`) are passed via Python arguments to `run`, `run_to_file`, or `batches` instead of `Search`.

## Inputs

- **FASTA files and in-memory sequences**:
  - `Sequences.from_fasta(path)` (automatically decompresses `.gz`, `.bz2`, and `.zst` files)
  - `Sequences.from_text(text)`
  - `Sequences.from_pairs([(title, sequence), ...])`
- **Formatted BLAST databases (`-db`)**: `Database("core_nt")` resolves names using `BLASTDB` and `.ncbirc` exactly like the CLI:

```python
search = pylosat.Search("blastn", taxids="9606", evalue=1e-5, num_threads=4)
outcome = search.run(query, pylosat.Database("core_nt"), outfmt="6")
```

- **PSI-BLAST and PSSM Checkpoints (`-in_pssm`)**: Pass `Pssm.from_file("model.asn")` to `psiblast` or `tblastn`. For iterative `psiblast`, `outcome.rounds[q]` holds per-round `QueryStats` and `converged` (`Search has CONVERGED!`), while `outcome.queries()[q].rounds()` splits HSPs per iteration (`hsp.round`). `-out_pssm` and `-out_ascii_pssm` write checkpoint files to disk as on the CLI.

## Results and Arrow `HspTable`

`outcome.hsps` implements the Arrow PyCapsule C Data Interface (`__arrow_c_array__` / `__arrow_c_stream__`) with unrounded numeric values and zero-copy consumption by `polars.DataFrame(outcome.hsps)`, `pyarrow.table(outcome.hsps)`, or `pandas.DataFrame.from_arrow(outcome.hsps)`. No external Arrow Python dependency is required.

All 31 schema columns exposed by `HspTable`:

| Group | Columns |
| :--- | :--- |
| **Default `-outfmt 6`** | `qaccver`, `saccver`, `pident`, `length`, `mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, `bitscore` |
| **Extended alignment** | `score`, `nident`, `positive`, `gaps`, `qlen`, `slen`, `qframe`, `sframe`, `sstrand`, `qseq`, `sseq` |
| **Metadata & PSI-BLAST** | `query_strand`, `sum_n`, `comp_adjust_method`, `index`, `q_idx`, `s_idx`, `rank`, `round` |

Row access via `outcome.hsps[i]` returns an `Hsp` object with 1-based inclusive coordinates (`q_start`, `q_end`, `s_start`, `s_end`) plus 0-based half-open Python `slice` helpers (`q_range`/`query_range`, `subject_range`). All `Outcome`, `QueryHits`, `SubjectHits`, `QueryStats`, `QueryRounds`, `Report`, and `Search` objects are picklable.

To write `-outfmt 13` / `-outfmt 14` reports (with side files next to `-out`) or large outputs directly to disk:

```python
diagnostics = search.run_to_file(query, subject, outfmt="14", path="report.xml")
```

## Streaming Large Query Sets and Interruption

Iterate over NCBI query batches incrementally with `search.batches(query, subject)` to bound peak memory and track progress:

```python
batches = search.batches(query, subject)
for batch in batches:
    print(f"Processed up to query {batch.query_indices.stop} / {batches.query_count}")
```

The Python GIL is released while searches execute. Between query batches (or when leaving the loop / pressing Ctrl-C during `batches` or report-free `run`), execution stops cleanly.

## Database and Masking Helpers

Run `makeblastdb`, `blastdbcmd`, `dustmasker`, `segmasker`, and `convert2blastmask` with optional `input=` for standard input and receive a `ToolOutput` (`.bytes`, `.text`, `.diagnostics`):

```python
pylosat.makeblastdb("-dbtype", "prot", "-out", "proteins", input=fasta_text)
info = pylosat.blastdbcmd("-db", "proteins", "-info").text
masked = pylosat.dustmasker("-outfmt", "fasta", input=nucleotide_fasta).text
```

## Exceptions

- `LosatError`: Base exception for all `pylosat` errors
  - `ArgumentError`: Invalid CLI options or values
  - `InputError`: Invalid FASTA or input data
  - `SearchError`: Runtime search, database, or tool failure

## License

Released under the [MIT License](https://github.com/satoshikawato/LOSAT/blob/main/LICENSE).
