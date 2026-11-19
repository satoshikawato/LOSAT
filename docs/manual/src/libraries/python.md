# Python (`pylosat`)

`pylosat` embeds the LOSAT search engine directly inside Python 3.10+ via PyO3. Searches use the exact CLI option parser, run with the GIL released, yield zero-copy [Apache Arrow](https://arrow.apache.org/) HSP tables (`polars`, `pyarrow`, `pandas`), and return exact BLAST report bytes with zero subprocess overhead.

## Installation

Pre-built binary wheels for Linux (`x86_64`, `aarch64`), macOS (Apple silicon and Intel), and Windows (`x64`) are scheduled on PyPI immediately after v0.2.0:

```bash
pip install pylosat
```

Before PyPI publication (or to build from source), install from the `python/` directory of the repository using [maturin](https://www.maturin.rs/) (requires Rust 1.92+), or download wheels from the [Python wheels workflow](https://github.com/satoshikawato/LOSAT/actions/workflows/python-wheels.yml):

```bash
git clone https://github.com/satoshikawato/LOSAT.git
cd LOSAT/python
pip install maturin
maturin build --release
pip install target/wheels/pylosat-*.whl
```

## Quick start

```python
import polars as pl
import pylosat

search = pylosat.Search("blastn", task="megablast", evalue=1e-5, num_threads=4)
query = pylosat.Sequences.from_fasta("query.fna")
subject = pylosat.Sequences.from_fasta("subject.fna")

outcome = search.run(query, subject, outfmt=["6", "0"])

# Exact bytes of `losat blastn ... -outfmt 6`
print(outcome.reports[0].text)

# Zero-copy Arrow table into Polars, PyArrow, or Pandas
df = pl.DataFrame(outcome.hsps)

# Per-query Karlin-Altschul statistics
print(outcome.query_stats[0].eff_space)
```

## Supported programs

`pylosat.PROGRAMS` exposes:

- **Searches:** `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `psiblast`
- **Database and masking helpers:** `makeblastdb`, `blastdbcmd`, `dustmasker`, `segmasker`, `convert2blastmask`

For `rpsblast`, `rpstblastn`, and `rpsbproc`, use the [command-line binary](../installation.md).

## Options and introspection

`Search` accepts CLI flags as keyword arguments without the leading dash (`evalue=1e-5` corresponds to `-evalue 1e-5`, boolean flags take `True`/`False`). Input and output routing options (`-query`, `-subject`, `-db`, `-in_pssm`, `-out`, `-outfmt`) are passed directly to `run()`, `run_to_file()`, or `batches()`:

```python
# Inspect options programmatically directly from the engine
desc = pylosat.describe("blastp")
for param in desc.parameters[:5]:
    print(param.flag, param.value_kind, param.default)

# Or build directly from CLI tokens
search = pylosat.Search.from_argv("blastn", ["-task", "dc-megablast", "-evalue", "1e-5"])
```

Every argument is checked by LOSAT's CLI parser (`ArgumentError` on invalid options or values).

## Inputs: FASTA, BLAST databases, and PSSMs

- **`Sequences.from_fasta(path)`:** Reads FASTA from disk and automatically decompresses `.gz`, `.bz2`, and `.zst` inputs. Pass `name=...` to control how the input is named inside reports.
- **`Sequences.from_text(text)`** and **`Sequences.from_pairs([(title, seq), ...])`:** Construct FASTA inputs in memory without touching disk.
- **`Database(name)`:** Opens a local BLAST database (`-db`), resolved using the current working directory, `BLASTDB`, and `.ncbirc` (see [BLAST databases](../databases.md)). Filtering options such as `taxids="9606"` or `seqidlist=...` are passed on `Search`.

```python
search = pylosat.Search("blastn", taxids="9606", evalue=1e-5, num_threads=4)
outcome = search.run(query, pylosat.Database("core_nt"), outfmt="6")
```

- **`Pssm.from_file(path)`:** Loads an ASN.1 checkpoint (`-in_pssm`) for `psiblast` or `tblastn`.

## Iterative PSI-BLAST and results

`outcome.hsps` implements the Arrow PyCapsule C Data Interface (`__arrow_c_array__` and `__arrow_c_stream__`), zero-copy readable via `polars.DataFrame(hsps)`, `pyarrow.table(hsps)`, or `pandas.DataFrame.from_arrow(hsps)` without requiring Arrow libraries at runtime.

| Group | Arrow `HspTable` columns |
| :--- | :--- |
| **Default `-outfmt 6`** | `qaccver`, `saccver`, `pident`, `length`, `mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, `bitscore` |
| **Extended fields** | `score`, `nident`, `positive`, `gaps`, `qlen`, `slen`, `qframe`, `sframe`, `sstrand`, `qseq`, `sseq` |
| **Metadata & PSI-BLAST** | `query_strand`, `sum_n`, `comp_adjust_method`, `index`, `q_idx`, `s_idx`, `rank`, `round` |

All numeric values in `HspTable` and `Hsp` are full engine precision. Indexing `outcome.hsps[i]` yields an `Hsp` object with 1-based inclusive coordinates (`q_start`, `q_end`, `s_start`, `s_end`) plus 0-based `slice` helpers (`query_range`, `subject_range`). `qframe` and `sframe` are `None` when not applicable.

In iterative `psiblast` searches:
- `hsp.round` identifies which iteration found each hit (1 for non-iterative searches).
- `outcome.rounds[q]` gives per-round `QueryStats` and `converged` (`Search has CONVERGED!`).
- `outcome.queries()[q].rounds()` splits HSPs per round with each round's statistics.
- `-out_pssm` and `-out_ascii_pssm` write checkpoint files to disk identically to the CLI.

All result objects (`Outcome`, `QueryHits`, `SubjectHits`, `QueryStats`, `QueryRounds`, `Report`, `Search`) are picklable for use with `multiprocessing`.

To write large reports or `-outfmt 13` / `14` (which write companion files next to `-out`) straight to disk:

```python
diagnostics = search.run_to_file(query, subject, outfmt="14", path="report.xml")
```

## Streaming batches and cancellation

Use `search.batches(query, subject)` to iterate one NCBI query batch at a time and bound peak memory:

```python
batches = search.batches(query, subject)
for batch in batches:
    print(f"{batch.query_indices.stop} / {batches.query_count}")
```

The GIL is released during search execution. Leaving the `batches` loop or pressing Ctrl-C cleanly cancels execution between query batches (when calling `run` with formats requested, Ctrl-C takes effect after the engine call returns).

## Database and masking helpers

```python
pylosat.makeblastdb("-dbtype", "prot", "-out", "proteins", input=fasta_text)
info = pylosat.blastdbcmd("-db", "proteins", "-info").text
masked = pylosat.dustmasker("-outfmt", "fasta", input=nucleotide_fasta).text
```

Available helpers: `makeblastdb`, `blastdbcmd`, `dustmasker`, `segmasker`, `convert2blastmask`. All accept optional `input=` (`bytes` or `str`) as standard input and return a `ToolOutput` with `.bytes`, `.text`, and `.diagnostics`.

## Exceptions

- `LosatError`: base exception
  - `ArgumentError`: invalid options or values
  - `InputError`: FASTA parsing failure
  - `SearchError`: search, database, or tool error
