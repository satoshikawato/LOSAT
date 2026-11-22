# Python (`pylosat`)

`pylosat` runs LOSAT searches inside your Python process. You give it the command line's options and FASTA input; you get back the command line's reports, byte for byte, and the HSPs as a table that pandas, Polars and PyArrow read without a copy. No executable is started. It supports CPython 3.10 and later.

## Installation

```bash
pip install pylosat
```

[PyPI](https://pypi.org/project/pylosat/) has wheels for CPython 3.10 and later on Linux (x86_64, aarch64; glibc 2.17 or later), macOS (arm64, x86_64) and Windows (x64). On other platforms pip builds `pylosat` from its source distribution, which needs Rust 1.92 or later ([rustup](https://rustup.rs/)).

To build a wheel from the repository instead, use [maturin](https://www.maturin.rs/):

```bash
git clone https://github.com/satoshikawato/LOSAT.git
cd LOSAT/python
pip install maturin
maturin build --release
pip install target/wheels/pylosat-*.whl
```

For development, `maturin develop` builds the extension into the active virtual environment.

## Quick start

```python
import polars
import pylosat

search = pylosat.Search("blastn", task="megablast", evalue=1e-5, num_threads=4)
query = pylosat.Sequences.from_fasta("query.fna")
subject = pylosat.Sequences.from_fasta("subject.fna")

outcome = search.run(query, subject, outfmt=["6", "0"])
print(outcome.reports[0].text)          # the -outfmt 6 report, as the command line writes it
outcome.reports[1].bytes                # the bytes of -outfmt 0
hits = polars.DataFrame(outcome.hsps)   # one row per HSP, the numbers not rounded
outcome.query_stats[0].eff_space        # what fixed the E-values of the first query
outcome.diagnostics                     # what the command line writes to its standard error
```

## What it runs

`pylosat.PROGRAMS` names the six search programs: `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx` and `psiblast`. The tools `makeblastdb`, `blastdbcmd`, `dustmasker`, `segmasker`, `windowmasker`, `convert2blastmask` and `makembindex` are functions of their own ([below](#database-and-masking-tools)). `rpsblast`, `rpstblastn` and `rpsbproc` run only on the [command line](../installation.md).

## Options

`Search(program, **options)` takes the command line's options as keywords without the dash: `evalue=1e-5` is `-evalue 1e-5`, and `ungapped=True` is `-ungapped`. The package's type stubs list the options of each program, so an editor can complete them.

```python
desc = pylosat.describe("blastp")          # the options, from the engine's own table
for param in desc.parameters[:5]:
    print(param.flag, param.value_kind, param.default)

# or the words of a command line, after the program name
search = pylosat.Search.from_argv("blastn", ["-task", "dc-megablast", "-evalue", "1e-5"])
```

An unknown keyword, or a value of the wrong Python type, raises `TypeError`. The command line's parser checks everything else, so an option means, and fails, as it does there (`ArgumentError`, with the parser's message).

`-query`, `-subject`, `-db`, `-in_pssm`, `-in_msa`, `-out` and `-outfmt` are not options here: the inputs, the formats and the output file are arguments of `run`, `batches` and `run_to_file`. `psiblast` runs its rounds (`num_iterations`) as on the command line, and its PSSM files (`out_pssm` and `out_ascii_pssm`, with `save_pssm_after_last_round` and `save_each_pssm`) are written where the command line writes them. `psiblast` starts from a Clustal alignment as `-in_msa` does: give an `Alignment` as the query, and choose its master with `msa_master_idx` or `ignore_msa_master`, which are options of `Search` like the others (the command line takes them with `-in_msa` only: a search of another query fails with the parser's message).

## Inputs

- **`Sequences.from_fasta(path)`** reads a FASTA file, and inflates it when it is named `.gz`, `.bz2` or `.zst`, as the command line does. The reports show `name=...` where the command line shows the file name; by default, the path.
- **`Sequences.from_text(text)`** and **`Sequences.from_pairs([(title, sequence), ...])`** make FASTA input in memory.
- **`Database(name)`** is `-db name`. It is found as the command line finds it (`BLASTDB`, `.ncbirc`; see [BLAST databases](../databases.md)). The options of a database search (`taxids`, `seqidlist`, ...) go to `Search` with the others. In a database search, `hsp.s_idx` is the sequence's ordinal id, the same in every search of that database.
- **`Pssm.from_file(path)`** is `-in_pssm path` of `psiblast` and `tblastn`: a PSSM file as `psiblast -out_pssm` writes it. Its one query is the PSSM's sequence (query 0).
- **`Alignment.from_file(path)`** is `-in_msa path` of `psiblast`: a multiple sequence alignment, which the search reads as the command line reads it (a file it cannot read fails with the command line's message). Its one query is the alignment's master (query 0), the first row or the one `msa_master_idx` names; the PSSM is built from the alignment as NCBI BLAST+ builds it.

```python
search = pylosat.Search("blastn", taxids="9606", evalue=1e-5, num_threads=4)
outcome = search.run(query, pylosat.Database("core_nt"), outfmt="6")
```

## Results

`search.run(...)` returns an `Outcome`:

- **`outcome.reports`**: one `Report` per requested format, in the order you asked. `report.bytes` and `report.text` are the report; `report.spans[hsp.index]` is the slice of `report.bytes` that holds the row or section of an HSP.
- **`outcome.hsps`**: an `HspTable`, one row per HSP in the order of the final hit list. Any Arrow consumer reads it without a copy, through the Arrow PyCapsule interface: `polars.DataFrame(hsps)`, `pyarrow.table(hsps)`, `pandas.DataFrame.from_arrow(hsps)`. `pylosat` itself needs no Arrow library.
- **`outcome.query_stats[q]`**: the statistics of query `q`: `eff_space` (effective search space), `hsp_len` (length adjustment), `lambda_`, `kappa` and `entropy` (Karlin-Altschul lambda, K and H).
- **`outcome.diagnostics`**: what the command line writes to its standard error, as text with one U+FFFD per maximal ill-formed subpart of its bytes; `outcome.diagnostics_bytes` is the bytes.

The table has 31 columns. The first twelve are those of `-outfmt 6`:

| Columns | Meaning |
|:---|:---|
| `qaccver`, `saccver`, `pident`, `length`, `mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, `bitscore` | The `-outfmt 6` columns |
| `score`, `nident`, `positive`, `gaps`, `qlen`, `slen`, `qframe`, `sframe`, `sstrand`, `qseq`, `sseq` | More `-outfmt` columns |
| `query_strand`, `sum_n`, `comp_adjust_method`, `index`, `q_idx`, `s_idx`, `rank`, `round` | `pylosat`'s own: the query's strand, the number of HSPs linked by sum statistics, the composition-based statistics method, the HSP's place in the hit list, the query and subject numbers, the HSP's place among the HSPs of its query, all 0-based, and the `psiblast` round |

The numbers are the engine's and are not rounded. For text to show, read the report. `qframe` and `sframe` are empty (null) where the program does not translate that side; the reports print a number there. The ids (`qaccver`, `saccver`, `hsp.query_id`, `hsp.subject_id`) are text, with one U+FFFD per maximal ill-formed subpart of their bytes, as `bytes.decode("utf-8", "replace")` gives; the report's row keeps the bytes (`report.bytes[report.spans[hsp.index]]`).

`outcome.hsps[i]` is an `Hsp`, with the field names of the Rust library (`bit_score`, `e_value`, `q_start`, ...). Its coordinates are 1-based, as in the reports; `hsp.query_range` and `hsp.subject_range` are 0-based Python slices on the plus strand of the record. `outcome.queries()` and `QueryHits.subjects()` group the rows by query and by subject.

In `psiblast`, the rounds of each query follow one another in the table, round 1 first. `hsp.round` is the round that found the HSP (1 in the other programs). `outcome.rounds[q]` holds the statistics of each round of query `q` and whether its rounds converged (`Search has CONVERGED!`), and `outcome.queries()[q].rounds()` gives each round's HSPs with its statistics.

Everything a search returns can be pickled, so a process pool can return it.

## Writing to a file

`run_to_file` is the command line's `-out`. Use it for large reports and for `-outfmt 13` and `14`, whose side files are written next to the file. It returns a `FileOutcome`: what the command line writes to its standard error, as `diagnostics_bytes` and as the text `diagnostics` (U+FFFD for a byte that is not UTF-8).

```python
written = search.run_to_file(query, subject, outfmt="14", path="report.xml")
print(written.diagnostics)
```

## Long searches

`search.batches(query, subject)` searches one of NCBI's query batches at a time, so your program holds the results of one batch. The rows and statistics of the batches, one after the other, are those of the whole search. Leaving the loop stops the search.

```python
batches = search.batches(query, subject)
for batch in batches:
    print(f"{batch.query_indices.stop} of {batches.query_count} queries")
    hits = polars.DataFrame(batch.hsps)
```

The GIL is released while a search runs, so other Python threads go on. Python cannot stop a thread inside its process, so Ctrl-C takes effect between two query batches, in `batches` and in a `run` without reports. A `run` with reports is one call of the engine, and Ctrl-C takes effect when it returns.

## Database and masking tools

`makeblastdb`, `blastdbcmd`, `dustmasker`, `segmasker`, `windowmasker`, `convert2blastmask` and `makembindex` take the command line's words and return a `ToolOutput`: what the tool writes to standard output (`.bytes`, `.text`) and to standard error (`.diagnostics_bytes`, `.diagnostics`). Each text has one U+FFFD per maximal ill-formed subpart of its bytes. `input=` (text or bytes) is what the tool reads as standard input, so you can build a database from FASTA in memory, give ids to `blastdbcmd -entry_batch -`, or mask sequences in memory.

```python
pylosat.makeblastdb("-dbtype", "prot", "-title", "proteins", "-out", "proteins", input=fasta_text)
info = pylosat.blastdbcmd("-db", "proteins", "-info").text
masked = pylosat.dustmasker("-outfmt", "fasta", input=nucleotide_fasta).text
```

As on the command line, `makeblastdb` needs `-title` when it reads standard input.

## Errors

All errors derive from `LosatError` and carry the command line's message for the same failure:

- `ArgumentError`: an option or value that the parser refuses.
- `InputError`: a FASTA input that does not open, or a subject record that the reader refuses (the command line reads every subject before the first query).
- `SearchError`: a search, database or tool that fails. A query record that the reader refuses is one too: the command line reads its queries one batch at a time, so the search ends in the batch of that record, after the batches before it. `batches` gives those batches, and `run_to_file` writes their reports, as the command line's `-out` does.

`str(error)` is the message as text, with one U+FFFD per maximal ill-formed subpart of its bytes; `error.message_bytes` is the bytes that the command line writes.

## What `pylosat` promises, and what it does not

- It never changes your sequences, options or statistics, and never splits one search into several searches whose results are joined.
- Reports are the command line's bytes. `pylosat` does not rewrite, filter or regenerate them.
- The Arrow table, pickling and the Python types are `pylosat` formats. They are not NCBI output.
- A search stops only between query batches.
- Sequences, headers, file names and results are never sent over the network.
- The extension calls no NCBI program and starts no process.

The boundary is set by [the decision record on the Python binding boundary](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-PYTHON-BINDING-BOUNDARY.md).
