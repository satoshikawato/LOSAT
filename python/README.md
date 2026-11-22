# pylosat

The Python package of [LOSAT](https://github.com/satoshikawato/LOSAT), a pure-Rust reimplementation of NCBI BLAST+ 2.17.0. A search runs inside your Python process, with the command line's options, and returns the command line's reports byte for byte. The HSPs come as a table that pandas, Polars and PyArrow read without a copy. No executable is started. It supports CPython 3.10 and later.

- **Documentation:** [satoshikawato.github.io/LOSAT/libraries/python.html](https://satoshikawato.github.io/LOSAT/libraries/python.html)
- **Repository:** [github.com/satoshikawato/LOSAT](https://github.com/satoshikawato/LOSAT)
- **Issues:** [github.com/satoshikawato/LOSAT/issues](https://github.com/satoshikawato/LOSAT/issues)

```python
import polars
import pylosat

search = pylosat.Search("blastn", task="megablast", evalue=1e-5, num_threads=4)
query = pylosat.Sequences.from_fasta("query.fna")
subject = pylosat.Sequences.from_fasta("subject.fna")

outcome = search.run(query, subject, outfmt=["6", "0"])
print(outcome.reports[0].text)          # the -outfmt 6 report, as `losat blastn ... -outfmt 6` writes it
outcome.reports[1].bytes                # the bytes of -outfmt 0
hits = polars.DataFrame(outcome.hsps)   # one row per HSP, the numbers not rounded
outcome.query_stats[0].eff_space        # what fixed the E-values of the first query
outcome.diagnostics                     # what the command line writes to its standard error (bytes: diagnostics_bytes)
```

## Installation

```bash
pip install pylosat
```

PyPI has wheels for CPython 3.10 and later on Linux (x86_64, aarch64; glibc 2.17 or later), macOS (arm64, x86_64) and Windows (x64). On other platforms pip builds pylosat from its source distribution, which needs Rust 1.92 or later ([rustup](https://rustup.rs/)).

To build a wheel from the repository instead, use [maturin](https://www.maturin.rs/):

```bash
git clone https://github.com/satoshikawato/LOSAT.git
cd LOSAT/python
pip install maturin
maturin build --release
pip install target/wheels/pylosat-*.whl
```

For development, `maturin develop` builds the extension into the active virtual environment.

## Programs

`pylosat.PROGRAMS` names the search programs: `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx` and `psiblast`. The tools `makeblastdb`, `blastdbcmd`, `dustmasker`, `segmasker`, `windowmasker`, `convert2blastmask` and `makembindex` are functions of their own (below). `rpsblast`, `rpstblastn` and `rpsbproc` run only on the command line.

## Options

`Search` takes the program and the command line's options, each as a keyword without its dash: `evalue=1e-5` is `-evalue 1e-5`, and `ungapped=True` is `-ungapped`. `pylosat.describe("blastn")` lists them from the engine's own table, and the type stubs list them by program. `Search.from_argv("blastn", ["-evalue", "1e-5"])` takes the words themselves. An unknown keyword, or a value of the wrong Python type, raises `TypeError`; the command line's parser checks the rest, so an option means, and fails, as it does there (`ArgumentError`).

`-query`, `-subject`, `-db`, `-in_pssm`, `-in_msa`, `-out` and `-outfmt` are not options here: the inputs and the formats are arguments of `run`. psiblast runs its rounds (`num_iterations`) as on the command line, and its `-out_pssm` and `-out_ascii_pssm` (with `-save_pssm_after_last_round` and `-save_each_pssm`) write the command line's PSSM files under their names.

## Inputs

A `Sequences` is one FASTA input, and the search reads it with the command line's reader:

- `Sequences.from_fasta(path)`: a file, inflated when it is named `.gz`, `.bz2` or `.zst`.
- `Sequences.from_text(text)`: FASTA text.
- `Sequences.from_pairs([(title, sequence), ...])`: the FASTA file of those records.

A PSSM takes the place of the query of psiblast and tblastn: `Pssm.from_file(path)` is `-in_pssm path` (a PSSM file as `-out_pssm` writes it). Its one query is the PSSM's sequence: `q_idx` 0.

An alignment takes the place of the query of psiblast: `Alignment.from_file(path)` is `-in_msa path` (a multiple sequence alignment; the search reads it as the command line does). Its one query is the alignment's master: `q_idx` 0. `msa_master_idx` and `ignore_msa_master` are options of `Search`, as on the command line, and choose the master.

A BLAST database takes the place of the subject: `Database("core_nt")` is `-db core_nt`, found as the command line finds it (`BLASTDB`, `.ncbirc`). The options of a database search (`taxids`, `seqidlist`, ...) go to `Search` with the others. Then `hsp.s_idx` is the sequence's ordinal id, the same in every search of that database.

```python
search = pylosat.Search("blastn", taxids="9606", evalue=1e-5)
outcome = search.run(query, pylosat.Database("core_nt"), outfmt="6")
```

## Results

`outcome.hsps` is an `HspTable` of 31 columns. Any Arrow consumer reads it without a copy, through the Arrow PyCapsule interface: `polars.DataFrame(hsps)`, `pyarrow.table(hsps)`, `pandas.DataFrame.from_arrow(hsps)`. pylosat itself needs no Arrow library. The columns have the names of `-outfmt`, the default report's first (`qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore`), then `score nident positive gaps qlen slen qframe sframe sstrand qseq sseq`, then pylosat's own `query_strand sum_n comp_adjust_method index q_idx s_idx rank round`.

The numbers are the engine's, not rounded. Text to show comes from a report: `report.spans[hsp.index]` is the slice of `report.bytes` that holds the row or section of an HSP. `qframe` and `sframe` are null where the program does not translate that side; the reports print a number there. The ids (`qaccver`, `saccver`, `hsp.query_id`, `hsp.subject_id`) are text, with one U+FFFD per maximal ill-formed subpart of their bytes, as `bytes.decode("utf-8", "replace")` gives; the report's row keeps the bytes (`report.bytes[report.spans[hsp.index]]`).

`hsps[i]` is an `Hsp` with the field names of LOSAT's Rust API (`bit_score`, `e_value`, `q_start`). Its coordinates are 1-based, as in the reports; `query_range` and `subject_range` are 0-based slices on the record's plus strand. `outcome.query_stats[q]` holds the statistics of query `q` (`eff_space`, `hsp_len`, `lambda_`, `kappa`, `entropy`). `outcome.queries()` and `QueryHits.subjects()` group the rows by query and by subject. Everything that a search returns can be pickled, so a process pool can return it.

psiblast's rounds of each query follow one another in the table, round 1 first: `hsp.round` (the `round` column) is the round that found the HSP, 1 in the other searches. `outcome.rounds[q]` holds the statistics of each round of query `q` and whether its rounds converged (`Search has CONVERGED!`), and `QueryHits.rounds()` gives each round's HSPs with its statistics.

## Files

`search.run_to_file(query, subject, outfmt="14", path="report.xml")` is the command line's `-out`, and the way to outfmt 13 and 14: their side files lie next to `path`. It returns a `FileOutcome`: what the command line writes to its standard error, as `diagnostics_bytes` and as the text `diagnostics`.

## Long searches

```python
batches = search.batches(query, subject)
for batch in batches:
    print(f"{batch.query_indices.stop} of {batches.query_count} queries")
    hits = polars.DataFrame(batch.hsps)
```

`batches` searches one of NCBI's query batches at a time, so a run holds the results of one batch. Their rows and statistics, one batch after the other, are those of the whole search. Leaving the loop stops the search.

The GIL is released while a search runs, so other Python threads go on. Python cannot stop a thread inside its process, so Ctrl-C takes effect between two query batches, in `batches` and in a `run` without reports. A `run` with reports is one call of the engine, and Ctrl-C takes effect when it returns.

## Databases and masks

```python
pylosat.makeblastdb("-dbtype", "prot", "-title", "proteins", "-out", "proteins", input=fasta_text)
print(pylosat.blastdbcmd("-db", "proteins", "-info").text)
print(pylosat.dustmasker("-outfmt", "fasta", input=nucleotide_fasta).text)
print(pylosat.segmasker("-outfmt", "fasta", input=protein_fasta).text)
print(pylosat.convert2blastmask("-masking_algorithm", "other", "-masking_options", "mine",
                                input=lowercase_masked_fasta).text)
```

`makeblastdb`, `blastdbcmd`, `dustmasker`, `segmasker`, `windowmasker`, `convert2blastmask` and `makembindex` take the command line's words and return a `ToolOutput`: what the tool writes to standard output (`bytes`, `text`) and to standard error (`diagnostics_bytes`, `diagnostics`). `input` is what the command line reads as standard input, so a database can be built from FASTA in memory (with `-title`, which makeblastdb needs for standard input), ids given to `blastdbcmd -entry_batch -`, sequences in memory masked, and the lower case letters of FASTA in memory turned into masking information.

## Errors

`ArgumentError` (the options), `InputError` (a FASTA input that does not open, or a subject record that the reader refuses) and `SearchError` (the search, a query record that the reader refuses, in its query batch as on the command line, or a tool) derive from `LosatError`, each with the command line's message for the same failure: `message_bytes` holds its bytes, and `str(error)` is the text. A search's message of NCBI's keeps its newlines; LOSAT's own refusals have no `Error: ` and no last newline; an input's error, a tool's and that of a start that ends the program are the command line's text without its last newlines; text that a search writes before the failure, such as a warning, is not in the message.

Every text that pylosat gives of what the command line writes (`str(error)`, `outcome.diagnostics`, `batch.diagnostics`, `ToolOutput.diagnostics` and `.text`, `FileOutcome.diagnostics` from `run_to_file`) has one U+FFFD per maximal ill-formed subpart of its bytes, as `bytes.decode("utf-8", "replace")` gives; the bytes are `message_bytes`, `diagnostics_bytes` and `ToolOutput.bytes`.

## Development

See [`CONTRIBUTING.md`](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/CONTRIBUTING.md). The boundary of this package is set by [`PD-LOSAT-PYTHON-BINDING-BOUNDARY`](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-PYTHON-BINDING-BOUNDARY.md).

## License

LOSAT is a Rust translation of the search programs and database tools of NCBI BLAST+ 2.17.0. The code written by NCBI is in the public domain in the United States. LOSAT's own work (the changes and additions of the translation, threading, the WebAssembly builds, the tests, the Web application and the Python package) is under the [MIT License](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/LICENSE). Parts that come from other projects (LMDB, Boost, Arm Optimized Routines, Sun's fdlibm, Microsoft's C runtime, BitMagic, the Rust crates bzip2 and ruzstd) keep their own notices: see [THIRD_PARTY_NOTICES](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/THIRD_PARTY_NOTICES).
