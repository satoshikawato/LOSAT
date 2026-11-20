# Using LOSAT as a library

You can run LOSAT searches from your own program instead of the command line, in Rust (the `losat` crate) or in Python (`pylosat`). A search then runs inside your process, with the command line's options, and returns the command line's reports byte for byte. No `losat` executable is started.

| Interface | Status in v0.2.0 |
|:---|:---|
| Rust crate `losat`, module `losat::api` | Not yet supported as a download: not on crates.io. Use it from the repository, by path or git. |
| Python package `pylosat` | Not yet supported as a download: not on PyPI. Build it from the repository. |

Both live in the repository. Neither is one of the release files.

## Rust

The library is the crate `losat`. It is not on crates.io yet. Until then, depend on the repository's `LOSAT/` directory by path or by git. It needs Rust 1.92 or later.

The public surface is the module `losat::api`. Every other module is the command line or the engine; they are hidden from the documentation and outside the compatibility promise. `cargo doc --no-deps --open` in the `LOSAT/` directory shows the API.

### What `api` offers

- **Searches in memory.** `read_fasta` reads FASTA bytes into records with the command line's reader. `validate` takes a program and the command line's options (without `-query`, `-subject`, `-db`, `-in_pssm`, `-out` and `-outfmt`) and checks them with the command line's parser. `run_collect` runs the search and returns the reports, plus one typed record per HSP and the statistics of each query. The programs are `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx` and `psiblast`.
- **Inputs.** Records in memory (`Inputs::new`), a BLAST database by name (`Inputs::database`, the command line's `-db`), or the PSSM of `psiblast` and `tblastn` by file (`Inputs::pssm` and `Inputs::pssm_database`, the command line's `-in_pssm`).
- **Batches.** `run_batches` gives the same records one of NCBI's query batches at a time. A caller keeps one batch in memory, reports progress between batches and stops by not asking for the next.
- **Files.** `run_to_file` is the command line's `-out`, and the way to `-outfmt 13` and `14` with their side files. `run` writes each report to writers you give it.
- **Tools.** `run_makeblastdb`, `run_blastdbcmd`, `run_dustmasker`, `run_segmasker` and `run_convert2blastmask` run the tools with your reader as their standard input and your writers as their output.
- **Errors and output.** Errors carry the command line's messages. The library writes nothing to standard output or standard error; what the command line writes to its standard error comes back as diagnostics.

The reports are the command line's. The HSP records give the numbers of the `-outfmt` columns (scores, coordinates, strands, frames, counts), and the statistics give each query's effective search space, length adjustment, lambda, K and H.

### An example

```rust
use losat::api::{read_fasta, run_collect, validate, Inputs, Program, Role};

let query = read_fasta(Program::Blastp, Role::Query, &std::fs::read("query.faa")?)?;
let subject = read_fasta(Program::Blastp, Role::Subject, &std::fs::read("subject.faa")?)?;
let search = validate(Program::Blastp, &["-evalue", "1e-5", "-num_threads", "4"])?;
let outcome = run_collect(&search, Inputs::new(&query, &subject), &["6", "0"])?;

// outcome.reports[0].bytes   the bytes of `-outfmt 6`
// outcome.reports[1].bytes   the bytes of `-outfmt 0`
// outcome.hsps[i]            the numbers of one HSP
// outcome.query_stats[q]     the statistics of query q
```

The repository has a runnable program that does this for any search program. From the `LOSAT/` directory:

```bash
cargo run --release --example library_search -- blastp query.faa subject.faa -evalue 1e-5
```

It prints the `-outfmt 6` report.

## Python

`pylosat` is a Python binding of `losat::api`, built with PyO3. A search runs in the Python process, with the command line's options and report bytes. The HSPs come as an Arrow table. Ctrl-C stops a search between query batches. It supports CPython 3.10 and later.

`pylosat` is not on PyPI. Build a wheel from the repository's `python/` directory with [maturin](https://www.maturin.rs/) (`maturin build --release`, then install the wheel with `pip`), or for development `maturin develop`. The repository's `Python wheels` workflow also builds a wheel for each platform of the release binaries.

```python
import polars
import pylosat

search = pylosat.Search("blastn", task="megablast", evalue=1e-5, num_threads=4)
query = pylosat.Sequences.from_fasta("query.fna")
subject = pylosat.Sequences.from_fasta("subject.fna")

outcome = search.run(query, subject, outfmt=["6", "0"])
outcome.reports[0].text           # the bytes of `-outfmt 6`, as text
outcome.reports[1].bytes          # the bytes of `-outfmt 0`
polars.DataFrame(outcome.hsps)    # one row per HSP, the numbers not rounded
outcome.query_stats[0].eff_space  # what fixed the E-values of the first query
outcome.diagnostics               # what the command line writes to its standard error
```

What `pylosat` offers:

- **Options.** `Search(program, **options)` takes the command line's options as keywords without the dash (`evalue=1e-5` is `-evalue 1e-5`). `pylosat.describe("blastn")` lists them. `Search.from_argv` takes the words themselves. The command line's parser checks every option.
- **Inputs.** `Sequences.from_fasta(path)`, `Sequences.from_text(text)` and `Sequences.from_pairs([(title, sequence), ...])` for FASTA; `Pssm.from_file(path)` for `-in_pssm`; `Database("name")` for `-db`.
- **Results.** `outcome.hsps` is a table with the names of the `-outfmt` columns. Any Arrow consumer reads it without a copy, for example `polars.DataFrame(hsps)` or `pyarrow.table(hsps)`. `outcome.reports` holds the command line's reports. `outcome.query_stats` holds the statistics of each query.
- **Long searches.** `search.batches(query, subject)` searches one query batch at a time. `search.run_to_file(...)` is the command line's `-out`.
- **Tools.** `pylosat.makeblastdb`, `blastdbcmd`, `dustmasker`, `segmasker` and `convert2blastmask` take the command line's words and return what the tool writes.
- **Errors.** `ArgumentError`, `InputError` and `SearchError` derive from `LosatError`, each with the command line's message.

The package README, `python/README.md` in the repository, has the full description.

### What `pylosat` promises, and what it does not

- The binding never changes your sequences, options or statistics, and never splits one search into several searches whose results are joined.
- Reports are the command line's bytes. The binding does not rewrite, filter or regenerate them.
- The Arrow table, pickling and the Python types are `pylosat` formats. They are not NCBI output, and they never claim to be.
- Python cannot stop a thread inside its process, so a search stops only between query batches.
- Sequences, headers, file names and results are never sent over the network.
- The extension calls no NCBI program and starts no process.

The boundary is set by [the decision record on the Python binding boundary](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-PYTHON-BINDING-BOUNDARY.md).

## In the browser

LOSAT also builds to WebAssembly interfaces for use from a web page, as [gbdraw](https://github.com/satoshikawato/gbdraw) does. These interfaces are for integrations and are not documented here. For the browser application, see [Platforms and WebAssembly](platforms.md).
