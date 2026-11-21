# Using LOSAT as a library

You can run LOSAT searches from your own program instead of the command line. From Rust or Python, the search runs inside your process, takes the command line's options, and returns the reports that the command line writes, byte for byte. You also get one typed record per HSP and the statistics of each query, so you do not have to parse the reports. No `losat` executable is started, and no NCBI program is called. A web page can run the same engine in the browser through the npm package `losat`.

## Choose an interface

| Interface | What it runs | Page |
|:---|:---|:---|
| Python package `pylosat` | `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`, `psiblast`, and the tools `makeblastdb`, `blastdbcmd`, `dustmasker`, `segmasker`, `windowmasker`, `convert2blastmask`, `makembindex` | [Python (`pylosat`)](libraries/python.md) |
| Rust crate `losat`, module `losat::api` | The same searches and tools | [Rust (`losat`)](libraries/rust.md) |
| npm package `losat`, for web pages | `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx` on FASTA sequences held by the page | [WebAssembly and Browser](libraries/webassembly.md) |

`rpsblast`, `rpstblastn` and `rpsbproc` run only on the command line. The WASI builds of the command line (`losat.wasm`) are on [Platforms and WebAssembly](platforms.md).

## Availability in v0.2.0

The v0.2.0 release files on [GitHub Releases](https://github.com/satoshikawato/LOSAT/releases) and the [Bioconda](https://anaconda.org/bioconda/losat) package hold the `losat` command line. The Rust crate and the Python package are on their package registries; the npm package is built from the [repository](https://github.com/satoshikawato/LOSAT):

| Interface | How to get it |
|:---|:---|
| Rust crate `losat` | From [crates.io](https://crates.io/crates/losat): `cargo add losat`. |
| Python package `pylosat` | From [PyPI](https://pypi.org/project/pylosat/): `pip install pylosat`. |
| npm package `losat` | Not on npm. Build it from the repository's `web/packages/losat/` directory. |

Each page says how.
