# Installation

LOSAT is one executable, `losat` (`losat.exe` on Windows). Install it from Bioconda, download a pre-built file from the GitHub release page, or build it from source. Searching and building databases need nothing else: no NCBI BLAST+ and no extra libraries. In v0.2.0 the domain annotation of `rpsbproc` reads an archive written by NCBI's `rpsblast` or `rpstblastn`, and NCBI's prebuilt databases are downloaded from NCBI.

## Where LOSAT is published

| Channel | What you get |
|:---|:---|
| [Bioconda](https://anaconda.org/bioconda/losat) | `losat` for Linux (x86_64, aarch64) and macOS (Intel, Apple silicon) |
| [GitHub Releases](https://github.com/satoshikawato/LOSAT/releases) | `losat` for Linux, macOS and Windows, the two WASI modules, and the source crate |
| This repository | Build `losat` from source; the Rust library, the Python binding `pylosat` and the WebAssembly engine for web pages |

LOSAT is not published on crates.io, Homebrew, PyPI, Docker Hub or npm. Bioconda is the only package-manager channel.

## Bioconda

[Bioconda](https://anaconda.org/bioconda/losat) carries LOSAT for Linux (x86_64 and aarch64) and macOS (Intel and Apple silicon). Windows is not part of Bioconda; use the release file below.

```bash
# conda
conda install -c bioconda -c conda-forge losat

# mamba
mamba install -c bioconda losat

# pixi
pixi add losat
```

For `pixi`, the project needs the `bioconda` and `conda-forge` channels. Bioconda gets a release after it is published on GitHub, so it can lag behind. Run `losat --version` to see which one you have; this manual describes v0.2.0.

<!-- example: installation-1 -->
```bash
losat --version
```

This prints `losat 0.2.0` for v0.2.0.

## Pre-built release files

Each release from v0.2.0 on lists eight items on the [Releases](https://github.com/satoshikawato/LOSAT/releases) page. Earlier releases list their own items on their release pages.

| Item | Target | File type | Program inside |
|:---|:---|:---|:---|
| Linux x86_64 | `x86_64-unknown-linux-gnu` | `.tar.gz` | `losat` |
| Linux aarch64 | `aarch64-unknown-linux-gnu` | `.tar.gz` | `losat` |
| macOS arm64 (Apple silicon) | `aarch64-apple-darwin` | `.tar.gz` | `losat` |
| macOS x86_64 (Intel) | `x86_64-apple-darwin` | `.tar.gz` | `losat` |
| Windows x86_64 | `x86_64-pc-windows-msvc` | `.zip` | `losat.exe` |
| WASI serial | `wasm32-wasip1` | `.tar.gz` | `losat.wasm` |
| WASI threads | `wasm32-wasip1-threads` | `.tar.gz` | `losat.wasm` |
| Source crate | none | `.crate` | none |

The archives keep the name `LOSAT` of v0.1.0, for example `LOSAT-0.2.0-x86_64-unknown-linux-gnu.tar.gz`. The program inside is `losat`. The source crate is `losat-0.2.0.crate`. An archive unpacks into a folder of its own name that holds the program, `LICENSE`, `README.md` and `RELEASE-METADATA.json`.

The release also carries `SHA256SUMS`, a software bill of materials (`LOSAT-0.2.0.cdx.json`) and GitHub artifact attestations for the eight files.

### Which file for which platform

| Platform | File | Status in v0.2.0 |
|:---|:---|:---|
| Linux x86_64 | `LOSAT-0.2.0-x86_64-unknown-linux-gnu.tar.gz` | Supported |
| Linux aarch64 | `LOSAT-0.2.0-aarch64-unknown-linux-gnu.tar.gz` | Supported |
| Windows x86_64 | `LOSAT-0.2.0-x86_64-pc-windows-msvc.zip` | Supported |
| macOS arm64 | `LOSAT-0.2.0-aarch64-apple-darwin.tar.gz` | Supported |
| macOS x86_64 | `LOSAT-0.2.0-x86_64-apple-darwin.tar.gz` | Supported |
| WASI serial | `LOSAT-0.2.0-wasm32-wasip1.tar.gz` | Supported |
| WASI threads | `LOSAT-0.2.0-wasm32-wasip1-threads.tar.gz` | Supported: tested in CI, outside the release's NCBI parity claim |

"Supported" without a note means the output matches NCBI BLAST+ 2.17.0 on the release test set; on macOS, the last digit of a full-precision E value can differ by 1. [Platforms and WebAssembly](platforms.md) gives the details for each platform, including the macOS and WASI notes.

### Install a native file

On Linux or macOS, unpack the archive and put `losat` on your `PATH`:

```bash
tar -xzf LOSAT-0.2.0-x86_64-unknown-linux-gnu.tar.gz
mkdir -p ~/.local/bin
cp LOSAT-0.2.0-x86_64-unknown-linux-gnu/losat ~/.local/bin/
~/.local/bin/losat --version
```

If `~/.local/bin` is not on your `PATH` yet, add it (for example `export PATH="$HOME/.local/bin:$PATH"` in your shell's start-up file).

On Windows, unzip the `.zip` file and put `losat.exe` in a folder that is on your `PATH`. The binaries are not code-signed or notarized; see [Troubleshooting](troubleshooting.md) if macOS or Windows stops a downloaded file.

### Check a download

`SHA256SUMS` lists the SHA-256 checksum of each release file, in the format that `sha256sum` reads. Download it next to the files, then check the file you have:

```bash
sha256sum --check --ignore-missing SHA256SUMS
```

On macOS, use `shasum -a 256 --check --ignore-missing SHA256SUMS` (or compare the output of `shasum -a 256 FILE` with the line in `SHA256SUMS`). On Windows PowerShell, `Get-FileHash FILE -Algorithm SHA256` prints the value to compare.

If you use the GitHub command line tool, you can also check the attestation of a file. It shows that the file was built by the repository's release workflow:

```bash
gh attestation verify LOSAT-0.2.0-x86_64-unknown-linux-gnu.tar.gz --repo satoshikawato/LOSAT
```

A second attestation covers the software bill of materials of the same files. Check it with the SBOM predicate type:

```bash
gh attestation verify LOSAT-0.2.0-x86_64-unknown-linux-gnu.tar.gz --repo satoshikawato/LOSAT --predicate-type https://cyclonedx.org/bom
```

## Build from source

You need the [Rust toolchain](https://rustup.rs/) (Rust 1.92 or later; the repository pins 1.92.0). LOSAT is not on crates.io, so `cargo` builds it from the repository:

```bash
git clone https://github.com/satoshikawato/LOSAT.git
cd LOSAT/LOSAT
cargo build --release
```

The program is `target/release/losat` (`target\release\losat.exe` on Windows). It needs no C or C++ compiler.

## The libraries

The libraries are used from this repository; none of them is published on a package registry.

- **Rust.** Depend on the crate `losat` in the repository's `LOSAT/` folder, by path or by git. See [Rust (`losat`)](libraries/rust.md).
- **Python.** Build `pylosat` from the repository's `python/` folder (Python 3.10 or later). See [Python (`pylosat`)](libraries/python.md).
- **Web pages.** The WebAssembly engine for the browser is in `web/packages/losat`. See [WebAssembly and Browser](libraries/webassembly.md).

## The WASI files

The two WASI files are WebAssembly command modules, `losat.wasm`, for hosts such as Node.js. The serial module (`wasm32-wasip1`) runs on one thread. The threads module (`wasm32-wasip1-threads`) needs a host that provides shared memory and threads. The repository's Node.js runners (Node.js 24 or later) handle what the host must pass in; how to run them is on the [Platforms and WebAssembly](platforms.md) page.

To build them yourself, add the two targets with `rustup target add wasm32-wasip1 wasm32-wasip1-threads`, then, from the `LOSAT` folder of the repository:

```bash
# serial WASI command
cargo build --release --bin losat --target wasm32-wasip1 --no-default-features

# threaded WASI command (needs a shared-memory host)
cargo build --release --bin losat --target wasm32-wasip1-threads --features wasm-threads
```

## Next steps

Try the [Quick start](quick-start.md), then the recipes in [Common workflows](workflows.md). If you use NCBI BLAST+ today, read [Coming from NCBI BLAST+](from-ncbi-blast.md) to make `blastn`, `blastp` and the other names run LOSAT.
