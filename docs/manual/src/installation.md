# Installation

LOSAT ships as a single self-contained executable (`losat`, or `losat.exe` on Windows), plus embeddable Rust and Python packages and WebAssembly modules. Searching and building databases require no C/C++ runtime libraries and no local NCBI BLAST+ installation. (For conserved-domain annotation in v0.2.0, `losat rpsbproc` reads `-outfmt 11` archives written by NCBI's `rpsblast` or `rpstblastn`.)

## Release channels at a glance

| Channel | Target | Availability |
| :--- | :--- | :--- |
| **[Bioconda](https://anaconda.org/bioconda/losat)** | Linux (`x86_64`, `aarch64`), macOS (Apple silicon & Intel) | Available with v0.2.0 |
| **[GitHub Releases](https://github.com/satoshikawato/LOSAT/releases)** | Linux, macOS, Windows binaries, WASI modules, `.crate`, SBOM & attestations | Available with v0.2.0 |
| **crates.io (`losat`)** | Rust crate (`losat::api`) & CLI via `cargo install losat` | Publishing right after v0.2.0 |
| **PyPI (`pylosat`)** | Python 3.10+ wheels (`pip install pylosat`) | Publishing right after v0.2.0 |
| **From source** | Any platform with Rust 1.92+ | Repository `HEAD` / release tag |

## Bioconda

[Bioconda](https://anaconda.org/bioconda/losat) packages LOSAT for Linux (`x86_64` and `aarch64`) and macOS (Intel and Apple silicon). For Windows, use the GitHub Release archive below.

```bash
# conda
conda install -c bioconda -c conda-forge losat

# mamba
mamba install -c bioconda losat

# pixi (requires bioconda and conda-forge channels)
pixi add losat
```

Confirm the installed version (Bioconda updates shortly after each GitHub Release):

<!-- example: installation-1 -->
```bash
losat --version
```

## Pre-built release files

Each release on [GitHub Releases](https://github.com/satoshikawato/LOSAT/releases) provides eight archives plus checksums, a CycloneDX SBOM (`LOSAT-0.2.0.cdx.json`), and GitHub artifact attestations:

| Item | Target | Archive | Binary | Status in v0.2.0 |
| :--- | :--- | :--- | :--- | :--- |
| Linux x86_64 | `x86_64-unknown-linux-gnu` | `.tar.gz` | `losat` | Supported |
| Linux aarch64 | `aarch64-unknown-linux-gnu` | `.tar.gz` | `losat` | Supported |
| macOS arm64 (Apple silicon) | `aarch64-apple-darwin` | `.tar.gz` | `losat` | Supported |
| macOS x86_64 (Intel) | `x86_64-apple-darwin` | `.tar.gz` | `losat` | Supported |
| Windows x86_64 | `x86_64-pc-windows-msvc` | `.zip` | `losat.exe` | Supported |
| WASI serial | `wasm32-wasip1` | `.tar.gz` | `losat.wasm` | Supported |
| WASI threads | `wasm32-wasip1-threads` | `.tar.gz` | `losat.wasm` | Supported (tested in CI, outside NCBI parity scope) |
| Source crate | — | `.crate` | — | Source package |

"Supported" denotes matching NCBI BLAST+ 2.17.0 outputs across the release test suite; see [Platforms and WebAssembly](platforms.md) for details.

### Install a native binary

On Linux or macOS, unpack the archive and place `losat` on your `PATH`:

```bash
tar -xzf LOSAT-0.2.0-x86_64-unknown-linux-gnu.tar.gz
mkdir -p ~/.local/bin
cp LOSAT-0.2.0-x86_64-unknown-linux-gnu/losat ~/.local/bin/
losat --version
```

On Windows, extract the `.zip` archive and add its directory containing `losat.exe` to `PATH`. Pre-built binaries are not Apple-notarized or Authenticode-signed; see [Troubleshooting](troubleshooting.md) if macOS Gatekeeper or Windows SmartScreen warns on first launch.

### Check a download

Compare your downloaded archive against `SHA256SUMS`:

```bash
# Linux
sha256sum --check --ignore-missing SHA256SUMS

# macOS
shasum -a 256 --check --ignore-missing SHA256SUMS
```

On Windows PowerShell: `Get-FileHash FILE -Algorithm SHA256`.

With the GitHub CLI (`gh`), you can also check build provenance and SBOM attestations:

```bash
gh attestation verify LOSAT-0.2.0-x86_64-unknown-linux-gnu.tar.gz --repo satoshikawato/LOSAT
gh attestation verify LOSAT-0.2.0-x86_64-unknown-linux-gnu.tar.gz --repo satoshikawato/LOSAT --predicate-type https://cyclonedx.org/bom
```

## Rust crate and Cargo (`losat`)

Once `losat` is published on crates.io immediately after v0.2.0 (Rust **1.92+** required):

```bash
cargo install losat
```

Or build from a local repository clone:

```bash
git clone https://github.com/satoshikawato/LOSAT.git
cd LOSAT/LOSAT
cargo build --release
```

The built binary is `target/release/losat` (`target\release\losat.exe` on Windows). See [Rust (`losat`)](libraries/rust.md) to use `losat::api` in your own Rust crate.

## Python library (`pylosat`)

Pre-built wheels for `pylosat` are scheduled for PyPI right after v0.2.0 (`pip install pylosat`). Until then, build with `maturin` from `python/` or fetch workflow wheels as described in [Python (`pylosat`)](libraries/python.md).

## WebAssembly (WASI modules)

The two WASI modules (`losat.wasm`) run in WebAssembly hosts such as Node.js 24+:

```bash
rustup target add wasm32-wasip1 wasm32-wasip1-threads
cd LOSAT

# serial WASI command
cargo build --release --bin losat --target wasm32-wasip1 --no-default-features

# threaded WASI command (requires host thread/shared-memory support)
cargo build --release --bin losat --target wasm32-wasip1-threads --features wasm-threads
```

See [Platforms and WebAssembly](platforms.md) and [WebAssembly and Browser](libraries/webassembly.md) for host requirements and Node.js launcher scripts.

## Next steps

Follow the five-minute [Quick start](quick-start.md), browse [Common workflows](workflows.md), or read [Coming from NCBI BLAST+](from-ncbi-blast.md) to alias `blastn`, `blastp`, and the other NCBI command names to `losat`.
