# Bioconda packaging for LOSAT

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/losat)
[![Conda Version](https://img.shields.io/conda/vn/bioconda/losat.svg?style=flat)](https://anaconda.org/bioconda/losat)

LOSAT is officially available on [Bioconda](https://anaconda.org/bioconda/losat) for Linux (`linux-64`, `linux-aarch64`) and macOS (`osx-64`, `osx-arm64`).

```bash
conda install -c bioconda -c conda-forge losat
# or
mamba install -c bioconda losat
```

This directory maintains the upstream packaging template and build script for `bioconda-recipes`.

## Intended Bioconda platform matrix

Bioconda builds the default x86_64 Linux and macOS platforms and supports the
additional ARM platforms through `extra.additional-platforms`.

LOSAT v0.1.0 is intended to be submitted for:

| Bioconda subdir | Architecture |
| --- | --- |
| `linux-64` | Linux x86_64 |
| `linux-aarch64` | Linux ARM64 |
| `osx-64` | macOS Intel x86_64 |
| `osx-arm64` | macOS Apple Silicon ARM64 |

Windows is not part of the Bioconda platform matrix. Windows distribution must
use a separate channel or GitHub release artifact if added later.

## Files

- `meta.yaml.template` — proposed Bioconda recipe. It contains one deliberate
  placeholder for the SHA-256 of the tagged `v0.1.0` GitHub source archive.
- `build.sh` — native Rust build/install script used by the proposed recipe.

The recipe builds LOSAT from source on each Bioconda worker. It does not embed,
link, invoke, or fall back to NCBI BLAST+.

## After the v0.1.0 tag

1. Download the immutable GitHub tag archive:

   ```bash
   curl -L \
     -o LOSAT-v0.1.0.tar.gz \
     https://github.com/satoshikawato/LOSAT/archive/refs/tags/v0.1.0.tar.gz
   sha256sum LOSAT-v0.1.0.tar.gz
   ```

2. Copy `meta.yaml.template` to `meta.yaml` in a new
   `bioconda-recipes/recipes/losat/` directory.
3. Replace only:

   ```text
   __REPLACE_WITH_V0_1_0_TAG_TARBALL_SHA256__
   ```

   with the observed archive SHA-256.
4. Copy `build.sh` alongside `meta.yaml`.
5. Run the current Bioconda lint/build workflow without changing LOSAT's
   certified program scope.
6. Submit the recipe to `bioconda/bioconda-recipes`.

## Recipe design constraints

The proposed recipe deliberately stays small:

- build from the tagged LOSAT source archive;
- use the Bioconda/conda-forge Rust compiler toolchain;
- install the CLI with `cargo install --locked`;
- generate third-party license metadata with `cargo-bundle-licenses`;
- smoke-test `LOSAT --version`, top-level help, and the three certified program
  subcommands;
- opt in to `linux-aarch64` and `osx-arm64` in addition to Bioconda's default
  x86_64 Linux/macOS builds.

The recipe does not claim generic BLAST compatibility, database/remote search,
threaded Wasm support, or a stable public Rust/Wasm API.

## Pre-tag audit boundary

This upstream repository PR prepares the recipe source-of-truth only. It does
not claim that Bioconda has already built or published all four platforms.
Four-platform publication becomes true only after the downstream Bioconda PR
passes its platform builds and is merged.
