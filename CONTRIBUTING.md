# Contributing

LOSAT targets bit-perfect compatibility with NCBI BLAST behavior. Contributions
are welcome when they preserve that goal and keep unsupported behavior explicit.

## Behavioral Authority

NCBI BLAST C/C++ source code is the only behavioral authority for LOSAT.

- Read the corresponding NCBI source before changing behavior.
- Include NCBI source file paths and line numbers in code comments for behavior
  ports or behavior changes.
- Do not use NCBI BLAST+ executables, libraries, bindings, or subprocess calls
  as LOSAT runtime, build, feature, fallback, or unsupported-feature
  implementation paths.
- If LOSAT has not ported a feature yet, return an explicit unsupported or
  unimplemented error rather than delegating to NCBI.

NCBI BLAST+ may be used only as a validation oracle in tests, diagnostics, and
release parity evidence.

## Scope And Unsupported Behavior

Before implementing a feature, confirm that it exists in NCBI BLAST for the same
program and task. Do not add behavior that has no NCBI equivalent.

For v0.1.0, the release scope is documented in:

- [docs/v0.1.0_scope.md](docs/v0.1.0_scope.md)
- [docs/release/v0.1.0.md](docs/release/v0.1.0.md)

Experimental areas must remain documented as experimental until current fixtures
prove the exact parity being claimed.

## Development Checks

Run focused checks for the area you changed, then run the release gate before a
release candidate:

```bash
cd LOSAT
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
cargo build --release
```

For crates.io packaging changes:

```bash
cd LOSAT
cargo package --list
cargo publish --dry-run --locked \
  --config 'build.target-dir="/tmp/losat-cargo-publish-target"'
```

## Parity Evidence

Behavioral changes need current comparison evidence against NCBI BLAST+.
Record:

- NCBI BLAST+ version.
- LOSAT commit.
- Command lines and input files.
- LOSAT and NCBI output paths.
- Diff summary covering coordinates, raw scores, bit scores, E-values, ordering,
  and formatting.

Use the comparison entry points documented in [README.md](README.md) and
[docs/release/v0.1.0.md](docs/release/v0.1.0.md).

## Repository Hygiene

Do not commit generated comparison output, debug traces, package artifacts,
Python cache files, or temporary files unless they are deliberately promoted to
documented canonical fixtures.

The release source should not include `target/`, `tests/*_out/`, BLAST database
index files, plots, `.tmp/`, `__pycache__/`, `*.pyc`, `*.tmp`, or `*.orig`.

## Pull Request Checklist

- The change is inside the documented release scope, or it documents an
  experimental/unsupported area honestly.
- NCBI source references are included for behavior changes.
- Unsupported behavior fails explicitly.
- No runtime/build/fallback path invokes NCBI BLAST+.
- Focused tests or comparison evidence were run and recorded.
- Generated output was not added to the release tree.
