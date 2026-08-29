# LOSAT

LOSAT (LOcal Sequence Alignment Tool) is a Rust implementation of
BLAST-compatible local sequence alignment. The project is built for native CLI
use and WebAssembly-oriented embedding, without delegating runtime behavior to
NCBI BLAST+ executables or libraries.

NCBI BLAST C/C++ source code is the behavioral authority for LOSAT. NCBI BLAST+
may be used as a validation oracle in tests and release checks, but it is not a
runtime dependency, build dependency, feature fallback, or implementation
component.

## v0.1.0 Scope

The v0.1.0 release candidate is parity-gated. A final tag must include current
NCBI comparison evidence for every supported area.

| Area | v0.1.0 status | Notes |
| --- | --- | --- |
| TBLASTX local `-query`/`-subject` | Supported candidate | Must pass native NCBI parity fixtures, except for the approved local-subject `--db-gencode` behavior below. |
| BLASTN `--task megablast` | Experimental | Promotion requires exact source-defined parity and the narrow source-compatible equal-HSP contract documented in the Product Decision. |
| BLASTN `--task blastn` | Experimental | Promotion requires exact source-defined parity and the narrow source-compatible equal-HSP contract documented in the Product Decision. |
| BLASTP | Experimental | The CLI subcommand exists, but BLASTP remains secondary to TBLASTX/BLASTN for v0.1.0. |
| Native CLI | Supported candidate | Release artifacts must include build metadata and checksums. |
| `wasm32-wasip1` serial command build | Experimental | Must match the corresponding native output before being promoted. |
| `wasm32-wasip1-threads` | Experimental | Requires the `wasm-threads` feature and a WASI runtime with thread support. |
| Rust library API | Internal only | No semver-stable API commitment in v0.1.0. |
| Web or embeddable Wasm API | Internal only | Public ABI and memory ownership are not yet release-stable. |

TBLASTX local `-subject` searches intentionally honor `--db-gencode` for subject
translation/search/reporting for every non-default genetic code. This is the
only approved v0.1.0 behavior difference from NCBI BLAST+ local `-subject`
semantics; all other timing, ordering, scoring, filtering, statistics, pruning,
and formatting behavior remains NCBI-parity gated.

## CLI Usage

Build the native CLI:

```bash
cd LOSAT
cargo build --release
```

Run local TBLASTX:

```bash
target/release/LOSAT tblastx \
  -q tests/fasta/LC738874.fasta \
  -s tests/fasta/LC738875.fasta \
  --outfmt 6
```

Run local BLASTN:

```bash
target/release/LOSAT blastn \
  -q tests/fasta/EDL933.fna \
  -s tests/fasta/Sakai.fna \
  --task megablast \
  --outfmt 6
```

Run BLASTP experimental mode:

```bash
target/release/LOSAT blastp \
  -q tests/fasta/AP027078.faa \
  -s tests/fasta/AP027131.faa \
  --outfmt 6
```

Supported output formats currently exposed by the CLI are:

- `0`: pairwise alignment view
- `6`: tabular output
- `7`: tabular output with comment lines

Custom field specifications are supported for formats `6` and `7` where the
selected program implements the requested fields.

## Verification

Release validation compares LOSAT output against NCBI BLAST+ output using local
fixtures. The release gate records the NCBI BLAST+ version, LOSAT commit,
commands, inputs, output paths, checksums, and diff summaries.

Primary comparison entry points:

```bash
cd LOSAT/tests
bash run_comparison.sh
python3 compare_blastn_parity.py
bash run_blastp_comparison.sh

cd ../..
bash tests/compare_tblastx_native_ncbi_parity.sh
bash tests/compare_self_tblastx.sh
bash tests/compare_long_sequences_debug.sh
```

NCBI BLAST+ is allowed only in these comparison and diagnostic workflows. LOSAT
runtime code must fail explicitly for unsupported behavior rather than invoking
NCBI tools as a fallback.

## WebAssembly

The crate can be built for command-style WASI targets. Plain `wasm32-wasip1`
builds are serial. Threaded Wasm builds require `wasm32-wasip1-threads`, the
`wasm-threads` feature, and a compatible runtime.

```bash
cd LOSAT
cargo build --release --target wasm32-wasip1 --no-default-features
cargo build --release --target wasm32-wasip1-threads --features wasm-threads
```

Browser-facing or embeddable Wasm APIs are not release-stable in v0.1.0.

## Limitations

- LOSAT is focused on local `-query`/`-subject` comparisons, not large database
  searches against external BLAST databases.
- Unsupported options must be treated as unsupported, not silently delegated to
  NCBI BLAST+.
- BLASTN and BLASTP are experimental in the v0.1.0 scope. BLASTN promotion
  requires exact parity for source-defined behavior; its demonstrated
  source-underdetermined equal-HSP tie is governed by
  [Product Decision Version 1.1](docs/product_decisions/PD-BLASTN-HSP-CANONICALIZATION.md).
- Existing comparison outputs under `LOSAT/tests/*_out` are release hygiene
  targets; regenerated output should be treated as artifact or scratch data
  unless explicitly documented as canonical fixture data.

## Release Documents

- Scope: [docs/v0.1.0_scope.md](docs/v0.1.0_scope.md)
- Readiness plan: [docs/v0.1.0_release_readiness_plan.md](docs/v0.1.0_release_readiness_plan.md)
- Release note draft: [docs/release/v0.1.0.md](docs/release/v0.1.0.md)
- Release procedure: [RELEASE.md](RELEASE.md)
- Changelog: [CHANGELOG.md](CHANGELOG.md)
- Contributing: [CONTRIBUTING.md](CONTRIBUTING.md)
- Security policy: [SECURITY.md](SECURITY.md)

## References

- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [NCBI BLAST source](https://github.com/ncbi/ncbi-cxx-toolkit-public)

## License

[MIT License](LICENSE)
