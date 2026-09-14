# LOSAT

LOSAT (LOcal Sequence Alignment Tool) is a Rust implementation of NCBI BLAST
local-sequence-alignment behavior for the certified profiles described below.
The project is built for native CLI use and WebAssembly-oriented embedding,
without delegating runtime behavior to NCBI BLAST+ executables or libraries.

NCBI BLAST C/C++ source code is the behavioral authority for LOSAT. NCBI BLAST+
may be used as a validation oracle in tests and release checks, but it is not a
runtime dependency, build dependency, feature fallback, or implementation
component.

## v0.1.0 Scope

The v0.1.0 release candidate is certification-gated. Program-profile
certification is complete for the scopes below, but the release itself remains
unpublished. Exact-SHA artifact certification is emitted by the manual release
workflow and is not copied back into source documentation.

| Area | v0.1.0 status | Notes |
| --- | --- | --- |
| BLASTN `-task blastn` and `-task megablast` | Supported for the certified local profile | The [14-case certification](docs/release/blastn_v0.1.0_certification.md) covers 13 exact source-defined cases and one Version 1.2 source-underdetermined equal-HSP case. |
| BLASTP / LOSATP | Supported for the certified gbdraw local profiles | The [nine-case certification](docs/release/blastp_v0.1.0_certification.md) covers gbdraw P1-P3 local query/subject workflows with standard outfmt 6. |
| TBLASTX / TLOSATX | Supported for the certified gbdraw local profiles | The [20-case certification](docs/release/tblastx_v0.1.0_certification.md) covers gbdraw P1-P2 local query/subject workflows and the approved `-db_gencode` behavior below. |
| Native CLI | Supported candidate on Linux x64, Windows x64, macOS arm64, and macOS x64 | PR 6 certified all 43 declared native contracts against the frozen Linux LOSAT output on each non-Linux target. Exact-SHA release archives are produced and certified only by the final RC workflow. |
| `wasm32-wasip1` serial command build | Supported candidate for directly applicable certified rows | PR 5 established raw-byte equality with native LOSAT for all 41 directly applicable rows: BLASTN 14, BLASTP 7, and TBLASTX 20. This is not generic Wasm support. |
| `wasm32-wasip1-threads` | Experimental | Requires the `wasm-threads` feature and a WASI runtime with thread support. |
| Rust library API | Internal only | No semver-stable API commitment in v0.1.0. |
| Web or embeddable Wasm API | Internal only | Public ABI and memory ownership are not yet release-stable. |

TBLASTX local `-subject` searches intentionally honor `-db_gencode` for subject
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
  -query tests/fasta/LC738874.fasta \
  -subject tests/fasta/LC738875.fasta \
  -outfmt 6
```

Run local BLASTN:

```bash
target/release/LOSAT blastn \
  -query tests/fasta/EDL933.fna \
  -subject tests/fasta/Sakai.fna \
  -task megablast \
  -outfmt 6
```

Run local BLASTP:

```bash
target/release/LOSAT blastp \
  -query tests/fasta/AP027078.faa \
  -subject tests/fasta/AP027131.faa \
  -outfmt 6
```

CLI v2 accepts NCBI single-dash option names only. Old spellings such as
`--query`, `--num-threads`, and `-q` are rejected. Common defaults are
`-evalue 10`, `-num_threads 1`, `-max_target_seqs 500`, and `-outfmt 0`.
For v0.1.0, BLASTP publicly accepts only `-task blastp` (also the default).
Other task values are rejected at CLI parsing. The ordinary E-value default
remains 10; an explicit `-evalue` overrides it.
Thread count 0 is invalid; callers must resolve AUTO before invoking LOSAT.

| Program | Implemented output formats |
| --- | --- |
| BLASTN | 6 and 7, standard fields only |
| BLASTP | 0, 6 and 7; implemented custom fields in 6/7 |
| TBLASTX | 6, standard fields only |

BLASTN and TBLASTX reject the default format 0 until it is implemented; pass
`-outfmt 6` explicitly for the examples above. Unsupported formats or fields
fail instead of silently selecting tabular output.

Protein filtering uses `-seg no`, `-seg yes`, or `-seg "12 2.2 2.5"`.
BLASTN uses `-dust no`, `-dust yes`, or `-dust "20 64 1"`.
`-max_hsps` is optional and must be positive when supplied. Run
`losat blastp -help` (or `--help`) for canonical options and task defaults.
See the [CLI v2 implementation record](docs/cli_v2_migration.md) for capability
boundaries, regression evidence, and validation commands.

## Verification

The program records above, the
[integrated native/serial-Wasm record](docs/release/pure_rust_runtime_v0.1.0_certification.md),
and the PR 6 cross-platform certificate are the support authorities. The
integrated result is 43/43 policy-accepted native contracts and 41/41 directly
applicable serial-Wasm/native byte equalities. PR 6 run `33625511701` produced
`CROSS_PLATFORM_NATIVE_CERTIFIED` for Windows x64, macOS arm64, and macOS x64.
The committed program gate entry points are:

```bash
cd LOSAT
cargo build --release --locked
python3 tests/compare_blastn_parity.py \
  --manifest tests/blastn_parity_manifest.tsv \
  --fresh-paired \
  --paired-output-dir /tmp/losat-blastn-v010-certification/paired-base \
  --losat-bin target/release/LOSAT \
  --ncbi-bin /home/kawato/tools/ncbi-blast-oracle/ncbi-blast-2.17.0+/bin/blastn
python3 tests/certify_blastn_v010.py \
  --manifest tests/blastn_parity_manifest.tsv \
  --paired-output-dir /tmp/losat-blastn-v010-certification/paired-base \
  --exceptions tests/blastn_v010_source_exceptions.tsv

cd ..
python LOSAT/tests/audit_blastp_v010.py \
  --output-dir /tmp/losat-blastp-v010-audit/final-native
python LOSAT/tests/audit_tblastx_v010.py \
  --output-dir /tmp/losat-tblastx-v010-certification/final-native
```

Older broad comparison scripts remain useful diagnostics, but they are not the
v0.1.0 support authorities.

NCBI BLAST+ is allowed only in these comparison and diagnostic workflows. LOSAT
runtime code must fail explicitly for unsupported behavior rather than invoking
NCBI tools as a fallback.

## Benchmark

### Hit Distribution

![v0.1.0 hit distribution](benchmarks/v0.1.0/hit_distribution.png)

The current v0.1.0 certified alignment snapshot covers all 43 declared
contracts and compares LOSAT with NCBI BLAST+. See the
[complete methodology and provenance](benchmarks/v0.1.0/README.md).

### Execution Time

![v0.1.0 execution time](benchmarks/v0.1.0/execution_time.png)

The execution-time figure uses exact v0.1.0 benchmark SHA
`af3e2ea837afdb8a00cf19920f68be4f0bf3bfb5`, NCBI BLAST+ 2.17.0, and six
representative cases across NCBI/native n1/n8 plus serial and requested-n8
threaded Wasm modes. It shows the median and min–max range of all five retained
timed repetitions after one protocol warmup. This was a same-machine controlled
benchmark collected in two execution segments, with binary/toolchain identities
revalidated after restart; six restart-conditioning warmups were excluded from
timed statistics. Results are machine-specific and are not a cross-platform
performance guarantee, and this characterization does not expand the certified
feature or target scope above.

## WebAssembly

The standard Wasm builds use `wasm32-wasip1-threads` with the `wasm-threads`
feature. The same command artifact supports both `-num_threads 1` and multiple
threads. A compatible runtime with shared-memory/thread support is required.
This build default does not expand the historical certification scope above.

```bash
cd LOSAT
cargo build-wasi-command
node tests/run_losat_wasi_threads.js target/threaded-command/wasm32-wasip1-threads/release/LOSAT.wasm --help
```

Use `cargo build-web` (also named `build-web-threaded-reactor`) for the separate
threaded reactor. `_start` commands and `_initialize` reactors are not
interchangeable. The reactor links the selected Rust toolchain's
`crt1-reactor.o`; initialize it once with `WASI.initialize`, then call the direct
API. Build and inspect both standard artifacts with
`python tests/build_wasi_artifacts.py --target-dir target --output-dir ../.tmp/wasi-artifacts`.
The helper records kind/imports/exports, hashes and build metadata. Use a fresh
output directory when changing the selected artifact kinds.

Serial Wasm remains an explicit compatibility option:
`cargo build-wasi-command-serial`, `cargo build-web-serial`, or
`build_wasi_artifacts.py --include-serial` (all four artifact kinds).
It uses the separate `run_losat_wasi.js` host and rejects multiple threads.
The current gbdraw browser integration still uses this serial path for some AUTO
workloads and fallback; switching its deployment to threaded-only is a separate
consumer migration. Existing frozen v0.1.0 release/certification workflows retain
their historical serial contract.

`-num_threads 1` runs on the caller without a compute pool. A supported
`-num_threads N` creates exactly N dedicated compute workers for that search and
joins them before returning. Unsupported targets, excessive requests, malformed
caps, and spawn failures return errors. Explicit `LOSAT_WASI_THREAD_CAP` is a
rejection limit; it never silently reduces N. Input size does not override N.

Threaded command/reactor tests use Rust 1.92.0 and Node 24.21.0. See
[threading tests](LOSAT/tests/README.md#threading-contract-gates) and the
[remediation plan](docs/wasm_threading_remediation_plan_20260913.md).
Browser-facing or embeddable Wasm APIs remain internal and are not release-stable
in v0.1.0; these gates do not expand the frozen release certification scope.

## Limitations

- LOSAT is focused on local `-query`/`-subject` comparisons, not large database
  searches against external BLAST databases.
- Unsupported options must be treated as unsupported, not silently delegated to
  NCBI BLAST+.
- BLASTN support is limited to the local query/subject `megablast` and `blastn`
  profile in the committed certification manifest. `dc-megablast`, database
  search, and threaded-Wasm BLASTN certification remain outside this claim.
  The one demonstrated source-underdetermined equal-HSP tie is governed by
  [Product Decision Version 1.2](docs/product_decisions/PD-BLASTN-HSP-CANONICALIZATION.md)
  and the [durable certification record](docs/release/blastn_v0.1.0_certification.md).
- BLASTP support is limited to the certified gbdraw P1-P3 local query/subject
  profiles with standard outfmt 6. Database/remote search, alternate tasks and
  options, other output formats, and threaded Wasm remain outside this claim.
- TBLASTX support is limited to the certified gbdraw P1-P2 local query/subject
  profiles with standard outfmt 6 and one thread per job. Other search modes,
  output formats, unexercised options, and threaded Wasm remain outside this
  claim.
- Serial Wasm evidence covers the 41 directly applicable rows in the declared
  program manifests; it does not promote unlisted options, threaded Wasm, or a
  browser/embeddable ABI.
- Existing comparison outputs under `LOSAT/tests/*_out` are release hygiene
  targets; regenerated output should be treated as artifact or scratch data
  unless explicitly documented as canonical fixture data.

## Release Documents

- Scope: [docs/v0.1.0_scope.md](docs/v0.1.0_scope.md)
- Readiness plan: [docs/v0.1.0_release_readiness_plan.md](docs/v0.1.0_release_readiness_plan.md)
- Release note draft: [docs/release/v0.1.0.md](docs/release/v0.1.0.md)
- Exact-SHA RC contract: [docs/release/v0.1.0_rc_contract.json](docs/release/v0.1.0_rc_contract.json)
- Release procedure: [RELEASE.md](RELEASE.md)
- Changelog: [CHANGELOG.md](CHANGELOG.md)
- Contributing: [CONTRIBUTING.md](CONTRIBUTING.md)
- Security policy: [SECURITY.md](SECURITY.md)

## References

- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [NCBI BLAST source](https://github.com/ncbi/ncbi-cxx-toolkit-public)

## License

[MIT License](LICENSE)
