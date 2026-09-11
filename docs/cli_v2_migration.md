# CLI v2 migration and local qualification

CLI v2 uses NCBI single-dash, snake_case options under `LOSAT blastn`,
`LOSAT blastp`, and `LOSAT tblastx`. Legacy search aliases are rejected.
The active, unreleased v0.1.0 RC contract and invocation-producing certification
scripts now emit this interface. Historical recorded argv remains readable in
the PR 5 evidence identity reader; it is never accepted by the public parser.

## Canonical source authority (2026-09-11)

The clean parent-derived CLI v2 source is the active implementation authority:
`CANONICAL_CLI_V2_QUALIFIED_SOURCE`. It reconstructs the accepted semantics on
parent `e2abf5848b09309c761c44a373d3bd188daafd80` and is identified by this
commit's Git tree plus the separately retained canonical source manifests.

Both earlier dirty-worktree frozen manifests are
`HISTORICAL_QUALIFICATION_SOURCE`. Their qualification remains historical
semantic evidence; their incidental CRLF bytes are not the final Git target.
The final dirty-worktree source manifest SHA-256 was
`ab6b0fb7b9c62470deedf434fc585c174d1ffa680cc4fc0101c248ebffea82e9`.

```text
old dirty-worktree qualification
    -> semantic reconstruction
clean parent-derived source
    -> fresh qualification
CANONICAL_CLI_V2_QUALIFIED_SOURCE
```

All 13 authorized production files match the accepted non-EOL content. The
TBLASTX `run_impl.rs` change is only the typed SEG adapter; its 3,779 unchanged
lines retain parent LF bytes and its four inserted lines also use LF. The nine
scope-external CRLF-only production files remain byte-identical to the parent.
The retained 41-path inventory below is historical: the canonical commit has
40 changed paths because `packaging/bioconda/README.md` retains the parent's
curl `-o` invocation. Changing curl to `-out` is unrelated to LOSAT CLI grammar.
No non-EOL production hunk was excluded and no algorithm was changed.

Fresh qualification from this clean source passed 607 Cargo tests (3 existing
ignored tests), all-target Clippy, focused rustfmt, 104 Python harness tests,
and shell syntax checks. The accepted biological matrix passed nine native
raw-output equalities, four exact NCBI 2.17.0+ oracle checks, and the three
removed-public-capability parser rejections. This remains the bounded CLI
qualification, not a new cross-platform release certification or performance
claim. No FASTA provenance or release packaging work was performed.

Production manifest SHA-256 (153 Rust files):
`a892a1f91f286155f2fbe5ffb9e20e8cc75d957bb874481b22807f1c0265296a`.

| Clean artifact | Bytes | SHA-256 |
| --- | ---: | --- |
| native | 2909136 | `ea137912d80e439289a8bed47821b75992bd94a557b5759320dbf21fe96f6367` |
| serial | 2051371 | `309b8d2890be345cc58f736ebdd3870dc311ff718445ef10a5ea1657147590e3` |
| threaded | 2288108 | `8ab2cdf8916a6b6e765b0ebe00f9859b7e47f87b8e8911356cd9cf820f15e873` |

`WASM_REPRODUCIBILITY: EXACT`: both clean Wasm builds match the currently
adopted gbdraw hashes. The previous adopted-browser qualification remains
applicable; no new browser campaign or gbdraw adoption was performed.

The canonical full-tree/production manifests, semantic hunk audit, build
records, raw outputs and test logs are retained as separate qualification
evidence. The earlier evidence directories remain untouched. The ordinary
dirty checkout remains at its original branch, HEAD, index and file bytes.

## Current v0.1.0 public surface

The public BLASTP task is `blastp` only, also selected when `-task` is omitted.
Other BLASTP task values are rejected during CLI parsing. Internal task
resolution remains available for future implementation work. The ordinary
E-value default remains 10 and explicit overrides are preserved.

The historical dirty-worktree restriction and adopted Wasm identities are recorded in
[`artifacts/cli-v2-v010-surface-20260911`](../artifacts/cli-v2-v010-surface-20260911/).
The initial migration artifact identities in the following section are historical.

Historical final restriction qualification passed: Cargo 607 passed, 0 failed (3 ignored), Clippy,
9 native byte equalities, and 18/18 candidate plus 18/18 adopted browser checks.
No rebuild followed asset adoption.

| Artifact | Bytes | SHA-256 |
| --- | ---: | --- |
| native | 2909136 | `ea137912d80e439289a8bed47821b75992bd94a557b5759320dbf21fe96f6367` |
| serial | 2051371 | `309b8d2890be345cc58f736ebdd3870dc311ff718445ef10a5ea1657147590e3` |
| threaded | 2288108 | `8ab2cdf8916a6b6e765b0ebe00f9859b7e47f87b8e8911356cd9cf820f15e873` |


## Initial migration qualification (before the v0.1.0 surface restriction)

CLI v2 is accepted for the declared migration scope. Full Cargo tests: 607
passed, 3 ignored. Clippy and focused rustfmt passed. Native regression: 11
migration checks and four NCBI checks passed. Candidate and adopted browser
runs: 18/18 each, zero page/console errors and zero external requests.
Independent source review found no remaining CLI defects. No rebuild followed
asset adoption; source and artifact hashes were rechecked at handoff.

| Artifact | Bytes | SHA-256 |
| --- | ---: | --- |
| native | 2909792 | `69bd921733902a445837fda3bc2b24221d4384c32eb138e63e0ff0537a606109` |
| serial | 2052384 | `c4ea93cb8a3ca63079ce1c95102fd4e9a6b8320855350e0564494f37ce188b5a` |
| threaded | 2288993 | `d4cf0510f4c0588928b3ac0f1cca8a5544628a3c55558f6f4ad19586c589724d` |

Machine-readable final fields: [`acceptance.json`](../artifacts/cli-v2-20260911/acceptance.json).

## Configuration boundary

Both `-query` and `-subject` require file paths. Stdin (`-`) and `-db` are not
implemented. Values are consumed as values, including negative numbers and
paths beginning with a dash. `-help`/`--help` display canonical names; top-level
`--version` remains available.

| Setting | CLI v2 behavior |
| --- | --- |
| `-num_threads` | Integer 1..2147483647; default 1; product AUTO resolves an explicit count |
| `-evalue` | Finite, nonnegative; ordinary default 10 |
| BLASTP task and omitted E-value | Only `blastp`; parser retains `None`; resolve gives 10 |
| `-max_target_seqs` | Positive integer; default 500 |
| `-max_hsps` | BLASTN/BLASTP optional positive integer; omission becomes internal unlimited sentinel 0 |
| `-seg` | BLASTP/TBLASTX shared `no`, `yes`, or one quoted `WINDOW LOCUT HICUT` value |
| SEG defaults | BLASTP no; TBLASTX 12/2.2/2.5; finite cutoffs retain NCBI normalization |
| `-dust` | BLASTN `no`, `yes`, or one quoted `LEVEL WINDOW LINKER`; default 20/64/1 |
| `-comp_based_stats` | Entire-token grammar 0/1/2/3/F/f/D/d/T/t; lowercase `u` only with enabled modes |
| BLASTN tasks | megablast, blastn; engine task-specific scoring resolution preserved |
| BLASTN word size | 4..100 |
| TBLASTX word size | 3 only; other sizes were never implemented by this engine |
| Genetic codes | NCBI valid IDs; query_gencode and db_gencode default 1 |

BLASTN `hitlist_size`, `scan_step`, and `min_diag_separation` are internal fields,
not exposed ineffective controls. Other retained LOSAT engine controls are
identified as such in help. No search algorithm or product scheduling policy
was changed by this migration.

| Program | Output capability |
| --- | --- |
| BLASTN | Formats 6/7, standard fields |
| BLASTP | Formats 0/6/7, implemented custom fields for 6/7 |
| TBLASTX | Format 6, standard fields |

The represented default is format 0. BLASTN/TBLASTX explicitly reject it until
implemented; callers must pass `-outfmt 6` or a supported alternative.

The public task restriction does not remove internal task resolution or its
focused tests. It changes only which configurations the CLI advertises and
accepts for v0.1.0.

## Qualification evidence

Initial migration evidence and logs (before the public restriction) are in
[`artifacts/cli-v2-20260911`](../artifacts/cli-v2-20260911/).
The 153 production Rust files are frozen in `frozen-source.json`.
The original checkout was commit `e2abf5848b09309c761c44a373d3bd188daafd80` plus
pre-existing working-tree changes. This is a local qualification of that
working tree, not a release or a new cross-platform certification.

An environment restart removed the initial `/tmp/losat-cli-v2` raw baseline
and logs. The initial six output hashes were
recorded before the migration and independently reviewed. The current public regression retains the four ordinary-program hashes and
checks the hidden capabilities as parser rejections. It does not claim
a newly rerun pre-migration baseline. All post-restart evidence is durable here.

The ordinary regression matrix comprises BLASTN formats 6/7, BLASTP format 6,
and TBLASTX format 6 at one and two threads: eight byte-equality checks, plus
four exact NCBI 2.17.0+ oracle checks. The two historical additional BLASTP task searches
were migration-equivalence checks before they were hidden from the public CLI. Their original outputs already differ from
NCBI: short-supported has a bit-score whitespace difference; fast-default
has different retained rows. CLI acceptance does not certify those tasks as
universally NCBI-identical. These existing algorithm/format gaps were not changed.

Browser QA uses actual gbdraw command, serial direct API, and threaded service/
worker paths for all three programs, twice each (18 comparisons). Each raw
result must equal the final native binary. The command path also verifies an
invalid thread-count failure followed by success. Fresh contexts block external
network requests, record page/console errors, and require threaded product-path progress.
BLASTN and BLASTP each spawned one additional worker. The short TBLASTX fixture
requested two threads through the threaded module but spawned no additional
worker; this is path qualification, not a parallel speedup claim.
Adopted mode serves the ordinary gbdraw Wasm paths with no candidate substitution.
This tests product runtime services, not an end-to-end SVG generation journey.

## Commands

From `LOSAT/`:

```sh
cargo test --locked --offline
cargo clippy --locked --offline --all-targets
cargo build --release --locked --offline --bin LOSAT
cargo build --release --locked --offline --bin LOSAT --target wasm32-wasip1 --no-default-features
cargo build --release --locked --offline --bin LOSAT --target wasm32-wasip1-threads --features wasm-threads
```

From the repository root, after building:

```sh
python LOSAT/tests/check_cli_v2_regression.py --output-dir artifacts/cli-v2-v010-surface-20260911/regression
python artifacts/cli-v2-v010-surface-20260911/browser_qa.py
# Copy the qualified artifacts, compare SHA-256, and do not rebuild afterward.
python artifacts/cli-v2-v010-surface-20260911/browser_qa.py --adopted
```

Focused rustfmt checks every CLI-modified Rust file listed in `rust-files.txt`
with `--edition 2021 --config skip_children=true`. Whole-crate format retains
the unrelated existing `LOSAT/src/algorithm/blastp/gapalign.rs` formatting issue.
The historical ordinary-worktree v0.1.0 Python run had one pre-existing
native-input identity failure: AP027152/AP027202 had CRLF byte identities
inconsistent with the registered fingerprints. That run preserved authority
and FASTA bytes without bypassing the release gate. The canonical checkout
retains exact parent fixture bytes and its 82-test v0.1.0 harness suite passes.
The ordinary-worktree fixture provenance remains outside this session.

## Exact change inventory

The paths below are this CLI migration's change footprint. Most already had
pre-existing changes; those overlapping changes were preserved. This list does
not assign all `git diff` content in a listed file to CLI v2.

### LOSAT (paths relative to its repository root)

production:

- `LOSAT/src/algorithm/blastn/args.rs`
- `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- `LOSAT/src/algorithm/blastn/coordination.rs`
- `LOSAT/src/algorithm/blastp/args.rs`
- `LOSAT/src/algorithm/blastp/blast_engine.rs`
- `LOSAT/src/algorithm/tblastx/args.rs`
- `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`
- `LOSAT/src/blastinput/mod.rs`
- `LOSAT/src/blastinput/value_parsers.rs`
- `LOSAT/src/cli.rs`
- `LOSAT/src/lib.rs`
- `LOSAT/src/main.rs`
- `LOSAT/src/web_api.rs`

test/harness:

- `LOSAT/tests/audit_blastp_v010.py`
- `LOSAT/tests/audit_tblastx_v010.py`
- `LOSAT/tests/certify_integrated_runtime_v010.py`
- `LOSAT/tests/certify_platform_native_v010.py`
- `LOSAT/tests/check_cli_v2_regression.py`
- `LOSAT/tests/cli_smoke.rs`
- `LOSAT/tests/cli_v2.rs`
- `LOSAT/tests/compare_blastn_parity.py`
- `LOSAT/tests/run_comparison.sh`
- `LOSAT/tests/run_wasm_comparison.sh`
- `LOSAT/tests/test_certify_integrated_runtime_v010.py`
- `LOSAT/tests/test_certify_platform_native_v010.py`
- `LOSAT/tests/test_compare_blastn_parity.py`
- `LOSAT/tests/test_prepare_release_candidate_v010.py`
- `LOSAT/tests/test_wasm_performance.py`
- `LOSAT/tests/unit/blastn/args.rs`
- `LOSAT/tests/unit/tblastx/args.rs`
- `LOSAT/tests/unit/tblastx/run_engine.rs`
- `LOSAT/tests/wasm_performance.py`
- `tests/compare_long_sequences_debug.sh`
- `tests/compare_self_tblastx.sh`
- `tests/compare_tblastx_native_ncbi_parity.sh`
- `tests/compare_tblastx_wasm_parity.sh`

documentation/contract:

- `README.md`
- `docs/cli_v2_migration.md`
- `docs/release/blastp_v0.1.0_certification.md`
- `docs/release/v0.1.0_rc_contract.json`
- `packaging/bioconda/README.md`

### gbdraw (paths relative to its repository root)

production:

- `gbdraw/analysis/protein_colinearity.py`
- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/services/losat-runtime.js`
- `gbdraw/web/js/workers/losat-threaded-worker.js`

test/harness:

- `tests/test_protein_colinearity.py`
- `tests/web/contracts/vibrio-full-generation.serial.spec.js`
- `tests/web/losat-cli-v2.test.mjs`
- `tests/web/run-analysis-derived-cache.test.mjs`
- `tests/web/run-analysis-simple-path.test.mjs`

documentation/contract:

- `gbdraw/web/wasm/losat/README.md`

generated asset:

- `gbdraw/web/wasm/losat/losat-threaded.wasm`
- `gbdraw/web/wasm/losat/losat.wasm`

Generated qualification files are listed individually in the final evidence
initial manifest in `artifacts/cli-v2-20260911/evidence-files.json` and the final
restriction manifest in `artifacts/cli-v2-v010-surface-20260911/evidence-files.json`;
build outputs are
`LOSAT/target/release/LOSAT`, `LOSAT/target/wasm32-wasip1/release/LOSAT.wasm`, and
`LOSAT/target/wasm32-wasip1-threads/release/LOSAT.wasm`.

Unrelated modified LOSAT paths are separately enumerated in
`artifacts/cli-v2-20260911/preexisting-losat-paths.txt`; the gbdraw counterpart is
`preexisting-gbdraw-paths.txt`. Existing algorithm optimizations, host lifecycle,
AUTO selection, T01/T04/T09, FASTA line endings, and older experiment reports
remain outside this migration. The historical migration session performed no
commit or push. Its historical handoff title was `Migrate LOSAT and gbdraw to
CLI v2`, covering canonical grammar, migrated callers/contracts and the then
adopted Wasm assets. The subsequent canonical LOSAT-only commit and authority
are described at the top of this document; no gbdraw adoption or push is part
of the canonicalization session.
