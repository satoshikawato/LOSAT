# LOSATP v0.1.0 gbdraw BLASTP profile certification

## Outcome and scope

Primary outcome: `CERTIFIABLE_AS_IS`.

LOSATP is supported only for the declared gbdraw v0.1.0 local query/subject
profiles in this document. This certification does not claim broad NCBI
BLASTP compatibility.

The audit made no BLASTP or shared production-code change. It did not change a
protein fixture or generate a replacement expected-output file.

## Frozen inputs

The audit started only after GitHub reported PR #85, `Certify the BLASTN v0.1.0
supported profile`, merged. The remote `fix_blastn` SHA and clean audit base
were both:

```text
31de524827cbb2c0150918118f29d75537fc2c56
```

The gbdraw product profile was read from the clean current checkout at:

```text
3f7328255300f611d448a6fd070b4a252187ce02
```

The fixed official oracle was:

```text
Path: /home/kawato/tools/ncbi-blast-oracle/ncbi-blast-2.17.0+/bin/blastp
Version: blastp 2.17.0+
Package: blast 2.17.0
Build: Jul 1 2025 08:59:18
SHA-256: 5ce267c04e4988c265357bfbedc64e809545b6fcfae7ff6775266fabbee8ba0e
```

No PATH, Conda, self-built, or alternate BLASTP binary was used. NCBI BLAST
source was not built.

## gbdraw-required profiles

| Profile | gbdraw use case | Task | Search topology | `max_hsps` | `max_target_seqs` | Output | Other search options | Required |
|---|---|---|---|---:|---|---|---|---|
| P1 | Pairwise adjacent protein similarity | `blastp` | Local query/subject, display direction | 1 | Omitted (BLAST default 500) or an explicit positive `candidateLimit` | Standard outfmt 6 | Default BLOSUM62, gap 11/1, word size 3, composition mode 2, SEG off, search e-value 10; positive native thread count or serial Wasm | Yes |
| P2 | Similarity grouping / orthogroup | `blastp` | Local self searches and both directions for every record pair | Omitted/default 0 (all HSPs) | Omitted (BLAST default 500) or an explicit positive `candidateLimit` | Standard outfmt 6 | Same defaults as P1 | Yes |
| P3 | Collinear workflow | `blastp` | Local self searches and both directions for adjacent or all record pairs | Omitted/default 0 (all HSPs) | Omitted (BLAST default 500) or an explicit positive `candidateLimit` | Standard outfmt 6 | Same defaults as P1 | Yes |

gbdraw applies its e-value 1e-5, bit-score 50, identity 70, alignment-length 0,
pairwise display-hit, grouping, and collinearity rules after parsing the raw
BLASTP outfmt-6 rows. Those consumer rules are not BLASTP command options and
are outside this LOSATP search certification.

`candidateLimit=None` means gbdraw does not send a target-limit override. Both
the official oracle and LOSATP then use the BLASTP default of 500 targets. This
document does not describe the omitted setting as an unlimited target search.

The option-by-option inventory is in
`GBDRAW_BLASTP_PROFILE_COVERAGE.tsv`. Every BLASTP search option needed by P1,
P2, or P3 is `COVERED`. gbdraw-owned post-search filters are `OUT_OF_SCOPE`.

## Parity matrix

The committed manifest is `LOSAT/tests/blastp_v010_parity_manifest.tsv`. It has
nine cases and reuses only committed protein fixtures as opaque software test
inputs. The matrix covers:

- P1 default target handling, explicit target limits 1 and 5, one HSP per
  subject, multiple target ordering, query contexts with no retained row, and
  native 1- and 4-thread execution;
- P2/P3 all-HSP forward, reverse, self, default-target, explicit-target, and
  native 1- and 4-thread execution.

The certification command is:

```bash
cd LOSAT
cargo build --release --locked
cd ..
python LOSAT/tests/audit_blastp_v010.py \
  --output-dir /tmp/losat-blastp-v010-audit/final-native
```

The script rejects an oracle path, version, or SHA-256 change; runs the official
oracle and current LOSATP fresh; validates the standard 12 fields; compares row
count, HSP key, coordinates, identity, alignment length, mismatch, gap-open,
E-value, bit score, ordering, and raw bytes; and runs LOSATP three times per
case for repeatability evidence.

## Fresh native result

Classification totals:

```text
EXACT_TEXT  9
EXACT_DATA  0
ORDER_ONLY  0
VALUE_DIFF  0
HSP_SET_DIFF  0
MISSING  0
EXECUTION_ERROR  0
```

| Case | Profiles | Rows, NCBI/LOSATP | No-hit queries, NCBI/LOSATP | SHA-256 | Classification |
|---|---|---:|---:|---|---|
| `pairwise_default_serial` | P1 | 475/475 | 0/0 | `fd4b010800e32ce6c823cb38b42a10b7845f3342edae892acccc8f554f9edf34` | `EXACT_TEXT` |
| `pairwise_limit1_serial` | P1 | 89/89 | 0/0 | `2825ce099d7858d7d23a811bff98f4c42fee96c9519003ee936f25855648d3af` | `EXACT_TEXT` |
| `pairwise_limit5_serial` | P1 | 373/373 | 0/0 | `cd02727a3d7bd2548e1db2759a96b6be72ab0709e7bb0e2179312c01d42419c7` | `EXACT_TEXT` |
| `pairwise_nohit_serial` | P1 | 334/334 | 2/2 | `8c39d7c6c27963e78c01fac416bd86542b64501301919bd74df280e31aa6f4b8` | `EXACT_TEXT` |
| `pairwise_default_thread4` | P1 | 475/475 | 0/0 | `fd4b010800e32ce6c823cb38b42a10b7845f3342edae892acccc8f554f9edf34` | `EXACT_TEXT` |
| `all_hsp_forward_serial` | P2, P3 | 530/530 | 0/0 | `a6d38a3200c9570726550a12858101d1eaa414dca236f925da5651e95d819b5d` | `EXACT_TEXT` |
| `all_hsp_reverse_limit5_serial` | P2, P3 | 445/445 | 0/0 | `2db7a6604c7d44731e69f75546340b127037e63fe56ad9370e7d64849ad7162d` | `EXACT_TEXT` |
| `all_hsp_self_serial` | P2, P3 | 568/568 | 0/0 | `67911f4a93327e9622afe10d15a66d9e52ae8e1c3d905db9ca5fb018d6612e46` | `EXACT_TEXT` |
| `all_hsp_forward_thread4` | P2, P3 | 530/530 | 0/0 | `a6d38a3200c9570726550a12858101d1eaa414dca236f925da5651e95d819b5d` | `EXACT_TEXT` |

The complete generated evidence is outside the repository at
`/tmp/losat-blastp-v010-audit/final-native/`.

There was no non-exact case to classify by software stage. No BLASTP-specific
NCBI source audit or source exception was needed. No BLASTN exception was
carried into BLASTP.

## Repeatability

Every manifest case ran three times. Each case produced one SHA-256 value across
all three runs. The 1-thread and 4-thread variants also matched one another for
both the P1 and P2/P3 raw search profiles. Native repeatability passed.

## Wasm

This required build passed:

```bash
cargo build --release --target wasm32-wasip1 --no-default-features
```

The command-Wasm artifact SHA-256 was:

```text
780d59f6b18fd50b5d8fce7e7c56b378dde4961b2c7ea64d870c27cd56491733
```

The existing Node WASI preview1 command-harness contract ran P1 and P2/P3. The
Wasm, native, and NCBI output SHA-256 values were respectively
`fd4b010800e32ce6c823cb38b42a10b7845f3342edae892acccc8f554f9edf34`
and `a6d38a3200c9570726550a12858101d1eaa414dca236f925da5651e95d819b5d`.

The current gbdraw `runLosatPairDirect` path was also exercised with the current
vendored `browser_wasi_shim`. Its returned P1 and P2/P3 strings had the same two
SHA-256 values. The shim printed `wasi: 0 0` to the Node probe console; that line
was not part of `runLosatPairDirect`'s returned BLASTP text.

These checks certify the serial `wasm32-wasip1` fallback for the declared raw
profiles. They do not certify threaded command Wasm, a browser packaging bundle,
or unrelated BLASTP modes.

## Exclusions

The following are not certified by this record:

- database, remote, or non-local-subject searches;
- `blastp-short`, `blastp-fast`, custom matrices, non-default gaps, SEG-enabled
  searches, other composition modes, ungapped search, or Smith-Waterman
  traceback;
- outfmt 0, 7, custom tabular fields, or any format other than standard outfmt
  6;
- generic BLASTP compatibility outside P1, P2, and P3;
- threaded Wasm runtime behavior.

## Quality and independent review

The audit passed every required gate:

| Gate | Result |
|---|---|
| `cargo fmt -- --check` | Passed |
| `cargo test --release --locked blastp` | Passed |
| `cargo test --release --locked` | Passed |
| `cargo clippy --release --all-targets --locked -- -D warnings` | Passed |
| `cargo build --release --locked` | Passed |
| `cargo build --release --target wasm32-wasip1 --no-default-features` | Passed; existing lib/bin output-name collision warnings only |
| `python -m unittest LOSAT/tests/test_audit_blastp_v010.py` | Passed, 6 tests |
| strict nine-case certification script | Passed, 9 `EXACT_TEXT` and 9 `REPEATABLE` |

The independent read-only audit is stored outside the repository at:

```text
/tmp/losat-blastp-v010-audit-review.md
```
