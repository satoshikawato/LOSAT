# BLASTN v0.1.0 Certification

Status: native certified on 2026-08-29. Serial Wasm build verified with two
required runtime cases; full-manifest Wasm runtime certification remains
pending.

## Candidate and oracle identity

| Item | Certified value |
| --- | --- |
| LOSAT production commit under test | `458b87dc0bf6e05f8390775f819590de097dce5e` |
| Certification branch | `release/blastn-v010-certification` |
| Native target | `x86_64-unknown-linux-gnu` |
| Rust | `rustc 1.92.0 (ded5c06cf 2025-12-08)` |
| Native LOSAT binary SHA-256 | `f418b85583f3ba5809d578d266203227890079139f904297c890d7aa536a9711` |
| Official oracle path | `/home/kawato/tools/ncbi-blast-oracle/ncbi-blast-2.17.0+/bin/blastn` |
| Official oracle version | `blastn: 2.17.0+` (package `blast 2.17.0`, build 2025-07-01) |
| Official oracle SHA-256 | `33b64bc67d3149cee2459b2f7766b363323df632cf12c099546de00aea9698b5` |
| Official release source | `/home/kawato/tools/ncbi-blast-source/ncbi-blast-2.17.0+-src` |
| Official source archive SHA-256 | `502057a88e9990e34e62758be21ea474cc0ad68d6a63a2e37b2372af1e5ea147` |

The certification PR changes no Rust or other production file. The tested
binary therefore corresponds to the production tree at the recorded base
commit; the PR adds only certification tooling, policy data, and documentation.

## Supported profile

The certified native scope is the 14 cases in
[`LOSAT/tests/blastn_parity_manifest.tsv`](../../LOSAT/tests/blastn_parity_manifest.tsv).
Fixture contents remain opaque test inputs; the manifest is the authority for
their paths and command options.

The supported profile is local query/subject BLASTN with:

- `--task blastn` using reward `2`, penalty `-3`, gap open `5`, gap extend `2`,
  and word size `11`;
- `--task megablast` using reward `1`, penalty `-2`, gap open `0`, gap extend
  `0`, and word size `28`;
- DUST and soft masking enabled, both strands, E-value `10` for the ordinary
  search cases (with `1e-100` only for the compact no-hit probe), maximum target
  sequences `500`, unlimited per-subject HSPs, and native serial execution;
- default 12-field outfmt 6 and 7 output, including multi-query and no-hit
  outfmt 7 behavior represented by the compact cases.

This is not generic BLASTN compatibility. `dc-megablast`, database search,
options absent from the manifest, and threaded native/Wasm output are not
certified by this record.

## Native diagnostic and certification result

Fresh paired outputs were written outside the repository under
`/tmp/losat-blastn-v010-certification/paired-base/`.

The unchanged diagnostic comparator reported:

| Classification | Cases |
| --- | ---: |
| `EXACT_TEXT` | 13 |
| `VALUE_DIFF` | 1 |
| `HSP_SET_DIFF` | 0 |
| `ORDER_ONLY` | 0 |
| `MISSING` | 0 |

The sole `VALUE_DIFF` is `Sakai.MG1655.megablast`. It has 6476 rows in both
outputs, identical coordinate/HSP-key multisets, no NCBI-only or LOSAT-only
key, and identical E-values and bit scores for every common key. Differences
are confined to the four reviewed edit-script-derived fields.

The separate certification gate consumes the committed manifest, fresh paired
outputs, and the one-row exception specification at
[`LOSAT/tests/blastn_v010_source_exceptions.tsv`](../../LOSAT/tests/blastn_v010_source_exceptions.tsv).
It passed with the frozen footprint:

| Invariant | Observed |
| --- | ---: |
| Differing coordinate keys | 5 |
| `pident` differences | 2 |
| Alignment-length differences | 2 |
| Mismatch differences | 2 |
| Gap-open differences | 5 |

Every non-exception case defaults to exact byte equality. The exception fails
on a row-count, coordinate/HSP-set, E-value, bit-score, header, delimiter,
newline, unexpected-field, missing-output, new-case, or expanded-footprint
difference. The raw comparator remains truthful and continues to expose the
actual `VALUE_DIFF`.

## NCBI source provenance

Read-only confirmation used official BLAST+ 2.17.0 release source
`c++/src/algo/blast/core/blast_hits.c:2268-2535`:

- `s_QueryOffsetCompareHSPs` and `s_QueryEndCompareHSPs` omit `gap_info` and can
  return equality for distinct edit scripts;
- `Blast_HSPListPurgeHSPsWithCommonEndpoints` applies C `qsort` with those
  comparators and removes later common-endpoint entries, retaining the first
  survivor.

This certification session did not configure or compile NCBI source. No NCBI
source was linked or used as a LOSAT runtime dependency.

## Native repeatability

The same native candidate binary ran all 14 LOSAT manifest cases three times
with identical inputs and options. Each row below is the SHA-256 observed in
runs 1, 2, and 3.

| Case | Run 1 = run 2 = run 3 SHA-256 |
| --- | --- |
| `PesePMNV.MjPMNV.task_blastn` | `cadccd5ba65c5632193b2866713328ca3689e684ab8718a773df8decf4d8284c` |
| `PmeNMV.MjPMNV.task_blastn` | `30717748c1eac9b994512c1117eb440d1c5a2ec7f05f45dcfe3d50bd3ca13177` |
| `PmeNMV.PesePMNV.task_blastn` | `aa2df5afde8454f4565b1c832c9fdc0b078b3f505a99f8a06112a8e7ada263c1` |
| `PeseMJNV.PemoMJNVB.task_blastn` | `01235a77f9b666c2fba84618d170ff70fe759870c704936480ab348e69a2feb8` |
| `MelaMJNV.PemoMJNVA.task_blastn` | `338eb1fa9f090fbdc4068691080608cdd41af1ee85c57e2c3adaa4d1aaff3e7b` |
| `PemoMJNVA.PeseMJNV.task_blastn` | `19f2bc9c15512dd9ac03947dac6c01868c38acff5e6fcb672601d24dce833a31` |
| `MjeNMV.MelaMJNV.task_blastn` | `780e6fc526f86d258a6830f90f5c593d5136fa5cdea4d4cb750ad4c5506e7bf1` |
| `MjPMNV.MlPMNV.task_blastn` | `ea96f79824e8ecd5a704a6df1950da122eeaee352646ec5cf92594d71ee58c96` |
| `NZ_CP006932.NZ_CP006932.task_blastn` | `f68b187559ec9119b535a0d9e57910061a649ee7ae6d78a7c0b3dd74d990df01` |
| `EDL933.Sakai.megablast` | `4b4b2501f1d87303c3060c5b0c100fbd381604a72e4df130e23701ff81e49fac` |
| `Sakai.MG1655.megablast` | `9edd4883881316976c82c3e3674cab6d4d863e17ac89e3db2052b2ed7996df74` |
| `compact.multi_query.outfmt6` | `ac4a20b64c8da3eaea725ed09dc89880294f58e971b7207c943489faacdbe03e` |
| `compact.multi_query.outfmt7` | `58b0673dc9379427777880c709b9edb35647806e78de74197899680fac402547` |
| `compact.multi_query.no_hit.outfmt7` | `a685460a43cb83c416bd0b3e1f10e4736d12b2f37d44da382907427ef9c6283b` |

The full hash table is retained outside the repository at
`/tmp/losat-blastn-v010-certification/native-repeatability-sha256.tsv`.

## Wasm evidence

`cargo build --release --target wasm32-wasip1 --no-default-features` passed.
The command-Wasm SHA-256 was
`780d59f6b18fd50b5d8fce7e7c56b378dde4961b2c7ea64d870c27cd56491733`.

The existing Node WASI preview1 harness contract was exercised with Node
`v18.19.1`:

| Runtime case | Result |
| --- | --- |
| Source-defined `PesePMNV.MjPMNV.task_blastn` | Official NCBI, native LOSAT, and serial Wasm were byte-identical; SHA-256 `07edabbe23bc1fe31cd9b3b556d2183e8d9df1a0460a82c01de24642ee0a81cc`. |
| Source-underdetermined `Sakai.MG1655.megablast` | Native LOSAT and serial Wasm were byte-identical; SHA-256 `03505c8a5d176bd11157b7bb4dbd3bac5ed2c6931279477b8e2bc343c6ecd066`. The native/NCBI pair retained the certified five-key footprint. |

No native/serial-Wasm incompatibility was observed for the tested cases. This
is targeted runtime evidence, not full 14-case Wasm output certification.
Threaded Wasm BLASTN remains pending, and no canonicalization was added.

## gbdraw profile coverage

Read-only inspection used the local gbdraw checkout at
`64ba172f4c6ba9ffc4749983b10f721f4ee2cdcd`. That checkout had six uncommitted
Web files; their diff changed BLASTP-related option handling and did not alter
the BLASTN option surface described here. No gbdraw file was changed.

| gbdraw-used option/profile | Covered? | Evidence |
| --- | --- | --- |
| Local in-memory query/subject pairs | Yes | The LOSAT Web runtime supplies `query.fa` and `subject.fa`; every certification case is local query/subject. |
| `--task megablast` (default) | Yes | `EDL933.Sakai.megablast` and the narrow `Sakai.MG1655.megablast` case. |
| `--task blastn` | Yes | Nine task-blastn manifest cases plus the compact multi-query cases. |
| `--task dc-megablast` (optional UI choice) | No; pending | Not in the v0.1.0 manifest or declared supported profile. It is not required for the default gbdraw integration. |
| Default 12-field outfmt 6 | Yes | `compact.multi_query.outfmt6`; the same fields feed gbdraw comparison parsing and filters. |
| Outfmt 7 input compatibility | Yes | Manifest task-blastn and megablast cases plus compact multi-query/no-hit cases. gbdraw-generated LOSAT output defaults to outfmt 6. |
| Search E-value `10`, DUST on, both strands, default target/HSP limits | Yes | All ordinary manifest rows use E-value `10`, DUST/soft masking, both strands, target limit `500`, and unlimited per-subject HSPs. Only the compact no-hit probe overrides E-value to `1e-100`; gbdraw does not override the BLASTN search defaults. |
| Downstream E-value, bit-score, identity, and alignment-length filters | Yes as output inputs | The certified default fields include all four inputs. gbdraw filters after LOSAT; it does not pass identity or length as BLASTN search filters. |
| Circular serial fallback (`threadsPerJob: 1`) | Yes | Native three-run evidence and targeted serial-Wasm runtime evidence. |
| Linear threaded per-job Wasm execution | No; pending | The v0.1.0 manifest is native serial and this session did not certify the threaded-Wasm BLASTN path. gbdraw retains a serial fallback. |

The optional `dc-megablast` and threaded paths remain explicitly outside this
certification rather than broadening LOSAT support.

## Commands and quality gates

The principal commands were:

```bash
/home/kawato/tools/ncbi-blast-oracle/ncbi-blast-2.17.0+/bin/blastn -version
sha256sum /home/kawato/tools/ncbi-blast-oracle/ncbi-blast-2.17.0+/bin/blastn
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
```

Quality results:

| Gate | Result |
| --- | --- |
| `cargo fmt -- --check` | Passed |
| `cargo test --release --locked` | Passed: 585 tests, 3 ignored, 0 failed |
| `python3 -m unittest tests/test_compare_blastn_parity.py` | Passed: 5 tests |
| `python3 -m unittest tests/test_certify_blastn_v010.py` | Passed: 10 tests |
| `cargo clippy --release --all-targets --locked -- -D warnings` | Passed |
| `cargo build --release --locked` | Passed |
| `cargo build --release --target wasm32-wasip1 --no-default-features` | Passed; Cargo emitted the existing lib/bin output-name collision warning |
| Fresh raw diagnostic | Passed expected classification: 13 `EXACT_TEXT`, 1 real `VALUE_DIFF`, no other class |
| Fresh certification command | Passed: 14 cases, 13 exact, 1 source-underdetermined |

## Known limitations

- Certification is limited to the committed v0.1.0 local query/subject
  manifest and the exact options described above.
- The accepted equal-HSP residual is source-underdetermined and can affect
  percent identity and alignment length near downstream thresholds. The
  exception does not authorize any coordinate, membership, E-value, bit-score,
  formatter, or unrelated algorithm difference.
- Full-manifest serial Wasm, threaded Wasm, `dc-megablast`, and database search
  certification remain pending.
- Large paired outputs and logs are retained only under
  `/tmp/losat-blastn-v010-certification/` and are not committed.
