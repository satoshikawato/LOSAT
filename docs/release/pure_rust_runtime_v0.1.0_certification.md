# Pure-Rust runtime v0.1.0 integrated certification

## Decision

`INTEGRATED_RUNTIME_CERTIFIED`

LOSAT commit `5845d22ed9842449628a647f29b8c6762511ca59` completed the declared Linux native and directly applicable serial `wasm32-wasip1` v0.1.0 evidence matrix, artifact reproduction, and repeatability resume. Authority reconciliation and fresh independent reconsideration returned disposition `A — NO POLICY CONFLICT`; all 43 declared native contracts are policy-accepted under their exact classifications.

This PR makes no BLAST production-algorithm change. Its integration base is `1f5dd1bcd215f52e7e23008f74c2298dfac5a8fb`.

## Frozen toolchain and artifacts

| Item | Frozen identity |
| --- | --- |
| Rust | `rustc 1.92.0 (ded5c06cf 2025-12-08)`, LLVM `21.1.3` |
| Cargo | `cargo 1.92.0 (344c4567c 2025-10-21)` |
| Node | `v18.19.1` |
| Host/target | `x86_64-unknown-linux-gnu` |
| Native build | `cargo build --release --locked`; default features; empty output-affecting environment flags |
| Serial-Wasm build | `cargo build --release --target wasm32-wasip1 --no-default-features --locked`; `+simd128` |
| Cargo config | `LOSAT/.cargo/config.toml`, SHA-256 `dfc8cebcb5b37d2f9863387934d74a8e09bdbea13d28db7ae3846070f7dc2d2e` |
| Native artifact | SHA-256 `058587655963bf1eb253f232fdb5553819105646c0d1c50bbece96e410446819` |
| Serial-Wasm artifact | SHA-256 `000691b99169ddca8072954fdfa64348ad94b75546432655f6a05daea1d5eb47` |

The exact-SHA `/tmp` build worktree disappeared after shutdown. A fresh detached worktree was created at the same certification SHA, the retained `CERT_TOOLCHAIN` and build settings were verified first, and only the two frozen release artifacts were rebuilt. Both reproduced hashes match the original hashes byte-for-byte. The runtime certification SHA was not changed and the completed matrix evidence was not invalidated.

## Oracle provenance

The retained oracle identity record uses official NCBI BLAST+ 2.17.0 artifacts only:

| Oracle | SHA-256 |
| --- | --- |
| `blastn` | `33b64bc67d3149cee2459b2f7766b363323df632cf12c099546de00aea9698b5` |
| `blastp` | `5ce267c04e4988c265357bfbedc64e809545b6fcfae7ff6775266fabbee8ba0e` |
| `tblastx` | `583e5d60bbd444ac455d20e0956c5aa0aeef675da8daee8204d8f9376ddb8804` |
| `ncbi-blast-2.17.0+-src.tar.gz` | `502057a88e9990e34e62758be21ea474cc0ad68d6a63a2e37b2372af1e5ea147` |

No NCBI executable was rerun during resume. NCBI remains a comparison oracle only and is not a LOSAT runtime, build, fallback, FFI, or subprocess dependency.

## Zero-delegation boundary

The exact-SHA boundary checker reports:

- temporary project-owned production findings: `0`;
- reviewed dependency `links` metadata signals: `2` (`rayon-core` and `wasm-bindgen-shared`);
- unexpected, stale, or invalid findings: `0`.

The two dependency signals are reviewed non-algorithm metadata, not BLAST implementation delegation.

## Integrated matrix

The completed native classifier evidence was preserved and re-aggregated without rerunning any matrix or oracle command:

| Program | Native contracts | Classification |
| --- | ---: | --- |
| BLASTN | 14/14 policy-accepted | 13 `EXACT_TEXT`; 1 unchanged narrow `SOURCE_UNDETERMINED_ACCEPTED` contract |
| BLASTP | 9/9 | 9 `EXACT_TEXT` |
| TBLASTX | 20/20 | 14 `EXACT_TEXT`; exactly 6 approved local-subject non-default `db-gencode` `HSP_SET_DIFF` cases |
| Total | 43/43 policy-accepted | PASS |

`Sakai.MG1655.megablast` is not byte-exact to the frozen NCBI binary and is not an intentional LOSAT product divergence. It is the existing frozen BLASTN source-underdetermined comparator-tie contract. Its retained five-row (five coordinate-key) difference remains fully inside the authorized footprint: both outputs retain the same 6476 rows, HSP membership and coordinate keys, E-values, and bit scores; only `pident`, `length`, `mismatch`, and/or `gapopen` differ. No classifier, footprint, expected output, or contract was broadened or rebaselined.

The authority reconciliation distinguishes this source-underdetermined NCBI behavior from the narrow TBLASTX intentional product deviation. The committed BLASTN Product Decision, frozen source-exception row, and explicit PR 5 requirement authorize 13 exact BLASTN rows plus this one unchanged narrow source-underdetermined row. The six TBLASTX local-subject non-default `db-gencode` rows remain the only approved intentional LOSAT deviations.

The serial-Wasm evidence is raw-byte equality against the corresponding native LOSAT output:

| Program | Directly applicable serial-Wasm rows | Result |
| --- | ---: | --- |
| BLASTN | 14/14 | `EXACT_BYTES` |
| BLASTP | 7/7 | `EXACT_BYTES` |
| TBLASTX | 20/20 | `EXACT_BYTES` |
| Total | 41/41 | PASS |

The two native BLASTP four-thread rows are each byte-identical to their paired native serial row: 2/2 PASS. Threaded Wasm is not claimed.

The TBLASTX implicit/explicit default-code probe remains PASS with common SHA-256 `86c05a04efb50e4026720e2d44fe2db2e6446f9594e174f3fde56931d09d5b49`.

## Repeatability and checkpoint resume

The persistent checkpoint authority is:

```text
/mnt/c/Users/genom/LOSAT-certification-evidence/
losat-pr5-integrated-5845d22-20260831-final
```

Before resume, `CHECKPOINT.txt` matched SHA-256 `7cdac0ea675de9c0ab248a065321aca76589198268e63568adfed941b92f4225`, and the complete restored tree matched pre-resume digest `33be791e47af207d0cba7280bb85840d13618cf95272a2da6938682e730a28f6`.

Completed run-1 and earlier repeatability evidence was not regenerated. Exactly the ten missing executions were run: p11 Wasm runs 2–3; p14 native and Wasm runs 2–3; and d06 native and Wasm runs 2–3. No representative run 1, NCBI oracle command, completed matrix row, completed repeatability run, or benchmark was rerun.

All six semantic representatives are byte-repeatable across three runs on both supported targets:

| Representative | Native runs 1–3 | Serial-Wasm runs 1–3 |
| --- | --- | --- |
| BLASTN ordinary exact | `fcd9236c7b02488d868feaf35166a472c73e9f9e6cc45fc6d7449a9edc7d4120` | same |
| BLASTN equal-HSP | `3f27c1f1396b59ecb7e78827b5da390e9c97a276c05cbc7a6c82f56bd0a03262` | same |
| BLASTP representative | `fd4b010800e32ce6c823cb38b42a10b7845f3342edae892acccc8f554f9edf34` | same |
| TBLASTX p11 linking-heavy exact | `1eb11f5caa4d1030a016e67bafc90acaecb0a1fe76f9b11256fe1343f0571fe4` | same |
| TBLASTX p14 query-gencode exact | `ef33958a8c31146eea7a2b7fea4e7ff1c690e451058c007f1a0cb60047a43de6` | same |
| TBLASTX d06 approved db-gencode deviation | `719339294d0f22d0f415dd33106a0fd93c25ac372b1456d19d92390fce081df1` | same |

## Historical performance lineage

PR 5 introduces no production algorithm change, so the historical PR 2–4 performance statements were retained without mechanical reruns:

| Migration | Candidate | Retained conclusion |
| --- | --- | --- |
| [PR 91](https://github.com/satoshikawato/LOSAT/pull/91), BLASTN stable sorting | `0a3cf9ad7d6ca9f221e3c9e55e31c008e60cd41c` | measured performance showed no median regression |
| [PR 92](https://github.com/satoshikawato/LOSAT/pull/92), shared output stable ordering | `0652593f113f2813d926415a28b54aaa9ef6efd7` | performance showed no regression |
| [PR 93](https://github.com/satoshikawato/LOSAT/pull/93), TBLASTX search/linking stable sorting | `d92905569317ff9c6bde4e9dabcc66f6e4c15f81` | no material regression: p11 `639.21 s` baseline vs `636.70 s` candidate; d06 `60.26 s` vs `60.65 s` |

These previously accepted PR-level conclusions and their independent read-only approvals remain the performance lineage authority. PR 91 and PR 92 numeric samples are not duplicated in the PR 5 tree; PR 93 retains only the scalar values quoted by its merged PR. That retention limit does not globally invalidate the accepted lineage, and PR 5 introduced no separate production/build change or concrete lineage invalidation. No broader numerical claim was invented and no benchmark was rerun.

## Lightweight quality gates

| Gate | Result |
| --- | --- |
| `cargo fmt --check` | PASS |
| `cargo clippy --all-targets --all-features --locked -- -D warnings` | PASS |
| `cargo test --all-features --locked` | PASS: 177 passed, 0 failed, 2 ignored |
| runtime-boundary checker | PASS: 0 temporary findings, 2 reviewed dependency signals |
| runtime-boundary unit tests | PASS: 7 tests |
| integrated-driver unit tests | PASS: 5 tests |
| `cargo package --list` | PASS |
| `cargo publish --dry-run --locked` | PASS: 320 files, 26.7 MiB, 7.9 MiB compressed; upload aborted by dry-run as expected |

The initial sandboxed package dry-run attempt could not resolve `index.crates.io`; the recorded network-enabled retry passed. No tag, upload, publication, or deployment was performed. The reproduced native and serial-Wasm hashes remained unchanged after the gates.

## Evidence records

The completed resume evidence tree contains the original per-command outputs plus:

- `summary.json` — reconciled certification decision and exact classifications;
- `native_contracts.tsv` — all 43 native classifications and output hashes;
- `native_wasm_equalities.tsv` — all 41 raw-byte equalities;
- `blastp_native_thread_equalities.tsv` — both native thread pairs;
- `repeatability.tsv` — all 12 target rows and run hashes;
- `resume_history.json` — checkpoint restoration and binary-reproduction history;
- `performance_lineage.json` — retained PR 2–4 performance provenance;
- `quality_gates.json` and `quality/` — lightweight gate results and logs.
- `independent_audit.md` — preserved initial `REQUEST_CHANGES`, authority reconciliation, and fresh independent reconsideration disposition A.

The immutable pre-resume persistent checkpoint remains unchanged. The premature audit-blocked tree is also preserved unchanged at:

```text
/mnt/c/Users/genom/LOSAT-certification-evidence/
losat-pr5-integrated-5845d22-20260831-resumed-audit-blocked
```

It is a historical intermediate produced under the superseded policy interpretation and is not the final evidence authority. The corrected final evidence is persisted separately at:

```text
/mnt/c/Users/genom/LOSAT-certification-evidence/
losat-pr5-integrated-5845d22-20260831-resumed-final
```

Its `FINAL_EVIDENCE_MANIFEST.sha256` contains 642 entries. The manifest file SHA-256 is `b9fc98a376d2849274c86b4e4769d2ee38b76025adbb4d63d3da4a5e3e7cdb5c`.

## Independent review and authority reconciliation

A fresh read-only `ncbi_parity_auditor` initially returned `REQUEST_CHANGES`. It independently verified the exact SHA and toolchain, byte-identical reproduced binaries, immutable checkpoint, exactly ten resumed executions, all repeatability hashes, all 41 native/Wasm equalities, both BLASTP thread pairs, zero-delegation result, and lightweight quality logs.

It blocked acceptance because:

1. `Sakai.MG1655.megablast` was interpreted as non-byte-exact and lacking an exception authorized by `AGENTS.md`;
2. the completed evidence needed a distinct persistent path;
3. retained PR 2–4 performance statements are historical and unauditable from the surviving samples.

The initial audit is retained as history. A supervisor-requested reconciliation then reviewed the complete authority chain and found that the first interpretation conflated an intentional product exception with source-underdetermined NCBI behavior. The PR 5 instruction explicitly requires the unchanged narrow BLASTN contract, and the committed BLASTN Product Decision makes one precompiled binary's comparator-equal edit-script choice diagnostic rather than a unique source contract. A fresh independent auditor reconsideration agreed with disposition `A — NO POLICY CONFLICT` and retained the independently accepted PR 2-4 performance lineage.

The final accepted native decision is therefore BLASTN 13 exact + 1 source-underdetermined accepted, BLASTP 9 exact, and TBLASTX 14 exact + 6 approved intentional deviations: 43/43 policy accepted. Serial-Wasm/native equality remains 41/41, BLASTP thread equality remains 2/2, all six repeatability representatives remain PASS, checkpoint/resume integrity and deterministic artifact reproduction remain verified, and no contract was weakened or rebaselined.

## Exclusions and next boundary

This certification covers the declared Linux native and directly applicable serial `wasm32-wasip1` v0.1.0 profiles only. It does not certify threaded Wasm, Windows, macOS, unlisted program/task/options, generic BLAST compatibility, release packaging on every platform, signing, notarization, publication, or biological conclusions.

PR 6 has not begun. Tagging and publication remain separately authorized maintainer actions.
