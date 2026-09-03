# TBLASTX v0.1.0 supported-profile certification

## Decision

`CERTIFIABLE_AS_IS`

At `BASE_SHA` `4d7608cc7d044932b07c97000dfbd36db4d96a95`, LOSAT's TBLASTX behavior is certified only for the local query/subject profiles currently required by gbdraw v0.1.0. This is not a claim of generic TBLASTX compatibility.

The machine-checkable certification gate was added by commit `d6071f5408858bebdcced83c9faa0748628e8a9d` (`Add TBLASTX v0.1.0 certification gate`). No production implementation or biological fixture changed.

## Frozen references

| Reference | Identity |
| --- | --- |
| LOSAT integration base | `fix_blastn` at `4d7608cc7d044932b07c97000dfbd36db4d96a95`; merged PR #86 is at the same commit |
| Official executable oracle | `/home/kawato/tools/ncbi-blast-oracle/ncbi-blast-2.17.0+/bin/tblastx` |
| Oracle version | `tblastx: 2.17.0+`; `Package: blast 2.17.0, build Jul 1 2025 08:59:18` |
| Oracle SHA-256 | `583e5d60bbd444ac455d20e0956c5aa0aeef675da8daee8204d8f9376ddb8804` |
| Read-only official-source checkout | `/mnt/c/Users/genom/GitHub/ncbi-blast`, commit `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4` (`ncbi-blast-2.17.0+-src`) |
| Audited gbdraw checkout | `/mnt/c/Users/genom/GitHub/gbdraw`, commit `cc1ac4f49ba93bcee288072682bccea2e0ff2551` |

The NCBI checkout was inspected read-only for provenance. NCBI was not built from source. Because the fresh source-defined gate had no mismatch, no mismatch-specific source audit was needed.

## Exact supported profile

Current executable gbdraw code defines two required TBLASTX profiles:

| Profile | gbdraw use case | topology and endpoint roles | query gencode | db gencode | outfmt | thread mode | other LOSAT options | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| P1 | Browser-generated circular similarity rings | Comparison FASTA is the TBLASTX query; the displayed circular reference record is the local subject | Per-comparison-record code | Displayed-reference code | Standard `6` | One thread per job | Only `--query-gencode` and `--db-gencode`; remaining search options use LOSAT defaults | `COVERED` |
| P2 | Browser-generated linear pairwise TBLASTX links | The current comparison plan selects local query/subject record endpoints and direction | Query-record code | Subject-record code | Standard `6` | One thread per job | Only `--query-gencode` and `--db-gencode`; remaining search options use LOSAT defaults | `COVERED` |

gbdraw may schedule separate pair jobs concurrently, but each TBLASTX invocation is serial (`--num_threads 1`). It does not request a TBLASTX matrix, E-value, word-size, SEG, strand, maximum-target, or maximum-HSP override. For this profile, LOSAT therefore uses E-value `10`, BLOSUM62, word size `3`, query SEG defaults `12/2.2/2.5`, window size `40`, and maximum targets `500`. gbdraw's downstream identity, E-value, bit-score, and alignment-length display filters do not alter the TBLASTX invocation and are outside this software-output certification.

Coverage includes local query/subject input, forward and reverse pair directions, a self comparison, implicit and explicit default genetic codes, non-default query code with default subject code, the approved non-default subject-code behavior, standard 12-column outfmt 6 values, raw row ordering, IDs, native repeatability, and representative serial Wasm execution. Remote/database search, other output formats, CLI options not issued by gbdraw, native thread-count equivalence, threaded Wasm, performance, public API stability, and biological conclusions are excluded.

## Manifest and fresh native result

The manifest is `LOSAT/tests/tblastx_v010_parity_manifest.tsv`: 20 committed cases covering P1 and P2, comprising 14 source-defined parity cases and 6 approved local-subject `db-gencode` deviation cases. The runner is `LOSAT/tests/audit_tblastx_v010.py`; its classifier/contract tests are in `LOSAT/tests/test_audit_tblastx_v010.py`.

Fresh results were written outside the repository under `/tmp/losat-tblastx-v010-certification/final-native/`:

| Classification | Count | Contract result |
| --- | ---: | --- |
| `EXACT_TEXT` | 14 | All source-defined cases passed as raw byte-for-byte, order-preserving matches |
| `HSP_SET_DIFF` | 6 | All and only the declared non-default local-subject `db-gencode` cases passed the approved-deviation contract |

There were no `EXACT_DATA`, `ORDER_ONLY`, `VALUE_DIFF`, `MISSING_HSP`, execution-error, or other unexpected non-exact cases. All standard outfmt 6 rows had 12 fields and the expected query/subject IDs. Sorting was not used to accept parity.

The targeted `LC738874.fasta` query versus `LC738875.fasta` subject case produced 4,319 rows. NCBI and LOSAT were raw-byte identical with SHA-256 `86c05a04efb50e4026720e2d44fe2db2e6446f9594e174f3fde56931d09d5b49`. Omitting both genetic-code arguments and explicitly passing default code `1/1` produced this same hash for both executables.

## Approved local-subject `db-gencode` deviation

The only accepted product deviation is the repository policy that local `-s/--subject` searches honor an explicit non-default `--db-gencode` during subject translation, search, and reporting. The evidence isolates that boundary:

- all six declared cases with non-default subject code were classified `HSP_SET_DIFF`, not silently promoted to exact parity;
- `p14_ap027131_ap027133_query4` (`query-gencode 4`, `db-gencode 1`) was raw-byte identical to NCBI with 14,871 rows, isolating non-default query-code behavior from the deviation;
- `d06_ap027131_ap027133_db4` (`query-gencode 1`, `db-gencode 4`) was the subject-only deviation probe: NCBI produced 12,672 rows and LOSAT 13,644, classified `HSP_SET_DIFF` as declared;
- `d06` was raw-byte repeatable across three native runs with LOSAT SHA-256 `719339294d0f22d0f415dd33106a0fd93c25ac372b1456d19d92390fce081df1`;
- the existing engine test `db_gencode_controls_local_subject_search_translation` passed, and the production code at the certified base constructs subject frames from `db_code` and reuses that selected frame set for scoring, reevaluation, identity, and output.

No second implicit deviation was observed or permitted by the gate. The default subject code remains in the source-defined parity contract.

## Repeatability and Wasm

Every one of the 14 parity cases plus the subject-only deviation probe `d06` ran three times natively. All 15 cases were raw-byte stable across the three runs. Native thread-count comparison is not applicable because both required gbdraw profiles invoke each TBLASTX job with one thread.

The serial command-Wasm build succeeded with:

```text
cargo build --release --target wasm32-wasip1 --no-default-features --locked
```

Artifact `target/wasm32-wasip1/release/LOSAT.wasm` had SHA-256 `780d59f6b18fd50b5d8fce7e7c56b378dde4961b2c7ea64d870c27cd56491733`. Node's WASI preview1 runtime executed three required probes, and each Wasm output was raw-byte identical to its native LOSAT output:

| Probe | Rows | Native/Wasm SHA-256 | Oracle relation |
| --- | ---: | --- | --- |
| `p03_mela_pemojnva` | 4,794 | `bbef153094025db89dd4a69b2a3c4246b0b8853396fe740bc7e1bd654db105b6` | Raw-byte identical to NCBI |
| `p12_lc738874_lc738875_default` | 4,319 | `86c05a04efb50e4026720e2d44fe2db2e6446f9594e174f3fde56931d09d5b49` | Raw-byte identical to NCBI |
| `d06_ap027131_ap027133_db4` | 13,644 | `719339294d0f22d0f415dd33106a0fd93c25ac372b1456d19d92390fce081df1` | Matches native approved-deviation output |

Threaded Wasm is not part of this certification.

## Quality gates

All gates were rerun on the certification branch after adding the machine-checkable gate:

| Gate | Result |
| --- | --- |
| `cargo fmt -- --check` | PASS |
| `cargo test --release --locked tblastx` | PASS: 65 focused unit, 1 CLI smoke, and 76 focused integration tests; 1 existing integration test ignored |
| `cargo test --release --locked` | PASS: 399 library tests, 5 CLI tests, 4 integration-smoke tests, and 177 integration tests; 3 existing tests ignored in total |
| `cargo clippy --release --all-targets --locked -- -D warnings` | PASS |
| `cargo build --release --locked` | PASS |
| serial `wasm32-wasip1` build | PASS; existing lib/bin output-name collision warning only |
| `python -m unittest LOSAT/tests/test_audit_tblastx_v010.py` | PASS: 9 tests |
| fresh native certification runner | PASS: all 20 manifest contracts |
| serial Wasm runtime comparisons | PASS: all 3 probes |

Generated outputs are intentionally uncommitted and remain under `/tmp/losat-tblastx-v010-certification/`.

## Independent audit

The required publication audit passed. Its read-only checklist and evidence are stored outside the repository at `/tmp/losat-tblastx-v010-certification-audit.md`.

## Limits of this certification

This record certifies only the two gbdraw v0.1.0 local query/subject TBLASTX profiles above at the frozen base and oracle identities. It does not certify generic TBLASTX behavior, database or remote searches, non-standard output formats, unexercised CLI options, multi-thread native equivalence, threaded Wasm, performance, biological-input correctness, or downstream gbdraw filtering and visualization semantics.
