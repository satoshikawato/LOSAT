# PR01 follow-up evidence — 2026-09-05

Decision: **NOT_PROVEN / NOT_ACCEPTED**. PR02 implementation remains blocked by its unaccepted PR01 prerequisite. The declared follow-up execution is complete; no further timeout increase or later-PR implementation was performed.

## Reconciled correctness

The combined 43 native / 41 directly applicable serial catalog has **74 PASS and 10 TIMEOUT** (native 38/43; serial 36/41). All 74 completed LOSAT outputs match frozen raw bytes/hashes and pass their existing behavioral or exception contracts. No completed-output mismatch was found. The ten timeouts cover both native and serial for p11, d01, d02, d03 and d05.

Compared with the original 68 PASS, 15 TIMEOUT and one oracle failure, six rows are now validated: p06 native/serial, p10 native/serial, and d04 native/serial. The p10 native output is reused; its prior ORACLE_FAIL was due to the missing oracle. This is completed validation, not an algorithm defect fixed.

| Retried case | Oracle process | Native raw + contract | Serial raw + contract |
| --- | --- | --- | --- |
| p06_pemojnva_lvmjnv | PASS | PASS | PASS |
| p10_meenmjnv_mejomjnv | PASS | PASS | PASS |
| p11_avclpv_psclpv | TIMEOUT | TIMEOUT | TIMEOUT |
| d01_nz_self_code4 | PASS | TIMEOUT | TIMEOUT |
| d02_ap027132_nz_code4 | PASS | TIMEOUT | TIMEOUT |
| d03_ap027078_ap027131_code4 | PASS | TIMEOUT | TIMEOUT |
| d04_ap027131_ap027133_code4 | PASS | PASS | PASS |
| d05_ap027133_ap027132_code4 | PASS | TIMEOUT | TIMEOUT |

All catalog oracle executions are now complete except p11 (42/43). An oracle process PASS does not certify an unfinished LOSAT output: eight remaining rows have behavioral validation NOT_RUN, and two p11 rows retain oracle TIMEOUT. d04 uses the existing local-subject db-gencode HSP_SET_DIFF contract, not EXACT_NCBI. Its query-code4/db-code4 command cannot reuse the earlier query-code1/db-code4 protected profile (d06).

Per-row process paths, commands and gates: [combined-audit.json](combined-audit.json). The join rechecks exact selector sets without duplicates, ordered argv/cwd, SHA and byte lengths, process results and the current exception validators. Original rows are retained unchanged.

## Identity, execution budget and preservation

- All 457 pre-existing tracked-file hashes match the starting snapshot. Native, serial and threaded artifacts retain their original SHA; Wasm imports/exports match. Current source, inputs, command runners, canonical/exception authorities and official oracle executables pass preflight. See [identity-check.json](identity-check.json) and [preservation-check.json](preservation-check.json).
- The existing release binaries are reused, with prior build logs and source identity; no rebuild or production edit was needed. The official source archive and four additional references match the inspected source text (only checkout CRLF ignored): [source-reference-check.json](source-reference-check.json). Output bytes are never normalized.
- Original 120-second manifests, raw evidence, samples and canonical expected remain unchanged. A separate [600-second manifest](manifest.json) declared the initial retries. After p11 oracle timed out, the sequential driver was stopped during p11 native; [sequential-interruption.json](sequential-interruption.json) records only owned processes. Partial files remain INTERRUPTED and are not accepted results.
- A further [concurrent manifest](manifest-concurrent.json) declared 21 remaining attempts, each at most 600 seconds, at most three concurrent searches. It reuses the completed p11 oracle timeout without another attempt, the successful p10 native and d04 oracle. All 21 attempts finished with saved process results. The earlier p11 oracle is the 22nd completed follow-up audit execution; the interrupted native attempt is recorded separately.
- Audit concurrency and changing contention invalidate performance comparisons between these retry durations. No retry duration enters a median. No new expected output, exception, sensitivity change or fallback was introduced.

## Reused performance evidence and additional checks

The retained-data exporter passed in a new destination: [retained-revalidation](retained-revalidation/RESULTS.md). All 90 completed cold series (450 timed samples, 630 successful outputs) and 60 successful warm outputs are reused; the nine cold p11 configurations and its warm series remain NOT_RUN because the oracle still times out. Protected profiles remain eight PASS. The original top-level samples, medians and raw files were not regenerated. No eligible case was removed and no whole-eligible-set attainment is claimed.

Focused harness negative tests: **17 PASS** ([log](focused-tests.log)); new orchestration scripts passed syntax inspection and independent read-only review. No production algorithm or existing harness changed, so another broad algorithm execution was not justified beyond the incomplete catalog retries.

The added raw-gated serial BLASTP CPU profile passed: 35,004 bytes, frozen SHA 67911f4a93327e9622afe10d15a66d9e52ae8e1c3d905db9ca5fb018d6612e46. Of 1,475 samples, 997 (67.6%) include the redo workspace function on the stack; this includes 685 alignment leaf samples. It supports the scoped Kappa/redo investigation priority. Exact engine wall-time shares and allocation/reset dominance remain unproved. See [profile summary](blastp-cpu-profile/summary.json) and the appended [HOTSPOTS record](../HOTSPOTS.md). This diagnostic is excluded from latency medians.

## Acceptance and PR02 handoff

Required before complete PR01 acceptance: finish the ten remaining catalog executions and p11 oracle under a new declared budget/environment, then fill p11 cold native1/2/4/8, serial, threaded1/2/4/8 (warmup1 + timed5 + diagnostic1) and serial warm (warmup1 + timed5), revalidate the joined evidence and obtain acceptance review. Existing successful data can continue to be reused only while identity and exact command contracts remain valid.

Native effective-N/active workers stay null where not observed. The protocol permits unavailable values with a reason; these are not an independent universal PR01 or PR02-start blocker, but they prohibit same-effective-N and efficiency claims. Source-derived pool capacity does not fill runtime observations. PR02 must measure jobs/effective/pool/spawn/caps and its own threshold and 0/1/2/many-job fixtures before accepting a policy. Lifecycle/resident work belongs to PR03, finer attribution to PR07/10, real browsers to PR12 and registered platform Gate B to release/PR13. Full classification: [REQUIREMENTS.md](REQUIREMENTS.md).

Independent audit: [auditor-review.json](auditor-review.json). It supports only the stated local follow-up scope, not PR01 acceptance or release certification. No optimization, PR02 implementation, commit, push, publication or external message was performed.

English commit title: `Revalidate PR01 baselines with preserved evidence`. English PR handoff: [PR_DESCRIPTION.md](PR_DESCRIPTION.md).
