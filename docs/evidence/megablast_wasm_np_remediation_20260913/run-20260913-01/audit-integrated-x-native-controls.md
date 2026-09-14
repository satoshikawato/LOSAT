# Independent audit: integrated native TLOSATX controls

Read-only reviewer: `/root/audit_megablast` (`ncbi_parity_auditor`).
Disposition: **FAIL — four of six native time nonregression conditions fail**.
The existing records are retained; no sample is excluded or corrected.

| Fixture | Threads | B1 median s | Integrated median s | Reduction | Time gate |
|---|---:|---:|---:|---:|---|
| MelaMJNV / PemoMJNVA | 1 | 3.750522890 | 3.877465233 | -3.38466% | PASS |
| MelaMJNV / PemoMJNVA | 8 | 2.191553932 | 2.297203809 | -4.82077% | PASS |
| MjeNMV / MelaMJNV | 1 | 16.985607133 | 18.459421946 | -8.67685% | FAIL |
| MjeNMV / MelaMJNV | 8 | 10.206621052 | 11.049535520 | -8.25851% | FAIL |
| AP027280 / AP027280 | 1 | 30.007486731 | 32.648401755 | -8.80085% | FAIL |
| AP027280 / AP027280 | 8 | 18.298717971 | 19.813536844 | -8.27828% | FAIL |

All four failed conditions have disjoint five-sample baseline/candidate ranges.
CPU user+system medians also rise: Mje n1 17.37 to18.89 s, n8 18.05 to19.65 s;
AP n1 30.80 to33.46 s, n8 32.13 to34.95 s. Mela n8 has only3.928 ms margin.
All six RSS comparisons pass the declared max(10%,16MiB) allowance.

The reviewer read actual output bytes for both this group and the older
C-X2-native-controls group: each87 records/151,222,240 bytes, all matching the
corresponding official NCBI BLAST+2.17.0+ oracle. Each group contains3 oracle,
12 diagnostic,12 warmup and60 measured searches. All use code1 local subjects;
the non-default subject-genetic-code exception does not apply.

Command/result/usage records, artifact hashes, ordered argv, source identities,
AB/BA order, non-overlapping monotonic intervals and native pool sizes were
checked. Old/new argv differ only in executable/output paths; these are direct
native executions and do not use Node or a JS runner.

B1 native SHA256: ae494b43e50ad15ae004d23ded13b325cfe8c3bc3b036928ca3302d2d2adf922.
C-X2 native SHA256: 994837f553163851bb812cdaebbf00851e28d001c0d721decf762f1b38a82b5d.
Integrated native SHA256: f6e273760561a034ede0edb1a938b8fb63f619c07a40a613ea9a629438aa56aa.

Root and snapshot589-file bindings and both candidate154 production-source
bindings match. C-X2/integrated native Cargo bin/lib fingerprints match:
Rust1.92.0, LLVM21.1.3, default/parallel/rayon, empty rustflags, release opt3,
LTO, codegen-units1, panic abort. Cargo.lock/build.rs are identical; manifest
changes are comments and config changes are Wasm aliases. Production source
differences are N4 and P1; X2 differs only in a test reference comment.
NCBI link_hsps.c:990-994,1080-1085 and linking.rs:790-811 correspond to the
unchanged reviewed X2 final sort/replay. Parallel group order remains at717-723.
This establishes a regression of the integrated native artifact; it does not
identify X2, N4, P1 or code generation as the cause.

There are12 realtime clock disagreements,11 measured; the largest is
+1.777796547 s. Every accepted wall value equals its stored monotonic endpoint
difference; BOOTTIME disagreement is at most40.657 microseconds. No correction
or exclusion is warranted. Continuous CPU frequency/core-placement data were
not collected. This limitation does not invalidate the retained regression.

The separate continuation stopped with x-native-controls exit1. The original
seven-step ledger and N/P reuse failure remain unchanged; eight later groups
are unexecuted. Prior standalone X2 PASS and the pending memory-policy question
do not waive this integrated native time failure.
