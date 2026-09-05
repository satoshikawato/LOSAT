# PR01 measured results

Status: NOT_PROVEN. These are baseline measurements, not optimization or release certification.

| Run | Status | Recorded rows |
| --- | --- | ---: |
| initial-audit | FAIL | 84 |
| cold | FAIL | 639 |
| warm | FAIL | 61 |
| protected-profile | PASS_EXECUTED_SCOPE | 8 |

Initial audit raw Gate A revalidation: {'PASS': 69, 'TIMEOUT': 15}.
Initial audit oracle-contract revalidation: {'PASS': 68, 'TIMEOUT': 14, 'NOT_RUN': 2}.

Cold medians below each require one successful warmup, five successful timed outputs, and the separate diagnostic gate. Missing cells are NOT_RUN. Times are seconds.

| Case | Native 1 | Plain serial Wasm | Threaded Wasm 1 | Serial/native | Requested 4 speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| PesePMNV.MjPMNV.task_blastn | 0.230 | 0.732 | 0.741 | 3.180 | 1.008 |
| EDL933.Sakai.megablast | 1.061 | 1.718 | 1.735 | 1.619 | 0.942 |
| compact.multi_query.outfmt6 | 0.037 | 0.167 | 0.163 | 4.506 | 0.827 |
| compact.multi_query.no_hit.outfmt7 | 0.047 | 0.170 | 0.177 | 3.607 | 1.099 |
| pairwise_default_serial | 0.792 | 1.316 | 1.362 | 1.662 | 1.877 |
| all_hsp_self_serial | 1.046 | 1.725 | 1.756 | 1.649 | 2.032 |
| p03_mela_pemojnva | 4.564 | 8.110 | 8.183 | 1.777 | 1.023 |
| p11_avclpv_psclpv | NOT_RUN | NOT_RUN | NOT_RUN | NOT_RUN | NOT_RUN |
| p12_lc738874_lc738875_default | 4.817 | 10.128 | 10.166 | 2.103 | 0.992 |
| pr01_blastp_no_hit | 0.338 | 0.620 | 0.625 | 1.832 | 1.085 |
| pr01_tblastx_no_hit | 0.046 | 0.189 | 0.185 | 4.138 | 0.823 |

Complete median ranges and RSS: [cold/summary.json](cold/summary.json). All 2/4/8 requested-N comparisons: [comparisons.json](comparisons.json). Warm serial results: [warm/summary.json](warm/summary.json).

[samples.tsv](samples.tsv) retains failures, warmups and diagnostics as separate rows. [stage-profiles.json](stage-profiles.json) links diagnostic output and logs. Smoke runs are excluded from these tables. Stage timers are nested and may sum work across threads; they are not an additive wall-time partition.

No whole-eligible-set geometric mean or performance attainment is accepted while an eligible case is missing. Requested-N speedup is not evidence that N computation threads were active. Platform Gate B, real browsers, resident API, threaded warm lifecycle, and native active-worker measurements remain NOT_RUN.
