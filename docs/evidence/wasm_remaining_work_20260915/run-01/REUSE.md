# Bounded serial reuse verification

Same artifacts, ordinary Node, two independent AB/BA sessions. Each case has one warmup and three measured calls per version/session. All192 invocations returned exact fresh-oracle bytes. The time and RSS guards pass **in each session independently**; linear-memory sequences are identical. This does not resolve existing long-window reactor-memory issues.

| Artifact / mode | Case | Baseline seconds, session0 /1 medians | H1 seconds, session0 /1 medians | Peak RSS baseline / H1 (MiB) | Memory / guards |
|---|---|---:|---:|---:|---|
| serial-command / compiled-module | valid-nohit | 0.005643 / 0.005006 | 0.006197 / 0.005376 | 83.41 / 78.80 | equal / PASS |
| serial-command / compiled-module | single32 | 0.034201 / 0.030948 | 0.034307 / 0.032021 | 108.31 / 116.38 | equal / PASS |
| serial-command / compiled-module | word7 | 0.033383 / 0.030206 | 0.033432 / 0.031307 | 126.16 / 134.02 | equal / PASS |
| serial-command / compiled-module | EDL933.Sakai.losatn.megablast | 1.115191 / 1.073575 | 1.079538 / 1.126610 | 1025.97 / 1041.93 | equal / PASS |
| serial-reactor / compiled-module | valid-nohit | 0.005175 / 0.005329 | 0.005200 / 0.005309 | 95.96 / 97.57 | equal / PASS |
| serial-reactor / compiled-module | single32 | 0.029523 / 0.029428 | 0.030364 / 0.031008 | 121.11 / 125.74 | equal / PASS |
| serial-reactor / compiled-module | word7 | 0.032434 / 0.033341 | 0.033445 / 0.034155 | 127.22 / 129.62 | equal / PASS |
| serial-reactor / compiled-module | EDL933.Sakai.losatn.megablast | 1.140727 / 1.109751 | 1.126597 / 1.119825 | 1116.32 / 1114.45 | equal / PASS |
| serial-reactor / same-instance | valid-nohit | 0.003213 / 0.003110 | 0.003189 / 0.003165 | 89.12 / 89.52 | equal / PASS |
| serial-reactor / same-instance | single32 | 0.028625 / 0.027795 | 0.028120 / 0.027159 | 112.79 / 113.52 | equal / PASS |
| serial-reactor / same-instance | word7 | 0.031484 / 0.031180 | 0.031981 / 0.030801 | 119.51 / 114.54 | equal / PASS |
| serial-reactor / same-instance | EDL933.Sakai.losatn.megablast | 1.049385 / 1.025191 | 1.070519 / 1.027094 | 404.03 / 398.58 | equal / PASS |

The measured boundary includes instance preparation (when fresh), input copy, search, result copy and input release. Compiled-module samples also include per-job close. Same-instance final per-case close is outside invocation samples and is recorded separately as `timings.final_instance_close_seconds`; it is not included in the setup table below. Serial calls create no workers. Same-instance mode uses one instance per case. Source input text is loaded outside invocation timing; its resident bytes are recorded. Compile/validation setup is separate:

| Artifact / mode / session | Baseline setup (ms) | H1 setup (ms) |
|---|---:|---:|
| serial-command / compiled-module / 0 | 12.936 | 4.913 |
| serial-command / compiled-module / 1 | 10.873 | 4.615 |
| serial-reactor / compiled-module / 0 | 9.116 | 3.974 |
| serial-reactor / compiled-module / 1 | 9.770 | 4.101 |
| serial-reactor / same-instance / 0 | 10.273 | 3.670 |
| serial-reactor / same-instance / 1 | 10.196 | 4.090 |

All raw calls, warmups, setup components, CPU/RSS samples and exact linear-memory sequences are in `measurement.tar.gz` under `H1-runtime/`. RSS is process-lifetime peak, and different cases can contribute to it; it is not an allocation counter. No steady-state speedup or unlimited-repeat memory stability is inferred.
