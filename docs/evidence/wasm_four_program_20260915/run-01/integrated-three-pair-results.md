# Integrated production cold results: first three pairs

C0 versus integrated; one warmup, then exactly repeats 1–3 at every condition. Default Node 26.8.2, affinity CPUs 0–7. Positive change means faster. The user stopped measurement during repeat 4; later records remain preserved and are excluded by count.

| Case | Runtime | n | C0 median (s) | Integrated median (s) | Faster (%) | Time guard | RSS guard |
|---|---|---:|---:|---:|---:|---|---|
| AP027078.AP027131.losatp | native | 1 | 29.694982 | 29.349279 | +1.16 | PASS | PASS |
| AP027078.AP027131.losatp | native | 8 | 6.595026 | 5.994218 | +9.11 | PASS | PASS |
| AP027078.AP027131.losatp | serial | 1 | 31.158873 | 30.312198 | +2.72 | PASS | PASS |
| AP027078.AP027131.losatp | threaded | 1 | 30.715644 | 30.479851 | +0.77 | PASS | PASS |
| AP027078.AP027131.losatp | threaded | 8 | 7.947678 | 8.532982 | -7.36 | FAIL | PASS |
| AP027132.NZ_CP006932.losatp | native | 1 | 43.934452 | 45.826301 | -4.31 | PASS | PASS |
| AP027132.NZ_CP006932.losatp | native | 8 | 8.632675 | 8.996885 | -4.22 | PASS | PASS |
| AP027132.NZ_CP006932.losatp | serial | 1 | 45.836657 | 44.380952 | +3.18 | PASS | PASS |
| AP027132.NZ_CP006932.losatp | threaded | 1 | 47.329470 | 45.466225 | +3.94 | PASS | PASS |
| AP027132.NZ_CP006932.losatp | threaded | 8 | 10.915194 | 10.779925 | +1.24 | PASS | PASS |
| AP027280.AP027280.tlosatx | native | 1 | 29.727952 | 26.451539 | +11.02 | PASS | PASS |
| AP027280.AP027280.tlosatx | native | 8 | 18.419366 | 16.062567 | +12.80 | PASS | PASS |
| AP027280.AP027280.tlosatx | serial | 1 | 43.917071 | 31.927045 | +27.30 | PASS | PASS |
| AP027280.AP027280.tlosatx | threaded | 1 | 44.263204 | 39.361184 | +11.07 | PASS | PASS |
| AP027280.AP027280.tlosatx | threaded | 8 | 74.281582 | 63.584923 | +14.40 | PASS | PASS |
| EDL933.Sakai.losatn.megablast | native | 1 | 1.562946 | 1.401015 | +10.36 | PASS | PASS |
| EDL933.Sakai.losatn.megablast | native | 8 | 1.446128 | 1.452088 | -0.41 | PASS | PASS |
| EDL933.Sakai.losatn.megablast | serial | 1 | 2.178150 | 1.914365 | +12.11 | PASS | PASS |
| EDL933.Sakai.losatn.megablast | threaded | 1 | 2.419483 | 2.046134 | +15.43 | PASS | PASS |
| EDL933.Sakai.losatn.megablast | threaded | 8 | 3.014250 | 3.071963 | -1.91 | PASS | PASS |
| MelaMJNV.PemoMJNVA.losatn.blastn | native | 1 | 0.477898 | 0.488086 | -2.13 | PASS | PASS |
| MelaMJNV.PemoMJNVA.losatn.blastn | native | 8 | 0.460249 | 0.401463 | +12.77 | PASS | PASS |
| MelaMJNV.PemoMJNVA.losatn.blastn | serial | 1 | 1.549531 | 1.668804 | -7.70 | FAIL | PASS |
| MelaMJNV.PemoMJNVA.losatn.blastn | threaded | 1 | 0.932750 | 0.867999 | +6.94 | PASS | PASS |
| MelaMJNV.PemoMJNVA.losatn.blastn | threaded | 8 | 1.560742 | 1.657703 | -6.21 | FAIL | PASS |
| MjPMNV.MlPMNV.losatn.blastn | native | 1 | 8.841056 | 8.921009 | -0.90 | PASS | PASS |
| MjPMNV.MlPMNV.losatn.blastn | native | 8 | 4.665271 | 5.178615 | -11.00 | FAIL | PASS |
| MjPMNV.MlPMNV.losatn.blastn | serial | 1 | 37.617360 | 37.458220 | +0.42 | PASS | PASS |
| MjPMNV.MlPMNV.losatn.blastn | threaded | 1 | 13.695291 | 13.632410 | +0.46 | PASS | PASS |
| MjPMNV.MlPMNV.losatn.blastn | threaded | 8 | 7.529271 | 7.755592 | -3.01 | PASS | PASS |
| MjeNMV.MelaMJNV.tlosatx | native | 1 | 16.747868 | 14.441021 | +13.77 | PASS | PASS |
| MjeNMV.MelaMJNV.tlosatx | native | 8 | 10.241742 | 8.645299 | +15.59 | PASS | PASS |
| MjeNMV.MelaMJNV.tlosatx | serial | 1 | 25.525346 | 19.824949 | +22.33 | PASS | PASS |
| MjeNMV.MelaMJNV.tlosatx | threaded | 1 | 25.500857 | 22.015548 | +13.67 | PASS | PASS |
| MjeNMV.MelaMJNV.tlosatx | threaded | 8 | 42.591743 | 35.640448 | +16.32 | PASS | PASS |
| NZ_CP006932.NZ_CP006932.losatn.megablast | native | 1 | 0.251358 | 0.414610 | -64.95 | FAIL | PASS |
| NZ_CP006932.NZ_CP006932.losatn.megablast | native | 8 | 0.277749 | 0.231909 | +16.50 | PASS | PASS |
| NZ_CP006932.NZ_CP006932.losatn.megablast | serial | 1 | 0.673436 | 0.581243 | +13.69 | PASS | PASS |
| NZ_CP006932.NZ_CP006932.losatn.megablast | threaded | 1 | 0.595579 | 0.587433 | +1.37 | PASS | PASS |
| NZ_CP006932.NZ_CP006932.losatn.megablast | threaded | 8 | 1.592003 | 1.566445 | +1.61 | PASS | PASS |
| single-query-32-matches | native | 1 | 0.031462 | 0.032015 | -1.76 | PASS | PASS |
| single-query-32-matches | native | 8 | 0.012920 | 0.012905 | +0.12 | PASS | PASS |
| single-query-32-matches | serial | 1 | 0.139560 | 0.159847 | -14.54 | PASS | PASS |
| single-query-32-matches | threaded | 1 | 0.187019 | 0.217974 | -16.55 | PASS | PASS |
| single-query-32-matches | threaded | 8 | 0.745989 | 0.848402 | -13.73 | FAIL | PASS |

Time guard: regression ≤ max(5%, 50 ms). RSS guard: increase ≤ max(10%, 16 MiB). All 270 retained timed outputs were rehashed and equal their recorded raw hashes; each condition has one shared C0/integrated output hash. All samples and same-version ratios are in the JSON.

## Same-version Wasm / native cold ratios

| Case | Runtime | n | C0 | Integrated |
|---|---|---:|---:|---:|
| AP027078.AP027131.losatp | serial | 1 | 1.049 | 1.033 |
| AP027078.AP027131.losatp | threaded | 1 | 1.034 | 1.039 |
| AP027078.AP027131.losatp | threaded | 8 | 1.205 | 1.424 |
| AP027132.NZ_CP006932.losatp | serial | 1 | 1.043 | 0.968 |
| AP027132.NZ_CP006932.losatp | threaded | 1 | 1.077 | 0.992 |
| AP027132.NZ_CP006932.losatp | threaded | 8 | 1.264 | 1.198 |
| AP027280.AP027280.tlosatx | serial | 1 | 1.477 | 1.207 |
| AP027280.AP027280.tlosatx | threaded | 1 | 1.489 | 1.488 |
| AP027280.AP027280.tlosatx | threaded | 8 | 4.033 | 3.959 |
| EDL933.Sakai.losatn.megablast | serial | 1 | 1.394 | 1.366 |
| EDL933.Sakai.losatn.megablast | threaded | 1 | 1.548 | 1.460 |
| EDL933.Sakai.losatn.megablast | threaded | 8 | 2.084 | 2.116 |
| MelaMJNV.PemoMJNVA.losatn.blastn | serial | 1 | 3.242 | 3.419 |
| MelaMJNV.PemoMJNVA.losatn.blastn | threaded | 1 | 1.952 | 1.778 |
| MelaMJNV.PemoMJNVA.losatn.blastn | threaded | 8 | 3.391 | 4.129 |
| MjPMNV.MlPMNV.losatn.blastn | serial | 1 | 4.255 | 4.199 |
| MjPMNV.MlPMNV.losatn.blastn | threaded | 1 | 1.549 | 1.528 |
| MjPMNV.MlPMNV.losatn.blastn | threaded | 8 | 1.614 | 1.498 |
| MjeNMV.MelaMJNV.tlosatx | serial | 1 | 1.524 | 1.373 |
| MjeNMV.MelaMJNV.tlosatx | threaded | 1 | 1.523 | 1.525 |
| MjeNMV.MelaMJNV.tlosatx | threaded | 8 | 4.159 | 4.123 |
| NZ_CP006932.NZ_CP006932.losatn.megablast | serial | 1 | 2.679 | 1.402 |
| NZ_CP006932.NZ_CP006932.losatn.megablast | threaded | 1 | 2.369 | 1.417 |
| NZ_CP006932.NZ_CP006932.losatn.megablast | threaded | 8 | 5.732 | 6.755 |
| single-query-32-matches | serial | 1 | 4.436 | 4.993 |
| single-query-32-matches | threaded | 1 | 5.944 | 6.809 |
| single-query-32-matches | threaded | 8 | 57.738 | 65.742 |
