# Integrated performance progress

The frozen integrated source passed the build and correctness gates. The main BLASTN/BLASTP measurement group is complete: one warmup and five measured AB/BA pairs, shared B1 runner paths, unchanged Node settings and exact raw-output admission. Independent review of every sample, source/tool binding and recalculated gate passes. The remaining control/reuse/capacity groups are separate gates.

| Input | Baseline n8 cold median (s) | Integrated median (s) | Reduction | Time/RSS |
|---|---:|---:|---:|---|
| MelaMJNV.PemoMJNVA.losatn.blastn | 0.802802 | 0.803456 | -0.08% | PASS/PASS |
| MjPMNV.MlPMNV.losatn.blastn | 11.574859 | 6.207426 | 46.37% | PASS/PASS |
| AP027078.AP027131.losatp | 7.449706 | 6.618647 | 11.16% | PASS/PASS |
| AP027132.NZ_CP006932.losatp | 10.786954 | 9.556909 | 11.40% | PASS/PASS |

The large BLASTN input and both BLASTP inputs pass the required 5% floor. The short BLASTN input is 0.65 ms slower in this cohort, within the unchanged 50 ms nonregression bound; its speed floor is waived by the explicit user instruction. Earlier component cohorts remain preserved and are not pooled into this fresh integrated result.

## TBLASTX main comparison

The predeclared acceptance scope is final-sort data movement reduction, with exact raw output and unchanged nonregression/RSS gates. There is no TBLASTX speed floor and no whole-search speedup claim. Independent actual-byte/sample/tool-binding review and recalculation of every gate pass.

| Input | Baseline n8 cold median (s) | Integrated median (s) | Reduction | Time/RSS |
|---|---:|---:|---:|---|
| MelaMJNV.PemoMJNVA.tlosatx | 4.120223 | 4.130872 | -0.26% | PASS/PASS |
| MjeNMV.MelaMJNV.tlosatx | 13.576849 | 13.923771 | -2.56% | PASS/PASS |
| AP027280.AP027280.tlosatx | 25.556511 | 25.562424 | -0.02% | PASS/PASS |

## Native BLASTN/BLASTP controls

All eight native conditions pass the declared time and RSS gates. Independent review of all 116 records and actual identities passes. BLASTP is slower in these native medians; these data are not a native speedup claim. The short BLASTN n8 percentage exceeds 5% but its absolute change is below the unchanged 50 ms tolerance. Long native n1 BLASTP has a narrow margin to the 5% bound, retained below and in the full sample ranges.

| Input | Threads | Baseline median (s) | Integrated median (s) | Reduction | Delta (s) | Allowed increase (s) |
|---|---:|---:|---:|---:|---:|---:|
| MelaMJNV.PemoMJNVA.losatn.blastn | 1 | 0.491106 | 0.497652 | -1.33% | 0.006546 | 0.050000 |
| MelaMJNV.PemoMJNVA.losatn.blastn | 8 | 0.350307 | 0.369669 | -5.53% | 0.019363 | 0.050000 |
| MjPMNV.MlPMNV.losatn.blastn | 1 | 11.473342 | 9.075481 | 20.90% | -2.397860 | 0.573667 |
| MjPMNV.MlPMNV.losatn.blastn | 8 | 6.859844 | 4.463187 | 34.94% | -2.396656 | 0.342992 |
| AP027078.AP027131.losatp | 1 | 28.694368 | 30.000106 | -4.55% | 1.305738 | 1.434718 |
| AP027078.AP027131.losatp | 8 | 5.403868 | 5.571826 | -3.11% | 0.167958 | 0.270193 |
| AP027132.NZ_CP006932.losatp | 1 | 43.468891 | 45.382233 | -4.40% | 1.913342 | 2.173445 |
| AP027132.NZ_CP006932.losatp | 8 | 8.358210 | 8.516648 | -1.90% | 0.158438 | 0.417911 |

## Threaded-Wasm n1 BLASTN/BLASTP controls

All four conditions pass raw output, nonregression and RSS gates. Both BLASTP n1 main medians exceed the separate 5% minimum for a single-threaded-Wasm speedup claim; independent review of all 60 actual records and recalculated gates passes. This uses the same threaded artifact at n1, with zero worker spawns, and is separate from plain serial compatibility Wasm.

| Input | Baseline median (s) | Integrated median (s) | Reduction | Time/RSS |
|---|---:|---:|---:|---|
| MelaMJNV.PemoMJNVA.losatn.blastn | 0.819908 | 0.819672 | 0.03% | PASS/PASS |
| MjPMNV.MlPMNV.losatn.blastn | 18.905531 | 13.631000 | 27.90% | PASS/PASS |
| AP027078.AP027131.losatp | 35.381984 | 30.093954 | 14.95% | PASS/PASS |
| AP027132.NZ_CP006932.losatp | 52.543935 | 44.853448 | 14.64% | PASS/PASS |

## Additional BLASTN, megablast and BLASTP controls

All nine threaded-n8 control conditions pass actual raw output, median nonregression and RSS checks. These are control inputs; the 5% main-input minimum is not imposed on them. Independent actual-file review and recalculation of all 135 records pass. Seven adjustable-clock diagnostic disagreements are documented in the audit; monotonic/boottime agree and all samples remain included.

| Input | Baseline median (s) | Integrated median (s) | Reduction | Time/RSS |
|---|---:|---:|---:|---|
| PesePMNV.MjPMNV.losatn.blastn | 0.649265 | 0.654976 | -0.88% | PASS/PASS |
| NZ_CP006932.NZ_CP006932.losatn.blastn | 6.177420 | 6.053328 | 2.01% | PASS/PASS |
| PeseMJNV.PemoMJNVB.losatn.blastn | 1.831472 | 1.780202 | 2.80% | PASS/PASS |
| NZ_CP006932.NZ_CP006932.losatn.megablast | 0.705832 | 0.706675 | -0.12% | PASS/PASS |
| EDL933.Sakai.losatn.megablast | 1.892816 | 1.888755 | 0.21% | PASS/PASS |
| Sakai.MG1655.losatn.megablast | 1.849298 | 1.865824 | -0.89% | PASS/PASS |
| AP027131.NZ_CP006932.losatp | 11.103862 | 9.691967 | 12.72% | PASS/PASS |
| WSSV.PajaWSV.losatp | 1.220243 | 1.170931 | 4.04% | PASS/PASS |
| SicyWSV.CoBV.losatp | 0.724288 | 0.692525 | 4.39% | PASS/PASS |
