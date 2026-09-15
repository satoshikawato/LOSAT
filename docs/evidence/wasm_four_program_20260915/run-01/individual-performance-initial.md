# Initial individual cold measurements

Node26.8.2 default flags, n8, CPU0–7; one warmup plus five alternating pairs.
These are initial single-candidate results, not final integration acceptance.
Positive percentages mean shorter elapsed time. All shown time/RSS guards pass.

| Candidate | Case | Native seconds A → B | Threaded Wasm seconds A → B | Wasm reduction |
|---|---|---:|---:|---:|
| M1 | EDL933.Sakai.losatn.megablast | 0.970 → 0.935 | 1.786 → 1.722 | +3.63% |
| M1 | NZ_CP006932.NZ_CP006932.losatn.megablast | 0.239 → 0.221 | 0.709 → 0.693 | +2.22% |
| M1 | Sakai.MG1655.losatn.megablast | 0.940 → 0.882 | 1.715 → 1.656 | +3.45% |
| N1 | MelaMJNV.PemoMJNVA.losatn.blastn | 0.315 → 0.306 | 0.764 → 0.754 | +1.36% |
| N1 | MjPMNV.MlPMNV.losatn.blastn | 4.197 → 4.269 | 5.520 → 5.734 | -3.87% |
| N2 | MelaMJNV.PemoMJNVA.losatn.blastn | 0.302 → 0.306 | 0.754 → 0.764 | -1.41% |
| N2 | MjPMNV.MlPMNV.losatn.blastn | 4.065 → 3.967 | 5.426 → 5.369 | +1.05% |
| N3 | MelaMJNV.PemoMJNVA.losatn.blastn | 0.313 → 0.309 | 0.763 → 0.741 | +2.94% |
| N3 | MjPMNV.MlPMNV.losatn.blastn | 4.033 → 4.026 | 5.521 → 5.396 | +2.27% |
| N3 | NZ_CP006932.NZ_CP006932.losatn.megablast | 0.237 → 0.234 | 0.725 → 0.723 | +0.24% |
| P1 | AP027078.AP027131.losatp | 5.153 → 5.295 | 6.293 → 6.329 | -0.57% |
| P1 | AP027132.NZ_CP006932.losatp | 8.015 → 7.726 | 9.044 → 8.700 | +3.80% |
| P1 | SicyWSV.CoBV.losatp | 0.205 → 0.198 | 0.683 → 0.679 | +0.69% |
| P2 | single-query-32-matches | 0.011 → 0.011 | 0.459 → 0.436 | +5.13% |
| P2 | AP027078.AP027131.losatp | 4.981 → 4.847 | 5.896 → 5.652 | +4.14% |

| X2 | MjeNMV.MelaMJNV.tlosatx | 9.003 → 8.034 | 36.891 → 32.292 | +12.47% |
| X2 | AP027280.AP027280.tlosatx | 16.547 → 14.316 | 63.898 → 61.603 | +3.59% |
| X2 | MelaMJNV.PemoMJNVA.tlosatx | 1.878 → 1.805 | 6.320 → 5.582 | +11.69% |

N1/N2 speed effects are inconclusive. N3/M1 show small positive Wasm medians below the5% goal. P1 has mixed effects. P2 reaches a5.13% initial median reduction on its predeclared single-query fixture; this does not certify integration, reuse or memory.
M1 uses the corrected M0 baseline; other candidates use the frozen input baseline.
The AP027078 multi-query input is an out-of-scope control for P2: its observed timing change cannot be attributed to single-query workspace reuse.

Full samples and per-condition metadata are retained in the corresponding temporary performance directories and will be archived with the completed run.
