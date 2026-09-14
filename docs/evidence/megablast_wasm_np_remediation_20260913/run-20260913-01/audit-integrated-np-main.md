# Independent integrated BLASTN/BLASTP main audit

The required read-only `ncbi_parity_auditor` supports the completed fresh threaded-n8 main comparison. All 60 records and 40 measured observations were independently checked against actual raw output bytes. Medians, ranges, maximum process RSS, AB/BA order and non-overlap from exact monotonic endpoints were independently recalculated.

| Input | B1 median (s) | Integrated median (s) | Reduction |
|---|---:|---:|---:|
| BLASTN LC738874/LC738870 | 0.802802 | 0.803456 | −0.08% |
| BLASTN AP027202/LC738875 | 11.574859 | 6.207426 | 46.37% |
| BLASTP AP027078/AP027131 | 7.449706 | 6.618647 | 11.16% |
| BLASTP AP027132/NZ_CP006932 | 10.786954 | 9.556909 | 11.40% |

All raw outputs match their official oracle. Every time nonregression and RSS condition passes, as do the required 5% floors for the large BLASTN input and both BLASTP inputs. Only the explicitly user-approved short-LC speed-floor waiver is applied; no other gate is relaxed.

The auditor verified the frozen source, actual integrated and B1 artifacts, Node 24.21.0, NCBI 2.17.0+, identical actual runner paths/hashes and absence of TBLASTX-only Node flags in NP runs. Integrated threaded command SHA-256 starts `5b4f4283`; full hashes remain in group metadata and artifact records.

This audit accepts the main comparison scope, not the still-pending controls, measured reuse or scratch-capacity claims. Earlier component medians and failed/excluded cohorts remain separate.
