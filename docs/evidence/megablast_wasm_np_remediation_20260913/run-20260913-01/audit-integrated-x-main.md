# Independent integrated TBLASTX main audit

The required read-only `ncbi_parity_auditor` supports the completed nonregression and raw-output scope. All 45 actual raw records (three oracle, six diagnostic, six warmup and 30 measured) were independently verified. Saved argv, usage/result records, 154 source hashes, Cargo configuration, Node/runner/artifact identities, AB/BA order and non-overlap at exact monotonic endpoints match.

| Input | B1 median (s) | Integrated median (s) | Reduction | B1 maximum RSS (bytes) | Integrated maximum RSS (bytes) |
|---|---:|---:|---:|---:|---:|
| MelaMJNV/PemoMJNVA | 4.120223462 | 4.130871705 | −0.25844% | 323469312 | 331804672 |
| MjeNMV/MelaMJNV | 13.576848569 | 13.923771091 | −2.55525% | 346382336 | 335249408 |
| AP027280 self | 25.556511331 | 25.562424232 | −0.02314% | 350535680 | 358797312 |

All three recalculated time and RSS gates pass. All three raw-output hashes match the previous official oracle bytes. The result is limited to the predeclared 4N final-sort copy-removal contract; it does not demonstrate whole-search speedup.

The reviewer also examined the unexecuted summary-table script: 37 cold conditions, 28 reused-invocation conditions and 16 separate preparation records, with first warmup and subsequent measured medians kept separate. A missing uniqueness/completeness assertion was identified in that aggregator, not in actual measured data. The script was strengthened before execution to compare exact cold check/summary key sets and exact reused invocation/observation key sets, rejecting duplicates before dictionary overwrites could hide a missing condition. Measurement scripts, samples and frozen production source were unchanged.

Follow-up static review confirms the aggregator correction resolves the duplicate/missing-key issue. Syntax parsing passes; the aggregator was not executed during review. Reviewed script SHA-256: `6d915083af012380fa710ea534d92b9fdcff4852bcc15c38412d65f03c94bcb6`.
