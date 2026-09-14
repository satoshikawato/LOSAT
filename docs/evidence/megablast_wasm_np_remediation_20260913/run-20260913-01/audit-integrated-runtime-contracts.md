# Independent integrated runtime contract audit

The required read-only `ncbi_parity_auditor` independently examined saved actual bytes, argv, worker events and artifacts. Result: PASS for the completed command/reactor and threshold scope.

| Command checker record category | Standard threaded | Explicit compatibility |
|---|---:|---:|
| Official oracle records | 23 | 23 |
| Successful raw comparisons | 166 | 212 |
| Expected failures | 57 | 124 |
| stdout cases | 6 | 9 |
| API harness records | 4 | 5 |
| Total | 256 | 373 |

Every successful raw comparison matches the actual saved oracle bytes and recorded argv; the expected negative-condition sets match exactly. These totals are not all NCBI comparison counts.

Each reactor suite contains 144 records (108 successful, 36 expected failures). The auditor checked every record, including the 1→2→4→8→2→1 sequence, n2 recovery after spawn rejection, format bytes, and clearing results on invalid input. Each host has 294 spawn attempts with no thread-ID reuse; recorded events exactly match actual stderr events. In the 24-repeat exercise, the final 12 linear-memory observations plateau at 12,976,128 bytes for the default build and 12,910,592 bytes for compatibility. This is neither scratch-capacity nor speed evidence.

The capacity-limit suites (36 each plus 12 serial records) have the expected statuses, empty results on rejection and exact recovery bytes. The threshold suites (nine default and 12 compatibility records) pass actual raw and worker checks. Actual artifact separation and all 13 successful gate exits are consistent.

API formatter checks use their recorded API n1/reference and default/compatibility byte contracts. No CLI-oracle label normalization is performed. This audit does not accept pending integrated performance, measured reuse or diagnostic scratch-capacity claims.
