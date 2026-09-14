# Independent C-P1 review

Reviewer: `ncbi_parity_auditor` (`/root/audit_megablast`), read-only, 2026-09-13.

The final `gapalign.rs` SHA256 is `ae3d1cc128aa2049f56c72947f527920b27f2574f0ffb6dd9aff1504a0cc6fbd`. Source, patch and identity agree; other production/build inputs match B1. NCBI `blast_gapalign.c:563–575,634–639,666–669` supports the same fence decision, cell recurrence, write order and final used-row span. The row allocator clears length, actual-band reserve guarantees capacity, DP and trace slices belong to separate vectors, and every visited non-fence cell initializes one byte, including X-drop failures. The only early break saves the fence index before the write. `set_len` exposes exactly the initialized prefix after the borrow ends. Unwinding before that point leaves length zero; no uninitialized byte is exposed. Existing band-growth and final row resize remain unchanged.

The reviewer checked the added forward/reverse fence tests at the first, middle and last residue, including the reused capacity8/len2=8/gap_extend=0 boundary where the old allocator's relative reserve alone cannot cover the nine-cell band. Focused final release tests report20 passed,1 ignored; the ignored diagnostic is not counted as coverage.

The primary five-pair confirmation contains30 actual output files: oracle2, diagnostic4, cold24 including warmup. Every raw output matches the official NCBI2.17.0+ oracle, with no exclusions. Warmup exclusion, alternating order, process non-overlap, identical runner/Node/search options, source/artifact/input hashes, medians/ranges and maximum RSS were independently recomputed.

| Main input | B1 median [range] s | C-P1 median [range] s | Improvement | Maximum RSS B1→C-P1 bytes |
|---|---:|---:|---:|---:|
| AP027078/AP027131 | 8.235 [8.066–8.270] | 7.267 [7.205–7.440] | 11.76% | 312467456→309039104 |
| AP027132/NZ_CP006932 | 11.848 [11.600–11.957] | 10.372 [10.089–11.123] | 12.45% | 327868416→322727936 |

This supports the primary threaded n8 performance condition only. Native controls were running; threaded1/task controls, formats, reuse/linear memory and final integration were pending. No full adoption claim is made. Serial is outside the revised normal adoption axis and remains explicit compatibility scope.
