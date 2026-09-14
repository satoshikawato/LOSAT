# Independent native CPU sampling audit

Reviewer: /root/audit_megablast (ncbi_parity_auditor), read-only.
Disposition: COMPLETE_DIAGNOSTIC_ONLY; native time gate FAIL unchanged.

All1192 bound file hashes, five binary/source/fixture bindings, argv and exact
outputs pass. Five actual outputs total8,643,295 bytes, each SHA256
b309f77fee038f0559d870737ebeba732c11eba4f7e495f045acbfb831817fd7.
The fixed order was B1,N4,P1,X2,integrated, one run each, all timed=false.
All process intervals are disjoint and record/report commands exit0.

The reviewer parsed actual perf.data: software CPU-clock, user-space only,
frequency199, no callchain/stack dump, no lost/throttle/unthrottle records.
Each profile has one PID/TID. All16–27 unresolved samples per file are vdso.
No IPC or hardware-counter evidence exists.

| Artifact | All samples | Link function self samples | Percentage |
|---|---:|---:|---:|
| B1 | 3282 | 2909 | 88.63% |
| C-N4 | 3638 | 3243 | 89.14% |
| C-P1 | 3617 | 3234 | 89.41% |
| C-X2 | 3298 | 2920 | 88.54% |
| integrated | 3739 | 3344 | 89.44% |

The link_hsp_group_ncbi function dominates this fixture's CPU samples.
The fixed-order single runs cannot establish individual N4/P1 causality,
regression rates, instruction/cache causes or caller-specific costs.
Inline work is included in the self-symbol aggregation.
The stderr sum_stats_linking timer at run_impl.rs2450–2498 measures only
postliminary linking; the preliminary call at2305 is outside it. Its roughly
1.1-second value therefore does not contradict the sampled link hotspot.
