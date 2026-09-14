# Independent TBLASTX movement-cost audit

Reviewer: `ncbi_parity_auditor` (`/root/audit_megablast`), read-only, 2026-09-14. Recalculation errors:0. This is diagnostic evidence, not candidate adoption or a predicted speedup.

All15 actual outputs match corresponding current NCBI TBLASTX2.17.0+ oracle bytes, genetic codes1/1 and outfmt6. The9 Wasm invocations use the same `--no-liftoff --no-wasm-tier-up` flags and B1 runner path. The155-source instrumentation hash map, artifact `23e019ad…`, other artifacts, fixtures, oracle and runner hashes match metadata.

| Input | Linking-sort logical clone bytes | Aggregate sort seconds |
|---|---:|---:|
| MelaMJNV/PemoMJNVA |20725584|0.021482657|
| MjeNMV/MelaMJNV |48783504|0.054337746|
| AP027280 self |57000240|0.067010365|

These are this Wasm target's88-byte HSP payload times clone occurrences, not unique HSP counts or memory-bus traffic. Each linking call's3 sorts copy2N records each. The recorded6 sort calls across2 linking calls agree with source.

The local linking-sort timers are contained in both global sort timers and link-total timers. Global sort durations are0.040918766/0.085574993/0.102785829s; link totals are1.003705976/10.092151677/22.432198155s. They must not be added or equated with removable wall time. Clone byte categories also overlap.

The initial parser rejected AP027280 because stderr70–72 contains two complete host worker-done lines interleaved between `[X_COST] name=` and `x_sort_original_clone`. Reassembling only that exact transport span reproduces all7 unique expected diagnostic records and the saved fragment list. Other logs need no reassembly. Original stderr and the initial parser remain available; raw search outputs were never normalized. The reviewer independently confirmed this reconstruction.
