# Independent integrated N/P n2/n4 boundary audit

The required read-only `ncbi_parity_auditor` reviewed all 20 records: four
official oracles and 16 diagnostics. All are non-measured, with zero exclusions;
the cold summary and engineering checks are correctly empty.

All four main inputs, both versions, and n2/n4 have matching actual argv, raw
bytes, hashes, result records, and usage. N explicitly selects `-task blastn`;
all searches use outfmt 6. Each requested pool has exactly the requested workers,
with unique IDs and `spawn_attempt -> spawned -> ready -> exited(code=0)`.
Saved stage records match actual stderr. Exact monotonic endpoints do not
overlap within the group or with the preceding task controls.

All raw results equal the same inputs' already audited native n1/n8 and Wasm
n1/n8 oracle bytes. The auditor checked NCBI 2.17.0+, Node 24.21.0 without extra
flags, the common `/tmp/.../B1-artifacts` runner, actual B1 `20d604bb...` and
integrated `5b4f4283...` Wasm bytes, 154 frozen source hashes, Cargo inputs,
49 runner hashes, and eight fixture hashes. No genetic-code exception applies.

Two non-measured diagnostics retain adjustable-clock disagreement:

| Relative directory | Realtime minus monotonic |
|---|---:|
| diagnostic/AP027078.AP027131.losatp/baseline-threaded-n4 | +567.770999 ms |
| diagnostic/AP027078.AP027131.losatp/candidate-threaded-n4 | +580.818863 ms |

All 20 records pass monotonic/boottime agreement; maximum difference is
6.209 microseconds. Nothing was excluded or corrected. This audit supports
correctness and worker lifecycle, not a speed claim.

For subsequent reuse evaluation, the final-three plateau applies only to the
same reactor instance. Command's `tail_three_memory_stable=true` is a
not-applicable sentinel because each command job creates a new instance.
