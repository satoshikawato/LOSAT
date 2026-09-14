# Independent integrated threaded-Wasm n1 controls

The required read-only `ncbi_parity_auditor` supports all four raw, median nonregression and RSS gates. Both BLASTP inputs support a greater-than-5% improvement specifically for cold threaded command-Wasm at n1.

All 60 records (four official oracle, eight diagnostic, eight warmup, 40 measured) were independently checked. There are zero exclusions and every actual raw byte comparison passes. The five-pair medians, ranges and maximum RSS match the saved calculations.

| Input | B1 median (s) | Integrated median (s) | Reduction | B1 maximum RSS (bytes) | Integrated maximum RSS (bytes) |
|---|---:|---:|---:|---:|---:|
| LC738874/LC738870 | 0.819908 | 0.819672 | 0.03% | 138588160 | 137228288 |
| AP027202/LC738875 | 18.905531 | 13.631000 | 27.90% | 190246912 | 183951360 |
| AP027078/AP027131 | 35.381984 | 30.093954 | 14.95% | 165564416 | 166494208 |
| AP027132/NZ_CP006932 | 52.543935 | 44.853448 | 14.64% | 174620672 | 180764672 |

P078 ranges are 34.865–35.851→29.851–30.613s; P132 ranges are 52.468–53.652→44.539–45.513s. The short LC difference is approximately 0.24 ms and is not treated as a meaningful speedup.

The auditor verified actual blastn/task-blastn and blastp outfmt 6 / n1 argv, eight fixture hashes and official NCBI 2.17.0+ binary hashes. There is no genetic-code exception in this group. Both versions use the same actual B1-artifacts runner path and helper hashes, Node 24.21.0 and no extra Node flags. Actual B1 (`20d604bb…`) and integrated (`5b4f4283…`) Wasm hashes match the metadata, as do the frozen snapshot, HEAD, 154 production-source hashes, Cargo inputs and 49 runner hashes.

Every diagnostic has pool_threads=0, effective_compute_threads=1, nonparallel stages and zero spawn_attempt/spawned/ready/exited events. No worker events occur in stderr. Recorded argv/result/usage/raw hashes, per-repeat AB/BA and exact monotonic endpoints establish non-overlap, including with the preceding native group. Raw bytes also match the already audited native n1/n8 and threaded n8 oracles for these inputs.

The speed boundary includes command startup, validation, guard, compile, search, output and termination. It does not imply native BLASTP improvement or completion of remaining task controls, measured reuse, capacity or formal certification.
