# C-X2 long gencode4 independent review

Recorded from the read-only `ncbi_parity_auditor` review by `/root/audit_megablast` on 2026-09-13 UTC.

All nine records pass actual raw-byte comparison: the fresh official NCBI database oracle and B1/C-X2 native/threaded n1/n8. Every output is 2,270,680 bytes with SHA-256 `1da3bfe7d96d50a96726a44e2df98f5e81c05bee9cd655ad51bb8d44aa377963`.

The reviewer checked the ten database files, makeblastdb executable/input hashes, ordered search argv, 154 C-X2 Rust sources, 49 runner snapshots, and all four measured LOSAT artifacts. At n1 there are no host workers. Both threaded n8 runs have unique IDs 1–8 matched through spawn, readiness and exit code zero; four linking jobs select the parallel path.

These are correctness diagnostics, with zero adoption samples. No gencode4 timing or speedup claim follows from their durations. The local-subject non-default subject genetic code uses the approved database-oracle contract.

The final candidate-only replay design is also independently reviewed: it keeps all four integrated target/thread conditions while reusing this exact saved oracle, with input/tool/output identity checks, full source guards, and Node/live runner/import/artifact hashes before and after each invocation. The updated driver matches the predeclared amendment hash. Its execution remains pending.

The same reviewer completed the LC738874/LC738875 threshold sweep: all nine oracle/native n1/threaded n4 records match at evalue 10/100/10000 (321,589 / 508,077 / 3,582,531 bytes). Input/artifact/Node/runner identities and unique workers through exit zero are verified. The reviewer supports selecting C-X2 for integration under the 4N copy-reduction contract; integrated acceptance remains separate.
