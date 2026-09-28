# Independent read-only TBLASTN Stage F audit — 2026-09-25

Verdict: **SUPPORTED / PASS for the stated fixture-scoped native and command-WASI Stage F parity gate**. This review covered the final uncommitted candidate from Stage E pass `911d91f5cfc31eae3be48cf5ca9efac9970474bd` and fixed NCBI source `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`. It does not certify a release or a speedup.

The independent `ncbi_parity_auditor` read the Rust scheduling diff, pinned NCBI source, comparison scripts and refreshed evidence. Two concrete findings were resolved before this verdict:

1. Subject work now finishes preliminary search, link/direct E-value and preliminary reap within each subject job, as `c++/src/algo/blast/core/blast_engine.c:804-905,1469-1475` requires. The source-order collector and later Kappa, containment, heap and hitlist operations remain serial.
2. Parallel work is collected in batches capped at `pool.threads()`, then reduced in subject OID order. This bounds retained intermediate subject HSPs by worker count instead of total subject count.

The auditor found no remaining concrete parity or order defect in this diff. The independent row and hash review found 3,000/3,000 native complete-byte matches (1,794 physical, 396 single-subject genetic-code, 810 eight-subject genetic-code), 564/564 command-WASI matches (102 plain serial, 462 threaded), the full 27-code and 0/6/7 grids, and repeated 2/4/8 runs on all eight-subject code cases. All generated eight-subject input checksums were rehashed; the recorded native/WASI executable hashes match the built binaries. Actual simultaneous subject work peaked at 2/4/8. Code 1 was calibrated to local NCBI CLI bytes; code 32 used the comparison-only selected-code `FindGeneticCode(32)` C++ API oracle. The only accepted subject-code exception is `PD-TLOSAN-LOCAL-GENCODE-32`.

Focused TBLASTN tests passed 106/106, search-pool tests 3/3, Stage D oracle replay 130/130, Stage D checksum validation 615/615, and clippy/fmt passed. The final evidence manifest must pass `sha256sum -c evidence.sha256` after this document and the gate report are finalized.

Stage G retains the existing BLASTN/BLASTP/TBLASTX regression and full certification matrices, formal speed/resource protocol, broader platform/option certification, and final release-facing audit.
