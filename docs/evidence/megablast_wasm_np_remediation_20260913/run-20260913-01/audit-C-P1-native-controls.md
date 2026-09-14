# Independent C-P1 native control audit

Reviewer: `ncbi_parity_auditor` (`/root/audit_megablast`), read-only, 2026-09-14. This audit supports the declared native nonregression controls; it does not support a native speedup claim.

All58 actual output files (2 oracle,8 diagnostic,48 cold) equal their current official oracle and prior C-P1 primary outputs. Four configurations each retain1 warmup and5 measured baseline/candidate pairs:40 measured samples. All24 pair orders, reversed configuration order on odd repetitions and recorded process non-overlap agree. Minimum inter-process gap is0.031778s; main→native group gap is3.094371s. Exclusions are empty and run status is complete. Recorded argv, environment, sinks, result JSON and GNU time entries agree; diagnostic n1/n8 thread counts agree.

| Input | Threads | Baseline median s | Candidate median s | Improvement | Baseline max RSS bytes | Candidate max RSS bytes |
|---|---:|---:|---:|---:|---:|---:|
| AP027078/AP027131 |1|30.494676|31.885836|−4.56198%|68149248|68222976|
| AP027078/AP027131 |8|5.719275|5.605488|1.98953%|100089856|99151872|
| AP027132/NZ_CP006932 |1|48.296780|49.951307|−3.42575%|91078656|91074560|
| AP027132/NZ_CP006932 |8|8.562051|8.818120|−2.99075%|106692608|107753472|

All ranges and maximum RSS values match saved summaries. All conditions satisfy time growth≤max(5% baseline median,50ms) and RSS growth≤max(10% baseline maximum,16MiB). Both n1 medians are slower; AP027078 n1 has only0.134s of remaining tolerance. These results must not be described as native improvement.

The154-file candidate snapshot, build inputs, binaries,4 fixtures,49 saved harness sources,5 shared JS runners and actual oracle identities match metadata and the prior main audit. Candidate `gapalign.rs` SHA is `ae3d1cc128aa2049f56c72947f527920b27f2574f0ffb6dd9aff1504a0cc6fbd`; native baseline SHA is `ae494b43e50ad15ae004d23ded13b325cfe8c3bc3b036928ca3302d2d2adf922`, candidate `569e6430ecaa1738889118b8af0f57ded851c91e8f86a12ce1870f0eb1f5b6ad`. The actual BLASTP2.17.0+ oracle SHA is `fdcbaa25cfdee5359231ae178713d2c28272e694614b847c5b2b7407c3d8a08b`.

Frozen certification, reuse/linear memory, additional tasks and final integration are outside this audit. The unchanged traceback source proof remains in `audit-C-P1.md`. No biological exception applies.
