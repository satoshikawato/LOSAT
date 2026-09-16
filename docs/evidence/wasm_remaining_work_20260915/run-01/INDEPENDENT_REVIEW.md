# Independent read-only review

Reviewer: `ncbi_parity_auditor`, 2026-09-15.

The reviewer inspected source, records and actual output files. It did not run
new tests or benchmarks during exclusive measurements.

Final reconciliation found no material H1 implementation defect:

- Five applied host-file hashes match the tested candidate.
- All 158 Rust/build/configuration inputs match baseline.
- All 15 final fixture outputs match actual oracle files and recorded hashes.
- All 192 reuse output files match actual oracle hashes.
- All 24 session/case timing and RSS guards pass after independent recomputation.
- Four artifact-kind inspector JSON comparisons and serial reactor recovery
  output/error files match between old and new hosts.

The reviewer supports rejecting X1 because all four default body medians
regress. H1 adoption is supported only as a measured Node/WASI serial preparation
improvement: 7.540686 to 4.699592 ms. Its primary complete-process difference is
51.205032 to 50.398820 ms; this is not a 5% overall improvement claim.

The requested documentation corrections were applied: retained threaded timings
are labelled baseline references, same-instance final close is distinguished
from per-job close, and STATUS identifies the body boundary clocks used for
measurements. The initial stronger compiler-tier attribution was narrowed, then
updated with the symmetric TurboFan diagnostic evidence.

Archive integrity is verified separately by the primary agent in
`archive-verification.json`; it is not attributed to this reviewer. This review
is not full-input, platform-release or browser certification and does not resolve
the existing reactor-memory or invalid-TBLASTX-query status issues.
