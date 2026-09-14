# Independent C-NA64 N/P native audit

The read-only ncbi_parity_auditor supports raw parity, time nonregression and RSS
PASS for all eight native N/P conditions. Combined with the separately audited
X group, the new native configuration passes all fourteen conditions.

| Input | n | B1 median (s) | Candidate median (s) | Reduction |
|---|---:|---:|---:|---:|
| LC738874/LC738870 | 1 | 0.469533 | 0.471080 | -0.33% |
| LC738874/LC738870 | 8 | 0.356120 | 0.334276 | +6.13% |
| AP027202/LC738875 | 1 | 10.528382 | 8.438493 | +19.85% |
| AP027202/LC738875 | 8 | 6.162598 | 4.088946 | +33.65% |
| AP027078/AP027131 | 1 | 27.095781 | 28.126090 | -3.80% |
| AP027078/AP027131 | 8 | 5.001962 | 5.165810 | -3.28% |
| AP027132/NZ_CP006932 | 1 | 41.578226 | 42.969299 | -3.35% |
| AP027132/NZ_CP006932 | 8 | 7.927434 | 7.992196 | -0.82% |

All 116 records (4 oracles, 16 diagnostics, 16 warmups, 80 measured runs) have
actual raw output equal to the corresponding fresh official oracle: 155,326,784
bytes checked. All commands, results and usage records, five-pair arrays,
medians, ranges and maximum RSS agree. Zero exclusions, prescribed alternating
order and no monotonic overlap; minimum process gap 47.541304 ms, with 7.000680
seconds after the preceding X group.

RSS uses max(10% of baseline maximum, 16 MiB); every condition passes. Native
BLASTP is not faster: AP027078 n1 ranges are disjoint (B1 26.572–27.576 seconds,
candidate 27.774–29.304), although the regression stays within 5%. Short BLASTN
n8 improves by 21.844 ms with overlapping ranges.

Fifteen records (eleven timed) have realtime disagreement. The largest difference
is oracle P132 +2.194024 seconds; the largest timed difference is repeat2 P132
baseline n1 +2.162081 seconds. Maximum absolute boottime difference is 8.193
microseconds. No sample was corrected or excluded.

The auditor verified 589 source paths/hashes at all three source locations,
154 production metadata hashes, actual fixtures, tools and binaries. B1 is
`ae494b43...`; candidate is `a1ab2845...`; new manifest is `08521898...`. These are
NCBI 2.17.0+ native x86_64 outfmt6 n1/n8 comparisons without a parity exception.
The mutable current-fingerprint caveat remains. Wasm byte identity and later
gates require separate evidence; this native audit does not qualify them.
