# I3: limit the real call boundary to threaded WASI

I2b's serial Mje screen numerically failed: body medians 13.807070 →
15.123494 s (+9.53%). Its final candidate call overlaps an external benchmark
for 10.254939 s, so the complete study is invalid for performance acceptance.
All samples, the numerical failure and overlap remain in `I2b-serial-controls`.
Earlier pairs also show slower candidate times but are not a replacement study.
AP was not reached. No I2b serial regression is certified from this study.

I2b forces a real scan call on every Wasm target, although the diagnosed tier
problem is in threaded WASI. I3 changes only the two inlining attributes and
their explanatory comment: `inline(never)` on the existing threaded-WASI cfg,
`inline(always)` otherwise. Selection, prefix/indexed traversal, all state and
FP statements, candidate set, trace and scheduling stay unchanged. No macro
refactor or algorithm duplication is introduced. Native behavior of these
attributes is unchanged; compiler output must still be checked.

Before timing, inspect actual serial emitted code to verify the helper call is
absent, and run the focused raw/state gates on the new artifact. First fixed
performance screen is Mje serial n1: one warmup plus three AB/BA/AB pairs under
ordinary Node. The original time and memory guards remain. Stop on failure;
do not extend samples. If it passes, complete AP serial, current threaded
n8/n1, Native, reuse and browser conditions. Previous I2b measurements remain
bound to I2b unless rebuilt binary equality is actually established.

NCBI authority: `c++/src/algo/blast/core/link_hsps.c:827–863`. Inlining is a
Rust/compiler implementation choice; no NCBI search behavior changes.

## Decision

Rejected after the fixed Mje serial screen: body 13.540885 → 16.847026 s
(+24.416%); process 13.637316 → 16.979467 s (+24.507%). RSS passes.
All eight raw outputs match. Pre/during/post process snapshots showed no
competing benchmark/build/test; this is sampled observation, not continuous
monitoring. No samples were added and no subsequent condition was run.
