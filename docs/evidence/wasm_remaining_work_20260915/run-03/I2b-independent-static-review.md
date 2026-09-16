# I2b independent read-only review

Recorded by the parent from `ncbi_parity_auditor` responses. The auditor
modified no files and ran no performance measurements.

NCBI authority: `c++/src/algo/blast/core/link_hsps.c:675-684,827-863,876-885`.
The I2 prefix chooses original index `length+1`; assigning its tail is the
original decrement. On rejection, `next_larger.saturating_sub(1)` retains
exactly the intended backward range. Targets 0/1 end the scan. The no-visit
guard preserves an empty helper slice for h_idx=2. Macro hygiene and continue
targets are consistent. Restoring the two cursor statements produces the same
249 body tokens as I1; all helper-external source also matches.

The first review found that I2's atomics selector is false on this Rust target.
I2b corrects only the two cfg attributes to the project's existing
`all(losat_wasi_threads, feature="wasm-threads")` and its negation. The auditor
independently checked this exact I2→I2b relation. No new static blocker was found.

The linking hash `fe3ced55866481fe33f42aee36a880fc0d0237657186768b3ffd5a9f1a3e8da5`
matches source and differential manifests. The auditor checked the harness
hash, real spawn/join, explicit build.rs-equivalent cfg, 802 cases, and all
12 successful execution records. Full stdout hashes agree; trace stderr hashes
agree. The I1 full-group comparison summary was checked, without repeating its
large raw comparisons in this short review. Dynamic I2b full-group state and
performance remain separate acceptance work.

All auditor tools closed at **2026-09-16 01:49:43 UTC**; no background work.

## Emitted-code follow-up

The auditor independently confirmed that helper 1368 retains r8 from
+0x191/+0x195 through the contender branch to coordinate loads at
+0x47f/+0x489/+0x49b. I1's second multiply/add is absent there. Initial
slice, b0 prefix/saturating-subtraction, and payload bounds remain.
The raw artifact hash was recomputed and matches the metadata. This supports
the bounded code hypothesis, not a speed or universal-tier claim. All tools
closed at **2026-09-16 01:51:55 UTC**, with no background work.
