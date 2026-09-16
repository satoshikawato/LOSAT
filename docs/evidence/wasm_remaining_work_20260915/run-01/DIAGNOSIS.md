# Diagnosis and candidate selection

## Current reproducible behavior

Both primary TBLASTX inputs have the same per-group HSP counts, helper visits,
rejection counts, predecessor updates and chain counts at requested n1/n2/n4/n8.
Total large-gap helper visits (both linking stages):

| Input | Visits at each thread count |
|---|---:|
| MjeNMV / MelaMJNV | 4,047,645,435 |
| AP027280 self | 9,907,891,074 |

The initial linker dominates this diagnostic. The existing `LOSAT_TIMING`
linking counter measures only the second linker; it cannot explain total work.
Counter-run durations overlapped correctness work and are **not adoption timing**.

The actual default-Node code-print run returned the exact oracle bytes. It
recorded four Liftoff compilations of group function 1340 followed by TurboFan.
In the emitted Wasm, the helper stream uses a 28-byte stride, with sum1 at offset 4
and next_larger at offset 20. V8 documents
[no on-stack replacement for Wasm](https://v8.dev/docs/wasm-compilation-pipeline):
a call that enters Liftoff keeps that code even if optimized code becomes available.
The equal-work result and generated tiers support a compilation-tier hypothesis;
the code-print log alone does not map actual invocations to tiers or isolate
contention. The symmetric TurboFan condition below tests this hypothesis directly. This evidence is scoped to Node/WASI, not browser certification.

## Candidate decisions before adoption timing

- **X1:** implement one separate experimental layout. Move sum1 and next_larger
  into one dense stream; preserve every NCBI update, sentinel, compact index,
  state and reduction boundary. All reservation sizes remain unchanged.
- **X2:** defer. Large capacities are observed, but allocation/initialization
  was not demonstrated to dominate. Helper initialization is measured at group
  boundaries; capacity alone does not prove physical-memory or wall-time waste.
- **M1:** defer at the eligibility gate. Only 0.070–0.134% of uncompressed match
  calls reach 16 matching bases. Means are 0.776–1.418 matching bases. This does
  not support adding 16-lane setup to the dominant short calls. No claim that an
  untested SIMD implementation is slower is made.
- **H1a:** isolated probes measured a serial second compile API request of
  2.5–2.7 ms. Evaluate retaining the inspected Module, with host preparation as
  its primary interval and end-to-end time/RSS as controls. This does not prove
  duplicate machine-code generation.
- **H1b:** defer. Threaded raw inspection takes about 4 ms, guard transformation
  about 47 ms, and guarded compile about 2.5 ms. Execution must retain the guard;
  original and guarded bytes differ. A new validation parser is outside this round.

## Local scheduling proposal considered after D0

A separate C1 proposal would complete the real first frame group on the linker's
calling thread, then dispatch remaining groups through the same requested pool
and indexed reduction. The source permits this scheduling choice: NCBI
`link_hsps.c:510–531,553–563,955–994` separates groups and sorts after all groups.
It introduces no synthetic search or extra linking call. The calling thread can
be any Rayon slot in a multi-subject search.

This remains **unimplemented and unaccepted**. It serializes a prefix on every
linking invocation, including already-optimized module/instance reuse. Cold
improvement alone would not meet the plan's reuse and control guards. The
independent reviewer requires actual child-isolate tier evidence and the full
scheduling/lifecycle/reuse proof before accepting this candidate. The X1 layout
experiment does not include this scheduling change.

## Supplementary compiler condition

One diagnostic search per version/input/n used the same
`--no-liftoff --no-wasm-tier-up` flags. All eight raw outputs matched. These are
single diagnostic observations, not three-pair adoption estimates.

| Input | Default baseline body n1 / n8, three-pair medians | TurboFan baseline body n1 / n8, single observations |
|---|---:|---:|
| MjeNMV / MelaMJNV | 22.406 / 35.077 s | 16.954 / 10.715 s |
| AP027280 self | 37.147 / 61.850 s | 32.210 / 19.821 s |

The n8 reversal disappears for both inputs under forced optimized compilation.
Together with equal work and actual emitted tiers, this strongly supports the
compilation tier as a major contributor under this Node/V8 version. The log does
not establish an exact group-to-tier mapping, and it does not separately quantify
all scheduler/contention effects. V8 documents that already-running Wasm calls
cannot switch from Liftoff to newly compiled TurboFan code mid-call.

X1 also remains slower than baseline in all four supplementary observations.
It is rejected on the primary **default-Node** results: body medians regress
1.26%, 7.83%, 13.13%, and 5.63%. Its exact state/output preservation does not make
this layout a performance improvement. The experimental Wasm hot stream is
8 bytes per helper, with a separate 20-byte cold stream; total helper field
storage remains 28 bytes before capacity/allocation overhead. No production Rust
layout, pruning, floating-point expression, pool size or scheduling is changed.
