# I1 compiled-module reuse and optimized-loop audit

## Verdict and scope

The first fixed Mje n8 compiled-module session **fails its prospective time guard**. No arithmetic, raw-output, thread-budget, or selected-clock error was found that would justify replacing or extending its samples. It is one completed AB session; the planned BA session did not run after failure. The magnitude is not a general performance estimate. I1 adoption is unsupported by this failed guard.

The optimized code supplies a concrete follow-up hypothesis: I1 recomputes a helper address after the score test where baseline retains it. Both versions already have the predecessor bounds check. Neither the address recomputation, function size, nor call overhead has been shown to cause the measured regression. This audit ran no build, search, benchmark, profiler, or long transcript diff. Production source was not changed.

## Identity and authority

- LOSAT HEAD: `8f23f774b44d6812d0943149877d93835eab1d52`, plus the frozen dirty sources in `baseline-source-manifest.json` and `I1-source-manifest.json`.
- Fixture: `MjeNMV.fasta` against `MelaMJNV.fasta`, TBLASTX local subject, `-outfmt 6 -query_gencode 1 -db_gencode 1 -num_threads 8`; query/subject hashes are `f91873bd2957bfaefe377aa221ffd1e0b57e657e5d26ba86c52fce891ac0cedf` / `491f80b482acc96600773c3493358c0262dbf714cb5d76f8f5560a92ecccac1a`.
- Fresh oracle: `/home/kawato/micromamba/bin/tblastx`, BLAST+ 2.17.0+, build Aug 11 2025; recorded in `environment.json`. No genetic-code exception applies to these default-code runs.
- Target: release `wasm32-wasip1-threads`, `wasm-threads`; exact Cargo argv retained in `work/{baseline,I1}/build.json`.
- Runtime: Node 26.8.2 / V8 14.6.202.34-node.28, rustc 1.92.0 / LLVM 21.1.3, WSL2 i9-14900HX, affinity 0–31; no Node compiler flags in reuse measurements.
- Raw command SHA256, independently recomputed: baseline `71a8b1b851c108eebac33414156108a1b2cc61b390ff3cb8562315626a5d70b9`; I1 `86ab290585570df5f4cf1dd295e8c460638c441a3239e3b0a855598115977374`. Both agree with the actual host reports.

NCBI owner is `c++/src/algo/blast/core/link_hsps.c`, `s_BlastEvenGapLinkHSPs`: previous-best minus one at 812–823; predecessor loop at 827–861; caller sum/xsum/state updates at 863–894. `next_larger` is constructed at 675–684 and updated at 876–885. LOSAT ownership/timing, frame-relative trimmed amino-acid coordinates, immutable scan inputs, strict comparisons, f64 copying, unchanged downstream arithmetic, group ordering and reduction were established in `I1-independent-actual-diff-audit.md`; nothing in this code diagnosis changes those findings. NCBI executables remain comparison oracles only.

## Fixed reuse evidence

`measure_reuse.py` writes the policy before launching either process, runs one warmup and three timed searches per version, compares raw bytes, and stops immediately after any time/RSS/linear-memory guard fails. The planned sessions are AB then BA. Only session 0 AB exists, consistent with that rule. Each process prepares one module, then creates and closes a fresh command instance and seven child workers per search. This is compiled-module reuse, not reuse of an exited command instance.

| Metric | Baseline samples | I1 samples | Result |
|---|---|---|---|
| Monotonic wall seconds | 9.681366307, 9.537426704, 9.849780190 | 9.977141889, 10.311214258, 10.370315738 | medians 9.681366307 → 10.311214258; +6.505775%, fail |
| Process-lifetime peak RSS bytes | 396886016, 483508224, 564199424 | 382828544, 461496320, 542601216 | medians 483508224 → 461496320, pass |
| Linear-memory bytes | 94371840 each | 94306304 each | pass |

The wall allowance is `max(0.05 * 9.681366307, 0.050) = 0.484068315 s`; actual increase is 0.629847951 s. RSS allowance is 48,350,822.4 bytes; linear-memory allowance is 16,777,216 bytes. The memory measurement is cumulative process peak, correctly named; it does not establish a long-running memory plateau.

All eight retained raw searches, including warmups, independently hash to the oracle `b309f77fee038f0559d870737ebeba732c11eba4f7e495f045acbfb831817fd7`. Host records show seven distinct successful spawn/ready/exit lifecycles per search. Baseline process monotonic interval is 2599.824391834–2653.347704926; I1 is 2653.562178442–2696.599910698. They do not overlap. Logs alone cannot prove absence of unrelated machine activity; the parent declared an exclusive measurement window.

The common `benchmark_wasi_reuse.js`, `wasi_thread_host.js`, and `wasm_performance.py` hashes still exactly match the frozen baseline manifest. Their measured boundary includes instantiation, command execution, output read, worker join, and close; evidence serialization and output hashing occur outside that elapsed interval. Both versions use the same host and environment policy.

## What the lifecycle intervals establish

`wasi_thread_host.js:171` records **performance.now()**, not Date.now(), for every event. The saved monotonic event partitions therefore remain usable despite adjustable-clock disagreements.

| Version/sample | First spawn to all ready (s) | Last ready to last exit event (s) |
|---|---:|---:|
| Baseline 0 | 0.363954346 | 9.314648636 |
| Baseline 1 | 0.372668339 | 9.159812596 |
| Baseline 2 | 0.388923177 | 9.458174696 |
| I1 0 | 0.358074115 | 9.616353729 |
| I1 1 | 0.399067522 | 9.909809867 |
| I1 2 | 0.767507994 | 9.596454467 |

I1 sample 2 has two larger spawned-to-ready waits, 196.684 and 293.689 ms. Instantiation is 1.85–5.88 ms, exit wait 0.49–1.12 ms, and close 0.012–0.029 ms. Startup delay affects one sample, but startup alone does not explain the other slower invocations. The post-ready interval includes all remaining search, formatting/I/O, scheduling and host receipt of exit events; it is **not scan CPU time**. Do not subtract it or the startup intervals to override the fixed failed metric.

Adjustable-clock disagreements are retained: baseline warmup and timed sample 2, and I1 timed sample 1. For I1 sample 1, Date.now reports 11.850 s against monotonic 10.311214258 s. Process CLOCK_MONOTONIC and CLOCK_BOOTTIME agree. The selected metric was prospectively monotonic; substituting realtime would be unjustified. This audit does not identify the cause of the realtime adjustments.

## Concrete emitted-code findings

Both saved code prints identify actual TurboFan functions by name: baseline group 1340 and I1 scan 1368. Diagnostic argv contains `--prof` and the print filter; these runs are separate from acceptance timing. `tier-pc-evidence.json` reports I1 helper TurboFan samples in child isolates (7,225 total) and baseline group Liftoff samples (23,058). Its named ranges match the inspected code headers. This review did not recompute every profile tick; code allocation alone would not prove execution, nor do these profiles establish the tier of every reused invocation.

1. **The predecessor bounds check is shared, not newly introduced.** Baseline TurboFan `stdout.txt:71575–71578`, offsets `+0x4b4b/+0x4b4e/+0x4b54`, compares helper length against current index, branches to bounds failure, and multiplies the index by 28. I1 `stdout.txt:5065–5069`, `+0x1c4/+0x1c7/+0x1cd`, does the same. Both unroll successive b0 jumps with corresponding checks. Baseline WAT 8105–8120 and I1 WAT 435–452 confirm the Rust-level relation. It is incorrect to explain +6.5% as a bounds check added by extraction.
2. **I1 repeats address arithmetic for score contenders.** I1 WAT 525–533 redoes `base + current_idx * 28` after leaving the b0 scan; actual TurboFan `stdout.txt:5160–5161`, `+0x34b/+0x34e`, emits another multiply/add. Baseline WAT 8175–8177 reuses the helper address from 8115; baseline TurboFan `+0x4e46` loads coordinates through the retained `r8` address. This is a concrete additional operation on this I1 path, without a measured frequency or attributable time.
3. **Extraction creates real per-call costs.** I1 WAT 22–41 allocates a 192-byte guest stack frame and stores scalar inputs; 48–187 prepares formatting argument address packs even before the first trace-condition branch. The TurboFan prologue `+0x36..+0x108` retains this work. Result fields are stored to the caller return area at WAT 492–506, then loaded by group WAT 7842–7860. These costs occur per eligible HSP scan call, including empty scans in current emission, rather than a new group call. They do not prove the regression's cause or net cost relative to baseline register pressure.
4. Both versions retain diagnostic temporary stores and selected-HSP bounds checks. The helper takes direct slice data/length values: there is no evidence here of a new double-Vec-owner indirection. Some baseline accesses reload owner fields; instruction counts or code size cannot substitute for a path-specific cost measurement.

## Bounded predecessor-prefix follow-up

A borrowed prefix can encode the existing valid index range without changing the selected set or algorithm. For a nonempty scan, use original helpers `2..h_lh_idx`. If its length is L, `split_last` selects original index L+1. Compute `current_idx` before replacing the prefix by its tail. Assigning the tail at the original decrement point leaves indices 2 through current_idx−1. If b0 is true, truncating that tail to `next_larger.saturating_sub(1)` leaves exactly original indices 2 through next_larger. Targets 0 and 1 both produce an empty cursor.

NCBI construction/update guarantees `next_larger < current_idx`, so the truncation is valid for production state. Each split/truncation is O(1), borrows existing memory, and needs no allocation, unsafe indexing, candidate pruning or duplicated owner. Keep sum/next/b0 reads, decrement, jump and existing trace fallthrough in their current order; the borrowed current helper must still be used for traced predicates after a jump.

**Empty-boundary trap:** current unit test `large_gap_scan_preserves_initial_state_and_zero_one_jumps` deliberately calls `h_lh_idx=2` with an empty helper slice. Unconditional `[2..2]` would panic. Preserve the legitimate no-visit return before forming the prefix. Do not replace malformed state with a silent empty fallback.

Such a candidate is justified as removal of an existing range check and potential address recomputation, not as a proven fix for the failed reuse result. Require fresh emission confirming the intended change, actual nonadjacent jumps, 0/1, negative sums, strict ties/mappings, exact trace/FP differential states and raw outputs before a frozen performance screen. Preserve Native/serial performance controls and all remaining adoption gates. Lost round-02 source or timings provide no evidence for this candidate.

## Evidence bindings

SHA256 values independently read during this audit:

- baseline group WAT: `92cf16e4dcd7bcd5172250cb89001dc8210b59eff818bd6df3fe6d40a2f336b5`
- I1 helper WAT: `81057941dcfe67c9b4c9a9341185d4df96a82b220904ea34f16286b24ce5e47a`
- baseline code print: `76cfaf078d080a74947122255d36333839c246b5432228ed8fd5c8d79d977e94`
- I1 code print: `78a6ccb65b336a2abf9251e587a337d5240712e88b81de0b910f1bb3132e5500`
- reuse policy: `13fe6dad823d8ff103d0fff24f536d0766c0c7e8fe2b7f17ffe53421124e587c`
- reuse runs: `90f244423d58befa5dd2453570fe7c6a082e2290dfc234b10652d21161fcdeda`
- reuse summary: `03879551b4fe103177e7c49aee3a1c66fbc69f1252114aab21afd172db55732f`
- measurement script: `3258aa0d27a92825f89970e5e597318476fd0bb265d8b36e9353d2e4789a12f4`
