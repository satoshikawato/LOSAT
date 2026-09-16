# Independent preflight: I1 large-gap predecessor function boundary

Status: the proposed extraction has a sound, bounded ownership and control-flow boundary, subject to the constraints below. **An n8 improvement is unproven. The same large-gap extraction was previously tried as C-X3.** This revision must be described using its concrete differences and new diagnostic question, not as an unexplored idea or an adopted fix.

This audit is read-only except for this note. No implementation edit, test, build, benchmark, browser run or archive extraction was performed. Repository AGENTS.md and the verify-ncbi-parity-and-speed skill/references were read. The working tree contains extensive user changes and was preserved.

## Current source and evidence authority

- LOSAT HEAD: `8f23f774b44d6812d0943149877d93835eab1d52`, plus existing working-tree changes. HEAD alone is not the source identity.
- Current `LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs` SHA256: `e9ad6da8849ac2c18be1f05889446302780178238bf55011c24a6d66a8e0647d`. Independently verified equal to the corresponding frozen `run-02/BASELINE.json` entry. Parent separately verified all 158 Rust/config/build inputs.
- NCBI source: `/mnt/c/users/genom/github/ncbi-blast/c++/src/algo/blast/core/link_hsps.c`, SHA256 `22c749d1501403d2ab0926e4110617a2f02a9ccf8e613bdb24a19715b1e4c695`, `s_BlastEvenGapLinkHSPs` beginning at line 415; selected predecessor loop 827–861, call-context block 774–823, downstream updates 863–894.
- Existing durable run-01 comparison records use NCBI BLAST+ 2.17.0+ and Node 26.8.2. Run-03 must freeze the actual current toolchain/runtime/executable identities before its measurements; old versions are not assumed to describe the resumed environment.
- Raw round-02 scratch artifacts, transcripts and measurements are unavailable after the environment transition. `run-02/PROGRESS.md` is historical narrative, not reconstructable adoption evidence. This note does not carry C8d/C9 outputs or timings. Current source and durable older records support only the preflight described here.

Primary follow-up fixtures remain TBLASTX local subject, outfmt 6, query/db genetic code 1:

| Fixture | Query SHA256 | Subject SHA256 |
|---|---|---|
| MjeNMV.fasta / MelaMJNV.fasta | `f91873bd2957bfaefe377aa221ffd1e0b57e657e5d26ba86c52fce891ac0cedf` | `491f80b482acc96600773c3493358c0262dbf714cb5d76f8f5560a92ecccac1a` |
| AP027280.fasta self | `14c68a026e685896afdc5366902999d5960589cf301a76b3f1a732a8f4708b6a` | same |

These input hashes were freshly recomputed. n1 and n8 are separate required normal-runtime conditions; n2/n4 remain useful for work/order diagnosis. Plain wasm32-wasip1 is serial; actual threaded Wasm requires wasm32-wasip1-threads plus wasm-threads. No genetic-code exception applies to these main gencode-1 fixtures. The existing AP027131/AP027133 code-4 local-subject regression requires the approved database-oracle treatment and does not authorize other differences. NCBI binaries remain comparison oracles only.

## Exact call path and extraction boundary

NCBI `blast_engine.c:870–899` links/reaps preliminary HSPs; the later translated-subject path reevaluates and links again at 1515–1520. `BLAST_LinkHsps` (`link_hsps.c:1761`, even-gap dispatch around 1780) invokes `s_BlastEvenGapLinkHSPs`. That function sorts/groups HSPs, constructs the helper stream, and repeatedly runs the existing small-/large-gap DP and chain selection.

LOSAT follows `run_impl.rs` preliminary `apply_sum_stats_even_gap_linking_with_parallel` call around 2305, reevaluation, then the final call around 2465. `linking.rs` sorts with the NCBI reverse comparator, partitions by translated context/sign, invokes `link_hsp_group_ncbi`, and reduces indexed parallel results back to original group order. All those call sites and scheduling/reduction owners must remain unchanged. The existing LOSAT_TIMING linking timer surrounds only the later call; it is insufficient to diagnose total linker cost.

Within the current group kernel:

1. `linking.rs:1566–1617`: each existing large-gap DP HSP initializes local selection state; checks the NCBI unchanged-previous-choice fast path; retains long-subject diagnostic counters and strict score cutoff.
2. `1617–1654`: only an HSP requiring a real scan and satisfying `i_score > cutoff_big` obtains trim ends, emits the existing target header, and possibly sets initial `h_sum = previous_sum - 1`. NCBI 812–823 deliberately changes only that sum to retain original tie selection.
3. **Extract exactly `1656–1718`: cursor initialization plus the entire predecessor while loop.** Call once per scan-eligible HSP in each existing DP pass. This is not once per HSP globally: the chain algorithm may revisit an HSP in later passes. Do not call the helper on the unchanged-choice fast path or filtered scores merely to create entries/warmup.
4. `1725–1800`: the caller retains new_sum and the left-associative f64 new_xsum expression, all HSP/helper writes, maxsum/next_larger rebuilding, best selection, xsum assignment and linked_to increments. None moves into the scan helper or across its call.

The scan has one local `continue`, targeting its own predecessor loop. It has no return, outer-loop break, or mutation of pool/helper/hit arrays. Therefore an ordinary private function (no FnMut closure or shared state) can own the cursor and four scalar selection values for the duration of this one scan. Immutable slices borrow the existing vectors; the borrow ends before caller updates. No heap allocation, payload cloning, extra traversal or asymptotic change is needed.

## Required state and invariant preservation

The function requires the original `h_lh_idx`, helper/link/hit slices, `h_qe`, `h_se`, `is_target_hsp`, and initial selection values. Return final values of exactly the current types: h_sum i32, h_num i16, h_xsum f64, h_link usize. A small by-value tuple is sufficient; no duplicate persistent state owner is needed. Carrying all four initial values is safer and more explicit than the previous C-X3 shortcut that initialized three values inside its helper.

Preserve these details together:

- Start at `h_lh_idx - 1`; test cursor >1; save current_idx and bind the current helper before decrement or jump.
- Read sum[1] and next_larger; compute `sum <= h_sum`; decrement once and replace the cursor with next_larger on b0. Jumps to **0 and 1** both terminate. Backward-jump construction/update remain in their original owners (NCBI 675–684,876–885).
- Preserve the current adopted non-traced b0 early continue. Do not reintroduce unconditional coordinate loads from old C-X3. With tracing enabled, b0 still falls through to compute all predicates and emit the same second skip diagnostic using the original current helper/index, not the jump destination.
- Keep coordinate comparisons `qo <= h_qe`, `so <= h_se` and selection `!(b0 || b1 || b2)`. Do not change tie operators or choose a different equal-scoring predecessor.
- `helper.hsp_idx` maps the helper to the link/hit arrays. It is not generally the helper cursor. Use the mapped index for every num/sum/xsum/link and trace read.
- Copy selected num, sum, xsum and link in the existing order. The scan introduces no f64 arithmetic. Return/copy the f64 bits intact; do not recompute a normalized score or combine the caller's `(h_xsum + score * lambda) - logK` expression.
- Inputs are already frame-relative trimmed amino-acid coordinates constructed from `frame_relative_coords` around 1110–1130. Do not apply nucleotide-coordinate or frame/context adjustments inside this helper.
- Preserve every trace condition, format string and argument-evaluation location inside the scan. Keep the caller's initial-best and final-update diagnostics outside it. Do not introduce a per-call clock or shared atomic in the performance artifact.

No changes are warranted to nucleotide/NCBISTDAA encoding, SEG, cutoffs, candidate construction, chain member filtering, statistics or output formatting. Output filtering of linked non-head members stays in reporting.

## Durable prior attempts: do not erase the precedent

1. `docs/evidence/wasm_performance_20260913/run-20260913-01/P1.md` documents a **small-gap** predecessor helper (NCBI 703–746) with Wasm-only non-inline boundary. On its MelaMJNV/PemoMJNVA fixture, the recorded local candidate improved n8 only 1.78% and was rejected under that run's threshold. Different loop, fixture and Node 24.21.0 mean that result neither proves nor disproves this large-gap revision. The invalid earlier build missing build.rs is explicitly excluded by that record.
2. `docs/evidence/megablast_wasm_np_remediation_20260913/run-20260913-01/C-X3-exploration.patch` already defines **scan_large_gap_predecessors** over NCBI 827–861, called once per eligible HSP, with three immutable slices and returned selection tuple. Its attribute is unconditional `#[inline(never)]`; only initial h_sum is passed, while num/xsum/link are initialized inside. It predates the currently adopted deferred-coordinate/early-continue path. `C-X3-exploration-plan.json` names the integrated baseline and a native Mje diagnostic; `REPORT.md:181–185` records rejection. No new quantitative rejection percentage is inferred here from missing individual measurement files.
3. `audit-C-X3-static-exploration.md` found no semantic/borrowing defect but identified a meaningful test gap: its 1050-combination test did not cover the production jump-to-0 sentinel. Retain that lesson in the new differential cases.
4. The durable four-program `execution-evidence-manifest.json` also records X0loop source/build and native/threaded-n8 profile attempts for both main fixtures. Its archive members were not extracted in this bounded audit, so no exact source/timing claim is made from that manifest alone.
5. `run-02/PROGRESS.md` describes whole-round C2/C4 boundaries and C6 diagnostic outlining, which have a different granularity. Because their raw scratch evidence is lost, they cannot provide fresh acceptance evidence. Their history still warns against equating an emitted call boundary with a speedup.

Concrete new dimensions are: Native inline(always) instead of C-X3's unconditional no-inline; preserve the current deferred-coordinate scan; pass the complete initial scalar state; and investigate repeated actual callee entry/tier within the original long-lived group call. A new investigation is justified by these differences plus durable run-01 equal-work and compiler-tier observations. It is not permission to re-adopt an old rejected patch without fresh proof.

## Compiler and performance hypotheses

A Wasm non-inline function creates an intended repeated real entry point per eligible HSP. A helper entering after optimized code becomes available may execute that code even while its caller remains in an older compiled activation. This is a hypothesis to verify in the actual normal runtime, not guaranteed by a Rust attribute. Native inline(always) requests optimizer expansion; it does not guarantee identical native machine code or nonregression.

Minimum useful generated-code evidence:

- Map the emitted function index to the actual complete large-gap helper and its real group call site; include the host shared-memory transform's index mapping. Check that no wrapper/closure contains the hot loop elsewhere and that native expansion occurs as intended.
- Inspect real code size, argument/result lowering, stack spills and bounds checks. The tuple may lower to stack/return storage; three slices carry pointer/length pairs. More calls can create cost even with identical algorithmic work. Existing diagnostic formatting may keep the helper larger than its source loop suggests.
- Show **executed** default-runtime helper PCs/tier in actual child isolates, along with compilation records. A TurboFan code-allocation message alone does not prove later calls entered it. Distinguish compilation, completed calls and work. Profiling/code-print measurements stay outside the adoption samples.
- Confirm identical helper visits, rejection/jump counts, selected predecessors and per-group work at n1/n8. Include both preliminary and final linker invocations. Preserve requested n and n−1 children; do not shrink pools or schedule a serial prefix as part of this candidate.
- Compare cold calls and compiled-module/reactor reuse. Once the original group is already optimized, per-HSP call overhead can dominate. If evidence does not show the proposed execution-tier change or n8 benefit, report that boundary rather than layering further micro-candidates onto it.

Durable run-01 DIAGNOSIS.md supports the earlier version's tier hypothesis: equal group work at n1/n2/n4/n8 and reversal disappearance under symmetric forced optimized compilation. It explicitly lacks exact invocation-to-tier mapping and does not quantify all contention/scheduler causes. It is not a blanket explanation of the resumed environment or browser performance. Do not set forced-tier flags as the default fix.

## Focused differential proof before timing

Freeze the original inline loop as test reference and compare the extracted helper's final tuple and full visitation/selection traces, not merely hit counts. Required nontrivial cases include empty scan (h_lh_idx=2), one predecessor, backward jumps to 0/1/interior entries, strict sum ties and trim equalities, coordinate rejection after passing b0, negative/zero/positive initial sums, nonzero initial state preserved when nothing improves, nonidentity helper→HSP mappings, and distinct xsum bit patterns. Exercise trace off and b0=true trace fallthrough; retain exact f64 bits within each target.

Full group differential states must still cover helper/link/active-chain ownership across multiple passes and unchanged fast-path behavior. Compare Native, plain serial Wasm and real threaded Wasm against their own frozen baseline; preserve any already characterized cross-target f64 delta instead of normalizing it away. Current fixture oracle outputs, threshold10/100/10000 and target trace, primary n1/n8 raw output/order, and thread/reuse contracts are the relevant integration checks. These are requested future checks, not tests run by this auditor.

The resumed run needs new source/artifact/runtime bindings and fixed primary diagnostic/measurement policy. Normal-runtime n8 versus n1 is the active user concern; do not substitute a small host-startup benefit or best individual sample. No output-changing speedup is acceptable. Any supported benchmark must retain all declared samples and output hashes, use monotonic boundaries, and keep other builds/tests/archives/audits outside timing. The user's 1.20 ratio remains an aspiration, not a waiver of output, timing or memory guards.

## Handoff conclusion

No static blocker was found for the proposed narrow extraction. The prior C-X3 precedent is material and must be acknowledged. Keep one helper with one algorithm body, the existing outer callers/guards and all pool/state updates. Phrase the source comment as providing a repeated real function-entry boundary; do not promise optimized execution before observing it. No adoption or n8-resolution claim is currently supported.

Audit completed: 2026-09-15T22:58:08.816827+00:00 UTC. All auditor tool commands have exited; no background audit work remains.
