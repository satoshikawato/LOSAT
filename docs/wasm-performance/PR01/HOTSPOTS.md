# PR01 diagnostic priorities

These are baseline investigation priorities. No later PR is implemented, and a
large parent-stage timer does not prove a particular inner loop or allocation
is responsible. All ratios below use the separate requested-one-thread
**diagnostic engine total**, not the cold-process median. Stage timers are
rounded, nested, and not an additive partition. Raw logs are retained in
`cold/<program>/<case>/<target>-1/diagnostic/stderr.txt`; the extracted values
are in [stage-profiles.json](stage-profiles.json).

| Follow-up | Current fixture and measurement | Supported investigation |
| --- | --- | --- |
| PR04 greedy matching | EDL933/Sakai megablast: native greedy traceback 0.377 / total 1.029 s (36.6%); serial Wasm 0.531 / 1.736 s (30.6%) | Profile the greedy traceback inner loops, including first-mismatch scanning. No isolated match-loop share or SIMD benefit is proven. |
| PR05 traceback/pruning | PesePMNV/MjPMNV task-blastn: native DP traceback 0.089 / 0.220 s (40.5%); serial Wasm 0.426 / 0.614 s (69.4%). EDL serial traceback/prune parent is 0.603 s, including 0.531 s greedy alignment. | Prioritize DP traceback on task-blastn. EDL tree checks and score sorting individually measure 0.001–0.004 s; these do not support prioritizing generic pruning over alignment. |
| PR07 Kappa workspaces | BLASTP pairwise_default native Kappa 0.523 / 0.789 s (66.3%); all_hsp_self 0.724 / 1.046 s (69.2%). Prelim gapped is 0.101/0.141 s, scan/ungapped 0.079/0.089 s. | Kappa is the leading measured native stage. Allocation/reset share and benefit of workspace reuse remain unproved. Existing Wasm BLASTP stage timers are disabled, so Wasm shares are NOT_RUN. |
| PR10 TBLASTX | p03 Mela/PeMoJNV: native scan 0.083, extension 0.397, linking 0.009, search 4.526, total 4.717 s; serial Wasm scan 0.148, extension 1.291, linking 0.011, search 8.491, total 8.794 s. | The named counters leave substantial search work unaccounted for. Extension is the largest named counter, but its 14.7% of serial engine total does not establish the dominant hotspot. CPU profile evidence is recorded separately below. |

The TBLASTX scan counter wraps only `s_blast_aa_scan_subject`; the following
hit loop and diagonal/two-hit bookkeeping are outside that counter
(`LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:1796–1848`). This identifies
an instrumentation gap, not an allocation of the residual time to that loop.
The p11 AvCLPV/PsCLPV linking fixture timed out at the oracle gate. Its missing
profile prevents a general conclusion that linking is insignificant.

NCBI 2.17.0 source owners inspected for these boundaries:

- `c++/src/algo/blast/core/greedy_align.c:313–377`, `s_FindFirstMismatch`;
  callers at 457, 571, 877, 1071 preserve direction, ambiguity, packed-base and
  fence behavior.
- `blast_traceback.c:358–379,633–692`: traceback and post-traceback ordering,
  reevaluation, and containment.
- `blast_kappa.c:3429–3459`: per-match scheduling and thread-local state.
- `aa_ungapped.c:344–400`: returned hit iteration, diagonal state, window checks,
  and `s_BlastAaExtendTwoHit` admission.

All are source references for later investigation, not authorization to change
candidate admission, pruning, logical order, scores, or statistics.

## Separate TBLASTX CPU profile

The serial p03 command produced the expected 351,051 raw bytes and SHA
`bbef153094025db89dd4a69b2a3c4246b0b8853396fe740bc7e1bd654db105b6`.
The V8 profile collected 7,112 samples with a requested 1 ms interval:

- `sum_stats_linking::linking::link_hsp_group_ncbi`: 4,153 leaf samples
  (**58.4% of all samples**).
- `blast_engine::run_impl::run_internal` closure: 1,950 (27.4%). Release
  inlining prevents assigning this parent function to a finer source region.
- `getenv` and `strncmp`: 247 and 239 (6.8% combined); these identify runtime
  lookup overhead but do not establish a permissible change to diagnostics.

**PR10 should investigate linking first for this fixture.** The source explains
why the existing 0.011 s `sum_stats_linking` timer is insufficient: initial
linking at `run_impl.rs:2556–2570` has no `t_linking` timer; the timer starts only
before the later relinking call at `run_impl.rs:2708–2729`. NCBI likewise calls
`BLAST_LinkHsps` in `blast_engine.c:870–874` before preliminary E-value reaping,
then relinks after ambiguity reevaluation at `blast_engine.c:1509–1520`.
`link_hsps.c:415` owns `s_BlastEvenGapLinkHSPs`.

The sampled function contains work from both linking calls. The profile does
not isolate an exact per-call or inner-loop percentage, and no native or other
fixture share is inferred. Neither linking call may be removed or moved on
the basis of this measurement. The earlier untimed scan-hit-loop observation
remains an instrumentation gap, but cannot explain away this measured linking
hotspot.

Raw profile: [tblastx.cpuprofile](cpu-profile/tblastx.cpuprofile).
Counts, function names and SHA: [CPU summary](cpu-profile/summary.json).
Exact invocation and raw gate: [command](cpu-profile/search/command.json),
[gate](cpu-profile/gate.json). This diagnostic is excluded from every median.

## Follow-up serial BLASTP CPU profile

The unchanged serial artifact ran `all_hsp_self_serial` and passed the frozen
raw gate: 35,004 bytes, SHA
`67911f4a93327e9622afe10d15a66d9e52ae8e1c3d905db9ca5fb018d6612e46`.
Of 1,475 V8 samples (requested 1 ms interval), 997 (67.6%) have
`blast_redo_one_match_with_workspace` on their stack. This inclusive count
contains 685 leaf samples in `align_ex_protein`; the redo function itself has
253 leaf samples. `align_ex_protein_score_only` has 186 leaf samples (12.6%)
on a separate engine closure stack. Do not add inclusive and leaf counts.

This supports the PR07 Kappa/redo investigation priority on this serial-Wasm
fixture as well as the earlier native evidence. It does not establish that
allocation or workspace reset dominates, or that workspace reuse will help.
Exact Wasm engine-stage timer shares remain NOT_RUN. Release inlining limits
finer attribution. The diagnostic could overlap the oracle audit and is
excluded from every latency median; these are sampled CPU proportions.

NCBI owners read for this follow-up: `blast_kappa.c:3623–3642` calls
`Blast_RedoOneMatch` with thread-local alignment and composition state;
`composition_adjustment/redo_alignment.c:1101–1111` receives the composition
workspace. Corresponding LOSAT entry:
`LOSAT/src/core/composition_adjustment/redo_alignment.rs:1694`.

Raw [profile](followup-20260905/blastp-cpu-profile/blastp.cpuprofile),
[sample counts](followup-20260905/blastp-cpu-profile/summary.json),
[command](followup-20260905/blastp-cpu-profile/search/command.json), and
[gate](followup-20260905/blastp-cpu-profile/gate.json) are retained separately.
