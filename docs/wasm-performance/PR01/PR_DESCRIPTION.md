Title: Freeze Wasm performance baselines and exact-output gates

Follow-up handoff: [2026-09-05 description](followup-20260905/PR_DESCRIPTION.md).
The original description below retains the 120-second results; the follow-up
reconciles correctness to 74 PASS / 10 TIMEOUT and keeps PR01 NOT_ACCEPTED.

Add a dedicated LOSAT performance harness that binds explicit native, serial
command-Wasm, and threaded command-Wasm artifacts to the existing frozen
output contracts. It records ordered arguments, controlled input paths,
source/input/oracle/artifact hashes, repeated process measurements, and raw
output without normalizing headers, newlines, or result order.

Implementation:
- Reuse existing program command builders and canonical authority for the
  43 native and 41 directly applicable serial-Wasm contracts. Validate ordered
  arguments and selectors again before execution; preserve omitted defaults.
- Preserve NCBI 2.17.0 `blast_format.cpp:794–808,828–832` and
  `tabular.cpp:1098–1108` output contracts. NCBI executables are test oracles
  only. Production Rust and existing runners are unchanged.
- Separate cold process measurements from serial cached-Module/new-instance
  measurements. A watchdog, subprocess status, exact raw hashes, and an RSS
  budget determine whether each sample is usable. Diagnostic timings and
  thread configuration logs are collected separately.

Validation:
- All three explicit release command artifacts built successfully with
  `--locked --bin LOSAT`; Wasm exports/imports distinguish serial and threaded
  commands. The existing pure-Rust runtime boundary check passed.
- Seventeen focused harness tests passed, including wrong hashes, missing
  artifacts/oracles, bad exits, timeouts, header/path/newline/order differences,
  modified canonical arguments, and wrong thread-artifact labels.
- Serial warm-instance smoke: six exact outputs. Threaded BLASTP smoke: four
  exact outputs across requested 1/4; diagnostic runs reported effective 1/4
  and 0/3 helper spawns. Smoke timings are excluded from baseline conclusions.
- Initial 43/41 audit: 68 PASS, 15 TIMEOUT, one NCBI-oracle timeout with
  canonical LOSAT output. All 69 completed LOSAT outputs match frozen hashes.
- Cold matrix: 630 PASS outputs (90 series, five timed repeats each), nine
  NOT_RUN configurations due to the p11 oracle timeout. Warm serial: 60 PASS
  outputs and one NOT_RUN case. Protected native1 profiles: eight PASS.
- Cold native1/serial medians: BLASTN task-blastn 0.230/0.732 s,
  megablast 1.061/1.718 s, BLASTP self 1.046/1.725 s, TBLASTX Mela/PeMoJNV
  4.564/8.110 s. Peak timed-process RSS is 399,867,904 bytes.
- A separate raw-gated serial CPU profile attributes 58.4% of samples to
  TBLASTX linking. Existing timers omit initial linking; this is an
  investigation priority, not a measured optimization.
- The read-only NCBI parity auditor verified harness/exporter corrections.
  Evidence export rechecks raw hashes, full series status and median inputs.
  Full local results, ranges, profiles and limitations are in EVIDENCE.md,
  RESULTS.md and HOTSPOTS.md. Acceptance remains NOT_PROVEN.

Remaining limits:
- Timed-out searches are unresolved validation, not parity failures or passes.
  The declared per-search timeout is 120 seconds. Long TBLASTX cases require
  further runs with an explicitly revised execution budget.
- Native effective-N and active-worker utilization are not measured. Wasm
  BLASTP internal timing is disabled in the existing engine. Threaded warm
  lifecycle, resident direct API, real browsers, and platform-native Gate B
  remain NOT_RUN. No release certification, optimization speedup, or complete
  PR01 acceptance is claimed.

Publishing actions: none.
