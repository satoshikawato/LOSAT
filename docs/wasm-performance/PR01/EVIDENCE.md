# PR01 evidence

Latest follow-up (2026-09-05): **NOT_PROVEN / NOT_ACCEPTED**. Reusing the
completed evidence and retrying incomplete executions under separate
600-second manifests gives **74 PASS / 10 TIMEOUT** across the 43/41 catalog.
The remaining rows are native/serial p11, d01, d02, d03 and d05; only p11's
oracle remains incomplete. PR02 implementation remains blocked by PR01's
unaccepted prerequisite. Native effective-N, lifecycle and browser evidence
have distinct requirements and are not bundled into that blocker. See the
[follow-up evidence](followup-20260905/EVIDENCE.md),
[dependency boundaries](followup-20260905/REQUIREMENTS.md) and
[new independent audit](followup-20260905/auditor-review.json).
The original 120-second evidence record below is retained as history.

Status: NOT_PROVEN / NOT_ACCEPTED. Local baseline and diagnostic collection is complete within the declared 120-second per-search budget. Missing coverage remains explicit.

## Identity

- PR / goal / scope: PR01; freeze a local performance harness and exact-output gates. Production algorithms are unchanged.
- Baseline HEAD: `af3e2ea837afdb8a00cf19920f68be4f0bf3bfb5`. Dirty source identity: `a8362900aee39f2b563563f24c321966b9ed12d8d4c49007fe4d5d7f72536c97`. All 32 selected source files matched the supplied plan snapshot before work. The full tracked-file snapshot is [source-identity.json](source-identity.json).
- Candidate: none; no optimization is implemented or claimed.
- Rust: 1.92.0; release O3, LTO, codegen-units=1, panic abort. `--locked --bin LOSAT` avoids lib/bin artifact collisions. Native default features; serial `--no-default-features`; threaded `--features wasm-threads`; both Wasm targets use simd128.
- Native SHA-256: `058587655963bf1eb253f232fdb5553819105646c0d1c50bbece96e410446819`.
- Serial Wasm SHA-256: `000691b99169ddca8072954fdfa64348ad94b75546432655f6a05daea1d5eb47`.
- Threaded Wasm SHA-256: `557f7cb99aa0a48dc64207e5efe7833b753cf343b6f3550cbd3b83040d3d25f1`.
- Artifact paths, full build argv, exports/imports, compiler/lockfile hashes, Node v18.19.1, runtime flags, CPU and OS identity: [contract/manifest.json](contract/manifest.json).
- Hardware has 32 logical CPUs. Host cap is 8. Power mode, temperature and active-worker utilization are not independently measured.
- Official NCBI BLAST+ 2.17.0 executable hashes are checked against the existing integrated-certifier authority. The official source archive SHA is `502057a88e9990e34e62758be21ea474cc0ad68d6a63a2e37b2372af1e5ea147`.
- [oracle-archive-identity.json](oracle-archive-identity.json) verifies the published binary-archive checksum and equality of extracted versus installed executable bytes. [ncbi-source-identity.json](ncbi-source-identity.json) binds the inspected sources to the official source archive. Only source-checkout CRLF is disregarded in that source-text comparison; output gates never normalize bytes.
- Controlled lexical root: `/tmp/losat-pr5-runtime-cert-5845d22/LOSAT`. Cwd and input byte identities are frozen in the manifest. Manifest omission remains omission, not numeric zero.

## NCBI mapping

- `c++/src/algo/blast/core/blast_engine.c:1410–1455`: subject iteration and per-subject parameter updates. The harness executes complete searches, without changing subject/frame/chunk units.
- `blast_traceback.c:358–379,633–692`: score-ordered traceback, endpoint purge, reevaluation and containment. No pruning, comparator, floating-point, context, strand or frame code changed.
- `blast_kappa.c:3429–3459`: match scheduling and thread-local search state. No Kappa workspace or admission behavior changed.
- `blast_format.cpp:794–808,828–832` and `c++/src/objtools/align_format/tabular.cpp:1098–1108`: header path identity, retained-alignment order, fields, delimiters and line endings. These are the owners of the raw-output contract enforced by the harness.
- Existing Sakai source-underdetermined classification and the six local-subject TBLASTX genetic-code deviations are preserved. Their validators and frozen LOSAT hashes remain separate from platform-native fingerprints.
- First current divergence: no completed-output raw mismatch observed. Timeouts are incomplete executions, not resolved or diagnosed parity defects.
- Mutable search state is not reused. The serial warm runner shares a compiled Module only and creates fresh WASI, instance and memory for every command. It never calls `_start` twice on one instance.

## Change and correctness

- PR01 production diff: none relative to the starting dirty snapshot. The existing Git diff against HEAD contains prior user changes; those bytes were preserved. Existing command runners and frozen canonical/authority files are unchanged.
- Tests: dedicated Python harness, focused fixture policy, two exact-empty extra contracts, warm serial host, negative tests, and evidence summarizer under `LOSAT/tests/`.
- Documentation: this record, the runbook, result tables, and an English PR description.
- Generated evidence: manifests, build logs, source/oracle identities, raw output/stdout/stderr, commands, resource usage, JSONL samples and diagnostic profiles.
- Initial full-catalog integrity: 43 native / 41 directly applicable serial rows. Initial execution recorded 68 PASS, 15 TIMEOUT and 1 ORACLE_FAIL. The latter is p10 native: LOSAT matched its frozen hash but the NCBI oracle timed out. Thus 69 completed LOSAT outputs matched their canonical SHA; no completed-output mismatch was found. Initial execution and retained-evidence revalidation are distinct operations.
- Focused matrix: native 1/2/4/8, plain serial, threaded 1/2/4/8. Ten completed cases cover small, complete no-hit, and representative computational inputs for all three programs. Cold recorded 630 PASS outputs across 90 complete series (450 timed samples), plus 9 NOT_RUN configurations for the p11 oracle timeout. Warm serial recorded 60 PASS outputs (50 timed samples), plus p11 NOT_RUN. Exact output, not timing, controls acceptance.
- Harness negative tests: 17 passed; lint/format and Node syntax checks passed. Wrong hashes, missing artifacts/oracles, failed processes, timeouts, path/header/newline/order differences, canonical argument tampering and thread-artifact mislabeling are tested.
- Build validation: all three release command builds passed. Existing pure-Rust runtime boundary audit passed.
- Smoke validation: serial warm has six exact outputs; threaded BLASTP has four exact outputs at requested 1/4, with configured effective 1/4 and helpers 0/3 in diagnostic runs. These concurrent smoke timings are excluded from baseline results.
- Protected native1 profile: 8 PASS outputs across NZ self, Sakai/MG1655, AP027280 self and the approved AP027131/AP027133 db4 fixture. Separate serial p03 CPU diagnostic: raw gate PASS; linking accounts for 58.4% of 7,112 V8 leaf samples, while the existing timer covers only the later linking call.
- Final broader algorithm re-execution: NOT_RUN; no production algorithm changed. Retained initial results are rehashed and checked with current existing exception validators by the summarizer. This does not convert timeouts into passes.
- Direct API, real browsers, threaded warm lifecycle and Windows/macOS platform Gate B: NOT_RUN. No unknown platform is accepted and no native fingerprint becomes LOSAT expected output.

## Performance

- Eligible and protected case policy was frozen before candidate optimization; no candidate exists. `benchmark=false` means a protection/profile case, not an excluded correctness contract. The long AvCLPV/PsCLPV eligible case remains eligible even when execution times out.
- Cold: one warmup, five fresh process samples, then a separate diagnostic run. Blocking waitpid plus watchdog avoids timeout-polling bias. Timing series are scheduled after the initial audit to avoid overlapping searches.
- Warm serial: one cached Module, one warmup plus five new instances; compile cost and command-return latency are separate. No persistent engine or resident throughput claim.
- B/C interleaving: not applicable; baseline only. Noise is reported through complete sample ranges. No best-single-sample or whole-set attainment claim is made.
- Resource budget: 4 GiB peak process RSS rejection threshold. Threaded linear memory maximum is 1 GiB; serial host limit is null/module-defined. Missing measurements stay null.
- Existing native stage timers are reused. Wasm BLASTP timers are disabled in the existing engine. Nested or multi-thread timer totals are not an additive wall-time partition.
- Raw samples and medians: [RESULTS.md](RESULTS.md), [samples.tsv](samples.tsv), [RUN_MANIFEST.jsonl](RUN_MANIFEST.jsonl), and [comparisons.json](comparisons.json). Cold requested-N comparisons are not proof of the same effective N; native effective-N and active worker counts remain unknown. Peak RSS over successful cold timed samples was 399,867,904 bytes, below the 4 GiB rejection budget. Allocation counts and native active workers are NOT_RUN because the existing measurement path does not expose them.
- Representative cold native1 / serial Wasm medians: task-blastn 0.230 / 0.732 s; EDL megablast 1.061 / 1.718 s; BLASTP self 1.046 / 1.725 s; TBLASTX Mela/PeMoJNV 4.564 / 8.110 s. Complete min/max ranges and all requested thread counts are retained. There is no candidate/baseline improvement comparison.
- PR04/05/07/10 investigation priorities and the limits of stage attribution are documented in [HOTSPOTS.md](HOTSPOTS.md). The separate V8 CPU profile is excluded from latency medians. Startup/compile/instantiate allocation bytes, resident latency, and resident independent-job throughput remain NOT_RUN; the new-instance boundary does not establish a resident engine.

## Decision and handoff

- Implementation: NOT_PROVEN; harness and local evidence are implemented, but incomplete gates do not permit ACCEPTED.
- Correctness: NOT_RUN for complete certification; completed executed scope passed its raw gates, while timeouts remain unresolved coverage.
- Performance: INCONCLUSIVE for improvement and target attainment; baseline only, with a missing eligible case. No whole-eligible-set geometric mean or thread-efficiency claim.
- Independent reviewer: `ncbi_parity_auditor`, read-only. Command binding, watchdog, missing-data handling, warm-host and exporter series-status/median issues found during review were fixed and independently rechecked. The final read-only review independently rehashed the retained outputs, checked all 99 series statuses, medians/RSS, and CPU sample counts, and supports the scoped hotspot conclusions. It does not approve full PR01 acceptance with missing evidence. See [auditor-review.json](auditor-review.json).
- Remaining goals: complete long-case validation under a separately declared larger execution budget, native thread utilization, Wasm BLASTP engine timing, finer TBLASTX hotspot attribution, and later lifecycle/browser goals. No later PR is executed.
- Preservation: all 457 pre-existing tracked-file hashes equal the starting snapshot; see [preservation-check.json](preservation-check.json).
- English commit title: `Freeze Wasm performance baselines and exact-output gates`.
- English PR description: [PR_DESCRIPTION.md](PR_DESCRIPTION.md).
- Publishing actions: none.
