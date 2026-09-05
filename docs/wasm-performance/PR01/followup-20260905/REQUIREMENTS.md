# PR01 acceptance and PR02 dependency boundaries

The PR01 baseline must be accepted before PR02 implementation relies on it.
An unfinished search is missing coverage, not a byte mismatch. A larger
timeout does not relax the output contract or remove an eligible fixture.

| Evidence | PR01 / PR02 boundary | Required next action |
| --- | --- | --- |
| Current 43 native / 41 applicable serial catalog, including the existing exception validators | Required for complete PR01 correctness; unresolved affected paths cannot enter optimization | Reconcile retained successful executions with the separately budgeted retries; retain every unresolved row |
| p11 AvCLPV/PsCLPV eligible cold matrix | Required PR01 baseline; cannot drop this case because it is slow | After its oracle passes, obtain native1/2/4/8, serial, threaded1/2/4/8 with warmup1 + timed5 + diagnostic1 |
| p11 serial warm Module / new instance | Required to complete the already declared PR01 warm matrix | Obtain warmup1 + timed5 after oracle completion; this is not a resident engine |
| Native effective computation threads | Protocol requests observation; null with a reason is permitted when unavailable | Preserve requested-N labels. Runtime-observed effective N is required before claiming same-effective-N ratios or efficiency; source-derived pool capacity is not active utilization |
| Native active workers, allocations and reset/copy cost | Not a universal PR01 acceptance or PR02-start blocker; required for claims that depend on these measurements | Collect during the corresponding scheduler/workspace optimization, outside baseline timing samples |
| Wasm BLASTP stage attribution | PR01 needs evidence for hotspot selection | New raw-gated CPU profile supplies scoped redo attribution. Exact engine timer shares remain unavailable; instrumentation/reuse work belongs to PR07 |
| Finer TBLASTX per-call linking attribution | Existing raw-gated p03 CPU profile supports the PR10 investigation priority | Isolate initial linking, relinking and inner loops in PR10 before claiming a particular optimization |
| Useful jobs, pool size, configured effective threads, helper spawns, caps and serial reason | PR02 must measure these before accepting its policy | PR02's own 0/1/2/many jobs, requests1/2/4/8, caps1/2/8 and threshold boundaries; do not implement them in PR01 follow-up |
| Threaded warm, resident API, pool lifetime | PR03 dependency and lifecycle proof | Keep NOT_RUN; serial new-instance success does not transfer certification |
| Real browser, transfer/worker-ready/teardown timing | PR12, with host source and real browsers | Keep separate from Node results |
| Windows/macOS registered native-NCBI Gate B | Release/PR13 platform gate, outside this Linux baseline | Run on each registered platform; never change LOSAT expected bytes to native fingerprints |

NCBI references read during follow-up: `blast_args.cpp:3205–3228` caps the
requested value and then applies the local-subject rule;
`prelim_stage.cpp:82–88` configures preliminary threads;
`blast_engine.c:1410–1455` owns subject iteration and per-subject setup.
These references explain thread configuration, not the observed number of
busy LOSAT workers.

Current LOSAT source gives a limited configuration inference: native N/X use
the nonzero request as the local pool capacity when the parallel path is used
(`blastn/blast_engine/run.rs:4638–4652,4887–4903`,
`tblastx/blast_engine/run_impl.rs:879–898,1230–1235,138–157`). Native BLASTP
uses `min(requested, num_cpus::get().max(1))` and per-stage admission
(`blastp/blast_engine.rs:1008–1051,1310–1356`). Native pool builders do not
use the Wasm-only `use_current_thread` configuration. None of these source
observations retroactively fills the runtime-observed fields in the original
manifest or proves concurrent useful work. The harness's cap of eight limits
the selected matrix; it is not evidence of an identical native engine cap.

No performance improvement or complete eligible-set target attainment is
required of a baseline. Complete, reproducible correctness and the declared
baseline coverage are required. PR01 NOT_ACCEPTED therefore blocks automatic
PR02 implementation, while read-only PR02 planning can use the scoped evidence.
