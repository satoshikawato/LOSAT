# Independent read-only audit

Auditor: `/root/audit_measurement`, custom `ncbi_parity_auditor`. The auditor inspected source, actual saved outputs, manifests, argv, timing/usage and lifecycle logs. It performed no edits, builds, tests or benchmarks. This is a primary-agent transcription of the findings delivered in this conversation.

## Accepted TBLASTX profile

The auditor supports the limited **cold Node/WASI comparison launch profile** claim. All 43 fresh records bind actual raw bytes, command, result, GNU CPU/RSS and checksums. Every output equals the current oracle's 351,051 bytes. There are no excluded records. The 6 diagnostics prove requested pool counts and each tid's spawn → ready → exit(0) sequence. Detailed worker events are present in diagnostic runs; ordinary cold timing rows do not independently record those events.

Warmup order is AB; measured repetitions are BA/AB/BA/AB/BA on each axis. There are 6 warmup and 30 measured rows. The independently recomputed median reductions are serial n1 15.0886%, threaded n1 14.5351%, threaded n8 32.3135%. The differences between each condition's maximum measured RSS are −14.0000, +10.6914, +27.9492 MiB. All satisfy the declared policy.

All 43 fresh monotonic/boottime pairs agree; the largest difference is 5.16 µs. Seven realtime disagreements remain diagnostics under the already declared CLOCK_MONOTONIC policy. They were not corrected into physical wall time.

For the two preserved p1-controls series, the auditor rechecked all 30 records' raw bytes and usage plus four diagnostics. MjeNMV/MelaMJNV is 39.767963 → 14.364419 s, RSS +47.328125 MiB; AP027280 self is 70.531921 → 25.571660 s, RSS +36.8515625 MiB. Both pass the newly authorized max(10%,48 MiB) increase. Mje has only 0.671875 MiB observed allowance remaining. This policy is an acceptance criterion, not a runtime memory cap.

Legacy rows have no boottime crosscheck. The original CLOCK_MONOTONIC measurement boundary and prior clock-audit limitations remain; these are reassessed old measurements, not newly repeated controls. p1-controls as a whole remains PARTIAL because of its separate failed input.

Node executable / 24.21.0 / V8, flags, search artifact, production hashes, JS execution runners and input identities agree across old/fresh runs. The old metadata's absent Cargo.toml hash is supplemented by the frozen baseline source archive, whose hash matches the fresh recording. NCBI is 2.17.0; these TBLASTX cases use code 1, with no parity exception applied. All six source bindings in tblastx-evaluation.json also agree.

All 12 real simple searches bind actual argv, executable/log/output hashes and exact current oracle bytes. Only the two Wasm invocations in the TBLASTX-default case receive the dedicated flags. Opt-out, BLASTP and BLASTN receive none. The plot contains the three successful BLASTP conditions. Final Python log: 39 tests / OK.

## Argument selection and source provenance

`comparison_data.node_args(program)` appends TBLASTX-only arguments after common arguments. An explicit empty TBLASTX array removes only dedicated additions. Invalid JSON, non-string arguments and NUL fail explicitly. The manifest and NUL files use the same function, and shell quoted arrays preserve argument boundaries. Each search selects its array afresh; megablast maps to blastn and cannot receive the TBLASTX profile by default. Native, oracle, search argv and existing JS runners are unchanged.

The original run-02 start archive omitted shell and documentation. The shell baseline is independently bound to run-01/final-task-source.tar.gz, SHA-256 `193bd1857acb8e1f29fb808cf342aaf9fe95436264d679e6babfc7046b1db94c`, also matching the run-01 real comparison's runner hash. Applying the four exact known tblastx_profile.py replacements produces the implementation shell byte-for-byte. A separately labelled supplement is appropriate; the original archive must not be relabelled as having contained it. snapshot-v2 records this limitation.

## Rejected Rust candidates

### Subject SEG cache folder reuse

The source audit supports the candidate's bounded semantics on the covered BLASTP query-parallel CLI path: cache keys select fixed subject OID/range, the query-dependent near-identical/bias/CBS checks remain before lookup on each call, and cached bytes/bias are cloned. No query-adjusted matrix is cached. Full-subject windows have frame 0; query heaps replay in query index order. Cache lifetime is search-folder-local, not a guaranteed worker lifetime. The pointer used by synchronous postprocessing is restored even on errors.

NCBI reference: blast_kappa.c:1414–1454 (fixed SEG parameters 10, 1.8, 2.1), 1626–1643 (query-dependent masking gate), 3493–3503 (thread-local state); redo_alignment.c:800–806 (whole subject window). Current CLI supports BLOSUM62, gaps 11/1 and CBS mode 2 on this path; unsupported alternate combinations remain fail-fast. Query `-seg no` does not disable Kappa subject SEG.

The measured cache hits increased but whole-process time regressed. The candidate is rejected, so additional adoption unit tests (masked/near-identical transitions, false bias, empty SEG, clone isolation, different OIDs) and the broader native/serial matrix are not claimed as run.

### Wasm DP compilation boundary

Only two wasm32 `inline(never)` attributes differ at score-only and traceback implementations. Arguments, recurrence, tie-breaks and loops do not change. NCBI blast_gapalign.c:543–578 and 841–866 describe the corresponding DP work. The source attribute does not prove final V8 inlining, machine-code size or a sampled CPU share. This affects Wasm targets but not native; only isolated threaded exploration was measured. Rejected after time regressions.

### Adjusted-matrix row specialization

The auditor checked all 18 matrix/direction dispatch arms: mode 0 BLOSUM62, mode 1 adjusted, mode 2 standard fallback. Reverse flags and arguments match the original code. Loop bodies and cell semantics remain the same; adjusted mode's `unreachable` is not reachable from the checked wrappers. No matrix copy or lifetime extension is added. The candidate is rejected at at most 5.45% improvement.

Reference-comment errata in the archived rejected version: actual NCBI blast_gapalign.c cell sites are 578, 866, 1115 and 1180, not 578/859/1099. The restricted wrapper should reference 1078–1115 and 1180. Correct these in a new source version before reviving the candidate; the measured source and manifests remain immutable.

If revived, useful missing tests include wrapper versus generic adjusted/standard dispatch with asymmetric matrices, non-palindromic unequal sequence lengths, non-zero bases, restricted boundaries 9/10/11 and 19/20/21, fresh/reused scratch, and traceback script identity. Those broader gates were not run for this rejected candidate.

The audit does not establish browser speed, warm-reuse gains, a faster Rust search kernel, or complete release parity. The existing megablast and formal certification limitations remain open.
