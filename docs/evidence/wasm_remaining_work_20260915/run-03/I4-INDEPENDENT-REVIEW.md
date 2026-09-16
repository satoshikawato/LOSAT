# Independent read-only review of I4

Final audit closed: 2026-09-16 06:08:35 UTC. Reviewer: ncbi_parity_auditor.
No material blocker was found in the requested scope. No production changes,
builds, tests or timing runs were performed by the reviewer.

- Rehashed all 79 distinct command raw outputs plus four long-code4 outputs;
  compared oracle bytes, recorded output hashes and current I4 artifact hashes.
- Confirmed long oracle argv uses -db, -query_gencode 4, -db_gencode 4,
  with no -subject argument. This is the saved run-01 DB oracle, SHA256
  1da3bfe7d96d50a96726a44e2df98f5e81c05bee9cd655ad51bb8d44aa377963.
- Recomputed 58 additional metrics: Native12, serial3, serial AP3,
  reactor24, browser16. All original guards pass and sample arrays agree with
  saved run records. The previous bounded review recomputed cold/module36
  metrics and confirmed the n8/n1 reversal on both inputs.
- All seven environment records bind their policy, monitor SHA and referring
  study; return code and observed environment status pass.
- Matched artifacts for all 16 case/thread/version/session records in each of
  module, reactor and browser studies. Rehashed one timed raw output per record,
  48 spot checks total. All matched. Module/reactor metadata for all 128 samples
  reports raw/thread PASS. Both browser thread counts used threaded-command.

Limits: the final pass did not rehash every reuse/browser raw output, redo the
source/state/served-assets review, or validate root integration. Earlier bounded
reviews covered macro hygiene, scan binding, group-state manifests, selected
emission and browser manifests. Root integration is separately recorded in
I4-integration.json and I4-integrated-checks.json, with all 158 build inputs
matching the tested source and focused tests/clippy/format passing.

The evidence supports the declared fixtures/runtimes under observed normal
desktop background. CPU accounting can miss new/exited processes, workload-name
detection is heuristic, and I/O/cache/frequency effects are not excluded.
It is not proof of strict isolation or broad release certification.
