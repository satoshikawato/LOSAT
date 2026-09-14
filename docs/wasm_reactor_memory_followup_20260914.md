# Repeated reactor memory follow-up

Status: open shared-runtime investigation, separated from the search performance
remediation by explicit user instruction on 2026-09-14. The original memory
failures remain failures under their original criteria. See the
[scope decision](evidence/megablast_wasm_np_remediation_20260913/run-20260913-01/runtime-memory-scope-decision.json).

The fixed observation window ran four N/P inputs for 16 calls each, in two
baseline/candidate sessions (256 untimed jobs). Six of sixteen per-version
series still grew during their last eight calls, including both baseline and
candidate. All raw outputs, worker lifecycles and memory budgets passed. Seven
of eight paired candidate maximum linear-memory observations were smaller than
the baseline; P132 session0 was 655,360 bytes larger. Allocated linear memory is
not the same measure as retained live allocations or process RSS.

A separately instrumented, validated allocator diagnostic ran 96 untimed jobs.
All sixteen series reached their observed post-call live maximum by call2; the
maximum and final usable bytes/block counts were equal between versions and
sessions. Fifteen series nevertheless had later linear-memory growth.

| Input | Observed retained usable bytes | Blocks |
|---|---:|---:|
| Short BLASTN | 1,419,328 | 40 |
| Large BLASTN | 6,334,528 | 40 |
| BLASTP AP027078/AP027131 | 1,616,940 | 41 |
| BLASTP AP027132/NZ_CP006932 | 2,141,132 | 41 |

These measurements are after input deallocation and worker exit and include the
retained result. Usable allocation size excludes allocator headers and free
chunks. They do not establish a bound for unlimited repetitions and are not
performance samples. The diagnostic wrapper's initial validation failures are
also preserved.

Source inspection records delayed crossbeam epoch retirement as one possible
contributor, but does not attribute the observed byte counts to a specific
allocation type. No production allocator flush, early result drop, permanent
worker pool or runtime behavior change was introduced as a workaround.

Further work should characterize allocation/retirement and fragmentation at the
same validated boundary, bind worker/TLS cleanup to actual allocation records,
and evaluate any fix against both versions with a predeclared finite observation
window. A future fix must preserve raw output, failure recovery, search-scoped
worker termination, direct-API ownership and existing time/RSS/budget limits.
Repeatedly extending a window until it happens to flatten is not a resolution.

Evidence and reproduction scripts are indexed in
[REPRODUCE.md](evidence/megablast_wasm_np_remediation_20260913/run-20260913-01/REPRODUCE.md),
with the fixed-window results under `integrated-np-reactor-highwater/`, allocator
results under `allocator-diagnostic/`, and source analysis in
`reuse-epoch-retirement-source.md` within the same run directory.
