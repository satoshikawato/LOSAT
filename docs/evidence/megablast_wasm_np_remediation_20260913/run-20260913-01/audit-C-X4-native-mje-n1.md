# Independent C-X4 Mje native n1 audit

Reviewer: /root/audit_megablast (ncbi_parity_auditor), read-only.
Disposition: raw/time/RSS PASS for this condition, no native speedup claim.

Five-sample median B1 16.70943805697607s / X4 17.252233978011645s:
3.248439% slower, delta0.542795921s within allowance0.835471903s.
Ranges are disjoint: B16.272627–16.828095s, C17.182476–17.438217s.
CPU user+system medians17.08s/17.78s; maximum RSS94,982,144/96,870,400 bytes,
increase1,888,256 bytes within16MiB.

All15 actual raw outputs (25,929,885 bytes) match the official oracle.
Records comprise1 oracle,2 diagnostic,2 warmup and10 timed searches. Inputs,
argv, artifacts, source hashes, runner snapshots, usage/results, AB/BA order,
n1 pool0, monotonic non-overlap and zero exclusions were independently checked.
Reatime diagnostic disagreements:8 total/5 timed, largest+1.045917486s.
All acceptance wall times match monotonic endpoints; BOOTTIME differs at most
7.565 microseconds. No correction or exclusion is justified.

Candidate source linking.rs SHA256:
8a503fde9a7e03e770be91ef25d79c8b4cfa99f099a2e2e69bcbc2cd2bb9883f.
Candidate native binary SHA256:
4c5cfa371e703d81dc10892740c0e2e8703f789c9480c711f95a04667793a692.
The original integrated native failure remains. Other native cells, Wasm,
reuse and whole-change acceptance are not established by this one-condition audit.
