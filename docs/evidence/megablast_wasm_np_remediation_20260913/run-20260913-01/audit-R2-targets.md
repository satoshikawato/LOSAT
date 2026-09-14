# Independent read-only R2 follow-up

Auditor: ncbi_parity_auditor `/root/audit_megablast`. No builds or benchmarks were run during this review.

Supported scope: all 30 saved target records have matching raw file size, SHA256 and equality flags. All 18 B1 target outputs (three cases × native8, serial1, threaded1/2/4/8) match the official oracle. Actual hashes of six B0/B1 artifacts match records. All 14 existing task raw pairs also match, including outfmt6/7, multiple queries and no-hit. In all 18 threaded logs, per-worker spawn attempt → spawned → ready → exited order and counts are consistent. Worker lifecycle is not utilization or performance evidence. B0 PASS means execution success; its NZ/Sakai oracle mismatches remain explicit.

The frozen comparison has 6,482 lines on both sides and differs in one lexical input-path header and five data rows. Hashes match the record. This diagnostic is not formal Gate A/B certification; fixing the path cannot remove the five data differences. Three existing release parity manifests retain their initial hashes.

The new integration test has valid static public CLI → run wiring and is auto-discovered by Cargo. Its new 38,401-byte golden equals the official NZ oracle. Compilation/execution remained pending at the time of review. The initial option/default source rationale is supported. No new implementation blocker was found. Claims do not extend to outfmt0/custom, all boundaries, formal release certification or speed targets.
