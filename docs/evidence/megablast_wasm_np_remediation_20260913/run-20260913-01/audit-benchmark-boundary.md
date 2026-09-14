# Independent read-only benchmark boundary audit

Reviewer: `ncbi_parity_auditor` (`/root/audit_megablast`), 2026-09-13.

The reviewer independently rehashed all five shared runner JS files, inspected actual cold argv, filesystem types, worker imports and launch paths, recalculated the path probe medians, and checked the exclusion records. Identical B1 Wasm with ext4 runners had median 0.762769 s versus 1.265441 s with DrvFS runners (one warmup plus five alternating pairs, LC738874/LC738870, n8); all twelve outputs were byte-identical. This proves a path-related timing effect for that probe, not a transferable correction factor.

The 24 initial C-N1 cold commands confirm the runner filesystem mismatch. Their comparative timings are excluded. All 16 fair-exploration cold commands use the same ext4 runner path, Node flags and normalized search argv, with both artifacts on ext4. Actual artifact hashes agree with metadata. The new `--candidate-runners` option is correctly wired; reuse still uses the existing common host and is separately labeled. The three new NCBI source-comment ranges were corrected from 177–188 to 171–188 to include `Run()` at line 173 and `Join()` at line 180.

The reviewer also checked every saved cold command for historical run-01 `p2-exploration` (24 commands), `p3-exploration` (12), and run-02 `cache-exploration`, `dp-exploration`, `matrix-exploration` (12 each). All compared baseline runners **and Wasm** on `/tmp` against candidate runners and Wasm in the repository. The recorded five JS content hashes agree. These historical comparisons therefore do not isolate candidate performance. Their original decisions remain historical records; any reconsideration requires a new B1-based candidate and fair measurements. The current probe difference must not be subtracted from old samples.

The reviewer supports the exclusions and corrected boundary. No adoption or 5% achievement claim was accepted: fair exploration has only three measured pairs, and controls/final gates remain outstanding. No source edits, builds, or measurements were performed by the reviewer.
