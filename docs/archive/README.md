# docs/archive

Plan documents of closed lines, moved here from the `docs/` root on 2026-10-23 (slim-20261022 P2b; the
Owner's decision of 2026-10-22: old plans move to `docs/archive/`, not deleted). The files are unchanged except
that their relative links were rewritten for the new location. Run records under `docs/evidence/` that name the
old paths were left as they are.

| Document | Line |
| --- | --- |
| [blastn_bit_parity_fix_plan.md](blastn_bit_parity_fix_plan.md) | BLASTN parity, 2026-05 |
| [blastn_bitwise_parity_implementation_plan.md](blastn_bitwise_parity_implementation_plan.md) | BLASTN parity, 2026-05 to 2026-08 |
| [blastn_edl933_sakai_megablast_extra_hsps_fix_plan.md](blastn_edl933_sakai_megablast_extra_hsps_fix_plan.md) | BLASTN parity, 2026-05 |
| [blastn_hit_score_divergence_fix_plan.md](blastn_hit_score_divergence_fix_plan.md) | BLASTN parity, 2026-05 |
| [blastn_traceback_prune_performance_plan.md](blastn_traceback_prune_performance_plan.md) | BLASTN speed, 2026-05 |
| [blastn_v010_supported_parity_plan.md](blastn_v010_supported_parity_plan.md) | v0.1.0 BLASTN scope, 2026-08 |
| [blastn_visual_divergence_fix_plan.md](blastn_visual_divergence_fix_plan.md) | BLASTN parity, 2026-05 |
| [blastp_bit_parity_implementation_plan.md](blastp_bit_parity_implementation_plan.md) | BLASTP parity, 2026-05 |
| [blastp_performance_analysis.md](blastp_performance_analysis.md) | BLASTP speed analysis, 2026-05 |
| [blastp_wasm_current_thread_pool_plan.md](blastp_wasm_current_thread_pool_plan.md) | Wasm threads (superseded 2026-09-13) |
| [megablast_parity_wasm_blastn_blastp_remediation_plan_20260913.md](megablast_parity_wasm_blastn_blastp_remediation_plan_20260913.md) | megablast + Wasm, 2026-09-13 |
| [tblastx_native_ncbi_parity_plan.md](tblastx_native_ncbi_parity_plan.md) | TBLASTX parity, 2026-05 |
| [tblastx_wasm_parity_plan.md](tblastx_wasm_parity_plan.md) | TBLASTX Wasm parity, 2026-05 |
| [v0.1.0_release_readiness_plan.md](v0.1.0_release_readiness_plan.md) | v0.1.0 release (published 2026-09-17) |
| [v0.1.0_scope.md](v0.1.0_scope.md) | v0.1.0 release scope |
| [wasm_four_program_optimization_plan_20260915.md](wasm_four_program_optimization_plan_20260915.md) | Wasm speed, 2026-09-15 |
| [wasm_multithreading_improvement_plan.md](wasm_multithreading_improvement_plan.md) | Wasm threads (superseded 2026-09-13) |
| [wasm_performance_optimization_plan_20260913.md](wasm_performance_optimization_plan_20260913.md) | Wasm speed, 2026-09-13 |
| [wasm_remaining_work_implementation_plan_20260915.md](wasm_remaining_work_implementation_plan_20260915.md) | Wasm speed, 2026-09-15 |
| [wasm_threading_audit_20260913.md](wasm_threading_audit_20260913.md) | Wasm threads audit, 2026-09-13 |
| [wasm_threading_remediation_plan_20260913.md](wasm_threading_remediation_plan_20260913.md) | Wasm threads, 2026-09-13 |
| [wasm_threading_remediation_results_20260913.md](wasm_threading_remediation_results_20260913.md) | Wasm threads, 2026-09-13 |
| [wasm_total_threads_20260914.md](wasm_total_threads_20260914.md) | Wasm thread count (N includes the caller), 2026-09-14 |
| [losat_wasm_np_plan_20260906/](losat_wasm_np_plan_20260906/README.md) | Wasm BLASTN/BLASTP speed sessions, 2026-09-06 |

Plans that an open line still follows stay at the `docs/` root: the LOSAT Web plan and its sessions, the reactor
memory follow-up, the LOSATX (BLASTX) and TLOSAN (TBLASTN) v0.2.0 plans with their sessions (v0.2.0 is not
released), and the user documents `cli_v2_migration.md` and `losatx_blastx_native_cli.md`. When such a line
closes, its plan moves here the same way: `git mv`, then rewrite the links to and from it.
