//! The `blastp` test binary: blastp tests. Each file listed here is one module of it, named after
//! the file, so `cargo test --test blastp <file>::` runs one file's tests. A new file of
//! LOSAT/tests/ gets its two lines in its program's suite (`tests/suite_layout.rs` checks that
//! every file is in exactly one test binary).

#[path = "../blastp_culling.rs"]
mod blastp_culling;
#[path = "../blastp_db_cut.rs"]
mod blastp_db_cut;
#[path = "../blastp_diag_batch.rs"]
mod blastp_diag_batch;
#[path = "../blastp_diag_defer.rs"]
mod blastp_diag_defer;
#[path = "../blastp_fast_fill_buckets.rs"]
mod blastp_fast_fill_buckets;
#[path = "../blastp_fast_scan_bf.rs"]
mod blastp_fast_scan_bf;
#[path = "../blastp_fast_seed_bf.rs"]
mod blastp_fast_seed_bf;
#[path = "../blastp_fast_stage.rs"]
mod blastp_fast_stage;
#[path = "../blastp_hitbf.rs"]
mod blastp_hitbf;
#[path = "../blastp_hsp_tie_order.rs"]
mod blastp_hsp_tie_order;
#[path = "../blastp_ih_pool.rs"]
mod blastp_ih_pool;
#[path = "../blastp_listed_lean.rs"]
mod blastp_listed_lean;
#[path = "../blastp_lookup_ahead.rs"]
mod blastp_lookup_ahead;
#[path = "../blastp_lookup_flat.rs"]
mod blastp_lookup_flat;
#[path = "../blastp_prelim_compact.rs"]
mod blastp_prelim_compact;
#[path = "../blastp_prelim_free.rs"]
mod blastp_prelim_free;
#[path = "../blastp_query_blocks.rs"]
mod blastp_query_blocks;
#[path = "../blastp_query_blocks_fmt.rs"]
mod blastp_query_blocks_fmt;
#[path = "../blastp_redo_db_subjects.rs"]
mod blastp_redo_db_subjects;
#[path = "../blastp_redo_fit.rs"]
mod blastp_redo_fit;
#[path = "../blastp_redo_memory.rs"]
mod blastp_redo_memory;
#[path = "../blastp_redo_par.rs"]
mod blastp_redo_par;
#[path = "../blastp_restricted_gapped.rs"]
mod blastp_restricted_gapped;
#[path = "../blastp_subj_meta.rs"]
mod blastp_subj_meta;
#[path = "../blastp_subject_runs.rs"]
mod blastp_subject_runs;
#[path = "../blastp_ungap_sent.rs"]
mod blastp_ungap_sent;
#[path = "../run_local_blastp.rs"]
mod run_local_blastp;
#[path = "../run_local_blastp_blocks.rs"]
mod run_local_blastp_blocks;
