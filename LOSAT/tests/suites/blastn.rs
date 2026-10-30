//! The `blastn` test binary: blastn tests. Each file listed here is one module of it, named after
//! the file, so `cargo test --test blastn <file>::` runs one file's tests. A new file of
//! LOSAT/tests/ gets its two lines in its program's suite (`tests/suite_layout.rs` checks that
//! every file is in exactly one test binary).

#[path = "../blastn_bn_csr.rs"]
mod blastn_bn_csr;
#[path = "../blastn_bn_seed_lists.rs"]
mod blastn_bn_seed_lists;
#[path = "../blastn_db_cut.rs"]
mod blastn_db_cut;
#[path = "../blastn_db_search.rs"]
mod blastn_db_search;
#[path = "../blastn_db_stream.rs"]
mod blastn_db_stream;
#[path = "../blastn_dc_seed_lean.rs"]
mod blastn_dc_seed_lean;
#[path = "../blastn_diag_lv.rs"]
mod blastn_diag_lv;
#[path = "../blastn_disc_lookup_par.rs"]
mod blastn_disc_lookup_par;
#[path = "../blastn_few_subj.rs"]
mod blastn_few_subj;
#[path = "../blastn_greedy_row_fuse.rs"]
mod blastn_greedy_row_fuse;
#[path = "../blastn_hit_group.rs"]
mod blastn_hit_group;
#[path = "../blastn_hl_grow.rs"]
mod blastn_hl_grow;
#[path = "../blastn_lut_drop.rs"]
mod blastn_lut_drop;
#[path = "../blastn_mb_csr.rs"]
mod blastn_mb_csr;
#[path = "../blastn_mb_fill.rs"]
mod blastn_mb_fill;
#[path = "../blastn_mb_lookup_once.rs"]
mod blastn_mb_lookup_once;
#[path = "../blastn_mb_lut_keep.rs"]
mod blastn_mb_lut_keep;
#[path = "../blastn_pairpar_disc.rs"]
mod blastn_pairpar_disc;
#[path = "../blastn_pairpar_disc_arr.rs"]
mod blastn_pairpar_disc_arr;
#[path = "../blastn_pairpar_fold.rs"]
mod blastn_pairpar_fold;
#[path = "../blastn_pairpar_merge.rs"]
mod blastn_pairpar_merge;
#[path = "../blastn_pairpar_na.rs"]
mod blastn_pairpar_na;
#[path = "../blastn_qctx_sentinel.rs"]
mod blastn_qctx_sentinel;
#[path = "../blastn_qgroup_par.rs"]
mod blastn_qgroup_par;
#[path = "../blastn_query_drop.rs"]
mod blastn_query_drop;
#[path = "../blastn_query_dust_par.rs"]
mod blastn_query_dust_par;
#[path = "../blastn_query_prep_par.rs"]
mod blastn_query_prep_par;
#[path = "../blastn_query_scan_par.rs"]
mod blastn_query_scan_par;
#[path = "../blastn_query_stream.rs"]
mod blastn_query_stream;
#[path = "../blastn_report_stream.rs"]
mod blastn_report_stream;
#[path = "../blastn_strand.rs"]
mod blastn_strand;
#[path = "../blastn_subject_store.rs"]
mod blastn_subject_store;
#[path = "../blastn_subject_stream.rs"]
mod blastn_subject_stream;
#[path = "../blastn_tab_par.rs"]
mod blastn_tab_par;
#[path = "../blastn_ungap_count.rs"]
mod blastn_ungap_count;
#[path = "../blastn_window_masker.rs"]
mod blastn_window_masker;
#[path = "../megablast_xdrop.rs"]
mod megablast_xdrop;
#[path = "../run_local_blastn.rs"]
mod run_local_blastn;
