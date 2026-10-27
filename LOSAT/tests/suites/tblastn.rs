//! The `tblastn` test binary: tblastn tests. Each file listed here is one module of it, named after
//! the file, so `cargo test --test tblastn <file>::` runs one file's tests. A new file of
//! LOSAT/tests/ gets its two lines in its program's suite (`tests/suite_layout.rs` checks that
//! every file is in exactly one test binary).

#[path = "../run_local_tblastn.rs"]
mod run_local_tblastn;
#[path = "../tblastn_batch_rule.rs"]
mod tblastn_batch_rule;
#[path = "../tblastn_culling.rs"]
mod tblastn_culling;
#[path = "../tblastn_db_search.rs"]
mod tblastn_db_search;
#[path = "../tblastn_db_stream.rs"]
mod tblastn_db_stream;
#[path = "../tblastn_link_window.rs"]
mod tblastn_link_window;
#[path = "../tblastn_oid_chunk.rs"]
mod tblastn_oid_chunk;
#[path = "../tblastn_reap_idx.rs"]
mod tblastn_reap_idx;
#[path = "../tblastn_redo_lpt.rs"]
mod tblastn_redo_lpt;
#[path = "../tblastn_scan_bf.rs"]
mod tblastn_scan_bf;
#[path = "../tblastn_scan_byte.rs"]
mod tblastn_scan_byte;
#[path = "../tblastn_subject_store.rs"]
mod tblastn_subject_store;
#[path = "../tblastn_ungap.rs"]
mod tblastn_ungap;
#[path = "../tblastn_unit_xlate.rs"]
mod tblastn_unit_xlate;
#[path = "../tblastn_xlate_simd.rs"]
mod tblastn_xlate_simd;
