//! The `tblastx` test binary: tblastx tests. Each file listed here is one module of it, named after
//! the file, so `cargo test --test tblastx <file>::` runs one file's tests. A new file of
//! LOSAT/tests/ gets its two lines in its program's suite (`tests/suite_layout.rs` checks that
//! every file is in exactly one test binary).

#[path = "../run_local_tblastx.rs"]
mod run_local_tblastx;
#[path = "../tblastx_band.rs"]
mod tblastx_band;
#[path = "../tblastx_batch_pipe.rs"]
mod tblastx_batch_pipe;
#[path = "../tblastx_db_search.rs"]
mod tblastx_db_search;
#[path = "../tblastx_extbatch.rs"]
mod tblastx_extbatch;
#[path = "../tblastx_hitbf.rs"]
mod tblastx_hitbf;
#[path = "../tblastx_kmerge.rs"]
mod tblastx_kmerge;
#[path = "../tblastx_link_fast.rs"]
mod tblastx_link_fast;
#[path = "../tblastx_link_slices.rs"]
mod tblastx_link_slices;
#[path = "../tblastx_link_sort.rs"]
mod tblastx_link_sort;
#[path = "../tblastx_matrix.rs"]
mod tblastx_matrix;
#[path = "../tblastx_prefetch.rs"]
mod tblastx_prefetch;
#[path = "../tblastx_qslice.rs"]
mod tblastx_qslice;
#[path = "../tblastx_query_frames.rs"]
mod tblastx_query_frames;
#[path = "../tblastx_ref_lut.rs"]
mod tblastx_ref_lut;
#[path = "../tblastx_report_batch.rs"]
mod tblastx_report_batch;
#[path = "../tblastx_report_par.rs"]
mod tblastx_report_par;
#[path = "../tblastx_resolve.rs"]
mod tblastx_resolve;
#[path = "../tblastx_search_timing.rs"]
mod tblastx_search_timing;
#[path = "../tblastx_seed_bucket.rs"]
mod tblastx_seed_bucket;
#[path = "../tblastx_seed_compact.rs"]
mod tblastx_seed_compact;
#[path = "../tblastx_seed_release.rs"]
mod tblastx_seed_release;
#[path = "../tblastx_seed_unmap.rs"]
mod tblastx_seed_unmap;
#[path = "../tblastx_store.rs"]
mod tblastx_store;
#[path = "../tblastx_subject_next.rs"]
mod tblastx_subject_next;
#[path = "../tblastx_ungap_sent.rs"]
mod tblastx_ungap_sent;
#[path = "../tblastx_unit_gate.rs"]
mod tblastx_unit_gate;
