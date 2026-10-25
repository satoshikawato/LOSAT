//! The `blastx` test binary: blastx tests. Each file listed here is one module of it, named after
//! the file, so `cargo test --test blastx <file>::` runs one file's tests. A new file of
//! LOSAT/tests/ gets its two lines in its program's suite (`tests/suite_layout.rs` checks that
//! every file is in exactly one test binary).

#[path = "../blastx_db_search.rs"]
mod blastx_db_search;
#[path = "../blastx_db_stream.rs"]
mod blastx_db_stream;
#[path = "../blastx_diag.rs"]
mod blastx_diag;
#[path = "../blastx_diag_bitfield.rs"]
mod blastx_diag_bitfield;
#[path = "../blastx_flow.rs"]
mod blastx_flow;
#[path = "../blastx_hitbf.rs"]
mod blastx_hitbf;
#[path = "../blastx_scan_byte.rs"]
mod blastx_scan_byte;
#[path = "../blastx_subject_flow.rs"]
mod blastx_subject_flow;
#[path = "../blastx_ungap_sent.rs"]
mod blastx_ungap_sent;
#[path = "../blastx_unit_cap.rs"]
mod blastx_unit_cap;
