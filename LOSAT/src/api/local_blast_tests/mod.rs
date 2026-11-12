//! The shared entries `run_local_*` of `api/local_blast.rs` in process: what an entry does by
//! itself (its sinks, its observer, its refusals of hand-set arguments, its errors), and the
//! hit records and ranges that it gives. The tests that compare an entry's reports with the
//! `losat` binary are integration tests and go through `LOSAT::api` (`LOSAT/tests/run_local_*.rs`).

mod blastn;
mod blastp;
mod blastx;
mod outfmt13_14;
mod support;
mod tblastn;
mod tblastx;
