//! The `common` test binary: tests of the command line, the output formats and the parts several
//! programs share. Each file listed here is one module of it, named after the file, so
//! `cargo test --test common <file>::` runs one file's tests. A new file of LOSAT/tests/ gets its
//! two lines in its suite (`common` for the command line, the output formats and shared parts;
//! `tests/suite_layout.rs` checks that every file is in exactly one test binary).

#[path = "../blastdbcmd_info.rs"]
mod blastdbcmd_info;
#[path = "../cli_smoke.rs"]
mod cli_smoke;
#[path = "../cli_v2.rs"]
mod cli_v2;
#[path = "../compressed_lookup.rs"]
mod compressed_lookup;
#[path = "../default_bundle.rs"]
mod default_bundle;
#[path = "../newton_lanes8_default.rs"]
mod newton_lanes8_default;
#[path = "../outfmt11_blastn.rs"]
mod outfmt11_blastn;
#[path = "../outfmt11_blastp.rs"]
mod outfmt11_blastp;
#[path = "../outfmt11_blastx.rs"]
mod outfmt11_blastx;
#[path = "../outfmt11_tblastn.rs"]
mod outfmt11_tblastn;
#[path = "../outfmt11_tblastx.rs"]
mod outfmt11_tblastx;
#[path = "../outfmt13_14_files.rs"]
mod outfmt13_14_files;
#[path = "../outfmt16_blastn.rs"]
mod outfmt16_blastn;
#[path = "../outfmt16_blastp.rs"]
mod outfmt16_blastp;
#[path = "../outfmt16_blastx.rs"]
mod outfmt16_blastx;
#[path = "../outfmt16_tblastn.rs"]
mod outfmt16_tblastn;
#[path = "../outfmt16_tblastx.rs"]
mod outfmt16_tblastx;
#[path = "../outfmt17_blastn.rs"]
mod outfmt17_blastn;
#[path = "../outfmt5_blastx.rs"]
mod outfmt5_blastx;
#[path = "../outfmt5_no_query.rs"]
mod outfmt5_no_query;
#[path = "../outfmt5_refusals.rs"]
mod outfmt5_refusals;
#[path = "../outfmt5_tblastn.rs"]
mod outfmt5_tblastn;
#[path = "../outfmt5_tblastx.rs"]
mod outfmt5_tblastx;
#[path = "../pairwise_par.rs"]
mod pairwise_par;
#[path = "../redo_contain_idx.rs"]
mod redo_contain_idx;
#[path = "../report_counts_cli.rs"]
mod report_counts_cli;
#[path = "../seqalign_samples.rs"]
mod seqalign_samples;
#[path = "../suite_layout.rs"]
mod suite_layout;
