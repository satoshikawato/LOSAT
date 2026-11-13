//! The `common` test binary: tests of the command line, the output formats and the parts several
//! programs share. Each file listed here is one module of it, named after the file, so
//! `cargo test --test common <file>::` runs one file's tests. A new file of LOSAT/tests/ gets its
//! two lines in its suite (`common` for the command line, the output formats and shared parts;
//! `tests/suite_layout.rs` checks that every file is in exactly one test binary).

#[path = "../blastdbcmd_alias_keys.rs"]
mod blastdbcmd_alias_keys;
#[path = "../blastdbcmd_batch_line_ends.rs"]
mod blastdbcmd_batch_line_ends;
#[path = "../blastdbcmd_discovery.rs"]
mod blastdbcmd_discovery;
#[path = "../blastdbcmd_info.rs"]
mod blastdbcmd_info;
#[path = "../blastdbcmd_masks.rs"]
mod blastdbcmd_masks;
#[path = "../blastdbcmd_retrieve.rs"]
mod blastdbcmd_retrieve;
#[path = "../blastdbcmd_several_names.rs"]
mod blastdbcmd_several_names;
#[path = "../blastx_library_db_hold.rs"]
mod blastx_library_db_hold;
#[path = "../cli_smoke.rs"]
mod cli_smoke;
#[path = "../db_empty_title.rs"]
mod db_empty_title;
#[path = "../db_seqidlist.rs"]
mod db_seqidlist;
#[path = "../db_taxid_root.rs"]
mod db_taxid_root;
#[path = "../default_bundle.rs"]
mod default_bundle;
#[path = "../entry_names.rs"]
mod entry_names;
#[path = "../export_search_strategy.rs"]
mod export_search_strategy;
#[path = "../library_api.rs"]
mod library_api;
#[path = "../library_examples.rs"]
mod library_examples;
#[path = "../makeblastdb_ncbi_exceptions.rs"]
mod makeblastdb_ncbi_exceptions;
#[path = "../makeblastdb_output_streams.rs"]
mod makeblastdb_output_streams;
#[path = "../makeblastdb_past_32_bits.rs"]
mod makeblastdb_past_32_bits;
#[path = "../makeblastdb_rename_lock.rs"]
mod makeblastdb_rename_lock;
#[path = "../masker_fixtures.rs"]
mod masker_fixtures;
#[path = "../ncbi_lookup_paths.rs"]
mod ncbi_lookup_paths;
#[path = "../ncbi_toolkit_dryrun.rs"]
mod ncbi_toolkit_dryrun;
#[path = "../ncbi_toolkit_file_options.rs"]
mod ncbi_toolkit_file_options;
#[path = "../ncbi_toolkit_help.rs"]
mod ncbi_toolkit_help;
#[path = "../ncbi_toolkit_version.rs"]
mod ncbi_toolkit_version;
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
#[path = "../outfmt11_titles.rs"]
mod outfmt11_titles;
#[path = "../outfmt13_14_files.rs"]
mod outfmt13_14_files;
#[path = "../outfmt13_16_titles.rs"]
mod outfmt13_16_titles;
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
#[path = "../parse_deflines_subject.rs"]
mod parse_deflines_subject;
#[path = "../redo_contain_idx.rs"]
mod redo_contain_idx;
#[path = "../report_counts_cli.rs"]
mod report_counts_cli;
#[path = "../rpsbproc_fixtures.rs"]
mod rpsbproc_fixtures;
#[path = "../search_environment.rs"]
mod search_environment;
#[path = "../seq_id_first_line.rs"]
mod seq_id_first_line;
#[path = "../seqalign_samples.rs"]
mod seqalign_samples;
#[path = "../suite_layout.rs"]
mod suite_layout;
#[path = "../tax_report_streams.rs"]
mod tax_report_streams;
#[path = "../taxonomy_read_order.rs"]
mod taxonomy_read_order;
#[path = "../web_v1_search_space.rs"]
mod web_v1_search_space;
#[path = "../x_reach_directory.rs"]
mod x_reach_directory;
