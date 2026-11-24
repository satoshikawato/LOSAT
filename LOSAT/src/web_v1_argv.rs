//! The words of a request of LOSAT Web's ABI v1 (frozen, plan TD-1), as its BLASTX and TBLASTN
//! read them: with the reader that LOSAT's command line had before it read its arguments with
//! NCBI's argument reader (`cli::reader`). v1 keeps its own copy of that reader, of its command
//! (a builder of clap's command that the command line's derive made, `frozen_command`) and of its
//! texts (line RG DECISIONS K4i-2 D5, K4i-3-1 Q4), so that its frozen messages stay
//! (`tests/fixtures/web_v1_tblastn`); the words that it takes are read into the program's options
//! as the command line reads them (`cli::reader::options`). No other caller reads an argv here.
//!
//! No NCBI counterpart: NCBI has no Web engine.

use std::ffi::OsString;
use std::path::PathBuf;

use crate::blastinput::value_parsers::{genetic_code_number, ncbi_double, ncbi_int8, ncbi_integer};
use crate::cli::Commands;
use crate::utils::genetic_code::GeneticCode;
use clap::builder::{PossibleValuesParser, ValueParser};
use clap::error::ErrorKind;

/// The command of v1's `program` (BLASTX or TBLASTN): the command line's when v1 was frozen
/// (`frozen_command`), with the value parsers and the exclusions of the arguments that v1's reader
/// had (`v1_parsers`, `V1_EXCLUSIONS`): the command line's arguments are NCBI's reader's since line
/// RG item K4i-2b.
pub(crate) fn command(program: &str) -> clap::Command {
    v1_command(frozen_command(program), program)
}

/// `command` of v1's program `program` on the program's command `command`.
fn v1_command(command: clap::Command, program: &str) -> clap::Command {
    command.mut_subcommand(program, |mut scope| {
        for (name, parser) in v1_parsers(program) {
            scope = scope.mut_arg(name, |arg| arg.value_parser(parser));
        }
        for (name, excluded) in V1_EXCLUSIONS {
            if scope.get_arguments().any(|arg| arg.get_id() == *name) {
                scope = scope.mut_arg(*name, |arg| arg.conflicts_with_all(*excluded));
            }
        }
        if program == "tblastn" {
            scope = scope.mut_arg("in_pssm", |arg| {
                arg.conflicts_with_all(["query", "query_loc"])
            });
        }
        scope.mut_arg("mt_mode", |arg| arg.requires("num_threads"))
    })
}

/// How an argument of v1's frozen command reads its value: the value parser that clap's derive
/// gave its field when v1 was frozen (`blastinput::value_parsers`).
#[derive(Debug, Clone, Copy)]
enum Reads {
    /// A flag (`ArgAction::SetTrue`).
    Flag,
    /// Text (clap's `String`).
    Text,
    /// A file name's bytes (`ncbi_path`).
    Path,
    /// `ncbi_integer`.
    Integer,
    /// `ncbi_int8`.
    Int8,
    /// `ncbi_double`.
    Double,
    /// `ncbi_boolean`.
    Boolean,
    /// `ncbi_integer_as::<usize>`.
    Count,
    /// `ncbi_integer_as::<u8>`.
    Small,
}

use Reads::{Boolean, Count, Double, Flag, Int8, Integer, Path, Small, Text};

/// An argument of v1's frozen command: its id, its long name, its value name, how it reads its
/// value, and its default.
type Frozen = (
    &'static str,
    &'static str,
    &'static str,
    Reads,
    Option<&'static str>,
);

/// BLASTX's arguments in v1's command, in the order of clap's derive of `BlastxArgs` when v1 was
/// frozen (line RG item K4i-3c; the derive and its help texts left the command line there).
#[rustfmt::skip]
const BLASTX: &[Frozen] = &[
    ("query", "query", "PATH", Path, Some("-")),
    ("subject", "subject", "PATH", Path, None),
    ("subject_loc", "subject_loc", "RANGE", Text, None),
    ("db", "db", "NAME", Text, None),
    ("seqidlist", "seqidlist", "filename", Text, None),
    ("negative_seqidlist", "negative_seqidlist", "filename", Text, None),
    ("taxids", "taxids", "taxids", Text, None),
    ("taxidlist", "taxidlist", "filename", Text, None),
    ("negative_taxids", "negative_taxids", "taxids", Text, None),
    ("negative_taxidlist", "negative_taxidlist", "filename", Text, None),
    ("no_taxid_expansion", "no_taxid_expansion", "NO_TAXID_EXPANSION", Flag, None),
    ("db_soft_mask", "db_soft_mask", "filtering_algorithm", Text, None),
    ("db_hard_mask", "db_hard_mask", "filtering_algorithm", Text, None),
    ("dbsize", "dbsize", "LETTERS", Int8, None),
    ("show_gis", "show_gis", "SHOW_GIS", Flag, Some("false")),
    ("html", "html", "HTML", Flag, Some("false")),
    ("parse_deflines", "parse_deflines", "PARSE_DEFLINES", Flag, Some("false")),
    ("searchsp", "searchsp", "INT_VALUE", Int8, None),
    ("out", "out", "PATH", Path, None),
    ("export_search_strategy", "export_search_strategy", "filename", Path, None),
    ("task", "task", "TASK", Text, Some("blastx")),
    ("evalue", "evalue", "EVALUE", Double, Some("10")),
    ("word_size", "word_size", "WORD_SIZE", Integer, None),
    ("max_hsps", "max_hsps", "MAX_HSPS", Integer, None),
    ("qcov_hsp_perc", "qcov_hsp_perc", "QCOV_HSP_PERC", Double, None),
    ("sum_stats", "sum_stats", "SUM_STATS", Boolean, None),
    ("gap_open", "gapopen", "GAP_OPEN", Integer, None),
    ("gap_extend", "gapextend", "GAP_EXTEND", Integer, None),
    ("max_intron_length", "max_intron_length", "MAX_INTRON_LENGTH", Integer, Some("0")),
    ("seg", "seg", "SEG", Text, Some("12 2.2 2.5")),
    ("soft_masking", "soft_masking", "SOFT_MASKING", Boolean, Some("false")),
    ("matrix", "matrix", "MATRIX", Text, None),
    ("threshold", "threshold", "THRESHOLD", Double, None),
    ("window_size", "window_size", "WINDOW_SIZE", Integer, None),
    ("xdrop_ungap", "xdrop_ungap", "float_value", Double, None),
    ("xdrop_gap", "xdrop_gap", "float_value", Double, None),
    ("xdrop_gap_final", "xdrop_gap_final", "float_value", Double, None),
    ("culling_limit", "culling_limit", "CULLING_LIMIT", Integer, None),
    ("best_hit_overhang", "best_hit_overhang", "float_value", Double, None),
    ("best_hit_score_edge", "best_hit_score_edge", "float_value", Double, None),
    ("subject_besthit", "subject_besthit", "SUBJECT_BESTHIT", Flag, None),
    ("ungapped", "ungapped", "UNGAPPED", Flag, None),
    ("lcase_masking", "lcase_masking", "LCASE_MASKING", Flag, None),
    ("use_sw_tback", "use_sw_tback", "USE_SW_TBACK", Flag, None),
    ("strand", "strand", "STRAND", Text, Some("both")),
    ("query_loc", "query_loc", "RANGE", Text, None),
    ("query_gencode", "query_gencode", "QUERY_GENCODE", Small, Some("1")),
    ("outfmt", "outfmt", "OUTFMT", Text, Some("0")),
    ("max_target_seqs", "max_target_seqs", "MAX_TARGET_SEQS", Integer, None),
    ("num_descriptions", "num_descriptions", "int_value", Count, None),
    ("num_alignments", "num_alignments", "int_value", Count, None),
    ("line_length", "line_length", "line_length", Count, None),
    ("sort_hits", "sorthits", "sort_hits", Small, None),
    ("sort_hsps", "sorthsps", "sort_hsps", Small, None),
    ("num_threads", "num_threads", "NUM_THREADS", Integer, Some("1")),
    ("mt_mode", "mt_mode", "int_value", Small, None),
    ("comp_based_stats", "comp_based_stats", "COMP_BASED_STATS", Text, Some("2")),
];

/// TBLASTN's arguments in v1's command, as `BLASTX`.
#[rustfmt::skip]
const TBLASTN: &[Frozen] = &[
    ("query", "query", "PATH", Path, Some("-")),
    ("subject", "subject", "PATH", Path, None),
    ("db", "db", "NAME", Text, None),
    ("seqidlist", "seqidlist", "filename", Text, None),
    ("negative_seqidlist", "negative_seqidlist", "filename", Text, None),
    ("taxids", "taxids", "taxids", Text, None),
    ("taxidlist", "taxidlist", "filename", Text, None),
    ("negative_taxids", "negative_taxids", "taxids", Text, None),
    ("negative_taxidlist", "negative_taxidlist", "filename", Text, None),
    ("no_taxid_expansion", "no_taxid_expansion", "NO_TAXID_EXPANSION", Flag, None),
    ("db_soft_mask", "db_soft_mask", "filtering_algorithm", Text, None),
    ("db_hard_mask", "db_hard_mask", "filtering_algorithm", Text, None),
    ("dbsize", "dbsize", "LETTERS", Int8, None),
    ("show_gis", "show_gis", "SHOW_GIS", Flag, Some("false")),
    ("html", "html", "HTML", Flag, Some("false")),
    ("parse_deflines", "parse_deflines", "PARSE_DEFLINES", Flag, Some("false")),
    ("searchsp", "searchsp", "INT_VALUE", Int8, None),
    ("task", "task", "TASK", Text, Some("tblastn")),
    ("out", "out", "PATH", Path, None),
    ("export_search_strategy", "export_search_strategy", "filename", Path, None),
    ("evalue", "evalue", "EVALUE", Double, None),
    ("word_size", "word_size", "WORD_SIZE", Integer, None),
    ("gap_open", "gapopen", "GAP_OPEN", Integer, None),
    ("gap_extend", "gapextend", "GAP_EXTEND", Integer, None),
    ("db_gencode", "db_gencode", "DB_GENCODE", Small, Some("1")),
    ("max_intron_length", "max_intron_length", "MAX_INTRON_LENGTH", Integer, None),
    ("matrix", "matrix", "MATRIX", Text, None),
    ("threshold", "threshold", "THRESHOLD", Double, None),
    ("xdrop_ungap", "xdrop_ungap", "XDROP_UNGAP", Double, None),
    ("xdrop_gap", "xdrop_gap", "XDROP_GAP", Double, None),
    ("xdrop_gap_final", "xdrop_gap_final", "XDROP_GAP_FINAL", Double, None),
    ("comp_based_stats", "comp_based_stats", "COMP_BASED_STATS", Text, Some("2")),
    ("outfmt", "outfmt", "SPEC", Text, Some("0")),
    ("query_loc", "query_loc", "RANGE", Text, None),
    ("subject_loc", "subject_loc", "RANGE", Text, None),
    ("seg", "seg", "SEG", Text, None),
    ("lcase_masking", "lcase_masking", "LCASE_MASKING", Flag, None),
    ("use_sw_tback", "use_sw_tback", "USE_SW_TBACK", Flag, None),
    ("soft_masking", "soft_masking", "SOFT_MASKING", Boolean, None),
    ("sum_stats", "sum_stats", "SUM_STATS", Boolean, None),
    ("window_size", "window_size", "WINDOW_SIZE", Integer, None),
    ("max_target_seqs", "max_target_seqs", "MAX_TARGET_SEQS", Count, None),
    ("num_descriptions", "num_descriptions", "int_value", Count, None),
    ("num_alignments", "num_alignments", "int_value", Count, None),
    ("line_length", "line_length", "line_length", Count, None),
    ("sort_hits", "sorthits", "sort_hits", Small, None),
    ("sort_hsps", "sorthsps", "sort_hsps", Small, None),
    ("max_hsps", "max_hsps", "MAX_HSPS", Count, None),
    ("qcov_hsp_perc", "qcov_hsp_perc", "QCOV_HSP_PERC", Double, None),
    ("culling_limit", "culling_limit", "CULLING_LIMIT", Integer, None),
    ("best_hit_overhang", "best_hit_overhang", "float_value", Double, None),
    ("best_hit_score_edge", "best_hit_score_edge", "float_value", Double, None),
    ("subject_besthit", "subject_besthit", "SUBJECT_BESTHIT", Flag, None),
    ("num_threads", "num_threads", "NUM_THREADS", Count, Some("1")),
    ("mt_mode", "mt_mode", "int_value", Small, None),
    ("ungapped", "ungapped", "UNGAPPED", Flag, None),
    ("in_pssm", "in_pssm", "psi_chkpt_file", Path, None),
];

/// The argument `frozen` of v1's command, as clap's derive built it (without its help text, which
/// v1 never prints: its callers refuse the help words).
fn frozen_arg(&(id, long, value_name, reads, default): &Frozen) -> clap::Arg {
    use crate::blastinput::value_parsers::{ncbi_boolean, ncbi_integer_as, ncbi_path};
    let arg = clap::Arg::new(id).long(long).value_name(value_name);
    let arg = match reads {
        Flag => arg.action(clap::ArgAction::SetTrue),
        Text => arg.value_parser(ValueParser::string()),
        Path => arg.value_parser(ncbi_path()),
        Integer => arg.value_parser(ncbi_integer),
        Int8 => arg.value_parser(ncbi_int8),
        Double => arg.value_parser(ncbi_double),
        Boolean => arg.value_parser(ncbi_boolean),
        Count => arg.value_parser(ncbi_integer_as::<usize>),
        Small => arg.value_parser(ncbi_integer_as::<u8>),
    };
    let arg = match reads {
        Flag => arg,
        _ => arg.action(clap::ArgAction::Set),
    };
    match default {
        Some(default) => arg.default_value(default),
        None => arg,
    }
}

/// v1's command of `program` (BLASTX or TBLASTN) as LOSAT's command line had it when v1 was frozen
/// (`Cli` of that program: `losat`, its version, and the program's arguments).
fn frozen_command(program: &str) -> clap::Command {
    let (name, args) = match program {
        "blastx" => ("blastx", BLASTX),
        "tblastn" => ("tblastn", TBLASTN),
        _ => unreachable!("v1 reads BLASTX and TBLASTN"),
    };
    let scope = clap::Command::new(name).args(args.iter().map(frozen_arg));
    clap::Command::new("losat")
        .version(env!("CARGO_PKG_VERSION"))
        .subcommand_required(true)
        .subcommand(scope)
}

/// The long names of the programs' arguments of LOSAT's command line when v1 was frozen, in the
/// order of its command (each program's, from blastx to rpsbproc), which v1's messages write with
/// one dash (`render`).
#[rustfmt::skip]
const LONGS: &[&str] = &[
    "query", "subject", "subject_loc", "db", "seqidlist", "negative_seqidlist", "taxids",
    "taxidlist", "negative_taxids", "negative_taxidlist", "no_taxid_expansion", "db_soft_mask",
    "db_hard_mask", "dbsize", "show_gis", "html", "parse_deflines", "searchsp", "out",
    "export_search_strategy", "task", "evalue", "word_size", "max_hsps", "qcov_hsp_perc",
    "sum_stats", "gapopen", "gapextend", "max_intron_length", "seg", "soft_masking", "matrix",
    "threshold", "window_size", "xdrop_ungap", "xdrop_gap", "xdrop_gap_final", "culling_limit",
    "best_hit_overhang", "best_hit_score_edge", "subject_besthit", "ungapped", "lcase_masking",
    "use_sw_tback", "strand", "query_loc", "query_gencode", "outfmt", "max_target_seqs",
    "num_descriptions", "num_alignments", "line_length", "sorthits", "sorthsps", "num_threads",
    "mt_mode", "comp_based_stats", "query", "subject", "db", "seqidlist", "negative_seqidlist",
    "taxids", "taxidlist", "negative_taxids", "negative_taxidlist", "no_taxid_expansion",
    "db_soft_mask", "db_hard_mask", "dbsize", "show_gis", "html", "parse_deflines", "searchsp",
    "task", "template_type", "template_length", "word_size", "num_threads", "mt_mode", "use_index",
    "index_name", "evalue", "perc_identity", "max_target_seqs", "num_descriptions",
    "num_alignments", "line_length", "sorthits", "sorthsps", "max_hsps", "qcov_hsp_perc", "out",
    "export_search_strategy", "reward", "penalty", "gapopen", "gapextend", "dust", "lcase_masking",
    "soft_masking", "filtering_db", "window_masker_taxid", "window_masker_db", "strand",
    "query_loc", "subject_loc", "subject_besthit", "culling_limit", "best_hit_overhang",
    "best_hit_score_edge", "ungapped", "xdrop_ungap", "xdrop_gap", "xdrop_gap_final", "no_greedy",
    "min_raw_gapped_score", "window_size", "off_diagonal_range", "outfmt", "query", "subject",
    "db", "seqidlist", "negative_seqidlist", "taxids", "taxidlist", "negative_taxids",
    "negative_taxidlist", "no_taxid_expansion", "db_soft_mask", "db_hard_mask", "dbsize",
    "show_gis", "html", "parse_deflines", "searchsp", "task", "evalue", "threshold", "word_size",
    "num_threads", "mt_mode", "out", "export_search_strategy", "max_target_seqs",
    "num_descriptions", "num_alignments", "line_length", "sorthits", "sorthsps", "max_hsps",
    "qcov_hsp_perc", "culling_limit", "best_hit_overhang", "best_hit_score_edge",
    "subject_besthit", "ungapped", "xdrop_ungap", "xdrop_gap", "xdrop_gap_final", "window_size",
    "matrix", "gapopen", "gapextend", "comp_based_stats", "seg", "soft_masking", "lcase_masking",
    "use_sw_tback", "query_loc", "subject_loc", "outfmt", "query", "subject", "db", "seqidlist",
    "negative_seqidlist", "taxids", "taxidlist", "negative_taxids", "negative_taxidlist",
    "no_taxid_expansion", "db_soft_mask", "db_hard_mask", "dbsize", "show_gis", "html",
    "parse_deflines", "searchsp", "evalue", "threshold", "matrix", "word_size", "num_threads",
    "out", "export_search_strategy", "query_gencode", "db_gencode", "max_intron_length",
    "max_target_seqs", "num_descriptions", "num_alignments", "line_length", "sorthits", "sorthsps",
    "max_hsps", "qcov_hsp_perc", "seg", "soft_masking", "sum_stats", "lcase_masking", "strand",
    "xdrop_ungap", "window_size", "outfmt", "culling_limit", "best_hit_overhang",
    "best_hit_score_edge", "subject_besthit", "query_loc", "subject_loc", "query", "subject", "db",
    "seqidlist", "negative_seqidlist", "taxids", "taxidlist", "negative_taxids",
    "negative_taxidlist", "no_taxid_expansion", "db_soft_mask", "db_hard_mask", "dbsize",
    "show_gis", "html", "parse_deflines", "searchsp", "task", "out", "export_search_strategy",
    "evalue", "word_size", "gapopen", "gapextend", "db_gencode", "max_intron_length", "matrix",
    "threshold", "xdrop_ungap", "xdrop_gap", "xdrop_gap_final", "comp_based_stats", "outfmt",
    "query_loc", "subject_loc", "seg", "lcase_masking", "use_sw_tback", "soft_masking",
    "sum_stats", "window_size", "max_target_seqs", "num_descriptions", "num_alignments",
    "line_length", "sorthits", "sorthsps", "max_hsps", "qcov_hsp_perc", "culling_limit",
    "best_hit_overhang", "best_hit_score_edge", "subject_besthit", "num_threads", "mt_mode",
    "ungapped", "in_pssm", "query", "subject", "db", "seqidlist", "negative_seqidlist", "taxids",
    "taxidlist", "negative_taxids", "negative_taxidlist", "no_taxid_expansion", "dbsize",
    "show_gis", "parse_deflines", "searchsp", "sum_stats", "evalue", "threshold", "word_size",
    "num_threads", "out", "export_search_strategy", "max_target_seqs", "num_descriptions",
    "num_alignments", "line_length", "sorthits", "sorthsps", "max_hsps", "qcov_hsp_perc",
    "culling_limit", "best_hit_overhang", "best_hit_score_edge", "subject_besthit", "xdrop_ungap",
    "xdrop_gap", "xdrop_gap_final", "gap_trigger", "window_size", "matrix", "gapopen", "gapextend",
    "comp_based_stats", "seg", "soft_masking", "lcase_masking", "use_sw_tback", "query_loc",
    "subject_loc", "outfmt", "num_iterations", "pseudocount", "inclusion_ethresh", "out_pssm",
    "out_ascii_pssm", "save_pssm_after_last_round", "save_each_pssm", "in_pssm", "in_msa",
    "msa_master_idx", "ignore_msa_master", "query", "db", "out", "export_search_strategy",
    "evalue", "dbsize", "searchsp", "comp_based_stats", "outfmt", "show_gis", "num_descriptions",
    "num_alignments", "line_length", "sorthits", "sorthsps", "max_target_seqs", "max_hsps",
    "qcov_hsp_perc", "xdrop_ungap", "xdrop_gap", "xdrop_gap_final", "sum_stats", "use_sw_tback",
    "seg", "soft_masking", "window_size", "lcase_masking", "query_loc", "parse_deflines",
    "num_threads", "mt_mode", "query", "query_gencode", "query_loc", "strand", "db", "out",
    "export_search_strategy", "evalue", "comp_based_stats", "outfmt", "show_gis",
    "num_descriptions", "num_alignments", "line_length", "sorthits", "sorthsps", "seg",
    "soft_masking", "lcase_masking", "qcov_hsp_perc", "max_hsps", "max_target_seqs", "dbsize",
    "searchsp", "sum_stats", "xdrop_ungap", "xdrop_gap", "xdrop_gap_final", "window_size",
    "ungapped", "parse_deflines", "num_threads", "mt_mode", "use_sw_tback", "in", "input_type",
    "dbtype", "title", "out", "blastdb_version", "max_file_sz", "metadata_output_prefix",
    "logfile", "hash_index", "taxid", "taxid_map", "mask_data", "mask_id", "mask_desc", "gi_mask",
    "gi_mask_name", "parse_seqids", "oid_masks", "randomize", "i", "o", "e", "m", "t", "d", "f",
    "q",
];

/// The exclusions of v1's BLASTX and TBLASTN (the clap attributes of their arguments when v1 was
/// frozen): an argument and the arguments it excludes.
const V1_EXCLUSIONS: &[(&str, &[&str])] = &[
    (
        "subject_loc",
        &crate::blastinput::db_filter::DB_FILTER_OPTIONS,
    ),
    ("db", &["subject", "subject_loc"]),
    ("best_hit_overhang", &["culling_limit"]),
    ("best_hit_score_edge", &["culling_limit"]),
    ("max_target_seqs", &["num_descriptions", "num_alignments"]),
    (
        "seqidlist",
        &[
            "negative_seqidlist",
            "taxids",
            "taxidlist",
            "negative_taxids",
            "negative_taxidlist",
            "no_taxid_expansion",
            "subject",
        ],
    ),
    (
        "negative_seqidlist",
        &[
            "taxids",
            "taxidlist",
            "negative_taxids",
            "negative_taxidlist",
            "no_taxid_expansion",
            "subject",
        ],
    ),
    (
        "taxids",
        &[
            "taxidlist",
            "negative_taxids",
            "negative_taxidlist",
            "subject",
        ],
    ),
    (
        "taxidlist",
        &["negative_taxids", "negative_taxidlist", "subject"],
    ),
    ("negative_taxids", &["negative_taxidlist", "subject"]),
    ("negative_taxidlist", &["subject"]),
    ("no_taxid_expansion", &["subject"]),
    ("db_soft_mask", &["subject", "subject_loc", "db_hard_mask"]),
    ("db_hard_mask", &["subject", "subject_loc"]),
];

/// The value parsers of v1's `program` (BLASTX or TBLASTN) where they are not the command line's
/// conversion: NCBI's constraints, checked with v1's own messages.
fn v1_parsers(program: &str) -> Vec<(&'static str, ValueParser)> {
    let mut parsers = vec![
        ("searchsp", ValueParser::new(nonnegative_ncbi_int8)),
        ("out", ValueParser::new(output_path())),
        ("export_search_strategy", ValueParser::new(output_path())),
        ("word_size", ValueParser::new(protein_word_size)),
        ("qcov_hsp_perc", ValueParser::new(percentage)),
        (
            "max_intron_length",
            ValueParser::new(nonnegative_ncbi_integer),
        ),
        ("window_size", ValueParser::new(nonnegative_ncbi_integer)),
        ("culling_limit", ValueParser::new(nonnegative_ncbi_integer)),
        ("threshold", ValueParser::new(ncbi_threshold)),
        ("best_hit_overhang", ValueParser::new(best_hit_value)),
        ("best_hit_score_edge", ValueParser::new(best_hit_value)),
        ("num_descriptions", ValueParser::new(report_count)),
        ("num_alignments", ValueParser::new(report_count)),
        ("mt_mode", ValueParser::new(ncbi_mt_mode)),
        ("line_length", ValueParser::new(formatting_line_length)),
        ("sort_hits", ValueParser::new(formatting_sort_option)),
        ("sort_hsps", ValueParser::new(formatting_sort_option)),
    ];
    if program == "blastx" {
        parsers.extend([
            (
                "task",
                ValueParser::new(PossibleValuesParser::new(["blastx", "blastx-fast"])),
            ),
            ("max_hsps", ValueParser::new(ncbi_positive_integer)),
            (
                "strand",
                ValueParser::new(PossibleValuesParser::new(["both", "plus", "minus"])),
            ),
            ("query_gencode", ValueParser::new(genetic_code)),
            ("max_target_seqs", ValueParser::new(ncbi_positive_integer)),
            ("num_threads", ValueParser::new(ncbi_positive_integer)),
        ]);
    } else {
        parsers.extend([
            (
                "task",
                ValueParser::new(PossibleValuesParser::new(["tblastn", "tblastn-fast"])),
            ),
            ("db_gencode", ValueParser::new(tblastn_genetic_code)),
            ("max_target_seqs", ValueParser::new(blastn_count)),
            ("max_hsps", ValueParser::new(blastn_count)),
            ("num_threads", ValueParser::new(blastn_count)),
        ]);
    }
    parsers
}

/// The arguments of the request `argv` (`losat <program> <words>`, BLASTX or TBLASTN) read with
/// `command` (the program's command, `cli::command_of`, with v1's own value parsers), or the
/// message of their error as v1 gives it (`render`). ABI v1 runs no NCBI application start: the
/// versions and the toolkit's words that the start reads (`-logfile`, `-conffile`, `-dryrun`) are
/// refused; its callers refuse the help words before.
pub(crate) fn parse(command: clap::Command, argv: Vec<OsString>) -> Result<Commands, String> {
    read(command, argv).map_err(|error| render(&error))
}

/// `parse` without the rendering of the error: v1's command reads the words, then the program's
/// options are read from them (`checked_options`).
fn read(command: clap::Command, argv: Vec<OsString>) -> Result<Commands, clap::Error> {
    let (program, given) = checked_words(&command, argv)?;
    checked_options(&program, &given).map_err(raw)
}

/// The program of `argv` and the words that v1's `command` took (each long name and its value,
/// `None` for a flag), or v1's error.
fn checked_words(
    command: &clap::Command,
    argv: Vec<OsString>,
) -> Result<(String, Vec<(String, Option<OsString>)>), clap::Error> {
    use crate::blastinput::ncbi_environment::{preparse, Preparsed};
    let program = argv
        .get(1)
        .and_then(|word| word.to_str())
        .unwrap_or_default()
        .to_string();
    let words = argv.get(2..).unwrap_or(&[]);
    // The refusals name the program as its other refusals do (`Tool::refusal_name`).
    let tool =
        crate::cli::reader::tool(&program).expect("v1 reads the words of one of its programs");
    let words = match preparse(words) {
        Preparsed::Version(request) => {
            return Err(raw(crate::cli::version_word_refusal(tool, request.word())))
        }
        Preparsed::Words { words: read, .. } => {
            if let Some(word) = crate::cli::toolkit_file_word(words, &program) {
                return Err(raw(crate::cli::toolkit_file_word_refusal(tool, word)));
            }
            read
        }
    };
    let scope = command
        .find_subcommand(&program)
        .expect("v1 reads the words of one of its programs");
    let mut translated: Vec<OsString> = vec!["losat".into(), program.clone().into()];
    let mut given: Vec<(String, Option<OsString>)> = Vec::new();
    let mut input = words.into_iter();
    while let Some(token) = input.next() {
        let text = token.to_str().unwrap_or("").to_string();
        // A last `--` changes nothing in TBLASTN, as NCBI's reader; BLASTX's reader of v1 did not
        // take it.
        if text == "--" && input.len() == 0 && program == "tblastn" {
            continue;
        }
        let key = text.strip_prefix('-').filter(|key| !key.starts_with('-'));
        let (name, inline) = key
            .unwrap_or("")
            .split_once('=')
            .map_or((key.unwrap_or(""), None), |(k, v)| (k, Some(v)));
        let Some(arg) = scope
            .get_arguments()
            .find(|arg| arg.get_long() == Some(name))
        else {
            return Err(unknown_option(&program, name, &text));
        };
        if arg.get_action().takes_values() {
            let value = match inline {
                Some(value) => OsString::from(value),
                None => input.next().ok_or_else(|| {
                    clap::Error::raw(
                        ErrorKind::InvalidValue,
                        format!("-{name} requires one value"),
                    )
                })?,
            };
            let mut word = OsString::from(format!("--{name}="));
            word.push(&value);
            translated.push(word);
            given.push((name.to_string(), Some(value)));
        } else {
            if inline.is_some() {
                return Err(clap::Error::raw(
                    ErrorKind::TooManyValues,
                    format!("-{name} is a flag and does not take a value"),
                ));
            }
            translated.push(format!("--{name}").into());
            given.push((name.to_string(), None));
        }
    }
    command.clone().try_get_matches_from(&translated)?;
    Ok((program, given))
}

/// The options of v1's `program` of the words `given` that its command took (`checked_words`):
/// NCBI's arguments of them (`ncbi_app::args::checked`), read into the program's options as the
/// command line reads its own (`cli::reader::options`).
fn checked_options(
    program: &str,
    given: &[(String, Option<OsString>)],
) -> Result<Commands, String> {
    let tool = crate::cli::reader::tool(program).expect("v1's programs are on the reader");
    let args = crate::ncbi_app::args::checked(tool, given).ok_or_else(|| {
        format!(
            "LOSAT's {} does not convert the arguments that its reader read",
            tool.refusal_name
        )
    })?;
    crate::cli::reader::options(tool, &args, false)
}

/// A refusal of v1's reader, as its reader wrote them (a clap error of LOSAT's text).
fn raw(message: String) -> clap::Error {
    clap::Error::raw(ErrorKind::InvalidValue, message)
}

/// v1's rejection of an option that the program does not have: NCBI's options that LOSAT
/// refuses by name, and the rest.
fn unknown_option(program: &str, name: &str, text: &str) -> clap::Error {
    let tool = crate::cli::reader::tool(program).expect("v1's programs are on the reader");
    match crate::cli::reader::named_refusal(tool, name) {
        Some(refusal) => raw(refusal),
        None => clap::Error::raw(
            ErrorKind::UnknownArgument,
            format!("unknown option or argument '{text}'; use -help to list the options"),
        ),
    }
}

/// The text of v1's `error`: clap's, with each program's long options written with one dash, as
/// v1's messages have them.
fn render(error: &clap::Error) -> String {
    let mut message = error
        .to_string()
        .replace("-h, --help", "-help")
        .replace("-V, --version", "--version");
    let toolkit = ["logfile", "conffile"].into_iter();
    let longs = LONGS.iter().copied().chain(toolkit);
    for name in longs.filter(|name| !matches!(*name, "help" | "version") && name.len() > 1) {
        message = message.replace(&format!("--{name}"), &format!("-{name}"));
    }
    message
        .replace("'-h'", "'-help'")
        .replace("-h, --help", "-help")
        .replace("-V, --version", "--version")
        .replace("see a summary with '-help'", "see a summary with '-h'")
        .replace("see more with '--help'", "see more with '-help'")
}

// v1's copies of the value parsers of the command line before line RG item K4i-2b (D5).

/// An integer argument with a constraint: after the conversion, NCBI's constraint reads the
/// value again with `NStr::StringToDouble` (default flags, `ncbi_str::string_to_double`), which
/// does not read a `0x` prefix without digits (the only form that `ncbi_integer` reads and
/// `strtod` does not end at), so the value is illegal.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:110-113
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1231-1243
fn ncbi_constrained_integer(value: &str) -> Result<i32, String> {
    let n = ncbi_integer(value)?;
    if !reads_as_double(value) {
        return Err("Illegal value".into());
    }
    Ok(n)
}

/// Whether the constraint's `NStr::StringToDouble` of `value` (default flags) reads `value`
/// (`ncbi_constrained_integer`).
fn reads_as_double(value: &str) -> bool {
    crate::ncbi_str::string_to_double(value.as_bytes()).is_ok()
}

/// An integer argument of `CArg_Integer` with NCBI's lower bound (`at_least`,
/// `CArgAllowValuesGreaterThanOrEqual`).
fn ncbi_integer_at_least(value: &str, at_least: i32) -> Result<i32, String> {
    let n = ncbi_constrained_integer(value)?;
    if n < at_least {
        return Err(format!("expected an integer >= {at_least}"));
    }
    Ok(n)
}

/// A non-negative `Int8` argument (`CArg_Int8` with `CArgAllowValuesGreaterThanOrEqual` of 0):
/// `-searchsp`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:234-238
/// The constraint reads the value as `ncbi_constrained_integer`'s does.
pub(crate) fn nonnegative_ncbi_int8(value: &str) -> Result<i64, String> {
    let n = ncbi_int8(value)?;
    if !reads_as_double(value) {
        return Err("Illegal value".into());
    }
    if n < 0 {
        return Err("expected an integer >= 0".into());
    }
    Ok(n)
}

/// A non-negative integer argument (`CArg_Integer` with `CArgAllowValuesGreaterThanOrEqual` of 0),
/// such as `-culling_limit`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3296-3302
fn nonnegative_ncbi_integer(value: &str) -> Result<i32, String> {
    ncbi_integer_at_least(value, 0)
}

/// A `-num_descriptions` or `-num_alignments` value: 0 or more, an `Int4` in NCBI's
/// integer forms (the sign and hexadecimal forms of `nonnegative_ncbi_integer`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2665-2685
fn report_count(value: &str) -> Result<usize, String> {
    nonnegative_ncbi_integer(value).map(|n| n as usize)
}

/// A BLASTN count of 1 or more (`-num_threads`, `-max_target_seqs`, `-max_hsps`), as
/// NCBI's arguments (c++/src/algo/blast/blastinput/blast_args.cpp:203-207,2731-2732,3162-3163).
fn blastn_count(value: &str) -> Result<usize, String> {
    ncbi_positive_integer(value).map(|n| n as usize)
}

/// `blastn_count` as an `int`, the type of BLASTX's fields.
fn ncbi_positive_integer(value: &str) -> Result<i32, String> {
    ncbi_integer_at_least(value, 1)
}

/// A percentage, in 0..=100 (NaN is outside): `-perc_identity`, and
/// `-qcov_hsp_perc` of every program (the `*_query_coverage` parsers).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:196-201
/// The constraint reads the value with the argument's own conversion (`ncbi_double`), so it
/// tests the argument's number (`1e` is 1, `0x10` is 16).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4821-4830
fn percentage(value: &str) -> Result<f64, String> {
    let n = ncbi_double(value)?;
    if !(0.0..=100.0).contains(&n) {
        return Err("expected a percentage in 0..=100".into());
    }
    Ok(n)
}

/// A `-best_hit_overhang` or `-best_hit_score_edge` value: above 0 and below 0.5 (NaN is
/// outside), of every program that has the options.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3304-3328
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/hspfilter_besthit.h:80-91
/// The constraint reads the value with the argument's own conversion (`ncbi_double`); NCBI
/// BLAST+ 2.17.0 refuses 0 and 0.5 with "Illegal value, expected (>0 and <0.5)" (line NP O6a's
/// probe).
fn best_hit_value(value: &str) -> Result<f64, String> {
    let n = ncbi_double(value)?;
    if !(n > 0.0 && n < 0.5) {
        return Err("expected a number above 0 and below 0.5".into());
    }
    Ok(n)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:995-1006
fn genetic_code(value: &str) -> Result<u8, String> {
    let n = genetic_code_number(value).ok_or("invalid genetic code")?;
    if !crate::ncbi_app::args::is_genetic_code(i32::from(n)) {
        return Err("unsupported genetic code".into());
    }
    Ok(n)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-51
/// A BLASTN output file: a name shorter than 256 bytes, as NCBI's constraint; `-` is
/// standard output, and a file that cannot be created is reported when it is opened
/// (`blastn/blast_engine/run.rs`).
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:79-87
fn output_path() -> impl clap::builder::TypedValueParser<Value = PathBuf> {
    use clap::builder::TypedValueParser;
    clap::builder::OsStringValueParser::new().try_map(|value| {
        if crate::ncbi_app::args::file_name_length(value.as_encoded_bytes()) >= 256 {
            return Err("Illegal value, expected file name length < 256".to_string());
        }
        Ok(PathBuf::from(value))
    })
}

/// A `-threshold` value: NCBI's constraint reads it again with `NStr::StringToDouble` and its
/// default flags (`ncbi_str::string_to_double`: glibc's `strtod` alone, so a value that only
/// `StringToDoublePosix` reads, such as `1e`, makes it throw), and requires at least 0 (a NaN
/// fails it); a value that the constraint reads is the same number as the argument's.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:578-583
/// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:110-113
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1231-1243
fn ncbi_threshold(value: &str) -> Result<f64, String> {
    let threshold = ncbi_double(value)?;
    // The check that `value` is at least `m_MinValue` (a NaN fails it).
    let allowed =
        crate::ncbi_str::string_to_double(value.as_bytes()).is_ok_and(|checked| checked >= 0.0);
    if !allowed {
        return Err("Illegal value, expected greater or equal to 0".into());
    }
    Ok(threshold)
}

/// A `-mt_mode` value: an integer from 0 to 2 (`CArgAllowValuesBetween` from 0 to 2, with `true`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3176-3177
fn ncbi_mt_mode(value: &str) -> Result<u8, String> {
    let mode = ncbi_integer_at_least(value, 0)?;
    u8::try_from(mode)
        .ok()
        .filter(|mode| *mode <= 2)
        .ok_or_else(|| "Illegal value: expected an integer from 0 to 2".to_string())
}

/// A protein `-word_size` (`CArgAllowValuesGreaterThanOrEqual` of 2 for a protein query).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:166-170
fn protein_word_size(value: &str) -> Result<i32, String> {
    ncbi_integer_at_least(value, 2)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:997-1056
// Stage A product decision admits 32 from gc.prt:340-347 for TBLASTN only.
/// `NStr::StringToInt` of `value` in `CArgAllowGeneticCodeInteger` is
/// `value_parsers::genetic_code_number` (LOSAT's "invalid" message for the values outside 0..=255).
fn tblastn_genetic_code(value: &str) -> std::result::Result<u8, String> {
    let id = genetic_code_number(value).ok_or("invalid genetic code ID")?;
    GeneticCode::try_from_id(id).map(|_| id)
}

/// A `-sorthits` or `-sorthsps` value: an integer from 0 to 4 (`CArgAllowValuesBetween` from 0 to 4,
/// with `true`, cited at `FormattingArgs`).
fn formatting_sort_option(value: &str) -> Result<u8, String> {
    let option = ncbi_integer_at_least(value, 0)?;
    u8::try_from(option)
        .ok()
        .filter(|option| *option <= 4)
        .ok_or_else(|| "Illegal value: expected an integer from 0 to 4".to_string())
}

/// A `-line_length` value: an integer of 1 or more (`CArgAllowValuesGreaterThanOrEqual` of 1,
/// cited at `FormattingArgs`).
fn formatting_line_length(value: &str) -> Result<usize, String> {
    let length = ncbi_integer_at_least(value, 1)?;
    usize::try_from(length).map_err(|_| "Illegal value: expected an integer >= 1".to_string())
}

#[cfg(test)]
#[path = "web_v1_argv_tests.rs"]
mod tests;
