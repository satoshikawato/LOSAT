//! TBLASTN's options (NCBI c++/src/algo/blast/blastinput/tblastn_args.cpp and blast_args.cpp)
//! and LOSAT's limits.
use anyhow::Result;
use clap::Args;
use std::io::Write;
use std::path::PathBuf;

use crate::blastinput::value_parsers::*;

use crate::utils::genetic_code::GeneticCode;

use super::stage_d_pipeline::LocalStageDScoring;
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::redo_alignment::BlastCompoAdjustMode;

// NCBI reference: c++/src/algo/blast/blastinput/tblastn_args.cpp:54-61
// ```c++
//     static const char kDefaultTask[] = "tblastn";
//     SetTask(kDefaultTask);
//     set<string> tasks;
//     tasks.insert(kDefaultTask);
//     tasks.insert("tblastn-fast");
//     arg.Reset(new CTaskCmdLineArgs(tasks, kDefaultTask));
//     m_Args.push_back(arg);
// ```
/// TBLASTN's `-task`: the command line's parser, and the check of `run_local` (`TaskParser::check`).
pub(crate) const TBLASTN_TASK: TaskParser =
    task_parser("TBLASTN", &["tblastn", "tblastn-fast"], &[]);

// NCBI reference: c++/src/algo/blast/api/tblastn_options.cpp:54-85
// ```c++
// SetWordThreshold(BLAST_WORD_THRESHOLD_TBLASTN);
// m_Opts->SetSumStatisticsMode();
// m_Opts->SetCompositionBasedStats(eCompositionMatrixAdjust);
// SetDbGeneticCode(BLAST_GENETIC_CODE);
// ```
// The arguments are read with NCBI's grammar (`blastinput/value_parsers.rs`); the values
// that NCBI checks after parsing are checked by `check_options`, in NCBI's order. The
// NCBI options that LOSAT's TBLASTN does not implement are rejected by the parser
// (`cli.rs`).
#[derive(Args, Debug, Clone)]
#[command(rename_all = "snake_case")]
pub struct TblastnArgs {
    // NCBI blast_args.cpp:3425-3427:
    // arg_desc.AddDefaultKey(kArgQuery, "input_file", "Input file name",
    //                        CArgDescriptions::eInputFile, kDfltArgQuery);
    // The default `-` is standard input.
    #[arg(long, value_parser = ncbi_input_path(), value_name = "PATH", default_value = "-")]
    pub query: PathBuf,
    #[arg(long, value_parser = ncbi_input_path(), value_name = "PATH")]
    pub subject: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2364-2371
    // ```c++
    //         arg_desc.AddOptionalKey(kArgSubject, "subject_input_file",
    //                                 "Subject sequence(s) to search",
    //                                 CArgDescriptions::eInputFile);
    //         ITERATE(vector<string>, dbarg, database_args) {
    //             arg_desc.SetDependency(kArgSubject, CArgDescriptions::eExcludes,
    //                                    *dbarg);
    //         }
    // ```
    // (and the same exclusion for `-subject_loc`, blast_args.cpp:2377-2381). The exclusion
    // is an error of the argument parser: LOSAT's parser reports it with its own text
    // (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES, exception 1). The database is searched by
    // `tblastn/database.rs`.
    #[arg(
        long,
        value_name = "NAME",
        conflicts_with_all = ["subject", "subject_loc"],
        help = "BLAST database name (a nucleotide database, with or without parsed Seq-ids)"
    )]
    pub db: Option<String>,
    /// `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist` and
    /// `-no_taxid_expansion` (line DI-8).
    #[command(flatten)]
    pub db_filter: crate::blastinput::db_filter::DbFilterArgs,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2231-2237
    // ```c++
    //     // DB size
    //     if (!m_IsMapper) {
    //         arg_desc.SetCurrentGroup("Statistical options");
    //         arg_desc.AddOptionalKey(kArgDbSize, "num_letters",
    //                                 "Effective length of the database ",
    //                                 CArgDescriptions::eInt8);
    //     }
    // ```
    // No constraint and no dependency: read with `-db` and with `-subject`
    // (`stage_d_pipeline::SubjectSetNumbers`).
    #[arg(long, value_parser = ncbi_int8, value_name = "LETTERS", help = "Effective length of the database")]
    pub dbsize: Option<i64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2663
    // ```c++
    //     arg_desc.AddFlag(kArgShowGIs, "Show NCBI GIs in deflines?", true);
    // ```
    #[arg(
        long = "show_gis",
        default_value_t = false,
        help = "Show NCBI GIs in deflines?"
    )]
    pub show_gis: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2050-2052
    // ```c++
    //     arg_desc.AddDefaultKey(kArgParseDeflines, "TF", "Should the query and "
    // ```
    // `-parse_deflines` is a flag of the programs (blast_args.cpp:1961).
    #[arg(
        long = "parse_deflines",
        default_value_t = false,
        help = "Should the query and subject defline(s) be parsed?"
    )]
    pub parse_deflines: bool,
    /// The ids of the search's records under `-parse_deflines` (`run`).
    #[arg(skip)]
    pub id_scope: Option<std::sync::Arc<crate::seqid::id_scope::IdScope>>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:234-238
    // ```c++
    //     arg_desc.AddOptionalKey(kArgEffSearchSpace, "int_value",
    //                             "Effective length of the search space",
    //                             CArgDescriptions::eInt8);
    // ```
    #[arg(long, value_parser = crate::blastinput::value_parsers::nonnegative_ncbi_int8, value_name = "INT_VALUE", help = "Effective length of the search space")]
    pub searchsp: Option<i64>,
    // `tblastn-fast`: word size 5, the compressed lookup table and threshold 20 (`resolve`).
    #[arg(long, default_value = "tblastn", value_parser = TBLASTN_TASK)]
    pub task: String,
    #[arg(long, value_name = "PATH", value_parser = ncbi_output_path())]
    pub out: Option<PathBuf>,
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3532-3535
    /// (quoted at `BlastpArgs::export_search_strategy`). The search strategy file
    /// (`report::archive::strategy`), created by `run` before the other options are read.
    #[arg(long, value_name = "filename", value_parser = ncbi_output_path(), help = "File name to record the search strategy used")]
    pub export_search_strategy: Option<PathBuf>,
    // An omitted value keeps the default of the task (10).
    #[arg(long, value_parser = tblastn_real, help = "Expectation value (E) threshold for saving hits [default: 10]")]
    pub evalue: Option<f64>,
    #[arg(long, value_parser = protein_word_size, help = "Word size for wordfinder algorithm [default: 3; 5 for tblastn-fast]")]
    pub word_size: Option<i32>,
    // NCBI blast_args.cpp:258-273: when omitted, each cost keeps the selected matrix's
    // BLAST_MATRIX_BEST value.
    #[arg(long = "gapopen", value_parser = ncbi_integer, help = "Cost to open a gap [default: 11, or the matrix's best value]")]
    pub gap_open: Option<i32>,
    #[arg(long = "gapextend", value_parser = ncbi_integer, help = "Cost to extend a gap [default: 1, or the matrix's best value]")]
    pub gap_extend: Option<i32>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1029-1056
    // ```c++
    // arg_desc.AddDefaultKey(kArgDbGeneticCode, "int_value",
    //                        "Genetic code to use to translate "
    //                        "database/subjects (see user manual for details)\n",
    //                        CArgDescriptions::eInteger,
    //                        NStr::IntToString(BLAST_GENETIC_CODE));
    // opt.SetDbGeneticCode(args[kArgDbGeneticCode].AsInteger());
    // ```
    // PD-TLOSAN-LOCAL-GENCODE-32 adds gc.prt ID 32 only for TBLASTN.
    #[arg(long, default_value_t = 1, value_parser = tblastn_genetic_code)]
    pub db_gencode: u8,
    #[arg(long, value_parser = nonnegative_ncbi_integer, help = "Length of the largest intron allowed in a translated nucleotide sequence when linking multiple distinct alignments [default: 0]")]
    pub max_intron_length: Option<i32>,
    #[arg(long, help = "Scoring matrix name [default: BLOSUM62]")]
    pub matrix: Option<String>,
    #[arg(long, value_parser = tblastn_threshold_value, help = "Minimum word score such that the word is added to the BLAST lookup table [default: 13, or the matrix's suggestion; 19.3, 21 or 20.25 for word sizes 5, 6 and 7; 20 for tblastn-fast]")]
    pub threshold: Option<f64>,
    // NCBI c++/src/algo/blast/blastinput/blast_args.cpp:220-229,280-286:
    // arg_desc.AddOptionalKey(kArgGappedXDropoff, "float_value", ...);
    // arg_desc.AddOptionalKey(kArgFinalGappedXDropoff, "float_value", ...);
    // opt.SetGapXDropoff(args[kArgGappedXDropoff].AsDouble());
    // opt.SetGapXDropoffFinal(args[kArgFinalGappedXDropoff].AsDouble());
    // NCBI c++/src/algo/blast/blastinput/blast_args.cpp:209-214,276-278:
    // arg_desc.AddOptionalKey(kArgUngappedXDropoff, "float_value", ...);
    // opt.SetXDropoff(args[kArgUngappedXDropoff].AsDouble());
    #[arg(long, value_parser = tblastn_real, help = "X-dropoff value (in bits) for ungapped extensions [default: 7]")]
    pub xdrop_ungap: Option<f64>,
    #[arg(long, value_parser = tblastn_real, help = "X-dropoff value (in bits) for preliminary gapped extensions [default: 15]")]
    pub xdrop_gap: Option<f64>,
    #[arg(long, value_parser = tblastn_real, help = "X-dropoff value (in bits) for final gapped alignment [default: 25]")]
    pub xdrop_gap_final: Option<f64>,
    // Read when NCBI's handler reads it (`blastinput/app.rs` `parse_comp_based_stats`).
    #[arg(
        long,
        default_value = "2",
        help = "Use composition-based statistics: 0 (F, f) none, 1 composition-based statistics, 2 (D, d, T, t) conditional compositional score matrix adjustment, 3 unconditional compositional score matrix adjustment"
    )]
    pub comp_based_stats: String,
    #[arg(long, default_value = "0")]
    pub outfmt: String,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1945-1949
    // ```c++
    //     // query location
    //     arg_desc.AddOptionalKey(kArgQueryLocation, "range",
    //                             "Location on the query sequence in 1-based offsets "
    //                             "(Format: start-stop)",
    //                             CArgDescriptions::eString);
    // ```
    // Read by the query options handler (`check_options`).
    #[arg(
        long = "query_loc",
        value_name = "RANGE",
        help = "Location on the query sequence in 1-based offsets (Format: start-stop)"
    )]
    pub query_loc: Option<String>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2372-2376
    // ```c++
    //         // subject location
    //         arg_desc.AddOptionalKey(kArgSubjectLocation, "range",
    //                         "Location on the subject sequence in 1-based offsets "
    //                         "(Format: start-stop)",
    //                         CArgDescriptions::eString);
    // ```
    // Read by the database arguments handler when it reads the subjects (`run`).
    #[arg(
        long = "subject_loc",
        value_name = "RANGE",
        conflicts_with_all = crate::blastinput::db_filter::DB_FILTER_OPTIONS,
        help = "Location on the subject sequence in 1-based offsets (Format: start-stop)"
    )]
    pub subject_loc: Option<String>,
    // Read when NCBI's filtering handler reads it (`blastinput/app.rs` `parse_seg_option`);
    // an omitted value is NCBI's default for tblastn, "12 2.2 2.5".
    #[arg(
        long,
        help = "SEG: no, yes, or \"WINDOW LOCUT HICUT\" [default: 12 2.2 2.5]"
    )]
    pub seg: Option<String>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1938-1942,2549-2557
    // arg_desc.AddFlag(kArgUseLCaseMasking,
    //     "Use lower case filtering in query and subject sequence(s)?", true);
    // ReadSequencesToBlast(..., use_lcase_masks, subjects, ...);
    #[arg(long)]
    pub lcase_masking: bool,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:339-342
    // ```c++
    //     arg_desc.AddDefaultKey(kArgLookupTableMaskingOnly, "soft_masking",
    //                            "Apply filtering locations as soft masks",
    //                            CArgDescriptions::eBoolean,
    //                            kDfltArgLookupTableMaskingOnlyProt);
    // ```
    #[arg(long, value_parser = ncbi_boolean, help = "Apply filtering locations as soft masks [default: false]")]
    pub soft_masking: Option<bool>,
    #[arg(long, value_parser = ncbi_boolean, help = "Use sum statistics [default: true]")]
    pub sum_stats: Option<bool>,
    #[arg(long, value_parser = nonnegative_ncbi_integer, help = "Multiple hits window size, use 0 to specify 1-hit algorithm [default: 40, or the matrix's suggestion]")]
    pub window_size: Option<i32>,
    // An omitted value keeps the default hit list size (500).
    #[arg(long, value_parser = blastn_count, help = "Maximum number of aligned sequences to keep [default: 500]", conflicts_with_all = ["num_descriptions", "num_alignments"])]
    pub max_target_seqs: Option<usize>,
    // NCBI blast_args.cpp:2665-2685,2732-2738: 0 or more; -max_target_seqs excludes both.
    #[arg(long = "num_descriptions", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show one-line descriptions for. Not applicable for outfmt > 4 [default: 500]")]
    pub num_descriptions: Option<usize>,
    #[arg(long = "num_alignments", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show alignments for [default: 250]")]
    pub num_alignments: Option<usize>,
    /// `-line_length`, `-sorthits` and `-sorthsps` (NCBI's formatting handler,
    /// `blastinput::report_counts::FormattingArgs`).
    #[command(flatten)]
    pub formatting: crate::blastinput::report_counts::FormattingArgs,
    // NCBI blast_args.cpp:196-206: -qcov_hsp_perc in 0..100, -max_hsps >= 1.
    #[arg(long = "max_hsps", value_parser = blastn_count, help = "Set maximum number of HSPs per subject sequence to save for each query")]
    pub max_hsps: Option<usize>,
    #[arg(long = "qcov_hsp_perc", value_parser = tblastn_query_coverage, help = "Percent query coverage per hsp")]
    pub qcov_hsp_perc: Option<f64>,
    // NCBI blast_args.cpp:3335-3340 (CHspFilteringArgs::ExtractAlgorithmOptions):
    // opts.SetCullingLimit(args[kArgCullingLimit].AsInteger());
    /// If the query range of a hit is enveloped by that of at least this many
    /// higher-scoring hits, delete the hit (default: 0, no culling).
    #[arg(long, value_parser = nonnegative_ncbi_integer)]
    pub culling_limit: Option<i32>,
    #[arg(long, default_value_t = 1, value_parser = blastn_count)]
    pub num_threads: usize,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3169-3181
    // ```c++
    //     if (m_MTMode >= 0) {
    //     	arg_desc.AddDefaultKey(kArgMTMode, "int_value",
    //                                "Multi-thread mode to use in BLAST search:\n "
    //                                "0 auto split by database or queries \n "
    //                                "1 split by queries\n "
    //                                "2 split by database",
    //                                CArgDescriptions::eInteger,
    //                                NStr::IntToString(0));
    //     	arg_desc.SetConstraint(kArgMTMode,
    //                                new CArgAllowValuesBetween(0, 2, true));
    //     	arg_desc.SetDependency(kArgMTMode,
    //                                CArgDescriptions::eRequires,
    //                                kArgNumThreads);
    //     }
    // ```
    /// `-mt_mode` is read and has no effect: LOSAT's output is the same at every thread count
    /// and for every value, and equals NCBI's `-num_threads 1` output (Owner 2026-10-13;
    /// CLI DECISIONS 5). NCBI's `-mt_mode 1` (and `0` when it chooses to split by queries)
    /// runs the queries in batches and writes the end of an outfmt 0 report or each XML
    /// document per batch (`blastn_app.cpp:134-170`, `blastx_app.cpp` `x_RunMTBySplitQuery`);
    /// that framing is not reproduced. The value must be 0 to 2 and `-num_threads` must be
    /// given with it, as in NCBI (clap's `requires` ignores a default value).
    #[arg(long, value_name = "int_value", value_parser = crate::blastinput::value_parsers::ncbi_mt_mode, requires = "num_threads", help = "Multi-thread mode: 0 auto, 1 split by queries, 2 split by database (accepted; the output does not depend on it)")]
    pub mt_mode: Option<u8>,
    #[arg(long)]
    pub ungapped: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1133-1143,1241-1247
    // ```c++
    //     if (m_DbTarget == eNucleotideDb) {
    //         arg_desc.SetCurrentGroup("PSI-TBLASTN options");
    //
    //         // PSI-tblastn checkpoint
    //         arg_desc.AddOptionalKey(kArgPSIInputChkPntFile, "psi_chkpt_file",
    //                                 "PSI-TBLASTN checkpoint file",
    //                                 CArgDescriptions::eInputFile);
    //         arg_desc.SetDependency(kArgPSIInputChkPntFile,
    //                                 CArgDescriptions::eExcludes,
    //                                 kArgRemote);
    // ...
    //     if (!m_IsDeltaBlast) {
    //         arg_desc.SetDependency(kArgPSIInputChkPntFile,
    //                                CArgDescriptions::eExcludes,
    //                                kArgQuery);
    //         arg_desc.SetDependency(kArgPSIInputChkPntFile,
    //                                CArgDescriptions::eExcludes,
    //                                kArgQueryLocation);
    //     }
    // ```
    // PSI-TBLASTN (line PS, item PS7, DECISIONS 35): the query is a PSSM. The exclusions are
    // errors of the argument parser, reported with LOSAT's text
    // (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES, exception 1; `-remote` is refused by the parser,
    // `cli.rs`). The file is read with the options (`check_options`, `in_pssm.rs`); the default
    // `-query -` is not read.
    #[arg(
        long = "in_pssm",
        value_name = "psi_chkpt_file",
        value_parser = ncbi_input_path(),
        conflicts_with_all = ["query", "query_loc"],
        help = "PSI-TBLASTN checkpoint file"
    )]
    pub in_pssm: Option<PathBuf>,
    /// The run's BLAST archive (outfmt 11), made by `run` and `run_local` before the search
    /// reads the arguments (`archive::handle_for`).
    #[arg(skip)]
    pub archive: Option<crate::report::archive::ArchiveHandle>,
    /// The run's search strategy export (`-export_search_strategy`), made by `run`
    /// (`strategy::StrategyExport::for_run`).
    #[arg(skip)]
    pub strategy: Option<crate::report::archive::strategy::StrategyExport>,
    /// Whether `taxdb.bti` is on NCBI's search path (`SeqDB_ResolveDbPath("taxdb.bti")`), asked
    /// by the command line before the search where a tabular field list shows a name of the
    /// taxonomy database (`algorithm::tblastn::tabular`); `None` elsewhere (`run_local`'s callers
    /// give no search path: not found, DECISIONS 7 of line CLI C2). Not a public option.
    #[arg(skip)]
    pub tabular_taxdb: Option<bool>,
    /// The taxonomy database of a `-subject` search's organism report (outfmt 18), made by the
    /// command line with `tabular_taxdb` (`common::database::SubjectTaxdb`); `None` elsewhere.
    /// Not a public option.
    #[arg(skip)]
    pub subject_taxdb: Option<crate::algorithm::common::database::SubjectTaxdb>,
    // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input.hpp:180-182
    // ```c++
    //     int GetLocalIdCounterInitValue() const { return m_LocalIdCounter; }
    //     ...
    //     void SetLocalIdCounterInitValue(int val) { m_LocalIdCounter = val; }
    // ```
    // The position in the whole input of the first of these records, as a reader's counter
    // starts there (`api::run_batches`): a query's `Query_N` in the messages. 0 on the command
    // line, whose records are the whole input.
    #[arg(skip)]
    pub first_query_ordinal: usize,
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:997-1056
// ```c++
// static int gcs[] = {1,2,3,4,5,6,9,10,11,12,13,14,15,16,
//                     21,22,23,24,25,26,27,28,29,30,31,33};
// return (genetic_codes.find(val) != genetic_codes.end());
// ```
// Stage A product decision admits 32 from gc.prt:340-347 for TBLASTN only.
/// `NStr::StringToInt(value)` of `CArgAllowGeneticCodeInteger` is
/// `value_parsers::genetic_code_number` (LOSAT's "invalid" message for the values outside 0..=255).
pub fn tblastn_genetic_code(value: &str) -> std::result::Result<u8, String> {
    let id = crate::blastinput::value_parsers::genetic_code_number(value)
        .ok_or("invalid genetic code ID")?;
    GeneticCode::try_from_id(id).map(|_| id)
}

/// The options of a TBLASTN search as NCBI's option handlers set them.
#[derive(Debug, Clone)]
pub struct ResolvedTblastnArgs {
    pub task: String,
    pub evalue: f64,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub xdrop_gap: f64,
    pub xdrop_gap_final: f64,
    /// `-xdrop_ungap` or 7 bits (`CBlastProteinOptionsHandle::SetInitialWordOptionsDefaults`,
    /// blast_prot_options.cpp:100-104).
    pub xdrop_ungap: f64,
    pub word_size: i32,
    /// Whether the compressed-alphabet lookup table is used (word sizes 5-7, tblastn-fast).
    pub compressed_lookup: bool,
    pub sum_stats: bool,
    pub db_gencode: u8,
    pub gapped: bool,
    pub max_intron_length: i32,
    pub seg: SegSpec,
    pub soft_masking: bool,
    pub lcase_masking: bool,
    /// The matrix name as typed (the report prints it so).
    pub matrix_name: String,
    pub threshold: f64,
    pub window_size: i32,
    /// The hit list size (`ReportCounts::hitlist_size`): the -max_target_seqs value, the
    /// larger of the numbers of descriptions and alignments of the pairwise report when
    /// one of them is given, the -num_alignments value of a tabular format, or 500.
    pub hitlist_size: usize,
    /// The numbers of descriptions and alignments of the pairwise report: the
    /// -max_target_seqs value, -num_descriptions and -num_alignments, or 500 and 250.
    pub num_descriptions: usize,
    pub num_alignments: usize,
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:263-265
    // ```c++
    //         formatter.SetLineLength(fmt_args->GetLineLength());
    //         formatter.SetHitsSortOption(fmt_args->GetHitsSortOption());
    //         formatter.SetHspsSortOption(fmt_args->GetHspsSortOption());
    // ```
    /// `-line_length`, `-sorthits` and `-sorthsps` of the pairwise report.
    pub formatting: crate::blastinput::report_counts::FormattingArgs,
    /// Whether -max_target_seqs was given (the text of LOSAT's limit on the hit list size).
    pub max_target_seqs_given: bool,
    /// The -max_hsps value, or 0 (every HSP).
    pub max_hsps: usize,
    /// The -qcov_hsp_perc value, or 0 (no filter).
    pub query_cov_hsp_perc: f64,
    /// The -culling_limit value (0: no culling).
    pub culling_limit: i32,
    pub composition_mode: crate::blastinput::app::CompositionMode,
    /// The `-comp_based_stats` value as typed.
    pub comp_based_stats: String,
    pub num_threads: usize,
    /// The `-query_loc` range, as given.
    pub query_range: Option<crate::blastinput::seq_range::SequenceRange>,
    /// `-in_pssm`: the query's PSSM (PSI-TBLASTN, `in_pssm.rs`).
    pub pssm: Option<std::sync::Arc<crate::algorithm::pssm_search::PsiSearch>>,
}

impl TblastnArgs {
    /// `-culling_limit`, 0 (no culling) when it is not given. The BLAST archive records
    /// whether it was given (`culling_limit`).
    pub fn culling_limit_value(&self) -> i32 {
        self.culling_limit.unwrap_or(0)
    }

    /// The `-db` value as NCBI keeps it (white space at its ends removed), `None` for a
    /// search of `-subject`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
    /// ```c++
    ///     if (args.Exist(kArgDb) && args[kArgDb]) {
    /// 	std::string local_dblist = NStr::TruncateSpaces( args[kArgDb].AsString() );
    ///
    ///         m_SearchDb.Reset(new CSearchDatabase( local_dblist,
    ///                                              mol_type));
    /// ```
    pub fn database_name(&self) -> Option<&str> {
        self.db.as_deref().map(crate::seqdb::db_argument_name)
    }

    /// The options as `check_options` sets them, without its warnings.
    pub fn resolve(&self) -> Result<ResolvedTblastnArgs> {
        self.check_options(&mut std::io::sink())
    }

    /// `check_options` for the output formats of a `run_local` call (their numbers).
    pub fn check_options_for_formats(
        &self,
        outfmts: &[i32],
        diagnostics: &mut dyn Write,
    ) -> Result<ResolvedTblastnArgs> {
        self.check_options_with(Some(outfmts), diagnostics)
    }

    /// NCBI's `CBlastAppArgs::SetOptions` for tblastn after the files are opened: the task's
    /// options, each argument handler in the order of `CTblastnAppArgs` (with the warnings
    /// that they post to `diagnostics`), and `Validate`, with NCBI's errors.
    ///
    /// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3631-3639
    /// ```c
    ///     NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args) {
    ///         (*arg)->ExtractAlgorithmOptions(args, opts);
    ///     }
    ///
    ///     m_IsUngapped = !opts.GetGappedMode();
    ///     try { retval->Validate(); }
    ///     catch (const CBlastException& e) {
    ///         NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
    ///     }
    /// ```
    pub fn check_options(&self, diagnostics: &mut dyn Write) -> Result<ResolvedTblastnArgs> {
        self.check_options_with(None, diagnostics)
    }

    fn check_options_with(
        &self,
        outfmts: Option<&[i32]>,
        diagnostics: &mut dyn Write,
    ) -> Result<ResolvedTblastnArgs> {
        use crate::blastinput::app::{
            formatting_handler_check, parse_comp_based_stats, parse_formatting_string,
            parse_seg_option,
        };
        use crate::blastinput::protein_options_check::{
            validate_protein_options, ProteinOptionsCheck,
        };
        use crate::stats::protein_options::{
            protein_gap_existence_extend_params, suggested_threshold, suggested_window_size,
            SuggestionProgram,
        };
        // NCBI reference: c++/src/algo/blast/api/blast_options_handle.cpp:436-447
        // ```c
        //     else if (!NStr::CompareNocase(task, "tblastn") ||
        //     		 !NStr::CompareNocase(task, "tblastn-fast"))
        //     {
        //     	CTBlastnOptionsHandle* opts =
        //     	            dynamic_cast<CTBlastnOptionsHandle*>
        //          (CBlastOptionsFactory::Create(eTblastn, locality));
        //     	if(task == "tblastn-fast") {
        //     		opts->SetWordSize(5);
        //             opts->SetOptions().SetLookupTableType(eCompressedAaLookupTable);
        //             opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP_FAST);
        //     	}
        //     	retval = opts;
        // ```
        // NCBI reference: c++/src/algo/blast/api/tblastn_options.cpp:54-85 (the handle's
        // defaults: threshold 13, word size 3, window 40, X-drops 15 and 25, BLOSUM62
        // 11/1, e-value 10, sum statistics, composition mode 2, genetic code 1)
        // ```c++
        // SetWordThreshold(BLAST_WORD_THRESHOLD_TBLASTN);
        // m_Opts->SetSumStatisticsMode();
        // m_Opts->SetCompositionBasedStats(eCompositionMatrixAdjust);
        // SetDbGeneticCode(BLAST_GENETIC_CODE);
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:138-149
        // ```c++
        //     if (args.Exist(kArgPSIInputChkPntFile) && args[kArgPSIInputChkPntFile]) {
        //        CPSIBlastOptionsHandle * rv(new CPSIBlastOptionsHandle(locality));
        //        rv->SetPSITblastnDefaults();
        //        if (args[kTask].AsString() == "blastx-fast") {
        //            rv->SetWordSize(6);
        // ...
        //        return CRef<CBlastOptionsHandle>(rv);
        //     }
        //     else {
        //         return x_CreateOptionsHandleWithTask(locality, args[kTask].AsString());
        //     }
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_options.cpp:69-75
        // ```c++
        // void CPSIBlastOptionsHandle::SetPSITblastnDefaults(void)
        // {
        //     m_Opts->SetProgram(ePSITblastn);
        //     m_Opts->SetSmithWatermanMode(false);
        //     m_Opts->SetDbGeneticCode(BLAST_GENETIC_CODE);
        //     m_Opts->SetWordThreshold(BLAST_WORD_THRESHOLD_TBLASTN);
        // }
        // ```
        // With `-in_pssm` (PSI-TBLASTN, DECISIONS 35 of line PS) the options are psiblast's handle
        // with tblastn's program, threshold 13 and genetic code: the defaults of blastp's protein
        // handle (word size 3, window 40, X-drops 15 and 25, BLOSUM62 11/1, e-value 10, the same
        // values as tblastn's) without sum statistics (the handle never calls
        // `SetSumStatisticsMode`). `-task` is checked and has no effect (the test names
        // `blastx-fast`, which is not a tblastn task).
        let psi = self.in_pssm.is_some();
        let fast = !psi && self.task == "tblastn-fast";
        let mut word_size: i32 = if fast { 5 } else { 3 };
        let mut compressed_lookup = fast;
        let mut threshold: f64 = if fast { 20.0 } else { 13.0 };
        let mut evalue = 10.0;
        let (mut gap_open, mut gap_extend) = (11, 1);
        let (mut xdrop_gap, mut xdrop_gap_final) = (15.0, 25.0);
        let mut sum_stats = !psi;
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:251-322 (the
        // generic search handler)
        // ```c
        //     if (args.Exist(kArgEvalue) && args[kArgEvalue]) {
        //         opt.SetEvalueThreshold(args[kArgEvalue].AsDouble());
        //     }
        //
        //     int gap_open=0, gap_extend=0;
        //     if (args.Exist(kArgMatrixName) && args[kArgMatrixName])
        //          BLAST_GetProteinGapExistenceExtendParams
        //              (args[kArgMatrixName].AsString().c_str(), &gap_open, &gap_extend);
        // ...
        //     if ( args.Exist(kArgWordSize) && args[kArgWordSize]) {
        //         if (m_QueryIsProtein && args[kArgWordSize].AsInteger() > 4){
        //            opt.SetLookupTableType(eCompressedAaLookupTable);
        //            opt.SetWordThreshold(19.3);
        // ...
        //     if (args.Exist(kArgSumStats) && args[kArgSumStats]) {
        //         opt.SetSumStatisticsMode(args[kArgSumStats].AsBoolean());
        //     }
        // ```
        if let Some(value) = self.evalue {
            evalue = value;
        }
        let (matrix_gap_open, matrix_gap_extend) = self
            .matrix
            .as_deref()
            .and_then(protein_gap_existence_extend_params)
            .unwrap_or((0, 0));
        match (self.gap_open, &self.matrix) {
            (Some(value), _) => gap_open = value,
            (None, Some(_)) => gap_open = matrix_gap_open,
            (None, None) => {}
        }
        match (self.gap_extend, &self.matrix) {
            (Some(value), _) => gap_extend = value,
            (None, Some(_)) => gap_extend = matrix_gap_extend,
            (None, None) => {}
        }
        if let Some(value) = self.xdrop_gap {
            xdrop_gap = value;
        }
        if let Some(value) = self.xdrop_gap_final {
            xdrop_gap_final = value;
        }
        if let Some(value) = self.word_size {
            if value > 4 {
                compressed_lookup = true;
                threshold = 19.3;
                if value > 5 {
                    threshold = 21.0;
                }
                if value > 6 {
                    threshold = 20.25;
                }
            }
            word_size = value;
            // NCBI reference: c++/src/algo/blast/api/blast_options_local_priv.hpp:612-619
            // ```c
            // CBlastOptionsLocal::SetWordSize(int ws)
            // {
            //     m_LutOpts->word_size = ws;
            //     if (m_LutOpts->lut_type == eCompressedAaLookupTable && ws <= 4)
            // 	m_LutOpts->lut_type = eAaLookupTable;
            //     else if (m_LutOpts->lut_type == eAaLookupTable && ws > 4)
            // 	m_LutOpts->lut_type = eCompressedAaLookupTable;
            // }
            // ```
            // (tblastn-fast with -word_size 3 or 4 uses the plain lookup table.)
            if compressed_lookup && value <= 4 {
                compressed_lookup = false;
            }
        }
        if let Some(value) = self.sum_stats {
            sum_stats = value;
        }
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:387-408 (the
        // filtering handler; the default -seg of tblastn is "12 2.2 2.5")
        // ```c
        //     if (args[kArgLookupTableMaskingOnly]) {
        //         opt.SetMaskAtHash(args[kArgLookupTableMaskingOnly].AsBoolean());
        //     }
        //
        //     vector<string> tokens;
        //
        //     try {
        //         if (m_QueryIsProtein && args[kArgSegFiltering]) {
        //             const string& seg_opts = args[kArgSegFiltering].AsString();
        // ```
        let soft_masking = self.soft_masking.unwrap_or(false);
        let seg = parse_seg_option(self.seg.as_deref().unwrap_or("12 2.2 2.5"))?;
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:586-623 (the word
        // threshold handler: the matrix's suggestion replaces only tblastn's default 13)
        // ```c
        //     } else if (program == eTblastn &&
        //                word_threshold != BLAST_WORD_THRESHOLD_TBLASTN) {
        //         retval = false;
        // ...
        //     if (args[kArgWordScoreThreshold]) {
        //         opt.SetWordThreshold(args[kArgWordScoreThreshold].AsDouble());
        //     } else if (s_IsDefaultWordThreshold(opt.GetProgram(),
        //                                         opt.GetWordThreshold())) {
        // ```
        let matrix_name = self
            .matrix
            .clone()
            .unwrap_or_else(|| "BLOSUM62".to_string());
        if let Some(value) = self.threshold {
            threshold = value;
        } else if threshold as i32 == 13 {
            threshold = suggested_threshold(SuggestionProgram::TranslatedSubject, &matrix_name);
        }
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:485-497
        // ```c
        // CWindowSizeArg::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& opt)
        // {
        //     if (args[kArgWindowSize]) {
        //         opt.SetWindowSize(args[kArgWindowSize].AsInteger());
        //     } else {
        //         int window = -1;
        //         BLAST_GetSuggestedWindowSize(opt.GetProgramType(),
        //                                      opt.GetMatrixName(),
        //                                      &window);
        // ```
        let window_size = self
            .window_size
            .unwrap_or_else(|| suggested_window_size(&matrix_name));
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1995-1999
        // ```c++
        //     // set the sequence range
        //     if (args.Exist(kArgQueryLocation) && args[kArgQueryLocation]) {
        //         m_Range = ParseSequenceRange(args[kArgQueryLocation].AsString(),
        //                                      "Invalid specification of query location");
        //     }
        // ```
        // The query options handler comes after the window size handler and before the
        // formatting handler (tblastn_args.cpp:44-135).
        let query_range = crate::blastinput::seq_range::parse_optional_range(
            self.query_loc.as_deref(),
            crate::blastinput::seq_range::RangeRole::Query,
        )?;
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2874-2886 and
        // 2960-2977 (the formatting handler: other programs' formats, the hit list size
        // and its warning)
        // ```c
        //     if(hitlist_size < 5){
        //    		ERR_POST(Warning << "Examining 5 or more matches is recommended");
        //     }
        // ```
        let choice = parse_formatting_string(&self.outfmt)?;
        formatting_handler_check(&choice, crate::blastinput::app::Program::Tblastn)?;
        // The hit list size and the counts of the pairwise report
        // (`blastinput::report_counts`), with the handler's warnings. The formats of a
        // `run_local` call are those of `outfmts`; the command line has the one `-outfmt`.
        let own_outfmt = [choice.number];
        let (counts, warnings) = crate::blastinput::report_counts::run_report_counts_with(
            "tblastn",
            outfmts.unwrap_or(&own_outfmt),
            self.max_target_seqs,
            self.num_descriptions,
            self.num_alignments,
            &self.formatting,
        )?;
        diagnostics.write_all(&warnings)?;
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:893-903 (the
        // composition-based statistics handler; `blastinput/app.rs`)
        // ```c
        //     if (args[kArgCompBasedStats]) {
        //         unique_ptr<bool> ungapped(args.Exist(kArgUngapped)
        //             ? new bool(args[kArgUngapped]) : 0);
        //         s_SetCompositionBasedStats(opt,
        //                                    args[kArgCompBasedStats].AsString(),
        //                                    args[kArgUseSWTraceback],
        //                                    ungapped.get());
        //     }
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:866-875
        // ```c
        //         if(program == ePSITblastn) {
        //             compo_mode = eNoCompositionBasedStats;
        //         }
        //
        //         if (ungapped && *ungapped && compo_mode != eNoCompositionBasedStats) {
        //             NCBI_THROW(CInputException, eInvalidInput,
        // ```
        // PSI-TBLASTN takes any value as mode 0, before the test of an ungapped search.
        let (composition_mode, _) =
            parse_comp_based_stats(&self.comp_based_stats, false, self.ungapped && !psi)?;
        let composition_mode = if psi {
            crate::blastinput::app::CompositionMode::NoCompositionBasedStats
        } else {
            composition_mode
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1338-1340
        // ```c++
        //     if (args.Exist(kArgPSIInputChkPntFile) && args[kArgPSIInputChkPntFile]) {
        //         CNcbiIstream& in = args[kArgPSIInputChkPntFile].AsInputFile();
        //         _ASSERT(m_Pssm.Empty());
        // ```
        // The PSSM is read by the last handler of tblastn's options (`CPsiBlastArgs`, after the
        // composition-based statistics, tblastn_args.cpp:128-130), before `Validate`.
        let pssm = match self.in_pssm.as_deref() {
            Some(path) => Some(std::sync::Arc::new(super::in_pssm::read(
                path,
                &matrix_name,
            )?)),
            None => None,
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1050-1055
        // ```c
        //     if (m_Target == eDatabase && args[kArgDbGeneticCode] &&
        //         (program == eTblastn || program == eTblastx) ) {
        //         opt.SetDbGeneticCode(args[kArgDbGeneticCode].AsInteger());
        //     }
        // ```
        // PSI-TBLASTN keeps the genetic code 1 of its handle; the value is only checked.
        let db_gencode = if psi {
            super::in_pssm::check_db_gencode(self.db_gencode)?;
            1
        } else {
            self.db_gencode
        };
        validate_protein_options(&ProteinOptionsCheck {
            gapped: !self.ungapped,
            gap_open,
            gap_extend,
            matrix_name: &matrix_name,
            threshold,
            word_size,
            compressed_lookup,
            hitlist_size: counts.hitlist_size,
            evalue,
            psiblast: false,
            psi_tblastn: psi,
            x_dropoff: self.xdrop_ungap,
        })?;
        // The gapped X-drops that NCBI's kernels cannot size are refused before the search
        // (`core::x_drop::check_gapped_x_dropoffs`), at the search's gapped Lambda: the matrix's,
        // or a PSSM query's `kbp_gap_psi` (`PsiSearch::gapped_karlin`, whose Lambda NCBI's X-drops
        // read, blast_parameters.c:455-458).
        // (A name that LOSAT has no matrix for is refused when the search starts.)
        let matrix = matrix_name.parse::<ScoringMatrix>();
        if let (false, Ok(matrix)) = (self.ungapped, matrix) {
            let gapped = crate::stats::lookup_protein_params(&crate::config::ProteinScoringSpec {
                matrix,
                gap_open,
                gap_extend,
            });
            // The final X-drop reaches the unscaled `ALIGN_EX` only without the
            // composition-based redo: mode 0, or IDENTITY's reset (`search_settings`).
            crate::core::x_drop::check_gapped_x_dropoffs(
                "TBLASTN",
                (xdrop_gap, xdrop_gap_final),
                crate::core::x_drop::gapped_x_dropoffs(
                    xdrop_gap,
                    xdrop_gap_final,
                    pssm.as_deref()
                        .map_or(gapped.lambda, |psi| psi.gapped_karlin(&gapped).lambda),
                ),
                (gap_open, gap_extend),
                composition_mode
                    == crate::blastinput::app::CompositionMode::NoCompositionBasedStats
                    || matrix_name.eq_ignore_ascii_case("identity"),
            )?;
        }
        Ok(ResolvedTblastnArgs {
            task: self.task.clone(),
            evalue,
            gap_open,
            gap_extend,
            xdrop_gap,
            xdrop_gap_final,
            xdrop_ungap: self.xdrop_ungap.unwrap_or(7.0),
            word_size,
            compressed_lookup,
            sum_stats,
            db_gencode,
            gapped: !self.ungapped,
            max_intron_length: self.max_intron_length.unwrap_or(0),
            seg,
            soft_masking,
            lcase_masking: self.lcase_masking,
            matrix_name,
            threshold,
            window_size,
            hitlist_size: counts.hitlist_size,
            num_descriptions: counts.num_descriptions,
            formatting: self.formatting,
            num_alignments: counts.num_alignments,
            max_target_seqs_given: self.max_target_seqs.is_some(),
            // NCBI blast_args.cpp:313-319: SetQueryCovHspPerc, SetMaxHspsPerSubject.
            max_hsps: self.max_hsps.unwrap_or(0),
            query_cov_hsp_perc: self.qcov_hsp_perc.unwrap_or(0.0),
            culling_limit: self.culling_limit_value(),
            composition_mode,
            comp_based_stats: self.comp_based_stats.clone(),
            num_threads: self.num_threads,
            query_range,
            pssm,
        })
    }
}

/// LOSAT's limits, checked where NCBI starts the search: the options that NCBI runs and
/// LOSAT's TBLASTN does not implement are rejected, and the others select the search.
pub(super) fn search_settings(args: &ResolvedTblastnArgs) -> Result<SearchSettings> {
    use crate::blastinput::app::CompositionMode;
    // IDENTITY with the compressed lookup table is refused where the search reaches its first
    // query batch (`check_identity_compressed_lookup_of`).
    if !args.gapped {
        anyhow::bail!("{}", crate::cli::run_refusal("TBLASTN", "ungapped"));
    }
    // NCBI reference: c++/src/algo/blast/core/blast_parameters.c:774-815 (uneven gap
    // linking of the HSPs, `longest_intron`)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:781-787
    // ```c
    //    if (params->do_sum_stats) {
    //       BlastLinkHSPParametersNew(program_number, gapped_calculation,
    //                                 &params->link_hsp_params);
    //
    //       if((Blast_QueryIsTranslated(program_number) ||
    // 	  Blast_SubjectIsTranslated(program_number)) &&
    // 	 program_number != eBlastTypeTblastx) {
    // ```
    // The largest intron is read by the linking of the HSPs only: without sum statistics (the
    // default of `-in_pssm`, or `-sum_stats false`) NCBI ignores it (DECISIONS 35 Q3 of line PS).
    if args.max_intron_length != 0 && args.sum_stats {
        anyhow::bail!(
            "-max_intron_length {} (uneven gap linking) is not supported by LOSAT's TBLASTN",
            args.max_intron_length
        );
    }
    // NCBI reference: c++/src/algo/blast/core/blast_traceback.c:1486-1499
    // ```c
    //     } else if (ext_params->options->compositionBasedStats > 0 ||
    //                ext_params->options->eTbackExt == eSmithWatermanTbck) {
    // ...
    //         retval =
    //                 Blast_RedoAlignmentCore_MT(program_number,
    // ...
    //     } else {
    // ```
    // The composition mode selects the ordinary traceback (mode 0) or the composition
    // redo, whose matrix adjustment rule modes 1 to 3 pick (`blast_adjust_scores`; item SC5 of
    // line NP).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:143-155
    // ```c++
    //     if (opts_memento->m_ScoringOpts->matrix &&
    //         strcasecmp("identity", opts_memento->m_ScoringOpts->matrix) == 0 &&
    //         opts_memento->m_ExtnOpts->compositionBasedStats > 0) {
    //
    //         opts_memento->m_ExtnOpts->compositionBasedStats = 0;
    // ...
    //         search_messages.AddMessageAllQueries(eBlastSevWarning, 0, os.str());
    //     }
    // ```
    // IDENTITY has no composition data: NCBI searches it in mode 0, whatever mode was asked,
    // with a warning for each query of the first query batch (`IDENTITY_CBS_RESET_MESSAGE`; the
    // later batches start in mode 0, `tblastn_app::search`). The reset is made with the score
    // block, before the HSP writer takes the preliminary hit list size
    // (`blast_aux_priv.cpp:175-229`), so that size is mode 0's too.
    let identity_cbs_reset = args.matrix_name.eq_ignore_ascii_case("identity")
        && args.composition_mode != CompositionMode::NoCompositionBasedStats;
    let compo_adjust_mode = match args.composition_mode {
        _ if identity_cbs_reset => BlastCompoAdjustMode::NoCompositionBasedStats,
        CompositionMode::NoCompositionBasedStats => BlastCompoAdjustMode::NoCompositionBasedStats,
        CompositionMode::CompositionBasedStats => BlastCompoAdjustMode::CompositionBasedStats,
        CompositionMode::CompositionMatrixAdjust => BlastCompoAdjustMode::CompositionMatrixAdjust,
        CompositionMode::CompoForceFullMatrixAdjust => BlastCompoAdjustMode::ForceFullMatrixAdjust,
    };
    // NCBI c++/src/algo/blast/blastinput/blast_args.cpp:3152-3187:
    // arg_desc.SetConstraint(kArgNumThreads, new CArgAllowValuesGreaterThanOrEqual(1));
    // NCBI c++/src/algo/blast/api/prelim_stage.cpp:145-147:
    // TBlastThreads the_threads(GetNumberOfThreads());
    crate::utils::threading::validate_threads(args.num_threads)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:409,3687
    // ```c
    //     *pbestEvalue = DBL_MAX;
    //                 if (best_evalue <= hitParams->options->expect_value) {
    // ```
    // A list whose HSPs all fail keeps the best e-value DBL_MAX, so an -evalue of DBL_MAX
    // or more (+inf, 1e999, 1.7976931348623157e308) puts an empty list into the hit list,
    // and NCBI's tblastn dies of SIGSEGV reading its first HSP (blast_hits.c:3266) when
    // the search keeps enough HSPs, and runs otherwise; LOSAT cannot tell beforehand, so
    // it rejects the value (decision D12 of docs/evidence/losat_web_e2e/AUTHORITY.md). A
    // smaller -evalue such as 1e308 runs.
    if args.evalue >= f64::MAX {
        anyhow::bail!(
            "an -evalue of DBL_MAX or more ({}), with which NCBI BLAST+'s tblastn crashes on some inputs, is not supported by LOSAT's TBLASTN",
            args.evalue
        );
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:43-70 (the
    // preliminary hit list size through NCBI's `int`; `blastn/hsp.rs`)
    // A hit list size whose preliminary size is not positive crashes NCBI.
    crate::algorithm::blastn::hsp::check_prelim_hitlist_size(
        args.hitlist_size,
        args.max_target_seqs_given,
        compo_adjust_mode.is_enabled(),
        true,
        "TBLASTN",
    )?;
    let matrix: Option<ScoringMatrix> = args.matrix_name.parse().ok();
    let threshold = args.threshold;
    let window = args.window_size;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:137 and
    // c++/src/algo/blast/core/aa_ungapped.c:214-224
    // ```c
    //         Boolean multiple_hits = (word_params->options->window_size > 0);
    // ...
    //     if (ewp->diag_table->multiple_hits) {
    //         status = s_BlastAaWordFinder_TwoHit(subject, query,
    // ...
    //     } else {
    //         status = s_BlastAaWordFinder_OneHit(subject, query,
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:910-936
    // ```c
    //       if (options->gapped_calculation && !Blast_ProgramIsRpsBlast(program_number)) {
    //           Int2 status=0;
    //           if ((status=Blast_KarlinBlkGappedLoadFromTables(NULL, options->gap_open,
    //              options->gap_extend, options->matrix, std_matrix_only)) != 0)
    // ```
    // Every matrix of NCBI's tables (the nine of blastp and tblastn, IDENTITY with them) and
    // every gap cost pair of its table, as `validate_protein_options` has checked them; the
    // threshold and window are the options' (the matrix's suggested values by default). The
    // standard table takes word sizes 2 to 4 and the compressed one 5 to 7; the word finder is
    // the two-hit one (a window above 0) or the one-hit one (a window of 0,
    // `XTwoHit::x_pairs_one_hit`).
    let profile_supported = (2..=7).contains(&args.word_size);
    let Some(matrix) = matrix.filter(|_| profile_supported) else {
        anyhow::bail!(
            "the matrix {} with word size {} is not supported by LOSAT's TBLASTN; it implements word sizes 2 to 7",
            args.matrix_name,
            args.word_size
        );
    };
    Ok(SearchSettings {
        compo_adjust_mode,
        identity_cbs_reset,
        scoring: LocalStageDScoring {
            matrix,
            gap_open: args.gap_open,
            gap_extend: args.gap_extend,
            word_size: args.word_size as usize,
            threshold,
            window,
            gap_xdrop_bits: args.xdrop_gap,
            final_xdrop_bits: args.xdrop_gap_final,
            word_xdrop_bits: args.xdrop_ungap,
        },
    })
}

/// The checks of the options alone, without the inputs (the `validate` of web ABI v2):
/// NCBI's processing of the options (`TblastnArgs::check_options`) and LOSAT's limits.
pub fn check_options(args: &TblastnArgs) -> Result<()> {
    // The host has no database files (`run_local`).
    if args.db.is_some() {
        anyhow::bail!("-db (a database search) is not supported in the browser yet");
    }
    // The host has no PSSM files either (line PS, item PSW1).
    if args.in_pssm.is_some() {
        anyhow::bail!("-in_pssm (a PSSM query) is not supported in the browser yet");
    }
    let resolved = args.resolve()?;
    search_settings(&resolved)?;
    // Without the inputs: IDENTITY with the compressed lookup table, which the search refuses where
    // it searches its first query batch.
    check_identity_compressed_lookup_of(&resolved)
}

/// `protein_options_check::check_identity_compressed_lookup` of a TBLASTN search.
pub(super) fn check_identity_compressed_lookup_of(args: &ResolvedTblastnArgs) -> Result<()> {
    crate::blastinput::protein_options_check::check_identity_compressed_lookup(
        &args.matrix_name,
        args.compressed_lookup,
        args.word_size,
        "tblastn",
        "TBLASTN",
    )
}

// NCBI reference: c++/src/algo/blast/api/tblastn_options.cpp:54-58
// ```c++
// CTBlastnOptionsHandle::SetLookupTableDefaults()
// {
//     CBlastProteinOptionsHandle::SetLookupTableDefaults();
//     SetWordThreshold(BLAST_WORD_THRESHOLD_TBLASTN);
// }
// ```
// The options handle keeps the composition mode and the scoring choices that the
// search reads; these are the ones the Stage D port uses.
/// The search settings that the options select.
pub(super) struct SearchSettings {
    /// The composition mode of the search (mode 0 for IDENTITY, `identity_cbs_reset`).
    pub(super) compo_adjust_mode: BlastCompoAdjustMode,
    /// IDENTITY with a composition mode above 0, which NCBI searches in mode 0 with a warning
    /// for each query of the first query batch (`IDENTITY_CBS_RESET_MESSAGE`).
    pub(super) identity_cbs_reset: bool,
    pub(super) scoring: LocalStageDScoring,
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::cli::{try_parse_from, Cli, Commands};

    // TBLASTN's `-query_gencode` and `-db_gencode`: `genetic_code_number`, then the codes LOSAT
    // builds (32 included, Stage A); "invalid genetic code ID" for what is not a number of 0 to
    // 255 written without a minus sign.
    #[test]
    fn tblastn_genetic_codes_are_read_as_ncbi_reads_them() {
        for value in [
            "-0",
            "-5",
            "256",
            "2147483648",
            "",
            " 5",
            "5 ",
            "5.0",
            "0x",
            "0x5",
        ] {
            assert_eq!(
                tblastn_genetic_code(value),
                Err("invalid genetic code ID".to_string()),
                "{value:?}"
            );
        }
        for (value, expected) in [("+5", 5), ("005", 5), ("11", 11), ("33", 33)] {
            assert_eq!(tblastn_genetic_code(value), Ok(expected), "{value}");
        }
        for value in ["0", "7", "8", "255"] {
            assert!(tblastn_genetic_code(value).is_err(), "{value}");
        }
    }

    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1938-1942
    // arg_desc.AddFlag(kArgUseLCaseMasking,
    //     "Use lower case filtering in query and subject sequence(s)?", true);
    #[test]
    fn lcase_masking_flag_is_explicit() {
        assert!(!parse(&[]).unwrap().lcase_masking);
        assert!(parse(&["-lcase_masking"]).unwrap().lcase_masking);
    }

    fn parse(extra: &[&str]) -> std::result::Result<TblastnArgs, clap::Error> {
        let mut argv = vec!["losat", "tblastn", "-query", "q.faa", "-subject", "s.fna"];
        argv.extend_from_slice(extra);
        let cli: Cli = try_parse_from(argv)?;
        let Commands::Tblastn(args) = cli.command else {
            unreachable!()
        };
        Ok(args)
    }

    fn error(extra: &[&str]) -> String {
        parse(extra).unwrap().resolve().unwrap_err().to_string()
    }

    fn limit(extra: &[&str]) -> String {
        check_options(&parse(extra).unwrap())
            .unwrap_err()
            .to_string()
    }

    // NCBI blast_args.cpp:2894-2973 (`blastinput::report_counts`).
    #[test]
    fn num_descriptions_and_num_alignments_set_the_counts() {
        let counts = |extra: &[&str]| {
            let resolved = parse(extra).unwrap().resolve().unwrap();
            (
                resolved.hitlist_size,
                resolved.num_descriptions,
                resolved.num_alignments,
            )
        };
        assert_eq!(counts(&[]), (500, 500, 250));
        assert_eq!(counts(&["-num_descriptions", "3"]), (250, 3, 250));
        assert_eq!(counts(&["-num_alignments", "2"]), (500, 500, 2));
        assert_eq!(counts(&["-max_target_seqs", "7"]), (7, 7, 7));
        assert_eq!(counts(&["-outfmt", "6", "-num_alignments", "2"]), (2, 2, 2));
        assert_eq!(
            error(&["-num_descriptions", "0", "-num_alignments", "0"]),
            "BLAST query/options error: No hits are being saved\nPlease refer to the BLAST+ user manual.\n"
        );
    }

    // NCBI tblastn_options.cpp:54-85 and blast_prot_options.cpp:85-141:
    // SetWordThreshold(BLAST_WORD_THRESHOLD_TBLASTN);
    // m_Opts->SetSumStatisticsMode(); SetDbGeneticCode(BLAST_GENETIC_CODE);
    #[test]
    fn ordinary_defaults_and_code32_are_accepted() {
        let args = parse(&[]).unwrap();
        assert_eq!(args.outfmt, "0");
        assert_eq!(args.comp_based_stats, "2");
        let resolved = args.resolve().unwrap();
        assert_eq!(resolved.task, "tblastn");
        assert_eq!(resolved.db_gencode, 1);
        assert_eq!(resolved.word_size, 3);
        assert_eq!(resolved.threshold, 13.0);
        assert_eq!(resolved.window_size, 40);
        assert_eq!((resolved.gap_open, resolved.gap_extend), (11, 1));
        assert_eq!(resolved.evalue, 10.0);
        assert_eq!(resolved.matrix_name, "BLOSUM62");
        assert!(resolved.sum_stats && !resolved.soft_masking);
        assert_eq!(resolved.max_intron_length, 0);
        assert_eq!(resolved.hitlist_size, 500);
        check_options(&args).unwrap();
        let resolved = parse(&["-soft_masking", "true", "-sum_stats", "off"])
            .unwrap()
            .resolve()
            .unwrap();
        assert!(resolved.soft_masking && !resolved.sum_stats);
        assert!(parse(&["-sum_stats", "maybe"]).is_err());
        // NCBI blast_args.cpp:586-623: the matrix's suggestion (+2 for a translated
        // subject) replaces only the default 13.
        let resolved = parse(&["-matrix", "BLOSUM45"]).unwrap().resolve().unwrap();
        assert_eq!(
            (
                resolved.threshold,
                resolved.window_size,
                resolved.gap_open,
                resolved.gap_extend
            ),
            (16.0, 60, 14, 2)
        );
        for id in [
            1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29,
            30, 31, 32, 33,
        ] {
            let id_string = id.to_string();
            assert_eq!(parse(&["-db_gencode", &id_string]).unwrap().db_gencode, id);
        }
    }

    // NCBI blast_args.cpp:997-1056 has 26 CLI values; gc.prt:340-347
    // defines ID 32 for the TBLASTN-only product decision.
    #[test]
    fn invalid_codes_are_explicit_errors() {
        for id in ["0", "7", "8", "17", "20", "34", "255", "256", "abc", "-1"] {
            let error = parse(&["-db_gencode", id]).unwrap_err().to_string();
            assert!(error.contains("genetic code"), "{id}: {error}");
        }
    }

    // NCBI's checks of the options (blast_options.c BLAST_ValidateOptions, the option
    // handlers of blast_args.cpp) and LOSAT's limits, where NCBI starts the search.
    #[test]
    fn ncbi_checks_and_losat_limits_are_explicit() {
        assert!(error(&["-matrix", "BAD"]).contains("BAD is not a supported matrix"));
        assert!(error(&["-gapopen", "1", "-gapextend", "1"])
            .contains("Gap existence and extension values of 1 and 1 not supported for BLOSUM62"));
        assert!(error(&["-matrix", "IDENTITY", "-word_size", "6"])
            .contains("Word size larger than 5 is not supported for the identity scoring matrix"));
        assert!(error(&["-evalue", "0"])
            .contains("expect value or cutoff score must be greater than zero"));
        assert!(error(&["-word_size", "8"]).contains("Word-size must be less than 8"));
        assert!(error(&["-threshold", "0"]).contains("Non-zero threshold required"));
        assert!(error(&["-ungapped"]).contains("Composition-adjusted searched"));
        assert!(error(&["-seg", "12 2.2"]).contains("Invalid number of arguments"));
        assert!(parse(&["-word_size", "1"]).is_err());
        // NCBI blast_args.cpp:255-275 and blast_options.c:1174-1235: a matrix without gap costs
        // takes its best pair, and its suggested threshold (+2) and window.
        let resolved = parse(&["-matrix", "PAM30"]).unwrap().resolve().unwrap();
        assert_eq!(
            (
                resolved.gap_open,
                resolved.gap_extend,
                resolved.threshold,
                resolved.window_size
            ),
            (9, 1, 18.0, 15)
        );
        for extra in [
            &["-ungapped", "-comp_based_stats", "F"][..],
            &["-max_intron_length", "100", "-sum_stats", "true"],
            &["-max_target_seqs", "1073741799"],
            &["-comp_based_stats", "0", "-max_target_seqs", "2147483598"],
        ] {
            assert!(
                limit(extra).contains("not supported by LOSAT's TBLASTN"),
                "{extra:?}"
            );
        }
        // NCBI's tblastn crashes with IDENTITY's compressed lookup table.
        for extra in [
            &["-task", "tblastn-fast", "-matrix", "IDENTITY"][..],
            &["-matrix", "identity", "-word_size", "5"],
        ] {
            assert!(
                limit(extra).contains("with which NCBI BLAST+'s tblastn crashes"),
                "{extra:?}"
            );
        }
        // Every matrix of NCBI's tables with every gap cost pair of its table, word sizes 2 to 4,
        // and composition modes 0 and 2 (series SC, SC1) and 1 and 3 (SC5).
        for extra in [
            &["-matrix", "PAM30"][..],
            &["-comp_based_stats", "1"],
            &["-comp_based_stats", "3"],
            &["-matrix", "PAM30", "-comp_based_stats", "1"],
            &[
                "-task",
                "tblastn-fast",
                "-matrix",
                "BLOSUM45",
                "-comp_based_stats",
                "3",
            ],
            &["-matrix", "pam250", "-comp_based_stats", "0"],
            &["-matrix", "BLOSUM45", "-gapopen", "13", "-gapextend", "3"],
            &["-matrix", "BLOSUM45", "-word_size", "3"],
            &["-matrix", "BLOSUM80", "-gapopen", "25", "-gapextend", "2"],
            &["-gapopen", "13", "-gapextend", "1"],
            &["-word_size", "2"],
            &["-word_size", "4", "-matrix", "PAM70"],
            &[
                "-matrix",
                "BLOSUM90",
                "-threshold",
                "11",
                "-window_size",
                "60",
            ],
            // The compressed lookup table with every matrix but IDENTITY (series SC, SC11).
            &["-task", "tblastn-fast", "-matrix", "BLOSUM80"],
            &[
                "-matrix",
                "pam30",
                "-word_size",
                "5",
                "-comp_based_stats",
                "0",
            ],
            &["-matrix", "BLOSUM45", "-word_size", "6", "-threshold", "15"],
            &[
                "-task",
                "tblastn-fast",
                "-matrix",
                "PAM250",
                "-word_size",
                "7",
            ],
            // The one-hit word finder (line RG, X2; aa_ungapped.c:214-224).
            &["-window_size", "0"],
            &["-word_size", "5", "-window_size", "0"],
            &["-task", "tblastn-fast", "-window_size", "0"],
        ] {
            check_options(&parse(extra).unwrap()).unwrap();
        }
        // NCBI blast_args.cpp:606-620: the matrix's suggested threshold replaces tblastn's 13
        // only, not the compressed table's 19.3, 21 or 20.25 (nor tblastn-fast's 20); the
        // window is the matrix's (BLOSUM45's 60).
        let settings = search_settings(
            &parse(&["-matrix", "BLOSUM45", "-word_size", "6"])
                .unwrap()
                .resolve()
                .unwrap(),
        )
        .unwrap();
        assert_eq!(
            (
                settings.scoring.word_size,
                settings.scoring.threshold,
                settings.scoring.window
            ),
            (6, 21.0, 60)
        );
        // NCBI blast_args.cpp:836-879: the typed mode is the search's (`compo_adjust_mode`).
        for (value, mode) in [
            ("0", BlastCompoAdjustMode::NoCompositionBasedStats),
            ("F", BlastCompoAdjustMode::NoCompositionBasedStats),
            ("1", BlastCompoAdjustMode::CompositionBasedStats),
            ("2", BlastCompoAdjustMode::CompositionMatrixAdjust),
            ("T", BlastCompoAdjustMode::CompositionMatrixAdjust),
            ("3", BlastCompoAdjustMode::ForceFullMatrixAdjust),
        ] {
            let settings = search_settings(
                &parse(&["-comp_based_stats", value])
                    .unwrap()
                    .resolve()
                    .unwrap(),
            )
            .unwrap();
            assert_eq!(settings.compo_adjust_mode, mode, "{value}");
        }
        // NCBI setup_factory.cpp:143-155: IDENTITY is searched in mode 0 with a warning (for the
        // first query batch), whatever mode is asked.
        for (extra, reset) in [
            (&["-matrix", "IDENTITY"][..], true),
            (&["-matrix", "identity", "-comp_based_stats", "1"], true),
            (&["-matrix", "IDENTITY", "-comp_based_stats", "3"], true),
            (&["-matrix", "IDENTITY", "-comp_based_stats", "0"], false),
            (&["-matrix", "BLOSUM62"], false),
        ] {
            let settings = search_settings(&parse(extra).unwrap().resolve().unwrap()).unwrap();
            assert_eq!(
                (
                    settings.identity_cbs_reset,
                    settings.compo_adjust_mode.is_enabled()
                ),
                (reset, !reset && extra.len() == 2),
                "{extra:?}"
            );
        }
        // NCBI api/blast_options_handle.cpp:436-447: tblastn-fast sets word size 5, the
        // compressed lookup table and BLAST_WORD_THRESHOLD_BLASTP_FAST (20); with BLOSUM62 11/1
        // the table of word size 3, 5, 6 or 7 takes any threshold and window size.
        let fast = parse(&["-task", "tblastn-fast"]).unwrap();
        let resolved = fast.resolve().unwrap();
        assert_eq!(
            (
                resolved.word_size,
                resolved.compressed_lookup,
                resolved.threshold
            ),
            (5, true, 20.0)
        );
        check_options(&fast).unwrap();
        for extra in [
            &[
                "-task",
                "tblastn-fast",
                "-word_size",
                "5",
                "-window_size",
                "15",
                "-threshold",
                "100",
            ][..],
            &["-word_size", "6"],
            &["-word_size", "7", "-threshold", "15"],
            &["-task", "tblastn-fast", "-word_size", "3"],
            &["-word_size", "3", "-window_size", "15", "-threshold", "100"],
            &["-threshold", "13.5", "-window_size", "1"],
        ] {
            check_options(&parse(extra).unwrap()).unwrap();
        }
        let settings =
            search_settings(&parse(&["-word_size", "7"]).unwrap().resolve().unwrap()).unwrap();
        assert_eq!(
            (settings.scoring.word_size, settings.scoring.threshold),
            (7, 20.25)
        );
        // NCBI blast_args.cpp:834-892: the first character chooses the mode (any other is
        // mode 0); a `u` asks for unified P-values only for blastp.
        for value in ["bogus", "0", "F", "2u", "D", "t"] {
            check_options(&parse(&["-comp_based_stats", value]).unwrap()).unwrap();
        }
        // NCBI blast_hits.c:43-70 through an int: 2^30 without composition statistics
        // keeps a preliminary hit list of 10 subjects.
        check_options(
            &parse(&["-comp_based_stats", "0", "-max_target_seqs", "1073741824"]).unwrap(),
        )
        .unwrap();
        // NCBI blast_args.cpp:313-319: -qcov_hsp_perc and -max_hsps.
        let resolved = parse(&["-max_hsps", "2", "-qcov_hsp_perc", "50.5"])
            .unwrap()
            .resolve()
            .unwrap();
        assert_eq!((resolved.max_hsps, resolved.query_cov_hsp_perc), (2, 50.5));
        let resolved = parse(&[]).unwrap().resolve().unwrap();
        assert_eq!((resolved.max_hsps, resolved.query_cov_hsp_perc), (0, 0.0));
        for extra in [
            &["-max_hsps", "0"][..],
            &["-qcov_hsp_perc", "100.5"],
            &["-qcov_hsp_perc", "nan"],
        ] {
            assert!(parse(extra).is_err(), "{extra:?}");
        }
        for extra in [
            &["-db_soft_mask", "1"][..],
            &["-remote"],
            &["-use_sw_tback"],
        ] {
            let error = crate::cli::render_message(&parse(extra).unwrap_err());
            assert!(
                error.contains("is not supported by LOSAT's TBLASTN"),
                "{extra:?}: {error}"
            );
        }
    }

    /// The tblastn arguments `-in_pssm <psiblast's frozen PSSM> -subject s.fna <extra>`.
    fn parse_pssm(extra: &[&str]) -> std::result::Result<TblastnArgs, clap::Error> {
        let pssm = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/db_search/db.ps.asn.rnd.6.asn"
        );
        let mut argv = vec!["losat", "tblastn", "-in_pssm", pssm, "-subject", "s.fna"];
        argv.extend_from_slice(extra);
        let cli: Cli = try_parse_from(argv)?;
        let Commands::Tblastn(args) = cli.command else {
            unreachable!()
        };
        Ok(args)
    }

    // NCBI tblastn_args.cpp:138-149, psiblast_options.cpp:69-75 and blast_args.cpp:866-868,
    // 1050-1055: PSI-TBLASTN is psiblast's handle with tblastn's threshold and genetic code 1, no
    // sum statistics, composition-based statistics 0 whatever `-comp_based_stats` says, and
    // `-task tblastn-fast` without effect (line PS, DECISIONS 35).
    #[test]
    fn in_pssm_takes_psi_tblastns_options() {
        let resolved = parse_pssm(&[]).unwrap().resolve().unwrap();
        assert!(resolved.pssm.is_some());
        assert!(!resolved.sum_stats && !resolved.compressed_lookup);
        assert_eq!((resolved.word_size, resolved.threshold), (3, 13.0));
        assert_eq!(resolved.db_gencode, 1);
        assert_eq!(
            resolved.composition_mode,
            crate::blastinput::app::CompositionMode::NoCompositionBasedStats
        );
        for value in ["0", "1", "2", "3", "D", "t", "x", "", "2u"] {
            let resolved = parse_pssm(&["-comp_based_stats", value])
                .unwrap()
                .resolve()
                .unwrap();
            assert_eq!(
                resolved.composition_mode,
                crate::blastinput::app::CompositionMode::NoCompositionBasedStats,
                "{value:?}"
            );
        }
        let fast = parse_pssm(&["-task", "tblastn-fast"])
            .unwrap()
            .resolve()
            .unwrap();
        assert_eq!(
            (fast.word_size, fast.threshold, fast.compressed_lookup),
            (3, 13.0, false)
        );
        let resolved = parse_pssm(&["-db_gencode", "4", "-sum_stats", "true"])
            .unwrap()
            .resolve()
            .unwrap();
        assert_eq!((resolved.db_gencode, resolved.sum_stats), (1, true));
        // The matrix's suggestion replaces the default threshold 13 (+2 for a translated subject).
        let resolved = parse_pssm(&["-matrix", "BLOSUM45"])
            .unwrap()
            .resolve()
            .unwrap();
        assert_eq!((resolved.threshold, resolved.window_size), (16.0, 60));
        // A matrix query keeps tblastn's handle.
        assert!(parse(&[]).unwrap().resolve().unwrap().pssm.is_none());
    }

    // NCBI blast_options.c:910-936 (`std_matrix_only` for a program that is neither blastp nor
    // tblastn) and 1358-1363 (the word size of PSI-TBLASTN); the code 32 is not NCBI's
    // (blast_args.cpp:1011-1015); `-query` and `-query_loc` exclude `-in_pssm`
    // (blast_args.cpp:1241-1247).
    #[test]
    fn in_pssm_errors_are_ncbis_and_losats() {
        let error = |extra: &[&str]| {
            parse_pssm(extra)
                .unwrap()
                .resolve()
                .unwrap_err()
                .to_string()
        };
        for size in ["5", "6", "7", "8"] {
            assert_eq!(
                error(&["-word_size", size]),
                "BLAST query/options error: Word-size must be less than 6 for protein comparison\nPlease refer to the BLAST+ user manual.\n",
                "{size}"
            );
        }
        assert!(error(&["-task", "tblastn-fast", "-word_size", "6"])
            .contains("Word-size must be less than 6 for protein comparison"));
        let identity = error(&["-matrix", "IDENTITY"]);
        assert!(
            identity.starts_with(
                "BLAST query/options error: IDENTITY is not a supported matrix, supported matrices are:\nBLOSUM80 \n"
            ),
            "{identity}"
        );
        assert!(!identity[40..].contains("IDENTITY"), "{identity}");
        assert!(error(&["-db_gencode", "32"]).contains("-db_gencode 32 with -in_pssm"));
        assert!(error(&["-threshold", "0"]).contains("Non-zero threshold required"));
        let pssm = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/db_search/db.ps.asn.rnd.6.asn"
        );
        for extra in [&["-query", "q.faa"][..], &["-query_loc", "1-10"]] {
            let mut argv = vec!["losat", "tblastn", "-in_pssm", pssm, "-subject", "s.fna"];
            argv.extend_from_slice(extra);
            assert!(try_parse_from::<Cli, _, _>(argv).is_err(), "{extra:?}");
        }
    }

    // NCBI blast_parameters.c:781-787: the largest intron is read by the linking only (DECISIONS
    // 35 Q3 of line PS: accepted without `-sum_stats true`, refused with it as before).
    #[test]
    fn the_largest_intron_is_refused_with_sum_statistics_only() {
        let settings = |args: TblastnArgs| search_settings(&args.resolve().unwrap()).map(|_| ());
        assert!(settings(parse_pssm(&["-max_intron_length", "100"]).unwrap()).is_ok());
        assert!(
            settings(parse(&["-max_intron_length", "100", "-sum_stats", "false"]).unwrap()).is_ok()
        );
        for args in [
            parse_pssm(&["-max_intron_length", "100", "-sum_stats", "true"]).unwrap(),
            parse(&["-max_intron_length", "100"]).unwrap(),
        ] {
            assert!(settings(args)
                .unwrap_err()
                .to_string()
                .contains("-max_intron_length 100 (uneven gap linking) is not supported"));
        }
    }

    // blast_args.cpp:276-278 and blast_options.c:608-618: `-xdrop_ungap` reaches the search's
    // options (7 bits by default), one of 0 or below is NCBI's options error, and a gapped X-drop
    // whose raw score NCBI's kernels cannot size is refused (`core::x_drop`).
    #[test]
    fn the_x_drops_reach_the_search_or_ncbis_errors() {
        let resolved = parse(&[]).unwrap().resolve().unwrap();
        assert_eq!(
            (
                resolved.xdrop_ungap,
                resolved.xdrop_gap,
                resolved.xdrop_gap_final
            ),
            (7.0, 15.0, 25.0)
        );
        let resolved = parse(&["-xdrop_ungap", "3.5", "-xdrop_gap", "1e300"])
            .unwrap()
            .resolve()
            .unwrap();
        assert_eq!((resolved.xdrop_ungap, resolved.xdrop_gap), (3.5, 1e300));
        assert_eq!(
            search_settings(&resolved).unwrap().scoring.word_xdrop_bits,
            3.5
        );
        for value in ["0", "-1", "-0.0", "-inf"] {
            let error = parse(&["-xdrop_ungap", value])
                .unwrap()
                .resolve()
                .unwrap_err();
            assert!(
                format!("{error:#}").contains("x_dropoff must be greater than zero"),
                "{value}: {error:#}"
            );
        }
        // 4.2e8 bits are 1.09e9 raw at BLOSUM62 11/1's Lambda, past the kernels' 2^30 - 103: the
        // preliminary X-drop at every mode, the final one only without the redo (mode 0, or
        // IDENTITY's reset).
        for extra in [
            &["-xdrop_gap", "4.2e8"][..],
            &["-xdrop_gap_final", "4.2e8", "-comp_based_stats", "0"],
        ] {
            let error = parse(extra).unwrap().resolve().unwrap_err();
            assert!(
                format!("{error:#}").starts_with(&format!("an {} of 420000000.0 bits", extra[0])),
                "{extra:?}: {error:#}"
            );
        }
        for cbs in ["1", "2", "3"] {
            assert!(
                parse(&["-xdrop_gap_final", "4.2e8", "-comp_based_stats", cbs])
                    .unwrap()
                    .resolve()
                    .is_ok()
            );
        }
    }
}

#[cfg(test)]
#[path = "tblastn_culling_tests.rs"]
mod tblastn_culling_tests;
