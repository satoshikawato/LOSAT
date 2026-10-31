//! Command-line arguments for BLASTP

use anyhow::{bail, Result};
use clap::Args;
use std::path::PathBuf;

use crate::config::{ProteinScoringSpec, ScoringMatrix};
use crate::utils::seg::SegParams;

pub use crate::blastinput::value_parsers::SegSpec as BlastpSegSpec;
use crate::blastinput::value_parsers::*;

/// The message for `-task blastp-short`, which the search entrance printed for the task's
/// word size (`blast_engine::check_options`) and which is now printed when the arguments
/// are parsed.
const BLASTP_SHORT_UNSUPPORTED: &str =
    "-word_size 2 (the lookup table) is not supported by LOSAT's BLASTP";

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:825-876
// ```c
// switch (comp_stat_string[0]) {
// case '0': case 'F': case 'f':
//     compo_mode = eNoCompositionBasedStats;
//     break;
// case '1':
//     compo_mode = eCompositionBasedStats;
//     break;
// case 'D': case 'd':
//     ...
//     compo_mode = eCompositionMatrixAdjust;
//     break;
// case '2':
//     compo_mode = eCompositionMatrixAdjust;
//     break;
// case '3':
//     compo_mode = eCompoForceFullMatrixAdjust;
//     break;
// case 'T': case 't':
//     compo_mode = ... eCompositionMatrixAdjust;
//     break;
// }
// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlastpCompositionMode {
    NoCompositionBasedStats,
    CompositionBasedStats,
    CompositionMatrixAdjust,
    ForceFullMatrixAdjust,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:877-886
// ```c
// opt.SetCompositionBasedStats(compo_mode);
// if (program == eBlastp &&
//     compo_mode != eNoCompositionBasedStats &&
//     tolower(comp_stat_string[1]) == 'u') {
//     opt.SetUnifiedP(1);
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:44-56
// ```c++
// CPsiBlastAppArgs::CPsiBlastAppArgs()
// {
//     bool const kQueryIsProtein = true;
//     bool const kFilterByDefault = false;
//     static const string kProgram("psiblast");
//     CRef<IBlastCmdLineArgs> arg;
//     arg.Reset(new CProgramDescriptionArgs(kProgram,
//                                           "Position-Specific Iterated BLAST"));
//     m_Args.push_back(arg);
//     m_ClientId = kProgram + " " + CBlastVersion().Print();
//
//     static const char kDefaultTask[] = "psiblast";
//     SetTask(kDefaultTask);
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:418-421
// ```c++
//         // Allow the fasta reader to complain on invalid sequence input
//         SetDiagPostLevel(eDiag_Warning);
//         SetDiagPostPrefix("psiblast");
//         SetDiagHandler(&m_bah, false);
// ```
/// The NCBI application whose search the blastp engine runs: blastp, or the first round of
/// psiblast (`algorithm/psiblast`), whose program is `ePSIBlast` (`eBlastTypePsiBlast` in the
/// core) with the same search as blastp until a PSSM exists.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BlastpApp {
    #[default]
    Blastp,
    Psiblast,
}

impl BlastpApp {
    /// The program's name in NCBI's messages (`Warning: [psiblast] ...`) and reports
    /// (`# PSIBLAST 2.17.0+` in upper case).
    pub fn ncbi_name(self) -> &'static str {
        match self {
            BlastpApp::Blastp => "blastp",
            BlastpApp::Psiblast => "psiblast",
        }
    }

    /// The name in LOSAT's rejections ("... is not supported by LOSAT's PSIBLAST").
    pub fn name(self) -> &'static str {
        match self {
            BlastpApp::Blastp => "BLASTP",
            BlastpApp::Psiblast => "PSIBLAST",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlastpCompBasedStats {
    pub mode: BlastpCompositionMode,
    pub unified_p: bool,
}

impl BlastpCompBasedStats {
    #[inline]
    pub fn is_enabled(&self) -> bool {
        self.mode != BlastpCompositionMode::NoCompositionBasedStats
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:825-886
    // ```c
    // switch (comp_stat_string[0]) {
    // case '0': ... case '1': ... case '2': ... case '3': ...
    // }
    // if (program == eBlastp && compo_mode != eNoCompositionBasedStats &&
    //     tolower(comp_stat_string[1]) == 'u') {
    //     opt.SetUnifiedP(1);
    // }
    // ```
    pub fn to_ncbi_cli_string(&self) -> String {
        let mode = match self.mode {
            BlastpCompositionMode::NoCompositionBasedStats => '0',
            BlastpCompositionMode::CompositionBasedStats => '1',
            BlastpCompositionMode::CompositionMatrixAdjust => '2',
            BlastpCompositionMode::ForceFullMatrixAdjust => '3',
        };
        if self.unified_p && self.mode != BlastpCompositionMode::NoCompositionBasedStats {
            format!("{mode}u")
        } else {
            mode.to_string()
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:838-886
// ```c
// switch (comp_stat_string[0]) {
// ...
// }
// ...
// if (program == eBlastp &&
//     compo_mode != eNoCompositionBasedStats &&
//     tolower(comp_stat_string[1]) == 'u') {
//     opt.SetUnifiedP(1);
// }
// ```
pub fn parse_comp_based_stats(value: &str) -> Result<BlastpCompBasedStats, String> {
    let mut chars = value.chars();
    let Some(first) = chars.next() else {
        return Err("composition-based statistics option cannot be empty".to_string());
    };

    let mode = match first {
        '0' | 'F' | 'f' => BlastpCompositionMode::NoCompositionBasedStats,
        '1' => BlastpCompositionMode::CompositionBasedStats,
        '2' | 'D' | 'd' | 'T' | 't' => BlastpCompositionMode::CompositionMatrixAdjust,
        '3' => BlastpCompositionMode::ForceFullMatrixAdjust,
        _ => {
            return Err(format!(
                "invalid composition-based statistics mode '{value}'"
            ))
        }
    };

    // CLI v2 requires the complete mode token. NCBI blast_args.cpp:879-883:
    // if (program == eBlastp && compo_mode != eNoCompositionBasedStats &&
    //     tolower(comp_stat_string[1]) == 'u') { opt.SetUnifiedP(1); }
    // The sole suffix is lowercase u, and only on enabled modes.
    let unified_p = match (chars.next(), chars.next()) {
        (None, None) => false,
        (Some('u'), None) if mode != BlastpCompositionMode::NoCompositionBasedStats => true,
        _ => {
            return Err(format!(
                "invalid composition-based statistics token '{value}'"
            ))
        }
    };

    Ok(BlastpCompBasedStats { mode, unified_p })
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_prot_options.cpp:55-83
// ```c
// void CBlastProteinOptionsHandle::SetWordSize(int ws) {
//     m_Opts->SetWordSize(ws);
//     switch (ws) {
//     case 3: m_Opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP); break;
//     case 5: m_Opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP_FAST); break;
//     case 6: m_Opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP_WD_SZ_6); break;
//     case 7: m_Opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP_WD_SZ_7); break;
//     }
//     if (ws > 4) m_Opts->SetLookupTableType(eCompressedAaLookupTable);
// }
// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlastpLookupTableType {
    AaLookupTable,
    CompressedAaLookupTable,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_prot_options.cpp:85-148
// ```c
// SetWordSize(BLAST_WORDSIZE_PROT);
// SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP);
// SetWindowSize(BLAST_WINDOW_SIZE_PROT);
// SetMatrixName(BLAST_DEFAULT_MATRIX);
// SetGapOpeningCost(BLAST_GAP_OPEN_PROT);
// SetGapExtensionCost(BLAST_GAP_EXTN_PROT);
// SetHitlistSize(500);
// SetEvalueThreshold(BLAST_EXPECT_VALUE);
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_advprot_options.cpp:48-66
// ```c
// CBlastProteinOptionsHandle::SetGappedExtensionDefaults();
// m_Opts->SetCompositionBasedStats(eCompositionMatrixAdjust);
// m_Opts->SetSmithWatermanMode(false);
// ...
// CBlastProteinOptionsHandle::SetQueryOptionDefaults();
// SetSegFiltering(false);
// ```
// CLI v2 names/defaults: NCBI c++/src/algo/blast/blastinput/blast_args.cpp:
// 166-170,203-207,332-349,2657-2660,3158-3163; cmdline_flags.cpp:46-94.
// arg_desc.AddOptionalKey(kArgMaxHSPsPerSubject, "int_value", ..., eInteger);
// arg_desc.SetConstraint(kArgMaxHSPsPerSubject, new CArgAllowValuesGreaterThanOrEqual(1));
// arg_desc.AddDefaultKey(kArgNumThreads, "int_value", ..., NStr::IntToString(kDfltValue));
// The single-dash lexical translation is owned by crate::cli.
// NCBI reference: c++/src/algo/blast/blastinput/blastp_args.cpp:44-60
// ```c
// CBlastpAppArgs::CBlastpAppArgs()
// {
//     CRef<IBlastCmdLineArgs> arg;
//     static const string kProgram("blastp");
//     arg.Reset(new CProgramDescriptionArgs(kProgram, "Protein-Protein BLAST"));
//     const bool kQueryIsProtein = true;
//     bool const kFilterByDefault = false;
//     m_Args.push_back(arg);
//     m_ClientId = kProgram + " " + CBlastVersion().Print();
//
//     static const char kDefaultTask[] = "blastp";
//     SetTask(kDefaultTask);
//     set<string> tasks
//         (CBlastOptionsFactory::GetTasks(CBlastOptionsFactory::eProtProt));
//     arg.Reset(new CTaskCmdLineArgs(tasks, kDefaultTask));
//     m_Args.push_back(arg);
// ```
// The values are read as NCBI's argument parser reads them (`value_parsers.rs`), with
// NCBI's declared constraints only; the checks that NCBI makes after reading the
// arguments are in `check_options`, in NCBI's order. The single-dash lexical translation
// is owned by crate::cli.
#[derive(Args, Debug, Clone)]
#[command(rename_all = "snake_case")]
pub struct BlastpArgs {
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427
    // ```c
    //     arg_desc.AddDefaultKey(kArgQuery, "input_file",
    //                      "Input file name",
    //                      CArgDescriptions::eInputFile, kDfltArgQuery);
    // ```
    // `-` (the default) is standard input.
    #[arg(long, default_value = "-", value_parser = ncbi_input_path(), value_name = "PATH")]
    pub query: PathBuf,
    // Optional for the parser: a missing `-subject` is NCBI's error after the parsing
    // (`blastinput/app.rs` `missing_subject_error`).
    #[arg(long, value_parser = ncbi_input_path(), value_name = "PATH")]
    pub subject: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2191-2199
    // ```c++
    //     arg_desc.SetCurrentGroup("General search options");
    //     // database filename
    //     if (m_IsIgBlast){
    //         arg_desc.AddOptionalKey(kArgDb, "database_name", "Optional additional database name",
    //                                 CArgDescriptions::eString);
    //     } else {
    //         arg_desc.AddOptionalKey(kArgDb, "database_name", "BLAST database name",
    //                                 CArgDescriptions::eString);
    //     }
    // ```
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
    // (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES, exception 1).
    #[arg(
        long,
        value_name = "NAME",
        conflicts_with_all = ["subject", "subject_loc"],
        help = "BLAST database name (a protein database, with or without parsed Seq-ids)"
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
    // No constraint and no dependency: read with `-db` and with `-subject`.
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:234-238
    // ```c++
    //     arg_desc.AddOptionalKey(kArgEffSearchSpace, "int_value",
    //                             "Effective length of the search space",
    //                             CArgDescriptions::eInt8);
    // ```
    #[arg(long, value_parser = crate::blastinput::value_parsers::nonnegative_ncbi_int8, value_name = "INT_VALUE", help = "Effective length of the search space")]
    pub searchsp: Option<i64>,
    // NCBI reference: c++/src/algo/blast/api/blast_options_handle.cpp:225-228
    // ```c
    //     if (choice == eProtProt || choice == eAll) {
    //         retval.insert("blastp");
    //         retval.insert("blastp-short");
    //         retval.insert("blastp-fast");
    // ```
    // `blastp-short` (word size 2, PAM30 9/1) is an NCBI task that LOSAT's BLASTP does not
    // implement: it fails when the arguments are parsed, as the search entrance did
    // (`blast_engine::check_options`).
    #[arg(long, default_value = "blastp", value_parser = task_parser(&["blastp", "blastp-fast"], &[("blastp-short", BLASTP_SHORT_UNSUPPORTED)]))]
    pub task: String,
    // An omitted value keeps the default of the task (10).
    #[arg(long, value_parser = blastp_real, help = "Expectation value (E) threshold for saving hits [default: 10]")]
    pub evalue: Option<f64>,
    #[arg(long, value_parser = blastp_threshold, help = "Minimum word score such that the word is added to the BLAST lookup table [default: 11, or the matrix's suggestion; 19.3, 21 or 20.25 for word sizes 5, 6 and 7; 20 for blastp-fast]")]
    pub threshold: Option<f64>,
    #[arg(long, value_parser = protein_word_size, help = "Word size for wordfinder algorithm [default: 3; 5 for blastp-fast]")]
    pub word_size: Option<i32>,
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
    #[arg(long, value_name = "PATH", value_parser = ncbi_output_path())]
    pub out: Option<PathBuf>,
    // An omitted value keeps the default hit list size (500) and NCBI's numbers of
    // descriptions and alignments of the pairwise report (500, 250).
    #[arg(long, value_parser = blastn_count, help = "Maximum number of aligned sequences to keep [default: 500]", conflicts_with_all = ["num_descriptions", "num_alignments"])]
    pub max_target_seqs: Option<usize>,
    // NCBI blast_args.cpp:2665-2685,2732-2738: 0 or more; -max_target_seqs excludes both.
    #[arg(long = "num_descriptions", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show one-line descriptions for. Not applicable for outfmt > 4 [default: 500]")]
    pub num_descriptions: Option<usize>,
    #[arg(long = "num_alignments", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show alignments for [default: 250]")]
    pub num_alignments: Option<usize>,
    #[arg(long = "max_hsps", value_parser = blastn_count)]
    pub max_hsps_per_subject: Option<usize>,
    // NCBI blast_args.cpp:196-201,313-315: no filter when omitted or 0.
    #[arg(long = "qcov_hsp_perc", value_parser = blastp_query_coverage, help = "Percent query coverage per hsp")]
    pub qcov_hsp_perc: Option<f64>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3335-3340
    // ```c
    // CHspFilteringArgs::ExtractAlgorithmOptions(const CArgs& args,
    //                                            CBlastOptions& opts)
    // {
    //     if (args[kArgCullingLimit]) {
    //         opts.SetCullingLimit(args[kArgCullingLimit].AsInteger());
    //     }
    // ```
    /// If the query range of a hit is enveloped by that of at least this many
    /// higher-scoring hits, delete the hit (default: 0, no culling).
    #[arg(long, value_parser = nonnegative_ncbi_integer)]
    pub culling_limit: Option<i32>,
    #[arg(long)]
    pub ungapped: bool,
    #[arg(long, value_parser = nonnegative_ncbi_integer, help = "Multiple hits window size, use 0 to specify 1-hit algorithm [default: 40, or the matrix's suggestion]")]
    pub window_size: Option<i32>,
    #[arg(long, help = "Scoring matrix name [default: BLOSUM62]")]
    pub matrix: Option<String>,
    #[arg(long = "gapopen", value_parser = ncbi_integer, help = "Cost to open a gap [default: 11, or the matrix's best value]")]
    pub gap_open: Option<i32>,
    #[arg(long = "gapextend", value_parser = ncbi_integer, help = "Cost to extend a gap [default: 1, or the matrix's best value]")]
    pub gap_extend: Option<i32>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:796-797
    // ```c
    //     arg_desc.AddDefaultKey(kArgCompBasedStats, "compo", legend,
    //                            CArgDescriptions::eString, m_DefaultOpt);
    // ```
    // Read when NCBI's handler reads it (`blastinput/app.rs` `parse_comp_based_stats`).
    #[arg(
        long,
        default_value = "2",
        help = "Use composition-based statistics: 0 (F, f) none, 1 composition-based statistics, 2 (D, d, T, t) conditional compositional score matrix adjustment, 3 unconditional compositional score matrix adjustment"
    )]
    pub comp_based_stats: String,
    // Read when NCBI's filtering handler reads it (`blastinput/app.rs` `parse_seg_option`);
    // an omitted value keeps the default of the task (no filtering).
    #[arg(long, help = "SEG: no, yes, or \"WINDOW LOCUT HICUT\" [default: no]")]
    pub seg: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:339-342
    // ```c++
    //         arg_desc.AddDefaultKey(kArgLookupTableMaskingOnly, "soft_masking",
    //                         "Apply filtering locations as soft masks",
    //                         CArgDescriptions::eBoolean,
    //                         kDfltArgLookupTableMaskingOnlyProt);
    // ```
    #[arg(long, value_parser = ncbi_boolean, help = "Apply filtering locations as soft masks [default: false]")]
    pub soft_masking: Option<bool>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:803-806
    // ```c
    //     arg_desc.AddFlag(kArgUseSWTraceback,
    //                      "Compute locally optimal Smith-Waterman alignments?",
    //                      true);
    // ```
    #[arg(long)]
    pub use_sw_tback: bool,
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
    // NCBI's argument is a string (blast_args.cpp:2657-2660), parsed before the option
    // handlers (`blastinput/app.rs` `parse_formatting_string`).
    #[arg(long, default_value = "0", value_name = "SPEC")]
    pub outfmt: String,
    /// The run's BLAST archive (outfmt 11), made by `run` and `run_local` before the search
    /// reads the arguments (`archive::handle_for`).
    #[arg(skip)]
    pub archive: Option<crate::report::archive::ArchiveHandle>,
    /// The application whose search this is: psiblast's command line sets it
    /// (`algorithm/psiblast`).
    #[arg(skip)]
    pub app: BlastpApp,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1117-1119
    // ```c++
    //     if (args[kArgPSIInclusionEThreshold]) {
    //         opt.SetInclusionThreshold(args[kArgPSIInclusionEThreshold].AsDouble());
    //     }
    // ```
    /// psiblast's `-inclusion_ethresh`; `None` keeps the options' default
    /// (`PSI_INCLUSION_ETHRESH`, which blastp's options have too).
    #[arg(skip)]
    pub inclusion_ethresh: Option<f64>,
}

// NCBI reference: c++/src/algo/blast/api/blast_options_handle.cpp:381-402
// ```c
//     else if (!NStr::CompareNocase(task, "blastp") ||
//              !NStr::CompareNocase(task, "blastp-short") ||
//              !NStr::CompareNocase(task, "blastp-fast"))
//     {
//          CBlastAdvancedProteinOptionsHandle* opts =
//                dynamic_cast<CBlastAdvancedProteinOptionsHandle*>
//                 (CBlastOptionsFactory::Create(eBlastp, locality));
//          if (task == "blastp-short") {
//             opts->SetMatrixName("PAM30");
//             opts->SetGapOpeningCost(9);
//             opts->SetGapExtensionCost(1);
//             opts->SetEvalueThreshold(20000);
//             opts->SetWordSize(2);
//             opts->ClearFilterOptions();
//          } else if (task == "blastp-fast") {
//             opts->SetWordSize(5);
//             opts->SetOptions().SetLookupTableType(eCompressedAaLookupTable);
//             opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP_FAST);
//             opts->SetChaining(true);
//          }
//          retval = opts;
//     }
// ```
// NCBI reference: c++/src/algo/blast/api/blast_prot_options.cpp:56-82
// ```c
// void CBlastProteinOptionsHandle::SetWordSize(int ws) {
//
//    	m_Opts->SetWordSize(ws);
//    	switch (ws) {
//    		case 3:
//    		m_Opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP);
//    		break;
// ...
//    		default:
//    		m_Opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP);
//    		break;
//    	}
//
//    	if (ws > 4) {
//    		m_Opts->SetLookupTableType(eCompressedAaLookupTable);
//    	}
//    	else {
//    		m_Opts->SetLookupTableType(eAaLookupTable);
//    	}
// }
// ```
// The options of a task before the command line's values (BLAST_EXPECT_VALUE 10,
// BLAST_WORDSIZE_PROT 3, BLAST_WORD_THRESHOLD_BLASTP 11, BLOSUM62 11/1, window 40, SEG
// off for the advanced protein handle).
#[derive(Debug, Clone)]
struct BlastpTaskOptions {
    evalue: f64,
    matrix_name: String,
    gap_open: i32,
    gap_extend: i32,
    word_size: i32,
    threshold: f64,
    lookup_table_type: BlastpLookupTableType,
    window_size: i32,
    seg: BlastpSegSpec,
    chaining: bool,
}

impl BlastpTaskOptions {
    fn create(task: &str) -> Result<Self> {
        let mut options = Self {
            evalue: 10.0,
            matrix_name: "BLOSUM62".to_string(),
            gap_open: 11,
            gap_extend: 1,
            word_size: 3,
            threshold: 11.0,
            lookup_table_type: BlastpLookupTableType::AaLookupTable,
            window_size: 40,
            seg: BlastpSegSpec::No,
            chaining: false,
        };
        match task {
            "blastp" => {}
            "blastp-short" => {
                options.matrix_name = "PAM30".to_string();
                options.gap_open = 9;
                options.gap_extend = 1;
                options.evalue = 20000.0;
                options.word_size = 2;
                options.threshold = 11.0;
                options.lookup_table_type = BlastpLookupTableType::AaLookupTable;
            }
            "blastp-fast" => {
                options.word_size = 5;
                options.threshold = 20.0;
                options.lookup_table_type = BlastpLookupTableType::CompressedAaLookupTable;
                options.chaining = true;
            }
            _ => bail!("unsupported blastp task '{task}': expected one of blastp, blastp-fast"),
        }
        Ok(options)
    }

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
    fn set_word_size(&mut self, ws: i32) {
        self.word_size = ws;
        if self.lookup_table_type == BlastpLookupTableType::CompressedAaLookupTable && ws <= 4 {
            self.lookup_table_type = BlastpLookupTableType::AaLookupTable;
        } else if self.lookup_table_type == BlastpLookupTableType::AaLookupTable && ws > 4 {
            self.lookup_table_type = BlastpLookupTableType::CompressedAaLookupTable;
        }
    }
}

#[derive(Debug, Clone)]
pub struct ResolvedBlastpArgs {
    pub query: PathBuf,
    pub subject: PathBuf,
    /// The `-db` value as NCBI keeps it (`db_name`), when the search is of a database.
    pub db: Option<String>,
    /// `-dbsize` (`eff_len_options->db_length`; 0 when it is not given).
    pub dbsize: i64,
    /// `-searchsp`: the effective search space when above 0, and at any value no Gumbel block
    /// (blast_args.cpp:303-307, blast_stat.c:921-923).
    pub searchsp: Option<i64>,
    pub task: String,
    pub evalue: f64,
    pub threshold: f64,
    pub word_size: usize,
    pub lookup_table_type: BlastpLookupTableType,
    pub num_threads: usize,
    pub out: Option<PathBuf>,
    /// The hit list size (`ReportCounts::hitlist_size`): the -max_target_seqs value, the
    /// larger of the numbers of descriptions and alignments of the pairwise report when
    /// one of them is given, the -num_alignments value of a tabular format, or 500.
    pub max_target_seqs: usize,
    /// The numbers of descriptions and alignments of the pairwise report: the
    /// -max_target_seqs value, -num_descriptions and -num_alignments, or 500 and 250.
    pub num_descriptions: usize,
    pub num_alignments: usize,
    /// Whether -max_target_seqs was given (the text of LOSAT's limit on the hit list size).
    pub max_target_seqs_given: bool,
    pub max_hsps_per_subject: usize,
    /// The -qcov_hsp_perc value, or 0 (no filter).
    pub query_cov_hsp_perc: f64,
    /// The -culling_limit value (0: no culling, as NCBI's `SetCullingLimit(0)`).
    pub culling_limit: i32,
    pub ungapped: bool,
    pub window_size: usize,
    /// The matrix as the command line names it (NCBI keeps the spelling); `scoring.matrix`
    /// is its LOSAT matrix when LOSAT has it, otherwise BLOSUM62 (and the search rejects
    /// the name, `blast_engine.rs` `validate_requested_blastp_support`).
    pub matrix_name: String,
    pub scoring: ProteinScoringSpec,
    pub comp_based_stats: BlastpCompBasedStats,
    pub seg: BlastpSegSpec,
    /// `-soft_masking` (`mask_at_hash`): the SEG intervals keep words out of the lookup table
    /// but the searched query keeps its letters (`encoding.rs` `encode_query_frame`).
    pub soft_masking: bool,
    pub use_sw_tback: bool,
    pub chaining: bool,
    /// The `-query_loc` range, as given.
    pub query_range: Option<crate::blastinput::seq_range::SequenceRange>,
    pub outfmt: String,
    /// The run's BLAST archive (outfmt 11), shared by the blocks and batches of the run.
    pub archive: Option<crate::report::archive::ArchiveHandle>,
    /// The application whose search this is (`BlastpArgs::app`).
    pub app: BlastpApp,
    /// The inclusion threshold (`psi_options->inclusion_ethresh`): psiblast's
    /// `-inclusion_ethresh`, or `PSI_INCLUSION_ETHRESH`.
    pub inclusion_ethresh: f64,
    /// `traceback_hitlist_size`, computed once (psiblast's reads `ADAPTIVE_CBS`).
    traceback_hitlist: usize,
    /// Whether `taxdb.bti` is on NCBI's search path for a `-subject` search whose custom tabular
    /// report checks it (`x_CheckTaxDB`, `blast_engine.rs` `search_cli`); false otherwise, and
    /// for `run_local`, which has no search path (C2 DECISIONS 7). A `-db` search asks its
    /// database (`BlastpDatabase::taxdb_found`).
    pub subject_taxdb_found: bool,
}

impl ResolvedBlastpArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:213-224
    // ```c++
    //     // When dealing with PSI-BLAST iterations, we need to keep a larger
    //     // alignment for the PSSM engine as to replicate blastpgp's behavior
    //     int hitlist_size_backup = m_OptsMemento->m_HitSaveOpts->hitlist_size;
    //     if (m_OptsMemento->m_ProgramType == eBlastTypePsiBlast ) {
    //         SBlastHitsParameters* bhp = NULL;
    //         SBlastHitsParametersNew(m_OptsMemento->m_HitSaveOpts,
    //                                 m_OptsMemento->m_ExtnOpts,
    //                                 m_OptsMemento->m_ScoringOpts,
    //                                 &bhp);
    //         m_OptsMemento->m_HitSaveOpts->hitlist_size = bhp->prelim_hitlist_size;
    //         SBlastHitsParametersFree(bhp);
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:92-93
    // ```c
    //        (*retval)->prelim_hitlist_size = GetPrelimHitlistSize(hit_options->hitlist_size,
    //     		   ext_options->compositionBasedStats, scoring_options->gapped_calculation);
    // ```
    /// The hit list size of the traceback (`hitParams->options->hitlist_size`: the heaps and the
    /// hit lists of the composition-based redo): the hit list size, or psiblast's preliminary hit
    /// list size. psiblast's formatter then keeps the first `max_target_seqs` subjects
    /// (`PruneSeqalign` with `m_HitlistSize`, blast_format.cpp:813). The preliminary size is
    /// positive (`validate_requested_blastp_support`).
    pub fn traceback_hitlist_size(&self) -> usize {
        self.traceback_hitlist
    }

    fn compute_traceback_hitlist_size(&self) -> usize {
        match self.app {
            BlastpApp::Blastp => self.max_target_seqs,
            BlastpApp::Psiblast => crate::algorithm::blastn::hsp::get_prelim_hitlist_size(
                self.max_target_seqs,
                self.comp_based_stats.is_enabled(),
                true,
            )
            .max(1) as usize,
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:465-468
    // ```c++
    //             CRef<CBlastFastaInputSource> fasta(new CBlastFastaInputSource(
    //                                          m_CmdLineArgs->GetInputStream(),
    //                                          iconfig));
    //             input.Reset(new CBlastInput(&*fasta, 1));
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:215-216
    // ```c++
    //         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
    //         CBlastInput input(&fasta, m_CmdLineArgs->GetQueryBatchSize());
    // ```
    /// The query batch size: blastp's (`BATCH_SIZE` or 10000 letters), or psiblast's 1, which
    /// reads no variable (a batch of one query, and of the queries without letters before it).
    pub fn query_batch_size(&self) -> Result<u32> {
        match self.app {
            BlastpApp::Blastp => crate::blastinput::app::query_batch_size("BLASTP", 10000),
            BlastpApp::Psiblast => Ok(1),
        }
    }
}

impl BlastpArgs {
    /// `-culling_limit`, 0 (no culling) when it is not given. The BLAST archive records
    /// whether it was given (`culling_limit`).
    pub fn culling_limit_value(&self) -> i32 {
        self.culling_limit.unwrap_or(0)
    }

    /// The options as `check_options` resolves them, without NCBI's warnings.
    pub fn resolve(&self) -> Result<ResolvedBlastpArgs> {
        self.check_options(&mut std::io::sink())
    }

    /// NCBI's `CBlastAppArgs::SetOptions` for blastp after the files are opened: the task's
    /// options, each argument handler in the order of `CBlastpAppArgs` (with the warnings
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
    pub fn check_options(
        &self,
        diagnostics: &mut dyn std::io::Write,
    ) -> Result<ResolvedBlastpArgs> {
        use crate::blastinput::app::{options_error, CompositionMode};
        use crate::stats::protein_options::{
            protein_gap_existence_extend_params, suggested_threshold, suggested_window_size,
            ProteinOptionsCheck, SuggestionProgram,
        };
        let mut options = BlastpTaskOptions::create(&self.task)?;
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:251-301
        // ```c
        //     if (args.Exist(kArgEvalue) && args[kArgEvalue]) {
        //         opt.SetEvalueThreshold(args[kArgEvalue].AsDouble());
        //     }
        //
        //     int gap_open=0, gap_extend=0;
        //     if (args.Exist(kArgMatrixName) && args[kArgMatrixName])
        //          BLAST_GetProteinGapExistenceExtendParams
        //              (args[kArgMatrixName].AsString().c_str(), &gap_open, &gap_extend);
        //
        //     if (args.Exist(kArgGapOpen) && args[kArgGapOpen]) {
        //         opt.SetGapOpeningCost(args[kArgGapOpen].AsInteger());
        //     }
        //     else if (args.Exist(kArgMatrixName) && args[kArgMatrixName]) {
        //         opt.SetGapOpeningCost(gap_open);
        //     }
        // ...
        //     if ( args.Exist(kArgWordSize) && args[kArgWordSize]) {
        //         if (m_QueryIsProtein && args[kArgWordSize].AsInteger() > 4){
        //            opt.SetLookupTableType(eCompressedAaLookupTable);
        //            opt.SetWordThreshold(19.3);
        //            if (args[kArgWordSize].AsInteger() > 5) {
        //                opt.SetWordThreshold(21.0);
        //            }
        //            if (args[kArgWordSize].AsInteger() > 6) {
        //                opt.SetWordThreshold(20.25);
        //            }
        //         }
        //         opt.SetWordSize(args[kArgWordSize].AsInteger());
        //
        //     }
        // ```
        if let Some(evalue) = self.evalue {
            options.evalue = evalue;
        }
        let (matrix_gap_open, matrix_gap_extend) = self
            .matrix
            .as_deref()
            .and_then(protein_gap_existence_extend_params)
            .unwrap_or((0, 0));
        match (self.gap_open, &self.matrix) {
            (Some(gap_open), _) => options.gap_open = gap_open,
            (None, Some(_)) => options.gap_open = matrix_gap_open,
            (None, None) => {}
        }
        match (self.gap_extend, &self.matrix) {
            (Some(gap_extend), _) => options.gap_extend = gap_extend,
            (None, Some(_)) => options.gap_extend = matrix_gap_extend,
            (None, None) => {}
        }
        if let Some(word_size) = self.word_size {
            if word_size > 4 {
                options.lookup_table_type = BlastpLookupTableType::CompressedAaLookupTable;
                options.threshold = 19.3;
                if word_size > 5 {
                    options.threshold = 21.0;
                }
                if word_size > 6 {
                    options.threshold = 20.25;
                }
            }
            options.set_word_size(word_size);
        }
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:396-406 (the
        // filtering handler; `blastinput/app.rs` `parse_seg_option`)
        // ```c
        //         if (m_QueryIsProtein && args[kArgSegFiltering]) {
        //             const string& seg_opts = args[kArgSegFiltering].AsString();
        // ```
        if let Some(seg) = &self.seg {
            options.seg = crate::blastinput::app::parse_seg_option(seg, self.app.name())?;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:387-391
        // ```c++
        // CFilteringArgs::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& opt)
        // {
        //     if (args[kArgLookupTableMaskingOnly]) {
        //         opt.SetMaskAtHash(args[kArgLookupTableMaskingOnly].AsBoolean());
        //     }
        // ```
        // The argument always has its default (`kDfltArgLookupTableMaskingOnlyProt`, "false"),
        // and no other handler of blastp sets `mask_at_hash`.
        let soft_masking = self.soft_masking.unwrap_or(false);
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:634-639
        // ```c
        // CMatrixNameArg::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& opt)
        // {
        //     if (args[kArgMatrixName]) {
        //         opt.SetMatrixName(args[kArgMatrixName].AsString().c_str());
        //     }
        // }
        // ```
        if let Some(matrix) = &self.matrix {
            options.matrix_name = matrix.clone();
        }
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:586-623
        // ```c
        // s_IsDefaultWordThreshold(EProgram program, double threshold)
        // {
        //     int word_threshold = static_cast<int>(threshold);
        //     bool retval = true;
        //     if (program == eBlastp &&
        //         word_threshold != BLAST_WORD_THRESHOLD_BLASTP) {
        //         retval = false;
        // ...
        //     if (args[kArgWordScoreThreshold]) {
        //         opt.SetWordThreshold(args[kArgWordScoreThreshold].AsDouble());
        //     } else if (s_IsDefaultWordThreshold(opt.GetProgram(),
        //                                         opt.GetWordThreshold())) {
        //         double threshold = -1;
        //         BLAST_GetSuggestedThreshold(opt.GetProgramType(),
        //                                     opt.GetMatrixName(),
        //                                     &threshold);
        // ```
        if let Some(threshold) = self.threshold {
            options.threshold = threshold;
        } else if options.threshold as i32 == 11 {
            options.threshold =
                suggested_threshold(SuggestionProgram::Protein, &options.matrix_name);
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
        options.window_size = match self.window_size {
            Some(window) => window,
            None => suggested_window_size(&options.matrix_name),
        };
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1995-1999
        // ```c++
        //     // set the sequence range
        //     if (args.Exist(kArgQueryLocation) && args[kArgQueryLocation]) {
        //         m_Range = ParseSequenceRange(args[kArgQueryLocation].AsString(),
        //                                      "Invalid specification of query location");
        //     }
        // ```
        // The query options handler comes after the window size handler and before the
        // formatting handler (blastp_args.cpp:44-120).
        let query_range = crate::blastinput::seq_range::parse_optional_range(
            self.query_loc.as_deref(),
            crate::blastinput::seq_range::RangeRole::Query,
            self.app.name(),
        )?;
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2874-2886 and
        // 2960-2977 (the formatting handler: other programs' formats, the hit list size
        // and its warning)
        // ```c
        //     if(hitlist_size < 5){
        //    		ERR_POST(Warning << "Examining 5 or more matches is recommended");
        //     }
        // ```
        let choice = crate::blastinput::app::parse_formatting_string(&self.outfmt)?;
        crate::blastinput::app::formatting_handler_check(
            &choice,
            crate::blastinput::app::Program::Blastp,
        )?;
        // The hit list size and the counts of the pairwise report
        // (`blastinput::report_counts`), with the handler's warnings.
        let (counts, warnings) = crate::blastinput::report_counts::run_report_counts(
            self.app.ncbi_name(),
            &[choice.number],
            self.max_target_seqs,
            self.num_descriptions,
            self.num_alignments,
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
        // The unified P-values are blastp's only (`program == eBlastp`, blast_args.cpp:880-883).
        let (mode, unified_p) = crate::blastinput::app::parse_comp_based_stats(
            &self.comp_based_stats,
            self.app == BlastpApp::Blastp,
            self.ungapped,
        )?;
        // NCBI's `Validate` (`BLAST_ValidateOptions`), shared with TBLASTN.
        crate::stats::protein_options::validate_protein_options(&ProteinOptionsCheck {
            gapped: !self.ungapped,
            gap_open: options.gap_open,
            gap_extend: options.gap_extend,
            matrix_name: &options.matrix_name,
            threshold: options.threshold,
            word_size: options.word_size,
            compressed_lookup: options.lookup_table_type
                == BlastpLookupTableType::CompressedAaLookupTable,
            hitlist_size: counts.hitlist_size,
            evalue: options.evalue,
            psiblast: self.app == BlastpApp::Psiblast,
        })?;
        let matrix = options
            .matrix_name
            .parse::<ScoringMatrix>()
            .unwrap_or(ScoringMatrix::Blosum62);
        let mut resolved = ResolvedBlastpArgs {
            query: self.query.clone(),
            subject: self.subject.clone().unwrap_or_default(),
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
            // ```c++
            //     if (args.Exist(kArgDb) && args[kArgDb]) {
            // 	std::string local_dblist = NStr::TruncateSpaces( args[kArgDb].AsString() );
            //
            //         m_SearchDb.Reset(new CSearchDatabase( local_dblist,
            //                                              mol_type));
            // ```
            // The name without the white space at its ends: what the reader resolves and
            // what the reports show.
            db: self
                .db
                .as_deref()
                .map(crate::seqdb::db_argument_name)
                .map(str::to_string),
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2570-2572
            // ```c++
            //     if (args.Exist(kArgDbSize) && args[kArgDbSize]) {
            //         opts.SetDbLength(args[kArgDbSize].AsInt8());
            //     }
            // ```
            // The options start with a database length of 0
            // (blast_prot_options.cpp:150-156, `SetDbLength(0)`).
            dbsize: self.dbsize.unwrap_or(0),
            searchsp: self.searchsp,
            task: self.task.clone(),
            evalue: options.evalue,
            threshold: options.threshold,
            word_size: options.word_size as usize,
            lookup_table_type: options.lookup_table_type,
            num_threads: self.num_threads,
            out: self.out.clone(),
            max_target_seqs: counts.hitlist_size,
            num_descriptions: counts.num_descriptions,
            num_alignments: counts.num_alignments,
            max_target_seqs_given: self.max_target_seqs.is_some(),
            // NCBI blast_args.cpp:317-318: if (args[kArgMaxHSPsPerSubject])
            // opt.SetMaxHspsPerSubject(args[kArgMaxHSPsPerSubject].AsInteger());
            max_hsps_per_subject: self.max_hsps_per_subject.unwrap_or(0),
            // NCBI blast_args.cpp:313-315: if (args[kArgQueryCovHspPerc])
            // opt.SetQueryCovHspPerc(args[kArgQueryCovHspPerc].AsDouble());
            query_cov_hsp_perc: self.qcov_hsp_perc.unwrap_or(0.0),
            culling_limit: self.culling_limit_value(),
            ungapped: self.ungapped,
            window_size: usize::try_from(options.window_size).unwrap_or(0),
            matrix_name: options.matrix_name,
            scoring: ProteinScoringSpec {
                matrix,
                gap_open: options.gap_open,
                gap_extend: options.gap_extend,
            },
            comp_based_stats: BlastpCompBasedStats {
                mode: match mode {
                    CompositionMode::NoCompositionBasedStats => {
                        BlastpCompositionMode::NoCompositionBasedStats
                    }
                    CompositionMode::CompositionBasedStats => {
                        BlastpCompositionMode::CompositionBasedStats
                    }
                    CompositionMode::CompositionMatrixAdjust => {
                        BlastpCompositionMode::CompositionMatrixAdjust
                    }
                    CompositionMode::CompoForceFullMatrixAdjust => {
                        BlastpCompositionMode::ForceFullMatrixAdjust
                    }
                },
                unified_p,
            },
            seg: options.seg,
            soft_masking,
            use_sw_tback: self.use_sw_tback,
            chaining: options.chaining,
            query_range,
            outfmt: self.outfmt.clone(),
            archive: self.archive.clone(),
            app: self.app,
            inclusion_ethresh: self
                .inclusion_ethresh
                .unwrap_or(super::blast_engine::PSI_INCLUSION_ETHRESH),
            traceback_hitlist: 0,
            subject_taxdb_found: false,
        };
        resolved.traceback_hitlist = resolved.compute_traceback_hitlist_size();
        Ok(resolved)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3152-3171 (`s_TruncateSpaces`
    // with `eTrunc_Both`: `isspace` at both ends; `seqdb::db_argument_name`).
    #[test]
    fn the_database_name_loses_the_white_space_at_its_ends() {
        let resolved_name = |given: &str| {
            let cli = parse_cli(["losat", "blastp", "-query", "q.faa", "-db", given]);
            let TestCommand::Blastp(args) = cli.command;
            args.resolve().unwrap().db
        };
        for (given, kept) in [
            ("db", "db"),
            (" db ", "db"),
            ("\tdir/db\n", "dir/db"),
            ("\r\x0b\x0c db", "db"),
            ("a b", "a b"),
            (" a  b ", "a  b"),
            ("", ""),
            ("   ", ""),
            ("\u{a0}db", "\u{a0}db"),
        ] {
            assert_eq!(resolved_name(given).as_deref(), Some(kept), "{given:?}");
        }
    }

    // NCBI blast_args.cpp:332-349: parse the same single string as the public CLI.
    fn parse_cli<'a>(args: impl IntoIterator<Item = &'a str>) -> TestCli {
        crate::cli::try_parse_from(args).unwrap()
    }

    #[derive(clap::Parser, Debug)]
    struct TestCli {
        #[command(subcommand)]
        command: TestCommand,
    }

    #[derive(clap::Subcommand, Debug)]
    enum TestCommand {
        Blastp(BlastpArgs),
    }

    #[test]
    fn test_default_blastp_options_resolve_to_ncbi_defaults() {
        let cli = parse_cli(["losat", "blastp", "-query", "q.faa", "-subject", "s.faa"]);
        let TestCommand::Blastp(args) = cli.command;
        assert_eq!(args.evalue, None);
        assert_eq!(args.threshold, None);
        assert_eq!(args.word_size, None);
        assert_eq!(args.window_size, None);
        assert_eq!(args.matrix, None);
        assert_eq!(args.gap_open, None);
        assert_eq!(args.gap_extend, None);
        assert_eq!(args.comp_based_stats, "2");
        assert_eq!(args.seg, None);
        assert!(!args.ungapped);
        assert!(!args.use_sw_tback);

        let resolved = args.resolve().expect("resolved blastp args");
        assert_eq!(resolved.task, "blastp");
        assert_eq!(resolved.evalue, 10.0);
        assert_eq!(resolved.word_size, 3);
        assert_eq!(
            resolved.lookup_table_type,
            BlastpLookupTableType::AaLookupTable
        );
        assert_eq!(resolved.threshold, 11.0);
        assert_eq!(resolved.window_size, 40);
        assert_eq!(resolved.scoring.matrix, ScoringMatrix::Blosum62);
        assert_eq!(resolved.scoring.gap_open, 11);
        assert_eq!(resolved.scoring.gap_extend, 1);
        assert_eq!(
            resolved.comp_based_stats,
            BlastpCompBasedStats {
                mode: BlastpCompositionMode::CompositionMatrixAdjust,
                unified_p: false,
            }
        );
        assert_eq!(resolved.seg, BlastpSegSpec::No);
        assert_eq!(resolved.max_target_seqs, 500);
        assert_eq!(resolved.max_hsps_per_subject, 0);
        assert!(!resolved.ungapped);
        assert!(!resolved.use_sw_tback);
        assert!(!resolved.chaining);
        assert!(!resolved.soft_masking);
        assert_eq!(resolved.outfmt, "0");
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:387-391
    // ```c++
    //     if (args[kArgLookupTableMaskingOnly]) {
    //         opt.SetMaskAtHash(args[kArgLookupTableMaskingOnly].AsBoolean());
    //     }
    // ```
    // `-soft_masking` takes NCBI's boolean spellings (`NStr::StringToBool`); a task does not
    // change it (`blastp-short`'s `ClearFilterOptions` comes before the handler).
    #[test]
    fn soft_masking_sets_mask_at_hash() {
        let resolve = |extra: &[&str]| {
            let mut argv = vec!["losat", "blastp", "-query", "q.faa", "-subject", "s.faa"];
            argv.extend_from_slice(extra);
            let TestCommand::Blastp(args) = parse_cli(argv).command;
            args.resolve().expect("resolved blastp args").soft_masking
        };
        for value in ["true", "T", "1", "YES", "on"] {
            assert!(resolve(&["-soft_masking", value]), "{value}");
        }
        for value in ["false", "F", "0", "no", "off"] {
            assert!(!resolve(&["-soft_masking", value]), "{value}");
        }
        // `-task blastp-short` stops at the command line, so the task is set by hand.
        let TestCommand::Blastp(mut args) = parse_cli([
            "losat",
            "blastp",
            "-query",
            "q.faa",
            "-subject",
            "s.faa",
            "-soft_masking",
            "true",
        ])
        .command;
        args.task = "blastp-short".into();
        assert!(args.resolve().expect("resolved blastp args").soft_masking);
        assert!(resolve(&[
            "-task",
            "blastp-fast",
            "-soft_masking",
            "true",
            "-seg",
            "yes"
        ]));
    }

    // NCBI blast_args.cpp:2894-2973 (`blastinput::report_counts`): the hit list size and the
    // counts of outfmt 0 from -num_descriptions / -num_alignments; "No hits are being saved".
    #[test]
    fn num_descriptions_and_num_alignments_set_the_counts() {
        let resolve = |extra: &[&str]| {
            let mut argv = vec!["losat", "blastp", "-query", "q.faa", "-subject", "s.faa"];
            argv.extend_from_slice(extra);
            let cli = parse_cli(argv);
            let TestCommand::Blastp(args) = cli.command;
            args.resolve()
        };
        let counts = |extra: &[&str]| {
            let resolved = resolve(extra).unwrap();
            (
                resolved.max_target_seqs,
                resolved.num_descriptions,
                resolved.num_alignments,
            )
        };
        assert_eq!(counts(&[]), (500, 500, 250));
        assert_eq!(counts(&["-num_descriptions", "3"]), (250, 3, 250));
        assert_eq!(counts(&["-num_alignments", "2"]), (500, 500, 2));
        assert_eq!(
            counts(&["-num_descriptions", "2", "-num_alignments", "4"]),
            (4, 2, 4)
        );
        assert_eq!(counts(&["-max_target_seqs", "7"]), (7, 7, 7));
        assert_eq!(counts(&["-outfmt", "6", "-num_alignments", "2"]), (2, 2, 2));
        assert_eq!(
            counts(&["-outfmt", "6", "-num_descriptions", "2"]),
            (500, 500, 500)
        );
        let error = resolve(&["-num_descriptions", "0", "-num_alignments", "0"]).unwrap_err();
        assert_eq!(
            error.to_string(),
            "BLAST query/options error: No hits are being saved\nPlease refer to the BLAST+ user manual.\n"
        );
        assert!(resolve(&["-outfmt", "7", "-num_alignments", "0"]).is_err());
        assert!(resolve(&["-num_descriptions", "0"]).is_ok());
    }

    #[test]
    fn test_blastp_short_task_resolves_ncbi_defaults() {
        // NCBI api/blast_options_handle.cpp:395-399: blastp-fast sets word size 5, the
        // compressed lookup and BLAST_WORD_THRESHOLD_BLASTP_FAST (20), which is kept
        // because (int)20 is not BLASTP's default threshold (blast_args.cpp:586-600).
        let cli = parse_cli(["losat", "blastp", "-query", "q.faa", "-subject", "s.faa"]);
        let TestCommand::Blastp(mut args) = cli.command;
        args.task = "blastp-short".into();
        let resolved = args.resolve().expect("resolved blastp-short args");
        assert_eq!(resolved.task, "blastp-short");
        assert_eq!(resolved.evalue, 20000.0);
        assert_eq!(resolved.scoring.matrix, ScoringMatrix::Pam30);
        assert_eq!(resolved.scoring.gap_open, 9);
        assert_eq!(resolved.scoring.gap_extend, 1);
        assert_eq!(resolved.word_size, 2);
        assert_eq!(resolved.threshold, 16.0);
        assert_eq!(resolved.window_size, 15);
        assert_eq!(resolved.seg, BlastpSegSpec::No);
        assert!(!resolved.chaining);
        args.evalue = Some(42.0);
        assert_eq!(args.resolve().unwrap().evalue, 42.0);
    }

    #[test]
    fn test_blastp_fast_task_resolves_ncbi_defaults() {
        // NCBI api/blast_options_handle.cpp:395-399: blastp-fast sets word size 5, the
        // compressed lookup and BLAST_WORD_THRESHOLD_BLASTP_FAST (20), which is kept
        // because (int)20 is not BLASTP's default threshold (blast_args.cpp:586-600).
        let cli = parse_cli(["losat", "blastp", "-query", "q.faa", "-subject", "s.faa"]);
        let TestCommand::Blastp(mut args) = cli.command;
        args.task = "blastp-fast".into();
        let resolved = args.resolve().expect("resolved blastp-fast args");
        assert_eq!(resolved.task, "blastp-fast");
        assert_eq!(resolved.evalue, 10.0);
        assert_eq!(resolved.word_size, 5);
        assert_eq!(resolved.threshold, 20.0);
        assert_eq!(
            resolved.lookup_table_type,
            BlastpLookupTableType::CompressedAaLookupTable
        );
        assert!(resolved.chaining);
        args.evalue = Some(42.0);
        assert_eq!(args.resolve().unwrap().evalue, 42.0);
    }

    #[test]
    fn test_seg_yes_no_and_custom_parser() {
        let cli = parse_cli([
            "losat", "blastp", "-query", "q.faa", "-subject", "s.faa", "-seg", "no",
        ]);
        let TestCommand::Blastp(args) = cli.command;
        assert_eq!(args.resolve().unwrap().seg, BlastpSegSpec::No);

        let cli = parse_cli([
            "losat", "blastp", "-query", "q.faa", "-subject", "s.faa", "-seg", "yes",
        ]);
        let TestCommand::Blastp(args) = cli.command;
        assert_eq!(args.resolve().unwrap().seg, BlastpSegSpec::Yes);

        let cli = parse_cli([
            "losat",
            "blastp",
            "-query",
            "q.faa",
            "-subject",
            "s.faa",
            "-seg",
            "15 2.5 3.0",
        ]);
        let TestCommand::Blastp(args) = cli.command;
        assert_eq!(
            args.resolve().unwrap().seg,
            BlastpSegSpec::WindowLocutHicut {
                window: 15,
                locut: 2.5,
                hicut: 3.0,
            }
        );
    }

    #[test]
    fn test_comp_based_stats_parser_supports_ncbi_modes() {
        assert_eq!(
            parse_comp_based_stats("0").unwrap(),
            BlastpCompBasedStats {
                mode: BlastpCompositionMode::NoCompositionBasedStats,
                unified_p: false,
            }
        );
        assert_eq!(
            parse_comp_based_stats("1").unwrap(),
            BlastpCompBasedStats {
                mode: BlastpCompositionMode::CompositionBasedStats,
                unified_p: false,
            }
        );
        assert_eq!(
            parse_comp_based_stats("2u").unwrap(),
            BlastpCompBasedStats {
                mode: BlastpCompositionMode::CompositionMatrixAdjust,
                unified_p: true,
            }
        );
        assert_eq!(
            parse_comp_based_stats("3").unwrap(),
            BlastpCompBasedStats {
                mode: BlastpCompositionMode::ForceFullMatrixAdjust,
                unified_p: false,
            }
        );
        assert_eq!(
            parse_comp_based_stats("T").unwrap(),
            BlastpCompBasedStats {
                mode: BlastpCompositionMode::CompositionMatrixAdjust,
                unified_p: false,
            }
        );
    }

    #[test]
    fn test_matrix_resolution_updates_gap_threshold_and_window_defaults() {
        let cli = parse_cli([
            "losat", "blastp", "-query", "q.faa", "-subject", "s.faa", "-matrix", "BLOSUM45",
        ]);
        let TestCommand::Blastp(args) = cli.command;
        let resolved = args.resolve().expect("resolved blastp args");
        assert_eq!(resolved.scoring.matrix, ScoringMatrix::Blosum45);
        assert_eq!(resolved.scoring.gap_open, 14);
        assert_eq!(resolved.scoring.gap_extend, 2);
        assert_eq!(resolved.threshold, 14.0);
        assert_eq!(resolved.window_size, 60);
    }

    #[test]
    fn test_explicit_gap_values_override_matrix_defaults_independently() {
        let cli = parse_cli([
            "losat", "blastp", "-query", "q.faa", "-subject", "s.faa", "-matrix", "PAM70",
            "-gapopen", "11",
        ]);
        let TestCommand::Blastp(args) = cli.command;
        let resolved = args.resolve().expect("resolved blastp args");
        assert_eq!(resolved.scoring.matrix, ScoringMatrix::Pam70);
        assert_eq!(resolved.scoring.gap_open, 11);
        assert_eq!(resolved.scoring.gap_extend, 1);
    }

    #[test]
    fn test_word_size_above_four_selects_compressed_lookup_thresholds() {
        let cli = parse_cli([
            "losat",
            "blastp",
            "-query",
            "q.faa",
            "-subject",
            "s.faa",
            "-word_size",
            "5",
        ]);
        let TestCommand::Blastp(args) = cli.command;
        let resolved = args.resolve().expect("resolved blastp args");
        assert_eq!(
            resolved.lookup_table_type,
            BlastpLookupTableType::CompressedAaLookupTable
        );
        assert_eq!(resolved.threshold, 19.3);

        let cli = parse_cli([
            "losat",
            "blastp",
            "-query",
            "q.faa",
            "-subject",
            "s.faa",
            "-word_size",
            "7",
        ]);
        let TestCommand::Blastp(args) = cli.command;
        let resolved = args.resolve().expect("resolved blastp args");
        assert_eq!(resolved.threshold, 20.25);
    }

    #[test]
    fn test_ungapped_with_composition_based_stats_is_rejected() {
        let cli = parse_cli([
            "losat",
            "blastp",
            "-query",
            "q.faa",
            "-subject",
            "s.faa",
            "-ungapped",
            "-comp_based_stats",
            "2",
        ]);
        let TestCommand::Blastp(args) = cli.command;
        let err = args
            .resolve()
            .expect_err("ungapped comp-based stats should fail");
        assert!(err
            .to_string()
            .contains("Composition-adjusted searched are not supported with an ungapped search"));
    }

    #[test]
    fn test_comp_based_stats_cli_string_matches_ncbi_modes() {
        assert_eq!(
            BlastpCompBasedStats {
                mode: BlastpCompositionMode::NoCompositionBasedStats,
                unified_p: false,
            }
            .to_ncbi_cli_string(),
            "0"
        );
        assert_eq!(
            BlastpCompBasedStats {
                mode: BlastpCompositionMode::CompositionBasedStats,
                unified_p: false,
            }
            .to_ncbi_cli_string(),
            "1"
        );
        assert_eq!(
            BlastpCompBasedStats {
                mode: BlastpCompositionMode::CompositionMatrixAdjust,
                unified_p: true,
            }
            .to_ncbi_cli_string(),
            "2u"
        );
        assert_eq!(
            BlastpCompBasedStats {
                mode: BlastpCompositionMode::ForceFullMatrixAdjust,
                unified_p: false,
            }
            .to_ncbi_cli_string(),
            "3"
        );
    }

    #[test]
    fn test_seg_cli_string_matches_ncbi_forms() {
        assert_eq!(BlastpSegSpec::No.to_ncbi_cli_string(), "no");
        assert_eq!(BlastpSegSpec::Yes.to_ncbi_cli_string(), "yes");
        assert_eq!(
            BlastpSegSpec::WindowLocutHicut {
                window: 12,
                locut: 2.2,
                hicut: 2.5,
            }
            .to_ncbi_cli_string(),
            "12 2.2 2.5"
        );
    }
}
