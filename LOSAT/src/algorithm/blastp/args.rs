//! Command-line arguments for BLASTP

use anyhow::{bail, Result};
use clap::Args;
use std::path::PathBuf;

use crate::config::{ProteinScoringSpec, ScoringMatrix};
use crate::utils::seg::SegParams;

pub use crate::blastinput::value_parsers::SegSpec as BlastpSegSpec;
use crate::blastinput::value_parsers::*;

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

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:845-886
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
    // NCBI reference: c++/src/algo/blast/api/blast_options_handle.cpp:225-228
    // ```c
    //     if (choice == eProtProt || choice == eAll) {
    //         retval.insert("blastp");
    //         retval.insert("blastp-short");
    //         retval.insert("blastp-fast");
    // ```
    #[arg(long, default_value = "blastp", value_parser = ["blastp", "blastp-fast", "blastp-short"])]
    pub task: String,
    // An omitted value keeps the default of the task (10; 20000 for blastp-short).
    #[arg(long, value_parser = blastp_real, help = "Expectation value (E) threshold for saving hits [default: 10; 20000 for blastp-short]")]
    pub evalue: Option<f64>,
    #[arg(long, value_parser = blastp_threshold, help = "Minimum word score such that the word is added to the BLAST lookup table [default: 11, or the matrix's suggestion; 19.3, 21 or 20.25 for word sizes 5, 6 and 7; 20 for blastp-fast]")]
    pub threshold: Option<f64>,
    #[arg(long, value_parser = protein_word_size, help = "Word size for wordfinder algorithm [default: 3; 5 for blastp-fast, 2 for blastp-short]")]
    pub word_size: Option<i32>,
    #[arg(long, default_value_t = 1, value_parser = blastn_count)]
    pub num_threads: usize,
    #[arg(long, value_name = "PATH", value_parser = ncbi_output_path())]
    pub out: Option<PathBuf>,
    // An omitted value keeps the default hit list size (500) and NCBI's numbers of
    // descriptions and alignments of the pairwise report (500, 250).
    #[arg(long, value_parser = blastn_count, help = "Maximum number of aligned sequences to keep [default: 500]")]
    pub max_target_seqs: Option<usize>,
    #[arg(long = "max_hsps", value_parser = blastn_count)]
    pub max_hsps_per_subject: Option<usize>,
    #[arg(long)]
    pub ungapped: bool,
    #[arg(long, value_parser = nonnegative_ncbi_integer, help = "Multiple hits window size, use 0 to specify 1-hit algorithm [default: 40, or the matrix's suggestion]")]
    pub window_size: Option<i32>,
    #[arg(
        long,
        help = "Scoring matrix name [default: BLOSUM62; PAM30 for blastp-short]"
    )]
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
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:803-806
    // ```c
    //     arg_desc.AddFlag(kArgUseSWTraceback,
    //                      "Compute locally optimal Smith-Waterman alignments?",
    //                      true);
    // ```
    #[arg(long)]
    pub use_sw_tback: bool,
    // NCBI's argument is a string (blast_args.cpp:2657-2660), parsed before the option
    // handlers (`blastinput/app.rs` `parse_formatting_string`).
    #[arg(long, default_value = "0", value_name = "SPEC")]
    pub outfmt: String,
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
            _ => bail!("unsupported blastp task '{task}': expected one of blastp, blastp-short, blastp-fast"),
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
    pub task: String,
    pub evalue: f64,
    pub threshold: f64,
    pub word_size: usize,
    pub lookup_table_type: BlastpLookupTableType,
    pub num_threads: usize,
    pub out: Option<PathBuf>,
    /// The hit list size: the -max_target_seqs value, or 500.
    pub max_target_seqs: usize,
    /// The -max_target_seqs value, which also sets the numbers of descriptions and
    /// alignments of the pairwise report.
    pub max_target_seqs_given: Option<usize>,
    pub max_hsps_per_subject: usize,
    pub ungapped: bool,
    pub window_size: usize,
    /// The matrix as the command line names it (NCBI keeps the spelling); `scoring.matrix`
    /// is its LOSAT matrix when LOSAT has it, otherwise BLOSUM62 (and the search rejects
    /// the name, `blast_engine.rs` `validate_requested_blastp_support`).
    pub matrix_name: String,
    pub scoring: ProteinScoringSpec,
    pub comp_based_stats: BlastpCompBasedStats,
    pub seg: BlastpSegSpec,
    pub use_sw_tback: bool,
    pub chaining: bool,
    pub outfmt: String,
}

impl BlastpArgs {
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
            options.seg = crate::blastinput::app::parse_seg_option(seg, "BLASTP")?;
        }
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
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2874-2886 and
        // 2960-2977 (the formatting handler: other programs' formats, the hit list size
        // and its warning)
        // ```c
        //     if(hitlist_size < 5){
        //    		ERR_POST(Warning << "Examining 5 or more matches is recommended");
        //     }
        // ```
        let choice = crate::blastinput::app::parse_formatting_string(&self.outfmt)?;
        crate::blastinput::app::formatting_handler_check(&choice, false)?;
        if self.max_target_seqs.is_some_and(|size| size < 5) {
            diagnostics.write_all(&crate::report::query_warnings::few_matches_warning(
                "blastp",
            ))?;
        }
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
        let (mode, unified_p) = crate::blastinput::app::parse_comp_based_stats(
            &self.comp_based_stats,
            true,
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
            evalue: options.evalue,
        })?;
        let matrix = options
            .matrix_name
            .parse::<ScoringMatrix>()
            .unwrap_or(ScoringMatrix::Blosum62);
        Ok(ResolvedBlastpArgs {
            query: self.query.clone(),
            subject: self.subject.clone().unwrap_or_default(),
            task: self.task.clone(),
            evalue: options.evalue,
            threshold: options.threshold,
            word_size: options.word_size as usize,
            lookup_table_type: options.lookup_table_type,
            num_threads: self.num_threads,
            out: self.out.clone(),
            max_target_seqs: self.max_target_seqs.unwrap_or(500),
            max_target_seqs_given: self.max_target_seqs,
            // NCBI blast_args.cpp:317-318: if (args[kArgMaxHSPsPerSubject])
            // opt.SetMaxHspsPerSubject(args[kArgMaxHSPsPerSubject].AsInteger());
            max_hsps_per_subject: self.max_hsps_per_subject.unwrap_or(0),
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
            use_sw_tback: self.use_sw_tback,
            chaining: options.chaining,
            outfmt: self.outfmt.clone(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    // NCBI blast_args.cpp:332-349: parse the same single string as the public CLI.
    fn parse_cli<const N: usize>(args: [&str; N]) -> TestCli {
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
        assert_eq!(resolved.outfmt, "0");
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
