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
#[derive(Args, Debug, Clone)]
#[command(rename_all = "snake_case")]
pub struct BlastpArgs {
    #[arg(long, value_parser = file_path(), value_name = "PATH")]
    pub query: PathBuf,
    #[arg(long, value_parser = file_path(), value_name = "PATH")]
    pub subject: PathBuf,
    // NCBI api/blast_options_handle.cpp:381-387:
    // CBlastOptionsFactory::Create(eBlastp, locality);
    // v0.1.0 exposes only the comparison-supported ordinary task.
    #[arg(long, default_value = "blastp", value_parser = ["blastp"])]
    pub task: String,
    // NCBI c++/src/algo/blast/api/blast_options_handle.cpp:388-394:
    // if (task == "blastp-short") { ... opts->SetEvalueThreshold(20000); }
    // Keep omission distinct so resolve() retains the selected task's default.
    #[arg(long, value_parser = nonnegative_f64, help = "[default: 10]")]
    pub evalue: Option<f64>,
    #[arg(long, value_parser = positive_f64, help = "[default: resolved from configuration]", long_help = "Default blastp configuration: BLOSUM62, word size 3, threshold 11, window 40, gaps 11/1. Explicit matrix/word size overrides resolve their associated defaults.")]
    pub threshold: Option<f64>,
    #[arg(long, value_parser = blastp_word_size, help = "[default: resolved from configuration]", long_help = "Default blastp configuration: BLOSUM62, word size 3, threshold 11, window 40, gaps 11/1. Explicit matrix/word size overrides resolve their associated defaults.")]
    pub word_size: Option<usize>,
    #[arg(long, default_value_t = 1, value_parser = positive_usize)]
    pub num_threads: usize,
    #[arg(long, value_name = "PATH")]
    pub out: Option<PathBuf>,
    #[arg(long, default_value_t = 500, value_parser = positive_usize)]
    pub max_target_seqs: usize,
    #[arg(long = "max_hsps", value_parser = positive_usize)]
    pub max_hsps_per_subject: Option<usize>,
    #[arg(long)]
    pub ungapped: bool,
    #[arg(long, help = "[default: resolved from configuration]", long_help = "Default blastp configuration: BLOSUM62, word size 3, threshold 11, window 40, gaps 11/1. Explicit matrix/word size overrides resolve their associated defaults.", value_parser = nonnegative_usize)]
    pub window_size: Option<usize>,
    #[arg(
        long,
        help = "[default: resolved from configuration]",
        long_help = "Default blastp configuration: BLOSUM62, word size 3, threshold 11, window 40, gaps 11/1. Explicit matrix/word size overrides resolve their associated defaults."
    )]
    pub matrix: Option<String>,
    #[arg(long = "gapopen", value_parser = nonnegative_i32, help = "[default: resolved from configuration]", long_help = "Default blastp configuration: BLOSUM62, word size 3, threshold 11, window 40, gaps 11/1. Explicit matrix/word size overrides resolve their associated defaults.")]
    pub gap_open: Option<i32>,
    #[arg(long = "gapextend", value_parser = nonnegative_i32, help = "[default: resolved from configuration]", long_help = "Default blastp configuration: BLOSUM62, word size 3, threshold 11, window 40, gaps 11/1. Explicit matrix/word size overrides resolve their associated defaults.")]
    pub gap_extend: Option<i32>,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
    // ```c
    // const string kArgCompBasedStats("comp_based_stats");
    // ```
    #[arg(long, value_parser = parse_comp_based_stats)]
    pub comp_based_stats: Option<BlastpCompBasedStats>,
    #[arg(
        long,
        value_parser = parse_seg_filtering,
        action = clap::ArgAction::Set,
        num_args = 1
    , help = "SEG: no, yes, or \"WINDOW LOCUT HICUT\" [default: no]")]
    pub seg: Option<BlastpSegSpec>,
    // NCBI api/blast_advprot_options.cpp:58:
    // m_Opts->SetSmithWatermanMode(false);
    // Retain internal configuration; this unported path is not public in v0.1.0.
    #[arg(skip)]
    pub use_sw_tback: bool,
    #[arg(long, default_value = "0", value_name = "SPEC", value_parser = super::blast_engine::validate_cli_outfmt)]
    pub outfmt: String,
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_options_handle.cpp:388-400
// ```c
// else if (!NStr::CompareNocase(task, "blastp") ||
//          !NStr::CompareNocase(task, "blastp-short") ||
//          !NStr::CompareNocase(task, "blastp-fast"))
// {
//      CBlastAdvancedProteinOptionsHandle* opts =
//            dynamic_cast<CBlastAdvancedProteinOptionsHandle*>
//             (CBlastOptionsFactory::Create(eBlastp, locality));
//      if (task == "blastp-short") {
//         opts->SetMatrixName("PAM30");
//         opts->SetGapOpeningCost(9);
//         opts->SetGapExtensionCost(1);
//         opts->SetEvalueThreshold(20000);
//         opts->SetWordSize(2);
//         opts->ClearFilterOptions();
//      } else if (task == "blastp-fast") {
//         opts->SetWordSize(5);
//         opts->SetOptions().SetLookupTableType(eCompressedAaLookupTable);
//         opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP_FAST);
//         opts->SetChaining(true);
//      }
//      retval = opts;
// }
// ```
#[derive(Debug, Clone, Copy)]
struct BlastpTaskDefaults {
    task: &'static str,
    evalue: f64,
    matrix: ScoringMatrix,
    word_size: usize,
    chaining: bool,
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_options_handle.cpp:388-400
// ```c
// if (task == "blastp-short") {
//    opts->SetMatrixName("PAM30");
//    opts->SetGapOpeningCost(9);
//    opts->SetGapExtensionCost(1);
//    opts->SetEvalueThreshold(20000);
//    opts->SetWordSize(2);
//    opts->ClearFilterOptions();
// } else if (task == "blastp-fast") {
//    opts->SetWordSize(5);
//    opts->SetOptions().SetLookupTableType(eCompressedAaLookupTable);
//    opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP_FAST);
//    opts->SetChaining(true);
// }
// ```
fn blastp_task_defaults(task: &str) -> Result<BlastpTaskDefaults> {
    match task {
        "blastp" => Ok(BlastpTaskDefaults {
            task: "blastp",
            evalue: 10.0,
            matrix: ScoringMatrix::Blosum62,
            word_size: 3,
            chaining: false,
        }),
        "blastp-short" => Ok(BlastpTaskDefaults {
            task: "blastp-short",
            evalue: 20000.0,
            matrix: ScoringMatrix::Pam30,
            word_size: 2,
            chaining: false,
        }),
        "blastp-fast" => Ok(BlastpTaskDefaults {
            task: "blastp-fast",
            evalue: 10.0,
            matrix: ScoringMatrix::Blosum62,
            word_size: 5,
            chaining: true,
        }),
        _ => bail!(
            "unsupported blastp task '{}': expected one of blastp, blastp-short, blastp-fast",
            task
        ),
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_prot_options.cpp:55-83
// ```c
// if (ws > 4) {
//     m_Opts->SetLookupTableType(eCompressedAaLookupTable);
// }
// else {
//     m_Opts->SetLookupTableType(eAaLookupTable);
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:275-285
// ```c
// if (m_QueryIsProtein && args[kArgWordSize].AsInteger() > 4){
//     opt.SetLookupTableType(eCompressedAaLookupTable);
//     opt.SetWordThreshold(19.3);
//     if (args[kArgWordSize].AsInteger() > 5) {
//         opt.SetWordThreshold(21.0);
//     }
//     if (args[kArgWordSize].AsInteger() > 6) {
//         opt.SetWordThreshold(20.25);
//     }
// }
// ```
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
    pub max_target_seqs: usize,
    pub max_hsps_per_subject: usize,
    pub ungapped: bool,
    pub window_size: usize,
    pub scoring: ProteinScoringSpec,
    pub comp_based_stats: BlastpCompBasedStats,
    pub seg: BlastpSegSpec,
    pub use_sw_tback: bool,
    pub chaining: bool,
    pub outfmt: String,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_options.c:1174-1198
// ```c
// if(strcasecmp(matrixName, "BLOSUM62") == 0)
//     *threshold = kB62_threshold;
// else if(strcasecmp(matrixName, "BLOSUM45") == 0)
//     *threshold = 14;
// else if(strcasecmp(matrixName, "BLOSUM80") == 0)
//     *threshold = 12;
// else if(strcasecmp(matrixName, "PAM30") == 0)
//     *threshold = 16;
// else if(strcasecmp(matrixName, "PAM70") == 0)
//     *threshold = 14;
// else
//     *threshold = kB62_threshold;
// ```
#[inline]
fn suggested_threshold(matrix: ScoringMatrix) -> f64 {
    match matrix {
        ScoringMatrix::Blosum45 => 14.0,
        ScoringMatrix::Blosum80 => 12.0,
        ScoringMatrix::Pam30 => 16.0,
        ScoringMatrix::Pam70 => 14.0,
        _ => 11.0,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_options.c:1200-1227
// ```c
// if(strcasecmp(matrixName, "BLOSUM62") == 0)
//     *window_size = kB62_windowsize;
// else if(strcasecmp(matrixName, "BLOSUM45") == 0)
//     *window_size = 60;
// else if(strcasecmp(matrixName, "BLOSUM80") == 0)
//     *window_size = 25;
// else if(strcasecmp(matrixName, "PAM30") == 0)
//     *window_size = 15;
// else if(strcasecmp(matrixName, "PAM70") == 0)
//     *window_size = 20;
// else
//     *window_size = kB62_windowsize;
// ```
#[inline]
fn suggested_window_size(matrix: ScoringMatrix) -> usize {
    match matrix {
        ScoringMatrix::Blosum45 => 60,
        ScoringMatrix::Blosum80 => 25,
        ScoringMatrix::Pam30 => 15,
        ScoringMatrix::Pam70 => 20,
        _ => 40,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:200-217,238-255,273-287,303-314,328-337,359-376,393-404,421-430
// ```c
// static Int4 blosum45_prefs[BLOSUM45_VALUES_MAX] = { ..., BLAST_MATRIX_BEST, ... };
// static Int4 blosum50_prefs[BLOSUM50_VALUES_MAX] = { ..., BLAST_MATRIX_BEST, ... };
// static Int4 blosum62_prefs[BLOSUM62_VALUES_MAX] = { ..., BLAST_MATRIX_BEST, ... };
// static Int4 blosum80_prefs[BLOSUM80_VALUES_MAX] = { ..., BLAST_MATRIX_BEST, ... };
// static Int4 blosum90_prefs[BLOSUM90_VALUES_MAX] = { ..., BLAST_MATRIX_BEST, ... };
// static Int4 pam250_prefs[PAM250_VALUES_MAX] = { ..., BLAST_MATRIX_BEST, ... };
// static Int4 pam30_prefs[PAM30_VALUES_MAX] = { ..., BLAST_MATRIX_BEST, ... };
// static Int4 pam70_prefs[PAM70_VALUES_MAX] = { ..., BLAST_MATRIX_BEST, ... };
// ```
#[inline]
fn preferred_gap_costs(matrix: ScoringMatrix) -> (i32, i32) {
    match matrix {
        ScoringMatrix::Blosum45 => (14, 2),
        ScoringMatrix::Blosum50 => (13, 2),
        ScoringMatrix::Blosum62 => (11, 1),
        ScoringMatrix::Blosum80 => (10, 1),
        ScoringMatrix::Blosum90 => (10, 1),
        ScoringMatrix::Pam30 => (9, 1),
        ScoringMatrix::Pam70 => (10, 1),
        ScoringMatrix::Pam250 => (15, 2),
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_prot_options.cpp:57-83
// ```c
// switch (ws) {
// case 3: m_Opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP); break;
// case 5: m_Opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP_FAST); break;
// case 6: m_Opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP_WD_SZ_6); break;
// case 7: m_Opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP_WD_SZ_7); break;
// default: m_Opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP); break;
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:275-285
// ```c
// if (m_QueryIsProtein && args[kArgWordSize].AsInteger() > 4){
//     opt.SetLookupTableType(eCompressedAaLookupTable);
//     opt.SetWordThreshold(19.3);
//     if (args[kArgWordSize].AsInteger() > 5) {
//         opt.SetWordThreshold(21.0);
//     }
//     if (args[kArgWordSize].AsInteger() > 6) {
//         opt.SetWordThreshold(20.25);
//     }
// }
// ```
#[inline]
fn word_size_threshold_override(word_size: usize) -> Option<f64> {
    if word_size > 6 {
        Some(20.25)
    } else if word_size > 5 {
        Some(21.0)
    } else if word_size > 4 {
        Some(19.3)
    } else {
        None
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:257-273,275-285,491-499,607-620,893-899
// ```c
// if (args.Exist(kArgMatrixName) && args[kArgMatrixName])
//     BLAST_GetProteinGapExistenceExtendParams(args[kArgMatrixName].AsString().c_str(), &gap_open, &gap_extend);
// ...
// else if (args.Exist(kArgMatrixName) && args[kArgMatrixName]) {
//     opt.SetGapOpeningCost(gap_open);
// }
// ...
// if (args.Exist(kArgWordSize) && args[kArgWordSize]) {
//     if (m_QueryIsProtein && args[kArgWordSize].AsInteger() > 4){
//         opt.SetLookupTableType(eCompressedAaLookupTable);
//         opt.SetWordThreshold(19.3);
//         ...
//     }
//     opt.SetWordSize(args[kArgWordSize].AsInteger());
// }
// ...
// BLAST_GetSuggestedWindowSize(opt.GetProgramType(), opt.GetMatrixName(), &window);
// ...
// BLAST_GetSuggestedThreshold(opt.GetProgramType(), opt.GetMatrixName(), &threshold);
// ...
// s_SetCompositionBasedStats(opt, args[kArgCompBasedStats].AsString(), ...);
// ```
impl BlastpArgs {
    pub fn resolve(&self) -> Result<ResolvedBlastpArgs> {
        let canonical_task = self.task.to_ascii_lowercase();
        let task_defaults = blastp_task_defaults(&canonical_task)?;

        let matrix = match self.matrix.as_deref() {
            Some(name) => name.parse::<ScoringMatrix>().map_err(anyhow::Error::msg)?,
            None => task_defaults.matrix,
        };
        let word_size = self.word_size.unwrap_or(task_defaults.word_size);
        let lookup_table_type = if word_size > 4 {
            BlastpLookupTableType::CompressedAaLookupTable
        } else {
            BlastpLookupTableType::AaLookupTable
        };

        let (preferred_gap_open, preferred_gap_extend) = preferred_gap_costs(matrix);
        let scoring = ProteinScoringSpec {
            matrix,
            gap_open: match (self.gap_open, self.matrix.as_ref()) {
                (Some(value), _) => value,
                (None, Some(_)) => preferred_gap_open,
                (None, None) => preferred_gap_open,
            },
            gap_extend: match (self.gap_extend, self.matrix.as_ref()) {
                (Some(value), _) => value,
                (None, Some(_)) => preferred_gap_extend,
                (None, None) => preferred_gap_extend,
            },
        };

        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:867-876
        // ```c
        // if (ungapped && *ungapped && compo_mode != eNoCompositionBasedStats) {
        //     NCBI_THROW(CInputException, eInvalidInput,
        //                "Composition-adjusted searched are not supported with "
        //                "an ungapped search, please add -comp_based_stats F or "
        //                "do a gapped search");
        // }
        // ```
        let comp_based_stats = self.comp_based_stats.unwrap_or(BlastpCompBasedStats {
            mode: BlastpCompositionMode::CompositionMatrixAdjust,
            unified_p: false,
        });
        if self.ungapped && comp_based_stats.is_enabled() {
            bail!(
                "Composition-adjusted searched are not supported with an ungapped search, please add -comp_based_stats F or do a gapped search"
            );
        }

        let threshold = if let Some(value) = self.threshold {
            value
        } else if let Some(value) = word_size_threshold_override(word_size) {
            value
        } else {
            suggested_threshold(matrix)
        };

        let window_size = self
            .window_size
            .unwrap_or_else(|| suggested_window_size(matrix));
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_advprot_options.cpp:64-67
        // ```c
        // void
        // CBlastAdvancedProteinOptionsHandle::SetQueryOptionDefaults()
        // {
        //     CBlastProteinOptionsHandle::SetQueryOptionDefaults();
        //     SetSegFiltering(false); // disable SEG filtering because of eCompositionMatrixAdjust mode
        // }
        // ```
        let seg = self.seg.clone().unwrap_or(BlastpSegSpec::No);

        Ok(ResolvedBlastpArgs {
            query: self.query.clone(),
            subject: self.subject.clone(),
            task: task_defaults.task.to_string(),
            evalue: self.evalue.unwrap_or(task_defaults.evalue),
            threshold,
            word_size,
            lookup_table_type,
            num_threads: self.num_threads,
            out: self.out.clone(),
            max_target_seqs: self.max_target_seqs,
            // NCBI blast_args.cpp:317-318: if (args[kArgMaxHSPsPerSubject])
            // opt.SetMaxHspsPerSubject(args[kArgMaxHSPsPerSubject].AsInteger());
            max_hsps_per_subject: self.max_hsps_per_subject.unwrap_or(0),
            ungapped: self.ungapped,
            window_size,
            scoring,
            comp_based_stats,
            seg,
            use_sw_tback: self.use_sw_tback,
            chaining: task_defaults.chaining,
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
        assert_eq!(args.comp_based_stats, None);
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
        // NCBI api/blast_options_handle.cpp:388-400 preserves future internal task defaults.
        // opts->SetEvalueThreshold(20000); opts->SetWordSize(5);
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
        // NCBI api/blast_options_handle.cpp:388-400 preserves future internal task defaults.
        // opts->SetEvalueThreshold(20000); opts->SetWordSize(5);
        let cli = parse_cli(["losat", "blastp", "-query", "q.faa", "-subject", "s.faa"]);
        let TestCommand::Blastp(mut args) = cli.command;
        args.task = "blastp-fast".into();
        let resolved = args.resolve().expect("resolved blastp-fast args");
        assert_eq!(resolved.task, "blastp-fast");
        assert_eq!(resolved.evalue, 10.0);
        assert_eq!(resolved.word_size, 5);
        assert_eq!(resolved.threshold, 19.3);
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
        assert_eq!(args.seg, Some(BlastpSegSpec::No));

        let cli = parse_cli([
            "losat", "blastp", "-query", "q.faa", "-subject", "s.faa", "-seg", "yes",
        ]);
        let TestCommand::Blastp(args) = cli.command;
        assert_eq!(args.seg, Some(BlastpSegSpec::Yes));

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
            args.seg,
            Some(BlastpSegSpec::WindowLocutHicut {
                window: 15,
                locut: 2.5,
                hicut: 3.0,
            })
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
            "-gapopen", "12",
        ]);
        let TestCommand::Blastp(args) = cli.command;
        let resolved = args.resolve().expect("resolved blastp args");
        assert_eq!(resolved.scoring.matrix, ScoringMatrix::Pam70);
        assert_eq!(resolved.scoring.gap_open, 12);
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
