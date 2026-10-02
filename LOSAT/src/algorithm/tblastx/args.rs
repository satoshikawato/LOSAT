//! Command-line arguments for TBLASTX

use crate::blastinput::value_parsers::*;
use clap::Args;
use std::path::PathBuf;

/// Command-line arguments for TBLASTX (translated DNA vs translated DNA search)
///
/// NCBI reference (CTblastxAppArgs argument set):
/// ncbi-blast/c++/src/algo/blast/blastinput/tblastx_args.cpp:66-118
/// ```c
/// arg.Reset(new CGenericSearchArgs( !kQueryIsProtein, false, false, true));
/// ...
/// m_HspFilteringArgs.Reset(new CHspFilteringArgs);
/// ...
/// arg.Reset(new CWindowSizeArg);
/// ...
/// m_QueryOptsArgs.Reset(new CQueryOptionsArgs(kQueryIsProtein));
/// ...
/// arg.Reset(new CGeneticCodeArgs(CGeneticCodeArgs::eQuery));
/// arg.Reset(new CGeneticCodeArgs(CGeneticCodeArgs::eDatabase));
/// ...
/// arg.Reset(new CFormattingArgs);
/// arg.Reset(new CMTArgs);
/// arg.Reset(new CRemoteArgs);
/// arg.Reset(new CDebugArgs);
/// ```
// CLI v2 names/defaults: NCBI c++/src/algo/blast/blastinput/blast_args.cpp:
// 166-170,203-207,332-349,2657-2660,3158-3163; cmdline_flags.cpp:46-94.
// arg_desc.AddOptionalKey(kArgMaxHSPsPerSubject, "int_value", ..., eInteger);
// arg_desc.SetConstraint(kArgMaxHSPsPerSubject, new CArgAllowValuesGreaterThanOrEqual(1));
// arg_desc.AddDefaultKey(kArgNumThreads, "int_value", ..., NStr::IntToString(kDfltValue));
// The single-dash lexical translation is owned by crate::cli.
#[derive(Args, Debug)]
#[command(rename_all = "snake_case")]
pub struct TblastxArgs {
    #[arg(long, value_parser = file_path(), value_name = "PATH")]
    pub query: PathBuf,
    #[arg(long, value_parser = file_path(), value_name = "PATH")]
    pub subject: PathBuf,
    #[arg(long, default_value_t = 10.0, value_parser = nonnegative_f64)]
    pub evalue: f64,
    #[arg(long, default_value_t = 13, value_parser = positive_i32)]
    pub threshold: i32,
    #[arg(long, default_value_t = 3, value_parser = tblastx_word_size)]
    pub word_size: usize,
    #[arg(long, default_value_t = 1, value_parser = positive_usize)]
    pub num_threads: usize,

    #[arg(long, value_name = "PATH")]
    pub out: Option<PathBuf>,
    #[arg(long, default_value_t = 1, value_parser = genetic_code)]
    pub query_gencode: u8,
    #[arg(long, default_value_t = 1, value_parser = genetic_code)]
    pub db_gencode: u8,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2910-2927
    // ```c
    //     m_NumDescriptions = m_DfltNumDescriptions;
    //     m_NumAlignments = m_DfltNumAlignments;
    //     ...
    //     if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
    //         m_NumDescriptions = args[kArgMaxTargetSequences].AsInteger();
    //         m_NumAlignments = args[kArgMaxTargetSequences].AsInteger();
    //         hitlist_size = m_NumAlignments;
    //     }
    // ```
    // NCBI reference: c++/src/objtools/align_format/format_flags.cpp:219,221
    // ```c
    // const size_t kDfltArgNumDescriptions = 500;
    // const size_t kDfltArgNumAlignments = 250;
    // ```
    // An omitted value keeps the default hit list size (500) but the pairwise report then
    // shows 250 alignments, so the option has no clap default.
    #[arg(long, value_parser = positive_usize, help = "Maximum number of aligned sequences to keep (default: 500)")]
    pub max_target_seqs: Option<usize>,
    // NCBI low-complexity filtering selection:
    // - dust is used only for blastn (and mapping)
    // - otherwise seg is used
    //
    // NCBI reference (verbatim):
    //   else if (*ptr == 'L' || *ptr == 'T')
    //   { /* do low-complexity filtering; dust for blastn, otherwise seg.*/
    //       if (program_number == eBlastTypeBlastn
    //           || program_number == eBlastTypeMapping)
    //           SDustOptionsNew(&dustOptions);
    //       else
    //           SSegOptionsNew(&segOptions);
    //       ptr++;
    //   }
    // Source: ncbi-blast/c++/src/algo/blast/core/blast_filter.c:572-580
    //
    // Therefore, for tblastx we do NOT apply nucleotide-level DUST masking.

    // NCBI blast_args.cpp:396-406: opt.SetSegFiltering(false/true);
    // opt.SetSegFilteringWindow(...); opt.SetSegFilteringLocut(...);
    // opt.SetSegFilteringHicut(...);
    #[arg(long, default_value = "12 2.2 2.5", value_parser = parse_seg_filtering, help = "SEG: no, yes, or WINDOW LOCUT HICUT")]
    pub seg: SegSpec,

    /// Two-hit window size for triggering ungapped extension (default: 40).
    /// Smaller values are more strict, larger values are more sensitive.
    /// 0 (NCBI's one-hit word finder) is not supported by LOSAT's TBLASTX.
    #[arg(long, default_value_t = 40, value_parser = nonnegative_usize)]
    pub window_size: usize,

    /// Output format: 0 (pairwise), 6 or 7 (tabular), without custom fields.
    #[arg(long, default_value = "0", value_name = "SPEC", value_parser = tblastx_outfmt)]
    pub outfmt: String,

    // NCBI reference: cmdline_flags.cpp:127-128 (kDfltArgCullingLimit = 0)
    /// HSP culling limit (default: 0, no culling). A limit above 0 is not supported by
    /// LOSAT's TBLASTX (its HSP culling differs from NCBI's hspfilter_culling.c).
    #[arg(long, default_value_t = 0)]
    pub culling_limit: u32,
}
