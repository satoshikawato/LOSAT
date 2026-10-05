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
// The arguments are read with NCBI's grammar (`blastinput/value_parsers.rs`); the values
// that NCBI checks after parsing (-seg, -threshold, -word_size, -evalue, -outfmt) are
// checked in NCBI's order by `blast_engine::check_options`.
#[derive(Args, Debug, Clone)]
#[command(rename_all = "snake_case")]
pub struct TblastxArgs {
    // NCBI blast_args.cpp:3425-3427:
    // arg_desc.AddDefaultKey(kArgQuery, "input_file", "Input file name",
    //                        CArgDescriptions::eInputFile, kDfltArgQuery);
    // The default `-` is standard input.
    #[arg(long, value_parser = ncbi_input_path(), value_name = "PATH", default_value = "-")]
    pub query: PathBuf,
    #[arg(long, value_parser = ncbi_input_path(), value_name = "PATH")]
    pub subject: Option<PathBuf>,
    #[arg(long, default_value_t = 10.0, value_parser = tblastx_real)]
    pub evalue: f64,
    // NCBI blast_args.cpp:578-583: a double of at least 0; the lookup table takes its
    // (Int4) value and the report prints the double.
    #[arg(long, default_value_t = 13.0, value_parser = tblastx_threshold)]
    pub threshold: f64,
    #[arg(long, default_value_t = 3, value_parser = protein_word_size)]
    pub word_size: i32,
    #[arg(long, default_value_t = 1, value_parser = blastn_count)]
    pub num_threads: usize,

    #[arg(long, value_name = "PATH", value_parser = ncbi_output_path())]
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
    #[arg(long, value_parser = blastn_count, help = "Maximum number of aligned sequences to keep (default: 500)")]
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
    // Read when NCBI's filtering handler reads it (`blastinput/app.rs` `parse_seg_option`).
    #[arg(
        long,
        default_value = "12 2.2 2.5",
        help = "SEG: no, yes, or WINDOW LOCUT HICUT"
    )]
    pub seg: String,

    /// Two-hit window size for triggering ungapped extension (default: 40).
    /// Smaller values are more strict, larger values are more sensitive.
    /// 0 (NCBI's one-hit word finder) is not supported by LOSAT's TBLASTX.
    #[arg(long, default_value_t = 40, value_parser = nonnegative_ncbi_integer)]
    pub window_size: i32,

    /// Output format: 0 (pairwise), 6 or 7 (tabular), without custom fields.
    // Read when NCBI reads it (`blastinput/app.rs` `parse_formatting_string`).
    #[arg(long, default_value = "0", value_name = "SPEC")]
    pub outfmt: String,

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
    #[arg(long, default_value_t = 0, value_parser = nonnegative_ncbi_integer)]
    pub culling_limit: i32,
}

impl TblastxArgs {
    /// The -seg value as NCBI's filtering handler reads it (`check_options` reports its
    /// errors first).
    pub fn seg_spec(&self) -> anyhow::Result<SegSpec> {
        crate::blastinput::app::parse_seg_option(&self.seg, "TBLASTX")
    }

    /// The subject file name as the reports show it.
    pub fn subject_label(&self) -> std::path::Display<'_> {
        self.subject
            .as_deref()
            .unwrap_or(std::path::Path::new(""))
            .display()
    }
}
