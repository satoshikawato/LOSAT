//! Command-line arguments for TBLASTX

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
#[derive(Args, Debug)]
pub struct TblastxArgs {
    #[arg(short, long)]
    pub query: PathBuf,
    #[arg(short, long)]
    pub subject: PathBuf,
    #[arg(short, long, default_value_t = 10.0)]
    pub evalue: f64,
    #[arg(short, long, default_value_t = 13)]
    pub threshold: i32,
    #[arg(short, long, default_value_t = 3)]
    pub word_size: usize,
    #[arg(short = 'n', long, default_value_t = 0)]
    pub num_threads: usize,

    #[arg(short, long)]
    pub out: Option<PathBuf>,
    #[arg(long, default_value_t = 1)]
    pub query_gencode: u8,
    #[arg(long, default_value_t = 1)]
    pub db_gencode: u8,
    #[arg(long, default_value_t = 500)]
    pub max_target_seqs: usize,
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

    // SEG filter options for masking low-complexity regions in amino acid sequences
    // NCBI BLAST default: enabled for translated/protein searches.
    // To disable (for debugging/perf experiments), pass `--seg=false`.
    #[arg(
        long,
        default_value_t = true,
        action = clap::ArgAction::Set,
        num_args = 0..=1,
        default_missing_value = "true"
    )]
    pub seg: bool,
    #[arg(long, default_value_t = 12)]
    pub seg_window: usize,
    #[arg(long, default_value_t = 2.2)]
    pub seg_locut: f64,
    #[arg(long, default_value_t = 2.5)]
    pub seg_hicut: f64,

    /// Two-hit window size for triggering ungapped extension (default: 40)
    /// Smaller values are more strict, larger values are more sensitive
    /// Use 0 to enable one-hit mode (like NCBI BLAST's -window_size 0)
    #[arg(long, default_value_t = 40)]
    pub window_size: usize,

    /// Output format (NCBI BLAST compatible).
    ///
    /// Supported formats:
    ///   0 = Pairwise alignment view (traditional BLAST output)
    ///   6 = Tabular (tab-separated values, default)
    ///   7 = Tabular with comment lines (headers)
    ///
    /// Custom field specification is supported for formats 6 and 7:
    ///   -outfmt "6 qaccver saccver pident length"
    ///
    /// Default fields: qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
    #[arg(long, default_value = "6")]
    pub outfmt: String,

    /// HSP culling limit (number of HSPs allowed per query region).
    ///
    /// When > 0, applies NCBI's interval tree-based HSP culling algorithm to remove
    /// dominated HSPs based on score/length tradeoff. Default: 0 (disabled, matches NCBI tblastx default).
    ///
    /// NCBI reference: hspfilter_culling.c, cmdline_flags.cpp:127-128 (kDfltArgCullingLimit = 0)
    #[arg(long, default_value_t = 0)]
    pub culling_limit: u32,
}


