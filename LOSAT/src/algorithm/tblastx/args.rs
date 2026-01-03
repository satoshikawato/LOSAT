//! Command-line arguments for TBLASTX

use clap::Args;
use std::path::PathBuf;

/// Command-line arguments for TBLASTX (translated DNA vs translated DNA search)
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

    /// Include stop codon (`*`) in 3-mer seed lookup indexing (experimental).
    ///
    /// Default behavior (NCBI parity) includes stop codons in lookup indexing.
    /// To opt out (legacy LOSAT safety behavior), pass `--include-stop-seeds=false`.
    ///
    /// NOTE: This affects seeding and can change the hit set.
    #[arg(
        long,
        default_value_t = true,
        action = clap::ArgAction::Set,
        num_args = 0..=1,
        default_missing_value = "true"
    )]
    pub include_stop_seeds: bool,

    /// Use NCBI BLOSUM62 stop-stop score (`*-* = +1`) during seeding/extension (experimental).
    ///
    /// Default behavior (NCBI parity) uses `*-* = +1`.
    /// To opt out (legacy LOSAT safety behavior), pass `--ncbi-stop-stop-score=false`.
    ///
    /// NOTE: This can change extensions; always re-check for overlong-tail regressions.
    #[arg(
        long,
        default_value_t = true,
        action = clap::ArgAction::Set,
        num_args = 0..=1,
        default_missing_value = "true"
    )]
    pub ncbi_stop_stop_score: bool,

    /// Maximum number of lookup hits per k-mer before suppressing that k-mer (high-frequency word suppression).
    ///
    /// NCBI's standard protein lookup construction does not apply this suppression.
    /// Default: disabled (usize::MAX) to match NCBI BLAST+ behavior exactly.
    #[arg(long, default_value_t = usize::MAX)]
    pub max_hits_per_kmer: usize,

    /// Use pre-computed neighbor map for fast scanning.
    ///
    /// When enabled, neighbor relationships are pre-computed once and reused during
    /// every scan. This drastically reduces lookup table size (from ~8M to ~567K entries)
    /// while maintaining NCBI-compatible output.
    ///
    /// The neighbor map computation takes a few seconds but enables much faster scanning.
    #[arg(long, default_value_t = false)]
    pub neighbor_map: bool,

    /// Use NCBI BLAST compatible parameters (TWO_HIT_WINDOW=16, X_DROP_UNGAPPED=7)
    /// When enabled, uses stricter parameters that match NCBI BLAST+ defaults
    #[arg(long, default_value_t = false)]
    pub ncbi_compat: bool,

    /// Two-hit window size for triggering ungapped extension (default: 40)
    /// Smaller values are more strict, larger values are more sensitive
    /// Use 0 to enable one-hit mode (like NCBI BLAST's -window_size 0)
    #[arg(long, default_value_t = 40)]
    pub window_size: usize,

    /// Restrict search to a single query translation frame (1,2,3,-1,-2,-3).
    /// Useful for debugging slow one-hit (`--window-size 0`) runs.
    #[arg(long, allow_hyphen_values = true)]
    pub only_qframe: Option<i8>,

    /// Restrict search to query strand only.
    /// Values: 1 (plus strand; frames 1,2,3) or -1 (minus strand; frames -1,-2,-3).
    ///
    /// NOTE: If `--only-qframe` is provided, it takes precedence.
    #[arg(long, allow_hyphen_values = true)]
    pub only_qstrand: Option<i8>,

    /// Restrict search to a single subject translation frame (1,2,3,-1,-2,-3).
    /// Useful for debugging slow one-hit (`--window-size 0`) runs.
    #[arg(long, allow_hyphen_values = true)]
    pub only_sframe: Option<i8>,

    /// Restrict search to subject strand only.
    /// Values: 1 (plus strand; frames 1,2,3) or -1 (minus strand; frames -1,-2,-3).
    ///
    /// NOTE: If `--only-sframe` is provided, it takes precedence.
    #[arg(long, allow_hyphen_values = true)]
    pub only_sstrand: Option<i8>,
}


