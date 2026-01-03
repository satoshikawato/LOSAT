use clap::Args;
use std::path::PathBuf;

#[derive(Args, Debug)]
pub struct BlastnArgs {
    #[arg(short, long)]
    pub query: PathBuf,
    #[arg(short, long)]
    pub subject: PathBuf,
    #[arg(long, default_value = "megablast")]
    pub task: String,
    #[arg(short, long, default_value_t = 28)]
    pub word_size: usize,
    #[arg(short = 'n', long, default_value_t = 0)]
    pub num_threads: usize,
    #[arg(long, default_value_t = 10.0)]
    pub evalue: f64,
    #[arg(long, default_value_t = 500)]
    pub max_target_seqs: usize,
    #[arg(short, long)]
    pub out: Option<PathBuf>,
    // Scoring parameters - defaults are for megablast task
    // For blastn task, these are overridden in run() based on --task
    #[arg(long, default_value_t = 1)]
    pub reward: i32,
    #[arg(long, default_value_t = -2)]
    pub penalty: i32,
    #[arg(long, default_value_t = 0)]
    pub gap_open: i32,
    #[arg(long, default_value_t = 0)]
    pub gap_extend: i32,
    // DUST filter options for masking low-complexity regions
    #[arg(long, default_value_t = true)]
    pub dust: bool,
    #[arg(long, default_value_t = 20)]
    pub dust_level: u32,
    #[arg(long, default_value_t = 64)]
    pub dust_window: usize,
    #[arg(long, default_value_t = 1)]
    pub dust_linker: usize,
    #[arg(long, short = 'v', default_value_t = false)]
    pub verbose: bool,
    /// Enable HSP chaining to merge nearby HSPs into longer alignments.
    /// By default, chaining is disabled for BLAST-compatible output (individual HSPs).
    /// Enable this for visualization tools like gbdraw that benefit from longer continuous alignments.
    #[arg(long, default_value_t = false)]
    pub chain: bool,
    /// Scan stride for subject sequence scanning (NCBI BLAST optimization).
    /// Higher values skip more positions, reducing k-mer lookups but potentially missing some seeds.
    /// Default: 0 (auto-calculate based on word_size: 1 for word_size < 16, 4 for word_size >= 16).
    /// For megablast (word_size=28), scan_step=4 reduces lookups by ~4x with minimal sensitivity loss.
    #[arg(long, default_value_t = 0)]
    pub scan_step: usize,

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
}

