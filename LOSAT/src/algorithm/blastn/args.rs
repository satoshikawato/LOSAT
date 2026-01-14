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
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2960-2968
    // ```c
    // if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
    //    m_NumDescriptions = args[kArgMaxTargetSequences].AsInteger();
    //    m_NumAlignments = args[kArgMaxTargetSequences].AsInteger();
    //    hitlist_size = m_NumAlignments;
    // }
    // ```
    #[arg(long)]
    pub max_target_seqs: Option<usize>,
    /// Maximum number of hits to save (NCBI BLAST hitlist_size)
    /// Reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:231-270
    #[arg(long, default_value_t = 500)]
    pub hitlist_size: usize,
    /// Remove word seeds with high frequency in the searched database.
    /// Reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:257 (limit_lookup)
    #[arg(long = "limit_lookup", default_value_t = false)]
    pub limit_lookup: bool,
    /// Maximum database word count for lookup filtering.
    /// Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:172-174
    /// #define MAX_DB_WORD_COUNT_MAPPER 30
    #[arg(long = "max_db_word_count", default_value_t = 30)]
    pub max_db_word_count: u8,
    /// Maximum number of HSPs per subject (0 = unlimited)
    /// Reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:231-270
    #[arg(long, default_value_t = 0)]
    pub max_hsps_per_subject: usize,
    /// Minimum diagonal separation between HSPs on the same subject (0 = auto, task-specific)
    /// Reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:231-270
    /// Default: 50 for blastn, 6 for megablast
    #[arg(long, default_value_t = 0)]
    pub min_diag_separation: usize,
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
    /// Apply subject best hit filtering (disabled by default).
    /// Reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:135
    #[arg(long = "subject_besthit", default_value_t = false)]
    pub subject_besthit: bool,
    #[arg(long, short = 'v', default_value_t = false)]
    pub verbose: bool,
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

