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
    #[arg(long, default_value_t = true)]
    pub dust: bool,
    #[arg(long, default_value_t = 20)]
    pub dust_level: u32,
    #[arg(long, default_value_t = 64)]
    pub dust_window: usize,
    #[arg(long, default_value_t = 1)]
    pub dust_linker: usize,
    
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

