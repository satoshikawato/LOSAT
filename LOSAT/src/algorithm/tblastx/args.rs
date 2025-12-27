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
}

