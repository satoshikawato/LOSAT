pub mod algorithm;
pub mod common;
pub mod sequence;
pub mod utils;

pub mod align;
pub mod config;
pub mod post;
pub mod report;
pub mod seed;
pub mod stats;

// New NCBI-style modules (matching NCBI BLAST structure)
pub mod api;
pub mod blastinput;
pub mod core;
pub mod format;

#[cfg(target_arch = "wasm32")]
pub mod web_api;
