//! Coordination module for BLASTN pipeline setup
//! 
//! This module handles:
//! - Task-specific parameter configuration (word size, scoring, etc.)
//! - Sequence reading and preprocessing
//! - DUST masking
//! - Lookup table building

use anyhow::Result;
use bio::io::fasta;
use crate::utils::dust::{DustMasker, MaskedInterval};
use super::args::BlastnArgs;
use super::constants::{MIN_UNGAPPED_SCORE_MEGABLAST, MIN_UNGAPPED_SCORE_BLASTN, MAX_DIRECT_LOOKUP_WORD_SIZE};
use super::lookup::{build_lookup, build_pv_direct_lookup, build_two_stage_lookup, TwoStageLookup, PvDirectLookup, KmerLookup};

/// Task-specific configuration for BLASTN
pub struct TaskConfig {
    pub effective_word_size: usize,
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub min_ungapped_score: i32,
    pub use_dp: bool,
    pub scan_step: usize,
    pub use_direct_lookup: bool,
    pub use_two_stage: bool,
    pub lut_word_length: usize,
}

/// Lookup tables for seed finding
pub struct LookupTables {
    pub two_stage_lookup: Option<TwoStageLookup>,
    pub pv_direct_lookup: Option<PvDirectLookup>,
    pub hash_lookup: Option<KmerLookup>,
}

/// Sequence data and metadata
pub struct SequenceData {
    pub queries: Vec<fasta::Record>,
    pub query_ids: Vec<String>,
    pub subjects: Vec<fasta::Record>,
    pub query_masks: Vec<Vec<MaskedInterval>>,
    pub db_len_total: usize,
    pub db_num_seqs: usize,
}

/// Determine effective word size based on task
pub fn determine_effective_word_size(args: &BlastnArgs) -> usize {
    match args.task.as_str() {
        "megablast" => args.word_size,
        "blastn" | "dc-megablast" => {
            if args.word_size == 28 {
                11
            } else {
                args.word_size
            }
        }
        _ => args.word_size,
    }
}

/// Determine effective scoring parameters based on task
pub fn determine_scoring_params(args: &BlastnArgs) -> (i32, i32, i32, i32) {
    match args.task.as_str() {
        "megablast" => {
            let r = if args.reward == 1 { 1 } else { args.reward };
            let p = if args.penalty == -2 { -2 } else { args.penalty };
            let go = if args.gap_open == 0 { 0 } else { args.gap_open };
            let ge = if args.gap_extend == 0 { 0 } else { args.gap_extend };
            (r, p, go, ge)
        }
        "blastn" | "dc-megablast" => {
            let r = if args.reward == 1 { 2 } else { args.reward };
            let p = if args.penalty == -2 { -3 } else { args.penalty };
            let go = if args.gap_open == 0 { -5 } else { args.gap_open };
            let ge = if args.gap_extend == 0 { -2 } else { args.gap_extend };
            (r, p, go, ge)
        }
        "blastn-short" => {
            let r = if args.reward == 1 { 1 } else { args.reward };
            let p = if args.penalty == -2 { -3 } else { args.penalty };
            let go = if args.gap_open == 0 { -5 } else { args.gap_open };
            let ge = if args.gap_extend == 0 { -2 } else { args.gap_extend };
            (r, p, go, ge)
        }
        _ => (args.reward, args.penalty, args.gap_open, args.gap_extend),
    }
}

/// Calculate initial scan step based on word size
pub fn calculate_initial_scan_step(effective_word_size: usize, user_scan_step: usize) -> usize {
    if user_scan_step > 0 {
        user_scan_step
    } else {
        if effective_word_size >= 16 {
            4
        } else if effective_word_size >= 11 {
            2
        } else {
            1
        }
    }
}

/// Configure task-specific parameters
pub fn configure_task(args: &BlastnArgs) -> TaskConfig {
    let effective_word_size = determine_effective_word_size(args);
    let (reward, penalty, gap_open, gap_extend) = determine_scoring_params(args);
    
    let min_ungapped_score = match args.task.as_str() {
        "megablast" => MIN_UNGAPPED_SCORE_MEGABLAST,
        _ => MIN_UNGAPPED_SCORE_BLASTN,
    };
    
    let use_dp = match args.task.as_str() {
        "megablast" => false,
        _ => true,
    };
    
    let scan_step = calculate_initial_scan_step(effective_word_size, args.scan_step);
    
    let use_direct_lookup = effective_word_size <= MAX_DIRECT_LOOKUP_WORD_SIZE;
    let use_two_stage = effective_word_size >= 11;
    let lut_word_length = if use_two_stage {
        if effective_word_size >= 16 {
            8
        } else if effective_word_size == 11 {
            8
        } else {
            8
        }
    } else {
        effective_word_size.min(MAX_DIRECT_LOOKUP_WORD_SIZE)
    };
    
    TaskConfig {
        effective_word_size,
        reward,
        penalty,
        gap_open,
        gap_extend,
        min_ungapped_score,
        use_dp,
        scan_step,
        use_direct_lookup,
        use_two_stage,
        lut_word_length,
    }
}

/// Read query and subject sequences
pub fn read_sequences(args: &BlastnArgs) -> Result<(Vec<fasta::Record>, Vec<String>, Vec<fasta::Record>)> {
    eprintln!("Reading query & subject...");
    let query_reader = fasta::Reader::from_file(&args.query)?;
    let queries: Vec<fasta::Record> = query_reader.records().filter_map(|r| r.ok()).collect();
    let query_ids: Vec<String> = queries
        .iter()
        .map(|r| {
            r.id()
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string()
        })
        .collect();
    
    let subject_reader = fasta::Reader::from_file(&args.subject)?;
    let subjects: Vec<fasta::Record> = subject_reader.records().filter_map(|r| r.ok()).collect();
    
    Ok((queries, query_ids, subjects))
}

/// Apply DUST filter to query sequences
pub fn apply_dust_masking(
    args: &BlastnArgs,
    queries: &[fasta::Record],
) -> Vec<Vec<MaskedInterval>> {
    if args.dust {
        eprintln!(
            "Applying DUST filter (level={}, window={}, linker={})...",
            args.dust_level, args.dust_window, args.dust_linker
        );
        let masker = DustMasker::new(args.dust_level, args.dust_window, args.dust_linker);
        let masks: Vec<Vec<MaskedInterval>> = queries
            .iter()
            .map(|record| masker.mask_sequence(record.seq()))
            .collect();
        
        let total_masked: usize = masks.iter().map(|m| m.iter().map(|i| i.end - i.start).sum::<usize>()).sum();
        let total_bases: usize = queries.iter().map(|r| r.seq().len()).sum();
        if total_bases > 0 {
            eprintln!(
                "DUST masked {} bases ({:.2}%) across {} sequences",
                total_masked,
                100.0 * total_masked as f64 / total_bases as f64,
                queries.len()
            );
        }
        masks
    } else {
        vec![Vec::new(); queries.len()]
    }
}

/// Build lookup tables based on configuration
pub fn build_lookup_tables(
    config: &TaskConfig,
    args: &BlastnArgs,
    queries: &[fasta::Record],
    query_masks: &[Vec<MaskedInterval>],
) -> (LookupTables, usize) {
    eprintln!(
        "Building lookup (Task: {}, Word: {}, TwoStage: {}, LUTWord: {}, Direct: {}, DUST: {})...",
        args.task, config.effective_word_size, config.use_two_stage, config.lut_word_length, config.use_direct_lookup, args.dust
    );
    
    let two_stage_lookup: Option<TwoStageLookup> = if config.use_two_stage {
        Some(build_two_stage_lookup(queries, config.effective_word_size, config.lut_word_length, query_masks))
    } else {
        None
    };
    
    let pv_direct_lookup: Option<PvDirectLookup> = if !config.use_two_stage && config.use_direct_lookup {
        Some(build_pv_direct_lookup(queries, config.effective_word_size, query_masks))
    } else {
        None
    };
    
    let hash_lookup: Option<KmerLookup> = if !config.use_two_stage && !config.use_direct_lookup {
        Some(build_lookup(queries, config.effective_word_size, query_masks))
    } else {
        None
    };
    
    // Update scan_step for two-stage lookup
    let scan_step = if config.use_two_stage {
        (config.effective_word_size as isize - config.lut_word_length as isize + 1).max(1) as usize
    } else {
        config.scan_step
    };
    
    if args.verbose && config.use_two_stage {
        eprintln!("[INFO] Using two-stage lookup: lut_word_length={}, word_length={}, scan_step={}", 
                  config.lut_word_length, config.effective_word_size, scan_step);
    }
    
    (
        LookupTables {
            two_stage_lookup,
            pv_direct_lookup,
            hash_lookup,
        },
        scan_step,
    )
}

/// Prepare all sequence data and configuration
pub fn prepare_sequence_data(
    args: &BlastnArgs,
    queries: Vec<fasta::Record>,
    query_ids: Vec<String>,
    subjects: Vec<fasta::Record>,
) -> SequenceData {
    let query_masks = apply_dust_masking(args, &queries);
    let db_len_total: usize = subjects.iter().map(|r| r.seq().len()).sum();
    let db_num_seqs: usize = subjects.len();
    
    SequenceData {
        queries,
        query_ids,
        subjects,
        query_masks,
        db_len_total,
        db_num_seqs,
    }
}

