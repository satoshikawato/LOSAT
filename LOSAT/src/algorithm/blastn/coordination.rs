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
use super::constants::{MIN_UNGAPPED_SCORE_MEGABLAST, MIN_UNGAPPED_SCORE_BLASTN, MAX_DIRECT_LOOKUP_WORD_SIZE, X_DROP_GAPPED_NUCL, X_DROP_GAPPED_GREEDY, X_DROP_GAPPED_FINAL, SCAN_RANGE_BLASTN, SCAN_RANGE_MEGABLAST, LUT_WIDTH_11_THRESHOLD_8, LUT_WIDTH_11_THRESHOLD_10, MIN_DIAG_SEPARATION_BLASTN, MIN_DIAG_SEPARATION_MEGABLAST};
use super::lookup::{
    build_db_word_counts, build_na_lookup, build_pv_direct_lookup, build_two_stage_lookup,
    NaLookupTable, PvDirectLookup, TwoStageLookup,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum LookupTableKind {
    Small,
    Mb,
    Na,
}

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
    pub x_drop_gapped: i32, // Task-specific gapped X-dropoff (blastn: 30, megablast: 25)
    pub x_drop_final: i32, // Final traceback X-dropoff (100 for all nucleotide tasks)
    pub scan_range: usize, // Scan range for off-diagonal hit detection (blastn: 4, megablast: 0)
    pub min_diag_separation: i32, // NCBI reference: blast_nucl_options.cpp:239,259 (blastn: 50, megablast: 6)
}

/// Lookup tables for seed finding
pub struct LookupTables {
    pub two_stage_lookup: Option<TwoStageLookup>,
    pub pv_direct_lookup: Option<PvDirectLookup>,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:131-156
    // ```c
    // typedef struct BlastNaLookupTable {
    //     ...
    // } BlastNaLookupTable;
    // ```
    pub na_lookup: Option<NaLookupTable>,
}

/// Subject metadata needed for search setup.
pub struct SubjectMetadata {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1409
    // ```c
    // db_length = BlastSeqSrcGetTotLen(seq_src);
    // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    // ```
    pub db_len_total: usize,
    pub db_num_seqs: usize,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-155
    // ```c
    // Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    // ```
    pub subject_ids: Vec<String>,
}

/// Sequence data and metadata
pub struct SequenceData {
    pub queries: Vec<fasta::Record>,
    pub query_ids: Vec<String>,
    pub query_masks: Vec<Vec<MaskedInterval>>,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1409
    // ```c
    // db_length = BlastSeqSrcGetTotLen(seq_src);
    // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    // ```
    pub db_len_total: usize,
    pub db_num_seqs: usize,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-155
    // ```c
    // Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    // ```
    pub subject_ids: Vec<String>,
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
            // NCBI BLAST: gap penalties are specified as positive values (cost)
            // Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:84-96
            let go = if args.gap_open == 0 { 5 } else { args.gap_open };
            let ge = if args.gap_extend == 0 { 2 } else { args.gap_extend };
            (r, p, go, ge)
        }
        "blastn-short" => {
            let r = if args.reward == 1 { 1 } else { args.reward };
            let p = if args.penalty == -2 { -3 } else { args.penalty };
            // NCBI BLAST: gap penalties are specified as positive values (cost)
            // Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:84-96
            let go = if args.gap_open == 0 { 5 } else { args.gap_open };
            let ge = if args.gap_extend == 0 { 2 } else { args.gap_extend };
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
    
    // NCBI BLAST algorithm selection:
    // - megablast: eGreedyScoreOnly (greedy alignment)
    // - blastn: eDynProgScoreOnly (dynamic programming)
    // Reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:182, 192
    let use_dp = match args.task.as_str() {
        "megablast" => false,
        _ => true,
    };
    
    // NCBI BLAST: Task-specific gapped X-dropoff
    // Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:122-148
    // blastn (non-greedy): 30, megablast (greedy): 25
    // Reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:177-194
    let x_drop_gapped = match args.task.as_str() {
        "megablast" => X_DROP_GAPPED_GREEDY, // 25
        _ => X_DROP_GAPPED_NUCL, // 30
    };

    // NCBI BLAST: Final traceback X-dropoff (100 for all nucleotide tasks)
    // Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:146
    // BLAST_GAP_X_DROPOFF_FINAL_NUCL = 100
    // Used in traceback phase to extend alignments further than preliminary extension
    let x_drop_final = X_DROP_GAPPED_FINAL; // 100
    
    let scan_step = calculate_initial_scan_step(effective_word_size, args.scan_step);
    
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:163-168
    // ```c
    // SetWindowSize(BLAST_WINDOW_SIZE_NUCL);
    // SetOffDiagonalRange(BLAST_SCAN_RANGE_NUCL);
    // ```
    let scan_range = match args.task.as_str() {
        "megablast" => SCAN_RANGE_MEGABLAST, // 0
        _ => SCAN_RANGE_BLASTN, // 0
    };

    // NCBI reference: blast_nucl_options.cpp:239, 259
    // Minimum diagonal separation for HSP containment checking
    // Used in MB_HSP_CLOSE macro (blast_gapalign_priv.h:123-124)
    let min_diag_separation = match args.task.as_str() {
        "megablast" => MIN_DIAG_SEPARATION_MEGABLAST, // 6
        _ => MIN_DIAG_SEPARATION_BLASTN, // 50
    };

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:45-185
    // ```c
    // lut_type = BlastChooseNaLookupTable(lookup_options,
    //                                     approx_table_entries,
    //                                     max_q_off,
    //                                     &lut_width);
    // ```
    let use_two_stage = false;
    let use_direct_lookup = effective_word_size <= MAX_DIRECT_LOOKUP_WORD_SIZE;
    let lut_word_length = effective_word_size;

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
        x_drop_gapped,
        x_drop_final,
        scan_range,
        min_diag_separation,
    }
}

/// Choose lookup table kind and LUT width using NCBI's BlastChooseNaLookupTable logic.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:45-185
/// ```c
/// if (lookup_options->mb_template_length > 0) {
///    *lut_width = lookup_options->word_size;
///    return eMBLookupTable;
/// }
/// switch(lookup_options->word_size) {
/// case 4: case 5: case 6:
///    lut_type = eSmallNaLookupTable;
///    *lut_width = lookup_options->word_size;
///    break;
/// case 7:
///    lut_type = eSmallNaLookupTable;
///    if (approx_table_entries < 250) *lut_width = 6;
///    else *lut_width = 7;
///    break;
/// case 8:
///    lut_type = eSmallNaLookupTable;
///    if (approx_table_entries < 8500) *lut_width = 7;
///    else *lut_width = 8;
///    break;
/// case 9:
///    if (approx_table_entries < 1250) { *lut_width = 7; lut_type = eSmallNaLookupTable; }
///    else if (approx_table_entries < 21000) { *lut_width = 8; lut_type = eSmallNaLookupTable; }
///    else { *lut_width = 9; lut_type = eMBLookupTable; }
///    break;
/// case 10:
///    if (approx_table_entries < 1250) { *lut_width = 7; lut_type = eSmallNaLookupTable; }
///    else if (approx_table_entries < 8500) { *lut_width = 8; lut_type = eSmallNaLookupTable; }
///    else if (approx_table_entries < 18000) { *lut_width = 9; lut_type = eMBLookupTable; }
///    else { *lut_width = 10; lut_type = eMBLookupTable; }
///    break;
/// case 11:
///    if (approx_table_entries < 12000) { *lut_width = 8; lut_type = eSmallNaLookupTable; }
///    else if (approx_table_entries < 180000) { *lut_width = 10; lut_type = eMBLookupTable; }
///    else { *lut_width = 11; lut_type = eMBLookupTable; }
///    break;
/// case 12:
///    if (approx_table_entries < 8500) { *lut_width = 8; lut_type = eSmallNaLookupTable; }
///    else if (approx_table_entries < 18000) { *lut_width = 9; lut_type = eMBLookupTable; }
///    else if (approx_table_entries < 60000) { *lut_width = 10; lut_type = eMBLookupTable; }
///    else if (approx_table_entries < 900000) { *lut_width = 11; lut_type = eMBLookupTable; }
///    else { *lut_width = 12; lut_type = eMBLookupTable; }
///    break;
/// default:
///    if (approx_table_entries < 8500) { *lut_width = 8; lut_type = eSmallNaLookupTable; }
///    else if (approx_table_entries < 300000) { *lut_width = 11; lut_type = eMBLookupTable; }
///    else { *lut_width = 12; lut_type = eMBLookupTable; }
///    break;
/// }
/// if (lut_type == eSmallNaLookupTable &&
///     (approx_table_entries >= 32767 || max_q_off >= 32768)) {
///    lut_type = eNaLookupTable;
/// }
/// ```
fn choose_na_lookup_table(
    word_size: usize,
    approx_table_entries: usize,
    max_q_off: usize,
    discontig_template: bool,
) -> (LookupTableKind, usize) {
    debug_assert!(word_size >= 4);

    if discontig_template {
        return (LookupTableKind::Mb, word_size);
    }

    let (mut lut_kind, mut lut_width) = match word_size {
        4 | 5 | 6 => (LookupTableKind::Small, word_size),
        7 => {
            if approx_table_entries < 250 {
                (LookupTableKind::Small, 6)
            } else {
                (LookupTableKind::Small, 7)
            }
        }
        8 => {
            if approx_table_entries < 8_500 {
                (LookupTableKind::Small, 7)
            } else {
                (LookupTableKind::Small, 8)
            }
        }
        9 => {
            if approx_table_entries < 1_250 {
                (LookupTableKind::Small, 7)
            } else if approx_table_entries < 21_000 {
                (LookupTableKind::Small, 8)
            } else {
                (LookupTableKind::Mb, 9)
            }
        }
        10 => {
            if approx_table_entries < 1_250 {
                (LookupTableKind::Small, 7)
            } else if approx_table_entries < 8_500 {
                (LookupTableKind::Small, 8)
            } else if approx_table_entries < 18_000 {
                (LookupTableKind::Mb, 9)
            } else {
                (LookupTableKind::Mb, 10)
            }
        }
        11 => {
            if approx_table_entries < LUT_WIDTH_11_THRESHOLD_8 {
                (LookupTableKind::Small, 8)
            } else if approx_table_entries < LUT_WIDTH_11_THRESHOLD_10 {
                (LookupTableKind::Mb, 10)
            } else {
                (LookupTableKind::Mb, 11)
            }
        }
        12 => {
            if approx_table_entries < 8_500 {
                (LookupTableKind::Small, 8)
            } else if approx_table_entries < 18_000 {
                (LookupTableKind::Mb, 9)
            } else if approx_table_entries < 60_000 {
                (LookupTableKind::Mb, 10)
            } else if approx_table_entries < 900_000 {
                (LookupTableKind::Mb, 11)
            } else {
                (LookupTableKind::Mb, 12)
            }
        }
        _ => {
            if approx_table_entries < 8_500 {
                (LookupTableKind::Small, 8)
            } else if approx_table_entries < 300_000 {
                (LookupTableKind::Mb, 11)
            } else {
                (LookupTableKind::Mb, 12)
            }
        }
    };

    if lut_kind == LookupTableKind::Small
        && (approx_table_entries >= 32_767 || max_q_off >= 32_768)
    {
        lut_kind = LookupTableKind::Na;
    }

    (lut_kind, lut_width)
}

/// Finalize task configuration with query-dependent parameters.
/// Must be called after queries are loaded to enable adaptive lookup table selection.
///
/// NCBI reference: blast_nalookup.c (BlastChooseNaLookupTable)
/// - Uses total query length (approx_table_entries) and max query offset
/// - Handles discontiguous template forcing MB lookup
pub fn finalize_task_config(
    config: &mut TaskConfig,
    total_query_length: usize,
    max_query_length: usize,
    discontig_template: bool,
) {
    let max_q_off = max_query_length.saturating_sub(1);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:45-185
    // ```c
    // lut_type = BlastChooseNaLookupTable(lookup_options,
    //                                     approx_table_entries,
    //                                     max_q_off,
    //                                     &lut_width);
    // ```
    let (lut_kind, lut_width) = choose_na_lookup_table(
        config.effective_word_size,
        total_query_length,
        max_q_off,
        discontig_template,
    );

    config.lut_word_length = lut_width;
    config.use_two_stage =
        lut_kind == LookupTableKind::Mb || config.lut_word_length < config.effective_word_size;
    config.use_direct_lookup =
        !config.use_two_stage && config.lut_word_length <= MAX_DIRECT_LOOKUP_WORD_SIZE;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:397,1334
    // ```c
    // lookup->scan_step = lookup->word_length - lookup->lut_word_length + 1;
    // mb_lt->scan_step = mb_lt->word_length - mb_lt->lut_word_length + 1;
    // ```
    config.scan_step = (config.effective_word_size - config.lut_word_length + 1).max(1);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input_aux.cpp:221-246
// ```c
// ReadSequencesToBlast(CNcbiIstream& in,
//                      bool read_proteins,
//                      const TSeqRange& range,
//                      bool parse_deflines,
//                      bool use_lcase_masking,
//                      CRef<CBlastQueryVector>& sequences,
//                      bool gaps_to_Ns)
// {
//     ...
//     sequences = input->GetAllSeqs(*scope);
//     return scope;
// }
// ```
/// Read query sequences
pub fn read_queries(args: &BlastnArgs) -> Result<(Vec<fasta::Record>, Vec<String>)> {
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

    Ok((queries, query_ids))
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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1413
// ```c
// db_length = BlastSeqSrcGetTotLen(seq_src);
// itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
// while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//        != BLAST_SEQSRC_EOF) {
// ```
pub fn subject_metadata_from_records(records: &[fasta::Record]) -> SubjectMetadata {
    let mut subject_ids = Vec::with_capacity(records.len());
    let mut db_len_total: usize = 0;
    for record in records {
        subject_ids.push(
            record
                .id()
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string(),
        );
        db_len_total += record.seq().len();
    }
    SubjectMetadata {
        db_len_total,
        db_num_seqs: records.len(),
        subject_ids,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
// ```c
// db_length = BlastSeqSrcGetTotLen(seq_src);
// itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
// while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//        != BLAST_SEQSRC_EOF) {
//     if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
//         continue;
//     }
// }
// ```
pub fn scan_subjects_metadata(args: &BlastnArgs) -> Result<SubjectMetadata> {
    let subject_reader = fasta::Reader::from_file(&args.subject)?;
    let mut subject_ids: Vec<String> = Vec::new();
    let mut db_len_total: usize = 0;
    let mut db_num_seqs: usize = 0;

    for record in subject_reader.records().filter_map(|r| r.ok()) {
        subject_ids.push(
            record
                .id()
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string(),
        );
        db_len_total += record.seq().len();
        db_num_seqs += 1;
    }

    Ok(SubjectMetadata {
        db_len_total,
        db_num_seqs,
        subject_ids,
    })
}

// NCBI reference: ncbi-blast/c++/src/objtools/readers/fasta.cpp:856-874
// ```c
// case 'a': case 'b': case 'c': case 'd':
// case 'g': case 'h':
// ...
//     char_type = eCharType_MaskedNonGap;
//     break;
// ```
// NCBI reference: ncbi-blast/c++/src/objtools/readers/fasta.cpp:1079-1089
// ```c
// void CFastaReader::x_OpenMask(void)
// {
//     m_MaskRangeStart = GetCurrentPos(ePosWithGapsAndSegs);
// }
// void CFastaReader::x_CloseMask(void)
// {
//     m_CurrentMask->SetPacked_int().AddInterval(...);
// }
// ```
pub fn collect_lowercase_masks(seq: &[u8]) -> Vec<MaskedInterval> {
    let mut masks = Vec::new();
    let mut current_start: Option<usize> = None;

    for (idx, &base) in seq.iter().enumerate() {
        let is_lower = base.is_ascii_lowercase();
        match (current_start, is_lower) {
            (None, true) => current_start = Some(idx),
            (Some(start), false) => {
                if start < idx {
                    masks.push(MaskedInterval::new(start, idx));
                }
                current_start = None;
            }
            _ => {}
        }
    }

    if let Some(start) = current_start {
        if start < seq.len() {
            masks.push(MaskedInterval::new(start, seq.len()));
        }
    }

    masks
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2547-2556
// ```c
// const bool use_lcase_masks = args.Exist(kArgUseLCaseMasking) ? ... : kDfltArgUseLCaseMasking;
// ReadSequencesToBlast(... use_lcase_masks, subjects, ...);
// ```
fn collect_lowercase_masks_for_records(
    records: &[fasta::Record],
) -> Vec<Vec<MaskedInterval>> {
    records
        .iter()
        .map(|record| collect_lowercase_masks(record.seq()))
        .collect()
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/api/dust_filter.cpp:92-128
// ```c
// const int kTopFlags = CSeq_loc::fStrand_Ignore|CSeq_loc::fMerge_All|CSeq_loc::fSort;
// if (orig_query_mask.NotEmpty()) {
//     orig_query_mask->Add(*query_masks,  kTopFlags, 0);
// } else {
//     query_masks->Merge(kTopFlags, 0);
// }
// ```
fn merge_mask_intervals(mut intervals: Vec<MaskedInterval>) -> Vec<MaskedInterval> {
    if intervals.is_empty() {
        return intervals;
    }

    intervals.sort_by_key(|interval| interval.start);
    let mut merged: Vec<MaskedInterval> = Vec::with_capacity(intervals.len());
    let mut current = intervals[0].clone();

    for interval in intervals.into_iter().skip(1) {
        if interval.start <= current.end {
            current.end = current.end.max(interval.end);
        } else {
            merged.push(current);
            current = interval;
        }
    }
    merged.push(current);
    merged
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
    query_offsets: &[i32],
    subjects: &[fasta::Record],
    subjects_packed: Option<&[Vec<u8>]>,
    approx_table_entries: usize,
) -> (LookupTables, usize) {
    eprintln!(
        "Building lookup (Task: {}, Word: {}, TwoStage: {}, LUTWord: {}, Direct: {}, DUST: {})...",
        args.task, config.effective_word_size, config.use_two_stage, config.lut_word_length, config.use_direct_lookup, args.dust
    );
    
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1312-1325
    // ```c
    // if (lookup_options->db_filter) {
    //    counts = (Uint1*)calloc(mb_lt->hashsize / 2, sizeof(Uint1));
    // }
    // if (lookup_options->db_filter) {
    //    s_FillPV(query, location, mb_lt, lookup_options);
    //    s_ScanSubjectForWordCounts(seqsrc, mb_lt, counts,
    //                               lookup_options->max_db_word_count);
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1270-1306
    // ```c
    // pv_size = (Int4)(mb_lt->hashsize >> PV_ARRAY_BTS);
    // mb_lt->pv_array_bts = ilog2(mb_lt->hashsize / pv_size);
    // ```
    let db_word_counts = if args.limit_lookup {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
        // ```c
        // if (subj_is_na) {
        //     BlastSeqBlkSetSequence(subj, sequence.data.release(),
        //        ((sentinels == eSentinels) ? sequence.length - 2 :
        //         sequence.length));
        //     ...
        //     SBlastSequence compressed_seq =
        //         subjects.GetBlastSequence(i, eBlastEncodingNcbi2na,
        //                                   eNa_strand_plus, eNoSentinels);
        //     BlastSeqBlkSetCompressedSequence(subj,
        //                               compressed_seq.data.release());
        // }
        // ```
        Some(build_db_word_counts(
            queries,
            query_masks,
            subjects,
            config.lut_word_length,
            args.max_db_word_count,
            approx_table_entries,
            subjects_packed,
        ))
    } else {
        None
    };
    let db_word_counts_ref = db_word_counts.as_deref();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1027-1034
    // ```c
    // /* Also add 1 to all indices, because lookup table indices count
    //    from 1. */
    // mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
    // mb_lt->hashtable[ecode] = index;
    // ```
    let two_stage_lookup: Option<TwoStageLookup> = if config.use_two_stage {
        Some(build_two_stage_lookup(
            queries,
            query_offsets,
            config.effective_word_size,
            config.lut_word_length,
            query_masks,
            db_word_counts_ref,
            args.max_db_word_count,
            approx_table_entries,
        ))
    } else {
        None
    };
    
    let pv_direct_lookup: Option<PvDirectLookup> = if !config.use_two_stage && config.use_direct_lookup {
        Some(build_pv_direct_lookup(
            queries,
            query_offsets,
            config.effective_word_size,
            config.effective_word_size,
            query_masks,
            db_word_counts_ref,
            args.max_db_word_count,
            approx_table_entries,
            false,
        ))
    } else {
        None
    };
    
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:548-583
    // ```c
    // Int4 BlastNaLookupTableNew(BLAST_SequenceBlk* query,
    //                            BlastSeqLoc* locations,
    //                            BlastNaLookupTable * *lut,
    //                            const LookupTableOptions * opt,
    //                            const QuerySetUpOptions* query_options,
    //                            Int4 lut_width)
    // ```
    let na_lookup: Option<NaLookupTable> = if !config.use_two_stage && !config.use_direct_lookup {
        Some(build_na_lookup(
            queries,
            query_offsets,
            config.effective_word_size,
            config.lut_word_length,
            query_masks,
            db_word_counts_ref,
            args.max_db_word_count,
        ))
    } else {
        None
    };
    
    // Use scan_step from config - if user specified a value > 0, use it; otherwise use the configured default
    // For two-stage lookup, the default was already calculated in configure_task
    let scan_step = config.scan_step;
    
    if args.verbose && config.use_two_stage {
        eprintln!("[INFO] Using two-stage lookup: lut_word_length={}, word_length={}, scan_step={}", 
                  config.lut_word_length, config.effective_word_size, scan_step);
    }
    
    (
        LookupTables {
            two_stage_lookup,
            pv_direct_lookup,
            na_lookup,
        },
        scan_step,
    )
}

/// Prepare all sequence data and configuration
pub fn prepare_sequence_data(
    args: &BlastnArgs,
    queries: Vec<fasta::Record>,
    query_ids: Vec<String>,
    subject_metadata: SubjectMetadata,
) -> SequenceData {
    let mut query_masks = apply_dust_masking(args, &queries);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2547-2556
    // ```c
    // const bool use_lcase_masks = args.Exist(kArgUseLCaseMasking) ? ... : kDfltArgUseLCaseMasking;
    // ReadSequencesToBlast(... use_lcase_masks, subjects, ...);
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/dust_filter.cpp:92-128
    // ```c
    // orig_query_mask->Add(*query_masks,  kTopFlags, 0);
    // query_masks->Merge(kTopFlags, 0);
    // ```
    if args.lcase_masking {
        let lcase_masks = collect_lowercase_masks_for_records(&queries);
        for (dust_masks, lcase) in query_masks.iter_mut().zip(lcase_masks) {
            if !lcase.is_empty() {
                dust_masks.extend(lcase);
                let merged = merge_mask_intervals(std::mem::take(dust_masks));
                *dust_masks = merged;
            }
        }
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1409
    // ```c
    // db_length = BlastSeqSrcGetTotLen(seq_src);
    // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    // ```
    let db_len_total = subject_metadata.db_len_total;
    let db_num_seqs = subject_metadata.db_num_seqs;
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-155
    // ```c
    // Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    // ```
    let subject_ids = subject_metadata.subject_ids;
    
    SequenceData {
        queries,
        query_ids,
        query_masks,
        db_len_total,
        db_num_seqs,
        subject_ids,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference: ncbi-blast/c++/src/objtools/readers/fasta.cpp:867-949
    // ```c
    // case 'a': case 'b': case 'c': case 'd': ...
    //     char_type = eCharType_MaskedNonGap;
    // ...
    // case eCharType_MaskedNonGap:
    //     OpenMask();
    // ```
    #[test]
    fn test_collect_lowercase_masks() {
        let masks = collect_lowercase_masks(b"AAaaBBbC");
        assert_eq!(
            masks,
            vec![MaskedInterval::new(2, 4), MaskedInterval::new(6, 7)]
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/dust_filter.cpp:121-127
    // ```c
    // const int kTopFlags = ... | CSeq_loc::fMerge_All | ...;
    // query_masks->Merge(kTopFlags, 0);
    // ```
    #[test]
    fn test_merge_mask_intervals() {
        let merged = merge_mask_intervals(vec![
            MaskedInterval::new(0, 2),
            MaskedInterval::new(2, 5),
            MaskedInterval::new(6, 7),
        ]);
        assert_eq!(
            merged,
            vec![MaskedInterval::new(0, 5), MaskedInterval::new(6, 7)]
        );
    }
}
