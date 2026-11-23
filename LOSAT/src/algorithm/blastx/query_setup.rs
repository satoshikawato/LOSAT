//! BLASTX context/translation/mask preparation before score/statistics setup.
use super::{args::ResolvedOptions, input::FastaRecord};
use crate::algorithm::common::masks::{combine_masks, dna_to_protein_masks, protein_to_dna_masks};
use crate::algorithm::tblastx::translation::generate_frames;
use crate::core::blast_seg::SegMasker;
use crate::utils::genetic_code::GeneticCode;
use anyhow::{bail, Result};

#[derive(Debug, Clone)]
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:86-94
pub struct QueryContext {
    pub query_index: usize,
    pub frame: i8,
    pub offset: usize,
    pub length: usize,
    pub is_valid: bool,
    pub lowercase_masks: Vec<(i32, i32)>,
    pub masks: Vec<(i32, i32)>,
    pub dna_masks: Vec<(i32, i32)>,
}
#[derive(Debug, Clone)]
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:264-266
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:499-501
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:654-659
pub struct PreparedQueryBatch {
    pub original_lengths: Vec<usize>,
    pub contexts: Vec<QueryContext>,
    pub first_context: usize,
    pub max_length: usize,
    pub min_length: usize,
    pub sequence_start: Vec<u8>,
    pub sequence_start_nomask: Vec<u8>,
    pub lookup_segments: Vec<(i32, i32)>,
    /// The values NCBI takes of the query batch (reader batch) that this unit is a part of, when
    /// the batch is searched in several units (`reader_batch.rs`); `None` when the unit is the
    /// whole batch.
    pub(crate) reader_shape: Option<std::sync::Arc<super::reader_batch::ReaderBatchShape>>,
    /// EXPERIMENT (LOSAT_X_BXUNGAPSENT): whether `sequence_start` is on the path of the seed
    /// loops' sentinel walks (`blastx/x_ungap_sent.rs`, `query_fits`, another matrix than BLOSUM62
    /// `query_fits_of`; false with the switch off).
    pub(crate) x_ungap_sent: bool,
    /// The options' search spaces that NCBI's preliminary stage sets for the contexts of the query
    /// batch (`SplitQuery_SetEffectiveSearchSpace`, the whole query's, by context index; a query
    /// chunk reads them by its own index), which the statistics of each subject read
    /// (`parameters::subject_update`); `None` where the search has none (`search::x_search_core`).
    pub(crate) subject_option_spaces: Option<std::sync::Arc<[i64]>>,
}
impl PreparedQueryBatch {
    /// EXPERIMENT (LOSAT_X_BXENDMEM): what the results of a searched batch keep: the batch without
    /// its residues and lookup segments, which the report does not read.
    pub(crate) fn x_for_results(&self) -> Self {
        if !super::kappa::x_end_mem() {
            return self.clone();
        }
        Self {
            original_lengths: self.original_lengths.clone(),
            contexts: self.contexts.clone(),
            first_context: self.first_context,
            max_length: self.max_length,
            min_length: self.min_length,
            sequence_start: Vec::new(),
            sequence_start_nomask: Vec::new(),
            lookup_segments: Vec::new(),
            reader_shape: self.reader_shape.clone(),
            x_ungap_sent: self.x_ungap_sent,
            subject_option_spaces: self.subject_option_spaces.clone(),
        }
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:123-136
pub const BATCH_SIZE: usize = 10002;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:134-139
/// The letters of RPSTBLASTN's query batch (`GetQueryBatchSize` of `eRPSTblastn`, the `default:`
/// case; DECISIONS 38 D10): its batches are read with it, never split (`rps_search::translated`).
pub const RPSTBLASTN_BATCH_SIZE: usize = 10000;
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:139-165
/// The search units of a reader batch (`records`, NCBI's query batch, searched against the
/// database `db`, `None` for `-subject`): parts that end at the first query that brings them to
/// `BATCH_SIZE` letters (`batch_ranges`), or the whole batch where a part cannot give NCBI's
/// result of the batch alone (`searches_batch_whole`). No NCBI counterpart of its own: LOSAT's
/// units are its own scheduling (`native.rs`, `XUnit`).
pub(crate) fn unit_ranges(
    records: &[FastaRecord],
    options: &ResolvedOptions,
    db: Option<&crate::seqdb::SeqDb>,
) -> Vec<std::ops::Range<usize>> {
    if !records.is_empty() && searches_batch_whole(options, db) {
        vec![0..records.len()]
    } else {
        batch_ranges(records, BATCH_SIZE)
    }
}

/// The places a unit of `records` takes (`unit_ranges`; no NCBI counterpart, scheduling only): for
/// a reader batch searched whole (`searches_batch_whole`), the parts `batch_ranges` would make of
/// it, whose letters its search holds at once (a lookup table, a preparation and the workers'
/// diagonal tables of all of them); 1 for every other unit. The batch flow keeps the places
/// searched at a time to the threads (`native.rs`, `x_flow`), and such a batch is left the threads
/// of its places (`search::x_bx_mem::whole_batch_threads`).
#[cfg(feature = "parallel")]
pub(crate) fn unit_places(
    records: &[FastaRecord],
    options: &ResolvedOptions,
    db: Option<&crate::seqdb::SeqDb>,
) -> usize {
    if !records.is_empty() && searches_batch_whole(options, db) {
        batch_ranges(records, BATCH_SIZE).len()
    } else {
        1
    }
}

/// Whether BLASTX searches each reader batch of a search of `db` (`None`: `-subject`) as one
/// unit (`unit_ranges`): with `-use_sw_tback`, and where the `-db_hard_mask` walk keeps a database
/// sequence apart (`keeps_sequences_apart`). NCBI's Smith-Waterman redo keeps state over the whole
/// query batch, which a unit cannot see alone: the status of the batch's last redone match, the
/// heap of its first query, and, in e-value mode, the heaps that a window's query context indexes
/// among the batch's queries (`runtime::kappa_batch`).
pub(crate) fn searches_batch_whole(
    options: &ResolvedOptions,
    db: Option<&crate::seqdb::SeqDb>,
) -> bool {
    options.use_sw_tback || keeps_sequences_apart(options, db)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-291
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:46-65
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:519-573
// NCBI searches each query batch as one: one lookup table, and one diagonal table per thread of
// the power of two just at least the batch's length plus the window, whose hits count their query
// offsets from the batch's start. A unit of the batch gives the batch's rows where every entry
// within reach of a hit was written by the same chunk's scan at a lower query offset of the hit's
// diagonal (the context test then never wraps and the ring's size changes nothing:
// `common::hard_chunk_walk::DiagonalWalk::apart`), which holds for every sequence of a search
// without `-db_hard_mask`. A sequence that the walk keeps apart (one on a new table, whose
// entries no sequence wrote, or one whose walk reads back to an earlier chunk's entries) meets
// entries at offsets that depend on its query's place in the batch, and an earlier chunk's hits
// of every query of the batch: such a database is searched with the batch as one unit (line NP,
// M4b-3b2; decision 57 of series M).
/// Whether the search of the database `db` (`None`: `-subject`) with `-db_hard_mask` holds a
/// sequence that the walk keeps apart (`blastinput::db_mask::ResolvedMask::hard_starts`,
/// `HardStarts::keeps_apart`: the walk that `check_losat_limits` ran, or the same walk once where
/// it did not run).
pub(crate) fn keeps_sequences_apart(
    options: &ResolvedOptions,
    db: Option<&crate::seqdb::SeqDb>,
) -> bool {
    let Some(mask) = options
        .db_mask
        .as_ref()
        .filter(|mask| mask.kind == crate::blastinput::db_mask::MaskKind::Hard)
    else {
        return false;
    };
    db.is_some_and(|db| {
        mask.hard_starts(db, options.word_size, options.window_size)
            .keeps_apart()
    })
}

pub fn batch_ranges(records: &[FastaRecord], batch_size: usize) -> Vec<std::ops::Range<usize>> {
    let mut result = Vec::new();
    let mut start = 0;
    while start < records.len() {
        let mut end = start;
        let mut length = 0;
        while end < records.len() && length < batch_size {
            length += records[end].sequence.len();
            end += 1;
        }
        result.push(start..end);
        start = end;
    }
    result
}
// No NCBI counterpart of its own: context `i` (frames 1, 2, 3, -1, -2, -3) of a query of `letters`
// nucleotides as `prepare_queries` lays the contexts out (blast_setup_cxx.cpp:191-219, cited
// there): the frame's protein length, and whether the context's strand is searched (a context
// of a strand not searched has length 0).
pub(crate) fn context_frame(letters: usize, i: usize, strand: &str) -> (usize, bool) {
    let len = letters.saturating_sub(i % 3) / 3;
    let enabled = match strand {
        "plus" => i < 3,
        "minus" => i >= 3,
        _ => true,
    };
    (len, enabled)
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:71-88
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:191-219
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:577-587
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:620-641
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:997-1003
/// The error of a query batch with a `-matrix` that NCBI's tables do not have
/// (`ResolvedOptions::unknown_matrix`): NCBI's (`unknown_matrix_message`), or LOSAT's refusal when
/// NCBI reads a file for the name (`ResolvedOptions::matrix_file`): NCBI 2.17.0 reads it where it
/// sets up the batch, and aborts after a matrix or fails on another file (`matrix_file_refusal`).
/// A run without a batch reads no file.
pub(crate) fn matrix_error_message(records: &[FastaRecord], options: &ResolvedOptions) -> String {
    match &options.matrix_file {
        Some(file) => crate::blastinput::ncbi_environment::matrix_file_refusal(
            &options.matrix,
            file,
            "BLASTX",
        ),
        None => unknown_matrix_message(records),
    }
}

/// NCBI's error of a query batch whose score block cannot load the `-matrix` of an ungapped search
/// (a name its matrix loader does not know, `args.rs` `unknown_matrix`): the message of each query
/// of the batch, a query without letters with its own after it.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:158-182 (the score block of
/// the batch: a name that the matrix loader does not know fails with -1, blast_stat.c:1679-1680,
/// posted for all queries, blast_setup.c:498-501)
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:246 (the text of -1)
/// One message per query (`TSearchMessages::ToString`), `Error: ` first (`GetMessage`); a query
/// without letters keeps its `Sequence contains no data` (blast_setup_cxx.cpp:631-638) after the
/// error, which `RemoveDuplicates` sorts first (`CSearchMessage::operator<`,
/// blast_types.hpp:295-304), as TBLASTX's (`tblastx/blast_engine/run_impl.rs`) and TBLASTN's
/// (`tblastn/stage_d_pipeline.rs`): the text of `common::unknown_matrix::unknown_matrix_message_of`.
pub(crate) fn unknown_matrix_message(records: &[FastaRecord]) -> String {
    crate::algorithm::common::unknown_matrix::unknown_matrix_message_of(records.iter().map(
        |record| {
            if record.sequence.is_empty() {
                vec![crate::blastinput::query_batch::NO_DATA_MESSAGE.to_string()]
            } else {
                Vec::new()
            }
        },
    ))
}

pub fn prepare_queries(
    records: &[FastaRecord],
    options: &ResolvedOptions,
) -> Result<PreparedQueryBatch> {
    prepare_unit_queries(records, None, options)
}

/// `prepare_queries` for `records`, one unit of the reader batch of `batch` when it is given
/// (`reader_batch.rs`): whether the setup fails for want of data is the reader batch's.
pub(crate) fn prepare_unit_queries(
    records: &[FastaRecord],
    batch: Option<&super::reader_batch::ReaderBatchShape>,
    options: &ResolvedOptions,
) -> Result<PreparedQueryBatch> {
    if records.is_empty() {
        bail!("cannot prepare an empty BLASTX query batch");
    }
    crate::blastinput::value_parsers::genetic_code(&options.query_gencode.to_string())
        .map_err(anyhow::Error::msg)?;
    let code = GeneticCode::try_from_id(options.query_gencode).map_err(anyhow::Error::msg)?;
    // Diagnostic (LOSAT_TIMING): the contexts' layout and the six-frame translation for the
    // common line (utils/xtiming.rs). No NCBI counterpart.
    let t_translate = crate::utils::xtiming::clock();
    let (mut contexts, mut max_length, mut min_length) = (Vec::new(), 0, usize::MAX);
    let mut offset = 0;
    for (q, record) in records.iter().enumerate() {
        for i in 0..6 {
            let frame = if i < 3 {
                (i + 1) as i8
            } else {
                -((i - 2) as i8)
            };
            let (len, enabled) = context_frame(record.sequence.len(), i, &options.strand);
            max_length = max_length.max(len);
            min_length = min_length.min(len);
            let length = if enabled { len } else { 0 };
            let lowercase_masks = if enabled {
                dna_to_protein_masks(&record.lowercase_masks, record.sequence.len(), frame)
            } else {
                Vec::new()
            };
            contexts.push(QueryContext {
                query_index: q,
                frame,
                offset,
                length,
                is_valid: length > 0,
                lowercase_masks,
                masks: Vec::new(),
                dna_masks: Vec::new(),
            });
            if length > 0 {
                offset += length + 1;
            }
        }
    }
    if options.strand != "both" {
        min_length = 0;
    }
    let last = contexts.last().unwrap();
    let buffer_length = last.offset + last.length + if last.length > 0 { 2 } else { 1 };
    let mut buffer = vec![0; buffer_length];
    for (q, record) in records.iter().enumerate() {
        for frame in generate_frames(&record.sequence, &code) {
            let i = if frame.frame > 0 {
                frame.frame as usize - 1
            } else {
                3 + (-frame.frame) as usize - 1
            };
            let context = &contexts[q * 6 + i];
            if context.length > 0 {
                buffer[context.offset..context.offset + context.length + 2]
                    .copy_from_slice(&frame.aa_seq);
            }
        }
    }
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::Translate, t_translate);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup.hpp:192-197
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:634-639
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:650-652
    // `qinfo` is the query batch's: a unit of a reader batch (`batch`) fails only when no
    // context of the reader batch is valid, with the messages of all its empty queries.
    let (fails, empty) = match batch {
        Some(batch) => (batch.set_up_fails(), batch.empty_queries),
        None => {
            let empty = records.iter().filter(|q| q.sequence.is_empty()).count();
            (!contexts.iter().any(|c| c.is_valid) && empty > 0, empty)
        }
    };
    if fails {
        anyhow::bail!("{}", "Warning: Sequence contains no data ".repeat(empty));
    }
    // The matrix error comes after the set-up's errors above, with the messages of the whole
    // query batch: a unit of a reader batch raises its batch's (`ReaderBatchShape::matrix_error`).
    if options.unknown_matrix {
        let message = match batch {
            Some(batch) => batch.matrix_error.as_deref().map(str::to_string),
            None => None,
        };
        anyhow::bail!(
            "{}",
            message.unwrap_or_else(|| matrix_error_message(records, options))
        );
    }
    let nomask = buffer.clone();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1144-1158
    // Same SEG parameters as the block above; the masker is built once for all contexts instead of
    // once per context (placement only).
    let seg_masker = options.seg.filter_options().map(|(window, locut, hicut)| {
        SegMasker::new(
            if window > 0 { window as usize } else { 12 },
            if locut > 0.0 { locut } else { 2.2 },
            if hicut > 0.0 { hicut } else { 2.5 },
        )
        // LOSAT_X_SEGTRIMMEMO (graduated 2026-10-26): the chunks of a split query (split.rs) are masked
        // again after the whole query, mostly over the same segments, so this masker takes the
        // trim of a segment trimmed before (on any thread) from the memo.
        .x_with_trim_memo()
    });
    // SEG reads only the residues of its own context and the hard masking below writes only there,
    // so the contexts of a large batch (the six frames of the full query) are filtered on the
    // search pool before the loop consumes the intervals in context order.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1280-1290
    // NCBI filters the contexts one after another in this loop. SEG of one context reads only that
    // context, so with a pool of 2+ threads the contexts of a large batch are filtered in parallel
    // and the intervals are kept in seg_ready; the loop below still consumes them in context order.
    #[allow(unused_mut)]
    let mut seg_ready: Vec<Option<Vec<(i32, i32)>>> = Vec::new();
    #[cfg(feature = "parallel")]
    if let Some(masker) = seg_masker.as_ref() {
        if super::runtime::x_bx_par_rule(|| {
            buffer.len() >= 60_000
                && !super::runtime::x_inner_serial()
                && rayon::current_thread_index().is_some()
                && rayon::current_num_threads() > 1
        }) {
            use rayon::prelude::*;
            let buffer = &buffer;
            // Diagnostic (LOSAT_TIMING): the parallel SEG's wall time for the common line
            // (utils/xtiming.rs). No NCBI counterpart.
            let t_seg = crate::utils::xtiming::clock();
            seg_ready = contexts
                .par_iter()
                .map(|context| {
                    context.is_valid.then(|| {
                        masker
                            .mask_sequence(
                                &buffer[context.offset + 1..context.offset + 1 + context.length],
                            )
                            .into_iter()
                            .map(|r| (r.start as i32, r.end as i32 - 1))
                            .collect()
                    })
                })
                .collect();
            crate::utils::xtiming::add(crate::utils::xtiming::Phase::Seg, t_seg);
        }
    }
    for (context_index, context) in contexts.iter_mut().enumerate() {
        if !context.is_valid {
            continue;
        }
        let mut masks = Vec::new();
        if let Some(masker) = seg_masker.as_ref() {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1280-1290
            // Dispatch point: a context filtered ahead of the loop uses its stored intervals; any
            // other context is filtered here exactly as before.
            if let Some(ready) = seg_ready.get_mut(context_index).and_then(Option::take) {
                masks.extend(ready);
            } else {
                // Diagnostic (LOSAT_TIMING): this context's SEG for the common line
                // (utils/xtiming.rs). No NCBI counterpart.
                let t_seg = crate::utils::xtiming::clock();
                masks.extend(
                    masker
                        .mask_sequence(
                            &buffer[context.offset + 1..context.offset + 1 + context.length],
                        )
                        .into_iter()
                        .map(|r| (r.start as i32, r.end as i32 - 1)),
                );
                crate::utils::xtiming::add(crate::utils::xtiming::Phase::Seg, t_seg);
            }
        }
        // SEG's intervals first, as `BlastSeqLocAppend` of `lcase_mask_slp` to `filter_out`: the sort of
        // `BlastSeqLocCombine` compares only the left ends, and a one-letter run on the first letter
        // of a SEG interval is merged only in this order (the results' masks).
        masks.extend_from_slice(&context.lowercase_masks);
        context.masks = combine_masks(masks);
        if !options.soft_masking {
            for &(from, to) in &context.masks {
                for index in from..=to {
                    buffer[context.offset + 1 + index as usize] = 21;
                }
            }
        }
        context.dna_masks = protein_to_dna_masks(
            &context.masks,
            core_dna_length(&contexts_lengths(
                records[context.query_index].sequence.len(),
            )),
            context.frame,
        );
    }
    let lookup_segments = complement_masks(&contexts);
    // The sentinel walks take BLOSUM62's keys, or another matrix's where it has them (item SC6b,
    // `x_ungap_sent.rs` `query_fits_of`).
    let x_ungap_sent = if options.scoring == crate::config::ScoringMatrix::Blosum62 {
        super::seed::x_ungap_sent::query_fits(&buffer, options.word_size, options.window_size)
    } else {
        super::seed::x_ungap_sent::query_fits_of(
            options.scoring,
            &buffer,
            options.word_size,
            options.window_size,
        )
    };
    Ok(PreparedQueryBatch {
        original_lengths: records.iter().map(|r| r.sequence.len()).collect(),
        contexts,
        first_context: if options.strand == "minus" { 3 } else { 0 },
        max_length,
        min_length,
        sequence_start: buffer,
        sequence_start_nomask: nomask,
        lookup_segments,
        reader_shape: None,
        x_ungap_sent,
        subject_option_spaces: None,
    })
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:923-929
fn contexts_lengths(length: usize) -> [usize; 3] {
    [
        length / 3,
        length.saturating_sub(1) / 3,
        length.saturating_sub(2) / 3,
    ]
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:143-160
fn core_dna_length(lengths: &[usize; 3]) -> usize {
    lengths.iter().sum::<usize>() + 2
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1073-1115
fn complement_masks(contexts: &[QueryContext]) -> Vec<(i32, i32)> {
    let mut result = Vec::new();
    for c in contexts.iter().filter(|c| c.is_valid) {
        let start = c.offset as i32;
        let end = start + c.length as i32 - 1;
        if c.masks.is_empty() {
            result.push((start, end));
            continue;
        }
        let mut first = true;
        let mut last_interval_open = true;
        let mut left = 0;
        for &(a, b) in &c.masks {
            let filter_start = start + a;
            let filter_end = start + b;
            if first {
                last_interval_open = true;
                first = false;
                if filter_start > start {
                    left = start;
                } else {
                    left = filter_end + 1;
                    continue;
                }
            }
            result.push((left, filter_start - 1));
            if filter_end >= end {
                last_interval_open = false;
                break;
            } else {
                left = filter_end + 1;
            }
        }
        if last_interval_open {
            result.push((left, end));
        }
    }
    result
}

#[cfg(test)]
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:86-94
#[path = "query_setup_tests.rs"]
mod tests;
