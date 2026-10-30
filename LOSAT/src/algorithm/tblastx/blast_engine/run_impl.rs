//! Main run() function for backbone mode
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c

use super::*;
use crate::algorithm::tblastx::database::DatabaseTitles;
use crate::algorithm::tblastx::hit_collector::HitCollector;
use crate::algorithm::tblastx::subjects::{
    SearchSpaceOverflow, SubjectSetNumbers, Subjects, SubjectsUnreadable,
};
use crate::blastinput::report_counts::ReportCounts;
use std::sync::Arc;

use crate::algorithm::blastn::blast_engine::{QueryReading, NO_DATA_MESSAGE};
use crate::blastinput::fasta_reader::{
    read_queries, read_subjects, FastaInputSource, FastaRecord, QueryEnd, QueryRecords,
    ReaderConfig,
};
use crate::blastinput::seq_range::{
    cut_queries, cut_subjects, next_ranged_batch, parse_optional_range, subjects_read, Placements,
    QueryInput, RangeRole, SequenceRange,
};

// EXPERIMENT (PAIRPAR, graduated on 2026-10-22: no variable): the seed stage of a subject as units
// of (subject frame, chunk, query context) on the search pool; a child module of this file so that
// it uses the chunk helpers below unchanged (see tblastx/x_pair_par.rs).
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-516
// ```c
// while (scan_range[1] <= scan_range[2]) {
//     /* scan the subject sequence for hits */
//     hits = scansub(lookup_wrap, subject,
//                               offset_pairs, array_size, scan_range);
// ...
//     for (i = 0; i < hits; ++i) {
// ```
// The module ports the per-hit body of this loop and runs it per query context; it is used
// where its rule (`x_pair_par::decide`) takes the batch.
#[path = "../x_pair_par.rs"]
mod x_pair_par;
pub use x_pair_par::print_summary as x_pair_par_print_summary;

// The seed-stage tables of a query batch are held by two
// per-batch states (the lookup table; the unit path with its group tables) and reached only
// through one stage of each per subject; they are freed when the batch's
// last stage of that kind ends, before the last subject's HSP list is merged and linked (see
// tblastx/x_seed_release.rs). A child module of this file, as x_pair_par is.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1604-1605
// ```c
// word_params = BlastInitialWordParametersFree(word_params);
// s_BlastCoreAuxStructFree(aux_struct);
// ```
// NCBI frees its diagonal table after the subject loop and its lookup table with the search's
// internal data; when LOSAT frees its copies is not a value NCBI computes.
#[path = "../x_seed_release.rs"]
mod x_seed_release;
pub use x_seed_release::print_summary as x_seed_release_print_summary;

// Diagnostic (LOSAT_TIMING): the stage times of a batch search and of its subjects, two
// `[TIMING]` lines after the ones this file prints (see tblastx/x_search_timing.rs). A child
// module of this file, as x_seed_release is.
// No NCBI counterpart: a timer; it computes no value NCBI computes and none the search reads.
#[path = "../x_search_timing.rs"]
mod x_search_timing;
// EXPERIMENT (LOSAT_X_TBXQFRAMES): the SEG of a query batch's frames on the search pool, one task
// per frame (see tblastx/x_query_frames.rs). A child module of this file, as x_seed_release is.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1280-1281
// ```c
//     for (context = query_info->first_context;
//          context <= query_info->last_context; ++context) {
// ```
// NCBI filters the contexts one after the other; each context's filtering reads and writes only
// that context, so where it runs is scheduling only.
#[path = "../x_query_frames.rs"]
mod x_query_frames;
// EXPERIMENT (LOSAT_X_TBXNEXT): the subjects of the parallel subject loop handed to the threads one
// at a time, each thread taking the next subject when it is free, with one worker state per
// thread (see tblastx/x_subject_next.rs). A child module of this file, as x_seed_release is.
// No NCBI counterpart: scheduling only; NCBI's threads take the next subjects from the database
// when they are free (blast_engine.c:1411-1413), which is what the module does.
#[path = "../x_subject_next.rs"]
#[cfg_attr(
    not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )),
    allow(dead_code)
)]
mod x_subject_next;
// EXPERIMENT (LOSAT_X_TBXPIPE): the next query batch starts when the batch before has ended its
// seed stage, while that batch still links its subjects; the batch loop takes the results in its
// order (see tblastx/x_batch_pipe.rs). A child module of this file, as x_subject_next is.
// No NCBI counterpart: scheduling only; NCBI searches its query batches one after the other
// (tblastx_app.cpp:176-195), and each batch's result is the same whenever it is searched.
#[path = "../x_batch_pipe.rs"]
#[cfg_attr(
    not(all(feature = "parallel", not(target_arch = "wasm32"))),
    allow(dead_code)
)]
mod x_batch_pipe;
// The reports of each query batch written when the batch loop takes it (line DBID, ID4t fix 2b;
// see tblastx/report_stream.rs). A child module of this file, as x_batch_pipe is.
#[path = "../report_stream.rs"]
mod report_stream;
use x_search_timing::XStage;
// Diagnostic (LOSAT_TIMING): the stage times of the two linking calls, two more `[TIMING]`
// lines after those (see tblastx/x_link_timing.rs).
use crate::algorithm::tblastx::x_link_timing::{LinkCall, XLinkTimer};

/// The queries of a search: the input records, the searched queries (each cut to its
/// `-query_loc` interval), the query input (where each searched query lies in its record,
/// its input position, the skipped records), and where each subject's searched letters
/// lie in its record.
pub(crate) struct QueryRuns<'a> {
    pub input_records: &'a [FastaRecord],
    pub records: &'a [FastaRecord],
    pub input: &'a QueryInput,
    pub subject_placements: &'a Placements,
}

/// The subjects of the batches of a search, the numbers its statistics read, and what its
/// reports show of them.
struct RunSubjects<'a> {
    subjects: Subjects<'a>,
    numbers: SubjectSetNumbers,
    shown: ShownSubjects<'a>,
}

/// The subject set as the reports show it.
enum ShownSubjects<'a> {
    /// `-subject`: `User specified sequence set (Input: PATH)`; the records are the whole set.
    /// EXPERIMENT (LOSAT_X_TBXSTORE): with `x_stored`, the records hold no letters; the stores
    /// hold them, by record (the ordinals `0..n`, which the reports' `SubjectLetters::Stored`
    /// reads).
    Fasta {
        records: &'a [FastaRecord],
        x_stored: Option<(
            &'a [crate::sequence::packed_subject::PackedNucleotide],
            &'a [u32],
        )>,
    },
    /// `-db`: the `-db` value (outfmt 7's `# Database:` line), the database's `Database:` lines
    /// (outfmt 0), and the sequences by ordinal id, whose reported ones the reports show by
    /// rank (`database::reported_subjects`).
    Database {
        database: &'a super::super::database::TblastxDatabase,
        pairwise: crate::report::pairwise::PairwiseDatabase,
        stores: &'a [crate::sequence::packed_subject::PackedNucleotide],
    },
}

/// Where the subjects of a search come from.
enum SubjectSource<'a> {
    /// The records of `-subject`, and where `-subject_loc` cut them.
    Fasta {
        records: &'a [FastaRecord],
        placements: &'a Placements,
        /// EXPERIMENT (LOSAT_X_TBXSTORE): the records' letters, packed as read (the records
        /// hold none).
        x_stores: Option<&'a [crate::sequence::packed_subject::PackedNucleotide]>,
    },
    /// A `-db` database (`tblastx/database.rs`).
    Database(&'a super::super::database::TblastxDatabase),
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:54
// ```c
// #define MAX_DBSEQ_LEN 5000000
// ```
const MAX_DBSEQ_LEN: usize = 5_000_000;

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_gapalign.h:54
// ```c
// #define MAX_DBSEQ_LEN 5000000
// ```
#[inline]
fn tblastx_max_dbseq_len_for_run() -> usize {
    #[cfg(debug_assertions)]
    {
        if let Ok(value) = std::env::var("LOSAT_TBLASTX_TEST_CHUNK_SIZE") {
            if let Ok(parsed) = value.parse::<usize>() {
                return parsed.max(DBSEQ_CHUNK_OVERLAP + 3);
            }
        }
    }
    MAX_DBSEQ_LEN
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:492-505
// ```c
// while (scan_range[1] <= scan_range[2]) {
//     hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
// ```
#[inline]
fn tblastx_scan_chunk_size_for_run(search_unit_len: usize, num_threads: usize) -> usize {
    #[cfg(debug_assertions)]
    {
        if let Ok(value) = std::env::var("LOSAT_TBLASTX_TEST_SCAN_CHUNK_SIZE") {
            if let Ok(parsed) = value.parse::<usize>() {
                return parsed.max(1);
            }
        }
    }

    if num_threads <= 1 {
        search_unit_len.max(1)
    } else {
        search_unit_len.div_ceil(num_threads).max(1)
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:221-251
// ```c
// if (backup->offset + MAX_DBSEQ_LEN <
//     backup->hard_ranges[backup->hm_index].right) {
//     subject->length = MAX_DBSEQ_LEN;
//     backup->next = backup->offset + MAX_DBSEQ_LEN - dbseq_chunk_overlap;
// } else {
//     subject->length = backup->hard_ranges[backup->hm_index].right
//                     - backup->offset;
// }
// ```
fn tblastx_estimated_aa_chunks(aa_len: usize, max_dbseq_len: usize) -> usize {
    if aa_len == 0 {
        return 0;
    }
    if aa_len <= max_dbseq_len {
        return 1;
    }
    let stride = max_dbseq_len.saturating_sub(DBSEQ_CHUNK_OVERLAP).max(1);
    1 + aa_len.saturating_sub(max_dbseq_len).div_ceil(stride)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:1296-1308
// ```c
// for (context = 0; context < num_frames; ++context) {
//    int frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
//    BLAST_GetTranslation(subject_blk->sequence_start, nucl_seq_rev,
//       subject_blk->length, frame, retval->translations[context], gen_code_string);
// }
// ```
fn tblastx_estimated_subject_chunk_work_items(
    subject_nucl_len: usize,
    max_dbseq_len: usize,
) -> usize {
    let mut max_frame_chunks = 0usize;
    for frame_offset in 0..3usize {
        let aa_len = subject_nucl_len.saturating_sub(frame_offset) / 3;
        let frame_chunks = tblastx_estimated_aa_chunks(aa_len, max_dbseq_len);
        max_frame_chunks = max_frame_chunks.max(frame_chunks);
    }
    max_frame_chunks
}
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:75-127
// ```c
// while (s_DetermineScanningOffsets(subject, word_length, word_length, s_range)) {
//     ...
//     for (s = s_first; s <= s_last; s++) {
// ```
fn tblastx_scan_interiors(
    search_unit_len: usize,
    scan_chunk_size: Option<usize>,
) -> Vec<(usize, usize)> {
    let chunk_size = scan_chunk_size.unwrap_or(search_unit_len).max(1);
    let mut interiors = Vec::new();
    let mut start = 0usize;
    while start < search_unit_len {
        let end = start.saturating_add(chunk_size).min(search_unit_len);
        interiors.push((start, end));
        start = end;
    }
    interiors
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:221-310
// ```c
// subject->seq_ranges = subject->seq_ranges_allocated;
// subject->num_seq_ranges = 0;
// ...
// subject->seq_ranges[subject->num_seq_ranges].left = MAX(...);
// subject->seq_ranges[subject->num_seq_ranges].right = MIN(...);
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:83-123
// ```c
// for (s = s_first; s <= s_last; s++) {
//     ...
//     offset_pairs[i + totalhits].qs_offsets.s_off = s_off;
// }
// ```
fn clip_tblastx_seq_ranges_for_scan_interior(
    seq_ranges: &[(i32, i32)],
    interior_start: usize,
    interior_end: usize,
    wordsize: usize,
    subject_len: usize,
) -> Vec<(i32, i32)> {
    if interior_start >= interior_end || subject_len < wordsize {
        return Vec::new();
    }

    let emit_left = interior_start.min(subject_len);
    let emit_right = interior_end.min(subject_len);
    let read_right = emit_right
        .saturating_add(wordsize.saturating_sub(1))
        .min(subject_len);
    let mut clipped = Vec::with_capacity(seq_ranges.len());
    for &(left, right) in seq_ranges {
        let range_left = (left.max(emit_left as i32)) as usize;
        let range_right = (right.min(read_right as i32)) as usize;
        if range_right.saturating_sub(range_left) >= wordsize {
            clipped.push((range_left as i32, range_right as i32));
        }
    }
    clipped
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:190-192
// ```c
// /** Size of overlap in splitting query or database sequence */
// #define DBSEQ_CHUNK_OVERLAP 100
// ```
const DBSEQ_CHUNK_OVERLAP: usize = 100;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1534-1539
// ```c
// /** Maximal diagonal distance between HSP starting offsets, within which HSPs
//  * from search of different chunks of subject sequence are considered for
//  * merging.
//  */
// #define OVERLAP_DIAG_CLOSE 10
// ```
const OVERLAP_DIAG_CLOSE: i32 = 10;

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:122-143
// ```c
// typedef struct SubjectSplitStruct {
//    Uint1* sequence;
//    SSeqRange  full_range;
//    SSeqRange* hard_ranges;
//    Int4 num_hard_ranges;
//    Int4 hm_index;
//    Int4 offset;
//    Int4 next;
// } SubjectSplitStruct;
// ```
struct SubjectSplitState {
    full_right: i32,
    hard_ranges: [(i32, i32); 1],
    hm_index: usize,
    offset: i32,
    next: i32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct SubjectChunk {
    offset: usize,
    length: usize,
    overlap: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SubjectChunkStatus {
    Done,
    Ok(SubjectChunk),
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:452-584
// ```c
// BlastInitHitListReset(init_hitlist);
// aux_struct->WordFinder(..., init_hitlist, ...);
// BLAST_GetUngappedHSPList(..., &hsp_list);
// Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
// status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list, ...);
// ```
#[derive(Default)]
struct TblastxChunkScanStats {
    hsp_saved: usize,
    hsp_filtered_by_cutoff: usize,
    score_distribution: Vec<i32>,
}

struct TblastxChunkScanResult {
    chunk: SubjectChunk,
    hits: Vec<UngappedHit>,
    stats: TblastxChunkScanStats,
}

impl SubjectSplitState {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:146-184
    // ```c
    // backup->full_range.left = 0;
    // backup->full_range.right = subject->length;
    // backup->hard_ranges = &(backup->full_range);
    // backup->num_hard_ranges = 1;
    // backup->hm_index = 0;
    // backup->offset = backup->hard_ranges[0].left;
    // backup->next = backup->offset;
    // ```
    fn new(subject_len: usize) -> Self {
        let full_right = subject_len as i32;
        Self {
            full_right,
            hard_ranges: [(0, full_right)],
            hm_index: 0,
            offset: 0,
            next: 0,
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:221-310
    // ```c
    // if (backup->next >= backup->full_range.right) return SUBJECT_SPLIT_DONE;
    // residual = is_nucleotide ?  backup->next % COMPRESSION_RATIO : 0;
    // backup->offset = backup->next - residual;
    // if (backup->offset + MAX_DBSEQ_LEN <
    //     backup->hard_ranges[backup->hm_index].right) {
    //     subject->length = MAX_DBSEQ_LEN;
    //     backup->next = backup->offset + MAX_DBSEQ_LEN - dbseq_chunk_overlap;
    // } else {
    //     subject->length = backup->hard_ranges[backup->hm_index].right
    //                     - backup->offset;
    //     backup->hm_index++;
    //     backup->next = (backup->hm_index < backup->num_hard_ranges) ?
    //                     backup->hard_ranges[backup->hm_index].left :
    //                     backup->full_range.right;
    // }
    // ```
    fn next_chunk(&mut self, max_dbseq_len: usize, chunk_overlap: usize) -> SubjectChunkStatus {
        if self.next >= self.full_right {
            return SubjectChunkStatus::Done;
        }

        let residual = 0;
        self.offset = self.next - residual;
        let offset = self.offset;
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:580-584
        // ```c
        // Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
        // overlap = (backup.offset == backup.hard_ranges[backup.hm_index].left) ?
        //           0 : dbseq_chunk_overlap;
        // ```
        //
        // Preserve the current hard-range start before `hm_index` advances in the
        // final-chunk branch below; the first chunk in a hard range has no prior
        // chunk to overlap.
        let hard_left = self.hard_ranges[self.hm_index].0;
        let hard_right = self.hard_ranges[self.hm_index].1;
        let length = if offset + (max_dbseq_len as i32) < hard_right {
            self.next = offset + max_dbseq_len as i32 - chunk_overlap as i32;
            max_dbseq_len
        } else {
            let length = (hard_right - offset).max(0) as usize;
            self.hm_index = self.hm_index.saturating_add(1);
            self.next = if self.hm_index < self.hard_ranges.len() {
                self.hard_ranges[self.hm_index].0
            } else {
                self.full_right
            };
            length
        };

        let overlap = if offset == hard_left {
            0
        } else {
            chunk_overlap
        };

        SubjectChunkStatus::Ok(SubjectChunk {
            offset: offset as usize,
            length,
            overlap,
        })
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:162-173
// ```c
// Blast_ExtendWordExit(Blast_ExtendWord * ewp, Int4 subject_length)
// {
//    if (ewp->diag_table->offset >= INT4_MAX / 4) {
//       ewp->diag_table->offset = ewp->diag_table->window;
//       s_BlastDiagClear(ewp->diag_table);
//    } else {
//       ewp->diag_table->offset += subject_length + ewp->diag_table->window;
//    }
// }
// ```
#[inline]
fn advance_tblastx_diag_offset(
    diag_offset: &mut i32,
    diag_array: &mut [DiagStruct],
    window: i32,
    subject_len: usize,
) {
    if *diag_offset >= i32::MAX / 4 {
        *diag_offset = window;
        for d in diag_array.iter_mut() {
            *d = DiagStruct::clear(window);
        }
    } else {
        *diag_offset += subject_len as i32 + window;
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits_priv.h:66-72
// ```c
// #define CONTAINED_IN_HSP(a,b,c,d,e,f) \
//     (((a <= c && b >= c) && (d <= f && e >= f)) ? TRUE : FALSE)
// ```
#[inline]
fn contained_in_hsp(a: usize, b: usize, c: usize, d: usize, e: usize, f: usize) -> bool {
    a <= c && b >= c && d <= f && e >= f
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1488-1533
// ```c
// static Boolean
// s_BlastMergeTwoHSPs(BlastHSP* hsp1, BlastHSP* hsp2, Boolean allow_gap)
// {
//    if (!allow_gap &&
//        hsp1->subject.offset - hsp2->subject.offset
//        - hsp1->query.offset + hsp2->query.offset) return FALSE;
//    if(hsp1->subject.frame != hsp2->subject.frame) return FALSE;
//    if (CONTAINED_IN_HSP(...) || CONTAINED_IN_HSP(...)) {
//       ...
//       return TRUE;
//    }
//    return FALSE;
// }
// ```
fn merge_two_tblastx_chunk_hsps(hsp1: &mut UngappedHit, hsp2: &UngappedHit) -> bool {
    if hsp1.s_aa_start as isize - hsp2.s_aa_start as isize - hsp1.q_aa_start as isize
        + hsp2.q_aa_start as isize
        != 0
    {
        return false;
    }
    if hsp1.s_frame != hsp2.s_frame {
        return false;
    }

    if contained_in_hsp(
        hsp1.q_aa_start as usize,
        hsp1.q_aa_end as usize,
        hsp2.q_aa_start as usize,
        hsp1.s_aa_start as usize,
        hsp1.s_aa_end as usize,
        hsp2.s_aa_start as usize,
    ) || contained_in_hsp(
        hsp1.q_aa_start as usize,
        hsp1.q_aa_end as usize,
        hsp2.q_aa_end as usize,
        hsp1.s_aa_start as usize,
        hsp1.s_aa_end as usize,
        hsp2.s_aa_end as usize,
    ) {
        let len1 = (hsp1.q_aa_end as usize).saturating_sub(hsp1.q_aa_start as usize) as f64;
        let len2 = (hsp2.q_aa_end as usize).saturating_sub(hsp2.q_aa_start as usize) as f64;
        let score_density = (hsp1.raw_score as f64 + hsp2.raw_score as f64) / (len1 + len2);

        hsp1.q_aa_start = hsp1.q_aa_start.min(hsp2.q_aa_start);
        hsp1.s_aa_start = hsp1.s_aa_start.min(hsp2.s_aa_start);
        hsp1.q_aa_end = hsp1.q_aa_end.max(hsp2.q_aa_end);
        hsp1.s_aa_end = hsp1.s_aa_end.max(hsp2.s_aa_end);

        if hsp2.raw_score > hsp1.raw_score {
            hsp1.q_seed_off = hsp2.q_seed_off;
            hsp1.s_seed_off = hsp2.s_seed_off;
            hsp1.raw_score = hsp2.raw_score;
        }

        let new_len = (hsp1.q_aa_end as usize).saturating_sub(hsp1.q_aa_start as usize) as f64;
        hsp1.raw_score = hsp1.raw_score.max((score_density * new_len) as i32);
        return true;
    }

    false
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2857-2995
// ```c
// if (contexts_per_query < 0) {      /* subject seq is split */
//    if (hsp1->subject.end > split_offsets[0]) { ... }
//    if (hsp2->subject.offset < split_offsets[0] + chunk_overlap_size) { ... }
// }
// ...
// if (!hsp2 || hsp1->context != hsp2->context) continue;
// end_diag = s_HSPEndDiag(hsp1);
// start_diag = s_HSPStartDiag(hsp2);
// if (ABS(end_diag - start_diag) < OVERLAP_DIAG_CLOSE) {
//    if (s_BlastMergeTwoHSPs(hsp1, hsp2, allow_gap)) { ... }
// }
// ```
fn merge_tblastx_subject_chunk_hits(
    combined: &mut Vec<UngappedHit>,
    mut incoming: Vec<UngappedHit>,
    split_offset: usize,
    chunk_overlap_size: usize,
) {
    if incoming.is_empty() {
        return;
    }
    if combined.is_empty() {
        combined.append(&mut incoming);
        return;
    }

    let mut combined_overlap = Vec::new();
    for (idx, hsp) in combined.iter().enumerate() {
        if hsp.s_aa_end as usize > split_offset {
            combined_overlap.push(idx);
        }
    }

    let mut incoming_overlap = Vec::new();
    for (idx, hsp) in incoming.iter().enumerate() {
        if (hsp.s_aa_start as usize) < split_offset.saturating_add(chunk_overlap_size) {
            incoming_overlap.push(idx);
        }
    }

    let mut deleted = vec![false; incoming.len()];
    for &i in combined_overlap.iter() {
        let hsp1_context = combined[i].ctx_idx;
        for &j in incoming_overlap.iter() {
            if deleted[j] || incoming[j].ctx_idx != hsp1_context {
                continue;
            }
            let end_diag = combined[i].q_aa_end as i32 - combined[i].s_aa_end as i32;
            let start_diag = incoming[j].q_aa_start as i32 - incoming[j].s_aa_start as i32;
            if (end_diag - start_diag).abs() < OVERLAP_DIAG_CLOSE
                && merge_two_tblastx_chunk_hsps(&mut combined[i], &incoming[j])
            {
                deleted[j] = true;
            }
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3000-3035
    // ```c
    // Blast_HSPListPurgeNullHSPs(hsp_list);
    // ...
    // s_BlastHSPListsCombineByScore(hsp_list, combined_hsp_list, new_hspcnt);
    // hsp_list = Blast_HSPListFree(hsp_list);
    // ```
    let mut delete_idx = 0usize;
    incoming.retain(|_| {
        let keep = !deleted[delete_idx];
        delete_idx += 1;
        keep
    });
    combined.append(&mut incoming);
    if !ungapped_hits_is_sorted_by_score_ncbi(combined) {
        sort_ungapped_hits_by_score_ncbi(combined);
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3038-3051
// ```c
// void Blast_HSPListAdjustOffsets(BlastHSPList* hsp_list, Int4 offset)
// {
//    if (offset == 0) return;
//    for (index=0; index<hsp_list->hspcnt; index++) {
//       hsp->subject.offset += offset;
//       hsp->subject.end += offset;
//       hsp->subject.gapped_start += offset;
//    }
// }
// ```
fn adjust_tblastx_chunk_subject_offsets(hits: &mut [UngappedHit], offset: usize) {
    if offset == 0 {
        return;
    }
    for hit in hits {
        hit.s_aa_start = u32::try_from(hit.s_aa_start as usize + offset)
            .expect("s_aa_start: a subject offset fits u32 (NCBI's BlastSeg offset is Int4)");
        hit.s_aa_end = u32::try_from(hit.s_aa_end as usize + offset)
            .expect("s_aa_end: a subject offset fits u32 (NCBI's BlastSeg end is Int4)");
        hit.s_seed_off = u32::try_from(hit.s_seed_off as usize + offset)
            .expect("s_seed_off: a subject offset fits u32 (NCBI's gapped_start is Int4)");
    }
}

/// `settings` are the NCBI application settings that LOSAT reproduces
/// (`ncbi_environment::check_ncbi_application_settings`): the input readers use them once
/// the program reads with `fasta_reader` (steps S3-S8 of the port plan).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:67-72
/// ```c++
///     CNcbiApplication* app = CNcbiApplication::Instance();
///     if (app) {
///         const CNcbiRegistry& registry = app->GetConfig();
///         x_LoadDataLoadersConfig(registry);
///         x_LoadBlastDbDataLoaderConfig(registry);
///     }
/// ```
pub fn run(
    mut args: TblastxArgs,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:106-111
    // ```c
    // if(RecoverSearchStrategy(args, m_CmdLineArgs)) {
    // 	opts_hndl.Reset(&*m_CmdLineArgs->SetOptionsForSavedStrategy(args));
    // }
    // else {
    // 	opts_hndl.Reset(&*m_CmdLineArgs->SetOptions(args));
    // }
    // ```
    // NCBI's option handlers read the subjects first, then open the query and the output,
    // then process the formatting options (blast_args.cpp:3631-3636). The thread count is
    // checked before the inputs are read, as before. LOSAT's own debug and timing output
    // (LOSAT_WASI_THREADS_DEBUG, LOSAT_TIMING, LOSAT_DEBUG_SCAN_SOFF) follows the reads and
    // does not count their time.
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
    // ```c
    //     if (GetExportSearchStrategyStream(args) ||
    //            m_FormattingArgs->ArchiveFormatRequested(args)) {
    //         locality = CBlastOptions::eBoth;
    //     }
    // ```
    // `ArchiveFormatRequested` parses `-outfmt` (blast_args.cpp:2745-2748) before the
    // option handlers run; a format that NCBI runs and LOSAT does not write is rejected
    // there too.
    let choice = crate::blastinput::app::parse_formatting_string(&args.outfmt)?;
    let format = crate::blastinput::app::report_format(
        &choice,
        crate::blastinput::app::Program::Tblastx,
        args.out.as_deref(),
    )?;
    // The BLAST archive (outfmt 11): its request, and the refusals of what LOSAT cannot write,
    // before anything is read or written.
    // The error archives follow the `-outfmt` text (`writes_error_archive`).
    args.archive = super::super::archive::handle_for(
        &args,
        format == Some(crate::blastinput::app::ReportFormat::ArchiveFormat),
    )?
    .map(|handle| handle.with_outfmt_text(&args.outfmt));
    let archive = args.archive.clone();
    crate::utils::threading::validate_threads(args.num_threads)?;
    use crate::blastinput::input_files;
    // The diagnostics stream: with the archive it also saves what it writes as the messages of
    // the next archive.
    let mut stderr = std::io::stderr();
    let mut archive_tee;
    let diagnostics: &mut (dyn std::io::Write + Send) = match &archive {
        Some(handle) => {
            archive_tee = crate::report::archive::DiagTee::new(&mut stderr, handle.clone());
            &mut archive_tee
        }
        None => &mut stderr,
    };
    // `-db`: no subject file; the database is opened after the options are processed
    // (`search_cli`), where NCBI's `InitializeSubject` opens it (tblastx_app.cpp:113-119).
    let subjects_read = if args.db.is_some() {
        Ok((Vec::new(), Placements::default(), None))
    } else {
        read_subject_set(&args, settings.clone(), &mut *diagnostics)
    };
    // NCBI reads the `-subject` file inside `Run`, after the diagnostics handler that saves the
    // messages is set: an error archive follows the error.
    let (subjects, subject_placements, x_stores) = match subjects_read {
        Ok(read) => read,
        Err(error) => {
            write_run_error_archive(&args, archive.as_ref(), &error);
            return Err(error);
        }
    };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3456-3481
    // ```c
    //     if (args.Exist(kArgQuery) && args[kArgQuery].HasValue() &&
    //         m_InputStream == NULL) {
    // ...
    //         else {
    //             m_InputStream = &args[kArgQuery].AsInputFile();
    //         }
    //     }
    // ...
    //     else {
    //         m_OutputStream = &args[kArgOutput].AsOutputFile();
    //     }
    // ```
    // The query is opened, then the output file is created (`-` is standard output); the
    // query is read after both, so an output file that is the query file reads empty.
    input_files::check_utf8_file_name(&args.query, "query", "TBLASTX")?;
    // NCBI opens the query file inside `Run` too: its error is saved and an error archive
    // follows.
    let query_file = match input_files::open_input(&args.query, "query", "TBLASTX") {
        Ok(file) => file,
        Err(error) => {
            write_run_error_archive(&args, archive.as_ref(), &error);
            return Err(error);
        }
    };
    // Whether this run creates `-out` (a Q3 refusal removes a file it created for nothing).
    let out_is_new = args
        .out
        .as_deref()
        .is_some_and(|path| path.as_os_str() != "-" && !path.exists());
    let out_file = open_out_file(&args)?;
    let out_path = args.out.clone();
    let pairwise = format == Some(crate::blastinput::app::ReportFormat::Pairwise);
    // LOSAT's reports read the format number that NCBI reads: the number and the custom fields,
    // or the value as given when it has a delimiter, which `normalized` drops and the custom
    // tabular report reads (`report::output_format`, which parses it again as NCBI does).
    let outfmt = if choice.delimiter.is_empty() {
        choice.normalized()
    } else {
        args.outfmt.clone()
    };
    let mut stream = cli_report_stream(out_file);
    let result = {
        let mut outputs =
            ReportOutputs::single(&outfmt, OutputSink::Writer(&mut stream), &mut *diagnostics);
        search_cli(
            args,
            query_file,
            &subjects,
            &subject_placements,
            x_stores.as_deref(),
            &mut outputs,
            settings,
        )
    };
    // NCBI writes an error archive after the archives of the batches before when the
    // application stops with one of its own errors, and when the run ends with messages that
    // no archive carried (`PrintErrorArchive` after its `CATCH_ALL`, whatever the status).
    if let Some(handle) = &archive {
        match &result {
            Err(error) => {
                if let Some(lines) = super::super::archive::error_lines(error, false) {
                    let _ =
                        crate::report::archive::write_error_archive(handle, &lines, &mut stream);
                }
            }
            Ok(()) => {
                if handle.has_messages() {
                    let _ = crate::report::archive::write_error_archive(handle, &[], &mut stream);
                }
            }
        }
    }
    // What was written before an error (such as the outfmt 0 prolog) stays in the file.
    let flushed = std::io::Write::flush(&mut stream);
    // A title that the archive cannot hold is refused where NCBI would write its Bioseq, which
    // is known only after the search. BLASTX/TBLASTN refuse before they create `-out`, and so
    // does this run when the file was not there and nothing was written to it: no empty file is
    // left behind.
    let write_failed = stream.failed;
    drop(stream);
    if out_is_new
        && matches!(&result, Err(error) if error.downcast_ref::<crate::report::archive::TitleRefusal>().is_some())
    {
        if let Some(path) = out_path.as_deref() {
            if std::fs::metadata(path).is_ok_and(|meta| meta.len() == 0) {
                let _ = std::fs::remove_file(path);
            }
        }
    }
    // NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.hpp:252-255
    // ```c
    //     catch (const std::ios::failure&) {                                      \
    //         LOG_POST(Error << "BLAST failed to write output");                  \
    //         exit_code = BLAST_OUTPUT_ERROR;                                     \
    //     }                                                                       \
    // ```
    // The outfmt 0 formatter's stream throws when a write fails (blast_format.cpp:118-119).
    // For outfmt 6/7 NCBI aborts instead, and LOSAT reports the error
    // (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES).
    if write_failed && pairwise {
        return Err(crate::cli::NativeError {
            exit: 6,
            message: "BLAST failed to write output\n".to_string(),
        }
        .into());
    }
    result?;
    flushed.context("failed to write the output")
}

/// The error archive of a run that stopped in `Run` with `error` before its search wrote
/// anything, when NCBI writes one (`PrintErrorArchive`): into `-out`, or nothing when that
/// cannot be opened (NCBI's call swallows that).
fn write_run_error_archive(
    args: &TblastxArgs,
    archive: Option<&crate::report::archive::ArchiveHandle>,
    error: &anyhow::Error,
) {
    let Some(handle) = archive else {
        return;
    };
    let Some(lines) = super::super::archive::error_lines(error, true) else {
        return;
    };
    if let Ok(out_file) = open_out_file(args) {
        let mut report = cli_report_stream(out_file);
        let _ = crate::report::archive::write_error_archive(handle, &lines, &mut report);
        let _ = std::io::Write::flush(&mut report);
    }
}

/// The `-out` file, created (`-`: standard output).
fn open_out_file(args: &TblastxArgs) -> Result<Option<std::io::BufWriter<std::fs::File>>> {
    if let Some(path) = args.out.as_deref() {
        crate::blastinput::input_files::check_utf8_file_name(path, "out", "TBLASTX")?;
    }
    Ok(
        match args.out.as_deref().filter(|path| path.as_os_str() != "-") {
            Some(path) => Some(std::io::BufWriter::new(
                std::fs::File::create(path).map_err(|_| crate::cli::inaccessible("out", path))?,
            )),
            None => None,
        },
    )
}

/// The report stream of the command line: the `-out` file, or standard output.
fn cli_report_stream(
    out_file: Option<std::io::BufWriter<std::fs::File>>,
) -> crate::cli::ReportStream {
    crate::cli::ReportStream {
        inner: match out_file {
            Some(file) => Box::new(file) as Box<dyn std::io::Write + Send>,
            None => crate::cli::report_standard_output(),
        },
        failed: false,
    }
}

/// The subjects of `-subject` (after `-subject_loc`) as NCBI's option handlers read them, with
/// their messages.
fn read_subject_set(
    args: &TblastxArgs,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    diagnostics: &mut dyn std::io::Write,
) -> Result<(
    Vec<FastaRecord>,
    Placements,
    Option<Vec<crate::sequence::packed_subject::PackedNucleotide>>,
)> {
    use crate::blastinput::input_files;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2553-2557
    // ```c
    //         CRef<blast::CBlastQueryVector> subjects;
    //         m_Scope = ReadSequencesToBlast(*subj_input_stream, IsProtein(),
    //                                        subj_range, parse_deflines,
    //                                        use_lcase_masks, subjects, m_IsMapper);
    //         m_Subjects.Reset(new blast::CObjMgr_QueryFactory(*subjects));
    // ```
    // The first handler opens and reads the subjects (an empty subject set fails there),
    // as BLASTN's `run`.
    let Some(subject_path) = args.subject.clone() else {
        return Err(crate::blastinput::app::missing_subject_error());
    };
    input_files::check_utf8_file_name(&subject_path, "subject", "TBLASTX")?;
    let subject_file = input_files::open_input(&subject_path, "subject", "TBLASTX")?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2537-2545
    // ```c++
    //             subj_input_stream = &args[kArgSubject].AsInputFile();
    //         }
    //
    //         TSeqRange subj_range;
    //         if (args.Exist(kArgSubjectLocation) && args[kArgSubjectLocation]) {
    //             subj_range =
    //                 ParseSequenceRange(args[kArgSubjectLocation].AsString(),
    //                             "Invalid specification of subject location");
    //         }
    // ```
    // The subject range is read after the subject file is opened, before it is read.
    let subject_range =
        parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject, "TBLASTX")?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:242-246
    // ```c++
    //     CRef<CBlastFastaInputSource> fasta(new CBlastFastaInputSource(in, iconfig));
    //     CRef<CBlastInput> input(new CBlastInput(fasta));
    //     CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
    //     sequences = input->GetAllSeqs(*scope);
    //     return scope;
    // ```
    // NCBI reads the records one at a time, writes each one's messages as it reads it, and
    // checks each one's range after reading it: a range that starts past the end of a
    // record stops the reading there (`read_subjects`). Records without residues (and
    // intervals without letters) are read without a message and reported when the search
    // is set up (`search`).
    let mut subject_source = FastaInputSource::from_argument(
        &subject_path,
        subject_file,
        ReaderConfig::subject("TBLASTX", false, settings.data_loaders)
            .with_defline(args.archive.is_some()),
    );
    // Diagnostic (LOSAT_TIMING): the subject FASTA for the common line (utils/xtiming.rs); with
    // LOSAT_X_TBXSTORE the packing of each record is in this interval too. No NCBI counterpart.
    let t_read_fasta = crate::utils::xtiming::clock();
    // EXPERIMENT (LOSAT_X_TBXSTORE): each record's letters packed as it is read
    // (`x_read_packed_subjects`, cut to `-subject_loc` first); the search reads the stores as
    // a database's (`Subjects::Stored`), the reports by record. A set with a byte without an
    // ncbi4na code (the FASTA reader gives none: it reads `U` as `T`) or a record over `u32::MAX`
    // bases keeps every record's letters. The archive of `-subject` writes every subject's whole
    // letters (and its deflines) in every batch, which the packed stores do not keep (they hold
    // the interval of a `-subject_loc` and no letters), so a run that writes it reads the records
    // as the switch's `=0` side does (as TBLASTN's `x_tbn_store`).
    let mut warn = |message: &[u8]| std::io::Write::write_all(&mut *diagnostics, message);
    let (subjects, subject_placements, stores) =
        if super::super::subjects::x_tbx_store() && args.archive.is_none() {
            let (subjects, stores, placements) =
                crate::blastinput::x_packed_subjects::x_read_packed_subjects(
                    &mut subject_source,
                    subject_range.as_ref(),
                    &mut warn,
                    args.num_threads,
                    "TBLASTX",
                )?;
            crate::utils::xreach::x_reach!("LOSAT_X_TBXSTORE", stores.is_some());
            (subjects, placements, stores)
        } else {
            let read_subject_records =
                read_subjects(&mut subject_source, subject_range.as_ref(), &mut warn)?;
            match cut_subjects(&read_subject_records, subject_range.as_ref())? {
                Some((cut, placements)) => (cut, placements, None),
                None => (read_subject_records, Placements::default(), None),
            }
        };
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
    drop(subject_source);
    check_subjects_not_empty(&subjects)?;
    Ok((subjects, subject_placements, stores))
}

/// NCBI's error for a subject file without records, raised when it reads the subjects.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
/// ```c
/// CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
///     : m_QueryVector(& queries)
/// {
///     if (queries.Empty()) {
///         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
///     }
/// ```
/// NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.hpp:225-227
/// ```c
///             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
///             exit_code = BLAST_ENGINE_ERROR;                                 \
/// ```
fn check_subjects_not_empty(subjects: &[FastaRecord]) -> Result<()> {
    if subjects.is_empty() {
        return Err(crate::cli::NativeError {
            exit: 3,
            message: "BLAST engine error: Empty CBlastQueryVector\n".to_string(),
        }
        .into());
    }
    Ok(())
}

/// NCBI's processing of the options of tblastx after the files are opened: the handlers in
/// the order of `CTblastxAppArgs` (the -seg value; the formats of other programs and the
/// warning for a hit list size below 5), then `Validate` (the threshold and the word size
/// of the lookup table, then the e-value), with NCBI's errors.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2975-2977
/// ```c
///     if(hitlist_size < 5){
///    		ERR_POST(Warning << "Examining 5 or more matches is recommended");
///     }
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_options.c:1303-1310,1358-1364
/// ```c
///     if (program_number != eBlastTypeBlastn &&
///         program_number != eBlastTypeMapping &&
///         (!Blast_ProgramIsRpsBlast(program_number)) &&
///         options->threshold <= 0)
///     {
///         Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
///                          "Non-zero threshold required");
/// ...
///         else {
///             Blast_MessageWrite(blast_msg, eBlastSevError,
///                                kBlastMessageNoContext,
///                                "Word-size must be less "
///                                "than 6 for protein comparison");
///             return BLASTERR_OPTION_VALUE_INVALID;
///         }
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_options.c:1518-1523
/// ```c
///     if (options->expect_value <= 0.0 && options->cutoff_score <= 0)
///     {
///         Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
///          "expect value or cutoff score must be greater than zero");
///         return BLASTERR_OPTION_VALUE_INVALID;
///     }
/// ```
/// TBLASTX has no cutoff score option, so an `-evalue` of 0 (or one that reads as 0) fails.
fn check_ncbi_options(
    args: &TblastxArgs,
    outfmts: &[i32],
    diagnostics: &mut dyn std::io::Write,
) -> Result<(Option<SequenceRange>, ReportCounts)> {
    use crate::blastinput::app::{
        formatting_handler_check, options_error, parse_formatting_string, Program,
    };
    args.seg_spec()?;
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1995-1999
    // ```c++
    //     // set the sequence range
    //     if (args.Exist(kArgQueryLocation) && args[kArgQueryLocation]) {
    //         m_Range = ParseSequenceRange(args[kArgQueryLocation].AsString(),
    //                                      "Invalid specification of query location");
    //     }
    // ```
    // The query options handler comes after the filtering handler and before the
    // formatting handler (tblastx_args.cpp:44-125).
    let query_range = crate::blastinput::seq_range::parse_optional_range(
        args.query_loc.as_deref(),
        crate::blastinput::seq_range::RangeRole::Query,
        "TBLASTX",
    )?;
    formatting_handler_check(&parse_formatting_string(&args.outfmt)?, Program::Tblastx)?;
    // The hit list size and the counts of the pairwise report
    // (`blastinput::report_counts`), with the handler's warnings.
    let (counts, warnings) = crate::blastinput::report_counts::run_report_counts(
        "tblastx",
        outfmts,
        args.max_target_seqs,
        args.num_descriptions,
        args.num_alignments,
    )?;
    diagnostics.write_all(&warnings)?;
    // (`threshold <= 0` is false for NaN, as in C.)
    if args.word_threshold() <= 0.0 {
        return Err(options_error("Non-zero threshold required"));
    }
    if args.word_size() > 4 {
        return Err(options_error(
            "Word-size must be less than 6 for protein comparison",
        ));
    }
    // The hit saving options' check of the hit list size comes before the e-value's.
    counts.validate()?;
    if args.evalue <= 0.0 {
        return Err(options_error(
            "expect value or cutoff score must be greater than zero",
        ));
    }
    Ok((query_range, counts))
}

/// The number of `-outfmt` of the command line (`check_ncbi_options`).
fn cli_outfmts(args: &TblastxArgs) -> Result<[i32; 1]> {
    Ok([crate::blastinput::app::parse_formatting_string(&args.outfmt)?.number])
}

/// The options that LOSAT's TBLASTX rejects, where NCBI would start the search.
pub(crate) fn check_losat_limits(args: &TblastxArgs) -> Result<()> {
    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:214-224
    // ```c
    //     if (ewp->diag_table->multiple_hits) {
    //         status = s_BlastAaWordFinder_TwoHit(subject, query,
    //     ...
    //     } else {
    //         status = s_BlastAaWordFinder_OneHit(subject, query,
    // ```
    // LOSAT's TBLASTX scans with the two-hit word finder only; the one-hit word finder of
    // a window of 0 is not ported.
    if args.window() == 0 {
        anyhow::bail!(
            "-window_size 0 (the one-hit word finder) is not supported by LOSAT's TBLASTX"
        );
    }
    // NCBI reference: c++/src/algo/blast/core/blast_aalookup.c:237
    // ```c
    //     lookup->word_length = opt->word_size;
    // ```
    // LOSAT's TBLASTX builds and scans 3-residue words only.
    if args.word_size() != 3 {
        anyhow::bail!(
            "-word_size {} is not supported by LOSAT's TBLASTX (it implements word size 3)",
            args.word_size()
        );
    }
    Ok(())
}

/// The checks of the options alone, without the inputs (the `validate` of web ABI v2):
/// NCBI's processing of the options and LOSAT's limits.
pub fn check_options(args: &TblastxArgs) -> Result<()> {
    // The host has no database files (`run_local_with`).
    if args.db.is_some() {
        anyhow::bail!("-db (a database search) is not supported by LOSAT Web yet");
    }
    check_ncbi_options(args, &cli_outfmts(args)?, &mut std::io::sink())?;
    check_losat_limits(args)?;
    // The host validates the argv it runs, with its `-num_threads` (docs/web/abi_v2.md).
    crate::utils::threading::validate_threads(args.num_threads)
}

/// The part of `run` after the output is opened: the options, `Query is Empty!`, the
/// reading of the queries and the search.
fn search_cli(
    mut args: TblastxArgs,
    query_file: std::fs::File,
    subjects: &[FastaRecord],
    subject_placements: &Placements,
    x_stores: Option<&[crate::sequence::packed_subject::PackedNucleotide]>,
    outputs: &mut ReportOutputs<'_>,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    let (query_range, counts) =
        check_ncbi_options(&args, &cli_outfmts(&args)?, outputs.diagnostics)?;
    args.report_counts = counts;
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:113-119
    // ```c++
    //         /*** Initialize the database/subject ***/
    //         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
    //         CRef<CLocalDbAdapter> db_adapter;
    //         CRef<CScope> scope;
    //         InitializeSubject(db_args, opts_hndl, m_CmdLineArgs->ExecuteRemotely(),
    //                          db_adapter, scope);
    // ```
    // `-db`: the database opens after the options are processed, before `Query is Empty!`
    // (`TblastxDatabase::open`).
    // The taxid list (line DI-8) is read and expanded first: NCBI's handler of the database
    // arguments reads it (`s_GetTaxIDList`) before the database opens.
    let database = match args.database_name() {
        Some(name) => {
            let taxids = args.taxid_filter.resolve(&settings, outputs.diagnostics)?;
            Some(super::super::database::TblastxDatabase::open_restricted(
                name,
                &settings,
                taxids.as_ref(),
                outputs.diagnostics,
            )?
            .showing_gis(args.show_gis))
        }
        None => None,
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1093-1104
    // ```c++
    // void CBlastTabularInfo::x_CheckTaxDB()
    // {
    //    	if( x_IsFieldRequested(eSubjectSciNames) ||
    //     ...
    //    		string resolved = SeqDB_ResolveDbPath("taxdb.bti");
    // ```
    // A tabular field list that shows a name of the taxonomy database asks whether `taxdb.bti`
    // is on the search path: its formatters warn without it (the plural names), and LOSAT
    // refuses a name of a subject with a taxid with it (`tabular`), as BLASTN and TBLASTN.
    if outputs.formats.iter().any(|format| {
        matches!(
            report::output_format(format.outfmt),
            report::TblastxOutputFormat::CustomTabular(spec) if spec.needs().taxonomy_names
        )
    }) {
        args.tabular_taxdb = Some(match database.as_ref() {
            Some(database) => database.taxdb_found(),
            None => crate::algorithm::common::database::subject_taxdb_found(&settings, "TBLASTX")?,
        });
    }
    // LOSAT's refusals of a database for the Seq-annot formats, before the search
    // (`report::serial::seqalign::check_database`).
    if let Some(database) = database.as_ref() {
        for format in &outputs.formats {
            if let report::TblastxOutputFormat::Seqalign(serial) =
                report::output_format(format.outfmt)
            {
                crate::report::serial::seqalign::check_database(
                    serial,
                    database.name(),
                    &database.pairwise_database().definition,
                    crate::blastinput::app::Program::Tblastx,
                )?;
            }
        }
    }
    // The BLAST XML (outfmt 5, and the XML2 and JSON documents of 16 and 15) cannot write a
    // `-db` name outside US-ASCII: refused before the search (`names::check_database_name`).
    if let Some(database) = database.as_ref() {
        if outputs
            .formats
            .iter()
            .any(|format| report::output_format(format.outfmt).is_structured_xml())
        {
            crate::report::blast_xml::names::check_database_name("TBLASTX", database.name())?;
        }
    }
    let source = match &database {
        Some(database) => SubjectSource::Database(database),
        None => SubjectSource::Fasta {
            records: subjects,
            placements: subject_placements,
            x_stores,
        },
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:125-135
    // ```c++
    //         SDataLoaderConfig dlconfig =
    //             InitializeQueryDataLoaderConfiguration(query_opts->QueryIsProtein(),
    //                                                    db_adapter);
    //         CBlastInputSourceConfig iconfig(dlconfig, query_opts->GetStrand(),
    //                                      query_opts->UseLowercaseMasks(),
    //                                      query_opts->GetParseDeflines(),
    //                                      query_opts->GetRange());
    //         if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())){
    //            	ERR_POST(Warning << "Query is Empty!");
    //            	return BLAST_EXIT_SUCCESS;
    //         }
    // ```
    //
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-860
    // ```c++
    // 	char c;
    // 	CNcbiStreampos orig_p = in.tellg();
    // 	// Piped input
    // 	if(orig_p < 0)
    // 		return false;
    // ```
    // The query's data loaders are those of the subjects (`SDataLoaderConfig`). The position
    // is taken on the opened file, before it is read (`stream_is_empty`), as BLASTN's
    // `search_cli`. When the subjects were read from standard input too, `cin` has reached
    // its end (a failed stream), so NCBI gets no position for the query either, and its
    // reader is at its end: no query batch is read.
    let shared_standard_input = args.query.as_os_str() == "-"
        && args
            .subject
            .as_deref()
            .is_some_and(|path| path.as_os_str() == "-");
    let mut query_source = FastaInputSource::from_argument(
        &args.query,
        query_file,
        ReaderConfig::query("TBLASTX", false, settings.data_loaders)
            .with_defline(args.archive.is_some()),
    );
    if !shared_standard_input && query_source.stream_is_empty() {
        outputs
            .diagnostics
            .write_all(b"Warning: [tblastx] Query is Empty!\n")?;
        // NCBI returns from `Run` here, before its `PrintErrorArchive`: the message is not an
        // archive's.
        if let Some(handle) = &args.archive {
            handle.take_messages();
        }
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:176-178
    // ```c++
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
    //
    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
    // ```
    // The records are read here, and their messages are written batch by batch
    // (`run_in_pool`). What ends the reading (the end of the input, an `eEOF`, a reader
    // error or LOSAT's rejection of a Seq-id line) comes in the batch that NCBI reads it in,
    // after the reports of the batches before (`run_in_pool`).
    // Diagnostic (LOSAT_TIMING): the query FASTA for the common line (utils/xtiming.rs). No NCBI
    // counterpart.
    let t_read_fasta = crate::utils::xtiming::clock();
    let QueryRecords { records, end } = read_queries(&mut query_source);
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
    drop(query_source);
    let reading = match end {
        QueryEnd::Input => QueryReading::End,
        QueryEnd::EofError => QueryReading::BlankLines,
        QueryEnd::Error { error, warnings } => QueryReading::Error { error, warnings },
    };
    search(args, &records, reading, query_range, source, outputs, None)
}

// NCBI reference: ncbi-blast/c++/src/app/blast/blast_formatter.cpp:429-467
// ```c
// CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
// formatter.PrintProlog();
// ...
// ITERATE(CSearchResultSet, result, *results) {
//     ...
//         formatter.PrintOneResultSet(**result, queries);
//     ...
// }
// ```
// NCBI formats one result set without searching again; several requested formats
// are several CBlastFormat printers over the same result set.
/// Runs one TBLASTX search over already parsed records and writes every requested
/// output format from the same result (the shared entry of the CLI, web ABI v1 and
/// v2), and gives the caller the final hit list (`ReportOutputs::hits`).
///
/// Each requested `-outfmt` is validated as on the command line, before the search
/// starts. The `-query` and `-subject` values of `args` are used only as display names.
/// The records are those of NCBI's reader (`blastinput::fasta_reader`; their messages
/// are written where NCBI reads them): queries without a record are NCBI's
/// `Query is Empty!`.
pub fn run_local(
    args: TblastxArgs,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    run_local_with(args, query_records, subject_records, outputs, None)
}

/// `run_local`, with the checks that ABI v1 makes of its `bio` records (`v1`,
/// `web_api::v1_tblastx::run_web_pair`).
pub(crate) fn run_local_with(
    mut args: TblastxArgs,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
    v1: Option<&dyn V1Checks>,
) -> Result<()> {
    // The host has no database files; LOSAT Web searches the records it gives.
    if args.db.is_some() {
        anyhow::bail!("-db (a database search) is not supported by LOSAT Web yet");
    }
    // NCBI parses -outfmt before its option handlers read the subjects (BLASTN's
    // `run_local`); `search` checks the formats again.
    check_report_formats(outputs, None)?;
    // No NCBI counterpart: ABI v1 is frozen (plan TD-1; line OF, outfmt-20261019 DECISIONS
    // 10): it adds no format, so it refuses the Seq-annot formats (outfmt 8, 9, 12), and outfmt
    // 10 and 20 and a custom specification (fields or a delimiter) of the tabular formats (line
    // CLI C2, PR 2), with the messages that `check_report_formats` gave them before TBLASTX wrote
    // them. (Its `-outfmt` reader takes 6 alone, `tblastx_v1_outfmt`.)
    if v1.is_some() {
        for format in &outputs.formats {
            let choice = crate::blastinput::app::parse_formatting_string(format.outfmt)?;
            let number = choice.number;
            if matches!(number, 8 | 9 | 10 | 11 | 12 | 20) {
                anyhow::bail!("output format {number} is not supported by LOSAT's TBLASTX");
            }
            if choice.custom() {
                anyhow::bail!(
                    "a custom output format specification (fields or a delimiter) is not supported by LOSAT's TBLASTX"
                );
            }
        }
    }
    // The BLAST archive (outfmt 11) among the formats: its request, the refusals of what LOSAT
    // cannot write, and a diagnostics stream that saves the messages. The records' deflines are
    // those the caller's reader kept; records without one get `>` and their title.
    let archive_requested = outputs
        .formats
        .iter()
        .any(|format| report::output_format(format.outfmt) == report::TblastxOutputFormat::Archive);
    if archive_requested && args.archive.is_none() {
        args.archive = super::super::archive::handle_for(&args, true)?;
        if let Some(handle) = args.archive.clone() {
            return outputs.with_archive_diagnostics(handle, |view| {
                run_local_with(args, query_records, subject_records, view, v1)
            });
        }
    }
    // The subject range and its record checks come where NCBI reads the subjects (`run`):
    // the messages of each record read, up to a record whose range starts past its end.
    let subject_range =
        parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject, "TBLASTX")?;
    for record in &subject_records[..subjects_read(subject_records, subject_range.as_ref())] {
        outputs.diagnostics.write_all(&record.warnings)?;
    }
    let ranged_subjects = cut_subjects(subject_records, subject_range.as_ref())?;
    check_subjects_not_empty(subject_records)?;
    let mut outfmts = Vec::with_capacity(outputs.formats.len());
    for format in &outputs.formats {
        outfmts.push(crate::blastinput::app::parse_formatting_string(format.outfmt)?.number);
    }
    let (query_range, counts) = check_ncbi_options(&args, &outfmts, outputs.diagnostics)?;
    args.report_counts = counts;
    let whole_subjects = Placements::default();
    let (searched_subjects, subject_placements) = match &ranged_subjects {
        Some((cut, placements)) => (cut.as_slice(), placements),
        None => (subject_records, &whole_subjects),
    };
    search(
        args,
        query_records,
        QueryReading::Records,
        query_range,
        SubjectSource::Fasta {
            records: searched_subjects,
            placements: subject_placements,
            x_stores: None,
        },
        outputs,
        v1,
    )
}

/// Each requested format as NCBI reads `-outfmt`: NCBI's errors, and a rejection of a
/// format that LOSAT's TBLASTX does not write.
///
/// `out` is the run's `-out`: outfmt 13 and 14 write side files next to it, and fail on standard
/// output (`xinclude_check_for`). `run_local` has no output file.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2887
/// ```c++
///     s_ValidateCustomDelim(m_CustomOutputFormatSpec,m_CustomDelim);
/// ```
/// The formatting handler checks the custom delimiter of each format that LOSAT writes too
/// (`formatting_handler_check`; the command line's comes before, in `check_ncbi_options`).
fn check_report_formats(outputs: &ReportOutputs<'_>, out: Option<&std::path::Path>) -> Result<()> {
    use crate::blastinput::app;
    for format in &outputs.formats {
        let choice = app::parse_formatting_string(format.outfmt)?;
        if app::report_format(&choice, app::Program::Tblastx, out)?.is_none() {
            app::formatting_handler_check(&choice, app::Program::Tblastx)?;
            app::xinclude_check_for(&choice, out)?;
        } else {
            app::formatting_handler_check(&choice, app::Program::Tblastx)?;
        }
    }
    Ok(())
}

/// The search of `run_local` and of the CLI, after NCBI's checks of the options and the
/// reading of the subjects: `Query is Empty!` for records without a query, the query batch
/// size, the subjects without letters, and the query batches.
#[allow(clippy::too_many_arguments)]
fn search(
    args: TblastxArgs,
    query_records: &[FastaRecord],
    reading: QueryReading,
    query_range: Option<SequenceRange>,
    source: SubjectSource<'_>,
    outputs: &mut ReportOutputs<'_>,
    v1: Option<&dyn V1Checks>,
) -> Result<()> {
    // NCBI reads the subjects before the queries (`run`); the records are those of NCBI's
    // reader (`fasta_reader`; `run_local`'s callers read them so). The CLI has already
    // reported an empty query file and read its queries; the records of the other callers
    // get NCBI's warning for an input without records.
    if query_records.is_empty() && matches!(reading, QueryReading::Records) {
        outputs
            .diagnostics
            .write_all(b"Warning: [tblastx] Query is Empty!\n")?;
        return Ok(());
    }
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:136-137
    // ```c
    //         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
    //         CBlastInput input(&fasta, m_CmdLineArgs->GetQueryBatchSize());
    // ```
    // The batch size is read before the formatter is made and before the queries, which
    // NCBI reads one batch at a time after the outfmt 0 prolog (`run_in_pool`).
    let batch_size = crate::blastinput::app::query_batch_size("TBLASTX", 10002)?;
    // ABI v1's checks of its `bio` records come where they came before the records of
    // NCBI's reader (plan TD-1).
    if let Some(v1) = v1 {
        v1.check_search(&args, query_range.as_ref(), outputs)?;
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:129-136
    // ```c
    //     if(m_IsDbScan) {
    // 	int num_seqs=0;
    //         int total_length=0;
    // 	if (!is_remote_search)
    //         {
    //                 BlastSeqSrc* seqsrc = db_adapter.MakeSeqSrc();
    //                 num_seqs=BlastSeqSrcGetNumSeqs(seqsrc);
    //                 total_length=static_cast<int>(BlastSeqSrcGetTotLen(seqsrc));
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:773-789
    // ```c
    //         catch(CBlastException & e ) {
    //         	// Skip bad subject sequence
    //         	if(e.GetErrCode() == CBlastException::eInvalidArgument) {
    //         		seqblk_vec->push_back(subj);
    //         ...
    //         		warning += "Subject sequence contains no data";
    //         		ERR_POST(Warning << warning);
    //         		continue;
    // ```
    // The formatter, made after the query batch size and before the outfmt 0 prolog, sets
    // up the subjects (as BLASTN's `search`): a subject without letters (a record without
    // residues, or an interval that starts just past its record's end) gets its warning,
    // stays in the database statistics with no letters, and is never searched
    // (`process_subject` skips subjects shorter than a word).
    // A database has no such warnings (its sequences are not set up by the formatter).
    match &source {
        SubjectSource::Fasta {
            records,
            x_stores: Some(stores),
            ..
        } => crate::blastinput::input_files::write_empty_subject_warnings_by_length(
            records,
            |index| stores[index].len(),
            "tblastx",
            outputs.diagnostics,
        )?,
        SubjectSource::Fasta { records, .. } => {
            crate::blastinput::input_files::write_empty_subject_warnings(
                records,
                "tblastx",
                outputs.diagnostics,
            )?
        }
        SubjectSource::Database(_) => {}
    }
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:170-173
    // ```c
    // 	if(UseXInclude(*fmt_args, args[kArgOutput].AsString())) {
    //         	formatter.SetBaseFile(args[kArgOutput].AsString());
    //         }
    //         formatter.PrintProlog();
    // ```
    // The XInclude formats written to standard output fail after the formatter is made
    // (`check_report_formats`: `xinclude_check`), as the other format checks of each
    // requested format.
    check_report_formats(outputs, args.out.as_deref())?;
    // LOSAT's limits come where NCBI starts the search, after its checks.
    check_losat_limits(&args)?;
    crate::blastinput::app::check_unsupported_environment("TBLASTX")?;
    // The query range applies to every query record as NCBI's batch reader reads it
    // (`seq_range::cut_queries`): a record is searched cut, and one whose interval starts
    // just past its end is a query without data (`run_in_pool`).
    let input_records = query_records;
    let (query_records, query_input) = match query_range.as_ref() {
        Some(range) => {
            let ranged = cut_queries(query_records, range);
            (std::borrow::Cow::Owned(ranged.records), ranged.input)
        }
        None => (
            std::borrow::Cow::Borrowed(query_records),
            QueryInput::whole(query_records),
        ),
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1426
    // ```c
    //     itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    //
    //     /* iterate over all subject sequences */
    //     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //            != BLAST_SEQSRC_EOF) {
    //        ...
    //        if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
    // ```
    // `-db`: the database's sequences, read once for every batch (`TblastxDatabase::read_stores`)
    // where NCBI reads them first, in the preliminary search of the first batch with a valid
    // context (a batch without one is not searched, `CheckInternalData`). A database without
    // sequences (BLASTN's `Database memory map file error`) or a sequence file that cannot be
    // read fails that batch, after the reports of the batches before (`run_in_pool`), and
    // LOSAT's refusal of a damaged file stops the run there; a run without such a batch never
    // reads the file.
    let whole_subjects = Placements::default();
    let (stores, stores_error) = match &source {
        SubjectSource::Database(database) => {
            let read = if database.numbers(None).count == 0 {
                Err(crate::algorithm::common::database::PreliminaryReadError(
                    "Database memory map file error",
                )
                .into())
            } else {
                database.read_stores()
            };
            match read {
                Ok(stores) => (stores, None),
                Err(error) => (Vec::new(), Some(error)),
            }
        }
        SubjectSource::Fasta { .. } => (Vec::new(), None),
    };
    // EXPERIMENT (LOSAT_X_TBXSTORE): the reports read the stores by record.
    let x_ordinals: Vec<u32> = match &source {
        SubjectSource::Fasta {
            x_stores: Some(stores),
            ..
        } => (0..stores.len() as u32).collect(),
        _ => Vec::new(),
    };
    let (run_subjects, subject_placements) = match source {
        SubjectSource::Fasta {
            records,
            placements,
            x_stores: Some(stores),
        } => (
            RunSubjects {
                subjects: Subjects::Stored(stores),
                numbers: SubjectSetNumbers {
                    total: stores.iter().map(|store| store.len() as i64).sum(),
                    count: records.len() as i64,
                    dbsize: args.dbsize,
                },
                shown: ShownSubjects::Fasta {
                    records,
                    x_stored: Some((stores, &x_ordinals)),
                },
            },
            placements,
        ),
        SubjectSource::Fasta {
            records,
            placements,
            x_stores: None,
        } => (
            RunSubjects {
                subjects: Subjects::Fasta(records),
                numbers: SubjectSetNumbers::of_records(records, args.dbsize),
                shown: ShownSubjects::Fasta {
                    records,
                    x_stored: None,
                },
            },
            placements,
        ),
        SubjectSource::Database(database) => (
            RunSubjects {
                subjects: if stores_error.is_some() {
                    Subjects::Unreadable
                } else {
                    Subjects::Stored(&stores)
                },
                numbers: database.numbers(args.dbsize),
                shown: ShownSubjects::Database {
                    database,
                    pairwise: database.pairwise_database(),
                    stores: &stores,
                },
            },
            &whole_subjects,
        ),
    };
    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // TBlastThreads the_threads(GetNumberOfThreads());
    // (*thread)->Run(); (*thread)->Join(&result);
    crate::utils::threading::with_search_pool(args.num_threads, "tblastx", |pool| {
        run_in_pool(
            args,
            &QueryRuns {
                input_records,
                records: &query_records,
                input: &query_input,
                subject_placements,
            },
            reading,
            &run_subjects,
            stores_error,
            outputs,
            pool,
            batch_size,
            v1.is_some(),
        )
    })
}

/// ABI v1's checks of its `bio` records (plan TD-1), which the v1 layer makes on its own
/// records (`web_api::v1_tblastx`): the search calls them where ABI v1 made them when it
/// searched those records (and `check_shown_subject_titles` before the reports), and
/// searches the records of NCBI's reader made from them, which those checks guarantee to be
/// read alike.
pub(crate) trait V1Checks: Sync {
    /// ABI v1's checks where the search starts, in their order: the outfmt 0 and 7 titles,
    /// LOSAT's limits and environment, the residues, the records without residues and the
    /// intervals without letters (the subjects cut to `-subject_loc`, the queries to
    /// `-query_loc`).
    fn check_search(
        &self,
        args: &TblastxArgs,
        query_range: Option<&SequenceRange>,
        outputs: &ReportOutputs<'_>,
    ) -> Result<()>;
}

/// ABI v1's rejection (`V1Checks`) of the outfmt 0 titles that it did not write as NCBI,
/// of the subjects that the reports show: those with hits in the final hit lists, which
/// NCBI's description tables list (blast_format.cpp:1540, `report/defline.rs`). The check
/// comes after the outfmt 0 prolog, where NCBI writes the titles of the first query's
/// report. ABI v1's titles are its `bio` deflines (`web_api::v1_bio::from_bio`), which
/// `check_report_titles` keeps ASCII.
fn check_shown_subject_titles(hits: &[TblastxHsp], subjects: &[FastaRecord]) -> Result<()> {
    let shown: std::collections::BTreeSet<usize> =
        hits.iter().map(|hsp| hsp.hit.s_idx as usize).collect();
    for index in shown {
        crate::report::defline::check_shown_subject_title(
            &String::from_utf8_lossy(&subjects[index].title),
            index + 1,
            "TBLASTX",
        )?;
    }
    Ok(())
}

/// What every report of a run needs besides its HSPs.
pub(crate) struct TblastxReportRun<'a> {
    pub args: &'a TblastxArgs,
    /// The queries reported (those of the batches searched so far) and the subjects.
    pub queries: &'a [FastaRecord],
    pub subjects: &'a [FastaRecord],
    /// The subjects' letters that the displayed rows read (`-db`: the stored codes of the
    /// reported sequences; `-subject` with LOSAT_X_TBXSTORE: the stores of the records; their
    /// records hold none).
    pub subject_letters: report::SubjectLetters<'a>,
    /// `-db`: the database and its reported sequences (`subjects` are their records, by rank).
    pub database: Option<ReportDatabase<'a>>,
    /// The tabular query IDs (`FastaRecord::shown_id`), by query.
    pub query_ids: &'a [Arc<[u8]>],
    /// The per-query statistics and whether the query's batch was searched.
    pub query_stats: &'a [TblastxQueryStats],
    pub searched: &'a [bool],
    /// The warnings written before each query's report (`QueryWarnings`): the reader's
    /// messages of the query's batch (before its first query) and the query's own.
    pub warnings: &'a [Vec<u8>],
    /// The queries' input and ranges.
    pub ranges: &'a QueryRuns<'a>,
    /// Whether the reports end with the epilog (not when the search stops at a batch).
    pub epilog: bool,
    /// The queries whose reports earlier calls wrote, which outfmt 7's epilog counts too
    /// (NCBI's `m_QueriesFormatted`): the batches that the loop wrote when it took them
    /// (`report_stream`); 0 otherwise.
    pub queries_before: usize,
    /// The subject set's numbers (the BLAST XML computes a query's length adjustment again
    /// from them).
    pub numbers: &'a SubjectSetNumbers,
    /// `-db`, XML2 and JSON: the taxid of each reported sequence by rank and what the taxonomy
    /// database says about it (empty for the other formats).
    pub xml_taxa: &'a super::super::database::XmlTaxa,
    /// `-db`, a custom tabular report that shows a taxonomy field: the taxonomy of each reported
    /// sequence by rank (`TblastxDatabase::tabular_taxonomy`; empty otherwise).
    pub tabular_taxonomy: &'a [crate::report::tabular::TabularTaxonomy],
    /// `-db`: the names of the taxonomy database that the rows of the run's custom report (the
    /// command line's one `-outfmt`) look up (`TblastxDatabase::tabular_names`); `None`: every
    /// lookup fails (no name field, no `taxdb.bti`, `-subject`).
    pub tabular_names: Option<&'a crate::report::tabular::TabularNames>,
    /// Whether ABI v1's rejection of outfmt 0 titles applies (`V1Checks`).
    pub v1_titles: bool,
    /// The queries of each batch searched, in order (the documents of the BLAST archive).
    pub archive_spans: &'a [std::ops::Range<usize>],
    /// EXPERIMENT (LOSAT_X_TBXRPTPAR): the pool the report stage runs on (tblastx/x_rpt_par.rs),
    /// `None` for the serial stage.
    pub x_pool: Option<&'a crate::utils::threading::SearchPool<'a>>,
    /// EXPERIMENT (LOSAT_X_PAIRWISEPAR): the pool that formats the outfmt 0 text
    /// (`report::x_pairwise_par`), `None` for the serial writer.
    pub x_pairwise_pool: Option<&'a crate::utils::threading::SearchPool<'a>>,
}

/// The database of a `-db` search as its reports show it.
#[derive(Clone, Copy)]
pub(crate) struct ReportDatabase<'a> {
    /// The `-db` value (outfmt 7's `# Database:` line).
    pub name: &'a str,
    /// The `Database:` lines of outfmt 0.
    pub pairwise: &'a crate::report::pairwise::PairwiseDatabase,
    /// The reported sequences by rank (their tabular ids).
    pub reported: &'a super::super::database::ReportedSubjects,
}

/// The `Database:` text of the reports of a local `-subject` search.
///
/// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:797-799
/// ```c++
/// 	    string dbname;
/// 	    if (m_IsDbScan)
/// 		dbname = string("User specified sequence set (Input: ") + m_SubjectTag + string(")");
/// ```
fn tblastx_database_name(args: &TblastxArgs) -> String {
    format!(
        "User specified sequence set (Input: {})",
        args.subject_label()
    )
}

/// Whether the outfmt 0 prolog of `sink` is written before the search, as NCBI writes it
/// (`write_pairwise_prologs`). A file sink is created when its report is written
/// (`OutputSink::open`), so its prolog stays at the start of the report.
fn prolog_before_search(sink: &OutputSink<'_>) -> bool {
    !matches!(sink, OutputSink::File(_))
}

// NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:173-178
// ```c
//         formatter.PrintProlog();
//
//         /*** Process the input ***/
//         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:118-119
// ```c
// {
//     m_Outfile.exceptions(NcbiBadbit);
// ```
/// Writes and flushes the outfmt 0 prolog of every outfmt 0 sink whose prolog comes before
/// the search (`before_search`), before NCBI reads the first query batch (a failed write
/// stops the run there), or, for a run that stops before any report, of the others.
fn write_pairwise_prologs(
    outputs: &mut ReportOutputs<'_>,
    args: &TblastxArgs,
    subjects: &ShownSubjects<'_>,
    before_search: bool,
) -> Result<()> {
    for format in outputs.formats.iter_mut() {
        if report::output_format(format.outfmt) == report::TblastxOutputFormat::Pairwise
            && prolog_before_search(&format.sink) == before_search
        {
            let mut writer = format.sink.open()?;
            match subjects {
                ShownSubjects::Fasta { records, x_stored } => {
                    crate::report::pairwise::write_tblastx_pairwise_prolog(
                        &mut writer,
                        report::NCBI_TBLASTX_VERSION,
                        &tblastx_database_name(args),
                        records.len(),
                        match x_stored {
                            Some((stores, _)) => stores.iter().map(|store| store.len()).sum(),
                            None => records.iter().map(|r| r.seq().len()).sum(),
                        },
                        None,
                    )?
                }
                ShownSubjects::Database { pairwise, .. } => {
                    crate::report::pairwise::write_tblastx_pairwise_prolog(
                        &mut writer,
                        report::NCBI_TBLASTX_VERSION,
                        "",
                        0,
                        0,
                        Some(pairwise),
                    )?
                }
            }
            writer.flush()?;
        }
    }
    Ok(())
}

// NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:203-209
// ```c
//                 ITERATE(CSearchResultSet, result, *results) {
//                     formatter.PrintOneResultSet(**result, query_batch);
//                 }
//             }
//         }
//
//         formatter.PrintEpilog(opt);
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_seqalign.cpp:1574-1577
// ```c
// // Sort HSPs with e-values as first priority and scores as
// // tie-breakers, since that is the order we want to see them in
// // in Seq-aligns.
// Blast_HSPListSortByEvalue(hsp_list);
// ```
/// Writes the final hits (already filtered by e-value) to every requested format, and
/// gives the caller the same final list (`ReportOutputs::hits`).
fn write_tblastx_outputs(
    mut hits: Vec<TblastxHsp>,
    outputs: &mut ReportOutputs<'_>,
    run: &TblastxReportRun<'_>,
) -> Result<()> {
    let args = run.args;
    let query_code = GeneticCode::from_id(args.query_gencode);
    let db_code = GeneticCode::from_id(args.db_gencode);
    let output_formats: Vec<report::TblastxOutputFormat> = outputs
        .formats
        .iter()
        .map(|format| report::output_format(format.outfmt))
        .collect();
    if run.v1_titles && output_formats.contains(&report::TblastxOutputFormat::Pairwise) {
        check_shown_subject_titles(&hits, run.subjects)?;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1199-1220
    // ```c++
    //     // Set the identity score
    //     if (hsp->num_ident >= 0) {
    //         static const string kNumIdent("num_ident");
    //         scores.push_back(s_MakeScore(kNumIdent, 0.0, hsp->num_ident, true));
    //     }
    //     ...
    //     if (hsp->num_positives > 0) {
    //             static const string kNumPositives("num_positives");
    //             scores.push_back(s_MakeScore(kNumPositives, 0.0, hsp->num_positives, true));
    //     }
    // ```
    // The Seq-aligns (outfmt 8, 12) show the engine's identities and positives
    // (`report::engine_counts`), read before the display counts replace the identities.
    // The positives of the reports are the search matrix's (a name without a table stops the
    // search before any report has a hit, `run_batch_loop`); BLOSUM62's counts keep their
    // constant matrix (`Blosum62Score`).
    let matrix = match run.args.scoring_matrix() {
        TblastxMatrix::Table(matrix) => matrix,
        TblastxMatrix::Unknown | TblastxMatrix::Aborts => ScoringMatrix::Blosum62,
    };
    let (engine_counts, pairwise) = if matrix == ScoringMatrix::Blosum62 {
        report_counts(
            Blosum62Score,
            &mut hits,
            outputs,
            run,
            &output_formats,
            &query_code,
            &db_code,
        )
    } else {
        report_counts(
            matrix,
            &mut hits,
            outputs,
            run,
            &output_formats,
            &query_code,
            &db_code,
        )
    };
    // The tabular rows print record coordinates (`report::shift_to_records`); the outfmt 0
    // rows above were read from the searched letters first.
    report::shift_to_records(
        &mut hits,
        &run.ranges.input.placements,
        run.ranges.subject_placements,
    );
    if let (Some(hits_sink), Some(pairwise)) = (outputs.hits.as_mut(), pairwise.as_ref()) {
        hits_sink(pairwise);
    }
    // Built before any format is written, so that a report that cannot be made fails the
    // run without output (`tabular_subject_ids`).
    let subject_ids = if output_formats.iter().any(|format| format.is_tabular()) {
        match &run.database {
            Some(database) => database.reported.tabular_subject_ids()?,
            None => tabular_subject_ids(&hits, run.subjects)?,
        }
    } else {
        Vec::new()
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:497-498
    // ```c++
    //     //get defline
    //     sdl->defline = CDeflineGenerator().GenerateDefline(m_ScopeRef->GetBioseqHandle(*(sdl->id)), sequence::CDeflineGenerator::fLeavePrefixSuffix);
    // ```
    //
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2273
    // ```c++
    // 	alnDispParams->title = CDeflineGenerator().GenerateDefline(bsp_handle);
    // ```
    // The outfmt 0 descriptions and headings are made from the subjects' title bytes
    // (`report::defline::generate_defline`), by subject index, as BLASTN's.
    let subject_titles: Vec<Arc<[u8]>> =
        if output_formats.contains(&report::TblastxOutputFormat::Pairwise) {
            run.subjects
                .iter()
                .map(|subject| Arc::from(subject.title.as_slice()))
                .collect()
        } else {
            Vec::new()
        };
    let database = match &run.database {
        Some(database) => database.name.to_string(),
        None => tblastx_database_name(args),
    };
    // The `Query=` and `# Query:` text of each query: its title bytes.
    let query_titles: Vec<&[u8]> = run.queries.iter().map(|query| query.title()).collect();
    let unsearched: Vec<bool> = run.searched.iter().map(|&searched| !searched).collect();
    // The warnings of each query are written before its report, once, with the first format
    // (`QueryWarnings`); the taxonomy database's messages of XML2 and JSON go to the
    // diagnostics at every format.
    let mut warning_lines = Some(run.warnings);
    let observer = &mut outputs.observer;
    for (format_index, (format, &output_format)) in
        outputs.formats.iter_mut().zip(&output_formats).enumerate()
    {
        // The warnings that a format other than the archive writes are not the archive's
        // messages (a call with several formats, `run_local`).
        if let Some(handle) = &run.args.archive {
            handle.set_muted(output_format != report::TblastxOutputFormat::Archive);
        }
        let mut probe = observer
            .as_deref_mut()
            .map(|observer| FormatProbe::new(observer, format_index));
        // So do the warnings of the formatters of a tabular field list that shows a taxonomy
        // name (`x_CheckTaxDB`, `tabular`).
        let mut warnings = (warning_lines.is_some()
            || matches!(output_format, report::TblastxOutputFormat::Xml2(_))
            || matches!(
                output_format,
                report::TblastxOutputFormat::CustomTabular(spec) if spec.checks_taxdb()
            ))
        .then(|| crate::report::query_warnings::QueryWarnings {
            before: warning_lines.take().unwrap_or(&[]),
            sink: &mut *outputs.diagnostics,
        });
        // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:68-93
        // ```c
        // CBlastFormat::CBlastFormat(..., CNcbiOstream& outfile, ...)
        //     : m_FormatType(format_type), ..., m_Outfile(outfile),
        // ```
        let prolog = !prolog_before_search(&format.sink);
        let mut writer = format.sink.open()?;
        match output_format {
            report::TblastxOutputFormat::Pairwise => {
                let queries: Vec<crate::report::pairwise::TblastxPairwiseQuery> = run
                    .queries
                    .iter()
                    .zip(run.query_stats)
                    .zip(&query_titles)
                    .enumerate()
                    .map(|(q_idx, ((query, stats), title))| {
                        crate::report::pairwise::TblastxPairwiseQuery {
                            query_name: title.to_vec(),
                            // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:742-744
                            // ```c++
                            //         if(cbs.IsSetInst() && cbs.GetInst().CanGetLength()){
                            //             out << "\nLength=";
                            //             out << cbs.GetInst().GetLength() <<"\n";
                            // ```
                            // The query record's length, also with `-query_loc`.
                            query_length: run
                                .ranges
                                .input
                                .placements
                                .length(q_idx, query.seq().len()),
                            karlin: stats.karlin,
                            effective_search_space: stats.eff_searchsp,
                        }
                    })
                    .collect();
                let (num_descriptions, num_alignments) = (
                    args.report_counts.num_descriptions,
                    args.report_counts.num_alignments,
                );
                let pairwise_report = crate::report::pairwise::TblastxPairwiseReport {
                    version: report::NCBI_TBLASTX_VERSION.to_string(),
                    database_name: database.clone(),
                    // `-db`: the `Database:` lines are the database's (`database` below).
                    database_num_sequences: match &run.database {
                        Some(database) => database.pairwise.num_sequences as usize,
                        None => run.subjects.len(),
                    },
                    database_total_letters: match &run.database {
                        Some(database) => database.pairwise.total_letters as usize,
                        None => (0..run.subjects.len())
                            .map(|s_idx| run.subject_letters.length(s_idx))
                            .sum(),
                    },
                    matrix_name: args.matrix_name().to_string(),
                    word_threshold: args.word_threshold(),
                    window_size: args.window() as usize,
                    num_descriptions,
                    num_alignments,
                    unsearched: unsearched.clone(),
                    // A search that stops at a query batch (`Empty CBlastQueryVector` with
                    // `-query_loc`) writes no epilog (tblastx_app.cpp:209, `PrintEpilog`).
                    epilog: run.epilog,
                    prolog,
                    parsed_subjects: match &run.database {
                        Some(database) => database.reported.parsed.clone(),
                        None => Vec::new(),
                    },
                };
                let config = crate::report::pairwise::PairwiseConfig {
                    program: "tblastx".to_string(),
                    protein_matrix: matrix,
                    ..crate::report::pairwise::PairwiseConfig::default()
                };
                let mut write = |rows: &[crate::report::pairwise::PairwiseHit],
                                 range: std::ops::Range<usize>| {
                    crate::report::pairwise::write_tblastx_pairwise_report(
                        rows,
                        &mut writer,
                        &config,
                        &queries,
                        range,
                        &subject_titles,
                        &pairwise_report,
                        run.database.map(|database| database.pairwise),
                        probe.as_mut(),
                        warnings.as_mut(),
                        run.x_pairwise_pool,
                    )
                };
                match pairwise.as_deref() {
                    Some(rows) => write(rows, 0..queries.len())?,
                    None => {
                        // The rows of a window of queries (`report::pairwise_windows`), written
                        // and dropped before the next window's (line DBID, ID4t): the same rows
                        // as `report_counts` makes for all queries, from the same hits
                        // (`shift_to_records` leaves them as they are here, `rows_up_front`).
                        // Making a row cannot fail, so no window's rows stop the report after
                        // the windows before were written.
                        let shown = report::PairwiseRows::new(
                            run.queries,
                            run.subjects,
                            run.subject_letters,
                            run.query_stats,
                            &query_code,
                            &db_code,
                            &run.ranges.input.placements,
                            run.ranges.subject_placements,
                        );
                        // The final list is in report order, queries in input order: each batch's
                        // lists come in query order (`HitCollector::finish`, `culled_hit_order`)
                        // and the batches follow each other (`run_batch_loop`); `-db` renumbers
                        // the subjects only (`reported_subjects`).
                        debug_assert!(hits.is_sorted_by_key(|hsp| hsp.hit.q_idx));
                        let windows =
                            report::pairwise_windows(&hits, queries.len(), run.x_pairwise_pool);
                        if windows.is_empty() {
                            // No query: the prolog and the epilog alone.
                            write(&[], 0..0)?;
                        }
                        for range in windows {
                            let first =
                                hits.partition_point(|hsp| (hsp.hit.q_idx as usize) < range.start);
                            let end =
                                hits.partition_point(|hsp| (hsp.hit.q_idx as usize) < range.end);
                            let window = &hits[first..end];
                            // EXPERIMENT (LOSAT_X_TBXRPTPAR): the rows as an indexed map on
                            // the pool.
                            let mut rows = match run.x_pool {
                                Some(pool) => crate::algorithm::tblastx::x_rpt_par::pairwise_hits(
                                    pool, window, &shown,
                                ),
                                None => {
                                    let mut buffer = Vec::new();
                                    window
                                        .iter()
                                        .map(|hsp| shown.hit(hsp, &mut buffer))
                                        .collect()
                                }
                            };
                            // The rows count BLOSUM62's positives (`report::set_row_positives`).
                            if matrix != ScoringMatrix::Blosum62 {
                                report::set_row_positives(&mut rows, matrix);
                            }
                            write(&rows, range)?;
                        }
                    }
                }
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:690-698
            // ```c++
            //     if (m_FormatType == CFormattingArgs::eAsnText || m_FormatType == CFormattingArgs::eJsonSeqalign) {
            //         if (results.HasAlignments()) {
            //     	    CRef<CSeq_align_set> aln_set (new CSeq_align_set);
            //             CBlastFormatUtil::PruneSeqalign(*(results.GetSeqAlign()), *aln_set, m_HitlistSize);
            //             if(m_FormatType == CFormattingArgs::eAsnText)
            //             	m_Outfile << MSerial_AsnText << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
            //             else
            //             	m_Outfile << MSerial_Json << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
            //         }
            // ```
            // The Seq-annot of each query with alignments.
            report::TblastxOutputFormat::Seqalign(serial) => {
                crate::algorithm::tblastx::seqalign::write_tblastx_seqalign(
                    serial,
                    &mut writer,
                    &hits,
                    &engine_counts,
                    run.queries,
                    &run.ranges.input.placements,
                    &match &run.database {
                        Some(database) => {
                            crate::algorithm::tblastx::seqalign::SeqalignSubjects::Database {
                                ordinals: &database.reported.ordinals,
                                title: &database.pairwise.definition,
                                parsed: &database.reported.parsed,
                            }
                        }
                        None => crate::algorithm::tblastx::seqalign::SeqalignSubjects::Records {
                            records: run.subjects,
                            placements: run.ranges.subject_placements,
                        },
                    },
                    run.subject_letters,
                    warnings.as_mut(),
                )?
            }
            // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:198
            // ```c++
            //                 if (fmt_args->ArchiveFormatRequested(args)) {
            //                     formatter.WriteArchive(*queries, *opts_hndl, *results, 0, m_Bah.GetMessages());
            //                     m_Bah.ResetMessages();
            // ```
            // One BLAST archive per query batch (`archive`), from the Seq-aligns of outfmt 8.
            report::TblastxOutputFormat::Archive => {
                let handle = run
                    .args
                    .archive
                    .as_ref()
                    .expect("a BLAST archive request has its handle");
                let subjects = match &run.database {
                    Some(database) => {
                        crate::algorithm::tblastx::seqalign::SeqalignSubjects::Database {
                            ordinals: &database.reported.ordinals,
                            title: &database.pairwise.definition,
                            parsed: &database.reported.parsed,
                        }
                    }
                    None => crate::algorithm::tblastx::seqalign::SeqalignSubjects::Records {
                        records: run.subjects,
                        placements: run.ranges.subject_placements,
                    },
                };
                crate::algorithm::tblastx::archive::write_tblastx_archives(
                    handle,
                    &crate::algorithm::tblastx::archive::ArchiveData {
                        queries: run.queries,
                        query_stats: run.query_stats,
                        searched: run.searched,
                        subject_records: match &run.database {
                            Some(_) => &[],
                            None => run.subjects,
                        },
                        hits: &hits,
                        engine: &engine_counts,
                        subjects: &subjects,
                        subject_letters: run.subject_letters,
                        query_placements: &run.ranges.input.placements,
                        numbers: run.numbers,
                        spans: run.archive_spans,
                    },
                    warnings.as_mut(),
                    &mut writer,
                )?
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:708-744
            // ```c++
            //     } else if (m_FormatType == CFormattingArgs::eXml) {
            //         ...
            //         	BlastXML_FormatReport(xml_output, &report_data, &m_Outfile,
            //                                   m_BlastXMLIncremental.GetPointer());
            // ```
            // The document of every query; the hits are the final list after `shift_to_records`
            // (record coordinates), the rows are read from the searched letters.
            report::TblastxOutputFormat::Xml => super::super::blast_xml::write_tblastx_xml(
                &mut crate::report::blast_xml::BlastXmlWriter::new(0),
                false,
                &mut writer,
                &hits,
                run,
                warnings.as_mut(),
            )?,
            // The XML2 and JSON documents (outfmt 16 and 15, and 14 and 13 as side files next to
            // `-out` with this output as the master file): the same records, written by
            // `Xml2Writer`.
            report::TblastxOutputFormat::Xml2(kind) => super::super::blast_xml::write_tblastx_xml(
                &mut crate::report::blast_xml2::Xml2Writer::for_run(
                    kind,
                    super::super::blast_xml::xml2_run(run),
                    run.args.out.as_deref(),
                ),
                true,
                &mut writer,
                &hits,
                run,
                warnings.as_mut(),
            )?,
            // The tabular reports with a field list or a delimiter, and outfmt 10 and 20
            // (`tabular`); never the default 6 and 7 below.
            report::TblastxOutputFormat::CustomTabular(spec) => {
                super::super::tabular::write_tblastx_tabular(
                    &super::super::tabular::TblastxTabular {
                        spec,
                        hits: &hits,
                        run,
                        query_titles: &query_titles,
                        database: &database,
                        subject_ids: &subject_ids,
                        unsearched: &unsearched,
                    },
                    &mut writer,
                    warnings.as_mut(),
                    probe.as_mut(),
                )?
            }
            report::TblastxOutputFormat::Tabular
            | report::TblastxOutputFormat::TabularWithComments => match run.x_pool {
                // EXPERIMENT (LOSAT_X_TBXRPTPAR): the rows formatted on the pool; with the format
                // observer (which flushes around every row) the serial writer.
                Some(pool) if probe.is_none() => {
                    crate::algorithm::tblastx::x_rpt_par::write_tabular(
                        pool,
                        &hits,
                        &mut writer,
                        output_format == report::TblastxOutputFormat::TabularWithComments,
                        &query_titles,
                        &database,
                        run.query_ids,
                        &subject_ids,
                        &unsearched,
                        run.epilog.then_some(run.queries_before),
                        warnings.as_mut(),
                    )?
                }
                _ => report::write_tabular(
                    &hits,
                    &mut writer,
                    output_format == report::TblastxOutputFormat::TabularWithComments,
                    &query_titles,
                    &database,
                    run.query_ids,
                    &subject_ids,
                    &unsearched,
                    run.epilog.then_some(run.queries_before),
                    probe.as_mut(),
                    warnings.as_mut(),
                )?,
            },
        }
        writer.flush()?;
    }
    if let Some(handle) = &run.args.archive {
        handle.set_muted(false);
    }
    // Without an output format the warnings are written all at once.
    if let Some(lines) = warning_lines {
        for query_warnings in lines {
            outputs.diagnostics.write_all(query_warnings)?;
        }
    }
    Ok(())
}

/// Whether the outfmt 0 rows of every query are made before any format is written
/// (`report_counts`): for the caller's hit records, for a format observer with outfmt 0 (its HSP
/// numbers count over the whole report), and for outfmt 0 where the hits move into the records
/// (`-query_loc`, `-subject_loc`: `report::shift_to_records` moves the hits that the rows are
/// read from). Otherwise outfmt 0 makes the rows of a window of queries, writes them and drops
/// them before the next window's (`write_tblastx_outputs`).
///
/// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1411-1585
/// ```c++
/// CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
/// ...
///     display.DisplaySeqalign(m_Outfile);
/// ```
/// NCBI makes and prints the alignments of one query's result set at a time.
fn rows_up_front(
    outputs: &ReportOutputs<'_>,
    run: &TblastxReportRun<'_>,
    output_formats: &[report::TblastxOutputFormat],
) -> bool {
    outputs.hits.is_some()
        || (output_formats.contains(&report::TblastxOutputFormat::Pairwise)
            && (outputs.observer.is_some()
                || run.ranges.input.placements.ranged()
                || run.ranges.subject_placements.ranged()))
}

/// The Seq-aligns' engine counts (outfmt 8, 12), the display counts of `hits` and the outfmt 0
/// rows made up front (`rows_up_front`), with the positives of the search's matrix `sc`
/// (`write_tblastx_outputs`).
#[allow(clippy::too_many_arguments)]
fn report_counts<S: TbxScore>(
    sc: S,
    hits: &mut Vec<TblastxHsp>,
    outputs: &ReportOutputs<'_>,
    run: &TblastxReportRun<'_>,
    output_formats: &[report::TblastxOutputFormat],
    query_code: &GeneticCode,
    db_code: &GeneticCode,
) -> (
    Vec<(u32, u32)>,
    Option<Vec<crate::report::pairwise::PairwiseHit>>,
) {
    // The BLAST archive (outfmt 11) writes the same Seq-aligns.
    let seqalign = output_formats.iter().any(|format| {
        matches!(
            format,
            report::TblastxOutputFormat::Seqalign(_) | report::TblastxOutputFormat::Archive
        )
    });
    let engine_counts = if seqalign {
        // EXPERIMENT (LOSAT_X_TBXRPTPAR): the counts over chunks of the list on the pool.
        match run.x_pool {
            Some(pool) => crate::algorithm::tblastx::x_rpt_par::engine_counts(
                pool,
                hits,
                run.queries,
                run.subject_letters,
                query_code,
                db_code,
                sc,
            ),
            None => report::engine_counts(
                hits,
                run.queries,
                run.subject_letters,
                query_code,
                db_code,
                sc,
            ),
        }
    } else {
        Vec::new()
    };
    // The display counts are read by the outfmt 0, 6 and 7 reports and the caller's hit
    // records only.
    if outputs.hits.is_some()
        || output_formats
            .iter()
            .any(|format| format.reads_displayed_counts())
    {
        // EXPERIMENT (LOSAT_X_TBXRPTPAR): the identities over chunks of the list on the pool.
        match run.x_pool {
            Some(pool) => crate::algorithm::tblastx::x_rpt_par::set_displayed_identities(
                pool,
                hits,
                run.queries,
                run.subject_letters,
                query_code,
                db_code,
                sc,
            ),
            None => report::set_displayed_identities(
                hits,
                run.queries,
                run.subject_letters,
                query_code,
                db_code,
                sc,
            ),
        }
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1411
    // ```c
    // CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
    // ```
    // The caller receives the same final result that every formatter prints, with the
    // displayed rows of outfmt 0. Otherwise outfmt 0 makes its rows a window of queries at a
    // time (`rows_up_front`).
    let mut pairwise = rows_up_front(outputs, run, output_formats)
        // EXPERIMENT (LOSAT_X_TBXRPTPAR): the rows as an indexed map on the pool.
        .then(|| match run.x_pool {
            Some(pool) => crate::algorithm::tblastx::x_rpt_par::pairwise_hits(
                pool,
                hits,
                &report::PairwiseRows::new(
                    run.queries,
                    run.subjects,
                    run.subject_letters,
                    run.query_stats,
                    query_code,
                    db_code,
                    &run.ranges.input.placements,
                    run.ranges.subject_placements,
                ),
            ),
            None => report::pairwise_hits(
                hits,
                run.queries,
                run.subjects,
                run.subject_letters,
                run.query_stats,
                query_code,
                db_code,
                &run.ranges.input.placements,
                run.ranges.subject_placements,
            ),
        });
    // The rows count BLOSUM62's positives (`report::set_row_positives`).
    if !S::BLOSUM62 {
        if let Some(rows) = pairwise.as_mut() {
            report::set_row_positives(rows, sc.matrix());
        }
    }
    (engine_counts, pairwise)
}

/// The subject ID of the tabular formats of the subject record `subject`, `None` where
/// `tabular_subject_ids` refuses it.
fn subject_tabular_id(subject: &FastaRecord) -> Option<Arc<[u8]>> {
    crate::report::defline::tabular_subject_id(subject.local_id.as_bytes(), &subject.title, false)
        .ok()
        .map(Arc::from)
}

/// The subject ID of the tabular formats of every subject, by subject index
/// (`report::defline::tabular_subject_id`: the title's first word, or the local ID, and
/// for `lcl|Subject_...` the first word of `GenerateDefline`), as BLASTN's.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:863-870
/// ```c++
///         } catch (const CException&) {
///             list<CRef<CSeq_id> > subject_ids;
///             CRef<CSeq_id> id(new CSeq_id());
///             id->Assign(align.GetSeq_id(1));
///             subject_ids.push_back(id);
///             SetSubjectId(subject_ids);
///             bioseqs_found = false;
///         }
/// ```
/// A subject with a row whose `GenerateDefline` throws (a `Subject_` first word and a title
/// in no encoding that `CUtf8::GuessEncoding` recognises) makes NCBI 2.17.0 write part of the
/// row and stop with a `CCoreException` that names its build's source files (exit 255;
/// `AUTHORITY.md` §J-6, §K-1, RP-20). No result can be defined for it, so LOSAT rejects such
/// a subject explicitly before it writes any format; a subject without rows is never
/// named.
fn tabular_subject_ids(hits: &[TblastxHsp], subjects: &[FastaRecord]) -> Result<Vec<Arc<[u8]>>> {
    let ids: Vec<Option<Arc<[u8]>>> = subjects.iter().map(subject_tabular_id).collect();
    for hsp in hits {
        let s_idx = hsp.hit.s_idx as usize;
        if ids.get(s_idx).is_some_and(Option::is_none) {
            anyhow::bail!(
                "subject record {} ({}) has a title whose first word starts with 'Subject_' and whose non-UTF-8 bytes are in no encoding that NCBI BLAST+ recognises; NCBI BLAST+ writes part of a tabular row for such a subject and stops with an exception that names its build's source files (exit 255), which is not supported by LOSAT's TBLASTX",
                s_idx + 1,
                String::from_utf8_lossy(subjects[s_idx].shown_id())
            );
        }
    }
    Ok(ids
        .into_iter()
        .map(|id| id.unwrap_or_else(|| Arc::from(&b""[..])))
        .collect())
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input_aux.cpp:85-91
// ```c
//     // used for experimentation purposes
//     char* batch_sz_str = getenv("BATCH_SIZE");
//     if (batch_sz_str) {
//         retval = NStr::StringToInt(batch_sz_str);
//         _TRACE("DEBUG: Using query batch size " << retval);
//         return retval;
//     }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input_aux.cpp:130-134
// ```c
//     case eTblastx:
//         // N.B.: the splitting is done on the nucleotide query sequences, then
//         // each of these chunks is translated
//         retval = 10002;
//         break;
// ```

// NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:136-137
// ```c
//         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
//         CBlastInput input(&fasta, m_CmdLineArgs->GetQueryBatchSize());
// ```
// NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:176-207
// ```c
//         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
//     ...
//                 CLocalBlast lcl_blast(queries, opts_hndl, db_adapter);
//                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
//                 results = lcl_blast.Run();
//     ...
//                 ITERATE(CSearchResultSet, result, *results) {
//                     formatter.PrintOneResultSet(**result, query_batch);
//                 }
//         }
// ```
// Each query batch is searched on its own: the linking cutoffs use the average length
// and the smallest Lambda of the batch's contexts (blast_parameters.c:1023-1026), so the
// batches change the HSPs of a run with several queries. NCBI writes the reports of a
// batch after its search; LOSAT searches the batches one at a time and writes the reports
// of each batch when the loop takes it where no refusal of LOSAT's can come after the first
// batch (`report_stream`), else of all of them once after the loop (`write_tblastx_outputs`),
// so that what stops the run at a batch comes after the reports of the batches before
// (without the epilog).
#[allow(clippy::too_many_arguments)]
fn run_in_pool(
    args: TblastxArgs,
    query_runs: &QueryRuns<'_>,
    reading: QueryReading,
    subjects: &RunSubjects<'_>,
    stores_error: Option<anyhow::Error>,
    outputs: &mut ReportOutputs<'_>,
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
    batch_size: u32,
    v1_titles: bool,
) -> Result<()> {
    check_report_formats(outputs, args.out.as_deref())?;
    let query_records = query_runs.records;
    let input = query_runs.input;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:474-504
    // (`s_ReplaceLocalId`)
    // The query ID of the tabular formats is the title's first word, or the local ID
    // (`FastaRecord::shown_id`), as BLASTN's.
    let query_ids: Vec<Arc<[u8]>> = query_records
        .iter()
        .map(|record| Arc::from(record.shown_id()))
        .collect();
    let reports = BatchReports {
        args: &args,
        query_runs,
        query_ids: &query_ids,
        subjects: &subjects.shown,
        numbers: &subjects.numbers,
        titles: std::cell::RefCell::new(DatabaseTitles::new()),
        v1_titles,
        x_search_pool: parallel_pool,
        queries_written: std::cell::Cell::new(0),
    };
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:176-178
    // ```c
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
    //
    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
    // ```
    // A reader at the end of its input reads no batch (`!input.End()`): the report of no
    // query, the prolog and the epilog (an empty pipe, or `-query -` after `-subject -`).
    let no_batch = query_runs.input_records.is_empty() && matches!(reading, QueryReading::End);
    // NCBI reference: c++/src/algo/blast/core/blast_stat.c:3103-3105, 3136 (the length
    // adjustment's tables of a name that the matrix loader read as IDENTITY,
    // raw_scoremat.c:157-161)
    // ```c
    //    num_values = Blast_GetMatrixValues(matrixName, &gapOpen_arr,
    //      &gapExtend_arr, NULL, NULL, NULL,  &alpha_arr, &beta_arr,
    //      &pref_flags);
    // ...
    //    sfree(gapOpen_arr);
    // ```
    // NCBI 2.17.0 frees arrays that it never set and aborts (`free(): invalid size`, rc 134)
    // when it sets up a batch; a run without a batch reads no matrix. LOSAT refuses a run with
    // a batch before its first byte.
    if !no_batch && args.scoring_matrix() == TblastxMatrix::Aborts {
        anyhow::bail!(
            "-matrix {} is read as IDENTITY without IDENTITY's statistics table (NCBI BLAST+ 2.17.0 aborts); this is not supported by LOSAT's TBLASTX",
            args.matrix_name()
        );
    }
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:173
    // ```c
    //         formatter.PrintProlog();
    // ```
    write_pairwise_prologs(outputs, &args, &subjects.shown, true)?;
    if no_batch {
        return reports.write(outputs, SearchedBatches::default(), &[], true);
    }
    // EXPERIMENT (LOSAT_X_TBXPIPE), dispatch: the batches the loop below will search, planned
    // before it (`x_pipe_plan`), when the pipe may run (`x_pipe_wanted`); `None` keeps the loop's
    // own searches.
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let x_plan = x_pipe_wanted(parallel_pool)
        .then(|| x_pipe_plan(input, &reading, batch_size, query_records));
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:193-195
    // ```c
    //                 results = lcl_blast.Run();
    // ```
    // The search of one batch, the queries `range`; with LOSAT_X_TBXPIPE its seed-stage end
    // (`seed_end`) may start the next batch (tblastx/x_batch_pipe.rs).
    let search_subjects = subjects.subjects;
    let numbers = &subjects.numbers;
    // The search's matrix, chosen once (`scoring.rs`): BLOSUM62's instance of the search, or
    // the other matrices' with the table of `-matrix` (a name without a table stops at the first
    // batch with letters, `run_batch_loop`; one that NCBI aborts on runs no batch).
    let matrix_scores = match args.scoring_matrix() {
        TblastxMatrix::Table(ScoringMatrix::Blosum62)
        | TblastxMatrix::Unknown
        | TblastxMatrix::Aborts => None,
        TblastxMatrix::Table(matrix) => Some(MatrixScores::new(matrix)),
    };
    let search: &x_batch_pipe::XSearch<'_, TblastxBatch> = &|range, seed_end| match &matrix_scores {
        None => search_query_batch(
            Blosum62Score,
            &args,
            &query_records[range],
            search_subjects,
            numbers,
            parallel_pool,
            seed_end,
        ),
        Some(scores) => search_query_batch(
            scores.scorer(),
            &args,
            &query_records[range],
            search_subjects,
            numbers,
            parallel_pool,
            seed_end,
        ),
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:203-209
    // The loop writes each batch's reports when it takes the batch where no refusal of LOSAT's
    // can come after the first batch (`report_stream`); otherwise after the loop.
    let x_check = Instant::now();
    let stream = report_stream::check(
        &args,
        outputs,
        v1_titles,
        query_runs,
        subjects,
        stores_error.is_some(),
        batch_size,
        parallel_pool,
    );
    // Diagnostic (LOSAT_X_STATS): whether the loop writes each batch's reports, or why not. No NCBI
    // counterpart.
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        eprintln!(
            "[X_STATS] TBXRPTSTREAM off={} check_s={:.6}",
            stream.map_or_else(report_stream::Off::name, |()| "none"),
            x_check.elapsed().as_secs_f64()
        );
    }
    let stream = stream.is_ok();
    // The loop's result of the batch of queries `start..end`: the search above, or with
    // LOSAT_X_TBXPIPE the pipe's (`XPipe::take`).
    let batch_loop = |take: &mut dyn FnMut(std::ops::Range<usize>) -> Result<TblastxBatch>| {
        run_batch_loop(
            &args,
            query_runs,
            reading,
            subjects,
            stores_error,
            outputs,
            batch_size,
            &reports,
            stream,
            take,
        )
    };
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    if let Some(plan) = x_plan {
        if let Some(pipe) = x_batch_pipe::XPipe::new(plan, search_subjects.count(), search) {
            let result = rayon::in_place_scope(|scope| {
                let _closer = pipe.closer();
                batch_loop(&mut |range| pipe.take(scope, range))
            });
            for line in pipe.lines() {
                eprintln!("{line}");
            }
            return result;
        }
    }
    batch_loop(&mut |range| search(range, None))
}

// EXPERIMENT (LOSAT_X_TBXPIPE): whether the pipe may run in this search (tblastx/x_batch_pipe.rs,
// "Off"). No NCBI counterpart.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn x_pipe_wanted(parallel_pool: &crate::utils::threading::SearchPool<'_>) -> bool {
    if !x_batch_pipe::on() {
        return false;
    }
    let wanted = parallel_pool.enabled()
        && parallel_pool.threads() > 1
        && rayon::current_thread_index().is_some()
        && !writes_debug_output();
    // The batches that start early are marked where they start (`XPipe::seed_ended`).
    if !wanted {
        crate::utils::xreach::x_reach!("LOSAT_X_TBXPIPE", false);
    }
    wanted
}

// No NCBI counterpart: LOSAT's debugging output, which keeps the batches searched one at a time
// (LOSAT_X_TBXPIPE) and the reports after the batch loop (`report_stream`).
/// Whether the search writes debugging output while it searches a batch (`LOSAT_DEBUG*`, the HSP
/// trace, the diagnostics, the stage dumps).
fn writes_debug_output() -> bool {
    std::env::vars_os().any(|(name, _)| {
        name.to_str()
            .is_some_and(|name| name.starts_with("LOSAT_DEBUG"))
    }) || trace_hsp_target().is_some()
        || diagnostics_enabled()
        || crate::utils::threading::diagnostics_enabled()
        || stage_dump::enabled()
}

// EXPERIMENT (LOSAT_X_TBXPIPE): the batches `run_batch_loop` will search unless a search fails, in
// its order: its batch boundaries and the checks of a batch that need no search (the reader's
// error where the reading ends, a batch without a query, a batch whose queries all lack a codon
// with one without letters). The check of the query splitter's environment concerns the first
// batch only, which the loop starts itself after it. No NCBI counterpart (scheduling only; the
// statements are the loop's, see there).
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn x_pipe_plan(
    input: &QueryInput,
    reading: &QueryReading,
    batch_size: u32,
    query_records: &[FastaRecord],
) -> Vec<std::ops::Range<usize>> {
    let mut plan = Vec::new();
    let reads_past_records = reading.reads_past_records();
    let reading_fails = matches!(reading, QueryReading::Error { .. });
    let mut reading_ended = false;
    let mut input_start = 0;
    while input_start < input.input_lengths.len() || (reads_past_records && !reading_ended) {
        let (input_end, reached_size) = next_ranged_batch(
            &input.input_lengths,
            &input.skipped,
            input_start,
            batch_size,
        );
        let std::ops::Range { start, end } = input.searched_in(input_start..input_end);
        if input_end == input.input_lengths.len() && !reached_size && reads_past_records {
            reading_ended = true;
            if reading_fails {
                break;
            }
        }
        if start == end {
            break;
        }
        let batch_queries = &query_records[start..end];
        if batch_queries.iter().any(|record| record.seq().is_empty())
            && batch_queries
                .iter()
                .all(|record| record.seq().len() < CODON_LENGTH)
        {
            break;
        }
        plan.push(start..end);
        input_start = input_end;
    }
    plan
}

/// The batch loop of `run_in_pool`: reads the batches, checks them, takes each one's search
/// result from `take`, and writes the reports: each batch's when the loop takes it with `stream`
/// (`report_stream`), else all of them after the loop.
#[allow(clippy::too_many_arguments)]
fn run_batch_loop(
    args: &TblastxArgs,
    query_runs: &QueryRuns<'_>,
    reading: QueryReading,
    subjects: &RunSubjects<'_>,
    mut stores_error: Option<anyhow::Error>,
    outputs: &mut ReportOutputs<'_>,
    batch_size: u32,
    reports: &BatchReports<'_>,
    stream: bool,
    take: &mut dyn FnMut(std::ops::Range<usize>) -> Result<TblastxBatch>,
) -> Result<()> {
    let query_records = query_runs.records;
    let input = query_runs.input;
    let mut done = SearchedBatches::default();
    // The warnings written before each query's report (`QueryWarnings`).
    let mut warnings: Vec<Vec<u8>> = vec![Vec::new(); query_records.len()];
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1431-1452
    // ```c++
    //     if (m_FormatType == CFormattingArgs::eAsnText
    //       ...
    //       || m_FormatType == CFormattingArgs::eJsonSeqalign
    //       || m_FormatType == CFormattingArgs::eSAM)
    //     {
    //         x_PrintStructuredReport(results, queries);
    //         return;
    //     }
    //     ...
    //     if (results.HasWarnings()) {
    //         ERR_POST(Warning << results.GetWarningStrings());
    //     }
    // ```
    // A structured format (outfmt 5, 8, 9, 12) writes none of a query's own warnings
    // (`ReportFormat::writes_query_warnings`; outfmt 5 puts them in the Iteration's message,
    // `blast_xml`); the warnings go with the first format.
    let result_warnings = outputs
        .formats
        .first()
        .is_none_or(|format| report::output_format(format.outfmt).writes_query_warnings());
    // Whether the batches' prefetch of a `-db` search writes the taxonomy database's messages
    // with the query warnings: the first format is outfmt 0, 5, 6 or 7
    // (`TblastxOutputFormat::prefetches_subjects`; see the batch's prefetch below).
    let prefetches = outputs
        .formats
        .first()
        .is_some_and(|format| report::output_format(format.outfmt).prefetches_subjects());
    // A custom tabular report with a plural name looks every taxid of its rows up, taxid 0 too,
    // which the prefetch does not: the rows of a batch with hits may open the taxonomy database
    // before a later batch's prefetch (`TaxdbRun::open_for_rows`).
    let plural_names = outputs.formats.first().is_some_and(|format| {
        matches!(
            report::output_format(format.outfmt),
            report::TblastxOutputFormat::CustomTabular(spec) if spec.checks_taxdb()
        )
    });
    // A batch is read from the input records (`input_start..input_end`); its searched
    // queries are `start..end` (all of them without `-query_loc`). A batch size counts the
    // whole records, and a record whose range starts past its end is skipped
    // (`next_ranged_batch`).
    let mut reading = reading;
    // Whether what ended the reading after the records (`reading`) has come.
    let mut reading_ended = false;
    let mut input_start = 0;
    // The reader is at the end of its input after the last record, unless an `eEOF` or an
    // error ends the reading there; that comes in the batch being read when the records
    // run out before the batch reaches its size, and otherwise in a batch of its own.
    while input_start < input.input_lengths.len()
        || (reading.reads_past_records() && !reading_ended)
    {
        let (input_end, reached_size) = next_ranged_batch(
            &input.input_lengths,
            &input.skipped,
            input_start,
            batch_size,
        );
        let std::ops::Range { start, end } = input.searched_in(input_start..input_end);
        let reading_ends_here =
            input_end == input.input_lengths.len() && !reached_size && reading.reads_past_records();
        // `CFastaReader` writes its messages about a batch's queries (their lines and
        // titles) when it reads them, after the report of the batch before: before the
        // report of the batch's first query (`QueryWarnings`). Skipped records are read too.
        let batch_warnings: Vec<u8> = query_runs.input_records[input_start..input_end]
            .iter()
            .flat_map(|record| record.warnings.iter().copied())
            .collect();
        // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input.cpp:146-152
        // ```c
        //         try { q.Reset(m_Source->GetNextSequence(scope)); }
        //         catch (const CObjReaderParseException& e) {
        //             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
        //                 break;
        //             }
        //             throw;
        //         }
        // ```
        // An `eEOF` ends the batch with the records read so far. A reader error (or LOSAT's
        // rejection of a Seq-id line) stops the run while the batch is read: after the
        // reports of the batches before (without the epilog) and the messages of the
        // batch's records, which are neither searched nor reported.
        if reading_ends_here {
            reading_ended = true;
            if let QueryReading::Error {
                error,
                warnings: error_warnings,
            } = std::mem::replace(&mut reading, QueryReading::End)
            {
                reports.write_before_failed_batch(outputs, done, &warnings)?;
                outputs.diagnostics.write_all(&batch_warnings)?;
                outputs.diagnostics.write_all(&error_warnings)?;
                return Err(error.into_app_error());
            }
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:378-380
        // ```c
        //     if (queries.Empty()) {
        //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
        //     }
        // ```
        // A batch without a query (a batch size of 0, a batch of skipped records only with
        // `-query_loc`, or a batch that an `eEOF` ends before its first record) fails after
        // the reports of the batches before (without the epilog), when NCBI has read that
        // batch (its records' messages).
        if start == end {
            reports.write_before_failed_batch(outputs, done, &warnings)?;
            outputs.diagnostics.write_all(&batch_warnings)?;
            return Err(crate::cli::NativeError {
                exit: 3,
                message: "BLAST engine error: Empty CBlastQueryVector\n".to_string(),
            }
            .into());
        }
        // The query splitter reads its environment variables for every batch, before the
        // queries are set up (`check_query_split_environment`); the first batch meets an
        // error there.
        if input_start == 0 {
            if let Err(error) =
                crate::blastinput::app::check_query_split_environment("TBLASTX", true)?
            {
                reports.write_before_failed_batch(outputs, done, &warnings)?;
                outputs.diagnostics.write_all(&batch_warnings)?;
                return Err(error);
            }
        }
        let batch_queries = &query_records[start..end];
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652
        // ```c
        //         } catch (const CException& e) {
        //             ...
        //             CRef<CSearchMessage> m
        //                 (new CSearchMessage(eBlastSevWarning, index, e.GetMsg()));
        //             messages[index].push_back(m);
        //             s_InvalidateQueryContexts(qinfo, index);
        //         }
        //     ...
        //     // Validate that at least one query context is valid
        //     if (BlastSetup_Validate(qinfo, NULL) != 0 && messages.HasMessages()) {
        //         NCBI_THROW(CBlastException, eSetup, messages.ToString());
        //     }
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_aux.cpp:1013-1025
        // ```c
        // TSearchMessages::ToString() const
        // {
        //     string retval;
        //     ITERATE(vector<TQueryMessages>, qm, *this) {
        //         if (qm->empty()) {
        //             continue;
        //         }
        //         ITERATE(TQueryMessages, msg, *qm) {
        //             retval += (*msg)->GetMessage() + " ";
        //         }
        //     }
        //     return retval;
        // }
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_util.c:923-929
        // ```c
        // BLAST_GetTranslatedProteinLength(size_t nucleotide_length, unsigned int context)
        // {
        //     if (nucleotide_length == 0 || nucleotide_length <= context % CODON_LENGTH) {
        //         return 0;
        //     }
        //     return (nucleotide_length - context % CODON_LENGTH) / CODON_LENGTH;
        // }
        // ```
        // A query's contexts are its six frames, of length 0 (invalid,
        // `s_QueryInfo_SetContext`) below one codon. A query without letters
        // (`Sequence contains no data`, blast_setup.hpp:190-198) gets a message and its
        // contexts are invalidated; a query of one or two letters gets none. A batch whose
        // queries all lack a valid context stops NCBI only when one of them has a message:
        // one message per query without letters (`GetMessage` puts the severity first),
        // after the reports of the batches before (`CATCH_ALL`: `BLAST engine error: `,
        // exit 3). Any other batch is set up.
        let no_data_queries = batch_queries
            .iter()
            .filter(|record| record.seq().is_empty())
            .count();
        if no_data_queries > 0
            && batch_queries
                .iter()
                .all(|record| record.seq().len() < CODON_LENGTH)
        {
            reports.write_before_failed_batch(outputs, done, &warnings)?;
            outputs.diagnostics.write_all(&batch_warnings)?;
            let mut message = String::from("BLAST engine error: ");
            for _ in 0..no_data_queries {
                message.push_str("Warning: ");
                message.push_str(NO_DATA_MESSAGE);
                message.push(' ');
            }
            message.push('\n');
            return Err(crate::cli::NativeError { exit: 3, message }.into());
        }
        // NCBI reference: c++/src/algo/blast/api/setup_factory.cpp:158-182 (the score block of
        // the batch: a name that the matrix loader does not know fails with -1,
        // blast_stat.c:1679-1680, posted for all queries, blast_setup.c:498-501)
        // ```c
        //     Int2 status = BLAST_MainSetUp(opts_memento->m_ProgramType,
        // ...
        //     Blast_Message2TSearchMessages(blast_msg.Get(), query_info, search_messages);
        //     if (status != 0 &&
        //         (!(blast_msg.Get()) || (blast_msg.Get() && blast_msg.Get()->severity == eBlastSevError)))
        //     {
        //         	string msg;
        //         	if (search_messages.HasMessages()) {
        //             		msg = search_messages.ToString();
        // ```
        // NCBI reference: c++/src/algo/blast/core/blast_message.c:246 (the text of -1)
        // ```c
        //             snprintf(buf, sizeof(buf) - 1, "Unknown error code %d", error_code);
        // ```
        // One message per query of the batch (`TSearchMessages::ToString` above, after
        // `RemoveDuplicates` in each query's list), `Error: ` first (`GetMessage`).
        // A query without letters keeps its own message (`Sequence contains no data`,
        // blast_setup_cxx.cpp:631-638, id = the query's index) after the error (id
        // kBlastMessageNoContext = -1), which `RemoveDuplicates` sorts first
        // (`CSearchMessage::operator<`, blast_types.hpp:295-304).
        if args.scoring_matrix() == TblastxMatrix::Unknown {
            reports.write_before_failed_batch(outputs, done, &warnings)?;
            outputs.diagnostics.write_all(&batch_warnings)?;
            let mut message = String::from("BLAST engine error: ");
            for record in batch_queries {
                message.push_str("Error: Unknown error code -1 ");
                if record.seq().is_empty() {
                    message.push_str("Warning: ");
                    message.push_str(NO_DATA_MESSAGE);
                    message.push(' ');
                }
            }
            message.push('\n');
            return Err(crate::cli::NativeError { exit: 3, message }.into());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:185-218 (the
        // BLAST XML report data of the first query that NCBI formats: the matrix of its
        // formatting table, found by the name as typed with `strcmp`)
        // ```c++
        //     else {
        //         string prog_name = Blast_ProgramNameFromType(
        //                                            m_Options.GetProgramType());
        //         if (prog_name != "blastn" && prog_name != "megablast") {
        //             NCBI_THROW(blast::CBlastException, eInvalidArgument,
        //                         "unsupported score matrix");
        // ```
        // NCBI searches the first batch and writes nothing of it (the XML2 and JSON documents
        // only their start, `PrintProlog`); the search does not change the output.
        if outputs
            .formats
            .iter()
            .any(|format| report::output_format(format.outfmt).is_structured_xml())
            && crate::report::blast_xml::matrix::FormattingMatrix::new(args.matrix_name()).is_none()
        {
            reports.write_before_failed_batch(outputs, done, &warnings)?;
            outputs.diagnostics.write_all(&batch_warnings)?;
            return Err(crate::cli::NativeError {
                exit: 3,
                message: "BLAST engine error: unsupported score matrix\n".to_string(),
            }
            .into());
        }
        warnings[start].extend_from_slice(&batch_warnings);
        let batch = match take(start..end) {
            Ok(batch) => {
                // LOSAT_TIMING: the batch's block, kept in its result by a search of the pipe
                // (LOSAT_X_TBXPIPE), in the batch order.
                for line in &batch.x_timing {
                    eprintln!("{line}");
                }
                batch
            }
            Err(error) if error.is::<SubjectsUnreadable>() => {
                // The first batch with a valid context of a `-db` search whose sequences could
                // not be read (`search`): NCBI's read error is an engine error after the
                // reports of the batches before; LOSAT's refusal (a truncated sequence file,
                // unreadable ambiguity data) stops the run, as any other refusal does.
                let error = stores_error
                    .take()
                    .expect("unreadable subjects keep their read error");
                let read_error =
                    error.downcast::<crate::algorithm::common::database::PreliminaryReadError>()?;
                reports.write_before_failed_batch(outputs, done, &warnings)?;
                outputs.diagnostics.write_all(&batch_warnings)?;
                return Err(crate::blastinput::app::engine_error(read_error.0));
            }
            Err(error) => match error.downcast::<ErrorAfterBatchReports>() {
                // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:203-209
                // ```c
                //                 ITERATE(CSearchResultSet, result, *results) {
                //                     formatter.PrintOneResultSet(**result, query_batch);
                //                 }
                //             }
                //         }
                //
                //         formatter.PrintEpilog(opt);
                // ```
                // The reports of the batches before are written (for the first batch, the
                // outfmt 0 prolog), then the messages of the failing batch's records (read
                // before its search), then the error, which skips `PrintEpilog`.
                Ok(ErrorAfterBatchReports(error)) => {
                    reports.write_before_failed_batch(outputs, done, &warnings)?;
                    outputs.diagnostics.write_all(&batch_warnings)?;
                    return Err(error);
                }
                Err(error) => return Err(error),
            },
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:283-293
        // ```c++
        //     m_TbackSearch.Reset(new CBlastTracebackSearch(m_QueryFactory,
        //                                                   m_InternalData,
        //                                                   m_Opts,
        //                                                   seqinfo_src,
        //                                                   search_msgs));
        // ```
        // `-db`: the definition lines of the batch's reported sequences are read where NCBI
        // makes the batch's Seq-aligns (`database::read_titles`, in report order); a header file
        // that NCBI cannot read fails this batch, after the reports of the batches before.
        if let ShownSubjects::Database { database, .. } = &subjects.shown {
            let ordered = batch.hits.as_slice();
            let read = super::super::database::read_titles(
                database,
                ordered,
                start,
                &mut reports.titles.borrow_mut(),
            );
            if let Err(error) = read {
                match error.downcast::<crate::algorithm::common::database::HeaderReadError>() {
                    Ok(failure) => {
                        reports.write_before_failed_batch(outputs, done, &warnings)?;
                        outputs.diagnostics.write_all(&batch_warnings)?;
                        return Err(failure.error);
                    }
                    Err(error) => return Err(error),
                }
            }
            // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:197-205
            // ```c++
            //             if (fmt_args->ArchiveFormatRequested(args)) {
            //                 formatter.WriteArchive(*queries, *opts_hndl, *results, 0, bah.GetMessages());
            //                 bah.ResetMessages();
            //             } else {
            //                 BlastFormatter_PreFetchSequenceData(*results, scope,
            //                 		                            fmt_args->GetFormattedOutputChoice());
            //                 ITERATE(CSearchResultSet, result, *results) {
            //                     formatter.PrintOneResultSet(**result, query_batch);
            //                 }
            // ```
            // Before the reports of a batch NCBI fetches every subject that the batch reports
            // (the sequences of `ordered`), and with `taxdb.bti` on the search path the taxonomy
            // database writes its messages there (`database::prefetch_messages`: the texts of
            // opening it at the run's first lookup, a volume's lookup of a taxid once per run):
            // after the reader's messages of the batch and before the warnings of its first
            // query, which are added below. They go with the warnings, where the first format is
            // outfmt 0, 5, 6 or 7 (`prefetches`), so every writer of the reports writes them
            // there (the reports after the loop, the error paths that write the reports of the
            // batches before). Without `taxdb.bti` nothing is read. A refusal (a taxonomy
            // database or a taxid that LOSAT cannot read as NCBI does) stops the run before any
            // report, as the refusals of the titles above and XML2's lookups (`read_taxids`) do.
            if prefetches && database.taxdb_found() {
                let messages = super::super::database::prefetch_messages(
                    database,
                    ordered,
                    &reports.titles.borrow(),
                )?;
                if plural_names && !ordered.is_empty() {
                    database.open_for_tabular_rows()?;
                }
                warnings[start].extend(messages);
            }
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1450-1452
        // ```c
        //     if (results.HasWarnings()) {
        //         ERR_POST(Warning << results.GetWarningStrings());
        //     }
        // ```
        // NCBI reference: c++/src/algo/blast/core/blast_stat.c:2815-2821
        // ```c
        //    if (valid_context == FALSE)
        //    {   /* No valid contexts were found. */
        //        /* Message for non-translated search issued above. */
        //        if (Blast_QueryIsTranslated(program) ) {
        //             Blast_MessageWrite(blast_message, eBlastSevWarning, kBlastMessageNoContext,
        //             kBlastErrMsg_CantCalculateUngappedKAParams);
        //        }
        // ```
        // NCBI reference: c++/src/algo/blast/api/local_blast.cpp:197-201,221
        // ```c
        //          for (index=0; index<local_query_data->GetNumQueries(); index++)
        //          {
        //               CConstRef<objects::CSeq_id> query_id(local_query_data->GetSeq_loc(index)->GetId());
        //               TQueryMessages q_msg;
        //               local_query_data->GetQueryMessages(index, q_msg);
        //     ...
        //          msg_vec.Combine(m_PrelimSearch->GetSearchMessages());
        // ```
        // The warnings of a query come with its report (`PrintOneResultSet`), before its
        // preamble (`QueryWarnings`): a query without letters has the message of its set-up
        // (above); a translated query has no Karlin-Altschul message of its own, and the
        // queries of a batch that is not searched (no valid context) all get the message of
        // the batch.
        for (offset, query) in batch_queries.iter().enumerate().filter(|_| result_warnings) {
            let index = start + offset;
            let mut messages = Vec::new();
            if !batch.searched {
                messages.push(crate::report::query_warnings::INVALID_QUERY_MESSAGE.to_string());
            }
            if query.seq().is_empty() {
                messages.push(NO_DATA_MESSAGE.to_string());
            }
            warnings[index].extend(crate::report::query_warnings::query_warning(
                "tblastx",
                input.ordinal(index),
                query,
                &messages,
            ));
        }
        if stream {
            // The batch's reports, written now, as NCBI's loop prints a batch's result sets
            // (`report_stream`), from its lists numbered from its first query as it made them.
            // Its lists and its queries' warnings are not kept.
            debug_assert_eq!(reports.queries_written.get(), start);
            let count = end - start;
            let batch_done = SearchedBatches {
                hits: batch.hits,
                query_stats: batch.queries,
                searched: vec![batch.searched; count],
                spans: vec![0..count],
            };
            reports.write(outputs, batch_done, &warnings, false)?;
            for query_warnings in &mut warnings[start..end] {
                *query_warnings = Vec::new();
            }
        } else {
            // The batch numbers its queries from 0.
            done.hits.extend(batch.hits.into_iter().map(|mut hsp| {
                hsp.hit.q_idx += start as u32;
                hsp
            }));
            done.query_stats.extend(batch.queries);
            done.searched
                .extend(std::iter::repeat_n(batch.searched, end - start));
            done.spans.push(start..end);
        }
        input_start = input_end;
    }
    reports.write(outputs, done, &warnings, true)
}

/// NCBI's `CODON_LENGTH`: a query shorter than one codon has no frame of positive length.
///
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:62-64
/// ```c
/// #ifndef CODON_LENGTH
/// #define CODON_LENGTH 3
/// #endif
/// ```
const CODON_LENGTH: usize = 3;

/// The results of the query batches searched so far, by query.
#[derive(Default)]
struct SearchedBatches {
    hits: Vec<TblastxHsp>,
    query_stats: Vec<TblastxQueryStats>,
    searched: Vec<bool>,
    /// The queries of each batch searched, in order: the batches of the BLAST archive (outfmt 11).
    spans: Vec<std::ops::Range<usize>>,
}

/// What the reports of the batches searched so far need besides their results.
struct BatchReports<'a> {
    args: &'a TblastxArgs,
    query_runs: &'a QueryRuns<'a>,
    query_ids: &'a [Arc<[u8]>],
    subjects: &'a ShownSubjects<'a>,
    numbers: &'a SubjectSetNumbers,
    /// `-db`: the stored titles of the batches' reported sequences (`run_in_pool`).
    titles: std::cell::RefCell<DatabaseTitles>,
    v1_titles: bool,
    /// EXPERIMENT (LOSAT_X_TBXRPTPAR): the search's pool, for the report stage
    /// (`TblastxReportRun::x_pool`).
    x_search_pool: &'a crate::utils::threading::SearchPool<'a>,
    /// The queries whose reports were written: those of the batches that the loop wrote when it
    /// took them (`report_stream`), which the next call follows.
    queries_written: std::cell::Cell<usize>,
}

impl BatchReports<'_> {
    /// Writes the reports of the queries of the batches searched (`done`), with NCBI's
    /// epilog (`PrintEpilog`) when `epilog`. `done` numbers its queries from the first query
    /// not yet written (`queries_written`: 0 unless the loop writes each batch's reports,
    /// `report_stream`). `warnings` are by query (at least those of `done`).
    fn write(
        &self,
        outputs: &mut ReportOutputs<'_>,
        done: SearchedBatches,
        warnings: &[Vec<u8>],
        epilog: bool,
    ) -> Result<()> {
        // Diagnostic (LOSAT_TIMING): the report for the common line (utils/xtiming.rs); written
        // after the batches are searched. No NCBI counterpart.
        let t_output = crate::utils::xtiming::clock();
        let first = self.queries_written.get();
        let reported = first..first + done.searched.len();
        let queries = &self.query_runs.records[reported.clone()];
        // The queries' input and ranges, numbered from `first` as `done` (`report_stream`).
        let batch_input;
        let batch_runs;
        let ranges = if first == 0 {
            self.query_runs
        } else {
            batch_input = report_stream::batch_input(self.query_runs.input, reported.clone());
            batch_runs = QueryRuns {
                input_records: self.query_runs.input_records,
                records: queries,
                input: &batch_input,
                subject_placements: self.query_runs.subject_placements,
            };
            &batch_runs
        };
        // The BLAST archive (outfmt 11) writes the Bioseq of every query of the batches NCBI
        // completes and, with `-subject`, of every subject in each archive: a title it cannot
        // write is refused before the first byte, exactly where NCBI would write that Bioseq.
        // NCBI writes none of them when no batch is completed (an error archive or nothing).
        if self.args.archive.is_some() && !reported.is_empty() {
            super::super::archive::check_titles(queries, "query")?;
            if let ShownSubjects::Fasta { records, .. } = self.subjects {
                super::super::archive::check_titles(records, "subject")?;
            }
        }
        let mut hits = done.hits;
        // `-db`: the reported sequences by rank (their records hold the stored titles, their
        // letters are read from the stored codes).
        let reported_subjects;
        let (subjects, subject_letters, database) = match self.subjects {
            ShownSubjects::Fasta { records, x_stored } => (
                *records,
                match x_stored {
                    Some((stores, ordinals)) => report::SubjectLetters::Stored { stores, ordinals },
                    None => report::SubjectLetters::Records(records),
                },
                None,
            ),
            ShownSubjects::Database {
                database,
                pairwise,
                stores,
            } => {
                let formats: Vec<report::TblastxOutputFormat> = outputs
                    .formats
                    .iter()
                    .map(|format| report::output_format(format.outfmt))
                    .collect();
                let pairwise_shown = formats.contains(&report::TblastxOutputFormat::Pairwise);
                let shown = crate::algorithm::common::database::ParsedShown {
                    tabular: formats.iter().any(|format| format.is_tabular()),
                    pairwise: pairwise_shown,
                    joined: pairwise_shown
                        || formats.iter().any(|format| format.is_structured_xml()),
                    sam: false,
                };
                reported_subjects = super::super::database::reported_subjects(
                    database,
                    &self.titles.borrow(),
                    &mut hits,
                    stores,
                    shown,
                )?;
                (
                    reported_subjects.records.as_slice(),
                    report::SubjectLetters::Stored {
                        stores,
                        ordinals: &reported_subjects.ordinals,
                    },
                    Some(ReportDatabase {
                        name: database.name(),
                        pairwise,
                        reported: &reported_subjects,
                    }),
                )
            }
        };
        // `-db`, XML2 and JSON: the taxid of each reported sequence, and what the taxonomy
        // database says about it (the refusals of both come before any byte).
        let xml_taxa = match (&self.subjects, &database) {
            (ShownSubjects::Database { database, .. }, Some(shown))
                if outputs.formats.iter().any(|format| {
                    matches!(
                        report::output_format(format.outfmt),
                        report::TblastxOutputFormat::Xml2(_)
                    )
                }) =>
            {
                database.read_taxids(shown.reported)?
            }
            _ => Default::default(),
        };
        // `-db`, a custom tabular report that shows `staxid`, `staxids` or a taxonomy name: the
        // taxonomy of each reported sequence (a header that does not decode is refused before
        // any byte, `tabular`).
        let tabular_taxonomy = match (&self.subjects, &database) {
            (ShownSubjects::Database { database, .. }, Some(shown))
                if outputs.formats.iter().any(|format| {
                    matches!(
                        report::output_format(format.outfmt),
                        report::TblastxOutputFormat::CustomTabular(spec) if spec.needs().taxonomy
                    )
                }) =>
            {
                database.tabular_taxonomy(shown.reported)?
            }
            _ => Vec::new(),
        };
        // And the names of the taxonomy database that its rows look up (with their refusals,
        // before any byte).
        let custom_tabular = outputs
            .formats
            .iter()
            .find_map(|format| match report::output_format(format.outfmt) {
                report::TblastxOutputFormat::CustomTabular(spec) => Some(spec),
                _ => None,
            });
        let tabular_names = match (&self.subjects, custom_tabular) {
            (ShownSubjects::Database { database, .. }, Some(spec)) => {
                database.tabular_names(spec, &tabular_taxonomy)?
            }
            _ => None,
        };
        let written = write_tblastx_outputs(
            hits,
            outputs,
            &TblastxReportRun {
                args: self.args,
                queries,
                subjects,
                subject_letters,
                database,
                query_ids: &self.query_ids[reported.clone()],
                query_stats: &done.query_stats,
                searched: &done.searched,
                warnings: &warnings[reported.clone()],
                ranges,
                epilog,
                queries_before: first,
                numbers: self.numbers,
                xml_taxa: &xml_taxa,
                tabular_taxonomy: &tabular_taxonomy,
                tabular_names: tabular_names.as_ref(),
                v1_titles: self.v1_titles,
                archive_spans: &done.spans,
                x_pool: crate::algorithm::tblastx::x_rpt_par::pool(self.x_search_pool),
                x_pairwise_pool: crate::report::x_pairwise_par::pool(self.x_search_pool),
            },
        );
        crate::utils::xtiming::add(crate::utils::xtiming::Phase::Output, t_output);
        self.queries_written.set(reported.end);
        written
    }

    /// Writes what NCBI has written when a query batch fails, before its error: the reports
    /// of the batches before (`done`; without the epilog), or, for the first batch, the
    /// outfmt 0 prolog of the sinks that write it with the reports (the others wrote it
    /// before the search, `write_pairwise_prologs`). Where the loop wrote the batches before
    /// when it took them (`report_stream`), nothing remains.
    fn write_before_failed_batch(
        &self,
        outputs: &mut ReportOutputs<'_>,
        done: SearchedBatches,
        warnings: &[Vec<u8>],
    ) -> Result<()> {
        if self.queries_written.get() > 0 {
            Ok(())
        } else if done.searched.is_empty() {
            write_pairwise_prologs(outputs, self.args, self.subjects, false)?;
            // The XML2 and JSON documents start with their prolog too, and outfmt 20 with its line
            // of the field names (and its formatter's warning, `tabular::write_tabular_prolog`).
            for format in outputs.formats.iter_mut() {
                match report::output_format(format.outfmt) {
                    report::TblastxOutputFormat::Xml2(kind) => {
                        let mut writer = format.sink.open()?;
                        writer
                            .write_all(&crate::report::blast_xml2::framing::prolog_bytes(kind))?;
                        writer.flush()?;
                    }
                    report::TblastxOutputFormat::CustomTabular(spec) => {
                        let mut writer = format.sink.open()?;
                        super::super::tabular::write_tabular_prolog(
                            spec,
                            self.args.tabular_taxdb.unwrap_or(false),
                            &mut writer,
                            &mut *outputs.diagnostics,
                        )?;
                        writer.flush()?;
                    }
                    _ => {}
                }
            }
            Ok(())
        } else {
            self.write(outputs, done, warnings, false)
        }
    }
}

/// An error of a query batch's search that NCBI raises after the reports of the batches
/// before (`run_in_pool`): the average subject length (`search_query_batch`).
#[derive(Debug)]
struct ErrorAfterBatchReports(anyhow::Error);

impl std::fmt::Display for ErrorAfterBatchReports {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for ErrorAfterBatchReports {}

fn search_query_batch<S: TbxScore>(
    // The scores of the search's matrix (`scoring.rs`).
    sc: S,
    args: &TblastxArgs,
    query_records: &[FastaRecord],
    subjects: Subjects<'_>,
    numbers: &SubjectSetNumbers,
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
    // EXPERIMENT (LOSAT_X_TBXPIPE): the batch pipe's view of this search's seed stage, `None`
    // without the pipe (tblastx/x_batch_pipe.rs).
    x_seed_end: Option<&x_batch_pipe::XSeedEnd<'_>>,
) -> Result<TblastxBatch> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:575-582
    // ```c
    // score = s_BlastAaExtendTwoHit(matrix, subject, query,
    //     last_hit + wordsize, subject_offset, query_offset,
    //     cutoffs->x_dropoff, &hsp_q, &hsp_s, &hsp_len, use_pssm,
    //     wordsize, &right_extend, &s_last_off);
    // ```
    // LOSAT-only existing diagnostics: capture configuration once per search;
    // the NCBI extension inputs and all biological decisions stay unchanged.
    let extension_debug_enabled = crate::utils::xdiag::var("LOSAT_DEBUG_EXTENSION").is_ok();
    // Optional timing breakdown (disabled by default to preserve output/parity logs)
    let timing_enabled = std::env::var_os("LOSAT_TIMING").is_some();
    let t_total = Instant::now();
    let mut t_read_queries = Duration::ZERO;
    let mut t_build_lookup = Duration::ZERO;
    let mut t_read_subjects = Duration::ZERO;

    // NCBI reference: ncbi-blast/c++/include/algo/blast/blastinput/blast_args.hpp:1290-1296
    // ```c
    // CMTArgs(...)
    // {
    // #ifdef NCBI_NO_THREADS
    //     m_NumThreads = CThreadable::kMinNumThreads;
    //     m_MTMode = eNotSupported;
    // #endif
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
    // ```c
    // const int kMaxValue = static_cast<int>(CSystemInfo::GetCpuCount());
    // ...
    // int num_threads = args[kArgNumThreads].AsInteger();
    // if (num_threads > kMaxValue) {
    //     m_NumThreads = kMaxValue;
    // } else {
    //     m_NumThreads = num_threads;
    // }
    // ```
    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // TBlastThreads the_threads(GetNumberOfThreads());
    // (*thread)->Run(); (*thread)->Join(&result);
    let num_threads = parallel_pool.threads();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:82-88
    // ```c
    // if (num_threads > 1) {
    //     SetNumberOfThreads(num_threads);
    // }
    // ```
    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // TBlastThreads the_threads(GetNumberOfThreads());
    // (*thread)->Run(); (*thread)->Join(&result);
    let requested_parallel = parallel_pool.enabled();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1002-1003
    // ```c
    // status = s_BlastSetUpAuxStructures(..., aux_struct);
    // ```
    //
    // The scan-chunk experiment keeps the canonical two-hit state serial inside
    // each NCBI WordFinder search unit, so subject-level parallelism is disabled
    // while this gate is active.
    let use_serial_scan_chunks = std::env::var_os("LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS").is_some();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:82-88
    // ```c
    // if (num_threads > 1) {
    //     SetNumberOfThreads(num_threads);
    // }
    // ```
    let query_code = GeneticCode::from_id(args.query_gencode);
    let db_code = GeneticCode::from_id(args.db_gencode);
    // LOSAT intentionally treats `--db-gencode` as the subject genetic code for
    // local `-s/--subject` searches. Unlike BLAST+, search/scoring and reporting
    // both honor the explicit subject code.

    // [C] window = diag->window;
    let window = args.window();
    // [C] wordsize = lookup->word_length;
    let wordsize: i32 = 3;

    // NCBI BLAST computes x_dropoff per-context using kbp[context]->Lambda:
    //   p->cutoffs[context].x_dropoff_init =
    //       (Int4)(sbp->scale_factor * ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:219-221
    //
    // For translated queries, NCBI computes kbp_std per context and applies check_ideal:
    //   if (check_ideal && kbp->Lambda >= sbp->kbp_ideal->Lambda)
    //      Blast_KarlinBlkCopy(kbp, sbp->kbp_ideal);
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2797
    // We still maintain per-context structure for parity.
    //
    // x_dropoff_per_context is populated after build_ncbi_lookup() creates the contexts.
    let ungapped_params_for_xdrop = lookup_protein_params_ungapped(sc.matrix());

    let diag_enabled = diagnostics_enabled();
    let debug_cutoffs_all = crate::utils::xdiag::var_os("LOSAT_DEBUG_CUTOFFS_ALL").is_some();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3265-3273
    // ```c
    // #if _BLAST_DEBUG
    // arg_desc.AddFlag("verbose", "Produce verbose output (show BLAST options)",
    //                  true);
    // arg_desc.AddFlag("remote_verbose",
    //                  "Produce verbose output for remote searches", true);
    // #endif /* _BLAST_DEBUG */
    // ```
    let debug_output_filter = std::env::var_os("LOSAT_DEBUG_OUTPUT_FILTER").is_some();
    let debug_hsp_saving = crate::utils::xdiag::var_os("LOSAT_DEBUG_HSP_SAVING").is_some();
    let diagnostics = std::sync::Arc::new(DiagnosticCounters::default());

    // Debug: optional scan output dump around a specific subject offset.
    // The offset_pairs correspond to s_BlastAaScanSubject copy loop.
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:99-115
    let scan_debug_center = crate::utils::xdiag::var("LOSAT_DEBUG_SCAN_SOFF")
        .ok()
        .and_then(|v| v.parse::<i32>().ok());
    let scan_debug_window = crate::utils::xdiag::var("LOSAT_DEBUG_SCAN_WINDOW")
        .ok()
        .and_then(|v| v.parse::<i32>().ok())
        .unwrap_or(0);
    let scan_debug_range =
        scan_debug_center.map(|center| (center - scan_debug_window, center + scan_debug_window));
    if let Some((lo, hi)) = scan_debug_range {
        eprintln!(
            "[DEBUG SCAN_OFF] enabled s_off_range=[{},{}] (center={} window={})",
            lo,
            hi,
            scan_debug_center.unwrap_or(0),
            scan_debug_window
        );
    }

    let t_phase_read_queries = Instant::now();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:486-651
    // ```c
    // SetupQueries_OMF(IBlastQuerySource& queries,
    //                  BlastQueryInfo* qinfo,
    //                  BLAST_SequenceBlk** seqblk,
    //                  EBlastProgramType prog,
    //                  ...)
    // ```
    let queries_raw: &[FastaRecord] = query_records;
    // NCBI tblastx low-complexity filtering uses SEG on translated protein sequences.
    // No nucleotide-level DUST masking is applied.
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

    // NCBI reference (translate all 6 frames for tblastx queries):
    // ncbi-blast/c++/src/algo/blast/core/blast_util.c:1076-1101
    // ```c
    // frame_offsets = (Uint4*) malloc((NUM_FRAMES+1)*sizeof(Uint4));
    // frame_offsets[0] = 0;
    // for (context = 0; context < NUM_FRAMES; ++context) {
    //    frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
    //    length = BLAST_GetTranslation(nucl_seq, nucl_seq_rev,
    //       nucl_length, frame, translation_buffer+offset, genetic_code);
    //    offset += length + 1;
    //    frame_offsets[context+1] = offset;
    // }
    // ```
    // Diagnostic (LOSAT_TIMING): the queries' six-frame translation for the common line
    // (utils/xtiming.rs). No NCBI counterpart.
    let t_translate = crate::utils::xtiming::clock();
    let mut query_frames: Vec<Vec<QueryFrame>> = queries_raw
        .iter()
        .map(|r| generate_frames(r.seq(), &query_code))
        .collect();
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::Translate, t_translate);

    // NCBI blast_args.cpp:404-406: opt.SetSegFilteringWindow(...);
    // opt.SetSegFilteringLocut(...); opt.SetSegFilteringHicut(...);
    if let Some(params) = args.seg_spec()?.params() {
        // Diagnostic (LOSAT_TIMING): the queries' SEG for the common line (utils/xtiming.rs). No
        // NCBI counterpart.
        let t_seg = crate::utils::xtiming::clock();
        let seg = SegMasker::new(params.window, params.locut, params.hicut);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1280-1290
        // ```c
        //     for (context = query_info->first_context;
        //          context <= query_info->last_context; ++context) {
        // ...
        //         status = s_GetFilteringLocationsForOneContext(query_blk,
        // ```
        // Dispatch point (LOSAT_X_TBXQFRAMES, strict set): with two or more threads the frames
        // are masked on the search pool, one task per frame; otherwise, and with the switch off,
        // the frames are masked one after the other as before (tblastx/x_query_frames.rs, which
        // holds the per-frame bodies `mask_frame` and `mask_frame_soft`).
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:621-624
        // ```c
        //     if (!mask_at_hash) {
        //         BlastSetUp_MaskQuery(query_blk, query_info, filter_maskloc,
        //                              program_number);
        //     }
        // ```
        // `-soft_masking true` keeps the intervals for the lookup table and the frames' letters.
        if args.soft_masking() {
            x_query_frames::mask_soft(&mut query_frames, &seg, parallel_pool);
        } else {
            x_query_frames::mask(&mut query_frames, &seg, parallel_pool);
        }
        crate::utils::xtiming::add(crate::utils::xtiming::Phase::Seg, t_seg);
    }

    t_read_queries = t_phase_read_queries.elapsed();

    let t_phase_build_lookup = Instant::now();
    // Note: karlin_params argument is now unused - computed per context in build_ncbi_lookup()
    // We still pass it for x_dropoff calculation (which uses ideal params for all contexts in tblastx)
    // NCBI lookup build always precomputes neighbors (no lazy mode).
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:446-543
    // NCBI reference (exact match indexing before neighbor expansion):
    // ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:454-456
    // ```c
    // BlastLookupIndexQueryExactMatches(exact_backbone, lookup->word_length,
    //                                   lookup->charsize, lookup->word_length,
    //                                   query, location);
    // ```
    // NCBI reference: c++/src/algo/blast/core/blast_aalookup.c:245
    // ```c
    //     lookup->threshold = (Int4)opt->threshold;
    // ```
    // The x86-64 conversion: `INT_MIN` for a threshold beyond `Int4` (every word).
    //
    // EXPERIMENT (PAIRLUT, graduated on 2026-10-22), dispatch point 0 (per batch, at the
    // lookup build): when PAIRPAR will run the seed stage of every subject of this batch
    // as units, the units' tables (one per query-context group) are built directly from the
    // query, and `lookup` becomes a stub with the full table's header and no entries
    // (x_pair_par.rs `lut_build`; lookup/x_pair_lut.rs for why the tables are the same). The
    // inputs are the ones dispatch point 1 below gives `prepare`, by the same expressions;
    // `prepare` decides again and stops the run if it disagrees. When the unit path will not run,
    // the reference builder below runs.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:93-97
    // ```c
    // BlastAaLookupTableNew(lookup_options, (BlastAaLookupTable* *)
    //                       &lookup_wrap->lut);
    // ((BlastAaLookupTable*)lookup_wrap->lut)->use_pssm = has_pssm;
    // BlastAaLookupIndexQuery( (BlastAaLookupTable*) lookup_wrap->lut, matrix,
    //                           query, lookup_segments, 0);
    // ```
    // The tables of TBXLUT, PAIRLUT and PAIRPAR hold BLOSUM62's scores: another matrix's
    // search builds NCBI's lookup table with its scores (`build_ncbi_lookup_for_profile`) and
    // runs the frame loop.
    let mut x_pair_lut_batch = if !S::BLOSUM62 {
        None
    } else {
        x_pair_par::lut_build(x_pair_par::XLutEarly {
            query_frames: &query_frames,
            threshold: crate::core::blast_util::ncbi_int4_from_double(args.word_threshold()),
            window,
            subjects,
            debugging: trace_hsp_target().is_some()
                || scan_debug_range.is_some()
                || extension_debug_enabled
                || diag_enabled
                || debug_hsp_saving,
            requested_parallel,
            use_serial_scan_chunks,
            pool: parallel_pool,
        })
    };
    let x_pair_lut_stub = x_pair_lut_batch
        .as_mut()
        .and_then(x_pair_par::XLutBatch::take_stub);
    // True when `lookup` is the stub: the reference frame loop must then run zero times
    // (asserted at dispatch point 2).
    let x_pair_lut_is_stub = x_pair_lut_stub.is_some();
    // Dispatch point (LOSAT_X_TBXLUT, on unless `LOSAT_X_TBXLUT=0`): the full table of the batch
    // from the chain builder (`build_ncbi_lookup`) or filled in place from per-cell counts; the
    // same table.
    let (lookup, contexts) = match x_pair_lut_stub {
        Some(stub) => {
            crate::algorithm::tblastx::lookup::x_tbxlut_reach_stub();
            stub
        }
        None if !S::BLOSUM62 => {
            // NCBI reference: c++/src/algo/blast/core/blast_setup.c:812-822 (the alpha and beta
            // of the length adjustment: the matrix's ungapped row, `BLAST_GetAlphaBeta`, for
            // the Karlin block of every context)
            // ```c
            //          BLAST_GetAlphaBeta(sbp->name, &alpha, &beta,
            //                             scoring_options->gapped_calculation,
            //                             scoring_options->gap_open,
            //                             scoring_options->gap_extend,
            //                             sbp->kbp_std[index]);
            //       }
            //       BLAST_ComputeLengthAdjustment(kbp->K, kbp->logK,
            //                                     alpha/kbp->Lambda, beta,
            // ```
            // (`compute_karlin_params_ungapped` gives every block BLOSUM62's.)
            let (lookup, mut contexts) =
                crate::algorithm::tblastx::lookup::build_ncbi_lookup_for_profile(
                    &query_frames,
                    crate::core::blast_util::ncbi_int4_from_double(args.word_threshold()),
                    &ungapped_params_for_xdrop,
                    true,
                    sc.matrix(),
                    3,
                );
            for ctx in contexts.iter_mut().filter(|ctx| ctx.is_valid) {
                ctx.karlin_params.alpha = ungapped_params_for_xdrop.alpha;
                ctx.karlin_params.beta = ungapped_params_for_xdrop.beta;
            }
            (lookup, contexts)
        }
        None => crate::algorithm::tblastx::lookup::x_tbx_build_ncbi_lookup(
            &query_frames,
            crate::core::blast_util::ncbi_int4_from_double(args.word_threshold()),
            &ungapped_params_for_xdrop, // Used for x_dropoff calculation only
            true,
        ),
    };
    t_build_lookup = t_phase_build_lookup.elapsed();
    // Diagnostic (LOSAT_TIMING): the same interval for the common line (utils/xtiming.rs). No
    // NCBI counterpart.
    crate::utils::xtiming::add_duration(crate::utils::xtiming::Phase::BuildLookup, t_build_lookup);

    // NCBI BLAST: word_params->cutoffs[context].x_dropoff_init
    // Compute per-context x_dropoff using kbp[context]->Lambda.
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:219-221
    //
    // For translated queries, NCBI computes kbp_std per context and applies check_ideal:
    //   if (check_ideal && kbp->Lambda >= kbp_ideal->Lambda) Blast_KarlinBlkCopy(kbp, kbp_ideal);
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2797
    //
    // NCBI BLAST dynamic x_dropoff (ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:380-383):
    //   if (curr_cutoffs->x_dropoff_init == 0)
    //      curr_cutoffs->x_dropoff = new_cutoff;  // x_dropoff = cutoff_score
    //   else
    //      curr_cutoffs->x_dropoff = curr_cutoffs->x_dropoff_init;
    //
    // For TBLASTX, x_dropoff_init is non-zero, so this dynamic update is not triggered.
    // However, we store x_dropoff_init here and apply the logic during subject processing
    // where cutoff_score is available.
    let x_dropoff_per_context: Vec<i32> = contexts
        .iter()
        .map(|ctx| x_drop_raw_score(X_DROP_UNGAPPED_BITS, &ctx.karlin_params, 1.0))
        .collect();

    let t_phase_read_subjects = Instant::now();
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
    // An unreadable database is not empty: its first batch with a valid context fails
    // (`SubjectsUnreadable`, below).
    if queries_raw.is_empty()
        || (subjects.count() == 0 && !matches!(subjects, Subjects::Unreadable))
    {
        return Ok(TblastxBatch {
            hits: Vec::new(),
            searched: true,
            queries: vec![TblastxQueryStats::default(); queries_raw.len()],
            x_timing: Vec::new(),
        });
    }
    t_read_subjects = t_phase_read_subjects.elapsed();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1475
    // ```c
    // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //        != BLAST_SEQSRC_EOF) {
    //    status = s_BlastSearchEngineCore(...);
    // }
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:452-584
    // ```c
    // while (TRUE) {
    //    status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
    //                                   dbseq_chunk_overlap);
    //    if (status == SUBJECT_SPLIT_DONE) break;
    //    ...
    //    status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list, ...);
    // }
    // ```
    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // TBlastThreads the_threads(GetNumberOfThreads());
    // (*thread)->Run(); (*thread)->Join(&result);
    let use_parallel = requested_parallel;
    let use_parallel_chunks = requested_parallel
        && !use_serial_scan_chunks
        && (std::env::var_os("LOSAT_TBLASTX_PARALLEL_CHUNKS").is_some()
            || (cfg!(all(target_arch = "wasm32", feature = "wasm-threads"))
                && subjects.lengths().any(|length| {
                    tblastx_estimated_subject_chunk_work_items(
                        length,
                        tblastx_max_dbseq_len_for_run(),
                    ) > 1
                })));
    crate::utils::threading::report_stage(
        "tblastx",
        "subjects",
        subjects.count(),
        use_parallel && !use_serial_scan_chunks && subjects.count() > 1,
    );

    // NCBI BLAST Karlin params for TBLASTX (ungapped-only algorithm):
    //
    // TBLASTX is explicitly ungapped-only (blast_options.c line 869-873):
    //   "Gapped search is not allowed for tblastx"
    //
    // For bit score and E-value calculation, NCBI uses sbp->kbp (ungapped):
    //   blast_hits.c line 1833: kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
    //   blast_hits.c line 1918: same pattern in Blast_HSPListGetBitScores
    //
    // NCBI reference (ncbi-blast/c++/src/algo/blast/core/blast_setup.c:768):
    //   kbp_ptr = (scoring_options->gapped_calculation ? sbp->kbp_gap_std : sbp->kbp);
    // For tblastx, gapped_calculation = FALSE, so kbp_ptr = sbp->kbp (ungapped params)
    // Therefore, ALL calculations (eff_searchsp, cutoff, bit score, E-value) use UNGAPPED params
    //
    // (the contexts' blocks, `QueryContext::karlin_params`).

    // Compute NCBI-style average query length for linking cutoffs
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:998-1082
    // NCBI uses average over ALL contexts including zero-length (frame restriction via strand)
    let query_nucl_lengths: Vec<usize> = queries_raw.iter().map(|r| r.seq().len()).collect();
    let avg_query_length = compute_avg_query_length_ncbi(&query_nucl_lengths);

    let t_search_start = Instant::now();
    let scan_ns = AtomicU64::new(0);
    let scan_calls = AtomicU64::new(0);
    let ungapped_ns = AtomicU64::new(0);
    let ungapped_calls = AtomicU64::new(0);
    let reeval_ns = AtomicU64::new(0);
    let reeval_calls = AtomicU64::new(0);
    let linking_ns = AtomicU64::new(0);
    let linking_calls = AtomicU64::new(0);
    let identity_ns = AtomicU64::new(0);
    let identity_calls = AtomicU64::new(0);
    // No NCBI counterpart: the stage timer of LOSAT_TIMING (x_search_timing.rs). It reads the
    // clock only when `timing_enabled`, and the search reads nothing from it.
    let x_times = x_search_timing::XSearchTimes::new(timing_enabled);
    // No NCBI counterpart: the linking timer of LOSAT_TIMING (x_link_timing.rs). It arms a
    // linking call only when `timing_enabled`, and the search reads nothing from it.
    let x_link_timer = XLinkTimer::new(timing_enabled);

    // NCBI verbose CLI flag is only present under _BLAST_DEBUG builds.
    // ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3265-3273
    // ```c
    // #if _BLAST_DEBUG
    // arg_desc.AddFlag("verbose", "Produce verbose output (show BLAST options)",
    //                  true);
    // arg_desc.AddFlag("remote_verbose",
    //                  "Produce verbose output for remote searches", true);
    // #endif /* _BLAST_DEBUG */
    // ```
    let bar = ProgressBar::hidden();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1475
    // ```c
    // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //        != BLAST_SEQSRC_EOF) {
    //    ...
    //    status = s_BlastSearchEngineCore(..., &hsp_list, ...);
    // }
    // ```
    // NCBI iterates subjects before reducing results. Only the parallel subject
    // traversal sends to this queue; a single-subject caller retains its hits
    // directly and must not keep a sender alive while joining the writer.
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let use_channel = use_parallel && !use_serial_scan_chunks && subjects.count() > 1;
    #[cfg(any(not(feature = "parallel"), target_arch = "wasm32"))]
    let use_channel = false;

    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let (tx_opt, mut rx_opt) = if use_channel {
        let (tx, rx) = channel::<(usize, Vec<TblastxHsp>)>();
        (Some(tx), Some(rx))
    } else {
        (None, None)
    };
    #[cfg(not(all(feature = "parallel", not(target_arch = "wasm32"))))]
    let tx_opt: Option<std::sync::mpsc::Sender<(usize, Vec<TblastxHsp>)>> = None;
    let evalue_threshold = args.evalue;
    // The hit lists of the batch's queries (tblastx/hit_collector.rs). Each subject's final HSPs
    // reach them as one message at the end of the subject (`emit_subject_hits`): through the
    // channel to the writer thread, or, without one (a serial search, the WASI targets), into
    // the collector the subject loop keeps here.
    let local_collector: std::sync::Mutex<Option<HitCollector>> =
        std::sync::Mutex::new(if use_channel {
            None
        } else {
            Some(HitCollector::new(args, query_records))
        });

    // Diagonal array sizing MUST match NCBI's `s_BlastDiagTableNew`:
    // it depends only on (query_length + window_size), not on subject length.
    //
    // NCBI reference (verbatim):
    //   diag_array_length = 1;
    //   while (diag_array_length < (qlen+window_size))
    //       diag_array_length = diag_array_length << 1;
    //   diag_table->diag_array_length = diag_array_length;
    //   diag_table->diag_mask = diag_array_length-1;
    // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-61
    //
    // `query_length` here must match BLAST_SequenceBlk->length, computed as
    // last_context.query_offset + last_context.query_length (excludes the final trailing NULLB).
    // References: ncbi-blast/c++/src/algo/blast/core/blast_query_info.c:311-315, 378-381
    // NCBI BLAST diag array sizing:
    // diag_array_length = next_power_of_2(qlen + window_size)
    // For tblastx, qlen = total concatenated query buffer (all 6 frames)
    // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-61
    let query_length: i32 = contexts
        .last()
        .map(|c| c.frame_base + c.aa_len as i32)
        .unwrap_or(0);
    crate::blastinput::app::check_diag_table_window(query_length, window, "TBLASTX")?;
    let mut diag_array_size: i32 = 1;
    while diag_array_size < (query_length + window) {
        diag_array_size <<= 1;
    }
    let diag_mask: i32 = diag_array_size - 1;
    if trace_hsp_target().is_some() {
        // NCBI: diag_array_length = next_power_of_2(qlen + window_size).
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-61
        eprintln!(
            "[TRACE_HSP] diag_table query_length={} window={} diag_array_size={} diag_mask={}",
            query_length, window, diag_array_size, diag_mask
        );
    }

    // [C] array_size for offset_pairs
    // NCBI: GetOffsetArraySize() = OFFSET_ARRAY_SIZE (4096) + lookup->longest_chain
    // Reference: ncbi-blast/c++/include/algo/blast/core/lookup_wrap.h + lookup_wrap.c
    const OFFSET_ARRAY_SIZE: i32 = 4096;
    let offset_array_size: i32 = OFFSET_ARRAY_SIZE + lookup.longest_chain.max(0);
    // EXPERIMENT (LOSAT_X_TBXSCANBF): with the switch, the offset arrays have X_SCAN_BF_SLACK
    // entries after `offset_array_size`.
    let offset_pairs_len =
        crate::algorithm::tblastx::x_scan_bf::x_offset_pairs_len(offset_array_size);
    // EXPERIMENT (LOSAT_X_TBXSCANBF): whether the reference frame loop scans this batch's table with
    // the branch-free loop (decided once per table).
    let x_scan_bf = crate::algorithm::tblastx::x_scan_bf::x_scan_bf_choice(&lookup);

    // EXPERIMENT (PAIRPAR), dispatch point 1 (per batch): the gates, the exactness guards (the
    // first one only without LOSAT_X_TBXBAND) and the cost rule, then one lookup table per query
    // context (each backbone chain filtered to the context, order kept). None (a gate, a guard
    // or the cost rule) keeps the reference loop for every subject of the batch.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-63
    // ```c
    // while (diag_array_length < (qlen+window_size))
    // {
    //         diag_array_length = diag_array_length << 1;
    // }
    // ...
    // diag_table->diag_mask = diag_array_length-1;
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:337-343
    // ```c
    // if (lookup->thin_backbone[i] ) {
    //     Int4 * dest = NULL;
    //     /* set the corresponding bit in the pv_array */
    //     PV_SET(pv, i, PV_ARRAY_BTS);
    //     bbc[i].num_used = lookup->thin_backbone[i][1];
    // ```
    // The guards compare the subject lengths (and, without LOSAT_X_TBXBAND, the context lengths)
    // with the slack of this table (size above minus query_length); the per-context tables use the
    // layout above.
    let x_pair_par_batch = if !S::BLOSUM62 {
        None
    } else {
        x_pair_par::XPairPar::prepare(x_pair_par::XBatch {
            lookup: &lookup,
            contexts: &contexts,
            query_length,
            diag_array_size,
            subjects,
            debugging: trace_hsp_target().is_some()
                || scan_debug_range.is_some()
                || extension_debug_enabled
                || diag_enabled
                || debug_hsp_saving,
            chunk_switches: use_serial_scan_chunks || use_parallel_chunks,
            // PAIRLUT: the groups built at dispatch point 0.
            lut: x_pair_lut_batch,
            // LOSAT_X_PAIRPARGATE: the thread count and threshold dispatch point 0 decided with.
            threads: parallel_pool.threads(),
            threshold: crate::core::blast_util::ncbi_int4_from_double(args.word_threshold()),
        })
    };
    // The reference frame loop scans the batch's table (LOSAT_X_TBXSCANBF's choice above).
    if x_pair_par_batch.is_none() {
        crate::algorithm::tblastx::x_scan_bf::x_reach_scan_bf(x_scan_bf);
    }
    // Dispatch point 1 of the seed release (per batch): the batch's tables move
    // into two holders, `x_seed_batch.lookup` (the lookup table) and `x_seed_batch.units` (the
    // unit path with its group tables). From here on they are reached only through the stages
    // that each call of `process_subject` opens (dispatch point 2), and the
    // batch's last units stage frees the unit path (dispatch point 3) and its last lookup stage
    // the lookup table (dispatch point 5). `x_worker_diag_cells` is the size the workers'
    // diagonal tables are created with: `diag_array_size`, or 0 when the
    // unit path takes every subject, so that the reference frame loop, the
    // only reader of that table, runs for none.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1002-1004
    // ```c
    // if ((status = BlastExtendWordNew(query->length, word_params,
    //                                 &aux_struct->ewp)) != 0)
    //   return status;
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1604-1605
    // ```c
    // word_params = BlastInitialWordParametersFree(word_params);
    // s_BlastCoreAuxStructFree(aux_struct);
    // ```
    // NCBI keeps its lookup table and diagonal table through the subject loop; when LOSAT
    // frees its copies is not a value NCBI computes (x_seed_release.rs).
    let x_seed_batch = x_seed_release::XSeedBatch::new(lookup, x_pair_par_batch, subjects.count());
    let x_worker_diag_cells = x_seed_batch.worker_table_cells(diag_array_size as usize);

    let contexts_ref = &contexts;
    // EXPERIMENT (LOSAT_X_TBXUNGAPSENT): the reference frame loops' extensions may take the walks
    // the sentinels stop (not with the debug output of `extend_hit_two_hit`). The sentinel walks
    // hold BLOSUM62's byte table: another matrix takes `extend_hit_two_hit_with` with its scores.
    let tbx_sent = S::BLOSUM62
        && !extension_debug_enabled
        && x_tbx_sent_contexts(contexts.iter().map(|c| c.aa_seq.as_slice()));

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1475
    // ```c
    // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //        != BLAST_SEQSRC_EOF) {
    //    ...
    //    status = s_BlastSearchEngineCore(...);
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1475
    // ```c
    // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //        != BLAST_SEQSRC_EOF) {
    //    ...
    //    status = s_BlastSearchEngineCore(...);
    // }
    // ```
    fn for_each_subjects<FInit, FBody, S>(subjects: usize, init: FInit, mut body: FBody)
    where
        FInit: FnOnce() -> WorkerState<S>,
        FBody: FnMut(&mut WorkerState<S>, usize),
    {
        let mut state = init();
        for s_idx in 0..subjects {
            body(&mut state, s_idx);
        }
    }

    // NCBI reference: c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-180,261-264;
    // c++/src/algo/blast/core/blast_engine.c:1372-1374,1434-1443
    // m_iTotalLength += (Int8) (*iter)->length;
    // avg_length = (Uint4) (total_length / num_seqs);
    // BlastInitialWordParametersUpdate(..., BlastSeqSrcGetAvgSeqLen(seq_src), ...);
    // if (db_length == 0) { BLAST_OneSubjectUpdateParameters(...); }
    // The set's numbers are read once per run (`SubjectSetNumbers`: the records' sums, or the
    // database index's numbers).
    let db_length_nucl = numbers.total;
    let db_num_seqs = numbers.count;
    // NCBI reference: c++/src/algo/blast/core/blast_setup.c:535-560
    // if (valid_context_found) { return 0; } else { return 1; }
    // NCBI reference: c++/src/algo/blast/api/local_blast.cpp:177-178,206-208
    // ```c
    //     int status = m_PrelimSearch->CheckInternalData();
    //     if (status != 0)
    //     ...
    //             pair<double, double> tmp_pair(-1.0, -1.0);
    //             CRef<CBlastAncillaryData>  tmp_ancillary_data(new CBlastAncillaryData(tmp_pair, tmp_pair, tmp_pair, 0));
    // ```
    // A batch without a valid context is not searched; each of its queries is reported
    // with the `-1` statistics.
    if !contexts_ref.iter().any(|ctx| ctx.is_valid) {
        return Ok(TblastxBatch {
            hits: Vec::new(),
            searched: false,
            queries: vec![TblastxQueryStats::default(); queries_raw.len()],
            x_timing: Vec::new(),
        });
    }
    // `-db` whose sequences could not be read: this first batch with a valid context fails
    // (`run_in_pool` raises the kept read error).
    if matches!(subjects, Subjects::Unreadable) {
        return Err(SubjectsUnreadable.into());
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:521-528
    // ```c
    //        } else {
    //           ASSERT(sbp->kbp_gap == NULL);
    //           /* for ungapped cases we do not have gbp filled */
    //           if (sbp->gbp) {
    //               sfree(sbp->gbp);
    //               sbp->gbp=NULL;
    //           }
    //        }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:908-939
    // ```c
    //    if (seq_src) {
    //       total_length = BlastSeqSrcGetTotLenStats(seq_src);
    //       if (total_length <= 0)
    //           total_length = BlastSeqSrcGetTotLen(seq_src);
    //       ...
    //       if (total_length > 0) {
    //           num_seqs = BlastSeqSrcGetNumSeqsStats(seq_src);
    //           if (num_seqs <= 0)
    //               num_seqs = BlastSeqSrcGetNumSeqs(seq_src);
    //       } else {
    //           /* Not a database search; each subject sequence is considered
    //              individually */
    //           Int4 oid = 0;  /* Get length of first sequence. */
    //           if ( (total_length = BlastSeqSrcGetSeqLen(seq_src, (void*) &oid)) < 0) {
    //               total_length = -1;
    //               num_seqs = -1;
    //           }
    //           num_seqs = 1;
    //       }
    //    }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:969-980
    // ```c
    //    if (sbp->gbp) {
    //        min_subject_length = BlastSeqSrcGetMinSeqLen(seq_src);
    //        if (Blast_SubjectIsTranslated(program_number)) {
    //            min_subject_length/=3;
    //        }
    //    } else {
    //        min_subject_length = (Int4) (total_length/num_seqs);
    //    }
    //
    //    if(min_subject_length <=0) {
    // 	   return BLASTERR_SUBJECT_LENGTH_INVALID;
    //    }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_message.c:215-219
    // ```c
    //     case BLASTERR_SUBJECT_LENGTH_INVALID:
    //         new_msg->message = strdup("The average subject length is too short");
    //         new_msg->severity = eBlastSevFatal;
    //         new_msg->context = context;
    //         break;
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:310-314
    // ```c
    //             retval = CPrelimSearchRunner(*m_InternalData, opts_memento.get())();
    //             if (retval) {
    //                 NCBI_THROW(CBlastException, eCoreBlastError,
    //                            BlastErrorCode2String(retval));
    //             }
    // ```
    // The set-up of the search comes after `CheckInternalData` (above): a batch with a valid
    // context stops NCBI when the subjects have fewer letters than subjects (subjects
    // without letters count with none; `AUTHORITY.md` §E3, BI-44). TBLASTX is ungapped, so
    // its Gumbel block is freed and the average length decides. The error comes after the
    // reports of the batches before (`ErrorAfterBatchReports`, `run_in_pool`).
    let (total_length, num_seqs) = if db_length_nucl > 0 {
        (db_length_nucl, db_num_seqs)
    } else {
        (
            if subjects.count() > 0 {
                subjects.length(0) as i64
            } else {
                0
            },
            1,
        )
    };
    if total_length / num_seqs <= 0 {
        return Err(ErrorAfterBatchReports(
            crate::cli::NativeError {
                exit: 3,
                message: "BLAST engine error: The average subject length is too short\n"
                    .to_string(),
            }
            .into(),
        )
        .into());
    }
    // Precompute per-context cutoff scores using NCBI BLAST algorithm.
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:280-419
    //
    // NCBI cutoff calculation for tblastx ungapped path:
    // 1. gap_trigger from ungapped params (kbp_std)
    // 2. cutoff_score_max from BlastHitSavingParametersNew (uses user's E-value)
    // 3. Initial word cutoff: cutoff_score_for_update_tblastx with CUTOFF_E_TBLASTX=1e-300
    // 4. Final cutoff = MIN(update_cutoff, gap_trigger, cutoff_score_max)
    //
    // The word cutoff uses the average subject nucleotide length. With E = 1e-300 (step 3) its
    // score exceeds gap_trigger whatever the lengths (above 2,000 against about 40), so step 4
    // never takes it and no output depends on this average; it is computed as NCBI computes it.
    let subject_len_nucl = db_length_nucl / db_num_seqs;
    // NCBI: cutoff scores are stored per query context (no subject-frame dimension).
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:320-324
    // ```c
    // BlastUngappedCutoffs *curr_cutoffs = parameters->cutoffs + context;
    // ...
    // curr_cutoffs->cutoff_score = new_cutoff;
    // ```
    let mut cutoff_scores: Vec<i32> = vec![0; contexts_ref.len()];
    // NCBI word_params->cutoff_score_min = min of cutoffs across all contexts
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:401-403
    let mut cutoff_score_min = i32::MAX;

    // =======================================================================
    // NCBI Parity: Pre-compute length_adjustment and eff_searchsp per context
    // =======================================================================
    // NCBI stores these in query_info->contexts[ctx].length_adjustment and
    // query_info->contexts[ctx].eff_searchsp via BLAST_CalcEffLengths
    // (ncbi-blast/c++/src/algo/blast/core/blast_setup.c:700-850).
    // Compute once for the complete subject set and share with all subject jobs.
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:846-847
    //   query_info->contexts[index].eff_searchsp = effective_search_space;
    //   query_info->contexts[index].length_adjustment = length_adjustment;
    let mut length_adj_per_context: Vec<i64> = Vec::with_capacity(contexts_ref.len());
    let mut eff_searchsp_per_context: Vec<i64> = Vec::with_capacity(contexts_ref.len());

    for (ctx_idx, ctx) in contexts_ref.iter().enumerate() {
        // NCBI reference: c++/src/algo/blast/core/blast_parameters.c:324-331;
        // blast_setup.c:774-847: invalid contexts keep zero effective lengths.
        // if (!query_info->contexts[context].is_valid) {
        //     curr_cutoffs->cutoff_score = INT4_MAX; continue;
        // }
        if !ctx.is_valid || ctx.aa_len == 0 {
            cutoff_scores[ctx_idx] = i32::MAX;
            length_adj_per_context.push(0);
            eff_searchsp_per_context.push(0);
            continue;
        }
        // NCBI: per-context kbp_std[context] is used throughout cutoff/length calcs.
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2797
        let ctx_params = &ctx.karlin_params;

        // NCBI: query_length = query_info->contexts[context].query_length
        let query_len_aa = ctx.aa_len as i64;

        // NCBI: gap_trigger uses kbp_std[context]->Lambda/logK.
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:340-345
        let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE, ctx_params);

        // =======================================================================
        // NCBI Parity: Use compute_eff_lengths_tblastx to get BOTH
        // length_adjustment and eff_searchsp from a single source of truth.
        // This mirrors BLAST_CalcEffLengths which computes and stores both values.
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:821-847
        // =======================================================================
        // `-dbsize` replaces the database length here only (`SubjectSetNumbers`).
        let eff_db_length_nucl = numbers.effective_length();
        let eff_lengths = crate::algorithm::tblastx::ncbi_cutoffs::try_compute_eff_lengths_tblastx(
            query_len_aa,
            eff_db_length_nucl,
            db_num_seqs,
            ctx_params, // tblastx uses per-context ungapped params (kbp_gap is NULL)
        )
        .ok_or(SearchSpaceOverflow {
            db_length: eff_db_length_nucl,
            query_length: ctx.aa_len,
        })?;
        let eff_searchsp = eff_lengths.eff_searchsp;
        length_adj_per_context.push(eff_lengths.length_adjustment);
        eff_searchsp_per_context.push(eff_searchsp);

        // Step 1: Compute cutoff_score_max from BlastHitSavingParametersNew
        // This uses the effective search space WITH length adjustment
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:942-946
        let cutoff_score_max = cutoff_score_max_for_tblastx(
            eff_searchsp,
            evalue_threshold, // User's E-value (typically 10.0)
            ctx_params,
        );

        // Step 2: Compute per-subject cutoff using BlastInitialWordParametersUpdate
        // This uses CUTOFF_E_TBLASTX=1e-300 and a simple searchsp formula
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:348-374
        let cutoff = cutoff_score_for_update_tblastx(
            query_len_aa,
            subject_len_nucl, // NUCLEOTIDE length, NOT divided by 3!
            gap_trigger,
            cutoff_score_max,
            BLAST_GAP_DECAY_RATE, // 0.5
            ctx_params,
            1.0, // scale_factor (standard BLOSUM62)
        );

        // DEBUG: Print cutoff values for first context
        static PRINTED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
        if diag_enabled && ctx_idx == 0 && !PRINTED.swap(true, std::sync::atomic::Ordering::Relaxed)
        {
            eprintln!(
                "[DEBUG CUTOFF] query_len_aa={}, subject_len_nucl={}",
                query_len_aa, subject_len_nucl
            );
            eprintln!("[DEBUG CUTOFF] eff_searchsp={}", eff_searchsp);
            eprintln!(
                "[DEBUG CUTOFF] length_adjustment={}",
                eff_lengths.length_adjustment
            );
            eprintln!("[DEBUG CUTOFF] cutoff_score_max={}", cutoff_score_max);
            eprintln!("[DEBUG CUTOFF] gap_trigger={}", gap_trigger);
            eprintln!("[DEBUG CUTOFF] final cutoff={}", cutoff);
        }
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:943-946
        // ```c
        // BLAST_Cutoffs(&new_cutoff, &evalue, kbp, searchsp, FALSE, 0);
        // params->cutoffs[context].cutoff_score = new_cutoff;
        // params->cutoffs[context].cutoff_score_max = new_cutoff;
        // ```
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:360-374
        // ```c
        // BLAST_Cutoffs(&new_cutoff, &cutoff_e, kbp,
        //               MIN((Uint8)subj_length, (Uint8)query_length)*((Uint8)subj_length),
        //               TRUE, gap_decay_rate);
        // new_cutoff = MIN(new_cutoff, gap_trigger);
        // new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
        // ```
        // Diagnostic-only dump of the same per-context values that feed
        // word_params->cutoff_score_min and CalculateLinkHSPCutoffs.
        if debug_cutoffs_all {
            eprintln!(
                "[DEBUG CUTOFF_ALL] ctx_idx={} q_frame={} query_len_aa={} subject_len_nucl={} eff_searchsp={} length_adjustment={} lambda={:.12e} k={:.12e} h={:.12e} cutoff_score_max={} gap_trigger={} word_cutoff={}",
                ctx_idx,
                ctx.frame,
                query_len_aa,
                subject_len_nucl,
                eff_searchsp,
                eff_lengths.length_adjustment,
                ctx_params.lambda,
                ctx_params.k,
                ctx_params.h,
                cutoff_score_max,
                gap_trigger,
                cutoff
            );
        }

        // Track minimum cutoff for linking
        cutoff_score_min = cutoff_score_min.min(cutoff);

        // All subject frames use the same cutoff (NCBI: per-context cutoffs only).
        cutoff_scores[ctx_idx] = cutoff;
    }
    // NCBI reference: c++/src/algo/blast/api/blast_results.cpp:82-100
    // ```c++
    //     // find the first valid context corresponding to this query
    //     for (i = 0; i < context_per_query; i++) {
    //         BlastContextInfo *ctx = query_info->contexts +
    //                                 query_number * context_per_query + i;
    //         if (ctx->is_valid) {
    //             m_SearchSpace = ctx->eff_searchsp;
    // 	    m_LengthAdjustment = ctx->length_adjustment;
    //             break;
    //     ...
    //     const int ctx_index = query_number * context_per_query + i;
    //     if (sbp->kbp_std) {
    //         s_InitializeKarlinBlk(sbp->kbp_std[ctx_index], &m_UngappedKarlinBlk);
    // ```
    // The per-query footer of the reports shows the first valid context's block and
    // search space. The SEG intervals of each frame are the query masks of the report.
    let mut query_stats: Vec<TblastxQueryStats> = query_frames
        .iter()
        .map(|frames| TblastxQueryStats {
            karlin: None,
            eff_searchsp: 0,
            aa_len: 0,
            context0_search_space: 0,
            seg_masks: frames
                .iter()
                .map(|frame| (frame.frame, frame.seg_masks.clone()))
                .collect(),
        })
        .collect();
    for (ctx_idx, ctx) in contexts_ref.iter().enumerate() {
        let stats = &mut query_stats[ctx.q_idx as usize];
        // Context 0 is the query's frame +1 (invalid contexts hold search space 0).
        if ctx.frame == 1 {
            stats.context0_search_space = eff_searchsp_per_context[ctx_idx];
        }
        if ctx.is_valid && stats.karlin.is_none() {
            stats.karlin = Some(ctx.karlin_params);
            stats.eff_searchsp = eff_searchsp_per_context[ctx_idx];
            stats.aa_len = ctx.aa_len;
        }
    }
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:401-416
    // ```c
    // if (new_cutoff < cutoff_min) {
    //    cutoff_min = new_cutoff;
    // }
    // parameters->cutoff_score_min = cutoff_min;
    // ```
    if debug_cutoffs_all {
        eprintln!(
            "[DEBUG CUTOFF_ALL] word_params_cutoff_score_min={}",
            cutoff_score_min
        );
    }

    // NCBI reference: c++/src/algo/blast/core/blast_parameters.c:101-110
    // if (kbp[index] && query_info->contexts[index].is_valid &&
    //     kbp[index]->Lambda > 0.0 && kbp[index]->Lambda < min_lambda) { ... }
    let context_params: Vec<KarlinParams> = contexts_ref
        .iter()
        .filter(|ctx| ctx.is_valid)
        .map(|ctx| ctx.karlin_params)
        .collect();
    let linking_params_for_cutoff = find_smallest_lambda_params(&context_params)
        .context("no valid query statistical parameters")?;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:82-88
    // ```c
    // if (num_threads > 1) {
    //     SetNumberOfThreads(num_threads);
    // }
    // ```
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let writer = if use_channel {
        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
        // ```c
        // typedef struct BlastHSPList {
        //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
        //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
        //                       Set to 0 if not applicable */
        // } BlastHSPList;
        // ```
        // The writer thread is NCBI's collector (tblastx/hit_collector.rs): it takes each
        // subject's message and returns the batch's hit lists in report order after the join,
        // for the one place that writes every output.
        let rx = rx_opt.take().expect("rx must be available for writer");
        let mut collector = HitCollector::new(args, query_records);
        Some(std::thread::spawn(move || -> Vec<TblastxHsp> {
            while let Ok((s_idx, hits)) = rx.recv() {
                collector.push(s_idx, hits);
            }
            collector.finish()
        }))
    } else {
        None
    };

    // One message per subject of the loop, with its final HSPs or none, so that the collector
    // reaches every OID (`HitCollector::push`): `st.hits` holds the subject's HSPs until here.
    let emit_subject_hits = |st: &mut WorkerState<S>, s_idx: usize| {
        let hits = std::mem::take(&mut st.hits);
        if let Some(tx) = &st.tx {
            tx.send((s_idx, hits)).unwrap();
        } else {
            local_collector
                .lock()
                .unwrap()
                .as_mut()
                .expect("a subject loop without a writer thread keeps the collector")
                .push(s_idx, hits);
        }
    };

    let process_subject = |st: &mut WorkerState<S>, s_idx: usize| {
        // No NCBI counterpart: the clock of this subject for LOSAT_TIMING. The first local of
        // the closure, so it is dropped last and its last stage includes the drops of the rest.
        let mut x_lap = x_times.subject();
        // The scores of the search's matrix, an argument of the closures below.
        let sc = st.sc;
        // EXPERIMENT (LOSAT_X_TBXPIPE): this subject's part of the batch's seed-stage end, ended at
        // dispatch point 5 below with the subject's HSPs, or with none when the closure returns
        // before it (tblastx/x_batch_pipe.rs). Nothing without the pipe.
        let mut x_pipe_seed = x_batch_pipe::XSeedGuard::new(x_seed_end);
        // Dispatch point 2 of the seed release (per subject): the seed stage of
        // this subject holds the batch's lookup table through `x_seed_stage` and the unit path
        // through `x_units_stage`. `lookup_ref` borrows from the first and `x_pair_par_ref`
        // from the second, so neither can be used after its stage has ended (`x_units_stage`
        // at dispatch point 3, `x_seed_stage` at dispatch point 5). A stage that is dropped
        // (the return just below) ends there. The run stops here if tables were already
        // released.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:503-505
        // ```c
        //     /* scan the subject sequence for hits */
        //     hits = scansub(lookup_wrap, subject,
        //                               offset_pairs, array_size, scan_range);
        // ```
        // `lookup_ref` is NCBI's `lookup_wrap`, the one table every subject is scanned with.
        let x_seed_stage = x_seed_batch.lookup.enter(s_idx);
        let x_units_stage = x_seed_batch.units.enter(s_idx);
        let lookup_ref: &crate::algorithm::tblastx::lookup::BlastAaLookupTable =
            x_seed_stage.tables();
        let x_pair_par_ref: Option<&x_pair_par::XPairPar> = if !S::BLOSUM62 {
            None
        } else {
            x_units_stage.tables().as_ref()
        };
        // NCBI reference: c++/src/algo/blast/core/blast_engine.c:1318-1330,1429-1431
        // return word_length * 3 + 2;
        // if (subject->length < min_subj_seq_length) { ... continue; }
        // Statistics include all records, including ones too short to search.
        let s_len = subjects.length(s_idx);
        if s_len < (wordsize * 3 + 2) as usize {
            emit_subject_hits(st, s_idx);
            return;
        }
        // NCBI: Creates ewp (diagonal table) ONCE per SUBJECT via BlastExtendWordNew
        // (blast_engine.c:1002). Each subject sequence gets a fresh diagonal array.
        // Reference: blast_extend.c:109-180 (BlastExtendWordNew) allocates with calloc,
        // which zeros all entries. The diag_offset is initialized to window.
        //
        // Reset diagonal state for each subject to match NCBI behavior:
        for d in st.diag_array.iter_mut() {
            *d = DiagStruct::default();
        }
        st.diag_offset = window;
        x_lap.lap(XStage::Enter);

        // NCBI reference (translate all frames for subject sequences):
        // ncbi-blast/c++/src/algo/blast/core/blast_util.c:1296-1308
        // ```c
        // for (context = 0; context < num_frames; ++context) {
        //    int frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
        //    retval->translations[context] = (Uint1*) malloc(...);
        //    BLAST_GetTranslation(subject_blk->sequence_start, nucl_seq_rev,
        //       subject_blk->length, frame, retval->translations[context], gen_code_string);
        // }
        // ```
        // Only the frames the preliminary stage translates are made here, in windows
        // (`Subjects::x_preliminary_frames`). For a subject with ambiguous bases the frames of
        // its letters (`s_frames`) are made after the seed stage, once the preliminary frames
        // are freed, and only when HSPs are left for the re-evaluation (`x_letter_frames_later`,
        // below); until then `s_frames` is empty and the HSP lists index the preliminary frames,
        // which have the same number and lengths (`s_frames_listed`). For a subject without
        // ambiguous bases both are the same frames.
        let mut x_letter_frames_later = false;
        let (s_frames, s_frames_preliminary_owned) = {
            let (frames, of_letters) = subjects.x_preliminary_frames(s_idx, &db_code);
            if of_letters {
                (frames, None)
            } else {
                x_letter_frames_later = true;
                (Vec::new(), Some(frames))
            }
        };
        let s_frames_report = &s_frames;
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:772-775
        // ```c
        //             BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
        //                                      backup.full_range.right,
        //                                      subject->gen_code_string, &translation_buffer,
        //                                      &frame_offsets, NULL);
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1492-1497
        // ```c
        //                status =
        //                   Blast_HSPListReevaluateUngapped(
        //                             program_number, hsp_list, query,
        //                             seq_arg.seq, word_params, hit_params,
        //                             query_info, sbp, score_params, seq_src,
        //                             seq_arg.seq->gen_code_string);
        // ```
        // The preliminary search scans and extends the translation of the compressed
        // subject, whose ambiguous bases NCBI replaced by compatible bases drawn from
        // `CRandom` (`preliminary_subject_bases`); the re-evaluation and the identities
        // read the frames of the ncbi4na subject (`s_frames`, ambiguous codons are X).
        let s_frames_preliminary: &[QueryFrame] =
            s_frames_preliminary_owned.as_deref().unwrap_or(&s_frames);
        // The frames the HSP lists of the seed stage index (`get_ungapped_hsp_list`, the frame
        // bases of the first linking call): only their number and lengths are read, which the
        // frames of the letters and the preliminary frames share.
        let s_frames_listed: &[QueryFrame] = if x_letter_frames_later {
            s_frames_preliminary
        } else {
            &s_frames
        };
        x_lap.lap(XStage::Translate);

        // [C] BlastOffsetPair *offset_pairs
        let offset_pairs = &mut st.offset_pairs;

        // [C] DiagStruct *diag_array = diag->hit_level_array;
        let diag_array = &mut st.diag_array;

        // [C] diag_offset = diag->offset;  (reset to window per-subject)
        let mut diag_offset: i32 = st.diag_offset;

        // NCBI reference: c++/src/algo/blast/core/blast_engine.c:1448-1455
        // CalculateLinkHSPCutoffs(program_number, query_info, gap_align->sbp,
        //     hit_params->link_hsp_params, word_params, db_length, subject->length);
        // Compute once before scanning; both linking passes read the same state.
        let subject_len_nucl = s_len as i64;
        let linking_params = LinkingParams {
            subject_len_nucl,
            gap_decay_rate: BLAST_GAP_DECAY_RATE,
            cutoffs: calculate_link_hsp_cutoffs_ncbi(
                avg_query_length,
                subject_len_nucl,
                db_length_nucl,
                cutoff_score_min,
                1.0,
                BLAST_GAP_DECAY_RATE,
                &linking_params_for_cutoff,
            ),
        };

        // NCBI: Each subject frame gets its own init_hitlist that is reset between frames.
        // Reference: blast_engine.c:491 BlastInitHitListReset(init_hitlist)
        // After each frame, init_hsps are converted and merged into combined_ungapped_hits.

        // Statistics for HSP saving analysis (long sequences only)
        let is_long_sequence = subject_len_nucl > 600_000;
        let collect_hsp_saving_stats = is_long_sequence && (debug_hsp_saving || diag_enabled);
        let mut stats_hsp_saved = 0usize;
        let mut stats_hsp_filtered_by_cutoff = 0usize;
        let mut stats_hsp_filtered_by_reeval = 0usize;
        let stats_hsp_filtered_by_hsp_test = 0usize;
        let mut stats_score_distribution: Vec<i32> = Vec::new();

        // NCBI: Combined HSP list across all subject frames
        // Reference: blast_engine.c:438 BlastHSPList* combined_hsp_list
        let mut combined_ungapped_hits: Vec<UngappedHit> = Vec::new();

        // EXPERIMENT (PAIRPAR), dispatch point 2 (per subject): the unit path fills
        // `combined_ungapped_hits` and the frame slice below becomes empty, so the reference
        // frame loop (unchanged) runs zero times.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805,835-844
        // ```c
        // for (context=first_context; context<=last_context; context++) {
        // ...
        //     status = s_BlastSearchEngineOneContext(program_number, query, query_info,
        // ...
        //     if (Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax)) {
        // ```
        // The unit path returns the list this frame loop builds (same calls, same order of
        // the appended frame lists); see x_pair_par.rs for why it is the same list.
        //
        // Dispatch point 3 of the seed release (per subject): the unit path runs
        // in two calls, the units (`run_subject_units`, the only reader of the group tables)
        // and the merge of their lists (`merge_subject`, a function without `self`), and the
        // units stage of this subject ends between them. What is still needed of the unit path
        // is read before that: `x_pair_par_stub` for the assertion below. After
        // `x_units_stage.end()` neither `x` nor `x_pair_par_ref` can be used. The batch's last
        // units stage frees the unit path there, before the merge, where the seed stage's
        // memory peaks.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:561-562
        // ```c
        // BLAST_GetUngappedHSPList(init_hitlist, query_info, subject,
        //         hit_params->options, &hsp_list);
        // ```
        // NCBI converts and merges its init HSPs with the lookup table still allocated (freed
        // with the search's internal data); that conversion reads no lookup table.
        let x_pair_par_stub = x_pair_par_ref.is_some_and(|x| x.stub());
        let s_frames_preliminary: &[QueryFrame] = match x_pair_par_ref {
            Some(x) => {
                let x_subject = x_pair_par::XSubject {
                    s_frames_preliminary,
                    s_frames: s_frames_listed,
                    contexts: contexts_ref,
                    cutoff_scores: &cutoff_scores,
                    x_dropoff_per_context: &x_dropoff_per_context,
                    window,
                    wordsize,
                    s_idx,
                    s_len,
                };
                x_lap.lap(XStage::Prepare);
                let x_units = x.run_subject_units(&x_subject, parallel_pool);
                x_lap.lap(XStage::Units);
                x_units_stage.end();
                x_lap.lap(XStage::UnitsRelease);
                let x_hits = x_pair_par::XPairPar::merge_subject(x_units, &x_subject);
                x_lap.lap(XStage::Merge);
                combined_ungapped_hits = x_hits;
                &[]
            }
            None => {
                x_lap.lap(XStage::Prepare);
                x_units_stage.end();
                x_lap.lap(XStage::UnitsRelease);
                s_frames_preliminary
            }
        };
        // EXPERIMENT (PAIRLUT): with the direct group tables the batch's `lookup` is a
        // stub without entries, so a scan of it by the frame loop below would find no hit. The
        // unit path above must have taken the subject (empty frame slice); anything else stops
        // the run here.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:503-505
        // ```c
        //     /* scan the subject sequence for hits */
        //     hits = scansub(lookup_wrap, subject,
        //                               offset_pairs, array_size, scan_range);
        // ```
        if x_pair_lut_is_stub {
            assert!(
                s_frames_preliminary.is_empty() && x_pair_par_stub,
                "PAIRLUT: the reference frame loop would scan the stub lookup table (subject {s_idx})"
            );
        }
        // Dispatch point 4 of the seed release (per subject): the frame loop below
        // indexes the worker's diagonal table through a raw pointer with `diag_mask`. A worker
        // of a batch that the unit path takes has a table of zero cells, and the worker that ended
        // the batch's last lookup stage has dropped its table. The loop must then run zero times;
        // anything else stops the run here.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-516
        // ```c
        // for (i = 0; i < hits; ++i) {
        //     Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
        //     Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
        // ...
        //     diag_coord = (query_offset - subject_offset) & diag_mask;
        // ```
        assert!(
            s_frames_preliminary.is_empty() || diag_array.len() == diag_array_size as usize,
            "the reference frame loop would run without the worker's diagonal table (subject {s_idx})"
        );

        for (s_f_idx, s_frame) in s_frames_preliminary.iter().enumerate() {
            // NCBI: Diagonal state is NOT reset between subject frames.
            // NCBI shares ewp (diag_table) across all 6 subject frame iterations.
            // Reference: blast_engine.c:805-855
            // The per-subject reset (above, line ~1646) matches NCBI's Blast_ExtendWordNew calloc.
            // Blast_ExtendWordExit at the end of each subject chunk increments diag_offset.

            let subject_full = &s_frame.aa_seq;
            // NCBI subject->sequence points past the leading NULLB sentinel, so offsets
            // are 0-based from the first residue; see blast_engine.c:811-812 and
            // blast_util.c:112-116.
            let subject_all = &subject_full[1..subject_full.len() - 1];
            let s_aa_len = s_frame.aa_len;

            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:804-855
            // ```c
            // for (context=first_context; context<=last_context; context++) {
            //    subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
            //    subject->sequence = translation_buffer + frame_offsets[context] + 1;
            //    subject->length = frame_offsets[context+1] - frame_offsets[context] - 1;
            //    status = s_BlastSearchEngineOneContext(..., &hsp_list_for_chunks, ...);
            //    Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax);
            // }
            // ```
            let mut frame_ungapped_hits: Vec<UngappedHit> = Vec::new();
            let mut split_state = SubjectSplitState::new(s_aa_len);
            let max_dbseq_len = tblastx_max_dbseq_len_for_run();

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-516
            // ```c
            // while (scan_range[1] <= scan_range[2]) {
            // ...
            //     for (i = 0; i < hits; ++i) {
            //         Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
            //         Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
            // ...
            //         diag_coord = (query_offset - subject_offset) & diag_mask;
            // ```
            // This closure is the reference loop (only renamed from scan_chunk): the port of the C loop above.
            // The dispatcher `scan_chunk` below chooses between it and `scan_chunk_bucketed`.
            let scan_chunk_plain = |chunk: SubjectChunk,
                                    offset_pairs: &mut [OffsetPair],
                                    diag_array: &mut [DiagStruct],
                                    diag_offset: &mut i32,
                                    scan_chunk_size: Option<usize>,
                                    sc: S|
             -> TblastxChunkScanResult {
                let mut init_hsps: Vec<InitHSP> = Vec::new();
                let mut stats = TblastxChunkScanStats::default();
                let chunk_end = chunk.offset.saturating_add(chunk.length);
                let subject = &subject_all[chunk.offset..chunk_end];
                // EXPERIMENT (LOSAT_X_TBXUNGAPSENT): the chunk's extensions take the walks the
                // sentinels stop where it is the whole frame; never on wasm32.
                let subject_sent = cfg!(not(target_arch = "wasm32"))
                    && x_tbx_sent_subject(tbx_sent, subject_full, chunk.offset, chunk.length);

                if subject.len() < wordsize as usize {
                    // NCBI still advances the diagonal table offset even when no hits can be found.
                    // References: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:445-446,
                    // ncbi-blast/c++/src/algo/blast/core/blast_extend.c:167-173
                    advance_tblastx_diag_offset(diag_offset, diag_array, window, chunk.length);
                    return TblastxChunkScanResult {
                        chunk,
                        hits: Vec::new(),
                        stats,
                    };
                }

                // NCBI subject seq_ranges are used by s_DetermineScanningOffsets (masksubj.inl).
                // With no subject masking, the range is [0, subject->length].
                // References: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:509-511,
                // ncbi-blast/c++/src/algo/blast/core/masksubj.inl:43-58
                let base_seq_ranges: [(i32, i32); 1] = [(0, chunk.length as i32)];
                let scan_interiors = tblastx_scan_interiors(chunk.length, scan_chunk_size);
                let mut previous_seed_s_off: Option<u32> = None;
                for (_scan_chunk_index, (interior_start, interior_end)) in
                    scan_interiors.into_iter().enumerate()
                {
                    let seq_ranges = clip_tblastx_seq_ranges_for_scan_interior(
                        &base_seq_ranges,
                        interior_start,
                        interior_end,
                        wordsize as usize,
                        subject.len(),
                    );
                    if seq_ranges.is_empty() {
                        continue;
                    }
                    // [C] scan_range[0] = 0;
                    // [C] scan_range[1] = subject->seq_ranges[0].left;
                    // [C] scan_range[2] = subject->seq_ranges[0].right - wordsize;
                    let mut scan_range: [i32; 3] = [0, seq_ranges[0].0, seq_ranges[0].1 - wordsize];

                    // [C] while (scan_range[1] <= scan_range[2])
                    while scan_range[1] <= scan_range[2] {
                        let prev_scan_left = scan_range[1];
                        // [C] hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
                        let t0 = if timing_enabled {
                            Some(Instant::now())
                        } else {
                            None
                        };
                        // EXPERIMENT (LOSAT_X_TBXSCANBF): the branch-free scan with the switch.
                        let hits = crate::algorithm::tblastx::x_scan_bf::x_scan_subject(
                            x_scan_bf,
                            lookup_ref,
                            subject,
                            &seq_ranges,
                            offset_pairs,
                            offset_array_size,
                            &mut scan_range,
                        );
                        if let Some(t0) = t0 {
                            scan_ns
                                .fetch_add(t0.elapsed().as_nanos() as u64, AtomicOrdering::Relaxed);
                            scan_calls.fetch_add(1, AtomicOrdering::Relaxed);
                        }

                        if diag_enabled && hits > 0 {
                            diagnostics
                                .base
                                .kmer_matches
                                .fetch_add(hits as usize, AtomicOrdering::Relaxed);

                            // DEBUG: Check for duplicate offset pairs in scan output
                            if is_long_sequence {
                                // NCBI BlastOffsetPair uses Uint4 offsets.
                                // Reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                                let mut seen: HashSet<(u32, u32)> =
                                    HashSet::with_capacity(hits as usize);
                                let mut duplicate_count = 0usize;
                                for i in 0..hits as usize {
                                    let pair = unsafe { &*offset_pairs.as_ptr().add(i) };
                                    if !seen.insert((pair.q_off, pair.s_off)) {
                                        duplicate_count += 1;
                                    }
                                }
                                if duplicate_count > 0 {
                                    eprintln!("[DEBUG SCAN_DUPES] s_f_idx={} scan_range=[{},{}] hits={} duplicates={} ({:.2}%)",
                                        s_f_idx, prev_scan_left, scan_range[1], hits, duplicate_count,
                                        (duplicate_count as f64 / hits as f64) * 100.0);
                                }
                            }
                        }

                        if hits == 0 && scan_range[1] == prev_scan_left {
                            // Safety guard: with correct NCBI-sized offset arrays, this should not happen.
                            // If it does, breaking avoids an infinite loop.
                            break;
                        }

                        // [C] for (i = 0; i < hits; ++i)
                        // OPTIMIZATION: Use raw pointers to eliminate bounds checking in hot loop
                        let diag_ptr = diag_array.as_mut_ptr();
                        let offset_pairs_ptr = offset_pairs.as_ptr();

                        for i in 0..hits as usize {
                            // SAFETY: i < hits, and hits <= offset_array_size (checked by scan)
                            let pair = unsafe { &*offset_pairs_ptr.add(i) };
                            let query_offset = pair.q_off;
                            let subject_offset = pair.s_off;
                            debug_assert!(
                            previous_seed_s_off
                                .map(|prev| subject_offset >= prev)
                                .unwrap_or(true),
                            "TBLASTX scan chunks must replay seeds in nondecreasing subject offset order"
                        );
                            previous_seed_s_off = Some(subject_offset);

                            // Debug: dump scan output around a target subject offset.
                            // The scan output is produced by s_BlastAaScanSubject.
                            // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:83-123
                            if let Some((lo, hi)) = scan_debug_range {
                                // NCBI offsets are Uint4; cast for debug range comparison only.
                                // Reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                                let subject_offset_i32 =
                                    subject_offset as i32 + chunk.offset as i32;
                                if subject_offset_i32 >= lo && subject_offset_i32 <= hi {
                                    eprintln!(
                                        "[DEBUG SCAN_OFF] s_f_idx={} s_off={} local_s_off={} q_off={} scan_range=[{},{}] chunk_offset={} diag_offset={}",
                                        s_f_idx,
                                        subject_offset_i32,
                                        subject_offset,
                                        query_offset,
                                        prev_scan_left,
                                        scan_range[1],
                                        chunk.offset,
                                        diag_offset
                                    );
                                }
                            }

                            // [C] diag_coord = (query_offset - subject_offset) & diag_mask;
                            // NCBI uses Uint4 for offsets; apply unsigned wrapping.
                            // References: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                            //             ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:534
                            let diag_coord = (query_offset.wrapping_sub(subject_offset)
                                & (diag_mask as u32))
                                as usize;

                            // SAFETY: diag_coord is masked by diag_mask, which is < diag_array.len()
                            let diag_entry = unsafe { &mut *diag_ptr.add(diag_coord) };

                            // [C] if (diag_array[diag_coord].flag)
                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:536-553
                            if diag_entry.flag() != 0 {
                                // [C] if ((Int4)(subject_offset + diag_offset) < diag_array[diag_coord].last_hit)
                                let subject_plus_offset =
                                    subject_offset.wrapping_add(*diag_offset as u32);
                                if subject_plus_offset < diag_entry.last_hit() as u32 {
                                    if diag_enabled {
                                        diagnostics
                                            .base
                                            .seeds_masked
                                            .fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    continue;
                                }
                                // [C] diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                                // [C] diag_array[diag_coord].flag = 0;
                                diag_entry.set_last_hit(subject_plus_offset as i32);
                                diag_entry.set_flag(0);
                                // Track flag reset (hit after previous extension zone)
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .seeds_flag_reset
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                            }
                            // [C] else
                            else {
                                // [C] last_hit = diag_array[diag_coord].last_hit - diag_offset;
                                let last_hit = diag_entry.last_hit() - *diag_offset;
                                // [C] diff = subject_offset - last_hit;
                                // NCBI uses Uint4 for subject_offset; compute with unsigned wrap.
                                // References: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                                //             ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:559-560
                                let diff = subject_offset.wrapping_sub(last_hit as u32) as i32;

                                // [C] if (diff >= window)
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:562-569
                                if diff >= window {
                                    if diag_enabled {
                                        diagnostics
                                            .base
                                            .seeds_second_hit_too_far
                                            .fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    diag_entry.set_last_hit(
                                        subject_offset.wrapping_add(*diag_offset as u32) as i32,
                                    );
                                    continue;
                                }

                                // [C] if (diff < wordsize)
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:573-580
                                if diff < wordsize {
                                    if diag_enabled {
                                        diagnostics
                                            .base
                                            .seeds_second_hit_overlap
                                            .fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    continue;
                                }

                                // [C] curr_context = BSearchContextInfo(query_offset, query_info);
                                // NCBI passes Uint4 query_offset into BSearchContextInfo (Int4).
                                // References: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                                //             ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:590
                                let ctx_idx = lookup_ref.get_context_idx(query_offset as i32);
                                let ctx = unsafe { contexts_ref.get_unchecked(ctx_idx) };
                                let q_raw =
                                    query_offset.wrapping_sub(ctx.frame_base as u32) as usize;
                                // NCBI uses masked sequence for extension; query->sequence is
                                // sequence_start + 1, so offsets are 0-based in that buffer.
                                // Reference: blast_query_info.c:311-315, blast_util.c:112-116.
                                let query_full = &ctx.aa_seq;
                                let query = &query_full[1..query_full.len() - 1];

                                // [C] if (query_offset - diff < query_info->contexts[curr_context].query_offset)
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:592-606
                                let q_minus_diff = query_offset.wrapping_sub(diff as u32);
                                if q_minus_diff < ctx.frame_base as u32 {
                                    if diag_enabled {
                                        diagnostics
                                            .base
                                            .seeds_ctx_boundary_fail
                                            .fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    diag_entry.set_last_hit(
                                        subject_offset.wrapping_add(*diag_offset as u32) as i32,
                                    );
                                    continue;
                                }

                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .seeds_second_hit_window
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                    diagnostics
                                        .base
                                        .seeds_passed
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }

                                // [C] cutoffs = word_params->cutoffs + curr_context;
                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:686-688
                                // ```c
                                // Int4 cutoff_score = word_params->cutoffs[hsp->context].cutoff_score;
                                // ```
                                let cutoff = unsafe { *cutoff_scores.get_unchecked(ctx_idx) };
                                // [C] cutoffs->x_dropoff (per-context x_dropoff)
                                // Reference: aa_ungapped.c:579
                                let x_dropoff =
                                    unsafe { *x_dropoff_per_context.get_unchecked(ctx_idx) };

                                // [C] score = s_BlastAaExtendTwoHit(matrix, subject, query,
                                //                                   last_hit + wordsize, subject_offset, query_offset, ...)
                                // Two-hit ungapped extension (NCBI `s_BlastAaExtendTwoHit`)
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1089-1158
                                let t0 = if timing_enabled {
                                    Some(Instant::now())
                                } else {
                                    None
                                };
                                let (
                                    hsp_q_u,
                                    hsp_qe_u,
                                    hsp_s_u,
                                    _hsp_se_u,
                                    score,
                                    right_extend,
                                    s_last_off_u,
                                ) = if !subject_sent {
                                    extend_hit_two_hit_with(
                                        sc,
                                        query,
                                        subject,
                                        (last_hit + wordsize) as usize,
                                        subject_offset as usize,
                                        q_raw as usize,
                                        x_dropoff,
                                        // NCBI aa_ungapped.c:576-582 (call above):
                                        // score = s_BlastAaExtendTwoHit(...);
                                        // Search-local LOSAT diagnostic flag only.
                                        extension_debug_enabled,
                                    )
                                } else {
                                    // SAFETY: the contexts passed `x_tbx_sent_contexts` (tbx_sent,
                                    // BLOSUM62 only) and `x_tbx_sent_subject` accepted the chunk
                                    // (`subject_sent`).
                                    unsafe {
                                        extend_hit_two_hit_sent(
                                            query_full,
                                            subject_full,
                                            (last_hit + wordsize) as usize,
                                            subject_offset as usize,
                                            q_raw as usize,
                                            x_dropoff,
                                        )
                                    }
                                };
                                if let Some(t0) = t0 {
                                    ungapped_ns.fetch_add(
                                        t0.elapsed().as_nanos() as u64,
                                        AtomicOrdering::Relaxed,
                                    );
                                    ungapped_calls.fetch_add(1, AtomicOrdering::Relaxed);
                                }

                                let hsp_q: i32 = hsp_q_u as i32;
                                let hsp_s: i32 = hsp_s_u as i32;
                                let hsp_len: i32 = (hsp_qe_u - hsp_q_u) as i32;
                                let s_last_off: i32 = s_last_off_u as i32;

                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .ungapped_extensions
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                    if right_extend {
                                        diagnostics
                                            .base
                                            .ungapped_two_hit_extensions
                                            .fetch_add(1, AtomicOrdering::Relaxed);
                                    } else {
                                        diagnostics
                                            .base
                                            .ungapped_one_hit_extensions
                                            .fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    if hsp_len > 0 {
                                        diagnostics
                                            .base
                                            .extension_total_length
                                            .fetch_add(hsp_len as usize, AtomicOrdering::Relaxed);
                                        atomic_max_usize(
                                            &diagnostics.base.extension_max_length,
                                            hsp_len as usize,
                                        );
                                    }
                                }

                                // NCBI: Update diagonal state based on right_extend
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:636-648
                                // if (right_extend) {
                                //     diag_array[diag_coord].flag = 1;
                                //     diag_array[diag_coord].last_hit = s_last_off - (wordsize - 1) + diag_offset;
                                // } else {
                                //     diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                                // }
                                if right_extend {
                                    diag_entry.set_flag(1);
                                    diag_entry
                                        .set_last_hit(s_last_off - (wordsize - 1) + *diag_offset);
                                    if diag_enabled {
                                        diagnostics
                                            .base
                                            .mask_updates
                                            .fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                } else {
                                    diag_entry.set_last_hit(
                                        subject_offset.wrapping_add(*diag_offset as u32) as i32,
                                    );
                                }

                                // [C] if (score >= cutoffs->cutoff_score)
                                // NCBI reference: aa_ungapped.c:575-591 (Extension後のcutoffチェック)
                                if collect_hsp_saving_stats {
                                    if score >= cutoff {
                                        stats.score_distribution.push(score);
                                        stats.hsp_saved += 1;
                                    } else {
                                        stats.hsp_filtered_by_cutoff += 1;
                                    }
                                }
                                if score >= cutoff {
                                    if diag_enabled {
                                        diagnostics
                                            .ungapped_only_hits
                                            .fetch_add(1, AtomicOrdering::Relaxed);
                                    }

                                    // Extra debug for a traced HSP: print seed/extension inputs and cutoffs.
                                    if let Some(target) = trace_hsp_target() {
                                        // Compute outfmt coords for this candidate init-hsp (same logic as trace_init_hsp_if_match).
                                        // NCBI offsets are 0-based in query/subject->sequence buffers.
                                        // Reference: blast_gapalign.c:4756-4768, blast_aascan.c:110-113.
                                        let q_aa_start = hsp_q_u as usize;
                                        let q_aa_end = hsp_qe_u as usize;
                                        let s_aa_start = hsp_s_u as usize + chunk.offset;
                                        let s_aa_end = _hsp_se_u as usize + chunk.offset;
                                        let (q_start_dna, q_end_dna) = convert_coords(
                                            q_aa_start,
                                            q_aa_end,
                                            ctx.frame,
                                            ctx.orig_len,
                                        );
                                        let (s_start_dna, s_end_dna) = convert_coords(
                                            s_aa_start,
                                            s_aa_end,
                                            s_frame.frame,
                                            s_len,
                                        );
                                        if trace_match_target(
                                            target,
                                            q_start_dna,
                                            q_end_dna,
                                            s_start_dna,
                                            s_end_dna,
                                        ) {
                                            eprintln!(
                                                "[TRACE_HSP] seed/extend ctx_idx={} s_f_idx={} q_frame={} s_frame={} score={} cutoff={} x_dropoff={} last_hit={} subject_offset={} chunk_offset={} diff={} q_raw={} query_offset={} diag_coord={} right_extend={} s_last_off={}",
                                                ctx_idx,
                                                s_f_idx,
                                                ctx.frame,
                                                s_frame.frame,
                                                score,
                                                cutoff,
                                                x_dropoff,
                                                last_hit,
                                                subject_offset,
                                                chunk.offset,
                                                diff,
                                                q_raw,
                                                query_offset,
                                                diag_coord,
                                                right_extend,
                                                s_last_off,
                                            );
                                            // NCBI two-hit gating checks (diff/window/wordsize/context).
                                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:531-606
                                            // NCBI uses Uint4 offsets; apply unsigned wrap here as well.
                                            // Reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                                            let q_minus_diff =
                                                query_offset.wrapping_sub(diff as u32);
                                            eprintln!(
                                                "[TRACE_HSP] two_hit_pass diff={} window={} wordsize={} diff>=window={} diff<wordsize={} q_minus_diff={} ctx_frame_base={} q_minus_diff<base={} diag_offset={} diag_mask={} diag_array_size={}",
                                                diff,
                                                window,
                                                wordsize,
                                                diff >= window,
                                                diff < wordsize,
                                                q_minus_diff,
                                                ctx.frame_base,
                                                q_minus_diff < ctx.frame_base as u32,
                                                diag_offset,
                                                diag_mask,
                                                diag_array_size
                                            );
                                        }
                                    }
                                    // NCBI: BlastSaveInitHsp equivalent
                                    // Reference: blast_extend.c:360-375 BlastSaveInitHsp
                                    // Store HSP with absolute coordinates (before coordinate conversion)
                                    //
                                    // hsp_q is frame-relative coordinate in query->sequence (0-based),
                                    // frame_base is the context query_offset in the concatenated buffer.
                                    // NCBI: ungapped_data->q_start is absolute query offset.
                                    // Reference: blast_gapalign.c:4756-4768, blast_query_info.c:311-315.
                                    let hsp_q_absolute = ctx.frame_base + hsp_q;
                                    let hsp_qe_absolute = ctx.frame_base + (hsp_q + hsp_len);

                                    let init = InitHSP {
                                        q_start_absolute: hsp_q_absolute,
                                        q_end_absolute: hsp_qe_absolute,
                                        s_start: hsp_s,
                                        s_end: hsp_s + hsp_len,
                                        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:589-592
                                        // ```c
                                        // BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                                        //                  query_offset, subject_offset, hsp_len,
                                        //                  score);
                                        // ```
                                        q_seed_absolute: query_offset as i32,
                                        s_seed: subject_offset as i32,
                                        score,
                                        ctx_idx,
                                        s_f_idx,
                                        q_idx: ctx.q_idx,
                                        s_idx: s_idx as u32,
                                        q_frame: ctx.frame,
                                        s_frame: s_frame.frame,
                                        q_orig_len: ctx.orig_len,
                                        s_orig_len: s_len,
                                    };
                                    trace_init_hsp_if_match("init_hsp_saved", &init, contexts_ref);
                                    init_hsps.push(init);
                                } else if diag_enabled {
                                    diagnostics
                                        .ungapped_cutoff_failed
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                    atomic_min_i32(
                                        &diagnostics.ungapped_cutoff_failed_min_score,
                                        score,
                                    );
                                    atomic_max_i32(
                                        &diagnostics.ungapped_cutoff_failed_max_score,
                                        score,
                                    );
                                }
                            }
                        }
                    }
                }

                // [C] Blast_ExtendWordExit(ewp, subject->length);
                //
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:614
                // ```c
                // /* increment the offset in the diagonal array */
                // Blast_ExtendWordExit(ewp, subject->length);
                // ```
                advance_tblastx_diag_offset(diag_offset, diag_array, window, chunk.length);

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:234-235
                // ```c
                // Blast_InitHitListSortByScore(init_hitlist);
                // return status;
                // ```
                // BlastAaWordFinder sorts the chunk's init hit list by
                // score_compare_match (blast_extend.c:274-310) before
                // BLAST_GetUngappedHSPList sees it.
                sort_init_hsps_by_score_ncbi(&mut init_hsps);
                let mut hits = if init_hsps.is_empty() {
                    Vec::new()
                } else {
                    // NCBI: BLAST_GetUngappedHSPList equivalent - per chunk conversion,
                    // then Blast_HSPListAdjustOffsets and Blast_HSPListsMerge.
                    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:561-584
                    // ```c
                    // BLAST_GetUngappedHSPList(init_hitlist, query_info, subject,
                    //         hit_params->options, &hsp_list);
                    // Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
                    // status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list,
                    //      kHspNumMax, &(backup.offset), INT4_MIN, overlap, ...);
                    // ```
                    get_ungapped_hsp_list(init_hsps, contexts_ref, s_frames_listed)
                };
                adjust_tblastx_chunk_subject_offsets(&mut hits, chunk.offset);
                TblastxChunkScanResult { chunk, hits, stats }
            };

            // EXPERIMENT (SEEDBUCKET, graduated on 2026-10-22: no variable): `scan_chunk_plain`
            // above is the reference loop, untouched.  This variant appends the hits of every scan
            // call to per-diagonal-range buckets and runs the same per-hit body (the macro
            // below, a copy of the loop body of `scan_chunk_plain`) bucket by bucket,
            // then puts the saved HSPs back into scan order.  See `x_seed_bucket`.
            //
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-516
            // ```c
            // while (scan_range[1] <= scan_range[2]) {
            // ...
            //     for (i = 0; i < hits; ++i) {
            //         Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
            //         Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
            // ...
            //         diag_coord = (query_offset - subject_offset) & diag_mask;
            // ```
            // This variant ports the same C loop with one change of order: the hits of a scan call are
            // grouped by range of diag_coord and processed range by range (see x_seed_bucket.rs).
            // The order of the hits within one diagonal, and so every diagonal cell state, is as in C.
            let scan_chunk_bucketed = |chunk: SubjectChunk,
                                       offset_pairs: &mut [OffsetPair],
                                       diag_array: &mut [DiagStruct],
                                       diag_offset: &mut i32,
                                       scan_chunk_size: Option<usize>,
                                       sc: S|
             -> TblastxChunkScanResult {
                let mut init_hsps: Vec<InitHSP> = Vec::new();
                let mut stats = TblastxChunkScanStats::default();
                let chunk_end = chunk.offset.saturating_add(chunk.length);
                let subject = &subject_all[chunk.offset..chunk_end];
                // EXPERIMENT (LOSAT_X_TBXUNGAPSENT): the chunk's extensions take the walks the
                // sentinels stop where it is the whole frame; never on wasm32.
                let subject_sent = cfg!(not(target_arch = "wasm32"))
                    && x_tbx_sent_subject(tbx_sent, subject_full, chunk.offset, chunk.length);

                if subject.len() < wordsize as usize {
                    // NCBI still advances the diagonal table offset even when no hits can be found.
                    // References: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:445-446,
                    // ncbi-blast/c++/src/algo/blast/core/blast_extend.c:167-173
                    advance_tblastx_diag_offset(diag_offset, diag_array, window, chunk.length);
                    return TblastxChunkScanResult {
                        chunk,
                        hits: Vec::new(),
                        stats,
                    };
                }

                // NCBI subject seq_ranges are used by s_DetermineScanningOffsets (masksubj.inl).
                // With no subject masking, the range is [0, subject->length].
                // References: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:509-511,
                // ncbi-blast/c++/src/algo/blast/core/masksubj.inl:43-58
                let base_seq_ranges: [(i32, i32); 1] = [(0, chunk.length as i32)];
                let scan_interiors = tblastx_scan_interiors(chunk.length, scan_chunk_size);
                let mut buckets = crate::algorithm::tblastx::x_seed_bucket::SeedBuckets::new(
                    diag_array_size as u32,
                    diag_mask as u32,
                );
                let x_bucket_budget = crate::algorithm::tblastx::x_seed_bucket::budget();
                // Per-flush positions of the HSPs the current flush has saved.
                let mut seq_keys: Vec<u32> = Vec::new();
                let diag_offset_value: i32 = *diag_offset;
                // The NCBI two-hit loop body for one hit (aa_ungapped.c:531-606): a copy
                // of the loop body of `scan_chunk_plain`, as a macro expanded in place in
                // the bucket flush callbacks.  Local names resolve at this definition site
                // (macro hygiene), so the body reads `init_hsps`, `seq_keys`, `stats` and
                // `diag_array` of this closure directly.  `break 'hit` is the `continue`
                // of the NCBI loop; `seq` is the hit's position within the flush.
                //
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-606
                // ```c
                // diag_coord = (query_offset - subject_offset) & diag_mask;
                // ...
                // if (diag_array[diag_coord].flag) {
                // ...
                //     last_hit = diag_array[diag_coord].last_hit - diag_offset;
                //     diff = subject_offset - last_hit;
                // ...
                //     if (diff >= window) {
                // ...
                //     score = s_BlastAaExtendTwoHit(matrix, subject, query,
                // ...
                //     if (score >= cutoffs->cutoff_score)
                // ...
                //     if (right_extend) {
                //         diag_array[diag_coord].flag = 1;
                //         diag_array[diag_coord].last_hit =
                //             s_last_off - (wordsize - 1) + diag_offset;
                // ```
                // The macro body is the loop body of scan_chunk_plain, which ports the C body above statement
                // for statement. The body reads and writes only diag_array[diag_coord], the sequences and the
                // matrix, so applying it bucket by bucket gives each diagonal the same sequence of states.
                macro_rules! x_two_hit_body {
                    ($qo:expr, $so:expr, $sq:expr) => {{
                    let query_offset: u32 = $qo;
                    let subject_offset: u32 = $so;
                    let seq: u32 = $sq;
                    'hit: {
                    let diag_ptr = diag_array.as_mut_ptr();
                    let diag_offset = &diag_offset_value;

                            // [C] diag_coord = (query_offset - subject_offset) & diag_mask;
                            // NCBI uses Uint4 for offsets; apply unsigned wrapping.
                            // References: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                            //             ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:534
                            let diag_coord = (query_offset.wrapping_sub(subject_offset)
                                & (diag_mask as u32))
                                as usize;

                            // SAFETY: diag_coord is masked by diag_mask, which is < diag_array.len()
                            let diag_entry = unsafe { &mut *diag_ptr.add(diag_coord) };

                            // [C] if (diag_array[diag_coord].flag)
                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:536-553
                            if diag_entry.flag() != 0 {
                                // [C] if ((Int4)(subject_offset + diag_offset) < diag_array[diag_coord].last_hit)
                                let subject_plus_offset =
                                    subject_offset.wrapping_add(*diag_offset as u32);
                                if subject_plus_offset < diag_entry.last_hit() as u32 {
                                    if diag_enabled {
                                        diagnostics
                                            .base
                                            .seeds_masked
                                            .fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    break 'hit;
                                }
                                // [C] diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                                // [C] diag_array[diag_coord].flag = 0;
                                diag_entry.set_last_hit(subject_plus_offset as i32);
                                diag_entry.set_flag(0);
                                // Track flag reset (hit after previous extension zone)
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .seeds_flag_reset
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                            }
                            // [C] else
                            else {
                                // [C] last_hit = diag_array[diag_coord].last_hit - diag_offset;
                                let last_hit = diag_entry.last_hit() - *diag_offset;
                                // [C] diff = subject_offset - last_hit;
                                // NCBI uses Uint4 for subject_offset; compute with unsigned wrap.
                                // References: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                                //             ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:559-560
                                let diff = subject_offset.wrapping_sub(last_hit as u32) as i32;

                                // [C] if (diff >= window)
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:562-569
                                if diff >= window {
                                    if diag_enabled {
                                        diagnostics
                                            .base
                                            .seeds_second_hit_too_far
                                            .fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    diag_entry.set_last_hit(
                                        subject_offset.wrapping_add(*diag_offset as u32) as i32,
                                    );
                                    break 'hit;
                                }

                                // [C] if (diff < wordsize)
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:573-580
                                if diff < wordsize {
                                    if diag_enabled {
                                        diagnostics
                                            .base
                                            .seeds_second_hit_overlap
                                            .fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    break 'hit;
                                }

                                // [C] curr_context = BSearchContextInfo(query_offset, query_info);
                                // NCBI passes Uint4 query_offset into BSearchContextInfo (Int4).
                                // References: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                                //             ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:590
                                let ctx_idx = lookup_ref.get_context_idx(query_offset as i32);
                                let ctx = unsafe { contexts_ref.get_unchecked(ctx_idx) };
                                let q_raw =
                                    query_offset.wrapping_sub(ctx.frame_base as u32) as usize;
                                // NCBI uses masked sequence for extension; query->sequence is
                                // sequence_start + 1, so offsets are 0-based in that buffer.
                                // Reference: blast_query_info.c:311-315, blast_util.c:112-116.
                                let query_full = &ctx.aa_seq;
                                let query = &query_full[1..query_full.len() - 1];

                                // [C] if (query_offset - diff < query_info->contexts[curr_context].query_offset)
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:592-606
                                let q_minus_diff = query_offset.wrapping_sub(diff as u32);
                                if q_minus_diff < ctx.frame_base as u32 {
                                    if diag_enabled {
                                        diagnostics
                                            .base
                                            .seeds_ctx_boundary_fail
                                            .fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    diag_entry.set_last_hit(
                                        subject_offset.wrapping_add(*diag_offset as u32) as i32,
                                    );
                                    break 'hit;
                                }

                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .seeds_second_hit_window
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                    diagnostics
                                        .base
                                        .seeds_passed
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }

                                // [C] cutoffs = word_params->cutoffs + curr_context;
                                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:686-688
                                // ```c
                                // Int4 cutoff_score = word_params->cutoffs[hsp->context].cutoff_score;
                                // ```
                                let cutoff = unsafe { *cutoff_scores.get_unchecked(ctx_idx) };
                                // [C] cutoffs->x_dropoff (per-context x_dropoff)
                                // Reference: aa_ungapped.c:579
                                let x_dropoff =
                                    unsafe { *x_dropoff_per_context.get_unchecked(ctx_idx) };

                                // [C] score = s_BlastAaExtendTwoHit(matrix, subject, query,
                                //                                   last_hit + wordsize, subject_offset, query_offset, ...)
                                // Two-hit ungapped extension (NCBI `s_BlastAaExtendTwoHit`)
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1089-1158
                                let t0 = if timing_enabled {
                                    Some(Instant::now())
                                } else {
                                    None
                                };
                                let (
                                    hsp_q_u,
                                    hsp_qe_u,
                                    hsp_s_u,
                                    _hsp_se_u,
                                    score,
                                    right_extend,
                                    s_last_off_u,
                                ) = if !subject_sent {
                                    extend_hit_two_hit_with(
                                        sc,
                                        query,
                                        subject,
                                        (last_hit + wordsize) as usize,
                                        subject_offset as usize,
                                        q_raw as usize,
                                        x_dropoff,
                                        // NCBI aa_ungapped.c:576-582 (call above):
                                        // score = s_BlastAaExtendTwoHit(...);
                                        // Search-local LOSAT diagnostic flag only.
                                        extension_debug_enabled,
                                    )
                                } else {
                                    // SAFETY: the contexts passed `x_tbx_sent_contexts` (tbx_sent,
                                    // BLOSUM62 only) and `x_tbx_sent_subject` accepted the chunk
                                    // (`subject_sent`).
                                    unsafe {
                                        extend_hit_two_hit_sent(
                                            query_full,
                                            subject_full,
                                            (last_hit + wordsize) as usize,
                                            subject_offset as usize,
                                            q_raw as usize,
                                            x_dropoff,
                                        )
                                    }
                                };
                                if let Some(t0) = t0 {
                                    ungapped_ns.fetch_add(
                                        t0.elapsed().as_nanos() as u64,
                                        AtomicOrdering::Relaxed,
                                    );
                                    ungapped_calls.fetch_add(1, AtomicOrdering::Relaxed);
                                }

                                let hsp_q: i32 = hsp_q_u as i32;
                                let hsp_s: i32 = hsp_s_u as i32;
                                let hsp_len: i32 = (hsp_qe_u - hsp_q_u) as i32;
                                let s_last_off: i32 = s_last_off_u as i32;

                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .ungapped_extensions
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                    if right_extend {
                                        diagnostics
                                            .base
                                            .ungapped_two_hit_extensions
                                            .fetch_add(1, AtomicOrdering::Relaxed);
                                    } else {
                                        diagnostics
                                            .base
                                            .ungapped_one_hit_extensions
                                            .fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    if hsp_len > 0 {
                                        diagnostics
                                            .base
                                            .extension_total_length
                                            .fetch_add(hsp_len as usize, AtomicOrdering::Relaxed);
                                        atomic_max_usize(
                                            &diagnostics.base.extension_max_length,
                                            hsp_len as usize,
                                        );
                                    }
                                }

                                // NCBI: Update diagonal state based on right_extend
                                // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:636-648
                                // if (right_extend) {
                                //     diag_array[diag_coord].flag = 1;
                                //     diag_array[diag_coord].last_hit = s_last_off - (wordsize - 1) + diag_offset;
                                // } else {
                                //     diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                                // }
                                if right_extend {
                                    diag_entry.set_flag(1);
                                    diag_entry
                                        .set_last_hit(s_last_off - (wordsize - 1) + *diag_offset);
                                    if diag_enabled {
                                        diagnostics
                                            .base
                                            .mask_updates
                                            .fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                } else {
                                    diag_entry.set_last_hit(
                                        subject_offset.wrapping_add(*diag_offset as u32) as i32,
                                    );
                                }

                                // [C] if (score >= cutoffs->cutoff_score)
                                // NCBI reference: aa_ungapped.c:575-591 (Extension後のcutoffチェック)
                                if collect_hsp_saving_stats {
                                    if score >= cutoff {
                                        stats.score_distribution.push(score);
                                        stats.hsp_saved += 1;
                                    } else {
                                        stats.hsp_filtered_by_cutoff += 1;
                                    }
                                }
                                if score >= cutoff {
                                    if diag_enabled {
                                        diagnostics
                                            .ungapped_only_hits
                                            .fetch_add(1, AtomicOrdering::Relaxed);
                                    }

                                    // Extra debug for a traced HSP: print seed/extension inputs and cutoffs.
                                    if let Some(target) = trace_hsp_target() {
                                        // Compute outfmt coords for this candidate init-hsp (same logic as trace_init_hsp_if_match).
                                        // NCBI offsets are 0-based in query/subject->sequence buffers.
                                        // Reference: blast_gapalign.c:4756-4768, blast_aascan.c:110-113.
                                        let q_aa_start = hsp_q_u as usize;
                                        let q_aa_end = hsp_qe_u as usize;
                                        let s_aa_start = hsp_s_u as usize + chunk.offset;
                                        let s_aa_end = _hsp_se_u as usize + chunk.offset;
                                        let (q_start_dna, q_end_dna) = convert_coords(
                                            q_aa_start,
                                            q_aa_end,
                                            ctx.frame,
                                            ctx.orig_len,
                                        );
                                        let (s_start_dna, s_end_dna) = convert_coords(
                                            s_aa_start,
                                            s_aa_end,
                                            s_frame.frame,
                                            s_len,
                                        );
                                        if trace_match_target(
                                            target,
                                            q_start_dna,
                                            q_end_dna,
                                            s_start_dna,
                                            s_end_dna,
                                        ) {
                                            eprintln!(
                                                "[TRACE_HSP] seed/extend ctx_idx={} s_f_idx={} q_frame={} s_frame={} score={} cutoff={} x_dropoff={} last_hit={} subject_offset={} chunk_offset={} diff={} q_raw={} query_offset={} diag_coord={} right_extend={} s_last_off={}",
                                                ctx_idx,
                                                s_f_idx,
                                                ctx.frame,
                                                s_frame.frame,
                                                score,
                                                cutoff,
                                                x_dropoff,
                                                last_hit,
                                                subject_offset,
                                                chunk.offset,
                                                diff,
                                                q_raw,
                                                query_offset,
                                                diag_coord,
                                                right_extend,
                                                s_last_off,
                                            );
                                            // NCBI two-hit gating checks (diff/window/wordsize/context).
                                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:531-606
                                            // NCBI uses Uint4 offsets; apply unsigned wrap here as well.
                                            // Reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                                            let q_minus_diff =
                                                query_offset.wrapping_sub(diff as u32);
                                            eprintln!(
                                                "[TRACE_HSP] two_hit_pass diff={} window={} wordsize={} diff>=window={} diff<wordsize={} q_minus_diff={} ctx_frame_base={} q_minus_diff<base={} diag_offset={} diag_mask={} diag_array_size={}",
                                                diff,
                                                window,
                                                wordsize,
                                                diff >= window,
                                                diff < wordsize,
                                                q_minus_diff,
                                                ctx.frame_base,
                                                q_minus_diff < ctx.frame_base as u32,
                                                diag_offset,
                                                diag_mask,
                                                diag_array_size
                                            );
                                        }
                                    }
                                    // NCBI: BlastSaveInitHsp equivalent
                                    // Reference: blast_extend.c:360-375 BlastSaveInitHsp
                                    // Store HSP with absolute coordinates (before coordinate conversion)
                                    //
                                    // hsp_q is frame-relative coordinate in query->sequence (0-based),
                                    // frame_base is the context query_offset in the concatenated buffer.
                                    // NCBI: ungapped_data->q_start is absolute query offset.
                                    // Reference: blast_gapalign.c:4756-4768, blast_query_info.c:311-315.
                                    let hsp_q_absolute = ctx.frame_base + hsp_q;
                                    let hsp_qe_absolute = ctx.frame_base + (hsp_q + hsp_len);

                                    let init = InitHSP {
                                        q_start_absolute: hsp_q_absolute,
                                        q_end_absolute: hsp_qe_absolute,
                                        s_start: hsp_s,
                                        s_end: hsp_s + hsp_len,
                                        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:589-592
                                        // ```c
                                        // BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                                        //                  query_offset, subject_offset, hsp_len,
                                        //                  score);
                                        // ```
                                        q_seed_absolute: query_offset as i32,
                                        s_seed: subject_offset as i32,
                                        score,
                                        ctx_idx,
                                        s_f_idx,
                                        q_idx: ctx.q_idx,
                                        s_idx: s_idx as u32,
                                        q_frame: ctx.frame,
                                        s_frame: s_frame.frame,
                                        q_orig_len: ctx.orig_len,
                                        s_orig_len: s_len,
                                    };
                                    trace_init_hsp_if_match("init_hsp_saved", &init, contexts_ref);
                                    init_hsps.push(init);
                                    seq_keys.push(seq);
                                } else if diag_enabled {
                                    diagnostics
                                        .ungapped_cutoff_failed
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                    atomic_min_i32(
                                        &diagnostics.ungapped_cutoff_failed_min_score,
                                        score,
                                    );
                                    atomic_max_i32(
                                        &diagnostics.ungapped_cutoff_failed_max_score,
                                        score,
                                    );
                                }
                            }
                    }
                    }};
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-516,588-591
                // ```c
                // for (i = 0; i < hits; ++i) {
                //     Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
                //     Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
                // ...
                //     diag_coord = (query_offset - subject_offset) & diag_mask;
                // ...
                //     if (score >= cutoffs->cutoff_score)
                //         BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                //                          query_offset, subject_offset, hsp_len,
                //                          score);
                // ```
                // One flush: the buffered hits go through the loop body bucket by bucket, then the HSPs this
                // flush saved (the tail of `init_hsps`) are put back into scan order, the order in which the C
                // loop saves them. The flushes follow each other in scan order, so `init_hsps` is in scan order
                // after every flush.
                macro_rules! x_flush {
                    () => {{
                        let x_flush_start = init_hsps.len();
                        buckets.flush(|q, s, seq| x_two_hit_body!(q, s, seq));
                        crate::algorithm::tblastx::x_seed_bucket::restore_flush_order(
                            &mut init_hsps,
                            x_flush_start,
                            &mut seq_keys,
                        );
                    }};
                }
                for (_scan_chunk_index, (interior_start, interior_end)) in
                    scan_interiors.into_iter().enumerate()
                {
                    let seq_ranges = clip_tblastx_seq_ranges_for_scan_interior(
                        &base_seq_ranges,
                        interior_start,
                        interior_end,
                        wordsize as usize,
                        subject.len(),
                    );
                    if seq_ranges.is_empty() {
                        continue;
                    }
                    // [C] scan_range[0] = 0;
                    // [C] scan_range[1] = subject->seq_ranges[0].left;
                    // [C] scan_range[2] = subject->seq_ranges[0].right - wordsize;
                    let mut scan_range: [i32; 3] = [0, seq_ranges[0].0, seq_ranges[0].1 - wordsize];

                    // [C] while (scan_range[1] <= scan_range[2])
                    while scan_range[1] <= scan_range[2] {
                        let prev_scan_left = scan_range[1];
                        // [C] hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
                        let t0 = if timing_enabled {
                            Some(Instant::now())
                        } else {
                            None
                        };
                        // EXPERIMENT (LOSAT_X_TBXSCANBF): the branch-free scan with the switch.
                        let hits = crate::algorithm::tblastx::x_scan_bf::x_scan_subject(
                            x_scan_bf,
                            lookup_ref,
                            subject,
                            &seq_ranges,
                            offset_pairs,
                            offset_array_size,
                            &mut scan_range,
                        );
                        if let Some(t0) = t0 {
                            scan_ns
                                .fetch_add(t0.elapsed().as_nanos() as u64, AtomicOrdering::Relaxed);
                            scan_calls.fetch_add(1, AtomicOrdering::Relaxed);
                        }

                        if diag_enabled && hits > 0 {
                            diagnostics
                                .base
                                .kmer_matches
                                .fetch_add(hits as usize, AtomicOrdering::Relaxed);

                            // DEBUG: Check for duplicate offset pairs in scan output
                            if is_long_sequence {
                                // NCBI BlastOffsetPair uses Uint4 offsets.
                                // Reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:141-150
                                let mut seen: HashSet<(u32, u32)> =
                                    HashSet::with_capacity(hits as usize);
                                let mut duplicate_count = 0usize;
                                for i in 0..hits as usize {
                                    let pair = unsafe { &*offset_pairs.as_ptr().add(i) };
                                    if !seen.insert((pair.q_off, pair.s_off)) {
                                        duplicate_count += 1;
                                    }
                                }
                                if duplicate_count > 0 {
                                    eprintln!("[DEBUG SCAN_DUPES] s_f_idx={} scan_range=[{},{}] hits={} duplicates={} ({:.2}%)",
                                        s_f_idx, prev_scan_left, scan_range[1], hits, duplicate_count,
                                        (duplicate_count as f64 / hits as f64) * 100.0);
                                }
                            }
                        }

                        if hits == 0 && scan_range[1] == prev_scan_left {
                            // Safety guard: with correct NCBI-sized offset arrays, this should not happen.
                            // If it does, breaking avoids an infinite loop.
                            break;
                        }

                        // [C] for (i = 0; i < hits; ++i)
                        // The hits of this scan call are appended to the buckets and
                        // processed, bucket by bucket, when the buffer is full.
                        //
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-516
                        // ```c
                        // for (i = 0; i < hits; ++i) {
                        //     Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
                        //     Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
                        // ...
                        //     diag_coord = (query_offset - subject_offset) & diag_mask;
                        // ```
                        // The C loop processes each offset pair here. This loop only appends the pairs to buckets,
                        // in scan order (a buffered hit keeps its position, which is used to restore that order).
                        let offset_pairs_ptr = offset_pairs.as_ptr();
                        for i in 0..hits as usize {
                            // SAFETY: i < hits, and hits <= offset_array_size (checked by scan)
                            let pair = unsafe { &*offset_pairs_ptr.add(i) };
                            buckets.push(pair.q_off, pair.s_off);
                        }
                        if buckets.flush_due(x_bucket_budget) {
                            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-516
                            // ```c
                            // Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
                            // Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
                            // ...
                            // diag_coord = (query_offset - subject_offset) & diag_mask;
                            // ```
                            // The buffered hits are processed when the budget is reached (a bucket at a time).
                            x_flush!();
                        }
                    }
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-516
                // ```c
                // Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
                // Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
                // ...
                // diag_coord = (query_offset - subject_offset) & diag_mask;
                // ```
                // The remaining hits are processed at the end of the subject chunk.
                x_flush!();
                // The saved HSPs are in scan order here, as in C: every flush restored its own. The
                // score sort below (sort_init_hsps_by_score_ncbi, with its NCBI reference) is stable,
                // so the order of comparator-equal HSPs matters.

                // [C] Blast_ExtendWordExit(ewp, subject->length);
                //
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:614
                // ```c
                // /* increment the offset in the diagonal array */
                // Blast_ExtendWordExit(ewp, subject->length);
                // ```
                advance_tblastx_diag_offset(diag_offset, diag_array, window, chunk.length);

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:234-235
                // ```c
                // Blast_InitHitListSortByScore(init_hitlist);
                // return status;
                // ```
                // BlastAaWordFinder sorts the chunk's init hit list by
                // score_compare_match (blast_extend.c:274-310) before
                // BLAST_GetUngappedHSPList sees it.
                sort_init_hsps_by_score_ncbi(&mut init_hsps);
                let mut hits = if init_hsps.is_empty() {
                    Vec::new()
                } else {
                    // NCBI: BLAST_GetUngappedHSPList equivalent - per chunk conversion,
                    // then Blast_HSPListAdjustOffsets and Blast_HSPListsMerge.
                    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:561-584
                    // ```c
                    // BLAST_GetUngappedHSPList(init_hitlist, query_info, subject,
                    //         hit_params->options, &hsp_list);
                    // Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
                    // status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list,
                    //      kHspNumMax, &(backup.offset), INT4_MIN, overlap, ...);
                    // ```
                    get_ungapped_hsp_list(init_hsps, contexts_ref, s_frames_listed)
                };
                adjust_tblastx_chunk_subject_offsets(&mut hits, chunk.offset);
                TblastxChunkScanResult { chunk, hits, stats }
            };

            // EXPERIMENT (SEEDBUCKET): the scan order (`scan_chunk_plain`) or the bucketed order
            // (`scan_chunk_bucketed`, only for tables of at least LOSAT_X_SEEDBUCKET_MIN_CELLS
            // cells).
            //
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-516
            // ```c
            // while (scan_range[1] <= scan_range[2]) {
            // ...
            //     for (i = 0; i < hits; ++i) {
            //         Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
            //         Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
            // ...
            //         diag_coord = (query_offset - subject_offset) & diag_mask;
            // ```
            // Dispatch point of SEEDBUCKET: `scan_chunk_plain` (the port of the C loop), or
            // `scan_chunk_bucketed` for a table of at least the minimum cells. A traced or debugged
            // run always uses `scan_chunk_plain`.
            let scan_chunk = |chunk: SubjectChunk,
                              offset_pairs: &mut [OffsetPair],
                              diag_array: &mut [DiagStruct],
                              diag_offset: &mut i32,
                              scan_chunk_size: Option<usize>,
                              sc: S|
             -> TblastxChunkScanResult {
                // The debug/trace paths exist only in the reference loop.
                let x_debugging = scan_debug_range.is_some() || trace_hsp_target().is_some();
                if crate::algorithm::tblastx::x_seed_bucket::bucketed(diag_array.len())
                    && !x_debugging
                {
                    scan_chunk_bucketed(
                        chunk,
                        offset_pairs,
                        diag_array,
                        diag_offset,
                        scan_chunk_size,
                        sc,
                    )
                } else {
                    scan_chunk_plain(
                        chunk,
                        offset_pairs,
                        diag_array,
                        diag_offset,
                        scan_chunk_size,
                        sc,
                    )
                }
            };

            // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:492-505
            // while (scan_range[1] <= scan_range[2]) { hits = scansub(...); }
            if use_serial_scan_chunks {
                loop {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:452-584
                    // ```c
                    // while (TRUE) {
                    //    status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
                    //                                   dbseq_chunk_overlap);
                    //    if (status == SUBJECT_SPLIT_DONE) break;
                    //    BlastInitHitListReset(init_hitlist);
                    //    aux_struct->WordFinder(..., init_hitlist, ...);
                    //    BLAST_GetUngappedHSPList(..., &hsp_list);
                    //    Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
                    //    status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list, ...);
                    // }
                    // ```
                    let chunk = match split_state.next_chunk(max_dbseq_len, DBSEQ_CHUNK_OVERLAP) {
                        SubjectChunkStatus::Done => break,
                        SubjectChunkStatus::Ok(chunk) => chunk,
                    };

                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:83-123
                    // ```c
                    // for (s = s_first; s <= s_last; s++) {
                    //     ...
                    //     offset_pairs[i + totalhits].qs_offsets.s_off = s_off;
                    // }
                    // ```
                    //
                    // Experimental scan chunks split only the scan walk inside this
                    // NCBI search unit. The reducer still extends against `subject`,
                    // the full real chunk, and `Blast_ExtendWordExit` runs once below.
                    let scan_size = tblastx_scan_chunk_size_for_run(chunk.length, num_threads);
                    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:492-505
                    // while (scan_range[1] <= scan_range[2]) { hits = scansub(...); }
                    crate::utils::threading::report_stage(
                        "tblastx",
                        "serial_scan_chunks",
                        chunk.length.div_ceil(scan_size),
                        false,
                    );
                    let result = scan_chunk(
                        chunk,
                        offset_pairs,
                        diag_array,
                        &mut diag_offset,
                        Some(scan_size),
                        sc,
                    );
                    stats_hsp_saved += result.stats.hsp_saved;
                    stats_hsp_filtered_by_cutoff += result.stats.hsp_filtered_by_cutoff;
                    stats_score_distribution.extend(result.stats.score_distribution);
                    if !result.hits.is_empty() {
                        merge_tblastx_subject_chunk_hits(
                            &mut frame_ungapped_hits,
                            result.hits,
                            result.chunk.offset,
                            result.chunk.overlap,
                        );
                    }
                }
            } else if use_parallel_chunks {
                #[cfg(all(
                    feature = "parallel",
                    any(not(target_arch = "wasm32"), feature = "wasm-threads")
                ))]
                {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:452-584
                    // ```c
                    // while (TRUE) {
                    //    status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
                    //                                   dbseq_chunk_overlap);
                    //    if (status == SUBJECT_SPLIT_DONE) break;
                    //    ...
                    //    status = Blast_HSPListsMerge(...);
                    // }
                    // ```
                    let mut chunks = Vec::new();
                    loop {
                        match split_state.next_chunk(max_dbseq_len, DBSEQ_CHUNK_OVERLAP) {
                            SubjectChunkStatus::Done => break,
                            SubjectChunkStatus::Ok(chunk) => chunks.push(chunk),
                        }
                    }
                    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:452-584
                    // status = s_GetNextSubjectChunk(subject, &backup, ...);
                    crate::utils::threading::report_stage(
                        "tblastx",
                        "subject_chunks",
                        chunks.len(),
                        chunks.len() > 1,
                    );
                    if chunks.len() <= 1 {
                        for chunk in chunks {
                            let result = scan_chunk(
                                chunk,
                                offset_pairs,
                                diag_array,
                                &mut diag_offset,
                                None,
                                sc,
                            );
                            stats_hsp_saved += result.stats.hsp_saved;
                            stats_hsp_filtered_by_cutoff += result.stats.hsp_filtered_by_cutoff;
                            stats_score_distribution.extend(result.stats.score_distribution);
                            if !result.hits.is_empty() {
                                merge_tblastx_subject_chunk_hits(
                                    &mut frame_ungapped_hits,
                                    result.hits,
                                    result.chunk.offset,
                                    result.chunk.overlap,
                                );
                            }
                        }
                    } else {
                        let mut chunk_results: Vec<TblastxChunkScanResult> =
                            Vec::with_capacity(chunks.len());
                        chunks
                            .par_iter()
                            .map_init(
                                || {
                                    (
                                        vec![OffsetPair::default(); offset_pairs_len],
                                        vec![DiagStruct::default(); diag_array_size as usize],
                                        sc,
                                    )
                                },
                                |state, chunk| {
                                    let (offset_pairs, diag_array, sc) = state;
                                    for diag in diag_array.iter_mut() {
                                        *diag = DiagStruct::default();
                                    }
                                    let mut chunk_diag_offset = window;
                                    scan_chunk(
                                        *chunk,
                                        offset_pairs,
                                        diag_array,
                                        &mut chunk_diag_offset,
                                        None,
                                        *sc,
                                    )
                                },
                            )
                            .collect_into_vec(&mut chunk_results);
                        chunk_results.sort_by_key(|result| result.chunk.offset);
                        for result in chunk_results {
                            // Keep the canonical per-subject diagonal offset moving in the
                            // same chunk order as NCBI before the next subject frame starts.
                            // The experimental workers above use local diagonal tables; the
                            // ordered merge below remains the NCBI `Blast_HSPListsMerge` path.
                            // References: blast_engine.c:561-584, blast_extend.c:167-173.
                            advance_tblastx_diag_offset(
                                &mut diag_offset,
                                diag_array,
                                window,
                                result.chunk.length,
                            );
                            stats_hsp_saved += result.stats.hsp_saved;
                            stats_hsp_filtered_by_cutoff += result.stats.hsp_filtered_by_cutoff;
                            stats_score_distribution.extend(result.stats.score_distribution);
                            if !result.hits.is_empty() {
                                merge_tblastx_subject_chunk_hits(
                                    &mut frame_ungapped_hits,
                                    result.hits,
                                    result.chunk.offset,
                                    result.chunk.overlap,
                                );
                            }
                        }
                    }
                }
                #[cfg(any(
                    not(feature = "parallel"),
                    all(target_arch = "wasm32", not(feature = "wasm-threads"))
                ))]
                {
                    unreachable!("parallel chunk mode is disabled for this target");
                }
            } else {
                loop {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:452-584
                    // ```c
                    // while (TRUE) {
                    //    status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
                    //                                   dbseq_chunk_overlap);
                    //    if (status == SUBJECT_SPLIT_DONE) break;
                    //    BlastInitHitListReset(init_hitlist);
                    //    aux_struct->WordFinder(..., init_hitlist, ...);
                    //    BLAST_GetUngappedHSPList(..., &hsp_list);
                    //    Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
                    //    status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list, ...);
                    // }
                    // ```
                    let chunk = match split_state.next_chunk(max_dbseq_len, DBSEQ_CHUNK_OVERLAP) {
                        SubjectChunkStatus::Done => break,
                        SubjectChunkStatus::Ok(chunk) => chunk,
                    };

                    // NCBI: BlastInitHitListReset(init_hitlist) - reset per chunk.
                    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:488-491
                    let result =
                        scan_chunk(chunk, offset_pairs, diag_array, &mut diag_offset, None, sc);
                    stats_hsp_saved += result.stats.hsp_saved;
                    stats_hsp_filtered_by_cutoff += result.stats.hsp_filtered_by_cutoff;
                    stats_score_distribution.extend(result.stats.score_distribution);
                    if !result.hits.is_empty() {
                        merge_tblastx_subject_chunk_hits(
                            &mut frame_ungapped_hits,
                            result.hits,
                            result.chunk.offset,
                            result.chunk.overlap,
                        );
                    }
                }
            }

            if !frame_ungapped_hits.is_empty() {
                // NCBI: Blast_HSPListAppend merges per-frame subject HSP lists
                // into the combined translated-subject list, then sorts the
                // combined list by score through s_BlastHSPListsCombineByScore.
                //
                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:804-845
                // ```c
                // for (context=first_context; context<=last_context; context++) {
                //     subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
                //     ...
                //     Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax);
                // }
                // ```
                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2758-2766
                // ```c
                // for (index=combined_hsp_list->hspcnt, index1=0;
                //      index1<hsp_list->hspcnt; index1++) {
                //    combined_hsp_list->hsp_array[index++] = hsp_list->hsp_array[index1];
                // }
                // combined_hsp_list->hspcnt = new_hspcnt;
                // Blast_HSPListSortByScore(combined_hsp_list);
                // ```
                combined_ungapped_hits.extend(frame_ungapped_hits);
                if !ungapped_hits_is_sorted_by_score_ncbi(&combined_ungapped_hits) {
                    sort_ungapped_hits_by_score_ncbi(&mut combined_ungapped_hits);
                }
            }
        } // End of subject frame loop
        x_lap.lap(XStage::FrameLoop);

        // Dispatch point 5 of the seed release (per subject): the seed stage of
        // this subject ends here; nothing below reads the batch's lookup table (`end` consumes
        // the stage that `lookup_ref` borrows from; the unit path's stage ended at dispatch
        // point 3) or the worker's diagonal table. The batch's last lookup
        // stage frees the lookup table and returns true, and this worker, which gets no
        // further subject, drops its diagonal table before the linking call below, where the
        // heap of the search peaks.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-873
        // ```c
        // if (hit_params->link_hsp_params) {
        //     status = BLAST_LinkHsps(program_number, hsp_list_out, query_info,
        //               subject->length, gap_align->sbp, hit_params->link_hsp_params,
        //               score_options->gapped_calculation);
        // ```
        // NCBI links with its lookup table and diagonal table still allocated (freed after the
        // subject loop, blast_engine.c:1605); `BLAST_LinkHsps` reads neither.
        if x_seed_stage.end() {
            x_seed_release::drop_worker_table(&mut st.diag_array);
        }
        // EXPERIMENT (LOSAT_X_TBXPIPE): the seed stage of this subject has ended with the HSP list
        // that the linking below takes; the batch's last one may start the next batch. While a
        // later batch may run, the two linking calls below do not use the pool (`x_link_on_pool`,
        // tblastx/x_batch_pipe.rs item 6); the lists are the same with or without it.
        let x_link_on_pool = use_parallel && x_pipe_seed.end(combined_ungapped_hits.len());
        x_lap.lap(XStage::SeedRelease);

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:561-584
        // ```c
        // BLAST_GetUngappedHSPList(init_hitlist, query_info, subject,
        //         hit_params->options, &hsp_list);
        // Blast_HSPListsMerge(...);
        // ```
        // This snapshot records internal frame-relative HSPs immediately after
        // the ungapped extension list is materialized. TBLASTX has no active
        // common-endpoint purge in this ungapped path, so the post-purge
        // diagnostic intentionally records the same NCBI-stage boundary.
        if stage_dump::enabled() {
            stage_dump::dump_ungapped_hits(
                "after_initial_ungapped_extension",
                &combined_ungapped_hits,
            );
            stage_dump::dump_ungapped_hits("after_common_endpoint_purge", &combined_ungapped_hits);
        }

        // Build NCBI-style subject frame base offsets for sum-statistics linking.
        // In NCBI, HSP coords live in a concatenated translation buffer with sentinels.
        // LOSAT uses per-frame sequences; for linking we emulate absolute offsets by
        // concatenating frames in the same order as `generate_frames()`.
        let mut subject_frame_bases: Vec<i32> = Vec::with_capacity(s_frames_listed.len());
        let mut base: i32 = 0;
        for f in s_frames_listed {
            subject_frame_bases.push(base);
            // NCBI concatenation shares the trailing NULLB sentinel between frames:
            //   offset += length + 1;
            // where `length` is the number of residues (excluding sentinels).
            // Source: ncbi-blast/c++/src/algo/blast/core/blast_util.c:1098-1101
            base += f.aa_seq.len() as i32 - 1;
        }
        // Nothing below reads the preliminary frames (the seed stage has ended); they are freed
        // before the linking call, where the heap peaks.
        drop(s_frames_preliminary_owned);

        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:870-899
        // ```c
        // if (hit_params->link_hsp_params) {
        //     status = BLAST_LinkHsps(program_number, hsp_list_out, query_info,
        //               subject->length, gap_align->sbp, hit_params->link_hsp_params,
        //               score_options->gapped_calculation);
        // }
        // ...
        // status = s_Blast_HSPListReapByPrelimEvalue(hsp_list_out, hit_params);
        // ```
        // NCBI links and prelim-evalue-reaps the raw ungapped HSP list inside
        // s_BlastSearchEngineCore before the outer translated-subject path
        // reevaluates ambiguities and calls BLAST_LinkHsps again.
        x_link_timer.arm(LinkCall::First);
        let mut prelinked_ungapped_hits = if combined_ungapped_hits.is_empty() {
            combined_ungapped_hits
        } else {
            apply_sum_stats_even_gap_linking_with_parallel(
                combined_ungapped_hits,
                &linking_params_for_cutoff,
                &linking_params,
                contexts_ref,
                &subject_frame_bases,
                &length_adj_per_context,
                &eff_searchsp_per_context,
                x_link_on_pool,
            )
        };
        x_link_timer.disarm();
        x_lap.lap(XStage::Link1);
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c:1802-1803
        // ```c
        //     /* Sort the HSP array by score */
        //     Blast_HSPListSortByScore(hsp_list);
        // ```
        // BLAST_LinkHsps ends by sorting the list it linked, here with the scores of the
        // preliminary (compressed-subject) search. The re-evaluation below may change or
        // trim those scores (random bases for ambiguous subject letters, two-hit
        // extension tails), after which Blast_HSPListReevaluateUngapped sorts again
        // (blast_hits.c:2733-2734): HSPs that tie after re-evaluation keep the order
        // of this first sort, not the order of the linking chains. That order is the
        // input order of the second BLAST_LinkHsps (blast_engine.c:1515-1520), whose
        // `>=` chain choices (link_hsps.c:613-622, 759, 887) give the better chain to
        // the later of two comparator-equal HSPs of different frames.
        // (The sort at the end of the second BLAST_LinkHsps is the score sort of each
        // subject's list in `hit_collector::QueryLists::add`.)
        if !ungapped_hits_is_sorted_by_score_ncbi(&prelinked_ungapped_hits) {
            sort_ungapped_hits_by_score_ncbi(&mut prelinked_ungapped_hits);
        }
        if stage_dump::enabled() {
            stage_dump::dump_ungapped_hits(
                "after_prelim_link_hsps_before_reevaluate",
                &prelinked_ungapped_hits,
            );
        }
        let prelim_linked_count = prelinked_ungapped_hits.len();
        // NCBI reference: c++/src/algo/blast/core/blast_hits.c:1988-1996
        // ```c
        //    cutoff = hit_options->expect_value;
        // ...
        //       if (hsp->evalue > cutoff) {
        // ```
        // An HSP is removed when its e-value is greater than the cutoff, so a NaN cutoff
        // (`-evalue -nan`) keeps every HSP.
        prelinked_ungapped_hits.retain(|h| !(h.e_value > evalue_threshold));
        x_lap.lap(XStage::Reap);
        if diag_enabled && prelim_linked_count != prelinked_ungapped_hits.len() {
            eprintln!(
                "[DEBUG PRELIM_REAP] linked_before={} kept={} filtered_by_prelim_evalue={} threshold={}",
                prelim_linked_count,
                prelinked_ungapped_hits.len(),
                prelim_linked_count - prelinked_ungapped_hits.len(),
                evalue_threshold
            );
        }
        if stage_dump::enabled() {
            stage_dump::dump_ungapped_hits(
                "after_prelim_evalue_reap_before_reevaluate",
                &prelinked_ungapped_hits,
            );
        }

        // NCBI: Blast_HSPListReevaluateUngapped equivalent
        // Reference: blast_engine.c:1492-1497, blast_hits.c:2609-2737
        // Perform batch reevaluation on all HSPs after merging all frames
        //
        // A subject with ambiguous bases gets the frames of its letters here, after its
        // preliminary frames were freed, when HSPs are left for the re-evaluation (none are read
        // otherwise).
        let x_letter_frames: Vec<QueryFrame>;
        let (s_frames, s_frames_report): (&[QueryFrame], &[QueryFrame]) = if x_letter_frames_later {
            x_letter_frames = if prelinked_ungapped_hits.is_empty() {
                Vec::new()
            } else {
                subjects.x_letter_frames(s_idx, &db_code)
            };
            x_lap.lap(XStage::Translate);
            (&x_letter_frames, &x_letter_frames)
        } else {
            (&s_frames, s_frames_report)
        };
        let mut ungapped_hits = reevaluate_ungapped_hsp_list(
            sc,
            prelinked_ungapped_hits,
            contexts_ref,
            s_frames,
            s_frames_report,
            &cutoff_scores,
            timing_enabled,
            &reeval_ns,
            &reeval_calls,
            collect_hsp_saving_stats,
            &mut stats_hsp_filtered_by_reeval,
        );
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2733-2734
        // ```c
        // /* Sort the HSP array by score (scores may have changed!) */
        // Blast_HSPListSortByScore(hsp_list);
        // ```
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1374-1381
        // ```c
        // if (!Blast_HSPListIsSortedByScore(hsp_list)) {
        //     qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
        //           ScoreCompareHSPs);
        // }
        // ```
        // This sort happens before the second BLAST_LinkHsps call in
        // blast_engine.c:1515-1520 and controls tie order for link_hsps.c
        // comparator-equal short HSPs.
        if !ungapped_hits_is_sorted_by_score_ncbi(&ungapped_hits) {
            sort_ungapped_hits_by_score_ncbi(&mut ungapped_hits);
        }
        x_lap.lap(XStage::Reevaluate);
        // Record the post-reevaluation list and the exact input to link_hsps.
        if stage_dump::enabled() {
            stage_dump::dump_ungapped_hits("after_reevaluate", &ungapped_hits);
            stage_dump::dump_ungapped_hits("before_link_hsps", &ungapped_hits);
        }

        if !ungapped_hits.is_empty() {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3265-3273
            // ```c
            // #if _BLAST_DEBUG
            // arg_desc.AddFlag("verbose", "Produce verbose output (show BLAST options)",
            //                  true);
            // #endif /* _BLAST_DEBUG */
            // ```
            // DEBUG: Print HSP statistics for long sequences only when requested.
            if collect_hsp_saving_stats {
                eprintln!(
                    "[DEBUG HSP_STATS] After reevaluation: {} HSPs",
                    ungapped_hits.len()
                );
                eprintln!("[DEBUG HSP_STATS] Saved by cutoff: {}", stats_hsp_saved);
                eprintln!(
                    "[DEBUG HSP_STATS] Filtered by cutoff: {}",
                    stats_hsp_filtered_by_cutoff
                );
                eprintln!(
                    "[DEBUG HSP_STATS] Filtered by reeval: {}",
                    stats_hsp_filtered_by_reeval
                );
                if !stats_score_distribution.is_empty() {
                    let min_score = stats_score_distribution.iter().min().unwrap();
                    let max_score = stats_score_distribution.iter().max().unwrap();
                    let avg_score = stats_score_distribution.iter().sum::<i32>() as f64
                        / stats_score_distribution.len() as f64;
                    eprintln!(
                        "[DEBUG HSP_STATS] Score range: {} - {} (avg: {:.2})",
                        min_score, max_score, avg_score
                    );
                }
            }
            if diag_enabled {
                diagnostics
                    .base
                    .hsps_before_chain
                    .fetch_add(ungapped_hits.len(), AtomicOrdering::Relaxed);
            }
            // Output HSP saving statistics for long sequences
            if collect_hsp_saving_stats
                && (stats_hsp_saved > 0
                    || stats_hsp_filtered_by_cutoff > 0
                    || stats_hsp_filtered_by_reeval > 0
                    || stats_hsp_filtered_by_hsp_test > 0)
            {
                let total_attempted = stats_hsp_saved
                    + stats_hsp_filtered_by_cutoff
                    + stats_hsp_filtered_by_reeval
                    + stats_hsp_filtered_by_hsp_test;
                eprintln!(
                    "[DEBUG HSP_SAVING] subject_len_nucl={}, cutoff={}",
                    subject_len_nucl, cutoff_score_min
                );
                eprintln!("[DEBUG HSP_SAVING] total_attempted={}, saved={}, filtered_by_cutoff={}, filtered_by_reeval={}, filtered_by_hsp_test={}", 
                        total_attempted, stats_hsp_saved, stats_hsp_filtered_by_cutoff, stats_hsp_filtered_by_reeval, stats_hsp_filtered_by_hsp_test);
                if !stats_score_distribution.is_empty() {
                    stats_score_distribution.sort();
                    let min_score = stats_score_distribution[0];
                    let max_score = stats_score_distribution[stats_score_distribution.len() - 1];
                    let median_score = stats_score_distribution[stats_score_distribution.len() / 2];
                    let low_score_count =
                        stats_score_distribution.iter().filter(|&&s| s < 30).count();
                    eprintln!("[DEBUG HSP_SAVING] score_distribution: min={}, max={}, median={}, low_score(<30)={}/{} ({:.2}%)", 
                            min_score, max_score, median_score, low_score_count, stats_score_distribution.len(),
                            if stats_score_distribution.len() > 0 { (low_score_count as f64 / stats_score_distribution.len() as f64) * 100.0 } else { 0.0 });
                }
            }

            // NCBI Parity: Pass pre-computed length_adjustment and eff_searchsp per context
            // These values are stored in query_info->contexts[ctx] in NCBI and referenced
            // by link_hsps.c for BLAST_SmallGapSumE/BLAST_LargeGapSumE calculations.
            let t_linking = if timing_enabled {
                Some(Instant::now())
            } else {
                None
            };
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:553-558
            // ```c
            // for (frame_index=0; frame_index<num_query_frames; frame_index++)
            // {
            //    hp_start->next = hp_frame_start[frame_index];
            //    number_of_hsps = hp_frame_number[frame_index];
            // }
            // ```
            x_link_timer.arm(LinkCall::Second);
            let linked = apply_sum_stats_even_gap_linking_with_parallel(
                ungapped_hits,
                &linking_params_for_cutoff,
                &linking_params,
                contexts_ref,
                &subject_frame_bases,
                &length_adj_per_context,
                &eff_searchsp_per_context,
                x_link_on_pool,
            );
            x_link_timer.disarm();
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:959-982
            // ```c
            // H->hsp->evalue = (best_evalue == -1) ? H->hsp->evalue :
            //                  MIN(H->hsp->evalue, best_evalue);
            // H->ordering_method = ordering_method;
            // ```
            // Record linked_set/start_of_chain/order/e-value before output
            // coordinate conversion can obscure frame-relative HSP identity.
            if stage_dump::enabled() {
                stage_dump::dump_ungapped_hits("after_link_hsps_before_output_conversion", &linked);
            }
            if trace_hsp_target().is_some() {
                for h in &linked {
                    trace_ungapped_hit_if_match("after_linking", h);
                }
            }
            if let Some(t) = t_linking {
                let elapsed = t.elapsed();
                linking_ns.fetch_add(elapsed.as_nanos() as u64, AtomicOrdering::Relaxed);
                linking_calls.fetch_add(1, AtomicOrdering::Relaxed);
            }
            x_lap.lap(XStage::Link2);
            if diag_enabled {
                diagnostics
                    .base
                    .hsps_after_chain
                    .fetch_add(linked.len(), AtomicOrdering::Relaxed);
            }

            // NCBI's culling (`-culling_limit`) runs on the final HSPs of all subjects, after
            // the e-value reap (`report::culled_hit_order`).

            let total_linked = linked.len();
            let mut stats_single_hsps = 0usize;
            let mut stats_chain_heads = 0usize;
            let mut stats_chain_members = 0usize;
            if debug_output_filter {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3265-3273
                // ```c
                // #if _BLAST_DEBUG
                // arg_desc.AddFlag("verbose", "Produce verbose output (show BLAST options)",
                //                  true);
                // #endif /* _BLAST_DEBUG */
                // ```
                // DEBUG: Collect statistics before filtering only when requested.
                for h in &linked {
                    if h.linked_set && !h.start_of_chain {
                        stats_chain_members += 1;
                    } else if !h.linked_set {
                        stats_single_hsps += 1;
                    } else if h.start_of_chain {
                        stats_chain_heads += 1;
                    }
                }
            }

            let mut final_hits: Vec<TblastxHsp> = Vec::new();
            let dump_output_stage = stage_dump::enabled();
            let mut output_snapshot_hits: Vec<UngappedHit> = Vec::new();
            let mut output_snapshot_pairs: Vec<(Hit, UngappedHit)> = Vec::new();
            let mut filtered_by_evalue = 0usize;
            for h in linked {
                // NCBI reference (verbatim, link_hsps.c:1018-1020):
                //   /* If this is not a single piece or the start of a chain, then Skip it. */
                //   if (H->linked_set == TRUE && H->start_of_chain == FALSE)
                //       continue;
                // NOTE: This "continue" in NCBI is NOT a filter - NCBI then walks the link
                // pointer from each chain head to include all chain members (lines 1047-1059).
                // LOSAT's linking already produces a flat list of ALL HSPs (including chain
                // members) so we output everything without this skip logic.

                // NCBI reference: E-value filtering is applied during output conversion
                // The exact timing may differ, but the threshold check is standard
                if h.e_value > evalue_threshold {
                    filtered_by_evalue += 1;
                    if diag_enabled {
                        diagnostics
                            .ungapped_evalue_failed
                            .fetch_add(1, AtomicOrdering::Relaxed);
                    }
                    continue;
                }

                if diag_enabled {
                    diagnostics
                        .ungapped_evalue_passed
                        .fetch_add(1, AtomicOrdering::Relaxed);
                }
                if dump_output_stage {
                    output_snapshot_hits.push(h.clone());
                }

                let ctx = &contexts_ref[h.ctx_idx as usize];
                let s_score_frame = &s_frames[usize::from(h.s_f_idx)];
                let s_frame = &s_frames_report[usize::from(h.s_f_idx)];

                let len = (h.q_aa_end as usize).saturating_sub(h.q_aa_start as usize);
                let q0 = h.q_aa_start as usize;
                let s0 = h.s_aa_start as usize;

                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2708-2720
                // ```c
                // Blast_HSPGetNumIdentitiesAndPositives(query_nomask,
                //     subject_start, hsp, score_params->options, &align_length, sbp);
                // delete_hsp = Blast_HSPTest(hsp, ...);
                // ```
                // `reevaluate_ungapped_hsp_list` already performs this NCBI step
                // using the unmasked query and reporting subject buffers, then stores
                // the identity count on the HSP. Reuse it here instead of rescanning
                // every final hit.
                let matches = h.num_ident as usize;
                let mismatch = len.saturating_sub(matches);
                let identity = if len > 0 {
                    (matches as f64 / len as f64) * 100.0
                } else {
                    0.0
                };

                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1918-1928
                // ```c
                // kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
                // hsp->bit_score =
                //    (hsp->score*kbp[hsp->context]->Lambda - kbp[hsp->context]->logK) /
                //    NCBIMATH_LN2;
                // ```
                let bit_params = &contexts_ref[h.ctx_idx as usize].karlin_params;
                let bit = calc_bit_score(h.raw_score, bit_params);
                let (q_start, q_end) = convert_coords(
                    h.q_aa_start as usize,
                    h.q_aa_end as usize,
                    ctx.frame,
                    ctx.orig_len,
                );
                let (s_start, s_end) = convert_coords(
                    h.s_aa_start as usize,
                    h.s_aa_end as usize,
                    s_frame.frame,
                    s_len,
                );

                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_filter.c:1379-1404
                // ```c
                // query_blk->sequence_start_nomask = BlastMemDup(query_blk->sequence_start, total_length);
                // query_blk->sequence_nomask = query_blk->sequence_start_nomask + 1;
                // Blast_MaskTheResidues(buffer, query_length, kIsNucl, ...);
                // ```
                // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:699-700
                // ```c
                // sum += matrix[*query & kResidueMask][*subject];
                // query++;
                // ```
                // Diagnostic only: compare the traced HSP's masked working-query
                // score against the preserved unmasked query copy used for identity
                // reporting. Normal runtime behavior is unchanged unless both
                // LOSAT_TRACE_HSP and LOSAT_TRACE_HSP_MASKS are set.
                if crate::utils::xdiag::var_os("LOSAT_TRACE_HSP_MASKS").is_some() {
                    if let Some(target) = trace_hsp_target() {
                        if trace_match_target(target, q_start, q_end, s_start, s_end) {
                            // NCBI uses query_nomask = query_blk->sequence_nomask + query_offset,
                            // with sequence_nomask pointing past the leading NULLB.
                            // Reference: blast_filter.c:1381-1382, blast_util.c:112-116.
                            let q_seq_nomask_full: &[u8] =
                                ctx.aa_seq_nomask.as_deref().unwrap_or(&ctx.aa_seq);
                            let q_seq_nomask = &q_seq_nomask_full[1..q_seq_nomask_full.len() - 1];
                            let s_seq = &s_frame.aa_seq[1..s_frame.aa_seq.len() - 1];
                            let q_seq_masked = &ctx.aa_seq[1..ctx.aa_seq.len() - 1];
                            let s_seq_scoring =
                                &s_score_frame.aa_seq[1..s_score_frame.aa_seq.len() - 1];
                            let mut masked_residues = 0usize;
                            let mut masked_score = 0i32;
                            let mut unmasked_score = 0i32;
                            let mut report_unmasked_score = 0i32;
                            let mut mask_runs = Vec::new();
                            let mut current_mask_start: Option<usize> = None;
                            for rel in 0..len {
                                let q_pos = q0 + rel;
                                let s_pos = s0 + rel;
                                if q_pos >= q_seq_masked.len()
                                    || q_pos >= q_seq_nomask.len()
                                    || s_pos >= s_seq.len()
                                    || s_pos >= s_seq_scoring.len()
                                {
                                    break;
                                }
                                let q_masked = q_seq_masked[q_pos];
                                let q_unmasked = q_seq_nomask[q_pos];
                                let subject_scoring = s_seq_scoring[s_pos];
                                let subject_report = s_seq[s_pos];
                                if q_masked != q_unmasked {
                                    masked_residues += 1;
                                    if current_mask_start.is_none() {
                                        current_mask_start = Some(rel);
                                    }
                                } else if let Some(start) = current_mask_start.take() {
                                    mask_runs.push(format!("{}..{}", start, rel));
                                }
                                masked_score += sc.score(q_masked, subject_scoring);
                                unmasked_score += sc.score(q_unmasked, subject_scoring);
                                report_unmasked_score += sc.score(q_unmasked, subject_report);
                            }
                            if let Some(start) = current_mask_start.take() {
                                mask_runs.push(format!("{}..{}", start, len));
                            }
                            eprintln!(
                                "[TRACE_HSP_MASKS] q={}-{} s={}-{} ctx_idx={} q_frame={} s_frame={} len={} raw_score={} masked_score={} unmasked_score={} report_unmasked_score={} masked_residues={} mask_runs={}",
                                q_start,
                                q_end,
                                s_start,
                                s_end,
                                h.ctx_idx,
                                ctx.frame,
                                s_frame.frame,
                                len,
                                h.raw_score,
                                masked_score,
                                unmasked_score,
                                report_unmasked_score,
                                masked_residues,
                                if mask_runs.is_empty() {
                                    "none".to_string()
                                } else {
                                    mask_runs.join(",")
                                }
                            );
                        }
                    }
                }

                // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
                // ```c
                // typedef struct BlastHSPList {
                //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
                //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
                //                       Set to 0 if not applicable */
                // } BlastHSPList;
                // ```
                let out_hit = Hit {
                    identity,
                    length: len,
                    mismatch: mismatch as i64,
                    gapopen: 0,
                    q_start,
                    q_end,
                    s_start,
                    s_end,
                    e_value: h.e_value,
                    bit_score: bit,
                    num_ident: matches,
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
                    // ```c
                    // if (hsp->query.frame != hsp->subject.frame) {
                    //    *q_end = query_length - hsp->query.offset;
                    //    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
                    // }
                    // ```
                    query_frame: ctx.frame as i32,
                    query_length: 0,
                    q_idx: ctx.q_idx,
                    s_idx: h.s_idx,
                    raw_score: h.raw_score,
                    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4754-4767
                    // ```c
                    // s_AdjustInitialHSPOffsets(init_hsp,
                    //                           query_info->contexts[context].query_offset);
                    // Blast_HSPInit(ungapped_data->q_start, ...,
                    //               ungapped_data->s_start, ...,
                    //               context, query_info->contexts[context].frame,
                    //               subject->frame, ungapped_data->score, NULL, &new_hsp);
                    // ```
                    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:808-812
                    // ```c
                    // subject->sequence = translation_buffer + frame_offsets[context] + 1;
                    // subject->length = frame_offsets[context+1] - frame_offsets[context] - 1;
                    // ```
                    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1347-1353
                    // ```c
                    // BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)
                    // BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)
                    // BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset)
                    // BLAST_CMP(hsp2->query.end, hsp1->query.end)
                    // ```
                    // TBLASTX ScoreCompareHSPs uses the context/frame-relative
                    // BlastSeg offsets stored in Blast_HSPInit, not formatted
                    // nucleotide coordinates or LOSAT's concatenated frame bases.
                    sort_query_offset: h.q_aa_start as usize,
                    sort_query_end: h.q_aa_end as usize,
                    sort_subject_offset: h.s_aa_start as usize,
                    sort_subject_end: h.s_aa_end as usize,
                    has_sort_offsets: true,
                    gap_info: None,
                    num_positives: matches,
                };
                trace_final_hit_if_match("output_hit", &out_hit);
                if dump_output_stage {
                    output_snapshot_pairs.push((out_hit.clone(), h.clone()));
                }
                // NCBI reference: c++/src/algo/blast/core/blast_gapalign.c:4760-4767
                // ```c
                //       Blast_HSPInit(ungapped_data->q_start,
                //                     ungapped_data->length+ungapped_data->q_start,
                //     ...
                //                     context, query_info->contexts[context].frame,
                //                     subject->frame, ungapped_data->score, NULL, &new_hsp);
                // ```
                // The HSP keeps its subject frame and its linked-set size (`num`) for the
                // pairwise report.
                final_hits.push(TblastxHsp {
                    hit: out_hit,
                    subject_frame: s_frame.frame,
                    num: h.num,
                });
            }
            x_lap.lap(XStage::Convert);
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:1018-1059
            // ```c
            // if (H->linked_set == TRUE && H->start_of_chain == FALSE)
            //     continue;
            // while (H->hsp_link.link[ordering_method]) { ... }
            // ```
            // The converted output hits are represented here by their original
            // internal HSP rows after e-value filtering and before common.rs
            // applies final HSP-list sorting for the selected output format.
            if dump_output_stage {
                stage_dump::dump_ungapped_hits(
                    "after_output_conversion_before_final_sort",
                    &output_snapshot_hits,
                );
                stage_dump::dump_final_output_order(
                    "after_final_output_sort",
                    &output_snapshot_pairs,
                );
            }

            if debug_output_filter {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3265-3273
                // ```c
                // #if _BLAST_DEBUG
                // arg_desc.AddFlag("verbose", "Produce verbose output (show BLAST options)",
                //                  true);
                // #endif /* _BLAST_DEBUG */
                // ```
                // DEBUG: Print output filtering statistics.
                // NCBI reference: link_hsps.c:1018-1020 - continue is NOT an output filter
                // Chain members are included in output via link pointer traversal.
                eprintln!("[DEBUG OUTPUT_FILTER] Total linked HSPs: {}", total_linked);
                eprintln!(
                    "[DEBUG OUTPUT_FILTER] Single HSPs (linked_set=false): {}",
                    stats_single_hsps
                );
                eprintln!(
                    "[DEBUG OUTPUT_FILTER] Chain heads (linked_set=true, start_of_chain=true): {}",
                    stats_chain_heads
                );
                eprintln!("[DEBUG OUTPUT_FILTER] Chain members (linked_set=true, start_of_chain=false): {} (included in output)", stats_chain_members);
                eprintln!(
                    "[DEBUG OUTPUT_FILTER] Filtered by E-value (threshold={}): {}",
                    evalue_threshold, filtered_by_evalue
                );
                eprintln!("[DEBUG OUTPUT_FILTER] Expected after E-value filter: {} (all HSPs - E-value filtered)", total_linked - filtered_by_evalue);
                eprintln!(
                    "[DEBUG OUTPUT_FILTER] Final hits after filtering: {}",
                    final_hits.len()
                );
            }

            if !final_hits.is_empty() {
                if diag_enabled {
                    diagnostics
                        .base
                        .hsps_after_overlap_filter
                        .fetch_add(final_hits.len(), AtomicOrdering::Relaxed);
                    diagnostics
                        .output_from_ungapped
                        .fetch_add(final_hits.len(), AtomicOrdering::Relaxed);
                }
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1497
                // ```c
                // while (...) {
                //    ...
                //    status = s_BlastSearchEngineCore(..., &hsp_list, ...);
                // }
                // ```
                // The subject's final HSPs, sent at the end of the subject (`emit_subject_hits`).
                debug_assert!(st.hits.is_empty());
                st.hits = final_hits;
            }
            x_lap.lap(XStage::Emit);
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1475
        // ```c
        // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
        //        != BLAST_SEQSRC_EOF) {
        //    ...
        //    status = s_BlastSearchEngineCore(...);
        // }
        // ```
        bar.inc(1);
        st.diag_offset = diag_offset;
        emit_subject_hits(st, s_idx);
    };
    // No NCBI counterpart: LOSAT_TIMING, the clock before the subject loop (x_search_timing.rs).
    let x_t_subjects = x_times.now();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:82-88
    // ```c
    // if (num_threads > 1) {
    //     SetNumberOfThreads(num_threads);
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1475
    // ```c
    // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //        != BLAST_SEQSRC_EOF) {
    //    ...
    //    status = s_BlastSearchEngineCore(...);
    // }
    // ```
    // LOSAT_X_TBXNEXT: its loop takes the subjects of a parallel search of several subjects.
    if x_subject_next::on() {
        crate::utils::xreach::x_reach!(
            "LOSAT_X_TBXNEXT",
            use_parallel && !use_serial_scan_chunks && subjects.count() > 1
        );
    }

    #[cfg(all(feature = "parallel", target_arch = "wasm32", feature = "wasm-threads"))]
    if use_parallel && !use_serial_scan_chunks && subjects.count() > 1 {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1427
        // ```c
        // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
        // /* iterate over all subject sequences */
        // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
        //        != BLAST_SEQSRC_EOF) {
        //    if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
        //        continue;
        //    }
        // ```
        // The threads' subjects go to the collector the subject loop keeps (`emit_subject_hits`;
        // no writer thread on this target), which takes them in OID order.
        let parallel_pool = parallel_pool;
        if x_subject_next::on() {
            // EXPERIMENT (LOSAT_X_TBXNEXT): as in native builds below, each thread of the pool
            // takes the next subject when it is free, with one worker state per thread made when
            // it takes its first subject (tblastx/x_subject_next.rs). Same `process_subject`.
            parallel_pool.install(|| {
                x_subject_next::for_each_subject(
                    subjects.count(),
                    &|| WorkerState {
                        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
                        // ```c
                        //     /* iterate over all subject sequences */
                        //     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
                        //            != BLAST_SEQSRC_EOF) {
                        // ```
                        // The state of one thread's iterations of this loop (as below).
                        tx: None,
                        hits: Vec::new(),
                        offset_pairs: vec![OffsetPair::default(); offset_pairs_len],
                        diag_array: vec![DiagStruct::default(); x_worker_diag_cells],
                        // NCBI: diag_table->offset = window_size;
                        // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:63
                        diag_offset: window,
                        sc,
                    },
                    &process_subject,
                );
            });
        } else {
            parallel_pool.install(|| {
                (0..subjects.count()).into_par_iter().for_each_init(
                    || WorkerState {
                        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1475
                        // ```c
                        // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
                        //        != BLAST_SEQSRC_EOF) {
                        //    ...
                        //    status = s_BlastSearchEngineCore(...);
                        // }
                        // ```
                        tx: None,
                        hits: Vec::new(),
                        offset_pairs: vec![OffsetPair::default(); offset_pairs_len],
                        diag_array: vec![DiagStruct::default(); x_worker_diag_cells],
                        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-63
                        // ```c
                        // diag_table->diag_array_length = diag_array_length;
                        // diag_table->diag_mask = diag_array_length-1;
                        // diag_table->offset = window_size;
                        // ```
                        diag_offset: window,
                        sc,
                    },
                    &process_subject,
                );
            });
        }
    } else {
        for_each_subjects(
            subjects.count(),
            || WorkerState {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1475
                // ```c
                // while (...) {
                //    ...
                //    status = s_BlastSearchEngineCore(..., &hsp_list, ...);
                // }
                // ```
                tx: tx_opt.clone(),
                hits: Vec::new(),
                offset_pairs: vec![OffsetPair::default(); offset_pairs_len],
                diag_array: vec![DiagStruct::default(); x_worker_diag_cells],
                // NCBI: diag_table->offset = window_size;
                // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:63
                diag_offset: window,
                sc,
            },
            &process_subject,
        );
    }

    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    if use_parallel && !use_serial_scan_chunks && subjects.count() > 1 && x_subject_next::on() {
        // EXPERIMENT (LOSAT_X_TBXNEXT): each thread of the pool takes the next subject when it is
        // free, with one worker state per thread, made by the same initializer as below when the
        // thread takes its first subject (tblastx/x_subject_next.rs). Same `process_subject`.
        let parallel_pool = parallel_pool;
        parallel_pool.install(|| {
            x_subject_next::for_each_subject(
                subjects.count(),
                &|| WorkerState {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
                    // ```c
                    //     /* iterate over all subject sequences */
                    //     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
                    //            != BLAST_SEQSRC_EOF) {
                    // ```
                    // The state of one thread's iterations of this loop (as below).
                    tx: tx_opt.clone(),
                    hits: Vec::new(),
                    offset_pairs: vec![OffsetPair::default(); offset_pairs_len],
                    diag_array: vec![DiagStruct::default(); x_worker_diag_cells],
                    // NCBI: diag_table->offset = window_size;
                    // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:63
                    diag_offset: window,
                    sc,
                },
                &process_subject,
            );
        });
    } else if use_parallel && !use_serial_scan_chunks && subjects.count() > 1 {
        let parallel_pool = parallel_pool;
        parallel_pool.install(|| {
            (0..subjects.count()).into_par_iter().for_each_init(
                || WorkerState {
                    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1475
                    // ```c
                    // while (...) {
                    //    ...
                    //    status = s_BlastSearchEngineCore(..., &hsp_list, ...);
                    // }
                    // ```
                    tx: tx_opt.clone(),
                    hits: Vec::new(),
                    offset_pairs: vec![OffsetPair::default(); offset_pairs_len],
                    diag_array: vec![DiagStruct::default(); x_worker_diag_cells],
                    // NCBI: diag_table->offset = window_size;
                    // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:63
                    diag_offset: window,
                    sc,
                },
                &process_subject,
            );
        });
    } else {
        for_each_subjects(
            subjects.count(),
            || WorkerState {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1475
                // ```c
                // while (...) {
                //    ...
                //    status = s_BlastSearchEngineCore(..., &hsp_list, ...);
                // }
                // ```
                tx: tx_opt.clone(),
                hits: Vec::new(),
                offset_pairs: vec![OffsetPair::default(); offset_pairs_len],
                diag_array: vec![DiagStruct::default(); x_worker_diag_cells],
                // NCBI: diag_table->offset = window_size;
                // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:63
                diag_offset: window,
                sc,
            },
            &process_subject,
        );
    }

    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    {
        for_each_subjects(
            subjects.count(),
            || WorkerState {
                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1475
                // ```c
                // while (...) {
                //    ...
                //    status = s_BlastSearchEngineCore(..., &hsp_list, ...);
                // }
                // ```
                tx: tx_opt.clone(),
                hits: Vec::new(),
                offset_pairs: vec![OffsetPair::default(); offset_pairs_len],
                diag_array: vec![DiagStruct::default(); x_worker_diag_cells],
                // NCBI: diag_table->offset = window_size;
                // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:63
                diag_offset: window,
                sc,
            },
            &process_subject,
        );
    }

    // No NCBI counterpart: LOSAT_TIMING, the clock after the subject loop (x_search_timing.rs).
    let x_t_collect = x_times.now();
    // Close channel so the writer can exit.
    if let Some(tx) = tx_opt {
        drop(tx);
    }

    bar.finish();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_seqalign.cpp:1569-1577
    // ```c
    // for (int index = 0; index < hit_list->hsplist_count; index++) {
    //     BlastHSPList* hsp_list = hit_list->hsplist_array[index];
    //     if (!hsp_list)
    //         continue;
    //     Blast_HSPListSortByEvalue(hsp_list);
    // }
    // ```
    // Each traversal ends with the batch's hit lists in report order (tblastx/hit_collector.rs):
    // from the writer thread, or from the collector the subject loop fed.
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let collected_hits: Option<Vec<TblastxHsp>> = writer.map(|writer| writer.join().unwrap());
    #[cfg(not(all(feature = "parallel", not(target_arch = "wasm32"))))]
    let collected_hits: Option<Vec<TblastxHsp>> = None;
    let hits = match collected_hits {
        Some(hits) => hits,
        None => local_collector
            .into_inner()
            .unwrap()
            .expect("a subject loop without a writer thread keeps the collector")
            .finish(),
    };
    if diag_enabled {
        print_diagnostics_summary(&diagnostics);
    }

    // EXPERIMENT (LOSAT_X_TBXPIPE): the block's lines, printed here, or kept in the batch's result
    // for the batch loop when the pipe runs this search (`x_seed_end`).
    let mut x_timing: Vec<String> = Vec::new();
    if timing_enabled {
        let t_search = t_search_start.elapsed();
        // Diagnostic (LOSAT_TIMING): the same interval for the common line (utils/xtiming.rs).
        // No NCBI counterpart.
        crate::utils::xtiming::add_duration(crate::utils::xtiming::Phase::Search, t_search);
        let scan_s = scan_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9;
        let scan_n = scan_calls.load(AtomicOrdering::Relaxed);
        let ungapped_s = ungapped_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9;
        let ungapped_n = ungapped_calls.load(AtomicOrdering::Relaxed);
        let reeval_s = reeval_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9;
        let reeval_n = reeval_calls.load(AtomicOrdering::Relaxed);
        let linking_s = linking_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9;
        let linking_n = linking_calls.load(AtomicOrdering::Relaxed);
        let identity_s = identity_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9;
        let identity_n = identity_calls.load(AtomicOrdering::Relaxed);

        x_timing.push(format!(
            "[TIMING] read_queries: {:.3}s",
            t_read_queries.as_secs_f64()
        ));
        x_timing.push(format!(
            "[TIMING] build_lookup: {:.3}s",
            t_build_lookup.as_secs_f64()
        ));
        x_timing.push(format!(
            "[TIMING] read_subjects: {:.3}s",
            t_read_subjects.as_secs_f64()
        ));
        x_timing.push(format!(
            "[TIMING] scan_subject: {:.3}s (calls={})",
            scan_s, scan_n
        ));
        x_timing.push(format!(
            "[TIMING] ungapped_extend: {:.3}s (calls={})",
            ungapped_s, ungapped_n
        ));
        x_timing.push(format!(
            "[TIMING] reevaluate: {:.3}s (calls={})",
            reeval_s, reeval_n
        ));
        x_timing.push(format!(
            "[TIMING] sum_stats_linking: {:.3}s (calls={})",
            linking_s, linking_n
        ));
        x_timing.push(format!(
            "[TIMING] identity_calc: {:.3}s (calls={})",
            identity_s, identity_n
        ));
        x_timing.push(format!(
            "[TIMING] search_total: {:.3}s",
            t_search.as_secs_f64()
        ));
        x_timing.push(format!(
            "[TIMING] total: {:.3}s",
            t_total.elapsed().as_secs_f64()
        ));
        // No NCBI counterpart: the stage lines of x_search_timing.rs, after the lines above.
        if let Some(lines) = x_times.lines(t_search_start, x_t_subjects, x_t_collect, t_search) {
            x_timing.extend(lines);
        }
        // No NCBI counterpart: the linking lines of x_link_timing.rs, after the stage lines.
        if let Some(lines) = x_link_timer.lines() {
            x_timing.extend(lines);
        }
    }
    if x_seed_end.is_none() {
        for line in x_timing.drain(..) {
            eprintln!("{line}");
        }
    }
    Ok(TblastxBatch {
        hits,
        searched: true,
        queries: query_stats,
        x_timing,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_chunk_hit(
        ctx_idx: usize,
        q_start: usize,
        q_end: usize,
        s_start: usize,
        s_end: usize,
        score: i32,
    ) -> UngappedHit {
        UngappedHit {
            q_idx: 0,
            s_idx: 0,
            ctx_idx: u32::try_from(ctx_idx).unwrap(),
            s_f_idx: 0,
            q_frame: 1,
            s_frame: 1,
            q_aa_start: u32::try_from(q_start).unwrap(),
            q_aa_end: u32::try_from(q_end).unwrap(),
            s_aa_start: u32::try_from(s_start).unwrap(),
            s_aa_end: u32::try_from(s_end).unwrap(),
            q_seed_off: u32::try_from(q_start).unwrap(),
            s_seed_off: u32::try_from(s_start).unwrap(),
            q_orig_len: 900,
            s_orig_len: 900,
            raw_score: score,
            e_value: f64::INFINITY,
            num_ident: 0,
            hsp_list_order: 0,
            ordering_method: 0,
            linked_set: false,
            start_of_chain: false,
            link_id: 0,
            chain_next_link_id: None,
            hsp_link_num: 0,
            num: 0,
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:221-310
    // ```c
    // if (backup->offset + MAX_DBSEQ_LEN <
    //     backup->hard_ranges[backup->hm_index].right) {
    //     subject->length = MAX_DBSEQ_LEN;
    //     backup->next = backup->offset + MAX_DBSEQ_LEN - dbseq_chunk_overlap;
    // } else {
    //     subject->length = backup->hard_ranges[backup->hm_index].right
    //                     - backup->offset;
    // }
    // ```
    #[test]
    fn subject_split_state_uses_ncbi_max_len_and_overlap() {
        let mut split = SubjectSplitState::new(MAX_DBSEQ_LEN + 250);

        assert_eq!(
            split.next_chunk(MAX_DBSEQ_LEN, DBSEQ_CHUNK_OVERLAP),
            SubjectChunkStatus::Ok(SubjectChunk {
                offset: 0,
                length: MAX_DBSEQ_LEN,
                overlap: 0
            })
        );
        assert_eq!(
            split.next_chunk(MAX_DBSEQ_LEN, DBSEQ_CHUNK_OVERLAP),
            SubjectChunkStatus::Ok(SubjectChunk {
                offset: MAX_DBSEQ_LEN - DBSEQ_CHUNK_OVERLAP,
                length: DBSEQ_CHUNK_OVERLAP + 250,
                overlap: DBSEQ_CHUNK_OVERLAP
            })
        );
        assert_eq!(
            split.next_chunk(MAX_DBSEQ_LEN, DBSEQ_CHUNK_OVERLAP),
            SubjectChunkStatus::Done
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2857-2995
    // ```c
    // if (hsp1->subject.end > split_offsets[0]) { ... }
    // if (hsp2->subject.offset < split_offsets[0] + chunk_overlap_size) { ... }
    // if (ABS(end_diag - start_diag) < OVERLAP_DIAG_CLOSE) {
    //    if (s_BlastMergeTwoHSPs(hsp1, hsp2, allow_gap)) { ... }
    // }
    // ```
    #[test]
    fn merge_tblastx_subject_chunk_hits_merges_overlap_strip_same_context() {
        let mut combined = vec![make_chunk_hit(0, 0, 120, 0, 120, 100)];
        let incoming = vec![make_chunk_hit(0, 100, 220, 100, 220, 120)];

        merge_tblastx_subject_chunk_hits(&mut combined, incoming, 100, DBSEQ_CHUNK_OVERLAP);

        assert_eq!(combined.len(), 1);
        assert_eq!(combined[0].q_aa_start, 0);
        assert_eq!(combined[0].q_aa_end, 220);
        assert_eq!(combined[0].s_aa_start, 0);
        assert_eq!(combined[0].s_aa_end, 220);
        assert_eq!(combined[0].q_seed_off, 100);
        assert_eq!(combined[0].s_seed_off, 100);
        assert_eq!(combined[0].raw_score, 201);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2968-2970
    // ```c
    // /* Skip already deleted HSPs, or HSPs from different contexts */
    // if (!hsp2 || hsp1->context != hsp2->context)
    //    continue;
    // ```
    #[test]
    fn merge_tblastx_subject_chunk_hits_keeps_different_context_hits() {
        let mut combined = vec![make_chunk_hit(0, 0, 120, 0, 120, 100)];
        let incoming = vec![make_chunk_hit(1, 100, 220, 100, 220, 120)];

        merge_tblastx_subject_chunk_hits(&mut combined, incoming, 100, DBSEQ_CHUNK_OVERLAP);

        assert_eq!(combined.len(), 2);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:221-310
    // ```c
    // if (backup->offset + MAX_DBSEQ_LEN <
    //     backup->hard_ranges[backup->hm_index].right) {
    //     subject->length = MAX_DBSEQ_LEN;
    //     backup->next = backup->offset + MAX_DBSEQ_LEN - dbseq_chunk_overlap;
    // } else {
    //     subject->length = backup->hard_ranges[backup->hm_index].right
    //                     - backup->offset;
    // }
    // ```
    #[test]
    fn subject_split_state_keeps_short_subjects_as_single_chunk() {
        let mut split = SubjectSplitState::new(MAX_DBSEQ_LEN - 1);

        assert_eq!(
            split.next_chunk(MAX_DBSEQ_LEN, DBSEQ_CHUNK_OVERLAP),
            SubjectChunkStatus::Ok(SubjectChunk {
                offset: 0,
                length: MAX_DBSEQ_LEN - 1,
                overlap: 0
            })
        );
        assert_eq!(
            split.next_chunk(MAX_DBSEQ_LEN, DBSEQ_CHUNK_OVERLAP),
            SubjectChunkStatus::Done
        );
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2857-2995
    // ```c
    // if (contexts_per_query < 0) {      /* subject seq is split */
    //    if (hsp1->subject.end > split_offsets[0]) { ... }
    //    if (hsp2->subject.offset < split_offsets[0] + chunk_overlap_size) { ... }
    // }
    // ```
    #[test]
    fn merge_tblastx_subject_chunk_hits_keeps_previous_hits_before_overlap_strip() {
        let mut combined = vec![make_chunk_hit(0, 0, 90, 0, 90, 100)];
        let incoming = vec![make_chunk_hit(0, 100, 180, 100, 180, 120)];

        merge_tblastx_subject_chunk_hits(&mut combined, incoming, 100, DBSEQ_CHUNK_OVERLAP);

        assert_eq!(combined.len(), 2);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2857-2995
    // ```c
    // if (contexts_per_query < 0) {      /* subject seq is split */
    //    if (hsp2->subject.offset < split_offsets[0] + chunk_overlap_size) { ... }
    // }
    // ```
    #[test]
    fn merge_tblastx_subject_chunk_hits_keeps_current_hits_after_overlap_strip() {
        let mut combined = vec![make_chunk_hit(0, 0, 140, 0, 140, 100)];
        let incoming = vec![make_chunk_hit(0, 200, 260, 200, 260, 120)];

        merge_tblastx_subject_chunk_hits(&mut combined, incoming, 100, DBSEQ_CHUNK_OVERLAP);

        assert_eq!(combined.len(), 2);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1488-1533
    // ```c
    // if(hsp1->subject.frame != hsp2->subject.frame) return FALSE;
    // if (CONTAINED_IN_HSP(...) || CONTAINED_IN_HSP(...)) {
    //    ...
    //    return TRUE;
    // }
    // ```
    #[test]
    fn merge_tblastx_subject_chunk_hits_merges_hsp_fully_inside_overlap_strip() {
        let mut combined = vec![make_chunk_hit(0, 0, 190, 0, 190, 190)];
        let incoming = vec![make_chunk_hit(0, 120, 180, 120, 180, 60)];

        merge_tblastx_subject_chunk_hits(&mut combined, incoming, 100, DBSEQ_CHUNK_OVERLAP);

        assert_eq!(combined.len(), 1);
        assert_eq!(combined[0].q_aa_start, 0);
        assert_eq!(combined[0].q_aa_end, 190);
        assert_eq!(combined[0].s_aa_start, 0);
        assert_eq!(combined[0].s_aa_end, 190);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1488-1533
    // ```c
    // if(hsp1->subject.frame != hsp2->subject.frame) return FALSE;
    // ```
    #[test]
    fn merge_tblastx_subject_chunk_hits_merges_negative_subject_frame_overlap() {
        let mut combined = vec![make_chunk_hit(0, 0, 120, 0, 120, 100)];
        combined[0].s_frame = -1;
        let mut incoming_hit = make_chunk_hit(0, 100, 220, 100, 220, 120);
        incoming_hit.s_frame = -1;

        merge_tblastx_subject_chunk_hits(
            &mut combined,
            vec![incoming_hit],
            100,
            DBSEQ_CHUNK_OVERLAP,
        );

        assert_eq!(combined.len(), 1);
        assert_eq!(combined[0].s_frame, -1);
        assert_eq!(combined[0].q_aa_end, 220);
        assert_eq!(combined[0].s_aa_end, 220);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3038-3051
    // ```c
    // hsp->subject.offset += offset;
    // hsp->subject.end += offset;
    // hsp->subject.gapped_start += offset;
    // ```
    #[test]
    fn adjust_tblastx_chunk_subject_offsets_adjusts_subject_coordinates() {
        let mut hits = vec![make_chunk_hit(0, 10, 40, 20, 50, 100)];

        adjust_tblastx_chunk_subject_offsets(&mut hits, 500);

        assert_eq!(hits[0].s_aa_start, 520);
        assert_eq!(hits[0].s_aa_end, 550);
        assert_eq!(hits[0].s_seed_off, 520);
        assert_eq!(hits[0].q_aa_start, 10);
        assert_eq!(hits[0].q_seed_off, 10);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:83-127
    // ```c
    // for (s = s_first; s <= s_last; s++) { ... }
    // s_range[1] = (Int4)(s - subject->sequence);
    // ```
    #[test]
    fn scan_interior_ranges_keep_right_lookahead_without_duplicate_emission() {
        let wordsize = 3usize;
        let subject_len = 8usize;
        let base_ranges = [(0, subject_len as i32)];

        let left =
            clip_tblastx_seq_ranges_for_scan_interior(&base_ranges, 0, 4, wordsize, subject_len);
        let right =
            clip_tblastx_seq_ranges_for_scan_interior(&base_ranges, 4, 8, wordsize, subject_len);

        assert_eq!(left, vec![(0, 6)]);
        assert_eq!(right, vec![(4, 8)]);

        let left_emitted: Vec<i32> = (left[0].0..=left[0].1 - wordsize as i32).collect();
        let right_emitted: Vec<i32> = (right[0].0..=right[0].1 - wordsize as i32).collect();
        assert_eq!(left_emitted, vec![0, 1, 2, 3]);
        assert_eq!(right_emitted, vec![4, 5]);
    }
}
