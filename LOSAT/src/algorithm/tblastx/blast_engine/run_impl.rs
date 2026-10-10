//! Main run() function for backbone mode
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c

use super::*;
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
/// A zeroed diagonal table (NCBI calloc), with huge pages advised (LOSAT_X_THP).
fn x_new_diag_array(size: usize) -> Vec<DiagStruct> {
    let mut v = vec![DiagStruct::default(); size];
    crate::utils::x_hugepage::advise(&mut v);
    v
}

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
        hsp1.q_aa_start,
        hsp1.q_aa_end,
        hsp2.q_aa_start,
        hsp1.s_aa_start,
        hsp1.s_aa_end,
        hsp2.s_aa_start,
    ) || contained_in_hsp(
        hsp1.q_aa_start,
        hsp1.q_aa_end,
        hsp2.q_aa_end,
        hsp1.s_aa_start,
        hsp1.s_aa_end,
        hsp2.s_aa_end,
    ) {
        let len1 = hsp1.q_aa_end.saturating_sub(hsp1.q_aa_start) as f64;
        let len2 = hsp2.q_aa_end.saturating_sub(hsp2.q_aa_start) as f64;
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

        let new_len = hsp1.q_aa_end.saturating_sub(hsp1.q_aa_start) as f64;
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
        if hsp.s_aa_end > split_offset {
            combined_overlap.push(idx);
        }
    }

    let mut incoming_overlap = Vec::new();
    for (idx, hsp) in incoming.iter().enumerate() {
        if hsp.s_aa_start < split_offset.saturating_add(chunk_overlap_size) {
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
        hit.s_aa_start += offset;
        hit.s_aa_end += offset;
        hit.s_seed_off += offset;
    }
}

/// The bases of a subject as the preliminary search reads them: NCBI's compressed (ncbi2na)
/// subject, in which each ambiguous base is a compatible base drawn from `CRandom` seeded
/// with the subject's length (`resolve_ncbi4na_to_ncbi2na`). `None` for a subject of
/// A, C, G and T only, whose compressed bases are its own.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520
/// ```c
///     virtual SBlastSequence GetCompressedPlusStrand() {
///         SBlastSequence retval(size());
///         string ncbi4na = kEmptyStr;
///         m_SeqVector.GetSeqData(m_SeqVector.begin(), m_SeqVector.end(), ncbi4na);
///         s_Ncbi4naToNcbi2na(ncbi4na, size(), retval.data.get());
///         return retval;
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:94-103
/// ```c
/// const Uint1 IUPACNA_TO_NCBI4NA[128]={
///  0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0,
///  0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0,
/// ```
/// The ncbi4na code of a letter is that of `GeneticCode::get` (both cases, `U` as `T`);
/// another byte is a gap (0).
fn preliminary_subject_bases(subject: &[u8]) -> Option<Vec<u8>> {
    if subject
        .iter()
        .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
    {
        return None;
    }
    let ncbi4na: Vec<u8> = subject
        .iter()
        .map(|&base| match base.to_ascii_uppercase() {
            b'U' => 8,
            upper => crate::core::blast_encoding::iupacna_to_ncbi4na(upper).unwrap_or(0),
        })
        .collect();
    Some(
        crate::core::blast_encoding::resolve_ncbi4na_to_ncbi2na(&ncbi4na)
            .into_iter()
            .map(|code| b"ACGT"[code as usize])
            .collect(),
    )
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
    args: TblastxArgs,
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
        "TBLASTX",
        false,
        args.out.as_deref(),
        false,
    )?;
    crate::utils::threading::validate_threads(args.num_threads)?;
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
    use crate::blastinput::input_files;
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
        ReaderConfig::subject("TBLASTX", false, settings.data_loaders),
    );
    let read_subject_records = {
        let mut stderr = std::io::stderr();
        read_subjects(
            &mut subject_source,
            subject_range.as_ref(),
            &mut |message: &[u8]| std::io::Write::write_all(&mut stderr, message),
        )?
    };
    drop(subject_source);
    let (subjects, subject_placements) =
        match cut_subjects(&read_subject_records, subject_range.as_ref())? {
            Some((cut, placements)) => (cut, placements),
            None => (read_subject_records, Placements::default()),
        };
    check_subjects_not_empty(&subjects)?;
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
    let query_file = input_files::open_input(&args.query, "query", "TBLASTX")?;
    if let Some(path) = args.out.as_deref() {
        input_files::check_utf8_file_name(path, "out", "TBLASTX")?;
    }
    let out_file = match args.out.as_deref().filter(|path| path.as_os_str() != "-") {
        Some(path) => Some(std::io::BufWriter::new(
            std::fs::File::create(path).map_err(|_| crate::cli::inaccessible("out", path))?,
        )),
        None => None,
    };
    let pairwise = format == Some(crate::blastinput::app::ReportFormat::Pairwise);
    // LOSAT's reports read the format number that NCBI reads.
    let outfmt = choice.normalized();
    let mut stderr = std::io::stderr();
    let mut stream = crate::cli::ReportStream {
        inner: match out_file {
            Some(file) => Box::new(file) as Box<dyn std::io::Write + Send>,
            None => crate::cli::report_standard_output(),
        },
        failed: false,
    };
    let result = {
        let mut outputs =
            ReportOutputs::single(&outfmt, OutputSink::Writer(&mut stream), &mut stderr);
        search_cli(
            args,
            query_file,
            &subjects,
            &subject_placements,
            &mut outputs,
            settings,
        )
    };
    // What was written before an error (such as the outfmt 0 prolog) stays in the file.
    let flushed = std::io::Write::flush(&mut stream);
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
    if stream.failed && pairwise {
        return Err(crate::cli::NativeError {
            exit: 6,
            message: "BLAST failed to write output\n".to_string(),
        }
        .into());
    }
    result?;
    flushed.context("failed to write the output")
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
    diagnostics: &mut dyn std::io::Write,
) -> Result<Option<SequenceRange>> {
    use crate::blastinput::app::{
        formatting_handler_check, options_error, parse_formatting_string,
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
    formatting_handler_check(&parse_formatting_string(&args.outfmt)?, false)?;
    if args
        .max_target_seqs
        .is_some_and(|max_target_seqs| max_target_seqs < 5)
    {
        diagnostics.write_all(&crate::report::query_warnings::few_matches_warning(
            "tblastx",
        ))?;
    }
    // (`threshold <= 0` is false for NaN, as in C.)
    if args.threshold <= 0.0 {
        return Err(options_error("Non-zero threshold required"));
    }
    if args.word_size > 4 {
        return Err(options_error(
            "Word-size must be less than 6 for protein comparison",
        ));
    }
    if args.evalue <= 0.0 {
        return Err(options_error(
            "expect value or cutoff score must be greater than zero",
        ));
    }
    Ok(query_range)
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
    if args.window_size == 0 {
        anyhow::bail!(
            "-window_size 0 (the one-hit word finder) is not supported by LOSAT's TBLASTX"
        );
    }
    // NCBI reference: c++/src/algo/blast/core/blast_aalookup.c:237
    // ```c
    //     lookup->word_length = opt->word_size;
    // ```
    // LOSAT's TBLASTX builds and scans 3-residue words only.
    if args.word_size != 3 {
        anyhow::bail!(
            "-word_size {} is not supported by LOSAT's TBLASTX (it implements word size 3)",
            args.word_size
        );
    }
    Ok(())
}

/// The checks of the options alone, without the inputs (the `validate` of web ABI v2):
/// NCBI's processing of the options and LOSAT's limits.
pub fn check_options(args: &TblastxArgs) -> Result<()> {
    check_ncbi_options(args, &mut std::io::sink())?;
    check_losat_limits(args)?;
    // The host validates the argv it runs, with its `-num_threads` (docs/web/abi_v2.md).
    crate::utils::threading::validate_threads(args.num_threads)
}

/// The part of `run` after the output is opened: the options, `Query is Empty!`, the
/// reading of the queries and the search.
fn search_cli(
    args: TblastxArgs,
    query_file: std::fs::File,
    subjects: &[FastaRecord],
    subject_placements: &Placements,
    outputs: &mut ReportOutputs<'_>,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    let query_range = check_ncbi_options(&args, outputs.diagnostics)?;
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
        ReaderConfig::query("TBLASTX", false, settings.data_loaders),
    );
    if !shared_standard_input && query_source.stream_is_empty() {
        outputs
            .diagnostics
            .write_all(b"Warning: [tblastx] Query is Empty!\n")?;
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
    let QueryRecords { records, end } = read_queries(&mut query_source);
    drop(query_source);
    let reading = match end {
        QueryEnd::Input => QueryReading::End,
        QueryEnd::EofError => QueryReading::BlankLines,
        QueryEnd::Error { error, warnings } => QueryReading::Error { error, warnings },
    };
    search(
        args,
        &records,
        reading,
        query_range,
        subjects,
        subject_placements,
        outputs,
        None,
    )
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
    args: TblastxArgs,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
    v1: Option<&dyn V1Checks>,
) -> Result<()> {
    // NCBI parses -outfmt before its option handlers read the subjects (BLASTN's
    // `run_local`); `search` checks the formats again.
    check_report_formats(outputs)?;
    // The subject range and its record checks come where NCBI reads the subjects (`run`):
    // the messages of each record read, up to a record whose range starts past its end.
    let subject_range =
        parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject, "TBLASTX")?;
    for record in &subject_records[..subjects_read(subject_records, subject_range.as_ref())] {
        outputs.diagnostics.write_all(&record.warnings)?;
    }
    let ranged_subjects = cut_subjects(subject_records, subject_range.as_ref())?;
    check_subjects_not_empty(subject_records)?;
    let query_range = check_ncbi_options(&args, outputs.diagnostics)?;
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
        searched_subjects,
        subject_placements,
        outputs,
        v1,
    )
}

/// Each requested format as NCBI reads `-outfmt`: NCBI's errors, and a rejection of a
/// format that LOSAT's TBLASTX does not write.
fn check_report_formats(outputs: &ReportOutputs<'_>) -> Result<()> {
    use crate::blastinput::app;
    for format in &outputs.formats {
        let choice = app::parse_formatting_string(format.outfmt)?;
        if app::report_format(&choice, "TBLASTX", false, None, false)?.is_none() {
            app::formatting_handler_check(&choice, false)?;
            app::xinclude_check(&choice)?;
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
    subject_records: &[FastaRecord],
    subject_placements: &Placements,
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
    crate::blastinput::input_files::write_empty_subject_warnings(
        subject_records,
        "tblastx",
        outputs.diagnostics,
    )?;
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
    check_report_formats(outputs)?;
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
            subject_records,
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
    /// Whether ABI v1's rejection of outfmt 0 titles applies (`V1Checks`).
    pub v1_titles: bool,
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
    subject_records: &[FastaRecord],
    before_search: bool,
) -> Result<()> {
    for format in outputs.formats.iter_mut() {
        if report::output_format(format.outfmt) == report::TblastxOutputFormat::Pairwise
            && prolog_before_search(&format.sink) == before_search
        {
            let mut writer = format.sink.open()?;
            crate::report::pairwise::write_tblastx_pairwise_prolog(
                &mut writer,
                report::NCBI_TBLASTX_VERSION,
                &tblastx_database_name(args),
                subject_records.len(),
                subject_records.iter().map(|r| r.seq().len()).sum(),
            )?;
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
    hits: Vec<TblastxHsp>,
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
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2960-2962,2978
    // ```c
    //     	if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
    //     		hitlist_size = args[kArgMaxTargetSequences].AsInteger();
    //     	}
    //     ...
    //     opt.SetHitlistSize(hitlist_size);
    // ```
    // The hit list size is the -max_target_seqs value, or 500 when it is omitted.
    // NCBI reference: c++/src/algo/blast/api/blast_options_local_priv.hpp:1320-1324
    // ```c
    // CBlastOptionsLocal::SetCullingLimit(int s)
    // {
    //     if (s <= 0) {
    //         return;
    //     }
    // ```
    // A culling limit above 0 installs NCBI's culling writer and pipe
    // (`report::culled_hit_order`).
    let hitlist_size = args.max_target_seqs.unwrap_or(500);
    let mut hits = if args.culling_limit > 0 {
        let query_lengths: Vec<usize> = run.queries.iter().map(|query| query.seq().len()).collect();
        report::culled_hit_order(hits, hitlist_size, args.culling_limit, &query_lengths)
    } else {
        report::final_hit_order(hits, hitlist_size)
    };
    if run.v1_titles && output_formats.contains(&report::TblastxOutputFormat::Pairwise) {
        check_shown_subject_titles(&hits, run.subjects)?;
    }
    report::set_displayed_identities(&mut hits, run.queries, run.subjects, &query_code, &db_code);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1411
    // ```c
    // CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
    // ```
    // The caller receives the same final result that every formatter prints, with the
    // displayed rows of outfmt 0.
    let pairwise = (outputs.hits.is_some()
        || output_formats.contains(&report::TblastxOutputFormat::Pairwise))
    .then(|| {
        report::pairwise_hits(
            &hits,
            run.queries,
            run.subjects,
            run.query_stats,
            &query_code,
            &db_code,
            &run.ranges.input.placements,
            run.ranges.subject_placements,
        )
    });
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
    let subject_ids = if output_formats
        .iter()
        .any(|&format| format != report::TblastxOutputFormat::Pairwise)
    {
        tabular_subject_ids(&hits, run.subjects)?
    } else {
        Vec::new()
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:497-498
    // ```c++
    //     //get defline
    //     sdl->defline = CDeflineGenerator().GenerateDefline(m_ScopeRef->GetBioseqHandle(*(sdl->id)), sequence::CDeflineGenerator::fLeavePrefixSuffix);
    // ```
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
    let database = tblastx_database_name(args);
    // The `Query=` and `# Query:` text of each query: its title bytes.
    let query_titles: Vec<&[u8]> = run.queries.iter().map(|query| query.title()).collect();
    let unsearched: Vec<bool> = run.searched.iter().map(|&searched| !searched).collect();
    let mut format_warnings = Some(crate::report::query_warnings::QueryWarnings {
        before: run.warnings,
        sink: &mut *outputs.diagnostics,
    });
    let observer = &mut outputs.observer;
    for (format_index, (format, &output_format)) in
        outputs.formats.iter_mut().zip(&output_formats).enumerate()
    {
        let mut probe = observer
            .as_deref_mut()
            .map(|observer| FormatProbe::new(observer, format_index));
        // The warnings are written once, with the first format (`QueryWarnings`).
        let mut warnings = format_warnings.take();
        // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:68-93
        // ```c
        // CBlastFormat::CBlastFormat(..., CNcbiOstream& outfile, ...)
        //     : m_FormatType(format_type), ..., m_Outfile(outfile),
        // ```
        let prolog = !prolog_before_search(&format.sink);
        let mut writer = format.sink.open()?;
        match output_format {
            report::TblastxOutputFormat::Pairwise => {
                let pairwise = pairwise
                    .as_deref()
                    .expect("pairwise hits are built for outfmt 0");
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
                let (num_descriptions, num_alignments) = match args.max_target_seqs {
                    Some(max_target_seqs) => (max_target_seqs, max_target_seqs),
                    None => (500, 250),
                };
                let pairwise_report = crate::report::pairwise::TblastxPairwiseReport {
                    version: report::NCBI_TBLASTX_VERSION.to_string(),
                    database_name: database.clone(),
                    database_num_sequences: run.subjects.len(),
                    database_total_letters: run.subjects.iter().map(|r| r.seq().len()).sum(),
                    word_threshold: args.threshold,
                    window_size: args.window_size as usize,
                    num_descriptions,
                    num_alignments,
                    unsearched: unsearched.clone(),
                    // A search that stops at a query batch (`Empty CBlastQueryVector` with
                    // `-query_loc`) writes no epilog (tblastx_app.cpp:209, `PrintEpilog`).
                    epilog: run.epilog,
                    prolog,
                };
                crate::report::pairwise::write_tblastx_pairwise_report(
                    pairwise,
                    &mut writer,
                    &crate::report::pairwise::PairwiseConfig {
                        program: "tblastx".to_string(),
                        ..crate::report::pairwise::PairwiseConfig::default()
                    },
                    &queries,
                    &subject_titles,
                    &pairwise_report,
                    probe.as_mut(),
                    warnings.as_mut(),
                )?;
            }
            report::TblastxOutputFormat::Tabular
            | report::TblastxOutputFormat::TabularWithComments => {
                report::write_tabular(
                    &hits,
                    &mut writer,
                    output_format == report::TblastxOutputFormat::TabularWithComments,
                    &query_titles,
                    &database,
                    run.query_ids,
                    &subject_ids,
                    &unsearched,
                    run.epilog,
                    probe.as_mut(),
                    warnings.as_mut(),
                )?;
            }
        }
        writer.flush()?;
    }
    // Without an output format the warnings are written all at once.
    if let Some(warnings) = format_warnings {
        for query_warnings in warnings.before {
            warnings.sink.write_all(query_warnings)?;
        }
    }
    Ok(())
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
    let ids: Vec<Option<Arc<[u8]>>> = subjects
        .iter()
        .map(|subject| {
            crate::report::defline::tabular_subject_id(
                subject.local_id.as_bytes(),
                &subject.title,
                false,
            )
            .ok()
            .map(Arc::from)
        })
        .collect();
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
// of all of them once (`write_tblastx_outputs`), so that what stops the run at a batch
// comes after the reports of the batches before (without the epilog).
#[allow(clippy::too_many_arguments)]
fn run_in_pool(
    args: TblastxArgs,
    query_runs: &QueryRuns<'_>,
    reading: QueryReading,
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
    batch_size: u32,
    v1_titles: bool,
) -> Result<()> {
    check_report_formats(outputs)?;
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
        subject_records,
        v1_titles,
    };
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:173
    // ```c
    //         formatter.PrintProlog();
    // ```
    write_pairwise_prologs(outputs, &args, subject_records, true)?;
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:176-178
    // ```c
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
    //
    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
    // ```
    // A reader at the end of its input reads no batch (`!input.End()`): the report of no
    // query, the prolog and the epilog (an empty pipe, or `-query -` after `-subject -`).
    if query_runs.input_records.is_empty() && matches!(reading, QueryReading::End) {
        return reports.write(outputs, SearchedBatches::default(), &[], true);
    }
    let mut done = SearchedBatches::default();
    // The warnings written before each query's report (`QueryWarnings`).
    let mut warnings: Vec<Vec<u8>> = vec![Vec::new(); query_records.len()];
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
        warnings[start].extend_from_slice(&batch_warnings);
        let batch = match search_query_batch(&args, batch_queries, subject_records, parallel_pool) {
            Ok(batch) => batch,
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
        for (offset, query) in batch_queries.iter().enumerate() {
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
        // The batch numbers its queries from 0.
        done.hits.extend(batch.hits.into_iter().map(|mut hsp| {
            hsp.hit.q_idx += start as u32;
            hsp
        }));
        done.query_stats.extend(batch.queries);
        done.searched
            .extend(std::iter::repeat_n(batch.searched, end - start));
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
}

/// What the reports of the batches searched so far need besides their results.
struct BatchReports<'a> {
    args: &'a TblastxArgs,
    query_runs: &'a QueryRuns<'a>,
    query_ids: &'a [Arc<[u8]>],
    subject_records: &'a [FastaRecord],
    v1_titles: bool,
}

impl BatchReports<'_> {
    /// Writes the reports of the queries of the batches searched (`done`), with NCBI's
    /// epilog (`PrintEpilog`) when `epilog`. `warnings` are by query (at least those of
    /// `done`).
    fn write(
        &self,
        outputs: &mut ReportOutputs<'_>,
        done: SearchedBatches,
        warnings: &[Vec<u8>],
        epilog: bool,
    ) -> Result<()> {
        let reported = done.searched.len();
        write_tblastx_outputs(
            done.hits,
            outputs,
            &TblastxReportRun {
                args: self.args,
                queries: &self.query_runs.records[..reported],
                subjects: self.subject_records,
                query_ids: &self.query_ids[..reported],
                query_stats: &done.query_stats,
                searched: &done.searched,
                warnings: &warnings[..reported],
                ranges: self.query_runs,
                epilog,
                v1_titles: self.v1_titles,
            },
        )
    }

    /// Writes what NCBI has written when a query batch fails, before its error: the reports
    /// of the batches before (`done`; without the epilog), or, for the first batch, the
    /// outfmt 0 prolog of the sinks that write it with the reports (the others wrote it
    /// before the search, `write_pairwise_prologs`).
    fn write_before_failed_batch(
        &self,
        outputs: &mut ReportOutputs<'_>,
        done: SearchedBatches,
        warnings: &[Vec<u8>],
    ) -> Result<()> {
        if done.searched.is_empty() {
            write_pairwise_prologs(outputs, self.args, self.subject_records, false)
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

fn search_query_batch(
    args: &TblastxArgs,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
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
    let extension_debug_enabled = std::env::var("LOSAT_DEBUG_EXTENSION").is_ok();
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
    let window = args.window_size as i32;
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
    let ungapped_params_for_xdrop = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);

    let diag_enabled = diagnostics_enabled();
    let debug_cutoffs_all = std::env::var_os("LOSAT_DEBUG_CUTOFFS_ALL").is_some();
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
    let debug_hsp_saving = std::env::var_os("LOSAT_DEBUG_HSP_SAVING").is_some();
    let diagnostics = std::sync::Arc::new(DiagnosticCounters::default());

    // Debug: optional scan output dump around a specific subject offset.
    // The offset_pairs correspond to s_BlastAaScanSubject copy loop.
    // Reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:99-115
    let scan_debug_center = std::env::var("LOSAT_DEBUG_SCAN_SOFF")
        .ok()
        .and_then(|v| v.parse::<i32>().ok());
    let scan_debug_window = std::env::var("LOSAT_DEBUG_SCAN_WINDOW")
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
    let mut query_frames: Vec<Vec<QueryFrame>> = queries_raw
        .iter()
        .map(|r| generate_frames(r.seq(), &query_code))
        .collect();

    // NCBI blast_args.cpp:404-406: opt.SetSegFilteringWindow(...);
    // opt.SetSegFilteringLocut(...); opt.SetSegFilteringHicut(...);
    if let Some(params) = args.seg_spec()?.params() {
        let seg = SegMasker::new(params.window, params.locut, params.hicut);
        for frames in &mut query_frames {
            for frame in frames {
                if frame.aa_seq.len() >= 3 {
                    for m in seg.mask_sequence(&frame.aa_seq[1..frame.aa_seq.len() - 1]) {
                        frame.seg_masks.push((m.start, m.end));
                    }

                    // NCBI BLAST query masking semantics (SEG, etc.): keep an unmasked copy
                    // (`sequence_nomask`), then overwrite masked residues in the working
                    // query sequence buffer with X.
                    //
                    // NCBI reference (verbatim):
                    //   const Uint1 kProtMask = 21;     /* X in NCBISTDAA */
                    //   query_blk->sequence_start_nomask = BlastMemDup(query_blk->sequence_start, total_length);
                    //   query_blk->sequence_nomask = query_blk->sequence_start_nomask + 1;
                    //   buffer[index] = kMaskingLetter;
                    //
                    // LOSAT uses NCBISTDAA encoding where X = 21.
                    if !frame.seg_masks.is_empty() {
                        if frame.aa_seq_nomask.is_none() {
                            frame.aa_seq_nomask = Some(frame.aa_seq.clone());
                        }
                        const X_MASK_NCBISTDAA: u8 = 21; // NCBI: kProtMask = 21
                        let raw_end_exclusive = frame.aa_seq.len().saturating_sub(1); // keep last sentinel untouched
                        for &(s, e) in &frame.seg_masks {
                            let raw_s = 1usize.saturating_add(s);
                            let raw_e = 1usize.saturating_add(e).min(raw_end_exclusive);
                            for pos in raw_s..raw_e {
                                frame.aa_seq[pos] = X_MASK_NCBISTDAA;
                            }
                        }
                    }
                }
            }
        }
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
    let (lookup, contexts) = build_ncbi_lookup(
        &query_frames,
        crate::core::blast_util::ncbi_int4_from_double(args.threshold),
        &ungapped_params_for_xdrop, // Used for x_dropoff calculation only
        true,
    );
    t_build_lookup = t_phase_build_lookup.elapsed();

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
    let subjects_raw: &[FastaRecord] = subject_records;
    if queries_raw.is_empty() || subjects_raw.is_empty() {
        return Ok(TblastxBatch {
            hits: Vec::new(),
            searched: true,
            queries: vec![TblastxQueryStats::default(); queries_raw.len()],
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
                && subjects_raw.iter().any(|r| {
                    tblastx_estimated_subject_chunk_work_items(
                        r.seq().len(),
                        tblastx_max_dbseq_len_for_run(),
                    ) > 1
                })));
    crate::utils::threading::report_stage(
        "tblastx",
        "subjects",
        subjects_raw.len(),
        use_parallel && !use_serial_scan_chunks && subjects_raw.len() > 1,
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
    // BLOSUM62 ungapped: lambda=0.3176, K=0.134 (used for tblastx)
    // BLOSUM62 gapped:   lambda=0.267,  K=0.041 (NOT used for tblastx)
    let ungapped_params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    // Note: gapped_params is kept for API compatibility but NOT used for tblastx
    let gapped_params = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    // Use UNGAPPED params for all calculations (NCBI parity for tblastx)
    let params = ungapped_params.clone();

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
    let use_channel = use_parallel && !use_serial_scan_chunks && subjects_raw.len() > 1;
    #[cfg(any(not(feature = "parallel"), target_arch = "wasm32"))]
    let use_channel = false;

    let (tx_opt, mut rx_opt) = if use_channel {
        let (tx, rx) = channel::<Vec<TblastxHsp>>();
        (Some(tx), Some(rx))
    } else {
        (None, None)
    };
    let evalue_threshold = args.evalue;

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

    let lookup_ref = &lookup;
    let contexts_ref = &contexts;
    let _gapped_params_ref = &gapped_params; // Unused - tblastx uses ungapped params

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
    fn for_each_subjects<FInit, FBody>(
        subjects: &[FastaRecord],
        init: FInit,
        mut body: FBody,
    ) -> WorkerState
    where
        FInit: FnOnce() -> WorkerState,
        FBody: FnMut(&mut WorkerState, (usize, &FastaRecord)),
    {
        let mut state = init();
        for (s_idx, s_rec) in subjects.iter().enumerate() {
            body(&mut state, (s_idx, s_rec));
        }
        state
    }

    // NCBI reference: c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-180,261-264;
    // c++/src/algo/blast/core/blast_engine.c:1372-1374,1434-1443
    // m_iTotalLength += (Int8) (*iter)->length;
    // avg_length = (Uint4) (total_length / num_seqs);
    // BlastInitialWordParametersUpdate(..., BlastSeqSrcGetAvgSeqLen(seq_src), ...);
    // if (db_length == 0) { BLAST_OneSubjectUpdateParameters(...); }
    let db_length_nucl = subjects_raw
        .iter()
        .map(|r| r.seq().len() as i64)
        .sum::<i64>();
    let db_num_seqs = subjects_raw.len() as i64;
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
        });
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
            subjects_raw.first().map_or(0, |record| record.seq().len()) as i64,
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
    // The word cutoff uses the average subject nucleotide length.
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
        let eff_lengths = compute_eff_lengths_tblastx(
            query_len_aa,
            db_length_nucl,
            db_num_seqs,
            ctx_params, // tblastx uses per-context ungapped params (kbp_gap is NULL)
        );
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
            seg_masks: frames
                .iter()
                .map(|frame| (frame.frame, frame.seg_masks.clone()))
                .collect(),
        })
        .collect();
    for (ctx_idx, ctx) in contexts_ref.iter().enumerate() {
        let stats = &mut query_stats[ctx.q_idx as usize];
        if ctx.is_valid && stats.karlin.is_none() {
            stats.karlin = Some(ctx.karlin_params);
            stats.eff_searchsp = eff_searchsp_per_context[ctx_idx];
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
        // The collector returns the final hits; they are written after the join, at
        // the one place that writes every output.
        let rx = rx_opt.take().expect("rx must be available for writer");
        Some(std::thread::spawn(move || -> Vec<TblastxHsp> {
            let mut all: Vec<TblastxHsp> = Vec::new();
            while let Ok(h) = rx.recv() {
                all.extend(h);
            }
            // NCBI reference: c++/src/algo/blast/core/blast_hits.c:1988-1996
            // ```c
            //    cutoff = hit_options->expect_value;
            // ...
            //       if (hsp->evalue > cutoff) {
            // ```
            // An HSP is removed when its e-value is greater than the cutoff, so a NaN cutoff
            // (`-evalue -nan`) keeps every HSP.
            all.retain(|h| !(h.hit.e_value > evalue_threshold));
            all
        }))
    } else {
        None
    };

    let process_subject = |st: &mut WorkerState, (s_idx, s_rec): (usize, &FastaRecord)| {
        // NCBI reference: c++/src/algo/blast/core/blast_engine.c:1318-1330,1429-1431
        // return word_length * 3 + 2;
        // if (subject->length < min_subj_seq_length) { ... continue; }
        // Statistics include all records, including ones too short to search.
        if s_rec.seq().len() < (wordsize * 3 + 2) as usize {
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
        let s_frames = generate_frames(s_rec.seq(), &db_code);
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
        let s_frames_preliminary_owned =
            preliminary_subject_bases(s_rec.seq()).map(|bases| generate_frames(&bases, &db_code));
        let s_frames_preliminary: &[QueryFrame] =
            s_frames_preliminary_owned.as_deref().unwrap_or(&s_frames);

        let s_len = s_rec.seq().len();

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

            let scan_chunk_plain = |chunk: SubjectChunk,
                                    offset_pairs: &mut [OffsetPair],
                                    diag_array: &mut [DiagStruct],
                                    diag_offset: &mut i32,
                                    scan_chunk_size: Option<usize>|
             -> TblastxChunkScanResult {
                let mut init_hsps: Vec<InitHSP> = Vec::new();
                let mut stats = TblastxChunkScanStats::default();
                let chunk_end = chunk.offset.saturating_add(chunk.length);
                let subject = &subject_all[chunk.offset..chunk_end];

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
                        let hits = s_blast_aa_scan_subject(
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
                                ) = extend_hit_two_hit(
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
                                );
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
                    get_ungapped_hsp_list(init_hsps, contexts_ref, &s_frames)
                };
                adjust_tblastx_chunk_subject_offsets(&mut hits, chunk.offset);
                TblastxChunkScanResult { chunk, hits, stats }
            };

            // EXPERIMENT (LOSAT_X_SEEDBUCKET): `scan_chunk_plain` above is the reference
            // loop, untouched.  This variant appends the hits of every scan call to
            // per-diagonal-range buckets and runs the same per-hit body (the macro
            // below, a copy of the loop body of `scan_chunk_plain`) bucket by bucket,
            // then puts the saved HSPs back into scan order.  See `x_seed_bucket`.
            let scan_chunk_bucketed = |chunk: SubjectChunk,
                                       offset_pairs: &mut [OffsetPair],
                                       diag_array: &mut [DiagStruct],
                                       diag_offset: &mut i32,
                                       scan_chunk_size: Option<usize>|
             -> TblastxChunkScanResult {
                let mut init_hsps: Vec<InitHSP> = Vec::new();
                let mut stats = TblastxChunkScanStats::default();
                let chunk_end = chunk.offset.saturating_add(chunk.length);
                let subject = &subject_all[chunk.offset..chunk_end];

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
                let mut x_seq_counter: u32 = 0;
                let mut seq_keys: Vec<u32> = Vec::new();
                let diag_offset_value: i32 = *diag_offset;
                // The NCBI two-hit loop body for one hit (aa_ungapped.c:531-606): a copy
                // of the loop body of `scan_chunk_plain`, as a macro expanded in place in
                // the bucket flush callbacks.  Local names resolve at this definition site
                // (macro hygiene), so the body reads `init_hsps`, `seq_keys`, `stats` and
                // `diag_array` of this closure directly.  `break 'hit` is the `continue`
                // of the NCBI loop; `seq` is the hit's position in the scan stream.
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
                                ) = extend_hit_two_hit(
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
                                );
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
                        let hits = s_blast_aa_scan_subject(
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
                        let offset_pairs_ptr = offset_pairs.as_ptr();
                        for i in 0..hits as usize {
                            // SAFETY: i < hits, and hits <= offset_array_size (checked by scan)
                            let pair = unsafe { &*offset_pairs_ptr.add(i) };
                            buckets.push(pair.q_off, pair.s_off, x_seq_counter + i as u32);
                        }
                        x_seq_counter += hits as u32;
                        if buckets.flush_due(x_bucket_budget) {
                            buckets.flush(|q, s, seq| x_two_hit_body!(q, s, seq));
                        }
                    }
                }
                buckets.flush(|q, s, seq| x_two_hit_body!(q, s, seq));
                // Restore the scan order of the saved HSPs (the score sort below is
                // stable, so the order of comparator-equal HSPs matters).
                crate::algorithm::tblastx::x_seed_bucket::restore_scan_order(
                    &mut init_hsps,
                    &seq_keys,
                );

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
                    get_ungapped_hsp_list(init_hsps, contexts_ref, &s_frames)
                };
                adjust_tblastx_chunk_subject_offsets(&mut hits, chunk.offset);
                TblastxChunkScanResult { chunk, hits, stats }
            };

            // EXPERIMENT (LOSAT_X_SEEDBUCKET): 0 = scan order, 1 = bucketed (only for
            // tables of at least LOSAT_X_SEEDBUCKET_MIN_CELLS cells), 2 = both (the
            // bucketed one on a copy of the diagonal table), compared.
            let scan_chunk = |chunk: SubjectChunk,
                              offset_pairs: &mut [OffsetPair],
                              diag_array: &mut [DiagStruct],
                              diag_offset: &mut i32,
                              scan_chunk_size: Option<usize>|
             -> TblastxChunkScanResult {
                // The debug/trace paths exist only in the reference loop.
                let x_debugging = scan_debug_range.is_some() || trace_hsp_target().is_some();
                match crate::algorithm::tblastx::x_seed_bucket::mode_for(diag_array.len()) {
                    1 if !x_debugging => scan_chunk_bucketed(
                        chunk,
                        offset_pairs,
                        diag_array,
                        diag_offset,
                        scan_chunk_size,
                    ),
                    2 => {
                        let mut diag_copy = diag_array.to_vec();
                        let mut offset_copy = *diag_offset;
                        let shadow = scan_chunk_bucketed(
                            chunk,
                            offset_pairs,
                            &mut diag_copy,
                            &mut offset_copy,
                            scan_chunk_size,
                        );
                        let reference = scan_chunk_plain(
                            chunk,
                            offset_pairs,
                            diag_array,
                            diag_offset,
                            scan_chunk_size,
                        );
                        assert_eq!(
                            offset_copy, *diag_offset,
                            "LOSAT_X_SEEDBUCKETSHADOW: diag offset differs"
                        );
                        assert!(
                            diag_copy
                                .iter()
                                .zip(diag_array.iter())
                                .all(|(a, b)| a.raw_bits() == b.raw_bits()),
                            "LOSAT_X_SEEDBUCKETSHADOW: diagonal table differs after the chunk"
                        );
                        assert_eq!(
                            shadow.hits.len(),
                            reference.hits.len(),
                            "LOSAT_X_SEEDBUCKETSHADOW: hit count differs"
                        );
                        for (k, (a, b)) in shadow.hits.iter().zip(reference.hits.iter()).enumerate()
                        {
                            assert!(
                                a.q_idx == b.q_idx
                                    && a.s_idx == b.s_idx
                                    && a.ctx_idx == b.ctx_idx
                                    && a.s_f_idx == b.s_f_idx
                                    && a.q_frame == b.q_frame
                                    && a.s_frame == b.s_frame
                                    && a.q_aa_start == b.q_aa_start
                                    && a.q_aa_end == b.q_aa_end
                                    && a.s_aa_start == b.s_aa_start
                                    && a.s_aa_end == b.s_aa_end
                                    && a.q_seed_off == b.q_seed_off
                                    && a.s_seed_off == b.s_seed_off
                                    && a.q_orig_len == b.q_orig_len
                                    && a.s_orig_len == b.s_orig_len
                                    && a.raw_score == b.raw_score
                                    && a.e_value.to_bits() == b.e_value.to_bits()
                                    && a.num_ident == b.num_ident
                                    && a.hsp_list_order == b.hsp_list_order,
                                "LOSAT_X_SEEDBUCKETSHADOW: hit {k} differs: {a:?} vs {b:?}"
                            );
                        }
                        crate::algorithm::tblastx::x_seed_bucket::SHADOW_CHUNKS
                            .fetch_add(1, AtomicOrdering::Relaxed);
                        crate::algorithm::tblastx::x_seed_bucket::SHADOW_HITS
                            .fetch_add(reference.hits.len() as u64, AtomicOrdering::Relaxed);
                        reference
                    }
                    _ => scan_chunk_plain(
                        chunk,
                        offset_pairs,
                        diag_array,
                        diag_offset,
                        scan_chunk_size,
                    ),
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
                            let result =
                                scan_chunk(chunk, offset_pairs, diag_array, &mut diag_offset, None);
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
                                        vec![OffsetPair::default(); offset_array_size as usize],
                                        x_new_diag_array(diag_array_size as usize),
                                    )
                                },
                                |state, chunk| {
                                    let (offset_pairs, diag_array) = state;
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
                        scan_chunk(chunk, offset_pairs, diag_array, &mut diag_offset, None);
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
        let mut subject_frame_bases: Vec<i32> = Vec::with_capacity(s_frames.len());
        let mut base: i32 = 0;
        for f in &s_frames {
            subject_frame_bases.push(base);
            // NCBI concatenation shares the trailing NULLB sentinel between frames:
            //   offset += length + 1;
            // where `length` is the number of residues (excluding sentinels).
            // Source: ncbi-blast/c++/src/algo/blast/core/blast_util.c:1098-1101
            base += f.aa_seq.len() as i32 - 1;
        }

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
                use_parallel,
            )
        };
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
        // (The sort at the end of the second BLAST_LinkHsps is `report::final_hit_order`.)
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
        let mut ungapped_hits = reevaluate_ungapped_hsp_list(
            prelinked_ungapped_hits,
            contexts_ref,
            &s_frames,
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
            let linked = apply_sum_stats_even_gap_linking_with_parallel(
                ungapped_hits,
                &linking_params_for_cutoff,
                &linking_params,
                contexts_ref,
                &subject_frame_bases,
                &length_adj_per_context,
                &eff_searchsp_per_context,
                use_parallel,
            );
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

                let ctx = &contexts_ref[h.ctx_idx];
                let s_score_frame = &s_frames[h.s_f_idx];
                let s_frame = &s_frames_report[h.s_f_idx];

                let len = h.q_aa_end.saturating_sub(h.q_aa_start);
                let q0 = h.q_aa_start;
                let s0 = h.s_aa_start;

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
                let matches = h.num_ident;
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
                let bit_params = &contexts_ref[h.ctx_idx].karlin_params;
                let bit = calc_bit_score(h.raw_score, bit_params);
                let (q_start, q_end) =
                    convert_coords(h.q_aa_start, h.q_aa_end, ctx.frame, ctx.orig_len);
                let (s_start, s_end) =
                    convert_coords(h.s_aa_start, h.s_aa_end, s_frame.frame, s_len);

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
                if std::env::var_os("LOSAT_TRACE_HSP_MASKS").is_some() {
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
                                masked_score +=
                                    crate::utils::matrix::blosum62_score(q_masked, subject_scoring);
                                unmasked_score += crate::utils::matrix::blosum62_score(
                                    q_unmasked,
                                    subject_scoring,
                                );
                                report_unmasked_score += crate::utils::matrix::blosum62_score(
                                    q_unmasked,
                                    subject_report,
                                );
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
                    mismatch,
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
                    sort_query_offset: h.q_aa_start,
                    sort_query_end: h.q_aa_end,
                    sort_subject_offset: h.s_aa_start,
                    sort_subject_end: h.s_aa_end,
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
                // Single-threaded path accumulates hits directly, parallel path uses a channel.
                if let Some(tx) = &st.tx {
                    tx.send(final_hits).unwrap();
                } else {
                    st.hits.extend(final_hits);
                }
            }
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
    };

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
    let mut single_state: Option<WorkerState> = None;
    #[cfg(all(feature = "parallel", target_arch = "wasm32", feature = "wasm-threads"))]
    let mut threaded_wasi_subject_hit_batches: Option<Vec<(usize, Vec<TblastxHsp>)>> = None;

    #[cfg(all(feature = "parallel", target_arch = "wasm32", feature = "wasm-threads"))]
    if use_parallel && !use_serial_scan_chunks && subjects_raw.len() > 1 {
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
        //
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_seqalign.cpp:1569-1577
        // ```c
        // for (int index = 0; index < hit_list->hsplist_count; index++) {
        //     BlastHSPList* hsp_list = hit_list->hsplist_array[index];
        //     if (!hsp_list)
        //         continue;
        //     Blast_HSPListSortByEvalue(hsp_list);
        // }
        // ```
        let parallel_pool = parallel_pool;
        let subject_hit_batches: Vec<(usize, Vec<TblastxHsp>)> = parallel_pool.install(|| {
            subjects_raw
                .par_iter()
                .enumerate()
                .map_init(
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
                        offset_pairs: vec![OffsetPair::default(); offset_array_size as usize],
                        diag_array: x_new_diag_array(diag_array_size as usize),
                        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:52-63
                        // ```c
                        // diag_table->diag_array_length = diag_array_length;
                        // diag_table->diag_mask = diag_array_length-1;
                        // diag_table->offset = window_size;
                        // ```
                        diag_offset: window,
                    },
                    |st, (s_idx, s_rec)| {
                        process_subject(st, (s_idx, s_rec));
                        if st.hits.is_empty() {
                            None
                        } else {
                            Some((s_idx, std::mem::take(&mut st.hits)))
                        }
                    },
                )
                .filter_map(|batch| batch)
                .collect()
        });
        threaded_wasi_subject_hit_batches = Some(subject_hit_batches);
    } else {
        single_state = Some(for_each_subjects(
            &subjects_raw,
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
                offset_pairs: vec![OffsetPair::default(); offset_array_size as usize],
                diag_array: x_new_diag_array(diag_array_size as usize),
                // NCBI: diag_table->offset = window_size;
                // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:63
                diag_offset: window,
            },
            &process_subject,
        ));
    }

    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    if use_parallel && !use_serial_scan_chunks && subjects_raw.len() > 1 {
        let parallel_pool = parallel_pool;
        parallel_pool.install(|| {
            subjects_raw.par_iter().enumerate().for_each_init(
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
                    offset_pairs: vec![OffsetPair::default(); offset_array_size as usize],
                    diag_array: x_new_diag_array(diag_array_size as usize),
                    // NCBI: diag_table->offset = window_size;
                    // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:63
                    diag_offset: window,
                },
                &process_subject,
            );
        });
    } else {
        single_state = Some(for_each_subjects(
            &subjects_raw,
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
                offset_pairs: vec![OffsetPair::default(); offset_array_size as usize],
                diag_array: x_new_diag_array(diag_array_size as usize),
                // NCBI: diag_table->offset = window_size;
                // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:63
                diag_offset: window,
            },
            &process_subject,
        ));
    }

    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    {
        single_state = Some(for_each_subjects(
            &subjects_raw,
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
                offset_pairs: vec![OffsetPair::default(); offset_array_size as usize],
                diag_array: x_new_diag_array(diag_array_size as usize),
                // NCBI: diag_table->offset = window_size;
                // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:63
                diag_offset: window,
            },
            &process_subject,
        ));
    }

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
    // Each traversal ends with the final hits of the search, filtered by e-value;
    // they are sorted into NCBI's order and written at one place below.
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let collected_hits: Option<Vec<TblastxHsp>> = writer.map(|writer| writer.join().unwrap());
    #[cfg(not(all(feature = "parallel", not(target_arch = "wasm32"))))]
    let collected_hits: Option<Vec<TblastxHsp>> = None;
    #[cfg(all(feature = "parallel", target_arch = "wasm32", feature = "wasm-threads"))]
    let threaded_wasi_hits: Option<Vec<TblastxHsp>> =
        threaded_wasi_subject_hit_batches.map(|mut subject_hit_batches| {
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1427
            // ```c
            // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
            // /* iterate over all subject sequences */
            // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
            //        != BLAST_SEQSRC_EOF) {
            //    ...
            //    status = s_BlastSearchEngineCore(...);
            // }
            // ```
            subject_hit_batches.sort_by_key(|(s_idx, _)| *s_idx);
            let mut all: Vec<TblastxHsp> = Vec::new();
            for (_, hits) in subject_hit_batches {
                all.extend(hits);
            }
            // NCBI reference: c++/src/algo/blast/core/blast_hits.c:1988-1996
            // ```c
            //    cutoff = hit_options->expect_value;
            // ...
            //       if (hsp->evalue > cutoff) {
            // ```
            // An HSP is removed when its e-value is greater than the cutoff, so a NaN cutoff
            // (`-evalue -nan`) keeps every HSP.
            all.retain(|h| !(h.hit.e_value > evalue_threshold));
            all
        });
    #[cfg(not(all(feature = "parallel", target_arch = "wasm32", feature = "wasm-threads")))]
    let threaded_wasi_hits: Option<Vec<TblastxHsp>> = None;

    let final_hits = if let Some(all) = collected_hits.or(threaded_wasi_hits) {
        Some(all)
    } else if let Some(rx) = rx_opt.take() {
        let mut all: Vec<TblastxHsp> = Vec::new();
        for h in rx {
            all.extend(h);
        }
        // NCBI reference: c++/src/algo/blast/core/blast_hits.c:1988-1996
        // ```c
        //    cutoff = hit_options->expect_value;
        // ...
        //       if (hsp->evalue > cutoff) {
        // ```
        // An HSP is removed when its e-value is greater than the cutoff, so a NaN cutoff
        // (`-evalue -nan`) keeps every HSP.
        all.retain(|h| !(h.hit.e_value > evalue_threshold));
        Some(all)
    } else if let Some(state) = single_state {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1497
        // ```c
        // while (...) {
        //    ...
        //    status = s_BlastSearchEngineCore(..., &hsp_list, ...);
        // }
        // ```
        let mut all = state.hits;
        // NCBI reference: c++/src/algo/blast/core/blast_hits.c:1988-1996
        // ```c
        //    cutoff = hit_options->expect_value;
        // ...
        //       if (hsp->evalue > cutoff) {
        // ```
        // An HSP is removed when its e-value is greater than the cutoff, so a NaN cutoff
        // (`-evalue -nan`) keeps every HSP.
        all.retain(|h| !(h.hit.e_value > evalue_threshold));
        Some(all)
    } else {
        None
    };
    // The hit list of each query (blast_hits.c Blast_HitListUpdate and its final e-value
    // sort) is applied over all batches in `report::final_hit_order`.
    let hits = final_hits.unwrap_or_default();
    if diag_enabled {
        print_diagnostics_summary(&diagnostics);
    }

    if timing_enabled {
        let t_search = t_search_start.elapsed();
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

        eprintln!(
            "[TIMING] read_queries: {:.3}s",
            t_read_queries.as_secs_f64()
        );
        eprintln!(
            "[TIMING] build_lookup: {:.3}s",
            t_build_lookup.as_secs_f64()
        );
        eprintln!(
            "[TIMING] read_subjects: {:.3}s",
            t_read_subjects.as_secs_f64()
        );
        eprintln!("[TIMING] scan_subject: {:.3}s (calls={})", scan_s, scan_n);
        eprintln!(
            "[TIMING] ungapped_extend: {:.3}s (calls={})",
            ungapped_s, ungapped_n
        );
        eprintln!("[TIMING] reevaluate: {:.3}s (calls={})", reeval_s, reeval_n);
        eprintln!(
            "[TIMING] sum_stats_linking: {:.3}s (calls={})",
            linking_s, linking_n
        );
        eprintln!(
            "[TIMING] identity_calc: {:.3}s (calls={})",
            identity_s, identity_n
        );
        eprintln!("[TIMING] search_total: {:.3}s", t_search.as_secs_f64());
        eprintln!("[TIMING] total: {:.3}s", t_total.elapsed().as_secs_f64());
    }
    Ok(TblastxBatch {
        hits,
        searched: true,
        queries: query_stats,
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
            ctx_idx,
            s_f_idx: 0,
            q_frame: 1,
            s_frame: 1,
            q_aa_start: q_start,
            q_aa_end: q_end,
            s_aa_start: s_start,
            s_aa_end: s_end,
            q_seed_off: q_start,
            s_seed_off: s_start,
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
