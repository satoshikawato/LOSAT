//! Main run() function for backbone mode
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c

use super::*;
use std::sync::Arc;

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
        ncbi_qsort_ungapped_hits_by_score(combined);
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

struct TblastxInMemoryRun<'a> {
    queries: Vec<fasta::Record>,
    subjects: Vec<fasta::Record>,
    output: &'a mut Vec<u8>,
}

#[cfg(target_arch = "wasm32")]
fn fasta_records_from_bytes(bytes: &[u8]) -> Vec<fasta::Record> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:486-651
    // ```c
    // void
    // SetupQueries_OMF(IBlastQuerySource& queries,
    //                  BlastQueryInfo* qinfo,
    //                  BLAST_SequenceBlk** seqblk,
    //                  EBlastProgramType prog,
    //                  ...)
    // ```
    fasta::Reader::new(bytes)
        .records()
        .filter_map(|record| record.ok())
        .collect()
}

#[cfg(target_arch = "wasm32")]
pub fn run_web_pair(args: TblastxArgs, query_fasta: &str, subject_fasta: &str) -> Result<Vec<u8>> {
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
    let queries = fasta_records_from_bytes(query_fasta.as_bytes());
    let subjects = fasta_records_from_bytes(subject_fasta.as_bytes());
    let mut output = Vec::new();
    run_internal(
        args,
        Some(TblastxInMemoryRun {
            queries,
            subjects,
            output: &mut output,
        }),
    )?;
    Ok(output)
}

pub fn run(args: TblastxArgs) -> Result<()> {
    run_internal(args, None)
}

fn run_internal(args: TblastxArgs, mut in_memory: Option<TblastxInMemoryRun<'_>>) -> Result<()> {
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
    let num_threads = {
        #[cfg(target_arch = "wasm32")]
        {
            1
        }
        #[cfg(not(target_arch = "wasm32"))]
        {
            #[cfg(feature = "parallel")]
            {
                if args.num_threads == 0 {
                    num_cpus::get()
                } else {
                    args.num_threads
                }
            }
            #[cfg(not(feature = "parallel"))]
            {
                1
            }
        }
    };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:82-88
    // ```c
    // if (num_threads > 1) {
    //     SetNumberOfThreads(num_threads);
    // }
    // ```
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let use_parallel = num_threads > 1;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1002-1003
    // ```c
    // status = s_BlastSetUpAuxStructures(..., aux_struct);
    // ```
    //
    // The scan-chunk experiment keeps the canonical two-hit state serial inside
    // each NCBI WordFinder search unit, so subject-level parallelism is disabled
    // while this gate is active.
    let use_parallel_scan_chunks = std::env::var_os("LOSAT_TBLASTX_PARALLEL_SCAN_CHUNKS").is_some();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:82-88
    // ```c
    // if (num_threads > 1) {
    //     SetNumberOfThreads(num_threads);
    // }
    // ```
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let use_parallel_chunks = use_parallel
        && !use_parallel_scan_chunks
        && std::env::var_os("LOSAT_TBLASTX_PARALLEL_CHUNKS").is_some();
    #[cfg(any(not(feature = "parallel"), target_arch = "wasm32"))]
    let use_parallel_chunks = false;
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

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:82-88
    // ```c
    // if (num_threads > 1) {
    //     SetNumberOfThreads(num_threads);
    // }
    // ```
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    if use_parallel {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .context("Failed to build thread pool")?;
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
    let queries_raw: Vec<fasta::Record> = if let Some(input) = in_memory.as_mut() {
        std::mem::take(&mut input.queries)
    } else {
        let query_reader = fasta::Reader::from_file(&args.query)?;
        query_reader.records().filter_map(|r| r.ok()).collect()
    };
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
    //                       Set to 0 if not applicable */
    //    BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
    //    Int4 hspcnt; /**< Number of HSPs saved */
    //    ...
    // } BlastHSPList;
    // ```
    let query_ids: Vec<Arc<str>> = queries_raw
        .iter()
        .map(|r| Arc::<str>::from(r.id().split_whitespace().next().unwrap_or("unknown")))
        .collect();
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

    if args.seg {
        let seg = SegMasker::new(args.seg_window, args.seg_locut, args.seg_hicut);
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
    let (lookup, contexts) = build_ncbi_lookup(
        &query_frames,
        args.threshold,
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
    let subjects_raw: Vec<fasta::Record> = if let Some(input) = in_memory.as_mut() {
        std::mem::take(&mut input.subjects)
    } else {
        let subject_reader = fasta::Reader::from_file(&args.subject)?;
        subject_reader.records().filter_map(|r| r.ok()).collect()
    };
    if queries_raw.is_empty() || subjects_raw.is_empty() {
        return Ok(());
    }
    t_read_subjects = t_phase_read_subjects.elapsed();

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
    //                       Set to 0 if not applicable */
    //    BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
    //    Int4 hspcnt; /**< Number of HSPs saved */
    //    ...
    // } BlastHSPList;
    // ```
    let subject_ids: Vec<Arc<str>> = subjects_raw
        .iter()
        .map(|r| Arc::<str>::from(r.id().split_whitespace().next().unwrap_or("unknown")))
        .collect();

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
    // Single-threaded NCBI runs the subject loop in-process, so we can
    // accumulate hits directly without an mpsc queue.
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let use_channel = use_parallel && !use_parallel_scan_chunks;
    #[cfg(any(not(feature = "parallel"), target_arch = "wasm32"))]
    let use_channel = false;

    let (tx_opt, mut rx_opt) = if use_channel {
        let (tx, rx) = channel::<Vec<Hit>>();
        (Some(tx), Some(rx))
    } else {
        (None, None)
    };
    let out_path = args.out.clone();
    let evalue_threshold = args.evalue;

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
        let query_ids_out = query_ids.clone();
        let subject_ids_out = subject_ids.clone();
        let out_path = out_path.clone();
        let rx = rx_opt.take().expect("rx must be available for writer");
        Some(std::thread::spawn(move || -> Result<()> {
            let mut all: Vec<Hit> = Vec::new();
            while let Ok(h) = rx.recv() {
                all.extend(h);
            }
            all.retain(|h| h.e_value <= evalue_threshold);
            // NCBI-style output ordering: query (input order) → subject (best_evalue/score/oid) → HSP (score/coords)
            // Reference: BLAST_LinkHsps() + s_EvalueCompareHSPLists() + ScoreCompareHSPs()
            write_output_ncbi_order(all, out_path.as_ref(), &query_ids_out, &subject_ids_out)?;
            Ok(())
        }))
    } else {
        None
    };

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
        subjects: &[fasta::Record],
        init: FInit,
        mut body: FBody,
    ) -> WorkerState
    where
        FInit: FnOnce() -> WorkerState,
        FBody: FnMut(&mut WorkerState, (usize, &fasta::Record)),
    {
        let mut state = init();
        for (s_idx, s_rec) in subjects.iter().enumerate() {
            body(&mut state, (s_idx, s_rec));
        }
        state
    }

    let process_subject = |st: &mut WorkerState, (s_idx, s_rec): (usize, &fasta::Record)| {
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

        let s_len = s_rec.seq().len();

        // [C] BlastOffsetPair *offset_pairs
        let offset_pairs = &mut st.offset_pairs;

        // [C] DiagStruct *diag_array = diag->hit_level_array;
        let diag_array = &mut st.diag_array;

        // [C] diag_offset = diag->offset;  (reset to window per-subject)
        let mut diag_offset: i32 = st.diag_offset;

        // Precompute per-context cutoff scores using NCBI BLAST algorithm.
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:280-419
        //
        // NCBI cutoff calculation for tblastx ungapped path:
        // 1. gap_trigger from ungapped params (kbp_std)
        // 2. cutoff_score_max from BlastHitSavingParametersNew (uses user's E-value)
        // 3. Per-subject update: cutoff_score_for_update_tblastx with CUTOFF_E_TBLASTX=1e-300
        // 4. Final cutoff = MIN(update_cutoff, gap_trigger, cutoff_score_max)
        //
        // For -subject mode, subject_len_nucl is used (not per-frame AA length).
        let subject_len_nucl = s_len as i64;
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
        // We precompute them here and pass to sum_stats_linking for NCBI parity.
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:846-847
        //   query_info->contexts[index].eff_searchsp = effective_search_space;
        //   query_info->contexts[index].length_adjustment = length_adjustment;
        let mut length_adj_per_context: Vec<i64> = Vec::with_capacity(contexts_ref.len());
        let mut eff_searchsp_per_context: Vec<i64> = Vec::with_capacity(contexts_ref.len());

        for (ctx_idx, ctx) in contexts_ref.iter().enumerate() {
            // NCBI: per-context kbp_std[context] is used throughout cutoff/length calcs.
            // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2797
            let ctx_params = &ctx.karlin_params;

            // NCBI: query_length = query_info->contexts[context].query_length
            let query_len_aa = ctx.aa_len as i64;

            // NCBI: gap_trigger uses kbp_std[context]->Lambda/logK.
            // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:340-345
            let gap_trigger = gap_trigger_raw_score(GAP_TRIGGER_BIT_SCORE, ctx_params);

            // =======================================================================
            // NCBI Parity: Use compute_eff_lengths_subject_mode_tblastx to get BOTH
            // length_adjustment and eff_searchsp from a single source of truth.
            // This mirrors BLAST_CalcEffLengths which computes and stores both values.
            // Reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:821-847
            // =======================================================================
            let eff_lengths = compute_eff_lengths_subject_mode_tblastx(
                query_len_aa,
                subject_len_nucl,
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
            static PRINTED: std::sync::atomic::AtomicBool =
                std::sync::atomic::AtomicBool::new(false);
            if diag_enabled
                && ctx_idx == 0
                && !PRINTED.swap(true, std::sync::atomic::Ordering::Relaxed)
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
        // If no contexts, use 0 as fallback
        if cutoff_score_min == i32::MAX {
            cutoff_score_min = 0;
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

        for (s_f_idx, s_frame) in s_frames.iter().enumerate() {
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

            let scan_chunk = |chunk: SubjectChunk,
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

            if use_parallel_scan_chunks {
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
                #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
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
                                        vec![DiagStruct::default(); diag_array_size as usize],
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
                #[cfg(any(not(feature = "parallel"), target_arch = "wasm32"))]
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
                    ncbi_qsort_ungapped_hits_by_score(&mut combined_ungapped_hits);
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

        // NCBI CalculateLinkHSPCutoffs parameters. These are computed before
        // both BLAST_LinkHsps calls below, matching the single per-subject
        // CalculateLinkHSPCutoffs call in blast_engine.c:1448-1454.
        let linking_params = LinkingParams {
            avg_query_length,
            subject_len_nucl,
            cutoff_score_min,
            scale_factor: 1.0,   // Standard BLOSUM62
            gap_decay_rate: 0.5, // BLAST_GAP_DECAY_RATE
        };
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

        // NCBI parity: s_BlastFindSmallestLambda selects smallest lambda across contexts.
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:92-112
        // Per-context kbp_std is computed from query composition (with check_ideal).
        // Reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2778-2797
        let context_params: Vec<KarlinParams> =
            contexts_ref.iter().map(|ctx| ctx.karlin_params).collect();
        let linking_params_for_cutoff =
            find_smallest_lambda_params(&context_params).unwrap_or_else(|| params.clone());

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
            apply_sum_stats_even_gap_linking(
                combined_ungapped_hits,
                &linking_params_for_cutoff,
                &linking_params,
                contexts_ref,
                &subject_frame_bases,
                &length_adj_per_context,
                &eff_searchsp_per_context,
            )
        };
        if stage_dump::enabled() {
            stage_dump::dump_ungapped_hits(
                "after_prelim_link_hsps_before_reevaluate",
                &prelinked_ungapped_hits,
            );
        }
        let prelim_linked_count = prelinked_ungapped_hits.len();
        prelinked_ungapped_hits.retain(|h| h.e_value <= evalue_threshold);
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
            ncbi_qsort_ungapped_hits_by_score(&mut ungapped_hits);
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
            let linked = apply_sum_stats_even_gap_linking(
                ungapped_hits,
                &linking_params_for_cutoff,
                &linking_params,
                contexts_ref,
                &subject_frame_bases,
                &length_adj_per_context,
                &eff_searchsp_per_context,
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

            // NCBI HSP culling (conditional on culling_limit > 0)
            // Reference: hspfilter_culling.c - applied after linking, before output
            // Default for tblastx: culling_limit = 0 (disabled)
            let linked_before_cull = linked.len();
            let linked = if args.culling_limit > 0 {
                hsp_culling::apply_culling(linked, contexts_ref, args.culling_limit)
            } else {
                // Default: no culling (matches NCBI tblastx default)
                linked
            };
            if diag_enabled && args.culling_limit > 0 {
                diagnostics.base.hsps_culled_dominated.fetch_add(
                    linked_before_cull.saturating_sub(linked.len()),
                    AtomicOrdering::Relaxed,
                );
            }

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

            let mut final_hits: Vec<Hit> = Vec::new();
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
                    gap_info: None,
                    num_positives: matches,
                };
                trace_final_hit_if_match("output_hit", &out_hit);
                if dump_output_stage {
                    output_snapshot_pairs.push((out_hit.clone(), h.clone()));
                }
                final_hits.push(out_hit);
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

    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    if use_parallel && !use_parallel_scan_chunks {
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
                diag_array: vec![DiagStruct::default(); diag_array_size as usize],
                // NCBI: diag_table->offset = window_size;
                // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:63
                diag_offset: window,
            },
            &process_subject,
        );
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
                diag_array: vec![DiagStruct::default(); diag_array_size as usize],
                // NCBI: diag_table->offset = window_size;
                // Source: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:63
                diag_offset: window,
            },
            &process_subject,
        ));
    }

    #[cfg(any(not(feature = "parallel"), target_arch = "wasm32"))]
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
                diag_array: vec![DiagStruct::default(); diag_array_size as usize],
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
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    if let Some(writer) = writer {
        writer.join().unwrap()?;
    }
    if let Some(rx) = rx_opt.take() {
        let mut all: Vec<Hit> = Vec::new();
        for h in rx {
            all.extend(h);
        }
        all.retain(|h| h.e_value <= evalue_threshold);
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1700-1704
        // ```c
        // s_Heapify((char*)hsp_array, (char*)hsp_array,
        //          (char*)&hsp_array[hsp_list->hspcnt/2 - 1],
        //          (char*)&hsp_array[hsp_list->hspcnt-1],
        //          sizeof(BlastHSP*), ScoreCompareHSPs);
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3202-3205
        // ```c
        // s_Heapify((char*)hit_list->hsplist_array, (char*)hit_list->hsplist_array,
        //          (char*)&hit_list->hsplist_array[hit_list->hsplist_count/2 - 1],
        //          (char*)&hit_list->hsplist_array[hit_list->hsplist_count-1],
        //          sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
        // ```
        write_output_ncbi_order(all, out_path.as_ref(), &query_ids, &subject_ids)?;
    } else if let Some(state) = single_state {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1497
        // ```c
        // while (...) {
        //    ...
        //    status = s_BlastSearchEngineCore(..., &hsp_list, ...);
        // }
        // ```
        let mut all = state.hits;
        all.retain(|h| h.e_value <= evalue_threshold);
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1700-1704
        // ```c
        // s_Heapify((char*)hsp_array, (char*)hsp_array,
        //          (char*)&hsp_array[hsp_list->hspcnt/2 - 1],
        //          (char*)&hsp_array[hsp_list->hspcnt-1],
        //          sizeof(BlastHSP*), ScoreCompareHSPs);
        // ```
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3202-3205
        // ```c
        // s_Heapify((char*)hit_list->hsplist_array, (char*)hit_list->hsplist_array,
        //          (char*)&hit_list->hsplist_array[hit_list->hsplist_count/2 - 1],
        //          (char*)&hit_list->hsplist_array[hit_list->hsplist_count-1],
        //          sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
        // ```
        if let Some(input) = in_memory.as_mut() {
            write_output_ncbi_order_to_writer(all, input.output, &query_ids, &subject_ids)?;
        } else {
            write_output_ncbi_order(all, out_path.as_ref(), &query_ids, &subject_ids)?;
        }
    }
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
    Ok(())
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
