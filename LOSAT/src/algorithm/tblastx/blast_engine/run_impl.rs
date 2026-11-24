//! Main run() function for backbone mode
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c

use super::*;
use crate::algorithm::common::translated_chunks::{translated_chunk_scan_ranges, TranslatedChunk};
use crate::algorithm::common::x_db_stream::XReadError;
use crate::algorithm::tblastx::database::DatabaseTitles;
use crate::algorithm::tblastx::hit_collector::HitCollector;
use crate::algorithm::tblastx::subjects::{SubjectSetNumbers, Subjects, SubjectsUnreadable};
use crate::blastinput::report_counts::ReportCounts;
use crate::report::diag_log::Posts;
use crate::report::tabular::SubjectFetch;
use crate::wording::ENGINE_ERROR;
use std::sync::Arc;

use crate::algorithm::blastn::blast_engine::{QueryReading, NO_DATA_MESSAGE};
use crate::blastinput::fasta_reader::{
    read_queries, read_subjects, FastaInputSource, FastaRecord, QueryEnd, QueryRecords,
    ReaderConfig,
};
use crate::blastinput::query_batch::BatchCount;
use crate::blastinput::seq_range::{
    cut_queries, cut_subjects, parse_optional_range, subjects_read, Placements, QueryInput,
    RangeRole, SequenceRange,
};

// EXPERIMENT (PAIRPAR, graduated on 2026-10-22: no variable): the seed stage of a subject as units
// of (subject frame, chunk, query context) on the search pool; a child module of this file so that
// it uses the chunk helpers below unchanged (see tblastx/x_pair_par.rs).
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-516
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
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:484-658
// (`SetupQueries_OMF`) and c++/src/algo/blast/core/lookup_wrap.c:61-195 (`LookupTableWrapInit_MT`).
// The flags, the query frames and the lookup table of a query batch (a child module, so that it
// uses the private modules of this file).
#[path = "batch_setup.rs"]
mod batch_setup;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:279-419
// (`BlastInitialWordParametersUpdate`) and c++/src/algo/blast/core/blast_setup.c:699-850
// (`BLAST_CalcEffLengths`).
// The cutoffs, length adjustments and search spaces of a query batch's contexts (a child module,
// so that it uses the private error type of this file).
#[path = "context_cutoffs.rs"]
mod context_cutoffs;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:41-67 (`s_BlastDiagTableNew`)
// and c++/src/algo/blast/core/lookup_wrap.c:255-289 (`GetOffsetArraySize`).
// The diagonal table and offset-array sizes and the seed-stage tables of a query batch (a child
// module, so that it uses the private modules of this file).
#[path = "seed_tables.rs"]
mod seed_tables;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1584 (the subject loop of
// `BLAST_PreliminarySearchEngine`) and c++/src/algo/blast/api/prelim_stage.cpp:145-188 (its
// threads).
// The worker states, the writer thread and the subject loop of a query batch (a child module, so
// that it uses the private modules of this file).
#[path = "subject_loop.rs"]
mod subject_loop;
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1569-1577 (in
// `BlastHitList2SeqAlign_OMF`).
// The hit lists and the `LOSAT_TIMING` lines at the end of a query batch's search (a child module,
// so that it uses the private timer of this file).
#[path = "batch_finish.rs"]
mod batch_finish;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-899 and :1492-1520 (the
// linking calls, the reap and the re-evaluation of a subject's HSP list).
// The end of a subject's search up to its linked HSP list (a child module, so that it uses the
// private modules of this file).
#[path = "subject_link.rs"]
mod subject_link;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1535-1541 (the reap and the
// bit scores of a subject's final HSPs).
// The report hits of a subject (a child module, so that it uses the private timer of this file).
#[path = "subject_report.rs"]
mod subject_report;
// No NCBI counterpart: the diagnostics at the end of a subject's search (a child module, as
// subject_link is).
#[path = "subject_debug.rs"]
mod subject_debug;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775 and :1429-1432 (the
// translation of a subject and the test of its length).
// The preliminary frames and the unit path of a subject (a child module, so that it uses the
// private modules of this file).
#[path = "subject_prep.rs"]
mod subject_prep;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-593 and :844 (the chunk
// loop of a subject frame and the append of the frame's list).
// The chunks of a subject frame (a child module, so that it uses the private chunk helpers of
// this file).
#[path = "frame_chunks.rs"]
mod frame_chunks;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:440-619
// (`s_BlastAaWordFinder_TwoHit`).
// The word finder that the two chunk scans of the frame loop share (a child module, so that it
// uses the private chunk types and helpers of this file).
#[path = "word_finder.rs"]
mod word_finder;
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
    /// (outfmt 0), and the stream of its sequences (`None` when they could not be read: nothing
    /// is reported), from which the batch loop reads the reported ones (`database::read_reported`)
    /// that the reports show by rank (`database::reported_subjects`).
    Database {
        database: &'a super::super::database::TblastxDatabase,
        pairwise: crate::report::pairwise::PairwiseDatabase,
        stream: Option<&'a crate::algorithm::common::x_db_stream::XDbStream<'a>>,
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

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:54
// The one copy of the macro (`core::subject_split`).
use crate::core::subject_split::MAX_DBSEQ_LEN;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:54
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-505
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-255
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1296-1308
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:75-127
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-310
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:83-123
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:816-831 (the ranges of context 0 serve the
// positive frames, those of context 3, by `subject->length` of frame -1, the negative ones) and
// :259-310 (one list per chunk; SUBJECT_SPLIT_NO_RANGE skips the chunk's WordFinder). The same
// ranges as TBLASTN's (`translated_chunk_scan_ranges`): a nucleotide subject's lower-case runs in
// `seq_ranges`, entered for both strands (blast_aux_priv.cpp:436-454).
/// The chunks of one frame of subject `s_idx` that have ranges, with their ranges, when the subject
/// has a lower-case letter. Cold and out of line (no NCBI counterpart: code placement only).
#[cold]
#[inline(never)]
fn x_soft_subject_calls(
    subjects: Subjects<'_>,
    s_idx: usize,
    frame: &QueryFrame,
    frames: &[QueryFrame],
) -> Option<Vec<(TranslatedChunk, Vec<(i32, i32)>)>> {
    let Subjects::Fasta(records) = subjects else {
        return None;
    };
    let letters = records[s_idx].seq();
    if !letters.iter().any(u8::is_ascii_lowercase) {
        return None;
    }
    // `translated_chunk_scan_ranges` lists the chunks of NCBI's MAX_DBSEQ_LEN.
    debug_assert_eq!(tblastx_max_dbseq_len_for_run(), MAX_DBSEQ_LEN);
    let negative_first_length = frames
        .iter()
        .find(|frame| frame.frame == -1)
        .map_or(0, |frame| frame.aa_len);
    Some(
        translated_chunk_scan_ranges(
            letters,
            frame.frame,
            frame.aa_len,
            negative_first_length,
            true,
        )
        .expect("the ranges of a soft-masked subject have their sentinel"),
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-830 (the frames' ranges
// mapped once per strand from the sequence's ranges in bases) and :221-311 (each frame chunked as
// a protein subject; with soft masks the ranges cut per chunk of the whole frame).
/// A masked `-db` search (`-db_soft_mask`, `-db_hard_mask`; series NP, M4b PR 2): the chunks of
/// frame `frame` of a subject whose ranges in bases are `ranges`, with their scan ranges
/// (`translated_db_mask_chunks`). Cold and out of line (no NCBI counterpart: code placement only).
#[cold]
#[inline(never)]
fn x_db_mask_calls(
    ranges: &[(i32, i32)],
    hard: bool,
    frame: &QueryFrame,
    frames: &[QueryFrame],
) -> Vec<(TranslatedChunk, Vec<(i32, i32)>)> {
    debug_assert_eq!(tblastx_max_dbseq_len_for_run(), MAX_DBSEQ_LEN);
    let negative_first_length = frames
        .iter()
        .find(|frame| frame.frame == -1)
        .map_or(0, |frame| frame.aa_len);
    crate::algorithm::common::translated_chunks::translated_db_mask_chunks(
        ranges,
        hard,
        frame.frame,
        frame.aa_len,
        negative_first_length,
        &mut crate::core::subject_split::HardSplitState::default(),
    )
    .expect("the ranges of a masked subject have their sentinel")
}

// No NCBI counterpart: looks up what `x_soft_subject_calls` computed.
/// The ranges of `chunk`; `None` for a chunk without ranges (SUBJECT_SPLIT_NO_RANGE).
#[cold]
#[inline(never)]
fn x_soft_chunk<'a>(
    calls: &'a [(TranslatedChunk, Vec<(i32, i32)>)],
    chunk: &SubjectChunk,
) -> Option<&'a [(i32, i32)]> {
    calls
        .iter()
        .find(|(call, _)| call.offset == chunk.offset)
        .map(|(_, ranges)| ranges.as_slice())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:816-831. A shorter frame takes the ranges
// of the longest frame of its strand, so the last word NCBI scans can end on the NULLB after the
// frame (blast_util.c:438-453, BLAST_GetAllTranslations).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:300-306,
// c++/src/algo/blast/core/blast_aascan.c:76-77 and c++/src/algo/blast/core/aa_ungapped.c:499-500
// NCBI's scan reads the subject's buffer wherever a range points: before the chunk for a later
// range of masks stored out of order, and past its last letter for a first range raised to one
// word at its left (the frame's next letters there, on a chunk that is not the frame's last).
/// The letters the scan of a chunk reads, and how many of them lie before the chunk: `subject`
/// alone; with `x_soft` ranges the frame from `common::translated_chunks::scan_lead` letters
/// before the chunk, to `wordsize - 1` letters (at least 1) after it, or to the sentinel after the
/// frame. The extensions of the hits read the same letters (`word_finder::ChunkCx`).
#[inline(always)]
fn x_soft_scan_letters<'a>(
    x_soft: Option<&[(i32, i32)]>,
    subject: &'a [u8],
    subject_full: &'a [u8],
    chunk: &SubjectChunk,
    wordsize: usize,
) -> (&'a [u8], usize) {
    match x_soft {
        Some(ranges) => {
            let lead = crate::algorithm::common::translated_chunks::scan_lead(ranges);
            let end = (1 + chunk.offset + chunk.length + wordsize.saturating_sub(1).max(1))
                .min(subject_full.len());
            (&subject_full[1 + chunk.offset - lead..end], lead)
        }
        None => (subject, 0),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:76-77 and :110
// NCBI scans from a range's left, which can lie before the chunk, and keeps each hit's offset
// from the chunk's first letter.
/// The ranges of a chunk as the scan of its letters from `lead` letters before it reads them
/// (`x_soft_scan_letters`): each moved by `lead`.
#[inline(always)]
fn x_scan_ranges(ranges: &[(i32, i32)], lead: usize) -> Vec<(i32, i32)> {
    let lead = lead as i32;
    ranges
        .iter()
        .map(|&(left, right)| (left + lead, right + lead))
        .collect()
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:190-192
// The one copy of the macro (`core::subject_split`).
use crate::core::subject_split::DBSEQ_CHUNK_OVERLAP;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1533-1539
const OVERLAP_DIAG_CLOSE: i32 = 10;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:122-143
struct SubjectSplitState {
    full_right: i32,
    hard_ranges: [(i32, i32); 1],
    hm_index: usize,
    offset: i32,
    next: i32,
    /// `-db_hard_mask` (series NP, M4b PR 2): the frame's chunks, listed in order
    /// (`x_db_mask_calls`), given instead of the cut of the whole frame.
    x_listed: Option<std::vec::IntoIter<SubjectChunk>>,
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:452-584
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:146-184
    fn new(subject_len: usize) -> Self {
        let full_right = subject_len as i32;
        Self {
            full_right,
            hard_ranges: [(0, full_right)],
            hm_index: 0,
            offset: 0,
            next: 0,
            x_listed: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-311 (the hard chunks
    // of `s_GetNextSubjectChunk`, cut in `translated_db_mask_chunks`) and :578-584 (the merge's
    // overlap of each).
    /// The state of a frame of `subject_len` residues whose chunks are `calls`' (a hard-masked
    /// frame), each cut to the frame: a chunk that runs onto the frame's end sentinel (a shorter
    /// frame reads its strand's ranges, mapped by the longest frame) scans and extends as the chunk
    /// that ends at the frame's end, since a word or an extension stops at the sentinel
    /// (`BLAST_SCORE_MIN`).
    fn listed(subject_len: usize, calls: &[(TranslatedChunk, Vec<(i32, i32)>)]) -> Self {
        let chunks: Vec<SubjectChunk> = calls
            .iter()
            .map(|(chunk, _)| SubjectChunk {
                offset: chunk.offset,
                length: chunk.length.min(subject_len.saturating_sub(chunk.offset)),
                overlap: chunk.overlap,
            })
            .collect();
        let mut state = Self::new(subject_len);
        state.x_listed = Some(chunks.into_iter());
        state
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-310
    fn next_chunk(&mut self, max_dbseq_len: usize, chunk_overlap: usize) -> SubjectChunkStatus {
        if let Some(listed) = &mut self.x_listed {
            return listed
                .next()
                .map_or(SubjectChunkStatus::Done, SubjectChunkStatus::Ok);
        }
        if self.next >= self.full_right {
            return SubjectChunkStatus::Done;
        }

        let residual = 0;
        self.offset = self.next - residual;
        let offset = self.offset;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:578-584
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits_priv.h:68-69
// On NCBI's `Int4` offsets and ends.
#[inline]
fn contained_in_hsp(a: i32, b: i32, c: i32, d: i32, e: i32, f: i32) -> bool {
    a <= c && b >= c && d <= f && e >= f
}

/// A merge of two HSPs from two chunks of a subject frame whose result NCBI leaves undefined
/// (`merge_two_tblastx_chunk_hsps`): their lengths add up to 0, which only an HSP of the one-hit
/// word finder of length 0 or less can give (`TWO_HIT_CHUNK_MERGE_IS_DEFINED`).
#[derive(Debug)]
struct ChunkMergeUndefined;

/// Why the chunk merge of a list of the two-hit word finder never returns `ChunkMergeUndefined`.
/// A merge never shortens an HSP: the merged HSP takes the minimum start and the maximum end, so its
/// length is at least the longer input's. A two-hit HSP merged again in a third chunk therefore
/// still has a positive length, and the lengths of two HSPs of any chunk add up to at least 2.
const TWO_HIT_CHUNK_MERGE_IS_DEFINED: &str = concat!(
    "the two-hit word finder saves HSPs of positive length only (`convert_init_hsps_unsorted` ",
    "asserts a positive span), so the lengths of two of them add up to at least 2 in a chunk merge"
);

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1488-1531 (`s_BlastMergeTwoHSPs`)
// TBLASTX is ungapped (`allow_gap` is FALSE: the diagonal test runs). The offsets and ends are
// NCBI's `Int4` (`UngappedHit::q_end_int4`): an HSP of the one-hit word finder (`-window_size 0`)
// may have a length of 0 or less (aa_ungapped.c:1054,1083), and the lengths are NCBI's
// differences, so such an HSP lowers the density's divisor: an HSP of length 0 contains the HSPs
// that start or end at its point, one of negative length contains nothing, and a positive HSP
// that absorbs one of negative length gets a higher density. The diagonal test is computed in
// `i64`, the value NCBI's `Int4` arithmetic gives where it does not overflow. Lengths that add up
// to 0 (an HSP of length 0 found by both chunks, or a positive HSP of length 1 or 2 that contains
// the start or end of one of length -1 or -2) make the density `+inf` (both scores are at least
// 1), and its product with the merged length `+inf`, or NaN when that length is 0: `(int)` of
// either is undefined (C11 6.3.1.4; x86-64 gives INT_MIN), `Err`, which the search refuses
// (PD-LOSAT-NCBI-DEFECTS 1.21, rule 3). Otherwise `(int)` of the finite product truncates toward
// 0, as `as i32` does for a value in range (the lengths are at least -2, so the product stays
// near the scores).
fn merge_two_tblastx_chunk_hsps(
    hsp1: &mut UngappedHit,
    hsp2: &UngappedHit,
) -> Result<bool, ChunkMergeUndefined> {
    let (q_off1, s_off1) = (hsp1.q_aa_start as i32, hsp1.s_aa_start as i32);
    let (q_off2, s_off2) = (hsp2.q_aa_start as i32, hsp2.s_aa_start as i32);
    if i64::from(s_off1) - i64::from(s_off2) - i64::from(q_off1) + i64::from(q_off2) != 0 {
        return Ok(false);
    }
    if hsp1.s_frame != hsp2.s_frame {
        return Ok(false);
    }
    let (q_end1, s_end1) = (hsp1.q_end_int4(), hsp1.s_end_int4());
    let (q_end2, s_end2) = (hsp2.q_end_int4(), hsp2.s_end_int4());
    if contained_in_hsp(q_off1, q_end1, q_off2, s_off1, s_end1, s_off2)
        || contained_in_hsp(q_off1, q_end1, q_end2, s_off1, s_end1, s_end2)
    {
        let length_sum = (q_end1 - q_off1) + (q_end2 - q_off2);
        if length_sum == 0 {
            return Err(ChunkMergeUndefined);
        }
        let score_density = (hsp1.raw_score as f64 + hsp2.raw_score as f64) / f64::from(length_sum);

        hsp1.q_aa_start = hsp1.q_aa_start.min(hsp2.q_aa_start);
        hsp1.s_aa_start = hsp1.s_aa_start.min(hsp2.s_aa_start);
        let q_end = q_end1.max(q_end2);
        hsp1.set_int4_ends(q_end, s_end1.max(s_end2));

        if hsp2.raw_score > hsp1.raw_score {
            hsp1.q_seed_off = hsp2.q_seed_off;
            hsp1.s_seed_off = hsp2.s_seed_off;
            hsp1.raw_score = hsp2.raw_score;
        }

        let new_length = q_end - hsp1.q_aa_start as i32;
        hsp1.raw_score = hsp1
            .raw_score
            .max((score_density * f64::from(new_length)) as i32);
        return Ok(true);
    }

    Ok(false)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2876-2881,2890-2911
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2956-3029
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2758-2766 (in
// `s_BlastHSPListsCombineByScore`; `hsp_num_max` is `kHspNumMax`, which keeps every HSP)
// The merge of a chunk's list (`incoming`, its offsets adjusted) into the frame's list of the
// chunks before it (`combined`), for every TBLASTX list. The swaps move the HSPs of each list that
// reach the overlap strip to its front, in their order, and the HSPs they displace behind later
// ones ([N1, N2, S] becomes [S, N2, N1]); the lists keep that order in the append, so the stable
// score sort (`Blast_HSPListSortByScore`, glibc's merge sort) leaves HSPs that tie on all its
// keys in NCBI's order. The ends are NCBI's `Int4` (the strip test of the earlier chunk's HSPs,
// `s_HSPEndDiag`, the merge and the sort), as a list of the one-hit word finder needs; on a list of
// the two-hit word finder every end is 0 or more, where the `Int4` and the `u32` are the same
// number. A merge NCBI leaves undefined stops the merging (`Err`): the search is refused, and its
// lists are not reported.
fn merge_tblastx_subject_chunk_hits(
    combined: &mut Vec<UngappedHit>,
    mut incoming: Vec<UngappedHit>,
    split_offset: usize,
    chunk_overlap_size: usize,
) -> Result<(), ChunkMergeUndefined> {
    if incoming.is_empty() {
        return Ok(());
    }
    if combined.is_empty() {
        combined.append(&mut incoming);
        return Ok(());
    }
    let split = i64::try_from(split_offset)
        .expect("a subject chunk's offset fits i64 (NCBI's split offset is Int4)");
    let strip_end = split
        + i64::try_from(chunk_overlap_size).expect("the chunk overlap fits i64 (NCBI's is Int4)");

    let mut hspcnt1 = 0usize;
    for index1 in 0..combined.len() {
        if i64::from(combined[index1].s_end_int4()) > split {
            combined.swap(hspcnt1, index1);
            hspcnt1 += 1;
        }
    }
    let mut hspcnt2 = 0usize;
    for index2 in 0..incoming.len() {
        if i64::from(incoming[index2].s_aa_start) < strip_end {
            incoming.swap(hspcnt2, index2);
            hspcnt2 += 1;
        }
    }

    let mut result = Ok(());
    if hspcnt1 > 0 && hspcnt2 > 0 {
        let mut deleted = vec![false; hspcnt2];
        'merge: for index1 in 0..hspcnt1 {
            for index2 in 0..hspcnt2 {
                if deleted[index2] || incoming[index2].ctx_idx != combined[index1].ctx_idx {
                    continue;
                }
                let end_diag = i64::from(combined[index1].q_end_int4())
                    - i64::from(combined[index1].s_end_int4());
                let start_diag =
                    i64::from(incoming[index2].q_aa_start) - i64::from(incoming[index2].s_aa_start);
                if (end_diag - start_diag).abs() < i64::from(OVERLAP_DIAG_CLOSE) {
                    match merge_two_tblastx_chunk_hsps(&mut combined[index1], &incoming[index2]) {
                        Ok(true) => deleted[index2] = true,
                        Ok(false) => {}
                        Err(undefined) => {
                            result = Err(undefined);
                            break 'merge;
                        }
                    }
                }
            }
        }
        let mut index2 = 0usize;
        incoming.retain(|_| {
            let keep = index2 >= hspcnt2 || !deleted[index2];
            index2 += 1;
            keep
        });
    }
    combined.append(&mut incoming);
    sort_ungapped_hits_by_score_int4(combined);
    result
}

/// `-subject_besthit` on a subject frame's list after the merge of a chunk that had HSPs, in every
/// merge of the frame's chunks (the frame loop's `frame_chunks::add_chunk_result` and the unit
/// path's merges, `x_pair_par` and `x_kmerge`): before the frame's list joins the subject's and
/// before the linking. The contexts are the query frames; the end of an HSP of length 0 or less is
/// NCBI's `Int4`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:564,586-591
fn frame_subject_besthit(frame_hits: &mut Vec<UngappedHit>) {
    crate::core::blast_hits::subject_besthit(frame_hits, false, |hit| {
        crate::core::blast_hits::BestHitExtent {
            context: i64::from(hit.ctx_idx),
            query_frame: i32::from(hit.q_frame),
            offset: i64::from(hit.q_aa_start),
            end: i64::from(hit.q_end_int4()),
            query_length: 0,
        }
    });
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3037-3051 (`Blast_HSPListAdjustOffsets`)
// For every TBLASTX chunk list; the subject's offset, end and seed are NCBI's `Int4`
// (`UngappedHit::s_start_int4`, `s_end_int4`, `s_seed_int4`): an HSP or a seed before its chunk
// (masks stored out of order) starts below 0, and an HSP of the one-hit word finder of length 0 or
// less may end before its start, below 0 in a frame's first chunk. Such a start lies in the frame
// (a range starts at the frame's first letter or after it), so the moved start is 0 or more.
fn adjust_tblastx_chunk_subject_offsets(hits: &mut [UngappedHit], offset: usize) {
    if offset == 0 {
        return;
    }
    let offset_int4 = i32::try_from(offset)
        .expect("a subject chunk's offset fits i32 (NCBI's BlastSeg offset is Int4)");
    let moved = |chunk_offset: i32, what: &str| -> i32 {
        chunk_offset
            .checked_add(offset_int4)
            .unwrap_or_else(|| panic!("{what}: a subject offset fits i32 (NCBI's Int4)"))
    };
    for hit in hits {
        let q_end = hit.q_end_int4();
        let s_end = moved(hit.s_end_int4(), "s_aa_end");
        hit.s_aa_start = u32::try_from(moved(hit.s_start_int4(), "s_aa_start"))
            .expect("s_aa_start: an HSP starts in its frame (NCBI's BlastSeg offset)");
        hit.set_int4_ends(q_end, s_end);
        hit.s_seed_off = u32::try_from(moved(hit.s_seed_int4(), "s_seed_off"))
            .expect("s_seed_off: a seed lies in its frame (NCBI's gapped_start)");
    }
}

/// Where a chunk merge of the one-hit word finder (`-window_size 0`) reports a refusal: the
/// search's failure and the subject's ordinal id (`frame_chunks::ChunkSums::one_hit`).
#[derive(Clone, Copy)]
struct OneHitChunkMerge<'a> {
    failure: &'a crate::algorithm::tblastx::subjects::XReadFailure,
    s_idx: usize,
}

impl OneHitChunkMerge<'_> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:581-584
    // A chunk merge of a frame of the one-hit word finder whose result NCBI leaves undefined
    // (`merge_tblastx_subject_chunk_hits`) fails the search (PD-LOSAT-NCBI-DEFECTS 1.21, rule 3),
    // as the earliest such subject's error (`XReadFailure`).
    /// Refuses the search of this subject.
    #[cold]
    #[inline(never)]
    fn refuse(self, _undefined: ChunkMergeUndefined) {
        // After the reports of the batches before (`ErrorAfterBatchReports`).
        self.failure.fail(
            self.s_idx,
            anyhow::Error::new(ErrorAfterBatchReports(anyhow::anyhow!(
                "TBLASTX -window_size 0: a one-hit ungapped alignment whose length is 0 or less meets another in the overlap of two 5,000,000-letter parts of a subject frame, and their lengths add up to 0 (NCBI BLAST+ divides by that sum when it merges them; its result is undefined), which LOSAT does not support"
            ))),
        );
    }
}

/// `settings` are the NCBI application settings that LOSAT reproduces
/// (`ncbi_environment::ApplicationSettings`): the input readers use them once
/// the program reads with `fasta_reader` (steps S3-S8 of the port plan).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:67-72
pub fn run(
    mut args: TblastxArgs,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    args.id_scope = args.parse_deflines.then(Default::default);
    // NCBI's score block looks for the file of a `-matrix` that its tables do not have when it
    // sets up a batch (`TblastxArgs::find_matrix_file`); the run refuses it there, as the name
    // that NCBI aborts on (`run_in_pool`).
    args.find_matrix_file(&settings.matrix_search)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:106-111
    // NCBI's option handlers read the subjects first, then open the query and the output,
    // then process the formatting options (blast_args.cpp:3631-3636). The thread count is
    // checked before the inputs are read, as before. LOSAT's own debug and timing output
    // (LOSAT_WASI_THREADS_DEBUG, LOSAT_TIMING, LOSAT_DEBUG_SCAN_SOFF) follows the reads and
    // does not count their time.
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
    // `ArchiveFormatRequested` parses `-outfmt` (blast_args.cpp:2745-2748) before the
    // option handlers run; a format that NCBI runs and LOSAT does not write is rejected
    // there too. `GetExportSearchStrategyStream` comes first: it creates the search strategy
    // file (`strategy::open_strategy_file`), which an error after it leaves empty; with it a
    // value that does not parse is left to the formatting handler
    // (`app::early_format_choice`).
    let strategy_file = crate::report::archive::strategy::open_strategy_file(
        args.export_search_strategy.as_deref(),
        args.out.as_deref(),
        "TBLASTX",
    )?;
    let choice =
        crate::blastinput::app::early_format_choice(&args.outfmt, strategy_file.is_some())?;
    let format = choice.report_format(
        crate::blastinput::app::Program::Tblastx,
        args.out.as_deref(),
    )?;
    // The BLAST archive (outfmt 11): its request, before anything is read or written.
    // The error archives follow the `-outfmt` text (`writes_error_archive`).
    // The run's options handle, which its archive and its search strategy share.
    let options = crate::report::archive::OptionsHandle::default();
    args.archive = super::super::archive::handle_for(
        &args,
        format == Some(crate::blastinput::app::ReportFormat::ArchiveFormat),
        &options,
    )
    .map(|handle| handle.with_outfmt_text(&args.outfmt));
    // The search strategy's request, with the archive's options handle.
    args.strategy = strategy_file.map(|file| {
        crate::report::archive::strategy::StrategyExport::for_args(
            &args,
            file,
            choice.number(),
            &options,
        )
    });
    let archive = args.archive.clone();
    crate::utils::threading::validate_threads(args.num_threads)?;
    use crate::blastinput::input_files;
    // The diagnostics stream: with the archive it also saves what it writes as the messages of
    // the next archive.
    let mut stderr = crate::cli::CliStderr::process();
    let mut archive_tee;
    let diagnostics: &mut (dyn std::io::Write + Send) = match &archive {
        Some(handle) => {
            archive_tee = crate::report::archive::DiagTee::new(&mut stderr, handle.clone());
            &mut archive_tee
        }
        None => &mut stderr,
    };
    // The taxid or seqid list (line DI-8): NCBI's handler of the database arguments
    // (`s_GetTaxIDList`, `CSeqDBFileGiList`) comes before the other option handlers
    // (`m_BlastDbArgs` is among the first of `CBlastpAppArgs`' `m_Args`, blastp_args.cpp:61-65)
    // and before `CStdCmdLineArgs`, which opens the query and the output: its warning before
    // theirs (probe F04 of line DI-8), and a list error leaves no `-out`, only the error archive
    // that NCBI's `PrintErrorArchive` writes for outfmt 11 (probe E of line DI-8).
    let list = match args.database_name() {
        Some(_) => match args
            .db_filter
            .resolve(&settings, "tblastx", &mut *diagnostics)
        {
            Ok(list) => list,
            Err(error) => {
                write_run_error_archive(&args, archive.as_ref(), &error);
                return Err(error);
            }
        },
        None => None,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3456-3481
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
    let out_file = open_out_file(&args)?;
    // Outfmt 0 to 4 (the query-anchored formats write through the same formatter).
    let pairwise = format
        .and_then(crate::report::pairwise::PairwiseView::of)
        .is_some();
    // LOSAT's reports read the format number that NCBI reads: the number and the custom fields,
    // or the value as given when it has a delimiter, which `normalized` drops and the custom
    // tabular report reads (`report::output_format`, which parses it again as NCBI does).
    // A value not read yet (`EarlyFormat::Deferred`) fails at the options' handlers (`search_cli`).
    let outfmt = match choice.parsed() {
        Some(parsed) if parsed.delimiter.is_empty() => parsed.normalized(),
        _ => args.outfmt.clone(),
    };
    // The organism report's `-out` file is written as NCBI's buffer writes it.
    let tax_report = format == Some(crate::blastinput::app::ReportFormat::TaxFormat);
    let mut stream = crate::cli::ReportStream::of_cli(out_file, tax_report, args.strategy.as_ref());
    let strategy = args.strategy.clone();
    let result = {
        let mut outputs =
            ReportOutputs::single(&outfmt, OutputSink::Writer(&mut stream), &mut *diagnostics);
        search_cli(
            args,
            list,
            query_file,
            &subjects,
            &subject_placements,
            x_stores.as_deref(),
            &mut outputs,
            settings,
        )
    };
    // A post of the search strategy left unwritten is LOSAT's internal error (`end_run`).
    let result = crate::report::archive::strategy::end_run(strategy.as_ref(), result);
    // NCBI writes an error archive after the archives of the batches before when the
    // application stops with one of its own errors, and when the run ends with messages that
    // no archive carried (`PrintErrorArchive` after its `CATCH_ALL`, whatever the status).
    if let Some(handle) = &archive {
        match &result {
            Err(error) => {
                if let Some(lines) = super::super::archive::error_lines(error, false) {
                    let _ = crate::report::archive::write_error_archive_after_post(
                        handle,
                        error,
                        &lines,
                        &mut stream,
                    );
                }
            }
            Ok(()) => {
                if handle.has_messages() {
                    let _ = crate::report::archive::write_error_archive(handle, &[], &mut stream);
                }
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:252-255
    // The outfmt 0 formatter's stream throws when a write fails (blast_format.cpp:118-119;
    // `cli::ReportStream::finish`).
    stream.finish(result, pairwise || tax_report)
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
    // `-outfmt "11 qseqid"`: NCBI's `PrintErrorArchive` stops at the number before it opens `-out`.
    if !handle.error_archive_enabled() {
        return;
    }
    let Some(lines) = super::super::archive::error_lines(error, true) else {
        return;
    };
    if let Ok(out_file) = open_out_file(args) {
        let mut report = crate::cli::ReportStream::of_cli(out_file, false, None);
        let _ = crate::report::archive::write_error_archive_after_post(
            handle,
            error,
            &lines,
            &mut report,
        );
        let _ = std::io::Write::flush(&mut report);
    }
}

/// The `-out` file, created (`-`: standard output).
fn open_out_file(args: &TblastxArgs) -> Result<Option<std::io::BufWriter<std::fs::File>>> {
    args.out
        .as_deref()
        .filter(|path| path.as_os_str() != "-")
        .map(|path| crate::cli::create_out_file(path, Some("TBLASTX")).map(std::io::BufWriter::new))
        .transpose()
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2553-2557
    // The first handler opens and reads the subjects (an empty subject set fails there),
    // as BLASTN's `run`.
    let Some(subject_path) = args.subject.clone() else {
        return Err(crate::blastinput::app::missing_subject_error());
    };
    input_files::check_utf8_file_name(&subject_path, "subject", "TBLASTX")?;
    let subject_file = input_files::open_input(&subject_path, "subject", "TBLASTX")?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2537-2545
    // The subject range is read after the subject file is opened, before it is read.
    let subject_range = parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:242-246
    // NCBI reads the records one at a time, writes each one's messages as it reads it, and
    // checks each one's range after reading it: a range that starts past the end of a
    // record stops the reading there (`read_subjects`). Records without residues (and
    // intervals without letters) are read without a message and reported when the search
    // is set up (`search`).
    let mut subject_source = FastaInputSource::from_file(
        subject_file,
        ReaderConfig::subject("TBLASTX", false, settings.data_loaders)
            .with_defline(args.archive.is_some() || args.strategy.is_some())
            .with_parse_deflines(args.parse_deflines),
    );
    subject_source.set_id_scope(args.id_scope.clone());
    // Diagnostic (LOSAT_TIMING): the subject FASTA for the common line (utils/xtiming.rs); with
    // LOSAT_X_TBXSTORE the packing of each record is in this interval too. No NCBI counterpart.
    let t_read_fasta = crate::utils::xtiming::clock();
    // EXPERIMENT (LOSAT_X_TBXSTORE): each record's letters packed as it is read
    // (`x_read_packed_subjects`, cut to `-subject_loc` first); the search reads the stores as
    // a database's (`Subjects::Stored`), the reports by record. A set with a byte without an
    // ncbi4na code (the FASTA reader gives none: it reads `U` as `T`) or a record over `u32::MAX`
    // bases keeps every record's letters. The archive of `-subject` writes every subject's whole
    // letters (and its deflines) in every batch, which the packed stores do not keep (they hold
    // the interval of a `-subject_loc` and no letters), so a run that writes it (or the search
    // strategy's requests, which hold the same Bioseqs) reads the records as the switch's `=0`
    // side does (as TBLASTN's `x_tbn_store`). So does `-lcase_masking`: the
    // scan reads the case of the letters (`x_soft_subject_calls`).
    let mut warn = |message: &[u8]| std::io::Write::write_all(&mut *diagnostics, message);
    let (subjects, subject_placements, stores) = if super::super::subjects::x_tbx_store()
        && args.archive.is_none()
        && args.strategy.is_none()
        && !args.lcase_masking
    {
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
fn check_subjects_not_empty(subjects: &[FastaRecord]) -> Result<()> {
    if subjects.is_empty() {
        return Err(crate::cli::NativeError {
            exit: 3,
            message: format!("{ENGINE_ERROR} Empty CBlastQueryVector\n").into(),
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2975-2977
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1303-1310,1358-1364
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1518-1523
/// TBLASTX has no cutoff score option, so an `-evalue` of 0 (or one that reads as 0) fails.
fn check_ncbi_options(
    args: &TblastxArgs,
    choices: &[crate::blastinput::app::EarlyFormat],
    diagnostics: &mut dyn std::io::Write,
) -> Result<(
    Option<SequenceRange>,
    ReportCounts,
    Option<std::sync::Arc<crate::report::html_ids::HtmlLinks>>,
)> {
    use crate::blastinput::app::{formatting_handler_check, options_error, Program};
    args.seg_spec()?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1975-1991 (the
    // query options handler reads the strand before the range)
    // The command line's parser has checked the value, a caller of `run_local` has not: a value
    // other than `both`, `plus` and `minus` is the parser's argument error.
    crate::algorithm::blastn::args::QueryStrand::parse(&args.strand)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1995-1999
    // The query options handler comes after the filtering handler and before the
    // formatting handler (tblastx_args.cpp:44-125).
    let query_range = crate::blastinput::seq_range::parse_optional_range(
        args.query_loc.as_deref(),
        crate::blastinput::seq_range::RangeRole::Query,
    )?;
    // Each requested format's (the library's `-outfmt` values are its outputs', not `args.outfmt`).
    let mut outfmts = Vec::with_capacity(choices.len());
    for choice in choices {
        let choice = choice.at_formatting_handler(&args.outfmt)?;
        formatting_handler_check(&choice, Program::Tblastx)?;
        outfmts.push(choice.number);
    }
    // The hit list size and the counts of the pairwise report
    // (`blastinput::report_counts`), with the handler's warnings.
    let html = crate::blastinput::app::html_links(
        args.html,
        Program::Tblastx,
        args.database_name(),
        &outfmts,
    )?;
    let (counts, warnings) = crate::blastinput::report_counts::run_report_counts_with(
        "tblastx",
        &outfmts,
        args.max_target_seqs,
        args.num_descriptions,
        args.num_alignments,
        &args.formatting,
    )?;
    diagnostics.write_all(&warnings)?;
    // (`threshold <= 0` is false for NaN, as in C.)
    if args.word_threshold() <= 0.0 {
        return Err(options_error("Non-zero threshold required"));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1312-1318
    // The command line's `-word_size` is at least 2 (`protein_word_size`); a caller of
    // `run_local` sets the field itself.
    if args.word_size() <= 0 {
        return Err(options_error("Word-size must be greater than zero"));
    }
    if args.word_size() > 4 {
        return Err(options_error(
            "Word-size must be less than 6 for protein comparison",
        ));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:608-618 (the initial
    // word options' check, after the lookup table's and before the hit saving options')
    // (`<= 0.0` is false for NaN, as in C.)
    if args.x_drop_ungapped_bits() <= 0.0 {
        return Err(options_error("x_dropoff must be greater than zero"));
    }
    // The hit saving options' check of the hit list size comes before the e-value's.
    counts.validate()?;
    if args.evalue <= 0.0 {
        return Err(options_error(
            "expect value or cutoff score must be greater than zero",
        ));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1525-1534 (after the
    // e-value's check, before the culling limit's)
    // `-max_intron_length` (blast_args.cpp:949-958 sets it, 0 when omitted); the culling limit is
    // at least 0 (`nonnegative_ncbi_integer`).
    if args.max_intron_length.unwrap_or(0) != 0 {
        return Err(options_error(
            "Uneven gap linking of HSPs is allowed for blastx, tblastn, and psitblastn only",
        ));
    }
    Ok((query_range, counts, html))
}

/// The `-outfmt` of the command line (`check_ncbi_options`), as `app::early_format_choice` reads
/// it: where the argv names a search strategy file a value that does not parse is
/// `check_ncbi_options`' to reject, at its formatting handler (the command line's run has made the
/// file then; the checks of the options alone, `check_options`, read the argv as that run).
fn cli_outfmts(args: &TblastxArgs) -> Result<[crate::blastinput::app::EarlyFormat; 1]> {
    Ok([crate::blastinput::app::early_format_choice(
        &args.outfmt,
        args.export_search_strategy.is_some(),
    )?])
}

/// The requested formats as NCBI parses `-outfmt` (`check_ncbi_options`), as the command line's
/// (`cli_outfmts`): where the argv names a search strategy file a value that does not parse is
/// left to the formatting handler, which reads the arguments' `-outfmt` again (LOSAT Web gives
/// such a value as the run's one format and as the arguments' `-outfmt`).
fn output_choices(
    outputs: &ReportOutputs<'_>,
    args: &TblastxArgs,
) -> Result<Vec<crate::blastinput::app::EarlyFormat>> {
    outputs
        .formats
        .iter()
        .map(|format| {
            crate::blastinput::app::early_format_choice(
                format.outfmt,
                args.export_search_strategy.is_some(),
            )
        })
        .collect()
}

/// The checks of the options alone, without the inputs (the `validate` of web ABI v2):
/// NCBI's processing of the options and LOSAT's limits.
pub fn check_options(args: &TblastxArgs) -> Result<()> {
    // `-db` is checked as the command line checks it: the database is opened by the run
    // (`run_local_database`, which LOSAT Web's host gives the run's files; line WQA D12a).
    check_ncbi_options(args, &cli_outfmts(args)?, &mut std::io::sink())?;
    // The host validates the argv it runs, with its `-num_threads` (docs/web/abi_v2.md).
    crate::utils::threading::validate_threads(args.num_threads)
}

/// What `search_cli` and `run_local_database` do after the options are processed and before the
/// query: the opening of the database of `-db` (`list`: its taxid or seqid list, `None` for a
/// search without `-db`), and what the requested reports ask of it.
fn open_search_database(
    args: &mut TblastxArgs,
    list: Option<&crate::seqdb::DbList>,
    settings: &crate::blastinput::ncbi_environment::ApplicationSettings,
    outputs: &mut ReportOutputs<'_>,
) -> Result<Option<super::super::database::TblastxDatabase>> {
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:113-119
    // `-db`: the database opens after the options are processed, before `Query is Empty!`
    // (`TblastxDatabase::open`).
    let database = match args.database_name() {
        Some(name) => {
            let database = super::super::database::TblastxDatabase::open_restricted(
                name,
                settings,
                list,
                args.db_mask.option().as_ref(),
                outputs.diagnostics,
            )?
            .showing_gis(args.show_gis)
            .after_calls(args.first_query_ordinal > 0);
            crate::algorithm::common::database::note_archive_database_filter(
                args.archive.as_ref(),
                args.strategy.as_ref(),
                &args.db_filter,
                database.list_taxids(),
                &database.statistics(),
                true,
                args.searchsp.is_some_and(|space| space > 0),
                args.dbsize.is_some_and(|size| size > 0),
                "TBLASTX",
            )?;
            crate::algorithm::common::database::check_formats_report_sums(
                || database.pairwise_database(),
                outputs
                    .formats
                    .iter_mut()
                    .map(|format| (format.outfmt, &mut format.refusal)),
                "TBLASTX",
            )?;
            Some(database)
        }
        None => None,
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:126-133
    // A tabular field list that looks every taxid up (the plural names of the taxonomy database,
    // the organism report) asks whether `taxdb.bti` is on the search path: the plural names'
    // formatters warn without it (`tabular`), as BLASTN and TBLASTN. A `-db` search asks its
    // database, a `-subject` search the search path of NCBI's databases (`subject_taxonomy`).
    match database.as_ref() {
        Some(database) => {
            if outputs.formats.iter().any(|format| {
                matches!(
                    report::output_format(format.outfmt),
                    report::TblastxOutputFormat::CustomTabular(spec) if spec.looks_up_every_taxid()
                )
            }) {
                args.tabular_taxdb = Some(database.taxdb_found());
            }
        }
        None => subject_taxonomy(
            args,
            Some(&crate::algorithm::common::database::TaxdbSetting::of(
                settings,
            )),
            outputs,
        ),
    }
    // LOSAT's refusal of a database name for the Seq-annot formats, before the search
    // (`report::serial::seqalign::check_database`).
    if let Some(database) = database.as_ref() {
        for format in &outputs.formats {
            if let report::TblastxOutputFormat::Seqalign(serial) =
                report::output_format(format.outfmt)
            {
                crate::report::serial::seqalign::check_database(
                    serial,
                    database.name(),
                    crate::blastinput::app::Program::Tblastx,
                )?;
            }
        }
    }
    Ok(database)
}

/// What the reports of `outputs` of a `-subject` search ask of the taxonomy database, through
/// `setting` (`common::database::subject_taxonomy`): asked by the command line
/// (`open_search_database`) and by a search of records (`run_local`) at the same point, after the
/// options are processed.
fn subject_taxonomy(
    args: &mut TblastxArgs,
    setting: Option<&crate::algorithm::common::database::TaxdbSetting>,
    outputs: &ReportOutputs<'_>,
) {
    let taxonomy = crate::algorithm::common::database::subject_taxonomy(
        setting,
        outputs
            .formats
            .iter()
            .filter_map(|format| match report::output_format(format.outfmt) {
                report::TblastxOutputFormat::CustomTabular(spec) => Some(spec),
                _ => None,
            }),
    );
    args.subject_taxdb = taxonomy.taxdb;
    args.tabular_taxdb = taxonomy.found;
}

/// The part of `run` after the output is opened: the options, `Query is Empty!`, the
/// reading of the queries and the search.
fn search_cli(
    mut args: TblastxArgs,
    list: Option<crate::seqdb::DbList>,
    query_file: crate::blastinput::input_files::InputFile,
    subjects: &[FastaRecord],
    subject_placements: &Placements,
    x_stores: Option<&[crate::sequence::packed_subject::PackedNucleotide]>,
    outputs: &mut ReportOutputs<'_>,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    let (query_range, counts, html_links) =
        check_ncbi_options(&args, &cli_outfmts(&args)?, outputs.diagnostics)?;
    args.report_counts = counts;
    args.html_links = html_links;
    let database = open_search_database(&mut args, list.as_ref(), &settings, outputs)?;
    let source = match &database {
        Some(database) => SubjectSource::Database(database),
        None => SubjectSource::Fasta {
            records: subjects,
            placements: subject_placements,
            x_stores,
        },
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:125-135
    //
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-860
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
    let mut query_source = FastaInputSource::from_file(
        query_file,
        ReaderConfig::query("TBLASTX", false, settings.data_loaders)
            .with_defline(args.archive.is_some() || args.strategy.is_some())
            .with_parse_deflines(args.parse_deflines)
            .with_strand(&args.strand),
    );
    query_source.set_id_scope(args.id_scope.clone());
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
    search(args, &records, reading, query_range, source, outputs)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:429-467
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
    mut args: TblastxArgs,
    query_records: &[FastaRecord],
    query_end: &crate::api::QueryEnd,
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // The command line checks `-strand` at parsing (clap), before anything else; a caller of
    // `run_local` has not (blastn's `run_local` does the same).
    crate::algorithm::blastn::args::QueryStrand::parse(&args.strand)?;
    // `run_local` searches subject records; a `-db` search runs through `run_local_database`.
    if args.db.is_some() {
        anyhow::bail!(crate::api::local_blast::RUN_LOCAL_DB_REFUSAL);
    }
    // NCBI parses -outfmt before its option handlers read the subjects (BLASTN's
    // `run_local`), with LOSAT's refusals; `check_ncbi_options` and `search` check the formats
    // again, after `-out` is opened, as on the command line.
    report_formats(outputs, &args)?;
    // The BLAST archive (outfmt 11) among the formats: its request and a diagnostics stream that
    // saves the messages. The records' deflines are
    // those the caller's reader kept; records without one get `>` and their title.
    let archive_requested = outputs
        .formats
        .iter()
        .any(|format| report::output_format(format.outfmt) == report::TblastxOutputFormat::Archive);
    // The archive shares the options handle of the run's search strategy (LOSAT Web's
    // `-export_search_strategy`, `web_adapter::give_strategy`), as on the command line (`run`).
    if archive_requested && args.archive.is_none() {
        args.archive = super::super::archive::handle_for(
            &args,
            true,
            &crate::report::archive::strategy::options_of(args.strategy.as_ref()),
        )
        .map(|handle| outputs.with_archive_outfmt(handle));
    }
    let archive = args.archive.clone();
    let strategy = args.strategy.clone();
    let no_query = query_records.is_empty() && query_end.is_input();
    outputs.run_records(strategy.as_ref(), archive.as_ref(), no_query, |outputs| {
        run_local_formats(args, query_records, query_end, subject_records, outputs)
    })
}

/// `run_local` from the thread count on, with its archive and its search strategy set up.
fn run_local_formats(
    mut args: TblastxArgs,
    query_records: &[FastaRecord],
    query_end: &crate::api::QueryEnd,
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // The thread count is checked before the inputs are read, as on the command line (`run`).
    crate::utils::threading::validate_threads(args.num_threads)?;
    // The subject range and its record checks come where NCBI reads the subjects (`run`):
    // the messages of each record read, up to a record whose range starts past its end.
    let subject_range = parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject)?;
    for record in &subject_records[..subjects_read(subject_records, subject_range.as_ref())] {
        outputs.diagnostics.write_all(&record.warnings)?;
    }
    let ranged_subjects = cut_subjects(subject_records, subject_range.as_ref())?;
    check_subjects_not_empty(subject_records)?;
    // `-out` (`run_to_file`) is made where the command line opens it (`run`).
    outputs.open_out()?;
    let (query_range, counts, html_links) =
        check_ncbi_options(&args, &output_choices(outputs, &args)?, outputs.diagnostics)?;
    args.report_counts = counts;
    args.html_links = html_links;
    // The taxonomy database of the reports, where the command line asks (`open_search_database`).
    let setting = args.taxdb_setting.take();
    subject_taxonomy(&mut args, setting.as_ref(), outputs);
    let whole_subjects = Placements::default();
    let (searched_subjects, subject_placements) = match &ranged_subjects {
        Some((cut, placements)) => (cut.as_slice(), placements),
        None => (subject_records, &whole_subjects),
    };
    search(
        args,
        query_records,
        query_end.reading(),
        query_range,
        SubjectSource::Fasta {
            records: searched_subjects,
            placements: subject_placements,
            x_stores: None,
        },
        outputs,
    )
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:113-119
/// Runs one TBLASTX search of the records `query_records` against the database of `args` and
/// writes every requested output format from the same result (the library API's entry; the
/// command line's is `run`). `settings` are the NCBI application settings, which find the
/// database. The steps are `run`'s and `search_cli`'s, in their order; what differs is what they
/// read: the records are read already, and the outputs are the caller's. An archive (outfmt 11)
/// ends as the command line's does (`ReportOutputs::finish_archive`).
pub fn run_local_database(
    mut args: TblastxArgs,
    query_records: &[FastaRecord],
    query_end: &crate::api::QueryEnd,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // The command line checks `-strand` at parsing (clap), before anything else; a caller of
    // `run_local` has not (blastn's `run_local` does the same).
    crate::algorithm::blastn::args::QueryStrand::parse(&args.strand)?;
    // The file of an unknown `-matrix`, then, as the command line's `run` looks it up.
    args.find_matrix_file(&settings.matrix_search)?;
    // NCBI parses -outfmt before its option handlers run, with LOSAT's refusals;
    // `check_ncbi_options` and `search` check the formats again, after `-out` is opened.
    report_formats(outputs, &args)?;
    let archive_requested = outputs
        .formats
        .iter()
        .any(|format| report::output_format(format.outfmt) == report::TblastxOutputFormat::Archive);
    // The archive shares the options handle of the run's search strategy (LOSAT Web's
    // `-export_search_strategy`, `web_adapter::give_strategy`), as on the command line (`run`).
    args.archive = super::super::archive::handle_for(
        &args,
        archive_requested,
        &crate::report::archive::strategy::options_of(args.strategy.as_ref()),
    )
    .map(|handle| outputs.with_archive_outfmt(handle));
    crate::utils::threading::validate_threads(args.num_threads)?;
    let archive = args.archive.clone();
    let strategy = args.strategy.clone();
    let no_query = query_records.is_empty() && query_end.is_input();
    outputs.run_records(strategy.as_ref(), archive.as_ref(), no_query, |outputs| {
        search_database(args, query_records, query_end, settings, outputs)
    })
}

/// `run_local_database` from the taxid or seqid list on.
fn search_database(
    mut args: TblastxArgs,
    query_records: &[FastaRecord],
    query_end: &crate::api::QueryEnd,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // The handler of the database arguments reads the list before the other handlers (`run`).
    let list = args
        .db_filter
        .resolve(&settings, "tblastx", &mut *outputs.diagnostics)?;
    // `-out` (`run_to_file`) is made where the command line opens it (`run`).
    outputs.open_out()?;
    let (query_range, counts, html_links) =
        check_ncbi_options(&args, &output_choices(outputs, &args)?, outputs.diagnostics)?;
    args.report_counts = counts;
    args.html_links = html_links;
    let database = open_search_database(&mut args, list.as_ref(), &settings, outputs)?
        .ok_or_else(crate::blastinput::app::missing_subject_error)?;
    search(
        args,
        query_records,
        query_end.reading(),
        query_range,
        SubjectSource::Database(&database),
        outputs,
    )
}

/// Each requested format as NCBI parses `-outfmt` before its option handlers run, and LOSAT's
/// rejection of a format that its TBLASTX does not write (`ReportOutputs::early_choices`): the
/// command line's checks before it opens `-out` (`run`). NCBI's formatting handler comes after
/// (`check_ncbi_options`, `check_report_formats`), and so does the error of a value left to it
/// with a search strategy file.
fn report_formats(outputs: &ReportOutputs<'_>, args: &TblastxArgs) -> Result<()> {
    outputs
        .early_choices(
            crate::blastinput::app::Program::Tblastx,
            args.out.as_deref(),
            args.export_search_strategy.is_some(),
        )
        .map(drop)
}

/// Each requested format as NCBI reads `-outfmt`: NCBI's errors, and a rejection of a
/// format that LOSAT's TBLASTX does not write.
///
/// `out` is the run's `-out` (a format with the caller's receiver of side files takes that
/// receiver's, `FormatOutput::side_out`): outfmt 13 and 14 write side files next to it, and fail
/// on standard output (`xinclude_check_for`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2887
/// The formatting handler checks the custom delimiter of each format that LOSAT writes too
/// (`formatting_handler_check`; the command line's comes before, in `check_ncbi_options`).
fn check_report_formats(outputs: &ReportOutputs<'_>, out: Option<&std::path::Path>) -> Result<()> {
    use crate::blastinput::app;
    for format in &outputs.formats {
        let choice = app::parse_formatting_string(format.outfmt)?;
        app::report_format(&choice, app::Program::Tblastx, format.side_out(out))?;
        app::formatting_handler_check(&choice, app::Program::Tblastx)?;
    }
    // The XInclude formats without a file, the one check of every program (`ReportOutputs`).
    outputs.xinclude_check(out)
}

/// The search of `run_local` and of the CLI, after NCBI's checks of the options and the
/// reading of the subjects: `Query is Empty!` for records without a query, the query batch
/// size, the subjects without letters, and the query batches.
fn search(
    args: TblastxArgs,
    query_records: &[FastaRecord],
    reading: QueryReading,
    query_range: Option<SequenceRange>,
    source: SubjectSource<'_>,
    outputs: &mut ReportOutputs<'_>,
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
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:136-137
    // The batch size is read before the formatter is made and before the queries, which
    // NCBI reads one batch at a time after the outfmt 0 prolog (`run_in_pool`).
    let batch_size = crate::blastinput::app::query_batch_size(10002)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:129-136
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:773-789
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
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:170-173
    // The XInclude formats written to standard output fail after the formatter is made
    // (`check_report_formats`: `xinclude_check`), as the other format checks of each
    // requested format.
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:140-153 (cited at
    // `blastinput::db_mask::DbMaskOption::resolve`)
    // The formatter resolves the database masking (`-db_soft_mask`, `-db_hard_mask`; series NP,
    // M4b PR 2) before its other work: its warning (one message of the archive) or its database
    // error; the archive and the search strategy write the result's program options.
    let db_mask = match (&source, args.db_mask.option()) {
        (SubjectSource::Database(database), Some(option)) => database.resolve_mask_option(
            &option,
            args.archive.as_ref(),
            args.strategy.as_ref(),
            outputs.diagnostics,
        )?,
        _ => None,
    };
    check_report_formats(outputs, args.out.as_deref())?;
    // LOSAT's refusals of the masked search (`-db_hard_mask`) come where NCBI starts the search,
    // after its checks.
    if let (SubjectSource::Database(database), Some(mask)) = (&source, db_mask.as_ref()) {
        database.check_mask_limits(mask, usize::try_from(args.word_size())?)?;
    }
    crate::blastinput::app::check_search_environment(
        "TBLASTX",
        !matches!(source, SubjectSource::Database(_)),
    )?;
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
    // `-db`: the database checked once before the first batch (`TblastxDatabase::x_stream`), as a
    // copy of its sequences would be read, and its sequences read where each batch searches them
    // (`Subjects::Stream`); NCBI reads them first in the preliminary search of the first batch
    // with a valid context (a batch without one is not searched, `CheckInternalData`). A
    // database without sequences (BLASTN's `Database memory map file error`) or a sequence file
    // that cannot be read fails that batch, after the reports of the batches before
    // (`run_in_pool`), and LOSAT's refusal of a damaged file stops the run there; a run without
    // such a batch never searches the file. A read that fails during a batch's search (a database
    // changed after the check) fails that batch in the same way (`XReadFailure`).
    let whole_subjects = Placements::default();
    let (stream, stores_error) = match &source {
        SubjectSource::Database(database) => {
            let opened = if database.first_oid_read_fails() {
                Err(crate::algorithm::common::database::PreliminaryReadError(
                    "Database memory map file error",
                )
                .into())
            } else {
                database.x_stream(args.num_threads, db_mask.as_ref())
            };
            match opened {
                Ok(stream) => (Some(stream), None),
                Err(error) => (None, Some(error)),
            }
        }
        SubjectSource::Fasta { .. } => (None, None),
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
                    link_total: stores.iter().map(|store| store.len() as i64).sum(),
                    average: None,
                    dbsize: args.dbsize,
                    searchsp: args.searchsp,
                    statistics: None,
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
                numbers: SubjectSetNumbers::of_records(records, args.dbsize)
                    .with_searchsp(args.searchsp),
                shown: ShownSubjects::Fasta {
                    records,
                    x_stored: None,
                },
            },
            placements,
        ),
        SubjectSource::Database(database) => (
            RunSubjects {
                subjects: match &stream {
                    Some(stream) => Subjects::Stream(stream),
                    None => Subjects::Unreadable,
                },
                numbers: database.numbers(args.dbsize).with_searchsp(args.searchsp),
                shown: ShownSubjects::Database {
                    database,
                    pairwise: database.pairwise_database_masked(db_mask.as_ref()),
                    stream: stream.as_ref(),
                },
            },
            &whole_subjects,
        ),
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-188
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
        )
    })
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
    /// A list kept no sequence of the database: the messages of the queries of the searched
    /// batches have `EMPTY_DATABASE_MESSAGE` (`SeqDb::filtered_empty`).
    pub filtered_empty: bool,
    /// The posts written before each query's report (`QueryWarnings`): the reader's messages of
    /// the query's batch and NCBI's posts of the bytes replaced in its search strategy request
    /// (before its first query), and the query's own warnings.
    pub warnings: &'a [Posts],
    /// The texts of the subjects that each query's report reads one by one (its batch's
    /// prefetch off), by query, which the first format's writer places where it reads them
    /// (`QueryWarnings::subject_read`).
    pub query_reads: &'a [crate::report::query_warnings::QueryReads],
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
    /// The OID of each ordinal id, which names a sequence without parsed Seq-ids.
    pub oids: crate::seqdb::OidMap<'a>,
}

/// The `Database:` text of the reports of a local `-subject` search.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:797-799
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

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:173-178
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:118-119
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
        if report::output_format(format.outfmt).is_pairwise()
            && prolog_before_search(&format.sink) == before_search
        {
            let flush = crate::report::pairwise::flushes_after_prolog(
                args.html_links.is_some(),
                !format.sink.is_writer(),
            );
            let mut writer = format.sink.open()?;
            match subjects {
                ShownSubjects::Fasta { records, x_stored } => {
                    crate::report::pairwise::write_tblastx_pairwise_prolog(
                        &mut writer,
                        crate::wording::REPORT_VERSION,
                        tblastx_database_name(args).as_bytes(),
                        records.len() as i64,
                        // NCBI's `Int8` sum (c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-181: `m_iTotalLength` grows by
                        // the `length` of each element, cast to `Int8`).
                        match x_stored {
                            Some((stores, _)) => {
                                stores.iter().map(|store| store.len() as i64).sum()
                            }
                            None => records.iter().map(|r| r.seq().len() as i64).sum(),
                        },
                        None,
                        args.html_links.as_deref().map(|links| {
                            let posts = crate::report::query_warnings::html_posts_sink(
                                outputs.html_posts.as_deref_mut(),
                                &mut *outputs.diagnostics,
                            );
                            (links, posts as &mut dyn std::io::Write)
                        }),
                    )?
                }
                ShownSubjects::Database { pairwise, .. } => {
                    crate::report::pairwise::write_tblastx_pairwise_prolog(
                        &mut writer,
                        crate::wording::REPORT_VERSION,
                        b"",
                        0,
                        0,
                        Some(pairwise),
                        args.html_links.as_deref().map(|links| {
                            let posts = crate::report::query_warnings::html_posts_sink(
                                outputs.html_posts.as_deref_mut(),
                                &mut *outputs.diagnostics,
                            );
                            (links, posts as &mut dyn std::io::Write)
                        }),
                    )?
                }
            }
            if flush {
                writer.flush()?;
            }
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:203-209
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1574-1577
/// Writes the final hits (already filtered by e-value) to every requested format, and
/// gives the caller the same final list (`ReportOutputs::hits`).
/// The subjects of a call's Seq-aligns (the archive's).
fn seqalign_subjects<'a>(
    run: &'a TblastxReportRun<'_>,
) -> crate::algorithm::tblastx::seqalign::SeqalignSubjects<'a> {
    match &run.database {
        Some(database) => crate::algorithm::tblastx::seqalign::SeqalignSubjects::Database {
            ordinals: &database.reported.ordinals,
            oids: database.oids,
            title: &database.pairwise.definition,
            parsed: &database.reported.parsed,
        },
        None => crate::algorithm::tblastx::seqalign::SeqalignSubjects::Records {
            records: run.subjects,
            placements: run.ranges.subject_placements,
        },
    }
}

/// What the archive reads of a call's report data.
fn archive_data<'a>(
    run: &'a TblastxReportRun<'_>,
    hits: &'a [TblastxHsp],
    engine_counts: &'a [(u32, u32)],
    subjects: &'a crate::algorithm::tblastx::seqalign::SeqalignSubjects<'a>,
) -> crate::algorithm::tblastx::archive::ArchiveData<'a> {
    crate::algorithm::tblastx::archive::ArchiveData {
        queries: run.queries,
        query_stats: run.query_stats,
        searched: run.searched,
        subject_records: match &run.database {
            Some(_) => &[],
            None => run.subjects,
        },
        hits,
        engine: engine_counts,
        subjects,
        subject_letters: run.subject_letters,
        query_placements: &run.ranges.input.placements,
        numbers: run.numbers,
        spans: run.archive_spans,
        lcase_masking: run.args.lcase_masking,
    }
}

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
    // Outfmt 0's refusal, each such output's (`refuse_shared_masks_in`; formats 1 to 4 write a
    // query's own record in its section, `query_anchored::aln_mix`, and end with NCBI's exception
    // where NCBI's merge throws, `query_anchored::write_section`).
    if args.parse_deflines {
        crate::blastinput::fasta_reader::believed::refuse_shared_masks_in(
            outputs
                .formats
                .iter_mut()
                .zip(&output_formats)
                .filter(|(_, output_format)| {
                    **output_format
                        == report::TblastxOutputFormat::Pairwise(
                            crate::report::pairwise::PairwiseView::Pairwise,
                        )
                })
                .map(|(format, _)| &mut format.refusal),
            run.queries,
            |query| {
                run.query_stats
                    .get(query)
                    .is_some_and(|stats| stats.seg_masks.iter().any(|(_, masks)| !masks.is_empty()))
            },
            args.lcase_masking,
            "TBLASTX",
        )?;
    }
    // The records of a subject set whose ids `-parse_deflines` read, by subject index.
    let records_parsed = if args.parse_deflines && run.database.is_none() {
        crate::report::db_subject::parsed_records(
            run.subjects
                .iter()
                .map(|record| (record.believed(), record.title.as_slice())),
            args.show_gis,
        )
    } else {
        Vec::new()
    };
    // LOSAT's refusals of these records in outfmt 0 to 4, each such format's (`FormatRefusal`).
    if !records_parsed.is_empty() {
        for (format, output_format) in outputs.formats.iter_mut().zip(&output_formats) {
            if output_format.is_pairwise() {
                format
                    .refusal
                    .check(crate::report::db_subject::check_pairwise_records(
                        &records_parsed,
                        false,
                        "TBLASTX",
                    ))?;
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1199-1220
    // The Seq-aligns (outfmt 8, 12) show the engine's identities and positives
    // (`report::engine_counts`), read before the display counts replace the identities.
    // The positives of the reports are the search matrix's (a name without a table stops the
    // search before any report has a hit, `run_batch_loop`); BLOSUM62's counts keep their
    // constant matrix (`Blosum62Score`).
    let matrix = match run.args.scoring_matrix() {
        TblastxMatrix::Table(matrix) => matrix,
        TblastxMatrix::Unknown | TblastxMatrix::Aborts => ScoringMatrix::Blosum62,
    };
    let (engine_counts, seqalign_identities, pairwise, pairwise_masks) =
        if matrix == ScoringMatrix::Blosum62 {
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
    let pairwise = pairwise.map(std::sync::Arc::new);
    if let (Some(hits_sink), Some(pairwise)) = (outputs.hits.as_mut(), pairwise.as_ref()) {
        hits_sink(std::sync::Arc::clone(pairwise));
    }
    // Built before any format is written, so that a report that cannot be made fails the
    // run without output (`tabular_subject_ids`).
    let subject_ids = if output_formats.iter().any(|format| format.is_tabular()) {
        match &run.database {
            Some(database) => database.reported.tabular_subject_ids(
                database.oids,
                output_formats
                    .iter()
                    .any(|format| format.tabular_reads_subject()),
            )?,
            None => tabular_subject_ids(run.subjects),
        }
    } else {
        Vec::new()
    };
    // The rows of subjects whose lookup fails (no id) where NCBI has no result
    // (`report::tabular::check_failed_lookups`), before any format is written: each tabular
    // format's refusal is that format's (`FormatRefusal::check`).
    let parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>] = match &run.database {
        Some(database) => &database.reported.parsed,
        None => &records_parsed,
    };
    // A subject set's subjects are all loaded.
    let every_loaded = crate::report::db_subject::HitLoads::default();
    let loads = run
        .database
        .as_ref()
        .map_or(&every_loaded, |database| &database.reported.loads);
    for (output, format) in outputs
        .formats
        .iter_mut()
        .zip(&output_formats)
        .filter(|(_, format)| format.is_tabular())
    {
        let spec = format.tabular_spec();
        if crate::algorithm::common::tabular::checks_failed_lookups(
            spec,
            &subject_ids,
            parsed,
            loads,
        ) {
            let checked = crate::report::tabular::check_failed_lookups(
                spec,
                (0..run.queries.len()).map(|q_idx| {
                    hits.iter()
                        .filter(move |hsp| hsp.hit.q_idx as usize == q_idx)
                        .map(|hsp| hsp.hit.s_idx as usize)
                }),
                |s_idx| {
                    crate::algorithm::common::tabular::subject_lookup_failure(
                        spec,
                        &subject_ids,
                        parsed,
                        loads,
                        run.database.is_some(),
                        false,
                        s_idx,
                    )
                },
                |s_idx| match &run.database {
                    Some(database) => crate::report::tabular::failed_database_subject(
                        database.oids.oid(database.reported.ordinals[s_idx]),
                        loads.thrown(s_idx),
                    ),
                    None => crate::report::tabular::failed_record_title(
                        s_idx,
                        run.subjects[s_idx].shown_id(),
                    ),
                },
                "TBLASTX",
            );
            output.refusal.check(checked)?;
        }
    }
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:497-498
    //
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2273
    // The outfmt 0 descriptions and headings are made from the subjects' title bytes
    // (`report::defline::generate_defline`), by subject index, as BLASTN's.
    let subject_titles: Vec<Arc<[u8]>> = if output_formats.iter().any(|format| format.is_pairwise())
    {
        run.subjects
            .iter()
            .map(|subject| Arc::from(subject.title.as_slice()))
            .collect()
    } else {
        Vec::new()
    };
    // The query-anchored formats (outfmt 1 to 4, line RG T3c): their refusals before any format
    // of the call is written (`anchored::refuse_anchored_reports`), each such format's
    // (`FormatRefusal`), and the sequences of their sections, made once.
    let anchored_view = |format: &report::TblastxOutputFormat| match format {
        report::TblastxOutputFormat::Pairwise(view)
            if *view != crate::report::pairwise::PairwiseView::Pairwise =>
        {
            Some(*view)
        }
        _ => None,
    };
    for (format, output_format) in outputs.formats.iter_mut().zip(&output_formats) {
        if let Some(view) = anchored_view(output_format) {
            format
                .refusal
                .check(super::super::anchored::refuse_anchored_reports(
                    run,
                    &hits,
                    seqalign_identities.as_deref(),
                    view.report_sort(args.formatting.report_sort()),
                    args.report_counts.num_descriptions,
                ))?;
        }
    }
    let anchored = output_formats
        .iter()
        .any(|format| anchored_view(format).is_some())
        .then(|| super::super::anchored::AnchoredSequences {
            queries: run.queries,
            query_placements: &run.ranges.input.placements,
            subjects: run.subjects,
            subject_letters: run.subject_letters,
            subject_placements: run.ranges.subject_placements,
            database: run.database,
            show_gis: args.show_gis,
            query_masks: run
                .query_stats
                .iter()
                .zip(run.queries)
                .enumerate()
                .map(|(q_idx, (stats, query))| {
                    report::query_frame_masks(
                        stats,
                        query.seq().len(),
                        run.ranges.input.placements.offset(q_idx),
                    )
                })
                .collect(),
            query_code: report::DisplayCode::new(&query_code),
            db_code: report::DisplayCode::new(&db_code),
        });
    let database = match &run.database {
        Some(database) => database.name.to_string(),
        None => tblastx_database_name(args),
    };
    // The `Query=` and `# Query:` text of each query: its title bytes.
    let query_lines: Vec<std::borrow::Cow<'_, [u8]>> =
        run.queries.iter().map(|query| query.query_line()).collect();
    let query_titles: Vec<&[u8]> = query_lines.iter().map(|line| line.as_ref()).collect();
    let unsearched: Vec<bool> = run.searched.iter().map(|&searched| !searched).collect();
    // The warnings of each query are written before its report, once, with the first format
    // (`QueryWarnings`); the taxonomy database's messages of XML2 and JSON go to the
    // diagnostics at every format.
    let mut warning_lines = Some(run.warnings);
    // The queries' result warnings that a structured first format does not write are each later
    // format's that writes them, to the caller's receiver of that format's
    // (`ReportOutputs::later_warnings`).
    // NCBI reference (598d8ae6): blast_format.cpp:1450-1452 (at `QueryWarnings::with_later`).
    let keeps_later = outputs.keeps_later_query_warnings();
    let observer = &mut outputs.observer;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-89
    // The statistics of each query for the observer, the values of the XML reports.
    if let Some(observer) = observer.as_deref_mut() {
        for (q_idx, stats) in run.query_stats.iter().enumerate() {
            observer.query_stats(
                run.queries_before + q_idx,
                super::super::blast_xml::query_stats(stats, run.numbers),
            );
        }
    }
    // Every structured format of the call is XML2 or JSON: none writes the report that a
    // `-subject` document holds back (`blast_xml::write_tblastx_xml`).
    let xml_only_xml2 = !output_formats.contains(&report::TblastxOutputFormat::Xml);
    for (format_index, (format, &output_format)) in
        outputs.formats.iter_mut().zip(&output_formats).enumerate()
    {
        // A format that failed alone writes nothing more; the next takes the warnings.
        if format.refusal.refused() {
            continue;
        }
        // The warnings that a format other than the archive writes are not the archive's
        // messages (a call with several formats, `run_local`).
        if let Some(handle) = &run.args.archive {
            handle.set_muted(output_format != report::TblastxOutputFormat::Archive);
        }
        let mut probe = observer
            .as_deref_mut()
            .map(|observer| FormatProbe::new(observer, format_index));
        // So do the warnings of the formatters of a tabular field list that shows a taxonomy
        // name (`x_CheckTaxDB`, `tabular`), and NCBI's posts of the bytes that the archive's
        // and the Seq-annots' streams replace (`write_tblastx_archives`,
        // `seqalign::write_query_to`).
        // NCBI reference (598d8ae6): blast_format.cpp:1450-1452 (at `QueryWarnings::with_later`).
        let later = format.later_warnings(keeps_later, &mut outputs.later_warnings, format_index);
        let mut warnings = (warning_lines.is_some()
            || later.is_some()
            || crate::report::query_warnings::report_posts(
                format.report_format(),
                args.html_links.is_some(),
            )
            || matches!(
                output_format,
                report::TblastxOutputFormat::CustomTabular(spec) if spec.looks_up_every_taxid()
            ))
        .then(|| {
            let first = warning_lines.take();
            let warnings = crate::report::query_warnings::QueryWarnings::of_posts(
                first.unwrap_or(&[]),
                run.args.archive.as_ref(),
                &mut *outputs.diagnostics,
            )
            .with_later(run.warnings, later)
            .for_report(format_index)
            .with_html_posts(
                outputs
                    .html_posts
                    .as_deref_mut()
                    .map(|posts| posts as &mut dyn std::io::Write),
            );
            // The first format's writer also writes the texts of the subjects that the reports
            // read one by one (`TblastxReportRun::query_reads`).
            match first {
                Some(_) => warnings.with_subject_reads(run.query_reads),
                None => warnings,
            }
        });
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:68-93
        let prolog = !prolog_before_search(&format.sink);
        // The caller's receiver of the side files of outfmt 13 and 14, if it gives one.
        let side_sink = format.sink.side_files().cloned();
        // The format's report, whose error goes to its refusal (below).
        let written = (|| -> anyhow::Result<()> {
            let mut writer = crate::cli::NcbiOstream::sink(&mut format.sink)?;
            match output_format {
                report::TblastxOutputFormat::Pairwise(view) => {
                    let queries: Vec<crate::report::pairwise::TblastxPairwiseQuery> = run
                        .queries
                        .iter()
                        .zip(run.query_stats)
                        .zip(&query_titles)
                        .enumerate()
                        .map(|(q_idx, ((query, stats), title))| {
                            crate::report::pairwise::TblastxPairwiseQuery {
                                query_name: title.to_vec(),
                                // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:742-744
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
                        version: crate::wording::REPORT_VERSION.to_string(),
                        database_name: database.clone().into(),
                        // `-db`: the `Database:` lines are the database's (`database` below).
                        // NCBI's `Int8` counts (`PairwiseDatabase`); a FASTA set's records are
                        // held in memory.
                        database_num_sequences: match &run.database {
                            Some(database) => database.pairwise.num_sequences,
                            None => run.subjects.len() as i64,
                        },
                        database_total_letters: match &run.database {
                            Some(database) => database.pairwise.total_letters,
                            None => (0..run.subjects.len())
                                .map(|s_idx| run.subject_letters.length(s_idx) as i64)
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
                            None => records_parsed.clone(),
                        },
                        subject_sources: match &run.database {
                            Some(database) => database.reported.sources.clone(),
                            None => Vec::new(),
                        },
                        sum_stats: args.sum_stats(),
                    };
                    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:166-168
                    let config = crate::report::pairwise::PairwiseConfig {
                        program: "tblastx".to_string(),
                        protein_matrix: matrix,
                        line_length: args.formatting.line_length(),
                        // `-sorthits` for every view, `-sorthsps` for outfmt 0 only (NCBI's warning
                        // for the others, `PairwiseView::report_sort`).
                        sort: view.report_sort(args.formatting.report_sort()),
                        view,
                        html: args.html_links.as_ref().map(|links| {
                            crate::report::html_ids::HtmlReport::new(
                                links,
                                run.database.into_iter().flat_map(|database| {
                                    database
                                        .reported
                                        .ordinals
                                        .iter()
                                        .map(move |&o| database.oids.oid(o))
                                }),
                            )
                        }),
                        ..crate::report::pairwise::PairwiseConfig::default()
                    };
                    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251 (quoted
                    // at `cli::NativeError::exception`)
                    // An exception of NCBI's merge in a query-anchored section ends the report after
                    // what it wrote before (`query_anchored::write_section`) with `CATCH_ALL`'s report
                    // (`NativeError::of_writer`), without the epilog or the later batches.
                    // A subject set's subjects are all loaded.
                    let every_loaded = crate::report::db_subject::HitLoads::default();
                    let mut write = |rows: &[crate::report::pairwise::PairwiseHit],
                                     masks: Option<&[crate::report::query_anchored::HspMasks]>,
                                     range: std::ops::Range<usize>|
                     -> Result<()> {
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
                            anchored.as_ref().map(|sequences| {
                                sequences as &dyn crate::report::query_anchored::ReportSequences
                            }),
                            masks,
                            run.database
                                .map_or(&every_loaded, |database| &database.reported.loads),
                        )
                        .map_err(crate::cli::NativeError::of_writer)
                    };
                    match pairwise.as_deref() {
                        Some(rows) => write(rows, pairwise_masks.as_deref(), 0..queries.len())?,
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
                            let windows = report::pairwise_windows(
                                &hits,
                                queries.len(),
                                run.x_pairwise_pool,
                                config.line_length,
                            );
                            if windows.is_empty() {
                                // No query: the prolog and the epilog alone.
                                write(&[], None, 0..0)?;
                            }
                            for range in windows {
                                let first = hits
                                    .partition_point(|hsp| (hsp.hit.q_idx as usize) < range.start);
                                let end = hits
                                    .partition_point(|hsp| (hsp.hit.q_idx as usize) < range.end);
                                let window = &hits[first..end];
                                let identities = seqalign_identities
                                    .as_deref()
                                    .map(|identities| &identities[first..end]);
                                // Under `-html` the rows' masks (`HtmlRowMasks`).
                                let mut masks = config.html.is_some().then(Vec::new);
                                // EXPERIMENT (LOSAT_X_TBXRPTPAR): the rows as an indexed map on
                                // the pool.
                                let mut rows = match run.x_pool {
                                    Some(pool) => {
                                        crate::algorithm::tblastx::x_rpt_par::pairwise_hits(
                                            pool,
                                            window,
                                            identities,
                                            &shown,
                                            masks.as_mut(),
                                        )
                                    }
                                    None => {
                                        let mut buffer = Vec::new();
                                        let mut row_masks =
                                            crate::report::query_anchored::HspMasks::default();
                                        window
                                            .iter()
                                            .enumerate()
                                            .map(|(index, hsp)| {
                                                let hit = shown.hit(
                                                    hsp,
                                                    identities.map(|identities| identities[index]),
                                                    &mut buffer,
                                                    &mut row_masks,
                                                );
                                                if let Some(table) = masks.as_mut() {
                                                    table.push(row_masks.copy_rows());
                                                }
                                                hit
                                            })
                                            .collect()
                                    }
                                };
                                // The rows count BLOSUM62's positives (`report::set_row_positives`).
                                if matrix != ScoringMatrix::Blosum62 {
                                    report::set_row_positives(&mut rows, matrix);
                                }
                                write(&rows, masks.as_deref(), range)?;
                            }
                        }
                    }
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:690-698
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
                                    oids: database.oids,
                                    title: &database.pairwise.definition,
                                    parsed: &database.reported.parsed,
                                }
                            }
                            None => {
                                crate::algorithm::tblastx::seqalign::SeqalignSubjects::Records {
                                    records: run.subjects,
                                    placements: run.ranges.subject_placements,
                                }
                            }
                        },
                        run.subject_letters,
                        warnings
                            .as_mut()
                            .expect("a Seq-annot writes its posts to the diagnostics"),
                    )?
                }
                // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:198
                // One BLAST archive per query batch (`archive`), from the Seq-aligns of outfmt 8.
                report::TblastxOutputFormat::Archive => {
                    let handle = run
                        .args
                        .archive
                        .as_ref()
                        .expect("a BLAST archive request has its handle");
                    let subjects = seqalign_subjects(run);
                    crate::algorithm::tblastx::archive::write_tblastx_archives(
                        handle,
                        &archive_data(run, &hits, &engine_counts, &subjects),
                        warnings
                            .as_mut()
                            .expect("the BLAST archive writes its posts to the diagnostics"),
                        &mut writer,
                    )?
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:708-744
                // The document of every query; the hits are the final list after `shift_to_records`
                // (record coordinates), the rows are read from the searched letters.
                report::TblastxOutputFormat::Xml => super::super::blast_xml::write_tblastx_xml(
                    &mut crate::report::blast_xml::BlastXmlWriter::new(0),
                    false,
                    xml_only_xml2,
                    &mut writer,
                    (&hits, seqalign_identities.as_deref()),
                    run,
                    warnings.as_mut(),
                )?,
                // The XML2 and JSON documents (outfmt 16 and 15, and 14 and 13 as side files next to
                // `-out` with this output as the master file): the same records, written by
                // `Xml2Writer`.
                report::TblastxOutputFormat::Xml2(kind) => {
                    super::super::blast_xml::write_tblastx_xml(
                        &mut crate::report::blast_xml2::Xml2Writer::for_run(
                            kind,
                            super::super::blast_xml::xml2_run(run),
                            run.args.out.as_deref(),
                            side_sink.as_ref(),
                        ),
                        true,
                        xml_only_xml2,
                        &mut writer,
                        (&hits, None),
                        run,
                        warnings.as_mut(),
                    )?
                }
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
                // Refused before the search (`TblastxOutputFormat::Unsupported`).
                report::TblastxOutputFormat::Unsupported => {
                    anyhow::bail!(
                        "LOSAT internal error: an output format that TBLASTX does not write"
                    )
                }
                report::TblastxOutputFormat::Tabular
                | report::TblastxOutputFormat::TabularWithComments => match run.x_pool {
                    // EXPERIMENT (LOSAT_X_TBXRPTPAR): the rows formatted on the pool; with the format
                    // observer (which flushes around every row) the serial writer, which also
                    // writes the rows of subjects whose load throws (`report::write_tabular`).
                    Some(pool) if probe.is_none() && !loads.any_throws() => {
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
                        loads,
                        &unsearched,
                        run.epilog.then_some(run.queries_before),
                        probe.as_mut(),
                        warnings.as_mut(),
                    )?,
                },
            }
            std::io::Write::flush(&mut writer)?;
            Ok(())
        })();
        // An error of the format's report fails the call, or this format alone (`FormatRefusal`):
        // the formats after it are written then.
        format.refusal.check(written)?;
    }
    if let Some(handle) = &run.args.archive {
        handle.set_muted(false);
    }
    // Without an output format the warnings are written all at once.
    if let Some(posts) = warning_lines {
        for query_posts in posts {
            query_posts.write(&mut *outputs.diagnostics, run.args.archive.as_ref())?;
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411-1585
/// NCBI makes and prints the alignments of one query's result set at a time.
fn rows_up_front(
    outputs: &ReportOutputs<'_>,
    run: &TblastxReportRun<'_>,
    output_formats: &[report::TblastxOutputFormat],
) -> bool {
    outputs.hits.is_some()
        || (output_formats.iter().any(|format| format.is_pairwise())
            && (outputs.observer.is_some()
                || run.ranges.input.placements.ranged()
                || run.ranges.subject_placements.ranged()))
}

/// The Seq-aligns' engine counts (outfmt 8, 12), their identities for the outfmt 0 rows and BLAST
/// XML's thrown HSPs when the reports read them (`report::seqalign_identities`), the display
/// counts of `hits` and the
/// outfmt 0 rows made up front (`rows_up_front`), with the positives of the search's matrix
/// `sc` (`write_tblastx_outputs`).
#[allow(clippy::too_many_arguments)]
fn report_counts<S: TbxScore>(
    sc: S,
    hits: &mut [TblastxHsp],
    outputs: &ReportOutputs<'_>,
    run: &TblastxReportRun<'_>,
    output_formats: &[report::TblastxOutputFormat],
    query_code: &GeneticCode,
    db_code: &GeneticCode,
) -> (
    Vec<(u32, u32)>,
    Option<Vec<u32>>,
    Option<Vec<crate::report::pairwise::PairwiseHit>>,
    // Under `-html`, the rows' masks of those rows (`HtmlRowMasks`, the spans of outfmt 0).
    Option<Vec<crate::report::query_anchored::HspMasks>>,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:518-522,1522-1530
    // (cited at `ReportSort::reads_hsp_scores`): only `-sorthits` and a `-sorthsps` above 0
    // read the HSPs' Seq-align scores in the pairwise report.
    // The rows of outfmt 0 and of the caller's hit records keep the engine's identities, read
    // before the display counts replace them, when the report reads them.
    let rows_read = run.args.formatting.report_sort().reads_hsp_scores()
        && (outputs.hits.is_some() || output_formats.iter().any(|format| format.is_pairwise()));
    // BLAST XML's `catch` of an HSP whose subject's lookup throws writes them too
    // (`XmlHsp::alignment_identities`): where the run writes BLAST XML and the load of some
    // subject of its database throws (`HitLoads::any_throws`). The one producer of both (review 1
    // of line DBID pfoff-4, L2).
    let xml_thrown = output_formats.contains(&report::TblastxOutputFormat::Xml)
        && run
            .database
            .as_ref()
            .is_some_and(|database| database.reported.loads.any_throws());
    let seqalign_identities = (rows_read || xml_thrown).then(|| report::seqalign_identities(hits));
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411
    // The caller receives the same final result that every formatter prints, with the
    // displayed rows of outfmt 0. Otherwise outfmt 0 makes its rows a window of queries at a
    // time (`rows_up_front`).
    let up_front = rows_up_front(outputs, run, output_formats);
    // Under `-html` the rows' masks of those rows (`HtmlRowMasks`).
    let mut pairwise_masks = (up_front && run.args.html_links.is_some()).then(Vec::new);
    let mut pairwise = up_front
        // EXPERIMENT (LOSAT_X_TBXRPTPAR): the rows as an indexed map on the pool.
        .then(|| match run.x_pool {
            Some(pool) => crate::algorithm::tblastx::x_rpt_par::pairwise_hits(
                pool,
                hits,
                seqalign_identities.as_deref(),
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
                pairwise_masks.as_mut(),
            ),
            None => report::pairwise_hits(
                hits,
                seqalign_identities.as_deref(),
                run.queries,
                run.subjects,
                run.subject_letters,
                run.query_stats,
                query_code,
                db_code,
                &run.ranges.input.placements,
                run.ranges.subject_placements,
                pairwise_masks.as_mut(),
            ),
        });
    // The rows count BLOSUM62's positives (`report::set_row_positives`).
    if !S::BLOSUM62 {
        if let Some(rows) = pairwise.as_mut() {
            report::set_row_positives(rows, sc.matrix());
        }
    }
    (engine_counts, seqalign_identities, pairwise, pairwise_masks)
}

/// The subject ID of the tabular formats of the subject record `subject`, `None` where
/// `GenerateDefline` throws (`tabular_subject_ids`).
fn subject_tabular_id(subject: &FastaRecord) -> Option<Arc<[u8]>> {
    subject
        .subject_tabular_id(false)
        .ok()
        .flatten()
        .map(Arc::from)
}

/// The subject ID of the tabular formats of every subject, by subject index
/// (`report::defline::tabular_subject_id`: the title's first word, or the local ID, and
/// for `lcl|Subject_...` the first word of `GenerateDefline`), as BLASTN's; `None` where
/// `GenerateDefline` throws (`common::tabular::SubjectId`), whose rows NCBI writes with its catch's values
/// (`report::tabular::SubjectLookup::Failed`) or stops at (`check_failed_lookups`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:863-870
fn tabular_subject_ids(
    subjects: &[FastaRecord],
) -> Vec<crate::algorithm::common::tabular::SubjectId> {
    subjects.iter().map(subject_tabular_id).collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:85-91
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:130-134

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:136-137
// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:176-207
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
) -> Result<()> {
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:253-259
    // What the run writes before its query batches ends here.
    if let Some(observer) = outputs.observer.as_deref_mut() {
        observer.queries_begin();
    }
    let query_records = query_runs.records;
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let input = query_runs.input;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:474-504
    // (`s_ReplaceLocalId`)
    // The query ID of the tabular formats is the title's first word, or the local ID
    // (`FastaRecord::shown_id`), as BLASTN's.
    let query_ids: Vec<Arc<[u8]>> = query_records
        .iter()
        .map(|record| Arc::from(record.shown_id()))
        .collect();
    let filtered_empty = match &subjects.shown {
        ShownSubjects::Database { database, .. } => database.filtered_empty(),
        ShownSubjects::Fasta { .. } => false,
    };
    let reports = BatchReports {
        filtered_empty,
        args: &args,
        query_runs,
        query_ids: &query_ids,
        subjects: &subjects.shown,
        numbers: &subjects.numbers,
        titles: std::cell::RefCell::new(DatabaseTitles::new()),
        stores: std::cell::RefCell::default(),
        lengths: std::cell::RefCell::default(),
        reads_letters: outputs.hits.is_some()
            || outputs.observer.is_some()
            || outputs
                .formats
                .iter()
                .any(|format| match report::output_format(format.outfmt) {
                    report::TblastxOutputFormat::CustomTabular(spec) => {
                        super::super::tabular::reads_letters(&spec.fields)
                    }
                    _ => true,
                }),
        x_search_pool: parallel_pool,
        queries_written: std::cell::Cell::new(0),
        reads: std::cell::RefCell::default(),
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:176-178
    // A reader at the end of its input reads no batch (its `End()` is true): the report of no
    // query, the prolog and the epilog (an empty pipe, or `-query -` after `-subject -`).
    let no_batch = query_runs.input_records.is_empty() && matches!(reading, QueryReading::End);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3103-3105, 3136 (the length
    // adjustment's tables of a name that the matrix loader read as IDENTITY,
    // raw_scoremat.c:157-161)
    // NCBI 2.17.0 frees arrays that it never set and aborts (`free(): invalid size`, rc 134)
    // when it sets up a batch; a run without a batch reads no matrix. LOSAT refuses a run with
    // a batch before its first byte.
    if !no_batch && args.scoring_matrix() == TblastxMatrix::Aborts {
        return Err(
            crate::algorithm::tblastx::args::identity_named_matrix_refusal(
                args.matrix_name(),
                "TBLASTX",
            ),
        );
    }
    // The same abort after NCBI reads the file of a name that its tables do not have
    // (`ApplicationSettings::matrix_file`, NCBI BLAST+ 2.17.0 with `BLASTMAT`: rc 134).
    if let (false, Some(file)) = (no_batch, args.matrix_file.as_deref()) {
        anyhow::bail!(
            "{}",
            crate::blastinput::ncbi_environment::matrix_file_refusal(
                args.matrix_name(),
                file,
                "TBLASTX"
            )
        );
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:173
    write_pairwise_prologs(outputs, &args, &subjects.shown, true)?;
    if no_batch {
        return reports.write(outputs, SearchedBatches::default(), &[], true);
    }
    // EXPERIMENT (LOSAT_X_TBXPIPE), dispatch: the batches the loop below will search, planned
    // before it (`x_pipe_plan`), when the pipe may run (`x_pipe_wanted`); `None` keeps the loop's
    // own searches.
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let x_plan = x_pipe_wanted(parallel_pool).then(|| {
        x_pipe_plan(
            input,
            &reading,
            batch_size,
            query_records,
            query_runs.input_records,
            batch_count(&args, subjects, outputs),
        )
    });
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:193-195
    // The search of one batch, the queries `range`; with LOSAT_X_TBXPIPE its seed-stage end
    // (`seed_end`) may start the next batch (tblastx/x_batch_pipe.rs).
    let search_subjects = subjects.subjects;
    let numbers = &subjects.numbers;
    // The search's matrix, chosen once (`scoring.rs`): BLOSUM62's instance of the search, or
    // the other matrices' with the table of `-matrix` (a name without a table stops at the first
    // batch with letters, `run_batch_loop`; one that NCBI aborts on runs no batch).
    // A word length other than 3 (`-word_size` 2 or 4) takes the other instance, BLOSUM62 too
    // (`TbxScore::word`).
    let word = args.word_size() as usize;
    let matrix_scores = match args.scoring_matrix() {
        TblastxMatrix::Unknown | TblastxMatrix::Aborts => None,
        TblastxMatrix::Table(ScoringMatrix::Blosum62) if word == 3 => None,
        TblastxMatrix::Table(matrix) => Some(MatrixScore::with_word(matrix, word)),
    };
    // EXPERIMENT (LOSAT_X_TBXBOARD): one board for the search's batches (tblastx/x_work_board.rs).
    // Only when the batches use it (`search_query_batch`, `x_board_batch`): with the serial chunk
    // scan, or several subjects without TBXNEXT's loop, the batches and the pipe run as before.
    let x_board = (crate::algorithm::tblastx::x_work_board::on()
        && parallel_pool.enabled()
        && parallel_pool.threads() > 1
        && std::env::var_os("LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS").is_none()
        && (search_subjects.count() == 1 || x_subject_next::on()))
    .then(crate::algorithm::tblastx::x_work_board::XBoard::default);
    let search: &x_batch_pipe::XSearch<'_, TblastxBatch> = &|range, seed_end| match &matrix_scores {
        None => search_query_batch(
            Blosum62Score,
            &args,
            &query_records[range],
            search_subjects,
            numbers,
            parallel_pool,
            seed_end,
            x_board.as_ref(),
            false,
        ),
        Some(scores) => search_query_batch(
            *scores,
            &args,
            &query_records[range],
            search_subjects,
            numbers,
            parallel_pool,
            seed_end,
            x_board.as_ref(),
            false,
        ),
    };
    // The set-up alone of the records of a refused batch that are first in their scope
    // (`run_batch_loop`; line WQA D10 round 5, 2026-11-20).
    let set_up = |records: &[FastaRecord]| match &matrix_scores {
        None => search_query_batch(
            Blosum62Score,
            &args,
            records,
            search_subjects,
            numbers,
            parallel_pool,
            None,
            None,
            true,
        ),
        Some(scores) => search_query_batch(
            *scores,
            &args,
            records,
            search_subjects,
            numbers,
            parallel_pool,
            None,
            None,
            true,
        ),
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:203-209
    // The loop writes each batch's reports when it takes the batch where no refusal of LOSAT's
    // can come after the first batch (`report_stream`); otherwise after the loop.
    let x_check = Instant::now();
    let stream = report_stream::check(
        &args,
        outputs,
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
            &set_up,
        )
    };
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    if let Some(plan) = x_plan {
        if let Some(pipe) =
            x_batch_pipe::XPipe::new(plan, search_subjects.count(), search, x_board.as_ref())
        {
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

/// NCBI's reader's count of the query records into the batches of a run (`BatchCount`), with
/// the scope of `run_batch_loop`'s `QueryScopes` (`None` without `-parse_deflines`). The pipe's
/// plan (`x_pipe_plan`) reads it.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn batch_count(
    args: &TblastxArgs,
    subjects: &RunSubjects<'_>,
    outputs: &ReportOutputs<'_>,
) -> Option<BatchCount> {
    let subject_records = matches!(subjects.shown, ShownSubjects::Fasta { .. });
    let xml = crate::blastinput::query_batch::outfmts_have_xml(
        outputs.formats.iter().map(|format| format.outfmt),
    );
    let span = crate::blastinput::query_batch::QueryScopeSpan::of(subject_records, xml);
    args.parse_deflines
        .then(|| BatchCount::new(span, "TBLASTX"))
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
    input_records: &[FastaRecord],
    mut count: Option<BatchCount>,
) -> Vec<std::ops::Range<usize>> {
    let mut plan = Vec::new();
    let reads_past_records = reading.reads_past_records();
    let reading_fails = matches!(reading, QueryReading::Error { .. });
    let mut reading_ended = false;
    let mut input_start = 0;
    while input_start < input.input_lengths.len() || (reads_past_records && !reading_ended) {
        let (input_end, reached_size) = BatchCount::next_batch(
            count.as_mut(),
            &|record| input_records[record].believed_shared(),
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
    set_up: &dyn Fn(&[FastaRecord]) -> Result<TblastxBatch>,
) -> Result<()> {
    let query_records = query_runs.records;
    let input = query_runs.input;
    let mut done = SearchedBatches::default();
    // The posts written before each query's report (`QueryWarnings`).
    let mut warnings: Vec<Posts> = vec![Posts::new(); query_records.len()];
    // Under `-parse_deflines`, the scopes of the batches' queries (`QueryScopes`), added as the
    // loop reads the batches. A streamed run has none that refuses (`report_stream::check`).
    // With outfmt 5 against a database, every query's ids from the start (`with_whole_run_of`).
    let subject_records = matches!(subjects.shown, ShownSubjects::Fasta { .. });
    let xml = args.parse_deflines
        && crate::blastinput::query_batch::outfmts_have_xml(
            outputs.formats.iter().map(|format| format.outfmt),
        );
    let mut query_scopes = crate::blastinput::query_batch::QueryScopes::of_search(
        args.parse_deflines,
        subject_records,
        xml,
        "TBLASTX",
    )
    .map(|scopes| {
        if xml && !subject_records {
            scopes.with_whole_run_of(query_runs.input_records)
        } else {
            scopes
        }
        .with_refusal(args.scope_refusal.clone())
        .with_xml_matrix_failure(reports.xml_matrix_failure(outputs))
    });
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1431-1452
    // A structured format (outfmt 5, 8, 9, 12) writes none of a query's own warnings
    // (`ReportFormat::writes_query_warnings`; outfmt 5 puts them in the Iteration's message,
    // `blast_xml`); the warnings go with the first format, or, where it is structured, they are
    // kept for the later formats that write them, each to the caller's receiver of its own
    // (`Posts::later`, `ReportOutputs::later_warnings`).
    let result_warnings = outputs.first_writes_query_warnings();
    let keeps_later = outputs.keeps_later_query_warnings();
    // Whether the batches' prefetch of a `-db` search writes the taxonomy database's messages
    // with the query warnings: the first format is outfmt 0, 5, 6 or 7, or with
    // `PRE_FETCH_SEQS_LIMIT=2147483647` 8, 9 or 12 (`TblastxOutputFormat::prefetches_subjects`;
    // see the batch's prefetch below).
    let prefetches = outputs
        .formats
        .first()
        .is_some_and(|format| report::output_format(format.outfmt).prefetches_subjects());
    // The first format, whose prefetch each batch asks (`database::prefetch_format`).
    let prefetch_format = super::super::database::prefetch_format(&outputs.formats);
    // The query records that NCBI's scope holds at the prefetch of the batch of the input records
    // `input` (`shared::PrefetchQueries`, made only where the batch stops): those that NCBI's
    // reader added for the batch (the records that `-query_loc` drops included, numbered in the
    // run, `FastaRecord::local_id`), and in outfmt 5, whose scope the formatter keeps
    // (`ResetScopeHistory`, cited at `shared::stop_lookup_unsure`), those of the batches before.
    let xml_scope = crate::blastinput::query_batch::outfmts_have_xml(
        outputs.formats.iter().map(|format| format.outfmt),
    );
    let batch_scope = |input: std::ops::Range<usize>| {
        let input = if xml_scope { 0..input.end } else { input };
        crate::blastinput::query_batch::QueryScope::of_shown_records(
            &query_runs.input_records[input],
        )
    };
    // The XML2 and JSON formats, whose prefetch NCBI runs too but which write the taxonomy texts
    // at each `Hit`: of the prefetch only its stop (`database::batch_prefetch_stop`).
    let xml2_prefetch = outputs.formats.first().is_some_and(|format| {
        matches!(
            report::output_format(format.outfmt),
            report::TblastxOutputFormat::Xml2(_)
        )
    });
    // Whether the first format writes the texts of the subjects that a query's report reads one
    // by one where it reads each (`ReportFormat::places_subject_reads`).
    let first_report = outputs
        .formats
        .first()
        .and_then(crate::api::local_blast::FormatOutput::report_format);
    let places_reads =
        first_report.is_some_and(crate::blastinput::app::ReportFormat::places_subject_reads);
    // The pairwise reports (outfmt 0 to 4) read the subjects that they show as their writers read
    // them (for the loads); the first format places the texts in the order of its reads.
    let pairwise = super::super::database::pairwise_reading(
        &outputs.formats,
        args.formatting.report_sort(),
        args.report_counts,
    );
    let pairwise_first = first_report
        .and_then(crate::report::pairwise::PairwiseView::of)
        .is_some();
    // Where a query's report of a batch whose prefetch is off ends the run at a subject whose
    // load throws, no report and no batch's prefetch after it loads anything
    // (`database::report_loads`).
    let run_loads = crate::report::db_subject::RunLoads::new(super::super::database::after_throw(
        outputs.formats.iter().map(|format| format.outfmt),
        outputs.other_outfmts,
        args.html_links.is_some(),
    ));
    // Whether the reports write the same where every load of a subject with a positive taxid throws
    // (`db_subject::unchanged_by_load_throws`).
    let reports_unchanged_by_load_throws =
        crate::report::db_subject::unchanged_by_load_throws_outfmts(
            outputs.formats.iter().map(|format| format.outfmt),
            crate::blastinput::app::Program::Tblastx,
            crate::report::archive::options::RequestKind::Tblastx,
        );
    let (query_code, db_code) = (
        crate::utils::genetic_code::GeneticCode::from_id(args.query_gencode),
        crate::utils::genetic_code::GeneticCode::from_id(args.db_gencode),
    );
    // A custom tabular report with a plural name looks every taxid of its rows up, taxid 0 too,
    // which the prefetch does not: the rows of a batch with hits may open the taxonomy database
    // before a later batch's prefetch (`TaxdbRun::open_for_rows`).
    let plural_names = outputs.formats.first().is_some_and(|format| {
        matches!(
            report::output_format(format.outfmt),
            report::TblastxOutputFormat::CustomTabular(spec) if spec.looks_up_every_taxid()
        )
    });
    // Which rows of the first format's report read the subjects' Bioseqs, which NCBI's
    // formatters fetch one by one where `PRE_FETCH_SEQS_LIMIT` turns a batch's prefetch off: a
    // custom tabular report's only for fields that read them (`TabularSpec::subject_fetch`).
    let rows = outputs
        .formats
        .first()
        .map_or(SubjectFetch::Every, |format| {
            match report::output_format(format.outfmt) {
                report::TblastxOutputFormat::CustomTabular(spec) => {
                    spec.subject_fetch(crate::blastinput::app::Program::Tblastx)
                }
                _ => SubjectFetch::Every,
            }
        });
    // Where the run's BLAST XML family formats fail for a matrix name that NCBI's XML does not
    // know (`BatchReports::xml_matrix_failure`; `None` for a name it knows).
    let xml_matrix_failure = reports.xml_matrix_failure(outputs);
    // The search strategy's requests of the batches (`-export_search_strategy`), with every
    // `-subject` record.
    let requests = args.strategy.as_ref().map(|_| {
        super::super::archive::requests(
            match &subjects.shown {
                ShownSubjects::Fasta { records, .. } => Some(*records),
                ShownSubjects::Database { .. } => None,
            },
            args.lcase_masking,
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
    // NCBI's reader's count of the records into the batches (`BatchCount`).
    let mut count = BatchCount::of(query_scopes.as_ref());
    // The reader is at the end of its input after the last record, unless an `eEOF` or an
    // error ends the reading there; that comes in the batch being read when the records
    // run out before the batch reaches its size, and otherwise in a batch of its own.
    while input_start < input.input_lengths.len()
        || (reading.reads_past_records() && !reading_ended)
    {
        let (input_end, reached_size) = BatchCount::next_batch(
            count.as_mut(),
            &|record| query_runs.input_records[record].believed_shared(),
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
        // NCBI's posts of the bytes replaced in the batch's search strategy request follow them
        // (below), so every error of the batch writes both.
        let mut batch_warnings = Posts::of_lines(
            query_runs.input_records[input_start..input_end]
                .iter()
                .flat_map(|record| record.warnings.iter().copied())
                .collect(),
        );
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:146-152
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
                batch_warnings.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
                outputs.diagnostics.write_all(&error_warnings)?;
                return Err(error.into_app_error());
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:378-380
        // A batch without a query (a batch size of 0, a batch of skipped records only with
        // `-query_loc`, or a batch that an `eEOF` ends before its first record) fails after
        // the reports of the batches before (without the epilog), when NCBI has read that
        // batch (its records' messages).
        if start == end {
            reports.write_before_failed_batch(outputs, done, &warnings)?;
            batch_warnings.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
            return Err(crate::cli::NativeError {
                exit: 3,
                message: format!("{ENGINE_ERROR} Empty CBlastQueryVector\n").into(),
            }
            .into());
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:178-181
        // LOSAT's refusal of two queries of the batch's scope with one id, in the one form of every
        // program (`query_batch::refused_batch_error`; line WQA D10 round 5, 2026-11-20): with
        // `-export_search_strategy` before the batch's request, which NCBI writes with the earlier
        // record's letters for the repeated id (no request of the batch is written); else after the
        // batch's checks before its set-up (the query splitter's environment, queries without
        // letters, the matrix name), where the set-up of the batch's records that are first in the
        // scope gives NCBI's error of the subjects or the refusal comes, before the search. The
        // stop is written as the Owner's decision Q-X2c of 2026-11-16 writes TBLASTX's refusals of
        // the search: after the reports of the batches before (for the first batch, the prologs)
        // and the messages of the batch's records. (Round 4 refused after the search, and wrote the
        // batch's request.)
        let mut scope_refusal = None;
        if let Some(scopes) = query_scopes.as_mut() {
            if let Err(refusal) = scopes.add_batch(
                &query_runs.input_records[input_start..input_end],
                &input.skipped[input_start..input_end],
                start,
            ) {
                scope_refusal = Some(refusal);
            }
        }
        if let Some(refusal) = scope_refusal.take_if(|_| args.strategy.is_some()) {
            reports.write_before_failed_batch(outputs, done, &warnings)?;
            batch_warnings.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
            return Err(refusal);
        }
        // The batch's search strategy request (`-export_search_strategy`): once its query factory
        // is made, before its set-up and search, so every error of the batch from here on follows
        // it (the file is apart from the reports). Its `EffectiveSearchSpace` option is the one
        // the batch before left (`note_searched` below). NCBI's posts of the bytes replaced in it
        // follow the reader's messages of the batch (`batch_warnings`).
        if let (Some(strategy), Some(requests)) = (args.strategy.as_ref(), requests.as_ref()) {
            requests.write(
                strategy,
                start,
                &query_records[start..end],
                input.placements.get(start),
            )?;
            batch_warnings.append(&strategy.take_all_posts());
        }
        // The query splitter reads its environment variables for every batch, before the
        // queries are set up (`check_query_split_environment`); the first batch meets an
        // error there.
        if input_start == 0 {
            if let Err(error) = crate::blastinput::app::check_query_split_environment(true, true) {
                reports.write_before_failed_batch(outputs, done, &warnings)?;
                batch_warnings.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
                return Err(error);
            }
        }
        let batch_queries = &query_records[start..end];
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1013-1025
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:923-929
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
            batch_warnings.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
            let mut message = format!("{ENGINE_ERROR} ");
            for _ in 0..no_data_queries {
                message.push_str("Warning: ");
                message.push_str(NO_DATA_MESSAGE);
                message.push(' ');
            }
            message.push('\n');
            return Err(crate::cli::NativeError {
                exit: 3,
                message: message.into(),
            }
            .into());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:158-182 (the score block of
        // the batch: a name that the matrix loader does not know fails with -1,
        // blast_stat.c:1679-1680, posted for all queries, blast_setup.c:498-501)
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:246 (the text of -1)
        // One message per query of the batch (`TSearchMessages::ToString` above, after
        // `RemoveDuplicates` in each query's list), `Error: ` first (`GetMessage`).
        // A query without letters keeps its own message (`Sequence contains no data`,
        // blast_setup_cxx.cpp:631-638, id is the query's index) after the error (context id
        // `kBlastMessageNoContext`, -1), which `RemoveDuplicates` sorts first
        // (`CSearchMessage::operator<`, blast_types.hpp:295-304).
        if args.scoring_matrix() == TblastxMatrix::Unknown {
            reports.write_before_failed_batch(outputs, done, &warnings)?;
            batch_warnings.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
            let text = crate::algorithm::common::unknown_matrix::unknown_matrix_message_of(
                batch_queries.iter().map(|record| {
                    if record.seq().is_empty() {
                        vec![crate::blastinput::query_batch::NO_DATA_MESSAGE.to_string()]
                    } else {
                        Vec::new()
                    }
                }),
            );
            return Err(crate::blastinput::app::engine_error(&text));
        }
        // The batch of LOSAT's refusal (`scope_refusal` above): NCBI's error of the set-up of its
        // records that are first in the scope (`set_up`, which ends before the search; the
        // database's kept read error where its sequences could not be read, as `take` fails such a
        // batch below), else the batch's search where that set-up passes and the run's XML failure
        // for its matrix name takes the refusal (`ScopeRefusal::xml_takes`; the failure comes
        // below), else the refusal.
        if let Some(refusal) = scope_refusal.take() {
            use crate::blastinput::query_batch::BatchSetUp;
            let refused = crate::blastinput::query_batch::refused_batch_error(
                refusal,
                false,
                |first| {
                    let records: Vec<FastaRecord> = (start..end)
                        .filter(|&query| first[input.ordinal(query) - input_start])
                        .map(|query| query_records[query].clone())
                        .collect();
                    match set_up(&records) {
                        Ok(batch) if batch.searched => BatchSetUp::Searched,
                        Ok(_) => BatchSetUp::Unknown,
                        Err(error) if error.is::<SubjectsUnreadable>() => stores_error
                            .take()
                            .and_then(|error| {
                                error
                                    .downcast::<crate::algorithm::common::database::PreliminaryReadError>()
                                    .ok()
                            })
                            .map(|read_error| crate::blastinput::app::engine_error(read_error.0))
                            .into(),
                        Err(error) => error
                            .downcast::<ErrorAfterBatchReports>()
                            .ok()
                            .map(|ErrorAfterBatchReports(error)| error)
                            .into(),
                    }
                },
            );
            if let Err(error) = refused {
                reports.write_before_failed_batch(outputs, done, &warnings)?;
                batch_warnings.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
                return Err(error);
            }
        }
        warnings[start].append(&batch_warnings);
        let mut batch = match take(start..end) {
            Ok(batch) => {
                // LOSAT_TIMING: the batch's block, kept in its result by a search of the pipe
                // (LOSAT_X_TBXPIPE), in the batch order.
                for line in &batch.x_timing {
                    eprintln!("{line}");
                }
                // The batch's set-up has written its `EffectiveSearchSpace` option, which the next
                // batch's search strategy request holds.
                if let Some(strategy) = args.strategy.as_ref() {
                    strategy.note_searched(
                        start..end,
                        super::super::archive::batch_searched_space(
                            batch.searched,
                            &batch.queries[0],
                        ),
                    )?;
                }
                batch
            }
            Err(error) if error.is::<SubjectsUnreadable>() || error.is::<XReadError>() => {
                // The first batch with a valid context of a `-db` search whose sequences could
                // not be read (`search`), or a read of the database that failed during the
                // batch's search (`XReadFailure`, with the error that the check gives for the same
                // read): NCBI's read error is an engine error after the reports of the batches
                // before; LOSAT's refusal (a truncated sequence file, unreadable ambiguity data)
                // stops the run, as any other refusal does.
                let error = match error.downcast::<XReadError>() {
                    Ok(XReadError(error)) => error,
                    Err(_) => stores_error
                        .take()
                        .expect("unreadable subjects keep their read error"),
                };
                let read_error =
                    error.downcast::<crate::algorithm::common::database::PreliminaryReadError>()?;
                reports.write_before_failed_batch(outputs, done, &warnings)?;
                batch_warnings.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
                return Err(crate::blastinput::app::engine_error(read_error.0));
            }
            Err(error) => match error.downcast::<ErrorAfterBatchReports>() {
                // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:203-209
                // The reports of the batches before are written (for the first batch, the
                // outfmt 0 prolog), then the messages of the failing batch's records (read
                // before its search), then the error, which skips `PrintEpilog`.
                Ok(ErrorAfterBatchReports(error)) => {
                    reports.write_before_failed_batch(outputs, done, &warnings)?;
                    batch_warnings.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
                    return Err(error);
                }
                Err(error) => return Err(error),
            },
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:283-293
        // `-db`: the definition lines of the batch's reported sequences are read where NCBI
        // makes the batch's Seq-aligns (`database::read_titles`, in report order); a header file
        // that NCBI cannot read fails this batch, after the reports of the batches before.
        // The taxonomy database's messages of each query's own fetches, where the batch's
        // prefetch is off (below), by query of the batch.
        let mut query_fetches: Vec<Vec<u8>> = Vec::new();
        if let ShownSubjects::Database { database, .. } = &subjects.shown {
            let aligned = super::super::database::drop_unaligned(database, &mut batch.hits);
            let earlier = crate::algorithm::common::database::EarlierUtf8::of_batch(
                start,
                batch.fetched_utf8,
                aligned,
            );
            let ordered = batch.hits.as_slice();
            let read = super::super::database::read_titles(
                database,
                ordered,
                start,
                (
                    earlier,
                    crate::algorithm::common::database::prefetch_may_be_off(prefetch_format),
                ),
                &mut reports.titles.borrow_mut(),
            );
            match read {
                // NCBI's warning of its first read of a UTF8String header follows the batch's
                // reader messages and request posts (`warnings[start]`), before its prefetch.
                Ok(read) => {
                    if let Some(posts) = crate::algorithm::common::database::Utf8Read::post_of(
                        read.as_ref(),
                        &database.db,
                        "tblastx",
                    ) {
                        warnings[start].append(&posts);
                        batch_warnings.append(&posts);
                    }
                }
                Err(error) => {
                    match error.downcast::<crate::algorithm::common::database::HeaderReadError>() {
                        Ok(failure) => {
                            reports.write_before_failed_batch(outputs, done, &warnings)?;
                            batch_warnings
                                .write(&mut *outputs.diagnostics, args.archive.as_ref())?;
                            if let Some(posts) =
                                crate::algorithm::common::database::Utf8Read::post_of(
                                    failure.utf8.as_ref(),
                                    &database.db,
                                    "tblastx",
                                )
                            {
                                posts.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
                            }
                            return Err(failure.error);
                        }
                        Err(error) => return Err(error),
                    }
                }
            }
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:197-202
        // The prefetch before the batch's reports reads `PRE_FETCH_SEQS_LIMIT` before it reads any
        // sequence, after the batch's search (a database's titles are read in its traceback,
        // above): a value that NCBI cannot convert stops the first batch here, after the
        // messages of its records and before its query warnings
        // (`app::check_prefetch_environment`).
        if let Err(error) =
            crate::blastinput::app::check_prefetch_environment(outputs.only_archive())
        {
            reports.write_before_failed_batch(outputs, done, &warnings)?;
            batch_warnings.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
            return Err(error);
        }
        if let ShownSubjects::Database { database, .. } = &subjects.shown {
            let ordered = batch.hits.as_slice();
            if let Some(scopes) = query_scopes.as_ref() {
                crate::algorithm::common::database::refuse_reported_ids_of_queries(
                    &reports.titles.borrow(),
                    ordered
                        .iter()
                        .map(|hsp| (start + hsp.hit.q_idx as usize, hsp.hit.s_idx)),
                    database.oid_map(),
                    scopes,
                    "TBLASTX",
                )?;
            }
            // The sequences that the batch reports, read once each where NCBI fetches them before
            // the batch's reports (`BlastFormatter_PreFetchSequenceData`, below), when a report
            // reads their letters, else their lengths alone (the index and the last byte of each):
            // with the reports of each batch (`report_stream`) they go after
            // them, else they stay to the reports after the loop, as the titles do. A read that
            // fails fails this batch as a read of its search does (above).
            if let ShownSubjects::Database {
                stream: Some(stream),
                ..
            } = &subjects.shown
            {
                let read = if reports.reads_letters {
                    super::super::database::read_reported(
                        stream,
                        ordered,
                        &mut reports.stores.borrow_mut(),
                    )
                } else {
                    super::super::database::read_reported_lengths(
                        stream,
                        ordered,
                        &mut reports.lengths.borrow_mut(),
                    )
                };
                if let Err(error) = read {
                    let XReadError(error) = error.downcast::<XReadError>()?;
                    let read_error = error
                        .downcast::<crate::algorithm::common::database::PreliminaryReadError>(
                    )?;
                    reports.write_before_failed_batch(outputs, done, &warnings)?;
                    batch_warnings.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
                    return Err(crate::blastinput::app::engine_error(read_error.0));
                }
            }
            // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:197-205
            // Before the reports of a batch NCBI fetches every subject that the batch reports
            // (the sequences of `ordered`), and with `taxdb.bti` on the search path the taxonomy
            // database writes its messages there (`database::prefetch_messages`: the texts of
            // opening it at the run's first lookup, a volume's lookup of a taxid once per run):
            // after the reader's messages of the batch and before the warnings of its first
            // query, which are added below. They go with the warnings, where the first format is
            // outfmt 0, 5, 6 or 7, or with `PRE_FETCH_SEQS_LIMIT=2147483647` 8, 9 or 12, which write
            // the reader's messages only (`prefetches`), so every writer of the reports writes them
            // there (the reports after the loop, the error paths that write the reports of the
            // batches before). Without `taxdb.bti` nothing is read. A refusal (a taxonomy
            // database or a taxid that LOSAT cannot read as NCBI does) stops the run before any
            // report, as the refusals of the titles above and XML2's lookups (`read_taxids`) do.
            // Where the prefetch fails at a subject with parsed Seq-ids, NCBI writes its warning
            // after those messages (`database::batch_prefetch`), with or without `taxdb.bti`.
            let taxdb_found = database.taxdb_found();
            // A lookup that fails on a file stops the run in the prefetch
            // (`shared::batch_prefetch`): after the reports of the batches before and this
            // batch's reader messages, as a header file that NCBI cannot read.
            macro_rules! prefetch_failed {
                ($error:expr) => {{
                    match $error.downcast::<crate::algorithm::common::database::HeaderReadError>() {
                        Ok(failure) => {
                            reports.write_before_failed_batch(outputs, done, &warnings)?;
                            batch_warnings
                                .write(&mut *outputs.diagnostics, args.archive.as_ref())?;
                            return Err(failure.error);
                        }
                        Err(error) => return Err(error),
                    }
                }};
            }
            let off = crate::algorithm::common::database::prefetch_off(
                prefetch_format,
                super::super::database::seq_aligns_of(ordered),
            );
            if xml2_prefetch && !run_loads.ended() && off {
                let (titles, mut reads) = (reports.titles.borrow(), reports.reads.borrow_mut());
                super::super::database::scope_utf8_hit_reads(
                    database,
                    ordered,
                    &titles,
                    start,
                    run_loads.answer().run(),
                    &mut reads,
                );
            }
            if xml2_prefetch
                && !run_loads.ended()
                && database.parsed_seqids()
                && !ordered.is_empty()
                && !off
            {
                let scope = || batch_scope(input_start..input_end);
                // The stop, or the warning of a UTF8String header that a load's lookup reads. The
                // batches before: their queries (`done.spans`, in the numbering of `done.hits`),
                // for what their reports load (`database::batch_prefetch_stop`).
                match super::super::database::batch_prefetch_stop(
                    database,
                    &done.hits,
                    done.searched.len(),
                    &done.spans,
                    (run_loads.answer(), prefetch_format),
                    &|hsp| {
                        super::super::database::hit_identities(
                            &reports.stores.borrow(),
                            hsp,
                            query_records[reports.queries_written.get() + hsp.hit.q_idx as usize]
                                .seq(),
                            &query_code,
                            &db_code,
                        )
                    },
                    ordered,
                    &reports.titles.borrow(),
                    &crate::algorithm::common::database::PrefetchQueries::new(start, &scope),
                ) {
                    Ok(warning) => warnings[start].lines(&warning),
                    Err(error) => prefetch_failed!(error),
                }
            }
            if prefetches && !run_loads.ended() && (taxdb_found || database.parsed_seqids()) {
                // The batch's Seq-aligns (`database::seq_aligns_of`).
                let seq_aligns = super::super::database::seq_aligns_of(ordered);
                // Where `PRE_FETCH_SEQS_LIMIT` turns the batch's prefetch off
                // (`shared::prefetch_off`), each query's report fetches the subjects it shows,
                // after the query's own warnings (`query_fetches`, added after them below), where
                // the first format's rows read their Bioseqs (`rows`; for `nident` alone the rows
                // of HSPs with the engine's identities, counted from the batch's stores, which a
                // report with `nident` reads above). NCBI looks the taxonomy up only where it
                // fetches a subject: in such a report, or in the prefetch of a batch with
                // Seq-aligns (one without them loads no id, `BlastFormatter_PreFetchSequenceData`).
                // The run's refusal of NCBI's taxonomy cache (`shared::TaxdbRun::cache_collision`)
                // is asked only there: the report stream is off for a run with that refusal
                // (`report_stream::Off::TaxonomyCache`), so every report is then written after
                // the loop and a later batch that fetches still refuses before any report (outfmt
                // 0's prolog is written before the search).
                if crate::algorithm::common::database::prefetch_off(prefetch_format, seq_aligns) {
                    // No prefetch; without `taxdb.bti` the reports' fetches look nothing up, but
                    // a lookup that reads a UTF8String header posts a warning (`load_read_texts`).
                    if taxdb_found || database.db.scope_utf8_seen() {
                        debug_assert!(ordered.is_sorted_by_key(|hsp| hsp.hit.q_idx));
                        let stores = reports.stores.borrow();
                        for query in 0..(end - start) as u32 {
                            let from = ordered.partition_point(|hsp| hsp.hit.q_idx < query);
                            let to = ordered.partition_point(|hsp| hsp.hit.q_idx <= query);
                            let query_letters = query_records[start + query as usize].seq();
                            let identities = |hsp: &TblastxHsp| {
                                super::super::database::hit_identities(
                                    &stores,
                                    hsp,
                                    query_letters,
                                    &query_code,
                                    &db_code,
                                )
                            };
                            // The first HSP of each subject whose row reads it, in report order (a
                            // subject's HSPs follow each other): the engine's identities are counted
                            // only for `nident` alone, until one of the subject's rows reads it.
                            // Not a subject that no report reaches after one whose load throws
                            // (`HitLoads::reached`; `database::report_loads`, after the reports
                            // before it).
                            let loads = super::super::database::report_loads(
                                (database, &run_loads),
                                &ordered[from..to],
                                &reports.titles.borrow(),
                                (pairwise, start),
                                identities,
                            );
                            let shown: Vec<&TblastxHsp> = match rows {
                                SubjectFetch::Every => ordered[from..to]
                                    .iter()
                                    .filter(|hsp| loads.reached(hsp.hit.s_idx as usize))
                                    .collect(),
                                SubjectFetch::Never => Vec::new(),
                                SubjectFetch::WithIdentities => {
                                    let mut read = std::collections::HashSet::new();
                                    ordered[from..to]
                                        .iter()
                                        .filter(|hsp| {
                                            loads.reached(hsp.hit.s_idx as usize)
                                                && !read.contains(&hsp.hit.s_idx)
                                                && identities(hsp) > 0
                                                && read.insert(hsp.hit.s_idx)
                                        })
                                        .collect()
                                }
                            };
                            // The first format writes the texts where it reads each subject
                            // (`ReportFormat::places_subject_reads`), else before the
                            // report.
                            query_fetches.push(if shown.is_empty() {
                                Vec::new()
                            } else if places_reads {
                                let titles = reports.titles.borrow();
                                let ordinals = super::super::database::read_ordinals(
                                    &shown,
                                    &ordered[from..to],
                                    &titles,
                                    pairwise.filter(|_| pairwise_first),
                                    &loads,
                                );
                                let texts = super::super::database::subject_read_texts(
                                    database,
                                    ordinals,
                                    &titles,
                                    plural_names,
                                )?;
                                drop(titles);
                                let mut reads = reports.reads.borrow_mut();
                                let index = start + query as usize;
                                if reads.len() <= index {
                                    reads.resize_with(index + 1, Vec::new);
                                }
                                reads[index] = texts;
                                Vec::new()
                            } else {
                                super::super::database::prefetch_messages(
                                    database,
                                    shown,
                                    &reports.titles.borrow(),
                                )?
                            });
                            // Only a row whose subject a report loads looks a name up: a row
                            // whose load throws keeps the names of the rows before (review 1 of
                            // line DBID pfoff-2, H1).
                            if plural_names
                                && ordered[from..to]
                                    .iter()
                                    .any(|hsp| loads.loaded(hsp.hit.s_idx as usize))
                            {
                                database.open_for_tabular_rows()?;
                            }
                        }
                    }
                } else if !ordered.is_empty() {
                    // Outfmt 5 first: whether NCBI's XML stream cut a report of the batches
                    // before (`shared::XmlCuts`), asked only where the batch's prefetch stops, by
                    // the record of their loads (`LoadedBefore`: a cut report loads nothing) and
                    // the texts that the cut reports left (`shared::XmlReportsBefore`).
                    let xml_cut = |query: usize| {
                        query_records.get(query).is_some_and(|record| {
                            crate::report::blast_xml::names::iteration_cut_before_hits(
                                record.local_id.as_bytes(),
                                record.warning_id(),
                                record.believed(),
                                &record.title,
                                false,
                            )
                        })
                    };
                    let reads = reports.reads.borrow();
                    // The posts before each query's report: those of the prefetch of the batch
                    // it begins (`warnings`, the batches before this one).
                    let prefetched = |query: usize| warnings.get(query);
                    let xml_cuts = crate::algorithm::common::database::XmlCuts::of(
                        first_report == Some(crate::blastinput::app::ReportFormat::Xml),
                        &xml_cut,
                    );
                    let xml_before = xml_cuts.reports_before(&reads, &prefetched);
                    let scope = || batch_scope(input_start..input_end);
                    let prefetch = super::super::database::batch_prefetch(
                        database,
                        super::super::database::LoadedBefore {
                            hits: &done.hits,
                            spans: &done.spans,
                            format: prefetch_format,
                            rows,
                            cuts: xml_cuts,
                        },
                        ordered,
                        &reports.titles.borrow(),
                        &crate::algorithm::common::database::PrefetchQueries::new(start, &scope)
                            .after_xml_reports(xml_before.as_ref())
                            .where_reports_unchanged_by_load_throws(
                                reports_unchanged_by_load_throws,
                            ),
                    );
                    drop(reads);
                    let messages = match prefetch {
                        Ok(messages) => messages,
                        Err(error) => prefetch_failed!(error),
                    };
                    if plural_names && !ordered.is_empty() {
                        database.open_for_tabular_rows()?;
                    }
                    warnings[start].append(&messages);
                }
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:213-218 (the
        // BLAST XML report data of a query that NCBI formats: the matrix of its formatting table,
        // found by the name as typed with `strcmp`)
        // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:197-205 (cited above: the
        // batch's prefetch, then `PrintOneResultSet` of each of its results)
        // NCBI fails where it formats the batch of the run's first searched query, or where the
        // report is held of NCBI's next query of another Seq-id (`BatchReports::xml_matrix_failure`;
        // `XmlMatrixFailure::SecondQueryOrEnd`, of the batches' scopes so far:
        // `QueryScopes::xml_failing_query`): after the batch's
        // search, the read of its titles (and LOSAT's refusal of two queries of one id, above) and
        // its prefetch, whose errors come first. It writes nothing of the batch's reports or of
        // the batches' before (the XML2 and JSON documents only their start, `PrintProlog`): the
        // reader messages and posts of the batches before (`write_before_failed_batch`), then
        // those of this batch up to the failing query, with its prefetch's messages. A run whose
        // held report has no such query fails after the loop. TBLASTX writes each batch's
        // search strategy request at the batch's start (above), so no request follows the
        // failing batch's.
        let failing = xml_matrix_failure
            .and_then(|failure| {
                crate::blastinput::query_batch::QueryScopes::xml_failing_query(
                    query_scopes.as_ref(),
                    failure,
                )
                .query()
            })
            .filter(|query| (start..end).contains(query));
        if let Some(failing) = failing {
            reports.write_before_failed_batch(outputs, done, &warnings)?;
            for posts in &warnings[start..=failing] {
                posts.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
            }
            return Err(crate::report::blast_xml::matrix::unsupported_score_matrix());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2815-2821
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:197-201,221
        // The warnings of a query come with its report (`PrintOneResultSet`), before its
        // preamble (`QueryWarnings`): a query without letters has the message of its set-up
        // (above); a translated query has no Karlin-Altschul message of its own, and the
        // queries of a batch that is not searched (no valid context) all get the message of
        // the batch.
        for (offset, query) in batch_queries
            .iter()
            .enumerate()
            .filter(|_| result_warnings || keeps_later)
        {
            let index = start + offset;
            let mut messages = Vec::new();
            if !batch.searched {
                messages.push(crate::report::query_warnings::INVALID_QUERY_MESSAGE.to_string());
            }
            if query.seq().is_empty() {
                messages.push(NO_DATA_MESSAGE.to_string());
            }
            if reports.filtered_empty && batch.searched {
                messages.push(crate::report::query_warnings::EMPTY_DATABASE_MESSAGE.to_string());
            }
            let text = crate::report::query_warnings::query_warning(
                "tblastx",
                reports.args.first_query_ordinal + input.ordinal(index),
                query,
                &messages,
            );
            match result_warnings {
                true => warnings[index].lines(&text),
                false => warnings[index].later(&text),
            }
        }
        for (offset, messages) in query_fetches.into_iter().enumerate() {
            warnings[start + offset].lines(&messages);
        }
        if stream {
            // The batch's reports, written now, as NCBI's loop prints a batch's result sets
            // (`report_stream`), from its lists numbered from its first query as it made them.
            // Its lists and its queries' warnings are not kept: a later batch's stop does not find
            // its hits among those the loop holds (the prefetches' `before`), so the run records
            // what NCBI's data loader loaded for the batch (`database::note_loaded`; nothing for
            // a database that cannot stop).
            if let ShownSubjects::Database { database, .. } = &subjects.shown {
                super::super::database::note_loaded(
                    database,
                    &batch.hits,
                    (start, 0..end - start),
                    &reports.titles.borrow(),
                    (prefetch_format, rows),
                );
            }
            debug_assert_eq!(reports.queries_written.get(), start);
            let count = end - start;
            let batch_done = SearchedBatches {
                hits: batch.hits,
                query_stats: batch.queries,
                searched: vec![batch.searched; count],
                spans: vec![0..count],
            };
            reports.write(outputs, batch_done, &warnings, false)?;
            reports.stores.borrow_mut().clear();
            reports.lengths.borrow_mut().clear();
            for query_warnings in &mut warnings[start..end] {
                *query_warnings = Posts::new();
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2204-2212 (cited at
    // `blast_xml::matrix::XmlMatrixFailure`: `PrintEpilog` makes the held report of the run's one
    // query)
    // A run whose held report has no query of another Seq-id fails at its end, instead of the
    // epilog, after what the batch loop wrote before a failing batch (the prologs, the batches'
    // reader messages and posts).
    if xml_matrix_failure.is_some_and(|failure| {
        crate::blastinput::query_batch::QueryScopes::xml_failing_query(
            query_scopes.as_ref(),
            failure,
        )
        .fails_at_end(done.searched.len())
    }) {
        reports.write_before_failed_batch(outputs, done, &warnings)?;
        return Err(crate::report::blast_xml::matrix::unsupported_score_matrix());
    }
    reports.write(outputs, done, &warnings, true)
}

/// NCBI's `CODON_LENGTH`: a query shorter than one codon has no frame of positive length.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:62-64
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
    filtered_empty: bool,
    args: &'a TblastxArgs,
    query_runs: &'a QueryRuns<'a>,
    query_ids: &'a [Arc<[u8]>],
    subjects: &'a ShownSubjects<'a>,
    numbers: &'a SubjectSetNumbers,
    /// `-db`: the stored titles of the batches' reported sequences (`run_in_pool`).
    titles: std::cell::RefCell<DatabaseTitles>,
    /// `-db`: the stores of the reported sequences of the batches not yet reported
    /// (`database::read_reported`), read when a report reads letters (`reads_letters`).
    stores: std::cell::RefCell<super::super::database::ReportedStores>,
    /// `-db`: the lengths of the reported sequences of the batches not yet reported
    /// (`database::read_reported_lengths`), read when no report reads letters.
    lengths: std::cell::RefCell<super::super::database::ReportedLengths>,
    /// Whether a report reads the letters of the subjects: every format but a custom tabular
    /// report without a field that reads them (`tabular::reads_letters`), and the hits and the
    /// observer of a library caller (`report_counts`).
    reads_letters: bool,
    /// EXPERIMENT (LOSAT_X_TBXRPTPAR): the search's pool, for the report stage
    /// (`TblastxReportRun::x_pool`).
    x_search_pool: &'a crate::utils::threading::SearchPool<'a>,
    /// The queries whose reports were written: those of the batches that the loop wrote when it
    /// took them (`report_stream`), which the next call follows.
    queries_written: std::cell::Cell<usize>,
    /// The texts of the subjects that each query's report reads one by one (its batch's prefetch
    /// off), by query, where the first format places them (`TblastxReportRun::query_reads`).
    reads: std::cell::RefCell<Vec<crate::report::query_warnings::QueryReads>>,
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
        warnings: &[Posts],
        epilog: bool,
    ) -> Result<()> {
        // Diagnostic (LOSAT_TIMING): the report for the common line (utils/xtiming.rs); written
        // after the batches are searched. No NCBI counterpart.
        let t_output = crate::utils::xtiming::clock();
        let first = self.queries_written.get();
        let reported = first..first + done.searched.len();
        // The texts of the reads of the reported queries (`BatchReports::reads`; a query after
        // the last one that reads a subject has none).
        let mut query_reads: Vec<crate::report::query_warnings::QueryReads> = {
            let reads = self.reads.borrow();
            reported
                .clone()
                .map(|query| reads.get(query).cloned().unwrap_or_default())
                .collect()
        };
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
                database, pairwise, ..
            } => {
                let formats: Vec<report::TblastxOutputFormat> = outputs
                    .formats
                    .iter()
                    .map(|format| report::output_format(format.outfmt))
                    .collect();
                // The first format, whose prefetch the batches ask (`database::prefetch_format`).
                let prefetch_format = super::super::database::prefetch_format(&outputs.formats);
                let pairwise_shown = formats.iter().any(|format| format.is_pairwise());
                let shown = crate::algorithm::common::database::ParsedShown {
                    tabular: formats.iter().any(|format| format.tabular_reads_subject()),
                    stops_at_failed_lookup: formats.iter().any(|format| {
                        format.is_tabular()
                            && crate::report::tabular::may_stop_at_failed_lookup(
                                format.tabular_spec(),
                            )
                    }),
                    pairwise: pairwise_shown,
                    structured: formats.iter().any(|format| format.is_structured_xml()),
                    joined: pairwise_shown
                        || formats.iter().any(|format| {
                            format.is_structured_xml() || format.tabular_joins_subjects()
                        }),
                    sam: false,
                    taxonomy: formats.iter().any(|format| format.writes_taxonomy()),
                    words: formats.iter().fold(
                        crate::report::tabular::PrintedIds::default(),
                        |words, format| words.or(format.printed_ids()),
                    ),
                    prefetched: formats
                        .first()
                        .is_some_and(|format| format.prefetches_subjects())
                        && !crate::algorithm::common::database::prefetch_always_off(
                            prefetch_format,
                        ),
                    one_by_one: crate::algorithm::common::database::prefetch_may_be_off(
                        prefetch_format,
                    ),
                };
                // The queries of each batch (`done.spans`, in the numbering of `hits`), whose
                // prefetch decides the subjects that only the prefetch reads; the hits' engine
                // identities tell which rows of `nident` alone read a subject.
                let (query_code, db_code) = (
                    crate::utils::genetic_code::GeneticCode::from_id(self.args.query_gencode),
                    crate::utils::genetic_code::GeneticCode::from_id(self.args.db_gencode),
                );
                reported_subjects = super::super::database::reported_subjects(
                    database,
                    &self.titles.borrow(),
                    &mut hits,
                    self.reads_letters.then_some(&*self.stores.borrow()),
                    &self.lengths.borrow(),
                    (shown, prefetch_format),
                    (
                        super::super::database::after_throw(
                            outputs.formats.iter().map(|format| format.outfmt),
                            outputs.other_outfmts,
                            self.args.html_links.is_some(),
                        ),
                        // The pairwise reports' reads (as `run_batch_loop`'s).
                        super::super::database::pairwise_reading(
                            &outputs.formats,
                            self.args.formatting.report_sort(),
                            self.args.report_counts,
                        ),
                    ),
                    &|hsp| {
                        super::super::database::hit_identities(
                            &self.stores.borrow(),
                            hsp,
                            queries[hsp.hit.q_idx as usize].seq(),
                            &query_code,
                            &db_code,
                        )
                    },
                    done.searched.len(),
                    &done.spans,
                )?;
                // The texts of the reads, made before the hits named the subjects by rank, by
                // rank as the writers name them.
                query_reads = crate::report::query_warnings::reads_by_rank(
                    &query_reads,
                    &reported_subjects.ordinals,
                );
                // The library's records name a database's sequences through this event.
                if outputs.hits.is_some() {
                    if let Some(observer) = outputs.observer.as_deref_mut() {
                        observer.database_subjects(reported_subjects.library_subjects(database));
                    }
                }
                (
                    reported_subjects.records.as_slice(),
                    report::SubjectLetters::Database {
                        lengths: &reported_subjects.lengths,
                        stores: &reported_subjects.stores,
                    },
                    Some(ReportDatabase {
                        name: database.name(),
                        pairwise,
                        reported: &reported_subjects,
                        oids: database.oid_map(),
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
        // A `-subject` search's organism report (outfmt 18) looks names up in the run's
        // database.
        let tabular_names = match (&self.subjects, custom_tabular) {
            (ShownSubjects::Database { database, .. }, Some(spec)) => {
                database.tabular_names(spec, &tabular_taxonomy)?
            }
            (_, Some(_)) => self
                .args
                .subject_taxdb
                .as_ref()
                .map(|taxdb| taxdb.report_names("TBLASTX")),
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
                filtered_empty: self.filtered_empty,
                warnings: &warnings[reported.clone()],
                query_reads: &query_reads,
                ranges,
                epilog,
                queries_before: first,
                numbers: self.numbers,
                xml_taxa: &xml_taxa,
                tabular_taxonomy: &tabular_taxonomy,
                tabular_names: tabular_names.as_ref(),
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
    /// when it took them (`report_stream`), nothing remains. (The failing batch's search strategy
    /// request, when NCBI made its query factory, was written before its search,
    /// `run_batch_loop`.)
    fn write_before_failed_batch(
        &self,
        outputs: &mut ReportOutputs<'_>,
        done: SearchedBatches,
        warnings: &[Posts],
    ) -> Result<()> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2362-2389 (quoted
        // at `blast_xml::matrix::XmlMatrixFailure`)
        // With a matrix name that NCBI's XML does not know, the batches before a failing batch
        // wrote no report: the XML formats are not streamed (`report_stream`), and their queries
        // come before the failing query (`xml_matrix_failure`), the one query of a run whose
        // report is held, which NCBI never writes. What they wrote is their reader messages and
        // posts, after the prologs.
        let xml_matrix_failed = self.xml_matrix_failure(outputs).is_some();
        if self.queries_written.get() > 0 {
            Ok(())
        } else if done.searched.is_empty() || xml_matrix_failed {
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
            for posts in &warnings[..done.searched.len()] {
                posts.write(&mut *outputs.diagnostics, self.args.archive.as_ref())?;
            }
            Ok(())
        } else {
            self.write(outputs, done, warnings, false)
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:121 (cited at
    // `blast_xml2::reports_held`)
    /// Where the run of the formats of `outputs` stops for a matrix name that NCBI's BLAST XML
    /// report data does not know (`blast_xml::matrix::xml_formatting_failure`), with the facts of
    /// the search: subject records (`-subject`, `m_IsBl2Seq`); TBLASTX does not iterate. The one
    /// decision of the batch loop's error (`run_batch_loop`).
    fn xml_matrix_failure(
        &self,
        outputs: &ReportOutputs<'_>,
    ) -> Option<crate::report::blast_xml::matrix::XmlMatrixFailure> {
        crate::report::blast_xml::matrix::xml_formatting_failure(
            outputs.formats.iter().map(|format| format.outfmt),
            self.args.matrix_name(),
            matches!(self.subjects, ShownSubjects::Fasta { .. }),
            false,
        )
    }
}

/// TBLASTX's errors of a query batch's search that come after the reports of the batches before
/// (`blastinput::app::ErrorAfterBatchReports`, which TBLASTN shares): NCBI's (`run_in_pool`), the
/// average subject length (`context_cutoffs::context_cutoffs`); and LOSAT's refusals of what the
/// search meets, an ungapped extension that scores higher past the end of its query frame
/// (`extension::query_batch::refusal`, the Owner's decision of 2026-11-16 on line RG's Q-X2c: "the
/// results of the query batches before are written"), a one-hit merge of lengths that add up to 0
/// (`OneHitChunkMerge::refuse`), the ungapped re-evaluation of HSPs on a subject (or its
/// `-subject_loc` range) of fewer than 12 bases (`subject_loop::SubjectEmit::refuse_untranslated_frame`,
/// PD-LOSAT-NCBI-DEFECTS 1.51) and the re-evaluation's fetch of a subject whose header cannot be read
/// under a negative seqid list (`subject_loop::x_fetch_keeps`, NCBI's engine error there; the inner
/// error stays a `FetchRefusal`). The run then writes the same bytes whether its reports are
/// streamed or held (`report_stream`).
pub(crate) use crate::blastinput::app::ErrorAfterBatchReports;

// EXPERIMENT (LOSAT_X_TBXBOARD): closes subject `.1`'s units on the board when its search returns
// (a no-op once it put its units there). No NCBI counterpart.
struct XCloseUnits<'a>(&'a crate::algorithm::tblastx::x_work_board::XBatch, usize);

impl Drop for XCloseUnits<'_> {
    fn drop(&mut self) {
        self.0.close_without_units(self.1);
    }
}

// EXPERIMENT (LOSAT_X_TBXBOARD): the batch leaves the board on every return of its search. No NCBI
// counterpart.
struct XEndBatch<'a>(
    &'a crate::algorithm::tblastx::x_work_board::XBoard,
    &'a std::sync::Arc<crate::algorithm::tblastx::x_work_board::XBatch>,
);

impl Drop for XEndBatch<'_> {
    fn drop(&mut self) {
        self.0.end_batch(self.1);
    }
}

/// Searches one query batch against the subjects and returns its hits and query statistics: the
/// set-up of the batch (the flags, the queries' frames, the lookup table: `batch_setup.rs`; the
/// diagonal table, the offset arrays and the seed release: `seed_tables.rs`; the cutoffs, length
/// adjustments and search spaces of the contexts: `context_cutoffs.rs`), then `process_subject`
/// for every subject (`subject_loop.rs`, with the end of each subject, `SubjectEmit`; the
/// statistics each subject reads: `context_cutoffs.rs`; its linking parameters, the linking, the
/// reap and the re-evaluation of its HSP list: `subject_link.rs`; the stages after the
/// re-evaluation and its report hits: `subject_report.rs`, with the diagnostics of
/// `subject_debug.rs`) and the collection of their hits (`batch_finish.rs`).
/// These preliminary hit lists are the final ones: without gapped calculation NCBI's traceback
/// stage aligns nothing, it only computes bit scores
/// (c++/src/algo/blast/core/blast_traceback.c:1526, 1691-1696; `hit_collector.rs`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1334-1607
/// (`BLAST_PreliminarySearchEngine`), with the set-up of
/// c++/src/algo/blast/api/blast_setup_cxx.cpp:484-658 (`SetupQueries_OMF`),
/// c++/src/algo/blast/core/lookup_wrap.c:61-195 (`LookupTableWrapInit_MT`),
/// c++/src/algo/blast/core/blast_parameters.c:279-419 (`BlastInitialWordParametersUpdate`) and
/// c++/src/algo/blast/core/blast_setup.c:699-850 (`BLAST_CalcEffLengths`).
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
    // EXPERIMENT (LOSAT_X_TBXBOARD): the search's board of units and link groups, `None` without it
    // (tblastx/x_work_board.rs).
    x_board: Option<&crate::algorithm::tblastx::x_work_board::XBoard>,
    // The set-up of the batch alone, up to its checks of the subjects (`context_cutoffs`): a batch
    // that passes them with a valid context returns without hits, `searched` (NCBI searches it).
    // The set-up of a batch that LOSAT refuses (`query_batch::refused_batch_error`; line WQA D10
    // round 5, 2026-11-20).
    setup_only: bool,
) -> Result<TblastxBatch> {
    // The search's flags, threads, genetic codes, window and word size, and the ungapped block of
    // the X-dropoffs (`batch_setup::batch_flags`).
    let batch_setup::BatchFlags {
        extension_debug_enabled,
        timing_enabled,
        t_total,
        num_threads,
        requested_parallel,
        use_serial_scan_chunks,
        query_code,
        db_code,
        window,
        wordsize,
        ungapped_params_for_xdrop,
        diag_enabled,
        debug_cutoffs_all,
        debug_output_filter,
        debug_hsp_saving,
        diagnostics,
        scan_debug_range,
    } = batch_setup::batch_flags(sc, args, parallel_pool);

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:485-651
    let queries_raw: &[FastaRecord] = query_records;
    // The six frames of the batch's queries with their masks (`batch_setup::query_frames`).
    let batch_setup::QueryFrames {
        query_frames,
        t_read_queries,
    } = batch_setup::query_frames(args, queries_raw, &query_code, parallel_pool)?;

    // The batch's lookup table, its contexts and their X-dropoffs (`batch_setup::batch_lookup`).
    let batch_setup::BatchLookup {
        x_soft_subjects,
        x_pair_lut_batch,
        x_pair_lut_is_stub,
        lookup,
        contexts,
        t_build_lookup,
        x_dropoff_per_context,
    } = batch_setup::batch_lookup(&batch_setup::LookupInputs {
        args,
        subjects,
        query_frames: &query_frames,
        window,
        scan_debug_range,
        extension_debug_enabled,
        diag_enabled,
        debug_hsp_saving,
        requested_parallel,
        use_serial_scan_chunks,
        parallel_pool,
        sc,
        ungapped_params_for_xdrop,
    });

    let t_phase_read_subjects = Instant::now();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1427
    // An unreadable database is not empty: its first batch with a valid context fails
    // (`SubjectsUnreadable`, below). A database that a list left without sequences has
    // statistics (`SeqDb::statistics_numbers`) and its queries are reported with them.
    if queries_raw.is_empty()
        || (subjects.count() == 0
            && numbers.total == 0
            && !matches!(subjects, Subjects::Unreadable))
    {
        return Ok(TblastxBatch {
            hits: Vec::new(),
            searched: true,
            queries: vec![TblastxQueryStats::default(); queries_raw.len()],
            x_timing: Vec::new(),
            fetched_utf8: None,
        });
    }
    let t_read_subjects = t_phase_read_subjects.elapsed();

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1475
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:452-584
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-188
    let use_parallel = requested_parallel;
    let use_parallel_chunks = requested_parallel
        && !use_serial_scan_chunks
        && (std::env::var_os("LOSAT_TBLASTX_PARALLEL_CHUNKS").is_some()
            || (cfg!(all(target_arch = "wasm32", feature = "wasm-threads"))
                && tblastx_estimated_subject_chunk_work_items(
                    subjects.longest_and_total().0,
                    tblastx_max_dbseq_len_for_run(),
                ) > 1));
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
    // For bit score and E-value calculation, NCBI uses the ungapped `kbp` table on `sbp`:
    //   blast_hits.c line 1833 selects `kbp_gap` when gapped or `kbp` when ungapped;
    //   blast_hits.c line 1918 follows the same selection in `Blast_HSPListGetBitScores`.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:768
    // For tblastx, `gapped_calculation` is false, so `kbp_ptr` points to the ungapped `kbp` table on `sbp`.
    // Therefore, ALL calculations (eff_searchsp, cutoff, bit score, E-value) use UNGAPPED params
    //
    // (the contexts' blocks, `QueryContext::karlin_params`).

    // Compute NCBI-style average query length for linking cutoffs
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:998-1082
    // NCBI uses average over ALL contexts including zero-length (frame restriction via strand)
    let query_nucl_lengths: Vec<usize> = queries_raw.iter().map(|r| r.seq().len()).collect();
    let avg_query_length = compute_avg_query_length_ncbi(
        &query_nucl_lengths,
        crate::algorithm::blastn::args::QueryStrand::from_option(&args.strand).searched(),
    );

    let t_search_start = Instant::now();
    let scan_ns = AtomicU64::new(0);
    let scan_calls = AtomicU64::new(0);
    let ungapped_ns = AtomicU64::new(0);
    let ungapped_calls = AtomicU64::new(0);
    let reeval_ns = AtomicU64::new(0);
    let reeval_calls = AtomicU64::new(0);
    let linking_ns = AtomicU64::new(0);
    let linking_calls = AtomicU64::new(0);
    // No NCBI counterpart: the stage timer of LOSAT_TIMING (x_search_timing.rs). It reads the
    // clock only when `timing_enabled`, and the search reads nothing from it.
    let x_times = x_search_timing::XSearchTimes::new(timing_enabled);
    // No NCBI counterpart: the linking timer of LOSAT_TIMING (x_link_timing.rs). It arms a
    // linking call only when `timing_enabled`, and the search reads nothing from it.
    let x_link_timer = XLinkTimer::new(timing_enabled);

    // NCBI verbose CLI flag is only present under _BLAST_DEBUG builds.
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3265-3274
    let bar = ProgressBar::hidden();

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1475
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
    // reach them as one message at the end of the subject (`subject_loop::SubjectEmit::emit`):
    // through the channel to the writer thread, or, without one (a serial search, the WASI
    // targets), into the collector the subject loop keeps here.
    let local_collector: std::sync::Mutex<Option<HitCollector>> =
        std::sync::Mutex::new(if use_channel {
            None
        } else {
            Some(HitCollector::new(
                args,
                query_records,
                numbers.link_total > 0,
            ))
        });

    // The sizes of the diagonal table and the offset arrays, the batch's unit path (PAIRPAR,
    // dispatch point 1) and its batch on the board (`seed_tables::scan_tables`).
    let seed_tables::ScanTables {
        diag_array_size,
        diag_mask,
        offset_array_size,
        offset_pairs_len,
        x_scan_bf,
        x_pair_par_batch,
        x_board_batch,
    } = seed_tables::scan_tables(
        &seed_tables::ScanInputs {
            contexts: &contexts,
            lookup: &lookup,
            window,
            subjects,
            scan_debug_range,
            extension_debug_enabled,
            diag_enabled,
            debug_hsp_saving,
            x_soft_subjects,
            use_serial_scan_chunks,
            use_parallel_chunks,
            use_parallel,
            parallel_pool,
            args,
            sc,
        },
        x_pair_lut_batch,
        x_board,
    )?;
    // EXPERIMENT (LOSAT_X_TBXBOARD): the batch leaves the board on every return of this search
    // (`XEndBatch`).
    let _x_board_end = x_board_batch
        .as_ref()
        .map(|(board, batch)| XEndBatch(board, batch));
    // The holders of the seed release, which take the lookup table and the unit path, the size of
    // the workers' diagonal tables and the extensions' sentinel walks (`seed_tables::seed_batch`).
    let seed_tables::SeedBatch {
        #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
        x_board_unit_path,
        x_seed_batch,
        x_worker_diag_cells,
        tbx_sent,
    } = seed_tables::seed_batch(
        &seed_tables::SeedInputs {
            x_board_batch: &x_board_batch,
            subjects,
            diag_array_size,
            contexts: &contexts,
            sc,
            extension_debug_enabled,
        },
        lookup,
        x_pair_par_batch,
    );

    let contexts_ref = &contexts;

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-180,261-264;
    // c++/src/algo/blast/core/blast_engine.c:1372-1374,1434-1443
    // When total database length is zero, NCBI updates parameters per subject via `BLAST_OneSubjectUpdateParameters`.
    // The set's numbers are read once per run (`SubjectSetNumbers`: the records' sums, or the
    // database index's numbers).
    let db_length_nucl = numbers.total;
    let db_num_seqs = numbers.count;
    // The linking cutoffs read `BlastSeqSrcGetTotLen`, not the statistics' total
    // (`SubjectSetNumbers::link_total`).
    let link_db_length_nucl = numbers.link_total;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:535-560
    // NCBI returns success only when at least one valid query context exists.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-178,206-208
    // A batch without a valid context is not searched; each of its queries is reported
    // with the `-1` statistics.
    if !contexts_ref.iter().any(|ctx| ctx.is_valid) {
        return Ok(TblastxBatch {
            hits: Vec::new(),
            searched: false,
            queries: vec![TblastxQueryStats::default(); queries_raw.len()],
            x_timing: Vec::new(),
            fetched_utf8: None,
        });
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980
    // A count of 0 (of a positive total) is the division by 0 where NCBI's setup of this first
    // batch with a valid context crashes, before it reads a sequence: LOSAT's refusal
    // (`common::database::check_setup_count`, D1a-5 of line RG), after the reports of the batches
    // before, also where the sequences could not be read (NCBI 2.17.0: `NSEQ 0` over a volume
    // without its `.nsq`, SIGFPE; review A of D1a-5, F1). A run whose batches have no valid context
    // does not reach it and is reported as NCBI's.
    crate::algorithm::common::database::check_setup_count(
        db_length_nucl,
        i32::try_from(db_num_seqs).unwrap_or(i32::MAX),
        false,
        "TBLASTX",
    )
    .map_err(ErrorAfterBatchReports)?;
    // `-db` whose sequences could not be read: this first batch with a valid context fails
    // (`run_in_pool` raises the kept read error).
    if matches!(subjects, Subjects::Unreadable) {
        return Err(SubjectsUnreadable.into());
    }
    // `-db`: the read of the database that failed in this batch's search, if one does (a database
    // changed or damaged after the check, an error of the file system): the batch fails with it
    // once every subject before it is searched (`XReadFailure`).
    let x_read_failure = crate::algorithm::tblastx::subjects::XReadFailure::new();
    // The average subject length check, then the cutoffs, length adjustments and search spaces of
    // every context, the queries' report statistics and the smallest-Lambda block of the linking
    // cutoffs (`context_cutoffs::context_cutoffs`).
    let context_cutoffs::ContextCutoffs {
        cutoff_scores,
        cutoff_score_min,
        length_adj_per_context,
        eff_searchsp_per_context,
        query_stats,
        linking_params_for_cutoff,
        option_spaces,
    } = context_cutoffs::context_cutoffs(&context_cutoffs::CutoffInputs {
        contexts_ref,
        numbers,
        db_length_nucl,
        db_num_seqs,
        subjects,
        evalue_threshold,
        diag_enabled,
        debug_cutoffs_all,
        query_frames: &query_frames,
    })?;
    // The set-up passed with a valid context: NCBI searches the batch (`searched`).
    if setup_only {
        return Ok(TblastxBatch {
            hits: Vec::new(),
            searched: true,
            queries: vec![TblastxQueryStats::default(); queries_raw.len()],
            x_timing: Vec::new(),
            fetched_utf8: None,
        });
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:82-88
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let writer = if use_channel {
        let rx = rx_opt.take().expect("rx must be available for writer");
        // The writer thread (`subject_loop::spawn_writer`).
        Some(subject_loop::spawn_writer(
            rx,
            args,
            query_records,
            numbers.link_total > 0,
        ))
    } else {
        None
    };

    // The negative seqid list's check of the re-evaluation's fetch, one for the batch
    // (`FetchExcluded`), which each subject with HSPs left for the re-evaluation reads
    // (`subject_loop::SubjectEmit::reevaluation_fetch_keeps`).
    let fetch_excluded = subjects
        .seqdb()
        .and_then(crate::algorithm::common::database::FetchExcluded::of);
    // One message per subject of the loop, with its final HSPs or none
    // (`subject_loop::SubjectEmit::emit`).
    let subject_emit = subject_loop::SubjectEmit {
        fetch_excluded: fetch_excluded.as_ref(),
        x_read_failure: &x_read_failure,
        local_collector: &local_collector,
    };

    // What the word finder of the frame loop reads of this batch (`word_finder::WordFinderCx`),
    // in the batch's statistics below; `process_subject` makes a copy with the subject's
    // `cutoff_scores` where NCBI computes the statistics of each subject.
    let word_finder_cx = word_finder::WordFinderCx {
        window,
        wordsize,
        diag_mask,
        diag_array_size,
        diag_enabled,
        diagnostics: &diagnostics,
        timing_enabled,
        extension_debug_enabled,
        contexts_ref,
        cutoff_scores: &cutoff_scores,
        x_dropoff_per_context: &x_dropoff_per_context,
        scan_ns: &scan_ns,
        scan_calls: &scan_calls,
        ungapped_ns: &ungapped_ns,
        ungapped_calls: &ungapped_calls,
        x_scan_bf,
        offset_array_size,
        tbx_sent,
        scan_debug_range,
        x_read_failure: &x_read_failure,
        query_batch: QueryBatch::new(contexts_ref, &x_read_failure),
        subject_besthit: args.subject_besthit,
    };
    // The batch's statistics (`SubjectStatistics`); `per_subject`: NCBI computes the statistics
    // of each subject (a database length of 0) with the preliminary stage's option spaces.
    let batch_statistics = context_cutoffs::SubjectStatistics {
        cutoff_scores: &cutoff_scores,
        cutoff_score_min,
        length_adj_per_context: &length_adj_per_context,
        eff_searchsp_per_context: &eff_searchsp_per_context,
        word_finder_cx,
    };
    let subject_option_spaces = option_spaces.as_deref().filter(|_| numbers.per_subject());

    let process_subject = |st: &mut WorkerState<S>, s_idx: usize| {
        // No NCBI counterpart: the clock of this subject for LOSAT_TIMING. The first local of
        // the closure, so it is dropped last and its last stage includes the drops of the rest.
        let mut x_lap = x_times.subject();
        // EXPERIMENT (LOSAT_X_TBXBOARD): this thread searches a subject of the board's batch; a
        // subject that returns before its units puts none on the board.
        let _x_board_subject = x_board_batch.as_ref().map(|(board, batch)| {
            (
                crate::algorithm::tblastx::x_work_board::enter(board, batch),
                XCloseUnits(batch, s_idx),
            )
        });
        // The scores of the search's matrix, an argument of the frame loop below.
        let sc = st.sc;
        // EXPERIMENT (LOSAT_X_TBXPIPE): this subject's part of the batch's seed-stage end, ended at
        // dispatch point 5 (`subject_link::end_subject_seed_stage`) with the subject's HSPs, or
        // with none when the closure returns before it (tblastx/x_batch_pipe.rs). Nothing
        // without the pipe.
        let mut x_pipe_seed = x_batch_pipe::XSeedGuard::new(x_seed_end);
        // Dispatch point 2 of the seed release (per subject): the seed stage of
        // this subject holds the batch's lookup table through `x_seed_stage` and the unit path
        // through `x_units_stage`. `lookup_ref` borrows from the first and `x_pair_par_ref`
        // (`subject_prep::run_unit_path`, which takes the stage by value) from the second, so
        // neither can be used after its stage has ended (`x_units_stage` at dispatch point 3,
        // `x_seed_stage` at dispatch point 5). A stage that is dropped (the return just below)
        // ends there. The run stops here if tables were already released.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:503-505
        // `lookup_ref` is NCBI's `lookup_wrap`, the one table every subject is scanned with.
        let x_seed_stage = x_seed_batch.lookup.enter(s_idx);
        let x_units_stage = x_seed_batch.units.enter(s_idx);
        let lookup_ref: &crate::algorithm::tblastx::lookup::BlastAaLookupTable =
            x_seed_stage.tables();
        // A subject after a read of the database that failed is not searched; it leaves through
        // the stages above, as a subject too short to search does.
        if x_read_failure.skips(s_idx) {
            subject_emit.emit(st, s_idx);
            return;
        }
        // NCBI: Creates ewp (diagonal table) ONCE per SUBJECT via BlastExtendWordNew
        // (blast_engine.c:1002). Each subject sequence gets a fresh diagonal array.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:109-180 (BlastExtendWordNew) allocates with calloc,
        // which zeros all entries. The diag_offset is initialized to window.
        //
        // Reset diagonal state for each subject to match NCBI behavior:
        for d in st.diag_array.iter_mut() {
            *d = DiagStruct::default();
        }
        st.diag_offset = window;
        x_lap.lap(XStage::Enter);

        // The preliminary frames of the subject (`subject_prep::subject_frames`). A subject too
        // short to search gets none, and a subject whose read failed is not searched: either
        // leaves here.
        let subject_prep::SubjectFrames {
            x_letter_frames_later,
            mut x_kept,
            s_len,
            s_frames,
            s_frames_preliminary_owned,
            x_db_mask,
        } = match subject_prep::subject_frames(subjects, s_idx, &db_code, wordsize) {
            Ok(Some(made)) => made,
            Ok(None) => {
                subject_emit.emit(st, s_idx);
                return;
            }
            Err(error) => {
                x_read_failure.fail(s_idx, error);
                subject_emit.emit(st, s_idx);
                return;
            }
        };
        let s_frames_report = &s_frames;
        // The statistics the rest of this closure reads (`SubjectStatistics`), once per subject
        // that passed the length check above: the subject's own where NCBI computes the
        // statistics of each subject (`context_cutoffs::per_subject_cutoffs`), else the batch's.
        // An overflow of the search space fails the batch as a failed read does.
        let x_own_cutoffs = match context_cutoffs::per_subject_cutoffs(
            subject_option_spaces,
            contexts_ref,
            numbers,
            s_len as i64,
            evalue_threshold,
        ) {
            Ok(cutoffs) => cutoffs,
            Err(overflow) => {
                x_read_failure.fail(s_idx, overflow.into());
                subject_emit.emit(st, s_idx);
                return;
            }
        };
        let x_own_statistics = x_own_cutoffs
            .as_ref()
            .map(|own| batch_statistics.of_subject(own));
        let statistics = x_own_statistics.as_ref().unwrap_or(&batch_statistics);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1492-1497
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

        // Reusable buffer for scanned word offset pairs (NCBI `BlastOffsetPair`).
        let offset_pairs = &mut st.offset_pairs;

        // Diagonal tracking array for hit extension (NCBI `diag->hit_level_array`).
        let diag_array = &mut st.diag_array;

        // Diagonal offset for this subject (NCBI `diag->offset`, reset to `window` per subject).
        let mut diag_offset: i32 = st.diag_offset;

        // The subject's linking parameters (`subject_link::subject_linking_params`), computed once
        // before scanning; both linking passes read the same state.
        let subject_len_nucl = s_len as i64;
        let linking_params = subject_link::subject_linking_params(
            args.sum_stats(),
            avg_query_length,
            subject_len_nucl,
            link_db_length_nucl,
            statistics.cutoff_score_min,
            &linking_params_for_cutoff,
        );

        // NCBI: Each subject frame gets its own init_hitlist that is reset between frames.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:491
        // Resets `init_hitlist` via `BlastInitHitListReset` between frames.
        // After each frame, init_hsps are converted and merged into combined_ungapped_hits.

        // Statistics for HSP saving analysis (long sequences only)
        let is_long_sequence = subject_len_nucl > 600_000;
        let collect_hsp_saving_stats = is_long_sequence && (debug_hsp_saving || diag_enabled);
        let mut stats_hsp_saved = 0usize;
        let mut stats_hsp_filtered_by_cutoff = 0usize;
        let mut stats_hsp_filtered_by_reeval = 0usize;
        let stats_hsp_filtered_by_hsp_test = 0usize;
        let mut stats_score_distribution: Vec<i32> = Vec::new();

        // NCBI: Combined HSP list across all subject frames (`combined_hsp_list`).
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:438
        let mut combined_ungapped_hits: Vec<UngappedHit> = Vec::new();

        // EXPERIMENT (PAIRPAR), dispatch points 2 and 3 (per subject): the unit path, when the
        // batch has one, takes the subject and fills `combined_ungapped_hits`, and the frame loop
        // gets no frame; the units stage ends in it (`subject_prep::run_unit_path`).
        let s_frames_preliminary: &[QueryFrame] = subject_prep::run_unit_path(
            &subject_prep::UnitPathArgs {
                s_frames_preliminary,
                s_frames_listed,
                contexts_ref,
                cutoff_scores: statistics.cutoff_scores,
                x_dropoff_per_context: &x_dropoff_per_context,
                window,
                wordsize,
                s_idx,
                s_len,
                parallel_pool,
                x_board_batch: &x_board_batch,
                x_seed_end,
                x_pair_lut_is_stub,
                x_read_failure: &x_read_failure,
                subject_besthit: statistics.word_finder_cx.subject_besthit,
            },
            x_units_stage,
            &mut x_lap,
            &mut x_pipe_seed,
            &mut combined_ungapped_hits,
        );
        // Dispatch point 4 of the seed release (per subject): the frame loop below
        // indexes the worker's diagonal table through a raw pointer with `diag_mask`. A worker
        // of a batch that the unit path takes has a table of zero cells, and the worker that ended
        // the batch's last lookup stage has dropped its table. The loop must then run zero times;
        // anything else stops the run here.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-516
        assert!(
            s_frames_preliminary.is_empty() || diag_array.len() == diag_array_size as usize,
            "the reference frame loop would run without the worker's diagonal table (subject {s_idx})"
        );

        // The frame loop of the subject (`frame_chunks::search_subject_frames`): each preliminary
        // frame's chunks, scanned by the word finder (`word_finder::scan_chunk`) and merged into
        // the frame's list, which joins `combined_ungapped_hits`.
        frame_chunks::search_subject_frames(
            &statistics.word_finder_cx,
            &frame_chunks::SubjectScan {
                subjects,
                x_soft_subjects,
                x_db_mask: x_db_mask.as_deref(),
                use_serial_scan_chunks,
                use_parallel_chunks,
                num_threads,
                #[cfg(all(
                    feature = "parallel",
                    any(not(target_arch = "wasm32"), feature = "wasm-threads")
                ))]
                offset_pairs_len,
                lookup_ref,
                s_idx,
                s_len,
                s_frames_preliminary,
                s_frames_listed,
                is_long_sequence,
                collect_hsp_saving_stats,
            },
            offset_pairs,
            diag_array,
            &mut diag_offset,
            frame_chunks::SubjectSums {
                stats_hsp_saved: &mut stats_hsp_saved,
                stats_hsp_filtered_by_cutoff: &mut stats_hsp_filtered_by_cutoff,
                stats_score_distribution: &mut stats_score_distribution,
                combined_ungapped_hits: &mut combined_ungapped_hits,
            },
            sc,
        );
        x_lap.lap(XStage::FrameLoop);

        // Dispatch point 5 of the seed release, the batch pipe's seed-stage end and the helped
        // mode of the linking calls (`subject_link::end_subject_seed_stage`). The guard
        // `_x_helped` ends with the subject.
        let (x_link_on_pool, _x_helped) = subject_link::end_subject_seed_stage(
            x_seed_stage,
            &mut st.diag_array,
            &mut x_pipe_seed,
            combined_ungapped_hits.len(),
            use_parallel,
            x_board_batch.is_some(),
        );
        x_lap.lap(XStage::SeedRelease);

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:561-584
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

        // NCBI-style subject frame base offsets for sum-statistics linking
        // (`subject_link::subject_frame_bases`).
        let subject_frame_bases = subject_link::subject_frame_bases(s_frames_listed);
        // Nothing after this reads the preliminary frames (the seed stage has ended); they are
        // freed before the first linking call, where the heap peaks.
        drop(s_frames_preliminary_owned);

        // What the two linking calls and the stages between and after them read
        // (subject_link.rs).
        let link_args = subject_link::SubjectLinkArgs {
            linking_params_for_cutoff: &linking_params_for_cutoff,
            linking_params: linking_params.as_ref(),
            contexts_ref,
            subject_frame_bases: &subject_frame_bases,
            length_adj_per_context: statistics.length_adj_per_context,
            eff_searchsp_per_context: statistics.eff_searchsp_per_context,
            x_link_on_pool,
        };
        let stage_args = subject_link::SubjectStageArgs {
            contexts_ref,
            cutoff_scores: statistics.cutoff_scores,
            evalue_threshold,
            timing_enabled,
            reeval_ns: &reeval_ns,
            reeval_calls: &reeval_calls,
            linking_ns: &linking_ns,
            linking_calls: &linking_calls,
            diag_enabled,
            diagnostics: &diagnostics,
            debug_output_filter,
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-899
        // NCBI links and prelim-evalue-reaps the raw ungapped HSP list inside
        // s_BlastSearchEngineCore before the outer translated-subject path
        // reevaluates ambiguities and calls BLAST_LinkHsps again.
        //
        // A list of the one-hit word finder (`-window_size 0`) that holds an HSP of length 0 or
        // less links with NCBI's `Int4` ends, and those HSPs are dropped after it; without sum
        // statistics the HSPs get their own E-values (`subject_link::first_link_subject_hits`).
        x_link_timer.arm(LinkCall::First);
        let prelinked_ungapped_hits =
            subject_link::first_link_subject_hits(&link_args, combined_ungapped_hits, window);
        x_link_timer.disarm();
        x_lap.lap(XStage::Link1);
        // The sort after the first linking call and the reap by the preliminary E-value
        // (`subject_link::sort_and_reap_after_first_link`).
        let mut prelinked_ungapped_hits = subject_link::sort_and_reap_after_first_link(
            prelinked_ungapped_hits,
            &link_args,
            &stage_args,
            &mut x_lap,
        );
        // The re-evaluation's fetch of a subject with HSPs, which the negative seqid list may
        // exclude before its letters are read (`subject_loop::SubjectEmit::reevaluation_fetch_keeps`).
        if !prelinked_ungapped_hits.is_empty() && !subject_emit.reevaluation_fetch_keeps(s_idx) {
            prelinked_ungapped_hits.clear();
        }

        // A subject with ambiguous bases gets the frames of its letters here, after its
        // preliminary frames were freed, when HSPs are left for the re-evaluation (none are read
        // otherwise).
        let x_letter_frames: Vec<QueryFrame>;
        let (s_frames, s_frames_report): (&[QueryFrame], &[QueryFrame]) = if x_letter_frames_later {
            x_letter_frames = if prelinked_ungapped_hits.is_empty() {
                Vec::new()
            } else {
                match subjects.x_letter_frames(s_idx, &db_code, x_kept.take()) {
                    Ok(frames) => frames,
                    Err(error) => {
                        x_read_failure.fail(s_idx, error);
                        st.hits.clear();
                        subject_emit.emit(st, s_idx);
                        return;
                    }
                }
            };
            x_lap.lap(XStage::Translate);
            (&x_letter_frames, &x_letter_frames)
        } else {
            (&s_frames, s_frames_report)
        };
        // Without sum statistics, the HSPs of length 0 or less, which NCBI's re-evaluation
        // deletes without reading them (`subject_link::drop_unread_without_linking`).
        subject_link::drop_unread_without_linking(&link_args, window, &mut prelinked_ungapped_hits);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2699-2705
        // The re-evaluation of HSPs on a subject (its `-subject_loc` range) of fewer than 12 bases
        // reads a frame NCBI never translates (`subject_loop::SubjectEmit::refuse_untranslated_frame`):
        // the batch fails, and the subject goes on without HSPs to its end, which emits it.
        if !prelinked_ungapped_hits.is_empty()
            && crate::core::blast_hits::target_translation_never_made(s_len)
        {
            subject_emit.refuse_untranslated_frame(s_idx);
            prelinked_ungapped_hits.clear();
        }
        // The re-evaluation with the subject's letters (`subject_link::reevaluate_subject_hits`).
        let ungapped_hits = subject_link::reevaluate_subject_hits(
            sc,
            prelinked_ungapped_hits,
            (s_frames, s_frames_report),
            &stage_args,
            collect_hsp_saving_stats,
            &mut stats_hsp_filtered_by_reeval,
            &mut x_lap,
        );

        if !ungapped_hits.is_empty() {
            // The stages after the re-evaluation, up to the subject's final HSPs in `st.hits`
            // (`subject_report::report_subject_hits`).
            subject_report::report_subject_hits(
                sc,
                ungapped_hits,
                subject_debug::HspSavingStats {
                    collect_hsp_saving_stats,
                    stats_hsp_saved,
                    stats_hsp_filtered_by_cutoff,
                    stats_hsp_filtered_by_reeval,
                    stats_hsp_filtered_by_hsp_test,
                    stats_score_distribution: &mut stats_score_distribution,
                    subject_len_nucl,
                    cutoff_score_min: statistics.cutoff_score_min,
                },
                &subject_report::SubjectReportArgs {
                    link: &link_args,
                    stage: &stage_args,
                    x_link_timer: &x_link_timer,
                    s_frames,
                    s_frames_report,
                    s_len,
                },
                &mut x_lap,
                &mut st.hits,
            );
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1475
        bar.inc(1);
        st.diag_offset = diag_offset;
        subject_emit.emit(st, s_idx);
    };
    // No NCBI counterpart: LOSAT_TIMING, the clock before the subject loop (x_search_timing.rs).
    let x_t_subjects = x_times.now();

    // The subject loop: `process_subject` for every subject of the batch, on the pool's threads
    // or on this one (`subject_loop::run_subject_loop`).
    subject_loop::run_subject_loop(
        &subject_loop::SubjectLoop {
            use_parallel,
            use_serial_scan_chunks,
            subjects,
            #[cfg(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            ))]
            parallel_pool,
            offset_pairs_len,
            x_worker_diag_cells,
            window,
            sc,
            tx_opt: &tx_opt,
            #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
            x_board_batch: &x_board_batch,
            #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
            x_board_unit_path,
            x_seed_end,
            x_board,
        },
        process_subject,
    );

    // EXPERIMENT (LOSAT_X_TBXBOARD): the batch leaves the board once its subject loop has ended
    // (`XEndBatch` removes it again when the search returns, which then finds nothing).
    if let Some((board, batch)) = &x_board_batch {
        board.end_batch(batch);
    }
    // No NCBI counterpart: LOSAT_TIMING, the clock after the subject loop (x_search_timing.rs).
    let x_t_collect = x_times.now();
    // Close channel so the writer can exit.
    if let Some(tx) = tx_opt {
        drop(tx);
    }

    bar.finish();
    // The batch's hit lists in report order (`batch_finish::collect_hits`).
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let hits = batch_finish::collect_hits(writer, local_collector);
    #[cfg(not(all(feature = "parallel", not(target_arch = "wasm32"))))]
    let hits = batch_finish::collect_hits(local_collector);
    if diag_enabled {
        print_diagnostics_summary(&diagnostics);
    }

    // EXPERIMENT (LOSAT_X_TBXPIPE): the block's lines, printed here, or kept in the batch's result
    // for the batch loop when the pipe runs this search (`x_seed_end`).
    let mut x_timing = batch_finish::timing_lines(&batch_finish::TimingInputs {
        timing_enabled,
        t_search_start,
        scan_ns: &scan_ns,
        scan_calls: &scan_calls,
        ungapped_ns: &ungapped_ns,
        ungapped_calls: &ungapped_calls,
        reeval_ns: &reeval_ns,
        reeval_calls: &reeval_calls,
        linking_ns: &linking_ns,
        linking_calls: &linking_calls,
        t_read_queries,
        t_build_lookup,
        t_read_subjects,
        t_total,
        x_times: &x_times,
        x_t_subjects,
        x_t_collect,
        x_link_timer: &x_link_timer,
    });
    if let Some(error) = x_read_failure.take() {
        return Err(error);
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
        fetched_utf8: fetch_excluded
            .as_ref()
            .and_then(crate::algorithm::common::database::FetchExcluded::take_utf8),
    })
}

#[cfg(test)]
#[path = "run_impl_tests.rs"]
mod tests;
