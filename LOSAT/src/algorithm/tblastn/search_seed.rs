//! Internal TBLASTN protein lookup and translated-subject seed stage: the scan plan and the
//! per-unit scan of the preliminary search.

use crate::algorithm::blastp::encoding::{
    encode_protein_query_frame_with_seg, encode_protein_query_frame_with_seg_timed,
};
use crate::algorithm::common::translated_chunks::TranslatedChunk;
use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair;
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{XCompressedLookup, X_TBLASTN};
#[cfg(test)]
use crate::algorithm::tblastx::lookup::{build_ncbi_lookup, build_ncbi_lookup_for_profile};
#[cfg(test)]
use crate::algorithm::tblastx::translation::generate_frames;
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::config::ScoringMatrix;
#[cfg(test)]
use crate::stats::lookup_protein_params_ungapped;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::matrix::ncbistdaa;
use crate::utils::seg::SegParams;
use anyhow::Result;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:486-502;
// c++/src/algo/blast/api/blast_setup_cxx.cpp:655-657;
// c++/src/algo/blast/core/blast_filter.c:1241-1256;
// c++/src/algo/blast/core/blast_setup.c:614-638:
// query lowercase spans join the filtering locations; their complement
// defines lookup segments, and hard masking replaces those residues with X.
// Keep the unmasked copy for identity calculation after traceback.
pub(super) fn encode_tblastn_lookup_query(
    query: &[u8],
    seg: Option<&SegParams>,
    mask_lowercase: bool,
) -> QueryFrame {
    encode_tblastn_lookup_query_with(query, seg, mask_lowercase, false)
}

/// `encode_tblastn_lookup_query`; with LOSAT_TIMING the encoding and the SEG go to the common
/// line (utils/xtiming.rs). For the driver's query preparation on the calling thread only.
pub(super) fn encode_tblastn_lookup_query_timed(
    query: &[u8],
    seg: Option<&SegParams>,
    mask_lowercase: bool,
) -> QueryFrame {
    encode_tblastn_lookup_query_with(query, seg, mask_lowercase, true)
}

fn encode_tblastn_lookup_query_with(
    query: &[u8],
    seg: Option<&SegParams>,
    mask_lowercase: bool,
    timed: bool,
) -> QueryFrame {
    let mut frame = if timed {
        encode_protein_query_frame_with_seg_timed(query, seg)
    } else {
        encode_protein_query_frame_with_seg(query, seg)
    };
    if !mask_lowercase {
        return frame;
    }
    let mut start = None;
    for (index, &residue) in query.iter().enumerate() {
        if residue.is_ascii_lowercase() {
            start.get_or_insert(index);
        } else if let Some(left) = start.take() {
            frame.seg_masks.push((left, index));
        }
    }
    if let Some(left) = start {
        frame.seg_masks.push((left, query.len()));
    }
    if frame.seg_masks.is_empty() {
        return frame;
    }
    frame.seg_masks.sort_unstable_by_key(|&(left, _)| left);
    let mut merged: Vec<(usize, usize)> = Vec::new();
    for (left, right) in frame.seg_masks.drain(..) {
        if let Some(last) = merged.last_mut() {
            if left <= last.1 {
                last.1 = last.1.max(right);
                continue;
            }
        }
        merged.push((left, right));
    }
    frame.seg_masks = merged;
    if frame.aa_seq_nomask.is_none() {
        frame.aa_seq_nomask = Some(frame.aa_seq.clone());
    }
    for &(left, right) in &frame.seg_masks {
        frame.aa_seq[1 + left..1 + right].fill(ncbistdaa::X);
    }
    frame
}

/// The runs of lower-case letters of `query`, as half-open intervals, appended to `out`.
fn push_lowercase_runs(query: &[u8], out: &mut Vec<(usize, usize)>) {
    let mut start = None;
    for (index, &residue) in query.iter().enumerate() {
        if residue.is_ascii_lowercase() {
            start.get_or_insert(index);
        } else if let Some(left) = start.take() {
            out.push((left, index));
        }
    }
    if let Some(left) = start {
        out.push((left, query.len()));
    }
}

/// The SEG intervals of `query` and, with `mask_lowercase`, its lower-case runs (half-open,
/// unsorted, not joined), for the masks of the BLAST XML report, which joins them as NCBI's
/// `BlastSeqLocCombine` of `filter_out` with 0 does: overlapping intervals only. The search's list
/// (`encode_tblastn_lookup_query`) also joins adjacent ones; that changes no masked letter,
/// but NCBI's XML mask loop (`s_MaskQuerySeq`) leaves a letter unmasked where one interval
/// ends on the second-last column of an alignment segment and the next begins on its last.
pub(super) fn tblastn_report_mask_intervals(
    query: &[u8],
    seg: Option<&SegParams>,
    mask_lowercase: bool,
) -> Vec<(usize, usize)> {
    let mut intervals = encode_protein_query_frame_with_seg(query, seg).seg_masks;
    if mask_lowercase {
        push_lowercase_runs(query, &mut intervals);
    }
    intervals
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:141-150
// Pairs the query and subject seed offsets (`q_off`, `s_off`) of `qs_offsets`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:808-812:
#[cfg(test)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Seed {
    pub frame: i8,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:241-250:
    // backup.offset identifies this chunk before WordFinder returns local offsets.
    pub chunk_offset: u32,
    pub query_offset: u32,
    pub subject_offset: u32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775,805-812
// NCBI translates the subject into its six frames once, then calls WordFinder for each frame and
// chunk with a read-only lookup table. XScanPlan holds what those calls need (the lookup table and
// the translated frames with their chunks and scan ranges) so that every unit can be scanned on its
// own.
// The scan state of the preliminary search, built once.  NCBI c++/src/algo/blast/core/
// blast_engine.c:478-500, 804-841 calls WordFinder once per unit; the lookup table is read-only
// there.
pub(super) struct XScanPlan {
    pub lookup: XPlanLookup,
    pub subject: std::sync::Arc<XSubjectPlan>,
    pub pair_capacity: usize,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:93-108
/// The lookup table of a scan plan: the standard table (`eAaLookupTable`), or the compressed
/// table of word sizes 5 to 7 (`eCompressedAaLookupTable`, tblastn-fast).
pub(super) enum XPlanLookup {
    Aa(std::sync::Arc<crate::algorithm::tblastx::lookup::BlastAaLookupTable>),
    Compressed(std::sync::Arc<XCompressedLookup>),
}

/// The compressed lookup table of a query set (word sizes 5 to 7), kept for all its subjects as
/// `XQuerySide` keeps the standard table: the first subject that needs it builds it while the
/// others wait. A refusal is returned to that subject and leaves nothing behind.
#[derive(Default)]
pub(super) struct XCompressedSide {
    table: std::sync::OnceLock<std::sync::Arc<XCompressedLookup>>,
    /// Held while `table` is being built (on the builder's thread alone: no pool task waits on
    /// it while it runs others).
    in_progress: std::sync::Mutex<()>,
}

// No NCBI counterpart: builds the table once for the subjects of a query batch, as NCBI builds it
// once before them (blast_aux_priv.cpp:206-213); it does not change any value NCBI computes.
impl XCompressedSide {
    fn table(
        &self,
        build: impl FnOnce() -> Result<XCompressedLookup>,
    ) -> Result<std::sync::Arc<XCompressedLookup>> {
        if let Some(table) = self.table.get() {
            return Ok(std::sync::Arc::clone(table));
        }
        let _in_progress = self
            .in_progress
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(table) = self.table.get() {
            return Ok(std::sync::Arc::clone(table));
        }
        let built = std::sync::Arc::new(build()?);
        Ok(std::sync::Arc::clone(self.table.get_or_init(|| built)))
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:104-108
/// The compressed lookup table of the query set (word sizes 5 to 7, `threshold` in its units,
/// `stage_d_pipeline::lookup_threshold`): the masked query frames of the standard table
/// (`x_scan_plan`), their concatenation and lookup segments as the standard table takes them (the
/// coordinates of the word finder's query, `x_lookup_query_for_profile`), and
/// blastp-fast's builders of NCBI's one table of the query batch (`x_tbn_build_compressed_lookup`:
/// the flat table, or the reference one with `LOSAT_X_BPCLUT=0`) on this thread. A matrix other
/// than BLOSUM62 takes the same builders with its own words (`x_tbn_build_compressed_lookup_of`).
fn x_build_compressed(
    queries: &[&[u8]],
    seg: Option<&SegParams>,
    mask_lowercase: bool,
    threshold: i32,
    matrix: ScoringMatrix,
    word_length: usize,
    // The query set's frames as the table's query, made once for the batch (`x_scan_plan`).
    x_frames: Option<&[QueryFrame]>,
) -> Result<XCompressedLookup> {
    use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
        x_tbn_build_compressed_lookup, x_tbn_build_compressed_lookup_of,
    };
    let query_frames: Vec<_> = match x_frames {
        Some(frames) => frames.iter().map(|frame| vec![frame.clone()]).collect(),
        None => queries
            .iter()
            .map(|query| vec![encode_tblastn_lookup_query(query, seg, mask_lowercase)])
            .collect(),
    };
    let (query, locations) =
        crate::algorithm::tblastx::lookup::x_lookup_query_for_profile(&query_frames, word_length);
    let tables = [0..locations.len()];
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1303-1310
    // Dispatch point: the words of BLOSUM62 (the default) or of the search's other matrix.
    let built = if matrix == ScoringMatrix::Blosum62 {
        x_tbn_build_compressed_lookup(word_length, threshold, &query, &locations, &tables)
    } else {
        x_tbn_build_compressed_lookup_of(
            matrix,
            word_length,
            threshold,
            &query,
            &locations,
            &tables,
            1,
        )
    };
    built
        .ok_or_else(|| {
            anyhow::anyhow!(
                "-word_size {word_length} (the compressed lookup table) is not supported by LOSAT's TBLASTN"
            )
        })?
        .map_err(|_| {
            // NCBI builds the table of each query batch in its search (`CLocalBlast`), after it
            // has printed the reports of the batches before: the refusal of its heap's corruption
            // comes after them (`blastinput::app::ErrorAfterBatchReports`).
            anyhow::Error::new(crate::blastinput::app::ErrorAfterBatchReports::new(anyhow::anyhow!(
                "-threshold {} with the compressed lookup table of -word_size {} is not supported by LOSAT's TBLASTN: the neighboring words of a query batch fill NCBI BLAST+'s 1024 overflow banks, past which NCBI BLAST+ 2.17.0 writes and corrupts its heap",
                f64::from(threshold) / 100.0,
                word_length
            )))
        })
}

/// The compressed table of a scan plan (`x_scan_plan`, word sizes 5 to 7) and its longest chain,
/// built once for the query set (`XCompressedSide`). Cold and not inlined, so that the standard
/// table's plan keeps its code.
#[cold]
#[inline(never)]
fn x_compressed_plan_lookup(
    x_compressed: &XCompressedSide,
    queries: &[&[u8]],
    seg: Option<&SegParams>,
    mask_lowercase: bool,
    threshold: i32,
    matrix: ScoringMatrix,
    word_length: usize,
    x_frames: Option<&[QueryFrame]>,
) -> Result<(XPlanLookup, i32)> {
    let table = x_compressed.table(|| {
        x_build_compressed(
            queries,
            seg,
            mask_lowercase,
            threshold,
            matrix,
            word_length,
            x_frames,
        )
    })?;
    let longest_chain = match &*table {
        XCompressedLookup::Reference(table) => table.longest_chain,
        XCompressedLookup::Flat(table) => table.longest_chain,
    };
    Ok((XPlanLookup::Compressed(table), longest_chain))
}

/// The word length of a compressed table (`BlastCompressedAaLookupTable::word_length`).
fn x_compressed_word_length(table: &XCompressedLookup) -> i32 {
    match table {
        XCompressedLookup::Reference(table) => table.word_length,
        XCompressedLookup::Flat(table) => table.word_length,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775,805-812
// The six translated frames of one subject and the (frame, chunk, scan range) units in NCBI call
// order.
/// Graduated 2026-10-22 (was `LOSAT_X_TBNSSIDE`): the part of a scan plan that depends on the
/// subject only.
pub(super) struct XSubjectPlan {
    pub frames: Vec<QueryFrame>,
    /// (index into `frames`, chunk, scan ranges) in NCBI call order.
    pub units: Vec<(usize, TranslatedChunk, Vec<(i32, i32)>)>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:627-631
// NCBI frees a subject's translation after the subject; a dropped plan gives its frames to the
// next plan of the thread instead.
impl Drop for XSubjectPlan {
    fn drop(&mut self) {
        for frame in self.frames.drain(..) {
            super::x_subject_translate::x_reuse_frame(frame.aa_seq);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775,805-812
// NCBI translates a subject once per query batch (BLAST_GetAllTranslations in
// s_BlastSearchEngineCore). The translation depends only on the subject, the genetic code and the
// lowercase-masking option, so XSubjectSide keeps it for later batches. A reuse of the same value;
// no value is computed differently.
/// Graduated 2026-10-22 (was `LOSAT_X_TBNSSIDE`; its rule stays: plans are kept from
/// `KEEP_FROM_BATCHES` batches on, up to `KEEP_LIMIT`): the subject plans of one search, by
/// subject index.
///
/// A plan (the six translated frames and their chunks) is a function of the
/// subject, the genetic code and the lowercase-masking option. The search
/// consists of one preliminary stage per query batch, and each of them made
/// the plan of every subject again. The owner of the search creates one value
/// for its fixed list of subjects and options; the first batch that reaches a
/// subject computes its plan, as before, and the later ones read it. Plans are
/// kept only up to a total size (`keeps`: decided from the subject's length by
/// the first batch that reaches it), beyond which they are made per use as the
/// plans of a search without a side are: without residues when the subject is
/// long (LOSAT_X_TBNUNITXLATE, `x_scan_plan`), so that four batches over a
/// genome hold the frames of the chromosomes they are on no more than one does.
/// A search of fewer than `KEEP_FROM_BATCHES` batches gets no value
/// (`for_batches`): one batch has no later reader, and two or three give the
/// time their later readers would save for the memory of the plans. With
/// `LOSAT_X_TBNKEEPFEW` (`keep_few`) a search of more batches that searches them
/// one after another gets none either unless its plans hold no more than its
/// threads could have in flight (`kept_within_threads`): many subjects are
/// translated again for each batch, as NCBI translates them, and few long ones
/// are kept. A search that searches batches side by side keeps them
/// (`stage_d_pipeline::x_searched_side_by_side`).
///
/// Batches can reach a subject at the same time. Whether its plan is kept is
/// decided by one of them while the others wait, and its bytes are added to
/// `kept_bytes` once; the kept plan is then made by one of them while the
/// others wait, so a plan is made once: `kept_bytes` is always the size of the
/// plans the cells hold.
pub(super) struct XSubjectSide {
    gencode: u8,
    mask_lowercase: bool,
    /// One cell per subject.
    slots: Vec<XSubjectSlot>,
    kept_bytes: std::sync::atomic::AtomicUsize,
    /// `KEEP_LIMIT` (a smaller limit in tests).
    limit: usize,
}

/// The cell of one subject of `XSubjectSide`: unset until a batch reaches the
/// subject.
#[derive(Default)]
struct XSubjectSlot {
    /// Whether the plan is kept (its bytes reserved within the limit), decided
    /// by the first batch that reaches the subject (`XSubjectSide::keeps`), and the subject's
    /// length (`XSubjectSlot::KEPT` and the bits below it): the later batches do not read the
    /// length of a `-db` sequence (`XSubjectSide::length`).
    kept: std::sync::OnceLock<u32>,
    /// The kept plan, made by the first batch that reaches the subject, or
    /// `None` when it could not be prepared (every later use then makes its
    /// own).
    plan: std::sync::OnceLock<Option<std::sync::Arc<XSubjectPlan>>>,
}

impl XSubjectSlot {
    /// The bit of `kept` that keeps the plan; a length (NCBI's `Int4`) is below it.
    const KEPT: u32 = 1 << 31;
}

// No NCBI counterpart: the methods below keep, count and hand out the cached subject plans (a memo
// of x_build_subject); it does not change any value NCBI computes.
impl XSubjectSide {
    /// Translated residues kept at most (two bytes per subject base).
    pub(super) const KEEP_LIMIT: usize = 256 << 20;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:627-631
    // NCBI frees a subject's translation when the subject has been searched, in every query
    // batch. Below this number of batches the search does the same.
    /// Query batches (20,000 query residues each) from which a search keeps its
    /// subject plans. A kept plan saves the later batches the translation of
    /// its subject and holds the translated residues until the search ends.
    /// Timed on 28,036 subjects (40.6 Mb) at 4 and 8 threads: the plans hold
    /// 115 MiB, and a search of two, three or four batches without them took a
    /// median of 1.03 to 1.11 of the time with them. The saving does not start
    /// at this number; below it the search gives that time for the memory.
    pub(super) const KEEP_FROM_BATCHES: usize = 4;

    /// EXPERIMENT (LOSAT_X_TBNKEEPFEW): a switch of the strict set (`utils/xswitch.rs`), on
    /// unless it is `0`. With it a search of
    /// `KEEP_FROM_BATCHES` or more batches keeps its plans only when `kept_within_threads`
    /// holds or it searches batches side by side (`stage_d_pipeline::x_searched_side_by_side`).
    pub(super) fn keep_few() -> bool {
        use std::sync::OnceLock;
        static ON: OnceLock<bool> = OnceLock::new();
        *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNKEEPFEW"))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:627-631
    // NCBI frees each subject's translation after the subject, so a search holds the
    // translations of the subjects its threads are on. The kept plans are held where they are no
    // more than that.
    /// EXPERIMENT (LOSAT_X_TBNKEEPFEW): true when the plans of subjects of these lengths hold
    /// no more than `threads` plans of the longest one: no more than the threads of a search
    /// that keeps none can have in flight, each on one subject. A few long subjects (one
    /// genome, a few chromosomes) are kept, and batches searched at the same time share one
    /// plan of each instead of making one per batch. Many subjects searched one batch after
    /// another are not: without the plans the search held 115 MiB less on 28,036 16S subjects
    /// (40.6 Mb) and 309 MiB less on 57,775 zebrafish cDNAs (122 Mb), with five query batches,
    /// and took 0.980-0.994 of the time at 4, 8 and 16 threads (workstation, 2026-10-19:
    /// translating a packed subject again, `LOSAT_X_TBNXLATESIMD`, costs less than holding its
    /// plan; with the frame reuse 0.998 at one thread).
    pub(super) fn kept_within_threads(longest: usize, total: usize, threads: usize) -> bool {
        total <= longest.saturating_mul(threads.max(1))
    }

    pub(super) fn new(subjects: usize, gencode: u8, mask_lowercase: bool) -> Self {
        Self::with_limit(subjects, gencode, mask_lowercase, Self::KEEP_LIMIT)
    }

    /// `new`, with `limit` bytes of kept plans instead of `KEEP_LIMIT`.
    fn with_limit(subjects: usize, gencode: u8, mask_lowercase: bool, limit: usize) -> Self {
        Self {
            gencode,
            mask_lowercase,
            slots: (0..subjects).map(|_| XSubjectSlot::default()).collect(),
            kept_bytes: std::sync::atomic::AtomicUsize::new(0),
            limit,
        }
    }

    /// The subject plans of a search of `batches` query batches, or `None`
    /// when the search keeps none: `on` is false (`LOSAT_X_TBNKEEPFEW` keeps
    /// none, `kept_within_threads`), or there are
    /// fewer than `KEEP_FROM_BATCHES` batches. A kept plan is read by
    /// the batches after the first one that reaches its subject. A search of
    /// fewer batches makes its plans per use and drops them after the subject,
    /// like those `keeps` refuses.
    pub(super) fn for_batches(
        on: bool,
        batches: usize,
        subjects: usize,
        gencode: u8,
        mask_lowercase: bool,
    ) -> Option<Self> {
        // The own side is recorded where a plan is kept or refused (KEEP_LIMIT).
        if on && batches < Self::KEEP_FROM_BATCHES {
            crate::utils::xreach::x_reach!("LOSAT_X_TBNSSIDE", false);
        }
        (on && batches >= Self::KEEP_FROM_BATCHES)
            .then(|| Self::new(subjects, gencode, mask_lowercase))
    }

    /// The cell of subject `oid`, when the side is one of this genetic code
    /// and masking option.
    fn slot(&self, oid: usize, gencode: u8, mask_lowercase: bool) -> Option<&XSubjectSlot> {
        if self.gencode != gencode || self.mask_lowercase != mask_lowercase {
            return None;
        }
        self.slots.get(oid)
    }

    /// Translated residues of a plan: what `KEEP_LIMIT` counts
    /// (`plan_bytes_of`, from the length; this one checks it in tests).
    #[cfg(test)]
    fn plan_bytes(plan: &XSubjectPlan) -> usize {
        plan.frames.iter().map(|frame| frame.aa_seq.len()).sum()
    }

    /// `plan_bytes` of the plan of a subject of `length` bases, from the
    /// length alone: the frames `i + 1` and `-(i + 1)` (`i + 3 <= length`)
    /// hold `(length - i) / 3` residues and two sentinels each
    /// (`x_subject_translate::x_packed_frames`, `x_frame`).
    fn plan_bytes_of(length: usize) -> usize {
        (0..3)
            .filter(|&i| i + 3 <= length)
            .map(|i| ((length - i) / 3 + 2).saturating_mul(2))
            .fold(0usize, usize::saturating_add)
    }

    /// Whether the plan of the subject in `slot`, of `length` bases, is kept:
    /// decided by the first batch that reaches the subject, which reserves
    /// the plan's bytes within the limit (`plan_bytes_of`). A plan that does
    /// not fit is made per use, as the plans of a search without a side are
    /// (`x_scan_plan`; `x_lean_subject` reads the decision before a `-db`
    /// sequence is read for its units).
    fn keeps(&self, slot: &XSubjectSlot, length: usize) -> bool {
        let decided = slot.kept.get_or_init(|| {
            let kept = Self::reserve(&self.kept_bytes, Self::plan_bytes_of(length), self.limit);
            // The own side is recorded where a plan is kept or refused.
            crate::utils::xreach::x_reach!("LOSAT_X_TBNSSIDE", kept);
            // A `-db` sequence's length (NCBI's `Int4`) is below the bit; a longer record's is
            // kept saturated (`length` reads `-db` sequences only).
            let length = u32::try_from(length)
                .unwrap_or(u32::MAX)
                .min(XSubjectSlot::KEPT - 1);
            length | if kept { XSubjectSlot::KEPT } else { 0 }
        });
        decided & XSubjectSlot::KEPT != 0
    }

    /// The length of subject `oid` that an earlier search of it gave `keeps` (none before), for
    /// the search's genetic code and masking option: a `-db` sequence whose plan the side keeps
    /// is not read again, not even for its length.
    pub(super) fn length(&self, oid: usize, gencode: u8, mask_lowercase: bool) -> Option<usize> {
        let decided = *self.slot(oid, gencode, mask_lowercase)?.kept.get()?;
        Some((decided & !XSubjectSlot::KEPT) as usize)
    }

    /// True when the plans of subjects of these lengths are all kept: a plan
    /// has six frames of at most `length / 3` residues and two sentinels, so
    /// `2 * length + 12` bounds `plan_bytes`.
    #[cfg(any(feature = "parallel", test))]
    pub(super) fn all_fit(total: usize, subjects: usize) -> bool {
        total
            .saturating_mul(2)
            .saturating_add(subjects.saturating_mul(12))
            <= Self::KEEP_LIMIT
    }

    /// Size of the plans held, for tests.
    #[cfg(test)]
    fn kept_bytes(&self) -> usize {
        self.kept_bytes.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// Adds `bytes` to `kept` if the total stays within `limit`. The total
    /// only changes when the plan fits, so a plan that does not fit never
    /// makes another one look too large.
    fn reserve(kept: &std::sync::atomic::AtomicUsize, bytes: usize, limit: usize) -> bool {
        use std::sync::atomic::Ordering;
        kept.fetch_update(Ordering::Relaxed, Ordering::Relaxed, |held| {
            held.checked_add(bytes).filter(|&total| total <= limit)
        })
        .is_ok()
    }

    /// The kept plan of the subject in `slot` (`keeps`, of `length` bases):
    /// made by the first caller and read by the others. When it could not be
    /// made, the first caller gets the error, its bytes are given back, and
    /// every later caller makes its own.
    fn plan(
        &self,
        slot: &XSubjectSlot,
        length: usize,
        build: impl Fn() -> Result<std::sync::Arc<XSubjectPlan>>,
    ) -> Result<std::sync::Arc<XSubjectPlan>> {
        let mut failed = None;
        let kept = slot.plan.get_or_init(|| match build() {
            Ok(plan) => Some(plan),
            Err(error) => {
                self.kept_bytes.fetch_sub(
                    Self::plan_bytes_of(length),
                    std::sync::atomic::Ordering::Relaxed,
                );
                failed = Some(error);
                None
            }
        });
        if let Some(error) = failed {
            return Err(error);
        }
        match kept {
            Some(plan) => Ok(std::sync::Arc::clone(plan)),
            None => build(),
        }
    }
}

// No NCBI counterpart: the rule of LOSAT_X_TBNUNITXLATE in `x_scan_plan` (its packed arm) for a
// subject that is not read yet; it computes no value NCBI computes.
/// Whether `x_scan_plan` makes the plan of a packed subject of
/// `length` bases without residues (LOSAT_X_TBNUNITXLATE): `x_units_translate` and the plan is
/// made for its use alone (`x_own_plan`). The length is read first, so only a subject long enough
/// for a plan without residues decides here whether its plan is kept; the others decide in
/// `x_scan_plan`, after its check of the subject's length (a subject too short to be searched
/// charges nothing).
pub(super) fn x_lean_subject(
    x_subject: Option<(&XSubjectSide, usize)>,
    length: usize,
    db_gencode: u8,
    mask_lowercase: bool,
) -> bool {
    super::x_subject_translate::x_units_translate(length)
        && x_own_plan(x_subject, length, db_gencode, mask_lowercase)
}

// No NCBI counterpart: where a plan is made for one use (NCBI makes every plan so); it computes
// no value NCBI computes.
/// Whether the plan of subject `oid` of `x_subject`, of `length` bases, is made for its use
/// alone: the search keeps no plans (no side, or one of another genetic code or masking
/// option), or none of this subject (`XSubjectSide::keeps`, decided once per subject).
fn x_own_plan(
    x_subject: Option<(&XSubjectSide, usize)>,
    length: usize,
    db_gencode: u8,
    mask_lowercase: bool,
) -> bool {
    x_subject.is_none_or(|(side, oid)| {
        side.slot(oid, db_gencode, mask_lowercase)
            .is_none_or(|slot| !side.keeps(slot, length))
    })
}

/// No NCBI counterpart: the offset array of `LOSAT_X_TBNPAIRCAP` holds this many entries more
/// than the longest chain, eight times NCBI's OFFSET_ARRAY_SIZE (lookup_wrap.h:119). A scan call
/// that stops on a full array starts again in the next call at the word that did not fit
/// (`x_scan_byte.rs` builds that word's block of candidates again), which costs little per call
/// but adds up at NCBI's size: dmel_b4.faa x zebrafish cDNAs (`-db`) at one thread took 1.0041 of
/// the instructions with 4096 + the longest chain, 1.0016 with 8192 +, 1.0004 with 16384 + and
/// 1.0001 with 32768 + (against 1024 x the longest chain; D14 s1 c1, 2026-10-26). With
/// `LOSAT_X_TBNSCANRESUME` a stopped call keeps its block, and the units that take its scan have
/// NCBI's size (`OFFSET_ARRAY_SIZE`).
const X_PAIR_ROOM: usize = 8 * 4096;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/lookup_wrap.h:119
// NCBI reference (598d8ae6): c++/include/algo/blast/core/lookup_wrap.h:119-119
/// NCBI's OFFSET_ARRAY_SIZE: the offset array of a unit scanned with `LOSAT_X_TBNSCANRESUME` holds
/// this many entries more than the longest chain (GetOffsetArraySize, lookup_wrap.c:264-266).
const OFFSET_ARRAY_SIZE: usize = 4096;

// No NCBI counterpart: reads the switch once; it does not change any value NCBI computes.
/// Whether a plan's offset array holds `X_PAIR_ROOM` + the longest chain entries (at most
/// 1024 x the longest chain, the size without the switch) (`LOSAT_X_TBNPAIRCAP`, read once).
fn x_pair_cap_on() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNPAIRCAP"))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775,805-812
// x_scan_plan does the setup that precedes NCBI's WordFinder calls for one subject (translation,
// scan ranges) and attaches the lookup table of the query set.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:91-100
/// `None` when the subject is too short to be searched.
pub(super) fn x_scan_plan(
    queries: &[&[u8]],
    // EXPERIMENT (LOSAT_X_TBNSTORE): the subject's letters, or its packed store; or a `-db`
    // sequence, read here when its plan is made.
    subject: super::x_subject_store::XPrelimSubject<'_>,
    db_gencode: u8,
    seg: Option<&SegParams>,
    threshold: i32,
    mask_lowercase: bool,
    matrix: ScoringMatrix,
    word_length: usize,
    // Where the table of this query set is kept (`XQuerySide`).
    x_lookup: &std::sync::OnceLock<
        std::sync::Arc<crate::algorithm::tblastx::lookup::BlastAaLookupTable>,
    >,
    // Where the compressed table of this query set is kept (word sizes 5 to 7, `XQuerySide`).
    x_compressed: &XCompressedSide,
    // The search's subject plans and this subject's index among them.
    x_subject: Option<(&XSubjectSide, usize)>,
    // The query set's frames as the table's query (`encode_tblastn_lookup_query` of each query
    // with `seg` and `mask_lowercase`), made once for the batch (`XQuerySide`), when the caller
    // has them.
    x_frames: Option<&[QueryFrame]>,
) -> Result<Option<XScanPlan>> {
    // An unknown genetic code is refused before the subject length is looked at.
    GeneticCode::try_from_id(db_gencode).map_err(anyhow::Error::msg)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1318-1334,1429-1432
    let length = subject.len();
    if length < word_length * 3 + 2 {
        return Ok(None);
    }
    // Everything below that reads the subject.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775,805-812
    // Translation and chunking of the subject (the shared plan memoises them).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1071-1074,1080-1090
    // What the shared plan memoises: the subject becomes ncbi2na in one pass and each frame is read
    // from the codon table of its strand, as BLAST_GetAllTranslations does for the compressed
    // subject (tblastn/x_subject_translate.rs). Graduated 2026-10-22 (was LOSAT_X_TBNXLATE).
    // EXPERIMENT (LOSAT_X_TBNSTORE): a packed subject is translated from its packed bases
    // (x_packed_subject_plan); a stage that reads the case of a subject with lower-case letters
    // gets the letters instead. The letters of a subject searched from its letters (the records
    // `run_local` is given, or a set that keeps its letters for a byte without an ncbi4na code)
    // are packed for the plan (`x_plan_store`), and the plan is made from them as for a packed
    // subject.
    // The slot of the subject's kept plan when the search keeps them (the dispatch point below);
    // a plan made for this use only may have no residues (LOSAT_X_TBNUNITXLATE).
    let x_slot = x_subject.and_then(|(side, oid)| {
        side.slot(oid, db_gencode, mask_lowercase)
            .map(|slot| (side, slot))
    });
    // Whether the plan is made for this use alone: the search keeps no plans, or none of this
    // subject (`XSubjectSide::keeps`, decided once per subject; `x_lean_subject` read the same
    // decision before a `-db` sequence was read for its units).
    let x_own = x_slot.is_none_or(|(side, slot)| !side.keeps(slot, length));
    let x_build_subject = || -> Result<std::sync::Arc<XSubjectPlan>> {
        use super::x_subject_store::XPrelimSubject;
        use crate::sequence::packed_subject::SubjectNt;
        match subject {
            // The plan is made from the sequence's bytes, read
            // for it (a database sequence has no lower-case letters).
            // In a masked search (`-db_soft_mask`, `-db_hard_mask`; series NP, M4b PR 2) every
            // sequence's units follow from its ranges in bases (`XDbStream::mask_ranges`).
            XPrelimSubject::Stream(stream, oid, _) => stream.with_block(oid, |sequence| {
                let length = sequence.len();
                let ranges = stream.mask_ranges(oid, length)?;
                let db_mask =
                    ranges
                        .as_deref()
                        .map(|ranges| super::x_subject_translate::XDbMaskRanges {
                            ranges,
                            hard: stream.mask().is_some_and(|mask| mask.hard()),
                        });
                super::x_subject_translate::x_packed_bytes_plan(
                    sequence.packed(),
                    length,
                    db_gencode,
                    db_mask,
                )
            }),
            XPrelimSubject::Nt(SubjectNt::Letters(letters)) => {
                match super::x_subject_store::x_plan_store(letters, mask_lowercase) {
                    Some(store) => {
                        super::x_subject_translate::x_packed_subject_plan(&store, db_gencode)
                    }
                    None => super::x_subject_translate::x_subject_plan(
                        letters,
                        db_gencode,
                        mask_lowercase,
                    ),
                }
            }
            XPrelimSubject::Nt(SubjectNt::Packed(store)) => {
                anyhow::ensure!(
                    !mask_lowercase || store.lowercase().is_empty(),
                    "TBLASTN: a packed subject with lower-case letters searched with lowercase masking"
                );
                // EXPERIMENT (LOSAT_X_TBNUNITXLATE): the plan of a long subject made for this
                // use only (no side, or a plan the side does not keep) has no residues; each
                // unit translates those it reads.
                if x_own {
                    let lean = super::x_subject_translate::x_units_translate(store.len());
                    crate::utils::xreach::x_reach!("LOSAT_X_TBNUNITXLATE", lean);
                    if lean {
                        return super::x_subject_translate::x_packed_lean_plan(store, db_gencode);
                    }
                }
                super::x_subject_translate::x_packed_subject_plan(store, db_gencode)
            }
        }
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775,805-812
    // Dispatch point: with the search's subject plans the plan is the shared one (made by the
    // first batch that reaches the subject); a search of fewer than four batches makes it here,
    // and so does every use of a plan the side does not keep (`XSubjectSide::keeps`).
    let subject_plan = match x_slot {
        Some((side, slot)) if !x_own => side.plan(slot, length, x_build_subject)?,
        _ => x_build_subject()?,
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:91-100
    // The query-side lookup table: BlastAaLookupTableNew and BlastAaLookupIndexQuery on the masked
    // query. The table-only variant (`x_lookup_table_for_profile`, graduated 2026-10-22, was
    // LOSAT_X_LUTSPLIT) skips the per-context Karlin blocks, which the table does not read.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:227-245:
    let x_build_lookup = || {
        let query_frames: Vec<_> = match x_frames {
            Some(frames) => frames.iter().map(|frame| vec![frame.clone()]).collect(),
            None => queries
                .iter()
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1241-1255 merges lowercase query and SEG masks.
                .map(|query| vec![encode_tblastn_lookup_query(query, seg, mask_lowercase)])
                .collect(),
        };
        std::sync::Arc::new(
            crate::algorithm::tblastx::lookup::x_lookup_table_for_profile(
                &query_frames,
                threshold,
                matrix,
                word_length,
            ),
        )
    };
    // The table depends on the query set and the word-finder options only; the first subject that
    // gets here builds it (`XQuerySide`, was LOSAT_X_TBNQSIDE), unless the batch built it before
    // its subjects (`XQuerySide::x_build_lookup_on`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:91-100
    // NCBI builds the lookup table once per query batch: the first subject builds it and the
    // others read it.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:612-619
    // Dispatch point: word sizes 5 to 7 take the compressed table (`x_build_compressed`).
    let (lookup, longest_chain) = if word_length > 4 {
        x_compressed_plan_lookup(
            x_compressed,
            queries,
            seg,
            mask_lowercase,
            threshold,
            matrix,
            word_length,
            x_frames,
        )?
    } else {
        let lookup = std::sync::Arc::clone(x_lookup.get_or_init(x_build_lookup));
        let longest_chain = lookup.longest_chain;
        (XPlanLookup::Aa(lookup), longest_chain)
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:264-270
    // blast_engine.c:992,1040-1041 allocates each thread's offset array of that size and
    // aa_ungapped.c:504-505 passes it to scansub (OFFSET_ARRAY_SIZE 4096, lookup_wrap.h:119).
    // The size only decides where a scan call stops; the next call resumes there, so the pairs
    // reach the caller in the same order at any size. With `LOSAT_X_TBNPAIRCAP` the array has
    // NCBI's form with `X_PAIR_ROOM` for OFFSET_ARRAY_SIZE, never more entries than without it
    // (1024 x the longest chain, which a table whose longest chain is below 33 keeps).
    let longest_chain = longest_chain.max(1) as usize;
    let pair_capacity = if x_pair_cap_on() {
        (X_PAIR_ROOM + longest_chain).min(longest_chain.saturating_mul(1024))
    } else {
        longest_chain * 1024
    };
    i32::try_from(pair_capacity).expect("NCBI offset array fits Int4");
    Ok(Some(XScanPlan {
        lookup,
        subject: subject_plan,
        pair_capacity,
    }))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-515
// x_scan_unit runs the scansub() calls of one unit and passes every filled offset array to the
// caller in scan order, as the loop at line 506 consumes it.
/// Scan one unit and hand every filled offset
/// array to `on_pairs` in scan order, as NCBI c++/src/algo/blast/core/
/// aa_ungapped.c:492-516 consumes it (a loop of `i` from 0 while it is below `hits`).
/// `subject` holds the unit's frame from `lead` residues before the chunk's first one on (the
/// frame's `aa_seq` from `1 + chunk.offset - lead`, or the residues the unit translated,
/// `x_unit_residues`), `lead` being `common::translated_chunks::scan_lead` of the unit's ranges
/// (0 for ranges in order); the pairs are at the chunk's offsets.
#[inline(always)]
pub(super) fn x_scan_unit(
    plan: &XScanPlan,
    unit: usize,
    subject: &[u8],
    lead: usize,
    pairs: &mut Vec<BlastOffsetPair>,
    on_pairs: impl FnMut(&[BlastOffsetPair]),
) {
    x_scan_unit_of::<true>(plan, unit, subject, lead, pairs, on_pairs)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:494-502 and 760-764
/// `x_scan_unit` of the two-hit word finder (`CLAMP`: a first range shorter than a word scans
/// the word at its left end, aa_ungapped.c:499-500) or of the one-hit word finder (`-window_size
/// 0`, `s_BlastAaWordFinder_OneHit`: no such clamp, so a unit whose first range is shorter than
/// a word is not scanned at all).
pub(super) fn x_scan_unit_of<const CLAMP: bool>(
    plan: &XScanPlan,
    unit: usize,
    subject: &[u8],
    lead: usize,
    pairs: &mut Vec<BlastOffsetPair>,
    mut on_pairs: impl FnMut(&[BlastOffsetPair]),
) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:300-306 and
    // c++/src/algo/blast/core/blast_aascan.c:76-77, :110
    // A later range of masks stored out of order can start before the chunk: NCBI scans the
    // frame's letters there and keeps the hits at the chunk's offsets (below 0, as Uint4). The
    // scan reads `subject` from `lead` residues before the chunk on the ranges moved by `lead`,
    // and the pairs are moved back (`x_pairs_at_chunk`) before `on_pairs` takes them.
    debug_assert_eq!(
        lead,
        crate::algorithm::common::translated_chunks::scan_lead(&plan.subject.units[unit].2)
    );
    let lead_i32 = lead as i32;
    let mut on_pairs = |pairs: &mut [BlastOffsetPair]| {
        if lead > 0 {
            crate::algorithm::tblastx::blast_aascan::x_pairs_at_chunk(pairs, lead);
        }
        on_pairs(&*pairs)
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1018-1019
    // The scan variant of the unit is chosen once per table, as NCBI picks one scansub: the
    // compressed table's (`x_scan_compressed_call`, range by range as below), or one of the
    // standard table's below.
    let lookup = match &plan.lookup {
        XPlanLookup::Aa(lookup) => lookup,
        XPlanLookup::Compressed(table) => {
            let (_, _, ranges) = &plan.subject.units[unit];
            if pairs.len() < plan.pair_capacity {
                pairs.resize(plan.pair_capacity, BlastOffsetPair::default());
            }
            let pair_capacity = plan.pair_capacity as i32;
            let word_length = x_compressed_word_length(table);
            for (range_index, &(left, right)) in ranges.iter().enumerate() {
                let (left, right) = (left + lead_i32, right + lead_i32);
                let mut end = right - word_length;
                if range_index == 0 && end < left {
                    if !CLAMP {
                        break;
                    }
                    end = left;
                }
                let mut scan_range = [0, left, end];
                while scan_range[1] <= scan_range[2] {
                    let hits = x_scan_compressed_call(
                        table,
                        subject,
                        (left, right),
                        pairs,
                        pair_capacity,
                        &mut scan_range,
                    );
                    on_pairs(&mut pairs[..hits as usize]);
                }
            }
            return;
        }
    };
    let (_, _, ranges) = &plan.subject.units[unit];
    //
    // EXPERIMENT (LOSAT_X_TBNSCANBYTE): the scan on presence bytes, for a table of tblastn's shape
    // (without AVX2 a sparse one); it writes up to X_SCAN_BYTE_SLACK entries after the capacity.
    let scan_byte = crate::algorithm::tblastx::x_scan_byte::x_scan_byte();
    let presence = crate::algorithm::tblastx::x_scan_byte::x_presence(scan_byte, lookup);
    // EXPERIMENT (LOSAT_X_TBNSCANRESUME): a call of the scan on presence bytes that stops on a full
    // offset array leaves its block in `held`, and the next call copies on from there instead of
    // making the block's candidates again; the array then has NCBI's size (OFFSET_ARRAY_SIZE + the
    // longest chain, lookup_wrap.c:264-266), never more entries than the plan's.
    let resume = presence.is_some() && crate::algorithm::tblastx::x_scan_byte::x_scan_resume();
    let capacity = if resume {
        let longest_chain = lookup.longest_chain.max(1) as usize;
        (OFFSET_ARRAY_SIZE + longest_chain).min(plan.pair_capacity)
    } else {
        plan.pair_capacity
    };
    // The branch-free scan writes up to X_SCAN_BF_SLACK entries after the capacity.
    let pairs_len = if presence.is_some() {
        capacity + crate::algorithm::tblastx::x_scan_byte::X_SCAN_BYTE_SLACK
    } else {
        capacity + crate::algorithm::tblastx::x_aascan_bf::X_SCAN_BF_SLACK
    };
    if pairs.len() < pairs_len {
        pairs.resize(pairs_len, BlastOffsetPair::default());
    }
    let pair_capacity = capacity as i32;
    let mut held = crate::algorithm::tblastx::x_scan_byte::XScanHeld::new();
    // The scan calls of the unit, for LOSAT_X_STATS.
    let mut x_calls = 0u64;
    for (range_index, &(left, right)) in ranges.iter().enumerate() {
        let (left, right) = (left + lead_i32, right + lead_i32);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:496-501
        let mut end = right - lookup.word_length;
        if range_index == 0 && end < left {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:762-765 and
            // c++/src/algo/blast/core/masksubj.inl:58-59
            // The one-hit word finder's loop does not start on a first range shorter than a
            // word, so `scansub`, which would step to the next ranges, is never called: the unit
            // has no hit.
            if !CLAMP {
                break;
            }
            end = left;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1098-1101 and
            // c++/src/algo/blast/core/blast_stat.c:1559-1575
            // A raised first range on the last letter of a frame's last chunk (a soft-masked
            // subject split into chunks) gives a word past `subject`, which ends at the frame's
            // closing NULLB. NCBI's word there holds that NULLB and the next byte of its buffer of
            // the six frames. NULLB is NCBIstdaa '-', whose scores are BLAST_SCORE_MIN, so no
            // neighbor word holds it, and no query word does (NCBI's FASTA reader drops hyphens;
            // TBLASTN takes no PSSM): NCBI's scan finds no hit, and the word is not scanned here
            // (TBLASTX does the same, tblastx/blast_engine/word_finder.rs).
            if left as usize + lookup.word_length as usize > subject.len() {
                continue;
            }
        }
        let mut scan_range = [0, left, end];
        while scan_range[1] <= scan_range[2] {
            let hits = if let Some(presence) = &presence {
                if resume {
                    crate::algorithm::tblastx::x_scan_byte::x_scan_one_range_held(
                        scan_byte,
                        lookup,
                        presence,
                        subject,
                        pairs,
                        pair_capacity,
                        &mut scan_range,
                        &mut held,
                    )
                } else {
                    crate::algorithm::tblastx::x_scan_byte::x_scan_one_range_byte(
                        scan_byte,
                        lookup,
                        presence,
                        subject,
                        pairs,
                        pair_capacity,
                        &mut scan_range,
                    )
                }
            } else {
                crate::algorithm::tblastx::x_aascan_bf::x_scan_one_range_bf(
                    lookup,
                    subject,
                    pairs,
                    pair_capacity,
                    &mut scan_range,
                )
            };
            on_pairs(&mut pairs[..hits as usize]);
            x_calls += 1;
        }
    }
    if presence.is_some() {
        crate::algorithm::tblastx::x_scan_byte::x_count_unit(x_calls);
    } else if crate::algorithm::tblastx::x_aascan_bf::x_stats() {
        crate::algorithm::tblastx::x_aascan_bf::x_count_unit(x_calls);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:487-500 and
// c++/src/algo/blast/core/blast_aascan.c:564-568
/// One scansub call of `s_BlastCompressedAaScanSubject` (the table's `scan_subject`, which writes
/// no entry past `array_size`) on the range `(left, right)` from `scan_range`, for the compressed
/// arm of `x_scan_unit` (`s_DetermineScanningOffsets` of a list of one range ends after it, and a
/// call that stops at a full array resumes in it). Not inlined, and it takes no closure, so that
/// the standard table's scan keeps its code.
#[inline(never)]
fn x_scan_compressed_call(
    table: &XCompressedLookup,
    subject: &[u8],
    (left, right): (i32, i32),
    pairs: &mut [BlastOffsetPair],
    pair_capacity: i32,
    scan_range: &mut [i32; 3],
) -> i32 {
    let range = [(left, right)];
    match table {
        XCompressedLookup::Reference(table) => {
            table.scan_subject(subject, &range, pairs, pair_capacity, scan_range)
        }
        XCompressedLookup::Flat(table) => {
            table.scan_subject::<X_TBLASTN>(subject, &range, pairs, pair_capacity, scan_range)
        }
    }
}

#[cfg(test)]
#[path = "search_seed_tests.rs"]
mod tests;
