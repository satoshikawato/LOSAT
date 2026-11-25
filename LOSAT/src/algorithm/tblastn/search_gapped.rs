//! Internal TBLASTN protein gapped stage of the preliminary search: the gapped score of each
//! (frame, chunk) unit, the endpoint purge, the merge and the append, and the traceback helpers.

use super::search_init::InitHsp;
use super::search_seed::encode_tblastn_lookup_query;
#[cfg(test)]
use super::search_seed::Seed;
use crate::algorithm::blastn::interval_tree::{
    x_grid_cells, BlastIntervalTree, IndexMethod, TreeHsp,
};
use crate::algorithm::blastp::encoding::encode_protein_query_frame_with_seg;
use crate::algorithm::blastp::gapalign::{
    gapped_start_window, protein_identities_from_edit_ops, GapAlignScratch,
};
use crate::algorithm::common::translated_chunks::{frame_length, ChunkLetters, TranslatedChunk};
use crate::algorithm::tblastx::translation::{
    translate_frame_into, x_codon_residues, x_codons_into, QueryFrame,
};
use crate::common::GapEditOp;
use crate::config::ScoringMatrix;
use crate::sequence::packed_subject::SubjectNt;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::matrix::protein_score;
use crate::utils::seg::SegParams;
use anyhow::{bail, Context, Result};

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:96-103,126-143:
// Mirrors `BlastSeg` (frame, offset, end, and gapped start) and `BlastHSP` (score with
// query and subject segments).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct GappedHsp {
    pub frame: i8,
    pub score: i32,
    pub q_start: i32,
    pub q_end: i32,
    pub q_gapped_start: i32,
    pub s_start: i32,
    pub s_end: i32,
    pub s_gapped_start: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-536,583-605;
// c++/include/algo/blast/core/blast_hits.h:126-143:
// Keep the transferred script on the same Rust HSP through purge and pruning.
#[derive(Clone, Debug)]
pub(super) struct TracebackOwnedHsp {
    pub hsp: GappedHsp,
    pub retried: bool,
    pub edit_script: Vec<GapEditOp>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-552:
// Score only this chunk's sorted initial HSPs before the next chunk scan.
/// `subject_sequence` is the chunk's residues of `frame` and `subject_frame` the frame's residues
/// from the chunk's first one on (`subject_sequence` and the letters after it); `score` is the
/// source of the scores (the matrix, or a PSSM query's rows, `score::TbnScore`); `one_hit` the
/// buffers that the one-hit word finder's HSPs of negative length read (`-window_size 0`, the
/// instance `ONE_HIT`; the other one is the two-hit units' code, which reads neither).
#[allow(clippy::too_many_arguments)]
fn score_gapped_chunk<S: super::score::TbnScore, const ONE_HIT: bool>(
    queries: &[&[u8]],
    query_frames: &[QueryFrame],
    context_offsets: &[i32],
    total_query_length: i32,
    frame: &QueryFrame,
    chunk: TranslatedChunk,
    subject_sequence: &[u8],
    subject_frame: &[u8],
    initial: &[InitHsp],
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    cutoff_scores: &[i32],
    matrix: ScoringMatrix,
    scratch: &mut S::Scratch,
    score: S,
    one_hit: Option<&XOneHitBuffers>,
) -> Result<Vec<(usize, GappedHsp)>> {
    let mut results = Vec::new();
    debug_assert_eq!(subject_sequence.len(), chunk.length);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:223-235:
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:273-313:
    // compare score DESC, subject start ASC, length DESC, query start ASC.
    let mut frame_initial: Vec<_> = initial
        .iter()
        .filter(|hit| hit.frame == frame.frame && hit.chunk_offset == chunk.offset as u32)
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:306-310
    // Stable, as glibc's qsort under the pinned NCBI BLAST+ (TN-5).
    frame_initial.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.s_start.cmp(&b.s_start))
            .then(b.length.cmp(&a.length))
            .then(a.q_start.cmp(&b.q_start))
    });
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:500:
    // Skip `GetGappedScore` when the initial hit list is empty.
    if frame_initial.is_empty() {
        return Ok(Vec::new());
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3827-3834,3908-3919,
    // 4065-4091:
    let mut tree = BlastIntervalTree::new(
        0,
        total_query_length + 1,
        0,
        i32::try_from(chunk.length)? + 1,
    );
    tree.reset_with_bounds_cells(
        0,
        total_query_length + 1,
        0,
        i32::try_from(chunk.length)? + 1,
        x_grid_cells(frame_initial.len()),
    );
    for hit in frame_initial {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2371-2468:
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:212-229:
        // `s_GetQueryStrandOffset` returns the protein context's `query_offset`,
        // which is also used by containment and insertion.
        let context = context_offsets.partition_point(|&offset| offset <= hit.q_seed as i32) - 1;
        let context_offset = context_offsets[context];
        let query_sequence = &query_frames[context].aa_seq[1..1 + queries[context].len()];
        let initial_tree_hsp = TreeHsp {
            query_offset: hit.q_start - context_offset,
            query_end: hit.q_start - context_offset + hit.length,
            subject_offset: hit.s_start,
            subject_end: hit.s_start + hit.length,
            score: hit.score,
            query_frame: 0,
            query_length: i32::try_from(queries[context].len())?,
            query_context_offset: context_offset,
            subject_frame_sign: i32::from(frame.frame.signum()),
        };
        if tree.contains_hsp(&initial_tree_hsp, context_offset, 0) {
            continue;
        }
        let q_start = usize::try_from(hit.q_start - context_offset)
            .context("negative context-local initial query start")?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3393-3404,3936-3953:
        // Selects the gapped alignment start and runs `s_BlastProtGappedAlignment`.
        // An initial alignment starting before its chunk (`-db_soft_mask` with masks stored out of
        // order) reaches here only with length in 0 to HSP_MAX_WINDOW and its middle in the chunk
        // (`x_crashes_ncbi` refused the rest): NCBI sets its gapped start to its middle without
        // reading a letter. A negative Int4 length comes only from the one-hit word finder
        // (`-window_size 0`, aa_ungapped.c:1054,1083); NCBI reads it as Uint4
        // (`gapped_start_int4`), and its window may read past the context and the chunk
        // (`x_gapped_start_negative`).
        let q_gapped_start = if hit.s_start < 0 {
            x_gapped_start_before_chunk_middle(q_start, hit)?
        } else {
            let s_start = hit.s_start as usize;
            if ONE_HIT && hit.length < 0 {
                x_gapped_start_negative(
                    score,
                    one_hit,
                    query_sequence,
                    frame,
                    chunk,
                    subject_frame,
                    q_start,
                    s_start,
                    hit.length,
                    matrix,
                )?
            } else {
                let length = usize::try_from(hit.length).context("negative initial HSP length")?;
                score.gapped_start(
                    query_sequence,
                    subject_sequence,
                    q_start,
                    length,
                    s_start,
                    matrix,
                )
            }
        };
        let s_gapped_start = x_gapped_subject_start(
            hit.s_seed,
            i32::try_from(q_gapped_start)?,
            i32::try_from(hit.q_seed)? - context_offset,
        )?;
        let gapped = if ONE_HIT && s_gapped_start >= chunk.length {
            // Only the window of a negative length reaches past the chunk (`q_gapped_start` is
            // then `q_start + 5`: inside the frame's letters, or with a PSSM, for a chunk that ends
            // within two letters of the frame's end, on the frame's NULLB or up to two letters
            // after it).
            x_score_only_past_chunk(
                score,
                query_sequence,
                frame,
                chunk,
                subject_frame,
                one_hit,
                q_gapped_start,
                s_gapped_start,
                matrix,
                gap_open,
                gap_extend,
                x_drop,
                scratch,
            )?
        } else {
            score.score_only(
                query_sequence,
                subject_sequence,
                q_gapped_start,
                s_gapped_start,
                matrix,
                gap_open,
                gap_extend,
                x_drop,
                scratch,
            )
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3924-3927,4058-4091:
        if gapped.score >= cutoff_scores[context] {
            let saved = GappedHsp {
                frame: frame.frame,
                score: gapped.score,
                q_start: gapped.query_start,
                q_end: gapped.query_stop,
                q_gapped_start: i32::try_from(q_gapped_start)?,
                s_start: gapped.subject_start,
                s_end: gapped.subject_stop,
                s_gapped_start: i32::try_from(s_gapped_start)?,
            };
            let saved_tree_hsp = TreeHsp {
                query_offset: saved.q_start,
                query_end: saved.q_end,
                subject_offset: saved.s_start,
                subject_end: saved.s_end,
                score: saved.score,
                query_frame: 0,
                query_length: i32::try_from(queries[context].len())?,
                query_context_offset: context_offset,
                subject_frame_sign: i32::from(frame.frame.signum()),
            };
            if ONE_HIT {
                x_add_one_hit_hsp(&mut tree, saved_tree_hsp, context_offset)?;
            } else {
                tree.add_hsp(saved_tree_hsp, context_offset, IndexMethod::QueryAndSubject);
            }
            results.push((context, saved));
        }
    }
    Ok(results)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:500-501,654,
// c++/src/algo/blast/core/blast_query_info.c:246-250 and
// c++/src/algo/blast/core/blast_util.c:1058,1063-1064,1080,1098-1101
/// The subject's side of NCBI's buffers that the gapped start of a one-hit HSP of negative length
/// reads (`-window_size 0`; `BlastGetStartForGappedAlignment` with the length read as Uint4,
/// blast_gapalign.c:3405-3413). Its window of HSP_MAX_WINDOW letters can pass the end of its
/// frame into the next frames of the subject's translation (`BLAST_GetAllTranslations`: frames 1,
/// 2, 3, -1, -2, -3, each after a NULLB, and a NULLB after the last), whose letters a PSSM's '-'
/// column can make count (`TbnScore::GAPPED_START_PAST_FRAME`). The translation buffer is set up
/// to the NULLB after frame -3 (the byte after it is allocated and never written); beyond it NCBI
/// reads a byte it never set or memory that is not its own, and can crash, and the window reads
/// every letter there as the sentinel NULLB (PD-LOSAT-NCBI-DEFECTS exception 3, which version 1.48
/// extends from BLASTP to TBLASTN). The query's side is the context's letters and NULLBs past them
/// (`gapped_start_window`, as BLASTP and BLASTX read it): with the matrix a window that passes its
/// context holds the context's NULLB, which scores BLAST_SCORE_MIN, so the bytes after it never
/// change the start; a PSSM's window reads its rows, not the query.
pub(super) struct XOneHitBuffers {
    /// `frame_offsets` of the translation buffer: the NULLB before each frame, and the NULLB after
    /// the last (index 6).
    frame_offsets: [usize; 7],
    /// The unit's frame, as an index of `frame_offsets`.
    frame_index: usize,
    /// For a chunk whose windows can pass its frame's end, the bytes of the translation buffer
    /// after the frame's NULLB (the next frames' letters and NULLBs, up to the buffer's last
    /// NULLB), as many as a window reads; empty otherwise.
    after_frame: Vec<u8>,
    /// For a gapped start on the frame's NULLB or after it (`x_score_only_past_chunk`), the bytes of
    /// NCBI's buffer from the chunk's first letter on: the frame's letters, its NULLB and
    /// `after_frame`; made once for the unit, at its first such start.
    frame_end: std::cell::OnceCell<Vec<u8>>,
}

impl XOneHitBuffers {
    /// The buffers of the unit of `chunk` of frame `frame_index` of `frames` (the subject's six
    /// frames in NCBI's order, `subject` their source).
    fn new(
        frames: &[QueryFrame],
        frame_index: usize,
        chunk: TranslatedChunk,
        subject: super::x_subject_store::XPrelimSubject<'_>,
        db_gencode: u8,
    ) -> Result<Self> {
        const NCBI_FRAMES: [i8; 6] = [1, 2, 3, -1, -2, -3];
        // Every frame of a subject of 5 bases or more is made, in this order (`x_generate_frames`,
        // `x_packed_lean_plan`); a shorter subject has no word.
        anyhow::ensure!(
            frames.len() == 6 && frames.iter().zip(NCBI_FRAMES).all(|(f, n)| f.frame == n),
            "TBLASTN -window_size 0: the subject's frames are not NCBI's six (internal error)"
        );
        let mut frame_offsets = [0usize; 7];
        for (index, frame) in frames.iter().enumerate() {
            frame_offsets[index + 1] = frame_offsets[index] + frame.aa_len + 1;
        }
        let frame = &frames[frame_index];
        let mut after_frame = Vec::new();
        // A window starts inside the chunk and reads HSP_MAX_WINDOW letters: it passes the frame's
        // end only from a chunk that ends fewer than HSP_MAX_WINDOW letters before it (the frame's
        // last chunk, or one that the frame goes on past by a few letters).
        if chunk.offset + chunk.length + HSP_MAX_WINDOW - 1 > frame.aa_len {
            let mut held = Vec::new();
            for next in &frames[frame_index + 1..] {
                let wanted = (HSP_MAX_WINDOW - 1).saturating_sub(after_frame.len());
                if wanted == 0 {
                    break;
                }
                let n = next.aa_len.min(wanted);
                if next.aa_seq.is_empty() {
                    // A frame of a plan without residues (LOSAT_X_TBNUNITXLATE): its first
                    // letters, translated as its units translate theirs.
                    let super::x_subject_store::XPrelimSubject::Nt(SubjectNt::Packed(store)) =
                        subject
                    else {
                        bail!(
                            "TBLASTN -window_size 0: a frame without residues of a subject that is not packed (internal error)"
                        );
                    };
                    let range = super::x_subject_translate::x_unit_residues(
                        store,
                        db_gencode,
                        next,
                        TranslatedChunk::new(0, n),
                        &mut held,
                    )?;
                    after_frame.extend_from_slice(&held[range][..n]);
                } else {
                    after_frame.extend_from_slice(&next.aa_seq[1..1 + n]);
                }
                if n == next.aa_len {
                    // The NULLB after the frame.
                    after_frame.push(0);
                }
            }
        }
        Ok(Self {
            frame_offsets,
            frame_index,
            after_frame,
            frame_end: std::cell::OnceCell::new(),
        })
    }

    /// The byte of the translation buffer at letter `position` of the unit's `frame` (its NULLB at
    /// `aa_len`, then the bytes after it), where `subject_frame` holds the frame's letters from
    /// the first letter of `chunk` on.
    fn frame_byte(
        &self,
        frame: &QueryFrame,
        chunk: TranslatedChunk,
        subject_frame: &[u8],
        position: usize,
    ) -> Result<u8> {
        let byte = match position.cmp(&frame.aa_len) {
            std::cmp::Ordering::Less => position
                .checked_sub(chunk.offset)
                .and_then(|index| subject_frame.get(index)),
            std::cmp::Ordering::Equal => Some(&0),
            std::cmp::Ordering::Greater => self.after_frame.get(position - frame.aa_len - 1),
        };
        byte.copied().with_context(|| {
            format!(
                "TBLASTN -window_size 0: letter {} of a frame of {} letters is not held (internal error)",
                position, frame.aa_len
            )
        })
    }

    /// Whether letter `position` of the unit's frame lies past the NULLB after frame -3, the last
    /// byte that NCBI's translation buffer sets (`frame_offsets[6]`).
    fn past_translation(&self, position: usize) -> bool {
        self.frame_offsets[self.frame_index] + 1 + position > self.frame_offsets[6]
    }

    /// The bytes NCBI reads in the subject's window of a gapped start at `s_start` of `chunk` of
    /// the unit's `frame`, with every letter past the NULLB after frame -3 read as the sentinel
    /// NULLB (`past_translation`; PD-LOSAT-NCBI-DEFECTS exception 3).
    fn subject_window(
        &self,
        frame: &QueryFrame,
        chunk: TranslatedChunk,
        subject_frame: &[u8],
        s_start: usize,
    ) -> Result<[u8; HSP_MAX_WINDOW]> {
        let s_position = chunk.offset + s_start;
        let mut window = [0u8; HSP_MAX_WINDOW];
        for (index, byte) in window.iter_mut().enumerate() {
            if !self.past_translation(s_position + index) {
                *byte = self.frame_byte(frame, chunk, subject_frame, s_position + index)?;
            }
        }
        Ok(window)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:603-605,635-657,703-712 and
// :102-109 (s_IntervalNodeInit)
/// `BlastIntervalTreeAddHSP` of an HSP of the one-hit word finder's chunk (`-window_size 0`). The
/// gapped start of a one-hit HSP of negative length can lie up to 3 letters past its chunk
/// (`x_score_only_past_chunk`), and the HSP's subject range then starts past the tree's `s_max`
/// (the chunk's length + 1, blast_gapalign.c:3829-3830). Two such HSPs in one subject subtree
/// always go right and never part: NCBI splits the node `[s_max, s_max]` into `[s_max + 1, s_max]`,
/// whose middle is `s_max` again, and goes on splitting it while its node pool doubles, so the
/// search never returns (PD-LOSAT-NCBI-DEFECTS, rule 3). A gapped start at the chunk's last letter
/// or past it can also give an alignment of score 0 whose right extension did not run (its start
/// one past its end, kept by a cutoff of 0, blast_gapalign.c:4290-4295, blast_stat.c:5236-5266),
/// which never straddles a node's middle: two of them split `[p, p]` the same way. The search
/// stops where that first node would be made (`BlastIntervalTree::try_add_hsp`). Ranges in order
/// inside the tree's bounds never make such a node, and only this word finder's HSPs are out of
/// order or leave the bounds, so the other searches keep `add_hsp`. NCBI stops there in the search
/// of a query batch, after it has printed the reports of the batches before: the refusal comes
/// after them (`blastinput::app::ErrorAfterBatchReports`).
#[cold]
#[inline(never)]
fn x_add_one_hit_hsp(
    tree: &mut BlastIntervalTree,
    hsp: TreeHsp,
    context_offset: i32,
) -> Result<()> {
    if !tree.try_add_hsp(hsp, context_offset) {
        return Err(anyhow::Error::new(
            crate::blastinput::app::ErrorAfterBatchReports::new(anyhow::anyhow!(
                "TBLASTN -window_size 0: two one-hit alignments of negative length whose gapped starts lie at the end of their subject chunk or past it (alignments past the chunk, or of score 0 and no letters) make NCBI BLAST+ split a node of its interval tree without end (it does not finish the search there), which LOSAT does not support"
            )),
        ));
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3405-3413
/// The gapped start of a one-hit HSP of negative `length` at `q_start` of the context's letters
/// (`query_sequence`) and `s_start` of `chunk`: the query's window with the letters past the
/// context read as NULLB (`gapped_start_window`), the subject's window in NCBI's translation
/// buffer and past it the sentinel (`XOneHitBuffers::subject_window`). Out of line, so that the
/// gapped stage of the other HSPs keeps its code.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_gapped_start_negative<S: super::score::TbnScore>(
    score: S,
    one_hit: Option<&XOneHitBuffers>,
    query_sequence: &[u8],
    frame: &QueryFrame,
    chunk: TranslatedChunk,
    subject_frame: &[u8],
    q_start: usize,
    s_start: usize,
    length: i32,
    matrix: ScoringMatrix,
) -> Result<usize> {
    let Some(buffers) = one_hit else {
        bail!(
            "TBLASTN: an initial HSP of negative length outside the one-hit word finder (internal error)"
        );
    };
    let query_window = gapped_start_window(query_sequence, q_start);
    let subject_window = buffers.subject_window(frame, chunk, subject_frame, s_start)?;
    score.gapped_start_int4(&query_window, &subject_window, q_start, length, matrix)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4251-4258,4266-4281,4295,4342-4345
/// `s_BlastProtGappedAlignment` of a gapped start at or past the end of its chunk
/// (`s_gapped_start >= chunk.length`), which only the window of a one-hit HSP of negative length
/// gives: NCBI's `subject->sequence` goes on past the chunk with the frame's next letters or, past
/// the frame's last one, with its NULLB and the next frame's letters, and the left extension reads
/// back from the start over them; `subject->length` is the chunk's, so the right extension does
/// not run and the HSP ends at the gapped start. It then goes through the purge,
/// `Blast_HSPListAdjustOffsets` and the merge of the chunks as any other HSP.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_score_only_past_chunk<S: super::score::TbnScore>(
    score: S,
    query_sequence: &[u8],
    frame: &QueryFrame,
    chunk: TranslatedChunk,
    subject_frame: &[u8],
    one_hit: Option<&XOneHitBuffers>,
    q_gapped_start: usize,
    s_gapped_start: usize,
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    scratch: &mut S::Scratch,
) -> Result<crate::algorithm::blastp::gapalign::BlastpScoreOnlyGappedAlignment> {
    let Some(buffers) = one_hit else {
        bail!(
            "TBLASTN: a gapped start past its chunk outside the one-hit word finder (internal error)"
        );
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1094-1101
    // A start at or past the frame's end (a chunk near it) needs a window that sums above 0
    // over the frame's NULLB, which only a PSSM whose '-' column is not BLAST_SCORE_MIN gives
    // (`TbnScore::GAPPED_START_PAST_FRAME`). NCBI's `subject->sequence` is the frame inside the
    // six frames' buffer: the left extension reads back from the start over the next frame's
    // letters (a start past the NULLB), the frame's NULLB (scored with the PSSM's '-' column) and
    // the frame's letters, all set by `BLAST_GetAllTranslations`. The HSP may then be dropped
    // before its traceback, or reach it (`x_traceback_start_past_frame`). A window that passes the
    // NULLB after frame -3 reads the letters past NCBI's buffer as NULLBs (exception 3); a start
    // that then lies past that NULLB is refused (`x_gapped_start_past_translation`).
    if buffers.past_translation(chunk.offset + s_gapped_start) {
        return Err(x_gapped_start_past_translation());
    }
    let in_frame = frame.aa_len.saturating_sub(chunk.offset);
    let bytes = if s_gapped_start < in_frame {
        // The frame's letters from the chunk's first one to the gapped start.
        subject_frame.get(..=s_gapped_start)
    } else {
        // The frame's letters from the chunk's first one, its NULLB, and the bytes of the buffer
        // after it (empty where `subject_frame` does not hold the frame's letters: an internal
        // error below).
        let held = buffers.frame_end.get_or_init(|| {
            let Some(letters) = subject_frame.get(..in_frame) else {
                return Vec::new();
            };
            let mut bytes = Vec::with_capacity(in_frame + 1 + buffers.after_frame.len());
            bytes.extend_from_slice(letters);
            bytes.push(0);
            bytes.extend_from_slice(&buffers.after_frame);
            bytes
        });
        held.get(..=s_gapped_start)
    }
    .with_context(|| {
        format!(
            "TBLASTN -window_size 0: letter {} past a chunk is not held (internal error)",
            s_gapped_start
        )
    })?;
    score.score_only_past_chunk(
        query_sequence,
        bytes,
        chunk.length,
        q_gapped_start,
        s_gapped_start,
        matrix,
        gap_open,
        gap_extend,
        x_drop,
        scratch,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3393-3404,3942-3945,
// 4300-4333
// An initial alignment starting before its chunk crashes NCBI BLAST+ when:
// - its length read as Uint4 exceeds HSP_MAX_WINDOW (11), including a negative one-hit length
//   (`BlastGetStartForGappedAlignment` reads the subject at the start taken as Uint4, about 4 GiB
//   past the frame); or
// - its length is at most 11 (so the gapped start is the alignment's middle without reading a
//   letter), and the subject gapped start still lies before the chunk (`Blast_SemiGappedAlign`
//   reads the right extension from the subject at that start taken as Uint4).
/// Whether a gapped search of `hit` crashes NCBI BLAST+ because `hit` starts before its chunk and
/// either exceeds `HSP_MAX_WINDOW` (11, as `Uint4`, so negative lengths too) or has its middle
/// before the chunk (PD-LOSAT-NCBI-DEFECTS, rule 3; series NP, M4b-T3).
pub(super) fn x_crashes_ncbi(hit: &InitHsp) -> bool {
    hit.s_start < 0
        && (hit.length < 0
            || hit.length > HSP_MAX_WINDOW as i32
            || hit.s_start + hit.length / 2 < 0)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3401-3404
/// The gapped query start (middle of the initial alignment, read from no letter, matrix and PSSM
/// alike) of a non-crashing initial alignment that starts before its chunk. Out of line so the
/// common gapped path keeps its code.
#[cold]
#[inline(never)]
fn x_gapped_start_before_chunk_middle(q_start: usize, hit: &InitHsp) -> Result<usize> {
    if hit.length < 0 || hit.length > HSP_MAX_WINDOW as i32 {
        bail!(
            "LOSAT internal error: negative-start initial alignment of length {} reached gapped scoring",
            hit.length
        );
    }
    Ok(q_start + (hit.length as usize) / 2)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3393-3413,3936-3953,
// 4300-4333 and c++/src/algo/blast/core/blast_itree.c:822-831
/// The refusal (PD-LOSAT-NCBI-DEFECTS, rule 3; series NP, M4b-T3) of a gapped search of a unit
/// that saved an initial HSP meeting [`x_crashes_ncbi`]: starting before its chunk (`-db_soft_mask`
/// with masks stored out of order on a frame longer than 5,000,000 residues) and either longer
/// than `HSP_MAX_WINDOW` (11, as `Uint4`) or with its middle before the chunk. Every gapped HSP
/// saved in a chunk starts in the chunk, so none contains it; NCBI BLAST+ reads past its buffer
/// (`BlastGetStartForGappedAlignment` or `Blast_SemiGappedAlign`) and ends with SIGSEGV. LOSAT
/// stops the search there and names the sequence and frame.
#[cold]
#[inline(never)]
fn x_gapped_start_before_chunk(
    subject: super::x_subject_store::XPrelimSubject<'_>,
    frame: i8,
    chunk_offset: usize,
) -> anyhow::Error {
    let name = match subject {
        super::x_subject_store::XPrelimSubject::Stream(stream, ordinal, _) => {
            format!("database sequence {}", stream.seqdb().oid(ordinal as u32))
        }
        super::x_subject_store::XPrelimSubject::Nt(_) => "the subject".to_owned(),
    };
    anyhow::anyhow!(
        "TBLASTN: {name}, frame {frame}: an initial alignment of the chunk from residue {} starts before the chunk (-db_soft_mask with masks stored out of order on a frame longer than 5,000,000 residues) and is longer than 11 residues or has its middle before the chunk, where NCBI BLAST+ reads past its buffer (a segmentation fault), which LOSAT does not support",
        chunk_offset + 1
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3942-3945 and
// c++/src/algo/blast/core/blast_aascan.c:110-113
/// The subject's gapped start of an initial HSP: its seed (`s_seed`, the chunk's `Uint4` offset,
/// which holds a seed before the chunk wrapped below 0) moved by the gapped start's distance from
/// the query's seed, in NCBI's wrapping 32-bit sum, read as `Int4`. Both a seed at -1 whose HSP
/// starts in the chunk and a short HSP starting before the chunk whose middle lies in the chunk
/// (`x_gapped_start_before_chunk_middle`) give a start in the chunk (`>= 0`); a start below 0 is
/// an internal error (such alignments are refused before by [`x_crashes_ncbi`] and
/// [`x_gapped_start_before_chunk`]).
pub(super) fn x_gapped_subject_start(
    s_seed: u32,
    q_gapped_start: i32,
    q_seed: i32,
) -> Result<usize> {
    let start = s_seed.wrapping_add(q_gapped_start.wrapping_sub(q_seed) as u32) as i32;
    usize::try_from(start).context("negative gapped subject start")
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1058-1101 and
// c++/src/algo/blast/core/blast_gapalign.c:4251-4258,4266-4281
/// The refusal (PD-LOSAT-NCBI-DEFECTS, rule 3) of a gapped start past the NULLB after frame -3,
/// the last byte that NCBI's translation buffer of the subject sets: a one-hit HSP of negative
/// length whose window of HSP_MAX_WINDOW letters reaches past that NULLB (read as NULLBs, exception
/// 3) and sums above 0, which only a PSSM whose '-' column scores above BLAST_SCORE_MIN gives
/// (`TbnScore::GAPPED_START_PAST_FRAME`). NCBI's preliminary gapped alignment then reads the
/// subject back from that start: the byte after the NULLB, which the buffer allocates and never
/// sets, or memory past the buffer. Exception 3 covers the window's letters, not this alignment's.
/// NCBI reads there in the search, after its prolog: the refusal comes after it
/// (`blastinput::app::ErrorAfterBatchReports`, as every refusal of a step where NCBI's search
/// crashes).
#[cold]
#[inline(never)]
fn x_gapped_start_past_translation() -> anyhow::Error {
    anyhow::Error::new(crate::blastinput::app::ErrorAfterBatchReports::new(anyhow::anyhow!(
        "TBLASTN -in_pssm with -window_size 0: the gapped start of a one-hit alignment of negative length lies past the end of the subject's translation (NCBI BLAST+ reads past its buffer there), which LOSAT does not support"
    )))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4609 and :569-575,
// blast_hits.c:1175-1177,1189,1208,1225-1227 and blast_util.c:453
/// The refusal of the traceback of an HSP whose gapped start lies after its frame's NULLB (a one-hit
/// HSP of negative length, `-in_pssm` with `-window_size 0`, `x_score_only_past_chunk`): NCBI's
/// traceback reads the subject at the start first, a byte past the last one either translation of
/// the frame sets (past the fence of the partial translation; on the retry after that fence, the
/// translation of the whole frame (`subject.offset` -1: `1 + n/3 + 2` bytes, no fence) past its
/// NULLB, or past its allocation): PSI-TBLASTN's result there is not defined
/// (PD-LOSAT-NCBI-DEFECTS, rule 3). An HSP dropped before its traceback (below the cutoff,
/// purged, contained, or cut by the hit list) leaves the search to go on, as does a start on the
/// NULLB (defined: the fence, then the whole frame's NULLB, scored with the PSSM's '-' column).
/// NCBI reads there in the search, after its prolog: the refusal comes after it
/// (`blastinput::app::ErrorAfterBatchReports`).
#[cold]
#[inline(never)]
fn x_traceback_start_past_frame() -> anyhow::Error {
    anyhow::Error::new(crate::blastinput::app::ErrorAfterBatchReports::new(anyhow::anyhow!(
        "TBLASTN -in_pssm with -window_size 0: the traceback of a one-hit alignment of negative length starts after the NULLB after its frame (NCBI BLAST+ reads a byte its translation does not set there), which LOSAT does not support"
    )))
}

/// The refusal of `x_traceback_start_past_frame` for a gapped start at letter `s_start` of the window
/// whose first byte is subject position `subject_base`, in the frame `frame` of a subject of
/// `subject_length` bases, after its NULLB (`aa_len`). A start on the NULLB goes on: NCBI reads the
/// fence there and, on the whole frame, the NULLB (`stage_e_report::x_refuse_a_column_of_the_nullb`
/// refuses its alignment where NCBI reports it).
fn x_traceback_start_in_frame(
    subject_length: usize,
    frame: i8,
    s_start: usize,
    subject_base: isize,
) -> Result<()> {
    let aa_len = frame_length(subject_length, frame).context("missing translated subject frame")?;
    if s_start as isize + subject_base > aa_len as isize {
        return Err(x_traceback_start_past_frame());
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
// WordFinder -> GetGappedScore -> endpoint purge -> offset adjustment ->
// Blast_HSPListsMerge runs within each chunk; Blast_HSPListAppend ends a frame.
// The preliminary search keeps that function and input order for every unit.
#[derive(Clone, Copy)]
pub(super) struct PreliminaryProfile<'a> {
    pub(super) seg: Option<&'a SegParams>,
    pub(super) soft_masking: bool,
    pub(super) threshold: i32,
    pub(super) window: i32,
    pub(super) word_xdrop: &'a [i32],
    pub(super) word_cutoff: &'a [i32],
    pub(super) mask_lowercase: bool,
    pub(super) matrix: ScoringMatrix,
    pub(super) word_size: usize,
    pub(super) gap_open: i32,
    pub(super) gap_extend: i32,
    pub(super) gap_xdrop: i32,
    pub(super) gapped_cutoff: &'a [i32],
    pub(super) hsp_num_max: usize,
    /// `score_options->gapped_calculation`: false for `-ungapped`, whose units save the initial
    /// HSPs as they are (`BLAST_GetUngappedHSPList`, `ungapped_chunk_hsps`) and merge the chunks
    /// of a frame on their diagonals only.
    pub(super) gapped: bool,
    /// `-subject_besthit`: each frame's list is filtered after each merge of a chunk.
    pub(super) subject_besthit: bool,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-586,804-850:
// one chunk's WordFinder and gapped result precedes the next chunk's scan;
// each frame's append precedes the next frame's scan.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-584,844
// Dispatch point: this is the preliminary stage of s_BlastSearchEngineCore for one subject. The
// (frame, chunk) units run independently and are merged in NCBI order
// (x_preliminary_protein_hsps_by_unit). The unit tests that check the stage against NCBI traces run
// this same path (`preliminary_hsps_of_letters`, `preliminary_hsps_traced`); the test probe next to
// `PreliminaryEvent` records its events.
/// The preliminary stage of one subject, as independent (frame, chunk) units, with the scores of
/// `score` (the matrix, or a PSSM query's rows).
pub(super) fn preliminary_protein_hsps_for_search<S: super::score::TbnScore>(
    queries: &[&[u8]],
    // EXPERIMENT (LOSAT_X_TBNSTORE): the subject's letters, or its packed store (the unit path
    // translates the store).
    subject: super::x_subject_store::XPrelimSubject<'_>,
    db_gencode: u8,
    profile: PreliminaryProfile<'_>,
    pool: Option<&crate::utils::threading::SearchPool<'_>>,
    x_query_side: &XQuerySide,
    x_subject: Option<(&super::search_seed::XSubjectSide, usize)>,
    score: S,
) -> Result<Vec<(usize, GappedHsp)>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:137 and
    // aa_ungapped.c:214-224 (cited at the units' word finder)
    // Dispatch point: `-window_size 0` (no `multiple_hits`) takes the one-hit word finder, in an
    // instance of its own (`ONE_HIT`), so that the two-hit units keep their code.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:503,560-563
    // NCBI extends the initial HSPs with gaps only in a gapped search; an ungapped search keeps
    // them as the unit's HSPs (`ungapped_chunk_hsps`).
    // Dispatch point: `-ungapped` takes the units that save the initial HSPs, in instances of
    // their own (`UNGAPPED`), so that the gapped units keep their code.
    macro_rules! units {
        ($one_hit:literal, $ungapped:literal) => {
            x_preliminary_protein_hsps_by_unit::<S, $one_hit, $ungapped>(
                queries,
                subject,
                db_gencode,
                profile,
                pool,
                x_query_side,
                x_subject,
                score,
            )
        };
    }
    match (profile.window == 0, profile.gapped) {
        (true, true) => units!(true, false),
        (false, true) => units!(false, false),
        (true, false) => units!(true, true),
        (false, false) => units!(false, true),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1241-1256
// XQueryFrames is the query-side state the unit-wise stage derives from the query set alone: the
// query frames with the lowercase and SEG masks merged (as in x_query_frames below), and the
// offsets of the contexts.
/// What the unit-wise preliminary stage derives from the query set alone.
struct XQueryFrames {
    query_frames: Vec<crate::algorithm::tblastx::translation::QueryFrame>,
    context_offsets: Vec<i32>,
    total_query_length: i32,
    query_sequence: Vec<u8>,
    init_context_offsets: Vec<i32>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1241-1256
// The query frames of the preliminary stage, built once so that one value can serve all subjects
// of a query set.
fn x_query_frames(queries: &[&[u8]], profile: PreliminaryProfile<'_>) -> Result<XQueryFrames> {
    let query_frames: Vec<_> = queries
        .iter()
        .map(|query| {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:620-623:
            // Soft masking retains the unmasked WordFinder query while lookup still excludes the
            // SEG intervals used by the candidate scanner.
            if profile.soft_masking {
                encode_protein_query_frame_with_seg(query, None)
            } else {
                encode_tblastn_lookup_query(query, profile.seg, profile.mask_lowercase)
            }
        })
        .collect();
    x_query_frames_of(query_frames, queries)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:614-623 (`BLAST_MainSetUp`)
// NCBI masks the query of a batch once, in its setup; the lookup table, the word finder and the
// score block read that query. `x_query_frames` with the frames made already (by the batch's setup,
// `stage_d_pipeline::query_set_setup`, which masks each query as `x_query_frames` does).
/// The query-side state of `queries` from their frames `query_frames` (one per query).
fn x_query_frames_of(query_frames: Vec<QueryFrame>, queries: &[&[u8]]) -> Result<XQueryFrames> {
    let mut context_offsets = Vec::with_capacity(queries.len());
    let mut total_query_length = 0i32;
    for query in queries {
        context_offsets.push(total_query_length);
        total_query_length += i32::try_from(query.len())? + 1;
    }
    total_query_length -= 1;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:246-250:
    // last context's query_offset + query_length is the sequence length;
    // each preceding protein context ends with one NULLB separator.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:68-96:
    // The concatenated query of the word finder: each masked context followed by its sentinel.
    let mut query_sequence = Vec::new();
    let mut init_context_offsets = Vec::with_capacity(queries.len());
    for frame in &query_frames {
        init_context_offsets.push(i32::try_from(query_sequence.len())?);
        query_sequence.extend_from_slice(&frame.aa_seq[1..]);
    }
    Ok(XQueryFrames {
        query_frames,
        context_offsets,
        total_query_length,
        query_sequence,
        init_context_offsets,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:91-100
// NCBI builds the lookup table once per query batch. XQuerySide keeps the query frames and the
// lookup table for all subjects of one query set. The frames are the ones the batch's setup masked
// (`x_with_frames`); the table is built before the subjects on a pool of two or more threads
// (`x_build_lookup_on`), otherwise by the first subject that needs it.
/// The query-side state of the unit-wise preliminary stage for one query set, shared by all
/// its subjects (graduated 2026-10-22, was `LOSAT_X_TBNQSIDE`).
///
/// The masked query frames and the neighbourhood lookup table depend on the
/// queries and the word-finder options, not on the subject, yet the stage is
/// entered once per subject and rebuilt them every time (NCBI builds its
/// lookup table once per query batch, `lookup_wrap.c:91-100`). The frames come
/// from the batch's setup (`x_with_frames`); the table is built once, before the
/// subjects (`x_build_lookup_on`) or by the first subject that needs it, and the
/// others read it. One value serves one query set with one set of options only.
pub(super) struct XQuerySide {
    frames: std::sync::OnceLock<XQueryFrames>,
    /// Held while `frames` is being made, so that subjects searched at the
    /// same time wait for the first one instead of each making the frames.
    frames_in_progress: std::sync::Mutex<()>,
    lookup:
        std::sync::OnceLock<std::sync::Arc<crate::algorithm::tblastx::lookup::BlastAaLookupTable>>,
    /// The compressed table of word sizes 5 to 7 (tblastn-fast), in place of `lookup`.
    compressed: super::search_seed::XCompressedSide,
    /// EXPERIMENT (LOSAT_X_TBNUNITXLATE): the buffers of the units of the batch's plans without
    /// residues.
    unit_residues: XUnitResiduesPool,
}

// No NCBI counterpart: the methods below create the shared query-side state and build its frames
// once; it does not change any value NCBI computes.
impl XQuerySide {
    pub(super) fn new() -> Self {
        use std::sync::OnceLock;
        Self {
            frames: OnceLock::new(),
            frames_in_progress: std::sync::Mutex::new(()),
            lookup: OnceLock::new(),
            compressed: Default::default(),
            unit_residues: XUnitResiduesPool::new(),
        }
    }

    /// The query frames, made by the first caller. A failure is returned to
    /// that caller and leaves nothing behind, as when every caller made them.
    fn frames(&self, build: impl FnOnce() -> Result<XQueryFrames>) -> Result<&XQueryFrames> {
        if let Some(frames) = self.frames.get() {
            return Ok(frames);
        }
        let _in_progress = self
            .frames_in_progress
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(frames) = self.frames.get() {
            return Ok(frames);
        }
        let built = build()?;
        Ok(self.frames.get_or_init(|| built))
    }

    /// The side of a query set whose frames its setup made (`query_set_setup`, one per query, the
    /// frames `x_query_frames` makes): the query is masked once per batch, as NCBI masks it once
    /// in `BLAST_MainSetUp` (`x_query_frames_of`).
    pub(super) fn x_with_frames(
        setup_frames: &[Vec<QueryFrame>],
        queries: &[&[u8]],
    ) -> Result<Self> {
        anyhow::ensure!(
            setup_frames.len() == queries.len() && setup_frames.iter().all(|f| f.len() == 1),
            "one query frame per TBLASTN query"
        );
        let side = Self::new();
        let frames = setup_frames.iter().map(|f| f[0].clone()).collect();
        let _ = side.frames.set(x_query_frames_of(frames, queries)?);
        Ok(side)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:86-99
    /// The query set's lookup table made by its setup: a PSSM query's (`in_pssm::lookup_table`),
    /// which the subjects' plans then read (`x_scan_plan`) instead of the matrix's.
    pub(super) fn x_set_lookup(
        &self,
        table: std::sync::Arc<crate::algorithm::tblastx::lookup::BlastAaLookupTable>,
    ) {
        let _ = self.lookup.set(table);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:91-100 (cited at
    // `x_scan_plan`) and c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // NCBI builds the batch's table in its setup, before its threads search the subjects. The table
    // `x_scan_plan` builds at the first subject that needs it, built now, its lists on `pool`'s
    // threads (`x_build_lean_in`): the same table, from the frames the setup masked, which are the
    // table's query without soft masking. With soft masking, or the compressed table of word sizes
    // 5 to 7, the first subject builds it.
    /// Builds the query set's lookup table now, on `pool`, where the side holds the table's query.
    pub(super) fn x_build_lookup_on(
        &self,
        profile: PreliminaryProfile<'_>,
        pool: &crate::utils::threading::SearchPool<'_>,
    ) {
        if profile.soft_masking || profile.word_size > 4 {
            return;
        }
        let Some(frames) = self.frames.get() else {
            return;
        };
        self.lookup.get_or_init(|| {
            let query_frames: Vec<Vec<QueryFrame>> = frames
                .query_frames
                .iter()
                .map(|frame| vec![frame.clone()])
                .collect();
            std::sync::Arc::new(
                crate::algorithm::tblastx::lookup::x_lookup_table_for_profile_in(
                    &query_frames,
                    profile.threshold,
                    profile.matrix,
                    profile.word_size,
                    Some(pool),
                ),
            )
        });
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1002-1003,1377-1379,1412-1413
// NCBI's search thread sets up its word-finder state (diagonal table, offset pairs, initial hit
// list) and its gapped alignment memory once, before its loop over the subjects, and keeps them for
// every subject. These are the same buffers of the unit path, kept by a thread between subjects.
/// The scan state of the unit path that a thread keeps from one subject to the next
/// (`X_UNIT_STATE`; graduated 2026-10-22, was `LOSAT_X_TBNSTATE`). A fresh state per subject would
/// be, on a search of many short subjects, hundreds of thousands of tables of up to a megabyte,
/// which glibc's per-thread heaps create and delete again.
struct XUnitBuffers {
    diagonals: super::search_init::Diagonals,
    pairs: Vec<crate::algorithm::tblastx::blast_aascan::BlastOffsetPair>,
    hits: Vec<InitHsp>,
    scratch: GapAlignScratch,
}

thread_local! {
    /// This thread's `XUnitBuffers`, taken while a subject's units
    /// run on the thread and put back after them. A thread that finds it empty (its first subject,
    /// or a subject started while another one on the same thread holds them) makes new buffers.
    static X_UNIT_STATE: std::cell::Cell<Option<XUnitBuffers>> = const { std::cell::Cell::new(None) };
}

/// Frees this thread's `XUnitBuffers` at the end of its loop over the chunks of the subject flow
/// (`x_oid_flow::x_loop_free`, LOSAT_X_TBNLOOPFREE). A thread that searches a subject after it
/// makes new buffers, as for its first subject.
#[cfg(feature = "parallel")]
pub(super) fn x_free_unit_buffers() {
    let _ = X_UNIT_STATE.try_with(|state| drop(state.take()));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-812
// NCBI's units read the subject's translation buffer; a unit of a plan without residues translates
// the part it reads into a buffer of its query batch's preliminary stage (up to a chunk of
// 5,000,000 residues and the margin, about 5 MB), freed at the end of that stage.
/// EXPERIMENT (LOSAT_X_TBNUNITXLATE): the idle buffers for the residues of the units of plans
/// without residues (`x_packed_lean_plan`) of one query batch (`XQuerySide`), shared by its
/// threads: a unit takes one (or a new empty one) and gives it back when it ends, and they are
/// freed with the batch's query-side state at the end of its preliminary stage. Before, each
/// thread kept its buffer until it ended, so the buffers of every thread that had run such a unit
/// were still there at the peak, which comes after the units (16 threads: 16 x 4.77 MiB).
pub(super) struct XUnitResiduesPool {
    idle: std::sync::Mutex<Vec<Vec<u8>>>,
}

// No NCBI counterpart: the pool above; it computes no value NCBI computes.
impl XUnitResiduesPool {
    fn new() -> Self {
        Self {
            idle: std::sync::Mutex::new(Vec::new()),
        }
    }

    fn idle(&self) -> std::sync::MutexGuard<'_, Vec<Vec<u8>>> {
        self.idle
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    /// A buffer for one unit: an idle one, or a new empty one.
    fn take(&self) -> XUnitResidues<'_> {
        let residues = self.idle().pop();
        if residues.is_none() {
            super::x_subject_translate::x_count_unit_buffers(1, 0);
        }
        XUnitResidues(self, residues.unwrap_or_default())
    }
}

impl Drop for XUnitResiduesPool {
    fn drop(&mut self) {
        let idle = self
            .idle
            .get_mut()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if !idle.is_empty() {
            super::x_subject_translate::x_count_unit_buffers(0, idle.len() as u64);
        }
    }
}

// No NCBI counterpart: the hold of a unit on a buffer of the pool (above); it computes no value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNUNITXLATE): a unit's buffer, given back to the pool when the unit ends,
/// also when it fails. A buffer only grows (`x_unit_residues`), so it is filled with zeros once.
struct XUnitResidues<'a>(&'a XUnitResiduesPool, Vec<u8>);

impl Drop for XUnitResidues<'_> {
    fn drop(&mut self) {
        let residues = std::mem::take(&mut self.1);
        if residues.capacity() > 0 {
            self.0.idle().push(residues);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-584,844
// The units of the loop at line 478 run independently here: WordFinder (two-hit pass with its own
// diagonal array), the gapped score and the common-endpoint purge are computed per (frame, chunk)
// unit; the merge (581) and the per-frame append (844) then run in NCBI order on one thread. This
// reproduces the order-dependent steps of the serial loop; the units read nothing an earlier unit
// wrote.
// The preliminary stage as independent (frame, chunk) units (graduated 2026-10-22, was
// `LOSAT_X_TBNPAR`).  NCBI c++/src/algo/blast/core/blast_engine.c:478-552 runs
// WordFinder, GetGappedScore and the endpoint purge of one unit before it scans
// the next; none of the three reads anything an earlier unit produced (the
// diagonal array is the only carried state, see `XTwoHit`).  Only the merge
// (lines 572-586) and the per-frame append (840-850) depend on the order, and
// they run here afterwards, in that order, on one thread.  The offset pairs are
// consumed straight from the scan buffer, so no per-frame seed list is kept.
fn x_preliminary_protein_hsps_by_unit<
    S: super::score::TbnScore,
    const ONE_HIT: bool,
    const UNGAPPED: bool,
>(
    queries: &[&[u8]],
    subject: super::x_subject_store::XPrelimSubject<'_>,
    db_gencode: u8,
    profile: PreliminaryProfile<'_>,
    pool: Option<&crate::utils::threading::SearchPool<'_>>,
    x_query_side: &XQuerySide,
    x_subject: Option<(&super::search_seed::XSubjectSide, usize)>,
    score: S,
) -> Result<Vec<(usize, GappedHsp)>> {
    use super::search_init::XTwoHit;
    use super::search_seed::{x_scan_plan, x_scan_unit, x_scan_unit_of};
    use super::x_subject_store::XPrelimSubject;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3924-3927:
    // Every query context must supply its own measured cutoff.
    if profile.gapped_cutoff.len() != queries.len() {
        bail!("one preliminary gapped cutoff is required per query context");
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:547-583:
    // Require one x-drop and cutoff per query context for the indexed read.
    anyhow::ensure!(
        profile.word_xdrop.len() == queries.len() && profile.word_cutoff.len() == queries.len(),
        "one WordFinder x-drop and cutoff per query context"
    );
    let x_frames: &XQueryFrames = x_query_side.frames(|| x_query_frames(queries, profile))?;
    anyhow::ensure!(
        x_frames.query_frames.len() == queries.len(),
        "query-side state used for another query set"
    );
    let query_frames = &x_frames.query_frames;
    let context_offsets = &x_frames.context_offsets;
    let total_query_length = x_frames.total_query_length;
    let query_sequence = &x_frames.query_sequence;
    let init_context_offsets = &x_frames.init_context_offsets;
    // A database sequence whose plan has no residues
    // (LOSAT_X_TBNUNITXLATE, `x_lean_subject`) is read before its plan (`read_shared`) and
    // searched as the packed subject of the copy is: its units translate the store they read.
    let x_store;
    let subject = match subject {
        // A masked search (`-db_soft_mask`, `-db_hard_mask`) makes every plan from its ranges
        // (`x_scan_plan`), with residues.
        XPrelimSubject::Stream(stream, oid, length)
            if stream.mask().is_none()
                && super::search_seed::x_lean_subject(
                    x_subject,
                    length,
                    db_gencode,
                    profile.mask_lowercase,
                ) =>
        {
            x_store = stream.read_shared(oid)?;
            XPrelimSubject::Nt(SubjectNt::Packed(&x_store))
        }
        subject => subject,
    };
    let Some(plan) = x_scan_plan(
        queries,
        subject,
        db_gencode,
        profile.seg,
        profile.threshold,
        profile.mask_lowercase,
        profile.matrix,
        profile.word_size,
        &x_query_side.lookup,
        &x_query_side.compressed,
        x_subject,
        // Without soft masking the frames are the lookup query's (both
        // `encode_tblastn_lookup_query` with the profile's SEG and lowercase options), masked once
        // for the batch.
        (!profile.soft_masking).then_some(query_frames.as_slice()),
    )?
    else {
        return Ok(Vec::new());
    };
    struct UnitState<'a, M> {
        two_hit: XTwoHit<'a>,
        pairs: Vec<crate::algorithm::tblastx::blast_aascan::BlastOffsetPair>,
        hits: Vec<InitHsp>,
        scratch: M,
    }
    // The units end with NCBI's offset step instead of starting on a cleared table, and the
    // sequential path below runs on this thread's kept buffers. `x_in_order`: a new table of
    // `(last_hit, flag)` cells (below).
    let state_on = |kept: Option<XUnitBuffers>, x_in_order: bool| {
        let (table, pairs, hits, scratch) = match kept {
            Some(kept) => (Some(kept.diagonals), kept.pairs, kept.hits, kept.scratch),
            None => (None, Vec::new(), Vec::new(), GapAlignScratch::new()),
        };
        let two_hit = if x_in_order {
            XTwoHit::x_with_layout(
                false,
                None,
                query_sequence,
                init_context_offsets,
                profile.word_xdrop,
                profile.word_cutoff,
                profile.matrix,
                profile.window,
                profile.word_size,
            )
        } else {
            XTwoHit::with_table(
                table,
                query_sequence,
                init_context_offsets,
                profile.word_xdrop,
                profile.word_cutoff,
                profile.matrix,
                profile.window,
                profile.word_size,
            )
        };
        UnitState {
            two_hit: S::x_two_hit(two_hit),
            pairs,
            hits,
            scratch: S::unit_scratch(scratch),
        }
    };
    #[cfg(feature = "parallel")]
    let new_state = || state_on(None, false);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:300-306,478-498,1377 and
    // c++/src/algo/blast/core/aa_ungapped.c:516-547 and
    // c++/src/algo/blast/core/blast_extend.c:87-107,162-173
    // A unit whose scan reads letters before its chunk (`-db_soft_mask` masks stored out of order,
    // `scan_lead`) or a `-db_hard_mask` subject recorded by `translated_walk` (`hard_start`:
    // overlapping masks or a sequence that meets unwritten entries / clears inside it) runs all
    // six frames' units in order on a fresh `(last_hit, flag)` diagonal table initialized to its
    // `TableStart::New` offset.
    let hard_start = match subject {
        super::x_subject_store::XPrelimSubject::Stream(stream, ordinal, _) => {
            stream.hard_start(ordinal)
        }
        super::x_subject_store::XPrelimSubject::Nt(_) => None,
    };
    let x_in_order = hard_start.is_some()
        || plan.subject.units.iter().any(|(_, _, ranges)| {
            crate::algorithm::common::translated_chunks::scan_lead(ranges) > 0
        });
    let run_unit = |state: &mut UnitState<'_, S::Scratch>,
                    unit: usize|
     -> Result<Vec<(usize, GappedHsp)>> {
        let (frame_index, chunk, ranges) = &plan.subject.units[unit];
        let frame = &plan.subject.frames[*frame_index];
        let chunk = *chunk;
        // The letters before the chunk that the unit's scan reads (masks stored out of order).
        let x_lead = crate::algorithm::common::translated_chunks::scan_lead(ranges);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-812
        // EXPERIMENT (LOSAT_X_TBNUNITXLATE): the frame from `x_lead` residues before the chunk's
        // first one on: the plan's, or for a plan without residues those the unit reads,
        // translated into a buffer of the query batch (`XUnitResiduesPool`).
        let mut x_held = None;
        let x_residues: &[u8] = if frame.aa_seq.is_empty() {
            let XPrelimSubject::Nt(SubjectNt::Packed(store)) = subject else {
                bail!("TBLASTN: a plan without residues of a subject that is not packed");
            };
            let held = x_held.insert(x_query_side.unit_residues.take());
            let range = super::x_subject_translate::x_unit_residues(
                store,
                db_gencode,
                frame,
                TranslatedChunk {
                    offset: chunk.offset - x_lead,
                    length: chunk.length + x_lead,
                    ..chunk
                },
                &mut held.1,
            )?;
            &held.1[range]
        } else {
            &frame.aa_seq[1 + chunk.offset - x_lead..]
        };
        let subject_residues = &x_residues[x_lead..];
        let subject_sequence = &subject_residues[..chunk.length];
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:300-306
        // The letters around the chunk for the extensions of the hits before it (a word past the
        // chunk's last letter is never extended: `XTwoHit::pairs_on`).
        let x_around = (x_lead > 0).then(|| ChunkLetters {
            letters: &x_residues[..x_lead + chunk.length],
            lead: x_lead,
            length: chunk.length,
        });
        // The buffers that the one-hit word finder's HSPs of negative length read (`-window_size
        // 0` only, in the gapped score).
        let x_one_hit = if ONE_HIT && !UNGAPPED {
            Some(XOneHitBuffers::new(
                &plan.subject.frames,
                *frame_index,
                chunk,
                subject,
                db_gencode,
            )?)
        } else {
            None
        };
        let chunk_offset = u32::try_from(chunk.offset)?;
        let UnitState {
            two_hit,
            pairs,
            hits,
            scratch,
        } = state;
        hits.clear();
        if x_around.is_some() {
            two_hit.x_set_around(x_around);
        }
        #[cfg(test)]
        probe_begin_unit(unit, frame.frame, chunk.offset);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:214-224
        // Dispatch point: `-window_size 0` (no `multiple_hits`, blast_extend.c:137) takes the
        // one-hit word finder, with its own scan (no clamp of a short first range) and loop
        // (`XTwoHit::x_pairs_one_hit`); the two-hit loops below keep their code. Then a PSSM
        // query's walks (`score::XPssmWalks`, in a scan of their own so that the matrices keep
        // their code), else the matrix's (BLOSUM62's or another's table).
        // An extension that cannot run (`x_pairs_one_hit`), reported after the table's exit.
        let mut x_one_hit_failed = None;
        if ONE_HIT {
            x_scan_unit_of::<false>(&plan, unit, x_residues, x_lead, pairs, |filled| {
                #[cfg(test)]
                probe_seeds(unit, frame.frame, chunk_offset, filled);
                if let Err(error) = two_hit.x_pairs_one_hit(
                    score,
                    frame.frame,
                    chunk_offset,
                    subject_sequence,
                    filled,
                    hits,
                ) {
                    x_one_hit_failed.get_or_insert(error);
                }
            })
        } else {
            match (score.x_walks(), two_hit.x_ungap_matrix()) {
                (Some(walks), _) => x_scan_unit(&plan, unit, x_residues, x_lead, pairs, |filled| {
                    #[cfg(test)]
                    probe_seeds(unit, frame.frame, chunk_offset, filled);
                    two_hit.x_pairs_bf_walks(
                        walks,
                        frame.frame,
                        chunk_offset,
                        subject_sequence,
                        filled,
                        hits,
                    )
                }),
                (None, None) => x_scan_unit(&plan, unit, x_residues, x_lead, pairs, |filled| {
                    #[cfg(test)]
                    probe_seeds(unit, frame.frame, chunk_offset, filled);
                    // The two-hit loop of the units without a branch on the common outcomes
                    // (tblastn/x_tbn_hitbf.rs).
                    two_hit.x_pairs_bf(frame.frame, chunk_offset, subject_sequence, filled, hits)
                }),
                // Another matrix's byte-table walks (`XTwoHit::x_pairs_bf_of`), in a scan of its own
                // so that BLOSUM62's keeps its code.
                (None, Some(table)) => {
                    x_scan_unit(&plan, unit, x_residues, x_lead, pairs, |filled| {
                        #[cfg(test)]
                        probe_seeds(unit, frame.frame, chunk_offset, filled);
                        two_hit.x_pairs_bf_of(
                            table,
                            frame.frame,
                            chunk_offset,
                            subject_sequence,
                            filled,
                            hits,
                        )
                    })
                }
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:614,807:
        // Right after the scan, before anything that can fail, so that the table is ready for the
        // next unit whatever this one returns.
        two_hit.exit(chunk.ncbi_length);
        if x_around.is_some() {
            two_hit.x_set_around(None);
        }
        if ONE_HIT {
            if let Some(error) = x_one_hit_failed {
                return Err(error);
            }
        }
        XTwoHit::sort(hits);
        #[cfg(test)]
        probe_event(unit, || {
            PreliminaryEvent::Initial(frame.frame, chunk.offset, hits.to_vec())
        });
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:500:
        if hits.is_empty() {
            return Ok(Vec::new());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3393-3404,3942-3945,
        // 4300-4333 and c++/src/algo/blast/core/blast_itree.c:822-831
        // Every gapped alignment saved in a chunk starts in the chunk, so an initial alignment
        // starting before it is never skipped as contained; refuse exactly the alignments on which
        // NCBI BLAST+ crashes (`x_crashes_ncbi`), in unit order, after the scan with NCBI's table.
        if !UNGAPPED && x_lead > 0 && hits.iter().any(x_crashes_ncbi) {
            return Err(x_gapped_start_before_chunk(
                subject,
                frame.frame,
                chunk.offset,
            ));
        }
        if UNGAPPED {
            return Ok(ungapped_chunk_hsps(context_offsets, frame.frame, hits));
        }
        let gapped = score_gapped_chunk::<S, ONE_HIT>(
            queries,
            query_frames,
            context_offsets,
            total_query_length,
            frame,
            chunk,
            subject_sequence,
            subject_residues,
            hits,
            profile.gap_open,
            profile.gap_extend,
            profile.gap_xdrop,
            profile.gapped_cutoff,
            profile.matrix,
            scratch,
            score,
            x_one_hit.as_ref(),
        )?;
        #[cfg(test)]
        let gapped = probe_gapped(unit, frame.frame, chunk.offset, gapped)?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:539-552,555:
        let purged = purge_preliminary_common_endpoints(gapped);
        #[cfg(test)]
        probe_event(unit, || {
            PreliminaryEvent::Purged(frame.frame, chunk.offset, purged.clone())
        });
        Ok(purged)
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-482,805-812
    // Dispatch point: NCBI visits the units one after another. With a pool the units run in
    // parallel (map_init gives each thread its own state); without one the loop below runs them in
    // order. The results are collected by unit index.
    let unit_count = plan.subject.units.len();
    let parallel = pool.is_some_and(crate::utils::threading::SearchPool::enabled)
        && unit_count > 1
        && !x_in_order;
    // No NCBI counterpart: stage report for the thread-use statistics; it does not change any value
    // NCBI computes.
    crate::utils::threading::report_stage("tblastn", "prelim_units", unit_count, parallel);
    #[cfg(test)]
    probe_require_serial(parallel);
    let mut purged: Vec<Result<Vec<(usize, GappedHsp)>>> = Vec::with_capacity(unit_count);
    #[cfg(feature = "parallel")]
    if parallel {
        use rayon::prelude::*;
        let pool = pool.expect("parallel pool selected");
        purged = pool.install(|| {
            (0..unit_count)
                .into_par_iter()
                .with_max_len(1)
                .map_init(&new_state, |state, unit| run_unit(state, unit))
                .collect()
        });
    }
    if purged.is_empty() {
        let mut state = state_on(
            if x_in_order {
                None
            } else {
                X_UNIT_STATE.take()
            },
            x_in_order,
        );
        if let Some(start) = hard_start {
            use crate::algorithm::common::hard_chunk_walk::{ChunkWork, DiagonalWalk, TableStart};
            state.two_hit.apply_table_start(TableStart::New(start))?;
            let mut walk = DiagonalWalk::new(profile.word_size as i32, profile.window, 0, false);
            walk.start_on_ncbi_table(start);
            let mut walk_frame: Option<usize> = None;
            for unit in 0..unit_count {
                let (frame_index, chunk, _) = &plan.subject.units[unit];
                if walk_frame != Some(*frame_index) {
                    walk.frame(plan.subject.frames[*frame_index].aa_len);
                    walk_frame = Some(*frame_index);
                }
                let clears = state.two_hit.offset() >= i32::MAX / 4;
                let res = match walk.chunk(chunk.offset, chunk.ncbi_length) {
                    ChunkWork::Words { .. } => run_unit(&mut state, unit),
                    ChunkWork::OneWord => x_run_one_word_unit::<S, UNGAPPED>(
                        &plan,
                        unit,
                        subject,
                        queries,
                        query_frames,
                        query_sequence,
                        context_offsets,
                        total_query_length,
                        profile,
                        score,
                        &mut state.two_hit,
                        &mut state.pairs,
                        &mut state.hits,
                        &mut state.scratch,
                    ),
                    ChunkWork::Nothing => {
                        state.two_hit.exit(chunk.ncbi_length);
                        Ok(Vec::new())
                    }
                };
                if clears {
                    walk.clear();
                }
                purged.push(res);
            }
        } else {
            if x_in_order {
                if let super::x_subject_store::XPrelimSubject::Stream(stream, ordinal, _) = subject
                {
                    let start = stream.lead_start(ordinal, profile.word_size, profile.window)?;
                    state.two_hit.apply_table_start(
                        crate::algorithm::common::hard_chunk_walk::TableStart::New(start),
                    )?;
                }
            }
            for unit in 0..unit_count {
                purged.push(run_unit(&mut state, unit));
            }
        }
        if !x_in_order {
            let UnitState {
                two_hit,
                pairs,
                hits,
                scratch,
            } = state;
            X_UNIT_STATE.set(Some(XUnitBuffers {
                diagonals: two_hit.into_table(),
                pairs,
                hits,
                scratch: S::kept_scratch(scratch),
            }));
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:581-584,844
    // Replay in NCBI order: the units of a frame are merged one by one (Blast_HSPListsMerge), then
    // the frame is appended (Blast_HSPListAppend). Only this thread runs it.
    let mut purged = purged.into_iter();
    let mut combined = Vec::new();
    let mut unit = 0usize;
    for frame_index in 0..plan.subject.frames.len() {
        let mut per_frame = Vec::new();
        while unit < unit_count && plan.subject.units[unit].0 == frame_index {
            #[cfg(test)]
            probe_replay_unit(unit);
            let chunk = plan.subject.units[unit].1;
            let unit_hsps = purged.next().expect("one result per unit")?;
            unit += 1;
            if unit_hsps.is_empty() {
                continue;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:572-586:
            if hard_start.is_some() && chunk.overlap_past_list && !per_frame.is_empty() {
                x_check_hard_merge_overlap(subject, chunk, &per_frame, &unit_hsps)?;
            }
            per_frame =
                merge_adjusted_chunk(per_frame, unit_hsps, chunk, profile.hsp_num_max, !UNGAPPED)?;
            #[cfg(test)]
            probe_push(|| {
                PreliminaryEvent::Merged(
                    plan.subject.frames[frame_index].frame,
                    chunk.offset,
                    per_frame.clone(),
                )
            });
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:586-591
            // `-subject_besthit` on the frame's list after each merge, before the frame is
            // appended: the HSPs of other subject frames are not in it.
            if profile.subject_besthit {
                crate::core::blast_hits::subject_besthit(&mut per_frame, false, |(query, hsp)| {
                    crate::core::blast_hits::BestHitExtent {
                        context: *query as i64,
                        query_frame: 0,
                        offset: i64::from(hsp.q_start),
                        end: i64::from(hsp.q_end),
                        query_length: 0,
                    }
                });
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:840-850:
        combined = append_preliminary_hsps(combined, per_frame, profile.hsp_num_max);
        #[cfg(test)]
        probe_push(|| {
            PreliminaryEvent::Appended(plan.subject.frames[frame_index].frame, combined.clone())
        });
    }
    Ok(combined)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-505, 516-614
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3827-3953, 4246-4295
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:510-795
/// Cold two-hit scan and gapped/ungapped stage for a `-db_hard_mask` chunk shorter than the word
/// (`ChunkWork::OneWord`): scans the single word at the chunk's offset from the frame's stored
/// letters through `common::hard_chunk_scan::OneWord` on TBLASTN's `(last_hit, flag)` table,
/// advances the diagonal offset by `chunk.ncbi_length`, and runs either the ungapped list or the
/// short-chunk gapped stage (refusing with `short_chunk_tree_refusal` where NCBI's interval tree
/// splits without end).
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_run_one_word_unit<S: super::score::TbnScore, const UNGAPPED: bool>(
    plan: &super::search_seed::XScanPlan,
    unit: usize,
    subject: super::x_subject_store::XPrelimSubject<'_>,
    queries: &[&[u8]],
    query_frames: &[QueryFrame],
    query_sequence: &[u8],
    context_offsets: &[i32],
    total_query_length: i32,
    profile: PreliminaryProfile<'_>,
    score: S,
    two_hit: &mut super::search_init::XTwoHit<'_>,
    pairs: &mut Vec<crate::algorithm::tblastx::blast_aascan::BlastOffsetPair>,
    hits: &mut Vec<InitHsp>,
    scratch: &mut S::Scratch,
) -> Result<Vec<(usize, GappedHsp)>> {
    use crate::algorithm::common::hard_chunk_scan::{
        short_chunk_tree_refusal, OneWord, WordContext,
    };
    let (frame_index, chunk, _) = &plan.subject.units[unit];
    let frame = &plan.subject.frames[*frame_index];
    let chunk = *chunk;
    let chunk_offset = u32::try_from(chunk.offset)?;
    let word_size_i32 = i32::try_from(profile.word_size)?;
    anyhow::ensure!(
        !frame.aa_seq.is_empty(),
        "TBLASTN: one-word hard chunk reached a frame without residues"
    );
    let stored = &frame.aa_seq[1..1 + frame.aa_len];
    let mut word_letters = Vec::with_capacity(profile.word_size + 1);
    if chunk.offset < stored.len() {
        let end = (chunk.offset + profile.word_size).min(stored.len());
        word_letters.extend_from_slice(&stored[chunk.offset..end]);
    }
    word_letters.resize(profile.word_size + 1, 0);

    let lookup_hits = super::search_seed::x_one_word_hits(
        &plan.lookup,
        &word_letters,
        word_size_i32,
        pairs,
        plan.pair_capacity,
    );
    hits.clear();
    let (mut ring, offset) = two_hit
        .pair_ring()
        .context("TBLASTN: one-word hard chunk requires a Pair diagonal table")?;
    let test = OneWord {
        query: query_sequence,
        letters: &word_letters[..profile.word_size],
        length: chunk.ncbi_length,
        word: word_size_i32,
        window: profile.window,
        scores: crate::utils::matrix::protein_score_rows(profile.matrix),
    };
    test.test_on_with(
        &mut ring,
        &lookup_hits,
        offset,
        total_query_length,
        |q, s_letter| score.residue_score(query_sequence, q, s_letter, profile.matrix),
        |q| {
            let ctx = context_offsets.partition_point(|&o| o <= q) - 1;
            WordContext {
                start: context_offsets[ctx],
                x_drop: profile.word_xdrop[ctx],
                cutoff: profile.word_cutoff[ctx],
            }
        },
        |q, seed| {
            hits.push(InitHsp {
                frame: frame.frame,
                chunk_offset,
                q_seed: q as u32,
                s_seed: 0,
                q_start: seed.q_start,
                s_start: seed.s_start,
                length: seed.length,
                score: seed.score,
            });
        },
    );
    two_hit.exit(chunk.ncbi_length);
    super::search_init::XTwoHit::sort(hits);
    if hits.is_empty() {
        return Ok(Vec::new());
    }
    if UNGAPPED {
        return Ok(ungapped_chunk_hsps(context_offsets, frame.frame, hits));
    }
    let s_max = chunk.ncbi_length.saturating_add(1);
    let mut tree = BlastIntervalTree::new(0, total_query_length + 1, 0, s_max);
    tree.reset_with_bounds_cells(
        0,
        total_query_length + 1,
        0,
        s_max,
        x_grid_cells(hits.len()),
    );
    let mut results = Vec::new();
    let subject_frame_sign = i32::from(frame.frame.signum());
    for hit in hits.iter() {
        let context = context_offsets.partition_point(|&o| o <= hit.q_seed as i32) - 1;
        let context_offset = context_offsets[context];
        let query_length = i32::try_from(queries[context].len())?;
        let ctx_query = &query_frames[context].aa_seq[1..1 + queries[context].len()];
        let initial_tree_hsp = TreeHsp {
            query_offset: hit.q_start - context_offset,
            query_end: hit.q_start - context_offset + hit.length,
            subject_offset: hit.s_start,
            subject_end: hit.s_start + hit.length,
            score: hit.score,
            query_frame: 0,
            query_length,
            query_context_offset: context_offset,
            subject_frame_sign,
        };
        if tree.contains_hsp(&initial_tree_hsp, context_offset, 0) {
            continue;
        }
        let q_start = usize::try_from(hit.q_start - context_offset)
            .context("negative context-local initial query start")?;
        let s_start = usize::try_from(hit.s_start).context("negative one-word subject start")?;
        let length = usize::try_from(hit.length).context("negative one-word HSP length")?;
        let q_gapped_start = score.gapped_start(
            ctx_query,
            &word_letters[..profile.word_size],
            q_start,
            length,
            s_start,
            profile.matrix,
        );
        let s_gapped_start = x_gapped_subject_start(
            hit.s_seed,
            i32::try_from(q_gapped_start)?,
            i32::try_from(hit.q_seed)? - context_offset,
        )?;
        let letters = usize::try_from(chunk.ncbi_length)
            .unwrap_or(0)
            .max(s_gapped_start + 1)
            .min(word_letters.len());
        let gapped = score.score_only(
            ctx_query,
            &word_letters[..letters],
            q_gapped_start,
            s_gapped_start,
            profile.matrix,
            profile.gap_open,
            profile.gap_extend,
            profile.gap_xdrop,
            scratch,
        );
        if gapped.score >= profile.gapped_cutoff[context] {
            let saved = GappedHsp {
                frame: frame.frame,
                score: gapped.score,
                q_start: gapped.query_start,
                q_end: gapped.query_stop,
                q_gapped_start: i32::try_from(q_gapped_start)?,
                s_start: gapped.subject_start,
                s_end: gapped.subject_stop,
                s_gapped_start: i32::try_from(s_gapped_start)?,
            };
            let saved_tree_hsp = TreeHsp {
                query_offset: saved.q_start,
                query_end: saved.q_end,
                subject_offset: saved.s_start,
                subject_end: saved.s_end,
                score: saved.score,
                query_frame: 0,
                query_length,
                query_context_offset: context_offset,
                subject_frame_sign,
            };
            if !tree.try_add_hsp(saved_tree_hsp, context_offset) {
                return Err(short_chunk_tree_refusal(
                    x_subject_oid(subject),
                    "TBLASTN",
                    chunk.ncbi_length,
                ));
            }
            results.push((context, saved));
        }
    }
    Ok(purge_preliminary_common_endpoints(results))
}

/// Database OID label of `target` for refusal messages (`"0"` for in-memory subjects).
fn x_subject_oid(target: super::x_subject_store::XPrelimSubject<'_>) -> String {
    match target {
        super::x_subject_store::XPrelimSubject::Stream(stream, ordinal, _) => {
            stream.seqdb().oid(ordinal as u32).to_string()
        }
        super::x_subject_store::XPrelimSubject::Nt(_) => "0".to_string(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:572-586
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2873-2911
/// Refuses a hard-chunk merge whose overlap NCBI reads past its list of ranges
/// (`overlap_past_list`) when `merge_joins` would enter NCBI's join loop.
#[cold]
#[inline(never)]
fn x_check_hard_merge_overlap(
    prelim_subject: super::x_subject_store::XPrelimSubject<'_>,
    chunk: TranslatedChunk,
    per_frame: &[(usize, GappedHsp)],
    unit_hsps: &[(usize, GappedHsp)],
) -> Result<()> {
    let hard = crate::core::subject_split::HardChunk {
        offset: chunk.offset,
        length: chunk.ncbi_length,
        residual: 0,
        full_range: chunk.full_range,
        range: 0,
        overlap: chunk.overlap,
        overlap_past_list: true,
    };
    let offset = i32::try_from(chunk.offset)?;
    if crate::algorithm::common::hard_chunk_scan::merge_joins(
        &hard,
        per_frame.iter().map(|(_, h)| h.s_end),
        unit_hsps.iter().map(|(_, h)| h.s_start + offset),
    ) {
        return Err(
            crate::algorithm::common::hard_chunk_scan::short_chunk_merge_refusal(
                x_subject_oid(prelim_subject),
                "TBLASTN",
            ),
        );
    }
    Ok(())
}

#[cfg(test)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum PreliminaryEvent {
    Candidates(i8, usize, Vec<Seed>),
    Initial(i8, usize, Vec<InitHsp>),
    Gapped(i8, usize, Vec<(usize, GappedHsp)>),
    Purged(i8, usize, Vec<(usize, GappedHsp)>),
    Merged(i8, usize, Vec<(usize, GappedHsp)>),
    Appended(i8, Vec<(usize, GappedHsp)>),
}

// No NCBI counterpart: the test probe below records the values the unit path computes at the
// steps of blast_engine.c:478-586,804-850 (WordFinder, GetGappedScore, endpoint purge,
// Blast_HSPListsMerge, Blast_HSPListAppend), so that unit tests can compare them with NCBI traces
// stored under docs/evidence. It only reads (and, through `at_purge_input`, edits the purge input
// of) what the search computes; the search does not call it.
/// Edits the gapped list of one (frame, chunk) unit before the endpoint purge: frame, chunk
/// offset, list. For comparison-only artificial NCBI input states.
#[cfg(test)]
pub(super) type PurgeInputHook =
    Box<dyn FnMut(i8, usize, &mut Vec<(usize, GappedHsp)>) -> Result<()>>;

/// The probe of `x_preliminary_protein_hsps_by_unit` (test builds, this thread, serial loop only).
#[cfg(test)]
#[derive(Default)]
struct PreliminaryProbe {
    /// The events in the serial harness's order, filled by the replay.
    events: Vec<PreliminaryEvent>,
    /// The events of each unit (index = unit), recorded when the unit runs and emitted by the
    /// replay in unit order, before that unit's merge.
    unit_events: Vec<Vec<PreliminaryEvent>>,
    at_purge_input: Option<PurgeInputHook>,
}

#[cfg(test)]
thread_local! {
    static PRELIMINARY_PROBE: std::cell::RefCell<Option<PreliminaryProbe>> =
        const { std::cell::RefCell::new(None) };
}

/// Runs `f` on the probe when this thread has set one.
#[cfg(test)]
fn with_probe<R>(f: impl FnOnce(&mut PreliminaryProbe) -> R) -> Option<R> {
    PRELIMINARY_PROBE.with_borrow_mut(|probe| probe.as_mut().map(f))
}

/// A unit starts: its buffer holds an empty `Candidates` event that the filled scan slices extend.
#[cfg(test)]
fn probe_begin_unit(unit: usize, frame: i8, chunk_offset: usize) {
    with_probe(|probe| {
        if probe.unit_events.len() <= unit {
            probe.unit_events.resize_with(unit + 1, Vec::new);
        }
        probe.unit_events[unit] = vec![PreliminaryEvent::Candidates(
            frame,
            chunk_offset,
            Vec::new(),
        )];
    });
}

/// One filled offset-pair slice of the unit's scan, in scan order.
#[cfg(test)]
fn probe_seeds(
    unit: usize,
    frame: i8,
    chunk_offset: u32,
    filled: &[crate::algorithm::tblastx::blast_aascan::BlastOffsetPair],
) {
    with_probe(|probe| {
        if let Some(PreliminaryEvent::Candidates(_, _, seeds)) = probe.unit_events[unit].first_mut()
        {
            seeds.extend(filled.iter().map(|pair| Seed {
                frame,
                chunk_offset,
                query_offset: pair.q_off,
                subject_offset: pair.s_off,
            }));
        }
    });
}

/// Appends an event to the unit's buffer.
#[cfg(test)]
fn probe_event(unit: usize, event: impl FnOnce() -> PreliminaryEvent) {
    with_probe(|probe| probe.unit_events[unit].push(event()));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:525-552:
// `GetGappedScore` returns the HSP list; immediately afterward the same list is
// passed to `Blast_HSPListPurgeHSPsWithCommonEndpoints` with endpoint purging enabled.
// The hook edits that input for comparison-only artificial NCBI input states.
/// Records `Gapped`, applies the purge-input hook, and returns the (edited) gapped list for the
/// search's purge.
#[cfg(test)]
fn probe_gapped(
    unit: usize,
    frame: i8,
    chunk_offset: usize,
    mut gapped: Vec<(usize, GappedHsp)>,
) -> Result<Vec<(usize, GappedHsp)>> {
    // A hook can panic (a test assertion); the probe is taken out while it runs.
    let Some(mut probe) = PRELIMINARY_PROBE.with_borrow_mut(Option::take) else {
        return Ok(gapped);
    };
    probe.unit_events[unit].push(PreliminaryEvent::Gapped(
        frame,
        chunk_offset,
        gapped.clone(),
    ));
    let edited = match probe.at_purge_input.as_mut() {
        Some(hook) => hook(frame, chunk_offset, &mut gapped),
        None => Ok(()),
    };
    PRELIMINARY_PROBE.with_borrow_mut(|slot| *slot = Some(probe));
    edited.map(|()| gapped)
}

/// The replay reaches a unit: its buffer joins the events.
#[cfg(test)]
fn probe_replay_unit(unit: usize) {
    with_probe(|probe| {
        if let Some(buffer) = probe.unit_events.get_mut(unit) {
            let buffer = std::mem::take(buffer);
            probe.events.extend(buffer);
        }
    });
}

/// The replay records one event (a merge or a frame's append).
#[cfg(test)]
fn probe_push(event: impl FnOnce() -> PreliminaryEvent) {
    with_probe(|probe| probe.events.push(event()));
}

/// The probe needs the serial loop (units run in order on this thread).
#[cfg(test)]
fn probe_require_serial(parallel: bool) {
    if parallel && with_probe(|_| ()).is_some() {
        panic!("the preliminary-stage test probe needs the serial unit loop (pool = None)");
    }
}

/// Takes the probe away again when the traced call ends, also when it fails or panics.
#[cfg(test)]
struct ProbeGuard;

#[cfg(test)]
impl Drop for ProbeGuard {
    fn drop(&mut self) {
        // No NCBI counterpart: test-only cleanup. `try_with` and `try_borrow_mut`: no second panic.
        let _ = PRELIMINARY_PROBE.try_with(|probe| {
            if let Ok(mut probe) = probe.try_borrow_mut() {
                *probe = None;
            }
        });
    }
}

// No NCBI counterpart: test helper. The exact search path (the unit path, no probe, no pool) on
// letters; the tests check its list against NCBI traces stored under docs/evidence.
#[cfg(test)]
pub(super) fn preliminary_hsps_of_letters(
    queries: &[&[u8]],
    subject: &[u8],
    db_gencode: u8,
    profile: PreliminaryProfile<'_>,
) -> Result<Vec<(usize, GappedHsp)>> {
    preliminary_protein_hsps_for_search(
        queries,
        SubjectNt::Letters(subject).into(),
        db_gencode,
        profile,
        None,
        &XQuerySide::new(),
        None,
        super::score::MatrixScore,
    )
}

// No NCBI counterpart: test helper. The same search path with the probe set: the list and the
// events of blast_engine.c:478-586,804-850 in NCBI order (per unit Candidates, Initial, Gapped,
// Purged; Merged and Appended in the replay). `at_purge_input` may edit a unit's gapped list
// before the endpoint purge (an artificial NCBI input state); it is `'static`, so a test that
// counts calls shares an `Rc<Cell<usize>>`.
#[cfg(test)]
pub(super) fn preliminary_hsps_traced(
    queries: &[&[u8]],
    subject: &[u8],
    db_gencode: u8,
    profile: PreliminaryProfile<'_>,
    at_purge_input: Option<PurgeInputHook>,
) -> Result<(Vec<(usize, GappedHsp)>, Vec<PreliminaryEvent>)> {
    preliminary_hsps_traced_of(
        queries,
        SubjectNt::Letters(subject),
        db_gencode,
        profile,
        at_purge_input,
    )
}

// No NCBI counterpart: test helper. `preliminary_hsps_traced` for the letters or the packed store
// of a subject (a packed one is searched as every `-db` sequence is).
#[cfg(test)]
fn preliminary_hsps_traced_of(
    queries: &[&[u8]],
    subject: SubjectNt<'_>,
    db_gencode: u8,
    profile: PreliminaryProfile<'_>,
    at_purge_input: Option<PurgeInputHook>,
) -> Result<(Vec<(usize, GappedHsp)>, Vec<PreliminaryEvent>)> {
    let _guard = ProbeGuard;
    PRELIMINARY_PROBE.with_borrow_mut(|probe| {
        *probe = Some(PreliminaryProbe {
            at_purge_input,
            ..PreliminaryProbe::default()
        });
    });
    let hsps = preliminary_protein_hsps_for_search(
        queries,
        subject.into(),
        db_gencode,
        profile,
        None,
        &XQuerySide::new(),
        None,
        super::score::MatrixScore,
    )?;
    let events = with_probe(|probe| std::mem::take(&mut probe.events)).unwrap_or_default();
    Ok((hsps, events))
}

// No NCBI counterpart: test helper. The probe events of one `preliminary_hsps_traced` call, sorted
// by stage: the offset pairs of the scans (blast_engine.c:478-500, aa_ungapped.c:492-516), the
// initial HSPs of the word finder (aa_ungapped.c:200-234), the GetGappedScore lists
// (blast_engine.c:525-552) with their chunk offsets, and the append snapshots (blast_engine.c:
// 840-850), each in event (= unit) order. It computes nothing.
#[cfg(test)]
pub(super) struct PreliminaryStages {
    /// The final list of the search.
    pub list: Vec<(usize, GappedHsp)>,
    /// Every `Candidates` event, concatenated.
    pub seeds: Vec<Seed>,
    /// Every `Initial` event, concatenated.
    pub initial: Vec<InitHsp>,
    /// Every `Gapped` event as (context, hsp, chunk offset).
    pub gapped: Vec<(usize, GappedHsp, usize)>,
    /// One snapshot per `Appended` event (one per subject frame).
    pub appended: Vec<Vec<(usize, GappedHsp)>>,
}

// No NCBI counterpart: test helper (see `PreliminaryStages`).
#[cfg(test)]
pub(super) fn preliminary_stages(
    queries: &[&[u8]],
    subject: &[u8],
    db_gencode: u8,
    profile: PreliminaryProfile<'_>,
) -> Result<PreliminaryStages> {
    preliminary_stages_of(queries, SubjectNt::Letters(subject), db_gencode, profile)
}

// No NCBI counterpart: test helper (see `PreliminaryStages`), for the letters or the packed store
// of a subject.
#[cfg(test)]
pub(super) fn preliminary_stages_of(
    queries: &[&[u8]],
    subject: SubjectNt<'_>,
    db_gencode: u8,
    profile: PreliminaryProfile<'_>,
) -> Result<PreliminaryStages> {
    let (list, events) = preliminary_hsps_traced_of(queries, subject, db_gencode, profile, None)?;
    let mut stages = PreliminaryStages {
        list,
        seeds: Vec::new(),
        initial: Vec::new(),
        gapped: Vec::new(),
        appended: Vec::new(),
    };
    for event in events {
        match event {
            PreliminaryEvent::Candidates(_, _, seeds) => stages.seeds.extend(seeds),
            PreliminaryEvent::Initial(_, _, hits) => stages.initial.extend(hits),
            PreliminaryEvent::Gapped(_, offset, hsps) => stages.gapped.extend(
                hsps.into_iter()
                    .map(|(context, hsp)| (context, hsp, offset)),
            ),
            PreliminaryEvent::Appended(_, hsps) => stages.appended.push(hsps),
            PreliminaryEvent::Purged(..) | PreliminaryEvent::Merged(..) => {}
        }
    }
    Ok(stages)
}

// No NCBI counterpart: test statics. Per-context cutoffs of INT4_MAX (nothing is saved) and
// per-context x-drops of 0, for the test profiles below.
#[cfg(test)]
static TEST_NEVER: [i32; 4] = [i32::MAX; 4];
#[cfg(test)]
static TEST_ZERO: [i32; 4] = [0; 4];

// No NCBI counterpart: test helper. A profile for tests that read the word finder's initial HSPs
// only: gap 11/1/38, and a gapped cutoff of INT4_MAX, so that GetGappedScore (blast_engine.c:
// 489-552) runs on the saved initial HSPs and saves none; `Initial` is recorded before it.
#[cfg(test)]
pub(super) fn word_finder_profile<'a>(
    seg: Option<&'a SegParams>,
    threshold: i32,
    window: i32,
    word_xdrop: &'a [i32],
    word_cutoff: &'a [i32],
    mask_lowercase: bool,
    matrix: ScoringMatrix,
    word_size: usize,
) -> PreliminaryProfile<'a> {
    assert!(word_xdrop.len() <= TEST_NEVER.len());
    PreliminaryProfile {
        seg,
        soft_masking: false,
        threshold,
        window,
        word_xdrop,
        word_cutoff,
        mask_lowercase,
        matrix,
        word_size,
        gap_open: 11,
        gap_extend: 1,
        gap_xdrop: 38,
        gapped_cutoff: &TEST_NEVER[..word_xdrop.len()],
        hsp_num_max: i32::MAX as usize,
        gapped: true,
        subject_besthit: false,
    }
}

// No NCBI counterpart: test helper. A profile for tests that read the scan's offset pairs only: a
// word cutoff of INT4_MAX saves no initial HSP, so GetGappedScore is never called
// (blast_engine.c:500); `Candidates` does not depend on the window or the x-drop.
#[cfg(test)]
pub(super) fn scan_profile(
    queries: usize,
    seg: Option<&SegParams>,
    threshold: i32,
    mask_lowercase: bool,
    matrix: ScoringMatrix,
    word_size: usize,
) -> PreliminaryProfile<'_> {
    assert!(queries <= TEST_NEVER.len());
    word_finder_profile(
        seg,
        threshold,
        40,
        &TEST_ZERO[..queries],
        &TEST_NEVER[..queries],
        mask_lowercase,
        matrix,
        word_size,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4742-4772
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2371-2394
/// The HSPs of an ungapped search's (frame, chunk) unit (`BLAST_GetUngappedHSPList`): each initial
/// HSP of `hits` (sorted by `XTwoHit::sort`, as NCBI's word finder leaves them) in its query
/// context (the context of its word's query offset), with the query offsets made the context's
/// and its word's offsets as the gapped start, then sorted by score. The subject offsets stay the
/// chunk's (`merge_adjusted_chunk` adds the chunk's offset). The search's `hsp_num_max` is
/// `INT4_MAX` (`PreliminaryProfile::hsp_num_max`), so every HSP is saved.
fn ungapped_chunk_hsps(
    context_offsets: &[i32],
    frame: i8,
    hits: &[InitHsp],
) -> Vec<(usize, GappedHsp)> {
    let mut hsps: Vec<(usize, GappedHsp)> = hits
        .iter()
        .map(|hit| {
            let q_off = hit.q_seed as i32;
            let context = context_offsets.partition_point(|&offset| offset <= q_off) - 1;
            let query_start = context_offsets[context];
            let q_start = hit.q_start - query_start;
            (
                context,
                GappedHsp {
                    frame,
                    score: hit.score,
                    q_start,
                    q_end: hit.length + q_start,
                    q_gapped_start: q_off - query_start,
                    s_start: hit.s_start,
                    s_end: hit.length + hit.s_start,
                    s_gapped_start: hit.s_seed as i32,
                },
            )
        })
        .collect();
    sort_preliminary_by_score(&mut hsps);
    hsps
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1382:
// ScoreCompareHSPs sorts score DESC, then subject offset ASC,
// subject end DESC, query offset ASC, query end DESC.
fn sort_preliminary_by_score(hsps: &mut [(usize, GappedHsp)]) {
    if hsps
        .windows(2)
        .any(|pair| compare_gapped_score(&pair[0].1, &pair[1].1).is_gt())
    {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1382
        // The pinned NCBI BLAST+ runs on glibc, whose qsort is a stable merge sort.
        // ScoreCompareHSPs ignores the subject frame, so HSPs of other frames at the
        // same protein offsets keep their order (frames 1, 2, 3, -1, -2, -3 as
        // appended); the endpoint purge and the containment test keep the first (TN-5).
        hsps.sort_by(|a, b| compare_gapped_score(&a.1, &b.1));
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2268-2379,2455-2537:
// Sorts `hsp_array` by `s_QueryOffsetCompareHSPs` and deletes duplicate starts,
// then sorts by `s_QueryEndCompareHSPs` and deletes duplicate ends.
// Context and subject frame must also agree before deletion.
fn purge_preliminary_common_endpoints(
    mut incoming: Vec<(usize, GappedHsp)>,
) -> Vec<(usize, GappedHsp)> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2268-2323,2455-2492:
    // `s_QueryOffsetCompareHSPs` returns 0 for identical score,
    // query/subject starts and ends. The pinned NCBI oracle retains the first
    // score-355 HSP's gapped_start (192, 47944) on that exact comparator tie.
    incoming.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then(a.1.q_start.cmp(&b.1.q_start))
            .then(a.1.s_start.cmp(&b.1.s_start))
            .then(b.1.score.cmp(&a.1.score))
            .then(b.1.q_end.cmp(&a.1.q_end))
            .then(b.1.s_end.cmp(&a.1.s_end))
    });
    let mut i = 0;
    while i + 1 < incoming.len() {
        let a = incoming[i];
        let b = incoming[i + 1];
        if a.0 == b.0
            && a.1.frame == b.1.frame
            && a.1.q_start == b.1.q_start
            && a.1.s_start == b.1.s_start
        {
            incoming.remove(i + 1);
        } else {
            i += 1;
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2330-2390,2497-2537:
    // s_QueryEndCompareHSPs also returns 0 for these equal-endpoint HSPs;
    // retain their arrival order through the second purge pass.
    incoming.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then(a.1.q_end.cmp(&b.1.q_end))
            .then(a.1.s_end.cmp(&b.1.s_end))
            .then(b.1.score.cmp(&a.1.score))
            .then(b.1.q_start.cmp(&a.1.q_start))
            .then(b.1.s_start.cmp(&a.1.s_start))
    });
    i = 0;
    while i + 1 < incoming.len() {
        let a = incoming[i];
        let b = incoming[i + 1];
        if a.0 == b.0 && a.1.frame == b.1.frame && a.1.q_end == b.1.q_end && a.1.s_end == b.1.s_end
        {
            incoming.remove(i + 1);
        } else {
            i += 1;
        }
    }
    sort_preliminary_by_score(&mut incoming);
    incoming
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2762-2864:
// If there is no existing combined list, adopt incoming; otherwise combine by `ScoreCompareHSPs`
// up to hsp_num_max, preferring the existing list on a comparator tie.
fn append_preliminary_hsps(
    mut combined: Vec<(usize, GappedHsp)>,
    mut incoming: Vec<(usize, GappedHsp)>,
    hsp_num_max: usize,
) -> Vec<(usize, GappedHsp)> {
    if incoming.is_empty() {
        return combined;
    }
    if combined.is_empty() {
        return incoming;
    }
    let retained = combined
        .len()
        .saturating_add(incoming.len())
        .min(hsp_num_max);
    if retained == combined.len() + incoming.len() {
        combined.extend(incoming);
        sort_preliminary_by_score(&mut combined);
        return combined;
    }
    sort_preliminary_by_score(&mut combined);
    sort_preliminary_by_score(&mut incoming);
    let mut selected = Vec::with_capacity(retained);
    let mut old = 0;
    let mut new = 0;
    while selected.len() < retained {
        if old < combined.len()
            && (new >= incoming.len()
                || compare_gapped_score(&combined[old].1, &incoming[new].1).is_le())
        {
            selected.push(combined[old]);
            old += 1;
        } else {
            selected.push(incoming[new]);
            new += 1;
        }
    }
    selected
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:572-586:
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3037-3050:
// Adjust subject offset, end and gapped_start before merging.
// `allow_gap` is `score_options->gapped_calculation` (blast_engine.c:584): an ungapped search
// merges the HSPs of a diagonal only.
fn merge_adjusted_chunk(
    per_frame: Vec<(usize, GappedHsp)>,
    mut incoming: Vec<(usize, GappedHsp)>,
    chunk: TranslatedChunk,
    hsp_num_max: usize,
    allow_gap: bool,
) -> Result<Vec<(usize, GappedHsp)>> {
    let offset = i32::try_from(chunk.offset)?;
    for (_, hsp) in &mut incoming {
        hsp.s_start = hsp
            .s_start
            .checked_add(offset)
            .context("subject offset overflow")?;
        hsp.s_end = hsp
            .s_end
            .checked_add(offset)
            .context("subject end overflow")?;
        hsp.s_gapped_start = hsp
            .s_gapped_start
            .checked_add(offset)
            .context("subject gapped start overflow")?;
    }
    // The chunk's overlap (`TranslatedChunk::overlap`: 0 where it starts the hard range at
    // `hm_index`, else DBSEQ_CHUNK_OVERLAP).
    let overlap = i32::try_from(chunk.overlap)?;
    Ok(merge_subject_chunk_hsps(
        per_frame,
        incoming,
        offset,
        overlap,
        hsp_num_max,
        allow_gap,
    ))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1465-1537:
// Compares start and end diagonals (query minus subject at each endpoint).
// If either new endpoint falls inside the old HSP, computes the combined score density
// as the sum of scores divided by the sum of query lengths, and takes the larger of
// the truncated density-scaled merged query length and the old score.
pub(super) fn merge_two_chunk_hsps(old: &mut GappedHsp, new: &GappedHsp, allow_gap: bool) -> bool {
    if !allow_gap && old.s_start - new.s_start - old.q_start + new.q_start != 0 {
        return false;
    }
    if old.frame != new.frame {
        return false;
    }
    let contained =
        |q: i32, s: i32| old.q_start <= q && q <= old.q_end && old.s_start <= s && s <= old.s_end;
    if !contained(new.q_start, new.s_start) && !contained(new.q_end, new.s_end) {
        return false;
    }
    let density = f64::from(old.score + new.score)
        / f64::from((old.q_end - old.q_start) + (new.q_end - new.q_start));
    old.q_start = old.q_start.min(new.q_start);
    old.s_start = old.s_start.min(new.s_start);
    old.q_end = old.q_end.max(new.q_end);
    old.s_end = old.s_end.max(new.s_end);
    if new.score > old.score {
        old.q_gapped_start = new.q_gapped_start;
        old.s_gapped_start = new.s_gapped_start;
        old.score = new.score;
    }
    old.score = ((density * f64::from(old.q_end - old.q_start)) as i32).max(old.score);
    true
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:572-586:
// Uses zero overlap when the chunk starts at the hard left boundary, otherwise `DBSEQ_CHUNK_OVERLAP`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2857-3035:
// move overlap HSPs to the front, match equal contexts whose diagonals
// differ by less than 10, merge, purge the consumed new HSPs, combine by score.
// Input HSP coordinates must already include the chunk offset.
#[allow(dead_code)]
fn merge_subject_chunk_hsps(
    mut old: Vec<(usize, GappedHsp)>,
    mut incoming: Vec<(usize, GappedHsp)>,
    split_offset: i32,
    overlap: i32,
    hsp_num_max: usize,
    allow_gap: bool,
) -> Vec<(usize, GappedHsp)> {
    if incoming.is_empty() {
        return old;
    }
    if old.is_empty() {
        return incoming;
    }
    let mut old_overlap = 0;
    for index in 0..old.len() {
        if old[index].1.s_end > split_offset {
            old.swap(old_overlap, index);
            old_overlap += 1;
        }
    }
    let mut new_overlap = 0;
    for index in 0..incoming.len() {
        if incoming[index].1.s_start < split_offset + overlap {
            incoming.swap(new_overlap, index);
            new_overlap += 1;
        }
    }
    let mut incoming: Vec<Option<(usize, GappedHsp)>> = incoming.into_iter().map(Some).collect();
    for old_index in 0..old_overlap {
        for new_index in 0..new_overlap {
            let Some((context, candidate)) = incoming[new_index] else {
                continue;
            };
            if old[old_index].0 != context {
                continue;
            }
            let old_end_diag = old[old_index].1.q_end - old[old_index].1.s_end;
            let new_start_diag = candidate.q_start - candidate.s_start;
            if (old_end_diag - new_start_diag).abs() < 10
                && merge_two_chunk_hsps(&mut old[old_index].1, &candidate, allow_gap)
            {
                incoming[new_index] = None;
            }
        }
    }
    let mut incoming: Vec<_> = incoming.into_iter().flatten().collect();
    let retained = old.len().saturating_add(incoming.len()).min(hsp_num_max);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2762-2807, 1330-1382:
    // on no cap, append then SortByScore if out of order;
    // on cap, sort both lists and merge score order, preferring old on ties.
    let sorted = |hsps: &mut Vec<(usize, GappedHsp)>| {
        if hsps
            .windows(2)
            .any(|pair| compare_gapped_score(&pair[0].1, &pair[1].1).is_gt())
        {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1382
            // Stable, as glibc's qsort under the pinned NCBI BLAST+ (TN-5).
            hsps.sort_by(|a, b| compare_gapped_score(&a.1, &b.1));
        }
    };
    if retained == old.len() + incoming.len() {
        old.extend(incoming);
        sorted(&mut old);
        return old;
    }
    sorted(&mut old);
    sorted(&mut incoming);
    let mut selected = Vec::with_capacity(retained);
    let mut old_index = 0;
    let mut new_index = 0;
    while selected.len() < retained {
        if old_index < old.len()
            && (new_index >= incoming.len()
                || compare_gapped_score(&old[old_index].1, &incoming[new_index].1).is_le())
        {
            selected.push(old[old_index]);
            old_index += 1;
        } else {
            selected.push(incoming[new_index]);
            new_index += 1;
        }
    }
    selected
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1268-1317:
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1156-1228:
#[derive(Default)]
struct TargetFrameTranslation {
    start: usize,
    stop: usize,
    /// The subject position of `sequence[0]`: -1 for the left fence of a window that
    /// starts at the first residue.
    base: isize,
    sequence: Vec<u8>,
    /// The request (`get`'s window of bases and whether it is fenced) that made `sequence`. The
    /// same request again (the report counts after the search's counts of one HSP) reads it as
    /// it is: a translation depends only on its window. No NCBI counterpart.
    made_from: Option<(usize, usize, bool)>,
}

/// The index in a translated window (`TargetTranslation::get`, whose subject position of the
/// first byte is `base`) of subject position `position`, or `None` before the window.
pub(super) fn window_index(position: i32, base: isize) -> Option<usize> {
    usize::try_from(position as isize - base).ok()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1284;
// c++/src/algo/blast/core/blast_hits.c:1162-1164,1234-1235
// NCBI skips the translation of a subject (its searched range, `-subject_loc`) of fewer than 12
// bases (`stop < n / 3 - 3` is false, `core::blast_hits::target_translation_never_made`) and
// returns the frame's translation that was never made: NULL + 1. Every caller reads from there
// next: the traceback's gapped start and ALIGN_EX
// (blast_traceback.c:429-440,509-512), the identities after the composition-based redo
// (blast_kappa.c:527-528, after its E-value purge, :3677-3700), the ungapped re-evaluation
// (blast_hits.c:2700-2705; LOSAT's TBLASTN `-ungapped` and TBLASTX refuse there,
// PD-LOSAT-NCBI-DEFECTS 1.79 and 1.51). NCBI BLAST+ 2.17.0 dies there with
// SIGSEGV (probes on ws, 2026-11-17: `-subject`, `-subject_loc` and `-db`, frames + and -, a matrix
// query at `-comp_based_stats` 0 to 3, `-in_pssm`): no result is defined (PD-LOSAT-NCBI-DEFECTS
// 1.46, rule 3). A subject that leaves no HSP to that point (cut by a cutoff, the redo's purge or
// the hit list; at the default word size 3, any of 9 or 10 bases in the probes; with `-word_size 2`
// they reach it too) lets the search go on, in NCBI and LOSAT alike. NCBI dies in the search of a
// query batch, after it has printed the reports of the batches before (tblastn_app.cpp:275-302):
// the refusal comes after them (`blastinput::app::ErrorAfterBatchReports`, line RG item
// X-short-bat).
/// The refusal of a read of the translation of a frame of a subject of fewer than 12 bases
/// (`TargetTranslation::get`): NCBI BLAST+ never makes that translation and reads from address 1
/// on (NULL + 1). The message is TBLASTX's too (`core::blast_hits::x_untranslated_frame_refusal`).
#[cold]
#[inline(never)]
fn x_traceback_of_an_untranslated_frame() -> anyhow::Error {
    anyhow::Error::new(crate::blastinput::app::ErrorAfterBatchReports::new(
        crate::core::blast_hits::x_untranslated_frame_refusal("TBLASTN", "the traceback"),
    ))
}

pub(super) struct TargetTranslation<'a> {
    subject: SubjectNt<'a>,
    code: &'a GeneticCode,
    frames: [TargetFrameTranslation; 6],
    /// EXPERIMENT (LOSAT_X_TBNSTORE): the letters of the window being translated, when the
    /// subject is packed.
    window: Vec<u8>,
    /// The bases `[from, to)` whose letters `window` holds.
    window_range: Option<(usize, usize)>,
    /// The residue of every codon of the code (`x_codon_residues`), taken at the first
    /// translation.
    residues: Option<std::rc::Rc<[u8; 4096]>>,
}

impl<'a> TargetTranslation<'a> {
    #[cfg(test)]
    pub(super) fn new(subject: &'a [u8], code: &'a GeneticCode) -> Self {
        Self::of(SubjectNt::Letters(subject), code)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1197-1201,1216-1218
    // NCBI translates `nucl_length` nucleotides from `nucl_shift`; a packed subject (EXPERIMENT,
    // LOSAT_X_TBNSTORE) gives the letters of that window when it is translated.
    /// The translations of `subject`, letters or (EXPERIMENT, LOSAT_X_TBNSTORE) packed.
    pub(super) fn of(subject: SubjectNt<'a>, code: &'a GeneticCode) -> Self {
        Self {
            subject,
            code,
            frames: std::array::from_fn(|_| TargetFrameTranslation::default()),
            window: Vec::new(),
            window_range: None,
            residues: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1154-1228:
    // Translates the requested subject window via `BLAST_GetTranslation`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:428-455:
    pub(super) fn get(
        &mut self,
        frame: i8,
        offset: i32,
        end: i32,
    ) -> Result<(&[u8], usize, isize)> {
        let context = if frame > 0 { frame - 1 } else { 2 - frame };
        let context = usize::try_from(context).context("invalid subject frame")?;
        if context >= self.frames.len() || frame == 0 {
            bail!("invalid subject frame {frame}");
        }
        let target = &mut self.frames[context];
        let n = self.subject.len();
        if target.start != 0 || target.stop < (n / 3).saturating_sub(3) {
            let (nucl_start, mut nucl_end) = if offset < 0 {
                (0, n)
            } else {
                (
                    usize::try_from(offset)?
                        .saturating_mul(3)
                        .saturating_sub(99),
                    usize::try_from(end)?
                        .saturating_mul(3)
                        .saturating_add(99)
                        .min(n),
                )
            };
            if n - nucl_end <= 21 {
                nucl_end = n;
            }
            let nucl_length = nucl_end - nucl_start;
            let translation_length = 1 + nucl_length / 3;
            let start_shift = nucl_start / 3;
            let request = (nucl_start, nucl_end, offset >= 0);
            if (start_shift < target.start || start_shift + translation_length > target.stop)
                && target.made_from != Some(request)
            {
                let nucl_shift = if frame < 0 { n - nucl_end } else { nucl_start };
                let window = (nucl_shift, nucl_shift + nucl_length);
                let letters = match self.subject {
                    SubjectNt::Letters(letters) => &letters[window.0..window.1],
                    // A packed subject (every `-db` sequence; with EXPERIMENT LOSAT_X_TBNSTORE every
                    // `-subject` record): the letters of a window are made once for consecutive
                    // reads of the same bases (an HSP's frame read again, or after a fence every
                    // frame, which reads the whole subject).
                    SubjectNt::Packed(store) => {
                        if self.window_range != Some(window) {
                            store.letters_into(window.0, window.1, &mut self.window);
                            self.window_range = Some(window);
                        }
                        &self.window[..]
                    }
                };
                let code = self.code;
                let residues = self.residues.get_or_insert_with(|| x_codon_residues(code));
                // NCBI translates only the HSP's frame (c++/src/algo/blast/core/blast_hits.c:1216-1218,
                // `BLAST_GetTranslation`); `translate_frame_into` writes its residues between the
                // window's ends (below) in the frame's buffer.
                target.sequence.clear();
                if offset >= 0 {
                    target.sequence.push(201);
                }
                let Some(aa_len) =
                    translate_frame_into(letters, frame, residues, &mut target.sequence)
                else {
                    *target = TargetFrameTranslation::default();
                    bail!("missing translated subject frame");
                };
                target.sequence.push(if offset >= 0 { 201 } else { 0 });
                target.start = start_shift;
                target.stop = start_shift + aa_len;
                target.made_from = Some(request);
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1224-1228
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1235:
                // The Rust slice begins at the left fence and `base` restores NCBI's
                // absolute subject coordinates without allocating the unused prefix
                // preceding a late partial window. A window from the first residue is
                // fenced too (its fence is at position -1): a traceback that reaches it
                // is redone with the full translation.
                target.base = if offset >= 0 {
                    start_shift as isize - 1
                } else {
                    0
                };
            }
        } else if target.sequence.is_empty() {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1162-1164,1234-1235: the frame never translated, NULL + 1
            // (`x_traceback_of_an_untranslated_frame`). Only a subject of fewer than 12 bases skips
            // the first translation; every translation leaves at least its closing byte.
            debug_assert!(
                crate::core::blast_hits::target_translation_never_made(n),
                "an empty translation of a subject of {n} bases"
            );
            return Err(x_traceback_of_an_untranslated_frame());
        }
        Ok((&target.sequence, target.stop, target.base))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:911-918
    // The reports translate the subject's codons again with the display's table (the rows of
    // `stage_e_report::aligned_sequences`), where an ambiguous codon whose residues are D/N, E/Q or
    // I/L is B, Z or J; the search's translation (`get`) reads X there.
    /// The residue that the reports show at residue `offset` of frame `frame` of the subject (its
    /// codon translated by `DisplayCode`), as an NCBIstdaa code; `buffer` holds the codon's
    /// letters of a packed subject.
    pub(super) fn shown_residue(
        &self,
        frame: i8,
        offset: usize,
        buffer: &mut Vec<u8>,
    ) -> Result<u8> {
        let n = self.subject.len();
        let first = 3 * offset + usize::from(frame.unsigned_abs()) - 1;
        anyhow::ensure!(first + 3 <= n, "TBLASTN shown residue outside the subject");
        // The codon's three bases, read from its start on the plus strand and from its end on the
        // minus strand (`DisplayCode::residue` of offset 0 in a frame of the same sign).
        let start = if frame > 0 { first } else { n - first - 3 };
        let codon = self.subject.fetch(start, start + 3, buffer);
        let display = crate::algorithm::tblastx::report::DisplayCode::new(self.code);
        let letter = display.residue(codon, frame.signum(), 0);
        Ok(crate::utils::matrix::aa_char_to_ncbistdaa(letter))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1648-1677
    // After a fence NCBI aligns every HSP of the list again in the whole translation of its frame
    // (`get` with offset -1: the residues, no byte before the first, a 0 after the last). No NCBI
    // counterpart: LOSAT's speed (a whole frame of a chromosome is millions of residues, translated
    // again for every list retried). The traceback reads a window of that translation with a fence
    // where the window ends before the translation does, and grows the window when it reads one.
    // Why the window's traceback is the whole translation's when it reads no fence: each extension
    // reads the subject from the HSP outward, one column further at a time, and a column past the
    // window is read only after the column of the fence before it; the band never reaches past
    // the window's letters without reading that column (`b_size` passes the window's last letter
    // only to a column the next row reads). So a traceback that reads no fence reads only residues
    // of the whole translation, the same cells as in it, and its subject length (the window's end)
    // stops no band that the whole length would let go on. The gapped start
    // (`blast_get_offsets_for_gapped_alignment_protein`) reads the HSP's residues and tests its
    // window of `HSP_MAX_WINDOW` residues against the subject's ends, and with at least that many
    // residues on either side of the HSP its reads and tests are those of the whole translation;
    // the identities read the alignment's residues.
    /// The residues `[from, to)` of the whole translation of `frame` (`get(frame, -1, _)`) around
    /// the HSP `[offset, end)` (`margin` residues on either side, at least `HSP_MAX_WINDOW`, cut at
    /// the translation's ends), as (residues, the whole translation's length, the subject position
    /// of the first byte). The residues end as the whole translation does where the window
    /// reaches its ends (no byte before its first residue, a 0 after its last) and with a fence
    /// elsewhere. A subject of fewer than 12 letters gives what `get(frame, -1, end)` gives (its
    /// refusal, `x_traceback_of_an_untranslated_frame`).
    pub(super) fn x_whole_window(
        &mut self,
        frame: i8,
        offset: i32,
        end: i32,
        margin: usize,
    ) -> Result<(&[u8], usize, isize)> {
        debug_assert!(margin >= HSP_MAX_WINDOW, "margin {margin}");
        let n = self.subject.len();
        // `get` translates nothing for a subject of fewer than 12 letters (`stop < n / 3 - 3` is
        // false) and refuses the read, as after a fence.
        if crate::core::blast_hits::target_translation_never_made(n) {
            return self.get(frame, -1, end);
        }
        let context = if frame > 0 { frame - 1 } else { 2 - frame };
        let context = usize::try_from(context).context("invalid subject frame")?;
        if context >= self.frames.len() || frame == 0 {
            bail!("invalid subject frame {frame}");
        }
        let length = frame_length(n, frame).context("missing translated subject frame")?;
        let from = usize::try_from(offset)?.saturating_sub(margin);
        let to = usize::try_from(end)?.saturating_add(margin).min(length);
        let target = &mut self.frames[context];
        if target.sequence.is_empty() || from < target.start || to > target.stop {
            // The letters of codons `[from, to)`: those of frame `i + 1` start at letter
            // `i + 3 * from`; codon k of frame `-(i + 1)` is the reverse complement of the letters
            // `[n - i - 3k - 3, n - i - 3k)`.
            let i = usize::from(frame.unsigned_abs()) - 1;
            let bases = if frame > 0 {
                (i + 3 * from, i + 3 * to)
            } else {
                (n - i - 3 * to, n - i - 3 * from)
            };
            let letters = match self.subject {
                SubjectNt::Letters(letters) => &letters[bases.0..bases.1],
                SubjectNt::Packed(store) => {
                    if self.window_range != Some(bases) {
                        store.letters_into(bases.0, bases.1, &mut self.window);
                        self.window_range = Some(bases);
                    }
                    &self.window[..]
                }
            };
            let code = self.code;
            let residues = self.residues.get_or_insert_with(|| x_codon_residues(code));
            target.sequence.clear();
            if from > 0 {
                target.sequence.push(201);
            }
            x_codons_into(letters, frame < 0, residues, &mut target.sequence);
            target.sequence.push(if to == length { 0 } else { 201 });
            target.start = from;
            target.stop = to;
            target.made_from = None;
            target.base = if from > 0 { from as isize - 1 } else { 0 };
        }
        Ok((&target.sequence, length, target.base))
    }
}

/// The residues on either side of an HSP in the windows of the whole translation
/// (`x_whole_window`) that the traceback reads after a fence; when one is too narrow, the list is
/// traced again with windows 8 times as wide.
const X_WHOLE_MARGIN: usize = 1 << 10;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1382:
// compare score DESC, subject.offset ASC, subject.end DESC,
// query.offset ASC, query.end DESC; sort only if an adjacent pair is out of order.
fn compare_gapped_score(a: &GappedHsp, b: &GappedHsp) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then(a.s_start.cmp(&b.s_start))
        .then(b.s_end.cmp(&a.s_end))
        .then(a.q_start.cmp(&b.q_start))
        .then(b.q_end.cmp(&a.q_end))
}

fn sort_gapped_score_if_needed(hsps: &mut [GappedHsp]) {
    if hsps
        .windows(2)
        .any(|pair| compare_gapped_score(&pair[0], &pair[1]).is_gt())
    {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1382
        // Stable, as glibc's qsort under the pinned NCBI BLAST+ (TN-5).
        hsps.sort_by(compare_gapped_score);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:537,817-838,931-995:
// Delegates containment checks to `s_HSPIsContained` using the node's HSP and bounds.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:352-360:
fn traceback_tree_hsp(hsp: &GappedHsp, query_length: i32) -> TreeHsp {
    TreeHsp {
        query_offset: hsp.q_start,
        query_end: hsp.q_end,
        subject_offset: hsp.s_start,
        subject_end: hsp.s_end,
        score: hsp.score,
        query_frame: 0,
        query_length,
        query_context_offset: 0,
        subject_frame_sign: i32::from(hsp.frame.signum()),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2268-2379,2455-2537:
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1382:
// Orders by descending score, then ascending subject offset, descending subject end,
// ascending query offset, and descending query end.
// This internal Stage C list has one protein query context; multi-query
// context offsets remain a separate Stage C boundary.
fn purge_traceback_common_endpoints(mut hsps: Vec<TracebackOwnedHsp>) -> Vec<TracebackOwnedHsp> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2487
    // glibc's qsort under the pinned NCBI BLAST+ is a stable merge sort. The
    // comparator ignores the subject frame and only the next HSP is compared, so
    // the order of other frames at equal offsets decides what is purged (TN-5).
    hsps.sort_by(|a, b| {
        a.hsp
            .q_start
            .cmp(&b.hsp.q_start)
            .then(a.hsp.s_start.cmp(&b.hsp.s_start))
            .then(b.hsp.score.cmp(&a.hsp.score))
            .then(b.hsp.q_end.cmp(&a.hsp.q_end))
            .then(b.hsp.s_end.cmp(&a.hsp.s_end))
    });
    let mut i = 0;
    while i + 1 < hsps.len() {
        let a = hsps[i].hsp;
        let b = hsps[i + 1].hsp;
        if a.frame == b.frame && a.q_start == b.q_start && a.s_start == b.s_start {
            hsps.remove(i + 1);
        } else {
            i += 1;
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2504
    // Stable, as glibc's qsort (TN-5).
    hsps.sort_by(|a, b| {
        a.hsp
            .q_end
            .cmp(&b.hsp.q_end)
            .then(a.hsp.s_end.cmp(&b.hsp.s_end))
            .then(b.hsp.score.cmp(&a.hsp.score))
            .then(b.hsp.q_start.cmp(&a.hsp.q_start))
            .then(b.hsp.s_start.cmp(&a.hsp.s_start))
    });
    i = 0;
    while i + 1 < hsps.len() {
        let a = hsps[i].hsp;
        let b = hsps[i + 1].hsp;
        if a.frame == b.frame && a.q_end == b.q_end && a.s_end == b.s_end {
            hsps.remove(i + 1);
        } else {
            i += 1;
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1382:
    if hsps
        .windows(2)
        .any(|pair| compare_gapped_score(&pair[0].hsp, &pair[1].hsp).is_gt())
    {
        // Stable, as glibc's qsort under the pinned NCBI BLAST+ (TN-5).
        hsps.sort_by(|a, b| compare_gapped_score(&a.hsp, &b.hsp));
    }
    hsps
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:303-312,408-440,503-536,
// 709-721,1644-1684:
// When a fence sentinel is hit, refetches the whole subject and retries traceback.
// The NCBI endpoint purge and interval-tree pass follow traceback; ambiguity
// reevaluation has no HSPs in this non-greedy, `purge`-enabled profile.
#[allow(dead_code)]
pub(super) fn full_translation_traceback_blosum62(
    query: &[u8],
    subject: &[u8],
    db_gencode: u8,
    gapped: &[GappedHsp],
    gap_open: i32,
    gap_extend: i32,
    x_drop_final: i32,
) -> Result<Vec<(GappedHsp, bool)>> {
    full_translation_traceback_with_matrix(
        query,
        subject,
        db_gencode,
        gapped,
        ScoringMatrix::Blosum62,
        gap_open,
        gap_extend,
        x_drop_final,
        0.0,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign_priv.h:120
/// The window of residues in which the traceback's gapped start is looked for.
const HSP_MAX_WINDOW: usize = 11;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3248-3321:
// When the query span fits one window, both starts advance by `q_length / 2`; otherwise
// slide the window and keep the first position whose score is strictly greater than the
// best so far, then place the subject start on the same diagonal (`max_offset - q_start + s_start`).
pub(super) fn blast_get_offsets_for_gapped_alignment_protein(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
    matrix: ScoringMatrix,
) -> Option<(usize, usize)> {
    let q_length = q_end.checked_sub(q_start)?;
    let s_length = s_end.checked_sub(s_start)?;
    if q_end > query.len() || s_end > subject.len() {
        return None;
    }
    if q_length <= HSP_MAX_WINDOW {
        return Some((q_start + q_length / 2, s_start + q_length / 2));
    }
    // NCBI uses raw pointers; these are Rust bounds checks for the same reads.
    if s_start.checked_add(HSP_MAX_WINDOW)? > subject.len() || s_end < HSP_MAX_WINDOW {
        return None;
    }
    let mut score = 0i32;
    for index in 0..HSP_MAX_WINDOW {
        score += protein_score(matrix, query[q_start + index], subject[s_start + index]);
    }
    let mut max_score = score;
    let mut max_offset = q_start + HSP_MAX_WINDOW - 1;
    for index in (q_start + HSP_MAX_WINDOW)..(q_start + q_length.min(s_length)) {
        let subject_index = s_start + index - q_start;
        score -= protein_score(
            matrix,
            query[index - HSP_MAX_WINDOW],
            subject[subject_index - HSP_MAX_WINDOW],
        );
        score += protein_score(matrix, query[index], subject[subject_index]);
        if score > max_score {
            max_score = score;
            max_offset = index;
        }
    }
    if max_score > 0 {
        return Some((max_offset, max_offset - q_start + s_start));
    }
    score = 0;
    for index in 0..HSP_MAX_WINDOW {
        score += protein_score(
            matrix,
            query[q_end - HSP_MAX_WINDOW + index],
            subject[s_end - HSP_MAX_WINDOW + index],
        );
    }
    if score > 0 {
        Some((q_end - HSP_MAX_WINDOW / 2, s_end - HSP_MAX_WINDOW / 2))
    } else {
        None
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:508-513:
#[allow(dead_code)]
fn full_translation_traceback_with_matrix(
    query: &[u8],
    subject: &[u8],
    db_gencode: u8,
    gapped: &[GappedHsp],
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop_final: i32,
    percent_identity: f64,
) -> Result<Vec<(GappedHsp, bool)>> {
    full_translation_traceback_with_matrix_and_events(
        query,
        subject,
        db_gencode,
        gapped,
        matrix,
        gap_open,
        gap_extend,
        x_drop_final,
        percent_identity,
        None,
        None,
        None,
        None,
        None,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583-605:
// This diagnostic observer records the real call order and input offsets.
#[allow(dead_code)]
fn full_translation_traceback_with_matrix_and_events(
    query: &[u8],
    subject: &[u8],
    db_gencode: u8,
    gapped: &[GappedHsp],
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop_final: i32,
    percent_identity: f64,
    seg: Option<&SegParams>,
    test_events: Option<&mut Vec<(bool, usize, usize, usize, usize, usize)>>,
    start_events: Option<&mut Vec<(bool, i32, i32, i32, i32, i32, i32)>>,
    containment_events: Option<&mut Vec<(bool, i8, i32, i32, i32, i32)>>,
    identity_events: Option<&mut Vec<(usize, usize)>>,
) -> Result<Vec<(GappedHsp, bool)>> {
    full_translation_traceback_with_matrix_and_events_with_mask_mode(
        query,
        subject,
        db_gencode,
        gapped,
        matrix,
        gap_open,
        gap_extend,
        x_drop_final,
        percent_identity,
        seg,
        false,
        false,
        test_events,
        start_events,
        containment_events,
        identity_events,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:614-625;
// c++/src/algo/blast/core/blast_traceback.c:380-391,583-596:
// Soft masking keeps both traceback and identity inputs unmasked.
#[allow(dead_code)]
fn full_translation_traceback_with_matrix_and_events_with_mask_mode(
    query: &[u8],
    subject: &[u8],
    db_gencode: u8,
    gapped: &[GappedHsp],
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop_final: i32,
    percent_identity: f64,
    seg: Option<&SegParams>,
    soft_masking: bool,
    mask_lowercase: bool,
    test_events: Option<&mut Vec<(bool, usize, usize, usize, usize, usize)>>,
    start_events: Option<&mut Vec<(bool, i32, i32, i32, i32, i32, i32)>>,
    containment_events: Option<&mut Vec<(bool, i8, i32, i32, i32, i32)>>,
    identity_events: Option<&mut Vec<(usize, usize)>>,
) -> Result<Vec<(GappedHsp, bool)>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:555,844:
    // The Stage C observers pass one subject's preliminary HSPs, which reach the
    // traceback in this score order.
    let mut ordered = gapped.to_vec();
    sort_gapped_score_if_needed(&mut ordered);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583-605,675-693:
    // return the retained HSPs in source order to existing Stage C observers.
    Ok(
        full_translation_traceback_with_matrix_and_events_with_mask_mode_owned(
            query,
            SubjectNt::Letters(subject),
            db_gencode,
            &ordered,
            matrix,
            gap_open,
            gap_extend,
            x_drop_final,
            percent_identity,
            seg,
            soft_masking,
            mask_lowercase,
            test_events,
            start_events,
            containment_events,
            identity_events,
        )?
        .0
        .into_iter()
        .map(|item| (item.hsp, item.retried))
        .collect(),
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-536,583-605,675-693:
// Sorts the HSP list by score (`Blast_HSPListSortByScore`) and purges contained HSPs afterward.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:294,425-433,717-719;
// core/blast_hits.c:1231-1235: the last translated_length returned while
// visiting input HSPs becomes the posttraceback statistical subject length.
pub(super) fn full_translation_traceback_with_matrix_and_events_with_mask_mode_owned(
    query: &[u8],
    subject: SubjectNt<'_>,
    db_gencode: u8,
    gapped: &[GappedHsp],
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop_final: i32,
    percent_identity: f64,
    seg: Option<&SegParams>,
    soft_masking: bool,
    mask_lowercase: bool,
    test_events: Option<&mut Vec<(bool, usize, usize, usize, usize, usize)>>,
    start_events: Option<&mut Vec<(bool, i32, i32, i32, i32, i32, i32)>>,
    containment_events: Option<&mut Vec<(bool, i8, i32, i32, i32, i32)>>,
    identity_events: Option<&mut Vec<(usize, usize)>>,
) -> Result<(Vec<TracebackOwnedHsp>, i32)> {
    x_traceback_in_windows(
        query,
        subject,
        db_gencode,
        gapped,
        matrix,
        gap_open,
        gap_extend,
        x_drop_final,
        percent_identity,
        seg,
        soft_masking,
        mask_lowercase,
        test_events,
        start_events,
        containment_events,
        identity_events,
        X_WHOLE_MARGIN,
        super::score::MatrixScore,
    )
    .map(|(results, stat_length, _)| (results, stat_length))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:429-439,509-512
// (`gap_align->sbp` holds the PSSM: `positionBased` reads its rows.)
/// `full_translation_traceback_with_matrix_and_events_with_mask_mode_owned` of a PSSM query
/// (`-in_pssm`, PSI-TBLASTN): the same traceback with the scores of the PSSM's rows; the
/// identities and positives read the unmasked query letters with `matrix` (NCBI's `sbp->matrix`,
/// the `-matrix`).
#[allow(clippy::too_many_arguments)]
pub(super) fn x_pssm_traceback_owned(
    rows: &crate::algorithm::pssm_search::rows::QueryPssm,
    query: &[u8],
    subject: SubjectNt<'_>,
    db_gencode: u8,
    gapped: &[GappedHsp],
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop_final: i32,
    seg: Option<&SegParams>,
    soft_masking: bool,
    mask_lowercase: bool,
) -> Result<(Vec<TracebackOwnedHsp>, i32)> {
    x_traceback_in_windows(
        query,
        subject,
        db_gencode,
        gapped,
        matrix,
        gap_open,
        gap_extend,
        x_drop_final,
        0.0,
        seg,
        soft_masking,
        mask_lowercase,
        None,
        None,
        None,
        None,
        X_WHOLE_MARGIN,
        super::score::PssmScore(rows),
    )
    .map(|(results, stat_length, _)| (results, stat_length))
}

/// `full_translation_traceback_with_matrix_and_events_with_mask_mode_owned` whose first windows
/// of the whole translation after a fence hold `first_margin` residues (at least
/// `HSP_MAX_WINDOW`) on either side of an HSP, and the number of times the list was traced in
/// such windows (0 when no fence was read); `score` is the source of the scores (the matrix, or a
/// PSSM query's rows).
fn x_traceback_in_windows<S: super::score::TbnScore>(
    query: &[u8],
    subject: SubjectNt<'_>,
    db_gencode: u8,
    gapped: &[GappedHsp],
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_drop_final: i32,
    percent_identity: f64,
    seg: Option<&SegParams>,
    soft_masking: bool,
    mask_lowercase: bool,
    mut test_events: Option<&mut Vec<(bool, usize, usize, usize, usize, usize)>>,
    mut start_events: Option<&mut Vec<(bool, i32, i32, i32, i32, i32, i32)>>,
    mut containment_events: Option<&mut Vec<(bool, i8, i32, i32, i32, i32)>>,
    mut identity_events: Option<&mut Vec<(usize, usize)>>,
    first_margin: usize,
    score: S,
) -> Result<(Vec<TracebackOwnedHsp>, i32, usize)> {
    let code = GeneticCode::try_from_id(db_gencode).map_err(anyhow::Error::msg)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:380-391,583-596:
    // Trace the masked working query, then count identities on its
    // separately retained unmasked sequence.
    let query_frame = if soft_masking {
        encode_protein_query_frame_with_seg(query, None)
    } else {
        encode_tblastn_lookup_query(query, seg, mask_lowercase)
    };
    let query_sequence = &query_frame.aa_seq[1..query_frame.aa_seq.len() - 1];
    let query_nomask = query_frame
        .aa_seq_nomask
        .as_ref()
        .map_or(query_sequence, |bytes| &bytes[1..bytes.len() - 1]);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:358-365
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3269-3284
    // The release traceback keeps the order the preliminary hit list left: by
    // score, or by E-value once that hit list was full and heapified. The order
    // decides which HSP is translated last, i.e. `stat_length` (TN-2).
    let ordered = gapped.to_vec();
    let query_length = i32::try_from(query.len())?;
    let subject_limit = i32::try_from(subject.len() / 3)? + 1;
    let mut target = TargetTranslation::of(subject, &code);
    // After a fence (pass 1) each HSP is traced in a window of the whole translation, `margin`
    // residues on either side (`x_whole_window`).
    let mut margin = first_margin.max(HSP_MAX_WINDOW);
    let mut pass = 0;
    let mut attempts = 0;
    loop {
        // The events before this pass, kept when a window of pass 1 turns out too narrow.
        let marks = [
            test_events.as_deref().map_or(0, Vec::len),
            start_events.as_deref().map_or(0, Vec::len),
            containment_events.as_deref().map_or(0, Vec::len),
            identity_events.as_deref().map_or(0, Vec::len),
        ];
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:294,425-433: initialized from raw subject
        // length, then overwritten by each positive translated_length.
        let mut stat_length = i32::try_from(subject.len())?;
        let mut results = Vec::new();
        let mut retry = false;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:352-360,401-405,
        // 583-605:
        // Initializes the traceback interval tree (`Blast_IntervalTreeInit`).
        let mut tree = BlastIntervalTree::new(0, query_length + 1, 0, subject_limit);
        for hit in &ordered {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:403-405,607-609:
            // Checks interval-tree containment (`BlastIntervalTreeContainsHSP`) before
            // obtaining target translation or attempting traceback.
            let contained = tree.contains_hsp(&traceback_tree_hsp(hit, query_length), 0, 0);
            if let Some(events) = containment_events.as_deref_mut() {
                events.push((
                    contained,
                    hit.frame,
                    hit.q_start,
                    hit.q_end,
                    hit.s_start,
                    hit.s_end,
                ));
            }
            if contained {
                continue;
            }
            let (subject_sequence, translated_length, subject_base) = if pass == 0 {
                target.get(hit.frame, hit.s_start, hit.s_end)?
            } else {
                target.x_whole_window(hit.frame, hit.s_start, hit.s_end, margin)?
            };
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:425-433:
            if translated_length > 0 {
                stat_length = i32::try_from(translated_length)?;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:436-449:
            // Recomputes starting offsets via `BlastGetOffsetsForGappedAlignment`; when
            // none is found, frees the HSP (`Blast_HSPFree`) and skips to the next entry.
            let (q_start, s_start) = if hit.q_gapped_start == 0 && hit.s_gapped_start == 0 {
                let q_offset = usize::try_from(hit.q_start)?;
                let q_end = usize::try_from(hit.q_end)?;
                let s_offset = window_index(hit.s_start, subject_base)
                    .context("HSP subject start before translated window")?;
                let s_end = window_index(hit.s_end, subject_base)
                    .context("HSP subject end before translated window")?;
                let start = score.offsets(
                    query_sequence,
                    subject_sequence,
                    q_offset,
                    q_end,
                    s_offset,
                    s_end,
                    matrix,
                );
                if let Some(events) = start_events.as_deref_mut() {
                    events.push((
                        start.is_some(),
                        hit.q_start,
                        hit.q_end,
                        hit.s_start,
                        hit.s_end,
                        start.map_or(-1, |row| row.0 as i32),
                        start.map_or(-1, |row| (row.1 as isize + subject_base) as i32),
                    ));
                }
                let Some(start) = start else { continue };
                start
            } else {
                (
                    usize::try_from(hit.q_gapped_start)?,
                    window_index(hit.s_gapped_start, subject_base)
                        .context("gapped start before translated window")?,
                )
            };
            // Only a gapped start of a one-hit HSP of negative length with a PSSM lies after the
            // frame's NULLB (its letter `aa_len`, the last byte of the window there: the fence in
            // pass 0, the NULLB in pass 1).
            if S::GAPPED_START_PAST_FRAME {
                x_traceback_start_in_frame(subject.len(), hit.frame, s_start, subject_base)?;
            }
            let mut fence_hit = false;
            let mut scratch = S::traceback_scratch();
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:509-512
            // The window's letters end before its sentinel (`subject_sequence` ends with it),
            // which only the last column of the right extension reads (TN-1).
            let alignment = score.traceback(
                query_sequence,
                subject_sequence,
                subject_sequence.len() - 1,
                q_start,
                s_start,
                matrix,
                gap_open,
                gap_extend,
                x_drop_final,
                &mut scratch,
                &mut fence_hit,
            );
            if fence_hit {
                retry = true;
                break;
            }
            let alignment = alignment.context("NCBI protein traceback returned no alignment")?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:516-520,583-596
            // Array/report consumers run only after the caller's fence retry decision.
            let query_start = usize::try_from(alignment.query_start)?;
            let query_stop = usize::try_from(alignment.query_stop)?;
            let subject_start = usize::try_from(alignment.subject_start)?;
            let subject_stop = usize::try_from(alignment.subject_stop)?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:585-605:
            let align_length: usize = alignment
                .edit_script
                .iter()
                .map(|op| op.num() as usize)
                .sum();
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583-596:
            // Computes identity and alignment length from the adjusted subject, HSP, and edit script.
            let num_ident = protein_identities_from_edit_ops(
                query_nomask,
                subject_sequence,
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583-596: identities use the post-fence valid HSP offsets.
                query_start,
                subject_start,
                &alignment.edit_script,
                matrix,
            );
            if let Some(events) = identity_events.as_deref_mut() {
                events.push((num_ident, align_length));
            }
            // NCBI's `min_hit_length` is always 0 (`core::blast_hits::hsp_test`).
            let delete_hsp = crate::core::blast_hits::hsp_test(
                num_ident as i64,
                align_length as i64,
                percent_identity,
                0,
            );
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:585-605 records updated HSP offsets
            // before Blast_HSPTest and before Blast_HSPAdjustSubjectOffset.
            if let Some(events) = test_events.as_deref_mut() {
                events.push((
                    delete_hsp,
                    align_length,
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583-596: observer sees the post-fence updated HSP offsets.
                    query_start,
                    query_stop,
                    (subject_start as isize + subject_base) as usize,
                    (subject_stop as isize + subject_base) as usize,
                ));
            }
            if delete_hsp {
                continue;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:446-448:
            // The Rust subject slice is rebased; restore its global offset.
            let saved = GappedHsp {
                frame: hit.frame,
                score: alignment.score,
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:583
                q_start: alignment.query_start,
                q_end: alignment.query_stop,
                q_gapped_start: i32::try_from(q_start)?,
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:599
                s_start: alignment.subject_start + i32::try_from(subject_base)?,
                s_end: alignment.subject_stop + i32::try_from(subject_base)?,
                s_gapped_start: i32::try_from(s_start as isize + subject_base)?,
            };
            tree.add_hsp(
                traceback_tree_hsp(&saved, query_length),
                0,
                IndexMethod::QueryAndSubject,
            );
            results.push(TracebackOwnedHsp {
                hsp: saved,
                retried: pass != 0,
                edit_script: alignment.edit_script,
            });
        }
        if !retry {
            let mut results = purge_traceback_common_endpoints(results);
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:672-693:
            let mut final_tree = BlastIntervalTree::new(0, query_length + 1, 0, subject_limit);
            results.retain(|item| {
                let hsp = &item.hsp;
                let tree_hsp = traceback_tree_hsp(hsp, query_length);
                let contained = final_tree.contains_hsp(&tree_hsp, 0, 0);
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:675-693:
                // After score sort, test
                // interval-tree containment again before retaining HSPs.
                if let Some(events) = containment_events.as_deref_mut() {
                    events.push((
                        contained,
                        hsp.frame,
                        hsp.q_start,
                        hsp.q_end,
                        hsp.s_start,
                        hsp.s_end,
                    ));
                }
                if contained {
                    false
                } else {
                    final_tree.add_hsp(tree_hsp, 0, IndexMethod::QueryAndSubject);
                    true
                }
            });
            return Ok((results, stat_length, attempts));
        }
        if pass == 1 {
            // A window of the whole translation was too narrow: the list again with wider
            // windows, its events as if this attempt had not been made. With a margin of a
            // frame's length every window is the whole translation, which has no fence.
            if margin >= subject.len() / 3 {
                bail!("NCBI full-subject TBLASTN traceback reached a fence");
            }
            margin = margin.saturating_mul(8);
            if let Some(events) = test_events.as_deref_mut() {
                events.truncate(marks[0]);
            }
            if let Some(events) = start_events.as_deref_mut() {
                events.truncate(marks[1]);
            }
            if let Some(events) = containment_events.as_deref_mut() {
                events.truncate(marks[2]);
            }
            if let Some(events) = identity_events.as_deref_mut() {
                events.truncate(marks[3]);
            }
        }
        pass = 1;
        attempts += 1;
        target = TargetTranslation::of(subject, &code);
    }
}

#[cfg(test)]
#[path = "search_gapped_tests.rs"]
mod tests;
