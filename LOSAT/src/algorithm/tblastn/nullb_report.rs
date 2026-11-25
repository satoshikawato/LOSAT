//! Reporting of TBLASTN `-in_pssm -window_size 0` alignments whose one-hit traceback started on
//! their frame's trailing NULLB sentinel (`search_gapped::x_traceback_start_in_frame`).
//!
//! When a PSSM scores NCBIstdaa 0 (`-`) above `BLAST_SCORE_MIN`, a one-hit ungapped seed near the
//! end of a translated frame can place its gapped alignment start on the frame's trailing NULLB
//! (subject gapped start equal to the frame length). NCBI's `BLAST_GappedAlignmentWithTraceback`
//! runs the left extension only (leaving `found_end` false, so `hsp.q_end` and `hsp.s_end` stop
//! one residue before the gapped start), while the left edit script includes the NULLB column in
//! its final `GapEditOp::Sub(k)`.
//!
//! This module is the single owner (P3) of how NCBI's Seq-align conversion
//! (`s_BlastHSP2SeqAlign`), `-subject` remapping (`s_RemapToSubjectLoc` / `CSeq_loc_Mapper`),
//! and formatters (`CreateDensegFromStdseg`, `CAlnVec`, tabular, pairwise, XML/XML2/JSON,
//! Seq-annot) process such an alignment.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4650-4653
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:114-145,181-268,306-368,1630-1645
//! NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1125-1175,1440-1585,3083-3094
//! NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:788-803
//! NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:893-904

use anyhow::{ensure, Context, Result};

use super::search_gapped::TargetTranslation;
use super::stage_d_linking::LinkedHsp;
use super::stage_d_results::{KappaHspPayload, KappaResultHitList};
use crate::algorithm::common::tabular::COUNTS_BLOCK;
use crate::common::GapEditOp;
use crate::core::composition_adjustment::redo_alignment::EMatrixAdjustRule;
use crate::report::row_letters::QueryLetters;
use crate::report::tabular::{TabularField, TabularSpec};
use crate::sequence::packed_subject::SubjectNt;
use crate::utils::genetic_code::GeneticCode;

/// Report-time state for a TBLASTN `-in_pssm -window_size 0` alignment whose traceback started on
/// its frame's trailing NULLB (`None` for every ordinary HSP).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct NullbReport {
    /// Original `BlastHSP` query span (`q_end - q_start`, before Seq-align / mapper adjustment)
    /// used by `hsp_percent_coverage` and tabular `qcovhsp` when scores stay intact.
    pub blast_hsp_q_span: i32,
    /// Whether `-subject`'s Seq-loc mapper truncated the alignment and invalidated its scores
    /// (`x_InvalidateScores`): outfmt 8/9/12 omit the `score` block, while pairwise, tabular, and
    /// XML/XML2/JSON read `-1` / `-1.0` / `-1.0` (displayed as `0.0` for E-value and `0.0` in
    /// the description table when no other HSP of the subject has a positive bit score).
    pub scores_invalidated: bool,
    /// `-subject`, positive frame, `rem in {1, 2}`: last Std-seg carries
    /// `(subject_fuzz_to_in_searched_coords, query_to_decrement)` where `subject_fuzz_to` is
    /// `subject_len - 1` (caller adds `subject.offset`) and `query_to_decrement` is `1` when
    /// `k >= 2` or `0` when `k == 1`.
    pub last_seg_fuzz: Option<(u32, i32)>,
    /// `-subject`, positive frame, `rem in {1, 2}`, `k >= 2`: `CreateDensegFromStdseg` throws
    /// `Std-seg segment lengths not accurate!` (`Seq_align.cpp:801`).
    pub std_seg_lengths_inaccurate: bool,
    /// Positive frame, `rem in {1, 2}`, and (`-db` or `-subject` with `k == 1`):
    /// `CAlnVec::TranslateNAToAA` throws `NA size expected to be divisible by 3`
    /// (`alnvec.cpp:898`).
    pub translate_na_failure: bool,
    /// Positive frame, `-db`, `rem == 0`: the final column of the last `Sub(k)` reads past the
    /// subject's end and is padded with `'-'` by `CAlnVec`.
    pub trailing_nullb_dash: bool,
    /// Negative frame, `-db`, `rem in {0, 1}`: the final `Sub(k)` has a negative start (`-3` or
    /// `-2`), written as unsigned `TSeqPos` in Std-seg and rendered as `k` trailing spaces in the
    /// subject display row (for `rem == 2`, `s_from == -1 == GAP_VALUE` turns the segment directly
    /// into `GapEditOp::Ins(k)` in `edit_script`).
    pub wrapped_minus_tail: usize,
}

/// Whether `linked` / `payload` is a NULLB-start alignment: its edit script consumes more query
/// letters than the HSP query range (specifically two more, since `found_end` stays false and
/// stops `q_end` one residue before `q_gs`, while the left extension covers `q_gs`).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4650-4653
#[inline]
pub(super) fn is_nullb_start(linked: &LinkedHsp, payload: &KappaHspPayload) -> bool {
    let query_letters: i64 = payload
        .edit_script
        .iter()
        .map(|op| match *op {
            GapEditOp::Sub(n) | GapEditOp::Ins(n) => i64::from(n),
            GapEditOp::Del(_) => 0,
        })
        .sum();
    let range = i64::from(linked.hsp.q_end) - i64::from(linked.hsp.q_start);
    query_letters != range
}

/// Whether any HSP in `hitlists` is a NULLB-start alignment (before or after report preparation).
pub(super) fn hitlists_have_nullb(hitlists: &[KappaResultHitList]) -> bool {
    hitlists.iter().any(|hitlist| {
        hitlist.lists().iter().any(|list| {
            list.hsps
                .hsps
                .iter()
                .zip(&list.payloads)
                .any(|(linked, payload)| payload.nullb.is_some() || is_nullb_start(linked, payload))
        })
    })
}

/// Transforms a single NULLB-start HSP in place to reflect NCBI's Seq-align construction (and,
/// when `remaps` is true for `-subject`, `CSeq_loc_Mapper` truncation / dropping) once, inside
/// `KappaResultHitList::put_in_report_order`. Returns `true` if the HSP was a NULLB-start HSP.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:114-145,306-368,1630-1645
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1125-1175,1440-1585,3083-3094
pub(super) fn prepare_nullb_hsp(
    linked: &mut LinkedHsp,
    payload: &mut KappaHspPayload,
    subject_len: usize,
    remaps: bool,
) -> bool {
    if !is_nullb_start(linked, payload) {
        return false;
    }
    let Some(GapEditOp::Sub(k_u32)) = payload.edit_script.last().copied() else {
        return false;
    };
    let k = k_u32 as usize;
    let k_i32 = k_u32 as i32;
    let script_q: i32 = payload
        .edit_script
        .iter()
        .map(|op| match *op {
            GapEditOp::Sub(n) | GapEditOp::Ins(n) => n as i32,
            GapEditOp::Del(_) => 0,
        })
        .sum();
    let script_s: i32 = payload
        .edit_script
        .iter()
        .map(|op| match *op {
            GapEditOp::Sub(n) | GapEditOp::Del(n) => n as i32,
            GapEditOp::Ins(_) => 0,
        })
        .sum();

    let hsp = &mut linked.hsp;
    let blast_hsp_q_span = hsp.q_end - hsp.q_start;
    let frame = i32::from(hsp.frame);
    let f = frame.abs();
    let len_i32 = subject_len as i32;
    let rem = (len_i32 - (f - 1)).rem_euclid(3);

    let mut report = NullbReport {
        blast_hsp_q_span,
        scores_invalidated: false,
        last_seg_fuzz: None,
        std_seg_lengths_inaccurate: false,
        translate_na_failure: false,
        trailing_nullb_dash: false,
        wrapped_minus_tail: 0,
    };

    if frame > 0 {
        if !remaps {
            // Positive frame, -db: raw Std-seg extends 1 codon past aa_len (subject end L + 2 - rem).
            hsp.q_end = hsp.q_start + script_q;
            hsp.s_end = hsp.s_start + script_s;
            if rem == 0 {
                report.trailing_nullb_dash = true;
            } else {
                report.translate_na_failure = true;
            }
        } else {
            // Positive frame, -subject: CSeq_loc_Mapper truncates to [0, L - 1] and invalidates scores.
            // NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:3083-3094
            report.scores_invalidated = true;
            hsp.score = -1;
            payload.bit_score = -1.0;
            linked.evalue = -1.0;
            linked.num = 0;
            payload.num_ident = 0;
            payload.num_positives = 0;
            payload.matrix_adjust_rule = EMatrixAdjustRule::DontAdjustMatrix;

            if rem == 0 {
                // Truncated tail has 3*(k - 1) bases (exact codons, no fuzz).
                if k >= 2 {
                    *payload.edit_script.last_mut().expect("non-empty script") =
                        GapEditOp::Sub(k_u32 - 1);
                } else {
                    payload.edit_script.pop();
                }
                payload.align_length = payload.align_length.saturating_sub(1);
                hsp.q_end = hsp.q_start + script_q - 1;
                hsp.s_end = hsp.s_start + script_s - 1;
            } else {
                // Truncated tail leaves rem in {1, 2} bases -> fuzz-to alt { L - 1 }.
                // NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1538-1565
                let query_to_sub = if k >= 2 { 1 } else { 0 };
                report.last_seg_fuzz = Some((subject_len.saturating_sub(1) as u32, query_to_sub));
                if k >= 2 {
                    report.std_seg_lengths_inaccurate = true;
                    hsp.q_end = hsp.q_start + script_q - 1;
                    hsp.s_end = hsp.s_start + script_s;
                } else {
                    report.translate_na_failure = true;
                    hsp.q_end = hsp.q_start + script_q;
                    hsp.s_end = hsp.s_start + script_s;
                }
            }
        }
    } else if !remaps {
        // Negative frame, -db: s_from = rem - 3 in {-3, -2, -1}.
        // When rem == 2, s_from == -1 == GAP_VALUE, so s_CreateStdSegs writes an empty subject loc
        // (equivalent to GapEditOp::Ins(k)). When rem in {0, 1}, s_from is -3 or -2 (written as
        // unsigned TSeqPos in Std-seg, treated as end-gap spaces by CAlnVec and GetAlignLengths,
        // while XML Hsp_gaps excludes it because starts[seg] != -1).
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:344-355
        hsp.q_end = hsp.q_start + script_q;
        hsp.s_end = hsp.s_start + script_s - k_i32;
        payload.gap_letters += k;
        payload.gap_opens += 1;
        if rem == 2 {
            *payload.edit_script.last_mut().expect("non-empty script") = GapEditOp::Ins(k_u32);
        } else {
            report.wrapped_minus_tail = k;
        }
    } else {
        // Negative frame, -subject: CSeq_loc_Mapper drops the invalid final Std-seg without
        // invalidating scores (scores stay intact).
        // NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1568-1580
        payload.edit_script.pop();
        payload.align_length = payload.align_length.saturating_sub(k);
        hsp.q_end = hsp.q_start + script_q - k_i32;
        hsp.s_end = hsp.s_start + script_s - k_i32;
    }

    payload.nullb = Some(report);
    true
}

/// Whether NCBI's `CBlastTabularInfo::SetFields` enters the `CreateDensegFromStdseg` block for
/// `spec` and the HSP's Seq-align `num_ident` score (`engine_ident`).
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-908
pub(super) fn tabular_enters_denseg(spec: &TabularSpec, engine_ident: usize) -> bool {
    spec.fields.iter().any(|field| COUNTS_BLOCK.contains(field))
        || (spec.has(TabularField::NumIdentical) && engine_ident > 0)
}

/// Whether NCBI's `CBlastTabularInfo::SetFields` (for TBLASTN, where `kNoFetchSequence` is false)
/// enters the `GetWholeAlnSeqString` block that translates subject nucleotides via
/// `CAlnVec::TranslateNAToAA`.
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:973-986
pub(super) fn tabular_enters_translate(spec: &TabularSpec, engine_ident: usize) -> bool {
    tabular_enters_denseg(spec, engine_ident)
        && (spec.needs().rows
            || spec.has(TabularField::NumIdentical)
            || spec.has(TabularField::Mismatches)
            || spec.has(TabularField::PercentIdentical))
}

/// Renders `(query_string, subject_string)` for a NULLB-start alignment whose final `Sub` segment
/// has either `trailing_nullb_dash` (positive frame `-db` `rem == 0`: last column padded with
/// `'-'`) or `wrapped_minus_tail > 0` (negative frame `-db` `rem in {0, 1}`: last `k` columns
/// rendered as end-gap spaces `' '`).
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:279-297,893-923
#[allow(clippy::too_many_arguments)]
pub(super) fn aligned_sequences_nullb(
    query: QueryLetters<'_>,
    subject: SubjectNt<'_>,
    window_buffer: &mut Vec<u8>,
    code: &GeneticCode,
    target: &mut TargetTranslation<'_>,
    frame: i8,
    q_start: i32,
    s_start: i32,
    script: &[GapEditOp],
    nullb: NullbReport,
) -> Result<(String, String)> {
    ensure!(!script.is_empty(), "TBLASTN report HSP has no edit script");
    let unbacked_tail = if nullb.trailing_nullb_dash {
        1usize
    } else {
        nullb.wrapped_minus_tail
    };
    let tail_char = if nullb.trailing_nullb_dash { '-' } else { ' ' };
    let total_subject_res: usize = script
        .iter()
        .filter(|op| !matches!(op, GapEditOp::Ins(_)))
        .map(|op| op.num() as usize)
        .sum();
    let backed_subject_res = total_subject_res
        .checked_sub(unbacked_tail)
        .context("NULLB tail exceeds subject residues")?;
    let backed_s_end = s_start + backed_subject_res as i32;
    let (translated, _, base) = if backed_subject_res > 0 {
        target.get(frame, s_start, backed_s_end)?
    } else {
        (&[][..], 0, s_start as isize)
    };
    let mut q = usize::try_from(q_start)?;
    let mut s = if backed_subject_res > 0 {
        super::search_gapped::window_index(s_start, base)
            .context("TBLASTN report translation window starts after HSP")?
    } else {
        0
    };
    let (nucleotides, first) = match subject {
        SubjectNt::Letters(letters) => (letters, 0),
        SubjectNt::Packed(_) if backed_subject_res == 0 => (&[][..], 0),
        SubjectNt::Packed(_) => {
            let first = usize::try_from(s_start)?;
            let frame_shift = usize::from(frame.unsigned_abs());
            let length = 3 * backed_subject_res + frame_shift - 1;
            let start = if frame > 0 {
                Some(3 * first)
            } else {
                (subject.len() + 1).checked_sub(3 * (first + backed_subject_res) + frame_shift)
            }
            .filter(|start| start + length <= subject.len())
            .context("TBLASTN report HSP outside the subject")?;
            (subject.fetch(start, start + length, window_buffer), first)
        }
    };
    let display = crate::algorithm::tblastx::report::DisplayCode::new(code);
    let displayed_subject = |s: usize| -> u8 {
        display.residue(nucleotides, frame, (base + s as isize) as usize - first)
    };
    const NCBISTDAA_TO_AA: &[u8; 28] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
    let columns = script.iter().map(|op| op.num() as usize).sum();
    let mut query_string = String::with_capacity(columns);
    let mut subject_string = String::with_capacity(columns);
    let mut s_emitted = 0usize;
    for &op in script {
        for _ in 0..op.num() {
            match op {
                GapEditOp::Sub(_) => {
                    let qa = query
                        .get(q)
                        .context("TBLASTN report query offset overflow")?;
                    query_string.push(char::from(qa));
                    if s_emitted < backed_subject_res {
                        let &sa = translated
                            .get(s)
                            .context("TBLASTN report subject offset overflow")?;
                        ensure!(
                            NCBISTDAA_TO_AA.get(sa as usize).is_some(),
                            "invalid TBLASTN subject amino acid"
                        );
                        subject_string.push(char::from(displayed_subject(s)));
                        s += 1;
                    } else {
                        subject_string.push(tail_char);
                    }
                    q += 1;
                    s_emitted += 1;
                }
                GapEditOp::Del(_) => {
                    let &sa = translated
                        .get(s)
                        .context("TBLASTN report subject offset overflow")?;
                    query_string.push('-');
                    ensure!(
                        NCBISTDAA_TO_AA.get(sa as usize).is_some(),
                        "invalid TBLASTN subject amino acid"
                    );
                    subject_string.push(char::from(displayed_subject(s)));
                    s += 1;
                    s_emitted += 1;
                }
                GapEditOp::Ins(_) => {
                    let qa = query
                        .get(q)
                        .context("TBLASTN report query offset overflow")?;
                    query_string.push(char::from(qa));
                    subject_string.push('-');
                    q += 1;
                }
            }
        }
    }
    Ok((query_string, subject_string))
}

#[cfg(test)]
#[path = "nullb_report_tests.rs"]
mod tests;
