//! The BLASTN pairwise report (outfmt 0): the displayed alignment of every HSP of the final
//! hit list, and the per-query statistics that the report prints.

use std::iter::repeat_n;

use anyhow::{ensure, Context, Result};
use bio::io::fasta;

use crate::common::GapEditOp;
use crate::core::blast_encoding::encode_iupac_to_blastna;
use crate::report::pairwise::PairwiseHit;
use crate::stats::{
    compute_karlin_params_ungapped, score_freq_profile_from_probabilities, KarlinParams,
};
use crate::utils::dust::MaskedInterval;

use super::alignment::build_blastna_matrix;
use super::hsp::{BlastnHitList, BlastnHsp};
use super::lookup::reverse_complement;

/// The regions that the report prints in lowercase.
pub(crate) struct DisplayMasks<'a> {
    /// Per query, on its plus strand: DUST and, with `-lcase_masking`, the input lowercase.
    pub query: &'a [Vec<MaskedInterval>],
    /// Per subject, with `-lcase_masking`: the input lowercase.
    pub subject: Option<&'a [Vec<MaskedInterval>]>,
}

/// Every HSP of the final hit list, in the order of the tabular formats (queries, then the
/// subjects and HSPs of each query as ranked), with its displayed rows.
pub(crate) fn pairwise_hits(
    hit_lists: &[Option<BlastnHitList>],
    queries: &[fasta::Record],
    subjects: &[fasta::Record],
    masks: &DisplayMasks<'_>,
) -> Result<Vec<PairwiseHit>> {
    let mut hits = Vec::new();
    for hit_list in hit_lists.iter().flatten() {
        for hsp_list in &hit_list.hsplist_array {
            for hsp in &hsp_list.hsps {
                let query = queries
                    .get(hsp.q_idx as usize)
                    .context("BLASTN HSP of an unknown query")?;
                let subject = subjects
                    .get(hsp.s_idx as usize)
                    .context("BLASTN HSP of an unknown subject")?;
                let (mut query_row, mut subject_row) =
                    displayed_rows(hsp, query.seq(), subject.seq())?;
                let query_masks =
                    shown_query_masks(&masks.query[hsp.q_idx as usize], query.seq().len());
                lowercase_masked(&mut query_row, hsp.q_start, 1, &query_masks);
                if let Some(subject_masks) = masks.subject {
                    let step = if hsp.s_start > hsp.s_end { -1 } else { 1 };
                    lowercase_masked(
                        &mut subject_row,
                        hsp.s_start,
                        step,
                        &subject_masks[hsp.s_idx as usize],
                    );
                }
                let gaps = query_row
                    .iter()
                    .chain(&subject_row)
                    .filter(|&&residue| residue == b'-')
                    .count();
                hits.push(PairwiseHit {
                    hit: hsp.clone().into_hit(),
                    query_seq: Some(String::from_utf8(query_row).context("BLASTN query row")?),
                    subject_seq: Some(
                        String::from_utf8(subject_row).context("BLASTN subject row")?,
                    ),
                    query_frame: None,
                    subject_frame: None,
                    positives: None,
                    gaps: Some(gaps),
                    subject_length: Some(subject.seq().len()),
                    subject_title: subject.desc().map(str::to_string),
                    comp_adjust_method: None,
                    sum_n: None,
                });
            }
        }
    }
    Ok(hits)
}

// The rows of one HSP as NCBI shows them: the query on its plus strand.
//
// NCBI reference: c++/src/objtools/align_format/showalign.cpp:1852-1858
// ```c++
//         if((ds->IsSetStrands()
//             && ds->GetStrands().front()==eNa_strand_minus)
//            && !(ds->IsSetWidths() && ds->GetWidths()[0] == 3)){
//             //show plus strand if master is minus for non-translated case
//             finalDenseg->Reverse();
// ```
// NCBI reference: c++/include/algo/blast/core/gapinfo.h:43-61
// ```c
// typedef enum EGapAlignOpType {
//    eGapAlignDel = 0, /**< Deletion: a gap in query */
//    ...
//    eGapAlignSub = 3, /**< Substitution */
//    ...
//    eGapAlignIns = 6, /**< Insertion: a gap in subject */
// ```
// The HSP offsets and its edit script are those of the query context (the reverse
// complement of the query for a minus-strand HSP) against the plus-strand subject.
// Reversing the alignment of a minus-strand HSP reverse-complements both rows. The
// residues are the input letters in uppercase (NCBI's IUPAC sequence data); masking is
// applied afterwards.
fn displayed_rows(hsp: &BlastnHsp, query: &[u8], subject: &[u8]) -> Result<(Vec<u8>, Vec<u8>)> {
    let query = query.to_ascii_uppercase();
    let context = if hsp.query_frame < 0 {
        reverse_complement(&query)
    } else {
        query
    };
    let subject = subject.to_ascii_uppercase();
    let ungapped = [GapEditOp::Sub(
        (hsp.internal_q_end_0 - hsp.internal_q_offset_0) as u32,
    )];
    let ops = hsp.gap_info.as_deref().unwrap_or(&ungapped);
    let (mut q, mut s) = (hsp.internal_q_offset_0, hsp.internal_s_offset_0);
    let (mut query_row, mut subject_row) = (Vec::new(), Vec::new());
    for op in ops {
        let n = op.num() as usize;
        let (take_query, take_subject) = match op {
            GapEditOp::Sub(_) => (true, true),
            GapEditOp::Del(_) => (false, true),
            GapEditOp::Ins(_) => (true, false),
        };
        if take_query {
            query_row.extend_from_slice(
                context
                    .get(q..q + n)
                    .context("BLASTN edit script past the query")?,
            );
            q += n;
        } else {
            query_row.extend(repeat_n(b'-', n));
        }
        if take_subject {
            subject_row.extend_from_slice(
                subject
                    .get(s..s + n)
                    .context("BLASTN edit script past the subject")?,
            );
            s += n;
        } else {
            subject_row.extend(repeat_n(b'-', n));
        }
    }
    ensure!(
        q == hsp.internal_q_end_0 && s == hsp.internal_s_end_0,
        "BLASTN edit script does not span its HSP"
    );
    if hsp.query_frame < 0 {
        query_row = reverse_complement(&query_row);
        subject_row = reverse_complement(&subject_row);
    }
    Ok((query_row, subject_row))
}

// NCBI reference: c++/src/algo/blast/api/blast_aux.cpp:885-892
// ```c++
//             const TSeqRange kTarget(s_GetQueryRange...
//             ...
//             if (range.NotEmpty() && range != kTarget) {
//                 ...
//                     (new CSeqLocInfo(seqint, CSeqLocInfo::eFrameNotSet));
// ```
// A mask that covers the whole query is not reported, so it is not shown.
fn shown_query_masks(masks: &[MaskedInterval], query_length: usize) -> Vec<MaskedInterval> {
    masks
        .iter()
        .filter(|mask| !(mask.start == 0 && mask.end >= query_length))
        .cloned()
        .collect()
}

// NCBI reference: c++/src/objtools/align_format/showalign.cpp:2495-2521
// ```c++
//         if(id.Which() != CSeq_id::e_not_set){
//             /*only do this for sequence but not for others like middle line,
//               features*/
//             ...
//                 } else if (m_SeqLocChar==eLowerCase){
//                     actualSeq[i-start]=tolower((unsigned char) actualSeq[i-start]);
// ```
// The row's residue positions start at `first` (1-based) and move by `step`; a residue
// whose position lies in a mask is shown in lowercase.
fn lowercase_masked(row: &mut [u8], first: usize, step: isize, masks: &[MaskedInterval]) {
    if masks.is_empty() {
        return;
    }
    let mut position = first as isize;
    for residue in row.iter_mut().filter(|residue| **residue != b'-') {
        let zero_based = (position - 1) as usize;
        if masks
            .iter()
            .any(|mask| mask.start <= zero_based && zero_based < mask.end)
        {
            residue.make_ascii_lowercase();
        }
        position += step;
    }
}

/// The ungapped Karlin-Altschul block of a blastn query, from its composition; `None`
/// when NCBI cannot compute it and marks the query invalid.
///
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:2778-2792
/// ```c
///       Blast_ResFreqString(sbp, rfp, (char*)buffer, query_length);
///       sbp->sfp[context] = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
///       BlastScoreFreqCalc(sbp, sbp->sfp[context], rfp, stdrfp);
///       sbp->kbp_std[context] = kbp = Blast_KarlinBlkNew();
///       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
///       if (loop_status) {
///           contexts[context].is_valid = FALSE;
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:2004-2013
/// ```c
///    for (lp = str, lpmax = lp+length; lp < lpmax; lp++)
///    {
///       ++rcp->comp[(int)(*lp & mask)];
///    }
///
///    /* Don't count ambig. residues. */
///    for (index=0; index<sbp->ambig_occupy; index++)
///    {
///       rcp->comp[sbp->ambiguous_res[index]] = 0;
///    }
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_setup.c:346-350
/// ```c
///     if (program_number == eBlastTypeBlastn ||
///         program_number == eBlastTypeMapping) {
///
///         BLAST_ScoreSetAmbigRes(sbp, 'N');
///         BLAST_ScoreSetAmbigRes(sbp, '-');
/// ```
/// The residues are counted in BLASTNA without `N` and the gap, normalized, and paired
/// with the standard composition (A, C, G, T at 0.25; `Blast_ResFreqStdComp`) through the
/// BLASTNA matrix (`BlastScoreFreqCalc`, blast_stat.c:2151-2205). Both strands of a query
/// give the same block, so the plus strand is used. A query of `N` only has no counted
/// residue, so every score probability is 0 and the calculation fails.
pub(crate) fn query_ungapped_karlin(
    query: &[u8],
    reward: i32,
    penalty: i32,
) -> Option<KarlinParams> {
    const BLASTNA_N: usize = 14;
    const BLASTNA_GAP: usize = 15;
    let mut counts = [0f64; 16];
    for code in encode_iupac_to_blastna(&query.to_ascii_uppercase()) {
        counts[(code & 0x0f) as usize] += 1.0;
    }
    counts[BLASTNA_N] = 0.0;
    counts[BLASTNA_GAP] = 0.0;
    let sum: f64 = counts.iter().sum();
    if sum == 0.0 {
        return None;
    }
    let matrix = build_blastna_matrix(reward, penalty);
    let mut probabilities = Vec::new();
    for (residue, &count) in counts.iter().enumerate() {
        if count == 0.0 {
            continue;
        }
        for base in 0..4 {
            probabilities.push((matrix[residue * 16 + base], count / sum * 0.25));
        }
    }
    let sfp = score_freq_profile_from_probabilities(
        penalty.min(reward),
        reward.max(penalty),
        &probabilities,
    );
    compute_karlin_params_ungapped(&sfp).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(seq: &[u8]) -> fasta::Record {
        fasta::Record::with_attrs("q", None, seq)
    }

    fn hsp(
        frame: i32,
        q: (usize, usize),
        s: (usize, usize),
        ops: Option<Vec<GapEditOp>>,
    ) -> BlastnHsp {
        BlastnHsp {
            identity: 0.0,
            length: 0,
            mismatch: 0,
            gapopen: 0,
            q_start: 0,
            q_end: 0,
            s_start: 0,
            s_end: 0,
            e_value: 0.0,
            bit_score: 0.0,
            num_ident: 0,
            query_frame: frame,
            query_length: 0,
            q_idx: 0,
            s_idx: 0,
            raw_score: 0,
            internal_q_offset_0: q.0,
            internal_q_end_0: q.1,
            internal_s_offset_0: s.0,
            internal_s_end_0: s.1,
            internal_query_context_offset: 0,
            gap_info: ops,
            num_positives: 0,
        }
    }

    #[test]
    fn plus_strand_rows_follow_the_edit_script() {
        let query = record(b"ACGTTACG");
        let subject = record(b"ACGTACCG");
        let ops = vec![GapEditOp::Sub(4), GapEditOp::Ins(1), GapEditOp::Sub(3)];
        let (q, s) = displayed_rows(
            &hsp(1, (0, 8), (0, 7), Some(ops)),
            query.seq(),
            subject.seq(),
        )
        .unwrap();
        assert_eq!(
            (q.as_slice(), s.as_slice()),
            (&b"ACGTTACG"[..], &b"ACGT-ACC"[..])
        );
    }

    // A minus-strand HSP is shown with the query on its plus strand: both rows of the
    // context alignment are reverse-complemented, gaps included.
    #[test]
    fn minus_strand_rows_are_reverse_complemented() {
        let query = record(b"AACCGGTT");
        let subject = record(b"CCGGT");
        // Context (reverse complement of the query) AACCGGTT[1..6] = ACCGG against CCGG-.
        let ops = vec![GapEditOp::Ins(1), GapEditOp::Sub(4)];
        let (q, s) = displayed_rows(
            &hsp(-1, (1, 6), (0, 4), Some(ops)),
            query.seq(),
            subject.seq(),
        )
        .unwrap();
        assert_eq!((q.as_slice(), s.as_slice()), (&b"CCGGT"[..], &b"CCGG-"[..]));
    }

    #[test]
    fn masked_positions_are_lowercase_in_the_row_direction() {
        let mut row = b"AC-GT".to_vec();
        lowercase_masked(&mut row, 10, -1, &[MaskedInterval { start: 7, end: 9 }]);
        assert_eq!(row, b"Ac-gT");
    }

    #[test]
    fn an_all_n_query_has_no_ungapped_block() {
        assert!(query_ungapped_karlin(b"NNNN", 2, -3).is_none());
        let acgt = query_ungapped_karlin(b"ACGTACGTTTGA", 2, -3).unwrap();
        assert!((acgt.lambda - 0.634).abs() < 5e-4);
    }
}
