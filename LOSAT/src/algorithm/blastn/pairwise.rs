//! The BLASTN pairwise report (outfmt 0): the displayed alignment of every HSP of the final
//! hit list, and the per-query statistics that the report prints.

use std::iter::repeat_n;

use crate::blastinput::fasta_reader::FastaRecord;
use anyhow::{ensure, Context, Result};

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
    /// Per query, on its plus strand of the searched letters: DUST and, with
    /// `-lcase_masking`, the input lowercase.
    pub query: &'a [Vec<MaskedInterval>],
    /// Where each query's searched letters lie in its record (`-query_loc`).
    pub query_offsets: &'a crate::blastinput::seq_range::Placements,
    /// Per subject, with `-lcase_masking`: the input lowercase of the subject's record.
    pub subject: Option<&'a [Vec<MaskedInterval>]>,
}

/// Every HSP of the final hit list, in the order of the tabular formats (queries, then the
/// subjects and HSPs of each query as ranked), with its displayed rows.
// Kept out of line: this runs only for outfmt 0 or hit records, and inlining it into the
// shared post-processing makes every BLASTN run compile it in Wasm hosts.
#[inline(never)]
pub(crate) fn pairwise_hits(
    hit_lists: &[Option<BlastnHitList>],
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    masks: &DisplayMasks<'_>,
) -> Result<Vec<PairwiseHit>> {
    let mut hits = Vec::new();
    let query_masks: Vec<Vec<MaskedInterval>> = masks
        .query
        .iter()
        .zip(queries)
        .enumerate()
        .map(|(q_idx, (query_masks, query))| {
            shown_query_masks(
                query_masks,
                query.seq().len(),
                masks.query_offsets.offset(q_idx),
            )
        })
        .collect();
    let subject_masks: Option<Vec<Vec<MaskedInterval>>> = masks
        .subject
        .map(|masks| masks.iter().map(|masks| merged(masks)).collect());
    for hit_list in hit_lists.iter().flatten() {
        for hsp_list in &hit_list.hsplist_array {
            let list_subject_masks = subject_masks.as_ref().and_then(|subject_masks| {
                let s_idx = hsp_list.hsps.first()?.s_idx as usize;
                Some(shown_subject_masks(&subject_masks[s_idx], &hsp_list.hsps))
            });
            for hsp in &hsp_list.hsps {
                let query = queries
                    .get(hsp.q_idx as usize)
                    .context("BLASTN HSP of an unknown query")?;
                let subject = subjects
                    .get(hsp.s_idx as usize)
                    .context("BLASTN HSP of an unknown subject")?;
                let (mut query_row, mut subject_row) =
                    displayed_rows(hsp, query.seq(), subject.seq())?;
                lowercase_masked(
                    &mut query_row,
                    hsp.q_start,
                    1,
                    &query_masks[hsp.q_idx as usize],
                );
                if let Some(subject_masks) = &list_subject_masks {
                    let step = if hsp.s_start > hsp.s_end { -1 } else { 1 };
                    lowercase_masked(&mut subject_row, hsp.s_start, step, subject_masks);
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
                    // The title after its first word (the `bio` description of the
                    // inputs that `bio` read alike); the reports take the subject's title
                    // bytes by subject index.
                    subject_title: subject.title.iter().position(|&byte| byte == b' ').map(
                        |space| String::from_utf8_lossy(&subject.title[space + 1..]).into_owned(),
                    ),
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
// applied afterwards. Only the aligned segments are read: the segment of the minus-strand
// context `[qs, qe)` is the reverse complement of query `[len - qe, len - qs)`.
fn displayed_rows(hsp: &BlastnHsp, query: &[u8], subject: &[u8]) -> Result<(Vec<u8>, Vec<u8>)> {
    let (qs, qe) = (hsp.internal_q_offset_0, hsp.internal_q_end_0);
    let (ss, se) = (hsp.internal_s_offset_0, hsp.internal_s_end_0);
    ensure!(
        qs <= qe && qe <= query.len() && ss <= se && se <= subject.len(),
        "BLASTN HSP outside its sequences"
    );
    let context = if hsp.query_frame < 0 {
        reverse_complement(&query[query.len() - qe..query.len() - qs])
    } else {
        query[qs..qe].to_ascii_uppercase()
    };
    let subject = subject[ss..se].to_ascii_uppercase();
    let ungapped = [GapEditOp::Sub((qe - qs) as u32)];
    let ops = hsp.gap_info.as_deref().unwrap_or(&ungapped);
    // Offsets into the two segments.
    let (mut q, mut s) = (0, 0);
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
        q == context.len() && s == subject.len(),
        "BLASTN edit script does not span its HSP"
    );
    if hsp.query_frame < 0 {
        query_row = reverse_complement(&query_row);
        subject_row = reverse_complement(&subject_row);
    }
    Ok((query_row, subject_row))
}

// NCBI reference: c++/src/algo/blast/api/blast_aux.cpp:878-892
// ```c++
//         const TSeqRange kTarget((*query_interval)->GetFrom(),
//                                 (*query_interval)->GetTo());
//         ...
//             TSeqRange range(Map(kTarget, masked_range));
//             if (range.NotEmpty() && range != kTarget) {
//                 ...
//                 CRef<CSeqLocInfo> seqlocinfo
//                     (new CSeqLocInfo(seqint, CSeqLocInfo::eFrameNotSet));
// ```
// NCBI reference: c++/src/algo/blast/api/blast_aux.cpp:825-842
// ```c++
// template <class Position>
// CRange<Position> Map(const CRange<Position>& target,
//                      const CRange<Position>& range)
// {
//     if (target.Empty()) {
//         throw std::runtime_error("Target range is empty");
//     }
//
//     if (range.Empty() ||
//         (range.GetFrom() > target.GetTo()) ||
//         ((range.GetFrom() + target.GetFrom()) > target.GetTo())) {
//         return target;
//     }
//
//     CRange<Position> retval;
//     retval.SetFrom(max(target.GetFrom() + range.GetFrom(), target.GetFrom()));
//     retval.SetTo(min(target.GetFrom() + range.GetTo(), target.GetTo()));
//     return retval;
// }
// ```
// A mask that covers the whole query (its searched interval, `kTarget`) is not reported,
// so it is not shown. The masks of the search lie in the searched letters (`query_length`
// of them); `Map` moves them by the start of the query's interval (`offset`, 0 without
// `-query_loc`) into record coordinates, as the report's HSPs.
fn shown_query_masks(
    masks: &[MaskedInterval],
    query_length: usize,
    offset: usize,
) -> Vec<MaskedInterval> {
    merged(
        &masks
            .iter()
            .filter(|mask| !(mask.start == 0 && mask.end >= query_length))
            .map(|mask| MaskedInterval::new(mask.start + offset, mask.end + offset))
            .collect::<Vec<_>>(),
    )
}

// NCBI reference: c++/src/algo/blast/api/blast_seqalign.cpp:1588-1600
// ```c++
//         // Union subject sequence ranges
//         vector <TSeqRange> ranges;
//         for (int i=0; i<hsp_list->hspcnt; i++) {
//             const BlastHSP* hsp = hsp_list->hsp_array[i];
//             TSeqRange rg;
//             rg.SetFrom(hsp->subject.offset);
//             rg.SetTo(hsp->subject.end);
//             ranges.push_back(rg);
//         }
//
//         // Extract subject masks
//         TMaskedSubjRegions masks;
//         if (!ranges.empty() && seqinfo_src->GetMasks(kOid, ranges, masks)) {
// ```
// NCBI reference: c++/src/algo/blast/api/seqinfosrc_seqvec.cpp:109-125
// ```c++
// static void
// s_SeqIntervalToSeqLocInfo(CRef<CSeq_interval> interval,
//                           const vector <TSeqRange>& target_ranges,
//                           const CSeqLocInfo::ETranslationFrame frame,
//                           TMaskedSubjRegions& retval)
// {
//     TSeqRange loc(interval->GetFrom(), interval->GetTo());
//
//     for (size_t ir=0; ir< target_ranges.size(); ir++) {
//         if (target_ranges[ir] != TSeqRange::GetEmpty() &&
//            loc.IntersectingWith(target_ranges[ir])) {
//            CRef<CSeqLocInfo> sli(new CSeqLocInfo(interval, frame));
//            retval.push_back(sli);
//            return;
//         }
//     }
// }
// ```
// The subject masks shown with the HSPs of one subject: the lowercase regions of its record
// (record coordinates) that meet the closed range `[offset, end]` of an HSP in the searched
// letters (the end is one past the HSP). With `-subject_loc` the two coordinates differ by
// the interval's start, so a region of the record is kept by where it would lie in the
// interval, and shown where it lies in the record.
fn shown_subject_masks(masks: &[MaskedInterval], hsps: &[BlastnHsp]) -> Vec<MaskedInterval> {
    masks
        .iter()
        .filter(|mask| {
            hsps.iter()
                .any(|hsp| mask.start <= hsp.internal_s_end_0 && hsp.internal_s_offset_0 < mask.end)
        })
        .cloned()
        .collect()
}

/// The masked positions as disjoint intervals in order, for `lowercase_masked`.
fn merged(masks: &[MaskedInterval]) -> Vec<MaskedInterval> {
    let mut sorted = masks.to_vec();
    sorted.sort_by_key(|mask| mask.start);
    let mut merged: Vec<MaskedInterval> = Vec::with_capacity(sorted.len());
    for mask in sorted {
        match merged.last_mut() {
            Some(last) if mask.start <= last.end => last.end = last.end.max(mask.end),
            _ => merged.push(mask),
        }
    }
    merged
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
// whose position lies in a mask is shown in lowercase. `masks` are disjoint and in order
// (`merged`), so the masks that overlap the row's span are found by binary search.
fn lowercase_masked(row: &mut [u8], first: usize, step: isize, masks: &[MaskedInterval]) {
    let residues = row.iter().filter(|&&residue| residue != b'-').count();
    if masks.is_empty() || residues == 0 {
        return;
    }
    let last = first as isize + step * (residues as isize - 1);
    let (low, high) = (
        (first as isize).min(last) - 1,
        (first as isize).max(last) - 1,
    );
    let (low, high) = (low as usize, high as usize);
    let candidates = &masks[masks.partition_point(|mask| mask.end <= low)
        ..masks.partition_point(|mask| mask.start <= high)];
    if candidates.is_empty() {
        return;
    }
    let mut position = first as isize;
    for residue in row.iter_mut().filter(|residue| **residue != b'-') {
        let zero_based = (position - 1) as usize;
        if candidates
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
/// BLASTNA matrix (`BlastScoreFreqCalc`, blast_stat.c:2151-2205). The minus strand of a
/// query (its reverse complement) has the same residues in another order, so its block
/// can differ in the last bits; the search computes one per strand (`scoring.rs`). A
/// query of `N` only has no counted residue, so every score probability is 0 and the
/// calculation fails.
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

    fn record(seq: &[u8]) -> FastaRecord {
        FastaRecord::new("Query_1", b"q", seq)
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
