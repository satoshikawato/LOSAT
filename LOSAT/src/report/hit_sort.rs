//! The order of one query's hits and HSPs in the pairwise report (outfmt 0) with `-sorthits` and
//! `-sorthsps` (NCBI's `CAlignFormatUtil::SortSeqalignForSortableFormat`), shared by BLASTN,
//! BLASTP, BLASTX, TBLASTN and TBLASTX, and the values that the description table's columns of
//! `-sorthits` read.
//!
//! NCBI sorts `std::list`s, whose `sort` is stable (the C++ standard: elements that compare equal
//! keep their order). Each sort here is Rust's stable `sort_by` with a total order that gives the
//! order of NCBI's comparator, plus the explicit steps for the comparators that are not strict
//! weak orders or that reorder the compared hits' HSPs (`sort_query_hits`). The steps were taken
//! from the outputs of `std::list::sort` built with g++ 13 for comparators of the same shapes, on
//! lists of 0 to 300 elements with many ties (the order of equal totals with `>=`; which hits'
//! HSPs end sorted).

use std::cmp::Ordering;

use super::pairwise::PairwiseHit;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2938-2946,2979-2988
// ```c++
//     if(m_OutputFormat <= eFlatQueryAnchoredNoIdentities) {
//     	...
//     	if (args[kArgLineLength]) {
//     	    m_LineLength = args[kArgLineLength].AsInteger();
//     	}
//     	if(args.Exist(kArgSortHits) && args[kArgSortHits])
//     	{
//     		m_HitsSortOption = args[kArgSortHits].AsInteger();
//     	}
//     ...
//     if(args.Exist(kArgSortHSPs) && args[kArgSortHSPs])
//     {
//         int hspsSortOption = args[kArgSortHSPs].AsInteger();
//         if(m_OutputFormat == ePairwise) {
//        	    m_HspsSortOption = hspsSortOption;
//         }
// ```
/// `-sorthits` and `-sorthsps` as the pairwise report reads them (NCBI's `m_HitsSortOption` and
/// `m_HspsSortOption`; `None` is NCBI's -1, the option not given). psiblast and rpsblast keep the
/// default: their applications do not hand the options to the formatter
/// (`psiblast_app.cpp:502`, `rpsblast_app.cpp:187` set the line length only).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ReportSort {
    /// `-sorthits` (0 to 4).
    pub hits: Option<u8>,
    /// `-sorthsps` (0 to 4).
    pub hsps: Option<u8>,
}

impl ReportSort {
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1522-1530
    // ```c++
    //     //invoke sorting only for m_HitsSortOption > CAlignFormatUtil::eEvalue or m_HspsSortOption > CAlignFormatUtil::eHspEvalue
    //     if(m_HitsSortOption > 0 || m_HspsSortOption > 0) {
    //         aln_set = CBlastFormatUtil::SortSeqalignForSortableFormat(
    // ```
    /// Whether the report sorts each query's hits (`m_HitsSortOption > 0 || m_HspsSortOption >
    /// 0`): a value of 0 (by E-value, the order of the search's results) sorts nothing.
    pub fn sorts(&self) -> bool {
        self.hits.is_some_and(|hits| hits > 0) || self.hsps.is_some_and(|hsps| hsps > 0)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:518-522
    // ```c++
    //     if(m_HitsSortOption >= 0) {
    //         flags |= CShowBlastDefline::eShowPercentIdent;
    //         flags |= CShowBlastDefline::eShowTotalScore;
    //         flags |= CShowBlastDefline::eShowQueryCoverage;
    //     }
    // ```
    /// Whether the description table shows the total score, query cover and identity columns:
    /// `-sorthits` is given, at any value (0 too).
    pub fn shows_columns(&self) -> bool {
        self.hits.is_some()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:518-522,1522-1530
    // ```c++
    //     if(m_HitsSortOption >= 0) {
    //         flags |= CShowBlastDefline::eShowPercentIdent;
    //     ...
    //     if(m_HitsSortOption > 0 || m_HspsSortOption > 0) {
    //         aln_set = CBlastFormatUtil::SortSeqalignForSortableFormat(
    // ```
    /// Whether the report reads the HSPs' `HspScores` (the Seq-aligns' identities,
    /// `PairwiseHit::seqalign_num_ident`): it sorts (`sorts`) or shows the columns
    /// (`shows_columns`). Without either the report reads none of them, and TBLASTX does not
    /// copy its HSPs' identities for it (`tblastx::report::seqalign_identities`).
    pub(crate) fn reads_hsp_scores(&self) -> bool {
        self.sorts() || self.shows_columns()
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:813-848
// ```c++
//     bool hasScore = false;
//     score = -1;
//     bits = -1;
//     evalue = -1;
//     sum_n = -1;
//     num_ident = -1;
//     ...
//     //look for scores at seqalign level first
//     hasScore = s_GetBlastScore(aln.GetScore(), score, bits, evalue,
//                                sum_n, num_ident, use_this_gi, comp_adj_method);
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:252-264
// ```c++
// TSeqPos CSeq_align::GetSeqStart(TDim row) const
// {
//     switch (GetSegs().Which()) {
//     case C_Segs::e_Denseg:
//         return GetSegs().GetDenseg().GetSeqStart(row);
//     ...
//     case C_Segs::e_Dendiag:
//     case C_Segs::e_Std:
//     case C_Segs::e_Disc:
//         return GetSeqRange(row).GetFrom();
// ```
/// The values of one HSP that the sorts and the description table's columns read: its scores
/// (`GetAlnScores`: the bit score, the E-value, `e_value` or the linked set's `sum_e`, which the
/// report shows, and the Seq-align's `num_ident`, the engine's identities, which every HSP's
/// Seq-align carries, `blast_seqalign.cpp:1199-1203`, `PairwiseHit::seqalign_num_ident`, not the
/// identities of the displayed rows), its alignment length in columns of letters (TBLASTX's in
/// amino acids) and the 0-based ranges of the query (`from..=to`) and the subject (`from`) that
/// `GetSeqRange` gives (BLASTX's and TBLASTX's query and TBLASTN's subject in nucleotides).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct HspScores {
    pub(crate) bits: f64,
    pub(crate) evalue: f64,
    pub(crate) num_ident: i64,
    pub(crate) length: i64,
    pub(crate) query_from: u64,
    pub(crate) query_to: u64,
    pub(crate) subject_from: u64,
}

impl HspScores {
    /// The values of the HSP of `hit` (its 1-based coordinates as the report shows them, from
    /// <= to after `GetSeqRange`'s swap of a minus strand).
    pub(crate) fn of(hit: &PairwiseHit) -> Self {
        Self::of_hit(&hit.hit, hit.seqalign_num_ident)
    }

    /// The values of the HSP `h` of a report's hit list whose Seq-align carries the identities
    /// `seqalign_num_ident` (`PairwiseHit::seqalign_num_ident`): [`Self::of`] of the report's row
    /// of it, which the models of the reads read before the rows are made
    /// (`x_pairwise_par::QueryGroup::read_order`).
    pub(crate) fn of_hit(h: &crate::common::Hit, seqalign_num_ident: usize) -> Self {
        Self::of_values(
            h.bit_score,
            h.e_value,
            seqalign_num_ident,
            h.length,
            (h.q_start, h.q_end),
            (h.s_start, h.s_end),
        )
    }

    /// The values of an HSP of the report whose row has the bit score `bits`, the E-value
    /// `e_value`, the Seq-align identities `seqalign_num_ident`, `length` columns and the 1-based
    /// coordinates `query` and `subject` (start, end) as the row's `Hit` holds them
    /// ([`Self::of_hit`]).
    pub(crate) fn of_values(
        bits: f64,
        e_value: f64,
        seqalign_num_ident: usize,
        length: usize,
        query: (usize, usize),
        subject: (usize, usize),
    ) -> Self {
        let zero_based = |position: usize| (position as u64).saturating_sub(1);
        Self {
            bits,
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1185-1190
            // ```c++
            //     // Set the E-Value
            //     double evalue = (hsp->evalue < SMALLEST_EVALUE) ? 0.0 : hsp->evalue;
            //     if (evalue >= 0.0) {
            //         string score_type = (hsp->num <= 1) ? "e_value" : "sum_e";
            //         scores.push_back(s_MakeScore(score_type, evalue, 0, false));
            //     }
            // ```
            // The comparators read the Seq-align's E-value (`GetAlnScores`), 0 below 1e-180:
            // two HSPs below it tie.
            evalue: crate::report::serial::seqalign::seqalign_evalue(e_value),
            num_ident: seqalign_num_ident as i64,
            length: length as i64,
            query_from: zero_based(query.0.min(query.1)),
            query_to: zero_based(query.0.max(query.1)),
            subject_from: zero_based(subject.0.min(subject.1)),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1462-1489
    // ```c++
    // int CAlignFormatUtil::GetAlignmentLength(const CSeq_align& aln, bool do_translation)
    // {
    //     ...
    //     if (aln.GetSegs().IsStd()) {
    //         CRef<CSeq_align> denseg_aln = aln.CreateDensegFromStdseg();
    //         // When both query and subject are translated, i.e. tblastx, convert
    //         // to a special type of Dense-seg.
    //         if (do_translation) {
    //             final_aln = denseg_aln->CreateTranslatedDensegFromNADenseg();
    //         } else {
    //             final_aln = denseg_aln;
    //
    //         }
    //     } else if (aln.GetSegs().IsDendiag()) {
    //         final_aln = CreateDensegFromDendiag(aln);
    //     }
    //     ...
    //     CAlnMap alnmap(ds);
    //     return alnmap.GetAlnStop() + 1;
    // ```
    /// `GetAlignmentLength(aln, translation)`: the HSP's columns, except TBLASTX's (`tblastx`),
    /// whose Std-seg of two nucleotide rows becomes a Dense-seg of nucleotide lengths, three per
    /// amino acid, unless `translation` (`CreateTranslatedDensegFromNADenseg`, the lengths in
    /// amino acids).
    pub(crate) fn alignment_length(&self, tblastx: bool, translation: bool) -> i64 {
        if tblastx && !translation {
            3 * self.length
        } else {
            self.length
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1383-1400
// ```c++
// bool CAlignFormatUtil::
// SortHspByScoreDescending(const CRef<CSeq_align>& info1,
//                          const CRef<CSeq_align>& info2)
// {
//     ...
//     GetAlnScores(*info1, score1,  bits1, evalue1, sum_n1, num_ident1, use_this_gi1);
//     GetAlnScores(*info2, score2,  bits2, evalue2, sum_n2, num_ident2, use_this_gi2);
//     return bits1 > bits2;
//
// }
// ```
/// `SortHspByScoreDescending` as an order: `info1` before `info2` when its bit score is higher,
/// equal (kept in their order by the stable sort) when the scores are equal.
pub(crate) fn hsp_by_score_descending(info1: &HspScores, info2: &HspScores) -> Ordering {
    cmp_f64(info2.bits, info1.bits)
}

/// `a` against `b` as C++ compares two doubles with `<`, `>` and `==` (-0.0 equals 0.0), with a
/// NaN above every number and equal to another NaN so that the order is total, as `sort_by`
/// needs. The values that the sorts compare (bit scores, E-values, their sums, identity ratios and
/// percentages) are never NaN.
fn cmp_f64(a: f64, b: f64) -> Ordering {
    a.partial_cmp(&b)
        .unwrap_or_else(|| a.is_nan().cmp(&b.is_nan()))
}

/// Sorts `order` with Rust's stable sort by the order that NCBI's comparator `less` gives: `x`
/// before `y` when `less(x, y)`, after it when `less(y, x)`, else equal (kept in their order, as
/// the stable `std::list::sort` keeps them). Called only on lists of which `less` is a strict weak
/// order, the condition under which the C++ standard defines `list::sort`'s result and under which
/// this is a total order of the list's elements.
fn sort_by_comparator(order: &mut [usize], less: impl Fn(usize, usize) -> bool) {
    order.sort_by(|&x, &y| {
        if less(x, y) {
            Ordering::Less
        } else if less(y, x) {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    });
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1403-1431
// ```c++
//     start1 = min(info1->GetSeqStart(0), info1->GetSeqStop(0));
//     start2 = min(info2->GetSeqStart(0), info2->GetSeqStop(0)) ;
//
//     if (start1 == start2) {
//         //same start then arrange by bits score
//         ...
//         GetAlnScores(*info1, score1,  bits1, evalue1, sum_n1, num_ident1, use_this_gi1);
//         GetAlnScores(*info2, score2,  bits2, evalue2, sum_n2, num_ident2, use_this_gi2);
//         return evalue1 < evalue2;
//
//     } else {
//
//         return start1 < start2;
//     }
// ```
/// `SortHspByMasterStartAscending` as an order: the lower query start first, equal starts by the
/// lower E-value, equal starts and E-values equal.
pub(crate) fn hsp_by_query_start_ascending(info1: &HspScores, info2: &HspScores) -> Ordering {
    info1
        .query_from
        .cmp(&info2.query_from)
        .then_with(|| cmp_f64(info1.evalue, info2.evalue))
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1433-1458
// ```c++
//     start1 = min(info1->GetSeqStart(1), info1->GetSeqStop(1));
//     start2 = min(info2->GetSeqStart(1), info2->GetSeqStop(1)) ;
//
//     if (start1 == start2) {
//         ...
//         return evalue1 < evalue2;
//
//     } else {
//
//         return start1 < start2;
//     }
// ```
/// `SortHspBySubjectStartAscending` as an order: the lower subject start first, equal starts by
/// the lower E-value, equal starts and E-values equal.
pub(crate) fn hsp_by_subject_start_ascending(info1: &HspScores, info2: &HspScores) -> Ordering {
    info1
        .subject_from
        .cmp(&info2.subject_from)
        .then_with(|| cmp_f64(info1.evalue, info2.evalue))
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1261-1296
// ```c++
//     int length1 = GetAlignmentLength(*info1, kTranslation);
//     int length2 = GetAlignmentLength(*info2, kTranslation);
//     bool retval = false;
//
//
//     if(length1 > 0 && length2 > 0 && num_ident1 > 0 &&num_ident2 > 0 ) {
//         if (((double)num_ident1)/length1 == ((double)num_ident2)/length2) {
//
//             retval = evalue1 < evalue2;
//
//         } else {
//             retval = ((double)num_ident1)/length1 >= ((double)num_ident2)/length2;
//
//         }
//     } else {
//         retval = evalue1 < evalue2;
//     }
//     return retval;
// ```
/// `SortHspByPercentIdentityDescending` with NCBI's `kTranslation` (`translation`; the lengths
/// of `HspScores::alignment_length`).
pub(crate) fn hsp_by_percent_identity_descending(
    info1: &HspScores,
    info2: &HspScores,
    tblastx: bool,
    translation: bool,
) -> bool {
    let length1 = info1.alignment_length(tblastx, translation);
    let length2 = info2.alignment_length(tblastx, translation);
    if length1 > 0 && length2 > 0 && info1.num_ident > 0 && info2.num_ident > 0 {
        let ratio1 = info1.num_ident as f64 / length1 as f64;
        let ratio2 = info2.num_ident as f64 / length2 as f64;
        if ratio1 == ratio2 {
            info1.evalue < info2.evalue
        } else {
            ratio1 >= ratio2
        }
    } else {
        info1.evalue < info2.evalue
    }
}

/// `list.sort(SortHspByPercentIdentityDescending)` of the HSPs `order` of `hsps` (with NCBI's
/// `kTranslation`, `translation`).
///
/// Between two HSPs with identities the comparator is the higher ratio first, equal ratios by the
/// lower E-value (`>=` decides only unequal ratios, where it is `>`); an HSP without identities
/// compares with every HSP by E-value. That is a strict weak order of the list, and NCBI's order is
/// the stable one, unless the list holds an HSP without identities (`order_holds`): two HSPs with
/// identities whose ratio order goes against their E-values (the higher ratio, the larger or equal
/// E-value) and an HSP without identities whose E-value lies between theirs (ends included) are
/// ordered in a cycle (or tie with an HSP that they do not tie with). On such a list NCBI's order
/// is the one that its C++ library's merge procedure gives, which no order of the HSPs reproduces;
/// LOSAT puts the HSPs with identities first, by ratio and E-value as NCBI's comparator orders
/// them, then those without, by E-value (each group stable).
fn sort_hsps_by_percent_identity(
    hsps: &[HspScores],
    order: &mut [usize],
    tblastx: bool,
    translation: bool,
) {
    let identities = |k: usize| {
        let hsp = &hsps[k];
        let length = hsp.alignment_length(tblastx, translation);
        (hsp.num_ident > 0 && length > 0).then(|| hsp.num_ident as f64 / length as f64)
    };
    if percent_identity_order_holds(hsps, order, identities) {
        sort_by_comparator(order, |x, y| {
            hsp_by_percent_identity_descending(&hsps[x], &hsps[y], tblastx, translation)
        });
    } else {
        order.sort_by(|&x, &y| match (identities(x), identities(y)) {
            (Some(ratio1), Some(ratio2)) => {
                cmp_f64(ratio2, ratio1).then_with(|| cmp_f64(hsps[x].evalue, hsps[y].evalue))
            }
            (Some(_), None) => Ordering::Less,
            (None, Some(_)) => Ordering::Greater,
            (None, None) => cmp_f64(hsps[x].evalue, hsps[y].evalue),
        });
    }
}

/// Whether `SortHspByPercentIdentityDescending` is a strict weak order of the HSPs `order` of
/// `hsps`, whose identity ratios `identities` gives (`None` without identities): no E-value is
/// NaN, and no HSP without identities has an E-value from `e_low` to `e_high` (ends included)
/// where an HSP with identities and E-value `e_high` has a higher ratio than one with E-value
/// `e_low` (the largest ratio among the HSPs with identities whose E-value is at least the HSP's
/// is above the smallest ratio among those whose E-value is at most the HSP's).
fn percent_identity_order_holds(
    hsps: &[HspScores],
    order: &[usize],
    identities: impl Fn(usize) -> Option<f64>,
) -> bool {
    if order.iter().any(|&k| hsps[k].evalue.is_nan()) {
        return false;
    }
    let with: Vec<(f64, f64)> = order
        .iter()
        .filter_map(|&k| identities(k).map(|ratio| (ratio, hsps[k].evalue)))
        .collect();
    if with.len() < 2 {
        return true;
    }
    order
        .iter()
        .filter(|&&k| identities(k).is_none())
        .all(|&k| {
            let evalue = hsps[k].evalue;
            let mut highest_at_or_above = f64::NEG_INFINITY;
            let mut lowest_at_or_below = f64::INFINITY;
            for &(ratio, other) in &with {
                if other >= evalue {
                    highest_at_or_above = highest_at_or_above.max(ratio);
                }
                if other <= evalue {
                    lowest_at_or_below = lowest_at_or_below.min(ratio);
                }
            }
            highest_at_or_above <= lowest_at_or_below
        })
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3163-3171
// ```c++
// double CAlignFormatUtil::GetPercentIdentity(int numerator, int denominator)
// {
//      if (numerator == denominator)
//         return 100;
//      else {
//        double retval =100*(double)numerator/(double)denominator;
//        return retval;
//      }
// }
// ```
/// `GetPercentIdentity(numerator, denominator)`.
pub(crate) fn percent_identity(numerator: i64, denominator: i64) -> f64 {
    if numerator == denominator {
        100.0
    } else {
        100.0 * numerator as f64 / denominator as f64
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4317-4356
// ```c++
// double CAlignFormatUtil::GetSeqAlignSetCalcPercentIdent(const CSeq_align_set& aln, bool do_translation)
// {
//     ...
//     if(aln.Get().empty())
//       return -1;
//
//     double highest_bits = 0;
//     int highest_length = 1;
//     int highest_ident = 0;
//     ...
//     ITERATE(CSeq_align_set::Tdata, iter, aln.Get()) {
//         int align_length = CAlignFormatUtil::GetAlignmentLength(**iter, do_translation);
//
//         CAlignFormatUtil::GetAlnScores(**iter, score, bits, evalue, sum_n,
//                                    num_ident, use_this_gi);
//         ...
//         if (bits > highest_bits) {  // this is the replacement condition (WB-1175)
//             highest_length = align_length;
//             highest_ident = num_ident;
// ///         highest_identity = 100*num_ident/align_length;
//             highest_bits = bits;
//         }
//     }
//
//     double percent_identity = CAlignFormatUtil::GetPercentIdentity(highest_ident, highest_length);
//     return percent_identity;
// }
// ```
/// `GetSeqAlignSetCalcPercentIdent` of the HSPs `order` of `hsps`.
fn set_percent_identity(
    hsps: &[HspScores],
    order: &[usize],
    tblastx: bool,
    translation: bool,
) -> f64 {
    if order.is_empty() {
        return -1.0;
    }
    let (mut highest_bits, mut highest_length, mut highest_ident) = (0.0, 1, 0);
    for &k in order {
        let hsp = &hsps[k];
        if hsp.bits > highest_bits {
            highest_length = hsp.alignment_length(tblastx, translation);
            highest_ident = hsp.num_ident;
            highest_bits = hsp.bits;
        }
    }
    percent_identity(highest_ident, highest_length)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:2969-2992
// ```c++
// int CAlignFormatUtil::GetMasterCoverage(const CSeq_align_set& alnset)
// {
//
//     list<CRange<TSeqPos> > merge_list;
//
//     list<CRange<TSeqPos> > temp;
//     ITERATE(CSeq_align_set::Tdata, iter, alnset.Get()) {
//         CRange<TSeqPos> seq_range = (*iter)->GetSeqRange(0);
//         //for minus strand
//         if(seq_range.GetFrom() > seq_range.GetTo()){
//             seq_range.Set(seq_range.GetTo(), seq_range.GetFrom());
//         }
//         temp.push_back(seq_range);
//     }
//
//     temp.sort(FromRangeAscendingSort);
//
//     merge_list = s_MergeRangeList(temp);
//
//     int master_covered_lenghth = 0;
//     ITERATE(list<CRange<TSeqPos> >, iter, merge_list) {
//         master_covered_lenghth += iter->GetLength();
//     }
//     return master_covered_lenghth;
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:2940-2963
// ```c++
// static list<CRange<TSeqPos> > s_MergeRangeList(list<CRange<TSeqPos> > &source)
// {
//     ...
//     ITERATE(list<CRange<TSeqPos> >, iter, source) {
//
//         if (is_first) {
//             merge_list.push_back(*iter);
//             is_first= false;
//             prev_range = *iter;
//         } else {
//             if (prev_range.IntersectingWith(*iter)) {
//                 merge_list.pop_back();
//                 CRange<TSeqPos> temp_range = prev_range.CombinationWith(*iter);
//                 merge_list.push_back(temp_range);
//                 prev_range = temp_range;
//             } else {
//                 merge_list.push_back(*iter);
//                 prev_range = *iter;
//             }
//         }
//
//     }
//     return merge_list;
// }
// ```
/// `GetMasterCoverage` (and the `master_covered_length` of `GetSeqAlignCoverageParams`, which
/// merges the same ranges, align_format_util.cpp:2997-3010): the number of query positions in
/// the union of the HSPs' query ranges. The ranges sorted by their start and merged where they
/// intersect cover each position once, in any order of the HSPs. NCBI sums the lengths in an
/// `int`; the sum here is a `u64` (larger than any query).
pub(crate) fn query_coverage<'a>(hsps: impl IntoIterator<Item = &'a HspScores>) -> u64 {
    let mut ranges: Vec<(u64, u64)> = hsps
        .into_iter()
        .map(|hsp| (hsp.query_from, hsp.query_to))
        .collect();
    ranges.sort_unstable();
    let mut covered = 0u64;
    let mut merged: Option<(u64, u64)> = None;
    for (from, to) in ranges {
        merged = match merged {
            Some((first, last)) if from <= last => Some((first, last.max(to))),
            Some((first, last)) => {
                covered += last - first + 1;
                Some((from, to))
            }
            None => Some((from, to)),
        };
    }
    if let Some((first, last)) = merged {
        covered += last - first + 1;
    }
    covered
}

/// The order of one query's hits after `sort_query_hits`: the hits (indices of the input) in
/// their new order, and each hit's HSPs (indices of that hit's HSPs, by input hit).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SortedHits {
    pub(crate) hits: Vec<usize>,
    pub(crate) hsps: Vec<Vec<usize>>,
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3103-3118
// ```c++
// CRef<CSeq_align_set>
// CAlignFormatUtil::SortSeqalignForSortableFormat(CSeq_align_set& aln_set,
//                                              bool nuc_to_nuc_translation,
//                                              int hit_sort,
//                                              int hsp_sort) {
//
//     if (hit_sort <= eEvalue && hsp_sort <= eHspEvalue) {
//        return (CRef<CSeq_align_set>) &aln_set;
//     }
//
// //  seqalign_vec[0] = const_cast<CSeq_align_set*>(&aln_set);
//     list< CRef<CSeq_align_set> > seqalign_hit_total_list = SortOneSeqalignForSortableFormat(aln_set,
//                                                             nuc_to_nuc_translation,
//                                                             hit_sort,
//                                                             hsp_sort);
//     return HitListToHspList(seqalign_hit_total_list);
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3065-3100
// ```c++
//     list< CRef<CSeq_align_set> > seqalign_hit_total_list;
//     list< CRef<CSeq_align_set> > seqalign_hit_list;
//     HspListToHitList(seqalign_hit_list, source);
//
//     if (hit_sort == eTotalScore) {
//         seqalign_hit_list.sort(SortHitByTotalScoreDescending);
//     } else if (hit_sort == eHighestScore) {
//         seqalign_hit_list.sort(CAlignFormatUtil::SortHitByScoreDescending);
//     } else if (hit_sort == ePercentIdentity) {
//         SortHitByPercentIdentityDescending(seqalign_hit_list,
//                                                nuc_to_nuc_translation);
//     } else if (hit_sort == eQueryCoverage) {
//         seqalign_hit_list.sort(SortHitByMasterCoverageDescending);
//     }
//
//     ITERATE(list< CRef<CSeq_align_set> >, iter2, seqalign_hit_list) {
//         CRef<CSeq_align_set> temp(*iter2);
//         if (hsp_sort == eQueryStart) {
//             temp->Set().sort(SortHspByMasterStartAscending);
//         } else if (hsp_sort == eHspPercentIdentity) {
//             temp->Set().sort(SortHspByPercentIdentityDescending);
//         } else if (hsp_sort == eScore) {
//             temp->Set().sort(SortHspByScoreDescending);
//         } else if (hsp_sort == eSubjectStart) {
//             temp->Set().sort(SortHspBySubjectStartAscending);
//
//         }
//         seqalign_hit_total_list.push_back(temp);
//     }
//     return seqalign_hit_total_list;
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1249-1257
// ```c++
// SortHitByPercentIdentityDescending(list< CRef<CSeq_align_set> >&
//                                    seqalign_hit_list,
//                                    bool do_translation
//                                    )
// {
//
//     kTranslation = do_translation;
//     seqalign_hit_list.sort(SortHitByPercentIdentityDescendingEx);
// }
// ```
/// `SortSeqalignForSortableFormat(aln_set, tblastx, hit_sort, hsp_sort)` of one query whose
/// hits (`HspListToHitList`: one per subject, as the report groups them) are `hits`, each with its
/// HSPs in the report's order; `None` when `sort` sorts nothing (no work, the order stays).
///
/// NCBI's `kTranslation` is a static that only the hit sort by percent identity sets, to
/// `m_Program == "tblastx"`, before it sorts each query's hits; it is zero-initialised, so it is
/// `tblastx && -sorthits 3` for the whole run, and the HSP sort by percent identity reads it too.
pub(crate) fn sort_query_hits(
    sort: ReportSort,
    tblastx: bool,
    hits: &[Vec<HspScores>],
) -> Option<SortedHits> {
    if !sort.sorts() {
        return None;
    }
    let translation = tblastx && sort.hits == Some(3);
    let mut orders: Vec<Vec<usize>> = hits.iter().map(|hsps| (0..hsps.len()).collect()).collect();
    let mut order: Vec<usize> = (0..hits.len()).collect();
    match sort.hits {
        Some(2) => sort_hits_by_total_score(hits, &orders, &mut order),
        Some(1) => sort_hits_by_score(hits, &mut orders, &mut order),
        Some(3) => {
            sort_hits_by_percent_identity(hits, &mut orders, &mut order, tblastx, translation)
        }
        Some(4) => sort_hits_by_query_coverage(hits, &orders, &mut order),
        _ => {}
    }
    for (hsps, hsp_order) in hits.iter().zip(orders.iter_mut()) {
        match sort.hsps {
            Some(2) => hsp_order.sort_by(|&x, &y| hsp_by_query_start_ascending(&hsps[x], &hsps[y])),
            Some(3) => sort_hsps_by_percent_identity(hsps, hsp_order, tblastx, translation),
            Some(1) => hsp_order.sort_by(|&x, &y| hsp_by_score_descending(&hsps[x], &hsps[y])),
            Some(4) => {
                hsp_order.sort_by(|&x, &y| hsp_by_subject_start_ascending(&hsps[x], &hsps[y]))
            }
            _ => {}
        }
    }
    Some(SortedHits {
        hits: order,
        hsps: orders,
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1597-1622
// ```c++
//     double total_bits1 = 0, total_bits2 = 0;
//
//     ITERATE(CSeq_align_set::Tdata, iter, info1->Get()) {
//         CAlignFormatUtil::GetAlnScores(**iter, score1, bits, evalue,
//                                        sum_n, num_ident, use_this_gi);
//         total_bits1 += bits;
//     }
//     ...
//     return total_bits1 >= total_bits2;
// ```
/// The sum of the bit scores of the HSPs `order` of `hsps`, added in that order from 0
/// (`SortHitByTotalScoreDescending` compares two with `>=`).
pub(crate) fn total_score(hsps: &[HspScores], order: &[usize]) -> f64 {
    order.iter().fold(0.0, |total, &k| total + hsps[k].bits)
}

/// `seqalign_hit_list.sort(SortHitByTotalScoreDescending)` of the hits `order` (in the search's
/// order), each with its HSPs in `orders`: the higher total first. The comparator's `>=` holds
/// both ways for equal totals; `std::list::sort` with it ends equal totals in the reverse of their
/// order (seen in every list of 0 to 300 elements that g++ 13's `std::list::sort` sorted with a
/// `>=` comparator): the list reversed, then sorted stably by total.
fn sort_hits_by_total_score(hits: &[Vec<HspScores>], orders: &[Vec<usize>], order: &mut [usize]) {
    let totals: Vec<f64> = hits
        .iter()
        .zip(orders)
        .map(|(hsps, hsp_order)| total_score(hsps, hsp_order))
        .collect();
    order.reverse();
    order.sort_by(|&a, &b| cmp_f64(totals[b], totals[a]));
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1298-1316
// ```c++
// bool CAlignFormatUtil::
// SortHitByScoreDescending(const CRef<CSeq_align_set>& info1,
//                          const CRef<CSeq_align_set>& info2)
// {
//     CRef<CSeq_align_set> i1(info1), i2(info2);
//
//     i1->Set().sort(SortHspByScoreDescending);
//     i2->Set().sort(SortHspByScoreDescending);
//     ...
//     GetAlnScores(*(info1->Get().front()), score1,  bits1, evalue1, sum_n1, num_ident1, use_this_gi1);
//     GetAlnScores(*(info2->Get().front()), score2,  bits2, evalue2, sum_n2, num_ident2, use_this_gi2);
//     return bits1 > bits2;
// }
// ```
/// `seqalign_hit_list.sort(SortHitByScoreDescending)` of the hits `order`, each with its HSPs in
/// `orders`. Each comparison sorts both hits' HSPs by score, and the hits keep that order of their
/// HSPs. With two hits or more every hit takes part in a comparison (in every list of 2 to 300
/// elements that g++ 13's `std::list::sort` sorted), so every hit's HSPs end sorted by score; a
/// stable sort of HSPs already in that order leaves them as they are, so one sort per hit gives the
/// same order. One hit is not compared and keeps its order. The hits then go by their first HSP's
/// bit score, the highest first (`>`, a strict weak order: equal scores keep their order).
fn sort_hits_by_score(hits: &[Vec<HspScores>], orders: &mut [Vec<usize>], order: &mut [usize]) {
    if order.len() < 2 {
        return;
    }
    for &hit in order.iter() {
        let hsps = &hits[hit];
        orders[hit].sort_by(|&x, &y| hsp_by_score_descending(&hsps[x], &hsps[y]));
    }
    order.sort_by(|&a, &b| cmp_f64(hits[b][orders[b][0]].bits, hits[a][orders[a][0]].bits));
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1543-1595
// ```c++
// static void s_CalcAlnPercentIdent(const CRef<CSeq_align_set>& info1,
//                                const CRef<CSeq_align_set>& info2,
//                                double &percentIdent1,
//                                double &percentIdent2)
// {
//
//     CRef<CSeq_align_set> i1(info1), i2(info2);
//     percentIdent1 = -1;
//     percentIdent2 = -1;
//
//     i1->Set().sort(CAlignFormatUtil::SortHspByPercentIdentityDescending);
//     i2->Set().sort(CAlignFormatUtil::SortHspByPercentIdentityDescending);
//
//     percentIdent1 = CAlignFormatUtil::GetSeqAlignSetCalcPercentIdent(*info1, kTranslation);
//     percentIdent2 = CAlignFormatUtil::GetSeqAlignSetCalcPercentIdent(*info2, kTranslation);
//     return;
// }
// ...
//     unique_ptr<CAlignFormatUtil::SSeqAlignSetCalcParams> seqSetInfo1( CAlignFormatUtil::GetSeqAlignSetCalcParamsFromASN(*info1));
//     unique_ptr<CAlignFormatUtil::SSeqAlignSetCalcParams> seqSetInfo2( CAlignFormatUtil::GetSeqAlignSetCalcParamsFromASN(*info2));
//     double evalue1 = seqSetInfo1->evalue;
//     double evalue2 = seqSetInfo2->evalue;
//     double percentIdent1 = seqSetInfo1->percent_identity;
//     double  percentIdent2 = seqSetInfo2->percent_identity;
//
//     bool retval = false;
//     if(percentIdent1 < 0 || percentIdent2 < 0) {
//         s_CalcAlnPercentIdent(info1, info2,percentIdent1,percentIdent2);
//     }
//     if(percentIdent1 > 0 &&percentIdent2 > 0) {
//         if (percentIdent1 == percentIdent2) {
//             retval = evalue1 < evalue2;
//
//         } else {
//             retval = percentIdent1 >= percentIdent2;
//         }
//     } else {
//         retval = evalue1 < evalue2;
//     }
//     return retval;
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4474-4485
// ```c++
//     bool hasScore = false;
//     double evalue = -1;
//     double bitScore = -1;
//     double totalBitScore = -1;
//     int percentCoverage = -1;
//     double percentIdent = -1;
//     ...
//     hasScore = s_GetBlastScore(aln.GetScore(),evalue,bitScore, totalBitScore,percentCoverage,percentIdent,hspNum,totalLen,rawScore,sum_n,use_this_gi);
// ```
/// `SortHitByPercentIdentityDescendingEx` of two hits whose `GetSeqAlignSetCalcPercentIdent` are
/// `percent1` and `percent2`. The Seq-aligns of a BLAST+ search carry none of the "seq_*" scores
/// that `GetSeqAlignSetCalcParamsFromASN` reads, so both E-values and percentages start at -1 and
/// the comparison takes the percentages of `s_CalcAlnPercentIdent`; equal percentages, or one not
/// above 0, compare their E-values of -1 (`-1 < -1`, false).
fn hit_by_percent_identity_descending(percent1: f64, percent2: f64) -> bool {
    let (evalue1, evalue2) = (-1.0f64, -1.0f64);
    if percent1 > 0.0 && percent2 > 0.0 {
        if percent1 == percent2 {
            evalue1 < evalue2
        } else {
            percent1 >= percent2
        }
    } else {
        evalue1 < evalue2
    }
}

/// `SortHitByPercentIdentityDescending(seqalign_hit_list, kTranslation)` of the hits `order`, each
/// with its HSPs in `orders` (`translation`: NCBI's `kTranslation`).
///
/// Each comparison sorts both hits' HSPs by percent identity (`s_CalcAlnPercentIdent`), and the
/// hits keep that order. As with `sort_hits_by_score`, two hits or more mean that every hit's
/// HSPs end sorted, and one sort per hit gives that order (a second sort leaves a list that its
/// comparator orders as it is); one hit is not compared. The hits then go by their percentages, the
/// higher first: `>=` decides only unequal percentages, where it is `>`; equal percentages tie. A
/// percentage not above 0 (the first HSP of the highest score has no identities) ties with every
/// hit, so the comparator is a strict weak order of the list only when there is no such hit or the
/// other hits' percentages are all equal (`holds`), and NCBI's order is then the stable one.
/// Otherwise NCBI's order is the one that its C++ library's merge procedure gives, and LOSAT puts
/// the hits with a percentage above 0 first, the higher first, then the others (each group
/// stable). The same holds within a hit whose HSPs `sort_hsps_by_percent_identity` cannot order as
/// NCBI does: NCBI sorts them again at each comparison of the hit.
fn sort_hits_by_percent_identity(
    hits: &[Vec<HspScores>],
    orders: &mut [Vec<usize>],
    order: &mut [usize],
    tblastx: bool,
    translation: bool,
) {
    if order.len() < 2 {
        return;
    }
    for &hit in order.iter() {
        sort_hsps_by_percent_identity(&hits[hit], &mut orders[hit], tblastx, translation);
    }
    let percents: Vec<f64> = hits
        .iter()
        .zip(orders.iter())
        .map(|(hsps, hsp_order)| set_percent_identity(hsps, hsp_order, tblastx, translation))
        .collect();
    let above: Vec<f64> = order
        .iter()
        .map(|&hit| percents[hit])
        .filter(|&percent| percent > 0.0)
        .collect();
    let holds = above.len() == order.len() || above.iter().all(|&percent| percent == above[0]);
    if holds {
        sort_by_comparator(order, |a, b| {
            hit_by_percent_identity_descending(percents[a], percents[b])
        });
    } else {
        order.sort_by(|&a, &b| {
            let (percent1, percent2) = (percents[a], percents[b]);
            (percent2 > 0.0)
                .cmp(&(percent1 > 0.0))
                .then_with(|| cmp_f64(percent2, percent1))
        });
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1318-1344
// ```c++
//     int cov1 = GetMasterCoverage(*info1);
//     int cov2 = GetMasterCoverage(*info2);
//     bool retval = false;
//
//     if (cov1 > cov2) {
//         retval = cov1 > cov2;
//     } else if (cov1 == cov2) {
//         ...
//         GetAlnScores(*(info1->Get().front()), score1,  bits1, evalue1, sum_n1, num_ident1, use_this_gi1);
//         GetAlnScores(*(info2->Get().front()), score2,  bits2, evalue2, sum_n2, num_ident2, use_this_gi2);
//         retval = evalue1 < evalue2;
//     }
//
//     return retval;
// ```
/// `seqalign_hit_list.sort(SortHitByMasterCoverageDescending)` of the hits `order`, each with its
/// HSPs in `orders`: the higher query coverage first, equal coverages by their first HSP's lower
/// E-value (a strict weak order: equal coverages and E-values keep their order).
fn sort_hits_by_query_coverage(
    hits: &[Vec<HspScores>],
    orders: &[Vec<usize>],
    order: &mut [usize],
) {
    if order.len() < 2 {
        return;
    }
    let coverages: Vec<u64> = hits.iter().map(query_coverage).collect();
    order.sort_by(|&a, &b| {
        coverages[b]
            .cmp(&coverages[a])
            .then_with(|| cmp_f64(hits[a][orders[a][0]].evalue, hits[b][orders[b][0]].evalue))
    });
}

#[cfg(test)]
#[path = "hit_sort_tests.rs"]
mod hit_sort_tests;
