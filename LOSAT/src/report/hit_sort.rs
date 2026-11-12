//! The order of one query's hits and HSPs in the pairwise report (outfmt 0) with `-sorthits` and
//! `-sorthsps` (NCBI's `CAlignFormatUtil::SortSeqalignForSortableFormat`), shared by BLASTN,
//! BLASTP, BLASTX, TBLASTN and TBLASTX, and the values that the description table's columns of
//! `-sorthits` read.
//!
//! NCBI sorts `std::list`s with comparators of which some are not strict weak orders (`>=`) and
//! some sort the HSP lists of both hits they compare as a side effect, so the order is the one
//! that libstdc++'s merge sort gives with that comparator: `list_sort` is that algorithm, call
//! for call, and the comparators keep their side effects on the hits' HSP orders.

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
        let h = &hit.hit;
        let zero_based = |position: usize| (position as u64).saturating_sub(1);
        Self {
            bits: h.bit_score,
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
            evalue: crate::report::serial::seqalign::seqalign_evalue(h.e_value),
            num_ident: hit.seqalign_num_ident as i64,
            length: h.length as i64,
            query_from: zero_based(h.q_start.min(h.q_end)),
            query_to: zero_based(h.q_start.max(h.q_end)),
            subject_from: zero_based(h.s_start.min(h.s_end)),
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

// libstdc++ (GCC 13, the C++ library of NCBI BLAST+'s Linux builds): bits/list.tcc:614-662
// ```c++
//       sort(_StrictWeakOrdering __comp)
//       {
// 	// Do nothing if the list has length 0 or 1.
// 	if (this->_M_impl._M_node._M_next != &this->_M_impl._M_node
// 	    && this->_M_impl._M_node._M_next->_M_next != &this->_M_impl._M_node)
// 	{
// 	  using __detail::_Scratch_list;
// 	  _Scratch_list __carry;
// 	  _Scratch_list __tmp[64];
// 	  _Scratch_list* __fill = __tmp;
// 	  _Scratch_list* __counter;
// 	  ...
// 	      do
// 		{
// 		  __carry._M_take_one(begin()._M_node);
//
// 		  for(__counter = __tmp;
// 		      __counter != __fill && !__counter->empty();
// 		      ++__counter)
// 		    {
//
// 		      __counter->merge(__carry, __ptr_comp);
// 		      __carry.swap(*__counter);
// 		    }
// 		  __carry.swap(*__counter);
// 		  if (__counter == __fill)
// 		    ++__fill;
// 		}
// 	      while ( !empty() );
//
// 	      for (__counter = __tmp + 1; __counter != __fill; ++__counter)
// 		__counter->merge(__counter[-1], __ptr_comp);
// 	      __fill[-1].swap(this->_M_impl._M_node);
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3072-3075
// ```c++
//     HspListToHitList(seqalign_hit_list, source);
//
//     if (hit_sort == eTotalScore) {
//         seqalign_hit_list.sort(SortHitByTotalScoreDescending);
// ```
/// libstdc++'s `std::list::sort(comp)` of `list`: the bottom-up merge sort of 64 buckets, with
/// `comp(a, b)` (whether `a` goes before `b`) called on the same elements in the same order as
/// libstdc++ calls it (`comp` may change state that later calls read). A list of 0 or 1
/// elements is left as it is without a call. The Vec stands for the list: `merge` builds the
/// merged list as libstdc++ splices it.
pub(crate) fn list_sort<T>(list: &mut Vec<T>, mut comp: impl FnMut(&T, &T) -> bool) {
    if list.len() < 2 {
        return;
    }
    let mut tmp: Vec<Vec<T>> = (0..64).map(|_| Vec::new()).collect();
    let mut fill = 0usize;
    for item in std::mem::take(list) {
        let mut carry = vec![item];
        let mut counter = 0usize;
        while counter != fill && !tmp[counter].is_empty() {
            merge(&mut tmp[counter], &mut carry, &mut comp);
            std::mem::swap(&mut carry, &mut tmp[counter]);
            counter += 1;
        }
        std::mem::swap(&mut carry, &mut tmp[counter]);
        if counter == fill {
            fill += 1;
        }
    }
    // `__counter->merge(__counter[-1])`: bucket `counter` holds the older elements (`this`),
    // bucket `counter - 1` the newer ones (`x`).
    for counter in 1..fill {
        let (newer, older) = tmp.split_at_mut(counter);
        merge(&mut older[0], &mut newer[counter - 1], &mut comp);
    }
    *list = std::mem::take(&mut tmp[fill - 1]);
}

// libstdc++ (GCC 13): bits/stl_list.h:194-214
// ```c++
// 	merge(_List_node_base& __x, _Cmp __comp)
// 	{
// 	  _List_node_base* __first1 = _M_next;
// 	  _List_node_base* const __last1 = this;
// 	  _List_node_base* __first2 = __x._M_next;
// 	  _List_node_base* const __last2 = std::__addressof(__x);
//
// 	  while (__first1 != __last1 && __first2 != __last2)
// 	    {
// 	      if (__comp(__first2, __first1))
// 		{
// 		  _List_node_base* __next = __first2->_M_next;
// 		  __first1->_M_transfer(__first2, __next);
// 		  __first2 = __next;
// 		}
// 	      else
// 		__first1 = __first1->_M_next;
// 	    }
// 	  if (__first2 != __last2)
// 	    this->_M_transfer(__first2, __last2);
// 	}
// ```
/// `this.merge(x, comp)`: an element of `x` goes before the element of `this` it is compared
/// with only when `comp(x's, this's)`; the rest of `x` goes after `this`. `x` is left empty.
fn merge<T>(this: &mut Vec<T>, x: &mut Vec<T>, comp: &mut impl FnMut(&T, &T) -> bool) {
    if x.is_empty() {
        return;
    }
    let mut merged = Vec::with_capacity(this.len() + x.len());
    let mut first1 = std::mem::take(this).into_iter().peekable();
    let mut first2 = std::mem::take(x).into_iter().peekable();
    while let (Some(older), Some(newer)) = (first1.peek(), first2.peek()) {
        if comp(newer, older) {
            merged.extend(first2.next());
        } else {
            merged.extend(first1.next());
        }
    }
    merged.extend(first1);
    merged.extend(first2);
    *this = merged;
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
/// `SortHspByScoreDescending`.
pub(crate) fn hsp_by_score_descending(info1: &HspScores, info2: &HspScores) -> bool {
    info1.bits > info2.bits
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
/// `SortHspByMasterStartAscending`.
pub(crate) fn hsp_by_query_start_ascending(info1: &HspScores, info2: &HspScores) -> bool {
    if info1.query_from == info2.query_from {
        info1.evalue < info2.evalue
    } else {
        info1.query_from < info2.query_from
    }
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
/// `SortHspBySubjectStartAscending`.
pub(crate) fn hsp_by_subject_start_ascending(info1: &HspScores, info2: &HspScores) -> bool {
    if info1.subject_from == info2.subject_from {
        info1.evalue < info2.evalue
    } else {
        info1.subject_from < info2.subject_from
    }
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
        Some(2) => list_sort(&mut order, |&a, &b| {
            total_score(&hits[a], &orders[a]) >= total_score(&hits[b], &orders[b])
        }),
        Some(1) => list_sort(&mut order, |&a, &b| {
            hit_by_score_descending(hits, &mut orders, a, b)
        }),
        Some(3) => list_sort(&mut order, |&a, &b| {
            hit_by_percent_identity_descending(hits, &mut orders, a, b, tblastx, translation)
        }),
        Some(4) => list_sort(&mut order, |&a, &b| {
            hit_by_query_coverage_descending(hits, &orders, a, b)
        }),
        _ => {}
    }
    for (hsps, hsp_order) in hits.iter().zip(orders.iter_mut()) {
        match sort.hsps {
            Some(2) => list_sort(hsp_order, |&x, &y| {
                hsp_by_query_start_ascending(&hsps[x], &hsps[y])
            }),
            Some(3) => list_sort(hsp_order, |&x, &y| {
                hsp_by_percent_identity_descending(&hsps[x], &hsps[y], tblastx, translation)
            }),
            Some(1) => list_sort(hsp_order, |&x, &y| {
                hsp_by_score_descending(&hsps[x], &hsps[y])
            }),
            Some(4) => list_sort(hsp_order, |&x, &y| {
                hsp_by_subject_start_ascending(&hsps[x], &hsps[y])
            }),
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
/// `SortHitByScoreDescending` of hits `a` and `b`: sorts both hits' HSPs by score (kept: the
/// order of their HSPs afterwards), then compares their first HSPs' bit scores.
fn hit_by_score_descending(
    hits: &[Vec<HspScores>],
    orders: &mut [Vec<usize>],
    a: usize,
    b: usize,
) -> bool {
    for hit in [a, b] {
        let hsps = &hits[hit];
        list_sort(&mut orders[hit], |&x, &y| {
            hsp_by_score_descending(&hsps[x], &hsps[y])
        });
    }
    hits[a][orders[a][0]].bits > hits[b][orders[b][0]].bits
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
/// `SortHitByPercentIdentityDescendingEx` of hits `a` and `b`. The Seq-aligns of a BLAST+ search
/// carry none of the "seq_*" scores that `GetSeqAlignSetCalcParamsFromASN` reads, so both
/// E-values and percentages start at -1: the comparison sorts both hits' HSPs by percent
/// identity (kept) and takes their `GetSeqAlignSetCalcPercentIdent`; equal percentages, or one
/// not above 0, compare their E-values of -1 (`-1 < -1`, false).
fn hit_by_percent_identity_descending(
    hits: &[Vec<HspScores>],
    orders: &mut [Vec<usize>],
    a: usize,
    b: usize,
    tblastx: bool,
    translation: bool,
) -> bool {
    let (evalue1, evalue2) = (-1.0f64, -1.0f64);
    for hit in [a, b] {
        let hsps = &hits[hit];
        list_sort(&mut orders[hit], |&x, &y| {
            hsp_by_percent_identity_descending(&hsps[x], &hsps[y], tblastx, translation)
        });
    }
    let percent1 = set_percent_identity(&hits[a], &orders[a], tblastx, translation);
    let percent2 = set_percent_identity(&hits[b], &orders[b], tblastx, translation);
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
/// `SortHitByMasterCoverageDescending` of hits `a` and `b` (their first HSPs in their current
/// orders).
fn hit_by_query_coverage_descending(
    hits: &[Vec<HspScores>],
    orders: &[Vec<usize>],
    a: usize,
    b: usize,
) -> bool {
    let cov1 = query_coverage(&hits[a]);
    let cov2 = query_coverage(&hits[b]);
    if cov1 > cov2 {
        true
    } else if cov1 == cov2 {
        hits[a][orders[a][0]].evalue < hits[b][orders[b][0]].evalue
    } else {
        false
    }
}

#[cfg(test)]
#[path = "hit_sort_tests.rs"]
mod hit_sort_tests;
