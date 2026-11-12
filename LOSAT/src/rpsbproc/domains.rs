//! The domain hits of one protein query: the alignments of an archive turned into rpsbproc's
//! `TDomSeqAlignment`s, their classification into regions, representative, standard and full
//! sets (`CCdAlignProcessor::Calculate`), and the records that the DOMAINS output prints
//! (`TDomSeqAlignIndex::CreateRecordSets`).

use std::collections::HashMap;

use super::archive::{Align, DenseSeg};
use super::data::{dom_acc_type, CddData, TOTALSIGS};
use super::segset::SegSet;
use super::sites::ResiduePos;
use crate::utils::std_sort::std_sort;

/// `TDomSeqAlignment::GAP_THRESHOLD`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:1203-1205
// ```c++
// struct TDomSeqAlignment: public TSeqAlignment
// {
// 	static constexpr const SeqLen_t GAP_THRESHOLD = 35;
// ```
const GAP_THRESHOLD: i64 = 35;

/// `TDomSeqAlignment`, with the members that the domain output reads.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:1084-1123,1206-1221
// ```c++
// 	// -- scores
// 	SeqLen_t m_uiAlignedLen;
// 	double m_dAlignedPct;
// 	int m_iScore;
// 	double m_dEValue;
// 	double m_dBitScore;
// 	int m_iNumIdent;
// 	double m_dSeqIdentity;
// 	...
// 	SeqPos_t m_iFrom;
// 	SeqPos_t m_iTo;
// // -- processed align info
// 	std::vector<SeqPos_t> m_vecMStarts;	//regardless the align type, MStarts and SStarts always be protein-to-protein.
// 	std::vector<SeqPos_t> m_vecSStarts;
// 	std::vector<SeqLen_t> m_vecLens;
// 	...
// 	PssmId_t m_uiPssmId;
// 	int m_iRegionIdx;	//region on the query sequence.
//
// 	double m_dNMissing;
// 	double m_dCMissing;
// 	bool m_bSpecQualified;	//higher bitscore than threshold
// 	int m_iRepClass;	//single and multi -- sort separately
// 	bool m_bRep;
// 	...
// 	TDomSeqAlignment(void):
// 		m_uiPssmId(0), m_iRegionIdx(0), m_dNMissing(0.0), m_dCMissing(0.0), m_bSpecQualified(false), m_iRepClass(0), m_bRep(0), m_bLifted(0), m_bSuppressed(0)
// 	{};
// ```
#[derive(Clone, Debug)]
pub(crate) struct DomAlign {
    pub(crate) pssm_id: u32,
    pub(crate) from: i32,
    pub(crate) to: i32,
    pub(crate) m_starts: Vec<i32>,
    pub(crate) s_starts: Vec<i32>,
    pub(crate) lens: Vec<u32>,
    pub(crate) e_value: f64,
    pub(crate) bit_score: f64,
    pub(crate) region: i32,
    pub(crate) n_missing: f64,
    pub(crate) c_missing: f64,
    pub(crate) spec_qualified: bool,
    pub(crate) rep_class: i32,
    pub(crate) rep: bool,
    /// `m_ClipSet`: the query ranges (as positions on the model, `ori_from` side) that a
    /// structural-domain hit keeps for its motifs (set by `calculate`, empty otherwise).
    pub(crate) clip_set: SegSet,
}

/// `TDomSeqAlignIndex`: the indices of the alignments in each set.
#[derive(Clone, Debug, Default)]
pub(crate) struct AlignIndex {
    pub(crate) sorted: Vec<usize>,
    pub(crate) concise: Vec<usize>,
    pub(crate) standard: Vec<usize>,
    pub(crate) qualified_feat: Vec<usize>,
    /// `m_vecSDIndice`: the structural-domain hits that keep a part of the query.
    pub(crate) structural: Vec<usize>,
}

/// `TDataModes`: `rep`, `std`, `full`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum DataMode {
    Rep,
    Std,
    Full,
}

/// `TDomSeqAlignIndex::__TCdAlignRecord` of the DOMAINS output.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Record {
    /// The alignment (`pAlign`).
    pub(crate) align: usize,
    /// The domain (`pCdInfo`, index in `CddData::domains`).
    pub(crate) domain: usize,
    /// The superfamily (`pClst`, index in `CddData::clusters`).
    pub(crate) cluster: Option<usize>,
}

/// `CleanAlignment`: the segments where both rows are aligned.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/impl/basedata_impl.hpp:2-19
// ```c++
// template<typename NCBISeqPosType, typename NCBISeqLenType>
// void CleanAlignment(const std::vector<NCBISeqPosType> & rSrcStarts, const std::vector<NCBISeqLenType> & rSrcLens, std::vector<SeqPos_t>& rMStarts, std::vector<SeqPos_t>& rSStarts, std::vector<SeqLen_t>& rLens)
// {
// 	rMStarts.clear();
// 	rSStarts.clear();
// 	rLens.clear();
//
// 	for (size_t i = 0; i < rSrcLens.size(); ++i)
// 	{
// 		size_t ii = i + i;
// 		if (rSrcStarts[ii] >=0 && rSrcStarts[ii + 1] >= 0)
// 		{
// 			rMStarts.push_back(rSrcStarts[ii]);
// 			rSStarts.push_back(rSrcStarts[ii + 1]);
// 			rLens.push_back(rSrcLens[i]);
// 		}
// 	}
// }
// ```
// (`starts` shorter than twice `lens` is read past its end by NCBI: refused by the caller.)
fn clean_alignment(seg: &DenseSeg, align: &mut DomAlign) {
    for (i, &len) in seg.lens.iter().enumerate() {
        let (m, s) = (seg.starts[2 * i], seg.starts[2 * i + 1]);
        if m >= 0 && s >= 0 {
            align.m_starts.push(m);
            align.s_starts.push(s);
            align.lens.push(len);
        }
    }
}

/// `ParseAlignSegs` of a `Dense-seg` and `ParseAlignScores`: the alignment, without its
/// PSSM-ID (`ProcessCDAlign` sets it).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_blast.cpp:9-43
// ```c++
// CRef<CSeq_id> ParseAlignSegs(const CSeq_align::TSegs &c_segs, TSeqAlignment &dst)
// {
// 	CRef<CSeq_id> hit_id(nullptr);
//
// 	dst.Reset();
// 	CSeq_align::TSegs::E_Choice segs_type = c_segs.Which();
//
// 	dst.m_iFrom = k_INT_MAX;
// 	dst.m_iTo = 0;
// 	ENa_strand qstrand = eNa_strand_other;
//
// 	switch (segs_type)
// 	{
//
// 	case CSeq_align::TSegs::e_Denseg:	//dense seg is always pr::pr
//
// 		{
// 			const CDense_seg &denseg = c_segs.GetDenseg();
// 			const auto &ids = denseg.GetIds();
// 			if (!ids.empty())
// 			{
// 				hit_id = SerialClone(*(ids[1]));
// 				CleanAlignment(denseg.GetStarts(), denseg.GetLens(), dst.m_vecMStarts, dst.m_vecSStarts, dst.m_vecLens);
//
//
// 				size_t segs = dst.m_vecLens.size();
//
// 				if (segs > 0)
// 				{
// 					dst.m_iFrom = dst.m_vecMStarts[0];
// 					dst.m_iTo = dst.m_vecMStarts[segs - 1] + dst.m_vecLens[segs - 1] - 1;
// 				}
// 			}
// 		}
// 		break;
// ```
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_blast.cpp:214-262
// (`ParseAlignScores`, quoted at `archive::ArchiveStream::align_score`, which keeps the scores
// it reads.)
// (`m_vecMStarts[...] + m_vecLens[...] - 1` adds an `int` and an `unsigned int`: the sum wraps
// as an `unsigned int` and is stored back in an `int`.)
pub(crate) fn parse_align(seg: &DenseSeg, align: &Align, pssm_id: u32) -> DomAlign {
    let mut dom = DomAlign {
        pssm_id,
        from: i32::MAX,
        to: 0,
        m_starts: Vec::new(),
        s_starts: Vec::new(),
        lens: Vec::new(),
        e_value: align.e_value.unwrap_or(0.0),
        bit_score: align.bit_score.unwrap_or(0.0),
        region: 0,
        n_missing: 0.0,
        c_missing: 0.0,
        spec_qualified: false,
        rep_class: 0,
        rep: false,
        clip_set: SegSet::new(),
    };
    clean_alignment(seg, &mut dom);
    if let (Some(&first), Some(&last), Some(&last_len)) =
        (dom.m_starts.first(), dom.m_starts.last(), dom.lens.last())
    {
        dom.from = first;
        dom.to = (last as u32).wrapping_add(last_len).wrapping_sub(1) as i32;
    }
    dom
}

/// `TSeqAlignment::CalcMasterGaps(GAP_THRESHOLD)`: the query ranges between segments that are
/// at least `GAP_THRESHOLD` long.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:1324-1339
// ```c++
// void TSeqAlignment::CalcMasterGaps(SeqLen_t gapThreshold, CSegSet &segset) const
// {
// 	segset.Clear();
// 	if (!m_vecLens.empty())
// 	{
// 		size_t segs = m_vecLens.size();
// 		for (size_t i = 0; i < segs - 1; ++i)
// 		{
// 			SeqPos_t gapstart = m_vecMStarts[i] + (SeqPos_t)m_vecLens[i];
// 			SeqPos_t gaplen = m_vecMStarts[i + 1] - gapstart;
//
// 			if (gaplen >= (SeqPos_t)gapThreshold)	//consider a gap
// 				segset.AddSeg(gapstart, gapstart + gaplen - 1);
// 		}
// 	}
// }
// ```
// (`(SeqPos_t)m_vecLens[i]` converts an `unsigned int` to `int`, which wraps in NCBI's build;
// the `int` sums after it are computed in `i64`, `segset::Seg`.)
fn master_gaps(align: &DomAlign) -> SegSet {
    let mut gaps = SegSet::new();
    for i in 0..align.lens.len().saturating_sub(1) {
        let gap_start = i64::from(align.m_starts[i]) + i64::from(align.lens[i] as i32);
        let gap_len = i64::from(align.m_starts[i + 1]) - gap_start;
        if gap_len >= GAP_THRESHOLD {
            gaps.add(gap_start, gap_start + gap_len - 1);
        }
    }
    gaps
}

/// `TDomSrcCount`: how many hits of each source database a region took (one each).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:115,154-166
// ```c++
// array<int, TDomSrcCount::TOTALSRCS> TDomSrcCount::MAXCOUNTS = {1, 1, 1, 1, 1, 1, 1};
// ...
// bool TDomSrcCount::CountSrc(const string &acxn)
// {
// 	ESrcIdx srctype = DomAccType(acxn);
// 	if (TOTALSIGS == srctype) return true;	//unknown src always enter.
// 	if (srctype >= ePRK)
// 		srctype = ePRK;
//
// 	map<ESrcIdx, int> :: iterator iter = m_SrcCounter.emplace(srctype, 0).first;
//
// 	if (iter->second >= m_max_counts[srctype]) return false;
// 	++iter->second;
// 	return true;
// }
// ```
fn count_source(counts: &mut [i32; 7], accession: &[u8]) -> bool {
    const E_PRK: usize = 6;
    let source = dom_acc_type(accession);
    if source == TOTALSIGS {
        return true;
    }
    let source = source.min(E_PRK);
    if counts[source] >= 1 {
        return false;
    }
    counts[source] += 1;
    true
}

/// A hit as `Calculate` sorts it (`TDomAlignFacility`).
#[derive(Clone, Copy, Debug)]
struct Facility {
    align: usize,
    domain: usize,
    /// The gaps (`ulGapsIdx`, an index in the gap sets).
    gaps: Option<usize>,
}

/// `CCdAlignProcessor::Calculate` of a protein query: the alignments' scores against their
/// models, the regions, and the sets of `index`. Structural domains (`sd`) are set apart and
/// take no region (their clipping, `m_vecSDIndice` and `m_ClipSet`, serves the motif output).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_base.cpp:48-130
// ```c++
// void CCdAlignProcessor::Calculate(vector<TDomSeqAlignment> &aligns, const vector<size_t> &selIdx, TDomSeqAlignIndex &dst, vector<PssmId_t> &missed, const int * stdmax, size_t nmax) const
// {
//
// 	vector <CSegSet> vecGaps;
// 	vecGaps.emplace_back(CSegSet());	//sentinel
// 	size_t gapIdx = 0;
// 	vector <TDomAlignFacility> vecPrivileged, vecNonMulti, vecLongMultiDom, vecMultiDom, vecConcise, vecSDs;
// 	// -- added 12/6/2016: Concise view now should always include curated models (if available).
// 	if (!selIdx.empty())
// 	{
// 		...
// 		for (vector<size_t> :: const_iterator iterIdx = selIdx.begin(), iterIdxEnd = selIdx.end(); iterIdx != iterIdxEnd; ++iterIdx)
// 		{
// 			TDomSeqAlignment &rAlign = aligns[*iterIdx];
//
// 			if (rAlign.m_vecSStarts.empty()) continue;	//somehow live blast search contains empty aligns.
// 			const TDomain *pCdInfo = nullptr;
// 			...
// 			pCdInfo = m_pDomSrc->FindCdInfo(rAlign.m_uiPssmId);
// 			...
// 			if (nullptr == pCdInfo)	//unrecognized CD
// 			{
//
// 				missed.emplace_back(rAlign.m_uiPssmId);
// 				missedIdx.emplace_back(*iterIdx);
// 				//alignIdx.emplace_back(*iterIdx);
// 				continue;
// 			}
// 			rAlign.m_dSeqIdentity = (double)(rAlign.m_iNumIdent) / (double)(pCdInfo->m_uiLength) * 100.0;
//
// 			int iNMissing = *(rAlign.m_vecSStarts.begin());
// 			int iCMissing = pCdInfo->m_uiLength - (*(rAlign.m_vecSStarts.rbegin()) + *(rAlign.m_vecLens.rbegin()));
//
// 			rAlign.m_dNMissing = (double)(iNMissing) / (double)(pCdInfo->m_uiLength);
// 			rAlign.m_dCMissing = (double)(iCMissing) / (double)(pCdInfo->m_uiLength);
// 			rAlign.m_uiAlignedLen = pCdInfo->m_uiLength - iNMissing - iCMissing;
// 			rAlign.m_dAlignedPct = (double)(rAlign.m_uiAlignedLen) / (double)(pCdInfo->m_uiLength) * 100.0;
//
// 			rAlign.m_bSpecQualified = (pCdInfo->m_dMinBitScore > 0.0 && rAlign.m_dBitScore >= pCdInfo->m_dMinBitScore);
// 			// -- sorting
//
// 			pholder.pAlign = &(rAlign);
// 			pholder.ulIdx = *iterIdx;
// 			pholder.pCdInfo = pCdInfo;
// 			pholder.ulGapsIdx = string::npos;
//
// 			rAlign.CalcMasterGaps(TDomSeqAlignment::GAP_THRESHOLD, vecGaps[gapIdx]);
//
// 			if (!vecGaps[gapIdx].IsEmpty())
// 			{
// 				pholder.ulGapsIdx = gapIdx;
// 				++gapIdx;
// 				vecGaps.emplace_back(CSegSet());
// 			}
//
// 			if (pCdInfo->m_bIsStructDom)
// 			{
// 				rAlign.m_iRepClass = 0x2;
// 				vecSDs.emplace_back(pholder);
// 			}
// 			// -- from 4/4/2018 m_bMultiDom should be ignored.
// 			else// if (pCdInfo->m_bCurated || !pCdInfo->m_bMultiDom)	//curated or non-multi are mixed together
// 			{
// 				rAlign.m_iRepClass = 0;
// 				vecNonMulti.emplace_back(pholder);
// 			}
// 			// -- from 4/4/2018 m_bMultiDom should be ignored.
// 			//else
// 			//{
// 			//	rAlign.m_iRepClass = 1;
// 			//	vecMultiDom.emplace_back(pholder);
// 			//}
// 		}
// ```
// The models missing from the data are looked up again after `x_LoadMissingDomains`, which the
// offline processor does not override (`cdalignproc_base.hpp:123`, an empty body): they stay
// missing and their hits are dropped (cdalignproc_base.cpp:133-195). `Multidom` hits (`m_iRepClass`
// 1) are not made since 2018-04-04 (the commented lines above), so the output's `Multidom` type
// is never printed.
// (`m_uiLength - (sstart + len)` is `unsigned int` arithmetic stored in an `int`; the ratios
// divide doubles, so a model of length 0 gives an infinity or a NaN as NCBI's does.)
pub(crate) fn calculate(aligns: &mut [DomAlign], data: &CddData) -> AlignIndex {
    let mut index = AlignIndex::default();
    let mut gap_sets: Vec<SegSet> = Vec::new();
    let mut non_multi: Vec<Facility> = Vec::new();
    let mut structural: Vec<Facility> = Vec::new();
    for (i, align) in aligns.iter_mut().enumerate() {
        if align.s_starts.is_empty() {
            continue;
        }
        let Some(domain_index) = data.find_domain(align.pssm_id) else {
            continue;
        };
        let domain = &data.domains[domain_index];
        let length = domain.length;
        let n_missing = align.s_starts[0];
        let c_missing = length.wrapping_sub(
            (*align.s_starts.last().expect("not empty") as u32)
                .wrapping_add(*align.lens.last().expect("not empty")),
        ) as i32;
        align.n_missing = f64::from(n_missing) / f64::from(length);
        align.c_missing = f64::from(c_missing) / f64::from(length);
        align.spec_qualified =
            domain.min_bit_score > 0.0 && align.bit_score >= domain.min_bit_score;
        let gaps = master_gaps(align);
        let gaps = if gaps.is_empty() {
            None
        } else {
            gap_sets.push(gaps);
            Some(gap_sets.len() - 1)
        };
        let facility = Facility {
            align: i,
            domain: domain_index,
            gaps,
        };
        if domain.is_struct_dom {
            align.rep_class = 2;
            structural.push(facility);
        } else {
            align.rep_class = 0;
            non_multi.push(facility);
        }
    }
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_base.hpp:62-70
    // ```c++
    // 	struct TDomAlignOrderByEValue
    // 	{
    // 		bool operator () (const TDomAlignFacility &p1, const TDomAlignFacility &p2)
    // 		{
    // 			if (p1.pAlign->m_dEValue < p2.pAlign->m_dEValue) return true;
    // 			else if (p1.pAlign->m_dEValue > p2.pAlign->m_dEValue) return false;
    // 			else return (p1.pAlign->m_dBitScore > p2.pAlign->m_dBitScore);
    // 		}
    // 	};
    // ```
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_base.cpp:197-270
    // ```c++
    // 		map<int, TDomSrcCount> dimAccTypeCount;//, dimAccTypeCount_multi;	//one for single doms and one for multi-doms
    // 		pair<int, TDomSrcCount> srccnt(0, TDomSrcCount(stdmax, nmax));
    // 	  // -- sort
    // 		sort(vecNonMulti.begin(), vecNonMulti.end(), sortingEV);
    // 		for (vector <TDomAlignFacility>::const_iterator iter = vecNonMulti.begin(), iterEnd = vecNonMulti.end(); iterEnd != iter; ++iter)
    // 		{
    // 			CSegSet s_overlaps;	// to calculate combined overlap region
    //
    // 			CSegSet thisSegs;
    // 			thisSegs.AddSeg(iter->pAlign->m_iFrom, iter->pAlign->m_iTo);
    // 			if (string::npos != iter->ulGapsIdx)
    // 				thisSegs.Clip(vecGaps[iter->ulGapsIdx]);
    //
    // 			int iThisLength = (int)thisSegs.GetTotalResidues();
    // 			double dThisLength = (double)iThisLength;
    //
    // 			for (vector <TDomAlignFacility>::const_iterator iterRep = vecConcise.begin(), iterRepEnd = vecConcise.end(); iterRepEnd != iterRep; ++iterRep)
    // 			{
    // 				CSegSet repSegs;
    // 				repSegs.AddSeg(iterRep->pAlign->m_iFrom, iterRep->pAlign->m_iTo);
    //
    // 				if (string::npos != iterRep->ulGapsIdx)
    // 					repSegs.Clip(vecGaps[iterRep->ulGapsIdx]);
    //
    // 				int iRepLength = (int)repSegs.GetTotalResidues();
    // 				double dRepLength = (double)iRepLength;
    //
    // 				// -- find gaps. any gap >= half of the domain model length are considered a gap and excluded from overlapping
    //
    // 				CSegSet olSegs(thisSegs);
    // 				olSegs.Cross(repSegs);
    //
    // 				int iOverlapLen = (int)olSegs.GetTotalResidues();
    // 				double dOverlapLen = (double)iOverlapLen;
    //
    // 				if (dOverlapLen > 0)
    // 				{
    //
    // 					if (dOverlapLen / dThisLength > 0.5 || (dOverlapLen / dRepLength > 0.5 && iter->pCdInfo->m_iClusterId == iterRep->pCdInfo->m_iClusterId))	//mutually overlap > 0.5
    // 					//if (dOverlapLen / dThisLength + dOverlapLen / dRepLength > 1.0)	//mutually overlap > 0.5
    // 					{
    // 						iter->pAlign->m_iRegionIdx = iterRep->pAlign->m_iRegionIdx;
    // 						iter->pAlign->m_bRep = false;
    //
    // 						goto labelNextNonMulti;
    // 					}
    // 				}
    // 				s_overlaps.Merge(olSegs);
    //
    // 			}
    //
    // 			if ((iter->pAlign->m_bRep = (((double)(s_overlaps.GetTotalResidues()) / dThisLength) < 0.5)))	//new region
    // 			{
    // 				iter->pAlign->m_iRegionIdx = (int)vecConcise.size();
    //
    // 				vecConcise.emplace_back(*iter);
    //
    // 			}
    // 		labelNextNonMulti:
    // 			srccnt.first = iter->pAlign->m_iRegionIdx;
    // 			map<int, TDomSrcCount> :: iterator iterSrcCounter = dimAccTypeCount.emplace(srccnt).first;
    //
    // 			if (iterSrcCounter->second.CountSrc(iter->pCdInfo->m_strAccession))
    // 			{
    // 				dst.m_vecStdIndice.emplace_back(iter->ulIdx);
    // 				iter->pAlign->m_bRep = true;
    // 				if (iter->pCdInfo->m_bCurated)
    // 					dst.m_vecQualifiedFeatIndice.emplace_back(iter->ulIdx);
    //
    // 			}
    //
    //
    // 			dst.m_vecSortedIndice.emplace_back(iter->ulIdx);
    // 		}
    // ```
    std_sort(&mut non_multi, &mut |a: &Facility, b: &Facility| {
        let (a, b) = (&aligns[a.align], &aligns[b.align]);
        if a.e_value < b.e_value {
            true
        } else if a.e_value > b.e_value {
            false
        } else {
            a.bit_score > b.bit_score
        }
    });
    let segs_of = |facility: &Facility, aligns: &[DomAlign], gap_sets: &[SegSet]| {
        let mut segs = SegSet::new();
        segs.add(
            i64::from(aligns[facility.align].from),
            i64::from(aligns[facility.align].to),
        );
        if let Some(gaps) = facility.gaps {
            segs.clip(&gap_sets[gaps]);
        }
        segs
    };
    let mut concise: Vec<Facility> = Vec::new();
    let mut source_counts: HashMap<i32, [i32; 7]> = HashMap::new();
    for facility in &non_multi {
        let this_segs = segs_of(facility, aligns, &gap_sets);
        let this_length = this_segs.total_residues() as f64;
        let mut overlaps = SegSet::new();
        let mut joined = false;
        for rep in &concise {
            let rep_segs = segs_of(rep, aligns, &gap_sets);
            let rep_length = rep_segs.total_residues() as f64;
            let mut overlap = this_segs.clone();
            overlap.cross(&rep_segs);
            let overlap_length = overlap.total_residues() as f64;
            if overlap_length > 0.0
                && (overlap_length / this_length > 0.5
                    || (overlap_length / rep_length > 0.5
                        && data.domains[facility.domain].cluster_id
                            == data.domains[rep.domain].cluster_id))
            {
                let region = aligns[rep.align].region;
                let align = &mut aligns[facility.align];
                align.region = region;
                align.rep = false;
                joined = true;
                break;
            }
            overlaps.merge(&overlap);
        }
        if !joined {
            let rep = overlaps.total_residues() as f64 / this_length < 0.5;
            aligns[facility.align].rep = rep;
            if rep {
                aligns[facility.align].region = concise.len() as i32;
                concise.push(*facility);
            }
        }
        let domain = &data.domains[facility.domain];
        let counts = source_counts
            .entry(aligns[facility.align].region)
            .or_insert([0; 7]);
        if count_source(counts, data.text(domain.accession)) {
            index.standard.push(facility.align);
            aligns[facility.align].rep = true;
            if domain.curated {
                index.qualified_feat.push(facility.align);
            }
        }
        index.sorted.push(facility.align);
    }
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_base.cpp:275-303
    // ```c++
    // 		// -- When done, convert setup dst.m_vecConciseIndice
    // 		size_t ttlWinners = vecConcise.size();
    // 		if (ttlWinners > 0)
    // 		{
    // 			dst.m_vecConciseIndice.reserve(ttlWinners);
    // 			for (size_t i = 0; i < ttlWinners; ++i)
    // 				dst.m_vecConciseIndice.emplace_back(vecConcise[i].ulIdx);
    // 		}
    // 		vecPrivileged.clear();
    //
    //
    // 		if (!dst.m_vecSortedIndice.empty())
    // 		{
    // 			if (dst.m_vecConciseIndice.empty())
    // 			{
    // 				size_t iEnd = vecPrivileged.size();
    // 				...
    // 			}
    //
    // 			if (dst.m_vecStdIndice.empty())
    // 			{
    // 				dst.m_vecStdIndice = dst.m_vecSortedIndice;
    // 			}
    // 		}
    // ```
    // (`vecPrivileged` is always empty.)
    index.concise = concise.iter().map(|facility| facility.align).collect();
    if !index.sorted.is_empty() && index.standard.is_empty() {
        index.standard = index.sorted.clone();
    }
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_base.cpp:305-342
    // ```c++
    // 		// -- finally calculate SDs
    // 		// -- SDs: models are sorted according to evalues, each model reserve their regions on the query sequence. that would need to change the from/to and the aligned segments. Features are
    // 		// -- then mapped to the query from these regions -- to guarantee non-redundency. Proposed by Aron.
    // 		// -- implementation: calculate a restriction segset for SD models, which will be used later to trim the sites first
    // 		if (!vecSDs.empty())
    // 		{
    //
    // 			sort(vecSDs.begin(), vecSDs.end(), sortingEV);
    // 			//vecConcise.clear();
    // 			dst.m_vecSDIndice.clear();
    // 			dst.m_vecSDIndice.reserve(vecSDs.size());
    // 			CSegSet covered_regions;
    //
    // 			for (vector <TDomAlignFacility>::const_iterator iter = vecSDs.begin(), iterEnd = vecSDs.end(); iterEnd != iter; ++iter)
    // 			{
    //
    // 				CSegSet sdSegSet;
    // 				sdSegSet.AddSeg(0, iter->pCdInfo->m_uiLength - 1);
    //
    // 				iter->pAlign->MapSegSet(sdSegSet);
    //
    // 				sdSegSet.Clip(covered_regions);
    //
    // 				if (!sdSegSet.IsEmpty())
    // 				{
    // 					iter->pAlign->m_ClipSet.Clear();
    // 					const CSegSet::TSegs& sdSegs = sdSegSet.GetSegs();
    //
    // 					for (CSegSet::TSegs::const_iterator iterSeg = sdSegs.begin(), iterSegEnd = sdSegs.end(); iterSegEnd != iterSeg; ++iterSeg)
    // 					{
    // 						iter->pAlign->m_ClipSet.AddSeg(iterSeg->ori_from, sdSegSet.GetOriTo(iterSeg));
    // 					}
    // 					covered_regions.Merge(sdSegSet);
    // 					dst.m_vecSDIndice.emplace_back(iter->ulIdx);
    // 				}
    // 			}
    //
    // 		}
    // ```
    // (`m_uiLength - 1` is `unsigned int`, converted to `SeqPos_t` by `AddSeg`: a model of length
    // 0 gives -1, an empty range.)
    if !structural.is_empty() {
        std_sort(&mut structural, &mut |a: &Facility, b: &Facility| {
            let (a, b) = (&aligns[a.align], &aligns[b.align]);
            if a.e_value < b.e_value {
                true
            } else if a.e_value > b.e_value {
                false
            } else {
                a.bit_score > b.bit_score
            }
        });
        let mut covered = SegSet::new();
        for facility in &structural {
            let mut sd_segs = SegSet::new();
            sd_segs.add(
                0,
                i64::from(data.domains[facility.domain].length.wrapping_sub(1) as i32),
            );
            map_seg_set(&mut sd_segs, &aligns[facility.align]);
            sd_segs.clip(&covered);
            if !sd_segs.is_empty() {
                let mut clip = SegSet::new();
                for seg in sd_segs.segs() {
                    clip.add(seg.ori_from, SegSet::ori_to(seg));
                }
                aligns[facility.align].clip_set = clip;
                covered.merge(&sd_segs);
                index.structural.push(facility.align);
            }
        }
    }
    index
}

/// `TSeqAlignment::MapSegSet(segset)` of a protein alignment.
pub(crate) fn map_seg_set(set: &mut SegSet, align: &DomAlign) {
    set.map_through(&align.m_starts, &align.s_starts, &align.lens);
}

/// `TSeqAlignment::GetTranslatedPosMap` of a protein alignment: the residues of a mapped set,
/// each with its position on the query and its original position.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:1252-1294
// ```c++
// READINGFRAME::TFRAMEINDEX TSeqAlignment::GetTranslatedPosMap(const CSegSet &mappedAAsegs, SeqLen_t qLen, vector<TSeg_base::TResiduePos> &rMappedAAPos) const
// {
// 	rMappedAAPos.clear();
//
// 	const CSegSet::TSegs &segs = mappedAAsegs.GetSegs();
// 	size_t ttlres = mappedAAsegs.GetTotalResidues();
//
// 	rMappedAAPos.reserve(ttlres);
// 	...
// 	// -- either protein or plus strand translation, no adjust needed.
// 	for (CSegSet::TSegs::const_iterator iterSeg = segs.begin(), iterSegEnd = segs.end(); iterSegEnd != iterSeg; ++iterSeg)
// 	{
// 		for (SeqPos_t c = iterSeg->from, inc = 0; c <= iterSeg->to; ++c, ++inc)
// 		{
// 			rMappedAAPos.emplace_back(c, iterSeg->ori_from + inc);
// 		}
// 	}
// ```
pub(crate) fn translated_pos_map(set: &SegSet) -> Vec<ResiduePos> {
    let mut positions = Vec::with_capacity(set.total_residues().max(0) as usize);
    for seg in set.segs() {
        let mut c = seg.from;
        let mut inc = 0i64;
        while c <= seg.to {
            positions.push(ResiduePos {
                curr: c as i32,
                ori: (seg.ori_from + inc) as i32,
            });
            c += 1;
            inc += 1;
        }
    }
    positions
}

/// `TDomSeqAlignIndex::CreateRecordSets`, the domain records (`rDomAligns`) of `mode`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:1918-2010
// ```c++
// void TDomSeqAlignIndex::CreateRecordSets(const vector<TDomSeqAlignment> &rAlignments, const TDomClusterIndexIfx & rDomInfo, vector< TDomSeqAlignIndex::__TCdAlignRecord> &rDomAligns, vector<TDomSeqAlignIndex::__TCdAlignRecord> &rFeatAligns, int mode) const
// {
// 	...
// 	if (TDataModes::e_rep == mode)	//all rep
// 	{
//
// 		for (size_t iidx = 0, iidxEnd = m_vecConciseIndice.size(); iidx < iidxEnd; ++iidx)
// 		{
// 			__TCdAlignRecord alignRec;
// 			alignRec.pAlign = &(rAlignments[m_vecConciseIndice[iidx]]);
// 			alignRec.pCdInfo = rDomInfo.FindCdInfo(alignRec.pAlign->m_uiPssmId);
//
//
// 			if (alignRec.pCdInfo->m_iClusterId > 0) alignRec.pClst = rDomInfo.FindClusterInfo(alignRec.pCdInfo->m_iClusterId);
// 			...
// 			rDomAligns.push_back(alignRec);
// 		}
// 	}
// 	else	//non-concise
// 	{
// 		...
// 		const vector<size_t> &nonConciseIdx = (TDataModes::e_std == mode ? m_vecStdIndice : m_vecSortedIndice);
//
// 		for (size_t iidx = 0, iidxEnd = nonConciseIdx.size(); iidx < iidxEnd; ++iidx)
// 		{
// 			__TCdAlignRecord alignRec;
// 			alignRec.pAlign = &(rAlignments[nonConciseIdx[iidx]]);
// 			alignRec.pCdInfo = rDomInfo.FindCdInfo(alignRec.pAlign->m_uiPssmId);
// 			if (alignRec.pCdInfo->m_iClusterId > 0) alignRec.pClst = rDomInfo.FindClusterInfo(alignRec.pCdInfo->m_iClusterId);
// 			...
// 			rDomAligns.push_back(alignRec);
// 		}
// 	}
// ```
// The hierarchy roots and the site records (`rFeatAligns`, `ExtractFeatAligns`) serve the site
// output (`feature_records`).
pub(crate) fn domain_records(
    index: &AlignIndex,
    aligns: &[DomAlign],
    data: &CddData,
    mode: DataMode,
) -> Vec<Record> {
    let indices = match mode {
        DataMode::Rep => &index.concise,
        DataMode::Std => &index.standard,
        DataMode::Full => &index.sorted,
    };
    indices
        .iter()
        .map(|&align| {
            let domain = data
                .find_domain(aligns[align].pssm_id)
                .expect("a classified hit has its model");
            let cluster_id = data.domains[domain].cluster_id;
            let cluster = if cluster_id > 0 {
                data.find_cluster(cluster_id)
            } else {
                None
            };
            Record {
                align,
                domain,
                cluster,
            }
        })
        .collect()
}

/// `TDomSeqAlignIndex::__TCdAlignRecord` of the site output: a hit whose model's sites are
/// printed, with the model of its hierarchy root (`pRootCdInfo`).
#[derive(Clone, Copy, Debug)]
pub(crate) struct FeatureRecord {
    pub(crate) align: usize,
    pub(crate) domain: usize,
    /// The root of the model's hierarchy, `None` when the model has none (or it is not in the
    /// data).
    pub(crate) root: Option<usize>,
}

/// `TDomSeqAlignIndex::CreateRecordSets`, the feature records (`rFeatAligns`) of `mode`: for
/// each representative hit, the qualified curated hits of its region (`ExtractFeatAligns`),
/// then the structural-domain hits.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:1918-2032
// ```c++
// 			if (0 == alignRec.pAlign->m_iRepClass && alignRec.pAlign->m_bRep)	//non-multi
// 			{
// 				ExtractFeatAligns(alignRec, rAlignments, rDomInfo, rFeatAligns);
// 				...
// 			}
// 	...
// 			if (alignRec.pAlign->m_bRep)
// 			{
// 				alignRec.idxidx = repIdx++;
// 				if (0 == alignRec.pAlign->m_iRepClass)	//non-multi
// 				{
// 					ExtractFeatAligns(alignRec, rAlignments, rDomInfo, rFeatAligns);
// 	...
// 	// -- added 9/9/2014 handling structure motifs -- attach to
// 	for (size_t iidx = 0, iidxEnd = m_vecSDIndice.size(); iidx < iidxEnd; ++iidx)
// 	{
// 		__TCdAlignRecord rec;
//
// 		rec.idx = m_vecSDIndice[iidx];
//
// 		rec.pAlign =  &(rAlignments[rec.idx]);
//
// 		rec.pCdInfo = rDomInfo.FindCdInfo(rec.pAlign->m_uiPssmId);
//
// 		rec.pClst = nullptr;
//
// 		rec.pRootCdInfo = nullptr;
//
// 		rec.idxidx = ccsBase + amendCount++;
//
// 		rFeatAligns.push_back(rec);
// 	}
// ```
// (The `idxidx` members only order the JSON output.)
pub(crate) fn feature_records(
    index: &AlignIndex,
    aligns: &[DomAlign],
    data: &CddData,
    mode: DataMode,
) -> Vec<FeatureRecord> {
    let indices = match mode {
        DataMode::Rep => &index.concise,
        DataMode::Std => &index.standard,
        DataMode::Full => &index.sorted,
    };
    let mut records = Vec::new();
    for &align in indices {
        let a = &aligns[align];
        if a.rep && a.rep_class == 0 {
            let domain = data
                .find_domain(a.pssm_id)
                .expect("a classified hit has its model");
            extract_feat_aligns(align, domain, index, aligns, data, &mut records);
        }
    }
    for &align in &index.structural {
        let domain = data
            .find_domain(aligns[align].pssm_id)
            .expect("a classified hit has its model");
        records.push(FeatureRecord {
            align,
            domain,
            root: None,
        });
    }
    records
}

/// `TDomSeqAlignIndex::ExtractFeatAligns`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:2208-2272
// ```c++
// void TDomSeqAlignIndex :: ExtractFeatAligns(const TDomSeqAlignIndex::__TCdAlignRecord &rRepRec, const vector<TDomSeqAlignment> &rAlignments, const TDomClusterIndexIfx & rDomInfo, vector<TDomSeqAlignIndex::__TCdAlignRecord> &rResult) const
// {
// 	__TCdAlignRecord rec;
// 	if (rRepRec.pCdInfo->m_iClusterId > 0)	//has cluster, do cluster match
// 	{
// 		for (size_t fiidx = 0, fiidxEnd = m_vecQualifiedFeatIndice.size(); fiidx < fiidxEnd; ++fiidx)
// 		{
// 			rec.pAlign = &(rAlignments[m_vecQualifiedFeatIndice[fiidx]]);
// 			// -- check if duplicated, one domain may have multiple features
// 			for (const auto & r : rResult)
// 			{
// 				if (r.pAlign == rec.pAlign)
// 					goto labelSkip1;
// 			}
//
// 			if (rec.pAlign->m_iRegionIdx == rRepRec.pAlign->m_iRegionIdx)	//region match
// 			{
// 				rec.pCdInfo = rDomInfo.FindCdInfo(rec.pAlign->m_uiPssmId);
// 				if (rec.pCdInfo->m_iClusterId == rRepRec.pCdInfo->m_iClusterId)	//matched cluster
// 				{
// 					rec.pClst = rRepRec.pClst;
// 					rec.idx = m_vecQualifiedFeatIndice[fiidx];
// 					//rec.idxidx = FeatIdx2Iidx(rec.idx);
// 					if (rec.pCdInfo->m_uiHierarchyRoot > 0)
// 					{
// 						if (rec.pCdInfo->m_uiHierarchyRoot == rec.pAlign->m_uiPssmId) rec.pRootCdInfo = rec.pCdInfo;
// 						else rec.pRootCdInfo = rDomInfo.FindCdInfo(rec.pCdInfo->m_uiHierarchyRoot);
// 					}
// 					rResult.push_back(rec);
// 				}
// 			}
// 		labelSkip1:
// 			;
// 		}
// 	}
// 	else	//no cluster, just match region class -- should not happen
// 	{
// 		... (the same, without the cluster comparison; `rec.pClst = rDomInfo.FindClusterInfo(...)`)
// 	}
// }
// ```
// `rec` is one object for the whole call: a model without a hierarchy root (`m_uiHierarchyRoot`
// 0) keeps the `pRootCdInfo` of the record pushed before it (`root` below, set only when the
// model has a root).
fn extract_feat_aligns(
    rep_align: usize,
    rep_domain: usize,
    index: &AlignIndex,
    aligns: &[DomAlign],
    data: &CddData,
    result: &mut Vec<FeatureRecord>,
) {
    let rep_cluster = data.domains[rep_domain].cluster_id;
    let mut root: Option<usize> = None;
    for &candidate in &index.qualified_feat {
        // One model may have several hits; each hit is listed once.
        if result.iter().any(|record| record.align == candidate) {
            continue;
        }
        if aligns[candidate].region != aligns[rep_align].region {
            continue;
        }
        let domain = data
            .find_domain(aligns[candidate].pssm_id)
            .expect("a classified hit has its model");
        if rep_cluster > 0 && data.domains[domain].cluster_id != rep_cluster {
            continue;
        }
        let hierarchy_root = data.domains[domain].hierarchy_root;
        if hierarchy_root > 0 {
            root = if hierarchy_root == aligns[candidate].pssm_id {
                Some(domain)
            } else {
                data.find_domain(hierarchy_root)
            };
        }
        result.push(FeatureRecord {
            align: candidate,
            domain,
            root,
        });
    }
}

#[cfg(test)]
#[path = "domains_tests.rs"]
mod tests;
