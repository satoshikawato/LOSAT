//! The domain hits of one protein query: the alignments of an archive turned into rpsbproc's
//! `TDomSeqAlignment`s, their classification into regions, representative, standard and full
//! sets (`CCdAlignProcessor::Calculate`), and the records that the DOMAINS output prints
//! (`TDomSeqAlignIndex::CreateRecordSets`).

use std::collections::HashMap;

use super::archive::{Align, DenseSeg, StdLoc, StdPair};
use super::data::{dom_acc_type, CddData, TOTALSIGS};
use super::frames::{minus_pr_to_plus_na, plus_pr_to_plus_na, plus_rf_idx_to_minus, RF_SIZE};
use super::segset::SegSet;
use super::sites::ResiduePos;
use crate::utils::std_sort::{try_std_sort, OutOfRange};

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
// (`m_iFrom` and `m_iTo` are `i64`: the XML input's `GetQuery_from() - COORDSBASE` of
// `-2147483648` stores `INT_MIN - 1` there, the value of NCBI's overflowing `int` arithmetic,
// PD-LOSAT-NCBI-DEFECTS exception 11; an archive's are `int`s, converted from `TSeqPos`.)
#[derive(Clone, Debug)]
pub(crate) struct DomAlign {
    pub(crate) pssm_id: u32,
    pub(crate) from: i64,
    pub(crate) to: i64,
    pub(crate) m_starts: Vec<i32>,
    pub(crate) s_starts: SubjectStarts,
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
    /// `m_eAlignType == ePr2Na`: a hit of a translated query (a `Std-seg` whose query locations
    /// carry a strand); its `m_starts` are positions on a "pseudo-protein".
    pub(crate) pr2na: bool,
    /// `m_ReadingFrame` (`TFRAMEINDEX`, an `unsigned int`): the plus-strand frame index of
    /// `from` in its two low bits, and for a minus-strand hit the aligned length in codons
    /// shifted left by 2.
    pub(crate) reading_frame: u32,
}

impl DomAlign {
    /// A hit without segments (`TSeqAlignment::Reset` and the start of `ParseAlignSegs`).
    fn empty(align: &Align, pssm_id: u32) -> Self {
        Self::new(
            pssm_id,
            i64::from(i32::MAX),
            0,
            align.e_value.unwrap_or(0.0),
            align.bit_score.unwrap_or(0.0),
        )
    }

    /// A hit of the range `from`..`to` without segments, of a protein query (the alignment
    /// type and the reading frame of `TDomSeqAlignment`'s constructor).
    pub(crate) fn new(pssm_id: u32, from: i64, to: i64, e_value: f64, bit_score: f64) -> Self {
        DomAlign {
            pssm_id,
            from,
            to,
            m_starts: Vec::new(),
            s_starts: SubjectStarts::default(),
            lens: Vec::new(),
            e_value,
            bit_score,
            region: 0,
            n_missing: 0.0,
            c_missing: 0.0,
            spec_qualified: false,
            rep_class: 0,
            rep: false,
            clip_set: SegSet::new(),
            pr2na: false,
            reading_frame: 0,
        }
    }

    /// The aligned length in codons of a minus-strand hit (`m_ReadingFrame >> 2`), 0 for the
    /// others.
    pub(crate) fn aligned_len(&self) -> u32 {
        self.reading_frame >> 2
    }

    /// The plus-strand frame index (`m_ReadingFrame & READINGFRAME::RF_SIZE`).
    pub(crate) fn plus_rf_idx(&self) -> u32 {
        self.reading_frame & RF_SIZE
    }
}

/// `m_vecSStarts` (`std::vector<SeqPos_t>`): the starts of the segments on the model. NCBI keeps
/// `int`s; LOSAT keeps `i32`s, 4 bytes a segment, until a start past `int`'s range comes, and
/// `i64`s from then on: `ParseAlignStrings` of the XML input moves its `int` `hf` past the range
/// with `++hf`, `--hf` (and `- COORDSBASE` of an `Hsp_hit-from` of `-2147483648`), which LOSAT
/// computes as their values (PD-LOSAT-NCBI-DEFECTS exception 11). An archive's starts are `int`s
/// (converted from `TSeqPos`), so they stay narrow, and so does every start of a real document.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:1121
// ```c++
// 	std::vector<SeqPos_t> m_vecSStarts;
// ```
// (No equality of two of them: the storage is not part of the value. The tests compare the values
// with a list, `PartialEq<Vec<i64>>`.)
#[derive(Clone, Debug)]
pub(crate) enum SubjectStarts {
    Narrow(Vec<i32>),
    Wide(Vec<i64>),
}

impl Default for SubjectStarts {
    fn default() -> Self {
        SubjectStarts::Narrow(Vec::new())
    }
}

impl SubjectStarts {
    /// `push_back`.
    pub(crate) fn push(&mut self, start: i64) {
        match self {
            SubjectStarts::Narrow(starts) => match i32::try_from(start) {
                Ok(start) => starts.push(start),
                Err(_) => {
                    let mut wide: Vec<i64> = starts.iter().map(|&start| i64::from(start)).collect();
                    wide.push(start);
                    *self = SubjectStarts::Wide(wide);
                }
            },
            SubjectStarts::Wide(starts) => starts.push(start),
        }
    }

    pub(crate) fn len(&self) -> usize {
        match self {
            SubjectStarts::Narrow(starts) => starts.len(),
            SubjectStarts::Wide(starts) => starts.len(),
        }
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// `m_vecSStarts[i]` (`i` below `len()`).
    pub(crate) fn get(&self, i: usize) -> i64 {
        match self {
            SubjectStarts::Narrow(starts) => i64::from(starts[i]),
            SubjectStarts::Wide(starts) => starts[i],
        }
    }
}

/// The starts' values against a list, whatever the storage (for the tests).
#[cfg(test)]
impl PartialEq<Vec<i64>> for SubjectStarts {
    fn eq(&self, other: &Vec<i64>) -> bool {
        self.len() == other.len()
            && other
                .iter()
                .enumerate()
                .all(|(i, &start)| self.get(i) == start)
    }
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
            align.s_starts.push(i64::from(s));
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
    let mut dom = DomAlign::empty(align, pssm_id);
    clean_alignment(seg, &mut dom);
    if let (Some(&first), Some(&last), Some(&last_len)) =
        (dom.m_starts.first(), dom.m_starts.last(), dom.lens.last())
    {
        dom.from = i64::from(first);
        dom.to = i64::from((last as u32).wrapping_add(last_len).wrapping_sub(1) as i32);
    }
    dom
}

/// `ParseAlignSegs` of the `Std-seg`s of a translated search and `ParseAlignScores`: the
/// alignment, without its PSSM-ID (`ProcessCDAlign` sets it), or the refusal of what NCBI
/// reads outside its vectors.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_blast.cpp:44-133
// ```c++
// 	case CSeq_align::TSegs::e_Std:
// 		{
//
// 			const CSeq_align::TSegs::TStd &stdsegs = c_segs.GetStd();
//
// 			if (!stdsegs.empty())
// 			{
// 				CSeq_align::TSegs::TStd::const_iterator iterStdSeg = stdsegs.begin(), iterStdSegEnd = stdsegs.end();
//
// 				const CStd_seg::TIds & idsTbl = (*iterStdSeg)->GetIds();	//vector< CRef< CSeq_id > >
// 				CStd_seg::TIds::const_iterator iter_hit_id = idsTbl.begin();
// 				++iter_hit_id;
// 				hit_id = SerialClone(**iter_hit_id);
//
// 				for (; iterStdSegEnd != iterStdSeg; ++iterStdSeg)
// 				{
// 					const CStd_seg::TLoc &seqlocs = (*iterStdSeg)->GetLoc();	//vector< CRef< CSeq_loc > >
//
// 					if (seqlocs[0]->IsNull() || seqlocs[1]->IsNull()) continue;	//unmapped. ignore
//
// 					CSeq_loc::E_Choice qloc_type = seqlocs[0]->Which(), hloc_type = seqlocs[1]->Which();
// 					if (CSeq_loc::e_Empty == qloc_type || CSeq_loc::e_Empty == hloc_type) continue;
//
// 					bool error = false;
// 					// -- query
// 					switch (qloc_type)
// 					{
// 					case CSeq_loc::e_Int:	//interval
// 						{
// 							const CSeq_interval &intv = seqlocs[0]->GetInt();
// 							SeqPos_t f = intv.GetFrom(), t = intv.GetTo();
// 							if (f < dst.m_iFrom) dst.m_iFrom = f;
// 							if (t > dst.m_iTo) dst.m_iTo = t;
// 							dst.m_vecMStarts.push_back(f);
//
// 							if (intv.CanGetStrand())
// 								qstrand = intv.GetStrand();
// 						}
// 						break;
// 					case CSeq_loc::e_Pnt:	//seq point
// 						{
// 							const CSeq_point &pnt = seqlocs[0]->GetPnt();
// 							SeqPos_t p = pnt.GetPoint();
// 							if (p < dst.m_iFrom) dst.m_iFrom = p;
// 							if (p > dst.m_iTo) dst.m_iTo = p;
// 							dst.m_vecMStarts.push_back(p);
//
// 							if (pnt.CanGetStrand())
// 								qstrand = pnt.GetStrand();
//
// 						}
// 						break;
// 					default:
// 						error = true;
// 						break;
// 					}
// 					if (!error)
// 					{
//
// 						// -- target, assume always protein
// 						switch (hloc_type)
// 						{
// 						case CSeq_loc::e_Int:	//interval
// 							{
// 								const CSeq_interval &intv = seqlocs[1]->GetInt();
// 								SeqPos_t f = intv.GetFrom(), t = intv.GetTo();
//
// 								dst.m_vecSStarts.push_back(f);
// 								dst.m_vecLens.push_back(t - f + 1);
// 							}
// 							break;
// 						case CSeq_loc::e_Pnt:	//seq point
//
// 							{
// 								const CSeq_point &pnt = seqlocs[1]->GetPnt();
// 								SeqPos_t p = pnt.GetPoint();
//
// 								dst.m_vecSStarts.push_back(p);
// 								dst.m_vecLens.push_back(1);
//
// 							}
// 							break;
// 						default:
// 							error = true;
// 							break;
// 						}
// 					}
// 				}
//
// 			}
// 		}
// 		break;
// ```
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_blast.cpp:166-209
// ```c++
// 	if (!hit_id.IsNull() && eNa_strand_other != qstrand)	//na and pr
// 	{
// 		dst.m_eAlignType = TSeqAlignment::ePr2Na;
// 		dst.m_ReadingFrame = (dst.m_iFrom % READINGFRAME::RF_SIZE);	//always be the plus side rf
//
// 		if (eNa_strand_minus == qstrand || eNa_strand_both_rev == qstrand)
// 		{
// 			dst.m_bIsMinus = true;
// 			SeqLen_t alignedLen = (dst.m_iTo + 1) / READINGFRAME::RF_SIZE;
// 			dst.m_ReadingFrame |= (alignedLen << 2);	//keep this for align computation
// 			...
// 			vector<SeqPos_t> :: iterator iterMS = dst.m_vecMStarts.begin(), iterMSEnd = dst.m_vecMStarts.end();
// 			vector<SeqLen_t> :: const_iterator iterL = dst.m_vecLens.begin();
// 			while (iterMS != iterMSEnd)
// 			{
// 				*iterMS = alignedLen - (*iterMS / READINGFRAME::RF_SIZE + *iterL);
// 				++iterMS;
// 				++iterL;
// 			}
// 		}
// 		else
// 		{
// 			for (vector<SeqPos_t> :: iterator iterMS = dst.m_vecMStarts.begin(), iterMSEnd = dst.m_vecMStarts.end(); iterMS != iterMSEnd; ++iterMS)
// 				*iterMS /= READINGFRAME::RF_SIZE;
// 		}
//
// 	}
// ```
// The types decide the arithmetic: `GetFrom`/`GetTo`/`GetPoint` are `TSeqPos` (unsigned)
// stored in `SeqPos_t` (`int`, wraps in NCBI's build); `t - f + 1` is `int` arithmetic, which
// overflows for an interval of 2^31 positions or more (undefined): computed as its value
// (PD-LOSAT-NCBI-DEFECTS exception 11), which `push_back` converts to the `unsigned int` of
// `m_vecLens` (modulo 2^32); `m_iFrom % RF_SIZE`, `*iterMS / RF_SIZE` and the minus-strand sum
// are `unsigned int` arithmetic (wraps); `m_iTo + 1` is `int` (`m_iTo` is at least 0), which
// overflows for a hit that ends at 2^31 - 1: computed as its value (exception 11), which the
// division by the `unsigned int` `RF_SIZE` converts modulo 2^32. A segment whose subject location is
// neither `int` nor `pnt` leaves its query start without a subject start and length: the
// vectors differ in length, which the plus-strand code and the classification read within
// bounds (LOSAT keeps the three vectors as NCBI does), but the minus-strand loop above reads
// `m_vecLens` past its end (refused). `ids` with fewer than 2 ids in the first segment is read
// past its end (refused, by the caller). `seqlocs[0]->IsNull() || seqlocs[1]->IsNull()` reads
// `seqlocs[1]` only when the first location is not `null`: a segment without locations, or
// whose only location is not `null`, is read past its end (refused), and a segment whose first
// location is `null` is skipped whatever follows (PS16c review oracle `nullone`). The strand
// that counts is the last one set; `unknown`, `plus` and `both` make a plus-strand hit, `minus`
// and `both-rev` a minus-strand one, `other` (or none) an alignment of two proteins.
pub(crate) fn parse_std_align(
    segs: &[StdPair],
    align: &Align,
    pssm_id: u32,
) -> Result<DomAlign, String> {
    const STRAND_MINUS: i32 = 2;
    const STRAND_BOTH_REV: i32 = 4;
    const STRAND_OTHER: i32 = 255;
    let mut dom = DomAlign::empty(align, pssm_id);
    let mut qstrand = STRAND_OTHER;
    for &StdPair { query, hit } in segs {
        if query == StdLoc::Missing {
            return Err(
                "a Std-seg without locations (NCBI reads past them) is not supported by LOSAT's rpsbproc".to_string(),
            );
        }
        if query == StdLoc::Null {
            continue;
        }
        if hit == StdLoc::Missing {
            return Err(
                "a Std-seg whose only location is not null (NCBI reads past it) is not supported by LOSAT's rpsbproc".to_string(),
            );
        }
        if hit == StdLoc::Null {
            continue;
        }
        if query == StdLoc::Empty || hit == StdLoc::Empty {
            continue;
        }
        let (start, end, strand) = match query {
            StdLoc::Int { from, to, strand } => (from as i32, to as i32, strand),
            StdLoc::Pnt { point, strand } => (point as i32, point as i32, strand),
            _ => continue,
        };
        if i64::from(start) < dom.from {
            dom.from = i64::from(start);
        }
        if i64::from(end) > dom.to {
            dom.to = i64::from(end);
        }
        dom.m_starts.push(start);
        if let Some(strand) = strand {
            qstrand = strand;
        }
        match hit {
            StdLoc::Int { from, to, .. } => {
                let (f, t) = (i64::from(from as i32), i64::from(to as i32));
                dom.s_starts.push(f);
                dom.lens.push((t - f + 1) as u32);
            }
            StdLoc::Pnt { point, .. } => {
                dom.s_starts.push(i64::from(point as i32));
                dom.lens.push(1);
            }
            _ => {}
        }
    }
    if qstrand != STRAND_OTHER {
        dom.pr2na = true;
        dom.reading_frame = (dom.from as u32) % RF_SIZE;
        if qstrand == STRAND_MINUS || qstrand == STRAND_BOTH_REV {
            let aligned_len = (dom.to + 1) as u32 / RF_SIZE;
            dom.reading_frame |= aligned_len << 2;
            if dom.m_starts.len() > dom.lens.len() {
                return Err(
                    "a minus-strand Std-seg with a subject location other than int or pnt (NCBI reads past its segment lengths) is not supported by LOSAT's rpsbproc".to_string(),
                );
            }
            for (start, &len) in dom.m_starts.iter_mut().zip(&dom.lens) {
                *start =
                    aligned_len.wrapping_sub((*start as u32 / RF_SIZE).wrapping_add(len)) as i32;
            }
        } else {
            for start in &mut dom.m_starts {
                *start = (*start as u32 / RF_SIZE) as i32;
            }
        }
    }
    Ok(dom)
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

/// `CCdAlignProcessor::Calculate` of the alignments `selected` (all of a protein query's, or
/// those of one reading frame of a nucleotide query's): their scores against their models, the
/// regions, and the sets of `index`. Structural domains (`sd`) are set apart and take no region
/// (their clipping, `m_vecSDIndice` and `m_ClipSet`, serves the motif output).
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
// (`m_uiLength - (sstart + len)` is `unsigned int` arithmetic stored in an `int`, of the start's
// value modulo 2^32 where the XML input took it past `int` (`SubjectStarts`); `iNMissing` is
// that value (PD-LOSAT-NCBI-DEFECTS exception 11); the ratios
// divide doubles, so a model of length 0 gives an infinity or a NaN as NCBI's does. A NaN E-value
// or bit score makes `TDomAlignOrderByEValue` no strict weak order (`order_by_e_value`): libstdc++'s
// sort, which `try_std_sort` follows step by step, is then still defined while it reads inside the
// vector, and LOSAT gives its order. The comparison stays irreflexive and asymmetric with a NaN, so
// the median-of-three and partition sentinels of `__unguarded_partition` hold; the walk of the
// final insertion sort past the first 16 has a sentinel left of it in every part right of a
// partition, but not in the leftmost part when the depth limit hands it to the heap sort (no input
// with rpsbproc's comparison was found that reaches it: PS16d delta review S3-1). Where a loop would
// read outside the vector (`OutOfRange`), NCBI's result is its undefined read: refused, `Err` with
// the refusal. PS16d review S3-3; the permutations against libstdc++ 13 in
// `utils::std_sort` tests.)
pub(crate) fn calculate(
    aligns: &mut [DomAlign],
    selected: &[usize],
    data: &CddData,
) -> Result<AlignIndex, String> {
    let mut index = AlignIndex::default();
    let mut gap_sets: Vec<SegSet> = Vec::new();
    let mut non_multi: Vec<Facility> = Vec::new();
    let mut structural: Vec<Facility> = Vec::new();
    for &i in selected {
        let align = &mut aligns[i];
        if align.s_starts.is_empty() {
            continue;
        }
        let Some(domain_index) = data.find_domain(align.pssm_id) else {
            continue;
        };
        let domain = &data.domains[domain_index];
        let length = domain.length;
        let n_missing = align.s_starts.get(0);
        let c_missing = length.wrapping_sub(
            (align.s_starts.get(align.s_starts.len() - 1) as u32)
                .wrapping_add(*align.lens.last().expect("not empty")),
        ) as i32;
        align.n_missing = n_missing as f64 / f64::from(length);
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
    try_std_sort(&mut non_multi, &mut |a: &Facility, b: &Facility| {
        order_by_e_value(&aligns[a.align], &aligns[b.align])
    })
    .map_err(sort_out_of_range)?;
    let segs_of = |facility: &Facility, aligns: &[DomAlign], gap_sets: &[SegSet]| {
        let mut segs = SegSet::new();
        segs.add(aligns[facility.align].from, aligns[facility.align].to);
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
        try_std_sort(&mut structural, &mut |a: &Facility, b: &Facility| {
            order_by_e_value(&aligns[a.align], &aligns[b.align])
        })
        .map_err(sort_out_of_range)?;
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
    Ok(index)
}

/// `TDomAlignOrderByEValue`: the E-value ascending, then the bit score descending (with C++'s
/// comparisons of doubles: a NaN is neither less nor greater).
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
fn order_by_e_value(a: &DomAlign, b: &DomAlign) -> bool {
    if a.e_value < b.e_value {
        true
    } else if a.e_value > b.e_value {
        false
    } else {
        a.bit_score > b.bit_score
    }
}

/// The refusal of a sort of hits that would read outside NCBI's vector (`try_std_sort`).
fn sort_out_of_range(_: OutOfRange) -> String {
    "hits whose NaN E-values or bit scores make NCBI's std::sort of them read outside its vector (undefined) are not supported by LOSAT's rpsbproc".to_string()
}

/// `CCdAlignProcessor::SortReadingFrames`: the indices of the alignments of each reading frame
/// (`GetRFIdx` of the query's length), or NCBI's exception for a query of length 0.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_base.cpp:20-45
// ```c++
// SeqLen_t CCdAlignProcessor::SortReadingFrames(const vector<TDomSeqAlignment> &rSrc, vector<size_t> rfIndice[READINGFRAME::TOTAL_RFS], SeqLen_t qlen)
// {
// 	...
// 	if (qlen > 0)	//calculate
// 		for (size_t idx = 0, len = rSrc.size(); idx < len; ++idx)
// 		{
// 			...
// 			rfIndice[rSrc[idx].GetRFIdx(qlen)].emplace_back(idx);
// 		}
// 	else	//should not happen
// 		THROW_SIMPLE("Cannot sort reading frames - query sequence length unknown");
// 	return qlen;
// }
// ```
pub(crate) fn sort_reading_frames(
    aligns: &[DomAlign],
    query_len: u32,
) -> Result<[Vec<usize>; 6], &'static str> {
    if query_len == 0 {
        return Err("Cannot sort reading frames - query sequence length unknown");
    }
    let mut frames: [Vec<usize>; 6] = Default::default();
    for (i, align) in aligns.iter().enumerate() {
        frames[rf_index(align, query_len) as usize].push(i);
    }
    Ok(frames)
}

/// `TSeqAlignment::GetRFIdx(seqLen)`: the reading frame index (0 to 5) of an alignment.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:1342-1351
// ```c++
// READINGFRAME::TFRAMEINDEX TSeqAlignment::GetRFIdx(SeqLen_t seqLen) const
// {
// 	if (ePr2Na != m_eAlignType) return 0;
// 	SeqLen_t alignedLen = m_ReadingFrame >> 2;	//m_ReadingFrame contains the positive side reading frame, always.
// 	READINGFRAME::TFRAMEINDEX rfidx = m_ReadingFrame & READINGFRAME::RF_SIZE;
//
// 	if (alignedLen > 0)	//negative -- minus strand
// 		rfidx = READINGFRAME::PlusRFIdx2MinusRFIdx(rfidx, seqLen);
// 	return rfidx;
// }
// ```
pub(crate) fn rf_index(align: &DomAlign, query_len: u32) -> u32 {
    if !align.pr2na {
        return 0;
    }
    let rfidx = align.plus_rf_idx();
    if align.aligned_len() > 0 {
        plus_rf_idx_to_minus(rfidx, query_len)
    } else {
        rfidx
    }
}

/// `TSeqAlignment::MapSegSet(segset)` of a protein alignment.
pub(crate) fn map_seg_set(set: &mut SegSet, align: &DomAlign) {
    set.map_through(&align.m_starts, &align.s_starts, &align.lens);
}

/// `TSeqAlignment::GetTranslatedPosMap`: the residues of a mapped set, each with its position
/// on the query (on the translation of its frame for a hit of a translated query) and its
/// original position; `Err` with the text of NCBI's exception when a minus-strand hit reaches
/// past the query's length.
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
// 	SeqLen_t aaLen = qLen;	//assume it's aa
// 	READINGFRAME::TFRAMEINDEX rfidx = m_ReadingFrame & READINGFRAME::RF_SIZE;
//
// 	if (ePr2Na == m_eAlignType)
// 	{
// 		SeqLen_t alignedLen = m_ReadingFrame >> 2;
// 		aaLen = (qLen - rfidx) / READINGFRAME::RF_SIZE;
// 		if (alignedLen > 0)	//minus strand. Plus strand will have 0 == alignedLen
// 		{
// 			if (alignedLen > aaLen)	//error
// 				THROW_SIMPLE("Invalid protein length " << aaLen << ": shorter than aligned range " << alignedLen);
// 			rfidx = READINGFRAME::PlusRFIdx2MinusRFIdx(rfidx, qLen);
// 			SeqLen_t offset = aaLen - alignedLen;
//
// 			for (CSegSet::TSegs::const_iterator iterSeg = segs.begin(), iterSegEnd = segs.end(); iterSegEnd != iterSeg; ++iterSeg)
// 			{
// 				for (SeqPos_t c = iterSeg->from, inc = 0; c <= iterSeg->to; ++c, ++inc)
// 				{
// 					rMappedAAPos.emplace_back(c + offset, iterSeg->ori_from + inc);
// 				}
// 			}
// 			goto labelReturn;
// 		}
// 	}
// 	// -- either protein or plus strand translation, no adjust needed.
// 	for (CSegSet::TSegs::const_iterator iterSeg = segs.begin(), iterSegEnd = segs.end(); iterSegEnd != iterSeg; ++iterSeg)
// 	{
// 		for (SeqPos_t c = iterSeg->from, inc = 0; c <= iterSeg->to; ++c, ++inc)
// 		{
// 			rMappedAAPos.emplace_back(c, iterSeg->ori_from + inc);
// 		}
// 	}
// ```
// (`qLen - rfidx` and `c + offset` are `unsigned int` arithmetic, which wraps; `c + offset`
// is stored in an `int`. `++c` and `ori_from + inc` are `int` arithmetic, which overflows past
// `INT_MAX` (undefined; a mapped range at the limit, made by hand): computed as their values,
// which the plus strand stores as they are (PD-LOSAT-NCBI-DEFECTS exception 11); a position
// outside the translation is then read by the out-of-range rules of its reader.)
pub(crate) fn translated_pos_map(
    set: &SegSet,
    align: &DomAlign,
    query_len: u32,
) -> Result<Vec<ResiduePos>, String> {
    // `Some(offset)` on the minus strand, whose `c + offset` is an `unsigned int`.
    let mut minus_offset = None;
    if align.pr2na {
        let aligned_len = align.aligned_len();
        let aa_len = query_len.wrapping_sub(align.plus_rf_idx()) / RF_SIZE;
        if aligned_len > 0 {
            if aligned_len > aa_len {
                return Err(format!(
                    "Invalid protein length {aa_len}: shorter than aligned range {aligned_len}"
                ));
            }
            minus_offset = Some(aa_len - aligned_len);
        }
    }
    let mut positions = Vec::with_capacity(set.total_residues().max(0) as usize);
    for seg in set.segs() {
        let mut c = seg.from;
        let mut inc = 0i64;
        while c <= seg.to {
            let curr = match minus_offset {
                Some(offset) => i64::from((c as u32).wrapping_add(offset) as i32),
                None => c,
            };
            positions.push(ResiduePos {
                curr,
                ori: seg.ori_from + inc,
            });
            c += 1;
            inc += 1;
        }
    }
    Ok(positions)
}

/// The two ends of a mapped set after `TSeqAlignment::Pr2NaConvert` (`GetLeft`, `GetRight`):
/// positions on the query's plus strand for a hit of a translated query; the set's own ends
/// for the others (`None`: an empty set). A position that NCBI's `int` arithmetic took past
/// `int`'s range holds its value, and `pr - 1` of `INT_MIN` is computed as its value
/// (`frames::minus_pr_to_plus_na`; PD-LOSAT-NCBI-DEFECTS exception 11): the conversions after
/// them are `unsigned int`'s, modulo 2^32.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:1137-1177
// ```c++
// void TSeqAlignment::Pr2NaConvert(CSegSet &segset) const
// {
// 	if (ePr2Na == m_eAlignType)
// 	{
// 		...
// 		SeqLen_t alnLen = (m_ReadingFrame >> 2);
//
// 		READINGFRAME::TFRAMEINDEX rfidx = m_ReadingFrame & READINGFRAME::RF_SIZE;	//reading frame at positive strand
// 		segset.m_iFactor = READINGFRAME::RF_SIZE;
// 		if (alnLen > 0)	//is minus
// 		{
// 			// -- alnLen is aligned pr len. Convert to na len
// 			alnLen = alnLen * READINGFRAME::RF_SIZE + rfidx;
// 			...
// 			for (CSegSet::TSegs::iterator iterSeg = segset.m_lstContainer.begin(); iterSeg != segset.m_lstContainer.end(); ++iterSeg)
// 			{
// 				SeqPos_t newFrom = READINGFRAME::MinusPr2PlusNA(iterSeg->from, READINGFRAME::RF_SIZE, alnLen);
// 				iterSeg->from = READINGFRAME::MinusPr2PlusNA(iterSeg->to, READINGFRAME::RF_SIZE, alnLen) - READINGFRAME::RF_SIZE + 1;
// 				iterSeg->to = newFrom;
//
// 				// -- reverse ori_as well
// 				iterSeg->ori_from += (iterSeg->to - iterSeg->from + 1) / segset.m_iFactor - 1;
//
// 			}
// 			segset.m_iFactor = -segset.m_iFactor;	//change direction
// 			segset.m_lstContainer.reverse();
// 		}
// 		else	//plus strand
// 		{
// 			for (CSegSet::TSegs::iterator iter = segset.m_lstContainer.begin(); iter != segset.m_lstContainer.end(); ++iter)
// 			{
// 				iter->from = READINGFRAME::PlusPr2PlusNA(iter->from, rfidx, 0);
// 				iter->to = READINGFRAME::PlusPr2PlusNA(iter->to, rfidx, 0) + READINGFRAME::RF_SIZE - 1;
//
// 			}
// 		}
// 	}
// }
// ```
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/impl/segset_impl.hpp:281-290
// ```c++
// SeqPos_t CSegSetTmpl<SEGTYPE>::GetLeft(void) const
// {
// 	return IsEmpty() ? -1 : m_lstContainer.begin()->from;
// }
//
// template <typename SEGTYPE>
// SeqPos_t CSegSetTmpl<SEGTYPE>::GetRight(void) const
// {
// 	return IsEmpty() ? -1 : m_lstContainer.rbegin()->to;
// }
// ```
// Only the ends are printed (the motif lines), so only the first and the last ranges are
// converted; `ori_from` is not read after the conversion. The minus strand reverses the list:
// the left end comes from the last range's `to`, the right end from the first range's `from`.
// The minus strand's `ori_from += (to - from + 1) / m_iFactor - 1` is not computed: its value
// is never read, so no printed number depends on it (PS16c review S3-4).
pub(crate) fn pr2na_ends(set: &SegSet, align: &DomAlign) -> Option<(i64, i64)> {
    let (first, last) = (set.segs().first()?, set.segs().last()?);
    if !align.pr2na {
        return Some((first.from, last.to));
    }
    let rfidx = align.plus_rf_idx();
    let aligned_len = align.aligned_len();
    if aligned_len > 0 {
        let na_len = aligned_len.wrapping_mul(RF_SIZE).wrapping_add(rfidx);
        let left = (minus_pr_to_plus_na(last.to, RF_SIZE, na_len) as u32)
            .wrapping_sub(RF_SIZE)
            .wrapping_add(1) as i32;
        let right = minus_pr_to_plus_na(first.from, RF_SIZE, na_len);
        Some((i64::from(left), i64::from(right)))
    } else {
        let left = plus_pr_to_plus_na(first.from, rfidx);
        let right = (plus_pr_to_plus_na(last.to, rfidx) as u32)
            .wrapping_add(RF_SIZE)
            .wrapping_sub(1) as i32;
        Some((i64::from(left), i64::from(right)))
    }
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
