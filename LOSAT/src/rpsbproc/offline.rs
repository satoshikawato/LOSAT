//! The query of one `Iteration` of a BlastOutput XML document (rpsbproc's deprecated XML input):
//! `COfflCdAlignProcessor::ParseBlastOutput` with `_TOfflDomQuery`, which takes the query's
//! label, title and length from the Iteration, its hits from the HSPs whose E-value passes the
//! cutoff, the segments of each hit and the query's residues from the HSP's alignment strings,
//! and then classifies the hits (`Calculate`, per reading frame for a nucleotide query).
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_cd_align_proc.hpp:46-57
//! ```c++
//! struct _TOfflDomQuery: public TDomQuery
//! {
//! 	_TOfflDomQuery(void):
//! 		TDomQuery(),
//! 		m_strMessage(), m_dimTranslated() {};
//! 	void ParseAlignStrings(const std::string &qseq, const std::string &hseq, SeqPos_t qf, SeqPos_t hf, SeqPos_t start, SeqLen_t n, SeqLen_t qofs, TDomSeqAlignment &dst, READINGFRAME::TFRAMEINDEX rfidx = 0);
//!
//! 	std::string m_strMessage;
//! 	std::string m_dimTranslated[READINGFRAME::TOTAL_RFS];
//! 	// -- internals
//! 	void NACommit(void);
//! 	void PRCommit(void);
//! };
//! ```
//!
//! The HSPs are processed as the reader reads them (`OfflineQuery::hsp`), in NCBI's order, so a
//! query keeps its hits and its strings, not the text of its alignments. What `ParseBlastOutput`
//! throws is kept (`OfflineQuery::failure`) and thrown when NCBI would, after the document is read
//! and the queries before it printed (`finish`); so are LOSAT's refusals of what NCBI does not
//! define.

use super::data::CddData;
use super::domains::{calculate, sort_reading_frames, DomAlign};
use super::frames::{
    minus_id_to_plus_rf_idx, plus_id_to_plus_rf_idx, plus_na_to_minus_pr, plus_na_to_plus_pr,
    plus_rf_idx_to_minus, RF_SIZE,
};
use super::query::{unassigned_member, DomQuery, ProcessError};

/// `COfflCdAlignProcessor::COORDSBASE`.
const COORDSBASE: i64 = 1;

/// The memory from which LOSAT refuses NCBI's strings of a query (1 GiB at once, PS16b's limit
/// of DECISIONS 37 and 40: below it NCBI's allocation succeeds on the machines that run it).
const MEMORY_LIMIT: u64 = 1 << 30;

/// One HSP as `ParseBlastOutput` reads it: the members it gets, `None` for an OPTIONAL one that
/// the document lacks (its getter throws).
#[derive(Debug, Default)]
pub(crate) struct Hsp {
    pub(crate) bit_score: f64,
    pub(crate) evalue: f64,
    pub(crate) query_from: i32,
    pub(crate) query_to: i32,
    pub(crate) hit_from: i32,
    pub(crate) query_frame: Option<i32>,
    pub(crate) identity: Option<i32>,
    pub(crate) align_len: Option<i32>,
    pub(crate) qseq: Vec<u8>,
    pub(crate) hseq: Vec<u8>,
}

fn unsupported(what: &str) -> ProcessError {
    ProcessError::Unsupported(format!(
        "{what} in the XML input is not supported by LOSAT's rpsbproc"
    ))
}

/// `_TOfflDomQuery` while `ParseBlastOutput` reads its Iteration.
#[derive(Debug, Default)]
pub(crate) struct OfflineQuery {
    accession: Vec<u8>,
    title: Vec<u8>,
    /// `m_uiSeqLen`.
    length: u32,
    /// `m_bIsNa`.
    is_na: bool,
    /// `m_dimTranslated`.
    translated: [Vec<u8>; 6],
    /// `m_strSeqData`, when it is not `translated[0]` (`seq_data_is_translated`).
    seq_data: Vec<u8>,
    /// Whether `m_strSeqData` is a copy of `m_dimTranslated[0]` as it is now (assigned after each
    /// protein hit); LOSAT copies it only when `m_dimTranslated[0]` changes otherwise.
    seq_data_is_translated: bool,
    /// `m_vecAlignments`.
    aligns: Vec<DomAlign>,
    /// What `ParseBlastOutput` threw (or LOSAT refused) at an HSP; the later HSPs are not read.
    failure: Option<ProcessError>,
}

impl OfflineQuery {
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_cd_align_proc.cpp:125-142
    // ```c++
    // void COfflCdAlignProcessor::ParseBlastOutput(_TOfflDomQuery &dst, const CIteration &blastout, vector<PssmId_t> &missed, double evcut) const
    // {
    // 	dst.m_vecAlignments.clear();
    //
    // 	if (blastout.CanGetQuery_ID())
    // 		dst.m_strCleanedInput = dst.m_strAccession = blastout.GetQuery_ID();
    //
    // 	if (blastout.CanGetQuery_def())
    // 		dst.m_strTitle = blastout.GetQuery_def();	//length of original sequence, could be na
    //
    // 	if (blastout.CanGetQuery_len())
    // 		dst.m_uiSeqLen = blastout.GetQuery_len();	//length of original sequence, could be na
    //
    //
    // 	if (blastout.CanGetMessage())
    // 		dst.m_strMessage = blastout.GetMessage();	//length of original sequence, could be na
    // ```
    // (The Iteration's members come before its hits; `m_strMessage` is never printed.)
    /// `Iteration_query-ID`.
    pub(crate) fn set_accession(&mut self, id: &[u8]) {
        self.accession = id.to_vec();
    }

    /// `Iteration_query-def`.
    pub(crate) fn set_title(&mut self, def: &[u8]) {
        self.title = def.to_vec();
    }

    /// `Iteration_query-len`, an `int` stored in the `unsigned int` `m_uiSeqLen`.
    pub(crate) fn set_length(&mut self, len: i32) {
        self.length = len as u32;
    }

    /// Makes `seq_data` hold `m_strSeqData` before `m_dimTranslated[0]` changes.
    fn keep_seq_data(&mut self) {
        if self.seq_data_is_translated {
            self.seq_data = self.translated[0].clone();
            self.seq_data_is_translated = false;
        }
    }

    /// `NACommit`: a nucleotide query, the six frames' strings of `-` (the residues of each
    /// frame of `m_uiSeqLen` bases, in `unsigned int` arithmetic). A query of 0 or 1 base gives
    /// strings of more than a billion bytes, as does a query of half a billion bases: NCBI then
    /// holds 1 GiB and more, which LOSAT refuses (DECISIONS 37, 40).
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_cd_align_proc.cpp:17-41
    // ```c++
    // void _TOfflDomQuery::NACommit(void)
    // {
    // 	m_bIsNa = true;
    // 	//m_iMolType = CSeq_inst::eMol_na;
    // 	SeqLen_t aaLen = m_uiSeqLen / 3;
    //
    // 	m_dimTranslated[0].clear();
    // 	m_dimTranslated[0].append(aaLen, '-');
    //
    // 	m_dimTranslated[3] = m_dimTranslated[0];
    //
    // 	aaLen = (m_uiSeqLen - 1) / 3;
    //
    // 	m_dimTranslated[1].clear();
    // 	m_dimTranslated[1].append(aaLen, '-');
    //
    // 	m_dimTranslated[4] = m_dimTranslated[1];
    //
    // 	aaLen = (m_uiSeqLen - 2) / 3;
    //
    // 	m_dimTranslated[2].clear();
    // 	m_dimTranslated[2].append(aaLen, '-');
    // 	m_dimTranslated[5] = m_dimTranslated[2];
    // }
    // ```
    fn na_commit(&mut self) -> Result<(), ProcessError> {
        let len = self.length;
        let frames = [len / 3, len.wrapping_sub(1) / 3, len.wrapping_sub(2) / 3];
        let total: u64 = frames.iter().map(|&n| 2 * u64::from(n)).sum();
        if total >= MEMORY_LIMIT {
            return Err(unsupported(&format!(
                "a nucleotide query of length {len}, whose six frames NCBI's rpsbproc keeps as strings of {total} bytes (1 GiB or more; NCBI's result depends on the host's memory),"
            )));
        }
        self.keep_seq_data();
        self.is_na = true;
        for (frame, &n) in frames.iter().enumerate() {
            self.translated[frame].clear();
            self.translated[frame].resize(n as usize, b'-');
            self.translated[frame + 3] = self.translated[frame].clone();
        }
        Ok(())
    }

    /// `PRCommit`: a protein query, its residues `-` (NCBI then copies them to `m_strSeqData`:
    /// a query of half a billion residues or more is refused as `na_commit` refuses).
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_cd_align_proc.cpp:43-49
    // ```c++
    // void _TOfflDomQuery::PRCommit(void)
    // {
    // 	m_bIsNa = false;
    // 	//m_iMolType = CSeq_inst::eMol_aa;
    // 	m_dimTranslated[0].clear();
    // 	m_dimTranslated[0].append(m_uiSeqLen, '-');
    // }
    // ```
    fn pr_commit(&mut self) -> Result<(), ProcessError> {
        let len = self.length;
        let total = 2 * u64::from(len);
        if total >= MEMORY_LIMIT {
            return Err(unsupported(&format!(
                "a protein query of length {len}, whose residues and their copy NCBI's rpsbproc keeps as strings of {total} bytes (1 GiB or more; NCBI's result depends on the host's memory),"
            )));
        }
        self.keep_seq_data();
        self.is_na = false;
        self.translated[0].clear();
        self.translated[0].resize(len as usize, b'-');
        Ok(())
    }

    /// One HSP of a hit whose `Hit_accession` gives `pssm_id` (`(unsigned int)atoi`), read when
    /// its element ends; nothing after a failure.
    pub(crate) fn hsp(&mut self, pssm_id: u32, hsp: &Hsp, evcut: f64) {
        if self.failure.is_some() {
            return;
        }
        if let Err(error) = self.parse_hsp(pssm_id, hsp, evcut) {
            self.failure = Some(error);
        }
    }

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_cd_align_proc.cpp:143-280
    // ```c++
    // 	if (blastout.CanGetHits())
    // 	{
    // 		const CIteration::THits &hits = blastout.GetHits();	//list< CRef< CHit > >
    // 		if (!hits.empty())
    // 		{
    // 			dst.m_vecAlignments.reserve(hits.size() * 3);	//estimated
    //
    // 			TDomSeqAlignment alnVal;
    //
    // 			for (CIteration::THits::const_iterator iterHit = hits.begin(), iterHitEnd = hits.end(); iterHitEnd != iterHit; ++iterHit)
    // 			{
    // 				PssmId_t pssmid = (unsigned int)atoi((*iterHit)->GetAccession().c_str());
    //
    // 				SeqLen_t cdLen = (*iterHit)->GetLen();
    // 				const CHit::THsps &hsps = (*iterHit)->GetHsps();	//list< CRef< CHsp > > THsps
    // 				for (CHit::THsps::const_iterator iterHsp = hsps.begin(), iterHspEnd = hsps.end(); iterHspEnd != iterHsp; ++iterHsp)
    // 				{
    // 					double eval = (*iterHsp)->GetEvalue();
    // 					if (eval > evcut) continue;
    //
    // 					dst.m_vecAlignments.push_back(alnVal);
    // 					TDomSeqAlignment &currAlign = *dst.m_vecAlignments.rbegin();
    //
    // 					currAlign.m_uiPssmId = pssmid;
    // 					currAlign.m_uiAlignedLen = (*iterHsp)->GetAlign_len();
    // 					currAlign.m_iScore = (int)((*iterHsp)->GetScore() + 0.5);
    // 					currAlign.m_dEValue = eval;
    // 					currAlign.m_dBitScore = (*iterHsp)->GetBit_score();
    // 					currAlign.m_iNumIdent = (*iterHsp)->GetIdentity();
    // 					currAlign.m_dSeqIdentity = (double)currAlign.m_iNumIdent / (double)cdLen * 100.0;
    // 					//currAlign.m_ReadingFrame = (*iterHsp)->GetQuery_frame();
    //
    // 					currAlign.m_iFrom = (*iterHsp)->GetQuery_from() - COORDSBASE;
    // 					currAlign.m_iTo = (*iterHsp)->GetQuery_to() - COORDSBASE;
    //
    // 					int rf = (*iterHsp)->GetQuery_frame();
    //
    //
    // 					if (rf > 0)	//NA committed
    // 					{
    // 						currAlign.m_ReadingFrame = READINGFRAME::PlusId2PlusRFIdx(rf, dst.m_uiSeqLen);
    // 						if (dst.m_dimTranslated[0].empty())	//not filled yet
    // 							dst.NACommit();
    // 					}
    // 					else if (rf < 0)	//NA determined.
    // 					{
    // 						currAlign.m_ReadingFrame = READINGFRAME::MinusId2PlusRFIdx(rf, dst.m_uiSeqLen) | ((currAlign.m_iTo / 3) << 2);
    // 						if (dst.m_dimTranslated[0].empty())	//not filled yet
    // 							dst.NACommit();
    // 					}
    // 					else if (dst.m_dimTranslated[0].empty())	//not filled yet
    // 						dst.PRCommit();
    //
    //
    //
    // 					// -- start to parse sequence
    // 					const string &qseq = (*iterHsp)->GetQseq(), hseq = (*iterHsp)->GetHseq();
    // 					SeqLen_t alnLen = qseq.size();
    // 					SeqPos_t alnResIdx = 0;	//assume qseq.size() == hseq.size()
    //
    // 					// -- incase of blast problem...
    // 					// -- trim unaligned regions from beginning
    // 					while (alnResIdx < alnLen && ('-' == qseq[alnResIdx] || '-' == hseq[alnResIdx])) ++alnResIdx;
    //
    // 					// -- see if NA
    // 					if (dst.m_dimTranslated[0].empty())	//not filled yet
    // 					{
    // 						SeqLen_t aaRes = 0;
    // 						SeqPos_t t0 = alnResIdx;
    //
    // 						while (t0 < alnLen)
    // 						{
    // 							if ('-' == qseq[t0])
    // 								++t0;
    // 							else
    // 							{
    // 								size_t t1 = qseq.find('-', t0);
    // 								if (string::npos == t1)	//
    // 									t1 = alnLen;
    // 								aaRes += t1 - t0;
    // 								t0 = t1 + 1;
    // 							}
    // 						}
    // 						if (aaRes == currAlign.m_iTo - currAlign.m_iFrom + 1)
    // 							dst.PRCommit();
    // 						else if (aaRes * 3 == currAlign.m_iTo - currAlign.m_iFrom + 1)
    // 							dst.NACommit();
    // 					}
    //
    // 					//if (dst.m_iMolType == CSeq_inst::eMol_aa)
    // 					if (!dst.m_bIsNa)
    // 					{
    // 						currAlign.m_ReadingFrame = 0;	//0 is invalid for real reading frames, but a sign for protein sequence
    // 						currAlign.m_eAlignType = TSeqAlignment::eNormal;
    // 						//currAlign.m_eStrand = eNa_strand_unknown;
    // 						currAlign.m_bIsMinus = false;
    // 						SeqPos_t qfrom = currAlign.m_iFrom, hfrom = (*iterHsp)->GetHit_from() - COORDSBASE;
    //
    // 						dst.ParseAlignStrings(qseq, hseq, qfrom, hfrom, alnResIdx, alnLen, 0, currAlign, 0);
    // 						dst.m_strSeqData = dst.m_dimTranslated[0];
    // 					}
    // 					else// if (dst.m_iMolType == CSeq_inst::eMol_na)	//na
    // 					{
    //
    //
    // 						currAlign.m_eAlignType = TSeqAlignment::ePr2Na;
    // 						SeqLen_t alignedLen = currAlign.m_ReadingFrame >> 2;
    //
    //
    // 						// -- calculate
    // 						int rfidx = currAlign.m_ReadingFrame & READINGFRAME::RF_SIZE;
    //
    // 						if (0 == alignedLen)	//plus
    // 						{
    // 							//currAlign.m_eStrand = eNa_strand_plus;
    // 							currAlign.m_bIsMinus = false;
    // 							//currAlign.m_ReadingFrame -= 1;	//from 1, 2, 3 to 0, 1 ,2
    //
    // 							SeqPos_t qfrom = READINGFRAME::PlusNA2PlusPr(currAlign.m_iFrom, currAlign.m_ReadingFrame, dst.m_uiSeqLen),
    // 								hfrom = (*iterHsp)->GetHit_from() - COORDSBASE;
    //
    // 							//SeqPos_t qfrom = currAlign.NAPlus2Pr(currAlign.m_iFrom),
    //
    //
    // 							dst.ParseAlignStrings(qseq, hseq, qfrom, hfrom, alnResIdx, alnLen, 0, currAlign, rfidx);
    //
    //
    // 						}	//plus strand
    // 						else	//minus strand
    // 						{
    // 							rfidx = READINGFRAME::PlusRFIdx2MinusRFIdx(rfidx, dst.m_uiSeqLen);
    //
    // 							//currAlign.m_eStrand = eNa_strand_minus;
    // 							currAlign.m_bIsMinus = true;
    // 							// -- squeeze m_uiSeqLen in currAlign.m_ReadingFrame, rf -1, -2, -3 -> 0, 1, 2
    // 							//currAlign.m_ReadingFrame = (-currAlign.m_ReadingFrame - 1) | (m_uiSeqLen << 2);
    // 							int qfrom = READINGFRAME::PlusNA2MinusPr(currAlign.m_iTo, rfidx, dst.m_uiSeqLen),
    // 								hfrom = (*iterHsp)->GetHit_from() - COORDSBASE;
    //
    // 							dst.ParseAlignStrings(qseq, hseq, qfrom, hfrom, alnResIdx, alnLen, dst.m_dimTranslated[rfidx].size() - alignedLen, currAlign, rfidx);
    //
    // 						}
    // 					}
    // 					...
    // 				}
    // 			}
    // ```
    // The getters of OPTIONAL members without a default throw `CUnassignedMember` when the
    // document lacks them (`unassigned_member`, quoted there): `Hsp_align-len`,
    // `Hsp_identity` and `Hsp_query-frame`, in this order, after the E-value test. What the
    // types decide: `GetQuery_from() - COORDSBASE` (and `GetQuery_to()`'s, `GetHit_from()`'s)
    // overflows `int` for `-2147483648` (undefined): computed as its value, `INT_MIN - 1`, which
    // `m_iFrom`, `m_iTo` and `hfrom` keep (PD-LOSAT-NCBI-DEFECTS exception 11, as for an archive's
    // coordinates; PS16-Q1 oracles `xqfmin`, `xqtmin`, `xhfmin` and their nearby inputs);
    // `m_iTo / 3` is then an `int` again, and `(m_iTo / 3) << 2` is C++20's `int` shift (NCBI's
    // builds are C++20: modulo 2^32), or-ed into the `unsigned int` frame; `rf - 1`,
    // `MinusId2PlusRFIdx` and the frame conversions are `unsigned int` arithmetic (`frames`, of
    // the values modulo 2^32); `m_iTo - m_iFrom + 1` is `int` (computed as its value where it
    // overflows, exception 11) compared as `unsigned int` (modulo 2^32) with the `unsigned int`
    // `aaRes`, which `t1 - t0` (a `size_t`) adds to modulo 2^32;
    // `m_dimTranslated[rfidx].size() - alignedLen` is a `size_t` difference passed as the
    // `unsigned int` `qofs`. `hseq` is read at the indices of `qseq`: its terminating 0 at its
    // length, past it undefined (refused). The values that nothing prints or compares are not
    // kept: `m_uiAlignedLen`, `m_iScore` (whose `(int)` of a score out of `int`'s range or of
    // NaN is a conversion that C++ leaves undefined; it is stored and never read), `m_iNumIdent`,
    // `m_dSeqIdentity`, `m_bIsMinus` (the frame's bits decide the strand everywhere it is read,
    // `domains::rf_index`).
    fn parse_hsp(&mut self, pssm_id: u32, hsp: &Hsp, evcut: f64) -> Result<(), ProcessError> {
        let eval = hsp.evalue;
        if eval > evcut {
            return Ok(());
        }
        let member = |name| unassigned_member("NCBI-BlastOutput", "Hsp", name);
        hsp.align_len.ok_or_else(|| member("align-len"))?;
        hsp.identity.ok_or_else(|| member("identity"))?;
        let coordinate = |value: i32| i64::from(value) - COORDSBASE;
        let from = coordinate(hsp.query_from);
        let to = coordinate(hsp.query_to);
        let rf = hsp.query_frame.ok_or_else(|| member("query-frame"))?;
        let mut align = DomAlign::new(pssm_id, from, to, eval, hsp.bit_score);
        if rf > 0 {
            align.reading_frame = plus_id_to_plus_rf_idx(rf);
            if self.translated[0].is_empty() {
                self.na_commit()?;
            }
        } else if rf < 0 {
            align.reading_frame =
                minus_id_to_plus_rf_idx(rf, self.length) | (((to / 3) as u32) << 2);
            if self.translated[0].is_empty() {
                self.na_commit()?;
            }
        } else if self.translated[0].is_empty() {
            self.pr_commit()?;
        }
        let (qseq, hseq) = (&hsp.qseq[..], &hsp.hseq[..]);
        let aln_len = alignment_length(qseq.len())?;
        let mut aln_res_idx: i32 = 0;
        while (aln_res_idx as u32) < aln_len
            && (qseq[aln_res_idx as usize] == b'-' || hseq_at(hseq, aln_res_idx)? == b'-')
        {
            aln_res_idx += 1;
        }
        if self.translated[0].is_empty() {
            let mut aa_res: u32 = 0;
            let mut t0 = aln_res_idx;
            while (t0 as u32) < aln_len {
                if qseq[t0 as usize] == b'-' {
                    t0 += 1;
                } else {
                    let start = t0 as usize;
                    let t1 = qseq[start..]
                        .iter()
                        .position(|&c| c == b'-')
                        .map_or(aln_len as usize, |gap| start + gap);
                    aa_res = aa_res.wrapping_add((t1 - start) as u32);
                    t0 = (t1 + 1) as i32;
                }
            }
            let span = (to - from + 1) as u32;
            if aa_res == span {
                self.pr_commit()?;
            } else if aa_res.wrapping_mul(3) == span {
                self.na_commit()?;
            }
        }
        let hfrom = coordinate(hsp.hit_from);
        if !self.is_na {
            align.reading_frame = 0;
            self.parse_align_strings(
                qseq,
                hseq,
                from,
                hfrom,
                aln_res_idx,
                aln_len,
                0,
                &mut align,
                0,
            )?;
            self.seq_data_is_translated = true;
        } else {
            align.pr2na = true;
            let aligned_len = align.reading_frame >> 2;
            let rfidx = align.reading_frame & RF_SIZE;
            if aligned_len == 0 {
                let qfrom = i64::from(plus_na_to_plus_pr(from, align.reading_frame));
                if rfidx == 0 {
                    self.keep_seq_data();
                }
                self.parse_align_strings(
                    qseq,
                    hseq,
                    qfrom,
                    hfrom,
                    aln_res_idx,
                    aln_len,
                    0,
                    &mut align,
                    rfidx as usize,
                )?;
            } else {
                let rfidx = plus_rf_idx_to_minus(rfidx, self.length);
                let qfrom = i64::from(plus_na_to_minus_pr(to, rfidx, self.length));
                let qofs = (self.translated[rfidx as usize].len() as u64)
                    .wrapping_sub(u64::from(aligned_len)) as u32;
                self.parse_align_strings(
                    qseq,
                    hseq,
                    qfrom,
                    hfrom,
                    aln_res_idx,
                    aln_len,
                    qofs,
                    &mut align,
                    rfidx as usize,
                )?;
            }
        }
        self.aligns.push(align);
        Ok(())
    }

    /// `ParseAlignStrings`: the segments of an alignment from its strings (from `start`), on the
    /// query from residue `qf` (less `qofs`) and on the domain from `hf`, and the query's residues
    /// written into its string of frame `rfidx`.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_cd_align_proc.cpp:51-122
    // ```c++
    // void _TOfflDomQuery::ParseAlignStrings(const string &qseq, const string &hseq, SeqPos_t qf, SeqPos_t hf, SeqPos_t start, SeqLen_t n, SeqLen_t qofs, TDomSeqAlignment &dst, READINGFRAME::TFRAMEINDEX rfidx)
    // {
    //
    // 	SeqLen_t segLen = 0;
    // 	bool bEffSeg = true;
    //
    // 	while (start < n)
    // 	{
    //
    // 		if (bEffSeg)	//initially, both chains are in effect
    // 		{
    // 			if ('-' != qseq[start] && '-' != hseq[start])
    // 				++segLen;
    // 			else
    // 			{
    //
    // 				bEffSeg = false;
    //
    // 				dst.m_vecMStarts.push_back(qf - qofs);
    // 				dst.m_vecSStarts.push_back(hf);
    // 				dst.m_vecLens.push_back(segLen);
    //
    // 				// -- write in master sequence data
    // 				m_dimTranslated[rfidx].replace(qf, segLen, qseq, start - segLen, segLen);
    //
    // 				qf += segLen;
    // 				hf += segLen;
    // 				segLen = 0;
    // 				if ('-' == qseq[start])
    // 					--qf;
    // 				else
    // 				{
    // 					--hf;
    // 					//write in one more residue
    // 					m_dimTranslated[rfidx][qf] = qseq[start];
    // 				}
    // 			}
    // 		}
    // 		else	//not in effective seg
    // 		{
    // 			//!bEffSeg
    // 			bool segStart = true;
    // 			if ('-' != qseq[start])
    // 			{
    // 				++qf;
    // 				m_dimTranslated[rfidx][qf] = qseq[start];
    // 			}
    // 			else
    // 				segStart = false;
    //
    // 			if ('-' != hseq[start])
    // 				++hf;
    // 			else
    // 				segStart = false;
    //
    // 			if ((bEffSeg = segStart)) ++segLen;
    //
    // 		}
    //
    // 		++start;
    // 	}//while
    // 	// -- last segment
    // 	if (bEffSeg && segLen > 0)
    // 	{
    //
    // 		dst.m_vecMStarts.push_back(qf - qofs);
    // 		dst.m_vecSStarts.push_back(hf);
    // 		dst.m_vecLens.push_back(segLen);
    //
    // 		m_dimTranslated[rfidx].replace(qf, segLen, qseq, start - segLen, segLen);
    // 	}
    //
    // }
    // ```
    // (`qf - qofs`, `qf += segLen`, `hf += segLen` and `start - segLen` are `unsigned int`
    // arithmetic (of `qf`'s and `hf`'s values modulo 2^32), stored in an `int` or passed as a
    // `size_t`; `++qf`, `--qf`, `++hf` and `--hf` are `int`, undefined past its limits: computed
    // as their values, which `qf`, `hf` and `m_vecSStarts` keep (`domains::SubjectStarts`;
    // PD-LOSAT-NCBI-DEFECTS exception 11; PS16-Q1 oracles `xhfmax`, `xhfmin` and their nearby
    // inputs). `hf` starts from `Hsp_hit-from` and passes them with a hit near `INT_MAX` or
    // `INT_MIN`; `qf` starts from `m_iFrom` (`INT_MIN - 1` for an `Hsp_query-from` of
    // `-2147483648`), a `size_t` far past the string's length that the first `replace` throws at.
    // `replace` is libstdc++'s (`string_replace`); `operator[]` past the string's length is
    // undefined (refused), at its length only a NUL is defined (`write_residue`).)
    #[allow(clippy::too_many_arguments)]
    fn parse_align_strings(
        &mut self,
        qseq: &[u8],
        hseq: &[u8],
        mut qf: i64,
        mut hf: i64,
        mut start: i32,
        n: u32,
        qofs: u32,
        dst: &mut DomAlign,
        rfidx: usize,
    ) -> Result<(), ProcessError> {
        // `x += segLen` of an `int` `x`: an `unsigned int` sum, stored back in the `int`.
        let add_seg_len = |x: i64, seg_len: u32| i64::from((x as u32).wrapping_add(seg_len) as i32);
        let mut seg_len: u32 = 0;
        let mut eff_seg = true;
        let translated = &mut self.translated[rfidx];
        while (start as u32) < n {
            let q = qseq[start as usize];
            if eff_seg {
                if q != b'-' && hseq_at(hseq, start)? != b'-' {
                    seg_len = seg_len.wrapping_add(1);
                } else {
                    eff_seg = false;
                    dst.m_starts.push((qf as u32).wrapping_sub(qofs) as i32);
                    dst.s_starts.push(hf);
                    dst.lens.push(seg_len);
                    string_replace(
                        translated,
                        qf,
                        seg_len,
                        qseq,
                        (start as u32).wrapping_sub(seg_len),
                        seg_len,
                    )?;
                    qf = add_seg_len(qf, seg_len);
                    hf = add_seg_len(hf, seg_len);
                    seg_len = 0;
                    if q == b'-' {
                        qf -= 1;
                    } else {
                        hf -= 1;
                        write_residue(translated, qf, q)?;
                    }
                }
            } else {
                let mut seg_start = true;
                if q != b'-' {
                    qf += 1;
                    write_residue(translated, qf, q)?;
                } else {
                    seg_start = false;
                }
                if hseq_at(hseq, start)? != b'-' {
                    hf += 1;
                } else {
                    seg_start = false;
                }
                eff_seg = seg_start;
                if eff_seg {
                    seg_len = seg_len.wrapping_add(1);
                }
            }
            start += 1;
        }
        if eff_seg && seg_len > 0 {
            dst.m_starts.push((qf as u32).wrapping_sub(qofs) as i32);
            dst.s_starts.push(hf);
            dst.lens.push(seg_len);
            string_replace(
                translated,
                qf,
                seg_len,
                qseq,
                (start as u32).wrapping_sub(seg_len),
                seg_len,
            )?;
        }
        Ok(())
    }

    /// The query of the Iteration (`ParseBlastOutput`'s classification of its hits), or what
    /// it threw or LOSAT refuses.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_cd_align_proc.cpp:286-315
    // ```c++
    // 			// -- calculate
    // 			if (!dst.m_vecAlignments.empty())
    // 			{
    // 				size_t ttlAligns = dst.m_vecAlignments.size();
    // 				//if (CSeq_inst::eMol_aa == dst.m_iMolType)	//protein
    // 				if (!dst.m_bIsNa)	//protein
    // 				{
    // 					vector<size_t> ttlIndice;
    // 					ttlIndice.reserve(ttlAligns);
    // 					for (size_t i = 0; i < ttlAligns; ++i)
    // 						ttlIndice.push_back(i);
    // 					Calculate(dst.m_vecAlignments, ttlIndice, dst.m_dimSplitAligns[0], missed);
    // 				}
    // 				else
    // 				{
    // 					vector<size_t> rfIndice[READINGFRAME::TOTAL_RFS];
    // 					for (size_t i = 0; i < READINGFRAME::TOTAL_RFS; ++i)
    // 						rfIndice[i].reserve(ttlAligns);
    //
    // 					SortReadingFrames(dst.m_vecAlignments, rfIndice, dst.m_uiSeqLen);
    //
    // 					for (size_t i = 0; i < READINGFRAME::TOTAL_RFS; ++i)
    // 						Calculate(dst.m_vecAlignments, rfIndice[i], dst.m_dimSplitAligns[i], missed);
    // 				}
    // 			}
    // ```
    // (A NaN E-value or bit score is kept, as in an archive (`query::process_cd_query`):
    // `domains::calculate` sorts the hits as NCBI's `std::sort` does and refuses only where that
    // sort would read outside its vector. `Calculate` of a frame without hits gives an empty
    // classification. `missed`, the PSSM-IDs that the data lack, is never printed.)
    pub(crate) fn finish(self, data: &CddData) -> Result<DomQuery, ProcessError> {
        if let Some(failure) = self.failure {
            return Err(failure);
        }
        // The sites of a protein query read `m_strSeqData`, those of a nucleotide query
        // `m_dimTranslated` (rpsbproc_main.cpp:1370-1379, quoted at `nucleotide::print_features`).
        let mut translated = self.translated;
        let (residues, translated) = match (self.is_na, self.seq_data_is_translated) {
            (true, true) => (translated[0].clone(), Some(Box::new(translated))),
            (true, false) => (self.seq_data, Some(Box::new(translated))),
            (false, true) => (std::mem::take(&mut translated[0]), None),
            (false, false) => (self.seq_data, None),
        };
        let mut query = DomQuery {
            accession: self.accession,
            title: self.title,
            length: self.length,
            is_na: self.is_na,
            residues,
            aligns: self.aligns,
            frames: Default::default(),
            translated,
        };
        if query.aligns.is_empty() {
            return Ok(query);
        }
        if !query.is_na {
            let all: Vec<usize> = (0..query.aligns.len()).collect();
            query.frames[0] =
                calculate(&mut query.aligns, &all, data).map_err(ProcessError::Unsupported)?;
            return Ok(query);
        }
        let frames = sort_reading_frames(&query.aligns, query.length)
            .map_err(|what| ProcessError::Std(what.to_string()))?;
        for (frame, selected) in frames.iter().enumerate() {
            if !selected.is_empty() {
                query.frames[frame] = calculate(&mut query.aligns, selected, data)
                    .map_err(ProcessError::Unsupported)?;
            }
        }
        Ok(query)
    }
}

/// `SeqLen_t alnLen = qseq.size()` of an `Hsp_qseq` below 2^31 letters. One of 2^31 letters or
/// more is refused under the memory rule (DECISIONS 37, 40; PD-LOSAT-NCBI-DEFECTS 1.17 and 1.19):
/// NCBI's `CBlastOutput` holds it with its `Hsp_hseq` (as long, or one letter shorter: a shorter one
/// is read past its end), 4 GiB or more at once, so NCBI's result depends on the host's memory.
/// (Its `SeqPos_t` positions past `INT_MAX`, `++start` of `ParseAlignStrings` quoted at
/// `OfflineQuery::parse_hsp`, would overflow too.)
fn alignment_length(len: usize) -> Result<u32, ProcessError> {
    if len > i32::MAX as usize {
        return Err(unsupported(
            "an Hsp_qseq of 2^31 or more letters (NCBI's rpsbproc holds it and its Hsp_hseq, 4 GiB or more at once; NCBI's result depends on the host's memory)",
        ));
    }
    Ok(len as u32)
}

/// `hseq[i]` of a `const std::string` at an index of `qseq`: its terminating 0 at its length,
/// past it undefined (refused).
fn hseq_at(hseq: &[u8], i: i32) -> Result<u8, ProcessError> {
    let i = i as usize;
    match i.cmp(&hseq.len()) {
        std::cmp::Ordering::Less => Ok(hseq[i]),
        std::cmp::Ordering::Equal => Ok(0),
        std::cmp::Ordering::Greater => Err(unsupported(
            "an Hsp whose Hsp_hseq is shorter than its Hsp_qseq (NCBI reads past it)",
        )),
    }
}

/// `m_dimTranslated[rfidx][qf] = c`: a write within the string; at its length, the write of the
/// terminating `charT()` that the standard defines (a NUL residue: the string is unchanged); any
/// other write at or past its length is undefined (refused).
// C++17 and C++20 [string.access]/2 (LWG 2475; NCBI's 0.5.1 builds are C++20): "Returns:
// *(begin() + pos) if pos < size(). Otherwise, returns a reference to an object of type charT with
// value charT(), where modifying the object to any value other than charT() leads to undefined
// behavior."
// (A NUL residue comes only from a CDATA section of `Hsp_qseq`: `ReadTagData` drops a NUL
// character otherwise. PS16d review S3-2.)
fn write_residue(text: &mut [u8], at: i64, c: u8) -> Result<(), ProcessError> {
    // `qf` as the `size_t` of `operator[]` (of its value, modulo 2^64).
    let at = at as u64;
    if at == text.len() as u64 && c == 0 {
        return Ok(());
    }
    if at >= text.len() as u64 {
        return Err(unsupported(
            "an alignment whose residue falls outside NCBI's string of the query's frame (NCBI writes past it, undefined)",
        ));
    }
    text[at as usize] = c;
    Ok(())
}

/// libstdc++'s `basic_string::replace(pos, n1, str, pos2, n2)`: `n1` bytes from `pos` (fewer at
/// the end) replaced by `n2` bytes of `str` from `pos2` (fewer at its end); the string may grow.
/// Its `std::out_of_range` when a position is past its string's length, `pos2` checked first.
// NCBI reference (GCC 13 libstdc++): include/c++/13/bits/basic_string.h:395-402,2311-2316,2352-2358
// ```c++
//       size_type
//       _M_check(size_type __pos, const char* __s) const
//       {
// 	if (__pos > this->size())
// 	  __throw_out_of_range_fmt(__N("%s: __pos (which is %zu) > "
// 				       "this->size() (which is %zu)"),
// 				   __s, __pos, this->size());
// 	return __pos;
//       }
// ...
//       basic_string&
//       replace(size_type __pos1, size_type __n1, const basic_string& __str,
// 	      size_type __pos2, size_type __n2 = npos)
//       { return this->replace(__pos1, __n1, __str._M_data()
// 			     + __str._M_check(__pos2, "basic_string::replace"),
// 			     __str._M_limit(__pos2, __n2)); }
// ...
//       basic_string&
//       replace(size_type __pos, size_type __n1, const _CharT* __s,
// 	      size_type __n2)
//       {
// 	__glibcxx_requires_string_len(__s, __n2);
// 	return _M_replace(_M_check(__pos, "basic_string::replace"),
// 			  _M_limit(__pos, __n1), __s, __n2);
//       }
// ```
// (`_M_limit(pos, n)` is `min(n, size() - pos)`. `qf`, an `int` (its value past `int`'s range
// where exception 11 took it there), is the `size_t` `pos`: a negative one is beyond any length.)
fn string_replace(
    text: &mut Vec<u8>,
    pos: i64,
    n1: u32,
    src: &[u8],
    pos2: u32,
    n2: u32,
) -> Result<(), ProcessError> {
    let out_of_range = |pos: u64, size: usize| {
        ProcessError::Std(format!(
            "basic_string::replace: __pos (which is {pos}) > this->size() (which is {size})"
        ))
    };
    let pos2 = pos2 as usize;
    if pos2 > src.len() {
        return Err(out_of_range(pos2 as u64, src.len()));
    }
    let n2 = (n2 as usize).min(src.len() - pos2);
    let pos = pos as u64;
    if pos > text.len() as u64 {
        return Err(out_of_range(pos, text.len()));
    }
    let pos = pos as usize;
    let n1 = (n1 as usize).min(text.len() - pos);
    text.splice(pos..pos + n1, src[pos2..pos2 + n2].iter().copied());
    Ok(())
}

#[cfg(test)]
#[path = "offline_tests.rs"]
mod tests;
