//! The motif check of a site (`TDomSite::MotifCheck`) with NCBI's process-wide motif position
//! list.
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:669-740
//! ```c++
//! int TDomSite::MotifCheck(const vector<TSeg_base::TResiduePos> &rMappedRes, const string &aaSeq) const
//! {
//! 	static vector<SeqPos_t> st_motif_pos;
//!
//! 	if (!m_strMotif.empty())
//! 	{
//! 		if (st_motif_pos.empty())
//! 			GetMotifResPos(st_motif_pos);
//! 		...
//! ```
//! `st_motif_pos` is a function-static: the positions of the first motif-bearing site that is
//! checked in a process serve every later check of any site, domain and query of that process
//! (survey 8.1). In LOSAT a process is one `losat rpsbproc` run, which may share a process
//! with other runs (the library, the Web), so the list belongs to the run: `MotifState`.

use super::data::{CddData, Site};
use super::prosite::{MapEnd, MatchError, ProSite};

/// A mapped residue (`TSeg_base::TResiduePos`): its position on the query (`curr`) and the
/// position on the model it came from (`ori`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ResiduePos {
    pub(crate) curr: i32,
    pub(crate) ori: i32,
}

/// NCBI's `static vector<SeqPos_t> st_motif_pos` of `MotifCheck`, for one run.
#[derive(Debug, Default)]
pub(crate) struct MotifState {
    positions: Vec<i32>,
    /// Counts the reads of `st_motif_pos` past its end (the cases of the Owner's rule below),
    /// for the tests.
    #[cfg(test)]
    pub(crate) past_end_reads: usize,
}

/// Why a motif check does not give 0, 1 or 2.
#[derive(Debug, PartialEq, Eq)]
pub(crate) enum MotifError {
    /// `THROW_SIMPLE("Motif string parse error -- Motif = " << m_strMotif << ", error position: " << errPos)`.
    Parse { motif: Vec<u8>, position: usize },
    /// `GetMinimalXMap` (or the `push_back` of the text) throws, or holds what the host decides
    /// (`ProSite::minimal_map`).
    Map(MapEnd),
    /// The match of the residues (`ProSite::matches`).
    Match(MatchError),
}

impl MotifState {
    /// `TDomSite::MotifCheck`: 0 = the site passes (or has no motif), 1 = the essential
    /// positions are not complete, 2 = they are but the residues do not fit the motif.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:669-740
    // ```c++
    // 	if (!m_strMotif.empty())
    // 	{
    // 		if (st_motif_pos.empty())
    // 			GetMotifResPos(st_motif_pos);
    //
    // 		CProSite ps;
    // 		size_t errPos;
    // 		CProSite::EParseError err = ps.Parse(m_strMotif, errPos);
    // 		if (CProSite::eNoError != err)
    // 		{
    // 			THROW_SIMPLE("Motif string parse error -- Motif = " << m_strMotif << ", error position: " << errPos);
    // 		}
    // 		string minMap(k_strEmptyString);
    // 		ps.GetMinimalXMap(minMap);
    //
    // 		size_t mtfIdx = 0, mtfLen = minMap.size();
    // 		size_t mapIdx = 0, mappedLen = rMappedRes.size();
    //
    // 		while (mtfIdx < mtfLen)
    // 		{
    // 			while (mapIdx < mappedLen)
    // 			{
    // 				if (rMappedRes[mapIdx].ori > st_motif_pos[mtfIdx])	//failed
    // 					return 1;
    // 				else if (rMappedRes[mapIdx].ori == st_motif_pos[mtfIdx])
    // 				{
    // 					++mapIdx;
    // 					goto labelContinue;
    // 				}
    // 				++mapIdx;
    //
    // 			}
    // 			// -- exhausted, fail
    // 			return 1;
    // 		labelContinue:
    // 			++mtfIdx;
    // 		}
    //
    // 		size_t seqLen = aaSeq.size();
    // 		if (!aaSeq.empty())
    // 		{
    // 			minMap.clear();	//borrow this for other use
    //
    // 			for (size_t i = 0; i < mappedLen; ++i)
    // 				if ((size_t)rMappedRes[i].curr < seqLen)
    // 				{
    // 					minMap.push_back(aaSeq[rMappedRes[i].curr]);
    // 				}
    //
    // 			size_t endPos = ps.Match(minMap, seqLen, 0);
    // 			if (string::npos == endPos)
    // 			{
    // 				return 2;
    // 			}
    // 		}
    // 	}
    //
    // 	return 0;
    // ```
    // `GetMotifResPos` (basedata.cpp:605-625) lists the positions of the site's ranges.
    //
    // Owner's rule, DECISIONS 25 Q2 = A (2026-11-11): `st_motif_pos[mtfIdx]` is read past the
    // end of the list when the motif's minimal map is longer than the first motif site's
    // positions (`mtfIdx >= st_motif_pos.size()`, the compare at :704 reads heap bytes: C++
    // undefined behavior). LOSAT treats the check as failed (returns 1), as the loop does when
    // the mapped residues are exhausted. PD-LOSAT-NCBI-DEFECTS exception 10 (version 1.14).
    pub(crate) fn motif_check(
        &mut self,
        data: &CddData,
        site: &Site,
        mapped: &[ResiduePos],
        residues: &[u8],
    ) -> Result<u8, MotifError> {
        let motif = data.text(site.motif);
        if motif.is_empty() {
            return Ok(0);
        }
        if self.positions.is_empty() {
            for seg in data.site_segs(site) {
                // `for (SeqPos_t c = v->from; c <= v->to; ++c)`
                let mut c = i64::from(seg[0]);
                while c <= i64::from(seg[1]) {
                    self.positions.push(c as i32);
                    c += 1;
                }
            }
        }
        let pattern = ProSite::parse(motif).map_err(|(_, position)| MotifError::Parse {
            motif: motif.to_vec(),
            position,
        })?;
        let map = pattern.minimal_map().map_err(MotifError::Map)?;
        let motif_len = map.len;
        let (mut motif_idx, mut map_idx) = (0usize, 0usize);
        while (motif_idx as u64) < motif_len {
            loop {
                if map_idx >= mapped.len() {
                    // Exhausted, fail.
                    return Ok(1);
                }
                let Some(&position) = self.positions.get(motif_idx) else {
                    // NCBI reads past the end of `st_motif_pos` here: the Owner's rule.
                    #[cfg(test)]
                    {
                        self.past_end_reads += 1;
                    }
                    return Ok(1);
                };
                let ori = mapped[map_idx].ori;
                if ori > position {
                    return Ok(1);
                }
                map_idx += 1;
                if ori == position {
                    break;
                }
            }
            motif_idx += 1;
        }
        if !residues.is_empty() {
            // `minMap` holds the residues now: its capacity bounds the bytes beyond them that
            // are NCBI's own, and its heap buffer lives during the match (`ProSite::matches`).
            let mut buffer = map;
            buffer.clear();
            let mut text = Vec::new();
            for res in mapped {
                if res.curr >= 0 && (res.curr as usize) < residues.len() {
                    text.push(residues[res.curr as usize]);
                    buffer.push_back().map_err(MotifError::Map)?;
                }
            }
            return match pattern.matches(&text, residues.len(), 0, buffer) {
                Ok(false) => Ok(2),
                Ok(true) => Ok(0),
                Err(error) => Err(MotifError::Match(error)),
            };
        }
        Ok(0)
    }
}

#[cfg(test)]
#[path = "sites_tests.rs"]
mod tests;
