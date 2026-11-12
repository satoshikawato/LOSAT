//! `CSegSet` of rpsbproc's library: a set of residue ranges on one sequence, kept sorted and
//! merged (ranges that touch become one), with the original position of each range's start
//! (`ori_from`), and the set algebra that the domain classification uses.
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/segset.hpp:15-45,62-131
//! ```c++
//! struct TSeg_base: public Cacheable
//! {
//! 	...
//! 	SeqPos_t from;
//! 	SeqPos_t to;
//!
//! 	// -- The ori_from should always be pointing to the position on original sequence that corresponding to the "from" of the current mapped position. ...
//! 	SeqPos_t ori_from;
//! 	...
//! 	TSeg_base(SeqPos_t f = 0, SeqPos_t t = 0): from(f), to(t), ori_from(f) {};
//! 	bool IsValid(void) const;
//! ...
//! template <typename SEGTYPE>
//! class CSegSetTmpl: public Cacheable
//! {
//! 	...
//! 	typedef std::list<SEGTYPE> TSegs;
//! 	CSegSetTmpl(void): m_iFactor(1), m_lstContainer(), m_ulGapThreshold(k_SeqLen_Max) {};
//! ```
//!
//! NCBI keeps the ranges in a `std::list` and changes it in place; LOSAT keeps them in a `Vec`
//! and builds the changed set in one pass over the old ranges (`merge`, `cross`, `clip`), which
//! gives the same ranges in the same order with the same `ori_from`, in time linear in the sizes
//! of the two sets as NCBI's list operations. `m_iFactor` is 1 in every set that PS16a's
//! classification makes (it changes only for nucleotide queries, `Pr2NaConvert`).

/// `TSeg_base`. Its positions are `SeqPos_t` (`int`) in NCBI; LOSAT keeps them in `i64`, so
/// the sums of the set algebra (`from - 1`, `to + 1`, the total residues) do not overflow where
/// NCBI's do (positions near 2^31 of a crafted Dense-seg; undefined in C++, which LOSAT widens
/// under the rule of 2026-10-10, as PD-LOSAT-NCBI-DEFECTS exception 4 does). The positions that
/// come in are `int`s, so every result NCBI defines is the same.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Seg {
    pub(crate) from: i64,
    pub(crate) to: i64,
    pub(crate) ori_from: i64,
}

impl Seg {
    /// `TSeg_base(f, t)`: `ori_from` is `from`.
    pub(crate) fn new(from: i64, to: i64) -> Self {
        Seg {
            from,
            to,
            ori_from: from,
        }
    }

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/segset.cpp:4-7
    // ```c++
    // bool TSeg_base::IsValid(void) const
    // {
    // 	return (to >= from);
    // }
    // ```
    fn is_valid(&self) -> bool {
        self.to >= self.from
    }
}

/// `CSegSet` (`CSegSetTmpl<TSeg_base>`) with `m_iFactor` 1.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct SegSet {
    segs: Vec<Seg>,
}

impl SegSet {
    pub(crate) fn new() -> Self {
        SegSet { segs: Vec::new() }
    }

    /// `GetSegs`.
    #[cfg(test)]
    pub(crate) fn segs(&self) -> &[Seg] {
        &self.segs
    }

    /// A set of the ranges `segs` as they are (for the tests).
    #[cfg(test)]
    pub(crate) fn from_segs(segs: Vec<Seg>) -> Self {
        SegSet { segs }
    }

    /// `IsEmpty`.
    pub(crate) fn is_empty(&self) -> bool {
        self.segs.is_empty()
    }

    /// `AddSeg(f, t)`.
    pub(crate) fn add(&mut self, from: i64, to: i64) {
        self.add_seg(Seg::new(from, to));
    }

    /// `AddSeg(seg)`: inserts the range, merging it with the ranges it overlaps or touches.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/impl/segset_impl.hpp:63-108
    // ```c++
    // void CSegSetTmpl<SEGTYPE>::AddSeg(const SEGTYPE& seg)
    // {
    // 	if (!seg.IsValid()) return;
    // 	typename TSegs::iterator curr = m_lstContainer.begin(), dstEnd = m_lstContainer.end();
    //
    // 	while (dstEnd != curr && curr -> to < seg.from - 1) ++curr;
    // 	if (curr == m_lstContainer.end())	//add to end
    // 	{
    // 		m_lstContainer.emplace(dstEnd, seg);
    // 		return;
    // 	}
    // 	else if (curr->from > seg.to + 1)	// gap big enough
    // 	{
    // 		m_lstContainer.emplace(curr, seg);
    // 		return;
    // 	}
    // 	else	//merge the seg in
    // 	{
    // 		if (curr->from > seg.from)
    // 		{
    // 			curr->ori_from -= (curr->from - seg.from) / m_iFactor;
    // 			curr->from = seg.from;	//extended left
    //
    // 		}
    // 		if (curr->to < seg.to)
    // 		{
    // 			curr->to = seg.to;	//extended right
    // 		}
    //
    // 		typename TSegs::iterator nxt = curr;
    //
    // 		++nxt;
    // 		while (nxt != dstEnd && nxt->from <= curr->to + 1)
    // 		{
    // 			if (curr->to < nxt->to)
    // 			{
    // 				curr->to = nxt->to;	//merge next
    // 			}
    // 			m_lstContainer.erase(nxt);
    // 			nxt = curr;
    // 			++nxt;
    // 		}
    // 	}
    // }
    // ```
    pub(crate) fn add_seg(&mut self, seg: Seg) {
        if !seg.is_valid() {
            return;
        }
        let mut curr = 0;
        while curr < self.segs.len() && self.segs[curr].to < seg.from - 1 {
            curr += 1;
        }
        if curr == self.segs.len() {
            self.segs.push(seg);
            return;
        }
        if self.segs[curr].from > seg.to + 1 {
            self.segs.insert(curr, seg);
            return;
        }
        let c = &mut self.segs[curr];
        if c.from > seg.from {
            c.ori_from -= c.from - seg.from;
            c.from = seg.from;
        }
        if c.to < seg.to {
            c.to = seg.to;
        }
        let mut next = curr + 1;
        while next < self.segs.len() && self.segs[next].from <= self.segs[curr].to + 1 {
            if self.segs[curr].to < self.segs[next].to {
                self.segs[curr].to = self.segs[next].to;
            }
            next += 1;
        }
        self.segs.drain(curr + 1..next);
    }

    /// `Merge(src)`: the union, keeping this set's `ori_from` where ranges merge.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/impl/segset_impl.hpp:111-159
    // ```c++
    // void CSegSetTmpl<SEGTYPE>::Merge(const CSegSetTmpl<SEGTYPE> &src)
    // {
    // 	typename TSegs::iterator curr = m_lstContainer.begin();
    // 	typename TSegs::const_iterator scur = src.m_lstContainer.begin();
    //
    // 	while (curr!= m_lstContainer.end() && scur != src.m_lstContainer.end())
    // 	{
    // 		if (curr->to < scur->from - 1) ++curr;
    // 		else if (scur->to < curr->from - 1)
    // 		{
    // 			m_lstContainer.emplace(curr, *scur);
    // 			++scur;
    // 		}
    // 		else
    // 		{
    // 			// -- merge this src seg
    // 			if (curr->from > scur->from)
    // 			{
    // 				curr->ori_from -= (curr->from - scur->from) / m_iFactor;
    // 				curr->from = scur->from;
    // 			}
    // 			if (curr->to < scur->to)
    // 			{
    // 				curr->to = scur->to;
    // 			}
    // 			// -- to next src seg
    // 			++scur;
    //
    // 			// maintain this set
    // 			typename TSegs::iterator nxt = curr;
    // 			++nxt;
    // 			while (nxt != m_lstContainer.end() && nxt->from <= curr->to + 1)
    // 			{
    // 				if (curr->to < nxt->to)
    // 				{
    // 					curr->to = nxt->to;	//merge next
    // 				}
    // 				m_lstContainer.erase(nxt);
    // 				nxt = curr;
    // 				++nxt;
    // 			}
    // 		}
    // 	}
    // 	// -- add the rest segs in
    //
    // 	while (scur != src.m_lstContainer.end())
    // 		m_lstContainer.emplace_back(*scur++);
    // }
    // ```
    pub(crate) fn merge(&mut self, src: &SegSet) {
        let mut old = std::mem::take(&mut self.segs);
        let mut out = Vec::with_capacity(old.len() + src.segs.len());
        let (mut curr, mut scur) = (0, 0);
        while curr < old.len() && scur < src.segs.len() {
            let s = src.segs[scur];
            if old[curr].to < s.from - 1 {
                out.push(old[curr]);
                curr += 1;
            } else if s.to < old[curr].from - 1 {
                out.push(s);
                scur += 1;
            } else {
                let c = &mut old[curr];
                if c.from > s.from {
                    c.ori_from -= c.from - s.from;
                    c.from = s.from;
                }
                if c.to < s.to {
                    c.to = s.to;
                }
                scur += 1;
                // The ranges after `curr` that it now reaches are erased: `curr` moves to the
                // place of the last of them, which the loop never reads again.
                let mut next = curr + 1;
                while next < old.len() && old[next].from <= old[curr].to + 1 {
                    if old[curr].to < old[next].to {
                        old[curr].to = old[next].to;
                    }
                    next += 1;
                }
                old[next - 1] = old[curr];
                curr = next - 1;
            }
        }
        out.extend_from_slice(&old[curr..]);
        out.extend_from_slice(&src.segs[scur..]);
        self.segs = out;
    }

    /// `Cross(src)`: the intersection.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/impl/segset_impl.hpp:161-208
    // ```c++
    // void CSegSetTmpl<SEGTYPE>::Cross(const CSegSetTmpl<SEGTYPE> & src)
    // {
    // 	typename TSegs::iterator curr = m_lstContainer.begin();
    // 	typename TSegs::const_iterator scur = src.m_lstContainer.begin();
    //
    // 	while (curr!= m_lstContainer.end() && scur != src.m_lstContainer.end())
    // 	{
    // 		if (curr->to < scur->from)
    // 		{
    // 			typename TSegs::iterator del = curr;
    // 			++curr;
    // 			m_lstContainer.erase(del);
    // 		}
    // 		else if (scur->to < curr->from)
    // 		{
    // 			++scur;
    // 		}
    // 		else
    // 		{
    // 			if (curr->from < scur->from)
    // 			{
    // 				curr->ori_from += (scur->from - curr->from) / m_iFactor;
    // 				curr->from = scur->from;
    // 			}
    // 			if (curr->to > scur->to + 1)
    // 			{
    // 				SEGTYPE temp(curr->from, scur->to);
    // 				temp.ori_from = curr->ori_from;
    // 				m_lstContainer.emplace(curr, temp);
    //
    // 				curr->ori_from += (scur->to + 2 - curr->from) / m_iFactor;
    // 				curr->from = scur->to + 2;
    // 				++scur;
    // 			}
    // 			else
    // 			{
    // 				if (curr->to > scur->to)
    // 				{
    // 					curr->to = scur->to;
    // 					++scur;
    // 				}
    // 				++curr;
    // 			}
    // 		}
    // 	}
    // 	m_lstContainer.erase(curr, m_lstContainer.end());
    // }
    // ```
    pub(crate) fn cross(&mut self, src: &SegSet) {
        let mut old = std::mem::take(&mut self.segs);
        let mut out = Vec::with_capacity(old.len() + src.segs.len());
        let (mut curr, mut scur) = (0, 0);
        while curr < old.len() && scur < src.segs.len() {
            let s = src.segs[scur];
            if old[curr].to < s.from {
                curr += 1;
            } else if s.to < old[curr].from {
                scur += 1;
            } else {
                let c = &mut old[curr];
                if c.from < s.from {
                    c.ori_from += s.from - c.from;
                    c.from = s.from;
                }
                if c.to > s.to + 1 {
                    let mut temp = Seg::new(c.from, s.to);
                    temp.ori_from = c.ori_from;
                    out.push(temp);
                    c.ori_from += s.to + 2 - c.from;
                    c.from = s.to + 2;
                    scur += 1;
                } else {
                    if c.to > s.to {
                        c.to = s.to;
                        scur += 1;
                    }
                    out.push(*c);
                    curr += 1;
                }
            }
        }
        self.segs = out;
    }

    /// `Clip(src)`: this set less `src`.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/impl/segset_impl.hpp:210-248
    // ```c++
    // void CSegSetTmpl<SEGTYPE>::Clip(const CSegSetTmpl<SEGTYPE>& src)
    // {
    // 	typename TSegs::iterator curr = m_lstContainer.begin();
    // 	typename TSegs::const_iterator scur = src.m_lstContainer.begin();
    //
    // 	while (curr!= m_lstContainer.end() && scur != src.m_lstContainer.end())
    // 	{
    // 		if (curr->to < scur->from)
    // 		{
    // 			++curr;
    // 		}
    // 		else if (scur->to < curr->from)
    // 		{
    // 			++scur;
    // 		}
    // 		else
    // 		{
    // 			if (curr->from < scur->from)
    // 			{
    // 				SEGTYPE temp(curr->from, scur->from - 1);
    // 				temp.ori_from = curr->ori_from;
    // 				m_lstContainer.emplace(curr, temp);
    // 			}
    // 			if (curr->to <= scur->to)
    // 			{
    // 				typename TSegs::iterator del = curr;
    // 				++curr;
    // 				m_lstContainer.erase(del);
    // 			}
    // 			else
    // 			{
    // 				curr->ori_from += (scur->to + 1 - curr->from) / m_iFactor;
    // 				curr->from = scur->to + 1;
    // 				++scur;
    // 			}
    // 		}
    // 	}
    // }
    // ```
    pub(crate) fn clip(&mut self, src: &SegSet) {
        let mut old = std::mem::take(&mut self.segs);
        let mut out = Vec::with_capacity(old.len() + src.segs.len());
        let (mut curr, mut scur) = (0, 0);
        while curr < old.len() && scur < src.segs.len() {
            let s = src.segs[scur];
            if old[curr].to < s.from {
                out.push(old[curr]);
                curr += 1;
            } else if s.to < old[curr].from {
                scur += 1;
            } else {
                let c = &mut old[curr];
                if c.from < s.from {
                    let mut temp = Seg::new(c.from, s.from - 1);
                    temp.ori_from = c.ori_from;
                    out.push(temp);
                }
                if c.to <= s.to {
                    curr += 1;
                } else {
                    c.ori_from += s.to + 1 - c.from;
                    c.from = s.to + 1;
                    scur += 1;
                }
            }
        }
        out.extend_from_slice(&old[curr..]);
        self.segs = out;
    }

    /// `GetTotalResidues`.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/impl/segset_impl.hpp:293-302
    // ```c++
    // SeqPos_t CSegSetTmpl<SEGTYPE>::GetTotalResidues(void) const
    // {
    // 	SeqPos_t iResult = 0;
    // 	for (typename TSegs::const_iterator iter = m_lstContainer.begin(); iter != m_lstContainer.end(); ++iter)
    // 	{
    // 		iResult += iter->to - iter->from + 1;
    // 	}
    // 	return iResult/abs(m_iFactor);
    // }
    // ```
    pub(crate) fn total_residues(&self) -> i64 {
        self.segs.iter().map(|seg| seg.to - seg.from + 1).sum()
    }
}

#[cfg(test)]
#[path = "segset_tests.rs"]
mod tests;
