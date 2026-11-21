//! `CSegSet` of rpsbproc's library: a set of residue ranges on one sequence, kept sorted and
//! merged (ranges that touch become one), with the original position of each range's start
//! (`ori_from`), and the set algebra that the domain classification uses.
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/segset.hpp:15-45,62-131
//!
//! NCBI keeps the ranges in a `std::list` and changes it in place; LOSAT keeps them in a `Vec`
//! and builds the changed set in one pass over the old ranges (`merge`, `cross`, `clip`), which
//! gives the same ranges in the same order with the same `ori_from`, in time linear in the sizes
//! of the two sets as NCBI's list operations. `m_iFactor` is 1 in every set that PS16a's
//! classification makes (`Pr2NaConvert` changes it for a nucleotide query's motifs, whose ends
//! only are printed: `domains::pr2na_ends`).

use super::domains::SubjectStarts;

/// `TSeg_base`. Its positions are `SeqPos_t` (`int`) in NCBI; LOSAT keeps them in `i64`, so
/// the `int` sums of the set algebra (`from - 1`, `to + 1`, the total residues, the mapping's
/// `diff` and `ori_from` sums) do not overflow where NCBI's do (positions near 2^31 of a crafted
/// alignment; undefined in C++): they are computed as their values and kept, under
/// PD-LOSAT-NCBI-DEFECTS exception 12 (the Owner's rule of 2026-10-10 and Q1 = B of 2026-11-14).
/// The positions that come in are `int`s and the `unsigned int` arithmetic is done modulo 2^32,
/// so every result NCBI defines is the same.
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

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/impl/segset_impl.hpp:281-291
    pub(crate) fn left(&self) -> i64 {
        self.segs.first().map_or(-1, |seg| seg.from)
    }

    /// `GetRight`.
    pub(crate) fn right(&self) -> i64 {
        self.segs.last().map_or(-1, |seg| seg.to)
    }

    /// `GetOriTo(citer)` of a range: `ori_from + (to - from) / m_iFactor`.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/impl/segset_impl.hpp:332-337
    pub(crate) fn ori_to(seg: &Seg) -> i64 {
        seg.ori_from + (seg.to - seg.from)
    }

    /// `TSeqAlignment::MapSegSet(segset, false)` of a protein alignment: this set (positions on
    /// the model, the slave row) turned into positions on the query (the master row) along the
    /// aligned segments (`m_vecMStarts`, `m_vecSStarts`, `m_vecLens`), keeping `ori_from`.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:1180-1250
    // (rpsbproc calls it with `doConvert` false everywhere; the conversion of a nucleotide
    // query's motifs is `domains::pr2na_ends`. NCBI changes the list in place; LOSAT writes the
    // new ranges in order into a new list. `m_vecSStarts[idx] + m_vecLens[idx] - 1` adds an
    // `int` and an `unsigned int`: an `unsigned int` sum, modulo 2^32, stored in the `int` `end`;
    // the other sums are `int` arithmetic, computed as their values (`Seg`).)
    pub(crate) fn map_through(&mut self, m_starts: &[i32], s_starts: &SubjectStarts, lens: &[u32]) {
        if self.segs.is_empty() {
            return;
        }
        let old = std::mem::take(&mut self.segs);
        let mut out = Vec::with_capacity(old.len() + 1);
        let mut idx = 0usize;
        let mut source = old.into_iter();
        let mut current = source.next();
        while let Some(mut iter) = current {
            if idx >= lens.len() || s_starts.get(idx) > iter.to {
                current = source.next();
                continue;
            }
            let start = s_starts.get(idx);
            let diff = start - i64::from(m_starts[idx]);
            let end = i64::from((start as u32).wrapping_add(lens[idx]).wrapping_sub(1) as i32);
            if end >= iter.to {
                if iter.from < start {
                    iter.ori_from += start - iter.from;
                    iter.from = start;
                }
                if end == iter.to {
                    idx += 1;
                }
                iter.from -= diff;
                iter.to -= diff;
                out.push(iter);
                current = source.next();
            } else if end >= iter.from {
                let mut temp = Seg::new(iter.from, end);
                temp.ori_from = iter.ori_from;
                if temp.from < start {
                    temp.ori_from += start - temp.from;
                    temp.from = start;
                }
                temp.from -= diff;
                temp.to -= diff;
                out.push(temp);
                iter.ori_from += end - iter.from + 1;
                iter.from = end + 1;
                idx += 1;
                current = Some(iter);
            } else {
                idx += 1;
                current = Some(iter);
            }
        }
        self.segs = out;
    }

    /// `GetTotalResidues`.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/impl/segset_impl.hpp:293-302
    pub(crate) fn total_residues(&self) -> i64 {
        self.segs.iter().map(|seg| seg.to - seg.from + 1).sum()
    }
}

#[cfg(test)]
#[path = "segset_tests.rs"]
mod tests;
