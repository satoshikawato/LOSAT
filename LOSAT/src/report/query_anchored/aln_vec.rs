//! NCBI's `CAlnMap` and `CAlnVec` (`objtools/alnmgr/alnmap.cpp`, `alnvec.cpp`) over the merged
//! Dense-seg of a query's section (`aln_mix::merge`), as the query-anchored reports read it:
//! the anchor (`SetAnchor`), the alignment coordinates of the segments and of each row
//! (`GetAlnStart`, `GetAlnStop`, `GetSeqAlnRange`, `GetSeqStart`, `GetSeqStop`),
//! `GetAlnPosFromSeqPos` with `eBackwards` and `try_reverse_dir` (the masks, `x_FillLocList`)
//! and `GetWholeAlnSeqString` with its insertions and per-line coordinates.
//!
//! Every row is a plus-strand row of width 1 (`aln_mix` module doc): the minus-strand halves of
//! these functions (`plus ? ... : ...`, `last - cur`) and `TranslateNAToAA` are not ported.

use super::{Denseg, SeqKey};

/// NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnmap.hpp:220-233
/// ```c++
///     class CNumSegWithOffset
///     {
///     public:
///         CNumSegWithOffset(TNumseg seg, int offset = 0)
///             : m_AlnSeg(seg), m_Offset(offset) { }
/// ```
#[derive(Clone, Copy)]
struct SegWithOffset {
    aln_seg: i32,
    offset: i32,
}

/// The CAlnMap of a merged Dense-seg (`ds`), with the anchor when set.
pub(crate) struct AlnMap<'d> {
    ds: &'d Denseg,
    num_rows: usize,
    num_segs: usize,
    anchor: Option<usize>,
    /// `m_AlnSegIdx`: the raw segment of each alignment segment (anchored).
    aln_seg_idx: Vec<usize>,
    /// `m_NumSegWithOffsets`: the alignment segment of each raw segment (anchored).
    num_seg_with_offsets: Vec<SegWithOffset>,
    /// `m_AlnStarts`.
    aln_starts: Vec<u32>,
    /// `x_GetSeqLeftSeg` and `x_GetSeqRightSeg` of each row (raw segments).
    left_segs: Vec<usize>,
    right_segs: Vec<usize>,
}

/// The letters, insertions and per-line coordinates of one row (`GetWholeAlnSeqString`).
pub(crate) struct RowString {
    pub(crate) buffer: Vec<u8>,
    /// `insert_aln_starts`, `insert_starts` and `insert_lens`: each insertion's alignment
    /// position (the anchor's columns before it), its first letter and its length.
    pub(crate) inserts: std::collections::VecDeque<(u32, u32, u32)>,
    /// `scrn_lefts` and `scrn_rights`: the row's first and last position on each line, as
    /// `TSeqPos` (a row's `-1` is 4294967295 there).
    pub(crate) lefts: std::collections::VecDeque<u32>,
    pub(crate) rights: std::collections::VecDeque<u32>,
}

impl<'d> AlnMap<'d> {
    /// `CAlnMap(ds)`: no anchor (`x_CreateAlnStarts`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:49-60
    /// ```c++
    /// void CAlnMap::x_CreateAlnStarts(void)
    /// {
    ///     m_AlnStarts.clear();
    ///     m_AlnStarts.reserve(GetNumSegs());
    ///
    ///     int start = 0, len = 0;
    ///     for (int i = 0;  i < GetNumSegs();  ++i) {
    ///         start += len;
    ///         m_AlnStarts.push_back(start);
    ///         len = m_Lens[i];
    ///     }
    /// }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:716-755 (`x_GetSeqLeftSeg`,
    /// `x_GetSeqRightSeg`: the first and the last segment of the row with a start; "Invalid
    /// Dense-seg: Row ... contains gaps only." for a row of none, which `x_CreateDenseg` removed).
    pub(crate) fn new(ds: &'d Denseg) -> Self {
        let num_rows = ds.ids.len();
        let num_segs = ds.lens.len();
        let mut aln_starts = Vec::with_capacity(num_segs);
        let mut start = 0u32;
        let mut len = 0u32;
        for &seg_len in &ds.lens {
            start = start.wrapping_add(len);
            aln_starts.push(start);
            len = seg_len;
        }
        let mut left_segs = vec![0; num_rows];
        let mut right_segs = vec![0; num_rows];
        for row in 0..num_rows {
            let present = |seg: &usize| ds.starts[seg * num_rows + row] >= 0;
            left_segs[row] = (0..num_segs).find(present).unwrap_or(0);
            right_segs[row] = (0..num_segs).rev().find(present).unwrap_or(0);
        }
        Self {
            ds,
            num_rows,
            num_segs,
            anchor: None,
            aln_seg_idx: Vec::new(),
            num_seg_with_offsets: Vec::new(),
            aln_starts,
            left_segs,
            right_segs,
        }
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:79-121
    /// ```c++
    ///     int start = 0, len = 0, aln_seg = -1, offset = 0;
    ///
    ///     m_Anchor = anchor;
    ///     for (int i = 0, pos = m_Anchor;  i < m_NumSegs;
    ///          ++i, pos += m_NumRows) {
    ///         if (m_Starts[pos] != -1) {
    ///             ++aln_seg;
    ///             offset = 0;
    ///             m_AlnSegIdx.push_back(i);
    ///             m_NumSegWithOffsets.push_back(CNumSegWithOffset(aln_seg));
    ///             start += len;
    ///             m_AlnStarts.push_back(start);
    ///             len = m_Lens[i];
    ///         } else {
    ///             ++offset;
    ///             m_NumSegWithOffsets.push_back(CNumSegWithOffset(aln_seg, offset));
    ///         }
    ///     }
    ///     if (!m_AlnSegIdx.size()) {
    ///         NCBI_THROW(CAlnException, eInvalidDenseg,
    ///                    "CAlnVec::SetAnchor(): "
    ///                    "Invalid Dense-seg: No sequence on the anchor row");
    ///     }
    /// ```
    /// (The query row, the anchor, has the letters of every query segment.)
    pub(crate) fn set_anchor(&mut self, anchor: usize) {
        self.aln_seg_idx.clear();
        self.aln_starts.clear();
        self.num_seg_with_offsets.clear();
        let (mut start, mut len, mut aln_seg, mut offset) = (0u32, 0u32, -1i32, 0i32);
        self.anchor = Some(anchor);
        for seg in 0..self.num_segs {
            if self.ds.starts[seg * self.num_rows + anchor] != -1 {
                aln_seg += 1;
                offset = 0;
                self.aln_seg_idx.push(seg);
                self.num_seg_with_offsets
                    .push(SegWithOffset { aln_seg, offset: 0 });
                start = start.wrapping_add(len);
                self.aln_starts.push(start);
                len = self.ds.lens[seg];
            } else {
                offset += 1;
                self.num_seg_with_offsets
                    .push(SegWithOffset { aln_seg, offset });
            }
        }
    }

    pub(crate) fn num_rows(&self) -> usize {
        self.num_rows
    }

    pub(crate) fn seq_key(&self, row: usize) -> SeqKey {
        self.ds.ids[row]
    }

    /// `GetNumSegs`.
    fn aln_num_segs(&self) -> usize {
        if self.anchor.is_some() {
            self.aln_seg_idx.len()
        } else {
            self.num_segs
        }
    }

    /// `x_GetRawSegFromSeg`.
    fn raw_seg_from_seg(&self, seg: usize) -> usize {
        if self.anchor.is_some() {
            self.aln_seg_idx[seg]
        } else {
            seg
        }
    }

    /// `x_GetSegFromRawSeg`.
    fn seg_from_raw_seg(&self, raw_seg: usize) -> SegWithOffset {
        if self.anchor.is_some() {
            self.num_seg_with_offsets[raw_seg]
        } else {
            SegWithOffset {
                aln_seg: raw_seg as i32,
                offset: 0,
            }
        }
    }

    fn raw_start(&self, row: usize, seg: usize) -> i32 {
        self.ds.starts[seg * self.num_rows + row]
    }

    // NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnmap.hpp:481-498
    // ```c++
    // inline TSeqPos CAlnMap::GetAlnStart(TNumseg seg) const
    // {
    //     return m_AlnStarts[seg];
    // }
    //
    // inline
    // TSeqPos CAlnMap::GetAlnStop(TNumseg seg) const
    // {
    //     return m_AlnStarts[seg] + m_Lens[x_GetRawSegFromSeg(seg)] - 1;
    // }
    //
    // inline
    // TSeqPos CAlnMap::GetAlnStop(void) const
    // {
    //     return GetAlnStop(GetNumSegs() - 1);
    // }
    // ```
    fn aln_start(&self, seg: usize) -> u32 {
        self.aln_starts[seg]
    }

    fn seg_aln_stop(&self, seg: usize) -> u32 {
        self.aln_starts[seg]
            .wrapping_add(self.ds.lens[self.raw_seg_from_seg(seg)])
            .wrapping_sub(1)
    }

    /// `GetAlnStop()`.
    pub(crate) fn aln_stop(&self) -> u32 {
        self.seg_aln_stop(self.aln_num_segs() - 1)
    }

    // NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnmap.hpp:665-681
    // ```c++
    // inline
    // TSeqPos CAlnMap::GetSeqStart(TNumrow row) const
    // {
    //     return
    //         m_Starts[(IsPositiveStrand(row) ?
    //                   x_GetSeqLeftSeg(row) :
    //                   x_GetSeqRightSeg(row)) * m_NumRows + row];
    // }
    //
    // inline
    // TSeqPos CAlnMap::GetSeqStop(TNumrow row) const
    // {
    //     const TNumseg& seg = IsPositiveStrand(row) ?
    //         x_GetSeqRightSeg(row) : x_GetSeqLeftSeg(row);
    //     return m_Starts[seg * m_NumRows + row] + x_GetLen(row, seg) - 1;
    // }
    // ```
    /// `GetSeqStart` (the row's first position).
    pub(crate) fn seq_start(&self, row: usize) -> u32 {
        self.raw_start(row, self.left_segs[row]) as u32
    }

    /// `GetSeqStop` (the row's last position).
    pub(crate) fn seq_stop(&self, row: usize) -> u32 {
        let seg = self.right_segs[row];
        (self.raw_start(row, seg) as u32)
            .wrapping_add(self.ds.lens[seg])
            .wrapping_sub(1)
    }

    /// `GetSeqAlnRange`: the row's first and last alignment position.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:969-998
    /// ```c++
    /// TSignedSeqPos CAlnMap::GetSeqAlnStart(TNumrow row) const
    /// {
    ///     if (IsSetAnchor()) {
    ///         TNumseg seg = -1;
    ///         while (++seg < (TNumseg) m_AlnSegIdx.size()) {
    ///             if (m_Starts[m_AlnSegIdx[seg] * m_NumRows + row] >= 0) {
    ///                 return GetAlnStart(seg);
    ///             }
    ///         }
    ///         return -1;
    ///     } else {
    ///         return GetAlnStart(x_GetSeqLeftSeg(row));
    ///     }
    /// }
    ///
    ///
    /// TSignedSeqPos CAlnMap::GetSeqAlnStop(TNumrow row) const
    /// {
    ///     if (IsSetAnchor()) {
    ///         TNumseg seg = TNumseg(m_AlnSegIdx.size());
    ///         while (seg--) {
    ///             if (m_Starts[m_AlnSegIdx[seg] * m_NumRows + row] >= 0) {
    ///                 return GetAlnStop(seg);
    ///             }
    ///         }
    ///         return -1;
    ///     } else {
    ///         return GetAlnStop(x_GetSeqRightSeg(row));
    ///     }
    /// }
    /// ```
    pub(crate) fn seq_aln_range(&self, row: usize) -> (i64, i64) {
        if self.anchor.is_some() {
            let present = |seg: &usize| self.raw_start(row, self.aln_seg_idx[*seg]) >= 0;
            let start = (0..self.aln_seg_idx.len())
                .find(present)
                .map_or(-1, |seg| i64::from(self.aln_start(seg) as i32));
            let stop = (0..self.aln_seg_idx.len())
                .rev()
                .find(present)
                .map_or(-1, |seg| i64::from(self.seg_aln_stop(seg) as i32));
            (start, stop)
        } else {
            (
                i64::from(self.aln_start(self.left_segs[row]) as i32),
                i64::from(self.seg_aln_stop(self.right_segs[row]) as i32),
            )
        }
    }

    /// `GetRawSeg(row, seq_pos, eBackwards, true)` of a plus-strand row.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:400-524
    /// ```c++
    ///     TSignedSeqPos start = -1, sseq_pos = seq_pos;
    ///     TNumseg       btm, top, mid, cur, last, cur_top, cur_btm;
    ///     btm = cur_btm = 0; cur = top = last = cur_top = m_NumSegs - 1;
    ///
    ///     bool plus = IsPositiveStrand(row);
    ///
    ///     // if out-of-range, return either -1 or the closest seg in dir direction
    ///     if (sseq_pos < (TSignedSeqPos)GetSeqStart(row)) {
    ///         if (dir == eNone) {
    ///             return -1;
    ///         } else if (dir == eForward  ||
    ///                    dir == (plus ? eRight : eLeft)  ||
    ///                    try_reverse_dir) {
    ///             TNumseg seg;
    ///             if (plus) {
    ///                 seg = -1;
    ///                 while (++seg < m_NumSegs) {
    ///                     if (m_Starts[seg * m_NumRows + row] >= 0) {
    ///                         return seg;
    ///                     }
    ///                 }
    ///             } else {
    ///                 ...
    ///             }
    ///         }
    ///     } else if (sseq_pos > (TSignedSeqPos)GetSeqStop(row)) {
    ///         if (dir == eNone) {
    ///             return -1;
    ///         } else if (dir == eBackwards  ||
    ///                    dir == (plus ? eLeft : eRight)  ||
    ///                    try_reverse_dir) {
    ///             TNumseg seg;
    ///             if (plus) {
    ///                 seg = m_NumSegs;
    ///                 while (seg--) {
    ///                     if (m_Starts[seg * m_NumRows + row] >= 0) {
    ///                         return seg;
    ///                     }
    ///                 }
    ///             } else {
    ///                 ...
    ///             }
    ///         }
    ///     }
    ///
    ///     // main loop
    ///     while (btm <= top) {
    ///         cur = mid = (top + btm) / 2;
    ///
    ///         while (cur <= top
    ///                &&  (start = m_Starts[(plus ? cur : last - cur)
    ///                                              * m_NumRows + row]) < 0) {
    ///             ++cur;
    ///         }
    ///         if (cur <= top && start >= 0) {
    ///             if (sseq_pos >= start &&
    ///                 seq_pos < start + x_GetLen(row, plus ? cur : last - cur)) {
    ///                 return (plus ? cur : last - cur); // found
    ///             }
    ///             if (sseq_pos > start) {
    ///                 btm = cur + 1;
    ///                 cur_btm = cur;
    ///             } else {
    ///                 top = mid - 1;
    ///                 cur_top = cur;
    ///             }
    ///             continue;
    ///         }
    ///
    ///         cur = mid-1;
    ///         while (cur >= btm &&
    ///                (start = m_Starts[(plus ? cur : last - cur)
    ///                                          * m_NumRows + row]) < 0) {
    ///             --cur;
    ///         }
    ///         if (cur >= btm && start >= 0) {
    ///             if (sseq_pos >= start
    ///                 &&  seq_pos < start + x_GetLen(row, plus ? cur : last - cur)) {
    ///                 return (plus ? cur : last - cur); // found
    ///             }
    ///             if (sseq_pos > start) {
    ///                 btm = mid + 1;
    ///                 cur_btm = cur;
    ///             } else {
    ///                 top = cur - 1;
    ///                 cur_top = cur;
    ///             }
    ///             continue;
    ///         }
    ///
    ///         // if we get here, seq_pos falls into an unaligned region
    ///         // return either -1 or the closest segment in dir direction
    ///         if (dir == eNone) {
    ///             return -1;
    ///         } else if (dir == eBackwards  ||  dir == (plus ? eLeft : eRight)) {
    ///             return (plus ? cur_btm : last - cur_btm);
    ///         } else if (dir == eForward  ||  dir == (plus ? eRight : eLeft)) {
    ///             return (plus ? cur_top : last - cur_top);
    ///         }
    ///     }
    ///
    ///     // return either -1 or the closest segment in dir direction
    ///     if (dir == eNone) {
    ///         return -1;
    ///     } else if (dir == eBackwards  ||  dir == (plus ? eLeft : eRight)) {
    ///         return (plus ? cur_btm : last - cur_btm);
    ///     } else if (dir == eForward  ||  dir == (plus ? eRight : eLeft)) {
    ///         return (plus ? cur_top : last - cur_top);
    ///     }
    ///
    ///     return -1;
    /// ```
    fn raw_seg_backwards(&self, row: usize, seq_pos: u32) -> usize {
        let sseq_pos = seq_pos as i32;
        let present = |seg: i64| self.raw_start(row, seg as usize) >= 0;
        let num_segs = self.num_segs as i64;
        if sseq_pos < self.seq_start(row) as i32 {
            if let Some(seg) = (0..num_segs).find(|&seg| present(seg)) {
                return seg as usize;
            }
        } else if sseq_pos > self.seq_stop(row) as i32 {
            if let Some(seg) = (0..num_segs).rev().find(|&seg| present(seg)) {
                return seg as usize;
            }
        }
        let (mut btm, mut top, mut cur_btm) = (0i64, num_segs - 1, 0i64);
        while btm <= top {
            let mid = (top + btm) / 2;
            let mut cur = mid;
            let mut start = -1;
            while cur <= top && {
                start = self.raw_start(row, cur as usize);
                start < 0
            } {
                cur += 1;
            }
            if cur <= top && start >= 0 {
                if sseq_pos >= start
                    && i64::from(seq_pos) < i64::from(start) + i64::from(self.ds.lens[cur as usize])
                {
                    return cur as usize;
                }
                if sseq_pos > start {
                    btm = cur + 1;
                    cur_btm = cur;
                } else {
                    top = mid - 1;
                }
                continue;
            }
            cur = mid - 1;
            while cur >= btm && {
                start = self.raw_start(row, cur as usize);
                start < 0
            } {
                cur -= 1;
            }
            if cur >= btm && start >= 0 {
                if sseq_pos >= start
                    && i64::from(seq_pos) < i64::from(start) + i64::from(self.ds.lens[cur as usize])
                {
                    return cur as usize;
                }
                if sseq_pos > start {
                    btm = mid + 1;
                    cur_btm = cur;
                } else {
                    top = cur - 1;
                }
                continue;
            }
            return cur_btm as usize;
        }
        cur_btm as usize
    }

    /// `GetAlnPosFromSeqPos(row, seq_pos, eBackwards, true)` of a plus-strand row.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:527-594
    /// ```c++
    ///     TNumseg raw_seg = GetRawSeg(row, seq_pos, dir, try_reverse_dir);
    ///     if (raw_seg < 0) { // out of seq range
    ///         return -1;
    ///     }
    ///
    ///     TSeqPos start = m_Starts[raw_seg * m_NumRows + row];
    ///     TSeqPos len   = x_GetLen(row, raw_seg);
    ///     TSeqPos stop  = start + len -1;
    ///     bool    plus  = IsPositiveStrand(row);
    ///
    ///     CNumSegWithOffset seg = x_GetSegFromRawSeg(raw_seg);
    ///
    ///     if (dir == eNone) {
    ///         ...
    ///     } else {
    ///         // check if within unaligned region
    ///         // if seq_pos outside the segment returned by GetRawSeg
    ///         // then return the edge alnpos
    ///         if ((plus ? seq_pos < start : seq_pos > stop)) {
    ///             return GetAlnStart(seg.GetAlnSeg());
    ///         }
    ///         if ((plus ? seq_pos > stop : seq_pos < start)) {
    ///             return GetAlnStop(seg.GetAlnSeg());
    ///         }
    ///
    ///
    ///         // check if within an insert
    ///         if (seg.GetOffset()  &&
    ///             (dir == eRight  ||
    ///              dir == (plus ? eForward : eBackwards))) {
    ///             ...
    ///         }
    ///         if (seg.GetOffset()  &&
    ///             (dir == eLeft  ||
    ///              dir == (plus ? eBackwards : eForward))) {
    ///
    ///             // seek the nearest alnpos on left
    ///             if (seg.GetAlnSeg() >= 0) {
    ///                 return GetAlnStop(seg.GetAlnSeg());
    ///             } else if (try_reverse_dir) {
    ///                 return GetAlnStart(seg.GetAlnSeg() + 1);
    ///             } else {
    ///                 return -1;
    ///             }
    ///
    ///         }
    ///     }
    ///
    ///     // main case: seq_pos is within an alnseg
    ///     //assert(seq_pos >= start  &&  seq_pos <= stop);
    ///     TSeqPos delta = (seq_pos - start) / GetWidth(row);
    ///     return m_AlnStarts[seg.GetAlnSeg()]
    ///         + (plus ? delta : m_Lens[raw_seg] - 1 - delta);
    /// ```
    pub(crate) fn aln_pos_from_seq_pos(&self, row: usize, seq_pos: u32) -> i32 {
        let raw_seg = self.raw_seg_backwards(row, seq_pos);
        let start = self.raw_start(row, raw_seg) as u32;
        let stop = start.wrapping_add(self.ds.lens[raw_seg]).wrapping_sub(1);
        let seg = self.seg_from_raw_seg(raw_seg);
        let edge = |aln_seg: i32, stop: bool| -> i32 {
            let aln_seg = aln_seg as usize;
            (if stop {
                self.seg_aln_stop(aln_seg)
            } else {
                self.aln_start(aln_seg)
            }) as i32
        };
        if seq_pos < start {
            return edge(seg.aln_seg, false);
        }
        if seq_pos > stop {
            return edge(seg.aln_seg, true);
        }
        if seg.offset != 0 {
            return if seg.aln_seg >= 0 {
                edge(seg.aln_seg, true)
            } else {
                edge(seg.aln_seg + 1, false)
            };
        }
        self.aln_starts[seg.aln_seg as usize].wrapping_add(seq_pos - start) as i32
    }

    /// `GetWholeAlnSeqString(row, buffer, &insert_aln_starts, &insert_starts, &insert_lens,
    /// scrn_width, &scrn_lefts, &scrn_rights)` of a plus-strand row of width 1, the letters read
    /// by `letters` (`GetSeqString`). The gap character is `-` and the end character a space
    /// (`x_PrepareRowData`, showalign.cpp:1209-1215).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:199-405
    /// ```c++
    ///     TSeqPos       aln_pos = 0,
    ///         len = 0,
    ///         curr_pos = 0,
    ///         anchor_pos = 0,
    ///         scrn_pos = 0,
    ///         prev_len = 0,
    ///         ttl_len = 0;
    ///     TSignedSeqPos start = -1,
    ///         stop = -1,
    ///         scrn_lft_seq_pos = -1,
    ///         scrn_rgt_seq_pos = -1,
    ///         prev_aln_pos = -1,
    ///         prev_start = -1;
    ///     TNumseg       seg;
    ///     int           pos, nscrns, delta;
    ///
    ///     TSeqPos aln_len = GetAlnStop() + 1;
    ///
    ///     bool anchored = IsSetAnchor();
    ///     bool plus     = IsPositiveStrand(row);
    ///     int  width    = GetWidth(row);
    ///
    ///     scrn_width *= width;
    ///
    ///     const bool record_inserts = insert_starts && insert_lens;
    ///     const bool record_coords  = scrn_width && scrn_lefts && scrn_rights;
    ///
    ///     // allocate space for the row
    ///     buffer.clear();
    ///     buffer.reserve(aln_len);
    ///     string buff;
    ///
    ///     const TNumseg& left_seg = x_GetSeqLeftSeg(row);
    ///     const TNumseg& right_seg = x_GetSeqRightSeg(row);
    ///
    ///     // loop through all segments
    ///     for (seg = 0, pos = row, aln_pos = 0, anchor_pos = m_Anchor;
    ///          seg < m_NumSegs;
    ///          ++seg, pos += m_NumRows, anchor_pos += m_NumRows) {
    ///
    ///         const TSeqPos& seg_len = m_Lens[seg];
    ///         start = m_Starts[pos];
    ///         len = seg_len * width;
    ///
    ///         if (anchored  &&  m_Starts[anchor_pos] < 0) {
    ///             if (start >= 0) {
    ///                 // record the insert if requested
    ///                 if (record_inserts) {
    ///                     if (prev_aln_pos == (TSignedSeqPos)(aln_pos / width)  &&
    ///                         start == (TSignedSeqPos)(plus ? prev_start + prev_len :
    ///                                   prev_start - len)) {
    ///                         // consolidate the adjacent inserts
    ///                         ttl_len += len;
    ///                         insert_lens->pop_back();
    ///                         insert_lens->push_back(ttl_len);
    ///                         if (!plus) {
    ///                             insert_starts->pop_back();
    ///                             insert_starts->push_back(start);
    ///                         }
    ///                     } else {
    ///                         prev_aln_pos = aln_pos / width;
    ///                         ttl_len = len;
    ///                         insert_starts->push_back(start);
    ///                         insert_aln_starts->push_back(prev_aln_pos);
    ///                         insert_lens->push_back(len);
    ///                     }
    ///                     prev_start = start;
    ///                     prev_len = len;
    /// 		}
    ///             }
    ///         } else {
    ///             if (start >= 0) {
    ///                 stop = start + len - 1;
    ///
    ///                 // add regular sequence to buffer
    ///                 GetSeqString(buff, row, start, stop);
    ///                 TSeqPos buf_len = min<TSeqPos>(TSeqPos(buff.size()), seg_len);
    ///                 buffer += buff;
    ///                 if (buf_len < seg_len) {
    ///                     // Not enough chars in the sequence, add gap
    ///                     buf_len = seg_len - buf_len;
    ///                     char fill_ch;
    ///
    ///                     if (seg < left_seg  ||  seg > right_seg) {
    ///                         fill_ch = GetEndChar();
    ///                     } else {
    ///                         fill_ch = GetGapChar(row);
    ///                     }
    ///
    ///                     for (size_t i = 0; i < buf_len; ++i) {
    ///                         buffer += fill_ch;
    ///                     }
    ///                 }
    ///
    ///                 // take care of coords if necessary
    ///                 if (record_coords) {
    ///                     if (scrn_lft_seq_pos < 0) {
    ///                         scrn_lft_seq_pos = plus ? start : stop;
    ///                         if (scrn_rgt_seq_pos < 0) {
    ///                             scrn_rgt_seq_pos = scrn_lft_seq_pos;
    ///                         }
    ///                     }
    ///                     // previous scrns
    ///                     nscrns = (aln_pos - scrn_pos) / scrn_width;
    ///                     for (int i = 0; i < nscrns; i++) {
    ///                         scrn_lefts->push_back(scrn_lft_seq_pos);
    ///                         scrn_rights->push_back(scrn_rgt_seq_pos);
    ///                         if (i == 0) {
    ///                             scrn_lft_seq_pos = plus ? start : stop;
    ///                         }
    ///                         scrn_pos += scrn_width;
    ///                     }
    ///                     if (nscrns > 0) {
    ///                         scrn_lft_seq_pos = plus ? start : stop;
    ///                     }
    ///                     // current scrns
    ///                     nscrns = (aln_pos + len - scrn_pos) / scrn_width;
    ///                     curr_pos = aln_pos;
    ///                     for (int i = 0; i < nscrns; i++) {
    ///                         delta = (plus ?
    ///                                  scrn_width - (curr_pos - scrn_pos) :
    ///                                  curr_pos - scrn_pos - scrn_width);
    ///
    ///                         scrn_lefts->push_back(scrn_lft_seq_pos);
    ///                         if (plus ?
    ///                             scrn_lft_seq_pos < start :
    ///                             scrn_lft_seq_pos > stop) {
    ///                             scrn_lft_seq_pos = (plus ? start : stop) +
    ///                                 delta;
    ///                             scrn_rgt_seq_pos = scrn_lft_seq_pos +
    ///                                 (plus ? -1 : 1);
    ///                         } else {
    ///                             scrn_rgt_seq_pos = scrn_lft_seq_pos + (plus ? -1 : 1)
    ///                                 + delta;
    ///                             scrn_lft_seq_pos += delta;
    ///                         }
    ///                         if (seg == left_seg  &&
    ///                             scrn_lft_seq_pos == scrn_rgt_seq_pos) {
    ///                             if (plus) {
    ///                                 scrn_rgt_seq_pos--;
    ///                             } else {
    ///                                 scrn_rgt_seq_pos++;
    ///                             }
    ///                         }
    ///                         scrn_rights->push_back(scrn_rgt_seq_pos);
    ///                         curr_pos = scrn_pos += scrn_width;
    ///                     }
    ///                     if (aln_pos + len <= scrn_pos) {
    ///                         scrn_lft_seq_pos = -1; // reset
    ///                     }
    ///                     scrn_rgt_seq_pos = plus ? stop : start;
    ///                 }
    ///             } else {
    ///                 // add appropriate number of gap/end chars
    ///
    ///                 char fill_ch;
    ///
    ///                 if (seg < left_seg  ||  seg > right_seg) {
    ///                     fill_ch = GetEndChar();
    ///                 } else {
    ///                     fill_ch = GetGapChar(row);
    ///                 }
    ///
    ///                 for (size_t i = 0; i < seg_len; ++i) {
    ///                     buffer += fill_ch;
    ///                 }
    ///             }
    ///             aln_pos += len;
    ///         }
    ///
    ///     }
    ///
    ///     // take care of the remaining coords if necessary
    ///     if (record_coords) {
    ///         // previous scrns
    ///         TSeqPos pos_diff = aln_pos - scrn_pos;
    ///         if (pos_diff > 0) {
    ///             nscrns = pos_diff / scrn_width;
    ///             if (pos_diff % scrn_width) {
    ///                 nscrns++;
    ///             }
    ///             for (int i = 0; i < nscrns; i++) {
    ///                 scrn_lefts->push_back(scrn_lft_seq_pos);
    ///                 scrn_rights->push_back(scrn_rgt_seq_pos);
    ///                 if (i == 0) {
    ///                     scrn_lft_seq_pos = scrn_rgt_seq_pos;
    ///                 }
    ///                 scrn_pos += scrn_width;
    ///             }
    ///         }
    ///     }
    ///     return buffer;
    /// ```
    /// (The arithmetic is that of NCBI's types: `TSeqPos` (`u32`) wraps, `TSignedSeqPos` and
    /// `int` are `i32`, and a quotient stored in an `int` is cut to 32 bits.)
    pub(crate) fn whole_aln_seq_string(
        &self,
        row: usize,
        scrn_width: u32,
        letters: &mut dyn FnMut(SeqKey, u32, u32, &mut Vec<u8>),
    ) -> RowString {
        let mut out = RowString {
            buffer: Vec::new(),
            inserts: std::collections::VecDeque::new(),
            lefts: std::collections::VecDeque::new(),
            rights: std::collections::VecDeque::new(),
        };
        let (mut aln_pos, mut scrn_pos, mut prev_len, mut ttl_len) = (0u32, 0u32, 0u32, 0u32);
        let (mut scrn_lft_seq_pos, mut scrn_rgt_seq_pos) = (-1i32, -1i32);
        let (mut prev_aln_pos, mut prev_start) = (-1i32, -1i32);
        let width = 1u32;
        let scrn_width = scrn_width.wrapping_mul(width);
        let record_coords = scrn_width != 0;
        out.buffer.reserve(self.aln_stop().wrapping_add(1) as usize);
        let (left_seg, right_seg) = (self.left_segs[row], self.right_segs[row]);
        let key = self.ds.ids[row];
        let mut buff = Vec::new();
        for seg in 0..self.num_segs {
            let seg_len = self.ds.lens[seg];
            let start = self.raw_start(row, seg);
            let len = seg_len.wrapping_mul(width);
            if self
                .anchor
                .is_some_and(|anchor| self.raw_start(anchor, seg) < 0)
            {
                if start >= 0 {
                    if prev_aln_pos == (aln_pos / width) as i32
                        && start == prev_start.wrapping_add(prev_len as i32)
                    {
                        ttl_len = ttl_len.wrapping_add(len);
                        out.inserts.back_mut().expect("the insert before").2 = ttl_len;
                    } else {
                        prev_aln_pos = (aln_pos / width) as i32;
                        ttl_len = len;
                        out.inserts
                            .push_back((prev_aln_pos as u32, start as u32, len));
                    }
                    prev_start = start;
                    prev_len = len;
                }
                continue;
            }
            if start >= 0 {
                let stop = start.wrapping_add(len as i32).wrapping_sub(1);
                buff.clear();
                letters(key, start as u32, stop as u32, &mut buff);
                let buf_len = (buff.len() as u32).min(seg_len);
                out.buffer.extend_from_slice(&buff);
                if buf_len < seg_len {
                    let fill_ch = if seg < left_seg || seg > right_seg {
                        b' '
                    } else {
                        b'-'
                    };
                    out.buffer
                        .extend(std::iter::repeat_n(fill_ch, (seg_len - buf_len) as usize));
                }
                if record_coords {
                    if scrn_lft_seq_pos < 0 {
                        scrn_lft_seq_pos = start;
                        if scrn_rgt_seq_pos < 0 {
                            scrn_rgt_seq_pos = scrn_lft_seq_pos;
                        }
                    }
                    let nscrns = (aln_pos.wrapping_sub(scrn_pos) / scrn_width) as i32;
                    for i in 0..nscrns {
                        out.lefts.push_back(scrn_lft_seq_pos as u32);
                        out.rights.push_back(scrn_rgt_seq_pos as u32);
                        if i == 0 {
                            scrn_lft_seq_pos = start;
                        }
                        scrn_pos = scrn_pos.wrapping_add(scrn_width);
                    }
                    if nscrns > 0 {
                        scrn_lft_seq_pos = start;
                    }
                    let nscrns =
                        (aln_pos.wrapping_add(len).wrapping_sub(scrn_pos) / scrn_width) as i32;
                    let mut curr_pos = aln_pos;
                    for _ in 0..nscrns {
                        let delta = scrn_width.wrapping_sub(curr_pos.wrapping_sub(scrn_pos)) as i32;
                        out.lefts.push_back(scrn_lft_seq_pos as u32);
                        if scrn_lft_seq_pos < start {
                            scrn_lft_seq_pos = start.wrapping_add(delta);
                            scrn_rgt_seq_pos = scrn_lft_seq_pos.wrapping_sub(1);
                        } else {
                            scrn_rgt_seq_pos = scrn_lft_seq_pos.wrapping_sub(1).wrapping_add(delta);
                            scrn_lft_seq_pos = scrn_lft_seq_pos.wrapping_add(delta);
                        }
                        if seg == left_seg && scrn_lft_seq_pos == scrn_rgt_seq_pos {
                            scrn_rgt_seq_pos = scrn_rgt_seq_pos.wrapping_sub(1);
                        }
                        out.rights.push_back(scrn_rgt_seq_pos as u32);
                        scrn_pos = scrn_pos.wrapping_add(scrn_width);
                        curr_pos = scrn_pos;
                    }
                    if aln_pos.wrapping_add(len) <= scrn_pos {
                        scrn_lft_seq_pos = -1;
                    }
                    scrn_rgt_seq_pos = stop;
                }
            } else {
                let fill_ch = if seg < left_seg || seg > right_seg {
                    b' '
                } else {
                    b'-'
                };
                out.buffer
                    .extend(std::iter::repeat_n(fill_ch, seg_len as usize));
            }
            aln_pos = aln_pos.wrapping_add(len);
        }
        if record_coords {
            let pos_diff = aln_pos.wrapping_sub(scrn_pos);
            if pos_diff > 0 {
                let mut nscrns = (pos_diff / scrn_width) as i32;
                if pos_diff % scrn_width != 0 {
                    nscrns += 1;
                }
                for i in 0..nscrns {
                    out.lefts.push_back(scrn_lft_seq_pos as u32);
                    out.rights.push_back(scrn_rgt_seq_pos as u32);
                    if i == 0 {
                        scrn_lft_seq_pos = scrn_rgt_seq_pos;
                    }
                    scrn_pos = scrn_pos.wrapping_add(scrn_width);
                }
            }
        }
        out
    }
}
