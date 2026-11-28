//! NCBI's `CAlnMap` and `CAlnVec` (`c++/src/objtools/alnmgr/alnmap.cpp`, `c++/src/objtools/alnmgr/alnvec.cpp`) over the merged
//! Dense-seg of a query's section (`aln_mix::merge`), as the query-anchored reports read it:
//! the anchor (`SetAnchor`), the alignment coordinates of the segments and of each row
//! (`GetAlnStart`, `GetAlnStop`, `GetSeqAlnRange`, `GetSeqStart`, `GetSeqStop`),
//! `GetAlnPosFromSeqPos` with `eBackwards` and `try_reverse_dir` (the masks, `x_FillLocList`)
//! and `GetWholeAlnSeqString` with its insertions and per-line coordinates.
//!
//! A row is on the plus or the minus strand and has a width of 1, or of 3 (a TBLASTN subject's
//! nucleotides, translated: `GetSeqString`, `TranslateNAToAA`); `aln_mix` module doc.

use std::io;

use super::{internal, Denseg, SeqKey, Sequences};

/// NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnmap.hpp:321-325
#[derive(Clone, Copy)]
struct SegWithOffset {
    aln_seg: i32,
    offset: i32,
}

/// The CAlnMap of a Dense-seg (`ds`: a section's merged one, or one HSP's in outfmt 0,
/// `super::HspAlnVec`), with the anchor when set, in the vectors `b` (`MapBuffers`).
pub(crate) struct AlnMap<'d> {
    ds: &'d Denseg,
    num_rows: usize,
    num_segs: usize,
    anchor: Option<usize>,
    b: &'d mut MapBuffers,
}

/// The vectors of an `AlnMap`, which the caller keeps from map to map to reuse them
/// (`AlnMap::new_in`: outfmt 0 makes one map per HSP, `super::HspAlnVec`; the query-anchored
/// section one per group of frames). No NCBI counterpart: NCBI makes a new `CAlnVec` per HSP
/// (showalign.cpp:1857-1860).
#[derive(Default)]
pub(crate) struct MapBuffers {
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

/// What `GetWholeAlnSeqString` makes besides the per-line coordinates when asked for it
/// (`AlnMap::aln_seq_string`): the row's letters (`buffer`, read by `GetSeqString` from
/// `sequences`) and, in an anchored map, its insertions (`insert_aln_starts`, `insert_starts` and
/// `insert_lens`: `RowString::inserts`).
struct Letters<'a> {
    sequences: &'a dyn Sequences,
    buffer: &'a mut Vec<u8>,
    inserts: &'a mut std::collections::VecDeque<(u32, u32, u32)>,
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
    /// `CAlnMap` of `ds` in the vectors `b`: no anchor (`x_CreateAlnStarts`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:49-60
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:716-755 (`x_GetSeqLeftSeg`,
    /// `x_GetSeqRightSeg`: the first and the last segment of the row with a start; "Invalid
    /// Dense-seg: Row ... contains gaps only." for a row of none, which `x_CreateDenseg` removed,
    /// and which an HSP's own Dense-seg does not have).
    /// A Dense-seg without segments is an internal error (the merge of an HSP makes at least one,
    /// as does every producer of an HSP's segments, `report::serial::seqalign`, and `GetAlnStop`
    /// would read `m_AlnStarts[-1]`).
    pub(crate) fn new_in(ds: &'d Denseg, b: &'d mut MapBuffers) -> Result<Self, &'static str> {
        let num_rows = ds.ids.len();
        let num_segs = ds.lens.len();
        if num_segs == 0 {
            return Err("a Dense-seg without segments");
        }
        b.aln_seg_idx.clear();
        b.num_seg_with_offsets.clear();
        b.aln_starts.clear();
        b.aln_starts.reserve(num_segs);
        let mut start = 0u32;
        let mut len = 0u32;
        for &seg_len in &ds.lens {
            start = start.wrapping_add(len);
            b.aln_starts.push(start);
            len = seg_len;
        }
        b.left_segs.clear();
        b.right_segs.clear();
        let starts = &ds.starts[..num_segs * num_rows];
        for row in 0..num_rows {
            let present = |seg: &usize| starts[seg * num_rows + row] >= 0;
            b.left_segs.push((0..num_segs).find(present).unwrap_or(0));
            b.right_segs
                .push((0..num_segs).rev().find(present).unwrap_or(0));
        }
        Ok(Self {
            ds,
            num_rows,
            num_segs,
            anchor: None,
            b,
        })
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:79-121
    /// (The query row, the anchor, has the letters of every query segment; NCBI's throw is an
    /// internal error here.)
    pub(crate) fn set_anchor(&mut self, anchor: usize) -> Result<(), &'static str> {
        self.b.aln_seg_idx.clear();
        self.b.aln_starts.clear();
        self.b.num_seg_with_offsets.clear();
        let (mut start, mut len, mut aln_seg, mut offset) = (0u32, 0u32, -1i32, 0i32);
        self.anchor = Some(anchor);
        for seg in 0..self.num_segs {
            if self.ds.starts[seg * self.num_rows + anchor] != -1 {
                aln_seg += 1;
                offset = 0;
                self.b.aln_seg_idx.push(seg);
                self.b
                    .num_seg_with_offsets
                    .push(SegWithOffset { aln_seg, offset: 0 });
                start = start.wrapping_add(len);
                self.b.aln_starts.push(start);
                len = self.ds.lens[seg];
            } else {
                offset += 1;
                self.b
                    .num_seg_with_offsets
                    .push(SegWithOffset { aln_seg, offset });
            }
        }
        if self.b.aln_seg_idx.is_empty() {
            return Err("no sequence on the anchor row");
        }
        Ok(())
    }

    pub(crate) fn num_rows(&self) -> usize {
        self.num_rows
    }

    pub(crate) fn seq_key(&self, row: usize) -> SeqKey {
        self.ds.ids[row]
    }

    // NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnmap.hpp:559-581,599-603
    /// `GetWidth`.
    pub(crate) fn width(&self, row: usize) -> u32 {
        self.ds.widths.get(row).copied().unwrap_or(1)
    }

    /// `IsPositiveStrand`.
    pub(crate) fn is_positive_strand(&self, row: usize) -> bool {
        self.ds.plus[row]
    }

    /// `x_GetLen`: segment `seg`'s length in the row's letters.
    fn x_get_len(&self, row: usize, seg: usize) -> u32 {
        self.ds.lens[seg].wrapping_mul(self.width(row))
    }

    /// `GetNumSegs`.
    fn aln_num_segs(&self) -> usize {
        if self.anchor.is_some() {
            self.b.aln_seg_idx.len()
        } else {
            self.num_segs
        }
    }

    /// `x_GetRawSegFromSeg`.
    fn raw_seg_from_seg(&self, seg: usize) -> usize {
        if self.anchor.is_some() {
            self.b.aln_seg_idx[seg]
        } else {
            seg
        }
    }

    /// `x_GetSegFromRawSeg`.
    fn seg_from_raw_seg(&self, raw_seg: usize) -> SegWithOffset {
        if self.anchor.is_some() {
            self.b.num_seg_with_offsets[raw_seg]
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
    fn aln_start(&self, seg: usize) -> u32 {
        self.b.aln_starts[seg]
    }

    fn seg_aln_stop(&self, seg: usize) -> u32 {
        self.b.aln_starts[seg]
            .wrapping_add(self.ds.lens[self.raw_seg_from_seg(seg)])
            .wrapping_sub(1)
    }

    /// `GetAlnStop()`. There is a segment: `new` and `set_anchor` refuse a map of none.
    #[inline(always)]
    pub(crate) fn aln_stop(&self) -> u32 {
        self.seg_aln_stop(self.aln_num_segs() - 1)
    }

    // NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnmap.hpp:665-681
    /// `GetSeqStart` (the row's lowest position).
    pub(crate) fn seq_start(&self, row: usize) -> u32 {
        let seg = if self.is_positive_strand(row) {
            self.b.left_segs[row]
        } else {
            self.b.right_segs[row]
        };
        self.raw_start(row, seg) as u32
    }

    /// `GetSeqStop` (the row's highest position).
    pub(crate) fn seq_stop(&self, row: usize) -> u32 {
        let seg = if self.is_positive_strand(row) {
            self.b.right_segs[row]
        } else {
            self.b.left_segs[row]
        };
        (self.raw_start(row, seg) as u32)
            .wrapping_add(self.x_get_len(row, seg))
            .wrapping_sub(1)
    }

    /// `GetSeqAlnRange`: the row's first and last alignment position.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:969-998
    #[inline(always)]
    pub(crate) fn seq_aln_range(&self, row: usize) -> (i64, i64) {
        if self.anchor.is_some() {
            let present = |seg: &usize| self.raw_start(row, self.b.aln_seg_idx[*seg]) >= 0;
            let start = (0..self.b.aln_seg_idx.len())
                .find(present)
                .map_or(-1, |seg| i64::from(self.aln_start(seg) as i32));
            let stop = (0..self.b.aln_seg_idx.len())
                .rev()
                .find(present)
                .map_or(-1, |seg| i64::from(self.seg_aln_stop(seg) as i32));
            (start, stop)
        } else {
            (
                i64::from(self.aln_start(self.b.left_segs[row]) as i32),
                i64::from(self.seg_aln_stop(self.b.right_segs[row]) as i32),
            )
        }
    }

    /// `GetRawSeg` of `row` and `seq_pos` with `eBackwards` and true: on a minus-strand row the segments are searched
    /// from the last (`last - cur`), whose positions are the lowest.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:400-524
    fn raw_seg_backwards(&self, row: usize, seq_pos: u32) -> usize {
        let sseq_pos = seq_pos as i32;
        let plus = self.is_positive_strand(row);
        let present = |seg: i64| self.raw_start(row, seg as usize) >= 0;
        let num_segs = self.num_segs as i64;
        let last = num_segs - 1;
        let index = |cur: i64| if plus { cur } else { last - cur };
        if sseq_pos < self.seq_start(row) as i32 {
            let found = if plus {
                (0..num_segs).find(|&seg| present(seg))
            } else {
                (0..num_segs).rev().find(|&seg| present(seg))
            };
            if let Some(seg) = found {
                return seg as usize;
            }
        } else if sseq_pos > self.seq_stop(row) as i32 {
            let found = if plus {
                (0..num_segs).rev().find(|&seg| present(seg))
            } else {
                (0..num_segs).find(|&seg| present(seg))
            };
            if let Some(seg) = found {
                return seg as usize;
            }
        }
        let (mut btm, mut top, mut cur_btm) = (0i64, last, 0i64);
        while btm <= top {
            let mid = (top + btm) / 2;
            let mut cur = mid;
            let mut start = -1;
            while cur <= top && {
                start = self.raw_start(row, index(cur) as usize);
                start < 0
            } {
                cur += 1;
            }
            if cur <= top && start >= 0 {
                let len = self.x_get_len(row, index(cur) as usize);
                if sseq_pos >= start && i64::from(seq_pos) < i64::from(start) + i64::from(len) {
                    return index(cur) as usize;
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
                start = self.raw_start(row, index(cur) as usize);
                start < 0
            } {
                cur -= 1;
            }
            if cur >= btm && start >= 0 {
                let len = self.x_get_len(row, index(cur) as usize);
                if sseq_pos >= start && i64::from(seq_pos) < i64::from(start) + i64::from(len) {
                    return index(cur) as usize;
                }
                if sseq_pos > start {
                    btm = mid + 1;
                    cur_btm = cur;
                } else {
                    top = cur - 1;
                }
                continue;
            }
            return index(cur_btm) as usize;
        }
        index(cur_btm) as usize
    }

    /// `GetAlnPosFromSeqPos` of `row` and `seq_pos` with `eBackwards` and true. Within an insertion a plus-strand
    /// row takes the column on its left, a minus-strand row (`eBackwards` is its `eForward`) the
    /// column on its right.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:527-596
    pub(crate) fn aln_pos_from_seq_pos(&self, row: usize, seq_pos: u32) -> i32 {
        let raw_seg = self.raw_seg_backwards(row, seq_pos);
        let start = self.raw_start(row, raw_seg) as u32;
        let stop = start
            .wrapping_add(self.x_get_len(row, raw_seg))
            .wrapping_sub(1);
        let plus = self.is_positive_strand(row);
        let seg = self.seg_from_raw_seg(raw_seg);
        let edge = |aln_seg: i32, stop: bool| -> i32 {
            let aln_seg = aln_seg as usize;
            (if stop {
                self.seg_aln_stop(aln_seg)
            } else {
                self.aln_start(aln_seg)
            }) as i32
        };
        if if plus {
            seq_pos < start
        } else {
            seq_pos > stop
        } {
            return edge(seg.aln_seg, false);
        }
        if if plus {
            seq_pos > stop
        } else {
            seq_pos < start
        } {
            return edge(seg.aln_seg, true);
        }
        if seg.offset != 0 && !plus {
            return if (seg.aln_seg as i64) < self.aln_num_segs() as i64 - 1 {
                edge(seg.aln_seg + 1, false)
            } else {
                edge(seg.aln_seg, true)
            };
        }
        if seg.offset != 0 && plus {
            return if seg.aln_seg >= 0 {
                edge(seg.aln_seg, true)
            } else {
                edge(seg.aln_seg + 1, false)
            };
        }
        let delta = (seq_pos - start) / self.width(row);
        let delta = if plus {
            delta
        } else {
            self.ds.lens[raw_seg] - 1 - delta
        };
        self.b.aln_starts[seg.aln_seg as usize].wrapping_add(delta) as i32
    }

    /// `GetAlnPosFromSeqPos` of `row` and `seq_pos` with its defaults (`eNone`): the alignment
    /// position of the row's letter `seq_pos`, -1 where no segment of the row holds it (outside
    /// the row's letters, or between two of its segments). A map without an anchor has no
    /// insertions, so the letter's segment is its alignment segment.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:527-541
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:589-596
    /// `GetRawSeg` with `eNone` (alnmap.cpp:400-524) runs the search of `raw_seg_backwards` and
    /// returns -1 where that search gives a segment that does not hold `seq_pos` (out of the
    /// row's range, or the closest segment of an unaligned letter): the segment it gives where one
    /// holds the letter is the same, as both return from the same loop there.
    pub(crate) fn aln_pos_from_seq_pos_none(&self, row: usize, seq_pos: u32) -> i32 {
        debug_assert!(self.anchor.is_none(), "a map without an anchor");
        let sseq_pos = seq_pos as i32;
        if sseq_pos < self.seq_start(row) as i32 || sseq_pos > self.seq_stop(row) as i32 {
            return -1;
        }
        let raw_seg = self.raw_seg_backwards(row, seq_pos);
        let start = self.raw_start(row, raw_seg);
        let holds = start >= 0
            && sseq_pos >= start
            && i64::from(seq_pos) < i64::from(start) + i64::from(self.x_get_len(row, raw_seg));
        if !holds {
            return -1;
        }
        let delta = (seq_pos - start as u32) / self.width(row);
        let delta = if self.is_positive_strand(row) {
            delta
        } else {
            self.ds.lens[raw_seg].wrapping_sub(1).wrapping_sub(delta)
        };
        self.b.aln_starts[raw_seg].wrapping_add(delta) as i32
    }

    /// `GetSeg`: the alignment segment of the alignment position `aln_pos`, which is at most
    /// `aln_stop` (a map without an anchor: the raw segment).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:373-396
    fn seg_of(&self, aln_pos: u32) -> usize {
        let (mut btm, mut top) = (0usize, self.b.aln_starts.len() - 1);
        while btm < top {
            let mid = (top + btm) / 2;
            if self.b.aln_starts[mid] as i32 == aln_pos as i32 {
                return mid;
            }
            if self.b.aln_starts[mid + 1] as i32 <= aln_pos as i32 {
                btm = mid + 1;
            } else {
                top = mid;
            }
        }
        top
    }

    /// `GetSeqPosFromAlnPos` of `row` and `aln_pos` with its defaults (`eNone`): the row's letter
    /// in the column `aln_pos` (a position past the alignment's end taken as its last column),
    /// -1 where the row has a gap there.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:663-686
    pub(crate) fn seq_pos_from_aln_pos(&self, row: usize, aln_pos: u32) -> i32 {
        debug_assert!(self.anchor.is_none(), "a map without an anchor");
        let aln_pos = aln_pos.min(self.aln_stop());
        let seg = self.seg_of(aln_pos);
        let pos = self.raw_start(row, seg);
        if pos < 0 {
            return -1;
        }
        let delta = aln_pos
            .wrapping_sub(self.aln_start(seg))
            .wrapping_mul(self.width(row));
        let pos = pos as u32;
        (if self.is_positive_strand(row) {
            pos.wrapping_add(delta)
        } else {
            pos.wrapping_add(self.x_get_len(row, seg))
                .wrapping_sub(1)
                .wrapping_sub(delta)
        }) as i32
    }

    /// `CAlignFormatUtil::GetAlignLengths` of the map (no anchor): the alignment's length (the
    /// columns of row 0's chunks), its gap letters and its gap openings (the gap chunks of every
    /// row), over `GetAlnChunks` of the row in row 0's `GetSeqAlnRange`, with the default flags
    /// (`fAlnSegsOnly`). A chunk is a run of adjacent segments that all hold letters of the row
    /// or all hold none, cut where its first segment has unaligned letters of the row on its
    /// right or the next segment has some on its left (`x_SetRawSegTypes`); without an anchor no
    /// segment is skipped, and row 0's range starts and ends at whole segments.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1098-1121
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:1000-1031
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:1081-1208
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:234-367
    pub(crate) fn align_lengths(&self) -> (i64, i64, i64) {
        debug_assert!(self.anchor.is_none(), "a map without an anchor");
        const SEQ: u8 = 1;
        const UNALIGNED_LEFT: u8 = 2;
        const UNALIGNED_RIGHT: u8 = 4;
        let (left, right) = (self.b.left_segs[0], self.b.right_segs[0]);
        let (mut length, mut gaps, mut gap_opens) = (0i64, 0i64, 0i64);
        let mut types = vec![0u8; self.num_segs];
        for row in 0..self.num_rows {
            let plus = self.is_positive_strand(row);
            types.fill(0);
            let mut left_seq_pos: i32 = -1;
            for (seg, flags) in types.iter_mut().enumerate() {
                let start = self.raw_start(row, seg);
                if start >= 0 {
                    *flags |= SEQ;
                    let len = self.x_get_len(row, seg);
                    if left_seq_pos > 0
                        && if plus {
                            start > left_seq_pos
                        } else {
                            (start as u32).wrapping_add(len) < left_seq_pos as u32
                        }
                    {
                        *flags |= UNALIGNED_LEFT;
                    }
                    left_seq_pos = if plus {
                        (start as u32).wrapping_add(len) as i32
                    } else {
                        start
                    };
                }
            }
            let mut right_seq_pos: i32 = -1;
            for (seg, flags) in types.iter_mut().enumerate().rev() {
                let start = self.raw_start(row, seg);
                if start >= 0 {
                    let len = self.x_get_len(row, seg);
                    if right_seq_pos > 0
                        && if plus {
                            (start as u32).wrapping_add(len) < right_seq_pos as u32
                        } else {
                            start > right_seq_pos
                        }
                    {
                        *flags |= UNALIGNED_RIGHT;
                    }
                    right_seq_pos = if plus {
                        start
                    } else {
                        (start as u32).wrapping_add(len) as i32
                    };
                }
            }
            let mut seg = left;
            while seg <= right {
                let first = types[seg];
                let mut stop = seg;
                while stop < right {
                    let next = types[stop + 1];
                    if (first & SEQ) != (next & SEQ)
                        || first & UNALIGNED_RIGHT != 0
                        || next & UNALIGNED_LEFT != 0
                    {
                        break;
                    }
                    stop += 1;
                }
                let chunk = i64::from(self.seg_aln_stop(stop) as i32)
                    - i64::from(self.aln_start(seg) as i32)
                    + 1;
                let chunk = chunk.max(0);
                if first & SEQ == 0 {
                    gap_opens += 1;
                    gaps += chunk;
                }
                if row == 0 {
                    length += chunk;
                }
                seg = stop + 1;
            }
        }
        (length, gaps, gap_opens)
    }

    /// `GetWholeAlnSeqString` of `row` and `buffer`, with the insert alignment starts, starts and
    /// lengths, `scrn_width`, and the screen lefts and rights (`whole_aln_seq_string`, outfmt 1-4), the
    /// letters read by `seq_string` (`GetSeqString`), or, without `letters`, its per-line
    /// coordinates alone (`screen_coords`, outfmt 0): the one owner of the coordinates of a
    /// row's lines. The gap character is `-` and the end character a space (`x_PrepareRowData`,
    /// showalign.cpp:1209-1215). `emit(left, right)` gets each screen's entries of `scrn_lefts`
    /// and `scrn_rights`, in order.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:199-405
    /// (The arithmetic is that of NCBI's types: `TSeqPos` (`u32`) wraps, `TSignedSeqPos` and
    /// `int` are `i32`, and a quotient stored in an `int` is cut to 32 bits.) A row of width 3 counts
    /// its alignment positions and screens in its nucleotides (`scrn_width` multiplied by `width`). Without
    /// `letters` the insert pointers are null (`record_inserts` is false) and no letter is read:
    /// the coordinates do not depend on them. Inlined into its two callers, so that the one without
    /// letters (`screen_coords`) has no code of theirs.
    #[inline(always)]
    fn aln_seq_string(
        &self,
        row: usize,
        scrn_width: u32,
        mut letters: Option<Letters<'_>>,
        mut emit: impl FnMut(u32, u32),
    ) -> io::Result<()> {
        let (mut aln_pos, mut scrn_pos, mut prev_len, mut ttl_len) = (0u32, 0u32, 0u32, 0u32);
        let (mut scrn_lft_seq_pos, mut scrn_rgt_seq_pos) = (-1i32, -1i32);
        let (mut prev_aln_pos, mut prev_start) = (-1i32, -1i32);
        let plus = self.is_positive_strand(row);
        let width = self.width(row);
        let scrn_width = scrn_width.wrapping_mul(width);
        let record_coords = scrn_width != 0;
        if let Some(letters) = letters.as_mut() {
            letters.buffer.clear();
            letters
                .buffer
                .reserve(self.aln_stop().wrapping_add(1) as usize);
        }
        let (left_seg, right_seg) = (self.b.left_segs[row], self.b.right_segs[row]);
        let mut buff = Vec::new();
        for seg in 0..self.num_segs {
            let seg_len = self.ds.lens[seg];
            let start = self.raw_start(row, seg);
            let len = seg_len.wrapping_mul(width);
            if self
                .anchor
                .is_some_and(|anchor| self.raw_start(anchor, seg) < 0)
            {
                if let (true, Some(letters)) = (start >= 0, letters.as_mut()) {
                    let adjacent = if plus {
                        prev_start.wrapping_add(prev_len as i32)
                    } else {
                        prev_start.wrapping_sub(len as i32)
                    };
                    if prev_aln_pos == (aln_pos / width) as i32 && start == adjacent {
                        ttl_len = ttl_len.wrapping_add(len);
                        let last = letters.inserts.back_mut().expect("the insert before");
                        last.2 = ttl_len;
                        if !plus {
                            last.1 = start as u32;
                        }
                    } else {
                        prev_aln_pos = (aln_pos / width) as i32;
                        ttl_len = len;
                        letters
                            .inserts
                            .push_back((prev_aln_pos as u32, start as u32, len));
                    }
                    prev_start = start;
                    prev_len = len;
                }
                continue;
            }
            let fill_ch = if seg < left_seg || seg > right_seg {
                b' '
            } else {
                b'-'
            };
            if start >= 0 {
                let stop = start.wrapping_add(len as i32).wrapping_sub(1);
                if let Some(letters) = letters.as_mut() {
                    buff.clear();
                    self.seq_string(row, start as u32, stop as u32, letters.sequences, &mut buff)?;
                    let buf_len = (buff.len() as u32).min(seg_len);
                    letters.buffer.extend_from_slice(&buff);
                    if buf_len < seg_len {
                        letters
                            .buffer
                            .extend(std::iter::repeat_n(fill_ch, (seg_len - buf_len) as usize));
                    }
                }
                if record_coords {
                    let near = if plus { start } else { stop };
                    if scrn_lft_seq_pos < 0 {
                        scrn_lft_seq_pos = near;
                        if scrn_rgt_seq_pos < 0 {
                            scrn_rgt_seq_pos = scrn_lft_seq_pos;
                        }
                    }
                    let nscrns = (aln_pos.wrapping_sub(scrn_pos) / scrn_width) as i32;
                    for i in 0..nscrns {
                        emit(scrn_lft_seq_pos as u32, scrn_rgt_seq_pos as u32);
                        if i == 0 {
                            scrn_lft_seq_pos = near;
                        }
                        scrn_pos = scrn_pos.wrapping_add(scrn_width);
                    }
                    if nscrns > 0 {
                        scrn_lft_seq_pos = near;
                    }
                    let nscrns =
                        (aln_pos.wrapping_add(len).wrapping_sub(scrn_pos) / scrn_width) as i32;
                    let step = if plus { -1 } else { 1 };
                    let mut curr_pos = aln_pos;
                    for _ in 0..nscrns {
                        let delta = if plus {
                            scrn_width.wrapping_sub(curr_pos.wrapping_sub(scrn_pos))
                        } else {
                            curr_pos.wrapping_sub(scrn_pos).wrapping_sub(scrn_width)
                        } as i32;
                        let left = scrn_lft_seq_pos;
                        if if plus {
                            scrn_lft_seq_pos < start
                        } else {
                            scrn_lft_seq_pos > stop
                        } {
                            scrn_lft_seq_pos = near.wrapping_add(delta);
                            scrn_rgt_seq_pos = scrn_lft_seq_pos.wrapping_add(step);
                        } else {
                            scrn_rgt_seq_pos =
                                scrn_lft_seq_pos.wrapping_add(step).wrapping_add(delta);
                            scrn_lft_seq_pos = scrn_lft_seq_pos.wrapping_add(delta);
                        }
                        if seg == left_seg && scrn_lft_seq_pos == scrn_rgt_seq_pos {
                            scrn_rgt_seq_pos = scrn_rgt_seq_pos.wrapping_add(step);
                        }
                        emit(left as u32, scrn_rgt_seq_pos as u32);
                        scrn_pos = scrn_pos.wrapping_add(scrn_width);
                        curr_pos = scrn_pos;
                    }
                    if aln_pos.wrapping_add(len) <= scrn_pos {
                        scrn_lft_seq_pos = -1;
                    }
                    scrn_rgt_seq_pos = if plus { stop } else { start };
                }
            } else if let Some(letters) = letters.as_mut() {
                letters
                    .buffer
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
                    emit(scrn_lft_seq_pos as u32, scrn_rgt_seq_pos as u32);
                    if i == 0 {
                        scrn_lft_seq_pos = scrn_rgt_seq_pos;
                    }
                    scrn_pos = scrn_pos.wrapping_add(scrn_width);
                }
            }
        }
        Ok(())
    }

    /// `GetWholeAlnSeqString` of the query-anchored reports (`aln_seq_string` with the letters
    /// and the insertions): the row's `RowString`.
    pub(crate) fn whole_aln_seq_string(
        &self,
        row: usize,
        scrn_width: u32,
        sequences: &dyn Sequences,
    ) -> io::Result<RowString> {
        let mut out = RowString {
            buffer: Vec::new(),
            inserts: std::collections::VecDeque::new(),
            lefts: std::collections::VecDeque::new(),
            rights: std::collections::VecDeque::new(),
        };
        let letters = Letters {
            sequences,
            buffer: &mut out.buffer,
            inserts: &mut out.inserts,
        };
        let (lefts, rights) = (&mut out.lefts, &mut out.rights);
        self.aln_seq_string(row, scrn_width, Some(letters), |left, right| {
            lefts.push_back(left);
            rights.push_back(right);
        })?;
        Ok(out)
    }

    /// The reads of `GetWholeAlnSeqString` of row `row` of a map without an anchor (the tabular
    /// formats'): `GetSeqString` of each segment that the row is in, in order, from its start to
    /// its start plus its length less 1 (`seq_from - 1` for a segment of 0 letters, in the
    /// unsigned `TSeqPos`), given to `read` (`aln_seq_string`'s calls of `seq_string`).
    pub(crate) fn seq_string_reads<E>(
        &self,
        row: usize,
        mut read: impl FnMut(u32, u32) -> Result<(), E>,
    ) -> Result<(), E> {
        debug_assert!(self.anchor.is_none(), "a map without an anchor");
        let width = self.width(row);
        for seg in 0..self.num_segs {
            let start = self.raw_start(row, seg);
            if start >= 0 {
                let len = self.ds.lens[seg].wrapping_mul(width);
                let stop = start.wrapping_add(len as i32).wrapping_sub(1);
                read(start as u32, stop as u32)?;
            }
        }
        Ok(())
    }

    /// The per-line coordinates of row `row` alone (`aln_seq_string` without the letters): each
    /// screen's `scrn_lefts` and `scrn_rights` entries of `GetWholeAlnSeqString`, given in order to
    /// `emit(left, right)`. Outfmt 0 reads them (`x_PrepareRowData` calls `GetWholeAlnSeqString`
    /// with the same `m_LineLen` cast to `int`, showalign.cpp:1281-1283; `super::HspAlnVec`).
    pub(crate) fn screen_coords(
        &self,
        row: usize,
        scrn_width: u32,
        emit: impl FnMut(u32, u32),
    ) -> io::Result<()> {
        self.aln_seq_string(row, scrn_width, None, emit)
    }

    /// `GetSeqString` of `buffer`, `row`, `seq_from` and `seq_to`: the letters `seq_from..=seq_to` of the row's
    /// sequence (`Sequences::letters`, its plus strand), reverse-complemented on a minus-strand row
    /// (the `eStrand_Minus` sequence vector, `x_GetSeqVector`), translated on a row of width 3
    /// with the row's genetic code (`Sequences::genetic_code`; `SetGenCode`, showalign.cpp:
    /// 2100-2101).
    ///
    /// NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnvec.hpp:287-313
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:893-904
    /// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:554-557
    /// (A row's positions lie in its sequence and a row of width 3 has whole codons: the
    /// `TranslateNAToAA` throw is an internal error here.) A segment of length 0 (an HSP of length
    /// 0, `serial::seqalign::empty_row_segments`) reads `seq_from` to `seq_from - 1`, whose
    /// `GetSeqData` interval is empty (`seq_to + 1` is `seq_from` in the unsigned `TSeqPos`, also
    /// at position 0): no letter, and still a translation, which `GetTransTable` makes of the
    /// row's genetic code (`Sequences::genetic_code`, row 0's or the database's): without one
    /// (`None`, BLASTX's code 0) it throws (`UnknownGeneticCode`).
    pub(crate) fn seq_string(
        &self,
        row: usize,
        seq_from: u32,
        seq_to: u32,
        sequences: &dyn Sequences,
        out: &mut Vec<u8>,
    ) -> io::Result<()> {
        let key = self.ds.ids[row];
        let plus = self.is_positive_strand(row);
        let empty = seq_to.wrapping_add(1) == seq_from;
        if self.width(row) == 3 {
            let mut nucleotides = Vec::new();
            if !empty {
                sequences.letters(key, seq_from, seq_to, &mut nucleotides);
            }
            if nucleotides.len() % 3 != 0 {
                return Err(io::Error::other(TranslateNaToAa));
            }
            let Some(code) = sequences.genetic_code(row == 0) else {
                return Err(if row == 0 {
                    internal("a translated query row without a genetic code")
                } else {
                    io::Error::other(UnknownGeneticCode)
                });
            };
            code.translate(&nucleotides, !plus, out);
        } else if !empty {
            let first = out.len();
            sequences.letters(key, seq_from, seq_to, out);
            if !plus {
                crate::algorithm::blastn::lookup::reverse_complement_in_place(&mut out[first..]);
            }
        }
        Ok(())
    }
}

/// NCBI's `Unable to find genetic code number 0` (`CGen_code_table_imp::GetTransTable`), which a
/// translation of a BLASTX row after row 0 throws (the database's code 0,
/// `Sequences::genetic_code`): the error of `AlnMap::seq_string`, which `super::write_section`
/// catches as `DisplaySeqalign`'s `try` does.
#[derive(Debug)]
pub(crate) struct UnknownGeneticCode;

impl UnknownGeneticCode {
    /// Whether `error` is this one.
    pub(crate) fn is(error: &io::Error) -> bool {
        error.get_ref().is_some_and(|inner| inner.is::<Self>())
    }
}

impl std::fmt::Display for UnknownGeneticCode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Unable to find genetic code number 0")
    }
}

impl std::error::Error for UnknownGeneticCode {}

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:895-900
/// NCBI's `CAlnVec::TranslateNAToAA(): NA size expected to be divisible by 3`, thrown when a
/// translated row's clamped nucleotide slice is not a multiple of 3 (TBLASTN `-in_pssm
/// -window_size 0` alignment starting on a frame's NULLB when `rem in {1, 2}`): caught by
/// `super::write_section` as `DisplaySeqalign`'s `try` does.
#[derive(Debug)]
pub(crate) struct TranslateNaToAa;

impl TranslateNaToAa {
    /// Whether `error` is this one.
    pub(crate) fn is(error: &io::Error) -> bool {
        error.get_ref().is_some_and(|inner| inner.is::<Self>())
    }
}

impl std::fmt::Display for TranslateNaToAa {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("CAlnVec::TranslateNAToAA(): NA size expected to be divisible by 3")
    }
}

impl std::error::Error for TranslateNaToAa {}
