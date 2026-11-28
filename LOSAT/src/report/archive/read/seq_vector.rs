//! The sequence map of a request's Bioseq (`CSeqMap`) and the reads of a row of it by the
//! tabular formats' `CAlnVec` (`RowVector::read`, `GetSeqString`): the moves of the one
//! `CSeqVector_CI` of the row's `CSeqVector`, the fills of its cache of 1024 letters and their
//! checks, so that a Bioseq whose letters do not all read (data shorter than its length, segments
//! that end before it, a literal of `kInvalidSeqPos` letters, a coding that does not convert)
//! throws NCBI's exception at the row and the read where NCBI's does (probes `probe3.py`,
//! `probe4.py` of agent E3, 2026-11-27: `reference out of range of Seq-inst data: 1124` for a row
//! from letter 101, 2048 for a minus-strand row over the first 2048 letters, `cannot get
//! seq-data in range`, `cannot locate segment at`, `Invalid sequence length`). The letters of a
//! read that does not throw are the map's (`SeqMap::letters`, the plus strand), which
//! `aln_vec::AlnMap::seq_string` reverse-complements.
//!
//! The model holds what decides an outcome: the segments' positions as the map resolves them
//! (none past a literal of `kInvalidSeqPos` letters, whose length every resolution throws), the
//! segment iterator over the whole map on the vector's strand (`Seg`, `CSeqMap_CI`), the cache's
//! and the backup's positions and sizes, not their letters. The preload of the scanned range
//! (`x_CheckForward`, `x_CheckBackward`) calls `CanGetRange`, which catches every exception and
//! changes nothing that a later read depends on: it is left out.

use super::{
    coding, data_bytes, seq_length, sequence_map, Bioseq, NcbiException, Piece, RenderError,
    ThrowSite, INCOMPATIBLE_CODINGS, INVALID_LENGTH, INVALID_SEQ_POS, IUPACAA, IUPACNA,
    KCACHE_SIZE, MOL_AA, NA2_TO_IUPACNA, NA4_TO_IUPACNA, NCBI2NA, NCBI4NA, NCBIEAA, NCBISTDAA,
    NO_SEGMENT, OUT_OF_RANGE, REPR_MAP, REPR_VIRTUAL, V,
};
use crate::sequence::stdaa_letters::protein_fasta_letters;

/// `kCacheSize` as a `TSeqPos`.
const CACHE: u32 = KCACHE_SIZE as u32;

// NCBI reference (598d8ae6): c++/src/objmgr/seq_map.cpp:202-291
// NCBI reference (598d8ae6): c++/src/objmgr/seq_map.cpp:343-372
// NCBI reference (598d8ae6): c++/src/objmgr/seq_map.cpp:375-432
/// The sequence map of a Bioseq (`CSeqMap` of its `inst`): its segments
/// (`queries::sequence_map`) and their positions as the map resolves them, each from the one
/// before.
pub(crate) struct SeqMap<'a> {
    pieces: Vec<Piece<'a>>,
    /// `starts[i]`: the position of piece `i`, for each piece up to the first of
    /// `kInvalidSeqPos` letters (`first_invalid`, whose length every resolution throws: `Invalid
    /// sequence length`); without one, one more: the end of the last piece.
    starts: Vec<u32>,
    first_invalid: Option<usize>,
    /// `m_SeqLength` as the constructor sets it: the Bioseq's length, `kInvalidSeqPos` where it
    /// declares 4294967295 (`GetLength` then resolves every segment).
    declared: u32,
    /// Whether the Bioseq is a protein (`CAlnVec`'s vector reads it in ncbieaa, else IUPACna).
    protein: bool,
}

impl<'a> SeqMap<'a> {
    /// The map of `bioseq`, or the constructor's throw. A Bioseq without its length is LOSAT's
    /// refusal (NCBI's reader of the archive keeps none); so are segments whose positions reach
    /// 2^32 - 1 (`Sequence position overflow`, which LOSAT does not reproduce).
    pub(crate) fn of(bioseq: Bioseq<'a>) -> Result<Self, RenderError> {
        let declared = match bioseq.length() {
            Some(length) => seq_length(length)?,
            None => {
                return Err(RenderError::Unsupported(
                    "a Bioseq without its length".into(),
                ))
            }
        };
        let inst = bioseq.value().field("inst");
        // A Bioseq without data of the length `kInvalidSeqPos`: a `virtual` one is a gap of 0
        // letters, a `map` one throws (seq_map.cpp:244-248, 258-271).
        let unset = declared == INVALID_SEQ_POS && inst.get("seq-data").is_none();
        let ext = inst.get("ext").map(|ext| ext.choice().0);
        let (pieces, declared) = match (unset, ext, inst.field("repr").int()) {
            (true, Some("map"), _) => {
                return Err(map_throw(
                    246,
                    "CSeq_inst.ext.map with unset CSeq_inst.length",
                ))
            }
            (true, None, REPR_VIRTUAL) => (vec![Piece::Gap(0)], 0),
            (true, None, REPR_MAP) => {
                return Err(map_throw(
                    269,
                    "CSeq_inst.repr map with unset CSeq_inst.length",
                ))
            }
            _ => (sequence_map(inst, declared)?, declared),
        };
        let mut starts = Vec::with_capacity(pieces.len() + 1);
        let mut first_invalid = None;
        let mut at = 0u32;
        for (i, piece) in pieces.iter().enumerate() {
            starts.push(at);
            let length = piece.length();
            if length == INVALID_SEQ_POS {
                first_invalid = Some(i);
                break;
            }
            at = match at.checked_add(length) {
                Some(end) if end != INVALID_SEQ_POS => end,
                _ => {
                    return Err(RenderError::Unsupported(
                        "a Bioseq whose segments pass 2^32 - 2 letters".into(),
                    ))
                }
            };
        }
        if first_invalid.is_none() {
            starts.push(at);
        }
        Ok(SeqMap {
            pieces,
            starts,
            first_invalid,
            declared,
            protein: bioseq.molecule() == MOL_AA,
        })
    }

    /// `GetLength`: the declared length, else the end of the last segment, whose resolution
    /// throws at a segment of `kInvalidSeqPos` letters.
    pub(crate) fn length(&self) -> Result<u32, RenderError> {
        if self.declared != INVALID_SEQ_POS {
            Ok(self.declared)
        } else if self.first_invalid.is_some() {
            Err(invalid_length())
        } else {
            Ok(self.end())
        }
    }

    /// The end of the last segment (no segment of `kInvalidSeqPos` letters).
    fn end(&self) -> u32 {
        *self
            .starts
            .last()
            .expect("a start for each piece and the end")
    }

    /// `x_FindSegment`: the piece that holds `pos`, `None` past the last; the resolution of the
    /// lengths up to it throws at a segment of `kInvalidSeqPos` letters at `pos` or before.
    fn find(&self, pos: u32) -> Result<Option<usize>, RenderError> {
        if let Some(k) = self.first_invalid {
            if self.starts[k] <= pos {
                return Err(invalid_length());
            }
        }
        let known = self.first_invalid.unwrap_or(self.pieces.len());
        // The last piece that starts at `pos` or before (pieces of 0 letters are passed over).
        let at = self.starts[..=known].partition_point(|&start| start <= pos) - 1;
        Ok((at < known).then_some(at))
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_map_ci.cpp:106-127 (`x_Move`)
    // NCBI reference (598d8ae6): c++/src/objmgr/seq_map_ci.cpp:483-538 (`x_Push`)
    /// Whether an iteration of the plus strand's positions `start..stop` throws: the find of
    /// `start`, or the move from each piece that it visits to the next, which resolves the next
    /// one's length (that of the first piece that starts at `stop` or after too). An iteration
    /// that finds no piece throws nothing.
    pub(crate) fn plus_range_throws(&self, start: u32, stop: u32) -> bool {
        let mut at = match self.find(start) {
            Ok(Some(at)) => at,
            Ok(None) => return false,
            Err(_) => return true,
        };
        loop {
            let next = at + 1;
            if Some(next) == self.first_invalid {
                return true;
            }
            if next >= self.pieces.len() || self.starts[next] >= stop {
                return false;
            }
            at = next;
        }
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_map.cpp:986-1117
    /// `CanResolveRange` of the range of `length` letters from `start` on the strand `minus`:
    /// within the map's length, no throw, and segments over all of it. The range is in the
    /// map's (plus) positions whatever the strand (`SSeqMapSelector::SetRange`): a minus-strand
    /// iteration finds its first segment at the range's last letter and goes back from there,
    /// over the pieces that the find resolved.
    fn can_resolve_range(&self, start: u32, length: u32, minus: bool) -> bool {
        let Ok(map_length) = self.length() else {
            return false;
        };
        let Some(stop) = start.checked_add(length) else {
            return false;
        };
        if start >= map_length || stop > map_length {
            return false;
        }
        if minus {
            return matches!(self.find(stop - 1), Ok(Some(_)));
        }
        if self.plus_range_throws(start, stop) {
            return false;
        }
        // The segments reach `stop`: before a segment of `kInvalidSeqPos` letters, which the
        // iteration did not reach, or at the end of the last.
        self.first_invalid.is_some() || self.end() >= stop
    }

    /// The number of pieces.
    fn count(&self) -> usize {
        self.pieces.len()
    }

    /// The position of the segment `idx` of an iterator (`0` and `count() + 1` the end
    /// markers).
    fn seg_start(&self, idx: usize) -> u32 {
        match idx {
            0 => 0,
            _ if idx <= self.count() => self.starts[idx - 1],
            _ => self.end(),
        }
    }

    /// The length of the segment `idx` (0 for the end markers).
    fn seg_length(&self, idx: usize) -> u32 {
        match idx {
            0 => 0,
            _ if idx <= self.count() => self.pieces[idx - 1].length(),
            _ => 0,
        }
    }

    /// Whether the segment `idx` is one of `kInvalidSeqPos` letters, whose resolution throws.
    fn unresolved(&self, idx: usize) -> bool {
        idx >= 1 && Some(idx - 1) == self.first_invalid
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_vector.cpp:1074-1104 (`sx_GetGapChar`)
    /// Appends the letters `from..end` of the plus strand in the coding of `CAlnVec`'s vector
    /// (ncbieaa for a protein, IUPACna otherwise; a gap's letter is `-` in ncbieaa and `N` in
    /// IUPACna). Only for letters that a read gave (`RowVector::read` without a throw): the
    /// segments and their Seq-data hold them.
    pub(crate) fn letters(&self, from: u32, end: u32, out: &mut Vec<u8>) {
        if from >= end {
            return;
        }
        let Ok(Some(mut at)) = self.find(from) else {
            debug_assert!(false, "letters that a read gave");
            return;
        };
        out.reserve((end - from) as usize);
        let mut pos = from;
        while pos < end && at < self.count() {
            let start = self.starts[at];
            let stop = start.saturating_add(self.pieces[at].length()).min(end);
            match self.pieces[at] {
                Piece::Gap(_) => out.extend(std::iter::repeat_n(
                    if self.protein { b'-' } else { b'N' },
                    (stop - pos) as usize,
                )),
                Piece::Data(data, _) => {
                    data_range_letters(data, (pos - start) as usize, (stop - start) as usize, out)
                }
            }
            pos = stop;
            at += 1;
        }
    }
}

/// The throw of a segment of `kInvalidSeqPos` letters (`queries::INVALID_LENGTH`).
fn invalid_length() -> RenderError {
    RenderError::Exception(NcbiException::new(
        INVALID_LENGTH,
        "Invalid sequence length",
    ))
}

/// A throw of the `CSeqMap` constructor of a `CSeq_inst` at `line`.
fn map_throw(line: u32, message: &str) -> RenderError {
    RenderError::Exception(NcbiException::new(
        ThrowSite {
            file: "src/objmgr/seq_map.cpp",
            line,
            kind: "CSeqMapException::eSeqDataError",
            function: Some("ncbi::objects::CSeqMap::CSeqMap"),
        },
        message,
    ))
}

/// Whether `CSeqVector` converts the coding `source` to `target` (`sx_GetConvertTable`): the
/// nucleotide codings but NCBI8na and NCBIpna to IUPACna or NCBI4na, the protein codings but
/// NCBI8aa and NCBIpaa to ncbieaa or NCBIstdaa.
pub(super) fn convertible(source: u32, target: u32) -> bool {
    matches!(
        (source, target),
        (IUPACNA | NCBI2NA | NCBI4NA, IUPACNA | NCBI4NA)
            | (NCBIEAA | IUPACAA | NCBISTDAA, NCBIEAA | NCBISTDAA)
    )
}

/// The letters `from..to` of the Seq-data `data`, which holds them, converted as
/// `CSeqVector` converts them to IUPACna or ncbieaa (`queries::data_letters`).
pub(super) fn data_range_letters(data: V<'_>, from: usize, to: usize, out: &mut Vec<u8>) {
    let (bytes, _) = data_bytes(data);
    match coding(data.choice().0) {
        IUPACNA | IUPACAA | NCBIEAA => out.extend_from_slice(&bytes[from..to]),
        NCBISTDAA => out.extend_from_slice(&protein_fasta_letters(&bytes[from..to])),
        NCBI4NA => out.extend((from..to).map(|i| {
            let byte = bytes[i / 2];
            NA4_TO_IUPACNA[usize::from(if i % 2 == 0 { byte >> 4 } else { byte & 15 })]
        })),
        NCBI2NA => out.extend((from..to).map(|i| {
            let shift = 6 - 2 * (i % 4);
            NA2_TO_IUPACNA[usize::from((bytes[i / 4] >> shift) & 3)]
        })),
        other => unreachable!("a coding that converts, not {other}"),
    }
}

/// The segment iterator of a vector (`CSeqMap_CI` over the whole map, on its strand): the
/// segment (`idx`: 0 and `count() + 1` the end markers), its position and length in the vector's
/// positions (on the minus strand, from the map's end), clipped to the vector's size.
#[derive(Clone, Copy, Debug)]
struct Seg {
    idx: usize,
    position: u32,
    length: u32,
}

impl Seg {
    fn end(&self) -> u32 {
        self.position.wrapping_add(self.length)
    }
}

// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:56 (`kCacheSize`)
/// A cache of the iterator (or its backup): where it starts and how many letters it holds.
#[derive(Clone, Copy, Debug, Default)]
struct Cache {
    pos: u32,
    size: u32,
}

impl Cache {
    fn end(&self) -> u32 {
        self.pos.wrapping_add(self.size)
    }
}

/// The state of `CSeqVector_CI` that decides its reads.
#[derive(Debug, Default)]
struct Iter {
    /// `m_Seg`, `None` before its first `x_InitSeg`.
    seg: Option<Seg>,
    cache: Cache,
    /// The position of `m_Cache` in its block `m_CacheData`.
    offset: u32,
    backup: Cache,
    /// Whether the cache's buffers are made (`x_InitializeCache`: `m_Cache` no longer null).
    made: bool,
}

// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector.cpp:152-163
// NCBI reference (598d8ae6): c++/include/objmgr/seq_vector.hpp:249-259
/// A `CSeqVector` of a Bioseq's map on one strand, with its iterator (`m_Iterator`, made by the
/// first read at its first position, `x_CreateIterator`, then moved, `SetPos`).
struct SeqVector<'m, 'a> {
    map: &'m SeqMap<'a>,
    /// `m_Size`: the map's length when the vector is made (`GetBioseqLength`).
    size: u32,
    minus: bool,
    /// The coding of `CAlnVec`'s vector (`x_GetSeqVector`): ncbieaa for a protein, IUPACna.
    target: u32,
    iter: Option<Iter>,
}

impl<'m, 'a> SeqVector<'m, 'a> {
    /// The vector of `map` on the strand `minus`: its size throws where the map's length does.
    fn new(map: &'m SeqMap<'a>, minus: bool) -> Result<Self, RenderError> {
        Ok(SeqVector {
            map,
            size: map.length()?,
            minus,
            target: if map.protein { NCBIEAA } else { IUPACNA },
            iter: None,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_vector.cpp:304-308
    // NCBI reference (598d8ae6): c++/include/objmgr/seq_vector_ci.hpp:517-525
    /// `GetSeqData(start, stop, buffer)`: the iterator made at `start` or moved there, moved
    /// there again, then the letters from `start` to `stop` (exclusive) up to the size.
    fn get_seq_data(&mut self, start: u32, stop: u32) -> Result<(), RenderError> {
        if self.iter.is_none() {
            self.iter = Some(Iter::default());
            self.x_set_pos(start)?;
        } else {
            self.set_pos(start)?;
        }
        self.set_pos(start)?;
        if start > stop {
            return Ok(());
        }
        self.read(stop - start)
    }

    fn it(&mut self) -> &mut Iter {
        self.iter
            .as_mut()
            .expect("the iterator made by the first read")
    }

    /// `SetPos`: within the cache, the offset; else `x_SetPos`.
    fn set_pos(&mut self, pos: u32) -> Result<(), RenderError> {
        let it = self.it();
        let offset = pos.wrapping_sub(it.cache.pos);
        if offset >= it.cache.size {
            self.x_set_pos(pos)
        } else {
            it.offset = offset;
            Ok(())
        }
    }

    /// `x_SwapCache`.
    fn swap_cache(&mut self) {
        let it = self.it();
        std::mem::swap(&mut it.cache, &mut it.backup);
        it.offset = 0;
    }

    /// `x_ResetCache`.
    fn reset_cache(&mut self) {
        let it = self.it();
        it.cache.size = 0;
        it.offset = 0;
    }

    /// `x_InitializeCache`: the buffers made the first time (the backup empty), else the cache
    /// emptied.
    fn initialize_cache(&mut self) {
        let it = self.it();
        if !it.made {
            it.made = true;
            it.backup.size = 0;
        }
        it.cache.size = 0;
        it.offset = 0;
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:631-681
    fn x_set_pos(&mut self, pos: u32) -> Result<(), RenderError> {
        let size = self.size;
        if pos >= size {
            if self.it().cache.size > 0 {
                self.swap_cache();
                self.reset_cache();
            }
            let it = self.it();
            it.cache.pos = size;
            it.offset = 0;
            return Ok(());
        }
        self.update_seg(pos)?;
        self.swap_cache();
        let cache = self.it().cache;
        let cache_offset = pos.wrapping_sub(cache.pos);
        if cache_offset < cache.size {
            self.it().offset = cache_offset;
            return Ok(());
        }
        self.initialize_cache();
        let old_pos = self.it().backup.pos;
        let seg_end = self.it().seg.expect("a segment found").end();
        if pos < old_pos && pos >= old_pos.wrapping_sub(CACHE) && seg_end >= old_pos {
            self.update_cache_down(old_pos - 1)?;
            let it = self.it();
            it.offset = pos - it.cache.pos;
        } else {
            self.update_cache_up(pos)?;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:484-510
    fn update_cache_up(&mut self, pos: u32) -> Result<(), RenderError> {
        let seg = self.it().seg.expect("a segment found");
        self.fill_cache(pos, CACHE.min(seg.end() - pos))?;
        self.it().offset = 0;
        Ok(())
    }

    fn update_cache_down(&mut self, pos: u32) -> Result<(), RenderError> {
        let seg = self.it().seg.expect("a segment found");
        let cache_offset = (CACHE - 1).min(pos - seg.position);
        self.fill_cache(pos - cache_offset, cache_offset + 1)?;
        self.it().offset = cache_offset;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:512-629
    // NCBI reference (598d8ae6): c++/include/objmgr/impl/seq_vector_cvt.hpp:47-131
    /// `x_FillCache` of `count` letters from `start` in the current segment: a gap's letters, or
    /// the Seq-data's converted (`Incompatible sequence codings`: a coding of the other
    /// molecule, NCBI8na, NCBIpna, NCBI8aa, NCBIpaa) from its position `dataPos` (from the
    /// segment's start on the plus strand, back from its end on the minus strand, which reads it
    /// reversed), whose end throws where it, its half or its quarter passes the data's bytes
    /// (`copy_8bit_any`, `copy_4bit_any`, `copy_2bit_any`). A packed coding's fill that ends
    /// inside the byte past the data reads memory that the input does not decide: refused.
    fn fill_cache(&mut self, start: u32, count: u32) -> Result<(), RenderError> {
        let (minus, target, map) = (self.minus, self.target, self.map);
        let it = self.it();
        it.cache.size = count;
        let seg = it.seg.expect("a segment found");
        if let Piece::Data(data, _) = map.pieces[seg.idx - 1] {
            let source = coding(data.choice().0);
            if source != target && !convertible(source, target) {
                return Err(RenderError::Exception(NcbiException::new(
                    INCOMPATIBLE_CODINGS,
                    format!("Incompatible sequence codings: {source} -> {target}"),
                )));
            }
            let offset = start - seg.position;
            let data_pos = if minus {
                seg.length - offset - count
            } else {
                offset
            };
            let (bytes, per_byte) = data_bytes(data);
            let end = u64::from(data_pos) + u64::from(count);
            if end / per_byte > bytes.len() as u64 {
                return Err(RenderError::Exception(NcbiException::new(
                    OUT_OF_RANGE,
                    format!("reference out of range of Seq-inst data: {end}"),
                )));
            }
            if end > bytes.len() as u64 * per_byte {
                return Err(RenderError::Unsupported(
                    "a packed Seq-data that ends inside its last byte before the Bioseq's \
                     length (NCBI reads past the data)"
                        .into(),
                ));
            }
        }
        self.it().cache.pos = start;
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:683-718
    /// `x_UpdateSeg`: the segment that holds `pos`, by the iterator's moves from the current one
    /// (over segments of 0 letters) or a new one; `cannot locate segment at` where none does.
    fn update_seg(&mut self, pos: u32) -> Result<(), RenderError> {
        let seg = match self.it().seg {
            None => self.init_seg(pos)?,
            Some(mut seg) if seg.position > pos => {
                loop {
                    self.prev(&mut seg)?;
                    if !(self.valid(&seg) && seg.length == 0) {
                        break;
                    }
                }
                if !self.valid(&seg) || seg.position > pos {
                    self.init_seg(pos)?
                } else {
                    seg
                }
            }
            Some(mut seg) if seg.end() <= pos => {
                loop {
                    self.next(&mut seg)?;
                    if !(self.valid(&seg) && seg.length == 0) {
                        break;
                    }
                }
                if !self.valid(&seg) || seg.end() <= pos {
                    self.init_seg(pos)?
                } else {
                    seg
                }
            }
            Some(seg) => seg,
        };
        self.it().seg = Some(seg);
        let valid = self.valid(&seg);
        if !valid && pos == self.size {
            return Ok(());
        }
        if !valid || pos < seg.position || pos >= seg.end() {
            return Err(RenderError::Exception(NcbiException::new(
                NO_SEGMENT,
                format!("CSeqVector_CI: cannot locate segment at {pos}"),
            )));
        }
        Ok(())
    }

    /// The clipped length of the segment `idx` (`x_CalcLength` over the range `0..size`).
    fn calc_length(&self, idx: usize) -> u32 {
        let start = self.map.seg_start(idx);
        let end = start.saturating_add(self.map.seg_length(idx));
        end.min(self.size).saturating_sub(start)
    }

    /// The vector's position of the segment `idx` (`x_GetTopOffset`).
    fn top_offset(&self, idx: usize) -> u32 {
        let start = self.map.seg_start(idx);
        if self.minus {
            let end = start.saturating_add(self.map.seg_length(idx));
            self.size.saturating_sub(end)
        } else {
            start.min(self.size)
        }
    }

    /// Whether the iterator is on a segment (`CSeqMap_CI::IsValid`: in range, not an end
    /// marker).
    fn valid(&self, seg: &Seg) -> bool {
        let start = self.map.seg_start(seg.idx);
        (1..=self.map.count()).contains(&seg.idx)
            && start < self.size
            && start.saturating_add(self.map.seg_length(seg.idx)) > 0
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:257-272 (`x_InitSeg`)
    // NCBI reference (598d8ae6): c++/src/objmgr/seq_map_ci.cpp:237-262 (`x_Select`)
    /// A new iterator at `pos` (below the size): the segment that holds it (on the minus strand
    /// the one that holds the map's position `size - 1 - pos`), or an end marker.
    fn init_seg(&self, pos: u32) -> Result<Seg, RenderError> {
        let find = if self.minus {
            self.size.wrapping_sub(1).wrapping_sub(pos)
        } else {
            pos
        };
        let idx = match self.map.find(find)? {
            Some(at) => at + 1,
            None if self.minus => 0,
            None => self.map.count() + 1,
        };
        let marker = !(1..=self.map.count()).contains(&idx);
        Ok(Seg {
            idx,
            position: self.top_offset(idx),
            // At an end marker `x_SettleNext` moves no further.
            length: if marker { 0 } else { self.calc_length(idx) },
        })
    }

    /// The move of `seg` toward the map's end (`forward`) or its start; whether it lands within
    /// the range. Toward the end it resolves the next segment's length.
    fn step(&self, seg: &mut Seg, forward: bool) -> Result<bool, RenderError> {
        if forward {
            if self.map.seg_start(seg.idx) > self.size || seg.idx > self.map.count() {
                return Ok(false);
            }
            seg.idx += 1;
            if self.map.unresolved(seg.idx) {
                return Err(invalid_length());
            }
            Ok(self.map.seg_start(seg.idx) < self.size)
        } else {
            if seg.idx == 0 {
                return Ok(false);
            }
            let old = self.map.seg_start(seg.idx);
            seg.idx -= 1;
            Ok(old > 0)
        }
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_map_ci.cpp:558-570 (`x_TopNext`)
    // NCBI reference (598d8ae6): c++/src/objmgr/seq_map_ci.cpp:685-692 (`x_SettleNext`)
    /// The step of `m_Seg` to the next segment of the vector, past the end markers.
    fn next(&self, seg: &mut Seg) -> Result<(), RenderError> {
        loop {
            seg.position = seg.position.wrapping_add(seg.length);
            if !self.step(seg, !self.minus)? {
                seg.length = 0;
                return Ok(());
            }
            seg.length = self.calc_length(seg.idx);
            if (1..=self.map.count()).contains(&seg.idx) {
                return Ok(());
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_map_ci.cpp:573-585 (`x_TopPrev`)
    /// `--m_Seg`: the segment before in the vector.
    fn prev(&self, seg: &mut Seg) -> Result<(), RenderError> {
        loop {
            if !self.step(seg, self.minus)? {
                seg.length = 0;
                return Ok(());
            }
            seg.length = self.calc_length(seg.idx);
            seg.position = seg.position.wrapping_sub(seg.length);
            if (1..=self.map.count()).contains(&seg.idx) {
                return Ok(());
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:720-755
    // NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:206-232 (`CanGetRange`)
    /// `GetSeqData(buffer, count)` from the iterator's position: at most the letters up to the
    /// size, all of which `CanGetRange` must find (`cannot get seq-data in range`), read chunk by
    /// chunk from the cache; a chunk that ends at the cache's end fills the next cache
    /// (`x_NextCacheSeg`), after the last chunk too.
    fn read(&mut self, count: u32) -> Result<(), RenderError> {
        let pos = {
            let it = self.it();
            it.cache.pos.wrapping_add(it.offset)
        };
        let count = count.min(self.size.wrapping_sub(pos));
        if count == 0 {
            return Ok(());
        }
        if !self.map.can_resolve_range(pos, count, self.minus) {
            return Err(RenderError::Exception(NcbiException::new(
                ThrowSite {
                    file: "src/objmgr/seq_vector_ci.cpp",
                    line: 731,
                    kind: "CSeqVectorException::eSeqDataError",
                    function: Some("ncbi::objects::CSeqVector_CI::GetSeqData"),
                },
                format!(
                    "CSeqVector_CI::GetSeqData: cannot get seq-data in range: {pos}-{}",
                    pos.wrapping_add(count)
                ),
            )));
        }
        let mut left = count;
        while left > 0 {
            let it = self.it();
            let available = it.cache.size - it.offset;
            let chunk = left.min(available);
            debug_assert!(chunk > 0, "a cache with letters");
            left -= chunk;
            if chunk == available {
                self.next_cache_seg()?;
            } else {
                it.offset += chunk;
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:757-799
    /// `x_NextCacheSeg`: at the size, an empty cache there; else the backup where it holds the
    /// position, or a fill from it.
    fn next_cache_seg(&mut self) -> Result<(), RenderError> {
        let pos = self.it().cache.end();
        if pos >= self.size {
            // A chunk of letters ends here: the cache starts before its end.
            debug_assert!(self.it().cache.pos < pos);
            self.swap_cache();
            self.reset_cache();
            self.it().cache.pos = pos;
            return Ok(());
        }
        self.swap_cache();
        self.update_seg(pos)?;
        let cache = self.it().cache;
        if pos < cache.end() && pos >= cache.pos {
            self.it().offset = pos - cache.pos;
        } else {
            self.reset_cache();
            self.update_cache_up(pos)?;
        }
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:109-140 (`x_GetSeqVector`)
// NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnvec.hpp:287-313 (`GetSeqString`)
/// The reads of one row of a `CAlnVec` (the tabular formats make one for each row) of a
/// request's Bioseq on the row's strand: its vector, made by the first read (`x_GetSeqVector`,
/// whose size throws where the map's length does), and its iterator, which each read moves on.
pub(crate) struct RowVector<'m, 'a> {
    map: &'m SeqMap<'a>,
    minus: bool,
    vector: Option<SeqVector<'m, 'a>>,
}

impl<'m, 'a> RowVector<'m, 'a> {
    /// The reads of a row of the Bioseq of `map` on the plus strand (`plus`) or the minus strand.
    pub(crate) fn new(map: &'m SeqMap<'a>, plus: bool) -> Self {
        RowVector {
            map,
            minus: !plus,
            vector: None,
        }
    }

    /// `GetSeqString` of the row from `seq_from` to `seq_to`: the letters `seq_from..=seq_to` of
    /// the plus strand (`seq_to` is `seq_from - 1` for a segment of 0 letters: a read of none,
    /// which still moves the iterator), read on the row's strand: on the minus strand from
    /// `size - seq_to - 1` to `size - seq_from` of the vector, in the unsigned `TSeqPos` (none
    /// where `seq_to` is past the size, `RowLetters`). The exception of the read, or none.
    pub(crate) fn read(&mut self, seq_from: u32, seq_to: u32) -> Result<(), RenderError> {
        let vector = match &mut self.vector {
            Some(vector) => vector,
            vector @ None => vector.insert(SeqVector::new(self.map, self.minus)?),
        };
        if vector.minus {
            let size = vector.size;
            vector.get_seq_data(
                size.wrapping_sub(seq_to).wrapping_sub(1),
                size.wrapping_sub(seq_from),
            )
        } else {
            vector.get_seq_data(seq_from, seq_to.wrapping_add(1))
        }
    }
}

#[cfg(test)]
#[path = "seq_vector_tests.rs"]
mod tests;
