//! NCBI's index factory (`CDbIndex_Factory`, `dbindex_factory.cpp`): one index volume made from
//! the sequences of a stream (`CDbIndex::MakeIndex`), with its subject map
//! (`CSubjectMap_Factory`: the packed sequences, their chunks and the logical sequence ids) and
//! its offset lists (`COffsetData_Factory`: for each Nmer, the encoded positions where it ends).
//!
//! The offset lists keep NCBI's order and contents (each list in the order its entries were
//! added, truncated as NCBI truncates it); LOSAT stores them in blocks of 15 words chained by
//! index, with three words per Nmer for the heads (NCBI: a `COffsetList` of 48 bytes per Nmer
//! and units of 21 words), so the bytes written are NCBI's with less memory.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::rc::Rc;

use super::format::Header;
use super::stream::{SeqData, SequenceStream, StreamError};
use super::{code_bits, min_offset, Options, SeqNum, Word, CR};

/// `CDbIndex_Factory::MEGABYTE`.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1537
const MEGABYTE: u64 = 1024 * 1024;

/// `CSubjectMap_Factory_Base::CMaskHelper`: an iterator over the masked intervals of a sequence
/// (the packed locations of its `mask_locs_`, each a list of `(from, to)`), which `In` moves
/// forward and `Adjust` back, as NCBI's does (the chunks of a long sequence overlap).
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:169-241
#[derive(Default)]
struct MaskHelper {
    /// `c_locs_`.
    locs: Vec<Vec<(u32, u32)>>,
    /// `vit_`: an index of `locs`; `locs.len()` is its end.
    vit: usize,
    /// `it_`: an index of `locs[vit]`.
    it: usize,
    /// `start_` and `stop_` (one past the end) of the interval at `it`.
    start: u32,
    stop: u32,
}

impl MaskHelper {
    fn new(locs: Vec<Vec<(u32, u32)>>) -> Self {
        let mut helper = MaskHelper {
            locs,
            ..MaskHelper::default()
        };
        helper.init();
        helper
    }

    /// `Good()`.
    fn good(&self) -> bool {
        self.vit != self.locs.len()
    }

    /// `start_` and `stop_` of the interval at `it` (`stop_` is one more than `GetTo()`, a
    /// `TSeqPos`).
    fn set(&mut self) {
        let (from, to) = self.locs[self.vit][self.it];
        self.start = from;
        self.stop = to.wrapping_add(1);
    }

    /// `Init()`: the first interval.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:538-553
    fn init(&mut self) {
        self.vit = 0;
        while self.vit != self.locs.len() {
            self.it = 0;
            if !self.locs[self.vit].is_empty() {
                self.set();
                break;
            }
            self.vit += 1;
        }
    }

    /// `Advance()`: the next interval.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:556-568
    fn advance(&mut self) {
        while self.good() {
            self.it += 1;
            if self.it != self.locs[self.vit].len() {
                self.set();
                return;
            }
            self.vit += 1;
            if self.good() {
                self.it = 0;
            }
        }
    }

    /// `Adjust` of `pos`: back to the first interval to the left of `pos`, or to the start.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:571-578
    fn adjust(&mut self, pos: u32) {
        loop {
            let notdone = self.retreat();
            if !(notdone && pos < self.stop) {
                break;
            }
        }
    }

    /// `Retreat()`: the previous interval; false when there is none.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:581-632
    fn retreat(&mut self) -> bool {
        if self.locs.is_empty() {
            return false;
        }
        if !self.good() {
            self.vit -= 1;
            while self.vit != 0 && self.locs[self.vit].is_empty() {
                self.vit -= 1;
            }
            if !self.locs[self.vit].is_empty() {
                self.it = self.locs[self.vit].len() - 1;
                self.set();
                return true;
            }
            self.vit = self.locs.len();
            return false;
        }
        if self.it != 0 {
            self.it -= 1;
            self.set();
            return true;
        }
        if self.vit == 0 {
            self.init();
            return false;
        }
        self.vit -= 1;
        while self.vit != 0 && self.locs[self.vit].is_empty() {
            self.vit -= 1;
        }
        if !self.locs[self.vit].is_empty() {
            self.it = self.locs[self.vit].len() - 1;
            self.set();
            return true;
        }
        self.init();
        false
    }

    /// `In` of `pos`: whether `pos` is in the interval reached by moving forward.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:635-640
    fn contains(&mut self, pos: u32) -> bool {
        while self.good() && pos >= self.stop {
            self.advance();
        }
        if !self.good() {
            return false;
        }
        pos >= self.start
    }
}

/// `SSeqInfo`: a chunk of a sequence (a logical sequence).
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:368-389
struct SeqInfo {
    /// `seq_start_`: the chunk's first byte in the packed sequences.
    seq_start: Word,
    /// `segs_`: the valid intervals `(start, stop)` of the chunk (all but the masked and the
    /// ambiguous bases), relative to its start. Emptied once the offset lists have them.
    segs: Vec<(u32, u32)>,
}

/// `SLIdMapElement`: a local sequence id, the chunks combined under it.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:452-458
#[derive(Clone, Copy)]
struct LidEntry {
    start: SeqNum,
    end: SeqNum,
    seq_start: u32,
    seq_end: u32,
}

/// `CSubjectMap_Factory` (with its bases `CSubjectMap_Factory_TBase` and `_Base`).
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:118-535
struct SubjectMap {
    chunk_size: u64,
    chunk_overlap: u64,
    /// `last_chunk_`: the number of the last valid chunk.
    last_chunk: SeqNum,
    /// `c_chunk_`: the chunk of the current sequence.
    c_chunk: SeqNum,
    /// `c_seq_`: the current sequence (its letters); kept by a `next` that returned none.
    c_seq: Option<Rc<SeqData>>,
    /// `seq_store_`: the packed sequences (`stride` zero bytes first).
    seq_store: Vec<u8>,
    /// `subjects_`: for each sequence, the number of its first chunk (1-based; 0 for none).
    subjects: Vec<SeqNum>,
    mask_helper: MaskHelper,
    stride: u64,
    min_offset: u64,
    /// `chunks_`.
    chunks: Vec<SeqInfo>,
    /// `lengths_`: the length of each sequence.
    lengths: Vec<Word>,
    /// `lid_map_`.
    lid_map: Vec<LidEntry>,
    /// `cur_lid_len_`.
    cur_lid_len: u32,
    /// `offset_bits_`.
    offset_bits: u8,
}

impl SubjectMap {
    /// `CSubjectMap_Factory` of `options`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:273-285
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:804-811
    fn new(options: &Options) -> Self {
        let min_offset = min_offset(options.stride);
        let max_len = (1 + options.chunk_size / options.stride) + min_offset;
        let mut offset_bits = 16u8;
        while (max_len >> offset_bits) != 0 {
            offset_bits += 1;
        }
        SubjectMap {
            chunk_size: options.chunk_size,
            chunk_overlap: options.chunk_overlap,
            last_chunk: 0,
            c_chunk: 0,
            c_seq: None,
            seq_store: vec![0; options.stride as usize],
            subjects: Vec::new(),
            mask_helper: MaskHelper::default(),
            stride: options.stride,
            min_offset,
            chunks: Vec::new(),
            lengths: Vec::new(),
            lid_map: Vec::new(),
            cur_lid_len: 0,
            offset_bits,
        }
    }

    /// `c_seq_.size()`, a `TSeqPos`.
    fn seq_len(&self) -> u32 {
        self.c_seq.as_ref().map_or(0, |seq| seq.bases.len() as u32)
    }

    /// `NewSequenceInit` with `sd` and `start_chunk` (with the lengths table of
    /// `CSubjectMap_Factory`): a new sequence, or none at the end of the input, which keeps the
    /// last one's letters (its length goes into the table again). The first word of the title,
    /// `unknown` for none.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:665-686
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:480-485
    fn new_sequence_init(&mut self, sd: Option<Rc<SeqData>>, start_chunk: SeqNum) -> Vec<u8> {
        let mut result = b"unknown".to_vec();
        self.subjects.push(0);
        self.c_chunk = start_chunk;
        if let Some(sd) = sd {
            result = sd.id.clone().unwrap_or_default();
            self.mask_helper = MaskHelper::new(sd.masks.clone());
            self.c_seq = Some(sd);
        }
        self.lengths.push(self.seq_len());
        result
    }

    /// `CSubjectMap_Factory_TBase::AddSequenceChunk` of `seq_off`: the valid intervals of the
    /// next chunk of the current sequence; false when it has no more chunks.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:714-778
    fn add_chunk_base(&mut self, seq_off: usize) -> bool {
        let chunk_start = ((self.chunk_size - self.chunk_overlap) * u64::from(self.c_chunk)) as u32;
        self.c_chunk += 1;
        let seq_len = self.seq_len();
        if chunk_start >= seq_len {
            self.c_chunk -= 1;
            return false;
        }
        let chunk_end = ((u64::from(chunk_start) + self.chunk_size) as u32).min(seq_len);
        let chunk_len = chunk_end - chunk_start;
        let mut segs: Vec<(u32, u32)> = Vec::new();
        if chunk_len > 0 {
            let invalid = self.c_seq.as_ref().expect("a sequence").bases.invalid();
            let mut run = invalid.partition_point(|&(_, to)| to <= chunk_start);
            let mut inside = false;
            self.mask_helper.adjust(chunk_start);
            for pos in chunk_start..chunk_end {
                while run < invalid.len() && invalid[run].1 <= pos {
                    run += 1;
                }
                let valid = run == invalid.len() || pos < invalid[run].0;
                // The masks are asked of the letters A, C, G and T only, as NCBI asks them.
                let in1 = !valid || self.mask_helper.contains(pos);
                if in1 && !inside {
                    if segs.is_empty() {
                        segs.push((0, 0));
                    }
                    segs.last_mut().expect("a segment").1 = pos - chunk_start;
                    inside = true;
                } else if !in1 && inside {
                    segs.push((pos - chunk_start, 0));
                    inside = false;
                }
            }
            if !inside {
                if segs.is_empty() {
                    segs.push((0, 0));
                }
                segs.last_mut().expect("a segment").1 = chunk_end - chunk_start;
            }
        }
        self.chunks.push(SeqInfo {
            seq_start: seq_off as Word,
            segs,
        });
        if let Some(first) = self.subjects.last_mut() {
            if *first == 0 {
                *first = self.chunks.len() as SeqNum;
            }
        }
        self.last_chunk = self.chunks.len() as SeqNum;
        true
    }

    /// `CSubjectMap_Factory::AddSequenceChunk` of `overflow`: the next chunk, under the current
    /// local id or a new one; the whole sequence is packed with its first chunk. `(added,
    /// overflow)`: overflow when a new local id was needed and none is left.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:814-876
    fn add_sequence_chunk(&mut self) -> Result<(bool, bool), StreamError> {
        let starting = self.c_chunk == 0;
        let step = self.chunk_size - self.chunk_overlap;
        let chunk_start = (step * u64::from(self.c_chunk)) as u32;
        let seq_off = if starting {
            self.seq_store.len()
        } else {
            self.chunks.last().expect("a chunk").seq_start as usize + (step / CR) as usize
        };
        if !self.add_chunk_base(seq_off) {
            return Ok((false, false));
        }
        let seq_len = self.seq_len();
        let length_limit: u32 = 1 << (self.offset_bits - 1);
        let chunk_end = ((u64::from(chunk_start) + self.chunk_size) as u32).min(seq_len);
        let chunk_len = chunk_end - chunk_start;
        if self.lid_map.is_empty() || self.cur_lid_len.wrapping_add(chunk_len) > length_limit {
            let lid_bits = 32 - u32::from(self.offset_bits);
            let lid_limit = 1u64 << lid_bits;
            if self.lid_map.len() as u64 >= lid_limit {
                return Ok((true, true));
            }
            self.lid_map.push(LidEntry {
                start: self.chunks.len() as SeqNum - 1,
                end: 0,
                seq_start: seq_off as u32,
                seq_end: 0,
            });
            self.cur_lid_len = 0;
        }
        let chunks = self.chunks.len() as SeqNum;
        self.cur_lid_len = self.cur_lid_len.wrapping_add(chunk_len);
        let lid = self.lid_map.last_mut().expect("a local id");
        lid.end = chunks;
        lid.seq_end = lid.seq_start.wrapping_add(self.cur_lid_len);
        if starting && seq_len > 0 {
            let packed = self.c_seq.as_ref().expect("a sequence").bases.packed();
            // The packed sequences are indexed by NCBI's 32-bit words (`seq_start_`).
            if self.seq_store.len() + packed.len() > u32::MAX as usize {
                return Err(StreamError::Refused(super::unsupported(
                    "an index volume of more than 4 GiB of packed sequences (-volsize)",
                )));
            }
            self.seq_store.extend_from_slice(packed);
        }
        Ok((true, false))
    }

    /// The local id of the chunk whose data start at byte `soff` (the last one that starts at
    /// or before it) and the position of that byte in it.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:879-888
    fn lid_of(&self, soff: u32) -> (usize, u32) {
        let at = self
            .lid_map
            .iter()
            .rposition(|lid| lid.seq_start <= soff)
            .expect("a local id at or before the chunk");
        (
            at,
            (soff - self.lid_map[at].seq_start).wrapping_mul(CR as u32),
        )
    }

    /// `MakeOffset` (with `seq` and `off`) of the chunk whose data start at byte `soff`, by its
    /// local id `lid` and base `base` (`lid_of`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:891-903
    #[inline]
    fn make_offset(&self, (lid, base): (usize, u32), off: u32) -> Word {
        let mut off = u64::from(off.wrapping_add(base)) / self.stride;
        off = u64::from((off as u32).wrapping_add(self.min_offset as u32));
        (((lid as u64) << self.offset_bits) + off) as Word
    }

    /// `MakeOffset` with `seqnum` and `off`: of chunk `seqnum` (0-based).
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:906-912
    fn make_chunk_offset(&self, seqnum: SeqNum, off: u32) -> Word {
        let soff = self.chunks[seqnum as usize].seq_start;
        self.make_offset(self.lid_of(soff), off)
    }

    /// `GetSeqInfo` of `snum`: chunk `snum` (1-based) while it is valid.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:411-418
    fn has_seq_info(&self, snum: SeqNum) -> bool {
        snum <= self.last_chunk
    }

    /// `RollBack()`: the current sequence is dropped.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:781-788
    fn roll_back(&mut self) {
        if let Some(first) = self.subjects.last_mut() {
            self.last_chunk = first.wrapping_sub(1);
            self.c_chunk = 0;
            *first = 0;
        }
    }

    /// `Commit()`: the chunks past the last valid one and their data are dropped.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:791-801
    fn commit(&mut self) {
        if (self.last_chunk as usize) < self.chunks.len() {
            let newsize = self.chunks[self.last_chunk as usize].seq_start as usize;
            self.seq_store.truncate(newsize);
            self.chunks.truncate(self.last_chunk as usize);
        }
    }

    /// `total()`: the size of the packed sequences.
    fn total(&self) -> Word {
        self.seq_store.len() as Word
    }

    /// `Save` of `CSubjectMap_Factory` (to `os`), then of its base.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:915-938
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:689-711
    fn save(&self, os: &mut impl Write) -> std::io::Result<()> {
        write_word(os, (4 * self.lengths.len()) as Word)?;
        write_word(os, Word::from(self.offset_bits))?;
        for &length in &self.lengths {
            write_word(os, length)?;
        }
        write_word(os, (16 * self.lid_map.len()) as Word)?;
        for lid in &self.lid_map {
            for word in [lid.start, lid.end, lid.seq_start, lid.seq_end] {
                write_word(os, word)?;
            }
        }
        write_word(
            os,
            (4 * self.subjects.len() + 4 * self.chunks.len()) as Word,
        )?;
        for &first in &self.subjects {
            write_word(os, first)?;
        }
        for chunk in &self.chunks {
            write_word(os, chunk.seq_start)?;
        }
        write_word(os, self.seq_store.len() as Word)?;
        write_word(os, self.seq_store.len() as Word)?;
        os.write_all(&self.seq_store)
    }
}

/// `WriteWord` with `os` and `word`: the word's bytes (little-endian, as NCBI's on its targets).
#[inline]
fn write_word(os: &mut impl Write, word: Word) -> std::io::Result<()> {
    os.write_all(&word.to_le_bytes())
}

/// The words of a block of the offset lists (NCBI's units hold 21).
const UNIT: usize = 15;

/// A block of an offset list: its words and the number of the next block (0 for none).
#[derive(Clone, Copy)]
struct Block {
    data: [Word; UNIT],
    next: u32,
}

/// The blocks in a page (NCBI's pool takes blocks of a million units).
const PAGE: usize = 1 << 16;

/// The blocks of the offset lists (NCBI's `COffsetList::CDataPool`), numbered from 1, with the
/// list of the blocks that truncation freed.
#[derive(Default)]
struct Pool {
    pages: Vec<Box<[Block]>>,
    /// The blocks handed out from the pages.
    count: u32,
    /// The first free block, 0 for none.
    free: u32,
}

impl Pool {
    fn block(&self, id: u32) -> &Block {
        let at = (id - 1) as usize;
        &self.pages[at / PAGE][at % PAGE]
    }

    fn block_mut(&mut self, id: u32) -> &mut Block {
        let at = (id - 1) as usize;
        &mut self.pages[at / PAGE][at % PAGE]
    }

    /// `alloc()`: a block whose next is none.
    fn alloc(&mut self) -> u32 {
        let id = if self.free != 0 {
            let id = self.free;
            self.free = self.block(id).next;
            id
        } else {
            if self.count as usize == self.pages.len() * PAGE {
                self.pages.push(
                    vec![
                        Block {
                            data: [0; UNIT],
                            next: 0
                        };
                        PAGE
                    ]
                    .into_boxed_slice(),
                );
            }
            self.count += 1;
            self.count
        };
        self.block_mut(id).next = 0;
        id
    }

    /// `free` of a block: the chain from block `id` goes to the free blocks.
    fn free_chain(&mut self, id: u32) {
        if id == 0 {
            return;
        }
        let mut last = id;
        while self.block(last).next != 0 {
            last = self.block(last).next;
        }
        let free = self.free;
        self.block_mut(last).next = free;
        self.free = id;
    }
}

/// `COffsetData_Factory`: the offset lists, one for each Nmer (`hash_table_`), each a chain of
/// blocks (`heads`, `tails`, `sizes`: NCBI's `COffsetList` and its `CData`).
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:945-1219
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1292-1396
struct OffsetData {
    heads: Vec<u32>,
    tails: Vec<u32>,
    sizes: Vec<u32>,
    pool: Pool,
    /// `total_`: the words of the lists (a `TWord`).
    total: Word,
    hkey_width: u64,
    /// `last_seq_`: the last chunk the lists have.
    last_seq: SeqNum,
    stride: u64,
    code_bits: u64,
    min_offset: u64,
    /// `mult_` of the lists: `ws_hint` less `hkey_width`, plus 1, divided by `stride`.
    mult: u64,
}

impl OffsetData {
    /// `COffsetData_Factory` (with `subject_map`, `options` and `pool`) and each list's
    /// `SetIndexParams`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1303-1320
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:953-957
    fn new(options: &Options) -> Self {
        let lists = 1usize << (2 * options.hkey_width);
        OffsetData {
            heads: vec![0; lists],
            tails: vec![0; lists],
            sizes: vec![0; lists],
            pool: Pool::default(),
            total: 0,
            hkey_width: options.hkey_width,
            last_seq: 0,
            stride: options.stride,
            code_bits: code_bits(options.stride),
            min_offset: min_offset(options.stride),
            mult: (options.ws_hint - options.hkey_width + 1) / options.stride,
        }
    }

    /// `AddData` (with `item` and `total_`) of the list of `nmer`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1255-1259
    #[inline]
    fn add_data(&mut self, nmer: usize, item: Word) {
        let size = self.sizes[nmer] as usize;
        let slot = size % UNIT;
        if size == 0 {
            let id = self.pool.alloc();
            self.heads[nmer] = id;
            self.tails[nmer] = id;
        } else if slot == 0 {
            let id = self.pool.alloc();
            let tail = self.tails[nmer];
            self.pool.block_mut(tail).next = id;
            self.tails[nmer] = id;
        }
        let tail = self.tails[nmer];
        self.pool.block_mut(tail).data[slot] = item;
        self.sizes[nmer] += 1;
        self.total = self.total.wrapping_add(1);
    }

    /// The words of the list of `nmer`, in order.
    fn list(&self, nmer: usize) -> impl Iterator<Item = Word> + '_ {
        let size = self.sizes[nmer] as usize;
        let mut block = self.heads[nmer];
        (0..size).map(move |i| {
            if i > 0 && i % UNIT == 0 {
                block = self.pool.block(block).next;
            }
            self.pool.block(block).data[i % UNIT]
        })
    }

    /// `TruncateList` (with `offset` and `total_`) of the list of `nmer`: the entries from the
    /// first that is at least `offset` are dropped, with the code before it. Every list is read
    /// from its start, as NCBI reads it: a list is not always in increasing order (with `-stride 1`
    /// a local id of many short sequences has positions whose offsets reach into the local id's
    /// bits, above the offsets of a later local id).
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1262-1286
    fn truncate_list(&mut self, nmer: usize, offset: Word) {
        let size = self.sizes[nmer] as usize;
        if size == 0 {
            return;
        }
        let min_offset = self.min_offset;
        let mut flag = false;
        let mut cut = None;
        for (i, word) in self.list(nmer).enumerate() {
            if u64::from(word) < min_offset {
                flag = true;
                continue;
            }
            if word >= offset {
                cut = Some(if flag { i - 1 } else { i });
                break;
            }
            flag = false;
        }
        let Some(newsize) = cut else {
            return;
        };
        self.total = self.total.wrapping_sub((size - newsize) as Word);
        self.resize(nmer, newsize);
    }

    /// `CData::resize` to the smaller size `newsize`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1178-1201
    fn resize(&mut self, nmer: usize, newsize: usize) {
        if newsize == 0 {
            self.pool.free_chain(self.heads[nmer]);
            self.heads[nmer] = 0;
            self.tails[nmer] = 0;
            self.sizes[nmer] = 0;
            return;
        }
        let mut block = self.heads[nmer];
        for _ in 1..newsize.div_ceil(UNIT) {
            block = self.pool.block(block).next;
        }
        let rest = self.pool.block(block).next;
        self.pool.free_chain(rest);
        self.pool.block_mut(block).next = 0;
        self.tails[nmer] = block;
        self.sizes[nmer] = newsize as u32;
    }

    /// `EncodeAndAddOffset` with `nmer`, `start`, `stop`, `curr` and `offset`: the offset, after
    /// the code of its distances to the ends of its valid interval when one of them is at most the
    /// stride.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1450-1467
    #[inline]
    fn encode_and_add_offset(
        &mut self,
        nmer: Word,
        start: u32,
        stop: u32,
        curr: u32,
        offset: Word,
    ) {
        let mut start_diff = curr
            .wrapping_add(2)
            .wrapping_sub(self.hkey_width as u32)
            .wrapping_sub(start);
        let mut end_diff = stop.wrapping_sub(curr);
        if u64::from(start_diff) <= self.stride || u64::from(end_diff) <= self.stride {
            if u64::from(start_diff) > self.stride {
                start_diff = 0;
            }
            if u64::from(end_diff) > self.stride {
                end_diff = 0;
            }
            let code = (start_diff << self.code_bits).wrapping_add(end_diff);
            self.add_data(nmer as usize, code);
        }
        self.add_data(nmer as usize, offset);
    }

    /// `AddSeqSeg` with `seq`, `seqlen`, `start` and `stop`: the Nmers ending in the valid interval
    /// `start .. stop` of the chunk whose data start at byte `soff`, at the positions that are
    /// multiples of the stride in its local id (`CheckOffset`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1470-1490
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:879-888
    fn add_seq_seg(&mut self, map: &SubjectMap, soff: u32, start: u32, stop: u32) {
        let nmer_mask: Word = ((1u64 << (2 * self.hkey_width)) - 1) as Word;
        let lid = map.lid_of(soff);
        let seq = &map.seq_store[soff as usize..];
        let mut nmer: Word = 0;
        let mut count = 0u64;
        for curr in start..stop {
            let unit = seq[(curr / CR as u32) as usize];
            let letter = (unit >> (6 - 2 * (curr % CR as u32))) & 0x3;
            nmer = ((nmer << 2) & nmer_mask) + Word::from(letter);
            if count >= self.hkey_width - 1
                && u64::from(curr.wrapping_add(lid.1)) % self.stride == 0
            {
                let offset = map.make_offset(lid, curr);
                self.encode_and_add_offset(nmer, start, stop, curr, offset);
            }
            count += 1;
        }
    }

    /// `Truncate()`: the lists lose what they have of the chunks past the last valid one.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1504-1513
    fn truncate(&mut self, map: &SubjectMap) {
        self.last_seq = map.last_chunk;
        let offset = map.make_chunk_offset(self.last_seq, 0);
        for nmer in 0..self.sizes.len() {
            self.truncate_list(nmer, offset);
        }
    }

    /// `Update()`: the lists take the valid chunks of the subject map that they do not have.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1516-1528
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1493-1501
    fn update(&mut self, map: &mut SubjectMap) {
        if map.last_chunk < self.last_seq {
            self.truncate(map);
        }
        while map.has_seq_info(self.last_seq + 1) {
            let at = self.last_seq as usize;
            let segs = std::mem::take(&mut map.chunks[at].segs);
            let soff = map.chunks[at].seq_start;
            for &(start, stop) in &segs {
                self.add_seq_seg(map, soff, start, stop);
            }
            self.last_seq += 1;
        }
    }

    /// `COffsetList::Save` (to `os`): the list's codes, each with its offset, and the offsets that
    /// are multiples of `mult_`, then for each `m` from one less than `mult_` down to 1 the other
    /// offsets whose largest divisor up to `mult_` is `m`, then 0 (for a list that is not empty).
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1222-1252
    fn save_list(
        &self,
        nmer: usize,
        os: &mut impl Write,
        later: &mut Vec<(u64, Word)>,
    ) -> std::io::Result<()> {
        later.clear();
        let mut words = self.list(nmer);
        while let Some(word) = words.next() {
            if u64::from(word) < self.min_offset {
                write_word(os, word)?;
                write_word(os, words.next().expect("the offset after a code"))?;
            } else if u64::from(word) % self.mult == 0 {
                write_word(os, word)?;
            } else {
                // The largest divisor below `mult_` (`word` itself when `word` is below it).
                let class = (1..=(self.mult - 1).min(u64::from(word)))
                    .rev()
                    .find(|&m| u64::from(word) % m == 0)
                    .expect("1 divides every offset");
                later.push((class, word));
            }
        }
        later.sort_by(|a, b| b.0.cmp(&a.0));
        for &(_, word) in later.iter() {
            write_word(os, word)?;
        }
        if self.sizes[nmer] != 0 {
            write_word(os, 0)?;
        }
        Ok(())
    }

    /// `Save` (to `os`): the total of the words (with the 0 that ends each list and the 0 before
    /// the first), the hash table (each Nmer's first word in the lists, 0 for an empty list), then
    /// the lists; with `stat_file_name`, the size of each list that is not empty.
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1399-1447
    fn save(&mut self, os: &mut impl Write, options: &Options) -> std::io::Result<()> {
        self.total = self.total.wrapping_add(1);
        for &size in &self.sizes {
            if size > 0 {
                self.total = self.total.wrapping_add(1);
            }
        }
        let mut stats = options
            .stat_file_name
            .as_deref()
            .and_then(|name| File::create(crate::utils::os_bytes::path_of_bytes(name)).ok())
            .map(BufWriter::new);
        write_word(os, self.total)?;
        let mut tot: Word = 0;
        for (nmer, &size) in self.sizes.iter().enumerate() {
            if size != 0 {
                tot = tot.wrapping_add(1);
                write_word(os, tot)?;
            } else {
                write_word(os, 0)?;
            }
            tot = tot.wrapping_add(size);
            if size > 0 {
                if let Some(stats) = &mut stats {
                    let _ = writeln!(stats, "{nmer:>10x} {size}");
                }
            }
        }
        if let Some(mut stats) = stats {
            let _ = stats.flush();
        }
        write_word(os, self.total)?;
        write_word(os, 0)?;
        let mut later = Vec::new();
        for nmer in 0..self.sizes.len() {
            if self.sizes[nmer] != 0 {
                self.save_list(nmer, os, &mut later)?;
            }
        }
        Ok(())
    }
}

/// Why a volume was not made, after what NCBI's would have posted.
pub(crate) enum Stop {
    Stream(StreamError),
    /// A volume for which NCBI's makembindex gives no defined index: its `Save` does not end
    /// (a zero `mult_`), or its offsets come from a shift past an `int`'s width (a stride of
    /// 32768 or more, C++ undefined behaviour).
    NoDefinedIndex,
}

impl From<StreamError> for Stop {
    fn from(error: StreamError) -> Self {
        Stop::Stream(error)
    }
}

/// `CDbIndex::MakeIndex` with `input`, `oname`, `start`, `stop` and `options` (`do_create_1_2`):
/// the index volume `oname` of the sequences from `start` of `input`; `stop` becomes the number
/// after its last sequence (`start` when the input had none left, and then nothing is written). A
/// sequence that the volume cannot take whole goes back to the input. `post` gets the messages
/// of the input's sequences; `cerr` gets the lines NCBI writes to `std::cerr`.
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1657-1760
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/dbindex_factory.cpp:1797-1803
pub(crate) fn make_index(
    input: &mut dyn SequenceStream,
    oname: &[u8],
    start: SeqNum,
    stop: &mut SeqNum,
    options: &Options,
    post: &mut dyn FnMut(&[u8]),
    cerr: &mut dyn FnMut(&[u8]),
) -> Result<(), Stop> {
    let mut subject_map = SubjectMap::new(options);
    let mut offset_data = OffsetData::new(options);
    let mut i = start;
    if i >= *stop {
        *stop = start;
        return Ok(());
    }
    let mut idmap: Vec<Vec<u8>> = Vec::new();
    let mut start_chunk: SeqNum = 0;
    let mut stop_chunk: SeqNum = 0;
    while i < *stop {
        let sd = input.next();
        for line in input.take_warnings() {
            post(&line);
        }
        let sd = sd?;
        if let Some(sd) = &sd {
            for text in &sd.texts {
                post(text.as_bytes());
            }
        }
        let real = sd.is_some();
        let idstr = subject_map.new_sequence_init(sd, start_chunk);
        idmap.push(idstr);
        if !real {
            if i == start {
                *stop = start;
                return Ok(());
            }
            *stop = i;
            stop_chunk = 0;
            break;
        }
        // NCBI gives no defined volume with a sequence for a zero `mult_` (the loop over the
        // divisor classes of its `Save` starts below 1 and never ends) or a stride of 32768 or
        // more (`GetMinOffset` shifts an `int` by 32 bits or more, which C++ leaves undefined):
        // `makembindex` refuses them here, at the volume's first sequence.
        if offset_data.mult == 0 || options.stride >= super::MAX_STRIDE {
            return Err(Stop::NoDefinedIndex);
        }
        loop {
            let (added, overflow) = subject_map.add_sequence_chunk()?;
            if !added {
                break;
            }
            if !overflow {
                offset_data.update(&mut subject_map);
            } else {
                cerr(b"WARNING: logical sequence id overflow. Starting new volume.\n");
            }
            let total = u64::from(subject_map.total()) + 4 * u64::from(offset_data.total);
            if total > MEGABYTE * options.max_index_size || overflow {
                input.putback();
                subject_map.roll_back();
                offset_data.update(&mut subject_map);
                subject_map.commit();
                *stop = start
                    .wrapping_add(subject_map.subjects.len() as SeqNum)
                    .wrapping_sub(1);
                stop_chunk = subject_map.c_chunk;
                break;
            }
        }
        subject_map.commit();
        start_chunk = 0;
        i = i.wrapping_add(1);
    }
    // NCBI opens the volume as a binary output stream without exceptions: a file that does not
    // open, or a write that fails, is not reported.
    if let Ok(file) = File::create(crate::utils::os_bytes::path_of_bytes(oname)) {
        let mut os = BufWriter::with_capacity(1 << 20, file);
        let header = Header::of(
            options.legacy,
            options.hkey_width as Word,
            options.stride as Word,
            options.ws_hint as Word,
            [start, start_chunk, *stop, stop_chunk],
        );
        let _ = header
            .write(&mut os)
            .and_then(|()| offset_data.save(&mut os, options))
            .and_then(|()| subject_map.save(&mut os))
            .and_then(|()| os.flush());
    } else {
        let _ = offset_data.save(&mut std::io::sink(), options);
    }
    if options.idmap {
        let map = [oname, b".map"].concat();
        if let Ok(file) = File::create(crate::utils::os_bytes::path_of_bytes(&map)) {
            let mut maps = BufWriter::new(file);
            for id in &idmap {
                let _ = maps.write_all(id);
                let _ = maps.write_all(b"\n");
            }
            let _ = maps.flush();
        }
    }
    Ok(())
}
