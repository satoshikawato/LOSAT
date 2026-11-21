//! The scan of a user query with the prebuilt lookup table of a profile volume
//! (`BlastRPSScanSubject`): the hits of each word go to the bucket of 2,048 concatenated rows
//! that holds them, in the order of the scan, and a scan stops before the word whose hits would
//! pass 4,000,000; the next scan restarts at that word.
//!
//! The hits are made and read in passes of `RPS_GROUP_BUCKETS` buckets (series PR, item PS12m), so
//! that a worker holds the hits of 64 buckets, not of the whole volume: the scan itself makes the
//! hits of the first group and notes whether the scan has more (`scan_subject`); each later group
//! is made by a pass over the same words of the subject (`RpsBuckets::fill_group`), which keeps 8
//! bytes for a word whose cell has an overflow run (more than `RPS_HITS_PER_CELL` = 3 entries, the
//! table format's own limit of the entries in the cell) and hits past the first group, and nothing
//! for the other words; a pass is made only for the groups that can have hits. A volume of at most
//! 64 buckets is one group and has no further pass.

use crate::seqdb::rps::{
    entry_q_off, q_off_bucket, RPS_BUCKET_BITS, RPS_BUCKET_SIZE, RPS_GROUP_BUCKETS,
    RPS_HITS_PER_CELL, RPS_MAX_SCAN_HITS, RPS_WORDSIZE,
};

/// What a scan reads of a volume's table: its cells, overflow array and presence bits.
#[derive(Clone, Copy)]
pub(crate) struct RpsTable<'a> {
    pub(crate) cells: &'a [[i32; 4]],
    pub(crate) overflow: &'a [i32],
    pub(crate) pv: &'a [u64],
}

/// The bits of a concatenated position inside its bucket: `RPS_BUCKET_SIZE` is `2^BUCKET_BITS`, so
/// a position is its bucket's index times 2,048 plus a value below 2,048, for every position of
/// every volume.
const BUCKET_BITS: u32 = RPS_BUCKET_BITS;
const BUCKET_MASK: u32 = (1 << BUCKET_BITS) - 1;

/// The bits of a query offset above those an entry holds (`RpsBuckets`): an entry keeps the
/// offset's low `SEGMENT_SHIFT` bits, its block the rest (the offset's segment).
const SEGMENT_SHIFT: u32 = u32::BITS - BUCKET_BITS;

/// The entries of one block of a bucket (`RpsBuckets`): 4 KB.
pub(crate) const BLOCK_ENTRIES: usize = 1024;

/// One block of a bucket's entries.
type Block = Box<[u32; BLOCK_ENTRIES]>;

/// The blocks of a run's buckets that no bucket holds (series PR, item PS12s, M2): one pool for
/// every worker of the run, so that the blocks a bucket gives back as soon as its hits are read
/// (`RpsBuckets::next_block`) serve the scans of the other workers, and the run keeps the blocks
/// its workers' buckets hold at the same time, not each worker's largest scan. A worker takes one
/// block for 1,024 hits of a bucket and gives a bucket's blocks back together: a lock for a block,
/// none for a hit.
#[derive(Default)]
pub(crate) struct RpsBlockPool {
    free: std::sync::Mutex<Vec<Block>>,
}

impl RpsBlockPool {
    /// A block of the pool (a new one when the pool has none).
    fn take(&self) -> Block {
        let block = self
            .free
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .pop();
        block.unwrap_or_else(|| Box::new([0; BLOCK_ENTRIES]))
    }

    /// `blocks` back in the pool.
    fn give_back(&self, blocks: impl Iterator<Item = Block>) {
        self.free
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .extend(blocks);
    }

    /// The blocks in the pool.
    #[cfg(test)]
    fn len(&self) -> usize {
        self.free
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .len()
    }
}

/// A bucket's block closed before its current one: its entries (the first `len`) and the segment
/// of their query offsets.
struct ClosedBlock {
    block: Block,
    len: u32,
    segment: u32,
}

/// A bucket's current block: its entries (the first `fill`) and the segment of their query
/// offsets.
struct BucketHead {
    block: Block,
    fill: u32,
    segment: u32,
}

/// A word of the scan whose cell has an overflow run (more than `RPS_HITS_PER_CELL` entries: the
/// first is in the cell and the others in the overflow array; that number is NCBI's table format,
/// not a size chosen here) and has hits beyond the first group of buckets, with where
/// `RpsBuckets::fill_group` stands in its entries: 8 bytes, and such a cell is at least four of the
/// scan's 4,000,000 hits, so the words that a scan records take at most 8 MB. (The other words'
/// entries are all in the cell, which each pass reads again.)
#[derive(Clone, Copy)]
struct ScanWord {
    /// The query offset of the word's first letter (NCBI's `s_off`).
    s_off: u32,
    /// The entries of the word's cell already made into hits (a prefix of its entries).
    taken: u32,
}

/// The buckets of a volume (`lookup->bucket_array`, `rows / RPS_BUCKET_SIZE + 1` of them), of
/// which a worker holds those of one group at a time (`fill_group`), kept by the worker for its
/// run (series PR, item PS12s, D4): a bucket is a list of fixed-size blocks from the run's pool
/// (`RpsBlockPool`, M2), its hits in the order of the scan; the hits are read bucket by bucket,
/// each in the order of the scan (`next_block`), as NCBI reads its buckets, and each bucket's
/// blocks but its current one go back to the pool of free blocks as soon as they are read.
///
/// A hit takes 4 bytes (M1 of PS12s), where NCBI's `BlastOffsetPair` takes 8: the concatenated
/// position inside its bucket (`BUCKET_BITS` bits; the bucket gives the rest) and the low
/// `SEGMENT_SHIFT` bits of the query offset, whose high bits (its segment) the block holds once.
/// A segment holds 2^21 query offsets; when the scan reaches the next one, every bucket's current
/// block that holds hits is closed (`start_segment`), so a block's hits are in one segment and
/// every offset reads back whole, whatever the length of the query. A scan that crosses k segment
/// boundaries therefore holds at most min(hits, k x buckets) partly filled blocks (4 KB each) in a
/// group.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:197-208
pub(crate) struct RpsBuckets {
    /// The buckets of the volume (`lookup->num_buckets`).
    num_buckets: usize,
    /// The buckets of a group (`RPS_GROUP_BUCKETS`; a field so that tests can read other sizes).
    group_buckets: usize,
    /// The current group's buckets: `g0 .. g1`; `heads` and `chains` are indexed from `g0`.
    g0: usize,
    g1: usize,
    /// The words of the current scan: the subject's offsets `first..end` (what `scan_subject`
    /// scanned; every word of them with hits is in the scan's total).
    scanned: (i32, i32),
    /// The words of the current scan with an overflow run that have hits past the first group, in the
    /// order of the scan.
    words: Vec<ScanWord>,
    /// Whether the current scan has hits past the first group (`scan_subject`).
    beyond: bool,
    /// The groups that the entries of the current scan's inline cells (`RPS_HITS_PER_CELL` entries
    /// or fewer, in the cell itself) with hits past the first group lie in, one bit per group: a
    /// pass is needed for these, and for every group when `words` is not empty.
    reached: Vec<u64>,
    /// The first concatenated position past the current group (`u32::MAX` for the last group).
    g1_q: u32,
    heads: Vec<BucketHead>,
    /// Each bucket's closed blocks in the order they were filled.
    chains: Vec<Vec<ClosedBlock>>,
    /// The query offsets' segment of the scan's current word (`entry_bits`).
    segment: u32,
    /// The run's blocks that no bucket holds.
    pool: std::sync::Arc<RpsBlockPool>,
}

impl RpsBuckets {
    /// The buckets of a volume of `num_rows` concatenated rows, empty, their blocks from `pool`.
    pub(crate) fn new(num_rows: usize, pool: std::sync::Arc<RpsBlockPool>) -> Self {
        Self::with_group(num_rows, pool, RPS_GROUP_BUCKETS)
    }

    /// `new` with groups of `group_buckets` buckets (tests read other sizes than production).
    pub(crate) fn with_group(
        num_rows: usize,
        pool: std::sync::Arc<RpsBlockPool>,
        group_buckets: usize,
    ) -> Self {
        let mut buckets = Self {
            num_buckets: 0,
            group_buckets,
            g0: 0,
            g1: 0,
            scanned: (0, 0),
            words: Vec::new(),
            beyond: false,
            reached: Vec::new(),
            g1_q: u32::MAX,
            heads: Vec::new(),
            chains: Vec::new(),
            segment: 0,
            pool,
        };
        buckets.reset_for(num_rows);
        buckets
    }

    /// The buckets of a volume of `num_rows` concatenated rows, emptied; the blocks stay.
    pub(crate) fn reset_for(&mut self, num_rows: usize) {
        self.num_buckets = num_rows / RPS_BUCKET_SIZE + 1;
        let count = self.num_buckets.min(self.group_buckets);
        self.clear();
        if self.heads.len() > count {
            self.pool
                .give_back(self.heads.drain(count..).map(|head| head.block));
        }
        while self.heads.len() < count {
            let block = self.take();
            self.heads.push(BucketHead {
                block,
                fill: 0,
                segment: 0,
            });
        }
        self.chains.resize_with(count, Vec::new);
        self.beyond = false;
        self.words.clear();
    }

    /// The groups of the volume's buckets that hold hits of the current scan: the first one, and
    /// the others when the scan found hits past it (`scan_subject`).
    pub(crate) fn groups(&self) -> usize {
        if !self.beyond {
            1
        } else {
            self.num_buckets.div_ceil(self.group_buckets)
        }
    }

    /// Whether the scan has hits in `group` to make: the first group always (its hits are the
    /// scan's), a later one when a cell with an overflow run has hits past the first group (its
    /// entries are not known to lie in this group or not: `words`) or an inline cell has an entry
    /// in it (`reached`).
    pub(crate) fn needs_pass(&self, group: usize) -> bool {
        group == 0
            || !self.words.is_empty()
            || self
                .reached
                .get(group / 64)
                .is_some_and(|bits| bits & (1u64 << (group % 64)) != 0)
    }

    /// The groups of the entries of an inline cell that lie past the first group, noted in
    /// `reached`.
    fn reach_groups(&mut self, entries: &[i32]) {
        let g1_q = self.g1_q;
        for &entry in entries {
            let q_off = entry_q_off(entry);
            if q_off >= g1_q {
                let bucket = q_off_bucket(q_off);
                // A shift in production (64 buckets), a division for the odd sizes of the tests.
                let group = if self.group_buckets.is_power_of_two() {
                    bucket >> self.group_buckets.trailing_zeros()
                } else {
                    bucket / self.group_buckets
                };
                if self.reached.len() <= group / 64 {
                    self.reached.resize(group / 64 + 1, 0);
                }
                self.reached[group / 64] |= 1u64 << (group % 64);
            }
        }
    }

    /// Every bucket emptied, its closed blocks back in the pool.
    fn clear(&mut self) {
        for (head, chain) in self.heads.iter_mut().zip(&mut self.chains) {
            if !chain.is_empty() {
                self.pool
                    .give_back(chain.drain(..).map(|closed| closed.block));
            }
            head.fill = 0;
            head.segment = 0;
        }
        self.segment = 0;
    }

    /// The reading cursor at the first bucket of `group` with its buckets filled: group 0 holds
    /// the hits that `scan_subject` made, the other groups are made by a pass over the scan's
    /// words (`fill_group`; `subject` is the one `scan_subject` scanned).
    pub(crate) fn start_group(
        &mut self,
        table: &RpsTable<'_>,
        subject: &[u8],
        group: usize,
    ) -> RpsReadCursor {
        if group == 0 {
            return RpsReadCursor {
                bucket: 0,
                block: 0,
            };
        }
        self.fill_group(table, subject, group)
    }

    /// The buckets of `group` filled with the hits of the scan's words, the reading cursor at the
    /// group's first bucket. The pass takes the words of the subject again, in the scan's order
    /// (`scan_subject`'s loop over `scanned`: every word of it with hits was in the scan), and
    /// adds from each cell the entries that lie in the group. The entries of a cell are in groups
    /// in non-decreasing order (a table made by makeprofiledb lists them ascending, and
    /// `seqdb::rps` regroups the cells of a table that does not when it opens the volume), so a
    /// word's entries in the group follow the entries of the groups before it: for a cell with an
    /// overflow run the number of those is `taken` of its `ScanWord` (made by the scan for group 0 and kept
    /// by each pass), for a cell without one (up to `RPS_HITS_PER_CELL` entries, in the cell itself)
    /// they are the entries before the group's first position, found in the cell at once. So
    /// every bucket receives exactly the hits NCBI's `BlastRPSScanSubject` appends to it, in its
    /// order (the words in the order of the scan, a word's entries in the order of its cell), and
    /// the word finder reads the same hits in the same order.
    ///
    /// The state per word is the 8 bytes of a cell with an overflow run and hits past group 0 (such a
    /// cell is at least four of the scan's hits: at most 8 MB, where the dev scan's hits are 16 MB),
    /// and none for a cell without a run, so that a scan of millions of words with one hit each past the
    /// first group (a repeat of one rare word) costs no more memory than the dev scan. A pass
    /// costs a rescan of the subject's words (the presence bits and the cells, which are 512 KB)
    /// and no search in a cell's entries.
    pub(crate) fn fill_group(
        &mut self,
        table: &RpsTable<'_>,
        subject: &[u8],
        group: usize,
    ) -> RpsReadCursor {
        let cursor = self.open_group(group);
        let (first, end) = self.scanned;
        let g0_q = (self.g0 as u32) << BUCKET_BITS;
        let mut words = std::mem::take(&mut self.words);
        let mut next = 0;
        let mut index = word_index_start(subject, first);
        for s in first..end {
            index = next_word_index(index, subject, s);
            let word = index as usize;
            if table.pv[word / 64] & (1u64 << (word % 64)) == 0 {
                continue;
            }
            let cell = &table.cells[word];
            let (taken, entry) = if cell[0] > RPS_HITS_PER_CELL {
                // The words of cells with an overflow run and hits past group 0, in the order of the scan.
                if words.get(next).is_none_or(|w| w.s_off != s as u32) {
                    continue;
                }
                next += 1;
                let taken = words[next - 1].taken;
                let entry = match taken {
                    0 => Some(cell[1]),
                    _ if taken as i32 >= cell[0] => None,
                    _ => {
                        let start = cell[2] as usize / std::mem::size_of::<i32>();
                        Some(table.overflow[start + taken as usize - 1])
                    }
                };
                (taken, entry)
            } else {
                let entries = &cell[1..=cell[0] as usize];
                let taken = entries.partition_point(|&entry| entry_q_off(entry) < g0_q);
                (taken as u32, entries.get(taken).copied())
            };
            if entry.is_none_or(|entry| entry_q_off(entry) >= self.g1_q) {
                continue;
            }
            let s_bits = self.entry_bits(s as u32);
            let (taken, _) = self.add_cell_hits(table, cell, taken, s_bits, self.g0);
            if cell[0] > RPS_HITS_PER_CELL {
                words[next - 1].taken = taken;
            }
        }
        self.words = words;
        cursor
    }

    /// The hits of the entries of `cell` from the `taken`-th on that lie before `g1_q` (the end of
    /// the current group), into their buckets: the number of entries made so far and the position
    /// of the first one past the group (`u32::MAX` when there is none). The first entry is
    /// `cell[1]` whatever `cell[0]` is, the others the cell's own or the overflow run's
    /// (`BlastRPSScanSubject`).
    #[inline(always)]
    fn add_cell_hits(
        &mut self,
        table: &RpsTable<'_>,
        cell: &[i32; 4],
        taken: u32,
        s_bits: u32,
        g0: usize,
    ) -> (u32, u32) {
        let numhits = cell[0] as usize;
        if cell[0] <= RPS_HITS_PER_CELL {
            return self.add_list(&cell[1 + taken as usize..=numhits], taken, s_bits, g0);
        }
        let mut taken = taken;
        if taken == 0 {
            // The first hit is in the cell, the others in the overflow array.
            let (taken_now, past) = self.add_list(&cell[1..2], 0, s_bits, g0);
            if past != u32::MAX {
                return (taken_now, past);
            }
            taken = 1;
        }
        let start = cell[2] as usize / std::mem::size_of::<i32>();
        let run = &table.overflow[start + taken as usize - 1..start + numhits - 1];
        self.add_list(run, taken, s_bits, g0)
    }

    /// `s_AddToRPSBucket`: the hit of the concatenated position `q_off` whose query offset gives
    /// `s_bits` (`entry_bits`) at the end of the bucket of its row.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:422-436
    ///
    /// The hits of `list` (entries of a cell, `taken` of the cell's entries before it) into the
    /// buckets of the group (`g0` its first bucket; the lookup of the bucket is the one range test of
    /// a hit, where the dev loop has its bounds check): the number of the cell's entries made and
    /// the position of the first entry past the group (`u32::MAX` when there is none).
    #[inline(always)]
    fn add_list(&mut self, list: &[i32], taken: u32, s_bits: u32, g0: usize) -> (u32, u32) {
        for (made, &entry) in list.iter().enumerate() {
            let q_off = entry_q_off(entry);
            let bucket = q_off_bucket(q_off).wrapping_sub(g0);
            let Some(head) = self.heads.get_mut(bucket) else {
                return (taken + made as u32, self.past_the_group(q_off));
            };
            if head.fill as usize >= BLOCK_ENTRIES {
                self.close(bucket);
                let head = &mut self.heads[bucket];
                head.block[head.fill as usize] = s_bits | (q_off & BUCKET_MASK);
                head.fill += 1;
            } else {
                head.block[head.fill as usize] = s_bits | (q_off & BUCKET_MASK);
                head.fill += 1;
            }
        }
        (taken + list.len() as u32, u32::MAX)
    }

    /// The buckets of `group` emptied, the reading cursor at the first of them.
    fn open_group(&mut self, group: usize) -> RpsReadCursor {
        self.clear();
        self.g0 = group * self.group_buckets;
        self.g1 = (self.g0 + self.group_buckets).min(self.num_buckets);
        self.g1_q = if self.g1 == self.num_buckets {
            u32::MAX
        } else {
            (self.g1 as u32) << BUCKET_BITS
        };
        RpsReadCursor {
            bucket: self.g0,
            block: 0,
        }
    }

    /// A block of the run's pool.
    fn take(&mut self) -> Block {
        self.pool.take()
    }

    /// The bits that the entries of the hits of the word at the query offset `s_off` hold of it,
    /// after the buckets' blocks are closed where the offset starts a segment (`start_segment`).
    #[inline(always)]
    fn entry_bits(&mut self, s_off: u32) -> u32 {
        let segment = s_off >> SEGMENT_SHIFT;
        if segment != self.segment {
            self.start_segment(segment);
        }
        s_off << BUCKET_BITS
    }

    /// The scan reaches the query offsets of `segment`: every current block that holds hits is
    /// closed, and every bucket's next hits go to blocks of that segment.
    #[inline(never)]
    fn start_segment(&mut self, segment: u32) {
        for bucket in 0..self.heads.len() {
            if self.heads[bucket].fill > 0 {
                self.close(bucket);
            }
            self.heads[bucket].segment = segment;
        }
        self.segment = segment;
    }

    /// One hit into the buckets (`add_list`'s body), for the tests.
    #[cfg(test)]
    fn add(&mut self, q_off: u32, s_bits: u32) {
        let entry = q_off as i32 + (RPS_WORDSIZE as i32 - 1);
        let g0 = self.g0;
        let (_, past) = self.add_list(&[entry], 0, s_bits, g0);
        assert_eq!(past, u32::MAX, "a hit past the group");
    }

    /// `add_list` found `q_off` outside the buckets of the group (the one range test of a hit, where
    /// the dev loop has its bounds check): `q_off` itself when it lies past the group, the end of the
    /// entries of its cell that belong to the group (the cells are in groups in order). A position
    /// before the group means the table's cells are not: only a table that `seqdb::rps` did not
    /// regroup when the volume was opened reaches it, so the volume's files were not opened as the
    /// search needs.
    #[cold]
    #[inline(never)]
    fn past_the_group(&self, q_off: u32) -> u32 {
        if q_off >= self.g1_q {
            return q_off;
        }
        panic!(
            "RPS scan: a hit at position {q_off} lies before the buckets {}..{} of its group; \
             the profile table's cells must be in groups of buckets in order \
             (`seqdb::rps::regroup_cells`, applied when the volume is opened)",
            self.g0,
            self.g0 + self.heads.len()
        )
    }

    /// The current block of `bucket` closed at the end of the bucket's list, a block of the pool
    /// in its place.
    #[inline(never)]
    fn close(&mut self, bucket: usize) {
        let fresh = self.take();
        let head = &mut self.heads[bucket];
        let block = std::mem::replace(&mut head.block, fresh);
        self.chains[bucket].push(ClosedBlock {
            block,
            len: head.fill,
            segment: head.segment,
        });
        head.fill = 0;
    }

    /// The next block of hits of the current group (`cursor`: where the reading stands, from
    /// `fill_group`), bucket by bucket, each bucket's blocks in the order they were filled: every
    /// hit in the order NCBI reads its buckets. A bucket is emptied when the reading passes it, its
    /// closed blocks back in the pool.
    pub(crate) fn next_block(&mut self, cursor: &mut RpsReadCursor) -> Option<RpsHits<'_>> {
        while cursor.bucket < self.g1 {
            let bucket = cursor.bucket - self.g0;
            let q_base = (cursor.bucket as u32) << BUCKET_BITS;
            let closed = self.chains[bucket].len();
            if cursor.block < closed {
                cursor.block += 1;
                let block = &self.chains[bucket][cursor.block - 1];
                return Some(RpsHits {
                    entries: &block.block[..block.len as usize],
                    q_base,
                    s_base: block.segment << SEGMENT_SHIFT,
                });
            }
            if cursor.block == closed && self.heads[bucket].fill > 0 {
                cursor.block += 1;
                let head = &self.heads[bucket];
                return Some(RpsHits {
                    entries: &head.block[..head.fill as usize],
                    q_base,
                    s_base: head.segment << SEGMENT_SHIFT,
                });
            }
            // The bucket is read: emptied, its closed blocks back in the pool.
            let chain = &mut self.chains[bucket];
            if !chain.is_empty() {
                self.pool
                    .give_back(chain.drain(..).map(|closed| closed.block));
            }
            self.heads[bucket].fill = 0;
            cursor.bucket += 1;
            cursor.block = 0;
        }
        None
    }

    /// Every hit of the current scan as (`q_off`, `s_off`), group by group in the order of
    /// `next_block`; the buckets emptied.
    #[cfg(test)]
    pub(crate) fn drain(&mut self, table: &RpsTable<'_>, subject: &[u8]) -> Vec<(u32, u32)> {
        let mut hits = Vec::new();
        for group in 0..self.groups() {
            if !self.needs_pass(group) {
                continue;
            }
            let cursor = self.start_group(table, subject, group);
            hits.extend(self.read_group(cursor));
        }
        hits
    }

    /// The hits of the open group from `cursor` on, as (`q_off`, `s_off`); the buckets emptied.
    #[cfg(test)]
    fn read_group(&mut self, mut cursor: RpsReadCursor) -> Vec<(u32, u32)> {
        let mut hits = Vec::new();
        while let Some(block) = self.next_block(&mut cursor) {
            hits.extend(block.pairs());
        }
        hits
    }

    /// The hits of the bucket `bucket` (of the current group) as (`q_off`, `s_off`), in order
    /// (not emptied).
    #[cfg(test)]
    fn bucket_hits(&self, bucket: usize) -> Vec<(u32, u32)> {
        let q_base = (bucket as u32) << BUCKET_BITS;
        let head = &self.heads[bucket - self.g0];
        let blocks = self.chains[bucket - self.g0]
            .iter()
            .map(|closed| (&closed.block[..closed.len as usize], closed.segment))
            .chain(std::iter::once((
                &head.block[..head.fill as usize],
                head.segment,
            )));
        blocks
            .flat_map(|(entries, segment)| {
                RpsHits {
                    entries,
                    q_base,
                    s_base: segment << SEGMENT_SHIFT,
                }
                .pairs()
            })
            .collect()
    }

    /// The bytes of the blocks the buckets hold or the pool keeps.
    #[cfg(test)]
    pub(crate) fn kept_bytes(&self) -> usize {
        self.blocks() * BLOCK_ENTRIES * std::mem::size_of::<u32>()
    }

    /// The blocks the buckets hold or the pool keeps.
    #[cfg(test)]
    fn blocks(&self) -> usize {
        self.heads.len() + self.chains.iter().map(Vec::len).sum::<usize>() + self.pool.len()
    }
}

/// The hits of one block of a bucket (`RpsBuckets::next_block`): its entries, the bucket's first
/// concatenated position and the first query offset of the block's segment. A hit's concatenated
/// position is `q_base` plus its entry's `q_place`, its query offset `s_base` plus its `s_place`
/// (the places are below `2^BUCKET_BITS` and `2^SEGMENT_SHIFT`: the sums are the bits joined).
pub(crate) struct RpsHits<'a> {
    pub(crate) entries: &'a [u32],
    pub(crate) q_base: u32,
    pub(crate) s_base: u32,
}

impl<'a> RpsHits<'a> {
    /// The hits in the order of the scan: the concatenated position of the word's first letter
    /// and the query offset (NCBI's `q_off` and `s_off`).
    #[cfg(test)]
    pub(crate) fn pairs(&self) -> impl Iterator<Item = (u32, u32)> + 'a {
        let (q_base, s_base) = (self.q_base, self.s_base);
        self.entries
            .iter()
            .map(move |&entry| (q_base + q_place(entry), s_base + s_place(entry)))
    }
}

/// The place of an entry's hit in its bucket (`RpsHits`).
#[inline(always)]
pub(crate) fn q_place(entry: u32) -> u32 {
    entry & BUCKET_MASK
}

/// The place of an entry's query offset in its block's segment (`RpsHits`).
#[inline(always)]
pub(crate) fn s_place(entry: u32) -> u32 {
    entry >> BUCKET_BITS
}

/// Where the reading of the buckets stands (`RpsBuckets::next_block`): a bucket and the number of
/// its blocks read.
#[derive(Default)]
pub(crate) struct RpsReadCursor {
    bucket: usize,
    block: usize,
}

/// `ComputeTableIndex(wordsize - 1, charsize, s_first)`: the index of the first `wordsize - 1`
/// letters of the word at `s_first`, which `next_word_index` completes.
#[inline(always)]
fn word_index_start(subject: &[u8], s_first: i32) -> i32 {
    let charsize = crate::seqdb::rps::RPS_CHARSIZE as i32;
    let mut index = 0i32;
    for i in 0..RPS_WORDSIZE as i32 - 1 {
        index = (index << charsize) | i32::from(subject[(s_first + i) as usize]);
    }
    index
}

/// `ComputeTableIndexIncremental`: the index of the word at `s`, from the index before it (the
/// word's last letter in).
#[inline(always)]
fn next_word_index(index: i32, subject: &[u8], s: i32) -> i32 {
    let charsize = crate::seqdb::rps::RPS_CHARSIZE as i32;
    let mask = (crate::seqdb::rps::RPS_BACKBONE_SIZE - 1) as i32;
    ((index << charsize) | i32::from(subject[(s + RPS_WORDSIZE as i32 - 1) as usize])) & mask
}

/// `BlastRPSScanSubject`: the hits of the words of `subject` (its letters, `subject->sequence`)
/// from `*offset` into the emptied buckets of the first group (the code below is NCBI's, with the
/// hits past the group left to `RpsBuckets::fill_group`, which makes them from the same cells in
/// the same order, after the scan, group by group); `*offset` becomes the word the next scan starts
/// at. Returns the number of hits.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:472-550
/// The table's hits name rows inside the concatenation (checked when the volume is opened,
/// `seqdb::rps`), so every bucket index is in range.
pub(crate) fn scan_subject(
    table: &RpsTable<'_>,
    subject: &[u8],
    offset: &mut i32,
    buckets: &mut RpsBuckets,
    max_hits: i32,
) -> i32 {
    // The hits of group 0 go to its buckets as the scan finds them; a word with hits past the
    // group makes the scan's later passes (`fill_group`) necessary.
    buckets.open_group(0);
    buckets.beyond = false;
    buckets.words.clear();
    buckets.reached.iter_mut().for_each(|bits| *bits = 0);
    let s_first = *offset;
    let s_last = subject.len() as i32 - RPS_WORDSIZE as i32;
    let mut totalhits = 0i32;
    // `ComputeTableIndex(wordsize - 1, charsize, s_first)`.
    let mut index = word_index_start(subject, s_first);
    let mut s = s_first;
    while s <= s_last {
        index = next_word_index(index, subject, s);
        let word = index as usize;
        if table.pv[word / 64] & (1u64 << (word % 64)) != 0 {
            let cell = &table.cells[word];
            let numhits = cell[0];
            if numhits <= max_hits - totalhits {
                let s_bits = buckets.entry_bits(s as u32);
                let (taken, next_q) = buckets.add_cell_hits(table, cell, 0, s_bits, 0);
                if next_q != u32::MAX {
                    buckets.beyond = true;
                    if numhits > RPS_HITS_PER_CELL {
                        buckets.words.push(ScanWord {
                            s_off: s as u32,
                            taken,
                        });
                    } else {
                        buckets.reach_groups(&cell[1 + taken as usize..=numhits as usize]);
                    }
                }
                totalhits += numhits;
            } else {
                break;
            }
        }
        s += 1;
    }
    *offset = s;
    buckets.scanned = (s_first, s);
    totalhits
}

/// The hits of one scan as NCBI collects them (`RPS_MAX_SCAN_HITS`).
pub(crate) const MAX_SCAN_HITS: i32 = RPS_MAX_SCAN_HITS;

#[cfg(test)]
#[path = "scan_tests.rs"]
mod tests;
