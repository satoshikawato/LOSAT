//! The scan of a user query with the prebuilt lookup table of a profile volume
//! (`BlastRPSScanSubject`): the hits of each word go to the bucket of 2,048 concatenated rows
//! that holds them, in the order of the scan, and a scan stops before the word whose hits would
//! pass 4,000,000; the next scan restarts at that word.
//!
//! The hits are made and read in passes of `RPS_GROUP_BUCKETS` buckets (series PR, item PS12m), so
//! that a worker holds the hits of 64 buckets, not of the whole volume: the scan itself makes the
//! hits of the first group and notes whether the scan has more (`scan_subject`); each later group
//! is made by a pass over the same words of the subject (`RpsBuckets::fill_group`), which keeps
//! no state per word. A volume of at most 64 buckets is one group and has no further pass.

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
/// ```c
///     lookup->num_buckets = num_pssm_rows / RPS_BUCKET_SIZE + 1;
///     lookup->bucket_array = (RPSBucket *) malloc(lookup->num_buckets *
///                                                 sizeof(RPSBucket));
///     for (i = 0; i < lookup->num_buckets; i++) {
///         RPSBucket *bucket = lookup->bucket_array + i;
///         bucket->num_filled = 0;
///         bucket->num_alloc = 1000;
/// ```
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
    /// Whether the current scan has hits past the first group (`scan_subject`).
    beyond: bool,
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
            beyond: false,
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
    /// from each cell the entries that lie in the group, which are the entries from the first
    /// one at or past the group's first position up to the first past the group. The entries of a
    /// cell are in groups in non-decreasing order (a table made by makeprofiledb lists them
    /// ascending, and `seqdb::rps` regroups the cells of a table that does not when it opens the
    /// volume), so that first entry is found by a binary search and the group's entries follow
    /// it. So every bucket receives exactly the hits NCBI's `BlastRPSScanSubject` appends to it,
    /// in its order (the words in the order of the scan, a word's entries in the order of its
    /// cell), and the word finder reads the same hits in the same order.
    ///
    /// The pass keeps nothing per word (the worker's state is the group's hits, at most the
    /// scan's 4,000,000 of 4 bytes, as the dev scan's): a scan of millions of words that have
    /// one hit each past the first group costs no more memory than the dev scan, and a word
    /// without a hit in the group costs a binary search in its cell.
    pub(crate) fn fill_group(
        &mut self,
        table: &RpsTable<'_>,
        subject: &[u8],
        group: usize,
    ) -> RpsReadCursor {
        let cursor = self.open_group(group);
        let (first, end) = self.scanned;
        let g0_q = (self.g0 as u32) << BUCKET_BITS;
        let mut index = word_index_start(subject, first);
        for s in first..end {
            index = next_word_index(index, subject, s);
            let word = index as usize;
            if table.pv[word / 64] & (1u64 << (word % 64)) == 0 {
                continue;
            }
            let cell = &table.cells[word];
            let (from, entry) = cell_group_start(cell, table.overflow, g0_q);
            if entry.is_none_or(|entry| entry_q_off(entry) >= self.g1_q) {
                continue;
            }
            let s_bits = self.entry_bits(s as u32);
            self.add_cell_hits(table, cell, from, s_bits, self.g0);
        }
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
        let buckets = self.heads.len();
        if cell[0] <= RPS_HITS_PER_CELL {
            return self.add_list(
                &cell[1 + taken as usize..=numhits],
                taken,
                s_bits,
                g0,
                buckets,
            );
        }
        let mut taken = taken;
        if taken == 0 {
            // The first hit is in the cell, the others in the overflow array.
            let (taken_now, past) = self.add_list(&cell[1..2], 0, s_bits, g0, buckets);
            if past != u32::MAX {
                return (taken_now, past);
            }
            taken = 1;
        }
        let start = cell[2] as usize / std::mem::size_of::<i32>();
        let run = &table.overflow[start + taken as usize - 1..start + numhits - 1];
        self.add_list(run, taken, s_bits, g0, buckets)
    }

    /// `s_AddToRPSBucket`: the hit of the concatenated position `q_off` whose query offset gives
    /// `s_bits` (`entry_bits`) at the end of the bucket of its row.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:422-436
    /// ```c
    /// static void s_AddToRPSBucket(RPSBucket * b, Uint4 q_off, Uint4 s_off)
    /// {
    ///     ...
    ///     offset_pairs[i].qs_offsets.q_off = q_off;
    ///     offset_pairs[i].qs_offsets.s_off = s_off;
    ///     b->num_filled++;
    /// }
    /// ```
    ///
    /// The hits of `list` (entries of a cell, `taken` of the cell's entries before it) into the
    /// buckets of the group (`g0` its first bucket, `buckets` their number: the one range test of
    /// a hit, where the dev loop has its bounds check): the number of the cell's entries made and
    /// the position of the first entry past the group (`u32::MAX` when there is none).
    #[inline(always)]
    fn add_list(
        &mut self,
        list: &[i32],
        taken: u32,
        s_bits: u32,
        g0: usize,
        buckets: usize,
    ) -> (u32, u32) {
        for (made, &entry) in list.iter().enumerate() {
            let q_off = entry_q_off(entry);
            let bucket = q_off_bucket(q_off).wrapping_sub(g0);
            if bucket >= buckets {
                return (taken + made as u32, self.past_the_group(q_off));
            }
            if self.heads[bucket].fill as usize >= BLOCK_ENTRIES {
                self.close(bucket);
            }
            let head = &mut self.heads[bucket];
            head.block[head.fill as usize] = s_bits | (q_off & BUCKET_MASK);
            head.fill += 1;
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
        let (g0, buckets) = (self.g0, self.heads.len());
        let (_, past) = self.add_list(&[entry], 0, s_bits, g0, buckets);
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

/// Where the entries of a cell with `cell[0] > 0` that lie at or past the position `g0_q` start:
/// the number of the entries before them (the first one is `cell[1]` whatever `cell[0]` is, the
/// others the cell's own or the overflow run's) and the first of them. The entries are in groups
/// in non-decreasing order, so the entries before `g0_q` are a prefix: a binary search finds it.
#[inline(always)]
fn cell_group_start(cell: &[i32; 4], overflow: &[i32], g0_q: u32) -> (u32, Option<i32>) {
    let numhits = cell[0] as usize;
    let before = |entry: &i32| entry_q_off(*entry) < g0_q;
    if cell[0] <= RPS_HITS_PER_CELL {
        let entries = &cell[1..=numhits];
        let from = entries.partition_point(before);
        (from as u32, entries.get(from).copied())
    } else if !before(&cell[1]) {
        (0, Some(cell[1]))
    } else {
        let start = cell[2] as usize / std::mem::size_of::<i32>();
        let run = &overflow[start..start + numhits - 1];
        let from = run.partition_point(before);
        (from as u32 + 1, run.get(from).copied())
    }
}

/// `BlastRPSScanSubject`: the hits of the words of `subject` (its letters, `subject->sequence`)
/// from `*offset` into the emptied buckets of the first group (the code below is NCBI's, with the
/// hits past the group left to `RpsBuckets::fill_group`, which makes them from the same cells in
/// the same order, after the scan, group by group); `*offset` becomes the word the next scan starts
/// at. Returns the number of hits.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:472-547
/// ```c
///     const Int4 max_hits = 4000000;
///     ...
///     for (index = 0; index < lookup->num_buckets; index++)
///         bucket_array[index].num_filled = 0;
///
///     s_first = abs_start + *offset;
///     s_last = abs_start + sequence->length - lookup->wordsize;
///     pv = lookup->pv;
///     ...
///     table_correction = lookup->wordsize - 1;
///
///     /* prime the index */
///     index = ComputeTableIndex(lookup->wordsize - 1,
///                           lookup->charsize, s_first);
///
///     for (s = s_first; s <= s_last; s++) {
///         /* compute the index value */
///         index = ComputeTableIndexIncremental(lookup->wordsize,
///                                              lookup->charsize,
///                                              lookup->mask, s, index);
///
///         /* if there are hits... */
///         if (PV_TEST(pv, index, PV_ARRAY_BTS)) {
///             cell = &lookup->rps_backbone[index];
///             numhits = cell->num_used;
///
///             ASSERT(numhits != 0);
///
///             if (numhits <= (max_hits - totalhits)) {
///                 Int4 *src;
///                 Int4 i;
///                 Uint4 q_off;
///                 Uint4 s_off = (Uint4)(s - abs_start);
///                 if (numhits <= RPS_HITS_PER_CELL) {
///                     for (i = 0; i < numhits; i++) {
///                         q_off = cell->entries[i] - table_correction;
///                         s_AddToRPSBucket(bucket_array +
///                                          q_off / RPS_BUCKET_SIZE, q_off,
///                                          s_off);
///                     }
///                 } else {
///                     /* hits (past the first) live in overflow array */
///                     src =
///                         lookup->overflow + (cell->entries[1] / sizeof(Int4));
///                     q_off = cell->entries[0] - table_correction;
///                     s_AddToRPSBucket(bucket_array + q_off / RPS_BUCKET_SIZE,
///                                      q_off, s_off);
///                     for (i = 0; i < (numhits - 1); i++) {
///                         q_off = src[i] - table_correction;
///                         s_AddToRPSBucket(bucket_array +
///                                          q_off / RPS_BUCKET_SIZE, q_off,
///                                          s_off);
///                     }
///                 }
///
///                 totalhits += numhits;
///             } else
///                 /* not enough space in the destination array; return early */
///             {
///                 break;
///             }
///         }
///     }
///
///     /* if we get here, we fell off the end of the sequence */
///     *offset = (Int4)(s - abs_start);
///
///     return totalhits;
/// ```
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
                let (_, next_q) = buckets.add_cell_hits(table, cell, 0, s_bits, 0);
                buckets.beyond |= next_q != u32::MAX;
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
mod tests {
    use super::*;
    use crate::seqdb::rps::RPS_BACKBONE_SIZE;

    fn word(letters: [u8; 3]) -> usize {
        (usize::from(letters[0]) << 10) | (usize::from(letters[1]) << 5) | usize::from(letters[2])
    }

    struct Table {
        cells: Vec<[i32; 4]>,
        overflow: Vec<i32>,
        pv: Vec<u64>,
    }

    impl Table {
        fn new() -> Self {
            Self {
                cells: vec![[0; 4]; RPS_BACKBONE_SIZE + 1],
                overflow: Vec::new(),
                pv: vec![0; RPS_BACKBONE_SIZE / 64],
            }
        }

        fn set(&mut self, w: usize, hits: &[i32]) {
            let mut cell = [hits.len() as i32, 0, 0, 0];
            if hits.len() <= 3 {
                cell[1..=hits.len()].copy_from_slice(hits);
            } else {
                cell[1] = hits[0];
                cell[2] = (self.overflow.len() * 4) as i32;
                self.overflow.extend_from_slice(&hits[1..]);
            }
            self.cells[w] = cell;
            self.pv[w / 64] |= 1 << (w % 64);
        }

        fn view(&self) -> RpsTable<'_> {
            RpsTable {
                cells: &self.cells,
                overflow: &self.overflow,
                pv: &self.pv,
            }
        }
    }

    // Hits go to the bucket of their row in the order of the scan; a cell's overflow run keeps
    // its order.
    #[test]
    fn hits_fill_the_buckets_in_scan_order() {
        let mut table = Table::new();
        table.set(word([1, 2, 3]), &[4100, 10]);
        table.set(word([2, 3, 1]), &[20, 2060, 30, 2050, 5000]);
        let subject = [1u8, 2, 3, 1, 2, 3];
        let mut buckets = RpsBuckets::new(6000, Default::default());
        let mut offset = 0;
        let hits = scan_subject(
            &table.view(),
            &subject,
            &mut offset,
            &mut buckets,
            MAX_SCAN_HITS,
        );
        assert_eq!(hits, 2 + 5 + 2);
        assert_eq!(offset, 4);
        let pairs = |b: usize| -> Vec<(u32, u32)> { buckets.bucket_hits(b) };
        // Word 1 2 3 at s 0 and 3 (rows of the last letter less 2); word 2 3 1 at s 1.
        assert_eq!(pairs(0), vec![(8, 0), (18, 1), (28, 1), (8, 3)]);
        assert_eq!(pairs(1), vec![(2058, 1), (2048, 1)]);
        assert_eq!(pairs(2), vec![(4098, 0), (4998, 1), (4098, 3)]);
    }

    // A scan stops before the word whose hits do not fit, and the next one starts there.
    #[test]
    fn a_full_scan_restarts_at_the_word_that_did_not_fit() {
        let mut table = Table::new();
        table.set(word([1, 1, 1]), &[2, 3, 4]);
        let subject = [1u8; 6];
        let mut buckets = RpsBuckets::new(10, Default::default());
        let mut offset = 0;
        // Room for seven hits: two words (six hits), then the third does not fit.
        let first = scan_subject(&table.view(), &subject, &mut offset, &mut buckets, 7);
        assert_eq!((first, offset), (6, 2));
        assert_eq!(buckets.bucket_hits(0).len(), 6);
        let second = scan_subject(&table.view(), &subject, &mut offset, &mut buckets, 7);
        assert_eq!((second, offset), (6, 4));
        assert_eq!(buckets.bucket_hits(0)[0], (0, 2));
        let last = scan_subject(&table.view(), &subject, &mut offset, &mut buckets, 7);
        assert_eq!((last, offset), (0, 4));
    }

    /// A deterministic sequence of pseudo-random numbers (no NCBI counterpart: test data).
    fn lcg(state: &mut u64) -> u32 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (*state >> 33) as u32
    }

    // D4 and M1 of PS12s: the blocks of a bucket hold its hits in the order of the scan, and
    // `next_block` reads them bucket by bucket, as one `Vec` of NCBI's pairs per bucket holds and
    // reads them; a drained bucket's closed blocks go back to the pool, which keeps its blocks
    // from scan to scan and volume to volume (fewer or more buckets): what the largest scan's
    // buckets held.
    // Query offsets cross segments of 2^21 (also inside a bucket's block) and read back whole.
    #[test]
    fn blocks_read_the_hits_in_the_order_of_vectors_per_bucket() {
        let mut state = 12_345u64;
        // One group of every bucket (the hits are added to it directly).
        let mut buckets = RpsBuckets::with_group(0, Default::default(), usize::MAX / 2);
        let mut largest = 0usize;
        let segment = 1u32 << SEGMENT_SHIFT;
        for (scan, &(rows, hits, first_s)) in [
            (5_000usize, 10usize, 0u32),
            (300_000, 9_000, 0),
            (300_000, 0, 0),
            (40_000, 3 * BLOCK_ENTRIES + 1, segment - 2_000),
            (600_000, 25_000, 0),
            (2_047, 2 * BLOCK_ENTRIES, 3 * segment - 700),
            (2_047, 5, u32::MAX - 4),
            (70_000, 6_000, 5 * segment + 12),
        ]
        .iter()
        .enumerate()
        {
            buckets.reset_for(rows);
            let cursor = buckets.open_group(0);
            let count = rows / RPS_BUCKET_SIZE + 1;
            assert_eq!(buckets.heads.len(), count, "scan {scan}");
            let mut expected: Vec<Vec<(u32, u32)>> = vec![Vec::new(); count];
            for i in 0..hits as u32 {
                // Words of one to three hits at increasing query offsets.
                let s_off = first_s.wrapping_add(i / 2);
                let q_off = lcg(&mut state) % rows.max(1) as u32;
                let s_bits = buckets.entry_bits(s_off);
                buckets.add(q_off, s_bits);
                expected[q_off as usize / RPS_BUCKET_SIZE].push((q_off, s_off));
            }
            for (b, hits) in expected.iter().enumerate() {
                assert_eq!(&buckets.bucket_hits(b), hits, "scan {scan} bucket {b}");
            }
            // Every closed block is full but where a segment started after it.
            let crossed = first_s >> SEGMENT_SHIFT
                != first_s.wrapping_add((hits as u32).saturating_sub(1) / 2) >> SEGMENT_SHIFT;
            for chain in &buckets.chains {
                for closed in chain {
                    assert!(
                        crossed || closed.len as usize == BLOCK_ENTRIES,
                        "scan {scan}"
                    );
                }
            }
            let held = buckets.heads.len() + buckets.chains.iter().map(Vec::len).sum::<usize>();
            largest = largest.max(held);
            let read = buckets.read_group(cursor);
            assert_eq!(read, expected.concat(), "scan {scan}");
            assert!(buckets.chains.iter().all(Vec::is_empty), "scan {scan}");
            assert!(
                buckets.heads.iter().all(|head| head.fill == 0),
                "scan {scan}"
            );
            // The pool holds what the buckets held at most at one time.
            assert_eq!(buckets.blocks(), largest, "scan {scan}");
            assert_eq!(
                buckets.kept_bytes(),
                buckets.blocks() * BLOCK_ENTRIES * std::mem::size_of::<u32>()
            );
        }
    }

    // M1 of PS12s: a scan whose query offsets cross a segment of 2^21 gives every hit back whole,
    // bucket by bucket in the order of the scan, as NCBI's pairs, whatever the groups (PS12m; the
    // cells list their hits by bucket, as a table made by makeprofiledb does). Ten buckets and
    // groups of 1 to 4: a word whose hits lie in a later group is skipped in the earlier ones,
    // also after the segment changed.
    #[test]
    fn a_scan_across_a_segment_reads_the_query_offsets_whole() {
        let cells: [([u8; 3], &[i32]); 3] = [
            ([1, 2, 3], &[10, 4100, 12_000]),
            ([2, 3, 1], &[20, 30, 2050, 2060, 5000]),
            ([3, 1, 2], &[18_000, 18_100]),
        ];
        let mut table = Table::new();
        for (letters, entries) in cells {
            table.set(word(letters), entries);
        }
        let boundary = 1usize << SEGMENT_SHIFT;
        let mut subject = vec![20u8; boundary + 12];
        for start in [
            boundary - 7,
            boundary - 4,
            boundary - 1,
            boundary + 2,
            boundary + 6,
        ] {
            subject[start..start + 3].copy_from_slice(&[1, 2, 3]);
            subject[start + 3] = 1;
        }
        subject[boundary - 6..boundary - 3].copy_from_slice(&[3, 1, 2]);
        let mut expected: Vec<Vec<(u32, u32)>> = vec![Vec::new(); 10];
        for s in 0..=subject.len() - 3 {
            let w = (usize::from(subject[s]) << 10)
                | (usize::from(subject[s + 1]) << 5)
                | usize::from(subject[s + 2]);
            let entries = cells
                .iter()
                .find(|(letters, _)| word(*letters) == w)
                .map_or(&[][..], |(_, entries)| entries);
            for &entry in entries {
                let q_off = (entry - 2) as u32;
                expected[q_off as usize / RPS_BUCKET_SIZE].push((q_off, s as u32));
            }
        }
        assert!(expected[0].iter().any(|&(_, s)| s as usize >= boundary));
        assert!(expected[8].iter().any(|&(_, s)| s as usize >= boundary - 7));
        for group_buckets in [RPS_GROUP_BUCKETS, 1, 2, 3, 4] {
            let mut buckets = RpsBuckets::with_group(20_000, Default::default(), group_buckets);
            let mut offset = 0;
            let hits = scan_subject(
                &table.view(),
                &subject,
                &mut offset,
                &mut buckets,
                MAX_SCAN_HITS,
            );
            assert_eq!(hits as usize, expected.iter().map(Vec::len).sum::<usize>());
            let read = buckets.drain(&table.view(), &subject);
            assert_eq!(read, expected.concat(), "{group_buckets} buckets a group");
        }
    }

    // M2 of PS12s: the buckets of two workers share one pool: the blocks one worker's buckets
    // give back once read serve the other's scan, so the run holds the blocks its buckets hold
    // at the same time.
    #[test]
    fn workers_share_the_blocks_given_back_once_read() {
        let pool: std::sync::Arc<RpsBlockPool> = Default::default();
        let mut a = RpsBuckets::new(5_000, pool.clone());
        let mut b = RpsBuckets::new(5_000, pool.clone());
        let fill = |buckets: &mut RpsBuckets| {
            let cursor = buckets.open_group(0);
            for s_off in 0..5 * BLOCK_ENTRIES as u32 {
                let s_bits = buckets.entry_bits(s_off);
                buckets.add(s_off % 3_000, s_bits);
            }
            cursor
        };
        let closed = |buckets: &RpsBuckets| buckets.chains.iter().map(Vec::len).sum::<usize>();
        let cursor = fill(&mut a);
        let closed_a = closed(&a);
        assert!(closed_a > 0);
        assert_eq!(pool.len(), 0);
        assert_eq!(a.read_group(cursor).len(), 5 * BLOCK_ENTRIES);
        assert_eq!(pool.len(), closed_a);
        // The other worker's scan takes those blocks: no block more.
        let cursor = fill(&mut b);
        assert_eq!(closed(&b), closed_a);
        assert_eq!(pool.len(), 0);
        assert_eq!(b.read_group(cursor).len(), 5 * BLOCK_ENTRIES);
        assert_eq!(pool.len(), closed_a);
    }

    /// NCBI's model of one scan (`BlastRPSScanSubject` with `s_AddToRPSBucket`): for each bucket
    /// the pairs appended in order, the scan's total and the offset the next scan starts at. `lists`
    /// are the cells' hits by word (empty: no hit), in the order of the cell.
    fn reference_scan(
        lists: &[Vec<i32>],
        subject: &[u8],
        offset: &mut i32,
        max_hits: i32,
        num_buckets: usize,
    ) -> (i32, Vec<Vec<(u32, u32)>>) {
        let mut buckets = vec![Vec::new(); num_buckets];
        let mut total = 0i32;
        let mut s = *offset;
        while s as usize + 3 <= subject.len() {
            let at = s as usize;
            let w = word([subject[at], subject[at + 1], subject[at + 2]]);
            let list = &lists[w];
            if !list.is_empty() {
                if list.len() as i32 > max_hits - total {
                    break;
                }
                for &entry in list {
                    let q_off = (entry - 2) as u32;
                    buckets[q_off as usize / RPS_BUCKET_SIZE].push((q_off, s as u32));
                }
                total += list.len() as i32;
            }
            s += 1;
        }
        *offset = s;
        (total, buckets)
    }

    /// A random table of `rows` rows over the letters 1 to 4: a hit list for about half of the
    /// words, of 1, 2, 3, 4 (the form boundary), 50 or (twice at most) 6,000 entries in ascending
    /// order, with equal buckets and repeated entries too.
    fn random_lists(state: &mut u64, rows: usize) -> Vec<Vec<i32>> {
        let mut lists = vec![Vec::new(); crate::seqdb::rps::RPS_BACKBONE_SIZE];
        let mut huge = 0;
        for a in 1..=4u8 {
            for b in 1..=4u8 {
                for c in 1..=4u8 {
                    if lcg(state).is_multiple_of(2) {
                        continue;
                    }
                    let mut len = [1usize, 2, 3, 4, 50, 3, 4, 2, 6_000][lcg(state) as usize % 9];
                    if len == 6_000 {
                        huge += 1;
                        if huge > 2 {
                            len = 50;
                        }
                    }
                    let mut list: Vec<i32> = (0..len)
                        .map(|_| 2 + (lcg(state) as usize % (rows - 2)) as i32)
                        .collect();
                    list.sort_unstable();
                    lists[word([a, b, c])] = list;
                }
            }
        }
        lists
    }

    fn table_of(lists: &[Vec<i32>]) -> Table {
        let mut table = Table::new();
        for (w, list) in lists.iter().enumerate() {
            if !list.is_empty() {
                table.set(w, list);
            }
        }
        table
    }

    /// Scans `subject` completely (scan after scan) into buckets of groups of `group_buckets` and
    /// asserts, for each scan, NCBI's total, next offset and hits per bucket in order.
    fn assert_scans_equal_ncbi(
        table: &RpsTable<'_>,
        lists: &[Vec<i32>],
        subject: &[u8],
        rows: usize,
        max_hits: i32,
        group_buckets: usize,
    ) -> usize {
        let num_buckets = rows / RPS_BUCKET_SIZE + 1;
        let mut buckets = RpsBuckets::with_group(rows, Default::default(), group_buckets);
        let (mut offset, mut ref_offset, mut scans) = (0i32, 0i32, 0);
        while offset as usize + 3 <= subject.len() {
            let total = scan_subject(table, subject, &mut offset, &mut buckets, max_hits);
            let (ref_total, expected) =
                reference_scan(lists, subject, &mut ref_offset, max_hits, num_buckets);
            assert_eq!((total, offset), (ref_total, ref_offset));
            if group_buckets >= num_buckets {
                // One group: the scan made every hit, as the scan of dev did, and no pass follows.
                assert!(!buckets.beyond);
                assert_eq!(buckets.groups(), 1);
            }
            assert_eq!(
                buckets.drain(table, subject),
                expected.concat(),
                "scan {scans} of {group_buckets} buckets a group"
            );
            scans += 1;
        }
        scans
    }

    // PS12m (second review, S2): a scan whose every word has one hit past the first group keeps
    // nothing per word: the later group is made by a pass over the subject, and the hits are the
    // scan's, in its order (the dev scan's state is the same hits).
    #[test]
    fn a_repeated_word_with_one_hit_past_the_first_group_needs_no_state_per_word() {
        let mut table = Table::new();
        let q_off = 100 * RPS_BUCKET_SIZE + 77; // bucket 100: group 1
        table.set(word([1, 1, 1]), &[q_off as i32 + 2]);
        let subject = vec![1u8; 20_000];
        let mut buckets = RpsBuckets::new(300_000, Default::default());
        let mut offset = 0;
        let hits = scan_subject(
            &table.view(),
            &subject,
            &mut offset,
            &mut buckets,
            MAX_SCAN_HITS,
        );
        assert_eq!(hits, 20_000 - 2);
        assert_eq!(buckets.groups(), 3);
        let read = buckets.drain(&table.view(), &subject);
        let expected: Vec<(u32, u32)> = (0..20_000 - 2).map(|s| (q_off as u32, s)).collect();
        assert_eq!(read, expected);
    }

    // PS12m: whatever the group of buckets (1, 2, 3, 7, 64, or all), every scan gives the word
    // finder NCBI's hits in NCBI's order: the same total and next offset, and the concatenation
    // of the buckets' pairs. Random tables (cells of 1, 2, 3, 4, 50 and 6,000 hits, up to 342
    // buckets), random subjects, and `max_hits` from the longest cell to a few times it, so
    // that a subject takes several scans and a scan ends inside a run of heavy words.
    #[test]
    fn groups_of_buckets_read_the_hits_in_ncbis_order() {
        let mut state = 7u64;
        let mut cut_scans = 0;
        // 131,071 and 131,074 rows: volumes of exactly 64 and of 65 buckets, the production
        // group's edge (one group, and a second group of one bucket).
        for case in 0..28 {
            let rows = [
                1_000usize, 5_000, 40_000, 140_000, 700_000, 131_071, 131_074,
            ][case % 7];
            let lists = random_lists(&mut state, rows);
            let table = table_of(&lists);
            let longest = lists.iter().map(Vec::len).max().unwrap() as i32;
            let length = 20 + lcg(&mut state) as usize % 1_480;
            let subject: Vec<u8> = (0..length)
                .map(|_| 1 + (lcg(&mut state) % 4) as u8)
                .collect();
            let max_hits = longest + (lcg(&mut state) % 40_000) as i32;
            let buckets = rows / RPS_BUCKET_SIZE + 1;
            for group_buckets in [1, 2, 3, 7, RPS_GROUP_BUCKETS, buckets, buckets + 5] {
                let scans = assert_scans_equal_ncbi(
                    &table.view(),
                    &lists,
                    &subject,
                    rows,
                    max_hits,
                    group_buckets,
                );
                cut_scans += scans.saturating_sub(1);
            }
        }
        assert!(
            cut_scans > 20,
            "{cut_scans} scans ended before the subject did"
        );
    }

    // PS12m: cells whose hits are not in groups of buckets in order are regrouped when the volume
    // is opened (`seqdb::rps`); read from the regrouped cells, every group gives the hits that
    // NCBI's scan of the cells as they are gives (the model is built from the original, shuffled
    // lists). The regrouping is by the group size read (production: `RPS_GROUP_BUCKETS`).
    #[test]
    fn regrouped_cells_read_the_hits_in_the_order_of_the_original_cells() {
        let mut state = 99u64;
        let mut regrouped_cases = 0;
        for case in 0..10 {
            let rows = [5_000usize, 40_000, 140_000, 700_000, 300_000][case % 5];
            let mut lists = random_lists(&mut state, rows);
            for list in lists.iter_mut().filter(|list| list.len() > 1) {
                // A shuffle (Fisher-Yates), so that most long lists lose their order.
                for i in (1..list.len()).rev() {
                    list.swap(i, lcg(&mut state) as usize % (i + 1));
                }
            }
            let table = table_of(&lists);
            let length = 20 + lcg(&mut state) as usize % 1_000;
            let subject: Vec<u8> = (0..length)
                .map(|_| 1 + (lcg(&mut state) % 4) as u8)
                .collect();
            let longest = lists.iter().map(Vec::len).max().unwrap() as i32;
            let buckets = rows / RPS_BUCKET_SIZE + 1;
            for group_buckets in [1, 3, RPS_GROUP_BUCKETS, buckets] {
                let regrouped =
                    crate::seqdb::rps::regroup_cells(&table.cells, &table.overflow, group_buckets)
                        .unwrap();
                if group_buckets < buckets {
                    assert!(
                        regrouped.is_some(),
                        "shuffled lists in {group_buckets} of {buckets} buckets"
                    );
                    regrouped_cases += 1;
                }
                let view = match &regrouped {
                    Some((cells, overflow)) => RpsTable {
                        cells,
                        overflow,
                        pv: &table.pv,
                    },
                    None => table.view(),
                };
                assert_scans_equal_ncbi(
                    &view,
                    &lists,
                    &subject,
                    rows,
                    longest + 3_000,
                    group_buckets,
                );
            }
        }
        assert!(regrouped_cases >= 10, "{regrouped_cases} regrouped cases");
    }
}
