//! The scan of a user query with the prebuilt lookup table of a profile volume
//! (`BlastRPSScanSubject`): the hits of each word go to the bucket of 2,048 concatenated rows
//! that holds them, in the order of the scan, and a scan stops before the word whose hits would
//! pass 4,000,000; the next scan restarts at that word.

use crate::seqdb::rps::{RPS_HITS_PER_CELL, RPS_MAX_SCAN_HITS, RPS_WORDSIZE};

/// `RPS_BUCKET_SIZE`.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h:349-351
/// ```c
/// /** The number of regions into which the concatenated RPS blast
///     database is split via bucket sorting */
/// #define RPS_BUCKET_SIZE 2048
/// ```
pub(crate) const RPS_BUCKET_SIZE: usize = 2048;

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
const BUCKET_BITS: u32 = 11;
const _: () = assert!(RPS_BUCKET_SIZE == 1 << BUCKET_BITS);
const BUCKET_MASK: u32 = (1 << BUCKET_BITS) - 1;

/// The bits of a query offset above those an entry holds (`RpsBuckets`): an entry keeps the
/// offset's low `SEGMENT_SHIFT` bits, its block the rest (the offset's segment).
const SEGMENT_SHIFT: u32 = u32::BITS - BUCKET_BITS;

/// The entries of one block of a bucket (`RpsBuckets`): 4 KB.
pub(crate) const BLOCK_ENTRIES: usize = 1024;

/// One block of a bucket's entries.
type Block = Box<[u32; BLOCK_ENTRIES]>;

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

/// The buckets of a volume (`lookup->bucket_array`, `rows / RPS_BUCKET_SIZE + 1` of them), kept
/// by a worker for its run (series PR, item PS12s, D4): a bucket is a list of fixed-size blocks,
/// its hits in the order of the scan; the hits are read bucket by bucket, each in the order of
/// the scan (`next_block`), as NCBI reads its buckets, and each bucket's blocks but its current one go
/// back to the pool of free blocks as soon as they are read.
///
/// A hit takes 4 bytes (M1 of PS12s), where NCBI's `BlastOffsetPair` takes 8: the concatenated
/// position inside its bucket (`BUCKET_BITS` bits; the bucket gives the rest) and the low
/// `SEGMENT_SHIFT` bits of the query offset, whose high bits (its segment) the block holds once.
/// A segment holds 2^21 query offsets; when the scan reaches the next one, every bucket's current
/// block that holds hits is closed (`start_segment`), so a block's hits are in one segment and
/// every offset reads back whole, whatever the length of the query.
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
    heads: Vec<BucketHead>,
    /// Each bucket's closed blocks in the order they were filled.
    chains: Vec<Vec<ClosedBlock>>,
    /// The query offsets' segment of the scan's current word (`entry_bits`).
    segment: u32,
    /// The blocks that no bucket holds.
    free: Vec<Block>,
}

impl RpsBuckets {
    /// The buckets of a volume of `num_rows` concatenated rows, empty.
    pub(crate) fn new(num_rows: usize) -> Self {
        let mut buckets = Self {
            heads: Vec::new(),
            chains: Vec::new(),
            segment: 0,
            free: Vec::new(),
        };
        buckets.reset_for(num_rows);
        buckets
    }

    /// The buckets of a volume of `num_rows` concatenated rows, emptied; the blocks stay.
    pub(crate) fn reset_for(&mut self, num_rows: usize) {
        let count = num_rows / RPS_BUCKET_SIZE + 1;
        self.clear();
        while self.heads.len() > count {
            let head = self.heads.pop().expect("a bucket above the count");
            self.free.push(head.block);
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
    }

    /// Every bucket emptied, its closed blocks back in the pool.
    fn clear(&mut self) {
        for (head, chain) in self.heads.iter_mut().zip(&mut self.chains) {
            self.free.extend(chain.drain(..).map(|closed| closed.block));
            head.fill = 0;
            head.segment = 0;
        }
        self.segment = 0;
    }

    /// A block of the pool (a new one when the pool has none).
    fn take(&mut self) -> Block {
        self.free
            .pop()
            .unwrap_or_else(|| Box::new([0; BLOCK_ENTRIES]))
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
    #[inline(always)]
    fn add(&mut self, q_off: u32, s_bits: u32) {
        let bucket = (q_off >> BUCKET_BITS) as usize;
        if self.heads[bucket].fill as usize >= BLOCK_ENTRIES {
            self.close(bucket);
        }
        let head = &mut self.heads[bucket];
        head.block[head.fill as usize] = s_bits | (q_off & BUCKET_MASK);
        head.fill += 1;
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

    /// The next block of hits of the current scan (`cursor`: where the reading stands; it starts
    /// at the default), bucket by bucket, each bucket's blocks in the order they were filled: every
    /// hit in the order NCBI reads its buckets. A bucket is emptied when the reading passes it, its
    /// closed blocks back in the pool.
    pub(crate) fn next_block(&mut self, cursor: &mut RpsReadCursor) -> Option<RpsHits<'_>> {
        while cursor.bucket < self.heads.len() {
            let bucket = cursor.bucket;
            let q_base = (bucket as u32) << BUCKET_BITS;
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
                self.free.extend(chain.drain(..).map(|closed| closed.block));
            }
            self.heads[bucket].fill = 0;
            cursor.bucket += 1;
            cursor.block = 0;
        }
        None
    }

    /// Every hit of the current scan as (`q_off`, `s_off`), in the order of `next_block`; the
    /// buckets emptied.
    #[cfg(test)]
    fn drain(&mut self) -> Vec<(u32, u32)> {
        let mut hits = Vec::new();
        let mut cursor = RpsReadCursor::default();
        while let Some(block) = self.next_block(&mut cursor) {
            hits.extend(block.pairs());
        }
        hits
    }

    /// The hits of `bucket` as (`q_off`, `s_off`), in order (not emptied).
    #[cfg(test)]
    fn bucket_hits(&self, bucket: usize) -> Vec<(u32, u32)> {
        let q_base = (bucket as u32) << BUCKET_BITS;
        let head = &self.heads[bucket];
        let blocks = self.chains[bucket]
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
        self.heads.len() + self.chains.iter().map(Vec::len).sum::<usize>() + self.free.len()
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

/// `BlastRPSScanSubject`: the hits of the words of `subject` (its letters, `subject->sequence`)
/// from `*offset` into the emptied buckets; `*offset` becomes the word the next scan starts at.
/// Returns the number of hits.
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
    buckets.clear();
    let wordsize = RPS_WORDSIZE as i32;
    let charsize = crate::seqdb::rps::RPS_CHARSIZE as i32;
    let mask = (crate::seqdb::rps::RPS_BACKBONE_SIZE - 1) as i32;
    let s_first = *offset;
    let s_last = subject.len() as i32 - wordsize;
    let table_correction = wordsize - 1;
    let mut totalhits = 0i32;
    // `ComputeTableIndex(wordsize - 1, charsize, s_first)`.
    let mut index = 0i32;
    for i in 0..wordsize - 1 {
        index = (index << charsize) | i32::from(subject[(s_first + i) as usize]);
    }
    let mut s = s_first;
    while s <= s_last {
        // `ComputeTableIndexIncremental`: the word's last letter.
        index = ((index << charsize) | i32::from(subject[(s + wordsize - 1) as usize])) & mask;
        let word = index as usize;
        if table.pv[word / 64] & (1u64 << (word % 64)) != 0 {
            let cell = &table.cells[word];
            let numhits = cell[0];
            if numhits <= max_hits - totalhits {
                let s_bits = buckets.entry_bits(s as u32);
                if numhits <= RPS_HITS_PER_CELL {
                    for &entry in &cell[1..=numhits as usize] {
                        let q_off = (entry - table_correction) as u32;
                        buckets.add(q_off, s_bits);
                    }
                } else {
                    let q_off = (cell[1] - table_correction) as u32;
                    buckets.add(q_off, s_bits);
                    let start = cell[2] as usize / std::mem::size_of::<i32>();
                    for &entry in &table.overflow[start..start + numhits as usize - 1] {
                        let q_off = (entry - table_correction) as u32;
                        buckets.add(q_off, s_bits);
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
        let mut buckets = RpsBuckets::new(6000);
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
        let mut buckets = RpsBuckets::new(10);
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
    // `next_block` reads them bucket by bucket, as one `Vec` of NCBI's pairs per bucket holds and reads
    // them; a drained bucket's closed blocks go back to the pool, which keeps its blocks from scan
    // to scan and volume to volume (fewer or more buckets): what the largest scan's buckets held.
    // Query offsets cross segments of 2^21 (also inside a bucket's block) and read back whole.
    #[test]
    fn blocks_read_the_hits_in_the_order_of_vectors_per_bucket() {
        let mut state = 12_345u64;
        let mut buckets = RpsBuckets::new(0);
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
            let read = buckets.drain();
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
    // bucket by bucket in the order of the scan, as NCBI's pairs.
    #[test]
    fn a_scan_across_a_segment_reads_the_query_offsets_whole() {
        let mut table = Table::new();
        table.set(word([1, 2, 3]), &[4100, 10]);
        table.set(word([2, 3, 1]), &[20, 2060, 30, 2050, 5000]);
        let boundary = 1usize << SEGMENT_SHIFT;
        let mut subject = vec![20u8; boundary + 8];
        for start in [boundary - 4, boundary - 1, boundary + 2] {
            subject[start..start + 3].copy_from_slice(&[1, 2, 3]);
            subject[start + 3] = 1;
        }
        let mut buckets = RpsBuckets::new(6000);
        let mut offset = 0;
        let hits = scan_subject(
            &table.view(),
            &subject,
            &mut offset,
            &mut buckets,
            MAX_SCAN_HITS,
        );
        let mut expected: Vec<Vec<(u32, u32)>> = vec![Vec::new(); 3];
        for s in 0..=subject.len() - 3 {
            let w = (usize::from(subject[s]) << 10)
                | (usize::from(subject[s + 1]) << 5)
                | usize::from(subject[s + 2]);
            let entries: &[i32] = if w == word([1, 2, 3]) {
                &[4100, 10]
            } else if w == word([2, 3, 1]) {
                &[20, 2060, 30, 2050, 5000]
            } else {
                &[]
            };
            for &entry in entries {
                let q_off = (entry - 2) as u32;
                expected[q_off as usize / RPS_BUCKET_SIZE].push((q_off, s as u32));
            }
        }
        assert_eq!(hits as usize, expected.iter().map(Vec::len).sum::<usize>());
        assert!(expected[0].iter().any(|&(_, s)| s as usize >= boundary));
        let read = buckets.drain();
        assert_eq!(read, expected.concat());
    }
}
