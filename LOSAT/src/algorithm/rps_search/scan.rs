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

/// One word hit: the concatenated position of the word's first letter and the query offset
/// (`BlastOffsetPair`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct RpsOffsetPair {
    pub(crate) q_off: u32,
    pub(crate) s_off: u32,
}

/// The pairs of one block of a bucket (`RpsBuckets`): 8 KB.
pub(crate) const BLOCK_PAIRS: usize = 1024;

/// One block of a bucket's pairs.
type Block = Box<[RpsOffsetPair; BLOCK_PAIRS]>;

/// The last block of a bucket and the pairs it holds; `fill` is `BLOCK_PAIRS` for a bucket
/// without a block, so that its first pair takes one.
#[derive(Clone, Copy)]
struct BucketHead {
    block: u32,
    fill: u32,
}

const EMPTY_HEAD: BucketHead = BucketHead {
    block: 0,
    fill: BLOCK_PAIRS as u32,
};

/// The buckets of a volume (`lookup->bucket_array`, `rows / RPS_BUCKET_SIZE + 1` of them), kept
/// by a worker for its run (series PR, item PS12s, D4): a bucket is a list of fixed-size blocks
/// taken from the worker's pool of blocks, its pairs in the order of the scan; a scan gives every
/// block back to the pool. The pool keeps the blocks of the largest scan the worker made and at
/// most one partly filled block per bucket, where buckets of their own (NCBI's arrays, a `Vec`
/// each) each keep their largest size. The pairs are read bucket by bucket, each in the order of
/// the scan (`slices`), as NCBI reads its buckets.
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
    /// Each bucket's blocks (indexes into `blocks`) in the order they were filled.
    chains: Vec<Vec<u32>>,
    /// The pool: the first `used` blocks belong to the buckets of the current scan.
    blocks: Vec<Block>,
    used: usize,
}

impl RpsBuckets {
    /// The buckets of a volume of `num_rows` concatenated rows, without a block.
    pub(crate) fn new(num_rows: usize) -> Self {
        let mut buckets = Self {
            heads: Vec::new(),
            chains: Vec::new(),
            blocks: Vec::new(),
            used: 0,
        };
        buckets.reset_for(num_rows);
        buckets
    }

    /// The buckets of a volume of `num_rows` concatenated rows, emptied; the pool's blocks stay.
    pub(crate) fn reset_for(&mut self, num_rows: usize) {
        let count = num_rows / RPS_BUCKET_SIZE + 1;
        self.chains.resize_with(count, Vec::new);
        self.chains.truncate(count);
        self.heads.resize(count, EMPTY_HEAD);
        self.heads.truncate(count);
        self.clear();
    }

    /// Every bucket emptied, every block back in the pool.
    fn clear(&mut self) {
        for chain in &mut self.chains {
            chain.clear();
        }
        self.heads.fill(EMPTY_HEAD);
        self.used = 0;
    }

    /// `s_AddToRPSBucket`: the pair at the end of the bucket of its row.
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
    fn add(&mut self, q_off: u32, s_off: u32) {
        let bucket = q_off as usize / RPS_BUCKET_SIZE;
        let mut head = self.heads[bucket];
        if head.fill as usize >= BLOCK_PAIRS {
            head = self.take_block(bucket);
        }
        self.blocks[head.block as usize][head.fill as usize] = RpsOffsetPair { q_off, s_off };
        self.heads[bucket].fill = head.fill + 1;
    }

    /// A block from the pool (a new one when every block is in use) at the end of `bucket`.
    #[inline(never)]
    fn take_block(&mut self, bucket: usize) -> BucketHead {
        if self.used == self.blocks.len() {
            self.blocks.push(Box::new(
                [RpsOffsetPair { q_off: 0, s_off: 0 }; BLOCK_PAIRS],
            ));
        }
        let block = u32::try_from(self.used).expect("fewer than 2^32 blocks");
        self.used += 1;
        self.chains[bucket].push(block);
        let head = BucketHead { block, fill: 0 };
        self.heads[bucket] = head;
        head
    }

    /// The pairs of the current scan, bucket by bucket, each bucket's in the order of the scan,
    /// as slices of its blocks.
    pub(crate) fn slices(&self) -> impl Iterator<Item = &[RpsOffsetPair]> + '_ {
        self.chains
            .iter()
            .zip(&self.heads)
            .flat_map(move |(chain, head)| {
                let last = chain.len().wrapping_sub(1);
                chain.iter().enumerate().map(move |(i, &block)| {
                    let filled = if i == last {
                        head.fill as usize
                    } else {
                        BLOCK_PAIRS
                    };
                    &self.blocks[block as usize][..filled]
                })
            })
    }

    /// The bytes of the pool's blocks (what the buckets keep).
    #[cfg(test)]
    pub(crate) fn kept_bytes(&self) -> usize {
        self.blocks.len() * BLOCK_PAIRS * std::mem::size_of::<RpsOffsetPair>()
    }
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
                let s_off = s as u32;
                if numhits <= RPS_HITS_PER_CELL {
                    for &entry in &cell[1..=numhits as usize] {
                        let q_off = (entry - table_correction) as u32;
                        buckets.add(q_off, s_off);
                    }
                } else {
                    let q_off = (cell[1] - table_correction) as u32;
                    buckets.add(q_off, s_off);
                    let start = cell[2] as usize / std::mem::size_of::<i32>();
                    for &entry in &table.overflow[start..start + numhits as usize - 1] {
                        let q_off = (entry - table_correction) as u32;
                        buckets.add(q_off, s_off);
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
        let pairs = |b: usize| -> Vec<(u32, u32)> {
            bucket(&buckets, b)
                .iter()
                .map(|p| (p.q_off, p.s_off))
                .collect()
        };
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
        assert_eq!(bucket(&buckets, 0).len(), 6);
        let second = scan_subject(&table.view(), &subject, &mut offset, &mut buckets, 7);
        assert_eq!((second, offset), (6, 4));
        assert_eq!(bucket(&buckets, 0)[0], RpsOffsetPair { q_off: 0, s_off: 2 });
        let last = scan_subject(&table.view(), &subject, &mut offset, &mut buckets, 7);
        assert_eq!((last, offset), (0, 4));
    }

    /// The pairs of bucket `b` of the current scan, in order.
    fn bucket(buckets: &RpsBuckets, b: usize) -> Vec<RpsOffsetPair> {
        let head = buckets.heads[b];
        let chain = &buckets.chains[b];
        let mut pairs = Vec::new();
        for (i, &block) in chain.iter().enumerate() {
            let filled = if i + 1 == chain.len() {
                head.fill as usize
            } else {
                BLOCK_PAIRS
            };
            pairs.extend_from_slice(&buckets.blocks[block as usize][..filled]);
        }
        pairs
    }

    /// A deterministic sequence of pseudo-random numbers (no NCBI counterpart: test data).
    fn lcg(state: &mut u64) -> u32 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (*state >> 33) as u32
    }

    // D4 of PS12s: the blocks of a bucket hold its pairs in the order of the scan, and `slices`
    // reads them bucket by bucket, as one `Vec` per bucket (NCBI's arrays) holds and reads them;
    // the pool keeps its blocks from scan to scan and volume to volume (fewer or more buckets),
    // and holds what the largest scan took plus at most one partly filled block per bucket.
    #[test]
    fn blocks_read_the_pairs_in_the_order_of_vectors_per_bucket() {
        let mut state = 12_345u64;
        let mut buckets = RpsBuckets::new(0);
        let mut largest = 0usize;
        for (scan, &(rows, hits)) in [
            (5_000usize, 10usize),
            (300_000, 9_000),
            (300_000, 0),
            (40_000, 3 * BLOCK_PAIRS + 1),
            (600_000, 25_000),
            (2_047, 2 * BLOCK_PAIRS),
        ]
        .iter()
        .enumerate()
        {
            buckets.reset_for(rows);
            let count = rows / RPS_BUCKET_SIZE + 1;
            assert_eq!(buckets.heads.len(), count, "scan {scan}");
            let mut expected: Vec<Vec<RpsOffsetPair>> = vec![Vec::new(); count];
            for s_off in 0..hits as u32 {
                let q_off = lcg(&mut state) % rows.max(1) as u32;
                buckets.add(q_off, s_off);
                expected[q_off as usize / RPS_BUCKET_SIZE].push(RpsOffsetPair { q_off, s_off });
            }
            let read: Vec<RpsOffsetPair> = buckets.slices().flatten().copied().collect();
            assert_eq!(read, expected.concat(), "scan {scan}");
            for (b, pairs) in expected.iter().enumerate() {
                assert_eq!(&bucket(&buckets, b), pairs, "scan {scan} bucket {b}");
            }
            let partly_filled = expected
                .iter()
                .filter(|pairs| pairs.len() % BLOCK_PAIRS != 0)
                .count();
            let taken: usize = expected
                .iter()
                .map(|pairs| pairs.len().div_ceil(BLOCK_PAIRS))
                .sum();
            assert_eq!(buckets.used, taken, "scan {scan}");
            assert!(taken <= hits / BLOCK_PAIRS + partly_filled, "scan {scan}");
            largest = largest.max(taken);
            assert_eq!(buckets.blocks.len(), largest, "scan {scan}");
            assert_eq!(
                buckets.kept_bytes(),
                largest * BLOCK_PAIRS * std::mem::size_of::<RpsOffsetPair>()
            );
        }
    }
}
