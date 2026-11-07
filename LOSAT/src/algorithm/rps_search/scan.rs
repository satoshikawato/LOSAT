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

/// The buckets of a volume (`lookup->bucket_array`, `rows / RPS_BUCKET_SIZE + 1` of them).
pub(crate) struct RpsBuckets {
    pub(crate) buckets: Vec<Vec<RpsOffsetPair>>,
}

impl RpsBuckets {
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
    pub(crate) fn new(num_rows: usize) -> Self {
        Self {
            buckets: (0..num_rows / RPS_BUCKET_SIZE + 1)
                .map(|_| Vec::new())
                .collect(),
        }
    }
}

/// `s_AddToRPSBucket` (the bucket grows as NCBI's doubles its array).
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
#[inline]
fn add_to_bucket(buckets: &mut [Vec<RpsOffsetPair>], q_off: u32, s_off: u32) {
    buckets[q_off as usize / RPS_BUCKET_SIZE].push(RpsOffsetPair { q_off, s_off });
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
    for bucket in &mut buckets.buckets {
        bucket.clear();
    }
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
                        add_to_bucket(&mut buckets.buckets, q_off, s_off);
                    }
                } else {
                    let q_off = (cell[1] - table_correction) as u32;
                    add_to_bucket(&mut buckets.buckets, q_off, s_off);
                    let start = cell[2] as usize / std::mem::size_of::<i32>();
                    for &entry in &table.overflow[start..start + numhits as usize - 1] {
                        let q_off = (entry - table_correction) as u32;
                        add_to_bucket(&mut buckets.buckets, q_off, s_off);
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
            buckets.buckets[b]
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
        assert_eq!(buckets.buckets[0].len(), 6);
        let second = scan_subject(&table.view(), &subject, &mut offset, &mut buckets, 7);
        assert_eq!((second, offset), (6, 4));
        assert_eq!(buckets.buckets[0][0], RpsOffsetPair { q_off: 0, s_off: 2 });
        let last = scan_subject(&table.view(), &subject, &mut offset, &mut buckets, 7);
        assert_eq!((last, offset), (0, 4));
    }
}
