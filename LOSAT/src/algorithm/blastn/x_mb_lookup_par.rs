//! EXPERIMENT (LOSAT_X_PAIRPAR, shadow LOSAT_X_PAIRPARSHADOW): the megablast lookup table of
//! `build_mb_lookup` filled by the threads of the search pool, each thread owning a fixed set
//! of table cells.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1032-1033,1043-1045,1081-1091
//! ```c
//!       for (index = from; index <= last_offset; index++) {
//!          val = *++seq;
//! ...
//!          ecode = ((ecode << BITS_PER_NUC) & kLutMask) + val;
//!          if (seq < pos)
//!             continue;
//! ...
//!          if (mb_lt->hashtable[ecode] == 0) {
//! ...
//!             PV_SET(pv_array, ecode, pv_array_bts);
//!          }
//!          else {
//!             helper_array[ecode/kCompressionFactor]++;
//!          }
//!          mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
//!          mb_lt->hashtable[ecode] = index;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1104-1106
//! ```c
//!    longest_chain = 2;
//!    for (index = 0; index < mb_lt->hashsize / kCompressionFactor; index++)
//!        longest_chain = MAX(longest_chain, helper_array[index]);
//! ```
//!
//! Why the table is the same. The update above reads and writes only the cell `ecode` of
//! `hashtable`, the presence-vector word `ecode >> pv_array_bts`, the counter
//! `ecode / kCompressionFactor` and `next_pos[index]`, where `index` is the word's own
//! (unique) query offset. The cells are cut into blocks of `2^block_shift` cells, with
//! `2^block_shift` a multiple of both `kCompressionFactor` (2048) and `2^pv_array_bts`, so a
//! block owns whole presence-vector words and whole counters. Every worker reads the whole
//! query in the original order (the same loop, the same ambiguity resets and masks) and
//! applies the update only to the words whose cell lies in one of its blocks. Each cell
//! therefore receives the same words in the same order as in the serial loop, which fixes
//! its chain (`hashtable` and the `next_pos` links of its words), its presence bit and its
//! counter; `longest_chain` is then computed from the same counters by the original code.
//! The words of a worker are applied in reading order (a first-in first-out delay of 16
//! words on x86_64, for a prefetch, as LOSAT_X_MBDELAY does). Workers write disjoint
//! elements, so the number of workers changes the running time only.
//!
//! The re-threading of `ascending_cells` (blast_lookup.c:74-76 order, see `build_mb_lookup`)
//! touches the same cells and links of the same words, so each worker applies it to its own
//! words.

use super::{build_unmasked_ranges, MaskedInterval, PvArrayType, BLAST2NA_MASK, PV_ARRAY_MASK};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1081-1091
/// ```c
///          if (mb_lt->hashtable[ecode] == 0) {
///             PV_SET(pv_array, ecode, pv_array_bts);
///          }
///          else {
///             helper_array[ecode/kCompressionFactor]++;
///          }
/// ```
/// `kCompressionFactor`, the number of cells that share one counter.
const K_COMPRESSION_FACTOR_SHIFT: usize = 11;

/// No NCBI counterpart: reads the switches once; they do not change any value NCBI computes.
/// (LOSAT_X_PAIRPAR, LOSAT_X_PAIRPARSHADOW)
fn switches() -> (bool, bool) {
    static SWITCHES: OnceLock<(bool, bool)> = OnceLock::new();
    *SWITCHES.get_or_init(|| {
        (
            std::env::var_os("LOSAT_X_PAIRPAR").is_some(),
            std::env::var_os("LOSAT_X_PAIRPARSHADOW").is_some(),
        )
    })
}

/// No NCBI counterpart: shadow counters for the summary line (tables compared).
pub(crate) static SHADOW_TABLES: AtomicU64 = AtomicU64::new(0);
/// No NCBI counterpart: shadow counters for the summary line (cells compared).
pub(crate) static SHADOW_CELLS: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart: the number of threads of the search pool this thread runs in, when
/// there is one with two threads or more (scheduling only).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn pool_threads() -> Option<usize> {
    rayon::current_thread_index()?;
    let threads = rayon::current_num_threads();
    (threads > 1).then_some(threads)
}

#[cfg(any(
    not(feature = "parallel"),
    all(target_arch = "wasm32", not(feature = "wasm-threads"))
))]
fn pool_threads() -> Option<usize> {
    None
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1049-1051
/// ```c
///          if (kDbFilter) {
///              if (!(ecode & 1)) {
///                  if ((counts[ecode / 2] >> 4) >= max_word_count) {
/// ```
/// Dispatch test of `build_mb_lookup`: the number of workers when the switch is on, the
/// search runs on a pool of two threads or more, and the table is the plain one (`eligible`:
/// no database word counts, no debug counters). `None` runs the original loops.
pub(super) fn parts(eligible: bool) -> Option<usize> {
    let (on, shadow) = switches();
    if !on || shadow || !eligible {
        return None;
    }
    pool_threads()
}

/// No NCBI counterpart: LOSAT_X_PAIRPARSHADOW is set.
pub(super) fn shadow() -> bool {
    switches().1
}

/// The tables `build_mb_lookup` fills (NCBI `BlastMBLookupTable` fields and its local
/// `helper_array`).
pub(super) struct Tables<'a> {
    pub(super) hashtable: &'a mut [u32],
    pub(super) next_pos: &'a mut [u32],
    pub(super) pv_array: &'a mut [PvArrayType],
    pub(super) pv_array_bts: usize,
    pub(super) helper_array: &'a mut [u32],
}

/// The read-only inputs of the fill loop.
pub(super) struct Words<'a> {
    pub(super) queries_blastna: &'a [Vec<u8>],
    pub(super) query_offsets: &'a [i32],
    pub(super) query_masks: &'a [Vec<MaskedInterval>],
    pub(super) word_length: usize,
    pub(super) lut_word_length: usize,
    pub(super) kmer_mask: u64,
    pub(super) ascending_cells: bool,
}

// No NCBI counterpart: a slice shared by the workers, each of which reads and writes only the
// elements it owns (see the module documentation); the pointer outlives the scope.
#[derive(Clone, Copy)]
struct Shared<T> {
    ptr: *mut T,
    len: usize,
}

// SAFETY: the workers access disjoint elements only (ownership by cell block, and by the
// word's own query offset for `next_pos`), and the slices outlive the scope that uses them.
unsafe impl<T: Send> Send for Shared<T> {}
// SAFETY: as above.
unsafe impl<T: Send> Sync for Shared<T> {}

impl<T: Copy> Shared<T> {
    fn new(slice: &mut [T]) -> Self {
        Self {
            ptr: slice.as_mut_ptr(),
            len: slice.len(),
        }
    }

    /// SAFETY: `i < len`, and no other worker accesses element `i`.
    #[inline(always)]
    unsafe fn get(self, i: usize) -> T {
        debug_assert!(i < self.len);
        unsafe { *self.ptr.add(i) }
    }

    /// SAFETY: `i < len`, and no other worker accesses element `i`.
    #[inline(always)]
    unsafe fn set(self, i: usize, value: T) {
        debug_assert!(i < self.len);
        unsafe { *self.ptr.add(i) = value }
    }
}

#[derive(Clone, Copy)]
struct SharedTables {
    hashtable: Shared<u32>,
    next_pos: Shared<u32>,
    pv_array: Shared<PvArrayType>,
    pv_array_bts: usize,
    helper_array: Shared<u32>,
}

impl SharedTables {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1081-1091
    // ```c
    //          if (mb_lt->hashtable[ecode] == 0) {
    //             PV_SET(pv_array, ecode, pv_array_bts);
    //          }
    //          else {
    //             helper_array[ecode/kCompressionFactor]++;
    //          }
    //          mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
    //          mb_lt->hashtable[ecode] = index;
    // ```
    // The update of one word, as in the loop of `build_mb_lookup` (same statements);
    // the presence bit is set in the presence-vector word this worker owns.
    /// SAFETY: the worker owns the block of `bucket` and the word at `q_off_1`.
    #[inline(always)]
    unsafe fn apply(self, bucket: usize, q_off_1: u32) {
        unsafe {
            let head = self.hashtable.get(bucket);
            if head == 0 {
                // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:49-50
                // ```c
                // #define PV_SET(lookup, index, shift) \
                //     lookup[(index) >> (shift)] |= (PV_ARRAY_TYPE)1 << ((index) & PV_ARRAY_MASK)
                // ```
                // As `pv_set_shift`, on the presence-vector word this worker owns.
                let array_idx = bucket >> self.pv_array_bts;
                if array_idx < self.pv_array.len {
                    let bit = (1 as PvArrayType) << (bucket & PV_ARRAY_MASK);
                    self.pv_array
                        .set(array_idx, self.pv_array.get(array_idx) | bit);
                }
            } else {
                let h = bucket >> K_COMPRESSION_FACTOR_SHIFT;
                self.helper_array
                    .set(h, self.helper_array.get(h).saturating_add(1));
            }
            self.next_pos.set(q_off_1 as usize, head);
            self.hashtable.set(bucket, q_off_1);
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1032-1045
/// ```c
///       for (index = from; index <= last_offset; index++) {
///          val = *++seq;
///          /* if an ambiguity is encountered, do not add
///             any words that would contain it */
///          if ((val & BLAST2NA_MASK) != 0) {
///             ecode = 0;
///             pos = seq + kLutWordLength;
///             continue;
///          }
///
///          /* get next base */
///          ecode = ((ecode << BITS_PER_NUC) & kLutMask) + val;
///          if (seq < pos)
///             continue;
/// ```
/// One worker: the word loop of `build_mb_lookup` (same reading order, same resets), applying
/// the update only to the words whose cell block `owner` gives to `part`.
fn fill_part(part: u8, owner: &[u8], block_shift: usize, words: &Words<'_>, t: SharedTables) {
    let lut_word_length = words.lut_word_length;
    let mut own_words: Vec<(u32, u32)> = Vec::new();
    #[cfg(target_arch = "x86_64")]
    const X_DELAY: usize = 16;
    #[cfg(target_arch = "x86_64")]
    let mut x_ring = [(0u32, 0u32); X_DELAY];
    #[cfg(target_arch = "x86_64")]
    let mut x_ring_len = 0usize;
    for (q_idx, seq_blastna) in words.queries_blastna.iter().enumerate() {
        let seq = seq_blastna.as_slice();
        if seq.len() < lut_word_length {
            continue;
        }

        let masks = words
            .query_masks
            .get(q_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);
        let ranges = build_unmasked_ranges(seq.len(), masks);
        let query_offset = words.query_offsets[q_idx].max(0) as usize;

        for (range_start, range_end) in ranges {
            let range_len = range_end.saturating_sub(range_start);
            if words.word_length > range_len {
                continue;
            }

            let mut current_kmer: u64 = 0;
            let mut valid_bases = 0usize;

            for pos in range_start..range_end {
                let base = seq[pos];
                if (base & BLAST2NA_MASK) != 0 {
                    current_kmer = 0;
                    valid_bases = 0;
                    continue;
                }

                current_kmer = ((current_kmer << 2) | base as u64) & words.kmer_mask;
                valid_bases += 1;
                if valid_bases < lut_word_length {
                    continue;
                }

                let bucket = current_kmer as usize;
                if owner[bucket >> block_shift] != part {
                    continue;
                }
                let q_off_1 = (query_offset + (pos + 1 - lut_word_length) + 1) as u32;
                if words.ascending_cells {
                    own_words.push((bucket as u32, q_off_1));
                }
                #[cfg(target_arch = "x86_64")]
                {
                    // SAFETY: `bucket` < table size (masked k-mer); a prefetch does not fault.
                    unsafe {
                        core::arch::x86_64::_mm_prefetch(
                            t.hashtable.ptr.add(bucket) as *const i8,
                            core::arch::x86_64::_MM_HINT_T0,
                        );
                    }
                    let slot = x_ring_len % X_DELAY;
                    if x_ring_len >= X_DELAY {
                        let (b, q) = x_ring[slot];
                        // SAFETY: the cell block of `b` is owned by this part; `q` is the
                        // offset of one of this part's words.
                        unsafe { t.apply(b as usize, q) };
                    }
                    x_ring[slot] = (bucket as u32, q_off_1);
                    x_ring_len += 1;
                }
                #[cfg(not(target_arch = "x86_64"))]
                {
                    // SAFETY: as above.
                    unsafe { t.apply(bucket, q_off_1) };
                }
            }
        }
    }
    #[cfg(target_arch = "x86_64")]
    {
        let pending = x_ring_len.min(X_DELAY);
        for k in x_ring_len - pending..x_ring_len {
            let (b, q) = x_ring[k % X_DELAY];
            // SAFETY: as above.
            unsafe { t.apply(b as usize, q) };
        }
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_lookup.c:74-76
    // ```c
    // /* add the hit */
    // chain[chain[1] + 2] = query_offset;
    // chain[1]++;
    // ```
    // The re-threading of `build_mb_lookup` for `ascending_cells`, on this part's words.
    if words.ascending_cells {
        for &(bucket, _) in &own_words {
            // SAFETY: own cell.
            unsafe { t.hashtable.set(bucket as usize, 0) };
        }
        for &(bucket, q_off_1) in own_words.iter().rev() {
            // SAFETY: own cell and own word.
            unsafe {
                t.next_pos
                    .set(q_off_1 as usize, t.hashtable.get(bucket as usize));
                t.hashtable.set(bucket as usize, q_off_1);
            }
        }
    }
}

/// No NCBI counterpart: the cell blocks and their owners (scheduling only). Blocks of
/// `2^block_shift` cells, dealt round robin to `parts` workers.
fn plan(table_size: usize, pv_array_bts: usize, parts: usize) -> (usize, Vec<u8>) {
    let block_shift = K_COMPRESSION_FACTOR_SHIFT.max(pv_array_bts);
    let blocks = (table_size >> block_shift).max(1);
    let parts = parts.clamp(1, 255);
    let owner = (0..blocks).map(|block| (block % parts) as u8).collect();
    (block_shift, owner)
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1032-1091
/// ```c
///       for (index = from; index <= last_offset; index++) {
/// ...
///          mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
///          mb_lt->hashtable[ecode] = index;
/// ```
/// Fills the tables as the word loop of `build_mb_lookup` does, with `parts` workers on the
/// search pool (or one after another, `in_pool == false`, for the shadow without a pool).
pub(super) fn fill(parts: usize, words: &Words<'_>, tables: Tables<'_>, in_pool: bool) {
    let table_size = tables.hashtable.len();
    let (block_shift, owner) = plan(table_size, tables.pv_array_bts, parts);
    let parts = parts.clamp(1, 255);
    let shared = SharedTables {
        hashtable: Shared::new(tables.hashtable),
        next_pos: Shared::new(tables.next_pos),
        pv_array: Shared::new(tables.pv_array),
        pv_array_bts: tables.pv_array_bts,
        helper_array: Shared::new(tables.helper_array),
    };
    let owner = owner.as_slice();
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if in_pool {
        rayon::scope(|scope| {
            for part in 1..parts {
                scope.spawn(move |_| fill_part(part as u8, owner, block_shift, words, shared));
            }
            fill_part(0, owner, block_shift, words, shared);
        });
        return;
    }
    let _ = in_pool;
    for part in 0..parts {
        fill_part(part as u8, owner, block_shift, words, shared);
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1081-1091,1104-1106
/// ```c
///          mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
///          mb_lt->hashtable[ecode] = index;
/// ...
///    longest_chain = 2;
///    for (index = 0; index < mb_lt->hashsize / kCompressionFactor; index++)
///        longest_chain = MAX(longest_chain, helper_array[index]);
/// ```
/// LOSAT_X_PAIRPARSHADOW: builds the tables again with the partitioned fill (on the pool when
/// there is one, otherwise with four parts one after another) and stops at the first element
/// that differs from the tables the original loop built.
#[allow(clippy::too_many_arguments)]
pub(super) fn shadow_check(
    words: &Words<'_>,
    pv_array_bts: usize,
    hashtable: &[u32],
    next_pos: &[u32],
    pv_array: &[PvArrayType],
    helper_array: &[u32],
    longest_chain: usize,
) {
    let mut x_hashtable = vec![0u32; hashtable.len()];
    let mut x_next_pos = vec![0u32; next_pos.len()];
    let mut x_pv_array = vec![0 as PvArrayType; pv_array.len()];
    let mut x_helper_array = vec![0u32; helper_array.len()];
    let (parts, in_pool) = match pool_threads() {
        Some(threads) => (threads, true),
        None => (4, false),
    };
    fill(
        parts,
        words,
        Tables {
            hashtable: &mut x_hashtable,
            next_pos: &mut x_next_pos,
            pv_array: &mut x_pv_array,
            pv_array_bts,
            helper_array: &mut x_helper_array,
        },
        in_pool,
    );
    let mut x_longest_chain = 2usize;
    for &value in &x_helper_array {
        x_longest_chain = x_longest_chain.max(value as usize);
    }
    for (name, reference, partitioned) in [
        ("hashtable", hashtable, x_hashtable.as_slice()),
        ("next_pos", next_pos, x_next_pos.as_slice()),
        ("pv_array", pv_array, x_pv_array.as_slice()),
        ("helper_array", helper_array, x_helper_array.as_slice()),
    ] {
        if let Some(i) = (0..reference.len()).find(|&i| reference[i] != partitioned[i]) {
            panic!(
                "LOSAT_X_PAIRPARSHADOW: megablast lookup {name}[{i}] differs: reference {} partitioned {} (parts {parts}, lut_word_length {}, ascending_cells {})",
                reference[i], partitioned[i], words.lut_word_length, words.ascending_cells
            );
        }
    }
    assert_eq!(
        longest_chain, x_longest_chain,
        "LOSAT_X_PAIRPARSHADOW: megablast lookup longest_chain differs (parts {parts})"
    );
    SHADOW_TABLES.fetch_add(1, Ordering::Relaxed);
    SHADOW_CELLS.fetch_add(
        (hashtable.len() + next_pos.len() + pv_array.len() + helper_array.len()) as u64,
        Ordering::Relaxed,
    );
}

#[cfg(test)]
mod tests {
    use super::super::build_two_stage_lookup;
    use super::*;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1081-1091
    // ```c
    //          mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
    //          mb_lt->hashtable[ecode] = index;
    // ```
    // Random queries (with ambiguity codes and masks) through the original loop and through
    // the partitioned fill with 1..8 workers: every table element must be the same.
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n.max(1)
        }
    }

    fn random_case(
        rng: &mut Lcg,
    ) -> (
        Vec<Vec<u8>>,
        Vec<i32>,
        Vec<Vec<MaskedInterval>>,
        usize,
        usize,
    ) {
        let lut_word_length = [4usize, 8, 9, 11, 12][rng.below(5) as usize];
        let word_length = lut_word_length + rng.below(17) as usize;
        let n_queries = 1 + rng.below(3) as usize;
        let mut queries = Vec::new();
        let mut offsets = Vec::new();
        let mut masks = Vec::new();
        let mut offset = 0i32;
        for _ in 0..n_queries {
            let len = 1 + rng.below(3000) as usize;
            // a small alphabet of repeats makes long chains
            let motif: Vec<u8> = (0..1 + rng.below(40)).map(|_| rng.below(4) as u8).collect();
            let seq: Vec<u8> = (0..len)
                .map(|i| {
                    let r = rng.below(100);
                    if r < 2 {
                        14 // N
                    } else if r < 50 {
                        motif[i % motif.len()]
                    } else {
                        rng.below(4) as u8
                    }
                })
                .collect();
            let mut m = Vec::new();
            let mut cursor = 0usize;
            while cursor < len && rng.below(3) != 0 {
                let start = cursor + rng.below(400) as usize;
                let end = start + 1 + rng.below(200) as usize;
                if start >= len {
                    break;
                }
                m.push(MaskedInterval {
                    start,
                    end: end.min(len),
                });
                cursor = end;
            }
            offsets.push(offset);
            offset += len as i32 + 1;
            queries.push(seq);
            masks.push(m);
        }
        (queries, offsets, masks, word_length, lut_word_length)
    }

    #[test]
    fn partitioned_fill_builds_the_same_tables() {
        let mut rng = Lcg(0x5EED_1234);
        for case in 0..60 {
            let (queries, offsets, masks, word_length, lut_word_length) = random_case(&mut rng);
            let ascending_cells = case % 2 == 1;
            let approx = queries.iter().map(Vec::len).sum::<usize>();
            let reference = build_two_stage_lookup(
                &queries,
                &offsets,
                word_length,
                lut_word_length,
                &masks,
                None,
                0,
                approx,
                ascending_cells,
            );
            let mb = &reference.mb_lookup;
            for parts in 1..=8usize {
                let mut hashtable = vec![0u32; mb.hashtable.len()];
                let mut next_pos = vec![0u32; mb.next_pos.len()];
                let mut pv_array = vec![0 as PvArrayType; mb.pv_array.len()];
                let mut helper_array = vec![0u32; (mb.hashtable.len() >> 11).max(1)];
                let words = Words {
                    queries_blastna: &queries,
                    query_offsets: &offsets,
                    query_masks: &masks,
                    word_length,
                    lut_word_length,
                    kmer_mask: (1u64 << (2 * lut_word_length)) - 1,
                    ascending_cells,
                };
                fill(
                    parts,
                    &words,
                    Tables {
                        hashtable: &mut hashtable,
                        next_pos: &mut next_pos,
                        pv_array: &mut pv_array,
                        pv_array_bts: mb.pv_array_bts,
                        helper_array: &mut helper_array,
                    },
                    false,
                );
                let mut longest_chain = 2usize;
                for &value in &helper_array {
                    longest_chain = longest_chain.max(value as usize);
                }
                assert_eq!(hashtable, mb.hashtable, "case {case} parts {parts}");
                assert_eq!(next_pos, mb.next_pos, "case {case} parts {parts}");
                assert_eq!(pv_array, mb.pv_array, "case {case} parts {parts}");
                assert_eq!(longest_chain, mb.longest_chain, "case {case} parts {parts}");
            }
        }
    }

    // The same with real threads (a pool of 1..8 threads; the owner is inside the pool).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn partitioned_fill_on_a_pool_builds_the_same_tables() {
        let mut rng = Lcg(0xFACE_0FF5);
        for case in 0..12 {
            let (queries, offsets, masks, word_length, lut_word_length) = random_case(&mut rng);
            let ascending_cells = case % 2 == 0;
            let reference = build_two_stage_lookup(
                &queries,
                &offsets,
                word_length,
                lut_word_length,
                &masks,
                None,
                0,
                1_000_000,
                ascending_cells,
            );
            let mb = &reference.mb_lookup;
            for threads in 1..=8usize {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let (hashtable, next_pos, pv_array) = pool.install(|| {
                    let mut hashtable = vec![0u32; mb.hashtable.len()];
                    let mut next_pos = vec![0u32; mb.next_pos.len()];
                    let mut pv_array = vec![0 as PvArrayType; mb.pv_array.len()];
                    let mut helper_array = vec![0u32; (mb.hashtable.len() >> 11).max(1)];
                    let words = Words {
                        queries_blastna: &queries,
                        query_offsets: &offsets,
                        query_masks: &masks,
                        word_length,
                        lut_word_length,
                        kmer_mask: (1u64 << (2 * lut_word_length)) - 1,
                        ascending_cells,
                    };
                    fill(
                        threads,
                        &words,
                        Tables {
                            hashtable: &mut hashtable,
                            next_pos: &mut next_pos,
                            pv_array: &mut pv_array,
                            pv_array_bts: mb.pv_array_bts,
                            helper_array: &mut helper_array,
                        },
                        true,
                    );
                    (hashtable, next_pos, pv_array)
                });
                assert_eq!(hashtable, mb.hashtable, "case {case} threads {threads}");
                assert_eq!(next_pos, mb.next_pos, "case {case} threads {threads}");
                assert_eq!(pv_array, mb.pv_array, "case {case} threads {threads}");
            }
        }
    }
}
