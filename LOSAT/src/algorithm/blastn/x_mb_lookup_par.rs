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
//!
//! EXPERIMENT (LOSAT_X_MBLUTPAD; off unless set, not in the strict set): private counters and
//! presence-vector words. The blocks are dealt round robin, so one 64-byte line of
//! `helper_array` (16 counters, one per 2048-cell block when `pv_array_bts <= 11`, as for every
//! table of `lut_word_length` 12 or less) holds the counters of up to 16 different workers, and a
//! line of `pv_array` holds the presence-vector words of `16 >> (11 - pv_array_bts)` blocks (up to
//! 8 workers at `pv_array_bts` 10). Every word that finds its cell non-empty increments a counter
//! and every word that finds it empty sets a presence bit: each is a read and a write of a line
//! that other cores write at the same time (false sharing). Measured on the workstation (16
//! threads, megablast of 31,168 D. melanogaster cDNAs against 3 Mb, 37 tables), inside the fills:
//! 115 M demand reads found their line modified in another core (`ocr.demand_data_rd.l3_hit.
//! snoop_hitm`), against 0.07 M for the serial loop and 0.13 M when one thread ran every worker's
//! part one after another; with the private copies 0.8 M. The stores of `next_pos` still take
//! lines from other cores (each worker writes the offsets of its own words, which interleave
//! with the other workers'; 32 M such stores in that run), but stores wait in the store buffer;
//! writing `next_pos` afterwards by spans of offsets removed them and cost more than it saved
//! (lists of 8 bytes a word and their growth), so it is not done.
//! With the switch each worker counts into a private copy of `helper_array` and sets its bits in a
//! private copy of `pv_array`, each a whole number of 64-byte lines starting on a line boundary
//! (`Lines`), and at the end of its part adds its own blocks' counters into `helper_array`
//! (`saturating_add`) and ORs its own blocks' presence-vector words into `pv_array`.
//!
//! Why the counters and presence bits are the loop's. Counter `h` covers the cells
//! `h * 2048 .. (h + 1) * 2048` and presence-vector word `w` the cells
//! `w << pv_array_bts .. (w + 1) << pv_array_bts`; both ranges lie inside one block, so only the
//! block's owner sees words of those cells, and it sees all of them, in the loop's order (above).
//! Its test `hashtable[cell] == 0` therefore takes the loop's value for every word, so the private
//! copy, starting at 0, receives exactly the loop's increments of counter `h` (each a saturating
//! `+1`, so the copy ends at `min(n, u32::MAX)` for `n` increments) and exactly the loop's bits of
//! word `w`. The owner then writes `helper_array[h] = helper_array[h].saturating_add(copy[h])`,
//! which is `min(s + n, u32::MAX)` for a start value `s`, the loop's `n` saturating increments of
//! `s` (the start is 0 in every build), and `pv_array[w] |= copy[w]`, the loop's bits ORed into the
//! start value (OR does not depend on order). No other worker writes counter `h` or word `w`
//! during the fill (they write their own copies and their own blocks' entries), and nothing reads
//! `helper_array` or `pv_array` before the fill returns (`longest_chain` and the presence tests of
//! the scan come after). The re-threading of `ascending_cells` changes no counter and no bit. The
//! update of `hashtable` and `next_pos` is the original statement for statement. So every element
//! of the four arrays, and `longest_chain`, is the loop's; the switch changes the running time
//! only. Memory: two private copies per worker, `ceil(hashsize / 2048 / 16)` and
//! `ceil(pv_size / 16)` lines of 64 bytes (at most 32 KiB + 128 KiB per worker for a
//! `lut_word_length` 12 table, 8 KiB + 128 KiB at 11), in one zeroed allocation that the calling
//! thread makes before the fill and frees after it (`Copies`; with each worker allocating its own
//! copies, glibc's per-thread arenas kept 8-19 MiB more at the end of a search).

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
pub(super) const K_COMPRESSION_FACTOR_SHIFT: usize = 11;

/// No NCBI counterpart: reads the switches once; they do not change any value NCBI computes.
/// (LOSAT_X_PAIRPAR, LOSAT_X_PAIRPARSHADOW, LOSAT_X_MBLUTPAD; the last is off unless set to any
/// value and not in the strict set)
fn switches() -> (bool, bool, bool) {
    static SWITCHES: OnceLock<(bool, bool, bool)> = OnceLock::new();
    *SWITCHES.get_or_init(|| {
        (
            crate::utils::xswitch::on("LOSAT_X_PAIRPAR"),
            std::env::var_os("LOSAT_X_PAIRPARSHADOW").is_some(),
            std::env::var_os("LOSAT_X_MBLUTPAD").is_some(),
        )
    })
}

/// No NCBI counterpart: LOSAT_TIMING is set (read once).
fn timing_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_TIMING").is_some())
}

// No NCBI counterpart (LOSAT_TIMING only): the `[TIMING] mb_lookup_*` counters of every megablast
// table build of the process (not the shadow's second build).
static BUILD_NANOS: AtomicU64 = AtomicU64::new(0);
static ALLOC_NANOS: AtomicU64 = AtomicU64::new(0);
static LOOP_NANOS: AtomicU64 = AtomicU64::new(0);
static POOL_FILL_NANOS: AtomicU64 = AtomicU64::new(0);
static CHAIN_NANOS: AtomicU64 = AtomicU64::new(0);
static BUILDS: AtomicU64 = AtomicU64::new(0);
static BUILDS_ON_POOL: AtomicU64 = AtomicU64::new(0);
static BUILDS_PADDED: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart: the marks of one `build_mb_lookup` for the `[TIMING]` lines (only with
/// LOSAT_TIMING): the start of the build, of the fill and of the longest chain.
pub(super) struct BuildClock {
    start: Option<std::time::Instant>,
    fill: Option<std::time::Instant>,
    chain: Option<std::time::Instant>,
}

impl BuildClock {
    /// The start of a build (before the allocation of the tables).
    pub(super) fn start() -> Self {
        Self {
            start: timing_on().then(std::time::Instant::now),
            fill: None,
            chain: None,
        }
    }

    /// The start of the fill (the original loop or the pool fill).
    pub(super) fn fill(&mut self) {
        if self.start.is_some() {
            self.fill = Some(std::time::Instant::now());
        }
    }

    /// The start of the longest chain (the end of the fill and of the re-threading of
    /// `ascending_cells`).
    pub(super) fn chain(&mut self) {
        if self.start.is_some() {
            self.chain = Some(std::time::Instant::now());
        }
    }

    /// The end of the build: adds its intervals to the `[TIMING]` counters.
    pub(super) fn record(self, on_pool: bool) {
        let (Some(start), Some(fill), Some(chain)) = (self.start, self.fill, self.chain) else {
            return;
        };
        let end = std::time::Instant::now();
        let ns = |from: std::time::Instant, to: std::time::Instant| {
            to.saturating_duration_since(from).as_nanos() as u64
        };
        let r = Ordering::Relaxed;
        BUILD_NANOS.fetch_add(ns(start, end), r);
        ALLOC_NANOS.fetch_add(ns(start, fill), r);
        if on_pool {
            &POOL_FILL_NANOS
        } else {
            &LOOP_NANOS
        }
        .fetch_add(ns(fill, chain), r);
        CHAIN_NANOS.fetch_add(ns(chain, end), r);
        BUILDS.fetch_add(1, r);
        BUILDS_ON_POOL.fetch_add(u64::from(on_pool), r);
        BUILDS_PADDED.fetch_add(u64::from(on_pool && switches().2), r);
    }
}

/// No NCBI counterpart: shadow counters of LOSAT_X_MBLUTPAD (tables whose second build was the
/// padded fill, and their elements).
static SHADOW_PADDED_TABLES: AtomicU64 = AtomicU64::new(0);
static SHADOW_PADDED_CELLS: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart: printed by `main` at exit. With LOSAT_X_MBLUTPAD and LOSAT_X_PAIRPARSHADOW,
/// `[X_SHADOW] MBLUTPAD tables=<compared> cells=<elements compared> diffs=0 identical` (a
/// difference panics in `shadow_check`, so the line only ever says 0). With LOSAT_TIMING, when at
/// least one megablast table was built, the `[TIMING] mb_lookup_build` and `[TIMING]
/// mb_lookup_split` lines. They print counters and change no value NCBI computes.
pub(super) fn print_summary() {
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
    let (_, shadow, pad) = switches();
    if shadow && pad {
        eprintln!(
            "[X_SHADOW] MBLUTPAD tables={} cells={} diffs=0 identical",
            load(&SHADOW_PADDED_TABLES),
            load(&SHADOW_PADDED_CELLS)
        );
    }
    if !timing_on() || load(&BUILDS) == 0 {
        return;
    }
    let s = |c: &AtomicU64| load(c) as f64 / 1e9;
    eprintln!(
        "[TIMING] mb_lookup_build: {:.3}s (builds={} parallel={} padded={})",
        s(&BUILD_NANOS),
        load(&BUILDS),
        load(&BUILDS_ON_POOL),
        load(&BUILDS_PADDED)
    );
    eprintln!(
        "[TIMING] mb_lookup_split: alloc={:.4}s loop={:.4}s pool_fill={:.4}s chain={:.4}s",
        s(&ALLOC_NANOS),
        s(&LOOP_NANOS),
        s(&POOL_FILL_NANOS),
        s(&CHAIN_NANOS)
    );
}

/// No NCBI counterpart: shadow counters for the summary line (tables compared).
pub(crate) static SHADOW_TABLES: AtomicU64 = AtomicU64::new(0);
/// No NCBI counterpart: shadow counters for the summary line (cells compared).
pub(crate) static SHADOW_CELLS: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart: the number of threads of the search pool this thread runs in, when
/// there is one with two threads or more (scheduling only). Also the dispatch test of
/// LOSAT_X_DISCLUTPAR (`x_disc_lookup_par.rs`).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) fn pool_threads() -> Option<usize> {
    rayon::current_thread_index()?;
    let threads = rayon::current_num_threads();
    (threads > 1).then_some(threads)
}

#[cfg(any(
    not(feature = "parallel"),
    all(target_arch = "wasm32", not(feature = "wasm-threads"))
))]
pub(super) fn pool_threads() -> Option<usize> {
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
    let (on, shadow, _) = switches();
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
    pub(super) queries_blastna: &'a [&'a [u8]],
    pub(super) query_offsets: &'a [i32],
    pub(super) query_masks: &'a [Vec<MaskedInterval>],
    pub(super) word_length: usize,
    pub(super) lut_word_length: usize,
    pub(super) kmer_mask: u64,
    pub(super) ascending_cells: bool,
}

// No NCBI counterpart: a slice shared by the workers, each of which reads and writes only the
// elements it owns (see the module documentation); the pointer outlives the scope. Also used by
// LOSAT_X_DISCLUTPAR (`x_disc_lookup_par.rs`), whose workers own elements the same way.
#[derive(Clone, Copy)]
pub(super) struct Shared<T> {
    pub(super) ptr: *mut T,
    pub(super) len: usize,
}

// SAFETY: the workers access disjoint elements only (ownership by cell block, and by the
// word's own query offset for `next_pos`), and the slices outlive the scope that uses them.
unsafe impl<T: Send> Send for Shared<T> {}
// SAFETY: as above.
unsafe impl<T: Send> Sync for Shared<T> {}

impl<T: Copy> Shared<T> {
    pub(super) fn new(slice: &mut [T]) -> Self {
        Self {
            ptr: slice.as_mut_ptr(),
            len: slice.len(),
        }
    }

    /// SAFETY: `i < len`, and no other worker accesses element `i`.
    #[inline(always)]
    pub(super) unsafe fn get(self, i: usize) -> T {
        debug_assert!(i < self.len);
        unsafe { *self.ptr.add(i) }
    }

    /// SAFETY: `i < len`, and no other worker accesses element `i`.
    #[inline(always)]
    pub(super) unsafe fn set(self, i: usize, value: T) {
        debug_assert!(i < self.len);
        unsafe { *self.ptr.add(i) = value }
    }
}

/// No NCBI counterpart (LOSAT_X_MBLUTPAD): sixteen 32-bit words on a 64-byte line of their own.
#[derive(Clone, Copy)]
#[repr(C, align(64))]
struct Line([u32; 16]);

const _: () = assert!(std::mem::size_of::<Line>() == 64 && std::mem::align_of::<Line>() == 64);
const _: () = assert!(std::mem::size_of::<PvArrayType>() == std::mem::size_of::<u32>());

/// No NCBI counterpart (LOSAT_X_MBLUTPAD): `n` zeroed lines in one allocation (all-zero bytes are
/// a valid `Line`); a fresh allocation of this size comes from zero pages, so a line's page is
/// first touched by the worker that writes it.
fn zeroed_lines(n: usize) -> Vec<Line> {
    if n == 0 {
        return Vec::new();
    }
    let layout = std::alloc::Layout::array::<Line>(n).expect("private copies fit in memory");
    // SAFETY: the layout has a non-zero size; zeroed memory is a valid `[Line]`; the vector owns
    // the allocation with this layout (`Layout::array::<Line>(n)`, as `Vec<Line>` frees it).
    unsafe {
        let ptr = std::alloc::alloc_zeroed(layout) as *mut Line;
        if ptr.is_null() {
            std::alloc::handle_alloc_error(layout);
        }
        Vec::from_raw_parts(ptr, n, n)
    }
}

/// No NCBI counterpart (LOSAT_X_MBLUTPAD): the private copies of all workers in one zeroed
/// allocation made by the calling thread (one allocation per fill, so the workers' own malloc
/// arenas take nothing), `stride` lines per worker: first the copy of `helper_array`
/// (`counter_lines`), then that of `pv_array`. Each worker's lines start on a line boundary and
/// are a whole number of lines, so no two workers' copies share a line. Empty when the switch is
/// off.
struct Copies {
    lines: Vec<Line>,
    counter_lines: usize,
    stride: usize,
}

impl Copies {
    fn new(pad: bool, parts: usize, helper_len: usize, pv_len: usize) -> Self {
        let (counter_lines, pv_lines) = if pad {
            (helper_len.div_ceil(16), pv_len.div_ceil(16))
        } else {
            (0, 0)
        };
        let stride = counter_lines + pv_lines;
        Self {
            lines: zeroed_lines(stride * parts),
            counter_lines,
            stride,
        }
    }

    /// The workers' copies, in part order (each `stride` lines; empty slices when the switch is
    /// off).
    fn parts(&mut self, parts: usize) -> impl Iterator<Item = Private<'_>> {
        let (stride, counter_lines) = (self.stride, self.counter_lines);
        let mut rest: &mut [Line] = &mut self.lines;
        (0..parts).map(move |_| {
            let (mine, tail) = std::mem::take(&mut rest).split_at_mut(stride);
            rest = tail;
            let (counters, pv) = mine.split_at_mut(counter_lines);
            Private {
                counters: Lines { lines: counters },
                pv: Lines { lines: pv },
            }
        })
    }
}

/// No NCBI counterpart (LOSAT_X_MBLUTPAD): a worker's private copy of `helper_array` or
/// `pv_array` (element `i` in line `i / 16`).
struct Lines<'a> {
    lines: &'a mut [Line],
}

impl Lines<'_> {
    #[inline(always)]
    fn get(&self, i: usize) -> u32 {
        self.lines[i >> 4].0[i & 15]
    }

    #[inline(always)]
    fn slot(&mut self, i: usize) -> &mut u32 {
        &mut self.lines[i >> 4].0[i & 15]
    }
}

/// No NCBI counterpart (LOSAT_X_MBLUTPAD): the private copies of one worker (empty when the
/// switch is off).
struct Private<'a> {
    counters: Lines<'a>,
    pv: Lines<'a>,
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
    // the presence bit is set in the presence-vector word this worker owns. With `PAD`
    // (LOSAT_X_MBLUTPAD) the counter and the presence bit go to the worker's private copies
    // (module documentation); without it this is the original update.
    /// SAFETY: the worker owns the block of `bucket` and the word at `q_off_1`.
    #[inline(always)]
    unsafe fn apply<const PAD: bool>(self, bucket: usize, q_off_1: u32, x: &mut Private<'_>) {
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
                    if PAD {
                        *x.pv.slot(array_idx) |= bit;
                    } else {
                        self.pv_array
                            .set(array_idx, self.pv_array.get(array_idx) | bit);
                    }
                }
            } else {
                let h = bucket >> K_COMPRESSION_FACTOR_SHIFT;
                if PAD {
                    let c = x.counters.slot(h);
                    *c = c.saturating_add(1);
                } else {
                    self.helper_array
                        .set(h, self.helper_array.get(h).saturating_add(1));
                }
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
/// the update only to the words whose cell block `owner` gives to `part`. With `PAD`
/// (LOSAT_X_MBLUTPAD) the counters and presence bits go to private copies, which the worker adds
/// into `helper_array` and `pv_array` for its own blocks at the end (module documentation).
fn fill_part<const PAD: bool>(
    part: u8,
    owner: &[u8],
    block_shift: usize,
    words: &Words<'_>,
    t: SharedTables,
    mut x: Private<'_>,
) {
    let lut_word_length = words.lut_word_length;
    let mut own_words: Vec<(u32, u32)> = Vec::new();
    #[cfg(target_arch = "x86_64")]
    const X_DELAY: usize = 16;
    #[cfg(target_arch = "x86_64")]
    let mut x_ring = [(0u32, 0u32); X_DELAY];
    #[cfg(target_arch = "x86_64")]
    let mut x_ring_len = 0usize;
    for (q_idx, seq_blastna) in words.queries_blastna.iter().enumerate() {
        let seq = seq_blastna;
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
                        unsafe { t.apply::<PAD>(b as usize, q, &mut x) };
                    }
                    x_ring[slot] = (bucket as u32, q_off_1);
                    x_ring_len += 1;
                }
                #[cfg(not(target_arch = "x86_64"))]
                {
                    // SAFETY: as above.
                    unsafe { t.apply::<PAD>(bucket, q_off_1, &mut x) };
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
            unsafe { t.apply::<PAD>(b as usize, q, &mut x) };
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
    if PAD {
        flush_private(part, owner, block_shift, t, &x);
    }
}

/// No NCBI counterpart (LOSAT_X_MBLUTPAD): adds the counters and ORs the presence-vector words
/// of `part`'s own blocks from its private copies into `helper_array` and `pv_array` (module
/// documentation: these elements are written by no other worker, and their copies hold exactly
/// the loop's increments and bits).
fn flush_private(part: u8, owner: &[u8], block_shift: usize, t: SharedTables, x: &Private<'_>) {
    for (block, &o) in owner.iter().enumerate() {
        if o != part {
            continue;
        }
        let (lo, hi) = (block << block_shift, (block + 1) << block_shift);
        let counters = (lo >> K_COMPRESSION_FACTOR_SHIFT).min(t.helper_array.len)
            ..(hi >> K_COMPRESSION_FACTOR_SHIFT).min(t.helper_array.len);
        for h in counters {
            let c = x.counters.get(h);
            if c != 0 {
                // SAFETY: counter `h` lies in a block of this part.
                unsafe {
                    t.helper_array
                        .set(h, t.helper_array.get(h).saturating_add(c))
                };
            }
        }
        let words =
            (lo >> t.pv_array_bts).min(t.pv_array.len)..(hi >> t.pv_array_bts).min(t.pv_array.len);
        for w in words {
            let bits = x.pv.get(w);
            if bits != 0 {
                // SAFETY: presence-vector word `w` lies in a block of this part.
                unsafe { t.pv_array.set(w, t.pv_array.get(w) | bits) };
            }
        }
    }
}

/// No NCBI counterpart: the cell blocks and their owners (scheduling only). Blocks of
/// `2^block_shift` cells, dealt round robin to `parts` workers. Also the plan of
/// LOSAT_X_DISCLUTPAR (`x_disc_lookup_par.rs`), for both templates' tables.
pub(super) fn plan(table_size: usize, pv_array_bts: usize, parts: usize) -> (usize, Vec<u8>) {
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
    fill_mode(switches().2, parts, words, tables, in_pool);
}

/// No NCBI counterpart: `fill` with LOSAT_X_MBLUTPAD given as `pad` (the tests take both).
fn fill_mode(pad: bool, parts: usize, words: &Words<'_>, tables: Tables<'_>, in_pool: bool) {
    // Dispatch point of LOSAT_X_MBLUTPAD: private counters and presence-vector words (module
    // documentation); off, the original fill.
    if pad {
        fill_with::<true>(parts, words, tables, in_pool);
    } else {
        fill_with::<false>(parts, words, tables, in_pool);
    }
}

/// No NCBI counterpart: `fill` with the private copies of LOSAT_X_MBLUTPAD (`PAD`) or without.
fn fill_with<const PAD: bool>(parts: usize, words: &Words<'_>, tables: Tables<'_>, in_pool: bool) {
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
    let mut copies = Copies::new(PAD, parts, shared.helper_array.len, shared.pv_array.len);
    let mut copies = copies.parts(parts);
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if in_pool {
        let first = copies.next().expect("one part at least");
        rayon::scope(|scope| {
            for (part, x) in (1..parts).zip(copies) {
                scope.spawn(move |_| {
                    fill_part::<PAD>(part as u8, owner, block_shift, words, shared, x)
                });
            }
            fill_part::<PAD>(0, owner, block_shift, words, shared, first);
        });
        return;
    }
    let _ = in_pool;
    for (part, x) in (0..parts).zip(copies) {
        fill_part::<PAD>(part as u8, owner, block_shift, words, shared, x);
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
/// With LOSAT_X_MBLUTPAD the second build is the padded fill.
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
    shadow_check_mode(
        switches().2,
        words,
        pv_array_bts,
        hashtable,
        next_pos,
        pv_array,
        helper_array,
        longest_chain,
    );
}

/// No NCBI counterpart: `shadow_check` with LOSAT_X_MBLUTPAD given as `pad`.
#[allow(clippy::too_many_arguments)]
fn shadow_check_mode(
    pad: bool,
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
    fill_mode(
        pad,
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
                "LOSAT_X_PAIRPARSHADOW: megablast lookup {name}[{i}] differs: reference {} partitioned {} (parts {parts}, lut_word_length {}, ascending_cells {}, padded {pad})",
                reference[i], partitioned[i], words.lut_word_length, words.ascending_cells
            );
        }
    }
    assert_eq!(
        longest_chain, x_longest_chain,
        "LOSAT_X_PAIRPARSHADOW: megablast lookup longest_chain differs (parts {parts}, padded {pad})"
    );
    SHADOW_TABLES.fetch_add(1, Ordering::Relaxed);
    let cells = (hashtable.len() + next_pos.len() + pv_array.len() + helper_array.len()) as u64;
    SHADOW_CELLS.fetch_add(cells, Ordering::Relaxed);
    if pad {
        SHADOW_PADDED_TABLES.fetch_add(1, Ordering::Relaxed);
        SHADOW_PADDED_CELLS.fetch_add(cells, Ordering::Relaxed);
    }
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
            let query_refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
            let ascending_cells = case % 2 == 1;
            let approx = queries.iter().map(Vec::len).sum::<usize>();
            let reference = build_two_stage_lookup(
                &query_refs,
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
                    queries_blastna: &query_refs,
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
            let query_refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
            let ascending_cells = case % 2 == 0;
            let reference = build_two_stage_lookup(
                &query_refs,
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
                        queries_blastna: &query_refs,
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

    // LOSAT_X_MBLUTPAD: the tables of one case, built by the unpadded fill in one part (the loop's
    // statements in reading order; the tests above check it against `build_two_stage_lookup`),
    // from start values `helper0` / `pv0` of the counters and presence-vector words.
    struct Built {
        hashtable: Vec<u32>,
        next_pos: Vec<u32>,
        pv_array: Vec<PvArrayType>,
        helper_array: Vec<u32>,
    }

    fn build_with(
        pad: bool,
        parts: usize,
        in_pool: bool,
        words: &Words<'_>,
        shape: (usize, usize, usize, usize),
        helper0: &[u32],
        pv0: &[PvArrayType],
    ) -> Built {
        let (cells, next_len, pv_len, pv_array_bts) = shape;
        let mut b = Built {
            hashtable: vec![0u32; cells],
            next_pos: vec![0u32; next_len],
            pv_array: pv0.to_vec(),
            helper_array: helper0.to_vec(),
        };
        assert_eq!(b.pv_array.len(), pv_len);
        fill_mode(
            pad,
            parts,
            words,
            Tables {
                hashtable: &mut b.hashtable,
                next_pos: &mut b.next_pos,
                pv_array: &mut b.pv_array,
                pv_array_bts,
                helper_array: &mut b.helper_array,
            },
            in_pool,
        );
        b
    }

    fn assert_same(a: &Built, b: &Built, what: &str) {
        assert_eq!(a.hashtable, b.hashtable, "{what}: hashtable");
        assert_eq!(a.next_pos, b.next_pos, "{what}: next_pos");
        assert_eq!(a.pv_array, b.pv_array, "{what}: pv_array");
        assert_eq!(a.helper_array, b.helper_array, "{what}: helper_array");
    }

    // A query with long unmasked repeats (many words per cell, so counters grow) and random letters.
    fn repeat_case(rng: &mut Lcg, len: usize) -> Vec<u8> {
        (0..len)
            .map(|i| {
                if (i / 700) % 3 == 0 {
                    [0u8, 0, 1, 0, 2][i % 5]
                } else if rng.below(200) == 0 {
                    14
                } else {
                    rng.below(4) as u8
                }
            })
            .collect()
    }

    // LOSAT_X_MBLUTPAD on real pools of 2, 3, 5 and 16 threads (and off the pool with the same
    // parts): every element of the four arrays and the longest chain equal the unpadded fill's
    // and `build_two_stage_lookup`'s, for random cases (ambiguity codes, masks, several queries,
    // `lut_word_length` 4..12, both `ascending_cells`) and for `lut_word_length` 11 and 12 tables
    // whose presence-vector words share lines between blocks (`pv_array_bts` 8 to 10).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn padded_fill_gives_the_loops_tables_on_pools() {
        let mut rng = Lcg(0x0B1D_5EED);
        let mut cases = Vec::new();
        for case in 0..8 {
            let (queries, offsets, masks, word_length, lut_word_length) = random_case(&mut rng);
            cases.push((
                queries,
                offsets,
                masks,
                word_length,
                lut_word_length,
                case % 2 == 1,
            ));
        }
        for (lut, asc) in [(11usize, false), (12, false), (12, true)] {
            let q = repeat_case(&mut rng, 40_000);
            cases.push((vec![q], vec![0], vec![Vec::new()], lut + 5, lut, asc));
        }
        let mut seen_shared_pv_lines = 0;
        for (case, (queries, offsets, masks, word_length, lut_word_length, asc)) in
            cases.iter().enumerate()
        {
            let query_refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
            let reference = build_two_stage_lookup(
                &query_refs,
                offsets,
                *word_length,
                *lut_word_length,
                masks,
                None,
                0,
                1_000_000,
                *asc,
            );
            let mb = &reference.mb_lookup;
            let shape = (
                mb.hashtable.len(),
                mb.next_pos.len(),
                mb.pv_array.len(),
                mb.pv_array_bts,
            );
            if mb.hashtable.len() >= 1 << 22 && mb.pv_array_bts > 7 {
                seen_shared_pv_lines += 1;
            }
            let helper0 = vec![0u32; (mb.hashtable.len() >> 11).max(1)];
            let pv0 = vec![0 as PvArrayType; mb.pv_array.len()];
            let words = Words {
                queries_blastna: &query_refs,
                query_offsets: offsets,
                query_masks: masks,
                word_length: *word_length,
                lut_word_length: *lut_word_length,
                kmer_mask: (1u64 << (2 * lut_word_length)) - 1,
                ascending_cells: *asc,
            };
            let loop_built = build_with(false, 1, false, &words, shape, &helper0, &pv0);
            assert_eq!(loop_built.hashtable, mb.hashtable, "case {case}");
            assert_eq!(loop_built.next_pos, mb.next_pos, "case {case}");
            assert_eq!(loop_built.pv_array, mb.pv_array, "case {case}");
            let chain = loop_built
                .helper_array
                .iter()
                .fold(2usize, |m, &v| m.max(v as usize));
            assert_eq!(chain, mb.longest_chain, "case {case}");
            for threads in [2usize, 3, 5, 16] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let padded =
                    pool.install(|| build_with(true, threads, true, &words, shape, &helper0, &pv0));
                assert_same(&loop_built, &padded, &format!("case {case} pool {threads}"));
                let unpadded = pool
                    .install(|| build_with(false, threads, true, &words, shape, &helper0, &pv0));
                assert_same(
                    &loop_built,
                    &unpadded,
                    &format!("case {case} pool {threads} off"),
                );
                let serial = build_with(true, threads, false, &words, shape, &helper0, &pv0);
                assert_same(
                    &loop_built,
                    &serial,
                    &format!("case {case} parts {threads}"),
                );
            }
        }
        assert!(
            seen_shared_pv_lines >= 3,
            "the lut 11 / 12 cases ({seen_shared_pv_lines})"
        );
    }

    // LOSAT_X_MBLUTPAD from start values other than 0 (production starts at 0): the flush's
    // `saturating_add` gives the loop's saturating increments (counters starting at u32::MAX,
    // u32::MAX - 0..6, 0 or a hash value), and the OR keeps start bits; some counters saturate,
    // and zero-start counters with no repeat stay 0 (checked for the lut 11 and 12 tables).
    #[test]
    fn padded_flush_matches_the_loop_from_any_start() {
        let mut rng = Lcg(0x5A7_0F1D);
        for (lut, asc) in [(8usize, false), (11, true), (12, false)] {
            let q = repeat_case(&mut rng, 30_000);
            let queries = vec![q];
            let query_refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
            let offsets = vec![0i32];
            let masks = vec![Vec::new()];
            let reference = build_two_stage_lookup(
                &query_refs,
                &offsets,
                lut + 4,
                lut,
                &masks,
                None,
                0,
                30_000,
                asc,
            );
            let mb = &reference.mb_lookup;
            let shape = (
                mb.hashtable.len(),
                mb.next_pos.len(),
                mb.pv_array.len(),
                mb.pv_array_bts,
            );
            let helper_len = (mb.hashtable.len() >> 11).max(1);
            let helper0: Vec<u32> = (0..helper_len)
                .map(|h| match h % 4 {
                    0 => u32::MAX,
                    1 => u32::MAX - (h as u32 % 7),
                    2 => 0,
                    _ => (h as u32).wrapping_mul(2_654_435_761),
                })
                .collect();
            let pv0: Vec<PvArrayType> = (0..mb.pv_array.len())
                .map(|w| {
                    if w % 3 == 0 {
                        (w as u32).wrapping_mul(40_503)
                    } else {
                        0
                    }
                })
                .collect();
            let words = Words {
                queries_blastna: &query_refs,
                query_offsets: &offsets,
                query_masks: &masks,
                word_length: lut + 4,
                lut_word_length: lut,
                kmer_mask: (1u64 << (2 * lut)) - 1,
                ascending_cells: asc,
            };
            let loop_built = build_with(false, 1, false, &words, shape, &helper0, &pv0);
            let zero_built = build_with(
                false,
                1,
                false,
                &words,
                shape,
                &vec![0u32; helper_len],
                &vec![0 as PvArrayType; mb.pv_array.len()],
            );
            let mut saturated = 0;
            let mut zero_stays = 0;
            for h in 0..helper_len {
                let n = zero_built.helper_array[h];
                let (s, v) = (helper0[h], loop_built.helper_array[h]);
                assert_eq!(v, s.saturating_add(n), "lut {lut} counter {h}");
                saturated += usize::from(s != u32::MAX && s.checked_add(n).is_none());
                zero_stays += usize::from(s == 0 && n == 0);
            }
            if lut >= 11 {
                assert!(
                    saturated > 0 && zero_stays > 0,
                    "lut {lut}: {saturated} {zero_stays}"
                );
            }
            for parts in [2usize, 3, 5, 16] {
                let padded = build_with(true, parts, false, &words, shape, &helper0, &pv0);
                assert_same(&loop_built, &padded, &format!("lut {lut} parts {parts}"));
            }
            #[cfg(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            ))]
            for threads in [2usize, 3, 5, 16] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let padded =
                    pool.install(|| build_with(true, threads, true, &words, shape, &helper0, &pv0));
                assert_same(&loop_built, &padded, &format!("lut {lut} pool {threads}"));
            }
        }
    }

    // LOSAT_X_MBLUTPAD: each worker's copies are a whole number of 64-byte lines starting on a line
    // boundary, the workers' copies are disjoint and in part order, zeroed, and the switch off
    // allocates nothing.
    #[test]
    fn private_copies_have_lines_of_their_own() {
        assert_eq!(std::mem::size_of::<Line>(), 64);
        assert_eq!(std::mem::align_of::<Line>(), 64);
        for (parts, helper_len, pv_len) in [
            (1usize, 1usize, 1usize),
            (2, 2048, 16384),
            (3, 8192, 32768),
            (5, 15, 17),
            (16, 2048, 32768),
            (16, 8192, 16384),
        ] {
            let mut copies = Copies::new(true, parts, helper_len, pv_len);
            assert_eq!(copies.lines.as_ptr() as usize % 64, 0);
            let stride = helper_len.div_ceil(16) + pv_len.div_ceil(16);
            let base = copies.lines.as_ptr() as usize;
            let mut seen = 0;
            for (p, mut x) in copies.parts(parts).enumerate() {
                let c0 = x.counters.lines.as_ptr() as usize;
                let v0 = x.pv.lines.as_ptr() as usize;
                assert_eq!(c0, base + p * stride * 64, "parts {parts} part {p}");
                assert_eq!(
                    v0,
                    c0 + helper_len.div_ceil(16) * 64,
                    "parts {parts} part {p}"
                );
                assert_eq!(x.counters.lines.len(), helper_len.div_ceil(16));
                assert_eq!(x.pv.lines.len(), pv_len.div_ceil(16));
                assert!((0..helper_len).all(|i| x.counters.get(i) == 0));
                assert!((0..pv_len).all(|i| x.pv.get(i) == 0));
                *x.counters.slot(helper_len - 1) = 1;
                *x.pv.slot(pv_len - 1) = 1;
                seen += 1;
            }
            assert_eq!(seen, parts);
            assert_eq!(copies.lines.len(), stride * parts);
            assert_eq!(
                copies
                    .lines
                    .iter()
                    .map(|l| l.0.iter().sum::<u32>())
                    .sum::<u32>(),
                2 * parts as u32,
                "each part wrote its own two elements"
            );
        }
        let mut off = Copies::new(false, 16, 8192, 32768);
        assert_eq!(off.lines.capacity(), 0);
        assert_eq!(
            off.parts(16)
                .filter(|x| x.counters.lines.is_empty() && x.pv.lines.is_empty())
                .count(),
            16
        );
    }

    // The shadow (`LOSAT_X_PAIRPARSHADOW`) with the padded second build: equal tables pass, and a
    // difference in any one of the four arrays or the longest chain stops it.
    #[test]
    fn shadow_with_padding_compares_every_array() {
        let mut rng = Lcg(0xC0FF_EE11);
        let q = repeat_case(&mut rng, 20_000);
        let queries = vec![q];
        let query_refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        let offsets = vec![0i32];
        let masks = vec![Vec::new()];
        let lut = 11usize;
        let reference = build_two_stage_lookup(
            &query_refs,
            &offsets,
            16,
            lut,
            &masks,
            None,
            0,
            1_000_000,
            false,
        );
        let mb = &reference.mb_lookup;
        let shape = (
            mb.hashtable.len(),
            mb.next_pos.len(),
            mb.pv_array.len(),
            mb.pv_array_bts,
        );
        let words = Words {
            queries_blastna: &query_refs,
            query_offsets: &offsets,
            query_masks: &masks,
            word_length: 16,
            lut_word_length: lut,
            kmer_mask: (1u64 << (2 * lut)) - 1,
            ascending_cells: false,
        };
        let helper_len = (mb.hashtable.len() >> 11).max(1);
        let b = build_with(
            false,
            1,
            false,
            &words,
            shape,
            &vec![0u32; helper_len],
            &vec![0 as PvArrayType; mb.pv_array.len()],
        );
        let chain = b
            .helper_array
            .iter()
            .fold(2usize, |m, &v| m.max(v as usize));
        let check = |b: &Built, chain: usize| {
            shadow_check_mode(
                true,
                &words,
                mb.pv_array_bts,
                &b.hashtable,
                &b.next_pos,
                &b.pv_array,
                &b.helper_array,
                chain,
            )
        };
        check(&b, chain);
        let h = b.helper_array.iter().position(|&v| v > 0).unwrap();
        let w = b.pv_array.iter().position(|&v| v != 0).unwrap();
        let cell = b.hashtable.iter().position(|&v| v != 0).unwrap();
        let q = b.next_pos.iter().position(|&v| v != 0).unwrap();
        for which in 0..5 {
            let mut c = Built {
                hashtable: b.hashtable.clone(),
                next_pos: b.next_pos.clone(),
                pv_array: b.pv_array.clone(),
                helper_array: b.helper_array.clone(),
            };
            let mut c_chain = chain;
            match which {
                0 => c.hashtable[cell] += 1,
                1 => c.next_pos[q] += 1,
                2 => c.pv_array[w] ^= 1,
                3 => c.helper_array[h] -= 1,
                _ => c_chain += 1,
            }
            let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| check(&c, c_chain)));
            assert!(r.is_err(), "difference {which} not caught");
        }
    }
}
