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
//! Private counters and presence-vector words, which the read-once fill of LOSAT_X_MBLUTONCE
//! (`x_mb_lookup_once.rs`) counts into; the fill of this module updates the shared arrays itself.
//! The blocks are dealt round robin, so one 64-byte line of
//! `helper_array` (16 counters, one per 2048-cell block when `pv_array_bts <= 11`, as for every
//! table of `lut_word_length` 12 or less) holds the counters of up to 16 different workers, and a
//! line of `pv_array` holds the presence-vector words of `16 >> (11 - pv_array_bts)` blocks when
//! `pv_array_bts` is 7 to 10 (the `lut_word_length` 11 and 12 tables; up to 8 workers at 10). At
//! smaller tables (`pv_array_bts` 5 or 6, `lut_word_length` 8 to 10) a block's words span several
//! lines and only the lines at block boundaries are shared. Both arrays are only 16-byte aligned
//! (`vec!`), so a block's words may also straddle a line. Every word that finds its cell
//! non-empty increments a counter and every word that finds it empty sets a presence bit: each is
//! a read and a write of a line that other cores write at the same time (false sharing). Measured
//! on the workstation (16 threads, megablast of 31,168 D. melanogaster cDNAs against 3 Mb, 37
//! tables), inside the fills:
//! 115 M demand reads found their line modified in another core (`ocr.demand_data_rd.l3_hit.
//! snoop_hitm`), against 0.07 M for the serial loop and 0.13 M when one thread ran every worker's
//! part one after another; with the private copies 0.8 M. The stores of `next_pos` still take
//! lines from other cores (each worker writes the offsets of its own words, which interleave
//! with the other workers'; 32 M such stores in that run), but stores wait in the store buffer;
//! writing `next_pos` afterwards by spans of offsets removed them and cost more than it saved
//! (lists of 8 bytes a word and their growth), so it is not done.
//! Each worker of the read-once fill counts into a private copy of `helper_array` and sets its bits in a
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
//! of the four arrays, and `longest_chain`, is the loop's; the copies change the running time
//! only. Memory: two private copies per worker, `ceil(hashsize / 2048 / 16)` and
//! `ceil(pv_size / 16)` lines of 64 bytes: 32 KiB + 128 KiB (`pv_array_bts` 9) or 32 + 64 KiB (10)
//! per worker for a `lut_word_length` 12 table, 8 + 128 KiB (7) or 8 + 64 KiB (8) at 11. All
//! workers' copies are one buffer of the calling thread (`Copies`; with each worker allocating its
//! own copies, glibc's per-thread arenas kept 8-19 MiB more at the end of a search), not zeroed by
//! it: each worker zeroes its own lines at the start of its part, and every line is written. So a
//! fill of P workers writes P times that, all resident: at `lut_word_length` 12, 2.5 / 1.5 MiB at 16
//! threads, 20 / 12 MiB at 128 and 40 / 24 MiB at 255 (the parts are at most 255); at 11, 2.1 /
//! 1.1, 17 / 9 and 34 / 18 MiB. On Linux with glibc the buffer is kept for the next fill and never
//! freed (`SPARE_LINES`: freeing it raised glibc's mmap threshold and with it the search's max RSS),
//! with a capacity of at least 32 MiB + 4 KiB (`lines_capacity`), so that glibc maps it and a buffer
//! that is freed after all (one too small for a fill, or one whose fill found the lock taken at its
//! end) never moves the threshold; the lines no fill wrote are never touched and not resident. What
//! the process keeps, for `n` threads and `s` the largest per-worker size of its tables: resident,
//! at most `n x s` from the fills on the whole pool plus, with `LOSAT_X_QBATCHPAR`, `m x s` for each
//! slot of `m` threads that filled a table while another fill ran (the slots' threads add up to
//! `n`, so below `2 n x s` in all), and a fill that finds the lock taken can add one more kept
//! buffer; address space, 32 MiB + 4 KiB per kept buffer (`n x s` when that is larger: 255 threads
//! at 160 KiB). At 64 threads and 160 KiB that is 10 MiB resident in one mapping of 32 MiB + 4 KiB,
//! and with four slots of 16 threads 17.5 MiB resident in four. Elsewhere (macOS, Windows, musl,
//! wasm32) the buffer has the fill's size and is freed at the end of the fill: those allocators
//! return large blocks themselves, and there is no threshold to keep low.
//!
//! EXPERIMENT (LOSAT_X_MBLUTONCE; in the strict set since 2026-10-18: on unless it is `0`, or
//! `LOSAT_X_DEFAULT=0` and it is not set): `fill`
//! hands the table to the read-once fill of `x_mb_lookup_once.rs` (its first dispatch), in which
//! each query position is read by one worker and each word handed to the owner of its cell,
//! instead of every worker reading every word as above. That fill counts into the private
//! copies of this module (`Copies`, `Unready::zeroed`, `flush_private`).

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
/// (LOSAT_X_PAIRPAR, LOSAT_X_PAIRPARSHADOW, LOSAT_X_MBLUTONCE)
fn switches() -> (bool, bool, bool) {
    static SWITCHES: OnceLock<(bool, bool, bool)> = OnceLock::new();
    *SWITCHES.get_or_init(|| {
        (
            crate::utils::xswitch::on("LOSAT_X_PAIRPAR"),
            std::env::var_os("LOSAT_X_PAIRPARSHADOW").is_some(),
            crate::utils::xswitch::on("LOSAT_X_MBLUTONCE"),
        )
    })
}

/// No NCBI counterpart: LOSAT_TIMING is set (read once).
pub(super) fn timing_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_TIMING").is_some())
}

// No NCBI counterpart (LOSAT_TIMING only): the `[TIMING] x_mb_lookup_*` counters of every megablast
// table build of the process (not the shadow's second build).
static BUILD_NANOS: AtomicU64 = AtomicU64::new(0);
static ALLOC_NANOS: AtomicU64 = AtomicU64::new(0);
static LOOP_NANOS: AtomicU64 = AtomicU64::new(0);
static POOL_FILL_NANOS: AtomicU64 = AtomicU64::new(0);
static CHAIN_NANOS: AtomicU64 = AtomicU64::new(0);
static BUILDS: AtomicU64 = AtomicU64::new(0);
static BUILDS_ON_POOL: AtomicU64 = AtomicU64::new(0);
static BUILDS_ONCE: AtomicU64 = AtomicU64::new(0);

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

    /// The end of the build: adds its intervals to the `[TIMING]` counters. `filled`: what the pool
    /// fill did (the default when the loop filled the table).
    pub(super) fn record(self, on_pool: bool, filled: FillStats) {
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
        BUILDS_ONCE.fetch_add(u64::from(filled.once), r);
    }
}

/// No NCBI counterpart: shadow counters of LOSAT_X_MBLUTONCE (tables whose second build was the
/// read-once fill, and their elements).
static SHADOW_ONCE_TABLES: AtomicU64 = AtomicU64::new(0);
static SHADOW_ONCE_CELLS: AtomicU64 = AtomicU64::new(0);

/// No NCBI counterpart: printed by `main` at exit. With LOSAT_X_MBLUTONCE and
/// LOSAT_X_PAIRPARSHADOW, for the tables whose second build was the read-once fill, `[X_SHADOW]
/// MBLUTONCE tables=<compared> cells=<elements compared> diffs=0 identical` (a difference panics
/// in `shadow_check`, so the line only ever says 0). With LOSAT_TIMING, when at
/// least one megablast table was built, the `[TIMING] x_mb_lookup_build` and `[TIMING]
/// x_mb_lookup_split` lines (`x_`: end-of-run counters of an experiment, as `[TIMING]
/// x_pairpar_seed`), and `[TIMING] x_mb_lookup_once` when the read-once fill built one
/// (`x_mb_lookup_once::print_timing`). They print counters and change no value NCBI computes.
pub(super) fn print_summary() {
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
    let (_, shadow, once) = switches();
    if shadow && once {
        eprintln!(
            "[X_SHADOW] MBLUTONCE tables={} cells={} diffs=0 identical",
            load(&SHADOW_ONCE_TABLES),
            load(&SHADOW_ONCE_CELLS)
        );
    }
    if !timing_on() || load(&BUILDS) == 0 {
        return;
    }
    let s = |c: &AtomicU64| load(c) as f64 / 1e9;
    eprintln!(
        "[TIMING] x_mb_lookup_build: {:.3}s (builds={} parallel={} once={})",
        s(&BUILD_NANOS),
        load(&BUILDS),
        load(&BUILDS_ON_POOL),
        load(&BUILDS_ONCE)
    );
    eprintln!(
        "[TIMING] x_mb_lookup_split: alloc={:.4}s loop={:.4}s pool_fill={:.4}s chain={:.4}s",
        s(&ALLOC_NANOS),
        s(&LOOP_NANOS),
        s(&POOL_FILL_NANOS),
        s(&CHAIN_NANOS)
    );
    super::x_mb_lookup_once::print_timing();
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

/// No NCBI counterpart (LOSAT_X_MBLUTONCE): sixteen 32-bit words on a 64-byte line of their own.
#[derive(Clone, Copy)]
#[repr(C, align(64))]
pub(super) struct Line([u32; 16]);

const _: () = assert!(std::mem::size_of::<Line>() == 64 && std::mem::align_of::<Line>() == 64);
const _: () = assert!(std::mem::size_of::<PvArrayType>() == std::mem::size_of::<u32>());

/// No NCBI counterpart (LOSAT_X_MBLUTONCE): Linux with glibc, where the buffers of the copies are
/// kept from fill to fill (`SPARE_LINES`) and given at least 32 MiB + 4 KiB (`lines_capacity`).
const GLIBC: bool = cfg!(all(target_os = "linux", target_env = "gnu"));

/// No NCBI counterpart (LOSAT_X_MBLUTONCE): a list of kept buffers (`SPARE_LINES`; the tests use
/// lists of their own).
type SpareLines = std::sync::Mutex<Vec<Vec<Line>>>;

/// No NCBI counterpart (LOSAT_X_MBLUTONCE): the line buffers of finished fills, kept for the next
/// fills and never freed, on Linux with glibc only (`GLIBC`). glibc raises its mmap threshold to the
/// size of a mapped block of up to 32 MiB when the block is freed, after which the search's other
/// allocations of up to that size stay in malloc's arenas: on the workstation (megablast of 31,168
/// cDNAs at 64 threads, 6 MiB of copies per fill) freeing the buffer after each fill gave max RSS
/// 263-305 MiB against 221-225 without the switch, and 228-232 with the threshold fixed
/// (`MALLOC_MMAP_THRESHOLD_`). A fill takes the largest kept buffer when it is large enough and
/// gives it back at its end. A kept buffer is never grown: `Line` is aligned to 64 bytes, more than
/// malloc's 16, so std's `realloc` of it allocates a new block, copies the old one on the calling
/// thread and frees it (`realloc_fallback`; glibc's `mremap` is never reached), which is the
/// threshold effect again. A buffer too small for a fill is freed and a new one made instead, and
/// every buffer is made with at least 32 MiB + 4 KiB (`lines_capacity`), above glibc's largest
/// dynamic threshold: freeing it never moves the threshold, and a fill needs a larger buffer only
/// above 32 MiB of copies (205 workers or more at 160 KiB each). The lock is only ever tried, never
/// waited for (`try_lock`: no fill waits for another's bookkeeping): a fill that finds it taken
/// uses a new buffer, which it keeps only if the lock is free at its end, so each such fill can add
/// one kept buffer. Elsewhere nothing is kept: macOS, Windows and musl return large blocks to the
/// system when they are freed, and wasm32's dlmalloc has no such threshold, so a kept buffer would
/// only hold memory (in the Web engine for the life of a reactor reused across searches).
static SPARE_LINES: SpareLines = std::sync::Mutex::new(Vec::new());

/// No NCBI counterpart (LOSAT_X_MBLUTONCE): the capacity, in lines, of a new buffer for `n` lines. On
/// Linux with glibc at least glibc's largest dynamic mmap threshold (`DEFAULT_MMAP_THRESHOLD_MAX`, 32
/// MiB on 64-bit) and one page more, the rule of `x_frame_capacity` in
/// `tblastn/x_subject_translate.rs` (LOSAT_X_TBNSTORE): glibc maps such a block and unmaps it when it
/// is freed without moving its threshold, and the pages beyond the lines the fills write are never
/// touched, so not resident. Elsewhere `n`.
pub(super) fn lines_capacity(n: usize) -> usize {
    if GLIBC {
        n.max(((32 << 20) + 4096) / std::mem::size_of::<Line>())
    } else {
        n
    }
}

/// No NCBI counterpart (LOSAT_X_MBLUTONCE): an empty buffer with room for `n` lines. On Linux with
/// glibc the largest buffer kept in `spare` when its lock is free and the buffer is large enough;
/// one too small is freed here, never grown (`SPARE_LINES`). Otherwise a new buffer of
/// `lines_capacity(n)` lines.
fn take_lines(spare: &SpareLines, n: usize) -> Vec<Line> {
    if GLIBC {
        let kept = spare.try_lock().ok().and_then(|mut spare| {
            let largest = (0..spare.len()).max_by_key(|&i| spare[i].capacity())?;
            Some(spare.swap_remove(largest))
        });
        if let Some(mut lines) = kept {
            if lines.capacity() >= n {
                lines.clear();
                return lines;
            }
            // Too small: freed (after the lock is released), not grown.
            drop(lines);
        }
    }
    Vec::with_capacity(lines_capacity(n))
}

/// No NCBI counterpart (LOSAT_X_MBLUTONCE): the private copies of all workers in one buffer of the
/// calling thread, not zeroed by it (one buffer per fill, on Linux with glibc kept from fill to fill:
/// `SPARE_LINES`; so the workers' own malloc arenas take nothing), `stride` lines per worker: first
/// the copy of `helper_array` (`counter_lines`), then that of `pv_array`. Each worker's lines start
/// on a line boundary and are a whole number of lines, so no two workers' copies share a line; each
/// worker zeroes its own lines (`Unready::zeroed`). `spare`: where the buffer comes from and goes back to (`SPARE_LINES`).
pub(super) struct Copies {
    lines: Vec<Line>,
    counter_lines: usize,
    stride: usize,
    parts: usize,
    spare: &'static SpareLines,
}

impl Copies {
    pub(super) fn new(parts: usize, helper_len: usize, pv_len: usize) -> Self {
        Self::with_spare(&SPARE_LINES, parts, helper_len, pv_len)
    }

    /// `new` with the kept buffers of `spare` (the tests use lists of their own).
    fn with_spare(
        spare: &'static SpareLines,
        parts: usize,
        helper_len: usize,
        pv_len: usize,
    ) -> Self {
        let (counter_lines, pv_lines) = (helper_len.div_ceil(16), pv_len.div_ceil(16));
        let stride = counter_lines + pv_lines;
        let n = stride * parts;
        Self {
            lines: if n == 0 {
                Vec::new()
            } else {
                take_lines(spare, n)
            },
            counter_lines,
            stride,
            parts,
            spare,
        }
    }

    /// The workers' copies, not yet zeroed, in part order (each `stride` lines).
    pub(super) fn parts(&mut self) -> impl Iterator<Item = Unready<'_>> {
        let (stride, counter_lines, parts) = (self.stride, self.counter_lines, self.parts);
        let mut rest = &mut self.lines.spare_capacity_mut()[..stride * parts];
        (0..parts).map(move |_| {
            let (mine, tail) = std::mem::take(&mut rest).split_at_mut(stride);
            rest = tail;
            Unready {
                lines: mine,
                counter_lines,
            }
        })
    }
}

impl Drop for Copies {
    /// On Linux with glibc gives the buffer back for the next fill when the lock is free
    /// (`SPARE_LINES`); otherwise the buffer is freed with `lines`.
    fn drop(&mut self) {
        if GLIBC && self.lines.capacity() > 0 {
            if let Ok(mut spare) = self.spare.try_lock() {
                spare.push(std::mem::take(&mut self.lines));
            }
        }
    }
}

/// No NCBI counterpart (LOSAT_X_MBLUTONCE): one worker's copies before it zeroes them.
pub(super) struct Unready<'a> {
    lines: &'a mut [std::mem::MaybeUninit<Line>],
    counter_lines: usize,
}

impl<'a> Unready<'a> {
    /// Zeroes every line (on the calling thread, the worker's) and gives the two copies.
    pub(super) fn zeroed(self) -> Private<'a> {
        for line in self.lines.iter_mut() {
            line.write(Line([0; 16]));
        }
        // SAFETY: every element was written just above, and `MaybeUninit<Line>` has the layout of
        // `Line`; the slice keeps its lifetime and length.
        let lines: &'a mut [Line] =
            unsafe { &mut *(self.lines as *mut [std::mem::MaybeUninit<Line>] as *mut [Line]) };
        let (counters, pv) = lines.split_at_mut(self.counter_lines);
        Private {
            counters: Lines { lines: counters },
            pv: Lines { lines: pv },
        }
    }
}

/// No NCBI counterpart (LOSAT_X_MBLUTONCE): a worker's private copy of `helper_array` or
/// `pv_array` (element `i` in line `i / 16`).
pub(super) struct Lines<'a> {
    lines: &'a mut [Line],
}

impl Lines<'_> {
    #[inline(always)]
    pub(super) fn get(&self, i: usize) -> u32 {
        self.lines[i >> 4].0[i & 15]
    }

    #[inline(always)]
    pub(super) fn slot(&mut self, i: usize) -> &mut u32 {
        &mut self.lines[i >> 4].0[i & 15]
    }
}

/// No NCBI counterpart (LOSAT_X_MBLUTONCE): the private copies of one worker.
pub(super) struct Private<'a> {
    pub(super) counters: Lines<'a>,
    pub(super) pv: Lines<'a>,
}

#[cfg(test)]
impl Private<'_> {
    /// Whether either copy holds a non-zero element (the tests' check that the copies took values).
    pub(super) fn has_values(&self) -> bool {
        let any = |l: &Lines<'_>| l.lines.iter().any(|line| line.0.iter().any(|&v| v != 0));
        any(&self.counters) || any(&self.pv)
    }
}

/// No NCBI counterpart (scheduling only): what one fill did. `once`: the read-once fill
/// (LOSAT_X_MBLUTONCE, `x_mb_lookup_once.rs`) made the table; for `[TIMING] x_mb_lookup_build`'s
/// `once=` and the shadow's counters. In tests also
/// `parts_with_values`, the number of parts (owners for the read-once fill) whose private copies
/// held a non-zero element before their flush.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct FillStats {
    pub(super) once: bool,
    #[cfg(test)]
    pub(super) parts_with_values: usize,
}

#[derive(Clone, Copy)]
pub(super) struct SharedTables {
    pub(super) hashtable: Shared<u32>,
    pub(super) next_pos: Shared<u32>,
    pub(super) pv_array: Shared<PvArrayType>,
    pub(super) pv_array_bts: usize,
    pub(super) helper_array: Shared<u32>,
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

/// No NCBI counterpart (LOSAT_X_MBLUTONCE): adds the counters and ORs the presence-vector words
/// of `part`'s own blocks from its private copies into `helper_array` and `pv_array` (module
/// documentation: these elements are written by no other worker, and their copies hold exactly
/// the loop's increments and bits).
pub(super) fn flush_private(
    part: u8,
    owner: &[u8],
    block_shift: usize,
    t: SharedTables,
    x: &Private<'_>,
) {
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
/// Returns what the fill did (the read-once fill of LOSAT_X_MBLUTONCE or not). `resident`: `hashtable` is a kept table whose pages are mapped
/// (LOSAT_X_MBLUTKEEP), on which the read-once fill skips its prefault round.
pub(super) fn fill(
    parts: usize,
    words: &Words<'_>,
    tables: Tables<'_>,
    in_pool: bool,
    resident: bool,
) -> FillStats {
    let (_, _, once) = switches();
    fill_mode(once, true, parts, words, tables, in_pool, resident)
}

/// No NCBI counterpart: `fill` with LOSAT_X_MBLUTONCE given as `once` (the tests take each). `timed`: a build of the search (not the shadow's second build),
/// whose read-once phases LOSAT_TIMING's `[TIMING] x_mb_lookup_once` line counts. `resident`:
/// as `fill`'s.
fn fill_mode(
    once: bool,
    timed: bool,
    parts: usize,
    words: &Words<'_>,
    tables: Tables<'_>,
    in_pool: bool,
    resident: bool,
) -> FillStats {
    // Dispatch point of LOSAT_X_MBLUTONCE (first): each query position read once, the words
    // handed to the owners of their cells (`x_mb_lookup_once.rs`, which counts into private
    // copies); off, the fill of this module.
    if once {
        let _stats = super::x_mb_lookup_once::fill(parts, words, tables, in_pool, timed, resident);
        return FillStats {
            once: true,
            #[cfg(test)]
            parts_with_values: _stats.owners_with_values,
        };
    }
    fill_with(parts, words, tables, in_pool)
}

/// No NCBI counterpart: `fill` without the read-once fill (the shared arrays are updated by the
/// owners of the cells).
fn fill_with(parts: usize, words: &Words<'_>, tables: Tables<'_>, in_pool: bool) -> FillStats {
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
    let ran = in_pool && {
        rayon::scope(|scope| {
            for part in 1..parts {
                scope.spawn(move |_| fill_part(part as u8, owner, block_shift, words, shared));
            }
            fill_part(0, owner, block_shift, words, shared);
        });
        true
    };
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    let ran = {
        let _ = in_pool;
        false
    };
    if !ran {
        for part in 0..parts {
            fill_part(part as u8, owner, block_shift, words, shared);
        }
    }
    FillStats {
        once: false,
        #[cfg(test)]
        parts_with_values: 0,
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
/// With LOSAT_X_MBLUTONCE the second build is the read-once fill; otherwise the fill of this
/// module.
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
    let (_, _, once) = switches();
    shadow_check_mode(
        once,
        words,
        pv_array_bts,
        hashtable,
        next_pos,
        pv_array,
        helper_array,
        longest_chain,
    );
}

/// No NCBI counterpart: `shadow_check` with LOSAT_X_MBLUTONCE given as `once`.
#[allow(clippy::too_many_arguments)]
fn shadow_check_mode(
    once: bool,
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
    let x_stats = fill_mode(
        once,
        false,
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
        false,
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
                "LOSAT_X_PAIRPARSHADOW: megablast lookup {name}[{i}] differs: reference {} partitioned {} (parts {parts}, lut_word_length {}, ascending_cells {}, once {once})",
                reference[i], partitioned[i], words.lut_word_length, words.ascending_cells
            );
        }
    }
    assert_eq!(
        longest_chain, x_longest_chain,
        "LOSAT_X_PAIRPARSHADOW: megablast lookup longest_chain differs (parts {parts}, once {once})"
    );
    SHADOW_TABLES.fetch_add(1, Ordering::Relaxed);
    let cells = (hashtable.len() + next_pos.len() + pv_array.len() + helper_array.len()) as u64;
    SHADOW_CELLS.fetch_add(cells, Ordering::Relaxed);
    if x_stats.once {
        SHADOW_ONCE_TABLES.fetch_add(1, Ordering::Relaxed);
        SHADOW_ONCE_CELLS.fetch_add(cells, Ordering::Relaxed);
    }
}

#[cfg(test)]
pub(super) mod tests {
    use super::super::build_two_stage_lookup;
    use super::*;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1081-1091
    // ```c
    //          mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
    //          mb_lt->hashtable[ecode] = index;
    // ```
    // Random queries (with ambiguity codes and masks) through the original loop and through
    // the partitioned fill with 1..8 workers: every table element must be the same.
    pub(crate) struct Lcg(pub(crate) u64);
    impl Lcg {
        pub(crate) fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        pub(crate) fn below(&mut self, n: u64) -> u64 {
            self.next() % n.max(1)
        }
    }

    pub(crate) fn random_case(
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

    /// The fill of this module without the private copies and without the read-once fill (the
    /// loop's statements, part by part; checked against the loop below), for the tests of
    /// `x_mb_lookup_once.rs`.
    pub(crate) fn fill_old(parts: usize, words: &Words<'_>, tables: Tables<'_>) {
        fill_mode(false, false, parts, words, tables, false, false);
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
                // The fill of this module (the read-once fill has its own tests in
                // `x_mb_lookup_once.rs`).
                fill_mode(
                    false,
                    false,
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
                    fill_mode(
                        false,
                        false,
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
                        false,
                    );
                    (hashtable, next_pos, pv_array)
                });
                assert_eq!(hashtable, mb.hashtable, "case {case} threads {threads}");
                assert_eq!(next_pos, mb.next_pos, "case {case} threads {threads}");
                assert_eq!(pv_array, mb.pv_array, "case {case} threads {threads}");
            }
        }
    }

    // The tables of one case, built by the fill of this module in one part (the loop's
    // statements in reading order; the tests above check it against `build_two_stage_lookup`),
    // from start values `helper0` / `pv0` of the counters and presence-vector words.
    struct Built {
        hashtable: Vec<u32>,
        next_pos: Vec<u32>,
        pv_array: Vec<PvArrayType>,
        helper_array: Vec<u32>,
        stats: FillStats,
    }

    fn build_with(
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
            stats: FillStats::default(),
        };
        assert_eq!(b.pv_array.len(), pv_len);
        b.stats = fill_mode(
            false,
            false,
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
            false,
        );
        b
    }

    fn assert_same(a: &Built, b: &Built, what: &str) {
        assert_eq!(a.hashtable, b.hashtable, "{what}: hashtable");
        assert_eq!(a.next_pos, b.next_pos, "{what}: next_pos");
        assert_eq!(a.pv_array, b.pv_array, "{what}: pv_array");
        assert_eq!(a.helper_array, b.helper_array, "{what}: helper_array");
    }

    // A query with long unmasked repeats (many words per cell, so counters grow) and random bases.
    pub(crate) fn repeat_case(rng: &mut Lcg, len: usize) -> Vec<u8> {
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

    // The fill on real pools of 2, 3, 5 and 16 threads (and off the pool with the same
    // parts): every element of the four arrays and the longest chain equal the one-part fill's
    // and `build_two_stage_lookup`'s, for random cases (ambiguity codes, masks, several queries,
    // `lut_word_length` 4..12, both `ascending_cells`) and for `lut_word_length` 11 and 12 tables
    // (`pv_array_bts` 8 to 10).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn fill_gives_the_loops_tables_on_pools() {
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
        let mut seen_large_tables = 0;
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
                seen_large_tables += 1;
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
            let loop_built = build_with(1, false, &words, shape, &helper0, &pv0);
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
                let pooled =
                    pool.install(|| build_with(threads, true, &words, shape, &helper0, &pv0));
                assert_same(&loop_built, &pooled, &format!("case {case} pool {threads}"));
                let serial = build_with(threads, false, &words, shape, &helper0, &pv0);
                assert_same(
                    &loop_built,
                    &serial,
                    &format!("case {case} parts {threads}"),
                );
            }
        }
        assert!(
            seen_large_tables >= 3,
            "the lut 11 / 12 cases ({seen_large_tables})"
        );
    }

    // Each worker's copies are a whole number of 64-byte lines starting on a line boundary, the
    // workers' copies are disjoint and in part order, and zeroed.
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
            let mut copies = Copies::new(parts, helper_len, pv_len);
            let stride = helper_len.div_ceil(16) + pv_len.div_ceil(16);
            assert!(copies.lines.capacity() >= stride * parts);
            let base = copies.lines.as_ptr() as usize;
            assert_eq!(base % 64, 0);
            let mut seen = 0;
            let mut ends = Vec::new();
            for (p, x) in copies.parts().enumerate() {
                let mut x = x.zeroed();
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
                assert!(!x.has_values());
                *x.counters.slot(helper_len - 1) = 1;
                *x.pv.slot(pv_len - 1) = 1;
                assert!(x.has_values());
                ends.push((c0, v0 + x.pv.lines.len() * 64));
                seen += 1;
            }
            assert_eq!(seen, parts);
            for pair in ends.windows(2) {
                assert!(pair[0].1 <= pair[1].0, "parts {parts}: copies overlap");
            }
        }
    }

    // The kept buffers, on a list of this test's own (the module's other tests run at
    // the same time on `SPARE_LINES`). On Linux with glibc a fill takes the largest kept buffer when
    // it is large enough; one too small is freed and the fill gets a new buffer of
    // `lines_capacity(n)` lines (32 MiB + 4 KiB or more), not the old one grown (std's aligned
    // realloc frees the old block, which moves glibc's mmap threshold); `Drop` gives the buffer back,
    // and frees it without waiting when the lock is taken. Elsewhere nothing is kept.
    #[test]
    fn kept_buffers_are_taken_largest_first_never_grown_and_given_back() {
        static SPARE: SpareLines = std::sync::Mutex::new(Vec::new());
        let kept = || -> Vec<(usize, usize)> {
            SPARE
                .lock()
                .unwrap()
                .iter()
                .map(|l| (l.as_ptr() as usize, l.capacity()))
                .collect()
        };
        let at = |c: &Copies| (c.lines.as_ptr() as usize, c.lines.capacity());
        // 4 parts of 2 + 16 lines.
        let (parts, helper_len, pv_len) = (4usize, 32usize, 256usize);
        let n = parts * (helper_len / 16 + pv_len / 16);
        let mapped = ((32 << 20) + 4096) / 64;
        assert_eq!(lines_capacity(n), if GLIBC { mapped } else { n });
        assert_eq!(lines_capacity(mapped + 7), mapped + 7);

        // The first fill makes a new buffer of `lines_capacity(n)` lines.
        let copies = Copies::with_spare(&SPARE, parts, helper_len, pv_len);
        let first = at(&copies);
        assert!(first.1 >= lines_capacity(n));
        drop(copies);
        if !GLIBC {
            // Nothing is kept, and a buffer in the list is not taken.
            assert!(kept().is_empty(), "kept off glibc");
            SPARE.lock().unwrap().push(Vec::with_capacity(4 * n));
            let before = kept();
            drop(Copies::with_spare(&SPARE, parts, helper_len, pv_len));
            assert_eq!(kept(), before, "taken or kept off glibc");
            SPARE.lock().unwrap().clear();
            return;
        }
        assert_eq!(kept(), [first], "Drop gives the buffer back");
        // The next fill takes it, and gives it back.
        let copies = Copies::with_spare(&SPARE, parts, helper_len, pv_len);
        assert_eq!(at(&copies), first);
        assert!(kept().is_empty());
        drop(copies);
        assert_eq!(kept(), [first]);
        // With a smaller buffer kept too, the fill takes the largest.
        SPARE.lock().unwrap().push(Vec::with_capacity(n));
        let copies = Copies::with_spare(&SPARE, parts, helper_len, pv_len);
        assert_eq!(at(&copies), first, "not the largest kept buffer");
        let small = kept();
        assert_eq!(small.len(), 1);
        assert!(small[0].1 >= n && small[0].1 < first.1);
        // The lock taken at the end of a fill: the buffer is freed, not kept, and nothing waits.
        {
            let _taken = SPARE.lock().unwrap();
            drop(copies);
        }
        assert_eq!(kept(), small, "kept although the lock was taken");
        // A fill larger than the only kept buffer: that one is freed, and the fill gets a new
        // buffer of `lines_capacity` lines, not the old one grown to its size.
        let copies = Copies::with_spare(&SPARE, 2 * parts, helper_len, pv_len);
        assert!(
            copies.lines.capacity() >= lines_capacity(2 * n),
            "grown to {} lines instead of a new buffer of {}",
            copies.lines.capacity(),
            lines_capacity(2 * n)
        );
        assert!(kept().is_empty(), "the buffer too small is kept");
        drop(copies);
        assert_eq!(kept().len(), 1);
        // More than 32 MiB of copies (205 parts of 160 KiB, not touched here): the kept buffer is
        // too small, and the fill gets a new one of its own size.
        let copies = Copies::with_spare(&SPARE, 205, 8192, 32768);
        assert_eq!(copies.stride, 2560);
        assert!(copies.lines.capacity() >= 205 * 2560);
        assert!(kept().is_empty());
        drop(copies);
        assert_eq!(kept().len(), 1);
        SPARE.lock().unwrap().clear();
    }

    // Many parts: off the pool 64 / 254 / 255 / 256 / 300 parts (the last two clamp to 255), on a
    // 64-thread pool 64 and 255 parts (more parts than threads, and at lut 8 more parts than
    // blocks); lut 8 / 11 / 12; the fill equals the loop.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn fill_with_many_parts_gives_the_loops_tables() {
        let mut rng = Lcg(0x00C0_FFEE);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(64)
            .build()
            .unwrap();
        for (lut, asc, approx) in [
            (8usize, true, 30_000usize),
            (11, false, 1_000_000),
            (12, true, 400_000),
        ] {
            let q = repeat_case(&mut rng, 25_000);
            let queries = vec![q];
            let refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
            let offsets = vec![0i32];
            let masks = vec![Vec::new()];
            let reference =
                build_two_stage_lookup(&refs, &offsets, lut + 4, lut, &masks, None, 0, approx, asc);
            let mb = &reference.mb_lookup;
            let shape = (
                mb.hashtable.len(),
                mb.next_pos.len(),
                mb.pv_array.len(),
                mb.pv_array_bts,
            );
            let words = Words {
                queries_blastna: &refs,
                query_offsets: &offsets,
                query_masks: &masks,
                word_length: lut + 4,
                lut_word_length: lut,
                kmer_mask: (1u64 << (2 * lut)) - 1,
                ascending_cells: asc,
            };
            let helper0 = vec![0u32; (mb.hashtable.len() >> 11).max(1)];
            let pv0 = vec![0 as PvArrayType; mb.pv_array.len()];
            let loop_built = build_with(1, false, &words, shape, &helper0, &pv0);
            assert_eq!(loop_built.hashtable, mb.hashtable);
            assert_eq!(loop_built.next_pos, mb.next_pos);
            assert_eq!(loop_built.pv_array, mb.pv_array);
            for parts in [64usize, 254, 255, 256, 300] {
                let b = build_with(parts, false, &words, shape, &helper0, &pv0);
                assert_same(&loop_built, &b, &format!("lut {lut} parts {parts}"));
            }
            for parts in [64usize, 255] {
                let b = pool.install(|| build_with(parts, true, &words, shape, &helper0, &pv0));
                assert_same(&loop_built, &b, &format!("lut {lut} pool 64 parts {parts}"));
            }
        }
    }

    // Every counter and presence-vector word is flushed by exactly one part, its block's owner,
    // for every reachable table shape (lut 4..12; `approx_table_entries` below, between and above
    // the cutoffs that halve the presence vector; with and without the database filter) and
    // P = 1..=17, 31, 64, 127, 128, 254, 255, 256 and 300 (the last two clamp): each part's
    // copies hold `part + 1` everywhere, so after all flushes a counter
    // holds its owner + 1 (the flush adds to 0) and a word its owner + 1 (one OR). A double flush,
    // a missed one or a wrong owner fails. Also: no presence-vector word spans two blocks, and
    // the copies are `P x (ceil(hashsize / 32768) + ceil(pv_size / 16))` lines, at most 160 KiB per
    // part (the module documentation's memory sentence).
    #[test]
    fn every_counter_and_word_is_flushed_once_by_its_owner() {
        let mut shapes = 0;
        for lut in 4usize..=12 {
            let mut worst = 0usize;
            for approx in [1_000usize, 20_000, 1_000_000] {
                for db in [false, true] {
                    let table = 1usize << (2 * lut);
                    let (pv_len, bts) = super::super::compute_mb_pv_params(table, approx, db, lut);
                    let helper_len = table.div_ceil(2048).max(1);
                    let per_part = table.div_ceil(32768) + pv_len.div_ceil(16);
                    worst = worst.max(per_part * 64);
                    for parts_in in (1usize..=17).chain([31, 64, 127, 128, 254, 255, 256, 300]) {
                        let (block_shift, owner) = plan(table, bts, parts_in);
                        let parts = parts_in.clamp(1, 255);
                        assert!(block_shift >= bts && block_shift >= 11);
                        let mut helper = vec![0u32; helper_len];
                        let mut pv = vec![0 as PvArrayType; pv_len];
                        let mut hashtable = vec![0u32; 0];
                        let mut next_pos = vec![0u32; 0];
                        let t = SharedTables {
                            hashtable: Shared::new(&mut hashtable),
                            next_pos: Shared::new(&mut next_pos),
                            pv_array: Shared::new(&mut pv),
                            pv_array_bts: bts,
                            helper_array: Shared::new(&mut helper),
                        };
                        let mut copies = Copies::new(parts, helper_len, pv_len);
                        assert!(copies.lines.capacity() >= parts * per_part, "lut {lut}");
                        assert_eq!(copies.stride, per_part, "lut {lut}");
                        for (p, x) in copies.parts().enumerate() {
                            let x = x.zeroed();
                            for line in x.counters.lines.iter_mut().chain(x.pv.lines.iter_mut()) {
                                line.0 = [p as u32 + 1; 16];
                            }
                            flush_private(p as u8, &owner, block_shift, t, &x);
                        }
                        let what = format!("lut {lut} approx {approx} db {db} P {parts_in}");
                        for (h, &value) in helper.iter().enumerate() {
                            let block = (h << 11) >> block_shift;
                            let o = owner[block.min(owner.len() - 1)] as u32;
                            assert_eq!(value, o + 1, "{what} counter {h}");
                        }
                        for (w, &value) in pv.iter().enumerate() {
                            let block = (w << bts) >> block_shift;
                            let o = owner[block.min(owner.len() - 1)] as u32;
                            assert_eq!(value, o + 1, "{what} word {w}");
                            assert_eq!(block, (((w + 1) << bts) - 1) >> block_shift, "{what}");
                        }
                        shapes += 1;
                    }
                }
            }
            if lut >= 11 {
                assert!(worst <= 160 * 1024, "lut {lut}: {worst} bytes per part");
            }
        }
        assert!(shapes > 1000, "{shapes}");
    }

    // The shadow (`LOSAT_X_PAIRPARSHADOW`) with the second build of this module's fill: equal
    // tables pass, and a difference in any one of the four arrays or the longest chain stops it.
    #[test]
    fn shadow_compares_every_array() {
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
                false,
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
                stats: b.stats,
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
