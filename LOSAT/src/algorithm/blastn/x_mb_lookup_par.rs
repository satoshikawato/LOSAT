//! EXPERIMENT (PAIRPAR, graduated on 2026-10-22: no variable): the megablast lookup table of
//! `build_mb_lookup` filled by the threads of the search pool, each thread owning a fixed set
//! of table cells.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1032-1033,1043-1045,1081-1091
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1104-1106
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
//! words on x86_64, for a prefetch, as `x_fill_mb_queued` does). Workers write disjoint
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
//! part one after another; with per-thread private copies (a padded fill since removed) 0.8 M. The stores of `next_pos` still take
//! lines from other cores (each worker writes the offsets of its own words, which interleave
//! with the other workers'; 32 M such stores in that run), but stores wait in the store buffer;
//! writing `next_pos` afterwards by spans of offsets removed them and cost more than it saved
//! (lists of 8 bytes a word and their growth), so it is not done.
//! Each worker of the read-once fill counts into a private copy of `helper_array` and sets its bits in a
//! private copy of `pv_array`, each a whole number of 64-byte lines starting on a line boundary
//! (`Lines`), zeroed in its first job, and after its last chunk adds its own blocks' counters into `helper_array`
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
//! at most `n x s` from the fills on the whole pool, and a fill that finds the lock taken can add
//! one more kept buffer; address space, 32 MiB + 4 KiB per kept buffer (`n x s` when that is
//! larger: 255 threads at 160 KiB). At 64 threads and 160 KiB that is 10 MiB resident in one
//! mapping of 32 MiB + 4 KiB. Elsewhere (macOS, Windows, musl, wasm32) the buffer has the fill's
//! size and is freed at the end of the fill: those allocators return large blocks themselves, and
//! there is no threshold to keep low.
//!
//! EXPERIMENT (LOSAT_X_MBLUTONCE; in the strict set since 2026-10-18: on unless it is `0`): `fill`
//! hands the table to the read-once fill of `x_mb_lookup_once.rs` (its first dispatch), in which
//! each query position is read by one worker and each word handed to the owner of its cell,
//! instead of every worker reading every word as above. That fill counts into the private
//! copies of this module (`Copies`, `Unready::zeroed`, `flush_private`).
//!
//! EXPERIMENT (LOSAT_X_MB24): a table of 3-byte elements (`Tables::x_h24`, `x_mb_narrow.rs`) runs
//! its own instance of `fill_part` (`H24`), which reads and writes every head and link as exactly
//! its 3 bytes: the cell after the last of a block and the offset after a worker's word may be
//! another worker's, so no 4-byte access may cover them. The plan, the reading order and the
//! updates are the same, so the table is the same.

use super::{build_unmasked_ranges, x_range_words, MaskedInterval, PvArrayType, PV_ARRAY_MASK};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1081-1091
/// `kCompressionFactor`, the number of cells that share one counter.
pub(super) const K_COMPRESSION_FACTOR_SHIFT: usize = 11;

/// No NCBI counterpart: reads LOSAT_X_MBLUTONCE once; it does not change any value NCBI computes.
fn once() -> bool {
    static ONCE: OnceLock<bool> = OnceLock::new();
    *ONCE.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_MBLUTONCE"))
}

/// No NCBI counterpart: LOSAT_TIMING is set (read once).
pub(super) fn timing_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_TIMING").is_some())
}

// No NCBI counterpart (LOSAT_TIMING only): the `[TIMING] x_mb_lookup_*` counters of every megablast
// table build of the process.
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

/// No NCBI counterpart: printed by `main` at exit. With LOSAT_TIMING, when at
/// least one megablast table was built, the `[TIMING] x_mb_lookup_build` and `[TIMING]
/// x_mb_lookup_split` lines (`x_`: end-of-run counters of an experiment, as `[TIMING]
/// x_pairpar_seed`), and `[TIMING] x_mb_lookup_once` when the read-once fill built one
/// (`x_mb_lookup_once::print_timing`). They print counters and change no value NCBI computes.
pub(super) fn print_summary() {
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
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
/// Dispatch test of `build_mb_lookup`: the number of workers when the search runs on a pool of
/// two threads or more and the table is the plain one (`eligible`: no database word counts, no
/// debug counters). `None` runs the original loops.
pub(super) fn parts(eligible: bool) -> Option<usize> {
    let once = once();
    if !eligible {
        crate::utils::xreach::x_reach!("LOSAT_X_PAIRPAR", false);
        if once {
            crate::utils::xreach::x_reach!("LOSAT_X_MBLUTONCE", false);
        }
        return None;
    }
    let parts = pool_threads();
    crate::utils::xreach::x_reach!("LOSAT_X_PAIRPAR", parts.is_some());
    if once {
        crate::utils::xreach::x_reach!("LOSAT_X_MBLUTONCE", parts.is_some());
    }
    parts
}

/// The tables `build_mb_lookup` fills (NCBI `BlastMBLookupTable` fields and its local
/// `helper_array`).
pub(super) struct Tables<'a> {
    pub(super) hashtable: &'a mut [u32],
    pub(super) next_pos: &'a mut [u32],
    pub(super) pv_array: &'a mut [PvArrayType],
    pub(super) pv_array_bts: usize,
    pub(super) helper_array: &'a mut [u32],
    /// EXPERIMENT (LOSAT_X_MB24): `hashtable` (`kmer_mask + 1` cells of `Words`) and `next_pos`
    /// hold 3-byte elements in their bytes (`x_mb_narrow.rs`).
    pub(super) x_h24: bool,
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
/// `once=`.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct FillStats {
    pub(super) once: bool,
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
    // The update of one word, as in the loop of `build_mb_lookup` (same statements);
    // the presence bit is set in the presence-vector word this worker owns. EXPERIMENT
    // (LOSAT_X_MB24): `H24`, on 3-byte elements, each read and written as exactly its 3 bytes
    // (another worker may own the next cell or offset).
    /// SAFETY: the worker owns the block of `bucket` and the word at `q_off_1`.
    #[inline(always)]
    unsafe fn apply<const H24: bool>(self, bucket: usize, q_off_1: u32) {
        unsafe {
            let head = if H24 {
                debug_assert!(3 * bucket + 3 <= 4 * self.hashtable.len);
                super::x_mb_narrow::load_exact(self.hashtable.ptr, bucket)
            } else {
                self.hashtable.get(bucket)
            };
            if head == 0 {
                // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:49-50
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
            if H24 {
                debug_assert!(3 * q_off_1 as usize + 3 <= 4 * self.next_pos.len);
                super::x_mb_narrow::store_exact(self.next_pos.ptr, q_off_1 as usize, head);
                super::x_mb_narrow::store_exact(self.hashtable.ptr, bucket, q_off_1);
            } else {
                self.next_pos.set(q_off_1 as usize, head);
                self.hashtable.set(bucket, q_off_1);
            }
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1032-1045
/// One worker: the word loop of `build_mb_lookup` (same reading order, same resets), applying
/// the update only to the words whose cell block `owner` gives to `part`; `H24`: on 3-byte
/// elements (LOSAT_X_MB24, `SharedTables::apply`; newest first).
fn fill_part<const H24: bool>(
    part: u8,
    owner: &[u8],
    block_shift: usize,
    words: &Words<'_>,
    t: SharedTables,
) {
    debug_assert!(!(H24 && words.ascending_cells));
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

            x_range_words(
                &seq[range_start..range_end],
                lut_word_length,
                words.kmer_mask,
                || {},
                |current_kmer, i| {
                    let pos = range_start + i;
                    let bucket = current_kmer as usize;
                    if owner[bucket >> block_shift] != part {
                        return;
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
                                if H24 {
                                    super::x_mb_narrow::at(t.hashtable.ptr, bucket)
                                } else {
                                    t.hashtable.ptr.add(bucket) as *const i8
                                },
                                core::arch::x86_64::_MM_HINT_T0,
                            );
                        }
                        let slot = x_ring_len % X_DELAY;
                        if x_ring_len >= X_DELAY {
                            let (b, q) = x_ring[slot];
                            // SAFETY: the cell block of `b` is owned by this part; `q` is the
                            // offset of one of this part's words.
                            unsafe { t.apply::<H24>(b as usize, q) };
                        }
                        x_ring[slot] = (bucket as u32, q_off_1);
                        x_ring_len += 1;
                    }
                    #[cfg(not(target_arch = "x86_64"))]
                    {
                        // SAFETY: as above.
                        unsafe { t.apply::<H24>(bucket, q_off_1) };
                    }
                },
            );
        }
    }
    #[cfg(target_arch = "x86_64")]
    {
        let pending = x_ring_len.min(X_DELAY);
        for k in x_ring_len - pending..x_ring_len {
            let (b, q) = x_ring[k % X_DELAY];
            // SAFETY: as above.
            unsafe { t.apply::<H24>(b as usize, q) };
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:74-76
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
/// Fills the tables as the word loop of `build_mb_lookup` does, with `parts` workers on the
/// search pool (or one after another, `in_pool == false`, where the unit tests run it without a
/// pool). Returns what the fill did (the read-once fill of LOSAT_X_MBLUTONCE or not). `resident`:
/// `hashtable` is a kept table whose pages are mapped (LOSAT_X_MBLUTKEEP), on which the read-once
/// fill skips its prefault round.
pub(super) fn fill(
    parts: usize,
    words: &Words<'_>,
    tables: Tables<'_>,
    in_pool: bool,
    resident: bool,
) -> FillStats {
    fill_mode(once(), true, parts, words, tables, in_pool, resident)
}

/// No NCBI counterpart: `fill` with LOSAT_X_MBLUTONCE given as `once` (the tests take each).
/// `timed`: a build of the search, whose read-once phases LOSAT_TIMING's `[TIMING]
/// x_mb_lookup_once` line counts. `resident`: as `fill`'s.
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
        let _ = super::x_mb_lookup_once::fill(parts, words, tables, in_pool, timed, resident);
        return FillStats { once: true };
    }
    // LOSAT_X_MB24: the instance for 3-byte elements.
    if tables.x_h24 {
        fill_with::<true>(parts, words, tables, in_pool)
    } else {
        fill_with::<false>(parts, words, tables, in_pool)
    }
}

/// No NCBI counterpart: `fill` without the read-once fill (the shared arrays are updated by the
/// owners of the cells).
/// `H24` (LOSAT_X_MB24): `tables` holds 3-byte elements (`Tables::x_h24`).
fn fill_with<const H24: bool>(
    parts: usize,
    words: &Words<'_>,
    tables: Tables<'_>,
    in_pool: bool,
) -> FillStats {
    debug_assert_eq!(tables.x_h24, H24);
    let table_size = if H24 {
        words.kmer_mask as usize + 1
    } else {
        tables.hashtable.len()
    };
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
                scope.spawn(move |_| {
                    fill_part::<H24>(part as u8, owner, block_shift, words, shared)
                });
            }
            fill_part::<H24>(0, owner, block_shift, words, shared);
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
            fill_part::<H24>(part as u8, owner, block_shift, words, shared);
        }
    }
    FillStats { once: false }
}

#[cfg(test)]
#[path = "x_mb_lookup_par_tests.rs"]
pub(super) mod tests;
