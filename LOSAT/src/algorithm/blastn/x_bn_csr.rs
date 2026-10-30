//! EXPERIMENT (LOSAT_X_BNCSR, in the strict set: default ON since 2026-10-25, off with
//! `LOSAT_X_BNCSR=0`): the chains of a contiguous table of `build_mb_lookup` (the two-stage table
//! of `blastn`'s word sizes 9 to 11, `lut_word_length` 8 to 11) laid out as CSR (compressed sparse
//! rows) instead of NCBI's `hashtable` / `next_pos`, as LOSAT_X_DCCSR lays out a discontiguous
//! table (`x_dc_csr.rs`, whose `Csr`, sum pass and scan `Ring` this module uses): `offsets[cell] ..
//! offsets[cell + 1]` indexes `qs`, which holds the cell's 1-based query offsets in the order its
//! chain lists them. A scan reads a cell's run as one slice, so a word with hits costs a miss on
//! `offsets` and one on its run, both prefetched ahead by `Ring`, instead of one miss per chain
//! element on `next_pos` (indexed by query offset, so in no order).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1081-1091
//! ```c
//!          if (mb_lt->hashtable[ecode] == 0) {
//!             PV_SET(pv_array, ecode, pv_array_bts);
//!          }
//!          else {
//!             helper_array[ecode/kCompressionFactor]++;
//!          }
//!          mb_lt->next_pos[index] = mb_lt->hashtable[ecode];
//!          mb_lt->hashtable[ecode] = index;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1406-1418
//! ```c
//!     Int4 q_off = lookup->hashtable[index];
//!
//!     while (q_off) {
//!         offset_pairs[i].qs_offsets.q_off   = q_off - 1;
//!         offset_pairs[i++].qs_offsets.s_off = s_off;
//!         q_off = lookup->next_pos[q_off];
//!     }
//! ```
//!
//! Why the runs are the chains. The CSR is built from the words of the fill's loop in the loop's
//! order (`each_word`: the loop of `build_mb_lookup` and `x_mb_lookup_par::fill_part`, the same
//! queries, unmasked ranges, ambiguity resets, k-mer and 1-based offset `q_off_1`, each word once).
//! A count pass adds 1 to the count of the word's cell, a sum pass (`x_dc_csr::sum_cells`) turns
//! the counts into inclusive running sums (the end of each cell's run), and a place pass, in the
//! loop's order again, writes each word into its cell's run. The two orders of a chain:
//! - newest first (`ascending_cells == false`, the megablast-type table, `-task blastn`'s table of
//!   word size 11 and its tables of 4^9 to 4^11 cells): the fill puts each word at the head of its
//!   cell's chain, so a chain lists its words in reverse reading order. The count of `cell` is
//!   `offsets[cell]`; the place pass writes each word at `--offsets[cell]`, so the first word read
//!   lands at the end of the run and the last at its start, and `offsets[cell]` is left at the
//!   run's start; `offsets[cells]` (no cell's end) is set to the number of words.
//! - reading order (`ascending_cells`, the table that stands for NCBI's small or standard table,
//!   which `build_mb_lookup` re-threads so that each chain lists its words in reading order,
//!   blast_lookup.c:74-76): the count of `cell` is `offsets[cell + 2]`, so after the sum pass
//!   `offsets[cell + 2]` is the end of the cell's run and `offsets[cell + 1]` its start (the end of
//!   the cell before, 0 for cell 0); the place pass writes each word at `offsets[cell + 1]` and
//!   moves it on, which leaves `offsets[cell + 1]` at the end of the run, so `offsets[cell] ..
//!   offsets[cell + 1]` is the run, in reading order, and `offsets[cells]` the number of words.
//!   The array has `cells + 2` elements; the table is its first `cells + 1`.
//!
//! The sum pass also makes the fill's presence vector and helper counters from the counts: the fill
//! sets the bit of a cell at its first word (`hashtable[ecode] == 0`) and adds 1 to the cell's
//! counter at every further word, i.e. bit `ecode & 31` of word `ecode >> pv_array_bts` for every
//! non-empty cell and `count - 1` on counter `ecode / kCompressionFactor` (wrapping in the sum pass
//! where the pool fill saturates; both are the plain sum, which stays below the number of words,
//! itself below 2^32); `longest_chain` is then computed from the counters by the original code. A
//! scan that reads `qs[offsets[c] .. offsets[c + 1]]` after the same presence test visits the same
//! query offsets in the same order as the chain walk, so the offset pairs, and everything made from
//! them, are the same. `contains_hit` (the masked-word test of `s_TypeOfWord`) searches the run,
//! which is sorted (`contains`).
//!
//! Scope (`plan`): a table without database word counts (the fill's `kDbFilter` path) and without
//! BLEMIR_DEBUG (whose counters the serial loop makes), of `lut_word_length` 8 to 11 (`MIN_LUT`,
//! `MAX_LUT`; megablast's tables of 4^12 cells keep their chains, v4-R R11 found no gain with a
//! subject of 5 Mb, unless LOSAT_X_MBCSR's size rule takes them: `X_MB12_LETTERS_PER_POSITION`),
//! whose word loop reads at least a quarter as many positions as the table has cells (the size rule
//! of `x_dc_csr`, `wanted`; below it the chains fit the caches and the CSR saves no miss). With
//! NCBI's choice of table (`choose_na_lookup_table`) that is: word size 11 (the default) with the
//! 4^11-cell table of a query batch of 1,048,576 positions or more (about 524 kb of query, both
//! strands); word size 10 with its 4^10-cell table from 262,144 positions; word size 9 with its
//! 4^9-cell table from 65,536 and its 4^8-cell table (reading order) from 16,384. Word sizes 12 and
//! above build 4^11-cell tables only for fewer positions than that (NCBI moves to 4^12 cells
//! first), so their tables keep the chains unless LOSAT_X_MBCSR takes them. EXPERIMENT (LOSAT_X_BNCSRPACK): a packed head (every
//! table of the scope but the reading-order one) has no size rule, so these tables (and megablast's
//! of 4^11 cells or fewer, for small query batches) are packed too. Its build writes one byte per
//! cell and two words per 64 cells (5 MiB at 4^11 cells, in `qs` and its bits) where the chains
//! take a table of four bytes per cell (16 MiB) and `next_pos`, and the batch loop's kept table
//! then never holds chains: below the rule a query batch (the first or the last of a run, beside
//! batches of the rule) built chains on the kept table, whose pages past those a packed head wrote
//! were mapped there, at the end of the run. The table is built on the search pool where the fill
//! would be (`x_mb_lookup_once::csr`), otherwise here (`serial`). Its `hashtable` and `next_pos`
//! are then empty: the CSR replaces them (`offsets` has two elements more than `hashtable`, `qs`
//! one element per word where `next_pos` has one per query position), so the table takes no more
//! memory. The batch loop's keep (LOSAT_X_MBLUTKEEP) holds no CSR: `offsets` is allocated by each
//! build and freed with the table (LOSAT_X_BNLUTDROP frees it before the traceback), and a build
//! laid out as CSR frees the table the keep holds, so the loop never holds a kept table beside it.
//! Every reader of the chains of a contiguous table goes through `MbLookupTable::x_csr`: the strip
//! scan of PAIRPAR (`x_pair_seed.rs`, `scan_strip_csr`), the reference scan of `search_query_batch`
//! (`run.rs`, its `x_on_word_csr`; the chains' `on_word` stays as it was), both reading their words
//! with hits through `x_dc_csr::Ring` (with LOSAT_X_BNCSRPACK through the ring of the head's form,
//! chosen once per table: `x_dc_csr::XHead`), and `for_each_hit`, `contains_hit` and `count_hits`
//! where the table has no chains.
//!
//! With LOSAT_TIMING `main` prints `[TIMING] x_bn_csr: builds=<CSR tables> parallel=<on a pool>
//! words=<their words>` (the build times are in `[TIMING] x_mb_lookup_build`).

use super::x_dc_csr::{
    assemble, prefetch, sum_cells, x_assemble_packed, x_count8, x_presence8, x_rank,
    x_rank_counts8, x_sum_packed, Csr,
};
use super::x_mb_lookup_par::Words;
use super::{build_unmasked_ranges, PvArrayType, BLAST2NA_MASK};
use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// No NCBI counterpart: reads LOSAT_X_BNCSR once (a switch of the strict set: on unless
/// `LOSAT_X_BNCSR=0`). It chooses between two LOSAT layouts of the same chains and changes no
/// value NCBI computes.
pub(super) fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNCSR"))
}

/// No NCBI counterpart: EXPERIMENT (LOSAT_X_BNCSRPACK, in the strict set: on unless it is `0`) the
/// head of a CSR table whose chains list their words newest first (every table of `plan` but the
/// reading-order one of 4^8 cells) is packed: two words per 64 cells (their presence bits and the
/// number of non-empty cells before them) and one offset per non-empty cell and one more, instead
/// of one offset per cell. A table of 4^11 cells has 4,194,305 offsets (16 MiB), of which a query
/// batch of about 1 Mb of `-task blastn` fills 1.2-1.6 M; the packed head is 1 MiB of bits and
/// ranks and 4.6-6.1 MiB of offsets. The build counts the words into one byte per cell (4 MiB, in
/// the capacity of the run list before the place pass writes it; `x_count8`), takes the presence
/// bits and ranks from the bytes and sums the non-empty cells instead of 16 MiB, with two passes
/// over the words as the unpacked build. A scan reads a word's bits (1 MiB) and, for a non-empty
/// cell, its two offsets; a word of an empty cell, whose presence-vector bit is set for 8 cells of
/// a 4^11 table, is dropped there instead of reading the 16 MiB array. The runs, their order and
/// the presence vector and helper counters are those of the unpacked head (`x_sum_packed`), so
/// every value read from the table is the same.
pub(super) fn x_pack_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNCSRPACK"))
}

/// No NCBI counterpart: whether a CSR of this order will be packed (LOSAT_X_BNCSRPACK; newest
/// first only), without a reach record.
pub(super) fn x_will_pack(ascending_cells: bool) -> bool {
    !ascending_cells && x_pack_on()
}

/// No NCBI counterpart: whether a CSR of this order is packed (LOSAT_X_BNCSRPACK; newest first only).
pub(super) fn x_packed(ascending_cells: bool) -> bool {
    let packed = x_will_pack(ascending_cells);
    crate::utils::xreach::x_reach!("LOSAT_X_BNCSRPACK", packed);
    packed
}

// No NCBI counterpart: counters of the summary line.
static BUILDS: AtomicU64 = AtomicU64::new(0);
static BUILDS_ON_POOL: AtomicU64 = AtomicU64::new(0);
static BUILT_WORDS: AtomicU64 = AtomicU64::new(0);

/// The `[TIMING]` line (with LOSAT_TIMING, when a CSR was built), printed by `main` at exit
/// (through `x_mb_lookup_print_summary`).
// No NCBI counterpart: prints counters; it does not change any value NCBI computes.
pub(super) fn print_summary() {
    let load = |c: &AtomicU64| c.load(Ordering::Relaxed);
    if std::env::var_os("LOSAT_TIMING").is_some() && load(&BUILDS) > 0 {
        eprintln!(
            "[TIMING] x_bn_csr: builds={} parallel={} words={}",
            load(&BUILDS),
            load(&BUILDS_ON_POOL),
            load(&BUILT_WORDS)
        );
    }
}

/// No NCBI counterpart: the `lut_word_length` range of the tables laid out as CSR (module
/// comment, scope).
const MIN_LUT: usize = 8;
const MAX_LUT: usize = 11;

/// No NCBI counterpart: EXPERIMENT (LOSAT_X_MBCSR, in the strict set: on unless it is `0`) megablast's
/// tables of 4^12 cells are laid out as CSR too, with a packed head, when the search scans at least
/// `X_MB12_LETTERS_PER_POSITION` subject letters per position of the word loop. It chooses between
/// two LOSAT layouts of the same chains and changes no value NCBI computes.
pub(super) fn x_mb12_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_MBCSR"))
}

/// No NCBI counterpart: LOSAT_X_MBCSR's size rule. At one thread on the workstation the CSR's build
/// costs 13-22 ns more per position of the word loop than the chains' fill, and its scan saves
/// about 3.2 ns per subject letter (draft200 x one 4 Gnt volume of ref_prok_rep_genomes, x0.72;
/// 1.5 Gnt of plasmids, x0.74), so the two meet at about 6 letters per position.
const X_MB12_LETTERS_PER_POSITION: u64 = 8;

/// No NCBI counterpart: whether a table of `lut_word_length` is in the scope (`MIN_LUT`, `MAX_LUT`;
/// LOSAT_X_MBCSR: 12 too, when its head is packed and the search scans enough letters).
/// `scan_letters`: the letters of the search's subjects, `None` when the caller does not say.
fn x_in_scope(words: &Words<'_>, scan_letters: Option<u64>) -> bool {
    if (MIN_LUT..=MAX_LUT).contains(&words.lut_word_length) {
        return true;
    }
    if words.lut_word_length != 12 || !x_will_pack(words.ascending_cells) || !x_mb12_on() {
        return false;
    }
    let taken = scan_letters.is_some_and(|letters| {
        letters >= (positions(words) as u64).saturating_mul(X_MB12_LETTERS_PER_POSITION)
    });
    crate::utils::xreach::x_reach!("LOSAT_X_MBCSR", taken);
    taken
}

/// No NCBI counterpart: whether `build_mb_lookup` lays out a table of `cells` cells over `words`
/// as CSR (module comment, scope). `db_word_counts`: the table is built with database word counts;
/// `debug`: BLEMIR_DEBUG; `scan_letters`: as in `x_in_scope`.
pub(super) fn plan(
    words: &Words<'_>,
    cells: usize,
    db_word_counts: bool,
    debug: bool,
    scan_letters: Option<u64>,
) -> bool {
    if !on() {
        return false;
    }
    // EXPERIMENT (LOSAT_X_BNCSRPACK): a packed head has no size rule (module comment).
    let csr = !db_word_counts
        && !debug
        && x_in_scope(words, scan_letters)
        && (x_will_pack(words.ascending_cells) || wanted(positions(words), cells));
    crate::utils::xreach::x_reach!("LOSAT_X_BNCSR", csr);
    csr
}

/// No NCBI counterpart: the size rule (module comment): at least one position of the word loop
/// per four cells.
fn wanted(positions: usize, cells: usize) -> bool {
    positions.saturating_mul(4) >= cells
}

/// No NCBI counterpart: the positions the word loop reads (the unmasked ranges of `word_length`
/// positions or more of the queries of `lut_word_length` letters or more; every word ends at one
/// of them).
fn positions(words: &Words<'_>) -> usize {
    let mut total = 0usize;
    for (q_idx, seq) in words.queries_blastna.iter().enumerate() {
        if seq.len() < words.lut_word_length {
            continue;
        }
        let masks = words
            .query_masks
            .get(q_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);
        for (start, end) in build_unmasked_ranges(seq.len(), masks) {
            let len = end.saturating_sub(start);
            if words.word_length <= len {
                total += len;
            }
        }
    }
    total
}

/// The word loop of `build_mb_lookup` (`x_mb_lookup_par::fill_part` without the owner test and
/// the update): calls `f(cell, q_off_1)` for each word, in the loop's order.
#[inline(always)]
fn each_word(words: &Words<'_>, mut f: impl FnMut(usize, u32)) {
    let lut_word_length = words.lut_word_length;
    for (q_idx, seq) in words.queries_blastna.iter().enumerate() {
        if seq.len() < lut_word_length {
            continue;
        }
        let masks = words
            .query_masks
            .get(q_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);
        let query_offset = words.query_offsets[q_idx].max(0) as usize;
        for (range_start, range_end) in build_unmasked_ranges(seq.len(), masks) {
            if words.word_length > range_end.saturating_sub(range_start) {
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
                f(
                    current_kmer as usize,
                    (query_offset + (pos + 1 - lut_word_length) + 1) as u32,
                );
            }
        }
    }
}

/// No NCBI counterpart: the CSR of a table that `plan` lays out so, with the presence bits and
/// helper counters into `pv_array` and `helper_array` (zeroed, of the fill's lengths): on the
/// search pool where the fill would run (`x_mb_lookup_once::csr`), otherwise `serial`. Returns it
/// and whether the pool built it.
pub(super) fn build(
    words: &Words<'_>,
    cells: usize,
    pv_array: &mut [PvArrayType],
    pv_array_bts: usize,
    helper_array: &mut [u32],
    x_offsets: Vec<u32>,
) -> (Csr, bool) {
    let (csr, on_pool) = match super::x_mb_lookup_par::pool_threads() {
        Some(parts) => super::x_mb_lookup_once::csr(
            parts,
            words,
            cells,
            pv_array,
            pv_array_bts,
            helper_array,
            x_offsets,
        ),
        None => (
            serial(
                words,
                cells,
                pv_array,
                pv_array_bts,
                helper_array,
                x_offsets,
            ),
            false,
        ),
    };
    BUILDS.fetch_add(1, Ordering::Relaxed);
    BUILDS_ON_POOL.fetch_add(u64::from(on_pool), Ordering::Relaxed);
    BUILT_WORDS.fetch_add(csr.words() as u64, Ordering::Relaxed);
    (csr, on_pool)
}

/// No NCBI counterpart: the element of `offsets` that holds the count of cell 0 (module comment:
/// `offsets[cell]` newest first, `offsets[cell + 2]` in reading order).
pub(super) fn count_shift(ascending_cells: bool) -> usize {
    if ascending_cells {
        2
    } else {
        0
    }
}

/// The serial build (module comment): the count pass (each word's count prefetched and counted
/// 16 words later), the sum pass with the presence bits and helper counters into `pv_array` and
/// `helper_array` (zeroed, of the fill's lengths), and the place pass.
pub(super) fn serial(
    words: &Words<'_>,
    cells: usize,
    pv_array: &mut [PvArrayType],
    pv_array_bts: usize,
    helper_array: &mut [u32],
    x_offsets: Vec<u32>,
) -> Csr {
    // Once per table: the bounds the unchecked writes below rely on.
    assert!(cells % 2048 == 0 && pv_array_bts >= 5);
    assert_eq!(helper_array.len(), cells >> 11);
    let ascending = words.ascending_cells;
    if x_packed(ascending) {
        return x_serial_packed(
            words,
            cells,
            pv_array,
            pv_array_bts,
            helper_array,
            x_offsets,
        );
    }
    x_serial_unpacked(
        words,
        cells,
        pv_array,
        pv_array_bts,
        helper_array,
        x_offsets,
    )
}

/// `serial` with one offset per cell (either order).
pub(super) fn x_serial_unpacked(
    words: &Words<'_>,
    cells: usize,
    pv_array: &mut [PvArrayType],
    pv_array_bts: usize,
    helper_array: &mut [u32],
    x_offsets: Vec<u32>,
) -> Csr {
    assert!(cells % 2048 == 0 && pv_array_bts >= 5);
    assert_eq!(helper_array.len(), cells >> 11);
    let ascending = words.ascending_cells;
    // EXPERIMENT (LOSAT_X_MBLUTKEEP): `x_offsets` is the batch loop's kept table or an empty vector.
    let mut offsets = super::x_dc_csr::x_offsets_buffer(x_offsets, cells + 2);
    let ob = offsets.as_mut_ptr();
    // SAFETY: `cells + 2` elements of the capacity, written here before anything reads them.
    unsafe { std::ptr::write_bytes(ob, 0, cells + 2) };
    let cnt = ob.wrapping_add(count_shift(ascending));
    const AHEAD: usize = 16;
    let mut ring = [0usize; AHEAD];
    let mut seen = 0usize;
    each_word(words, |cell, _| {
        assert!(cell < cells);
        prefetch(cnt.wrapping_add(cell));
        let slot = seen % AHEAD;
        if seen >= AHEAD {
            // SAFETY: `ring[slot] < cells` (asserted when it entered); `cnt` has `cells` elements.
            unsafe { *cnt.add(ring[slot]) += 1 };
        }
        ring[slot] = cell;
        seen += 1;
    });
    for i in seen.saturating_sub(AHEAD)..seen {
        // SAFETY: as above.
        unsafe { *cnt.add(ring[i % AHEAD]) += 1 };
    }
    let n_words = seen;
    assert!(
        n_words <= u32::MAX as usize,
        "BNCSR: more words than u32 offsets"
    );
    // SAFETY: the counts of all cells are this build's; the presence vector and the counters are
    // of the fill's lengths (`c >> 11 < helper_array.len()` for every cell).
    let total = unsafe {
        sum_cells(
            pv_array.as_mut_ptr(),
            pv_array.len(),
            helper_array.as_mut_ptr(),
            cnt,
            0,
            0,
            cells,
            pv_array_bts,
        )
    };
    debug_assert_eq!(total as usize, n_words);
    let mut qs: Vec<u32> = Vec::with_capacity(n_words);
    let qb = qs.as_mut_ptr();
    if ascending {
        each_word(words, |cell, q_off_1| {
            // SAFETY: `cell < cells` (the same words as the count pass); `offsets[cell + 1]`
            // walks the cell's run from its start, one place per word of the cell, so the place
            // is below `n_words`, in the capacity.
            unsafe {
                let cursor = ob.add(cell + 1);
                *qb.add(*cursor as usize) = q_off_1;
                *cursor += 1;
            }
        });
    } else {
        // SAFETY: `cells < cells + 2`; `offsets[cells]` is no cell's end.
        unsafe { *ob.add(cells) = total };
        each_word(words, |cell, q_off_1| {
            // SAFETY: `cell < cells` (the same words as the count pass), `offsets[cell]` is the end
            // of the cell's non-empty run (>= 1), so the place `e < n_words` is in the capacity.
            unsafe {
                let e = *ob.add(cell) - 1;
                *ob.add(cell) = e;
                *qb.add(e as usize) = q_off_1;
            }
        });
    }
    // SAFETY: `offsets[..cells + 1]` was written in full; the runs partition `0 .. n_words`, and
    // the place pass wrote each place once.
    unsafe { assemble(offsets, qs, cells, n_words) }
}

/// EXPERIMENT (LOSAT_X_BNCSRPACK): `serial` with a packed head (newest first): the count pass into
/// one byte per cell (`x_count8`, each word's byte prefetched and counted 16 words later), the
/// presence bits of the non-zero counts (`x_presence8`), the ranks of the blocks and the counts by
/// rank (`x_rank_counts8`; the bytes are in `qs`'s capacity, which the place pass writes then), the
/// sum pass (`x_sum_packed`, with the presence bits and helper counters), and the place pass, as
/// `x_serial_unpacked` does on the offsets of all cells: two passes over the words, as there.
pub(super) fn x_serial_packed(
    words: &Words<'_>,
    cells: usize,
    pv_array: &mut [PvArrayType],
    pv_array_bts: usize,
    helper_array: &mut [u32],
    x_offsets: Vec<u32>,
) -> Csr {
    assert!(cells % 2048 == 0 && pv_array_bts >= 5 && !words.ascending_cells);
    assert_eq!(helper_array.len(), cells >> 11);
    let blocks = cells >> 6;
    // `qs` has room for every position the loop reads (each word ends at one) and for the counts,
    // which take its first `cells` bytes, so the build frees nothing in the middle of a batch (a
    // freed block there is trimmed and faulted in again; with threads each trim interrupts them
    // all).
    let mut qs: Vec<u32> = Vec::with_capacity(positions(words).max(cells.div_ceil(4)));
    let cb = qs.as_mut_ptr() as *mut u8;
    // SAFETY: `cells` bytes of the capacity, before the place pass writes any element.
    unsafe { std::ptr::write_bytes(cb, 0, cells) };
    let mut more = HashMap::new();
    const AHEAD: usize = 16;
    let mut ring = [0usize; AHEAD];
    let mut seen = 0usize;
    each_word(words, |cell, _| {
        debug_assert!(cell < cells);
        prefetch(cb.wrapping_add(cell));
        let slot = seen % AHEAD;
        if seen >= AHEAD {
            // SAFETY: `ring[slot] < cells`, a byte of the counts; one thread.
            unsafe { x_count8(cb, &mut more, ring[slot]) };
        }
        ring[slot] = cell;
        seen += 1;
    });
    for i in seen.saturating_sub(AHEAD)..seen {
        // SAFETY: as above.
        unsafe { x_count8(cb, &mut more, ring[i % AHEAD]) };
    }
    let mut bits = vec![0u64; 2 * blocks];
    // SAFETY: `blocks` blocks of 64 cells: the `cells` bytes of the counts, the `2 blocks` words of
    // `bits`.
    let present = unsafe { x_presence8(bits.as_mut_ptr(), cb, 0, blocks) };
    // EXPERIMENT (LOSAT_X_MBLUTKEEP): `x_offsets` is the batch loop's kept table or an empty vector.
    let mut starts = super::x_dc_csr::x_offsets_buffer(x_offsets, present + 1);
    let sb = starts.as_mut_ptr();
    // SAFETY: as above, and the ranks `0 .. present` of the capacity's `present + 1` elements.
    let ranked = unsafe { x_rank_counts8(bits.as_mut_ptr(), cb, &more, sb, 0, blocks, 0) };
    debug_assert_eq!(ranked, present);
    drop(more);
    let n_words = seen;
    assert!(
        n_words <= u32::MAX as usize,
        "BNCSR: more words than u32 offsets"
    );
    // SAFETY: the counts of all non-empty cells are this build's; the presence vector and the
    // counters are of the fill's lengths.
    let total = unsafe {
        x_sum_packed(
            pv_array.as_mut_ptr(),
            pv_array.len(),
            helper_array.as_mut_ptr(),
            &bits,
            sb,
            0,
            0,
            blocks,
            pv_array_bts,
        )
    };
    debug_assert_eq!(total as usize, n_words);
    debug_assert!(qs.capacity() >= n_words);
    let qb = qs.as_mut_ptr();
    // SAFETY: `present < present + 1`; `starts[present]` is no cell's end.
    unsafe { *sb.add(present) = total };
    each_word(words, |cell, q_off_1| {
        let r = x_rank(&bits, cell).expect("BNCSRPACK: a word of an empty cell");
        // SAFETY: `r < present`, `starts[r]` is the end of the cell's non-empty run (>= 1), so the
        // place `e < n_words` is in the capacity.
        unsafe {
            let e = *sb.add(r) - 1;
            *sb.add(r) = e;
            *qb.add(e as usize) = q_off_1;
        }
    });
    // SAFETY: `starts[..present + 1]` was written in full; the runs partition `0 .. n_words`, and
    // the place pass wrote each place once.
    unsafe { x_assemble_packed(starts, bits, qs, present, n_words) }
}

/// No NCBI counterpart: whether `run` (a cell's run, sorted: increasing in reading order,
/// decreasing newest first; the offsets of a cell are distinct) holds `q_off_1`.
#[inline]
pub(super) fn contains(run: &[u32], q_off_1: u32) -> bool {
    match run {
        [first, .., last] if first < last => run.binary_search(&q_off_1).is_ok(),
        _ => run.binary_search_by(|q| q_off_1.cmp(q)).is_ok(),
    }
}

/// No NCBI counterpart: the check of the unit tests: `csr` (its presence vector and counters in
/// `pv` and `helper`) against the chains, presence vector and counters of the fill (`hashtable`,
/// `next_pos`, `fill_pv`, `fill_helper`) and their longest chain, and `contains` on each run.
/// Panics on a difference.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
pub(super) fn x_assert_runs_are_the_chains(
    csr: &Csr,
    pv: &[PvArrayType],
    helper: &[u32],
    hashtable: &[u32],
    next_pos: &[u32],
    fill_pv: &[PvArrayType],
    fill_helper: &[u32],
    what: &str,
) {
    let mut diffs = Vec::new();
    let mut chain = Vec::new();
    let mut words = 0usize;
    for (cell, &head) in hashtable.iter().enumerate() {
        chain.clear();
        let mut q = head;
        while q != 0 {
            chain.push(q);
            q = next_pos[q as usize];
        }
        words += chain.len();
        let run = csr.x_run_any(cell);
        if run != chain.as_slice() {
            diffs.push(format!(
                "cell {cell} run {:?} chain {:?}",
                &run[..run.len().min(8)],
                &chain[..chain.len().min(8)]
            ));
            break;
        }
        let mut sorted = chain.clone();
        sorted.sort_unstable();
        for &q in &chain {
            if !contains(run, q) {
                diffs.push(format!("cell {cell} does not contain {q}"));
            }
            let other = q + 1;
            if sorted.binary_search(&other).is_err() && contains(run, other) {
                diffs.push(format!("cell {cell} contains {other}"));
            }
        }
        if contains(run, 0) {
            diffs.push(format!("cell {cell} contains 0"));
        }
        if !diffs.is_empty() {
            break;
        }
    }
    if csr.words() != words {
        diffs.push(format!("{} words, {words} in the chains", csr.words()));
    }
    if pv != fill_pv {
        diffs.push("presence vector".to_string());
    }
    if helper != fill_helper {
        diffs.push("helper counters".to_string());
    }
    let longest = |h: &[u32]| h.iter().fold(2usize, |m, &v| m.max(v as usize));
    if longest(helper) != longest(fill_helper) {
        diffs.push("longest chain".to_string());
    }
    assert!(
        diffs.is_empty(),
        "{what}: the CSR differs from the chains (cells {}, words {}): {}",
        hashtable.len(),
        csr.words(),
        diffs.join("; ")
    );
}

#[cfg(test)]
pub(super) mod tests {
    use super::super::x_mb_lookup_par::tests::{fill_old, Lcg};
    use super::super::x_mb_lookup_par::Tables;
    use super::super::{build_two_stage_lookup, compute_mb_pv_params, MaskedInterval};
    use super::*;

    /// Queries (BLASTNA), their offsets in the concatenation and their masks.
    pub(crate) struct Case {
        pub(crate) queries: Vec<Vec<u8>>,
        pub(crate) offsets: Vec<i32>,
        pub(crate) masks: Vec<Vec<MaskedInterval>>,
    }

    impl Case {
        pub(crate) fn refs(&self) -> Vec<&[u8]> {
            self.queries.iter().map(Vec::as_slice).collect()
        }

        pub(crate) fn words<'a>(
            &'a self,
            refs: &'a [&'a [u8]],
            word_length: usize,
            lut_word_length: usize,
            ascending_cells: bool,
        ) -> Words<'a> {
            Words {
                queries_blastna: refs,
                query_offsets: &self.offsets,
                query_masks: &self.masks,
                word_length,
                lut_word_length,
                kmer_mask: (1u64 << (2 * lut_word_length)) - 1,
                ascending_cells,
            }
        }
    }

    /// EXPERIMENT (LOSAT_X_BNCSRPACK): `case` with four more queries whose words fill a few cells
    /// with more words than a byte counts: 1,200 letters A (one cell), 800 letters ACAC... (two
    /// cells), and one letter repeated so that its cell gets exactly 255 (T) or 256 (G) words of a
    /// table of `lut` letters.
    pub(crate) fn x_with_repeats(case: &mut Case, lut: usize) {
        let mut at = match (case.offsets.last(), case.queries.last()) {
            (Some(&o), Some(q)) => o + q.len() as i32 + 1,
            _ => 0,
        };
        for seq in [
            vec![0u8; 1200],
            (0..800).map(|i| (i % 2) as u8).collect(),
            vec![3u8; 255 + lut - 1],
            vec![2u8; 256 + lut - 1],
        ] {
            case.offsets.push(at);
            at += seq.len() as i32 + 1;
            case.queries.push(seq);
            case.masks.push(Vec::new());
        }
    }

    /// `n` random queries of about `len` letters: few distinct letters in half of them (so cells
    /// get long chains), ambiguities (N), masked intervals, and an empty query, a query of N only
    /// and a query of 5 letters (shorter than every table's words) among them.
    pub(crate) fn random_case(rng: &mut Lcg, n: usize, len: usize) -> Case {
        let mut queries = Vec::new();
        let mut offsets = Vec::new();
        let mut masks = Vec::new();
        let mut at = 0i32;
        for q in 0..n + 3 {
            let seq: Vec<u8> = match q {
                0 => Vec::new(),
                1 => vec![14; 40],
                2 => vec![0, 1, 2, 3, 0],
                _ => {
                    let l = len / 2 + rng.below(len as u64 + 1) as usize;
                    let letters = if q % 2 == 0 { 2 } else { 4 };
                    (0..l)
                        .map(|_| {
                            if rng.below(150) == 0 {
                                14
                            } else {
                                rng.below(letters) as u8
                            }
                        })
                        .collect()
                }
            };
            let mut m = Vec::new();
            if q % 3 == 1 && seq.len() > 60 {
                let s = rng.below((seq.len() - 40) as u64) as usize;
                m.push(MaskedInterval {
                    start: s,
                    end: s + 1 + rng.below(30) as usize,
                });
            }
            offsets.push(at);
            at += seq.len() as i32 + 1;
            queries.push(seq);
            masks.push(m);
        }
        Case {
            queries,
            offsets,
            masks,
        }
    }

    /// The table of `cells` cells as the loop fills it (`fill_old` with one part: the loop's
    /// statements, re-threaded for `ascending_cells`): `(hashtable, next_pos, pv_array,
    /// helper_array)`.
    pub(crate) fn chains(
        case: &Case,
        words: &Words<'_>,
        cells: usize,
        pv_len: usize,
        pv_array_bts: usize,
    ) -> (Vec<u32>, Vec<u32>, Vec<PvArrayType>, Vec<u32>) {
        let query_length = case
            .queries
            .iter()
            .zip(&case.offsets)
            .map(|(seq, &off)| off.max(0) as usize + seq.len())
            .max()
            .unwrap_or(0);
        let mut hashtable = vec![0u32; cells];
        let mut next_pos = vec![0u32; query_length + 1];
        let mut pv = vec![0 as PvArrayType; pv_len];
        let mut helper = vec![0u32; (cells >> 11).max(1)];
        fill_old(
            1,
            words,
            Tables {
                hashtable: &mut hashtable,
                next_pos: &mut next_pos,
                pv_array: &mut pv,
                pv_array_bts,
                helper_array: &mut helper,
                x_h24: false,
            },
        );
        (hashtable, next_pos, pv, helper)
    }

    /// The presence-vector length and shift of a table of `lut_word_length` for `case`.
    pub(crate) fn pv_params(case: &Case, lut_word_length: usize) -> (usize, usize) {
        let approx = case.queries.iter().map(Vec::len).sum();
        compute_mb_pv_params(1 << (2 * lut_word_length), approx, false, lut_word_length)
    }

    // The serial build lists the loop's chains, newest first and in reading order, for every
    // `lut_word_length` of the rule, with the fill's presence vector, counters and longest chain
    // (`x_assert_runs_are_the_chains`: every cell's run, `contains` on each run).
    #[test]
    fn the_runs_are_the_chains_in_both_orders() {
        let mut rng = Lcg(0xB4C5_0001);
        for lut in MIN_LUT..=MAX_LUT {
            for ascending in [false, true] {
                for (n, len) in [(1, 6), (3, 300), (5, 6000)] {
                    let case = random_case(&mut rng, n, len);
                    let refs = case.refs();
                    let word_length = lut + rng.below(6) as usize;
                    let words = case.words(&refs, word_length, lut, ascending);
                    let cells = 1usize << (2 * lut);
                    let (pv_len, bts) = pv_params(&case, lut);
                    let (hashtable, next_pos, fill_pv, fill_helper) =
                        chains(&case, &words, cells, pv_len, bts);
                    let mut pv = vec![0 as PvArrayType; pv_len];
                    let mut helper = vec![0u32; cells >> 11];
                    let csr = serial(&words, cells, &mut pv, bts, &mut helper, Vec::new());
                    x_assert_runs_are_the_chains(
                        &csr,
                        &pv,
                        &helper,
                        &hashtable,
                        &next_pos,
                        &fill_pv,
                        &fill_helper,
                        &format!("lut {lut} ascending {ascending} n {n} len {len}"),
                    );
                }
            }
        }
    }

    // LOSAT_X_BNCSRPACK: the packed head and the head of one offset per cell list the loop's chains
    // (newest first), with the fill's presence vector, counters and longest chain, whatever the
    // switch reads; many queries make tables with most cells non-empty.
    #[test]
    fn the_packed_head_lists_the_chains_too() {
        let mut rng = Lcg(0xB4C5_0002);
        for lut in MIN_LUT..=MAX_LUT {
            for (n, len, repeats) in [
                (1, 6, false),
                (3, 300, false),
                (5, 6000, false),
                (40, 3000, false),
                (2, 500, true),
            ] {
                let mut case = random_case(&mut rng, n, len);
                if repeats {
                    x_with_repeats(&mut case, lut);
                }
                let refs = case.refs();
                let word_length = lut + rng.below(6) as usize;
                let words = case.words(&refs, word_length, lut, false);
                let cells = 1usize << (2 * lut);
                let (pv_len, bts) = pv_params(&case, lut);
                let (hashtable, next_pos, fill_pv, fill_helper) =
                    chains(&case, &words, cells, pv_len, bts);
                for packed in [false, true] {
                    let mut pv = vec![0 as PvArrayType; pv_len];
                    let mut helper = vec![0u32; cells >> 11];
                    let build = if packed {
                        x_serial_packed
                    } else {
                        x_serial_unpacked
                    };
                    let csr = build(&words, cells, &mut pv, bts, &mut helper, Vec::new());
                    let what = format!("lut {lut} packed {packed} n {n} len {len}");
                    if repeats {
                        // the cells of 255, 256 and more words
                        let mut runs: Vec<usize> =
                            (0..cells).map(|c| csr.x_run_any(c).len()).collect();
                        runs.sort_unstable();
                        assert!(runs.contains(&255) && runs.contains(&256), "{what}");
                        assert!(runs[cells - 1] >= 1200 - lut, "{what}");
                    }
                    x_assert_runs_are_the_chains(
                        &csr,
                        &pv,
                        &helper,
                        &hashtable,
                        &next_pos,
                        &fill_pv,
                        &fill_helper,
                        &what,
                    );
                }
            }
        }
    }

    #[test]
    fn contains_searches_runs_of_either_order() {
        for (run, q, want) in [
            (&[][..], 1, false),
            (&[7][..], 7, true),
            (&[7][..], 8, false),
            (&[3, 9, 40][..], 9, true),
            (&[3, 9, 40][..], 3, true),
            (&[3, 9, 40][..], 40, true),
            (&[3, 9, 40][..], 10, false),
            (&[40, 9, 3][..], 9, true),
            (&[40, 9, 3][..], 40, true),
            (&[40, 9, 3][..], 3, true),
            (&[40, 9, 3][..], 4, false),
            (&[40, 9][..], 0, false),
        ] {
            assert_eq!(contains(run, q), want, "{run:?} {q}");
        }
    }

    #[test]
    fn the_rule_takes_a_quarter_of_the_cells_and_lut_8_to_11() {
        assert!(!wanted(0, 1 << 22));
        assert!(!wanted((1 << 20) - 1, 1 << 22));
        assert!(wanted(1 << 20, 1 << 22));
        assert!(wanted(usize::MAX, 1 << 22));
        // The positions are those the loop reads: masked ranges, queries shorter than
        // `lut_word_length` and ranges shorter than `word_length` do not count.
        let case = Case {
            queries: vec![vec![0; 100], vec![1; 7], vec![2; 30]],
            offsets: vec![0, 101, 109],
            masks: vec![
                vec![MaskedInterval { start: 10, end: 20 }],
                Vec::new(),
                vec![MaskedInterval { start: 12, end: 13 }],
            ],
        };
        let refs = case.refs();
        // ranges 0..10, 20..100 (both >= 9), the 7-letter query (< lut 8), 0..12 and 13..30
        assert_eq!(
            positions(&case.words(&refs, 9, 8, false)),
            10 + 80 + 12 + 17
        );
        // word size 12: the range of 10 letters is dropped, the one of 12 is not
        assert_eq!(positions(&case.words(&refs, 12, 8, false)), 80 + 12 + 17);
        if on() {
            let big = Case {
                queries: vec![vec![0; 20_000]],
                offsets: vec![0],
                masks: vec![Vec::new()],
            };
            let refs = big.refs();
            assert!(plan(
                &big.words(&refs, 9, 8, true),
                1 << 16,
                false,
                false,
                None
            ));
            assert!(!plan(
                &big.words(&refs, 9, 8, true),
                1 << 16,
                true,
                false,
                None
            ));
            assert!(!plan(
                &big.words(&refs, 9, 8, true),
                1 << 16,
                false,
                true,
                None
            ));
            assert!(!plan(
                &big.words(&refs, 9, 7, true),
                1 << 14,
                false,
                false,
                None
            ));
            assert!(!plan(
                &big.words(&refs, 13, 12, false),
                1 << 24,
                false,
                false,
                None
            ));
            // LOSAT_X_MBCSR: 4^12 cells with a search of at least 8 letters per position.
            let positions = big.refs()[0].len() as u64;
            assert_eq!(
                plan(
                    &big.words(&refs, 13, 12, false),
                    1 << 24,
                    false,
                    false,
                    Some(positions * X_MB12_LETTERS_PER_POSITION)
                ),
                x_will_pack(false) && x_mb12_on()
            );
            assert!(!plan(
                &big.words(&refs, 13, 12, false),
                1 << 24,
                false,
                false,
                Some(positions * X_MB12_LETTERS_PER_POSITION - 1)
            ));
            // below the rule: chains, unless the head is packed (LOSAT_X_BNCSRPACK has no rule)
            assert_eq!(
                plan(
                    &big.words(&refs, 11, 11, false),
                    1 << 22,
                    false,
                    false,
                    None
                ),
                x_will_pack(false)
            );
        }
    }

    // A table above the rule, built by `build_two_stage_lookup` off and on a pool, is laid out as
    // CSR (with the switch on) and answers `for_each_hit`, `count_hits`, `contains_hit` and
    // `longest_chain` as the loop's chains, in both orders.
    #[test]
    fn a_table_above_the_rule_answers_from_its_csr() {
        let mut rng = Lcg(0xB4C5_0002);
        // lut 8 in reading order (16,384 positions for 4^8 cells) and lut 9 newest first (65,536
        // for 4^9)
        for (lut, ascending, len) in [(8usize, true, 12_000usize), (9, false, 40_000)] {
            let case = random_case(&mut rng, 3, len);
            let refs = case.refs();
            let word_length = lut + 1;
            let words = case.words(&refs, word_length, lut, ascending);
            assert!(wanted(positions(&words), 1 << (2 * lut)));
            let cells = 1usize << (2 * lut);
            let (pv_len, bts) = pv_params(&case, lut);
            let (hashtable, next_pos, _, fill_helper) = chains(&case, &words, cells, pv_len, bts);
            let approx = case.queries.iter().map(Vec::len).sum();
            let check = |what: &str| {
                let lookup = build_two_stage_lookup(
                    &refs,
                    &case.offsets,
                    word_length,
                    lut,
                    &case.masks,
                    None,
                    0,
                    approx,
                    ascending,
                );
                assert_eq!(lookup.x_mb_csr().is_some(), on(), "{what}");
                let longest = fill_helper.iter().fold(2usize, |m, &v| m.max(v as usize));
                assert_eq!(lookup.longest_chain(), longest, "{what}");
                let mut chain = Vec::new();
                for (cell, &head) in hashtable.iter().enumerate() {
                    chain.clear();
                    let mut q = head;
                    while q != 0 {
                        chain.push(q);
                        q = next_pos[q as usize];
                    }
                    let mut seen = Vec::new();
                    lookup.for_each_hit(cell as u64, |q| seen.push(q));
                    assert_eq!(seen, chain, "{what} cell {cell}");
                    assert_eq!(lookup.count_hits(cell as u64), chain.len(), "{what}");
                    for &q in chain.iter().take(3) {
                        assert!(lookup.contains_hit(cell as u64, q), "{what} cell {cell}");
                    }
                    assert!(!lookup.contains_hit(cell as u64, 1 << 30), "{what}");
                    if let Some(csr) = lookup.x_mb_csr() {
                        assert_eq!(csr.x_run_any(cell), chain.as_slice(), "{what} cell {cell}");
                    }
                }
            };
            check(&format!("lut {lut} off a pool"));
            #[cfg(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            ))]
            crate::utils::threading::with_search_pool(4, "blastn", |_pool| {
                check(&format!("lut {lut} on a pool of 4"));
                Ok(())
            })
            .unwrap();
        }
    }
}
