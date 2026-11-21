//! The BLASTX seed loop reads the scan buffer directly, and the buffer is kept per thread.
//!
//! NCBI allocates the offset-pair array once per search thread and its word finder reads the
//! pairs straight from that array:
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1040-1041
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-511
//!
//! The reference (`Lookup::scan` + `word_finder`, which the tests compare with) allocates and
//! zero-fills a buffer of `capacity()` pairs per (chunk, subject) and copies every batch into a
//! new `(i32, i32)` vector before the diagonal loop reads it. `x_word_finder_scan` is
//! `word_finder` with the scan loop of `Lookup::scan` inlined: the same `scan_range` state
//! machine, the same `array_size` (`capacity()`), the scanner writing into a per-thread buffer
//! (grown when a table needs more, never shrunk), and the diagonal statements reading the pairs in
//! place.
//!
//! Identity argument: the scanners (`s_blast_aa_scan_subject_one_range`, the compressed
//! `scan_subject`) only write the buffer and return the number of pairs written, and only
//! `pairs[0..hits)` is read afterwards, so stale contents beyond `hits` are never read. The
//! pairs are visited in the same order with the same values (`u32 as i32`, as the copy did),
//! and the diagonal statements below are those of `word_finder`. The searches take this path
//! for every subject but a `-subject` with scan ranges (`x_soft_subject.rs`).
//!
//! BXSCANBF (graduated 2026-10-26: no variable, always on): the scan of the standard
//! lookup table is tblastn's branch-free scan (`tblastx/x_aascan_bf.rs`, `x_scan_one_range_bf`,
//! the loop of the former TBNSCANBF): the presence bits of a block of words advance a candidate list
//! instead of a branch, and the listed words' hits are copied in position order with NCBI's stop
//! and resume when the offset array is full. Its pairs and `scan_range` values are those of
//! `s_blast_aa_scan_subject_one_range`, so the loop above reads the same pairs; it writes up to `X_SCAN_BF_SLACK` pairs past `array_size`, which the
//! per-thread buffer then holds and nothing reads. With LOSAT_X_BXDIAG the fused loop takes the
//! same candidate loop (`x_seed_diag.rs`). The compressed lookup table keeps its scan.
//!
//! EXPERIMENT (LOSAT_X_BXHITBF): the direct loop takes the two outcomes most two-hit tests end in,
//! too far and overlap, without a branch (tblastn's branch-free two-hit loop and TBLASTX's
//! LOSAT_X_TBXHITBF do the same in their loops):
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:538-551
//! Where the direct loop scans every subject (`LOSAT_X_BXDIAG=0`, and the default without a cell
//! test of batches across the cells' borders: `LOSAT_X_BXCROSS=0`, aarch64, x86-64 without SSSE3,
//! SSE4.1 or POPCNT, and `LOSAT_X_BXCROSS128=0` on x86-64 without AVX2 and in WebAssembly)
//! 79-81 % of the hits of the three blastx cases of v4-A10 probe bx1 are too
//! far and 14.5-17.4 % overlap, which of the two is hard to predict, and 1.2-1.4 % meet a flag-1
//! cell and 2.1-2.8 % are candidates for the extension.
//!
//! The test is blastp's (`blastp::x_seed_fused::XHitBf`, whose module header gives the arithmetic):
//! it reads the packed cell `raw` (`last_hit` in bits 0-30, the flag in bit 31) and works with
//! twice the values, `d2 = 2 (s + diag_offset) - (raw << 1)` in wrapping `u32` arithmetic read as
//! `i32`, which is twice the reference's `diff = s - (last_hit - diag_offset)` for a flag-0 cell
//! whose `last_hit` is in `[-window, diag_offset + length)` where `XHitBf::fits` holds (checked
//! once per subject): `0 <= diag_offset` and `diag_offset + length + 2 window < 2^30` (`length` the
//! subject's letters). This table holds no other flag-0 value: the cleared value, 0 in a new table,
//! or an `s + diag_offset` of this subject or of an earlier one. Every store with flag 0 is such an
//! `s + diag_offset` (the extension's `s_last_off - (wordsize - 1) + diag_offset` is stored with
//! flag 1, and the flag goes back to 0 only with a store of `s + diag_offset`); the offset only
//! grows between clears (`Diagonals::finish_subject`), so an earlier subject's stores are below
//! this subject's offset, and none of them wrapped in the 31-bit field, as they were below 2^30.
//! Then too far (`diff >= window`) is `d2 >= 2 window`, and a candidate (`wordsize <= diff <
//! window`) is `(d2 - 2 wordsize) as u32 < 2 (window - wordsize)`, one unsigned compare. The rule
//! also asks `1 <= wordsize < window`, a two-hit search (the one-hit search, window 0, keeps the
//! statements). A flag-1 cell (`raw as i32 < 0`) or a candidate goes on to the statements of
//! `word_finder`, unchanged; any other hit stores `s + diag_offset` (too far; the flag stays
//! 0) or what it read (overlap), selected without a branch. So the table holds after each hit what
//! the statements leave there, the statements see the same table for every hit they run, and the
//! init HSPs are the same, in the same order.
//!
//! The switch is read once; it is in the strict set (`utils/xswitch.rs`): on unless
//! `LOSAT_X_BXHITBF=0`. The loop with the test is its own
//! function (`x_word_finder_direct_bf`, never inlined), so the loop without it keeps its code. The
//! fused loop of LOSAT_X_BXDIAG keeps its own tests.
//!
//! EXPERIMENT (LOSAT_X_BXUNGAPSENT, `x_ungap_sent.rs`): for a subject on that path the direct
//! loop's two-hit extension is the sentinel walks, in instances of their own
//! (`x_word_finder_direct_sent`, never inlined); the loops without it keep their code.
//!
//! EXPERIMENT (LOSAT_X_BXSCANBYTE): the direct loop scans the standard table with tblastn's scan on
//! presence bytes (`tblastx/x_scan_byte.rs` `x_scan_one_range_byte`, LOSAT_X_TBNSCANBYTE's code,
//! whose header gives the exactness argument): the indexes of a block of 256 words made from the
//! letters at once, the candidates listed from one byte per index, the hits of a block copied with
//! one capacity test. It returns the pairs and `scan_range` values of
//! `s_blast_aa_scan_subject_one_range` on every call, so the loop above reads the same pairs; it
//! writes up to `X_SCAN_BYTE_SLACK` pairs past `array_size` (more than the branch-free scan, for
//! the chain runs of LOSAT_X_CHAINRUN), which the per-thread buffer holds and nothing reads. The presence bytes (32 KB) are built the first
//! time the direct loop scans with the table and kept with it (`XKeptPresence` in
//! `Lookup::Standard`), so they go with the table and no other table can see them. The switch
//! takes the CPU policy of LOSAT_X_TBNSCANBYTE (`x_mode`): on unless it is `0`, the AVX2 forms
//! where the CPU has AVX2, the plain loops for sparse tables (at most three tenths of the backbone
//! cells with hits) on x86-64 without AVX2 and in WebAssembly, nothing on aarch64 (not timed); `=2`
//! the plain loops for every table, `=3` for sparse tables, on any CPU. Where it does not act the
//! branch-free scan of BXSCANBF (graduated, always on) runs. The fused loop of LOSAT_X_BXDIAG keeps
//! its candidate loop. Where LOSAT_X_BXDIAG's fused loop runs by default (x86-64 with AVX2 and
//! POPCNT, LOSAT_X_BXCROSS; the WebAssembly build and x86-64 with SSSE3, SSE4.1 and POPCNT,
//! LOSAT_X_BXCROSS128) the direct loop scans the subjects of batches with few hits per used cell
//! (BXDIAGMIN: short queries); elsewhere every batch.

use super::x_ungap_sent::{by_word_size, XSent};
use super::*;
use crate::algorithm::blastp::extension::BlastpTwoHitUngappedResult;
use crate::algorithm::blastp::x_seed_fused::XHitBf;
use crate::algorithm::tblastx::x_aascan_bf::{x_scan_one_range_bf, X_SCAN_BF_SLACK};
use crate::algorithm::tblastx::x_scan_byte::{
    x_kept_presence, x_mode, x_scan_one_range_byte, x_table_presence, XPresence, XScanByte,
    X_SCAN_BYTE_SLACK,
};
use std::cell::RefCell;
use std::sync::OnceLock;

// No NCBI counterpart: reads the LOSAT_X_BXHITBF switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BXHITBF): whether the direct loop takes the too-far and overlap outcomes
/// without a branch where [`XHitBf::fits`] holds (module header); a switch of the strict set, on
/// unless it is `0`.
pub(crate) fn x_bxhitbf() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXHITBF"))
}

// No NCBI counterpart: reads the LOSAT_X_BXSCANBYTE switch once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_BXSCANBYTE): the direct loop's scan of the standard table (module header);
/// a switch of the strict set, on unless it is `0`.
pub(crate) fn x_bxscan_byte() -> XScanByte {
    static MODE: OnceLock<XScanByte> = OnceLock::new();
    *MODE.get_or_init(|| x_mode("LOSAT_X_BXSCANBYTE"))
}

// No NCBI counterpart: the presence vector in another layout, kept with the table; it does not
// change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXSCANBYTE): a standard table's presence bytes and whether the table is
/// sparse (`tblastx/x_scan_byte.rs` `x_table_presence`), built the first time the direct loop scans
/// with the table on presence bytes; `None` inside for a table they do not fit. It also carries the
/// table's build number (`build`), which names the table for the choice of the loop
/// (LOSAT_X_BXDIAGMIN, `x_seed_diag.rs`): a table built from pooled vectors
/// can sit at the address of the one dropped before it, with the same backbone vector; and whether
/// the table was built pooled (`pooled`), so that its vectors go back to the pool when it is
/// dropped (LOSAT_X_BXPOOLN1 builds a one-thread search's tables without the pool).
pub struct XKeptPresence {
    presence: OnceLock<Option<(XPresence, bool)>>,
    /// EXPERIMENT (LOSAT_X_BXCTXONCE): the context table of the query batch the table was built
    /// from (`x_seed_diag.rs` `x_with_context_table`), made by the first thread that needs it.
    pub(super) contexts: OnceLock<crate::algorithm::blastp::x_seed_fused::XContextTable>,
    pub(crate) build: u64,
    pub(crate) pooled: bool,
}

// No NCBI counterpart: numbers LOSAT's tables (speed only).
/// The next build number (from 1; 0 names no table).
static X_BUILDS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);

impl XKeptPresence {
    // No NCBI counterpart (see the struct).
    /// The companion of a table built now: no presence bytes yet and a new build number.
    pub fn new() -> Self {
        Self::x_built(false)
    }

    // No NCBI counterpart (see the struct).
    /// [`Self::new`] for a table built pooled or not (`pooled`).
    pub(crate) fn x_built(pooled: bool) -> Self {
        Self {
            presence: OnceLock::new(),
            contexts: OnceLock::new(),
            build: X_BUILDS.fetch_add(1, std::sync::atomic::Ordering::Relaxed),
            pooled,
        }
    }
}

impl Default for XKeptPresence {
    fn default() -> Self {
        Self::new()
    }
}

impl Lookup {
    // No NCBI counterpart: chooses LOSAT's scan; both give NCBI's pairs and scan ranges.
    /// EXPERIMENT (LOSAT_X_BXSCANBYTE): the presence bytes for `x_scan_one_range_byte` when `mode`
    /// scans this table that way (a standard table, the switch on, the table sparse if `mode`
    /// asks), built on the first call that needs them; else `None` (the scan before).
    fn x_scan_presence(&self, mode: XScanByte) -> Option<&XPresence> {
        match self {
            Self::Standard(table, kept) if mode.on() => {
                x_kept_presence(mode, kept.presence.get_or_init(|| x_table_presence(table)))
            }
            _ => None,
        }
    }
}

/// The `[X_SHADOW]` summaries of the BXDIAG, BXUNGAPSENT and BXSUBJFLOW switches (printed by `main`
/// at exit, one line per switch that ran in shadow mode).
// No NCBI counterpart: prints the shadow-mode counters; it does not change any value NCBI computes.
pub fn print_shadow_stats() {
    super::x_seed_diag::x_print_shadow_summary();
    super::x_ungap_sent::print_shadow_stats();
    #[cfg(feature = "parallel")]
    super::super::search::x_subject_flow::x_print_shadow_summary();
}

thread_local! {
    /// The per-thread offset-pair array (blast_engine.c:1040-1041).
    static X_PAIRS: RefCell<Vec<OffsetPair>> = const { RefCell::new(Vec::new()) };
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-607
/// `word_finder` without the per-call buffer and the batch copy (the seed loop of a search that
/// traces nothing).
/// EXPERIMENT (LOSAT_X_BXUNGAPSENT): with `sent` (the subject and its batch on that path,
/// `x_ungap_sent.rs`) the fused and the direct loops extend with the sentinel walks.
/// With the matrix of `scoring`: another matrix than BLOSUM62 takes [`x_word_finder_matrix`], its
/// instance out of line.
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_word_finder_scan<S: BxScoring>(
    scoring: S,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &BxLookup,
    diagonals: &mut Diagonals,
    sent: Option<XSent<'_, S>>,
) -> Result<Vec<InitHsp>> {
    if !S::BLOSUM62 {
        return x_word_finder_matrix(
            scoring, batch, parameters, options, subject, lookup, diagonals, sent,
        );
    }
    x_word_finder_scan_with(
        scoring, batch, parameters, options, subject, lookup, diagonals, sent,
    )
}

/// `x_word_finder_scan` with another matrix than BLOSUM62 (items SC6 and SC6b of line NP): the same
/// loops with the matrix's extensions (`MatrixBx`): the fused loop on the tables it pays for
/// (LOSAT_X_BXDIAG, x_seed_diag.rs), the sentinel walks with the matrix's keys where `sent` is
/// given (LOSAT_X_BXUNGAPSENT), the direct loop otherwise. Not inlined, so that BLOSUM62's loops
/// keep their code.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_word_finder_matrix<S: BxScoring>(
    scoring: S,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &BxLookup,
    diagonals: &mut Diagonals,
    sent: Option<XSent<'_, S>>,
) -> Result<Vec<InitHsp>> {
    x_word_finder_scan_with(
        scoring, batch, parameters, options, subject, lookup, diagonals, sent,
    )
}

/// The body of [`x_word_finder_scan`] with the matrix of `scoring`: BLOSUM62's direct loop without
/// `sent` is `x_word_finder_direct`, another matrix's `x_word_finder_matrix_direct`.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn x_word_finder_scan_with<S: BxScoring>(
    scoring: S,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &BxLookup,
    diagonals: &mut Diagonals,
    sent: Option<XSent<'_, S>>,
) -> Result<Vec<InitHsp>> {
    // No NCBI counterpart: records which extension this subject takes (LOSAT_X_REACH); every loop
    // below takes `sent`.
    if super::x_ungap_sent::on() {
        crate::utils::xreach::x_reach!("LOSAT_X_BXUNGAPSENT", sent.is_some());
    }
    // The table of 16-bit hits is scanned by the fused loop (x_seed_diag.rs), never by the loops of
    // `Int4` hits below.
    let lookup = match lookup {
        BxLookup::Wide(lookup) => lookup,
        BxLookup::Narrow(lookup) => {
            return super::x_seed_diag::x_word_finder_narrow(
                scoring, batch, parameters, options, subject, lookup, diagonals, sent,
            );
        }
    };
    // EXPERIMENT (LOSAT_X_BXDIAG): the fused loop where it acts (x_seed_diag.rs); where it does
    // not, the direct loop.
    if super::x_seed_diag::x_bxdiag_mode() != 0 {
        if let Some(hits) = super::x_seed_diag::x_word_finder_diag(
            scoring, batch, parameters, options, subject, lookup, diagonals, sent,
        )? {
            return Ok(hits);
        }
    }
    if let Some(sent) = sent {
        // The instance of the word size, chosen once for the subject (item SC9b-1).
        return by_word_size!(
            options.word_size,
            W => x_word_finder_direct_sent::<S, W>(
                batch,
                parameters,
                options,
                subject,
                lookup,
                diagonals,
                x_bxscan_byte(),
                x_bxhitbf(),
                sent,
            )
        );
    }
    if !S::BLOSUM62 {
        return x_word_finder_matrix_direct(
            scoring, batch, parameters, options, subject, lookup, diagonals,
        );
    }
    x_word_finder_direct(
        batch,
        parameters,
        options,
        subject,
        lookup,
        diagonals,
        x_bxscan_byte(),
        x_bxhitbf(),
    )
}

/// The direct loop of another matrix than BLOSUM62 (`x_word_finder_direct_with`, the port of NCBI's
/// loop) with the matrix's extensions (`MatrixBx`): a subject without the sentinel walks.
#[allow(clippy::too_many_arguments)]
fn x_word_finder_matrix_direct<S: BxScoring>(
    scoring: S,
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
) -> Result<Vec<InitHsp>> {
    let scan_byte = x_bxscan_byte();
    let hit_bf = x_bxhitbf();
    let acts = x_hitbf_acts(hit_bf, options, subject, diagonals);
    if hit_bf {
        crate::utils::xreach::x_reach!("LOSAT_X_BXHITBF", acts);
    }
    macro_rules! direct {
        ($hit_bf:literal) => {
            x_word_finder_direct_with::<$hit_bf, _, _>(
                batch,
                parameters,
                options,
                subject,
                lookup,
                diagonals,
                scan_byte,
                |q: &[u8], s: &[u8], s_left_off, s_right_off, q_right_off, dropoff, word| {
                    bx_two_hit(
                        scoring,
                        q,
                        s,
                        s_left_off,
                        s_right_off,
                        q_right_off,
                        dropoff,
                        word,
                    )
                },
                |q: &[u8], s: &[u8], q_off, s_off, dropoff, word| {
                    bx_one_hit(scoring, q, s, q_off, s_off, dropoff, word)
                },
            )
        };
    }
    if acts {
        direct!(true)
    } else {
        direct!(false)
    }
}

// No NCBI counterpart: chooses LOSAT's loop; both loops give the reference's hits and table.
/// EXPERIMENT (LOSAT_X_BXHITBF): whether the direct loop for `subject` is the branch-free one:
/// `hit_bf` (the switch) and [`XHitBf::fits`] at the table's offset before the subject.
fn x_hitbf_acts(
    hit_bf: bool,
    options: &ResolvedOptions,
    subject: &[u8],
    diagonals: &Diagonals,
) -> bool {
    hit_bf
        && XHitBf::fits(
            diagonals.offset,
            subject.len().saturating_sub(1),
            options.window_size,
            options.word_size,
        )
}

// No NCBI counterpart: chooses the instance of LOSAT's direct loop (`x_word_finder_direct_with`, which
// ports the NCBI code); both give the reference's hits and table.
// The direct loop: `Lookup::scan` and `word_finder` in one function (`x_word_finder_direct_with`).
// The standard table is scanned with `x_scan_one_range_bf` (module header; LOSAT_X_BXSCANBF
// graduated).
// EXPERIMENT (LOSAT_X_BXSCANBYTE): `scan_byte` scans it on presence bytes where that mode acts
// (`x_scan_one_range_byte`, before `x_scan_one_range_bf`).
// EXPERIMENT (LOSAT_X_BXHITBF): with `hit_bf`, where `XHitBf::fits` holds for the subject, the
// instance with the branch-free two-hit test, `x_word_finder_direct_bf` (module header).
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_word_finder_direct(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    scan_byte: XScanByte,
    hit_bf: bool,
) -> Result<Vec<InitHsp>> {
    let acts = x_hitbf_acts(hit_bf, options, subject, diagonals);
    if hit_bf {
        crate::utils::xreach::x_reach!("LOSAT_X_BXHITBF", acts);
    }
    if acts {
        return x_word_finder_direct_bf(
            batch, parameters, options, subject, lookup, diagonals, scan_byte,
        );
    }
    x_word_finder_direct_with::<false, _, _>(
        batch,
        parameters,
        options,
        subject,
        lookup,
        diagonals,
        scan_byte,
        extend_two_hit_blosum62,
        extend_one_hit_blosum62,
    )
}

// No NCBI counterpart: instances of `x_word_finder_direct_with`, kept out of line; they do not
// change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXUNGAPSENT): `x_word_finder_direct` with the extension of `sent`
/// (`x_ungap_sent.rs`), in a function of its own so that the loops without it keep their code.
/// With another matrix (item SC6b) the walks are its keys' and the one-hit argument is never
/// called (`sent` exists only with a two-hit window). `W` is `options.word_size` as a const (2 to
/// 4), the walks' instance: the caller chooses it once per subject (`by_word_size!`).
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_word_finder_direct_sent<S: BxScoring, const W: usize>(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    scan_byte: XScanByte,
    hit_bf: bool,
    sent: XSent<'_, S>,
) -> Result<Vec<InitHsp>> {
    // The one-hit argument below is BLOSUM62's; another matrix's `sent` exists only with a two-hit
    // window (`query_fits_of`), where it is never called.
    debug_assert!(S::BLOSUM62 || options.window_size > 0);
    let acts = x_hitbf_acts(hit_bf, options, subject, diagonals);
    if hit_bf {
        crate::utils::xreach::x_reach!("LOSAT_X_BXHITBF", acts);
    }
    macro_rules! direct {
        ($hit_bf:literal, $shadow:literal) => {
            x_word_finder_direct_with::<$hit_bf, _, _>(
                batch,
                parameters,
                options,
                subject,
                lookup,
                diagonals,
                scan_byte,
                move |_: &[u8], _: &[u8], s_left_off, s_right_off, q_right_off, dropoff, _| {
                    sent.two_hit::<$shadow, W>(s_left_off, s_right_off, q_right_off, dropoff)
                },
                extend_one_hit_blosum62,
            )
        };
    }
    match (acts, sent.shadow()) {
        (true, false) => direct!(true, false),
        (true, true) => direct!(true, true),
        (false, false) => direct!(false, false),
        (false, true) => direct!(false, true),
    }
}

// No NCBI counterpart: an instance of `x_word_finder_direct_with`, kept out of line; it does not
// change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXHITBF): the direct loop with the branch-free two-hit test, in a function
/// of its own so that the loop without the test keeps its code (module header).
#[inline(never)]
fn x_word_finder_direct_bf(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    scan_byte: XScanByte,
) -> Result<Vec<InitHsp>> {
    x_word_finder_direct_with::<true, _, _>(
        batch,
        parameters,
        options,
        subject,
        lookup,
        diagonals,
        scan_byte,
        extend_two_hit_blosum62,
        extend_one_hit_blosum62,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-511
// `Lookup::scan` (its range set-up and loop) and `word_finder` (the statements per pair)
// in one function; the pairs are read from the scanner's buffer. The standard table is scanned
// with `x_scan_one_range_bf` (module header).
// EXPERIMENT (LOSAT_X_BXSCANBYTE): where `scan_byte` acts for the table (`Lookup::x_scan_presence`)
// it scans the standard table with `x_scan_one_range_byte` instead (module header).
// EXPERIMENT (LOSAT_X_BXHITBF): with `HIT_BF` (only where the caller found `XHitBf::fits` to hold)
// each hit first goes through the branch-free two-hit test (`XHitBf::test`), and only a flag-1 cell
// or a candidate reaches the statements of `word_finder`.
// EXPERIMENT (LOSAT_X_BXUNGAPSENT): `extend` is the two-hit extension, `extend_two_hit_blosum62`'s
// arguments and result (the sentinel walks of `x_ungap_sent.rs` in `x_word_finder_direct_sent`).
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn x_word_finder_direct_with<const HIT_BF: bool, E, O>(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    scan_byte: XScanByte,
    extend: E,
    extend_one: O,
) -> Result<Vec<InitHsp>>
where
    E: Fn(&[u8], &[u8], usize, usize, usize, i32, usize) -> Option<BlastpTwoHitUngappedResult>,
    O: Fn(&[u8], &[u8], usize, usize, i32, usize) -> Option<BlastpOneHitUngappedResult>,
{
    let query = &batch.sequence_start[1..];
    let word = options.word_size;
    let mask = diagonals.entries.len() as i32 - 1;
    let mut hits = Vec::new();
    // Lookup::scan's set-up, with length = subject.len() - 1 as word_finder passes it.
    let length = subject.len().saturating_sub(1);
    let window = options.window_size;
    let lookup_word = match lookup {
        Lookup::Standard(l, _) => l.word_length,
        Lookup::Compressed(l) => x_compressed_word_length(l),
    };
    let ranges = [(0, length as i32)];
    let mut range = [0, 0, length as i32 - lookup_word];
    if window > 0 && range[2] < range[1] {
        range[2] = range[1];
    }
    // The scanners read `scanned`: the subject, or a padded copy when its one scan position's
    // word runs past the sentinel (seed.rs `x_scan_bytes`); the statements below read `subject`.
    let mut pad = Vec::new();
    let scanned = x_scan_bytes(subject, range[2], lookup_word, &mut pad);
    let capacity = lookup.capacity();
    // No NCBI counterpart: chooses LOSAT's scan and records the side taken per subject
    // (LOSAT_X_REACH); both scans give NCBI's pairs. EXPERIMENT (LOSAT_X_BXSCANBYTE): the presence
    // bytes where the scan on them runs for this table.
    let presence = lookup.x_scan_presence(scan_byte);
    if scan_byte.on() && matches!(lookup, Lookup::Standard(..)) {
        crate::utils::xreach::x_reach!("LOSAT_X_BXSCANBYTE", presence.is_some());
    }
    // Room for the pairs the branch-free scan or the scan on presence bytes (LOSAT_X_BXSCANBYTE)
    // writes past `array_size` (never read).
    let room = capacity
        + if presence.is_some() {
            X_SCAN_BYTE_SLACK
        } else {
            X_SCAN_BF_SLACK
        };
    // EXPERIMENT (LOSAT_X_BXHITBF): the constants of the branch-free test (the caller checked
    // XHitBf::fits).
    let bf = XHitBf::in_range(diagonals.offset, window, word);
    X_PAIRS.with(|buffer| {
        let mut buffer = buffer.borrow_mut();
        if buffer.len() < room {
            buffer.resize(room, OffsetPair::default());
        }
        let pairs = &mut buffer[..room];
        let array_size = capacity as i32;
        while range[1] <= range[2] {
            let n = match (lookup, presence) {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:75-127
                // (`s_BlastAaScanSubject`; tblastx/x_scan_byte.rs gives the loop): the same pairs
                // and `scan_range` on every call.
                (Lookup::Standard(l, _), Some(presence)) => x_scan_one_range_byte(
                    scan_byte, l, presence, scanned, pairs, array_size, &mut range,
                ),
                (Lookup::Standard(l, _), None) => {
                    x_scan_one_range_bf(l, scanned, pairs, array_size, &mut range)
                }
                (Lookup::Compressed(l), _) => {
                    x_compressed_scan(l, scanned, &ranges, pairs, array_size, &mut range)
                }
            };
            for pair in &pairs[..n as usize] {
                let (q, s) = (pair.q_off as i32, pair.s_off as i32);
                // EXPERIMENT (LOSAT_X_BXHITBF): a flag-0 cell's too-far and overlap outcomes
                // without a branch; a flag-1 cell or a candidate goes on to the statements below.
                if HIT_BF {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516
                    let index = ((q - s) & mask) as usize;
                    debug_assert!(index < diagonals.entries.len());
                    // SAFETY: 0 <= index <= mask = entries.len() - 1 (a power of two less one,
                    // `Diagonals::new`).
                    let cell = unsafe { diagonals.entries.get_unchecked_mut(index) };
                    if bf.test(cell, bf.here(s)) {
                        continue;
                    }
                }
                // From here to the end of the loop body: the statements of word_finder.
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:560,783-784
                let context_of = |q: i32| {
                    batch
                        .contexts
                        .partition_point(|c| c.offset as i32 <= q)
                        .saturating_sub(1)
                };
                let context;
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516,776
                let index = if options.window_size == 0 {
                    (s - q) & mask
                } else {
                    (q - s) & mask
                } as usize;
                let entry = &mut diagonals.entries[index];
                let (u, end, extended) = if options.window_size == 0 {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:777-781
                    if s.wrapping_sub(entry.last_hit().wrapping_sub(diagonals.offset)) < 0 {
                        continue;
                    }
                    context = context_of(q);
                    let p = &parameters[context];
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:787-790
                    let Some(result) = extend_one(
                        query,
                        subject,
                        q as usize,
                        s as usize,
                        p.word_xdrop,
                        word as usize,
                    ) else {
                        continue;
                    };
                    (result.ungapped_data, result.s_last_off, true)
                } else {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:519-551
                    if entry.flag() != 0 {
                        if s.wrapping_add(diagonals.offset) < entry.last_hit() {
                            continue;
                        }
                        entry.set_last_hit(s.wrapping_add(diagonals.offset));
                        entry.set_flag(0);
                        continue;
                    }
                    let previous = entry.last_hit().wrapping_sub(diagonals.offset);
                    let diff = s.wrapping_sub(previous);
                    if diff >= options.window_size {
                        entry.set_last_hit(s.wrapping_add(diagonals.offset));
                        continue;
                    }
                    if diff < word {
                        continue;
                    }
                    context = context_of(q);
                    let p = &parameters[context];
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:566-573
                    if q - diff < batch.contexts[context].offset as i32 {
                        entry.set_last_hit(s.wrapping_add(diagonals.offset));
                        continue;
                    }
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:576-582
                    let Some(result) = extend(
                        query,
                        subject,
                        (previous + word) as usize,
                        s as usize,
                        q as usize,
                        p.word_xdrop,
                        word as usize,
                    ) else {
                        continue;
                    };
                    (result.ungapped_data, result.s_last_off, result.right_extend)
                };
                let p = &parameters[context];
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:588-606
                if u.score >= p.word_cutoff {
                    hits.push(InitHsp {
                        q_seed: q,
                        s_seed: s,
                        q_start: u.q_start,
                        s_start: u.s_start,
                        length: u.length,
                        score: u.score,
                    });
                }
                if options.window_size == 0 {
                    entry.set_last_hit((end - (word - 1)).wrapping_add(diagonals.offset));
                } else if extended {
                    entry.set_flag(1);
                    entry.set_last_hit((end - (word - 1)).wrapping_add(diagonals.offset));
                } else {
                    entry.set_last_hit(s.wrapping_add(diagonals.offset));
                }
            }
        }
    });
    // The rest of word_finder, unchanged.
    diagonals.finish_subject(subject.len().saturating_sub(1))?;
    hits.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.s_start.cmp(&b.s_start))
            .then(b.length.cmp(&a.length))
            .then(a.q_start.cmp(&b.q_start))
    });
    Ok(hits)
}

#[cfg(test)]
#[path = "x_seed_scan_tests.rs"]
mod tests;
