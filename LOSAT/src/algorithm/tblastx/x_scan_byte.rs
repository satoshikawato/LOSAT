//! EXPERIMENT (LOSAT_X_TBNSCANBYTE): tblastn's lookup scan with the table indexes of a block of
//! words made at once and their presence read from one byte per index, and the hits of a block
//! copied with one capacity test.
//!
//! The module lives in tblastx, the lowest program that scans with it (tblastx's dense unit-group
//! tables, LOSAT_X_TBXSCANBYTE; tblastn's units; blastx's direct loop), so that every caller reaches
//! it down the layering of `tests/check_program_edges.py`; it reads only tblastx's lookup table.
//!
//! NCBI's scan (c++/src/algo/blast/core/blast_aascan.c:75-127, `s_BlastAaScanSubject`, kept as
//! `tblastx/blast_aascan.rs` `s_blast_aa_scan_subject_one_range`) rolls the index of the next word
//! from the one before (`ComputeTableIndexIncremental`), tests its presence-vector bit and copies
//! the word's hits when the offset array has room. The branch-free scan (`tblastx/x_aascan_bf.rs`) takes the
//! branch on the bit away, in blocks of `X_BLOCK` words; on a genome searched with one protein
//! (HsG_p1 x hsap_dna `-db`, n1) the scan was still 55-58 % of the cycles, 5.38 cycles per word,
//! with the rolling index a chain from word to word (line F probe P3, v4-F/REPORT.md).
//!
//! For a table of word size 3 and 5-bit letters (tblastn's, `fits`) the index of the word at `s`
//! is `(subject[s] << 10 | subject[s + 1] << 5 | subject[s + 2]) & 0x7fff`: the rolled index is
//! `((index << 5) | letter) & mask`, and the letters before `s` are shifted past bit 14 and masked
//! away, so both give the same number for any bytes. Here each block's indexes are made from the
//! letters directly, 16 at a time with AVX2 (or by a plain loop), with no chain between words. Its
//! candidates are listed from `presence`, the presence vector with one byte per index (1 where
//! `pv` has the bit, 0 elsewhere, 32 KB for the 15-bit index space), so the test is one load
//! instead of a load, a shift and a mask. Indexes past the backbone read 0 and are never
//! candidates.
//!
//! The copy then takes the block's candidates in position order, as NCBI's loop does, each once.
//! While three hits for every candidate left in the block fit in the offset array (one test, made
//! again only after a cell of more hits), the hits of a cell of up to `AA_HITS_PER_CELL` are
//! written without the capacity test, as three pairs of one 8-byte store each (with AVX2 four pairs
//! in one 32-byte store: the cell broadcast, one `vpermd`, one `vpblendd` with the subject offset)
//! of which the first `num_used` count. A cell of more hits, and every candidate met while the
//! block could fill the array, takes NCBI's test; a larger cell's hits are then written four pairs
//! at a time from the overflow array. When a word's hits do not fit, `s_range[1]` is that word's
//! position and the scan returns, so the next call resumes there. The offset pairs and the
//! `s_range` values of every call are NCBI's (the words after the stop were looked at but not
//! used; the next call looks at them again, as with the branch-free scan).
//!
//! EXPERIMENT (LOSAT_X_TBNSCANRESUME): the unit path keeps the block where a call stops (its
//! indexes, its candidates and the word that did not fit, `XScanHeld`), and the next call of the
//! same range copies on from that candidate instead of making the block's candidates again
//! (`x_scan_one_range_held`; the copy loop is the one above). A stop then costs a call, not a
//! block of words, and the unit path gives the offset array NCBI's size, 4096 + the longest chain
//! (`search_seed.rs` `x_scan_unit`; `lookup_wrap.c:264-266`; at most the plan's size), instead of
//! `LOSAT_X_TBNPAIRCAP`'s 32768 + it. The pairs and the `s_range` of every call are still NCBI's.
//!
//! EXPERIMENT (LOSAT_X_CHAINRUN): with the AVX2 forms, a cell of more than `AA_HITS_PER_CELL` hits
//! (a chain in the overflow array) is written as one run of `X_CHAIN_RUN` pairs, six 32-byte stores
//! with no test between them, when the overflow array holds `X_CHAIN_RUN` query offsets from the
//! chain's start; the chain's hits after the run are written four pairs at a time, as before, and a
//! chain nearer the array's end takes the copy before. The copy before loops over every chain four
//! pairs at a time, and the end of that loop, set by the cell's hit count, which changes from
//! candidate to candidate, is a mispredicted branch for most chains: on BUSCO's eukaryota
//! ancestral variants x the yeast genome (`-db`, one thread; 1.26 G candidate words of 11.9 hits on
//! average, 84 % of them chains, 76 % of at most 16 hits) 69 % of the search's 2.48 G mispredicted
//! branches were in the copy, which took 31 % of the cycles at 1.1 instructions per cycle. The run
//! covers 88 % of the chains there. A cell of up to `AA_HITS_PER_CELL` hits keeps its one store, so
//! a table whose candidates are such cells (one protein x the human genome: 99.9 % of 1.23 G
//! candidates) runs the copy before. The pairs past a cell's hits are not part of the result, as
//! before, and every call returns the same pairs and `s_range` (unit test
//! `x_scan_byte_equals_ncbi_scan_on_random_subjects`). On a CPU without AVX2 and with the plain
//! loops nothing changes. `LOSAT_X_CHAINRUN=0` runs the copy before.
//!
//! The offset array must have `X_SCAN_BYTE_SLACK` entries after `array_size`; the entries past the
//! returned count are not part of the result. `presence` depends on the lookup table only; the
//! first thread that asks for a table builds it once and every thread reads it (`presence_bytes`).
//!
//! Probe and source of the form: line F s3 P3 V2 (candidates) and s5 P3b arm 7 (copy),
//! `exp/f2-tbn-scan` ff918e61: HsG_p1 x hsap_dna `-db` n1 whole-run cycles x0.703 of dev, with
//! the offset pairs and `s_range` of every call and the output equal.
//!
//! Both stages have an AVX2 form and a plain-loop form (the same values; the plain index loop is
//! left to the compiler's vectorizer). On the workstation (one core, the scan before this switch =
//! the branch-free scan) the AVX2 forms take x0.75 of the time of one protein x the human genome and
//! x0.98-0.99 where most words are candidates (many proteins x a genome, 55 proteins x 16S); the
//! plain loops take x0.79 and x0.83 (the LOSAT Web serial reactor, one protein x E. coli) of the
//! time where few words are candidates, but x1.01-1.06 where most are: their copy of a cell's hits
//! costs more than the copy before. So by default the switch takes the AVX2 forms where the CPU
//! has AVX2, and elsewhere the plain loops for sparse tables only: those where at most three tenths
//! of the backbone cells have hits (`X_SPARSE_TENTHS`; known when the presence bytes are built).
//! The first k proteins of AP027131 x E. coli (one core): the plain loops against the scan before
//! took x0.84-0.95 in the LOSAT Web serial reactor and cycles x0.84-0.90 on x86-64 at a fill of
//! 0.09-0.24 (k = 1-8), x1.00 and x0.95 at 0.34 (k = 16), and x1.02-1.03 and x0.99-1.00 at 0.42-0.46.
//!
//! The switch is read once. It is in the strict set (`utils/xswitch.rs`): on unless
//! `LOSAT_X_TBNSCANBYTE=0`; on, it takes the AVX2 forms
//! of both stages on x86-64 with AVX2, the plain loops for sparse tables on x86-64 without AVX2 and
//! in WebAssembly, and changes nothing elsewhere (aarch64: not timed). `LOSAT_X_TBNSCANBYTE=2`
//! takes the plain loops on any CPU for every table, `=3` the plain loops for sparse tables on any
//! CPU. A table of another shape, and without AVX2 a table that is not sparse, is scanned as before
//! (the branch-free scan of `tblastx/x_aascan_bf.rs`). With `LOSAT_X_STATS` set, `main` prints
//! `[X_STATS] TBNSCANBYTE units= calls= avx2=` (the units scanned this way, their scan calls, and
//! whether the AVX2 forms ran).
//!
//! EXPERIMENT (LOSAT_X_BXSCANBYTE): blastx's direct seed loop (`blastx/x_seed_scan.rs`) scans its
//! standard table (word size 3, 5-bit letters: every one `fits`) with `x_scan_one_range_byte` too,
//! under its own switch of the same values and CPU policy (`x_mode`). Its presence bytes are built
//! once per table and kept with it (`x_table_presence`; blastx's `Lookup`).

use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair;
use crate::algorithm::tblastx::lookup::{
    BackboneCell, BlastAaLookupTable, AA_HITS_PER_CELL, PV_ARRAY_BTS, PV_ARRAY_MASK,
};
use std::mem::MaybeUninit;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): subject words whose indexes are made before their hits are
/// copied.
const X_BLOCK: usize = 256;

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): the index mask of a table of word size 3 and 5-bit letters.
const X_MASK: usize = (1 << 15) - 1;

/// EXPERIMENT (LOSAT_X_CHAINRUN): the pairs written for a chain before its hit count is read.
const X_CHAIN_RUN: usize = 24;
const _: () = assert!(X_CHAIN_RUN.is_multiple_of(4) && X_CHAIN_RUN > AA_HITS_PER_CELL + 1);

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): offset-array entries needed after `array_size`. A cell of up
/// to `AA_HITS_PER_CELL` hits is written as `AA_HITS_PER_CELL + 1` pairs, a chain of `n` hits as
/// `n` rounded up to four, or as at least `X_CHAIN_RUN` (LOSAT_X_CHAINRUN): a chain is written only
/// where its `n > AA_HITS_PER_CELL` hits fit, so its run ends at most `X_CHAIN_RUN -
/// (AA_HITS_PER_CELL + 1)` entries after `array_size`.
pub(crate) const X_SCAN_BYTE_SLACK: usize = X_CHAIN_RUN - (AA_HITS_PER_CELL + 1);
const _: () = assert!(X_SCAN_BYTE_SLACK > AA_HITS_PER_CELL);

/// A table is sparse when at most `X_SPARSE_TENTHS` tenths of its backbone cells have hits
/// (presence bytes of 1): the plain loops then win (`XScanByte::Sparse`).
const X_SPARSE_TENTHS: usize = 3;

/// Whether CPUs without AVX2 take the plain loops for sparse tables by default: on the targets they
/// were timed on (x86-64 without AVX2, WebAssembly); elsewhere (aarch64) only with `=3`.
const X_SPARSE_DEFAULT: bool = cfg!(any(target_arch = "x86_64", target_arch = "wasm32"));

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): how the unit path scans a table that `fits`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum XScanByte {
    /// As before the switch (the branch-free scan of `tblastx/x_aascan_bf.rs`).
    Off,
    /// `x_scan_one_range_byte`; `avx2` = the AVX2 forms of both stages.
    On { avx2: bool },
    /// `x_scan_one_range_byte` with the plain loops for a sparse table (`presence_bytes`), as
    /// before for the others.
    Sparse,
}

impl XScanByte {
    /// Whether `x_scan_one_range_byte` runs.
    pub(crate) fn on(self) -> bool {
        self != XScanByte::Off
    }
}

// No NCBI counterpart: reads LOSAT_X_TBNSCANBYTE once and detects AVX2; it does not change any
// value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): the scan of the unit path. LOSAT_X_TBNSCANBYTE is a switch of
/// the strict set, on unless it is `0`. On, it takes the
/// AVX2 forms where the CPU has AVX2, and elsewhere the plain loops for sparse tables only (the
/// plain loops lose where most words are candidates, see the module text; on aarch64 nothing);
/// `2` takes the plain loops on any CPU and table, `3` the plain loops for sparse tables on any
/// CPU.
pub(crate) fn x_scan_byte() -> XScanByte {
    static MODE: std::sync::OnceLock<XScanByte> = std::sync::OnceLock::new();
    *MODE.get_or_init(|| x_mode("LOSAT_X_TBNSCANBYTE"))
}

// No NCBI counterpart: reads a switch and detects AVX2; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNSCANBYTE, LOSAT_X_BXSCANBYTE): the scan that the strict-set switch `name`
/// (one of these two) asks for on this CPU, read now (the callers keep it): `x_scan_byte`'s values.
pub(crate) fn x_mode(name: &str) -> XScanByte {
    let value = crate::utils::xswitch::value(name, 1..=3, 1);
    let avx2 = !matches!(value, Some(2 | 3)) && avx2_detected();
    if value == Some(2) || (value.is_some() && avx2) {
        XScanByte::On { avx2 }
    } else if value == Some(3) || (value == Some(1) && X_SPARSE_DEFAULT) {
        XScanByte::Sparse
    } else {
        XScanByte::Off
    }
}

/// Whether this machine has AVX2 (always false off x86-64).
fn avx2_detected() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): whether `x_scan_one_range_byte` can scan this table (word
/// size 3, 5-bit letters, a 15-bit mask).
pub(super) fn fits(lookup: &BlastAaLookupTable) -> bool {
    lookup.word_length == 3 && lookup.charsize == 5 && lookup.mask as usize == X_MASK
}

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): the presence bytes of a table (`build_presence`), one byte
/// per index of the 15-bit space, each 0 or 1. Only `build_presence` makes them, so the scan may
/// move its list of candidates on by a word's byte without a test (`list`).
#[derive(Clone)]
pub(crate) struct XPresence(Arc<[u8; X_MASK + 1]>);

impl std::ops::Deref for XPresence {
    type Target = [u8; X_MASK + 1];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:55-57
/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): the presence vector of `lookup` with one byte per index of
/// the 15-bit space: 1 where `pv` has the bit and the backbone has the cell, 0 elsewhere (a bit is
/// set only for a cell with hits, so this is `pv`; the second test keeps every candidate's cell
/// inside the backbone whatever `pv` holds).
fn build_presence(lookup: &BlastAaLookupTable) -> XPresence {
    let bytes: Arc<[u8]> = (0..=X_MASK)
        .map(|ix| {
            let bit = lookup
                .pv
                .get(ix >> PV_ARRAY_BTS)
                .map_or(0, |word| ((word >> (ix & PV_ARRAY_MASK)) & 1) as u8);
            if ix < lookup.backbone.len() {
                bit
            } else {
                0
            }
        })
        .collect();
    XPresence(bytes.try_into().expect("one byte per index"))
}

/// A table, its presence bytes and whether it is sparse. The `Weak` keeps the table's allocation,
/// so no other table can have its address while the entry is kept.
type PresenceEntry = (std::sync::Weak<BlastAaLookupTable>, XPresence, bool);

/// The entry `presence_bytes` gave last on this thread (its bytes are those of `SHARED`).
type KeptPresence = Option<PresenceEntry>;

thread_local! {
    /// The presence bytes this thread read last (`KeptPresence`).
    static PRESENCE: std::cell::RefCell<KeptPresence> = const { std::cell::RefCell::new(None) };
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h:95-140
// NCBI's presence vector is part of its table, one array that every search thread reads. These
// are the presence bytes of the live tables, one array per table that every thread reads: an
// entry whose table is gone goes when a thread next looks up a table. It does not change any value
// NCBI computes.
static SHARED: std::sync::Mutex<Vec<PresenceEntry>> = std::sync::Mutex::new(Vec::new());

// No NCBI counterpart: the presence vector in another layout, made once per table; it does not
// change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): the presence bytes of `lookup` (see `build_presence`) and
/// whether it is sparse, built by the first thread that asks for this table and read by every
/// thread after it (`SHARED`); a thread asks again only when it saw another table last. The table
/// must `fit`.
pub(super) fn presence_bytes(lookup: &Arc<BlastAaLookupTable>) -> (XPresence, bool) {
    PRESENCE.with(|kept| {
        let mut kept = kept.borrow_mut();
        if let Some((table, bytes, sparse)) = &*kept {
            if table.as_ptr() == Arc::as_ptr(lookup) {
                return (bytes.clone(), *sparse);
            }
        }
        let (bytes, sparse) = shared_presence(lookup);
        *kept = Some((Arc::downgrade(lookup), bytes.clone(), sparse));
        (bytes, sparse)
    })
}

// No NCBI counterpart: finds or makes the entry of `SHARED`; it does not change any value NCBI
// computes.
/// The `SHARED` entry of `lookup`, made now if there is none; the entries of tables that are gone
/// go first.
fn shared_presence(lookup: &Arc<BlastAaLookupTable>) -> (XPresence, bool) {
    let mut shared = SHARED
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    shared.retain(|(table, _, _)| table.strong_count() > 0);
    if let Some((_, bytes, sparse)) = shared
        .iter()
        .find(|(table, _, _)| table.as_ptr() == Arc::as_ptr(lookup))
    {
        return (bytes.clone(), *sparse);
    }
    let bytes = build_presence(lookup);
    let sparse = sparse(&bytes, lookup);
    shared.push((Arc::downgrade(lookup), bytes.clone(), sparse));
    (bytes, sparse)
}

// No NCBI counterpart: chooses LOSAT's scan; it does not change any value NCBI computes.
/// Whether at most `X_SPARSE_TENTHS` tenths of the backbone cells of `lookup` have hits.
fn sparse(presence: &[u8; X_MASK + 1], lookup: &BlastAaLookupTable) -> bool {
    let ones = presence.iter().filter(|&&b| b == 1).count();
    sparse_ones(ones, lookup)
}

// No NCBI counterpart: chooses LOSAT's scan; it does not change any value NCBI computes.
/// `sparse` for a table whose presence bytes hold `ones` ones.
fn sparse_ones(ones: usize, lookup: &BlastAaLookupTable) -> bool {
    ones * 10 <= lookup.backbone.len() * X_SPARSE_TENTHS
}

// No NCBI counterpart: counts the presence bits of `lookup`; it does not change any value NCBI
// computes.
/// The ones `build_presence` would write for `lookup`: the bits of `pv` at the indexes of the
/// 15-bit space inside the backbone.
fn presence_ones(lookup: &BlastAaLookupTable) -> usize {
    let end = lookup.backbone.len().min(X_MASK + 1);
    lookup
        .pv
        .iter()
        .enumerate()
        .map(|(w, &word)| {
            let first = w << PV_ARRAY_BTS;
            let bits = end.saturating_sub(first).min(PV_ARRAY_MASK + 1);
            let mask = if bits > PV_ARRAY_MASK {
                u32::MAX
            } else {
                (1u32 << bits) - 1
            };
            (word & mask).count_ones() as usize
        })
        .sum()
}

// No NCBI counterpart: chooses LOSAT's scan; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBXSCANBYTE): whether `mode` scans `lookup` on presence bytes (it is on, the
/// table `fits`, and it is sparse if `mode` asks), decided from `pv` without the bytes: the
/// same answer as `x_kept_presence(mode, &x_table_presence(lookup)).is_some()`, for a caller that
/// builds the bytes when a scan first needs them (`x_build_presence`; tblastx's unit-group tables).
pub(crate) fn x_table_takes(mode: XScanByte, lookup: &BlastAaLookupTable) -> bool {
    mode.on()
        && fits(lookup)
        && (mode != XScanByte::Sparse || sparse_ones(presence_ones(lookup), lookup))
}

// No NCBI counterpart: the presence vector in another layout; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_TBXSCANBYTE): the presence bytes of a table that `fits` (`build_presence`),
/// for a caller that decided with `x_table_takes`.
pub(crate) fn x_build_presence(lookup: &BlastAaLookupTable) -> XPresence {
    debug_assert!(fits(lookup), "presence bytes of a table that does not fit");
    build_presence(lookup)
}

// No NCBI counterpart: chooses LOSAT's scan; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): the presence bytes for `x_scan_one_range_byte` when `mode`
/// scans `lookup` this way (it is on, the table `fits`, and it is sparse if `mode` asks), else
/// `None` (the scan before).
pub(crate) fn x_presence(mode: XScanByte, lookup: &Arc<BlastAaLookupTable>) -> Option<XPresence> {
    if !mode.on() || !fits(lookup) {
        return None;
    }
    let (bytes, sparse) = presence_bytes(lookup);
    (mode != XScanByte::Sparse || sparse).then_some(bytes)
}

// No NCBI counterpart: the presence vector in another layout, made once per table; it does not
// change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXSCANBYTE): the presence bytes of `lookup` (see `build_presence`) and
/// whether it is sparse, for a caller that keeps them with the table (blastx's `Lookup`); `None`
/// for a table that does not `fit`.
pub(crate) fn x_table_presence(lookup: &BlastAaLookupTable) -> Option<(XPresence, bool)> {
    fits(lookup).then(|| {
        let bytes = build_presence(lookup);
        let sparse = sparse(&bytes, lookup);
        (bytes, sparse)
    })
}

// No NCBI counterpart: chooses LOSAT's scan; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXSCANBYTE): the bytes of `kept` (from `x_table_presence`) when `mode`
/// scans that table on presence bytes (it is on, and the table is sparse if `mode` asks), else
/// `None` (the scan before).
pub(crate) fn x_kept_presence(
    mode: XScanByte,
    kept: &Option<(XPresence, bool)>,
) -> Option<&XPresence> {
    match kept {
        Some((bytes, sparse)) if mode.on() && (mode != XScanByte::Sparse || *sparse) => Some(bytes),
        _ => None,
    }
}

// No NCBI counterpart: LOSAT_X_STATS counters; they do not change any value NCBI computes.
/// Units scanned by `x_scan_one_range_byte`, their scan calls, and the calls whose AVX2 copy wrote
/// chains as runs (`[0]`, `[1]`, `[2]`; LOSAT_X_CHAINRUN).
static COUNTS: [AtomicU64; 3] = [const { AtomicU64::new(0) }; 3];

/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): counts one unit and its `calls` for the
/// `[X_STATS] TBNSCANBYTE` line (`LOSAT_X_STATS` set).
pub(crate) fn x_count_unit(calls: u64) {
    if crate::algorithm::tblastx::x_aascan_bf::x_stats() {
        COUNTS[0].fetch_add(1, Ordering::Relaxed);
        COUNTS[1].fetch_add(calls, Ordering::Relaxed);
        if x_scan_byte() == (XScanByte::On { avx2: true }) && x_chain_run() {
            COUNTS[2].fetch_add(calls, Ordering::Relaxed);
        }
    }
}

/// The `[X_STATS] TBNSCANBYTE` line (printed by `main` at exit) with `LOSAT_X_STATS`, when the unit
/// path scanned any unit this way.
// No NCBI counterpart: prints the counters above; it does not change any value NCBI computes.
pub fn print_summary() {
    let [units, calls, runs] = COUNTS.each_ref().map(|c| c.load(Ordering::Relaxed));
    let mode = x_scan_byte_if_read();
    if units > 0 {
        let avx2 = matches!(mode, Some(XScanByte::On { avx2: true }));
        eprintln!("[X_STATS] TBNSCANBYTE units={units} calls={calls} avx2={avx2} runs={runs}");
    }
}

/// The mode if a search read it (the summary must not read the switch for programs that never
/// scanned).
fn x_scan_byte_if_read() -> Option<XScanByte> {
    (COUNTS[0].load(Ordering::Relaxed) > 0).then(x_scan_byte)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:75-127 (cited at `scan`).
/// EXPERIMENT (LOSAT_X_TBNSCANBYTE): `s_blast_aa_scan_subject_one_range` for a table that `fits`
/// (see the module text). Same arguments, same result, same `s_range`; `presence` is
/// `presence_bytes(lookup)` and `offset_pairs` must hold `array_size + X_SCAN_BYTE_SLACK` entries.
/// Not inlined: the switch-off scan of `x_scan_unit` keeps its code. Blastx's direct loop calls it
/// too (LOSAT_X_BXSCANBYTE).
#[inline(never)]
pub(crate) fn x_scan_one_range_byte(
    mode: XScanByte,
    lookup: &BlastAaLookupTable,
    presence: &XPresence,
    subject: &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
) -> i32 {
    x_scan_one_range_held(
        mode,
        lookup,
        presence,
        subject,
        offset_pairs,
        array_size,
        s_range,
        &mut XScanHeld::new(),
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:75-127 (cited at `scan`).
/// EXPERIMENT (LOSAT_X_TBNSCANRESUME): `x_scan_one_range_byte` that keeps in `held` the block of a
/// call that stops on a full offset array, and copies the rest of it when the next call is that
/// call's resumption (the same subject and `s_range[2]`, `s_range[1]` the word that did not fit).
/// Same result and `s_range` as NCBI's scan. A new `held` holds nothing.
#[allow(clippy::too_many_arguments)]
#[inline(never)]
pub(crate) fn x_scan_one_range_held(
    mode: XScanByte,
    lookup: &BlastAaLookupTable,
    presence: &XPresence,
    subject: &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
    held: &mut XScanHeld,
) -> i32 {
    match mode {
        XScanByte::On { .. } | XScanByte::Sparse => scan(
            mode == XScanByte::On { avx2: true },
            x_chain_run(),
            lookup,
            presence,
            subject,
            offset_pairs,
            array_size,
            s_range,
            held,
        ),
        XScanByte::Off => unreachable!("x_scan_one_range_byte with the switch off"),
    }
}

/// EXPERIMENT (LOSAT_X_TBNSCANRESUME): the block of the scan call that stopped on a full offset
/// array, for its resumption (`x_scan_one_range_held`): the indexes and candidates of the block,
/// the first candidate not copied, and the call that resumes it.
pub(crate) struct XScanHeld {
    index: [u16; X_BLOCK],
    cand: [MaybeUninit<u16>; X_BLOCK],
    /// The block's candidates (`cand[..n]`); a block is held while `next < n`.
    n: usize,
    /// The first candidate not copied.
    next: usize,
    /// The subject offset of the block's first word.
    block_first: usize,
    /// The resuming call: the subject's address and length, `s_range[1]` and `s_range[2]`.
    at: (usize, usize, i32, i32),
}

impl XScanHeld {
    /// One that holds nothing.
    pub(crate) fn new() -> Self {
        XScanHeld {
            index: [0; X_BLOCK],
            cand: [MaybeUninit::uninit(); X_BLOCK],
            n: 0,
            next: 0,
            block_first: 0,
            at: (0, 0, 0, 0),
        }
    }
}

// No NCBI counterpart: reads the switch once; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_TBNSCANRESUME): whether the unit path's scan on presence bytes resumes a
/// stopped call in its block (`x_scan_one_range_held`) and takes NCBI's offset-array size
/// (`x_scan_unit`). In the strict set: on unless `LOSAT_X_TBNSCANRESUME=0`; read once.
pub(crate) fn x_scan_resume() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNSCANRESUME"))
}

// No NCBI counterpart: reads LOSAT_X_CHAINRUN once; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_CHAINRUN): whether the AVX2 copy writes a chain as one run of
/// `X_CHAIN_RUN` pairs (module text). A switch of the strict set, on unless it is `0`.
fn x_chain_run() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_CHAINRUN"))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:83-127
/// The scan of one call (module text), with the AVX2 forms when `avx2` (and their chain runs when
/// `run`, LOSAT_X_CHAINRUN); the block where it stops is kept in `held`, and copied on from there
/// when this call resumes it.
#[allow(clippy::too_many_arguments)]
fn scan(
    avx2: bool,
    run: bool,
    lookup: &BlastAaLookupTable,
    presence: &XPresence,
    subject: &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
    held: &mut XScanHeld,
) -> i32 {
    if s_range[1] > s_range[2] {
        return 0;
    }
    assert!(
        fits(lookup),
        "LOSAT_X_TBNSCANBYTE: table of word size {}",
        lookup.word_length
    );
    assert!(
        array_size >= 0 && offset_pairs.len() >= array_size as usize + X_SCAN_BYTE_SLACK,
        "offset array of {} entries for array_size {array_size}",
        offset_pairs.len()
    );
    #[cfg(target_arch = "x86_64")]
    debug_assert!(!avx2 || avx2_detected(), "AVX2 forms without AVX2");
    let s_first = s_range[1] as usize;
    let s_last = s_range[2] as usize;
    // NCBI's loop reads `subject[s + 2]` for every `s <= s_last`. The callers keep every word in
    // `subject`: TBLASTN's `x_scan_unit` does not scan a raised first-range word past it
    // (aa_ungapped.c:499-500), and BLASTX pads its subject (`x_scan_bytes`).
    assert!(
        s_last + 2 < subject.len(),
        "LOSAT_X_TBNSCANBYTE: the word at {s_last} runs past a subject of {} letters",
        subject.len()
    );
    let table = Table {
        backbone: lookup.backbone.as_ptr(),
        overflow: lookup.overflow.as_ptr(),
        overflow_len: lookup.overflow.len(),
    };
    let dest = offset_pairs.as_mut_ptr();
    let mut totalhits: i32 = 0;
    let mut block_first = s_first;
    // The held block, from its first candidate not copied, when this call resumes the one that
    // stopped in it (its first candidate is the word at `s_range[1]`).
    let subject_at = (subject.as_ptr() as usize, subject.len());
    let mut first = 0;
    if held.next < held.n && held.at == (subject_at.0, subject_at.1, s_range[1], s_range[2]) {
        first = held.next;
        block_first = held.block_first;
    }
    while block_first <= s_last {
        let block_last = s_last.min(block_first + X_BLOCK - 1);
        if first == 0 {
            let letters = &subject[block_first..block_last + 3];
            held.n = candidates_any(avx2, letters, presence, &mut held.index, &mut held.cand);
        }
        // SAFETY: `presence` has an entry for every 15-bit index and is 1 only at indexes of the
        // backbone (`build_presence`), so every candidate's cell is inside it; `dest` has
        // `array_size + X_SCAN_BYTE_SLACK` entries; AVX2 only when the caller detected it.
        let (total, stop) = unsafe {
            copy_any(
                avx2,
                run,
                &table,
                dest,
                held,
                first,
                block_first,
                totalhits,
                array_size,
            )
        };
        if let Some(s_off) = stop {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:118-123
            s_range[1] = s_off as i32;
            // The word's place among the block's candidates (which ascend), for the next call.
            // SAFETY: `cand[..n]` was written by the candidate stage.
            let cand =
                unsafe { &*(&held.cand[..held.n] as *const [MaybeUninit<u16>] as *const [u16]) };
            held.next = cand.partition_point(|&k| (k as usize) < s_off as usize - block_first);
            held.block_first = block_first;
            held.at = (subject_at.0, subject_at.1, s_range[1], s_range[2]);
            return total;
        }
        totalhits = total;
        first = 0;
        block_first = block_last + 1;
    }
    held.n = 0;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:125-127
    s_range[1] = (s_last + 1) as i32;
    totalhits
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:96-107,121-127,55-57
// (cited at `candidates`).
/// The candidate stage of a block (`candidates`, or `candidates_avx2` when `avx2`).
#[inline(always)]
fn candidates_any(
    avx2: bool,
    letters: &[u8],
    presence: &[u8; X_MASK + 1],
    index: &mut [u16; X_BLOCK],
    cand: &mut [MaybeUninit<u16>; X_BLOCK],
) -> usize {
    #[cfg(target_arch = "x86_64")]
    if avx2 {
        // SAFETY: AVX2 only when the caller detected it.
        return unsafe { candidates_avx2(letters, presence, index, cand) };
    }
    let _ = avx2;
    candidates(letters, presence, index, cand)
}

/// The copy stage of `held`'s block from its candidate `first` (`copy`, or `copy_avx2` when
/// `avx2`, `copy_avx2_run` when `run` too).
///
/// # Safety
/// As `copy_block`; AVX2 must be there when `avx2`.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
unsafe fn copy_any(
    avx2: bool,
    run: bool,
    table: &Table,
    dest: *mut BlastOffsetPair,
    held: &XScanHeld,
    first: usize,
    block_first: usize,
    totalhits: i32,
    array_size: i32,
) -> (i32, Option<u32>) {
    #[cfg(target_arch = "x86_64")]
    if avx2 && run {
        return copy_avx2_run(
            table,
            dest,
            &held.index,
            &held.cand,
            held.n,
            first,
            block_first,
            totalhits,
            array_size,
        );
    }
    #[cfg(target_arch = "x86_64")]
    if avx2 {
        return copy_avx2(
            table,
            dest,
            &held.index,
            &held.cand,
            held.n,
            first,
            block_first,
            totalhits,
            array_size,
        );
    }
    let _ = (avx2, run);
    copy(
        table,
        dest,
        &held.index,
        &held.cand,
        held.n,
        first,
        block_first,
        totalhits,
        array_size,
    )
}

/// The lookup table's backbone and overflow, for the copy stage.
struct Table {
    backbone: *const BackboneCell,
    overflow: *const i32,
    overflow_len: usize,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:96-107,121-127,55-57
/// The plain-loop candidate stage: writes the index of every word of the block (word `k` has the
/// letters `letters[k..k + 3]`; `letters` holds the block's words and two letters more) to `index`
/// and lists the positions (from the block's first word) of the words whose presence byte is 1;
/// returns how many.
#[inline(never)]
fn candidates(
    letters: &[u8],
    presence: &[u8; X_MASK + 1],
    index: &mut [u16; X_BLOCK],
    cand: &mut [MaybeUninit<u16>; X_BLOCK],
) -> usize {
    let words = letters.len() - 2;
    let index = &mut index[..words];
    let (a, b, c) = (
        &letters[..words],
        &letters[1..words + 1],
        &letters[2..words + 2],
    );
    // The mask is read once per call: as a constant, the vectorized loop's mask vector is built
    // again at every use by V8 (the Wasm build), three instructions in each 8 words. This is a
    // code-generation hint for V8 (`black_box` returns its argument), to check again when the
    // toolchain or Node changes; the indexes do not depend on it.
    let mask = std::hint::black_box(X_MASK as u16);
    for k in 0..words {
        index[k] = ((a[k] as u16) << 10 | (b[k] as u16) << 5 | c[k] as u16) & mask;
    }
    // SAFETY: every index was masked with `mask`, which is `X_MASK`.
    unsafe { list(index, presence, cand) }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:83,90
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:55-57
/// NCBI's `PV_TEST` of every word of a block, without a branch: the positions of the words of
/// `index` whose presence byte is 1, in `cand`, and how many. Each word's position is written at
/// the next free entry, which moves on by the word's presence byte (one add per word; the Wasm
/// build has no scaled addressing to form the entry from a count).
///
/// # Safety
/// Every index of `index` is at most `X_MASK` (both candidate stages write them masked).
#[inline(always)]
unsafe fn list(
    index: &[u16],
    presence: &[u8; X_MASK + 1],
    cand: &mut [MaybeUninit<u16>; X_BLOCK],
) -> usize {
    let start = cand.as_mut_ptr().cast::<u16>();
    let mut out = start;
    for (k, &ix) in index.iter().enumerate() {
        debug_assert!(ix as usize <= X_MASK);
        // SAFETY: (caller) `ix <= X_MASK`, the last index of `presence` (its type has
        // `X_MASK + 1` bytes).
        let bit = unsafe { *presence.get_unchecked(ix as usize) } as usize;
        // SAFETY: every presence byte is 0 or 1 (`XPresence`, made only by `build_presence`), so
        // `out` is `start + n` with `n <= k < index.len() <= X_BLOCK` (`n` = the candidates
        // before word `k`), inside `cand`.
        unsafe {
            out.write(k as u16);
            out = out.add(bit);
        }
    }
    // SAFETY: `out` was moved on from `start` inside `cand`, never back.
    unsafe { out.offset_from_unsigned(start) }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:96-107,121-127,55-57
// (cited at `candidates`)
/// The AVX2 candidate stage: `candidates` with the indexes made 16 at a time (three shifted loads
/// of 16 letters each); the same values.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline(never)]
unsafe fn candidates_avx2(
    letters: &[u8],
    presence: &[u8; X_MASK + 1],
    index: &mut [u16; X_BLOCK],
    cand: &mut [MaybeUninit<u16>; X_BLOCK],
) -> usize {
    use std::arch::x86_64::*;
    let words = letters.len() - 2;
    let l = letters.as_ptr();
    let ixp = index.as_mut_ptr();
    let mask = _mm256_set1_epi16(X_MASK as i16);
    let mut k = 0usize;
    // Words k .. k + 16 read the letters k .. k + 18 (the last load ends at k + 18).
    while k + 16 <= words {
        let a = _mm256_cvtepu8_epi16(_mm_loadu_si128(l.add(k) as *const __m128i));
        let b = _mm256_cvtepu8_epi16(_mm_loadu_si128(l.add(k + 1) as *const __m128i));
        let c = _mm256_cvtepu8_epi16(_mm_loadu_si128(l.add(k + 2) as *const __m128i));
        let v = _mm256_or_si256(
            _mm256_or_si256(_mm256_slli_epi16::<10>(a), _mm256_slli_epi16::<5>(b)),
            c,
        );
        _mm256_storeu_si256(ixp.add(k) as *mut __m256i, _mm256_and_si256(v, mask));
        k += 16;
    }
    while k < words {
        let v = (*l.add(k) as u16) << 10 | (*l.add(k + 1) as u16) << 5 | *l.add(k + 2) as u16;
        *ixp.add(k) = v & X_MASK as u16;
        k += 1;
    }
    // SAFETY: every index was masked with `X_MASK` (the vector loop's `mask`, the tail's `&`).
    list(&index[..words], presence, cand)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:96,100
/// Whether three hits for each of the `left` candidates fit after `totalhits` (then no candidate
/// of up to `AA_HITS_PER_CELL` hits can stop the scan): NCBI's test of the copy stage made once
/// for the cells of up to `AA_HITS_PER_CELL` hits that follow. `0 <= totalhits <= array_size` (a
/// cell is written only after NCBI's test or under this room) and `left <= X_BLOCK`, so the sum is
/// exact in `usize` on every target (in the Wasm build's 32 bits too, which a 64-bit sum turned
/// into 64-bit loop counters).
#[inline(always)]
fn room_for_rest(totalhits: i32, left: usize, array_size: i32) -> bool {
    debug_assert!(0 <= totalhits && totalhits <= array_size && left <= X_BLOCK);
    totalhits as usize + AA_HITS_PER_CELL * left <= array_size as usize
}

/// The plain-loop copy stage for the candidates `cand[..n]` of the block that starts at
/// `block_first` (`copy_block`).
///
/// # Safety
/// As `copy_block`.
#[allow(clippy::too_many_arguments)]
#[inline(never)]
unsafe fn copy(
    table: &Table,
    dest: *mut BlastOffsetPair,
    index: &[u16; X_BLOCK],
    cand: &[MaybeUninit<u16>; X_BLOCK],
    n: usize,
    first: usize,
    block_first: usize,
    totalhits: i32,
    array_size: i32,
) -> (i32, Option<u32>) {
    copy_block::<false, false>(
        table,
        dest,
        index,
        cand,
        n,
        first,
        block_first,
        totalhits,
        array_size,
    )
}

/// The AVX2 copy stage (`copy_block` with the AVX2 writes).
///
/// # Safety
/// As `copy`; AVX2 must be there.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
#[inline(never)]
unsafe fn copy_avx2(
    table: &Table,
    dest: *mut BlastOffsetPair,
    index: &[u16; X_BLOCK],
    cand: &[MaybeUninit<u16>; X_BLOCK],
    n: usize,
    first: usize,
    block_first: usize,
    totalhits: i32,
    array_size: i32,
) -> (i32, Option<u32>) {
    copy_block::<true, false>(
        table,
        dest,
        index,
        cand,
        n,
        first,
        block_first,
        totalhits,
        array_size,
    )
}

/// The AVX2 copy stage with the chain runs (`copy_block` with `write_avx2`'s runs,
/// LOSAT_X_CHAINRUN).
///
/// # Safety
/// As `copy`; AVX2 must be there.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(clippy::too_many_arguments)]
#[inline(never)]
unsafe fn copy_avx2_run(
    table: &Table,
    dest: *mut BlastOffsetPair,
    index: &[u16; X_BLOCK],
    cand: &[MaybeUninit<u16>; X_BLOCK],
    n: usize,
    first: usize,
    block_first: usize,
    totalhits: i32,
    array_size: i32,
) -> (i32, Option<u32>) {
    copy_block::<true, true>(
        table,
        dest,
        index,
        cand,
        n,
        first,
        block_first,
        totalhits,
        array_size,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:91-123
/// The copy stage for the candidates `cand[first..n]` of the block that starts at `block_first`
/// (NCBI's statements above, in NCBI's order): the new total, and the subject offset of the word
/// at which the offset array has no room, if one has none. While three hits for every candidate
/// left fit, a cell of up to `AA_HITS_PER_CELL` hits is written without NCBI's test (after it they
/// still fit); a larger cell takes the test, and the room is computed again after it. A candidate
/// met without that room takes the test whatever its size. `AVX2` = the AVX2 writes (`write`),
/// `RUN` their chain runs (LOSAT_X_CHAINRUN).
///
/// # Safety
/// The candidates' indexes are those of backbone cells of `table`; `dest` has `array_size +
/// X_SCAN_BYTE_SLACK` entries; AVX2 must be there when `AVX2`.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
unsafe fn copy_block<const AVX2: bool, const RUN: bool>(
    table: &Table,
    dest: *mut BlastOffsetPair,
    index: &[u16; X_BLOCK],
    cand: &[MaybeUninit<u16>; X_BLOCK],
    n: usize,
    first: usize,
    block_first: usize,
    mut totalhits: i32,
    array_size: i32,
) -> (i32, Option<u32>) {
    let cand = cand.as_ptr().cast::<u16>();
    let index = index.as_ptr();
    let mut i = first;
    'block: while i < n {
        if room_for_rest(totalhits, n - i, array_size) {
            // The cells of up to `AA_HITS_PER_CELL` hits, up to the next larger cell (`i < n`
            // here). The larger cell leaves the loop with its position, cell and hit count, and
            // its test and copy are after the loop: its values (the overflow array, the array
            // size) are not kept through the loop (the Wasm build kept them in stack slots and
            // reloaded the subject offset from one for every candidate), and its copy does not
            // test the hit count again.
            let (k, cell, numhits) = loop {
                let k = *cand.add(i) as usize;
                let cell = table.backbone.add(*index.add(k) as usize);
                let numhits = (*cell).num_used;
                if numhits as usize > AA_HITS_PER_CELL {
                    break (k, cell, numhits);
                }
                write_cell::<AVX2>(dest.add(totalhits as usize), cell, (block_first + k) as u32);
                totalhits += numhits;
                i += 1;
                if i == n {
                    break 'block;
                }
            };
            let s_off = (block_first + k) as u32;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:96,118-123
            if numhits > array_size - totalhits {
                return (totalhits, Some(s_off));
            }
            let d = dest.add(totalhits as usize);
            write_chain::<AVX2, RUN>(table, d, cell, numhits as usize, s_off);
            totalhits += numhits;
        } else {
            let k = *cand.add(i) as usize;
            let cell = table.backbone.add(*index.add(k) as usize);
            let numhits = (*cell).num_used;
            let s_off = (block_first + k) as u32;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:96,118-123
            // (cited above)
            if numhits > array_size - totalhits {
                return (totalhits, Some(s_off));
            }
            let d = dest.add(totalhits as usize);
            write::<AVX2, RUN>(table, d, cell, numhits as usize, s_off);
            totalhits += numhits;
        }
        i += 1;
    }
    (totalhits, None)
}

/// One offset pair as the 8 bytes it is in memory (`q_off`, then `s_off`: `BlastOffsetPair` is
/// `repr(C)`), to write it with one store.
#[inline(always)]
fn pair_bits(q_off: u32, s_off: u32) -> u64 {
    if cfg!(target_endian = "little") {
        q_off as u64 | (s_off as u64) << 32
    } else {
        (q_off as u64) << 32 | s_off as u64
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:100-105
/// Writes the `numhits` hits of `cell` (its `num_used`; subject offset `s_off`) at `d`: a cell of
/// up to `AA_HITS_PER_CELL` hits by `write_cell`, a larger one by `write_chain`.
///
/// # Safety
/// As `write_chain` for a larger cell, as `write_cell` for the others.
#[inline(always)]
unsafe fn write<const AVX2: bool, const RUN: bool>(
    table: &Table,
    d: *mut BlastOffsetPair,
    cell: *const BackboneCell,
    numhits: usize,
    s_off: u32,
) {
    if numhits <= AA_HITS_PER_CELL {
        write_cell::<AVX2>(d, cell, s_off);
    } else {
        write_chain::<AVX2, RUN>(table, d, cell, numhits, s_off);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:100-102,110-114
/// Writes the hits of `cell`, a cell of up to `AA_HITS_PER_CELL` hits (subject offset `s_off`), at
/// `d` as `AA_HITS_PER_CELL` pairs (with AVX2 `AA_HITS_PER_CELL + 1` in one store), of which the
/// first `num_used` are the hits; the others are not part of the result.
///
/// # Safety
/// `cell` is a backbone cell of up to `AA_HITS_PER_CELL` hits; `d` has room for its hits and
/// `X_SCAN_BYTE_SLACK` entries more; AVX2 must be there when `AVX2`.
#[inline(always)]
unsafe fn write_cell<const AVX2: bool>(
    d: *mut BlastOffsetPair,
    cell: *const BackboneCell,
    s_off: u32,
) {
    #[cfg(target_arch = "x86_64")]
    if AVX2 {
        return write_cell_avx2(d, cell, s_off);
    }
    let d = d.cast::<u64>();
    for j in 0..AA_HITS_PER_CELL {
        d.add(j)
            .write_unaligned(pair_bits((*cell).entries[j] as u32, s_off));
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:103-105,110-114
/// Writes the `numhits` hits of `cell`, a cell of more than `AA_HITS_PER_CELL` (its `num_used`;
/// subject offset `s_off`), at `d` from the overflow array, four pairs at a time. The pairs past
/// the cell's hits are not part of the result; it writes at most three of them (with `RUN` and
/// AVX2, at most `X_CHAIN_RUN - 4` when it has fewer hits than that), and reads past its query
/// offsets only while they stay inside the overflow array (else its last pairs are written one by
/// one).
///
/// # Safety
/// `cell` is a backbone cell of `table` of `numhits > AA_HITS_PER_CELL` hits; `d` has room for
/// them and `X_SCAN_BYTE_SLACK` entries more; AVX2 must be there when `AVX2`.
#[inline(always)]
unsafe fn write_chain<const AVX2: bool, const RUN: bool>(
    table: &Table,
    d: *mut BlastOffsetPair,
    cell: *const BackboneCell,
    numhits: usize,
    s_off: u32,
) {
    debug_assert!(numhits == (*cell).num_used as usize && numhits > AA_HITS_PER_CELL);
    #[cfg(target_arch = "x86_64")]
    if AVX2 {
        return write_chain_avx2::<RUN>(table, d, cell, numhits, s_off);
    }
    let d = d.cast::<u64>();
    let cursor = (*cell).cursor();
    let src = table.overflow.add(cursor);
    let whole = chain_whole(table, cursor, numhits);
    // SAFETY: the first `whole` query offsets are inside the overflow array (`chain_whole`) and `d`
    // has room for `whole <= numhits + 3` pairs; the two do not overlap.
    chain_fours(
        std::slice::from_raw_parts_mut(d.cast::<u8>(), 8 * whole),
        std::slice::from_raw_parts(src, whole),
        s_off,
    );
    for j in whole..numhits {
        d.add(j)
            .write_unaligned(pair_bits(*src.add(j) as u32, s_off));
    }
}

/// Writes the pairs of the query offsets `src` (a multiple of four) with the subject offset
/// `s_off` to `dst` (8 bytes per pair), four at a time. The two slices are separate arguments so
/// that the compiler knows they do not overlap and can widen the four without a runtime test.
#[inline(always)]
fn chain_fours(dst: &mut [u8], src: &[i32], s_off: u32) {
    for (out, q) in dst.chunks_exact_mut(32).zip(src.chunks_exact(4)) {
        for l in 0..4 {
            out[8 * l..8 * l + 8].copy_from_slice(&pair_bits(q[l] as u32, s_off).to_ne_bytes());
        }
    }
}

/// The pairs of an overflow chain of `numhits` hits at `cursor` written four at a time: all of
/// them, rounded up to four, when those query offsets are inside the overflow array, else the
/// whole fours.
#[inline(always)]
fn chain_whole(table: &Table, cursor: usize, numhits: usize) -> usize {
    if cursor + numhits.next_multiple_of(4) <= table.overflow_len {
        numhits.next_multiple_of(4)
    } else {
        numhits & !3
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:100-102,110-114 (cited at
// `write_cell`)
/// `write_cell` with AVX2: one 32-byte store (`[e0, s, e1, s, e2, s, e2, s]` from the cell
/// `[num_used, e0, e1, e2]`).
///
/// # Safety
/// As `write_cell`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn write_cell_avx2(d: *mut BlastOffsetPair, cell: *const BackboneCell, s_off: u32) {
    use std::arch::x86_64::*;
    const _: () = assert!(std::mem::size_of::<BackboneCell>() == 16 && AA_HITS_PER_CELL == 3);
    const _: () =
        assert!(std::mem::size_of::<BlastOffsetPair>() == 8 && cfg!(target_endian = "little"));
    let order = _mm256_setr_epi32(1, 0, 2, 0, 3, 0, 3, 0);
    let c = _mm256_broadcastsi128_si256(_mm_loadu_si128(cell as *const __m128i));
    let s = _mm256_set1_epi32(s_off as i32);
    let v = _mm256_blend_epi32::<0xAA>(_mm256_permutevar8x32_epi32(c, order), s);
    _mm256_storeu_si256(d as *mut __m256i, v);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:103-105,110-114 (cited at
// `write_chain`)
/// `write_chain` with AVX2: 32-byte stores of four query offsets widened to 64 bits with `s_off`
/// in the high halves. With `RUN` (LOSAT_X_CHAINRUN), a chain that has `X_CHAIN_RUN` query
/// offsets inside the overflow array from its start writes its first `X_CHAIN_RUN` pairs without
/// reading its hit count, and loops only past them.
///
/// # Safety
/// As `write_chain`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn write_chain_avx2<const RUN: bool>(
    table: &Table,
    d: *mut BlastOffsetPair,
    cell: *const BackboneCell,
    numhits: usize,
    s_off: u32,
) {
    use std::arch::x86_64::*;
    const _: () =
        assert!(std::mem::size_of::<BlastOffsetPair>() == 8 && cfg!(target_endian = "little"));
    let cursor = (*cell).cursor();
    let src = table.overflow.add(cursor);
    let high = _mm256_set1_epi64x((s_off as i64) << 32);
    let mut j = 0;
    if RUN && cursor + X_CHAIN_RUN <= table.overflow_len {
        // A fixed count: the loop is unrolled and has no branch on `numhits`.
        while j < X_CHAIN_RUN {
            let q = _mm256_cvtepu32_epi64(_mm_loadu_si128(src.add(j) as *const __m128i));
            _mm256_storeu_si256(d.add(j) as *mut __m256i, _mm256_or_si256(q, high));
            j += 4;
        }
        if numhits <= X_CHAIN_RUN {
            return;
        }
    }
    let whole = chain_whole(table, cursor, numhits);
    while j < whole {
        let q = _mm256_cvtepu32_epi64(_mm_loadu_si128(src.add(j) as *const __m128i));
        _mm256_storeu_si256(d.add(j) as *mut __m256i, _mm256_or_si256(q, high));
        j += 4;
    }
    while j < numhits {
        *d.add(j) = BlastOffsetPair {
            q_off: *src.add(j) as u32,
            s_off,
        };
        j += 1;
    }
}

#[cfg(test)]
#[path = "x_scan_byte_tests.rs"]
mod tests;
