//! TBLASTX's two-hit ungapped extension
//! (`two_hit.rs::extend_hit_two_hit_fitting`: the word loop and the `s_BlastAaExtendLeft` /
//! `s_BlastAaExtendRight` walks of NCBI `s_BlastAaExtendTwoHit`) with the scores read from one
//! flat byte table.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:846-866
//! (right walk; the left walk, aa_ungapped.c:886-921, counts `i` down from `n`, has no
//! `score <= 0` test, and returns `n - best_i + 1`.)
//!
//! When these walks were written, the walks of `two_hit.rs` called `get_score` for every step,
//! which tests both codes against the matrix size, branches to the sentinel and the rare codes,
//! and multiplies for the row: 22 instructions and five branches a step on the left, about 30 and
//! six on the right (whose bounds are reloaded from the stack; release build of dev 7cdcd4a7).
//! Since line RG item X2c-s they read the matrix's rows (`TbxScore::score`, `ScoreRows`: a row and
//! its score, two loads), for the X-drops and the one-hit extension that these walks do not
//! cover: 16 to 19 instructions a step (the walks inlined into `word_finder` and `x_tbx_qslice`,
//! release builds of X2c-s). A step here loads the two codes, forms `q << 8 | s`, loads one signed
//! byte from the table and makes the reference's update and exit tests; the statements, their
//! order and the exit tests are those of the reference walks. Each walk is a function of its own, so that its loop has the registers to
//! itself (inlined into the caller, the X-drop and the sequence pointers were reloaded from the
//! stack at every step).
//!
//! wasm32 (plan v4 line W, s17): V8's code for the walks reads the table at a fixed address
//! (`FIXED_TABLE`, one add fewer a step), the left walk hands its pair back as one u64 and the right
//! walk is inlined into `two_hit`. A Wasm call keeps no register across it (every live value goes
//! to the stack and back) and a returned pair or triple goes through linear memory; inlined, the
//! right walk's loop keeps its 22 instructions a step, while an inlined left walk lost more cycles
//! than the call cost. The statements of the walks are the same on every target.
//!
//! Why the result is the same.
//! * Bounds. The reference reads `q_seq[q_r - 1 - i]` for `i < min(q_r, s_r)` on the left and
//!   `q_seq[q_r + j]` while both `q_r + j < q_seq.len()` and `s_r + j < s_seq.len()` on the right
//!   (NCBI's walks run on the query batch and read the NULLB past the context's ends, which ends a
//!   walk of an X-drop of 127 or less there: the reference's continuation into the batch,
//!   query_batch.rs, changes nothing for these X-drops).
//!   The walks here take the same counts and read the same positions, so every position is inside
//!   both slices; the word loop reads the word's positions (3 for BLOSUM62's `two_hit`, `W` for
//!   `two_hit_in_word::<W>`) after a test that they fit (the reference reads them unchecked; a
//!   word that does not fit runs the reference).
//! * Table. `T[q << 8 | s]` is `get_score(q, s)` (`blosum62_score_ncbistdaa_direct`, defined for
//!   every pair of bytes) for every two bytes, with `BLAST_SCORE_MIN` (-32768, the gap and
//!   sentinel code 0 against anything) stored as -128; the table exists only if every other score
//!   is within -127..=11 (`MAX_WORD_SCORE`, the sentinel argument below). Only the rows of the codes present are touched (28 rows of 256 bytes,
//!   and of each row the first cache line).
//! * Sentinel. Take `1 <= x_drop <= 127` (otherwise the reference runs). The score before a step
//!   is never above the best, so a step of -128 makes `best - score >= 128 > x_drop` as -32768
//!   does: both walks stop at that step, it does not raise the best in either, and the best, its
//!   length and the position of the stop are the same. The word loop adds three scores without a
//!   stop and only compares the sums with a best that starts at 0: a -128 is not made positive by
//!   the two other scores (at most `MAX_WORD_SCORE` = 11 each, which the table build checks), as
//!   -32768 is not, so every comparison has the same outcome.
//! * Overflow. Both paths add the same i32 scores in the same order, so they agree up to an
//!   overflow; a walk would need a best near 2^31 (some 195 M residues of positive diagonal) for a
//!   -128 and a -32768 step to wrap differently, where the reference (and NCBI's `Int4`) overflows
//!   too.
//! * Other matrices (`-matrix`, SC7b) and other word lengths (`-word_size` 2 and 4, any
//!   matrix). `two_hit_in_word::<W>` (`two_hit_in`, `two_hit_in2`, `two_hit_in4` for `W` = 3, 2,
//!   4) runs the same word loop over `W` letters (`TbxScore::word`, the reference's `k_size`) and
//!   the same walks on a byte table of the search's matrix (`matrix_walk_table`: a static of each
//!   matrix, built at the first search with it from the scores the reference walks read,
//!   `TbxScore::score`), for the calls `two_hit` covers. The table holds the 32 x 32 cells of the
//!   bytes below 32 at `q << 5 | s` (the frames hold the NCBISTDAA letters 0..=27, the sentinel 0
//!   among them); its 256 rows keep that index inside it for any two bytes. The arguments above
//!   hold with -127..=`MAX_MATRIX_SCORE` (42) in place of -127..=11 and `W` letters in place of 3:
//!   a word's sum of a -128 and `W - 1` (at most 3) scores of at most 42 stays at or below 0
//!   (3 x 42 = 126), as with -32768. NCBI's tables score -17..=17 (BLOSUM62
//!   -4..=11). wasm32's fixed table stays BLOSUM62's at 3 letters.
//!
//! The table walks act on every target (the former switch `LOSAT_X_TBXUNGAP` graduated on
//! 2026-10-26). Until line W of plan v4 measured them (item W7), wasm32 kept the 16-residue v128
//! walks of `two_hit.rs`, which take each lane out of the vector and call `get_score` for it
//! (scalar steps); with the table the two-hit extension of the LOSAT Web engine takes 0.77-0.89
//! of that time (`docs/agents/diagnostics.md`, entry `LOSAT_X_TBXUNGAP`). The calls the table
//! does not cover (the Why section's bounds and X-drop) run `extend_hit_two_hit_reference`.
//!
//! LOSAT_X_TBXUNGAPSENT (native builds; in the strict set, `=0` turns it off): `two_hit_sent`, the
//! same extension with walks that the sentinels stop instead of a bound, inlined with the word loop
//! into the callers' extension code (`extend_hit_two_hit_sent`: the unit path's `x_hit` and the
//! batched `resolve`, the reference frame loops).
//! * Instances. The unit path chooses for each unit of its branch-free hit loops
//!   (LOSAT_X_TBXHITBF; `XPairPar::x_unit_sent`, in the dispatches `scan_unit_or_bf`,
//!   `scan_unit_batched_or_bf` and `scan_unit_sliced_or_bf` of x_tbx_hitbf.rs) between instances
//!   of the loop's function (`const SENT`: `scan_unit_bf`, `scan_unit_batched_bf`,
//!   `scan_unit_sliced_bf`): without the walks (0), so that it keeps the machine code it had
//!   before the switch, and with the walks of the search's word length (3, or 2 for the other
//!   matrices' instance; `TbxScore::sent_two_hit_of`). The units of the other hit loop
//!   (`LOSAT_X_TBXHITBF=0`, the probes of `LOSAT_TIMING` and `LOSAT_X_STATS`) take `two_hit`.
//!   The reference frame loops (`blast_engine/word_finder.rs`: `begin_chunk_scan` sets
//!   `subject_sent` per chunk, `two_hit_one_hit` reads it) choose at every call.
//! * Where. A query context's buffer is `[0, residues.., 0]` and a subject frame's too
//!   (`translation.rs`, NCBI's NULLB at each end, blast_util.c:438-453); the callers' sequences are
//!   the residues. A subject frame of more than `MAX_DBSEQ_LEN` (5,000,000) residues is searched in
//!   overlapping chunks, and a chunk's ends are residues of the frame, so only a chunk that is the
//!   whole frame takes the path (`sent_subject`, per unit or chunk); the contexts' sentinels are
//!   checked once per batch (`sent_contexts`). The others take `two_hit`; the debug output of the
//!   reference frame loops keeps the reference.
//! * Sentinel. The reference's left walk takes `min(q_r, s_r)` steps and its right walk the room
//!   left in both sequences; the position one step further is just outside the sequence that set
//!   the bound, which holds a sentinel. A step on it scores -128 and stops the walk as the bound
//!   does, with the same best, length and stop position (the argument of the sentinel above: the
//!   step does not raise the best and makes `best - score >= 128 > x_drop`). Every step before it
//!   reads the positions the reference reads. The word loop reads the same positions as in
//!   `two_hit` (or `two_hit_in_word::<W>`, the word's `W`), after the same test that they fit.
//! * Overflow. A walk without a bound would step over a sentinel if its score wrapped. The left walk
//!   (from 0) and then the right walk (from the left walk's best) read each subject position once at
//!   most, at +`MAX_MATRIX_SCORE` (42) at most, so a score stays below `42 * (subject letters + 2)`, and
//!   it is never below `best - 255` (a step costs 128 at most, and the walk stops once `best - score`
//!   reaches `x_drop <= 127`); the path takes frames of at most `MAX_SENT_LEN` (`i32::MAX / 64`)
//!   bytes, above any whole frame (a longer one is searched in chunks of `MAX_DBSEQ_LEN`).
//! * Keys. The walks are the key walks of blastp and BLASTX (algorithm/common/x_key_walk.rs, item
//!   SC2e: score and step in one `i64`, so that one comparison keeps the best and its length), and
//!   the word loop runs on the same keys: a step adds `(t << 32) - 1` for the score `t` of the two
//!   bytes, read at `q << 5 | s` (`KeyTable`: BLOSUM62's from the cells of `build_table`, another
//!   matrix's beside its `MatrixWalkTable`). The keys hold the byte tables' scores of the bytes
//!   below `MATRIX_LETTERS`: the letters (NCBISTDAA 0..=27) and the sentinel 0, the only bytes a
//!   frame or a context holds (`translation.rs`), as `MatrixWalkTable` takes them; a context or a
//!   frame takes the path only if every byte of its buffer is below it (`letters` in
//!   `sent_contexts` and `sent_subject`), so every step reads a key of the table. A walk takes
//!   fewer than `MAX_SENT_LEN` steps, inside the keys' step part (`2^32 - 1`), and the scores stay
//!   as in the item above.
//! * Other matrices and word lengths (`TbxScore::sent_two_hit`, `two_hit_sent_in::<W>`): the same
//!   walks on the keys of the matrix's byte table after the word loop over `W` letters, for a table
//!   whose cells of code 0 all hold -128 (`MatrixWalkTable::x_sentinels_stop`); the arguments above
//!   hold with its scores, at most `MAX_MATRIX_SCORE`.
//!
//! No NCBI counterpart for the split itself and the switch: they choose between LOSAT code paths
//! with the same result.

use super::get_score;
#[cfg(test)]
use super::two_hit::extend_hit_two_hit_reference;
use crate::algorithm::common::x_key_walk::{best_of, walk, KeyTable, KEY_TABLE_LEN};
use crate::config::ScoringMatrix;
use crate::utils::matrix::BLAST_SCORE_MIN;
use std::cell::UnsafeCell;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Once, OnceLock};

/// The seven values of `extend_hit_two_hit_fitting`: (q_start, q_end, s_start, s_end, score,
/// right_extended, s_last_off).
pub(crate) type TwoHit = (usize, usize, usize, usize, i32, bool, usize);

/// `q << 8 | s` for every two bytes.
const TABLE_LEN: usize = 1 << 16;
/// Largest X-drop the byte table is exact for (module comment).
const MAX_X_DROP: i32 = 127;
/// Residues of BLOSUM62's word loop (`k_size` of `two_hit.rs`, `TbxScore::word` of
/// `Blosum62Score`).
const WORD: usize = 3;
/// The largest matrix score the word-loop argument allows (module comment).
const MAX_WORD_SCORE: i32 = 11;
/// The longest word of the walks on a matrix's table (`two_hit_in_word::<W>`, `two_hit_sent_in::<W>`):
/// NCBI's TBLASTX takes `-word_size` 2..=4 (`TbxScore::word`).
const MAX_MATRIX_WORD: usize = 4;
/// The largest score of a matrix's table (`matrix_walk_table`): the word loop's sums of a -128 and
/// up to `MAX_MATRIX_WORD - 1` scores stay at or below 0 while `(W - 1) * score <= 128`, as those
/// of -32768 do.
const MAX_MATRIX_SCORE: i32 = 42;
const _: () = assert!((MAX_MATRIX_WORD as i32 - 1) * MAX_MATRIX_SCORE <= 128);
/// The bytes another matrix's table holds the scores of (`MatrixWalkTable`): the frames' letters
/// (NCBISTDAA 0..=27) are below it.
const MATRIX_LETTERS: usize = 32;
const _: () = assert!(MATRIX_LETTERS * MATRIX_LETTERS == KEY_TABLE_LEN);
/// `q << 5 | s` for every two bytes (`MatrixWalkTable`).
const MATRIX_TABLE_LEN: usize = 256 * MATRIX_LETTERS;
#[expect(
    clippy::ineffective_bit_mask,
    reason = "the assertion spells the largest index `q << 5 | s` of two bytes as the walks compute it"
)]
const _: () = assert!((255 << 5 | 255) < MATRIX_TABLE_LEN);

type Table = [i8; TABLE_LEN];

/// wasm32: a copy of the table at a fixed address of linear memory, so that a load of the walks
/// carries the address as its offset instead of adding the table pointer at every step (the native
/// walks address `[table + index]` in one instruction). Zero until `state` copies the table in.
#[cfg(target_arch = "wasm32")]
struct FixedTable(std::cell::UnsafeCell<Table>);
// SAFETY: written once, inside `STATE.get_or_init`, before `State::table` points at it; read only
// through that reference.
#[cfg(target_arch = "wasm32")]
unsafe impl Sync for FixedTable {}
#[cfg(target_arch = "wasm32")]
static FIXED_TABLE: FixedTable = FixedTable(std::cell::UnsafeCell::new([0; TABLE_LEN]));

struct State {
    /// LOSAT_X_TBXUNGAPSENT (native builds, with the table): the walks the sentinels stop.
    sent: bool,
    /// The table of the walks; `None` (the reference runs) if a score does not fit it.
    #[cfg(not(target_arch = "wasm32"))]
    table: Option<Box<Table>>,
    /// The table of the walks (`FIXED_TABLE`); `None` (the reference runs) if a score does not fit it.
    #[cfg(target_arch = "wasm32")]
    table: Option<&'static Table>,
    /// LOSAT_X_TBXUNGAPSENT: the keys of `table`'s letters for the sentinel walks (`fill_keys`);
    /// `None` when the path is off.
    keys: Option<Box<KeyTable>>,
}

/// `tab[i]`. wasm32 reads `FIXED_TABLE` (where `tab` points) at its static address.
///
/// SAFETY: `i < TABLE_LEN`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:846-921 `matrix[q[i]][s[i]]` (the table of `build_table`).
#[inline(always)]
unsafe fn score_at(
    #[cfg_attr(
        target_arch = "wasm32",
        expect(
            unused_variables,
            reason = "wasm32 reads `FIXED_TABLE` at its static address"
        )
    )]
    tab: &Table,
    i: usize,
) -> i32 {
    #[cfg(target_arch = "wasm32")]
    let tab: &Table = &*FIXED_TABLE.0.get();
    *tab.get_unchecked(i) as i32
}

/// The scores the walks read: BLOSUM62's `Table` or another matrix's `MatrixWalkTable`.
trait WalkTable {
    /// The score of the query byte `q` against the subject byte `s`.
    fn at(&self, q: u8, s: u8) -> i32;
}

impl WalkTable for Table {
    #[inline(always)]
    fn at(&self, q: u8, s: u8) -> i32 {
        // SAFETY: q << 8 | s < TABLE_LEN for any two bytes.
        unsafe { score_at(self, (q as usize) << 8 | s as usize) }
    }
}

/// The byte table of a matrix other than BLOSUM62 (`two_hit_in`): the score of `q` against `s` at
/// `q << 5 | s` for the bytes below `MATRIX_LETTERS`. A byte of 32 or more, which no frame holds,
/// reads a cell that is not its score. Page-aligned: the cells of the letters (the first 1 KiB),
/// the only ones written (`fill_matrix_table`), take one page of memory. Beside it, the keys of
/// the sentinel walks (LOSAT_X_TBXUNGAPSENT, `fill_keys`; 8 KiB).
#[repr(align(4096))]
pub(crate) struct MatrixWalkTable([i8; MATRIX_TABLE_LEN], KeyTable);

/// LOSAT_X_TBXUNGAPSENT: writes into `k` the keys of the sentinel walks (algorithm/common/
/// x_key_walk.rs, item SC2e) for the letters of a byte table `t` read at `q << SHIFT | s`
/// (BLOSUM62's `Table`: 8; a `MatrixWalkTable`: 5): `(t << 32) - 1` at `q << 5 | s` for the score
/// `t` of every two bytes below `MATRIX_LETTERS`, the only bytes of a buffer that takes the walks
/// (`letters`).
// No NCBI counterpart: the scores the reference walks add, in the form of the key walks.
fn fill_keys<const SHIFT: u32>(k: &mut KeyTable, t: &[i8]) {
    for q in 0..MATRIX_LETTERS {
        for s in 0..MATRIX_LETTERS {
            k[q << 5 | s] = ((t[q << SHIFT | s] as i64) << 32) - 1;
        }
    }
}

impl MatrixWalkTable {
    /// Code 0 against anything reads -128, so that the sentinel walks stop at a sentinel.
    pub(crate) fn x_sentinels_stop(&self) -> bool {
        (0..MATRIX_LETTERS).all(|c| self.0[c] == -128 && self.0[c << 5] == -128)
    }
}

impl WalkTable for MatrixWalkTable {
    #[inline(always)]
    fn at(&self, q: u8, s: u8) -> i32 {
        debug_assert!(
            (q as usize) < MATRIX_LETTERS && (s as usize) < MATRIX_LETTERS,
            "a byte of the frames is a letter: {q} {s}"
        );
        // SAFETY: q << 5 | s < MATRIX_TABLE_LEN for any two bytes.
        unsafe { *self.0.get_unchecked((q as usize) << 5 | s as usize) as i32 }
    }
}

/// `walk_left`'s (best score, length). wasm32 hands them back as one u64 (the score's bits in the
/// low half, the length in the high half): a returned pair goes through linear memory there.
#[cfg(not(target_arch = "wasm32"))]
type LeftWalk = (i32, usize);
#[cfg(target_arch = "wasm32")]
type LeftWalk = u64;

// No NCBI counterpart: carries two values of `s_BlastAaExtendLeft` (aa_ungapped.c:886-921).
#[inline(always)]
fn left_walk(best: i32, len: usize) -> LeftWalk {
    #[cfg(not(target_arch = "wasm32"))]
    return (best, len);
    #[cfg(target_arch = "wasm32")]
    return best as u32 as u64 | (len as u64) << 32;
}

// No NCBI counterpart: the inverse of `left_walk`.
#[inline(always)]
fn left_walk_parts(w: LeftWalk) -> (i32, usize) {
    #[cfg(not(target_arch = "wasm32"))]
    return w;
    #[cfg(target_arch = "wasm32")]
    return (w as u32 as i32, (w >> 32) as usize);
}

/// The table of the walks, `None` if a score other than `BLAST_SCORE_MIN` is outside
/// -127..=`MAX_WORD_SCORE`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:846-921 `matrix[q[i]][s[i]]`; the scores are those of the function
// the reference walks call.
fn build_table() -> Option<Box<Table>> {
    let mut t: Box<Table> = vec![0i8; TABLE_LEN].into_boxed_slice().try_into().ok()?;
    for q in 0..=255u8 {
        for s in 0..=255u8 {
            let v = get_score(q, s);
            t[(q as usize) << 8 | s as usize] = if v == BLAST_SCORE_MIN {
                -128
            } else if (-127..=MAX_WORD_SCORE).contains(&v) {
                v as i8
            } else {
                return None;
            };
        }
    }
    Some(t)
}

/// Writes the scores `score(q, s)` of the bytes below `MATRIX_LETTERS` into another matrix's
/// table `t`; `false` if one other than `BLAST_SCORE_MIN` is outside -127..=`MAX_MATRIX_SCORE`.
fn fill_matrix_table(t: &mut [i8; MATRIX_TABLE_LEN], score: impl Fn(u8, u8) -> i32) -> bool {
    for q in 0..MATRIX_LETTERS as u8 {
        for s in 0..MATRIX_LETTERS as u8 {
            let v = score(q, s);
            t[(q as usize) << 5 | s as usize] = if v == BLAST_SCORE_MIN {
                -128
            } else if (-127..=MAX_MATRIX_SCORE).contains(&v) {
                v as i8
            } else {
                return false;
            };
        }
    }
    true
}

/// The state of one matrix's table of `matrix_walk_table` (`MATRIX_TABLES`).
struct MatrixTableSlot {
    /// Writes the table, at the first call for the matrix.
    once: Once,
    /// Whether the matrix's scores fit the table (`fill_matrix_table`), set inside `once`.
    fits: AtomicBool,
}

/// The tables of `matrix_walk_table`, one for each matrix, zero until a slot's `once` writes the
/// cells of the letters (the first 1 KiB). Apart from the slots: a `Once` does not start at zero,
/// and a static with one (or a `OnceLock`, with its uninitialized bytes) lands in .data, 72 KiB
/// of the binary; all zero, the tables stay in .bss.
struct MatrixTables([UnsafeCell<MatrixWalkTable>; ScoringMatrix::ALL.len()]);
// SAFETY: a table is written once, inside its slot's `once`, and read only through the references
// `matrix_walk_table` hands out after that `once` completed.
unsafe impl Sync for MatrixTables {}

static MATRIX_SLOTS: [MatrixTableSlot; ScoringMatrix::ALL.len()] = [const {
    MatrixTableSlot {
        once: Once::new(),
        fits: AtomicBool::new(false),
    }
}; ScoringMatrix::ALL.len()];
static MATRIX_TABLES: MatrixTables = MatrixTables(
    [const { UnsafeCell::new(MatrixWalkTable([0; MATRIX_TABLE_LEN], [0; KEY_TABLE_LEN])) };
        ScoringMatrix::ALL.len()],
);

/// The byte table of `matrix`, a matrix other than BLOSUM62 whose scores are `score(q, s)`: a
/// static written at the first call for the matrix, or `None` (the reference walks run) if a score
/// other than `BLAST_SCORE_MIN` is outside -127..=`MAX_MATRIX_SCORE`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:846-921 `matrix[q[i]][s[i]]`; the scores are those the reference
// walks read (`TbxScore::score` of the search's matrix).
pub(crate) fn matrix_walk_table(
    matrix: ScoringMatrix,
    score: impl Fn(u8, u8) -> i32,
) -> Option<&'static MatrixWalkTable> {
    let slot = &MATRIX_SLOTS[matrix as usize];
    let table = &MATRIX_TABLES.0[matrix as usize];
    slot.once.call_once(|| {
        // SAFETY: the only write to the table, inside `once`, before any reference to it
        // (MatrixTables).
        let t = unsafe { &mut *table.get() };
        let fits = fill_matrix_table(&mut t.0, score);
        fill_keys::<5>(&mut t.1, &t.0);
        // `once` orders this store before the loads that follow its completion.
        slot.fits.store(fits, Ordering::Relaxed);
    });
    // SAFETY: `once` has completed, so nothing writes the table any more (MatrixTables).
    slot.fits
        .load(Ordering::Relaxed)
        .then(|| unsafe { &*table.get() })
}

static STATE: OnceLock<State> = OnceLock::new();

fn state() -> &'static State {
    STATE.get_or_init(|| {
        let table = build_table();
        // SAFETY: the only write to FIXED_TABLE (FixedTable).
        #[cfg(target_arch = "wasm32")]
        let table = table.map(|t| unsafe {
            *FIXED_TABLE.0.get() = *t;
            &*FIXED_TABLE.0.get()
        });
        // LOSAT_X_TBXUNGAPSENT: a strict-set switch on the table walks; the wasm32 builds keep the
        // walks of line W.
        let sent = table.is_some()
            && cfg!(not(target_arch = "wasm32"))
            && crate::utils::xswitch::on("LOSAT_X_TBXUNGAPSENT");
        let keys = table.as_ref().filter(|_| sent).map(|t| {
            let mut k = Box::new([0i64; KEY_TABLE_LEN]);
            fill_keys::<8>(&mut k, &t[..]);
            k
        });
        State { sent, table, keys }
    })
}

/// `s_BlastAaExtendLeft` ending at `(qp, sp)` over `n` steps: step `i` (1..=n) reads `qp[-i]`,
/// `sp[-i]`. Returns the best score and its length.
///
/// SAFETY: the `n` positions before `qp` and before `sp` are readable.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:886-921 (module comment). Step `i` reads index `n - i` of the
// `n` positions before the pointers, so one counter runs the loop and the length is
// `n - (the index of the best step)`.
#[inline(never)]
unsafe fn walk_left<T: WalkTable>(
    tab: &T,
    qp: *const u8,
    sp: *const u8,
    n: usize,
    x_drop: i32,
) -> LeftWalk {
    let q_lo = qp.sub(n);
    let s_lo = sp.sub(n);
    let mut score = 0i32;
    let mut best = 0i32;
    let mut best_at = n;
    let mut at = n;
    while at > 0 {
        at -= 1;
        let qc = *q_lo.add(at);
        let sc = *s_lo.add(at);
        score += tab.at(qc, sc);
        if score > best {
            best = score;
            best_at = at;
        }
        if best - score >= x_drop {
            break;
        }
    }
    left_walk(best, n - best_at)
}

/// `s_BlastAaExtendRight` from `(qp, sp)` over `n` steps, from the score `start`: step `j` (0..n)
/// reads `qp[j]`, `sp[j]`. Returns the best score, its length and the step at which the loop left
/// (`n` if it ran out).
///
/// SAFETY: the `n` positions from `qp` and from `sp` are readable.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:846-866 (module comment). Inlined on wasm32 (module comment).
#[cfg_attr(not(target_arch = "wasm32"), inline(never))]
#[cfg_attr(target_arch = "wasm32", inline(always))]
unsafe fn walk_right<T: WalkTable>(
    tab: &T,
    qp: *const u8,
    sp: *const u8,
    n: usize,
    x_drop: i32,
    start: i32,
) -> (i32, usize, usize) {
    let mut score = start;
    let mut best = start;
    let mut best_len = 0usize;
    let mut taken = 0usize;
    while taken < n {
        let qc = *qp.add(taken);
        let sc = *sp.add(taken);
        score += tab.at(qc, sc);
        taken += 1;
        if score > best {
            best = score;
            best_len = taken;
        }
        if score <= 0 || best - score >= x_drop {
            return (best, best_len, taken - 1);
        }
    }
    (best, best_len, n)
}

/// `extend_hit_two_hit_fitting(q_seq, s_seq, s_left_off, s_right_off, q_right_off, x_drop, false)`, or
/// `None` when the reference has to run (an X-drop outside 1..=127, a word that does not fit).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit` (two_hit.rs follows it).
#[inline(never)]
pub(super) fn two_hit(
    q_seq: &[u8],
    s_seq: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> Option<TwoHit> {
    if !(1..=MAX_X_DROP).contains(&x_drop)
        || q_right_off >= q_seq.len()
        || q_seq.len() - q_right_off < WORD
        || s_right_off >= s_seq.len()
        || s_seq.len() - s_right_off < WORD
    {
        return None;
    }
    #[cfg(not(target_arch = "wasm32"))]
    let tab: &Table = state().table.as_deref()?;
    // SAFETY: `state().table` is set, so FIXED_TABLE holds the table (FixedTable). wasm32: the word
    // loop below reads it at its static address, as `score_at` does.
    #[cfg(target_arch = "wasm32")]
    let tab: &Table = {
        state().table?;
        unsafe { &*FIXED_TABLE.0.get() }
    };

    // The best end of the word at R (two_hit.rs, the `k_size` loop).
    let mut score = 0i32;
    let mut seed_score = 0i32;
    let mut right_d = 0usize;
    for i in 0..WORD {
        // SAFETY: the word fits in both sequences (tested above).
        let (qc, sc) = unsafe {
            (
                *q_seq.get_unchecked(q_right_off + i),
                *s_seq.get_unchecked(s_right_off + i),
            )
        };
        score += tab[(qc as usize) << 8 | sc as usize] as i32;
        if score > seed_score {
            seed_score = score;
            right_d = i + 1;
        }
    }
    let q_r = q_right_off + right_d;
    let s_r = s_right_off + right_d;

    // SAFETY: q_r <= q_seq.len() and s_r <= s_seq.len() (the word fits); the left walk reads the
    // min(q_r, s_r) positions before them, the right walk the room left after them in both.
    let (left_score, left_disp) = left_walk_parts(unsafe {
        walk_left(
            tab,
            q_seq.as_ptr().add(q_r),
            s_seq.as_ptr().add(s_r),
            q_r.min(s_r),
            x_drop,
        )
    });

    let reached_first_hit = left_disp >= s_r.saturating_sub(s_left_off);
    let (right_score, right_disp, s_last_off) = if reached_first_hit {
        let n = (q_seq.len() - q_r).min(s_seq.len() - s_r);
        // SAFETY: see above.
        let (best, len, end) = unsafe {
            walk_right(
                tab,
                q_seq.as_ptr().add(q_r),
                s_seq.as_ptr().add(s_r),
                n,
                x_drop,
                left_score,
            )
        };
        (best, len, s_r + end)
    } else {
        (0, 0, s_r)
    };

    Some((
        q_r - left_disp,
        q_r + right_disp,
        s_r - left_disp,
        s_r + right_disp,
        left_score.max(right_score),
        reached_first_hit,
        s_last_off,
    ))
}

/// `two_hit` with the byte table `t` of a matrix and words of 3 letters: the values of
/// `extend_hit_two_hit_reference` with that matrix's scores and `TbxScore::word` 3, or `None` when
/// the reference has to run. [`two_hit_in2`] and [`two_hit_in4`] are the same for words of 2 and 4
/// letters (`-word_size`).
///
/// Three functions around one generic body (`two_hit_in_word::<W>`) rather than one generic
/// function: `TbxScore::table_two_hit`, inlined into the callers, calls them, and a generic
/// function there makes the walks it calls (`walk_left`, `walk_right`) exported symbols of the
/// library, which BLOSUM62's `two_hit` then calls through the GOT (no-LTO build); out of line and
/// not generic, the walks stay local as before.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit` (two_hit.rs follows it).
#[inline(never)]
pub(crate) fn two_hit_in(
    t: &MatrixWalkTable,
    q_seq: &[u8],
    s_seq: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> Option<TwoHit> {
    two_hit_in_word::<3>(
        t,
        q_seq,
        s_seq,
        s_left_off,
        s_right_off,
        q_right_off,
        x_drop,
    )
}

/// [`two_hit_in`] with words of 2 letters (`-word_size 2`).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit` (`word_size` 2).
#[inline(never)]
pub(crate) fn two_hit_in2(
    t: &MatrixWalkTable,
    q_seq: &[u8],
    s_seq: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> Option<TwoHit> {
    two_hit_in_word::<2>(
        t,
        q_seq,
        s_seq,
        s_left_off,
        s_right_off,
        q_right_off,
        x_drop,
    )
}

/// [`two_hit_in`] with words of 4 letters (`-word_size 4`).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit` (`word_size` 4).
#[inline(never)]
pub(crate) fn two_hit_in4(
    t: &MatrixWalkTable,
    q_seq: &[u8],
    s_seq: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> Option<TwoHit> {
    two_hit_in_word::<4>(
        t,
        q_seq,
        s_seq,
        s_left_off,
        s_right_off,
        q_right_off,
        x_drop,
    )
}

/// The body of [`two_hit_in`], [`two_hit_in2`] and [`two_hit_in4`]: the table walks with words of
/// `W` letters, or `None` when the reference has to run.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit` (two_hit.rs follows it; `W` is
// its `word_size`).
#[inline(always)]
fn two_hit_in_word<const W: usize>(
    t: &MatrixWalkTable,
    q_seq: &[u8],
    s_seq: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> Option<TwoHit> {
    if !fits::<W>(q_seq, s_seq, s_right_off, q_right_off, x_drop) {
        return None;
    }
    Some(two_hit_on::<W>(
        t,
        q_seq,
        s_seq,
        s_left_off,
        s_right_off,
        q_right_off,
        x_drop,
    ))
}

/// Whether the table walks cover the call: an X-drop in 1..=127 and a word of `W` letters that
/// fits (the test of `two_hit`).
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
#[inline(always)]
fn fits<const W: usize>(
    q_seq: &[u8],
    s_seq: &[u8],
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> bool {
    (1..=MAX_X_DROP).contains(&x_drop)
        && q_right_off < q_seq.len()
        && q_seq.len() - q_right_off >= W
        && s_right_off < s_seq.len()
        && s_seq.len() - s_right_off >= W
}

/// The word loop (over `W` letters) and the two walks of `two_hit` (the same statements) on a
/// matrix's table `tab`. `two_hit` keeps its own copy, so that BLOSUM62's machine code stays as it
/// was; the tests compare both with the reference walks.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit` (two_hit.rs follows it).
#[inline(always)]
fn two_hit_on<const W: usize>(
    tab: &MatrixWalkTable,
    q_seq: &[u8],
    s_seq: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> TwoHit {
    // The word loop's argument (module comment) holds for `W` letters.
    const { assert!(W >= 1 && W <= MAX_MATRIX_WORD) };
    // The best end of the word at R (two_hit.rs, the `k_size` loop).
    let mut score = 0i32;
    let mut seed_score = 0i32;
    let mut right_d = 0usize;
    for i in 0..W {
        // SAFETY: the word fits in both sequences (tested above).
        let (qc, sc) = unsafe {
            (
                *q_seq.get_unchecked(q_right_off + i),
                *s_seq.get_unchecked(s_right_off + i),
            )
        };
        score += tab.at(qc, sc);
        if score > seed_score {
            seed_score = score;
            right_d = i + 1;
        }
    }
    let q_r = q_right_off + right_d;
    let s_r = s_right_off + right_d;

    // SAFETY: q_r <= q_seq.len() and s_r <= s_seq.len() (the word fits); the left walk reads the
    // min(q_r, s_r) positions before them, the right walk the room left after them in both.
    let (left_score, left_disp) = left_walk_parts(unsafe {
        walk_left(
            tab,
            q_seq.as_ptr().add(q_r),
            s_seq.as_ptr().add(s_r),
            q_r.min(s_r),
            x_drop,
        )
    });

    let reached_first_hit = left_disp >= s_r.saturating_sub(s_left_off);
    let (right_score, right_disp, s_last_off) = if reached_first_hit {
        let n = (q_seq.len() - q_r).min(s_seq.len() - s_r);
        // SAFETY: see above.
        let (best, len, end) = unsafe {
            walk_right(
                tab,
                q_seq.as_ptr().add(q_r),
                s_seq.as_ptr().add(s_r),
                n,
                x_drop,
                left_score,
            )
        };
        (best, len, s_r + end)
    } else {
        (0, 0, s_r)
    };

    (
        q_r - left_disp,
        q_r + right_disp,
        s_r - left_disp,
        s_r + right_disp,
        left_score.max(right_score),
        reached_first_hit,
        s_last_off,
    )
}

/// Longest subject frame of the sentinel walks (`[0, frame.., 0]` bytes): no score wraps (module
/// comment, section LOSAT_X_TBXUNGAPSENT).
const MAX_SENT_LEN: usize = (i32::MAX / 64) as usize;
const _: () = assert!(MAX_MATRIX_SCORE as i64 * (MAX_SENT_LEN as i64 + 3) < i32::MAX as i64);

/// LOSAT_X_TBXUNGAPSENT: whether a batch whose query contexts have the buffers `contexts` may take
/// the sentinel walks: the switch is on and every buffer is `[0, residues.., 0]` with every byte
/// below `MATRIX_LETTERS` (the bytes of the walks' keys, `KeyTable`).
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
pub(crate) fn sent_contexts<'a>(mut contexts: impl Iterator<Item = &'a [u8]>) -> bool {
    if !state().sent {
        return false;
    }
    let ok = contexts.all(|b| b.len() >= 2 && b[0] == 0 && b[b.len() - 1] == 0 && letters(b));
    if !ok {
        crate::utils::xreach::x_reach!("LOSAT_X_TBXUNGAPSENT", false);
    }
    ok
}

/// LOSAT_X_TBXUNGAPSENT: whether the extensions of the chunk `offset..offset + len` of
/// `frame[1..frame.len() - 1]` (`frame`: a subject frame's buffer, `[0, residues.., 0]`) take the
/// sentinel walks: the batch does (`batch`, from [`sent_contexts`]), the chunk is the whole frame
/// and the buffer has its sentinels, at most `MAX_SENT_LEN` bytes and every byte below
/// `MATRIX_LETTERS` (the bytes of the walks' keys, `KeyTable`). Otherwise they take [`two_hit`].
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
#[inline]
pub(crate) fn sent_subject(batch: bool, frame: &[u8], offset: usize, len: usize) -> bool {
    if !batch {
        return false;
    }
    let whole = offset == 0
        && len.checked_add(2) == Some(frame.len())
        && frame.len() <= MAX_SENT_LEN
        && frame[0] == 0
        && frame[len + 1] == 0
        && letters(frame);
    crate::utils::xreach::x_reach!("LOSAT_X_TBXUNGAPSENT", whole);
    whole
}

/// Whether every byte of `buf` is below `MATRIX_LETTERS`: the frames' letters (NCBISTDAA 0..=27)
/// and the sentinel 0 are, so the check refuses no frame of a search; it keeps the walks' keys
/// exact for every buffer that takes them (module comment, LOSAT_X_TBXUNGAPSENT).
// No NCBI counterpart: a condition of a LOSAT code path.
#[inline]
fn letters(buf: &[u8]) -> bool {
    buf.iter().fold(0u8, |a, &c| a | c) < MATRIX_LETTERS as u8
}

/// LOSAT_X_TBXUNGAPSENT: `two_hit(&q_full[1..q_full.len() - 1], &s_full[1..s_full.len() - 1],
/// s_left_off, s_right_off, q_right_off, x_drop)` with walks that the sentinels stop instead of a
/// bound, inlined with the word loop into the callers' extension code. `None` where `two_hit`
/// returns `None`.
///
/// Inlined, the walks reload the table pointer and the X-drop from the stack (one or two loads
/// every two steps), and a call costs more: on AP027131 x AP027133 (code 4) at 1 thread the
/// inlined form ran x0.983 of the same function behind a call (`#[inline(never)]`), with 1.6 %
/// fewer instructions (workstation, 8 alternating rounds).
///
/// SAFETY: `q_full` is `[0, residues.., 0]` (a context of a batch [`sent_contexts`] accepted) and
/// `s_full` the buffer of a subject frame [`sent_subject`] accepted (`[0, residues.., 0]`, at most
/// `MAX_SENT_LEN` bytes).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit` (two_hit.rs follows it).
#[inline(always)]
pub(super) unsafe fn two_hit_sent(
    q_full: &[u8],
    s_full: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> Option<TwoHit> {
    debug_assert!(q_full.len() >= 2 && q_full[0] == 0 && q_full[q_full.len() - 1] == 0);
    debug_assert!(s_full.len() >= 2 && s_full[0] == 0 && s_full[s_full.len() - 1] == 0);
    debug_assert!(s_full.len() <= MAX_SENT_LEN);
    two_hit_sent_on::<WORD>(
        || Some(state().keys.as_deref()?.as_ptr()),
        q_full,
        s_full,
        s_left_off,
        s_right_off,
        q_right_off,
        x_drop,
    )
}

/// LOSAT_X_TBXUNGAPSENT: `two_hit_in_word::<W>(t, &q_full[1..q_full.len() - 1],
/// &s_full[1..s_full.len() - 1], s_left_off, s_right_off, q_right_off, x_drop)` with the walks of
/// `two_hit_sent`.
///
/// SAFETY: those of `two_hit_sent`, and `t.x_sentinels_stop()`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit` (two_hit.rs follows it; `W` is
// its `word_size`).
#[inline(always)]
pub(crate) unsafe fn two_hit_sent_in<const W: usize>(
    t: &MatrixWalkTable,
    q_full: &[u8],
    s_full: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> Option<TwoHit> {
    two_hit_sent_on::<W>(
        || Some(t.1.as_ptr()),
        q_full,
        s_full,
        s_left_off,
        s_right_off,
        q_right_off,
        x_drop,
    )
}

/// `two_hit_sent` on the keys that `keys` gives (read after the calls the path does not cover
/// have returned; `KeyTable`, read at `q << 5 | s`), with words of `W` letters (`WORD` for
/// BLOSUM62's keys) and the key walks (`x_key_walk::walk`).
///
/// SAFETY: those of `two_hit_sent` and `x_key_walk::walk` (module comment, LOSAT_X_TBXUNGAPSENT).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit` (two_hit.rs follows it; `W` is
// its `word_size`).
#[inline(always)]
unsafe fn two_hit_sent_on<const W: usize>(
    keys: impl FnOnce() -> Option<*const i64>,
    q_full: &[u8],
    s_full: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    x_drop: i32,
) -> Option<TwoHit> {
    // The word loop's argument (module comment) holds for `W` letters.
    const { assert!(W >= 1 && W <= MAX_MATRIX_WORD) };
    let q_len = q_full.len() - 2;
    let s_len = s_full.len() - 2;
    if !(1..=MAX_X_DROP).contains(&x_drop) || q_right_off + W > q_len || s_right_off + W > s_len {
        return None;
    }
    let tab = keys()?;
    // offset 0 of the reference's sequences; one sentinel before each
    let qb = q_full.as_ptr().add(1);
    let sb = s_full.as_ptr().add(1);

    // The best end of the word at R (two_hit.rs, the `k_size` loop; inside both sequences).
    let qw = qb.add(q_right_off);
    let sw = sb.add(s_right_off);
    let mut key = 0i64;
    let mut seed = 0i64;
    for i in 0..W {
        key += *tab.add((*qw.add(i) as usize) << 5 | *sw.add(i) as usize);
        if key > seed {
            seed = key;
        }
    }
    let right_d = best_of(seed).1 as usize;
    let q_r = q_right_off + right_d;
    let s_r = s_right_off + right_d;

    // left: from (q_r - 1, s_r - 1) down to the first sentinel (index -1 of the shorter side)
    let (left_score, left_disp, _) =
        walk::<false, u8>(tab, qb.add(q_r).sub(1), sb.add(s_r).sub(1), x_drop, 0);
    let left_disp = left_disp as usize;

    let reached_first_hit = left_disp >= s_r.saturating_sub(s_left_off);
    let (right_score, right_disp, s_last_off) = if reached_first_hit {
        // right: from (q_r, s_r) up to the first sentinel (index q_len / s_len)
        let (best, len, end) = walk::<true, u8>(tab, qb.add(q_r), sb.add(s_r), x_drop, left_score);
        (best, len as usize, s_r + end)
    } else {
        (0, 0, s_r)
    };

    Some((
        q_r - left_disp,
        q_r + right_disp,
        s_r - left_disp,
        s_r + right_disp,
        left_score.max(right_score),
        reached_first_hit,
        s_last_off,
    ))
}

#[cfg(test)]
#[path = "x_two_hit_fast_tests.rs"]
mod tests;
