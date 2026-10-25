//! TBLASTX's two-hit ungapped extension
//! (`two_hit.rs::extend_hit_two_hit`: the word loop and the `s_BlastAaExtendLeft` /
//! `s_BlastAaExtendRight` walks of NCBI `s_BlastAaExtendTwoHit`) with the scores read from one
//! flat byte table.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:846-866
//! ```c
//! for (i = 0; i < n; i++) {
//!     score += matrix[q[i]][s[i]];
//!     if (score > maxscore) { maxscore = score; best_i = i; }
//!     if (score <= 0 || (maxscore - score) >= dropoff) break;
//! }
//! *length = best_i + 1;
//! *s_last_off = s_off + i;
//! ```
//! (right walk; the left walk, aa_ungapped.c:886-921, counts `i` down from `n`, has no
//! `score <= 0` test, and returns `n - best_i + 1`.)
//!
//! The walks of `two_hit.rs` call `get_score` for every step, which tests both codes against the
//! matrix size, branches to the sentinel and the rare codes, and multiplies for the row: 22
//! instructions and five branches a step on the left, about 30 and six on the right (whose
//! bounds are reloaded from the stack; release build of dev 7cdcd4a7). A step here loads the two
//! codes, forms `q << 8 | s`, loads one signed byte from the table and makes the reference's
//! update and exit tests; the statements, their order and the exit tests are those of the
//! reference walks. Each walk is a function of its own, so that its loop has the registers to
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
//!   `q_seq[q_r + j]` while both `q_r + j < q_seq.len()` and `s_r + j < s_seq.len()` on the right.
//!   The walks here take the same counts and read the same positions, so every position is inside
//!   both slices; the word loop reads its three positions after a test that they fit (the
//!   reference reads them unchecked; a word that does not fit runs the reference).
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
//!   `scan_unit_batched_or_bf` and `scan_unit_sliced_or_bf` of x_tbx_hitbf.rs) between two
//!   instances of the loop's function (`const SENT`: `scan_unit_bf`, `scan_unit_batched_bf`,
//!   `scan_unit_sliced_bf`), so that the instances without the walks keep the machine code they had
//!   before the switch. The units of the other hit loop (`LOSAT_X_TBXHITBF=0`, the probes of
//!   `LOSAT_TIMING` and `LOSAT_X_STATS`) take `two_hit`. The reference frame loops (run_impl.rs)
//!   choose at every call.
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
//!   reads the positions the reference reads. The word loop reads the same three positions as in
//!   `two_hit`, after the same test that they fit.
//! * Overflow. A walk without a bound would step over a sentinel if its score wrapped. The left walk
//!   (from 0) and then the right walk (from the left walk's best) read each subject position once at
//!   most, at +11 at most (`MAX_WORD_SCORE`), so a score stays below `11 * (subject letters + 2)`, and
//!   it is never below `best - 255` (a step costs 128 at most, and the walk stops once `best - score`
//!   reaches `x_drop <= 127`); the path takes frames of at most `MAX_SENT_LEN` (`i32::MAX / 16`)
//!   bytes.
//!
//! No NCBI counterpart for the split itself and the switch: they choose between LOSAT code paths
//! with the same result.

use super::get_score;
#[cfg(test)]
use super::two_hit::extend_hit_two_hit_reference;
use crate::utils::matrix::BLAST_SCORE_MIN;
use std::sync::OnceLock;

/// The seven values of `extend_hit_two_hit`: (q_start, q_end, s_start, s_end, score,
/// right_extended, s_last_off).
pub(super) type TwoHit = (usize, usize, usize, usize, i32, bool, usize);

/// `q << 8 | s` for every two bytes.
const TABLE_LEN: usize = 1 << 16;
/// Largest X-drop the byte table is exact for (module comment).
const MAX_X_DROP: i32 = 127;
/// Residues of the word loop (`k_size` of `two_hit.rs`).
const WORD: usize = 3;
/// The largest matrix score the word-loop argument allows (module comment).
const MAX_WORD_SCORE: i32 = 11;

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
}

/// `tab[i]`. wasm32 reads `FIXED_TABLE` (where `tab` points) at its static address.
///
/// SAFETY: `i < TABLE_LEN`.
// NCBI reference: aa_ungapped.c:846-921 `matrix[q[i]][s[i]]` (the table of `build_table`).
#[inline(always)]
unsafe fn score_at(tab: &Table, i: usize) -> i32 {
    #[cfg(target_arch = "wasm32")]
    let tab: &Table = &*FIXED_TABLE.0.get();
    *tab.get_unchecked(i) as i32
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
// NCBI reference: aa_ungapped.c:846-921 `matrix[q[i]][s[i]]`; the scores are those of the function
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
        State { sent, table }
    })
}

/// `s_BlastAaExtendLeft` ending at `(qp, sp)` over `n` steps: step `i` (1..=n) reads `qp[-i]`,
/// `sp[-i]`. Returns the best score and its length.
///
/// SAFETY: the `n` positions before `qp` and before `sp` are readable.
// NCBI reference: aa_ungapped.c:886-921 (module comment). Step `i` reads index `n - i` of the
// `n` positions before the pointers, so one counter runs the loop and the length is
// `n - (the index of the best step)`.
#[inline(never)]
unsafe fn walk_left(tab: &Table, qp: *const u8, sp: *const u8, n: usize, x_drop: i32) -> LeftWalk {
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
        // SAFETY: q << 8 | s < TABLE_LEN for any two bytes.
        score += score_at(tab, (qc as usize) << 8 | sc as usize);
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
// NCBI reference: aa_ungapped.c:846-866 (module comment). Inlined on wasm32 (module comment).
#[cfg_attr(not(target_arch = "wasm32"), inline(never))]
#[cfg_attr(target_arch = "wasm32", inline(always))]
unsafe fn walk_right(
    tab: &Table,
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
        // SAFETY: as in `walk_left`.
        score += score_at(tab, (qc as usize) << 8 | sc as usize);
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

/// `extend_hit_two_hit(q_seq, s_seq, s_left_off, s_right_off, q_right_off, x_drop, false)`, or
/// `None` when the reference has to run (an X-drop outside 1..=127, a word that does not fit).
// NCBI reference: aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit` (two_hit.rs follows it).
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
    let tab: &Table = state().table.as_deref()?;
    // SAFETY: `state().table` is set, so FIXED_TABLE holds the table (FixedTable). wasm32: the word
    // loop below reads it at its static address, as `score_at` does.
    #[cfg(target_arch = "wasm32")]
    let tab: &Table = unsafe { &*FIXED_TABLE.0.get() };

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

/// Longest subject frame of the sentinel walks (`[0, frame.., 0]` bytes): no score wraps (module
/// comment, section LOSAT_X_TBXUNGAPSENT).
const MAX_SENT_LEN: usize = (i32::MAX / 16) as usize;
const _: () = assert!(11 * (MAX_SENT_LEN as i64 + 3) < i32::MAX as i64);

/// LOSAT_X_TBXUNGAPSENT: whether a batch whose query contexts have the buffers `contexts` may take
/// the sentinel walks: the switch is on and every buffer is `[0, residues.., 0]`.
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
pub(crate) fn sent_contexts<'a>(mut contexts: impl Iterator<Item = &'a [u8]>) -> bool {
    if !state().sent {
        return false;
    }
    let ok = contexts.all(|b| b.len() >= 2 && b[0] == 0 && b[b.len() - 1] == 0);
    if !ok {
        crate::utils::xreach::x_reach!("LOSAT_X_TBXUNGAPSENT", false);
    }
    ok
}

/// LOSAT_X_TBXUNGAPSENT: whether the extensions of the chunk `offset..offset + len` of
/// `frame[1..frame.len() - 1]` (`frame`: a subject frame's buffer, `[0, residues.., 0]`) take the
/// sentinel walks: the batch does (`batch`, from [`sent_contexts`]), the chunk is the whole frame
/// and the buffer has its sentinels and at most `MAX_SENT_LEN` bytes. Otherwise they take
/// [`two_hit`].
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
        && frame[len + 1] == 0;
    crate::utils::xreach::x_reach!("LOSAT_X_TBXUNGAPSENT", whole);
    whole
}

/// LOSAT_X_TBXUNGAPSENT: one walk of `s_BlastAaExtendRight` (`RIGHT`, from `start`) or
/// `s_BlastAaExtendLeft` (from 0) from `qp`, `sp` (the first step's positions; the left walk
/// goes down), two steps a loop, until the X-drop test (or, right, `score <= 0`) holds: (best,
/// length, steps before the stop).
///
/// SAFETY: walking from `qp`, `sp` in the direction, both sequences reach a sentinel (code 0) with
/// every position up to it readable; `tab` is the table of `build_table` (code 0 against anything
/// -128, every other score at most `MAX_WORD_SCORE`); `1 <= x_drop <= 127`; and the score cannot
/// wrap before the nearer sentinel (`start` plus `MAX_WORD_SCORE` for every step up to it stays
/// below `i32::MAX`: `two_hit_sent` takes frames of at most `MAX_SENT_LEN` bytes).
// NCBI reference: aa_ungapped.c:846-866 and 886-921 (module comment).
#[inline(always)]
unsafe fn walk_sent<const RIGHT: bool>(
    tab: *const i8,
    qp: *const u8,
    sp: *const u8,
    x_drop: i32,
    start: i32,
) -> (i32, usize, usize) {
    let mut score = start;
    let mut best = start;
    // the offset of the best step; the start value means "no best yet" (length 0)
    let none: isize = if RIGHT { -1 } else { 1 };
    let mut best_off = none;
    let mut off: isize = 0;
    let dir: isize = if RIGHT { 1 } else { -1 };
    macro_rules! step {
        ($l:lifetime) => {{
            let qc = *qp.offset(off) as usize;
            let sc = *sp.offset(off) as usize;
            // SAFETY: q << 8 | s < TABLE_LEN for any two bytes.
            score += *tab.add(qc << 8 | sc) as i32;
            if score > best {
                best = score;
                best_off = off;
            }
            let stop = if RIGHT {
                score <= 0 || best - score >= x_drop
            } else {
                best - score >= x_drop
            };
            if stop {
                break $l;
            }
            off += dir;
        }};
    }
    'w: loop {
        step!('w);
        step!('w);
    }
    let len = if best_off == none {
        0
    } else {
        best_off.unsigned_abs() + 1
    };
    (best, len, off.unsigned_abs())
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
// NCBI reference: aa_ungapped.c:1089-1158 `s_BlastAaExtendTwoHit` (two_hit.rs follows it).
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
    let q_len = q_full.len() - 2;
    let s_len = s_full.len() - 2;
    if !(1..=MAX_X_DROP).contains(&x_drop)
        || q_right_off + WORD > q_len
        || s_right_off + WORD > s_len
    {
        return None;
    }
    let tab = state().table.as_deref()?.as_ptr();
    // offset 0 of the reference's sequences; one sentinel before each
    let qb = q_full.as_ptr().add(1);
    let sb = s_full.as_ptr().add(1);

    // The best end of the word at R (two_hit.rs, the `k_size` loop; inside both sequences).
    let qw = qb.add(q_right_off);
    let sw = sb.add(s_right_off);
    let mut score = 0i32;
    let mut seed_score = 0i32;
    let mut right_d = 0usize;
    for i in 0..WORD {
        score += *tab.add((*qw.add(i) as usize) << 8 | *sw.add(i) as usize) as i32;
        if score > seed_score {
            seed_score = score;
            right_d = i + 1;
        }
    }
    let q_r = q_right_off + right_d;
    let s_r = s_right_off + right_d;

    // left: from (q_r - 1, s_r - 1) down to the first sentinel (index -1 of the shorter side)
    let (left_score, left_disp, _) =
        walk_sent::<false>(tab, qb.add(q_r).sub(1), sb.add(s_r).sub(1), x_drop, 0);

    let reached_first_hit = left_disp >= s_r.saturating_sub(s_left_off);
    let (right_score, right_disp, s_last_off) = if reached_first_hit {
        // right: from (q_r, s_r) up to the first sentinel (index q_len / s_len)
        let (best, len, end) = walk_sent::<true>(tab, qb.add(q_r), sb.add(s_r), x_drop, left_score);
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::tblastx::scoring::Blosum62Score;

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    /// `len` residues from `lo..hi`.
    fn random_seq(rng: &mut Rng, len: usize, lo: u8, hi: u8) -> Vec<u8> {
        (0..len)
            .map(|_| lo + rng.below((hi - lo) as u64) as u8)
            .collect()
    }

    /// A copy of `base` with substitutions of probability `1/rate`.
    fn mutate(rng: &mut Rng, base: &[u8], rate: u64) -> Vec<u8> {
        base.iter()
            .map(|&c| {
                if rng.below(rate) != 0 {
                    c
                } else {
                    1 + rng.below(24) as u8
                }
            })
            .collect()
    }

    /// One call both ways; the fast path must handle it.
    fn check_call(q: &[u8], s: &[u8], sl: usize, sr: usize, qr: usize, d: i32) {
        let expect = extend_hit_two_hit_reference(Blosum62Score, q, s, sl, sr, qr, d, false);
        let got = two_hit(q, s, sl, sr, qr, d).expect("the fast path handles the call");
        assert_eq!(got, expect, "sl={sl} sr={sr} qr={qr} d={d}");
    }

    fn check_pair(q: &[u8], s: &[u8], rng: &mut Rng, calls: usize, dmax: i32) -> u64 {
        if q.len() < WORD || s.len() < WORD {
            return 0;
        }
        for _ in 0..calls {
            let d = 1 + rng.below(dmax as u64) as i32;
            let qr = rng.below((q.len() - WORD + 1) as u64) as usize;
            let sr = rng.below((s.len() - WORD + 1) as u64) as usize;
            // mostly at or before R; sometimes after it, or wrapped (`last_hit + wordsize` as usize)
            let sl = match rng.below(16) {
                0 => sr + 1 + rng.below(3) as usize,
                1 => usize::MAX - rng.below(3) as usize,
                _ => sr.saturating_sub(rng.below(60) as usize),
            };
            check_call(q, s, sl, sr, qr, d);
        }
        calls as u64
    }

    // No NCBI counterpart (two LOSAT paths). Random and related sequences (long walks), every
    // code 1..=27, short and long, X-drops up to the largest the table is exact for.
    #[test]
    fn fast_two_hit_equals_the_reference_on_random_inputs() {
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        let mut compared = 0u64;
        for round in 0..400 {
            let len = match round % 4 {
                0 => WORD + rng.below(40) as usize,
                1 => 20 + rng.below(120) as usize,
                2 => 100 + rng.below(600) as usize,
                _ => 600 + rng.below(3000) as usize,
            };
            let q = random_seq(&mut rng, len, 1, 28);
            let rate = 2 + rng.below(12);
            let s_related = mutate(&mut rng, &q, rate);
            let extra = rng.below(30) as usize;
            let s_random = random_seq(&mut rng, len + extra, 1, 28);
            let dmax = [4, 16, 40, 127][(round / 4 % 4) as usize];
            compared += check_pair(&q, &s_related, &mut rng, 60, dmax);
            compared += check_pair(&q, &s_random, &mut rng, 60, dmax);
            compared += check_pair(
                &s_related[..(len / 2 + WORD).min(len)],
                &q,
                &mut rng,
                20,
                dmax,
            );
        }
        assert!(compared > 50_000);
    }

    // No NCBI counterpart. Two-letter alphabets give ties for the best and drops equal to the
    // X-drop; a gap or sentinel (code 0) inside a sequence stops a walk as BLAST_SCORE_MIN does;
    // codes 28..=255 score as the reference scores them.
    #[test]
    fn fast_two_hit_equals_the_reference_at_ties_sentinels_and_odd_codes() {
        let mut rng = Rng(0xD1B5_4A32_D192_ED03);
        for (a, b) in [
            (1u8, 2u8),
            (1, 20),
            (4, 5),
            (11, 12),
            (19, 25),
            (5, 40),
            (27, 255),
        ] {
            for _ in 0..300 {
                let len = WORD + rng.below(200) as usize;
                let pick = |rng: &mut Rng| if rng.below(2) == 0 { a } else { b };
                let q: Vec<u8> = (0..len).map(|_| pick(&mut rng)).collect();
                let s: Vec<u8> = (0..len + rng.below(10) as usize)
                    .map(|_| pick(&mut rng))
                    .collect();
                let d = 1 + rng.below(12) as i32;
                check_pair(&q, &s, &mut rng, 20, d);
                let mut inner = s.clone();
                let at = rng.below(inner.len() as u64) as usize;
                inner[at] = 0;
                check_pair(&q, &inner, &mut rng, 20, 127);
            }
        }
    }

    // No NCBI counterpart. The calls the fast path leaves to the reference: X-drops outside
    // 1..=127 and a word that does not fit.
    #[test]
    fn fast_two_hit_refuses_what_the_table_cannot_do() {
        let q = vec![5u8; 50];
        let s = vec![5u8; 50];
        assert!(two_hit(&q, &s, 10, 20, 20, 0).is_none());
        assert!(two_hit(&q, &s, 10, 20, 20, 128).is_none());
        check_call(&q, &s, 10, 20, 20, 127);
        assert!(two_hit(&q, &s, 10, 20, 48, 16).is_none());
        assert!(two_hit(&q, &s, 10, 48, 20, 16).is_none());
        check_call(&q, &s, 10, 47, 47, 16);
    }

    // No NCBI counterpart. Every cell of the table is `get_score` of its two bytes, with
    // BLAST_SCORE_MIN as -128.
    #[test]
    fn the_table_holds_get_score_for_every_two_bytes() {
        let tab = build_table().expect("BLOSUM62 fits the byte table");
        for q in 0..=255u8 {
            for s in 0..=255u8 {
                let v = get_score(q, s);
                let want = if v == BLAST_SCORE_MIN { -128 } else { v };
                assert_eq!(
                    tab[(q as usize) << 8 | s as usize] as i32,
                    want,
                    "q={q} s={s}"
                );
            }
        }
    }

    /// `[W.., 0, seq.., 0, W..]` with `PAD` letters W around, and the range of `[0, seq.., 0]`:
    /// a step of the sentinel walks past a sentinel would read a W (+11) and change the result.
    const PAD: usize = 8;
    fn with_sentinels(seq: &[u8]) -> (Vec<u8>, std::ops::Range<usize>) {
        let mut v = vec![20u8; PAD];
        v.push(0);
        v.extend_from_slice(seq);
        v.push(0);
        let r = PAD..v.len();
        v.extend(std::iter::repeat_n(20u8, PAD));
        (v, r)
    }

    /// One call of the sentinel walks against the reference; they must handle it.
    fn check_sent(q: &[u8], s: &[u8], sl: usize, sr: usize, qr: usize, d: i32) {
        let (qb, qr_range) = with_sentinels(q);
        let (sb, sr_range) = with_sentinels(s);
        let expect = extend_hit_two_hit_reference(Blosum62Score, q, s, sl, sr, qr, d, false);
        // SAFETY: both buffers are [0, residues.., 0] and short.
        let got = unsafe { two_hit_sent(&qb[qr_range], &sb[sr_range], sl, sr, qr, d) }
            .expect("the sentinel walks handle the call");
        assert_eq!(got, expect, "sl={sl} sr={sr} qr={qr} d={d}");
    }

    fn check_sent_pair(q: &[u8], s: &[u8], rng: &mut Rng, calls: usize, dmax: i32) -> u64 {
        if q.len() < WORD || s.len() < WORD {
            return 0;
        }
        for _ in 0..calls {
            let d = 1 + rng.below(dmax as u64) as i32;
            // the word anywhere, and often at either end of either sequence
            let qr = match rng.below(8) {
                0 => 0,
                1 => q.len() - WORD,
                _ => rng.below((q.len() - WORD + 1) as u64) as usize,
            };
            let sr = match rng.below(8) {
                0 => 0,
                1 => s.len() - WORD,
                _ => rng.below((s.len() - WORD + 1) as u64) as usize,
            };
            let sl = match rng.below(16) {
                0 => sr + 1 + rng.below(3) as usize,
                1 => usize::MAX - rng.below(3) as usize,
                _ => sr.saturating_sub(rng.below(60) as usize),
            };
            check_sent(q, s, sl, sr, qr, d);
        }
        calls as u64
    }

    // No NCBI counterpart (two LOSAT paths): LOSAT_X_TBXUNGAPSENT against the reference on random
    // and related sequences (walks that run to the sentinels at both ends), short and long, every
    // code 1..=27, X-drops up to 127.
    #[test]
    fn sentinel_walks_equal_the_reference_on_random_inputs() {
        let mut rng = Rng(0x2545_F491_4F6C_DD1D);
        let mut compared = 0u64;
        for round in 0..400 {
            let len = match round % 4 {
                0 => WORD + rng.below(40) as usize,
                1 => 20 + rng.below(120) as usize,
                2 => 100 + rng.below(600) as usize,
                _ => 600 + rng.below(3000) as usize,
            };
            let q = random_seq(&mut rng, len, 1, 28);
            let rate = 2 + rng.below(40);
            let s_related = mutate(&mut rng, &q, rate);
            let extra = rng.below(30) as usize;
            let s_random = random_seq(&mut rng, len + extra, 1, 28);
            let dmax = [4, 16, 40, 127][(round / 4 % 4) as usize];
            compared += check_sent_pair(&q, &s_related, &mut rng, 60, dmax);
            compared += check_sent_pair(&q, &s_random, &mut rng, 60, dmax);
            // the subject shorter than the query, and the other way round
            compared += check_sent_pair(
                &s_related[..(len / 2 + WORD).min(len)],
                &q,
                &mut rng,
                20,
                dmax,
            );
            compared += check_sent_pair(
                &q,
                &s_related[(len / 3).min(len - WORD)..],
                &mut rng,
                20,
                dmax,
            );
        }
        assert!(compared > 60_000);
    }

    // No NCBI counterpart. Identical sequences of W (+11 a step): both walks run out at the
    // sentinels (the reference's bound), from every word position, with the subject offset from
    // the query on both sides.
    #[test]
    fn sentinel_walks_stop_at_both_ends_as_the_bound_does() {
        for (q_len, s_len) in [(3, 3), (3, 9), (9, 3), (17, 17), (40, 25), (25, 40)] {
            let q = vec![20u8; q_len];
            let s = vec![20u8; s_len];
            for qr in 0..=q_len - WORD {
                for sr in 0..=s_len - WORD {
                    for sl in [sr.saturating_sub(5), sr, sr + 1, usize::MAX] {
                        for d in [1, 11, 127] {
                            check_sent(&q, &s, sl, sr, qr, d);
                        }
                    }
                }
            }
        }
    }

    // No NCBI counterpart. Two-letter alphabets (ties for the best, drops equal to the X-drop)
    // and a code 0 inside a sequence (it stops a walk as BLAST_SCORE_MIN does), codes 28..=255.
    #[test]
    fn sentinel_walks_equal_the_reference_at_ties_inner_sentinels_and_odd_codes() {
        let mut rng = Rng(0x9E6C_63D0_676A_9A99);
        for (a, b) in [
            (1u8, 2u8),
            (1, 20),
            (4, 5),
            (11, 12),
            (19, 25),
            (5, 40),
            (27, 255),
        ] {
            for _ in 0..300 {
                let len = WORD + rng.below(200) as usize;
                let pick = |rng: &mut Rng| if rng.below(2) == 0 { a } else { b };
                let q: Vec<u8> = (0..len).map(|_| pick(&mut rng)).collect();
                let s: Vec<u8> = (0..len + rng.below(10) as usize)
                    .map(|_| pick(&mut rng))
                    .collect();
                let d = 1 + rng.below(12) as i32;
                check_sent_pair(&q, &s, &mut rng, 20, d);
                let mut inner = s.clone();
                let at = rng.below(inner.len() as u64) as usize;
                inner[at] = 0;
                check_sent_pair(&q, &inner, &mut rng, 20, 127);
                check_sent_pair(&inner, &q, &mut rng, 20, 127);
            }
        }
    }

    // No NCBI counterpart. The calls the sentinel walks leave to the reference (as `two_hit`):
    // X-drops outside 1..=127 and a word that does not fit.
    #[test]
    fn sentinel_walks_refuse_what_the_table_cannot_do() {
        let (qb, qr) = with_sentinels(&[5u8; 50]);
        let (sb, sr) = with_sentinels(&[5u8; 50]);
        let (q, s) = (&qb[qr.clone()], &sb[sr.clone()]);
        // SAFETY: both buffers are [0, residues.., 0].
        unsafe {
            assert!(two_hit_sent(q, s, 10, 20, 20, 0).is_none());
            assert!(two_hit_sent(q, s, 10, 20, 20, 128).is_none());
            assert!(two_hit_sent(q, s, 10, 20, 48, 16).is_none());
            assert!(two_hit_sent(q, s, 10, 48, 20, 16).is_none());
        }
        check_sent(&[5u8; 50], &[5u8; 50], 10, 20, 20, 127);
        check_sent(&[5u8; 50], &[5u8; 50], 10, 47, 47, 16);
        // `extend_hit_two_hit_sent` runs the reference for the X-drops (the callers' words fit).
        let (qv, sv) = (&qb[qr.start + 1..qr.end - 1], &sb[sr.start + 1..sr.end - 1]);
        for (sl, srr, qrr, d) in [(10, 20, 20, 0), (10, 20, 20, 128), (10, 47, 47, 0)] {
            let expect =
                extend_hit_two_hit_reference(Blosum62Score, qv, sv, sl, srr, qrr, d, false);
            // SAFETY: both buffers are [0, residues.., 0].
            let got =
                unsafe { super::super::two_hit::extend_hit_two_hit_sent(q, s, sl, srr, qrr, d) };
            assert_eq!(got, expect, "sl={sl} sr={srr} qr={qrr} d={d}");
        }
    }

    // No NCBI counterpart. A chunk takes the sentinel walks only when the batch does, the chunk
    // is the whole frame and the buffer has its two sentinels; otherwise the walks of `two_hit`.
    #[test]
    fn only_a_whole_frame_with_its_sentinels_takes_the_sentinel_walks() {
        let frame = [0u8, 4, 5, 6, 7, 0];
        assert!(sent_subject(true, &frame, 0, 4));
        assert!(!sent_subject(false, &frame, 0, 4));
        // a chunk of the frame (MAX_DBSEQ_LEN splits frames longer than 5 M letters)
        assert!(!sent_subject(true, &frame, 1, 3));
        assert!(!sent_subject(true, &frame, 0, 3));
        assert!(!sent_subject(true, &frame, 0, usize::MAX));
        assert!(!sent_subject(true, &[0u8, 4, 5, 6, 7, 8], 0, 4));
        assert!(!sent_subject(true, &[1u8, 4, 5, 6, 7, 0], 0, 4));
        // contexts: every buffer needs both sentinels (and the switch)
        let good: [&[u8]; 2] = [&[0, 4, 5, 6, 0], &[0, 0]];
        assert_eq!(sent_contexts(good.into_iter()), state().sent);
        let bad: [&[u8]; 2] = [&[0, 4, 5, 6, 0], &[0, 4, 5]];
        assert!(!sent_contexts(bad.into_iter()));
        assert!(!sent_contexts([&[0u8][..]].into_iter()));
    }

    // No NCBI counterpart. A gap or sentinel code (0) at every position of the query and of the
    // subject in turn, around a word of W (the largest score, 11, which bounds the word-loop
    // argument): in the word, at the first left step (R - 1), at the first right step (R), and
    // further out; with the first hit at, before and after R, and X-drops 1, 11, 127.
    #[test]
    fn a_sentinel_anywhere_stops_the_walks_as_the_reference_does() {
        const W: u8 = 20;
        assert_eq!(get_score(W, W), MAX_WORD_SCORE);
        let len = 48;
        for at in 0..len {
            for in_query in [true, false] {
                let mut q = vec![W; len];
                let mut s = vec![W; len];
                if in_query {
                    q[at] = 0;
                } else {
                    s[at] = 0;
                }
                for (qr, sr) in [
                    (20, 20),
                    (20, 24),
                    (0, 0),
                    (len - WORD, len - WORD),
                    (3, 40),
                ] {
                    for sl in [sr.saturating_sub(10), sr, sr + 2, usize::MAX] {
                        for d in [1, 11, 127] {
                            check_call(&q, &s, sl, sr, qr, d);
                        }
                    }
                }
            }
        }
    }
}
