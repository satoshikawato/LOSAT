//! EXPERIMENT (LOSAT_X_BPUNGAP / LOSAT_X_BPUNGAPSHADOW, off unless set): the BLOSUM62 protein
//! ungapped extension (`s_BlastAaExtendOneHit`, `s_BlastAaExtendTwoHit` and the
//! `s_BlastAaExtendLeft` / `s_BlastAaExtendRight` walks they share) with a tighter inner loop.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:804-883
//! ```c
//! for (i = 0; i < n; i++) {
//!     score += matrix[q[i]][s[i]];
//!     if (score > maxscore) { maxscore = score; best_i = i; }
//!     if (score <= 0 || (maxscore - score) >= dropoff) break;
//! }
//! *length = best_i + 1;
//! *s_last_off = s_off + i;
//! ```
//! (right walk; the left walk counts `i` down from `n`, has no `score <= 0` test, and returns
//! `n - best_i + 1`.)
//!
//! The reference walks cost about 22 instructions a step: two bounds checks, a range test and a
//! branch for each of the two residue codes, a multiply for the matrix row, and the update. A
//! step here is a load of each residue, `q << 5 | s`, one signed-byte load from a flat table, the
//! add, and the two compare-and-select updates of NCBI's loop; the statements, their order and
//! the exit tests are those of NCBI. (A 16-lane AVX2 version of the walk, with the scores
//! gathered from the same table and the stop found by prefix sums and prefix maxima, was
//! measured too and was not faster: the two `vpgatherdd` of a block and the dependent scans cost
//! as much as 16 scalar steps, whose exit branch is predictable.)
//!
//! Why the result is the same.
//! * Bounds. The reference indexes with `[]`. Here each walk reads `n` positions, with `n`
//!   derived from the slice lengths exactly as the reference derives it (right: the room left in
//!   both slices; left: `min(q_off, s_off) + 1`, after an `assert!` that the start is inside),
//!   so every position is inside both slices.
//! * Table. `T[q << 5 | s]` is `blosum62_score_ncbistdaa_direct(q, s)` (the function the
//!   reference calls) for `q, s < 32`, with the one change that `BLAST_SCORE_MIN` (-32768, the
//!   score of the sentinel code 0 against anything) is stored as -128 so that an entry is one
//!   byte; the table is used only if every other entry is within -127..=127. A residue code
//!   >= 32 would index a part of the 8 KiB table that holds zeros, so each extension ORs
//!   together every code it reads and is handed to the reference when the OR is >= 32.
//! * Sentinel. Take `1 <= dropoff <= 127` (otherwise the reference runs). The score before a
//!   step is never above the best (`best` is raised whenever `score` exceeds it, and both start
//!   equal), so a sentinel step makes `best - score >= 128 > dropoff` with the stored -128 as it
//!   does with -32768: the walk stops at that step in both, the step does not raise the best in
//!   either, and the rest of the result (`best`, `best_i`, the position of the stop) is the
//!   same. The word loops add scores without a stop: for a word of at most 11 residues a
//!   sentinel (-128) is not made positive again by the others (<= 11 each), as -32768 is not,
//!   so every comparison with a score >= 0 has the same outcome; longer words run the
//!   reference.
//!
//! LOSAT_X_BPWALK (off unless set, to any value; it acts only where LOSAT_X_BPUNGAP or its
//! shadow runs this module) replaces the walk `walk` with `drop_walk_right` and
//! `drop_walk_left` / `drop_walk_down`, which keep the drop below the best,
//! `drop = maxscore - score`, in place of `maxscore`, make no code test per step, count with one
//! variable and run in functions of their own (the loop then has the registers to itself):
//! about 40 % fewer instructions a step.
//! Why the result is the same.
//! * Recurrence. Both values start at the `maxscore` passed in, so the drop starts at 0, and
//!   `score` is never above `maxscore` before a step (`drop >= 0`). A step adds `t`. NCBI
//!   raises `maxscore` (and moves `best_i`) exactly when `score + t > maxscore`, that is when
//!   `d = drop - t < 0`, and the drop is then 0; otherwise `maxscore` stays and the drop is
//!   `d`. So per step `best_i` moves iff `d < 0` (strictly, as NCBI's `>`) and the drop becomes
//!   `max(d, 0)`; NCBI's `(maxscore - score) >= dropoff` is `max(d, 0) >= dropoff`, that is
//!   `d >= dropoff` (`dropoff >= 1`); and at the end (the stop step or the last step)
//!   `maxscore = score + drop`. The loops keep `slack = dropoff - drop` (1..=dropoff before a
//!   step): a step adds `t` to it (giving `dropoff - d`), stops the walk when it is <= 0, and
//!   when it is above `dropoff` moves `best_i` and sets it to `dropoff`; so the stop is the sign
//!   of one add, and `maxscore = score + dropoff - slack` at the end.
//! * Order. A step that stops on the drop (`d >= dropoff >= 1`) does not raise the best, so
//!   the test can come before the update; so can the right walk's `score <= 0` when the walk
//!   starts at `maxscore >= 0` (a step with `score <= 0 <= maxscore` does not raise it). The
//!   step's result is NCBI's in either order.
//! * Bounds. With the byte table `t` is in -128..=127 and `1 <= slack <= dropoff <= 127` before
//!   a step, so `slack + t` is in -127..=254: no overflow. The right walk's `score` is above 0
//!   before every step after the first and grows by at most 127 a step, so it stays in i32
//!   when `maxscore + 127 n <= i32::MAX`. The loops run only for `1 <= dropoff <= 127`, the
//!   right one only for `maxscore >= 0` and that sum in range (`drop_walk_proven`; `two_hit`
//!   and `one_hit` always pass it), otherwise `walk` runs. A sentinel step (`t = -128`) gives
//!   `slack + t <= -1`: the walk stops at that step and `best_i` does not move, as in `walk`
//!   and in the reference.
//! * Codes. Instead of `codes |= q | s` at every step, the bytes the walk read (positions 0
//!   to the stop step, or all `n`, a contiguous range of each sequence) are ORed once after the
//!   walk, in 8-byte words. The walk reads only positions inside both slices; a code >= 32
//!   reads a zero from the table (any two bytes index inside it) and makes the OR >= 32, and
//!   the extension is handed to the reference as before. When every code read is below 32 the
//!   scores are exact, so the walk stops where NCBI stops and has read what NCBI read.

use super::extension::{
    BlastpOneHitUngappedResult, BlastpTwoHitUngappedResult, BlastpUngappedData,
};
use crate::utils::matrix::{blosum62_score_ncbistdaa_direct, BLAST_SCORE_MIN};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// Codes below this index the table; `q << 5 | s` is the index.
const CODES: usize = 32;
/// `q` and `s` are bytes, so every `q << 5 | s` is below this: a code >= 32 reads zeros.
const TABLE_LEN: usize = 255 * CODES + 256;
/// Largest X-drop the byte table is exact for (see the module comment).
const MAX_DROPOFF: i32 = 127;
/// Longest word whose sentinel sum stays below zero (11 * 11 < 128).
const MAX_WORD_SIZE: usize = 11;

type Table = [i8; TABLE_LEN];

struct State {
    /// 0 = off, 1 = LOSAT_X_BPUNGAP, 2 = LOSAT_X_BPUNGAPSHADOW (both paths, compared).
    mode: u8,
    /// LOSAT_X_BPWALK: the walks of the fast path are `drop_walk_left` / `drop_walk_right`.
    drop_walk: bool,
    table: Box<Table>,
}

/// The table of the walks, `None` if a score does not fit a byte.
// NCBI reference: aa_ungapped.c:804-883 `matrix[q[i]][s[i]]`; the scores are those of the
// function the reference walks call.
fn build_table() -> Option<Box<Table>> {
    let mut t = Box::new([0i8; TABLE_LEN]);
    for q in 0..CODES {
        for s in 0..CODES {
            let v = blosum62_score_ncbistdaa_direct(q as u8, s as u8);
            t[q << 5 | s] = if v == BLAST_SCORE_MIN {
                -128
            } else if (-127..=127).contains(&v) {
                v as i8
            } else {
                return None;
            };
        }
    }
    Some(t)
}

fn state() -> &'static State {
    static STATE: OnceLock<State> = OnceLock::new();
    STATE.get_or_init(|| {
        // Not in `xswitch::STRICT_SET`, so not read with `xswitch::on`.
        let requested = if std::env::var_os("LOSAT_X_BPUNGAPSHADOW").is_some() {
            2
        } else if std::env::var_os("LOSAT_X_BPUNGAP").is_some() {
            1
        } else {
            0
        };
        // Not in `xswitch::STRICT_SET` either; it changes nothing while `mode` is 0.
        let drop_walk = std::env::var_os("LOSAT_X_BPWALK").is_some();
        // The table is built in every mode (the tests use it); a score that does not fit a byte
        // keeps the switch off.
        match build_table() {
            Some(table) => State {
                mode: requested,
                drop_walk,
                table,
            },
            None => State {
                mode: 0,
                drop_walk: false,
                table: Box::new([0i8; TABLE_LEN]),
            },
        }
    })
}

/// 0 = off, 1 = LOSAT_X_BPUNGAP, 2 = LOSAT_X_BPUNGAPSHADOW.
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
#[inline]
pub(crate) fn mode() -> u8 {
    state().mode
}

static SHADOW_COMPARED: AtomicU64 = AtomicU64::new(0);

/// Counts one compared extension (shadow mode).
// No NCBI counterpart: a counter.
#[inline]
pub(crate) fn count_shadow() {
    SHADOW_COMPARED.fetch_add(1, Ordering::Relaxed);
}

/// The summary line (printed by `main` at exit), only in shadow mode.
// No NCBI counterpart: prints a counter; it does not change any value NCBI computes.
pub fn print_shadow_stats() {
    if mode() == 2 {
        // With LOSAT_X_BPWALK the fast side of every comparison was the drop walk.
        let walk = if state().drop_walk { " walk=drop" } else { "" };
        eprintln!(
            "[X_SHADOW] BPUNGAP compared={}{walk} (ungapped extensions identical)",
            SHADOW_COMPARED.load(Ordering::Relaxed)
        );
    }
}

/// What one walk returns: best score, length, step at which the loop left, and the OR of every
/// residue code it read.
struct Walk {
    best: i32,
    len: i32,
    end: usize,
    codes: u8,
}

/// One walk of `s_BlastAaExtendRight` (`RIGHT`) or `s_BlastAaExtendLeft` over `n` steps: step `i`
/// reads `qp[i]`, `sp[i]` (right) or `qp[-i]`, `sp[-i]` (left).
///
/// SAFETY: the `n` positions of both sequences are readable.
// NCBI reference: aa_ungapped.c:804-883 (see the module comment).
#[inline(always)]
unsafe fn walk<const RIGHT: bool>(
    tab: &Table,
    qp: *const u8,
    sp: *const u8,
    n: usize,
    dropoff: i32,
    maxscore: i32,
) -> Walk {
    let mut score = maxscore;
    let mut best = maxscore;
    let mut best_len = 0i32;
    let mut codes = 0u8;
    let mut i = 0usize;
    while i < n {
        let (qc, sc) = if RIGHT {
            (*qp.add(i), *sp.add(i))
        } else {
            (*qp.sub(i), *sp.sub(i))
        };
        codes |= qc | sc;
        // SAFETY: q << 5 | s < 255 * 32 + 256 = TABLE_LEN for any two bytes.
        score += *tab.get_unchecked((qc as usize) << 5 | sc as usize) as i32;
        if score > best {
            best = score;
            best_len = i as i32 + 1;
        }
        let stop = if RIGHT {
            score <= 0 || (best - score) >= dropoff
        } else {
            (best - score) >= dropoff
        };
        if stop {
            break;
        }
        i += 1;
    }
    Walk {
        best,
        len: best_len,
        end: i,
        codes,
    }
}

/// The OR of the `len` bytes at `qp` and the `len` bytes at `sp` (the codes a walk read), in
/// 8-byte words; a byte may be read more than once (OR does not change).
///
/// SAFETY: the `len` bytes at `qp` and at `sp` are readable.
// No NCBI counterpart: NCBI indexes `matrix[q[i]][s[i]]` with whatever code it reads; the drop
// walks check once per walk, after it, that every code read indexes the matrix part of the
// table (module comment).
#[inline(always)]
unsafe fn or_codes(qp: *const u8, sp: *const u8, len: usize) -> u8 {
    let word = |p: *const u8, at: usize| p.add(at).cast::<u64>().read_unaligned();
    let acc = if len >= 8 {
        // words at 0, 8, 16, .. below len - 8, and the last 8 bytes
        let mut acc = word(qp, len - 8) | word(sp, len - 8);
        let mut at = 0;
        while at + 8 < len {
            acc |= word(qp, at) | word(sp, at);
            at += 8;
        }
        acc
    } else if len >= 4 {
        let half = |p: *const u8, at: usize| u64::from(p.add(at).cast::<u32>().read_unaligned());
        half(qp, 0) | half(qp, len - 4) | half(sp, 0) | half(sp, len - 4)
    } else if len > 0 {
        // 0, len / 2 and len - 1 are every position of 1..=3 bytes
        let byte = |p: *const u8, at: usize| u64::from(*p.add(at));
        byte(qp, 0)
            | byte(qp, len / 2)
            | byte(qp, len - 1)
            | byte(sp, 0)
            | byte(sp, len / 2)
            | byte(sp, len - 1)
    } else {
        0
    };
    let acc = acc | acc >> 32;
    let acc = acc | acc >> 16;
    (acc | acc >> 8) as u8
}

/// Whether the arithmetic of the drop walks is proven for these inputs (module comment,
/// Bounds): `1 <= dropoff <= 127`, and for the right walk `maxscore >= 0` and
/// `maxscore + 127 n <= i32::MAX`. `two_hit` and `one_hit` always pass them (`supported`; the
/// right walk starts from a left walk's best, which is never below its start of 0 or of the
/// word's best >= 0; a walk of 16.9 million steps is longer than any sequence it runs on);
/// otherwise the walks run `walk`.
// No NCBI counterpart: a range check of LOSAT's own walk.
#[inline(always)]
fn drop_walk_proven(right: bool, n: usize, dropoff: i32, maxscore: i32) -> bool {
    (1..=MAX_DROPOFF).contains(&dropoff)
        && (!right || (maxscore >= 0 && n <= ((i32::MAX - maxscore) / 127) as usize))
}

/// The walk of `s_BlastAaExtendRight` over `n` steps from `qp`, `sp` (step `i` reads `qp[i]`,
/// `sp[i]`), keeping the drop below the best (LOSAT_X_BPWALK, module comment); the same
/// `Walk` as `walk::<true>`. Not inlined: the loop gets the registers to itself.
///
/// SAFETY: the `n` positions of both sequences are readable.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:846-864
// ```c
// for (i = 0; i < n; i++) {
//     score += matrix[q[i]][s[i]];
//
//     if (score > maxscore) {
//         maxscore = score;
//         best_i = i;
//     }
//     ...
//     if (score <= 0 || (maxscore - score) >= dropoff)
//         break;
// }
//
// *length = best_i + 1;
// *s_last_off = s_off + i;
// ```
#[inline(never)]
unsafe fn drop_walk_right(
    tab: &Table,
    qp: *const u8,
    sp: *const u8,
    n: usize,
    dropoff: i32,
    maxscore: i32,
) -> Walk {
    if !drop_walk_proven(true, n, dropoff, maxscore) {
        return walk::<true>(tab, qp, sp, n, dropoff, maxscore);
    }
    let mut score = maxscore;
    // dropoff - (maxscore - score), 1..=dropoff before a step
    let mut slack = dropoff;
    let mut best_i = -1isize;
    let mut i = 0usize;
    while i < n {
        let (qc, sc) = (*qp.add(i), *sp.add(i));
        // SAFETY: q << 5 | s < 255 * 32 + 256 = TABLE_LEN for any two bytes.
        let t = *tab.get_unchecked((qc as usize) << 5 | sc as usize) as i32;
        // SAFETY: 0 <= score <= maxscore + 127 i before step i (score > 0 after every step
        // that did not stop the walk), and maxscore + 127 n <= i32::MAX (`drop_walk_proven`).
        score = score.unchecked_add(t);
        // NCBI tests `score <= 0` after the update of the best; with maxscore >= 0 such a step
        // does not raise the best, so the test can come first (the slack of the step is still
        // added: `best` below is `score + (dropoff - slack)`).
        if score <= 0 {
            slack += t;
            break;
        }
        // SAFETY: 1 <= slack <= dropoff <= 127 and -128 <= t <= 127.
        slack = slack.unchecked_add(t);
        // (maxscore - score) >= dropoff; such a step does not raise the best either
        if slack <= 0 {
            break;
        }
        // score > maxscore: the best moves here and the drop is 0
        if slack > dropoff {
            best_i = i as isize;
            slack = dropoff;
        }
        i += 1;
    }
    // positions 0..=i were read (all n when the loop ran out)
    let read = if i < n { i + 1 } else { n };
    Walk {
        best: score + (dropoff - slack),
        len: (best_i + 1) as i32,
        end: i,
        codes: or_codes(qp, sp, read),
    }
}

/// The walk of `s_BlastAaExtendLeft` over `n` steps down from `qp`, `sp` (step `k` reads
/// `qp[-k]`, `sp[-k]`), keeping the drop below the best (LOSAT_X_BPWALK, module comment); the
/// same `Walk` as `walk::<false>`.
///
/// SAFETY: the `n` positions of both sequences are readable.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:897-901
// ```c
// n = MIN(s_off, q_off);
// best_i = n + 1;
//
// s = subject->sequence + s_off - n;
// q = query->sequence + q_off - n;
// ```
#[inline(always)]
unsafe fn drop_walk_left(
    tab: &Table,
    qp: *const u8,
    sp: *const u8,
    n: usize,
    dropoff: i32,
    maxscore: i32,
) -> Walk {
    if n == 0 {
        return Walk {
            best: maxscore,
            len: 0,
            end: 0,
            codes: 0,
        };
    }
    // NCBI's n, and its q and s: the lowest positions of the walk
    drop_walk_down(tab, qp.sub(n - 1), sp.sub(n - 1), n - 1, dropoff, maxscore)
}

/// The loop of [`drop_walk_left`], from position `last` of `q`, `s` down to position 0, in
/// NCBI's indexing. Not inlined: the loop gets the registers to itself.
///
/// SAFETY: positions 0..=last of both sequences are readable.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:903-919
// ```c
// for (i = n; i >= 0; i--) {
//     score += matrix[q[i]][s[i]];
//
//     if (score > maxscore) {
//         maxscore = score;
//         best_i = i;
//     }
//     ...
//     if ((maxscore - score) >= dropoff)
//         break;
// }
//
// *length = n - best_i + 1;
// ```
#[inline(never)]
unsafe fn drop_walk_down(
    tab: &Table,
    q: *const u8,
    s: *const u8,
    last: usize,
    dropoff: i32,
    maxscore: i32,
) -> Walk {
    if !drop_walk_proven(false, last + 1, dropoff, maxscore) {
        return walk::<false>(tab, q.add(last), s.add(last), last + 1, dropoff, maxscore);
    }
    let mut score = maxscore;
    // dropoff - (maxscore - score), 1..=dropoff before a step
    let mut slack = dropoff;
    let mut best_i = last as isize + 1;
    let mut i = last as isize;
    while i >= 0 {
        let (qc, sc) = (*q.offset(i), *s.offset(i));
        // SAFETY: q << 5 | s < 255 * 32 + 256 = TABLE_LEN for any two bytes.
        let t = *tab.get_unchecked((qc as usize) << 5 | sc as usize) as i32;
        score += t;
        // SAFETY: 1 <= slack <= dropoff <= 127 and -128 <= t <= 127.
        slack = slack.unchecked_add(t);
        // (maxscore - score) >= dropoff; such a step does not raise the best, so the test can
        // come before the update
        if slack <= 0 {
            break;
        }
        // score > maxscore: the best moves here and the drop is 0
        if slack > dropoff {
            best_i = i;
            slack = dropoff;
        }
        i -= 1;
    }
    // positions i..=last were read (0..=last when the loop ran out, i = -1)
    let first = i.max(0) as usize;
    Walk {
        best: score + (dropoff - slack),
        len: (last as isize - best_i + 1) as i32,
        end: (last as isize - i) as usize,
        codes: or_codes(q.add(first), s.add(first), last + 1 - first),
    }
}

/// `s_BlastAaExtendLeft(.., s_off, q_off, ..)`: (score, length, codes read); the drop walk with
/// `DROP` (LOSAT_X_BPWALK).
#[inline(always)]
fn extend_left<const DROP: bool>(
    tab: &Table,
    query: &[u8],
    subject: &[u8],
    q_off: i32,
    s_off: i32,
    dropoff: i32,
    maxscore: i32,
) -> (i32, i32, u8) {
    let n = s_off.min(q_off);
    if n < 0 {
        // NCBI reference: aa_ungapped.c:853-883, `best_i = n + 1` with no step taken.
        return (maxscore, 0, 0);
    }
    // n + 1 steps down from (q_off, s_off); both offsets are >= n.
    assert!((q_off as usize) < query.len() && (s_off as usize) < subject.len());
    let (qp, sp) = (query.as_ptr(), subject.as_ptr());
    // SAFETY: the assert and n <= min(q_off, s_off) keep the n + 1 positions inside.
    let w = unsafe {
        let (qp, sp) = (qp.add(q_off as usize), sp.add(s_off as usize));
        if DROP {
            drop_walk_left(tab, qp, sp, n as usize + 1, dropoff, maxscore)
        } else {
            walk::<false>(tab, qp, sp, n as usize + 1, dropoff, maxscore)
        }
    };
    (w.best, w.len, w.codes)
}

/// `s_BlastAaExtendRight(.., s_off, q_off, ..)`: (score, length, s_last_off, codes read); the
/// drop walk with `DROP` (LOSAT_X_BPWALK).
#[inline(always)]
fn extend_right<const DROP: bool>(
    tab: &Table,
    query: &[u8],
    subject: &[u8],
    q_off: i32,
    s_off: i32,
    dropoff: i32,
    maxscore: i32,
) -> (i32, i32, i32, u8) {
    let n = (subject.len() as i32 - s_off).min(query.len() as i32 - q_off);
    if n <= 0 {
        // NCBI reference: aa_ungapped.c:804-835, the loop does not run: i = 0.
        return (maxscore, 0, s_off, 0);
    }
    // n steps up from (q_off, s_off); n is the room left in both slices. A negative offset was
    // an index panic in the reference.
    assert!(q_off >= 0 && s_off >= 0);
    let (qp, sp) = (query.as_ptr(), subject.as_ptr());
    // SAFETY: the n positions are inside both slices.
    let w = unsafe {
        let (qp, sp) = (qp.add(q_off as usize), sp.add(s_off as usize));
        if DROP {
            drop_walk_right(tab, qp, sp, n as usize, dropoff, maxscore)
        } else {
            walk::<true>(tab, qp, sp, n as usize, dropoff, maxscore)
        }
    };
    (w.best, w.len, s_off + w.end as i32, w.codes)
}

/// Whether the fast path covers these parameters.
#[inline(always)]
fn supported(query: &[u8], subject: &[u8], dropoff: i32, word_size: usize) -> bool {
    (1..=MAX_DROPOFF).contains(&dropoff)
        && (1..=MAX_WORD_SIZE).contains(&word_size)
        && query.len() <= i32::MAX as usize
        && subject.len() <= i32::MAX as usize
}

/// `Some(result)` when the fast path handled the call (`result` may be `None`, the reference's
/// early return), `None` when the reference has to run. The walks are those of LOSAT_X_BPWALK
/// when it is set.
// No NCBI counterpart: chooses between two LOSAT walks with the same result; the extension is
// `two_hit_with` (aa_ungapped.c:1089-1155 `s_BlastAaExtendTwoHit`).
#[inline]
pub(crate) fn two_hit(
    query: &[u8],
    subject: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<Option<BlastpTwoHitUngappedResult>> {
    if state().drop_walk {
        two_hit_with::<true>(
            query,
            subject,
            s_left_off,
            s_right_off,
            q_right_off,
            dropoff,
            word_size,
        )
    } else {
        two_hit_with::<false>(
            query,
            subject,
            s_left_off,
            s_right_off,
            q_right_off,
            dropoff,
            word_size,
        )
    }
}

/// [`two_hit`] with the walks `walk` or, with `DROP`, `drop_walk_left` / `drop_walk_right`.
// NCBI reference: aa_ungapped.c:1089-1155 `s_BlastAaExtendTwoHit`.
#[inline(never)]
fn two_hit_with<const DROP: bool>(
    query: &[u8],
    subject: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<Option<BlastpTwoHitUngappedResult>> {
    // The early return of the reference wrapper (a word that does not fit).
    if q_right_off + word_size > query.len() || s_right_off + word_size > subject.len() {
        return Some(None);
    }
    if !supported(query, subject, dropoff, word_size) {
        return None;
    }
    let tab: &Table = &state().table;
    let mut q_right_off = q_right_off as i32;
    let s_left_off = s_left_off as i32;
    let mut s_right_off = s_right_off as i32;
    let word_size = word_size as i32;

    let mut right_d = 0i32;
    let mut left_score = 0i32;
    let mut score = 0i32;
    let mut codes = 0u8;
    for i in 0..word_size {
        // the word fits in both sequences (checked above)
        let qc = query[(q_right_off + i) as usize];
        let sc = subject[(s_right_off + i) as usize];
        codes |= qc | sc;
        score += tab[(qc as usize) << 5 | sc as usize] as i32;
        if score > left_score {
            left_score = score;
            right_d = i + 1;
        }
    }
    q_right_off += right_d;
    s_right_off += right_d;

    let (left_score, left_d, c) = extend_left::<DROP>(
        tab,
        query,
        subject,
        q_right_off - 1,
        s_right_off - 1,
        dropoff,
        0,
    );
    codes |= c;

    let mut right_extend = false;
    let mut s_last_off = s_right_off;
    let mut right_score = 0i32;
    let mut right_d = 0i32;
    if left_d >= (s_right_off - s_left_off) {
        right_extend = true;
        let (rs, rd, last, c) = extend_right::<DROP>(
            tab,
            query,
            subject,
            q_right_off,
            s_right_off,
            dropoff,
            left_score,
        );
        (right_score, right_d, s_last_off) = (rs, rd, last);
        codes |= c;
    }
    if codes >= CODES as u8 {
        return None;
    }

    Some(Some(BlastpTwoHitUngappedResult {
        ungapped_data: BlastpUngappedData {
            q_start: q_right_off - left_d,
            s_start: s_right_off - left_d,
            length: left_d + right_d,
            score: left_score.max(right_score),
        },
        right_extend,
        s_last_off,
    }))
}

/// As [`two_hit`], for `s_BlastAaExtendOneHit`.
// No NCBI counterpart: chooses between two LOSAT walks with the same result; the extension is
// `one_hit_with` (aa_ungapped.c:1020-1083 `s_BlastAaExtendOneHit`).
#[inline]
pub(crate) fn one_hit(
    query: &[u8],
    subject: &[u8],
    q_off: usize,
    s_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<Option<BlastpOneHitUngappedResult>> {
    if state().drop_walk {
        one_hit_with::<true>(query, subject, q_off, s_off, dropoff, word_size)
    } else {
        one_hit_with::<false>(query, subject, q_off, s_off, dropoff, word_size)
    }
}

/// [`one_hit`] with the walks `walk` or, with `DROP`, `drop_walk_left` / `drop_walk_right`.
// NCBI reference: aa_ungapped.c:1020-1083 `s_BlastAaExtendOneHit`.
#[inline(never)]
fn one_hit_with<const DROP: bool>(
    query: &[u8],
    subject: &[u8],
    q_off: usize,
    s_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<Option<BlastpOneHitUngappedResult>> {
    if q_off + word_size > query.len() || s_off + word_size > subject.len() {
        return Some(None);
    }
    if !supported(query, subject, dropoff, word_size) {
        return None;
    }
    let tab: &Table = &state().table;
    let q_off = q_off as i32;
    let s_off = s_off as i32;
    let word_size = word_size as i32;

    let mut score = 0i32;
    let mut sum = 0i32;
    let mut q_left_off = q_off;
    let mut q_right_off = q_off + word_size;
    let mut q_best_left_off = q_off;
    let mut codes = 0u8;
    for i in 0..word_size {
        // the word fits in both sequences (checked above)
        let qc = query[(q_off + i) as usize];
        let sc = subject[(s_off + i) as usize];
        codes |= qc | sc;
        sum += tab[(qc as usize) << 5 | sc as usize] as i32;
        if sum > score {
            score = sum;
            q_best_left_off = q_left_off;
            q_right_off = q_off + i;
        } else if sum <= 0 {
            sum = 0;
            q_left_off = q_off + i + 1;
        }
    }

    let init_hit_width = q_right_off - q_left_off + 1;
    q_left_off = q_best_left_off;
    let s_left_off = q_left_off + (s_off - q_off);
    let s_right_off = q_right_off + (s_off - q_off);

    let (left_score, left_disp, c) = extend_left::<DROP>(
        tab,
        query,
        subject,
        q_left_off - 1,
        s_left_off - 1,
        dropoff,
        score,
    );
    codes |= c;
    let (total_score, right_disp, s_last_off, c) = extend_right::<DROP>(
        tab,
        query,
        subject,
        q_right_off + 1,
        s_right_off + 1,
        dropoff,
        left_score,
    );
    codes |= c;
    if codes >= CODES as u8 {
        return None;
    }

    Some(Some(BlastpOneHitUngappedResult {
        ungapped_data: BlastpUngappedData {
            q_start: q_left_off - left_disp,
            s_start: s_left_off - left_disp,
            length: left_disp + right_disp + init_hit_width,
            score: total_score,
        },
        s_last_off,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastp::extension::{
        extend_left as ref_left, extend_one_hit, extend_right as ref_right, extend_two_hit,
    };
    use crate::config::ScoringMatrix;

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

    /// A sequence of `len` residues from `lo..hi` with the sentinel 0 at both ends.
    fn random_seq(rng: &mut Rng, len: usize, lo: u8, hi: u8) -> Vec<u8> {
        let mut v = vec![0u8];
        for _ in 0..len {
            v.push(lo + rng.below((hi - lo) as u64) as u8);
        }
        v.push(0);
        v
    }

    /// A mutated copy of `base` (substitutions with probability `1/rate`).
    fn mutate(rng: &mut Rng, base: &[u8], rate: u64) -> Vec<u8> {
        base.iter()
            .enumerate()
            .map(|(i, &c)| {
                if i == 0 || i + 1 == base.len() || rng.below(rate) != 0 {
                    c
                } else {
                    1 + rng.below(24) as u8
                }
            })
            .collect()
    }

    /// Both paths on one call (the reference through the matrix-generic entry points, which the
    /// switch does not touch); the fast path must handle it, with either walk (`walk` and the
    /// drop walks of LOSAT_X_BPWALK) and through the entry point that chooses between them.
    fn check_call(q: &[u8], s: &[u8], sl: usize, sr: usize, qr: usize, d: i32, w: usize) {
        let expect = extend_two_hit(ScoringMatrix::Blosum62, q, s, sl, sr, qr, d, w);
        for (got, walk) in [
            (two_hit_with::<false>(q, s, sl, sr, qr, d, w), "walk"),
            (two_hit_with::<true>(q, s, sl, sr, qr, d, w), "drop walk"),
            (two_hit(q, s, sl, sr, qr, d, w), "chosen walk"),
        ] {
            let got = got.expect("fast path handles the call");
            assert_eq!(
                got, expect,
                "two_hit ({walk}) sl={sl} sr={sr} qr={qr} d={d} w={w}"
            );
        }
        let expect = extend_one_hit(ScoringMatrix::Blosum62, q, s, qr, sr, d, w);
        for (got, walk) in [
            (one_hit_with::<false>(q, s, qr, sr, d, w), "walk"),
            (one_hit_with::<true>(q, s, qr, sr, d, w), "drop walk"),
            (one_hit(q, s, qr, sr, d, w), "chosen walk"),
        ] {
            let got = got.expect("fast path handles the call");
            assert_eq!(got, expect, "one_hit ({walk}) qr={qr} sr={sr} d={d} w={w}");
        }
    }

    fn check_pair(q: &[u8], s: &[u8], rng: &mut Rng, calls: usize, dmax: i32) -> u64 {
        for _ in 0..calls {
            let w = 1 + rng.below(4) as usize;
            let d = 1 + rng.below(dmax as u64) as i32;
            let qr = rng.below(q.len() as u64) as usize;
            let sr = rng.below(s.len() as u64) as usize;
            let sl = if sr > 0 {
                rng.below(sr as u64) as usize
            } else {
                0
            };
            check_call(q, s, sl, sr, qr, d, w);
        }
        2 * calls as u64
    }

    /// Random sequences, related sequences (long walks), all codes, short and long.
    #[test]
    fn fast_extensions_equal_the_reference_on_random_inputs() {
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        let mut compared = 0u64;
        for round in 0..400 {
            let len = match round % 4 {
                0 => 1 + rng.below(40) as usize,
                1 => 20 + rng.below(120) as usize,
                2 => 100 + rng.below(600) as usize,
                _ => 600 + rng.below(3000) as usize,
            };
            // residues 1..=27, so every code except the sentinel appears
            let q = random_seq(&mut rng, len, 1, 28);
            let rate = 2 + rng.below(12);
            let s_related = mutate(&mut rng, &q, rate);
            let extra = rng.below(30) as usize;
            let s_random = random_seq(&mut rng, len + extra, 1, 28);
            let dmax = [4, 16, 40, 127][(round / 4 % 4) as usize];
            compared += check_pair(&q, &s_related, &mut rng, 60, dmax);
            compared += check_pair(&q, &s_random, &mut rng, 60, dmax);
            // a sentinel inside the sequence (between contexts)
            let mut inner = q.clone();
            let at = 1 + rng.below(len as u64) as usize;
            inner[at] = 0;
            compared += check_pair(&inner, &s_related, &mut rng, 40, dmax);
        }
        assert!(compared > 100_000);
    }

    /// X-drop at the boundary and ties for the best position: walks whose score returns to the
    /// best value exactly, drops that equal the dropoff, repeated patterns.
    #[test]
    fn fast_extensions_equal_the_reference_at_boundaries_and_ties() {
        let mut rng = Rng(0xD1B5_4A32_D192_ED03);
        let mut compared = 0u64;
        // Two-letter alphabets give plus, zero and minus steps in equal measure, so the score
        // returns to the same value over and over (ties) and meets best - score == dropoff.
        for pairs in [
            (1u8, 2u8),
            (1, 20),
            (20, 23),
            (4, 5),
            (13, 14),
            (1, 4),
            (17, 20),
        ] {
            for round in 0..300 {
                let len = 1 + rng.below(if round % 3 == 0 { 50 } else { 400 }) as usize;
                let mut draw = |n: usize| -> Vec<u8> {
                    std::iter::once(0)
                        .chain((0..n).map(|_| if rng.below(2) == 0 { pairs.0 } else { pairs.1 }))
                        .chain(std::iter::once(0))
                        .collect()
                };
                let q = draw(len);
                let s = draw(len + 3);
                compared += check_pair(&q, &s, &mut rng, 30, 25);
            }
        }
        // Every dropoff 1..=60 on a periodic pattern: +11 (W) among +4 (A-A) steps.
        let w_code = 20u8;
        let q: Vec<u8> = [0u8]
            .into_iter()
            .chain((0..300).map(|i| if i % 7 == 3 { w_code } else { 1 }))
            .chain([0u8])
            .collect();
        let s: Vec<u8> = [0u8]
            .into_iter()
            .chain((0..300).map(|i| if i % 5 == 1 { w_code } else { 1 }))
            .chain([0u8])
            .collect();
        for d in 1..=60 {
            for sr in 1..q.len() - 4 {
                for w in [1usize, 3] {
                    check_call(&q, &s, sr.saturating_sub(5), sr, sr, d, w);
                    compared += 2;
                }
            }
        }
        assert!(compared > 50_000);
    }

    /// Very long walks: identical sequences whose score passes 16 bits, with the sentinels at both
    /// ends, and walks that start next to either end.
    #[test]
    fn fast_extensions_equal_the_reference_on_very_long_walks() {
        let mut rng = Rng(0x2545_F491_4F6C_DD1D);
        // 9000 tryptophans: 99,000 > 2^16 (W-W = 11)
        let mut q = vec![0u8];
        q.extend(std::iter::repeat_n(20u8, 9000));
        q.push(0);
        let mut s = q.clone();
        for d in [1, 5, 16, 100, 127] {
            for (qr, sr) in [
                (1usize, 1usize),
                (4000, 4000),
                (8990, 8990),
                (5000, 20),
                (20, 5000),
            ] {
                check_call(&q, &s, sr.saturating_sub(10), sr, qr, d, 3);
            }
        }
        // a long walk with mismatches that never reach the dropoff
        for i in (0..s.len()).step_by(9) {
            if s[i] != 0 {
                s[i] = 1;
            }
        }
        for _ in 0..200 {
            let qr = 1 + rng.below(8990) as usize;
            let sr = 1 + rng.below(8990) as usize;
            let sl = sr.saturating_sub(rng.below(30) as usize);
            let d = 20 + rng.below(100) as i32;
            check_call(&q, &s, sl, sr, qr, d, 3);
        }
    }

    /// Words of 1..=11 residues with sentinels and the extreme scores inside them, and words next
    /// to either end of a sequence.
    #[test]
    fn fast_extensions_equal_the_reference_for_sentinels_in_words() {
        let mut rng = Rng(0x1234_5678_9ABC_DEF1);
        for _ in 0..20_000 {
            let len = 12 + rng.below(40) as usize;
            let mut q = random_seq(&mut rng, len, 1, 28);
            let mut s = random_seq(&mut rng, len, 1, 28);
            // plenty of W/W (11) and sentinels (0) in the sequences
            for v in [&mut q, &mut s] {
                for _ in 0..rng.below(6) {
                    let at = rng.below(v.len() as u64) as usize;
                    v[at] = if rng.below(2) == 0 { 20 } else { 0 };
                }
            }
            let w = 1 + rng.below(11) as usize;
            let d = 1 + rng.below(127) as i32;
            let qr = rng.below(q.len() as u64) as usize;
            let sr = rng.below(s.len() as u64) as usize;
            let sl = sr.saturating_sub(rng.below(8) as usize);
            check_call(&q, &s, sl, sr, qr, d, w);
        }
    }

    /// Parameters outside the proven range, and a code the table does not cover, go to the
    /// reference.
    #[test]
    fn fast_path_hands_unsupported_inputs_to_the_reference() {
        let q = vec![
            0u8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 0,
        ];
        let s = q.clone();
        for d in [i32::MIN, -5, 0, 128, 8193, i32::MAX] {
            assert_eq!(two_hit(&q, &s, 2, 6, 6, d, 3), None, "d={d}");
            assert_eq!(one_hit(&q, &s, 6, 6, d, 3), None, "d={d}");
            assert_eq!(two_hit_with::<true>(&q, &s, 2, 6, 6, d, 3), None, "d={d}");
            assert_eq!(one_hit_with::<true>(&q, &s, 6, 6, d, 3), None, "d={d}");
        }
        for w in [0usize, 12, 20] {
            assert_eq!(two_hit(&q, &s, 0, 2, 2, 16, w), None, "w={w}");
            assert_eq!(one_hit(&q, &s, 2, 2, 16, w), None, "w={w}");
            assert_eq!(two_hit_with::<true>(&q, &s, 0, 2, 2, 16, w), None, "w={w}");
            assert_eq!(one_hit_with::<true>(&q, &s, 2, 2, 16, w), None, "w={w}");
        }
        let mut wide = vec![0u8; 80];
        let mut other = vec![0u8; 80];
        for i in 1..79 {
            wide[i] = 1 + (i % 20) as u8;
            other[i] = 1 + ((i * 3) % 20) as u8;
        }
        wide[36] = 200;
        assert_eq!(two_hit(&wide, &other, 30, 35, 35, 16, 3), None);
        assert_eq!(one_hit(&wide, &other, 35, 35, 16, 3), None);
        assert_eq!(two_hit_with::<true>(&wide, &other, 30, 35, 35, 16, 3), None);
        assert_eq!(one_hit_with::<true>(&wide, &other, 35, 35, 16, 3), None);
        // a code >= 32 only in the walk, not in the word
        let mut deep = other.clone();
        deep[20] = 77;
        assert_eq!(two_hit(&other, &deep, 28, 30, 30, 100, 3), None);
        assert_eq!(
            two_hit_with::<true>(&other, &deep, 28, 30, 30, 100, 3),
            None
        );
        // the early return for a word that does not fit is the reference's
        assert_eq!(two_hit(&q, &s, 0, 21, 3, 16, 3), Some(None));
        assert_eq!(one_hit(&q, &s, 3, 21, 16, 3), Some(None));
    }

    /// The table the fast path scores with is the matrix function the reference calls (the
    /// sentinel is stored as -128, a stop in the same step).
    #[test]
    fn score_table_equals_the_reference_matrix_function() {
        let t = build_table().expect("BLOSUM62 fits bytes");
        for q in 0..32usize {
            for s in 0..32usize {
                let v = blosum62_score_ncbistdaa_direct(q as u8, s as u8);
                let expect = if v == BLAST_SCORE_MIN { -128 } else { v };
                assert_eq!(t[q << 5 | s] as i32, expect, "q={q} s={s}");
            }
        }
        assert!(t[CODES * CODES..].iter().all(|&b| b == 0));
    }

    // LOSAT_X_BPWALK: the drop walks against `walk`, the reference functions of extension.rs and
    // a transcription of NCBI's loops.

    /// Why one of NCBI's walks stopped.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum Exit {
        /// after all `n` steps
        Ran,
        /// right walk: `score <= 0`
        NotPositive,
        /// `maxscore - score == dropoff`
        DropEqual,
        /// `maxscore - score > dropoff`
        DropAbove,
    }

    /// NCBI's walk (aa_ungapped.c:846-864 right, 903-919 left) over `steps` (residue pairs in
    /// step order) in i64 with the matrix function's scores (the sentinel's -32768): (maxscore,
    /// steps to the best, the step it stopped at or `n`, why, whether it passed a drop of
    /// `dropoff - 1` without stopping).
    fn ncbi_walk(
        steps: &[(u8, u8)],
        right: bool,
        dropoff: i32,
        maxscore: i32,
    ) -> (i64, i64, usize, Exit, bool) {
        let dropoff = dropoff as i64;
        let mut score = maxscore as i64;
        let mut best = score;
        let mut best_len = 0i64;
        let mut near = false;
        for (i, &(q, s)) in steps.iter().enumerate() {
            score += blosum62_score_ncbistdaa_direct(q, s) as i64;
            if score > best {
                best = score;
                best_len = i as i64 + 1;
            }
            if right && score <= 0 {
                return (best, best_len, i, Exit::NotPositive, near);
            }
            if best - score >= dropoff {
                let why = if best - score == dropoff {
                    Exit::DropEqual
                } else {
                    Exit::DropAbove
                };
                return (best, best_len, i, why, near);
            }
            near |= best - score == dropoff - 1;
        }
        (best, best_len, steps.len(), Exit::Ran, near)
    }

    fn fields(w: Walk) -> (i32, i32, usize, u8) {
        (w.best, w.len, w.end, w.codes)
    }

    /// The walks over `steps` (codes below 32) both ways: right from position 0 of two sequences
    /// of exactly `steps.len()` residues, left from the last position of the reversed sequences
    /// down to position 0. The drop walk must return `walk`'s `Walk`, the extension functions
    /// with either walk the reference functions' results, and the walk NCBI's values. Returns
    /// NCBI's (stop step, exit, near) of the right and of the left walk.
    fn check_walks(
        tab: &Table,
        steps: &[(u8, u8)],
        dropoff: i32,
        maxscore: i32,
    ) -> [(usize, Exit, bool); 2] {
        let n = steps.len();
        let q: Vec<u8> = steps.iter().map(|p| p.0).collect();
        let s: Vec<u8> = steps.iter().map(|p| p.1).collect();
        let ql: Vec<u8> = q.iter().rev().copied().collect();
        let sl: Vec<u8> = s.iter().rev().copied().collect();
        let what = format!("n={n} dropoff={dropoff} maxscore={maxscore} steps={steps:?}");
        let last = n.saturating_sub(1);
        // SAFETY: the n positions up from 0 and down from n - 1 are inside the vectors.
        let (right_old, right_new, left_old, left_new) = unsafe {
            (
                fields(walk::<true>(
                    tab,
                    q.as_ptr(),
                    s.as_ptr(),
                    n,
                    dropoff,
                    maxscore,
                )),
                fields(drop_walk_right(
                    tab,
                    q.as_ptr(),
                    s.as_ptr(),
                    n,
                    dropoff,
                    maxscore,
                )),
                fields(walk::<false>(
                    tab,
                    ql.as_ptr().add(last),
                    sl.as_ptr().add(last),
                    n,
                    dropoff,
                    maxscore,
                )),
                fields(drop_walk_left(
                    tab,
                    ql.as_ptr().add(last),
                    sl.as_ptr().add(last),
                    n,
                    dropoff,
                    maxscore,
                )),
            )
        };
        assert_eq!(right_new, right_old, "right walk {what}");
        assert_eq!(left_new, left_old, "left walk {what}");

        let expect = ref_right::<true>(ScoringMatrix::Blosum62, &q, &s, 0, 0, dropoff, maxscore);
        for got in [
            extend_right::<false>(tab, &q, &s, 0, 0, dropoff, maxscore),
            extend_right::<true>(tab, &q, &s, 0, 0, dropoff, maxscore),
        ] {
            assert_eq!((got.0, got.1, got.2), expect, "extend_right {what}");
            assert!(got.3 < CODES as u8);
        }
        let off = n as i32 - 1;
        let expect = ref_left::<true>(
            ScoringMatrix::Blosum62,
            &ql,
            &sl,
            off,
            off,
            dropoff,
            maxscore,
        );
        for got in [
            extend_left::<false>(tab, &ql, &sl, off, off, dropoff, maxscore),
            extend_left::<true>(tab, &ql, &sl, off, off, dropoff, maxscore),
        ] {
            assert_eq!((got.0, got.1), expect, "extend_left {what}");
            assert!(got.2 < CODES as u8);
        }

        let mut out = [(0, Exit::Ran, false); 2];
        for (k, (right, got)) in [(true, right_new), (false, left_new)]
            .into_iter()
            .enumerate()
        {
            let (best, len, end, why, near) = ncbi_walk(steps, right, dropoff, maxscore);
            assert_eq!(
                (got.0 as i64, got.1 as i64, got.2),
                (best, len, end),
                "NCBI's walk (right {right}) {what}"
            );
            out[k] = (end, why, near);
        }
        out
    }

    /// The exits a set of walks at one dropoff reached ([right, left]).
    #[derive(Debug, Default)]
    struct Seen {
        ran: [bool; 2],
        not_positive: bool,
        drop_equal: [bool; 2],
        drop_above: [bool; 2],
        near: [bool; 2],
        first_step: [bool; 2],
        n0: bool,
        n1: bool,
    }

    impl Seen {
        fn add(&mut self, n: usize, exits: [(usize, Exit, bool); 2]) {
            self.n0 |= n == 0;
            self.n1 |= n == 1;
            for (k, (end, why, near)) in exits.into_iter().enumerate() {
                match why {
                    Exit::Ran => self.ran[k] = true,
                    Exit::NotPositive => self.not_positive = true,
                    Exit::DropEqual => self.drop_equal[k] = true,
                    Exit::DropAbove => self.drop_above[k] = true,
                }
                self.near[k] |= near;
                self.first_step[k] |= why != Exit::Ran && end == 0;
            }
        }

        fn all(&self) -> bool {
            let both = |v: [bool; 2]| v[0] && v[1];
            both(self.ran)
                && self.not_positive
                && both(self.drop_equal)
                && both(self.drop_above)
                && both(self.near)
                && both(self.first_step)
                && self.n0
                && self.n1
        }
    }

    /// A residue pair (standard amino acids) with BLOSUM62 score `v`.
    fn pair_scoring(v: i32) -> (u8, u8) {
        for q in 1..=22u8 {
            for s in 1..=22u8 {
                if blosum62_score_ncbistdaa_direct(q, s) == v {
                    return (q, s);
                }
            }
        }
        panic!("no pair scores {v}");
    }

    /// Every dropoff 1..=127: walks built to stop by a drop of exactly the dropoff, to pass a
    /// drop of the dropoff less one, to step past the dropoff, to stop on `score <= 0` (right),
    /// to stop at the first step, to run out (n = 0, 1 and more), with ties (two-letter
    /// alphabets) and with random codes (all 28, the sentinel among them); each of these exits
    /// must have occurred at each dropoff.
    #[test]
    fn drop_walks_equal_walk_at_every_dropoff_and_exit() {
        let tab = build_table().expect("BLOSUM62 fits bytes");
        let (w11, a4, m1, m4) = (
            pair_scoring(11),
            pair_scoring(4),
            pair_scoring(-1),
            pair_scoring(-4),
        );
        let tie_alphabets = [
            (a4, m4),
            (pair_scoring(2), pair_scoring(-2)),
            (pair_scoring(1), m1),
            (pair_scoring(3), pair_scoring(-3)),
            (w11, m1),
            (pair_scoring(0), m1),
        ];
        let sentinel = (0u8, 1u8);
        let mut rng = Rng(0x6A09_E667_F3BC_C908);
        let mut walks = 0u64;
        for d in 1..=127i32 {
            let du = d as usize;
            let mut seen = Seen::default();
            let mut run = |steps: &[(u8, u8)], m: i32, seen: &mut Seen| {
                seen.add(steps.len(), check_walks(&tab, steps, d, m));
                walks += 2;
            };
            // a rise that keeps the right walk's score above 0 below
            let rise = vec![w11; du / 11 + 2];
            run(
                &[&rise[..], &vec![m1; du][..], &[a4; 3][..]].concat(),
                0,
                &mut seen,
            );
            run(
                &[&rise[..], &vec![m1; du - 1][..], &[a4; 2][..]].concat(),
                0,
                &mut seen,
            );
            run(
                &[&rise[..], &vec![m1; du - 1][..], &[m4][..]].concat(),
                0,
                &mut seen,
            );
            run(
                &[&vec![m1; du - 1][..], &[a4, m1, m1][..]].concat(),
                0,
                &mut seen,
            );
            // the right walk's score reaches 0 before the drop reaches the dropoff
            run(&vec![m1; du], d - 1, &mut seen);
            run(&[m1, a4], 0, &mut seen);
            // the first step stops (a sentinel; a step of -4 at dropoffs up to 4)
            run(&[sentinel, a4, a4], 0, &mut seen);
            run(&[sentinel], 50, &mut seen);
            run(&[m4, a4], 10, &mut seen);
            // n = 0 and 1, and walks that run out
            run(&[], 0, &mut seen);
            run(&[], d, &mut seen);
            run(&[a4], 0, &mut seen);
            run(&[m1], 3, &mut seen);
            run(&[sentinel], 0, &mut seen);
            let up = vec![a4; 1 + rng.below(30) as usize];
            run(&up, rng.below(100) as i32, &mut seen);
            // ties: the score comes back to the best exactly
            for _ in 0..60 {
                let (x, y) = tie_alphabets[rng.below(tie_alphabets.len() as u64) as usize];
                let len = rng.below(90) as usize;
                let steps: Vec<(u8, u8)> = (0..len)
                    .map(|_| if rng.below(2) == 0 { x } else { y })
                    .collect();
                let m = rng.below(d as u64 + 20) as i32;
                run(&steps, m, &mut seen);
            }
            // a start below 0 (the right walk then runs `walk`, the left walk does not care)
            run(&[&rise[..], &vec![m1; du][..]].concat(), -5, &mut seen);
            run(&[m1, a4, w11], -1, &mut seen);
            // random codes, the sentinel among them
            for _ in 0..40 {
                let len = rng.below(70) as usize;
                let steps: Vec<(u8, u8)> = (0..len)
                    .map(|_| (rng.below(28) as u8, rng.below(28) as u8))
                    .collect();
                let m = [0, rng.below(30) as i32, rng.below(3000) as i32][rng.below(3) as usize];
                run(&steps, m, &mut seen);
            }
            assert!(seen.all(), "dropoff {d}: not every exit reached: {seen:?}");
        }
        assert!(walks > 25_000);
    }

    /// The range in which the drop walks run their own loop (`drop_walk_proven`) and its edges:
    /// the walks there and next to it equal `walk` and the reference.
    #[test]
    fn drop_walks_at_the_edges_of_their_proven_range() {
        let tab = build_table().expect("BLOSUM62 fits bytes");
        for d in [0, 1, 127, 128] {
            assert_eq!(drop_walk_proven(false, 9, d, 0), (1..=127).contains(&d));
            assert_eq!(drop_walk_proven(true, 9, d, 0), (1..=127).contains(&d));
        }
        assert!(!drop_walk_proven(true, 9, 16, -1));
        assert!(drop_walk_proven(false, 9, 16, -1));
        let m = i32::MAX - 127 * 5;
        assert!(drop_walk_proven(true, 5, 16, m));
        assert!(!drop_walk_proven(true, 6, 16, m));
        let (w11, m1) = (pair_scoring(11), pair_scoring(-1));
        for d in [1, 16, 127] {
            for steps in [
                vec![w11; 5],
                vec![m1; 5],
                vec![m1; 6],
                vec![w11, m1, m1, w11, m1],
            ] {
                check_walks(&tab, &steps, d, m);
            }
        }
    }

    /// Random and related sequences of 0..3600 residues over all 28 codes (with and without
    /// sentinels inside), random offsets (left walks also from -1, right walks also from the
    /// end), random dropoffs 1..=127 and start scores: the extension functions with the drop
    /// walk return what they return with `walk` and what the reference functions return.
    #[test]
    fn drop_walks_equal_walk_and_the_reference_on_random_sequences() {
        let tab = build_table().expect("BLOSUM62 fits bytes");
        let mut rng = Rng(0xBB67_AE85_84CA_A73B);
        let mut walks = 0u64;
        for round in 0..300 {
            let len = match round % 5 {
                0 => rng.below(4) as usize,
                1 => 1 + rng.below(40) as usize,
                2 => 20 + rng.below(200) as usize,
                3 => 200 + rng.below(1000) as usize,
                _ => 1000 + rng.below(2600) as usize,
            };
            let lo = if round % 3 == 0 { 0u8 } else { 1 };
            let q: Vec<u8> = (0..len)
                .map(|_| lo + rng.below(28 - lo as u64) as u8)
                .collect();
            let s: Vec<u8> = if round % 2 == 0 {
                let rate = 2 + rng.below(12);
                mutate(&mut rng, &q, rate)
            } else {
                (0..len + rng.below(30) as usize)
                    .map(|_| lo + rng.below(28 - lo as u64) as u8)
                    .collect()
            };
            for _ in 0..60 {
                let d = 1 + rng.below(127) as i32;
                let m = [0, rng.below(20) as i32, rng.below(2000) as i32][rng.below(3) as usize];
                let qo = rng.below(q.len() as u64 + 1) as i32 - 1;
                let so = rng.below(s.len() as u64 + 1) as i32 - 1;
                let old = extend_left::<false>(&tab, &q, &s, qo, so, d, m);
                let new = extend_left::<true>(&tab, &q, &s, qo, so, d, m);
                let what = format!("len={len} q_off={qo} s_off={so} d={d} m={m}");
                assert_eq!(new, old, "extend_left {what}");
                let expect = ref_left::<true>(ScoringMatrix::Blosum62, &q, &s, qo, so, d, m);
                assert_eq!((new.0, new.1), expect, "extend_left {what}");
                let qo = rng.below(q.len() as u64 + 1) as i32;
                let so = rng.below(s.len() as u64 + 1) as i32;
                let old = extend_right::<false>(&tab, &q, &s, qo, so, d, m);
                let new = extend_right::<true>(&tab, &q, &s, qo, so, d, m);
                let what = format!("len={len} q_off={qo} s_off={so} d={d} m={m}");
                assert_eq!(new, old, "extend_right {what}");
                let expect = ref_right::<true>(ScoringMatrix::Blosum62, &q, &s, qo, so, d, m);
                assert_eq!((new.0, new.1, new.2), expect, "extend_right {what}");
                walks += 2;
            }
        }
        assert_eq!(walks, 36_000);
    }

    /// Two-hit and one-hit extensions with both walks at every dropoff 1..=127 on two-letter
    /// alphabets (ties) and on related sequences over all codes.
    #[test]
    fn two_hit_and_one_hit_with_both_walks_at_every_dropoff() {
        let mut rng = Rng(0x3C6E_F372_FE94_F82B);
        for d in 1..=127 {
            let (x, y) =
                [(1u8, 2u8), (1, 20), (20, 23), (4, 5), (13, 14), (17, 20)][rng.below(6) as usize];
            let len = 1 + rng.below(200) as usize;
            let mut draw = |n: usize| -> Vec<u8> {
                std::iter::once(0)
                    .chain((0..n).map(|_| if rng.below(2) == 0 { x } else { y }))
                    .chain(std::iter::once(0))
                    .collect()
            };
            let (q, s) = (draw(len), draw(len + 2));
            for _ in 0..40 {
                let w = 1 + rng.below(4) as usize;
                let qr = rng.below(q.len() as u64) as usize;
                let sr = rng.below(s.len() as u64) as usize;
                let sl = sr.saturating_sub(rng.below(6) as usize);
                check_call(&q, &s, sl, sr, qr, d, w);
            }
            let len = 50 + rng.below(400) as usize;
            let q = random_seq(&mut rng, len, 1, 28);
            let rate = 3 + rng.below(8);
            let s = mutate(&mut rng, &q, rate);
            for _ in 0..40 {
                let w = 1 + rng.below(11) as usize;
                let qr = rng.below(q.len() as u64) as usize;
                let sr = (qr + rng.below(3) as usize).min(s.len() - 1);
                let sl = sr.saturating_sub(rng.below(40) as usize);
                check_call(&q, &s, sl, sr, qr, d, w);
            }
        }
    }

    /// The OR of the codes a walk read: `or_codes` against a byte loop for every length 0..=40
    /// with a code >= 32 at every position of either sequence; then walks over sequences with
    /// such codes inside: the drop walk returns `walk`'s `Walk` (codes included), so two-hit and
    /// one-hit extensions go to the reference (`None`) exactly when they do with `walk`, and
    /// otherwise equal the reference.
    #[test]
    fn drop_walks_read_the_codes_walk_reads() {
        let tab = build_table().expect("BLOSUM62 fits bytes");
        let mut rng = Rng(0xA54F_F53A_5F1D_36F1);
        for len in 0..=40usize {
            for at in 0..len.max(1) {
                for high in [0u8, 32, 77, 128, 255] {
                    let mut q: Vec<u8> = (0..len).map(|_| rng.below(28) as u8).collect();
                    let mut s: Vec<u8> = (0..len).map(|_| rng.below(28) as u8).collect();
                    if len > 0 {
                        if rng.below(2) == 0 {
                            q[at] = high;
                        } else {
                            s[at] = high;
                        }
                    }
                    let expect = q.iter().chain(&s).fold(0u8, |a, &b| a | b);
                    // SAFETY: both vectors hold len bytes.
                    let got = unsafe { or_codes(q.as_ptr(), s.as_ptr(), len) };
                    assert_eq!(got, expect, "len={len} at={at} high={high}");
                }
            }
        }
        let (mut handled, mut fallbacks) = (0u64, 0u64);
        for _ in 0..3000 {
            let len = 1 + rng.below(80) as usize;
            let mut q = random_seq(&mut rng, len, 1, 28);
            let mut s = if rng.below(2) == 0 {
                mutate(&mut rng, &q, 4)
            } else {
                random_seq(&mut rng, len, 1, 28)
            };
            for _ in 0..1 + rng.below(3) {
                let v = if rng.below(2) == 0 { &mut q } else { &mut s };
                let at = rng.below(v.len() as u64) as usize;
                v[at] = 32 + rng.below(224) as u8;
            }
            let d = 1 + rng.below(127) as i32;
            let m = rng.below(60) as i32;
            // walks from a random start in both directions, as far as the sequences go
            let (qo, so) = (
                rng.below(q.len() as u64) as usize,
                rng.below(s.len() as u64) as usize,
            );
            let right = (q.len() - qo).min(s.len() - so);
            let left = qo.min(so) + 1;
            // SAFETY: `right` positions up and `left` positions down from (qo, so) are inside.
            unsafe {
                let (qp, sp) = (q.as_ptr().add(qo), s.as_ptr().add(so));
                assert_eq!(
                    fields(drop_walk_right(&tab, qp, sp, right, d, m)),
                    fields(walk::<true>(&tab, qp, sp, right, d, m)),
                    "right qo={qo} so={so} d={d} m={m}"
                );
                assert_eq!(
                    fields(drop_walk_left(&tab, qp, sp, left, d, m)),
                    fields(walk::<false>(&tab, qp, sp, left, d, m)),
                    "left qo={qo} so={so} d={d} m={m}"
                );
            }
            let w = 1 + rng.below(4) as usize;
            let sl = so.saturating_sub(rng.below(6) as usize);
            let two = two_hit_with::<true>(&q, &s, sl, so, qo, d, w);
            assert_eq!(two, two_hit_with::<false>(&q, &s, sl, so, qo, d, w));
            let one = one_hit_with::<true>(&q, &s, qo, so, d, w);
            assert_eq!(one, one_hit_with::<false>(&q, &s, qo, so, d, w));
            match two {
                Some(got) => {
                    let expect = extend_two_hit(ScoringMatrix::Blosum62, &q, &s, sl, so, qo, d, w);
                    assert_eq!(got, expect, "two_hit qo={qo} so={so} sl={sl} d={d} w={w}");
                    handled += 1;
                }
                None => fallbacks += 1,
            }
            match one {
                Some(got) => {
                    let expect = extend_one_hit(ScoringMatrix::Blosum62, &q, &s, qo, so, d, w);
                    assert_eq!(got, expect, "one_hit qo={qo} so={so} d={d} w={w}");
                    handled += 1;
                }
                None => fallbacks += 1,
            }
        }
        assert!(handled > 500 && fallbacks > 500, "{handled} {fallbacks}");
    }
}
