//! The walks of `s_BlastAaExtendLeft` / `s_BlastAaExtendRight` on keys that carry the score and
//! the step in one `i64` (item SC2e). blastp's and BLASTX's sentinel two-hit extension
//! (`blastp/x_ungapped_sent.rs` `two_hit_in`, which `blastx/x_ungap_sent.rs` calls) and TBLASTX's
//! (`tblastx/extension/x_two_hit_fast.rs` `two_hit_sent_on`, LOSAT_X_TBXUNGAPSENT) take them. They
//! live here, below every program, because they read only two code buffers and a table of keys.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:804-883 (the right walk) and
//! 886-921 (the left walk). A step adds the matrix score of the two codes, a score above the best
//! becomes the best (a tie keeps the earlier step), and the walk stops once the best minus the
//! score reaches the X-drop; the right walk also stops at a score of 0 or less.
//!
//! The keys: after step `i` (from 0) of a walk the key is `score << 32 - (i + 1)`, before the
//! first step `maxscore << 32`, and a step adds the matrix's key `(t << 32) - 1` for the score `t`
//! of the two codes. With `0 < i + 1 < 2^32` (a walk ends at a sentinel of a buffer of fewer than
//! `i32::MAX` bytes): `key > best` exactly when the score passes the best score (on a tie the best
//! stays at the first step, the reference's `>`); `best - key >= dropoff << 32` exactly when
//! `best_score - score >= dropoff` (the step part of the difference is in `0..2^32`); `key <= 0`
//! exactly when `score <= 0`; and the best's length is `-best & 0xffff_ffff` (`i + 1` for a best at
//! step `i`, 0 for none), its score `(best + length) >> 32` ([`best_of`]). The walk keeps no offset
//! of the best beside the key. The walks have no bound test: the callers' buffers end in sentinels
//! whose keys stop them (each caller's module comment gives the argument for its buffers).

/// The number of keys: one for each pair of codes below 32, the key of query code `q` against
/// subject code `s` at `q << 5 | s`.
pub(crate) const KEY_TABLE_LEN: usize = 32 * 32;

/// The key of a step for every pair of codes below 32 (`q << 5 | s`): `(t << 32) - 1`, `t` the
/// score of the two codes (-128 for a sentinel). The walks add it (module comment).
pub(crate) type KeyTable = [i64; KEY_TABLE_LEN];

/// The score and the length of a best key (module comment).
// No NCBI counterpart: decodes the key of the module comment.
#[inline(always)]
pub(crate) fn best_of(best: i64) -> (i32, i32) {
    let len = (best.wrapping_neg() as u64 & 0xffff_ffff) as i64;
    (((best + len) >> 32) as i32, len as i32)
}

/// A query code as the walks read it: the row of the keys (`code << 5`). blastp's block is kept
/// shifted (`u16`), BLASTX's and TBLASTX's buffers are bytes.
pub(crate) trait XQueryCode: Copy {
    /// `code << 5`.
    fn row(self) -> usize;
}

impl XQueryCode for u8 {
    #[inline(always)]
    fn row(self) -> usize {
        (self as usize) << 5
    }
}

impl XQueryCode for u16 {
    #[inline(always)]
    fn row(self) -> usize {
        self as usize
    }
}

/// One walk of `s_BlastAaExtendRight` (`RIGHT`) or `s_BlastAaExtendLeft` from `qp`, `sp`, until
/// the X-drop test (or, right, `score <= 0`) holds: (best, length, steps before the stop), on the
/// keys of the module comment.
///
/// SAFETY: walking from `qp`, `sp` in the direction, both sequences reach a sentinel (code 0)
/// with every position up to it readable and fewer than `2^32 - 1` positions before it; `keys`
/// holds the key of every pair the walk reads at `row | s` (`row` the query code's
/// [`XQueryCode::row`], `s` the subject code), the sentinel's key against every code that of a
/// score of -128; `1 <= dropoff <= 127`; and no score wraps up to the nearer sentinel.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:804-883 and 886-921 (module comment).
#[inline(always)]
pub(crate) unsafe fn walk<const RIGHT: bool, Q: XQueryCode>(
    keys: *const i64,
    qp: *const Q,
    sp: *const u8,
    dropoff: i32,
    maxscore: i32,
) -> (i32, i32, usize) {
    let mut key = (maxscore as i64) << 32;
    let mut best = key;
    let drop = (dropoff as i64) << 32;
    let mut off: isize = 0;
    let dir: isize = if RIGHT { 1 } else { -1 };
    macro_rules! step {
        ($l:lifetime) => {{
            let row = (*qp.offset(off)).row();
            let sc = *sp.offset(off) as usize;
            key += *keys.add(row | sc);
            if key > best {
                best = key;
            }
            let stop = if RIGHT {
                key <= 0 || key + drop <= best
            } else {
                key + drop <= best
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
    let (best, len) = best_of(best);
    (best, len, off.unsigned_abs())
}
