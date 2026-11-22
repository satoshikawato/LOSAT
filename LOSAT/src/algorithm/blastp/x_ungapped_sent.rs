//! EXPERIMENT (LOSAT_X_BPUNGAPSENT): blastp's two-hit ungapped extension (`s_BlastAaExtendTwoHit`
//! with the `s_BlastAaExtendLeft` / `s_BlastAaExtendRight` walks) for word size 3, inlined into the
//! candidate closure of the fused seed loop (`blast_engine.rs`), with walks that the sentinels stop
//! instead of a bound.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:804-883 (the walks) and
//! 1089-1155 (`s_BlastAaExtendTwoHit`), as in `x_ungapped_fast.rs`, whose byte
//! table this path scores with and whose statements it keeps in the same order. What changes:
//! * The checks the fast path makes on every call are made once: residue codes below 32 and the
//!   X-drops of every context in `1..=127` for a query block ([`XSentQuery::new`]), codes below 32
//!   and a sentinel at each end for a subject ([`XSentSubject::new`]), word size 3. A block or a
//!   subject that fails a check runs the reference closure.
//! * The walks have no bound test. The reference walks `n` steps at most, where `n` is the room
//!   left in both slices (right) or `min(q_off, s_off) + 1` (left); the position `n` steps away
//!   is just past the end of the slice that sets `n`, and that position holds a sentinel (code 0):
//!   the subject buffer carries one at each end (`aa_seq`), and the query block is copied with one
//!   added at its start and its end. A step on a sentinel scores -128 in the byte table; the score
//!   before a step is never above the best, so the step leaves the best as it is and makes
//!   `best - score >= 128 > dropoff`: the walk stops at step `n` with the best, its length and the
//!   stop position of the reference, which left its loop at `i = n`. Steps before `n` read the
//!   same bytes as the reference (a sentinel among them stops both, as in `x_ungapped_fast.rs`).
//! * The word loop is unrolled for three letters (the same statements and order) and the walks
//!   take two steps per loop.
//! * The score and the step are one `i64` key (item SC2e): after step `i` (from 0) of a walk the
//!   key is `score << 32 - (i + 1)`, before the first step `maxscore << 32`, and a step adds the
//!   matrix's key `(t << 32) - 1` ([`KeyTable`]). With `0 < i + 1 < 2^32` (a walk ends at a
//!   sentinel of a buffer of fewer than `i32::MAX` bytes): `key > best` exactly when the score
//!   passes the best score (on a tie the best stays at the first step, the reference's `>`);
//!   `best - key >= dropoff << 32` exactly when `best_score - score >= dropoff` (the step part of
//!   the difference is in `0..2^32`); `key <= 0` exactly when `score <= 0`; and the best's length
//!   is `-best & 0xffff_ffff` (`i + 1` for a best at step `i`, 0 for none), its score
//!   `(best + length) >> 32`. The walk keeps no offset of the best beside the key. The word loop
//!   is three such steps from 0 without a stop (its `right_d` is the length of its best).
//! * blastp's query block is kept as `u16` codes shifted by 5 (`q << 5`, [`XQueryCode`]), so that a
//!   step's index into the keys is one OR with the subject's code; BLASTX's buffers stay bytes.
//! * No score wraps (a wrapped score could step over a sentinel without stopping): a block on the
//!   path has at most `i32::MAX / 16` letters, the word (`W` letters, at most 4) and the two walks
//!   of one extension read distinct query positions at +11 a step at most (the walks' positions
//!   are distinct, the word's score only decides where they start), so a score stays below
//!   `11 * (len + 4)`, and
//!   a score is never below `best - 255` (a step costs at most 128 and the walk stops once
//!   `best - score` reaches `dropoff <= 127`).
//!
//! BLASTX's seed loops call the same extension, [`two_hit_in`] with `TAIL` 1, on buffers laid out
//! differently (`blastx/x_ungap_sent.rs`, whose header gives the argument for them).
//!
//! A matrix other than BLOSUM62 (item SC2c, [`XSentQuery::new_of`]) takes the same walks with its
//! byte table (`x_ungapped_fast::matrix_table`). The argument holds with the matrix's largest score
//! `M` for 11 and the table's -128 for the sentinel: the path takes a matrix whose table gives the
//! sentinel -128 against every code and whose `M` is at most 31, and a block of at most
//! `i32::MAX / 32` letters, so a score stays below `31 * (len + 4)`.
//!
//! In a replay of the closure on 4.3 million captured candidates of an OrthoFinder job (line F,
//! `v4-F/results/s1/p2`, arm V6) every result equalled the reference's, at 204.7 instead of 264.7
//! cycles per candidate. The key walks (item SC2e, arm K16) replayed on that capture and on the
//! candidates of BLOSUM45, BLOSUM50, BLOSUM90 and PAM250 (39 million in all) gave every result of
//! the reference at 0.89 to 0.93 of the time of the walks before them (178 instead of 199 cycles
//! per BLOSUM50 candidate, 472 instead of 648 instructions).

use super::extension::{BlastpTwoHitUngappedResult, BlastpUngappedData};
use super::x_ungapped_fast::Table;
// The key walks (`walk`, `best_of`), their keys (`KeyTable`) and the query codes they read
// (`XQueryCode`) are in algorithm/common (`x_key_walk.rs`), where TBLASTX's sentinel walks take
// them too.
use crate::algorithm::common::x_key_walk::{best_of, walk, KEY_TABLE_LEN};
pub(crate) use crate::algorithm::common::x_key_walk::{KeyTable, XQueryCode};
use std::sync::OnceLock;

// The keys (`KeyTable`, in algorithm/common) cover the codes below `CODES`.
const _: () = assert!(KEY_TABLE_LEN == (CODES as usize) * (CODES as usize));

/// The keys of the byte table `table` (its first 32 x 32 codes).
// No NCBI counterpart: the reference's scores in LOSAT's form (module comment).
pub(crate) fn key_table(table: &Table) -> Box<KeyTable> {
    let mut keys = Box::new([0i64; KEY_TABLE_LEN]);
    for (k, &t) in keys.iter_mut().zip(table.iter()) {
        *k = ((t as i64) << 32) - 1;
    }
    keys
}

/// BLOSUM62's keys, from its scores as `x_ungapped_fast::build_table` takes them (-128 for
/// `BLAST_SCORE_MIN`), built in place (no byte table beside them); `None` if a score does not fit
/// a byte.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:804-883 `matrix[q[i]][s[i]]` (the scores the reference walks add).
pub(crate) fn blosum62_keys() -> Option<Box<KeyTable>> {
    let mut keys = Box::new([0i64; KEY_TABLE_LEN]);
    for q in 0..CODES {
        for c in 0..CODES {
            let v = crate::utils::matrix::blosum62_score_ncbistdaa_direct(q, c);
            let t = if v == crate::utils::matrix::BLAST_SCORE_MIN {
                -128
            } else if (-127..=127).contains(&v) {
                v
            } else {
                return None;
            };
            keys[(q as usize) << 5 | c as usize] = ((t as i64) << 32) - 1;
        }
    }
    Some(keys)
}

/// `[0, query.., 0]` with each code shifted by 5 (module comment).
fn padded_rows(query: &[u8]) -> Vec<u16> {
    let mut padded = Vec::with_capacity(query.len() + 2);
    padded.push(0);
    padded.extend(query.iter().map(|&c| u16::from(c) << 5));
    padded.push(0);
    padded
}

/// Codes below this index the table (`q << 5 | s`).
const CODES: u8 = 32;
/// The only word size blastp's path takes ([`XSentQuery::new`]); BLASTX's takes 2 to 4, as the
/// instance `W` of [`two_hit_in`] (`blastx/x_ungap_sent.rs`).
const WORD_SIZE: usize = 3;
/// Largest X-drop the sentinel stop is exact for (see the module comment).
const MAX_DROPOFF: i32 = 127;
/// Longest query block on the path: no score of a walk wraps (see the module comment).
const MAX_BLOCK_LEN: usize = (i32::MAX / 16) as usize;
const _: () = assert!(11 * (MAX_BLOCK_LEN as i64 + 4) < i32::MAX as i64);
/// Longest query block on the path for a matrix other than BLOSUM62, whose largest score is at
/// most 31 (`MatrixTable::sent_table`): no score of a walk wraps. BLASTX's query batches of
/// another matrix take the same bound (`blastx/x_ungap_sent.rs`, item SC6b).
pub(crate) const MAX_BLOCK_LEN_OF: usize = (i32::MAX / 32) as usize;
const _: () = assert!(31 * (MAX_BLOCK_LEN_OF as i64 + 4) < i32::MAX as i64);

struct State {
    /// `LOSAT_X_BPUNGAPSENT`.
    on: bool,
    /// BLOSUM62's keys.
    keys: Option<Box<KeyTable>>,
}

fn state() -> &'static State {
    static STATE: OnceLock<State> = OnceLock::new();
    STATE.get_or_init(|| State {
        // `LOSAT_X_BPUNGAPSENT` is a strict-set switch.
        on: crate::utils::xswitch::on("LOSAT_X_BPUNGAPSENT"),
        keys: blosum62_keys(),
    })
}

/// Whether the switch is on.
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
pub(crate) fn on() -> bool {
    state().on
}

/// A query block made ready for [`two_hit`]: a copy with a sentinel added at each end.
pub(crate) struct XSentQuery {
    /// `[0, query.., 0]`, each code shifted by 5 ([`XQueryCode`]).
    padded: Vec<u16>,
    keys: &'static KeyTable,
}

impl XSentQuery {
    /// `None` when the switch is off or the block is outside the path: a code >= 32, an X-drop
    /// outside `1..=127`, a word size other than 3 (`lookup_word_size` is the lookup table's,
    /// `extension_word_size` the one the extension is called with), or more than `i32::MAX / 16`
    /// letters (no score wraps; module comment).
    // No NCBI counterpart: chooses between two LOSAT code paths with the same result.
    pub(crate) fn new(
        query: &[u8],
        x_drops: &[i32],
        lookup_word_size: i32,
        extension_word_size: usize,
    ) -> Option<Self> {
        let state = state();
        if !state.on {
            return None;
        }
        let keys = state.keys.as_deref()?;
        let fits = lookup_word_size == WORD_SIZE as i32
            && extension_word_size == WORD_SIZE
            && x_drops.iter().all(|d| (1..=MAX_DROPOFF).contains(d))
            && query.len() <= MAX_BLOCK_LEN
            && query.iter().fold(0u8, |a, &c| a | c) < CODES;
        if !fits {
            return None;
        }
        Some(Self {
            padded: padded_rows(query),
            keys,
        })
    }
}

impl XSentQuery {
    /// [`XSentQuery::new`] for a matrix other than BLOSUM62 (item SC2c): the walks score with the
    /// matrix's byte table (`x_ungapped_fast::matrix_table`), on the path only where that table
    /// gives the sentinel -128 against every code and the matrix's largest score is at most 31
    /// (`MatrixTable::sent_table`), with the same checks and at most `i32::MAX / 32` letters (no
    /// score wraps: the module comment with 31 for 11).
    // No NCBI counterpart: chooses between two LOSAT code paths with the same result.
    pub(crate) fn new_of(
        matrix: crate::config::ScoringMatrix,
        query: &[u8],
        x_drops: &[i32],
        lookup_word_size: i32,
        extension_word_size: usize,
    ) -> Option<Self> {
        if !state().on {
            return None;
        }
        let keys = keys_of(matrix)?;
        if !fits_of(
            query,
            x_drops,
            lookup_word_size,
            extension_word_size,
            MAX_BLOCK_LEN_OF,
        ) {
            return None;
        }
        Some(Self {
            padded: padded_rows(query),
            keys,
        })
    }
}

/// The keys of another matrix's table for the walks (`MatrixTable::sent_table`), built at the first
/// call for the matrix; `None` where that table is not on the path. BLASTX's walks take them too
/// (`blastx/x_ungap_sent.rs` `XSent::new_of`, item SC6b).
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
pub(crate) fn keys_of(matrix: crate::config::ScoringMatrix) -> Option<&'static KeyTable> {
    static KEYS: [OnceLock<Option<Box<KeyTable>>>; crate::config::ScoringMatrix::ALL.len()] =
        [const { OnceLock::new() }; crate::config::ScoringMatrix::ALL.len()];
    KEYS[matrix as usize]
        .get_or_init(|| {
            super::x_ungapped_fast::matrix_table(matrix)
                .and_then(|table| table.sent_table())
                .map(key_table)
        })
        .as_deref()
}

/// The checks of [`XSentQuery::new_of`] on a block of at most `max_len` letters.
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
fn fits_of(
    query: &[u8],
    x_drops: &[i32],
    lookup_word_size: i32,
    extension_word_size: usize,
    max_len: usize,
) -> bool {
    lookup_word_size == WORD_SIZE as i32
        && extension_word_size == WORD_SIZE
        && x_drops.iter().all(|d| (1..=MAX_DROPOFF).contains(d))
        && query.len() <= max_len
        && query.iter().fold(0u8, |a, &c| a | c) < CODES
}

/// A subject made ready for [`two_hit`]: its buffer with the sentinels, checked.
#[derive(Clone, Copy)]
pub(crate) struct XSentSubject<'a> {
    /// `[0, subject.., 0]`.
    aa_seq: &'a [u8],
}

impl<'a> XSentSubject<'a> {
    /// `None` unless `aa_seq` is a subject buffer with a sentinel at each end, no code >= 32
    /// inside and a length within `i32`.
    // No NCBI counterpart: chooses between two LOSAT code paths with the same result.
    #[inline]
    pub(crate) fn new(aa_seq: &'a [u8]) -> Option<Self> {
        let ok = aa_seq.len() >= 2
            && aa_seq.len() < i32::MAX as usize
            && aa_seq[0] == 0
            && aa_seq[aa_seq.len() - 1] == 0
            && aa_seq.iter().fold(0u8, |a, &c| a | c) < CODES;
        ok.then_some(Self { aa_seq })
    }

    /// The buffer, `[0, subject.., 0]`.
    #[inline]
    pub(crate) fn with_sentinels(self) -> &'a [u8] {
        self.aa_seq
    }
}

/// `s_BlastAaExtendTwoHit` for word size 3 on a checked query block and subject (`subject_raw`
/// offsets, as `extend_two_hit_blosum62` takes them): the reference's result, `None` for its
/// early return (a word that does not fit).
///
/// SAFETY: `1 <= dropoff <= 127` (an X-drop that [`XSentQuery::new`] or
/// [`XSentQuery::new_of`] checked).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1155 `s_BlastAaExtendTwoHit`.
#[inline(always)]
pub(crate) unsafe fn two_hit(
    query: &XSentQuery,
    subject: XSentSubject<'_>,
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    dropoff: i32,
) -> Option<BlastpTwoHitUngappedResult> {
    two_hit_in::<2, u16, WORD_SIZE>(
        query.keys,
        &query.padded,
        subject.aa_seq,
        s_left_off,
        s_right_off,
        q_right_off,
        dropoff,
    )
}

/// [`two_hit`] on two buffers that start with a sentinel and end with one: the reference's
/// sequences are `query[1..query.len() - TAIL + 1]` and the same of `subject` (`TAIL` 2: blastp's
/// query block and `subject_raw`, without the last sentinel; 1: BLASTX's six-frame query and
/// subject, which keep it), for the word size `W` (the reference's `word_size` argument: blastp 3,
/// BLASTX 2 to 4). A const per caller keeps each caller's instance its own, and a const `W` makes
/// the word loop's length the instance's, chosen once per subject by the caller (never per
/// candidate); nothing else in the function reads `W` (the module comment of
/// `blastx/x_ungap_sent.rs` gives the argument for `W` other than 3).
///
/// SAFETY: `1 <= dropoff <= 127`; both buffers have a sentinel (code 0) at index 0 and at their
/// last index; `table` scores the sentinel -128 against every code, both ways, and with its
/// largest score `M`, `M * (len + 4) < i32::MAX` for the `len` letters of `query` (BLOSUM62's 11
/// and at most `i32::MAX / 16` letters, or `M <= 31` and at most `i32::MAX / 32`: module
/// comment); fewer than `i32::MAX` bytes in `subject`; `TAIL` is 1 or 2; and `W` is 1 to 4.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1155 `s_BlastAaExtendTwoHit`.
#[inline(always)]
pub(crate) unsafe fn two_hit_in<const TAIL: usize, Q: XQueryCode, const W: usize>(
    keys: &KeyTable,
    query: &[Q],
    subject: &[u8],
    s_left_off: usize,
    s_right_off: usize,
    q_right_off: usize,
    dropoff: i32,
) -> Option<BlastpTwoHitUngappedResult> {
    debug_assert!((1..=MAX_DROPOFF).contains(&dropoff));
    debug_assert!(TAIL == 1 || TAIL == 2);
    debug_assert!((1..=4).contains(&W));
    let qlen = query.len() - TAIL;
    let slen = subject.len() - TAIL;
    // The early return of the reference wrapper (a word that does not fit).
    if q_right_off + W > qlen || s_right_off + W > slen {
        return None;
    }
    let tab = keys.as_ptr();
    // offset 0 of the reference's sequences; one sentinel before each
    let qb = query.as_ptr().add(1);
    let sb = subject.as_ptr().add(1);
    let mut q_right_off = q_right_off as i32;
    let s_left_off = s_left_off as i32;
    let mut s_right_off = s_right_off as i32;

    // the word loop of `W` letters (inside both sequences: checked above), on keys from 0: its
    // best is `left_score` and the best's length `right_d` (module comment)
    let mut key = 0i64;
    let mut best = 0i64;
    let qw = qb.add(q_right_off as usize);
    let sw = sb.add(s_right_off as usize);
    for i in 0..W {
        key += *tab.add((*qw.add(i)).row() | *sw.add(i) as usize);
        if key > best {
            best = key;
        }
    }
    let right_d = best_of(best).1 as i32;
    q_right_off += right_d;
    s_right_off += right_d;

    // left: from (q_right_off - 1, s_right_off - 1) down; no step when either is -1
    let (left_score, left_d) = if q_right_off == 0 || s_right_off == 0 {
        (0, 0)
    } else {
        // SAFETY: the query block and the subject have a sentinel before offset 0.
        let (best, len, _) = walk::<false, Q>(
            tab,
            qb.offset(q_right_off as isize - 1),
            sb.offset(s_right_off as isize - 1),
            dropoff,
            0,
        );
        (best, len as i32)
    };

    let mut right_extend = false;
    let mut s_last_off = s_right_off;
    let mut right_score = 0i32;
    let mut right_d = 0i32;
    if left_d >= (s_right_off - s_left_off) {
        right_extend = true;
        // SAFETY: with TAIL 2, q_right_off <= qlen and s_right_off <= slen, and both have a
        // sentinel at offset qlen / slen. With TAIL 1 the last sentinel is at offset qlen - 1 /
        // slen - 1, and the walk starts at or before it: the word ends at that offset at most, and
        // a word whose last letter is a sentinel has right_d < W (that letter's step scores -128,
        // so its key is below the key before it, which is at most the best: a sentinel step is
        // never the best, whatever `W` and the other letters' scores, and the table gives the
        // sentinel -128 against every code), so the start is past it only when it is not in the
        // word.
        let (best, len, end) = walk::<true, Q>(
            tab,
            qb.add(q_right_off as usize),
            sb.add(s_right_off as usize),
            dropoff,
            left_score,
        );
        (right_score, right_d, s_last_off) = (best, len as i32, s_right_off + end as i32);
    }

    Some(BlastpTwoHitUngappedResult {
        ungapped_data: BlastpUngappedData {
            q_start: q_right_off - left_d,
            s_start: s_right_off - left_d,
            length: left_d + right_d,
            score: left_score.max(right_score),
        },
        right_extend,
        s_last_off,
    })
}

#[cfg(test)]
#[path = "x_ungapped_sent_tests.rs"]
mod tests;
