//! EXPERIMENT (LOSAT_X_BXUNGAPSENT / LOSAT_X_BXUNGAPSENTSHADOW): the BLASTX two-hit ungapped
//! extension of the seed loops (the fused loop of `x_seed_diag.rs` and the direct loop of
//! `x_seed_scan.rs`) as blastp's LOSAT_X_BPUNGAPSENT makes it (`blastp/x_ungapped_sent.rs`): the
//! walks stopped by the sentinels instead of a bound, the checks made once per query batch and
//! subject, and the extension inlined into the loop.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:804-883 (the walks) and
//! 1089-1155 (`s_BlastAaExtendTwoHit`), which the seed loops reach through
//! `extend_two_hit_blosum62` (and `x_ungapped_fast::two_hit`), another matrix's through
//! `extend_two_hit_table` (`seed.rs` `bx_two_hit`).
//!
//! Why the result is the reference's. The extension is `x_ungapped_sent::two_hit_in::<1, u8, W>`,
//! whose statements are blastp's; what differs is where the sentinels are, and `W`, the word size
//! (2, 3 or 4; below).
//! * The query is the batch's buffer (`PreparedQueryBatch::sequence_start`): the six frames of each
//!   query, each between two sentinels (code 0), the reference's query being the buffer without its
//!   first byte (it keeps the last sentinel). The subject is its residues between two sentinels,
//!   the reference's subject being the buffer without its first byte (it keeps the last sentinel,
//!   as the seed loops take it).
//! * Right walk: the reference walks `n` steps at most, `n` the room left in both of its
//!   sequences; each of them ends with a sentinel, which the walk reaches at step `n - 1` at the
//!   latest, and a step on a sentinel stops a walk (`x_ungapped_sent.rs`), so the bound never
//!   acts and the walk without it is the same walk.
//! * Left walk: the reference walks `min(q_off, s_off) + 1` steps; the position one step further is
//!   the first byte of a buffer, a sentinel, which stops the walk without changing its result (the
//!   argument of `x_ungapped_sent.rs`).
//! * Checks: the word size 2 to 4, a two-hit window, every code of the query batch below 32, at most
//!   `i32::MAX / 16` letters and the sentinels at both ends, once per batch ([`query_fits`]); the
//!   subject's codes and sentinels once per batch and subject (`XSentSubject::new`); the X-drop of
//!   each candidate's context in `1..=127` (on each candidate: the X-drops are the contexts'
//!   parameters, which the batch does not hold; a candidate outside takes the reference, out of
//!   line). A batch or a subject that fails a check runs the reference.
//!
//! The word size (item SC9b-1 of line NP). `W` is a const of the walk's instance (`two_hit::<_, W>`,
//! `two_hit_in::<1, u8, W>`) for the word sizes 2, 3 and 4, and the seed loops (`x_seed_diag.rs`,
//! `x_seed_scan.rs`) take the instance `options.word_size` names once per subject, never per
//! candidate; the loops' own statements still read the run-time word size, so the instance for 3
//! is the loop of the word size 3 before this item. `W` appears in four places of the extension,
//! and in none of them does the argument above or its arguments for the sentinel stop depend on
//! its value:
//! * The reference's word loop (aa_ungapped.c:1106-1118, `for (i = 0; i < word_size; i++)`) is
//!   `W` steps from 0 without a stop, with `score > left_score` for its best (the first of equal
//!   scores): the walk's keys give that on `W` steps (`key > best`; a tie keeps the earlier step
//!   for any `W < 2^32`, as for the walks), and its `right_d` is the best's length, 0 when no
//!   step scores above 0. The word is inside both sequences by the reference wrapper's early
//!   return (`q_right_off + W > qlen`, the same test).
//! * The right walk starts in the buffer. With the subject and the query ending in the sentinel
//!   that the reference keeps (`TAIL` 1), the word may end on it, and then the right walk starts
//!   at `q_right_off + right_d`, which must not pass it. The word's last letter being a sentinel,
//!   its step adds the key `(-128 << 32) - 1`, which is below the key before it, and the best is
//!   at least that key (it is the largest key so far, 0 before the first step): a sentinel step
//!   is never the best, so `right_d < W` and the start is at the sentinel at most. That does not
//!   depend on `W` or on the other letters of the word (sentinels inside it included; the keys
//!   are blastp's, with -128 for a sentinel against every code, both ways).
//! * The two walks read distinct positions: the left walk those before `q_right_off` (already
//!   moved by `right_d`) and the right walk those at and after it; the word loop's key only
//!   decides `right_d`, its score is not added to a walk's. So the walks' scores stay below
//!   `M * len` for the largest score `M` (11, or at most 31 for another matrix), the word loop's
//!   below `M * W`, and `M * (len + 4)` bounds all of them for `W <= 4` (the bound of
//!   `x_ungapped_sent.rs`, whose `+ 3` was the word of 3).
//! * The reference (the candidate outside `1..=127`, the shadow) is called with the same `W`
//!   ([`XSent::two_hit`] passes it on to `extend_two_hit`'s `word_size`).
//!
//! Another matrix than BLOSUM62 (item SC6b of line NP, `XSent<MatrixBx>`, [`query_fits_of`]) takes
//! the same walks with the keys of its byte table (blastp's `x_ungapped_sent::keys_of`, item SC2c),
//! which exist only where that table gives the sentinel -128 against every code and the matrix's
//! largest score is at most 31; the argument above then holds with 31 for BLOSUM62's 11 (blastp's
//! module comment), so such a batch is on the path with at most `i32::MAX / 32` letters. Its
//! reference (a candidate whose X-drop is outside `1..=127`, and the shadow's comparison) is the
//! matrix's (`extend_two_hit`).

use super::BxScoring;
use crate::algorithm::blastp::extension::{
    extend_two_hit, extend_two_hit_blosum62, BlastpTwoHitUngappedResult,
};
use crate::algorithm::blastp::x_ungapped_sent::{
    blosum62_keys, keys_of, two_hit_in, KeyTable, XSentSubject,
};
use crate::config::ScoringMatrix;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// Codes below this index the table (`q << 5 | s`).
const CODES: u8 = 32;
/// The smallest and the largest word size this path takes (module header, item SC9b-1); the
/// word sizes above them, the compressed tables of 5 to 7, run the reference.
const MIN_WORD_SIZE: i32 = 2;
const MAX_WORD_SIZE: i32 = 4;
/// Largest X-drop the sentinel stop is exact for (`x_ungapped_sent.rs`).
const MAX_DROPOFF: i32 = 127;
/// Longest query batch on the path (no score of a walk wraps, `x_ungapped_sent.rs`).
const MAX_QUERY_LEN: usize = (i32::MAX / 16) as usize;
/// Longest query batch on the path for another matrix than BLOSUM62, whose largest score is at most
/// 31 (module header): no score of a walk wraps (blastp's bound for its blocks).
const MAX_QUERY_LEN_OF: usize = crate::algorithm::blastp::x_ungapped_sent::MAX_BLOCK_LEN_OF;

struct State {
    /// 0 = off, 1 = LOSAT_X_BXUNGAPSENT, 2 = LOSAT_X_BXUNGAPSENTSHADOW (both, compared).
    mode: u8,
    /// BLOSUM62's keys for the walks (`x_ungapped_sent::KeyTable`).
    keys: Option<Box<KeyTable>>,
}

fn state() -> &'static State {
    static STATE: OnceLock<State> = OnceLock::new();
    STATE.get_or_init(|| {
        // `LOSAT_X_BXUNGAPSENT` is a strict-set switch; the shadow variable is read first, as
        // "set to any value", and implies the switch.
        let mode = if std::env::var_os("LOSAT_X_BXUNGAPSENTSHADOW").is_some() {
            2
        } else if crate::utils::xswitch::on("LOSAT_X_BXUNGAPSENT") {
            1
        } else {
            0
        };
        State {
            mode,
            keys: if mode != 0 { blosum62_keys() } else { None },
        }
    })
}

/// Whether the switch (or its shadow) is on.
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
#[inline]
pub(crate) fn on() -> bool {
    state().mode != 0
}

/// Whether a query batch is on the path (computed once per batch, `query_setup.rs`): the switch
/// is on, the word size is 2, 3 or 4 (not the compressed table of 5 to 7), the search has a two-hit
/// window (not `-window_size 0`'s one-hit extension), and the buffer has a sentinel at each end,
/// at most `i32::MAX / 16` letters and no code >= 32.
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
pub(crate) fn query_fits(buffer: &[u8], word_size: i32, window_size: i32) -> bool {
    on() && word_size_on_path(word_size) && window_size > 0 && query_on_path(buffer)
}

/// Whether the walks have an instance for the word size (`two_hit::<_, W>`, module header).
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
#[inline(always)]
fn word_size_on_path(word_size: i32) -> bool {
    (MIN_WORD_SIZE..=MAX_WORD_SIZE).contains(&word_size)
}

/// [`query_fits`] for another matrix than BLOSUM62 (module header): also its keys exist, and the
/// buffer has at most `i32::MAX / 32` letters.
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
pub(crate) fn query_fits_of(
    matrix: ScoringMatrix,
    buffer: &[u8],
    word_size: i32,
    window_size: i32,
) -> bool {
    on() && word_size_on_path(word_size)
        && window_size > 0
        && keys_of(matrix).is_some()
        && query_on_path_within(buffer, MAX_QUERY_LEN_OF)
}

/// [`query_fits`] without the switch.
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
fn query_on_path(buffer: &[u8]) -> bool {
    query_on_path_within(buffer, MAX_QUERY_LEN)
}

/// [`query_on_path`] with at most `max_len` letters.
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
#[inline(always)]
fn query_on_path_within(buffer: &[u8], max_len: usize) -> bool {
    buffer.len() >= 2
        && buffer.len() - 2 <= max_len
        && buffer[0] == 0
        && buffer[buffer.len() - 1] == 0
        && buffer.iter().fold(0u8, |a, &c| a | c) < CODES
}

/// Runs `$body` with the const `$w` set to the word size `$word` (2, 3 or 4), the instance of the
/// walks (`two_hit::<_, W>`) that a seed loop takes for a subject: the seed loops use it once per
/// subject, outside their candidate loops (module header, item SC9b-1). A word size that
/// [`query_fits`] does not take cannot reach it (`XSent` exists only for a batch that passed).
// No NCBI counterpart: chooses between instances of one LOSAT code path with the same result.
macro_rules! by_word_size {
    ($word:expr, $w:ident => $body:expr) => {
        match $word {
            2 => {
                const $w: usize = 2;
                $body
            }
            3 => {
                const $w: usize = 3;
                $body
            }
            4 => {
                const $w: usize = 4;
                $body
            }
            word => unreachable!("word size {word} is not on the sentinel walks' path"),
        }
    };
}
pub(crate) use by_word_size;

static SHADOW_SUBJECTS: AtomicU64 = AtomicU64::new(0);
static SHADOW_COMPARED: AtomicU64 = AtomicU64::new(0);

/// The summary line (printed by `main` at exit), only in shadow mode and after a subject took the
/// path: `subjects=` counts the pairs of a query batch and a subject made ready for the path,
/// `compared=` the extensions compared.
// No NCBI counterpart: prints a counter; it does not change any value NCBI computes.
pub(crate) fn print_shadow_stats() {
    let subjects = SHADOW_SUBJECTS.load(Ordering::Relaxed);
    if subjects > 0 && state().mode == 2 {
        eprintln!(
            "[X_SHADOW] BXUNGAPSENT subjects={subjects} compared={} (two-hit ungapped extensions identical)",
            SHADOW_COMPARED.load(Ordering::Relaxed)
        );
    }
}

/// One subject of one query batch made ready for [`XSent::two_hit`], with the matrix of the walks
/// as a type (`BxScoring`): BLOSUM62's (`Blosum62Bx`, no bytes) or another matrix's (`MatrixBx`,
/// item SC6b), whose reference walks a candidate off the path ([`XSent::reference`]).
#[derive(Clone, Copy)]
pub(crate) struct XSent<'a, S: BxScoring> {
    /// The batch's buffer, `[0, frames.., 0]`.
    query: &'a [u8],
    /// `[0, residues.., 0]`.
    subject: XSentSubject<'a>,
    keys: &'static KeyTable,
    shadow: bool,
    scoring: S,
}

impl<'a, S: BxScoring> XSent<'a, S> {
    /// The pair of a query batch that [`query_fits`] (another matrix: [`query_fits_of`]) passed (its
    /// buffer `query`) and a subject buffer (`with_sentinels`, the residues between two sentinels),
    /// with the keys of `scoring`'s matrix; `None` when the switch is off, the matrix has no keys or
    /// `XSentSubject::new` refuses the subject.
    // No NCBI counterpart: chooses between two LOSAT code paths with the same result.
    #[inline]
    pub(crate) fn new(scoring: S, query: &'a [u8], with_sentinels: &'a [u8]) -> Option<Self> {
        let state = state();
        let keys = if S::BLOSUM62 {
            state.keys.as_deref()?
        } else if state.mode == 0 {
            return None;
        } else {
            keys_of(scoring.matrix())?
        };
        let subject = XSentSubject::new(with_sentinels)?;
        if state.mode == 2 {
            SHADOW_SUBJECTS.fetch_add(1, Ordering::Relaxed);
        }
        Some(Self {
            query,
            subject,
            keys,
            shadow: state.mode == 2,
            scoring,
        })
    }

    /// A pair for the tests of the seed loops, with or without the shadow's comparison, with the
    /// keys of `scoring`'s matrix (one that has them).
    #[cfg(test)]
    pub(crate) fn for_test(
        scoring: S,
        query: &'a [u8],
        with_sentinels: &'a [u8],
        shadow: bool,
    ) -> Self {
        static KEYS: OnceLock<Box<KeyTable>> = OnceLock::new();
        let keys: &'static KeyTable = if S::BLOSUM62 {
            assert!(query_on_path(query));
            KEYS.get_or_init(|| blosum62_keys().expect("BLOSUM62 fits a byte table"))
        } else {
            assert!(query_on_path_within(query, MAX_QUERY_LEN_OF));
            keys_of(scoring.matrix()).expect("the matrix's table is on the path")
        };
        Self {
            query,
            subject: XSentSubject::new(with_sentinels).expect("a subject on the path"),
            keys,
            shadow,
            scoring,
        }
    }

    /// LOSAT_X_BXUNGAPSENTSHADOW: run the reference extension too and compare (the callers
    /// take `two_hit::<true>`).
    #[inline]
    pub(crate) fn shadow(&self) -> bool {
        self.shadow
    }

    /// `extend_two_hit_blosum62(query, subject, ..)` (with another matrix, `extend_two_hit` with it)
    /// for the reference's sequences (the buffers without their first byte) and the word size `W`
    /// (2 to 4; the caller's instance for `options.word_size`, chosen once per subject): the
    /// sentinel walks for an X-drop in `1..=127`, the reference otherwise; with `SHADOW` both,
    /// compared (a difference stops the run).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1155 `s_BlastAaExtendTwoHit` (module header).
    #[inline(always)]
    pub(crate) fn two_hit<const SHADOW: bool, const W: usize>(
        &self,
        s_left_off: usize,
        s_right_off: usize,
        q_right_off: usize,
        dropoff: i32,
    ) -> Option<BlastpTwoHitUngappedResult> {
        if !(1..=MAX_DROPOFF).contains(&dropoff) {
            return self.reference::<W>(s_left_off, s_right_off, q_right_off, dropoff);
        }
        // SAFETY: the dropoff is in 1..=127 (above); the query buffer passed `query_fits` and
        // the subject buffer `XSentSubject::new` (sentinels at both ends, the lengths); `W` is
        // the word size that `query_fits` took, 2 to 4.
        let got = unsafe {
            two_hit_in::<1, u8, W>(
                self.keys,
                self.query,
                self.subject.with_sentinels(),
                s_left_off,
                s_right_off,
                q_right_off,
                dropoff,
            )
        };
        if SHADOW {
            let expected = self.reference::<W>(s_left_off, s_right_off, q_right_off, dropoff);
            assert_eq!(
                got, expected,
                "[X_SHADOW] BXUNGAPSENT two-hit mismatch: s_left_off={s_left_off} s_right_off={s_right_off} q_right_off={q_right_off} dropoff={dropoff}"
            );
            SHADOW_COMPARED.fetch_add(1, Ordering::Relaxed);
        }
        got
    }

    /// `extend_two_hit_blosum62` (another matrix: `extend_two_hit` with it) on the reference's
    /// sequences and the word size `W`, out of line (a candidate whose X-drop is above 127, and the
    /// shadow's comparison).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1155 `s_BlastAaExtendTwoHit`.
    #[cold]
    #[inline(never)]
    fn reference<const W: usize>(
        &self,
        s_left_off: usize,
        s_right_off: usize,
        q_right_off: usize,
        dropoff: i32,
    ) -> Option<BlastpTwoHitUngappedResult> {
        if !S::BLOSUM62 {
            return extend_two_hit(
                self.scoring.matrix(),
                &self.query[1..],
                self.subject_raw(),
                s_left_off,
                s_right_off,
                q_right_off,
                dropoff,
                W,
            );
        }
        extend_two_hit_blosum62(
            &self.query[1..],
            self.subject_raw(),
            s_left_off,
            s_right_off,
            q_right_off,
            dropoff,
            W,
        )
    }

    /// The subject as the seed loops take it: the residues and the last sentinel.
    #[inline]
    fn subject_raw(&self) -> &'a [u8] {
        &self.subject.with_sentinels()[1..]
    }
}

#[cfg(test)]
#[path = "x_ungap_sent_tests.rs"]
mod tests;
