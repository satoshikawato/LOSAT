//! EXPERIMENT (LOSAT_X_BXUNGAPSENT / LOSAT_X_BXUNGAPSENTSHADOW): the BLASTX two-hit ungapped
//! extension of the seed loops (the fused loop of `x_seed_diag.rs` and the direct loop of
//! `x_seed_scan.rs`) as blastp's LOSAT_X_BPUNGAPSENT makes it (`blastp/x_ungapped_sent.rs`): the
//! walks stopped by the sentinels instead of a bound, the checks made once per query batch and
//! subject, and the extension inlined into the loop.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:804-883 (the walks) and
//! 1089-1155 (`s_BlastAaExtendTwoHit`), which the seed loops reach through
//! `extend_two_hit_blosum62` (and LOSAT_X_BPUNGAP's `x_ungapped_fast::two_hit`).
//!
//! Why the result is the reference's. The extension is `x_ungapped_sent::two_hit_in::<1>`, whose
//! statements are blastp's; what differs is where the sentinels are.
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
//! * Checks: the word size 3, a two-hit window, every code of the query batch below 32, at most
//!   `i32::MAX / 16` letters and the sentinels at both ends, once per batch ([`query_fits`]); the
//!   subject's codes and sentinels once per batch and subject (`XSentSubject::new`); the X-drop of
//!   each candidate's context in `1..=127` (on each candidate: the X-drops are the contexts'
//!   parameters, which the batch does not hold; a candidate outside takes the reference, out of
//!   line). A batch or a subject that fails a check runs the reference.

use crate::algorithm::blastp::extension::{extend_two_hit_blosum62, BlastpTwoHitUngappedResult};
use crate::algorithm::blastp::x_ungapped_fast::{build_table, Table};
use crate::algorithm::blastp::x_ungapped_sent::{two_hit_in, XSentSubject};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

/// Codes below this index the table (`q << 5 | s`).
const CODES: u8 = 32;
/// The only word size this path takes.
const WORD_SIZE: i32 = 3;
/// Largest X-drop the sentinel stop is exact for (`x_ungapped_sent.rs`).
const MAX_DROPOFF: i32 = 127;
/// Longest query batch on the path (no score of a walk wraps, `x_ungapped_sent.rs`).
const MAX_QUERY_LEN: usize = (i32::MAX / 16) as usize;

struct State {
    /// 0 = off, 1 = LOSAT_X_BXUNGAPSENT, 2 = LOSAT_X_BXUNGAPSENTSHADOW (both, compared).
    mode: u8,
    table: Option<Box<Table>>,
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
            table: if mode != 0 { build_table() } else { None },
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
/// is on, the word size is 3 (not `-word_size 5`'s compressed table), the search has a two-hit
/// window (not `-window_size 0`'s one-hit extension), and the buffer has a sentinel at each end,
/// at most `i32::MAX / 16` letters and no code >= 32.
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
pub(crate) fn query_fits(buffer: &[u8], word_size: i32, window_size: i32) -> bool {
    on() && word_size == WORD_SIZE && window_size > 0 && query_on_path(buffer)
}

/// [`query_fits`] without the switch.
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
fn query_on_path(buffer: &[u8]) -> bool {
    buffer.len() >= 2
        && buffer.len() - 2 <= MAX_QUERY_LEN
        && buffer[0] == 0
        && buffer[buffer.len() - 1] == 0
        && buffer.iter().fold(0u8, |a, &c| a | c) < CODES
}

static SHADOW_SUBJECTS: AtomicU64 = AtomicU64::new(0);
static SHADOW_COMPARED: AtomicU64 = AtomicU64::new(0);

/// The summary line (printed by `main` at exit), only in shadow mode and after a subject took the
/// path: `subjects=` counts the pairs of a query batch and a subject made ready for the path (the
/// seed loop then extends with it, or, with `LOSAT_X_BXSCAN=0` where the fused loop does not act,
/// runs `word_finder`), `compared=` the extensions compared.
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

/// One subject of one query batch made ready for [`XSent::two_hit`].
#[derive(Clone, Copy)]
pub(crate) struct XSent<'a> {
    /// The batch's buffer, `[0, frames.., 0]`.
    query: &'a [u8],
    /// `[0, residues.., 0]`.
    subject: XSentSubject<'a>,
    table: &'static Table,
    shadow: bool,
}

impl<'a> XSent<'a> {
    /// The pair of a query batch that [`query_fits`] passed (its buffer `query`) and a subject
    /// buffer (`with_sentinels`, the residues between two sentinels); `None` when the switch is off
    /// or `XSentSubject::new` refuses the subject.
    // No NCBI counterpart: chooses between two LOSAT code paths with the same result.
    #[inline]
    pub(crate) fn new(query: &'a [u8], with_sentinels: &'a [u8]) -> Option<Self> {
        let state = state();
        let table = state.table.as_deref()?;
        let subject = XSentSubject::new(with_sentinels)?;
        if state.mode == 2 {
            SHADOW_SUBJECTS.fetch_add(1, Ordering::Relaxed);
        }
        Some(Self {
            query,
            subject,
            table,
            shadow: state.mode == 2,
        })
    }

    /// A pair for the tests of the seed loops, with or without the shadow's comparison.
    #[cfg(test)]
    pub(crate) fn for_test(query: &'a [u8], with_sentinels: &'a [u8], shadow: bool) -> Self {
        static TABLE: OnceLock<Box<Table>> = OnceLock::new();
        assert!(query_on_path(query));
        Self {
            query,
            subject: XSentSubject::new(with_sentinels).expect("a subject on the path"),
            table: TABLE.get_or_init(|| build_table().expect("BLOSUM62 fits a byte table")),
            shadow,
        }
    }

    /// LOSAT_X_BXUNGAPSENTSHADOW: run the reference extension too and compare (the callers
    /// take `two_hit::<true>`).
    #[inline]
    pub(crate) fn shadow(&self) -> bool {
        self.shadow
    }

    /// `extend_two_hit_blosum62(query, subject, ..)` for the reference's sequences (the buffers
    /// without their first byte) and word size 3: the sentinel walks for an X-drop in `1..=127`,
    /// the reference otherwise; with `SHADOW` both, compared (a difference stops the run).
    // NCBI reference: aa_ungapped.c:1089-1155 `s_BlastAaExtendTwoHit` (module header).
    #[inline(always)]
    pub(crate) fn two_hit<const SHADOW: bool>(
        &self,
        s_left_off: usize,
        s_right_off: usize,
        q_right_off: usize,
        dropoff: i32,
    ) -> Option<BlastpTwoHitUngappedResult> {
        if !(1..=MAX_DROPOFF).contains(&dropoff) {
            return self.reference(s_left_off, s_right_off, q_right_off, dropoff);
        }
        // SAFETY: the dropoff is in 1..=127 (above); the query buffer passed `query_fits` and
        // the subject buffer `XSentSubject::new` (sentinels at both ends, the lengths).
        let got = unsafe {
            two_hit_in::<1>(
                self.table,
                self.query,
                self.subject.with_sentinels(),
                s_left_off,
                s_right_off,
                q_right_off,
                dropoff,
            )
        };
        if SHADOW {
            let expected = self.reference(s_left_off, s_right_off, q_right_off, dropoff);
            assert_eq!(
                got, expected,
                "[X_SHADOW] BXUNGAPSENT two-hit mismatch: s_left_off={s_left_off} s_right_off={s_right_off} q_right_off={q_right_off} dropoff={dropoff}"
            );
            SHADOW_COMPARED.fetch_add(1, Ordering::Relaxed);
        }
        got
    }

    /// `extend_two_hit_blosum62` on the reference's sequences, out of line (a candidate whose
    /// X-drop is above 127, and the shadow's comparison).
    // NCBI reference: aa_ungapped.c:1089-1155 `s_BlastAaExtendTwoHit`.
    #[cold]
    #[inline(never)]
    fn reference(
        &self,
        s_left_off: usize,
        s_right_off: usize,
        q_right_off: usize,
        dropoff: i32,
    ) -> Option<BlastpTwoHitUngappedResult> {
        extend_two_hit_blosum62(
            &self.query[1..],
            self.subject_raw(),
            s_left_off,
            s_right_off,
            q_right_off,
            dropoff,
            WORD_SIZE as usize,
        )
    }

    /// The subject as the seed loops take it: the residues and the last sentinel.
    #[inline]
    fn subject_raw(&self) -> &'a [u8] {
        &self.subject.with_sentinels()[1..]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastp::extension::extend_two_hit;
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

    /// A residue code: mostly 1..=24, sometimes the stop codon 25, a code 26..=31, or (with
    /// `zero`) an inner 0.
    fn code(rng: &mut Rng, zero: bool) -> u8 {
        match rng.below(40) {
            0 => 25,
            1 => 26 + rng.below(6) as u8,
            2 if zero => 0,
            _ => 1 + rng.below(24) as u8,
        }
    }

    /// A six-frame query buffer as `query_setup.rs` lays it out: `[0, frame, 0, frame, .., 0]`,
    /// frames of 1 letter and more, some copied from `subject` with substitutions so that walks
    /// run long.
    fn query_buffer(rng: &mut Rng, subject: &[u8]) -> Vec<u8> {
        let mut buffer = vec![0u8];
        for _ in 0..6 {
            let most = if rng.below(4) == 0 { 3 } else { 90 };
            let len = 1 + rng.below(most) as usize;
            if rng.below(2) == 0 && subject.len() > 2 {
                let from = 1 + rng.below((subject.len() - 2) as u64) as usize;
                for i in 0..len {
                    let c = subject[1 + (from + i) % (subject.len() - 2)];
                    let c = if rng.below(6) == 0 {
                        code(rng, false)
                    } else {
                        c
                    };
                    buffer.push(c);
                }
            } else {
                for _ in 0..len {
                    buffer.push(code(rng, false));
                }
            }
            buffer.push(0);
        }
        buffer
    }

    /// The path gives NCBI's extension (`extend_two_hit`, the reference without LOSAT_X_BPUNGAP)
    /// on the reference's sequences (the buffers without their first byte) for every offset of a
    /// word, the ends of both sequences and the frame borders included, subjects of 0 to 3
    /// letters, inner sentinels, codes up to 31, and X-drops on and off the path (above 127: the
    /// reference).
    #[test]
    fn two_hit_is_the_reference_on_blastx_buffers() {
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        let mut compared = 0u64;
        for round in 0..400 {
            let slen = rng.below(if round % 8 == 0 { 4 } else { 120 }) as usize;
            let mut subject = vec![0u8];
            subject.extend((0..slen).map(|_| code(&mut rng, true)));
            subject.push(0);
            let query = query_buffer(&mut rng, &subject);
            let sent = XSent::for_test(&query, &subject, false);
            let (q_ref, s_ref) = (&query[1..], &subject[1..]);
            for _ in 0..200 {
                let q = rng.below(q_ref.len() as u64 + 2) as usize;
                let s = rng.below(s_ref.len() as u64 + 2) as usize;
                let s_left = rng.below(s as u64 + 4) as usize;
                let dropoff = [1, 7, 16, 40, 127, 128, 300][rng.below(7) as usize];
                let want = extend_two_hit(
                    ScoringMatrix::Blosum62,
                    q_ref,
                    s_ref,
                    s_left,
                    s,
                    q,
                    dropoff,
                    3,
                );
                let got = sent.two_hit::<false>(s_left, s, q, dropoff);
                assert_eq!(
                    got, want,
                    "round {round}: s_left={s_left} s={s} q={q} dropoff={dropoff}"
                );
                compared += u64::from(want.is_some());
            }
        }
        assert!(compared > 10_000, "only {compared} extensions");
    }

    /// The batch check: sentinels at both ends and codes below 32 (`query_fits` adds the switch,
    /// the word size and the window).
    #[test]
    fn query_on_path_checks_the_buffer() {
        assert!(query_on_path(&[0, 3, 4, 0, 5, 0]));
        assert!(query_on_path(&[0, 0]));
        assert!(!query_on_path(&[0]));
        assert!(!query_on_path(&[]));
        assert!(!query_on_path(&[1, 3, 0]));
        assert!(!query_on_path(&[0, 3, 1]));
        assert!(!query_on_path(&[0, 32, 0]));
        assert!(query_on_path(&[0, 31, 0]));
    }
}
