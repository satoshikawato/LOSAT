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
//! * Checks: every code of the query batch below 32, at most `i32::MAX / 16` letters, and the
//!   sentinels at both ends, once per batch ([`query_fits`]); the subject's codes and sentinels
//!   once per subject (`XSentSubject::new`); the word size 3; the X-drop of each candidate's
//!   context in `1..=127` (on each candidate, which otherwise takes the reference). A batch or a
//!   subject that fails a check runs the reference.

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

/// Whether a query batch's buffer is on the path (computed once per batch, `query_setup.rs`): the
/// switch is on, the buffer has a sentinel at each end, at most `i32::MAX / 16` letters and no
/// code >= 32.
// No NCBI counterpart: chooses between two LOSAT code paths with the same result.
pub(crate) fn query_fits(buffer: &[u8]) -> bool {
    on() && buffer.len() >= 2
        && buffer.len() - 2 <= MAX_QUERY_LEN
        && buffer[0] == 0
        && buffer[buffer.len() - 1] == 0
        && buffer.iter().fold(0u8, |a, &c| a | c) < CODES
}

static SHADOW_SUBJECTS: AtomicU64 = AtomicU64::new(0);
static SHADOW_COMPARED: AtomicU64 = AtomicU64::new(0);

/// The summary line (printed by `main` at exit), only in shadow mode and after a subject took the
/// path.
// No NCBI counterpart: prints a counter; it does not change any value NCBI computes.
pub(crate) fn print_shadow_stats() {
    let subjects = SHADOW_SUBJECTS.load(Ordering::Relaxed);
    if state().mode == 2 && subjects > 0 {
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
    /// `None` when the switch is off or the pair is outside the path: a batch that
    /// [`query_fits`] did not pass (`fits`), a word size other than 3, or a subject buffer
    /// (`with_sentinels`, the residues between two sentinels) that `XSentSubject::new` refuses.
    // No NCBI counterpart: chooses between two LOSAT code paths with the same result.
    #[inline]
    pub(crate) fn new(
        query: &'a [u8],
        fits: bool,
        word_size: i32,
        with_sentinels: &'a [u8],
    ) -> Option<Self> {
        if !fits || word_size != WORD_SIZE {
            return None;
        }
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
        let reference = || {
            extend_two_hit_blosum62(
                &self.query[1..],
                self.subject_raw(),
                s_left_off,
                s_right_off,
                q_right_off,
                dropoff,
                WORD_SIZE as usize,
            )
        };
        if !(1..=MAX_DROPOFF).contains(&dropoff) {
            return reference();
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
            let expected = reference();
            assert_eq!(
                got, expected,
                "[X_SHADOW] BXUNGAPSENT two-hit mismatch: s_left_off={s_left_off} s_right_off={s_right_off} q_right_off={q_right_off} dropoff={dropoff}"
            );
            SHADOW_COMPARED.fetch_add(1, Ordering::Relaxed);
        }
        got
    }

    /// The subject as the seed loops take it: the residues and the last sentinel.
    #[inline]
    fn subject_raw(&self) -> &'a [u8] {
        &self.subject.with_sentinels()[1..]
    }
}
