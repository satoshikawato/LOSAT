//! The protein ungapped extension (`s_BlastAaExtendOneHit`, `s_BlastAaExtendTwoHit` and the
//! `s_BlastAaExtendLeft` / `s_BlastAaExtendRight` walks they share) with a tighter inner loop:
//! both extensions with BLOSUM62 (`one_hit`, `two_hit`), the two-hit one with another matrix
//! (`two_hit_of`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:804-883
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
//!   byte; the table is used only if every other entry is within -127..=127. A residue code >= 32
//!   would index a part of the 8 KiB table that holds zeros, so each extension ORs together every
//!   code it reads and is handed to the reference when the OR is >= 32.
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
//! Another matrix (`matrix_table`, `two_hit_of`). The same walks with a table of the matrix's
//! scores (`protein_score`, the function the reference calls for a matrix other than BLOSUM62),
//! when every score other than `BLAST_SCORE_MIN` is within -127..=127. The arguments above hold
//! with the matrix's largest score `M` in place of BLOSUM62's 11: the bounds, the codes and the
//! walks' sentinel steps do not depend on the scores, and in a word of `w` residues with
//! `w * M <= 127` the scores after a sentinel add up to at most `(w - 1) * M < 128`. In the two-hit
//! word loop a score after a sentinel is then at most the score before it (itself at most
//! `left_score`), with -128 as with -32768, so `left_score` is not raised in either. A longer
//! word (`MatrixTable::max_word_size`) runs the reference. (There is no one-hit extension with
//! another matrix's table.)

use super::extension::{
    BlastpOneHitUngappedResult, BlastpTwoHitUngappedResult, BlastpUngappedData,
};
use crate::config::ScoringMatrix;
use crate::utils::matrix::{blosum62_score_ncbistdaa_direct, protein_score, BLAST_SCORE_MIN};
use std::sync::OnceLock;

/// Codes below this index the table; `q << 5 | s` is the index.
const CODES: usize = 32;
/// `q` and `s` are bytes, so every `q << 5 | s` is below this: a code >= 32 reads zeros.
const TABLE_LEN: usize = 255 * CODES + 256;
/// Largest X-drop the byte table is exact for (see the module comment).
const MAX_DROPOFF: i32 = 127;
/// Longest word whose sentinel sum stays below zero (11 * 11 < 128).
const MAX_WORD_SIZE: usize = 11;

pub(crate) type Table = [i8; TABLE_LEN];

/// The table of the walks, `None` if a score does not fit a byte (`x_ungapped_sent.rs` scores
/// with it too).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:804-883 `matrix[q[i]][s[i]]`; the scores are those of the
// function the reference walks call.
pub(crate) fn build_table() -> Option<Box<Table>> {
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

/// The table of the walks, built once; `None` if a score does not fit a byte (the walks then hand
/// every call to the reference).
fn table() -> Option<&'static Table> {
    static TABLE: OnceLock<Option<Box<Table>>> = OnceLock::new();
    TABLE.get_or_init(build_table).as_deref()
}

/// Sum that the scores of a word may reach after a sentinel (module comment, another matrix).
const WORD_SCORE_ROOM: i32 = 127;

/// Another matrix's table of the walks (`two_hit_of`) and the longest word the walks take with it.
pub(crate) struct MatrixTable {
    table: Box<Table>,
    /// `WORD_SCORE_ROOM / M`, `M` the matrix's largest score (module comment).
    max_word_size: usize,
}

/// The table of the walks for `matrix`, `None` if a score other than `BLAST_SCORE_MIN` does not
/// fit a byte.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:804-883 `matrix[q[i]][s[i]]`; the scores are those of the
// function the reference walks call for a matrix other than BLOSUM62 (`protein_score`).
fn build_matrix_table(matrix: ScoringMatrix) -> Option<MatrixTable> {
    let mut table = Box::new([0i8; TABLE_LEN]);
    let mut max_score = 1i32;
    for q in 0..CODES {
        for s in 0..CODES {
            let v = protein_score(matrix, q as u8, s as u8);
            table[q << 5 | s] = if v == BLAST_SCORE_MIN {
                -128
            } else if (-127..=127).contains(&v) {
                max_score = max_score.max(v);
                v as i8
            } else {
                return None;
            };
        }
    }
    Some(MatrixTable {
        table,
        max_word_size: (WORD_SCORE_ROOM / max_score) as usize,
    })
}

impl MatrixTable {
    // No NCBI counterpart: chooses between two LOSAT code paths with the same result.
    /// The table for the sentinel walks of `x_ungapped_sent.rs` (item SC2c): `None` unless the
    /// sentinel (code 0) scores -128 against every code below 32, both ways, and the matrix's
    /// largest score is at most 31 (`max_word_size` at least 4).
    pub(crate) fn sent_table(&self) -> Option<&Table> {
        let sentinel = (0..CODES).all(|c| self.table[c] == -128 && self.table[c << 5] == -128);
        (sentinel && self.max_word_size >= 4).then_some(&*self.table)
    }
}

/// The table of the walks for `matrix`, built at the first call for the matrix; `None` if a score
/// does not fit a byte (the reference then runs).
pub(crate) fn matrix_table(matrix: ScoringMatrix) -> Option<&'static MatrixTable> {
    static TABLES: [OnceLock<Option<MatrixTable>>; ScoringMatrix::ALL.len()] =
        [const { OnceLock::new() }; ScoringMatrix::ALL.len()];
    TABLES[matrix as usize]
        .get_or_init(|| build_matrix_table(matrix))
        .as_ref()
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:804-883 (see the module comment).
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

/// `s_BlastAaExtendLeft` at `s_off` and `q_off`: (score, length, codes read).
#[inline(always)]
fn extend_left(
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:853-883
        // `best_i` is one more than `n`, with no step taken.
        return (maxscore, 0, 0);
    }
    // n + 1 steps down from (q_off, s_off); both offsets are >= n.
    assert!((q_off as usize) < query.len() && (s_off as usize) < subject.len());
    // SAFETY: the assert and n <= min(q_off, s_off) keep the n + 1 positions inside.
    let w = unsafe {
        walk::<false>(
            tab,
            query.as_ptr().add(q_off as usize),
            subject.as_ptr().add(s_off as usize),
            n as usize + 1,
            dropoff,
            maxscore,
        )
    };
    (w.best, w.len, w.codes)
}

/// `s_BlastAaExtendRight` at `s_off` and `q_off`: (score, length, s_last_off, codes read).
#[inline(always)]
fn extend_right(
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:804-835, the loop does not run: i = 0.
        return (maxscore, 0, s_off, 0);
    }
    // n steps up from (q_off, s_off); n is the room left in both slices. A negative offset was
    // an index panic in the reference.
    assert!(q_off >= 0 && s_off >= 0);
    // SAFETY: the n positions are inside both slices.
    let w = unsafe {
        walk::<true>(
            tab,
            query.as_ptr().add(q_off as usize),
            subject.as_ptr().add(s_off as usize),
            n as usize,
            dropoff,
            maxscore,
        )
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

/// The steps of [`two_hit`] and [`two_hit_of`] after their checks, on the table `$tab`: one text
/// for both, expanded in place, so that BLOSUM62's `two_hit` compiles as it did before `two_hit_of`
/// (an `#[inline(always)]` function in its place changed `two_hit`'s register allocation).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1155 `s_BlastAaExtendTwoHit`.
macro_rules! two_hit_steps {
    ($tab:ident, $query:ident, $subject:ident, $s_left_off:ident, $s_right_off:ident,
     $q_right_off:ident, $dropoff:ident, $word_size:ident) => {{
        let mut q_right_off = $q_right_off as i32;
        let s_left_off: i32 = $s_left_off;
        let mut s_right_off = $s_right_off as i32;
        let word_size = $word_size as i32;

        let mut right_d = 0i32;
        let mut left_score = 0i32;
        let mut score = 0i32;
        let mut codes = 0u8;
        for i in 0..word_size {
            // the word fits in both sequences (checked above)
            let qc = $query[(q_right_off + i) as usize];
            let sc = $subject[(s_right_off + i) as usize];
            codes |= qc | sc;
            score += $tab[(qc as usize) << 5 | sc as usize] as i32;
            if score > left_score {
                left_score = score;
                right_d = i + 1;
            }
        }
        q_right_off += right_d;
        s_right_off += right_d;

        let (left_score, left_d, c) = extend_left(
            $tab,
            $query,
            $subject,
            q_right_off - 1,
            s_right_off - 1,
            $dropoff,
            0,
        );
        codes |= c;

        let mut right_extend = false;
        let mut s_last_off = s_right_off;
        let mut right_score = 0i32;
        let mut right_d = 0i32;
        if left_d >= (s_right_off - s_left_off) {
            right_extend = true;
            let (rs, rd, last, c) = extend_right(
                $tab,
                $query,
                $subject,
                q_right_off,
                s_right_off,
                $dropoff,
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
    }};
}

/// `Some(result)` when the fast path handled the call (`result` may be `None`, the reference's
/// early return), `None` when the reference has to run.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1155 `s_BlastAaExtendTwoHit`.
pub(crate) fn two_hit(
    query: &[u8],
    subject: &[u8],
    s_left_off: i32,
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
    let tab = table()?;
    two_hit_steps!(
        tab,
        query,
        subject,
        s_left_off,
        s_right_off,
        q_right_off,
        dropoff,
        word_size
    )
}

/// [`two_hit`] with another matrix's table ([`matrix_table`]), for words of up to
/// `matrix.max_word_size` residues.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1089-1155 `s_BlastAaExtendTwoHit`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn two_hit_of(
    matrix: &MatrixTable,
    query: &[u8],
    subject: &[u8],
    s_left_off: i32,
    s_right_off: usize,
    q_right_off: usize,
    dropoff: i32,
    word_size: usize,
) -> Option<Option<BlastpTwoHitUngappedResult>> {
    // The early return of the reference wrapper (a word that does not fit).
    if q_right_off + word_size > query.len() || s_right_off + word_size > subject.len() {
        return Some(None);
    }
    if !((1..=MAX_DROPOFF).contains(&dropoff)
        && (1..=matrix.max_word_size).contains(&word_size)
        && query.len() <= i32::MAX as usize
        && subject.len() <= i32::MAX as usize)
    {
        return None;
    }
    let tab: &Table = &matrix.table;
    two_hit_steps!(
        tab,
        query,
        subject,
        s_left_off,
        s_right_off,
        q_right_off,
        dropoff,
        word_size
    )
}

/// As [`two_hit`], for `s_BlastAaExtendOneHit`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1020-1083 `s_BlastAaExtendOneHit`.
pub(crate) fn one_hit(
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
    let tab = table()?;
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

    let (left_score, left_disp, c) = extend_left(
        tab,
        query,
        subject,
        q_left_off - 1,
        s_left_off - 1,
        dropoff,
        score,
    );
    codes |= c;
    let (total_score, right_disp, s_last_off, c) = extend_right(
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
#[path = "x_ungapped_fast_tests.rs"]
mod tests;
