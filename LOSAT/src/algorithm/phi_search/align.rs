//! PHI-BLAST's gapped alignments around a pattern hit: the preliminary score-only extensions to
//! the left of the two pattern starts and to the right of the two pattern ends (the pattern itself
//! is not scored), and the traceback: `ALIGN_EX` to the left, the pattern region aligned position
//! by position (equal wildcard runs as substitutions, unequal ones by NCBI's banded aligner) and
//! `ALIGN_EX` to the right. The extensions use the query's matrix rows (`QueryPssm`).

use super::pattern::{
    refusal_because, PatternKind, PhiError, PhiPattern, MASK_AA_ALPHABET_BITS,
    PHI_BITS_PACKED_PER_WORD, PHI_MAX_HIT,
};
use super::scan::{find_short, right_one_bits, words_left_shift, PatternHit};
use crate::algorithm::pssm_search::gapalign::{
    phi_align_ex, phi_edit_script, phi_score_only, PssmGapScratch,
};
use crate::algorithm::pssm_search::rows::{PssmRows, QueryPssm};
use crate::common::GapEditOp;
use crate::core::blast_util::BLASTAA_SIZE;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:500-501
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:246-250
/// The zero bytes after the last letter of NCBI's query buffer: the query's sentinel and the byte
/// that the buffer's allocation adds after it. A long subject's merged HSP whose pattern region
/// passes the query's end reads them in its pattern alignment (`align_patterns`); past them NCBI
/// reads outside its buffer.
pub(crate) const QUERY_ZEROS_PAST_END: usize = 2;

/// The ends and score of an alignment (`BlastGapAlignStruct`'s): starts inclusive, stops as NCBI
/// sets them for each stage.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PhiAlignment {
    pub(crate) query_start: i32,
    pub(crate) query_stop: i32,
    pub(crate) subject_start: i32,
    pub(crate) subject_stop: i32,
    pub(crate) score: i32,
    pub(crate) edit_script: Vec<GapEditOp>,
}

/// The gap costs and X-dropoff of a stage.
#[derive(Debug, Clone, Copy)]
pub(crate) struct PhiGapCosts {
    pub(crate) gap_open: i32,
    pub(crate) gap_extend: i32,
    pub(crate) x_drop: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_gapalign.c:671-737
/// `s_PHIGappedAlignment`: the left extension is skipped when either pattern start is 0. `subject`
/// is the subject's letters (then its sentinel).
#[allow(clippy::too_many_arguments)]
pub(crate) fn prelim_alignment(
    rows: &QueryPssm,
    subject: &[u8],
    subject_length: usize,
    q_off: usize,
    s_off: usize,
    q_pat_length: usize,
    s_pat_length: usize,
    costs: PhiGapCosts,
    scratch: &mut PssmGapScratch,
) -> PhiAlignment {
    let query_length = rows.len();
    let mut score_left = 0;
    let (mut query_start, mut subject_start) = (q_off as i32, s_off as i32);
    if q_off != 0 && s_off != 0 {
        let (a, b, score) = phi_score_only(
            rows,
            q_off,
            q_off,
            subject,
            0,
            s_off,
            costs.gap_open,
            costs.gap_extend,
            costs.x_drop,
            true,
            scratch,
        );
        score_left = score;
        query_start = q_off as i32 - a as i32 + 1;
        subject_start = s_off as i32 - b as i32 + 1;
    }
    let q_last = q_off + q_pat_length - 1;
    let s_last = s_off + s_pat_length - 1;
    let mut score_right = 0;
    let (mut query_stop, mut subject_stop) =
        ((q_off + q_pat_length) as i32, (s_off + s_pat_length) as i32);
    if q_last < query_length && s_last < subject_length {
        let (a, b, score) = phi_score_only(
            rows,
            query_length - q_last - 1,
            q_last,
            subject,
            s_last,
            subject_length - s_last - 1,
            costs.gap_open,
            costs.gap_extend,
            costs.x_drop,
            false,
            scratch,
        );
        score_right = score;
        query_stop = (a + q_last) as i32;
        subject_stop = (b + s_last) as i32;
    }
    PhiAlignment {
        query_start,
        query_stop,
        subject_start,
        subject_stop,
        score: score_left + score_right,
        edit_script: Vec::new(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/gapinfo.c:175-185
/// `GapPrelimEditBlockAdd`: nothing for a count of 0, else merged into the last operation of its
/// kind.
fn add_op(ops: &mut Vec<GapEditOp>, op: GapEditOp) {
    let count = match op {
        GapEditOp::Sub(n) | GapEditOp::Del(n) | GapEditOp::Ins(n) => n,
    };
    if count == 0 {
        return;
    }
    match (ops.last_mut(), op) {
        (Some(GapEditOp::Sub(last)), GapEditOp::Sub(n))
        | (Some(GapEditOp::Del(last)), GapEditOp::Del(n))
        | (Some(GapEditOp::Ins(last)), GapEditOp::Ins(n)) => *last += n,
        _ => ops.push(op),
    }
}

/// `s_Align`'s and `s_BandedAlign`'s view of the two segments: `seq1[i]` (from 1) is the query
/// letter at `query_base + i - 1`, scored through its matrix row; `seq2[j]` the subject letter at
/// `subject_base + j - 1`.
struct Segments<'a> {
    rows: &'a QueryPssm,
    /// The scores of the zero byte after the query's letters (`QUERY_ZEROS_PAST_END`).
    gap_row: &'a [i32; BLASTAA_SIZE],
    subject: &'a [u8],
    query_base: usize,
    subject_base: usize,
    gap_open: i32,
    gap_extend: i32,
}

impl Segments<'_> {
    fn score(&self, i: i32, j: i32) -> i32 {
        let letter = self.subject[self.subject_base + j as usize - 1];
        let row = self
            .rows
            .rows
            .get(self.query_base + i as usize - 1)
            .unwrap_or(self.gap_row);
        row[letter as usize]
    }
}

const DIAGONAL_INSERT: i8 = 1;
const DIAGONAL_DELETE: i8 = 2;
const INSERT_CODE: i8 = 10;
const DELETE_CODE: i8 = 20;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_gapalign.c:85-251
/// `s_Align`: the banded global alignment of `seq1[1..=end1]` and `seq2[1..=end2]` between the
/// diagonals `low_diag` and `high_diag`, its operations added to `ops`; NCBI's caller drops the
/// score.
#[allow(clippy::too_many_arguments)]
fn banded_dp(
    seg: &Segments<'_>,
    end1: i32,
    end2: i32,
    low_diag: i32,
    high_diag: i32,
    gap_open: i32,
    gap_extend: i32,
    ops: &mut Vec<GapEditOp>,
) {
    #[derive(Clone, Copy, Default)]
    struct Cell {
        best: i32,
        best_gap: i32,
    }
    let min_score = i32::MIN / 2;
    let gap_cost = gap_open + gap_extend;
    let band = (high_diag - low_diag + 1) as usize;
    let width = band + 2;
    let mut cells = vec![Cell::default(); width];
    let mut state = vec![0i8; (end1 as usize + 1) * width];
    let at = |row: i32, d: i32| row as usize * width + d as usize;
    let mut leftd = 1 - low_diag;
    let mut rightd = high_diag - low_diag + 1;
    cells[leftd as usize].best = 0;
    state[at(0, leftd)] = -1;
    let mut initial = -gap_open;
    for d in leftd + 1..=rightd {
        initial -= gap_extend;
        cells[d as usize].best = initial;
        cells[d as usize - 1].best_gap = initial - gap_cost;
        state[at(0, d)] = DIAGONAL_INSERT;
    }
    cells[rightd as usize + 1].best = min_score;
    cells[rightd as usize].best_gap = min_score;
    cells[leftd as usize - 1].best_gap = -gap_cost;
    cells[leftd as usize - 1].best = min_score;
    let mut temp_sub = 0i32;
    for i in 1..=end1 {
        if i > end2 - high_diag {
            rightd -= 1;
        }
        if leftd > 1 {
            leftd -= 1;
        }
        let mut temp_indel = cells[leftd as usize].best_gap;
        let mut next_state = 0i8;
        let index2 = leftd + low_diag - 1 + i;
        if index2 > 0 {
            temp_sub = cells[leftd as usize].best + seg.score(i, index2);
        }
        if temp_indel > temp_sub || index2 <= 0 {
            temp_sub = temp_indel;
            next_state = DIAGONAL_DELETE;
        }
        let mut temp_hor = temp_sub - gap_cost;
        if leftd >= 1 {
            temp_indel -= gap_extend;
            if temp_indel >= temp_hor {
                cells[leftd as usize - 1].best_gap = temp_indel;
                next_state += DELETE_CODE;
            } else {
                cells[leftd as usize - 1].best_gap = temp_hor;
            }
        }
        state[at(i, leftd)] = next_state;
        cells[leftd as usize].best = temp_sub;
        for curd in leftd + 1..=rightd {
            let d = curd as usize;
            temp_sub = cells[d].best + seg.score(i, curd + low_diag - 1 + i);
            temp_indel = cells[d].best_gap;
            let cell_state;
            if temp_indel > temp_sub {
                if temp_indel > temp_hor {
                    cells[d].best = temp_indel;
                    temp_hor -= gap_extend;
                    cells[d - 1].best_gap = temp_indel - gap_extend;
                    cell_state = DELETE_CODE + INSERT_CODE + DIAGONAL_DELETE;
                } else {
                    cells[d].best = temp_hor;
                    temp_hor -= gap_extend;
                    cells[d - 1].best_gap = temp_indel - gap_extend;
                    cell_state = DELETE_CODE + INSERT_CODE + DIAGONAL_INSERT;
                }
            } else if temp_hor > temp_sub {
                cells[d].best = temp_hor;
                temp_hor -= gap_extend;
                cells[d - 1].best_gap = temp_indel - gap_extend;
                cell_state = DELETE_CODE + INSERT_CODE + DIAGONAL_INSERT;
            } else {
                cells[d].best = temp_sub;
                temp_sub -= gap_cost;
                temp_hor -= gap_extend;
                if temp_sub > temp_hor {
                    temp_hor = temp_sub;
                    next_state = 0;
                } else {
                    next_state = INSERT_CODE;
                }
                temp_indel -= gap_extend;
                if temp_sub > temp_indel {
                    cell_state = next_state;
                    cells[d - 1].best_gap = temp_sub;
                } else {
                    cell_state = next_state + DELETE_CODE;
                    cells[d - 1].best_gap = temp_indel;
                }
            }
            state[at(i, curd)] = cell_state;
        }
    }
    let mut instructions = Vec::new();
    let mut index1 = end1;
    let mut diag = rightd;
    let mut edit_step = 0i8;
    while index1 >= 0 {
        let decoder = state[at(index1, diag)];
        if decoder == -1 {
            break;
        }
        let mut next_step = decoder % INSERT_CODE;
        if edit_step == DIAGONAL_INSERT && (decoder / INSERT_CODE) % 2 == 1 {
            next_step = DIAGONAL_INSERT;
        }
        if edit_step == DIAGONAL_DELETE && decoder / DELETE_CODE == 1 {
            next_step = DIAGONAL_DELETE;
        }
        if next_step == DIAGONAL_INSERT {
            diag -= 1;
            index1 += 1;
        } else if next_step == DIAGONAL_DELETE {
            diag += 1;
        }
        edit_step = next_step;
        instructions.push(edit_step);
        index1 -= 1;
    }
    for &step in instructions.iter().rev() {
        match step {
            0 => add_op(ops, GapEditOp::Sub(1)),
            DIAGONAL_INSERT => add_op(ops, GapEditOp::Del(1)),
            DIAGONAL_DELETE => add_op(ops, GapEditOp::Ins(1)),
            _ => {}
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_gapalign.c:261-313
/// `s_BandedAlign`: a run of `start1` query and `start2` subject letters inside the pattern.
fn banded_align(seg: &Segments<'_>, start1: i32, start2: i32, ops: &mut Vec<GapEditOp>) {
    const BAND_LOW: i32 = -5;
    const BAND_HIGH: i32 = 5;
    let (gap_open, gap_extend) = (seg.gap_open, seg.gap_extend);
    let low_diag = BAND_LOW.max(-start1).min((start2 - start1).min(0));
    let high_diag = start2.min(BAND_HIGH).max((start2 - start1).max(0));
    if start2 <= 0 {
        if start1 > 0 {
            add_op(ops, GapEditOp::Ins(start1 as u32));
        }
        return;
    }
    if start1 <= 0 {
        add_op(ops, GapEditOp::Del(start2 as u32));
        return;
    }
    if high_diag <= low_diag {
        for _ in 1..=start1 {
            add_op(ops, GapEditOp::Sub(1));
        }
        return;
    }
    banded_dp(
        seg, start1, start2, low_diag, high_diag, gap_open, gap_extend, ops,
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_gapalign.c:325-355
/// `s_PHIGetShortPattern`: the expanded positions (first, last) of the alternative that `seq`
/// (the pattern region) matches.
fn short_pattern_ends(pattern: &PhiPattern, seq: &[u8]) -> (i32, i32) {
    let mask = pattern.match_mask;
    let mask_shift_plus_1 = (mask << 1) + 1;
    let mut matched = 0i32;
    for &letter in seq {
        matched = ((matched << 1) | mask_shift_plus_1) & pattern.which_positions[letter as usize];
    }
    let (right_one, right_mask_only) = right_one_bits(matched, mask);
    (right_mask_only + 1, right_one)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_gapalign.c:367-417
/// `s_PHIGetLongPattern`: as `short_pattern_ends` for a multi-word pattern.
fn long_pattern_ends(pattern: &PhiPattern, seq: &[u8]) -> (i32, i32) {
    let words = pattern.num_words;
    let match_mask = &pattern.match_mask_l[..words];
    let mut mask = match_mask.to_vec();
    words_left_shift(&mut mask, 1);
    let mut matched = vec![0i32; words];
    for &letter in seq {
        words_left_shift(&mut matched, 0);
        let letter_words = &pattern.bit_pattern_by_letter[letter as usize];
        for (w, word) in matched.iter_mut().enumerate() {
            *word = (*word | mask[w]) & letter_words[w];
        }
    }
    for (word, &mask_word) in matched.iter_mut().zip(match_mask) {
        *word &= mask_word;
    }
    let mut right_mask_only = -1;
    let mut word = 0usize;
    let mut bit = 0i32;
    'words: while word < words {
        bit = 0;
        while bit < PHI_BITS_PACKED_PER_WORD {
            if (matched[word] >> bit) % 2 == 1 {
                break 'words;
            } else if (match_mask[word] >> bit) % 2 == 1 {
                right_mask_only = word as i32 * PHI_BITS_PACKED_PER_WORD + bit;
            }
            bit += 1;
        }
        word += 1;
    }
    (
        right_mask_only + 1,
        word as i32 * PHI_BITS_PACKED_PER_WORD + bit,
    )
}

/// The slots of NCBI's array for one word's hits in `s_PHIGetExtraLongPattern`.
const MAX_HITS_IN_WORD: usize = 64;

fn traceback_array_refusal(what: &str) -> PhiError {
    refusal_because(what, "writes past its traceback arrays")
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_gapalign.c:430-483
/// `s_PHIGetExtraLongPattern`: the end offsets of each word of a very long pattern in `seq` (the
/// pattern region), the first decomposition that ends at the region's end moved to the front;
/// NCBI's array of 20,000 slots (the rest as the search left it).
fn extra_long_offsets(pattern: &PhiPattern, seq: &[u8]) -> Result<Vec<i32>, PhiError> {
    let len = seq.len() as i32;
    let words = pattern.num_words;
    let mut hit_array = vec![0i32; PHI_MAX_HIT];
    let mut next = vec![0i32; PHI_MAX_HIT];
    hit_array[0] = pattern.num_places_in_word[0];
    let mut tuples_end = 1usize;
    let mut pos = 0usize;
    for word in 1..words {
        let which = &pattern.sll[word];
        let mask = pattern.match_mask_l[word];
        pos = 0;
        let mut j = 0usize;
        while j < tuples_end {
            let last = hit_array[j + word - 1];
            let window = (len - last)
                .min(pattern.spacing[word - 1] + pattern.num_places_in_word[word])
                .max(0);
            let mut found: Vec<PatternHit> = Vec::new();
            find_short(
                which,
                mask,
                &seq[last as usize..(last + window) as usize],
                &mut found,
            );
            if found.len() * 2 > MAX_HITS_IN_WORD {
                return Err(traceback_array_refusal(
                    "more than 32 hits of one word of a very long pattern in an alignment",
                ));
            }
            for hit in &found {
                if pos + word + 1 > PHI_MAX_HIT {
                    return Err(traceback_array_refusal(
                        "more than 20,000 word positions of a very long pattern in an alignment",
                    ));
                }
                next[pos..pos + word].copy_from_slice(&hit_array[j..j + word]);
                next[pos + word] = next[pos + word - 1] + hit.end + 1;
                pos += word + 1;
            }
            j += word;
        }
        tuples_end = pos;
        hit_array[..pos].copy_from_slice(&next[..pos]);
    }
    let mut j = 0usize;
    while j < pos {
        if hit_array[j + words - 1] == len {
            hit_array.copy_within(j..j + words, 0);
            break;
        }
        j += words;
    }
    Ok(hit_array)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_gapalign.c:497-656
/// `s_PHIBlastAlignPatterns`: the operations of the pattern region `query[q_start..][..q_len]`
/// against `subject[s_start..][..s_len]`. `query` and `subject` are the letters as NCBI's buffers
/// hold them, with the bytes NCBI reads after them (`traceback`); `gap_row` scores those after
/// the query's letters.
#[allow(clippy::too_many_arguments)]
fn align_patterns(
    rows: &QueryPssm,
    gap_row: &[i32; BLASTAA_SIZE],
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    s_start: usize,
    q_len: usize,
    s_len: usize,
    pattern: &PhiPattern,
    gap_open: i32,
    gap_extend: i32,
) -> Result<Vec<GapEditOp>, PhiError> {
    let mut ops = Vec::new();
    // A long subject's merged HSP aligns the pattern region from one occurrence's start with the
    // length of another (`phi_hooks::long_subject_prelim`), which can pass the last letter of
    // the query or subject: NCBI reads the zero bytes after the letters in its buffers, and past
    // them outside the buffer or the next sequence of a database.
    let (Some(query_region), Some(subject_region)) = (
        query.get(q_start..q_start + q_len),
        subject.get(s_start..s_start + s_len),
    ) else {
        return Err(refusal_because(
            "a pattern region past the bytes after the query or subject",
            "reads past the sequence's buffer",
        ));
    };
    let segments = |query_base: usize, subject_base: usize| Segments {
        rows,
        gap_row,
        subject,
        query_base,
        subject_base,
        gap_open,
        gap_extend,
    };
    let (mut q_cur, mut s_cur) = (q_start, s_start);
    let masked = &pattern.masked;
    let (start_q, end_q, start_d, end_d) = match pattern.kind {
        PatternKind::OneWord => {
            let (sq, eq) = short_pattern_ends(pattern, query_region);
            let (sd, ed) = short_pattern_ends(pattern, subject_region);
            (sq, eq, sd, ed)
        }
        PatternKind::MultiWord => {
            let (sq, eq) = long_pattern_ends(pattern, query_region);
            let (sd, ed) = long_pattern_ends(pattern, subject_region);
            (sq, eq, sd, ed)
        }
        PatternKind::VeryLong => {
            let query_words = extra_long_offsets(pattern, query_region)?;
            let subject_words = extra_long_offsets(pattern, subject_region)?;
            let (mut q_matched, mut s_matched) = (0i32, 0i32);
            let (mut q_word, mut s_word) = (0usize, 0usize);
            let (mut q_var, mut s_var) = (0i32, 0i32);
            for &entry in &masked[..pattern.highest_place] {
                if entry < 0 {
                    // The leading variable region is shifted out, so a word comes first.
                    if q_word == 0 || s_word == 0 {
                        return Err(traceback_array_refusal(
                            "a variable wildcard before the first word",
                        ));
                    }
                    q_var += query_words[q_word]
                        - query_words[q_word - 1]
                        - pattern.num_places_in_word[q_word];
                    s_var += subject_words[s_word]
                        - subject_words[s_word - 1]
                        - pattern.num_places_in_word[s_word];
                } else if entry == MASK_AA_ALPHABET_BITS {
                    q_var += 1;
                    s_var += 1;
                } else {
                    if q_var != 0 || s_var != 0 {
                        if q_var == s_var {
                            add_op(&mut ops, GapEditOp::Sub(q_var as u32));
                        } else {
                            banded_align(&segments(q_cur, s_cur), q_var, s_var, &mut ops);
                        }
                        q_matched += q_var;
                        q_cur += q_var as usize;
                        s_matched += s_var;
                        s_cur += s_var as usize;
                    }
                    add_op(&mut ops, GapEditOp::Sub(1));
                    q_matched += 1;
                    q_cur += 1;
                    s_matched += 1;
                    s_cur += 1;
                    q_var = 0;
                    s_var = 0;
                }
                if q_matched + q_var >= query_words[q_word] {
                    q_word += 1;
                }
                if s_matched + s_var >= subject_words[s_word] {
                    s_word += 1;
                }
            }
            return Ok(ops);
        }
    };
    let (mut pos_q, mut pos_d) = (start_q, start_d);
    while pos_q <= end_q || pos_d <= end_d {
        if masked[pos_q as usize] != MASK_AA_ALPHABET_BITS
            && masked[pos_d as usize] != MASK_AA_ALPHABET_BITS
        {
            add_op(&mut ops, GapEditOp::Sub(1));
            pos_q += 1;
            pos_d += 1;
            q_cur += 1;
            s_cur += 1;
        } else {
            let mut q_run = 0i32;
            while masked[pos_q as usize] == MASK_AA_ALPHABET_BITS && pos_q <= end_q {
                pos_q += 1;
                q_run += 1;
            }
            let mut s_run = 0i32;
            while masked[pos_d as usize] == MASK_AA_ALPHABET_BITS && pos_d <= end_d {
                pos_d += 1;
                s_run += 1;
            }
            if q_run == s_run {
                if q_run == 0 {
                    // NCBI adds a substitution and repeats this step without end.
                    return Err(refusal_because(
                        "a pattern region whose wildcards do not pair up",
                        "does not end",
                    ));
                }
                for _ in 0..q_run {
                    add_op(&mut ops, GapEditOp::Sub(1));
                }
                q_cur += q_run as usize;
                s_cur += q_run as usize;
            } else {
                banded_align(&segments(q_cur, s_cur), q_run, s_run, &mut ops);
                q_cur += q_run as usize;
                s_cur += s_run as usize;
            }
        }
    }
    Ok(ops)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_gapalign.c:837-906
/// `PHIGappedAlignmentWithTraceback`: `ALIGN_EX` left of the pattern starts (also at 0), the
/// pattern region, `ALIGN_EX` right of the pattern ends; the score is the two extensions'.
/// `query` is the query's letters followed by the zero bytes after them in NCBI's buffer
/// (`QUERY_ZEROS_PAST_END`), `gap_row` their scores; `subject` the subject's letters followed by
/// its sentinel, both read by a pattern region that passes their ends (`align_patterns`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn traceback(
    rows: &QueryPssm,
    gap_row: &[i32; BLASTAA_SIZE],
    query: &[u8],
    subject: &[u8],
    subject_length: usize,
    q_start: usize,
    s_start: usize,
    q_pat_length: usize,
    s_pat_length: usize,
    pattern: &PhiPattern,
    costs: PhiGapCosts,
    scratch: &mut PssmGapScratch,
) -> Result<PhiAlignment, PhiError> {
    let query_length = rows.len();
    let (a, b, score_left, mut rev) = phi_align_ex(
        rows,
        q_start,
        q_start,
        subject,
        0,
        s_start,
        costs.gap_open,
        costs.gap_extend,
        costs.x_drop,
        true,
        scratch,
    );
    let query_start = q_start as i32 - a as i32;
    let subject_start = s_start as i32 - b as i32;
    for op in align_patterns(
        rows,
        gap_row,
        query,
        subject,
        q_start,
        s_start,
        q_pat_length,
        s_pat_length,
        pattern,
        costs.gap_open,
        costs.gap_extend,
    )? {
        add_op(&mut rev, op);
    }
    let q_end = q_start + q_pat_length - 1;
    let s_end = s_start + s_pat_length - 1;
    let (query_stop, subject_stop, score_right, fwd) =
        if q_end < query_length && s_end < subject_length {
            let (a, b, score, fwd) = phi_align_ex(
                rows,
                query_length - q_end - 1,
                q_end,
                subject,
                s_end,
                subject_length - s_end - 1,
                costs.gap_open,
                costs.gap_extend,
                costs.x_drop,
                false,
                scratch,
            );
            ((q_end + a + 1) as i32, (s_end + b + 1) as i32, score, fwd)
        } else {
            (q_end as i32, s_end as i32, 0, Vec::new())
        };
    Ok(PhiAlignment {
        query_start,
        query_stop,
        subject_start,
        subject_stop,
        score: score_left + score_right,
        edit_script: phi_edit_script(rev, fwd),
    })
}

#[cfg(test)]
#[path = "align_tests.rs"]
mod tests;
