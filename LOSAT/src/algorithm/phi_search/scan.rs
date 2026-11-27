//! The pattern's hits in a sequence (`FindPatternHits`): the one-word finder (at most 10,000 hits,
//! NCBI's cap), the multi-word finder and the very long finder (the most specific word, then the
//! words to its right and to its left inside their spacing), the query's occurrences inside its
//! lookup segments and the subject's hits. A hit is NCBI's pair of slots: the end, then the start,
//! both 0-based in the scanned sequence. Where NCBI writes past one of its hit arrays, the scan
//! stops with LOSAT's refusal.

use super::pattern::{
    refusal_because, PatternKind, PhiError, PhiPattern, PHI_BITS_PACKED_PER_WORD, PHI_MAX_HIT,
};
use crate::core::blast_util::BLASTAA_SIZE;

/// A hit: the 0-based end and start of the pattern in the scanned sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct PatternHit {
    pub(crate) end: i32,
    pub(crate) start: i32,
}

fn hit_array_refusal(what: &str) -> PhiError {
    refusal_because(what, "writes past its hit array")
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/pattern.c:61-101
/// `_PHIGetRightOneBits`: the lowest bit set in both `s` and `mask` (0 when none) and the highest
/// bit of `mask` below it (-1 when none).
pub(crate) fn right_one_bits(s: i32, mask: i32) -> (i32, i32) {
    let checking = s & mask;
    let mut left = -1;
    let mut right = 0;
    while right < PHI_BITS_PACKED_PER_WORD {
        if (checking >> right) % 2 == 1 {
            break;
        }
        if (mask >> right) % 2 == 1 {
            left = right;
        }
        right += 1;
    }
    if right == PHI_BITS_PACKED_PER_WORD {
        right = 0;
    }
    (right, left)
}

/// `s_LenOf`.
fn len_of(s: i32, mask: i32) -> i32 {
    let (right, left) = right_one_bits(s, mask);
    right - left
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/pattern.c:104-141
/// `_PHIBlastFindHitsShort` of the word `which`/`mask` over `seq`, appended to `hits`: NCBI stops
/// at its 10,000th hit of the call.
pub(crate) fn find_short(
    which: &[i32; BLASTAA_SIZE],
    mask: i32,
    seq: &[u8],
    hits: &mut Vec<PatternHit>,
) {
    let mask_shift_plus_1 = (mask << 1) + 1;
    let mut matched = 0i32;
    let mut count = 0usize;
    for (i, &letter) in seq.iter().enumerate() {
        matched = ((matched << 1) | mask_shift_plus_1) & which[letter as usize];
        if matched & mask != 0 {
            let end = i as i32;
            hits.push(PatternHit {
                end,
                start: end - len_of(matched, mask) + 1,
            });
            count += 2;
            if count == PHI_MAX_HIT {
                break;
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/pattern.c:236-254
/// `_PHIPatternWordsLeftShift`: the words shifted left by one with `carry` in.
pub(crate) fn words_left_shift(words: &mut [i32], mut carry: i32) {
    let overflow = 1i32 << PHI_BITS_PACKED_PER_WORD;
    for word in words {
        let x = (*word << 1) + carry;
        if x >= overflow {
            *word = x - overflow;
            carry = 1;
        } else {
            *word = x;
            carry = 0;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/pattern.c:286-303
/// `s_LenOfL`.
fn len_of_long(s: &[i32], mask: &[i32]) -> i32 {
    let mut first_one_in_mask = -1;
    for (word, (&s_word, &mask_word)) in s.iter().zip(mask).enumerate() {
        for bit in 0..PHI_BITS_PACKED_PER_WORD {
            let place = word as i32 * PHI_BITS_PACKED_PER_WORD + bit;
            if (s_word >> bit) % 2 == 1 {
                return place - first_one_in_mask;
            }
            if (mask_word >> bit) % 2 == 1 {
                first_one_in_mask = place;
            }
        }
    }
    -1
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/pattern.c:315-356
/// `s_FindHitsLong`: no cap; past `capacity` hits NCBI writes past the caller's array.
fn find_long(
    pattern: &PhiPattern,
    seq: &[u8],
    capacity: usize,
) -> Result<Vec<PatternHit>, PhiError> {
    let words = pattern.num_words;
    let match_mask = &pattern.match_mask_l[..words];
    let mut mask = match_mask.to_vec();
    words_left_shift(&mut mask, 1);
    let mut matched = vec![0i32; words];
    let mut result = vec![0i32; words];
    let mut hits = Vec::new();
    for (i, &letter) in seq.iter().enumerate() {
        words_left_shift(&mut matched, 0);
        let letter_words = &pattern.bit_pattern_by_letter[letter as usize];
        for (w, word) in matched.iter_mut().enumerate() {
            *word = (*word | mask[w]) & letter_words[w];
        }
        let mut any = false;
        for (w, word) in result.iter_mut().enumerate() {
            *word = matched[w] & match_mask[w];
            any |= *word != 0;
        }
        if any {
            if hits.len() == capacity {
                return Err(hit_array_refusal(
                    "more than 10,000 hits of a multi-word pattern in one sequence",
                ));
            }
            let end = i as i32;
            hits.push(PatternHit {
                end,
                start: end - len_of_long(&result, match_mask) + 1,
            });
        }
    }
    Ok(hits)
}

/// The one-word items of word `word` of a very long pattern.
fn very_long_word(pattern: &PhiPattern, word: usize) -> (&[i32; BLASTAA_SIZE], i32) {
    (&pattern.sll[word], pattern.match_mask_l[word])
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/pattern.c:368-466
/// `s_FindHitsVeryLong`: the most specific word's hits, then each word to the right found after
/// the hit within the spacing before it, then each word to the left before the hit within its
/// spacing. Past 20,000 slots of NCBI's local array, or past `capacity` hits of the caller's,
/// NCBI writes out of bounds.
fn find_very_long(
    pattern: &PhiPattern,
    seq: &[u8],
    capacity: usize,
) -> Result<Vec<PatternHit>, PhiError> {
    let len = seq.len() as i32;
    let first = pattern.most_specific;
    let (which, mask) = very_long_word(pattern, first);
    let mut hits = Vec::new();
    find_short(which, mask, seq, &mut hits);
    if hits.is_empty() {
        return Ok(hits);
    }
    let local_limit = PHI_MAX_HIT / 2;
    let too_many = || {
        hit_array_refusal("more than 10,000 partial hits of a very long pattern in one sequence")
    };
    for word in first + 1..pattern.num_words {
        let (which, mask) = very_long_word(pattern, word);
        let mut next = Vec::new();
        for hit in &hits {
            let from = hit.end + 1;
            let n = (len - hit.end - 1)
                .min(pattern.spacing[word - 1] + pattern.num_places_in_word[word]);
            let before = next.len();
            find_short(
                which,
                mask,
                &seq[from as usize..(from + n) as usize],
                &mut next,
            );
            if next.len() > local_limit {
                return Err(too_many());
            }
            for found in &mut next[before..] {
                found.end += hit.end + 1;
                found.start = hit.start;
            }
        }
        if next.is_empty() {
            return Ok(next);
        }
        if next.len() > capacity {
            return Err(too_many());
        }
        hits = next;
    }
    for word in (0..first).rev() {
        let (which, mask) = very_long_word(pattern, word);
        let mut next = Vec::new();
        for hit in &hits {
            let start =
                (hit.start - pattern.spacing[word] - pattern.num_places_in_word[word]).max(0);
            let before = next.len();
            find_short(
                which,
                mask,
                &seq[start as usize..hit.start as usize],
                &mut next,
            );
            if next.len() > local_limit {
                return Err(too_many());
            }
            for found in &mut next[before..] {
                found.end = hit.end;
                found.start += start;
            }
        }
        if next.is_empty() {
            return Ok(next);
        }
        if next.len() > capacity {
            return Err(too_many());
        }
        hits = next;
    }
    Ok(hits)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/pattern.c:468-476
/// `FindPatternHits` of a protein pattern over `seq` (NCBIstdaa letters), whose caller's array
/// holds `capacity` hits.
pub(crate) fn find_pattern_hits(
    pattern: &PhiPattern,
    seq: &[u8],
    capacity: usize,
) -> Result<Vec<PatternHit>, PhiError> {
    match pattern.kind {
        PatternKind::OneWord => {
            let mut hits = Vec::new();
            find_short(&pattern.which_positions, pattern.match_mask, seq, &mut hits);
            Ok(hits)
        }
        PatternKind::MultiWord => find_long(pattern, seq, capacity),
        PatternKind::VeryLong => find_very_long(pattern, seq, capacity),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_lookup.c:725-756
/// The subject's hits (`PHIBlastScanSubject`, into NCBI's array of `PHI_MAX_HIT` slots).
pub(crate) fn subject_hits(pattern: &PhiPattern, seq: &[u8]) -> Result<Vec<PatternHit>, PhiError> {
    find_pattern_hits(pattern, seq, PHI_MAX_HIT / 2)
}

/// A pattern occurrence in the query: its 0-based offset and its length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Occurrence {
    pub(crate) offset: i32,
    pub(crate) length: i32,
}

/// The query's occurrences, or why there are none to search with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Occurrences {
    Found(Vec<Occurrence>),
    /// One hit covers the whole query (NCBI's `INT4_MAX`).
    WholeQuery,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/pattern.c:553-596
/// `PHIGetPatternOccurrences`: the hits inside each lookup segment (`segments`, inclusive 0-based
/// ranges of `query`) in order; a hit at 0 as long as the query ends the search at once.
pub(crate) fn query_occurrences(
    pattern: &PhiPattern,
    query: &[u8],
    segments: &[(usize, usize)],
) -> Result<Occurrences, PhiError> {
    let capacity = query.len();
    let mut occurrences = Vec::new();
    for &(from, to) in segments {
        let hits = find_pattern_hits(pattern, &query[from..=to], capacity)?;
        for hit in hits {
            let offset = hit.start + from as i32;
            let length = hit.end - hit.start + 1;
            if offset == 0 && length as usize == query.len() {
                return Ok(Occurrences::WholeQuery);
            }
            occurrences.push(Occurrence { offset, length });
        }
    }
    Ok(Occurrences::Found(occurrences))
}

#[cfg(test)]
#[path = "scan_tests.rs"]
mod tests;
