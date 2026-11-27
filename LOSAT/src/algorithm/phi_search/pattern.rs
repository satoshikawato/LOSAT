//! PHI-BLAST's pattern compiler (`SPHIPatternSearchBlkNew`): the pattern's positions as residue
//! masks, its probability and minimum match length, the expansion of the variable wildcards into
//! their alternatives and the three representations the scanners read (one word, several words,
//! very long). The arrays are NCBI's, with NCBI's sizes, so that the positions past the written
//! part keep what NCBI's keep (the stale "last position" flags of the expansion, the unexpanded
//! entries of a very long pattern). A write that NCBI makes past one of its arrays, and an input
//! whose result NCBI leaves undefined, is LOSAT's refusal (`PhiError::Refusal`), raised
//! here, before the prolog; NCBI's own messages are `PhiError::Ncbi`, which the search
//! reports where NCBI's engine raises them.

use crate::core::blast_util::BLASTAA_SIZE;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/pattern.h:50-58
/// The positions of a pattern (`PHI_MAX_PATTERN_LENGTH`).
pub(crate) const PHI_MAX_PATTERN_LENGTH: usize = 330;
/// The bits of one word of the multi-word representations (`PHI_BITS_PACKED_PER_WORD`).
pub(crate) const PHI_BITS_PACKED_PER_WORD: i32 = 30;
/// The words of a multi-word pattern's arrays (`PHI_MAX_WORD_SIZE`).
const PHI_MAX_WORD_SIZE: usize = 11;
/// The words of a very long pattern's arrays (`PHI_BUF_SIZE`, `PHI_MAX_WORDS_IN_PATTERN`).
pub(crate) const PHI_MAX_WORDS_IN_PATTERN: usize = 100;
/// The slots of NCBI's hit arrays (`PHI_MAX_HIT`): two per hit.
pub(crate) const PHI_MAX_HIT: usize = 20000;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_lookup.c:41
/// Every residue of the alphabet (`kMaskAaAlphabetBits`).
pub(crate) const MASK_AA_ALPHABET_BITS: i32 = (1 << BLASTAA_SIZE) - 1;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:105-113
/// NCBI's `AMINOACID_TO_NCBISTDAA` (128 entries): what the compiler indexes with a pattern byte.
pub(crate) const AMINOACID_TO_NCBISTDAA: [u8; 128] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 27, 10, 11, 12, 13, 26, 14, 15, 16, 17, 18, 24, 19, 20, 21,
    22, 23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
];

/// The three representations (`EPatternType`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PatternKind {
    /// The expanded pattern fits in one word of 30 bits (`eOneWord`).
    OneWord,
    /// The expanded pattern fits in 330 positions (`eMultiWord`).
    MultiWord,
    /// The words of the unexpanded pattern, scanned from the most specific (`eVeryLong`).
    VeryLong,
}

/// A compiled pattern (`SPHIPatternSearchBlk` with its one-word, multi-word and very-long items).
#[derive(Debug, Clone)]
pub(crate) struct PhiPattern {
    pub(crate) kind: PatternKind,
    /// `patternProbability`.
    pub(crate) probability: f64,
    /// `minPatternMatchLength`, in 64 bits: NCBI's `Int4` sum of the repeat counts overflows for
    /// counts near the limits of `int` (signed overflow, widened: Owner 2026-10-10, 10-21).
    pub(crate) min_match_length: i64,
    /// The pattern as NCBI saves it (`pattern_blk->pattern`): the raw text, upper-cased.
    pub(crate) text: Vec<u8>,
    /// One word: `match_mask` and `whichPositionPtr` (for each residue, the positions it matches).
    pub(crate) match_mask: i32,
    pub(crate) which_positions: [i32; BLASTAA_SIZE],
    /// Multi-word and very long: `numWords`, `match_maskL`, `bitPatternByLetter` and `SLL`.
    pub(crate) num_words: usize,
    pub(crate) match_mask_l: [i32; PHI_MAX_WORDS_IN_PATTERN],
    pub(crate) bit_pattern_by_letter: [[i32; PHI_MAX_WORD_SIZE]; BLASTAA_SIZE],
    pub(crate) sll: Vec<[i32; BLASTAA_SIZE]>,
    /// `inputPatternMasked`: the expanded positions (one word, multi-word) or the unexpanded ones
    /// (very long), then what the expansion left past them.
    pub(crate) masked: [i32; PHI_MAX_PATTERN_LENGTH],
    /// Very long: `numPlacesInWord`, `spacing`, `highestPlace`, `whichMostSpecific`.
    pub(crate) num_places_in_word: [i32; PHI_MAX_WORDS_IN_PATTERN],
    pub(crate) spacing: [i32; PHI_MAX_WORDS_IN_PATTERN],
    pub(crate) highest_place: usize,
    pub(crate) most_specific: usize,
}

/// Why a pattern does not compile or a PHI search stops.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum PhiError {
    /// NCBI's message (`Blast_MessageWrite` of `SPHIPatternSearchBlkNew`), which NCBI's engine
    /// raises after the prolog: `Warning: <message> ` with exit 3.
    Ncbi(String),
    /// LOSAT's refusal of a pattern where NCBI reads or writes past its arrays or reads an
    /// undefined value (rule 2; DESIGN_ps8 §3.1).
    Refusal(String),
}

/// LOSAT's refusal of a pattern with `what`, where NCBI BLAST+ 2.17.0 `why`.
pub(crate) fn refusal_because(what: &str, why: &str) -> PhiError {
    PhiError::Refusal(format!(
        "a -phi_pattern with {what} (NCBI BLAST+ 2.17.0 {why}) is not supported by LOSAT's PSIBLAST"
    ))
}

/// The refusal of `what` by the reason NCBI's compiler has for it.
fn refusal(what: &str) -> PhiError {
    let why = match what {
        "x(n,m) without two numbers" => "reads uninitialised values",
        "a count outside 32-bit integers" => "overflows a 32-bit integer",
        "a non-ASCII byte in the pattern" => "reads past its residue table",
        "a pattern longer than NCBI's 330 positions"
        | "330 wildcard positions"
        | "more than 100 words in a very long pattern" => "writes or reads past its arrays",
        "two variable wildcards in a row in a very long pattern" => "shifts by a negative count",
        _ => "crashes",
    };
    refusal_because(what, why)
}

/// The pattern's bytes as NCBI's parser reads them: `text[i]` for `i` below the length, NUL at
/// the length (`calloc` of one more byte), and past it a read beyond the string (refused).
struct PatternBytes<'a>(&'a [u8]);

impl PatternBytes<'_> {
    fn at(&self, index: usize, what: &str) -> Result<u8, PhiError> {
        match index.cmp(&self.0.len()) {
            std::cmp::Ordering::Less => Ok(self.0[index]),
            std::cmp::Ordering::Equal => Ok(0),
            std::cmp::Ordering::Greater => Err(refusal(what)),
        }
    }
}

/// `isspace` of the C locale.
fn c_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// A decimal integer as C's `%d` and `atoi` read it from `at`: white space, a sign, digits. `None`
/// without a digit; LOSAT's refusal for a value outside `int` (undefined in C). Returns the value
/// and the index after it.
fn read_int(bytes: &PatternBytes<'_>, at: usize) -> Result<Option<(i32, usize)>, PhiError> {
    let mut index = at;
    let mut byte = bytes.at(index, "an unterminated count")?;
    while c_space(byte) {
        index += 1;
        byte = bytes.at(index, "an unterminated count")?;
    }
    let negative = byte == b'-';
    if byte == b'-' || byte == b'+' {
        index += 1;
        byte = bytes.at(index, "an unterminated count")?;
    }
    if !byte.is_ascii_digit() {
        return Ok(None);
    }
    let mut value: i64 = 0;
    while byte.is_ascii_digit() {
        value = value * 10 + i64::from(byte - b'0');
        if value > i64::from(i32::MAX) + 1 {
            return Err(refusal("a count outside 32-bit integers"));
        }
        index += 1;
        byte = bytes.at(index, "an unterminated count")?;
    }
    let value = if negative { -value } else { value };
    match i32::try_from(value) {
        Ok(value) => Ok(Some((value, index))),
        Err(_) => Err(refusal("a count outside 32-bit integers")),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_lookup.c:141-202
/// `s_ExpandPattern` on the arrays from `base`: each variable entry replaced by its alternatives
/// (the one-wildcard form in place, then 0, then 2 to the entry's count), the copies' flags taken
/// from the snapshot and the inserted wildcards' flags left as they were. -1 when the expanded
/// pattern reaches `max_length`.
fn expand(
    masked: &mut [i32; PHI_MAX_PATTERN_LENGTH],
    flags: &mut [u8; PHI_MAX_PATTERN_LENGTH],
    base: usize,
    length: usize,
    max_length: usize,
) -> i32 {
    for i in 0..length {
        let this_place = -i64::from(masked[base + i]);
        if this_place > 0 {
            masked[base + i] = MASK_AA_ALPHABET_BITS;
            let temp_mask: Vec<i32> = masked[base..base + length].to_vec();
            let temp_flags: Vec<u8> = flags[base..base + length].to_vec();
            let first = expand(masked, flags, base, length, max_length);
            if first == -1 {
                return -1;
            }
            let mut rec1 = first as usize;
            let mut rec2 = first as usize;
            for num_pos in 0..=this_place as usize {
                if num_pos == 1 {
                    continue;
                }
                for (k, (&mask, &flag)) in temp_mask.iter().zip(&temp_flags).enumerate() {
                    if k == i {
                        for _ in 0..num_pos {
                            masked[base + rec1] = MASK_AA_ALPHABET_BITS;
                            rec1 += 1;
                            if rec1 >= max_length {
                                return -1;
                            }
                        }
                    } else {
                        masked[base + rec1] = mask;
                        flags[base + rec1] = flag;
                        rec1 += 1;
                        if rec1 >= max_length {
                            return -1;
                        }
                    }
                    if rec1 >= max_length {
                        return -1;
                    }
                }
                let added = expand(
                    masked,
                    flags,
                    base + rec2,
                    length + num_pos - 1,
                    max_length - rec2,
                );
                if added == -1 {
                    return -1;
                }
                rec2 += added as usize;
                rec1 = rec2;
            }
            return rec1 as i32;
        }
    }
    length as i32
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_lookup.c:268-278
/// `s_NumOfOne`: the set bits of a positive mask (0 for a negative one).
fn num_of_one(mut a: i32) -> i32 {
    let mut count = 0;
    while a > 0 {
        if a % 2 == 1 {
            count += 1;
        }
        a >>= 1;
    }
    count
}

impl PhiPattern {
    /// `s_PatternSearchItemsInit`'s blank items.
    fn blank(text: Vec<u8>) -> Self {
        Self {
            kind: PatternKind::OneWord,
            probability: 1.0,
            min_match_length: 0,
            text,
            match_mask: 0,
            which_positions: [0; BLASTAA_SIZE],
            num_words: 0,
            match_mask_l: [0; PHI_MAX_WORDS_IN_PATTERN],
            bit_pattern_by_letter: [[0; PHI_MAX_WORD_SIZE]; BLASTAA_SIZE],
            sll: vec![[0; BLASTAA_SIZE]; PHI_MAX_WORDS_IN_PATTERN],
            masked: [0; PHI_MAX_PATTERN_LENGTH],
            num_places_in_word: [0; PHI_MAX_WORDS_IN_PATTERN],
            spacing: [0; PHI_MAX_WORDS_IN_PATTERN],
            highest_place: 0,
            most_specific: 0,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_lookup.c:388-688
    /// `SPHIPatternSearchBlkNew` of a protein pattern (`raw`: the text as read, without the NUL
    /// that ends NCBI's string) with the standard residue probabilities (`Blast_ResFreqStdComp`).
    pub(crate) fn compile(raw: &[u8]) -> Result<Self, PhiError> {
        if raw.len() >= PHI_MAX_PATTERN_LENGTH {
            return Err(PhiError::Ncbi(format!(
                "Pattern is too long ({} but only {} supported)",
                raw.len(),
                PHI_MAX_PATTERN_LENGTH
            )));
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_lookup.c:370-385
        // `s_MakePatternUpperCase`: `a` to `z` only.
        let text: Vec<u8> = raw.iter().map(u8::to_ascii_uppercase).collect();
        let prob = crate::stats::karlin_calc::compute_std_aa_composition();
        let mut blk = Self::blank(text.clone());
        let mut flags = [0u8; PHI_MAX_PATTERN_LENGTH];
        let pattern = PatternBytes(&text);
        let length = text.len();
        let mut pos: usize = 0;
        // A write of the parser at `pos`: past the array is LOSAT's refusal (NCBI's repeat loop and
        // the `continue` after `x(n,m)` write without a check).
        fn write(
            masked: &mut [i32; PHI_MAX_PATTERN_LENGTH],
            pos: &mut usize,
            value: i32,
        ) -> Result<(), PhiError> {
            if *pos >= PHI_MAX_PATTERN_LENGTH {
                return Err(refusal("a pattern longer than NCBI's 330 positions"));
            }
            masked[*pos] = value;
            *pos += 1;
            Ok(())
        }
        let order = |byte: u8, what: &str| -> Result<usize, PhiError> {
            AMINOACID_TO_NCBISTDAA
                .get(byte as usize)
                .map(|&code| code as usize)
                .ok_or_else(|| refusal(what))
        };
        let mut index = 0usize;
        while index < length {
            let next = text[index];
            if next == 0 || next == b'\r' || next == b'\n' {
                break;
            }
            if matches!(next, b'-' | b'.' | b'>' | b' ' | b'<') {
                index += 1;
                continue;
            }
            let char_set_mask: i32;
            let position_probability: f64;
            if next != b'[' && next != b'{' {
                if next == b'X' {
                    if pattern.at(index + 1, "an unterminated x(")? == b'(' {
                        index += 1;
                        let mut second = index;
                        loop {
                            let byte = pattern.at(second, "an unterminated x(")?;
                            if byte == b',' || byte == b')' {
                                break;
                            }
                            if second >= length {
                                return Err(refusal("an unterminated x("));
                            }
                            second += 1;
                        }
                        if pattern.at(second, "an unterminated x(")? == b')' {
                            index -= 1;
                            char_set_mask = MASK_AA_ALPHABET_BITS;
                            position_probability = 1.0;
                        } else {
                            index += 1;
                            // `sscanf` of `%d,%d`: both values or LOSAT's refusal (NCBI reads
                            // uninitialised variables).
                            let Some((min_wildcard, after_first)) = read_int(&pattern, index)?
                            else {
                                return Err(refusal("x(n,m) without two numbers"));
                            };
                            if pattern.at(after_first, "x(n,m) without two numbers")? != b',' {
                                return Err(refusal("x(n,m) without two numbers"));
                            }
                            let Some((max_wildcard, _)) = read_int(&pattern, after_first + 1)?
                            else {
                                return Err(refusal("x(n,m) without two numbers"));
                            };
                            let Some(variable) = max_wildcard.checked_sub(min_wildcard) else {
                                return Err(refusal("a count outside 32-bit integers"));
                            };
                            if variable == i32::MIN {
                                return Err(refusal("a count outside 32-bit integers"));
                            }
                            blk.min_match_length += i64::from(min_wildcard);
                            let mut remaining = min_wildcard;
                            while remaining > 0 {
                                remaining -= 1;
                                write(&mut blk.masked, &mut pos, MASK_AA_ALPHABET_BITS)?;
                                if pos >= PHI_MAX_PATTERN_LENGTH {
                                    return Err(PhiError::Ncbi("Pattern too long".to_string()));
                                }
                            }
                            if variable != 0 {
                                write(&mut blk.masked, &mut pos, -variable)?;
                                blk.probability *= f64::from(variable);
                            }
                            // The skip to `)` starts after the first byte of the numbers.
                            loop {
                                index += 1;
                                if pattern.at(index, "an unterminated (")? == b')' {
                                    break;
                                }
                                if index >= length {
                                    return Err(refusal("an unterminated ("));
                                }
                            }
                            index += 1;
                            continue;
                        }
                    } else {
                        char_set_mask = MASK_AA_ALPHABET_BITS;
                        position_probability = 1.0;
                    }
                } else if next == b'U' {
                    char_set_mask = MASK_AA_ALPHABET_BITS * 2 + 1;
                    position_probability = 1.0;
                } else {
                    let code = order(next, "a non-ASCII byte in the pattern")?;
                    char_set_mask = 1 << code;
                    position_probability = prob[code];
                }
            } else if next == b'[' {
                let mut mask = 0i32;
                let mut probability = 0.0;
                loop {
                    index += 1;
                    let byte = pattern.at(index, "an unterminated [")?;
                    if byte == b']' {
                        break;
                    }
                    // A byte of 0x80 or more is a negative `char`: below `A`.
                    if !byte.is_ascii_uppercase() {
                        return Err(PhiError::Ncbi(
                            "pattern description has a non-alphabeticcharacter inside a bracket"
                                .to_string(),
                        ));
                    }
                    let code = order(byte, "a non-ASCII byte in the pattern")?;
                    mask |= 1 << code;
                    probability += prob[code];
                }
                char_set_mask = mask;
                position_probability = probability;
            } else {
                let mut mask = MASK_AA_ALPHABET_BITS;
                let mut probability = 1.0;
                loop {
                    index += 1;
                    let byte = pattern.at(index, "an unterminated {")?;
                    if byte == b'}' {
                        break;
                    }
                    if index >= length {
                        return Err(refusal("an unterminated {"));
                    }
                    let code = order(byte, "a non-ASCII byte in the pattern")?;
                    mask -= mask & (1 << code);
                    probability -= prob[code];
                }
                char_set_mask = mask;
                position_probability = probability;
            }
            if pattern.at(index + 1, "an unterminated (")? == b'(' {
                index += 1;
                index += 1;
                let num_identical = read_int(&pattern, index)?.map_or(0, |(value, _)| value);
                blk.min_match_length += i64::from(num_identical);
                loop {
                    index += 1;
                    if pattern.at(index, "an unterminated (")? == b')' {
                        break;
                    }
                    if index >= length {
                        return Err(refusal("an unterminated ("));
                    }
                }
                let mut remaining = num_identical;
                while remaining > 0 {
                    remaining -= 1;
                    write(&mut blk.masked, &mut pos, char_set_mask)?;
                    blk.probability *= position_probability;
                }
            } else {
                write(&mut blk.masked, &mut pos, char_set_mask)?;
                blk.min_match_length += 1;
                blk.probability *= position_probability;
            }
            if pos >= PHI_MAX_PATTERN_LENGTH {
                return Err(PhiError::Ncbi("Pattern is too long".to_string()));
            }
            index += 1;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_lookup.c:613-635
        // Trailing variable entries are dropped (the entry before the first, NCBI's read at -1, is
        // 0 in NCBI's zeroed structure); a leading variable region after the first non-wildcard
        // entry's search is shifted out.
        while pos > 0 && blk.masked[pos - 1] < 0 {
            pos -= 1;
        }
        let mut first = 0usize;
        while first < pos && blk.masked[first] == MASK_AA_ALPHABET_BITS {
            first += 1;
        }
        // 330 wildcard positions leave NCBI reading the entry past its array.
        let Some(&first_entry) = blk.masked.get(first) else {
            return Err(refusal("330 wildcard positions"));
        };
        if first_entry < 0 {
            let mut second = first + 1;
            while second < pos && blk.masked[second] <= 0 {
                second += 1;
            }
            while second < pos {
                blk.masked[first] = blk.masked[second];
                second += 1;
                first += 1;
            }
            pos = first;
        }
        // A pattern of no position marks nothing (NCBI writes its flag before its array).
        if pos > 0 {
            flags[pos - 1] = 1;
        }
        if blk.probability > 1.0 {
            blk.probability = 1.0;
        }
        let temp_masked = blk.masked;
        let temp_pos = pos;
        let expanded = expand(&mut blk.masked, &mut flags, 0, pos, PHI_MAX_PATTERN_LENGTH);
        if expanded == -1 {
            blk.kind = PatternKind::VeryLong;
            blk.pack_very_long(&temp_masked, temp_pos)?;
            blk.masked[..temp_pos].copy_from_slice(&temp_masked[..temp_pos]);
            blk.highest_place = temp_pos;
            return Ok(blk);
        }
        let expanded = expanded as usize;
        if expanded > PHI_BITS_PACKED_PER_WORD as usize {
            blk.kind = PatternKind::MultiWord;
            blk.pack_long(expanded, &flags);
            return Ok(blk);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_lookup.c:211-220
        // `s_PackPattern`: the flags of the expanded positions as a bit mask.
        blk.match_mask = flags[..expanded]
            .iter()
            .enumerate()
            .filter(|(_, &flag)| flag != 0)
            .fold(0, |mask, (i, _)| mask + (1 << i));
        for (residue, positions) in blk.which_positions.iter_mut().enumerate() {
            let mut this_mask = 0i32;
            for (place, &mask) in blk.masked[..expanded].iter().enumerate() {
                if (1 << residue) & mask != 0 {
                    this_mask |= 1 << place;
                }
            }
            *positions = this_mask;
        }
        Ok(blk)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_lookup.c:231-262
    /// `s_PackLongPattern`: the words of the flags and each residue's positions per word.
    fn pack_long(&mut self, num_places: usize, flags: &[u8; PHI_MAX_PATTERN_LENGTH]) {
        let bits = PHI_BITS_PACKED_PER_WORD as usize;
        self.num_words = (num_places - 1) / bits + 1;
        for word in 0..self.num_words {
            let mut bit_pattern = 0i32;
            for i in 0..bits {
                if flags[word * bits + i] != 0 {
                    bit_pattern += 1 << i;
                }
            }
            self.match_mask_l[word] = bit_pattern;
        }
        for residue in 0..BLASTAA_SIZE {
            for word in 0..self.num_words {
                let mut bit_pattern = 0i32;
                for i in 0..bits {
                    if (1 << residue) & self.masked[word * bits + i] != 0 {
                        bit_pattern |= 1 << i;
                    }
                }
                self.bit_pattern_by_letter[residue][word] = bit_pattern;
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/phi_lookup.c:286-344
    /// `s_PackVeryLongPattern`: the unexpanded positions cut into words at each variable entry and
    /// every 30 positions, each word's residue positions, the spacing after it and the most
    /// specific word (the lowest product of the positions' residue fractions). A word of no
    /// position (two variable entries in a row: NCBI shifts by -1) and a 101st word (past NCBI's
    /// arrays) are LOSAT's refusals.
    fn pack_very_long(
        &mut self,
        masked: &[i32; PHI_MAX_PATTERN_LENGTH],
        num_places: usize,
    ) -> Result<(), PhiError> {
        let bits = PHI_BITS_PACKED_PER_WORD as usize;
        let mut most_specific = 1.0f64;
        self.most_specific = 0;
        let mut word_probability = 1.0f64;
        let mut word = 0usize;
        let mut place_in_word: i64 = 0;
        let mut place: i64 = 0;
        while place <= num_places as i64 {
            let p = place as usize;
            if p == num_places || masked[p] < 0 || place_in_word == bits as i64 {
                if place_in_word == 0 {
                    return Err(refusal(
                        "two variable wildcards in a row in a very long pattern",
                    ));
                }
                if word >= PHI_MAX_WORDS_IN_PATTERN {
                    return Err(refusal("more than 100 words in a very long pattern"));
                }
                let in_word = place_in_word as usize;
                self.match_mask_l[word] = 1 << (in_word - 1);
                for residue in 0..BLASTAA_SIZE {
                    let mut one_word_mask = 0i32;
                    for place_in_word2 in 0..in_word {
                        if (1 << residue) & masked[p - in_word + place_in_word2] != 0 {
                            one_word_mask |= 1 << place_in_word2;
                        }
                    }
                    self.sll[word][residue] = one_word_mask;
                }
                self.num_places_in_word[word] = in_word as i32;
                if word_probability < most_specific {
                    most_specific = word_probability;
                    self.most_specific = word;
                }
                if p == num_places {
                    self.spacing[word] = 0;
                } else if masked[p] < 0 {
                    self.spacing[word] = -masked[p];
                } else {
                    place -= 1;
                    self.spacing[word] = 0;
                }
                word += 1;
                place_in_word = -1;
                word_probability = 1.0;
            } else {
                word_probability *= f64::from(num_of_one(masked[p])) / BLASTAA_SIZE as f64;
            }
            place += 1;
            place_in_word += 1;
        }
        self.num_words = word;
        Ok(())
    }
}

#[cfg(test)]
#[path = "pattern_tests.rs"]
mod tests;
