//! Discontiguous megablast words: the templates, the word index of each template, and the
//! subject scans that read discontiguous words (NCBI `blast_nalookup.h`,
//! `blast_nalookup.c` `s_FillDiscMBTable`, `blast_nascan.c` `s_MB_DiscWordScanSubject_*`).
//!
//! The lookup table itself is the megablast table (`lookup.rs` `TwoStageLookup`), built
//! with these indices (`lookup.rs` `build_disc_mb_lookup`).

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:184-189
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DiscWordType {
    Coding = 0,
    Optimal = 1,
    TwoTemplates = 2,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:191-235
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DiscTemplateType {
    Contiguous = 0,
    T11_16Coding = 1,
    T11_16Optimal = 2,
    T12_16Coding = 3,
    T12_16Optimal = 4,
    T11_18Coding = 5,
    T11_18Optimal = 6,
    T12_18Coding = 7,
    T12_18Optimal = 8,
    T11_21Coding = 9,
    T11_21Optimal = 10,
    T12_21Coding = 11,
    T12_21Optimal = 12,
}

impl DiscTemplateType {
    /// The template whose enumerated value is `value` (`(EDiscTemplateType) temp_int`).
    pub(crate) const fn from_value(value: i32) -> Self {
        match value {
            1 => Self::T11_16Coding,
            2 => Self::T11_16Optimal,
            3 => Self::T12_16Coding,
            4 => Self::T12_16Optimal,
            5 => Self::T11_18Coding,
            6 => Self::T11_18Optimal,
            7 => Self::T12_18Coding,
            8 => Self::T12_18Optimal,
            9 => Self::T11_21Coding,
            10 => Self::T11_21Optimal,
            11 => Self::T12_21Coding,
            12 => Self::T12_21Optimal,
            _ => Self::Contiguous,
        }
    }

    /// The second template of a two-template search: the next enumerated value.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:707-715
    pub fn next(self) -> Self {
        Self::from_value(self as i32 + 1)
    }
}

/// The template of a word size (weight), template length and template type; every
/// unsupported combination is the contiguous template.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:600-643
pub fn get_disc_template_type(
    weight: i32,
    length: u8,
    word_type: DiscWordType,
) -> DiscTemplateType {
    let coding = matches!(word_type, DiscWordType::Coding | DiscWordType::TwoTemplates);
    let optimal = word_type == DiscWordType::Optimal;
    let pick = |c: DiscTemplateType, o: DiscTemplateType| {
        if coding {
            c
        } else if optimal {
            o
        } else {
            DiscTemplateType::Contiguous
        }
    };
    match (weight, length) {
        (11, 16) => pick(
            DiscTemplateType::T11_16Coding,
            DiscTemplateType::T11_16Optimal,
        ),
        (11, 18) => pick(
            DiscTemplateType::T11_18Coding,
            DiscTemplateType::T11_18Optimal,
        ),
        (11, 21) => pick(
            DiscTemplateType::T11_21Coding,
            DiscTemplateType::T11_21Optimal,
        ),
        (12, 16) => pick(
            DiscTemplateType::T12_16Coding,
            DiscTemplateType::T12_16Optimal,
        ),
        (12, 18) => pick(
            DiscTemplateType::T12_18Coding,
            DiscTemplateType::T12_18Optimal,
        ),
        (12, 21) => pick(
            DiscTemplateType::T12_21Coding,
            DiscTemplateType::T12_21Optimal,
        ),
        _ => DiscTemplateType::Contiguous,
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:300-316
#[inline(always)]
fn discontig_index_11_16_coding(accum: u64) -> u32 {
    let lo = accum as u32;
    (lo & 0x00000003)
        | ((lo & 0x000000f0) >> 2)
        | ((lo & 0x00003c00) >> 4)
        | ((lo & 0x000f0000) >> 6)
        | ((lo & 0x03c00000) >> 8)
        | ((lo & 0xf0000000) >> 10)
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:318-333
#[inline(always)]
fn discontig_index_11_16_optimal(accum: u64) -> u32 {
    let lo = accum as u32;
    (lo & 0x0000003f)
        | ((lo & 0x00000f00) >> 2)
        | ((lo & 0x0003c000) >> 4)
        | ((lo & 0x00300000) >> 6)
        | ((lo & 0xfc000000) >> 10)
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:335-353
#[inline(always)]
fn discontig_index_11_18_coding(accum: u64) -> u32 {
    let lo = accum as u32;
    let hi = (accum >> 32) as u32;
    (lo & 0x00000003)
        | ((lo & 0x000000f0) >> 2)
        | ((lo & 0x00003c00) >> 4)
        | ((lo & 0x00030000) >> 6)
        | ((lo & 0x03c00000) >> 10)
        | ((lo & 0xf0000000) >> 12)
        | ((hi & 0x0000000c) << 18)
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:355-373
#[inline(always)]
fn discontig_index_11_18_optimal(accum: u64) -> u32 {
    let lo = accum as u32;
    let hi = (accum >> 32) as u32;
    (lo & 0x0000003f)
        | ((lo & 0x00000300) >> 2)
        | ((lo & 0x0003c000) >> 6)
        | ((lo & 0x00300000) >> 8)
        | ((lo & 0x0c000000) >> 12)
        | ((lo & 0xc0000000) >> 14)
        | ((hi & 0x0000000f) << 18)
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:375-394
#[inline(always)]
fn discontig_index_11_21_coding(accum: u64) -> u32 {
    let lo = accum as u32;
    let hi = (accum >> 32) as u32;
    (lo & 0x00000003)
        | ((lo & 0x000000f0) >> 2)
        | ((lo & 0x00000c00) >> 4)
        | ((lo & 0x000f0000) >> 8)
        | ((lo & 0x00c00000) >> 10)
        | ((lo & 0xf0000000) >> 14)
        | ((hi & 0x0000000c) << 16)
        | ((hi & 0x00000300) << 12)
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:396-414
#[inline(always)]
fn discontig_index_11_21_optimal(accum: u64) -> u32 {
    let lo = accum as u32;
    let hi = (accum >> 32) as u32;
    (lo & 0x0000003f)
        | ((lo & 0x00000300) >> 2)
        | ((lo & 0x0000c000) >> 6)
        | ((lo & 0x00c00000) >> 12)
        | ((lo & 0x0c000000) >> 14)
        | ((hi & 0x00000003) << 14)
        | ((hi & 0x000003f0) << 12)
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:416-431
#[inline(always)]
fn discontig_index_12_16_coding(accum: u64) -> u32 {
    let lo = accum as u32;
    (lo & 0x00000003)
        | ((lo & 0x000000f0) >> 2)
        | ((lo & 0x00003c00) >> 4)
        | ((lo & 0x000f0000) >> 6)
        | ((lo & 0xffc00000) >> 8)
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:433-448
#[inline(always)]
fn discontig_index_12_16_optimal(accum: u64) -> u32 {
    let lo = accum as u32;
    (lo & 0x0000003f)
        | ((lo & 0x00000f00) >> 2)
        | ((lo & 0x0003c000) >> 4)
        | ((lo & 0x00f00000) >> 6)
        | ((lo & 0xfc000000) >> 8)
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:450-468
#[inline(always)]
fn discontig_index_12_18_coding(accum: u64) -> u32 {
    let lo = accum as u32;
    let hi = (accum >> 32) as u32;
    (lo & 0x00000003)
        | ((lo & 0x000000f0) >> 2)
        | ((lo & 0x00003c00) >> 4)
        | ((lo & 0x000f0000) >> 6)
        | ((lo & 0x03c00000) >> 8)
        | ((lo & 0xf0000000) >> 10)
        | ((hi & 0x0000000c) << 20)
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:470-488
#[inline(always)]
fn discontig_index_12_18_optimal(accum: u64) -> u32 {
    let lo = accum as u32;
    let hi = (accum >> 32) as u32;
    (lo & 0x0000003f)
        | ((lo & 0x00000f00) >> 2)
        | ((lo & 0x0000c000) >> 4)
        | ((lo & 0x00f00000) >> 8)
        | ((lo & 0x0c000000) >> 10)
        | ((lo & 0xc0000000) >> 12)
        | ((hi & 0x0000000f) << 20)
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:490-509
#[inline(always)]
fn discontig_index_12_21_coding(accum: u64) -> u32 {
    let lo = accum as u32;
    let hi = (accum >> 32) as u32;
    (lo & 0x00000003)
        | ((lo & 0x000000f0) >> 2)
        | ((lo & 0x00000c00) >> 4)
        | ((lo & 0x000f0000) >> 8)
        | ((lo & 0x03c00000) >> 10)
        | ((lo & 0xf0000000) >> 12)
        | ((hi & 0x0000000c) << 18)
        | ((hi & 0x00000300) << 14)
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:511-529
#[inline(always)]
fn discontig_index_12_21_optimal(accum: u64) -> u32 {
    let lo = accum as u32;
    let hi = (accum >> 32) as u32;
    (lo & 0x0000003f)
        | ((lo & 0x00000300) >> 2)
        | ((lo & 0x0000c000) >> 6)
        | ((lo & 0x00f00000) >> 10)
        | ((lo & 0x0c000000) >> 12)
        | ((hi & 0x00000003) << 16)
        | ((hi & 0x000003f0) << 14)
}

/// The lookup table index of the discontiguous word that ends with the base in the
/// low-order bits of `accum`.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_nalookup.h:540-588
#[inline(always)]
pub fn compute_discontiguous_index(accum: u64, template_type: DiscTemplateType) -> u32 {
    match template_type {
        DiscTemplateType::T11_16Coding => discontig_index_11_16_coding(accum),
        DiscTemplateType::T12_16Coding => discontig_index_12_16_coding(accum),
        DiscTemplateType::T11_16Optimal => discontig_index_11_16_optimal(accum),
        DiscTemplateType::T12_16Optimal => discontig_index_12_16_optimal(accum),
        DiscTemplateType::T11_18Coding => discontig_index_11_18_coding(accum),
        DiscTemplateType::T12_18Coding => discontig_index_12_18_coding(accum),
        DiscTemplateType::T11_18Optimal => discontig_index_11_18_optimal(accum),
        DiscTemplateType::T12_18Optimal => discontig_index_12_18_optimal(accum),
        DiscTemplateType::T11_21Coding => discontig_index_11_21_coding(accum),
        DiscTemplateType::T12_21Coding => discontig_index_12_21_coding(accum),
        DiscTemplateType::T11_21Optimal => discontig_index_11_21_optimal(accum),
        DiscTemplateType::T12_21Optimal => discontig_index_12_21_optimal(accum),
        DiscTemplateType::Contiguous => 0,
    }
}

/// The scan routine NCBI chooses for a discontiguous table.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:2617-2627
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DiscScanSubjectKind {
    General,
    TwoTemplates,
    Scan11_18,
    Scan11_21,
}

pub fn choose_disc_scan_subject(
    template_type: DiscTemplateType,
    two_templates: bool,
) -> DiscScanSubjectKind {
    if two_templates {
        DiscScanSubjectKind::TwoTemplates
    } else if template_type == DiscTemplateType::T11_18Coding {
        DiscScanSubjectKind::Scan11_18
    } else if template_type == DiscTemplateType::T11_21Coding {
        DiscScanSubjectKind::Scan11_21
    } else {
        DiscScanSubjectKind::General
    }
}

/// Scans the packed (ncbi2na) subject from `scan_range[0]` to `scan_range[1]` (the start
/// offsets of the words) with stride 1, calling `on_word(s_off, index, index2)` for every
/// word: `index` is the word's index for the first template and `index2` for the second
/// one (two templates only). NCBI looks the words up in the scan (`MB_ACCESS_HITS`,
/// `MB_ACCESS_HITS2`) and stops when its offset array is full, to resume at the same
/// offset; the caller here looks each word up in the same order.
pub fn disc_word_scan_subject(
    kind: DiscScanSubjectKind,
    packed: &[u8],
    template_type: DiscTemplateType,
    second_template_type: DiscTemplateType,
    template_length: i32,
    scan_range: [i32; 2],
    on_word: &mut impl FnMut(i32, u32, Option<u32>),
) {
    disc_word_scan_subject_until(
        kind,
        packed,
        template_type,
        second_template_type,
        template_length,
        scan_range,
        &mut |s_off, index, index2| {
            on_word(s_off, index, index2);
            true
        },
    );
}

/// `disc_word_scan_subject` with a stop: the scan returns after the word for which `on_word`
/// returns `false`, as NCBI's scan returns when its offset array is full (`total_hits >=
/// max_hits`, blast_nascan.c:1447-1454) with the next offset in `scan_range[0]`; a scan started
/// again from any word's offset gives the words from that one on, as NCBI's next call does.
/// EXPERIMENT (LOSAT_X_SEEDPIECE): the strips paths' seed pieces (`x_pair_seed_piece.rs`).
pub(crate) fn disc_word_scan_subject_until(
    kind: DiscScanSubjectKind,
    packed: &[u8],
    template_type: DiscTemplateType,
    second_template_type: DiscTemplateType,
    template_length: i32,
    scan_range: [i32; 2],
    on_word: &mut impl FnMut(i32, u32, Option<u32>) -> bool,
) {
    match kind {
        DiscScanSubjectKind::General => {
            mb_disc_word_scan_subject_1(packed, template_type, template_length, scan_range, on_word)
        }
        DiscScanSubjectKind::TwoTemplates => mb_disc_word_scan_subject_two_templates_1(
            packed,
            template_type,
            second_template_type,
            template_length,
            scan_range,
            on_word,
        ),
        DiscScanSubjectKind::Scan11_18 => {
            mb_disc_word_scan_subject_11_18_1(packed, scan_range, on_word)
        }
        DiscScanSubjectKind::Scan11_21 => {
            mb_disc_word_scan_subject_11_21_1(packed, scan_range, on_word)
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:2202-2273
fn mb_disc_word_scan_subject_1(
    packed: &[u8],
    template_type: DiscTemplateType,
    template_length: i32,
    mut scan_range: [i32; 2],
    on_word: &mut impl FnMut(i32, u32, Option<u32>) -> bool,
) {
    let mut s = (scan_range[0] / COMPRESSION_RATIO_I32) as usize;
    let mut accum: u64 = 0;
    let mut index = scan_range[0] - (scan_range[0] % COMPRESSION_RATIO_I32);
    while index < scan_range[0] + template_length {
        accum = (accum << 8) | packed[s] as u64;
        s += 1;
        index += COMPRESSION_RATIO_I32;
    }
    let mut phase = match index - (scan_range[0] + template_length) {
        1 => 3,
        2 => 2,
        3 => {
            accum >>= 8;
            s -= 1;
            1
        }
        _ => 0,
    };
    let mut rounds_done = false;
    loop {
        match phase {
            0 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                if !rounds_done {
                    rounds_done = true;
                    let go = x_disc_rounds(
                        packed,
                        (&mut s, &mut accum, &mut scan_range),
                        on_word,
                        |accum, b| (accum << 8) | u64::from(b),
                        |accum, k| {
                            let shifted = accum >> X_DISC_PLACE_SHIFT[k];
                            (compute_discontiguous_index(shifted, template_type), None)
                        },
                    );
                    if !go {
                        return;
                    }
                    continue;
                }
                if !on_word(
                    scan_range[0],
                    compute_discontiguous_index(accum, template_type),
                    None,
                ) {
                    return;
                }
                scan_range[0] += 1;
                phase = 1;
            }
            1 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                accum = (accum << 8) | packed[s] as u64;
                s += 1;
                if !on_word(
                    scan_range[0],
                    compute_discontiguous_index(accum >> 6, template_type),
                    None,
                ) {
                    return;
                }
                scan_range[0] += 1;
                phase = 2;
            }
            2 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                if !on_word(
                    scan_range[0],
                    compute_discontiguous_index(accum >> 4, template_type),
                    None,
                ) {
                    return;
                }
                scan_range[0] += 1;
                phase = 3;
            }
            _ => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                if !on_word(
                    scan_range[0],
                    compute_discontiguous_index(accum >> 2, template_type),
                    None,
                ) {
                    return;
                }
                scan_range[0] += 1;
                phase = 0;
            }
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:2289-2367
fn mb_disc_word_scan_subject_two_templates_1(
    packed: &[u8],
    template_type: DiscTemplateType,
    second_template_type: DiscTemplateType,
    template_length: i32,
    mut scan_range: [i32; 2],
    on_word: &mut impl FnMut(i32, u32, Option<u32>) -> bool,
) {
    let mut s = (scan_range[0] / COMPRESSION_RATIO_I32) as usize;
    let mut accum: u64 = 0;
    let mut index = scan_range[0] - (scan_range[0] % COMPRESSION_RATIO_I32);
    while index < scan_range[0] + template_length {
        accum = (accum << 8) | packed[s] as u64;
        s += 1;
        index += COMPRESSION_RATIO_I32;
    }
    let mut phase = match index - (scan_range[0] + template_length) {
        1 => 3,
        2 => 2,
        3 => {
            accum >>= 8;
            s -= 1;
            1
        }
        _ => 0,
    };
    // (a macro, not a closure: the rounds below hand `on_word` on as well)
    macro_rules! word {
        ($s_off:expr, $shifted:expr) => {{
            let shifted: u64 = $shifted;
            on_word(
                $s_off,
                compute_discontiguous_index(shifted, template_type),
                Some(compute_discontiguous_index(shifted, second_template_type)),
            )
        }};
    }
    let mut rounds_done = false;
    loop {
        match phase {
            0 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                if !rounds_done {
                    rounds_done = true;
                    let go = x_disc_rounds(
                        packed,
                        (&mut s, &mut accum, &mut scan_range),
                        on_word,
                        |accum, b| (accum << 8) | u64::from(b),
                        |accum, k| {
                            let shifted = accum >> X_DISC_PLACE_SHIFT[k];
                            (
                                compute_discontiguous_index(shifted, template_type),
                                Some(compute_discontiguous_index(shifted, second_template_type)),
                            )
                        },
                    );
                    if !go {
                        return;
                    }
                    continue;
                }
                if !word!(scan_range[0], accum) {
                    return;
                }
                scan_range[0] += 1;
                phase = 1;
            }
            1 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                accum = (accum << 8) | packed[s] as u64;
                s += 1;
                if !word!(scan_range[0], accum >> 6) {
                    return;
                }
                scan_range[0] += 1;
                phase = 2;
            }
            2 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                if !word!(scan_range[0], accum >> 4) {
                    return;
                }
                scan_range[0] += 1;
                phase = 3;
            }
            _ => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                if !word!(scan_range[0], accum >> 2) {
                    return;
                }
                scan_range[0] += 1;
                phase = 0;
            }
        }
    }
}

/// No NCBI counterpart (bounds checks only): the whole rounds of a discontiguous scan from its
/// first place (`phase` 0): four words per round, the first from the accumulator `acc`, then the
/// next byte shifted in (`shift_in`) and the other three (`word(acc, k)`: the word or words at the
/// round's place `k`, as the scan's phases compute them). The rounds are counted once, those whose
/// four positions are at most `scan_range[1]` and whose byte is in `packed`, and read with no test
/// of the length or the range per word, as NCBI's loop body (whose subject has a sentinel byte);
/// the scan's own loop takes the words after them. Returns `false` when `on_word` stopped the scan.
#[inline(always)]
fn x_disc_rounds<A: Copy>(
    packed: &[u8],
    (s, acc, scan_range): (&mut usize, &mut A, &mut [i32; 2]),
    on_word: &mut impl FnMut(i32, u32, Option<u32>) -> bool,
    shift_in: impl Fn(A, u8) -> A,
    word: impl Fn(A, usize) -> (u32, Option<u32>),
) -> bool {
    let by_pos = match scan_range[1].checked_sub(scan_range[0] + 3) {
        Some(d) if d >= 0 => d as usize / 4 + 1,
        _ => 0,
    };
    let rounds = by_pos.min(packed.len().saturating_sub(*s));
    for _ in 0..rounds {
        let pos = scan_range[0];
        let (index, index2) = word(*acc, 0);
        if !on_word(pos, index, index2) {
            return false;
        }
        // SAFETY: `rounds` is at most `packed.len() - s` at the first round, and `s` moves one
        // byte per round, so `s < packed.len()` here.
        *acc = shift_in(*acc, unsafe { *packed.get_unchecked(*s) });
        *s += 1;
        for k in 1..4 {
            let (index, index2) = word(*acc, k);
            if !on_word(pos + k as i32, index, index2) {
                return false;
            }
        }
        scan_range[0] = pos + 4;
    }
    true
}

/// The shifts of a discontiguous scan's 64-bit accumulator at the four places of a byte
/// (`mb_disc_word_scan_subject_1`'s phases 0 to 3).
const X_DISC_PLACE_SHIFT: [u32; 4] = [0, 6, 4, 2];

/// The words of template 11 of 18 (coding) at the four places of a byte, in the scan's order:
/// the extractions of `mb_disc_word_scan_subject_11_18_1` below (NCBI's, in its reference).
#[inline(always)]
fn x_index_11_18_0(lo: u32, hi: u32) -> u32 {
    (lo & 0x00000003)
        | ((lo & 0x000000f0) >> 2)
        | ((lo & 0x00003c00) >> 4)
        | ((lo & 0x00030000) >> 6)
        | ((lo & 0x03c00000) >> 10)
        | ((lo & 0xf0000000) >> 12)
        | ((hi & 0x0000000c) << 18)
}

#[inline(always)]
fn x_index_11_18_1(lo: u32, hi: u32) -> u32 {
    ((lo & 0x000000c0) >> 6)
        | ((lo & 0x00003c00) >> 8)
        | ((lo & 0x000f0000) >> 10)
        | ((lo & 0x00c00000) >> 12)
        | ((lo & 0xf0000000) >> 16)
        | ((hi & 0x0000003c) << 14)
        | ((hi & 0x00000300) << 12)
}

#[inline(always)]
fn x_index_11_18_2(lo: u32, hi: u32) -> u32 {
    ((lo & 0x00000030) >> 4)
        | ((lo & 0x00000f00) >> 6)
        | ((lo & 0x0003c000) >> 8)
        | ((lo & 0x00300000) >> 10)
        | ((lo & 0x3c000000) >> 14)
        | ((hi & 0x0000000f) << 16)
        | ((hi & 0x000000c0) << 14)
}

#[inline(always)]
fn x_index_11_18_3(lo: u32, hi: u32) -> u32 {
    ((lo & 0x0000000c) >> 2)
        | ((lo & 0x000003c0) >> 4)
        | ((lo & 0x0000f000) >> 6)
        | ((lo & 0x000c0000) >> 8)
        | ((lo & 0x0f000000) >> 12)
        | ((lo & 0xc0000000) >> 14)
        | ((hi & 0x00000003) << 18)
        | ((hi & 0x00000030) << 16)
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:2366-2479
fn mb_disc_word_scan_subject_11_18_1(
    packed: &[u8],
    mut scan_range: [i32; 2],
    on_word: &mut impl FnMut(i32, u32, Option<u32>) -> bool,
) {
    const K_TEMPLATE_LENGTH: i32 = 18;
    let mut s = (scan_range[0] / COMPRESSION_RATIO_I32) as usize;
    let mut lo: u32 = 0;
    let mut hi: u32 = 0;
    let mut index = scan_range[0] - (scan_range[0] % COMPRESSION_RATIO_I32);
    while index < scan_range[0] + K_TEMPLATE_LENGTH {
        hi = (hi << 8) | (lo >> 24);
        lo = (lo << 8) | packed[s] as u32;
        s += 1;
        index += COMPRESSION_RATIO_I32;
    }
    let mut phase = match index - (scan_range[0] + K_TEMPLATE_LENGTH) {
        1 => 3,
        2 => 2,
        3 => {
            s -= 1;
            lo = (lo >> 8) | (hi << 24);
            hi >>= 8;
            1
        }
        _ => 0,
    };
    let mut rounds_done = false;
    loop {
        match phase {
            0 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                if !rounds_done {
                    rounds_done = true;
                    let mut acc = (lo, hi);
                    let go = x_disc_rounds(
                        packed,
                        (&mut s, &mut acc, &mut scan_range),
                        on_word,
                        |(lo, hi), b| ((lo << 8) | u32::from(b), (hi << 8) | (lo >> 24)),
                        |(lo, hi), k| {
                            let index = match k {
                                0 => x_index_11_18_0(lo, hi),
                                1 => x_index_11_18_1(lo, hi),
                                2 => x_index_11_18_2(lo, hi),
                                _ => x_index_11_18_3(lo, hi),
                            };
                            (index, None)
                        },
                    );
                    (lo, hi) = acc;
                    if !go {
                        return;
                    }
                    continue;
                }
                let index = x_index_11_18_0(lo, hi);
                if !on_word(scan_range[0], index, None) {
                    return;
                }
                scan_range[0] += 1;
                phase = 1;
            }
            1 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                hi = (hi << 8) | (lo >> 24);
                lo = (lo << 8) | packed[s] as u32;
                s += 1;
                let index = x_index_11_18_1(lo, hi);
                if !on_word(scan_range[0], index, None) {
                    return;
                }
                scan_range[0] += 1;
                phase = 2;
            }
            2 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                let index = x_index_11_18_2(lo, hi);
                if !on_word(scan_range[0], index, None) {
                    return;
                }
                scan_range[0] += 1;
                phase = 3;
            }
            _ => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                let index = x_index_11_18_3(lo, hi);
                if !on_word(scan_range[0], index, None) {
                    return;
                }
                scan_range[0] += 1;
                phase = 0;
            }
        }
    }
}

/// The words of template 11 of 21 (coding) at the four places of a byte, in the scan's order:
/// the extractions of `mb_disc_word_scan_subject_11_21_1` below (NCBI's, in its reference).
#[inline(always)]
fn x_index_11_21_0(lo: u32, hi: u32) -> u32 {
    (lo & 0x00000003)
        | ((lo & 0x000000f0) >> 2)
        | ((lo & 0x00000c00) >> 4)
        | ((lo & 0x000f0000) >> 8)
        | ((lo & 0x00c00000) >> 10)
        | ((lo & 0xf0000000) >> 14)
        | ((hi & 0x0000000c) << 16)
        | ((hi & 0x00000300) << 12)
}

#[inline(always)]
fn x_index_11_21_1(lo: u32, hi: u32) -> u32 {
    ((lo & 0x000000c0) >> 6)
        | ((lo & 0x00003c00) >> 8)
        | ((lo & 0x00030000) >> 10)
        | ((lo & 0x03c00000) >> 14)
        | ((lo & 0x30000000) >> 16)
        | ((hi & 0x0000003c) << 12)
        | ((hi & 0x00000300) << 10)
        | ((hi & 0x0000c000) << 6)
}

#[inline(always)]
fn x_index_11_21_2(lo: u32, hi: u32) -> u32 {
    ((lo & 0x00000030) >> 4)
        | ((lo & 0x00000f00) >> 6)
        | ((lo & 0x0000c000) >> 8)
        | ((lo & 0x00f00000) >> 12)
        | ((lo & 0x0c000000) >> 14)
        | ((hi & 0x0000000f) << 14)
        | ((hi & 0x000000c0) << 12)
        | ((hi & 0x00003000) << 8)
}

#[inline(always)]
fn x_index_11_21_3(lo: u32, hi: u32) -> u32 {
    ((lo & 0x0000000c) >> 2)
        | ((lo & 0x000003c0) >> 4)
        | ((lo & 0x00003000) >> 6)
        | ((lo & 0x003c0000) >> 10)
        | ((lo & 0x03000000) >> 12)
        | ((lo & 0xc0000000) >> 16)
        | ((hi & 0x00000003) << 16)
        | ((hi & 0x00000030) << 14)
        | ((hi & 0x00000c00) << 10)
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:2481-2598
fn mb_disc_word_scan_subject_11_21_1(
    packed: &[u8],
    mut scan_range: [i32; 2],
    on_word: &mut impl FnMut(i32, u32, Option<u32>) -> bool,
) {
    const K_TEMPLATE_LENGTH: i32 = 21;
    let mut s = (scan_range[0] / COMPRESSION_RATIO_I32) as usize;
    let mut lo: u32 = 0;
    let mut hi: u32 = 0;
    let mut index = scan_range[0] - (scan_range[0] % COMPRESSION_RATIO_I32);
    while index < scan_range[0] + K_TEMPLATE_LENGTH {
        hi = (hi << 8) | (lo >> 24);
        lo = (lo << 8) | packed[s] as u32;
        s += 1;
        index += COMPRESSION_RATIO_I32;
    }
    let mut phase = match index - (scan_range[0] + K_TEMPLATE_LENGTH) {
        1 => 3,
        2 => 2,
        3 => {
            s -= 1;
            lo = (lo >> 8) | (hi << 24);
            hi >>= 8;
            1
        }
        _ => 0,
    };
    let mut rounds_done = false;
    loop {
        match phase {
            0 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                if !rounds_done {
                    rounds_done = true;
                    let mut acc = (lo, hi);
                    let go = x_disc_rounds(
                        packed,
                        (&mut s, &mut acc, &mut scan_range),
                        on_word,
                        |(lo, hi), b| ((lo << 8) | u32::from(b), (hi << 8) | (lo >> 24)),
                        |(lo, hi), k| {
                            let index = match k {
                                0 => x_index_11_21_0(lo, hi),
                                1 => x_index_11_21_1(lo, hi),
                                2 => x_index_11_21_2(lo, hi),
                                _ => x_index_11_21_3(lo, hi),
                            };
                            (index, None)
                        },
                    );
                    (lo, hi) = acc;
                    if !go {
                        return;
                    }
                    continue;
                }
                let index = x_index_11_21_0(lo, hi);
                if !on_word(scan_range[0], index, None) {
                    return;
                }
                scan_range[0] += 1;
                phase = 1;
            }
            1 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                hi = (hi << 8) | (lo >> 24);
                lo = (lo << 8) | packed[s] as u32;
                s += 1;
                let index = x_index_11_21_1(lo, hi);
                if !on_word(scan_range[0], index, None) {
                    return;
                }
                scan_range[0] += 1;
                phase = 2;
            }
            2 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                let index = x_index_11_21_2(lo, hi);
                if !on_word(scan_range[0], index, None) {
                    return;
                }
                scan_range[0] += 1;
                phase = 3;
            }
            _ => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                let index = x_index_11_21_3(lo, hi);
                if !on_word(scan_range[0], index, None) {
                    return;
                }
                scan_range[0] += 1;
                phase = 0;
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:83
const COMPRESSION_RATIO_I32: i32 = 4;

#[cfg(test)]
#[path = "disc_lookup_tests.rs"]
mod tests;
