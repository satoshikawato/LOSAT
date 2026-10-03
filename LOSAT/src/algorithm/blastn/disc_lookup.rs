//! Discontiguous megablast words: the templates, the word index of each template, and the
//! subject scans that read discontiguous words (NCBI `blast_nalookup.h`,
//! `blast_nalookup.c` `s_FillDiscMBTable`, `blast_nascan.c` `s_MB_DiscWordScanSubject_*`).
//!
//! The lookup table itself is the megablast table (`lookup.rs` `TwoStageLookup`), built
//! with these indices (`lookup.rs` `build_disc_mb_lookup`).

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:184-189
// ```c
// /** General types of discontiguous word templates */
// typedef enum {
//    eMBWordCoding = 0,
//    eMBWordOptimal = 1,
//    eMBWordTwoTemplates = 2
// } EDiscWordType;
// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DiscWordType {
    Coding = 0,
    Optimal = 1,
    TwoTemplates = 2,
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:191-235
// ```c
// /** Enumeration of all discontiguous word templates; the enumerated values
//  * encode the weight, template length and type information
//  *
//  * <PRE>
//  *  Optimal word templates:
//  * Number of 1's in a template is word size (weight);
//  * total number of 1's and 0's - template length.
//  *   1,110,110,110,110,111      - 12 of 16
//  *   1,110,010,110,110,111      - 11 of 16
//  * 111,010,110,010,110,111      - 12 of 18
//  * 111,010,010,110,010,111      - 11 of 18
//  * 111,010,010,110,010,010,111  - 12 of 21
//  * 111,010,010,100,010,010,111  - 11 of 21
//  *  Coding word templates:
//  *    111,110,110,110,110,1     - 12 of 16
//  *    110,110,110,110,110,1     - 11 of 16
//  * 10,110,110,110,110,110,1     - 12 of 18
//  * 10,110,110,010,110,110,1     - 11 of 18
//  * 10,010,110,110,110,010,110,1 - 12 of 21
//  * 10,010,110,010,110,010,110,1 - 11 of 21
//  * </PRE>
//  ...
//  */
// typedef enum {
//    eDiscTemplateContiguous = 0,
//    eDiscTemplate_11_16_Coding = 1,
//    eDiscTemplate_11_16_Optimal = 2,
//    eDiscTemplate_12_16_Coding = 3,
//    eDiscTemplate_12_16_Optimal = 4,
//    eDiscTemplate_11_18_Coding = 5,
//    eDiscTemplate_11_18_Optimal = 6,
//    eDiscTemplate_12_18_Coding = 7,
//    eDiscTemplate_12_18_Optimal = 8,
//    eDiscTemplate_11_21_Coding = 9,
//    eDiscTemplate_11_21_Optimal = 10,
//    eDiscTemplate_12_21_Coding = 11,
//    eDiscTemplate_12_21_Optimal = 12
// } EDiscTemplateType;
// ```
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
    fn from_value(value: i32) -> Self {
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
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:707-715
    /// ```c
    ///    /* For now leave only one possibility for the second template.
    ///       Note that the intention here is to select both the coding
    ///       and the optimal templates for one combination of word size
    ///       and template length. */
    ///    if (kTwoTemplates) {
    ///       /* Use the temporaray to avoid annoying ICC warning. */
    ///       int temp_int = template_type + 1;
    ///       second_template_type =
    ///            mb_lt->second_template_type = (EDiscTemplateType) temp_int;
    /// ```
    pub fn next(self) -> Self {
        Self::from_value(self as i32 + 1)
    }
}

/// The template of a word size (weight), template length and template type; every
/// unsupported combination is the contiguous template.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:600-643
/// ```c
/// static EDiscTemplateType
/// s_GetDiscTemplateType(Int4 weight, Uint1 length,
///                       EDiscWordType type)
/// {
///    if (weight == 11) {
///       if (length == 16) {
///          if (type == eMBWordCoding || type == eMBWordTwoTemplates)
///             return eDiscTemplate_11_16_Coding;
///          else if (type == eMBWordOptimal)
///             return eDiscTemplate_11_16_Optimal;
///       } else if (length == 18) {
///          if (type == eMBWordCoding || type == eMBWordTwoTemplates)
///             return eDiscTemplate_11_18_Coding;
///          else if (type == eMBWordOptimal)
///             return eDiscTemplate_11_18_Optimal;
///       } else if (length == 21) {
///          if (type == eMBWordCoding || type == eMBWordTwoTemplates)
///             return eDiscTemplate_11_21_Coding;
///          else if (type == eMBWordOptimal)
///             return eDiscTemplate_11_21_Optimal;
///       }
///    } else if (weight == 12) {
///       if (length == 16) {
///          if (type == eMBWordCoding || type == eMBWordTwoTemplates)
///             return eDiscTemplate_12_16_Coding;
///          else if (type == eMBWordOptimal)
///             return eDiscTemplate_12_16_Optimal;
///       } else if (length == 18) {
///          if (type == eMBWordCoding || type == eMBWordTwoTemplates)
///             return eDiscTemplate_12_18_Coding;
///          else if (type == eMBWordOptimal)
///             return eDiscTemplate_12_18_Optimal;
///       } else if (length == 21) {
///          if (type == eMBWordCoding || type == eMBWordTwoTemplates)
///             return eDiscTemplate_12_21_Coding;
///          else if (type == eMBWordOptimal)
///             return eDiscTemplate_12_21_Optimal;
///       }
///    }
///    return eDiscTemplateContiguous; /* All unsupported cases default to 0 */
/// }
/// ```
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

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:300-316
// ```c
// static NCBI_INLINE Int4 DiscontigIndex_11_16_Coding(Uint8 accum)
// {
//     Uint4 lo = (Uint4)accum;
//     return ((lo & 0x00000003)      ) |
//            ((lo & 0x000000f0) >>  2) |
//            ((lo & 0x00003c00) >>  4) |
//            ((lo & 0x000f0000) >>  6) |
//            ((lo & 0x03c00000) >>  8) |
//            ((lo & 0xf0000000) >> 10);
// }
// ```
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

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:318-333
// ```c
// static NCBI_INLINE Int4 DiscontigIndex_11_16_Optimal(Uint8 accum)
// {
//     Uint4 lo = (Uint4)accum;
//     return ((lo & 0x0000003f)      ) |
//            ((lo & 0x00000f00) >>  2) |
//            ((lo & 0x0003c000) >>  4) |
//            ((lo & 0x00300000) >>  6) |
//            ((lo & 0xfc000000) >> 10);
// }
// ```
#[inline(always)]
fn discontig_index_11_16_optimal(accum: u64) -> u32 {
    let lo = accum as u32;
    (lo & 0x0000003f)
        | ((lo & 0x00000f00) >> 2)
        | ((lo & 0x0003c000) >> 4)
        | ((lo & 0x00300000) >> 6)
        | ((lo & 0xfc000000) >> 10)
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:335-353
// ```c
// static NCBI_INLINE Int4 DiscontigIndex_11_18_Coding(Uint8 accum)
// {
//     Uint4 lo = (Uint4)accum;
//     Uint4 hi = (Uint4)(accum >> 32);
//     return ((lo & 0x00000003)      ) |
//            ((lo & 0x000000f0) >>  2) |
//            ((lo & 0x00003c00) >>  4) |
//            ((lo & 0x00030000) >>  6) |
//            ((lo & 0x03c00000) >> 10) |
//            ((lo & 0xf0000000) >> 12) |
//            ((hi & 0x0000000c) << 18);
// }
// ```
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

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:355-373
// ```c
// static NCBI_INLINE Int4 DiscontigIndex_11_18_Optimal(Uint8 accum)
// {
//     Uint4 lo = (Uint4)accum;
//     Uint4 hi = (Uint4)(accum >> 32);
//     return ((lo & 0x0000003f)      ) |
//            ((lo & 0x00000300) >>  2) |
//            ((lo & 0x0003c000) >>  6) |
//            ((lo & 0x00300000) >>  8) |
//            ((lo & 0x0c000000) >> 12) |
//            ((lo & 0xc0000000) >> 14) |
//            ((hi & 0x0000000f) << 18);
// }
// ```
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

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:375-394
// ```c
// static NCBI_INLINE Int4 DiscontigIndex_11_21_Coding(Uint8 accum)
// {
//     Uint4 lo = (Uint4)accum;
//     Uint4 hi = (Uint4)(accum >> 32);
//     return ((lo & 0x00000003)      ) |
//            ((lo & 0x000000f0) >>  2) |
//            ((lo & 0x00000c00) >>  4) |
//            ((lo & 0x000f0000) >>  8) |
//            ((lo & 0x00c00000) >> 10) |
//            ((lo & 0xf0000000) >> 14) |
//            ((hi & 0x0000000c) << 16) |
//            ((hi & 0x00000300) << 12);
// }
// ```
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

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:396-414
// ```c
// static NCBI_INLINE Int4 DiscontigIndex_11_21_Optimal(Uint8 accum)
// {
//     Uint4 lo = (Uint4)accum;
//     Uint4 hi = (Uint4)(accum >> 32);
//     return ((lo & 0x0000003f)      ) |
//            ((lo & 0x00000300) >>  2) |
//            ((lo & 0x0000c000) >>  6) |
//            ((lo & 0x00c00000) >> 12) |
//            ((lo & 0x0c000000) >> 14) |
//            ((hi & 0x00000003) << 14) |
//            ((hi & 0x000003f0) << 12);
// }
// ```
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

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:416-431
// ```c
// static NCBI_INLINE Int4 DiscontigIndex_12_16_Coding(Uint8 accum)
// {
//     Uint4 lo = (Uint4)accum;
//     return ((lo & 0x00000003)      ) |
//            ((lo & 0x000000f0) >>  2) |
//            ((lo & 0x00003c00) >>  4) |
//            ((lo & 0x000f0000) >>  6) |
//            ((lo & 0xffc00000) >>  8);
// }
// ```
#[inline(always)]
fn discontig_index_12_16_coding(accum: u64) -> u32 {
    let lo = accum as u32;
    (lo & 0x00000003)
        | ((lo & 0x000000f0) >> 2)
        | ((lo & 0x00003c00) >> 4)
        | ((lo & 0x000f0000) >> 6)
        | ((lo & 0xffc00000) >> 8)
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:433-448
// ```c
// static NCBI_INLINE Int4 DiscontigIndex_12_16_Optimal(Uint8 accum)
// {
//     Uint4 lo = (Uint4)accum;
//     return ((lo & 0x0000003f)     ) |
//            ((lo & 0x00000f00) >> 2) |
//            ((lo & 0x0003c000) >> 4) |
//            ((lo & 0x00f00000) >> 6) |
//            ((lo & 0xfc000000) >> 8);
// }
// ```
#[inline(always)]
fn discontig_index_12_16_optimal(accum: u64) -> u32 {
    let lo = accum as u32;
    (lo & 0x0000003f)
        | ((lo & 0x00000f00) >> 2)
        | ((lo & 0x0003c000) >> 4)
        | ((lo & 0x00f00000) >> 6)
        | ((lo & 0xfc000000) >> 8)
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:450-468
// ```c
// static NCBI_INLINE Int4 DiscontigIndex_12_18_Coding(Uint8 accum)
// {
//     Uint4 lo = (Uint4)accum;
//     Uint4 hi = (Uint4)(accum >> 32);
//     return ((lo & 0x00000003)      ) |
//            ((lo & 0x000000f0) >>  2) |
//            ((lo & 0x00003c00) >>  4) |
//            ((lo & 0x000f0000) >>  6) |
//            ((lo & 0x03c00000) >>  8) |
//            ((lo & 0xf0000000) >> 10) |
//            ((hi & 0x0000000c) << 20);
// }
// ```
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

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:470-488
// ```c
// static NCBI_INLINE Int4 DiscontigIndex_12_18_Optimal(Uint8 accum)
// {
//     Uint4 lo = (Uint4)accum;
//     Uint4 hi = (Uint4)(accum >> 32);
//     return ((lo & 0x0000003f)      ) |
//            ((lo & 0x00000f00) >>  2) |
//            ((lo & 0x0000c000) >>  4) |
//            ((lo & 0x00f00000) >>  8) |
//            ((lo & 0x0c000000) >> 10) |
//            ((lo & 0xc0000000) >> 12) |
//            ((hi & 0x0000000f) << 20);
// }
// ```
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

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:490-509
// ```c
// static NCBI_INLINE Int4 DiscontigIndex_12_21_Coding(Uint8 accum)
// {
//     Uint4 lo = (Uint4)accum;
//     Uint4 hi = (Uint4)(accum >> 32);
//     return ((lo & 0x00000003)      ) |
//            ((lo & 0x000000f0) >>  2) |
//            ((lo & 0x00000c00) >>  4) |
//            ((lo & 0x000f0000) >>  8) |
//            ((lo & 0x03c00000) >> 10) |
//            ((lo & 0xf0000000) >> 12) |
//            ((hi & 0x0000000c) << 18) |
//            ((hi & 0x00000300) << 14);
// }
// ```
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

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:511-529
// ```c
// static NCBI_INLINE Int4 DiscontigIndex_12_21_Optimal(Uint8 accum)
// {
//     Uint4 lo = (Uint4)accum;
//     Uint4 hi = (Uint4)(accum >> 32);
//     return ((lo & 0x0000003f)      ) |
//            ((lo & 0x00000300) >>  2) |
//            ((lo & 0x0000c000) >>  6) |
//            ((lo & 0x00f00000) >> 10) |
//            ((lo & 0x0c000000) >> 12) |
//            ((hi & 0x00000003) << 16) |
//            ((hi & 0x000003f0) << 14);
// }
// ```
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
/// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_nalookup.h:540-588
/// ```c
/// static NCBI_INLINE Int4 ComputeDiscontiguousIndex(Uint8 accum,
///                                     EDiscTemplateType template_type)
/// {
///    Int4 index;
///
///    switch (template_type) {
///    case eDiscTemplate_11_16_Coding:
///       index = DiscontigIndex_11_16_Coding(accum);
///       break;
///    ...
///    case eDiscTemplate_12_21_Optimal:
///       index = DiscontigIndex_12_21_Optimal(accum);
///       break;
///    default:
///       index = 0;
///       break;
///    }
///
///    return index;
/// }
/// ```
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
/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:2617-2627
/// ```c
///     if (mb_lt->discontiguous) {
///         if (mb_lt->two_templates)
///             mb_lt->scansub_callback =
///                        (void *)s_MB_DiscWordScanSubject_TwoTemplates_1;
///         else if (mb_lt->template_type == eDiscTemplate_11_18_Coding)
///             mb_lt->scansub_callback = (void *)s_MB_DiscWordScanSubject_11_18_1;
///         else if (mb_lt->template_type == eDiscTemplate_11_21_Coding)
///             mb_lt->scansub_callback = (void *)s_MB_DiscWordScanSubject_11_21_1;
///         else
///             mb_lt->scansub_callback = (void *)s_MB_DiscWordScanSubject_1;
///     }
/// ```
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

/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:2202-2273
/// ```c
/// static Int4 s_MB_DiscWordScanSubject_1(const LookupTableWrap* lookup_wrap,
///        const BLAST_SequenceBlk* subject,
///        BlastOffsetPair* NCBI_RESTRICT offset_pairs, Int4 max_hits,
///        Int4* scan_range)
/// {
///    BlastMBLookupTable* mb_lt = (BlastMBLookupTable*) lookup_wrap->lut;
///    Uint1* s = subject->sequence + scan_range[0] / COMPRESSION_RATIO;
///    Int4 total_hits = 0;
///    Int4 index;
///    Uint8 accum = 0;
///    EDiscTemplateType template_type = mb_lt->template_type;
///    Int4 template_length = mb_lt->template_length;
///
///    ASSERT(lookup_wrap->lut_type == eMBLookupTable);
///    max_hits -= mb_lt->longest_chain;
///
///    /* fill the accumulator */
///    index = scan_range[0] - (scan_range[0] % COMPRESSION_RATIO);
///    while(index < scan_range[0] + template_length) {
///       accum = accum << 8 | *s++;
///       index += COMPRESSION_RATIO;
///    }
///
///    /* note that the part of the loop we jump to will
///       depend on the number of extra bases (0-3) in the
///       accumulator, and not on the value of scan_range[0] */
///    switch (index - (scan_range[0] + template_length)) {
///    case 1:
///        goto base_3;
///    case 2:
///        goto base_2;
///    case 3:
///        /* this branch of the main loop adds another
///           byte from s[], which is not needed in the
///           initial value of the accumulator */
///        accum = accum >> 8;
///        s--;
///        goto base_1;
///    }
///
///    while (scan_range[0] <= scan_range[1]) {
///
///       index = ComputeDiscontiguousIndex(accum, template_type);
///       MB_ACCESS_HITS();
///       scan_range[0]++;
///
/// base_1:
///       if (scan_range[0] > scan_range[1])
///          break;
///
///       accum = accum << 8 | *s++;
///       index = ComputeDiscontiguousIndex(accum >> 6, template_type);
///       MB_ACCESS_HITS();
///       scan_range[0]++;
///
/// base_2:
///       if (scan_range[0] > scan_range[1])
///          break;
///
///       index = ComputeDiscontiguousIndex(accum >> 4, template_type);
///       MB_ACCESS_HITS();
///       scan_range[0]++;
///
/// base_3:
///       if (scan_range[0] > scan_range[1])
///          break;
///
///       index = ComputeDiscontiguousIndex(accum >> 2, template_type);
///       MB_ACCESS_HITS();
///       scan_range[0]++;
///    }
///    return total_hits;
/// }
/// ```
fn mb_disc_word_scan_subject_1(
    packed: &[u8],
    template_type: DiscTemplateType,
    template_length: i32,
    mut scan_range: [i32; 2],
    on_word: &mut impl FnMut(i32, u32, Option<u32>),
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
    loop {
        match phase {
            0 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                on_word(
                    scan_range[0],
                    compute_discontiguous_index(accum, template_type),
                    None,
                );
                scan_range[0] += 1;
                phase = 1;
            }
            1 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                accum = (accum << 8) | packed[s] as u64;
                s += 1;
                on_word(
                    scan_range[0],
                    compute_discontiguous_index(accum >> 6, template_type),
                    None,
                );
                scan_range[0] += 1;
                phase = 2;
            }
            2 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                on_word(
                    scan_range[0],
                    compute_discontiguous_index(accum >> 4, template_type),
                    None,
                );
                scan_range[0] += 1;
                phase = 3;
            }
            _ => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                on_word(
                    scan_range[0],
                    compute_discontiguous_index(accum >> 2, template_type),
                    None,
                );
                scan_range[0] += 1;
                phase = 0;
            }
        }
    }
}

/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:2289-2367
/// ```c
/// static Int4 s_MB_DiscWordScanSubject_TwoTemplates_1(
///        const LookupTableWrap* lookup_wrap,
///        const BLAST_SequenceBlk* subject,
///        BlastOffsetPair* NCBI_RESTRICT offset_pairs, Int4 max_hits,
///        Int4* scan_range)
/// {
///    ...
///    /* fill the accumulator */
///    index = scan_range[0] - (scan_range[0] % COMPRESSION_RATIO);
///    while(index < scan_range[0] + template_length) {
///       accum = accum << 8 | *s++;
///       index += COMPRESSION_RATIO;
///    }
///    ...
///    switch (index - (scan_range[0] + template_length)) {
///    case 1:
///        goto base_3;
///    case 2:
///        goto base_2;
///    case 3:
///        ...
///        accum = accum >> 8;
///        s--;
///        goto base_1;
///    }
///
///    while (scan_range[0] <= scan_range[1]) {
///
///       index = ComputeDiscontiguousIndex(accum, template_type);
///       index2 = ComputeDiscontiguousIndex(accum, second_template_type);
///       MB_ACCESS_HITS2();
///       scan_range[0]++;
///
/// base_1:
///       if (scan_range[0] > scan_range[1])
///          break;
///
///       accum = accum << 8 | *s++;
///       index = ComputeDiscontiguousIndex(accum >> 6, template_type);
///       index2 = ComputeDiscontiguousIndex(accum >> 6, second_template_type);
///       MB_ACCESS_HITS2();
///       scan_range[0]++;
///
/// base_2:
///       if (scan_range[0] > scan_range[1])
///          break;
///
///       index = ComputeDiscontiguousIndex(accum >> 4, template_type);
///       index2 = ComputeDiscontiguousIndex(accum >> 4, second_template_type);
///       MB_ACCESS_HITS2();
///       scan_range[0]++;
///
/// base_3:
///       if (scan_range[0] > scan_range[1])
///          break;
///
///       index = ComputeDiscontiguousIndex(accum >> 2, template_type);
///       index2 = ComputeDiscontiguousIndex(accum >> 2, second_template_type);
///       MB_ACCESS_HITS2();
///       scan_range[0]++;
///    }
///    return total_hits;
/// }
/// ```
fn mb_disc_word_scan_subject_two_templates_1(
    packed: &[u8],
    template_type: DiscTemplateType,
    second_template_type: DiscTemplateType,
    template_length: i32,
    mut scan_range: [i32; 2],
    on_word: &mut impl FnMut(i32, u32, Option<u32>),
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
    let mut word = |s_off: i32, shifted: u64| {
        on_word(
            s_off,
            compute_discontiguous_index(shifted, template_type),
            Some(compute_discontiguous_index(shifted, second_template_type)),
        );
    };
    loop {
        match phase {
            0 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                word(scan_range[0], accum);
                scan_range[0] += 1;
                phase = 1;
            }
            1 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                accum = (accum << 8) | packed[s] as u64;
                s += 1;
                word(scan_range[0], accum >> 6);
                scan_range[0] += 1;
                phase = 2;
            }
            2 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                word(scan_range[0], accum >> 4);
                scan_range[0] += 1;
                phase = 3;
            }
            _ => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                word(scan_range[0], accum >> 2);
                scan_range[0] += 1;
                phase = 0;
            }
        }
    }
}

/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:2366-2479
/// ```c
/// static Int4 s_MB_DiscWordScanSubject_11_18_1(
///        const LookupTableWrap* lookup_wrap,
///        const BLAST_SequenceBlk* subject,
///        BlastOffsetPair* NCBI_RESTRICT offset_pairs, Int4 max_hits,
///        Int4* scan_range)
/// {
///    BlastMBLookupTable* mb_lt = (BlastMBLookupTable*) lookup_wrap->lut;
///    Uint1* s = subject->sequence + scan_range[0] / COMPRESSION_RATIO;
///    Int4 total_hits = 0;
///    Int4 index;
///    const Int4 kTemplateLength = 18;
///    Uint4 lo = 0;
///    Uint4 hi = 0;
///    ...
///    /* fill the accumulator */
///    index = scan_range[0] - (scan_range[0] % COMPRESSION_RATIO);
///    while(index < scan_range[0] + kTemplateLength) {
///       hi = (hi << 8) | (lo >> 24);
///       lo = lo << 8 | *s++;
///       index += COMPRESSION_RATIO;
///    }
///
///    switch (index - (scan_range[0] + kTemplateLength)) {
///    case 1:
///        goto base_3;
///    case 2:
///        goto base_2;
///    case 3:
///        s--;
///        lo = (lo >> 8) | (hi << 24);
///        hi = hi >> 8;
///        goto base_1;
///    }
///
///    while (scan_range[0] <= scan_range[1]) {
///
///       index = ((lo & 0x00000003)      ) |
///               ((lo & 0x000000f0) >>  2) |
///               ((lo & 0x00003c00) >>  4) |
///               ((lo & 0x00030000) >>  6) |
///               ((lo & 0x03c00000) >> 10) |
///               ((lo & 0xf0000000) >> 12) |
///               ((hi & 0x0000000c) << 18);
///       MB_ACCESS_HITS();
///       scan_range[0]++;
///
/// base_1:
///       if (scan_range[0] > scan_range[1])
///          break;
///
///       hi = (hi << 8) | (lo >> 24);
///       lo = lo << 8 | *s++;
///
///       index = ((lo & 0x000000c0) >>  6) |
///               ((lo & 0x00003c00) >>  8) |
///               ((lo & 0x000f0000) >> 10) |
///               ((lo & 0x00c00000) >> 12) |
///               ((lo & 0xf0000000) >> 16) |
///               ((hi & 0x0000003c) << 14) |
///               ((hi & 0x00000300) << 12);
///       MB_ACCESS_HITS();
///       scan_range[0]++;
///
/// base_2:
///       if (scan_range[0] > scan_range[1])
///          break;
///
///       index = ((lo & 0x00000030) >>  4) |
///               ((lo & 0x00000f00) >>  6) |
///               ((lo & 0x0003c000) >>  8) |
///               ((lo & 0x00300000) >> 10) |
///               ((lo & 0x3c000000) >> 14) |
///               ((hi & 0x0000000f) << 16) |
///               ((hi & 0x000000c0) << 14);
///       MB_ACCESS_HITS();
///       scan_range[0]++;
///
/// base_3:
///       if (scan_range[0] > scan_range[1])
///          break;
///
///       index = ((lo & 0x0000000c) >>  2) |
///               ((lo & 0x000003c0) >>  4) |
///               ((lo & 0x0000f000) >>  6) |
///               ((lo & 0x000c0000) >>  8) |
///               ((lo & 0x0f000000) >> 12) |
///               ((lo & 0xc0000000) >> 14) |
///               ((hi & 0x00000003) << 18) |
///               ((hi & 0x00000030) << 16);
///       MB_ACCESS_HITS();
///       scan_range[0]++;
///    }
///    return total_hits;
/// }
/// ```
fn mb_disc_word_scan_subject_11_18_1(
    packed: &[u8],
    mut scan_range: [i32; 2],
    on_word: &mut impl FnMut(i32, u32, Option<u32>),
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
    loop {
        match phase {
            0 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                let index = (lo & 0x00000003)
                    | ((lo & 0x000000f0) >> 2)
                    | ((lo & 0x00003c00) >> 4)
                    | ((lo & 0x00030000) >> 6)
                    | ((lo & 0x03c00000) >> 10)
                    | ((lo & 0xf0000000) >> 12)
                    | ((hi & 0x0000000c) << 18);
                on_word(scan_range[0], index, None);
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
                let index = ((lo & 0x000000c0) >> 6)
                    | ((lo & 0x00003c00) >> 8)
                    | ((lo & 0x000f0000) >> 10)
                    | ((lo & 0x00c00000) >> 12)
                    | ((lo & 0xf0000000) >> 16)
                    | ((hi & 0x0000003c) << 14)
                    | ((hi & 0x00000300) << 12);
                on_word(scan_range[0], index, None);
                scan_range[0] += 1;
                phase = 2;
            }
            2 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                let index = ((lo & 0x00000030) >> 4)
                    | ((lo & 0x00000f00) >> 6)
                    | ((lo & 0x0003c000) >> 8)
                    | ((lo & 0x00300000) >> 10)
                    | ((lo & 0x3c000000) >> 14)
                    | ((hi & 0x0000000f) << 16)
                    | ((hi & 0x000000c0) << 14);
                on_word(scan_range[0], index, None);
                scan_range[0] += 1;
                phase = 3;
            }
            _ => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                let index = ((lo & 0x0000000c) >> 2)
                    | ((lo & 0x000003c0) >> 4)
                    | ((lo & 0x0000f000) >> 6)
                    | ((lo & 0x000c0000) >> 8)
                    | ((lo & 0x0f000000) >> 12)
                    | ((lo & 0xc0000000) >> 14)
                    | ((hi & 0x00000003) << 18)
                    | ((hi & 0x00000030) << 16);
                on_word(scan_range[0], index, None);
                scan_range[0] += 1;
                phase = 0;
            }
        }
    }
}

/// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:2481-2598
/// ```c
/// static Int4 s_MB_DiscWordScanSubject_11_21_1(
///        const LookupTableWrap* lookup_wrap,
///        const BLAST_SequenceBlk* subject,
///        BlastOffsetPair* NCBI_RESTRICT offset_pairs, Int4 max_hits,
///        Int4* scan_range)
/// {
///    ...
///    const Int4 kTemplateLength = 21;
///    Uint4 lo = 0;
///    Uint4 hi = 0;
///    ...
///    index = scan_range[0] - (scan_range[0] % COMPRESSION_RATIO);
///    while(index < scan_range[0] + kTemplateLength) {
///       hi = (hi << 8) | (lo >> 24);
///       lo = lo << 8 | *s++;
///       index += COMPRESSION_RATIO;
///    }
///
///    switch (index - (scan_range[0] + kTemplateLength)) {
///    case 1:
///        goto base_3;
///    case 2:
///        goto base_2;
///    case 3:
///        s--;
///        lo = (lo >> 8) | (hi << 24);
///        hi = hi >> 8;
///        goto base_1;
///    }
///
///    while (scan_range[0] <= scan_range[1]) {
///
///       index = ((lo & 0x00000003)      ) |
///               ((lo & 0x000000f0) >>  2) |
///               ((lo & 0x00000c00) >>  4) |
///               ((lo & 0x000f0000) >>  8) |
///               ((lo & 0x00c00000) >> 10) |
///               ((lo & 0xf0000000) >> 14) |
///               ((hi & 0x0000000c) << 16) |
///               ((hi & 0x00000300) << 12);
///       MB_ACCESS_HITS();
///       scan_range[0]++;
///
/// base_1:
///       if (scan_range[0] > scan_range[1])
///          break;
///
///       hi = (hi << 8) | (lo >> 24);
///       lo = lo << 8 | *s++;
///
///       index = ((lo & 0x000000c0) >>  6) |
///               ((lo & 0x00003c00) >>  8) |
///               ((lo & 0x00030000) >> 10) |
///               ((lo & 0x03c00000) >> 14) |
///               ((lo & 0x30000000) >> 16) |
///               ((hi & 0x0000003c) << 12) |
///               ((hi & 0x00000300) << 10) |
///               ((hi & 0x0000c000) <<  6);
///       MB_ACCESS_HITS();
///       scan_range[0]++;
///
/// base_2:
///       if (scan_range[0] > scan_range[1])
///          break;
///
///       index = ((lo & 0x00000030) >>  4) |
///               ((lo & 0x00000f00) >>  6) |
///               ((lo & 0x0000c000) >>  8) |
///               ((lo & 0x00f00000) >> 12) |
///               ((lo & 0x0c000000) >> 14) |
///               ((hi & 0x0000000f) << 14) |
///               ((hi & 0x000000c0) << 12) |
///               ((hi & 0x00003000) <<  8);
///       MB_ACCESS_HITS();
///       scan_range[0]++;
///
/// base_3:
///       if (scan_range[0] > scan_range[1])
///          break;
///
///       index = ((lo & 0x0000000c) >>  2) |
///               ((lo & 0x000003c0) >>  4) |
///               ((lo & 0x00003000) >>  6) |
///               ((lo & 0x003c0000) >> 10) |
///               ((lo & 0x03000000) >> 12) |
///               ((lo & 0xc0000000) >> 16) |
///               ((hi & 0x00000003) << 16) |
///               ((hi & 0x00000030) << 14) |
///               ((hi & 0x00000c00) << 10);
///       MB_ACCESS_HITS();
///       scan_range[0]++;
///    }
///    return total_hits;
/// }
/// ```
fn mb_disc_word_scan_subject_11_21_1(
    packed: &[u8],
    mut scan_range: [i32; 2],
    on_word: &mut impl FnMut(i32, u32, Option<u32>),
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
    loop {
        match phase {
            0 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                let index = (lo & 0x00000003)
                    | ((lo & 0x000000f0) >> 2)
                    | ((lo & 0x00000c00) >> 4)
                    | ((lo & 0x000f0000) >> 8)
                    | ((lo & 0x00c00000) >> 10)
                    | ((lo & 0xf0000000) >> 14)
                    | ((hi & 0x0000000c) << 16)
                    | ((hi & 0x00000300) << 12);
                on_word(scan_range[0], index, None);
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
                let index = ((lo & 0x000000c0) >> 6)
                    | ((lo & 0x00003c00) >> 8)
                    | ((lo & 0x00030000) >> 10)
                    | ((lo & 0x03c00000) >> 14)
                    | ((lo & 0x30000000) >> 16)
                    | ((hi & 0x0000003c) << 12)
                    | ((hi & 0x00000300) << 10)
                    | ((hi & 0x0000c000) << 6);
                on_word(scan_range[0], index, None);
                scan_range[0] += 1;
                phase = 2;
            }
            2 => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                let index = ((lo & 0x00000030) >> 4)
                    | ((lo & 0x00000f00) >> 6)
                    | ((lo & 0x0000c000) >> 8)
                    | ((lo & 0x00f00000) >> 12)
                    | ((lo & 0x0c000000) >> 14)
                    | ((hi & 0x0000000f) << 14)
                    | ((hi & 0x000000c0) << 12)
                    | ((hi & 0x00003000) << 8);
                on_word(scan_range[0], index, None);
                scan_range[0] += 1;
                phase = 3;
            }
            _ => {
                if scan_range[0] > scan_range[1] {
                    break;
                }
                let index = ((lo & 0x0000000c) >> 2)
                    | ((lo & 0x000003c0) >> 4)
                    | ((lo & 0x00003000) >> 6)
                    | ((lo & 0x003c0000) >> 10)
                    | ((lo & 0x03000000) >> 12)
                    | ((lo & 0xc0000000) >> 16)
                    | ((hi & 0x00000003) << 16)
                    | ((hi & 0x00000030) << 14)
                    | ((hi & 0x00000c00) << 10);
                on_word(scan_range[0], index, None);
                scan_range[0] += 1;
                phase = 0;
            }
        }
    }
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_def.h:83
// ```c
// #define COMPRESSION_RATIO 4
// ```
const COMPRESSION_RATIO_I32: i32 = 4;

#[cfg(test)]
mod tests {
    use super::*;

    const ALL_TEMPLATES: [DiscTemplateType; 12] = [
        DiscTemplateType::T11_16Coding,
        DiscTemplateType::T11_16Optimal,
        DiscTemplateType::T12_16Coding,
        DiscTemplateType::T12_16Optimal,
        DiscTemplateType::T11_18Coding,
        DiscTemplateType::T11_18Optimal,
        DiscTemplateType::T12_18Coding,
        DiscTemplateType::T12_18Optimal,
        DiscTemplateType::T11_21Coding,
        DiscTemplateType::T11_21Optimal,
        DiscTemplateType::T12_21Coding,
        DiscTemplateType::T12_21Optimal,
    ];

    fn template_length(t: DiscTemplateType) -> usize {
        match t as i32 {
            1..=4 => 16,
            5..=8 => 18,
            _ => 21,
        }
    }

    fn pack(bases: &[u8]) -> Vec<u8> {
        let mut packed = vec![0u8; bases.len().div_ceil(4)];
        for (i, &b) in bases.iter().enumerate() {
            packed[i / 4] |= b << (2 * (3 - (i % 4)));
        }
        packed
    }

    fn bases(len: usize, mut seed: u64) -> Vec<u8> {
        (0..len)
            .map(|_| {
                seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((seed >> 33) & 3) as u8
            })
            .collect()
    }

    /// The index of the word at `start` read base by base, as `s_FillDiscMBTable` builds
    /// the accumulator (blast_nalookup.c:752-775).
    fn index_from_bases(seq: &[u8], start: usize, t: DiscTemplateType) -> u32 {
        let mut accum = 0u64;
        for &b in &seq[..start + template_length(t)] {
            accum = (accum << 2) | b as u64;
        }
        compute_discontiguous_index(accum, t)
    }

    // The NCBI template strings (blast_nalookup.h:195-210) read as masks over the template
    // positions: an index uses exactly the weight's bases, in order.
    #[test]
    fn templates_select_the_documented_positions() {
        let patterns: [(DiscTemplateType, &str); 12] = [
            (DiscTemplateType::T12_16Optimal, "1110110110110111"),
            (DiscTemplateType::T11_16Optimal, "1110010110110111"),
            (DiscTemplateType::T12_18Optimal, "111010110010110111"),
            (DiscTemplateType::T11_18Optimal, "111010010110010111"),
            (DiscTemplateType::T12_21Optimal, "111010010110010010111"),
            (DiscTemplateType::T11_21Optimal, "111010010100010010111"),
            (DiscTemplateType::T12_16Coding, "1111101101101101"),
            (DiscTemplateType::T11_16Coding, "1101101101101101"),
            (DiscTemplateType::T12_18Coding, "101101101101101101"),
            (DiscTemplateType::T11_18Coding, "101101100101101101"),
            (DiscTemplateType::T12_21Coding, "100101101101100101101"),
            (DiscTemplateType::T11_21Coding, "100101100101100101101"),
        ];
        for (t, pattern) in patterns {
            let len = pattern.len();
            assert_eq!(len, template_length(t));
            // Each weighted position alone, as a 3 base, gives the index bits of that
            // position; unweighted positions give no bits.
            let weighted: Vec<usize> = pattern
                .bytes()
                .enumerate()
                .filter(|&(_, c)| c == b'1')
                .map(|(i, _)| i)
                .collect();
            for i in 0..len {
                let mut seq = vec![0u8; len];
                seq[i] = 3;
                let index = index_from_bases(&seq, 0, t);
                match weighted.iter().position(|&w| w == i) {
                    Some(rank) => {
                        let shift = 2 * (weighted.len() - 1 - rank);
                        assert_eq!(index, 3 << shift, "{t:?} position {i}");
                    }
                    None => assert_eq!(index, 0, "{t:?} position {i}"),
                }
            }
        }
    }

    #[test]
    fn scans_give_the_index_of_every_word_from_any_start() {
        let seq = bases(300, 7);
        let packed = pack(&seq);
        for t in ALL_TEMPLATES {
            let tl = template_length(t) as i32;
            for kind in [
                DiscScanSubjectKind::General,
                DiscScanSubjectKind::TwoTemplates,
                DiscScanSubjectKind::Scan11_18,
                DiscScanSubjectKind::Scan11_21,
            ] {
                if (kind == DiscScanSubjectKind::Scan11_18 && t != DiscTemplateType::T11_18Coding)
                    || (kind == DiscScanSubjectKind::Scan11_21
                        && t != DiscTemplateType::T11_21Coding)
                {
                    continue;
                }
                let second = if kind == DiscScanSubjectKind::TwoTemplates {
                    if t as i32 % 2 == 0 {
                        continue;
                    }
                    t.next()
                } else {
                    DiscTemplateType::Contiguous
                };
                for start in 0..9 {
                    let end = seq.len() as i32 - tl;
                    let mut seen = Vec::new();
                    disc_word_scan_subject(
                        kind,
                        &packed,
                        t,
                        second,
                        tl,
                        [start, end],
                        &mut |s_off, index, index2| seen.push((s_off, index, index2)),
                    );
                    assert_eq!(seen.len() as i32, end - start + 1);
                    for (s_off, index, index2) in seen {
                        assert_eq!(index, index_from_bases(&seq, s_off as usize, t));
                        if kind == DiscScanSubjectKind::TwoTemplates {
                            assert_eq!(
                                index2,
                                Some(index_from_bases(&seq, s_off as usize, second))
                            );
                        } else {
                            assert_eq!(index2, None);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn template_choice_follows_weight_length_and_type() {
        assert_eq!(
            get_disc_template_type(11, 18, DiscWordType::Coding),
            DiscTemplateType::T11_18Coding
        );
        assert_eq!(
            get_disc_template_type(11, 18, DiscWordType::TwoTemplates),
            DiscTemplateType::T11_18Coding
        );
        assert_eq!(
            get_disc_template_type(12, 21, DiscWordType::Optimal),
            DiscTemplateType::T12_21Optimal
        );
        assert_eq!(
            get_disc_template_type(13, 18, DiscWordType::Coding),
            DiscTemplateType::Contiguous
        );
        assert_eq!(
            get_disc_template_type(11, 17, DiscWordType::Coding),
            DiscTemplateType::Contiguous
        );
        assert_eq!(
            DiscTemplateType::T11_18Coding.next(),
            DiscTemplateType::T11_18Optimal
        );
        assert_eq!(
            DiscTemplateType::T12_21Coding.next(),
            DiscTemplateType::T12_21Optimal
        );
    }
}
