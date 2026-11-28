//! NCBI's accession guide (`SAccGuide` of `c++/src/objects/seqloc/Seq_id.cpp`) and
//! `CSeq_id::IdentifyAccession`: the type of a bare accession such as `XP_011520000.1` (a
//! RefSeq protein) or `P69905` (Swiss-Prot), as `makeblastdb -parse_seqids` reads it.
//!
//! The rules are NCBI's built-in guide (`kBuiltInGuide`, `accguide2.inc` of 2.17.0, which
//! makeblastdb uses: no `accguide2.txt` data file is looked for, so the result does not depend
//! on the machine). They are kept in two files beside this module, made by
//! `LOSAT/tests/accguide/make_accguide.py`: `rules.txt` (every line that `AddRule` reads,
//! without its comment) and `ranges.bin` (the range lines of `special2` blocks, four bits per
//! character). The prefix rules are read into the guide once, at the first lookup; the specials
//! of a format (its `special` lines and `special2` blocks, almost all of the text) once, at the
//! first lookup that reads them, which few accessions do (see `Part`).
//!
//! What is not ported: the branches of `AddRule` for a type name that this build does not
//! know (the fallbacks: every name of the built-in guide is known, which the tests check),
//! the diagnostics of the loading (`ERR_POST_X` 2-21: the built-in guide raises none of
//! Warning severity, also checked, by the build of every part that makes them; the Info ones
//! are not shown by makeblastdb), and the fallback warning of `IdentifyAccession` (no value of
//! the built-in guide is a fallback).

mod acc_info;

use std::collections::BTreeMap;
use std::ops::Bound;
use std::sync::OnceLock;

pub use acc_info::*;

use super::model::ParsedSeqId;
// The masks of the guide are letters, `_` and `?` (`NStr::MatchesMask`); the loading refuses a mask
// with `*`, `[` or `\`, and a NUL character of an accession matches nothing.
use crate::ncbi_str::matches_mask;
use crate::utils::byte_search::find_byte;
use crate::utils::byte_search::find_first;

/// A `CSeq_id::EAccessionInfo` value.
pub type AccInfo = u32;

/// `CSeq_id::GetAccType`: the Seq-id type of an accession's value (an `E_Choice`, 0 for
/// `e_not_set`).
///
/// NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_id.hpp:562-563
pub fn acc_type(info: AccInfo) -> u32 {
    info & TYPE_MASK
}

/// The parse flags of `CSeq_id` (`EParseFlags`) that the accession functions read.
///
/// NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_id.hpp:75-103
pub mod parse_flags {
    pub const PARTIAL_OK: u32 = 0x01;
    pub const RAW_TEXT: u32 = 0x02;
    pub const RAW_GI: u32 = 0x04;
    pub const ANY_RAW: u32 = RAW_TEXT | RAW_GI;
    pub const VALID_LOCAL: u32 = 0x08;
    pub const ANY_LOCAL: u32 = 0x18;
    pub const NO_FASTA: u32 = 0x20;
    pub const FALLBACK_OK: u32 = 0x40;
    pub const DEFAULT: u32 = RAW_TEXT | VALID_LOCAL;
}
use parse_flags::RAW_GI;
use parse_flags::RAW_TEXT;

/// A diagnostic that NCBI's Seq-id code posts (`ERR_POST_X` with a subcode, as a warning carrying
/// the message, in Seq_id.cpp, whose error code is `NCBI_ERRCODE_X` 1306): the caller writes it as
/// `Warning: (1306.<subcode>) [<application>] <message>`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeqIdWarning {
    pub subcode: u32,
    /// The warning's text, as bytes (an id's text keeps its bytes: probe R3 of line LIB-B).
    pub message: Vec<u8>,
}

/// A set of bit positions as `bm::bvector<>` holds them for a small special: sorted,
/// disjoint, non-adjacent closed ranges (the guide's bit vectors are up to 10^9 bits long
/// and mostly empty).
#[derive(Default, Debug)]
struct BitRanges(Vec<(u32, u32)>);

/// `partition_point` of `items` and `before` (true for a prefix of `items`, false for the
/// rest), searched from the end: the items 1, 2, 4, ... places before the end are tested until
/// one is before the point, and the span after it is then searched in halves. The cost grows
/// with the distance of the point from the end, which is 0 when a range comes after the ones of
/// its set (the guide lists its ranges in increasing order).
fn partition_point_from_end<T>(items: &[T], before: impl Fn(&T) -> bool) -> usize {
    let mut end = items.len();
    let mut step = 1;
    while end > 0 {
        let probe = end.saturating_sub(step);
        if before(&items[probe]) {
            return probe + 1 + items[probe + 1..end].partition_point(&before);
        }
        end = probe;
        step *= 2;
    }
    0
}

impl BitRanges {
    /// `set_range(left, right)`: the bits from `left` to `right`, or from `right` to `left`.
    ///
    /// BitMagic 8.0.1 (NCBI 598d8ae6 copy): c++/include/util/bitset/bm.h:2377-2378
    fn set_range(&mut self, left: u32, right: u32) {
        let (left, right) = (left.min(right), left.max(right));
        let ranges = &mut self.0;
        // A range after every range is appended (most ranges: the guide lists the ranges of a
        // key in increasing order).
        if ranges
            .last()
            .is_none_or(|&(_, end)| end.saturating_add(1) < left)
        {
            ranges.push((left, right));
            return;
        }
        // The ranges that touch or overlap [left, right].
        let first = partition_point_from_end(ranges, |&(_, end)| end.saturating_add(1) < left);
        let last = partition_point_from_end(ranges, |&(start, _)| start <= right.saturating_add(1));
        if first == last {
            ranges.insert(first, (left, right));
        } else {
            let start = left.min(ranges[first].0);
            let end = right.max(ranges[last - 1].1);
            ranges[first] = (start, end);
            ranges.drain(first + 1..last);
        }
    }

    /// `clear_range(left, right)`, which is `set_range` to 0 (so with `left` and `right` in
    /// either order).
    ///
    /// BitMagic 8.0.1 (NCBI 598d8ae6 copy): c++/include/util/bitset/bm.h:1216-1217,2377-2378
    fn clear_range(&mut self, left: u32, right: u32) {
        let (left, right) = (left.min(right), left.max(right));
        let ranges = &mut self.0;
        let first = partition_point_from_end(ranges, |&(_, end)| end < left);
        let last = partition_point_from_end(ranges, |&(start, _)| start <= right);
        if first == last {
            return;
        }
        let mut kept = Vec::with_capacity(2);
        if ranges[first].0 < left {
            kept.push((ranges[first].0, left - 1));
        }
        if ranges[last - 1].1 > right {
            kept.push((right + 1, ranges[last - 1].1));
        }
        ranges.splice(first..last, kept);
    }

    /// `operator[]` (`test`).
    fn test(&self, bit: u32) -> bool {
        let index = self.0.partition_point(|&(_, end)| end < bit);
        index < self.0.len() && self.0[index].0 <= bit
    }

    /// `any()`.
    fn any(&self) -> bool {
        !self.0.is_empty()
    }
}

/// `SAccGuide::SSubMap`: the rules of one format (`<letters>+<digits>`).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:889-911
/// The specials are made at the first lookup that reads them (see [`Part`]); a build that
/// makes them sets them at once.
#[derive(Debug)]
struct SubMap {
    prefixes: BTreeMap<Vec<u8>, AccInfo>,
    wildcards: Vec<(Vec<u8>, AccInfo)>,
    specials: OnceLock<Specials>,
}

impl SubMap {
    /// The rules of a format that a build of `part` adds: with empty specials when the build
    /// makes specials (of this format or another: a lookup during the build then never makes
    /// them), without them for [`Part::Rules`].
    fn new(part: Part) -> Self {
        SubMap {
            prefixes: BTreeMap::new(),
            wildcards: Vec::new(),
            specials: match part {
                Part::Rules => OnceLock::new(),
                Part::Specials(_) => OnceLock::from(Specials::default()),
                #[cfg(test)]
                Part::All => OnceLock::from(Specials::default()),
            },
        }
    }
}

/// The specials of one format: `big_specials` and `small_specials` of `SSubMap`.
///
/// The multimap of the small specials is a map of each key to its options in the multimap's
/// order (an `emplace` adds an option at the end of its key's equal range), kept as the place
/// of the options in `small_options`, where the key is kept too: `SHints` holds that place.
/// The cleanup at the end of a build removes from `small_specials` a key whose options all
/// became empty; its place in `small_options` stays, and no lookup reaches it.
#[derive(Default, Debug)]
struct Specials {
    big_specials: BTreeMap<Vec<u8>, (Vec<u8>, AccInfo)>,
    small_specials: BTreeMap<Vec<u8>, usize>,
    small_options: Vec<(Vec<u8>, Vec<(BitRanges, AccInfo)>)>,
}

impl Specials {
    /// The options of the small special `key`: its equal range of the multimap.
    fn options(&self, key: &[u8]) -> &[(BitRanges, AccInfo)] {
        match self.small_specials.get(key) {
            Some(&slot) => &self.small_options[slot].1,
            None => &[],
        }
    }

    /// The place of the options of `key`, made empty for a new key.
    fn slot(&mut self, key: &[u8]) -> usize {
        if let Some(&slot) = self.small_specials.get(key) {
            return slot;
        }
        let slot = self.small_options.len();
        self.small_options.push((key.to_vec(), Vec::new()));
        self.small_specials.insert(key.to_vec(), slot);
        slot
    }

    /// The specials of `fmt` in the built-in guide (a build of [`Part::Specials`]).
    fn built_in(fmt: u32) -> Specials {
        AccGuide::built_in(Part::Specials(fmt))
            .rules
            .remove(&fmt)
            .and_then(|submap| submap.specials.into_inner())
            .expect("a build of the specials of a format with a prefix rule has its specials")
    }
}

/// The part of the guide that a build reads the rules for.
///
/// The prefixes, the wildcards and the general databases are read at the first lookup
/// ([`guide`]); the specials of a format (the `special` lines and the `special2` blocks, almost
/// all of the guide's text) at the first lookup that reads them (`Find` reads them for an
/// accession only when its prefix or wildcard has `fAcc_specials`). A build of the specials of a
/// format reads every line in order, so that `SHints` holds what it holds in NCBI's loading when
/// a line of that format is added; it adds no specials of other formats, which never change
/// those of this one: every line that adds a special first makes its format the current one
/// (`FindSubMap`), which forgets the small special of another format. It keeps no prefix
/// either, and it passes over the range lines of the other formats' blocks.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Part {
    /// The rules without the specials (the `special` and `special2` lines are not read).
    Rules,
    /// Every rule, with the specials of one format.
    Specials(u32),
    /// Every rule, with every special, and the diagnostics of the loading: the build of the
    /// tests, which check that the built-in guide raises none.
    #[cfg(test)]
    All,
}

impl Part {
    /// Whether a build of this part adds the specials of `fmt`.
    fn has_specials_of(self, fmt: u32) -> bool {
        match self {
            Part::Rules => false,
            Part::Specials(only) => only == fmt,
            #[cfg(test)]
            Part::All => true,
        }
    }

    /// Whether a build of this part keeps the prefixes and wildcards (a build of the specials of
    /// a format keeps only them).
    fn keeps_prefixes(self) -> bool {
        !matches!(self, Part::Specials(_))
    }
}

/// `SAccGuide`.
#[derive(Debug)]
pub struct AccGuide {
    rules: BTreeMap<u32, SubMap>,
    general: BTreeMap<Vec<u8>, AccInfo>,
    /// The diagnostics of Warning severity that NCBI's loading would post (none for the
    /// built-in guide).
    load_warnings: Vec<String>,
    /// What this guide was built with.
    part: Part,
}

/// `SAccGuide::SHints`: the state that `AddRule` keeps from one line to the next.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:907-945
/// `prev_small_special` is the place of the key's options in the current format's
/// `small_options` and the position of the option in its equal range (`None` for `end()`);
/// `prev_big_special` is only an insertion hint and is not kept. `prev_special2_acc`, the
/// running accession of a `special2` block, is kept as its key and its bit (`special2_key`,
/// `special2_bit`): every range of the built-in guide replaces digits of the bit only (see
/// `add_special2_ranges`). The members that serve unknown type names only (the fallbacks:
/// `prev_special_type` and its name) are not kept; those that only `FindSpecial` reads are
/// kept by the build of the tests (`SpecialHints`).
struct Hints {
    prev_type: AccInfo,
    prev_type_name: Vec<u8>,
    prev_submap: Option<u32>,
    prev_small_special: Option<(usize, usize)>,
    special2_key: Vec<u8>,
    special2_bit: u32,
    special2_submap: Option<u32>,
    special2_type: AccInfo,
    version: u32,
    #[cfg(test)]
    special: SpecialHints,
}

impl Default for Hints {
    fn default() -> Self {
        Hints {
            prev_type: UNKNOWN,
            prev_type_name: Vec::new(),
            prev_submap: None,
            prev_small_special: None,
            special2_key: Vec::new(),
            special2_bit: 0,
            special2_submap: None,
            special2_type: UNKNOWN,
            version: 1,
            #[cfg(test)]
            special: SpecialHints::default(),
        }
    }
}

/// The members of `SHints` that only `FindSpecial` reads, for the diagnostics of the loading
/// (`prev_special_format`, `prev_special_key`, `prev_special_base_key` and
/// `prev_special_base_type`).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:907-931
#[cfg(test)]
struct SpecialHints {
    format: u32,
    key: Vec<u8>,
    base_key: Vec<u8>,
    base_type: AccInfo,
}

#[cfg(test)]
impl Default for SpecialHints {
    fn default() -> Self {
        SpecialHints {
            format: 0,
            key: Vec::new(),
            base_key: Vec::new(),
            base_type: UNKNOWN,
        }
    }
}

/// `kMaxSmallSpecialDigits` and `kBVSizes`.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:866-884
const MAX_SMALL_SPECIAL_DIGITS: u32 = 9;
const BV_SIZES: [u32; 10] = [
    1,
    10,
    100,
    1_000,
    10_000,
    100_000,
    1_000_000,
    10_000_000,
    100_000_000,
    1_000_000_000,
];

/// `SAccGuide::s_Key`, given the number of `letters` and of `digits`.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:952-953
fn format_key(letters: usize, digits: usize) -> u32 {
    ((letters as u32 & 0xffff) << 16) | (digits as u32 & 0xffff)
}

/// `NStr::StringToUInt` with `NStr::fConvErr_NoThrow` (`ncbi_str::string_to_uint`, 0 for an
/// error) of the guide's numbers (`version`, the `<letters>+<digits>` of a rule) and of a bit
/// position (`NStr::StringToNumeric` into a value, `x_SplitSpecial`). These are ASCII digits
/// alone: the guide is compiled in (`rules.txt`), and the bit positions of an accession are cut
/// from the digits after its letters (`x_IdentifyAccession` keys the rules by their counts); on
/// such text NCBI's throwing `StringToNumeric` forms read the same values.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1048
fn bit_number(digits: &[u8]) -> u32 {
    crate::ncbi_str::string_to_uint(digits).unwrap_or(0)
}

/// The size of an accession of the format `fmt` (`x_SplitSpecial`'s normal size).
fn normal_size(fmt: u32) -> usize {
    (fmt >> 16) as usize + (fmt & 0xffff) as usize
}

/// The number of last digits that `x_SplitSpecial` reads as the bit of an accession of the
/// format's size, at most `kMaxSmallSpecialDigits`; the letters and digits before them key the
/// small special.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1599-1605
fn bit_digits(fmt: u32) -> u32 {
    let digits = fmt & 0xffff;
    if digits == MAX_SMALL_SPECIAL_DIGITS + 1 {
        digits - 2
    } else {
        digits.min(MAX_SMALL_SPECIAL_DIGITS)
    }
}

/// `SAccGuide::x_SplitSpecial` of `acc` and `fmt`: the length of the prefix that keys the special
/// and the number that is its bit.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1596-1625
fn split_special(acc: &[u8], fmt: u32) -> (usize, u32) {
    let raw_digits = fmt & 0xffff;
    let digits = bit_digits(fmt);
    if acc.len() == normal_size(fmt) {
        let pos = acc.len() - digits as usize;
        (pos, bit_number(&acc[pos..]))
    } else {
        let mut pos = (fmt >> 16) as usize + 2;
        let mut result = bit_number(acc.get(pos + 1..).unwrap_or_default());
        if digits == raw_digits {
            pos -= 2;
            let high = bit_number(&acc[pos..(pos + 2).min(acc.len())]);
            result = result.wrapping_add(high.wrapping_mul(BV_SIZES[digits as usize - 2]));
        }
        (pos, result)
    }
}

/// The digits of a range `<low>[-<high>]` as numbers with their counts (an empty high when
/// there is no `-`), when both are ASCII digits alone, at most `max` (at most 9) of them.
///
/// Follows the first-hyphen split of a `special2` range token into low and optional high
/// suffixes, parsing each side directly as base-10 digits.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1310
fn range_numbers(range: &[u8], max: usize) -> Option<[(u32, usize); 2]> {
    let mut parts = [(0, 0); 2];
    let mut part = 0;
    for &c in range {
        if c == b'-' && part == 0 {
            part = 1;
            continue;
        }
        let digit = c.wrapping_sub(b'0');
        let (number, count) = &mut parts[part];
        if digit > 9 || *count == max {
            return None;
        }
        *number = *number * 10 + u32::from(digit);
        *count += 1;
    }
    Some(parts)
}

/// `number` (of at least `count` digits) with its last `count` digits (at most 9) replaced by
/// those of `digits`.
///
/// Computes numerically the trailing-digit overwrite at the end of the running `special2`
/// accession for each low and high range endpoint.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1308,1311,1317
fn with_last_digits(number: u32, (digits, count): (u32, usize)) -> u32 {
    // Divisions by constants (multiplications in the machine code).
    let cleared = match count {
        0 => number,
        1 => number / 10 * 10,
        2 => number / 100 * 100,
        3 => number / 1_000 * 1_000,
        4 => number / 10_000 * 10_000,
        5 => number / 100_000 * 100_000,
        6 => number / 1_000_000 * 1_000_000,
        7 => number / 10_000_000 * 10_000_000,
        8 => number / 100_000_000 * 100_000_000,
        9 => number / 1_000_000_000 * 1_000_000_000,
        _ => 0,
    };
    cleared + digits
}

impl AccGuide {
    /// `SAccGuide::Find` (with `fmt`, `acc_or_pfx` and `key_used`): the value of an accession (or of a
    /// prefix) of the format `fmt`. `key_used` gets the wildcard or the prefix that gave
    /// it, as NCBI's (for the loading's diagnostics).
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1370-1441
    fn find(&self, fmt: u32, acc_or_pfx: &[u8], mut key_used: Option<&mut Vec<u8>>) -> AccInfo {
        let Some(submap) = self.rules.get(&fmt) else {
            return UNKNOWN;
        };
        let letters = ((fmt >> 16) as usize).min(acc_or_pfx.len());
        let pfx = &acc_or_pfx[..letters];
        let mut result = UNKNOWN;
        if let Some(&value) = submap.prefixes.get(pfx) {
            result = value;
        } else {
            for (mask, value) in &submap.wildcards {
                if !matches_mask(pfx, mask) {
                    continue;
                }
                let bad_match = mask.iter().enumerate().any(|(pos, &m)| {
                    m == b'?' && {
                        let c = pfx.get(pos).copied().unwrap_or(0);
                        !c.is_ascii_alphanumeric() && c != b'?'
                    }
                });
                if bad_match {
                    continue;
                }
                if let Some(key_used) = key_used.as_deref_mut() {
                    if acc_or_pfx != mask.as_slice() {
                        *key_used = mask.clone();
                    }
                }
                result = *value;
                break;
            }
        }
        if acc_or_pfx != pfx && result & F_SPECIALS != 0 {
            let specials = submap.specials.get_or_init(|| Specials::built_in(fmt));
            let (prefix_len, n) = split_special(acc_or_pfx, fmt);
            let prefix = &acc_or_pfx[..prefix_len];
            for (bits, value) in specials.options(prefix) {
                if bits.test(n) {
                    if let Some(key_used) = key_used {
                        key_used.clear();
                    }
                    return *value;
                }
            }
            if let Some((_, (first, value))) = specials
                .big_specials
                .range::<[u8], _>((Bound::Included(acc_or_pfx), Bound::Unbounded))
                .next()
            {
                if acc_or_pfx >= first.as_slice() {
                    if let Some(key_used) = key_used {
                        key_used.clear();
                    }
                    return *value;
                }
            }
            if let Some(key_used) = key_used {
                if key_used.is_empty() {
                    *key_used = prefix[..letters.min(prefix.len())].to_vec();
                }
            }
        }
        result
    }

    /// `SHints::FindAccInfo` of `name`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:974-989
    /// `None` is `kUnrecognized`.
    fn find_acc_info_hinted(hints: &mut Hints, name: &[u8]) -> Option<AccInfo> {
        if !hints.prev_type_name.is_empty() && name == hints.prev_type_name.as_slice() {
            return Some(hints.prev_type);
        }
        let value = find_acc_info(name)?;
        #[cfg(test)]
        hints.special.key.clear();
        hints.prev_type_name = name.to_ascii_lowercase();
        hints.prev_type = value;
        Some(value)
    }

    /// `SHints::FindSpecial` of the guide, `fmt` and an accession or range.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:991-1010
    /// (The members of the fallbacks that it sets are not kept: see [`Hints`].)
    #[cfg(test)]
    fn find_special(&self, hints: &mut Hints, fmt: u32, acc_or_range: &[u8]) -> AccInfo {
        let pfx = &acc_or_range[..((fmt >> 16) as usize).min(acc_or_range.len())];
        let special = &mut hints.special;
        if fmt == special.format {
            if acc_or_range == special.key.as_slice() {
                return hints.prev_type;
            } else if pfx == special.base_key.as_slice() {
                return special.base_type;
            }
        }
        special.format = fmt;
        special.base_key = pfx.to_vec();
        special.base_type = self.find(fmt, pfx, None);
        special.base_type
    }

    /// The diagnostic of `AddRule` for a special (`acc`, of the format `fmt`) whose prefix rule
    /// does not list specials (`FindSpecial` of `acc`, then `Find` of `first`), shown as
    /// `shown`. Only the build of every part makes it (the tests): it changes no rule, and
    /// the built-in guide raises none.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1139,1141-1152,1225-1238
    #[cfg(test)]
    fn check_specials_listed(
        &mut self,
        hints: &mut Hints,
        fmt: u32,
        acc: &[u8],
        first: &[u8],
        shown: &[u8],
    ) {
        if self.part != Part::All {
            return;
        }
        // `Find` reads the prefixes, which this build keeps.
        debug_assert!(self.part.keeps_prefixes());
        if self.find_special(hints, fmt, acc) & F_SPECIALS != 0 {
            return;
        }
        let mut key_used = Vec::new();
        self.find(fmt, first, Some(&mut key_used));
        if !key_used.is_empty() {
            self.load_warnings.push(format!(
                "Main listing for special {} doesn't indicate that specials are present.",
                String::from_utf8_lossy(shown)
            ));
        }
    }

    /// `SHints::FindSubMap` of the rules and `fmt`: the rules of `fmt`, made if absent.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1012-1028
    fn find_sub_map(&mut self, hints: &mut Hints, fmt: u32) -> u32 {
        if hints.prev_submap != Some(fmt) {
            let part = self.part;
            self.rules.entry(fmt).or_insert_with(|| SubMap::new(part));
            hints.prev_submap = Some(fmt);
            hints.prev_small_special = None;
        }
        fmt
    }

    /// The specials of `fmt` that a build adds to (made by `FindSubMap` of a build that makes
    /// them).
    fn specials_mut(&mut self, fmt: u32) -> &mut Specials {
        self.rules
            .get_mut(&fmt)
            .and_then(|submap| submap.specials.get_mut())
            .expect("FindSubMap made the rules and their specials")
    }

    /// `SAccGuide::x_AddSpecial` with the submap of `fmt` (its `specials`), hints, `fmt`, `from`,
    /// `to`, `value`, `old_name` and `new_name`, for a known type name (`old_name` null).
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1513-1594
    fn add_special(
        specials: &mut Specials,
        hints: &mut Hints,
        fmt: u32,
        from: &[u8],
        to: &[u8],
        value: AccInfo,
    ) {
        let (from_len, left) = split_special(from, fmt);
        let (to_len, right) = split_special(to, fmt);
        let (from_pfx, to_pfx) = (&from[..from_len], &to[..to_len]);
        if from_pfx != to_pfx {
            let entry = specials
                .big_specials
                .entry(to.to_vec())
                .or_insert_with(|| (from.to_vec(), value));
            entry.1 = value;
            return;
        }
        Self::add_small_special(specials, hints, from_pfx, left, right, value);
    }

    /// The small special of `x_AddSpecial`: the bits `left` to `right` of the key `from_pfx`
    /// get `value`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1541-1587
    /// A bit vector needs no size here (`BitRanges`); `x_SplitSpecial` never gives a bit past
    /// it.
    fn add_small_special(
        specials: &mut Specials,
        hints: &mut Hints,
        from_pfx: &[u8],
        left: u32,
        right: u32,
        value: AccInfo,
    ) {
        let hint = hints
            .prev_small_special
            .filter(|&(slot, _)| specials.small_options[slot].0 == from_pfx);
        let (slot, prev) = match hint {
            Some((slot, prev)) => (slot, Some(prev)),
            None => (specials.slot(from_pfx), None),
        };
        let options = &mut specials.small_options[slot].1;
        let index = Self::add_to_options(options, prev, left, right, value);
        hints.prev_small_special = Some((slot, index));
    }

    /// The bits `left` to `right` of a small special's key get `value` among the key's
    /// `options`, given the option of the hint (`prev`, when the hint is of this key); the
    /// option set, which the hint becomes.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1546-1583
    /// When the option of the hint has the value, NCBI clears the bits there and sets them
    /// again in the same option: the bits are those of setting them alone.
    fn add_to_options(
        options: &mut Vec<(BitRanges, AccInfo)>,
        prev: Option<usize>,
        left: u32,
        right: u32,
        value: AccInfo,
    ) -> usize {
        let same = |option: &(BitRanges, AccInfo)| option.1 & !F_FALLBACK == value & !F_FALLBACK;
        let mut index = match prev {
            Some(prev) if same(&options[prev]) => prev,
            Some(prev) => {
                options[prev].0.clear_range(left, right);
                let mut it = prev;
                loop {
                    if same(&options[it]) {
                        break it;
                    }
                    if it == 0 {
                        break prev + 1;
                    }
                    it -= 1;
                }
            }
            None => 0,
        };
        while index < options.len() && !same(&options[index]) {
            index += 1;
        }
        if index == options.len() {
            options.push((BitRanges::default(), value));
        }
        options[index].0.set_range(left, right);
        options[index].1 = value;
        index
    }

    /// `SAccGuide::AddRule` (with a rule and `hints`) for the lines of `rules.txt`; `ranges` gives the
    /// tokens of each `:` line of a `special2` block.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1030-1053
    fn add_rule(&mut self, tokens: &[&[u8]], hints: &mut Hints) {
        let Some(&first) = tokens.first() else {
            return;
        };
        if tokens.len() == 2 && first.eq_ignore_ascii_case(b"version") {
            hints.special2_submap = None;
            hints.version = bit_number(tokens[1]);
            if !(1..=2).contains(&hints.version) {
                self.load_warnings.push(format!(
                    "Unsupported version {}",
                    String::from_utf8_lossy(tokens[1])
                ));
            }
            return;
        }
        if let Some(plus) = first.iter().position(|&c| c == b'+') {
            if tokens.len() == 3 || (tokens.len() == 4 && tokens[3] == b"*") {
                self.add_prefix_rule(tokens, plus, hints);
                return;
            }
        }
        if tokens.len() == 3 && first.eq_ignore_ascii_case(b"special") {
            self.add_special_rule(tokens, hints);
        } else if tokens.len() >= 3 && first.eq_ignore_ascii_case(b"special2") {
            self.add_special2_rule(tokens, hints);
        } else if tokens.len() >= 2 && first == b":" {
            self.add_special2_ranges(tokens[1..].iter().copied(), hints);
        } else if tokens.len() == 3 && first.eq_ignore_ascii_case(b"gnl") {
            // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1324-1329,1349-1351
            hints.special2_submap = None;
            match Self::find_acc_info_hinted(hints, tokens[2]) {
                Some(value) => {
                    self.general.insert(tokens[1].to_ascii_uppercase(), value);
                }
                None => self.unknown_name(tokens[2]),
            }
        } else if tokens.len() == 3 && first.eq_ignore_ascii_case(b"fallback") {
            // Default fallbacks serve unknown type names only (Seq_id.cpp:1348-1363).
            hints.special2_submap = None;
        } else {
            self.load_warnings.push(format!(
                "ignoring invalid line: {}",
                String::from_utf8_lossy(&tokens.join(&b' '))
            ));
        }
    }

    /// A type name this build does not know: NCBI's fallback branches, not ported (see the
    /// module).
    fn unknown_name(&mut self, name: &[u8]) {
        self.load_warnings.push(format!(
            "unrecognized accession type {}",
            String::from_utf8_lossy(name)
        ));
    }

    /// The `<letters>+<digits> <prefix> <type> [*]` line of `AddRule`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1054-1066,1102-1131
    fn add_prefix_rule(&mut self, tokens: &[&[u8]], plus: usize, hints: &mut Hints) {
        hints.special2_submap = None;
        let fmt = format_key(
            bit_number(&tokens[0][..plus]) as usize,
            bit_number(&tokens[0][plus + 1..]) as usize,
        );
        let Some(mut value) = Self::find_acc_info_hinted(hints, tokens[2]) else {
            self.unknown_name(tokens[2]);
            return;
        };
        self.find_sub_map(hints, fmt);
        if !self.part.keeps_prefixes() {
            return;
        }
        let submap = self.rules.get_mut(&fmt).expect("FindSubMap made the rules");
        if tokens.len() == 4 {
            value |= F_SPECIALS;
        }
        let prefix = tokens[1];
        if !prefix.iter().any(|&c| c == b'?' || c == b'*') {
            submap.prefixes.insert(prefix.to_vec(), value);
        } else if let Some(wildcard) = submap.wildcards.iter_mut().find(|(mask, _)| mask == prefix)
        {
            wildcard.1 = value;
        } else {
            if prefix.iter().any(|&c| matches!(c, b'*' | b'[' | b'\\')) {
                self.load_warnings.push(format!(
                    "mask {} is not one of the guide's forms",
                    String::from_utf8_lossy(prefix)
                ));
            }
            submap.wildcards.push((prefix.to_vec(), value));
        }
    }

    /// The `special <acc>[-<acc>] <type>` line of `AddRule`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1133-1154,1193-1209
    fn add_special_rule(&mut self, tokens: &[&[u8]], hints: &mut Hints) {
        if self.part == Part::Rules {
            return;
        }
        hints.special2_submap = None;
        let range = tokens[1];
        let pos = range
            .iter()
            .position(u8::is_ascii_digit)
            .unwrap_or(range.len());
        let pos2 = range[pos..]
            .iter()
            .position(|&c| c == b'-')
            .map(|offset| pos + offset);
        let fmt = format_key(pos, pos2.unwrap_or(range.len()) - pos);
        let first = &range[..pos2.unwrap_or(range.len())];
        #[cfg(test)]
        self.check_specials_listed(hints, fmt, range, first, range);
        match Self::find_acc_info_hinted(hints, tokens[2]) {
            None => self.unknown_name(tokens[2]),
            Some(value) => {
                #[cfg(test)]
                {
                    hints.special.key = range.to_vec();
                }
                self.find_sub_map(hints, fmt);
                if self.part.has_specials_of(fmt) {
                    let last = pos2.map_or(range, |pos2| &range[pos2 + 1..]);
                    Self::add_special(self.specials_mut(fmt), hints, fmt, first, last, value);
                }
            }
        }
    }

    /// The `special2 <letters>+<digits> <type>` line of `AddRule`, which starts a block of
    /// range lines.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1210-1247
    /// A known name in the third token is the only case without a fallback (the elided
    /// lines, 1248-1294).
    fn add_special2_rule(&mut self, tokens: &[&[u8]], hints: &mut Hints) {
        if self.part == Part::Rules {
            return;
        }
        if hints.version < 2 {
            self.load_warnings
                .push("special2 valid only in version 2+ guides".to_string());
        }
        let spec = tokens[1];
        let plus = spec.iter().position(|&c| c == b'+').unwrap_or(spec.len());
        let letters = &spec[..plus];
        let digits = bit_number(spec.get(plus + 1..).unwrap_or_default()) as usize;
        let fmt = format_key(letters.len(), digits);
        #[cfg(test)]
        {
            hints.special.format = fmt;
        }
        hints.special2_submap = Some(self.find_sub_map(hints, fmt));
        // The block's accession, the letters and `digits` zeros: its key, and its bit 0.
        let key_len = normal_size(fmt) - bit_digits(fmt) as usize;
        hints.special2_key.clear();
        hints.special2_key.extend_from_slice(letters);
        hints.special2_key.resize(key_len, b'0');
        hints.special2_bit = 0;
        #[cfg(test)]
        {
            let mut acc = letters.to_vec();
            acc.resize(letters.len() + digits, b'0');
            self.check_specials_listed(hints, fmt, &acc, &acc, spec);
        }
        match Self::find_acc_info_hinted(hints, tokens[2]) {
            Some(value) => hints.special2_type = value,
            None => self.unknown_name(tokens[2]),
        }
    }

    /// The `: <range> ...` lines of a `special2` block: each range replaces the last digits
    /// of the block's running accession.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1292-1323
    /// Every range of the built-in guide is `<digits>[-<digits>]`, each side at most
    /// `bit_digits` digits (the digits that `x_SplitSpecial` reads as the bit; `make_accguide.py`
    /// refuses any other, and a test checks the files). Such a range keeps the accession's size
    /// and key, so `from` and `to` are a small special of that key, from the bit after the low
    /// digits to the bit after the high ones, and the bit is kept as a number.
    fn add_special2_ranges<'t>(
        &mut self,
        ranges: impl Iterator<Item = &'t [u8]>,
        hints: &mut Hints,
    ) {
        if hints.version < 2 {
            self.load_warnings
                .push("special2 continuation lines valid only in version 2+ guides".to_string());
        }
        let Some(fmt) = hints.special2_submap else {
            self.load_warnings
                .push("ignoring misplaced special2 ranges line.".to_string());
            return;
        };
        if !self.part.has_specials_of(fmt) {
            return;
        }
        let value = hints.special2_type;
        let digits = bit_digits(fmt) as usize;
        let specials = self.specials_mut(fmt);
        let key = std::mem::take(&mut hints.special2_key);
        let mut bit = hints.special2_bit;
        // The hint while it is known to be of this key (set by the range before).
        let mut known: Option<(usize, usize)> = None;
        for range in ranges {
            let [low, high] = range_numbers(range, digits)
                .expect("every special2 range of the built-in guide is of the bit's digits");
            let left = with_last_digits(bit, low);
            let right = with_last_digits(left, high);
            bit = right;
            match known {
                Some((slot, prev)) => {
                    let options = &mut specials.small_options[slot].1;
                    let index = Self::add_to_options(options, Some(prev), left, right, value);
                    hints.prev_small_special = Some((slot, index));
                }
                None => Self::add_small_special(specials, hints, &key, left, right, value),
            }
            known = hints.prev_small_special;
        }
        hints.special2_key = key;
        hints.special2_bit = bit;
    }

    /// `SAccGuide::SAccGuide()` with the built-in rules, then the cleanup of empty small
    /// specials and `x_InitGeneral`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1443-1498
    /// The comments and blank lines of `kBuiltInGuide` are not in `rules.txt` (`AddRule`
    /// returns at once for them), and no data file is looked for (`g_FindDataFile`):
    /// makeblastdb 2.17.0 reads the built-in guide when no newer `accguide2.txt` is on its
    /// data path, and LOSAT always reads it. A build of `part` passes over the range lines of
    /// the `special2` blocks it adds nothing for.
    fn built_in(part: Part) -> AccGuide {
        let mut guide = AccGuide {
            rules: BTreeMap::new(),
            general: BTreeMap::new(),
            load_warnings: Vec::new(),
            part,
        };
        let mut hints = Hints::default();
        let mut ranges = RangeLines::new(RANGES);
        let mut words: Vec<&[u8]> = Vec::new();
        for line in rule_lines(RULES.as_bytes()) {
            if line == b":" {
                if part == Part::Rules {
                    // No block is read, so the place in the range lines is not kept.
                    continue;
                }
                if hints
                    .special2_submap
                    .is_some_and(|fmt| !part.has_specials_of(fmt))
                {
                    ranges.skip_block();
                    continue;
                }
                while let Some(range_line) = ranges.next_line() {
                    guide.add_special2_ranges(range_line.split(|&c| c == b' '), &mut hints);
                }
                continue;
            }
            words.clear();
            words.extend(line.split(|&c| c == b' '));
            guide.add_rule(&words, &mut hints);
        }
        for specials in guide
            .rules
            .values_mut()
            .filter_map(|submap| submap.specials.get_mut())
        {
            let Specials {
                small_specials,
                small_options,
                ..
            } = specials;
            small_specials.retain(|_, slot| {
                let options = &mut small_options[*slot].1;
                options.retain(|(bits, _)| bits.any());
                !options.is_empty()
            });
        }
        if guide.general.is_empty() {
            for db in ["SRA", "TI", "TR_ASSM_CH", "TRACE_ASSM", "TRACE_CHGR"] {
                guide.general.insert(db.as_bytes().to_vec(), GENERAL_NUC);
            }
        }
        guide
    }
}

/// The lines of `kBuiltInGuide` that `AddRule` reads, except the range lines of `special2`
/// blocks; a line `:` stands for one block of them (`make_accguide.py`).
const RULES: &str = include_str!("rules.txt");

/// The lines of `text` as `str::lines` gives them (without the `\n` or `\r\n` that ends each).
fn rule_lines(text: &[u8]) -> impl Iterator<Item = &[u8]> {
    let mut rest = text;
    std::iter::from_fn(move || {
        if rest.is_empty() {
            return None;
        }
        let line = match find_byte(rest, b'\n') {
            Some(end) => {
                let line = &rest[..end];
                rest = &rest[end + 1..];
                line.strip_suffix(b"\r").unwrap_or(line)
            }
            None => std::mem::take(&mut rest),
        };
        Some(line)
    })
}

/// The range lines, four bits per character (`make_accguide.py`).
const RANGES: &[u8] = include_bytes!("ranges.bin");

/// The reader of `ranges.bin`: one line of a block at a time (its tokens joined by
/// spaces), `None` at the end of the block.
struct RangeLines<'a> {
    data: &'a [u8],
    nibble: usize,
    line: Vec<u8>,
}

impl<'a> RangeLines<'a> {
    fn new(data: &'a [u8]) -> Self {
        RangeLines {
            data,
            nibble: 0,
            line: Vec::new(),
        }
    }

    fn next_nibble(&mut self) -> u8 {
        let byte = self.data.get(self.nibble / 2).copied().unwrap_or(0xdd);
        let value = if self.nibble.is_multiple_of(2) {
            byte >> 4
        } else {
            byte & 15
        };
        self.nibble += 1;
        value
    }

    /// Passes over the rest of the block, up to the nibble that ends it for `next_line` (13 or
    /// more), a byte at a time.
    fn skip_block(&mut self) {
        if !self.nibble.is_multiple_of(2) && self.next_nibble() >= 13 {
            return;
        }
        let start = self.nibble / 2;
        let rest = self.data.get(start..).unwrap_or_default();
        self.nibble = match find_first(rest, |byte| byte >> 4 >= 13 || byte & 15 >= 13) {
            Some(offset) if rest[offset] >> 4 >= 13 => 2 * (start + offset) + 1,
            Some(offset) => 2 * (start + offset) + 2,
            None => 2 * self.data.len() + 1,
        };
    }

    fn next_line(&mut self) -> Option<&[u8]> {
        // The characters of the values below 12; 12 ends a line, 13 and more the block.
        const CHARS: [u8; 12] = *b"0123456789- ";
        self.line.clear();
        loop {
            let value = self.next_nibble();
            match CHARS.get(usize::from(value)) {
                Some(&c) => self.line.push(c),
                None => return (value == 12).then_some(self.line.as_slice()),
            }
        }
    }
}

/// The guide, read at its first use (`s_Guide`, a `CSafeStatic` made at the first call), with the
/// specials of each format read at the first lookup that reads them.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1627-1632
pub fn guide() -> &'static AccGuide {
    static GUIDE: OnceLock<AccGuide> = OnceLock::new();
    GUIDE.get_or_init(|| AccGuide::built_in(Part::Rules))
}

/// `CSeq_id::IdentifyAccession` of `acc` with `flags`: the type of a bare accession. `warn` gets the
/// diagnostic of an accession with a protein flag and an unexpected division.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1634-1663
/// (For an empty accession, `acc.size()` less 1 wraps to `NPOS`, which no `find` result passes.)
pub fn identify_accession(acc: &[u8], flags: u32, warn: &mut dyn FnMut(SeqIdWarning)) -> AccInfo {
    let (main_size, has_version) = match acc.iter().position(|&c| c == b'.') {
        None => (acc.len(), false),
        Some(dot) => {
            if dot + 1 >= acc.len() || !acc[dot + 1..].iter().all(u8::is_ascii_digit) {
                return UNKNOWN;
            }
            (dot, true)
        }
    };
    const MAIN_ACC_BUF_SIZE: usize = 32;
    if main_size <= MAIN_ACC_BUF_SIZE {
        let mut main_acc_buf = [0u8; MAIN_ACC_BUF_SIZE];
        let main_acc = &mut main_acc_buf[..main_size];
        main_acc.copy_from_slice(&acc[..main_size]);
        main_acc.make_ascii_uppercase();
        x_identify_accession(main_acc, flags, has_version, warn)
    } else {
        let main = acc[..main_size].to_ascii_uppercase();
        x_identify_accession(&main, flags, has_version, warn)
    }
}

/// `CSeq_id::x_IdentifyAccession` of `main_acc`, with `flags` and `has_version`.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1665-1772
/// In NCBI, cases 7, 6 and 5 of this `switch` are commented out, so only case 4 is live.
/// `strchr` also finds the terminating NUL (a NUL at index 4 is a delimiter). The value is
/// never a fallback (see the module), so the fallback branch (Seq_id.cpp:1770-1792) is not
/// ported.
pub fn x_identify_accession(
    main: &[u8],
    flags: u32,
    has_version: bool,
    warn: &mut dyn FnMut(SeqIdWarning),
) -> AccInfo {
    let alnum = |c: u8| c.is_ascii_alphanumeric();
    let alpha = |c: u8| c.is_ascii_alphabetic();
    let digit = |c: u8| c.is_ascii_digit();
    // `main_acc[i]` past the end is NUL (`CTempString::operator[]`).
    let at = |i: usize| main.get(i).copied().unwrap_or(0);
    let Some(digit_pos) = main.iter().position(|&c| digit(c)) else {
        return UNKNOWN;
    };
    let main_size = main.len();
    let mut flag_char = 0u8;
    let non_dig_pos = main[digit_pos..]
        .iter()
        .position(|&c| !digit(c))
        .map(|offset| digit_pos + offset);
    if let Some(non_dig_pos) = non_dig_pos {
        if flags & RAW_TEXT != 0 {
            if !has_version
                && digit_pos == 0
                && main_size >= 4
                && non_dig_pos < 5
                && alnum(at(1))
                && alnum(at(2))
                && alnum(at(3))
            {
                #[expect(
                    clippy::if_same_then_else,
                    reason = "NCBI returns eAcc_pdb both for a delimited ID and in its switch's case 4 (Seq_id.cpp:1683-1710)"
                )]
                if main_size > 4
                    && main_size <= 17
                    && matches!(at(4), b'|' | b'-' | b'_' | 0)
                    && (main_size <= 6 || alnum(at(5)))
                {
                    return acc(SEQID_PDB, OTHER, SEQ);
                } else if main_size == 4 {
                    return acc(SEQID_PDB, OTHER, SEQ);
                }
            }
            let oqp = matches!(at(0), b'O' | b'P' | b'Q');
            #[expect(
                clippy::if_same_then_else,
                reason = "NCBI has the two UniProt patterns as two branches that both return eAcc_swissprot (Seq_id.cpp:1713-1727)"
            )]
            if digit_pos == 1
                && main_size == 6
                && (oqp || alpha(at(2)))
                && digit(at(1))
                && alnum(at(2))
                && alnum(at(3))
                && alnum(at(4))
                && digit(at(5))
            {
                return acc(SEQID_SWISSPROT, OTHER, PROT);
            } else if digit_pos == 1
                && main_size == 10
                && !oqp
                && alpha(at(2))
                && alnum(at(3))
                && alnum(at(4))
                && digit(at(5))
                && alpha(at(6))
                && alnum(at(7))
                && alnum(at(8))
                && digit(at(9))
            {
                return acc(SEQID_SWISSPROT, OTHER, PROT);
            } else if !has_version
                && digit_pos == 0
                && (non_dig_pos == 6 || non_dig_pos == 7)
                && (main_size == non_dig_pos + 1
                    || at(non_dig_pos + 1) == b':'
                    || (alpha(at(non_dig_pos + 1))
                        && (main_size == non_dig_pos + 2 || at(non_dig_pos + 2) == b':')))
            {
                return acc(SEQID_PRF, OTHER, PROT);
            } else if digit_pos >= 4
                && non_dig_pos == digit_pos + 2
                && main_size - non_dig_pos >= 6
                && at(3) != b'_'
                && matches!(at(non_dig_pos), b'S' | b'P')
                && main[non_dig_pos + 1..].iter().all(|&c| digit(c))
            {
                flag_char = at(non_dig_pos);
            } else {
                return UNKNOWN;
            }
        }
    }

    if digit_pos == 0 {
        return if flags & RAW_GI != 0
            && !has_version
            && at(0) != b'0'
            && main.iter().all(|&c| digit(c))
        {
            GI
        } else {
            UNKNOWN
        };
    } else if flags & RAW_TEXT == 0 {
        return UNKNOWN;
    }

    let flag_len = usize::from(flag_char != 0);
    let digit_count = main_size - digit_pos - flag_len;
    let mut ai = guide().find(format_key(digit_pos, digit_count), main, None);
    ai &= !F_SPECIALS;
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1793-1829
    if flag_char == b'P' {
        match ai & DIVISION_MASK {
            TARGETED | WGS => ai = (ai & TYPE_MASK) | WGS | F_PROT,
            TSA => ai = (ai & TYPE_MASK) | TSA | F_PROT,
            division => {
                warn(SeqIdWarning {
                    subcode: 11,
                    message: [
                        main,
                        format!(
                            ": Protein flag found with unexpected division {}",
                            division >> 8
                        )
                        .as_bytes(),
                    ]
                    .concat(),
                });
                ai = (ai & (TYPE_MASK | DIVISION_MASK)) | F_PROT;
            }
        }
    }
    if matches!(ai & DIVISION_MASK, TARGETED | TSA | WGS | WGS_INTERMED)
        && digit_pos >= 4
        && main
            .get(digit_pos + 2..)
            .unwrap_or_default()
            .iter()
            .all(|&c| c == b'0')
    {
        return ai | F_MASTER;
    }
    ai
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1094-1108
/// `IsLocalId`: whether `id` is a local Seq-id, or one whose `IdentifyAccession` (the member
/// function's default flags) is `eAcc_local` or `eAcc_general`. The archive's and the search
/// strategy's request then holds the queries' Bioseqs (`CExportStrategy::x_Process_Query`).
pub fn is_local_id(id: &ParsedSeqId) -> bool {
    let info = identify_seq_id(
        id,
        parse_flags::ANY_RAW | parse_flags::ANY_LOCAL,
        &mut |_| {},
    );
    matches!(id, ParsedSeqId::Local(_)) || info == LOCAL || info == GENERAL
}

/// `CSeq_id::IdentifyAccession` with `flags`, on a Seq-id (the member function; its default
/// flags are `fParse_AnyRaw` and `fParse_AnyLocal` together).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1832-1876
pub fn identify_seq_id(
    id: &ParsedSeqId,
    flags: u32,
    warn: &mut dyn FnMut(SeqIdWarning),
) -> AccInfo {
    let seq_type = choice(id);
    match id {
        ParsedSeqId::Pir(_) | ParsedSeqId::Swissprot(_) | ParsedSeqId::Prf(_) => seq_type | F_PROT,
        ParsedSeqId::General(dbtag) => guide()
            .general
            .get(&dbtag.db.to_ascii_uppercase())
            .copied()
            .unwrap_or(GENERAL),
        _ => match id.textseq() {
            Some(text) => match &text.accession {
                Some(accession) => {
                    let ai = identify_accession(accession, flags, warn);
                    let type2 = acc_type(ai);
                    let div2 = ai & DIVISION_MASK;
                    if type2 == SEQID_NOT_SET {
                        (ai & FLAG_MASK) | seq_type
                    } else if type2 == seq_type {
                        ai
                    } else if seq_type == SEQID_TPE
                        && type2 == SEQID_EMBL
                        && (div2 == OTHER || div2 == WGS)
                    {
                        (ai & !TYPE_MASK) | seq_type
                    } else {
                        seq_type
                    }
                }
                None => seq_type,
            },
            None => seq_type,
        },
    }
}

/// `CSeq_id::Which()` of a Seq-id: its `E_Choice`.
pub fn choice(id: &ParsedSeqId) -> u32 {
    u32::from(id.alternative()) + 1
}

#[cfg(test)]
mod tests;
