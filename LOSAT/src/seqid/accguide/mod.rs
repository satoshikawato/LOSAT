//! NCBI's accession guide (`SAccGuide` of `c++/src/objects/seqloc/Seq_id.cpp`) and
//! `CSeq_id::IdentifyAccession`: the type of a bare accession such as `XP_011520000.1` (a
//! RefSeq protein) or `P69905` (Swiss-Prot), as `makeblastdb -parse_seqids` reads it.
//!
//! The rules are NCBI's built-in guide (`kBuiltInGuide`, `accguide2.inc` of 2.17.0, which
//! makeblastdb uses: no `accguide2.txt` data file is looked for, so the result does not depend
//! on the machine). They are kept in two files beside this module, made by
//! `LOSAT/tests/accguide/make_accguide.py`: `rules.txt` (every line that `AddRule` reads,
//! without its comment) and `ranges.bin` (the range lines of `special2` blocks, four bits per
//! character). They are read into the guide once, at the first lookup.
//!
//! What is not ported: the branches of `AddRule` for a type name that this build does not
//! know (the fallbacks: every name of the built-in guide is known, which the tests check),
//! the diagnostics of the loading (`ERR_POST_X` 2-21: the built-in guide raises none of
//! Warning severity, also checked; the Info ones are not shown by makeblastdb), and the
//! fallback warning of `IdentifyAccession` (no value of the built-in guide is a fallback).

mod acc_info;

use std::collections::BTreeMap;
use std::sync::OnceLock;

pub use acc_info::*;

use super::model::ParsedSeqId;

/// A `CSeq_id::EAccessionInfo` value.
pub type AccInfo = u32;

/// `CSeq_id::GetAccType`: the Seq-id type of an accession's value (an `E_Choice`, 0 for
/// `e_not_set`).
///
/// NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_id.hpp:562-563
/// ```c++
///     static E_Choice GetAccType(EAccessionInfo info)
///         { return static_cast<E_Choice>(info & eAcc_type_mask); }
/// ```
pub fn acc_type(info: AccInfo) -> u32 {
    info & TYPE_MASK
}

/// The parse flags of `CSeq_id` (`EParseFlags`) that the accession functions read.
///
/// NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_id.hpp:75-103
/// ```c++
///     enum EParseFlags {
///         fParse_PartialOK  = 0x01,
///         fParse_RawText    = 0x02, ///< Try to ID raw non-numeric accessions
///         fParse_RawGI      = 0x04, ///< Treat raw numbers as GIs, not local IDs
///         fParse_AnyRaw     = fParse_RawText | fParse_RawGI,
///         ...
///         fParse_ValidLocal = 0x08,
///         ...
///         fParse_AnyLocal   = 0x18,
///         fParse_NoFASTA    = 0x20, ///< Don't bother checking for a tag
///         ...
///         fParse_FallbackOK = 0x40,
///         ...
///         fParse_Default    = fParse_RawText | fParse_ValidLocal
///     };
/// ```
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

/// A diagnostic that NCBI's Seq-id code posts (`ERR_POST_X(subcode, Warning << message)`
/// of Seq_id.cpp, whose error code is `NCBI_ERRCODE_X` 1306): the caller writes it as
/// `Warning: (1306.<subcode>) [<application>] <message>`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeqIdWarning {
    pub subcode: u32,
    pub message: String,
}

/// A set of bit positions as `bm::bvector<>` holds them for a small special: sorted,
/// disjoint, non-adjacent closed ranges (the guide's bit vectors are up to 10^9 bits long
/// and mostly empty).
#[derive(Default, Debug)]
struct BitRanges(Vec<(u32, u32)>);

impl BitRanges {
    /// `set_range(left, right)`.
    fn set_range(&mut self, left: u32, right: u32) {
        let ranges = &mut self.0;
        // The ranges that touch or overlap [left, right].
        let first = ranges.partition_point(|&(_, end)| end.saturating_add(1) < left);
        let last = ranges.partition_point(|&(start, _)| start <= right.saturating_add(1));
        if first == last {
            ranges.insert(first, (left, right));
        } else {
            let start = left.min(ranges[first].0);
            let end = right.max(ranges[last - 1].1);
            ranges.splice(first..last, [(start, end)]);
        }
    }

    /// `clear_range(left, right)`.
    fn clear_range(&mut self, left: u32, right: u32) {
        let ranges = &mut self.0;
        let first = ranges.partition_point(|&(_, end)| end < left);
        let last = ranges.partition_point(|&(start, _)| start <= right);
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
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:891-911
/// ```c++
///     typedef map<string, TAccInfo>   TPrefixes;
///     typedef pair<string, TAccInfo>  TPair;
///     typedef list<TPair>             TPairs; // not vector -- need stable ptrs
///     typedef map<string, TPair>      TBigSpecialMap; // last -> first -> value
///     typedef pair<bm::bvector<>, TAccInfo> TSmallSpecialOption;
///     typedef multimap<string, TSmallSpecialOption> TSmallSpecialMap;
///     ...
///     struct SSubMap {
///         TPrefixes         prefixes;
///         TPairs            wildcards;
///         TBigSpecialMap    big_specials;
///         TSmallSpecialMap  small_specials;
///     };
///     typedef map<TFormatCode, SSubMap> TMainMap;
/// ```
/// The multimap is a map of each key to its options in the multimap's order (an
/// `emplace` adds an option at the end of its key's equal range).
#[derive(Default, Debug)]
struct SubMap {
    prefixes: BTreeMap<Vec<u8>, AccInfo>,
    wildcards: Vec<(Vec<u8>, AccInfo)>,
    big_specials: BTreeMap<Vec<u8>, (Vec<u8>, AccInfo)>,
    small_specials: BTreeMap<Vec<u8>, Vec<(BitRanges, AccInfo)>>,
}

/// `SAccGuide`.
#[derive(Debug)]
pub struct AccGuide {
    rules: BTreeMap<u32, SubMap>,
    general: BTreeMap<Vec<u8>, AccInfo>,
    /// The diagnostics of Warning severity that NCBI's loading would post (none for the
    /// built-in guide).
    load_warnings: Vec<String>,
}

/// `SAccGuide::SHints`: the state that `AddRule` keeps from one line to the next.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:913-945
/// ```c++
///     struct SHints {
///         SHints()
///             : prev_type(CSeq_id::eAcc_unknown), prev_submap(NULL),
///               prev_special_format(0),
///               prev_special_type(CSeq_id::eAcc_unknown),
///               prev_special_base_type(CSeq_id::eAcc_unknown),
///               special2_submap(nullptr), version(1)
///             {}
///         ...
///         TAccInfo              prev_type;
///         CTempString           prev_type_name;
///         TMainMap::value_type* prev_submap;
///         TBigSpecialMap::iterator   prev_big_special;
///         TSmallSpecialMap::iterator prev_small_special;
///         TFormatCode           prev_special_format;
///         string                prev_special_key;
///         CTempString           prev_special_type_name;
///         string                prev_special_base_key;
///         TAccInfo              prev_special_type;
///         TAccInfo              prev_special_base_type;
///         string                prev_special2_acc;
///         map<string, CTempString> default_fallbacks;
///         string                special2_name;
///         unique_ptr<string>    special2_old_name;
///         SSubMap*              special2_submap;
///         TAccInfo              special2_type;
///         unsigned int          version;
///     };
/// ```
/// `prev_small_special` is the key and the position of the option in its equal range
/// (`None` for `end()`); `prev_big_special` is only an insertion hint and is not kept.
/// The fallback members serve unknown type names only.
struct Hints {
    prev_type: AccInfo,
    prev_type_name: Vec<u8>,
    prev_submap: Option<u32>,
    prev_small_special: Option<(Vec<u8>, usize)>,
    prev_special_format: u32,
    prev_special_key: Vec<u8>,
    prev_special_type_name: Vec<u8>,
    prev_special_base_key: Vec<u8>,
    prev_special_type: AccInfo,
    prev_special_base_type: AccInfo,
    prev_special2_acc: Vec<u8>,
    special2_submap: Option<u32>,
    special2_type: AccInfo,
    version: u32,
}

impl Default for Hints {
    fn default() -> Self {
        Hints {
            prev_type: UNKNOWN,
            prev_type_name: Vec::new(),
            prev_submap: None,
            prev_small_special: None,
            prev_special_format: 0,
            prev_special_key: Vec::new(),
            prev_special_type_name: Vec::new(),
            prev_special_base_key: Vec::new(),
            prev_special_type: UNKNOWN,
            prev_special_base_type: UNKNOWN,
            prev_special2_acc: Vec::new(),
            special2_submap: None,
            special2_type: UNKNOWN,
            version: 1,
        }
    }
}

/// `kMaxSmallSpecialDigits` and `kBVSizes`.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:866-884
/// ```c++
/// static const unsigned int kMaxSmallSpecialDigits = 9;
/// static const bm::bvector<>::size_type kBVSizes[kMaxSmallSpecialDigits + 1] = {
///     1,
///     10,
///     ...
///     1000000000
/// };
/// ```
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

/// `SAccGuide::s_Key(letters, digits)`.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:962-963
/// ```c++
///     static TFormatCode s_Key(unsigned short letters, unsigned short digits)
///         { return TFormatCode(letters) << 16 | digits; }
/// ```
fn format_key(letters: usize, digits: usize) -> u32 {
    ((letters as u32 & 0xffff) << 16) | (digits as u32 & 0xffff)
}

/// `NStr::StringToNumeric(str, &value)` of a bit position: 0 when it is not a number that
/// fits.
fn bit_number(digits: &[u8]) -> u32 {
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        return 0;
    }
    std::str::from_utf8(digits)
        .ok()
        .and_then(|text| text.parse().ok())
        .unwrap_or(0)
}

/// `SAccGuide::x_SplitSpecial(acc, fmt)`: the length of the prefix that keys the special
/// and the number that is its bit.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1599-1629
/// ```c++
/// bm::bvector_size_type SAccGuide::x_SplitSpecial(CTempString& acc,
///                                                 TFormatCode fmt)
/// {
///     auto raw_digits = fmt & 0xffff, digits = raw_digits;
///     auto normal_size = (fmt >> 16) + digits;
///     if (digits == kMaxSmallSpecialDigits + 1) {
///         digits -= 2;
///     } else if (digits > kMaxSmallSpecialDigits) {
///         digits = kMaxSmallSpecialDigits;
///     }
///     SIZE_TYPE pos;
///     bm::bvector_size_type result;
///     if (acc.size() == normal_size) {
///         pos = acc.size() - digits;
///         NStr::StringToNumeric(acc.substr(pos), &result);
///     } else {
///         _ASSERT(acc.size() == normal_size + 1);
///         _ASSERT(digits >= 3);
///         pos = (fmt >> 16) + 2;
///         _ASSERT(isalpha(static_cast<unsigned char>(acc[pos])));
///         NStr::StringToNumeric(acc.substr(pos + 1), &result);
///         if (digits == raw_digits) {
///             pos -= 2;
///             result += (NStr::StringToNumeric<Uint1>(acc.substr(pos, 2))
///                        * kBVSizes[digits - 2]);
///         }
///     }
///     acc.erase(pos);
///     return result;
/// }
/// ```
fn split_special(acc: &[u8], fmt: u32) -> (usize, u32) {
    let raw_digits = fmt & 0xffff;
    let mut digits = raw_digits;
    let normal_size = (fmt >> 16) as usize + digits as usize;
    if digits == MAX_SMALL_SPECIAL_DIGITS + 1 {
        digits -= 2;
    } else if digits > MAX_SMALL_SPECIAL_DIGITS {
        digits = MAX_SMALL_SPECIAL_DIGITS;
    }
    if acc.len() == normal_size {
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

/// `NStr::MatchesMask(str, mask)` (case-sensitive) for the masks of the guide, which are
/// letters, `_` and `?` (any one character); a NUL character of `str` matches nothing.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:291-307,359-386
/// ```c++
/// static EMatchesMaskResult s_MatchesMask(CTempString str, CTempString mask, bool ignore_case)
/// {
///     char s, m;
///     size_t str_pos = 0, mask_pos = 0;
///
///     for ( ; (m = mask[mask_pos]); ++str_pos, ++mask_pos) {
///
///         s = str[str_pos];
///
///         if (!s  &&  m != '*') {
///             return eMismatch;
///         }
///         // Analyze mask symbol
///         switch ( m ) {
///         case '?':
///             _ASSERT(s);
///             continue;
///         ...
///         default:
///             // Compare non pattern character in mask and name
///             _ASSERT(s  &&  m);
///             if (ignore_case) {
///                 ...
///             } else {
///                 if (s != m)
///                     return eNoMatch;
///             }
///             continue;
///         }
///     }
///     // Matches if we reach the end of the string and mask at the same time only
///     if ( str[str_pos] ) {
///         return eNoMatch;
///     }
///     return eMatch;
/// }
/// ```
/// (`CTempString::operator[]` gives NUL past the end.) The loading refuses a mask with `*`,
/// `[` or `\`.
fn matches_mask(text: &[u8], mask: &[u8]) -> bool {
    for (pos, &m) in mask.iter().enumerate() {
        if m == 0 {
            break;
        }
        let s = text.get(pos).copied().unwrap_or(0);
        if s == 0 {
            return false;
        }
        if m != b'?' && s != m {
            return false;
        }
    }
    let end = mask.iter().position(|&m| m == 0).unwrap_or(mask.len());
    text.get(end).copied().unwrap_or(0) == 0
}

impl AccGuide {
    /// `SAccGuide::Find(fmt, acc_or_pfx, key_used)`: the value of an accession (or of a
    /// prefix) of the format `fmt`. `key_used` gets the wildcard or the prefix that gave
    /// it, as NCBI's (for the loading's diagnostics).
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1370-1441
    /// ```c++
    /// const SAccGuide::TAccInfo& SAccGuide::Find(TFormatCode fmt,
    ///                                            const CTempString& acc_or_pfx,
    ///                                            string* key_used) const
    /// {
    ///     static const TAccInfo kUnknown = CSeq_id::eAcc_unknown;
    ///     TMainMap::const_iterator it = rules.find(fmt);
    ///     if (it == rules.end()) {
    ///         return kUnknown;
    ///     }
    ///
    ///     const SSubMap&            submap = it->second;
    ///     const TAccInfo*           result = &kUnknown;
    ///     CTempString               pfx     (acc_or_pfx, 0, fmt >> 16);
    ///     TPrefixes::const_iterator pit    = submap.prefixes.find(pfx);
    ///     if (pit != submap.prefixes.end()) {
    ///         result = &pit->second;
    ///     } else {
    ///         ITERATE (TPairs, wit, submap.wildcards) {
    ///             if (NStr::MatchesMask(pfx, wit->first)) {
    ///                 bool bad_match = false; // Limit ? to matching letters
    ///                 SIZE_TYPE pos = wit->first.find('?');
    ///                 while (pos != NPOS) {
    ///                     if ( !isalnum(pfx[pos])  &&  pfx[pos] != '?' ) {
    ///                         bad_match = true;
    ///                         break;
    ///                     } else {
    ///                         pos = wit->first.find('?', pos + 1);
    ///                     }
    ///                 }
    ///                 if (bad_match) {
    ///                     continue;
    ///                 }
    ///                 if (key_used  &&  acc_or_pfx != wit->first) {
    ///                     *key_used = wit->first;
    ///                 }
    ///                 result = &wit->second;
    ///                 break;
    ///             }
    ///         }
    ///     }
    ///     if (acc_or_pfx != pfx  &&  (*result & CSeq_id::fAcc_specials) != 0) {
    ///         pfx = acc_or_pfx;
    ///         auto n = x_SplitSpecial(pfx, fmt);
    ///         for (auto ssit = submap.small_specials.lower_bound(pfx);
    ///              ssit != submap.small_specials.end()  &&  ssit->first == pfx;
    ///              ++ssit) {
    ///             if (ssit->second.first[n]) {
    ///                 if (key_used) {
    ///                     key_used->erase();
    ///                 }
    ///                 return ssit->second.second;
    ///             }
    ///         }
    ///         TBigSpecialMap::const_iterator bsit
    ///             = submap.big_specials.lower_bound(acc_or_pfx);
    ///         if (bsit != submap.big_specials.end()
    ///             &&  !(acc_or_pfx < bsit->second.first) ) {
    ///             if (key_used) {
    ///                 key_used->erase();
    ///             }
    ///             return bsit->second.second;
    ///         } else {
    ///             if (key_used  &&  key_used->empty()) {
    ///                 *key_used = pfx.substr(0, fmt >> 16);
    ///             }
    ///             return *result;
    ///         }
    ///     } else /* if (*result != CSeq_id::eAcc_unknown) */ {
    ///         return *result;
    ///     }
    /// }
    /// ```
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
            let (prefix_len, n) = split_special(acc_or_pfx, fmt);
            let prefix = &acc_or_pfx[..prefix_len];
            if let Some(options) = submap.small_specials.get(prefix) {
                for (bits, value) in options {
                    if bits.test(n) {
                        if let Some(key_used) = key_used {
                            key_used.clear();
                        }
                        return *value;
                    }
                }
            }
            if let Some((_, (first, value))) =
                submap.big_specials.range(acc_or_pfx.to_vec()..).next()
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

    /// `SHints::FindAccInfo(name)`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:981-996
    /// ```c++
    /// inline
    /// SAccGuide::TAccInfo SAccGuide::SHints::FindAccInfo(CTempString name)
    /// {
    ///     if (name == prev_type_name) {
    ///         return prev_type;
    ///     } else {
    ///         TAccInfoMap::const_iterator it = sc_AccInfoMap.find(name);
    ///         if (it == sc_AccInfoMap.end()) {
    ///             return kUnrecognized;
    ///         } else {
    ///             prev_special_key.clear();
    ///             prev_type_name = it->first;
    ///             return prev_type = it->second;
    ///         }
    ///     }
    /// }
    /// ```
    /// `None` is `kUnrecognized`.
    fn find_acc_info_hinted(hints: &mut Hints, name: &[u8]) -> Option<AccInfo> {
        if !hints.prev_type_name.is_empty() && name == hints.prev_type_name.as_slice() {
            return Some(hints.prev_type);
        }
        let value = find_acc_info(name)?;
        hints.prev_special_key.clear();
        hints.prev_type_name = name.to_ascii_lowercase();
        hints.prev_type = value;
        Some(value)
    }

    /// `SHints::FindSpecial(guide, fmt, acc_or_range)`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:998-1010
    /// ```c++
    /// inline
    /// SAccGuide::TAccInfo SAccGuide::SHints::FindSpecial(const SAccGuide& guide,
    ///                                                    TFormatCode fmt,
    ///                                                    CTempString acc_or_range)
    /// {
    ///     CTempString pfx(acc_or_range, 0, fmt >> 16);
    ///     if (fmt == prev_special_format) {
    ///         if (acc_or_range == prev_special_key) {
    ///             prev_special_type = prev_type;
    ///             prev_special_type_name = prev_type_name;
    ///             return prev_type;
    ///         } else if (pfx == prev_special_base_key) {
    ///             return prev_special_base_type;
    ///         }
    ///     }
    ///     prev_special_format    = fmt;
    ///     prev_special_base_key  = pfx;
    ///     prev_special_base_type = guide.Find(fmt, pfx);
    ///     return prev_special_base_type;
    /// }
    /// ```
    fn find_special(&self, hints: &mut Hints, fmt: u32, acc_or_range: &[u8]) -> AccInfo {
        let pfx = &acc_or_range[..((fmt >> 16) as usize).min(acc_or_range.len())];
        if fmt == hints.prev_special_format {
            if acc_or_range == hints.prev_special_key.as_slice() {
                hints.prev_special_type = hints.prev_type;
                hints.prev_special_type_name = hints.prev_type_name.clone();
                return hints.prev_type;
            } else if pfx == hints.prev_special_base_key.as_slice() {
                return hints.prev_special_base_type;
            }
        }
        hints.prev_special_format = fmt;
        hints.prev_special_base_key = pfx.to_vec();
        hints.prev_special_base_type = self.find(fmt, pfx, None);
        hints.prev_special_base_type
    }

    /// `SHints::FindSubMap(rules, fmt)`: the rules of `fmt`, made if absent.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1012-1028
    /// ```c++
    /// inline
    /// SAccGuide::SSubMap& SAccGuide::SHints::FindSubMap(SAccGuide::TMainMap& rules,
    ///                                                   SAccGuide::TFormatCode fmt)
    /// {
    ///     if (prev_submap != NULL  &&  prev_submap->first == fmt) {
    ///         return prev_submap->second;
    ///     } else {
    ///         SAccGuide::TMainMap::iterator it = rules.lower_bound(fmt);
    ///         if (it == rules.end() || it->first != fmt) {
    ///             it = rules.insert(it, make_pair(fmt, SAccGuide::SSubMap()));
    ///         }
    ///         prev_submap = &*it;
    ///         prev_big_special   = it->second.big_specials.end();
    ///         prev_small_special = it->second.small_specials.end();
    ///         return it->second;
    ///     }
    /// }
    /// ```
    fn find_sub_map(&mut self, hints: &mut Hints, fmt: u32) -> u32 {
        if hints.prev_submap != Some(fmt) {
            self.rules.entry(fmt).or_default();
            hints.prev_submap = Some(fmt);
            hints.prev_small_special = None;
        }
        fmt
    }

    /// `SAccGuide::x_AddSpecial(submap, hints, fmt, from, to, value, old_name, new_name)`
    /// for a known type name (`old_name` null).
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1513-1597
    /// ```c++
    /// void SAccGuide::x_AddSpecial(SSubMap& submap, SHints& hints, TFormatCode fmt,
    ///                              CTempString from, CTempString to, TAccInfo value,
    ///                              const string* old_name,
    ///                              const CTempString& new_name)
    /// {
    ///     CTempString from_pfx = from, to_pfx = to;
    ///     auto left  = x_SplitSpecial(from_pfx, fmt),
    ///          right = x_SplitSpecial(to_pfx, fmt);
    ///     const TAccInfo* value_ptr = nullptr;
    ///     if (from_pfx != to_pfx) {
    ///         hints.prev_big_special
    ///             = submap.big_specials.insert(hints.prev_big_special,
    ///                                          make_pair(to, TPair(from, value)));
    ///         // Account for possible refinement.
    ///         hints.prev_big_special->second.second = value;
    ///         ...
    ///     } else {
    ///         TSmallSpecialMap::iterator it = submap.small_specials.end();
    ///         if (hints.prev_small_special != submap.small_specials.end()
    ///             &&  hints.prev_small_special->first == from_pfx) {
    ///             it = hints.prev_small_special;
    ///             it->second.first.clear_range(left, right);
    ///             while ((it->second.second & ~CSeq_id::fAcc_fallback)
    ///                    != (value & ~CSeq_id::fAcc_fallback)) {
    ///                 if (it == submap.small_specials.begin()
    ///                     ||  (--it)->first != from_pfx) {
    ///                     it = hints.prev_small_special;
    ///                     ++it;
    ///                     break;
    ///                 }
    ///             }
    ///         } else {
    ///             it = submap.small_specials.lower_bound(from_pfx);
    ///         }
    ///         while (it != submap.small_specials.end()) {
    ///             if (it->first != from_pfx) {
    ///                 it = submap.small_specials.end();
    ///                 break;
    ///             } else if ((it->second.second & ~CSeq_id::fAcc_fallback)
    ///                        == (value & ~CSeq_id::fAcc_fallback)) {
    ///                 break;
    ///             } else {
    ///                 ++it;
    ///             }
    ///         }
    ///         if (it != submap.small_specials.end()) {
    ///             ...
    ///             hints.prev_small_special = it;
    ///         } else {
    ///             auto size = kBVSizes[min(fmt & 0xffff, kMaxSmallSpecialDigits)];
    ///             hints.prev_small_special =
    ///                 submap.small_specials.emplace(
    ///                     from_pfx, make_pair(bm::bvector<>(size), value));
    ///         }
    ///         hints.prev_small_special->second.first.set_range(left, right);
    ///         // Account for possible refinement.
    ///         hints.prev_small_special->second.second = value;
    ///         ...
    ///     }
    ///     ...
    /// }
    /// ```
    /// A bit vector needs no size here (`BitRanges`); `x_SplitSpecial` never gives a bit
    /// past it.
    fn add_special(&mut self, hints: &mut Hints, fmt: u32, from: &[u8], to: &[u8], value: AccInfo) {
        let (from_len, left) = split_special(from, fmt);
        let (to_len, right) = split_special(to, fmt);
        let submap = self.rules.get_mut(&fmt).expect("FindSubMap made the rules");
        let (from_pfx, to_pfx) = (&from[..from_len], &to[..to_len]);
        if from_pfx != to_pfx {
            let entry = submap
                .big_specials
                .entry(to.to_vec())
                .or_insert_with(|| (from.to_vec(), value));
            entry.1 = value;
            return;
        }
        let same = |option: &(BitRanges, AccInfo)| option.1 & !F_FALLBACK == value & !F_FALLBACK;
        let options = submap.small_specials.entry(from_pfx.to_vec()).or_default();
        let mut index = match &hints.prev_small_special {
            Some((key, prev)) if key.as_slice() == from_pfx => {
                let prev = *prev;
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
            _ => 0,
        };
        while index < options.len() && !same(&options[index]) {
            index += 1;
        }
        if index == options.len() {
            options.push((BitRanges::default(), value));
        }
        options[index].0.set_range(left, right);
        options[index].1 = value;
        hints.prev_small_special = Some((from_pfx.to_vec(), index));
    }

    /// `SAccGuide::AddRule(rule, hints)` for the lines of `rules.txt`; `ranges` gives the
    /// tokens of each `:` line of a `special2` block.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1030-1051
    /// ```c++
    /// void SAccGuide::AddRule(const CTempString& rule, SHints& hints)
    /// {
    ///     CTempString         tmp1, tmp2;
    ///     vector<CTempStringEx> tokens;
    ///     SIZE_TYPE           pos, pos2;
    ///
    ///     ++count;
    ///     tmp1.assign(rule, 0, rule.find('#')); // strip comment
    ///     if (tmp1.empty())
    ///         return;
    ///     tokens.reserve(3);
    ///     NStr::Split(tmp1, " \t", tokens,
    ///                 NStr::fSplit_MergeDelimiters | NStr::fSplit_Truncate);
    ///     if (tokens.empty()) {
    ///         return;
    ///     } else if (tokens.size() == 2
    ///                &&  NStr::EqualNocase(tokens[0], "version")) {
    ///         hints.special2_submap = nullptr;
    ///         hints.version = NStr::StringToUInt(tokens[1], NStr::fConvErr_NoThrow);
    ///         if (hints.version > 2  ||  hints.version < 1) {
    ///             ERR_POST_X(2, "SAccGuide::AddRule: " << count
    ///                           << ": Unsupported version " << tokens[1]);
    ///             return;
    ///         }
    /// ```
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
            self.add_special2_ranges(&tokens[1..], hints);
        } else if tokens.len() == 3 && first.eq_ignore_ascii_case(b"gnl") {
            // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1320-1328,1345-1347
            // ```c++
            //     } else if (tokens.size() == 3 && NStr::EqualNocase(tokens[0], "gnl")) {
            //         hints.special2_submap = nullptr;
            //         string key(tokens[1]);
            //         NStr::ToUpper(key);
            //         TAccInfo value = hints.FindAccInfo(tokens[2]);
            //         if (value == kUnrecognized) {
            //         ...
            //         } else {
            //             general[key] = value;
            //         }
            // ```
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
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1052-1064,1102-1131
    /// ```c++
    ///     } else if ((pos = tokens[0].find('+')) != NPOS
    ///                &&  (tokens.size() == 3
    ///                     ||  (tokens.size() == 4  &&  tokens[3] == "*"))) {
    ///         hints.special2_submap = nullptr;
    ///         // _VERIFY(NStr::SplitInTwo(tokens[0], "+", tmp1, tmp2));
    ///         tmp1.assign(tokens[0], 0, pos);
    ///         tmp2.assign(tokens[0], pos + 1, NPOS);
    ///         TFormatCode fmt
    ///             = s_Key(NStr::StringToUInt(tmp1, NStr::fConvErr_NoThrow),
    ///                     NStr::StringToUInt(tmp2, NStr::fConvErr_NoThrow));
    ///         TAccInfo value = hints.FindAccInfo(tokens[2]);
    ///         unique_ptr<string> old_name;
    ///         if (value == kUnrecognized) {
    ///         ...
    ///         if (value != kUnrecognized) {
    ///             SSubMap& submap = hints.FindSubMap(rules, fmt);
    ///             if (tokens.size() == 4) {
    ///                 value = TAccInfo(value | CSeq_id::fAcc_specials);
    ///             }
    ///             const TAccInfo* value_ptr = NULL;
    ///             if (tokens[1].find_first_of("?*") == NPOS) {
    ///                 value_ptr = &(submap.prefixes[tokens[1]] = value);
    ///             } else {
    ///                 // Account for possible refinements of fallback definitions
    ///                 NON_CONST_ITERATE (TPairs, wit, submap.wildcards) {
    ///                     if (wit->first == tokens[1]) {
    ///                         wit->second = value;
    ///                         value_ptr   = &wit->second;
    ///                         break;
    ///                     }
    ///                 }
    ///                 if (value_ptr == NULL) {
    ///                     submap.wildcards.push_back(TPair(tokens[1], value));
    ///                     value_ptr = &submap.wildcards.back().second;
    ///                 }
    ///             }
    /// ```
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
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1132-1149,1193-1210
    /// ```c++
    ///     } else if (tokens.size() == 3 && NStr::EqualNocase(tokens[0], "special")) {
    ///         hints.special2_submap = nullptr;
    ///         pos  = tokens[1].find_first_of(kDigits);
    ///         pos2 = tokens[1].find('-', pos);
    ///         TFormatCode fmt
    ///             = s_Key(pos, ((pos2 == NPOS) ? tokens[1].size() : pos2) - pos);
    ///         TAccInfo old   = hints.FindSpecial(*this, fmt, tokens[1]);
    ///         TAccInfo value = hints.FindAccInfo(tokens[2]);
    ///         if ((old & CSeq_id::fAcc_specials) != 0) {
    ///             old = TAccInfo(old & ~CSeq_id::fAcc_specials);
    ///         } else {
    ///             string key_used;
    ///             Find(fmt, tokens[1].substr(0, pos2), &key_used);
    ///             if ( !key_used.empty() ) {
    ///                 ERR_POST_X(13, Warning
    ///                            << "SAccGuide::AddRule: Main listing for special "
    ///                            << tokens[1]
    ///                            << " doesn't indicate that specials are present.");
    ///             }
    ///         }
    ///         unique_ptr<string> old_name;
    ///         if (value == kUnrecognized) {
    ///         ...
    ///         } else {
    ///             _ASSERT(hints.prev_type == value);
    ///             hints.prev_special_key = tokens[1];
    ///         }
    ///         if (value != kUnrecognized) {
    ///             SSubMap& submap = hints.FindSubMap(rules, fmt);
    ///             if (pos2 == NPOS) {
    ///                 tmp1 = tmp2 = tokens[1];
    ///             } else {
    ///                 tmp1.assign(tokens[1], 0, pos2);
    ///                 tmp2.assign(tokens[1], pos2 + 1, NPOS);
    ///             }
    ///             x_AddSpecial(submap, hints, fmt, tmp1, tmp2, value, old_name.get(),
    ///                          tokens[2]);
    ///         }
    ///         hints.prev_special_type_name.clear();
    ///         hints.prev_special_type = CSeq_id::eAcc_unknown;
    /// ```
    fn add_special_rule(&mut self, tokens: &[&[u8]], hints: &mut Hints) {
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
        let old = self.find_special(hints, fmt, range);
        let value = Self::find_acc_info_hinted(hints, tokens[2]);
        let first = &range[..pos2.unwrap_or(range.len())];
        if old & F_SPECIALS == 0 {
            let mut key_used = Vec::new();
            self.find(fmt, first, Some(&mut key_used));
            if !key_used.is_empty() {
                self.load_warnings.push(format!(
                    "Main listing for special {} doesn't indicate that specials are present.",
                    String::from_utf8_lossy(range)
                ));
            }
        }
        match value {
            None => self.unknown_name(tokens[2]),
            Some(value) => {
                hints.prev_special_key = range.to_vec();
                self.find_sub_map(hints, fmt);
                let last = pos2.map_or(range, |pos2| &range[pos2 + 1..]);
                self.add_special(hints, fmt, first, last, value);
            }
        }
        hints.prev_special_type_name.clear();
        hints.prev_special_type = UNKNOWN;
    }

    /// The `special2 <letters>+<digits> <type>` line of `AddRule`, which starts a block of
    /// range lines.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1213-1247
    /// ```c++
    ///     } else if (tokens.size() >= 3
    ///                &&  NStr::EqualNocase(tokens[0], "special2")) {
    ///         if (hints.version < 2) {
    ///             ERR_POST_X(18,
    ///                        Warning << "SAccGuide::AddRule: " << count
    ///                        << ": special2 valid only in version 2+ guides");
    ///         }
    ///         NStr::SplitInTwo(tokens[1], "+", tmp1, tmp2);
    ///         auto digits = NStr::StringToNumeric<unsigned short>(tmp2);
    ///         hints.prev_special_format = s_Key(tmp1.size(), digits);
    ///         hints.special2_name = tokens[2];
    ///         hints.special2_old_name.reset();
    ///         hints.special2_submap
    ///             = &hints.FindSubMap(rules, hints.prev_special_format);
    ///         hints.prev_special2_acc = tmp1 + string(digits, '0');
    ///         TAccInfo old = hints.FindSpecial(*this, hints.prev_special_format,
    ///                                          hints.prev_special2_acc);
    ///         string why;
    ///         if ((old & CSeq_id::fAcc_specials) != 0) {
    ///             old = TAccInfo(old & ~CSeq_id::fAcc_specials);
    ///         } else {
    ///             Find(hints.prev_special_format, hints.prev_special2_acc, &why);
    ///             if ( !why.empty() ) {
    ///                 ERR_POST_X(13, Warning
    ///                            << "SAccGuide::AddRule: Main listing for special "
    ///                            << tokens[1]
    ///                            << " doesn't indicate that specials are present.");
    ///             }
    ///         }
    ///         for (size_t i = 2;  i < tokens.size();  ++i) {
    ///             hints.special2_type = hints.FindAccInfo(tokens[i]);
    ///             if (hints.special2_type != kUnrecognized) {
    ///                 if (i > 2) {
    ///                     hints.special2_old_name.reset(new string(tokens[i]));
    ///                 }
    ///                 break;
    ///             }
    ///         }
    /// ```
    /// A known name in the third token is the only case without a fallback (the elided
    /// lines, 1248-1294).
    fn add_special2_rule(&mut self, tokens: &[&[u8]], hints: &mut Hints) {
        if hints.version < 2 {
            self.load_warnings
                .push("special2 valid only in version 2+ guides".to_string());
        }
        let spec = tokens[1];
        let plus = spec.iter().position(|&c| c == b'+').unwrap_or(spec.len());
        let letters = &spec[..plus];
        let digits = bit_number(spec.get(plus + 1..).unwrap_or_default()) as usize;
        let fmt = format_key(letters.len(), digits);
        hints.prev_special_format = fmt;
        hints.special2_submap = Some(self.find_sub_map(hints, fmt));
        let mut acc = letters.to_vec();
        acc.resize(letters.len() + digits, b'0');
        hints.prev_special2_acc = acc.clone();
        let old = self.find_special(hints, fmt, &acc);
        if old & F_SPECIALS == 0 {
            let mut why = Vec::new();
            self.find(fmt, &acc, Some(&mut why));
            if !why.is_empty() {
                self.load_warnings.push(format!(
                    "Main listing for special {} doesn't indicate that specials are present.",
                    String::from_utf8_lossy(spec)
                ));
            }
        }
        match Self::find_acc_info_hinted(hints, tokens[2]) {
            Some(value) => hints.special2_type = value,
            None => self.unknown_name(tokens[2]),
        }
        hints.prev_special_type_name.clear();
        hints.prev_special_type = UNKNOWN;
    }

    /// The `: <range> ...` lines of a `special2` block: each range replaces the last digits
    /// of the block's running accession.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1297-1319
    /// ```c++
    ///     } else if (tokens.size() >= 2  &&  tokens[0] == ":") {
    ///         ...
    ///         if (hints.special2_submap == nullptr) {
    ///             ERR_POST_X(20,
    ///                        Warning <<
    ///                        "SAccGuide::AddRule: " << count
    ///                        << ": ignoring misplaced special2 ranges line.");
    ///             return;
    ///         }
    ///         string s;
    ///         CTempString from;
    ///         char *p = &hints.prev_special2_acc[hints.prev_special2_acc.size()];
    ///         for (size_t i = 1;  i < tokens.size();  ++i) {
    ///             NStr::SplitInTwo(tokens[i], "-", tmp1, tmp2);
    ///             memcpy(p - tmp1.size(), tmp1.data(), tmp1.size());
    ///             if (tmp2.empty()) {
    ///                 from = hints.prev_special2_acc;
    ///             } else {
    ///                 s = hints.prev_special2_acc;
    ///                 from = s;
    ///                 memcpy(p - tmp2.size(), tmp2.data(), tmp2.size());
    ///             }
    ///             x_AddSpecial(*hints.special2_submap, hints,
    ///                          hints.prev_special_format, from,
    ///                          hints.prev_special2_acc, hints.special2_type,
    ///                          hints.special2_old_name.get(), hints.special2_name);
    ///         }
    /// ```
    fn add_special2_ranges(&mut self, ranges: &[&[u8]], hints: &mut Hints) {
        if hints.version < 2 {
            self.load_warnings
                .push("special2 continuation lines valid only in version 2+ guides".to_string());
        }
        let Some(fmt) = hints.special2_submap else {
            self.load_warnings
                .push("ignoring misplaced special2 ranges line.".to_string());
            return;
        };
        let value = hints.special2_type;
        let mut acc = std::mem::take(&mut hints.prev_special2_acc);
        for range in ranges {
            let (low, high) = match range.iter().position(|&c| c == b'-') {
                Some(dash) => (&range[..dash], &range[dash + 1..]),
                None => (&range[..], &b""[..]),
            };
            let end = acc.len();
            acc[end - low.len()..].copy_from_slice(low);
            if high.is_empty() {
                let from = acc.clone();
                self.add_special(hints, fmt, &from, &acc, value);
            } else {
                let from = acc.clone();
                acc[end - high.len()..].copy_from_slice(high);
                self.add_special(hints, fmt, &from, &acc, value);
            }
        }
        hints.prev_special2_acc = acc;
    }

    /// `SAccGuide::SAccGuide()` with the built-in rules, then the cleanup of empty small
    /// specials and `x_InitGeneral`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1443-1498
    /// ```c++
    ///         static const unsigned int kNumBuiltInRules
    ///             = sizeof(kBuiltInGuide) / sizeof(*kBuiltInGuide);
    ///         SHints hints;
    ///         for (unsigned int i = 0;  i < kNumBuiltInRules;  ++i) {
    ///             AddRule(kBuiltInGuide[i], hints);
    ///         }
    ///     }
    ///     for (auto &rit : rules) {
    ///         ERASE_ITERATE(TSmallSpecialMap, sit, rit.second.small_specials) {
    ///             if (sit->second.first.any()) {
    ///                 sit->second.first.optimize();
    ///             } else {
    ///                 rit.second.small_specials.erase(sit);
    ///             }
    ///         }
    ///     }
    ///     x_InitGeneral();
    /// }
    ///
    /// void SAccGuide::x_InitGeneral(void)
    /// {
    ///     if (general.empty()) {
    ///         // Populate with a hard-coded list by default; there are only
    ///         // a few tags to worry about, but listing them in accguide.txt
    ///         // right away would yield warnings from old Toolkit versions.
    ///         static const char* const kNucDBs[] = {
    ///             "SRA", "TI", "TR_ASSM_CH", "TRACE_ASSM", "TRACE_CHGR", NULL
    ///         };
    ///         for (const char* const* p = kNucDBs;  *p;  ++p) {
    ///             general[*p] = CSeq_id::eAcc_general_nuc;
    ///         }
    ///     }
    /// }
    /// ```
    /// The comments and blank lines of `kBuiltInGuide` are not in `rules.txt` (`AddRule`
    /// returns at once for them), and no data file is looked for (`g_FindDataFile`):
    /// makeblastdb 2.17.0 reads the built-in guide when no newer `accguide2.txt` is on its
    /// data path, and LOSAT always reads it.
    fn built_in() -> AccGuide {
        let mut guide = AccGuide {
            rules: BTreeMap::new(),
            general: BTreeMap::new(),
            load_warnings: Vec::new(),
        };
        let mut hints = Hints::default();
        let mut ranges = RangeLines::new(RANGES);
        for line in RULES.lines() {
            if line == ":" {
                while let Some(tokens) = ranges.next_line() {
                    let mut words: Vec<&[u8]> = vec![b":"];
                    words.extend(tokens.split(|&c| c == b' '));
                    guide.add_rule(&words, &mut hints);
                }
                continue;
            }
            let words: Vec<&[u8]> = line.split(' ').map(str::as_bytes).collect();
            guide.add_rule(&words, &mut hints);
        }
        for submap in guide.rules.values_mut() {
            submap.small_specials.retain(|_, options| {
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
        let value = if self.nibble % 2 == 0 {
            byte >> 4
        } else {
            byte & 15
        };
        self.nibble += 1;
        value
    }

    fn next_line(&mut self) -> Option<&[u8]> {
        self.line.clear();
        loop {
            match self.next_nibble() {
                digit @ 0..=9 => self.line.push(b'0' + digit),
                10 => self.line.push(b'-'),
                11 => self.line.push(b' '),
                12 => return Some(&self.line),
                _ => return None,
            }
        }
    }
}

/// The guide, read at its first use (`s_Guide`, a `CSafeStatic` made at the first call).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1631-1634
/// ```c++
/// static CRef<SAccGuide>* s_CreateGuide(void)
/// {
///     return new CRef<SAccGuide>(new SAccGuide);
/// }
///
/// static CSafeStatic<CRef<SAccGuide> > s_Guide(s_CreateGuide, NULL);
/// ```
pub fn guide() -> &'static AccGuide {
    static GUIDE: OnceLock<AccGuide> = OnceLock::new();
    GUIDE.get_or_init(AccGuide::built_in)
}

/// `CSeq_id::IdentifyAccession(acc, flags)`: the type of a bare accession. `warn` gets the
/// diagnostic of an accession with a protein flag and an unexpected division.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1634-1663
/// ```c++
/// CSeq_id::EAccessionInfo CSeq_id::IdentifyAccession(const CTempString& acc,
///                                                    TParseFlags flags)
/// {
///     SIZE_TYPE main_size = acc.find('.');
///     bool has_version = true;
///     if (main_size == NPOS) {
///         has_version = false;
///         main_size = acc.size();
///     } else if (main_size >= acc.size() - 1
///                ||  acc.find_first_not_of(kDigits, main_size + 1) != NPOS) {
///         return eAcc_unknown; // non-numeric "version"
///     }
///
///     static const SIZE_TYPE kMainAccBufSize = 32;
///     if (main_size <= kMainAccBufSize) {
///         const unsigned char* ucdata = (const unsigned char*)acc.data();
///         char main_acc_buf[kMainAccBufSize];
///         for (SIZE_TYPE i = 0;  i < main_size;  ++i) {
///             main_acc_buf[i] = toupper(ucdata[i]);
///         }
///         CTempString main_acc(main_acc_buf, main_size);
///         return x_IdentifyAccession(main_acc, flags, has_version);
///     } else {
///         // Unlikely to prove recognizable (far too long for any standard
///         // format as of January 2016), but try anyway.
///         string main_acc(acc, 0, main_size);
///         NStr::ToUpper(main_acc);
///         return x_IdentifyAccession(main_acc, flags, has_version);
///     }
/// }
/// ```
/// (`acc.size() - 1` of an empty accession wraps to `NPOS`, which no `find` result passes.)
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
    let main = acc[..main_size].to_ascii_uppercase();
    x_identify_accession(&main, flags, has_version, warn)
}

/// `CSeq_id::x_IdentifyAccession(main_acc, flags, has_version)`.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1665-1769
/// ```c++
/// CSeq_id::EAccessionInfo
/// CSeq_id::x_IdentifyAccession(const CTempString& main_acc, TParseFlags flags,
///                              bool has_version)
/// {
///     SIZE_TYPE digit_pos = main_acc.find_first_of(kDigits),
///         main_size = main_acc.size();
///     char flag_char = '\0';
///     if (digit_pos == NPOS) {
///         return eAcc_unknown;
///     } else {
///         SIZE_TYPE non_dig_pos = main_acc.find_first_not_of(kDigits, digit_pos);
///         const unsigned char* ucdata = (const unsigned char*)main_acc.data();
///         if (non_dig_pos != NPOS  &&  (flags & fParse_RawText) != 0) {
///             if ( !has_version  &&  digit_pos == 0  &&  main_size >= 4
///                 &&  non_dig_pos < 5  &&  isalnum(ucdata[1])
///                 &&  isalnum(ucdata[2])  &&  isalnum(ucdata[3])) {
///                 // Possible PDB (always unversioned); examine further
///                 // to avoid false positives.
///                 if (main_size > 4  &&  main_size <= 17
///                     &&  strchr("|-_", main_acc[4])
///                     &&  (main_size <= 6  ||  isalnum(ucdata[5]))) {
///                     // Conventionally delimited
///                     return eAcc_pdb;
///                 } else switch (main_size) {
///                 /* (cases 7, 6 and 5 commented out) */
///                 case 4:
///                     return eAcc_pdb;
///                 }
///             }
///             if (digit_pos == 1  &&  main_size == 6
///                 &&  (main_acc[0] == 'O'  ||  main_acc[0] == 'P'
///                      ||  main_acc[0] == 'Q' ||  isalpha(ucdata[2]))
///                 &&  isdigit(ucdata[1])  &&  isalnum(ucdata[2])
///                 &&  isalnum(ucdata[3])  &&  isalnum(ucdata[4])
///                 &&  isdigit(ucdata[5])) {
///                 return eAcc_swissprot;
///             } else if (digit_pos == 1  &&  main_size == 10
///                        &&  main_acc[0] != 'O'  &&  main_acc[0] != 'P'
///                        &&  main_acc[0] != 'Q'
///                        &&  isalpha(ucdata[2])  &&  isalnum(ucdata[3])
///                        &&  isalnum(ucdata[4])  &&  isdigit(ucdata[5])
///                        &&  isalpha(ucdata[6])  &&  isalnum(ucdata[7])
///                        &&  isalnum(ucdata[8])  &&  isdigit(ucdata[9])) {
///                 return eAcc_swissprot;
///             } else if ( !has_version  &&  digit_pos == 0
///                        &&  (non_dig_pos == 6  ||  non_dig_pos == 7)
///                        &&  (main_size == non_dig_pos + 1
///                             ||  main_acc[non_dig_pos + 1] == ':'
///                             ||  (isalpha(ucdata[non_dig_pos + 1])
///                                  &&  (main_size == non_dig_pos + 2
///                                       ||  main_acc[non_dig_pos + 2] == ':')))) {
///                 // A formal spec appears to be elusive, but all examples in ID
///                 // contain six or seven digits followed by one or two letters,
///                 // followed in some rare cases by a tag such as :PDB=...
///                 return eAcc_prf;
///             } else if (digit_pos >= 4  &&  non_dig_pos == digit_pos + 2
///                        &&  main_size - non_dig_pos >= 6  &&  main_acc[3] != '_'
///                        &&  (main_acc[non_dig_pos] == 'S'
///                             ||  main_acc[non_dig_pos] == 'P')
///                        &&  (main_acc.find_first_not_of
///                             (kDigits, non_dig_pos + 1) == NPOS)) {
///                 flag_char = main_acc[non_dig_pos];
///             } else {
///                 return eAcc_unknown;
///             }
///         }
///     }
///
///     if (digit_pos == 0) {
///         if ((flags & fParse_RawGI) != 0  &&  !has_version
///             &&  main_acc[0] != '0'
///             &&  main_acc.find_first_not_of(kDigits) == NPOS) {
///             return eAcc_gi; // just digits
///         } else {
///             return eAcc_unknown; // PDB already handled
///         }
///     } else if ((flags & fParse_RawText) == 0) {
///         return eAcc_unknown;
///     }
///
///     SIZE_TYPE flag_len = (flag_char == '\0') ? 0 : 1;
///     SIZE_TYPE digit_count = main_size - digit_pos - flag_len;
///     auto& guide = *s_Guide;
///     const EAccessionInfo& found_ai
///         = guide->Find(SAccGuide::s_Key(digit_pos, digit_count), main_acc);
///     EAccessionInfo ai = found_ai;
///     if ((ai & fAcc_specials) != 0) {
///         ai = EAccessionInfo(ai & ~fAcc_specials);
///     }
/// ```
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
    // ```c++
    //     if (flag_char == 'P') {
    //         switch (ai & eAcc_division_mask) {
    //         case eAcc_targeted:
    //         case eAcc_wgs:
    //         // case eAcc_wgs_intermed:
    //             ai = EAccessionInfo((ai & eAcc_type_mask) | eAcc_wgs | fAcc_prot);
    //             break;
    //         case eAcc_tsa:
    //             ai = EAccessionInfo((ai & eAcc_type_mask) | eAcc_tsa | fAcc_prot);
    //             break;
    //         default:
    //             ERR_POST_X(11,
    //                        Warning << main_acc
    //                        << ": Protein flag found with unexpected division "
    //                        << ((ai & eAcc_division_mask) >> 8));
    //             ai = EAccessionInfo((ai & (eAcc_type_mask | eAcc_division_mask))
    //                                 | fAcc_prot);
    //             break;
    //         }
    //     }
    //     switch (ai & eAcc_division_mask) {
    //     case eAcc_targeted:
    //     case eAcc_tsa:
    //     case eAcc_wgs:
    //     case eAcc_wgs_intermed:
    //         if (digit_pos >= 4
    //             &&  (main_acc.find_first_not_of
    //                  ("0", digit_pos /* + flag_len */ + 2) == NPOS)) {
    //             return EAccessionInfo(ai | fAcc_master);
    //         }
    //     default:
    //         break;
    //     }
    //     return ai;
    // }
    // ```
    if flag_char == b'P' {
        match ai & DIVISION_MASK {
            TARGETED | WGS => ai = (ai & TYPE_MASK) | WGS | F_PROT,
            TSA => ai = (ai & TYPE_MASK) | TSA | F_PROT,
            division => {
                warn(SeqIdWarning {
                    subcode: 11,
                    message: format!(
                        "{}: Protein flag found with unexpected division {}",
                        String::from_utf8_lossy(main),
                        division >> 8
                    ),
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

/// `CSeq_id::IdentifyAccession(flags)` of a Seq-id (the member function; its default flags
/// are `fParse_AnyRaw | fParse_AnyLocal`).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1832-1876
/// ```c++
/// CSeq_id::EAccessionInfo CSeq_id::IdentifyAccession(TParseFlags flags) const
/// {
///     E_Choice type = Which();
///     switch (type) {
///     case e_Pir: case e_Swissprot: case e_Prf: // but *NOT* e_Pdb
///         // always just protein
///         return (EAccessionInfo)((EAccessionInfo)type | fAcc_prot);
///
///     case e_Genbank: case e_Embl:  case e_Ddbj:
///     case e_Tpg:     case e_Tpe:   case e_Tpd:
///     case e_Other:   case e_Gpipe: case e_Named_annot_track:
///     {
///         const CTextseq_id* tsid = GetTextseq_Id();
///         if (tsid->IsSetAccession()) {
///             // Can't necessarily go straight to x_IdentifyAccession, as
///             // the accession may contain lowercase letters.
///             EAccessionInfo ai = IdentifyAccession(tsid->GetAccession(), flags);
///             E_Choice type2 = GetAccType(ai);
///             auto div2 = ai & eAcc_division_mask;
///             if (type2 == e_not_set) {
///                 // We *know* what the type should be....
///                 return (EAccessionInfo)((ai & eAcc_flag_mask) | type);
///             } else if (type2 == type) {
///                 return ai;
///             } else if (type == e_Tpe  &&  type2 == e_Embl
///                        &&  (div2 == eAcc_other  ||  div2 == eAcc_wgs)) {
///                 return (EAccessionInfo)((ai & ~eAcc_type_mask) | type);
///             } else { // misidentified or mislabeled; assume the former
///                 return static_cast<EAccessionInfo>(type);
///             }
///         } else {
///             return static_cast<EAccessionInfo>(type);
///         }
///     }
///
///     case e_General:
///     {
///         string db = GetGeneral().GetDb();
///         NStr::ToUpper(db);
///         SAccGuide::TPrefixes::const_iterator it = (*s_Guide)->general.find(db);
///         return it == (*s_Guide)->general.end() ? eAcc_general : it->second;
///     }
///
///     default:
///         return static_cast<EAccessionInfo>(type);
///     }
/// }
/// ```
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
