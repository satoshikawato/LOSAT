//! dev's builder of the accession guide before line NP M (f) (2cd5a1906, comments removed): the
//! oracle of `tests.rs`, which compares the guide it builds, and its lookups, with the new
//! builder's.

use std::collections::BTreeMap;
use std::ops::Bound;

use super::super::{
    find_acc_info, AccInfo, F_FALLBACK, F_SPECIALS, GENERAL_NUC, RANGES, RULES, UNKNOWN,
};

#[derive(Default, Debug)]
pub(super) struct BitRanges(pub(super) Vec<(u32, u32)>);

impl BitRanges {
    fn set_range(&mut self, left: u32, right: u32) {
        let ranges = &mut self.0;
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

    fn test(&self, bit: u32) -> bool {
        let index = self.0.partition_point(|&(_, end)| end < bit);
        index < self.0.len() && self.0[index].0 <= bit
    }

    fn any(&self) -> bool {
        !self.0.is_empty()
    }
}

#[derive(Default, Debug)]
pub(super) struct SubMap {
    pub(super) prefixes: BTreeMap<Vec<u8>, AccInfo>,
    pub(super) wildcards: Vec<(Vec<u8>, AccInfo)>,
    pub(super) big_specials: BTreeMap<Vec<u8>, (Vec<u8>, AccInfo)>,
    pub(super) small_specials: BTreeMap<Vec<u8>, Vec<(BitRanges, AccInfo)>>,
}

#[derive(Debug)]
pub(super) struct AccGuide {
    pub(super) rules: BTreeMap<u32, SubMap>,
    pub(super) general: BTreeMap<Vec<u8>, AccInfo>,
    pub(super) load_warnings: Vec<String>,
}

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

fn format_key(letters: usize, digits: usize) -> u32 {
    ((letters as u32 & 0xffff) << 16) | (digits as u32 & 0xffff)
}

fn bit_number(digits: &[u8]) -> u32 {
    crate::ncbi_str::string_to_uint(digits).unwrap_or(0)
}

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
    pub(super) fn find(
        &self,
        fmt: u32,
        acc_or_pfx: &[u8],
        mut key_used: Option<&mut Vec<u8>>,
    ) -> AccInfo {
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
            if let Some((_, (first, value))) = submap
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

    fn find_sub_map(&mut self, hints: &mut Hints, fmt: u32) -> u32 {
        if hints.prev_submap != Some(fmt) {
            self.rules.entry(fmt).or_default();
            hints.prev_submap = Some(fmt);
            hints.prev_small_special = None;
        }
        fmt
    }

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
            hints.special2_submap = None;
            match Self::find_acc_info_hinted(hints, tokens[2]) {
                Some(value) => {
                    self.general.insert(tokens[1].to_ascii_uppercase(), value);
                }
                None => self.unknown_name(tokens[2]),
            }
        } else if tokens.len() == 3 && first.eq_ignore_ascii_case(b"fallback") {
            hints.special2_submap = None;
        } else {
            self.load_warnings.push(format!(
                "ignoring invalid line: {}",
                String::from_utf8_lossy(&tokens.join(&b' '))
            ));
        }
    }

    fn unknown_name(&mut self, name: &[u8]) {
        self.load_warnings.push(format!(
            "unrecognized accession type {}",
            String::from_utf8_lossy(name)
        ));
    }

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

    pub(super) fn built_in() -> AccGuide {
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
