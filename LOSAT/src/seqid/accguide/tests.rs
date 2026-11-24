use super::*;

/// FNV-1a 64 of the compacted guide lines, as `make_accguide.py` printed it for
/// `accguide2.inc` of 2.17.0 (598d8ae6): 25,549 lines, 1,737,234 bytes.
const GUIDE_FNV1A64: u64 = 0x0f89_a3ba_8888_0aef;

fn fnv1a64(data: &[u8]) -> u64 {
    data.iter().fold(0xcbf2_9ce4_8422_2325, |hash, &byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3)
    })
}

#[test]
fn the_two_files_hold_the_lines_of_kbuiltinguide() {
    let mut text = Vec::new();
    let mut lines = 0usize;
    let mut tokens = 0usize;
    let mut ranges = RangeLines::new(RANGES);
    for line in RULES.lines() {
        if line == ":" {
            while let Some(range_line) = ranges.next_line() {
                text.extend_from_slice(b": ");
                text.extend_from_slice(range_line);
                text.push(b'\n');
                lines += 1;
                tokens += range_line.split(|&c| c == b' ').count();
            }
        } else {
            text.extend_from_slice(line.as_bytes());
            text.push(b'\n');
            lines += 1;
        }
    }
    assert_eq!((lines, tokens, text.len()), (25_549, 456_251, 1_737_234));
    assert_eq!(fnv1a64(&text), GUIDE_FNV1A64);
}

#[test]
fn the_built_in_guide_loads_without_diagnostics() {
    assert_eq!(guide().load_warnings, Vec::<String>::new());
    let all = AccGuide::built_in(Part::All);
    assert_eq!(all.load_warnings, Vec::<String>::new());
    // No value is a fallback (every type name is known).
    for submap in all.rules.values() {
        assert!(submap.prefixes.values().all(|v| v & F_FALLBACK == 0));
        assert!(submap.wildcards.iter().all(|(mask, v)| v & F_FALLBACK == 0
            && mask
                .iter()
                .all(|&c| c.is_ascii_uppercase() || c == b'_' || c == b'?')));
        let specials = submap
            .specials
            .get()
            .expect("a build of every part makes the specials");
        assert!(specials
            .big_specials
            .values()
            .all(|(_, v)| v & F_FALLBACK == 0));
        assert!(specials
            .small_specials
            .values()
            .flat_map(|&slot| &specials.small_options[slot].1)
            .all(|(bits, v)| v & F_FALLBACK == 0 && bits.any()));
    }
    assert_eq!(guide().general.get(&b"SRA"[..]), Some(&GENERAL_NUC));
}

fn ident(acc: &str) -> AccInfo {
    identify_accession(
        acc.as_bytes(),
        parse_flags::RAW_TEXT | parse_flags::FALLBACK_OK,
        &mut |w| panic!("unexpected warning {w:?}"),
    )
}

fn name_of(info: AccInfo) -> &'static str {
    ACC_INFO_NAMES
        .iter()
        .find(|(_, value)| *value == info)
        .map(|(name, _)| *name)
        .unwrap_or("?")
}

/// Accessions and the types the guide's lines give them (the rule or special quoted).
#[test]
fn accessions_get_the_types_of_the_guide() {
    for (acc, expected) in [
        ("XP_011520000.1", "refseq_prot_predicted"), // 3+9 XP_ refseq_prot_predicted
        ("NP_000549", "refseq_prot"),                // 3+6 NP_ refseq_prot
        ("NC_000913.3", "refseq_chromosome"),        // 3+6 NC_ refseq_chromosome
        ("WP_000000001.1", "refseq_unique_prot"),    // 3+9 WP_ refseq_unique_prot
        ("P69906", "swissprot"),                     // x_IdentifyAccession
        ("A0A024R161", "swissprot"),
        ("1ABC", "pdb"),
        ("1ABC_A", "pdb"),
        ("1234567A", "prf"),
        ("N00001", "gb_embl"),     // special N00001 gb_embl
        ("N00006", "gb_ddbj"),     // special N00006 gb_ddbj
        ("F38006", "embl_dirsub"), // special F38006 embl_dirsub
        ("L48736", "gb_gss"),      // special L48736-L48796 gb_gss
        ("L48796", "gb_gss"),
        ("AAA12345", "gb_prot"),   // 3+5 A?? gb_prot
        ("CAA12345", "embl_prot"), // 3+5 C?? embl_prot * (no special)
        ("BAA12345", "ddbj_prot"), // 3+5 B?? ddbj_prot
        ("UNIREF50_A0A5A9P0L4", "unknown"),
        ("CHR1", "unknown"),
        ("123456", "unknown"), // no fParse_RawGI
        ("NC_000914.1:1-50", "unknown"),
    ] {
        assert_eq!(name_of(ident(acc)), expected, "{acc}");
    }
    assert_eq!(
        identify_accession(b"123456", parse_flags::ANY_RAW, &mut |_| {}),
        GI
    );
}

#[test]
fn bit_ranges_set_and_clear_as_a_bit_vector() {
    let mut bits = BitRanges::default();
    bits.set_range(10, 20);
    bits.set_range(30, 40);
    bits.set_range(21, 29);
    assert_eq!(bits.0, [(10, 40)]);
    bits.clear_range(15, 35);
    assert_eq!(bits.0, [(10, 14), (36, 40)]);
    bits.set_range(0, 0);
    bits.set_range(5, 50);
    assert_eq!(bits.0, [(0, 0), (5, 50)]);
    bits.clear_range(0, 100);
    assert!(!bits.any());
    bits.set_range(7, 7);
    assert!(bits.test(7) && !bits.test(6) && !bits.test(8));
}

#[test]
fn masks_match_as_ncbi_s() {
    assert!(matches_mask(b"ABC", b"A??"));
    assert!(!matches_mask(b"ABC", b"B??"));
    assert!(!matches_mask(b"AB\0", b"A??"));
    assert!(!matches_mask(b"ABCD", b"A??"));
}

/// `set_range` and `clear_range` in any order of ranges give the bits of a plain bit vector, in
/// sorted, disjoint, non-adjacent ranges.
#[test]
fn bit_ranges_match_a_bit_vector_in_any_order() {
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next = |bound: u32| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % u64::from(bound)) as u32
    };
    for round in 0..300 {
        let mut bits = BitRanges::default();
        let mut plain = [false; 160];
        for _ in 0..(round % 40 + 1) {
            let left = next(150);
            let right = left + next(10);
            let set = next(3) != 0;
            // BitMagic takes the two ends in either order.
            let (first, second) = if next(4) == 0 {
                (right, left)
            } else {
                (left, right)
            };
            if set {
                bits.set_range(first, second);
            } else {
                bits.clear_range(first, second);
            }
            plain[left as usize..=right as usize].fill(set);
            let mut expected = Vec::new();
            let mut bit = 0;
            while bit < plain.len() {
                if plain[bit] {
                    let start = bit;
                    while bit + 1 < plain.len() && plain[bit + 1] {
                        bit += 1;
                    }
                    expected.push((start as u32, bit as u32));
                }
                bit += 1;
            }
            assert_eq!(bits.0, expected, "round {round}");
        }
    }
}

/// `skip_block` stops where `next_line` returns `None`, from any line of any block.
#[test]
fn skip_block_ends_where_the_lines_end() {
    let blocks = RULES.lines().filter(|&line| line == ":").count();
    assert_eq!(blocks, 448);
    let mut lines = RangeLines::new(RANGES);
    for block in 0..blocks {
        let mut skipped = RangeLines::new(RANGES);
        skipped.nibble = lines.nibble;
        for _ in 0..block % 3 {
            if skipped.next_line().is_none() {
                skipped.nibble -= 1;
                break;
            }
        }
        skipped.skip_block();
        while lines.next_line().is_some() {}
        assert_eq!(skipped.nibble, lines.nibble, "block {block}");
    }
    assert!(lines.nibble + 1 >= 2 * RANGES.len());
}

/// The guide in a form that both builders give: for each format, its prefixes, wildcards, big
/// specials, and the options of each small special with their ranges.
type Canonical = Vec<(
    u32,
    Vec<(Vec<u8>, AccInfo)>,
    Vec<(Vec<u8>, AccInfo)>,
    Vec<(Vec<u8>, Vec<u8>, AccInfo)>,
    Vec<(Vec<u8>, Vec<(Vec<(u32, u32)>, AccInfo)>)>,
)>;

/// The canonical form of `guide`, with the specials of each format that `specials` gives.
fn canonical_new<'a>(
    guide: &'a AccGuide,
    specials: impl Fn(u32, &'a SubMap) -> &'a Specials,
) -> Canonical {
    guide
        .rules
        .iter()
        .map(|(&fmt, submap)| {
            let specials = specials(fmt, submap);
            (
                fmt,
                submap.prefixes.clone().into_iter().collect(),
                submap.wildcards.clone(),
                specials
                    .big_specials
                    .iter()
                    .map(|(to, (from, value))| (to.clone(), from.clone(), *value))
                    .collect(),
                specials
                    .small_specials
                    .iter()
                    .map(|(key, &slot)| {
                        assert_eq!(&specials.small_options[slot].0, key);
                        let options = specials.small_options[slot]
                            .1
                            .iter()
                            .map(|(bits, value)| (bits.0.clone(), *value))
                            .collect();
                        (key.clone(), options)
                    })
                    .collect(),
            )
        })
        .collect()
}

/// FNV-1a 64 of the canonical form: every format, key, value and range of the guide, each
/// list after its length, each item after its length in bytes (`make_accguide.py` writes the
/// same stream).
fn fingerprint(canonical: &Canonical) -> u64 {
    let mut bytes = Vec::new();
    let mut put = |data: &[u8]| {
        bytes.extend_from_slice(&(data.len() as u32).to_le_bytes());
        bytes.extend_from_slice(data);
    };
    let count = |n: usize| (n as u32).to_le_bytes();
    put(&count(canonical.len()));
    for (fmt, prefixes, wildcards, big, small) in canonical {
        put(&fmt.to_le_bytes());
        for list in [prefixes, wildcards] {
            put(&count(list.len()));
            for (key, value) in list {
                put(key);
                put(&value.to_le_bytes());
            }
        }
        put(&count(big.len()));
        for (to, from, value) in big {
            put(to);
            put(from);
            put(&value.to_le_bytes());
        }
        put(&count(small.len()));
        for (key, options) in small {
            put(key);
            put(&count(options.len()));
            for (ranges, value) in options {
                put(&value.to_le_bytes());
                put(&count(ranges.len()));
                for (left, right) in ranges {
                    put(&left.to_le_bytes());
                    put(&right.to_le_bytes());
                }
            }
        }
    }
    fnv1a64(&bytes)
}

/// The fingerprint of NCBI's built-in guide, as `make_accguide.py` prints it (its own model of
/// `SAccGuide`'s loading, from `accguide2.inc` and the values of `acc_info.rs`): run it again
/// when the guide or a value changes. The guide of dev 2cd5a1906's builder had it too (the
/// tests of b259ca129 to af9d5619e compared that builder's guide with this one, whole and by
/// part, and every lookup at the edges of the ranges, 2026-11-24).
const GUIDE_FINGERPRINT: u64 = 0xe9e2_88a9_4ded_e627;

/// The guide has the fingerprint of NCBI's, whole, and part by part: the rules without the
/// specials, and the specials of each format built alone.
#[test]
fn the_guide_has_ncbi_s_fingerprint_whole_and_by_part() {
    let all = AccGuide::built_in(Part::All);
    let whole = canonical_new(&all, |_, submap| {
        submap
            .specials
            .get()
            .expect("made by a build of every part")
    });
    assert_eq!(fingerprint(&whole), GUIDE_FINGERPRINT);
    assert_eq!(all.general.len(), 5);
    let rules = AccGuide::built_in(Part::Rules);
    assert!(rules.load_warnings.is_empty());
    assert_eq!(rules.general, all.general);
    let specials: BTreeMap<u32, Specials> = rules
        .rules
        .keys()
        .map(|&fmt| (fmt, Specials::built_in(fmt)))
        .collect();
    let parts = canonical_new(&rules, |fmt, _| &specials[&fmt]);
    let mut parts = parts.iter().peekable();
    for entry in &whole {
        if parts.peek().is_some_and(|part| part.0 == entry.0) {
            assert!(parts.next() == Some(entry), "format {:#x} differs", entry.0);
        } else {
            // A format of `special` lines alone, which no lookup reaches (every format of the
            // built-in guide has a prefix rule or a wildcard: this serves another guide).
            assert!(
                entry.1.is_empty() && entry.2.is_empty(),
                "format {:#x}",
                entry.0
            );
        }
    }
    assert!(parts.next().is_none());
}

/// The lookups of the guide (`Find`, with the key it reports), which reads the specials of a
/// format at need, answer as the guide built whole does: at both ends of every range of every
/// small special and next to them, at both ends of every big special and next to them, and for
/// random accessions of every prefix and wildcard.
#[test]
fn lookups_answer_as_with_the_whole_guide() {
    let whole = AccGuide::built_in(Part::All);
    let new = guide();
    let mut lookups = 0usize;
    let mut check = |fmt: u32, acc: &[u8]| {
        let (mut whole_key, mut new_key) = (Vec::new(), Vec::new());
        let expected = whole.find(fmt, acc, Some(&mut whole_key));
        let got = new.find(fmt, acc, Some(&mut new_key));
        assert_eq!(
            (got, &new_key),
            (expected, &whole_key),
            "{fmt:#x} {}",
            String::from_utf8_lossy(acc)
        );
        assert_eq!(new.find(fmt, acc, None), expected);
        lookups += 1;
    };
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut next = |bound: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % bound
    };
    for (&fmt, submap) in &whole.rules {
        let letters = (fmt >> 16) as usize;
        let digits = (fmt & 0xffff) as usize;
        let size = letters + digits;
        let specials = submap
            .specials
            .get()
            .expect("made by a build of every part");
        for (key, &slot) in &specials.small_specials {
            // The digits that the bit is cut from (`x_SplitSpecial` of an accession of the
            // format's size).
            let Some(width) = size
                .checked_sub(key.len())
                .filter(|&w| (1..=9).contains(&w))
            else {
                continue;
            };
            let limit = 10u64.pow(width as u32);
            for (bits, _) in &specials.small_options[slot].1 {
                for &(left, right) in &bits.0 {
                    let (left, right) = (u64::from(left), u64::from(right));
                    for bit in [left.wrapping_sub(1), left, right, right + 1] {
                        if bit < limit {
                            let mut acc = key.clone();
                            acc.extend_from_slice(format!("{bit:0width$}").as_bytes());
                            check(fmt, &acc);
                        }
                    }
                }
            }
        }
        for (to, (from, _)) in &specials.big_specials {
            for acc in [from, to] {
                check(fmt, acc);
                for step in [b'0', b'9'] {
                    let mut near = acc.clone();
                    if let Some(last) = near.last_mut() {
                        if last.is_ascii_digit() && *last != step {
                            *last = if step == b'0' { *last - 1 } else { *last + 1 };
                            check(fmt, &near);
                        }
                    }
                }
            }
        }
        let keys = submap
            .prefixes
            .iter()
            .map(|(key, value)| (key, *value))
            .chain(submap.wildcards.iter().map(|(mask, value)| (mask, *value)));
        for (key, value) in keys {
            let tries = if value & F_SPECIALS != 0 { 3000 } else { 20 };
            for _ in 0..tries {
                let mut acc: Vec<u8> = key
                    .iter()
                    .map(|&c| if c == b'?' { b'A' + next(26) as u8 } else { c })
                    .collect();
                acc.truncate(letters);
                while acc.len() < size {
                    acc.push(b'0' + next(10) as u8);
                }
                check(fmt, &acc);
                check(fmt, &acc[..letters.min(acc.len())]);
            }
        }
    }
    assert!(lookups > 1_500_000, "{lookups}");
}

/// `rule_lines` splits as `str::lines`.
#[test]
fn rule_lines_split_as_str_lines() {
    for text in [
        RULES, "", "a", "a\n", "a\r\n", "a\n\nb", "a\r", "\r\n", "x\ry\n", "\n\n", "a\r\r\nb",
    ] {
        let expected: Vec<&[u8]> = text.lines().map(str::as_bytes).collect();
        let got: Vec<&[u8]> = rule_lines(text.as_bytes()).collect();
        assert!(got == expected, "{text:?}");
    }
}

/// Every range of every `special2` block of the built-in guide is `<digits>[-<digits>]` with
/// at most `bit_digits` digits a side, the one form that `add_special2_ranges` reads, and each
/// block follows its `special2` line.
#[test]
fn every_special2_range_is_of_the_bit_s_digits() {
    let mut ranges = RangeLines::new(RANGES);
    let mut previous: Option<&str> = None;
    let mut tokens = 0usize;
    for line in RULES.lines() {
        if line == ":" {
            let spec = previous
                .and_then(|rule| rule.strip_prefix("special2 "))
                .and_then(|rest| rest.split(' ').next())
                .expect("a block follows its special2 line");
            let (letters, digits) = spec.split_once('+').expect("<letters>+<digits>");
            let fmt = format_key(letters.len(), digits.parse().expect("digits"));
            let bit = bit_digits(fmt) as usize;
            while let Some(range_line) = ranges.next_line() {
                for range in range_line.split(|&c| c == b' ') {
                    assert!(
                        range_numbers(range, bit).is_some(),
                        "{spec} {}",
                        String::from_utf8_lossy(range)
                    );
                    tokens += 1;
                }
            }
        }
        previous = Some(line);
    }
    assert_eq!(tokens, 456_251);
}

/// The ranges of a `special2` block added from the text of the accession `acc`, range by range
/// (NCBI's loop).
fn add_ranges_from_text(guide: &mut AccGuide, hints: &mut Hints, acc: &mut [u8], ranges: &[&[u8]]) {
    let fmt = hints.special2_submap.expect("a special2 block");
    let value = hints.special2_type;
    for range in ranges {
        let (low, high) = match range.iter().position(|&c| c == b'-') {
            Some(dash) => (&range[..dash], &range[dash + 1..]),
            None => (&range[..], &b""[..]),
        };
        let end = acc.len();
        acc[end - low.len()..].copy_from_slice(low);
        let from = acc.to_vec();
        if !high.is_empty() {
            acc[end - high.len()..].copy_from_slice(high);
        }
        AccGuide::add_special(guide.specials_mut(fmt), hints, fmt, &from, acc, value);
    }
}

/// Random digits for one side of a range, at most `bit` of them.
fn random_range_part(next: &mut impl FnMut(u64) -> u64, bit: usize) -> Vec<u8> {
    let len = next(1 + bit as u64) as usize;
    (0..len).map(|_| b'0' + next(10) as u8).collect()
}

/// `add_special2_ranges` adds random ranges as NCBI's loop adds them from the text, in formats
/// of every split of `x_SplitSpecial` (no digits, fewer than 10, 10, more than 10), with blocks
/// of the same key and of others one after the other.
#[test]
fn special2_ranges_add_as_from_the_text() {
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next = |bound: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % bound
    };
    let empty = || AccGuide {
        rules: BTreeMap::new(),
        general: BTreeMap::new(),
        load_warnings: Vec::new(),
        part: Part::All,
    };
    for (letters, digits) in [(2usize, 0usize), (2, 6), (3, 7), (1, 9), (2, 10), (4, 12)] {
        let (mut numeric, mut text) = (empty(), empty());
        let (mut numeric_hints, mut text_hints) = (Hints::default(), Hints::default());
        let size = letters + digits;
        let fmt = format_key(letters, digits);
        let bit = bit_digits(fmt) as usize;
        let spec = format!("{}+{digits}", &"ABCD"[..letters]);
        let types = ["embl_prot", "ddbj_prot", "gb_prot"];
        let mut acc = Vec::new();
        for block in 0..60 {
            let mut letters_now = spec.clone().into_bytes();
            if block % 4 == 3 {
                letters_now[0] = b'Q';
            }
            acc.clear();
            acc.extend_from_slice(&letters_now[..letters]);
            acc.resize(size, b'0');
            let type_name = types[next(3) as usize].as_bytes();
            let header: [&[u8]; 3] = [b"special2", &letters_now, type_name];
            for (guide, hints) in [
                (&mut numeric, &mut numeric_hints),
                (&mut text, &mut text_hints),
            ] {
                guide.add_rule(&[&b"version"[..], &b"2"[..]], hints);
                guide.add_rule(&header, hints);
            }
            for _ in 0..1 + next(4) {
                let mut tokens: Vec<Vec<u8>> = Vec::new();
                for _ in 0..1 + next(30) {
                    let mut token = random_range_part(&mut next, bit);
                    if next(4) != 0 {
                        token.push(b'-');
                        token.extend(random_range_part(&mut next, bit));
                    }
                    tokens.push(token);
                }
                let ranges: Vec<&[u8]> = tokens.iter().map(Vec::as_slice).collect();
                numeric.add_special2_ranges(ranges.iter().copied(), &mut numeric_hints);
                add_ranges_from_text(&mut text, &mut text_hints, &mut acc, &ranges);
                let mut running = numeric_hints.special2_key.clone();
                if bit > 0 {
                    running.extend(format!("{:0bit$}", numeric_hints.special2_bit).bytes());
                }
                assert_eq!(
                    String::from_utf8_lossy(&running),
                    String::from_utf8_lossy(&acc)
                );
                assert_eq!(
                    numeric_hints.prev_small_special,
                    text_hints.prev_small_special
                );
            }
        }
        let specials = |guide: &AccGuide| {
            canonical_new(guide, |_, submap| {
                submap.specials.get().expect("a build of every part")
            })
        };
        assert!(specials(&numeric) == specials(&text), "{letters}+{digits}");
        assert!(!specials(&numeric)[0].4.is_empty());
    }
}

/// The emplace and set of `x_AddSpecial` for a small special as NCBI does it: the bits are
/// cleared in the hinted option before the option of the value is searched.
fn add_to_options_as_ncbi(
    options: &mut Vec<(BitRanges, AccInfo)>,
    prev: Option<usize>,
    left: u32,
    right: u32,
    value: AccInfo,
) -> usize {
    let same = |option: &(BitRanges, AccInfo)| option.1 & !F_FALLBACK == value & !F_FALLBACK;
    let mut index = match prev {
        Some(prev) => {
            options[prev].0.clear_range(left, right);
            (0..=prev)
                .rev()
                .find(|&it| same(&options[it]))
                .unwrap_or(prev + 1)
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

/// `add_to_options` (which skips the clearing when the hinted option has the value) leaves
/// the options and the hint as NCBI's order does, for random ranges, values and hints.
#[test]
fn small_special_options_add_as_ncbi_s() {
    let mut state = 0x1234_5678_9abc_def1_u64;
    let mut next = |bound: u32| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % u64::from(bound)) as u32
    };
    let values = [
        acc(SEQID_EMBL, OTHER, PROT),
        acc(SEQID_DDBJ, OTHER, PROT),
        acc(SEQID_EMBL, OTHER, PROT) | F_FALLBACK,
        acc(SEQID_GENBANK, OTHER, PROT),
    ];
    for round in 0..400 {
        let (mut ours, mut ncbi) = (Vec::new(), Vec::new());
        let mut hint = None;
        for _ in 0..(round % 60 + 1) {
            let left = next(200);
            let right = left + next(12);
            let value = values[next(4) as usize];
            let prev = if next(5) == 0 { None } else { hint };
            let a = AccGuide::add_to_options(&mut ours, prev, left, right, value);
            let b = add_to_options_as_ncbi(&mut ncbi, prev, left, right, value);
            assert_eq!(a, b, "round {round}");
            let plain = |options: &Vec<(BitRanges, AccInfo)>| -> Vec<(Vec<(u32, u32)>, AccInfo)> {
                options
                    .iter()
                    .map(|(bits, v)| (bits.0.clone(), *v))
                    .collect()
            };
            assert_eq!(plain(&ours), plain(&ncbi), "round {round}");
            hint = Some(a);
        }
    }
}
