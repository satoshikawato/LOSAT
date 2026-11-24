use super::*;

#[path = "tests_old.rs"]
mod old;

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
            if set {
                bits.set_range(left, right);
            } else {
                bits.clear_range(left, right);
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

fn canonical_old(guide: &old::AccGuide) -> Canonical {
    guide
        .rules
        .iter()
        .map(|(&fmt, submap)| {
            (
                fmt,
                submap.prefixes.clone().into_iter().collect(),
                submap.wildcards.clone(),
                submap
                    .big_specials
                    .iter()
                    .map(|(to, (from, value))| (to.clone(), from.clone(), *value))
                    .collect(),
                submap
                    .small_specials
                    .iter()
                    .map(|(key, options)| {
                        let options = options
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

/// The canonical form of the new guide `guide`, with the specials of each format that
/// `specials` gives.
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

/// FNV-1a 64 of the canonical form: every format, key, value and range of the guide.
fn fingerprint(canonical: &Canonical) -> u64 {
    let mut bytes = Vec::new();
    let mut put = |data: &[u8]| {
        bytes.extend_from_slice(&(data.len() as u32).to_le_bytes());
        bytes.extend_from_slice(data);
    };
    for (fmt, prefixes, wildcards, big, small) in canonical {
        put(&fmt.to_le_bytes());
        for (key, value) in prefixes.iter().chain(wildcards) {
            put(key);
            put(&value.to_le_bytes());
        }
        for (to, from, value) in big {
            put(to);
            put(from);
            put(&value.to_le_bytes());
        }
        for (key, options) in small {
            put(key);
            for (ranges, value) in options {
                put(&value.to_le_bytes());
                for (left, right) in ranges {
                    put(&left.to_le_bytes());
                    put(&right.to_le_bytes());
                }
            }
        }
    }
    fnv1a64(&bytes)
}

/// The new builder makes the guide of dev's builder (`tests_old.rs`), whole, and part by part:
/// the rules without the specials, and the specials of each format built alone.
#[test]
fn the_builder_makes_the_guide_of_the_old_builder() {
    let old = old::AccGuide::built_in();
    let all = AccGuide::built_in(Part::All);
    let expected = canonical_old(&old);
    let whole = canonical_new(&all, |_, submap| {
        submap
            .specials
            .get()
            .expect("made by a build of every part")
    });
    assert_eq!(old.load_warnings, all.load_warnings);
    assert_eq!(old.general, all.general);
    assert_eq!(whole.len(), expected.len());
    for (new, old) in whole.iter().zip(&expected) {
        assert!(new == old, "format {:#x} differs", old.0);
    }
    println!("guide fingerprint {:#018x}", fingerprint(&expected));
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
            // A format of `special` lines alone: no lookup reaches its specials.
            assert!(
                entry.1.is_empty() && entry.2.is_empty(),
                "format {:#x}",
                entry.0
            );
        }
    }
    assert!(parts.next().is_none());
}

/// The lookups of the guide (`Find`, with the key it reports) answer as dev's builder's guide
/// does: at both ends of every range of every small special and next to them, at both ends of
/// every big special and next to them, and for random accessions of every prefix and wildcard.
#[test]
fn lookups_answer_as_with_the_old_builder() {
    let old = old::AccGuide::built_in();
    let new = guide();
    let mut lookups = 0usize;
    let mut check = |fmt: u32, acc: &[u8]| {
        let (mut old_key, mut new_key) = (Vec::new(), Vec::new());
        let expected = old.find(fmt, acc, Some(&mut old_key));
        let got = new.find(fmt, acc, Some(&mut new_key));
        assert_eq!(
            (got, &new_key),
            (expected, &old_key),
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
    for (&fmt, submap) in &old.rules {
        let letters = (fmt >> 16) as usize;
        let digits = (fmt & 0xffff) as usize;
        let size = letters + digits;
        for (key, options) in &submap.small_specials {
            // The digits that the bit is cut from (`x_SplitSpecial` of an accession of the
            // format's size).
            let Some(width) = size
                .checked_sub(key.len())
                .filter(|&w| (1..=9).contains(&w))
            else {
                continue;
            };
            let limit = 10u64.pow(width as u32);
            for (bits, _) in options {
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
        for (to, (from, _)) in &submap.big_specials {
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

/// The ranges of a `special2` block added from the text of the accession, range by range
/// (NCBI's loop, which `add_special2_ranges` follows for a range that changes the key).
fn add_ranges_from_text(guide: &mut AccGuide, hints: &mut Hints, ranges: &[&[u8]]) {
    let fmt = hints.special2_submap.expect("a special2 block");
    let value = hints.special2_type;
    let mut acc = std::mem::take(&mut hints.prev_special2_acc);
    for range in ranges {
        let (low, high) = match range.iter().position(|&c| c == b'-') {
            Some(dash) => (&range[..dash], &range[dash + 1..]),
            None => (&range[..], &b""[..]),
        };
        let end = acc.len();
        acc[end - low.len()..].copy_from_slice(low);
        let from = acc.clone();
        if !high.is_empty() {
            acc[end - high.len()..].copy_from_slice(high);
        }
        AccGuide::add_special(guide.specials_mut(fmt), hints, fmt, &from, &acc, value);
    }
    hints.prev_special2_acc = acc;
}

/// Random digits for one side of a range: mostly at most `bit` of them, sometimes up to `size`
/// with a letter first.
fn random_range_part(next: &mut impl FnMut(u64) -> u64, bit: usize, size: usize) -> Vec<u8> {
    let long = next(10) == 0;
    let max = if long { size } else { bit };
    let len = next(1 + max as u64) as usize;
    let mut out: Vec<u8> = (0..len).map(|_| b'0' + next(10) as u8).collect();
    if long && !out.is_empty() && next(3) == 0 {
        out[0] = b'A' + next(26) as u8;
    }
    out
}

/// `add_special2_ranges` adds random ranges as they are added from the text, in formats of
/// every split of `x_SplitSpecial` (no digits, fewer than 10, 10, more than 10), with ranges
/// that keep the key and ranges that do not (too long, with a second `-`, with letters), and
/// blocks of the same key and of others one after the other.
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
        for block in 0..60 {
            let mut letters_now = spec.clone().into_bytes();
            if block % 4 == 3 {
                letters_now[0] = b'Q';
            }
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
                    let mut token = random_range_part(&mut next, bit, size);
                    if next(4) != 0 {
                        let mut high = random_range_part(&mut next, bit, size);
                        if next(40) == 0 {
                            high.push(b'-');
                            high.extend(random_range_part(&mut next, bit, size));
                        }
                        high.truncate(size);
                        token.push(b'-');
                        token.extend(high);
                    }
                    tokens.push(token);
                }
                let ranges: Vec<&[u8]> = tokens.iter().map(Vec::as_slice).collect();
                numeric.add_special2_ranges(ranges.iter().copied(), &mut numeric_hints);
                add_ranges_from_text(&mut text, &mut text_hints, &ranges);
                assert_eq!(
                    String::from_utf8_lossy(&numeric_hints.prev_special2_acc),
                    String::from_utf8_lossy(&text_hints.prev_special2_acc)
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
