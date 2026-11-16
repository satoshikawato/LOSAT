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
    let guide = guide();
    assert_eq!(guide.load_warnings, Vec::<String>::new());
    // No value is a fallback (every type name is known).
    for submap in guide.rules.values() {
        assert!(submap.prefixes.values().all(|v| v & F_FALLBACK == 0));
        assert!(submap.wildcards.iter().all(|(mask, v)| v & F_FALLBACK == 0
            && mask
                .iter()
                .all(|&c| c.is_ascii_uppercase() || c == b'_' || c == b'?')));
        assert!(submap
            .big_specials
            .values()
            .all(|(_, v)| v & F_FALLBACK == 0));
        assert!(submap
            .small_specials
            .values()
            .flatten()
            .all(|(bits, v)| v & F_FALLBACK == 0 && bits.any()));
    }
    assert_eq!(guide.general.get(&b"SRA"[..]), Some(&GENERAL_NUC));
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
