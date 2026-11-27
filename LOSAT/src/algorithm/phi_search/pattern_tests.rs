use super::*;

const ALL: i32 = MASK_AA_ALPHABET_BITS;
const A: i32 = 1 << 1;
const G: i32 = 1 << 7;
const K: i32 = 1 << 10;
const L: i32 = 1 << 11;
const S: i32 = 1 << 17;
const T: i32 = 1 << 18;

fn compiled(text: &str) -> PhiPattern {
    PhiPattern::compile(text.as_bytes()).unwrap()
}

fn ncbi(text: &[u8]) -> String {
    match PhiPattern::compile(text) {
        Err(PhiError::Ncbi(message)) => message,
        other => panic!("{text:?}: {other:?}"),
    }
}

fn refused(text: &[u8]) -> String {
    match PhiPattern::compile(text) {
        Err(PhiError::Refusal(message)) => message,
        other => panic!("{text:?}: {other:?}"),
    }
}

#[test]
fn ploop_is_one_word_with_ncbi_probability() {
    let blk = compiled("[ag]-x(4)-G-K-[ST]");
    assert_eq!(blk.kind, PatternKind::OneWord);
    assert_eq!(blk.text, b"[AG]-X(4)-G-K-[ST]".to_vec());
    assert_eq!(blk.min_match_length, 8);
    assert_eq!(&blk.masked[..8], &[A | G, ALL, ALL, ALL, ALL, G, K, S | T]);
    assert_eq!(blk.match_mask, 1 << 7);
    assert_eq!(blk.which_positions[1], 0b1_1111);
    assert_eq!(blk.which_positions[7], 0b11_1111);
    assert_eq!(blk.which_positions[10], 0b101_1110);
    assert_eq!(blk.which_positions[17], 0b1001_1110);
    assert_eq!(blk.which_positions[0], 0b1_1110);
    // NCBI prints `pattern probability=8.338e-05` (probe rd.base.tk.db.0).
    assert!(
        (blk.probability - 8.338e-5).abs() < 6e-10,
        "{}",
        blk.probability
    );
}

#[test]
fn variable_wildcard_expands_in_ncbi_order_with_stale_flags() {
    let blk = compiled("G-x(1,3)-K");
    assert_eq!(blk.kind, PatternKind::OneWord);
    assert_eq!(blk.min_match_length, 3);
    // The one-wildcard form in place, then no wildcard, then two.
    assert_eq!(
        &blk.masked[..12],
        &[G, ALL, ALL, K, G, ALL, K, G, ALL, ALL, ALL, K]
    );
    assert_eq!(blk.match_mask, (1 << 3) | (1 << 6) | (1 << 11));
    assert_eq!(blk.which_positions[7], 0b111_1011_0111);
    assert_eq!(blk.which_positions[10], 0b1111_0110_1110);
    let prob = crate::stats::karlin_calc::compute_std_aa_composition();
    assert_eq!(blk.probability, prob[7] * 2.0 * prob[10]);
}

#[test]
fn leading_variable_region_is_shifted_out() {
    let blk = compiled("x(2)-x(0,3)-G-K-[ST]");
    assert_eq!(&blk.masked[..5], &[ALL, ALL, G, K, S | T]);
    assert_eq!(blk.match_mask, 1 << 4);
    assert_eq!(blk.min_match_length, 5);
}

#[test]
fn reversed_and_negative_counts() {
    // `x(5,3)`: the entry -(-2) = 2 is the residue A.
    let blk = compiled("G-x(5,3)-K");
    assert_eq!(&blk.masked[..8], &[G, ALL, ALL, ALL, ALL, ALL, A, K]);
    // NCBI prints `pattern probability=-0.0045032` (el.negprob).
    let blk = compiled("G-x(7,0)-[AG]-K");
    assert!(
        (blk.probability + 0.0045032).abs() < 6e-9,
        "{}",
        blk.probability
    );
    // A negative repeat count writes nothing and lowers the minimum length.
    let blk = compiled("G-A(-5)");
    assert_eq!(blk.min_match_length, -4);
    assert_eq!(blk.match_mask, 1);
}

#[test]
fn patterns_without_positions_match_nothing() {
    for text in ["-", "...", "x(0)", "A(0)", "x(0,3)", "\r", "\nG"] {
        let blk = compiled(text);
        assert_eq!(blk.kind, PatternKind::OneWord, "{text:?}");
        assert_eq!(blk.match_mask, 0, "{text:?}");
    }
    // `[]`: one position that no residue matches.
    let blk = compiled("[]");
    assert_eq!(blk.match_mask, 1);
    assert_eq!(blk.which_positions, [0; BLASTAA_SIZE]);
    assert_eq!(blk.probability, 0.0);
}

#[test]
fn digits_star_and_u_are_residues() {
    let blk = compiled("1-*-U");
    assert_eq!(&blk.masked[..3], &[1, 1 << 25, ALL * 2 + 1]);
}

#[test]
fn multi_word_and_very_long_kinds() {
    let blk = compiled("[AG]-x(4)-G-K-[ST]-x(23)-[LD]");
    assert_eq!(blk.kind, PatternKind::MultiWord);
    assert_eq!(blk.num_words, 2);
    assert_eq!(&blk.match_mask_l[..2], &[0, 1 << 1]);
    assert_eq!(blk.bit_pattern_by_letter[11][1], 0b11);
    assert_eq!(
        blk.bit_pattern_by_letter[7][0],
        (1 << 30) - 1 - (1 << 6) - (1 << 7)
    );

    let blk = compiled("[AG]-x(4)-G-K-[ST]-x(0,25)-L");
    assert_eq!(blk.kind, PatternKind::VeryLong);
    assert_eq!(blk.num_words, 2);
    assert_eq!(blk.highest_place, 10);
    assert_eq!(
        &blk.masked[..10],
        &[A | G, ALL, ALL, ALL, ALL, G, K, S | T, -25, L]
    );
    assert_eq!(&blk.num_places_in_word[..2], &[8, 1]);
    assert_eq!(&blk.spacing[..2], &[25, 0]);
    assert_eq!(&blk.match_mask_l[..2], &[1 << 7, 1]);
    assert_eq!(blk.most_specific, 0);
    assert_eq!(blk.sll[1][11], 1);
    assert_eq!(blk.sll[0][7], 0b11_1111);
}

#[test]
fn ncbi_messages() {
    assert_eq!(
        ncbi(&[b'A'; 330]),
        "Pattern is too long (330 but only 330 supported)"
    );
    assert_eq!(ncbi(b"A(330)"), "Pattern is too long");
    assert_eq!(ncbi(b"M-x(328)-K"), "Pattern is too long");
    assert_eq!(ncbi(b"x(330,331)"), "Pattern too long");
    let bracket = "pattern description has a non-alphabeticcharacter inside a bracket";
    assert_eq!(ncbi(b"[A"), bracket);
    assert_eq!(ncbi(b"[A1]"), bracket);
    assert_eq!(ncbi(b"[A\xc3\xa9]"), bracket);
}

#[test]
fn refusals() {
    let says = |text: &[u8], what: &str, why: &str| {
        assert_eq!(
            refused(text),
            format!(
                "a -phi_pattern with {what} (NCBI BLAST+ 2.17.0 {why}) is not supported by LOSAT's PSIBLAST"
            ),
            "{text:?}"
        );
    };
    for text in [&b"x("[..], b"x(3"] {
        says(text, "an unterminated x(", "crashes");
    }
    for text in [&b"x(3,5"[..], b"A(3", b"A()"] {
        says(text, "an unterminated (", "crashes");
    }
    says(b"{AG", "an unterminated {", "crashes");
    for text in [&b"x(a,b)"[..], b"x(3,)", b"x(,3)", b"x(3 ,5)"] {
        says(
            text,
            "x(n,m) without two numbers",
            "reads uninitialised values",
        );
    }
    for text in [
        &b"A(2147483648)"[..],
        b"x(-2147483648,1)",
        b"x(0,-2147483648)",
    ] {
        says(
            text,
            "a count outside 32-bit integers",
            "overflows a 32-bit integer",
        );
    }
    for text in [&b"A(331)"[..], b"A(400)", b"M-x(330)", b"x(329,330)-K"] {
        says(
            text,
            "a pattern longer than NCBI's 330 positions",
            "writes or reads past its arrays",
        );
    }
    says(
        b"x(329)x(0,-268435455)",
        "330 wildcard positions",
        "writes or reads past its arrays",
    );
    for text in [&b"G-\xc3"[..], b"{A\xc3}", b"\xc3(2)"] {
        says(
            text,
            "a non-ASCII byte in the pattern",
            "reads past its residue table",
        );
    }
    says(
        b"A-x(0,25)-x(0,25)-L",
        "two variable wildcards in a row in a very long pattern",
        "shifts by a negative count",
    );
}
