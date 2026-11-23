use super::*;

#[test]
fn the_pattern_block_lists_the_positions_and_the_probability() {
    let info = PhiQueryInfo {
        text: b"[AG]-X(4)-G-K-[ST]",
        probability: 8.338e-05,
        occurrences: vec![(19, 8), (208, 8)],
    };
    let mut out = Vec::new();
    write_phi_info(&mut out, &info).unwrap();
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "2 occurrence(s) of pattern: \n[AG]-X(4)-G-K-[ST] at position(s) 20, 209 of query \
         sequence\npattern probability=8.338e-05\n"
    );
}

#[test]
fn the_pattern_text_ends_at_a_nul_byte() {
    let info = PhiQueryInfo {
        text: b"G-K\0-S",
        probability: 0.0153169,
        occurrences: vec![(12, 3)],
    };
    let mut out = Vec::new();
    write_phi_info(&mut out, &info).unwrap();
    assert!(String::from_utf8(out)
        .unwrap()
        .starts_with("1 occurrence(s) of pattern: \nG-K at position(s) 13 of"));
}

#[test]
fn an_occurrence_without_alignments_keeps_the_table_header() {
    let mut out = Vec::new();
    write_empty_table(&mut out, 500).unwrap();
    let text = String::from_utf8(out).unwrap();
    assert_eq!(
        text,
        format!(
            "{}Score     E\n{:<69}(Bits)  Value\n\n\n",
            " ".repeat(70),
            "Sequences producing significant alignments:"
        )
    );
    let mut out = Vec::new();
    write_empty_table(&mut out, 0).unwrap();
    assert_eq!(out, b"\n");
}
