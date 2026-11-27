use super::*;

fn s(bytes: Vec<u8>) -> String {
    String::from_utf8(bytes).expect("ASCII")
}

#[test]
fn strings_observed_from_ncbi() {
    // blastn megablast / blastn / dc-megablast: DUST at its defaults, mask at hash.
    assert_eq!(
        s(filter_string(Some((20, 64, 1)), None, None, None, true)),
        "L;m;"
    );
    // blastn-short and `-dust no`: only mask at hash.
    assert_eq!(s(filter_string(None, None, None, None, true)), "m;");
    // `-dust no -soft_masking false` and blastp without SEG: empty.
    assert_eq!(s(filter_string(None, None, None, None, false)), "F");
    // blastp `-seg yes -soft_masking true`; blastx / tblastn / tblastx default SEG.
    assert_eq!(
        s(filter_string(None, Some((12, 2.2, 2.5)), None, None, true)),
        "L;m;"
    );
    assert_eq!(
        s(filter_string(None, Some((12, 2.2, 2.5)), None, None, false)),
        "L;"
    );
    // `-seg '10 1.8 2.1'`
    assert_eq!(
        s(filter_string(None, Some((10, 1.8, 2.1)), None, None, false)),
        "S 10 1.8 2.1;"
    );
}

#[test]
fn custom_dust_and_seg_and_masker() {
    assert_eq!(
        s(filter_string(Some((25, 64, 1)), None, None, None, true)),
        "D 25 64 1;m;"
    );
    assert_eq!(
        s(filter_string(
            Some((20, 64, 1)),
            Some((12, 2.2, 2.5)),
            None,
            None,
            false
        )),
        "L;L;"
    );
    assert_eq!(
        s(filter_string(
            None,
            None,
            None,
            Some(WindowMasker::Taxid(9606)),
            true
        )),
        "W -t 9606;m;"
    );
    assert_eq!(
        s(filter_string(
            None,
            None,
            None,
            Some(WindowMasker::Database(b"db.wm")),
            false
        )),
        "W -d db.wm;"
    );
    assert_eq!(
        s(filter_string(
            None,
            None,
            None,
            Some(WindowMasker::Taxid(0)),
            false
        )),
        "F"
    );
}

#[test]
fn seg_cutoffs_are_printed_as_glibc_rounds_them() {
    // %1.1f of an exact tie goes to the even digit: 2.25 -> 2.2, 0.35 (below the tie in
    // binary) -> 0.3, 0.45 (above) -> 0.5.
    assert_eq!(
        s(filter_string(
            None,
            Some((10, 2.25, 0.35)),
            None,
            None,
            false
        )),
        "S 10 2.2 0.3;"
    );
    assert_eq!(
        s(filter_string(
            None,
            Some((10, 0.45, 0.125)),
            None,
            None,
            false
        )),
        "S 10 0.5 0.1;"
    );
}

// NCBI 2.17.0 blastp `-outfmt 5` with `-seg` cuts that are not finite (sw8/nfix/probe/seg2).
#[test]
fn infinite_and_nan_seg_cutoffs_are_printed_as_glibc_prints_them() {
    let seg = |locut: f64, hicut: f64| {
        s(filter_string(
            None,
            Some((12, locut, hicut)),
            None,
            None,
            false,
        ))
    };
    assert_eq!(seg(f64::INFINITY, 2.5), "S 12 inf 2.5;");
    assert_eq!(seg(2.2, -f64::NAN), "S 12 2.2 -nan;");
    assert_eq!(
        seg(f64::from_bits(0x7ff8_0000_0000_0001), 2.5),
        "S 12 nan 2.5;"
    );
    assert_eq!(seg(f64::NEG_INFINITY, 2.5), "S 12 -inf 2.5;");
    assert_eq!(seg(-f64::NAN, -f64::NAN), "S 12 -nan -nan;");
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/blastfilter_unit_test.cpp:1683-1697 FilterOptionsToStringLargeData
#[test]
fn a_long_repeat_database_is_written_whole() {
    // The database of blastn's `-filtering_db` is appended to the text whole, not through a
    // buffer of fixed size, after a custom DUST.
    let database = "X".repeat(4096);
    let text = s(filter_string(
        Some((20, 128, 1)),
        None,
        Some(database.as_bytes()),
        None,
        false,
    ));
    assert!(text.starts_with("D 20 128 1;R -d XXXXXXXXXXXXXXXXXXXX"));
    assert_eq!(text, format!("D 20 128 1;R -d {database};"));
}

#[test]
fn a_long_custom_string_is_cut_like_the_24_byte_buffer() {
    let text = s(filter_string(
        Some((2147483647, 2147483647, 2147483647)),
        None,
        None,
        None,
        false,
    ));
    assert_eq!(text, "D 2147483647 2147483647");
    assert_eq!(text.len(), 23);
}
