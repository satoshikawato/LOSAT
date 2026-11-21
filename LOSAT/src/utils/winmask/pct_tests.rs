use super::*;

fn lines(text: &[&str]) -> Vec<Vec<u8>> {
    text.iter().map(|line| line.as_bytes().to_vec()).collect()
}

#[test]
fn the_first_value_at_or_above_the_percentage_gives_the_count() {
    let map = [0.0, 50.0, 90.0, 99.5, 100.0];
    assert_eq!(count_by_pct(0.0, &map), 0);
    assert_eq!(count_by_pct(50.0, &map), 1);
    assert_eq!(count_by_pct(90.1, &map), 3);
    assert_eq!(count_by_pct(100.0, &map), 4);
    assert_eq!(count_by_pct(100.5, &[0.0, 99.0]), 1, "none: the last index");
}

#[test]
fn text_formats_read_every_line_and_binary_ones_the_first() {
    let md = lines(&[
        "##note:x",
        "##pct: 0 0",
        "##pct: 1 89.7027",
        "##pct: 2 1e-05 extra",
    ]);
    assert_eq!(
        count_map(&md, StatType::Ascii).unwrap(),
        vec![0.0, 89.7027, 1e-05]
    );
    let md = lines(&["##pct: 2 0 50.5 100", "##pct: 1 7 8"]);
    assert_eq!(
        count_map(&md, StatType::OBinary).unwrap(),
        vec![0.0, 50.5, 100.0]
    );
    assert!(count_map(&lines(&["##note:x"]), StatType::Binary)
        .unwrap()
        .is_empty());
}

#[test]
fn failed_and_missing_values_are_0_as_the_stream_reads_them() {
    let md = lines(&[
        "##pct: 1",
        "##pct: abc 50",
        "##pct: 3 xyz",
        "##pct: -4 12.5",
        "##pct: 99999999999 13",
        "##pct: 6 99.5 extra words",
        "##pct:7 100",
        "##pct: 8 5g",
        "##pct: 9 7,5",
        "##pct: 10 1e",
        "##pct: 11 .5",
        "##pct: 0x10 5",
        "##pct: 13 -0",
        "##pct: 14 +1.5e+1",
    ]);
    let map = count_map(&md, StatType::Ascii).unwrap();
    let want = [
        0.0, 0.0, 0.0, 12.5, 0.0, 99.5, 0.0, 5.0, 7.0, 0.0, 0.5, 0.0, -0.0, 15.0,
    ];
    assert_eq!(map, want);
    assert!(map[12].is_sign_negative());
    let binary = |line: &str| count_map(&lines(&[line]), StatType::Binary).unwrap();
    assert_eq!(binary("##pct: 3 0 1 2"), vec![0.0, 1.0, 2.0, 0.0]);
    assert_eq!(
        binary("##pct: 5 0 0.1 abc 0.3"),
        vec![0.0, 0.1, 0.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(binary("##pct: x 1"), vec![0.0]);
    assert_eq!(binary("##pct: 2"), vec![0.0; 3]);
    assert_eq!(binary("##pct: 1 5 6 7"), vec![5.0, 6.0]);
}

#[test]
fn words_that_c_plus_plus_libraries_read_differently_are_refused() {
    for line in [
        "##pct: 1 -nan",
        "##pct: 1 INF",
        "##pct: 1 infinity",
        "##pct: 1 0x10",
        "##pct: 1 5x",
        "##pct: 1 2.5f",
        "##pct: 1 1e5e",
        "##pct: 1 1e999",
        "##pct: 1 1e-320",
    ] {
        assert!(
            matches!(
                count_map(&lines(&[line]), StatType::Ascii),
                Err(WinMaskError::PctLine { .. })
            ),
            "{line}"
        );
    }
    // Not reached: the count failed, so the percentage is not read.
    assert_eq!(
        count_map(&lines(&["##pct: x inf"]), StatType::Ascii).unwrap(),
        vec![0.0]
    );
    for line in ["##pct: 4294967295 1", "##pct: 99999999999"] {
        assert!(
            matches!(
                count_map(&lines(&[line]), StatType::OBinary),
                Err(WinMaskError::PctLine { .. })
            ),
            "{line}"
        );
    }
}

#[test]
fn a_file_without_percentages_names_the_first_setting_given() {
    let mut posts = Vec::new();
    let pct = PctSettings {
        extend_pct: Some(99.0),
        thres_pct: Some(99.5),
        ..PctSettings::default()
    };
    let params = MaskerParams::default();
    let got = apply_pct(&params, &[], StatType::Ascii, &pct, &mut |w| posts.push(w)).unwrap();
    assert_eq!(got, params);
    assert_eq!(
        posts,
        vec![
            "counts file contains no percentage information; extend_pct settings will not be applied"
                .to_string()
        ]
    );
}
