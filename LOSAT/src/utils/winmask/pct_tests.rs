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
fn words_the_runtime_would_parse_are_refused() {
    for line in [
        "##pct: 1 -nan",
        "##pct: -1 5",
        "##pct: 1 5x",
        "##pct: 1 1e999",
        "##pct: 1 1e-320",
        "##pct: 1",
        "##pct: 99999999999 1",
    ] {
        assert!(
            matches!(
                count_map(&lines(&[line]), StatType::Ascii),
                Err(WinMaskError::PctLine { .. })
            ),
            "{line}"
        );
    }
    assert!(count_map(&lines(&["##pct: 3 0 1 2"]), StatType::Binary).is_err());
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
