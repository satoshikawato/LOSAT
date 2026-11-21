use super::*;

fn os(text: &str) -> Option<&OsStr> {
    Some(OsStr::new(text))
}

#[test]
fn a_switch_follows_its_variable() {
    for off in ["0", "0 ", " 0\t"] {
        assert!(!state(os(off)), "{off:?}");
    }
    for on in ["1", "16", "on", "00", "false", "", "  "] {
        assert!(state(os(on)), "{on:?}");
    }
    assert!(state(None));
}

#[test]
fn a_value_switch_takes_a_number_in_range_only() {
    let range = 2..=65_536;
    assert_eq!(given(os("32"), &range), Some(32));
    assert_eq!(given(os("1"), &range), None);
    assert_eq!(given(os("65537"), &range), None);
    assert_eq!(given(os(" 32 "), &range), Some(32));
    assert_eq!(given(os("3 2"), &range), None);
    assert_eq!(given(os("x"), &range), None);
    assert_eq!(given(os(""), &range), None);
    assert_eq!(given(None, &range), None);
}

#[test]
fn the_strict_set_names_are_distinct_and_sorted() {
    assert!(STRICT_SET.windows(2).all(|pair| pair[0] < pair[1]));
    assert!(STRICT_SET.iter().all(|name| name.starts_with("LOSAT_X_")));
}

#[test]
fn the_graduated_names_are_distinct_sorted_and_outside_the_set() {
    assert!(GRADUATED.windows(2).all(|pair| pair[0] < pair[1]));
    assert!(GRADUATED.iter().all(|name| name.starts_with("LOSAT_X_")));
    assert!(GRADUATED.iter().all(|name| !STRICT_SET.contains(name)));
}

/// The block of `docs/agents/diagnostics.md` that writes the set as `export` lines names the
/// switches of `STRICT_SET`, one per line, in the same order.
#[test]
fn the_documented_block_lists_the_strict_set() {
    let doc = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/agents/diagnostics.md"
    ));
    let anchor = doc
        .find("The strict set (the variables that are on;")
        .expect("the block's heading in diagnostics.md");
    let block = &doc[anchor..];
    let start = block.find("```bash\n").expect("the block") + "```bash\n".len();
    let end = start + block[start..].find("```").expect("the end of the block");
    let documented: Vec<&str> = block[start..end]
        .lines()
        .map(|line| {
            let name = line
                .strip_prefix("export ")
                .and_then(|rest| rest.split_once('='))
                .map(|(name, _)| name);
            name.unwrap_or_else(|| panic!("not one `export NAME=value` line: {line:?}"))
        })
        .collect();
    assert_eq!(documented, STRICT_SET);
}
