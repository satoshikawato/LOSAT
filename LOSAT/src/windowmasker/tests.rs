use super::*;

fn run_words(line: &[&str]) -> (i32, String, String) {
    let words: Vec<OsString> = line.iter().map(OsString::from).collect();
    let mut out = Vec::new();
    let mut err = Vec::new();
    let status = run_with(&words, &mut std::io::empty(), &mut out, &mut err);
    (
        status,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

// NCBI BLAST+ 2.17.0 `windowmasker -help` (5740 bytes) and `-h` (994 bytes, the usage, the
// description and the pointer to `-help`; probe b1 of line NP M5).
#[test]
fn the_short_usage_is_the_head_of_the_help() {
    assert_eq!(HELP.len(), 5740);
    let usage = crate::ncbi_app::short_usage(HELP);
    assert_eq!(usage.len(), 994);
    assert!(usage.starts_with("USAGE\n  windowmasker [-h] [-help] [-xmlhelp] [-ustat unit_counts]"));
}

#[test]
fn the_version_is_losat_s() {
    // LOSAT's version where NCBI's windowmasker states its own `1.0.0` (AGENTS.md rule 8).
    assert_eq!(
        run_words(&["-version"]),
        (
            0,
            format!(
                "windowmasker: LOSAT {}\n Package: LOSAT {}, output compatible with blast {}\n",
                env!("CARGO_PKG_VERSION"),
                env!("CARGO_PKG_VERSION"),
                crate::NCBI_BLAST_RELEASE
            ),
            String::new()
        )
    );
}

#[test]
fn the_library_refuses_xmlhelp_as_blastdbcmd_does() {
    assert_eq!(
        run_words(&["-xmlhelp"]),
        (
            1,
            String::new(),
            format!("Error: {}\n", crate::cli::xml_help_refusal("WINDOWMASKER"))
        )
    );
}

// The making and the conversion of counts are refused after the arguments are read.
#[test]
fn mk_counts_and_convert_are_refused() {
    for (word, what) in [
        (
            "-mk_counts",
            "-mk_counts (the making of a unit counts file)",
        ),
        (
            "-convert",
            "-convert (the conversion of a unit counts file)",
        ),
    ] {
        assert_eq!(
            run_words(&[word]),
            (
                1,
                String::new(),
                format!("Error: {what} is not supported by LOSAT's WINDOWMASKER\n")
            ),
            "{word}"
        );
    }
}

// The reader's view of the keys: the eBoolean defaults (`F`, `T`, `false`), the doubles and the
// integers that a constraint allows.
#[test]
fn the_keys_read_as_ncbi_s() {
    let words =
        |line: &[&str]| -> Vec<String> { line.iter().map(|word| word.to_string()).collect() };
    let args = crate::ncbi_app::args::parse(
        &TOOL,
        &words(&[
            "-ustat",
            "u",
            "-t_thres_pct",
            "99.5",
            "-window",
            "0x14",
            "-dust",
            "yes",
        ]),
    )
    .unwrap();
    assert_eq!(args.boolean("dust"), Some(true));
    assert_eq!(args.boolean("text_match"), Some(true));
    // `-ustat` excludes `-checkdup`: NCBI adds a `CArg_ExcludedValue`, not the default.
    assert_eq!(args.boolean("checkdup"), None);
    assert_eq!(args.double("t_thres_pct"), Some(99.5));
    assert_eq!(args.double("t_low_pct"), None);
    assert_eq!(args.integer("window"), Some(20));
    assert_eq!(args.integer("dust_level"), Some(20));
    assert_eq!(args.value("outfmt"), Some("interval"));
    assert!(crate::ncbi_app::args::parse(&TOOL, &words(&["-window", "0"])).is_err());
    assert!(crate::ncbi_app::args::parse(&TOOL, &words(&["-t_low_pct", "100.5"])).is_err());
    assert!(crate::ncbi_app::args::parse(&TOOL, &words(&["-dust", "maybe"])).is_err());
}

// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:213-254
// The dependencies are NCBI's 38 `SetDependency` calls, all exclusions, as its `multimap` holds
// them (the XML's): each in both directions, ordered by the first name. Which of two exclusions
// NCBI reports is checked against NCBI's outputs (`tests/windowmasker/ncbi/x_*.err`).
#[test]
fn the_dependencies_are_ncbi_s() {
    use crate::ncbi_app::args::Dep;
    let list = crate::ncbi_app::args::dependencies("windowmasker");
    assert_eq!(list.len(), 2 * 38);
    for &(first, dep, second) in &list {
        assert_eq!(dep, Dep::Excludes);
        assert!(list.contains(&(second, dep, first)), "{first} {second}");
    }
    assert!(list.windows(2).all(|pair| pair[0].0 <= pair[1].0));
    assert_eq!(list[0], ("checkdup", Dep::Excludes, "ustat"));
}
