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
fn the_version_is_windowmasker_s() {
    assert_eq!(
        run_words(&["-version"]),
        (
            0,
            format!(
                "windowmasker: 1.0.0\n Package: LOSAT {}, output compatible with blast {}\n",
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

// What the counts modes refuse (line NP item M5b): the percentages of an input without a unit
// of unambiguous letters (NCBI prints 0/0 as its C runtime does), `-infmt seqids` (sequences
// fetched over the network), and the optimizer's search without memory for units of up to 3
// letters (undefined behaviour in C++). Each after NCBI's messages of what came before.
#[test]
fn the_counts_modes_refuse_what_ncbi_gives_no_one_result_for() {
    let refusal = |what: &str| format!("Error: {what} is not supported by LOSAT's WINDOWMASKER\n");
    let nan = refusal(
        "-mk_counts of an input without a unit of unambiguous letters, whose percentages NCBI \
         BLAST+ prints as the C runtime prints 0/0 (-nan on x86-64 Linux)",
    );
    assert_eq!(
        run_words(&["-mk_counts", "-in", "tests/windowmasker/in/alln.fa"]),
        (
            1,
            String::new(),
            format!("computing the genome length\npass 1\n{nan}")
        )
    );
    assert_eq!(
        run_words(&[
            "-mk_counts",
            "-in",
            "tests/windowmasker/in/empty.fa",
            "-unit",
            "3",
            "-sformat",
            "binary"
        ]),
        (1, String::new(), format!("pass 1\n{nan}"))
    );
    assert_eq!(
        run_words(&["-mk_counts", "-infmt", "seqids", "-in", "ids.txt"]),
        (
            1,
            String::new(),
            refusal(
                "-mk_counts with -infmt seqids (NCBI BLAST+ fetches each sequence from NCBI over \
                 the network)"
            )
        )
    );
    let smem = refusal(
        "-smem of a multiple of 65536 (whose low 16 bits, 0, NCBI BLAST+ keeps) for units of 3 \
         letters, where NCBI BLAST+ sizes its hash table for a key of 253 bits, allocates 8 GiB \
         and shifts past the width of its integers (undefined behaviour in C++)",
    );
    for sformat in ["obinary", "oascii"] {
        assert_eq!(
            run_words(&[
                "-mk_counts",
                "-in",
                "tests/winmask/fasta/dm2L_3k.fa",
                "-unit",
                "3",
                "-smem",
                "65536",
                "-sformat",
                sformat
            ]),
            (
                1,
                String::new(),
                format!("pass 1\npass 2\noptimizing the data structure\n{smem}")
            ),
            "{sformat}"
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

// The dependencies are NCBI's 38 `SetDependency` calls in order.
#[test]
fn the_dependencies_are_ncbi_s() {
    let list = dependencies();
    assert_eq!(list.len(), 38);
    assert_eq!(list[0], ("mk_counts", Dep::Excludes, "outfmt"));
    assert_eq!(list[10], ("ustat", Dep::Excludes, "checkdup"));
    assert_eq!(list[37], ("convert", Dep::Excludes, "text_match"));
}
