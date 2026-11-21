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

// What the counts modes refuse (line NP item M5b): the percentages of an input without a unit
// of unambiguous letters (NCBI prints 0/0 as its C runtime does), `-infmt seqids` (sequences
// fetched over the network), and the optimizer's search without memory for units of up to 3
// letters (undefined behaviour in C++). Each after NCBI's messages of what came before.
#[test]
fn the_counts_modes_refuse_what_ncbi_gives_no_one_result_for() {
    let refusal = |what: &str| format!("Error: {what} is not supported by LOSAT's WINDOWMASKER\n");
    let warning =
        "FASTA-Reader: Start of first data line in seq is about 100% ambiguous nucleotides \
                   (shouldn't be over 40%)\n";
    let nan = refusal(
        "-mk_counts of an input without a unit of unambiguous letters, whose percentages NCBI \
         BLAST+ prints as the C runtime prints 0/0 (-nan on x86-64 Linux)",
    );
    assert_eq!(
        run_words(&["-mk_counts", "-in", "tests/windowmasker/in/alln.fa"]),
        (
            1,
            String::new(),
            format!("computing the genome length\n{warning}pass 1\n{warning}{nan}")
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
        // The binary format's extra word of the standard output, written when its writer is
        // made, as NCBI's is.
        (1, "\0\0\0\0".to_string(), format!("pass 1\n{nan}"))
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

// What -convert refuses (review of line NP item M5b, 2026-11-21): the optimized formats of units
// of 0 letters (NCBI's search for its table size never ends), units of 17 letters or more (NCBI's
// reverse complement shifts past 32 bits), a ##pct: word that C++ libraries read differently; and
// -mk_counts -infmt blastdb with a -fa_list line that is not UTF-8 (LOSAT opens databases by
// names of text).
#[test]
fn the_convert_mode_refuses_what_ncbi_gives_no_one_result_for() {
    let refusal = |what: &str| format!("Error: {what} is not supported by LOSAT's WINDOWMASKER\n");
    let optimized = refusal(
        "the optimized formats (oascii, obinary) of units of 0 letters, where NCBI BLAST+'s \
         search for the size of its hash table never ends (undefined behaviour in C++)",
    );
    for sformat in ["oascii", "obinary"] {
        let (status, _, err) = run_words(&[
            "-convert",
            "-in",
            "tests/windowmasker/counts/u0.ascii",
            "-sformat",
            sformat,
        ]);
        assert_eq!(status, 1, "{sformat}");
        assert!(
            err.ends_with(&format!("optimizing the data structure\n{optimized}")),
            "{sformat}: {err}"
        );
    }
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-windowmasker-convert-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let body = "0 5\n\n>t_threshold 6\n>t_extend 5\n>t_low 1\n>t_high 100\n";
    let u17 = dir.join("u17.ascii");
    std::fs::write(&u17, format!("17\n{body}")).unwrap();
    let u17 = u17.to_str().unwrap();
    assert_eq!(
        run_words(&["-convert", "-in", u17, "-sformat", "ascii"]),
        (
            1,
            String::new(),
            format!(
                "reading counts...\n{}",
                refusal(&format!(
                    "-convert of -in {u17}, whose units of 17 letters NCBI BLAST+ converts \
                     without an error while its reverse complement of a unit shifts the 32-bit \
                     unit by 32 bits or more (undefined behaviour in C++),"
                ))
            )
        )
    );
    let inf = dir.join("inf.ascii");
    std::fs::write(&inf, format!("##pct: 0 inf\n4\n{body}")).unwrap();
    let inf = inf.to_str().unwrap();
    assert_eq!(
        run_words(&["-convert", "-in", inf, "-sformat", "ascii"]),
        (
            1,
            String::new(),
            format!(
                "reading counts...\n{}",
                refusal(&format!(
                    "-convert of -in {inf}: the percentage line '##pct: 0 inf' of the counts \
                     file holds a number that C++ libraries read differently (an infinity, a \
                     NaN, a hexadecimal number, a number with letters after it, or one beyond \
                     the range of a normal double), so LOSAT does not read it;"
                ))
            )
        )
    );
    let list = dir.join("list.txt");
    std::fs::write(&list, b"tests/windowmasker/db/wm\n\xffwm\n").unwrap();
    assert_eq!(
        run_words(&[
            "-mk_counts",
            "-infmt",
            "blastdb",
            "-fa_list",
            "T",
            "-in",
            list.to_str().unwrap()
        ]),
        (
            1,
            String::new(),
            refusal("-infmt blastdb with a -fa_list line (a database name) that is not UTF-8 text")
        )
    );
    let _ = std::fs::remove_dir_all(&dir);
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
