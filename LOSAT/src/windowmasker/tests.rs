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
    let usage = crate::ncbi_app::short_usage(&TOOL);
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
// reverse complement shifts past 32 bits) and a ##pct: word that C++ libraries read differently.
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
    let _ = std::fs::remove_dir_all(&dir);
}

// A -fa_list line of -mk_counts -infmt blastdb names a database by its bytes, as NCBI's string
// does (line RG item K10a): a copy of `db_search/xb_n` under the name `\xffxb` counts as the
// original. Not on macOS, whose APFS refuses such a file name (as `tests/non_utf8_names.rs`).
#[cfg(all(unix, not(target_os = "macos")))]
#[test]
fn a_fa_list_database_name_is_read_by_its_bytes() {
    use std::os::unix::ffi::OsStrExt;
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-windowmasker-fa-list-bytes-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    for extension in ["nhr", "nin", "nsq"] {
        let mut name = b"\xffxb.".to_vec();
        name.extend_from_slice(extension.as_bytes());
        std::fs::copy(
            format!("tests/db_search/xb_n.{extension}"),
            dir.join(std::ffi::OsStr::from_bytes(&name)),
        )
        .unwrap();
    }
    let list = dir.join("list.txt");
    let count = |lines: &[u8]| {
        std::fs::write(&list, lines).unwrap();
        run_words(&[
            "-mk_counts",
            "-infmt",
            "blastdb",
            "-fa_list",
            "T",
            "-in",
            list.to_str().unwrap(),
        ])
    };
    let mut by_bytes = b"tests/db_search/xb_n\n".to_vec();
    by_bytes.extend_from_slice(dir.as_os_str().as_bytes());
    by_bytes.extend_from_slice(b"/\xffxb\n");
    let by_bytes = count(&by_bytes);
    assert_eq!(by_bytes.0, 0, "{by_bytes:?}");
    assert_eq!(
        by_bytes,
        count(b"tests/db_search/xb_n\ntests/db_search/xb_n\n")
    );
    let _ = std::fs::remove_dir_all(&dir);
}

// NCBI's counts writers open `-out` by its bytes (`CSeqMaskerOstatFactory::create` of the
// argument's string, line RG item K4i-2a): a name that is not UTF-8 is written as given, in the
// counts mode and in the convert mode, and the counts are those of a name of text. Not on macOS,
// whose APFS refuses a file name that is not UTF-8 (as `tests/non_utf8_names.rs`).
#[cfg(all(unix, not(target_os = "macos")))]
#[test]
fn the_counts_modes_write_an_out_name_that_is_not_utf8() {
    use std::os::unix::ffi::OsStringExt;
    let run = |line: Vec<OsString>| {
        let mut out = Vec::new();
        let mut err = Vec::new();
        let status = run_with(&line, &mut std::io::empty(), &mut out, &mut err);
        (status, out, String::from_utf8(err).unwrap())
    };
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-windowmasker-out-bytes-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let named = |name: &[u8]| {
        let mut path = dir.clone().into_os_string().into_vec();
        path.push(b'/');
        path.extend_from_slice(name);
        OsString::from_vec(path)
    };
    let words = |line: &[&str], out: OsString| -> Vec<OsString> {
        let mut words: Vec<OsString> = line.iter().map(OsString::from).collect();
        words.extend(["-out".into(), out]);
        words
    };
    let counts = ["-mk_counts", "-in", "tests/windowmasker/in/cd.fa"];
    let (status, out, err) = run(words(&counts, named(b"c\xff")));
    assert_eq!((status, out), (0, Vec::new()), "{err}");
    let (status, _, err) = run(words(&counts, named(b"c")));
    assert_eq!(status, 0, "{err}");
    let written = std::fs::read(named(b"c")).unwrap();
    assert!(!written.is_empty());
    assert_eq!(std::fs::read(named(b"c\xff")).unwrap(), written);
    let convert = [
        "-convert",
        "-in",
        "tests/windowmasker/counts/u4.ascii",
        "-sformat",
        "binary",
    ];
    let (status, _, err) = run(words(&convert, named(b"b\xff")));
    assert_eq!(status, 0, "{err}");
    let (status, _, err) = run(words(&convert, named(b"b")));
    assert_eq!(status, 0, "{err}");
    assert_eq!(
        std::fs::read(named(b"b\xff")).unwrap(),
        std::fs::read(named(b"b")).unwrap()
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

// The reader's view of the keys: the eBoolean defaults (`F`, `T`, `false`), the doubles and the
// integers that a constraint allows.
#[test]
fn the_keys_read_as_ncbi_s() {
    let words = |line: &[&str]| -> Vec<OsString> { line.iter().map(OsString::from).collect() };
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
    assert_eq!(
        args.bytes("outfmt").as_deref(),
        Some(b"interval".as_slice())
    );
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

/// A run of LOSAT Web (`run_in`, `RunFiles::Receiver`) gives the bytes of the command line's run in
/// the same directory (here the crate's, `run_with` on the file system): the exit status, the
/// standard output and error, and each file that the run writes (`-out FILE`, `-logfile FILE`),
/// which the receiver gets under its name in place of the disk; the files of `-mk_counts` (binary,
/// with `-t_low_pct`, which a named file alone takes), `-convert` and the masking; a folder that the
/// directory does not have opens no counts file (NCBI's writer then loses the counts) and no masking
/// output (NCBI's `File is not accessible`); `-xmlhelp` prints NCBI's XML, as the command line's run
/// does (line NP WW3a, 2026-11-22).
#[test]
fn a_run_of_losat_web_writes_the_files_of_the_command_line() {
    use crate::report::run_files::{FileReceiver, RunFiles};
    use std::sync::{Arc, Mutex};

    #[derive(Default)]
    struct Kept(Mutex<Vec<(Vec<u8>, Vec<u8>)>>);
    impl FileReceiver for Kept {
        fn receive(&self, name: &[u8], bytes: &[u8]) -> std::io::Result<()> {
            self.0.lock().unwrap().push((name.to_vec(), bytes.to_vec()));
            Ok(())
        }
    }
    let disk = crate::utils::test_dir::temp_dir()
        .join(format!("ww3a-wm-{}", crate::utils::test_dir::process_id()));
    std::fs::create_dir_all(&disk).unwrap();
    let fasta = "tests/winmask/fasta/dm2L_3k.fa";
    let counts = "tests/winmask/counts/dm2L_u9.binary";
    let rows: [&[&str]; 5] = [
        &[
            "-mk_counts",
            "-in",
            fasta,
            "-unit",
            "5",
            "-sformat",
            "binary",
            "-t_low_pct",
            "50",
            "-out",
            "@/c.counts",
            "-logfile",
            "@/run.log",
        ],
        &[
            "-convert",
            "-in",
            "tests/windowmasker/counts/u4.ascii",
            "-sformat",
            "obinary",
            "-out",
            "@/c.obinary",
        ],
        &[
            "-ustat",
            counts,
            "-in",
            fasta,
            "-outfmt",
            "fasta",
            "-out",
            "@/masked.fa",
        ],
        &[
            "-mk_counts",
            "-in",
            fasta,
            "-unit",
            "5",
            "-out",
            "ww3a-no-such-folder/c.counts",
        ],
        &[
            "-ustat",
            counts,
            "-in",
            fasta,
            "-out",
            "ww3a-no-such-folder/masked.fa",
        ],
    ];
    for row in rows {
        // The command line's names under a directory of the disk; LOSAT Web's in the working
        // directory, as the run's own names.
        let named = |prefix: &str| -> Vec<String> {
            row.iter().map(|word| word.replace('@', prefix)).collect()
        };
        let cli = named(disk.to_str().unwrap());
        let web = named(".");
        let cli_words: Vec<OsString> = cli.iter().map(OsString::from).collect();
        let web_words: Vec<OsString> = web.iter().map(OsString::from).collect();
        let (mut cli_out, mut cli_err) = (Vec::new(), Vec::new());
        let cli_status = run_with(
            &cli_words,
            &mut std::io::empty(),
            &mut cli_out,
            &mut cli_err,
        );
        let kept = Arc::new(Kept::default());
        let (mut web_out, mut web_err) = (Vec::new(), Vec::new());
        let web_status = run_in(
            &web_words,
            &RunFiles::Receiver(kept.clone()),
            &mut std::io::empty(),
            &mut web_out,
            &mut web_err,
        );
        assert_eq!(web_status, cli_status, "{row:?}");
        assert_eq!(web_out, cli_out, "{row:?}");
        let disk_prefix = format!("{}/", disk.to_str().unwrap());
        assert_eq!(
            String::from_utf8(web_err).unwrap(),
            String::from_utf8(cli_err)
                .unwrap()
                .replace(&disk_prefix, "./"),
            "{row:?}"
        );
        let mut made: Vec<(Vec<u8>, Vec<u8>)> = std::fs::read_dir(&disk)
            .unwrap()
            .map(|entry| {
                let path = entry.unwrap().path();
                let name = path.file_name().unwrap().to_str().unwrap().to_string();
                let bytes = std::fs::read(&path).unwrap();
                std::fs::remove_file(&path).unwrap();
                (name.into_bytes(), bytes)
            })
            .collect();
        let mut given = kept.0.lock().unwrap().clone();
        made.sort();
        given.sort();
        assert_eq!(given, made, "{row:?}");
        assert_eq!(
            row.iter().filter(|word| word.starts_with('@')).count(),
            made.len(),
            "{row:?}"
        );
    }
    std::fs::remove_dir_all(&disk).unwrap();
    let mut out = Vec::new();
    let status = run_in(
        &[OsString::from("-xmlhelp")],
        &RunFiles::Disk,
        &mut std::io::empty(),
        &mut out,
        &mut Vec::new(),
    );
    assert_eq!(
        (status, out),
        (
            0,
            crate::ncbi_xmlhelp::print("windowmasker", |_| false)
                .unwrap()
                .into_bytes()
        )
    );
}

/// The checks of LOSAT Web (`check_in`, `app::check`) read the values as the command line does in
/// each mode, list and file names that are not UTF-8 included: a name read as text fails here in a
/// debug build (`Args::value`'s check; review F1 of item K10a of line RG). A given list is read by
/// the run before `-infmt seqids` is refused, so the check passes it; without one the check refuses.
#[cfg(unix)]
#[test]
fn the_checks_of_each_mode_read_names_by_their_bytes() {
    use std::os::unix::ffi::OsStringExt;
    let check = |line: &[&str], names: &[(&str, &[u8])]| {
        let mut words: Vec<OsString> = line.iter().map(OsString::from).collect();
        for (key, value) in names {
            words.push(OsString::from(key));
            words.push(OsString::from_vec(value.to_vec()));
        }
        let mut out = Vec::new();
        let mut err = Vec::new();
        let status = check_in(
            &words,
            &crate::report::run_files::RunFiles::Disk,
            &mut out,
            &mut err,
        );
        (status, String::from_utf8_lossy(&err).into_owned())
    };
    let name: &[u8] = b"l\xe9";
    for list in ["-ids", "-exclude_ids"] {
        for line in [
            &["-mk_counts", "-in", "missing.fa"][..],
            &["-mk_counts", "-in", "missing.fa", "-infmt", "seqids"][..],
            &["-ustat", "missing.counts", "-in", "missing.fa"][..],
        ] {
            assert_eq!(check(line, &[(list, name)]).0, 0, "{line:?} {list}");
        }
    }
    assert_ne!(
        check(
            &["-mk_counts", "-in", "missing.fa", "-infmt", "seqids"],
            &[]
        )
        .0,
        0
    );
    assert_eq!(check(&["-mk_counts"], &[("-in", name)]).0, 0);
    assert_eq!(check(&["-convert", "-out", "o"], &[("-in", name)]).0, 0);
    assert_eq!(
        check(&["-in", "missing.fa"], &[("-ustat", name), ("-meta", name)]).0,
        0
    );
}

/// A `windowmasker` receiver run reads `-logfile` against the run's real files (D19): missing
/// folders and existing folders keep diagnostics on stderr with no file; an existing file is
/// refused (exit 1); `/dev/stdout` writes diagnostic messages into `out` in chronological order;
/// and normalized file names are delivered (` run.log` -> `run.log`, `logs/` -> `logs`).
#[test]
fn a_receiver_windowmasker_reads_logfile_against_the_runs_files() {
    use crate::report::run_files::{FileReceiver, RunFiles};
    use std::sync::{Arc, Mutex};

    #[derive(Default)]
    struct Kept(Mutex<Vec<(Vec<u8>, Vec<u8>)>>);
    impl FileReceiver for Kept {
        fn receive(&self, name: &[u8], bytes: &[u8]) -> std::io::Result<()> {
            self.0.lock().unwrap().push((name.to_vec(), bytes.to_vec()));
            Ok(())
        }
    }
    let fasta = "tests/winmask/fasta/dm2L_3k.fa";
    let run_web = |words: &[&str]| {
        let kept = Arc::new(Kept::default());
        let args: Vec<OsString> = words.iter().map(OsString::from).collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let status = run_in(
            &args,
            &RunFiles::Receiver(kept.clone()),
            &mut std::io::empty(),
            &mut out,
            &mut err,
        );
        let files = kept.0.lock().unwrap().clone();
        (status, out, err, files)
    };
    // Case 0 & 1: missing folder and existing folder (`src`) keep messages on stderr, no file.
    for log in ["d19-no-such-folder/run.log", "src", "src/"] {
        let (status, out, err, files) = run_web(&[
            "-mk_counts",
            "-in",
            fasta,
            "-unit",
            "5",
            "-out",
            "c.out",
            "-logfile",
            log,
        ]);
        assert_eq!(status, 0, "{log}");
        assert!(out.is_empty(), "{log}");
        assert_eq!(String::from_utf8(err).unwrap(), "pass 1\npass 2\n", "{log}");
        assert_eq!(files.len(), 1, "{log}: only -out file");
        assert_eq!(files[0].0, b"c.out", "{log}");
    }
    // Case 2: existing file of the run's directory is refused, exit status 1.
    let (status, out, err, files) = run_web(&["-mk_counts", "-in", fasta, "-logfile", fasta]);
    assert_eq!(status, 1);
    assert!(out.is_empty());
    assert!(files.is_empty());
    assert_eq!(
        String::from_utf8(err).unwrap(),
        format!(
            "Error: a -logfile file '{fasta}' that is a file of the run's directory (the command \
             line appends its messages to it; LOSAT Web reads the run's files and writes none of \
             them) is not supported by LOSAT's WINDOWMASKER\n"
        )
    );
    // Case 3: `-logfile /dev/stdout` writes messages into `out` in chronological order before counts.
    let (status_plain, out_plain, err_plain, _) =
        run_web(&["-mk_counts", "-in", fasta, "-unit", "5"]);
    assert_eq!(status_plain, 0);
    assert_eq!(err_plain, b"pass 1\npass 2\n");
    for extra in [&[][..], &["-out", "/dev/stdout"][..]] {
        let mut words = vec!["-mk_counts", "-in", fasta, "-unit", "5"];
        words.extend_from_slice(extra);
        words.extend(["-logfile", "/dev/stdout"]);
        let (status, out, err, files) = run_web(&words);
        assert_eq!(status, 0, "{extra:?}");
        assert!(err.is_empty(), "{extra:?}");
        assert!(files.is_empty(), "{extra:?}");
        let mut want = b"pass 1\npass 2\n".to_vec();
        want.extend_from_slice(&out_plain);
        assert_eq!(out, want, "{extra:?}");
    }
    // Case 3 with binary `-out -` prefix before `pass 1\n`: prefix comes first, then messages, then counts.
    let (s_bin, out_bin, _, _) = run_web(&[
        "-mk_counts",
        "-in",
        fasta,
        "-unit",
        "5",
        "-sformat",
        "binary",
        "-out",
        "-",
    ]);
    assert_eq!(s_bin, 0);
    let (s_both, out_both, err_both, files_both) = run_web(&[
        "-mk_counts",
        "-in",
        fasta,
        "-unit",
        "5",
        "-sformat",
        "binary",
        "-out",
        "-",
        "-logfile",
        "/dev/stdout",
    ]);
    assert_eq!(s_both, 0);
    assert!(err_both.is_empty());
    assert!(files_both.is_empty());
    let prefix = crate::utils::winmask::CountsFormat::Binary.stream_prefix();
    assert!(!prefix.is_empty());
    assert!(out_both.starts_with(prefix));
    assert!(out_both[prefix.len()..].starts_with(b"pass 1\npass 2\n"));
    let mut expected_both = prefix.to_vec();
    expected_both.extend_from_slice(b"pass 1\npass 2\n");
    expected_both.extend_from_slice(&out_bin[prefix.len()..]);
    assert_eq!(out_both, expected_both);
    // Case 4: normalized file names for both `-logfile` and `-out`.
    let (status, _, err, mut files) = run_web(&[
        "-mk_counts",
        "-in",
        fasta,
        "-unit",
        "5",
        "-out",
        "./x.counts",
        "-logfile",
        " run.log",
    ]);
    assert_eq!(status, 0);
    assert!(err.is_empty());
    files.sort_by(|a, b| a.0.cmp(&b.0));
    assert_eq!(
        files
            .iter()
            .map(|(name, _)| name.as_slice())
            .collect::<Vec<_>>(),
        vec![b"run.log".as_slice(), b"x.counts".as_slice()]
    );
    let (status, _, err, files) = run_web(&[
        "-mk_counts",
        "-in",
        fasta,
        "-unit",
        "5",
        "-out",
        "x.counts",
        "-logfile",
        "d19-new-log/",
    ]);
    assert_eq!(status, 0);
    assert!(err.is_empty());
    assert!(files
        .iter()
        .any(|(name, bytes)| name == b"d19-new-log" && bytes == b"pass 1\npass 2\n"));
}
