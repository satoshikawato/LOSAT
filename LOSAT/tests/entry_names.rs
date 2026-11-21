//! `losat --install-shims` (`cli::install_shims`) and the executables it links, called by a
//! program's name (`cli::entry_argv`; its rule is checked in `src/utils/entry_names_tests.rs`).
//! No NCBI counterpart: NCBI has one executable per program, and these
//! names stand for them on `PATH`.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn losat() -> &'static str {
    env!("CARGO_BIN_EXE_losat")
}

/// A new directory of this test.
fn temp_dir(tag: &str) -> PathBuf {
    // In the target directory: the hard link of `the_default_directory_...` needs the binary's
    // file system.
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("losat-entry-{tag}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// One run with no `LOSAT_*` variable that chooses a code path, in LOSAT/.
fn run(program: impl AsRef<std::ffi::OsStr>, words: &[&str]) -> Output {
    let mut command = Command::new(program);
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .args(words)
        .output()
        .unwrap()
}

fn shim(dir: &Path, name: &str) -> PathBuf {
    dir.join(format!("{name}{}", std::env::consts::EXE_SUFFIX))
}

fn install(dir: &Path) -> Output {
    run(losat(), &["--install-shims", dir.to_str().unwrap()])
}

/// The programs that the binary links, in the order it lists them.
const PROGRAMS: [&str; 16] = [
    "blastn",
    "blastp",
    "blastx",
    "tblastn",
    "tblastx",
    "psiblast",
    "rpsblast",
    "rpstblastn",
    "makeblastdb",
    "blastdbcmd",
    "dustmasker",
    "segmasker",
    "windowmasker",
    "convert2blastmask",
    "makembindex",
    "rpsbproc",
];

/// The program names that `losat --install-shims` listed on standard output (one path per
/// line): the file names without the executable suffix. The list must hold `PROGRAMS`.
fn installed_names(installed: &Output) -> Vec<String> {
    let names: Vec<String> = String::from_utf8(installed.stdout.clone())
        .unwrap()
        .lines()
        .map(|line| {
            let file = Path::new(line).file_name().unwrap().to_str().unwrap();
            file.strip_suffix(std::env::consts::EXE_SUFFIX)
                .unwrap_or(file)
                .to_string()
        })
        .collect();
    assert!(!names.is_empty(), "{installed:?}");
    for program in PROGRAMS {
        assert!(
            names.iter().any(|name| name == program),
            "{program}: {names:?}"
        );
    }
    names
}

#[test]
fn a_shim_prints_what_the_subcommand_prints() {
    let dir = temp_dir("texts");
    let installed = install(&dir);
    assert!(installed.status.success(), "{installed:?}");
    for name in installed_names(&installed).iter().map(String::as_str) {
        // `blastn -version` is what tools read the version from (`blastn: LOSAT 0.2.0`, where NCBI's
        // is `blastn: 2.17.0+`; rpsbproc and the maskers write LOSAT's too, where NCBI's are
        // `rpsbproc: 0.5.1`, `dustmasker: 1.0.0`, `segmasker: 1.0.0`, `windowmasker: 1.0.0`,
        // `makembindex: 2.17.0`).
        let version = run(shim(&dir, name), &["-version"]);
        assert!(version.status.success(), "{name}: {version:?}");
        let first_line = format!("{name}: LOSAT {}\n", env!("CARGO_PKG_VERSION"));
        assert!(
            version.stdout.starts_with(first_line.as_bytes()),
            "{name}: {version:?}"
        );
        for words in [
            &["-version"][..],
            &["-h"],
            &["-help"],
            &[],
            &["-no_such_option"],
        ] {
            let mut through_losat = vec![name];
            through_losat.extend(words);
            let expected = run(losat(), &through_losat);
            let output = run(shim(&dir, name), words);
            assert_eq!(
                output.status.code(),
                expected.status.code(),
                "{name} {words:?}"
            );
            assert_eq!(output.stdout, expected.stdout, "{name} {words:?}");
            assert_eq!(output.stderr, expected.stderr, "{name} {words:?}");
        }
    }
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_shim_searches_and_builds_as_the_subcommand_does() {
    let dir = temp_dir("runs");
    assert!(install(&dir).status.success());
    let db = dir.join("db");
    let db = db.to_str().unwrap();
    let fasta = "tests/db_search/taxfilter/q.fna";
    let built = run(
        shim(&dir, "makeblastdb"),
        &["-in", fasta, "-dbtype", "nucl", "-out", db],
    );
    assert!(built.status.success(), "{built:?}");
    let searches: [(&str, &[&str]); 3] = [
        ("blastn", &["-query", fasta, "-db", db, "-outfmt", "6"]),
        (
            "tblastx",
            &[
                "-query",
                fasta,
                "-subject",
                fasta,
                "-outfmt",
                "6",
                "-max_target_seqs",
                "1",
            ],
        ),
        ("blastdbcmd", &["-db", db, "-info"]),
    ];
    for (name, words) in searches {
        let mut through_losat = vec![name];
        through_losat.extend(words);
        let expected = run(losat(), &through_losat);
        let output = run(shim(&dir, name), words);
        assert!(expected.status.success(), "{name}: {expected:?}");
        assert!(!expected.stdout.is_empty(), "{name}");
        assert_eq!(output.stdout, expected.stdout, "{name}");
        assert_eq!(output.stderr, expected.stderr, "{name}");
    }
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn installing_twice_is_fine_and_another_file_is_left_alone() {
    let dir = temp_dir("taken");
    let first = install(&dir);
    assert!(first.status.success(), "{first:?}");
    let listed: Vec<PathBuf> = String::from_utf8(first.stdout.clone())
        .unwrap()
        .lines()
        .map(PathBuf::from)
        .collect();
    let dir = fs::canonicalize(&dir).unwrap();
    let names: Vec<PathBuf> = PROGRAMS.iter().map(|name| shim(&dir, name)).collect();
    assert_eq!(listed, names);
    let again = install(&dir);
    assert!(again.status.success(), "{again:?}");
    assert_eq!(again.stdout, first.stdout);

    // NCBI's program in the directory: not replaced, and no link is made beside it.
    let other = temp_dir("taken2");
    fs::write(shim(&other, "blastp"), "NCBI\n").unwrap();
    let refused = install(&other);
    assert_eq!(refused.status.code(), Some(1), "{refused:?}");
    assert!(
        String::from_utf8_lossy(&refused.stderr).contains("blastp"),
        "{refused:?}"
    );
    assert_eq!(fs::read(shim(&other, "blastp")).unwrap(), b"NCBI\n");
    assert!(!shim(&other, "blastn").exists());

    for words in [
        &["--install-shims", "a", "b"][..],
        &["--install-shims", "-h"],
    ] {
        let usage = run(losat(), words);
        assert_eq!(usage.status.code(), Some(1), "{words:?} {usage:?}");
    }
    assert!(!Path::new(env!("CARGO_MANIFEST_DIR")).join("-h").exists());

    // A link left by a `losat` that moved is replaced; the single-dash form of `-help` works.
    let moved = temp_dir("moved");
    #[cfg(unix)]
    std::os::unix::fs::symlink("/no/such/dir/losat", shim(&moved, "blastx")).unwrap();
    #[cfg(not(unix))]
    fs::copy(losat(), shim(&moved, "blastx")).unwrap();
    let replaced = run(losat(), &["-install-shims", moved.to_str().unwrap()]);
    assert!(replaced.status.success(), "{replaced:?}");
    let version = run(shim(&moved, "blastx"), &["-version"]);
    assert!(
        version
            .stdout
            .starts_with(concat!("blastx: LOSAT ", env!("CARGO_PKG_VERSION"), "\n").as_bytes()),
        "{version:?}"
    );
    let help = run(losat(), &["-h"]);
    assert!(
        String::from_utf8_lossy(&help.stdout).contains("losat -install-shims [DIR]"),
        "{help:?}"
    );
    fs::remove_dir_all(&moved).unwrap();
    fs::remove_dir_all(&dir).unwrap();
    fs::remove_dir_all(&other).unwrap();
}

#[test]
fn the_default_directory_is_the_executables_and_the_links_are_relative() {
    let dir = temp_dir("default");
    let copy = dir.join(format!("losat{}", std::env::consts::EXE_SUFFIX));
    // A hard link: a copy just written can be busy for the exec of another test's fork.
    fs::hard_link(losat(), &copy).unwrap();
    let installed = run(&copy, &["--install-shims"]);
    assert!(installed.status.success(), "{installed:?}");
    let version = run(shim(&dir, "blastn"), &["-version"]);
    assert!(
        version
            .stdout
            .starts_with(concat!("blastn: LOSAT ", env!("CARGO_PKG_VERSION"), "\n").as_bytes()),
        "{version:?}"
    );
    #[cfg(unix)]
    assert_eq!(
        fs::read_link(shim(&dir, "blastn")).unwrap(),
        Path::new("losat")
    );
    fs::remove_dir_all(&dir).unwrap();
}

/// The top level's `-help` names every program that `losat` runs (the `Commands:` list and the
/// "Other programs" list, derived from the entry names), and none of the notes that went stale.
#[test]
fn top_level_help_names_every_program() {
    for words in [&["-help"][..], &["-h"], &["--help"]] {
        let help = run(losat(), words);
        let text = String::from_utf8_lossy(&help.stdout).into_owned();
        for program in PROGRAMS {
            assert!(
                text.lines()
                    .any(|line| line.split_whitespace().next() == Some(program)),
                "{words:?}: no line for {program}\n{text}"
            );
        }
        for stale in [
            "format version 4\n",
            "its first round",
            "(local blast",
            "(local tblast",
        ] {
            assert!(!text.contains(stale), "{words:?}: {stale:?}\n{text}");
        }
        // NCBI's headings stay (BRIG reads them), and so does the shims note.
        assert!(text.starts_with("USAGE\n"), "{text}");
        assert!(text.contains("\nDESCRIPTION\n"), "{text}");
        assert!(text.contains("losat -install-shims [DIR]"), "{text}");
        // An "Other programs" line has a description.
        for program in [
            "blastdbcmd",
            "dustmasker",
            "segmasker",
            "windowmasker",
            "convert2blastmask",
            "makembindex",
        ] {
            let line = text
                .lines()
                .find(|line| line.split_whitespace().next() == Some(program))
                .unwrap();
            assert!(line.trim().len() > program.len() + 2, "{line:?}");
        }
    }
}

/// Every search program's `-help` says which `-outfmt` numbers it writes and which take a field
/// list (README's "Output formats" tables; `cli::outfmt_help`), and none of the old notes.
#[test]
fn search_program_help_says_what_outfmt_writes() {
    for (program, written, fields) in [
        (
            "blastn",
            "Written: 0 to 18 and 20.",
            "Formats 6, 7, 10 and 20 take",
        ),
        (
            "blastp",
            "Written: 0 to 16, 18 and 20.",
            "Formats 6, 7, 10 and 20 take",
        ),
        (
            "blastx",
            "Written: 0 to 16, 18 and 20.",
            "Formats 6, 7, 10 and 20 take",
        ),
        (
            "tblastn",
            "Written: 0 to 16, 18 and 20; with -in_pssm: 0 to 4, 6, 7 and 10.",
            "Formats 6, 7, 10 and 20 take",
        ),
        (
            "tblastx",
            "Written: 0 to 16, 18 and 20.",
            "Formats 6, 7, 10 and 20 take",
        ),
        (
            "psiblast",
            "Written: 0, 6, 7, 10 and 20; with -in_pssm or -in_msa: 0, 6 and 7.",
            "Formats 6, 7, 10 and 20 take",
        ),
        (
            "rpsblast",
            "Written: 0, 6, 7 and 10.",
            "Formats 6, 7 and 10 take",
        ),
    ] {
        for flag in ["-help", "-h"] {
            let help = run(losat(), &[program, flag]);
            let text = String::from_utf8_lossy(&help.stdout).into_owned();
            let text = text.split_whitespace().collect::<Vec<_>>().join(" ");
            // blastx keeps the value name `<OUTFMT>` of its frozen Web ABI v1 error text.
            let after = text
                .split_once("-outfmt <SPEC> ")
                .or_else(|| text.split_once("-outfmt <OUTFMT> "))
                .unwrap_or_else(|| panic!("{program} {flag}: no -outfmt line\n{text}"))
                .1;
            assert!(
                after.starts_with("Output format. "),
                "{program} {flag}: {after}"
            );
            assert!(after.contains(written), "{program} {flag}: {after}");
            assert!(after.contains(fields), "{program} {flag}: {after}");
            assert!(
                after.contains("a field list, for example -outfmt \"6 qseqid sseqid evalue\"")
                    && after.contains("qaccver saccver pident length mismatch"),
                "{program} {flag}: {after}"
            );
            for stale in [
                "not yet ported",
                "without custom fields",
                "Supported formats",
            ] {
                assert!(!text.contains(stale), "{program} {flag}: {stale:?}");
            }
        }
    }
}

/// LOSAT's own refusals and syntax errors (not NCBI's text) end with a line end, as clap's own
/// errors do, so that the shell prompt starts a line of its own (`cli::command_line_message`).
#[test]
fn argument_refusals_end_with_a_line_end() {
    for words in [
        &["blastn", "-query", "q.fa", "-subject", "q.fa", "-remote"][..],
        &["blastn", "-no_such_option"],
        &["blastn", "-query"],
        &["makeblastdb", "-no_such_option"],
        &["rpsbproc", "-no_such_option"],
        &["no_such_program"],
    ] {
        let output = run(losat(), words);
        assert_eq!(output.status.code(), Some(2), "{words:?} {output:?}");
        let text = String::from_utf8_lossy(&output.stderr).into_owned();
        assert!(text.starts_with("error: "), "{words:?}: {text:?}");
        assert!(
            text.ends_with('\n') && !text.ends_with("\n\n"),
            "{words:?}: {text:?}"
        );
    }
}

/// No program's `-help` or `-h` text shows a source reference (they belong in `//` comments),
/// and the hint at the foot of the help is true: `-h` prints the summary and `-help` the full
/// text, so a hint names the other one in NCBI's single-dash form.
#[test]
fn help_has_no_source_references_and_a_true_footer() {
    let mut programs = vec![""];
    programs.extend([
        "blastn",
        "blastp",
        "blastx",
        "tblastn",
        "tblastx",
        "psiblast",
        "rpsblast",
        "rpstblastn",
        "makeblastdb",
        "rpsbproc",
    ]);
    for program in programs {
        let words = |flag: &'static str| -> Vec<&'static str> {
            if program.is_empty() {
                vec![flag]
            } else {
                vec![program, flag]
            }
        };
        let mut texts = Vec::new();
        for flag in ["-h", "-help", "--help"] {
            let output = run(losat(), &words(flag));
            assert!(output.status.success(), "{program} {flag}: {output:?}");
            let text = String::from_utf8_lossy(&output.stdout).into_owned();
            for leak in ["ncbi-blast/", "```", "Reference:", "NCBI reference"] {
                assert!(!text.contains(leak), "{program} {flag}: {leak:?}\n{text}");
            }
            // No hint names `--help`, which the programs do not show.
            assert!(!text.contains("'--help'"), "{program} {flag}\n{text}");
            texts.push(text);
        }
        let (short, full) = (&texts[0], &texts[1]);
        assert_eq!(full, &texts[2], "{program}: --help is -help");
        // A hint at the foot of the short text promises a longer `-help`; one at the foot of
        // the full text promises a shorter `-h`; no hint, no difference.
        let short_hint = short.contains("(see more with '-help')");
        let full_hint = full.contains("(see a summary with '-h')");
        assert_eq!(short_hint, full_hint, "{program}");
        if short_hint {
            assert!(full.len() > short.len(), "{program}");
        } else if program != "rpstblastn" {
            assert_eq!(short, full, "{program}: -h and -help differ without a hint");
        }
    }
    // The one search program whose long help differs from its summary.
    let short = String::from_utf8_lossy(&run(losat(), &["blastn", "-h"]).stdout).into_owned();
    let full = String::from_utf8_lossy(&run(losat(), &["blastn", "-help"]).stdout).into_owned();
    assert!(
        short
            .trim_end()
            .ends_with("Print help (see more with '-help')"),
        "{short}"
    );
    assert!(
        full.trim_end()
            .ends_with("Print help (see a summary with '-h')"),
        "{full}"
    );
}

/// One form for LOSAT's refusal of an NCBI option by name, whatever the program and the other
/// arguments are (blastn, blastp and blastx, with and without `-subject`: the same bytes but the
/// program's name, exit status 2), and the hint of an argument error names `-help`, NCBI's form.
#[test]
fn refusals_by_name_have_one_form_and_the_hint_says_help() {
    let option = "-remote";
    let mut forms = Vec::new();
    for program in ["blastn", "blastp", "blastx"] {
        for with_subject in [false, true] {
            let mut words = vec![program, "-query", "q.fa"];
            if with_subject {
                words.extend(["-subject", "s.fa"]);
            } else {
                words.extend(["-db", "d"]);
            }
            words.push(option);
            let output = run(losat(), &words);
            assert_eq!(output.status.code(), Some(2), "{words:?} {output:?}");
            assert!(output.stdout.is_empty(), "{words:?}");
            let text = String::from_utf8(output.stderr).unwrap();
            forms.push(text.replace(&program.to_uppercase(), "PROGRAM"));
        }
    }
    let expected =
        format!("error: the NCBI BLAST+ option {option} is not supported by LOSAT's PROGRAM\n");
    assert!(forms.iter().all(|form| *form == expected), "{forms:?}");
    // blastx's other refusals by name (`-gilist`, `-use_sw_tback`) have the same form.
    for option in ["-gilist", "-use_sw_tback"] {
        let output = run(losat(), &["blastx", "-query", "q.fa", option]);
        assert_eq!(output.status.code(), Some(2), "{option} {output:?}");
        assert_eq!(
            String::from_utf8(output.stderr).unwrap(),
            format!("error: the NCBI BLAST+ option {option} is not supported by LOSAT's BLASTX\n")
        );
    }
}

/// blastx refuses a protein subject of more than 5,000,000 residues (NCBI searches it in
/// overlapping parts) in one form, from a `-subject` record and from a database sequence: the
/// same wording but the name of the subject, exit status 1, no output.
#[test]
fn blastx_refuses_a_subject_over_5_000_000_residues_in_one_form() {
    let dir = temp_dir("longsubject");
    let query = dir.join("q.fna");
    fs::write(&query, format!(">q\n{}\n", "ATGGCTGCTAAAGAAGGT".repeat(20))).unwrap();
    let subject = dir.join("s.faa");
    fs::write(&subject, format!(">long\n{}\n", "W".repeat(5_000_001))).unwrap();
    let wording = |what: &str| {
        format!(
            "Error: {what} has 5000001 residues; a subject longer than 5000000 residues, which NCBI BLAST+ searches in overlapping parts, is not supported by LOSAT's BLASTX\n"
        )
    };
    let q = query.to_str().unwrap();
    let s = subject.to_str().unwrap();
    let out = dir.join("out.txt");
    let o = out.to_str().unwrap();
    let by_subject = run(
        losat(),
        &[
            "blastx", "-query", q, "-subject", s, "-outfmt", "6", "-out", o,
        ],
    );
    assert_eq!(by_subject.status.code(), Some(1), "{by_subject:?}");
    assert!(by_subject.stdout.is_empty(), "{by_subject:?}");
    assert_eq!(
        String::from_utf8(by_subject.stderr).unwrap(),
        wording("subject sequence 1")
    );
    assert!(fs::read(&out).map_or(true, |bytes| bytes.is_empty()));

    let db = dir.join("db");
    let made = run(
        losat(),
        &[
            "makeblastdb",
            "-in",
            s,
            "-dbtype",
            "prot",
            "-out",
            db.to_str().unwrap(),
        ],
    );
    assert!(made.status.success(), "{made:?}");
    let by_db = run(
        losat(),
        &[
            "blastx",
            "-query",
            q,
            "-db",
            db.to_str().unwrap(),
            "-outfmt",
            "6",
            "-out",
            o,
        ],
    );
    assert_eq!(by_db.status.code(), Some(1), "{by_db:?}");
    assert!(by_db.stdout.is_empty(), "{by_db:?}");
    assert_eq!(
        String::from_utf8(by_db.stderr).unwrap(),
        wording("database sequence 0")
    );
    assert!(fs::read(&out).map_or(true, |bytes| bytes.is_empty()));
    fs::remove_dir_all(&dir).unwrap();
}

// The bioconda recipe's tests run `<program> -version` for every program the package links, and
// its description names them: a program added to `PROGRAMS` goes into the recipe too (makembindex
// had been left out of it). Not in a crate built from its package, which has no recipe.
#[test]
fn the_bioconda_recipe_tests_every_program() {
    let recipe =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../packaging/bioconda/meta.yaml.template");
    let Ok(recipe) = fs::read_to_string(&recipe) else {
        assert!(
            !Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../.git")
                .exists(),
            "{} is missing",
            recipe.display()
        );
        return;
    };
    let description = recipe.split("description:").nth(1).unwrap_or_default();
    for name in PROGRAMS {
        assert!(
            recipe.contains(&format!("    - {name} -version\n")),
            "the recipe's tests have no `{name} -version`"
        );
        assert!(
            description.contains(name),
            "the recipe's description does not name {name}"
        );
    }
}
