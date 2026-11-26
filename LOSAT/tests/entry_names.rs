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
const PROGRAMS: [&str; 21] = [
    "blastn",
    "blastp",
    "blastx",
    "tblastn",
    "tblastx",
    "psiblast",
    "deltablast",
    "rpsblast",
    "rpstblastn",
    "makeblastdb",
    "blastdbcmd",
    "dustmasker",
    "segmasker",
    "windowmasker",
    "convert2blastmask",
    "makembindex",
    "makeprofiledb",
    "rpsbproc",
    "blastdb_aliastool",
    "blastdbcheck",
    "cleanup-blastdb-volumes.py",
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
        let first_line = if name == "cleanup-blastdb-volumes.py" {
            format!("{name} LOSAT {}\n", env!("CARGO_PKG_VERSION"))
        } else {
            format!("{name}: LOSAT {}\n", env!("CARGO_PKG_VERSION"))
        };
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
            "makeprofiledb",
            "blastdb_aliastool",
            "blastdbcheck",
            "cleanup-blastdb-volumes.py",
        ] {
            let line = text
                .lines()
                .find(|line| line.split_whitespace().next() == Some(program))
                .unwrap();
            assert!(line.trim().len() > program.len() + 2, "{line:?}");
        }
    }
}

/// The programs' argument errors are NCBI's reader's (its usage, the line of `=` and its two
/// posts, exit status 1; line RG item K4i-2b) and LOSAT's refusals by name are one line after
/// `Error: ` (status 1); both end with one line end, so that the shell prompt starts a line of its
/// own. The top level's error (an unknown program) stays LOSAT's own, status 2.
#[test]
fn argument_errors_and_refusals_end_with_a_line_end() {
    for (words, ncbi) in [
        (
            &["blastn", "-query", "q.fa", "-subject", "q.fa", "-remote"][..],
            false,
        ),
        (&["blastn", "-no_such_option"], true),
        (&["blastn", "-query"], true),
        (&["makeblastdb", "-no_such_option"], true),
        (&["rpsbproc", "-no_such_option"], true),
    ] {
        let output = run(losat(), words);
        assert_eq!(output.status.code(), Some(1), "{words:?} {output:?}");
        let text = String::from_utf8_lossy(&output.stderr).into_owned();
        if ncbi {
            assert!(text.starts_with("USAGE\n"), "{words:?}: {text:?}");
            assert!(
                text.contains(&format!("\n{}\n\nError: ", "=".repeat(72))),
                "{words:?}: {text:?}"
            );
        } else {
            assert!(
                text.starts_with("Error: the NCBI BLAST+ option -"),
                "{words:?}: {text:?}"
            );
            assert_eq!(text.lines().count(), 1, "{words:?}: {text:?}");
        }
        assert!(
            text.ends_with('\n') && !text.ends_with("\n\n"),
            "{words:?}: {text:?}"
        );
        for clap in ["error: ", "For more information"] {
            assert!(!text.contains(clap), "{words:?}: {clap:?} in {text:?}");
        }
    }
    let output = run(losat(), &["no_such_program"]);
    assert_eq!(output.status.code(), Some(2), "{output:?}");
    let text = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(text.starts_with("error: "), "{text:?}");
    assert!(text.ends_with('\n') && !text.ends_with("\n\n"), "{text:?}");
}

/// No `-help` or `-h` text shows a source reference (they belong in `//` comments). The top
/// level's help is LOSAT's, whose hint names the other form in NCBI's single-dash form. The
/// programs' help is NCBI BLAST+ 2.17.0's text in LOSAT's words (line RG item K4i-2b): `-h` the
/// short usage with NCBI's hint at its foot, `-help` the full text, `-help-full` the same but
/// makeblastdb's (which also lists the hidden `-randomize`); `--help` is no word of theirs
/// (NCBI's `Unknown argument`, status 1).
#[test]
fn help_has_no_source_references_and_a_true_footer() {
    let leaks = ["ncbi-blast/", "```", "Reference:", "NCBI reference"];
    let mut texts = Vec::new();
    for flag in ["-h", "-help", "--help"] {
        let output = run(losat(), &[flag]);
        assert!(output.status.success(), "{flag}: {output:?}");
        let text = String::from_utf8_lossy(&output.stdout).into_owned();
        for leak in leaks {
            assert!(!text.contains(leak), "{flag}: {leak:?}\n{text}");
        }
        assert!(!text.contains("'--help'"), "{flag}\n{text}");
        texts.push(text);
    }
    assert_eq!(texts[1], texts[2], "--help is -help");
    assert_eq!(
        texts[0].contains("(see more with '-help')"),
        texts[1].contains("(see a summary with '-h')")
    );
    for program in [
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
    ] {
        let mut texts = Vec::new();
        for flag in ["-h", "-help", "-help-full"] {
            let output = run(losat(), &[program, flag]);
            assert!(output.status.success(), "{program} {flag}: {output:?}");
            assert!(output.stderr.is_empty(), "{program} {flag}: {output:?}");
            let text = String::from_utf8_lossy(&output.stdout).into_owned();
            for leak in leaks {
                assert!(!text.contains(leak), "{program} {flag}: {leak:?}\n{text}");
            }
            assert!(text.starts_with("USAGE\n"), "{program} {flag}\n{text}");
            texts.push(text);
        }
        let (short, full, help_full) = (&texts[0], &texts[1], &texts[2]);
        assert!(
            short.ends_with(
                "\nUse '-help' to print detailed descriptions of command line arguments\n"
            ),
            "{program}: {short}"
        );
        assert!(full.len() > short.len(), "{program}");
        if program == "makeblastdb" {
            assert!(!full.contains(" -randomize\n"), "{full}");
            assert!(help_full.contains(" -randomize\n"), "{help_full}");
        } else {
            assert_eq!(full, help_full, "{program}: -help-full is -help");
        }
        let output = run(losat(), &[program, "--help"]);
        assert_eq!(output.status.code(), Some(1), "{program}: {output:?}");
        assert!(
            String::from_utf8_lossy(&output.stderr)
                .ends_with("\nError:  Unknown argument: \"-help\"\n"),
            "{program}: {output:?}"
        );
    }
}

/// One form for LOSAT's refusal of an NCBI option by name, whatever the program and the other
/// arguments are (blastn, blastp and blastx, with and without `-subject`: the same bytes but the
/// program's name, exit status 1), after NCBI's reader took the command line (line RG item K4i-2b).
#[test]
fn refusals_by_name_have_one_form() {
    for option in ["-remote"] {
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
                assert_eq!(output.status.code(), Some(1), "{words:?} {output:?}");
                assert!(output.stdout.is_empty(), "{words:?}");
                let text = String::from_utf8(output.stderr).unwrap();
                forms.push(text.replace(&program.to_uppercase(), "PROGRAM"));
            }
        }
        let expected =
            format!("Error: the NCBI BLAST+ option {option} is not supported by LOSAT's PROGRAM\n");
        assert!(forms.iter().all(|form| *form == expected), "{forms:?}");
    }
    // blastx's remaining refusal by name (`-entrez_query` with `-remote`) has the same form
    // (`-gilist`, `-ipglist`, and `-import_search_strategy` run now).
    {
        let option = "-entrez_query";
        let argv = [
            "blastx",
            "-query",
            "q.fa",
            "-db",
            "d",
            "-entrez_query",
            "x",
            "-remote",
        ];
        let output = run(losat(), &argv);
        assert_eq!(output.status.code(), Some(1), "{option} {output:?}");
        assert_eq!(
            String::from_utf8(output.stderr).unwrap(),
            format!("Error: the NCBI BLAST+ option {option} is not supported by LOSAT's BLASTX\n")
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
