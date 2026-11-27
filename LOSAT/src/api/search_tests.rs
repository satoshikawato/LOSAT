use super::*;
use crate::api::Role;

#[test]
fn the_argv_is_parsed_as_on_the_command_line() {
    let search = validate(Program::Blastp, &["-evalue", "0.1", "-num_threads", "2"]).unwrap();
    assert_eq!(search.program(), Program::Blastp);
    assert_eq!(search.argv(), ["-evalue", "0.1", "-num_threads", "2"]);
    let inputs = Inputs::new(&[], &[]).named("q.faa", "s.faa");
    match search.command(inputs, &DiagTarget::Stderr).unwrap() {
        Commands::Blastp(args) => {
            assert_eq!(args.evalue, Some(0.1));
            assert_eq!(args.num_threads, 2);
            assert_eq!(args.query, std::path::PathBuf::from("q.faa"));
            assert_eq!(args.subject.as_deref(), Some(std::path::Path::new("s.faa")));
        }
        other => panic!("{other:?}"),
    }
    // An empty argv is the command line's defaults.
    validate(Program::Tblastx, &[]).unwrap();
    // A parser error is the command line's message.
    let error = validate(Program::Blastp, &["-evalue"]).unwrap_err();
    let cli = crate::cli::try_parse_from([
        "losat", "blastp", "-query", "query", "-subject", "subject", "-evalue",
    ])
    .unwrap_err();
    assert_eq!(error, Error::of_argv(cli));
    let error = validate(Program::Blastn, &["-nosuchoption", "1"]).unwrap_err();
    assert!(matches!(error, Error::Argv(_)), "{error}");
}

#[test]
fn the_options_the_api_owns_are_refused() {
    for word in [
        "-query",
        "--query",
        "-subject",
        "-db",
        "-out",
        "-outfmt",
        "-outfmt=6",
    ] {
        let error = validate(Program::Blastp, &[word, "x"]).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("not accepted by the library API"),
            "{word}: {error}"
        );
        assert!(error.to_string().starts_with(key(word)), "{word}: {error}");
    }
    for program in [Program::Blastp, Program::Blastx, Program::Psiblast] {
        assert!(validate(program, &["-parse_deflines"])
            .unwrap()
            .parse_deflines());
    }
    // A value that looks like an option is a value.
    validate(Program::Blastp, &["-matrix", "BLOSUM62"]).unwrap();
}

#[test]
fn the_search_strategy_file_is_refused_by_validate() {
    for program in [
        Program::Blastn,
        Program::Blastp,
        Program::Blastx,
        Program::Tblastn,
        Program::Tblastx,
        Program::Psiblast,
    ] {
        for argv in [
            &["-export_search_strategy", "s.asn"][..],
            &["-export_search_strategy=s.asn"],
            &["-evalue", "1", "--export_search_strategy", "s.asn"],
        ] {
            let key = argv.iter().find(|word| word.contains("export")).unwrap();
            assert_eq!(
                validate(program, argv),
                Err(Error::Argv(format!(
                    "{} is not accepted by the library API: the library API writes no search strategy file",
                    super::key(key)
                ).into())),
                "{program} {argv:?}"
            );
        }
    }
}

// Line LIB item 1.6: psiblast's alignment query is an input (`Inputs::msa`), which the search
// gives the parser as `-in_msa` in place of `-query`; in an argv the option is the API's, and
// only for psiblast (another program's parser rejects it as an unknown option).
#[test]
fn the_alignment_option_is_refused_by_validate() {
    for argv in [
        &["-in_msa", "m.aln"][..],
        &["-in_msa=m.aln"],
        &["-evalue", "1", "--in_msa", "m.aln"],
        &["-evalue", "1", "-in_msa", "m.aln"],
    ] {
        let word = argv
            .iter()
            .find(|word| word.contains("msa"))
            .copied()
            .unwrap();
        assert_eq!(
            validate(Program::Psiblast, argv),
            Err(Error::Argv(format!(
                "{} is not accepted by the library API: a query alignment is given as the inputs (Inputs::msa)",
                key(word)
            ).into())),
            "{argv:?}"
        );
    }
    for program in [Program::Blastp, Program::Tblastn] {
        let error = validate(program, &["-in_msa", "m.aln"]).unwrap_err();
        assert!(
            !error.to_string_lossy().contains("library API"),
            "{program}: {error}"
        );
    }
}

// Line LIB item 1.6: `-msa_master_idx` and `-ignore_msa_master` are options of the argv that
// `validate` takes, although the parser takes them with `-in_msa` only; the search gives the
// parser the alignment as the query (`Inputs::msa`), and any other inputs fail with the
// parser's own message, as the command line without `-in_msa`.
#[test]
fn the_alignment_master_options_go_with_the_alignment_input() {
    let path = std::path::Path::new("m.aln");
    let records = crate::api::read_fasta(Program::Psiblast, Role::Query, b">q\nMKT\n").unwrap();
    for argv in [
        &["-msa_master_idx", "2"][..],
        &["-msa_master_idx=2", "-evalue", "1"],
        &["-msa_master_idx", "3", "-num_iterations", "2"],
        &["-ignore_msa_master"],
    ] {
        let search = validate(Program::Psiblast, argv).unwrap();
        // The reader's options come from the argv before the inputs are known.
        assert!(!search.parse_deflines());
        for inputs in [
            Inputs::msa(path, &records),
            Inputs::msa_database(path, "db"),
        ] {
            assert_eq!(inputs.query_count(), 1);
            match search
                .command(inputs.named("q.faa", "s.faa"), &DiagTarget::Stderr)
                .unwrap()
            {
                Commands::Psiblast(args) => {
                    assert_eq!(args.in_msa.as_deref(), Some(path), "{argv:?}");
                    assert_eq!(args.msa_master_idx.is_some(), argv[0].contains("idx"));
                    assert_eq!(args.ignore_msa_master, argv[0].contains("ignore"));
                }
                other => panic!("{other:?}"),
            }
        }
        // Without the alignment the parser's message: the options need `-in_msa`.
        let error = search
            .command(Inputs::new(&records, &records), &DiagTarget::Stderr)
            .unwrap_err();
        assert!(matches!(error, Error::Argv(_)), "{argv:?}: {error:?}");
        assert!(!error.to_string().contains("library API"), "{error}");
        let error = search
            .command(Inputs::pssm(path, &records), &DiagTarget::Stderr)
            .unwrap_err();
        assert!(matches!(error, Error::Argv(_)), "{argv:?}: {error:?}");
    }
    // A value the parser refuses is its message, with the alignment as the query.
    let error = validate(Program::Psiblast, &["-msa_master_idx", "0"]).unwrap_err();
    assert!(matches!(error, Error::Argv(_)), "{error:?}");
    assert!(!error.to_string().contains("library API"), "{error}");
    // The two options exclude each other; `-query_loc` excludes the alignment.
    assert!(validate(
        Program::Psiblast,
        &["-msa_master_idx", "2", "-ignore_msa_master"]
    )
    .is_err());
    let ranged = validate(Program::Psiblast, &["-query_loc", "1-10"]).unwrap();
    let error = ranged
        .command(Inputs::msa(path, &[]), &DiagTarget::Stderr)
        .unwrap_err();
    assert!(
        matches!(&error, Error::Argv(_))
            && error.to_string().contains("Incompatible with argument:  `"),
        "{error}"
    );
    // The other programs have no such options: the parser's message.
    for program in [Program::Blastn, Program::Blastp, Program::Tblastn] {
        let error = validate(program, &["-msa_master_idx", "2"]).unwrap_err();
        assert!(
            !error.to_string().contains("library API"),
            "{program}: {error}"
        );
    }
}

// Line LIB item 1.6: an alignment query is psiblast's (a PSSM query is psiblast's and
// tblastn's), and a file, not the standard input.
#[test]
fn an_alignment_query_is_psiblasts_and_a_file() {
    let path = std::path::Path::new("m.aln");
    let search = validate(Program::Psiblast, &["-evalue", "1"]).unwrap();
    for inputs in [Inputs::msa(path, &[]), Inputs::msa_database(path, "db")] {
        let in_msa = match search
            .command(inputs.named("q.faa", "s.faa"), &DiagTarget::Stderr)
            .unwrap()
        {
            Commands::Psiblast(args) => args.in_msa,
            other => panic!("{other:?}"),
        };
        assert_eq!(in_msa.as_deref(), Some(path));
    }
    let standard_input = std::path::Path::new("-");
    let error = search
        .command(Inputs::msa(standard_input, &[]), &DiagTarget::Stderr)
        .unwrap_err();
    assert!(matches!(error, Error::Argv(_)), "{error:?}");
    assert!(error.to_string().contains("standard input"), "{error}");
    for program in [
        Program::Blastn,
        Program::Blastp,
        Program::Blastx,
        Program::Tblastn,
        Program::Tblastx,
    ] {
        let search = validate(program, &[]).unwrap();
        assert_eq!(
            search
                .command(Inputs::msa(path, &[]), &DiagTarget::Stderr)
                .unwrap_err(),
            Error::Argv(
                format!(
                    "a query alignment (Inputs::msa) is searched by psiblast, not by {}",
                    program.name()
                )
                .into()
            ),
            "{program}"
        );
    }
}

// Series Q item Q6a: psiblast's PSSM files go to the parser as on the command line, and the
// search writes them (`tests/library_api.rs`).
#[test]
fn psiblast_pssm_files_are_the_command_lines_options() {
    for argv in [
        &["-out_ascii_pssm", "p.txt"][..],
        &["-out_ascii_pssm=p.txt"],
        &["-out_pssm", "p.asn"],
        &[
            "-evalue",
            "1",
            "-save_pssm_after_last_round",
            "-out_pssm",
            "p.asn",
        ],
        &["-save_each_pssm", "-out_ascii_pssm", "p.txt"],
    ] {
        let search = validate(Program::Psiblast, argv).unwrap();
        let query = crate::api::read_fasta(Program::Psiblast, Role::Query, b">q\nMKT\n").unwrap();
        match search
            .command(Inputs::new(&query, &query), &DiagTarget::Stderr)
            .unwrap()
        {
            Commands::Psiblast(args) => assert!(
                args.out_pssm.is_some() || args.out_ascii_pssm.is_some(),
                "{argv:?}"
            ),
            other => panic!("{other:?}"),
        }
    }
    // NCBI's check of the options comes with the search, as on the command line
    // (`PsiblastPssmOptions::check`).
    validate(Program::Psiblast, &["-save_each_pssm"]).unwrap();
    // The parser takes one dash only.
    let error = validate(Program::Psiblast, &["--save_each_pssm"]).unwrap_err();
    assert!(!error.to_string().contains("library API"), "{error}");
    // blastp has no such option: the parser's message.
    let error = validate(Program::Blastp, &["-save_each_pssm"]).unwrap_err();
    assert!(!error.to_string().contains("PSSM"), "{error}");
}

// Series Q item Q6b: a PSSM query is an input (`Inputs::pssm`), which the search gives the
// parser as `-in_pssm` in place of `-query`; in an argv the option is the API's.
#[test]
fn a_pssm_query_is_an_input_that_the_search_gives_the_parser() {
    let refused = |key: &str| {
        Error::Argv(format!(
            "{key} is not accepted by the library API: a PSSM query is given as the inputs (Inputs::pssm)"
        ).into())
    };
    for program in [Program::Psiblast, Program::Tblastn] {
        assert_eq!(
            validate(program, &["-in_pssm", "p.asn"]),
            Err(refused("-in_pssm"))
        );
        assert_eq!(
            validate(program, &["--in_pssm=p.asn"]),
            Err(refused("--in_pssm"))
        );
    }
    // The other programs have no such option: the parser's message.
    for program in [
        Program::Blastn,
        Program::Blastp,
        Program::Blastx,
        Program::Tblastx,
    ] {
        let error = validate(program, &["-in_pssm", "p.asn"]).unwrap_err();
        assert!(
            !error.to_string().contains("library API"),
            "{program}: {error}"
        );
    }
    let path = std::path::Path::new("p.asn");
    for program in [Program::Psiblast, Program::Tblastn] {
        let search = validate(program, &["-evalue", "1"]).unwrap();
        for inputs in [Inputs::pssm(path, &[]), Inputs::pssm_database(path, "db")] {
            assert_eq!(inputs.query_count(), 1);
            let in_pssm = match search
                .command(inputs.named("q.faa", "s.faa"), &DiagTarget::Stderr)
                .unwrap()
            {
                Commands::Psiblast(args) => args.in_pssm,
                Commands::Tblastn(args) => args.in_pssm,
                other => panic!("{other:?}"),
            };
            assert_eq!(in_pssm.as_deref(), Some(path), "{program}");
        }
        // `-in_pssm` excludes `-query_loc`: the parser's error of the search, as `-db` and
        // `-subject_loc`.
        let ranged = validate(program, &["-query_loc", "1-10"]).unwrap();
        let error = ranged
            .command(Inputs::pssm(path, &[]), &DiagTarget::Stderr)
            .unwrap_err();
        assert!(
            matches!(&error, Error::Argv(_))
                && error.to_string().contains("Incompatible with argument:  `"),
            "{program}: {error}"
        );
    }
    let search = validate(Program::Psiblast, &[]).unwrap();
    let standard_input = std::path::Path::new("-");
    let error = search
        .command(Inputs::pssm(standard_input, &[]), &DiagTarget::Stderr)
        .unwrap_err();
    assert!(matches!(error, Error::Argv(_)), "{error:?}");
    assert!(error.to_string().contains("standard input"), "{error}");
    let search = validate(Program::Blastp, &[]).unwrap();
    assert_eq!(
        search
            .command(Inputs::pssm(path, &[]), &DiagTarget::Stderr)
            .unwrap_err(),
        Error::Argv(
            "a PSSM query (Inputs::pssm) is searched by psiblast and tblastn, not by blastp".into()
        )
    );
}

// Line RG item Q5: psiblast's rounds (`-num_iterations` other than 1) run in the library too;
// the records carry their round (`HspRecord::round`, `Outcome::rounds`).
#[test]
fn psiblast_rounds_are_accepted_by_validate_and_every_search() {
    let query = crate::api::read_fasta(
        Program::Psiblast,
        Role::Query,
        // Test protein: N-terminal part of E. coli AsnA (RefSeq NP_418200.1), Cys51 changed to Ala.
        b">q\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ\n",
    )
    .unwrap();
    let inputs = Inputs::new(&query, &query).named("q.faa", "s.faa");
    for (rounds, number) in [("0", 0), ("1", 1), ("3", 3)] {
        let search = validate(Program::Psiblast, &["-num_iterations", rounds]).unwrap();
        match search.command(inputs, &DiagTarget::Stderr).unwrap() {
            Commands::Psiblast(args) => assert_eq!(args.num_iterations, number),
            other => panic!("{other:?}"),
        }
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:79-87
// (cited at `ncbi_app::constraint::Constraint::MaxFileNameLength`)
/// `run_to_file`'s file is the parser's `-out` (`Search::writing_to`): any name, one that looks
/// like a toolkit word (`-dryrun`, `--`) or is not UTF-8 among them, is the file and leaves the
/// argv as it is; a name of 256 bytes or more is NCBI's error.
#[test]
fn the_file_of_run_to_file_is_the_parsers_out() {
    let search = validate(Program::Blastn, &["-ungapped"]).unwrap();
    let inputs = Inputs::new(&[], &[]).named("q.fna", "s.fna");
    let names = [
        std::ffi::OsString::from("-dryrun"),
        std::ffi::OsString::from("--"),
        std::ffi::OsString::from("-version"),
        std::ffi::OsString::from("dir/o.out"),
    ]
    .into_iter()
    // A name that is not UTF-8, which only Unix can make.
    .chain({
        #[cfg(unix)]
        let not_utf8 = Some(
            <std::ffi::OsStr as std::os::unix::ffi::OsStrExt>::from_bytes(b"o\xff.out")
                .to_os_string(),
        );
        #[cfg(not(unix))]
        let not_utf8 = None;
        not_utf8
    });
    for name in names {
        let writing = search.writing_to(std::path::Path::new(&name));
        assert!(!writing.dry_run(), "{name:?}");
        match writing.command(inputs, &DiagTarget::Stderr).unwrap() {
            Commands::Blastn(args) => {
                assert_eq!(args.out.as_deref(), Some(std::path::Path::new(&name)));
                assert!(args.ungapped, "{name:?}");
            }
            other => panic!("{other:?}"),
        }
    }
    let long = "a".repeat(256);
    let error = search
        .writing_to(std::path::Path::new(&long))
        .command(inputs, &DiagTarget::Stderr)
        .unwrap_err();
    assert!(matches!(error, Error::Argv(_)), "{error:?}");
    assert!(
        error
            .to_string()
            .contains("Illegal value, expected file name length < 256"),
        "{error}"
    );
}
