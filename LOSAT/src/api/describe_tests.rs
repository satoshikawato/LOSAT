use super::*;

#[test]
fn api_owned_options_are_not_described_and_search_options_are() {
    for program in Program::ALL {
        let description = describe(program);
        assert_eq!(description.program, program);
        let flags: Vec<&str> = description
            .parameters
            .iter()
            .map(|parameter| parameter.flag.as_str())
            .collect();
        for owned in API_OWNED
            .iter()
            .chain(file_options(program))
            .chain(alignment_owned(program))
        {
            assert!(
                !flags.contains(&format!("-{owned}").as_str()),
                "{program} lists -{owned}"
            );
        }
        for expected in ["-evalue", "-num_threads", "-max_target_seqs"] {
            assert!(flags.contains(&expected), "{program} lacks {expected}");
        }
        let evalue = description
            .parameters
            .iter()
            .find(|parameter| parameter.flag == "-evalue")
            .unwrap();
        assert!(evalue.takes_value);
        // The help texts are NCBI's comments of the options (line RG DECISIONS K4i-3-1).
        assert!(description
            .parameters
            .iter()
            .any(|parameter| !parameter.help.is_empty()));
    }
}

#[test]
fn task_choices_are_the_tasks_the_parser_accepts() {
    let parses = |program: Program, task: &str| {
        let argv = [
            "losat",
            program.name(),
            "-query",
            "q",
            "-subject",
            "s",
            "-task",
            task,
        ];
        crate::cli::try_parse_from(argv).is_ok()
    };
    // Every -task lists its tasks: the default among them, each one accepted.
    for program in Program::ALL {
        let description = describe(program);
        let Some(task) = description
            .parameters
            .iter()
            .find(|parameter| parameter.flag == "-task")
        else {
            continue;
        };
        let default = task.default.as_deref().unwrap();
        assert!(
            task.choices.iter().any(|choice| choice == default),
            "{program}"
        );
        for choice in &task.choices {
            assert!(parses(program, choice), "{program} -task {choice}");
        }
    }
    let blastn = describe(Program::Blastn);
    let task = blastn
        .parameters
        .iter()
        .find(|parameter| parameter.flag == "-task")
        .unwrap();
    assert_eq!(
        task.choices,
        [
            "megablast",
            "dc-megablast",
            "blastn",
            "blastn-short",
            "rmblastn"
        ]
    );
    // NCBI's rmblastn (blastn with reward 0 and penalty 0) is offered and parses.
    assert!(parses(Program::Blastn, "rmblastn"));
}

#[test]
fn genetic_codes_are_listed_for_the_translated_programs_only() {
    let tblastx = describe(Program::Tblastx);
    assert!(tblastx.query_gencodes.contains(&1));
    assert!(tblastx.query_gencodes.contains(&11));
    assert!(!tblastx.query_gencodes.contains(&0));
    assert!(!tblastx.subject_gencodes.is_empty());
    let blastx = describe(Program::Blastx);
    assert!(!blastx.query_gencodes.is_empty());
    assert!(blastx.subject_gencodes.is_empty());
    let tblastn = describe(Program::Tblastn);
    assert!(tblastn.query_gencodes.is_empty());
    assert!(!tblastn.subject_gencodes.is_empty());
    for program in [Program::Blastn, Program::Blastp] {
        let description = describe(program);
        assert!(description.query_gencodes.is_empty(), "{program}");
        assert!(description.subject_gencodes.is_empty(), "{program}");
    }
}

// NCBI's reader decides which genetic codes a run takes (line RG item K4i-2b): describe's codes,
// which come from the same argument tables, are those that the command line takes (the probe of
// every code below 64 that describe ran before K4i-2a).
#[test]
fn the_genetic_codes_are_those_that_the_parser_takes() {
    for program in Program::ALL {
        let description = describe(program);
        for (name, listed) in [
            ("-query_gencode", &description.query_gencodes),
            ("-db_gencode", &description.subject_gencodes),
        ] {
            let taken: Vec<u32> = (0..64)
                .filter(|code| {
                    let code = code.to_string();
                    let argv = [
                        "losat",
                        program.name(),
                        "-query",
                        "q",
                        "-subject",
                        "s",
                        name,
                        &code,
                    ];
                    crate::cli::try_parse_from(argv).is_ok()
                })
                .collect();
            assert_eq!(listed, &taken, "{program} {name}");
        }
    }
}

// Line RG item K4i-3a: describe reads NCBI's argument tables, and the options it lists are those
// that the command line's options read (`from_args`, line RG item K4i-3c): the names of each
// program's XML but the standard arguments and the arguments refused by name. The value kinds are
// the tables' types; psiblast's `-out_pssm` and `-out_ascii_pssm`, NCBI's output files, are paths.
#[test]
fn the_described_options_are_those_that_the_command_line_reads() {
    for program in Program::ALL {
        let tool = crate::cli::reader::tool(program.name())
            .expect("every API program is a program of the CLI");
        let mut read: Vec<String> = crate::cli::reader::options_tests::names_read(tool)
            .into_iter()
            .map(|name| format!("-{name}"))
            .collect();
        let mut described: Vec<String> = parameters(program)
            .map(|parameter| parameter.flag)
            .collect();
        read.sort();
        described.sort();
        assert_eq!(described, read, "{program}");
    }
    let kind = |flag: &str| {
        parameters(Program::Psiblast)
            .find(|parameter| parameter.flag == flag)
            .map(|parameter| parameter.value_kind)
    };
    assert_eq!(kind("-out_pssm"), Some(ValueKind::Path));
    assert_eq!(kind("-out_ascii_pssm"), Some(ValueKind::Path));
}

// The help of an option is NCBI's comment of it (the `<description>` of its XML, which `-help`
// prints), in LOSAT's words (`wording::losat_text`), and the options come in the order of the usage
// line of `-help`.
#[test]
fn the_help_is_ncbis_comment_in_the_order_of_the_help() {
    let blastn = describe(Program::Blastn);
    let help = |flag: &str| {
        let parameter = blastn
            .parameters
            .iter()
            .find(|parameter| parameter.flag == flag);
        parameter.map(|parameter| parameter.help.as_str())
    };
    assert_eq!(
        help("-evalue"),
        Some(
            "Expectation value (E) threshold for saving hits. Default = 10 (1000 for blastn-short)"
        )
    );
    assert_eq!(
        help("-num_alignments"),
        Some("Number of database sequences to show alignments for\nDefault = `250'")
    );
    let flags: Vec<&str> = blastn
        .parameters
        .iter()
        .map(|parameter| parameter.flag.as_str())
        .collect();
    // NCBI's usage of blastn: `-task` first (after the search strategy files, which the API does
    // not take), `-evalue`, `-word_size` and `-gapopen` in turn, `-mt_mode` last (before `-remote`,
    // which LOSAT refuses).
    let at = |flag: &str| flags.iter().position(|listed| *listed == flag).unwrap();
    assert_eq!(flags.first(), Some(&"-task"));
    assert!(at("-evalue") < at("-word_size") && at("-word_size") < at("-gapopen"));
    assert_eq!(flags.last(), Some(&"-mt_mode"));
}
