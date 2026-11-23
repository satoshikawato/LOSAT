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
            .chain(pattern_owned(program))
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
        // The help texts are the command line's; some options have none.
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
        crate::cli::try_parse_from::<Cli, _, _>(argv).is_ok()
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
