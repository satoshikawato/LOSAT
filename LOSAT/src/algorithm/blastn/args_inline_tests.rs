use super::{BlastnArgs, QueryStrand};
use clap::{CommandFactory, Parser};

#[derive(Parser)]
struct TestCli {
    #[command(flatten)]
    args: BlastnArgs,
}

// blastn_args.cpp:57-59: -task takes NCBI's nucleotide tasks (NCBI's reader checks the set);
// describe lists them in the order of NCBI's BLASTN page, and the help texts name every one.
#[test]
fn task_takes_ncbis_tasks_and_the_help_names_them() {
    let tasks = crate::ncbi_app::args::choices("blastn", "task").unwrap();
    assert_eq!(
        tasks,
        [
            "megablast",
            "dc-megablast",
            "blastn",
            "blastn-short",
            "rmblastn"
        ]
    );
    for task in &tasks {
        let argv = [
            "losat",
            "--query",
            "q.fa",
            "--subject",
            "s.fa",
            "--task",
            task,
        ];
        assert_eq!(TestCli::try_parse_from(argv).unwrap().args.task, *task);
    }
    let command = TestCli::command();
    let arg = command
        .get_arguments()
        .find(|arg| arg.get_long() == Some("task"))
        .unwrap();
    for help in [arg.get_help(), arg.get_long_help()] {
        let help = help.unwrap().to_string();
        for task in &tasks {
            assert!(help.contains(task.as_str()), "{help}");
        }
    }
}

// `-soft_masking` takes NCBI's boolean spellings (`NStr::StringToBool`); its default is true,
// also for blastn-short (the handler comes after `ClearFilterOptions`).
#[test]
fn soft_masking_sets_mask_at_hash() {
    let soft = |extra: &[&str]| {
        let mut argv = vec!["losat", "--query", "q.fa", "--subject", "s.fa"];
        argv.extend_from_slice(extra);
        TestCli::parse_from(argv).args.soft_masking()
    };
    assert!(soft(&[]));
    assert!(soft(&["--task", "blastn-short"]));
    for value in ["true", "T", "1", "YES"] {
        assert!(soft(&["--soft_masking", value]), "{value}");
    }
    for value in ["false", "F", "0", "no"] {
        assert!(!soft(&["--soft_masking", value]), "{value}");
    }
}

#[test]
fn the_strand_option_is_read_once() {
    assert_eq!(QueryStrand::from_option("both").searched(), (true, true));
    assert_eq!(QueryStrand::from_option("plus").searched(), (true, false));
    assert_eq!(QueryStrand::from_option("minus").searched(), (false, true));
}

#[test]
fn an_unknown_strand_is_an_argument_error() {
    for value in ["foo", "", "Plus", "PLUS", "both "] {
        let error = QueryStrand::parse(value).unwrap_err().to_string();
        assert_eq!(
            error,
            format!(
                "Argument \"strand\". Illegal value, expected `both', `minus', `plus':  `{value}'"
            )
        );
    }
}
