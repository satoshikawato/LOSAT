use super::{BlastnArgs, QueryStrand, BLASTN_TASKS};
use clap::{CommandFactory, Parser};

#[derive(Parser)]
struct TestCli {
    #[command(flatten)]
    args: BlastnArgs,
}

// blastn_args.cpp:57-59: -task takes NCBI's nucleotide tasks; LOSAT runs those of
// `BLASTN_TASKS` (rmblastn among them) and lists its tasks for any other name, as for the
// other programs' -task (`task_parser`).
#[test]
fn task_accepts_the_task_table_and_lists_it() {
    let parse = |task: &str| {
        TestCli::try_parse_from([
            "losat",
            "--query",
            "q.fa",
            "--subject",
            "s.fa",
            "--task",
            task,
        ])
    };
    for task in BLASTN_TASKS {
        assert_eq!(parse(task).unwrap().args.task, *task);
    }
    let command = TestCli::command();
    let arg = command
        .get_arguments()
        .find(|arg| arg.get_long() == Some("task"))
        .unwrap();
    let listed: Vec<String> = arg
        .get_possible_values()
        .iter()
        .map(|value| value.get_name().to_string())
        .collect();
    assert_eq!(listed, BLASTN_TASKS);
    // The help texts name every task in place of clap's list.
    assert!(arg.is_hide_possible_values_set());
    for help in [arg.get_help(), arg.get_long_help()] {
        let help = help.unwrap().to_string();
        for task in BLASTN_TASKS {
            assert!(help.contains(task), "{help}");
        }
    }
    // rmblastn: NCBI's matrix-only task (blast_options_handle.cpp:365-370).
    assert_eq!(parse("rmblastn").unwrap().args.task, "rmblastn");
    for other in ["vecscreen", "MEGABLAST", "blastn-fast"] {
        let error = parse(other).err().unwrap();
        assert_eq!(
            error.kind(),
            clap::error::ErrorKind::InvalidValue,
            "{other}"
        );
        assert!(
            error.to_string().contains(
                "[possible values: megablast, dc-megablast, blastn, blastn-short, rmblastn]"
            ),
            "{error}"
        );
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
