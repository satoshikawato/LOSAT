use super::{BlastnArgs, QueryStrand};

/// BLASTN's options of the command line `losat blastn -query q.fa -subject s.fa <words>`.
fn blastn(words: &[&str]) -> BlastnArgs {
    let mut argv = vec!["losat", "blastn", "-query", "q.fa", "-subject", "s.fa"];
    argv.extend_from_slice(words);
    crate::cli::try_parse_from(argv).unwrap().blastn()
}

// blastn_args.cpp:57-59: -task takes NCBI's nucleotide tasks (NCBI's reader checks the set);
// describe lists them in the order of NCBI's BLASTN page, and `-help` names every one.
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
        assert_eq!(blastn(&["-task", task]).task, *task);
    }
    let Err(crate::cli::ArgvError::Printed(help)) =
        crate::cli::try_parse_from(["losat", "blastn", "-help"])
    else {
        panic!("blastn -help prints no help");
    };
    for task in &tasks {
        assert!(help.contains(&format!("'{task}'")), "{help}");
    }
    // NCBI BLAST+ 2.17.0 blastn `-task <name>` of a name that is not NCBI's (frozen on the
    // workstation with `env -i`, line RG item K4i-2b round 2): the constraint's error, NCBI's posts in LOSAT's words; the case of the name counts.
    for other in ["vecscreen", "MEGABLAST", "blastn-fast"] {
        let argv = [
            "losat", "blastn", "-query", "q.fa", "-subject", "s.fa", "-task", other,
        ];
        let Err(crate::cli::ArgvError::Report(bytes)) = crate::cli::try_parse_from(argv) else {
            panic!("{other}: not NCBI's argument error");
        };
        let text = String::from_utf8(bytes).unwrap();
        let line = format!(
            "Argument \"task\". Illegal value, expected Permissible values: 'blastn' 'blastn-short' 'dc-megablast' 'megablast' 'rmblastn' :  `{other}'"
        );
        let ncbi = format!("Error: {line}\nError:  (CArgException::eConstraint) {line}\n");
        assert_eq!(
            text.split_once(&format!("{}\n\n", "=".repeat(72)))
                .unwrap()
                .1,
            crate::wording::losat_text(&ncbi),
            "{other}"
        );
    }
}

// `-soft_masking` takes NCBI's boolean spellings (`NStr::StringToBool`); its default is true,
// also for blastn-short (the handler comes after `ClearFilterOptions`).
#[test]
fn soft_masking_sets_mask_at_hash() {
    let soft = |extra: &[&str]| blastn(extra).soft_masking();
    assert!(soft(&[]));
    assert!(soft(&["-task", "blastn-short"]));
    for value in ["true", "T", "1", "YES"] {
        assert!(soft(&["-soft_masking", value]), "{value}");
    }
    for value in ["false", "F", "0", "no"] {
        assert!(!soft(&["-soft_masking", value]), "{value}");
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
