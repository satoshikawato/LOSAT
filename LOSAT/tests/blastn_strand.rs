//! `blastn -strand plus|minus|both` on whole searches: the option is accepted (an unknown value is
//! a clap error), `both` is the search without the option, and the searches of the two strands
//! together are the search of both (NCBI BLAST+ 2.17.0 `-num_threads 1`, checked over
//! 543 + 91 probe cases: the lines of `plus` are the plus/plus lines of `both` and the lines of
//! `minus` the plus/minus lines, with the same e-values; only the order of the subjects can
//! change). That the bytes equal NCBI's is the business of the probe (`ncbi-parity` O2a) and the
//! fixtures.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta/outfmt0")
        .join(name)
}

fn run(query: &str, subject: &str, extra: &[&str]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    command
        .arg("blastn")
        .args(["-query"])
        .arg(fasta(query))
        .args(["-subject"])
        .arg(fasta(subject))
        .args(["-outfmt", "6"])
        .args(extra);
    for (name, _) in std::env::vars() {
        if name.starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command.output().expect("LOSAT starts")
}

fn lines(output: &Output) -> Vec<String> {
    assert!(
        output.status.success() && output.stderr.is_empty(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout.clone())
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect()
}

/// The subject strand of an outfmt 6 line: `true` for the plus strand (sstart < send).
fn subject_plus(line: &str) -> bool {
    let fields: Vec<&str> = line.split('\t').collect();
    fields[8].parse::<u64>().unwrap() < fields[9].parse::<u64>().unwrap()
}

const CASES: [(&str, &str); 2] = [
    ("strand_query.fasta", "strand_subject.fasta"),
    ("multi_query.fasta", "multi_subject.fasta"),
];

#[test]
fn an_unknown_strand_is_a_clap_error() {
    for value in ["foo", "PLUS", "Both", "", "1"] {
        let output = run(CASES[0].0, CASES[0].1, &["-strand", value]);
        assert_eq!(output.status.code(), Some(2), "-strand {value:?}");
        assert!(output.stdout.is_empty(), "-strand {value:?}");
        assert!(!output.stderr.is_empty(), "-strand {value:?}");
    }
}

#[test]
fn both_is_the_search_without_the_option_and_the_strands_add_up() {
    for (query, subject) in CASES {
        for task in ["megablast", "blastn", "blastn-short"] {
            let base = ["-task", task];
            let none = run(query, subject, &base);
            let both = run(query, subject, &[&base[..], &["-strand", "both"]].concat());
            assert_eq!(none.stdout, both.stdout, "{query} {task}");
            let none = lines(&none);
            for threads in ["1", "4"] {
                let one_strand = |strand: &str| {
                    lines(&run(
                        query,
                        subject,
                        &[&base[..], &["-strand", strand, "-num_threads", threads]].concat(),
                    ))
                };
                let (plus, minus) = (one_strand("plus"), one_strand("minus"));
                assert!(plus.iter().all(|line| subject_plus(line)), "{query} {task}");
                assert!(
                    minus.iter().all(|line| !subject_plus(line)),
                    "{query} {task}"
                );
                let mut together: Vec<String> = plus.into_iter().chain(minus).collect();
                let mut want = none.clone();
                together.sort();
                want.sort();
                assert_eq!(together, want, "{query} {task} threads {threads}");
            }
        }
    }
}

#[test]
fn the_equals_form_is_accepted() {
    let (query, subject) = CASES[0];
    let spaced = run(query, subject, &["-strand", "minus"]);
    let joined = run(query, subject, &["-strand=minus"]);
    assert_eq!(lines(&spaced), lines(&joined));
    assert!(!lines(&joined).is_empty());
}
