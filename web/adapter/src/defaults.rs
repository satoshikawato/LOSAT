//! The defaults of the *describe* response (docs/web/abi_v2.md §9): the value that the engine
//! resolves for each option a search leaves out, by task, by the options that change other
//! options' values, and which options BLAST+ tells given from not given. Every value comes
//! from the engine's own resolution of an argv (`resolved`), so the application holds no
//! table and no rule of its own.

use std::collections::{BTreeMap, BTreeSet};

use clap::CommandFactory;
use LOSAT::web_adapter::{Cli, Commands};

use crate::json;
use crate::run::Program;

/// The options' values of one argv, by flag.
type Values = BTreeMap<&'static str, String>;

/// The hit list options: BLAST+ tells each given from not given (abi_v2.md §7), and
/// -max_target_seqs excludes the other two, so a default holds only while the others are not
/// given. Giving -num_descriptions at its default changes the tabular formats' hit list
/// size, which `resolved` (format 0) does not show, so they are named here.
const HIT_LIST: [(&str, &[&str]); 3] = [
    (
        "-max_target_seqs",
        &["-num_descriptions", "-num_alignments"],
    ),
    ("-num_descriptions", &["-max_target_seqs"]),
    ("-num_alignments", &["-max_target_seqs"]),
];

/// The word sizes tried for the values that a given word size changes (`by`): the
/// engine's own resolution decides which it accepts and what changes.
const WORD_SIZE_PROBES: std::ops::RangeInclusive<u32> = 1..=32;

/// The value of each option of `words` (program name first, as `run::parse` reads it) as the
/// engine resolves it before a search, by flag: the options without a value of their own
/// take their defaults, by task and from the other options.
pub fn resolved(words: &[&str]) -> Result<Values, String> {
    let (_, command) = crate::run::parse(words)?;
    let mut values = Values::new();
    let mut put = |flag: &'static str, value: String| {
        values.insert(flag, value);
    };
    let bool_text = |value: bool| value.to_string();
    match command {
        Commands::Blastn(mut args) => {
            use LOSAT::web_adapter::{
                configure_task, determine_effective_word_size, determine_evalue,
                determine_scoring_params, determine_window_size,
            };
            args.resolve_dust().map_err(|error| format!("{error:#}"))?;
            let (reward, penalty, gap_open, gap_extend) = determine_scoring_params(&args);
            let task = configure_task(&args);
            put("-evalue", determine_evalue(&args).to_string());
            put(
                "-word_size",
                determine_effective_word_size(&args).to_string(),
            );
            put("-reward", reward.to_string());
            put("-penalty", penalty.to_string());
            put("-gapopen", gap_open.to_string());
            put("-gapextend", gap_extend.to_string());
            put(
                "-dust",
                args.dust.params().map_or_else(
                    || "no".to_string(),
                    |(level, window, linker)| format!("{level} {window} {linker}"),
                ),
            );
            put("-soft_masking", bool_text(args.soft_masking()));
            put("-window_size", determine_window_size(&args).to_string());
            put("-xdrop_gap", task.x_drop_gapped.to_string());
            put("-xdrop_gap_final", task.x_drop_final.to_string());
            put("-culling_limit", args.culling_limit_value().to_string());
            hit_list(
                &mut put,
                args.max_target_seqs,
                args.num_descriptions,
                args.num_alignments,
            );
        }
        Commands::Blastp(args) => {
            let resolved = args.resolve().map_err(|error| format!("{error:#}"))?;
            put("-evalue", resolved.evalue.to_string());
            put("-word_size", resolved.word_size.to_string());
            put("-matrix", resolved.matrix_name.clone());
            put("-gapopen", resolved.scoring.gap_open.to_string());
            put("-gapextend", resolved.scoring.gap_extend.to_string());
            put("-threshold", resolved.threshold.to_string());
            put("-window_size", resolved.window_size.to_string());
            put("-seg", resolved.seg.to_ncbi_cli_string());
            put("-soft_masking", bool_text(resolved.soft_masking));
            put("-culling_limit", resolved.culling_limit.to_string());
            put("-max_target_seqs", resolved.max_target_seqs.to_string());
            put("-num_descriptions", resolved.num_descriptions.to_string());
            put("-num_alignments", resolved.num_alignments.to_string());
        }
        Commands::Tblastn(args) => {
            let resolved = args.resolve().map_err(|error| format!("{error:#}"))?;
            put("-evalue", resolved.evalue.to_string());
            put("-word_size", resolved.word_size.to_string());
            put("-matrix", resolved.matrix_name.clone());
            put("-gapopen", resolved.gap_open.to_string());
            put("-gapextend", resolved.gap_extend.to_string());
            put("-threshold", resolved.threshold.to_string());
            put("-window_size", resolved.window_size.to_string());
            put("-xdrop_gap", resolved.xdrop_gap.to_string());
            put("-xdrop_gap_final", resolved.xdrop_gap_final.to_string());
            put("-sum_stats", bool_text(resolved.sum_stats));
            put("-seg", resolved.seg.to_ncbi_cli_string());
            put("-soft_masking", bool_text(resolved.soft_masking));
            put("-culling_limit", resolved.culling_limit.to_string());
            put("-max_target_seqs", resolved.hitlist_size.to_string());
            put("-num_descriptions", resolved.num_descriptions.to_string());
            put("-num_alignments", resolved.num_alignments.to_string());
        }
        Commands::Tblastx(args) => {
            put("-evalue", args.evalue.to_string());
            put("-word_size", args.word_size().to_string());
            put("-matrix", args.matrix_name().to_string());
            put("-threshold", args.word_threshold().to_string());
            put("-window_size", args.window().to_string());
            put(
                "-seg",
                args.seg_spec()
                    .map_err(|error| format!("{error:#}"))?
                    .to_ncbi_cli_string(),
            );
            put("-soft_masking", bool_text(args.soft_masking()));
            put("-culling_limit", args.culling_limit().to_string());
            hit_list(
                &mut put,
                args.max_target_seqs,
                args.num_descriptions,
                args.num_alignments,
            );
        }
        Commands::Blastx(args) => {
            let resolved = args.resolve().map_err(|error| format!("{error:#}"))?;
            put("-evalue", resolved.evalue.to_string());
            put("-word_size", resolved.word_size.to_string());
            put("-matrix", resolved.matrix.clone());
            put("-gapopen", resolved.gap_open.to_string());
            put("-gapextend", resolved.gap_extend.to_string());
            put("-threshold", resolved.threshold.to_string());
            put("-window_size", resolved.window_size.to_string());
            put("-sum_stats", bool_text(resolved.sum_stats));
            put("-soft_masking", bool_text(resolved.soft_masking));
            put("-culling_limit", resolved.culling_limit.to_string());
            put("-max_target_seqs", resolved.hitlist_size.to_string());
            put("-num_descriptions", resolved.num_descriptions.to_string());
            put("-num_alignments", resolved.num_alignments.to_string());
        }
        // `run::parse` refuses the other commands (`Program::parse`).
        _ => return Err("no defaults for this command".to_string()),
    }
    Ok(values)
}

/// The hit list options of format 0 (`report_counts`), for the programs whose resolution
/// does not report them.
fn hit_list(
    put: &mut impl FnMut(&'static str, String),
    max_target_seqs: Option<usize>,
    num_descriptions: Option<usize>,
    num_alignments: Option<usize>,
) {
    let counts =
        LOSAT::web_adapter::report_counts(0, max_target_seqs, num_descriptions, num_alignments);
    put("-max_target_seqs", counts.hitlist_size.to_string());
    put("-num_descriptions", counts.num_descriptions.to_string());
    put("-num_alignments", counts.num_alignments.to_string());
}

/// The defaults of one program, as `describe` reports them.
pub struct Defaults {
    /// The values of each task with no other option given ("" for a program without -task).
    pub tasks: BTreeMap<String, Values>,
    /// For an option that changes other options' values (-matrix, -word_size): for each value
    /// that the engine accepts, the changed options' values by task. Only the options whose
    /// value differs from the task's for some value are listed.
    pub by: BTreeMap<&'static str, BTreeMap<String, BTreeMap<String, Values>>>,
    /// The options whose default holds only while these other options are not given.
    pub unless: BTreeMap<&'static str, Vec<&'static str>>,
    /// The options that BLAST+ tells given from not given, with the tasks where it does: given at
    /// the value shown as the default, they change the search.
    pub given: BTreeMap<&'static str, BTreeSet<String>>,
}

/// The argv of a program and task with `extra` options, for `resolved`.
fn words<'a>(program: Program, task: &'a str, extra: &[&'a str]) -> Vec<&'a str> {
    let mut words = vec![program.name(), "-query", "q", "-subject", "s"];
    if !task.is_empty() {
        words.extend(["-task", task]);
    }
    words.extend_from_slice(extra);
    words
}

/// The tasks of a program that the engine resolves, in the order of its parser ("" for a
/// program without -task).
fn tasks(program: Program) -> Vec<String> {
    let root = Cli::command();
    let Some(task) = root.find_subcommand(program.name()).and_then(|command| {
        command
            .get_arguments()
            .find(|arg| arg.get_long() == Some("task"))
    }) else {
        return vec![String::new()];
    };
    let mut names: Vec<String> = task
        .get_possible_values()
        .iter()
        .map(|value| value.get_name().to_string())
        .collect();
    names.retain(|name| resolved(&words(program, name, &[])).is_ok());
    names
}

/// Marks each option of `context` (the values that `extra`, no option or one -matrix or
/// -word_size value, resolves for `task`) given for the task where giving it at its value there
/// changes the resolution; `only` limits the options tried.
fn mark_given(
    given: &mut BTreeMap<&'static str, BTreeSet<String>>,
    program: Program,
    task: &str,
    extra: &[&str],
    context: &Values,
    only: Option<&BTreeSet<&'static str>>,
) {
    for (flag, value) in context {
        if extra.contains(flag)
            || only.is_some_and(|only| !only.contains(flag))
            || given.get(flag).is_some_and(|tasks| tasks.contains(task))
        {
            continue;
        }
        let mut words = words(program, task, extra);
        words.extend([*flag, value.as_str()]);
        if resolved(&words).is_ok_and(|values| &values != context) {
            given.entry(*flag).or_default().insert(task.to_string());
        }
    }
}

pub fn defaults(program: Program) -> Result<Defaults, String> {
    let mut out = Defaults {
        tasks: BTreeMap::new(),
        by: BTreeMap::new(),
        unless: BTreeMap::new(),
        given: BTreeMap::new(),
    };
    for task in tasks(program) {
        let values = resolved(&words(program, &task, &[]))?;
        out.tasks.insert(task, values);
    }
    let described = |flag: &str| {
        Cli::command()
            .find_subcommand(program.name())
            .is_some_and(|command| {
                command
                    .get_arguments()
                    .any(|arg| arg.get_long() == flag.strip_prefix('-'))
            })
    };
    let matrices: Vec<String> = LOSAT::web_adapter::ScoringMatrix::ALL
        .iter()
        .map(|matrix| matrix.ncbi_name().to_string())
        .collect();
    let word_sizes: Vec<String> = WORD_SIZE_PROBES.map(|size| size.to_string()).collect();
    for (decider, probes) in [("-matrix", matrices), ("-word_size", word_sizes)] {
        if !described(decider) {
            continue;
        }
        // The values of every accepted probe, then the options that any of them changes.
        let mut tried: BTreeMap<String, BTreeMap<String, Values>> = BTreeMap::new();
        let mut changed: BTreeSet<&'static str> = BTreeSet::new();
        for (task, base) in &out.tasks {
            for probe in &probes {
                let Ok(values) = resolved(&words(program, task, &[decider, probe.as_str()])) else {
                    continue;
                };
                for (flag, value) in &values {
                    if *flag != decider && base.get(flag) != Some(value) {
                        changed.insert(*flag);
                    }
                }
                tried
                    .entry(probe.clone())
                    .or_default()
                    .insert(task.clone(), values);
            }
        }
        if changed.is_empty() {
            continue;
        }
        // The values that the form shows where the decider is given: an option given at one
        // of them must give the same search, or it is told given for the task.
        for (probe, by_task) in &tried {
            for (task, values) in by_task {
                mark_given(
                    &mut out.given,
                    program,
                    task,
                    &[decider, probe],
                    values,
                    Some(&changed),
                );
            }
        }
        let table = tried
            .into_iter()
            .map(|(probe, by_task)| {
                let by_task = by_task
                    .into_iter()
                    .map(|(task, values)| {
                        let kept = values
                            .into_iter()
                            .filter(|(flag, _)| changed.contains(flag))
                            .collect();
                        (task, kept)
                    })
                    .collect();
                (probe, by_task)
            })
            .collect();
        out.by.insert(decider, table);
    }
    for (flag, others) in HIT_LIST {
        out.unless.insert(flag, others.to_vec());
        out.given.insert(flag, out.tasks.keys().cloned().collect());
    }
    // An option given at its task's default that changes the resolution (the word size of
    // the fast tasks sets NCBI's threshold for it).
    for (task, base) in &out.tasks {
        mark_given(&mut out.given, program, task, &[], base, None);
    }
    Ok(out)
}

/// Appends `defaults` as the *describe* JSON's `"defaults"` value.
pub fn write(out: &mut String, defaults: &Defaults) {
    let values = |out: &mut String, values: &Values| {
        out.push('{');
        for (index, (flag, value)) in values.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            json::string(out, flag);
            out.push(':');
            json::string(out, value);
        }
        out.push('}');
    };
    let by_task = |out: &mut String, by_task: &BTreeMap<String, Values>| {
        out.push('{');
        for (index, (task, task_values)) in by_task.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            json::string(out, task);
            out.push(':');
            values(out, task_values);
        }
        out.push('}');
    };
    out.push_str("{\"tasks\":");
    by_task(out, &defaults.tasks);
    out.push_str(",\"by\":{");
    for (index, (decider, table)) in defaults.by.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        json::string(out, decider);
        out.push_str(":{");
        for (probe_index, (probe, tasks)) in table.iter().enumerate() {
            if probe_index > 0 {
                out.push(',');
            }
            json::string(out, probe);
            out.push(':');
            by_task(out, tasks);
        }
        out.push('}');
    }
    out.push_str("},\"unless\":{");
    for (index, (flag, others)) in defaults.unless.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        json::string(out, flag);
        out.push_str(":[");
        for (other_index, other) in others.iter().enumerate() {
            if other_index > 0 {
                out.push(',');
            }
            json::string(out, other);
        }
        out.push(']');
    }
    out.push_str("},\"given\":{");
    for (index, (flag, tasks)) in defaults.given.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        json::string(out, flag);
        out.push_str(":[");
        for (task_index, task) in tasks.iter().enumerate() {
            if task_index > 0 {
                out.push(',');
            }
            json::string(out, task);
        }
        out.push(']');
    }
    out.push_str("}}");
}

#[cfg(test)]
mod tests {
    use super::*;

    const PROGRAMS: [Program; 5] = [
        Program::Blastn,
        Program::Blastp,
        Program::Tblastn,
        Program::Tblastx,
        Program::Blastx,
    ];

    /// The help text and the clap default of a program's option.
    fn help(program: Program, flag: &str) -> (String, Option<String>) {
        let root = Cli::command();
        let arg = root
            .find_subcommand(program.name())
            .and_then(|command| {
                command
                    .get_arguments()
                    .find(|arg| arg.get_long() == flag.strip_prefix('-'))
            })
            .unwrap_or_else(|| panic!("{} has no {flag}", program.name()));
        let defaults: Vec<String> = arg
            .get_default_values()
            .iter()
            .map(|value| value.to_string_lossy().into_owned())
            .collect();
        (
            arg.get_help().map(ToString::to_string).unwrap_or_default(),
            (!defaults.is_empty()).then(|| defaults.join(" ")),
        )
    }

    /// Whether `text` holds `value` as a word of its own ("3", not the 3 of "BLOSUM62" or "13").
    fn names(text: &str, value: &str) -> bool {
        text.match_indices(value).any(|(at, _)| {
            let before = text[..at].chars().next_back();
            let after = text[at + value.len()..].chars().next();
            !before.is_some_and(|c| c.is_alphanumeric() || matches!(c, '_' | '.' | '-'))
                && !after.is_some_and(|c| c.is_alphanumeric() || matches!(c, '_' | '.'))
        })
    }

    // The CLI states each default that the engine resolves for a task, as the clap default or
    // in the help: the help and the resolution agree (before, programs.ts held a copy of the
    // help, checked the same way).
    #[test]
    fn every_task_default_is_a_value_that_the_help_names() {
        let mut unnamed = Vec::new();
        for program in PROGRAMS {
            let defaults = defaults(program).unwrap();
            for (task, values) in &defaults.tasks {
                for (flag, value) in values {
                    let (text, clap_default) = help(program, flag);
                    let stated = match &clap_default {
                        Some(clap_default) => clap_default == value,
                        None => text.to_lowercase().contains("default") && names(&text, value),
                    };
                    if !stated {
                        unnamed.push(format!("{} {task} {flag} {value}: {text}", program.name()));
                    }
                }
            }
        }
        assert!(unnamed.is_empty(), "{unnamed:#?}");
    }

    #[test]
    fn tasks_are_the_engines() {
        let names =
            |program| -> Vec<String> { defaults(program).unwrap().tasks.into_keys().collect() };
        assert_eq!(
            names(Program::Blastn),
            ["blastn", "blastn-short", "dc-megablast", "megablast"]
        );
        assert_eq!(
            names(Program::Blastp),
            ["blastp", "blastp-fast", "blastp-short"]
        );
        assert_eq!(names(Program::Tblastn), ["tblastn", "tblastn-fast"]);
        assert_eq!(names(Program::Blastx), ["blastx", "blastx-fast"]);
        // TBLASTX has no -task.
        assert_eq!(names(Program::Tblastx), [""]);
    }

    #[test]
    fn blastn_resolves_each_tasks_preset() {
        let blastn = defaults(Program::Blastn).unwrap();
        let row = |flag: &str| -> Vec<&str> {
            ["megablast", "dc-megablast", "blastn", "blastn-short"]
                .iter()
                .map(|task| blastn.tasks[*task][flag].as_str())
                .collect()
        };
        assert_eq!(row("-word_size"), ["28", "11", "11", "7"]);
        assert_eq!(row("-reward"), ["1", "2", "2", "1"]);
        assert_eq!(row("-penalty"), ["-2", "-3", "-3", "-3"]);
        assert_eq!(row("-gapopen"), ["0", "5", "5", "5"]);
        assert_eq!(row("-gapextend"), ["0", "2", "2", "2"]);
        assert_eq!(row("-evalue"), ["10", "10", "10", "1000"]);
        assert_eq!(row("-dust"), ["20 64 1", "20 64 1", "20 64 1", "no"]);
        assert_eq!(row("-xdrop_gap"), ["25", "30", "30", "30"]);
        assert_eq!(row("-window_size"), ["0", "40", "0", "0"]);
    }

    #[test]
    fn a_matrix_decides_the_gap_costs_threshold_and_window() {
        // The engine's own values (blastp/args.rs tests): BLOSUM45 gives 14/2, threshold 14
        // and window 60; blastp-fast keeps its threshold 20.
        let blastp = defaults(Program::Blastp).unwrap();
        let blosum45 = &blastp.by["-matrix"]["BLOSUM45"];
        assert_eq!(blosum45["blastp"]["-gapopen"], "14");
        assert_eq!(blosum45["blastp"]["-gapextend"], "2");
        assert_eq!(blosum45["blastp"]["-threshold"], "14");
        assert_eq!(blosum45["blastp"]["-window_size"], "60");
        assert_eq!(blosum45["blastp-fast"]["-threshold"], "20");
        // Only the options that a matrix changes are listed.
        assert!(!blosum45["blastp"].contains_key("-evalue"));
        // TBLASTX: the threshold and the window, no gap costs.
        let tblastx = defaults(Program::Tblastx).unwrap();
        let pam30 = &tblastx.by["-matrix"]["PAM30"][""];
        assert_eq!(
            pam30.keys().copied().collect::<Vec<_>>(),
            ["-threshold", "-window_size"]
        );
        // BLASTX reads BLOSUM62 only, which changes nothing; BLASTN has no -matrix.
        assert!(!defaults(Program::Blastx)
            .unwrap()
            .by
            .contains_key("-matrix"));
        assert!(defaults(Program::Blastn).unwrap().by.is_empty());
    }

    #[test]
    fn a_word_size_given_sets_its_own_threshold() {
        for (program, fast) in [
            (Program::Blastp, "blastp-fast"),
            (Program::Tblastn, "tblastn-fast"),
            (Program::Blastx, "blastx-fast"),
        ] {
            let defaults = defaults(program).unwrap();
            // NCBI's 19.3 for a given word size 5, also at the fast task's own 5 (20 when
            // the word size is not given): the word size is told given from not given.
            assert_eq!(
                defaults.tasks[fast]["-threshold"],
                "20",
                "{}",
                program.name()
            );
            assert_eq!(
                defaults.by["-word_size"]["5"][fast]["-threshold"],
                "19.3",
                "{}",
                program.name()
            );
            // Told given only where it changes the search: the plain task's word size 3
            // leaves its threshold.
            assert_eq!(
                defaults.given["-word_size"],
                BTreeSet::from([fast.to_string()]),
                "{}",
                program.name()
            );
        }
        // BLASTX's parser reads word sizes 3 and 5 only.
        let blastx = defaults(Program::Blastx).unwrap();
        assert_eq!(
            blastx.by["-word_size"].keys().collect::<Vec<_>>(),
            ["3", "5"]
        );
        // TBLASTX's threshold does not follow the word size.
        assert!(!defaults(Program::Tblastx)
            .unwrap()
            .by
            .contains_key("-word_size"));
    }

    #[test]
    fn the_hit_list_options_are_given_and_exclude_each_other() {
        for program in PROGRAMS {
            let defaults = defaults(program).unwrap();
            let tasks: BTreeSet<String> = defaults.tasks.keys().cloned().collect();
            for (flag, others) in HIT_LIST {
                assert_eq!(defaults.given[flag], tasks, "{} {flag}", program.name());
                assert_eq!(defaults.unless[flag], others, "{} {flag}", program.name());
            }
            // No other option is told given, except the fast tasks' word size.
            let given: Vec<&str> = defaults.given.keys().copied().collect();
            let expected: &[&str] = match program {
                Program::Blastn | Program::Tblastx => {
                    &["-max_target_seqs", "-num_alignments", "-num_descriptions"]
                }
                _ => &[
                    "-max_target_seqs",
                    "-num_alignments",
                    "-num_descriptions",
                    "-word_size",
                ],
            };
            assert_eq!(given, expected, "{}", program.name());
            for values in defaults.tasks.values() {
                assert_eq!(values["-max_target_seqs"], "500");
                assert_eq!(values["-num_descriptions"], "500");
                assert_eq!(values["-num_alignments"], "250");
            }
        }
    }

    #[test]
    fn the_json_reads_as_the_tables() {
        let mut out = String::new();
        write(&mut out, &defaults(Program::Tblastx).unwrap());
        assert!(
            out.starts_with("{\"tasks\":{\"\":{\"-culling_limit\":\"0\","),
            "{out}"
        );
        assert!(out.contains(",\"by\":{\"-matrix\":{\"BLOSUM45\":{\"\":{\"-threshold\":\"16\",\"-window_size\":\"60\"}}"), "{out}");
        assert!(
            out.ends_with(
                ",\"given\":{\"-max_target_seqs\":[\"\"],\"-num_alignments\":[\"\"],\"-num_descriptions\":[\"\"]}}"
            ),
            "{out}"
        );
    }
}
