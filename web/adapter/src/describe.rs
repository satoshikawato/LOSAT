//! The *describe* response (docs/web/abi_v2.md §9), generated from the engine's clap
//! definitions and value parsers, so that the application shows exactly the options,
//! defaults and allowed values of the command line.

use std::fmt::Write as _;
use std::sync::{Mutex, PoisonError};

use clap::CommandFactory;
use LOSAT::web_adapter::Cli;

use crate::json;
use crate::run::Program;

/// Options that the adapter or the host sets itself or refuses (`run::parse`): the output
/// file and format (a run writes every format of the program; a host passes the words of a
/// `-outfmt` of its own, whose other formats and 6 or 7 with fields add an output,
/// `run::parse_run`) and the thread count, and the
/// options of NCBI's command line that the application has no control for (`-mt_mode`,
/// which needs `-num_threads` and has no effect; `-use_index`, which a search of records ignores
/// as NCBI ignores it with `-subject`;
/// `-parse_deflines`, since a record is read before the options of a run are known).
/// `-num_descriptions` and `-num_alignments` are described: a run searches once for each
/// distinct hit list size of its formats (`run::search_groups`).
const ADAPTER_OWNED: [&str; 6] = [
    "out",
    "outfmt",
    "num_threads",
    "mt_mode",
    "use_index",
    "parse_deflines",
];

/// Genetic codes are small integers; every accepted one is below this bound.
const GENETIC_CODE_PROBE: u32 = 64;

/// The options of `program` that a search of records takes none of (the PSSM files: tblastn's
/// `-in_pssm`, which the browser refuses until line PS's PSW1): the library API's list
/// (`LOSAT::web_adapter::records_refused_options`), so that both `describe`s leave out the same
/// options.
fn records_refused(program: Program) -> &'static [&'static str] {
    let program = match program {
        Program::Blastp => LOSAT::api::Program::Blastp,
        Program::Tblastn => LOSAT::api::Program::Tblastn,
        Program::Blastn => LOSAT::api::Program::Blastn,
        Program::Tblastx => LOSAT::api::Program::Tblastx,
        Program::Blastx => LOSAT::api::Program::Blastx,
    };
    LOSAT::web_adapter::records_refused_options(program)
}

/// Each program's response after its first `describe`: the response is pure, and its
/// `defaults` resolve some hundred probe argvs (`defaults::defaults`), while the application
/// asks again each time it opens a result.
static DESCRIBED: Mutex<Vec<(Program, String)>> = Mutex::new(Vec::new());

pub fn describe(program: &str) -> Result<String, String> {
    let program = Program::parse(program)?;
    let mut described = DESCRIBED.lock().unwrap_or_else(PoisonError::into_inner);
    if let Some((_, response)) = described.iter().find(|(known, _)| *known == program) {
        return Ok(response.clone());
    }
    let response = describe_program(program)?;
    described.push((program, response.clone()));
    Ok(response)
}

fn describe_program(program: Program) -> Result<String, String> {
    let root = Cli::command();
    let command = root
        .find_subcommand(program.name())
        .ok_or_else(|| format!("the engine has no {} command", program.name()))?;
    let mut out = String::from("{\"program\":");
    json::string(&mut out, program.name());
    out.push_str(",\"formats\":[");
    for (index, format) in program.formats().iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(out, "{format}");
    }
    out.push_str("],\"parameters\":[");
    let mut first = true;
    let refused = records_refused(program);
    for arg in command.get_arguments() {
        let Some(long) = arg.get_long() else { continue };
        if matches!(long, "help" | "version")
            || ADAPTER_OWNED.contains(&long)
            || refused.contains(&long)
        {
            continue;
        }
        if !first {
            out.push(',');
        }
        first = false;
        out.push_str("{\"flag\":");
        json::string(&mut out, &format!("-{long}"));
        out.push_str(",\"help\":");
        json::string(
            &mut out,
            &arg.get_help().map(ToString::to_string).unwrap_or_default(),
        );
        let _ = write!(out, ",\"takes_value\":{}", arg.get_action().takes_values());
        let defaults: Vec<String> = arg
            .get_default_values()
            .iter()
            .map(|value| value.to_string_lossy().into_owned())
            .collect();
        if !defaults.is_empty() {
            out.push_str(",\"default\":");
            json::string(&mut out, &defaults.join(" "));
        }
        let choices: Vec<String> = arg
            .get_possible_values()
            .iter()
            .map(|value| value.get_name().to_string())
            .collect();
        if !choices.is_empty() {
            out.push_str(",\"choices\":[");
            for (index, choice) in choices.iter().enumerate() {
                if index > 0 {
                    out.push(',');
                }
                json::string(&mut out, choice);
            }
            out.push(']');
        }
        out.push('}');
    }
    out.push(']');
    for (long, key) in [
        ("query_gencode", "query_gencodes"),
        ("db_gencode", "subject_gencodes"),
    ] {
        if !command
            .get_arguments()
            .any(|arg| arg.get_long() == Some(long))
        {
            continue;
        }
        // The engine's own parser decides which codes it accepts.
        let flag = format!("-{long}");
        let accepted: Vec<u32> = (0..GENETIC_CODE_PROBE)
            .filter(|code| {
                let code = code.to_string();
                crate::run::parse(&[program.name(), "-query", "q", "-subject", "s", &flag, &code])
                    .is_ok()
            })
            .collect();
        let _ = write!(out, ",\"{key}\":[");
        for (index, code) in accepted.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            let _ = write!(out, "{code}");
        }
        out.push(']');
    }
    out.push_str(",\"defaults\":");
    crate::defaults::write(&mut out, &crate::defaults::defaults(program)?);
    out.push('}');
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn options_the_adapter_owns_or_refuses_are_not_described() {
        for program in ["blastp", "tblastn", "blastn", "tblastx", "blastx"] {
            let described = describe(program).unwrap();
            for flag in ADAPTER_OWNED {
                assert!(
                    !described.contains(&format!("\"-{flag}\"")),
                    "{program} lists -{flag}"
                );
            }
            // The hit list options are described with the CLI's help text, which carries
            // the defaults (they are `Option`s: "given or not" is NCBI's).
            for (flag, default) in [
                ("num_descriptions", "[default: 500]"),
                ("num_alignments", "[default: 250]"),
            ] {
                let entry = described
                    .split("{\"flag\":")
                    .find(|entry| entry.starts_with(&format!("\"-{flag}\"")))
                    .unwrap_or_else(|| panic!("{program} does not list -{flag}"));
                assert!(entry.contains(default), "{program} -{flag}: {entry}");
                assert!(entry.contains("\"takes_value\":true"), "{program} -{flag}");
            }
            assert!(described.contains("\"-max_target_seqs\""), "{program}");
        }
    }

    // tblastn's `-in_pssm` (a PSSM query file): the browser refuses it until line PS's PSW1, and
    // neither `describe` lists it (the library API's list, `records_refused`).
    #[test]
    fn options_a_search_of_records_does_not_take_are_not_described() {
        assert_eq!(records_refused(Program::Tblastn), ["in_pssm"]);
        for program in ["blastp", "tblastn", "blastn", "tblastx", "blastx"] {
            let described = describe(program).unwrap();
            for flag in records_refused(Program::parse(program).unwrap()) {
                assert!(
                    !described.contains(&format!("\"-{flag}\"")),
                    "{program} lists -{flag}"
                );
            }
        }
    }
}
