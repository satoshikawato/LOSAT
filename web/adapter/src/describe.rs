//! The *describe* response (docs/web/abi_v2.md §9), generated from the engine's clap
//! definitions and value parsers, so that the application shows exactly the options,
//! defaults and allowed values of the command line.

use std::fmt::Write as _;

use clap::CommandFactory;
use LOSAT::cli::Cli;

use crate::json;
use crate::run::Program;

/// Options that the adapter sets itself.
const ADAPTER_OWNED: [&str; 3] = ["out", "outfmt", "num_threads"];

/// Genetic codes are small integers; every accepted one is below this bound.
const GENETIC_CODE_PROBE: u32 = 64;

pub fn describe(program: &str) -> Result<String, String> {
    let program = Program::parse(program)?;
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
    for arg in command.get_arguments() {
        let Some(long) = arg.get_long() else { continue };
        if matches!(long, "help" | "version") || ADAPTER_OWNED.contains(&long) {
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
    out.push('}');
    Ok(out)
}
