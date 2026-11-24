//! The *describe* response (docs/web/abi_v2.md §9), generated from NCBI's argument tables that the
//! command line's reader reads (the options in the order of the usage line of its `-help`, their
//! help, defaults and allowed values: the library API's rows,
//! `losat::web_adapter::described_parameters`), so that the application shows exactly the options,
//! defaults and allowed values of the command line.

use std::fmt::Write as _;
use std::sync::{Mutex, PoisonError};

use crate::json;
use crate::run::Program;

/// Options that the adapter or the host sets itself or refuses (`run::parse`): the output
/// file and format (a run writes every format of the program; a host passes the words of a
/// `-outfmt` of its own, whose other formats and 6 or 7 with fields add an output,
/// `run::parse_run`) and the thread count, and the
/// options of NCBI's command line that the application has no control for (`-mt_mode`,
/// which needs `-num_threads` and has no effect; `-use_index` and `-index_name`, which a search
/// of records ignores as NCBI ignores them with `-subject`).
/// `-num_descriptions` and `-num_alignments` are described: a run searches once for each
/// distinct hit list size of its formats (`run::search_groups`).
const ADAPTER_OWNED: [&str; 6] = [
    "out",
    "outfmt",
    "num_threads",
    "mt_mode",
    "use_index",
    "index_name",
];

/// The options of `program` that the browser's searches of records take none of (tblastn's
/// `-in_pssm`, which the browser refuses until line PS's PSW1): the browser's list
/// (`losat::web_adapter::records_refused_options`), which the adapter's `describe` leaves out. The
/// search strategy file `-export_search_strategy` is a file of the run (line RG's F1w), described.
fn records_refused(program: Program) -> &'static [&'static str] {
    losat::web_adapter::records_refused_options(program.api())
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
    for parameter in losat::web_adapter::described_parameters(program.api()) {
        let long = parameter.flag.trim_start_matches('-');
        if ADAPTER_OWNED.contains(&long) || refused.contains(&long) {
            continue;
        }
        if !first {
            out.push(',');
        }
        first = false;
        out.push_str("{\"flag\":");
        json::string(&mut out, &parameter.flag);
        out.push_str(",\"help\":");
        json::string(&mut out, &parameter.help);
        let _ = write!(out, ",\"takes_value\":{}", parameter.takes_value);
        // The default and the values that NCBI's argument reader allows (line RG DECISIONS
        // K4i-2b-1): none for a flag.
        if let Some(default) = &parameter.default {
            out.push_str(",\"default\":");
            json::string(&mut out, default);
        }
        if !parameter.choices.is_empty() {
            out.push_str(",\"choices\":[");
            for (index, choice) in parameter.choices.iter().enumerate() {
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
        // The codes of the key's constraint in the engine's argument tables (NCBI's, and LOSAT's
        // 32 of tblastn's `-db_gencode`); no list for a program without the key.
        let Some(accepted) = losat::web_adapter::genetic_codes(program.name(), long) else {
            continue;
        };
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
            // The hit list options are described with NCBI's help text, which carries the
            // defaults (NCBI's table has none: "given or not" is NCBI's).
            for (flag, default) in [
                ("num_descriptions", "Default = `500'"),
                ("num_alignments", "Default = `250'"),
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
    // the adapter's `describe` does not list it (the browser's list, `records_refused`). The search
    // strategy file of the five programs is a file of the run (line RG's F1w): listed.
    #[test]
    fn options_a_search_of_records_does_not_take_are_not_described() {
        assert_eq!(records_refused(Program::Tblastn), ["in_pssm"]);
        for program in [
            Program::Blastn,
            Program::Blastp,
            Program::Blastx,
            Program::Tblastx,
        ] {
            assert!(records_refused(program).is_empty());
        }
        for program in ["blastp", "tblastn", "blastn", "tblastx", "blastx"] {
            let described = describe(program).unwrap();
            for flag in records_refused(Program::parse(program).unwrap()) {
                assert!(
                    !described.contains(&format!("\"-{flag}\"")),
                    "{program} lists -{flag}"
                );
            }
            assert!(
                described.contains("\"-export_search_strategy\""),
                "{program} does not list -export_search_strategy"
            );
        }
    }
}
