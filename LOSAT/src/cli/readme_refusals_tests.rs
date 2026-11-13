//! No NCBI counterpart: a guard that README.md's section "What LOSAT refuses" names every option
//! that this parser refuses by name, per program, and every NCBI BLAST+ 2.17.0 executable that
//! LOSAT does not provide.
//!
//! The README is read at compile time. The lists under test are the ones that
//! `unknown_option_error` consults (`UNPORTED_*_ARGS`, `NCBI_TOOLKIT_ARGS`) and
//! `crate::writedb::makeblastdb::UNSUPPORTED_OPTIONS`. Refusals of a value (`-word_size 4`) or
//! of an input (a database sequence) are in no list and are not checked here.
//!
//! How the README is read (kept simple on purpose):
//! - The section is the text from `### What LOSAT refuses` to the next `---` or `## ` line.
//! - A line that starts with `**` opens a heading: the text between the first two pairs of `**`
//!   (`**psiblast** runs ...` is the heading `psiblast`). The lines that follow, up to the next
//!   heading, belong to it.
//! - The headings `psiblast`, `rpsblast`, `rpsbproc` and `makeblastdb` are the programs' own
//!   paragraphs. All other lines are general bullets.
//! - An option `-name` is named on a line when the line holds `-name` that is neither preceded
//!   nor followed by a letter, digit, `_` or `-` (so `-version` is not found in `-version-full`).
//! - A refusal of `-name` by a program is documented when its own paragraph names the option,
//!   or a general line names the option and the program (the program's name as a word, the
//!   same rule), or, for the five search programs, says `all five programs` or `all programs`.
//!   psiblast also takes what is documented for blastp when blastp's list holds the name, as its
//!   paragraph says ("It also refuses the `blastp` entries above that apply to it").
//! - The NCBI executables are covered when `ENTRY_NAMES` holds them or a backticked name in the
//!   headings `Programs LOSAT does not provide` and `Out of scope` is theirs.

use super::*;

const README: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../README.md"));

/// The programs that "all five programs" stands for in the README.
const FIVE_PROGRAMS: [&str; 5] = ["blastn", "blastp", "blastx", "tblastn", "tblastx"];

/// The programs whose README paragraph is their own.
const OWN_PARAGRAPHS: [&str; 4] = ["psiblast", "rpsblast", "rpsbproc", "makeblastdb"];

/// The programs that refuse the NCBI C++ Toolkit options by name (`unknown_option_error`).
const TOOLKIT_PROGRAMS: [&str; 10] = [
    "blastn",
    "blastp",
    "blastx",
    "tblastn",
    "tblastx",
    "psiblast",
    "rpsblast",
    "rpstblastn",
    "rpsbproc",
    "makeblastdb",
];

/// The executables in the `bin` directory of NCBI BLAST+ 2.17.0 (`ncbi-blast-2.17.0+-x64-linux`,
/// checked 2026-11-12).
const NCBI_EXECUTABLES: [&str; 28] = [
    "blast_formatter",
    "blast_formatter_vdb",
    "blast_vdb_cmd",
    "blastdb_aliastool",
    "blastdbcheck",
    "blastdbcmd",
    "blastn",
    "blastn_vdb",
    "blastp",
    "blastx",
    "cleanup-blastdb-volumes.py",
    "convert2blastmask",
    "deltablast",
    "dustmasker",
    "get_species_taxids.sh",
    "legacy_blast.pl",
    "makeblastdb",
    "makembindex",
    "makeprofiledb",
    "psiblast",
    "rpsblast",
    "rpstblastn",
    "segmasker",
    "tblastn",
    "tblastn_vdb",
    "tblastx",
    "update_blastdb.pl",
    "windowmasker",
];

/// The lines of the README's refusal section with the heading each belongs to ("" before the
/// first one).
fn refusal_lines() -> Vec<(&'static str, &'static str)> {
    let start = README
        .find("### What LOSAT refuses")
        .expect("README.md has the section `### What LOSAT refuses`");
    let mut heading = "";
    let mut lines = Vec::new();
    for line in README[start..].lines().skip(1) {
        if line.starts_with("---") || line.starts_with("## ") {
            break;
        }
        if let Some(rest) = line.strip_prefix("**") {
            heading = rest.split("**").next().unwrap_or("");
        }
        lines.push((heading, line));
    }
    assert!(
        !lines.is_empty(),
        "the refusal section of README.md is empty"
    );
    lines
}

fn is_name_char(c: char) -> bool {
    c.is_ascii_alphanumeric() || c == '_' || c == '-'
}

/// Whether `text` holds `word` with no name character on either side.
fn has_word(text: &str, word: &str) -> bool {
    text.match_indices(word).any(|(at, _)| {
        let before = text[..at].chars().next_back();
        let after = text[at + word.len()..].chars().next();
        !before.is_some_and(is_name_char) && !after.is_some_and(is_name_char)
    })
}

/// Whether `line` names the option `-name`.
fn names_option(line: &str, name: &str) -> bool {
    has_word(line, &format!("-{name}"))
}

/// Whether the README documents the refusal of `-name` by `program` (see the module comment).
fn documented(lines: &[(&str, &str)], program: &str, name: &str) -> bool {
    let own = OWN_PARAGRAPHS.contains(&program);
    let general = |line: &str| {
        has_word(line, program)
            || (FIVE_PROGRAMS.contains(&program)
                && (line.contains("all five programs") || line.contains("all programs")))
    };
    let direct = lines.iter().any(|&(heading, line)| {
        names_option(line, name)
            && ((own && heading == program)
                || (!OWN_PARAGRAPHS.contains(&heading) && general(line)))
    });
    direct
        || (program == "psiblast"
            && is_unported_blastp_arg(name)
            && documented(lines, "blastp", name))
}

#[test]
fn every_refused_option_is_in_the_readme_for_its_program() {
    let lines = refusal_lines();
    let psiblast: Vec<&str> = UNPORTED_PSIBLAST_OWN_ARGS
        .iter()
        .chain(UNPORTED_BLASTP_ARGS)
        .copied()
        .filter(|name| is_unported_psiblast_arg(name))
        .collect();
    let mut lists: Vec<(&str, Vec<&str>)> = vec![
        ("blastn", UNPORTED_BLASTN_ARGS.to_vec()),
        ("blastp", UNPORTED_BLASTP_ARGS.to_vec()),
        ("blastx", UNPORTED_BLASTX_ARGS.to_vec()),
        ("tblastn", UNPORTED_TBLASTN_ARGS.to_vec()),
        ("tblastx", UNPORTED_TBLASTX_ARGS.to_vec()),
        ("psiblast", psiblast),
        ("rpsblast", UNPORTED_RPSBLAST_ARGS.to_vec()),
        (
            "makeblastdb",
            crate::writedb::makeblastdb::UNSUPPORTED_OPTIONS.to_vec(),
        ),
    ];
    for program in TOOLKIT_PROGRAMS {
        lists.push((program, NCBI_TOOLKIT_ARGS.to_vec()));
    }
    let mut missing = Vec::new();
    for (program, names) in &lists {
        for name in names {
            if !documented(&lines, program, name) {
                missing.push(format!(
                    "{program}: the option -{name} is refused by LOSAT but README.md's \"What LOSAT refuses\" does not name it for {program}"
                ));
            }
        }
    }
    assert!(missing.is_empty(), "\n{}\n", missing.join("\n"));
}

#[test]
fn the_documented_toolkit_programs_refuse_the_toolkit_options() {
    let command = Cli::command();
    let mut wrong = Vec::new();
    for program in TOOLKIT_PROGRAMS {
        let scope = command
            .find_subcommand(program)
            .unwrap_or_else(|| panic!("{program} is a subcommand"));
        for name in NCBI_TOOLKIT_ARGS {
            let message = unknown_option_error(scope, name, &format!("-{name}")).to_string();
            if !message.contains("NCBI C++ Toolkit option") {
                wrong.push(format!("{program}: -{name} gives {message:?}"));
            }
        }
    }
    assert!(
        wrong.is_empty(),
        "the README says these programs refuse the NCBI C++ Toolkit options by name:\n{}",
        wrong.join("\n")
    );
}

// `-version` and the full versions are answered on the command line before the parser
// (`ncbi_version_request`; blastdbcmd, dustmasker and segmasker in `ncbi_app::run_with`), so
// they are no refusal and the README lists none of them (line RG item E3).
#[test]
fn the_versions_are_answered_by_every_program() {
    for program in TOOLKIT_PROGRAMS {
        for word in [
            "-version",
            "-version-full",
            "-version-full-xml",
            "-version-full-json",
        ] {
            assert!(
                ncbi_version_request(["losat", program, word]).is_some(),
                "{program} {word} is not answered, so it is a refusal that README.md must list"
            );
        }
    }
}

#[test]
fn every_ncbi_executable_is_provided_or_listed_as_not_provided() {
    let lines = refusal_lines();
    let listed: Vec<&str> = lines
        .iter()
        .filter(|(heading, _)| {
            matches!(*heading, "Programs LOSAT does not provide" | "Out of scope")
        })
        .map(|&(_, line)| line)
        .collect();
    let backticked = |line: &str| -> Vec<String> {
        line.split('`')
            .skip(1)
            .step_by(2)
            .map(str::to_string)
            .collect()
    };
    let not_provided: Vec<String> = lines
        .iter()
        .filter(|(heading, _)| *heading == "Programs LOSAT does not provide")
        .flat_map(|&(_, line)| backticked(line))
        .collect();
    assert!(
        !not_provided.is_empty(),
        "README.md's \"Programs LOSAT does not provide\" lists nothing"
    );
    let mut wrong = Vec::new();
    for executable in NCBI_EXECUTABLES {
        let provided = crate::cli::ENTRY_NAMES.contains(&executable);
        let named = listed
            .iter()
            .any(|line| backticked(line).iter().any(|name| name == executable));
        if provided && not_provided.iter().any(|name| name == executable) {
            wrong.push(format!(
                "{executable}: LOSAT provides it but README.md lists it under \"Programs LOSAT does not provide\""
            ));
        }
        if !provided && !named {
            wrong.push(format!(
                "{executable}: an executable of NCBI BLAST+ 2.17.0 that LOSAT does not provide and README.md does not list under \"Programs LOSAT does not provide\""
            ));
        }
    }
    for name in &not_provided {
        if !NCBI_EXECUTABLES.contains(&name.as_str()) {
            wrong.push(format!(
                "{name}: listed under \"Programs LOSAT does not provide\" but not an executable of NCBI BLAST+ 2.17.0 (NCBI_EXECUTABLES)"
            ));
        }
    }
    assert!(wrong.is_empty(), "\n{}\n", wrong.join("\n"));
}
