//! No NCBI counterpart: a guard that the user manual's page "What LOSAT refuses"
//! (docs/manual/src/refusals.md) names every option and run-time refusal per program and every
//! NCBI BLAST+ 2.17.0 executable that LOSAT does not provide, while rare and pathological inputs
//! live on the linked page docs/manual/src/refusals-edge-cases.md (not counted toward CLI option
//! coverage).
//!
//! The page is read at compile time. The list under test is the one list of the arguments that
//! LOSAT refuses by name or at the run (`ncbi_app::args::refusals`, the rows of `REFUSED` in
//! `build/ncbi_args.rs`). Refusals of a value (`-word_size 4`) or of an input (a database
//! sequence) are in no list and are not checked here.
//!
//! How the page is read (kept simple on purpose):
//! - A line that starts with `## ` opens a section; its heading is the rest of the line. The lines
//!   that follow, up to the next `## ` line, belong to it.
//! - A section is a program's own when the first word of its heading, without backticks, is the
//!   program (`## tblastn with `-in_pssm`` is tblastn's).
//! - An option `-name` is named on a line when the line holds `-name` that is neither preceded
//!   nor followed by a letter, digit, `_` or `-` (so `-version` is not found in `-version-full`).
//! - A refusal of `-name` by a program is documented when a line of the program's own section
//!   names the option, or, for the five search programs, a line of the section `All five search
//!   programs` does, or any line names the option and the program (the program's name as a word,
//!   the same rule) or says `all programs` (or, for the five search programs, `all five search
//!   programs`).
//! - The NCBI executables are covered when `ENTRY_NAMES` holds them or a backticked name in the
//!   sections `Programs LOSAT does not provide`, `On the list for a later release` and `Out of
//!   scope...` is theirs.
//! - blast_formatter's refusals (`crate::blast_formatter::refusals`, one list) are documented when
//!   the lines of the section `blast_formatter` hold the words of each one.

use super::*;

const PAGE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../docs/manual/src/refusals.md"
));

/// The programs of the section `All five search programs`.
const FIVE_PROGRAMS: [&str; 5] = ["blastn", "blastp", "blastx", "tblastn", "tblastx"];

/// The heading of the section that holds the refusals of the five search programs.
const FIVE_SECTION: &str = "All five search programs";

/// The programs that `try_parse_from` reads, which answer the NCBI C++ Toolkit options.
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

/// The lines of one manual page with the heading of the section each belongs to
/// ("" before the first `## ` heading).
fn page_lines(page: &'static str) -> Vec<(&'static str, &'static str)> {
    let mut lines = Vec::new();
    let mut heading = "";
    for line in page.lines() {
        if let Some(rest) = line.strip_prefix("## ") {
            heading = rest.trim();
        }
        lines.push((heading, line));
    }
    lines
}

/// The lines of `docs/manual/src/refusals.md` (primary refused options, run refusals, and
/// non-provided executable list; `refusals-edge-cases.md` is intentionally excluded).
fn primary_refusal_lines() -> Vec<(&'static str, &'static str)> {
    let lines = page_lines(PAGE);
    assert!(
        lines.iter().any(|&(heading, _)| heading == FIVE_SECTION),
        "refusals.md has no section `## {FIVE_SECTION}`"
    );
    lines
}

/// The program whose own section `heading` opens: the heading's first word without backticks.
fn section_program(heading: &str) -> &str {
    heading
        .split_whitespace()
        .next()
        .unwrap_or("")
        .trim_matches('`')
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

/// Whether the page documents the refusal of `-name` by `program` (see the module comment).
fn documented(lines: &[(&str, &str)], program: &str, name: &str) -> bool {
    let five = FIVE_PROGRAMS.contains(&program);
    lines.iter().any(|&(heading, line)| {
        names_option(line, name)
            && (section_program(heading) == program
                || (five && heading == FIVE_SECTION)
                || has_word(line, program)
                || line.contains("all programs")
                || (five && line.contains("all five search programs")))
    })
}

#[test]
fn every_refused_option_is_in_the_manual_for_its_program() {
    let lines = primary_refusal_lines();
    let mut missing = Vec::new();
    for (program, name, refusal) in crate::ncbi_app::args::refusals() {
        if matches!(refusal, Refusal::Read(_)) && !documented(&lines, program, name) {
            missing.push(format!(
                "{program}: the option -{name} is refused by LOSAT but docs/manual/src/refusals.md does not name it for {program} (rare inputs in refusals-edge-cases.md do not count)"
            ));
        }
    }
    assert!(missing.is_empty(), "\n{}\n", missing.join("\n"));
}

#[test]
fn every_run_refusal_is_in_the_manual_for_its_program() {
    let lines = primary_refusal_lines();
    let mut missing = Vec::new();
    for (program, name, refusal) in crate::ncbi_app::args::refusals() {
        if matches!(refusal, Refusal::Run(_)) && !documented(&lines, program, name) {
            missing.push(format!(
                "{program}: the run refusal -{name} is refused by LOSAT's run but docs/manual/src/refusals.md does not name it for {program} (rare inputs in refusals-edge-cases.md do not count)"
            ));
        }
    }
    assert!(missing.is_empty(), "\n{}\n", missing.join("\n"));
}

// `-version` and the full versions are answered on the command line before the parser
// (`ncbi_version_request`; blastdbcmd, dustmasker and segmasker in `ncbi_app::run_with`), so
// they are no refusal and the manual lists none of them (line RG item E3).
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
                "{program} {word} is not answered, so it is a refusal that refusals.md must list"
            );
        }
    }
}

#[test]
fn every_ncbi_executable_is_provided_or_listed_as_not_provided() {
    let lines = primary_refusal_lines();
    let listed: Vec<&str> = lines
        .iter()
        .filter(|(heading, _)| {
            *heading == "Programs LOSAT does not provide"
                || *heading == "On the list for a later release"
                || heading.starts_with("Out of scope")
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
        "docs/manual/src/refusals.md's \"Programs LOSAT does not provide\" lists nothing"
    );
    let mut wrong = Vec::new();
    for executable in NCBI_EXECUTABLES {
        let provided = crate::cli::ENTRY_NAMES.contains(&executable);
        let named = listed
            .iter()
            .any(|line| backticked(line).iter().any(|name| name == executable));
        if provided && not_provided.iter().any(|name| name == executable) {
            wrong.push(format!(
                "{executable}: LOSAT provides it but refusals.md lists it under \"Programs LOSAT does not provide\""
            ));
        }
        if !provided && !named {
            wrong.push(format!(
                "{executable}: an executable of NCBI BLAST+ 2.17.0 that LOSAT does not provide and refusals.md does not list under \"Programs LOSAT does not provide\""
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

#[test]
fn every_blast_formatter_refusal_is_in_its_manual_section() {
    let paragraph: Vec<&str> = primary_refusal_lines()
        .into_iter()
        .filter(|(heading, _)| *heading == "blast_formatter")
        .map(|(_, line)| line)
        .collect();
    assert!(
        !paragraph.is_empty(),
        "the manual's \"What LOSAT refuses\" has no section `## blast_formatter`"
    );
    let mut missing = Vec::new();
    for refusal in crate::blast_formatter::refusals::documented() {
        let words = refusal.manual_words();
        if !paragraph.iter().any(|line| line.contains(&words)) {
            missing.push(format!(
                "blast_formatter: {refusal:?} ({words}) is refused by LOSAT but the manual's section `## blast_formatter` does not name it"
            ));
        }
    }
    assert!(missing.is_empty(), "\n{}\n", missing.join("\n"));
}

/// The `-in_msa` refusals of an input that the code names by a cause (each `pub(crate)` constant
/// is the text after "is not supported by LOSAT's PSIBLAST:"), and a phrase of the row of
/// refusals.md that tells the user of it. The page's rows of `-in_msa` input refusals are those of
/// this list and the three that spell their cause where they refuse (a master sequence that is the
/// start of a longer one above it, one of X or O alone, `-matrix IDENTITY`): a refusal added to the
/// code needs its row, and a row needs its refusal.
#[test]
fn the_in_msa_input_refusals_of_the_code_are_listed_in_the_manual() {
    let refusals: [(&str, &str); 5] = [
        (
            crate::blastinput::aln_reader::phylip::ZERO_COUNT,
            "a PHYLIP header of 0 sequences",
        ),
        (
            crate::blastinput::aln_reader::nexus::EMPTY_COMMAND,
            "an empty command (a lone `;`) in a NEXUS file",
        ),
        (
            crate::algorithm::psiblast::in_msa::SHORT_ROW_READ,
            "a FASTA alignment whose rows have different lengths",
        ),
        (
            crate::blastinput::aln_reader::mods::BACKWARD_RANGE,
            "`secondary-accession` range that runs backwards",
        ),
        (
            crate::blastinput::aln_reader::mods::RANGE_OVERFLOW,
            "a number of more than 18 digits",
        ),
    ];
    let mut wrong = Vec::new();
    for (cause, phrase) in refusals {
        let row = PAGE
            .lines()
            .find(|line| line.starts_with('|') && line.contains(phrase));
        match row {
            None => wrong.push(format!("refusals.md has no row with \"{phrase}\"")),
            Some(row) if !row.contains("Refused") => {
                wrong.push(format!("the row with \"{phrase}\" does not say Refused"))
            }
            Some(_) if cause.is_empty() => wrong.push(format!("\"{phrase}\": the cause is empty")),
            Some(_) => {}
        }
    }
    assert!(wrong.is_empty(), "\n{}\n", wrong.join("\n"));
}
