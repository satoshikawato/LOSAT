//! Web ABI v1's TBLASTN (`web_api::run_pair`, `run_pair_handles`), added on the Owner's request
//! of 2026-11-04 (gbdraw's serial module searches a protein query against its records); the
//! programs ABI v1 had before stay frozen (plan TD-1). The CLI's parser reads the words, and
//! the CLI's run (`TblastnArgs::run_in_memory`) reads the call's bytes as the files with the
//! same bytes: the report, the diagnostics and the errors are the command line's.
use super::TblastnArgs;
use anyhow::{bail, Result};

// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:114-117
/// The words of a request (`extra_args`, then `-outfmt` and the `outfmt` parameter) as the CLI's
/// parser reads them after `-query query -subject subject`: its defaults, value checks and errors.
/// The inputs and the output are the call's, so `-query`, `-subject`, `-out`, `-db` and the
/// help and version words are refused.
pub(crate) fn parse_args(extra: &[&str]) -> Result<TblastnArgs> {
    // The host has no database files; LOSAT Web searches the records it gives (`run_local`).
    if extra
        .iter()
        .any(|key| key.split('=').next().unwrap_or(key) == "-db")
    {
        bail!("-db (a database search) is not supported in the browser yet");
    }
    // The host has no PSSM files either (`-in_pssm`, line PS item PSW1).
    if extra
        .iter()
        .any(|key| key.split('=').next().unwrap_or(key) == "-in_pssm")
    {
        bail!("-in_pssm (a PSSM query) is not supported in the browser yet");
    }
    if extra.iter().any(|key| {
        matches!(
            key.split('=').next().unwrap_or(key),
            "-query"
                | "-subject"
                | "-out"
                | "-h"
                | "-help"
                | "-help-full"
                | "-xmlhelp"
                | "--help"
                | "-version"
                | "--version"
        )
    }) {
        bail!("file/help arguments are unsupported for the TBLASTN web API");
    }
    let mut argv = vec!["losat", "tblastn", "-query", "query", "-subject", "subject"];
    argv.extend_from_slice(extra);
    match crate::cli::try_parse_from::<crate::cli::Cli, _, _>(argv) {
        Ok(cli) => match cli.command {
            crate::cli::Commands::Tblastn(args) => Ok(args),
            _ => unreachable!("fixed TBLASTN subcommand"),
        },
        // The CLI's message for the words (`main` writes it to standard error).
        Err(error) => Err(anyhow::Error::msg(crate::cli::render_message(&error))),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:213-219
/// One search of `query` (protein FASTA) against `subject` (nucleotide FASTA) with `args`
/// (`parse_args`), as `losat tblastn` searches the files `query_label` and `subject_label` with
/// the same bytes, without NCBI's registry (`.ncbirc`, `tblastn.ini`, `NCBI_CONFIG__*`); the other
/// variables the CLI reads are read alike. The labels are only the names of the inputs (outfmt 0
/// and 7 show `Input: <subject_label>`; no TBLASTN report shows the query's).
pub(crate) fn run_pair(
    mut args: TblastnArgs,
    query: &[u8],
    subject: &[u8],
    query_label: &str,
    subject_label: &str,
) -> Result<Vec<u8>> {
    args.query = query_label.into();
    args.subject = Some(subject_label.into());
    // The report is returned, never written to a file (`parse_args` refuses `-out`).
    args.out = None;
    args.run_in_memory(Default::default(), query, subject)
}

/// The shared vector of ABI v1's TBLASTN (`tests/fixtures/web_v1_tblastn/cases.tsv`): the bytes
/// `losat tblastn` writes for each case (`tests/web_v1_tblastn.rs` checks them against the CLI),
/// which `run_pair` (the tests below) and the ABI's exports (`web_api.rs`) must return.
#[cfg(test)]
#[path = "web_vector_tests.rs"]
pub(crate) mod vector;

#[cfg(test)]
#[path = "web_tests.rs"]
mod tests;
