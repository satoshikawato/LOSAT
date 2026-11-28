//! dustmasker's `Run` (`CDustMaskApplication::Run`) on the arguments that NCBI's reader read
//! (`ncbi_app::args::Args`): the streams of `-out` and `-in`, then `work::run`.

use std::io::Write;

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::ncbi_app::args::Args;
use crate::ncbi_app::files::Cin;
use crate::ncbi_app::Diag;

use super::work::{self, report, Options, Report};

// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:275-276
/// The first step of `run`: `-out` opened (`AsOutputFile`, before `-in`), with its error as `run`
/// reports it (`Err` with the exit status): NCBI's `CArgException` of a name that does not open
/// (`File is not accessible`), which the application reports as an argument error, and LOSAT Web's
/// refusal of an output that it does not write (`OpenFiles::create`, windowmasker's place and
/// words). LOSAT Web's check takes it too (`check`).
fn output<'o>(
    args: &Args,
    out: &'o mut dyn Write,
    diag: &mut Diag,
) -> Result<crate::ncbi_app::files::OutputFile<'o>, i32> {
    let output = args
        .output_file("out", out)
        .map_err(|error| diag.run_arg_exception(error))?;
    if let Some(reason) = output.refusal() {
        return Err(diag.refuse(&work::unsupported(reason)));
    }
    Ok(output)
}

// No NCBI counterpart: LOSAT Web's check of the words (line NP WW4c; line NP WW4f, 2026-11-28).
/// `run` up to where it reads the input past its start: its first steps, `-out` opened (`output`),
/// then the writer made and `-in` opened and its reader made (`work::check`), whose failure is the
/// check's, as `run` reports it; 0 where they pass, as `run` then goes on to read the input. The
/// check's files are LOSAT Web's check's (`web_adapter::check_tool_in`), which make no file.
pub(super) fn check(
    args: &Args,
    settings: &ApplicationSettings,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    let output = match output(args, out, diag) {
        Ok(output) => output,
        Err(status) => return status,
    };
    let checked = work::check(&options(args, settings), stdin);
    drop(output);
    match checked {
        Ok(()) => 0,
        Err(failure) => reported(failure, diag),
    }
}

/// The options of `work::run` and `work::check` from `args` and `settings`.
fn options<'a>(args: &'a Args, settings: &'a ApplicationSettings) -> Options<'a> {
    let integer = |name| args.integer(name).expect("an integer key with a default");
    Options {
        input: args.os_value("in").expect("-in has a default"),
        infmt: args.bytes("infmt").expect("-infmt has a default"),
        outfmt: args.value("outfmt").expect("-outfmt has a default"),
        window: integer("window"),
        level: integer("level"),
        linker: integer("linker"),
        parse_seqids: args.flag("parse_seqids"),
        hard_masking: args.flag("hard_masking"),
        database: settings.database(),
    }
}

/// The exit status of the run's `failure`, with what the application writes for it (`report`).
fn reported(failure: work::Failure, diag: &mut Diag) -> i32 {
    match report(failure) {
        Report::Posts { lines, status } => {
            diag.write_text(&lines);
            status
        }
        Report::Refused(reason) => diag.refuse(&reason),
        Report::Arg(error) => diag.run_arg_exception(error),
    }
}

/// `CDustMaskApplication::Run` on `args`; the exit status (3 after an exception of `Run` that is
/// not an argument error, ncbiapp.cpp:711-734).
///
/// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:273-285
/// `stdin` is the standard input, which `-in -` reads (`work::run` opens `-in`).
pub(super) fn run(
    args: &Args,
    settings: &ApplicationSettings,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    let mut output = match output(args, out, diag) {
        Ok(output) => output,
        Err(status) => return status,
    };
    match work::run(&options(args, settings), stdin, &mut output, diag) {
        Ok(()) => 0,
        Err(failure) => reported(failure, diag),
    }
}
