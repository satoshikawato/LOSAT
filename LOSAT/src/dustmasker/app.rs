//! dustmasker's `Run` (`CDustMaskApplication::Run`) on the arguments that NCBI's reader read
//! (`ncbi_app::args::Args`): the streams of `-out` and `-in`, then `work::run`.

use std::io::Write;

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::ncbi_app::args::Args;
use crate::ncbi_app::files::Cin;
use crate::ncbi_app::Diag;

use super::work::{self, report, Options, Report};

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
    let mut output = match args.output_file("out", out) {
        Ok(output) => output,
        Err(error) => return diag.run_arg_exception(error),
    };
    let integer = |name| args.integer(name).expect("an integer key with a default");
    let options = Options {
        input: args.value("in").expect("-in has a default"),
        infmt: args.value("infmt").expect("-infmt has a default"),
        outfmt: args.value("outfmt").expect("-outfmt has a default"),
        window: integer("window"),
        level: integer("level"),
        linker: integer("linker"),
        parse_seqids: args.flag("parse_seqids"),
        hard_masking: args.flag("hard_masking"),
        blastdb_registry: match &settings.blastdb {
            Ok(value) => Ok(value.as_deref()),
            Err(reason) => Err(reason.as_str()),
        },
    };
    match work::run(&options, stdin, &mut output, diag) {
        Ok(()) => 0,
        Err(failure) => match report(failure) {
            Report::Posts { lines, status } => {
                diag.write_text(&lines);
                status
            }
            Report::Refused(reason) => diag.refuse(&reason),
            Report::Arg(error) => diag.run_arg_exception(error),
        },
    }
}
