//! dustmasker's `Run` (`CDustMaskApplication::Run`): the work on the arguments that NCBI's reader
//! read (`ncbi_app::args::Args`).
//!
//! STUB of the command-line package (line RG, item K1): it does the steps of NCBI's `Run` that
//! end in an argument error or in an error before any sequence is read (the streams of `-out`
//! and `-in`, the check of `-hard_masking`), then returns 0 and writes nothing. The writer, the
//! reader (`-infmt blastdb` and an unknown `-infmt` included) and the masking come with the
//! body.
//!
//! What the body gets (`ncbi_app::Body`):
//! - `args`: `flag(name)` (`-parse_seqids`, `-hard_masking`), `integer(name)` (the `i32` of
//!   `-window`, `-level`, `-linker`, their defaults applied: NCBI's `AsInteger()`), `value(name)`
//!   (the text of `-infmt`, `-outfmt`, `-in`, `-out`, defaults applied), `given(name)`, and the
//!   streams `output_file("out", out)` and `input_file("in")` (`AsOutputFile()`,
//!   `AsInputFile()`: `-` is the standard output or input; `Err` is NCBI's `CArgException`, which
//!   goes to `diag.run_arg_exception`). With `-infmt blastdb`, `-in` is a database name
//!   (`args.value("in")`), not opened as a file.
//! - `settings`: the application's settings (the `BLASTDB` of the registry for `-infmt blastdb`).
//! - `out`: the standard output (`-out -`).
//! - `diag`: the diagnostics, on the target of `-logfile`: `write_text` for NCBI's reports made
//!   by `crate::ncbi_exception` (`err_post_x(106, 7, ...)` of a `std::exception`,
//!   `NcbiException::report_default_ex(106, 16, ...)` of a `CException`), `post` for one line.

use std::io::Write;

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::ncbi_app::args::Args;
use crate::ncbi_app::Diag;

/// `CDustMaskApplication::Run` on `args`; the exit status (3 after an exception of `Run` that is
/// not an argument error, ncbiapp.cpp:711-734).
pub(super) fn run(
    args: &Args,
    _settings: &ApplicationSettings,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    // NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:273-285
    // ```c++
    // int CDustMaskApplication::Run (void)
    // {
    //     // Set up the input and output streams.
    //     CNcbiOstream& output_stream = GetArgs()[kOutput].AsOutputFile();
    //
    //     // Set up the object manager.
    //     CRef<CObjectManager> om(CObjectManager::GetInstance());
    //
    //     // Set up the duster object.
    //     Uint4 level = GetArgs()["level"].AsInteger();
    //     duster_type::size_type window = GetArgs()["window"].AsInteger();
    //     duster_type::size_type linker = GetArgs()["linker"].AsInteger();
    // ```
    let _output = match args.output_file("out", out) {
        Ok(output) => output,
        Err(error) => return diag.run_arg_exception(error),
    };
    let _level = args.integer("level");
    let _window = args.integer("window");
    let _linker = args.integer("linker");
    // NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:125-131
    // ```c++
    //     const CArgs& args = GetArgs();
    //     const string& format(args[kOutputFormat].AsString());
    //     CMaskWriter* retval = NULL;
    //
    //     if (args["hard_masking"].AsBoolean() && (format != "fasta")) {
    //         throw runtime_error("Hard masking can only be applied for fasta output");
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:730-734
    // ```c++
    //             catch (const exception& e) {
    //                 ERR_POST_X(7, "Application's execution failed: " << e.what());
    //                 *got_exception = true;
    //                 *exit_code = 3;
    //             }
    // ```
    if args.flag("hard_masking") && args.value("outfmt") != Some("fasta") {
        diag.write_text(
            crate::ncbi_exception::err_post_x(
                106,
                7,
                "Application's execution failed: Hard masking can only be applied for fasta output",
            )
            .as_bytes(),
        );
        return 3;
    }
    // NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp:170-177
    // ```c++
    //     const CArgs & args( GetArgs() );
    //     const string & format( args[kInputFormat].AsString() );
    //
    //     if( format == "fasta" ) {
    //         CNcbiIstream& input_stream  = GetArgs()[kInput].AsInputFile();
    //         return new CMaskFastaReader(
    //                 input_stream, true, args["parse_seqids"] );
    //     }
    // ```
    if args.value("infmt") == Some("fasta") {
        let _input = match args.input_file("in") {
            Ok(input) => input,
            Err(error) => return diag.run_arg_exception(error),
        };
    }
    0
}
