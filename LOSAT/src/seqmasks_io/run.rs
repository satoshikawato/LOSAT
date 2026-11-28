//! The parts of a masking program's `Run` that dustmasker and segmasker share around their
//! reader and writer: what a failed write, a failed read and a closed pipe become, and the end
//! of the run (the writer's object, the serial stream's report of a failed write). Each program keeps its own failure type ([`RunFailure`]), as
//! NCBI's programs report their exceptions differently (dustmasker through `AppMain`,
//! segmasker in the `catch` of its `Run`).

use crate::ncbi_app::files::{ends_like_sigpipe, OutputFile};
use crate::ncbi_app::Diag;

use super::fasta_reader::ios_failure_what;
use super::writer::{writes_serial_objects, MaskWriter, WriteError};

/// The failures of a masking program's run that the shared parts make.
pub(crate) trait RunFailure: Sized {
    /// LOSAT's refusal: the whole message.
    fn refused(message: String) -> Self;
    /// A write found the pipe closed where SIGPIPE ends NCBI's process (`OutputFile::closed`):
    /// the run ends at once, exit status 0, no message.
    fn closed_pipe() -> Self;
    /// A `std::exception` that leaves NCBI's `Run` (`AppMain`'s error post `ERR_POST_X` with
    /// subcode 7, exit status 3), with its `what()`.
    fn std_exception(what: String) -> Self;
    /// Whether the run stops writing at this failure: LOSAT's refusal comes before NCBI's run
    /// would end, and a closed pipe ends it where SIGPIPE ends NCBI's.
    fn stops_writing(&self) -> bool;
}

/// LOSAT's refusal of `what` (a noun phrase) in `program` (its name in refusals).
pub(crate) fn unsupported(what: &str, program: &str) -> String {
    format!("{what} is not supported by LOSAT's {program}")
}

// No NCBI counterpart: LOSAT Web's validation of a masking program's words (line NP WW4c,
// 2026-11-28, as windowmasker's `app::check`, D-WW3-7 and D-WW3-10).
/// The check of `-out` in a run of dustmasker, segmasker or convert2blastmask (`program`, its name
/// in refusals) that reaches `-out`: an `-out` that names a file which LOSAT Web does not write
/// (`RunFiles::refused_output`) is the refusal that the run makes when it opens `-out`. 0
/// otherwise; the standard output `-` is not refused. dustmasker opens `-out` first; segmasker and
/// convert2blastmask open `-in` and make its reader first, and call this only when those steps pass
/// (review of WW4c, F1).
pub(crate) fn check_output(
    args: &crate::ncbi_app::args::Args,
    diag: &mut Diag,
    program: &str,
) -> i32 {
    let out = args.os_value("out").unwrap_or(std::ffi::OsStr::new("-"));
    if out != "-" {
        if let Some(reason) = args.files().refused_output(out) {
            return diag.refuse(&unsupported(&reason, program));
        }
    }
    0
}

/// A failed write of the writer: a closed pipe, LOSAT's refusal of what it cannot write.
pub(crate) fn write_failure<F: RunFailure>(error: WriteError, program: &str) -> F {
    match error {
        WriteError::Io(error) if ends_like_sigpipe(&error) => F::closed_pipe(),
        WriteError::Io(error) => F::refused(format!("writing the output failed: {error}")),
        WriteError::Refused(what) => F::refused(unsupported(&what, program)),
    }
}

/// A failed read of `-in`: the `std::ios_base::failure` of NCBI's stream
/// (`fasta_reader::ios_failure_what`), which leaves NCBI's `Run`.
pub(crate) fn read_failure<F: RunFailure>(error: std::io::Error, program: &str) -> F {
    match ios_failure_what(&error) {
        Some(what) => F::std_exception(what),
        None => F::refused(unsupported(
            &format!("an input that cannot be read ({error})"),
            program,
        )),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_bdb_reader.cpp:44-49
// `GetBioseq` posts the taxonomy database's texts (`ERR_POST`, severity Error, no prefix in
// dustmasker and segmasker) while it reads a sequence, before the writer prints it; `cerr` is tied
// to `cout`, so each text comes after NCBI's writes of the output so far, and on a closed pipe
// SIGPIPE ends NCBI's process at that write, before the text.
/// Posts `texts`, the taxonomy database's texts of reading a database sequence
/// (`bdb_reader::BdbRecord::texts`), before the sequence is written: `Error: <text>` to the
/// diagnostics (`-logfile`), NCBI's writes of the output made first (`Diag`'s tie); a closed pipe
/// ends the run there, without the text.
pub(crate) fn post_read_texts<F: RunFailure>(texts: &[String], diag: &mut Diag) -> Result<(), F> {
    for text in texts {
        diag.err_post(text);
        if diag.closed() {
            return Err(F::closed_pipe());
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:338-353
/// LOSAT's refusal when the output `out` of a `seqloc_*` or `maskinfo_*` writer failed. Those
/// write through NCBI's serial stream (`CObjectOStream`), which checks its writes: NCBI posts
/// `cannot close output stream` for each object, then stops on `write fault` at the next one,
/// or aborts (134) for a `maskinfo` object larger than its 8 KiB buffer (dustmasker exits with
/// 3, 0 or 134 by the case, segmasker with 1 or 134). The `interval`, `acclist` and `fasta`
/// writers write to the `ostream` itself and go on silently, as here.
pub(crate) fn serial_write_failure<F: RunFailure>(
    out: &OutputFile<'_>,
    outfmt: &str,
    program: &str,
) -> Result<(), F> {
    if out.closed() {
        return Err(F::closed_pipe());
    }
    if !writes_serial_objects(outfmt) {
        return Ok(());
    }
    // NCBI's writes so far are made: whether one failed.
    let bad = match out.bad() {
        Ok(bad) => bad,
        Err(error) if ends_like_sigpipe(&error) => return Err(F::closed_pipe()),
        Err(_) => true,
    };
    if bad {
        return Err(F::refused(unsupported(
            &format!(
                "a failed write of the -outfmt {outfmt} output (NCBI's serial stream reports it)"
            ),
            program,
        )));
    }
    Ok(())
}

/// The end of the run whose masking ended with `result`: the writer writes its object (a
/// `maskinfo_*` writer's destructor writes it also when an exception leaves NCBI's `Run`)
/// unless the run stopped writing; then the failure, or the serial stream's report of a failed
/// write. The diagnostics of a failure come after NCBI's writes of the output (`Diag`'s tie; a
/// closed pipe ends the run there, without them), and the stream ends after the run
/// (`OpenFiles::finish_streams`).
pub(crate) fn finish<F: RunFailure>(
    writer: MaskWriter,
    out: &mut OutputFile<'_>,
    diag: &mut Diag,
    mut result: Result<(), F>,
    outfmt: &str,
    program: &str,
) -> Result<(), F> {
    if !matches!(&result, Err(failure) if failure.stops_writing()) {
        if let Err(error) = writer.finish(out, diag) {
            result = Err(write_failure(error, program));
        }
    }
    if out.closed() {
        return Err(F::closed_pipe());
    }
    result?;
    serial_write_failure(out, outfmt, program)
}
