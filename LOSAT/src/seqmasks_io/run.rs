//! The parts of a masking program's `Run` that dustmasker and segmasker share around their
//! reader and writer: what a failed write, a failed read and a closed pipe become, and the end
//! of the run (the writer's object, the flush before the diagnostics, the serial stream's
//! report of a failed write). Each program keeps its own failure type ([`RunFailure`]), as
//! NCBI's programs report their exceptions differently (dustmasker through `AppMain`,
//! segmasker in the `catch` of its `Run`).

use std::io::Write;

use crate::ncbi_app::files::{ends_like_sigpipe, Ostream};
use crate::ncbi_app::Diag;

use super::fasta_reader::ios_failure_what;
use super::writer::{writes_serial_objects, MaskWriter, WriteError};

/// The failures of a masking program's run that the shared parts make.
pub(crate) trait RunFailure: Sized {
    /// LOSAT's refusal: the whole message.
    fn refused(message: String) -> Self;
    /// A write found the pipe closed where SIGPIPE ends NCBI's process (`Ostream::closed`): the
    /// run ends at once, exit status 0, no message.
    fn closed_pipe() -> Self;
    /// A `std::exception` that leaves NCBI's `Run` (`AppMain`'s `ERR_POST_X(7, ...)`, exit
    /// status 3), with its `what()`.
    fn std_exception(what: String) -> Self;
    /// Whether the run stops writing at this failure: LOSAT's refusal comes before NCBI's run
    /// would end, and a closed pipe ends it where SIGPIPE ends NCBI's.
    fn stops_writing(&self) -> bool;
}

/// LOSAT's refusal of `what` (a noun phrase) in `program` (its name in refusals).
pub(crate) fn unsupported(what: &str, program: &str) -> String {
    format!("{what} is not supported by LOSAT's {program}")
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

// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:338-353
// ```c++
// void CObjectOStream::Close(void)
// {
//     if (m_Fail != fNotOpen) {
//         try {
//             DefaultFlush();
//             ...
//             m_Output.Close();
//         }
//         catch (CException& exc) {
//             if ( InGoodState() )
//                 RethrowError(fWriteError, "cannot close output stream",exc);
//         }
// ```
/// LOSAT's refusal when the output `out` of a `seqloc_*` or `maskinfo_*` writer failed. Those
/// write through NCBI's serial stream (`CObjectOStream`), which checks its writes: NCBI posts
/// `cannot close output stream` for each object, then stops on `write fault` at the next one,
/// or aborts (134) for a `maskinfo` object larger than its 8 KiB buffer (dustmasker exits with
/// 3, 0 or 134 by the case, segmasker with 1 or 134). The `interval`, `acclist` and `fasta`
/// writers write to the `ostream` itself and go on silently, as here.
pub(crate) fn serial_write_failure<F: RunFailure>(
    out: &Ostream<'_>,
    outfmt: &str,
    program: &str,
) -> Result<(), F> {
    if out.closed() {
        return Err(F::closed_pipe());
    }
    if out.bad() && writes_serial_objects(outfmt) {
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
/// unless the run stopped writing; the output is flushed (NCBI's `cerr` is tied to `cout`, so
/// the output comes before the diagnostics of a failure, and a closed pipe ends the run there,
/// without them); then the failure, or the serial stream's report of a failed write.
pub(crate) fn finish<F: RunFailure>(
    writer: MaskWriter,
    out: &mut Ostream<'_>,
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
    let _ = out.flush();
    if out.closed() {
        return Err(F::closed_pipe());
    }
    result?;
    serial_write_failure(out, outfmt, program)
}
