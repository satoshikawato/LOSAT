//! The process's standard output as LOSAT's command line writes it: with no buffer of Rust's.
//!
//! Rust's `std::io::Stdout` is line-buffered. A write that fails leaves up to 1023 bytes of it in
//! that buffer, and a later write or the flush at the exit writes them once the failure has
//! passed (`EAGAIN` on a non-blocking standard output, a device that gets space back). NCBI's
//! programs lose those bytes: a failed write leaves glibc's buffer of `cout` empty
//! (`crate::writedb::cout`). [`standard_output`] passes each write to the device at once and
//! keeps nothing (line RG item KC-11b-2). A write's error is the device's (a full device, a
//! descriptor opened only for reading), as NCBI's write gets it.
//!
//! What a failed write loses is up to the writer above. The stream states of NCBI
//! (`crate::writedb::cout`, `crate::ncbi_app::files::Ostream`, rpsbproc's report) drop the
//! bytes of the failed write and write nothing after it. Where that cuts the output follows
//! each writer's own write sizes; NCBI's for the tools' records are item T8a's (line RG-T5).
//! A `std::io::BufWriter` keeps the bytes. The search programs' report stream
//! (`crate::report::ostream::ReportStream`) holds NCBI's buffer of `cout` over this device.
//!
//! The library leaves the process's standard output to `Stdout`, which its caller shares.

/// The process's standard output, written with no buffer ([`standard_output`]).
pub(crate) struct StandardOutput(Device);

/// Where [`StandardOutput`] writes.
enum Device {
    /// Descriptor 1 (`STD_OUTPUT_HANDLE` on Windows) viewed as a file and never closed.
    #[cfg(any(unix, windows, target_os = "wasi"))]
    Raw(std::mem::ManuallyDrop<std::fs::File>),
    /// A Windows console (`GetConsoleMode` succeeds on the handle). Rust's `Stdout` writes to it
    /// with the console's wide-character call: UTF-8 is converted to UTF-16, and bytes that are
    /// not UTF-8 are refused (`InvalidData`). This flushes it after each write, so its line
    /// buffer holds nothing between writes.
    #[cfg(windows)]
    Console(std::io::Stdout),
    /// No standard output: a Windows process started without the handle, or a target that has
    /// no standard output. Writes succeed and go nowhere, as they do through Rust's `Stdout`.
    #[cfg(not(any(unix, target_os = "wasi")))]
    Missing,
}

/// The process's standard output, written with no buffer ([`StandardOutput`]). Rust's
/// standard library writes its `Stdout` the same way under its line buffer: descriptor 1 in a
/// `ManuallyDrop`, or on Windows the handle, with a console written as UTF-16.
/// `blastinput::input_files` reads the standard input the same way.
pub(crate) fn standard_output() -> StandardOutput {
    #[cfg(any(unix, target_os = "wasi"))]
    {
        use std::os::fd::{AsRawFd, FromRawFd};
        // SAFETY: descriptor 1 stays open for the whole process. On Unix the Rust runtime opens
        // `/dev/null` there before `main` if it was closed; on WASI it is the host's standard
        // output. LOSAT never closes it. The `File` sits in a `ManuallyDrop` that is never
        // dropped, so it only borrows the descriptor and never closes it.
        let file = unsafe { std::fs::File::from_raw_fd(std::io::stdout().as_raw_fd()) };
        StandardOutput(Device::Raw(std::mem::ManuallyDrop::new(file)))
    }
    #[cfg(windows)]
    {
        use std::os::windows::io::{AsRawHandle, FromRawHandle};
        let stdout = std::io::stdout();
        let handle = stdout.as_raw_handle();
        if handle.is_null() {
            return StandardOutput(Device::Missing);
        }
        let mut mode = 0u32;
        // SAFETY: `GetConsoleMode` only reads the handle it is given and writes the `u32` it is
        // given. The handle is the process's standard output, open for the call. The call
        // fails (returns 0) on a handle that is not a console, including the pipes of an
        // msys or cygwin terminal.
        let console =
            unsafe { windows_sys::Win32::System::Console::GetConsoleMode(handle, &mut mode) != 0 };
        if console {
            return StandardOutput(Device::Console(stdout));
        }
        // SAFETY: the standard output handle stays open for the whole process and LOSAT never
        // closes it. The `File` sits in a `ManuallyDrop` that is never dropped, so it only
        // borrows the handle and never closes it.
        let file = unsafe { std::fs::File::from_raw_handle(handle) };
        StandardOutput(Device::Raw(std::mem::ManuallyDrop::new(file)))
    }
    #[cfg(not(any(unix, windows, target_os = "wasi")))]
    {
        StandardOutput(Device::Missing)
    }
}

impl std::io::Write for StandardOutput {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match &mut self.0 {
            #[cfg(any(unix, windows, target_os = "wasi"))]
            Device::Raw(file) => std::io::Write::write(&mut **file, buf),
            #[cfg(windows)]
            Device::Console(stdout) => {
                std::io::Write::write_all(stdout, buf)?;
                std::io::Write::flush(stdout)?;
                Ok(buf.len())
            }
            #[cfg(not(any(unix, target_os = "wasi")))]
            Device::Missing => Ok(buf.len()),
        }
    }

    /// The device's own gathered write (`writev` of the descriptor), so that the bytes that NCBI's
    /// buffer holds and the bytes after them reach it as one write where NCBI's stream writes them
    /// so (`report::ostream`'s device writes). A console writes each part.
    fn write_vectored(&mut self, bufs: &[std::io::IoSlice<'_>]) -> std::io::Result<usize> {
        match &mut self.0 {
            #[cfg(any(unix, windows, target_os = "wasi"))]
            Device::Raw(file) => std::io::Write::write_vectored(&mut **file, bufs),
            #[cfg(windows)]
            Device::Console(stdout) => {
                let Some(first) = bufs.iter().find(|buf| !buf.is_empty()) else {
                    return Ok(0);
                };
                std::io::Write::write_all(stdout, first)?;
                std::io::Write::flush(stdout)?;
                Ok(first.len())
            }
            #[cfg(not(any(unix, target_os = "wasi")))]
            Device::Missing => Ok(bufs.iter().map(|buf| buf.len()).sum()),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

#[cfg(test)]
#[path = "standard_output_tests.rs"]
mod tests;
