//! The input files of LOSAT's command line as NCBI's arguments open them (`-query`,
//! `-subject`, `-out`), and NCBI's messages about the subjects that its formatter sets up.

use anyhow::Result;

/// Whether a FASTA file has no character but white space: NCBI reads it as a file without
/// records (an empty query, or no subject).
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:862-866
/// `in >> c` skips the characters of C's `isspace`: the test is `NStr::IsBlank`'s
/// (`ncbi_str::is_blank`); the library API's `read_query_records` reads it too.
pub fn is_blank(bytes: &[u8]) -> bool {
    crate::ncbi_str::is_blank(bytes)
}

/// Opens a FASTA input of `program` as NCBI's argument does when a handler asks for its stream: `-` is
/// standard input, and a file that does not open gets NCBI's error
/// (`crate::cli::inaccessible`). A file named `.gz` is read through gzip (`open_file_input`).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:717-735
/// `-` is the standard input itself, which this opening and every other `-` read in place and
/// leave open (`InputFile`, as `m_DeleteFlag` stays false for `cin`).
pub fn open_input(path: &std::path::Path, role: &str, program: &str) -> Result<InputFile> {
    if path.as_os_str() == "-" {
        return standard_input().map_err(|_| {
            anyhow::anyhow!(
                "reading the {role} from standard input ('-') on this platform is not supported by LOSAT's {program}"
            )
        });
    }
    open_file_input(path, role).map(InputFile::from)
}

/// An opened input as the FASTA readers read it (`open_input`): a file of its own, closed when
/// this is dropped, or the standard input, which every `-` reads in place and which stays open.
/// It reads, seeks and is a `std::fs::File` (`Deref`) as either.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:609-613
/// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:715-721
/// Every `-` argument is the one `cin`, never closed: its position is the process's standard
/// input's (`-query -` with `-subject -` reads what the subjects left, nothing). The standard
/// input is not copied (`dup` on Unix, `DuplicateHandle` on Windows), which the WASI targets
/// cannot do: on every target it is the process's own descriptor or handle, viewed as a file.
pub struct InputFile(InputHandle);

/// The handle of an `InputFile`.
enum InputHandle {
    /// A file that the input opened (or the inflated copy of a compressed one).
    Own(std::fs::File),
    /// The process's standard input (`standard_input`), never closed by LOSAT.
    StandardInput(std::mem::ManuallyDrop<std::fs::File>),
    /// A file read as the standard input is read (`InputFile::as_standard_input`), closed when
    /// dropped: a test's, or the copy of the standard input that a run reads twice
    /// (`InputFile::into_seekable_copy`).
    CopiedStandardInput(std::fs::File),
}

impl InputFile {
    /// `file` read in place and never closed: the standard input (`standard_input`).
    fn in_place(file: std::fs::File) -> Self {
        let file = std::mem::ManuallyDrop::new(file);
        Self(InputHandle::StandardInput(file))
    }

    /// A `file` read as `-` reads the standard input (`is_standard_input`), but closed when
    /// dropped.
    pub(crate) fn as_standard_input(file: std::fs::File) -> Self {
        Self(InputHandle::CopiedStandardInput(file))
    }

    /// The bytes of this input from its position to its end, copied into an unnamed temporary
    /// file under `temp_dir()` and read from there from the start, as this input is read (the
    /// standard input stays `cin`'s kind of stream, `is_standard_input`; a path stays a file).
    /// For a run that reads its query twice, which a pipe or a stream cannot be read for.
    /// `None` where there is no temporary directory (WASI) or the copy fails; the read
    /// consumes the input as the run's own read would.
    pub(crate) fn into_seekable_copy(mut self) -> Option<Self> {
        use std::io::{Seek, Write};
        let directory = temporary_directory()?;
        let mut copy = unnamed_temporary_file(&directory, "query").ok()?;
        let mut sink = std::io::BufWriter::with_capacity(1 << 16, &mut copy);
        let mut buffer = vec![0u8; 1 << 16];
        loop {
            match std::io::Read::read(&mut self, &mut buffer) {
                Ok(0) => break,
                Ok(read) => sink.write_all(&buffer[..read]).ok()?,
                Err(error) if error.kind() == std::io::ErrorKind::Interrupted => {}
                Err(_) => return None,
            }
        }
        sink.flush().ok()?;
        drop(sink);
        copy.seek(std::io::SeekFrom::Start(0)).ok()?;
        Some(if self.is_standard_input() {
            Self::as_standard_input(copy)
        } else {
            Self::from(copy)
        })
    }

    /// The file to read, seek and ask about.
    fn file(&self) -> &std::fs::File {
        match &self.0 {
            InputHandle::Own(file) => file,
            InputHandle::StandardInput(file) => file,
            InputHandle::CopiedStandardInput(file) => file,
        }
    }

    /// Whether this is the standard input that `-` names, which NCBI reads through `cin`
    /// (`m_Ios` points to `cin`, cited at `InputFile`), not through a `CNcbiIfstream` of a path (a path
    /// such as `/dev/stdin` is a file of its own). The FASTA readers read the two as NCBI's two
    /// streams answer a refill (`fasta_reader::stream::InputBackend::of`).
    pub(crate) fn is_standard_input(&self) -> bool {
        match self.0 {
            InputHandle::Own(_) => false,
            InputHandle::StandardInput(_) => true,
            InputHandle::CopiedStandardInput(_) => true,
        }
    }
}

impl From<std::fs::File> for InputFile {
    fn from(file: std::fs::File) -> Self {
        Self(InputHandle::Own(file))
    }
}

impl std::ops::Deref for InputFile {
    type Target = std::fs::File;

    fn deref(&self) -> &std::fs::File {
        self.file()
    }
}

// Reads and seeks go to the file, as `&std::fs::File` reads and seeks (one position per open
// file, shared by every handle of it).
impl std::io::Read for InputFile {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.file().read(buf)
    }
    fn read_vectored(&mut self, bufs: &mut [std::io::IoSliceMut<'_>]) -> std::io::Result<usize> {
        self.file().read_vectored(bufs)
    }
    fn read_to_end(&mut self, buf: &mut Vec<u8>) -> std::io::Result<usize> {
        self.file().read_to_end(buf)
    }
    fn read_to_string(&mut self, buf: &mut String) -> std::io::Result<usize> {
        self.file().read_to_string(buf)
    }
}

impl std::io::Read for &InputFile {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.file().read(buf)
    }
    fn read_vectored(&mut self, bufs: &mut [std::io::IoSliceMut<'_>]) -> std::io::Result<usize> {
        self.file().read_vectored(bufs)
    }
    fn read_to_end(&mut self, buf: &mut Vec<u8>) -> std::io::Result<usize> {
        self.file().read_to_end(buf)
    }
    fn read_to_string(&mut self, buf: &mut String) -> std::io::Result<usize> {
        self.file().read_to_string(buf)
    }
}

impl std::io::Seek for InputFile {
    fn seek(&mut self, pos: std::io::SeekFrom) -> std::io::Result<u64> {
        self.file().seek(pos)
    }
    fn stream_position(&mut self) -> std::io::Result<u64> {
        self.file().stream_position()
    }
}

impl std::io::Seek for &InputFile {
    fn seek(&mut self, pos: std::io::SeekFrom) -> std::io::Result<u64> {
        self.file().seek(pos)
    }
    fn stream_position(&mut self) -> std::io::Result<u64> {
        self.file().stream_position()
    }
}

// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
/// The bytes of an opened `-query` or `-subject` input (`open_input`, `open_file_input`) from its
/// position to its end, as the FASTA readers read it (`fasta_reader::stream`): a read that fails
/// ends the input as its end does, with the bytes read before it (NCBI's line reader is at its end
/// once the stream has failed). A directory opens and fails its first read: no byte, so the
/// command line's query of a directory is empty (`Query is Empty!`) and its subject set has no
/// record. The library's `api::read_input_file` reads its files here.
pub(crate) fn read_to_failure(input: &mut impl std::io::Read) -> Vec<u8> {
    let mut bytes = Vec::new();
    // `read_to_end` keeps the bytes before a failed read; the failure itself is the end.
    let _ = input.read_to_end(&mut bytes);
    bytes
}

/// Opens the file of a `-query` or `-subject` argument (`path` is not `-`): a file that does not
/// open gets NCBI's error (`crate::cli::inaccessible`), and a file named `.gz`, `.bz2` or `.zst`
/// (`compressed::compression_of_name`) is inflated into an unnamed temporary file first, which is
/// returned in its place: the program then reads it exactly as it reads the plain file (the same
/// lines and messages, the same readable-byte hint of `FastaStream::from_file`, the same second
/// handle of blastn's query stream and second reads of BLASTX's query), while the decoder holds
/// its window. Concatenated members, streams and frames (bgzip, pbzip2, pzstd, `cat a.gz b.gz`)
/// are read whole, and a `.gz` file as makeblastdb reads it (`compressed::GzipMembers`: a member's
/// trailer is not checked). A file that is not what its name says, or is cut short (for `.gz`:
/// where makeblastdb's decompressor fails), is refused here, before the search (line GZ,
/// DECISIONS 3, 5 and 9; for `.gz` narrowed by swebok-improve SW8 DECISIONS 38 and 41).
///
/// No NCBI counterpart: the NCBI search programs read such a file as text and stop at its first
/// line (`CFastaReader: Near line 1, there's a line that doesn't look like plausible data`);
/// NCBI's makeblastdb reads `-in` names with these suffixes through a decompressing stream
/// (`writedb::makeblastdb`), and LOSAT applies the same rule, by the name alone, to the search
/// inputs (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`, exception 8).
pub fn open_file_input(path: &std::path::Path, role: &str) -> Result<std::fs::File> {
    let file = std::fs::File::open(path).map_err(|_| crate::cli::inaccessible(role, path))?;
    let Some(compression) = crate::blastinput::compressed::compression_of_name(path) else {
        return Ok(file);
    };
    // The inflated copy goes beside the input, whose directory holds the data (`temp_dir()` is
    // memory where `/tmp` is a tmpfs); when that directory refuses the file or has no room for
    // it, the copy goes under `temp_dir()` (DECISIONS 12).
    let mut source = Some(file);
    let mut refusals = Vec::new();
    for directory in inflation_directories(path) {
        let input = match source.take() {
            Some(file) => file,
            None => std::fs::File::open(path).map_err(|_| crate::cli::inaccessible(role, path))?,
        };
        match inflate_to_temporary_file(compression.decoder(input), &directory, role) {
            Ok(inflated) => return Ok(inflated),
            Err(InflateFault::Data(error)) => {
                return Err(anyhow::anyhow!(
                    "-{role}: {}",
                    compression.not_data(&path.display().to_string(), &error.to_string())
                ));
            }
            Err(InflateFault::Copy(error)) => {
                refusals.push(format!("{error} in '{}'", directory.display()));
            }
        }
    }
    Err(anyhow::anyhow!(
        "-{role}: no room for the inflated copy of '{}' ({})",
        path.display(),
        refusals.join("; ")
    ))
}

/// Where the inflated copy of the input `path` may go: the input's own directory first, then
/// the temporary directory (`temporary_directory`) when that is another directory. The second
/// is looked up only when the caller asks for it, after the first refused the copy.
fn inflation_directories(path: &std::path::Path) -> impl Iterator<Item = std::path::PathBuf> {
    let beside = match path.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => parent.to_path_buf(),
        _ => std::path::PathBuf::from("."),
    };
    let fallback = {
        let beside = beside.clone();
        std::iter::once_with(move || temporary_directory().filter(|fallback| *fallback != beside))
            .flatten()
    };
    std::iter::once(beside).chain(fallback)
}

/// The directory of temporary files: `temp_dir()`.
#[cfg(not(target_os = "wasi"))]
fn temporary_directory() -> Option<std::path::PathBuf> {
    #[allow(clippy::disallowed_methods)] // not WASI: `temp_dir()` is only called off WASI
    Some(std::env::temp_dir())
}

/// None on WASI: Rust's `temp_dir()` panics on the WASI targets ("no filesystem on wasm",
/// `library/std/src/sys/pal/wasip1/os.rs:119` of 1.92), where a host gives a module directories
/// only by its preopens. The WASI command puts the inflated copy beside the input or refuses it.
#[cfg(target_os = "wasi")]
fn temporary_directory() -> Option<std::path::PathBuf> {
    None
}

/// What stopped an inflation: the compressed data (`Data`: not that format, cut short, empty)
/// or the copy (`Copy`: the directory refused the temporary file or ran out of room).
enum InflateFault {
    Data(std::io::Error),
    Copy(std::io::Error),
}

/// Copies the bytes of `inflater` (a decoder of `compressed::Compression`) into an unnamed
/// temporary file in `directory` (`unnamed_temporary_file`) and returns it, positioned at its
/// start.
fn inflate_to_temporary_file(
    mut inflater: impl std::io::Read,
    directory: &std::path::Path,
    role: &str,
) -> std::result::Result<std::fs::File, InflateFault> {
    use std::io::{Seek, Write};
    let mut temporary = unnamed_temporary_file(directory, role).map_err(InflateFault::Copy)?;
    {
        let mut sink = std::io::BufWriter::with_capacity(1 << 16, &mut temporary);
        let mut buffer = vec![0u8; 1 << 16];
        loop {
            let read = match inflater.read(&mut buffer) {
                Ok(0) => break,
                Ok(read) => read,
                Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(error) => return Err(InflateFault::Data(error)),
            };
            sink.write_all(&buffer[..read])
                .map_err(InflateFault::Copy)?;
        }
        sink.flush().map_err(InflateFault::Copy)?;
    }
    temporary
        .seek(std::io::SeekFrom::Start(0))
        .map_err(InflateFault::Copy)?;
    Ok(temporary)
}

/// A new, empty file in `directory` that no name reaches: on Unix and WASI its name is removed
/// as soon as it is open, on Windows it is hidden and deleted when its last handle closes
/// (`FILE_FLAG_DELETE_ON_CLOSE`, shared for reading, writing and deletion so that `try_clone`
/// and a second `open` of the handle keep working). It is readable and writable.
fn unnamed_temporary_file(
    directory: &std::path::Path,
    role: &str,
) -> std::io::Result<std::fs::File> {
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|since| since.as_nanos())
        .unwrap_or(0);
    for attempt in 0u32..64 {
        let path = directory.join(format!(
            ".losat-{role}-{}-{stamp}-{attempt}.inflated",
            crate::utils::process_id::process_id()
        ));
        let mut options = std::fs::OpenOptions::new();
        options.read(true).write(true).create_new(true);
        #[cfg(windows)]
        {
            use std::os::windows::fs::OpenOptionsExt;
            use windows_sys::Win32::Storage::FileSystem::{
                FILE_ATTRIBUTE_HIDDEN, FILE_ATTRIBUTE_TEMPORARY, FILE_FLAG_DELETE_ON_CLOSE,
                FILE_SHARE_DELETE, FILE_SHARE_READ, FILE_SHARE_WRITE,
            };
            options
                .share_mode(FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE)
                .custom_flags(
                    FILE_ATTRIBUTE_HIDDEN | FILE_ATTRIBUTE_TEMPORARY | FILE_FLAG_DELETE_ON_CLOSE,
                );
        }
        match options.open(&path) {
            Ok(file) => {
                // The name goes at once; the open handle keeps the file. On WASI a host that
                // cannot remove an open file leaves the name until the process ends.
                #[cfg(not(windows))]
                let _ = std::fs::remove_file(&path);
                return Ok(file);
            }
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(error),
        }
    }
    Err(std::io::Error::new(
        std::io::ErrorKind::AlreadyExists,
        format!("no free temporary file name in {}", directory.display()),
    ))
}

/// The standard input as an `InputFile` (`open_input`'s `-`, NCBI's `cin`): the process's own
/// descriptor 0 (handle `STD_INPUT_HANDLE` on Windows) viewed as a file and never closed, so every
/// `-` shares its position. The Rust standard library reads its `Stdin` the same way (descriptor
/// 0 in a `ManuallyDrop`, `library/std/src/sys/stdio/unix.rs` and `wasip1.rs` of 1.92). Before,
/// the descriptor was copied (`try_clone_to_owned`, `dup`), which the WASI targets refuse.
fn standard_input() -> std::io::Result<InputFile> {
    #[cfg(any(unix, target_os = "wasi"))]
    {
        use std::os::fd::{AsRawFd, FromRawFd};
        // SAFETY: descriptor 0 is open for the whole process (on Unix the Rust runtime opens
        // `/dev/null` there before `main` when it was closed; on WASI it is the host's standard
        // input) and LOSAT never closes it. The `File` is held in a `ManuallyDrop` that is never
        // dropped (`InputFile::in_place`), so it only borrows the descriptor: it is never closed
        // through this `File`.
        let file = unsafe { std::fs::File::from_raw_fd(std::io::stdin().as_raw_fd()) };
        Ok(InputFile::in_place(file))
    }
    #[cfg(windows)]
    {
        use std::os::windows::io::{AsRawHandle, FromRawHandle};
        // A process started without a standard input handle (Rust gives a null handle then)
        // has no standard input to read: the copy of the handle failed there before.
        let handle = std::io::stdin().as_raw_handle();
        if handle.is_null() {
            return Err(std::io::ErrorKind::NotFound.into());
        }
        // SAFETY: the standard input handle is open for the whole process and LOSAT never
        // closes it. The `File` is held in a `ManuallyDrop` that is never dropped
        // (`InputFile::in_place`), so it only borrows the handle: it is never closed through
        // this `File`.
        let file = unsafe { std::fs::File::from_raw_handle(handle) };
        Ok(InputFile::in_place(file))
    }
    #[cfg(not(any(unix, windows, target_os = "wasi")))]
    {
        Err(std::io::ErrorKind::Unsupported.into())
    }
}

/// Whether an opened input has a stream position the way a regular file on Linux has one:
/// `tellg()` of NCBI's `IsIStreamEmpty` is not negative (a pipe has none, a terminal has none,
/// while a character device such as `/dev/null` has one).
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-860
/// On Unix and WASI the answer is whether the position can be read, as before
/// (`lseek` fails on a pipe and on a terminal, not on `/dev/null`; the position is not moved).
/// On Windows `SetFilePointerEx` succeeds on a pipe or a console handle (its result is
/// undefined there), so an empty pipe would look like an empty file. Follow Linux by the type
/// of the handle (`GetFileType`): a disk file has a position; a character device has one
/// unless it is a console (`GetConsoleMode` succeeds on a console only; `NUL` has a position,
/// as `/dev/null` has); a pipe or an unknown type has none. LOSAT takes the Linux branch of
/// `IsIStreamEmpty` on every platform (`AUTHORITY.md` §C1).
pub fn has_stream_position(file: &std::fs::File) -> bool {
    #[cfg(windows)]
    {
        use std::os::windows::io::AsRawHandle;
        use windows_sys::Win32::Foundation::HANDLE;
        use windows_sys::Win32::Storage::FileSystem::{
            GetFileType, FILE_TYPE_CHAR, FILE_TYPE_DISK,
        };
        use windows_sys::Win32::System::Console::GetConsoleMode;
        let handle: HANDLE = file.as_raw_handle().cast();
        // SAFETY: `GetFileType` only reads the type of the handle it is given; the handle comes
        // from `file`, which is open for the call. A failure returns FILE_TYPE_UNKNOWN (0),
        // which has no position.
        match unsafe { GetFileType(handle) } {
            FILE_TYPE_DISK => true,
            FILE_TYPE_CHAR => {
                let mut mode = 0u32;
                // SAFETY: `handle` is open as above and `mode` is a valid `u32` to write; the
                // call fails (returns 0) on a handle that is not a console.
                unsafe { GetConsoleMode(handle, &mut mode) == 0 }
            }
            _ => false,
        }
    }
    #[cfg(not(windows))]
    {
        std::io::Seek::stream_position(&mut &*file).is_ok()
    }
}

/// NCBI's warning for a subject without residues, which it skips when it sets up the
/// subjects for the search (`CBlastFormat` is made after `Query is Empty!` and before the
/// outfmt 0 prolog; `AUTHORITY.md` §B row 7, BI-43, BI-49). The ID of a subject is
/// `Subject_<n>` and its title is the record's title (`FastaRecord::title`), written
/// as bytes; a subject without a title gives `Subject_<n> : ...`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:773-788
pub fn write_empty_subject_warnings(
    records: &[crate::blastinput::fasta_reader::FastaRecord],
    program: &str,
    diagnostics: &mut dyn std::io::Write,
) -> std::io::Result<()> {
    write_empty_subject_warnings_by_length(
        records,
        |index| records[index].seq().len(),
        program,
        diagnostics,
    )
}

/// `write_empty_subject_warnings` for records whose letters are held elsewhere
/// (`sequence::packed_subject`): `length(index)` is the number of letters of record `index`.
pub fn write_empty_subject_warnings_by_length(
    records: &[crate::blastinput::fasta_reader::FastaRecord],
    length: impl Fn(usize) -> usize,
    program: &str,
    diagnostics: &mut dyn std::io::Write,
) -> std::io::Result<()> {
    write_empty_subject_warnings_of(
        records
            .iter()
            .enumerate()
            .map(|(index, record)| (record.names(), length(index))),
        program,
        diagnostics,
    )
}

/// `write_empty_subject_warnings` of the subjects whose names (`fasta_reader::RecordNames`) and
/// numbers of letters `records` gives, in their order: the shared reader's records and BLASTX's
/// (`algorithm::blastx::input::FastaRecord::names`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:773-788 (cited at
/// `write_empty_subject_warnings`)
pub fn write_empty_subject_warnings_of<'a>(
    records: impl IntoIterator<Item = (crate::blastinput::fasta_reader::RecordNames<'a>, usize)>,
    program: &str,
    diagnostics: &mut dyn std::io::Write,
) -> std::io::Result<()> {
    for (index, (record, length)) in records.into_iter().enumerate() {
        if length == 0 {
            let mut text = match record.believed {
                Some(ids) => {
                    let mut text = ids.warning_id.clone();
                    text.push(b' ');
                    text
                }
                None => format!("Subject_{} ", index + 1).into_bytes(),
            };
            text.extend_from_slice(record.title);
            text.extend_from_slice(b": Subject sequence contains no data");
            // A post of NCBI's main thread (`ncbi_exception::post_line`; tblastn's formatter, which
            // makes it, comes after the PSSM query's annotations are indexed).
            diagnostics.write_all(&crate::ncbi_exception::post_line(
                crate::ncbi_exception::Severity::Warning,
                None,
                Some(program),
                &text,
            ))?;
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "input_files_tests.rs"]
mod tests;
