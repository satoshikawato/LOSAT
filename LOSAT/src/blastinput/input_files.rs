//! The input files of LOSAT's command line as NCBI's arguments open them (`-query`,
//! `-subject`, `-out`), and NCBI's messages about the subjects that its formatter sets up.

use anyhow::Result;

/// Whether a FASTA file has no character but white space: NCBI reads it as a file without
/// records (an empty query, or no subject).
///
/// NCBI reference: c++/src/app/blast/blast_app_util.cpp:862-866
/// ```c
/// 	IOS_BASE::iostate orig_state = in.rdstate();
/// 	IOS_BASE::fmtflags orig_flags = in.setf(ios::skipws);
///
/// 	if(! (in >> c))
/// 		return true;
/// ```
/// `in >> c` skips the characters of C's `isspace`.
pub fn is_blank(bytes: &[u8]) -> bool {
    bytes
        .iter()
        .all(|byte| matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r'))
}

/// Opens a FASTA input of `program` as NCBI's argument does when a handler asks for its stream: `-` is
/// standard input, and a file that does not open gets NCBI's error
/// (`crate::cli::inaccessible`). A file named `.gz` is read through gzip (`open_file_input`).
///
/// NCBI reference: ncbi-blast/c++/src/corelib/ncbiargs.cpp:717-735
/// ```c
///     if (AsString() == "-") {
/// #if defined(NCBI_OS_MSWIN)
///         NcbiSys_setmode(NcbiSys_fileno(stdin), (mode & IOS_BASE::binary) ? O_BINARY : O_TEXT);
/// #endif
///         m_Ios  = &cin;
///     } else if ( !AsString().empty() ) {
///         if (!fstrm) {
///             fstrm = new CNcbiIfstream;
///         }
///         if (fstrm) {
///             fstrm->open(AsString().c_str(),IOS_BASE::in | mode);
///             if ( !fstrm->is_open() ) {
///                 delete fstrm;
///                 fstrm = NULL;
///             } else {
///                 m_DeleteFlag = true;
///             }
///         }
///         m_Ios = fstrm;
///     }
/// ```
pub fn open_input(path: &std::path::Path, role: &str, program: &str) -> Result<std::fs::File> {
    if path.as_os_str() == "-" {
        return standard_input().map_err(|_| {
            anyhow::anyhow!(
                "reading the {role} from standard input ('-') on this platform is not supported by LOSAT's {program}"
            )
        });
    }
    open_file_input(path, role)
}

/// Opens the file of a `-query` or `-subject` argument (`path` is not `-`): a file that does not
/// open gets NCBI's error (`crate::cli::inaccessible`), and a file named `.gz`, `.bz2` or `.zst`
/// (`compressed::compression_of_name`) is inflated into an unnamed temporary file first, which is
/// returned in its place: the program then reads it exactly as it reads the plain file (the same
/// lines and messages, the same readable-byte hint of `FastaStream::from_file`, the same second
/// handle of blastn's query stream and second reads of BLASTX's query), while the decoder holds
/// its window. Concatenated members, streams and frames (bgzip, pbzip2, pzstd, `cat a.gz b.gz`)
/// are read whole. A file that is not what its name says, or is cut short, is refused here,
/// before the search (line GZ, DECISIONS 3, 5 and 9).
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
                    "-{role}: '{}' is named {suffix} but is not {format} data ({error}); LOSAT reads a file named {suffix} through {format}",
                    path.display(),
                    suffix = compression.suffix(),
                    format = compression.format(),
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
            process_id()
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

/// The process id in the temporary file's name: `std::process::id()`.
#[cfg(not(target_os = "wasi"))]
fn process_id() -> u32 {
    std::process::id()
}

/// 0 on WASI, which has no process id (Rust's `std::process::id()` panics there,
/// "unsupported", `library/std/src/sys/pal/wasip1/os.rs:131` of 1.92); `create_new` and the
/// time stamp keep the names apart.
#[cfg(target_os = "wasi")]
fn process_id() -> u32 {
    0
}

/// Standard input as a file that shares its position, as `cin` does.
fn standard_input() -> std::io::Result<std::fs::File> {
    #[cfg(any(unix, target_os = "wasi"))]
    {
        use std::os::fd::AsFd;
        std::io::stdin()
            .as_fd()
            .try_clone_to_owned()
            .map(std::fs::File::from)
    }
    #[cfg(windows)]
    {
        use std::os::windows::io::AsHandle;
        std::io::stdin()
            .as_handle()
            .try_clone_to_owned()
            .map(std::fs::File::from)
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
/// NCBI reference: c++/src/app/blast/blast_app_util.cpp:856-859
/// ```c
/// 	char c;
/// 	CNcbiStreampos orig_p = in.tellg();
/// 	// Piped input
/// 	if(orig_p < 0)
/// 		return false;
/// ```
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

/// Rejects a file name that is not UTF-8, where NCBI opens the file (the subjects, then the
/// queries, then the output: `CBlastDatabaseArgs` comes before `CStdCmdLineArgs`).
///
/// NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.cpp:903-911
/// ```c
/// GetSubjectFile(const CArgs& args)
/// {
/// 	string filename="";
///
/// 	if (args.Exist(kArgSubject) && args[kArgSubject].HasValue())
/// 		filename = args[kArgSubject].AsString();
///
/// 	return filename;
/// }
/// ```
/// NCBI takes a file name as bytes: it writes the `-subject` name into the outfmt 0 and 7
/// reports (`Database: User specified sequence set (Input: ...)`) and every name into its
/// error messages as they are, which LOSAT's UTF-8 strings do not reproduce (plan DW-13).
pub fn check_utf8_file_name(path: &std::path::Path, role: &str, program: &str) -> Result<()> {
    if path.to_str().is_none() {
        anyhow::bail!(
            "the -{role} file name {:?} is not UTF-8; NCBI BLAST+ writes the bytes of file names as they are, which is not supported by LOSAT's {program}",
            path.to_string_lossy()
        );
    }
    Ok(())
}

/// NCBI's warning for a subject without residues, which it skips when it sets up the
/// subjects for the search (`CBlastFormat` is made after `Query is Empty!` and before the
/// outfmt 0 prolog; `AUTHORITY.md` §B row 7, BI-43, BI-49). The ID of a subject is
/// `Subject_<n>` and its title is the record's title (`FastaRecord::title`), written
/// as bytes; a subject without a title gives `Subject_<n> : ...`.
///
/// NCBI reference: c++/src/algo/blast/api/blast_setup_cxx.cpp:773-788
/// ```c
///         catch(CBlastException & e ) {
///         	// Skip bad subject sequence
///         	if(e.GetErrCode() == CBlastException::eInvalidArgument) {
///         		seqblk_vec->push_back(subj);
///         		string warning = kEmptyStr;
///         		const CSeq_id *  id = subjects.GetSeqId(i);
///         		string title = subjects.GetTitle(i);
///         		if(id != NULL) {
///         			warning = id->GetSeqIdString() + " ";
///         		}
///         		warning += subjects.GetTitle(i);
///         		if(warning != kEmptyStr){
///         			warning += ": ";
///         		}
///         		warning += "Subject sequence contains no data";
///         		ERR_POST(Warning << warning);
/// ```
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
    for (index, record) in records.iter().enumerate() {
        if length(index) == 0 {
            let mut warning = match record.believed() {
                Some(ids) => {
                    let mut warning = format!("Warning: [{program}] ").into_bytes();
                    warning.extend_from_slice(&ids.warning_id);
                    warning.push(b' ');
                    warning
                }
                None => format!("Warning: [{program}] Subject_{} ", index + 1).into_bytes(),
            };
            warning.extend_from_slice(&record.title);
            warning.extend_from_slice(b": Subject sequence contains no data\n");
            diagnostics.write_all(&warning)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{has_stream_position, open_file_input};
    use crate::blastinput::compressed::test_data::{bzip2_streams, gzip_members, zstd_frames};
    use crate::blastinput::compressed::Compression;
    use std::io::Read;

    #[cfg(unix)]
    fn owned(reader: std::io::PipeReader) -> std::os::fd::OwnedFd {
        reader.into()
    }

    #[cfg(windows)]
    fn owned(reader: std::io::PipeReader) -> std::os::windows::io::OwnedHandle {
        reader.into()
    }

    fn compressed(compression: Compression, parts: &[&[u8]]) -> Vec<u8> {
        match compression {
            Compression::Gzip => gzip_members(parts),
            Compression::Bzip2 => bzip2_streams(parts),
            Compression::Zstd => zstd_frames(parts),
        }
    }

    fn temporary(stem: &str, bytes: &[u8]) -> std::path::PathBuf {
        // The stem ends the name: the suffix is what `open_file_input` reads.
        let path = std::env::temp_dir().join(format!(
            "losat-input-files-{}-{}-{stem}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&path, bytes).unwrap();
        path
    }

    fn read_all(mut file: std::fs::File) -> Vec<u8> {
        let mut bytes = Vec::new();
        file.read_to_end(&mut bytes).unwrap();
        bytes
    }

    /// Names left in `directory` by this process's inflated inputs (`unnamed_temporary_file`):
    /// none on Unix, where the name goes as soon as the file is open.
    #[cfg(unix)]
    fn inflated_names(directory: &std::path::Path) -> Vec<std::path::PathBuf> {
        let prefix = format!(".losat-query-{}-", std::process::id());
        std::fs::read_dir(directory)
            .unwrap()
            .filter_map(|entry| entry.ok().map(|entry| entry.path()))
            .filter(|path| {
                path.file_name()
                    .and_then(|name| name.to_str())
                    .is_some_and(|name| name.starts_with(&prefix))
            })
            .collect()
    }

    // No NCBI counterpart (line GZ, DECISIONS 3, 4 and 9): a compressed input opens as the
    // inflated bytes, every member / stream / frame of a concatenated file, as a regular file
    // with a position, and leaves no name behind; a plain name opens as it is.
    #[test]
    fn a_compressed_input_opens_as_the_inflated_file_with_every_part() {
        let plain = temporary("plain.fa", b">q1\nACGT\n>q2\nTTTT\n");
        assert_eq!(
            read_all(open_file_input(&plain, "query").unwrap()),
            b">q1\nACGT\n>q2\nTTTT\n"
        );
        std::fs::remove_file(&plain).unwrap();
        for (compression, stem) in [
            (Compression::Gzip, "members.fa.GZ"),
            (Compression::Bzip2, "streams.fa.Bz2"),
            (Compression::Zstd, "frames.fa.zst"),
        ] {
            let path = temporary(
                stem,
                &compressed(compression, &[b">q1\nACGT\n", b">q2\nTTTT\n", b""]),
            );
            let file = open_file_input(&path, "query").unwrap();
            assert!(has_stream_position(&file), "{stem}");
            assert!(file.metadata().unwrap().is_file(), "{stem}");
            assert_eq!(read_all(file), b">q1\nACGT\n>q2\nTTTT\n", "{stem}");
            std::fs::remove_file(&path).unwrap();
        }
        #[cfg(unix)]
        assert_eq!(
            inflated_names(&std::env::temp_dir()),
            Vec::<std::path::PathBuf>::new()
        );
    }

    /// The directory of an open, unlinked file: `/proc/self/fd/N` reads
    /// `<directory>/<name> (deleted)`.
    #[cfg(target_os = "linux")]
    fn directory_of_open_file(file: &std::fs::File) -> std::path::PathBuf {
        use std::os::fd::AsRawFd;
        std::fs::read_link(format!("/proc/self/fd/{}", file.as_raw_fd()))
            .unwrap()
            .parent()
            .unwrap()
            .to_path_buf()
    }

    // No NCBI counterpart (DECISIONS 12): the inflated copy is made in the input's directory,
    // where the data lives (`temp_dir()` may be memory), and under `temp_dir()` when that
    // directory refuses it; either way no name is left behind.
    #[cfg(unix)]
    #[test]
    fn the_inflated_copy_is_made_beside_the_input_or_else_under_temp_dir() {
        use std::os::unix::fs::PermissionsExt;
        let directory = std::env::temp_dir().join(format!(
            "losat-input-files-{}-{}-beside",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir(&directory).unwrap();
        let path = directory.join("q.fa.gz");
        std::fs::write(&path, gzip_members(&[b">q1\nACGT\n"])).unwrap();
        let file = open_file_input(&path, "query").unwrap();
        #[cfg(target_os = "linux")]
        assert_eq!(
            directory_of_open_file(&file),
            std::fs::canonicalize(&directory).unwrap()
        );
        assert_eq!(inflated_names(&directory), Vec::<std::path::PathBuf>::new());
        assert_eq!(read_all(file), b">q1\nACGT\n");

        // The directory refuses the copy (read-only; not as root, which writes anyway): the
        // copy is made under `temp_dir()`.
        std::fs::set_permissions(&directory, std::fs::Permissions::from_mode(0o555)).unwrap();
        let probe = directory.join(".probe");
        if std::fs::File::create(&probe).is_ok() {
            std::fs::remove_file(&probe).unwrap();
        } else {
            let file = open_file_input(&path, "query").unwrap();
            #[cfg(target_os = "linux")]
            assert_eq!(
                directory_of_open_file(&file),
                std::fs::canonicalize(std::env::temp_dir()).unwrap()
            );
            assert_eq!(read_all(file), b">q1\nACGT\n");
        }
        std::fs::set_permissions(&directory, std::fs::Permissions::from_mode(0o755)).unwrap();
        std::fs::remove_dir_all(&directory).unwrap();
    }

    #[test]
    fn the_inflated_copy_goes_beside_the_input_first_and_then_under_temp_dir() {
        let temp = std::env::temp_dir();
        let nested = std::path::Path::new("data").join("q.fa.gz");
        let directories =
            |path: &std::path::Path| super::inflation_directories(path).collect::<Vec<_>>();
        assert_eq!(
            directories(&nested),
            vec![std::path::PathBuf::from("data"), temp.clone()]
        );
        assert_eq!(
            directories(std::path::Path::new("q.fa.gz")),
            vec![std::path::PathBuf::from("."), temp.clone()]
        );
        assert_eq!(directories(&temp.join("q.fa.gz")), vec![temp]);
    }

    // No NCBI counterpart (line GZ, DECISIONS 5): a compressed name whose bytes are not that
    // format, or are cut short, is refused with the file's name and the format; a file that does
    // not open keeps NCBI's error.
    #[test]
    fn a_compressed_name_without_that_data_is_refused() {
        for compression in Compression::ALL {
            let stem = format!("text.fa{}", compression.suffix());
            let plain = temporary(&stem, b">q1\nACGT\n");
            let error = open_file_input(&plain, "query").unwrap_err().to_string();
            assert!(error.starts_with("-query: '"), "{error}");
            let expected = format!(
                "' is named {suffix} but is not {format} data (",
                suffix = compression.suffix(),
                format = compression.format()
            );
            assert!(error.contains(&expected), "{error}");
            assert!(
                error.ends_with(&format!(
                    "); LOSAT reads a file named {} through {}",
                    compression.suffix(),
                    compression.format()
                )),
                "{error}"
            );
            std::fs::remove_file(&plain).unwrap();
            let whole = compressed(
                compression,
                &[b">q1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"],
            );
            let cut = temporary(
                &format!("cut.fa{}", compression.suffix()),
                &whole[..whole.len() - 12],
            );
            let error = open_file_input(&cut, "subject").unwrap_err().to_string();
            assert!(error.starts_with("-subject: '"), "{error}");
            std::fs::remove_file(&cut).unwrap();
        }
        let missing = std::env::temp_dir().join("losat-input-files-missing.fa.gz");
        let error = open_file_input(&missing, "query").unwrap_err().to_string();
        assert!(
            error.starts_with(
                "Command line argument error: Argument \"query\". File is not accessible:"
            ),
            "{error}"
        );
    }

    // NCBI reference: c++/src/app/blast/blast_app_util.cpp:856-859 (`tellg() < 0` is piped
    // input). A regular file has a position; a pipe has none, on Linux because `lseek`
    // fails and on Windows because `GetFileType` is not `FILE_TYPE_DISK`.
    #[test]
    fn a_regular_file_has_a_position_and_a_pipe_has_none() {
        let path = std::env::temp_dir().join(format!(
            "losat-input-files-has-stream-position-{}",
            std::process::id()
        ));
        std::fs::write(&path, b">q\nACGT\n").unwrap();
        let file = std::fs::File::open(&path).unwrap();
        assert!(has_stream_position(&file));
        drop(file);
        std::fs::remove_file(&path).unwrap();
        // WASI has no pipe to make here.
        #[cfg(any(unix, windows))]
        {
            let (reader, mut writer) = std::io::pipe().unwrap();
            std::io::Write::write_all(&mut writer, b">q\nACGT\n").unwrap();
            drop(writer);
            let pipe = std::fs::File::from(owned(reader));
            assert!(!has_stream_position(&pipe));
        }
        // `NUL` is a character device that is not a console: it has a position, as
        // `/dev/null` has on Linux (`Stdio::null()` with `-query -` is an empty file).
        #[cfg(windows)]
        assert!(has_stream_position(&std::fs::File::open("NUL").unwrap()));
        #[cfg(unix)]
        assert!(has_stream_position(
            &std::fs::File::open("/dev/null").unwrap()
        ));
    }
}
