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
/// (`crate::cli::inaccessible`).
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
    std::fs::File::open(path).map_err(|_| crate::cli::inaccessible(role, path))
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
/// `tellg()` of NCBI's `IsIStreamEmpty` is not negative (a pipe has none, a console has none).
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
/// (`lseek` fails on a pipe; the position is not moved). On Windows `SetFilePointerEx` succeeds on a pipe or a console
/// handle (its result is undefined there), so an empty pipe would look like an empty file:
/// only a handle of `FILE_TYPE_DISK` has a position. LOSAT takes the Linux branch of
/// `IsIStreamEmpty` on every platform (`AUTHORITY.md` §C1).
pub fn has_stream_position(file: &std::fs::File) -> bool {
    #[cfg(windows)]
    {
        use std::os::windows::io::AsRawHandle;
        const FILE_TYPE_DISK: u32 = 1;
        #[link(name = "kernel32")]
        extern "system" {
            fn GetFileType(handle: *mut std::ffi::c_void) -> u32;
        }
        // SAFETY: `GetFileType` only reads the type of the handle it is given; the handle comes
        // from `file`, which is open for the call. A failure returns FILE_TYPE_UNKNOWN (0),
        // which is not a disk: no position.
        unsafe { GetFileType(file.as_raw_handle().cast()) == FILE_TYPE_DISK }
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
            let mut warning = format!("Warning: [{program}] Subject_{} ", index + 1).into_bytes();
            warning.extend_from_slice(&record.title);
            warning.extend_from_slice(b": Subject sequence contains no data\n");
            diagnostics.write_all(&warning)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::has_stream_position;

    #[cfg(unix)]
    fn owned(reader: std::io::PipeReader) -> std::os::fd::OwnedFd {
        reader.into()
    }

    #[cfg(windows)]
    fn owned(reader: std::io::PipeReader) -> std::os::windows::io::OwnedHandle {
        reader.into()
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
    }
}
