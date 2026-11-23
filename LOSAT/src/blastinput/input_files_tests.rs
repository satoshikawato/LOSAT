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
    let path = crate::utils::test_dir::temp_dir().join(format!(
        "losat-input-files-{}-{}-{stem}",
        crate::utils::test_dir::process_id(),
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
    let prefix = format!(".losat-query-{}-", crate::utils::test_dir::process_id());
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
// with a position, and leaves no name behind; a plain name opens as it is. The inputs lie in a
// directory of this test's own, where their inflated copies are made (DECISIONS 12), so the check
// for names left behind reads only this test's: other tests of the process open compressed
// inputs under `temp_dir()` at the same time, whose names exist there for a moment (the race that
// line PS's gate of PS15b-1 met, 2026-11-23).
#[test]
fn a_compressed_input_opens_as_the_inflated_file_with_every_part() {
    let directory = crate::utils::test_dir::temp_dir().join(format!(
        "losat-input-files-{}-{}-parts",
        crate::utils::test_dir::process_id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir(&directory).unwrap();
    let plain = directory.join("plain.fa");
    std::fs::write(&plain, b">q1\nACGT\n>q2\nTTTT\n").unwrap();
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
        let path = directory.join(stem);
        std::fs::write(
            &path,
            compressed(compression, &[b">q1\nACGT\n", b">q2\nTTTT\n", b""]),
        )
        .unwrap();
        let file = open_file_input(&path, "query").unwrap();
        assert!(has_stream_position(&file), "{stem}");
        assert!(file.metadata().unwrap().is_file(), "{stem}");
        assert_eq!(read_all(file), b">q1\nACGT\n>q2\nTTTT\n", "{stem}");
        std::fs::remove_file(&path).unwrap();
    }
    #[cfg(unix)]
    assert_eq!(inflated_names(&directory), Vec::<std::path::PathBuf>::new());
    std::fs::remove_dir(&directory).unwrap();
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
    let directory = crate::utils::test_dir::temp_dir().join(format!(
        "losat-input-files-{}-{}-beside",
        crate::utils::test_dir::process_id(),
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
            std::fs::canonicalize(crate::utils::test_dir::temp_dir()).unwrap()
        );
        assert_eq!(read_all(file), b">q1\nACGT\n");
    }
    std::fs::set_permissions(&directory, std::fs::Permissions::from_mode(0o755)).unwrap();
    std::fs::remove_dir_all(&directory).unwrap();
}

#[test]
fn the_inflated_copy_goes_beside_the_input_first_and_then_under_temp_dir() {
    let temp = crate::utils::test_dir::temp_dir();
    let nested = std::path::Path::new("data").join("q.fa.gz");
    let directories =
        |path: &std::path::Path| super::inflation_directories(path).collect::<Vec<_>>();
    // WASI has no directory of temporary files (`temporary_directory`): beside the input only.
    let with_temp = |beside: &str| {
        let mut all = vec![std::path::PathBuf::from(beside)];
        if cfg!(not(target_os = "wasi")) {
            all.push(temp.clone());
        }
        all
    };
    assert_eq!(directories(&nested), with_temp("data"));
    assert_eq!(directories(std::path::Path::new("q.fa.gz")), with_temp("."));
    assert_eq!(directories(&temp.join("q.fa.gz")), vec![temp]);
}

// No NCBI counterpart (line GZ, DECISIONS 5): a compressed name whose bytes are not that
// format, or are cut short (for `.gz`: where makeblastdb's decompressor fails), is refused
// with the file's name and the format; a file that does not open keeps NCBI's error.
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
    let missing = crate::utils::test_dir::temp_dir().join("losat-input-files-missing.fa.gz");
    let error = open_file_input(&missing, "query").unwrap_err().to_string();
    assert!(
        error.starts_with(
            "Command line argument error: Argument \"query\". File is not accessible:"
        ),
        "{error}"
    );
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:715-721
// Every `-` is the process's descriptor 0 itself (not a copy, which the WASI targets cannot
// make), and dropping the input leaves it open (`m_DeleteFlag` false: `cin` is never deleted).
#[cfg(unix)]
#[test]
fn standard_input_is_descriptor_0_itself_and_stays_open() {
    use std::os::fd::{AsFd, AsRawFd};
    for _ in 0..2 {
        let input = super::open_input(std::path::Path::new("-"), "query", "BLASTN").unwrap();
        assert_eq!(input.as_raw_fd(), 0);
        drop(input);
        assert!(std::io::stdin().as_fd().try_clone_to_owned().is_ok());
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:609-613
// Two inputs read in place (`InputFile::in_place`, as two `-` arguments read `cin`) share one
// position, and dropping one leaves the descriptor open for the other: here a pipe, whose
// bytes the first takes in part and the second reads on from there.
#[cfg(unix)]
#[test]
fn inputs_read_in_place_share_the_position_and_leave_the_descriptor_open() {
    use std::os::fd::{FromRawFd, IntoRawFd, OwnedFd};
    let (reader, mut writer) = std::io::pipe().unwrap();
    std::io::Write::write_all(&mut writer, b"abcdef").unwrap();
    drop(writer);
    let raw = owned(reader).into_raw_fd();
    // SAFETY: `raw` is the open read end of the pipe, which this test owns; the two inputs
    // only borrow it (`in_place` never closes), and it is closed once, by the `OwnedFd`
    // made of it at the end.
    let view = || super::InputFile::in_place(unsafe { std::fs::File::from_raw_fd(raw) });
    let (mut first, mut second) = (view(), view());
    let mut head = [0u8; 3];
    first.read_exact(&mut head).unwrap();
    assert_eq!(&head, b"abc");
    drop(first);
    let mut rest = Vec::new();
    second.read_to_end(&mut rest).unwrap();
    assert_eq!(rest, b"def");
    drop(second);
    // SAFETY: `raw` is still open (neither input closed it) and nothing else owns it.
    drop(unsafe { OwnedFd::from_raw_fd(raw) });
}

// No NCBI counterpart (the copy only makes a stream a file that a run can read twice): the
// bytes of a pipe from the input's position are read from the start of the copy, which has a
// position and is read as `cin` when the input was the standard input.
#[test]
fn a_stream_copied_for_a_second_read_has_a_position_and_keeps_its_kind() {
    use std::io::{Read, Seek, SeekFrom};
    for standard in [false, true] {
        let path = crate::utils::test_dir::temp_dir().join(format!(
            "losat-input-files-seekable-copy-{standard}-{}",
            crate::utils::test_dir::process_id()
        ));
        std::fs::write(&path, b">q\nACGT\n").unwrap();
        let mut file = std::fs::File::open(&path).unwrap();
        file.seek(SeekFrom::Start(3)).unwrap();
        let input = if standard {
            super::InputFile::as_standard_input(file)
        } else {
            super::InputFile::from(file)
        };
        let mut copy = input.into_seekable_copy().unwrap();
        assert_eq!(copy.is_standard_input(), standard);
        assert!(super::has_stream_position(&copy));
        for _ in 0..2 {
            copy.seek(SeekFrom::Start(0)).unwrap();
            let mut bytes = Vec::new();
            copy.read_to_end(&mut bytes).unwrap();
            assert_eq!(bytes, b"ACGT\n");
        }
        std::fs::remove_file(&path).unwrap();
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-859 (a negative `tellg()` is
// piped input). A regular file has a position; a pipe has none, on Linux because `lseek`
// fails and on Windows because `GetFileType` is not `FILE_TYPE_DISK`.
#[test]
fn a_regular_file_has_a_position_and_a_pipe_has_none() {
    let path = crate::utils::test_dir::temp_dir().join(format!(
        "losat-input-files-has-stream-position-{}",
        crate::utils::test_dir::process_id()
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
