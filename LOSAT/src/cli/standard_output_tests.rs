use super::*;
use std::io::Write;

/// A file of this test, made empty (tests of one binary run in parallel).
fn scratch_file(name: &str) -> std::path::PathBuf {
    let path = crate::utils::test_dir::temp_dir().join(format!(
        "losat-standard-output-{name}-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::write(&path, b"").unwrap();
    path
}

#[cfg(any(unix, windows, target_os = "wasi"))]
#[test]
fn each_write_reaches_the_device_at_once_and_nothing_is_held() {
    let path = scratch_file("at-once");
    let file = std::fs::OpenOptions::new()
        .append(true)
        .open(&path)
        .unwrap();
    let mut out = StandardOutput(Device::Raw(std::mem::ManuallyDrop::new(file)));
    // A write with no newline: `Stdout`'s line buffer would hold it until a newline or a flush.
    assert_eq!(out.write(b"no newline").unwrap(), 10);
    assert_eq!(std::fs::read(&path).unwrap(), b"no newline");
    out.write_all(b" and the rest\n").unwrap();
    assert_eq!(std::fs::read(&path).unwrap(), b"no newline and the rest\n");
    // The test owns this file: close it.
    // Windows has other devices; on Unix and WASI the pattern always matches.
    #[allow(irrefutable_let_patterns)]
    if let StandardOutput(Device::Raw(file)) = out {
        drop(std::mem::ManuallyDrop::into_inner(file));
    }
    std::fs::remove_file(&path).unwrap();
}

#[cfg(unix)]
#[test]
fn a_descriptor_that_is_not_open_counts_as_written_and_other_errors_do_not() {
    let closed = std::io::Error::from_raw_os_error(rustix::io::Errno::BADF.raw_os_error());
    assert!(is_closed_descriptor(&closed));
    for other in [
        rustix::io::Errno::AGAIN,
        rustix::io::Errno::NOSPC,
        rustix::io::Errno::PIPE,
    ] {
        assert!(!is_closed_descriptor(&std::io::Error::from_raw_os_error(
            other.raw_os_error()
        )));
    }
    assert!(!is_closed_descriptor(
        &std::io::ErrorKind::WouldBlock.into()
    ));
}

#[cfg(windows)]
#[test]
fn a_handle_that_is_not_open_counts_as_written_and_other_errors_do_not() {
    use windows_sys::Win32::Foundation::{ERROR_DISK_FULL, ERROR_INVALID_HANDLE, ERROR_NO_DATA};
    assert!(is_closed_descriptor(&std::io::Error::from_raw_os_error(
        ERROR_INVALID_HANDLE as i32
    )));
    for other in [ERROR_DISK_FULL, ERROR_NO_DATA] {
        assert!(!is_closed_descriptor(&std::io::Error::from_raw_os_error(
            other as i32
        )));
    }
}
