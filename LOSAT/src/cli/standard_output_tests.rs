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

/// A gathered write is the device's own (`writev`): both parts land in one call, as the report
/// stream writes NCBI's held bytes and the bytes after them.
#[cfg(any(unix, target_os = "wasi"))]
#[test]
fn a_gathered_write_is_one_write_of_the_device() {
    let path = scratch_file("gathered");
    let file = std::fs::OpenOptions::new()
        .append(true)
        .open(&path)
        .unwrap();
    let mut out = StandardOutput(Device::Raw(std::mem::ManuallyDrop::new(file)));
    let parts = [
        std::io::IoSlice::new(b"held "),
        std::io::IoSlice::new(b"and new\n"),
    ];
    assert_eq!(out.write_vectored(&parts).unwrap(), 13);
    assert_eq!(std::fs::read(&path).unwrap(), b"held and new\n");
    let StandardOutput(Device::Raw(file)) = out;
    drop(std::mem::ManuallyDrop::into_inner(file));
    std::fs::remove_file(&path).unwrap();
}

/// A write's error is the device's: on a descriptor opened only for reading (`1<file`) the
/// write fails, as NCBI's does (it reports the failure). Rust's `Stdout` counts it as done.
#[cfg(unix)]
#[test]
fn a_write_to_a_descriptor_opened_for_reading_fails() {
    let path = scratch_file("read-only");
    let file = std::fs::File::open(&path).unwrap();
    let mut out = StandardOutput(Device::Raw(std::mem::ManuallyDrop::new(file)));
    let error = out.write(b"text").unwrap_err();
    assert_eq!(
        error.raw_os_error(),
        Some(rustix::io::Errno::BADF.raw_os_error())
    );
    #[allow(irrefutable_let_patterns)]
    if let StandardOutput(Device::Raw(file)) = out {
        drop(std::mem::ManuallyDrop::into_inner(file));
    }
    std::fs::remove_file(&path).unwrap();
}
