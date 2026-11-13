//! `OfStream` against libstdc++ 13.2's `ofstream`: what reaches the file after each call.

use super::{IoState, OfStream, BUFSIZ};
use std::path::{Path, PathBuf};

/// A directory of its own under the system's temporary directory.
struct TempDir(PathBuf);

impl TempDir {
    fn new(name: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat-ofstream-{name}-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn on_disk(path: &Path) -> u64 {
    std::fs::metadata(path).unwrap().len()
}

/// A stream whose file cannot be opened (a directory) is failed from the start: nothing is
/// written, every output does nothing and the state only gains `failbit`.
#[test]
fn a_file_that_cannot_be_opened_fails_every_output() {
    let tmp = TempDir::new("open");
    let mut stream = OfStream::create(&tmp.0);
    assert_eq!(stream.rdstate(), IoState::FAIL);
    stream.write(b"abc");
    stream.put(b'x');
    stream.flush();
    stream.seekp(8);
    stream.endl();
    assert_eq!(stream.rdstate(), IoState::FAIL);
    assert!(stream.fail() && !stream.bad());
    stream.close();
    assert_eq!(stream.rdstate(), IoState::FAIL);
    assert!(tmp.0.is_dir());
    assert_eq!(std::fs::read_dir(&tmp.0).unwrap().count(), 0);
}

/// With its state cleared, a stream without a file takes no byte: `write` sets `badbit`.
#[test]
fn a_write_without_a_file_sets_badbit() {
    let mut stream = OfStream::new();
    assert!(stream.good());
    stream.write(b"abc");
    assert_eq!(stream.rdstate(), IoState::BAD);
    // The next output's sentry adds `failbit`.
    stream.write(b"d");
    assert_eq!(stream.rdstate(), IoState::BAD | IoState::FAIL);
    stream.clear(IoState::GOOD);
    stream.put(b'x');
    assert_eq!(stream.rdstate(), IoState::BAD);
}

/// The put area holds `BUFSIZ - 1` bytes. A write smaller than the room left is buffered; one
/// that is not goes out at once with what the buffer holds, and leaves the put area empty.
#[test]
fn writes_go_through_a_put_area_of_8191_bytes() {
    let tmp = TempDir::new("area");
    let path = tmp.0.join("f");
    let mut stream = OfStream::create(&path);
    assert_eq!(on_disk(&path), 0);
    // Uncommitted: the room is the whole put area.
    stream.write(&vec![1u8; BUFSIZ - 2]);
    assert_eq!((on_disk(&path), stream.pending()), (0, BUFSIZ - 2));
    // One byte left: a write of one byte is not smaller than the room, so it goes out with the
    // buffer.
    stream.write(&[2]);
    assert_eq!((on_disk(&path), stream.pending()), (BUFSIZ as u64 - 1, 0));
    // 4096 + 4095 bytes fit; the next 4095 do not (4095 left).
    stream.write(&[3u8; 4096]);
    stream.write(&[4u8; 4094]);
    assert_eq!((on_disk(&path), stream.pending()), (8191, 8190));
    stream.write(&[5u8; 1]);
    assert_eq!((on_disk(&path), stream.pending()), (8191 + 8191, 0));
    // A first write of the whole put area goes out directly.
    let other = tmp.0.join("g");
    let mut direct = OfStream::create(&other);
    direct.write(&vec![6u8; BUFSIZ - 1]);
    assert_eq!((on_disk(&other), direct.pending()), (8191, 0));
    direct.write(&vec![7u8; BUFSIZ - 2]);
    assert_eq!((on_disk(&other), direct.pending()), (8191, 8190));
    drop(direct);
    assert_eq!(on_disk(&other), 8191 + 8190);
    stream.close();
    assert!(stream.good());
    let bytes = std::fs::read(&path).unwrap();
    assert_eq!(bytes.len(), 8191 + 8191);
    assert_eq!(&bytes[8190..8192], &[2, 3]);
    assert_eq!(bytes[bytes.len() - 1], 5);
}

/// `put` into a full put area appends its byte in the buffer's last byte and writes all 8192.
#[test]
fn a_put_into_a_full_area_writes_8192_bytes() {
    let tmp = TempDir::new("put");
    let path = tmp.0.join("f");
    let mut stream = OfStream::create(&path);
    for _ in 0..BUFSIZ - 1 {
        stream.put(b'a');
    }
    assert_eq!((on_disk(&path), stream.pending()), (0, BUFSIZ - 1));
    stream.put(b'b');
    assert_eq!((on_disk(&path), stream.pending()), (BUFSIZ as u64, 0));
    stream.endl();
    assert_eq!(on_disk(&path), BUFSIZ as u64 + 1);
}

/// `flush`, `seekp` and the close write what is buffered; `seekp` moves the file position.
#[test]
fn flush_seekp_and_close_write_the_buffer() {
    let tmp = TempDir::new("seek");
    let path = tmp.0.join("f");
    let mut stream = OfStream::create(&path);
    stream.write(&[0u8; 24]);
    stream.flush();
    assert_eq!((on_disk(&path), stream.pending()), (24, 0));
    stream.write(b"tail");
    stream.seekp(8);
    assert_eq!(on_disk(&path), 28);
    stream.write(b"12345678");
    assert_eq!(stream.pending(), 8);
    stream.close();
    assert!(stream.good());
    let mut expected = vec![0u8; 24];
    expected[8..16].copy_from_slice(b"12345678");
    expected.extend_from_slice(b"tail");
    assert_eq!(std::fs::read(&path).unwrap(), expected);
    // A closed stream fails its next close.
    stream.close();
    assert_eq!(stream.rdstate(), IoState::FAIL);
}

/// A stream dropped while it holds bytes writes them (`~basic_filebuf` closes).
#[test]
fn the_destructor_writes_the_buffer() {
    let tmp = TempDir::new("drop");
    let path = tmp.0.join("f");
    let mut stream = OfStream::create(&path);
    stream.insert(b"#\n# Alias file created: ");
    assert_eq!(on_disk(&path), 0);
    drop(stream);
    assert_eq!(std::fs::read(&path).unwrap(), b"#\n# Alias file created: ");
}

/// A device that takes no byte (`/dev/full`): a buffered write succeeds and the flush fails;
/// the bytes stay buffered and every later output does nothing.
#[cfg(target_os = "linux")]
#[test]
fn a_write_that_fails_sets_badbit_and_later_output_does_nothing() {
    let mut stream = OfStream::create(Path::new("/dev/full"));
    assert!(stream.good());
    stream.write(&[1u8; 479]);
    assert!(stream.good());
    assert_eq!(stream.pending(), 479);
    stream.flush();
    assert_eq!(stream.rdstate(), IoState::BAD);
    assert_eq!(stream.pending(), 479);
    stream.write(&[2u8; 10]);
    assert_eq!(stream.rdstate(), IoState::BAD | IoState::FAIL);
    assert_eq!(stream.pending(), 479);
    // With its state cleared, a write that does not fit fails at once and keeps the buffer.
    stream.clear(IoState::GOOD);
    stream.write(&[3u8; BUFSIZ - 479]);
    assert_eq!(stream.rdstate(), IoState::BAD);
    assert_eq!(stream.pending(), 479);
    // A first write of 4096 bytes is buffered; the second of 4096 does not fit the 4095 left.
    let mut full = OfStream::create(Path::new("/dev/full"));
    full.write(&[0u8; 4096]);
    assert!(full.good());
    full.write(&[0u8; 4096]);
    assert_eq!(full.rdstate(), IoState::BAD);
}

/// `CWriteDB_File::Write` advances `m_Offset` by what it was given, whether the stream took it or
/// not (writedb_files.cpp:218-221): the offsets of a volume whose file could not be created are
/// those of the bytes it would hold.
#[test]
fn a_volume_file_counts_its_offset_when_its_stream_failed() {
    let tmp = TempDir::new("offset");
    let base = tmp.0.join("db").to_string_lossy().into_owned();
    std::fs::create_dir(format!("{base}.00.phr")).unwrap();
    let mut file = crate::writedb::files::WriteDbFile::lazy(&base, "phr", 0, 0);
    file.create_now();
    assert!(file.created());
    assert_eq!(file.write(b"abc").unwrap(), 3);
    assert_eq!(file.write_int4(7).unwrap(), 7);
    file.close();
    assert!(Path::new(&format!("{base}.00.phr")).is_dir());
}
