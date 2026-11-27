//! The stream against the observation tables of makeblastdb 2.17.0's C++ runtime
//! (`tests/fixtures/makeblastdb_ofstream/`): scripted operations on a `std::ofstream` whose file
//! failed past a byte limit, with the system calls made and the stream state after each
//! operation. Each table row is replayed on the stream over a file that logs its calls the same
//! way and fails the same way, and the log and the state must be the ones observed; at the end
//! of a sequence the file's bytes are those observed (size and checksum).
//!
//! Row format (see the fixture's README): `SEQ <n> limit=<L> mode=<M>` starts a sequence; a line
//! `<k> <op> | <calls> | st=<state>` is one operation and what it did; `DESTROY | <calls> |
//! FILE <size> <hash>` ends it: the calls of the destructor and the output file's size and
//! FNV-1a checksum. Operations: `O` open the output file, `OB` open a path that cannot be
//! opened, `OF` open `/dev/full` (every write fails, no call logged), `OM` open nothing,
//! `W n`/`I n` write/insert n bytes, `P` put a byte, `F` flush, `E` endl, `S pos` seekp, `C`
//! close, `K errno` make `close(2)` of the output file fail with `errno` from now on (0: stop),
//! `X` clear, `T bits` setstate. The bytes written come from one running counter per
//! sequence: byte k is `(k * 37 + (k >> 8) * 11 + 1) & 255`. A call is logged as `w <bytes> ->
//! <taken>` (`write`), `v ...` (`writev`), `l <pos> 0 -> <pos>` (`lseek`), `close`, `OPEN`;
//! a refused call ends `-1 ENOSPC`. Limit `L` (-1: none) is the number of bytes the output file
//! takes; mode 0: the call that reaches it takes up to it and the later calls fail; mode 1: a
//! call that would pass it fails and takes nothing.

use std::cell::RefCell;
use std::io;
use std::rc::Rc;

use super::{IoState, RawFile, Stream};

type Log = Rc<RefCell<Vec<String>>>;
type Disk = Rc<RefCell<Vec<u8>>>;

// The `errno` that `close(2)` fails with in the sequence being replayed (0: it succeeds).
// Set by the `K` operation; the tests of this module run one replay per thread.
thread_local! {
    static CLOSE_ERRNO: std::cell::Cell<i32> = const { std::cell::Cell::new(0) };
}

/// A file that follows the preloaded library's rules, logging its calls.
struct MockFile {
    limit: i64,
    mode: u8,
    total: i64,
    pos: usize,
    /// `None` for a file whose calls are not logged (`/dev/full`).
    log: Option<Log>,
    disk: Disk,
}

impl MockFile {
    fn event(&self, text: String) {
        if let Some(log) = &self.log {
            log.borrow_mut().push(text);
        }
    }
}

impl RawFile for MockFile {
    fn open(_path: &std::path::Path) -> io::Result<Self> {
        unreachable!("the tests open with `open_with`")
    }

    fn write(&mut self, segments: &[&[u8]]) -> io::Result<usize> {
        let count: usize = segments.iter().map(|s| s.len()).sum();
        let kind = if segments.len() == 1 { 'w' } else { 'v' };
        let take = if self.limit < 0 {
            Some(count)
        } else {
            let room = (self.limit - self.total) as usize;
            if room >= count {
                Some(count)
            } else if self.mode == 0 && room > 0 {
                Some(room)
            } else {
                None
            }
        };
        let Some(take) = take else {
            self.event(format!("{kind} {count} -> -1 ENOSPC"));
            return Err(io::Error::from_raw_os_error(28));
        };
        let mut disk = self.disk.borrow_mut();
        let mut left = take;
        for segment in segments {
            let part = &segment[..segment.len().min(left)];
            if self.pos > disk.len() {
                disk.resize(self.pos, 0);
            }
            let end = self.pos + part.len();
            if end > disk.len() {
                disk.resize(end, 0);
            }
            disk[self.pos..end].copy_from_slice(part);
            self.pos = end;
            left -= part.len();
        }
        self.total += take as i64;
        drop(disk);
        self.event(format!("{kind} {count} -> {take}"));
        Ok(take)
    }

    fn seek(&mut self, pos: u64) -> io::Result<()> {
        self.pos = pos as usize;
        self.event(format!("l {pos} 0 -> {pos}"));
        Ok(())
    }

    fn close(self) -> io::Result<()> {
        // `/dev/full` (no log) is not subject to the injected close error.
        let errno = if self.log.is_some() {
            CLOSE_ERRNO.with(|e| e.get())
        } else {
            0
        };
        if errno == 0 {
            self.event("close".to_string());
            Ok(())
        } else {
            self.event(format!("close -> -1 E{errno}"));
            Err(io::Error::from_raw_os_error(errno))
        }
    }
}

fn fnv(bytes: &[u8]) -> u64 {
    let mut h: u64 = 1469598103934665603;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(1099511628211);
    }
    h
}

fn pattern(counter: &mut u64, n: usize) -> Vec<u8> {
    (0..n)
        .map(|_| {
            let k = *counter;
            *counter += 1;
            ((k.wrapping_mul(37) + (k >> 8) * 11 + 1) & 255) as u8
        })
        .collect()
}

/// Replays the tables of `text`; returns the rows that differ (at most `max` of them).
pub(super) fn replay(text: &str, max: usize) -> (usize, Vec<String>) {
    let mut sequences = 0;
    let mut failures = Vec::new();
    let mut lines = text.lines().peekable();
    // The output file lives on between sequences: a sequence that opens nothing sees the
    // earlier one's bytes.
    let out_disk: Disk = Rc::new(RefCell::new(Vec::new()));
    while let Some(head) = lines.next() {
        let mut it = head.split_whitespace();
        assert_eq!(it.next(), Some("SEQ"), "row {head}");
        it.next();
        let limit: i64 = it
            .next()
            .unwrap()
            .trim_start_matches("limit=")
            .parse()
            .unwrap();
        let mode: u8 = it
            .next()
            .unwrap()
            .trim_start_matches("mode=")
            .parse()
            .unwrap();
        sequences += 1;
        let log: Log = Rc::new(RefCell::new(Vec::new()));
        let mut stream: Stream<MockFile> = Stream::new();
        let mut counter = 0u64;
        CLOSE_ERRNO.with(|e| e.set(0));
        let take_log = |log: &Log| log.borrow_mut().drain(..).collect::<Vec<_>>().join(" ; ");
        while let Some(&row) = lines.peek() {
            if row.starts_with("SEQ") {
                break;
            }
            lines.next();
            let mut cols = row.split(" | ");
            let op_text = cols.next().unwrap();
            let events = cols.next().unwrap().to_string();
            let tail = cols.next().unwrap().to_string();
            let got = if op_text == "DESTROY" {
                drop(std::mem::replace(&mut stream, Stream::new()));
                let disk = out_disk.borrow();
                format!(
                    "{} | FILE {} {:016x}",
                    take_log(&log),
                    disk.len(),
                    fnv(&disk)
                )
            } else {
                let mut words = op_text.split_whitespace();
                words.next();
                let op = words.next().unwrap();
                let arg = words.next().map(|w| w.parse::<u64>().unwrap());
                match (op, arg) {
                    ("O", _) => {
                        let (l, d) = (log.clone(), out_disk.clone());
                        stream.open_with(|| {
                            d.borrow_mut().clear();
                            l.borrow_mut().push("OPEN".to_string());
                            Ok(MockFile {
                                limit,
                                mode,
                                total: 0,
                                pos: 0,
                                log: Some(l),
                                disk: d,
                            })
                        });
                    }
                    ("OB", _) => {
                        stream.open_with(|| Err(io::Error::from_raw_os_error(21)));
                    }
                    ("OF", _) => {
                        stream.open_with(|| {
                            Ok(MockFile {
                                limit: 0,
                                mode: 1,
                                total: 0,
                                pos: 0,
                                log: None,
                                disk: Rc::new(RefCell::new(Vec::new())),
                            })
                        });
                    }
                    ("OM", _) => {}
                    ("W", Some(n)) => {
                        let data = pattern(&mut counter, n as usize);
                        stream.write(&data);
                    }
                    ("I", Some(n)) => {
                        let data = pattern(&mut counter, n as usize);
                        stream.insert(&data);
                    }
                    ("P", _) => {
                        let data = pattern(&mut counter, 1);
                        stream.put(data[0]);
                    }
                    ("F", _) => stream.flush(),
                    ("E", _) => stream.endl(),
                    ("S", Some(pos)) => stream.seekp(pos),
                    ("C", _) => stream.close(),
                    ("K", Some(errno)) => CLOSE_ERRNO.with(|e| e.set(errno as i32)),
                    ("X", _) => stream.clear(IoState::GOOD),
                    ("T", Some(bits)) => stream.setstate(IoState(bits as u8)),
                    other => panic!("unknown operation {other:?}"),
                }
                format!("{} | st={}", take_log(&log), stream.rdstate().0)
            };
            let want = format!("{events} | {tail}");
            if got != want && failures.len() < max {
                failures.push(format!("{head}\n  {op_text}\n  got  {got}\n  want {want}"));
            }
        }
    }
    (sequences, failures)
}

macro_rules! table_test {
    ($name:ident, $file:literal, $expect_seqs:expr) => {
        #[test]
        fn $name() {
            let text = include_str!(concat!("../../tests/fixtures/makeblastdb_ofstream/", $file));
            let (sequences, failures) = replay(text, 5);
            assert_eq!(sequences, $expect_seqs, "sequences in {}", $file);
            assert!(
                failures.is_empty(),
                "{} differ:\n{}",
                $file,
                failures.join("\n")
            );
        }
    };
}

table_test!(buffer_boundaries_table, "rows_boundaries.txt", 1274);
table_test!(failure_table, "rows_failures.txt", 108);
table_test!(state_table, "rows_states.txt", 416);
table_test!(random_table, "rows_random.txt", 800);
table_test!(close_table, "rows_close.txt", 807);

/// A file whose calls the test scripts: every write is taken in full, `lseek` and `close(2)`
/// fail as told, and the bytes taken are kept.
#[derive(Default)]
struct Scripted {
    seek_fails: bool,
    close_fails: bool,
    bytes: Rc<RefCell<Vec<u8>>>,
}

impl RawFile for Scripted {
    fn open(_path: &std::path::Path) -> io::Result<Self> {
        unreachable!("the tests open with `open_with`")
    }

    fn write(&mut self, segments: &[&[u8]]) -> io::Result<usize> {
        let mut n = 0;
        for segment in segments {
            self.bytes.borrow_mut().extend_from_slice(segment);
            n += segment.len();
        }
        Ok(n)
    }

    fn seek(&mut self, _pos: u64) -> io::Result<()> {
        if self.seek_fails {
            Err(io::Error::from_raw_os_error(29))
        } else {
            Ok(())
        }
    }

    fn close(self) -> io::Result<()> {
        if self.close_fails {
            Err(io::Error::from_raw_os_error(5))
        } else {
            Ok(())
        }
    }
}

fn scripted(seek_fails: bool, close_fails: bool) -> (Stream<Scripted>, Rc<RefCell<Vec<u8>>>) {
    let bytes = Rc::new(RefCell::new(Vec::new()));
    let mut stream = Stream::new();
    let b = bytes.clone();
    stream.open_with(move || {
        Ok(Scripted {
            seek_fails,
            close_fails,
            bytes: b,
        })
    });
    (stream, bytes)
}

// The standard: `seekp` sets `failbit` when the file buffer's seek fails (no observation of
// a failing `lseek` was made: the tables above have none).
#[test]
fn a_failed_lseek_sets_failbit_after_the_held_bytes_went_out() {
    let (mut stream, bytes) = scripted(true, false);
    stream.write(b"abc");
    stream.seekp(0);
    assert!(stream.fail() && !stream.bad());
    assert_eq!(*bytes.borrow(), b"abc");
}

// The standard: a failed `close()` of the file buffer sets `failbit`, not `badbit`; the
// held bytes are written first and the destructor of a stream that failed to close is silent.
#[test]
fn a_failed_close_sets_failbit_only_and_the_destructor_ignores_it() {
    let (mut stream, bytes) = scripted(false, true);
    stream.write(b"abc");
    stream.close();
    assert!(stream.fail() && !stream.bad());
    assert_eq!(*bytes.borrow(), b"abc");
    stream.clear(IoState::GOOD);
    stream.write(b"x");
    assert!(stream.bad(), "a stream with no file takes no bytes");
    let (mut stream, _) = scripted(false, true);
    stream.write(b"abc");
    drop(stream);
}

/// Bytes through a file on disk are those given, however the calls split them.
#[test]
fn the_disk_file_holds_the_bytes_written() {
    use super::{OfStream, BUFSIZ};
    let path = crate::utils::test_dir::temp_dir().join(format!(
        "ofstream-disk-{}",
        crate::utils::test_dir::process_id()
    ));
    let mut want = Vec::new();
    {
        let mut stream = OfStream::create(&path);
        let mut counter = 0u8;
        for n in [
            0usize,
            1,
            BUFSIZ - 2,
            BUFSIZ - 1,
            BUFSIZ,
            BUFSIZ + 1,
            3 * BUFSIZ + 5,
            7,
        ] {
            let data: Vec<u8> = (0..n)
                .map(|_| {
                    counter = counter.wrapping_add(37);
                    counter
                })
                .collect();
            stream.write(&data);
            want.extend_from_slice(&data);
            stream.put(counter);
            want.push(counter);
            stream.endl();
            want.push(b'\n');
        }
        stream.seekp(2);
        stream.write(b"XY");
        want[2] = b'X';
        want[3] = b'Y';
        assert!(stream.good());
        stream.close();
        assert!(stream.good());
    }
    assert_eq!(std::fs::read(&path).unwrap(), want);
    let _ = std::fs::remove_file(&path);
}

/// A device that takes no byte: the failing flush sets `badbit` and the later `close` sets
/// `failbit` (the held bytes cannot go out).
#[cfg(target_os = "linux")]
#[test]
fn the_disk_file_on_a_full_device_fails_at_the_flush() {
    let mut stream = super::OfStream::create(std::path::Path::new("/dev/full"));
    stream.write(b"abc");
    assert!(stream.good());
    stream.flush();
    assert!(stream.bad());
    stream.close();
    assert!(stream.fail());
}

#[test]
fn discard_unflushed_drops_held_bytes_without_writing_them() {
    let (mut stream, bytes) = scripted(false, false);
    stream.write(b"flushed\n");
    stream.flush();
    stream.write(b"buffered_only\n");
    stream.discard_unflushed();
    drop(stream);
    assert_eq!(*bytes.borrow(), b"flushed\n");
}
