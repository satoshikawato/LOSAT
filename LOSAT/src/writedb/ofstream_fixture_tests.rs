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
//! close, `X` clear, `T bits` setstate. The bytes written come from one running counter per
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

    fn close(self) {
        self.event("close".to_string());
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
        let mut bad = false;
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
                    ("X", _) => stream.clear(IoState::GOOD),
                    ("T", Some(bits)) => stream.setstate(IoState(bits as u8)),
                    other => panic!("unknown operation {other:?}"),
                }
                format!("{} | st={}", take_log(&log), stream.rdstate().0)
            };
            let want = format!("{events} | {tail}");
            if got != want {
                bad = true;
                if failures.len() < max {
                    failures.push(format!("{head}\n  {op_text}\n  got  {got}\n  want {want}"));
                }
            }
        }
        let _ = bad;
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
