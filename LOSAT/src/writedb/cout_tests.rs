use std::io::{self, Write};

use super::cout;
use crate::cli::NcbiStdio;
use crate::writedb::ofstream::StreamBuf;

/// A device that takes `limit` bytes in all: the write that reaches the limit is short and the
/// writes after it fail (a file under `RLIMIT_FSIZE` with SIGXFSZ ignored, a device that fills
/// up); `calls` are the sizes of the writes that reached it.
struct Limited {
    limit: usize,
    taken: Vec<u8>,
    calls: Vec<usize>,
    error: io::ErrorKind,
}

impl Limited {
    fn new(limit: usize) -> Self {
        Limited {
            limit,
            taken: Vec::new(),
            calls: Vec::new(),
            error: io::ErrorKind::FileTooLarge,
        }
    }
}

impl Write for Limited {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.calls.push(buf.len());
        let room = self.limit - self.taken.len();
        if room == 0 && !buf.is_empty() {
            return Err(io::Error::from(self.error));
        }
        let n = buf.len().min(room);
        self.taken.extend_from_slice(&buf[..n]);
        Ok(n)
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

const FILE: NcbiStdio = NcbiStdio::Cout {
    block: 4096,
    lines: false,
};

#[test]
fn bytes_reach_the_device_where_standard_output_s_buffer_writes_them() {
    let mut device = Limited::new(usize::MAX);
    let mut stream = cout(&mut device, FILE);
    // The first insertion is held but for its whole blocks; a later one that fits is held, one
    // that does not fit writes the buffer full and the whole blocks after it.
    stream.write(&[b'a'; 4096]);
    stream.write(&[b'b'; 100]);
    stream.write(&[b'c'; 4000]);
    stream.write(&[b'd'; 4096]);
    stream.flush();
    assert!(!stream.fail());
    drop(stream);
    assert_eq!(device.calls, [4096, 4096, 4096, 4]);
    assert_eq!(device.taken.len(), 4096 + 100 + 4000 + 4096);
}

#[test]
fn a_failed_write_sets_badbit_and_leaves_the_buffer_empty() {
    let mut device = Limited::new(10);
    let mut stream = cout(&mut device, FILE);
    stream.write(b"0123456789abc");
    assert!(!stream.fail(), "held");
    stream.flush();
    assert!(stream.fail(), "the flush failed");
    // The stream cleared: nothing is held, so a flush writes nothing and succeeds.
    stream.clear(super::super::ofstream::IoState::GOOD);
    stream.flush();
    assert!(!stream.fail());
    // A failed stream writes nothing (its sentry).
    stream.setstate(super::super::ofstream::IoState::BAD);
    stream.write(b"xyz");
    stream.flush();
    drop(stream);
    assert_eq!(device.taken, b"0123456789");
}

#[test]
fn a_closed_pipe_is_told_apart() {
    let mut device = Limited::new(0);
    device.error = io::ErrorKind::BrokenPipe;
    let mut stream = cout(&mut device, FILE);
    stream.put(b'\n');
    stream.flush();
    assert!(stream.fail());
    assert_eq!(stream.rdbuf().closed_pipe(), cfg!(not(windows)));
}
