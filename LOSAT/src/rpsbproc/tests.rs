use super::*;

// ncbifile.cpp:820-971 (cited at `ncbi_environment::normalize_path`).
#[test]
fn normalize_path_as_ncbi() {
    for (path, normalized) in [
        ("", ""),
        ("/", "/"),
        ("///", "/"),
        ("losat", "losat"),
        ("./losat", "losat"),
        ("bin/./losat", "bin/losat"),
        ("bin//losat", "bin/losat"),
        ("bin/../losat", "losat"),
        ("../bin/losat", "../bin/losat"),
        ("/../bin/losat", "/bin/losat"),
        ("bin/losat/", "bin/losat"),
        (".", "."),
        ("./", "."),
    ] {
        assert_eq!(
            normalize_path(path.as_bytes()),
            normalized.as_bytes(),
            "{path}"
        );
    }
}

// rpsbproc_main.cpp:1809,1836 (cited at `print_dom_line`, `print_cluster_line`): `m_iFrom` plus
// `COORDSBASE` for a hit that ends at 2^31 - 1 overflows `int` (both 0.5.1 builds print
// -2147483648: PS16c review oracle `plusmax`); its value, 2^31, is what NCBI prints for the hit
// moved 3 bases down plus 3 (PS16-Q1 oracle `ntplusmax` and its nearby input `ntplusmax3`;
// PD-LOSAT-NCBI-DEFECTS exception 12). The XML input's value one less than `INT_MIN` prints
// `INT_MIN`.
#[test]
fn printed_positions_past_int() {
    assert_eq!(printed_position(i32::MAX), 2_147_483_648);
    assert_eq!(printed_position(i32::MAX - 3), 2_147_483_645);
    assert_eq!(printed_position(i64::from(i32::MIN) - 1), -2_147_483_648);
    assert_eq!(printed_position(0), 1);
}

// Streaming a `double` of -0.0 to an `ostream` is glibc's `%g`, "-0".
#[test]
fn format_g_keeps_the_sign_of_zero() {
    assert_eq!(format_g(-0.0), "-0");
    assert_eq!(format_g(0.0), "0");
    assert_eq!(format_g(3.06789e-55), "3.06789e-55");
    assert_eq!(format_g(179.393), "179.393");
    assert_eq!(OVERLAPLEADING, 1_410_065_408);
}

/// A writer whose bytes, and the size of each of its writes (a gathered write is one, as
/// `writev`), the test reads while the report is written.
#[derive(Clone, Default)]
struct Shared(std::rc::Rc<std::cell::RefCell<(Vec<u8>, Vec<usize>)>>);

impl Write for Shared {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.write_vectored(&[std::io::IoSlice::new(buf)])
    }

    fn write_vectored(&mut self, bufs: &[std::io::IoSlice<'_>]) -> std::io::Result<usize> {
        let mut shared = self.0.borrow_mut();
        let (bytes, writes) = &mut *shared;
        let before = bytes.len();
        for buf in bufs {
            bytes.extend_from_slice(buf);
        }
        writes.push(bytes.len() - before);
        Ok(bytes.len() - before)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// The report's writes before line RG-KB 9a-2: every insertion copied into `pending`, which is
/// written once NCBI's stream has written a block of it.
struct CopyingReport {
    pending: Vec<u8>,
    buffer: crate::cli::NcbiBuffer,
    file: Vec<u8>,
    /// The size of each write (one `write_all` of a block).
    writes: Vec<usize>,
}

impl CopyingReport {
    fn write_out(&mut self, n: usize) {
        if n > 0 {
            self.file.extend_from_slice(&self.pending[..n]);
            self.writes.push(n);
        }
        self.pending.drain(..n);
    }

    fn write(&mut self, bytes: &[u8]) {
        self.pending.extend_from_slice(bytes);
        self.buffer.sputn(bytes.len(), &self.pending);
        if bytes.last() == Some(&b'\n') {
            self.buffer.flush();
        }
        if self.pending.len() - self.buffer.held() >= REPORT_BLOCK {
            self.write_out(self.pending.len() - self.buffer.held());
        }
    }

    fn put(&mut self, byte: u8) {
        self.pending.push(byte);
        self.buffer.sputc(&self.pending);
        if self.pending.len() - self.buffer.held() >= REPORT_BLOCK {
            self.write_out(self.pending.len() - self.buffer.held());
        }
    }
}

// No NCBI counterpart: how LOSAT gets NCBI's writes of the report (`cli::NcbiBuffer`, observed
// from NCBI's rpsbproc) to the file. The insertions that NCBI's stream has written go to the file
// from where they are (`Report::write_through`): the file, its writes (one per block, of the same
// size) and the held bytes after every insertion are those of the report that copied each
// insertion first, for every stream, and a string of a megabyte is not kept in the report's
// buffer.
#[test]
fn report_writes_long_insertions_where_they_are() {
    let streams = [
        crate::cli::NcbiStdio::File,
        crate::cli::NcbiStdio::Cout {
            block: 4096,
            lines: false,
        },
        crate::cli::NcbiStdio::Cout {
            block: 4096,
            lines: true,
        },
        crate::cli::NcbiStdio::Cout {
            block: 8192,
            lines: false,
        },
        crate::cli::NcbiStdio::Cout {
            block: 64,
            lines: true,
        },
    ];
    let mut seed = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = move |bound: u64| {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed % bound
    };
    for stream in streams {
        let file = Shared::default();
        let diag = Diagnostics::new(crate::blastinput::ncbi_environment::DiagTarget::Stderr);
        let mut report = Report::new(Box::new(file.clone()), stream, true, diag);
        let mut old = CopyingReport {
            pending: Vec::new(),
            buffer: crate::cli::NcbiBuffer::new(stream),
            file: Vec::new(),
            writes: Vec::new(),
        };
        for step in 0..2000u64 {
            let size = match next(20) {
                0 => 0,
                1 | 2 => REPORT_BLOCK as u64 - 1 + next(3),
                3..=5 => next(3 * REPORT_BLOCK as u64),
                6 => 1 << 20,
                _ => next(300),
            } as usize;
            let mut bytes: Vec<u8> = (0..size)
                .map(|i| b'a' + ((step as usize + i) % 26) as u8)
                .collect();
            for _ in 0..next(3) {
                if size > 0 {
                    let at = next(size as u64) as usize;
                    bytes[at] = b'\n';
                }
            }
            match next(8) {
                0 => {
                    report.put(b'\t');
                    old.put(b'\t');
                }
                1 => {
                    report.sync();
                    let n = old.pending.len() - old.buffer.held();
                    old.write_out(n);
                }
                _ => {
                    report.write(&bytes);
                    old.write(&bytes);
                }
            }
            // (The bytes and writes of each step, compared and dropped.)
            assert_eq!(file.0.borrow().0, old.file, "{stream:?} step {step}");
            assert_eq!(file.0.borrow().1, old.writes, "{stream:?} step {step}");
            file.0.borrow_mut().0.clear();
            file.0.borrow_mut().1.clear();
            old.file.clear();
            old.writes.clear();
            assert_eq!(report.pending, old.pending, "{stream:?} step {step}");
            assert_eq!(report.buffer, old.buffer, "{stream:?} step {step}");
            assert!(
                report.pending.capacity() <= 2 * REPORT_BLOCK + crate::cli::NCBI_BUFSIZ,
                "{stream:?} step {step}: {}",
                report.pending.capacity()
            );
        }
        report.flush();
        let n = old.pending.len();
        old.write_out(n);
        assert_eq!(file.0.borrow().0, old.file, "{stream:?}");
        assert_eq!(file.0.borrow().1, old.writes, "{stream:?}");
    }
}

/// A device that takes at most `limit` bytes of a write (of all the parts of a gathered one), is
/// interrupted at its first write, and takes no byte once it is `full`.
struct Device {
    taken: Vec<u8>,
    writes: usize,
    limit: usize,
    full: bool,
}

impl Write for Device {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.write_vectored(&[std::io::IoSlice::new(buf)])
    }

    fn write_vectored(&mut self, bufs: &[std::io::IoSlice<'_>]) -> std::io::Result<usize> {
        self.writes += 1;
        if self.writes == 1 {
            return Err(std::io::ErrorKind::Interrupted.into());
        }
        if self.full {
            return Ok(0);
        }
        let before = self.taken.len();
        for buf in bufs {
            let room = self.limit - (self.taken.len() - before);
            self.taken.extend_from_slice(&buf[..buf.len().min(room)]);
        }
        Ok(self.taken.len() - before)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

// No NCBI counterpart: the report's write of two pieces (`report::ostream::write_device2`) is
// `write_all` of the two in one write where the device takes them: an interrupted write is made
// again, a device that takes part of a write gets the rest (each write's bytes counted), and one
// that takes no byte is `write_all`'s error.
#[test]
fn report_writes_two_pieces_as_write_all() {
    let device = |limit, full| Device {
        taken: Vec::new(),
        writes: 0,
        limit,
        full,
    };
    let write = |device: &mut Device, first: &[u8], second: &[u8]| {
        let mut landed = 0;
        crate::report::ostream::write_device2(device, &mut landed, first, second).map(|()| landed)
    };
    let mut one = device(usize::MAX, false);
    assert_eq!(write(&mut one, b"abc", b"defg").unwrap(), 7);
    assert_eq!((one.taken.as_slice(), one.writes), (&b"abcdefg"[..], 2));
    let mut parts = device(2, false);
    assert_eq!(write(&mut parts, b"abc", b"defg").unwrap(), 7);
    assert_eq!((parts.taken.as_slice(), parts.writes), (&b"abcdefg"[..], 5));
    let mut second = device(usize::MAX, false);
    assert_eq!(write(&mut second, b"", b"defg").unwrap(), 4);
    assert_eq!((second.taken.as_slice(), second.writes), (&b"defg"[..], 2));
    let error = write(&mut device(usize::MAX, true), b"a", b"b").unwrap_err();
    let expected = device(usize::MAX, true).write_all(b"ab").unwrap_err();
    assert_eq!(
        (error.kind(), error.to_string()),
        (expected.kind(), expected.to_string())
    );
}
