//! The stream against the independent C++ oracle of the pinned `CStreamLineReader`,
//! `NcbiGetline`, `CPushback_Streambuf` and `IsIStreamEmpty`
//! (`LOSAT/tests/unit/blastx_stage_e_io_stream_expected.tsv`, the table that BLASTX's copy
//! of this stream is tested with: `algorithm/blastx/input.rs`,
//! `pinned_cpp_stream_state_and_line_oracle`).
//!
//! Columns: name, input bytes (hex, or `~byte:count` runs), seekable, run the empty
//! check, read-fault offset (-1 for none), unused, `IsIStreamEmpty`, the lines joined by
//! `|`, the number of lines. `memory_*` rows are `CMemoryLineReader`, which the BLAST input
//! source does not use; `file_*` rows read a real file (`FastaStream::from_file`).

use super::*;
use std::io::Cursor;

/// The input of a table row, with `CNcbiIfstream`'s readable-byte hint (`in_avail`).
trait NativeRead: Read + Seek {
    fn available(&mut self) -> std::io::Result<usize>;
}

impl NativeRead for std::fs::File {
    fn available(&mut self) -> std::io::Result<usize> {
        file_available(self)
    }
}

/// A stream that reads one byte at a time, may refuse to seek (a pipe) and fails or ends
/// at `fault` (the oracle's `filebuf` read failure).
struct Input {
    data: Cursor<Vec<u8>>,
    seekable: bool,
    fault: Option<u64>,
    error: bool,
}

impl Read for Input {
    fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
        if self
            .fault
            .is_some_and(|fault| self.data.position() >= fault)
        {
            return if self.error {
                Err(std::io::Error::other("filebuf read failure"))
            } else {
                Ok(0)
            };
        }
        let length = bytes.len().min(1);
        self.data.read(&mut bytes[..length])
    }
}

impl Seek for Input {
    fn seek(&mut self, position: SeekFrom) -> std::io::Result<u64> {
        if self.seekable {
            self.data.seek(position)
        } else {
            Err(std::io::Error::new(
                std::io::ErrorKind::Unsupported,
                "nonseekable stream",
            ))
        }
    }
}

impl NativeRead for Input {
    fn available(&mut self) -> std::io::Result<usize> {
        Ok(0)
    }
}

/// The lines of `CStreamLineReader` (`operator++` until `AtEOF`).
fn lines_of<R: Read>(stream: &mut FastaStream<R>) -> Vec<Vec<u8>> {
    let mut lines = Vec::new();
    let mut line = Vec::new();
    while !stream.at_eof() {
        stream.advance(&mut line);
        lines.push(line.clone());
    }
    lines
}

fn decode(encoded: &str) -> Vec<u8> {
    if let Some(runs) = encoded.strip_prefix('~') {
        let mut bytes = Vec::new();
        for run in runs.split(',') {
            let (byte, count) = run.split_once(':').unwrap();
            let byte = u8::from_str_radix(byte, 16).unwrap();
            bytes.extend(std::iter::repeat_n(byte, count.parse::<usize>().unwrap()));
        }
        bytes
    } else {
        encoded
            .as_bytes()
            .chunks(2)
            .map(|pair| u8::from_str_radix(std::str::from_utf8(pair).unwrap(), 16).unwrap())
            .collect()
    }
}

#[test]
fn stream_matches_the_pinned_cpp_oracle() {
    let table = include_str!("../../../tests/unit/blastx_stage_e_io_stream_expected.tsv");
    let (mut checked, mut memory_rows) = (0usize, 0usize);
    for row in table.lines() {
        let fields: Vec<_> = row.split('\t').collect();
        assert_eq!(fields.len(), 9, "{row}");
        if fields[0].starts_with("memory_") {
            memory_rows += 1;
            continue;
        }
        let data = decode(fields[1]);
        let count = fields[8].parse::<usize>().unwrap();
        let expected: Vec<Vec<u8>> = if count == 0 {
            Vec::new()
        } else {
            fields[7].split('|').map(decode).collect()
        };
        assert_eq!(expected.len(), count, "{}", fields[0]);
        let fault: i64 = fields[4].parse().unwrap();
        let file_backend = fields[0].starts_with("file_");
        let errors: &[bool] = if file_backend {
            &[false]
        } else {
            &[false, true]
        };
        // With and without LOSAT's bulk copy of lines (`FastaStream::bulk`).
        for (&error, bulk) in errors
            .iter()
            .flat_map(|error| [(error, true), (error, false)])
        {
            let temporary = std::env::temp_dir().join(format!(
                "losat-fasta-reader-stream-{}-{}",
                std::process::id(),
                fields[0]
            ));
            let input: Box<dyn NativeRead> = if file_backend {
                std::fs::write(&temporary, &data).unwrap();
                Box::new(std::fs::File::open(&temporary).unwrap())
            } else {
                Box::new(Input {
                    data: Cursor::new(data.clone()),
                    seekable: fields[2] == "1",
                    fault: u64::try_from(fault).ok(),
                    error,
                })
            };
            let mut stream = FastaStream::new(input);
            stream.backend_available = Some(|input| input.available());
            stream.bulk = bulk;
            let empty = fields[3] == "1" && stream_is_empty(&mut stream);
            assert_eq!(
                empty,
                fields[6] == "1",
                "{} empty error={error} bulk={bulk}",
                fields[0]
            );
            let lines = if empty {
                Vec::new()
            } else {
                lines_of(&mut stream)
            };
            if file_backend {
                std::fs::remove_file(&temporary).unwrap();
            }
            assert_eq!(
                lines, expected,
                "{} lines error={error} bulk={bulk}",
                fields[0]
            );
            checked += 1;
        }
    }
    // 9674 rows read through this stream (1413 of them real files, the others with and
    // without a read error), each with and without the bulk copy; 2506 CMemoryLineReader
    // rows skipped.
    assert_eq!(
        (checked, memory_rows),
        (2 * (1413 + 2 * (12180 - 2506 - 1413)), 2506)
    );
}

/// `FastaStream::from_file` (the CLI) and `FastaStream::from_bytes` (the ABI) on the
/// table's real-file rows: the same lines as the oracle's `ifstream`, without the test's
/// `in_avail` hook.
#[test]
fn file_rows_through_from_file_and_from_bytes() {
    let table = include_str!("../../../tests/unit/blastx_stage_e_io_stream_expected.tsv");
    for row in table.lines().filter(|row| row.starts_with("file_")) {
        let fields: Vec<_> = row.split('\t').collect();
        let data = decode(fields[1]);
        let count = fields[8].parse::<usize>().unwrap();
        let expected: Vec<Vec<u8>> = if count == 0 {
            Vec::new()
        } else {
            fields[7].split('|').map(decode).collect()
        };
        let temporary = std::env::temp_dir().join(format!(
            "losat-fasta-reader-from-file-{}-{}",
            std::process::id(),
            fields[0]
        ));
        std::fs::write(&temporary, &data).unwrap();
        let mut stream = FastaStream::from_file(std::fs::File::open(&temporary).unwrap());
        let lines = lines_of(&mut stream);
        std::fs::remove_file(&temporary).unwrap();
        assert_eq!(lines, expected, "{} from_file", fields[0]);
        let mut stream = FastaStream::from_bytes(&data);
        assert_eq!(lines_of(&mut stream), expected, "{} from_bytes", fields[0]);
    }
}
