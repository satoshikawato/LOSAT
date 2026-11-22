use super::*;
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
use std::io::Cursor;
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:727-735
trait NativeRead: Read + Seek {
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
    fn available(&mut self) -> std::io::Result<usize>;
}
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:727-735
impl NativeRead for crate::blastinput::input_files::InputFile {
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
    fn available(&mut self) -> std::io::Result<usize> {
        (crate::blastinput::fasta_reader::stream::InputBackend::of(self).available)(self)
    }
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
struct Input {
    data: Cursor<Vec<u8>>,
    seekable: bool,
    fault: Option<u64>,
    error: bool,
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
impl Read for Input {
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
        if self.fault.is_some_and(|p| self.data.position() >= p) {
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
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-871
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
impl NativeRead for Input {
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
    fn available(&mut self) -> std::io::Result<usize> {
        Ok(0)
    }
}
impl Seek for Input {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-871
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
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:154-174
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-873
#[test]
fn pinned_cpp_stream_state_and_line_oracle() {
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:353-381
    // Decode lossless independent C++ fixture bytes; no Rust-generated
    // expected lines or normalization of actual public report streams.
    let decode = |encoded: &str| -> Vec<u8> {
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
                .map(|b| u8::from_str_radix(std::str::from_utf8(b).unwrap(), 16).unwrap())
                .collect()
        }
    };
    for row in include_str!("../../../tests/unit/blastx_stage_e_io_stream_expected.tsv").lines() {
        let fields: Vec<_> = row.split('\t').collect();
        assert_eq!(fields.len(), 9);
        let data = decode(fields[1]);
        let count = fields[8].parse::<usize>().unwrap();
        let expected: Vec<Vec<u8>> = if count == 0 {
            Vec::new()
        } else {
            fields[7].split('|').map(decode).collect()
        };
        assert_eq!(expected.len(), count);
        // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:353-381
        if fields[0].starts_with("memory_") {
            let lines: Vec<_> = FastaMemoryLines { remaining: &data }.collect();
            assert_eq!(
                lines,
                expected.iter().map(Vec::as_slice).collect::<Vec<_>>(),
                "{} memory lines",
                fields[0]
            );
            continue;
        }
        let fault: i64 = fields[4].parse().unwrap();
        let file_backend = fields[0].starts_with("file_");
        for &error in if file_backend {
            &[false][..]
        } else {
            &[false, true][..]
        } {
            let temporary = crate::utils::test_dir::temp_dir().join(format!(
                "losatx-e-stream-unit-{}-{}",
                crate::utils::test_dir::process_id(),
                fields[0]
            ));
            let input: Box<dyn NativeRead> = if file_backend {
                assert!(!temporary.exists());
                std::fs::write(&temporary, &data).unwrap();
                Box::new(crate::blastinput::input_files::InputFile::from(
                    std::fs::File::open(&temporary).unwrap(),
                ))
            } else {
                Box::new(Input {
                    data: Cursor::new(data.clone()),
                    seekable: fields[2] == "1",
                    fault: if fault < 0 { None } else { Some(fault as u64) },
                    error,
                })
            };
            let mut stream = FastaStream::new(input);
            // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
            stream.backend_available = Some(|input| input.available());
            let empty = fields[3] == "1" && stream_is_empty(&mut stream);
            assert_eq!(empty, fields[6] == "1", "{} empty error={error}", fields[0]);
            let lines: Vec<Vec<u8>> = if empty { Vec::new() } else { stream.collect() };
            if file_backend {
                std::fs::remove_file(&temporary).unwrap();
            }
            assert_eq!(lines, expected, "{} lines error={error}", fields[0]);
        }
    }
}
