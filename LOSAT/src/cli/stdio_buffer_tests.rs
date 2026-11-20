//! The tables of `tests/fixtures/stdio_buffer/` replayed through `NcbiBuffer`: each row is a
//! list of operations on a stream and the number of bytes that had reached the file or device
//! after each of them (`strace` of the program that made them; the header of each table says
//! which program, binary and host).

use super::{stdio_buffer_size, NcbiBuffer, NcbiStdio, NCBI_BUFSIZ};

/// The operations of a row: `W<n>[:<k>...]` an insertion of `n` bytes (newlines at the offsets
/// `k`), `C` the insertion of a character, `L` of a newline character, `F` a flush, `E` a write
/// to `cerr` (which flushes `cout`, to which it is tied; the `-out` file's stream is not), `N`
/// `endl`.
fn replay(stream: NcbiStdio, ops: &str) -> Vec<u64> {
    let mut buffer = NcbiBuffer::new(stream);
    let mut pending: Vec<u8> = Vec::new();
    let mut written = 0u64;
    let mut after = Vec::new();
    for op in ops.split(' ') {
        let mut parts = op[1..].split(':');
        match op.as_bytes()[0] {
            b'W' => {
                let n: usize = parts.next().unwrap().parse().unwrap();
                let mut bytes = vec![b'x'; n];
                for at in parts {
                    let at: usize = at.parse().unwrap();
                    if at < n {
                        bytes[at] = b'\n';
                    }
                }
                pending.extend_from_slice(&bytes);
                buffer.sputn(n, &pending);
            }
            c @ (b'C' | b'L') => {
                pending.push(if c == b'C' { b'x' } else { b'\n' });
                buffer.sputc(&pending);
            }
            b'N' => {
                pending.push(b'\n');
                buffer.sputc(&pending);
                buffer.flush();
            }
            b'F' => buffer.flush(),
            b'E' => {
                if matches!(stream, NcbiStdio::Cout { .. }) {
                    buffer.flush();
                }
            }
            _ => panic!("operation {op}"),
        }
        let kept = buffer.held();
        assert!(kept <= pending.len(), "{ops}: held {kept} of {}", pending.len());
        let out = pending.len() - kept;
        written += out as u64;
        pending.drain(..out);
        after.push(written);
    }
    after
}

fn check(table: &str, stream: NcbiStdio) -> usize {
    let mut rows = 0;
    for line in table.lines().filter(|line| !line.starts_with('#')) {
        let mut fields = line.split('\t');
        let ops = fields.next().unwrap();
        let observed: Vec<u64> = fields
            .next()
            .unwrap()
            .split(',')
            .map(|n| n.parse().unwrap())
            .collect();
        assert_eq!(replay(stream, ops), observed, "{stream:?}: {ops}");
        rows += 1;
    }
    rows
}

#[test]
fn a_pipe_or_file_as_standard_output_writes_in_blocks_of_4096() {
    // `cout` (C's `stdout`) on a pipe and on a file, `st_blksize` 4096.
    let rows = check(
        include_str!("../../tests/fixtures/stdio_buffer/cout_pipe.tsv"),
        NcbiStdio::Cout {
            block: 4096,
            lines: false,
        },
    );
    assert!(rows > 1000);
}

#[test]
fn a_file_as_standard_output_writes_in_blocks_of_4096() {
    let rows = check(
        include_str!("../../tests/fixtures/stdio_buffer/cout_file.tsv"),
        NcbiStdio::Cout {
            block: 4096,
            lines: false,
        },
    );
    assert!(rows > 1000);
}

#[test]
fn a_terminal_as_standard_output_writes_in_blocks_of_1024_and_at_newlines() {
    // `cout` on a pty, `st_blksize` 1024, line-buffered.
    let rows = check(
        include_str!("../../tests/fixtures/stdio_buffer/cout_tty.tsv"),
        NcbiStdio::Cout {
            block: 1024,
            lines: true,
        },
    );
    assert!(rows > 1000);
}

#[test]
fn the_out_file_writes_when_it_reaches_8191_bytes() {
    // `std::ofstream` on a file.
    let rows = check(
        include_str!("../../tests/fixtures/stdio_buffer/ofstream_file.tsv"),
        NcbiStdio::File,
    );
    assert!(rows > 1000);
}

#[test]
fn a_new_buffer_has_written_nothing_and_a_flush_empties_it() {
    for stream in [
        NcbiStdio::File,
        NcbiStdio::Cout {
            block: 4096,
            lines: false,
        },
    ] {
        let mut buffer = NcbiBuffer::new(stream);
        assert_eq!((buffer.held(), buffer.stream()), (0, stream));
        buffer.sputn(10, &[b'a'; 10]);
        assert_eq!(buffer.held(), 10);
        buffer.flush();
        assert_eq!(buffer.held(), 0);
    }
}

#[test]
fn a_char_on_standard_output_is_held_in_a_buffer_that_it_makes() {
    // A first `put` makes the buffer and holds the char; a buffer of 4096 holds 4096 chars and
    // writes them with the next one.
    let mut buffer = NcbiBuffer::new(NcbiStdio::Cout {
        block: 4096,
        lines: false,
    });
    let mut stream = Vec::new();
    for count in 1..=4096usize {
        stream.push(b'x');
        buffer.sputc(&stream);
        assert_eq!(buffer.held(), count);
    }
    stream.push(b'x');
    buffer.sputc(&stream);
    assert_eq!(buffer.held(), 1);
}

#[test]
fn the_size_of_standard_outputs_buffer_is_below_bufsiz_or_bufsiz() {
    let size = stdio_buffer_size();
    assert!(size > 0 && size <= NCBI_BUFSIZ, "{size}");
    // A terminal's, a pipe's and a file's are what `st_blksize` says (1024, 4096).
    assert!(matches!(NcbiStdio::cout(), NcbiStdio::Cout { block, .. } if block == size));
}

#[test]
fn ncbis_tabular_row_is_written_to_the_out_file_where_blastn_writes_it() {
    // blastn's row `qseqid \t sseqid \n` with ids of L bytes: four insertions; the write calls
    // of NCBI BLAST+ 2.17.0 (official and conda builds) are the steps of the bytes written
    // after each insertion, and the held rest at the close.
    let table = include_str!("../../tests/fixtures/stdio_buffer/ncbi_blastn_out.tsv");
    let mut rows = 0;
    for line in table.lines().filter(|line| !line.starts_with('#')) {
        let (len, writes) = line.split_once('\t').unwrap();
        let len: usize = len.parse().unwrap();
        let mut buffer = NcbiBuffer::new(NcbiStdio::File);
        let mut pending = 0usize;
        let mut sizes = Vec::new();
        for n in [len, 1, len, 1] {
            pending += n;
            buffer.sputn(n, &vec![b'x'; pending]);
            let out = pending - buffer.held();
            if out > 0 {
                sizes.push(out);
                pending -= out;
            }
        }
        if pending > 0 {
            sizes.push(pending);
        }
        let sizes: Vec<String> = sizes.iter().map(usize::to_string).collect();
        assert_eq!(sizes.join(","), writes, "L = {len}");
        rows += 1;
    }
    assert_eq!(rows, 12);
}
