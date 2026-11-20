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
        assert!(
            kept <= pending.len(),
            "{ops}: held {kept} of {}",
            pending.len()
        );
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

/// The standard output writes (descriptor 1) of a row of `ncbi_rpsbproc_writes.tsv`.
fn rpsbproc_stdout_writes(case: &str, destination: &str) -> Vec<usize> {
    let table = include_str!("../../tests/fixtures/stdio_buffer/ncbi_rpsbproc_writes.tsv");
    let row = table
        .lines()
        .filter(|line| !line.starts_with('#'))
        .find(|line| line.starts_with(&format!("{case}\t{destination}\t")))
        .unwrap();
    row.split('\t')
        .nth(3)
        .unwrap()
        .split(' ')
        .filter_map(|write| write.strip_prefix("1:"))
        .map(|size| size.parse().unwrap())
        .collect()
}

#[test]
fn rpsbprocs_report_with_a_cut_line_is_written_to_a_pipe_or_file_where_ncbi_writes_it() {
    // `-logfile /dev/stdout` on a database name of 9,000 bytes (`longdb`): NCBI's stdout writes
    // 1127 bytes (the report so far, written by the `cerr` line that follows), then 4096 (a block
    // of the SESSION line's long string), then the rest (930) at exit; the two diagnostics of the
    // exception go to the other descriptor of the same pipe between them. `oracle` holds the
    // bytes of the pipe: the report with those 296 bytes at offset 1127 + 4096. The insertions
    // are the report's fields and delimiters (a field, a tab, ..., a line end), which the table
    // does not record: this replays the sizes of the writes, not NCBI's insertions.
    let oracle: &[u8] =
        include_bytes!("../../tests/fixtures/rpsbproc/oracle/linux_cout_longdb.out");
    let report = [&oracle[..1127 + 4096], &oracle[1127 + 4096 + 296..]].concat();
    for destination in ["pipe", "file"] {
        let observed = rpsbproc_stdout_writes("longdb", destination);
        assert_eq!(observed, [1127, 4096, 930], "{destination}");
        let mut buffer = NcbiBuffer::new(NcbiStdio::Cout {
            block: 4096,
            lines: false,
        });
        let mut pending: Vec<u8> = Vec::new();
        let mut writes = Vec::new();
        let insert = |buffer: &mut NcbiBuffer,
                      pending: &mut Vec<u8>,
                      writes: &mut Vec<usize>,
                      bytes: &[u8]| {
            pending.extend_from_slice(bytes);
            buffer.sputn(bytes.len(), pending);
            let out = pending.len() - buffer.held();
            if out > 0 {
                writes.push(out);
                pending.drain(..out);
            }
        };
        // The head, then the `cerr` line's flush.
        insert(&mut buffer, &mut pending, &mut writes, &report[..1127]);
        buffer.flush();
        writes.push(pending.len());
        pending.clear();
        // The rest, one insertion per field (the string of `D`s is one insertion).
        for piece in report[1127..].split_inclusive(|byte| *byte == b'\t' || *byte == b'\n') {
            insert(&mut buffer, &mut pending, &mut writes, piece);
        }
        buffer.flush();
        writes.push(pending.len());
        assert_eq!(writes, observed, "{destination}");
    }
}

#[test]
fn rpsbprocs_terminal_writes_are_the_lines_and_blocks_of_a_line_buffered_stream() {
    // The same run on a pty (block 1024, line-buffered): the table's writes cannot be replayed
    // through the entry points, because NCBI's insertions are not recorded (a line of several
    // fields is written, on a terminal, at its line end or where 1024 bytes are held, and the
    // observed sizes 1024, 18, 9, 71, 5 split the head at the points of its insertions). What
    // is checked here: the writes add up to the pipe's, the head's five writes, and that no
    // write is longer than the block of the pipe.
    let tty = rpsbproc_stdout_writes("longdb", "pty");
    let pipe = rpsbproc_stdout_writes("longdb", "pipe");
    assert_eq!(tty.iter().sum::<usize>(), pipe.iter().sum::<usize>());
    assert_eq!(tty[..5], [1024, 18, 9, 71, 5]);
    assert!(tty.iter().all(|&size| size <= 4096));
}
