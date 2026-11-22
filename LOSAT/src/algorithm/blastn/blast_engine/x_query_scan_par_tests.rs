use super::super::x_query_drop::Store;
use super::*;
use crate::blastinput::fasta_reader::{FastaInputSource, FastaRecord, ParseErrorCode, ReadError};
use crate::utils::threading::SearchPool;
use std::io::Write;

/// The first reading as `x_query_stream::scan` makes it (without the archive): the lengths,
/// and whether it ended at the end of the input.
fn serial(bytes: &[u8], config: ReaderConfig) -> (Vec<usize>, bool) {
    let mut source = FastaInputSource::from_bytes(bytes, config);
    let mut lengths = Vec::new();
    loop {
        if source.end() {
            return (lengths, true);
        }
        match source.next_sequence(&mut |_| Ok(())) {
            Ok(record) => lengths.push(record.sequence.len()),
            Err(ReadError::Parse {
                code: ParseErrorCode::Eof,
                ..
            }) => return (lengths, false),
            Err(_) => return (lengths, false),
        }
    }
}

fn file_of(bytes: &[u8], tag: &str) -> (std::path::PathBuf, File) {
    let path = crate::utils::test_dir::temp_dir().join(format!(
        "losat_qscanpar_{}_{tag}_{:?}.fa",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    let mut file = File::create(&path).expect("create the test file");
    file.write_all(bytes).expect("write the test file");
    drop(file);
    let file = File::open(&path).expect("open the test file");
    (path, file)
}

/// The parallel reading of `bytes` from byte `start`, without a minimum size.
fn parallel(
    bytes: &[u8],
    start: u64,
    config: ReaderConfig,
    threads: usize,
    limits: Limits,
) -> Option<Vec<usize>> {
    let (path, file) = file_of(bytes, "p");
    let lengths = lengths_with(&file, start, config, threads, 0, limits).map(|scan| scan.lengths);
    let _ = std::fs::remove_file(path);
    lengths
}

/// A small deterministic generator (xorshift).
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
}

/// A FASTA text with the reader's cases: ids and none, titles, lower case, IUPAC letters, `N`
/// runs, digits and spaces in data lines, blank and comment lines, deflines without letters,
/// data before the first defline, lines of any length and no final newline; in one text of four
/// `>?` gap lines and `>?_` deflines, in one of eight a line of bad letters (a reader error).
fn random_fasta(rng: &mut Rng, records: usize) -> Vec<u8> {
    let mut text = Vec::new();
    let gaps = rng.below(4) == 0;
    let bad = rng.below(8) == 0;
    if rng.below(4) == 0 {
        text.extend_from_slice(b"ACGTNACGT\n");
    }
    for index in 0..records {
        match rng.below(6) {
            0 => text.extend_from_slice(b">\n"),
            1 => text.extend_from_slice(format!(">lcl|q{index} title {index}\n").as_bytes()),
            2 => text.extend_from_slice(format!(">q{index}\n").as_bytes()),
            3 if gaps => text.extend_from_slice(b">?_x a title\n"),
            _ => text.extend_from_slice(format!(">gnl|db|{index} some words\n").as_bytes()),
        }
        let lines = rng.below(8);
        for _ in 0..lines {
            match rng.below(14) {
                0 => text.push(b'\n'),
                1 => text.extend_from_slice(b"; a comment\n"),
                2 if gaps => text.extend_from_slice(b">?10\n"),
                3 => text.extend_from_slice(b"  12 acgt acgt 13\n"),
                4 if bad && rng.below(4) == 0 => text.extend_from_slice(b"%%%%%%%%%%%%\n"),
                _ => {
                    let width = 1 + rng.below(90) as usize;
                    for _ in 0..width {
                        let letters = b"ACGTacgtNNnRYKMSWBDHV";
                        text.push(letters[rng.below(letters.len() as u64) as usize]);
                    }
                    text.push(b'\n');
                }
            }
        }
    }
    if !text.is_empty() && rng.below(3) == 0 {
        text.pop();
    }
    text
}

const SMALL: Limits = Limits {
    max_span: 97,
    max_reach: 400,
    parts_per_thread: 3,
};

#[test]
fn the_parts_give_the_first_readings_lengths() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let (mut taken, mut refused, mut failed) = (0, 0, 0);
    for round in 0..300 {
        let loaders = round % 2 == 0;
        let config = ReaderConfig::query("BLASTN", false, loaders);
        let text = random_fasta(&mut rng, 1 + round % 40);
        let (expected, clean) = serial(&text, config);
        for threads in [2, 3, 8] {
            for limits in [
                SMALL,
                Limits {
                    max_span: 7,
                    max_reach: 2000,
                    parts_per_thread: 1,
                },
                Limits::DEFAULT,
            ] {
                match parallel(&text, 0, config, threads, limits) {
                    Some(lengths) => {
                        assert!(clean, "round {round}: parts taken where the reading fails");
                        assert!(!has_gap_mark(&text), "round {round}: parts taken with `>?`");
                        assert_eq!(lengths, expected, "round {round} threads {threads}");
                        taken += 1;
                    }
                    None => refused += 1,
                }
                if !clean {
                    failed += 1;
                }
            }
        }
    }
    assert!(taken > 1000, "the parts were taken {taken} times");
    assert!(
        refused > 300 && failed > 30,
        "refused {refused}, failed readings {failed}"
    );
}

/// The second reading of each part (`read_again`) gives the records of the first reading's
/// parts as the file's reader reads them (their local IDs, titles, letters and messages, with
/// the file's line numbers), and the lines before each part are the line ends before it.
#[test]
fn the_parts_read_again_give_the_readers_records() {
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    let (mut taken, mut warned) = (0, 0);
    for round in 0..200 {
        let config = ReaderConfig::query("BLASTN", false, round % 2 == 0);
        let text = random_fasta(&mut rng, 1 + round % 40);
        let (expected, clean) = serial_records(&text, config);
        let (path, file) = file_of(&text, "again");
        let limits = Limits {
            max_span: 7,
            max_reach: 2000,
            parts_per_thread: 1,
        };
        for (threads, limits) in [(2, SMALL), (3, limits)] {
            let Some(scan) = lengths_with(&file, 0, config, threads, 0, limits) else {
                continue;
            };
            assert!(clean, "round {round}: parts taken where the reading fails");
            let mut got: Vec<FastaRecord> = Vec::new();
            for (index, span) in scan.spans.iter().enumerate() {
                let ends = text[..span.bytes.start as usize]
                    .iter()
                    .filter(|&&byte| byte == b'\n')
                    .count();
                assert_eq!(span.line, ends as u64, "round {round}");
                assert_eq!(span.first, got.len(), "round {round}");
                let last = index + 1 == scan.spans.len();
                let (records, failure) = read_again(&file, span, last, config, &scan.lengths);
                assert_eq!(failure, None, "round {round}");
                assert_eq!(records.len(), span.records, "round {round}");
                if index > 0 && records.iter().any(|record| !record.warnings.is_empty()) {
                    warned += 1;
                }
                got.extend(records);
            }
            assert_eq!(got, expected, "round {round} threads {threads}");
            taken += 1;
        }
        let _ = std::fs::remove_file(path);
    }
    assert!(taken > 100 && warned > 20, "taken {taken}, warned {warned}");
}

/// The records of `bytes` as `read_queries` reads them, and whether the reading ended at the
/// end of the input.
fn serial_records(bytes: &[u8], config: ReaderConfig) -> (Vec<FastaRecord>, bool) {
    use crate::blastinput::fasta_reader::{read_queries, QueryEnd};
    let records = read_queries(&mut FastaInputSource::from_bytes(bytes, config));
    (records.records, matches!(records.end, QueryEnd::Input))
}

/// The records of `bytes` read in parts from byte 0 into a store (without a minimum size), as
/// the loop takes them (in `batches` batches, by a pool of `threads` threads).
fn parallel_records(
    bytes: &[u8],
    config: ReaderConfig,
    threads: usize,
    limits: Limits,
    batches: usize,
) -> Option<Vec<FastaRecord>> {
    let (path, file) = file_of(bytes, "r");
    let parts = store_parts_with(&file, 0, config, threads, 0, limits);
    let _ = std::fs::remove_file(path);
    let store = Store::from_parts(parts?, config.id_prefix, true);
    let count = store.lengths().len();
    let step = count.div_ceil(batches.max(1)).max(1);
    let take = |pool: &SearchPool<'_>| {
        let mut records = Vec::new();
        for from in (0..count).step_by(step) {
            records.extend(store.take(from..(from + step).min(count), pool)?);
        }
        Ok(records)
    };
    Some(if threads > 1 {
        crate::utils::threading::with_search_pool(threads, "blastn", take).expect("take")
    } else {
        take(&SearchPool::x_serial()).expect("take")
    })
}

#[test]
fn the_parts_of_the_store_give_the_readings_records() {
    let mut rng = Rng(0x6a09_e667_f3bc_c908);
    let (mut taken, mut refused) = (0, 0);
    for round in 0..300 {
        let config = ReaderConfig::query("BLASTN", false, round % 2 == 0);
        let text = random_fasta(&mut rng, 1 + round % 40);
        let (expected, clean) = serial_records(&text, config);
        for threads in [2, 3, 8] {
            for limits in [
                SMALL,
                Limits {
                    max_span: 7,
                    max_reach: 2000,
                    parts_per_thread: 1,
                },
                Limits::DEFAULT,
            ] {
                let batches = 1 + round % 4;
                match parallel_records(&text, config, threads, limits, batches) {
                    Some(records) => {
                        assert!(clean, "round {round}: parts taken where the reading fails");
                        assert!(!has_gap_mark(&text), "round {round}: parts taken with `>?`");
                        assert_eq!(records, expected, "round {round} threads {threads}");
                        taken += 1;
                    }
                    None => refused += 1,
                }
            }
        }
    }
    // Most texts have a record with a reader message after the first part, a `>?` or a reader
    // error: the search's reader then reads them.
    assert!(taken > 150, "the parts were taken {taken} times");
    assert!(refused > 150, "refused {refused}");
}

#[test]
fn a_reader_message_after_the_first_part_leaves_the_records_to_the_reader() {
    let config = ReaderConfig::query("BLASTN", false, true);
    let mut good = Vec::new();
    for index in 0..40 {
        good.extend_from_slice(format!(">q{index} t\nACGTACGTAC\nACgtNNNNNNACGT\n").as_bytes());
    }
    let (expected, clean) = serial_records(&good, config);
    assert!(clean && expected.iter().all(|record| record.warnings.is_empty()));
    let got = parallel_records(&good, config, 4, SMALL, 3).expect("the parts are taken");
    assert_eq!(got, expected);
    assert_eq!(got[39].local_id, "Query_40");
    // A data line with a hyphen is read with a message that names its line: in the first part
    // the records are kept, after it the reader reads them.
    let mut late = good.clone();
    late.extend_from_slice(b">z\nACGT-ACGT\n");
    let (expected, clean) = serial_records(&late, config);
    assert!(clean && !expected.last().expect("a record").warnings.is_empty());
    assert_eq!(parallel_records(&late, config, 4, SMALL, 2), None);
    let mut early = b">z\nACGT-ACGT\n".to_vec();
    early.extend_from_slice(&good);
    let (expected, _) = serial_records(&early, config);
    assert!(!expected[0].warnings.is_empty());
    assert_eq!(
        parallel_records(&early, config, 4, SMALL, 5),
        Some(expected)
    );
    // One thread, or the minimum: the search's reader reads them.
    let (path, file) = file_of(&good, "rmin");
    assert!(store_parts_with(&file, 0, config, 1, 0, SMALL).is_none());
    assert!(store_parts_with(&file, 0, config, 4, 1 << 20, SMALL).is_none());
    let _ = std::fs::remove_file(path);
}

fn has_gap_mark(text: &[u8]) -> bool {
    text.windows(2).any(|pair| pair == b">?")
}

/// The record start of `record_start` found by looking at every byte.
fn record_start_by_bytes(text: &[u8], target: usize, reach: usize) -> Option<usize> {
    if target >= text.len() {
        return Some(text.len());
    }
    let limit = text.len().min(target + reach);
    (target..limit)
        .find(|&at| text[at] == b'>' && text[at - 1] == b'\n' && text.get(at + 1) != Some(&b'?'))
        .or((limit == text.len()).then_some(text.len()))
}

#[test]
fn record_starts_across_blocks_are_those_of_every_byte() {
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    for round in 0..6 {
        let mut text = random_fasta(&mut rng, 120 + 40 * round);
        // Record starts at and around the block edges of the search (8,192 bytes from a target).
        for edge in [8_190usize, 8_191, 8_192, 8_193, 16_383, 16_384] {
            if edge + 2 < text.len() {
                text[edge - 1] = b'\n';
                text[edge] = b'>';
            }
        }
        let (path, file) = file_of(&text, "rsb");
        let end = text.len() as u64;
        for target in (1..text.len())
            .step_by(37)
            .chain([1, 2, 8_190, 8_191, 8_192, 8_193])
        {
            for reach in [1usize, 50, 8_191, 8_192, 8_193, 20_000, 1 << 20] {
                assert_eq!(
                    record_start(&file, 0, target as u64, end, reach as u64),
                    record_start_by_bytes(&text, target, reach).map(|at| at as u64),
                    "round {round} target {target} reach {reach}"
                );
            }
        }
        let _ = std::fs::remove_file(path);
    }
}

#[test]
fn a_part_refuses_gap_marks_across_reads_and_bytes_it_cannot_read() {
    let config = ReaderConfig::query("BLASTN", false, true);
    let failed = std::sync::atomic::AtomicBool::new(false);
    // The reader's input window is 8,191 bytes: a `>` as its last byte, `?` the next read's first.
    let mut text = b">a\n".to_vec();
    while text.len() < 8_189 {
        text.extend_from_slice(b"ACGTACGTAC\n");
    }
    text.truncate(8_189);
    text.push(b'\n');
    let mut gap = text.clone();
    gap.extend_from_slice(b">?5\nACGT\n");
    assert_eq!(gap[8_190..8_192], *b">?");
    let mut plain = text.clone();
    plain.extend_from_slice(b">b\nACGT\n");
    for (bytes, taken) in [(&gap, false), (&plain, true)] {
        let (path, file) = file_of(bytes, "gap");
        let end = bytes.len() as u64;
        assert_eq!(read_part(&file, 0..end, config, &failed).is_some(), taken);
        // A part past the end of the file: the reader would take the early end for its own.
        assert_eq!(read_part(&file, 0..end + 100, config, &failed), None);
        let _ = std::fs::remove_file(path);
    }
}

#[test]
fn a_reading_that_fails_or_cannot_be_cut_is_left_to_the_first_reading() {
    let config = ReaderConfig::query("BLASTN", false, true);
    let good = b">a\nACGT\n>b\nACGTACGT\n>c\nAC\n>d\nACGTA\n".to_vec();
    assert_eq!(parallel(&good, 0, config, 4, SMALL), Some(vec![4, 8, 2, 5]));
    // A `\r` anywhere.
    let mut cr = good.clone();
    cr.splice(20..20, *b"\r");
    assert_eq!(parallel(&cr, 0, config, 4, SMALL), None);
    // A reader error in a late part (a line of mostly bad letters).
    let mut bad = good.clone();
    bad.extend_from_slice(b">e\n%%%%%%%%%%%%\n");
    assert!(!serial(&bad, config).1);
    assert_eq!(parallel(&bad, 0, config, 4, SMALL), None);
    // A record longer than the reach of a target.
    let mut long = b">a\n".to_vec();
    long.extend(std::iter::repeat_n(b'A', 5000));
    long.extend_from_slice(b"\n>b\nAC\n");
    let limits = Limits {
        max_span: 100,
        max_reach: 200,
        parts_per_thread: 1,
    };
    assert_eq!(parallel(&long, 0, config, 4, limits), None);
    // A Seq-id line first (the data loaders' check of the first line), and comment lines only
    // (an `eEOF`): the parts are taken only where the first reading ends at the end.
    for text in [
        b"NC_000913.3\n>a\nACGT\n>b\nAC\n".as_slice(),
        b"; one\n; two\n\n",
    ] {
        let (expected, clean) = serial(text, config);
        let parts = parallel(text, 0, config, 4, SMALL);
        assert_eq!(
            parts.is_some(),
            clean,
            "{:?}",
            String::from_utf8_lossy(text)
        );
        if let Some(lengths) = parts {
            assert_eq!(lengths, expected);
        }
    }
    assert!(!serial(b"; one\n; two\n\n", config).1);
    // One thread, no byte after the reader's start, or fewer bytes than the minimum.
    assert_eq!(parallel(&good, 0, config, 1, SMALL), None);
    assert_eq!(parallel(&good, good.len() as u64, config, 4, SMALL), None);
    let (path, file) = file_of(&good, "min");
    assert_eq!(lengths_with(&file, 0, config, 4, 1 << 20, SMALL), None);
    let _ = std::fs::remove_file(path);
}

#[test]
fn the_parts_start_where_the_reader_starts() {
    let config = ReaderConfig::query("BLASTN", false, true);
    let mut text = b"   \n\n".to_vec();
    let body = b">a\nACGT\n>b\nACGTACGT\n>c\nAC\n>d\nACGTA\n";
    text.extend_from_slice(body);
    let start = 5;
    assert_eq!(
        parallel(&text, start, config, 3, SMALL),
        Some(serial(&text[start as usize..], config).0)
    );
}

#[test]
fn record_starts_skip_gap_lines_and_lines_that_do_not_begin_with_the_mark() {
    let text = b">a\nAC\n>?5\n x>y\n>b\nAC\n".to_vec();
    let (path, file) = file_of(&text, "rs");
    let end = text.len() as u64;
    assert_eq!(record_start(&file, 0, 1, end, 100), Some(15));
    assert_eq!(record_start(&file, 0, 15, end, 100), Some(15));
    assert_eq!(record_start(&file, 0, 16, end, 100), Some(end));
    assert_eq!(record_start(&file, 0, 1, end, 5), None);
    assert_eq!(record_start(&file, 0, 1, end, 15), Some(15));
    let _ = std::fs::remove_file(path);
}
