use super::*;
use crate::blastinput::fasta_reader::{read_queries, QueryRecords};
use crate::blastinput::seq_range::{cut_queries, RangeRole};

/// The query input made from the lengths is `cut_queries`'s (with ranges that start inside,
/// at the end of, and past some records) and `QueryInput::whole`'s.
#[test]
fn query_input_from_lengths_is_that_of_the_records() {
    let lengths = [0usize, 1, 5, 99, 100, 101, 250, 3, 100, 7, 1000];
    let records: Vec<FastaRecord> = lengths
        .iter()
        .enumerate()
        .map(|(index, &length)| {
            let letters: Vec<u8> = (0..length)
                .map(|at| b"ACGTacgtN"[(at + index) % 9])
                .collect();
            FastaRecord::new(format!("Query_{}", index + 1), b"t", &letters)
        })
        .collect();
    let (whole, whole_lengths) = query_input_of(lengths.to_vec(), None);
    same_query_input(&whole, &QueryInput::whole(&records)).unwrap();
    assert_eq!(whole_lengths, lengths);
    for text in [
        "1-5",
        "3-200",
        "100-150",
        "101-5000",
        "6-7",
        "2-3",
        "1000-2000",
    ] {
        let range =
            crate::blastinput::seq_range::parse_sequence_range(text, RangeRole::Query).unwrap();
        let cut = cut_queries(&records, &range);
        let (input, searched) = query_input_of(lengths.to_vec(), Some(&range));
        same_query_input(&input, &cut.input).unwrap_or_else(|error| panic!("{text}: {error}"));
        let cut_lengths: Vec<usize> = cut.records.iter().map(|r| r.seq().len()).collect();
        assert_eq!(searched, cut_lengths, "{text}");
        // Each searched record is the input record cut at its placement and length.
        for (index, record) in cut.records.iter().enumerate() {
            let ordinal = input.ordinal(index);
            let from = input.placements.offset(index);
            assert_eq!(
                &records[ordinal].cut(from, from + searched[index]),
                record,
                "{text}"
            );
        }
    }
}

fn config() -> ReaderConfig {
    ReaderConfig::query("BLASTN", false, false)
}

/// A file of this test in the temporary directory, removed when dropped.
struct TempFile(std::path::PathBuf);

impl TempFile {
    fn new(tag: &str, bytes: &[u8]) -> Self {
        static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let path = crate::utils::test_dir::temp_dir().join(format!(
            "losat_qstream_{}_{}_{tag}.fa",
            crate::utils::test_dir::process_id(),
            NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        ));
        std::fs::write(&path, bytes).expect("write the test file");
        TempFile(path)
    }
}

impl Drop for TempFile {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// The store of `path` as `search_cli` makes it (`IsIStreamEmpty`, then the first reading)
/// with the switch on, and how the first reading ended (one thread: the reading in parts,
/// `x_query_scan_par.rs`, has its own tests).
fn store(path: &Path, range: Option<&SequenceRange>) -> (Store, QueryEnd) {
    let file = File::open(path).expect("open the test file");
    let opening = Opening::with_mode(true, path, &file, false, true, NO_RULE);
    assert!(opening.acts());
    let mut source = FastaInputSource::from_file(file, config());
    assert!(!source.stream_is_empty());
    Store::read(opening, &mut source, config(), range, 1).expect("the first reading")
}

/// The records as `search_cli` reads them without the switch.
fn old(path: &Path) -> QueryRecords {
    let mut source = FastaInputSource::from_file(File::open(path).unwrap(), config());
    assert!(!source.stream_is_empty());
    read_queries(&mut source)
}

/// xorshift64*.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

/// A FASTA text of `seed`: LF, CR LF or mixed line ends; records of 0 to about 20,000 letters
/// (one of 1.1 MiB in some texts: the reader's residue buffer is then taken, not copied), in
/// lines of 60 or 80 letters or one line (longer than the line reader's 8,191-byte window);
/// lower case runs, `N`s, IUPAC letters, digits (the reader's warning), records without a title
/// or without letters, blank and comment lines; letters before the first defline in some; and
/// the end with or without a line end, with blank lines, or with a record that the reader
/// refuses.
fn text(seed: u64) -> Vec<u8> {
    let mut rng = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
    let ends = seed % 3;
    let mut text = Vec::new();
    let eol = |text: &mut Vec<u8>, rng: &mut Rng| match ends {
        0 => text.push(b'\n'),
        1 => text.extend_from_slice(b"\r\n"),
        _ => {
            if rng.below(2) == 0 {
                text.push(b'\n')
            } else {
                text.extend_from_slice(b"\r\n")
            }
        }
    };
    if seed.is_multiple_of(7) {
        text.extend_from_slice(b"ACGTACGTTT");
        eol(&mut text, &mut rng);
    }
    let records = 1 + rng.below(40);
    for index in 0..records {
        match rng.below(10) {
            0 => {
                eol(&mut text, &mut rng);
            }
            1 => {
                text.extend_from_slice(b";a comment line");
                eol(&mut text, &mut rng);
            }
            _ => {}
        }
        text.push(b'>');
        if rng.below(8) != 0 {
            text.extend_from_slice(format!("rec{index} some words {seed}").as_bytes());
        }
        eol(&mut text, &mut rng);
        let length = if seed.is_multiple_of(5) && index == 2 {
            (1 << 20) + 100_000
        } else {
            match rng.below(20) {
                0 => 0,
                1 => 9_000 + rng.below(11_000) as usize,
                _ => rng.below(3_000) as usize,
            }
        };
        let width = match rng.below(3) {
            0 => 60,
            1 => 80,
            _ => length.max(1),
        };
        let mut letters: Vec<u8> = (0..length)
            .map(|_| b"ACGTACGTACGTNRY"[rng.below(15) as usize])
            .collect();
        if length > 200 && rng.below(3) == 0 {
            let at = rng.below(length as u64 - 150) as usize;
            letters[at..at + 120].make_ascii_lowercase();
        }
        for (line, chunk) in letters.chunks(width).enumerate() {
            text.extend_from_slice(chunk);
            if line == 1 && rng.below(6) == 0 {
                text.extend_from_slice(b"12");
            }
            eol(&mut text, &mut rng);
        }
    }
    match seed % 4 {
        0 => {
            // No line end after the last line.
            while matches!(text.last(), Some(b'\n' | b'\r')) {
                text.pop();
            }
        }
        1 => {
            eol(&mut text, &mut rng);
            eol(&mut text, &mut rng);
        }
        2 => text.extend_from_slice(b">bad\n@@@@\n"),
        _ => {}
    }
    text
}

/// `store` with the first reading in small parts on `threads` threads where the file can be cut
/// (`x_query_scan_par::small_lengths`), so that the second reading reads those parts.
#[cfg(all(unix, feature = "parallel"))]
fn store_in_parts(path: &Path, range: Option<&SequenceRange>, threads: usize) -> (Store, QueryEnd) {
    let file = File::open(path).expect("open the test file");
    let opening = Opening::with_mode(true, path, &file, false, true, NO_RULE);
    assert!(opening.acts());
    let mut source = FastaInputSource::from_file(file, config());
    assert!(!source.stream_is_empty());
    Store::read_with(opening, &mut source, config(), range, |file, start| {
        super::super::x_query_scan_par::small_lengths(file, start, config(), threads)
    })
    .expect("the first reading")
}

/// The batches of the second reading in `sizes` records (cycled): every record equals the
/// first's.
fn read_in_batches(
    path: &Path,
    old: &QueryRecords,
    range: Option<&SequenceRange>,
    sizes: &[usize],
) {
    let (store, end) = store(path, range);
    let pool = crate::utils::threading::SearchPool::x_serial();
    check_batches(&store, &end, path, old, range, sizes, &pool);
}

/// The batches of `store`'s second reading in `sizes` records (cycled), taken with `pool`:
/// every record equals the first's.
fn check_batches(
    store: &Store,
    end: &QueryEnd,
    path: &Path,
    old: &QueryRecords,
    range: Option<&SequenceRange>,
    sizes: &[usize],
    pool: &crate::utils::threading::SearchPool<'_>,
) {
    assert_eq!(end_text(end), end_text(&old.end), "{path:?}");
    let cut = range.map(|range| cut_queries(&old.records, range));
    let searched_old: &[FastaRecord] = cut.as_ref().map_or(&old.records, |cut| &cut.records);
    let queries = Queries::Streamed(store);
    let input = store.query_input();
    let count = input.input_lengths.len();
    assert_eq!(count, old.records.len());
    assert_eq!(queries.input_count(), count);
    let lengths: Vec<usize> = searched_old.iter().map(|r| r.seq().len()).collect();
    assert_eq!(queries.searched_lengths(), lengths);
    let mut start = 0;
    let mut turn = 0;
    while start < count {
        let end = (start + sizes[turn % sizes.len()]).min(count);
        turn += 1;
        let searched = input.searched_in(start..end);
        let batch = queries
            .batch_in(start..end, searched.clone(), pool)
            .unwrap();
        assert_eq!(batch.input(), &old.records[start..end], "{path:?}");
        assert_eq!(
            batch.searched(),
            &searched_old[searched.clone()],
            "{path:?}"
        );
        drop(batch);
        queries.done(end);
        start = end;
    }
}

/// Texts that exercise the reader's states between records (line ends, buffer windows, long
/// records, warnings, an error at the end): the second reading reads the records of the
/// first, in batches of every size, with and without `-query_loc`.
#[test]
#[cfg_attr(
    target_os = "wasi",
    ignore = "File::try_clone is unsupported on WASI, so LOSAT_X_QSTREAM never reads a file batch by batch there"
)]
fn the_second_reading_reads_the_records_of_the_first() {
    let range =
        crate::blastinput::seq_range::parse_sequence_range("50-1200", RangeRole::Query).unwrap();
    let (mut warned, mut failed, mut long) = (0, 0, 0);
    for seed in 1..=40u64 {
        let file = TempFile::new("records", &text(seed));
        let old = old(&file.0);
        warned += usize::from(old.records.iter().any(|r| !r.warnings.is_empty()));
        failed += usize::from(matches!(old.end, QueryEnd::Error { .. }));
        long += usize::from(old.records.iter().any(|r| r.seq().len() > 1 << 20));
        for sizes in [&[1usize][..], &[3, 1, 7], &[usize::MAX]] {
            read_in_batches(&file.0, &old, None, sizes);
            read_in_batches(&file.0, &old, Some(&range), sizes);
        }
    }
    assert!(
        warned >= 10 && failed >= 5 && long >= 5,
        "{warned} {failed} {long}"
    );
}

/// The second reading in parts (the first reading in parts: a file without `\r`, whose reading
/// ends at its end) reads the records of the first, in batches of every size, on a pool of 2 and
/// 4 threads, with and without `-query_loc`: their local IDs, titles, letters and messages (the
/// line numbers of the file); a file with a `\r` is read by one reader.
#[cfg(all(unix, feature = "parallel"))]
#[test]
fn the_second_reading_in_parts_reads_the_records_of_the_first() {
    let range =
        crate::blastinput::seq_range::parse_sequence_range("50-1200", RangeRole::Query).unwrap();
    let (mut in_parts, mut warned) = (0, 0);
    for seed in 1..=60u64 {
        let text = text(seed);
        let file = TempFile::new("parts", &text);
        let old = old(&file.0);
        for threads in [2, 4] {
            for ranged in [None, Some(&range)] {
                for sizes in [&[1usize][..], &[3, 1, 7], &[usize::MAX]] {
                    let (store, end) = store_in_parts(&file.0, ranged, threads);
                    let parts = store.acting().scan_par;
                    assert_eq!(
                        parts,
                        matches!(store.acting().lock().second, Second::Parts { .. })
                    );
                    if text.contains(&b'\r') {
                        assert!(!parts, "seed {seed}: parts with a `\\r`");
                    }
                    in_parts += usize::from(parts);
                    warned += usize::from(
                        parts && old.records.iter().skip(1).any(|r| !r.warnings.is_empty()),
                    );
                    crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
                        check_batches(&store, &end, &file.0, &old, ranged, sizes, pool);
                        Ok(())
                    })
                    .unwrap();
                }
            }
        }
    }
    assert!(in_parts >= 100 && warned >= 20, "{in_parts} {warned}");
}

/// Records of about 1,100 bytes each (a title and ten lines of 100 letters), `count` of them,
/// the letters of record `index` from `letter`.
#[cfg(all(unix, feature = "parallel"))]
fn even_records(count: usize, letter: impl Fn(usize) -> u8) -> Vec<u8> {
    let mut bytes = Vec::new();
    for index in 0..count {
        bytes.extend_from_slice(format!(">r{index} a title\n").as_bytes());
        for _ in 0..10 {
            bytes.extend(std::iter::repeat_n(letter(index), 100));
            bytes.push(b'\n');
        }
    }
    bytes
}

/// The batches of one record of the second reading in parts of the file `path` changed to
/// `changed` after the first reading: the records read up to the first error (in order) and
/// the error, where the reading failed. A second reading of the failed batch is the same error.
#[cfg(all(unix, feature = "parallel"))]
fn read_changed(path: &Path, changed: &[u8]) -> (Vec<FastaRecord>, Option<(usize, String)>) {
    let (store, _) = store_in_parts(path, None, 2);
    assert!(store.acting().scan_par, "the first reading in parts");
    std::fs::write(path, changed).unwrap();
    let queries = Queries::Streamed(&store);
    let count = queries.input_count();
    crate::utils::threading::with_search_pool(2, "blastn", |pool| {
        let mut records = Vec::new();
        for index in 0..count {
            match queries.batch_in(index..index + 1, index..index + 1, pool) {
                Ok(batch) => {
                    records.extend_from_slice(batch.input());
                    drop(batch);
                    queries.done(index + 1);
                }
                Err(error) => {
                    let again = queries
                        .batch_in(index..index + 1, index..index + 1, pool)
                        .err()
                        .expect("the error again");
                    assert_eq!(again.to_string(), error.to_string());
                    return Ok((records, Some((index, error.to_string()))));
                }
            }
        }
        Ok((records, None))
    })
    .unwrap()
}

/// The second reading in parts of a file changed after the first reading: a record of another
/// length is LOSAT's error at the batch that needs it (the parts are read whole, ahead of the
/// loop's batch: the error waits for that batch, and the records before it are read); a change
/// that keeps the bytes of every part and its lines is read as the new content; a change that
/// moves a part is an error; letters added to the last record are an error, a record added after
/// it is not read.
#[cfg(all(unix, feature = "parallel"))]
#[test]
fn a_file_changed_between_the_readings_in_parts() {
    let bytes = even_records(40, |_| b'A');
    let file = TempFile::new("changed_parts", &bytes);
    let read = |path: &Path| old(path).records;
    let old_records = read(&file.0);
    // Record 30 one letter shorter.
    let mut shorter = bytes.clone();
    let at = shorter
        .windows(5)
        .position(|window| window == b">r30 ")
        .unwrap();
    let line = at + b">r30 a title\n".len();
    shorter.remove(line);
    let (records, error) = read_changed(&file.0, &shorter);
    let (index, error) = error.expect("an error");
    assert_eq!(index, 30);
    assert!(
        error.contains("changed while LOSAT read it (record 31: 999 letters where the first reading read 1000)"),
        "{error}"
    );
    assert_eq!(records, old_records[..30]);
    // Record 25's letters changed, every byte in its place: the new records.
    std::fs::write(&file.0, &bytes).unwrap();
    let same = even_records(40, |index| if index == 25 { b'C' } else { b'A' });
    let (records, error) = read_changed(&file.0, &same);
    assert_eq!(error, None);
    let file_same = TempFile::new("changed_same", &same);
    assert_eq!(records, read(&file_same.0));
    assert_ne!(records, old_records);
    // A longer title in record 5 moves every later byte: its part, which still ends at the old
    // place, reads records of the old lengths (the new title), and the next part no longer
    // starts at a record start: an error there.
    std::fs::write(&file.0, &bytes).unwrap();
    let mut moved = bytes.clone();
    let at = moved
        .windows(5)
        .position(|window| window == b">r5 a")
        .unwrap();
    moved.insert(at + 4, b'x');
    let (records, error) = read_changed(&file.0, &moved);
    let (index, error) = error.expect("an error");
    assert!(
        index > 5 && error.contains("no longer starts at a record start"),
        "{index} {error}"
    );
    assert_eq!(records[..5], old_records[..5]);
    assert_eq!(records[5].title, b"r5 xa title");
    assert_eq!(records[6..], old_records[6..index]);
    // Letters added to the last record: an error at the last batch.
    std::fs::write(&file.0, &bytes).unwrap();
    let mut longer = bytes.clone();
    longer.extend_from_slice(b"ACGT\n");
    let (records, error) = read_changed(&file.0, &longer);
    let (index, error) = error.expect("an error");
    assert_eq!(index, 39);
    assert!(error.contains("record 40: 1004 letters"), "{error}");
    assert_eq!(records, old_records[..39]);
    // A record added after the last: not read.
    std::fs::write(&file.0, &bytes).unwrap();
    let mut added = bytes.clone();
    added.extend_from_slice(b">r40\nACGT\n");
    let (records, error) = read_changed(&file.0, &added);
    assert_eq!(error, None);
    assert_eq!(records, old_records);
}

/// The warnings of every record read again (a search that reads its batches and searches
/// none) are those of the first reading's records.
#[test]
#[cfg_attr(
    target_os = "wasi",
    ignore = "File::try_clone is unsupported on WASI, so LOSAT_X_QSTREAM never reads a file batch by batch there"
)]
fn the_input_warnings_are_those_of_the_records() {
    for seed in [3u64, 6, 12, 18] {
        let file = TempFile::new("warnings", &text(seed));
        let old = old(&file.0);
        let expected: Vec<u8> = old
            .records
            .iter()
            .flat_map(|record| record.warnings.iter().copied())
            .collect();
        let (store, _) = store(&file.0, None);
        let mut written = Vec::new();
        Queries::Streamed(&store)
            .write_input_warnings(&mut written)
            .unwrap();
        assert_eq!(written, expected, "seed {seed}");
        // The same where the second reading reads the parts of the first (the texts that end
        // with a record the reader refuses are read by one reader).
        #[cfg(all(unix, feature = "parallel"))]
        {
            let (store, _) = store_in_parts(&file.0, None, 4);
            assert_eq!(store.acting().scan_par, seed % 4 != 2, "seed {seed}");
            let mut written = Vec::new();
            Queries::Streamed(&store)
                .write_input_warnings(&mut written)
                .unwrap();
            assert_eq!(written, expected, "seed {seed} in parts");
        }
    }
}

/// A record whose length changed between the two readings ends the search with LOSAT's error.
#[test]
#[cfg_attr(
    target_os = "wasi",
    ignore = "File::try_clone is unsupported on WASI, so LOSAT_X_QSTREAM never reads a file batch by batch there"
)]
fn a_file_changed_between_the_readings_is_an_error() {
    let file = TempFile::new("changed", b">a\nACGTACGT\n>b\nACGT\n");
    let (store, _) = store(&file.0, None);
    std::fs::write(&file.0, b">a\nACGTACG\n>b\nACGT\n").unwrap();
    let error = Queries::Streamed(&store)
        .batch(0..2, 0..2)
        .err()
        .expect("an error")
        .to_string();
    assert!(error.contains("changed while LOSAT read it"), "{error}");
}

/// A read error (a record whose length changed) fails every later reading: the loop's next
/// call is LOSAT's error again, not a batch with records missing (the records before the
/// failed one are lost, and the reader is past it).
#[test]
#[cfg_attr(
    target_os = "wasi",
    ignore = "File::try_clone is unsupported on WASI, so LOSAT_X_QSTREAM never reads a file batch by batch there"
)]
fn a_read_error_fails_every_later_reading() {
    let mut bytes = Vec::new();
    for index in 0..6 {
        bytes.extend_from_slice(format!(">r{index}\nACGT\n").as_bytes());
    }
    let file = TempFile::new("read_error", &bytes);
    let (store, _) = store(&file.0, None);
    let mut changed = Vec::new();
    for index in 0..6 {
        let letters = if index == 3 { "ACG" } else { "ACGT" };
        changed.extend_from_slice(format!(">r{index}\n{letters}\n").as_bytes());
    }
    std::fs::write(&file.0, changed).unwrap();
    let queries = Queries::Streamed(&store);
    let first = queries.batch(0..1, 0..1).unwrap();
    assert_eq!(first.input().len(), 1);
    drop(first);
    queries.done(1);
    let failed = queries
        .batch(1..4, 1..4)
        .err()
        .expect("the reading's error")
        .to_string();
    assert!(failed.contains("record 4: 3 letters"), "{failed}");
    let error = queries
        .batch(1..4, 1..4)
        .err()
        .expect("the error again")
        .to_string();
    assert_eq!(error, failed);
}

/// The loop's batch is released when the loop is done with it: no chunk of the window and no
/// `current` holds it, and its records are freed with the loop's last handle.
#[test]
#[cfg_attr(
    target_os = "wasi",
    ignore = "File::try_clone is unsupported on WASI, so LOSAT_X_QSTREAM never reads a file batch by batch there"
)]
fn done_releases_the_loops_batch() {
    let mut bytes = Vec::new();
    for index in 0..6 {
        bytes.extend_from_slice(format!(">r{index}\nACGTACGT\n").as_bytes());
    }
    let file = TempFile::new("release", &bytes);
    let (store, _) = store(&file.0, None);
    let queries = Queries::Streamed(&store);
    let batch = queries.batch(0..3, 0..3).unwrap();
    let held = match &batch {
        BatchRecords::Read { input, .. } => Arc::clone(input),
        BatchRecords::Borrowed { .. } | BatchRecords::Taken(_) => panic!("records read again"),
    };
    // The batch's handle, `current` and this one.
    assert_eq!(Arc::strong_count(&held), 3);
    drop(batch);
    queries.done(3);
    {
        let state = store.acting().lock();
        assert!(state.current.is_none());
        assert!(state.chunks.iter().all(|(start, _)| *start >= 3));
    }
    assert_eq!(Arc::strong_count(&held), 1);
}

/// Where the switch does not act: off, the reports after the loop, standard input, a directory
/// (not a regular file).
#[test]
#[cfg_attr(
    target_os = "wasi",
    ignore = "File::try_clone is unsupported on WASI, so LOSAT_X_QSTREAM never reads a file batch by batch there"
)]
fn inputs_that_cannot_be_read_again() {
    let file = TempFile::new("off", b">a\nACGT\n");
    let handle = File::open(&file.0).unwrap();
    assert_eq!(
        Opening::with_mode(false, &file.0, &handle, false, true, NO_RULE).off,
        Some(Off::NotSet)
    );
    assert_eq!(
        Opening::with_mode(true, &file.0, &handle, false, false, NO_RULE).off,
        Some(Off::Reports)
    );
    assert_eq!(
        Opening::with_mode(true, Path::new("-"), &handle, false, true, NO_RULE).off,
        Some(Off::StandardInput)
    );
    assert_eq!(
        Opening::with_mode(true, &file.0, &handle, true, true, NO_RULE).off,
        Some(Off::StandardInput)
    );
    assert!(Opening::with_mode(true, &file.0, &handle, false, true, NO_RULE).acts());
    // The rule turns the switch off; it reads the file's length. The file is 8 bytes: off for a
    // batch of 5 letters, on for 4.
    for (letters, off) in [(5, Some(Off::Batch)), (4, None)] {
        let rule = Rule {
            batch_letters: letters,
        };
        assert_eq!(
            Opening::with_mode(true, &file.0, &handle, false, true, rule).off,
            off
        );
    }
    #[cfg(unix)]
    {
        let dir = File::open(crate::utils::test_dir::temp_dir()).unwrap();
        assert_eq!(
            Opening::with_mode(
                true,
                &crate::utils::test_dir::temp_dir(),
                &dir,
                false,
                true,
                NO_RULE
            )
            .off,
            Some(Off::NotAFile)
        );
    }
}

/// The first reading's lengths and end are those of `read_queries`.
#[test]
fn the_first_reading_is_that_of_read_queries() {
    for seed in 1..=40u64 {
        let file = TempFile::new("scan", &text(seed));
        let old = old(&file.0);
        let mut source = FastaInputSource::from_file(File::open(&file.0).unwrap(), config());
        assert!(!source.stream_is_empty());
        let (lengths, end) = scan(&mut source);
        let old_lengths: Vec<usize> = old.records.iter().map(|r| r.seq().len()).collect();
        assert_eq!(lengths, old_lengths, "seed {seed}");
        assert_eq!(end_text(&end), end_text(&old.end), "seed {seed}");
    }
}

/// A rule that never turns the switch off.
const NO_RULE: Rule = Rule { batch_letters: 0 };

/// The rule: off for a file shorter than twice the largest batch; the files of v4-D's ws windows
/// (1,775 cDNAs 5,398,904 bytes, 3,492 cDNAs 10,779,147, 5,000 cDNAs 15,440,261, 40 MB
/// 40,005,015, 31,168 cDNAs 96,191,719).
#[test]
fn the_rule_follows_the_batch() {
    const BLASTN: u64 = 999_000;
    const DC: u64 = 4_999_000;
    let rule = |batch_letters| Rule { batch_letters };
    assert_eq!(NO_RULE.off(0), None);
    // The batch: blastn off below 1,998,000 bytes, megablast and dc-megablast below 9,998,000,
    // BATCH_SIZE 4,000 below 8,000.
    assert_eq!(rule(BLASTN).off(1_997_999), Some(Off::Batch));
    assert_eq!(rule(BLASTN).off(1_998_000), None);
    for bytes in [5_398_904, 10_779_147, 15_440_261, 40_005_015, 96_191_719] {
        assert_eq!(rule(BLASTN).off(bytes), None);
    }
    assert_eq!(rule(DC).off(5_398_904), Some(Off::Batch));
    assert_eq!(rule(DC).off(9_997_999), Some(Off::Batch));
    assert_eq!(rule(DC).off(9_998_000), None);
    assert_eq!(rule(DC).off(96_191_719), None);
    assert_eq!(rule(4_000).off(7_999), Some(Off::Batch));
    assert_eq!(rule(4_000).off(8_000), None);
    // The ends saturate.
    assert_eq!(rule(u64::MAX).off(u64::MAX - 1), Some(Off::Batch));
    assert_eq!(rule(u64::MAX).off(u64::MAX), None);
    assert_eq!(rule(0).off(0), None);
}
