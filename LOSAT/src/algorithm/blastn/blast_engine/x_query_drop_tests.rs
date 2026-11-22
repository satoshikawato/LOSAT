use super::*;
use crate::blastinput::fasta_reader::{read_queries, ReaderConfig};

/// Nine records (three batches of three), every second one packed, as `search_cli` stores them.
fn store() -> (Vec<FastaRecord>, Store) {
    let records: Vec<FastaRecord> = (0..9)
        .map(|index| {
            // Mostly ACGT (packed), a run of N and a lower-case stretch (kept beside the
            // codes).
            let mut letters: Vec<u8> = b"ACGTTGCAAGCTTCGAGATCCTAGGCATTACG"
                .iter()
                .cycle()
                .skip(index)
                .take(400 + 37 * index)
                .copied()
                .collect();
            letters[100..104].copy_from_slice(b"NNNN");
            letters[200..240].make_ascii_lowercase();
            FastaRecord::new(format!("Query_{index}"), b"title", &letters)
        })
        .collect();
    let mut held = records.clone();
    let mut packed = Vec::new();
    for (index, record) in held
        .iter_mut()
        .enumerate()
        .filter(|(index, _)| index % 2 == 1)
    {
        if let Some(store) = PackedLetters::pack(&record.sequence) {
            record.sequence = Vec::new();
            packed.push((index, store));
        }
    }
    assert!(!packed.is_empty(), "no record packed");
    (records, Store::new(held, packed))
}

/// The records of `text` as `read_queries` reads them, and how the reading ended.
fn read_all(text: &[u8], config: ReaderConfig) -> (Vec<FastaRecord>, String) {
    let records = read_queries(&mut FastaInputSource::from_bytes(text, config));
    (records.records, format!("{:?}", records.end))
}

/// The records of `text` read into parts of `part_bytes` and taken back out of their store in
/// batches of `batch` records (by a pool of `threads` threads), how the reading ended and the
/// number of parts.
fn through_parts(
    text: &[u8],
    config: ReaderConfig,
    part_bytes: usize,
    batch: usize,
    threads: usize,
) -> (Vec<FastaRecord>, String, usize) {
    let mut source = FastaInputSource::from_bytes(text, config);
    let (parts, end) = read_parts_of(&mut source, config.id_prefix, part_bytes);
    let count = parts.len();
    let store = Store::from_parts(parts, config.id_prefix, false);
    let total = store.lengths().len();
    let take = |pool: &SearchPool<'_>| -> anyhow::Result<Vec<FastaRecord>> {
        let mut records = Vec::new();
        for from in (0..total).step_by(batch.max(1)) {
            records.extend(store.take(from..(from + batch).min(total), pool)?);
        }
        Ok(records)
    };
    let records = if threads > 1 {
        crate::utils::threading::with_search_pool(threads, "blastn", take).expect("take")
    } else {
        take(&SearchPool::x_serial()).expect("take")
    };
    (records, format!("{end:?}"), count)
}

/// FASTA of `records` records: titles, `N` runs, lower case, a record without letters, one
/// without an id, ids of several forms (for `-parse_deflines`), then `tail`.
fn text(records: usize, tail: &[u8]) -> Vec<u8> {
    let mut text = Vec::new();
    for index in 0..records {
        match index % 7 {
            0 => text.extend_from_slice(b">\n"),
            1 => text.extend_from_slice(format!(">gi|{index}|ref|NC_{index}.1| t\n").as_bytes()),
            2 => text.extend_from_slice(format!(">lcl|q{index}\n").as_bytes()),
            _ => text.extend_from_slice(format!(">read{index} a title {index}\n").as_bytes()),
        }
        if index % 11 == 5 {
            continue;
        }
        let letters: Vec<u8> = b"ACGTTGCAAGCTTCGAGATCCTAGGCATTACG"
            .iter()
            .cycle()
            .skip(index)
            .take(20 + (index * 37) % 300)
            .copied()
            .collect();
        for line in letters.chunks(60) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
        if index % 9 == 4 {
            text.extend_from_slice(b"NNNNNNNNNNacgtacgtacgtnnnnACGT\n");
        }
    }
    text.extend_from_slice(tail);
    text
}

// The records read into parts and taken back out of their store are the reader's, with how the
// reading ended (the end, a blank-line `eEOF`, a reader error) and the reader's messages, also
// with the ids `-parse_deflines` reads, across parts and batches of any size, by one thread and
// the pool.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn records_read_into_parts_are_the_readers_records() {
    let plain = ReaderConfig::query("BLASTN", false, false);
    let parsed = plain.with_parse_deflines(true);
    // A record that the reader does not pack (a run on almost every letter) and that is long
    // enough to keep its own vector.
    let mut dense = b">dense runs\n".to_vec();
    for line in b"ARCYGKTMNSAW".repeat(6_000).chunks(70) {
        dense.extend_from_slice(line);
        dense.push(b'\n');
    }
    dense.extend_from_slice(b">after\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
    for (tail, ends) in [
        (dense.as_slice(), "Input"),
        (b"".as_slice(), "Input"),
        (b"\n; a comment\n\n".as_slice(), ""),
        (b">hyphen\nACGT-ACGT\n".as_slice(), "Input"),
        (b">bad\n%%%%%%%%%%%%%%%%\n".as_slice(), "Error"),
    ] {
        for config in [plain, parsed] {
            let text = text(400, tail);
            let (expected, end) = read_all(&text, config);
            assert!(end.starts_with(ends), "{end}");
            for part_bytes in [1, 700, 5_000, 1 << 20] {
                for (batch, threads) in [(1, 1), (37, 1), (64, 4), (500, 3)] {
                    let (records, got_end, parts) =
                        through_parts(&text, config, part_bytes, batch, threads);
                    assert_eq!(got_end, end);
                    assert_eq!(records, expected, "part_bytes {part_bytes} batch {batch}");
                    if tail == dense.as_slice() {
                        let long = &records[records.len() - 2];
                        assert!(long.sequence.len() >= OWN_LETTERS);
                        assert!(PackedLetters::pack(&long.sequence).is_none());
                    }
                    assert!(part_bytes >= 5_000 || parts > 5, "{parts} parts");
                }
            }
        }
    }
}

// The loop's batches taken out of the store are the records read, the packed letters made
// again (one thread and the pool).
#[test]
fn the_loop_takes_its_batches_with_their_letters() {
    let (records, held) = store();
    let serial = SearchPool::x_serial();
    for range in [0..3, 3..6, 6..9] {
        let batch = held.take(range.clone(), &serial).expect("take");
        assert_eq!(&batch[..], &records[range]);
    }
    assert!(held.take(9..10, &serial).is_err());
}
