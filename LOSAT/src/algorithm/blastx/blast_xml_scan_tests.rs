//! The read of the queries before the search (`scan_queries`): a title the documents cannot
//! write is refused where a report that they write holds it, and a `-subject` XML2 or JSON
//! document never writes the last query before a batch that stops the run.
use super::*;

/// A record of `letters` pseudo-random nucleotides (valid contexts, as `reader_batch`'s tests
/// make them) with the title `title`.
fn record(title: &str, letters: usize) -> FastaRecord {
    let mut state = letters as u64 + 12345;
    let sequence = (0..letters)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 60) as usize & 3]
        })
        .collect();
    FastaRecord {
        internal_id: String::new(),
        title: title.as_bytes().to_vec(),
        sequence,
        lowercase_masks: Vec::new(),
        warnings: Vec::new(),
        source: None,
        around: None,
        skipped: None,
    }
}

/// Scans the reader batches `batches`, then ends the read with an error (`stops`: a reader
/// error) or at the end of the queries.
fn scan(batches: &[Vec<FastaRecord>], stops: bool, holds_back: bool) -> Result<QueryScan> {
    let options = super::super::stage_c_tests::options();
    scan_queries(&options, holds_back, true, &mut |deliver| {
        for batch in batches {
            deliver(batch)?;
        }
        if stops {
            anyhow::bail!("a reader error");
        }
        Ok(())
    })
}

#[test]
fn the_last_query_before_a_stop_is_checked_only_where_its_report_is_written() {
    let plain = || record("plain", 300);
    let bad = || record(", ,", 300);
    let refused =
        |scanned: Result<QueryScan>| scanned.err().is_some_and(|error| error.is::<XmlRefusal>());
    // The last query of the batch before a reader error: held back and never written.
    let held = [vec![plain(), bad()]];
    let found = scan(&held, true, true).unwrap();
    assert_eq!(found.stops_after, Some(2));
    // BLAST XML and `-db` write it.
    assert!(refused(scan(&held, true, false)));
    // At the end of the queries the epilog writes it.
    assert!(refused(scan(&held, false, true)));
    // A later batch that NCBI writes: written when its first query arrives.
    assert!(refused(scan(
        &[vec![plain(), bad()], vec![plain()]],
        true,
        true
    )));
    // A later batch that NCBI fails at its set-up (no valid context, a query without
    // letters) stops the run before it, as a reader error does; its records are not read.
    let empty = [vec![plain(), bad()], vec![record(", ,", 0)]];
    let found = scan(&empty, false, true).unwrap();
    assert_eq!(found.stops_after, Some(2));
    assert!(refused(scan(&empty, false, false)));
    // A title before the last of its batch is written when the next query arrives.
    assert!(refused(scan(&[vec![bad(), plain()]], true, true)));
    // A run that stops at its first batch writes nothing.
    let found = scan(&[], true, true).unwrap();
    assert_eq!(found.stops_after, None);
    // A run whose XML format fails for its matrix name writes no report and no title
    // (`xml_matrix_failure`, probe K4f of line RG): its batches are read, no title is checked.
    let options = super::super::stage_c_tests::options();
    let found = scan_queries(&options, false, false, &mut |deliver| {
        deliver(&[bad(), plain()])?;
        deliver(&[bad()])
    })
    .unwrap();
    assert_eq!(found.stops_after, None);
}

/// The (title, letters) of every record that the search's stream of `file` reads.
fn search_reads(
    settings: &super::super::input::QueryStreamSettings,
    file: crate::blastinput::input_files::InputFile,
    batch_size: usize,
) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut read = Vec::new();
    super::super::input::read_fasta_stream_batches(
        &mut settings.stream(file),
        false,
        false,
        batch_size,
        false,
        None,
        &mut |batch| {
            read.extend(
                batch
                    .iter()
                    .map(|record| (record.title.clone(), record.sequence.clone())),
            );
            Ok(())
        },
        &mut |_| Ok(()),
    )
    .unwrap();
    read
}

// The read before the search takes the open query file, reads it as the search does and gives
// the same open file back at the position it had (no second handle, which the WASI targets
// cannot make): the search then reads every record from there, as from a file nothing read.
#[test]
fn the_read_before_the_search_gives_the_open_file_back_where_it_was() {
    use std::io::{Seek, SeekFrom};
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-blastx-prescan-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("query.fna");
    let mut text = Vec::new();
    for (number, letters) in [(1, 300), (2, 450), (3, 120), (4, 300)] {
        let query = record("", letters);
        text.extend_from_slice(format!(">q{number} query {number}\n").as_bytes());
        for line in query.sequence.chunks(60) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    std::fs::write(&path, &text).unwrap();
    let second = text.windows(4).position(|w| w == b">q2 ").unwrap() as u64;
    let settings = super::super::input::QueryStreamSettings {
        program: super::super::args::ProgramName::BLASTX,
        parse_deflines: false,
        id_scope: None,
        strand: crate::report::archive::bioseq::SeqLocStrand::Both,
        data_loaders: false,
    };
    let options = super::super::stage_c_tests::options();
    let open = |start: u64| {
        let mut file =
            crate::blastinput::input_files::InputFile::from(std::fs::File::open(&path).unwrap());
        file.seek(SeekFrom::Start(start)).unwrap();
        file
    };
    // A batch of every query, a batch per query (the reader's refills), and no batch.
    for batch_size in [100_000, 100, 0] {
        for start in [0, second] {
            let unread = search_reads(&settings, open(start), batch_size.max(1));
            assert_eq!(unread.len(), if start == 0 { 4 } else { 3 });
            let file = open(start);
            #[cfg(unix)]
            let descriptor = std::os::fd::AsRawFd::as_raw_fd(&*file);
            let (scan, mut file) = prescan_queries(
                file, &path, &settings, &options, false, batch_size, None, false, true, None,
            )
            .unwrap();
            assert_eq!(scan.stops_after, None);
            // The same open file, at the position it had.
            #[cfg(unix)]
            assert_eq!(std::os::fd::AsRawFd::as_raw_fd(&*file), descriptor);
            assert_eq!(file.stream_position().unwrap(), start, "{batch_size}");
            assert_eq!(
                search_reads(&settings, file, batch_size.max(1)),
                unread,
                "{batch_size} {start}"
            );
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}
