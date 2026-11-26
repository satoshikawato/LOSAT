use super::*;
use crate::algorithm::blastx::args::ProgramName;
use crate::algorithm::blastx::input::{scope_records, shown_ids};
use crate::blastinput::input_files::InputFile;
use crate::blastinput::seq_range::RecordPlacement;
use std::path::Path;

/// The records of `data` (a query, or with `protein` a subject), read as the command line reads
/// them, the whole input as one batch (`BatchCutter`).
pub(crate) fn parse_fasta(data: &[u8], protein: bool, lowercase: bool) -> Result<Vec<FastaRecord>> {
    let config = if protein {
        ReaderConfig::subject("BLASTX", true, false)
    } else {
        ReaderConfig::query("BLASTX", false, false)
    };
    let mut source = FastaInputSource::from_bytes(data, config);
    let mut kept = Vec::new();
    if source.end() {
        return Ok(kept);
    }
    let mut cutter = BatchCutter::new(CutSettings {
        range: None,
        keep_whole: false,
        around: 0,
        keep_skipped: false,
        batch_size: usize::MAX,
        count: None,
        blank_batch: true,
    });
    while !source.end() {
        match source.next_sequence_upper(lowercase, &mut |_| Ok(())) {
            Ok((record, masks)) => {
                let record = from_upper(record, masks);
                kept.push(record.clone());
                cutter.add(record, &mut |_| Ok(()))?;
            }
            Err(ReadError::Parse {
                code: ParseErrorCode::Eof,
                ..
            }) => break,
            Err(error) => return Err(error.into_app_error()),
        }
    }
    cutter.finish(&mut |_| Ok(()))?;
    Ok(kept)
}

/// The records of the file `path`.
pub(crate) fn read_fasta(path: &Path, protein: bool, lowercase: bool) -> Result<Vec<FastaRecord>> {
    let data = anyhow::Context::with_context(std::fs::read(path), || {
        format!("cannot read BLASTX input {}", path.display())
    })?;
    parse_fasta(&data, protein, lowercase)
}

/// The directory of this test's files (removed at its end).
fn dir(tag: &str) -> std::path::PathBuf {
    crate::utils::test_dir::temp_dir().join(format!(
        "losat-blastx-input-source-{tag}-{}",
        crate::utils::test_dir::process_id()
    ))
}

/// A file of this test with the bytes `text`, open for reading.
fn file(dir: &std::path::Path, name: &str, text: &[u8]) -> std::fs::File {
    std::fs::create_dir_all(dir).unwrap();
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    std::fs::File::open(&path).unwrap()
}

fn settings(parse_deflines: bool, data_loaders: bool) -> QueryStreamSettings {
    QueryStreamSettings {
        program: ProgramName::BLASTX,
        parse_deflines,
        id_scope: parse_deflines.then(Default::default),
        strand: crate::report::archive::bioseq::SeqLocStrand::Both,
        data_loaders,
        batch_count: None,
    }
}

/// The batches of `stream`, each given as the records it holds (cloned).
fn read_all(
    stream: &mut QueryStream<InputFile>,
    batch_size: usize,
) -> Result<Vec<Vec<FastaRecord>>> {
    let mut batches = Vec::new();
    stream.read_batches(
        false,
        batch_size,
        None,
        &mut |records| {
            batches.push(records.to_vec());
            Ok(())
        },
        &mut |_| Ok(()),
    )?;
    Ok(batches)
}

/// A delivered batch: the ids, letters and placements of its records.
type Batch = Vec<(String, Vec<u8>, RecordPlacement)>;

fn summary(records: &[FastaRecord]) -> Batch {
    records
        .iter()
        .map(|record| {
            (
                record.internal_id.clone(),
                record.sequence.clone(),
                record.placement(),
            )
        })
        .collect()
}

/// The reader batches of the nucleotide query file `data` (`-lcase_masking`), read as the command
/// line reads it: the records of each batch, cloned.
fn batches_of(
    data: &[u8],
    range: Option<&SequenceRange>,
    batch_size: usize,
    keep_skipped: bool,
) -> Vec<Vec<FastaRecord>> {
    let mut source =
        FastaInputSource::from_bytes(data, ReaderConfig::query("BLASTX", false, false));
    let mut batches = Vec::new();
    read_source_batches(
        &mut source,
        true,
        CutSettings {
            range,
            keep_whole: false,
            around: 0,
            keep_skipped,
            batch_size,
            count: None,
            blank_batch: true,
        },
        &mut |records| {
            batches.push(records.to_vec());
            Ok(())
        },
        &mut |_| Ok(()),
        ReadError::into_app_error,
    )
    .unwrap();
    batches
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:144-168
// A batch counts the whole length of each searched record; a skipped record (its interval
// starts more than one letter past its end) counts nothing, is not delivered and keeps its
// number; the letters delivered are the interval's.
#[test]
fn a_range_cuts_the_records_and_skips_the_short_ones() {
    let range = SequenceRange { from: 4, to: 7 };
    let data = b">a\nACGTACGTAC\n>b\nACG\n>c\nTTTTGGGG\n>d\nACGT\n>e\nACGTA\n";
    let batches = batches_of(data, Some(&range), 15, true);
    let place = |offset, length| RecordPlacement { offset, length };
    assert_eq!(
        batches
            .iter()
            .map(|batch| summary(batch))
            .collect::<Vec<_>>(),
        vec![
            vec![
                ("Query_1".to_string(), b"ACGT".to_vec(), place(4, 10)),
                ("Query_3".to_string(), b"GGGG".to_vec(), place(4, 8)),
            ],
            // `d` starts just past its end: an interval without letters; `e` keeps one letter.
            vec![
                ("Query_4".to_string(), Vec::new(), place(4, 4)),
                ("Query_5".to_string(), b"A".to_vec(), place(4, 5)),
            ],
        ]
    );
    // The skipped `b` stays in its batch's scope, on the batch's first record (review 2 of
    // line DBID loadord, H2; blast_fasta_input.cpp:395 adds it to the scope first); the second
    // batch skipped none.
    let ids = |records: &[FastaRecord]| -> Vec<String> {
        shown_ids(records).map(|(id, _)| id.to_owned()).collect()
    };
    assert_eq!(ids(&batches[0]), ["Query_1", "Query_3", "Query_2"]);
    assert_eq!(ids(&batches[1]), ["Query_4", "Query_5"]);
    // A run whose prefetch cannot stop keeps none of them (review 5 of line DBID loadord, L2);
    // the numbers stay.
    let unkept = batches_of(data, Some(&range), 15, false);
    assert_eq!(ids(&unkept[0]), ["Query_1", "Query_3"]);
    assert_eq!(ids(&unkept[1]), ["Query_4", "Query_5"]);
    // The scopes of the queries take the records in the order read, the skipped `b` with its
    // whole letters between `a` and `c` (`scope_records`; line WQA D10 round 6).
    let order: Vec<(usize, Option<usize>)> = scope_records(&batches[0])
        .into_iter()
        .map(|(_, length, place)| (length, place))
        .collect();
    assert_eq!(order, [(10, Some(0)), (3, None), (8, Some(1))]);
}

// The batches of a read without a range, by their size: a batch ends at the record that brings
// it to `batch_size` letters, the rest is delivered at the end, and the numbers go on.
#[test]
fn the_batches_end_at_the_record_that_reaches_the_batch_size() {
    let data =
        b">a\nACGTACGTAC\n>b\nACG\n>c\nTTTTGGGG\n>d\nACGT\n>e\nACGTA\n>f\nAC\n>g\nACGTACGT\n";
    let ids = |batch_size| -> Vec<Vec<String>> {
        batches_of(data, None, batch_size, false)
            .iter()
            .map(|batch| batch.iter().map(|r| r.internal_id.clone()).collect())
            .collect()
    };
    assert_eq!(ids(usize::MAX).len(), 1);
    assert_eq!(ids(10).len(), 4);
    assert_eq!(ids(10)[3], ["Query_7"]);
    assert_eq!(ids(1).len(), 7);
}

// A batch of skipped records alone is delivered empty (NCBI's query factory refuses it), also
// at the end of the file.
#[test]
fn a_last_batch_of_skipped_records_is_delivered_empty() {
    let range = SequenceRange { from: 4, to: 7 };
    let delivered = batches_of(b">a\nACGTACGTAC\n>b\nACG\n>c\nA\n", Some(&range), 10, false);
    assert_eq!(delivered.len(), 2);
    assert_eq!(delivered[0].len(), 1);
    assert!(delivered[1].is_empty());
    let delivered = batches_of(b">a\nACG\n", Some(&range), 10, false);
    assert_eq!(delivered.len(), 1);
    assert!(delivered[0].is_empty());
}

// Blank lines alone, an input without a record, are one batch without a query; an input without
// a line gives no batch; ABI v1's read gives none for either.
#[test]
fn an_input_of_no_record_gives_one_empty_batch_or_none() {
    assert_eq!(batches_of(b"\n\n", None, 10, false).len(), 1);
    assert!(batches_of(b"\n\n", None, 10, false)[0].is_empty());
    assert!(batches_of(b"", None, 10, false).is_empty());
    let mut web = Vec::new();
    read_web_queries(
        b"\n\n",
        false,
        10,
        &mut |records| {
            web.push(records.len());
            Ok(())
        },
        &mut |_| Ok(()),
    )
    .unwrap();
    assert!(web.is_empty());
}

// The read delivers each batch and drops it: it holds one reader batch, not the file. Under
// `-parse_deflines` (without the search's ids) a record is the only owner of its Seq-ids, so a
// weak reference to them shows whether the record is still held.
#[test]
fn a_read_frees_each_batch_once_delivered() {
    use crate::blastinput::fasta_reader::believed::BelievedIds;
    use std::sync::{Arc, Weak};
    let mut source = FastaInputSource::from_bytes(
        b">lcl|q1\nACGT\n>lcl|q2\nACGT\n>lcl|q3\nACGT\n",
        ReaderConfig::query("BLASTX", false, false).with_parse_deflines(true),
    );
    let mut ids: Vec<Weak<BelievedIds>> = Vec::new();
    // At each delivery, the records of the earlier batches still held.
    let mut held_before = Vec::new();
    read_source_batches(
        &mut source,
        false,
        CutSettings {
            range: None,
            keep_whole: false,
            around: 0,
            keep_skipped: false,
            batch_size: 1,
            count: None,
            blank_batch: true,
        },
        &mut |records| {
            held_before.push(ids.iter().filter(|id| id.upgrade().is_some()).count());
            for record in records {
                let believed = record
                    .source
                    .as_ref()
                    .and_then(|source| source.believed.as_ref())
                    .expect("parsed Seq-ids");
                ids.push(Arc::downgrade(believed));
            }
            Ok(())
        },
        &mut |_| Ok(()),
        ReadError::into_app_error,
    )
    .unwrap();
    assert_eq!(held_before, vec![0, 0, 0]);
    assert_eq!(ids.iter().filter(|id| id.upgrade().is_some()).count(), 0);
}

// A record whose letter before the interval is lower case keeps its whole letters after the
// cut (that letter's run, which a split query's chunks restrict), whether the run ends there
// or goes on into the interval; a record with runs only elsewhere, or none, keeps only its
// interval.
#[test]
fn a_cut_record_keeps_its_whole_letters_for_a_lower_case_letter_before_it() {
    let range = SequenceRange { from: 4, to: 7 };
    let data = b">a\nACGtACGTAC\n>b\nAcgtaCGTAC\n>c\naCGTaCGTac\n>d\nACGTACGTAC\n";
    let kept = batches_of(data, Some(&range), 100, false).concat();
    assert_eq!(kept[0].sequence, b"ACGT");
    assert!(kept[0].lowercase_masks.is_empty());
    assert_eq!(kept[0].archive_letters(), b"ACGtACGTAC");
    assert_eq!(kept[0].archive_lowercase_runs(), vec![(3, 3)]);
    assert_eq!(kept[1].lowercase_masks, vec![(0, 0)]);
    assert_eq!(kept[1].archive_letters(), b"AcgtaCGTAC");
    assert_eq!(kept[2].lowercase_masks, vec![(0, 0)]);
    assert_eq!(kept[2].archive_letters(), b"ACGT");
    assert_eq!(kept[3].archive_letters(), b"ACGT");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1041-1046
#[test]
fn a_cut_clips_and_moves_the_lower_case_masks() {
    let mut record = parse_fasta(b">q\nacgTTTTaaaCCgg\n", false, true)
        .unwrap()
        .remove(0);
    assert_eq!(record.lowercase_masks, vec![(0, 2), (7, 9), (12, 13)]);
    let mut kept = record.clone();
    record.cut_in_place(2, 9, false);
    assert_eq!(record.sequence, b"GTTTTAA");
    assert_eq!(record.lowercase_masks, vec![(0, 0), (5, 6)]);
    assert_eq!(
        record.placement(),
        RecordPlacement {
            offset: 2,
            length: 14
        }
    );
    // Only a record that keeps its whole letters has them (and its own runs) after the cut.
    assert_eq!(record.archive_letters(), b"GTTTTAA");
    assert_eq!(record.archive_lowercase_runs(), vec![(0, 0), (5, 6)]);
    kept.cut_in_place(2, 9, true);
    assert_eq!(kept.sequence, record.sequence);
    assert_eq!(kept.lowercase_masks, record.lowercase_masks);
    assert_eq!(kept.placement(), record.placement());
    assert_eq!(kept.archive_letters(), b"acgTTTTaaaCCgg");
    assert_eq!(
        kept.archive_lowercase_runs(),
        vec![(0, 2), (7, 9), (12, 13)]
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1082-1090
// Both `next_sequence_upper` + `from_upper` and `next_sequence` + `from_shared_ref` produce
// identical BLASTX records (upper-case letters, identical mask intervals over `>?` gaps, U->T,
// protein `*`, trailing gaps, and empty masks without `-lcase_masking`).
#[test]
fn the_adapter_makes_upper_case_letters_and_masks() {
    let read_shared = |data: &[u8], protein: bool| {
        let cfg = if protein {
            ReaderConfig::subject("BLASTX", true, false)
        } else {
            ReaderConfig::query("BLASTX", false, false)
        };
        let mut source = FastaInputSource::from_bytes(data, cfg);
        source.next_sequence(&mut |_| Ok(())).unwrap()
    };
    let read_upper = |data: &[u8], protein: bool, keep_masks: bool| {
        let cfg = if protein {
            ReaderConfig::subject("BLASTX", true, false)
        } else {
            ReaderConfig::query("BLASTX", false, false)
        };
        let mut source = FastaInputSource::from_bytes(data, cfg);
        let (record, masks) = source
            .next_sequence_upper(keep_masks, &mut |_| Ok(()))
            .unwrap();
        from_upper(record, masks)
    };
    for (data, protein, expected_seq, expected_masks) in [
        // The gap `>?5` lies in the run that `gt` continues: one interval over the gap.
        (
            &b">r0 t\nac\n>?5\ngtAC\n"[..],
            false,
            &b"ACNNNNNGTAC"[..],
            vec![(0, 8)],
        ),
        // Trailing gap in a masked run.
        (b">r1\nac\n>?5\n", false, b"ACNNNNN", vec![(0, 6)]),
        // Gap between masked and unmasked, and between unmasked and masked.
        (b">r2\nac\n>?5\nAC\n", false, b"ACNNNNNAC", vec![(0, 6)]),
        (b">r3\nAC\n>?5\nac\n", false, b"ACNNNNNAC", vec![(7, 8)]),
        // U->T under masking.
        (b">r4\nuuAC\n", false, b"TTAC", vec![(0, 1)]),
        // Two runs, one at each end.
        (
            b">q\nacgTTTTaaaCCgg\n",
            false,
            b"ACGTTTTAAACCGG",
            vec![(0, 2), (7, 9), (12, 13)],
        ),
        // A record without lower case.
        (b">q\nACGT\n", false, b"ACGT", vec![]),
        // Protein with `*` splitting masked runs.
        (b">p\nac*gt\n", true, b"AC*GT", vec![(0, 1), (3, 4)]),
    ] {
        let shared = read_shared(data, protein);
        let masked_upper = read_upper(data, protein, true);
        let masked_ref = from_shared_ref(&shared, true).unwrap();
        assert_eq!(masked_upper.sequence, expected_seq);
        assert_eq!(masked_upper.lowercase_masks, expected_masks);
        assert_eq!(masked_ref.sequence, masked_upper.sequence);
        assert_eq!(masked_ref.lowercase_masks, masked_upper.lowercase_masks);
        assert_eq!(masked_ref.internal_id, masked_upper.internal_id);
        assert_eq!(masked_ref.title, masked_upper.title);

        let unmasked_upper = read_upper(data, protein, false);
        let unmasked_ref = from_shared_ref(&shared, false).unwrap();
        assert_eq!(unmasked_upper.sequence, expected_seq);
        assert!(unmasked_upper.lowercase_masks.is_empty());
        assert_eq!(unmasked_ref.sequence, unmasked_upper.sequence);
        assert!(unmasked_ref.lowercase_masks.is_empty());
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:130-133
// With the data loaders on, the stream refuses a first line that NCBI tries as a Seq-id with the
// reader's text, naming the run's program, also under `-parse_deflines`; with them off the line
// is data.
#[test]
fn a_first_line_ncbi_tries_as_a_seq_id_is_refused_with_the_data_loaders_on() {
    let dir = dir("seqid");
    for parse_deflines in [false, true] {
        let mut stream = settings(parse_deflines, true)
            .stream(file(&dir, "lcl", b"lcl|abc\n>q1 x\nACGTACGT\n"), false);
        let error = read_all(&mut stream, 100).unwrap_err().to_string();
        assert!(
            error.contains("is not a defline and may be a sequence identifier")
                && error.contains("not supported by LOSAT's BLASTX"),
            "{error}"
        );
    }
    let mut rps = settings(false, true);
    rps.program = ProgramName::RPSTBLASTN;
    let error = read_all(
        &mut rps.stream(file(&dir, "rps", b"12345\n>q1 x\nAC\n"), false),
        100,
    )
    .unwrap_err()
    .to_string();
    assert!(
        error.contains("not supported by LOSAT's RPSTBLASTN"),
        "{error}"
    );
    // With the loaders off the first line is data.
    let batches = read_all(
        &mut settings(false, false).stream(file(&dir, "off", b"ACGTACGTAC\n>q1 x\nACGT\n"), false),
        100,
    )
    .unwrap();
    assert_eq!(batches.concat().len(), 2);
    std::fs::remove_dir_all(dir).unwrap();
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:383-398
// Under -parse_deflines a gap line before any defline stops the read at its line, whatever its
// size or modifiers, with the reader's text as the application writes it; without it, the gap line
// starts a record without a defline.
#[test]
fn a_gap_line_before_any_defline_stops_a_parse_deflines_read() {
    let dir = dir("nodefline");
    for (data, line) in [
        (&b">?5\nACGT\n>r1\nACGT\n"[..], 1),
        (b">?abc\nACGT\n", 1),
        (b"  >?5\r\nACGT\r\n", 1),
        (b";c\n\n#x\n>?5\nACGT\n", 4),
    ] {
        let mut stream = settings(true, false).stream(file(&dir, "gap", data), false);
        let error = read_all(&mut stream, 100).unwrap_err().to_string();
        assert!(
            error.contains(&format!(
                "LOSAT query error: CFastaReader: Input doesn't start with a defline or comment around line {line}"
            )),
            "{}: {error}",
            String::from_utf8_lossy(data)
        );
    }
    let batches = read_all(
        &mut settings(false, false).stream(file(&dir, "plain", b">?5\nACGT\n>r1\nACGT\n"), false),
        100,
    )
    .unwrap();
    let records = batches.concat();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].sequence, b"NNNNNACGT");
    std::fs::remove_dir_all(dir).unwrap();
}

// The read before the search reads as the search does (under `-parse_deflines` a first line
// that is not a defline is the reader's error) and adds its ids to a copy of the search's: a
// subject with the id of a query that only that read has read is not refused, and that read
// stops at a query with the id of a subject the search has read, as the search does.
#[test]
fn a_read_before_the_search_reads_as_the_search_and_keeps_the_search_ids() {
    let dir = dir("settings");
    let settings = settings(true, false);
    let error = read_all(
        &mut settings.stream_before_search(file(&dir, "nodef", b"ACGT\n>q1 x\nACGT\n")),
        100,
    )
    .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("Input doesn't start with a defline or comment around line 1"),
        "{error}"
    );
    let queries = read_all(
        &mut settings.stream_before_search(file(&dir, "query", b">q1 x\nACGTACGT\n")),
        100,
    )
    .unwrap();
    assert_eq!(queries.concat().len(), 1);
    let mut subjects = FastaInputSource::from_file(
        file(&dir, "subject", b">q1 y\nMKVL\n"),
        ReaderConfig::subject("BLASTX", true, false).with_parse_deflines(true),
    );
    subjects.set_id_scope(settings.id_scope.clone());
    fasta_reader::read_subjects(&mut subjects, None, &mut |_| Ok(())).unwrap();
    let error = read_all(
        &mut settings.stream_before_search(file(&dir, "again", b">q1 x\nACGTACGT\n")),
        100,
    )
    .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("a query and a subject record have the Seq-id lcl|q1"),
        "{error}"
    );
    std::fs::remove_dir_all(dir).unwrap();
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:717-721 (`-` is `cin`) and
// c++/src/corelib/stream_utils.cpp:235-236 (a refill asks for `in_avail()` bytes), cited at
// `fasta_reader::stream::InputBackend::of`.
/// BLASTX's and RPSTBLASTN's query stream reads standard input as NCBI reads `cin`, by the
/// readers' one rule (`InputBackend::of`): the tail that the same bytes lose in a file named by
/// its path is kept (NCBI BLAST+ 2.17.0, blastx and rpstblastn: `Length=107` from a pipe and
/// from a redirected file, 85 from the path; fixture rows `blastx.stdin_q_lost_g55_pipe`,
/// `blastx.stdin_q_lost_g55_file`, `db.rpstbn.g55.*`).
#[test]
fn standard_input_keeps_the_tail_that_a_named_file_loses() {
    let text = b">first\r\nTC-TTGGCTCAATCCTAGGTGGGCATGTTTCCTAATGCCC\rTT-TTTAACGTGAGGGTTCGCGTTTTTATCCCACCTAGC\r\r\nACGTACGT\n-ACGTACGTACGTACGTACGTAC";
    let dir = dir("stdin");
    let path = {
        std::fs::create_dir_all(&dir).unwrap();
        dir.join("input")
    };
    std::fs::write(&path, text).unwrap();
    let lengths = |input: InputFile| {
        let mut stream = settings(false, false).stream(input, false);
        read_all(&mut stream, usize::MAX)
            .unwrap()
            .concat()
            .iter()
            .map(|record| record.sequence.len())
            .collect::<Vec<_>>()
    };
    let named = lengths(InputFile::from(std::fs::File::open(&path).unwrap()));
    let standard_input = lengths(InputFile::as_standard_input(
        std::fs::File::open(&path).unwrap(),
    ));
    std::fs::remove_dir_all(&dir).unwrap();
    assert_eq!(named, vec![85]);
    assert_eq!(standard_input, vec![107]);
}

// `xSetSeqMol` + `AssignMolType` + `CSeqportUtil::Pack` through BLASTX's unified reader
// (`input_source` -> `fasta_reader::FastaInputSource::next_sequence_upper`, line RG-XT item X7-5):
// a protein subject with a nucleotide-named id (`gnl|TI|`) drops non-IUPACna letters with NCBI's
// nucleotide warning; a nucleotide query with a protein-named id (`sp|`) keeps protein letters
// without warning, maps U->T and non-IUPACna letters to N, and preserves lowercase masks.
#[test]
fn blastx_other_molecule_accession_follows_ncbis_reading() {
    let read_pd = |input: &[u8], protein: bool| {
        let mut warnings = Vec::new();
        let config = if protein {
            ReaderConfig::subject("BLASTX", true, false)
        } else {
            ReaderConfig::query("BLASTX", false, false)
        }
        .with_parse_deflines(true);
        let mut source = FastaInputSource::from_bytes(input, config);
        let mut records = Vec::new();
        while !source.end() {
            match source.next_sequence_upper(true, &mut |warning| {
                crate::algorithm::blastx::native::reader_warning(warning, &mut warnings)
            }) {
                Ok((record, masks)) => records.push(from_upper(record, masks)),
                Err(ReadError::Parse {
                    code: ParseErrorCode::Eof,
                    ..
                }) => break,
                Err(error) => panic!("unexpected read error: {error}"),
            }
        }
        (records, String::from_utf8(warnings).unwrap())
    };
    let (records, warnings) = read_pd(b">gnl|TI|12347 trace\nMKLVEFu\n", true);
    assert_eq!(records[0].sequence, b"MKVU");
    assert_eq!(records[0].lowercase_masks, vec![(3, 3)]);
    assert_eq!(
        warnings,
        "FASTA-Reader: Ignoring invalid nucleotide residues at position(s): On line 2: 3, 5-6\n"
    );

    // Nucleotide queries under protein-named accessions (including ranged ids and `>?` gaps).
    let nuc_cases: [(&[u8], &[u8], &[(i32, i32)]); 2] = [
        (
            b">sp|P99995.1|NAME_W\nACGUefILPQXZ*\n",
            b"ACGTNNNNNNNNN",
            &[(4, 5)],
        ),
        (b">sp|P99995.1:1-20\nae\n>?2\nfgA\n", b"ANNNNGA", &[(0, 5)]),
    ];
    for (input, seq, masks) in nuc_cases {
        let (records, warnings) = read_pd(input, false);
        assert_eq!(records[0].sequence, seq);
        assert_eq!(records[0].lowercase_masks, masks);
        assert!(warnings.is_empty(), "{warnings}");
    }

    // Lowercase mask across a `>?` gap line on a ranged nucleotide-named protein subject.
    let (records, warnings) = read_pd(b">gnl|TI|12347:1-20 trace\naeu\n>?2\nfgcA\n", true);
    assert_eq!(records[0].sequence, b"AUXXGCA");
    assert_eq!(records[0].lowercase_masks, vec![(0, 5)]);
    assert!(
        warnings.contains("On line 2: 2") && warnings.contains("On line 4: 1"),
        "{warnings}"
    );
}
