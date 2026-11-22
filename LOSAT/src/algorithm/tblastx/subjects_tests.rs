use super::*;
use crate::algorithm::tblastx::translation::generate_frames;
use crate::seqdb::volume::AmbiguityRun;

fn assert_same(got: &[QueryFrame], expected: &[QueryFrame], what: &str) {
    assert_eq!(got.len(), expected.len(), "{what}");
    for (g, e) in got.iter().zip(expected) {
        assert_eq!(
            (g.frame, &g.aa_seq, g.aa_len, g.orig_len),
            (e.frame, &e.aa_seq, e.aa_len, e.orig_len),
            "{what}"
        );
    }
}

/// The two methods against the frames of `with_search_bases` (the frames `process_subject`
/// made before its frames were made in windows), for the subjects of both arms.
fn check(subjects: Subjects<'_>, table: &GeneticCode) {
    for oid in 0..subjects.count() {
        let (letters, preliminary) = subjects.with_search_bases(oid, |letters, bases| {
            (
                generate_frames(letters, table),
                bases.map(|bases| generate_frames(bases, table)),
            )
        });
        let (_, frames, of_letters, kept) = subjects
            .x_preliminary_frames(oid, table, 0)
            .unwrap()
            .unwrap();
        assert_eq!(of_letters, preliminary.is_none(), "oid {oid}");
        assert_same(
            &frames,
            preliminary.as_deref().unwrap_or(&letters),
            "preliminary",
        );
        assert_same(
            &subjects.x_letter_frames(oid, table, kept).unwrap(),
            &letters,
            "letters",
        );
    }
}

fn letters(seed: u64, len: usize, alphabet: &[u8]) -> Vec<u8> {
    let mut state = seed | 1;
    (0..len)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            alphabet[(state % alphabet.len() as u64) as usize]
        })
        .collect()
}

#[test]
fn frames_of_records_are_those_of_their_letters_and_bases() {
    let mut records = Vec::new();
    for (i, len) in [0usize, 2, 3, 4, 50, 1000, 200_003].into_iter().enumerate() {
        records.push(FastaRecord::new(
            format!("a{i}"),
            b"",
            &letters(i as u64, len, b"ACGTacgt"),
        ));
        records.push(FastaRecord::new(
            format!("n{i}"),
            b"",
            &letters(i as u64 + 99, len, b"ACGTACGTacgtNNNRYKMSWBDHVnU-"),
        ));
    }
    // One N block across a window border.
    let mut block = letters(7, 400_000, b"ACGT");
    block[196_000..197_000].fill(b'N');
    records.push(FastaRecord::new("block", b"", &block));
    for id in [1, 4, 11] {
        check(Subjects::Fasta(&records), &GeneticCode::from_id(id));
    }
}

/// LOSAT_X_TBXRESOLVE: every cell of `ONE_BASE` against `ncbi4na_of` and `resolve_ncbi4na_base`:
/// a byte of the table gets the resolve's base without a draw, and every other byte has an
/// ncbi4na value of no base or of more than one.
#[test]
fn one_base_is_the_resolve_without_a_draw() {
    use crate::core::blast_encoding::{resolve_ncbi4na_base, NcbiRandom};
    for byte in 0..=u8::MAX {
        let code = ncbi4na_of(byte);
        let known = ONE_BASE[usize::from(byte)];
        if known == 0 {
            assert_ne!(code.count_ones(), 1, "byte {byte}");
            continue;
        }
        let (mut drawn, mut fresh) = (NcbiRandom::new(7), NcbiRandom::new(7));
        let base = b"ACGT"[usize::from(resolve_ncbi4na_base(code, &mut drawn))];
        assert_eq!(known, base, "byte {byte}");
        assert_eq!(drawn.get_rand(), fresh.get_rand(), "byte {byte}: a draw");
    }
}

/// LOSAT_X_TBXRESOLVE: the table walk against the letter-by-letter walk (`=0`), every frame,
/// for records of every IUPAC letter in both cases, `U`, gaps and bytes without a code, with
/// one ambiguous letter, and with an N block across a window border.
#[test]
fn the_table_walk_gives_the_letter_walks_frames() {
    let mut records = vec![
        b"N".to_vec(),
        b"acgtN".to_vec(),
        b"uUuUaN".to_vec(),
        letters(5, 3, b"ACGTN"),
    ];
    for (i, len) in [1usize, 2, 3, 4, 7, 50, 1000, 200_003]
        .into_iter()
        .enumerate()
    {
        records.push(letters(
            i as u64 + 11,
            len,
            b"ACGTacgtNnRrYyKkMmSsWwBbDdHhVvUu-*0 ",
        ));
        let mut one = letters(i as u64 + 31, len, b"ACGTacgt");
        one[len / 2] = b'R';
        records.push(one);
    }
    let mut block = letters(17, 400_000, b"ACGTacgu");
    block[196_000..197_000].fill(b'N');
    records.push(block);
    for id in [1, 4, 11] {
        let table = GeneticCode::from_id(id);
        for (i, record) in records.iter().enumerate() {
            assert_same(
                &x_resolved_record_frames(record, &table, true),
                &x_resolved_record_frames(record, &table, false),
                &format!("record {i} code {id}"),
            );
        }
    }
}

/// LOSAT_X_TBXSTORE: a record searched from its store (`PackedNucleotide::from_subject_letters`,
/// as `x_read_packed_subjects` packs it) has the record's length, the frames the record's
/// preliminary stage translates (its letters, or the bases its ambiguous letters resolve to
/// from `CRandom` of its length) and the frames of its letters, for every IUPAC letter in both
/// cases and gaps, with one ambiguous letter, and with an N block across a window border. A
/// record with `U` has no store (`ncbi4na_of` reads it as `T`, the store would not).
#[test]
fn stored_records_give_the_frames_of_their_records() {
    let mut records = Vec::new();
    for (i, len) in [0usize, 1, 2, 3, 4, 7, 50, 1000, 200_003]
        .into_iter()
        .enumerate()
    {
        records.push(FastaRecord::new(
            format!("a{i}"),
            b"",
            &letters(i as u64, len, b"ACGTacgt"),
        ));
        records.push(FastaRecord::new(
            format!("n{i}"),
            b"",
            &letters(i as u64 + 41, len, b"ACGTacgtNnRrYyKkMmSsWwBbDdHhVv-"),
        ));
        let mut one = letters(i as u64 + 61, len, b"ACGTacgt");
        if len > 0 {
            one[len / 2] = b'r';
        }
        records.push(FastaRecord::new(format!("o{i}"), b"", &one));
    }
    let mut block = letters(17, 400_000, b"ACGTacgt");
    block[196_000..197_000].fill(b'N');
    records.push(FastaRecord::new("block", b"", &block));
    let stores: Vec<PackedNucleotide> = records
        .iter()
        .map(|record| PackedNucleotide::from_subject_letters(record.seq()).expect("IUPAC"))
        .collect();
    let (fasta, stored) = (Subjects::Fasta(&records), Subjects::Stored(&stores));
    for id in [1, 4, 11] {
        let table = GeneticCode::from_id(id);
        for oid in 0..records.len() {
            let what = format!("record {oid} code {id}");
            let (length, frames, of_letters, _) =
                fasta.x_preliminary_frames(oid, &table, 0).unwrap().unwrap();
            let (got_length, got, got_of_letters, _) = stored
                .x_preliminary_frames(oid, &table, 0)
                .unwrap()
                .unwrap();
            assert_eq!(length, records[oid].seq().len(), "{what}");
            assert_eq!(got_length, length, "{what}");
            assert_eq!(got_of_letters, of_letters, "{what}");
            assert_same(&got, &frames, &what);
            assert_same(
                &stored.x_letter_frames(oid, &table, None).unwrap(),
                &fasta.x_letter_frames(oid, &table, None).unwrap(),
                &what,
            );
        }
    }
    assert!(PackedNucleotide::from_subject_letters(b"ACGU").is_none());
    assert!(PackedNucleotide::from_subject_letters(b"ACGu").is_none());
}

/// A `-db` stream gives the frames of the copy of the database
/// (`common::database::read_nucleotide_stores`) for every fixture database of nucleotides
/// (ambiguity runs, stored bases flipped under them, several volumes, parsed ids), in blocks
/// of several sequences and alone, with genetic codes 1, 2, 4 and 11; the letters' frames
/// come from the kept block and from a read of the sequence alike.
#[test]
fn streamed_sequences_give_the_frames_of_the_copy() {
    use crate::algorithm::common::x_db_stream::XDbStream;
    use crate::seqdb::{SeqDb, SeqType};
    let tests = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let mut names: Vec<String> = ["seqdb/db/nucl_v4", "seqdb/db/nucl_v5", "seqdb/db/nucl_mv"]
        .map(String::from)
        .to_vec();
    for entry in std::fs::read_dir(tests.join("db_search")).unwrap() {
        let name = entry.unwrap().file_name().into_string().unwrap();
        if let Some(stem) = name.strip_suffix(".nsq") {
            names.push(format!("db_search/{stem}"));
        }
    }
    let mut ambiguous = 0;
    for name in &names {
        let path = tests.join(name);
        let db = SeqDb::open(path.to_str().unwrap(), SeqType::Nucleotide).unwrap();
        let copy =
            crate::algorithm::common::database::read_nucleotide_stores(&db, "TBLASTX").unwrap();
        ambiguous += copy
            .iter()
            .filter(|store| !store.ambiguity_runs().is_empty())
            .count();
        for max_bytes in [0, 64, u64::MAX] {
            let stream = XDbStream::new(&db, 1, "TBLASTX")
                .unwrap()
                .x_regrid(max_bytes);
            let (streamed, stored) = (Subjects::Stream(&stream), Subjects::Stored(&copy));
            assert_eq!(streamed.count(), stored.count(), "{name}");
            assert_eq!(streamed.longest_and_total(), stored.longest_and_total());
            for id in [1, 2, 4, 11] {
                let table = GeneticCode::from_id(id);
                check(streamed, &table);
                for oid in 0..copy.len() {
                    let (length, frames, of_letters, kept) = streamed
                        .x_preliminary_frames(oid, &table, 0)
                        .unwrap()
                        .unwrap();
                    let (want_length, want, want_of_letters, _) = stored
                        .x_preliminary_frames(oid, &table, 0)
                        .unwrap()
                        .unwrap();
                    assert_eq!(length, copy[oid].len(), "{name} {oid}");
                    assert_eq!(want_length, length, "{name} {oid}");
                    // A subject shorter than `shortest` gets no frame, with either side.
                    for subjects in [streamed, stored] {
                        let at = |shortest| {
                            subjects
                                .x_preliminary_frames(oid, &table, shortest)
                                .unwrap()
                                .is_some()
                        };
                        assert!(at(length) && !at(length + 1), "{name} {oid}");
                    }
                    assert_eq!(of_letters, want_of_letters, "{name} {oid}");
                    assert_same(&frames, &want, "preliminary");
                    let want = stored.x_letter_frames(oid, &table, None).unwrap();
                    assert_same(
                        &streamed.x_letter_frames(oid, &table, None).unwrap(),
                        &want,
                        "letters read again",
                    );
                    assert_same(
                        &streamed.x_letter_frames(oid, &table, kept).unwrap(),
                        &want,
                        "letters of the kept block",
                    );
                }
            }
        }
    }
    assert!(ambiguous > 0, "a fixture with ambiguity runs");
}

/// A read of a `-db` stream that fails during the search (`XDbStream::x_failing`) is the
/// `XReadError` of its subject's frames, preliminary and of the letters read again. A query
/// batch's failure (`XReadFailure`) keeps the error of the smallest ordinal id whose read
/// fails, whichever fails first, and skips the subjects after it only, so every subject
/// before it is searched at any number of threads.
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn a_failed_read_of_the_stream_fails_the_batch_at_its_smallest_subject() {
    use crate::algorithm::common::x_db_stream::{XDbStream, XReadError};
    use crate::seqdb::{SeqDb, SeqType};
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/db/nucl_v4");
    let db = SeqDb::open(path.to_str().unwrap(), SeqType::Nucleotide).unwrap();
    let table = GeneticCode::from_id(1);
    let count = db.num_seqs() as usize;
    let (smallest, last) = (count / 2, count - 1);
    let stream = XDbStream::new(&db, 4, "TBLASTX")
        .unwrap()
        .x_regrid(0)
        .x_failing(&[last, smallest], true);
    let subjects = Subjects::Stream(&stream);
    for oid in 0..count {
        let failing = oid == smallest || oid == last;
        let preliminary = subjects.x_preliminary_frames(oid, &table, 0);
        let letters = subjects.x_letter_frames(oid, &table, None);
        for error in [preliminary.err(), letters.err()] {
            assert_eq!(error.is_some(), failing, "{oid}");
            assert!(error.is_none_or(|error| error.is::<XReadError>()), "{oid}");
        }
    }
    let failure = XReadFailure::new();
    failure.fail(last, anyhow::anyhow!("last"));
    assert!(failure.skips(last + 1) && !failure.skips(last) && !failure.skips(smallest + 1));
    failure.fail(smallest, anyhow::anyhow!("smallest"));
    failure.fail(last, anyhow::anyhow!("last again"));
    assert!(failure.skips(smallest + 1) && !failure.skips(smallest));
    assert_eq!(failure.take().unwrap().to_string(), "smallest");
    assert!(failure.take().is_none());
    let message = format!("the read of database sequence {smallest} failed");
    for _ in 0..20 {
        let failure = XReadFailure::new();
        let next = AtomicUsize::new(0);
        let searched = Mutex::new(Vec::new());
        std::thread::scope(|scope| {
            for _ in 0..4 {
                scope.spawn(|| loop {
                    let oid = next.fetch_add(1, Ordering::Relaxed);
                    if oid >= count {
                        break;
                    }
                    if failure.skips(oid) {
                        continue;
                    }
                    match subjects.x_preliminary_frames(oid, &table, 0) {
                        Ok(_) => searched.lock().unwrap().push(oid),
                        Err(error) => failure.fail(oid, error),
                    }
                });
            }
        });
        let error = failure.take().expect("a read failed");
        assert!(error.to_string().contains(&message), "{error}");
        let searched = searched.into_inner().unwrap();
        assert!((0..smallest).all(|oid| searched.contains(&oid)));
    }
}

#[test]
fn frames_of_stored_sequences_are_those_of_their_codes_and_stored_bases() {
    let mut stores = Vec::new();
    for (i, len) in [3usize, 7, 64, 1001, 300_001].into_iter().enumerate() {
        let packed: Vec<u8> = letters(i as u64 + 5, len / 4 + 1, &(0..=255).collect::<Vec<u8>>());
        stores.push(PackedNucleotide::from_stored(&packed, len, &[]));
        // Runs of ambiguity codes, also of A/C/G/T codes other than the stored base, and
        // one across a window border; the stored bases under them are arbitrary (D2).
        let mut runs = vec![
            AmbiguityRun {
                offset: 0,
                length: 1,
                ncbi4na: 15,
            },
            AmbiguityRun {
                offset: 2,
                length: 1,
                ncbi4na: 5,
            },
        ];
        if len > 1000 {
            runs.push(AmbiguityRun {
                offset: 500,
                length: 3,
                ncbi4na: 1,
            });
            runs.push(AmbiguityRun {
                offset: 900,
                length: 40,
                ncbi4na: 14,
            });
        }
        if len > 200_000 {
            runs.push(AmbiguityRun {
                offset: 196_600,
                length: 50,
                ncbi4na: 15,
            });
        }
        stores.push(PackedNucleotide::from_stored(&packed, len, &runs));
    }
    for id in [1, 2, 11] {
        check(Subjects::Stored(&stores), &GeneticCode::from_id(id));
    }
}
