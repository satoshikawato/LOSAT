use super::*;
use crate::algorithm::blastn::hsp::{BlastnHsp, BlastnHspList};
#[cfg(unix)]
use crate::cli::NativeError;

fn fixture(name: &str) -> String {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/seqdb/db")
        .join(name)
        .to_str()
        .expect("fixture path is UTF-8")
        .to_string()
}

fn open(name: &str) -> BlastnDatabase {
    BlastnDatabase::open(&fixture(name), &ApplicationSettings::default())
        .unwrap_or_else(|error| panic!("{name}: {error:#}"))
}

fn open_dbmask(name: &str) -> BlastnDatabase {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/db_search/dbmask")
        .join(name);
    BlastnDatabase::open(path.to_str().unwrap(), &ApplicationSettings::default())
        .unwrap_or_else(|error| panic!("{name}: {error:#}"))
}

// blast_util.c:197-218 and seqdb.hpp:286-306 (cited at `seq_ranges_into`): the soft ranges
// of a masked database sequence are the complement of its masks in stored order, inverted
// members included; a sequence without masks has the one range of its letters.
#[test]
fn soft_ranges_are_the_complement_of_the_masks_in_stored_order() {
    let mut ranges = vec![(5, 6)];
    for (masks, expected) in [
        (&[][..], &[(0, 1200)][..]),
        (
            &[(600, 700), (100, 200)],
            &[(0, 600), (700, 100), (200, 1200)],
        ),
        (
            &[(100, 300), (200, 400)],
            &[(0, 100), (300, 200), (400, 1200)],
        ),
        (
            &[(0, 100), (1100, 1200)],
            &[(0, 0), (100, 1100), (1200, 1200)],
        ),
        (&[(0, 1200)], &[(0, 0), (1200, 1200)]),
    ] {
        crate::core::subject_split::complement_masks(masks, 1200, &mut ranges);
        assert_eq!(ranges, expected, "{masks:?}");
    }
}

// The masks of the M4a-2 edge database `edge1` (NCBI's makeblastdb of hand-made maskinfo; NCBI's
// blastdbcmd `-outfmt "%o %m11"` reads `600-700;100-200;` for sequence 4): read with the masks
// of the block that holds the sequence (one read for all its sequences, review 1 F3 of
// M4a-2a), the soft ranges in stored order, the display masks as stored.
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1518-1538,
// 1665-1673, 1694-1700 (cited at `mask_hard_letters`)
// `-db_hard_mask`'s traceback letters (series NP, M4a-2b decision H2): after the decode, every
// masked position of a decoded part is BLASTNA 'B' (10), as `s_SeqDBMaskSequence` writes
// NCBI4na 14 before `SeqDB_ncbina8_to_blastna8`, on the whole sequence and on each slice of
// the partial decode; the other letters (and the fences, the sentinel bytes) stay. Against a
// model of the C loop per slice, for random sorted masks and slices, then on a5's s1 (masks
// 1000-1100 and 2000-2024) through `traceback_listed` (whole) and the ranged decode.
#[test]
fn hard_masked_letters_are_b_in_each_decoded_part() {
    let model = |letters: &[u8], masks: &[(u32, u32)], slices: &[(usize, usize)]| {
        let mut seq = letters.to_vec();
        for &(begin, end) in slices {
            let mut i = 0;
            while i < masks.len() && masks[i].1 as usize <= begin {
                i += 1;
            }
            while i < masks.len() && (masks[i].0 as usize) < end {
                for j in (masks[i].0 as usize).max(begin)..(masks[i].1 as usize).min(end) {
                    seq[1 + j] = 10;
                }
                i += 1;
            }
        }
        seq
    };
    let mut seed = 0x2545_F491_4F6C_DD1Du64;
    let mut next = |bound: u64| {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed % bound) as usize
    };
    for _ in 0..500 {
        let length = 1 + next(400);
        let letters: Vec<u8> = (0..length + 2).map(|_| next(16) as u8).collect();
        let mut masks = Vec::new();
        let mut at = 0usize;
        while masks.len() < 5 {
            let start = at + next(80);
            let end = start + next(30);
            if end > length {
                break;
            }
            masks.push((start as u32, end as u32));
            at = end;
        }
        let mut slices = Vec::new();
        let mut at = 0usize;
        while at < length {
            let begin = at + next(60);
            let end = (begin + 1 + next(120)).min(length);
            if begin >= end {
                break;
            }
            slices.push((begin, end));
            at = end + 1;
        }
        for parts in [&[(0, length)][..], &slices] {
            let mut found = letters.clone();
            mask_hard_letters(&mut found, &masks, parts);
            assert_eq!(found, model(&letters, &masks, parts), "{masks:?} {parts:?}");
        }
    }

    let database = open_dbmask("a5");
    let stored = database.db.read_sequence(0).unwrap();
    let length = stored.nucleotide_len(0).unwrap();
    let whole = stored.nucleotide_blastna(0).unwrap();
    let masked =
        |position: usize| (1000..1100).contains(&position) || (2000..2024).contains(&position);
    let found = crate::utils::threading::with_search_pool(1, "blastn", |pool| {
        traceback_listed(
            &database,
            Some(11),
            pool,
            false,
            &[(0u32, ())],
            &|| (),
            &|_: &mut (), _: u32, subject: BlastnSubjectRef<'_>, _: &()| {
                Some(subject.blastna.expect("the traceback's letters").to_vec())
            },
            &|_: &()| Vec::new(),
        )
    })
    .unwrap();
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].len(), length);
    for (position, &code) in found[0].iter().enumerate() {
        let want = if masked(position) {
            10
        } else {
            whole[1 + position]
        };
        assert_eq!(code, want, "position {position}");
    }
    // The ranged decode: slices across both masks and between them, fences kept.
    let ranges = [(900usize, 1050usize), (1080, 1500), (1990, 2010)];
    let mut ranged = stored.x_nucleotide_blastna_ranges(0, &ranges).unwrap();
    let before = ranged.clone();
    mask_hard_letters(&mut ranged, &[(1000, 1100), (2000, 2024)], &ranges);
    for position in 0..length {
        let decoded = ranges
            .iter()
            .any(|&(from, to)| (from..to).contains(&position));
        let want = if decoded && masked(position) {
            10
        } else {
            before[1 + position]
        };
        assert_eq!(ranged[1 + position], want, "position {position}");
    }
}

#[test]
fn a_masked_search_reads_each_sequence_s_masks() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/db_search/dbmask/edge1")
        .to_str()
        .unwrap()
        .to_string();
    let database = BlastnDatabase::open(&path, &ApplicationSettings::default()).unwrap();
    let block = database.db.x_block_range(0, 5).unwrap();
    assert_eq!((block.first, block.end), (0, 5));
    let mut masks = MaskBlock::default();
    let mut ranges = Vec::new();
    database
        .seq_ranges_into(block, 4, 11, 1200, &mut masks, &mut ranges)
        .unwrap();
    assert_eq!(ranges, [(0, 600), (700, 100), (200, 1200)]);
    // Sequence 0 has no masks of algorithm 11.
    database
        .seq_ranges_into(block, 0, 11, 3000, &mut masks, &mut ranges)
        .unwrap();
    assert_eq!(ranges, [(0, 3000)]);
    // Another algorithm reads the block again.
    database
        .seq_ranges_into(block, 4, 15, 1200, &mut masks, &mut ranges)
        .unwrap();
    assert_eq!(ranges, [(0, 301), (302, 500), (503, 703), (706, 1200)]);
    // A block of one sequence gives its masks alike.
    let one = database.db.x_block_range(4, 5).unwrap();
    database
        .seq_ranges_into(one, 4, 11, 1200, &mut masks, &mut ranges)
        .unwrap();
    assert_eq!(ranges, [(0, 600), (700, 100), (200, 1200)]);
    let shown = database.display_masks(4, 12, &mut masks).unwrap();
    assert_eq!(
        shown,
        [
            crate::utils::dust::MaskedInterval::new(100, 300),
            crate::utils::dust::MaskedInterval::new(200, 400)
        ]
    );
}

// DECISIONS 3 Q8 of line DBID: NCBI's choice of the Seq-ids of a database with parsed
// Seq-ids depends on SEQ_ID_PREFER_ACCESSION_OVER_GI and SEQ_ID_AVOID_GI.
#[test]
fn a_seq_id_variable_refuses_a_database_with_parsed_seqids_only() {
    let path = |name: &str| {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join(name)
            .to_str()
            .expect("fixture path is UTF-8")
            .to_string()
    };
    let settings = ApplicationSettings {
        seq_id_rank_variable: Some("SEQ_ID_PREFER_ACCESSION_OVER_GI"),
        ..ApplicationSettings::default()
    };
    let parsed = path("tests/db_search/ids/nucl_v5");
    BlastnDatabase::open(&parsed, &ApplicationSettings::default())
        .expect("a database with parsed Seq-ids opens");
    let Err(error) = BlastnDatabase::open(&parsed, &settings) else {
        panic!("a database with parsed Seq-ids opened with SEQ_ID_PREFER_ACCESSION_OVER_GI")
    };
    assert_eq!(
        format!("{error:#}"),
        "the environment variable SEQ_ID_PREFER_ACCESSION_OVER_GI, which changes the Seq-ids that NCBI BLAST+ chooses for a database with parsed Seq-ids, is not supported by LOSAT's BLASTN"
    );
    BlastnDatabase::open(&path("tests/seqdb/db/nucl_v4"), &settings)
        .expect("a database without parsed Seq-ids opens with SEQ_ID_PREFER_ACCESSION_OVER_GI");
}

/// LOSAT/tests/seqdb/expected/nucl_mv.entries.tsv (blastdbcmd): the lengths, letters and
/// titles in ordinal order.
fn entries() -> Vec<(usize, String, String)> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/seqdb/expected/nucl_mv.entries.tsv");
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            (
                fields[1].parse().unwrap(),
                fields[3].to_string(),
                fields[4].to_string(),
            )
        })
        .collect()
}

/// A temporary directory of this test (removed by the test).
fn temp_dir(tag: &str) -> std::path::PathBuf {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-blastn-db-{tag}-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn test_hsp(s_idx: u32, raw_score: i32) -> BlastnHsp {
    BlastnHsp {
        identity: 100.0,
        length: 10,
        mismatch: 0,
        gapopen: 0,
        q_start: 1,
        q_end: 10,
        s_start: 1,
        s_end: 10,
        e_value: 1e-5,
        bit_score: 30.0,
        num_ident: 10,
        query_frame: 1,
        query_length: 20,
        q_idx: 0,
        s_idx,
        raw_score,
        internal_q_offset_0: 0,
        internal_q_end_0: 10,
        internal_s_offset_0: 0,
        internal_s_end_0: 10,
        internal_query_context_offset: 0,
        gap_info: None,
        num_positives: 10,
    }
}

fn hit_list(lists: Vec<(u32, f64)>, count: usize) -> Option<BlastnHitList> {
    let mut hit_list = BlastnHitList::new(10);
    hit_list.hsplist_array = lists
        .into_iter()
        .map(|(oid, best_evalue)| BlastnHspList {
            oid,
            query_index: 0,
            hsps: vec![test_hsp(oid, 60)],
            best_evalue,
        })
        .collect();
    hit_list.hsplist_count = count;
    Some(hit_list)
}

// Review 2 of line DBID pfoff-3, L5: BLASTN's model of what one query's pairwise report reads
// (`pairwise_read_order` of its hit list) is the writer's group of the rows that
// `pairwise::pairwise_hits` makes of that list (`x_pairwise_par::writer_reads`), for every
// sort and count; the rows name each subject by its rank, here its ordinal id.
#[test]
fn the_pairwise_reads_are_those_of_the_writers_rows() {
    use crate::algorithm::blastn::pairwise::{pairwise_hits, DisplayMasks, Rows};
    use crate::blastinput::fasta_reader::FastaRecord;
    let hsp = |s_idx: u32, bits: f64, e_value: f64, num_ident: usize, from: usize| {
        let mut hsp = test_hsp(s_idx, (bits * 2.0) as i32);
        (hsp.bit_score, hsp.e_value, hsp.num_ident) = (bits, e_value, num_ident);
        (hsp.q_start, hsp.q_end) = (from + 1, from + 10);
        (hsp.internal_q_offset_0, hsp.internal_q_end_0) = (from as i32, from as i32 + 10);
        hsp
    };
    let list = |oid: u32, hsps: Vec<BlastnHsp>| BlastnHspList {
        oid,
        query_index: 0,
        best_evalue: hsps[0].e_value,
        hsps,
    };
    let mut hit_list = BlastnHitList::new(10);
    hit_list.hsplist_array = vec![
        list(0, vec![hsp(0, 40.0, 1e-9, 4, 6), hsp(0, 30.0, 1e-7, 9, 0)]),
        list(1, vec![hsp(1, 60.0, 1e-12, 2, 3)]),
        list(2, vec![hsp(2, 35.0, 1e-8, 10, 1), hsp(2, 20.0, 1e-3, 8, 9)]),
    ];
    hit_list.hsplist_count = 3;
    let hit_lists = [Some(hit_list)];
    let letters = b"ACGTACGTTGCAACGTTGCAAGGTCC";
    let queries = [FastaRecord::new("q", b"q", letters)];
    let subjects: Vec<FastaRecord> = (0..3)
        .map(|_| FastaRecord::new("s", b"s", letters))
        .collect();
    let placements = crate::blastinput::seq_range::Placements::default();
    let masks = DisplayMasks {
        query: &[Vec::new()],
        query_offsets: &placements,
        subject: None,
        database: None,
    };
    let rows = Rows::new(
        &queries,
        crate::algorithm::blastn::x_subject_store::SubjectLetters::new(&subjects, None),
        &masks,
    );
    let rows = pairwise_hits(&hit_lists, &rows, None, None).unwrap();
    for hits in [None, Some(0), Some(1), Some(2), Some(3), Some(4)] {
        for hsps in [None, Some(0), Some(1), Some(2), Some(3), Some(4)] {
            for (descriptions, alignments) in [(0, 0), (1, 2), (2, 1), (500, 500)] {
                let reading = crate::report::pairwise::PairwiseReading {
                    sort: crate::report::hit_sort::ReportSort { hits, hsps },
                    shown: crate::report::db_subject::Shown {
                        descriptions,
                        alignments,
                    },
                };
                assert_eq!(
                    pairwise_read_order(&hit_lists, reading, |_| None),
                    crate::report::x_pairwise_par::writer_reads(
                        &rows,
                        |_, _| false,
                        reading,
                        false,
                        None
                    ),
                    "{reading:?}"
                );
            }
        }
    }
}

// NCBI BLAST+ 2.17.0 `blastdbcmd -info` of the fixture
// (LOSAT/tests/seqdb/expected/nucl_mv.info.tsv), and align_format_util.cpp:385-391.
#[test]
fn database_numbers_come_from_the_index_files() {
    let database = open("nucl_mv");
    assert_eq!(database.metadata().db_num_seqs, 12);
    assert_eq!(database.metadata().db_len_total, 850);
    let shown = database.pairwise_database();
    assert_eq!(shown.definition, b"nucleotide in volumes");
    assert_eq!(shown.databases[0].date, "Oct 12, 2026  4:12 PM");
    assert_eq!((shown.num_sequences, shown.total_letters), (12, 850));
    assert_eq!(database.name(), fixture("nucl_mv"));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:904-911
// An alias file's `LENGTH` of 2^32 letters or more over the 16 sequences of
// `tests/db_search/taxfilter/tfn` (`tests/db_search/dbcounts`, rows db.cnt.l4295.bn and
// db.cnt.avg2p32.bn) is NCBI's `Int8` on every target, wasm32's 32-bit `usize` included: the
// database opens and the statistics and the reports read the whole total.
#[test]
fn totals_of_2_to_the_32_letters_or_more_are_held_whole() {
    for (name, total) in [
        ("l4295", 4_294_967_296_000_000_000),
        ("l4294", 4_294_967_295_999_999_999),
        ("avg2p32", 68_719_476_736),
    ] {
        assert!(total > i64::from(u32::MAX), "{name}");
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/db_search/dbcounts")
            .join(name);
        let database =
            BlastnDatabase::open(path.to_str().unwrap(), &ApplicationSettings::default())
                .unwrap_or_else(|error| panic!("{name}: {error:#}"));
        let metadata = database.metadata();
        assert_eq!(
            (metadata.db_len_total, metadata.db_num_seqs),
            (total, 16),
            "{name}"
        );
        let shown = database.pairwise_database();
        assert_eq!(
            (shown.num_sequences, shown.total_letters),
            (16, total),
            "{name}"
        );
    }
}

// blast_app_util.hpp:185-203 and seqdb.cpp:155-165: NCBI's errors with exit status 2,
// LOSAT's refusals in LOSAT's wording, naming BLASTN.
// (A path in an error text: Unix spelling.)
#[cfg(unix)]
#[test]
fn open_errors_are_ncbis_database_errors_or_losats_refusals() {
    let native = |error: anyhow::Error| {
        let native = error
            .downcast_ref::<NativeError>()
            .unwrap_or_else(|| panic!("not a native error: {error:#}"));
        (
            native.exit,
            String::from_utf8(native.message.clone()).unwrap(),
        )
    };
    let settings = ApplicationSettings::default();
    let Err(error) = BlastnDatabase::open("", &settings) else {
        panic!("an empty name opened")
    };
    assert_eq!(
        native(error),
        (
            2,
            "LOSAT database error: Database name is required.\n".to_string()
        )
    );
    // A protein database is not found as a nucleotide one.
    for name in ["no_such_database", "prot_v4"] {
        let Err(error) = BlastnDatabase::open(&fixture(name), &settings) else {
            panic!("{name} opened")
        };
        let (exit, message) = native(error);
        assert_eq!(exit, 2);
        assert!(
            message.starts_with(
                "LOSAT database error: No alias or index file found for nucleotide database ["
            ),
            "{message}"
        );
    }
    // The reader's refusals name BLASTN.
    let error = database_error(SeqDbError::Unsupported {
        path: "/data/db.nal".into(),
        message: "the alias file filters sequences".to_string(),
    });
    assert!(error.downcast_ref::<NativeError>().is_none());
    assert!(
        format!("{error:#}").ends_with("; this database is not supported by LOSAT's BLASTN"),
        "{error:#}"
    );
    let read = search_read_error(SeqDbError::Open {
        path: "/data/db.nsq".into(),
        source: std::io::Error::from_raw_os_error(13),
    })
    .downcast::<PreliminaryReadError>()
    .expect("NCBI's error of the search stage");
    assert_eq!(
        read.to_string(),
        "LOSAT engine error: Database memory map file error"
    );
}

/// A kernel that records what it is given and keeps every third sequence.
fn recording_kernel(
    seen: &std::sync::Mutex<Vec<(u32, usize, Vec<u8>)>>,
) -> impl Fn(&mut u32, u32, BlastnSubjectRef<'_>) -> Option<u32> + Sync + Send + '_ {
    move |calls, ordinal, subject| {
        *calls += 1;
        assert!(subject.letters.is_none() && subject.blastna.is_none());
        seen.lock()
            .unwrap()
            .push((ordinal, subject.length, subject.packed.to_vec()));
        (ordinal % 3 == 0).then_some(ordinal)
    }
}

fn run_blocks(
    database: &BlastnDatabase,
    threads: usize,
    block_bytes: u64,
) -> (Vec<(u32, usize, Vec<u8>)>, Vec<u32>) {
    let seen = std::sync::Mutex::new(Vec::new());
    let mut merged = Vec::new();
    crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
        search_blocks_of(
            database,
            pool,
            threads > 1,
            block_bytes,
            2,
            &|| 0u32,
            &recording_kernel(&seen),
            &mut |ordinal, kept| {
                assert_eq!(ordinal, kept);
                merged.push(ordinal);
            },
        )
    })
    .unwrap();
    let mut seen = seen.into_inner().unwrap();
    seen.sort_by_key(|entry| entry.0);
    (seen, merged)
}

// blast_engine.c:1407-1426 (the sequences in the order of their ordinal ids) and
// blast_engine.c:1401-1405 (the stored 2-bit data): the block driver gives the kernel
// every sequence once with its stored data, and merges what it keeps in ordinal order,
// whatever the block size and the thread count, across the six volumes of the fixture.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn blocks_give_every_sequence_once_and_merge_in_ordinal_order() {
    let database = open("nucl_mv");
    let single = open("nucl_v4");
    let entries = entries();
    let (reference, merged) = run_blocks(&database, 1, u64::MAX);
    assert_eq!(reference.len(), 12);
    for (ordinal, (seen_ordinal, length, packed)) in reference.iter().enumerate() {
        assert_eq!(*seen_ordinal as usize, ordinal);
        assert_eq!(*length, entries[ordinal].0, "sequence {ordinal}");
        // The stored layout of the subject set's packed letters: length / 4 + 1 bytes,
        // the count of the last byte's bases in its two low bits.
        assert_eq!(packed.len(), length / 4 + 1, "sequence {ordinal}");
        assert_eq!(usize::from(packed[length / 4] & 3), length % 4);
    }
    assert_eq!(merged, [0, 3, 6, 9]);
    assert_eq!(run_blocks(&single, 1, u64::MAX).0, reference);
    for threads in [1usize, 2, 4] {
        for block_bytes in [1u64, 16, 64, 1 << 20] {
            for db in [&database, &single] {
                let (seen, order) = run_blocks(db, threads, block_bytes);
                assert_eq!(
                    seen, reference,
                    "{threads} threads, blocks of {block_bytes}"
                );
                assert_eq!(order, merged, "{threads} threads, blocks of {block_bytes}");
            }
        }
    }
}

// blast_traceback.c:819-828 and 1567-1573: the traceback reads each listed sequence
// again by ordinal id, in BLASTNA with its ambiguity codes; the results come in the
// order of the listed ordinal ids at every thread count.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn listed_sequences_are_read_again_with_their_ambiguity_codes() {
    let database = open("nucl_mv");
    let entries = entries();
    // n02, n06 (every ambiguity code), n08 (all N), n12 (ambiguity at both ends).
    let listed: Vec<(u32, u32)> = vec![(1, 10), (5, 50), (7, 70), (11, 110)];
    for threads in [1usize, 2, 4] {
        let found = crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
            traceback_listed(
                &database,
                None,
                pool,
                threads > 1,
                &listed,
                &|| (),
                &|_: &mut (), ordinal: u32, subject: BlastnSubjectRef<'_>, kept: &u32| {
                    assert_eq!(*kept, ordinal * 10);
                    let blastna = subject.blastna.expect("the traceback's letters");
                    assert_eq!(blastna.len(), subject.length);
                    assert_eq!(subject.packed.len(), subject.length / 4 + 1);
                    let letters: String = blastna
                        .iter()
                        .map(|&code| BLASTNA_TO_IUPACNA[usize::from(code)] as char)
                        .collect();
                    (ordinal != 7).then_some((ordinal, letters))
                },
                &|_: &u32| Vec::new(),
            )
        })
        .unwrap();
        let expected: Vec<(u32, String)> = [1u32, 5, 11]
            .iter()
            .map(|&ordinal| (ordinal, entries[ordinal as usize].1.clone()))
            .collect();
        assert_eq!(found, expected, "{threads} threads");
    }
}

// A traceback without nested parallel work makes a scratch per worker at most (a sequence
// takes an idle one; `map_init` made one per rayon task, at least twice the threads), and the
// results come in the order of `listed`: at 1 to 8 threads, uneven work per sequence so that
// the threads steal from each other.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn the_parallel_traceback_reuses_the_idle_scratches() {
    use std::sync::atomic::{AtomicU64, Ordering};
    let database = open("nucl_mv");
    let work = |entry: u32| -> u64 {
        let rounds = if entry % 37 == 5 {
            200_000
        } else {
            50 + (entry * 7919) % 3_000
        };
        (0..u64::from(rounds)).fold(u64::from(entry), |acc, x| {
            acc.wrapping_mul(31).wrapping_add(x)
        })
    };
    // Each sequence 50 times, in ascending ordinal order, numbered in `kept`.
    let listed: Vec<(u32, u32)> = (0..12u32 * 50).map(|i| (i / 50, i)).collect();
    let want: Vec<(u32, u64)> = listed
        .iter()
        .filter(|&&(ordinal, _)| ordinal % 3 != 2)
        .map(|&(_, entry)| (entry, work(entry)))
        .collect();
    for threads in 1..=8usize {
        let scratches = AtomicU64::new(0);
        let found = crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
            traceback_listed(
                &database,
                None,
                pool,
                threads > 1,
                &listed,
                &|| {
                    scratches.fetch_add(1, Ordering::Relaxed);
                },
                &|_: &mut (), ordinal: u32, _: BlastnSubjectRef<'_>, kept: &u32| {
                    (ordinal % 3 != 2).then(|| (*kept, work(*kept)))
                },
                &|_: &u32| Vec::new(),
            )
        })
        .unwrap();
        assert_eq!(found, want, "{threads} threads");
        let made = scratches.load(Ordering::Relaxed);
        assert!(
            (1..=threads as u64).contains(&made),
            "{threads} threads: {made} scratches"
        );
    }
}

// blast_seqsrc.c:467-519: each HSP widened by 1024 bases within the sequence, sorted by
// start, merged unless it starts more than 1024 bases after the range before.
#[test]
fn partial_ranges_are_ncbis() {
    assert_eq!(x_partial_ranges_padded([], 100, 1024), vec![]);
    assert_eq!(
        x_partial_ranges_padded([(10, 20)], 100, 1024),
        vec![(0, 100)]
    );
    // 3976..7024, 5476..8024, 8976..11124: each within 1024 of the one before.
    assert_eq!(
        x_partial_ranges_padded([(6500, 7000), (5000, 6000), (10000, 10100)], 20000, 1024),
        vec![(3976, 11124)]
    );
    // 0..1124 and 3172..6124: 3172 > 1124 + 1024.
    assert_eq!(
        x_partial_ranges_padded([(4196, 5100), (0, 100)], 5500, 1024),
        vec![(0, 1124), (3172, 5500)]
    );
    // 0..1124 and 2148..: exactly 1024 apart merges.
    assert_eq!(
        x_partial_ranges_padded([(0, 100), (3172, 3200)], 10000, 1024),
        vec![(0, 4224)]
    );
    // A range inside the one before.
    assert_eq!(
        x_partial_ranges_padded([(5000, 9000), (6000, 6100)], 20000, 1024),
        vec![(3976, 10024)]
    );
}

// seqdbvol.cpp:1565-1595 (`GetAmbigPartialSeq`): the bases of the ranges as the whole
// decode gives them, a fence before and after each range inside the sequence, the two
// sentinels, and nothing written elsewhere.
#[test]
fn blastna_in_ranges_is_the_whole_decode_inside_them() {
    let database = open("nucl_mv");
    let mut checked = 0;
    for ordinal in 0..12u32 {
        let block = database.db.read_sequence(ordinal).unwrap();
        let length = block.nucleotide_len(0).unwrap();
        let whole = block.nucleotide_blastna(0).unwrap();
        let mut lists: Vec<Vec<(usize, usize)>> = vec![
            vec![(0, length)],
            vec![(1, length.saturating_sub(1).max(1))],
        ];
        if length >= 3 {
            lists.push(vec![(0, 1), (length - 1, length)]);
        }
        if length >= 12 {
            lists.push(vec![(2, 5), (7, 9), (11, length)]);
            lists.push(vec![(0, 3), (length / 2, length / 2 + 2)]);
        }
        for ranges in lists {
            let ranges: Vec<(usize, usize)> = ranges.into_iter().filter(|r| r.0 < r.1).collect();
            let ranged = block.x_nucleotide_blastna_ranges(0, &ranges).unwrap();
            assert_eq!(ranged.len(), length + 2);
            let mut expected = vec![0u8; length + 2];
            for &(from, to) in &ranges {
                if from > 0 {
                    expected[from] = crate::seqdb::FENCE_SENTRY;
                }
                if to < length {
                    expected[to + 1] = crate::seqdb::FENCE_SENTRY;
                }
            }
            for &(from, to) in &ranges {
                expected[from + 1..=to].copy_from_slice(&whole[from + 1..=to]);
            }
            expected[0] = whole[0];
            expected[length + 1] = whole[length + 1];
            assert_eq!(ranged, expected, "sequence {ordinal}, ranges {ranges:?}");
            checked += 1;
        }
    }
    assert!(checked >= 36);
}

// blast_seqalign.cpp:1569-1582 and seqdbvol.cpp:2302-2310: NCBI opens a volume's
// header file when the Seq-aligns of a batch's final hit lists are made, query by query
// and in the order of each hit list; a list past the hit list's count is not read, and
// a title read before is not read again.
// (A path in the system's temporary directory: Unix spelling.)
#[cfg(unix)]
#[test]
fn titles_are_read_for_the_final_hit_lists_in_ncbis_order() {
    let source = open("nucl_mv");
    // The volume of each ordinal id.
    let mut volume_of = Vec::new();
    for (volume, vol) in source.db.volumes().iter().enumerate() {
        volume_of.extend(std::iter::repeat_n(volume, vol.num_oids() as usize));
    }
    assert_eq!(volume_of.len(), 12);
    let last = 11u32;
    let first_volume_ordinal = 0u32;
    assert_ne!(volume_of[last as usize], volume_of[0]);
    // A sequence of a third volume, read in the first batch.
    let middle = (0..12u32)
        .find(|&ordinal| {
            volume_of[ordinal as usize] != volume_of[0]
                && volume_of[ordinal as usize] != volume_of[last as usize]
        })
        .unwrap();
    let dir = temp_dir("nohdr");
    let missing = [
        format!("nucl_mv.{:02}.nhr", volume_of[0]),
        format!("nucl_mv.{:02}.nhr", volume_of[last as usize]),
    ];
    for entry in std::fs::read_dir(fixture("nucl_mv.nal").trim_end_matches("nucl_mv.nal")).unwrap()
    {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_str().unwrap().to_string();
        if name.starts_with("nucl_mv.") && !missing.contains(&name) {
            std::fs::copy(&path, dir.join(&name)).unwrap();
        }
    }
    let name = dir.join("nucl_mv").to_str().unwrap().to_string();
    let database = BlastnDatabase::open(&name, &ApplicationSettings::default()).unwrap();
    let mut titles = DatabaseTitles::new();
    // Query 1 reports `middle`; ordinal id 0 (no header file) is past the count.
    let batch = vec![
        None,
        hit_list(vec![(middle, 1e-9), (first_volume_ordinal, 1e-3)], 1),
    ];
    read_titles(&database, &batch, &mut titles).unwrap();
    assert_eq!(titles.keys().copied().collect::<Vec<_>>(), [middle]);
    assert_eq!(
        titles[&middle].stored().title,
        entries()[middle as usize].2.as_bytes()
    );
    // The next batch: query 0 names `middle` again (not read again), query 1 names the
    // last sequence before ordinal id 0: the first header file read is the last
    // volume's.
    let batch = vec![
        hit_list(vec![(middle, 1e-9)], 1),
        hit_list(vec![(last, 1e-9), (first_volume_ordinal, 1e-5)], 2),
    ];
    let error = read_titles(&database, &batch, &mut titles)
        .expect_err("a missing header file")
        .downcast::<HeaderReadError>()
        .expect("NCBI's error, with the query that reached it");
    assert_eq!(error.query, 1);
    let native = error
        .error
        .downcast_ref::<NativeError>()
        .expect("NCBI's error");
    assert_eq!(native.exit, 2);
    // The volume is recorded under its directory normalised with symbolic links
    // followed (seqdbalias.cpp:851-862), as in BLASTP's test.
    assert_eq!(
        String::from_utf8(native.message.clone()).unwrap(),
        format!(
            "LOSAT database error: Error: File ({}) not found.\n",
            dir.canonicalize().unwrap().join(&missing[1]).display()
        )
    );
    drop(database);
    std::fs::remove_dir_all(&dir).unwrap();
}

// The reports name the reported sequences by rank (the order of the ordinal ids), with
// the stored title, the id `gnl|BL_ORD_ID|n` and the stored data, whose letters are the
// stored letters in upper case with their ambiguity codes (seqdbvol.cpp:1406-1414,
// blast_encoding.c:80-83), whole and in any range, and the tabular ids are the first
// words of the titles.
#[test]
fn reported_subjects_are_renamed_to_ranks_with_their_stored_records() {
    let database = open("nucl_mv");
    let entries = entries();
    let mut hit_lists = vec![
        hit_list(vec![(9, 1e-9), (5, 1e-5)], 2),
        None,
        hit_list(vec![(5, 1e-9), (11, 1e-5)], 2),
    ];
    let mut titles = DatabaseTitles::new();
    read_titles(&database, &hit_lists, &mut titles).unwrap();
    // A list left without HSPs (its HSPs had score 0) is not reported.
    hit_lists[2].as_mut().unwrap().hsplist_array[1].hsps.clear();
    let reported = reported_subjects(
        &database,
        &titles,
        &mut hit_lists,
        true,
        ParsedShown::default(),
        (crate::report::db_subject::FormatThrow::Unread, None),
        &[],
        0,
        true,
        None,
    )
    .unwrap();
    assert_eq!(reported.ordinals, [5, 9]);
    let ranks = |hit_list: &Option<BlastnHitList>| -> Vec<(u32, Vec<u32>)> {
        hit_list
            .as_ref()
            .unwrap()
            .hsplist_array
            .iter()
            .map(|list| (list.oid, list.hsps.iter().map(|hsp| hsp.s_idx).collect()))
            .collect()
    };
    assert_eq!(ranks(&hit_lists[0]), [(1, vec![1]), (0, vec![0])]);
    assert_eq!(ranks(&hit_lists[2]), [(0, vec![0]), (11, vec![])]);
    let letters = crate::algorithm::blastn::x_subject_store::SubjectLetters::new(
        &reported.records,
        reported.stores.as_deref(),
    );
    for (rank, &ordinal) in reported.ordinals.iter().enumerate() {
        let record = &reported.records[rank];
        assert_eq!(record.local_id, format!("gnl|BL_ORD_ID|{ordinal}"));
        assert_eq!(record.title, entries[ordinal as usize].2.as_bytes());
        assert!(record.seq().is_empty(), "the records keep no letters");
        let stored = entries[ordinal as usize].1.as_bytes();
        assert_eq!(letters.len(rank), stored.len());
        assert_eq!(letters.upper_letters(rank, 0, stored.len()), stored);
        let (from, to) = (stored.len() / 3, stored.len() - stored.len() / 4);
        assert_eq!(letters.upper_letters(rank, from, to), &stored[from..to]);
        assert!(letters.lowercase_masks(rank).is_empty());
    }
    // n10 was stored from lower case letters: the database has upper case only.
    assert!(letters
        .upper_letters(1, 0, letters.len(1))
        .iter()
        .all(|letter| !letter.is_ascii_lowercase()));
    assert_eq!(
        reported
            .tabular_subject_ids(database.db.oid_map(), true)
            .unwrap()
            .iter()
            .map(|id| id.as_deref())
            .collect::<Vec<_>>(),
        [Some(&b"n06"[..]), Some(b"n10")]
    );
    // Without the letters (the tabular formats without alignment rows): no stored data.
    let mut hit_lists = vec![hit_list(vec![(9, 1e-9)], 1)];
    let reported = reported_subjects(
        &database,
        &titles,
        &mut hit_lists,
        false,
        ParsedShown::default(),
        (crate::report::db_subject::FormatThrow::Unread, None),
        &[],
        0,
        true,
        None,
    )
    .unwrap();
    assert!(reported.records[0].seq().is_empty());
    assert!(reported.stores.is_none());
}

// NCBI writes an error message for a stored title byte outside ' ' ..= '~' each time
// its formatter reads the definition line (`check_reported_title`).
// (A path in the system's temporary directory: Unix spelling.)
#[cfg(unix)]
#[test]
fn a_reported_title_with_a_byte_outside_the_visible_range_is_rejected() {
    let dir = temp_dir("badtitle");
    for ext in ["nin", "nsq"] {
        std::fs::copy(fixture("nucl_v4") + "." + ext, dir.join(format!("b.{ext}"))).unwrap();
    }
    let mut header = std::fs::read(fixture("nucl_v4") + ".nhr").unwrap();
    let at = header
        .windows(13)
        .position(|window| window == b"n02 length 41")
        .expect("the title of n02");
    header[at + 5] = 0xE9;
    std::fs::write(dir.join("b.nhr"), header).unwrap();
    let name = dir.join("b").to_str().unwrap().to_string();
    let database = BlastnDatabase::open(&name, &ApplicationSettings::default()).unwrap();
    let mut titles = DatabaseTitles::new();
    let mut clean = vec![hit_list(vec![(0, 1e-9)], 1)];
    read_titles(&database, &clean, &mut titles).unwrap();
    reported_subjects(
        &database,
        &titles,
        &mut clean,
        false,
        ParsedShown::default(),
        (crate::report::db_subject::FormatThrow::Unread, None),
        &[],
        0,
        true,
        None,
    )
    .unwrap();
    let mut bad = vec![hit_list(vec![(1, 1e-9)], 1)];
    read_titles(&database, &bad, &mut titles).unwrap();
    assert_eq!(titles[&1].stored().title, b"n02 l#ngth 41");
    let error = reported_subjects(
        &database,
        &titles,
        &mut bad,
        false,
        ParsedShown::default(),
        (crate::report::db_subject::FormatThrow::Unread, None),
        &[],
        0,
        true,
        None,
    )
    .err()
    .expect("the bad byte");
    assert!(
        format!("{error:#}")
            .starts_with("database sequence 1 has a stored title with a byte outside ' ' to '~'"),
        "{error:#}"
    );
    assert!(format!("{error:#}").ends_with("not supported by LOSAT's BLASTN"));
    drop(database);
    std::fs::remove_dir_all(&dir).unwrap();
}

// The records made on a pool are those made one at a time, and the error is the first one
// in the order of the ordinal ids.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn reported_records_on_a_pool_are_the_serial_ones() {
    let ordinals: Vec<u32> = (0..500).collect();
    let record = |&ordinal: &u32| -> Result<Reported> {
        if ordinal >= 300 && ordinal % 2 == 0 {
            anyhow::bail!("refused {ordinal}");
        }
        let record = FastaRecord {
            local_id: format!("gnl|BL_ORD_ID|{ordinal}"),
            title: Vec::new(),
            sequence: Vec::new(),
            warnings: Vec::new(),
            source: Default::default(),
        };
        Ok((record, None, None, Default::default()))
    };
    let ids = |made: Vec<Reported>| -> Vec<String> {
        made.into_iter()
            .map(|(record, _, _, _)| record.local_id)
            .collect()
    };
    let serial = x_reported_records(&ordinals[..300], None, &record).unwrap();
    crate::utils::threading::with_search_pool(4, "BLASTN", |pool| {
        let parallel = x_reported_records(&ordinals[..300], Some(pool), &record).unwrap();
        assert_eq!(ids(parallel), ids(serial));
        let error = x_reported_records(&ordinals, Some(pool), &record).expect_err("refused");
        assert_eq!(format!("{error}"), "refused 300");
        Ok(())
    })
    .unwrap();
    let error = x_reported_records(&ordinals, None, &record).expect_err("refused");
    assert_eq!(format!("{error}"), "refused 300");
}

/// The fixture database `name` of `tests/db_search`, with the taxonomy database variant
/// `taxdb` of `tests/seqdb/taxdb` (`None`: none) as `BLASTDB`, from the directory `cwd`.
fn open_with_taxdb(name: &str, taxdb: Option<&str>, cwd: &std::path::Path) -> BlastnDatabase {
    use crate::seqdb::alias::SearchPath;
    let tests = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let taxdb = taxdb.map(|variant| tests.join("seqdb/taxdb").join(variant));
    let search =
        SearchPath::new(cwd, taxdb.as_deref().map(std::path::Path::as_os_str), None).unwrap();
    let path = tests.join("db_search").join(name);
    let mut db = SeqDb::open_in(path.to_str().unwrap(), SeqType::Nucleotide, &search).unwrap();
    db.read_name_records(path.to_str().unwrap(), false).unwrap();
    let parsed_seqids = db.has_parsed_seqids().unwrap();
    BlastnDatabase {
        metadata: SubjectMetadata {
            db_len_total: db.total_length() as i64,
            db_num_seqs: db.num_seqs() as i32,
        },
        db,
        name: name.to_string(),
        taxdb: shared::TaxdbRun::default(),
        loaded_before_batch: Default::default(),
        parsed_seqids,
    }
}

// DBID item T: a batch's prefetch looks up the sequences with HSPs in its final hit lists
// (within the hit list size), once per volume and taxid over the batches of a run, and opens
// the taxonomy database at the run's first lookup; nothing without `taxdb.bti`.
// `x16_ntaxid_99990` (two sequences of taxid 99990) with `offset_last_past_end` (the last
// entry, 99990, is past the end of `taxdb.btd`), `x16_ntaxid` (taxid 9606) with `bad_magic`.
#[test]
fn the_prefetch_looks_up_the_sequences_with_hsps_once_per_run() {
    let cwd = temp_dir("prefetch");
    let text = |database: &BlastnDatabase, hit_lists: &[Option<BlastnHitList>]| -> String {
        let titles = DatabaseTitles::new();
        String::from_utf8(
            database
                .prefetch_messages(&titles, hit_lists, &HitLoads::default())
                .unwrap(),
        )
        .unwrap()
    };
    let database = open_with_taxdb("x16_ntaxid_99990", Some("offset_last_past_end"), &cwd);
    assert!(database.taxdb_found());
    // Sequence 0 without HSPs (its HSPs of score 0 were dropped), sequence 1 past the hit
    // list size: nothing is fetched, the taxonomy database is not opened.
    let mut none = hit_list(vec![(0, 1e-9), (1, 1e-8)], 1);
    none.as_mut().unwrap().hsplist_array[0].hsps.clear();
    assert_eq!(text(&database, &[None, none]), "");
    let offset = "Error: [blastn] Error: Offset error at end of taxdb file.\n";
    assert_eq!(text(&database, &[hit_list(vec![(1, 1e-9)], 1)]), offset);
    // The same volume and taxid in a later batch.
    assert_eq!(
        text(&database, &[hit_list(vec![(0, 1e-9), (1, 1e-8)], 2)]),
        ""
    );
    let database = open_with_taxdb("x16_ntaxid", Some("bad_magic"), &cwd);
    assert_eq!(
        text(&database, &[hit_list(vec![(0, 1e-9), (1, 1e-8)], 2)]),
        "Error: [blastn] Error: Tax database file has wrong magic number.\n"
    );
    let database = open_with_taxdb("x16_ntaxid", None, &cwd);
    assert!(!database.taxdb_found());
    assert_eq!(
        text(&database, &[hit_list(vec![(0, 1e-9), (1, 1e-8)], 2)]),
        ""
    );
    std::fs::remove_dir_all(&cwd).unwrap();
}
