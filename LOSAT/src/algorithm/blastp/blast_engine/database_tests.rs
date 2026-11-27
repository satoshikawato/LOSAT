use super::*;
use crate::algorithm::blastp::hsp::{BlastpHsp, BlastpHspList};
use crate::config::{ProteinScoringSpec, ScoringMatrix};
#[cfg(unix)]
use crate::seqdb::SeqDbError;
use crate::stats::lookup_protein_params;

/// NCBIstdaa code of pyrrolysine (O).
const NCBISTDAA_O: u8 = 26;

fn fixture(name: &str) -> String {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/seqdb/db")
        .join(name)
        .to_str()
        .expect("fixture path is UTF-8")
        .to_string()
}

fn open(name: &str) -> BlastpDatabase {
    BlastpDatabase::open(fixture(name), &ApplicationSettings::default(), "BLASTP")
        .unwrap_or_else(|error| panic!("{name}: {error:#}"))
}

// The byte-by-byte test that `holds_non_residue_byte` makes.
fn holds_non_residue_byte_by_byte(bytes: &[u8]) -> bool {
    bytes
        .iter()
        .any(|&residue| usize::from(residue) >= BLASTAA_SIZE)
}

// Lengths across vectors of 16 and 32 bytes, gap bytes from none to all, O, the last residue
// (27), the first byte that is none (28) and 255: the flag is that of the byte-by-byte test.
#[test]
fn holds_non_residue_byte_gives_the_value_of_the_byte_by_byte_test() {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut flagged = 0usize;
    for length in [0usize, 1, 2, 223, 224, 225, 447, 448, 449, 1_000, 4_096] {
        for gap_every in [0u64, 1, 2, 7, 300] {
            for wrong in [None, Some(27u8), Some(28), Some(255)] {
                let bytes: Vec<u8> = (0..length)
                    .map(|_| {
                        let r = next();
                        if gap_every != 0 && r % gap_every == 0 {
                            NCBISTDAA_GAP
                        } else if r % 11 == 0 {
                            NCBISTDAA_O
                        } else if r % 97 == 0 {
                            wrong.unwrap_or(1)
                        } else {
                            1 + (r % 25) as u8
                        }
                    })
                    .collect();
                let expected = holds_non_residue_byte_by_byte(&bytes);
                flagged += usize::from(expected);
                assert_eq!(
                    holds_non_residue_byte(&bytes),
                    expected,
                    "length {length}, a gap every {gap_every}, {wrong:?}"
                );
            }
        }
    }
    assert!(flagged > 20, "cases with a byte of 28 or more: {flagged}");
    // Every byte a gap.
    assert!(!holds_non_residue_byte(&[NCBISTDAA_GAP; 1_000]));
    // A byte of 28 alone: second, after the first vectors, last of the input.
    for (length, at) in [(2usize, 1usize), (300, 224), (511, 510)] {
        let mut bytes = vec![1u8; length];
        bytes[at] = 28;
        assert!(holds_non_residue_byte(&bytes), "28 at {at} of {length}");
    }
    // O and J, the last two residues.
    let mut bytes = vec![1u8; 300];
    bytes[223] = NCBISTDAA_O;
    bytes[224] = 27;
    assert!(!holds_non_residue_byte(&bytes));
}

fn blosum62_gapped() -> KarlinParams {
    lookup_protein_params(&ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: 11,
        gap_extend: 1,
    })
}

// NCBI BLAST+ 2.17.0, `blastp -db sp_1k -query q1.faa -num_threads 1` (1,000 sequences,
// 357,308 letters; a query of 379 letters): `Effective search space used: 84653016`,
// and with `-dbsize 1000000000`: `245967282000` (the task's
// logs/db3_probe_cli1.log, cases dbsize.db.fmt0_none and dbsize.db.fmt0).
#[test]
fn search_space_is_the_one_of_ncbis_footer() {
    let params = blosum62_gapped();
    let (space, exact) =
        blastp_database_search_space(379, 357_308, 1_000, &params, "BLASTP").unwrap();
    assert_eq!(exact, 84_653_016);
    assert_eq!(space.effective_space, 84_653_016.0);
    let adj = space.length_adjustment;
    assert_eq!(
        space.effective_space,
        ((379 - adj) * (357_308 - 1_000 * adj)) as f64
    );
    let (space, exact) =
        blastp_database_search_space(379, 1_000_000_000, 1_000, &params, "BLASTP").unwrap();
    assert_eq!(exact, 245_967_282_000);
    assert_eq!(space.effective_space, 245_967_282_000.0);
    let adj = space.length_adjustment;
    assert_eq!(
        space.effective_space,
        ((379 - adj) * (1_000_000_000 - 1_000 * adj)) as f64
    );
}

// blast_setup.c:821-844: the adjustment is taken off the query once and off the
// database once per sequence.
#[test]
fn search_space_of_a_small_database_length() {
    let params = blosum62_gapped();
    // `-dbsize 1`: no length adjustment is possible (blast_stat.c: c < 0), so the
    // search space is the query length.
    let (space, exact) = blastp_database_search_space(379, 1, 1_000, &params, "BLASTP").unwrap();
    assert_eq!((space.length_adjustment, exact), (0, 379));
    // A total of five letters per sequence.
    let (space, exact) =
        blastp_database_search_space(379, 5_000, 1_000, &params, "BLASTP").unwrap();
    assert_eq!(space.length_adjustment, 4);
    assert_eq!(space.effective_space, 1_000.0 * 375.0);
    assert_eq!(exact, 1_000 * 375);
}

// NCBI's `Int8` product overflows for such a `-dbsize`; LOSAT rejects it.
#[test]
fn search_space_beyond_int8_is_rejected() {
    let error = blastp_database_search_space(
        379,
        9_000_000_000_000_000_000,
        1_000,
        &blosum62_gapped(),
        "BLASTP",
    )
    .unwrap_err();
    assert!(
        format!("{error:#}").contains("overflows NCBI BLAST+'s 64-bit effective search space"),
        "{error:#}"
    );
}

#[test]
fn database_numbers_come_from_the_index_files() {
    // LOSAT/tests/seqdb/expected/prot_mv.info.tsv (blastdbcmd -info).
    let database = open("prot_mv");
    assert_eq!(database.num_sequences(), 10);
    assert_eq!(database.total_letters(), 801);
    let shown = database.pairwise_database();
    assert_eq!(shown.definition, b"protein in volumes");
    assert_eq!(shown.databases[0].date, "Oct 12, 2026  4:12 PM");
    assert_eq!((shown.num_sequences, shown.total_letters), (10, 801));
    // The name is the `-db` value as given.
    assert_eq!(database.name(), fixture("prot_mv").as_bytes());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:904-911
// An alias file's `LENGTH` of 2^32 letters or more over the 23 sequences of
// `tests/db_search/taxfilter/tfp` (`tests/db_search/dbcounts`, rows db.cnt.l4295.bp.*,
// db.cnt.l4294.bp.0 and db.cnt.avg2p32.bp_c0) is NCBI's `Int8` on every target, wasm32's
// 32-bit `usize` included.
#[test]
fn totals_of_2_to_the_32_letters_or_more_are_held_whole() {
    for (name, total) in [
        ("l4295", 4_294_967_296_000_000_000),
        ("l4294", 4_294_967_295_999_999_999),
        ("avg2p32", 98_784_247_808),
    ] {
        assert!(total > i64::from(u32::MAX), "{name}");
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/db_search/dbcounts")
            .join(name);
        let database = BlastpDatabase::open(
            path.to_str().unwrap(),
            &ApplicationSettings::default(),
            "BLASTP",
        )
        .unwrap_or_else(|error| panic!("{name}: {error:#}"));
        assert_eq!(
            (
                database.total_letters(),
                database.statistics_num_sequences()
            ),
            (total, 23),
            "{name}"
        );
        let shown = database.pairwise_database();
        assert_eq!(
            (shown.num_sequences, shown.total_letters),
            (23, total),
            "{name}"
        );
    }
}

// align_format_util.cpp:385-387: an empty title is shown as the name.
// (A path in the system's temporary directory or in an error text: Unix spelling.)
#[cfg(unix)]
#[test]
fn a_database_without_a_title_shows_its_name() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-db-notitle-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    // The fixture's index with an empty title: version, type, title length 0, then
    // the rest of the file after the stored title.
    let index = std::fs::read(fixture("prot_v4") + ".pin").unwrap();
    let title_len = u32::from_be_bytes([index[8], index[9], index[10], index[11]]) as usize;
    let mut untitled = index[..8].to_vec();
    untitled.extend_from_slice(&0u32.to_be_bytes());
    untitled.extend_from_slice(&index[12 + title_len..]);
    // The offsets of the index point into the header and sequence files, not into
    // the index itself, so they stay valid.
    std::fs::write(dir.join("u.pin"), untitled).unwrap();
    std::fs::copy(fixture("prot_v4") + ".psq", dir.join("u.psq")).unwrap();
    std::fs::copy(fixture("prot_v4") + ".phr", dir.join("u.phr")).unwrap();
    let name = dir.join("u").to_str().unwrap().to_string();
    let database = BlastpDatabase::open(&name, &ApplicationSettings::default(), "BLASTP").unwrap();
    assert_eq!(database.pairwise_database().definition, name.as_bytes());
    drop(database);
    std::fs::remove_dir_all(&dir).unwrap();
}

// The first header of each volume is read when the database is opened (`shared::
// check_seqids`), as before line DBID ID2a: a header that LOSAT does not read is refused
// there, and only one with parsed Seq-ids is let through (review 1 of ID2a, M2). A header
// that NCBI's read fails on is not refused: NCBI reads no header at the opening, and fails
// where a report reads it (RG-T5d G3); the volume is classified by its next header.
#[cfg(unix)]
#[test]
fn the_first_header_of_each_volume_is_read_at_open() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-db-badheader-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let header = std::fs::read(fixture("prot_v4") + ".phr").unwrap();
    std::fs::copy(fixture("prot_v4") + ".pin", dir.join("b.pin")).unwrap();
    std::fs::copy(fixture("prot_v4") + ".psq", dir.join("b.psq")).unwrap();
    let name = dir.join("b").to_str().unwrap().to_string();
    let settings = ApplicationSettings::default();
    // NCBI's `Unexpected tag: 17, expected: 16` at byte 0.
    let mut set_of = header.clone();
    set_of[0] = 0x31;
    std::fs::write(dir.join("b.phr"), &set_of).unwrap();
    for program in ["BLASTP", "PSIBLAST"] {
        assert!(
            BlastpDatabase::open(&name, &settings, program).is_ok(),
            "{program}"
        );
    }
    // A UTF8String title, which NCBI's reader takes with its warning once per process
    // (RG-T5d fix 1): the database opens, the header reads, and the reader names the member
    // of the warning, which the run posts once.
    let mut utf8 = header;
    assert_eq!(utf8[6], 0x1a);
    utf8[6] = 0x0c;
    std::fs::write(dir.join("b.phr"), &utf8).unwrap();
    for program in ["BLASTP", "PSIBLAST"] {
        let database = BlastpDatabase::open(&name, &settings, program)
            .unwrap_or_else(|error| panic!("{program}: {error:#}"));
        let headers = shared::HeaderReader::new(&database.db, false);
        shared::read_header(&headers, 0, 3, program).unwrap();
        let read = headers.first_utf8().expect("the UTF8String read").clone();
        assert_eq!(
            read,
            shared::Utf8Read {
                query: 3,
                path: "Blast-def-line-set.E.title".to_string()
            }
        );
        let posts = read.post(&database.db, "blastp").expect("the first post");
        assert_eq!(
            posts.text(),
            b"Warning: (803.10) [blastp] CObjectIStreamAsnBinary: UTF8String data for VisibleString member Blast-def-line-set.E.title, ASN.1 specification may need an update\n"
        );
        assert!(read.post(&database.db, "blastp").is_none());
        // Once posted, a reader names no member again.
        let again = shared::HeaderReader::new(&database.db, false);
        shared::read_header(&again, 0, 0, program).unwrap();
        assert!(again.first_utf8().is_none());
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// blast_app_util.hpp:185-203 and seqdb.cpp:155-165.
// (A path in the system's temporary directory or in an error text: Unix spelling.)
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
    let Err(error) = BlastpDatabase::open("", &settings, "BLASTP") else {
        panic!("an empty name opened")
    };
    assert_eq!(
        native(error),
        (
            2,
            "LOSAT database error: Database name is required.\n".to_string()
        )
    );
    let Err(error) = BlastpDatabase::open(fixture("no_such_database"), &settings, "BLASTP") else {
        panic!("a missing database opened")
    };
    let (exit, message) = native(error);
    assert_eq!(exit, 2);
    assert!(
        message.starts_with(
            "LOSAT database error: No alias or index file found for protein database ["
        ),
        "{message}"
    );
    // A database with parsed Seq-ids opens for BLASTP and psiblast.
    for program in ["BLASTP", "PSIBLAST"] {
        BlastpDatabase::open(fixture("prot_ids_v4"), &settings, program)
            .unwrap_or_else(|e| panic!("{program} opens a database with parsed Seq-ids: {e:#}"));
    }
}

// local_blast.cpp:240-251 and blast_app_util.hpp:185-203: a reader error of NCBI's own
// is the search stage's text during the preliminary search (`BLAST engine error`,
// exit 3) and the reader's text elsewhere (`BLAST Database error`, exit 2); what the
// reader does not read is LOSAT's refusal in both places.
#[cfg(unix)]
#[test]
fn reader_errors_of_the_search_and_of_the_other_reads() {
    let unreadable = || SeqDbError::Open {
        path: "/data/db.psq".into(),
        source: std::io::Error::from_raw_os_error(13),
    };
    let read = shared::search_read_error(unreadable(), "BLASTP")
        .downcast::<PreliminaryReadError>()
        .expect("NCBI's error of the search stage");
    assert_eq!(read.as_bytes(), b"Database memory map file error");
    assert_eq!(
        read.to_string(),
        "LOSAT engine error: Database memory map file error"
    );
    let native = shared::database_error(unreadable(), "BLASTP");
    let native = native.downcast_ref::<NativeError>().expect("NCBI's error");
    assert_eq!(native.exit, 2);
    assert_eq!(
        native.message,
        b"LOSAT database error: Database memory map file error - CMemoryFile: Cannot memory map file '/data/db.psq':Permission denied\n"
    );
    let too_many = || SeqDbError::Open {
        path: "/data/db.psq".into(),
        source: std::io::Error::from_raw_os_error(24),
    };
    let read = shared::search_read_error(too_many(), "BLASTP")
        .downcast::<PreliminaryReadError>()
        .expect("NCBI's error of the search stage");
    assert_eq!(
        read.as_bytes(),
        b"Too many open files, please raise the open file limit"
    );
    let missing = SeqDbError::NotFound("Error: File (/data/db.psq) not found.".to_string().into());
    let read = shared::search_read_error(missing, "BLASTP")
        .downcast::<PreliminaryReadError>()
        .expect("NCBI's error of the search stage");
    assert_eq!(read.as_bytes(), b"Database memory map file error");
    for refused in [shared::search_read_error, shared::database_error] {
        let error = refused(
            SeqDbError::Unsupported {
                path: "/data/db.pal".into(),
                message: "the alias file filters sequences".to_string(),
            },
            "BLASTP",
        );
        assert!(error.downcast_ref::<NativeError>().is_none());
        assert!(error.downcast_ref::<PreliminaryReadError>().is_none());
        assert!(
            format!("{error:#}").ends_with("; this database is not supported by LOSAT's BLASTP"),
            "{error:#}"
        );
    }
}

// blast_aalookup.c:494-508: NCBI searches a stored O as it is (the lookup meets it as a
// neighbour only, the matrix scores it as X). p03 of the fixture holds every NCBIstdaa
// letter but the gap (`ABCDEFGHIKLMNPQRSTVWXYZUOJMHCENQT`).
#[test]
fn a_search_block_holds_the_stored_letters_o_included() {
    let database = open("prot_v4");
    let range = database.db.next_block(0, u64::MAX).unwrap();
    assert_eq!((range.first_ordinal, range.len()), (0, 10));
    let stored = database.db.read_block(range).unwrap();
    let block = read_search_block(&database.db, range, "BLASTP").unwrap();
    for i in 0..block.len() {
        let expected = stored.protein(i);
        assert_eq!(block.protein(i), expected, "sequence {i}");
        let subject = block_subject(&block, i);
        assert_eq!(subject.aa_len, expected.len());
        assert_eq!(subject.aa_seq.len(), expected.len() + 2);
        assert_eq!(
            (subject.aa_seq[0], subject.aa_seq[expected.len() + 1]),
            (0, 0)
        );
        assert_eq!(&subject.aa_seq[1..=expected.len()], expected);
    }
    let letters: String = block
        .protein(2)
        .iter()
        .map(|&residue| ncbistdaa_to_ascii(residue))
        .collect();
    assert_eq!(letters, "ABCDEFGHIKLMNPQRSTVWXYZUOJMHCENQT");
    assert!(block.protein(2).contains(&NCBISTDAA_O));
}

#[test]
fn a_search_block_with_a_foreign_byte_is_rejected_and_one_with_a_gap_residue_is_not() {
    let database = open("prot_v4");
    let range = database.db.next_block(0, u64::MAX).unwrap();
    // The gap character in sequence 3 (its third residue), and as the whole of sequence 5:
    // residues that NCBI searches (row and column 0 of the matrix).
    let mut block = database.db.read_block(range).unwrap();
    let at = block.protein_bytes_mut().len() - block.protein_bytes_mut().len(); // 0: the leading sentinel
    assert_eq!(block.protein_bytes_mut()[at], 0);
    let offset_of = |block: &SequenceBlock, i: usize| -> usize {
        (0..i).map(|j| block.protein(j).len() + 1).sum::<usize>() + 1
    };
    let position = offset_of(&block, 3) + 2;
    block.protein_bytes_mut()[position] = 0;
    let fifth = offset_of(&block, 5);
    let fifth_len = block.protein(5).len();
    block.protein_bytes_mut()[fifth..fifth + fifth_len].fill(0);
    prepare_search_block(&block, database.db.oid_map(), "BLASTP").unwrap();
    assert_eq!(block.protein(3)[2], 0);
    assert!(block.protein(5).iter().all(|&residue| residue == 0));
    // A byte that is no NCBIstdaa residue in sequence 9.
    let mut block = database.db.read_block(range).unwrap();
    let position = offset_of(&block, 9);
    block.protein_bytes_mut()[position] = 28;
    let error = prepare_search_block(&block, database.db.oid_map(), "BLASTP").unwrap_err();
    assert!(
        format!("{error:#}").starts_with(
            "database sequence 9 holds a byte that is no NCBIstdaa residue at position 1,"
        ),
        "{error:#}"
    );
    // A sentinel that is not NUL (the byte after sequence 0), alone and with a NUL
    // inside sequence 9, so that the count of NUL bytes alone would not show it.
    for compensated in [false, true] {
        let mut block = database.db.read_block(range).unwrap();
        let after_first = offset_of(&block, 1) - 1;
        block.protein_bytes_mut()[after_first] = 1;
        if compensated {
            let inside = offset_of(&block, 9) + 5;
            block.protein_bytes_mut()[inside] = 0;
        }
        let error = prepare_search_block(&block, database.db.oid_map(), "BLASTP").unwrap_err();
        assert!(
            format!("{error:#}").starts_with("database sequence 0 is not between two NUL bytes"),
            "{error:#}"
        );
    }
}

#[test]
fn a_rejected_search_block_of_a_restricted_database_names_the_oid() {
    // prot_v4 through an alias whose range keeps OIDs 2 to 9 (line DBID ordnum): the block's
    // sequence 7 (rank 7) is OID 9.
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat_bp_ordnum_{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let alias = dir.join("prot_v4_r.pal");
    std::fs::write(
        &alias,
        format!(
            "TITLE prot_v4_r\nDBLIST {}\nFIRST_OID 3\nLAST_OID 10\n",
            fixture("prot_v4")
        ),
    )
    .unwrap();
    let database = BlastpDatabase::open(
        dir.join("prot_v4_r").to_str().unwrap(),
        &ApplicationSettings::default(),
        "BLASTP",
    )
    .unwrap();
    let range = database.db.next_block(0, u64::MAX).unwrap();
    assert_eq!((range.first_ordinal, range.len()), (0, 8));
    let mut block = database.db.read_block(range).unwrap();
    let position = (0..7).map(|j| block.protein(j).len() + 1).sum::<usize>() + 1;
    block.protein_bytes_mut()[position] = 28;
    let error = prepare_search_block(&block, database.db.oid_map(), "BLASTP").unwrap_err();
    assert!(
        format!("{error:#}").starts_with(
            "database sequence 9 holds a byte that is no NCBIstdaa residue at position 1,"
        ),
        "{error:#}"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

/// A kernel that records what it is given and "finds" every third sequence.
fn recording_kernel(
    window: i32,
    seen: &std::sync::Mutex<Vec<(usize, i32, Vec<u8>)>>,
) -> impl Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync + '_
{
    move |scratch, s_idx, subject, result| {
        result.reset(s_idx);
        assert_eq!(subject.aa_seq.len(), subject.aa_len + 2);
        seen.lock()
            .unwrap()
            .push((s_idx, scratch.diag_offset, subject.aa_seq.to_vec()));
        scratch.finish_subject(subject.aa_len, window);
        if s_idx % 3 == 0 {
            result.x_compact_hsp_lists.push((
                0,
                BlastpHspList {
                    oid: s_idx as u32,
                    query_index: 0,
                    hsps: vec![test_prelim_hsp(s_idx as u32, 50 + (s_idx as i32 % 2))],
                    best_evalue: 1e-5,
                },
            ));
        }
    }
}

fn test_prelim_hsp(s_idx: u32, raw_score: i32) -> XPrelimHsp {
    use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
    XPrelimHsp {
        hsp: BlastpPreliminaryHsp {
            query_start: 0,
            query_end: 10,
            subject_start: 0,
            subject_end: 10,
            gapped_query_start: 1,
            gapped_subject_start: 1,
            raw_score,
            query_context: 0,
            query_frame: 0,
            subject_frame: 0,
            query_length: 20,
            q_idx: 0,
            s_idx,
        },
        e_value: 1e-5,
    }
}

fn test_hsp(s_idx: u32, raw_score: i32) -> BlastpHsp {
    BlastpHsp {
        num: 1,
        identity: 100.0,
        length: 10,
        mismatch: 0,
        gapopen: 0,
        q_start: 1,
        q_end: 10,
        s_start: 1,
        s_end: 10,
        gapped_q_start: 0,
        gapped_s_start: 0,
        e_value: 1e-5,
        bit_score: 30.0,
        num_ident: 10,
        query_context: 0,
        query_frame: 0,
        subject_frame: 0,
        query_length: 20,
        q_idx: 0,
        s_idx,
        raw_score,
        gap_info: None,
        num_positives: 10,
        comp_adjustment_method: 0,
    }
}

fn run_blocks(
    database: &BlastpDatabase,
    threads: usize,
    block_bytes: u64,
    window: i32,
) -> (Vec<(usize, i32, Vec<u8>)>, Vec<u32>) {
    let seen = std::sync::Mutex::new(Vec::new());
    let mut hit_lists: Vec<Option<BlastpHitList<XPrelimHsp>>> = vec![None];
    crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
        search_blocks_of(
            database,
            pool,
            window,
            block_bytes,
            2,
            &|| BlastpSubjectScratch::new(16, 64, window, 21),
            &recording_kernel(window, &seen),
            &mut crate::algorithm::blastp::blast_engine::PrelimCollector {
                x_compact_hit_lists: &mut hit_lists,
                prelim_hitlist_size: 500,
            },
        )
    })
    .unwrap();
    let mut seen = seen.into_inner().unwrap();
    seen.sort_by_key(|entry| entry.0);
    let listed = hit_lists[0]
        .as_ref()
        .map(|list| {
            list.hsplist_array[..list.hsplist_count]
                .iter()
                .map(|l| l.oid)
                .collect()
        })
        .unwrap_or_default();
    (seen, listed)
}

// blast_extend.c:160-172 (the diagonal offset moves by the subject length and the
// window after every subject) and blast_engine.c:1411-1554 (the subjects in the order
// of their ordinal ids): the block driver gives the kernel every sequence once, with
// the diagonal offset of the one-thread order, whatever the block size and the thread
// count, across the six volumes of the fixture.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn blocks_give_every_sequence_the_diagonal_offset_of_the_one_thread_order() {
    let database = open("prot_mv");
    let single = open("prot_v4");
    let window = 40;
    // LOSAT/tests/seqdb/expected/prot_v4.entries.tsv: the lengths in ordinal order.
    let lengths = [61, 1, 33, 46, 130, 88, 45, 46, 51, 300];
    let mut expected_offsets = Vec::new();
    let mut offset = window;
    for length in lengths {
        expected_offsets.push(offset);
        offset += length + window;
    }
    let (reference, reference_listed) = run_blocks(&database, 1, u64::MAX, window);
    assert_eq!(reference.len(), 10);
    for (ordinal, (s_idx, diag_offset, aa_seq)) in reference.iter().enumerate() {
        assert_eq!(*s_idx, ordinal);
        assert_eq!(
            *diag_offset, expected_offsets[ordinal],
            "sequence {ordinal}"
        );
        assert_eq!(aa_seq.len() as i32, lengths[ordinal] + 2);
    }
    assert_eq!(reference_listed.len(), 4);
    // The same sequences from the single-volume fixture.
    assert_eq!(run_blocks(&single, 1, u64::MAX, window).0, reference);
    for threads in [1usize, 2, 4] {
        for block_bytes in [1u64, 64, 200, 1 << 20] {
            for db in [&database, &single] {
                let (seen, listed) = run_blocks(db, threads, block_bytes, window);
                assert_eq!(
                    seen, reference,
                    "{threads} threads, blocks of {block_bytes}"
                );
                assert_eq!(
                    listed, reference_listed,
                    "{threads} threads, blocks of {block_bytes}"
                );
            }
        }
    }
}

// blast_extend.c:168-170: at INT4_MAX / 4 the offset returns to the window.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn the_offset_recurrence_of_a_block_is_the_one_of_the_scratch() {
    let database = open("prot_v4");
    let range = database.db.next_block(0, u64::MAX).unwrap();
    let lengths = [61usize, 1, 33, 46, 130, 88, 45, 46, 51, 300];
    for (start, window) in [
        (40, 40),
        (0, 0),
        (i32::MAX / 4 - 100, 40),
        (i32::MAX / 4, 11),
    ] {
        let mut scratch = BlastpSubjectScratch::new(16, 64, window, 21);
        scratch.diag_offset = start;
        for length in lengths {
            scratch.finish_subject(length, window);
        }
        assert_eq!(
            diag_offset_after(&database.db, range, start, window),
            scratch.diag_offset,
            "from {start} with window {window}"
        );
    }
}

/// EXPERIMENT (LOSAT_X_BPDBCUT): the blocks cut by the database size and the thread count
/// (`database/x_db_cut.rs`) against the blocks of `BLOCK_BYTES`.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod x_db_cut_blocks {
    use super::super::x_db_cut::{
        x_count_blocks, x_cut_block_bytes, x_db_cut_block_bytes, x_search_blocks_cut,
        x_sequence_bytes,
    };
    use super::*;
    use std::path::{Path, PathBuf};

    type Seen = Vec<(usize, i32, Vec<u8>)>;

    /// What the recording kernel saw (sorted by ordinal id) and the listed ordinal ids of a
    /// search of `database` at `threads` threads in rounds of `blocks_per_thread` blocks per
    /// thread: `search_blocks_of` with blocks of `block_bytes` (`Some`), or the cut's driver
    /// (`None`).
    fn search(
        database: &BlastpDatabase,
        threads: usize,
        blocks_per_thread: usize,
        block_bytes: Option<u64>,
    ) -> Result<(Seen, Vec<u32>)> {
        let driver = match block_bytes {
            None => Driver::Cut,
            Some(block_bytes) => Driver::Of(block_bytes),
        };
        search_with(database, threads, blocks_per_thread, driver)
    }

    /// The drivers of the preliminary stage that the tests call.
    #[derive(Clone, Copy, Debug)]
    enum Driver {
        /// `x_search_blocks_cut`.
        Cut,
        /// `search_blocks_of` with blocks of this many bytes.
        Of(u64),
        /// `x_db_keep::x_search_blocks_kept` with blocks of this many bytes.
        Kept(u64),
        /// `search_blocks_made` (without LOSAT_X_BPDBKEEP) with blocks of this many bytes.
        Made(u64),
        /// `x_db_keep::x_search_window` with blocks of this many bytes, the window and the hook.
        Window(u64, super::super::x_db_keep::XWindowSize, Option<XHook>),
    }

    /// What the kernel does before it records a sequence (`Driver::Window`).
    #[derive(Clone, Copy, Debug)]
    enum XHook {
        /// Sleeps 50 ms at this ordinal id.
        SleepAt(usize),
        /// Panics at this ordinal id.
        PanicAt(usize),
    }

    /// `kernel` after the hook `hook`.
    fn hooked_kernel<'a>(
        kernel: impl Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult)
            + Sync
            + 'a,
        hook: Option<XHook>,
    ) -> impl Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult)
           + Sync
           + 'a {
        move |scratch, s_idx, subject, result| {
            match hook {
                Some(XHook::SleepAt(ordinal)) if ordinal == s_idx => {
                    std::thread::sleep(std::time::Duration::from_millis(50));
                }
                Some(XHook::PanicAt(ordinal)) if ordinal == s_idx => {
                    panic!("the kernel's hook at {ordinal}");
                }
                _ => {}
            }
            kernel(scratch, s_idx, subject, result);
        }
    }

    /// `search` with the driver `driver`.
    fn search_with(
        database: &BlastpDatabase,
        threads: usize,
        blocks_per_thread: usize,
        driver: Driver,
    ) -> Result<(Seen, Vec<u32>)> {
        let window = 40;
        let seen = std::sync::Mutex::new(Vec::new());
        let mut hit_lists: Vec<Option<BlastpHitList<XPrelimHsp>>> = vec![None];
        crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
            let new_scratch = || BlastpSubjectScratch::new(16, 64, window, 21);
            let hook = match driver {
                Driver::Window(_, _, hook) => hook,
                _ => None,
            };
            let kernel = hooked_kernel(recording_kernel(window, &seen), hook);
            let mut sink = crate::algorithm::blastp::blast_engine::PrelimCollector {
                x_compact_hit_lists: &mut hit_lists,
                prelim_hitlist_size: 500,
            };
            match driver {
                Driver::Cut => x_search_blocks_cut(
                    database,
                    pool,
                    window,
                    blocks_per_thread,
                    &new_scratch,
                    &kernel,
                    &mut sink,
                ),
                Driver::Of(block_bytes) => search_blocks_of(
                    database,
                    pool,
                    window,
                    block_bytes,
                    blocks_per_thread,
                    &new_scratch,
                    &kernel,
                    &mut sink,
                ),
                Driver::Kept(block_bytes) => super::super::x_db_keep::x_search_blocks_kept(
                    database,
                    pool,
                    window,
                    block_bytes,
                    blocks_per_thread,
                    &new_scratch,
                    &kernel,
                    &mut sink,
                ),
                Driver::Made(block_bytes) => search_blocks_made(
                    database,
                    pool,
                    window,
                    block_bytes,
                    blocks_per_thread,
                    &new_scratch,
                    &kernel,
                    &mut sink,
                ),
                Driver::Window(block_bytes, size, _) => super::super::x_db_keep::x_search_window(
                    database,
                    pool,
                    window,
                    block_bytes,
                    size,
                    &new_scratch,
                    &kernel,
                    &mut sink,
                ),
            }
        })?;
        let mut seen = seen.into_inner().unwrap();
        seen.sort_by_key(|entry| entry.0);
        let listed = hit_lists[0]
            .as_ref()
            .map(|list| {
                list.hsplist_array[..list.hsplist_count]
                    .iter()
                    .map(|l| l.oid)
                    .collect()
            })
            .unwrap_or_default();
        Ok((seen, listed))
    }

    /// The cut of `database` on a pool of `threads` threads.
    fn cut_of(database: &BlastpDatabase, threads: usize) -> u64 {
        crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
            Ok(x_db_cut_block_bytes(database, pool))
        })
        .unwrap()
    }

    /// The letters of the generated database: 300 sequences of 1 to 1,500 letters from a fixed
    /// generator, with sequences of 40,000 and 70,000 letters at ordinal ids 100 and 250.
    fn generated_lengths() -> Vec<usize> {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        (0..300)
            .map(|ordinal| match ordinal {
                100 => 40_000,
                250 => 70_000,
                _ => {
                    state = state
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1_442_695_040_888_963_407);
                    1 + (state >> 33) as usize % 1_500
                }
            })
            .collect()
    }

    /// The generated database written by `losat makeblastdb` into a new directory of the test
    /// (`label`), with volumes of at most `max_file_sz` bytes; returns the directory and the
    /// database name.
    fn generated_database(label: &str, max_file_sz: &str) -> (PathBuf, String) {
        const LETTERS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
        let mut fasta = Vec::new();
        let mut state = 7u64;
        for (ordinal, length) in generated_lengths().into_iter().enumerate() {
            fasta.extend_from_slice(format!(">s{ordinal} generated\n").as_bytes());
            for i in 0..length {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                fasta.push(LETTERS[(state >> 33) as usize % LETTERS.len()]);
                if i % 80 == 79 || i + 1 == length {
                    fasta.push(b'\n');
                }
            }
        }
        let dir = crate::utils::test_dir::temp_dir().join(format!(
            "losat-bpdbcut-{label}-{}",
            crate::utils::test_dir::process_id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let name = dir.join("gen").to_str().unwrap().to_string();
        let args = crate::writedb::makeblastdb::MakedbArgs {
            input: std::ffi::OsString::from("-"),
            input_type: "fasta".to_string(),
            dbtype: "prot".to_string(),
            title: Some(std::ffi::OsString::from("gen")),
            out: Some(std::ffi::OsString::from(&name)),
            blastdb_version: 4,
            max_file_sz: std::ffi::OsString::from(max_file_sz),
            metadata_output_prefix: None,
            oid_masks: None,
            logfile: None,
            parse_seqids: false,
            hash_index: false,
            taxid: None,
            taxid_map: None,
            mask_data: None,
            mask_id: None,
            mask_desc: None,
            gi_mask: false,
            gi_mask_name: None,
            randomize: false,
        };
        let now = crate::writedb::time::LocalTime {
            year: 2026,
            month: 10,
            day: 16,
            hour: 13,
            minute: 30,
            second: 0,
            unix_seconds: 1_792_157_400,
        };
        let (mut log, mut messages) = (Vec::new(), Vec::new());
        let status = crate::writedb::makeblastdb::run_with(
            &args,
            &crate::blastinput::ncbi_environment::ApplicationSettings::default(),
            Some(&now),
            crate::writedb::Platform::host(),
            &mut fasta.as_slice(),
            &mut log,
            &mut messages,
        );
        assert_eq!(
            status,
            0,
            "makeblastdb: {}",
            String::from_utf8_lossy(&messages)
        );
        (dir, name)
    }

    /// EXPERIMENT (LOSAT_X_BPDBKEEP): `cut` with the switch's largest block.
    fn x_keep_cut(cut: u64) -> u64 {
        if crate::seqdb::x_bp_db_keep() {
            cut.min(super::super::x_db_keep::X_KEEP_BLOCK_BYTES)
        } else {
            cut
        }
    }

    fn open_path(name: &str) -> BlastpDatabase {
        BlastpDatabase::open(name, &ApplicationSettings::default(), "BLASTP")
            .unwrap_or_else(|error| panic!("{name}: {error:#}"))
    }

    // The cut of the fixture databases (812 bytes of sequence, one volume or six): one thread
    // keeps 1 MiB, more threads the 4 KiB floor, and every volume is one block either way. The
    // sequence bytes are those of the one-volume sequence file.
    #[test]
    fn the_cut_of_tiny_and_multi_volume_databases() {
        let single = open("prot_v4");
        let volumes = open("prot_mv");
        let psq_bytes = |names: &[&str]| -> u64 {
            names
                .iter()
                .map(|name| std::fs::metadata(fixture(name)).unwrap().len())
                .sum()
        };
        assert_eq!(x_sequence_bytes(&single), psq_bytes(&["prot_v4.psq"]));
        assert_eq!(x_sequence_bytes(&single), 812);
        let mv: Vec<String> = (0..6).map(|v| format!("prot_mv.{v:02}.psq")).collect();
        let mv: Vec<&str> = mv.iter().map(String::as_str).collect();
        // Each further volume has its own leading sentinel byte.
        assert_eq!(x_sequence_bytes(&volumes) + 5, psq_bytes(&mv));
        assert_eq!(x_sequence_bytes(&volumes), x_sequence_bytes(&single));
        assert_eq!(cut_of(&single, 1), x_keep_cut(BLOCK_BYTES));
        assert_eq!(cut_of(&volumes, 1), x_keep_cut(BLOCK_BYTES));
        for threads in [2, 3, 4, 8] {
            assert_eq!(cut_of(&single, threads), x_cut_block_bytes(812, threads));
            assert_eq!(cut_of(&single, threads), 4096);
            assert_eq!(cut_of(&volumes, threads), 4096);
        }
        for block_bytes in [BLOCK_BYTES, 4096] {
            assert_eq!(x_count_blocks(&single.db, block_bytes), 1);
            assert_eq!(x_count_blocks(&volumes.db, block_bytes), 6);
        }
    }

    // blast_engine.c:1411-1554 and blast_extend.c:160-172: with the cut, the kernel gets every
    // sequence once with the one-thread order's diagonal offset, and the hit list gets the
    // one-thread order's updates, at 1, 2, 3 and 7 threads, in rounds of 2 and of
    // `BLOCKS_PER_THREAD` blocks per thread, on one volume and on several; the sequences of
    // 40,000 and 70,000 letters are larger than the cut at 2 threads and more (blocks of their
    // own).
    #[test]
    fn blocks_of_the_cut_give_the_one_thread_search() {
        let (single_dir, single_name) = generated_database("one", "3GB");
        let (multi_dir, multi_name) = generated_database("multi", "100000");
        let single = open_path(&single_name);
        let multi = open_path(&multi_name);
        let lengths = generated_lengths();
        let letters: usize = lengths.iter().sum();
        assert_eq!(single.total_letters(), letters as i64);
        assert_eq!(multi.total_letters(), letters as i64);
        assert_eq!(single.db.volumes().len(), 1);
        assert!(multi.db.volumes().len() > 2, "{}", multi.db.volumes().len());
        let (reference, reference_listed) = search(&single, 1, 2, Some(BLOCK_BYTES)).unwrap();
        assert_eq!(reference.len(), lengths.len());
        assert_eq!(reference_listed.len(), 100);
        assert_eq!(x_count_blocks(&single.db, BLOCK_BYTES), 1);
        for threads in [1usize, 2, 3, 7] {
            let cut = cut_of(&single, threads);
            if threads == 1 {
                assert_eq!(cut, x_keep_cut(BLOCK_BYTES));
            } else {
                assert_eq!(
                    cut,
                    x_cut_block_bytes(letters as u64 + lengths.len() as u64 + 1, threads)
                );
                // Both long sequences at 3 threads and more, the longer one at 2.
                assert!(cut < 70_000 + 2, "{threads} threads: {cut}");
                assert!(
                    x_count_blocks(&single.db, cut) > 2 * threads as u64,
                    "{threads} threads: {} blocks",
                    x_count_blocks(&single.db, cut)
                );
            }
            for blocks_per_thread in [2, BLOCKS_PER_THREAD] {
                for (what, db) in [("one volume", &single), ("volumes", &multi)] {
                    let (seen, listed) = search(db, threads, blocks_per_thread, None).unwrap();
                    let context =
                        format!("{what}, {threads} threads, {blocks_per_thread} blocks per thread");
                    assert!(seen == reference, "{context}");
                    assert_eq!(listed, reference_listed, "{context}");
                }
            }
        }
        let _ = std::fs::remove_dir_all(single_dir);
        let _ = std::fs::remove_dir_all(multi_dir);
    }

    // No NCBI counterpart: LOSAT_X_BPDBKEEP keeps a scratch and a block buffer per worker
    // (`x_db_keep.rs`); the kernel sees the same sequences from the same diagonal offsets, and
    // the hit lists receive the same updates, as with a scratch and a buffer per task.
    // At 1, 2, 3 and 7 threads, blocks of 1 MiB and of 4 KiB (the 40,000 and 70,000-letter
    // sequences in blocks of their own, so buffers grow and are refilled shorter), in rounds of 2 and of
    // `BLOCKS_PER_THREAD` blocks per thread, on one volume and on several.
    #[test]
    fn kept_scratches_and_buffers_give_the_search_without_them() {
        let (single_dir, single_name) = generated_database("keep1", "3GB");
        let (multi_dir, multi_name) = generated_database("keepn", "100000");
        let single = open_path(&single_name);
        let multi = open_path(&multi_name);
        let (reference, reference_listed) =
            search_with(&single, 1, 2, Driver::Made(BLOCK_BYTES)).unwrap();
        assert_eq!(reference.len(), generated_lengths().len());
        for threads in [1usize, 2, 3, 7] {
            for block_bytes in [BLOCK_BYTES, 4096] {
                for blocks_per_thread in [2, BLOCKS_PER_THREAD] {
                    for (what, db) in [("one volume", &single), ("volumes", &multi)] {
                        let context = format!(
                            "{what}, {threads} threads, {block_bytes}-byte blocks, {blocks_per_thread} blocks per thread"
                        );
                        let made =
                            search_with(db, threads, blocks_per_thread, Driver::Made(block_bytes))
                                .unwrap();
                        let kept =
                            search_with(db, threads, blocks_per_thread, Driver::Kept(block_bytes))
                                .unwrap();
                        assert!(made.0 == reference, "{context}");
                        assert_eq!(made.1, reference_listed, "{context}");
                        assert!(kept.0 == made.0, "{context}");
                        assert_eq!(kept.1, made.1, "{context}");
                    }
                }
            }
        }
        let _ = std::fs::remove_dir_all(single_dir);
        let _ = std::fs::remove_dir_all(multi_dir);
    }

    // blast_aalookup.c:494-508 (cited at `prepare_search_block`): NCBI scans the stored
    // letters. Every driver of the preliminary stage gives the kernel each sequence as it is
    // stored, O (26) of p03 included, between its sentinel bytes, on one volume and on six, at
    // 1 and 3 threads (the window of LOSAT_X_BPDBKEEP in its two waiting sizes).
    #[test]
    fn every_driver_gives_the_kernel_the_stored_letters_o_included() {
        let single = open("prot_v4");
        let range = single.db.next_block(0, u64::MAX).unwrap();
        let stored = single.db.read_block(range).unwrap();
        assert_eq!(stored.len(), 10);
        assert!(stored.protein(2).contains(&NCBISTDAA_O));
        let volumes = open("prot_mv");
        for (what, db) in [("one volume", &single), ("volumes", &volumes)] {
            for threads in [1usize, 3] {
                let [narrow, wide] = waiting_windows(threads);
                for driver in [
                    Driver::Cut,
                    Driver::Of(BLOCK_BYTES),
                    Driver::Made(BLOCK_BYTES),
                    Driver::Kept(BLOCK_BYTES),
                    Driver::Window(BLOCK_BYTES, narrow, None),
                    Driver::Window(4096, wide, None),
                ] {
                    let context = format!("{what}, {threads} threads, {driver:?}");
                    let (seen, _) = search_with(db, threads, 2, driver).unwrap();
                    assert_eq!(seen.len(), stored.len(), "{context}");
                    for (i, (ordinal, _, letters)) in seen.iter().enumerate() {
                        assert_eq!(*ordinal, i, "{context}");
                        assert_eq!(
                            letters.as_slice(),
                            stored.protein_with_sentinels(i),
                            "{context}, sequence {i}"
                        );
                    }
                }
            }
        }
    }

    /// The windows that wait on the database of the tests: of 2 blocks per thread, and wider
    /// with bytes that narrow it to those 2 per thread as soon as any result waits.
    fn waiting_windows(threads: usize) -> [super::super::x_db_keep::XWindowSize; 2] {
        let floor = 2 * threads;
        [
            super::super::x_db_keep::XWindowSize {
                floor,
                width: floor,
                max_waiting_bytes: usize::MAX,
            },
            super::super::x_db_keep::XWindowSize {
                floor,
                width: 64 * threads,
                max_waiting_bytes: 1,
            },
        ]
    }

    // No NCBI counterpart: the window of LOSAT_X_BPDBKEEP and LOSAT_X_BPDBWIDE (`x_db_keep.rs`)
    // with a slow block, the oldest one (the later blocks fill the window and their takers wait
    // for its merge) or a later one: the kernel sees the same sequences from the same diagonal
    // offsets and the hit lists receive the same updates as with a scratch and a buffer per
    // task. A panic of the kernel in either block ends the search with the panic (the waiting
    // workers stop). At 2, 3 and 7 threads, blocks of 4 KiB, windows that wait
    // (`waiting_windows`; the default window is wider than this database).
    #[test]
    fn the_window_waits_for_its_oldest_block_and_ends_on_a_panic() {
        let (dir, name) = generated_database("keepw", "3GB");
        let database = open_path(&name);
        let block_bytes = 4096;
        assert!(x_count_blocks(&database.db, block_bytes) > 40);
        for threads in [2usize, 3, 7] {
            let made = search_with(&database, threads, 2, Driver::Made(block_bytes)).unwrap();
            for size in waiting_windows(threads) {
                for ordinal in [0usize, 150] {
                    let context = format!("{threads} threads, ordinal {ordinal}, {size:?}");
                    let slow = search_with(
                        &database,
                        threads,
                        2,
                        Driver::Window(block_bytes, size, Some(XHook::SleepAt(ordinal))),
                    )
                    .unwrap();
                    assert!(slow.0 == made.0, "{context}");
                    assert_eq!(slow.1, made.1, "{context}");
                    let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        search_with(
                            &database,
                            threads,
                            2,
                            Driver::Window(block_bytes, size, Some(XHook::PanicAt(ordinal))),
                        )
                    }));
                    assert!(panicked.is_err(), "{context}");
                }
            }
        }
        let _ = std::fs::remove_dir_all(dir);
    }

    fn error_text(result: Result<(Seen, Vec<u32>)>) -> String {
        match result {
            Ok(_) => panic!("the search of the damaged database succeeded"),
            Err(error) => format!("{error:#}"),
        }
    }

    // A damaged database fails as with one thread at every thread count: a sequence file 500
    // bytes shorter than its index says (the byte range of the read names the block), and
    // then also a byte 28 in sequence 100 (a refusal of the search block that a smaller block
    // meets before the read that fails, while the one-thread order's 1 MiB block fails to read
    // first).
    #[test]
    fn a_damaged_database_gives_the_one_thread_error_with_the_cut() {
        let (dir, name) = generated_database("damaged", "3GB");
        let psq = Path::new(&name).with_extension("psq");
        let sequence_100 = open_path(&name).db.volumes()[0]
            .index()
            .sequence_start(100)
            .unwrap() as usize;
        let mut bytes = std::fs::read(&psq).unwrap();
        bytes.truncate(bytes.len() - 500);
        std::fs::write(&psq, &bytes).unwrap();
        for foreign in [false, true] {
            if foreign {
                bytes[sequence_100 + 1_000] = BLASTAA_SIZE as u8;
                std::fs::write(&psq, &bytes).unwrap();
            }
            let database = open_path(&name);
            let expected = error_text(search(&database, 1, 2, Some(BLOCK_BYTES)));
            assert!(
                expected.contains("passes the end of the file"),
                "{expected}"
            );
            assert_eq!(
                error_text(search(&database, 4, 2, Some(BLOCK_BYTES))),
                expected
            );
            for threads in [2usize, 4, 7] {
                let cut = cut_of(&database, threads);
                // Without the one-thread order's error, the cut's own error differs.
                let own = error_text(search(&database, threads, 2, Some(cut)));
                assert_ne!(own, expected, "{threads} threads");
                // EXPERIMENT (LOSAT_X_BPDBKEEP): the kept buffers give the same errors.
                for block_bytes in [BLOCK_BYTES, cut] {
                    assert_eq!(
                        error_text(search_with(
                            &database,
                            threads,
                            2,
                            Driver::Kept(block_bytes)
                        )),
                        error_text(search_with(
                            &database,
                            threads,
                            2,
                            Driver::Made(block_bytes)
                        )),
                        "{threads} threads, {block_bytes}-byte blocks"
                    );
                    // Windows whose workers wait when the failed block stops the search.
                    for size in waiting_windows(threads) {
                        assert_eq!(
                            error_text(search_with(
                                &database,
                                threads,
                                2,
                                Driver::Window(block_bytes, size, None)
                            )),
                            error_text(search_with(
                                &database,
                                threads,
                                2,
                                Driver::Made(block_bytes)
                            )),
                            "{threads} threads, {block_bytes}-byte blocks, {size:?}"
                        );
                    }
                }
                if foreign {
                    assert!(own.contains("is no NCBIstdaa residue"), "{own}");
                }
                assert_eq!(
                    error_text(search(&database, threads, 2, None)),
                    expected,
                    "{threads} threads"
                );
            }
        }
        let _ = std::fs::remove_dir_all(dir);
    }
}

// Line RG item X4-excl: under the negative seqid list of the rows `db.x4excl.*`
// (`tests/db_search/x4excl/fwd`), sx_frag (ordinal 3, a ti-only line) is kept without a line,
// so the redo's fetch excludes it (`exclude_fetched`): its lists leave the hit lists and leave
// marks of their places, by query, each with the number of kept lists before it in its
// query's order. Without the list, or without the fetches (`-ungapped`), nothing leaves and
// there is no mark.
#[test]
fn excluded_lists_leave_marks_of_their_places() {
    let name = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/db_search/x4excl/fwd");
    let mut warnings = Vec::new();
    let negative = crate::seqdb::DbList::Seqids(crate::seqdb::SeqidList {
        negative: true,
        file: crate::seqdb::seqidlist::read_si_list(b"XX999999.1\n", &mut warnings).unwrap(),
    });
    let settings = ApplicationSettings::default();
    // The hit list takes a list's best e-value from its HSPs (`update_best_evalue`).
    let list = |oid: u32, raw_score: i32, best_evalue: f64| BlastpHspList {
        oid,
        query_index: 0,
        hsps: vec![BlastpHsp {
            e_value: best_evalue,
            ..test_hsp(oid, raw_score)
        }],
        best_evalue,
    };
    // Query 0: sexact, sfail, sx_frag; query 2: sx_frag, sok1; query 1 has no list.
    let make = || {
        let mut first = BlastpHitList::new(10);
        for (oid, score, evalue) in [(0, 430, 1e-49), (1, 419, 1e-47), (3, 26, 0.04)] {
            first.update(list(oid, score, evalue));
        }
        let mut third = BlastpHitList::new(10);
        for (oid, score, evalue) in [(3, 90, 1e-9), (2, 60, 1e-5)] {
            third.update(list(oid, score, evalue));
        }
        vec![Some(first), None, Some(third)]
    };
    let oids = |lists: &[Option<BlastpHitList>]| -> Vec<Vec<u32>> {
        lists
            .iter()
            .flatten()
            .map(|hit_list| hit_list.hsplist_array.iter().map(|l| l.oid).collect())
            .collect()
    };
    let database = BlastpDatabase::open_restricted(
        name.as_bytes(),
        &settings,
        "BLASTP",
        Some(&negative),
        None,
        &mut std::io::sink(),
    )
    .unwrap();
    let mut lists = make();
    let listed = ListedSubjects::x_list(&database, &mut lists, &mut [], true).unwrap();
    assert_eq!(listed.ordinals, [0, 1, 2]);
    assert_eq!(oids(&lists), [vec![0, 1], vec![2]]);
    let mark = |query, position, ordinal, best_evalue| FetchExcludedMatch {
        query,
        position,
        ordinal,
        best_evalue,
    };
    assert_eq!(
        listed.fetch_excluded,
        [mark(0, 2, 3, 0.04), mark(2, 0, 3, 1e-9)]
    );
    // An ungapped search fetches no subject (`fetched` false): sx_frag stays in both lists, with
    // no mark (it leaves at the Seq-aligns, `drop_unaligned`).
    let mut lists = make();
    let listed = ListedSubjects::x_list(&database, &mut lists, &mut [], false).unwrap();
    assert_eq!(listed.ordinals, [0, 1, 2, 3]);
    assert_eq!(oids(&lists), [vec![0, 1, 3], vec![3, 2]]);
    assert!(listed.fetch_excluded.is_empty());
    // There it leaves at the Seq-aligns (ranks equal ordinals here); a gapped search's
    // `drop_unaligned` leaves the lists as they are (its fetches have excluded it).
    let mut kept = lists.clone();
    listed.drop_unaligned(&database, &mut kept, true);
    assert_eq!(oids(&kept), [vec![0, 1, 3], vec![3, 2]]);
    listed.drop_unaligned(&database, &mut lists, false);
    assert_eq!(oids(&lists), [vec![0, 1], vec![2]]);
    let database = BlastpDatabase::open(name, &settings, "BLASTP").unwrap();
    let mut lists = make();
    let listed = ListedSubjects::x_list(&database, &mut lists, &mut [], true).unwrap();
    assert_eq!(listed.ordinals, [0, 1, 2, 3]);
    assert_eq!(oids(&lists), [vec![0, 1, 3], vec![3, 2]]);
    assert!(listed.fetch_excluded.is_empty());
}

// `x_list` lists the subjects and renames them to their ranks with the stored lengths and no
// copy; the redo reads the stored residues by rank; after the redo `x_read_shown` keeps the
// residues that the reports read (p03 and p08 hold an O).
#[test]
fn listed_subjects_read_the_stored_residues_by_rank() {
    let database = open("prot_mv");
    let list = |oid: u32, raw_score: i32| BlastpHspList {
        oid,
        query_index: 0,
        hsps: vec![test_hsp(oid, raw_score), test_hsp(oid, raw_score - 1)],
        best_evalue: 1e-5,
    };
    let make = || {
        let mut first = BlastpHitList::new(10);
        for (oid, score) in [(7, 60), (2, 55), (9, 50)] {
            first.update(list(oid, score));
        }
        let mut second = BlastpHitList::new(10);
        second.update(list(2, 70));
        second.update(list(0, 65));
        vec![Some(first), None, Some(second)]
    };
    let mut lists = make();
    let listed = ListedSubjects::x_list(&database, &mut lists, &mut [], true).unwrap();
    assert_eq!(listed.ordinals, [0, 2, 7, 9]);
    assert!(listed.encoded.is_empty());
    // The subjects are renamed to their ranks: ordinals 7, 2, 9 and 2, 0 are ranks 2, 1, 3
    // and 1, 0.
    let ranks = |lists: &[Option<BlastpHitList>]| -> Vec<Vec<u32>> {
        lists
            .iter()
            .flatten()
            .map(|hit_list| hit_list.hsplist_array.iter().map(|l| l.oid).collect())
            .collect()
    };
    assert_eq!(ranks(&lists), [vec![2, 1, 3], vec![1, 0]]);
    // The stored residues of each rank (between the sentinels), straight from the database.
    let stored: Vec<Vec<u8>> = listed
        .ordinals
        .iter()
        .map(|&ordinal| {
            database
                .db
                .read_sequence(ordinal)
                .unwrap()
                .protein(0)
                .to_vec()
        })
        .collect();
    let lengths: Vec<usize> = stored.iter().map(Vec::len).collect();
    assert_eq!(lengths[3], 300);
    let subjects = listed.redo_subjects(&database, None);
    assert_eq!(subjects.len(), 4);
    for (rank, &length) in lengths.iter().enumerate() {
        assert_eq!(subjects.aa_len(rank), Some(length));
        let read = subjects.read(rank as u32).unwrap();
        assert_eq!(read.residues(), &stored[rank][..]);
    }
    assert_eq!(subjects.aa_len(4), None);
    // A final hit list of ranks 1 (p03, with an O) and 3 (without); rank 2 (p08) holds an O but
    // is in no final list.
    let mut final_list = BlastpHitList::new(10);
    final_list.update(list(1, 60));
    final_list.update(list(3, 50));
    let final_lists = vec![None, Some(final_list)];
    // With `WithO`, the HSPs of a query with an O (query 0 here) keep every subject's residues.
    for (residues, query_recounted, kept) in [
        (XShownResidues::WithO, &[][..], [false, true, false, false]),
        (
            XShownResidues::WithO,
            &[false, true],
            [false, true, false, false],
        ),
        (XShownResidues::WithO, &[true], [false, true, false, true]),
        (XShownResidues::Final, &[][..], [false, true, false, true]),
        (XShownResidues::Listed, &[][..], [true; 4]),
    ] {
        let mut shown = ListedSubjects::x_list(&database, &mut make(), &mut [], true).unwrap();
        // The redo reads the subjects of the final lists (and finds the O of p03).
        for rank in [1, 3] {
            shown.redo_subjects(&database, None).read(rank).unwrap();
        }
        shown
            .x_read_shown(&database, &final_lists, residues, query_recounted)
            .unwrap();
        assert_eq!(shown.encoded.len(), 4, "{residues:?}");
        for (rank, kept) in kept.into_iter().enumerate() {
            assert_eq!(
                shown.encoded[rank].aa_len, lengths[rank],
                "{residues:?} {rank}"
            );
            if kept {
                // Between two sentinel bytes.
                assert_eq!(
                    shown.encoded[rank].aa_seq.len(),
                    lengths[rank] + 2,
                    "{residues:?} {rank}"
                );
                assert_eq!(
                    &shown.encoded[rank].aa_seq[1..=lengths[rank]],
                    &stored[rank][..],
                    "{residues:?} {rank}"
                );
            } else {
                assert!(shown.encoded[rank].aa_seq.is_empty(), "{residues:?} {rank}");
            }
        }
    }
    // `read` (the unit tests' listed subjects) is `Listed`.
    let read = ListedSubjects::read(&database, &mut make(), &mut []).unwrap();
    assert_eq!(read.encoded[2].aa_len, lengths[2]);
    assert_eq!(read.encoded[3].aa_seq.len(), 302);
}

// The listed subjects stand for the subject set after the preliminary stage: ranks
// in ordinal order, with the stored letters (O is O; p03 and p08 of the fixture).
#[test]
fn listed_subjects_are_read_by_ordinal_and_renamed_to_ranks() {
    let database = open("prot_mv");
    let list = |oid: u32, raw_score: i32| BlastpHspList {
        oid,
        query_index: 0,
        hsps: vec![test_hsp(oid, raw_score), test_hsp(oid, raw_score - 1)],
        best_evalue: 1e-5,
    };
    let mut first = BlastpHitList::new(10);
    for (oid, score) in [(7, 60), (2, 55), (9, 50)] {
        first.update(list(oid, score));
    }
    let mut second = BlastpHitList::new(10);
    second.update(list(2, 70));
    second.update(list(0, 65));
    let mut hit_lists = vec![Some(first), None, Some(second)];
    let listed = ListedSubjects::read(&database, &mut hit_lists, &mut []).unwrap();
    assert_eq!(listed.ordinals, [0, 2, 7, 9]);
    let titles = listed.read_titles(&database, &hit_lists, &[0]).unwrap();
    let oids = |hit_list: &Option<BlastpHitList>| -> Vec<(u32, Vec<u32>)> {
        let hit_list = hit_list.as_ref().unwrap();
        hit_list.hsplist_array[..hit_list.hsplist_count]
            .iter()
            .map(|l| (l.oid, l.hsps.iter().map(|hsp| hsp.s_idx).collect()))
            .collect()
    };
    let mut renamed = oids(&hit_lists[0]);
    renamed.sort();
    assert_eq!(renamed, [(1, vec![1, 1]), (2, vec![2, 2]), (3, vec![3, 3])]);
    let mut renamed = oids(&hit_lists[2]);
    renamed.sort();
    assert_eq!(renamed, [(0, vec![0, 0]), (1, vec![1, 1])]);
    // p03: every letter; the redo and the report read the stored letters (O is O).
    assert_eq!(
        titles.records()[1].title(),
        b"p03 every letter of NCBIstdaa but the gap"
    );
    assert_eq!(
        titles.records()[1].seq(),
        b"ABCDEFGHIKLMNPQRSTVWXYZUOJMHCENQT"
    );
    let redone: String = listed.encoded[1].aa_seq[1..=listed.encoded[1].aa_len]
        .iter()
        .map(|&residue| ncbistdaa_to_ascii(residue))
        .collect();
    assert_eq!(redone, "ABCDEFGHIKLMNPQRSTVWXYZUOJMHCENQT");
    assert_eq!(listed.encoded[1].aa_seq.len(), 35);
    assert_eq!(
        titles.records()[2].seq(),
        b"MUGYDNEDMDNEPCNDCRACIGFRYKYMMVKWOFASVIHWKHVWDJ"
    );
    assert_eq!(listed.encoded[3].aa_len, 300);
    // The tabular ids of the subjects with rows, the stored titles for `stitle`.
    let (ids, accessions) = titles.tabular_subject_ids(&hit_lists, true).unwrap();
    assert!(accessions.is_empty());
    assert_eq!(
        ids.iter().map(|id| id.as_deref()).collect::<Vec<_>>(),
        [Some(&b"p01"[..]), Some(b"p03"), Some(b"p08"), Some(b"p10")]
    );
    assert_eq!(
        titles.stitles()[0].as_ref(),
        b"p01 first protein [Test organism]"
    );
    titles
        .check_shown_titles(
            &database,
            &[0usize, 1, 2, 3].into_iter().collect(),
            &HitLoads::default(),
        )
        .unwrap();
    titles
        .check_reported_titles(&database, &hit_lists, shared::TitleReads::Reports)
        .unwrap();
}

// EXPERIMENT (LOSAT_X_BPLISTEDLEAN): the lean titles are the records' titles, rank by rank,
// with one shared empty title for the subjects in no final hit list, and the stored residues
// give the records' letters (`x_listed_lean.rs`).
#[test]
fn lean_titles_are_the_records_titles_and_letters() {
    let database = open("prot_mv");
    let list = |oid: u32, raw_score: i32| BlastpHspList {
        oid,
        query_index: 0,
        hsps: vec![test_hsp(oid, raw_score)],
        best_evalue: 1e-5,
    };
    let mut first = BlastpHitList::new(10);
    for (oid, score) in [(7, 60), (2, 55), (9, 50)] {
        first.update(list(oid, score));
    }
    let mut second = BlastpHitList::new(10);
    second.update(list(2, 70));
    second.update(list(0, 65));
    let mut hit_lists = vec![Some(first), None, Some(second)];
    let listed = ListedSubjects::read(&database, &mut hit_lists, &mut []).unwrap();
    // Ordinal ids 7 and 9 (ranks 2 and 3) stay listed but leave the final hit lists.
    let first = hit_lists[0].as_mut().unwrap();
    first.hsplist_array.retain(|hsp_list| hsp_list.oid == 1);
    first.hsplist_count = first.hsplist_array.len();
    let records = listed.read_titles(&database, &hit_lists, &[0]).unwrap();
    let lean = listed
        .x_read_titles_lean(&database, &hit_lists, &[0])
        .unwrap();
    let ListedShown::Lean(lean_titles) = &lean.shown else {
        panic!("lean titles");
    };
    assert!(matches!(records.shown, ListedShown::Records(_)));
    assert_eq!(lean_titles.len(), 4);
    for rank in 0..4 {
        assert_eq!(
            lean.title(rank),
            records.records()[rank].title(),
            "rank {rank}"
        );
        assert_eq!(
            lean.replaced_title_bytes[rank],
            records.replaced_title_bytes[rank]
        );
    }
    assert_eq!(lean.title(1), b"p03 every letter of NCBIstdaa but the gap");
    assert!(lean.title(2).is_empty() && lean.title(3).is_empty());
    assert!(Arc::ptr_eq(&lean_titles[2], &lean_titles[3]));
    assert_eq!(
        lean.tabular_subject_ids(&hit_lists, true).unwrap(),
        records.tabular_subject_ids(&hit_lists, true).unwrap()
    );
    lean.check_shown_titles(
        &database,
        &[0usize, 1].into_iter().collect(),
        &HitLoads::default(),
    )
    .unwrap();
    lean.check_reported_titles(&database, &hit_lists, shared::TitleReads::Reports)
        .unwrap();
    assert_eq!(lean.stitles(), records.stitles());
    // Titles, flags and every letter against the stored residues.
    x_assert_same_titles_and_letters(&lean, &records);
}

// blast_seqalign.cpp:1569-1582 and seqdbvol.cpp:2302-2310: NCBI opens a volume's
// header file when the Seq-aligns of the final hit lists are made, query by query and
// in the order of each hit list; a subject that is in no final hit list is never read.
// prot_mv: volume 00 holds the ordinal ids 0-2, 01 holds 3, 05 holds 9.
// (A path in the system's temporary directory: Unix spelling.)
#[cfg(unix)]
#[test]
fn titles_are_read_for_the_final_hit_lists_in_ncbis_order() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-db-nohdr-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    for entry in std::fs::read_dir(fixture("prot_mv.pal").trim_end_matches("prot_mv.pal")).unwrap()
    {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_str().unwrap().to_string();
        if name.starts_with("prot_mv.") && name != "prot_mv.00.phr" && name != "prot_mv.05.phr" {
            std::fs::copy(&path, dir.join(&name)).unwrap();
        }
    }
    let name = dir.join("prot_mv").to_str().unwrap().to_string();
    let database = BlastpDatabase::open(&name, &ApplicationSettings::default(), "BLASTP").unwrap();
    let list = |oid: u32, best_evalue: f64| BlastpHspList {
        oid,
        query_index: 0,
        hsps: vec![test_hsp(oid, 60)],
        best_evalue,
    };
    let hit_list = |lists: Vec<BlastpHspList>, count: usize| {
        let mut hit_list = BlastpHitList::new(10);
        hit_list.hsplist_array = lists;
        hit_list.hsplist_count = count;
        Some(hit_list)
    };
    // Query 1 reports ordinal id 3 (volume 01); ordinal id 1 (volume 00, no header
    // file) is in its list past the count, and is not read.
    let mut hit_lists = vec![None, hit_list(vec![list(3, 1e-9), list(1, 1e-3)], 1)];
    let listed = ListedSubjects::read(&database, &mut hit_lists, &mut []).unwrap();
    assert_eq!(listed.ordinals, [1, 3]);
    let titles = listed.read_titles(&database, &hit_lists, &[0]).unwrap();
    assert_eq!(
        titles.records()[1].title(),
        b"p04 lower case and ambiguity codes"
    );
    assert_eq!(titles.records()[0].title(), b"");
    assert_eq!(titles.records()[0].seq(), b"I");
    // EXPERIMENT (LOSAT_X_BPLISTEDLEAN): the lean titles do not read it either.
    let lean = listed
        .x_read_titles_lean(&database, &hit_lists, &[0])
        .unwrap();
    assert_eq!(lean.title(1), b"p04 lower case and ambiguity codes");
    assert_eq!(lean.title(0), b"");
    // Query 2 reports ordinal id 9 (volume 05) before ordinal id 0 (volume 00): the
    // first header file read is that of volume 05, in query 2's batch.
    let mut hit_lists = vec![
        None,
        hit_list(vec![list(3, 1e-9), list(1, 1e-3)], 1),
        hit_list(vec![list(9, 1e-9), list(0, 1e-5)], 2),
    ];
    let listed = ListedSubjects::read(&database, &mut hit_lists, &mut []).unwrap();
    assert_eq!(listed.ordinals, [0, 1, 3, 9]);
    let error = listed
        .read_titles(&database, &hit_lists, &[0])
        .err()
        .expect("a missing header file")
        .downcast::<ErrorAfterBatchReports>()
        .expect("NCBI's error, with the query that reached it");
    assert_eq!(error.query, Some(2));
    let native = error
        .error
        .downcast_ref::<NativeError>()
        .expect("NCBI's error");
    assert_eq!(native.exit, 2);
    // EXPERIMENT (LOSAT_X_BPLISTEDLEAN): the lean titles are read in the same order, with
    // the same error.
    let lean_error = listed
        .x_read_titles_lean(&database, &hit_lists, &[0])
        .err()
        .expect("a missing header file")
        .downcast::<ErrorAfterBatchReports>()
        .expect("NCBI's error, with the query that reached it");
    assert_eq!(lean_error.query, Some(2));
    assert_eq!(
        lean_error
            .error
            .downcast_ref::<NativeError>()
            .unwrap()
            .message,
        native.message
    );
    // The volume is recorded under its directory normalised with symbolic links
    // followed (seqdbalias.cpp:851-862, `NormalizePath` of `dir_name` with `eFollowLinks`), so
    // the message names that directory: on macOS the temporary directory under `/var`
    // is shown under `/private/var`.
    assert_eq!(
        String::from_utf8(native.message.clone()).unwrap(),
        format!(
            "LOSAT database error: Error: File ({}) not found.\n",
            dir.canonicalize().unwrap().join("prot_mv.05.phr").display()
        )
    );
    drop(database);
    std::fs::remove_dir_all(&dir).unwrap();
}

// NCBI BLAST+ 2.17.0 writes `Error: (803.7) [blastp] Blast-def-line-set.E.title`,
// `Bad char [0xE9] in string at byte 24` and the stored title to standard error
// each time it reads such a definition line past SeqDB's header cache
// (`check_reported_titles`, `title_posts`).
// (A path in the system's temporary directory or in an error text: Unix spelling.)
#[cfg(unix)]
#[test]
fn a_reported_title_with_a_byte_outside_the_visible_range_is_rejected() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-db-badtitle-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    for ext in ["pin", "psq"] {
        std::fs::copy(fixture("prot_v4") + "." + ext, dir.join(format!("b.{ext}"))).unwrap();
    }
    // The title of p10 (`p10 last protein`) with a byte of Latin-1 in it.
    let mut header = std::fs::read(fixture("prot_v4") + ".phr").unwrap();
    let at = header
        .windows(16)
        .position(|window| window == b"p10 last protein")
        .expect("the title of p10");
    header[at + 5] = 0xE9;
    std::fs::write(dir.join("b.phr"), header).unwrap();
    let name = dir.join("b").to_str().unwrap().to_string();
    let database = BlastpDatabase::open(&name, &ApplicationSettings::default(), "BLASTP").unwrap();
    let list = |oid: u32| BlastpHspList {
        oid,
        query_index: 0,
        hsps: vec![test_hsp(oid, 60)],
        best_evalue: 1e-5,
    };
    // A subject whose title is clean: nothing to reject.
    let mut clean = BlastpHitList::new(10);
    clean.update(list(0));
    let mut hit_lists = vec![Some(clean)];
    let listed = ListedSubjects::read(&database, &mut hit_lists, &mut []).unwrap();
    let titles = listed.read_titles(&database, &hit_lists, &[0]).unwrap();
    titles
        .check_reported_titles(&database, &hit_lists, shared::TitleReads::Reports)
        .unwrap();
    // The subject with the bad byte: read as `#`, and rejected when it is reported.
    let mut bad = BlastpHitList::new(10);
    bad.update(list(9));
    bad.update(list(0));
    let mut hit_lists = vec![Some(bad)];
    let listed = ListedSubjects::read(&database, &mut hit_lists, &mut []).unwrap();
    let titles = listed.read_titles(&database, &hit_lists, &[0]).unwrap();
    assert_eq!(titles.records()[1].title(), b"p10 l#st protein");
    let error = titles
        .check_reported_titles(&database, &hit_lists, shared::TitleReads::Reports)
        .unwrap_err();
    assert!(
        format!("{error:#}")
            .starts_with("database sequence 9 has a stored title with a byte outside ' ' to '~'"),
        "{error:#}"
    );
    // EXPERIMENT (LOSAT_X_BPLISTEDLEAN): the lean titles hold the same title and flag, and
    // give the same refusal.
    let lean = listed
        .x_read_titles_lean(&database, &hit_lists, &[0])
        .unwrap();
    assert_eq!(lean.title(1), b"p10 l#st protein");
    let lean_error = lean
        .check_reported_titles(&database, &hit_lists, shared::TitleReads::Reports)
        .unwrap_err();
    assert_eq!(format!("{lean_error:#}"), format!("{error:#}"));
    // Where only the Seq-aligns read the definition lines (outfmt 8, 9, 11 and 12), NCBI's
    // posts are written instead: one per bad byte at the first read of the subject, none at a
    // read that SeqDB's header cache answers.
    titles
        .check_reported_titles(&database, &hit_lists, shared::TitleReads::SeqAligns)
        .unwrap();
    let posts = titles
        .title_posts(&database, &hit_lists, [&(0..1)])
        .unwrap();
    assert_eq!(posts.len(), 1);
    assert_eq!(posts[0].0, 0);
    assert_eq!(
        posts[0].1.text(),
        b"Error: (803.7) [blastp] Blast-def-line-set.E.title\nBad char [0xE9] in string at \
          byte 24\np10 l\xe9st protein\n"
    );
    assert!(titles
        .title_posts(&database, &hit_lists, [&(0..1)])
        .unwrap()
        .is_empty());
    drop(database);
    std::fs::remove_dir_all(&dir).unwrap();
}

// The compact preliminary lists (`x_prelim_compact.rs`) are listed and renamed to ranks with
// the lists of full HSPs, before they are expanded.
#[test]
fn compact_lists_are_listed_and_renamed_with_the_others() {
    use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
    let database = open("prot_mv");
    let compact = |oid: u32| BlastpHspList {
        oid,
        query_index: 1,
        hsps: vec![XPrelimHsp {
            hsp: BlastpPreliminaryHsp {
                query_start: 0,
                query_end: 10,
                subject_start: 0,
                subject_end: 10,
                gapped_query_start: 1,
                gapped_subject_start: 1,
                raw_score: 60,
                query_context: 1,
                query_frame: 0,
                subject_frame: 0,
                query_length: 20,
                q_idx: 1,
                s_idx: oid,
            },
            e_value: 1e-5,
        }],
        best_evalue: 1e-5,
    };
    let mut off = BlastpHitList::new(10);
    off.update(BlastpHspList {
        oid: 7,
        query_index: 0,
        hsps: vec![test_hsp(7, 60)],
        best_evalue: 1e-5,
    });
    let mut compact_list = BlastpHitList::new(10);
    compact_list.update(compact(9));
    compact_list.update(compact(2));
    let mut hit_lists = vec![Some(off), None];
    let mut compact_lists = vec![None, Some(compact_list)];
    let listed = ListedSubjects::read(&database, &mut hit_lists, &mut compact_lists).unwrap();
    assert_eq!(listed.ordinals, [2, 7, 9]);
    assert_eq!(listed.encoded.len(), 3);
    let off = &hit_lists[0].as_ref().unwrap().hsplist_array[0];
    assert_eq!((off.oid, off.hsps[0].s_idx), (1, 1));
    let mut renamed: Vec<(u32, u32)> = compact_lists[1]
        .as_ref()
        .unwrap()
        .hsplist_array
        .iter()
        .map(|list| (list.oid, list.hsps[0].hsp.s_idx))
        .collect();
    renamed.sort();
    assert_eq!(renamed, [(0, 0), (2, 2)]);
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): the one pass over the headers (`SeqDb::visit_deflines`)
// gives the answer of the two checks it replaces (`x_titles_read_allow_blocks`, then
// `x_taxids_allow_blocks` one sequence at a time), for a database whose headers read, one with
// taxids, an alias of volumes, and one whose header file is missing, at one thread and in a pool.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn titles_and_taxids_in_one_pass_give_the_answers_of_the_two_checks() {
    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let paja = manifest_dir.join("tests/db_search/paja_taxid");
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-db-qblock-onepass-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    for ext in ["pin", "psq"] {
        std::fs::copy(paja.with_extension(ext), dir.join(format!("nohdr.{ext}"))).unwrap();
    }
    let mut databases = vec![
        (open("prot_v4"), true),
        (open("prot_mv"), true),
        (
            BlastpDatabase::open(
                paja.to_str().unwrap(),
                &ApplicationSettings::default(),
                "BLASTP",
            )
            .unwrap(),
            true,
        ),
    ];
    // A database whose header file is missing: the search of all queries reports NCBI's error.
    if let Ok(missing) = BlastpDatabase::open(
        dir.join("nohdr").to_str().unwrap(),
        &ApplicationSettings::default(),
        "BLASTP",
    ) {
        databases.push((missing, false));
    }
    for (database, readable) in &databases {
        for threads in [1, 4, 8, 16] {
            crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                let two = x_titles_read_allow_blocks(database, pool)
                    && x_taxids_allow_blocks(database, pool);
                assert_eq!(two, *readable);
                assert_eq!(
                    x_titles_taxids_allow_blocks(database, pool),
                    two,
                    "threads {threads}"
                );
                Ok(())
            })
            .unwrap();
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): the block search of `-db` runs only when every definition
// line of the database reads and none is refused (`x_titles_allow_blocks`), at one thread and
// in a pool. (A path in the system's temporary directory: Unix spelling.)
#[cfg(unix)]
#[test]
fn blocks_are_allowed_only_when_every_title_reads_and_is_reported() {
    let allowed = |database: &BlastpDatabase| -> [bool; 2] {
        [1, 4].map(|threads| {
            crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                Ok(x_titles_allow_blocks(database, false, false, pool))
            })
            .unwrap()
        })
    };
    // EXPERIMENT (LOSAT_X_BPQBLOCK): the same databases for the
    // read-only check of outfmt 8, 9, 12 and 5.
    let read_allowed = |database: &BlastpDatabase| -> [bool; 2] {
        [1, 4].map(|threads| {
            crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                Ok(x_titles_read_allow_blocks(database, pool))
            })
            .unwrap()
        })
    };
    assert_eq!(allowed(&open("prot_v4")), [true, true]);
    assert_eq!(read_allowed(&open("prot_v4")), [true, true]);
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-db-qblock-titles-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    // A header file that cannot be read: the search of all queries reports NCBI's error.
    for entry in std::fs::read_dir(fixture("prot_mv.pal").trim_end_matches("prot_mv.pal")).unwrap()
    {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_str().unwrap().to_string();
        if name.starts_with("prot_mv.") && name != "prot_mv.05.phr" {
            std::fs::copy(&path, dir.join(&name)).unwrap();
        }
    }
    let missing = BlastpDatabase::open(
        dir.join("prot_mv").to_str().unwrap(),
        &Default::default(),
        "BLASTP",
    )
    .unwrap();
    assert_eq!(allowed(&missing), [false, false]);
    assert_eq!(read_allowed(&missing), [false, false]);
    // A stored title with a byte outside ' ' ..= '~' (LOSAT's refusal when it is reported).
    for ext in ["pin", "psq"] {
        std::fs::copy(fixture("prot_v4") + "." + ext, dir.join(format!("b.{ext}"))).unwrap();
    }
    let mut header = std::fs::read(fixture("prot_v4") + ".phr").unwrap();
    let at = header
        .windows(16)
        .position(|window| window == b"p10 last protein")
        .expect("the title of p10");
    header[at + 5] = 0xE9;
    std::fs::write(dir.join("b.phr"), header).unwrap();
    let bad = BlastpDatabase::open(
        dir.join("b").to_str().unwrap(),
        &Default::default(),
        "BLASTP",
    )
    .unwrap();
    assert_eq!(allowed(&bad), [false, false]);
    assert_eq!(read_allowed(&bad), [false, false]);
    // A title of punctuation only (`, ,`): the tabular id NCBI writes cannot be made from it
    // (outfmt 6 and 7 refuse it, `x_titles_allow_blocks`), but outfmt 8, 9, 12 and 5 write no
    // id, so `x_titles_read_allow_blocks` takes the blocks (the stricter function must not
    // replace it).
    for ext in ["pin", "psq"] {
        std::fs::copy(fixture("prot_v4") + "." + ext, dir.join(format!("c.{ext}"))).unwrap();
    }
    let mut header = std::fs::read(fixture("prot_v4") + ".phr").unwrap();
    let at = header
        .windows(16)
        .position(|window| window == b"p10 last protein")
        .expect("the title of p10");
    header[at..at + 16].copy_from_slice(b", , , , , , , , ");
    std::fs::write(dir.join("c.phr"), header).unwrap();
    let punctuation = BlastpDatabase::open(
        dir.join("c").to_str().unwrap(),
        &Default::default(),
        "BLASTP",
    )
    .unwrap();
    assert_eq!(allowed(&punctuation), [false, false]);
    assert_eq!(read_allowed(&punctuation), [true, true]);
    // An empty title (line DBID H3), which outfmt 0 shows as NCBI makes it from the name of
    // its taxid (`StoredTitle::source`; line RG T5): the blocks take it. Its tabular id is
    // `unnamed`.
    let bare = BlastpDatabase::open(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/db_search/ids/pbare_t")
            .to_str()
            .unwrap(),
        &Default::default(),
        "BLASTP",
    )
    .unwrap();
    let pairwise_allowed = [1, 4].map(|threads| {
        crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
            Ok(x_titles_allow_blocks(&bare, true, false, pool))
        })
        .unwrap()
    });
    assert_eq!(pairwise_allowed, [true, true]);
    assert_eq!(allowed(&bare), [true, true]);
    // A database with parsed Seq-ids (line DBID; RefSeq records among them, whose title made
    // from the source is the name alone and can read past its end): its definition lines are
    // not `SeqDb::deflines`, so the blocks never take it and the search of all queries checks
    // each title from its own source (`check_shown_parsed_subject`).
    let parsed = BlastpDatabase::open(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/db_search/ids/prot_v4")
            .to_str()
            .unwrap(),
        &Default::default(),
        "BLASTP",
    )
    .unwrap();
    for pairwise in [false, true] {
        let parsed_allowed = [1, 4].map(|threads| {
            crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                Ok(x_titles_allow_blocks(&parsed, pairwise, false, pool))
            })
            .unwrap()
        });
        assert_eq!(parsed_allowed, [false, false], "pairwise {pairwise}");
    }
    drop((missing, bad, punctuation, bare, parsed));
    std::fs::remove_dir_all(&dir).unwrap();
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): where the batches' prefetch reads the taxonomy database
// (`taxids`, DBID item T), the block search of `-db` runs only when no lookup can be refused:
// the one pass over the headers gives the answer of the titles' check followed by
// `x_taxids_allow_blocks`. `sprot_tax` of `tests/fixtures/outfmt16/db` (two sequences, taxid
// 9606) with no taxonomy database, ones whose lookups cannot be refused (a wrong magic number
// is a message; a scientific name of 9606 with bytes >= 0x80 is not asked for by the
// prefetch) and one whose record of 9606 is (an empty `taxdb.btd`, `shared::taxdb_may_refuse`),
// at one thread and in a pool.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn blocks_with_the_prefetch_need_lookups_that_cannot_be_refused() {
    use crate::seqdb::alias::SearchPath;
    let tests = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let cwd = crate::utils::test_dir::temp_dir().join(format!(
        "losat-db-qblock-prefetch-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&cwd).unwrap();
    let path = tests.join("fixtures/outfmt16/db/sprot_tax");
    for (taxdb, refused) in [
        (None, false),
        (Some("good"), false),
        (Some("bad_magic"), false),
        (Some("high_bytes_utf8"), false),
        (Some("btd_empty"), true),
    ] {
        let taxdb = taxdb.map(|name| tests.join("seqdb/taxdb").join(name));
        let search = SearchPath::new(&cwd, taxdb.as_deref().map(std::path::Path::as_os_str), None);
        let database = BlastpDatabase {
            db: SeqDb::open_in(path.to_str().unwrap().as_bytes(), SeqType::Protein, &search)
                .unwrap(),
            name: b"sprot_tax".to_vec(),
            taxdb: shared::TaxdbRun::default(),
            program: "BLASTP",
            rps_units: None,
        };
        assert_eq!(database.taxdb_found(), taxdb.is_some());
        for threads in [1, 4, 8, 16] {
            crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                for pairwise in [false, true] {
                    let titles = x_titles_allow_blocks(&database, pairwise, false, pool);
                    assert!(titles, "{taxdb:?} threads {threads}");
                    let two = titles && x_taxids_allow_blocks(&database, pool);
                    assert_eq!(two, !refused, "{taxdb:?} threads {threads}");
                    assert_eq!(
                        x_titles_allow_blocks(&database, pairwise, true, pool),
                        two,
                        "{taxdb:?} threads {threads}"
                    );
                }
                Ok(())
            })
            .unwrap();
        }
    }
    std::fs::remove_dir_all(&cwd).unwrap();
}
