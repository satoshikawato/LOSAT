//! NCBI's unit tests of the mask-data writer (`writedb_unit_test.cpp`, cases
//! `RegisterMaskingAlgorithms` to `MaskDataBoundsError`), ported with their assertions. NCBI's
//! cases read sequences from its test database `data/writedb_prot`; the ports write synthetic
//! protein sequences of the lengths the cases need (`MaskDataColumn` takes sequences of at
//! least 1024 residues, the other two cases use any sequence and its length) and, where NCBI
//! leaves its files untested ("the SeqDB test uses copies of these files"), read the database
//! back with LOSAT's reader.

use super::registry::{
    FILTER_PROGRAM_DUST, FILTER_PROGRAM_MAX, FILTER_PROGRAM_OTHER, FILTER_PROGRAM_REPEAT,
    FILTER_PROGRAM_SEG, FILTER_PROGRAM_WINDOWMASKER,
};
use super::tests::{bioseq, new_db};
use super::{MaskInfoRegistry, MaskedRanges};
use crate::seqdb::{MaskRanges, SeqDb, SeqType};
use crate::writedb::tests::{files_named, TempDir};
use crate::writedb::writer::DbVersion;
use crate::writedb::WriteDbError;

/// A `CWriteDBException`: `WriteDbError::Creation` in LOSAT.
fn is_write_db_exception(error: &WriteDbError) -> bool {
    matches!(error, WriteDbError::Creation(_))
}

// NCBI unit test (598d8ae6): c++/src/objtools/blast/seqdb_writer/unit_test/writedb_unit_test.cpp:1807-1824 RegisterMaskingAlgorithms
#[test]
fn register_masking_algorithms() {
    let mut registry = MaskInfoRegistry::default();
    let ids = [
        registry.add_program(FILTER_PROGRAM_SEG, b"", b"").unwrap(),
        registry.add_program(FILTER_PROGRAM_DUST, b"", b"").unwrap(),
        registry
            .add_program(FILTER_PROGRAM_WINDOWMASKER, b"", b"")
            .unwrap(),
        registry
            .add_program(FILTER_PROGRAM_REPEAT, b"9606", b"")
            .unwrap(),
        registry
            .add_program(FILTER_PROGRAM_OTHER, b"dummy1", b"")
            .unwrap(),
        registry
            .add_program(FILTER_PROGRAM_OTHER, b"dummy2", b"")
            .unwrap(),
        registry
            .add_program(FILTER_PROGRAM_OTHER, b"dummy2", b"funnyname")
            .unwrap(),
    ];
    for id in ids {
        assert!(registry.is_registered(id), "{id}");
    }
}

// NCBI unit test (598d8ae6): c++/src/objtools/blast/seqdb_writer/unit_test/writedb_unit_test.cpp:1826-1833 RegisterVariantsOfSameMaskingAlgorithm
#[test]
fn register_variants_of_same_masking_algorithm() {
    let mut registry = MaskInfoRegistry::default();
    let id1 = registry.add_program(FILTER_PROGRAM_SEG, b"", b"").unwrap();
    let id2 = registry
        .add_program(FILTER_PROGRAM_SEG, b"dummy", b"")
        .unwrap();
    assert_eq!(id1 + 1, id2);
}

/// `RegisterTooManyVariantsOfSameMaskingAlgorithm` (writedb_unit_test.cpp:1835-1882).
fn register_too_many_variants_of_same_masking_algorithm(masking_algo: i32, max_variants: usize) {
    let mut registry = MaskInfoRegistry::default();
    let mut algo_ids = Vec::new();
    for i in 0..max_variants * 2 {
        // For repeat and other there must be options, or the algorithm's value becomes its id.
        let options = if i == 0 && masking_algo < FILTER_PROGRAM_REPEAT {
            String::new()
        } else {
            i.to_string()
        };
        let result = registry.add_program(masking_algo, options.as_bytes(), b"");
        if i >= max_variants {
            assert!(is_write_db_exception(&result.unwrap_err()), "{i}");
        } else {
            algo_ids.push(result.unwrap());
        }
    }
    // The ids were assigned in increasing order.
    assert_eq!(algo_ids.len(), max_variants);
    for (i, id) in algo_ids.iter().enumerate() {
        assert_eq!(masking_algo + i as i32, *id);
    }
    // Only valid ids were assigned.
    for i in 0..max_variants * 2 {
        assert_eq!(
            registry.is_registered(masking_algo + i as i32),
            i < max_variants,
            "{i}"
        );
    }
}

// NCBI unit test (598d8ae6): c++/src/objtools/blast/seqdb_writer/unit_test/writedb_unit_test.cpp:1884-1889 RegisterTooManyVariantsOfDust
#[test]
fn register_too_many_variants_of_dust() {
    register_too_many_variants_of_same_masking_algorithm(
        FILTER_PROGRAM_DUST,
        (FILTER_PROGRAM_SEG - FILTER_PROGRAM_DUST) as usize,
    );
}

// NCBI unit test (598d8ae6): c++/src/objtools/blast/seqdb_writer/unit_test/writedb_unit_test.cpp:1891-1896 RegisterTooManyVariantsOfSeg
#[test]
fn register_too_many_variants_of_seg() {
    register_too_many_variants_of_same_masking_algorithm(
        FILTER_PROGRAM_SEG,
        (FILTER_PROGRAM_WINDOWMASKER - FILTER_PROGRAM_SEG) as usize,
    );
}

// NCBI unit test (598d8ae6): c++/src/objtools/blast/seqdb_writer/unit_test/writedb_unit_test.cpp:1898-1903 RegisterTooManyVariantsOfWindowMasker
#[test]
fn register_too_many_variants_of_window_masker() {
    register_too_many_variants_of_same_masking_algorithm(
        FILTER_PROGRAM_WINDOWMASKER,
        (FILTER_PROGRAM_REPEAT - FILTER_PROGRAM_WINDOWMASKER) as usize,
    );
}

// NCBI unit test (598d8ae6): c++/src/objtools/blast/seqdb_writer/unit_test/writedb_unit_test.cpp:1905-1910 RegisterTooManyVariantsOfRepeats
#[test]
fn register_too_many_variants_of_repeats() {
    register_too_many_variants_of_same_masking_algorithm(
        FILTER_PROGRAM_REPEAT,
        (FILTER_PROGRAM_OTHER - FILTER_PROGRAM_REPEAT) as usize,
    );
}

// NCBI unit test (598d8ae6): c++/src/objtools/blast/seqdb_writer/unit_test/writedb_unit_test.cpp:1912-1917 RegisterTooManyVariantsOfOther
#[test]
fn register_too_many_variants_of_other() {
    register_too_many_variants_of_same_masking_algorithm(
        FILTER_PROGRAM_OTHER,
        (FILTER_PROGRAM_MAX - FILTER_PROGRAM_OTHER) as usize,
    );
}

// NCBI unit test (598d8ae6): c++/src/objtools/blast/seqdb_writer/unit_test/writedb_unit_test.cpp:1919-2004 MaskDataColumn
#[test]
fn mask_data_column() {
    const NUM_SEQS: usize = 3;
    for version in [DbVersion::V4, DbVersion::V5] {
        let tmp = TempDir::new("ncbi-mask-data-column");
        let mut db = new_db(
            &tmp.base("mask-data-db"),
            b"Oct 12, 2026  3:22 PM",
            0,
            version,
        );
        let seg_id = db
            .register_mask_algorithm(FILTER_PROGRAM_SEG, b"", b"")
            .unwrap();
        let repeat_id = db
            .register_mask_algorithm(FILTER_PROGRAM_REPEAT, b"-species Desmodus_rotundus", b"")
            .unwrap();

        // The ranges set for each sequence, by algorithm.
        let mut expected: Vec<Vec<(i32, Vec<(u32, u32)>)>> = Vec::new();
        for i in 0..NUM_SEQS {
            // NCBI takes sequences of 1024 residues or more.
            db.add_sequence(bioseq(i + 1, 1100 + i)).unwrap();
            let mut ranges = Vec::new();
            if i & 1 != 0 {
                let offsets: Vec<(u32, u32)> = (0..i + 5)
                    .map(|j| {
                        let first = (i * 13 + j * 7 + 2) as u32;
                        (first, first + 3 + ((i + j) % 11) as u32)
                    })
                    .collect();
                ranges.push(MaskedRanges {
                    algorithm_id: seg_id,
                    offsets,
                });
            }
            if i & 2 != 0 {
                let offsets: Vec<(u32, u32)> = (0..i + 5)
                    .map(|j| {
                        let first = (i * 10 + j * 5 + 2) as u32;
                        (first, first + 20)
                    })
                    .collect();
                ranges.push(MaskedRanges {
                    algorithm_id: repeat_id,
                    offsets,
                });
            }
            // SetMaskData() is called when `i & 7` (not for the first sequence); calling it with
            // no ranges is the same as not calling it.
            if i & 7 != 0 {
                db.set_mask_data(&ranges).unwrap();
            }
            expected.push(
                ranges
                    .into_iter()
                    .map(|r| (r.algorithm_id, r.offsets))
                    .collect(),
            );
        }
        db.close().unwrap();

        // NCBI leaves the files untested; the column files are written, and the reader finds
        // what was set.
        let written = files_named(&tmp.0, "mask-data-db.");
        for extension in ["paa", "pab", "pac", "phr", "pin", "psq"] {
            assert!(
                written.contains(&format!("mask-data-db.{extension}")),
                "{version:?} {written:?}"
            );
        }
        let reader = SeqDb::open(&tmp.base("mask-data-db"), SeqType::Protein).unwrap();
        for (oid, sets) in expected.iter().enumerate() {
            for algo in [seg_id, repeat_id] {
                let want = sets
                    .iter()
                    .find(|(id, _)| *id == algo)
                    .map_or_else(Vec::new, |(_, offsets)| offsets.clone());
                assert_eq!(
                    reader
                        .mask_ranges_at(reader.locate_oid(oid as u32).unwrap(), algo)
                        .unwrap(),
                    MaskRanges::Ranges(want),
                    "{version:?} oid {oid} algorithm {algo}"
                );
            }
        }
    }
}

// NCBI unit test (598d8ae6): c++/src/objtools/blast/seqdb_writer/unit_test/writedb_unit_test.cpp:2006-2017 DuplicateAlgoId
#[test]
fn duplicate_algo_id() {
    let tmp = TempDir::new("ncbi-duplicate-algo-id");
    let mut db = new_db(
        &tmp.base("mask-data-db"),
        b"Oct 12, 2026  3:22 PM",
        0,
        DbVersion::V4,
    );
    db.register_mask_algorithm(FILTER_PROGRAM_SEG, b"", b"")
        .unwrap();
    let error = db
        .register_mask_algorithm(FILTER_PROGRAM_SEG, b"", b"")
        .unwrap_err();
    assert!(is_write_db_exception(&error), "{error}");
}

// NCBI unit test (598d8ae6): c++/src/objtools/blast/seqdb_writer/unit_test/writedb_unit_test.cpp:2019-2050 TooManyAlgoId
#[test]
fn too_many_algo_id() {
    let tmp = TempDir::new("ncbi-too-many-algo-id");
    let mut db = new_db(
        &tmp.base("mask-data-db"),
        b"Oct 12, 2026  3:22 PM",
        0,
        DbVersion::V4,
    );
    let masking_algorithm = FILTER_PROGRAM_SEG;
    let max_variants = (FILTER_PROGRAM_WINDOWMASKER - masking_algorithm) as usize;
    let mut algo_ids = Vec::new();
    // The last ones fail.
    for i in 0..max_variants * 2 {
        let options = if i == 0 { String::new() } else { i.to_string() };
        let result = db.register_mask_algorithm(masking_algorithm, options.as_bytes(), b"");
        if i >= max_variants {
            assert!(is_write_db_exception(&result.unwrap_err()), "{i}");
        } else {
            algo_ids.push(result.unwrap());
        }
    }
    // The ids were assigned in increasing order.
    assert_eq!(algo_ids.len(), max_variants);
    for (i, id) in algo_ids.iter().enumerate() {
        assert_eq!(masking_algorithm + i as i32, *id);
    }
}

// NCBI unit test (598d8ae6): c++/src/objtools/blast/seqdb_writer/unit_test/writedb_unit_test.cpp:2052-2089 UndefinedAlgoID
#[test]
fn undefined_algo_id() {
    let tmp = TempDir::new("ncbi-undefined-algo-id");
    let mut db = new_db(
        &tmp.base("mask-data-db"),
        b"Oct 12, 2026  3:22 PM",
        0,
        DbVersion::V4,
    );
    db.register_mask_algorithm(FILTER_PROGRAM_SEG, b"", b"")
        .unwrap();
    db.register_mask_algorithm(FILTER_PROGRAM_SEG, b"-species Aotus_vociferans", b"")
        .unwrap();
    db.register_mask_algorithm(FILTER_PROGRAM_REPEAT, b"-species Desmodus_rotundus", b"")
        .unwrap();

    let length = 300;
    db.add_sequence(bioseq(1, length)).unwrap();
    // `dust` was never registered.
    let ranges = [MaskedRanges {
        algorithm_id: FILTER_PROGRAM_DUST,
        offsets: vec![(length as u32 / 3, length as u32)],
    }];
    let error = db.set_mask_data(&ranges).unwrap_err();
    assert!(is_write_db_exception(&error), "{error}");
    assert!(
        error
            .to_string()
            .contains("Algorithm IDs must be registered before use"),
        "{error}"
    );

    db.close().unwrap();
}

// NCBI unit test (598d8ae6): c++/src/objtools/blast/seqdb_writer/unit_test/writedb_unit_test.cpp:2091-2126 MaskDataBoundsError
#[test]
fn mask_data_bounds_error() {
    let tmp = TempDir::new("ncbi-mask-bounds");
    let mut db = new_db(
        &tmp.base("mask-data-db"),
        b"Oct 12, 2026  3:22 PM",
        0,
        DbVersion::V4,
    );
    db.register_mask_algorithm(FILTER_PROGRAM_SEG, b"", b"")
        .unwrap();
    db.register_mask_algorithm(FILTER_PROGRAM_SEG, b"-species Aotus_vociferans", b"")
        .unwrap();
    db.register_mask_algorithm(FILTER_PROGRAM_REPEAT, b"-species Desmodus_rotundus", b"")
        .unwrap();

    let length = 300;
    db.add_sequence(bioseq(1, length)).unwrap();
    // The range ends one past the sequence (the algorithm `dust` is not registered either, as
    // in NCBI's case, which therefore stops at the unregistered id).
    let ranges = [MaskedRanges {
        algorithm_id: FILTER_PROGRAM_DUST,
        offsets: vec![(length as u32 / 3, length as u32 + 1)],
    }];
    let error = db.set_mask_data(&ranges).unwrap_err();
    assert!(is_write_db_exception(&error), "{error}");

    db.close().unwrap();
}
