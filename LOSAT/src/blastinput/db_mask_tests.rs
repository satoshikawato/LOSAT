//! Tests of the database masks' resolver (`db_mask.rs`), each case a probe of the M4a-2 survey
//! (`M/m4a2/SURVEY.md` §1, NCBI BLAST+ 2.17.0) on the databases of `tests/db_search/dbmask`, which
//! NCBI's makeblastdb wrote from the same inputs (README there): `a5` (dust 11, dust
//! `-level 40 -window 32` 12, windowmasker 30), `b5` (the second dust alone, its id 11), `k5`
//! (`-mask_id mydust,mywm`: 100, 101), `nm5` (no masks); `tests/seqdb/db/mask_nucl_colon_id`
//! (`-mask_id c:id`, stored with 0x01 for the colon).

use std::path::{Path, PathBuf};

use super::*;
use crate::seqdb::{SearchedPlace, SeqType};

fn dbmask(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/db_search/dbmask")
        .join(name)
}

fn open(names: &[PathBuf]) -> (SeqDb, String) {
    let name = names
        .iter()
        .map(|path| path.to_str().unwrap())
        .collect::<Vec<_>>()
        .join(" ");
    let mut db = SeqDb::open(&name, SeqType::Nucleotide).unwrap();
    db.read_name_records(&name, true).unwrap();
    (db, name)
}

fn soft(text: &str) -> DbMaskOption {
    DbMaskOption {
        kind: MaskKind::Soft,
        algorithm: FilteringAlgorithm::parse(text),
    }
}

/// NCBI's error of `result` as the command line writes it: (exit status, text).
fn native(result: Result<impl std::fmt::Debug>) -> (i32, String) {
    let error = result.unwrap_err();
    let native = error
        .downcast_ref::<NativeError>()
        .unwrap_or_else(|| panic!("not NCBI's error: {error:#}"));
    (
        native.exit,
        String::from_utf8(native.message.clone()).unwrap(),
    )
}

/// `resolve` with its warnings collected.
fn resolve(
    option: &DbMaskOption,
    db: &SeqDb,
    name: &str,
) -> (Result<Option<ResolvedMask>>, String) {
    let mut warnings = Vec::new();
    let resolved = option.resolve(db, name, "BLASTN", &mut |text| {
        warnings.extend_from_slice(text);
        Ok(())
    });
    (resolved, String::from_utf8(warnings).unwrap())
}

// Probes o01, o08_* (`+11`, `011`: the id 11) and the names ` 11`, `11 `, `11x`, `""` and
// `2147483648` (o08_*: NCBI's warning).
#[test]
fn an_int4_in_ncbi_s_form_is_an_id_and_any_other_text_a_name() {
    for (text, id) in [
        ("11", 11),
        ("+11", 11),
        ("011", 11),
        ("-1", -1),
        ("0", 0),
        ("-5", -5),
    ] {
        assert_eq!(
            FilteringAlgorithm::parse(text),
            FilteringAlgorithm::Id(id),
            "{text:?}"
        );
    }
    for text in [" 11", "11 ", "11x", "", "2147483648", "dust", "-"] {
        assert_eq!(
            FilteringAlgorithm::parse(text),
            FilteringAlgorithm::Name(text.to_string()),
            "{text:?}"
        );
    }
}

// blast_args.cpp:2517-2523: the soft option first (the parser refuses both together).
#[test]
fn the_soft_option_comes_first() {
    let args = DbMaskArgs {
        db_soft_mask: Some("dust".to_string()),
        db_hard_mask: Some("11".to_string()),
    };
    assert_eq!(
        args.option(),
        Some(DbMaskOption {
            kind: MaskKind::Soft,
            algorithm: FilteringAlgorithm::Name("dust".to_string()),
        })
    );
    let args = DbMaskArgs {
        db_soft_mask: None,
        db_hard_mask: Some("+30".to_string()),
    };
    assert_eq!(
        args.option(),
        Some(DbMaskOption {
            kind: MaskKind::Hard,
            algorithm: FilteringAlgorithm::Id(30),
        })
    );
    assert_eq!(DbMaskArgs::default().option(), None);
}

// Probes o02 (`99`), o07 (`11` on a database without masks): NCBI's options error when it opens
// the database; o01, o03 (a name), o08 (`0`, `-5`): no check there.
#[test]
fn an_id_above_0_must_be_the_database_s() {
    let (a5, a5_name) = open(&[dbmask("a5")]);
    for text in ["11", "+12", "30", "dust", "0", "-1", "-5", "99x"] {
        soft(text)
            .validate_at_open(&a5, &a5_name, true, "BLASTN")
            .unwrap_or_else(|error| panic!("{text}: {error:#}"));
    }
    for (text, id) in [("99", 99), ("10", 10)] {
        assert_eq!(
            native(soft(text).validate_at_open(&a5, "db/a5", true, "BLASTN")),
            (
                1,
                format!("LOSAT options error: Masking algorithm ID {id} is not supported in nucleotide 'db/a5' BLAST database\n")
            )
        );
    }
    let (nm5, _) = open(&[dbmask("nm5")]);
    assert_eq!(
        native(soft("11").validate_at_open(&nm5, "db/nm5", true, "BLASTN")),
        (
            1,
            "LOSAT options error: Masking algorithm ID 11 is not supported in nucleotide 'db/nm5' BLAST database\n"
                .to_string()
        )
    );
}

// Probes o03 (`dust`), k05 (`MYDUST`: names are case-sensitive), o08_v[-1] (`-1`), o07c (a name on
// `nm5`): NCBI's warning and a search without masks; o03b (the stored description as a name)
// and k01/k03 (`-mask_id`'s name, its id).
#[test]
fn a_name_is_the_stored_description_or_the_mask_id() {
    let (a5, a5_name) = open(&[dbmask("a5")]);
    for text in ["dust", "-1", "DUST", "", "11 "] {
        let (resolved, warning) = resolve(&soft(text), &a5, "db/a5");
        assert_eq!(resolved.unwrap(), None, "{text:?}");
        assert_eq!(
            warning,
            "Warning: [blastn] Subject mask not found in db/a5, proceeding without subject masking.\n",
            "{text:?}"
        );
    }
    let (resolved, warning) = resolve(&soft("10:window=64; level=20; linker=1"), &a5, &a5_name);
    let resolved = resolved.unwrap().unwrap();
    assert_eq!(warning, "");
    assert_eq!(
        (resolved.id, resolved.key.as_deref()),
        (11, Some("10:window=64; level=20; linker=1"))
    );
    assert_eq!(resolved.shown.len(), 1);
    assert_eq!(
        (
            resolved.shown[0].program_name.as_str(),
            resolved.shown[0].options.as_str()
        ),
        ("dust", "window=64; level=20; linker=1")
    );
    let (resolved, _) = resolve(&soft("+30"), &a5, &a5_name);
    let resolved = resolved.unwrap().unwrap();
    assert_eq!((resolved.id, resolved.key.as_deref()), (30, None));
    assert_eq!(
        (
            resolved.shown[0].program_name.as_str(),
            resolved.shown[0].options.as_str()
        ),
        ("windowmasker", "")
    );

    let (k5, k5_name) = open(&[dbmask("k5")]);
    for (text, id, key) in [
        ("mydust", 100, Some("mydust")),
        ("100", 100, None),
        ("mywm", 101, Some("mywm")),
    ] {
        let (resolved, warning) = resolve(&soft(text), &k5, &k5_name);
        let resolved = resolved.unwrap().unwrap();
        assert_eq!(warning, "");
        assert_eq!((resolved.id, resolved.key.as_deref()), (id, key));
    }
    let (resolved, warning) = resolve(&soft("MYDUST"), &k5, "db/k5");
    assert_eq!(resolved.unwrap(), None);
    assert!(warning.contains("Subject mask not found in db/k5"));
    let (nm5, _) = open(&[dbmask("nm5")]);
    let (resolved, warning) = resolve(&soft("dust"), &nm5, "db/nm5");
    assert_eq!(resolved.unwrap(), None);
    assert!(warning.contains("Subject mask not found in db/nm5"));
}

// Probes c01 (`c:id`: the stored id has 0x01 for the colon, NCBI's warning) and c06 (`c\x01id`
// finds it; the report shows the colon, `s_RestoreColon`).
#[test]
fn a_mask_id_with_a_colon_is_found_by_its_stored_form() {
    let colon = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/db/mask_nucl_colon_id");
    let (db, name) = open(&[colon]);
    let (resolved, warning) = resolve(&soft("c:id"), &db, &name);
    assert_eq!(resolved.unwrap(), None);
    assert!(warning.contains("Subject mask not found"));
    let (resolved, warning) = resolve(&soft("c\u{1}id"), &db, &name);
    let resolved = resolved.unwrap().unwrap();
    assert_eq!(warning, "");
    assert_eq!(resolved.id, 100);
    assert_eq!(resolved.shown[0].program_name, "c:id");
}

// Probes o08_v[0] (`0`), `-5`: the formatter's `GetMaskAlgorithmDetails` fails with the
// database's table, exit 2.
#[test]
fn ids_0_and_below_minus_1_are_the_formatter_s_database_error() {
    let (a5, a5_name) = open(&[dbmask("a5")]);
    let table = "\nAvailable filtering algorithms applied to database sequences:\n\n\
                 Algorithm ID Algorithm name                          Algorithm options                       \n\
                 11           dust                                    window=64; level=20; linker=1           \n\
                 12           dust                                    window=32; level=40; linker=1           \n\
                 30           windowmasker                            default options used                    \n";
    for id in [0, -5] {
        let (resolved, warning) = resolve(&soft(&id.to_string()), &a5, &a5_name);
        assert_eq!(warning, "");
        assert_eq!(
            native(resolved),
            (2, format!("LOSAT database error: Filtering algorithm ID {id} is not supported.\n{table}\n"))
        );
    }
}

// Probes m01, m02 (`a5 nm5`, `nm5 a5` with 11: `nm5` has no table), m03 (`a5 b5` with 12: `b5`'s
// table has its own 11 only), m04 (`b5 a5` with 11: each name's own algorithm 11).
#[test]
fn each_name_of_the_list_needs_the_id_in_its_own_database() {
    for names in [[dbmask("a5"), dbmask("nm5")], [dbmask("nm5"), dbmask("a5")]] {
        let (db, name) = open(&names);
        let (resolved, _) = resolve(&soft("11"), &db, &name);
        assert_eq!(
            native(resolved),
            (
                2,
                "LOSAT database error: Filtering algorithm ID 11 is not supported.\n\n".to_string()
            )
        );
    }
    let (db, name) = open(&[dbmask("a5"), dbmask("b5")]);
    let (resolved, _) = resolve(&soft("12"), &db, &name);
    assert_eq!(
        native(resolved).1,
        "LOSAT database error: Filtering algorithm ID 12 is not supported.\n\n\
         Available filtering algorithms applied to database sequences:\n\n\
         Algorithm ID Algorithm name                          Algorithm options                       \n\
         11           dust                                    window=32; level=40; linker=1           \n\n"
    );
    let (db, name) = open(&[dbmask("b5"), dbmask("a5")]);
    let (resolved, _) = resolve(&soft("11"), &db, &name);
    let resolved = resolved.unwrap().unwrap();
    let shown: Vec<_> = resolved
        .shown
        .iter()
        .map(|algorithm| (algorithm.program_name.as_str(), algorithm.options.as_str()))
        .collect();
    assert_eq!(
        shown,
        [
            ("dust", "window=32; level=40; linker=1"),
            ("dust", "window=64; level=20; linker=1")
        ]
    );
    // The search's table puts `a5` first (NCBI sorts the volumes by name): its 11 is `a5`'s
    // dust, which `b5`'s volume does not have (NCBI's stale volume-algorithm cache, M/DECISIONS
    // 9): LOSAT refuses the search (decision D5).
    let error = resolved
        .check_losat_limits(&db, "BLASTN", HardScan::Nucleotide { scanned_word: 11 })
        .unwrap_err();
    assert!(
        format!("{error:#}")
            .contains("the volume has mask data without the filtering algorithm 11"),
        "{error:#}"
    );
    // The alias of the same two (probe m05: `al_a_nm`; m06: `al_a_b` with 12) is one name.
    let (db, name) = open(&[dbmask("al_a_nm")]);
    let (resolved, _) = resolve(&soft("11"), &db, &name);
    let resolved = resolved.unwrap().unwrap();
    assert_eq!(resolved.shown.len(), 1);
    db.check_search_masks(resolved.id).unwrap();
}

// search_strategy.cpp:499-516: the key of a name (probes a04, k02), else the id (a03: `+11` is
// 11), then the masking type (blast_def.h:235-239).
#[test]
fn the_archive_writes_the_key_of_a_name_or_the_id() {
    let mut mask = ResolvedMask {
        kind: MaskKind::Soft,
        id: 11,
        key: None,
        shown: Vec::new(),
    };
    let params = |mask: &ResolvedMask| {
        mask.archive_params()
            .into_iter()
            .map(|param| (param.name, param.value))
            .collect::<Vec<_>>()
    };
    assert_eq!(
        params(&mask),
        [
            ("DbFilteringAlgorithmId", ParamValue::Integer(11)),
            ("SubjectMaskingType", ParamValue::Integer(1)),
        ]
    );
    mask.key = Some("mydust".to_string());
    mask.kind = MaskKind::Hard;
    assert_eq!(
        params(&mask),
        [
            (
                "DbFilteringAlgorithmKey",
                ParamValue::Str("mydust".to_string())
            ),
            ("SubjectMaskingType", ParamValue::Integer(2)),
        ]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:241-268 and 579-584 (quoted
// at `HardMaskShape::of_scan`)
// A search with `-db_hard_mask` is refused, before anything is written, at the first searched
// sequence whose masks NCBI's hard chunks do not handle (series NP, M4a-2b decisions H3-H5): a5's
// s4 (OID 3) is masked whole ([0, 300): the "no chunking" branch scans the full range of a chunk
// of no letters), refused at blastn-short's word 7 and searched at 11 and more (M4a-2c,
// `TAIL_MASK_SEARCHED_WORD`); edge1's algorithm 11 holds s5's masks out of order, 12 overlapping
// ones, at every word; 13 (adjacent masks: a hard range of no letters), 14 (masks at both ends)
// and 15 (one-letter masks) are searched; edge2's 14 (a tail mask) and 16 (s4 masked whole, OID 3
// before s5) are refused at 7 and searched at 11. `-db_soft_mask` keeps every shape.
#[test]
fn a_hard_masked_search_refuses_the_masks_that_ncbi_s_chunks_do_not_handle() {
    let mask = |kind, id| ResolvedMask {
        kind,
        id,
        key: None,
        shown: Vec::new(),
    };
    let refusal = |db: &SeqDb, id, word| {
        mask(MaskKind::Hard, id)
            .check_losat_limits(db, "BLASTN", HardScan::Nucleotide { scanned_word: word })
            .map_err(|error| error.to_string())
    };
    let (a5, _) = open(&[dbmask("a5")]);
    let error = refusal(&a5, 11, 7).unwrap_err();
    assert!(
        error.starts_with("-db_hard_mask: the sequence of OID 3 (volume ")
            && error.contains(
                "a5) has a mask, [0, 300), running to its end after the chunk from its start, for \
                 which NCBI BLAST+ 2.17.0 scans the masked letters past that chunk while it \
                 extends within the chunk (at a word of fewer than 11 letters: rows that differ \
                 from a near input's, or a crash)"
            )
            && error.ends_with("LOSAT's BLASTN does not search such a database with -db_hard_mask"),
        "{error}"
    );
    assert!(refusal(&a5, 11, 10).is_err());
    for word in [11, 16, 28] {
        refusal(&a5, 11, word).unwrap();
    }
    mask(MaskKind::Soft, 11)
        .check_losat_limits(&a5, "BLASTN", HardScan::Nucleotide { scanned_word: 7 })
        .unwrap();
    let (edge1, _) = open(&[dbmask("edge1")]);
    for word in [7, 11] {
        let error = refusal(&edge1, 11, word).unwrap_err();
        assert!(
            error.contains("OID 4 ")
                && error.contains("has masks out of order ([100, 200) after [600, 700))"),
            "{error}"
        );
        let error = refusal(&edge1, 12, word).unwrap_err();
        assert!(error.contains("([200, 400) after [100, 300))"), "{error}");
    }
    for id in [13, 14, 15] {
        refusal(&edge1, id, 11).unwrap();
    }
    let (edge2, _) = open(&[dbmask("edge2")]);
    let error = refusal(&edge2, 14, 7).unwrap_err();
    assert!(
        error.contains("OID 4 ") && error.contains("[1000, 1200), running to its end"),
        "{error}"
    );
    let error = refusal(&edge2, 16, 7).unwrap_err();
    assert!(
        error.contains("OID 3 ") && error.contains("[0, 300), running"),
        "{error}"
    );
    for id in [14, 15, 16] {
        refusal(&edge2, id, 11).unwrap();
    }
    for id in 11..=16 {
        mask(MaskKind::Soft, id)
            .check_losat_limits(&edge2, "BLASTN", HardScan::Nucleotide { scanned_word: 7 })
            .unwrap();
    }
}

// The shapes of `HardMaskShape::of_scan`, read off the chunks (blast_engine.c:221-276, 579-584, 1354,
// 1368-1370; na_ungapped.c:1647-1668): a mask that runs to the end after the chunk from 0 (empty
// masks at 0 or at the end around it; the first range at or below MAX_DBSEQ_LEN), in a subject of
// a word or more (`min_subj_seq_length` is 1 for BLASTN: 4 letters at word 4, review F2), at a
// word below `TAIL_MASK_SEARCHED_WORD` (series NP, M4a-2c: searched from 11 on); masks
// that go back; a last mask after which the final chunk starts before the end of an earlier chunk
// of another range, both chunks scanning a word (masksubj.inl:43-58; blast_engine.c:500, 565;
// blast_hits.c:2873-2911, 2956: NCBI's merge then takes the earlier chunk's HSPs by the overlap it
// reads past its list; review 2 F1: not where only a chunk without a word ends there, or where
// the final chunk has no word), not where the last range is longer than MAX_DBSEQ_LEN and its
// final chunk follows a full chunk of the same range (review F3). The length is read for every
// sequence with sorted masks (the chunks need it), never for masks out of order.
#[test]
fn the_hard_mask_shapes_follow_ncbi_s_chunks() {
    use std::cell::Cell;
    let shape_at = |masks: &[(u32, u32)], length: usize, word: usize| {
        let reads = Cell::new(0);
        let read = || {
            reads.set(reads.get() + 1);
            Some(length)
        };
        (shape_of(masks, None, &read, word).unwrap(), reads.get())
    };
    let shape = |masks: &[(u32, u32)], length: usize| shape_at(masks, length, 11);
    let tail = |start, end| Some(HardMaskShape::TailMask { start, end });
    // A mask that runs to the end after the chunk from 0 is refused below the word
    // `TAIL_MASK_SEARCHED_WORD` (blastn-short's 7) and searched from it on (series NP, M4a-2c);
    // the length is read all the same.
    let tails: [(&[(u32, u32)], usize, Option<HardMaskShape>); 7] = [
        (&[(100, 300)], 300, tail(100, 300)),
        (&[(0, 300)], 300, tail(0, 300)),
        (&[(100, 300), (300, 300)], 300, tail(100, 300)),
        // Empty masks at 0 before it (review F1): the chunk after them starts at 0 as well.
        (&[(0, 0), (1000, 1200)], 1200, tail(1000, 1200)),
        (&[(0, 0), (0, 0), (0, 1200)], 1200, tail(0, 1200)),
        (
            &[(5_000_000, 6_000_000)],
            6_000_000,
            tail(5_000_000, 6_000_000),
        ),
        // The shortest subject the scan of the full range reads a word of at 11.
        (&[(5, 11)], 11, tail(5, 11)),
    ];
    for (masks, length, refused) in tails {
        assert_eq!(shape_at(masks, length, 7), (refused, 1), "{masks:?}");
        assert_eq!(
            shape_at(masks, length, TAIL_MASK_SEARCHED_WORD - 1),
            (refused, 1),
            "{masks:?}"
        );
        for word in [TAIL_MASK_SEARCHED_WORD, 16, 28] {
            assert_eq!(shape_at(masks, length, word), (None, 1), "{masks:?} {word}");
        }
    }
    assert_eq!(shape_at(&[(0, 0)], 1200, 7), (None, 1));
    // Above MAX_DBSEQ_LEN the first chunk is a full chunk: NCBI chunks it as any other.
    assert_eq!(shape_at(&[(5_000_001, 6_000_000)], 6_000_000, 7), (None, 1));
    // An empty mask at the end, a mask before the end, two masks.
    assert_eq!(shape(&[(300, 300)], 300), (None, 1));
    assert_eq!(shape(&[(100, 299)], 300), (None, 1));
    assert_eq!(shape(&[(0, 100), (200, 300)], 300), (None, 1));
    let unsorted = |before, mask| Some(HardMaskShape::Unsorted { before, mask });
    assert_eq!(
        shape(&[(600, 700), (100, 200)], 1200),
        (unsorted((600, 700), (100, 200)), 0)
    );
    assert_eq!(
        shape(&[(100, 300), (200, 400)], 1200),
        (unsorted((100, 300), (200, 400)), 0)
    );
    assert_eq!(shape(&[(50, 40)], 1200), (unsorted((50, 40), (50, 40)), 0));
    assert_eq!(shape(&[(100, 200), (200, 300)], 1200), (None, 1));
    let short = |start, end| Some(HardMaskShape::ShortLastMask { start, end });
    // 1002 % 4 = 2 letters before the end in the final chunk, a mask of 1.
    assert_eq!(shape(&[(1001, 1002)], 1200), (short(1001, 1002), 1));
    assert_eq!(shape(&[(1000, 1002)], 1200), (None, 1));
    assert_eq!(shape(&[(1003, 1003)], 1200), (short(1003, 1003), 1));
    // The last range over MAX_DBSEQ_LEN letters (review F3): the chunk after the mask is a full
    // chunk with an overlap read in bounds, and the final chunk follows a chunk of its own range
    // (decision B-d); at exactly MAX_DBSEQ_LEN letters after the offset it is the final chunk.
    assert_eq!(shape(&[(1001, 1002)], 5_300_000), (None, 1));
    assert_eq!(shape(&[(1001, 1002)], 5_001_100), (None, 1));
    assert_eq!(shape(&[(1001, 1002)], 5_001_000), (short(1001, 1002), 1));
    // The last mask at the end: the last hard range is empty and the loop ends before it.
    assert_eq!(shape(&[(10, 20), (1199, 1199)], 1199), (None, 1));
    // A middle mask shorter than its residual merges with an overlap that NCBI reads in bounds.
    assert_eq!(shape(&[(301, 302), (500, 503)], 1200), (None, 1));
    // Review 2 F1: the chunk between two one-letter masks (101-102, scanned from 101: one letter)
    // and the chunk of an empty last mask's range (1000-1002, scanned from 1002) hold no word, and
    // the chunk before them ends at the final chunk's offset: no HSP of theirs reaches the merge.
    assert_eq!(shape(&[(100, 101), (102, 103)], 1200), (None, 1));
    assert_eq!(shape(&[(1000, 1002), (1002, 1002)], 1200), (None, 1));
    // The chunk from 0 ends at 101, after the final chunk's offset 100.
    assert_eq!(shape(&[(101, 102), (102, 103)], 1200), (short(102, 103), 1));
    // The final chunk (1000-1012, scanned from 1002) holds a word of 10 letters, not of 11.
    assert_eq!(
        shape_of(&[(1001, 1002)], None, &|| Some(1012), 11),
        Ok(None)
    );
    assert_eq!(
        shape_of(&[(1001, 1002)], None, &|| Some(1012), 10),
        Ok(short(1001, 1002))
    );
    // A subject that the scan reads no word of: shorter than the word.
    assert_eq!(shape_at(&[(0, 6)], 6, 7), (None, 1));
    assert_eq!(shape_at(&[(0, 7)], 7, 7), (tail(0, 7), 1));
    assert_eq!(shape_of(&[(0, 4)], None, &|| Some(4), 5), Ok(None));
    assert_eq!(shape_of(&[(2, 4)], None, &|| Some(4), 4), Ok(tail(2, 4)));
    assert_eq!(shape_of(&[(0, 4)], None, &|| Some(4), 4), Ok(tail(0, 4)));
    assert_eq!(shape_of(&[(1, 3)], None, &|| Some(3), 4), Ok(None));
    // The length that cannot be read ends the check.
    assert_eq!(shape_of(&[(100, 300)], None, &|| None, 11), Err(()));
    // The bounds of the length that the index gives (review 2 S4): the length is read only where
    // the shape differs between them (a mask that runs to the end, here at 300 or at 1201 of
    // 1200-1203), not where it is the same at all four (none, or a refusal read off the chunks).
    let bounded_at = |masks: &[(u32, u32)], bounds: (usize, usize), length: usize, word| {
        let reads = Cell::new(0);
        let read = || {
            reads.set(reads.get() + 1);
            Some(length)
        };
        (
            shape_of(masks, Some(bounds), &read, word).unwrap(),
            reads.get(),
        )
    };
    let bounded = |masks: &[(u32, u32)], bounds: (usize, usize), length: usize| {
        bounded_at(masks, bounds, length, 11)
    };
    assert_eq!(
        bounded_at(&[(100, 300)], (300, 303), 300, 7),
        (tail(100, 300), 1)
    );
    assert_eq!(bounded(&[(100, 300)], (300, 303), 300), (None, 1));
    assert_eq!(
        bounded_at(&[(1000, 1201)], (1200, 1203), 1201, 7),
        (tail(1000, 1201), 1)
    );
    assert_eq!(bounded(&[(1000, 1201)], (1200, 1203), 1201), (None, 1));
    assert_eq!(bounded(&[(100, 200)], (1200, 1203), 1201), (None, 0));
    assert_eq!(
        bounded(&[(1001, 1002)], (1200, 1203), 1202),
        (short(1001, 1002), 0)
    );
    // The final chunk 1000-1012 to 1000-1015 holds a word of 11 letters from 1013 on.
    assert_eq!(bounded(&[(1001, 1002)], (1012, 1015), 1012), (None, 1));
}

// blast_app_util.cpp:644-651 and local_db_adapter.cpp:54-58 (cited at `strategy_params`): after
// the formatter's warning the search strategy still writes a name as its key (probe of
// `-export_search_strategy` with `-db_soft_mask dust`; frozen as the strategy case
// blastn.db.dbmask.name), not the empty name or the id -1.
#[test]
fn the_strategy_keeps_a_name_the_database_does_not_know() {
    let names = |params: Vec<Param>| -> Vec<&'static str> {
        params.iter().map(|param| param.name).collect()
    };
    assert_eq!(
        names(soft("dust").strategy_params(None)),
        ["DbFilteringAlgorithmKey", "SubjectMaskingType"]
    );
    assert!(soft("").strategy_params(None).is_empty());
    assert!(soft("-1").strategy_params(None).is_empty());
    let mask = ResolvedMask {
        kind: MaskKind::Soft,
        id: 11,
        key: None,
        shown: Vec::new(),
    };
    assert_eq!(
        names(soft("+11").strategy_params(Some(&mask))),
        ["DbFilteringAlgorithmId", "SubjectMaskingType"]
    );
}

/// `open` of `names` restricted to the text seqid list `ids` (`-seqidlist`), as
/// `restrict_database` restricts a search's database before it reads the names' records.
fn open_listed(names: &[PathBuf], ids: &[&str]) -> (SeqDb, String) {
    let name = names
        .iter()
        .map(|path| path.to_str().unwrap())
        .collect::<Vec<_>>()
        .join(" ");
    let text: Vec<u8> = ids
        .iter()
        .flat_map(|id| format!("{id}\n").into_bytes())
        .collect();
    let file = crate::seqdb::seqidlist::read_si_list(&text, &mut Vec::new()).unwrap();
    let list = crate::seqdb::DbList::Seqids(crate::seqdb::SeqidList {
        negative: false,
        file,
    });
    let mut db = SeqDb::open(&name, SeqType::Nucleotide).unwrap();
    db.restrict(&list, &mut Vec::new()).unwrap();
    db.read_name_records(&name, true).unwrap();
    (db, name)
}

// Review 1 F2 of M4a-2a (probes d5_sid, d5_sid_n4, d5_nosid): `zz5` (the subjects as `t1`-`t5`,
// dust `-level 40 -window 32` alone: its own 11) lacks the list's algorithm 11 (`p5`'s dust; the
// list's 12 is zz5's). NCBI reads a sequence's masks only when it searches it (seqsrc_seqdb.cpp:
// 318-321) and translates the algorithm of a volume only for a blob that is not empty
// (seqdbimpl.cpp:2531-2546): a search that reads no blob of zz5 is NCBI's defined search (rc 0,
// 4 rows), one that reads one of them meets the stale volume-algorithm cache (decision D5).
#[test]
fn a_volume_without_the_algorithm_is_refused_only_where_the_search_reads_its_masks() {
    let names = [dbmask("p5"), dbmask("zz5")];
    let refused = |db: &SeqDb, name: &str| {
        let (resolved, _) = resolve(&soft("11"), db, name);
        let resolved = resolved.unwrap().unwrap();
        // Each name's own 11: p5's dust, zz5's dust40.
        let shown: Vec<_> = resolved
            .shown
            .iter()
            .map(|algorithm| algorithm.options.as_str())
            .collect();
        assert_eq!(
            shown,
            [
                "window=64; level=20; linker=1",
                "window=32; level=40; linker=1"
            ]
        );
        resolved
            .check_losat_limits(db, "BLASTN", HardScan::Nucleotide { scanned_word: 11 })
            .err()
            .map(|error| format!("{error:#}"))
    };
    // p5's sequences alone (probe d5_sid), and zz5's t2 and t5, which have no mask data.
    for ids in [&["s1", "s2", "s3"][..], &["s1", "t2", "t5"], &["t2"]] {
        let (db, name) = open_listed(&names, ids);
        assert_eq!(refused(&db, &name), None, "{ids:?}");
    }
    // A searched sequence of zz5 with mask data (t1, t3, t4), or the whole list (probe d5_nosid).
    for ids in [&["s1", "t1"][..], &["t4"]] {
        let (db, name) = open_listed(&names, ids);
        let error = refused(&db, &name).expect("refused");
        assert!(
            error.contains("zz5: unsupported database: the volume has mask data without the filtering algorithm 11"),
            "{ids:?}: {error}"
        );
    }
    let (db, name) = open(&names);
    assert!(refused(&db, &name).is_some());
}

// Review 1 S1 and T6 of M4a-2a: the formatter's algorithm of each name comes from the table that
// the name's record keeps (`read_name_records` opens each name once), one per record, and the
// pairwise report's `Masked using:` lines pair them with the names' `Database:` lines.
#[test]
fn each_name_s_masked_using_line_comes_from_its_own_record() {
    let (db, name) = open(&[dbmask("b5"), dbmask("a5"), dbmask("b5")]);
    assert!(db
        .name_records()
        .iter()
        .all(|record| matches!(record.mask_algorithms, Some(Ok(_)))));
    let (resolved, _) = resolve(&soft("11"), &db, &name);
    let resolved = resolved.unwrap().unwrap();
    let pairwise =
        crate::algorithm::common::database::pairwise_database_masked(&db, Some(&resolved.shown));
    let lines: Vec<_> = pairwise
        .databases
        .iter()
        .map(|database| {
            (
                String::from_utf8(database.definition.clone()).unwrap(),
                database.filt_algorithm_options.as_str(),
            )
        })
        .collect();
    assert_eq!(
        lines,
        [
            ("dust40 only".to_string(), "window=32; level=40; linker=1"),
            ("mask probe v5".to_string(), "window=64; level=20; linker=1"),
            ("dust40 only".to_string(), "window=32; level=40; linker=1"),
        ]
    );
    // Records read without the tables (another program's search) have none to give.
    let mut plain = SeqDb::open(&name, SeqType::Nucleotide).unwrap();
    plain.read_name_records(&name, false).unwrap();
    assert!(matches!(
        plain.name_mask_details(11),
        Err(crate::seqdb::NameMaskError::NotRead)
    ));
}

#[test]
#[should_panic(expected = "one masking algorithm per name of the -db list")]
fn a_masked_using_list_of_another_length_than_the_names_is_a_bug() {
    let (db, name) = open(&[dbmask("b5"), dbmask("a5")]);
    let (resolved, _) = resolve(&soft("11"), &db, &name);
    let resolved = resolved.unwrap().unwrap();
    crate::algorithm::common::database::pairwise_database_masked(&db, Some(&resolved.shown[..1]));
}

/// `HardMaskShape::of_scan` of BLASTN's scan, whose words have `scanned_word` letters.
fn shape_of(
    masks: &[(u32, u32)],
    bounds: Option<(usize, usize)>,
    length: &dyn Fn() -> Option<usize>,
    scanned_word: usize,
) -> Result<Option<HardMaskShape>, ()> {
    HardMaskShape::of_scan(
        masks,
        bounds,
        length,
        HardScan::Nucleotide { scanned_word },
        SearchedPlace::default(),
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-830 (the frames' ranges
// mapped from the bases' by x / 3, the minus frames reversed on frame -1's length), :221-276
// (the "no chunking" branch of a chunk from 0 whose next range starts at the frame's end) and
// :145-215 (the subject's list of ranges restored after the frame).
// A mask is refused where, after the mapping, a frame's last range that starts at 0 (not cut at
// `MAX_DBSEQ_LEN`, with no range after it or one at the frame's end) is not the frame's first range
// (`FrameStartMask`: NCBI frees memory that it did not allocate and crashes): a mask within the
// first two bases or within the one or two bases after the last whole codon (with no other mask in
// the frame). Where that range is the frame's first and shorter than the frame (a tail mask in
// frame +1, a head mask in frame -1, a mask that ends 2 or 3 bases before the end in frame +2:
// NCBI's "no chunking" branch) the chunk is searched alone since M4b-3 (the Owner's answers Q1 = B
// and Q2 = B of M4b, 2026-11-21; the assertions of those masks give `None`); a mask inside and a
// sequence shorter than the shortest searched one are searched. The bounds of the length are read
// where they agree, the length otherwise. The second block holds the cases of review F2 of M4b PR
// 2, where the manual's earlier words (an only mask, of 3 bases or more, within 4 bases of the end)
// were wrong both ways (`M/m4b/scripts3/rule_check4.py`: the manual's rule and a model of NCBI's
// chunks with the crash agree on every sorted set of 1 to 3 masks of sequences of 9 to 44 bases;
// NCBI BLAST+ 2.17.0 crashed on [0, 1), [1, 2), [0, 2) of 1,200 bases, [1200, 1202) and [1201,
// 1202) of 1,202, [0, 1) with [1199, 1200) of 1,200 and [0, 3) with [1200, 1201) of 1,201, and
// searched [0, 1) with [600, 700) and [300, 400) with [1200, 1202)).
#[test]
fn translated_hard_masks_are_refused_where_a_frame_restarts_at_its_start() {
    use std::cell::Cell;
    let scan = HardScan::Translated { min_length: 17 };
    let shape = |masks: &[(u32, u32)], length: usize| {
        HardMaskShape::of_scan(
            masks,
            None,
            &|| Some(length),
            scan,
            SearchedPlace::default(),
        )
        .unwrap()
    };
    let crash = |frame, start, end| Some(HardMaskShape::FrameStartMask { frame, start, end });
    assert_eq!(shape(&[(300, 390)], 1200), None);
    assert_eq!(shape(&[(300, 390), (900, 960)], 1200), None);
    assert_eq!(shape(&[(600, 700), (700, 800)], 1200), None);
    assert_eq!(shape(&[(1000, 1200)], 1200), None);
    assert_eq!(shape(&[(0, 1200)], 1200), None);
    assert_eq!(shape(&[(0, 100)], 1200), None);
    assert_eq!(shape(&[(0, 3)], 1200), None);
    assert_eq!(shape(&[(0, 2)], 1200), crash(1, 0, 2));
    assert_eq!(shape(&[(1000, 1198)], 1200), None);
    assert_eq!(shape(&[(1000, 1197)], 1200), None);
    assert_eq!(shape(&[(1000, 1196)], 1200), None);
    // A mask that runs to the end, or a head mask, after other masks: the chunk from the start
    // ends at the first mask (plus frames) or the last (minus frames).
    assert_eq!(shape(&[(300, 390), (900, 960)], 960), None);
    assert_eq!(shape(&[(0, 50), (600, 700)], 1200), None);
    // Empty masks before it: the chunk after them starts at 0 as well, the frame's second chunk
    // (an empty range at 0 before it: NCBI crashes, the mask between them named); in the minus
    // frames the empty ranges come last and the chunk from 0 is the first.
    assert_eq!(shape(&[(0, 0), (1000, 1200)], 1200), crash(1, 0, 0));
    assert_eq!(shape(&[(0, 0), (0, 100)], 1200), None);
    // A subject shorter than `min_length` bases is not searched.
    assert_eq!(shape(&[(0, 16)], 16), None);
    assert_eq!(shape(&[(0, 17)], 17), None);
    // Review F2: a mask of 1 or 2 bases at the start (frame -1); 3 bases within the first 6 that
    // leave the minus frames a letter after them; a mask at the end that leaves frame +1's last
    // codon; the minus frames one letter shorter than frame -1.
    assert_eq!(shape(&[(2, 3)], 1200), None);
    assert_eq!(shape(&[(1, 3)], 1200), None);
    assert_eq!(shape(&[(4, 7)], 1202), None);
    assert_eq!(shape(&[(4, 6)], 1200), None);
    assert_eq!(shape(&[(4, 6)], 1201), None);
    assert_eq!(shape(&[(1196, 1200)], 1200), None);
    // NCBI's crash: masks within the first two bases (the plus frames) or within the bases after
    // the last whole codon (the minus frames), alone or with a mask at the other end (frame +1
    // with [1199, 1200); frame -1 with [0, 3), which frame -1 also reads as its part from the
    // start ending before its end: the crash comes first); another mask in the frame searches it.
    assert_eq!(shape(&[(0, 1)], 1200), crash(1, 0, 1));
    assert_eq!(shape(&[(1, 2)], 1200), crash(1, 1, 2));
    assert_eq!(shape(&[(1200, 1202)], 1202), crash(-1, 1200, 1202));
    assert_eq!(shape(&[(1201, 1202)], 1202), crash(-1, 1201, 1202));
    assert_eq!(shape(&[(0, 1), (1199, 1200)], 1200), crash(1, 0, 1));
    assert_eq!(shape(&[(0, 3), (1200, 1201)], 1201), crash(-1, 1200, 1201));
    assert_eq!(shape(&[(0, 1), (600, 700)], 1200), None);
    assert_eq!(shape(&[(300, 400), (1200, 1202)], 1202), None);
    // The shortest searched subject: `min_length` (17 bases here) and one base less.
    assert_eq!(shape(&[(14, 17)], 17), None);
    assert_eq!(shape(&[(13, 16)], 16), None);
    // Masks out of order come first, as for the other programs.
    assert_eq!(
        shape(&[(600, 700), (300, 390)], 1200),
        Some(HardMaskShape::Unsorted {
            before: (600, 700),
            mask: (300, 390)
        })
    );
    // Bounds that give one shape at every length are taken without the length.
    let reads = Cell::new(0);
    let read = || {
        reads.set(reads.get() + 1);
        Some(1200)
    };
    let of = |masks: &[(u32, u32)], bounds| {
        HardMaskShape::of_scan(masks, Some(bounds), &read, scan, SearchedPlace::default())
    };
    assert_eq!(of(&[(300, 390)], (1197, 1200)), Ok(None));
    assert_eq!(reads.get(), 0);
    // Bounds whose lengths differed in shape before M4b-3 (at 1197 frame +2, of 398 residues, has
    // no letter after the mask, at 1200 every frame has): one shape now, the length is not read.
    assert_eq!(of(&[(1000, 1196)], (1197, 1200)), Ok(None));
    assert_eq!(reads.get(), 0);
    // A mask that ends past the shortest length: the length is read.
    assert_eq!(of(&[(1000, 1198)], (1197, 1200)), Ok(None));
    assert_eq!(reads.get(), 1);
    // A mask that ends past the sequence (not written by makeblastdb) is refused as for the
    // protein programs, before the frames.
    let past = Some(HardMaskShape::PastEnd {
        start: 300,
        end: 1250,
        length: 1200,
    });
    assert_eq!(of(&[(300, 1250)], (1197, 1200)), Ok(past));
    assert_eq!(reads.get(), 2);
    assert_eq!(
        shape(&[(0, 100), (300, 1250)], 1200),
        Some(HardMaskShape::PastEnd {
            start: 300,
            end: 1250,
            length: 1200,
        })
    );
    let text = HardMaskShape::FrameStartMask {
        frame: 1,
        start: 0,
        end: 1,
    }
    .text();
    assert!(
        text.starts_with(
            "has a mask, [0, 1), that leaves reading frame +1 an empty part at its start"
        ),
        "{text}"
    );
    assert!(!text.contains("/ 3") && !text.contains("chunk"), "{text}");
}

/// `HardMaskShape::of_scan` of the protein programs' scan (series NP, M4b): LOSAT's refusals that
/// NCBI BLAST+ cannot freeze (a range longer than its chunks; a mask past the sequence, review L1
/// of PR 1, which makeblastdb does not write), the ones that the rows db.pmask.*.rej.* also hold,
/// and masks that the search takes, chunks shorter than the word included (rows
/// db.pmask.*.pshort.*, db.pmask.*.pe_tiny.w0.*).
#[test]
fn protein_shapes_of_the_hard_masks() {
    let protein = HardScan::Protein {
        word: 3,
        window: 40,
        zero_words: false,
    };
    let of = |masks: &[(u32, u32)], length: usize, scan: HardScan| {
        HardMaskShape::of_scan(
            masks,
            Some((length, length)),
            &|| Some(length),
            scan,
            SearchedPlace::default(),
        )
        .unwrap()
    };
    // Masks within the sequence: searched, also where a chunk is shorter than the word (NCBI's
    // one-hit scan reads no word of it, its two-hit scan one word that extends nothing).
    assert_eq!(of(&[(10, 20), (50, 60)], 100, protein), None);
    assert_eq!(of(&[(10, 20), (22, 60)], 100, protein), None);
    assert_eq!(of(&[(10, 20), (20, 60)], 100, protein), None);
    assert_eq!(of(&[(1, 10)], 100, protein), None);
    // A mask that ends at the end after another chunk: the last chunk is not the first.
    assert_eq!(of(&[(10, 20), (50, 100)], 100, protein), None);
    // A mask that runs to the end after the chunk from the start (Q1 = B of M4b): searched over
    // that chunk, at every word; masked whole too.
    assert_eq!(of(&[(10, 100)], 100, protein), None);
    assert_eq!(of(&[(0, 100)], 100, protein), None);
    // The same chunk after an empty mask at 0 (a hand mask `from 0, to 4294967295`), as in
    // BLASTN: searched over that chunk (NCBI's "no chunking" branch after an allocated chunk
    // frees memory that it did not allocate; the near input's normal chunks give the result).
    assert_eq!(of(&[(0, 0), (10, 100)], 100, protein), None);
    assert_eq!(of(&[(0, 0), (0, 0), (0, 100)], 100, protein), None);
    // A mask that ends past the sequence, or starts past it: its own shape, not a long range.
    let past = of(&[(10, 20), (90, 120)], 100, protein);
    assert_eq!(
        past,
        Some(HardMaskShape::PastEnd {
            start: 90,
            end: 120,
            length: 100
        })
    );
    assert!(past
        .unwrap()
        .text()
        .starts_with("has a mask, [90, 120), ending past its 100 letters"));
    assert_eq!(
        of(&[(150, 160)], 100, protein),
        Some(HardMaskShape::PastEnd {
            start: 150,
            end: 160,
            length: 100
        })
    );
    // A hard range longer than NCBI's chunks.
    let long = crate::core::subject_split::MAX_DBSEQ_LEN + 100;
    let shape = of(&[(10, 20)], long, protein);
    assert_eq!(
        shape,
        Some(HardMaskShape::LongRange {
            left: 20,
            right: long as i32
        })
    );
    assert!(shape
        .unwrap()
        .text()
        .contains("longer than the 5,000,000 letters"));
    // At exactly MAX_DBSEQ_LEN letters the range is one chunk.
    assert_eq!(
        of(
            &[(10, 20)],
            crate::core::subject_split::MAX_DBSEQ_LEN + 20,
            protein
        ),
        None
    );
}

/// `HardMaskShape::of_scan` of the protein programs on masks that overlap or nest, or are out of
/// the order of their start (series NP, M4b-3b; the masks of probe 2 of M/m4b3b/NOTES_M4B3B.md, a
/// sequence of 240 letters, word 3): searched also where the walk of the chunks' diagonal offsets
/// leaves the sequence (each sequence on its own, PD-LOSAT-NCBI-DEFECTS exception 23: the search's
/// `DiagonalWalk::start`), `WordPastEnd` where a short chunk's one word read past the sentinel can
/// change a result; `PartsOutOfOrder` where a chunk that can hold
/// alignments starts before the end of an earlier one, `Inverted` for a mask that ends before it
/// starts; in the other programs `Unsorted` for every overlap, as before.
#[test]
fn protein_shapes_of_overlapping_masks() {
    let two_hit = HardScan::Protein {
        word: 3,
        window: 40,
        zero_words: false,
    };
    let one_hit = HardScan::Protein {
        word: 3,
        window: 0,
        zero_words: false,
    };
    // A query with the gap `-`, whose lookup table can hold a word with the residue 0 (review H2).
    let gap_query = HardScan::Protein {
        word: 3,
        window: 40,
        zero_words: true,
    };
    let of = |masks: &[(u32, u32)], length: usize, scan: HardScan| {
        HardMaskShape::of_scan(
            masks,
            Some((length, length)),
            &|| Some(length),
            scan,
            SearchedPlace::default(),
        )
        .unwrap()
    };
    // e3, nest and big: searched; in the one-hit search (no window) nest and big go back below the
    // sequence's start.
    assert_eq!(of(&[(40, 60), (55, 70)], 240, two_hit), None);
    assert_eq!(of(&[(40, 60), (55, 70)], 240, one_hit), None);
    assert_eq!(of(&[(40, 100), (50, 60)], 240, two_hit), None);
    assert_eq!(of(&[(40, 100), (50, 60)], 240, one_hit), None);
    assert_eq!(of(&[(40, 120), (60, 130)], 240, two_hit), None);
    assert_eq!(of(&[(40, 120), (60, 130)], 240, one_hit), None);
    // Masks with the same start.
    assert_eq!(of(&[(40, 100), (40, 60)], 240, two_hit), None);
    // lk: its entries reach 98 past its end in the two-hit search (137 in the one-hit search), and
    // with a query that can hit a word with the residue 0 its one word at 238, which covers the
    // sentinel, 20 below its start: searched (the walk's margins: `hard_chunk_walk` tests).
    assert_eq!(of(&[(40, 200), (60, 238)], 240, two_hit), None);
    assert_eq!(of(&[(40, 200), (60, 238)], 240, gap_query), None);
    assert_eq!(of(&[(40, 200), (60, 238)], 240, one_hit), None);
    // A short chunk at the last letter, whose one word reads past the sentinel, after a chunk that
    // went back within the window of its entries.
    // Only a query that can hit a word with the residue 0, the sentinel's, finds a hit there
    // (review H2).
    let past = of(&[(10, 60), (20, 99), (30, 40)], 100, gap_query);
    assert_eq!(past, Some(HardMaskShape::WordPastEnd { offset: 99 }));
    assert!(past
        .unwrap()
        .text()
        .starts_with("has a part between masks at 99 shorter than the word, whose one word"));
    assert_eq!(of(&[(10, 60), (20, 99), (30, 40)], 100, two_hit), None);
    assert_eq!(of(&[(10, 60), (20, 99), (30, 40)], 100, one_hit), None);
    // Out of the order of their start: searched where the chunks that can hold alignments do not
    // meet, as the seg masks of UniRef50_A0A835SUY3 (segmasker 2.17.0), whose third mask starts a
    // letter before the second (the chunk between them is inverted, -8 letters), and a mask
    // within the next one.
    let uniref = [(801, 819), (870, 877), (869, 916), (956, 994)];
    assert_eq!(of(&uniref, 1000, two_hit), None);
    assert_eq!(of(&uniref, 1000, one_hit), None);
    assert_eq!(of(&[(100, 130), (95, 140)], 240, two_hit), None);
    assert_eq!(of(&[(100, 130), (95, 140)], 240, one_hit), None);
    // A chunk shorter than the word holds no alignments ([100, 102), which ends past the last
    // chunk's offset 101); the one-hit search goes back below the sequence's start there.
    let short = [(50, 100), (102, 200), (60, 101)];
    assert_eq!(of(&short, 300, two_hit), None);
    assert_eq!(of(&short, 300, one_hit), None);
    // A chunk that can hold alignments before the end of an earlier one: the part after
    // [100, 130) starts before the end of the part before [250, 280) (the chunk between them is
    // inverted, so the masks name both); an empty mask at 0 after a mask, whose chunk then starts
    // at 0.
    let parts = |before, mask| Some(HardMaskShape::PartsOutOfOrder { before, mask });
    let shape = of(&[(250, 280), (100, 130)], 300, two_hit);
    assert_eq!(shape, parts((250, 280), (100, 130)));
    assert_eq!(
        shape.unwrap().text(),
        "has masks out of order ([100, 130) after [250, 280)): its part after [100, 130) starts \
         before the end of its part before [250, 280), and NCBI BLAST+ 2.17.0 joins the two parts' \
         alignments by a merge across that overlap, which LOSAT does not do yet"
    );
    assert_eq!(
        of(&[(250, 280), (100, 130)], 300, one_hit),
        parts((250, 280), (100, 130))
    );
    assert_eq!(
        of(&[(50, 60), (0, 0), (100, 300)], 300, two_hit),
        parts((50, 60), (0, 0))
    );
    // A mask that ends before it starts.
    let inverted = of(&[(50, 40)], 300, one_hit);
    assert_eq!(
        inverted,
        Some(HardMaskShape::Inverted { start: 50, end: 40 })
    );
    assert_eq!(
        inverted.unwrap().text(),
        "has a mask, [50, 40), that ends before it starts, which makeblastdb does not write"
    );
    // BLASTN and the translated programs refuse every overlap and every mask out of order.
    let unsorted = |before, mask| Some(HardMaskShape::Unsorted { before, mask });
    let blastn = HardScan::Nucleotide { scanned_word: 11 };
    let translated = HardScan::Translated { min_length: 11 };
    for scan in [blastn, translated] {
        assert_eq!(
            of(&[(40, 60), (55, 70)], 240, scan),
            unsorted((40, 60), (55, 70))
        );
        assert_eq!(
            of(&[(250, 280), (100, 130)], 300, scan),
            unsorted((250, 280), (100, 130))
        );
    }
}

/// The place of a sequence among the searched ones (`SearchedPlace`): the first searched
/// sequence's hits below its start meet the entries of NCBI's new table, a one word past the
/// sentinel with them; the margins below the start and past the end refuse no sequence at any
/// place (each sequence is searched on its own, PD-LOSAT-NCBI-DEFECTS exception 23).
#[test]
fn protein_shapes_at_every_searched_place() {
    let of = |masks: &[(u32, u32)], length: usize, window: i32, first: bool| {
        let scan = HardScan::Protein {
            word: 3,
            window,
            zero_words: true,
        };
        let place = SearchedPlace { first };
        HardMaskShape::of_scan(masks, Some((length, length)), &|| Some(length), scan, place)
            .unwrap()
    };
    for first in [false, true] {
        // nest, one-hit: 10 below the start; lk, two-hit with the residue 0: 20 below the start
        // and 98 past the end.
        assert_eq!(of(&[(40, 100), (50, 60)], 240, 0, first), None);
        assert_eq!(of(&[(40, 200), (60, 238)], 240, 40, first), None);
    }
}

/// A mask that runs to the end after the chunk from the start, with masks stored after it (review
/// L1, probe M/m4b3b/probe/fix1): refused in the protein programs, after a mask before it and
/// one inside it too; without the later masks searched as before (PD-LOSAT-NCBI-DEFECTS exception
/// 21).
#[test]
fn tail_masks_with_masks_stored_after_them() {
    let of = |masks: &[(u32, u32)], window: i32| {
        let scan = HardScan::Protein {
            word: 3,
            window,
            zero_words: false,
        };
        HardMaskShape::of_scan(
            masks,
            Some((240, 240)),
            &|| Some(240),
            scan,
            SearchedPlace::default(),
        )
        .unwrap()
    };
    let shape = HardMaskShape::TailMaskWithLaterMasks {
        start: 65,
        end: 240,
    };
    let refused = Some(shape);
    for window in [0, 40] {
        assert_eq!(of(&[(65, 240)], window), None);
        assert_eq!(of(&[(65, 240), (55, 64)], window), refused);
        assert_eq!(of(&[(65, 240), (100, 120)], window), refused);
        assert_eq!(of(&[(55, 64), (65, 240)], window), None);
    }
    assert_eq!(
        shape.text(),
        "has a mask, [65, 240), running to its end after the part from its start, and masks \
         stored after it, which NCBI BLAST+ 2.17.0 does not apply when it searches that part past \
         its end but applies to the letters of its alignments"
    );
}

/// Masks at random, apart in the order of their start, overlapping in that order, or out of it:
/// masks apart are never refused, and the shape is read at the sequence's exact length (bounds
/// that are not one length read it; `HardMaskShape::of_scan`). Masks out of the order of their
/// start are searched where their chunks allow it, and many are.
#[test]
fn protein_shapes_of_random_masks() {
    let mut state = 0x5851_f42d_4c95_7f2d_u64;
    let mut next = |below: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % below
    };
    let shape = |masks: &[(u32, u32)], length: usize, scan: HardScan| {
        HardMaskShape::of_scan(
            masks,
            Some((length, length)),
            &|| None,
            scan,
            SearchedPlace::default(),
        )
        .unwrap()
    };
    let (mut cases, mut refused) = (0, 0);
    let (mut out_of_order, mut searched, mut parts) = (0, 0, 0);
    for _ in 0..20_000 {
        let window = [0, 1, 40, 100][next(4) as usize];
        let scan = HardScan::Protein {
            word: 2 + next(5) as i32,
            window,
            zero_words: window == 100,
        };
        // 0: apart; 1: in the order of their start; 2: each start back by up to 29 letters.
        let order = next(3);
        let apart = order == 0;
        let mut masks = Vec::new();
        let mut start = 0u32;
        for _ in 0..1 + next(5) {
            start += next(60) as u32;
            if order == 2 {
                start = start.saturating_sub(next(30) as u32);
            }
            let end = start + next(50) as u32;
            masks.push((start, end));
            if apart {
                start = end;
            }
        }
        let last_end = masks.iter().map(|&(_, end)| end as usize).max().unwrap();
        let length = last_end + 1 + next(20) as usize;
        let high = length + 1 + next(12) as usize;
        let found = shape(&masks, length, scan);
        if apart {
            assert_eq!(found, None, "{masks:?} {length} {scan:?}");
        }
        // Bounds of more than one length read the length.
        assert_eq!(
            HardMaskShape::of_scan(
                &masks,
                Some((length, high)),
                &|| None,
                scan,
                SearchedPlace::default()
            ),
            Err(())
        );
        assert_eq!(
            HardMaskShape::of_scan(
                &masks,
                Some((length, high)),
                &|| Some(length),
                scan,
                SearchedPlace::default()
            ),
            Ok(found)
        );
        cases += 1;
        refused += usize::from(found.is_some());
        if masks.windows(2).any(|pair| pair[1].0 < pair[0].0) {
            out_of_order += 1;
            searched += usize::from(found.is_none());
            parts += usize::from(matches!(found, Some(HardMaskShape::PartsOutOfOrder { .. })));
        }
    }
    assert!(cases == 20_000 && refused > 50, "{cases} {refused}");
    assert!(
        out_of_order > 1_000 && searched > 500 && parts > 100,
        "{out_of_order} {searched} {parts}"
    );
}
