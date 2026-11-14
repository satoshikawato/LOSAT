//! Tests of the database masks' resolver (`db_mask.rs`), each case a probe of the M4a-2 survey
//! (`M/m4a2/SURVEY.md` §1, NCBI BLAST+ 2.17.0) on the databases of `tests/db_search/dbmask`, which
//! NCBI's makeblastdb wrote from the same inputs (README there): `a5` (dust 11, dust
//! `-level 40 -window 32` 12, windowmasker 30), `b5` (the second dust alone, its id 11), `k5`
//! (`-mask_id mydust,mywm`: 100, 101), `nm5` (no masks); `tests/seqdb/db/mask_nucl_colon_id`
//! (`-mask_id c:id`, stored with 0x01 for the colon).

use std::path::{Path, PathBuf};

use super::*;
use crate::seqdb::SeqType;

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
    (native.exit, native.message.clone())
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
                format!("BLAST options error: Masking algorithm ID {id} is not supported in nucleotide 'db/a5' BLAST database\n")
            )
        );
    }
    let (nm5, _) = open(&[dbmask("nm5")]);
    assert_eq!(
        native(soft("11").validate_at_open(&nm5, "db/nm5", true, "BLASTN")),
        (
            1,
            "BLAST options error: Masking algorithm ID 11 is not supported in nucleotide 'db/nm5' BLAST database\n"
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
            (2, format!("BLAST Database error: Filtering algorithm ID {id} is not supported.\n{table}\n"))
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
                "BLAST Database error: Filtering algorithm ID 11 is not supported.\n\n".to_string()
            )
        );
    }
    let (db, name) = open(&[dbmask("a5"), dbmask("b5")]);
    let (resolved, _) = resolve(&soft("12"), &db, &name);
    assert_eq!(
        native(resolved).1,
        "BLAST Database error: Filtering algorithm ID 12 is not supported.\n\n\
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
    let error = resolved.check_losat_limits(&db, "BLASTN", 11).unwrap_err();
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
// at `HardMaskShape::of`)
// A search with `-db_hard_mask` is refused, before anything is written, at the first searched
// sequence whose masks NCBI's hard chunks do not handle (series NP, M4a-2b decisions H3-H5): a5's
// s4 (OID 3) is masked whole ([0, 300): the "no chunking" branch scans the full range of a chunk
// of no letters); edge1's algorithm 11 holds s5's masks out of order, 12 overlapping ones; 13
// (adjacent masks: a hard range of no letters), 14 (masks at both ends) and 15 (one-letter masks)
// are searched; edge2's 14 (a tail mask) and 16 (s4 masked whole, OID 3 before s5) are refused.
// `-db_soft_mask` keeps every shape.
#[test]
fn a_hard_masked_search_refuses_the_masks_that_ncbi_s_chunks_do_not_handle() {
    let mask = |kind, id| ResolvedMask {
        kind,
        id,
        key: None,
        shown: Vec::new(),
    };
    let refusal = |db: &SeqDb, id| {
        mask(MaskKind::Hard, id)
            .check_losat_limits(db, "BLASTN", 11)
            .map_err(|error| error.to_string())
    };
    let (a5, _) = open(&[dbmask("a5")]);
    let error = refusal(&a5, 11).unwrap_err();
    assert!(
        error.starts_with("-db_hard_mask: the sequence of OID 3 (volume ")
            && error.contains(
                "a5) has its only mask, [0, 300), running to its end, for which NCBI BLAST+ 2.17.0 \
             scans the masked letters past the chunk before the mask"
            )
            && error.ends_with("LOSAT's BLASTN does not search such a database with -db_hard_mask"),
        "{error}"
    );
    mask(MaskKind::Soft, 11)
        .check_losat_limits(&a5, "BLASTN", 11)
        .unwrap();
    let (edge1, _) = open(&[dbmask("edge1")]);
    let error = refusal(&edge1, 11).unwrap_err();
    assert!(
        error.contains("OID 4 ")
            && error.contains("has masks out of order ([100, 200) after [600, 700))"),
        "{error}"
    );
    let error = refusal(&edge1, 12).unwrap_err();
    assert!(error.contains("([200, 400) after [100, 300))"), "{error}");
    for id in [13, 14, 15] {
        refusal(&edge1, id).unwrap();
    }
    let (edge2, _) = open(&[dbmask("edge2")]);
    let error = refusal(&edge2, 14).unwrap_err();
    assert!(
        error.contains("OID 4 ") && error.contains("[1000, 1200), running to its end"),
        "{error}"
    );
    let error = refusal(&edge2, 16).unwrap_err();
    assert!(
        error.contains("OID 3 ") && error.contains("[0, 300), running"),
        "{error}"
    );
    for id in 11..=16 {
        mask(MaskKind::Soft, id)
            .check_losat_limits(&edge2, "BLASTN", 11)
            .unwrap();
    }
}

// The shapes of `HardMaskShape::of` (blast_engine.c:241-268, 579-584, 1318-1331; na_ungapped.c:
// 1647-1668): the first mask's end at the subject's end with its start at or below MAX_DBSEQ_LEN, in
// a subject of a word or more and at least 5 letters (empty masks at the end may follow), masks
// that go back, and a last mask shorter than its end's residual before the subject's end. The
// length is read only for the first and the last.
#[test]
fn the_hard_mask_shapes_follow_ncbi_s_chunks() {
    use std::cell::Cell;
    let shape = |masks: &[(u32, u32)], length: usize| {
        let reads = Cell::new(0);
        let read = || {
            reads.set(reads.get() + 1);
            Some(length)
        };
        (HardMaskShape::of(masks, &read, 11).unwrap(), reads.get())
    };
    let tail = |start, end| Some(HardMaskShape::TailMask { start, end });
    assert_eq!(shape(&[(100, 300)], 300), (tail(100, 300), 1));
    assert_eq!(shape(&[(0, 300)], 300), (tail(0, 300), 1));
    assert_eq!(shape(&[(100, 300), (300, 300)], 300), (tail(100, 300), 1));
    assert_eq!(
        shape(&[(5_000_000, 6_000_000)], 6_000_000),
        (tail(5_000_000, 6_000_000), 1)
    );
    // Above MAX_DBSEQ_LEN the first chunk is a full chunk: NCBI chunks it as any other.
    assert_eq!(shape(&[(5_000_001, 6_000_000)], 6_000_000), (None, 0));
    // An empty mask at the end, a mask before the end, two masks.
    assert_eq!(shape(&[(300, 300)], 300), (None, 0));
    assert_eq!(shape(&[(100, 299)], 300), (None, 1));
    assert_eq!(shape(&[(0, 100), (200, 300)], 300), (None, 0));
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
    assert_eq!(shape(&[(100, 200), (200, 300)], 1200), (None, 0));
    let short = |start, end| Some(HardMaskShape::ShortLastMask { start, end });
    // 1002 % 4 = 2 letters before the end in the last chunk, a mask of 1.
    assert_eq!(shape(&[(1001, 1002)], 1200), (short(1001, 1002), 1));
    assert_eq!(shape(&[(1000, 1002)], 1200), (None, 1));
    assert_eq!(shape(&[(1003, 1003)], 1200), (short(1003, 1003), 1));
    // The last mask at the end: the last hard range is empty and the loop ends before it.
    assert_eq!(shape(&[(10, 20), (1199, 1199)], 1199), (None, 1));
    // A middle mask shorter than its residual merges with an overlap that NCBI reads in bounds.
    assert_eq!(shape(&[(301, 302), (500, 503)], 1200), (None, 0));
    // A subject that the scan reads no word of: shorter than the word, or than 5 letters (the
    // engine skips it).
    assert_eq!(shape(&[(0, 10)], 10), (None, 1));
    assert_eq!(shape(&[(5, 11)], 11), (tail(5, 11), 1));
    assert_eq!(HardMaskShape::of(&[(0, 4)], &|| Some(4), 4), Ok(None));
    assert_eq!(HardMaskShape::of(&[(0, 5)], &|| Some(5), 4), Ok(tail(0, 5)));
    // The length that cannot be read ends the check.
    assert_eq!(HardMaskShape::of(&[(100, 300)], &|| None, 11), Err(()));
}

// blast_app_util.cpp:644-651 and local_db_adapter.cpp:54-58 (quoted at `strategy_params`): after
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
            .check_losat_limits(db, "BLASTN", 11)
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
