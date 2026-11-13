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
    db.read_name_records(&name).unwrap();
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
    let error = resolved.check_losat_limits(&db, "BLASTN").unwrap_err();
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

// The search that `-db_hard_mask` masks is refused where NCBI starts it (series NP, M4a-2b).
#[test]
fn a_hard_masked_search_is_refused() {
    let (a5, _) = open(&[dbmask("a5")]);
    let mask = ResolvedMask {
        kind: MaskKind::Hard,
        id: 11,
        key: None,
        shown: Vec::new(),
    };
    let error = mask.check_losat_limits(&a5, "BLASTN").unwrap_err();
    assert_eq!(
        error.to_string(),
        "-db_hard_mask (hard masking of the database sequences) is not supported by LOSAT's BLASTN"
    );
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
