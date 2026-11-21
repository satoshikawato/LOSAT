use super::*;
use crate::ncbi_exception::NCBI_BUILD_ROOT;
use crate::utils::winmask::{FileCounts, FormatSite};

fn counts_path(name: &str) -> String {
    format!("{}/tests/winmask/counts/{name}", env!("CARGO_MANIFEST_DIR"))
}

// NCBI BLAST+ 2.17.0 `windowmasker -ustat mg1655_u9.obinary -outfmt interval` of the
// first record of tests/fixtures/blastn_regression/inputs/wm_query.fa (132 - 146, ...),
// and of its reverse complement (390 - 402, ...), which are the masks of blastn
// `-window_masker_db mg1655_u9.obinary` and of `-strand minus` (its archive's LCaseMask).
#[test]
fn the_masks_are_windowmaskers_of_the_strand_read() {
    let fasta = std::fs::read_to_string(format!(
        "{}/tests/fixtures/blastn_regression/inputs/wm_query.fa",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap();
    let first: String = fasta.split('>').nth(1).unwrap().lines().skip(1).collect();
    let masker = open(&counts_path("mg1655_u9.obinary"), &FileCounts).unwrap();
    let plus = masker.masks(first.as_bytes(), false).unwrap();
    assert_eq!(
        plus[..3],
        [
            MaskedInterval::new(132, 147),
            MaskedInterval::new(414, 429),
            MaskedInterval::new(489, 502)
        ]
    );
    // The lower-case letters of the record are letters like the others.
    assert_eq!(
        masker
            .masks(first.to_ascii_uppercase().as_bytes(), false)
            .unwrap(),
        plus
    );
    let minus = masker.masks(first.as_bytes(), true).unwrap();
    assert_eq!(
        minus[..2],
        [MaskedInterval::new(390, 403), MaskedInterval::new(559, 572)]
    );
}

#[test]
fn the_checks_of_the_file_are_ncbis() {
    let missing = open("nosuch.obinary", &FileCounts).unwrap_err();
    assert_eq!(missing.exit, 255);
    assert_eq!(
        missing.message,
        // NCBI's text in LOSAT's words (`crate::wording::losat_text`).
        crate::wording::losat_text(&format!("Error: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/algo/winmask/seq_masker_istat_factory.cpp\", line 171: Error: (Exception::open failed) ncbi::CSeqMaskerIstatFactory::DiscoverStatType() - could not open nosuch.obinary\n\n")).as_bytes()
    );
    // A directory opens and reads nothing on Unix and on WASI: the unknown format. Windows
    // cannot open a directory as a file, so the open fails first there: the "could not open"
    // exception.
    let directory = open(env!("CARGO_MANIFEST_DIR"), &FileCounts).unwrap_err();
    if cfg!(any(unix, target_os = "wasi")) {
        assert_eq!(directory.exit, 1);
        assert_eq!(
            directory.message,
            b"BLAST query/options error: Only optimized binary windowmasker stat files are supported\nPlease refer to the LOSAT user manual.\n"
        );
    } else {
        assert_eq!(directory.exit, 255);
        assert!(
            directory
                .message
                .ends_with(format!("could not open {}\n\n", env!("CARGO_MANIFEST_DIR")).as_bytes()),
            "{}",
            String::from_utf8_lossy(&directory.message)
        );
    }
    let binary = open(&counts_path("dm2L_u9.binary"), &FileCounts).unwrap();
    assert!(binary.build_error().is_none());
    assert_eq!(binary.name, counts_path("dm2L_u9.binary"));
}

// NCBI BLAST+ 2.17.0 blastn `-window_masker_taxid` (line NP M/STATE.md, "NCBI 2.17.0 の実機"):
// the file under the directory is used whatever its format is called; a path that is no file
// is the empty name, which fails each batch with `could not open ` and the factory's error;
// ascii counts are read (line NP M3d), their errors are NCBI's and their warnings come after
// the registry's posts at each load.
#[test]
fn the_file_of_a_taxid_is_the_masker_or_the_empty_name() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-blastn-window-masker-taxid-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    for (taxid, file) in [("562", "mg1655_u9.obinary"), ("7227", "dm2L_u9.binary")] {
        std::fs::create_dir_all(dir.join(taxid)).unwrap();
        std::fs::copy(counts_path(file), dir.join(taxid).join("wmasker.obinary")).unwrap();
    }
    std::fs::create_dir_all(dir.join("10")).unwrap();
    std::fs::write(
        dir.join("10/wmasker.obinary"),
        b"##windowmasker-statistics-format-version:ascii 1.1.0\n9\nAAAAAAAAA 5\n",
    )
    .unwrap();
    let root = WindowMaskerRegistry {
        path: Some(dir.to_str().unwrap().to_string()),
        ..WindowMaskerRegistry::default()
    };
    let root = &root;
    let found = open_taxid(562, root, &FileCounts).unwrap();
    assert_eq!(found.taxid, Some(562));
    assert!(found.name.ends_with("wmasker.obinary"), "{}", found.name);
    assert!(found.build_error().is_none());
    assert_eq!(found.filter_option(), FilterOption::Taxid(562));
    assert!(found.load_posts(true).is_empty() && found.load_posts(false).is_empty());
    // What each load of the registry writes, or only the first.
    let mut posted = Posts::new();
    posted.post(crate::report::diag_log::CRITICAL, "posted\n");
    for (once, later) in [(false, posted.clone()), (true, Posts::new())] {
        let registry = WindowMaskerRegistry {
            load_posts: posted.clone(),
            posts_once: once,
            ..root.clone()
        };
        let masker = open_taxid(562, &registry, &FileCounts).unwrap();
        assert_eq!(masker.load_posts(true), &posted);
        assert_eq!(masker.load_posts(false), &later);
    }
    // A binary counts file under the optimized binary's name.
    assert!(open_taxid(7227, root, &FileCounts)
        .unwrap()
        .build_error()
        .is_none());
    let missing = open_taxid(9606, root, &FileCounts).unwrap();
    assert_eq!(missing.name, "");
    let none = &Posts::new();
    let (warnings, error) = batch_error(missing.build_error().unwrap(), 2, none, none).unwrap();
    assert!(warnings.is_empty());
    // Each query's read of the file comes after a load of the registry.
    let (first, later) = (Posts::of_lines(b"first\n".to_vec()), posted);
    let (posts, _) = batch_error(missing.build_error().unwrap(), 3, &first, &later).unwrap();
    assert_eq!(posts.text(), b"first\nposted\nposted\n");
    let one = format!("Warning: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/algo/winmask/seq_masker_istat_factory.cpp\", line 171: Error: (Exception::open failed) ncbi::CSeqMaskerIstatFactory::DiscoverStatType() - could not open \n    T0 \"{NCBI_BUILD_ROOT}/src/algo/winmask/seq_masker_istat_factory.cpp\", line 358: Error: (Exception::creation failure) ncbi::CSeqMaskerIstatFactory::create() - could not create a unit counts container\n ");
    assert_eq!(error.exit, 3);
    // NCBI's text in LOSAT's words (`crate::wording::losat_text`).
    assert_eq!(
        error.message,
        crate::wording::losat_text(&format!("BLAST engine error: {one}{one}\n")).as_bytes()
    );
    // "AAAAAAAAA" is a unit beyond `Uint4` (base 16).
    let ascii = open_taxid(10, root, &FileCounts).unwrap();
    let report = ascii.build_error().unwrap().create_report().unwrap();
    let report = String::from_utf8(report).unwrap();
    assert!(
        report.contains("Cannot convert string 'AAAAAAAAA' to unsigned int, overflow"),
        "{report}"
    );
    // A second, higher t_low: its warning at every load, after the registry's posts.
    std::fs::create_dir_all(dir.join("15")).unwrap();
    std::fs::write(
        dir.join("15/wmasker.obinary"),
        b"9\n2 76\n>t_low 74\n>t_extend 147\n>t_threshold 173\n>t_high 210\n>t_low 80\n",
    )
    .unwrap();
    let warning = Posts::of_lines(b"Warning: [blastn] Requested value of t_low (74) is less than the one stored with the N-mer counts (80).The value 80 will be used.\n".to_vec());
    let mut both = Posts::of_lines(b"posted\n".to_vec());
    both.append(&warning);
    for (once, later) in [(false, &both), (true, &warning)] {
        let registry = WindowMaskerRegistry {
            load_posts: Posts::of_lines(b"posted\n".to_vec()),
            posts_once: once,
            ..root.clone()
        };
        let masker = open_taxid(15, &registry, &FileCounts).unwrap();
        assert!(masker.build_error().is_none());
        assert_eq!(masker.load_posts(true), &both);
        assert_eq!(masker.load_posts(false), later);
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// NCBI BLAST+ 2.17.0 (the official Linux x64 release) blastn with a counts file whose unit
// size is 17, and with one whose cache bit array is cut short (line NP M/STATE.md, "M2 の
// 実機確認"; the bytes are those of tests/blastn_window_masker.rs's fixtures).
#[test]
fn a_counts_file_that_fails_stops_the_batch_with_one_message_per_query() {
    let none = &Posts::new();
    let factory = format!("    T0 \"{NCBI_BUILD_ROOT}/src/algo/winmask/seq_masker_istat_factory.cpp\", line 358: Error: (Exception::creation failure) ncbi::CSeqMaskerIstatFactory::create() - could not create a unit counts container\n");
    let unit = format!("Warning: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/algo/winmask/seq_masker_istat_obinary.cpp\", line 114: Error: (Exception::bad parameter) ncbi::CSeqMaskerIstatOBinary::CSeqMaskerIstatOBinary() - unit size must be in [1,16]\n{factory} ");
    let (warnings, error) =
        batch_error(&WinMaskError::format(FormatSite::UnitSize), 2, none, none).unwrap();
    assert!(warnings.is_empty());
    assert_eq!(error.exit, 3);
    // NCBI's text in LOSAT's words (`crate::wording::losat_text`).
    assert_eq!(
        error.message,
        crate::wording::losat_text(&format!("BLAST engine error: {unit}{unit}\n")).as_bytes()
    );
    let short = WinMaskError::Format {
        site: FormatSite::HashTable,
        warnings: vec!["file read failed: bit array optimizations are not used.".into()],
    };
    let (warnings, _) = batch_error(&short, 2, none, none).unwrap();
    assert_eq!(
        warnings.text(),
        b"Warning: [blastn] file read failed: bit array optimizations are not used.\nWarning: [blastn] file read failed: bit array optimizations are not used.\n"
    );
    let (_, bad_index) = batch_error(
        &WinMaskError::BadIndex {
            key: 8764,
            value: 14779,
        },
        1,
        none,
        none,
    )
    .unwrap();
    assert_eq!(
        bad_index.message,
        b"BLAST engine error: Warning: bad index at key 8764 : 14779 \n"
    );
    assert!(batch_error(&WinMaskError::EmptyCacheBitArray { unit: 2 }, 1, none, none).is_none());
}
