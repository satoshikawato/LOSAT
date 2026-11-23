//! Unit tests of makembindex's command line (the index bytes are compared with NCBI's in
//! `tests/makembindex_fixtures.rs`).

use std::ffi::OsString;

use super::*;

/// `CArgAllow_Integers::GetUsage` of makembindex's constraints, as its `-help` prints them.
#[test]
fn integer_constraints_print_ncbi_usage() {
    let usage = |name: &str| crate::ncbi_app::args::constraint_usage("makembindex", name).unwrap();
    assert_eq!(usage("nmer"), "8..15");
    for name in ["ws_hint", "volsize", "stride"] {
        assert_eq!(usage(name), "greater or equal to 1");
        assert!(HELP.contains(&format!(" -{name} <Integer, greater or equal to 1>\n")));
    }
    assert!(HELP.contains(" -nmer <Integer, 8..15>\n"));
}

/// A run of `makembindex` with `words`, FASTA `input` on its standard input, in a directory of its
/// own (the old style writes `-output` relative to the current directory, so the words name it
/// whole): stdout, stderr, exit status and the files of the directory.
fn run_in_dir(tag: &str, words: &[&str], input: &[u8]) -> (String, String, i32, Vec<String>) {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-makembindex-{tag}-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let output = dir.join("o").to_string_lossy().into_owned();
    let mut words: Vec<OsString> = words.iter().map(OsString::from).collect();
    words.extend(["-output".into(), output.clone().into()]);
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let status = run_with(&words, &mut &input[..], &mut out, &mut err);
    let mut files: Vec<String> = std::fs::read_dir(&dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    files.sort();
    let err = String::from_utf8(err).unwrap().replace(&output, "@O@");
    let _ = std::fs::remove_dir_all(&dir);
    (String::from_utf8(out).unwrap(), err, status, files)
}

const FASTA: &[u8] = b">a\nACGTTGCAAGGCTTAACCGGTTAACCGGATCGACGTTGCAAGGCTTAACCGGTTAACCGGATCG\n";

/// `-legacy false` with a stride above `ws_hint - nmer + 1` (a zero `mult_`): NCBI's makembindex
/// does not finish the volume (its `Save` does not end); LOSAT refuses it when the volume has a
/// sequence, and writes no file. Without a sequence NCBI's run ends, and so does LOSAT's.
#[test]
fn a_stride_that_ncbi_does_not_finish_is_refused() {
    let words = [
        "-old_style_index",
        "true",
        "-legacy",
        "false",
        "-stride",
        "20",
    ];
    let (out, err, status, files) = run_in_dir("mult0", &words, FASTA);
    assert_eq!((out.as_str(), status, files.len()), ("", 1, 0));
    assert_eq!(
        err,
        "creating @O@.00.idx...\nError: -stride 20 above -ws_hint 28 minus -nmer 12 plus 1 \
         (NCBI's makembindex does not finish such an index volume) is not supported by LOSAT's \
         MAKEMBINDEX\n"
    );
    let (_, err, status, files) = run_in_dir("mult0empty", &words, b"");
    assert_eq!(
        (err.as_str(), status, files.len()),
        ("creating @O@.00.idx...removed (empty)\n", 0, 0)
    );
    // One less is NCBI's lowest multiple, 1: the volume is written.
    let words = [
        "-old_style_index",
        "true",
        "-legacy",
        "false",
        "-stride",
        "17",
    ];
    let (_, _, status, files) = run_in_dir("mult1", &words, FASTA);
    assert_eq!((status, files), (0, vec!["o.00.idx".to_string()]));
}

/// A stride of 32768 or more, where NCBI's `GetMinOffset` shifts an `int` by 32 bits or more
/// (undefined in C++): refused when the volume has a sequence.
#[test]
fn a_stride_past_ncbi_s_shift_is_refused() {
    let words = [
        "-old_style_index",
        "true",
        "-legacy",
        "false",
        "-stride",
        "32768",
        "-ws_hint",
        "40000",
    ];
    let (_, err, status, files) = run_in_dir("stride32768", &words, FASTA);
    assert_eq!((status, files.len()), (1, 0));
    assert!(
        err.ends_with(
            "Error: -stride 32768 (32768 or more, where NCBI's index has no defined result: \
             its offsets come from shifting an int past its width) is not supported by LOSAT's \
             MAKEMBINDEX\n"
        ),
        "{err}"
    );
    let words = [
        "-old_style_index",
        "true",
        "-legacy",
        "false",
        "-stride",
        "32767",
        "-ws_hint",
        "40000",
    ];
    let (_, _, status, files) = run_in_dir("stride32767", &words, FASTA);
    assert_eq!((status, files.len()), (0, 1));
}

/// `-old_style_index true -iformat blastdb` over an alias file whose `OID_MASK_TYPE` leaves a
/// sequence without a definition line (`tests/db_search/oidmask/ntall`, OID 1; line RG D1b-2):
/// NCBI's database reader gives no Bioseq for it and makembindex stops on the null pointer with a
/// report of the run's memory addresses (probe mk_old_ntall, exit status 3); LOSAT refuses it and
/// leaves no file. The volumes alone (`omn`) are indexed.
#[test]
fn an_old_style_index_over_a_sequence_without_a_line_is_refused() {
    let db = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/oidmask");
    let words = |name: &str| {
        let input = db.join(name).to_string_lossy().into_owned();
        let words = [
            "-old_style_index",
            "true",
            "-iformat",
            "blastdb",
            "-nmer",
            "8",
            "-input",
        ];
        let mut words: Vec<String> = words.iter().map(|word| word.to_string()).collect();
        words.push(input);
        words
    };
    let all = words("ntall");
    let all: Vec<&str> = all.iter().map(String::as_str).collect();
    let (out, err, status, files) = run_in_dir("oidmask-ntall", &all, b"");
    assert_eq!((out.as_str(), status, files.len()), ("", 1, 0));
    assert_eq!(
        err,
        "creating @O@.00.idx...\nError: sequence 1 of the database, which has no letters or no \
         definition line left by the alias files' OID_MASK_TYPE (NCBI's makembindex reads the null \
         sequence its database reader gives for it and stops with a report of the run's memory \
         addresses), is not supported by LOSAT's MAKEMBINDEX\n"
    );
    let volumes = words("omn");
    let volumes: Vec<&str> = volumes.iter().map(String::as_str).collect();
    let (_, _, status, files) = run_in_dir("oidmask-omn", &volumes, b"");
    assert_eq!((status, files), (0, vec!["o.00.idx".to_string()]));
}

/// The database `tests/db_search/t5/nhd_utf8` has a title stored as a UTF8String for a
/// VisibleString (line RG item T5d): NCBI 2.17.0's makembindex writes its warning 803.10, without
/// a prefix, at the read of that header, inside the line of the volume (probe of T5d merge 18,
/// 2026-11-23, Linux x86-64; the old style keeps the index out of the fixtures' directory).
#[test]
fn a_utf8_string_of_a_header_is_warned_in_the_volume_s_line() {
    let input =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/t5/nhd_utf8");
    let input = input.to_string_lossy().into_owned();
    let words = [
        "-old_style_index",
        "true",
        "-iformat",
        "blastdb",
        "-input",
        &input,
    ];
    let (out, err, status, files) = run_in_dir("utf8", &words, b"");
    assert_eq!((out.as_str(), status), ("", 0));
    assert_eq!(files, vec!["o.00.idx".to_string()]);
    assert_eq!(
        err,
        "creating @O@.00.idx...Warning: (803.10) CObjectIStreamAsnBinary: UTF8String data for \
         VisibleString member Blast-def-line-set.E.title, ASN.1 specification may need an \
         update\ndone\ncreating @O@.01.idx...removed (empty)\n"
    );
}
