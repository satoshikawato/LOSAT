//! Tests of the mask-file reader.
//!
//! The fixtures are in `tests/seqdb/mask_parse` (`make_mask_parse_fixtures.py`): files written by
//! NCBI's dustmasker, segmasker, windowmasker and convert2blastmask in both encodings, hand-written
//! text files and damaged binary files. `expected.tsv` records what makeblastdb 2.17.0 does with
//! each of them (exit status, standard error, the masks it stores); the tests below check that the
//! reader accepts what makeblastdb reads and refuses what its serial reader refuses.
//! `dump_dir` (ignored) prints one line per file of a directory, for comparing the reader with
//! makeblastdb on files made by a fuzzer.

use super::*;
use crate::writedb::WriteDbError;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn seqdb_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb")
}

fn parse_dir() -> PathBuf {
    seqdb_dir().join("mask_parse")
}

fn read_file(path: &Path) -> Vec<u8> {
    std::fs::read(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
}

/// What makeblastdb did with a file (`expected.tsv`).
struct Run {
    dbtype: String,
    fasta: String,
    rc: i32,
    stderr: String,
    no_db: bool,
    /// The options of the algorithm as `blastdbcmd -info` prints them; empty without a database.
    options: String,
    /// `MASK` rows: algorithm id -> oid -> ranges.
    masks: BTreeMap<u32, BTreeMap<u32, String>>,
}

fn unescape(field: &str) -> String {
    if field == "\\e" {
        return String::new();
    }
    let bytes = field.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] != b'\\' {
            out.push(bytes[i]);
            i += 1;
            continue;
        }
        match bytes[i + 1] {
            b'\\' => out.push(b'\\'),
            b't' => out.push(b'\t'),
            b'n' => out.push(b'\n'),
            b'r' => out.push(b'\r'),
            b's' => out.push(b' '),
            b'e' => {}
            b'x' => {
                out.push(u8::from_str_radix(&field[i + 2..i + 4], 16).unwrap());
                i += 2;
            }
            other => panic!("unknown escape \\{}", other as char),
        }
        i += 2;
    }
    // The standard error of NCBI's programs is bytes; the tests only compare it with ASCII.
    String::from_utf8_lossy(&out).into_owned()
}

fn read_expected() -> BTreeMap<String, Run> {
    let text = String::from_utf8(read_file(&parse_dir().join("expected.tsv"))).unwrap();
    let mut runs: BTreeMap<String, Run> = BTreeMap::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "RUN" => {
                runs.insert(
                    f[1].to_string(),
                    Run {
                        dbtype: f[2].to_string(),
                        fasta: f[3].to_string(),
                        rc: f[4].parse().unwrap(),
                        stderr: unescape(f[5]),
                        options: unescape(f[6]),
                        no_db: false,
                        masks: BTreeMap::new(),
                    },
                );
            }
            "NODB" => runs.get_mut(f[1]).unwrap().no_db = true,
            "MASK" => {
                runs.get_mut(f[1])
                    .unwrap()
                    .masks
                    .entry(f[2].parse().unwrap())
                    .or_default()
                    .insert(f[3].parse().unwrap(), unescape(f[4]));
            }
            other => panic!("unknown row {other}"),
        }
    }
    runs
}

/// `x_ProcessMaskData` for one file: the objects, each with the continuations that its `more`
/// announces.
fn read_all(reader: &mut MaskFileReader) -> Result<Vec<(MaskInfo, Vec<MaskList>)>, MaskFileError> {
    let mut all = Vec::new();
    while let Some(info) = reader.next_info()? {
        let mut more = info.masks.more;
        let mut lists = Vec::new();
        while more {
            let list = reader.next_list()?;
            more = list.more;
            lists.push(list);
        }
        all.push((info, lists));
    }
    Ok(all)
}

fn open_all(path: &Path) -> Result<Vec<(MaskInfo, Vec<MaskList>)>, MaskFileError> {
    read_all(&mut MaskFileReader::open(path)?)
}

fn infos(path: &Path) -> Vec<MaskInfo> {
    open_all(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .into_iter()
        .map(|(info, _)| info)
        .collect()
}

fn local(n: i32) -> Option<MaskSeqId> {
    Some(MaskSeqId::LocalId(n))
}

fn interval(from: u32, to: u32) -> MaskInterval {
    MaskInterval {
        from,
        to,
        strand: None,
    }
}

/// One line per file of `$MASK_DIFF_DIR` (sorted by name): what the reader makes of it. This is
/// the LOSAT side of the check against NCBI's makeblastdb that `make_mask_parse_fixtures.py`
/// and the notes of line NP/M4a describe; it is not a test.
#[test]
#[ignore = "helper for the comparison with makeblastdb"]
fn dump_dir() {
    let dir = PathBuf::from(std::env::var("MASK_DIFF_DIR").unwrap());
    let mut names: Vec<_> = std::fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().into_string().unwrap())
        .collect();
    names.sort();
    for name in names {
        println!("\n@@{name}\t{}", describe(&dir.join(&name)));
    }
}

fn describe(path: &std::path::Path) -> String {
    let mut reader = match MaskFileReader::open(path) {
        Ok(r) => r,
        Err(MaskFileError::Ncbi(e)) => return format!("UNKNOWN\t{e:?}"),
        Err(e) => return format!("OPENERR\t{e}"),
    };
    let mut out = String::new();
    let mut first = true;
    loop {
        match reader.next_info() {
            Ok(None) => return format!("OK\t{out}"),
            Ok(Some(info)) => {
                // makeblastdb registers the algorithm of the first object (it refuses
                // programs other than dust, seg, windowmasker, repeat, other) and stops at a
                // mask without a Seq-id before it reads on.
                let program_bad = first && ![10, 20, 30, 40, 100].contains(&info.algo_program);
                first = false;
                out += &format!(
                    "[a{} p{} o{:?} {}",
                    info.algo_id,
                    info.algo_program,
                    String::from_utf8_lossy(&info.algo_options),
                    dump_list(&info.masks)
                );
                let mut bad = program_bad || info.masks.masks.iter().any(|m| m.id.is_none());
                let mut more = info.masks.more;
                while more && !bad {
                    match reader.next_list() {
                        Ok(list) => {
                            out += &format!(" +{}", dump_list(&list));
                            bad = list.masks.iter().any(|m| m.id.is_none());
                            more = list.more;
                        }
                        Err(MaskFileError::EofInContinuation { .. }) => {
                            return format!("EOFCONT\t{out}");
                        }
                        Err(e) => return kind(&e),
                    }
                }
                if bad {
                    return format!("SEMERR\t{out}");
                }
                out += "]";
            }
            Err(e) => return kind(&e),
        }
    }
}

fn kind(e: &MaskFileError) -> String {
    match e {
        MaskFileError::Malformed { message, .. } => format!("MALFORMED\t{message}"),
        MaskFileError::Unsupported(m) => format!("UNSUPPORTED\t{m}"),
        other => format!("OTHER\t{other}"),
    }
}

fn dump_list(list: &MaskList) -> String {
    let mut out = format!("m{}", list.more as u8);
    for loc in &list.masks {
        out += " {";
        match &loc.id {
            Some(MaskSeqId::LocalId(n)) => out += &format!("L{n}"),
            Some(MaskSeqId::Other(b)) => out += &format!("O{}", b.len()),
            None => out += "-",
        }
        for iv in &loc.intervals {
            out += &format!(" {}-{}", iv.from, iv.to);
            if let Some(s) = iv.strand {
                out += &format!("s{s}");
            }
        }
        out += "}";
    }
    out
}

// ---------------------------------------------------------------------------------------------
// canonical encodings of Seq-ids

fn canon_str(s: &str) -> Vec<u8> {
    let mut v = vec![b's'];
    v.extend((s.len() as u64).to_le_bytes());
    v.extend(s.as_bytes());
    v
}

/// `local str "<s>"` as `MaskSeqId::Other`: Seq-id variant 0, Object-id variant 1, the string.
fn local_str(s: &str) -> MaskSeqId {
    let mut v = vec![0, 1];
    v.extend(canon_str(s));
    MaskSeqId::Other(v)
}

/// `general { db "<db>", tag str "<tag>" }`: Seq-id variant 10, a SEQUENCE of two members present,
/// the second a CHOICE whose variant 1 is a string.
fn general_str(db: &str, tag: &str) -> MaskSeqId {
    let mut v = vec![10, b'S', 1];
    v.extend(canon_str(db));
    v.extend([1, b'C', 1]);
    v.extend(canon_str(tag));
    MaskSeqId::Other(v)
}

// ---------------------------------------------------------------------------------------------
// the tool outputs

const TOOLS: [&str; 4] = ["dust", "seg", "wm", "c2m"];

#[test]
fn tool_outputs_text_and_binary_agree() {
    for tool in TOOLS {
        for suffix in ["", "_ps"] {
            let text_path = parse_dir().join(format!("{tool}{suffix}.asnt"));
            let bin_path = parse_dir().join(format!("{tool}{suffix}.asnb"));
            let mut text = MaskFileReader::open(&text_path).unwrap();
            let mut bin = MaskFileReader::open(&bin_path).unwrap();
            assert_eq!(text.format(), MaskFileFormat::TextAsn, "{tool}{suffix}");
            assert_eq!(bin.format(), MaskFileFormat::BinaryAsn, "{tool}{suffix}");
            let (text, bin) = (read_all(&mut text).unwrap(), read_all(&mut bin).unwrap());
            assert!(!text.is_empty(), "{tool}{suffix}");
            assert_eq!(text, bin, "{tool}{suffix}");
        }
    }
}

#[test]
fn tool_output_values() {
    // dustmasker on mask_nucl.fna: window=64; level=20; linker=1
    let dust = infos(&parse_dir().join("dust.asnt"));
    assert_eq!(dust.len(), 1);
    assert_eq!(dust[0].algo_program, 10);
    assert!(!dust[0].masks.more);
    let l1 = &dust[0].masks.masks[0];
    assert_eq!(l1.id, local(1));
    assert_eq!(l1.intervals, vec![interval(199, 280)]);
    // sequence 6 has two intervals, in the order of the file
    let l6 = dust[0]
        .masks
        .masks
        .iter()
        .find(|l| l.id == local(6))
        .unwrap();
    assert_eq!(l6.intervals.len(), 2);
    assert!(l6.intervals[0].from < l6.intervals[1].from);

    let seg = infos(&parse_dir().join("seg.asnb"));
    assert_eq!((seg[0].algo_program, seg[0].masks.masks.len()), (20, 3));

    // convert2blastmask: program `other` (100), the options as given, lower-case runs 30..60, 100..120
    let c2m = infos(&parse_dir().join("c2m.asnb"));
    assert_eq!(c2m[0].algo_program, 100);
    assert_eq!(c2m[0].masks.masks.len(), 8);
    assert_eq!(
        c2m[0].masks.masks[2].intervals,
        vec![interval(30, 59), interval(100, 119)]
    );
}

#[test]
fn tool_outputs_with_parsed_seqids() {
    for ext in ["asnt", "asnb"] {
        let dust = infos(&parse_dir().join(format!("dust_ps.{ext}")));
        let locs = &dust[0].masks.masks;
        assert_eq!(locs.len(), 16, "{ext}");
        // lcl|name1, gnl|mydb|tag1, ref|NC_000001.1|, ...
        assert_eq!(locs[0].id, Some(local_str("name1")), "{ext}");
        assert_eq!(locs[1].id, Some(general_str("mydb", "tag1")), "{ext}");
        // lcl|123 is `local id 123`; the other Seq-ids are kept whole
        assert_eq!(locs[10].id, local(123), "{ext}");
        let others: Vec<_> = locs
            .iter()
            .filter(|l| matches!(l.id, Some(MaskSeqId::Other(_))))
            .collect();
        assert_eq!(others.len(), 15, "{ext}");
        // Seq-ids of different kinds differ
        let mut distinct: Vec<_> = others.iter().map(|l| l.id.clone()).collect();
        distinct.sort_by(|a, b| format!("{a:?}").cmp(&format!("{b:?}")));
        distinct.dedup();
        assert_eq!(distinct.len(), 15, "{ext}");
        // pdb|1ABC|A: `pdb { mol "1ABC", chain 65, chain-id "A" }` is the 7th defline
        let Some(MaskSeqId::Other(pdb)) = &locs[6].id else {
            panic!("pdb id");
        };
        assert_eq!(pdb[0], 14, "Seq-id variant pdb");
    }
}

#[test]
fn masks_equal_what_blastdbcmd_prints() {
    let expected = read_expected();
    let mut compared = 0;
    for (file, run) in &expected {
        if run.masks.len() != 1 {
            continue;
        }
        let (algorithm, rows) = run.masks.iter().next().unwrap();
        let fasta = if parse_dir().join(&run.fasta).exists() {
            parse_dir().join(&run.fasta)
        } else {
            seqdb_dir().join(&run.fasta)
        };
        let sequences = read_file(&fasta)
            .split(|&c| c == b'\n')
            .filter(|l| l.starts_with(b">"))
            .count() as u32;
        let objects = match open_all(&parse_dir().join(file)) {
            // NCBI reads these; the rows are not checked
            Err(MaskFileError::Unsupported(_)) => continue,
            other => other.unwrap(),
        };
        let mut locs: Vec<&MaskLoc> = Vec::new();
        for (info, lists) in &objects {
            locs.extend(&info.masks.masks);
            for list in lists {
                locs.extend(&list.masks);
            }
        }
        for oid in 0..sequences {
            let mine: Vec<_> = locs
                .iter()
                .filter(|l| l.id == local(oid as i32 + 1))
                .collect();
            let expected_row = match mine.as_slice() {
                [] => "none".to_string(),
                // Several Seq-locs of one sequence are merged and sorted by NCBI (not here).
                [_, _, ..] => continue,
                [one] => one
                    .intervals
                    .iter()
                    .map(|iv| format!("{}-{};", iv.from, u64::from(iv.to) + 1))
                    .collect(),
            };
            assert_eq!(
                rows.get(&oid),
                Some(&expected_row),
                "{file}, algorithm {algorithm}, sequence {oid}"
            );
            compared += 1;
        }
    }
    assert!(compared > 100, "{compared} sequences compared");
}

// ---------------------------------------------------------------------------------------------
// what makeblastdb makes of each file

#[derive(Debug, PartialEq)]
enum Outcome {
    /// "Unknown encoding for mask data in ..." (rc 1)
    Unknown,
    /// the serial reader throws (rc 255)
    Serial,
    /// bad_alloc (rc 4)
    OutOfMemory,
    /// the build is abandoned without an error status (a missing continuation)
    NoDb,
    /// the objects are read (a later step may still fail: rc 1)
    Accepted,
}

fn outcome(run: &Run) -> Outcome {
    if run.rc == 1 && run.stderr.contains("Unknown encoding for mask data in '") {
        Outcome::Unknown
    } else if run.rc == 255 && run.stderr.starts_with("Error: NCBI C++ Exception:") {
        Outcome::Serial
    } else if run.rc == 4 {
        Outcome::OutOfMemory
    } else if run.rc == 0 && run.no_db {
        Outcome::NoDb
    } else {
        assert!(
            run.rc == 0 || run.rc == 1,
            "unexpected exit status {}: {}",
            run.rc,
            run.stderr
        );
        assert!(!run.stderr.contains("NCBI C++ Exception"), "{}", run.stderr);
        Outcome::Accepted
    }
}

/// Files that makeblastdb reads and LOSAT does not.
const UNSUPPORTED: [&str; 9] = [
    "h_unsupported_whole.asnt",
    "h_unsupported_pnt.asnt",
    "h_unsupported_mix.asnt",
    "h_unsupported_null.asnt",
    "h_unsupported_feat.asnt",
    "b_small_whole.asnb",
    // NCBI replaces a byte that is not printable in a VisibleString with '#' and logs it, and
    // warns about a UTF8String tag for a VisibleString; LOSAT does not write those lines
    "h_string_nonprint.asnt",
    "h_utf8_in_options.asnt",
    "b_utf8_options.asnb",
];

#[test]
fn every_file_is_read_as_makeblastdb_reads_it() {
    let expected = read_expected();
    assert!(expected.len() > 120);
    let mut seen = BTreeMap::new();
    for (file, run) in &expected {
        let path = parse_dir().join(file);
        let outcome = outcome(run);
        *seen.entry(format!("{outcome:?}")).or_insert(0) += 1;
        let result = open_all(&path);
        match outcome {
            Outcome::Unknown => {
                let Err(MaskFileError::Ncbi(WriteDbError::Options(text))) =
                    MaskFileReader::open(&path)
                else {
                    panic!("{file}: not rejected by the format guess");
                };
                // NCBI's stderr names the file as it was given; here the path is absolute.
                assert_eq!(
                    run.stderr,
                    format!(
                        "BLAST options error: {}\n",
                        text.replace(&path.display().to_string(), file)
                    ),
                    "{file}"
                );
            }
            Outcome::Serial | Outcome::OutOfMemory => {
                assert!(
                    matches!(result, Err(MaskFileError::Malformed { .. })),
                    "{file}: {result:?}"
                );
            }
            Outcome::NoDb => {
                assert!(
                    matches!(result, Err(MaskFileError::EofInContinuation { .. })),
                    "{file}: {result:?}"
                );
            }
            Outcome::Accepted
                if UNSUPPORTED.contains(&file.as_str()) || file == "h_different_other_ids.asnt" =>
            {
                assert!(
                    matches!(result, Err(MaskFileError::Unsupported(_))),
                    "{file}: {result:?}"
                );
            }
            Outcome::Accepted => {
                result.unwrap_or_else(|e| panic!("{file}: {e}"));
            }
        }
    }
    // every outcome is covered
    for key in ["Unknown", "Serial", "OutOfMemory", "NoDb", "Accepted"] {
        assert!(
            seen.contains_key(key),
            "no file with outcome {key}: {seen:?}"
        );
    }
}

#[test]
fn unsupported_files_name_what_is_not_read() {
    for (file, what) in [
        ("h_unsupported_whole.asnt", "`whole`"),
        ("h_unsupported_pnt.asnt", "`pnt`"),
        ("h_unsupported_mix.asnt", "`mix`"),
        ("h_unsupported_null.asnt", "`null`"),
        ("h_unsupported_feat.asnt", "`feat`"),
        ("b_small_whole.asnb", "`whole`"),
        ("h_different_other_ids.asnt", "different Seq-ids"),
        (
            "h_string_nonprint.asnt",
            "byte 0x09 in the VisibleString member `algo-options`",
        ),
        (
            "h_utf8_in_options.asnt",
            "byte 0xc3 in the VisibleString member `algo-options`",
        ),
        (
            "b_utf8_options.asnb",
            "a UTF8String tag for the VisibleString member `algo-options`",
        ),
    ] {
        let Err(MaskFileError::Unsupported(text)) = open_all(&parse_dir().join(file)) else {
            panic!("{file}");
        };
        assert!(text.contains(what), "{file}: {text}");
    }
}

#[test]
fn seq_ids_of_one_seq_loc() {
    let one = |file: &str| -> Vec<MaskLoc> { infos(&parse_dir().join(file)).remove(0).masks.masks };
    // two different `local id`: no Seq-id (makeblastdb: "Cannot get masked range Seq-id")
    assert_eq!(one("h_packed_two_ids.asnt")[0].id, None);
    assert_eq!(one("h_packed_two_ids.asnt")[0].intervals.len(), 2);
    assert_eq!(one("b_small_two_ids.asnb")[0].id, None);
    // an empty packed-int has none either
    assert_eq!(one("h_packed_empty.asnt")[0].id, None);
    assert!(one("h_packed_empty.asnt")[0].intervals.is_empty());
    // a `local id` and another kind
    assert_eq!(one("h_local_and_other_id.asnt")[0].id, None);
    // the same `gi 5`, written with other spacing and a comment: equal
    let same = one("h_same_other_id_spacing.asnt");
    assert!(matches!(same[0].id, Some(MaskSeqId::Other(_))));
    assert_eq!(same[0].intervals, vec![interval(5, 14), interval(20, 29)]);
    // `local str "1"` is not `local id 1`
    assert!(matches!(
        one("h_local_str_1.asnt")[0].id,
        Some(MaskSeqId::Other(_))
    ));
    // a negative `local id` is a number like any other
    assert_eq!(one("h_id_negative.asnt")[0].id, local(-3));
    // `int` and `packed-int` of one interval read alike
    assert_eq!(one("h_int_variant.asnt"), one("h_ok.asnt"));
    // the value of a stand-alone Seq-id does not depend on the encoding: gi 123456789012 in BER
    let Some(MaskSeqId::Other(gi)) = &one("b_small_gi.asnb")[0].id else {
        panic!("gi");
    };
    let mut expected = vec![11, b'i'];
    expected.extend(123_456_789_012i64.to_le_bytes());
    assert_eq!(gi, &expected);
}

#[test]
fn hand_written_text_values() {
    let read = |file: &str| infos(&parse_dir().join(file));
    let ok = read("h_ok.asnt");
    assert_eq!(
        ok,
        vec![MaskInfo {
            algo_id: 0,
            algo_program: 30,
            algo_options: b"x".to_vec(),
            masks: MaskList {
                masks: vec![MaskLoc {
                    id: local(3),
                    intervals: vec![interval(5, 14)],
                }],
                more: false,
            },
        }]
    );
    // comments, other white space, brackets around names, a vertical tab in front of a number
    // (`GetInt4` skips C white space) and a NUL in a comment do not change the value
    for file in [
        "h_comments.asnt",
        "h_crlf.asnt",
        "h_tab_and_cr.asnt",
        "h_bracket_member_names.asnt",
        "h_int_variant.asnt",
        "h_long_comment_before_object.asnt",
        "h_nul_in_comment.asnt",
        "h_vertical_tab_between_tokens.asnt",
    ] {
        assert_eq!(read(file), ok, "{file}");
    }
    // strings
    // enumerations
    assert_eq!(read("h_algo_program_number.asnt")[0].algo_program, 30);
    assert_eq!(read("h_algo_program_dust.asnt")[0].algo_program, 10);
    assert_eq!(read("h_algo_program_bad_number.asnt")[0].algo_program, 31);
    let strand = |file: &str| read(file)[0].masks.masks[0].intervals[0].strand;
    assert_eq!(strand("h_strand_plus.asnt"), Some(1));
    assert_eq!(strand("h_strand_number.asnt"), Some(3));
    assert_eq!(strand("h_ok.asnt"), None);
    let minus = &read("h_strand_minus_two.asnt")[0].masks.masks[0];
    assert_eq!(minus.id, local(3));
    assert_eq!(
        minus.intervals,
        vec![
            MaskInterval {
                from: 20,
                to: 29,
                strand: Some(2)
            },
            MaskInterval {
                from: 5,
                to: 14,
                strand: Some(2)
            },
        ]
    );
    // fuzz is read and dropped
    assert_eq!(read("h_fuzz_from.asnt"), ok);
    // numbers
    assert_eq!(read("h_algo_id_negative.asnt")[0].algo_id, i32::MIN);
    assert_eq!(
        read("h_from_after_to.asnt")[0].masks.masks[0].intervals,
        vec![interval(14, 5)]
    );
    assert_eq!(
        read("h_to_beyond_sequence.asnt")[0].masks.masks[0].intervals,
        vec![interval(5, 5000)]
    );
    // empty masks
    assert!(read("h_empty_masks.asnt")[0].masks.masks.is_empty());
}

#[test]
fn objects_continuations_and_the_end_of_the_file() {
    let path = |f: &str| parse_dir().join(f);
    // two objects: the second is read by a later call
    let two = open_all(&path("h_two_objects.asnt")).unwrap();
    assert_eq!(two.len(), 2);
    assert_eq!(
        (
            two[1].0.algo_id,
            two[1].0.algo_program,
            &two[1].0.algo_options[..]
        ),
        (1, 10, &b"y"[..])
    );
    // an object without a type name, as makeblastdb reads it
    assert_eq!(
        open_all(&path("h_two_objects_headless.asnt"))
            .unwrap()
            .len(),
        2
    );
    // more TRUE: the list that follows, then the list that it announces
    let cont = open_all(&path("h_more_true_two_continuations.asnt")).unwrap();
    assert_eq!(cont.len(), 1);
    assert!(cont[0].0.masks.more);
    assert_eq!(cont[0].1.len(), 2);
    assert_eq!(cont[0].1[0].masks[0].id, local(4));
    assert!(cont[0].1[0].more && !cont[0].1[1].more);
    // a missing continuation
    for file in [
        "h_more_true_no_continuation.asnt",
        "h_more_true_continuation_header_only.asnt",
        "b_small_more_true.asnb",
    ] {
        let mut r = MaskFileReader::open(&path(file)).unwrap();
        let info = r.next_info().unwrap().unwrap();
        assert!(info.masks.more, "{file}");
        assert!(
            matches!(r.next_list(), Err(MaskFileError::EofInContinuation { .. })),
            "{file}"
        );
        assert!(
            matches!(r.next_info(), Err(MaskFileError::EofInContinuation { .. })),
            "{file}: reading stops"
        );
    }
    // a continuation read as an object, and an object read as a continuation
    assert!(matches!(
        open_all(&path("h_more_true_wrong_type.asnt")),
        Err(MaskFileError::Malformed { .. })
    ));
    let mut r = MaskFileReader::open(&path("h_ok.asnt")).unwrap();
    assert!(
        matches!(r.next_list(), Err(MaskFileError::Malformed { .. })),
        "an object is not a list"
    );
    // the end of the file: nothing but comments and white space, or a file that ends in the type
    // name, after `::=`, or in the first bytes of a binary object, is the end
    for file in [
        "h_cut_after_header.asnt",
        "h_cut_in_type_name.asnt",
        "h_garbage_after_object.asnt",
        "h_comments.asnt",
        "b_cut_before_last_eoc_byte.asnb",
        "b_trailing_sequence_tag.asnb",
    ] {
        let mut r = MaskFileReader::open(&path(file)).unwrap();
        let first = r.next_info();
        if file == "b_cut_before_last_eoc_byte.asnb" {
            // the object being read is dropped
            assert!(first.unwrap().is_none(), "{file}");
        } else {
            assert!(first.unwrap().is_some(), "{file}");
            assert!(r.next_info().unwrap().is_none(), "{file}");
        }
        assert!(r.next_info().unwrap().is_none(), "{file}");
    }
    // ... but inside an object it is an error, and so is other text after the object
    for file in [
        "h_cut_after_brace.asnt",
        "h_cut_inside_object.asnt",
        "h_garbage_after_object_2.asnt",
        "h_trailing_brace.asnt",
        "b_trailing_zero.asnb",
        "b_trailing_sequence_start.asnb",
    ] {
        assert!(
            open_all(&path(file)).is_err_and(|e| matches!(e, MaskFileError::Malformed { .. })),
            "{file}"
        );
    }
}

#[test]
fn files_that_are_neither_encoding() {
    let expected = read_expected();
    for file in [
        "h_not_asn_fasta.fna",
        "h_not_asn_text.txt",
        "h_not_asn_empty.txt",
        "h_not_asn_comment_only.txt",
        "h_not_asn_comment_only_no_newline.txt",
        "h_not_asn_control_one.bin",
        "h_comment_mostly_punctuation.asnt",
    ] {
        assert_eq!(outcome(&expected[file]), Outcome::Unknown, "{file}");
        assert!(matches!(
            MaskFileReader::open(&parse_dir().join(file)),
            Err(MaskFileError::Ncbi(WriteDbError::Options(_)))
        ));
    }
    // a FASTA file given as a mask file, with the text that makeblastdb writes
    let fasta = seqdb_dir().join("mask_nucl.fna");
    let Err(error) = MaskFileReader::open(&fasta) else {
        panic!("a FASTA file is not a mask file");
    };
    assert_eq!(
        error.to_string(),
        format!(
            "BLAST options error: Unknown encoding for mask data in '{}'",
            fasta.display()
        )
    );
    // a byte-order mark and a control byte other than 0x01 make a file binary
    for file in ["h_bom_before_object.asnt", "h_binary_control_two.bin"] {
        let r = MaskFileReader::open(&parse_dir().join(file)).unwrap();
        assert_eq!(r.format(), MaskFileFormat::BinaryAsn, "{file}");
    }
    // a file that cannot be opened (or a directory, which cannot be read) is a failed stream for
    // the format guess: the same error
    for path in [parse_dir().join("no_such_file"), parse_dir()] {
        let Err(MaskFileError::Ncbi(WriteDbError::Options(text))) = MaskFileReader::open(&path)
        else {
            panic!("{}", path.display());
        };
        assert_eq!(
            text,
            format!("Unknown encoding for mask data in '{}'", path.display())
        );
    }
}

// ---------------------------------------------------------------------------------------------
// the fixtures of make_mask_fixtures.py

#[test]
fn windowmasker_fixtures() {
    let wm = |file: &str| infos(&seqdb_dir().join(file));
    let expect = |entries: &[(i32, &[(u32, u32)])]| -> Vec<MaskLoc> {
        entries
            .iter()
            .map(|&(id, ranges)| MaskLoc {
                id: local(id),
                intervals: ranges.iter().map(|&(a, b)| interval(a, b)).collect(),
            })
            .collect()
    };
    let m = wm("mask_wm.asnt");
    assert_eq!(m.len(), 1);
    assert_eq!(
        (m[0].algo_id, m[0].algo_program, &m[0].algo_options[..]),
        (0, 30, &b""[..])
    );
    assert_eq!(
        m[0].masks.masks,
        expect(&[
            (1, &[(0, 9), (300, 309), (500, 529)]),
            (2, &[(349, 349)]),
            (4, &[(10, 19), (20, 29), (380, 389)]),
            (5, &[(0, 0)]),
            (7, &[(0, 119), (200, 319)]),
        ])
    );
    assert!(!m[0].masks.more);
    let m = wm("mask_wm_opts.asnt");
    assert_eq!(&m[0].algo_options[..], b"-t_extend 3");
    assert_eq!(
        m[0].masks.masks,
        expect(&[(3, &[(5, 14)]), (6, &[(0, 99)])])
    );
    let m = wm("mask_wm_colon.asnt");
    assert_eq!(&m[0].algo_options[..], b"-x a:b");
    assert_eq!(
        m[0].masks.masks,
        expect(&[(1, &[(0, 9), (300, 309)]), (4, &[(10, 19)])])
    );
}

/// The rows `algorithm<TAB>oid<TAB>ranges` of an `expected/*.masks.tsv` for one algorithm.
fn blastdbcmd_rows(file: &str, algorithm: u32) -> BTreeMap<u32, String> {
    let text = String::from_utf8(read_file(&seqdb_dir().join("expected").join(file))).unwrap();
    text.lines()
        .filter_map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].parse::<u32>().unwrap() == algorithm)
                .then(|| (f[1].parse().unwrap(), f[2].to_string()))
        })
        .collect()
}

#[test]
fn dust_and_seg_fixtures() {
    for (file, expected, algorithm, program, sequences) in [
        ("mask_dust.asnb", "mask_nucl_v4.masks.tsv", 11, 10, 8),
        ("mask_seg.asnb", "mask_prot_v4.masks.tsv", 21, 20, 5),
    ] {
        let objects = infos(&seqdb_dir().join(file));
        assert_eq!(objects.len(), 1, "{file}");
        assert_eq!(objects[0].algo_program, program, "{file}");
        let rows = blastdbcmd_rows(expected, algorithm);
        assert_eq!(rows.len(), sequences, "{file}");
        for oid in 0..sequences as u32 {
            let locs: Vec<_> = objects[0]
                .masks
                .masks
                .iter()
                .filter(|l| l.id == local(oid as i32 + 1))
                .collect();
            assert!(locs.len() <= 1);
            let row = locs
                .first()
                .map(|l| {
                    l.intervals
                        .iter()
                        .map(|iv| format!("{}-{};", iv.from, iv.to + 1))
                        .collect()
                })
                .unwrap_or_else(|| "none".to_string());
            assert_eq!(rows[&oid], row, "{file}, sequence {oid}");
        }
    }
    // the text and the binary form of the same dustmasker output agree
    assert_eq!(
        infos(&seqdb_dir().join("mask_dust.asnb")),
        infos(&parse_dir().join("dust.asnb"))
    );
}

// ---------------------------------------------------------------------------------------------
// damaged files

fn read_bytes(label: &str, data: &[u8]) -> Result<Vec<(MaskInfo, Vec<MaskList>)>, MaskFileError> {
    read_all(&mut MaskFileReader::from_bytes(
        Path::new(label),
        data.to_vec(),
    )?)
}

#[test]
fn every_truncation_of_a_binary_file() {
    let dust = read_file(&seqdb_dir().join("mask_dust.asnb"));
    assert_eq!(dust.len(), 307);
    for n in 0..=dust.len() {
        let result = read_bytes("cut", &dust[..n]);
        match n {
            // no byte that is not a graphic character: neither format
            0 | 1 => assert!(
                matches!(result, Err(MaskFileError::Ncbi(_))),
                "{n}: {result:?}"
            ),
            // the file ends where the outermost SEQUENCE has no closing byte: the end of the file
            306 => assert_eq!(result.unwrap().len(), 0, "{n}"),
            307 => assert_eq!(result.unwrap().len(), 1, "{n}"),
            _ => assert!(
                matches!(result, Err(MaskFileError::Malformed { .. })),
                "{n}: {result:?}"
            ),
        }
    }
}

#[test]
fn damaged_binary_files() {
    let dust = read_file(&seqdb_dir().join("mask_dust.asnb"));
    let malformed = |label: &str, data: &[u8]| {
        let r = read_bytes(label, data);
        assert!(
            matches!(r, Err(MaskFileError::Malformed { .. })),
            "{label}: {r:?}"
        );
    };
    // a length larger than the file (an options string of 2^28 bytes), larger than a `size_t`
    // can hold the length of, and with the first length byte zero
    let with = |from: &[u8], to: &[u8]| -> Vec<u8> {
        let at = dust.windows(from.len()).position(|w| w == from).unwrap();
        [&dust[..at], to, &dust[at + from.len()..]].concat()
    };
    malformed(
        "beyond",
        &with(&[0x1a, 0x1d], &[0x1a, 0x84, 0x10, 0x00, 0x00, 0x00]),
    );
    malformed(
        "huge",
        &with(
            &[0x1a, 0x1d],
            &[0x1a, 0x88, 0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff],
        ),
    );
    malformed(
        "nine",
        &with(&[0x1a, 0x1d], &[0x1a, 0x89, 1, 0, 0, 0, 0, 0, 0, 0, 0]),
    );
    malformed(
        "zero first",
        &with(&[0x1a, 0x1d], &[0x1a, 0x82, 0x00, 0x1d]),
    );
    malformed("indefinite string", &with(&[0x1a, 0x1d], &[0x1a, 0x80]));
    // an indefinite length that is never closed (the closing bytes cut off)
    malformed("never closed", &dust[..dust.len() - 2]);
    // a ten-byte tag number
    malformed(
        "tag",
        &with(
            &[0xa0, 0x80, 0x02, 0x01, 0x02],
            &[
                0xbf, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x00, 0x80, 0x02, 0x01,
                0x02,
            ],
        ),
    );
    // a long-form tag number with the same value is read
    let ok = read_bytes(
        "tag",
        &with(
            &[0xa0, 0x80, 0x02, 0x01, 0x02],
            &[0xbf, 0x00, 0x80, 0x02, 0x01, 0x02],
        ),
    )
    .unwrap();
    assert_eq!(ok, read_bytes("dust", &dust).unwrap());
    // a length that is not the end of the item (a long-form definite length counts from the length
    // byte, so every such item is refused)
    let mut definite = vec![0x30, 0x82];
    definite.extend(((dust.len() - 4) as u16).to_be_bytes());
    definite.extend(&dust[2..dust.len() - 2]);
    malformed("definite long form", &definite);
    // wrong tags and bytes at the end
    malformed("first tag", &[&[0x31][..], &dust[1..]].concat());
    malformed("trailing zero", &[&dust[..], &[0u8][..]].concat());
    assert_eq!(
        read_bytes("two", &[&dust[..], &dust[..]].concat())
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn binary_edits_never_panic() {
    let dust = read_file(&seqdb_dir().join("mask_dust.asnb"));
    let mut results = (0, 0);
    for i in 0..dust.len() {
        for value in [0u8, 1, 0x1f, 0x30, 0x80, 0x81, 0xa0, 0xff, dust[i] ^ 1] {
            let mut data = dust.clone();
            data[i] = value;
            match read_bytes("flip", &data) {
                Ok(_) => results.0 += 1,
                Err(_) => results.1 += 1,
            }
        }
    }
    assert!(results.0 > 0 && results.1 > results.0, "{results:?}");
    // inserted and removed bytes
    for i in 0..dust.len() {
        let mut data = dust.clone();
        data.remove(i);
        let _ = read_bytes("remove", &data);
        let mut data = dust.clone();
        data.insert(i, 0x80);
        let _ = read_bytes("insert", &data);
    }
}

/// How `read_bytes` ends, for the truncation tables.
#[derive(Debug, PartialEq)]
enum Cut {
    Unknown,
    Malformed,
    Objects(usize),
    MissingContinuation,
}

fn cut_result(data: &[u8]) -> Cut {
    match read_bytes("cut", data) {
        Err(MaskFileError::Ncbi(_)) => Cut::Unknown,
        Err(MaskFileError::Malformed { .. }) => Cut::Malformed,
        Err(MaskFileError::EofInContinuation { .. }) => Cut::MissingContinuation,
        Err(other) => panic!("{other}"),
        Ok(objects) => Cut::Objects(objects.len()),
    }
}

/// Every prefix of a text file, against the table of what makeblastdb 2.17.0 does with it:
/// `unknown`, `ok` and `missing` are ranges of prefix lengths, the rest is refused.
fn check_text_prefixes(file: &str, unknown: Range, ok: &[(Range, usize)], missing: &[Range]) {
    let data = read_file(&parse_dir().join(file));
    for n in 0..=data.len() {
        let in_range = |r: &Range| (r.0..=r.1).contains(&n);
        let expected = if in_range(&unknown) {
            Cut::Unknown
        } else if let Some((_, count)) = ok.iter().find(|(r, _)| in_range(r)) {
            Cut::Objects(*count)
        } else if missing.iter().any(in_range) {
            Cut::MissingContinuation
        } else {
            Cut::Malformed
        };
        assert_eq!(cut_result(&data[..n]), expected, "{file}, {n} bytes");
    }
}

type Range = (usize, usize);

#[test]
fn every_truncation_of_text_files() {
    check_text_prefixes("h_ok.asnt", (0, 24), &[((163, 164), 1)], &[]);
    // inside the first comment, then inside the object, then after it the comment and white space
    check_text_prefixes(
        "h_comments.asnt",
        (0, 52),
        &[((53, 54), 0), ((232, 294), 1)],
        &[],
    );
    check_text_prefixes(
        "h_two_objects.asnt",
        (0, 24),
        &[((163, 183), 1), ((186, 187), 1), ((320, 321), 2)],
        &[],
    );
    check_text_prefixes(
        "h_more_true_continuation.asnt",
        (0, 24),
        &[((258, 259), 1)],
        &[(162, 179), (182, 183)],
    );
}

#[test]
fn text_edits_never_panic() {
    let data = read_file(&parse_dir().join("h_comments.asnt"));
    for i in 0..data.len() {
        for value in [
            0u8, b'"', b'-', b'{', b'}', b',', b'0', 0xff, b' ', b'\n', b'\r',
        ] {
            let mut edited = data.clone();
            edited[i] = value;
            let _ = read_bytes("edit", &edited);
        }
        let mut edited = data.clone();
        edited.remove(i);
        let _ = read_bytes("remove", &edited);
    }
    // nesting deeper than NCBI would need for any real file
    let mut deep = String::from("Blast-db-mask-info ::= { algo-id 0, algo-program windowmasker, algo-options \"\", masks { masks { ");
    for _ in 0..2000 {
        deep.push_str("mix { ");
    }
    assert!(matches!(
        read_bytes("deep", deep.as_bytes()),
        Err(MaskFileError::Unsupported(_)) | Err(MaskFileError::Malformed { .. })
    ));
}

#[test]
fn visible_strings_replace_what_is_not_printable() {
    let mut s = *b"a\tb\xe9 \x7f~\x00";
    super::reader::fix_visible_chars(&mut s);
    assert_eq!(&s, b"a#b# #~#");
}

/// The options of the first object, against what `blastdbcmd -info` prints for the database that
/// makeblastdb wrote from the file (`expected.tsv`, last field of `RUN`).
#[test]
fn algorithm_options_equal_what_blastdbcmd_prints() {
    let mut compared = 0;
    for (file, run) in read_expected() {
        if run.options.is_empty() {
            continue;
        }
        let objects = open_all(&parse_dir().join(&file));
        if matches!(objects, Err(MaskFileError::Unsupported(_))) {
            continue; // NCBI's text has '#' for bytes that LOSAT refuses
        }
        let options = String::from_utf8(objects.unwrap()[0].0.algo_options.clone()).unwrap();
        // The table prints an empty string as "default options used".
        let options = if options.is_empty() {
            "default options used".to_string()
        } else {
            options
        };
        assert_eq!(options, run.options, "{file}");
        compared += 1;
    }
    assert!(compared > 40, "{compared}");
}

/// After an error the reader does not go on, and does not report a clean end.
#[test]
fn an_error_is_final() {
    let mut r = MaskFileReader::open(&parse_dir().join("h_cut_inside_object.asnt")).unwrap();
    assert!(matches!(
        r.next_info(),
        Err(MaskFileError::Malformed { .. })
    ));
    for _ in 0..2 {
        assert!(matches!(
            r.next_info(),
            Err(MaskFileError::Malformed { .. })
        ));
        assert!(matches!(
            r.next_list(),
            Err(MaskFileError::Malformed { .. })
        ));
    }
    // a refusal as well
    let mut r = MaskFileReader::open(&parse_dir().join("h_unsupported_whole.asnt")).unwrap();
    assert!(matches!(r.next_info(), Err(MaskFileError::Unsupported(_))));
    assert!(matches!(
        r.next_info(),
        Err(MaskFileError::Malformed { .. })
    ));
    // a missing continuation
    let mut r =
        MaskFileReader::open(&parse_dir().join("h_more_true_no_continuation.asnt")).unwrap();
    r.next_info().unwrap().unwrap();
    for _ in 0..2 {
        assert!(matches!(
            r.next_list(),
            Err(MaskFileError::EofInContinuation { .. })
        ));
        assert!(matches!(
            r.next_info(),
            Err(MaskFileError::EofInContinuation { .. })
        ));
    }
    // a clean end stays one
    let mut r = MaskFileReader::open(&parse_dir().join("h_ok.asnt")).unwrap();
    assert!(r.next_info().unwrap().is_some());
    assert!(r.next_info().unwrap().is_none());
    assert!(r.next_info().unwrap().is_none());
}
