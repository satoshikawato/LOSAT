//! `losat blastdbcmd`: the retrieval of sequences (`-entry`, `-entry_batch`, `-outfmt`,
//! `-range`, `-strand`, `-line_length`, `-ctrl_a`, `-target_only`, `-get_dups`,
//! `-long_seqids`) against NCBI BLAST+ 2.17.0 blastdbcmd. The expected bytes are NCBI's,
//! frozen by line BDC's probe on ws (`tests/blastdbcmd/README.md`), or, for the approved
//! differences (line BDC DECISIONS 12), the correct output built from NCBI's neighbouring
//! runs; `{R}` stands for the directory of the probe's runs (here `tests/blastdbcmd`, the
//! current directory of each run). NCBI is not run here.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn tests_dir() -> PathBuf {
    fs::canonicalize(Path::new(env!("CARGO_MANIFEST_DIR")).join("tests")).unwrap()
}

/// A directory of its own for one test (the tests of a binary run in parallel), removed at
/// the end.
struct TestDir(PathBuf);

impl TestDir {
    fn new(stem: &str) -> Self {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "losat_blastdbcmd_{stem}_{}_{nanos}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).unwrap();
        TestDir(fs::canonicalize(dir).unwrap())
    }
}

impl Drop for TestDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// One run of `losat blastdbcmd words...` in `tests/blastdbcmd`, `{O}` in a word standing for
/// `out_dir`: stdout, stderr, exit status. No `BLASTDB`, no registry file (`NCBI_CONFIG_PATH`
/// names `out_dir` alone), no LOSAT variable.
fn blastdbcmd(words: &[&str], out_dir: &Path) -> (Vec<u8>, Vec<u8>, i32) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text == "BLASTDB" {
            command.env_remove(&name);
        }
    }
    let out_dir = out_dir.to_str().unwrap();
    let output = command
        .env("NCBI_CONFIG_PATH", out_dir)
        .current_dir(tests_dir().join("blastdbcmd"))
        .arg("blastdbcmd")
        .args(words.iter().map(|word| word.replace("{O}", out_dir)))
        .output()
        .unwrap();
    (output.stdout, output.stderr, output.status.code().unwrap())
}

/// The frozen file `name` of `tests/blastdbcmd/expected`, with `{R}` placed.
fn expected(name: &str) -> Vec<u8> {
    let tests = tests_dir();
    let text = fs::read(tests.join("blastdbcmd/expected").join(name)).unwrap();
    let text = String::from_utf8(text).unwrap();
    text.replace("{R}", tests.join("blastdbcmd").to_str().unwrap())
        .into_bytes()
}

fn check(id: &str, words: &[&str], status: i32, out_dir: &Path) {
    let (stdout, stderr, code) = blastdbcmd(words, out_dir);
    assert_eq!(
        String::from_utf8_lossy(&stdout),
        String::from_utf8_lossy(&expected(&format!("{id}.stdout"))),
        "{id} {words:?}: stdout"
    );
    assert_eq!(
        String::from_utf8_lossy(&stderr),
        String::from_utf8_lossy(&expected(&format!("{id}.stderr"))),
        "{id} {words:?}: stderr"
    );
    assert_eq!(code, status, "{id} {words:?}: exit status");
}

/// The probe's cases that LOSAT answers as NCBI does: id, words, NCBI's exit status. The
/// databases are the probe's (`db/`), without parsed Seq-ids (`p4n`, `p5n`, `n4n`, `n5n`,
/// `pL`, `nL`) or with them where no id is looked up (`-entry all`: `p4`, `p5`, `n4`, `n5`).
#[rustfmt::skip]
const MATCHES: &[(&str, &[&str], i32)] = &[
    ("F01", &["-db", "db/p5n", "-entry", "all"], 0),
    ("F02", &["-db", "db/p4n", "-entry", "all"], 0),
    ("F03", &["-db", "db/pL", "-entry", "all"], 0),
    ("F04", &["-db", "db/p5n", "-entry", "gnl|BL_ORD_ID|3"], 1),
    ("F05", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3"], 0),
    ("F06", &["-db", "db/p4n", "-entry", "3"], 1),
    ("F07", &["-db", "db/p5n", "-entry", "myprot9"], 1),
    ("F08", &["-db", "db/p4n", "-entry", "myprot9"], 1),
    ("F13", &["-db", "db/p4", "-entry", "all"], 0),
    ("F14", &["-db", "db/p5", "-entry", "all"], 0),
    ("F15", &["-db", "db/p4", "-entry", "all", "-outfmt", "%o|%a|%i|%g|%t|%T|%e|%l|%h|%P"], 0),
    ("F16", &["-db", "db/p5", "-entry", "all", "-outfmt", "%o|%a|%i|%g|%t|%T|%e|%l|%h|%P"], 0),
    ("F21", &["-db", "db/n4", "-entry", "all"], 0),
    ("F22", &["-db", "db/n5", "-entry", "all"], 0),
    ("F25", &["-db", "db/n4n", "-entry", "all", "-line_length", "70"], 0),
    ("F26", &["-db", "db/n5n", "-entry", "all", "-outfmt", "%o|%a|%i|%g|%t|%l"], 0),
    ("F27", &["-db", "db/n5n", "-entry", "gnl|BL_ORD_ID|0", "-range", "5-20", "-strand", "minus"], 1),
    ("F28", &["-db", "db/nL", "-entry", "all", "-range", "3-9", "-strand", "minus"], 0),
    ("S20", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-long_seqids"], 0),
    ("S44", &["-db", "db/p5n", "-entry", "all", "-outfmt", "%a|%i|%g|%T|%P|%e"], 0),
    ("S45", &["-db", "db/n5n", "-entry_batch", "in/ordid.txt"], 1),
    ("D01", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-outfmt", "%a %Z"], 1),
    ("D02", &["-db", "db/p4n", "-entry", "all", "-outfmt", "%f %b"], 1),
    ("D03", &["-db", "db/p4n", "-entry", "all", "-outfmt", "%a %m"], 1),
    ("D04", &["-db", "db/p4n", "-entry", "all", "-outfmt", "%a %m11"], 1),
    ("D05", &["-db", "db/p4n", "-entry", "all", "-mask_sequence_with", "11"], 1),
    ("D06", &["-db", "db/p4n", "-entry", "all", "-mask_sequence_with", "dust"], 1),
    ("D07", &["-db", "db/p4n", "-entry", "all", "-outfmt", ""], 1),
    ("D08", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-outfmt", "%"], 1),
    ("D09", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-outfmt", "a%"], 1),
    ("D10", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-outfmt", "%%"], 1),
    ("D11", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|0", "-outfmt", "%s"], 1),
    ("D12", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|0", "-outfmt", "%s", "-range", "5-20", "-strand", "minus"], 1),
    ("D13", &["-db", "db/n4n", "-entry", "all", "-outfmt", "%o %l %h"], 0),
    ("D14", &["-db", "db/p4n", "-entry", "all", "-outfmt", "%s %h"], 0),
    ("D15", &["-db", "db/pL", "-entry", "all", "-outfmt", "%o %s"], 0),
    ("D16", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|0", "-range", "5-20", "-strand", "minus", "-long_seqids"], 1),
    ("D17", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|0", "-range", "5-20", "-strand", "minus"], 1),
    ("D18", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|0", "-range", "5-20", "-long_seqids"], 1),
    ("D19", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|0,gnl|BL_ORD_ID|0", "-long_seqids"], 1),
    ("D20", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|0,gnl|BL_ORD_ID|0"], 1),
    ("D21", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-target_only"], 0),
    ("D22", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-target_only", "-outfmt", "%a %t"], 0),
    ("D23", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3,nosuch,gnl|BL_ORD_ID|1"], 1),
    ("D24", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|999"], 1),
    ("D25", &["-db", "db/p4n", "-entry", ""], 1),
    ("D26", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3,,gnl|BL_ORD_ID|1"], 1),
    ("D27", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3, gnl|BL_ORD_ID|1"], 1),
    ("D28", &["-db", "db/p4n", "-entry_batch", "in/ordid_d.txt"], 1),
    ("D29", &["-db", "db/p4n", "-entry_batch", "in/ordid_d.txt", "-get_dups"], 1),
    ("D30", &["-db", "db/n4n", "-entry_batch", "in/ordid_n.txt"], 1),
    ("D31", &["-db", "db/n4n", "-entry_batch", "in/ordid_n.txt", "-get_dups"], 1),
    ("D32", &["-db", "db/n5n", "-entry_batch", "in/ordid_n.txt"], 1),
    ("D33", &["-db", "db/n5n", "-entry_batch", "in/ordid_n.txt", "-get_dups"], 1),
    ("D34", &["-db", "db/p5n", "-entry", "3"], 1),
    ("D35", &["-db", "db/p5n", "-entry", ""], 1),
    ("D36", &["-db", "db/n4n", "-entry", "all", "-strand", "minus"], 0),
    ("D37", &["-db", "db/n4n", "-entry", "all", "-range", "100-"], 0),
    ("D38", &["-db", "db/n4n", "-entry", "all", "-range", "1000-2000"], 0),
    ("D39", &["-db", "db/n4n", "-entry", "all", "-line_length", "1", "-range", "1-3"], 0),
    ("D40", &["-db", "db/n4n", "-entry", "all", "-outfmt", "%a %l %s", "-range", "3-"], 0),
    ("D41", &["-db", "db/p4n", "-entry", "all", "-strand", "minus"], 0),
    ("D42", &["-db", "db/p4n", "-entry", "all", "-range", "5-10", "-strand", "minus"], 0),
    ("D43", &["-db", "db/n4n", "-entry", "all", "-ctrl_a"], 0),
    ("D44", &["-db", "db/n4", "-entry", "all", "-ctrl_a"], 0),
    ("D45", &["-db", "db/n4", "-entry", "all", "-long_seqids"], 0),
    ("D46", &["-db", "db/p4", "-entry", "all", "-long_seqids"], 0),
    ("D48", &["-db", "db/nL", "-entry", "all", "-outfmt", "%a|%i|%g|%t|%T|%e|%l|%h|%P|%X|%N|%C|%S|%L|%B|%K|%n"], 0),
    ("D49", &["-db", "db/p4", "-entry", "all", "-outfmt", "%a|%i|%g|%t|%T|%e|%l|%h|%P|%X|%N|%C|%S|%L|%B|%K|%n"], 0),
    ("D51", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|-1"], 1),
    ("D52", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|abc"], 1),
    ("D53", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|3|"], 0),
    ("D54", &["-db", "db/n4n", "-entry", "GNL|bl_ord_id|3"], 1),
    ("D55", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|03"], 1),
    ("D58", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-range", "0-"], 0),
    ("D59", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-range", "5-"], 0),
    ("D60", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-range", "5-x"], 1),
    ("D62", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-range", "3000000000-5"], 1),
    ("D63", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-range", "5--7"], 0),
    ("D64", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-range", "200-"], 0),
    ("D65", &["-db", "db/p4n", "-entry_batch", "in/crlf_d.txt"], 0),
    ("D66", &["-db", "db/p4n", "-entry_batch", "in/nosuch.txt"], 1),
    ("D67", &["-db", "db/n4n", "-entry", "all", "-outfmt", "%s", "-strand", "minus", "-range", "100-"], 0),
    ("D68", &["-db", "db/p4n", "-entry", "all", "-line_length", "1000"], 0),
    ("D69", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-outfmt", "%a %l", "-range", "5-10"], 0),
    ("D71", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-get_dups"], 0),
    ("D73", &["-db", "db/p4n", "-entry", "all", "-outfmt", "%o %a", "-dbtype", "prot"], 0),
    ("D74", &["-db", "db/p4n", "-entry", "all", "-dbtype", "nucl"], 2),
    ("D75", &["-db", "db/n4n", "-entry", "all", "-outfmt", "%t"], 0),
    ("D76", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-outfmt", "%g"], 0),
    ("D79", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-range", "99999999999999999999-5"], 1),
    ("D80", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-range", "+5-7"], 0),
    ("D81", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|0", "-range", "5-20", "-strand", "minus", "-long_seqids", "-ctrl_a"], 1),
    ("D83", &["-db", "db/n4", "-entry", "all", "-range", "3-9", "-strand", "minus"], 0),
    ("D85", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-outfmt", "%s"], 0),
    ("D86", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-outfmt", "%s", "-range", "5-20", "-strand", "minus"], 0),
    ("D88", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-range", "5-20", "-strand", "minus"], 0),
    ("D89", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-range", "5-20", "-long_seqids"], 0),
    ("D90", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1,gnl|BL_ORD_ID|1", "-long_seqids"], 1),
    ("D91", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1,gnl|BL_ORD_ID|1"], 0),
    ("D92", &["-db", "db/p4n", "-entry_batch", "in/ordid_d2.txt"], 1),
    ("D93", &["-db", "db/p4n", "-entry_batch", "in/ordid_d2.txt", "-get_dups"], 1),
    ("D94", &["-db", "db/n4n", "-entry_batch", "in/ordid_n2.txt"], 1),
    ("D95", &["-db", "db/n4n", "-entry_batch", "in/ordid_n2.txt", "-get_dups"], 1),
    ("D96", &["-db", "db/n5n", "-entry_batch", "in/ordid_n2.txt", "-get_dups"], 1),
    ("D98", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-strand", "minus", "-long_seqids"], 0),
    ("D99", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-target_only", "-outfmt", "%s %a"], 0),
    ("D100", &["-db", "db/p4n", "-entry_batch", "in/ordid_d2.txt", "-target_only"], 1),
    ("D102", &["-db", "db/n4", "-entry", "all", "-long_seqids", "-strand", "minus"], 0),
    ("D103", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-range", "500-", "-strand", "minus"], 0),
    ("D104", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-range", "500-", "-strand", "minus", "-long_seqids"], 0),
    ("D105", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|1", "-long_seqids", "-range", "5-20"], 0),
    ("D106", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-outfmt", "%s %h", "-strand", "minus"], 0),
    ("D107", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-outfmt", "%s", "-range", "5-5"], 0),
    ("D108", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-range", "5-5", "-strand", "minus"], 0),
    ("D109", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|1", "-outfmt", "%a", "-target_only", "-get_dups"], 1),
    ("D110", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-range", "488-"], 0),
    ("D111", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-range", "489-"], 0),
    ("D112", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-range", "488-600", "-strand", "minus"], 0),
];

#[test]
fn retrieval_is_ncbi_s_bytes() {
    let dir = TestDir::new("retrieve");
    for &(id, words, status) in MATCHES {
        check(id, words, status, &dir.0);
    }
}

/// The approved differences (line BDC DECISIONS 12): `-long_seqids` with `-range` on the
/// minus strand, where NCBI writes the forward letters under its `:c` header and LOSAT the
/// reverse complement (the letters NCBI writes for the same range with short ids: D83,
/// D88). The expected bytes are NCBI's header with those letters.
#[rustfmt::skip]
const APPROVED: &[(&str, &[&str], i32)] = &[
    ("D82", &["-db", "db/n4", "-entry", "all", "-long_seqids", "-range", "3-9", "-strand", "minus"], 0),
    ("D87", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-range", "5-20", "-strand", "minus", "-long_seqids"], 0),
    ("D97", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-range", "5-20", "-strand", "minus", "-long_seqids", "-ctrl_a"], 0),
    ("D101", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-long_seqids", "-range", "5-20", "-strand", "minus", "-line_length", "5"], 0),
];

#[test]
fn long_ids_on_the_minus_strand_write_the_reverse_complement() {
    let dir = TestDir::new("approved");
    for &(id, words, status) in APPROVED {
        check(id, words, status, &dir.0);
    }
}

/// `gnl|BL_ORD_ID|n` of a version 4 database in volumes is the one sequence `n` of the
/// whole database (line BDC DECISIONS 12: NCBI writes the sequence `n` of every volume
/// where `n` falls inside the database, probe S01): its FASTA is the `n`-th record of
/// `-entry all`, in every volume.
#[test]
fn ordinal_ids_of_a_database_in_volumes_are_database_wide() {
    let dir = TestDir::new("volumes");
    let (all, stderr, code) = blastdbcmd(&["-db", "../seqdb/db/prot_mv", "-entry", "all"], &dir.0);
    assert_eq!((stderr, code), (vec![], 0));
    let records: Vec<&[u8]> = all
        .split(|&c| c == b'>')
        .filter(|record| !record.is_empty())
        .collect();
    assert!(records.len() > 2, "{} records", records.len());
    for n in [0usize, 1, records.len() / 2, records.len() - 1] {
        let id = format!("gnl|BL_ORD_ID|{n}");
        let (one, stderr, code) =
            blastdbcmd(&["-db", "../seqdb/db/prot_mv", "-entry", &id], &dir.0);
        if n == 0 {
            // NCBI reads the tag `0` as a string and throws (probe D11).
            assert_eq!(code, 1, "{id}");
            assert_eq!(
                String::from_utf8_lossy(&stderr),
                "Error: [blastdbcmd] CObject_id::GetId(): Invalid choice selection: \
                 NCBI-General::Object-id.str\n"
            );
            continue;
        }
        assert_eq!((stderr, code), (vec![], 0), "{id}");
        assert_eq!(
            String::from_utf8_lossy(&one),
            format!(">{}", String::from_utf8_lossy(records[n])),
            "{id}"
        );
    }
}

// The files of `-out` and `-logfile` take what NCBI writes there: D70 the FASTA, D72 the
// diagnostic of the unknown specifier (its partial line stays on stdout), D84 that partial
// line in the file of `-out`.
#[test]
fn out_and_logfile_take_the_retrieval_and_its_diagnostics() {
    let dir = TestDir::new("files");
    let run = |words: &[&str]| blastdbcmd(words, &dir.0);
    let file = |name: &str| fs::read(dir.0.join(name)).unwrap();

    let (stdout, stderr, code) = run(&[
        "-db",
        "db/p4n",
        "-entry",
        "gnl|BL_ORD_ID|3",
        "-out",
        "{O}/D70.out",
    ]);
    assert_eq!((stdout, stderr, code), (vec![], vec![], 0), "D70");
    assert_eq!(file("D70.out"), expected("D70.out"), "D70 -out");

    let (stdout, stderr, code) = run(&[
        "-db",
        "db/p4n",
        "-entry",
        "gnl|BL_ORD_ID|3",
        "-outfmt",
        "%a %Z",
        "-logfile",
        "{O}/D72.logfile",
    ]);
    assert_eq!(
        (String::from_utf8_lossy(&stdout), stderr, code),
        (String::from_utf8_lossy(&expected("D72.stdout")), vec![], 1),
        "D72"
    );
    assert_eq!(file("D72.logfile"), expected("D72.logfile"), "D72 -logfile");

    let (stdout, stderr, code) = run(&[
        "-db",
        "db/p4n",
        "-entry",
        "gnl|BL_ORD_ID|3",
        "-outfmt",
        "%a %Z",
        "-out",
        "{O}/D84.out",
    ]);
    assert_eq!(
        (stdout, String::from_utf8_lossy(&stderr), code),
        (vec![], String::from_utf8_lossy(&expected("D84.stderr")), 1),
        "D84"
    );
    assert_eq!(file("D84.out"), expected("D84.out"), "D84 -out");
}

/// The probe's cases that LOSAT's BLASTDBCMD does not serve: the lookup of an accession, a
/// name or a GI in a database built with `-parse_seqids` (line BDC's BDC2b and BDC2c), an
/// alias of a database with masks and several names in `-db` (DECISIONS 6), `-range ""`
/// (NCBI crashes), and `-outfmt %d` / `%b` (DECISIONS 9).
#[rustfmt::skip]
const REFUSALS: &[(&str, &[&str])] = &[
    ("F09", &["-db", "db/p4", "-entry_batch", "in/pids.txt"]),
    ("F10", &["-db", "db/p5", "-entry_batch", "in/pids.txt"]),
    ("F11", &["-db", "db/p4", "-entry_batch", "in/pids.txt", "-get_dups"]),
    ("F12", &["-db", "db/p5", "-entry_batch", "in/pids.txt", "-get_dups"]),
    ("F17", &["-db", "db/p4", "-entry", "XP_000020.1", "-target_only"]),
    ("F18", &["-db", "db/p4", "-entry", "1819", "-target_only"]),
    ("F19", &["-db", "db/p5", "-entry", "XP_000020.1", "-ctrl_a"]),
    ("F20", &["-db", "db/p4", "-entry", "XP_000001.1,gi|123456,1ABC_A", "-long_seqids"]),
    ("F23", &["-db", "db/n4", "-entry_batch", "in/nids.txt"]),
    ("F32", &["-db", "db/palias", "-entry", "all", "-outfmt", "%o %a %t"]),
    ("F33", &["-db", "db/palias", "-entry", "P10001.1,XP_000001.1"]),
    ("F34", &["-db", "db/p5 db/pm5", "-entry", "all", "-outfmt", "%o %a"]),
    ("S13", &["-db", "db/p4", "-entry", "XP_000001.1,,P00012", "-outfmt", "%o %a"]),
    ("S14", &["-db", "db/p4", "-entry", ""]),
    ("S15", &["-db", "db/p4", "-entry", "emb|CAA00019.1|", "-target_only"]),
    ("S16", &["-db", "db/p4", "-entry", "gi|1819", "-target_only"]),
    ("S17", &["-db", "db/p4", "-entry", "CAA00019.1", "-target_only", "-outfmt", "%a %g"]),
    ("S18", &["-db", "db/n4", "-entry", "NM_000001.1", "-long_seqids", "-range", "5-20", "-strand", "minus"]),
    ("S19", &["-db", "db/n4", "-entry", "NM_000001.1,NM_000001.1", "-long_seqids"]),
    ("S46", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%o", "-get_dups"]),
    ("D61", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-range", ""]),
    ("asn-d", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-outfmt", "%d"]),
    ("asn-b", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-outfmt", "%b"]),
    ("mask-neg", &["-db", "db/p4n", "-entry", "all", "-mask_sequence_with", "-3"]),
];

#[test]
fn what_blastdbcmd_does_not_serve_is_refused_without_output() {
    let dir = TestDir::new("refuse");
    for &(id, words) in REFUSALS {
        let (stdout, stderr, code) = blastdbcmd(words, &dir.0);
        let stderr = String::from_utf8_lossy(&stderr);
        assert_eq!(stdout, b"", "{id} {words:?}");
        assert_eq!(code, 1, "{id} {words:?}: {stderr}");
        assert!(
            stderr.starts_with("Error: ")
                && stderr.ends_with("is not supported by LOSAT's BLASTDBCMD\n")
                && stderr.lines().count() == 1,
            "{id} {words:?}: {stderr}"
        );
    }
}
