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

/// `fs::canonicalize` without the `\\?\` prefix that Windows puts on its result: a verbatim path
/// takes `/` as a character of a name (`{O}/D70.out` below is not a file of the directory then),
/// and the product prints the paths it lists without the prefix (as `blastdbcmd_info.rs`).
fn canonical(path: &Path) -> PathBuf {
    let path = fs::canonicalize(path).unwrap();
    #[cfg(windows)]
    if let Some(plain) = path.to_str().and_then(|text| text.strip_prefix(r"\\?\")) {
        if plain.as_bytes().get(1) == Some(&b':') {
            return PathBuf::from(plain);
        }
    }
    path
}

fn tests_dir() -> PathBuf {
    canonical(&Path::new(env!("CARGO_MANIFEST_DIR")).join("tests"))
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
        TestDir(canonical(&dir))
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
    blastdbcmd_in(words, out_dir, None)
}

/// `blastdbcmd` with `BLASTDB` set to the directory `blastdb` of `tests/blastdbcmd` (where the
/// taxonomy database of the run is), or unset.
fn blastdbcmd_in(words: &[&str], out_dir: &Path, blastdb: Option<&str>) -> (Vec<u8>, Vec<u8>, i32) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text == "BLASTDB" {
            command.env_remove(&name);
        }
    }
    if let Some(dir) = blastdb {
        command.env("BLASTDB", tests_dir().join("blastdbcmd").join(dir));
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

/// The frozen file `name` of `tests/blastdbcmd/expected`, with `{R}` placed, in LOSAT's words
/// (NCBI's bytes through `losat::wording`, AGENTS.md rule 8).
fn expected(name: &str) -> Vec<u8> {
    let tests = tests_dir();
    let text = fs::read(tests.join("blastdbcmd/expected").join(name)).unwrap();
    let text = String::from_utf8(text).unwrap();
    // NCBI's files come from Linux; on Windows the name in "database [name]" has the platform's
    // separator and the search path is split at `;` (`SeqDB_ConvertOSPath`, as `blastdbcmd_info.rs`).
    #[cfg(windows)]
    let text = text
        .replace("database [db/", "database [db\\")
        .replace("{R}::]", "{R};;]");
    let text = text.replace("{R}", tests.join("blastdbcmd").to_str().unwrap());
    losat::wording::losat_output_bytes(text.as_bytes()).into_owned()
}

fn check(id: &str, words: &[&str], status: i32, out_dir: &Path) {
    check_in(id, words, status, out_dir, None);
}

/// `check` with `BLASTDB` set as `blastdbcmd_in` sets it.
fn check_in(id: &str, words: &[&str], status: i32, out_dir: &Path, blastdb: Option<&str>) {
    let (stdout, stderr, code) = blastdbcmd_in(words, out_dir, blastdb);
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
/// databases are the probe's (`db/`): without parsed Seq-ids (`p4n`, `p5n`, `n4n`, `n5n`,
/// `pL`, `nL`), version 4 with them (`p4`, `n4`; `d4` with accessions on two sequences, `i4`
/// with pdb ids of four digits, an EMBL and a GenBank id in one definition line, general ids
/// and versions, `it4` with trace ids, `dv4` in three volumes), and version 5 with them (`p5`,
/// `n5`, `d5` and `e5` with an accession on many sequences).
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
    ("F10", &["-db", "db/p5", "-entry_batch", "in/pids.txt"], 1),
    ("F12", &["-db", "db/p5", "-entry_batch", "in/pids.txt", "-get_dups"], 1),
    ("F19", &["-db", "db/p5", "-entry", "XP_000020.1", "-ctrl_a"], 0),
    ("F24", &["-db", "db/n5", "-entry_batch", "in/nids.txt"], 1),
    ("S03", &["-db", "db/d5", "-entry", "XP_800001", "-outfmt", "%o %a"], 0),
    ("S05", &["-db", "db/d5", "-entry", "XP_800002", "-outfmt", "%o %a"], 0),
    ("S07", &["-db", "db/d5", "-entry_batch", "in/dups.txt", "-outfmt", "%o %a"], 1),
    ("S09", &["-db", "db/d5", "-entry_batch", "in/dups.txt", "-get_dups", "-outfmt", "%o %a"], 1),
    ("S12", &["-db", "db/d5", "-entry", "XP_800001.3", "-outfmt", "%o %a"], 1),
    ("S46", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%o", "-get_dups"], 0),
    ("H01", &["-db", "db/p5", "-entry", ""], 1),
    ("H02", &["-db", "db/p5", "-entry", "XP_000001.1,,P00012", "-outfmt", "%o %a"], 1),
    ("H03", &["-db", "db/p5", "-entry", "emb|CAA00019.1|", "-target_only"], 0),
    ("H04", &["-db", "db/p5", "-entry", "gi|1819", "-target_only"], 0),
    ("H05", &["-db", "db/p5", "-entry", "1819", "-target_only"], 0),
    ("H06", &["-db", "db/p5", "-entry", "CAA00019.1", "-target_only", "-outfmt", "%a %g"], 0),
    ("H07", &["-db", "db/p5", "-entry", "XP_000020.1", "-target_only"], 0),
    ("H08", &["-db", "db/p5", "-entry", "CAA00019.1", "-outfmt", "%o %a %g"], 0),
    ("H11", &["-db", "db/p5", "-entry_batch", "in/pids.txt", "-target_only"], 1),
    ("H12", &["-db", "db/p5", "-entry_batch", "in/pids.txt", "-get_dups", "-outfmt", "%o %a"], 1),
    ("H13", &["-db", "db/p5", "-entry", "0"], 1),
    ("H14", &["-db", "db/p5", "-entry", "-5"], 1),
    ("H15", &["-db", "db/p5", "-entry", "+1819"], 0),
    ("H16", &["-db", "db/p5", "-entry", "99999999999"], 1),
    ("H17", &["-db", "db/p5", "-entry", "99999999999999999999"], 1),
    ("H18", &["-db", "db/p5", "-entry", "1818,1819,1820", "-outfmt", "%o %a"], 0),
    ("H22", &["-db", "db/p5", "-entry", "gnl|BL_ORD_ID|3"], 1),
    ("H23", &["-db", "db/p5", "-entry", "lcl|"], 1),
    ("H24", &["-db", "db/p5", "-entry", "gi|"], 1),
    ("H25", &["-db", "db/p5", "-entry", "|"], 1),
    ("H26", &["-db", "db/p5", "-entry", "_"], 1),
    ("H27", &["-db", "db/p5", "-entry", " XP_000001.1"], 0),
    ("H28", &["-db", "db/p5", "-entry", "XP_000001.1 "], 0),
    ("H29", &["-db", "db/p5", "-entry", "xp_000001.1", "-target_only", "-outfmt", "%a"], 0),
    ("H30", &["-db", "db/p5", "-entry", "XP_000001", "-outfmt", "%o %a"], 0),
    ("H31", &["-db", "db/p5", "-entry", "XP_000001.2", "-outfmt", "%o %a"], 1),
    ("H32", &["-db", "db/n5", "-entry", "654321", "-outfmt", "%o %a %g"], 0),
    ("H33", &["-db", "db/n5", "-entry", "gi|654321", "-target_only", "-outfmt", "%o %a %g"], 0),
    ("H34", &["-db", "db/n5", "-entry", "X00001", "-outfmt", "%o %a"], 0),
    ("H35", &["-db", "db/n5", "-entry_batch", "in/nids.txt", "-get_dups", "-outfmt", "%o %a"], 1),
    ("H36", &["-db", "db/n5", "-entry_batch", "in/nids.txt", "-target_only", "-outfmt", "%o %a"], 1),
    // An alias over two LMDB files with the same keys: `-entry` joins their OIDs, the bulk
    // lookup keeps the last file's.
    ("H120", &["-db", "db/h_dd", "-entry", "XP_800001", "-outfmt", "%o %a"], 0),
    ("H121", &["-db", "db/h_dd", "-entry", "XP_800002", "-outfmt", "%o %a"], 0),
    ("H122", &["-db", "db/h_dd", "-entry_batch", "in/dups.txt", "-outfmt", "%o %a"], 1),
    ("H123", &["-db", "db/h_dd", "-entry_batch", "in/dups.txt", "-get_dups", "-outfmt", "%o %a"], 1),
    // One LMDB file with `acc2oid` and one without: the exception stops `-entry` and `-get_dups`,
    // the bulk lookup swallows it and keeps the first file's answers.
    ("H58", &["-db", "db/h_p5p5n", "-entry", "XP_000001.1", "-outfmt", "%o %a"], 1),
    ("H59", &["-db", "db/h_p5p5n", "-entry_batch", "in/h_p.txt", "-outfmt", "%o %a"], 1),
    ("H60", &["-db", "db/h_p5p5n", "-entry_batch", "in/h_p.txt", "-get_dups", "-outfmt", "%o %a"], 1),
    ("H65", &["-db", "db/h_p5np5", "-entry_batch", "in/h_gi.txt", "-outfmt", "%o %a"], 1),
    // A database of two volumes (`pv`), and its second volume alone: accessions in the other
    // volume are not found, and the OIDs and GIs of the volume start from 0.
    ("H125", &["-db", "db/pv", "-entry", "XP_000001.1,DAA00016.1,1819", "-outfmt", "%o %a"], 0),
    ("H126", &["-db", "db/pv.01", "-entry", "DAA00016.1", "-outfmt", "%o %a"], 0),
    ("H127", &["-db", "db/pv.01", "-entry", "XP_000001.1", "-outfmt", "%o %a"], 1),
    ("H128", &["-db", "db/pv.01", "-entry", "1819", "-outfmt", "%o %a"], 0),
    ("H129", &["-db", "db/pv.01", "-entry_batch", "in/h_gi.txt", "-outfmt", "%o %a"], 1),
    ("H130", &["-db", "db/pv.00", "-entry", "1819", "-outfmt", "%o %a"], 1),
    ("H110", &["-db", "db/e5", "-entry", "XP_900002", "-outfmt", "%o %a"], 0),
    ("H111", &["-db", "db/e5", "-entry", "XP_900001.599", "-outfmt", "%o %a"], 0),
    // Version 4 databases with parsed ids (probe S and I): every spelling of an id in one batch
    // (I895-I907), the version check of an accession, the numbers that overflow, PIGs and trace
    // ids, the target of `-target_only`, keys on several sequences and volumes.
    ("F09", &["-db", "db/p4", "-entry_batch", "in/pids.txt"], 1),
    ("F11", &["-db", "db/p4", "-entry_batch", "in/pids.txt", "-get_dups"], 1),
    ("F17", &["-db", "db/p4", "-entry", "XP_000020.1", "-target_only"], 0),
    ("F18", &["-db", "db/p4", "-entry", "1819", "-target_only"], 0),
    ("F20", &["-db", "db/p4", "-entry", "XP_000001.1,gi|123456,1ABC_A", "-long_seqids"], 0),
    ("F23", &["-db", "db/n4", "-entry_batch", "in/nids.txt"], 1),
    ("S13", &["-db", "db/p4", "-entry", "XP_000001.1,,P00012", "-outfmt", "%o %a"], 1),
    ("S14", &["-db", "db/p4", "-entry", ""], 1),
    ("S15", &["-db", "db/p4", "-entry", "emb|CAA00019.1|", "-target_only"], 0),
    ("S16", &["-db", "db/p4", "-entry", "gi|1819", "-target_only"], 0),
    ("S17", &["-db", "db/p4", "-entry", "CAA00019.1", "-target_only", "-outfmt", "%a %g"], 0),
    ("S19", &["-db", "db/n4", "-entry", "NM_000001.1,NM_000001.1", "-long_seqids"], 1),
    ("S04", &["-db", "db/d4", "-entry", "XP_800001", "-outfmt", "%o %a"], 0),
    ("S06", &["-db", "db/d4", "-entry", "XP_800002", "-outfmt", "%o %a"], 0),
    ("S08", &["-db", "db/d4", "-entry_batch", "in/dups.txt", "-outfmt", "%o %a"], 1),
    ("S10", &["-db", "db/d4", "-entry_batch", "in/dups.txt", "-get_dups", "-outfmt", "%o %a"], 1),
    ("S11", &["-db", "db/d4", "-entry", "XP_800001.3", "-outfmt", "%o %a"], 1),
    ("I106", &["-db", "db/p4", "-entry", "NP_000002.2", "-outfmt", "%o %a"], 1),
    ("I107", &["-db", "db/p4", "-entry", "np_000002.2", "-outfmt", "%o %a"], 1),
    ("I117", &["-db", "db/p4", "-entry", "ref|NP_000002.2|", "-outfmt", "%o %a"], 1),
    ("I172", &["-db", "db/p4", "-entry", "gi|-5", "-outfmt", "%o %a"], 1),
    ("I173", &["-db", "db/p4", "-entry", "gi|99999999999", "-outfmt", "%o %a"], 1),
    ("I174", &["-db", "db/p4", "-entry", "gi|4294967296", "-outfmt", "%o %a"], 1),
    ("I183", &["-db", "db/p4", "-entry", "gnl|PIG|abc", "-outfmt", "%o %a"], 1),
    ("I187", &["-db", "db/p4", "-entry", "gnl|ti|abc", "-outfmt", "%o %a"], 1),
    ("I188", &["-db", "db/p4", "-entry", "gnl|ti|-5", "-outfmt", "%o %a"], 1),
    ("I895", &["-db", "db/p4", "-entry_batch", "in/i_p_all.txt", "-outfmt", "%o %a"], 1),
    ("I896", &["-db", "db/p4", "-entry_batch", "in/i_p_all.txt", "-get_dups", "-outfmt", "%o %a"], 1),
    ("I897", &["-db", "db/p4", "-entry_batch", "in/i_p_all.txt", "-target_only", "-outfmt", "%o %a"], 1),
    ("I898", &["-db", "db/p4", "-entry_batch", "in/i_p_safe.txt", "-outfmt", "%o %a"], 1),
    ("I899", &["-db", "db/p4", "-entry_batch", "in/i_p_safe.txt", "-get_dups", "-outfmt", "%o %a"], 1),
    ("I900", &["-db", "db/p4", "-entry_batch", "in/i_p_safe.txt", "-target_only", "-outfmt", "%o %a"], 1),
    ("I901", &["-db", "db/p4", "-entry_batch", "in/i_p_safe.txt"], 1),
    ("I902", &["-db", "db/p4", "-entry_batch", "in/i_p_safe.txt", "-target_only"], 1),
    ("I299", &["-db", "db/i4", "-entry", "1914", "-outfmt", "%o %a"], 0),
    ("I301", &["-db", "db/i4", "-entry", "pdb|1914|A", "-outfmt", "%o %a"], 0),
    ("I303", &["-db", "db/i4", "-entry", "pdb|1914|b", "-outfmt", "%o %a"], 0),
    ("I308", &["-db", "db/i4", "-entry", "1914A", "-outfmt", "%o %a"], 1),
    ("I313", &["-db", "db/i4", "-entry", "X10001.1", "-outfmt", "%o %a"], 0),
    ("I314", &["-db", "db/i4", "-entry", "U10002.1", "-outfmt", "%o %a"], 0),
    ("I322", &["-db", "db/i4", "-entry", "U30003.1", "-outfmt", "%o %a"], 1),
    ("I336", &["-db", "db/i4", "-entry", "gnl|PIG|77", "-outfmt", "%o %a"], 1),
    ("I337", &["-db", "db/i4", "-entry", "gnl|ti|88", "-outfmt", "%o %a"], 0),
    ("I338", &["-db", "db/i4", "-entry", "88", "-outfmt", "%o %a"], 1),
    ("I347", &["-db", "db/i4", "-entry", "NP_400001", "-outfmt", "%o %a"], 0),
    ("I350", &["-db", "db/i4", "-entry", "NP_400001.3", "-outfmt", "%o %a"], 1),
    ("I356", &["-db", "db/i4", "-entry", "ref|NP_400002.1|", "-outfmt", "%o %a"], 1),
    ("I365", &["-db", "db/i4", "-entry", "1914", "-target_only"], 0),
    ("I367", &["-db", "db/i4", "-entry", "pdb|1914|A", "-target_only"], 0),
    ("I368", &["-db", "db/i4", "-entry", "pdb|1914|a", "-target_only"], 1),
    ("I379", &["-db", "db/i4", "-entry", "X10001.1", "-target_only"], 0),
    ("I380", &["-db", "db/i4", "-entry", "U10002.1", "-target_only"], 1),
    ("I382", &["-db", "db/i4", "-entry", "gb|U10002.1|", "-target_only"], 1),
    ("I384", &["-db", "db/i4", "-entry", "X20002.1", "-target_only"], 1),
    ("I405", &["-db", "db/i4", "-entry", "gnl|BL_ORD_ID|2", "-target_only"], 1),
    ("I434", &["-db", "db/i4", "-entry", "pdb|1914|a", "-target_only", "-outfmt", "%o %a"], 1),
    ("I446", &["-db", "db/i4", "-entry", "U10002.1", "-target_only", "-outfmt", "%o %a"], 0),
    ("I471", &["-db", "db/i4", "-entry", "gnl|BL_ORD_ID|2", "-target_only", "-outfmt", "%o %a"], 1),
    ("I903", &["-db", "db/i4", "-entry_batch", "in/i_i_all.txt", "-outfmt", "%o %a"], 1),
    ("I904", &["-db", "db/i4", "-entry_batch", "in/i_i_all.txt", "-get_dups", "-outfmt", "%o %a"], 1),
    ("I905", &["-db", "db/i4", "-entry_batch", "in/i_i_all.txt", "-target_only", "-outfmt", "%o %a"], 1),
    ("I906", &["-db", "db/i4", "-entry_batch", "in/i_i_range.txt"], 1),
    ("I907", &["-db", "db/i4", "-entry_batch", "in/i_i_range.txt", "-get_dups"], 1),
    ("I695", &["-db", "db/it4", "-entry", "gnl|ti|12345", "-outfmt", "%o %a"], 0),
    ("I697", &["-db", "db/it4", "-entry", "gnl|ti|99", "-outfmt", "%o %a"], 1),
    ("I698", &["-db", "db/it4", "-entry", "12345", "-outfmt", "%o %a"], 1),
    ("I699", &["-db", "db/it4", "-entry", "gnl|ti|abc", "-outfmt", "%o %a"], 1),
    ("I704", &["-db", "db/it4", "-entry", "gnl|SRA|SRR1.1", "-outfmt", "%o %a"], 0),
    ("I705", &["-db", "db/it4", "-entry", "SRA:SRR1.1", "-outfmt", "%o %a"], 0),
    ("I707", &["-db", "db/it4", "-entry", "SRR1.1", "-outfmt", "%o %a"], 0),
    ("I709", &["-db", "db/it4", "-entry", "TSA:T1", "-outfmt", "%o %a"], 0),
    ("I711", &["-db", "db/n4", "-entry", "NM_000001", "-outfmt", "%o %a"], 0),
    ("I714", &["-db", "db/n4", "-entry", "NM_000001.2", "-outfmt", "%o %a"], 1),
    ("I724", &["-db", "db/n4", "-entry", "mydb:nuc7", "-outfmt", "%o %a"], 1),
    ("I734", &["-db", "db/n4", "-entry", "pat|US|1234567|14", "-outfmt", "%o %a"], 0),
    ("I741", &["-db", "db/n4", "-entry", "NM_000001", "-target_only"], 0),
    ("I761", &["-db", "db/n4", "-entry", "LOCUS13", "-target_only"], 1),
    ("I771", &["-db", "db/n4", "-entry", "NM_000001.1", "-range", "5-20", "-strand", "minus"], 0),
    ("I772", &["-db", "db/n4", "-entry", "NM_000001.1", "-range", "5-20", "-long_seqids"], 0),
    ("I773", &["-db", "db/n4", "-entry", "NM_000001.1", "-outfmt", "%s %h %l", "-strand", "minus", "-range", "3-9"], 0),
    ("I774", &["-db", "db/n4", "-entry", "gi|654321", "-outfmt", "%s %h %l", "-strand", "minus", "-range", "3-9"], 0),
    ("I781", &["-db", "db/d4", "-entry", "xp_800001", "-outfmt", "%o %a"], 0),
    ("I785", &["-db", "db/d4", "-entry", "XP_800001.256", "-outfmt", "%o %a"], 1),
    ("I786", &["-db", "db/d4", "-entry", "XP_800001", "-target_only", "-outfmt", "%o %a"], 0),
    ("I788", &["-db", "db/d4", "-entry", "XP_800001,XP_800001", "-outfmt", "%o %a"], 0),
    ("I789", &["-db", "db/d4", "-entry", "XP_800001,XP_800001.1", "-outfmt", "%o %a"], 0),
    ("I790", &["-db", "db/d4", "-entry_batch", "in/dups.txt", "-target_only", "-outfmt", "%o %a"], 1),
    ("I791", &["-db", "db/dv4", "-entry", "XP_800001", "-outfmt", "%o %a"], 0),
    ("I792", &["-db", "db/dv4", "-entry", "XP_800002", "-outfmt", "%o %a"], 0),
    ("I794", &["-db", "db/dv4", "-entry", "XP_800001.3", "-outfmt", "%o %a"], 1),
    ("I796", &["-db", "db/dv4", "-entry_batch", "in/dups.txt", "-outfmt", "%o %a"], 1),
    ("I797", &["-db", "db/dv4", "-entry_batch", "in/dups.txt", "-get_dups", "-outfmt", "%o %a"], 1),
    ("I801", &["-db", "db/dv4.00", "-entry", "XP_800001.3", "-outfmt", "%o %a"], 1),
    ("I805", &["-db", "db/dv4.01", "-entry", "XP_800001", "-outfmt", "%o %a"], 1),
    ("I810", &["-db", "db/dv4.01", "-entry_batch", "in/dups.txt", "-outfmt", "%o %a"], 1),
    ("I831", &["-db", "db/i_dvr", "-entry_batch", "in/dups.txt", "-outfmt", "%o %a"], 1),
    ("I832", &["-db", "db/i_dvr", "-entry_batch", "in/dups.txt", "-get_dups", "-outfmt", "%o %a"], 1),
    ("I838", &["-db", "db/i_dv1", "-entry_batch", "in/dups.txt", "-outfmt", "%o %a"], 1),
    ("I775", &["-db", "db/p4n", "-entry", "1", "-outfmt", "%o %a"], 1),
    ("I776", &["-db", "db/p4n", "-entry", "lcl|3", "-outfmt", "%o %a"], 1),
    ("I777", &["-db", "db/p4n", "-entry", "gi|5", "-outfmt", "%o %a"], 1),
    // Hexadecimal values of the keys with a constraint (`-line_length` `>=1`, `-ipg` `>=0`): NCBI's
    // reader converts them and its constraint reads their text (line RG item K4i-1, 2026-11-18;
    // G11 in `blastdbcmd_info.rs`).
    ("G01", &["-db", "db/n4n", "-entry", "all", "-line_length", "0x5"], 0),
    ("G02", &["-db", "db/p5n", "-entry", "all", "-line_length", "0X10"], 0),
    ("G03", &["-db", "db/n4n", "-entry", "all", "-line_length", "0x"], 1),
    ("G04", &["-db", "db/n4n", "-entry", "all", "-line_length", "0x0"], 1),
    ("G05", &["-db", "db/n4n", "-entry", "all", "-line_length", "0x1p3"], 1),
    ("G06", &["-db", "db/n4n", "-entry", "all", "-line_length", "0x80000000"], 1),
    ("G08", &["-db", "db/n4n", "-entry", "all", "-line_length", "0x5g"], 1),
    ("G09", &["-db", "db/ak/m_r", "-ipg", "0x3E9"], 0),
    ("G10", &["-db", "db/ak/m_r", "-ipg", "0x"], 1),
];

#[test]
fn retrieval_is_ncbi_s_bytes() {
    let dir = TestDir::new("retrieve");
    for &(id, words, status) in MATCHES {
        // D65's batch file ends its lines with `\r\n`, which NCBI reads apart on Linux alone
        // (`NcbiGetlineEOL`; the bytes of Windows and macOS: probe O27,
        // `blastdbcmd_batch_line_ends.rs`).
        if id == "D65" && (cfg!(windows) || cfg!(target_os = "macos")) {
            continue;
        }
        check(id, words, status, &dir.0);
    }
}

/// The probe's cases of `-outfmt %d` (the definition lines as ASN.1 text: `Q101`-`Q144`), `%b`
/// (the Bioseq as ASN.1 text: `Q201`-`Q244`) and both (`Q301`-`Q331`) that read no taxonomy
/// database (line RG, item Dbc): id, words, NCBI's exit status. `{O}/Q126.out` and `{O}/Q238.out`
/// are files of `-out`, whose expected bytes are the files of the same names. The databases are
/// those of `MATCHES`, the volumes of `j/lst` with `-taxids` and `-no_taxid_expansion`, `db/m/pmg`
/// (PIGs), several names in `-db` (`db/n5 db/n5n`, `db/p5 db/p5n`: the ordinal ids of the second
/// database's definition lines are database-wide, so the header in `%b` is written again) and the
/// Seq-id module's databases (`../seqid/db`: every Seq-id form that makeblastdb writes), `db/pl4`
/// (definition lines with memberships and other-info), and `db/pe4`, `db/pz4` and `db/pz5/pl5`:
/// `db/pl4` and `db/pl5` with a header made by hand, which makeblastdb never writes (a set of no
/// lines at OID 1; no bytes at OID 1, or at OID 0 of a database with taxids): NCBI's reader tells
/// a set of no lines from no set. `db/lc/lc4` and `db/lc/lc5` (makeblastdb `-parse_seqids`, version 4
/// and 5) hold local ids, read with `-target_only` in each format: `%b` reads the target with
/// `fParse_AnyRaw` (a bare local id is "Malformatted ID", a number is a GI that matches the local
/// id of its label), the other formats with the lookup's flags. `db/pzm`, `db/pzm2` and `db/pzm0` are
/// alias files over `db/pz4` with `OID_MASK_TYPE` 1, 2 and 0: NCBI's header read filters the set
/// under any type other than 0, and the filter of no set throws "Attempt to access NULL pointer.".
#[rustfmt::skip]
const ASN_TEXTS: &[(&str, &[&str], i32)] = &[
    ("Q101", &["-db", "db/p5", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q102", &["-db", "db/p4", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q103", &["-db", "db/n5", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q104", &["-db", "db/n4", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q105", &["-db", "db/p5n", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q106", &["-db", "db/n4n", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q107", &["-db", "db/i4", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q108", &["-db", "db/dv4", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q109", &["-db", "db/pv", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q110", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-outfmt", "%d"], 0),
    ("Q111", &["-db", "db/p5", "-entry", "AAA00018.1", "-outfmt", "%d"], 0),
    ("Q112", &["-db", "db/p5", "-entry", "CAA00019.1", "-target_only", "-outfmt", "%d"], 0),
    ("Q113", &["-db", "db/p5", "-entry", "AAA00018.1", "-target_only", "-outfmt", "%d"], 0),
    ("Q114", &["-db", "db/p5", "-entry", "all", "-outfmt", "%a %d"], 0),
    ("Q115", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "x%dy"], 0),
    ("Q116", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%d%d"], 0),
    ("Q117", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%%d"], 0),
    ("Q118", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%f %d"], 1),
    ("Q119", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%d %b"], 1),
    ("Q120", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%d %m11"], 0),
    ("Q121", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%d", "-range", "2-5"], 0),
    ("Q122", &["-db", "db/n5", "-entry", "all", "-outfmt", "%d", "-strand", "minus", "-range", "3-9"], 0),
    ("Q123", &["-db", "db/p5", "-entry_batch", "in/pids.txt", "-outfmt", "%d"], 1),
    ("Q124", &["-db", "db/p5", "-entry_batch", "in/pids.txt", "-outfmt", "%d", "-get_dups"], 1),
    ("Q126", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%d", "-out", "{O}/Q126.out"], 0),
    ("Q127", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%d"], 0),
    ("Q128", &["-db", "db/p5", "-entry", "AAA00018.1", "-outfmt", "%d", "-ctrl_a"], 0),
    ("Q129", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%d", "-mask_sequence_with", "11"], 0),
    ("Q132", &["-db", "db/p5", "-entry", "AAA00018.1", "-outfmt", "%d", "-long_seqids"], 0),
    ("Q133", &["-db", "db/p5", "-entry", "nosuch", "-outfmt", "%d"], 1),
    ("Q134", &["-db", "db/h_p5p5n", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q135", &["-db", "db/pL", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q136", &["-db", "db/nL", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q137", &["-db", "db/d5", "-entry", "XP_800001", "-outfmt", "%d"], 0),
    ("Q138", &["-db", "db/d5", "-entry", "XP_800001", "-target_only", "-outfmt", "%d"], 0),
    ("Q139", &["-db", "db/i4", "-entry_batch", "in/i_i_all.txt", "-target_only", "-outfmt", "%d"], 1),
    ("Q140", &["-db", "db/p5", "-entry", "1819", "-target_only", "-outfmt", "%d"], 0),
    ("Q141", &["-db", "db/it4", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q142", &["-db", "db/pm5", "-entry", "all", "-outfmt", "%d %m99"], 0),
    ("Q143", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%d %Z"], 0),
    ("Q144", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%Z %d"], 0),
    ("Q201", &["-db", "db/p5", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q202", &["-db", "db/p4", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q203", &["-db", "db/n5", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q204", &["-db", "db/n4", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q205", &["-db", "db/p5n", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q206", &["-db", "db/n4n", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q207", &["-db", "db/n5n", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q208", &["-db", "db/i4", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q210", &["-db", "db/pv", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q211", &["-db", "db/it4", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q212", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q213", &["-db", "db/nm4", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q214", &["-db", "db/pL", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q215", &["-db", "db/nL", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q216", &["-db", "db/pm5", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q217", &["-db", "db/p5", "-entry", "AAA00018.1", "-target_only", "-outfmt", "%b"], 0),
    ("Q218", &["-db", "db/p5", "-entry", "CAA00019.1", "-target_only", "-outfmt", "%b"], 0),
    ("Q219", &["-db", "db/p5", "-entry", "1819", "-target_only", "-outfmt", "%b"], 0),
    ("Q220", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%s %b"], 0),
    ("Q221", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%%b"], 0),
    ("Q222", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "x%by"], 0),
    ("Q223", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%b", "-range", "2-5"], 0),
    ("Q224", &["-db", "db/n5", "-entry", "all", "-outfmt", "%b", "-strand", "minus"], 0),
    ("Q225", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%b", "-mask_sequence_with", "11"], 0),
    ("Q226", &["-db", "db/p5", "-entry_batch", "in/pids.txt", "-outfmt", "%b", "-get_dups"], 1),
    ("Q228", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%b"], 0),
    ("Q229", &["-db", "db/p5", "-entry", "AAA00018.1", "-outfmt", "%b", "-ctrl_a"], 0),
    ("Q230", &["-db", "db/p5", "-entry", "AAA00018.1", "-outfmt", "%b", "-long_seqids"], 0),
    ("Q233", &["-db", "../db_search/t5/pmix", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q236", &["-db", "db/d5", "-entry", "XP_800001", "-outfmt", "%b"], 0),
    ("Q237", &["-db", "db/d5", "-entry", "XP_800001", "-target_only", "-outfmt", "%b"], 0),
    ("Q238", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%b", "-out", "{O}/Q238.out"], 0),
    ("Q239", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%f %b"], 1),
    ("Q240", &["-db", "db/h_p5p5n", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q241", &["-db", "db/i4", "-entry_batch", "in/i_i_all.txt", "-target_only", "-outfmt", "%b"], 1),
    ("Q242", &["-db", "db/p5", "-entry", "XP_000001.1", "-outfmt", "%b %m11"], 0),
    ("Q243", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-target_only", "-outfmt", "%b"], 0),
    ("Q244", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-target_only", "-outfmt", "%d"], 0),
    ("Q301", &["-db", "j/lst/p5", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%b"], 0),
    ("Q302", &["-db", "j/lst/p5", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%d"], 0),
    ("Q303", &["-db", "j/lst/p5 j/lst/p5n", "-taxids", "9606", "-no_taxid_expansion", "-target_only", "-outfmt", "%b"], 0),
    ("Q304", &["-db", "j/lst/p5 j/lst/p5n", "-taxids", "9606", "-no_taxid_expansion", "-target_only", "-outfmt", "%d"], 0),
    ("Q308", &["-db", "db/ak/or/A_lo3g_a", "-dbtype", "prot", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%b"], 0),
    ("Q309", &["-db", "db/d5", "-entry_batch", "in/dups.txt", "-get_dups", "-outfmt", "%b"], 1),
    ("Q310", &["-db", "db/d5", "-entry_batch", "in/dups.txt", "-get_dups", "-target_only", "-outfmt", "%d"], 1),
    ("Q311", &["-db", "db/i4", "-entry_batch", "in/i_p_safe.txt", "-target_only", "-outfmt", "%b"], 1),
    ("Q314", &["-db", "../seqid/db/nucl_ids", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q315", &["-db", "../seqid/db/nucl_ids", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q316", &["-db", "../seqid/db/prot_extra", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q317", &["-db", "../seqid/db/prot_extra", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q318", &["-db", "../seqid/db/variants", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q319", &["-db", "../seqid/db/variants", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q320", &["-db", "db/n5 db/n5n", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q321", &["-db", "db/p5 db/p5n", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q322", &["-db", "db/m/pmg", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q326", &["-db", "db/p5", "-entry", "all", "-outfmt", "%b", "-line_length", "10"], 0),
    ("Q327", &["-db", "db/m/pmg", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q328", &["-db", "db/n4n db/n5", "-entry", "all", "-outfmt", "%b"], 2),
    ("Q330", &["-db", "db/pl4", "-dbtype", "prot", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q331", &["-db", "db/pl4", "-dbtype", "prot", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q332", &["-db", "db/pe4", "-dbtype", "prot", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q333", &["-db", "db/pe4", "-dbtype", "prot", "-entry", "5003", "-outfmt", "%d"], 0),
    ("Q334", &["-db", "db/pe4", "-dbtype", "prot", "-entry", "5003", "-target_only", "-outfmt", "%d"], 1),
    ("Q335", &["-db", "db/pe4", "-dbtype", "prot", "-entry", "5003", "-target_only", "-outfmt", "%a"], 1),
    ("Q336", &["-db", "db/pe4", "-dbtype", "prot", "-entry", "5003", "-target_only", "-outfmt", "%b"], 1),
    ("Q337", &["-db", "db/pz4", "-dbtype", "prot", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q338", &["-db", "db/pz4", "-dbtype", "prot", "-entry", "5003", "-target_only", "-outfmt", "%b"], 1),
    ("Q339", &["-db", "db/pz4", "-dbtype", "prot", "-entry", "5003", "-target_only"], 1),
    ("Q340", &["-db", "db/pz4", "-dbtype", "prot", "-entry", "5003", "-target_only", "-outfmt", "%a"], 0),
    ("Q341", &["-db", "db/pz5/pl5", "-dbtype", "prot", "-taxids", "9606", "-target_only", "-outfmt", "%d"], 1),
    ("Q342", &["-db", "db/pz5/pl5", "-dbtype", "prot", "-taxids", "9606", "-target_only", "-outfmt", "%a"], 1),
    ("Q343", &["-db", "db/pz5/pl5", "-dbtype", "prot", "-taxids", "9606", "-outfmt", "%d"], 0),
    ("Q344", &["-db", "db/lc/lc4", "-entry", "myprot9", "-target_only", "-outfmt", "%b"], 1),
    ("Q345", &["-db", "db/lc/lc4", "-entry", "myprot9", "-target_only", "-outfmt", "%d"], 0),
    ("Q346", &["-db", "db/lc/lc4", "-entry", "myprot9", "-target_only", "-outfmt", "%f"], 0),
    ("Q347", &["-db", "db/lc/lc4", "-entry", "myprot9", "-target_only", "-outfmt", "%a"], 0),
    ("Q348", &["-db", "db/lc/lc4", "-entry", "lcl|myprot9", "-target_only", "-outfmt", "%b"], 0),
    ("Q349", &["-db", "db/lc/lc4", "-entry", "lcl|myprot9", "-target_only", "-outfmt", "%d"], 0),
    ("Q350", &["-db", "db/lc/lc4", "-entry", "lcl|myprot9", "-target_only", "-outfmt", "%f"], 0),
    ("Q351", &["-db", "db/lc/lc4", "-entry", "lcl|myprot9", "-target_only", "-outfmt", "%a"], 0),
    ("Q352", &["-db", "db/lc/lc4", "-entry", "myprot13", "-target_only", "-outfmt", "%b"], 1),
    ("Q353", &["-db", "db/lc/lc4", "-entry", "myprot13", "-target_only", "-outfmt", "%d"], 0),
    ("Q354", &["-db", "db/lc/lc4", "-entry", "myprot13", "-target_only", "-outfmt", "%f"], 0),
    ("Q355", &["-db", "db/lc/lc4", "-entry", "myprot13", "-target_only", "-outfmt", "%a"], 0),
    ("Q356", &["-db", "db/lc/lc4", "-entry", "12345", "-target_only", "-outfmt", "%b"], 0),
    ("Q357", &["-db", "db/lc/lc4", "-entry", "12345", "-target_only", "-outfmt", "%d"], 0),
    ("Q358", &["-db", "db/lc/lc4", "-entry", "12345", "-target_only", "-outfmt", "%f"], 0),
    ("Q359", &["-db", "db/lc/lc4", "-entry", "12345", "-target_only", "-outfmt", "%a"], 0),
    ("Q360", &["-db", "db/lc/lc4", "-entry", "lcl|12345", "-target_only", "-outfmt", "%b"], 0),
    ("Q361", &["-db", "db/lc/lc4", "-entry", "lcl|12345", "-target_only", "-outfmt", "%d"], 0),
    ("Q362", &["-db", "db/lc/lc4", "-entry", "lcl|12345", "-target_only", "-outfmt", "%f"], 0),
    ("Q363", &["-db", "db/lc/lc4", "-entry", "lcl|12345", "-target_only", "-outfmt", "%a"], 0),
    ("Q364", &["-db", "db/lc/lc5", "-entry", "myprot9", "-target_only", "-outfmt", "%b"], 1),
    ("Q365", &["-db", "db/lc/lc5", "-entry", "myprot9", "-target_only", "-outfmt", "%d"], 0),
    ("Q366", &["-db", "db/lc/lc5", "-entry", "myprot9", "-target_only", "-outfmt", "%f"], 0),
    ("Q367", &["-db", "db/lc/lc5", "-entry", "myprot9", "-target_only", "-outfmt", "%a"], 0),
    ("Q368", &["-db", "db/lc/lc5", "-entry", "lcl|myprot9", "-target_only", "-outfmt", "%b"], 0),
    ("Q369", &["-db", "db/lc/lc5", "-entry", "lcl|myprot9", "-target_only", "-outfmt", "%d"], 0),
    ("Q370", &["-db", "db/lc/lc5", "-entry", "lcl|myprot9", "-target_only", "-outfmt", "%f"], 0),
    ("Q371", &["-db", "db/lc/lc5", "-entry", "lcl|myprot9", "-target_only", "-outfmt", "%a"], 0),
    ("Q372", &["-db", "db/lc/lc5", "-entry", "myprot13", "-target_only", "-outfmt", "%b"], 1),
    ("Q373", &["-db", "db/lc/lc5", "-entry", "myprot13", "-target_only", "-outfmt", "%d"], 0),
    ("Q374", &["-db", "db/lc/lc5", "-entry", "myprot13", "-target_only", "-outfmt", "%f"], 0),
    ("Q375", &["-db", "db/lc/lc5", "-entry", "myprot13", "-target_only", "-outfmt", "%a"], 0),
    ("Q376", &["-db", "db/lc/lc5", "-entry", "12345", "-target_only", "-outfmt", "%b"], 0),
    ("Q377", &["-db", "db/lc/lc5", "-entry", "12345", "-target_only", "-outfmt", "%d"], 0),
    ("Q378", &["-db", "db/lc/lc5", "-entry", "12345", "-target_only", "-outfmt", "%f"], 0),
    ("Q379", &["-db", "db/lc/lc5", "-entry", "12345", "-target_only", "-outfmt", "%a"], 0),
    ("Q380", &["-db", "db/lc/lc5", "-entry", "lcl|12345", "-target_only", "-outfmt", "%b"], 0),
    ("Q381", &["-db", "db/lc/lc5", "-entry", "lcl|12345", "-target_only", "-outfmt", "%d"], 0),
    ("Q382", &["-db", "db/lc/lc5", "-entry", "lcl|12345", "-target_only", "-outfmt", "%f"], 0),
    ("Q383", &["-db", "db/lc/lc5", "-entry", "lcl|12345", "-target_only", "-outfmt", "%a"], 0),
    ("Q384", &["-db", "db/pzm", "-dbtype", "prot", "-entry", "all", "-outfmt", "%a"], 1),
    ("Q385", &["-db", "db/pzm", "-dbtype", "prot", "-entry", "all", "-outfmt", "%d"], 1),
    ("Q386", &["-db", "db/pzm", "-dbtype", "prot", "-entry", "all"], 1),
    ("Q387", &["-db", "db/pzm", "-dbtype", "prot", "-entry", "all", "-outfmt", "%b"], 1),
    ("Q388", &["-db", "db/pzm", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o"], 0),
    ("Q389", &["-db", "db/pzm", "-dbtype", "prot", "-entry", "XP_5001.1", "-outfmt", "%d"], 1),
    ("Q390", &["-db", "db/pzm", "-dbtype", "prot", "-entry", "XP_5001.1", "-target_only", "-outfmt", "%a"], 1),
    ("Q391", &["-db", "db/pzm2", "-dbtype", "prot", "-entry", "all", "-outfmt", "%a"], 1),
    ("Q392", &["-db", "db/pzm2", "-dbtype", "prot", "-entry", "all", "-outfmt", "%d"], 1),
    ("Q393", &["-db", "db/pzm2", "-dbtype", "prot", "-entry", "all", "-outfmt", "%b"], 1),
    ("Q394", &["-db", "db/pzm2", "-dbtype", "prot", "-entry", "all"], 1),
    ("Q395", &["-db", "db/pzm0", "-dbtype", "prot", "-entry", "all", "-outfmt", "%a"], 0),
    ("Q396", &["-db", "db/pzm0", "-dbtype", "prot", "-entry", "all", "-outfmt", "%d"], 0),
];

/// The cases of `%d` and `%b` that read a taxonomy database (`BLASTDB`, a directory of
/// `tests/blastdbcmd`): id, `BLASTDB`, words, NCBI's exit status. `Q235`'s names hold bytes
/// outside ' ' ..= '~', which NCBI's text stream writes as `#` and reports on the standard error.
#[rustfmt::skip]
const ASN_TEXTS_TAXDB: &[(&str, &str, &[&str], i32)] = &[
    ("Q130", "../seqdb/taxdb/t5", &["-db", "../db_search/t5/nmix", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q131", "../seqdb/taxdb/t5", &["-db", "../db_search/t5/pmix", "-entry", "all", "-outfmt", "%d"], 0),
    ("Q231", "../seqdb/taxdb/t5", &["-db", "../db_search/t5/nmix", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q232", "../seqdb/taxdb/t5", &["-db", "../db_search/t5/pmix", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q234", "../seqdb/taxdb/t5", &["-db", "../db_search/t5/nmix", "-entry", "AB311001.1", "-target_only", "-outfmt", "%b"], 0),
    ("Q235", "../seqdb/taxdb/k5", &["-db", "../db_search/t5/pk5", "-entry", "all", "-outfmt", "%b"], 0),
    ("Q305", "j/jtax5", &["-db", "j/lst/p5", "-taxids", "9604", "-outfmt", "%b"], 0),
    ("Q306", "j/jtax5", &["-db", "j/lst/p5", "-taxids", "10066", "-target_only", "-outfmt", "%b"], 0),
    ("Q307", "j/jtax5", &["-db", "j/lst/p5", "-taxids", "10066", "-target_only", "-outfmt", "%d"], 0),
];

#[test]
fn asn_texts_are_ncbi_s_bytes() {
    let dir = TestDir::new("asn");
    for &(id, words, status) in ASN_TEXTS {
        check(id, words, status, &dir.0);
        if let Some(out) = words.iter().find_map(|word| word.strip_prefix("{O}/")) {
            assert_eq!(
                fs::read(dir.0.join(out)).unwrap(),
                expected(out),
                "{id} -out"
            );
        }
    }
    for &(id, blastdb, words, status) in ASN_TEXTS_TAXDB {
        check_in(id, words, status, &dir.0, Some(blastdb));
    }
}

/// One case of NCBI's unit tests of `CSeqDB` on its databases `seqn` and `seqp` (copies in
/// `tests/seqdb/db`): blastdbcmd's `-outfmt fmt` of the sequence `entry` is `expected`, the text
/// that NCBI's case expects of `GetBioseq` or `GetHdr` (`MSerial_AsnText`).
fn ncbi_unit_case(db: &str, entry: &str, fmt: &str, expected: &str) {
    let dir = TestDir::new("unit");
    let db = format!("../seqdb/db/{db}");
    let words = ["-db", db.as_str(), "-entry", entry, "-outfmt", fmt];
    let (stdout, stderr, code) = blastdbcmd(&words, &dir.0);
    assert_eq!(String::from_utf8_lossy(&stdout), expected, "{words:?}");
    assert_eq!(
        (String::from_utf8_lossy(&stderr), code),
        ("".into(), 0),
        "{words:?}"
    );
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/seqdb_reader/seqdb_unit_test.cpp:534-575 GetBioseqN
/// `GetBioseq` of OID 2 on `seqn`: `ncbi4na` data (the sequence has an ambiguous base).
#[test]
fn get_bioseq_of_ncbi_s_seqn_is_its_unit_test_s_text() {
    ncbi_unit_case(
        "seqn",
        "46071107",
        "%b",
        r#"Bioseq ::= {
  id {
    gi 46071107,
    ddbj {
      accession "BP722514",
      version 1
    }
  },
  descr {
    title "Xenopus laevis NBRP cDNA clone:XL452f07ex, 3' end",
    user {
      type str "ASN1_BlastDefLine",
      data {
        {
          label str "ASN1_BlastDefLine",
          num 1,
          data oss {
            '30803080A0801A3158656E6F707573206C6165766973204E4252502063444E4120
636C6F6E653A584C34353266303765782C20332720656E640000A1803080AB80020402BEFD4300
00AC803080A1801A0842503732323531340000A38002010100000000000000000000A280020100
000000000000'H
          }
        }
      }
    }
  },
  inst {
    repr raw,
    mol na,
    length 165,
    seq-data ncbi4na '11428288218841844814141422818811214421121482118428221114
82211121141881228484211141128842148481121112222F882124422141148188842112118488
41114822882844214144144148281181'H
  }
}
"#,
    );
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/seqdb_reader/seqdb_unit_test.cpp:577-617 GetBioseqP
/// `GetBioseq` of OID 2 on `seqp`: `ncbistdaa` data.
#[test]
fn get_bioseq_of_ncbi_s_seqp_is_its_unit_test_s_text() {
    ncbi_unit_case(
        "seqp",
        "44268131",
        "%b",
        r#"Bioseq ::= {
  id {
    gi 44268131,
    genbank {
      accession "EAI08555",
      version 1
    }
  },
  descr {
    title "unknown [environmental sequence]",
    user {
      type str "ASN1_BlastDefLine",
      data {
        {
          label str "ASN1_BlastDefLine",
          num 1,
          data oss {
            '30803080A0801A20756E6B6E6F776E205B656E7669726F6E6D656E74616C207365
7175656E63655D0000A1803080AB80020402A37A630000A4803080A1801A084541493038353535
0000A38002010100000000000000000000A280020100000000000000'H
          }
        }
      }
    }
  },
  inst {
    repr raw,
    mol aa,
    length 64,
    seq-data ncbistdaa '0C0A0A0606090B0909060909060B09060909131004160F090A0A0A
0A0B0A0B0D0D010B0D110D0606090F12090D0A0B0904050D0A160D0B0B05051009100B1005'H
  }
}
"#,
    );
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/seqdb_reader/seqdb_unit_test.cpp:619-639 GetHdrN
/// `GetHdr` of OID 0 on `seqn`.
#[test]
fn get_hdr_of_ncbi_s_seqn_is_its_unit_test_s_text() {
    ncbi_unit_case(
        "seqn",
        "46071105",
        "%d",
        r#"Blast-def-line-set ::= {
  {
    title "Xenopus laevis NBRP cDNA clone:XL452f05ex, 3' end",
    seqid {
      gi 46071105,
      ddbj {
        accession "BP722512",
        version 1
      }
    },
    taxid 0
  }
}
"#,
    );
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/seqdb_reader/seqdb_unit_test.cpp:641-661 GetHdrP
/// `GetHdr` of OID 0 on `seqp`.
#[test]
fn get_hdr_of_ncbi_s_seqp_is_its_unit_test_s_text() {
    ncbi_unit_case(
        "seqp",
        "38083732",
        "%d",
        r#"Blast-def-line-set ::= {
  {
    title "similar to KIAA0960 protein [Mus musculus]",
    seqid {
      gi 38083732,
      other {
        accession "XP_357594",
        version 1
      }
    },
    taxid 0
  }
}
"#,
    );
}

/// The approved differences (line BDC DECISIONS 12): `-long_seqids` with `-range` on the
/// minus strand, where NCBI writes the forward letters under its `:c` header and LOSAT the
/// reverse complement (the letters NCBI writes for the same range with short ids: D83,
/// D88, I771). The expected bytes are NCBI's header with those letters. And
/// `gnl|BL_ORD_ID|n` in a version 4 database of several volumes (I795, I830): NCBI adds `n` to
/// the start of every volume, LOSAT writes the one sequence `n` (NCBI's first record).
#[rustfmt::skip]
const APPROVED: &[(&str, &[&str], i32)] = &[
    ("D82", &["-db", "db/n4", "-entry", "all", "-long_seqids", "-range", "3-9", "-strand", "minus"], 0),
    ("D87", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-range", "5-20", "-strand", "minus", "-long_seqids"], 0),
    ("D97", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-range", "5-20", "-strand", "minus", "-long_seqids", "-ctrl_a"], 0),
    ("D101", &["-db", "db/n4n", "-entry", "gnl|BL_ORD_ID|1", "-long_seqids", "-range", "5-20", "-strand", "minus", "-line_length", "5"], 0),
    ("H09", &["-db", "db/n5", "-entry", "NM_000001.1", "-long_seqids", "-range", "5-20", "-strand", "minus"], 0),
    ("S18", &["-db", "db/n4", "-entry", "NM_000001.1", "-long_seqids", "-range", "5-20", "-strand", "minus"], 0),
    ("I795", &["-db", "db/dv4", "-entry", "gnl|BL_ORD_ID|2", "-outfmt", "%o %a"], 0),
    ("I830", &["-db", "db/i_dvr", "-entry", "gnl|BL_ORD_ID|2", "-outfmt", "%o %a"], 0),
];

#[test]
fn long_ids_on_the_minus_strand_write_the_reverse_complement() {
    let dir = TestDir::new("approved");
    for &(id, words, status) in APPROVED {
        check(id, words, status, &dir.0);
    }
}

/// An accession on more than about 500 sequences (`PD-LOSAT-NCBI-DEFECTS` exception 6): NCBI
/// stops with `Accessions to Oids lookup error in e5` (probes H107-H109), LOSAT answers as NCBI
/// does for a key of fewer sequences (H110). `e5`: `XP_900001.i` on OID `i - 1`, 600 of them.
#[test]
fn an_accession_on_many_sequences_is_looked_up_like_any_other() {
    let dir = TestDir::new("many");
    check(
        "H107",
        &["-db", "db/e5", "-entry", "XP_900001", "-outfmt", "%o %a"],
        0,
        &dir.0,
    );
    check(
        "H108",
        &[
            "-db",
            "db/e5",
            "-entry_batch",
            "in/h_e.txt",
            "-outfmt",
            "%o %a",
        ],
        0,
        &dir.0,
    );
    check(
        "H109",
        &[
            "-db",
            "db/e5",
            "-entry_batch",
            "in/h_e.txt",
            "-get_dups",
            "-outfmt",
            "%o %a",
        ],
        0,
        &dir.0,
    );
}

/// Copies the files of the database `name` of `tests/blastdbcmd/db` into `dir`.
fn copy_db(name: &str, dir: &Path) {
    fs::create_dir_all(dir).unwrap();
    let db = tests_dir().join("blastdbcmd/db");
    for entry in fs::read_dir(&db).unwrap() {
        let path = entry.unwrap().path();
        if path.file_stem().unwrap().to_str() == Some(name) {
            fs::copy(&path, dir.join(path.file_name().unwrap())).unwrap();
        }
    }
}

/// A damaged lookup file (probes H78, H74, H86, H90-H106): NCBI's message at the opening of
/// the database, a silent miss of a GI without its ISAM, and refusals where NCBI crashes.
#[test]
fn damaged_lookup_files_are_reported_as_ncbi_s_or_refused() {
    let dir = TestDir::new("damaged");
    let damaged = |case: &str, change: &dyn Fn(&Path)| {
        let db = dir.0.join(case);
        copy_db("p5", &db);
        change(&db);
        // The platform's separator: NCBI's message writes the path with `\` on Windows.
        db.join("p5").display().to_string()
    };
    let open_error = |p5: &str, words: &[&str], message: String| {
        let mut all = vec!["-db", p5, "-dbtype", "prot"];
        all.extend_from_slice(words);
        let (stdout, stderr, code) = blastdbcmd(&all, &dir.0);
        assert_eq!(stdout, b"", "{all:?}");
        assert_eq!(
            String::from_utf8_lossy(&stderr),
            format!("{message}\n"),
            "{all:?}"
        );
        assert_eq!(code, 2, "{all:?}");
    };
    let nopdb = damaged("nopdb", &|db| fs::remove_file(db.join("p5.pdb")).unwrap());
    open_error(
        &nopdb,
        &["-entry", "all"],
        format!(
            "LOSAT database error: File {nopdb}.pdb not found. If you renamed any BLAST database \
             files, please use original file names, and makeblastdb to rename the database. If \
             you deleted any BLAST database files, you need to recreate the database."
        ),
    );
    let garbage = damaged("garbage", &|db| {
        fs::write(db.join("p5.pdb"), vec![0x5a; 8192]).unwrap();
    });
    for words in [&["-entry", "all"][..], &["-info"][..]] {
        open_error(
            &garbage,
            words,
            "LOSAT database error: LMDB runtime error: mdb_env_open: MDB_INVALID: File is not an \
             LMDB file"
                .to_string(),
        );
    }
    let empty = damaged("empty", &|db| fs::write(db.join("p5.pdb"), b"").unwrap());
    open_error(
        &empty,
        &["-entry", "all"],
        "LOSAT database error: LMDB runtime error: mdb_env_open: Bad file descriptor".to_string(),
    );
    let wrong = damaged("wrong", &|db| {
        let d5 = tests_dir().join("blastdbcmd/db/d5.pdb");
        fs::copy(d5, db.join("p5.pdb")).unwrap();
    });
    open_error(
        &wrong,
        &["-entry", "all"],
        "LOSAT database error: Input db vol does not match lmdb vol".to_string(),
    );

    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let unreadable = damaged("unreadable", &|db| {
            let mode = fs::Permissions::from_mode(0o000);
            fs::set_permissions(db.join("p5.pdb"), mode).unwrap();
        });
        // A superuser reads the file all the same.
        if fs::read(format!("{unreadable}.pdb")).is_err() {
            open_error(
                &unreadable,
                &["-entry", "all"],
                "LOSAT database error: LMDB runtime error: mdb_env_open: Permission denied"
                    .to_string(),
            );
        }
    }

    // NCBI stops with a bus error on a truncated LMDB file (H82).
    let trunc = damaged("trunc", &|db| {
        let bytes = fs::read(db.join("p5.pdb")).unwrap();
        fs::write(db.join("p5.pdb"), &bytes[..8192]).unwrap();
    });
    let (stdout, stderr, code) = blastdbcmd(&["-db", &trunc, "-entry", "all"], &dir.0);
    let stderr = String::from_utf8_lossy(&stderr);
    assert_eq!(stdout, b"");
    assert_eq!(code, 1, "{stderr}");
    assert!(
        stderr.ends_with("is not supported by LOSAT's BLASTDBCMD\n"),
        "{stderr}"
    );

    let missed = "Error: [blastdbcmd] Entry not found: 1819\nError: [blastdbcmd] Entry or entries \
                  not found in BLAST database\n";
    let nopnd = damaged("nopnd", &|db| fs::remove_file(db.join("p5.pnd")).unwrap());
    let nopni = damaged("nopni", &|db| fs::remove_file(db.join("p5.pni")).unwrap());
    let shortpnd = damaged("shortpnd", &|db| {
        let bytes = fs::read(db.join("p5.pnd")).unwrap();
        fs::write(db.join("p5.pnd"), &bytes[..bytes.len() - 8]).unwrap();
    });
    for (case, p5) in [("nopnd", nopnd), ("nopni", nopni), ("shortpnd", shortpnd)] {
        let (stdout, stderr, code) = blastdbcmd(&["-db", &p5, "-entry", "1819"], &dir.0);
        assert_eq!(stdout, b"", "{case}");
        assert_eq!(String::from_utf8_lossy(&stderr), missed, "{case}");
        assert_eq!(code, 1, "{case}");
    }
    // NCBI crashes on a GI index with a bad header (H98): only the lookups of a GI reach it.
    let badpni = damaged("badpni", &|db| {
        let mut bytes = fs::read(db.join("p5.pni")).unwrap();
        bytes[3] = 2;
        fs::write(db.join("p5.pni"), bytes).unwrap();
    });
    let (stdout, stderr, code) = blastdbcmd(&["-db", &badpni, "-entry", "1819"], &dir.0);
    let stderr = String::from_utf8_lossy(&stderr);
    assert_eq!(stdout, b"");
    assert_eq!(code, 1, "{stderr}");
    assert!(
        stderr.ends_with("is not supported by LOSAT's BLASTDBCMD\n"),
        "{stderr}"
    );
    let (stdout, _, code) = blastdbcmd(
        &["-db", &badpni, "-entry", "XP_000001.1", "-outfmt", "%o"],
        &dir.0,
    );
    assert_eq!(
        (String::from_utf8_lossy(&stdout).as_ref(), code),
        ("0\n", 0)
    );
}

/// A damaged ISAM of a version 4 database (probe I913-I1101): a missing file is a silent miss,
/// a string ISAM whose header does not fit its files is refused where NCBI crashes on it (the
/// lookups that open it: an accession, not a GI or an ordinal id), and so is a header with a
/// wrong count of lines, which NCBI does not read and LOSAT's reader refuses.
#[test]
fn damaged_isam_files_of_version_4_are_missed_or_refused() {
    let dir = TestDir::new("damaged4");
    let damaged = |case: &str, change: &dyn Fn(&Path)| {
        let db = dir.0.join(case);
        copy_db("p4", &db);
        change(&db);
        db.join("p4").display().to_string()
    };
    let set_byte = |db: &Path, file: &str, at: usize, value: u8| {
        let mut bytes = fs::read(db.join(file)).unwrap();
        bytes[at] = value;
        fs::write(db.join(file), bytes).unwrap();
    };
    let run = |p4: &str, id: &str| {
        let (stdout, stderr, code) =
            blastdbcmd(&["-db", p4, "-entry", id, "-outfmt", "%o %a"], &dir.0);
        (
            String::from_utf8(stdout).unwrap(),
            String::from_utf8(stderr).unwrap(),
            code,
        )
    };
    // OID 0 has three definition lines.
    let oid0 = || {
        (
            "0 XP_000001.1\n0 P99999.2\n0 AAA00001.1\n".to_string(),
            String::new(),
            0,
        )
    };
    let not_found = |id: &str| {
        (
            String::new(),
            format!(
                "Error: [blastdbcmd] Entry not found: {id}\nError: [blastdbcmd] Entry or entries not \
                 found in BLAST database\n"
            ),
            1,
        )
    };
    let refused = |(stdout, stderr, code): (String, String, i32)| {
        stdout.is_empty()
            && code == 1
            && stderr.lines().count() == 1
            && stderr.ends_with("is not supported by LOSAT's BLASTDBCMD\n")
    };

    for (case, file) in [("nopsi", "p4.psi"), ("nopsd", "p4.psd")] {
        let p4 = damaged(case, &|db| fs::remove_file(db.join(file)).unwrap());
        assert_eq!(run(&p4, "XP_000001.1"), not_found("XP_000001.1"), "{case}");
        assert_eq!(run(&p4, "1819").2, 0, "{case}");
    }
    let empty = damaged("emptypsi", &|db| fs::write(db.join("p4.psi"), b"").unwrap());
    let length = damaged("lenpsi", &|db| set_byte(db, "p4.psi", 11, 7));
    let version = damaged("verpsi", &|db| set_byte(db, "p4.psi", 3, 2));
    let terms = damaged("termspsi", &|db| set_byte(db, "p4.psi", 12, 16));
    for p4 in [&empty, &length, &version, &terms] {
        assert!(refused(run(p4, "XP_000001.1")), "{p4}");
        assert!(refused(run(p4, "NOSUCH")), "{p4}");
        let (stdout, stderr, code) = run(p4, "gnl|BL_ORD_ID|1");
        assert_eq!(
            (stdout.as_str(), stderr.as_str(), code),
            ("1 NP_000002.1\n", "", 0),
            "{p4}"
        );
        assert_eq!(run(p4, "gi|1819").2, 0, "{p4}");
    }
    let nopni = damaged("nopni", &|db| fs::remove_file(db.join("p4.pni")).unwrap());
    assert_eq!(run(&nopni, "gi|1819"), not_found("gi|1819"));
    assert_eq!(run(&nopni, "XP_000001.1"), oid0());
    let badpni = damaged("badpni", &|db| set_byte(db, "p4.pni", 3, 2));
    assert!(refused(run(&badpni, "1819")));
    assert_eq!(run(&badpni, "XP_000001.1"), oid0());
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

/// The probe's cases that LOSAT's BLASTDBCMD does not serve: `-range ""` (NCBI crashes) and a
/// `-mask_sequence_with` below -1.
#[rustfmt::skip]
const REFUSALS: &[(&str, &[&str])] = &[
    ("D61", &["-db", "db/p4n", "-entry", "gnl|BL_ORD_ID|3", "-range", ""]),
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
