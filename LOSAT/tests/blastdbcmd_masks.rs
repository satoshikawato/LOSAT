//! `losat blastdbcmd` on databases with masks (`makeblastdb -mask_data`): `-info`,
//! `-mask_sequence_with`, `-outfmt %m` and the mask id of a batch line against NCBI BLAST+
//! 2.17.0 blastdbcmd. The expected bytes are NCBI's, frozen by line BDC's probe M on ws
//! (`tests/blastdbcmd/README.md`), but for the approved differences 3 to 6 of
//! `PD-LOSAT-BLASTDBCMD-DEVIATIONS`. `{R}` stands for the directory of the probe's runs (here
//! `tests/blastdbcmd`, the current directory of each run). NCBI is not run here.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

/// `fs::canonicalize` without the `\\?\` prefix that Windows puts on its result (as
/// `blastdbcmd_retrieve.rs`).
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

fn run_dir() -> PathBuf {
    canonical(
        &Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("blastdbcmd"),
    )
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
            "losat_blastdbcmd_masks_{stem}_{}_{nanos}",
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

/// One run of `losat blastdbcmd words...` in `tests/blastdbcmd`: stdout, stderr, exit status.
/// No `BLASTDB`, no registry file (`NCBI_CONFIG_PATH` names `config` alone), no LOSAT variable,
/// nothing on the standard input.
fn blastdbcmd(words: &[&str], config: &Path) -> (Vec<u8>, Vec<u8>, i32) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text == "BLASTDB" {
            command.env_remove(&name);
        }
    }
    let output = command
        .env("NCBI_CONFIG_PATH", config)
        .stdin(std::process::Stdio::null())
        .current_dir(run_dir())
        .arg("blastdbcmd")
        .args(words)
        .output()
        .unwrap();
    (output.stdout, output.stderr, output.status.code().unwrap())
}

/// The frozen file `name` of `tests/blastdbcmd/expected`, with `{R}` placed.
fn expected(name: &str) -> Vec<u8> {
    let text = fs::read(run_dir().join("expected").join(name)).unwrap();
    let root = run_dir();
    let root = root.to_str().unwrap().as_bytes();
    let mut placed = Vec::with_capacity(text.len());
    let mut rest = &text[..];
    while let Some(at) = rest.windows(3).position(|window| window == b"{R}") {
        placed.extend_from_slice(&rest[..at]);
        placed.extend_from_slice(root);
        rest = &rest[at + 3..];
    }
    placed.extend_from_slice(rest);
    placed
}

/// The probe's cases: id, words, exit status. `db/nm5` and `db/nm4` have the dust masks 11 and
/// 12, `db/pm5` the seg masks 21, `db/palias` is an alias of `p5` and `pm5`, `db/m/nms4` has masks
/// 12 on one sequence, `db/m/pmg` an index of PIGs (a copy of its index of GIs). Approved
/// differences (expected bytes from NCBI's own neighbouring run, line BDC's
/// `probe/make_expected.py`): 3, a mask and `-strand minus` in FASTA (`M04`, `S31`, `M348`,
/// `M352`, `M372`, `M380`, `M392`, `M395`, `M397`, `M399`, `M401`, `M405`, `M406`, `M446`, `M449`,
/// `M452`, `M455`); 4, two `%m` (`M03`, `M130`, `M133`, `M137`, `M140`, `M142`); 5, a mask on a
/// protein in FASTA without a range (`M05`, `M285`, `M287`, `M288`, `M289`, `M339`, `M320`, `M323`,
/// `M324`, `M462`, `M465`); 6, the mask of the batch line before (`M425`, `M426`, `M428`).
#[rustfmt::skip]
const CASES: &[(&str, &[&str], i32)] = &[
    ("I16", &["-db", "db/nm5", "-info"], 0),
    ("I17", &["-db", "db/nm4", "-info"], 0),
    ("I18", &["-db", "db/pm5", "-info"], 0),
    ("I21", &["-db", "db/palias", "-info"], 0),
    ("I32", &["-db", "db/palias", "-info", "-exact_length"], 0),
    ("M123", &["-db", "db/nm5", "-info", "-dbtype", "nucl"], 0),
    ("M124", &["-db", "db/nm4", "-info", "-exact_length"], 0),
    ("M125", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%a %m11"], 0),
    ("M126", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%a %m12"], 0),
    ("M127", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%a %m0"], 1),
    ("M128", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%a %m99"], 1),
    ("M129", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%a %m"], 1),
    ("M131", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m11 %m"], 0),
    ("M132", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m %m11"], 1),
    ("M134", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m99 %m11"], 1),
    ("M135", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%%m11"], 1),
    ("M136", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%%m"], 1),
    ("M138", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%%m11 %m"], 0),
    ("M139", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m11 %%m12"], 0),
    ("M141", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m011"], 0),
    ("M146", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m11x"], 0),
    ("M147", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%s %m11 %h %l %o"], 0),
    ("M148", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%f %m11"], 0),
    ("M149", &["-db", "db/nm5", "-entry", "all", "-outfmt", "x%m11y"], 0),
    ("M150", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m11 %q"], 1),
    ("M151", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%q %m11"], 1),
    ("M152", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%a %m11", "-range", "100-200"], 0),
    ("M153", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%a %m11", "-strand", "minus"], 0),
    ("M154", &["-db", "db/nm5", "-entry", "NM_100001.1", "-outfmt", "%a %m11"], 0),
    ("M155", &["-db", "db/nm5", "-entry", "nomask3", "-outfmt", "%m11"], 0),
    ("M156", &["-db", "db/nm5", "-entry", "NM_100001.1,nomask3,AB100002.1", "-outfmt", "%a %m12"], 0),
    ("M157", &["-db", "db/nm4", "-entry", "all", "-outfmt", "%a %m11"], 0),
    ("M158", &["-db", "db/nm4", "-entry", "all", "-outfmt", "%a %m12 %s"], 0),
    ("M159", &["-db", "db/pm5", "-entry", "all", "-outfmt", "%a %m21"], 0),
    ("M160", &["-db", "db/pm5", "-entry", "all", "-outfmt", "%m21 %s"], 0),
    ("M163", &["-db", "db/m/nms4", "-entry", "all", "-outfmt", "%a %m12"], 0),
    ("M164", &["-db", "db/palias", "-entry", "all", "-outfmt", "%o %a %m21"], 0),
    ("M03", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%a %m11 %m12"], 0),
    ("M130", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m11 %m12"], 0),
    ("M133", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m11 %m99"], 1),
    ("M137", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%%m11 %m12"], 0),
    ("M140", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m11%m12"], 0),
    ("M142", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m0011 %m012"], 0),
    ("M01", &["-db", "db/nm5", "-entry", "all", "-mask_sequence_with", "11"], 0),
    ("M02", &["-db", "db/nm5", "-entry", "all", "-mask_sequence_with", "dust"], 1),
    ("M06", &["-db", "db/nm5", "-entry", "all", "-mask_sequence_with", "99"], 1),
    ("M199", &["-db", "db/nm5", "-entry", "all", "-mask_sequence_with=11"], 0),
    ("M200", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%a %s", "-mask_sequence_with=11"], 0),
    ("M203", &["-db", "db/nm5", "-entry", "all", "-mask_sequence_with=0"], 1),
    ("M204", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%a %s", "-mask_sequence_with=0"], 1),
    ("M205", &["-db", "db/nm5", "-entry", "all", "-mask_sequence_with="], 1),
    ("M207", &["-db", "db/nm5", "-entry", "all", "-mask_sequence_with=+11"], 0),
    ("M209", &["-db", "db/nm5", "-entry", "all", "-mask_sequence_with=011"], 0),
    ("M211", &["-db", "db/nm5", "-entry", "all", "-mask_sequence_with= 11"], 1),
    ("M215", &["-db", "db/nm5", "-entry", "all", "-mask_sequence_with=1e1"], 1),
    ("M219", &["-db", "db/nm5", "-entry", "all", "-mask_sequence_with=2147483648"], 1),
    ("M223", &["-db", "db/nm5", "-entry", "all", "-mask_sequence_with=-1"], 0),
    ("M224", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%a %s", "-mask_sequence_with=-1"], 0),
    ("M229", &["-db", "db/nm5", "-entry", "all", "-mask_sequence_with=dust"], 1),
    ("M235", &["-db", "db/nm5", "-entry", "all", "-mask_sequence_with=10"], 1),
    ("M271", &["-db", "db/nm4", "-entry", "all", "-mask_sequence_with", "11"], 0),
    ("M272", &["-db", "db/nm4", "-entry", "all", "-outfmt", "%s", "-mask_sequence_with", "11"], 0),
    ("M273", &["-db", "db/nm4", "-entry", "all", "-outfmt", "%f", "-mask_sequence_with", "11"], 0),
    ("M274", &["-db", "db/nm4", "-entry", "all", "-mask_sequence_with", "11", "-line_length", "30"], 0),
    ("M275", &["-db", "db/nm4", "-entry", "all", "-mask_sequence_with", "11", "-long_seqids"], 0),
    ("M276", &["-db", "db/nm4", "-entry", "all", "-mask_sequence_with", "11", "-range", "100-200"], 0),
    ("M277", &["-db", "db/nm4", "-entry", "all", "-outfmt", "%a %s", "-mask_sequence_with", "11", "-range", "100-200"], 0),
    ("M07", &["-db", "db/pm5", "-entry", "all", "-outfmt", "%a %s", "-mask_sequence_with", "21"], 0),
    ("M286", &["-db", "db/pm5", "-entry", "all", "-outfmt", "%s", "-mask_sequence_with", "21"], 0),
    ("M290", &["-db", "db/pm5", "-entry", "all", "-mask_sequence_with", "21", "-range", "100-200"], 0),
    ("M291", &["-db", "db/pm5", "-entry", "all", "-outfmt", "%a %s", "-mask_sequence_with", "21", "-range", "100-200"], 0),
    ("M403", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-strand", "plus"], 0),
    ("M404", &["-db", "db/nm5", "-entry", "NM_100001.1", "-outfmt", "%s", "-mask_sequence_with", "11", "-strand", "plus"], 0),
    ("M346", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "130-150", "-strand", "plus"], 0),
    ("M347", &["-db", "db/nm5", "-entry", "NM_100001.1", "-outfmt", "%a %s", "-mask_sequence_with", "11", "-range", "130-150", "-strand", "plus"], 0),
    ("M350", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "100-200", "-strand", "plus"], 0),
    ("M351", &["-db", "db/nm5", "-entry", "NM_100001.1", "-outfmt", "%a %s", "-mask_sequence_with", "11", "-range", "100-200", "-strand", "plus"], 0),
    ("M354", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "161-162", "-strand", "plus"], 0),
    ("M355", &["-db", "db/nm5", "-entry", "NM_100001.1", "-outfmt", "%a %s", "-mask_sequence_with", "11", "-range", "161-162", "-strand", "plus"], 0),
    ("M358", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "162-240", "-strand", "plus"], 0),
    ("M359", &["-db", "db/nm5", "-entry", "NM_100001.1", "-outfmt", "%a %s", "-mask_sequence_with", "11", "-range", "162-240", "-strand", "plus"], 0),
    ("M370", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "300-", "-strand", "plus"], 0),
    ("M371", &["-db", "db/nm5", "-entry", "NM_100001.1", "-outfmt", "%a %s", "-mask_sequence_with", "11", "-range", "300-", "-strand", "plus"], 0),
    ("M378", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "5000-6000", "-strand", "plus"], 0),
    ("M379", &["-db", "db/nm5", "-entry", "NM_100001.1", "-outfmt", "%a %s", "-mask_sequence_with", "11", "-range", "5000-6000", "-strand", "plus"], 0),
    ("M382", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "100-", "-strand", "plus"], 0),
    ("M383", &["-db", "db/nm5", "-entry", "NM_100001.1", "-outfmt", "%a %s", "-mask_sequence_with", "11", "-range", "100-", "-strand", "plus"], 0),
    ("M390", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "1-5000", "-strand", "plus"], 0),
    ("M391", &["-db", "db/nm5", "-entry", "NM_100001.1", "-outfmt", "%a %s", "-mask_sequence_with", "11", "-range", "1-5000", "-strand", "plus"], 0),
    ("M394", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "130-150", "-strand", "plus", "-long_seqids"], 0),
    ("M396", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "100-200", "-strand", "plus", "-long_seqids"], 0),
    ("M398", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "162-240", "-strand", "plus", "-long_seqids"], 0),
    ("M400", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "1-5000", "-strand", "plus", "-long_seqids"], 0),
    ("M402", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-strand", "plus", "-long_seqids"], 0),
    ("M04", &["-db", "db/nm4", "-entry", "NM_100001.1", "-mask_sequence_with", "12", "-range", "100-200", "-strand", "minus"], 0),
    ("S31", &["-db", "db/nm5", "-entry", "all", "-mask_sequence_with", "11", "-strand", "minus"], 0),
    ("S32", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%a %s", "-mask_sequence_with", "11", "-strand", "minus"], 0),
    ("M348", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "130-150", "-strand", "minus"], 0),
    ("M349", &["-db", "db/nm5", "-entry", "NM_100001.1", "-outfmt", "%a %s", "-mask_sequence_with", "11", "-range", "130-150", "-strand", "minus"], 0),
    ("M352", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "100-200", "-strand", "minus"], 0),
    ("M353", &["-db", "db/nm5", "-entry", "NM_100001.1", "-outfmt", "%a %s", "-mask_sequence_with", "11", "-range", "100-200", "-strand", "minus"], 0),
    ("M372", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "300-", "-strand", "minus"], 0),
    ("M373", &["-db", "db/nm5", "-entry", "NM_100001.1", "-outfmt", "%a %s", "-mask_sequence_with", "11", "-range", "300-", "-strand", "minus"], 0),
    ("M380", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "5000-6000", "-strand", "minus"], 0),
    ("M381", &["-db", "db/nm5", "-entry", "NM_100001.1", "-outfmt", "%a %s", "-mask_sequence_with", "11", "-range", "5000-6000", "-strand", "minus"], 0),
    ("M392", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "1-5000", "-strand", "minus"], 0),
    ("M393", &["-db", "db/nm5", "-entry", "NM_100001.1", "-outfmt", "%a %s", "-mask_sequence_with", "11", "-range", "1-5000", "-strand", "minus"], 0),
    ("M395", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "130-150", "-strand", "minus", "-long_seqids"], 0),
    ("M397", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "100-200", "-strand", "minus", "-long_seqids"], 0),
    ("M399", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "162-240", "-strand", "minus", "-long_seqids"], 0),
    ("M401", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-range", "1-5000", "-strand", "minus", "-long_seqids"], 0),
    ("M405", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-strand", "minus", "-long_seqids"], 0),
    ("M406", &["-db", "db/nm5", "-entry", "NM_100001.1", "-mask_sequence_with", "11", "-strand", "minus"], 0),
    ("M407", &["-db", "db/nm5", "-entry", "NM_100001.1", "-outfmt", "%s", "-mask_sequence_with", "11", "-strand", "minus"], 0),
    ("M05", &["-db", "db/pm5", "-entry", "all", "-mask_sequence_with", "21"], 0),
    ("M285", &["-db", "db/pm5", "-entry", "all", "-mask_sequence_with", "21"], 0),
    ("M287", &["-db", "db/pm5", "-entry", "all", "-outfmt", "%f", "-mask_sequence_with", "21"], 0),
    ("M288", &["-db", "db/pm5", "-entry", "all", "-mask_sequence_with", "21", "-line_length", "30"], 0),
    ("M289", &["-db", "db/pm5", "-entry", "all", "-mask_sequence_with", "21", "-long_seqids"], 0),
    ("M339", &["-db", "db/pm5", "-entry", "all", "-mask_sequence_with", "21", "-strand", "minus"], 0),
    ("M340", &["-db", "db/pm5", "-entry", "all", "-outfmt", "%s", "-mask_sequence_with", "21", "-strand", "minus"], 0),
    ("M320", &["-db", "db/palias", "-entry", "all", "-mask_sequence_with", "21"], 0),
    ("M321", &["-db", "db/palias", "-entry", "all", "-outfmt", "%s", "-mask_sequence_with", "21"], 0),
    ("M323", &["-db", "db/palias", "-entry", "all", "-mask_sequence_with", "21", "-line_length", "30"], 0),
    ("M324", &["-db", "db/palias", "-entry", "all", "-mask_sequence_with", "21", "-long_seqids"], 0),
    ("M420", &["-db", "db/m/nms4", "-entry", "all", "-mask_sequence_with", "12"], 0),
    ("M421", &["-db", "db/m/nms4", "-entry", "all", "-mask_sequence_with", "11"], 0),
    ("M422", &["-db", "db/m/nms4", "-entry", "all", "-outfmt", "%a %s", "-mask_sequence_with", "12"], 0),
    ("M423", &["-db", "db/m/nms4", "-entry", "AB100002.1,NM_100001.1", "-mask_sequence_with", "12"], 0),
    ("M424", &["-db", "db/m/nms4", "-entry", "NM_100001.1", "-mask_sequence_with", "12"], 0),
    ("M425", &["-db", "db/m/nms4", "-entry_batch", "in/m_b1.txt"], 0),
    ("M426", &["-db", "db/m/nms4", "-entry_batch", "in/m_b1.txt", "-get_dups"], 0),
    ("M427", &["-db", "db/m/nms4", "-entry_batch", "in/m_b1.txt", "-outfmt", "%a %s"], 0),
    ("M428", &["-db", "db/m/nms4", "-entry_batch", "in/m_b3.txt"], 0),
    ("M429", &["-db", "db/m/nms4", "-entry_batch", "in/m_b3.txt", "-long_seqids"], 1),
    ("M933", &["-db", "db/m/nms4", "-entry_batch", "in/m_b3.txt", "-outfmt", "%a %s"], 0),
    ("M446", &["-db", "db/nm5", "-entry_batch", "in/m_b2.txt"], 1),
    ("M447", &["-db", "db/nm5", "-entry_batch", "in/m_b2.txt", "-outfmt", "%a %s"], 1),
    ("M448", &["-db", "db/nm5", "-entry_batch", "in/m_b2.txt", "-outfmt", "%a %m11"], 1),
    ("M449", &["-db", "db/nm5", "-entry_batch", "in/m_b2.txt", "-get_dups"], 1),
    ("M450", &["-db", "db/nm5", "-entry_batch", "in/m_b2.txt", "-outfmt", "%a %s", "-get_dups"], 1),
    ("M451", &["-db", "db/nm5", "-entry_batch", "in/m_b2.txt", "-outfmt", "%a %m11", "-get_dups"], 1),
    ("M452", &["-db", "db/nm4", "-entry_batch", "in/m_b2.txt"], 1),
    ("M453", &["-db", "db/nm4", "-entry_batch", "in/m_b2.txt", "-outfmt", "%a %s"], 1),
    ("M454", &["-db", "db/nm4", "-entry_batch", "in/m_b2.txt", "-outfmt", "%a %m11"], 1),
    ("M455", &["-db", "db/nm4", "-entry_batch", "in/m_b2.txt", "-get_dups"], 1),
    ("M456", &["-db", "db/nm4", "-entry_batch", "in/m_b2.txt", "-outfmt", "%a %s", "-get_dups"], 1),
    ("M457", &["-db", "db/nm4", "-entry_batch", "in/m_b2.txt", "-outfmt", "%a %m11", "-get_dups"], 1),
    ("M458", &["-db", "db/m/pmg", "-ipg", "1001", "-outfmt", "%o %a %m21"], 0),
    ("M459", &["-db", "db/m/pmg", "-ipg", "9999", "-outfmt", "%o %m21"], 1),
    ("M460", &["-db", "db/m/pmg", "-ipg", "9999", "-outfmt", "x%o|%m21"], 1),
    ("M461", &["-db", "db/m/pmg", "-ipg", "9999", "-outfmt", "%m21|%o"], 1),
    ("M462", &["-db", "db/m/pmg", "-ipg", "1001", "-mask_sequence_with", "21"], 0),
    ("M463", &["-db", "db/m/pmg", "-ipg", "9999", "-mask_sequence_with", "21"], 1),
    ("M464", &["-db", "db/m/pmg", "-ipg", "1001", "-outfmt", "%s", "-mask_sequence_with", "21"], 0),
    ("M465", &["-db", "db/m/pmg", "-ipg_batch", "in/m_b4.txt"], 1),
    ("M466", &["-db", "db/m/pmg", "-ipg_batch", "in/m_b4.txt", "-outfmt", "%o %s"], 1),
    ("M467", &["-db", "db/m/pmg", "-ipg_batch", "in/m_b4.txt", "-outfmt", "%o %m21"], 1),
    ("M941", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m11%"], 1),
    ("M942", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%o %m11 %"], 1),
    ("M943", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%%m11%"], 1),
    ("M944", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m11%é"], 1),
    ("M945", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%é%m11"], 1),
    ("M946", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m11%%"], 0),
    ("M947", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%é %m11 %m12"], 1),
];

#[test]
fn masks_are_ncbi_s_bytes() {
    let config = TestDir::new("cases");
    for &(id, words, status) in CASES {
        // The lines of `-info` name a path with the platform's separator.
        if cfg!(windows) && words.contains(&"-info") {
            continue;
        }
        let (stdout, stderr, code) = blastdbcmd(words, &config.0);
        let (want_out, want_err) = (
            expected(&format!("{id}.stdout")),
            expected(&format!("{id}.stderr")),
        );
        assert_eq!(
            String::from_utf8_lossy(&stdout),
            String::from_utf8_lossy(&want_out),
            "{id} {words:?}: stdout"
        );
        assert_eq!(
            String::from_utf8_lossy(&stderr),
            String::from_utf8_lossy(&want_err),
            "{id} {words:?}: stderr"
        );
        // `-outfmt %é` puts one byte of the letter in NCBI's message.
        assert!(
            stdout == want_out && stderr == want_err,
            "{id} {words:?}: bytes"
        );
        assert_eq!(code, status, "{id} {words:?}: exit status");
    }
}

/// What NCBI does not answer: an algorithm id that overflows its `int`, a negative id it looks
/// up, an algorithm of the alias `mask_remap` that the sequence's volume does not know (NCBI
/// reads the masks of another algorithm).
#[rustfmt::skip]
const REFUSED: &[(&str, &[&str])] = &[
    ("M143", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m4294967307"]),
    ("M144", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m99999999999"]),
    ("M145", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%m2147483648"]),
    ("M225", &["-db", "db/nm5", "-entry", "all", "-mask_sequence_with=-2"]),
    ("M226", &["-db", "db/nm5", "-entry", "all", "-outfmt", "%a %s", "-mask_sequence_with=-2"]),
    ("M180", &["-db", "../seqdb/db/mask_remap", "-entry", "all", "-outfmt", "%o|%a|%l|%m11"]),
    ("M261", &["-db", "../seqdb/db/mask_remap", "-entry", "all", "-mask_sequence_with", "11"]),
];

#[test]
fn what_ncbi_does_not_answer_is_refused() {
    let config = TestDir::new("refused");
    for &(id, words) in REFUSED {
        let (stdout, stderr, code) = blastdbcmd(words, &config.0);
        let stderr = String::from_utf8_lossy(&stderr);
        assert!(stdout.is_empty(), "{id} {words:?}: stdout");
        assert_eq!(stderr.lines().count(), 1, "{id} {words:?}: {stderr}");
        assert!(
            stderr.contains("LOSAT's BLASTDBCMD"),
            "{id} {words:?}: {stderr}"
        );
        assert_eq!(code, 1, "{id} {words:?}: exit status");
    }
}

/// A copy of `db/nm4` named `b` in `dir`: its path, as `-db` takes it.
fn copy_of_nm4(dir: &Path) -> String {
    for ext in ["naa", "nab", "nac", "nhr", "nin", "nsq", "nsd", "nsi"] {
        fs::copy(
            run_dir().join("db").join(format!("nm4.{ext}")),
            dir.join(format!("b.{ext}")),
        )
        .unwrap();
    }
    dir.join("b").to_str().unwrap().to_string()
}

/// `bytes` written over `file` at the first place of `find`.
fn patch(file: &Path, find: &[u8], bytes: &[u8]) {
    let mut data = fs::read(file).unwrap();
    let at = data
        .windows(find.len())
        .position(|window| window == find)
        .unwrap();
    data[at..at + bytes.len()].copy_from_slice(bytes);
    fs::write(file, data).unwrap();
}

/// The little-endian offsets of a pair of ranges, as the data file `.nac` stores them.
fn offsets(numbers: &[u32]) -> Vec<u8> {
    numbers.iter().flat_map(|n| n.to_le_bytes()).collect()
}

/// `-info` stops at the filtering algorithms and `%m` at its algorithm id, with NCBI's text
/// (probes M722, M732, M742, M752, M762, M812 and M724, M764, M814), for a damaged index file.
#[test]
fn a_damaged_mask_column_is_ncbi_s_error() {
    let config = TestDir::new("damaged-config");
    let wrong: &[(&str, usize, &[u8], &str)] = &[
        (
            "version",
            0,
            &[0, 0, 0, 2],
            "Column file uses unknown format_version.",
        ),
        (
            "type",
            4,
            &[0, 0, 0, 2],
            "Column file uses unknown data type.",
        ),
        (
            "offset",
            8,
            &[0, 0, 0, 8],
            "Column file uses unsupported offset size.",
        ),
        (
            "oids",
            12,
            &[0, 0, 0, 5],
            "Error: column has wrong #oids [5 vs 4].",
        ),
        (
            "metastart",
            24,
            &[0, 0, 0, 77],
            "CSeqDBColumn: File format error.",
        ),
    ];
    for &(name, at, bytes, message) in wrong {
        let dir = TestDir::new(name);
        let db = copy_of_nm4(&dir.0);
        let index = dir.0.join("b.naa");
        let mut data = fs::read(&index).unwrap();
        data[at..at + bytes.len()].copy_from_slice(bytes);
        fs::write(&index, data).unwrap();

        let (stdout, stderr, code) = blastdbcmd(&["-db", &db, "-info"], &config.0);
        let stdout = String::from_utf8(stdout).unwrap();
        assert!(
            stdout.ends_with("\n\nBLASTDB Version: 4\n"),
            "{name}: {stdout}"
        );
        assert_eq!(
            String::from_utf8_lossy(&stderr),
            format!("BLAST Database error: {message}\n"),
            "{name}"
        );
        assert_eq!(code, 2, "{name}");

        let (stdout, stderr, code) = blastdbcmd(
            &["-db", &db, "-entry", "all", "-outfmt", "%o %m11"],
            &config.0,
        );
        assert!(stdout.is_empty(), "{name}");
        assert_eq!(
            String::from_utf8_lossy(&stderr),
            format!("Error: [blastdbcmd] {message}\n"),
            "{name}"
        );
        assert_eq!(code, 1, "{name}");

        // A request that does not read the masks does not open the column.
        let (stdout, stderr, code) = blastdbcmd(
            &["-db", &db, "-entry", "all", "-outfmt", "%o %l"],
            &config.0,
        );
        assert_eq!(stdout, b"0 400\n1 449\n2 300\n3 250\n", "{name}");
        assert!(stderr.is_empty(), "{name}");
        assert_eq!(code, 0, "{name}");
    }
}

/// The errors of what follows the fixed fields (M948-M954, M955-M961) and of a mask that runs
/// past its data (M992, M994), with NCBI's text.
#[test]
fn damaged_metadata_and_masks_are_ncbi_s_errors() {
    let config = TestDir::new("metadata-config");
    let dir = TestDir::new("dupkey");
    let db = copy_of_nm4(&dir.0);
    let index = dir.0.join("b.naa");
    let data = fs::read(&index).unwrap();
    let value = data
        .windows(12)
        .position(|window| window == b"10:window=32")
        .unwrap();
    let key = data[..value]
        .windows(2)
        .rposition(|window| window == b"12")
        .unwrap();
    let mut patched = data.clone();
    patched[key..key + 2].copy_from_slice(b"11");
    fs::write(&index, patched).unwrap();
    let (_, stderr, code) = blastdbcmd(&["-db", &db, "-info"], &config.0);
    assert_eq!(
        String::from_utf8_lossy(&stderr),
        "BLAST Database error: CSeqDBColumn: Error; duplicate metadata key.\n"
    );
    assert_eq!(code, 2);

    let mut patched = data;
    let same = b"10:window=64; level=20; linker=1";
    patched[value..value + same.len()].copy_from_slice(same);
    fs::write(&index, patched).unwrap();
    let (_, stderr, code) = blastdbcmd(
        &["-db", &db, "-entry", "all", "-outfmt", "%o %m11"],
        &config.0,
    );
    let stderr = String::from_utf8_lossy(&stderr);
    assert!(
        stderr.starts_with("Error: [blastdbcmd] Error: volume (")
            && stderr
                .ends_with(") mask data has duplicates value (10:window=64; level=20; linker=1)\n"),
        "{stderr}"
    );
    assert_eq!(code, 1);

    let dir = TestDir::new("npairs");
    let db = copy_of_nm4(&dir.0);
    let data_file = dir.0.join("b.nac");
    let data = fs::read(&data_file).unwrap();
    let pair = offsets(&[120, 161]);
    let at = data
        .windows(pair.len())
        .position(|window| window == pair.as_slice())
        .unwrap();
    let mut patched = data;
    patched[at - 4..at].copy_from_slice(&1000i32.to_be_bytes());
    fs::write(&data_file, patched).unwrap();
    let (stdout, stderr, code) = blastdbcmd(
        &["-db", &db, "-entry", "all", "-outfmt", "%o %m11"],
        &config.0,
    );
    assert!(stdout.is_empty());
    assert_eq!(
        String::from_utf8_lossy(&stderr),
        "Error: [blastdbcmd] CBlastDbBlob::x_ReadRaw: hit end of data\n"
    );
    assert_eq!(code, 1);
}

/// A column without its data files is no column (M772, M774); one with the big-endian data file
/// alone (M784), another column before the masks (M824) and, for `%s`, a stored mask outside the
/// sequence (M886, M896) are refused, where NCBI reads them as something else or crashes. The
/// FASTA path clips that mask (M887, M897).
#[test]
fn columns_that_ncbi_misreads_are_refused() {
    let config = TestDir::new("misread-config");
    let refused = |db: &str, words: &[&str]| {
        let mut all = vec!["-db", db, "-entry", "all"];
        all.extend_from_slice(words);
        let (stdout, stderr, code) = blastdbcmd(&all, &config.0);
        let stderr = String::from_utf8_lossy(&stderr);
        assert!(stdout.is_empty(), "{all:?}");
        assert!(stderr.contains("LOSAT's BLASTDBCMD"), "{all:?}: {stderr}");
        assert_eq!(code, 1, "{all:?}");
    };

    let dir = TestDir::new("nodata");
    let db = copy_of_nm4(&dir.0);
    fs::remove_file(dir.0.join("b.nab")).unwrap();
    fs::remove_file(dir.0.join("b.nac")).unwrap();
    let (stdout, _, code) = blastdbcmd(&["-db", &db, "-info"], &config.0);
    assert!(!String::from_utf8_lossy(&stdout).contains("Available filtering"));
    assert_eq!(code, 0);
    let (_, stderr, code) =
        blastdbcmd(&["-db", &db, "-entry", "all", "-outfmt", "%m11"], &config.0);
    assert_eq!(
        String::from_utf8_lossy(&stderr),
        "Error: [blastdbcmd] Invalid filtering algorithm ID: 11\n\
         Error: [blastdbcmd] Invalid filtering algorithm ID for outfmt %m.\n"
    );
    assert_eq!(code, 1);

    let dir = TestDir::new("big");
    let db = copy_of_nm4(&dir.0);
    fs::remove_file(dir.0.join("b.nac")).unwrap();
    refused(&db, &["-outfmt", "%o %m11"]);

    let dir = TestDir::new("other");
    let db = copy_of_nm4(&dir.0);
    for ext in ["a", "b", "c"] {
        fs::copy(
            dir.0.join(format!("b.na{ext}")),
            dir.0.join(format!("b.nb{ext}")),
        )
        .unwrap();
    }
    patch(
        &dir.0.join("b.naa"),
        b"BlastDb/MaskData",
        b"BlastDb/MaskDatX",
    );
    refused(&db, &["-outfmt", "%o %m11"]);

    let dir = TestDir::new("outside");
    let db = copy_of_nm4(&dir.0);
    patch(
        &dir.0.join("b.nac"),
        &offsets(&[120, 161]),
        &offsets(&[120, 100_000]),
    );
    refused(&db, &["-outfmt", "%o %s", "-mask_sequence_with", "11"]);
    let (stdout, _, code) = blastdbcmd(
        &[
            "-db",
            &db,
            "-entry",
            "NM_100001.1",
            "-mask_sequence_with",
            "11",
        ],
        &config.0,
    );
    let letters: Vec<u8> = stdout
        .split(|&byte| byte == b'\n')
        .skip(1)
        .flatten()
        .copied()
        .collect();
    assert_eq!(code, 0);
    assert_eq!(letters.len(), 400);
    assert!(letters[..120].iter().all(u8::is_ascii_uppercase));
    assert!(letters[120..].iter().all(u8::is_ascii_lowercase));
}
