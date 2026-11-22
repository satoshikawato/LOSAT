//! `losat blastdbcmd` on alias files with the keys of line RG, series D, item D1a
//! (`remaining-gaps-20261112/D/DESIGN_D1a.md`): `NSEQ` and `LENGTH` read as NCBI's walks of the alias
//! tree (values that are not the sums of the volumes included), `STATS_NSEQ`, `STATS_TOTLEN` and
//! `MIN_SEQ_LENGTH` (read, never printed), `FIRST_OID` / `LAST_OID` ranges (nested ranges ANDed, two
//! alias files over one volume ORed, the global OIDs of several volumes), `VDBLIST` (not read), a
//! volume that an alias file and a name both give, and NCBI's `CStringException` for a value that
//! does not convert, and (item D1a-2) `-ipg` / `-ipg_batch` for a PIG whose OID an alias file's
//! range leaves out (NCBI's `PigToOid` asks the volumes, and its formatters read that OID), and (D1a-2
//! part 2) `FIRST_OID` / `LAST_OID` of 2^31 or more and the iteration range of the only alias file
//! of a `-db`, as each sub-command reads NCBI's OID list (`make_d1a2a_fixtures.py`), and (D1b-1)
//! `OIDLIST` (`make_d1b_blastdbcmd.py`, 2026-11-21) and (D1b-2) `OID_MASK_TYPE`
//! (`make_d1b2_blastdbcmd.py`, 2026-11-21), and (D1c-1) `GILIST`, `TILIST` and `SEQIDLIST`
//! over the databases of `tests/db_search/idlist` and the lists over `db/pz4`'s header of no
//! bytes (`db/pzl*`, `make_d1c_blastdbcmd.py`, 2026-11-22; review 1 of D1c-1: string tags of
//! `PIG` and `BL_ORD_ID` under a list, TI and GI ISAMs NCBI does not initialise, a `SEQIDLIST`
//! over two volumes, A382-A389). The
//! expected bytes are NCBI BLAST+ 2.17.0 blastdbcmd's (the release build on the developer machine,
//! 2026-11-12 and 2026-11-13, `remaining-gaps-20261112/D/fixtures/make_d1a_blastdbcmd.py` and
//! `make_d1a2_blastdbcmd.py`); the alias files are in `tests/blastdbcmd/db/ak/`, the batch files
//! in `tests/blastdbcmd/in/`. `{R}` stands for the directory of the runs (`tests/blastdbcmd`), `{T}`
//! for the directory `tests`. NCBI is not run here.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

/// `fs::canonicalize` without the `\\\\?\\` prefix that Windows puts on its result (as
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
            "losat_blastdbcmd_alias_keys_{stem}_{}_{nanos}",
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

/// The frozen file `name` of `tests/blastdbcmd/expected`, as it is.
fn frozen(name: &str) -> Vec<u8> {
    fs::read(run_dir().join("expected").join(name)).unwrap()
}

/// `text` with `{R}` and `{T}` placed, in LOSAT's words (NCBI's bytes through `losat::wording`,
/// AGENTS.md rule 8).
fn placed(text: &[u8]) -> Vec<u8> {
    let root = run_dir();
    let tests = root.parent().unwrap().to_str().unwrap().as_bytes().to_vec();
    let root = root.to_str().unwrap().as_bytes().to_vec();
    let mut placed = Vec::with_capacity(text.len());
    let mut rest = text;
    while let Some(at) = rest
        .windows(3)
        .position(|window| window == b"{R}" || window == b"{T}")
    {
        placed.extend_from_slice(&rest[..at]);
        placed.extend_from_slice(if rest[at + 1] == b'R' { &root } else { &tests });
        rest = &rest[at + 3..];
    }
    placed.extend_from_slice(rest);
    losat::wording::losat_output_bytes(&placed).into_owned()
}

/// Whether `text` names a path under the run directory (`{R}`) or the tests directory (`{T}`).
fn holds_root(text: &[u8]) -> bool {
    text.windows(3)
        .any(|window| window == b"{R}" || window == b"{T}")
}

/// The cases: id, words, exit status.
#[rustfmt::skip]
const CASES: &[(&str, &[&str], i32)] = &[
    ("A01", &["-db", "db/ak/p_nseq", "-info"], 0),
    ("A02", &["-db", "db/ak/p_nseq", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A03", &["-db", "db/ak/j_nseq", "-metadata"], 0),
    ("A04", &["-db", "db/ak/p_range", "-info"], 0),
    ("A05", &["-db", "db/ak/p_range", "-info", "-exact_length"], 0),
    ("A06", &["-db", "db/ak/p_range", "-entry", "all", "-outfmt", "%o %a %l %T"], 0),
    ("A07", &["-db", "db/ak/p_range", "-entry", "XP_000001.1", "-outfmt", "%o %a"], 1),
    ("A08", &["-db", "db/ak/p_range", "-entry", "Q00003,A0A000", "-outfmt", "%o %a"], 0),
    ("A09", &["-db", "db/ak/p_range", "-entry", "all"], 0),
    ("A10", &["-db", "db/ak/p_union", "-info"], 0),
    ("A11", &["-db", "db/ak/p_union", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A12", &["-db", "db/ak/p_twice", "-info"], 0),
    ("A13", &["-db", "db/ak/p_twice", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A14", &["-db", "db/ak/p_stats", "-info"], 0),
    ("A15", &["-db", "db/ak/p_nseq0", "-info"], 0),
    ("A16", &["-db", "db/ak/d_glob", "-info"], 0),
    ("A17", &["-db", "db/ak/d_glob", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A18", &["-db", "db/ak/d_rev", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A19", &["-db", "db/ak/n_vdb", "-info"], 0),
    ("A20", &["-db", "db/ak/n_vdbonly", "-info", "-dbtype", "nucl"], 2),
    ("A21", &["-db", "db/ak/p_bad_nseq", "-info"], 255),
    ("A22", &["-db", "db/ak/p_bad_nseq", "-info", "-dbtype", "prot"], 255),
    ("A23", &["-db", "db/ak/p_bad_first", "-entry", "all", "-dbtype", "prot"], 255),
    ("A24", &["-db", "db/ak/p_nest", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A25", &["-db", "db/ak/p_empty", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A26", &["-db", "db/ak/p_empty", "-info"], 0),
    ("A27", &["-db", "db/ak/p_nseq4g", "-info"], 0),
    ("A28", &["-db", "db/ak/p_min", "-info"], 0),
    ("A29", &["-db", "db/ak/j_range", "-metadata"], 0),
    ("A30", &["-db", "db/palias db/p5", "-info"], 0),
    ("A31", &["-db", "db/palias db/p5", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A32", &["-db", "db/palias db/p5", "-entry", "P99999.2", "-outfmt", "%o %a"], 0),
    // D1a-2: `-ipg` / `-ipg_batch` read the OID of a PIG whether the alias files' ranges include
    // it or not (`PigToOid` asks the volumes; probe p4c_c 0-25 of D1a-2). A52 is approved difference
    // 5 of PD-LOSAT-BLASTDBCMD-DEVIATIONS (a mask on a protein in FASTA without a range): its
    // expected bytes are NCBI's header with the letters of NCBI's `-outfmt %s` run of the request.
    ("A33", &["-db", "db/ak/m_r", "-ipg", "1001"], 0),
    ("A34", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%o|%a|%g|%t|%l|%P|%T|%L|%e|%n|%h|%i"], 0),
    ("A35", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%s"], 0),
    ("A36", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%m21"], 0),
    ("A37", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%d"], 0),
    ("A38", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%b"], 0),
    ("A39", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%f"], 0),
    ("A40", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%X|%K|%S|%B|%C|%N"], 0),
    ("A41", &["-db", "db/ak/m_r", "-ipg", "1002", "-outfmt", "%o|%a|%l"], 0),
    ("A42", &["-db", "db/m/pmg", "-ipg", "1001", "-outfmt", "%o|%a|%g|%t|%l|%P|%T|%L|%e|%n|%h|%i"], 0),
    ("A43", &["-db", "db/m/pmg", "-ipg", "1001"], 0),
    ("A44", &["-db", "db/ak/m_s", "-ipg", "1003", "-outfmt", "%o|%a|%l|%s"], 0),
    ("A45", &["-db", "db/ak/m_s", "-ipg", "1003"], 0),
    ("A46", &["-db", "db/ak/m_t", "-ipg", "1001", "-outfmt", "%o|%a|%l"], 0),
    ("A47", &["-db", "db/ak/m_r", "-ipg", "9999"], 1),
    ("A48", &["-db", "db/ak/m_r", "-ipg", "9999", "-outfmt", "%o|%l"], 1),
    ("A49", &["-db", "db/ak/m_r", "-ipg", "9999", "-outfmt", "%o"], 0),
    ("A50", &["-db", "db/ak/m_r", "-ipg", "9999", "-outfmt", "%s"], 1),
    ("A51", &["-db", "db/ak/m_r", "-ipg", "9999", "-outfmt", "%a"], 1),
    ("A52", &["-db", "db/ak/m_r", "-ipg", "1001", "-mask_sequence_with", "21"], 0),
    ("A53", &["-db", "db/ak/m_r", "-ipg", "1001", "-range", "2-10"], 0),
    ("A54", &["-db", "db/ak/m_r", "-ipg", "1001", "-outfmt", "%s", "-range", "2-10", "-mask_sequence_with", "21"], 0),
    ("A55", &["-db", "db/ak/m_r", "-ipg_batch", "in/ak_pig.txt"], 1),
    ("A56", &["-db", "db/ak/m_r", "-ipg_batch", "in/ak_pig.txt", "-outfmt", "%o %a %l"], 1),
    ("A57", &["-db", "db/ak/m_r", "-ipg_batch", "in/ak_pigm.txt", "-outfmt", "%o %s"], 0),
    ("A58", &["-db", "db/ak/m_s", "-ipg_batch", "in/ak_pig.txt", "-outfmt", "%o %a"], 1),
    // D1a-2: `LENGTH` and `STATS_TOTLEN` over volumes of no OID (`db/ak/ep`, a copy of
    // `seqdb/db/empty`, and `db/ak/en`, a nucleotide index of no sequence): the alias file's
    // counts, no sequence (probes p4c_dp, p4c_dn of D1a-2).
    ("A59", &["-db", "db/ak/e_pB", "-dbtype", "prot", "-info"], 0),
    ("A60", &["-db", "db/ak/e_pB", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A61", &["-db", "db/ak/e_pA", "-dbtype", "prot", "-info"], 0),
    ("A62", &["-db", "db/ak/e_nB", "-dbtype", "nucl", "-info"], 0),
    ("A63", &["-db", "db/ak/e_nD", "-dbtype", "nucl", "-info"], 0),
    ("A64", &["-db", "db/ak/en", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A65", &["-db", "db/ak/e_nC", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o %a"], 0),
    // D1a-2 part 2 (survey A): the probe cases p4a_* over `db/ak/or/` (copies of taxfilter's `tfp` and
    // `tfn`, and `o`, whose OIDs come first in a `-db` with `tfp`): FIRST_OID / LAST_OID of 2^31 or
    // more, and the iteration range of the only alias file of a `-db`. `-entry all` and `-tax_info`
    // follow NCBI's loop (its first check before the OID list is set up can pass the first listed
    // OID, and the OID at the range's end passes when listed), `-entry` the included OIDs,
    // `-entry_batch` and `-taxids` the first check before the list is set up, `-target_only` every
    // OID of the list. A183-A190: GIs through `p_range` (`GiToOidwFilterCheck`: NCBI writes the OID
    // its check moved to when it does not keep the GI's).
    ("A66", &["-db", "db/ak/or/A_lo3g_a", "-dbtype", "prot", "-info"], 0),
    ("A67", &["-db", "db/ak/or/A_lo3g_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A68", &["-db", "db/ak/or/A_lo3g_a", "-dbtype", "prot", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%o %a"], 0),
    ("A69", &["-db", "db/ak/or/A_lo3g_a", "-dbtype", "prot", "-taxids", "9606", "-no_taxid_expansion", "-target_only", "-outfmt", "%o %a"], 0),
    ("A70", &["-db", "db/ak/or/A_lo3g_a", "-dbtype", "prot", "-taxids", "10090", "-no_taxid_expansion", "-outfmt", "%o %a"], 1),
    ("A71", &["-db", "db/ak/or/A_lo3g_a", "-dbtype", "prot", "-taxids", "10090", "-no_taxid_expansion", "-target_only", "-outfmt", "%o %a"], 0),
    ("A72", &["-db", "db/ak/or/A_lo3g_a", "-dbtype", "prot", "-tax_info"], 0),
    ("A73", &["-db", "db/ak/or/A_lo3g_a", "-dbtype", "prot", "-metadata"], 0),
    ("A74", &["-db", "db/ak/or/A_lo3g_a", "-dbtype", "prot", "-entry", "P09105.2", "-outfmt", "%o %a"], 1),
    ("A75", &["-db", "db/ak/or/A_lo3g_a", "-dbtype", "prot", "-entry", "P69905.2", "-outfmt", "%o %a"], 1),
    ("A76", &["-db", "db/ak/or/A_lo3g_a", "-dbtype", "prot", "-entry_batch", "in/or_A_lo3g_b1.txt", "-outfmt", "%o %a"], 1),
    ("A77", &["-db", "db/ak/or/A_lo3g_a", "-dbtype", "prot", "-entry_batch", "in/or_A_lo3g_b2.txt", "-outfmt", "%o %a"], 1),
    ("A78", &["-db", "db/ak/or/A_lo3g_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%f"], 0),
    ("A79", &["-db", "db/ak/or/A_fo3g_a", "-dbtype", "prot", "-info"], 0),
    ("A80", &["-db", "db/ak/or/A_fo3g_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A81", &["-db", "db/ak/or/A_fo3g_a", "-dbtype", "prot", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%o %a"], 1),
    ("A82", &["-db", "db/ak/or/A_fo3g_a", "-dbtype", "prot", "-taxids", "9606", "-no_taxid_expansion", "-target_only", "-outfmt", "%o %a"], 1),
    ("A83", &["-db", "db/ak/or/A_fo3_lo3g_a", "-dbtype", "prot", "-info"], 0),
    ("A84", &["-db", "db/ak/or/A_fo3_lo3g_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A85", &["-db", "db/ak/or/A_fo3_lo3g_a", "-dbtype", "prot", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%o %a"], 1),
    ("A86", &["-db", "db/ak/or/A_fo3_lo3g_a", "-dbtype", "prot", "-taxids", "9606", "-no_taxid_expansion", "-target_only", "-outfmt", "%o %a"], 0),
    ("A87", &["-db", "db/ak/or/A_fo1_lo3g_a", "-dbtype", "prot", "-info"], 0),
    ("A88", &["-db", "db/ak/or/A_fo1_lo3g_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A89", &["-db", "db/ak/or/A_fo3g_lo5_a", "-dbtype", "prot", "-info"], 0),
    ("A90", &["-db", "db/ak/or/A_fo3g_lo5_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A91", &["-db", "db/ak/or/A_fo2p31_a", "-dbtype", "prot", "-info"], 0),
    ("A92", &["-db", "db/ak/or/A_fo2p31_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A93", &["-db", "db/ak/or/A_fo2p31p1_a", "-dbtype", "prot", "-info"], 0),
    ("A94", &["-db", "db/ak/or/A_fo2p31p1_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A95", &["-db", "db/ak/or/A_fomax_a", "-dbtype", "prot", "-info"], 0),
    ("A96", &["-db", "db/ak/or/A_fomax_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A97", &["-db", "db/ak/or/A_lo2p31_a", "-dbtype", "prot", "-info"], 0),
    ("A98", &["-db", "db/ak/or/A_lo2p31_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A99", &["-db", "db/ak/or/A_lomax_a", "-dbtype", "prot", "-info"], 0),
    ("A100", &["-db", "db/ak/or/A_lomax_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A101", &["-db", "db/ak/or/A_fo2p31m1_a", "-dbtype", "prot", "-info"], 0),
    ("A102", &["-db", "db/ak/or/A_fo2p31m1_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A103", &["-db", "db/ak/or/A_lo2p31m1_a", "-dbtype", "prot", "-info"], 0),
    ("A104", &["-db", "db/ak/or/A_lo2p31m1_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A105", &["-db", "db/ak/or/A_nest_lo3g_t", "-dbtype", "prot", "-info"], 0),
    ("A106", &["-db", "db/ak/or/A_nest_lo3g_t", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A107", &["-db", "db/ak/or/A_nest_lo3g_t", "-dbtype", "prot", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%o %a"], 0),
    ("A108", &["-db", "db/ak/or/A_nest_lo3g_t", "-dbtype", "prot", "-taxids", "9606", "-no_taxid_expansion", "-target_only", "-outfmt", "%o %a"], 0),
    ("A109", &["-db", "db/ak/or/A_nest_fo3g_t", "-dbtype", "prot", "-info"], 0),
    ("A110", &["-db", "db/ak/or/A_nest_fo3g_t", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A111", &["-db", "db/ak/or/A_nest_fo3_lo3g_t", "-dbtype", "prot", "-info"], 0),
    ("A112", &["-db", "db/ak/or/A_nest_fo3_lo3g_t", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A113", &["-db", "db/ak/or/A_aa_lo3g_a db/ak/or/A_aa_lo3g_a", "-dbtype", "prot", "-info"], 0),
    ("A114", &["-db", "db/ak/or/A_aa_lo3g_a db/ak/or/A_aa_lo3g_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A115", &["-db", "db/ak/or/A_ab_lo3g_a db/ak/or/A_ab_lo3g_b", "-dbtype", "prot", "-info"], 0),
    ("A116", &["-db", "db/ak/or/A_ab_lo3g_a db/ak/or/A_ab_lo3g_b", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A117", &["-db", "db/ak/or/A_sub2_lo3g_a", "-dbtype", "prot", "-info"], 0),
    ("A118", &["-db", "db/ak/or/A_sub2_lo3g_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A119", &["-db", "db/ak/or/A_ao_lo3g_a db/ak/or/o", "-dbtype", "prot", "-info"], 0),
    ("A120", &["-db", "db/ak/or/A_ao_lo3g_a db/ak/or/o", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A121", &["-db", "db/ak/or/o db/ak/or/A_oa_lo3g_a", "-dbtype", "prot", "-info"], 0),
    ("A122", &["-db", "db/ak/or/o db/ak/or/A_oa_lo3g_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A123", &["-db", "db/ak/or/A_ao_fo3g_a db/ak/or/o", "-dbtype", "prot", "-info"], 0),
    ("A124", &["-db", "db/ak/or/A_ao_fo3g_a db/ak/or/o", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A125", &["-db", "db/ak/or/B_a_a", "-dbtype", "prot", "-info"], 0),
    ("A126", &["-db", "db/ak/or/B_a_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A127", &["-db", "db/ak/or/B_a_a", "-dbtype", "prot", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%o %a"], 1),
    ("A128", &["-db", "db/ak/or/B_a_a", "-dbtype", "prot", "-taxids", "9606", "-no_taxid_expansion", "-target_only", "-outfmt", "%o %a"], 1),
    ("A129", &["-db", "db/ak/or/B_aa_a db/ak/or/B_aa_a", "-dbtype", "prot", "-info"], 0),
    ("A130", &["-db", "db/ak/or/B_aa_a db/ak/or/B_aa_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A131", &["-db", "db/ak/or/B_aa_a db/ak/or/B_aa_a", "-dbtype", "prot", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%o %a"], 1),
    ("A132", &["-db", "db/ak/or/B_aa_a db/ak/or/B_aa_a", "-dbtype", "prot", "-taxids", "9606", "-no_taxid_expansion", "-target_only", "-outfmt", "%o %a"], 1),
    ("A133", &["-db", "db/ak/or/B_ao_a db/ak/or/o", "-dbtype", "prot", "-info"], 0),
    ("A134", &["-db", "db/ak/or/B_ao_a db/ak/or/o", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A135", &["-db", "db/ak/or/B_ao_a db/ak/or/o", "-dbtype", "prot", "-taxids", "10116", "-no_taxid_expansion", "-outfmt", "%o %a"], 0),
    ("A136", &["-db", "db/ak/or/B_ao_a db/ak/or/o", "-dbtype", "prot", "-taxids", "10116", "-no_taxid_expansion", "-target_only", "-outfmt", "%o %a"], 0),
    ("A137", &["-db", "db/ak/or/B_ao_a db/ak/or/o", "-dbtype", "prot", "-tax_info"], 0),
    ("A138", &["-db", "db/ak/or/B_ao_a db/ak/or/o", "-dbtype", "prot", "-entry", "o0", "-outfmt", "%o %a"], 1),
    ("A139", &["-db", "db/ak/or/o db/ak/or/B_oa_a", "-dbtype", "prot", "-info"], 0),
    ("A140", &["-db", "db/ak/or/o db/ak/or/B_oa_a", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A141", &["-db", "db/ak/or/o db/ak/or/B_oa_a", "-dbtype", "prot", "-taxids", "10116", "-no_taxid_expansion", "-outfmt", "%o %a"], 0),
    ("A142", &["-db", "db/ak/or/o db/ak/or/B_oa_a", "-dbtype", "prot", "-taxids", "10116", "-no_taxid_expansion", "-target_only", "-outfmt", "%o %a"], 0),
    ("A143", &["-db", "db/ak/or/o db/ak/or/B_oa_a", "-dbtype", "prot", "-entry_batch", "in/or_B_oa_b1.txt", "-outfmt", "%o %a"], 1),
    ("A144", &["-db", "db/ak/or/o db/ak/or/B_oa_a", "-dbtype", "prot", "-entry", "o0", "-outfmt", "%o %a"], 1),
    ("A145", &["-db", "db/ak/or/o db/ak/or/B_oa_a", "-dbtype", "prot", "-tax_info"], 0),
    ("A146", &["-db", "db/ak/or/B_ao_span_a db/ak/or/o", "-dbtype", "prot", "-info"], 0),
    ("A147", &["-db", "db/ak/or/B_ao_span_a db/ak/or/o", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A148", &["-db", "db/ak/or/B_ao_span_a db/ak/or/o", "-dbtype", "prot", "-taxids", "10116", "-no_taxid_expansion", "-outfmt", "%o %a"], 0),
    ("A149", &["-db", "db/ak/or/B_ao_span_a db/ak/or/o", "-dbtype", "prot", "-taxids", "10116", "-no_taxid_expansion", "-target_only", "-outfmt", "%o %a"], 0),
    ("A150", &["-db", "db/ak/or/B_ao_fo30_a db/ak/or/o", "-dbtype", "prot", "-info"], 0),
    ("A151", &["-db", "db/ak/or/B_ao_fo30_a db/ak/or/o", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A152", &["-db", "db/ak/or/B_ao_lo0_a db/ak/or/o", "-dbtype", "prot", "-info"], 0),
    ("A153", &["-db", "db/ak/or/B_ao_lo0_a db/ak/or/o", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A154", &["-db", "db/ak/or/B_ao_fo5_lo0_a db/ak/or/o", "-dbtype", "prot", "-info"], 0),
    ("A155", &["-db", "db/ak/or/B_ao_fo5_lo0_a db/ak/or/o", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A156", &["-db", "db/ak/or/B_ao_re_a db/ak/or/o", "-dbtype", "prot", "-info"], 0),
    ("A157", &["-db", "db/ak/or/B_ao_re_a db/ak/or/o", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A158", &["-db", "db/ak/or/B_ao_re_a db/ak/or/o", "-dbtype", "prot", "-taxids", "10116", "-no_taxid_expansion", "-outfmt", "%o %a"], 0),
    ("A159", &["-db", "db/ak/or/B_ao_re_a db/ak/or/o", "-dbtype", "prot", "-taxids", "10116", "-no_taxid_expansion", "-target_only", "-outfmt", "%o %a"], 0),
    ("A160", &["-db", "db/ak/or/B_ao_re_a db/ak/or/o", "-dbtype", "prot", "-tax_info"], 0),
    ("A161", &["-db", "db/ak/or/B_ao_re_a db/ak/or/o", "-dbtype", "prot", "-metadata"], 0),
    ("A162", &["-db", "db/ak/or/B_atfp_a db/ak/or/tfp", "-dbtype", "prot", "-info"], 0),
    ("A163", &["-db", "db/ak/or/B_atfp_a db/ak/or/tfp", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A164", &["-db", "db/ak/or/tfn", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A165", &["-db", "db/ak/or/Q3_B_a db/ak/or/tfn.01", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A166", &["-db", "db/ak/or/Q3_B_a db/ak/or/tfn.01", "-dbtype", "nucl", "-info"], 0),
    ("A167", &["-db", "db/ak/or/Q3_A_a", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A168", &["-db", "db/ak/or/Q3_A_a", "-dbtype", "nucl", "-info"], 0),
    ("A169", &["-db", "db/ak/or/Q3_Anest_t", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A170", &["-db", "db/ak/or/Q3_Anest_t", "-dbtype", "nucl", "-info"], 0),
    ("A171", &["-db", "db/ak/or/r2_re_a db/ak/or/tfp", "-dbtype", "prot", "-info"], 0),
    ("A172", &["-db", "db/ak/or/r2_re_a db/ak/or/tfp", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A173", &["-db", "db/ak/or/r2_re_a db/ak/or/tfp", "-dbtype", "prot", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%o %a"], 0),
    ("A174", &["-db", "db/ak/or/r2_re_a db/ak/or/tfp", "-dbtype", "prot", "-taxids", "9606", "-no_taxid_expansion", "-target_only", "-outfmt", "%o %a"], 0),
    ("A175", &["-db", "db/ak/or/r2_re_a db/ak/or/tfp", "-dbtype", "prot", "-tax_info"], 0),
    ("A176", &["-db", "db/ak/or/r2_re_a db/ak/or/tfp", "-dbtype", "prot", "-metadata"], 0),
    ("A177", &["-db", "db/ak/or/r2_re_a db/ak/or/tfp", "-dbtype", "prot", "-entry", "P69905.2", "-outfmt", "%o %a"], 1),
    ("A178", &["-db", "db/ak/or/r2_re_a db/ak/or/tfp", "-dbtype", "prot", "-entry_batch", "in/or_r2_re_b1.txt", "-outfmt", "%o %a"], 1),
    ("A179", &["-db", "db/ak/or/r2_ao_fo2147483648_a db/ak/or/o", "-dbtype", "prot", "-info"], 0),
    ("A180", &["-db", "db/ak/or/r2_ao_fo2147483648_a db/ak/or/o", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A181", &["-db", "db/ak/or/r2_ao_fo2147483649_a db/ak/or/o", "-dbtype", "prot", "-info"], 0),
    ("A182", &["-db", "db/ak/or/r2_ao_fo2147483649_a db/ak/or/o", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A183", &["-db", "db/ak/p_range", "-entry", "123456", "-outfmt", "%o %a"], 0),
    ("A184", &["-db", "db/ak/p_range", "-entry", "1818", "-outfmt", "%o %a"], 0),
    ("A185", &["-db", "db/ak/p_range", "-entry", "Q00003,123456,1818", "-outfmt", "%o %a"], 0),
    ("A186", &["-db", "db/ak/p_range", "-entry", "123456"], 0),
    ("A187", &["-db", "db/ak/p_range", "-entry_batch", "in/or_gi.txt", "-outfmt", "%o %a"], 0),
    ("A188", &["-db", "db/ak/p_range", "-entry_batch", "in/or_gi.txt", "-get_dups", "-outfmt", "%o %a"], 0),
    ("A189", &["-db", "db/ak/p_range", "-entry", "123456", "-outfmt", "%o %l"], 0),
    ("A190", &["-db", "db/p5", "-entry", "123456,1818", "-outfmt", "%o %a"], 0),
    // D1a-3 (`make_d1a3_fixtures.py`): `-info` of NCBI's `int` count and `Uint8` length with its
    // commas (`UInt8ToString` loses a digit from 4,294,967,296,000,000,000).
    ("A191", &["-db", "db/ak/cnt/negm1", "-dbtype", "prot", "-info"], 0),
    ("A192", &["-db", "db/ak/cnt/l4295", "-dbtype", "prot", "-info"], 0),
    ("A193", &["-db", "db/ak/cnt/nsm5", "-dbtype", "prot", "-info"], 0),
    ("A194", &["-db", "db/ak/cnt/neg1", "-dbtype", "prot", "-info"], 0),
    ("A195", &["-db", "db/ak/cnt/ns2p31", "-dbtype", "prot", "-info"], 0),
    // D1b-1: `OIDLIST` (`db/ak/ol/`, over `or/tfp` and `or/tfn`): the OID list file read relative
    // to each volume, NCBI's constructor scan of the OIDs it keeps unless a node above it has both
    // `NSEQ` and `LENGTH` (A202, A221: the walks), `-exact_length`'s exact scan (A211: the
    // approximate nucleotide lengths of A210), the bytes past the file (A215) and a header of
    // 0xFFFFFFFF (A213, A214), a range over the list in the only alias file of a `-db` (A216, A217),
    // and the errors of a missing file: at the constructor's scan (A207 and A208, which reports
    // the nucleotide type `-dbtype guess` tried first; A209 names the directory of the node's
    // volume), or at the first use of the OID list under `NSEQ` and `LENGTH` (A205, A206, A219,
    // A220; `-info` never reads it, A204).
    ("A196", &["-db", "db/ak/ol/p_k", "-info"], 0),
    ("A197", &["-db", "db/ak/ol/p_k", "-info", "-exact_length"], 0),
    ("A198", &["-db", "db/ak/ol/p_k", "-entry", "all", "-outfmt", "%o %a %l"], 0),
    ("A199", &["-db", "db/ak/ol/p_k", "-entry", "P01923.1", "-outfmt", "%o %a"], 1),
    ("A200", &["-db", "db/ak/ol/p_k", "-entry", "P69905.2,Q9TS35.2", "-outfmt", "%o %a"], 0),
    ("A201", &["-db", "db/ak/ol/p_k", "-metadata"], 0),
    ("A202", &["-db", "db/ak/ol/p_knl", "-info"], 0),
    ("A203", &["-db", "db/ak/ol/p_kr", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A204", &["-db", "db/ak/ol/p_missl", "-info"], 0),
    ("A205", &["-db", "db/ak/ol/p_missl", "-entry", "all", "-outfmt", "%o %a"], 1),
    ("A206", &["-db", "db/ak/ol/p_missl", "-entry", "P69905.2", "-outfmt", "%o %a"], 1),
    ("A207", &["-db", "db/ak/ol/p_miss", "-info", "-dbtype", "prot"], 2),
    ("A208", &["-db", "db/ak/ol/p_miss", "-info"], 2),
    ("A209", &["-db", "db/ak/ol/p_rel", "-info", "-dbtype", "prot"], 2),
    ("A210", &["-db", "db/ak/ol/n_n", "-info", "-dbtype", "nucl"], 0),
    ("A211", &["-db", "db/ak/ol/n_n", "-info", "-exact_length", "-dbtype", "nucl"], 0),
    ("A212", &["-db", "db/ak/ol/n_n", "-entry", "all", "-dbtype", "nucl", "-outfmt", "%o %a %l"], 0),
    ("A213", &["-db", "db/ak/ol/p_wrap", "-info"], 0),
    ("A214", &["-db", "db/ak/ol/p_wrap", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A215", &["-db", "db/ak/ol/p_h2", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A216", &["-db", "db/ak/ol/p_lr db/ak/or/tfp", "-info"], 0),
    ("A217", &["-db", "db/ak/ol/p_lr db/ak/or/tfp", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A218", &["-db", "db/ak/ol/p_k", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%o %a"], 0),
    ("A219", &["-db", "db/ak/ol/p_missl", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%o %a"], 1),
    ("A220", &["-db", "db/ak/ol/p_missl", "-metadata"], 2),
    ("A221", &["-db", "db/ak/ol/p_knl", "-info", "-exact_length"], 0),
    // Review A of line RG D1b-1 (2026-11-21, NCBI 2.17.0): an OID list over two range-only
    // sub-nodes keeps NCBI's OIDs (the union's bit array cut at the list's last OID where the
    // arrays differ in start or bytes, A222-A227; ANDed where they agree, A228); a range across
    // the volumes of tfn; a relative path through `../` as NCBI joins it; the list file beside
    // the first volume; `-ipg` reads a sequence without building the list (A232, A233) but not
    // its header (A234); a taxid list over a list file that fails at the scan.
    ("A222", &["-db", "db/ak/ol/p_cu2", "-entry", "all", "-outfmt", "%o|%a"], 0),
    ("A223", &["-db", "db/ak/ol/p_cu2", "-info"], 0),
    ("A224", &["-db", "db/ak/ol/p_cu8", "-entry", "all", "-outfmt", "%o|%a"], 0),
    ("A225", &["-db", "db/ak/ol/p_cu29", "-entry", "all", "-outfmt", "%o|%a"], 0),
    ("A226", &["-db", "db/ak/ol/p_cu0", "-info"], 0),
    ("A227", &["-db", "db/ak/ol/p_dq", "-entry", "all", "-outfmt", "%o|%a"], 0),
    ("A228", &["-db", "db/ak/ol/p_dr", "-entry", "all", "-outfmt", "%o|%a"], 0),
    ("A229", &["-db", "db/ak/ol/n_nr", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o|%a"], 0),
    ("A230", &["-db", "db/ak/ol/p_reln", "-info", "-dbtype", "prot"], 2),
    ("A231", &["-db", "db/ak/ol/p_fb", "-entry", "all", "-outfmt", "%o %a"], 0),
    ("A232", &["-db", "db/ak/ol/p_pl", "-dbtype", "prot", "-ipg", "1001", "-outfmt", "%s"], 0),
    ("A233", &["-db", "db/ak/ol/p_pl", "-dbtype", "prot", "-ipg", "1001", "-outfmt", "%h"], 0),
    ("A234", &["-db", "db/ak/ol/p_pl", "-dbtype", "prot", "-ipg", "1001", "-outfmt", "%s %a"], 1),
    ("A235", &["-db", "db/ak/ol/p_miss", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%o %a"], 2),
    // A range outside the bit array of the union below it: `-info` builds no list (A236); past
    // the array inside its last byte no OID is kept (A237). Before it and past its bytes,
    // `RANGE_REFUSALS`.
    ("A236", &["-db", "db/ak/ol/p_pb", "-info"], 0),
    ("A237", &["-db", "db/ak/ol/p_pa", "-entry", "all", "-outfmt", "%o|%a"], 0),
    // An OID list file of its header alone fails NCBI's assertion at the constructor's scan, a
    // `CSeqDBException` that NCBI's guess of the type catches to look for a nucleotide database.
    ("A239", &["-db", "db/ak/ol/p_h4", "-info"], 2),
    ("A240", &["-db", "db/ak/ol/p_h4", "-entry", "all", "-outfmt", "%o|%a"], 2),
    // OID_MASK_TYPE (D1b-2, `make_d1b2_blastdbcmd.py`) over ../db_search/oidmask: bit 0 erases
    // the lines with a predicted RefSeq id; a sequence left without a line has no line to print
    // (`%a`, `%t`, `%T`, FASTA), but its OID, length and letters are printed (`%o`, `%l`, `%s`);
    // an id erased still finds its sequence, whose other lines are printed, and `-target_only`
    // then finds no line of the id; the value converts as NCBI's unsigned int (A265, A269).
    ("A242", &["-db", "../db_search/oidmask/tall", "-entry", "all", "-outfmt", "%o|%a|%t"], 0),
    ("A243", &["-db", "../db_search/oidmask/ompx", "-entry", "all", "-outfmt", "%o|%a|%t"], 0),
    ("A244", &["-db", "../db_search/oidmask/ntall", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o|%a|%t"], 0),
    ("A245", &["-db", "../db_search/oidmask/omnx", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o|%a|%t"], 0),
    ("A246", &["-db", "../db_search/oidmask/tall", "-entry", "all"], 0),
    ("A247", &["-db", "../db_search/oidmask/omnx", "-dbtype", "nucl", "-entry", "all"], 0),
    ("A248", &["-db", "../db_search/oidmask/tall", "-entry", "all", "-outfmt", "%o"], 0),
    ("A249", &["-db", "../db_search/oidmask/tall", "-entry", "all", "-outfmt", "%o %l"], 0),
    ("A250", &["-db", "../db_search/oidmask/ompx", "-entry", "all", "-outfmt", "%o %s"], 0),
    ("A251", &["-db", "../db_search/oidmask/tall", "-entry", "all", "-outfmt", "%o %T"], 0),
    ("A252", &["-db", "../db_search/oidmask/tall", "-entry", "XP_000003.1", "-outfmt", "%o|%a|%t"], 0),
    ("A253", &["-db", "../db_search/oidmask/tall", "-entry", "XP_000003.1"], 0),
    ("A254", &["-db", "../db_search/oidmask/tall", "-entry", "XP_000002.1", "-outfmt", "%o|%a|%t"], 0),
    ("A255", &["-db", "../db_search/oidmask/ompx", "-entry", "XP_000002.1", "-outfmt", "%o|%a|%t"], 1),
    ("A256", &["-db", "../db_search/oidmask/tall", "-entry", "XP_000004.1", "-outfmt", "%o|%a|%t"], 0),
    ("A257", &["-db", "../db_search/oidmask/tall", "-entry", "XP_000004.1", "-target_only", "-outfmt", "%o|%a|%t"], 1),
    ("A258", &["-db", "../db_search/oidmask/tall", "-entry", "XP_000004.1", "-target_only"], 1),
    ("A259", &["-db", "../db_search/oidmask/tall", "-entry_batch", "in/om_batch.txt", "-outfmt", "%o|%a"], 0),
    ("A260", &["-db", "../db_search/oidmask/ompx", "-entry_batch", "in/om_batch.txt", "-outfmt", "%o|%a"], 1),
    ("A261", &["-db", "../db_search/oidmask/omp4x", "-entry", "all", "-outfmt", "%o|%a"], 0),
    ("A262", &["-db", "../db_search/oidmask/j5", "-entry", "all", "-outfmt", "%o|%a"], 0),
    ("A263", &["-db", "../db_search/oidmask/j2", "-entry", "all", "-outfmt", "%o|%a"], 0),
    ("A264", &["-db", "../db_search/oidmask/t2", "-entry", "all", "-outfmt", "%o|%a"], 0),
    ("A265", &["-db", "../db_search/oidmask/tabc", "-info"], 255),
    ("A266", &["-db", "../db_search/oidmask/tnol", "-entry", "all", "-outfmt", "%o|%a"], 0),
    ("A267", &["-db", "../db_search/oidmask/omnx", "-dbtype", "nucl", "-entry", "XM_000003.1,NM_000004.1", "-outfmt", "%o|%a|%t"], 0),
    ("A268", &["-db", "../db_search/oidmask/j1", "-entry", "all", "-outfmt", "%o|%a|%t"], 0),
    ("A269", &["-db", "../db_search/oidmask/tover", "-entry", "all", "-outfmt", "%o"], 255),
    ("A270", &["-db", "../db_search/oidmask/tall", "-entry", "all", "-outfmt", "%o|%a|%l|%s", "-line_length", "20"], 0),
    ("A272", &["-db", "../db_search/idlist/g_gi", "-info"], 0),
    ("A273", &["-db", "../db_search/idlist/g_gi", "-info", "-exact_length"], 0),
    ("A274", &["-db", "../db_search/idlist/g_gi", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A275", &["-db", "../db_search/idlist/g_gi", "-entry", "NP_000101.1", "-outfmt", "%o|%a|%t"], 0),
    ("A276", &["-db", "../db_search/idlist/g_gi", "-entry", "NP_000103.1", "-outfmt", "%o|%a"], 1),
    ("A277", &["-db", "../db_search/idlist/g_gi", "-entry", "P10102.1,AAA00108.1", "-outfmt", "%o|%a|%g"], 0),
    ("A278", &["-db", "../db_search/idlist/g_gil", "-info"], 0),
    ("A279", &["-db", "../db_search/idlist/g_gil", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A280", &["-db", "../db_search/idlist/g_gib", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A281", &["-db", "../db_search/idlist/g_nonl", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A282", &["-db", "../db_search/idlist/g_cr", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A283", &["-db", "../db_search/idlist/g_wrap", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A284", &["-db", "../db_search/idlist/g_unk", "-info"], 0),
    ("A285", &["-db", "../db_search/idlist/g_unk", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A286", &["-db", "../db_search/idlist/g_hash", "-info"], 2),
    ("A287", &["-db", "../db_search/idlist/g_hashl", "-info"], 0),
    ("A288", &["-db", "../db_search/idlist/g_hashl", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 1),
    ("A289", &["-db", "../db_search/idlist/g_alpha", "-info"], 2),
    ("A290", &["-db", "../db_search/idlist/g_nl1", "-info"], 2),
    ("A291", &["-db", "../db_search/idlist/g_empty", "-info"], 2),
    ("A292", &["-db", "../db_search/idlist/g_badcnt", "-info"], 2),
    ("A293", &["-db", "../db_search/idlist/g_short", "-info"], 2),
    ("A294", &["-db", "../db_search/idlist/g_miss", "-info"], 2),
    ("A295", &["-db", "../db_search/idlist/g_missl", "-info"], 0),
    ("A296", &["-db", "../db_search/idlist/g_missl", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 1),
    ("A297", &["-db", "../db_search/idlist/g_rel", "-info"], 2),
    ("A298", &["-db", "../db_search/idlist/g_miss", "-info", "-dbtype", "prot"], 2),
    ("A299", &["-db", "../db_search/idlist/g_multi", "-info"], 2),
    ("A300", &["-db", "../db_search/idlist/g_multil", "-info"], 2),
    ("A301", &["-db", "../db_search/idlist/g_multil", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 2),
    ("A302", &["-db", "../db_search/idlist/g_gr", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A303", &["-db", "../db_search/idlist/g_go", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A304", &["-db", "../db_search/idlist/g_nest", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A305", &["-db", "../db_search/idlist/g_nest", "-info"], 0),
    ("A306", &["-db", "../db_search/idlist/g_sib", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A307", &["-db", "../db_search/idlist/g_bare", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A308", &["-db", "../db_search/idlist/g_bare", "-info"], 0),
    ("A309", &["-db", "../db_search/idlist/g_bare", "-entry", "NP_000103.1", "-outfmt", "%o|%a"], 0),
    ("A310", &["-db", "../db_search/idlist/g_bare", "-entry", "NP_000101.1", "-outfmt", "%o|%a"], 0),
    ("A311", &["-db", "../db_search/idlist/g_ti", "-info"], 2),
    ("A312", &["-db", "../db_search/idlist/g_til", "-info"], 0),
    ("A313", &["-db", "../db_search/idlist/g_til", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 1),
    ("A314", &["-db", "../db_search/idlist/gt_ti", "-info"], 0),
    ("A315", &["-db", "../db_search/idlist/gt_ti", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A316", &["-db", "../db_search/idlist/gt_ti4", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A317", &["-db", "../db_search/idlist/gt_ti8", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A318", &["-db", "../db_search/idlist/gt_tib", "-info"], 2),
    ("A319", &["-db", "../db_search/idlist/gt_sib", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A320", &["-db", "../db_search/idlist/gt_sib", "-entry", "gnl|ti|1020", "-outfmt", "%o|%a"], 0),
    ("A321", &["-db", "../db_search/idlist/gu_gi", "-info"], 2),
    ("A322", &["-db", "../db_search/idlist/gu_gil", "-info"], 0),
    ("A323", &["-db", "../db_search/idlist/gu_gil", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 1),
    ("A324", &["-db", "../db_search/idlist/g4_sq", "-info"], 0),
    ("A325", &["-db", "../db_search/idlist/g4_sq", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A326", &["-db", "../db_search/idlist/g4_sqc", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A327", &["-db", "../db_search/idlist/g4_sqd", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A328", &["-db", "../db_search/idlist/g4_sqe", "-info"], 0),
    ("A329", &["-db", "../db_search/idlist/g4_sqe", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A330", &["-db", "../db_search/idlist/g4_sqm", "-info"], 2),
    ("A331", &["-db", "../db_search/idlist/g4_gs", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A332", &["-db", "../db_search/idlist/g4_gsu", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A333", &["-db", "../db_search/idlist/g4_gi", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A334", &["-db", "../db_search/idlist/g4_tx", "-info"], 0),
    ("A335", &["-db", "../db_search/idlist/g4_tx", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A336", &["-db", "../db_search/idlist/g4_txm", "-info"], 0),
    ("A337", &["-db", "../db_search/idlist/gm_gi", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A338", &["-db", "../db_search/idlist/gm_v0", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A339", &["-db", "../db_search/idlist/gm_v1", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A340", &["-db", "../db_search/idlist/gm_v01", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A341", &["-db", "../db_search/idlist/gm_mix", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A342", &["-db", "../db_search/idlist/gm_mix", "-info"], 0),
    ("A343", &["-db", "../db_search/idlist/gm_fr", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A344", &["-db", "../db_search/idlist/n_gi", "-dbtype", "nucl", "-info"], 0),
    ("A345", &["-db", "../db_search/idlist/n_gi", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A346", &["-db", "../db_search/idlist/n_gi", "-dbtype", "nucl", "-info", "-exact_length"], 0),
    ("A347", &["-db", "../db_search/idlist/n4_sq", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A348", &["-db", "../db_search/idlist/n4_gi", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 0),
    ("A349", &["-db", "../db_search/idlist/g_gi", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%o|%a|%T"], 0),
    ("A350", &["-db", "../db_search/idlist/g_gi", "-taxids", "10090", "-no_taxid_expansion", "-outfmt", "%o|%a|%T"], 0),
    ("A351", &["-db", "../db_search/idlist/g_bare", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%o|%a|%T"], 0),
    ("A352", &["-db", "../db_search/idlist/g_hash", "-dbtype", "prot", "-info"], 2),
    ("A353", &["-db", "../db_search/idlist/g_alpha", "-dbtype", "prot", "-info"], 2),
    ("A354", &["-db", "../db_search/idlist/g_nl1", "-dbtype", "prot", "-info"], 2),
    ("A355", &["-db", "../db_search/idlist/g_empty", "-dbtype", "prot", "-info"], 2),
    ("A356", &["-db", "../db_search/idlist/g_badcnt", "-dbtype", "prot", "-info"], 2),
    ("A357", &["-db", "../db_search/idlist/g_short", "-dbtype", "prot", "-info"], 2),
    ("A358", &["-db", "../db_search/idlist/g_rel", "-dbtype", "prot", "-info"], 2),
    ("A359", &["-db", "../db_search/idlist/g_multi", "-dbtype", "prot", "-info"], 2),
    ("A360", &["-db", "../db_search/idlist/g_multil", "-dbtype", "prot", "-info"], 2),
    ("A361", &["-db", "../db_search/idlist/g_ti", "-dbtype", "prot", "-info"], 2),
    ("A362", &["-db", "../db_search/idlist/gt_tib", "-dbtype", "prot", "-info"], 2),
    ("A363", &["-db", "../db_search/idlist/gu_gi", "-dbtype", "prot", "-info"], 2),
    ("A364", &["-db", "../db_search/idlist/g4_sqm", "-dbtype", "prot", "-info"], 2),
    ("A365", &["-db", "../db_search/idlist/g_multil", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 2),
    ("A366", &["-db", "../db_search/idlist/g_hashl", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%g|%T"], 1),
    ("A367", &["-db", "db/pzlg", "-dbtype", "prot", "-entry", "all", "-outfmt", "%a"], 0),
    ("A368", &["-db", "db/pzlgs", "-dbtype", "prot", "-entry", "all", "-outfmt", "%a"], 0),
    ("A369", &["-db", "db/pzlgs", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o"], 0),
    ("A370", &["-db", "db/pzlrg", "-dbtype", "prot", "-entry", "all", "-outfmt", "%a"], 1),
    ("A371", &["-db", "db/pzlrg", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o"], 0),
    ("A372", &["-db", "db/pzlrg", "-dbtype", "prot", "-entry", "all", "-outfmt", "%d"], 1),
    ("A373", &["-db", "db/pzlrg", "-dbtype", "prot", "-entry", "all"], 1),
    ("A374", &["-db", "db/pzlrg", "-dbtype", "prot", "-entry", "all", "-outfmt", "%b"], 1),
    ("A375", &["-db", "db/pzlrg", "-dbtype", "prot", "-entry", "5001", "-outfmt", "%o|%a"], 0),
    ("A376", &["-db", "db/pzlrn", "-dbtype", "prot", "-entry", "all", "-outfmt", "%a"], 1),
    ("A377", &["-db", "db/pzlrn", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o"], 0),
    ("A378", &["-db", "../db_search/idlist/tz_rg", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a"], 1),
    ("A379", &["-db", "../db_search/idlist/tz_rg", "-dbtype", "prot", "-entry", "all"], 1),
    ("A380", &["-db", "../db_search/idlist/tz_g", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a"], 0),
    ("A381", &["-db", "../db_search/idlist/tz_rg", "-dbtype", "prot", "-entry", "gnl|ti|77", "-outfmt", "%o|%a"], 0),
    ("A382", &["-db", "../db_search/idlist/pz_rg", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a"], 1),
    ("A383", &["-db", "../db_search/idlist/pz_og", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a"], 1),
    ("A384", &["-db", "../db_search/idlist/gtx_ti", "-dbtype", "prot", "-info"], 0),
    ("A385", &["-db", "../db_search/idlist/gtx_ti", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a"], 0),
    ("A386", &["-db", "../db_search/idlist/gd_gi", "-dbtype", "prot", "-info"], 2),
    ("A387", &["-db", "../db_search/idlist/gd_gi", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a"], 2),
    ("A388", &["-db", "../db_search/idlist/cv_sq", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a"], 0),
    ("A389", &["-db", "../db_search/idlist/cv_sq", "-dbtype", "prot", "-info"], 0),
];

/// A range outside the bit array of the union of the ranges below it, where NCBI clears bits
/// outside the array where it builds its OID list (line RG D1b-1), which LOSAT refuses there:
/// ending before the array (A238), NCBI crashes (SIGSEGV, probe D1a_P); starting past its bytes
/// (A241), NCBI writes over its heap, which crashes it unless the bytes stay in the unused end of
/// the array's allocation, as here (probes RB_pc and RB9 of review B of D1b-1).
#[rustfmt::skip]
const RANGE_REFUSALS: &[(&str, &[&str])] = &[
    ("A238", &["-db", "db/ak/ol/p_pb", "-entry", "all", "-outfmt", "%o|%a"]),
    ("A241", &["-db", "db/ak/ol/p_pc", "-entry", "all", "-outfmt", "%o|%a"]),
];

/// A definition line with bytes NCBI reads as '#' that an alias file's `OID_MASK_TYPE` erases
/// (`../db_search/oidmask/ombx`, line RG D1b-2, review A B1): NCBI decodes the whole header and
/// writes its error for each such byte before it erases the line (probe b1); LOSAT refuses it
/// as it refuses such bytes in a line it prints.
#[test]
fn removed_lines_with_bytes_ncbi_reads_as_hash_are_refused() {
    let config = TestDir::new("hash-refusals");
    let words = [
        "-db",
        "../db_search/oidmask/ombx",
        "-entry",
        "all",
        "-outfmt",
        "%o|%a|%t",
    ];
    let (stdout, stderr, code) = blastdbcmd(&words, &config.0);
    let stderr = String::from_utf8_lossy(&stderr);
    assert_eq!(
        (stdout.as_slice(), code),
        (b"".as_slice(), 1),
        "A271: {stderr}"
    );
    assert!(
        stderr.contains(
            "sequence 0 of the volume: a definition line that a list or the alias \
             files' OID_MASK_TYPE removes, with bytes NCBI reads as '#'"
        ),
        "A271: {stderr}"
    );
}

#[test]
fn ranges_outside_the_bit_array_below_them_are_refused() {
    let config = TestDir::new("range-refusals");
    for &(id, words) in RANGE_REFUSALS {
        let (stdout, stderr, code) = blastdbcmd(words, &config.0);
        let stderr = String::from_utf8_lossy(&stderr);
        assert_eq!(stdout, b"", "{id} {words:?}");
        assert_eq!(code, 1, "{id} {words:?}: {stderr}");
        assert!(
            stderr.contains("(NCBI crashes on it)"),
            "{id} {words:?}: {stderr}"
        );
    }
}

#[test]
fn alias_keys_are_ncbi_s_bytes() {
    let config = TestDir::new("cases");
    for &(id, words, status) in CASES {
        let (want_out, want_err) = (
            frozen(&format!("{id}.stdout")),
            frozen(&format!("{id}.stderr")),
        );
        // A path in the lines is written with the platform's separators, and `-metadata` cuts
        // a name at the platform's separator alone.
        if cfg!(windows)
            && (words.contains(&"-metadata") || holds_root(&want_out) || holds_root(&want_err))
        {
            continue;
        }
        let (stdout, stderr, code) = blastdbcmd(words, &config.0);
        assert_eq!(
            String::from_utf8_lossy(&stdout),
            String::from_utf8_lossy(&placed(&want_out)),
            "{id} {words:?}: stdout"
        );
        assert_eq!(
            String::from_utf8_lossy(&stderr),
            String::from_utf8_lossy(&placed(&want_err)),
            "{id} {words:?}: stderr"
        );
        assert_eq!(code, status, "{id} {words:?}: exit status");
    }
}
