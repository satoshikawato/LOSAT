//! `losat blastdbcmd`: the taxonomy (`-taxids`, `-taxidlist`, `-no_taxid_expansion`,
//! `-target_only`, `-tax_info`), the metadata (`-metadata`, `-metadata_output_prefix`), the
//! lists of databases (`-list`, `-list_outfmt`, `-recursive`, `-remove_redundant_dbs`),
//! `-show_blastdb_search_path` and the IPGs (`-ipg`) against NCBI BLAST+ 2.17.0 blastdbcmd. The
//! expected bytes are NCBI's, frozen by line BDC's probe J and K on ws (`tests/blastdbcmd/README.md`);
//! `{R}` stands for the directory of the probe's runs (here `tests/blastdbcmd`, the current
//! directory of each run). NCBI is not run here.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::ncbi_windows_paths::canonical;

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
            "losat_blastdbcmd_discovery_{stem}_{}_{nanos}",
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

/// One run of `losat blastdbcmd words...` in `cwd`: stdout, stderr, exit status. `BLASTDB` is
/// `blastdb` (a directory of `tests/blastdbcmd`) or unset; no registry file (`NCBI_CONFIG_PATH`
/// names `config` alone), no LOSAT variable, nothing on the standard input.
fn blastdbcmd(
    words: &[&str],
    cwd: &Path,
    blastdb: Option<&str>,
    config: &Path,
) -> (Vec<u8>, Vec<u8>, i32) {
    blastdbcmd_input(words, cwd, blastdb, config, None)
}

/// `blastdbcmd` with the file `stdin` (of `tests/blastdbcmd`) on the standard input.
fn blastdbcmd_input(
    words: &[&str],
    cwd: &Path,
    blastdb: Option<&str>,
    config: &Path,
    stdin: Option<&str>,
) -> (Vec<u8>, Vec<u8>, i32) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text == "BLASTDB" {
            command.env_remove(&name);
        }
    }
    if let Some(dir) = blastdb {
        command.env("BLASTDB", run_dir().join(dir));
    }
    let input = match stdin {
        Some(name) => std::process::Stdio::from(fs::File::open(run_dir().join(name)).unwrap()),
        None => std::process::Stdio::null(),
    };
    let output = command
        .env("NCBI_CONFIG_PATH", config)
        .stdin(input)
        .current_dir(cwd)
        .arg("blastdbcmd")
        .args(words)
        .output()
        .unwrap();
    (output.stdout, output.stderr, output.status.code().unwrap())
}

/// The frozen file `name` of `tests/blastdbcmd/expected`, with `{R}` placed.
fn expected(name: &str) -> Vec<u8> {
    let text = fs::read(run_dir().join("expected").join(name)).unwrap();
    let text = String::from_utf8(text).unwrap();
    // NCBI's footer in LOSAT's words (AGENTS.md rule 8).
    losat::wording::losat_text(&text.replace("{R}", run_dir().to_str().unwrap()))
        .into_owned()
        .into_bytes()
}

/// The probe's cases that LOSAT answers as NCBI does: id, the taxonomy directory of `BLASTDB`,
/// whether the output names a path (the lines of `-list`, the search path, the `dbname` of
/// `-metadata`, read with the platform's separator: Unix only), words, NCBI's exit status.
/// The databases are copies of the probe's (`j/lst`: version 4 and 5, without parsed ids, in
/// three volumes, an alias of two; `j/jd`: `p4` with other dates and a title to escape;
/// `j/jtax5`: a taxonomy over the taxids of `p5`). `K13` (`-taxids 1`, which never ends in NCBI:
/// taxid 1 is its own parent) is the approved difference of DI-8 DECISIONS 10, its expected
/// output NCBI's run of the explicit list of every taxid that taxid 1 reaches (K25).
/// A case: id, the taxonomy directory of `BLASTDB`, Unix only, words, NCBI's exit status.
type Case = (
    &'static str,
    Option<&'static str>,
    bool,
    &'static [&'static str],
    i32,
);

#[rustfmt::skip]
const CASES: &[Case] = &[
    ("J160", None, true, &["-list", "j/lst"], 0),
    ("J161", None, true, &["-list", "j/lst", "-list_outfmt", "%f|%p|%t|%d|%l|%n|%U|%v"], 0),
    ("J163", None, true, &["-list", "j/lst", "-dbtype", "prot", "-list_outfmt", "%f %t"], 0),
    ("J164", None, true, &["-list", "j/lst/", "-dbtype", "nucl"], 0),
    ("J169", None, false, &["-list", "j/lst", "-list_outfmt", "abc"], 1),
    ("J170", None, false, &["-list", "j/lst", "-list_outfmt", "%x"], 1),
    ("J172", None, false, &["-list", "j/lst", "-list_outfmt", "%f %"], 1),
    ("J173", None, false, &["-list", "nosuchdir", "-list_outfmt", "abc"], 1),
    ("J178", None, false, &["-list", "j/lst", "-list_outfmt", ""], 1),
    ("L06", None, false, &["-list", "nosuchdir"], 0),
    ("S38", None, false, &["-list", ""], 0),
    ("L07", None, true, &["-show_blastdb_search_path"], 0),
    ("J180", None, true, &["-db", "j/lst/p4", "-metadata"], 0),
    ("J181", None, true, &["-db", "j/lst/p5", "-metadata", "-metadata_output_prefix", "/x/y/"], 0),
    ("J184", None, true, &["-db", "j/lst/p5n", "-metadata"], 0),
    ("J185", None, true, &["-db", "j/jd/am12", "-metadata"], 0),
    ("J186", None, true, &["-db", "j/jd/pm12", "-metadata"], 0),
    ("J187", None, true, &["-db", "j/jd/d01", "-metadata"], 0),
    ("J188", None, true, &["-db", "j/jd/d1", "-metadata"], 0),
    ("J191", None, true, &["-db", "j/jd/esc", "-metadata"], 0),
    ("J192", None, true, &["-db", "j/lst/dv4", "-metadata"], 0),
    ("J197", None, true, &["-db", "j/lst/p5", "-metadata", "-metadata_output_prefix", ""], 0),
    ("J198", None, true, &["-db", "j/lst/pal", "-metadata"], 0),
    ("J146", None, false, &["-db", "j/lst/p5n", "-tax_info"], 0),
    ("J153", None, false, &["-db", "j/lst/p5", "-tax_info", "-outfmt", "%f"], 0),
    ("J155", None, false, &["-db", "j/lst/p5", "-tax_info", "-outfmt", "%n"], 0),
    ("J157", None, false, &["-db", "j/lst/n4", "-tax_info"], 2),
    ("J159", None, false, &["-db", "j/lst/p5", "-tax_info", "-outfmt", "%T %%"], 0),
    ("J143", None, false, &["-db", "j/lst/p5", "-tax_info", "-outfmt", "abc"], 1),
    ("T07", None, false, &["-db", "j/lst/p5", "-taxids", "9606", "-no_taxid_expansion", "-outfmt", "%o %a %T"], 0),
    ("T08", None, false, &["-db", "j/lst/p5n", "-taxids", "9606", "-no_taxid_expansion"], 1),
    ("S34", None, false, &["-db", "j/lst/p4", "-taxids", "9606", "-no_taxid_expansion"], 1),
    ("J106", None, false, &["-db", "j/lst/p5", "-taxids", "0", "-no_taxid_expansion", "-outfmt", "%o %a %T"], 0),
    ("J110", None, false, &["-db", "j/lst/p5", "-taxidlist", "in/nosuch.txt", "-no_taxid_expansion"], 1),
    ("J111", None, true, &["-db", "j/lst/p4", "-taxids", "9606", "-target_only"], 2),
    ("K09", None, false, &["-db", "j/lst/p5", "-taxids", "9606", "-no_taxid_expansion", "-dbtype", "prot", "-strand", "minus"], 0),
    ("K14", None, false, &["-db", "j/lst/p5", "-taxids", "9606", "-no_taxid_expansion", "-target_only"], 0),
    ("K17", None, false, &["-db", "j/lst/p5", "-taxids", "9606", "-target_only", "-dbtype", "prot"], 0),
    ("K10", Some("j/jtax5"), false, &["-db", "j/lst/p5", "-taxids", "9604", "-outfmt", "%o %a %T"], 0),
    ("K11", Some("j/jtax5"), false, &["-db", "j/lst/p5", "-taxids", "10066", "-target_only", "-outfmt", "%o %a %T"], 0),
    ("K12", Some("j/jtax5"), false, &["-db", "j/lst/p5", "-taxids", "9604,10066"], 0),
    ("K13", Some("j/jtax5"), false, &["-db", "j/lst/p5", "-taxids", "1", "-outfmt", "%o %T"], 0),
    // Line RG item T5: sequences without a title whose taxid `../seqdb/taxdb/t5` names
    // (`../db_search/t5/make_t5.py`): `-long_seqids` makes the title from the name (a protein's
    // `[name]`, without `unnamed protein product`), the short header keeps the empty title.
    ("T501", Some("../seqdb/taxdb/t5"), false, &["-db", "../db_search/t5/nclean", "-entry", "all", "-long_seqids"], 0),
    ("T502", Some("../seqdb/taxdb/t5"), false, &["-db", "../db_search/t5/pclean", "-entry", "all", "-long_seqids"], 0),
    ("T503", Some("../seqdb/taxdb/t5"), false, &["-db", "../db_search/t5/npall", "-entry", "all", "-long_seqids"], 0),
    ("T504", Some("../seqdb/taxdb/t5"), false, &["-db", "../db_search/t5/ppall", "-entry", "all", "-long_seqids"], 0),
    ("T505", Some("../seqdb/taxdb/t5"), false, &["-db", "../db_search/t5/pclean", "-entry", "all"], 0),
    ("T506", Some("../seqdb/taxdb/t5"), false, &["-db", "../db_search/t5/pclean", "-entry", "all", "-long_seqids", "-outfmt", "%f"], 0),
    ("T507", Some("../seqdb/taxdb/t5"), false, &["-db", "../db_search/t5/nclean", "-entry", "all", "-long_seqids", "-target_only"], 0),
    // Fix 1 of the review of T5 (`make_t5.py --mix-only`): several definition lines with other
    // taxids and titles, and every Seq-id kind; with `-target_only` the Bioseq keeps the target's
    // line for its ids and title, and its first source is the first positive taxid of every line,
    // the target's first (`x_GetTaxDefline`): of the other line where the target's has none
    // (T510, T511), the target's before the line stored first (T513, T514), and the suffix
    // `, partial [name]` rewritten from the other line's source (T512).
    ("T508", Some("../seqdb/taxdb/t5"), false, &["-db", "../db_search/t5/nmix", "-entry", "all", "-long_seqids"], 0),
    ("T509", Some("../seqdb/taxdb/t5"), false, &["-db", "../db_search/t5/pmix", "-entry", "all", "-long_seqids"], 0),
    ("T510", Some("../seqdb/taxdb/t5"), false, &["-db", "../db_search/t5/nmix", "-entry", "AB311001.1", "-target_only", "-long_seqids"], 0),
    ("T511", Some("../seqdb/taxdb/t5"), false, &["-db", "../db_search/t5/pmix", "-entry", "AB312001.1", "-target_only", "-long_seqids"], 0),
    ("T512", Some("../seqdb/taxdb/t5"), false, &["-db", "../db_search/t5/psuffix", "-entry", "AB362001.1", "-target_only", "-long_seqids"], 0),
    ("T513", Some("../seqdb/taxdb/t5"), false, &["-db", "../db_search/t5/nmix", "-entry", "AB301002.1", "-target_only", "-long_seqids"], 0),
    ("T514", Some("../seqdb/taxdb/t5"), false, &["-db", "../db_search/t5/pmix", "-entry", "AB302002.1", "-target_only", "-long_seqids"], 0),
    // Line RG item T5c (`make_t5.py --suffix-only`): the end of a protein title ending in `]`
    // rewritten from the Bioseq's first source (`x_AdjustProteinTitleSuffix`), with the taxonomy
    // database and without it (a source without a name drops a ` [ ]` end; ` []` stays), with
    // and without parsed Seq-ids (frozen from NCBI BLAST+ 2.17.0 blastdbcmd, 2026-11-14).
    ("T515", Some("../seqdb/taxdb/t5"), false, &["-db", "../db_search/t5/psfx", "-entry", "all", "-long_seqids"], 0),
    ("T516", None, false, &["-db", "../db_search/t5/psfx", "-entry", "all", "-long_seqids"], 0),
    ("T517", Some("../seqdb/taxdb/t5"), false, &["-db", "../db_search/t5/psfs", "-entry", "all", "-long_seqids"], 0),
    ("T518", None, false, &["-db", "../db_search/t5/psfs", "-entry", "all", "-long_seqids"], 0),
    // Fix 2 of line RG item T5c (`make_t5.py --taxid-only`): the taxid of the first record in an
    // encoding that NCBI's `ReadInt4` reads as makeblastdb's (the application INTEGER tag 0x42,
    // nine content bytes of sign extension; frozen from NCBI BLAST+ 2.17.0 blastdbcmd, 2026-11-15).
    ("T519", None, false, &["-db", "../db_search/t5/ptx_app42", "-entry", "all", "-outfmt", "%o %T %t"], 0),
    ("T520", None, false, &["-db", "../db_search/t5/ptx_sx9", "-entry", "all", "-outfmt", "%o %T %t"], 0),
    ("T521", None, false, &["-db", "../db_search/t5/ntx_app42", "-entry", "all", "-outfmt", "%o %T %t"], 0),
    ("T522", None, false, &["-db", "../db_search/t5/ntx_sx9", "-entry", "all", "-outfmt", "%o %T %t"], 0),
    // Merge 5 of line RG item T5c with item K4h (`make_t5.py --k5-only`): titles ending in `]` (`%t`
    // keeps the stored title) and scientific names with two-byte UTF-8 and a control byte, written
    // as they are (frozen from NCBI BLAST+ 2.17.0 blastdbcmd, 2026-11-17).
    ("T525", Some("../seqdb/taxdb/k5"), false, &["-db", "../db_search/t5/pk5", "-entry", "all", "-outfmt", "%t %T %S"], 0),
    ("T526", Some("../seqdb/taxdb/k5"), false, &["-db", "../db_search/t5/pk5s", "-entry", "all", "-outfmt", "%t %T %S"], 0),
    // Line RG item T5b: a damaged taxonomy database (`../seqdb/taxdb/<variant>`) on the search path
    // of `../db_search/t5b` (make_t5b.py there): NCBI posts the texts of opening it once, at the
    // first lookup, and the offset error of its last entry at each lookup that asks for taxid
    // 99990 (the FASTA formatter's two `GetBioseq` calls follow the volume's 256-slot cache;
    // `GetTaxInfo` of `%S` and `-tax_info` look up directly), and goes on (exit status 0). In
    // `btd_empty` the lookups of 99990 read no record (no crash). T5B501-T5B504: `-tax_info` of a
    // taxid found in the index without names writes empty names, `N/A` only for one not found.
    // NCBI BLAST+ 2.17.0 (local oracle), 2026-11-14, ad-RG-T5-4.
    ("T5B101", Some("../seqdb/taxdb/bad_magic"), false, &["-db", "../db_search/t5b/nmv", "-entry", "all"], 0),
    ("T5B102", Some("../seqdb/taxdb/bad_magic"), false, &["-db", "../db_search/t5b/nmv", "-entry", "all", "-long_seqids"], 0),
    ("T5B103", Some("../seqdb/taxdb/bad_magic"), false, &["-db", "../db_search/t5b/pmv", "-entry", "all", "-long_seqids"], 0),
    ("T5B104", Some("../seqdb/taxdb/bad_magic"), false, &["-db", "../db_search/t5b/nmv", "-entry", "all", "-outfmt", "%a %T %S"], 0),
    ("T5B105", Some("../seqdb/taxdb/bad_magic"), false, &["-db", "../db_search/t5b/nmv", "-entry", "all", "-outfmt", "%a %C %N"], 0),
    ("T5B106", Some("../seqdb/taxdb/bad_magic"), false, &["-db", "../db_search/t5b/nmv", "-entry", "AB410003.1,AB410002.1"], 0),
    ("T5B107", Some("../seqdb/taxdb/bad_magic"), false, &["-db", "../db_search/t5b/nmv", "-entry", "AB410008.1", "-target_only"], 0),
    ("T5B108", Some("../seqdb/taxdb/bad_magic"), false, &["-db", "../db_search/t5b/nmv", "-entry_batch", "../db_search/t5b/batch.txt"], 0),
    ("T5B109", Some("../seqdb/taxdb/bad_magic"), false, &["-db", "../db_search/t5b/nmv", "-taxids", "9606", "-no_taxid_expansion"], 0),
    ("T5B110", Some("../seqdb/taxdb/bad_magic"), false, &["-db", "../db_search/t5b/nmv", "-tax_info", "-outfmt", "%T %S"], 0),
    ("T5B201", Some("../seqdb/taxdb/count_zero"), false, &["-db", "../db_search/t5b/nmv", "-entry", "all"], 0),
    ("T5B202", Some("../seqdb/taxdb/count_zero"), false, &["-db", "../db_search/t5b/nmv", "-entry", "all", "-long_seqids"], 0),
    ("T5B203", Some("../seqdb/taxdb/count_zero"), false, &["-db", "../db_search/t5b/pmv", "-entry", "all", "-long_seqids"], 0),
    ("T5B204", Some("../seqdb/taxdb/count_zero"), false, &["-db", "../db_search/t5b/nmv", "-entry", "all", "-outfmt", "%a %T %S"], 0),
    ("T5B205", Some("../seqdb/taxdb/count_zero"), false, &["-db", "../db_search/t5b/nmv", "-entry", "all", "-outfmt", "%a %C %N"], 0),
    ("T5B206", Some("../seqdb/taxdb/count_zero"), false, &["-db", "../db_search/t5b/nmv", "-entry", "AB410003.1,AB410002.1"], 0),
    ("T5B207", Some("../seqdb/taxdb/count_zero"), false, &["-db", "../db_search/t5b/nmv", "-entry", "AB410008.1", "-target_only"], 0),
    ("T5B208", Some("../seqdb/taxdb/count_zero"), false, &["-db", "../db_search/t5b/nmv", "-entry_batch", "../db_search/t5b/batch.txt"], 0),
    ("T5B209", Some("../seqdb/taxdb/count_zero"), false, &["-db", "../db_search/t5b/nmv", "-taxids", "9606", "-no_taxid_expansion"], 0),
    ("T5B210", Some("../seqdb/taxdb/count_zero"), false, &["-db", "../db_search/t5b/nmv", "-tax_info", "-outfmt", "%T %S"], 0),
    ("T5B301", Some("../seqdb/taxdb/offset_last_past_end"), false, &["-db", "../db_search/t5b/nmv", "-entry", "all"], 0),
    ("T5B302", Some("../seqdb/taxdb/offset_last_past_end"), false, &["-db", "../db_search/t5b/nmv", "-entry", "all", "-long_seqids"], 0),
    ("T5B303", Some("../seqdb/taxdb/offset_last_past_end"), false, &["-db", "../db_search/t5b/pmv", "-entry", "all", "-long_seqids"], 0),
    ("T5B304", Some("../seqdb/taxdb/offset_last_past_end"), false, &["-db", "../db_search/t5b/nmv", "-entry", "all", "-outfmt", "%a %T %S"], 0),
    ("T5B305", Some("../seqdb/taxdb/offset_last_past_end"), false, &["-db", "../db_search/t5b/nmv", "-entry", "all", "-outfmt", "%a %C %N"], 0),
    ("T5B306", Some("../seqdb/taxdb/offset_last_past_end"), false, &["-db", "../db_search/t5b/nmv", "-entry", "AB410003.1,AB410002.1"], 0),
    ("T5B307", Some("../seqdb/taxdb/offset_last_past_end"), false, &["-db", "../db_search/t5b/nmv", "-entry", "AB410008.1", "-target_only"], 0),
    ("T5B308", Some("../seqdb/taxdb/offset_last_past_end"), false, &["-db", "../db_search/t5b/nmv", "-entry_batch", "../db_search/t5b/batch.txt"], 0),
    ("T5B309", Some("../seqdb/taxdb/offset_last_past_end"), false, &["-db", "../db_search/t5b/nmv", "-taxids", "9606", "-no_taxid_expansion"], 0),
    ("T5B310", Some("../seqdb/taxdb/offset_last_past_end"), false, &["-db", "../db_search/t5b/nmv", "-tax_info", "-outfmt", "%T %S"], 0),
    ("T5B401", Some("../seqdb/taxdb/btd_empty"), false, &["-db", "../db_search/t5b/nmv", "-entry", "AB410008.1", "-target_only"], 0),
    ("T5B402", Some("../seqdb/taxdb/btd_empty"), false, &["-db", "../db_search/t5b/nmv", "-entry_batch", "../db_search/t5b/batch.txt"], 0),
    ("T5B501", Some("../seqdb/taxdb/tabs0"), false, &["-db", "../db_search/t5b/nmv", "-tax_info", "-outfmt", "%T|%S|%L|%K|%B"], 0),
    ("T5B502", Some("../seqdb/taxdb/tabs1"), false, &["-db", "../db_search/t5b/nmv", "-tax_info", "-outfmt", "%T|%S|%L|%K|%B"], 0),
    ("T5B503", Some("../seqdb/taxdb/offset_last_at_end"), false, &["-db", "../db_search/t5b/nmv", "-tax_info", "-outfmt", "%T|%S|%L|%K|%B"], 0),
    ("T5B504", Some("../seqdb/taxdb/empty_fields"), false, &["-db", "../db_search/t5b/nmv", "-tax_info", "-outfmt", "%T|%S|%L|%K|%B"], 0),
    // Line RG item T5b fix 2: a record of `taxdb.btd` that reaches past NCBI's mapping of the file
    // (`../db_search/t5b/taxdb/<variant>`, make_t5b.py --taxdb). `kingdom_past_page`: 9606's
    // record runs past the mapping after its third tab, so its names are NCBI's and only its
    // super kingdom comes from memory past the mapping (`%K`, refused in taxdb_tool_texts.rs);
    // `names_past_page`: 9606's record lies wholly past it, so only the runs that write none of
    // its names go on (the loads of titled sequences, `%a %T`). NCBI BLAST+ 2.17.0 (local oracle),
    // 2026-11-15, the T5b fix 2 subagent of ad-RG-T5-5.
    ("T5B601", Some("../db_search/t5b/taxdb/kingdom_past_page"), false, &["-db", "../db_search/t5b/nmv", "-entry", "all"], 0),
    ("T5B602", Some("../db_search/t5b/taxdb/kingdom_past_page"), false, &["-db", "../db_search/t5b/pmv", "-entry", "all", "-long_seqids"], 0),
    ("T5B603", Some("../db_search/t5b/taxdb/kingdom_past_page"), false, &["-db", "../db_search/t5b/nmv", "-entry", "all", "-outfmt", "%a %T %S %L %B"], 0),
    ("T5B604", Some("../db_search/t5b/taxdb/kingdom_past_page"), false, &["-db", "../db_search/t5b/nmv", "-entry", "all", "-outfmt", "%a %C %N"], 0),
    ("T5B605", Some("../db_search/t5b/taxdb/kingdom_past_page"), false, &["-db", "../db_search/t5b/nmv", "-tax_info", "-outfmt", "%T|%S|%L|%B"], 0),
    ("T5B606", Some("../db_search/t5b/taxdb/names_past_page"), false, &["-db", "../db_search/t5b/nmv", "-entry", "AB410002.1,AB410011.1"], 0),
    ("T5B607", Some("../db_search/t5b/taxdb/names_past_page"), false, &["-db", "../db_search/t5b/nmv", "-entry", "all", "-outfmt", "%a %T"], 0),
    ("T5B608", Some("../db_search/t5b/taxdb/names_past_page"), false, &["-db", "../db_search/t5b/nmv", "-entry_batch", "../db_search/t5b/batch.txt"], 0),
    // Line RG item T5b fix 3: `zero_kingdom_past_page` (make_t5b.py --taxdb), whose record of taxid 0
    // runs past NCBI's mapping of taxdb.btd after its third tab: `-tax_info` writes `N/A` for taxid 0
    // whatever its lookup reads (blastdbcmd.cpp:711-718), also in a format with the super kingdom.
    // NCBI BLAST+ 2.17.0 (local oracle), 2026-11-15, the T5b fix 3 subagent of ad-RG-T5-5.
    ("T5B609", Some("../db_search/t5b/taxdb/zero_kingdom_past_page"), false, &["-db", "../db_search/t5b/nmv", "-tax_info", "-outfmt", "%T|%K"], 0),
    ("T5B610", Some("../db_search/t5b/taxdb/zero_kingdom_past_page"), false, &["-db", "../db_search/t5b/nmv", "-tax_info", "-outfmt", "%f"], 0),
];

/// The probe's cases that LOSAT refuses (one line naming LOSAT's BLASTDBCMD, exit status 1,
/// nothing written): NCBI's exception that names its source files (`-taxids` with an empty id,
/// a date `CTime` does not read), NCBI's protein letters written reversed (a taxid list on a
/// protein database read with `-dbtype guess` and `-strand minus`), the OID that NCBI's
/// `-ipg` never sets on a database without PIG files, and a first record whose taxid is outside
/// an `int` (`make_t5.py --taxid-only`: NCBI's `error while reading taxid`, exit 1, at that record).
#[rustfmt::skip]
const REFUSED: &[(&str, &[&str])] = &[
    ("J101", &["-db", "j/lst/p5", "-taxids", "4840,,4841", "-no_taxid_expansion", "-outfmt", "%a %T"]),
    ("J114", &["-db", "j/lst/p5", "-taxids", "9606", "-no_taxid_expansion", "-strand", "minus"]),
    ("J189", &["-db", "j/jd/feb30", "-metadata"]),
    ("J190", &["-db", "j/jd/bad", "-metadata"]),
    ("J204", &["-db", "j/lst/p5", "-ipg", "5"]),
    ("T523", &["-db", "../db_search/t5/ptx_big5", "-entry", "all", "-outfmt", "%o %T %t"]),
    ("T524", &["-db", "../db_search/t5/ntx_big5", "-entry", "all", "-outfmt", "%o %T %t"]),
];

#[test]
fn taxonomy_metadata_and_lists_are_ncbi_s_bytes() {
    let config = TestDir::new("cases");
    for &(id, blastdb, unix_only, words, status) in CASES {
        if unix_only && cfg!(windows) {
            continue;
        }
        let (stdout, stderr, code) = blastdbcmd(words, &run_dir(), blastdb, &config.0);
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
}

#[test]
fn what_ncbi_does_not_answer_is_refused() {
    let config = TestDir::new("refused");
    for &(id, words) in REFUSED {
        let (stdout, stderr, code) = blastdbcmd(words, &run_dir(), None, &config.0);
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

/// The entries of the probe's `db/j/lsx`: files the walk counts without opening them, a
/// directory named like a database, links to a file, to nowhere and to a directory, names that
/// sort apart from their bytes once quoted.
#[cfg(unix)]
fn make_lsx(dir: &Path) {
    use std::os::unix::fs::symlink;
    fs::create_dir_all(dir.join("sub").join("deep")).unwrap();
    fs::create_dir_all(dir.join("d.pin")).unwrap();
    for name in [
        "a b.pin", "a!.nin", ".h.pin", "X.PIN", "zz.pin", "n5.nin", "n5x.nin", "n55.nin",
    ] {
        fs::write(dir.join(name), b"").unwrap();
    }
    fs::write(dir.join("zz.pal"), b"DBLIST zz\n").unwrap();
    fs::write(dir.join("q.nal"), b"TITLE q\nDBLIST n5 \"n5x\"\n").unwrap();
    fs::write(dir.join("sub").join("y.nin"), b"").unwrap();
    fs::write(dir.join("sub").join("deep").join("z.pin"), b"").unwrap();
    fs::write(dir.join("d.pin").join("inner.pin"), b"").unwrap();
    symlink("zz.pin", dir.join("ln.pin")).unwrap();
    symlink("nowhere.nin", dir.join("br.nin")).unwrap();
    symlink("sub", dir.join("lsub")).unwrap();
}

// Probes J165-J167 and J176 of line BDC on a copy of `db/j/lsx` built here.
#[cfg(unix)]
#[test]
fn the_walk_lists_files_links_and_directories_as_ncbi_s() {
    let work = TestDir::new("lsx");
    let lsx = work.0.join("lsx");
    make_lsx(&lsx);
    let cases: [(&str, &[&str]); 4] = [
        ("J165", &["-list", "lsx"]),
        ("J166", &["-list", "lsx", "-recursive"]),
        ("J167", &["-list", "lsx", "-remove_redundant_dbs"]),
        (
            "J176",
            &[
                "-list",
                "lsx",
                "-recursive",
                "-dbtype",
                "nucl",
                "-remove_redundant_dbs",
                "-list_outfmt",
                "%p:%f",
            ],
        ),
    ];
    for (id, words) in cases {
        let (stdout, stderr, code) = blastdbcmd(words, &work.0, None, &work.0);
        let frozen = fs::read(run_dir().join("expected").join(format!("{id}.stdout"))).unwrap();
        let want = String::from_utf8(frozen)
            .unwrap()
            .replace("db/j/lsx", "lsx");
        assert_eq!(String::from_utf8_lossy(&stdout), want, "{id} {words:?}");
        assert!(
            stderr.is_empty(),
            "{id}: {}",
            String::from_utf8_lossy(&stderr)
        );
        assert_eq!(code, 0, "{id}");
    }
}

/// A run (stdout, stderr, exit status) against the frozen `id` of `tests/blastdbcmd/expected`.
fn check(id: &str, (stdout, stderr, code): (Vec<u8>, Vec<u8>, i32), status: i32) {
    assert_eq!(
        String::from_utf8_lossy(&stdout),
        String::from_utf8_lossy(&expected(&format!("{id}.stdout"))),
        "{id}: stdout"
    );
    assert_eq!(
        String::from_utf8_lossy(&stderr),
        String::from_utf8_lossy(&expected(&format!("{id}.stderr"))),
        "{id}: stderr"
    );
    assert_eq!(code, status, "{id}: exit status");
}

// Probe K3 of line BDC: an input file `-` is the standard input (`CArg_InputFile::x_Open`), an
// empty name opens nothing, a directory opens and reads nothing (Linux's stream), and `-list ""`
// lists the current directory without a prefix (`CDir::GetEntriesPtr`).
#[test]
fn input_files_and_the_current_directory_are_read_as_ncbi_s() {
    let config = TestDir::new("inputs");
    let words = ["-db", "db/p4", "-entry_batch", "-", "-outfmt", "%o %a"];
    let run = blastdbcmd_input(&words, &run_dir(), None, &config.0, Some("in/i_p_safe.txt"));
    check("K32", run, 1);
    let words = ["-db", "j/lst/p5", "-taxidlist", "", "-no_taxid_expansion"];
    check("K39", blastdbcmd(&words, &run_dir(), None, &config.0), 1);
    if cfg!(unix) {
        let words = ["-db", "db/p5", "-entry_batch", "in"];
        check("K35", blastdbcmd(&words, &run_dir(), None, &config.0), 0);
        let lst = run_dir().join("j").join("lst");
        check("K36", blastdbcmd(&["-list", ""], &lst, None, &config.0), 0);
    }
}
