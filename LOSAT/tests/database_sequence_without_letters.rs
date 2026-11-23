//! A nucleotide database whose first sequence has no letters (its offsets patched so, as no
//! makeblastdb writes it; line RG item D5a-t, 2026-11-21, review 3 of `remaining-gaps-20261112/D/d5`).
//! NCBI's `CSeqDB::GetBioseq` gives no Bioseq for it, before its taxonomy lookups. NCBI BLAST+
//! 2.17.0 windowmasker skips it with an id list matched as text: with an empty `taxdb.btd` whose
//! index holds the sequence's taxid it masks the other sequence and exits 0 (a lookup there would
//! read the empty file), and with a `taxdb.bti` of a wrong magic number it writes its error once,
//! at the other sequence's lookup. NCBI's makembindex dereferences the null Bioseq
//! (`CCoreException::eNullPtr`, exit status 3, a stack trace of the run's addresses): LOSAT refuses
//! (PD-LOSAT-NCBI-DEFECTS 1.61). NCBI's bytes: probe `D/d5/probe/zero` (local official release).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// NCBI windowmasker's standard output for both databases: the masks of `a1` (letters 20,000 to
/// 21,999 of `tests/winmask/fasta/dm2L_60k.fa`).
const MASKS: &str = ">titled\n20 - 74\n79 - 93\n111 - 125\n192 - 206\n345 - 382\n407 - 427\n\
464 - 485\n536 - 562\n596 - 610\n643 - 657\n886 - 900\n1117 - 1132\n1211 - 1225\n1338 - 1351\n\
1398 - 1412\n1486 - 1500\n1539 - 1562\n1821 - 1856\n1885 - 1913\n";

fn manifest() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

/// `losat <words>` in `dir` without LOSAT, NCBI or BLAST variables, with `BLASTDB` = `blastdb`.
fn losat(dir: &Path, blastdb: Option<&Path>, words: &[&str]) -> Output {
    let mut command = Command::new(BIN);
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_") || text.starts_with("NCBI_") || text.starts_with("BLAST") {
            command.env_remove(&name);
        }
    }
    if let Some(blastdb) = blastdb {
        command.env("BLASTDB", blastdb);
    }
    command
        .env("NCBI_CONFIG_PATH", dir)
        .current_dir(dir)
        .args(words)
        .output()
        .unwrap()
}

/// A directory with the database `z` (`lcl|e1` of 4 letters, taxid 9606, patched to 0 letters;
/// `lcl|a1 titled`, taxid `a1_taxid`), made by LOSAT's makeblastdb, and the id list `ids.txt`.
fn database(tag: &str, a1_taxid: u32) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("without-letters-{tag}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let dm: String = fs::read_to_string(manifest().join("tests/winmask/fasta/dm2L_60k.fa"))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('>'))
        .map(str::trim)
        .collect();
    let a1 = &dm[20000..22000];
    let mut fasta = String::from(">lcl|e1\nACGT\n>lcl|a1 titled\n");
    for line in a1.as_bytes().chunks(60) {
        fasta.push_str(std::str::from_utf8(line).unwrap());
        fasta.push('\n');
    }
    fs::write(dir.join("z.fa"), fasta).unwrap();
    fs::write(
        dir.join("tx.txt"),
        format!("lcl|e1\t9606\nlcl|a1\t{a1_taxid}\n"),
    )
    .unwrap();
    fs::write(dir.join("ids.txt"), "titled\n").unwrap();
    let made = losat(
        &dir,
        None,
        &[
            "makeblastdb",
            "-in",
            "z.fa",
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-taxid_map",
            "tx.txt",
            "-out",
            "z",
        ],
    );
    assert!(made.status.success(), "{made:?}");
    // The `.nin` ends with the header, sequence and ambiguity offsets of the 2 sequences (3 big
    // endian words each): the first sequence starts one byte later, at its last byte, whose
    // remainder is 0, so it has 0 letters.
    let nin = dir.join("z.nin");
    let mut bytes = fs::read(&nin).unwrap();
    let at = bytes.len() - 2 * 3 * 4;
    let start = u32::from_be_bytes(bytes[at..at + 4].try_into().unwrap());
    assert_eq!(start, 1);
    bytes[at..at + 4].copy_from_slice(&(start + 1).to_be_bytes());
    fs::write(&nin, bytes).unwrap();
    dir
}

fn taxdb(variant: &str) -> PathBuf {
    manifest().join("tests/seqdb/taxdb").join(variant)
}

fn windowmasker(dir: &Path, taxdb: &Path) -> Output {
    let ustat = manifest().join("tests/winmask/counts/dm2L_u9.binary");
    losat(
        dir,
        Some(taxdb),
        &[
            "windowmasker",
            "-ustat",
            ustat.to_str().unwrap(),
            "-infmt",
            "blastdb",
            "-in",
            "z",
            "-ids",
            "ids.txt",
        ],
    )
}

#[test]
fn windowmasker_skips_a_database_sequence_without_letters_before_its_lookup() {
    let dir = database("wm-empty-btd", 0);
    let run = windowmasker(&dir, &taxdb("btd_empty"));
    assert_eq!(
        (
            String::from_utf8_lossy(&run.stdout).as_ref(),
            String::from_utf8_lossy(&run.stderr).as_ref(),
            run.status.code()
        ),
        (MASKS, "", Some(0))
    );
    let dir = database("wm-bad-magic", 9606);
    let run = windowmasker(&dir, &taxdb("bad_magic"));
    assert_eq!(
        (
            String::from_utf8_lossy(&run.stdout).as_ref(),
            String::from_utf8_lossy(&run.stderr).as_ref(),
            run.status.code()
        ),
        (
            MASKS,
            "Error: Error: Tax database file has wrong magic number.\n",
            Some(0)
        )
    );
}

#[test]
fn makembindex_refuses_a_database_sequence_without_letters() {
    let dir = database("mk", 0);
    let expected = "creating @D@/z.00.idx...\nError: sequence 0 of the database, which has no \
                    letters or no definition line left by the alias files' OID_MASK_TYPE (NCBI's \
                    makembindex reads the null sequence its database reader gives for it and \
                    stops with a report of the run's memory addresses), is not supported by \
                    LOSAT's MAKEMBINDEX\n";
    // The index's name has the separator of the platform after the directory.
    let expected = if cfg!(windows) {
        crate::ncbi_windows_paths::ncbi_paths_on_windows(expected, &["@D@"])
    } else {
        expected.to_string()
    };
    for idmap in ["false", "true"] {
        let run = losat(&dir, None, &["makembindex", "-input", "z", "-idmap", idmap]);
        let err = String::from_utf8_lossy(&run.stderr).replace(&*dir.to_string_lossy(), "@D@");
        assert_eq!(
            (err.as_str(), run.status.code(), run.stdout.len()),
            (expected.as_str(), Some(1), 0),
            "-idmap {idmap}"
        );
    }
}
