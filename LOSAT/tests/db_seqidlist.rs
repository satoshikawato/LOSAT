//! `-seqidlist` with `-db` outside the manifest rows (line DI-8b, DI8 DECISIONS 17 and 20):
//! a list file that is not on the search path, which NCBI BLAST+ 2.17.0 reports with the path
//! of its build tree (exit status 255), and a list that keeps a sequence without a definition
//! line (a PDB molecule named without its chain), whose hits NCBI leaves out of its reports
//! without notice and which LOSAT refuses. NCBI is not run here; its texts are those of the
//! probe of dropin-survey-20261029/DI8 (S22, R01-R07).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn engine_dir() -> PathBuf {
    fs::canonicalize(Path::new(env!("CARGO_MANIFEST_DIR"))).unwrap()
}

/// One `losat PROGRAM -db DB -query QUERY words... -outfmt 6 -num_threads 1` in LOSAT/, with
/// no registry file, no `BLASTDB` and no other NCBI, BLAST or LOSAT variable. Returns stdout,
/// stderr and the exit status.
fn search(program: &str, db: &str, query: &str, words: &[&str]) -> (Vec<u8>, Vec<u8>, i32) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let name = name.to_string_lossy();
        if name.starts_with("LOSAT_")
            || name.starts_with("NCBI")
            || name.starts_with("BLAST")
            || name == "HOME"
        {
            command.env_remove(name.as_ref());
        }
    }
    let output = command
        .current_dir(engine_dir())
        .env("NCBI_CONFIG_PATH", engine_dir().join("tests/db_search/seqidlist"))
        .args([program, "-db", db, "-query", query])
        .args(words)
        .args(["-outfmt", "6", "-num_threads", "1"])
        .output()
        .unwrap();
    (output.stdout, output.stderr, output.status.code().unwrap_or(-1))
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:231-234
// ```c++
//     catch (const CException& e) {                                           \
//         LOG_POST(Error << "Error: " << e.what());                           \
//         exit_code = BLAST_UNKNOWN_ERROR;                                    \
//     }                                                                       \
// ```
// NCBI 2.17.0 (the conda build) writes `Error: NCBI C++ Exception:\n    T0 "/opt/conda/conda-bld/
// blast_1787399212531/work/c++/src/corelib/ncbifile.cpp", line 5754: Error:
// (CFileException::eMemoryMap) ncbi::CMemoryFileMap::CMemoryFileMap() - To be memory mapped the
// file must exist: ''\n\n` (probe S22); LOSAT names the build tree of the official Linux release
// and the test does not compare the path.
#[test]
fn a_list_file_off_the_search_path_is_ncbis_memory_map_error() {
    for program in ["blastp", "blastn"] {
        let (db, query) = match program {
            "blastp" => ("tests/db_search/taxfilter/tfp", "tests/db_search/taxfilter/q.faa"),
            _ => ("tests/db_search/taxfilter/tfn", "tests/db_search/taxfilter/q.fna"),
        };
        let (out, err, rc) = search(
            program,
            db,
            query,
            &["-seqidlist", "tests/db_search/seqidlist/missing.txt"],
        );
        let err = String::from_utf8_lossy(&err);
        assert_eq!(rc, 255, "{program}: {err}");
        assert!(out.is_empty(), "{program}");
        assert!(err.starts_with("Error: NCBI C++ Exception:\n    T0 \""), "{program}: {err}");
        assert!(
            err.ends_with(
                "/c++/src/corelib/ncbifile.cpp\", line 5754: Error: (CFileException::eMemoryMap) \
                 ncbi::CMemoryFileMap::CMemoryFileMap() - To be memory mapped the file must exist: \
                 ''\n\n"
            ),
            "{program}: {err}"
        );
        let (_, err, rc) = search(
            program,
            db,
            query,
            &["-negative_seqidlist", "tests/db_search/seqidlist/missing.txt"],
        );
        assert_eq!(rc, 255, "{program} negative: {}", String::from_utf8_lossy(&err));
    }
}

/// A list naming the PDB molecule `3ABC` selects the OID of `lcl|Subject_29|pdb|3ABC|C` of
/// `ids/prot_v5` through the lookup table but matches none of its definition lines (NCBI's
/// `FindId` asks for `3ABC_C`): NCBI counts the sequence and leaves its hits out of every
/// report; LOSAT refuses the search before it (DECISIONS 17).
#[test]
fn a_list_that_keeps_a_sequence_without_a_definition_line_is_refused() {
    let dir = std::env::temp_dir().join(format!("losat-di8b-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let list = dir.join("mol.txt");
    fs::write(&list, "3ABC\nDAA10001.1\n").unwrap();
    let (out, err, rc) = search(
        "blastp",
        "tests/db_search/ids/prot_v5",
        "tests/db_search/ids/ids_q.faa",
        &["-seqidlist", list.to_str().unwrap()],
    );
    fs::remove_dir_all(&dir).ok();
    let err = String::from_utf8_lossy(&err);
    assert_eq!(rc, 1, "{err}");
    assert!(out.is_empty());
    assert!(
        err.contains("but none of its definition lines, and NCBI leaves its hits out of the report without notice")
            && err.contains("this database is not supported by LOSAT's BLASTP"),
        "{err}"
    );
}
