//! `-seqidlist` with `-db` outside the manifest rows (line DI-8b, DI8 DECISIONS 17 and 20):
//! a list file that is not on the search path, which NCBI BLAST+ 2.17.0 reports with the path
//! of its build tree (exit status 255), and a list that keeps a sequence without a definition
//! line (a PDB molecule named without its chain), whose hits NCBI leaves out of its reports
//! without notice and which LOSAT refuses. Also the `-out` file of a list error, which the rows
//! cannot express (a taxid list too). NCBI is not run here; its texts are those of the probe of
//! dropin-survey-20261029/DI8 (S22, R01-R07, E*).

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
        .env(
            "NCBI_CONFIG_PATH",
            engine_dir().join("tests/db_search/seqidlist"),
        )
        .args([program, "-db", db, "-query", query])
        .args(words)
        .args(["-outfmt", "6", "-num_threads", "1"])
        .output()
        .unwrap();
    (
        output.stdout,
        output.stderr,
        output.status.code().unwrap_or(-1),
    )
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
            "blastp" => (
                "tests/db_search/taxfilter/tfp",
                "tests/db_search/taxfilter/q.faa",
            ),
            _ => (
                "tests/db_search/taxfilter/tfn",
                "tests/db_search/taxfilter/q.fna",
            ),
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
        assert!(
            err.starts_with("Error: NCBI C++ Exception:\n    T0 \""),
            "{program}: {err}"
        );
        assert!(
            err.ends_with(
                "/src/corelib/ncbifile.cpp\", line 5754: Error: (CFileException::eMemoryMap) \
                 ncbi::CMemoryFileMap::CMemoryFileMap() - To be memory mapped the file must exist: \
                 ''\n\n"
            ),
            "{program}: {err}"
        );
        let (_, err, rc) = search(
            program,
            db,
            query,
            &[
                "-negative_seqidlist",
                "tests/db_search/seqidlist/missing.txt",
            ],
        );
        assert_eq!(
            rc,
            255,
            "{program} negative: {}",
            String::from_utf8_lossy(&err)
        );
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

/// The programs with their taxonomy fixtures.
const LIST_ERROR_PROGRAMS: [(&str, &str, &str); 5] = [
    (
        "blastn",
        "tests/db_search/taxfilter/tfn",
        "tests/db_search/taxfilter/q.fna",
    ),
    (
        "blastp",
        "tests/db_search/taxfilter/tfp",
        "tests/db_search/taxfilter/q.faa",
    ),
    (
        "blastx",
        "tests/db_search/taxfilter/tfp",
        "tests/db_search/taxfilter/q.fna",
    ),
    (
        "tblastn",
        "tests/db_search/taxfilter/tfn",
        "tests/db_search/taxfilter/q.faa",
    ),
    (
        "tblastx",
        "tests/db_search/taxfilter/tfn",
        "tests/db_search/taxfilter/q.fna",
    ),
];

/// One `losat PROGRAM argv...` in LOSAT/ as `search` runs it, without its fixed format. Returns
/// stderr and the exit status.
fn run_with_out(program: &str, argv: &[&str]) -> (String, i32) {
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
        .env(
            "NCBI_CONFIG_PATH",
            engine_dir().join("tests/db_search/seqidlist"),
        )
        .arg(program)
        .args(argv)
        .output()
        .unwrap();
    (
        String::from_utf8_lossy(&output.stderr).into_owned(),
        output.status.code().unwrap_or(-1),
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3456-3481
// ```c++
//     if (args.Exist(kArgQuery) && args[kArgQuery].HasValue() &&
//         m_InputStream == NULL) {
// ...
//     else {
//         m_OutputStream = &args[kArgOutput].AsOutputFile();
//     }
// ```
// `CStdCmdLineArgs` comes after the handler of the database arguments, which reads the list: a
// list error stops the run before the query and the output are opened (probe E01-E20).
#[test]
fn list_error_leaves_no_out_file() {
    let dir = std::env::temp_dir().join(format!("losat-list-error-out-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let missing = "tests/db_search/taxfilter/no-such-list.txt";
    for (program, db, query) in LIST_ERROR_PROGRAMS {
        let out = dir.join(format!("{program}.out"));
        let out = out.to_str().unwrap();
        for words in [
            ["-taxidlist", missing],
            ["-negative_taxidlist", missing],
            ["-seqidlist", missing],
            ["-negative_seqidlist", missing],
            ["-taxids", "abc"],
        ] {
            let mut argv = vec!["-db", db, "-query", query, "-outfmt", "6", "-out", out];
            argv.extend(words);
            let (stderr, status) = run_with_out(program, &argv);
            assert_ne!(status, 0, "{program} {words:?}: {stderr}");
            assert!(
                !Path::new(out).exists(),
                "{program} {words:?} created -out: {stderr}"
            );
        }
        let taxid_list = ["-taxidlist", missing];
        let accessible = "File is not acessible";

        // The file of an earlier run stays as it was.
        fs::write(out, "OLD\n").unwrap();
        let mut argv = vec!["-db", db, "-query", query, "-outfmt", "6", "-out", out];
        argv.extend(taxid_list);
        let (stderr, status) = run_with_out(program, &argv);
        assert_eq!(status, 1, "{program}: {stderr}");
        assert_eq!(fs::read(out).unwrap(), b"OLD\n", "{program}");
        fs::remove_file(out).unwrap();

        // The list comes before the query and the output: its error is the one reported.
        let mut argv = vec!["-db", db, "-query", "no-such-query.fa", "-out", out];
        argv.extend(taxid_list);
        let (stderr, status) = run_with_out(program, &argv);
        assert_eq!(status, 1, "{program}: {stderr}");
        assert!(stderr.contains(accessible), "{program}: {stderr}");
        assert!(!Path::new(out).exists(), "{program}");
        let mut argv = vec!["-db", db, "-query", query, "-out", "no-such-dir/x.out"];
        argv.extend(taxid_list);
        let (stderr, status) = run_with_out(program, &argv);
        assert_eq!(status, 1, "{program}: {stderr}");
        assert!(stderr.contains(accessible), "{program}: {stderr}");

        // Outfmt 11: `PrintErrorArchive` opens the output after the error and writes the archive.
        let mut argv = vec!["-db", db, "-query", query, "-outfmt", "11", "-out", out];
        argv.extend(taxid_list);
        let (stderr, status) = run_with_out(program, &argv);
        assert_eq!(status, 1, "{program}: {stderr}");
        assert!(stderr.contains(accessible), "{program}: {stderr}");
        let archive = fs::read_to_string(out).unwrap();
        assert!(
            archive.starts_with("Blast4-archive"),
            "{program}: {archive}"
        );
        assert!(archive.contains(accessible), "{program}: {archive}");
        fs::remove_file(out).unwrap();

        // The database opens after the output: its error leaves the file empty (probe E15, E18).
        let argv = [
            "-db",
            "tests/db_search/taxfilter/no-such-db",
            "-query",
            query,
            "-taxids",
            "562",
            "-no_taxid_expansion",
            "-outfmt",
            "6",
            "-out",
            out,
        ];
        let (stderr, status) = run_with_out(program, &argv);
        assert_eq!(status, 2, "{program}: {stderr}");
        assert_eq!(fs::read(out).unwrap(), b"", "{program}");
        fs::remove_file(out).unwrap();
    }
    fs::remove_dir_all(&dir).unwrap();
}
