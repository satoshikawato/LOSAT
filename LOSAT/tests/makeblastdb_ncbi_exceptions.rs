//! `losat makeblastdb` with the environment variables of NCBI's LMDB writer set to values that
//! `NStr` cannot convert, against the bytes of makeblastdb 2.17.0 (the official Linux x64 release,
//! whose build tree the exception reports name; probe of 2026-11-08). The variables are read
//! when the first sequence is published: with one sequence that is in `Close`, whose exception
//! `EndBuild` turns into a creation error with the message only (exit status 1); with two it is
//! in `AddSequence`, and the exception reaches the general handler of makeblastdb, which prints
//! `Error: ` and its report (exit status 255). `LMDB_SPLIT_CHUNK_SIZE` and `LMDB_MIN_SPLIT_SIZE`
//! are read in the destructor of NCBI's LMDB writer, where makeblastdb aborts; LOSAT stops with
//! the creation error (PD-LOSAT-DATABASE-SEARCH, the inputs on which makeblastdb crashes).
//!
//! Each case runs the `losat` binary in a directory of its own, with the variable in the child's
//! environment only. The standard error stream and the exit status are compared as they are; the
//! standard output with the time, the run directory and the seconds of the log masked.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1136-1141
//! (`BLASTDB_LMDB_MAP_SIZE`, `NStr::StringToUInt8`), c++/src/objtools/blast/seqdb_writer/
//! writedb_lmdb.cpp:58-61 (`MAX_LMDB_TXN_ENTRY`, `NStr::StringToInt`), :281-290 (the split
//! variables, `NStr::StringToUInt`), c++/src/objtools/blast/seqdb_writer/build_db.cpp:1423-1473
//! (`EndBuild`).

use std::path::PathBuf;
use std::process::Command;

/// The build tree that makeblastdb 2.17.0's reports name (the official Linux x64 release).
const ROOT: &str = "/home/coremake/release_build/build/PrepareRelease_Linux64-Alma_JSID_01_170042_130.14.18.128_9008__PrepareRelease_Linux64-Alma_1751374548/c++/compilers/unix/../..";

const ONE: &[u8] = b">seq1 first\nMKVLAAGIVGLLLAAPAQA\n";
const TWO: &[u8] = b">seq1 first\nMKVLAAGIVGLLLAAPAQA\n>seq2 second\nMKKLLPTAAAGLLLLAAQPAMA\n";

/// The variables of the LMDB writer, removed from the child's environment before the case's.
const VARIABLES: [&str; 4] = [
    "BLASTDB_LMDB_MAP_SIZE",
    "MAX_LMDB_TXN_ENTRY",
    "LMDB_SPLIT_CHUNK_SIZE",
    "LMDB_MIN_SPLIT_SIZE",
];

/// The log of a build that stopped at its second sequence ("No volumes were created").
const LOG_TWO: &str = "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\n\nNo volumes were created.\n\n";
/// The log of a build of one sequence whose close failed.
const LOG_ONE: &str = "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 1 sequences in <S> seconds.\n\nNo volumes were created.\n\n";

struct Case {
    name: &'static str,
    variable: &'static str,
    value: &'static str,
    input: &'static [u8],
    parse_seqids: bool,
    status: i32,
    stdout: &'static str,
    /// The standard error stream; `{ROOT}` is the build tree.
    stderr: &'static str,
}

/// `Error: ` and the report of a `CStringException` of ncbistr.cpp.
const fn case(
    name: &'static str,
    variable: &'static str,
    value: &'static str,
    stderr: &'static str,
) -> Case {
    Case {
        name,
        variable,
        value,
        input: TWO,
        parse_seqids: false,
        status: 255,
        stdout: LOG_TWO,
        stderr,
    }
}

/// The creation error of one sequence.
const fn creation(
    name: &'static str,
    variable: &'static str,
    value: &'static str,
    stderr: &'static str,
) -> Case {
    Case {
        name,
        variable,
        value,
        input: ONE,
        parse_seqids: false,
        status: 1,
        stdout: LOG_ONE,
        stderr,
    }
}

const MAP: &str = "BLASTDB_LMDB_MAP_SIZE";
const TXN: &str = "MAX_LMDB_TXN_ENTRY";

const CASES: [Case; 15] = [
    creation("mapsz_abc_1seq", MAP, "abc", "BLAST Database creation error: Cannot convert string 'abc' to Uint8\n"),
    case("mapsz_abc_2seq", MAP, "abc", "Error: NCBI C++ Exception:\n    T0 \"{ROOT}/src/corelib/ncbistr.cpp\", line 908: Error: (CStringException::eConvert) ncbi::NStr::StringToUInt8() - Cannot convert string 'abc' to Uint8 (m_Pos = 0)\n\n"),
    creation("mapsz_plus_1seq", MAP, "+", "BLAST Database creation error: Cannot convert string '+' to Uint8\n"),
    case("mapsz_plus_2seq", MAP, "+", "Error: NCBI C++ Exception:\n    T0 \"{ROOT}/src/corelib/ncbistr.cpp\", line 895: Error: (CStringException::eConvert) ncbi::NStr::StringToUInt8() - Cannot convert string '+' to Uint8 (m_Pos = 1)\n\n"),
    creation("mapsz_1x_1seq", MAP, "1x", "BLAST Database creation error: Cannot convert string '1x' to Uint8\n"),
    case("mapsz_1x_2seq", MAP, "1x", "Error: NCBI C++ Exception:\n    T0 \"{ROOT}/src/corelib/ncbistr.cpp\", line 908: Error: (CStringException::eConvert) ncbi::NStr::StringToUInt8() - Cannot convert string '1x' to Uint8 (m_Pos = 1)\n\n"),
    // The fast conversion's overflow names no ", overflow".
    creation("mapsz_ovf_1seq", MAP, "99999999999999999999", "BLAST Database creation error: Cannot convert string '99999999999999999999' to Uint8\n"),
    case("mapsz_ovf_2seq", MAP, "99999999999999999999", "Error: NCBI C++ Exception:\n    T0 \"{ROOT}/src/corelib/ncbistr.cpp\", line 912: Error: (CStringException::eConvert) ncbi::NStr::StringToUInt8() - Cannot convert string '99999999999999999999' to Uint8 (m_Pos = 19)\n\n"),
    creation("txn_abc_1seq", TXN, "abc", "BLAST Database creation error: Cannot convert string 'abc' to Int8\n"),
    case("txn_abc_2seq", TXN, "abc", "Error: NCBI C++ Exception:\n    T0 \"{ROOT}/src/corelib/ncbistr.cpp\", line 861: Error: (CStringException::eConvert) ncbi::NStr::StringToInt8() - Cannot convert string 'abc' to Int8 (m_Pos = 0)\n\n"),
    case("txn_endptr_2seq", TXN, "-1x", "Error: NCBI C++ Exception:\n    T0 \"{ROOT}/src/corelib/ncbistr.cpp\", line 870: Error: (CStringException::eConvert) ncbi::NStr::StringToInt8() - Cannot convert string '-1x' to Int8 (m_Pos = 2)\n\n"),
    case("txn_intovf_2seq", TXN, "99999999999", "Error: NCBI C++ Exception:\n    T0 \"{ROOT}/src/corelib/ncbistr.cpp\", line 640: Error: (CStringException::eConvert) ncbi::NStr::StringToInt() - Cannot convert string '99999999999' to int, overflow (m_Pos = 0)\n\n"),
    case("txn_i8ovf_2seq", TXN, "99999999999999999999", "Error: NCBI C++ Exception:\n    T0 \"{ROOT}/src/corelib/ncbistr.cpp\", line 852: Error: (CStringException::eConvert) ncbi::NStr::StringToInt8() - Cannot convert string '99999999999999999999' to Int8, overflow (m_Pos = 18)\n\n"),
    // makeblastdb aborts in `~CWriteDB_LMDB` (SIGABRT) after the line "Adding sequences"; LOSAT
    // stops there with a creation error, its volumes removed.
    Case {
        parse_seqids: true,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 1 sequences in <S> seconds.\n\n\n",
        ..creation("split_chunk_abc_1seq", "LMDB_SPLIT_CHUNK_SIZE", "abc", "BLAST Database creation error: Cannot convert string 'abc' to Uint8\n")
    },
    Case {
        input: TWO,
        parse_seqids: true,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 2 sequences in <S> seconds.\n\n\n",
        ..creation("split_min_abc_2seq", "LMDB_MIN_SPLIT_SIZE", "abc", "BLAST Database creation error: Cannot convert string 'abc' to Uint8\n")
    },
];

/// A directory of its own for `name`, removed when the value is dropped.
struct TempDir(PathBuf);

impl TempDir {
    fn new(name: &str) -> Self {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos());
        let dir = std::env::temp_dir().join(format!(
            "losat-mkdb-exc-{name}-{}-{nanos}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        TempDir(dir)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// `text` with the line parts that change from run to run replaced: the time of the log head,
/// the directory in "New DB name" and the seconds of "Adding sequences".
fn masked(text: &str) -> String {
    text.lines()
        .map(|line| {
            if line.starts_with("Building a new DB, current time: ") {
                "Building a new DB, current time: <T>".to_string()
            } else if let (true, Some(name)) = (
                line.starts_with("New DB name:   "),
                line.rsplit_once(['/', '\\']),
            ) {
                format!("New DB name:   <D>/{}", name.1)
            } else if let (true, Some(at), Some(end)) = (
                line.starts_with("Adding sequences from FASTA; added "),
                line.find(" sequences in "),
                line.rfind(" seconds."),
            ) {
                format!("{} sequences in <S>{}", &line[..at], &line[end..])
            } else {
                line.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
        + if text.ends_with('\n') { "\n" } else { "" }
}

#[test]
fn unconvertible_lmdb_variables_end_as_makeblastdb_ends() {
    for case in &CASES {
        let dir = TempDir::new(case.name);
        std::fs::write(dir.0.join("in.fa"), case.input).unwrap();
        let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
        command.current_dir(&dir.0).args([
            "makeblastdb",
            "-in",
            "in.fa",
            "-dbtype",
            "prot",
            "-out",
            "db",
        ]);
        if case.parse_seqids {
            command.arg("-parse_seqids");
        }
        for (name, _) in std::env::vars_os() {
            if name.to_string_lossy().starts_with("LOSAT_") {
                command.env_remove(name);
            }
        }
        for name in VARIABLES {
            command.env_remove(name);
        }
        let output = command
            .env(case.variable, case.value)
            .output()
            .expect("run losat makeblastdb");
        assert_eq!(output.status.code(), Some(case.status), "{}", case.name);
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            case.stderr.replace("{ROOT}", ROOT),
            "{}",
            case.name
        );
        assert_eq!(
            masked(&String::from_utf8_lossy(&output.stdout)),
            case.stdout,
            "{}",
            case.name
        );
    }
}

#[test]
fn the_mask_keeps_everything_else() {
    assert_eq!(
        masked("\n\nBuilding a new DB, current time: 11/08/2026 19:07:40\nNew DB name:   /a b/c/db\nAdding sequences from FASTA; added 1 sequences in 0.000256197 seconds.\n\n"),
        "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nAdding sequences from FASTA; added 1 sequences in <S> seconds.\n\n"
    );
}
