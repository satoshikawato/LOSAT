//! `losat makeblastdb` on FASTA records of 4,294,967,296 letters or more (gaps included), whose
//! length NCBI's `CFastaReader` counts in a 32-bit `TSeqPos` (line RG DECISIONS K3-6, K3-8),
//! against the bytes of makeblastdb 2.17.0 (the official Linux x64 release, whose build tree the
//! exception reports name; probe of 2026-11-13, remaining-gaps-20261112/K/k3/f2probe/cli/gen.py).
//! None of these inputs makes NCBI allocate the letters of its gaps: `AddSequences` skips the
//! record (a gi above 4294967295, a length of 0 modulo 2^32, a delta that starts with a gap of
//! unknown length), stops at its wrapped length (above 2,147,483,647), or `s_FixBioseqDeltas`
//! throws for a protein whose gaps stayed a delta. A record whose data NCBI builds from more
//! letters than its length is LOSAT's refusal (`refused_where_ncbi_builds_the_letters`).
//!
//! Each case runs the `losat` binary in a directory of its own. The standard error stream and the
//! exit status are compared as they are; the standard output with the time, the run directory and
//! the seconds of the log masked.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:290-304
//! (`s_FixBioseqDeltas`), :805-876 (`AddSequences`); c++/src/objtools/readers/fasta.cpp:1074,
//! 1347, 1386 (`m_TotalGapLength`).

use std::path::PathBuf;
use std::process::Command;

/// The build tree that makeblastdb 2.17.0's reports name (the official Linux x64 release).
const ROOT: &str = "/home/coremake/release_build/build/PrepareRelease_Linux64-Alma_JSID_01_170042_130.14.18.128_9008__PrepareRelease_Linux64-Alma_1751374548/c++/compilers/unix/../..";

struct Case {
    name: &'static str,
    input: &'static [u8],
    dbtype: &'static str,
    parse_seqids: bool,
    status: i32,
    stdout: &'static str,
    /// The standard error stream; `{ROOT}` is the build tree.
    stderr: &'static str,
}

#[rustfmt::skip]
const CASES: [Case; 9] = [
    Case { name: "big_gi_nucl", input: ">gi|4294967296 big gi\nA\n>?4294967295\n>?1\nA\n>b\nACGTACGT\n".as_bytes(), dbtype: "nucl", parse_seqids: true, status: 0, stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Nucleotide\nKeep MBits: T\nMaximum file size: 3000000000B\nIgnoring gi '4294967296' as it has value larger then 4294967295\nAdding sequences from FASTA; added 1 sequences in <S> seconds.\n\n\n", stderr: "" },
    Case { name: "big_gi_prot", input: ">gi|4294967296 big gi\nM\n>?4294967295\n>?1\nK\n>b\nMKVLAAGG\n".as_bytes(), dbtype: "prot", parse_seqids: true, status: 0, stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nIgnoring gi '4294967296' as it has value larger then 4294967295\nAdding sequences from FASTA; added 1 sequences in <S> seconds.\n\n\n", stderr: "" },
    Case { name: "wrapped_to_0_nucl", input: ">a\nA\n>?4294967295\n>b\nACGTACGT\n".as_bytes(), dbtype: "nucl", parse_seqids: false, status: 0, stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Nucleotide\nKeep MBits: T\nMaximum file size: 3000000000B\nIgnoring sequence 'lcl|1' as it has no sequence data\nAdding sequences from FASTA; added 1 sequences in <S> seconds.\n\n\n", stderr: "" },
    Case { name: "wrapped_to_0_gaps", input: ">a\n>?4294967295\n>?1\n>b\nACGTACGT\n".as_bytes(), dbtype: "nucl", parse_seqids: false, status: 0, stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Nucleotide\nKeep MBits: T\nMaximum file size: 3000000000B\nIgnoring sequence 'lcl|1' as it has no sequence data\nAdding sequences from FASTA; added 1 sequences in <S> seconds.\n\n\n", stderr: "" },
    Case { name: "wrapped_to_0_prot", input: ">p\n>?4294967295\n>?1\n>q\nMKVLAAGG\n".as_bytes(), dbtype: "prot", parse_seqids: false, status: 0, stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nIgnoring sequence 'lcl|1' as it has no sequence data\nAdding sequences from FASTA; added 1 sequences in <S> seconds.\n\n\n", stderr: "" },
    Case { name: "lone_hyphen_first", input: ">a\n-\n>?4294967295\n>?1\nA\n>b\nACGTACGT\n".as_bytes(), dbtype: "nucl", parse_seqids: false, status: 0, stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Nucleotide\nKeep MBits: T\nMaximum file size: 3000000000B\nIgnoring sequence 'lcl|1' as it has no sequence data\nAdding sequences from FASTA; added 1 sequences in <S> seconds.\n\n\n", stderr: "" },
    Case { name: "too_long_after_wrap", input: ">a\n>?4294967295\n>?4294967295\n>b\nACGTACGT\n".as_bytes(), dbtype: "nucl", parse_seqids: false, status: 1, stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Nucleotide\nKeep MBits: T\nMaximum file size: 3000000000B\n\nNo volumes were created.\n\n", stderr: "BLAST Database creation error: Sequences longer than 2,147,483,647 bases are not supported.  Offending sequence is lcl|1\n" },
    Case { name: "protein_delta", input: ">p\nMKV\n>?4294967295\n>?1\nKKKK\n".as_bytes(), dbtype: "prot", parse_seqids: false, status: 255, stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\n\nNo volumes were created.\n\n", stderr: "Error: NCBI C++ Exception:\n    T0 \"{ROOT}/src/objtools/blast/seqdb_writer/build_db.cpp\", line 303: Error: (CMultisourceException::eArgErr) BLASTDB::ncbi::s_FixBioseqDeltas() - lcl|1: Protein delta sequences are not supported.\n\n" },
    Case { name: "protein_delta_parsed_id", input: ">q\nMKVL\n>sp|P12345|ABC_HUMAN x\nM\n>?2147483648\n>?2147483648\nK\n".as_bytes(), dbtype: "prot", parse_seqids: true, status: 255, stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\n\n\n", stderr: "Error: NCBI C++ Exception:\n    T0 \"{ROOT}/src/objtools/blast/seqdb_writer/build_db.cpp\", line 303: Error: (CMultisourceException::eArgErr) BLASTDB::ncbi::s_FixBioseqDeltas() - sp|P12345|ABC_HUMAN: Protein delta sequences are not supported.\n\n" },
];

/// A directory of its own for `name`, removed when the value is dropped.
struct TempDir(PathBuf);

impl TempDir {
    fn new(name: &str) -> Self {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos());
        let dir = std::env::temp_dir().join(format!(
            "losat-mkdb-32-{name}-{}-{nanos}",
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

/// `losat makeblastdb -in in.fa -dbtype <dbtype> -out db` (and `-parse_seqids`) on `input` in a
/// directory of its own, without LOSAT's variables.
fn makeblastdb(name: &str, input: &[u8], dbtype: &str, parse_seqids: bool) -> std::process::Output {
    let dir = TempDir::new(name);
    std::fs::write(dir.0.join("in.fa"), input).unwrap();
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    command.current_dir(&dir.0).args([
        "makeblastdb",
        "-in",
        "in.fa",
        "-dbtype",
        dbtype,
        "-out",
        "db",
    ]);
    if parse_seqids {
        command.arg("-parse_seqids");
    }
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") {
            command.env_remove(name);
        }
    }
    command.output().expect("run losat makeblastdb")
}

#[test]
fn records_past_32_bits_end_as_makeblastdb_ends() {
    for case in &CASES {
        let output = makeblastdb(case.name, case.input, case.dbtype, case.parse_seqids);
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

/// A nucleotide record whose delta starts with a literal, and a protein whose `>?` lines do not
/// sum to a multiple of 2^32: NCBI builds their data from the letters of the gaps (4 GiB or
/// more; "BLAST ran out of memory" under a cap), LOSAT refuses the record and allocates nothing.
#[test]
fn refused_where_ncbi_builds_the_letters() {
    for (name, input, dbtype, id, length) in [
        (
            "nucl",
            &b">a\nACGT\n>?4294967295\n>?4294967295\nAC\n"[..],
            "nucl",
            "lcl|1",
            "8,589,934,596",
        ),
        (
            "prot",
            b">p\nMKV\n>?4294967295\n>?2\nKKKK\n",
            "prot",
            "lcl|1",
            "4,294,967,304",
        ),
    ] {
        let output = makeblastdb(name, input, dbtype, false);
        assert_eq!(output.status.code(), Some(255), "{name}");
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            format!(
                "Error: a FASTA record of 4,294,967,296 letters or more, counting the letters of its \
                 gaps ({id}: {length} letters; NCBI counts the length of a record in 32 bits, which \
                 wrap there), is not supported by LOSAT's makeblastdb\n"
            ),
            "{name}"
        );
    }
}
