//! `losat makeblastdb` where an output file cannot be created or written, and where the last
//! sequence has a duplicate id, against the files and messages of makeblastdb 2.17.0 (the
//! official Linux x64 release; probes of 2026-11-12, `sw8/mkdb3` of the SWEBOK line).
//!
//! makeblastdb writes every output file through an `ofstream` it never tests: a file that
//! cannot be opened (a directory of its name that the removal of the old database cannot
//! delete, because it holds a directory that cannot be read) or written (a link to `/dev/full`)
//! is written to nowhere and the build goes on. The volume file under its numbered name is
//! renamed whatever it is (a directory, a link). The database is then opened: a missing index
//! file or alias file ends with the reader's message and exit status 2, any other missing file
//! is left out of the metadata. The metadata file is the one stream whose failure is reported,
//! by NCBI's serializer when it is destroyed; a flush that fails ends with exit status 6.
//! A duplicate id in the last sequence is found when the close publishes it: the files of the
//! build are removed, then the writer's destructors write the lookup files and the files a
//! volume makes when it closes.
//!
//! Each case runs the `losat` binary in a directory of its own, with the setup of the case made
//! before it. A case lists the top-level entries of the output directory with their sizes
//! (`name:size`, a directory `name/`, a link `name@`); the files that hold the build's date
//! (`.pin`, `.nin`, `.pjs`, `.njs`, `.pal`, `.nal`) have a size that depends on the time of day
//! (`*`). The cases of many volumes list the difference from the same build without the setup:
//! `+entry` added, `-name` gone, `~entry` changed.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:196-221
//! (`Create`, `Write`), c++/src/app/blastdb/makeblastdb.cpp:1296-1318 (the metadata file),
//! c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:100-108 (`~CWriteDB_Impl`).

#![cfg(unix)]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Proteins with a GI, a trace id and an accession (the probe's `RICH_P`).
const RICH_P: &str = ">gi|1001|ref|NP_000001.1| alpha\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV\n>gnl|ti|12345 trace one\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNA\n>sp|P55555|CCC_HUMAN Gamma\nMSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARR\n";
/// Three proteins whose last one repeats the id of the second.
const DUP_END: &str = ">a\nMKTAYIAKQ\n>b\nMKTAYIAKR\n>b\nMKTAYIAKS\n";

const LOG_3: &str = "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n";
const LOG_130: &str = "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 0B\nAdding sequences from FASTA; added 130 sequences in <S> seconds.\n\n\n";
const NO_INDEX: &str = "BLAST Database error: No alias or index file found for protein database [<D>/db] in search path [<C>::]. Please verify the spelling of the BLAST database and its molecule type.\n";

/// 130 proteins of 60 residues.
fn many() -> String {
    const AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let mut out = String::new();
    for i in 0..130usize {
        out.push_str(&format!(">sp|P{i:05}|S{i} desc {i}\n"));
        out.extend((0..60usize).map(|j| AA[(i * 7 + j * 3 + j * j % 11) % 20] as char));
        out.push('\n');
    }
    out
}

/// The kinds of build: the arguments and the input.
#[derive(Clone, Copy)]
enum Kind {
    /// Protein, `-parse_seqids -hash_index`, version 5, `RICH_P`.
    P5,
    /// `P5` with `-blastdb_version 4`.
    P4,
    /// `P5` with `DUP_END`.
    P5Dup,
    /// `P4` with `DUP_END`.
    P4Dup,
    /// 130 proteins, `-parse_seqids -max_file_sz 2000B`: 6 volumes.
    Six,
    /// 130 proteins, `-parse_seqids -max_file_sz 100B`: 130 volumes.
    Many,
}

#[derive(Clone, Copy)]
enum Setup {
    /// A directory that holds a directory that cannot be read (its removal fails).
    Blocked(&'static str),
    /// A symbolic link to `/dev/full`.
    Full(&'static str),
    None,
}

struct Case {
    name: &'static str,
    kind: Kind,
    setup: Setup,
    status: i32,
    stdout: &'static str,
    stderr: &'static str,
    /// The entries of the output directory, or their difference from a build without the setup.
    files: &'static str,
}

const CASES: &[Case] = &[
    // The index file under its numbered name cannot be created: the directory is renamed to
    // `db.pin` and the database cannot be opened.
    Case {
        name: "b_p5_00pin",
        kind: Kind::P5,
        setup: Setup::Blocked("db.00.pin"),
        status: 2,
        stdout: LOG_3,
        stderr: NO_INDEX,
        files: "db.pdb:32768 db.phd:37 db.phi:65 db.phr:221 db.pin/ db.pnd:8 db.pni:52 db.pog:44 db.pos:60 db.pot:44 db.psq:214 db.ptf:16384 db.pto:16",
    },
    // The header file cannot be created: the database opens and its metadata leaves it out.
    Case {
        name: "b_p5_00phr",
        kind: Kind::P5,
        setup: Setup::Blocked("db.00.phr"),
        status: 0,
        stdout: LOG_3,
        stderr: "",
        files: "db.pdb:32768 db.phd:37 db.phi:65 db.phr/ db.pin:* db.pjs:* db.pnd:8 db.pni:52 db.pog:44 db.pos:60 db.pot:44 db.psq:214 db.ptf:16384 db.pto:16",
    },
    Case {
        name: "b_p4_00psi",
        kind: Kind::P4,
        setup: Setup::Blocked("db.00.psi"),
        status: 0,
        stdout: LOG_3,
        stderr: "",
        files: "db.phd:37 db.phi:65 db.phr:221 db.pin:* db.pjs:* db.pnd:8 db.pni:52 db.pog:44 db.psd:62 db.psi/ db.psq:214",
    },
    Case {
        name: "b_p5_pos",
        kind: Kind::P5,
        setup: Setup::Blocked("db.pos"),
        status: 0,
        stdout: LOG_3,
        stderr: "",
        files: "db.pdb:32768 db.phd:37 db.phi:65 db.phr:221 db.pin:* db.pjs:* db.pnd:8 db.pni:52 db.pog:44 db.pos/ db.pot:44 db.psq:214 db.ptf:16384 db.pto:16",
    },
    // The metadata file cannot be created: the serializer reports the write fault at the end
    // of the text (479 bytes) and makeblastdb ends with exit status 0.
    Case {
        name: "b_p5_pjs",
        kind: Kind::P5,
        setup: Setup::Blocked("db.pjs"),
        status: 0,
        stdout: LOG_3,
        stderr: "Error: (802.5) [makeblastdb] CObjectOStream: error at byte 479: stack is empty: write fault\n",
        files: "db.pdb:32768 db.phd:37 db.phi:65 db.phr:221 db.pin:* db.pjs/ db.pnd:8 db.pni:52 db.pog:44 db.pos:60 db.pot:44 db.psq:214 db.ptf:16384 db.pto:16",
    },
    // The duplicate id of the last sequence: the files of the build are removed, then the
    // destructors write the `.pot`, `.pto` and `.pos` files, commit `.ptf`, and close the volume,
    // whose ISAM and OID to GI files are made under its numbered name.
    Case {
        name: "d_p5_end",
        kind: Kind::P5Dup,
        setup: Setup::None,
        status: 1,
        stdout: LOG_3,
        stderr: "BLAST Database creation error: Error: Duplicate seq_ids are found: \nLCL|B\n\n",
        files: "db.00.phd:26 db.00.phi:65 db.00.pog:40 db.pos:28 db.pot:32 db.ptf:16384 db.pto:12",
    },
    Case {
        name: "d_p4_end",
        kind: Kind::P4Dup,
        setup: Setup::None,
        status: 1,
        stdout: LOG_3,
        stderr: "BLAST Database creation error: Error: Duplicate seq_ids are found: \nLCL|B\n\n",
        files: "db.00.phd:26 db.00.phi:65 db.00.pog:40 db.00.psd:24 db.00.psi:56",
    },
    // The alias file of several volumes cannot be created.
    Case {
        name: "b_m2_pal",
        kind: Kind::Six,
        setup: Setup::Blocked("db.pal"),
        status: 2,
        stdout: LOG_130,
        stderr: NO_INDEX,
        files: "~db.pal/ -db.pjs",
    },
    Case {
        name: "b_m2_00pin",
        kind: Kind::Six,
        setup: Setup::Blocked("db.00.pin"),
        status: 2,
        stdout: LOG_130,
        stderr: "BLAST Database error: Could not find volume or alias file (db.00) referenced in alias file (<D>/db).\n",
        files: "~db.00.pin/ -db.pjs",
    },
    // A metadata text of more than 4096 bytes fails at the serializer's first write.
    Case {
        name: "b_m_pjs",
        kind: Kind::Many,
        setup: Setup::Blocked("db.pjs"),
        status: 0,
        stdout: LOG_130,
        stderr: "Error: (802.5) [makeblastdb] CObjectOStream: error at byte 4096: stack is empty: write fault\n",
        files: "~db.pjs/",
    },
];

/// The cases on a device that takes no byte (Linux only).
const FULL_CASES: &[Case] = &[
    // The flush of `COStreamBuffer::Flush` fails: `eFlush` ends makeblastdb with exit status 6
    // after the serializer's report.
    Case {
        name: "f_p5_pjs",
        kind: Kind::P5,
        setup: Setup::Full("db.pjs"),
        status: 6,
        stdout: LOG_3,
        stderr: "Error: (802.5) [makeblastdb] CObjectOStream: error at byte 479: stack is empty: cannot close output stream\nError: (CIOException::eFlush) [makeblastdb] COStreamBuffer::Flush: failed\nError: (CSerialException::eIoError) [makeblastdb] byte 479: cannot close output stream\nError: (802.4) [makeblastdb] Cannot close serializing output stream (CSerialException::eIoError) byte 479: cannot close output stream\nBLAST failed to write output: COStreamBuffer::Flush: failed\n",
        files: "db.pdb:32768 db.phd:37 db.phi:65 db.phr:221 db.pin:* db.pjs@ db.pnd:8 db.pni:52 db.pog:44 db.pos:60 db.pot:44 db.psq:214 db.ptf:16384 db.pto:16",
    },
    // The link under the numbered name is renamed as the link: `db.pin` is not a file.
    Case {
        name: "f_p5_00pin",
        kind: Kind::P5,
        setup: Setup::Full("db.00.pin"),
        status: 2,
        stdout: LOG_3,
        stderr: NO_INDEX,
        files: "db.pdb:32768 db.phd:37 db.phi:65 db.phr:221 db.pin@ db.pnd:8 db.pni:52 db.pog:44 db.pos:60 db.pot:44 db.psq:214 db.ptf:16384 db.pto:16",
    },
    // The first 4096 bytes fit the stream's buffer; the next 4096 do not, and their write fails.
    Case {
        name: "f_m_pjs",
        kind: Kind::Many,
        setup: Setup::Full("db.pjs"),
        status: 0,
        stdout: LOG_130,
        stderr: "Error: (802.5) [makeblastdb] CObjectOStream: error at byte 8192: stack is empty: write fault\n",
        files: "~db.pjs@",
    },
];

#[derive(Debug, PartialEq, Eq)]
enum Entry {
    File(u64),
    Dir,
    Link,
}

struct TempDir(PathBuf);

impl TempDir {
    fn new(name: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("losat-mkdb-streams-{name}-{}", std::process::id()));
        if dir.exists() {
            writable(&dir);
            std::fs::remove_dir_all(&dir).unwrap();
        }
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        writable(&self.0);
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// Gives every directory under `dir` its permissions back (the cases make some unreadable).
fn writable(dir: &Path) {
    use std::os::unix::fs::PermissionsExt;
    let _ = std::fs::set_permissions(dir, std::fs::Permissions::from_mode(0o755));
    if let Ok(entries) = std::fs::read_dir(dir) {
        for entry in entries.map_while(Result::ok) {
            if entry.file_type().is_ok_and(|t| t.is_dir()) {
                writable(&entry.path());
            }
        }
    }
}

/// A directory at `path` that holds a directory that cannot be read: removing it fails.
fn blocked(path: &Path) {
    use std::os::unix::fs::PermissionsExt;
    std::fs::create_dir_all(path.join("u").join("x")).unwrap();
    std::fs::write(path.join("u").join("x").join("f"), b"x").unwrap();
    std::fs::set_permissions(path.join("u"), std::fs::Permissions::from_mode(0o000)).unwrap();
}

struct Outcome {
    status: i32,
    stdout: String,
    stderr: String,
    entries: BTreeMap<String, Entry>,
}

fn run(case_name: &str, kind: Kind, setup: Setup) -> Outcome {
    let tmp = TempDir::new(case_name);
    let d = tmp.0.join("d");
    std::fs::create_dir(&d).unwrap();
    let mut args = vec!["-dbtype", "prot", "-parse_seqids"];
    let input = match kind {
        Kind::P5 | Kind::P4 => RICH_P.to_string(),
        Kind::P5Dup | Kind::P4Dup => DUP_END.to_string(),
        Kind::Six | Kind::Many => many(),
    };
    match kind {
        Kind::P5 | Kind::P5Dup => args.push("-hash_index"),
        Kind::P4 | Kind::P4Dup => args.extend(["-hash_index", "-blastdb_version", "4"]),
        Kind::Six => args.extend(["-max_file_sz", "2000B"]),
        Kind::Many => args.extend(["-max_file_sz", "100B"]),
    }
    std::fs::write(tmp.0.join("in.fa"), input).unwrap();
    match setup {
        Setup::Blocked(name) => blocked(&d.join(name)),
        Setup::Full(name) => std::os::unix::fs::symlink("/dev/full", d.join(name)).unwrap(),
        Setup::None => {}
    }
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    command
        .current_dir(&tmp.0)
        // `$HOME/.ncbirc` is read like the one of the current directory: it would change the
        // search path that the messages name.
        .env("HOME", &tmp.0)
        .arg("makeblastdb")
        .args(["-in", "in.fa"])
        .args(&args)
        .arg("-out")
        .arg(d.join("db"));
    for (name, _) in std::env::vars_os() {
        let name = name.to_string_lossy().into_owned();
        if name.starts_with("LOSAT_")
            || name.starts_with("NCBI_CONFIG")
            || name == "BLASTDB"
            || name == "BLASTDB_LMDB_MAP_SIZE"
            || name == "MAX_LMDB_TXN_ENTRY"
            || name.starts_with("LMDB_")
        {
            command.env_remove(name);
        }
    }
    let output = command.output().expect("run losat makeblastdb");
    let mask = |bytes: &[u8]| {
        masked(
            &String::from_utf8_lossy(bytes),
            d.to_str().unwrap(),
            tmp.0.to_str().unwrap(),
        )
    };
    let mut entries = BTreeMap::new();
    for entry in std::fs::read_dir(&d).unwrap().map_while(Result::ok) {
        let meta = std::fs::symlink_metadata(entry.path()).unwrap();
        let kind = if meta.file_type().is_symlink() {
            Entry::Link
        } else if meta.is_dir() {
            Entry::Dir
        } else {
            Entry::File(meta.len())
        };
        entries.insert(entry.file_name().to_string_lossy().into_owned(), kind);
    }
    Outcome {
        status: output.status.code().unwrap_or(-1),
        stdout: mask(&output.stdout),
        stderr: mask(&output.stderr),
        entries,
    }
}

/// The time, the seconds, the output directory and the run directory of a message masked.
fn masked(text: &str, d: &str, c: &str) -> String {
    let mut out = text.replace(d, "<D>").replace(c, "<C>");
    for (start, end, fill) in [
        ("current time: ", "\n", "<T>"),
        (" sequences in ", " seconds.", "<S>"),
    ] {
        let mut at = 0;
        while let Some(i) = out[at..].find(start) {
            let from = at + i + start.len();
            let Some(j) = out[from..].find(end) else {
                break;
            };
            out.replace_range(from..from + j, fill);
            at = from + fill.len();
        }
    }
    out
}

/// The name of a file that holds the build's date: its size is not compared.
fn dated(name: &str) -> bool {
    [".pin", ".nin", ".pjs", ".njs", ".pal", ".nal"]
        .iter()
        .any(|ext| name.ends_with(ext))
}

fn entry_text(name: &str, entry: &Entry) -> String {
    match entry {
        Entry::Dir => format!("{name}/"),
        Entry::Link => format!("{name}@"),
        Entry::File(_) if dated(name) => format!("{name}:*"),
        Entry::File(size) => format!("{name}:{size}"),
    }
}

fn listing(entries: &BTreeMap<String, Entry>) -> String {
    entries
        .iter()
        .map(|(name, entry)| entry_text(name, entry))
        .collect::<Vec<_>>()
        .join(" ")
}

fn delta(entries: &BTreeMap<String, Entry>, clean: &BTreeMap<String, Entry>) -> String {
    let mut out = Vec::new();
    for name in entries
        .keys()
        .chain(clean.keys())
        .collect::<std::collections::BTreeSet<_>>()
    {
        match (entries.get(name), clean.get(name)) {
            (None, _) => out.push(format!("-{name}")),
            (Some(entry), None) => out.push(format!("+{}", entry_text(name, entry))),
            (Some(a), Some(b)) if a != b => match (a, b) {
                (Entry::File(_), Entry::File(_)) if dated(name) => {}
                (Entry::File(a), Entry::File(b)) => out.push(format!("~{name}:{b}>{a}")),
                _ => out.push(format!("~{}", entry_text(name, a))),
            },
            _ => {}
        }
    }
    out.join(" ")
}

/// Root reads any directory: the cases that need an unreadable one say nothing then.
fn running_as_root() -> bool {
    let tmp = TempDir::new("root-check");
    blocked(&tmp.0.join("b"));
    std::fs::read_dir(tmp.0.join("b").join("u")).is_ok()
}

fn check(cases: &[Case]) {
    let mut failures = Vec::new();
    for case in cases {
        let outcome = run(case.name, case.kind, case.setup);
        let files = if matches!(case.kind, Kind::Six | Kind::Many) {
            let clean = run(&format!("{}-clean", case.name), case.kind, Setup::None);
            delta(&outcome.entries, &clean.entries)
        } else {
            listing(&outcome.entries)
        };
        let got = (outcome.status, outcome.stdout, outcome.stderr, files);
        let want = (
            case.status,
            case.stdout.to_string(),
            case.stderr.to_string(),
            case.files.to_string(),
        );
        if got != want {
            failures.push(format!("{}:\n  got  {got:?}\n  want {want:?}", case.name));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn makeblastdb_output_files_that_cannot_be_created_end_as_makeblastdb_ends() {
    if running_as_root() {
        eprintln!("skipped (running as root)");
        return;
    }
    check(CASES);
}

#[cfg(target_os = "linux")]
#[test]
fn makeblastdb_output_files_on_a_full_device_end_as_makeblastdb_ends() {
    check(FULL_CASES);
}
