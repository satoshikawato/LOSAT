//! `losat makeblastdb` where makeblastdb 2.17.0 ignores a volume file that cannot be renamed, and
//! where it creates and leaves LMDB's lock files, against the files and messages of makeblastdb
//! 2.17.0 (the official Linux x64 release; probes of 2026-11-11).
//!
//! makeblastdb gives a volume file its final name with `CDirEntry::Rename(fRF_Overwrite)` and does
//! not look at the result: a name that holds a directory the removal cannot delete (it holds a
//! directory that cannot be read) keeps the volume file under its numbered name, the build goes on
//! and ends where the database is read (`No alias or index file found ...`, `Could not find volume
//! or alias file ...`, exit status 2) or normally (exit status 0) when the missing file is not
//! one the reader opens. A file that is already there under the numbered name of a file the
//! writer never created is renamed as well. LMDB's `mdb_env_open` creates `NAME.pdb-lock` and
//! `NAME.ptf-lock`; the destructors of the writers remove them; one whose environment did not
//! open stays.
//!
//! Each case runs the `losat` binary in a directory of its own, with the setup of the case made
//! before it. A case lists the top-level entries of the output directory with their sizes
//! (`name:size`, a directory `name/`); the files that hold the build's date (`.pin`, `.nin`,
//! `.pjs`, `.njs`, `.pal`, `.nal`) have a size that depends on the time of day (`*`). The cases
//! of more than 100 volumes list the difference from the same build without the setup:
//! `+entry` added, `-name` gone, `~entry` or `~name:from>to` changed.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_files.cpp:267-304
//! (`RenameSingle`, `RenameFileIndex`), c++/src/corelib/ncbifile.cpp:2456-2590
//! (`CDirEntry::Rename`), c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:65-72, 437-444
//! (the destructors), c++/src/util/lmdb/mdb.c:4739-4915 (`mdb_env_setup_locks`).

#![cfg(unix)]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;

const PROT: &str = ">sp|P12345|AAA_HUMAN Alpha protein\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV\n>sp|Q99999|BBB_HUMAN Beta protein\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNA\n>sp|P55555|CCC_HUMAN Gamma\nMSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARR\n";
const NUC: &str = ">seq1 first\nACGTACGTACGTAGCTAGCTAGCTAGCATCGATCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCATGC\n>seq2 second\nTTGACCGGTTAACCGGTTAAGGCCTTAAGGCCTTAAACGTACGTACGTAGCTAGCTAGCATCGATCGATCG\n";

/// The kinds of build: the arguments and the input.
#[derive(Clone, Copy)]
enum Kind {
    /// Protein, `-parse_seqids -hash_index`, version 5.
    P5,
    /// Protein, `-parse_seqids -hash_index -blastdb_version 4`.
    P4,
    /// Nucleotide, `-parse_seqids -hash_index`, version 5.
    N5,
    /// Nucleotide, `-parse_seqids -hash_index -blastdb_version 4`.
    N4,
    /// Protein, no `-parse_seqids`.
    Plain,
    /// Protein, `-parse_seqids`, version 5.
    V5,
    /// 130 proteins, `-parse_seqids -max_file_sz 100B`: 130 volumes.
    Many,
    /// `Many` with `-hash_index`.
    ManyHash,
}

enum Setup {
    /// A directory that holds a directory that cannot be read (its removal fails) at each name.
    Blocked(&'static [&'static str]),
    /// A file `stale NAME` at each name.
    Stale(&'static [&'static str]),
    /// `Stale` for the first name, `Blocked` for the second.
    StaleAndBlocked(&'static str, &'static str),
}

struct Case {
    name: &'static str,
    kind: Kind,
    setup: Setup,
    status: i32,
    /// The standard output with the time, the seconds, the output directory (`<D>`) and the run
    /// directory (`<C>`) masked.
    stdout: &'static str,
    stderr: &'static str,
    /// The entries of the output directory, or their difference from a build without the setup.
    files: &'static str,
}

const CASES: &[Case] = &[
    Case {
        name: "s5_pin",
        kind: Kind::P5,
        setup: Setup::Blocked(&["db.pin"]),
        status: 2,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "BLAST Database error: No alias or index file found for protein database [<D>/db] in search path [<C>::]. Please verify the spelling of the BLAST database and its molecule type.\n",
        files: "db.00.pin:* db.pdb:32768 db.phd:37 db.phi:65 db.phr:264 db.pin/ db.pog:44 db.pos:53 db.pot:44 db.psq:214 db.ptf:16384 db.pto:16",
    },
    Case {
        name: "s5_pin_psq",
        kind: Kind::P5,
        setup: Setup::Blocked(&["db.pin", "db.psq"]),
        status: 2,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "BLAST Database error: No alias or index file found for protein database [<D>/db] in search path [<C>::]. Please verify the spelling of the BLAST database and its molecule type.\n",
        files: "db.00.pin:* db.00.psq:214 db.pdb:32768 db.phd:37 db.phi:65 db.phr:264 db.pin/ db.pog:44 db.pos:53 db.pot:44 db.psq/ db.ptf:16384 db.pto:16",
    },
    Case {
        name: "s5_pog",
        kind: Kind::P5,
        setup: Setup::Blocked(&["db.pog"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.00.pog:44 db.pdb:32768 db.phd:37 db.phi:65 db.phr:264 db.pin:* db.pjs:* db.pog/ db.pos:53 db.pot:44 db.psq:214 db.ptf:16384 db.pto:16",
    },
    Case {
        name: "s5_phd",
        kind: Kind::P5,
        setup: Setup::Blocked(&["db.phd"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.00.phd:37 db.pdb:32768 db.phd/ db.phi:65 db.phr:264 db.pin:* db.pjs:* db.pog:44 db.pos:53 db.pot:44 db.psq:214 db.ptf:16384 db.pto:16",
    },
    Case {
        name: "s5_phr",
        kind: Kind::P5,
        setup: Setup::Blocked(&["db.phr"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.00.phr:264 db.pdb:32768 db.phd:37 db.phi:65 db.phr/ db.pin:* db.pjs:* db.pog:44 db.pos:53 db.pot:44 db.psq:214 db.ptf:16384 db.pto:16",
    },
    Case {
        name: "s5_psq",
        kind: Kind::P5,
        setup: Setup::Blocked(&["db.psq"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.00.psq:214 db.pdb:32768 db.phd:37 db.phi:65 db.phr:264 db.pin:* db.pjs:* db.pog:44 db.pos:53 db.pot:44 db.psq/ db.ptf:16384 db.pto:16",
    },
    Case {
        name: "s4_pin",
        kind: Kind::P4,
        setup: Setup::Blocked(&["db.pin"]),
        status: 2,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "BLAST Database error: No alias or index file found for protein database [<D>/db] in search path [<C>::]. Please verify the spelling of the BLAST database and its molecule type.\n",
        files: "db.00.pin:* db.phd:37 db.phi:65 db.phr:264 db.pin/ db.pog:44 db.psd:63 db.psi:64 db.psq:214",
    },
    Case {
        name: "s4_psd",
        kind: Kind::P4,
        setup: Setup::Blocked(&["db.psd"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.00.psd:63 db.phd:37 db.phi:65 db.phr:264 db.pin:* db.pjs:* db.pog:44 db.psd/ db.psi:64 db.psq:214",
    },
    Case {
        name: "s4_psi",
        kind: Kind::P4,
        setup: Setup::Blocked(&["db.psi"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.00.psi:64 db.phd:37 db.phi:65 db.phr:264 db.pin:* db.pjs:* db.pog:44 db.psd:63 db.psi/ db.psq:214",
    },
    Case {
        name: "s4_phi",
        kind: Kind::P4,
        setup: Setup::Blocked(&["db.phi"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.00.phi:65 db.phd:37 db.phi/ db.phr:264 db.pin:* db.pjs:* db.pog:44 db.psd:63 db.psi:64 db.psq:214",
    },
    Case {
        name: "n5_nin",
        kind: Kind::N5,
        setup: Setup::Blocked(&["db.nin"]),
        status: 2,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Nucleotide\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 2 sequences in <S> seconds.\n\n\n",
        stderr: "BLAST Database error: No alias or index file found for nucleotide database [<D>/db] in search path [<C>::]. Please verify the spelling of the BLAST database and its molecule type.\n",
        files: "db.00.nin:* db.ndb:32768 db.nhd:25 db.nhi:65 db.nhr:97 db.nin/ db.nog:40 db.nos:34 db.not:32 db.nsq:37 db.ntf:16384 db.nto:12",
    },
    Case {
        name: "n5_nsq",
        kind: Kind::N5,
        setup: Setup::Blocked(&["db.nsq"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Nucleotide\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 2 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.00.nsq:37 db.ndb:32768 db.nhd:25 db.nhi:65 db.nhr:97 db.nin:* db.njs:* db.nog:40 db.nos:34 db.not:32 db.nsq/ db.ntf:16384 db.nto:12",
    },
    Case {
        name: "n5_nog",
        kind: Kind::N5,
        setup: Setup::Blocked(&["db.nog"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Nucleotide\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 2 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.00.nog:40 db.ndb:32768 db.nhd:25 db.nhi:65 db.nhr:97 db.nin:* db.njs:* db.nog/ db.nos:34 db.not:32 db.nsq:37 db.ntf:16384 db.nto:12",
    },
    Case {
        name: "n4_nin",
        kind: Kind::N4,
        setup: Setup::Blocked(&["db.nin"]),
        status: 2,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Nucleotide\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 2 sequences in <S> seconds.\n\n\n",
        stderr: "BLAST Database error: No alias or index file found for nucleotide database [<D>/db] in search path [<C>::]. Please verify the spelling of the BLAST database and its molecule type.\n",
        files: "db.00.nin:* db.nhd:25 db.nhi:65 db.nhr:97 db.nin/ db.nog:40 db.nsd:36 db.nsi:63 db.nsq:37",
    },
    Case {
        name: "n4_nsd",
        kind: Kind::N4,
        setup: Setup::Blocked(&["db.nsd"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Nucleotide\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 2 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.00.nsd:36 db.nhd:25 db.nhi:65 db.nhr:97 db.nin:* db.njs:* db.nog:40 db.nsd/ db.nsi:63 db.nsq:37",
    },
    Case {
        name: "s5plain_pin",
        kind: Kind::Plain,
        setup: Setup::Blocked(&["db.pin"]),
        status: 2,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "BLAST Database error: No alias or index file found for protein database [<D>/db] in search path [<C>::]. Please verify the spelling of the BLAST database and its molecule type.\n",
        files: "db.00.pin:* db.pdb:20480 db.phr:279 db.pin/ db.pot:44 db.psq:214 db.ptf:16384 db.pto:16",
    },
    Case {
        name: "s5plain_psq_phr",
        kind: Kind::Plain,
        setup: Setup::Blocked(&["db.psq", "db.phr"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.00.phr:279 db.00.psq:214 db.pdb:20480 db.phr/ db.pin:* db.pjs:* db.pot:44 db.psq/ db.ptf:16384 db.pto:16",
    },
    Case {
        name: "st_pnd_file",
        kind: Kind::P5,
        setup: Setup::Stale(&["db.00.pnd", "db.00.pni"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.pdb:32768 db.phd:37 db.phi:65 db.phr:264 db.pin:* db.pjs:* db.pnd:15 db.pni:15 db.pog:44 db.pos:53 db.pot:44 db.psq:214 db.ptf:16384 db.pto:16",
    },
    Case {
        name: "st_ppd_pi",
        kind: Kind::P5,
        setup: Setup::Stale(&["db.00.ppd", "db.00.ppi"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.pdb:32768 db.phd:37 db.phi:65 db.phr:264 db.pin:* db.pjs:* db.pog:44 db.pos:53 db.pot:44 db.ppd:15 db.ppi:15 db.psq:214 db.ptf:16384 db.pto:16",
    },
    Case {
        name: "st_ptd_ti",
        kind: Kind::P5,
        setup: Setup::Stale(&["db.00.ptd", "db.00.pti"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.00.ptd:15 db.00.pti:15 db.pdb:32768 db.phd:37 db.phi:65 db.phr:264 db.pin:* db.pjs:* db.pog:44 db.pos:53 db.pot:44 db.psq:214 db.ptf:16384 db.pto:16",
    },
    Case {
        name: "st_psd_v5",
        kind: Kind::P5,
        setup: Setup::Stale(&["db.00.psd", "db.00.psi"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.00.psd:15 db.00.psi:15 db.pdb:32768 db.phd:37 db.phi:65 db.phr:264 db.pin:* db.pjs:* db.pog:44 db.pos:53 db.pot:44 db.psq:214 db.ptf:16384 db.pto:16",
    },
    Case {
        name: "st_psd_v4",
        kind: Kind::P4,
        setup: Setup::Stale(&["db.00.psd", "db.00.psi"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.phd:37 db.phi:65 db.phr:264 db.pin:* db.pjs:* db.pog:44 db.psd:63 db.psi:64 db.psq:214",
    },
    Case {
        name: "st_nt",
        kind: Kind::N5,
        setup: Setup::Stale(&["db.00.ntd", "db.00.nti", "db.00.npd", "db.00.npi"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Nucleotide\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 2 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.00.npd:15 db.00.npi:15 db.00.ntd:15 db.00.nti:15 db.ndb:32768 db.nhd:25 db.nhi:65 db.nhr:97 db.nin:* db.njs:* db.nog:40 db.nos:34 db.not:32 db.nsq:37 db.ntf:16384 db.nto:12",
    },
    Case {
        name: "st_phd_hashless",
        kind: Kind::V5,
        setup: Setup::Stale(&["db.00.phd", "db.00.phi", "db.00.pnd"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.00.phd:15 db.00.phi:15 db.pdb:32768 db.phr:264 db.pin:* db.pjs:* db.pnd:15 db.pog:44 db.pos:53 db.pot:44 db.psq:214 db.ptf:16384 db.pto:16",
    },
    Case {
        name: "st_pog_file",
        kind: Kind::Plain,
        setup: Setup::Stale(&["db.00.pog"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from FASTA; added 3 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "db.00.pog:15 db.pdb:20480 db.phr:279 db.pin:* db.pjs:* db.pot:44 db.psq:214 db.ptf:16384 db.pto:16",
    },
    Case {
        name: "L_pdb_dir",
        kind: Kind::V5,
        setup: Setup::Blocked(&["db.pdb"]),
        status: 2,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\n\nNo volumes were created.\n\n",
        stderr: "BLAST Database error: LMDB runtime error: mdb_env_open: Is a directory\n",
        files: "db.pdb/ db.pdb-lock:8192",
    },
    Case {
        name: "L_ptf_dir",
        kind: Kind::V5,
        setup: Setup::Blocked(&["db.ptf"]),
        status: 2,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\n\nNo volumes were created.\n\n",
        stderr: "BLAST Database error: LMDB runtime error: mdb_env_open: Is a directory\n",
        files: "db.pdb:8192 db.ptf/ db.ptf-lock:8192",
    },
    Case {
        name: "L_pdblock_dir",
        kind: Kind::V5,
        setup: Setup::Blocked(&["db.pdb-lock"]),
        status: 2,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\n\nNo volumes were created.\n\n",
        stderr: "BLAST Database error: LMDB runtime error: mdb_env_open: Is a directory\n",
        files: "db.pdb-lock/",
    },
    Case {
        name: "L_ptflock_dir",
        kind: Kind::V5,
        setup: Setup::Blocked(&["db.ptf-lock"]),
        status: 2,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\n\nNo volumes were created.\n\n",
        stderr: "BLAST Database error: LMDB runtime error: mdb_env_open: Is a directory\n",
        files: "db.pdb:8192 db.ptf-lock/",
    },
    Case {
        name: "m_000pin",
        kind: Kind::Many,
        setup: Setup::Blocked(&["db.000.pin"]),
        status: 2,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 0B\nAdding sequences from FASTA; added 130 sequences in <S> seconds.\n\n\n",
        stderr: "BLAST Database error: Could not find volume or alias file (db.000) referenced in alias file (<D>/db).\n",
        files: "+db.00.pin:* ~db.000.pin/ -db.pjs",
    },
    Case {
        name: "m_000psq",
        kind: Kind::Many,
        setup: Setup::Blocked(&["db.000.psq"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 0B\nAdding sequences from FASTA; added 130 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "+db.00.psq:62 ~db.000.psq/",
    },
    Case {
        name: "m_000pog",
        kind: Kind::Many,
        setup: Setup::Blocked(&["db.000.pog"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 0B\nAdding sequences from FASTA; added 130 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "+db.00.pog:36 ~db.000.pog/",
    },
    Case {
        name: "m_050pin",
        kind: Kind::Many,
        setup: Setup::Blocked(&["db.050.pin"]),
        status: 2,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 0B\nAdding sequences from FASTA; added 130 sequences in <S> seconds.\n\n\n",
        stderr: "BLAST Database error: Could not find volume or alias file (db.050) referenced in alias file (<D>/db).\n",
        files: "~db.050.pin/ +db.50.pin:* -db.pjs",
    },
    Case {
        name: "m_099phr",
        kind: Kind::Many,
        setup: Setup::Blocked(&["db.099.phr"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 0B\nAdding sequences from FASTA; added 130 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "~db.099.phr/ +db.99.phr:79",
    },
    Case {
        name: "m_000pin_050psq",
        kind: Kind::Many,
        setup: Setup::Blocked(&["db.000.pin", "db.050.psq"]),
        status: 2,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 0B\nAdding sequences from FASTA; added 130 sequences in <S> seconds.\n\n\n",
        stderr: "BLAST Database error: Could not find volume or alias file (db.000) referenced in alias file (<D>/db).\n",
        files: "+db.00.pin:* ~db.000.pin/ ~db.050.psq/ +db.50.psq:62 -db.pjs",
    },
    Case {
        name: "mh_000phd",
        kind: Kind::ManyHash,
        setup: Setup::Blocked(&["db.000.phd"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 0B\nAdding sequences from FASTA; added 130 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "+db.00.phd:13 ~db.000.phd/",
    },
    Case {
        name: "mh_000phi",
        kind: Kind::ManyHash,
        setup: Setup::Blocked(&["db.000.phi"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 0B\nAdding sequences from FASTA; added 130 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "+db.00.phi:65 ~db.000.phi/",
    },
    Case {
        name: "st_m_pog",
        kind: Kind::Many,
        setup: Setup::Stale(&["db.00.pnd", "db.050.pnd", "db.100.pnd"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 0B\nAdding sequences from FASTA; added 130 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "+db.000.pnd:15 +db.050.pnd:16 +db.100.pnd:16",
    },
    Case {
        name: "st_m_ppd",
        kind: Kind::Many,
        setup: Setup::Stale(&["db.00.ppd", "db.00.ppi", "db.07.ppd", "db.130.ppd"]),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 0B\nAdding sequences from FASTA; added 130 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "+db.000.ppd:15 +db.000.ppi:15 +db.007.ppd:15 +db.130.ppd:16",
    },
    Case {
        name: "st_m_pnd_blocked",
        kind: Kind::Many,
        setup: Setup::StaleAndBlocked("db.00.pnd", "db.000.pnd"),
        status: 0,
        stdout: "\n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  in.fa\nSequence type: Protein\nKeep MBits: T\nMaximum file size: 0B\nAdding sequences from FASTA; added 130 sequences in <S> seconds.\n\n\n",
        stderr: "",
        files: "+db.00.pnd:15 +db.000.pnd/",
    },
];

struct TempDir(PathBuf);

impl TempDir {
    fn new(name: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("losat-mkdb-rename-{name}-{}", std::process::id()));
        make_writable(&dir);
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("d")).unwrap();
        TempDir(dir)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        make_writable(&self.0);
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn make_writable(path: &Path) {
    use std::os::unix::fs::PermissionsExt;
    if std::fs::symlink_metadata(path).is_ok_and(|meta| meta.is_dir()) {
        let _ = std::fs::set_permissions(path, std::fs::Permissions::from_mode(0o755));
        if let Ok(entries) = std::fs::read_dir(path) {
            for entry in entries.map_while(Result::ok) {
                make_writable(&entry.path());
            }
        }
    }
}

/// A directory at `path` that holds a directory that cannot be read.
fn blocked(path: &Path) {
    use std::os::unix::fs::PermissionsExt;
    std::fs::create_dir_all(path.join("u").join("x")).unwrap();
    std::fs::write(path.join("u").join("x").join("f"), b"x").unwrap();
    std::fs::set_permissions(path.join("u"), std::fs::Permissions::from_mode(0o000)).unwrap();
}

fn many() -> String {
    // 130 proteins of 60 letters; the letters do not matter, the volumes do.
    let letters = b"ACDEFGHIKLMNPQRSTVWY";
    let mut input = String::new();
    for i in 0..130usize {
        input.push_str(&format!(">sp|P{i:05}|S{i} desc {i}\n"));
        for j in 0..60usize {
            input.push(letters[(i * 7 + j * 3) % letters.len()] as char);
        }
        input.push('\n');
    }
    input
}

struct Outcome {
    status: i32,
    stdout: String,
    stderr: String,
    /// Top-level entries of the output directory.
    entries: BTreeMap<String, Entry>,
}

#[derive(Clone, PartialEq, Eq, Debug)]
enum Entry {
    File(u64),
    Dir,
}

fn run(case_name: &str, kind: Kind, setup: Option<&Setup>) -> Outcome {
    let tmp = TempDir::new(case_name);
    let d = tmp.0.join("d");
    let (input, args): (String, &[&str]) = match kind {
        Kind::P5 => (
            PROT.into(),
            &["-dbtype", "prot", "-parse_seqids", "-hash_index"],
        ),
        Kind::P4 => (
            PROT.into(),
            &[
                "-dbtype",
                "prot",
                "-parse_seqids",
                "-hash_index",
                "-blastdb_version",
                "4",
            ],
        ),
        Kind::N5 => (
            NUC.into(),
            &["-dbtype", "nucl", "-parse_seqids", "-hash_index"],
        ),
        Kind::N4 => (
            NUC.into(),
            &[
                "-dbtype",
                "nucl",
                "-parse_seqids",
                "-hash_index",
                "-blastdb_version",
                "4",
            ],
        ),
        Kind::Plain => (PROT.into(), &["-dbtype", "prot"]),
        Kind::V5 => (PROT.into(), &["-dbtype", "prot", "-parse_seqids"]),
        Kind::Many => (
            many(),
            &["-dbtype", "prot", "-parse_seqids", "-max_file_sz", "100B"],
        ),
        Kind::ManyHash => (
            many(),
            &[
                "-dbtype",
                "prot",
                "-parse_seqids",
                "-max_file_sz",
                "100B",
                "-hash_index",
            ],
        ),
    };
    std::fs::write(tmp.0.join("in.fa"), input).unwrap();
    match setup {
        None => {}
        Some(Setup::Blocked(names)) => names.iter().for_each(|n| blocked(&d.join(n))),
        Some(Setup::Stale(names)) => names
            .iter()
            .for_each(|n| std::fs::write(d.join(n), format!("stale {n}")).unwrap()),
        Some(Setup::StaleAndBlocked(stale, block)) => {
            std::fs::write(d.join(stale), format!("stale {stale}")).unwrap();
            blocked(&d.join(block));
        }
    }
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    command
        .current_dir(&tmp.0)
        .arg("makeblastdb")
        .args(["-in", "in.fa"])
        .args(args)
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
        entries.insert(
            entry.file_name().to_string_lossy().into_owned(),
            if meta.is_dir() {
                Entry::Dir
            } else {
                Entry::File(meta.len())
            },
        );
    }
    // The directory that was made unreadable has its write permission added by makeblastdb's
    // failed removal (`CDir::Remove`), and nothing else (probe of 2026-11-11).
    if let Some(Setup::Blocked(names)) = setup {
        for name in *names {
            use std::os::unix::fs::PermissionsExt;
            let mode = std::fs::metadata(d.join(name).join("u"))
                .map(|meta| meta.permissions().mode() & 0o777);
            assert_eq!(mode.ok(), Some(0o200), "{case_name}: {name}");
        }
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

fn is_many(kind: Kind) -> bool {
    matches!(kind, Kind::Many | Kind::ManyHash)
}

/// Root reads any directory: the cases that need an unreadable one say nothing then.
fn running_as_root(tag: &str) -> bool {
    let tmp = TempDir::new(&format!("root-check-{tag}"));
    blocked(&tmp.0);
    std::fs::read_dir(tmp.0.join("u")).is_ok()
}

#[test]
fn renames_that_fail_and_lock_files_end_as_makeblastdb_ends() {
    if running_as_root("table") {
        return;
    }
    let mut failures = Vec::new();
    for case in CASES {
        let outcome = run(case.name, case.kind, Some(&case.setup));
        let files = if is_many(case.kind) {
            let clean = run(case.name, case.kind, None);
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

/// The lock files makeblastdb leaves are the bytes of a lock region no transaction has used
/// (`mdb_env_setup_locks` on Linux x86-64), 8192 bytes, and the data files beside them are the
/// size makeblastdb leaves.
#[test]
fn the_lock_files_makeblastdb_leaves_hold_the_bytes_of_an_unused_lock_region() {
    if running_as_root("bytes") {
        return;
    }
    let mut expected = vec![0u8; 8192];
    expected[..8].copy_from_slice(&[0xDE, 0xC0, 0xEF, 0xBE, 0x01, 0x00, 0x01, 0x00]);
    expected[0x18] = 0x90;
    expected[0x50] = 0x90;
    for (name, lock) in [("L_pdb_dir", "db.pdb-lock"), ("L_ptf_dir", "db.ptf-lock")] {
        let case = CASES.iter().find(|case| case.name == name).unwrap();
        let tmp_name = format!("{name}-bytes");
        let tmp = TempDir::new(&tmp_name);
        // The same run as `run`, kept to read the file.
        let d = tmp.0.join("d");
        std::fs::write(tmp.0.join("in.fa"), PROT).unwrap();
        let Setup::Blocked(names) = &case.setup else {
            unreachable!()
        };
        names.iter().for_each(|n| blocked(&d.join(n)));
        let status = Command::new(env!("CARGO_BIN_EXE_losat"))
            .current_dir(&tmp.0)
            .args([
                "makeblastdb",
                "-in",
                "in.fa",
                "-dbtype",
                "prot",
                "-parse_seqids",
                "-out",
            ])
            .arg(d.join("db"))
            .env_remove("BLASTDB")
            .output()
            .unwrap()
            .status
            .code();
        assert_eq!(status, Some(2), "{name}");
        assert_eq!(std::fs::read(d.join(lock)).unwrap(), expected, "{name}");
    }
}
