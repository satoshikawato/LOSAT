//! The tools' posts to a `-logfile` that is a pipe whose reader has closed it, against NCBI BLAST+
//! 2.17.0 (the local official release, 2026-11-28, line RG item T8d, probe
//! `remaining-gaps-20261112/T/t8c/t8d`): `-logfile /dev/stdout` with the standard output such a
//! pipe (`Log::Stdout`, the standard error kept) and `-logfile /dev/stderr` with the standard error
//! such a pipe (`Log::Stderr`, the standard output kept); and, without `-logfile`, the posts and
//! the other texts on such a standard error (`Log::Cerr`; line RG item T8c round 3, probe
//! `remaining-gaps-20261112/T/t8c/r3`, the same release, 2026-11-28), blastdb_aliastool's and
//! blastdbcheck's included (the Owner's answer 「Q-T8d-W=A」 of 2026-11-28).
//!
//! NCBI's handler of the `-logfile` file writes each post at once and does not look at what the
//! write returns, so the first post that meets the closed pipe ends NCBI by SIGPIPE (exit status
//! 141 from the shell, no message): the messages NCBI collected before it opened the file (the
//! usage report's), the configuration's, and the run's warnings and errors. LOSAT ends there with
//! exit status 0 and no message: the output NCBI wrote before the post (its flushes: blastdbcmd's
//! `endl`, each ASN.1 object, the end of the FASTA formatter) is written, and what NCBI's buffers
//! still held is lost, as NCBI's ended process loses it (blastdb_aliastool's GI list file: empty)
//! (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 10). Without `-logfile` the posts go to `cerr`,
//! and so do the reports that segmasker and convert2blastmask write to it themselves, makembindex's
//! progress and the usage report's messages before the arguments: the first that meets the closed
//! pipe ends NCBI there, and LOSAT there. An argument
//! error keeps its exit status 1. windowmasker with `-infmt seqids` runs on to its exit status 1:
//! NCBI's GenBank loader, which it registers first, ignores SIGPIPE.
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:6683-6697
//! NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:144-147,202-204
//! NCBI reference (598d8ae6): c++/src/connect/ncbi_socket.c:866-873

#![cfg(unix)]

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

mod temp_support;
use temp_support::unique_dir;

/// Which stream `-logfile` names, the one that is the closed pipe.
#[derive(Clone, Copy, Debug)]
enum Log {
    /// `-logfile /dev/stdout`, the standard output closed; the standard error is kept.
    Stdout,
    /// `-logfile /dev/stderr`, the standard error closed; the standard output is kept.
    Stderr,
    /// No `-logfile`: the diagnostics go to `cerr`, the standard error, which is closed; the
    /// standard output is kept.
    Cerr,
}

/// Where a case runs.
#[derive(Clone, Copy, Debug)]
enum Place {
    /// The program's directory of fixtures (`tests/<name>`), which the run does not write.
    Fixtures(&'static str),
    /// A new directory holding a copy of `tests/<path>` as `in`.
    Fresh(&'static str),
    /// A new directory whose `.ncbirc` (`NCBI_CONFIG_PATH`) is malformed.
    MalformedNcbirc,
    /// The case directory `cases/<dir>` of the blastdbcheck fixtures, built in a new directory by
    /// their `setup.tsv` (`blastdbcheck_fixtures::build_case_dir`).
    CheckCase(&'static str),
}

/// One run: NCBI's exit status (141 is SIGPIPE, which LOSAT gives as 0), NCBI's kept stream, and
/// the files NCBI left (name, size, FNV-1a 64; LMDB's data files by name and size only, a
/// written `.pin` by name only: its size follows the length of the build's date).
struct Case {
    name: &'static str,
    program: &'static str,
    place: Place,
    words: &'static [&'static str],
    variables: &'static [(&'static str, &'static str)],
    log: Log,
    status: i32,
    kept: Kept,
    files: &'static [(&'static str, u64, u64)],
}

/// NCBI's kept stream.
enum Kept {
    Bytes(&'static [u8]),
    /// NCBI's text starts so (the usage of an argument error, which LOSAT writes in its words).
    StartsWith(&'static [u8]),
    /// The first bytes of NCBI's standard output of the row's run without `-logfile`, a file of
    /// `tests` (path, length).
    Prefix(&'static str, usize),
    /// The first lines of NCBI's standard output of the row's run without a closed pipe, a file of
    /// `tests` whose `<RUN>` is the run's directory (path, number of lines).
    Lines(&'static str, usize),
}

const EMPTY: u64 = 0xcbf2_9ce4_8422_2325;

/// The point where NCBI's run ends is named in each comment.
const CASES: &[Case] = &[
    // blastdbcmd: the error of a database that is not there.
    Case {
        name: "missing_db",
        program: "blastdbcmd",
        place: Place::Fixtures("blastdbcmd"),
        words: &["-db", "nonexistent_db", "-entry", "x"],
        variables: &[],
        log: Log::Stdout,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    Case {
        name: "missing_db",
        program: "blastdbcmd",
        place: Place::Fixtures("blastdbcmd"),
        words: &["-db", "nonexistent_db", "-entry", "x"],
        variables: &[],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    // blastdbcmd: the error of a duplicate Seq-id of the FASTA formatter, after the
    // destructor of its `CFastaOstream` flushed the record before it.
    Case {
        name: "D90",
        program: "blastdbcmd",
        place: Place::Fixtures("blastdbcmd"),
        words: &[
            "-db",
            "db/n4n",
            "-entry",
            "gnl|BL_ORD_ID|1,gnl|BL_ORD_ID|1",
            "-long_seqids",
        ],
        variables: &[],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Prefix("blastdbcmd/expected/D90.stdout", 559),
        files: &[],
    },
    // blastdbcmd `%d`: the error of an id that is not found, after the objects before it, each
    // flushed by the end of its object stream (3,658 of 7,879 bytes).
    Case {
        name: "Q123",
        program: "blastdbcmd",
        place: Place::Fixtures("blastdbcmd"),
        words: &[
            "-db",
            "db/p5",
            "-entry_batch",
            "in/pids.txt",
            "-outfmt",
            "%d",
        ],
        variables: &[],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Prefix("blastdbcmd/expected/Q123.stdout", 3658),
        files: &[],
    },
    // blastdbcmd: the error of an id that is not found, between found entries; the entries
    // before it are still in NCBI's buffer of the standard output, which the end loses.
    Case {
        name: "skip_mid",
        program: "blastdbcmd",
        place: Place::Fixtures("blastdbcmd"),
        words: &["-db", "db/p4", "-entry_batch", "in/pids.txt"],
        variables: &[],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    // blastdbcmd: the usage report's message of a malformed `.ncbirc`, written into the file
    // as NCBI opens it, before any entry.
    Case {
        name: "ncbirc_log_setup",
        program: "blastdbcmd",
        place: Place::MalformedNcbirc,
        words: &["-db", "{F}/blastdbcmd/db/nm5", "-entry", "all"],
        variables: &[],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    // blastdbcmd: the same message of the configuration's load (no usage report).
    Case {
        name: "ncbirc_config_load",
        program: "blastdbcmd",
        place: Place::MalformedNcbirc,
        words: &["-db", "{F}/blastdbcmd/db/nm5", "-entry", "all"],
        variables: &[("BLAST_USAGE_REPORT", "false")],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    // dustmasker: the reader's error of a record, after the record before it.
    Case {
        name: "db_t5d_seqid",
        program: "dustmasker",
        place: Place::Fixtures("dustmasker"),
        words: &["-infmt", "blastdb", "-in", "../db_search/t5/nhd_seqid"],
        variables: &[],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Bytes(b">record 0 of the T5d databases\n"),
        files: &[],
    },
    Case {
        name: "db_t5d_seqid",
        program: "dustmasker",
        place: Place::Fixtures("dustmasker"),
        words: &["-infmt", "blastdb", "-in", "../db_search/t5/nhd_seqid"],
        variables: &[],
        log: Log::Stdout,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    // dustmasker: a `-conffile` that NCBI cannot read (exit status 2 otherwise) or that is not
    // there (255).
    Case {
        name: "conf_malformed",
        program: "dustmasker",
        place: Place::Fixtures("dustmasker"),
        words: &["-in", "in/basic.fa", "-conffile", "conf/malformed.ini"],
        variables: &[],
        log: Log::Stdout,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    Case {
        name: "conf_missing",
        program: "dustmasker",
        place: Place::Fixtures("dustmasker"),
        words: &["-in", "in/basic.fa", "-conffile", "conf/missing.ini"],
        variables: &[],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    // segmasker: the taxonomy database's error of a record's title; the records before it are
    // written (762 of 2,695 bytes).
    Case {
        name: "db_t5b_offset_last_past_end_fasta",
        program: "segmasker",
        place: Place::Fixtures("segmasker"),
        words: &[
            "-infmt",
            "blastdb",
            "-in",
            "../db_search/t5b/pmv",
            "-outfmt",
            "fasta",
            "-conffile",
            "conf/t5b_offset_last_past_end.ini",
        ],
        variables: &[],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Prefix("segmasker/ncbi/db_t5b_offset_last_past_end_fasta.out", 762),
        files: &[],
    },
    // windowmasker: the configuration's error of an id of `-exclude_ids` that `CSeq_id` cannot
    // read, before any output.
    Case {
        name: "s_ps_ex_db",
        program: "windowmasker",
        place: Place::Fixtures("windowmasker"),
        words: &[
            "-ustat",
            "../winmask/counts/dm2L_u9.binary",
            "-infmt",
            "blastdb",
            "-in",
            "db/wm_ps",
            "-exclude_ids",
            "in/ids_db.txt",
            "-text_match",
            "F",
        ],
        variables: &[],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    // windowmasker: the counts' warning of a `-t_low` below the stored one, before any output.
    Case {
        name: "t_low20",
        program: "windowmasker",
        place: Place::Fixtures("windowmasker"),
        words: &[
            "-ustat",
            "../winmask/counts/dm2L_u9.binary",
            "-in",
            "../winmask/fasta/dm2L_60k.fa",
            "-t_low",
            "20",
        ],
        variables: &[],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    // windowmasker: the reader's warning of a record; the intervals before it are written.
    Case {
        name: "d_letters",
        program: "windowmasker",
        place: Place::Fixtures("windowmasker"),
        words: &[
            "-ustat",
            "../winmask/counts/dm2L_u9.binary",
            "-in",
            "in/letters.fa",
            "-dust",
            "T",
        ],
        variables: &[],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Bytes(
            b">t_lower some description\n248 - 261\n498 - 512\n1015 - 1029\n2431 - 2445\n\
              2773 - 2785\n3318 - 3331\n3358 - 3406\n3540 - 3560\n3850 - 3864\n4502 - 4516\n\
              4634 - 4654\n4871 - 4883\n4993 - 5068\n5218 - 5239\n5273 - 5312\n5357 - 5419\n\
              5848 - 5861\n5891 - 5904\n>t_nrun\n",
        ),
        files: &[],
    },
    // windowmasker `-infmt seqids`: SIGPIPE ignored, the error post lost, exit status 1.
    Case {
        name: "m_seqids",
        program: "windowmasker",
        place: Place::Fixtures("windowmasker"),
        words: &[
            "-ustat",
            "../winmask/counts/dm2L_u9.binary",
            "-infmt",
            "seqids",
            "-in",
            "in/ids_a.txt",
        ],
        variables: &[],
        log: Log::Stdout,
        status: 1,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    Case {
        name: "m_seqids_mi",
        program: "windowmasker",
        place: Place::Fixtures("windowmasker"),
        words: &[
            "-ustat",
            "../winmask/counts/dm2L_u9.binary",
            "-infmt",
            "seqids",
            "-in",
            "in/ids_a.txt",
            "-outfmt",
            "maskinfo_asn1_text",
        ],
        variables: &[],
        log: Log::Stderr,
        status: 1,
        kept: Kept::Bytes(
            b"Blast-db-mask-info ::= {\n  algo-id 3,\n  algo-program windowmasker,\n  \
              algo-options \"\",\n  masks {\n    masks {\n    },\n    more FALSE\n  }\n}\n",
        ),
        files: &[],
    },
    // windowmasker: an argument error keeps its exit status 1 and its usage.
    Case {
        name: "badopt",
        program: "windowmasker",
        place: Place::Fixtures("windowmasker"),
        words: &[
            "-ustat",
            "../winmask/counts/dm2L_u9.binary",
            "-nonexistent_option",
        ],
        variables: &[],
        log: Log::Stdout,
        status: 1,
        kept: Kept::StartsWith(b"USAGE\n  windowmasker "),
        files: &[],
    },
    // convert2blastmask: the serial stream's error of a byte of an id, in the first object
    // (27,622 bytes without `-logfile`).
    Case {
        name: "order_txt",
        program: "convert2blastmask",
        place: Place::Fixtures("convert2blastmask"),
        words: &[
            "-in",
            "in/order_del_lc.fa",
            "-masking_algorithm",
            "dust",
            "-masking_options",
            "o",
            "-parse_seqids",
        ],
        variables: &[],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    Case {
        name: "conf_malformed",
        program: "convert2blastmask",
        place: Place::Fixtures("convert2blastmask"),
        words: &[
            "-conffile",
            "conf/malformed.ini",
            "-in",
            "in/lc.fa",
            "-masking_algorithm",
            "dust",
            "-masking_options",
            "o",
        ],
        variables: &[],
        log: Log::Stdout,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    // makembindex: the error of a database that is not there.
    Case {
        name: "missing_db",
        program: "makembindex",
        place: Place::Fresh("makembindex/in"),
        words: &[
            "-input",
            "nonexistent_db",
            "-output",
            "idx",
            "-iformat",
            "blastdb",
        ],
        variables: &[],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    // makeprofiledb: the argument error of an `-in` that `Run` cannot open (`AsInputFile`), a
    // `CArgException` reported with the usage: NCBI's post meets the closed pipe after the usage,
    // and the exit status 1 of an argument error stays (version 1.19).
    Case {
        name: "H10",
        program: "makeprofiledb",
        place: Place::Fresh("fixtures/makeprofiledb/in"),
        words: &["-in", "in/nonexistent.lst", "-out", "db"],
        variables: &[],
        log: Log::Stdout,
        status: 1,
        kept: Kept::StartsWith(b"USAGE\n  makeprofiledb "),
        files: &[],
    },
    // blastdbcmd: `CATCH_ALL`'s argument error of an `-out` that `Run` cannot open, exit status 1.
    Case {
        name: "B58",
        program: "blastdbcmd",
        place: Place::Fixtures("blastdbcmd"),
        words: &["-db", "db/p5", "-info", "-out", ""],
        variables: &[],
        log: Log::Stderr,
        status: 1,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    // makeprofiledb: the error of the profile after the volume's close; the volume stays, the
    // RPS files hold the profile's rows, `.aux` what its buffer passed on (nothing).
    Case {
        name: "QD36",
        program: "makeprofiledb",
        place: Place::Fresh("fixtures/makeprofiledb/in"),
        words: &["-in", "in/l_QD36.lst", "-out", "db"],
        variables: &[],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Bytes(b"ib=1,n4sz=0\n"),
        files: &[
            ("db.aux", 0, EMPTY),
            ("db.freq", 16, 0x05c6_e401_c1af_5098),
            ("db.loo", 0, EMPTY),
            ("db.pdb", 32768, 0),
            ("db.phr", 642, 0xa75f_d20b_fa07_f754),
            ("db.pin", 0, 0),
            ("db.pog", 36, 0x4501_919c_e010_9275),
            ("db.pos", 27, 0xdb9a_1f64_7a43_b26d),
            ("db.pot", 20, 0x4192_1acc_144d_4237),
            ("db.psq", 3, 0xd796_ac18_6a9a_b80e),
            ("db.ptf", 16384, 0),
            ("db.pto", 8, 0x89cd_3129_1d2a_efa4),
            ("db.rps", 16, 0x05c6_e401_c1af_5098),
        ],
    },
    // Without `-logfile`, on a closed standard error (line RG item T8c round 3). segmasker: the
    // report that its `Run` writes to `cerr` itself, after the records before the error (986 of
    // 2,695 bytes).
    Case {
        name: "db_t5b_offset_last_past_end_fasta",
        program: "segmasker",
        place: Place::Fixtures("segmasker"),
        words: &[
            "-infmt",
            "blastdb",
            "-in",
            "../db_search/t5b/pmv",
            "-outfmt",
            "fasta",
            "-conffile",
            "conf/t5b_offset_last_past_end.ini",
        ],
        variables: &[],
        log: Log::Cerr,
        status: 141,
        kept: Kept::Prefix("segmasker/ncbi/db_t5b_offset_last_past_end_fasta.out", 986),
        files: &[],
    },
    // dustmasker: the error post of the record after the first.
    Case {
        name: "db_t5d_seqid",
        program: "dustmasker",
        place: Place::Fixtures("dustmasker"),
        words: &["-infmt", "blastdb", "-in", "../db_search/t5/nhd_seqid"],
        variables: &[],
        log: Log::Cerr,
        status: 141,
        kept: Kept::Bytes(b">record 0 of the T5d databases\n"),
        files: &[],
    },
    // windowmasker: the warning after the records before it (470 of 765 bytes).
    Case {
        name: "d_letters",
        program: "windowmasker",
        place: Place::Fixtures("windowmasker"),
        words: &[
            "-ustat",
            "../winmask/counts/dm2L_u9.binary",
            "-in",
            "in/letters.fa",
            "-dust",
            "T",
        ],
        variables: &[],
        log: Log::Cerr,
        status: 141,
        kept: Kept::Prefix("windowmasker/ncbi/d_letters.out", 470),
        files: &[],
    },
    // convert2blastmask: the report of the exception that its `Run` writes to `cerr` itself.
    Case {
        name: "order_txt",
        program: "convert2blastmask",
        place: Place::Fixtures("convert2blastmask"),
        words: &[
            "-in",
            "in/order_del_lc.fa",
            "-masking_algorithm",
            "dust",
            "-masking_options",
            "o",
            "-parse_seqids",
        ],
        variables: &[],
        log: Log::Cerr,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    // blastdbcmd `%d`: the error post of an id that is not found, after the objects before it.
    Case {
        name: "Q123",
        program: "blastdbcmd",
        place: Place::Fixtures("blastdbcmd"),
        words: &[
            "-db",
            "db/p5",
            "-entry_batch",
            "in/pids.txt",
            "-outfmt",
            "%d",
        ],
        variables: &[],
        log: Log::Cerr,
        status: 141,
        kept: Kept::Prefix("blastdbcmd/expected/Q123.stdout", 3658),
        files: &[],
    },
    // blastdbcmd: the usage report's message of a malformed `.ncbirc`, before `-version` is read.
    Case {
        name: "ncbirc_usage_report",
        program: "blastdbcmd",
        place: Place::MalformedNcbirc,
        words: &["-version"],
        variables: &[],
        log: Log::Cerr,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    // makembindex: its progress `creating idx.00.idx...`, before the volume's file is made.
    Case {
        name: "progress",
        program: "makembindex",
        place: Place::Fresh("makembindex/in"),
        words: &[
            "-input",
            "in/m.fa",
            "-output",
            "idx",
            "-iformat",
            "fasta",
            "-old_style_index",
            "true",
        ],
        variables: &[("BLAST_USAGE_REPORT", "false")],
        log: Log::Cerr,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    // makeprofiledb: the error post of the profile after the volume's close.
    Case {
        name: "QD36",
        program: "makeprofiledb",
        place: Place::Fresh("fixtures/makeprofiledb/in"),
        words: &["-in", "in/l_QD36.lst", "-out", "db"],
        variables: &[],
        log: Log::Cerr,
        status: 141,
        kept: Kept::Bytes(b"ib=1,n4sz=0\n"),
        files: &[
            ("db.aux", 0, EMPTY),
            ("db.freq", 16, 0x05c6_e401_c1af_5098),
            ("db.loo", 0, EMPTY),
            ("db.pdb", 32768, 0),
            ("db.phr", 642, 0xa75f_d20b_fa07_f754),
            ("db.pin", 0, 0),
            ("db.pog", 36, 0x4501_919c_e010_9275),
            ("db.pos", 27, 0xdb9a_1f64_7a43_b26d),
            ("db.pot", 20, 0x4192_1acc_144d_4237),
            ("db.psq", 3, 0xd796_ac18_6a9a_b80e),
            ("db.ptf", 16384, 0),
            ("db.pto", 8, 0x89cd_3129_1d2a_efa4),
            ("db.rps", 16, 0x05c6_e401_c1af_5098),
        ],
    },
    // blastdb_aliastool (the Owner's answer 「Q-T8d-W=A」): the post of the GI list's conversion,
    // before its file stream closes: the list file is empty, no alias file follows.
    Case {
        name: "A01",
        program: "blastdb_aliastool",
        place: Place::Fresh("fixtures/blastdb_aliastool/in"),
        words: &["-db", "in/p4", "-gilist", "in/gi_ok.txt", "-out", "ga"],
        variables: &[],
        log: Log::Stdout,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[("ga.p.gil", 0, EMPTY)],
    },
    Case {
        name: "A01",
        program: "blastdb_aliastool",
        place: Place::Fresh("fixtures/blastdb_aliastool/in"),
        words: &["-db", "in/p4", "-gilist", "in/gi_ok.txt", "-out", "ga"],
        variables: &[],
        log: Log::Cerr,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[("ga.p.gil", 0, EMPTY)],
    },
    // blastdb_aliastool `-gi_file_in`: the same post, the list file empty.
    Case {
        name: "G02",
        program: "blastdb_aliastool",
        place: Place::Fresh("fixtures/blastdb_aliastool/in"),
        words: &["-gi_file_in", "in/gi_ok.txt", "-gi_file_out", "ok.bgl"],
        variables: &[],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[("ok.bgl", 0, EMPTY)],
    },
    // blastdb_aliastool: the post of the alias file made (its date: by name).
    Case {
        name: "A03",
        program: "blastdb_aliastool",
        place: Place::Fresh("fixtures/blastdb_aliastool/in"),
        words: &["-db", "in/p4", "-gilist", "in/gi_ok.bgl", "-out", "gb"],
        variables: &[],
        log: Log::Cerr,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[("gb.pal", 0, 0)],
    },
    // blastdbcheck: the reader's posts (`Error: (803.7)`) of the first OID's header, after the
    // report's lines before them, each flushed (`endl`).
    Case {
        name: "d_p4_badchar",
        program: "blastdbcheck",
        place: Place::CheckCase("d_p4_badchar"),
        words: &["-db", "p4", "-full", "-verbosity", "4"],
        variables: &[],
        log: Log::Cerr,
        status: 141,
        kept: Kept::Lines("fixtures/blastdbcheck/ncbi/d_p4_badchar.out", 9),
        files: &[],
    },
    // blastdbcheck: the taxonomy database's error post of the MetaData test.
    Case {
        name: "t_bad_magic",
        program: "blastdbcheck",
        place: Place::CheckCase("t_bad_magic"),
        words: &["-db", "p5", "-full", "-verbosity", "3"],
        variables: &[],
        log: Log::Cerr,
        status: 141,
        kept: Kept::Lines("fixtures/blastdbcheck/ncbi/t_bad_magic.out", 8),
        files: &[],
    },
    Case {
        name: "t_bad_magic",
        program: "blastdbcheck",
        place: Place::CheckCase("t_bad_magic"),
        words: &["-db", "p5", "-full", "-verbosity", "3"],
        variables: &[],
        log: Log::Stderr,
        status: 141,
        kept: Kept::Bytes(b""),
        files: &[],
    },
    // blastdbcheck: an argument error keeps exit status 1 where its usage met the closed pipe.
    Case {
        name: "arg_dbtype",
        program: "blastdbcheck",
        place: Place::CheckCase("t"),
        words: &["-db", "../../in/p5", "-dbtype", "foo"],
        variables: &[],
        log: Log::Cerr,
        status: 1,
        kept: Kept::Bytes(b""),
        files: &[],
    },
];

/// FNV-1a 64 of `bytes`.
fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(EMPTY, |hash, &byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

fn tests_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests")
}

/// Copies the directory `from` to `to`, symbolic links as links.
fn copy_tree(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        let kind = entry.file_type().unwrap();
        if kind.is_dir() {
            copy_tree(&entry.path(), &target);
        } else if kind.is_symlink() {
            std::os::unix::fs::symlink(fs::read_link(entry.path()).unwrap(), &target).unwrap();
        } else {
            fs::copy(entry.path(), &target).unwrap();
        }
    }
}

/// A pipe whose reader has closed it.
fn closed_pipe() -> Stdio {
    let (reader, writer) = std::io::pipe().unwrap();
    drop(reader);
    writer.into()
}

/// The names of the files in `dir`, but `in` and `.ncbirc`.
fn left(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|name| name != "in" && name != ".ncbirc")
        .collect();
    names.sort();
    names
}

#[test]
fn a_text_on_a_closed_log_or_standard_error_ends_the_tool_where_ncbis_ends() {
    let tests = tests_dir();
    let fixtures = tests.to_string_lossy().into_owned();
    let mut failures = Vec::new();
    for (index, case) in CASES.iter().enumerate() {
        let label = format!("{} {} {:?}", case.program, case.name, case.log);
        let home = unique_dir(&format!("losat_logpipe_{index}"));
        let mut run_dir = home.to_path_buf();
        let dir = match case.place {
            Place::Fixtures(name) => tests.join(name),
            Place::Fresh(path) => {
                copy_tree(&tests.join(path), &home.join("in"));
                home.to_path_buf()
            }
            Place::MalformedNcbirc => {
                fs::write(home.join(".ncbirc"), "[BLAST]\nBLASTDB=db\ngarbage\n").unwrap();
                home.to_path_buf()
            }
            Place::CheckCase(name) => {
                run_dir = home.join("run");
                crate::blastdbcheck_fixtures::build_case_dir(&run_dir, name)
            }
        };
        let before = left(&dir);
        let mut run = Command::new(env!("CARGO_BIN_EXE_losat"));
        run.current_dir(&dir)
            .arg(case.program)
            .args(case.words.iter().map(|word| word.replace("{F}", &fixtures)))
            .args(match case.log {
                Log::Stdout => ["-logfile", "/dev/stdout"].as_slice(),
                Log::Stderr => ["-logfile", "/dev/stderr"].as_slice(),
                Log::Cerr => [].as_slice(),
            })
            .env("HOME", home.path())
            .env("NCBI_CONFIG_PATH", home.path())
            .env("LC_ALL", "C")
            .env("TZ", "UTC")
            .env("LOSAT_BLASTX_WORKER_LOG", "/dev/null")
            .env_remove("BLASTDB")
            .env_remove("BLAST_USAGE_REPORT")
            .stdin(Stdio::null());
        for (variable, value) in case.variables {
            run.env(variable, value);
        }
        match case.log {
            Log::Stdout => run.stdout(closed_pipe()).stderr(Stdio::piped()),
            Log::Stderr | Log::Cerr => run.stdout(Stdio::piped()).stderr(closed_pipe()),
        };
        let output = run.output().unwrap();
        let kept = match case.log {
            Log::Stdout => &output.stderr,
            Log::Stderr | Log::Cerr => &output.stdout,
        };
        let status = if case.status == 141 { 0 } else { case.status };
        if output.status.code() != Some(status) {
            failures.push(format!("{label}: exit status {:?}", output.status.code()));
        }
        let same = match case.kept {
            Kept::Bytes(want) => {
                kept.as_slice() == losat::wording::losat_output_bytes(want).as_ref()
            }
            Kept::StartsWith(want) => kept.starts_with(want),
            Kept::Prefix(path, length) => {
                let ncbi = fs::read(tests_dir().join(path)).unwrap();
                kept.as_slice() == losat::wording::losat_output_bytes(&ncbi[..length]).as_ref()
            }
            Kept::Lines(path, lines) => {
                let ncbi = fs::read_to_string(tests_dir().join(path)).unwrap();
                let run = run_dir.to_string_lossy();
                let want: String = ncbi
                    .split_inclusive('\n')
                    .take(lines)
                    .map(|line| line.replace("<RUN>", &run))
                    .collect();
                kept.as_slice() == losat::wording::losat_output_bytes(want.as_bytes()).as_ref()
            }
        };
        if !same {
            failures.push(format!(
                "{label}: kept stream {:?}",
                String::from_utf8_lossy(kept)
            ));
        }
        let made: Vec<String> = left(&dir)
            .into_iter()
            .filter(|name| !before.contains(name))
            .collect();
        let want: Vec<&str> = case.files.iter().map(|(name, _, _)| *name).collect();
        if made != want {
            failures.push(format!("{label}: files {made:?}"));
            continue;
        }
        for &(name, size, hash) in case.files {
            let bytes = fs::read(dir.join(name)).unwrap();
            let by_size = name.ends_with(".pdb") || name.ends_with(".ptf");
            let got = if name.ends_with(".pin") || name.ends_with(".pal") {
                (0, 0)
            } else {
                (
                    bytes.len() as u64,
                    if by_size { 0 } else { fnv1a64(&bytes) },
                )
            };
            if got != (size, hash) {
                failures.push(format!("{label}: {name} {got:?}"));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
