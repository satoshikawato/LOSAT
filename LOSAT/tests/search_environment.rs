//! The environment that NCBI BLAST+ reads around its search and report, in a search of
//! `-subject` sequences by blastp, psiblast, tblastn, tblastx and blastn (BLASTX's cases are in
//! `run_local_blastx.rs`; the `-db` cases are the `db.envshared.*` and `db.envconv.*` rows of
//! `db_search_manifest/`): `BL2SEQ_LEGACY` makes NCBI search each subject on its own and
//! write the bl2seq report; a `CHUNK_SIZE`, `OVERLAP_CHUNK_SIZE`, `BATCH_SIZE` or
//! `PRE_FETCH_SEQS_LIMIT` that `NStr::StringToInt` cannot convert stops NCBI where it reads it
//! with its `CStringException` (exit 255), except `PRE_FETCH_SEQS_LIMIT` for the archive of
//! outfmt 11, which NCBI's apps but psiblast write without the prefetch (line DBID N2e probe:
//! NCBI 2.17.0 exits 0 with the bytes it writes without the variable); and `OLD_FSC` (no Gumbel
//! block) is read as `-searchsp` is read (env parity G1; the `-db` cases are the `db.oldfsc.*`
//! rows).

mod run_local_support;

use std::process::{Command, Output};

use run_local_support::string_to_int_exception;

const BIN: &str = env!("CARGO_BIN_EXE_losat");
const INPUTS: &str = "tests/fasta/outfmt11";

/// The programs with their query and subject files and NCBI's name in LOSAT's messages.
const PROGRAMS: [(&str, &str, &str, &str); 5] = [
    ("blastp", "arch_bp_q.faa", "arch_bp_s.faa", "BLASTP"),
    ("psiblast", "arch_bp_q.faa", "arch_bp_s.faa", "PSIBLAST"),
    ("tblastn", "arch_tbn_q.faa", "arch_tbn_s.fna", "TBLASTN"),
    ("tblastx", "arch_tbx_q.fna", "arch_tbx_s.fna", "TBLASTX"),
    ("blastn", "arch_q.fna", "arch_s.fna", "BLASTN"),
];

/// One run of `losat <program> -query <query> -subject <subject> -outfmt <outfmt>` from
/// `LOSAT/` without the variables that change a search, plus `variable`.
fn run(
    (program, query, subject, _): (&str, &str, &str, &str),
    outfmt: &str,
    variable: Option<(&str, &str)>,
) -> Output {
    run_files(
        program,
        &format!("{INPUTS}/{query}"),
        &format!("{INPUTS}/{subject}"),
        outfmt,
        variable,
        &[],
    )
}

/// `run` with the paths of the query and subject files (relative to `LOSAT/`) and the options
/// `extra`.
fn run_files(
    program: &str,
    query: &str,
    subject: &str,
    outfmt: &str,
    variable: Option<(&str, &str)>,
    extra: &[&str],
) -> Output {
    run_search(
        program,
        &["-query", query, "-subject", subject],
        outfmt,
        variable,
        extra,
    )
}

/// One run of `losat <program> <inputs> -outfmt <outfmt> -num_threads 1 <extra>` from `LOSAT/`
/// without the variables that change a search, plus `variable`.
fn run_search(
    program: &str,
    inputs: &[&str],
    outfmt: &str,
    variable: Option<(&str, &str)>,
    extra: &[&str],
) -> Output {
    let mut command = Command::new(BIN);
    command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .arg(program)
        .args(inputs)
        .args(["-outfmt", outfmt, "-num_threads", "1"])
        .args(extra);
    for (name, _) in std::env::vars() {
        if name.starts_with("LOSAT_X_")
            || matches!(
                name.as_str(),
                "ARCHIVE_FORMAT"
                    | "BATCH_SIZE"
                    | "BL2SEQ_LEGACY"
                    | "BLASTDB"
                    | "CHUNK_SIZE"
                    | "OLD_FSC"
                    | "OVERLAP_CHUNK_SIZE"
                    | "PRE_FETCH_SEQS_LIMIT"
            )
        {
            command.env_remove(name);
        }
    }
    if let Some((name, value)) = variable {
        command.env(name, value);
    }
    command.output().expect("run LOSAT CLI")
}

// NCBI reference: c++/src/app/blast/blast_app_util.cpp:206-210,732-737, blastp_app.cpp:280-285,
// tblastn_app.cpp:294-299, tblastx_app.cpp:197-202, blastn_app.cpp:304-309 and
// psiblast_app.cpp:262-266 (psiblast fetches ahead before it writes the archive).
#[test]
fn bl2seq_legacy_is_refused_and_the_archive_is_written_without_the_prefetch() {
    for case in PROGRAMS {
        let (program, _, _, name) = case;
        // psiblast does not write the archive (outfmt 11 is refused before the environment).
        let archive = if program == "psiblast" { "6" } else { "11" };
        for (value, outfmt) in [("1", "6"), ("", archive)] {
            let output = run(case, outfmt, Some(("BL2SEQ_LEGACY", value)));
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert_eq!(
                output.status.code(),
                Some(1),
                "{program} BL2SEQ_LEGACY {outfmt}: {output:?}"
            );
            assert!(
                output.stdout.is_empty()
                    && stderr.contains("the environment variable BL2SEQ_LEGACY")
                    && stderr.contains(&format!("not supported by LOSAT's {name}")),
                "{program} BL2SEQ_LEGACY {outfmt}: {stderr}"
            );
        }
        // The archive (not psiblast's), and an integer in any format, are written as without
        // the variable.
        let mut written = vec![("5", "6")];
        if program != "psiblast" {
            written.extend([("abc", "11"), ("", "11")]);
        }
        for (value, outfmt) in written {
            let plain = run(case, outfmt, None);
            let output = run(case, outfmt, Some(("PRE_FETCH_SEQS_LIMIT", value)));
            assert!(
                plain.status.success() && !plain.stdout.is_empty(),
                "{program} {outfmt}: {plain:?}"
            );
            assert_eq!(
                output.status.code(),
                Some(0),
                "{program} {value} {outfmt}: {output:?}"
            );
            assert!(
                output.stdout == plain.stdout && output.stderr == plain.stderr,
                "{program} {value} {outfmt}"
            );
        }
    }
}

/// NCBI's error archive of outfmt 11 for that exception (`PrintErrorArchive` after `CATCH_ALL`):
/// the same 523 bytes for every input.
const ABC_ERROR_ARCHIVE: &str = "Blast4-archive ::= {\n  request {\n    body get-request-info {\n      request-id \"Error\"\n    }\n  },\n  results {\n  },\n  messages {\n    {\n      code 2,\n      message \"Error: NCBI C++ Exception: T0 \"\"/home/coremake/release_build/bu\nild/PrepareRelease_Linux64-Alma_JSID_01_170042_130.14.18.128_9008__PrepareRele\nase_Linux64-Alma_1751374548/c++/compilers/unix/../../src/corelib/ncbistr.cpp\n\"\", line 861: Error: (CStringException::eConvert) ncbi::NStr::StringToInt8() -\n Cannot convert string 'abc' to Int8 (m_Pos = 0)\"\n    }\n  }\n}\n";

/// Where NCBI stops a `-subject` search of `program` at a value of `variable` that it cannot
/// convert.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Stop {
    /// Before anything is written.
    Start,
    /// After the outfmt 0 prolog, before the first query's report.
    Prolog,
    /// Nowhere: the program does not read the variable.
    Never,
}

/// The outfmt 0 prolog of a report: its bytes up to the first query's report (psiblast's starts
/// with the line of its round).
fn prolog(report: &[u8]) -> &[u8] {
    let text = String::from_utf8_lossy(report);
    let at = ["\nResults from round ", "\n\n\nQuery= "]
        .iter()
        .filter_map(|start| text.find(start))
        .min()
        .expect("a query's report");
    &report[..at + 1]
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
// ```c++
//     catch (const CException& e) {                                           \
//         LOG_POST(Error << "Error: " << e.what());                           \
//         exit_code = BLAST_UNKNOWN_ERROR;                                    \
//     }                                                                       \
// ```
// NCBI 2.17.0 (the official Linux x64 release) with each variable set to `abc` (env parity probe
// h_probe of 2026-11-12, REPORT.md "Stop point"): `CHUNK_SIZE`, `OVERLAP_CHUNK_SIZE` (the query
// splitter of the first batch, local_blast.cpp:59-61, split_query_aux_priv.cpp:55-57) and
// `PRE_FETCH_SEQS_LIMIT` (the first batch's prefetch, blast_app_util.cpp:732-734) stop after the
// outfmt 0 prolog; `BATCH_SIZE` (blast_input_aux.cpp:86-88) before anything, in blastn after the
// prolog (blastn_app.cpp:258-262), and psiblast does not read it (psiblast_app.cpp:468). The
// archive of outfmt 11 is NCBI's error archive (not psiblast's, whose outfmt 11 LOSAT refuses;
// not for `PRE_FETCH_SEQS_LIMIT`, above). These bytes are NCBI's for these inputs.
#[test]
fn values_that_ncbi_cannot_convert_stop_where_ncbi_reads_them() {
    for case in PROGRAMS {
        let (program, ..) = case;
        for variable in [
            "CHUNK_SIZE",
            "OVERLAP_CHUNK_SIZE",
            "BATCH_SIZE",
            "PRE_FETCH_SEQS_LIMIT",
        ] {
            let stop = match (variable, program) {
                ("BATCH_SIZE", "psiblast") => Stop::Never,
                ("BATCH_SIZE", "blastn") => Stop::Prolog,
                ("BATCH_SIZE", _) => Stop::Start,
                _ => Stop::Prolog,
            };
            for outfmt in ["0", "6", "7", "11"] {
                if outfmt == "11" && (program == "psiblast" || variable == "PRE_FETCH_SEQS_LIMIT") {
                    continue;
                }
                let what = format!("{program} {variable}=abc -outfmt {outfmt}");
                let plain = run(case, outfmt, None);
                assert!(plain.status.success(), "{what}: {plain:?}");
                let output = run(case, outfmt, Some((variable, "abc")));
                if stop == Stop::Never {
                    assert!(
                        output.status.success()
                            && output.stdout == plain.stdout
                            && output.stderr == plain.stderr,
                        "{what}: {output:?}"
                    );
                    continue;
                }
                assert_eq!(output.status.code(), Some(255), "{what}: {output:?}");
                assert_eq!(
                    String::from_utf8_lossy(&output.stderr),
                    string_to_int_exception("abc", 861, 0),
                    "{what}"
                );
                let expected: &[u8] = match (outfmt, stop) {
                    ("11", _) => ABC_ERROR_ARCHIVE.as_bytes(),
                    ("0", Stop::Prolog) => prolog(&plain.stdout),
                    _ => b"",
                };
                assert!(
                    output.stdout == expected,
                    "{what}: {}",
                    String::from_utf8_lossy(&output.stdout)
                );
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:264-275
// ```c++
//     TSeqAlignVector aligns =
//         LocalBlastResults2SeqAlign(hsp_results,
//                                    *qdata,
//                                    *m_SeqInfoSrc,
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:262-264
// ```c++
//                    results = psiblast->Run();
//                    BlastFormatter_PreFetchSequenceData(*results, scope,
//                 		   m_CmdLineArgs->GetFormattingArgs()->GetFormattedOutputChoice());
// ```
// NCBI reads the titles of a batch's reported sequences in its traceback (the Seq-ids of the
// Seq-aligns; psiblast's in `Run`), before the batch's prefetch converts `PRE_FETCH_SEQS_LIMIT`:
// with a database whose header file is missing and a query that has hits, NCBI 2.17.0 stops with
// the header file's error (exit 2, after the outfmt 0 prolog) with `PRE_FETCH_SEQS_LIMIT=abc`
// (or an empty value) as without it (env parity H, probe h_probe/hdr of 2026-11-12: every
// program, outfmt 0, 6 and 7).
#[test]
fn a_header_file_that_ncbi_cannot_read_comes_before_the_prefetch() {
    for (program, db, header, query) in [
        ("blastp", "nophr", "phr", "tests/fasta/psiblast/tk.faa"),
        ("psiblast", "nophr", "phr", "tests/fasta/psiblast/tk.faa"),
        ("blastx", "bnophr", "phr", "tests/db_search/bslice_q.fna"),
        ("blastn", "nnophr", "nhr", "tests/seqdb/nucl.fna"),
        ("tblastx", "nnophr", "nhr", "tests/seqdb/nucl.fna"),
    ] {
        let run = |outfmt: &str, variable: Option<(&str, &str)>| {
            run_search(
                program,
                &["-db", &format!("tests/db_search/{db}"), "-query", query],
                outfmt,
                variable,
                &[],
            )
        };
        let error = format!(
            "BLAST Database error: Error: File ({}/tests/db_search/{db}.{header}) not found.\n",
            env!("CARGO_MANIFEST_DIR")
        );
        for outfmt in ["0", "6", "7"] {
            let what = format!("{program} -db {db} -outfmt {outfmt}");
            let plain = run(outfmt, None);
            assert_eq!(plain.status.code(), Some(2), "{what}: {plain:?}");
            assert_eq!(String::from_utf8_lossy(&plain.stderr), error, "{what}");
            for value in ["abc", ""] {
                let output = run(outfmt, Some(("PRE_FETCH_SEQS_LIMIT", value)));
                assert!(
                    output.status.code() == Some(2)
                        && output.stdout == plain.stdout
                        && output.stderr == plain.stderr,
                    "{what} PRE_FETCH_SEQS_LIMIT={value:?}: {output:?}"
                );
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:159-160
// ```c++
//     int from = NStr::StringToInt(tokens.front());
//     int to = NStr::StringToInt(tokens.back());
// ```
// NCBI 2.17.0 blastn and blastp `-query_loc` and `-subject_loc` (env parity probe h_probe of
// 2026-11-12, cases `LOC.*`): a part that `NStr::StringToInt` cannot convert ends the program
// with its `CStringException` before anything is written (exit 255), the start before the stop.
#[test]
fn a_location_part_that_ncbi_cannot_convert_is_its_exception() {
    for (program, query, subject) in [
        ("blastn", "arch_q.fna", "arch_s.fna"),
        ("blastp", "arch_bp_q.faa", "arch_bp_s.faa"),
    ] {
        for option in ["-query_loc", "-subject_loc"] {
            for (value, part, line, pos) in [
                ("a-10", "a", 861, 0),
                ("1-b", "b", 861, 0),
                ("1-1.5", "1.5", 870, 1),
            ] {
                for outfmt in ["0", "6"] {
                    let output = run_files(
                        program,
                        &format!("{INPUTS}/{query}"),
                        &format!("{INPUTS}/{subject}"),
                        outfmt,
                        None,
                        &[option, value],
                    );
                    let what = format!("{program} {option} {value} -outfmt {outfmt}");
                    assert_eq!(output.status.code(), Some(255), "{what}: {output:?}");
                    assert!(output.stdout.is_empty(), "{what}: {output:?}");
                    assert_eq!(
                        String::from_utf8_lossy(&output.stderr),
                        string_to_int_exception(part, line, pos),
                        "{what}"
                    );
                }
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:922-923
// ```c
//     use_old_fsc = getenv("OLD_FSC");
//     if (!use_old_fsc) sbp->gbp = s_BlastGumbelBlkNew();
// ```
// NCBI sets OLD_FSC for any -searchsp (blast_args.cpp:303-307), so BLASTP, BLASTX, TBLASTN and
// PSIBLAST with OLD_FSC (any value, the empty one too) write the frozen bytes of `-searchsp 0`
// (rows sp.<program>.zero.0 and ps.subj.sp.zero.tk.0 of outfmt0_manifest.tsv; NCBI 2.17.0 writes
// the same bytes with OLD_FSC=1, and psiblast with OLD_FSC=true, probe `ofsctrue.tk.subj.0` of line
// PS, item PS5c). RPSBLAST, which takes no -searchsp, refuses it (row db.oldfsc.rpsblast.rej).
#[test]
fn old_fsc_is_the_statistics_of_searchsp_0() {
    for (program, query, subject, frozen) in [
        (
            "blastp",
            "tests/db_search/bpqb_q.faa",
            "tests/db_search/bslice.faa",
            "sp.blastp.zero.0",
        ),
        (
            "blastx",
            "tests/db_search/bslice_q.fna",
            "tests/db_search/bslice.faa",
            "sp.blastx.zero.0",
        ),
        (
            "tblastn",
            "tests/db_search/tslice_q.faa",
            "tests/db_search/tslice.fna",
            "sp.tblastn.zero.0",
        ),
        (
            "psiblast",
            "tests/fasta/psiblast/tk.faa",
            "tests/fasta/PajaWSV.faa",
            "ps.subj.sp.zero.tk.0",
        ),
    ] {
        let expected = std::fs::read(format!(
            "{}/tests/fixtures/outfmt0/{frozen}.out",
            env!("CARGO_MANIFEST_DIR")
        ))
        .expect("frozen -searchsp 0 report");
        for value in ["1", ""] {
            let output = run_files(program, query, subject, "0", Some(("OLD_FSC", value)), &[]);
            assert_eq!(
                output.status.code(),
                Some(0),
                "{program} OLD_FSC={value:?}: {output:?}"
            );
            assert!(
                output.stdout == expected && output.stderr.is_empty(),
                "{program} OLD_FSC={value:?}: {}",
                String::from_utf8_lossy(&output.stderr)
            );
        }
    }
}
