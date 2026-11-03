//! The environment that NCBI BLAST+ reads around its search and report, in a search of
//! `-subject` sequences by blastp, psiblast, tblastn, tblastx and blastn (BLASTX's cases are in
//! `run_local_blastx.rs`; the `-db` cases are the `db.envshared.*` rows of
//! `db_search_manifest.tsv`): `BL2SEQ_LEGACY` makes NCBI search each subject on its own and
//! write the bl2seq report, a `PRE_FETCH_SEQS_LIMIT` that `NStr::StringToInt` cannot convert
//! (an empty one too) stops NCBI before each batch's report, except for the archive of outfmt 11,
//! which NCBI's apps but psiblast write without the prefetch (line DBID N2e probe: NCBI 2.17.0
//! exits 0 with the bytes it writes without the variable).

use std::process::{Command, Output};

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
    let mut command = Command::new(BIN);
    command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .arg(program)
        .args(["-query", &format!("{INPUTS}/{query}")])
        .args(["-subject", &format!("{INPUTS}/{subject}")])
        .args(["-outfmt", outfmt, "-num_threads", "1"]);
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
fn environment_that_ncbi_reads_around_the_report_is_refused_where_it_reads_it() {
    for case in PROGRAMS {
        let (program, _, _, name) = case;
        // psiblast does not write the archive (outfmt 11 is refused before the environment).
        let archive = if program == "psiblast" { "6" } else { "11" };
        for (variable, value, outfmt, text) in [
            (
                "BL2SEQ_LEGACY",
                "1",
                "6",
                "the environment variable BL2SEQ_LEGACY",
            ),
            (
                "BL2SEQ_LEGACY",
                "",
                archive,
                "the environment variable BL2SEQ_LEGACY",
            ),
            ("PRE_FETCH_SEQS_LIMIT", "abc", "6", "has the value \"abc\""),
            ("PRE_FETCH_SEQS_LIMIT", "", "0", "has the value \"\""),
        ] {
            let output = run(case, outfmt, Some((variable, value)));
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert_eq!(
                output.status.code(),
                Some(1),
                "{program} {variable} {outfmt}: {output:?}"
            );
            assert!(
                output.stdout.is_empty()
                    && stderr.contains(text)
                    && stderr.contains(&format!("not supported by LOSAT's {name}")),
                "{program} {variable} {outfmt}: {stderr}"
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
