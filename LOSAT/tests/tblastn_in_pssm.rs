//! tblastn `-in_pssm` (PSI-TBLASTN; line PS, item PS7, DECISIONS 35): what a fixture row cannot
//! hold. The searches and NCBI's errors are fixture rows (`db_search_manifest`, prefix
//! `db.tbn.inpssm.`, frozen from NCBI BLAST+ 2.17.0).

use std::path::Path;
use std::process::{Command, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// `losat tblastn -db tests/db_search/tslice -outfmt 6 -num_threads 1 <args>` from `LOSAT/`, without
/// the variables that change a search.
fn tblastn(args: &[&str], stdin: Stdio) -> std::process::Output {
    let mut command = Command::new(BIN);
    command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .arg("tblastn")
        .args([
            "-db",
            "tests/db_search/tslice",
            "-outfmt",
            "6",
            "-num_threads",
            "1",
        ])
        .args(args)
        .stdin(stdin);
    for (name, _) in std::env::vars() {
        if name.starts_with("LOSAT_X_")
            || name == "LOSAT_TIMING"
            || name.starts_with("NCBI")
            || name.starts_with("BLAST")
            || matches!(
                name.as_str(),
                "BATCH_SIZE" | "BL2SEQ_LEGACY" | "CHUNK_SIZE" | "PRE_FETCH_SEQS_LIMIT" | "OLD_FSC"
            )
        {
            command.env_remove(name);
        }
    }
    command.output().expect("run LOSAT CLI")
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:717-721
// ```c++
//     if (AsString() == "-") {
// #if defined(NCBI_OS_MSWIN)
//         NcbiSys_setmode(NcbiSys_fileno(stdin), (mode & IOS_BASE::binary) ? O_BINARY : O_TEXT);
// #endif
//         m_Ios  = &cin;
// ```
// NCBI 2.17.0 (step-1 case `pf.stdin.6` of line PS PS7: the bytes of the run that names the file):
// `-in_pssm -` reads the PSSM from the standard input, and the default `-query -` is not read.
#[test]
fn in_pssm_dash_reads_the_standard_input() {
    let pssm =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/db_search/db.ps.asn.rnd.6.asn");
    let named = tblastn(
        &["-in_pssm", pssm.to_str().unwrap(), "-evalue", "1000"],
        Stdio::null(),
    );
    let piped = tblastn(
        &["-in_pssm", "-", "-evalue", "1000"],
        std::fs::File::open(&pssm).unwrap().into(),
    );
    assert_eq!(named.status.code(), Some(0), "{named:?}");
    assert!(!named.stdout.is_empty());
    assert_eq!(
        (piped.status.code(), &piped.stdout, &piped.stderr),
        (named.status.code(), &named.stdout, &named.stderr)
    );
    let empty = tblastn(&["-in_pssm", "-"], Stdio::null());
    assert_eq!(empty.status.code(), Some(1));
    assert_eq!(
        empty.stderr,
        b"BLAST query/options error: Unsupported format for PSSM\nPlease refer to the BLAST+ user manual.\n"
    );
}
