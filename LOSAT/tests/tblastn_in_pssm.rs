//! tblastn `-in_pssm` (PSI-TBLASTN; line PS, item PS7, DECISIONS 35): what a fixture row cannot
//! hold. The searches and NCBI's errors are fixture rows (`db_search_manifest`, prefix
//! `db.tbn.inpssm.`, frozen from NCBI BLAST+ 2.17.0).

use std::path::Path;
use std::process::{Command, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

/// `losat tblastn -db tests/db_search/tslice -outfmt 6 -num_threads 1 <args>` from `LOSAT/`, without
/// the variables that change a search.
fn tblastn(args: &[&str], stdin: Stdio) -> std::process::Output {
    tblastn_on(&["-db", "tests/db_search/tslice"], "1", args, stdin)
}

/// `losat tblastn <subjects> -outfmt 6 -num_threads <threads> <args>` from `LOSAT/`, without the
/// variables that change a search.
fn tblastn_on(
    subjects: &[&str],
    threads: &str,
    args: &[&str],
    stdin: Stdio,
) -> std::process::Output {
    let mut command = Command::new(BIN);
    command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .arg("tblastn")
        .args(subjects)
        .args(["-outfmt", "6", "-num_threads", threads])
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

// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:159-162
// ```c++
// 	if ((mt_mode == CMTArgs::eSplitByQueries) && (num_threads > 1)){
// 		m_UsageReport.AddParam(CBlastUsageReport::eMTMode, CMTArgs::eSplitByQueries);
// 		return x_RunMTBySplitQuery();
// 	}
// ```
// NCBI 2.17.0 (step-1 cases `ts1.mt1.n2.6`, `ts1.mt1.n1.6`, `ts1.subj.mt1.n2.6` of line PS PS7): a PSSM
// search split by queries aborts with `-db` (LOSAT's refusal, DECISIONS 35 Q4), and runs as one
// thread with one thread or with `-subject` (whose thread warning LOSAT does not write,
// PD-LOSAT-CLI-NONSEARCH-DIFFERENCES). A fixture row cannot hold it: its runner sets `-num_threads`.
#[test]
fn mt_mode_1_with_threads_and_a_database_is_refused() {
    let pssm =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/db_search/db.ps.asn.rnd.6.asn");
    let pssm = pssm.to_str().unwrap();
    let db = ["-db", "tests/db_search/tslice"];
    let refused = tblastn_on(
        &db,
        "2",
        &["-in_pssm", pssm, "-mt_mode", "1"],
        Stdio::null(),
    );
    assert_eq!(refused.status.code(), Some(1), "{refused:?}");
    assert!(refused.stdout.is_empty());
    assert!(String::from_utf8_lossy(&refused.stderr)
        .contains("-mt_mode 1 with -num_threads above 1 and -db"));
    let args = ["-in_pssm", pssm, "-evalue", "1000"];
    let one = tblastn_on(&db, "1", &args, Stdio::null());
    let mut mt1 = args.to_vec();
    mt1.extend(["-mt_mode", "1"]);
    let one_mt1 = tblastn_on(&db, "1", &mt1, Stdio::null());
    assert_eq!(one.status.code(), Some(0));
    assert_eq!(
        (&one_mt1.stdout, one_mt1.status.code()),
        (&one.stdout, Some(0))
    );
    let subject = ["-subject", "tests/db_search/tslice.fna"];
    let subject_one = tblastn_on(&subject, "1", &args, Stdio::null());
    let subject_mt1 = tblastn_on(&subject, "2", &mt1, Stdio::null());
    assert_eq!(subject_one.status.code(), Some(0));
    assert_eq!(
        (&subject_mt1.stdout, subject_mt1.status.code()),
        (&subject_one.stdout, Some(0))
    );
}
