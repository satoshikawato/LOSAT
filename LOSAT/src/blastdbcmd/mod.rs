//! `losat blastdbcmd`: NCBI BLAST+ 2.17.0 blastdbcmd, of which LOSAT serves `-info` (the
//! information of a database: `-db`, `-dbtype`, `-exact_length`, `-out`, `-logfile`), the
//! retrieval of sequences (`retrieve`: `-entry all`, `-entry`, `-entry_batch`, `-ipg`,
//! `-ipg_batch`, `-taxids`, `-taxidlist`, every `-outfmt` but the ASN.1 texts, `-range`,
//! `-strand`, `-line_length`, `-ctrl_a`, `-target_only`, `-get_dups`, `-long_seqids`), the
//! taxids of a database (`taxonomy`: `-tax_info`), its metadata (`metadata`: `-metadata`,
//! `-metadata_output_prefix`), the databases of a directory (`list`: `-list`,
//! `-list_outfmt`, `-recursive`, `-remove_redundant_dbs`), `-show_blastdb_search_path`, `-h`,
//! `-help` and `-version`. Its command line is NCBI's (`args`, NCBI's argument reader), not
//! clap's: `main` hands it the words after `blastdbcmd`. The masks of `makeblastdb -mask_data`
//! are read (`-info`, `-mask_sequence_with`, `%m`). A `-db` of several names is read as NCBI's
//! reader reads it, as one database (`SeqDb::open_in`).

mod args;
mod fasta;
mod format;
mod letters;
mod list;
mod metadata;
mod range;
mod retrieve;
mod taxonomy;

use std::ffi::OsString;
use std::io::{Read, Write};

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::ncbi_app::args::Args;
use crate::ncbi_app::files::{Cin, Ostream};
use crate::ncbi_app::{Diag, Tool};
use crate::ncbi_str::{int_with_commas, with_commas};
use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType, TaxidList};

/// The name in LOSAT's refusals.
const PROGRAM: &str = "BLASTDBCMD";

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:907-909
// ```c++
//     arg_desc->SetUsageContext(GetArguments().GetProgramBasename(),
//                   "BLAST database client, version " + CBlastVersion().Print());
// ```
/// The text of `blastdbcmd -help` (NCBI 2.17.0's, `CArgDescriptions::PrintUsage` of the
/// descriptions in `args`). Its part before `OPTIONAL ARGUMENTS` is the usage that `-h` and the
/// argument errors print.
const HELP: &str = include_str!("help.txt");

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1175-1177
// ```c++
//     SetDiagPostLevel(eDiag_Warning);
//     SetDiagPostPrefix("blastdbcmd");
// ```
/// blastdbcmd as NCBI's application runs it (`ncbi_app`): its `-version` is that of BLAST+
/// (`CBlastDBCmdApp` sets `CBlastVersion`), its diagnostics carry the prefix `[blastdbcmd]`.
static TOOL: Tool = Tool {
    program: "blastdbcmd",
    refusal_name: PROGRAM,
    version: crate::NCBI_BLAST_VERSION,
    help: HELP,
    descs: args::DESCS,
    dependencies: args::dependencies,
    post_prefix: Some("blastdbcmd"),
};

/// Runs `blastdbcmd` with `words` (the command line after `losat blastdbcmd`) on the process's
/// standard streams; the exit status.
pub fn run(words: &[OsString]) -> i32 {
    crate::ncbi_app::run(&TOOL, words, run_body)
}

/// The run of one command line with `input` as its standard input (a batch file `-`), its
/// output on `out` and its errors on `err`; the exit status (NCBI's application around
/// blastdbcmd's `Run`: `ncbi_app::run_with`).
pub(crate) fn run_with(
    words: &[OsString],
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_with(&TOOL, words, input, out, err, run_body)
}

/// blastdbcmd's `Run` (`CBlastDBCmdApp::Run`) on the arguments that NCBI's reader read.
///
/// LOSAT's refusals come where NCBI would start the work they stand for, and before anything is
/// written: the program's other uses after the arguments are read, the database where it is
/// opened (after `-out` is opened, as NCBI's `Run` opens it first).
fn run_body(
    args: &Args,
    settings: &ApplicationSettings,
    stdin: &mut Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32 {
    let registry = match &settings.blastdb {
        Ok(value) => value.as_deref(),
        Err(reason) => return diag.refuse(&refusal(reason)),
    };
    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1179-1180
    // ```c++
    //     try {
    //         CNcbiOstream& out = args["out"].AsOutputFile();
    // ```
    let mut out = match args.output_file("out", out) {
        Ok(out) => out,
        Err(error) => {
            // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-180
            // ```c++
            //     catch (const CArgException& e) {                                        \
            //         LOG_POST(Error << "Command line argument error: " << e.GetMsg());   \
            //         exit_code = BLAST_INPUT_ERROR;                                      \
            //     }                                                                       \
            // ```
            diag.post(&format!("Command line argument error: {}", error.message()));
            return 1;
        }
    };
    // NCBI's `ostream` of `-out` for every mode: a failed write drops the later ones, and a
    // closed pipe ends the program where SIGPIPE ends NCBI's (`files::Ostream`).
    let mut out = Ostream::new(&mut out);
    let path = match SearchPath::from_environment(registry) {
        Ok(path) => path,
        Err(error) => return database_error(diag, error),
    };
    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1181-1184
    // ```c++
    //         if (args["show_blastdb_search_path"]) {
    //             out << CSeqDB::GenerateSearchPath() << NcbiEndl;
    //             return status;
    // ```
    if args.given("show_blastdb_search_path") {
        let _ = writeln!(out, "{}", path.as_str());
        let _ = out.flush();
        return 0;
    }
    if args.given("list") {
        return list::list(args, &path, diag, &mut out);
    }
    let name = args.value("db").unwrap_or("nr");
    let dbtype = args.value("dbtype").unwrap_or("guess");
    let selected = ["info", "metadata", "tax_info"]
        .into_iter()
        .any(|selector| args.given(selector));
    if !selected && (args.given("taxids") || args.given("taxidlist")) {
        // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:358-369
        // ```c++
        //     CSeqDB::ESeqType seqtype = ParseMoleculeTypeString(args[kArgDbType].AsString());
        //     m_DbIsProtein = static_cast<bool>(seqtype == CSeqDB::eProtein);
        //     m_TargetOnly = args["target_only"];
        //    	if(m_TargetOnly) {
        //     	CRef<CSeqDBGiList> taxid_list(new CSeqDBGiList());
        //     	taxid_list->AddTaxIds(m_TaxIdList);
        //    		m_BlastDb.Reset(new CSeqDBExpert(args[kArgDb].AsString(), seqtype, taxid_list.GetPointer()));
        //     }
        //    	else {
        //    		m_BlastDb.Reset(new CSeqDBExpert(args[kArgDb].AsString(), seqtype));
        //    	}
        // ```
        // With `-target_only` the list's errors are those of the open, inside NCBI's guess
        // of the type: a protein database the list fails on is then looked for as a
        // nucleotide one (probe J111 of line BDC: `No alias or index file found for
        // nucleotide database [db/p4]`).
        let tax_ids = match taxonomy::tax_id_list(args, &path, diag, stdin) {
            Ok(tax_ids) => tax_ids,
            Err(taxonomy::Ended(status)) => return status,
        };
        let target_only = args.given("target_only");
        let list = TaxidList::Positive(tax_ids.clone());
        let restrict = |db: &mut SeqDb| -> Result<(), SeqDbError> {
            if target_only {
                db.restrict_taxids(&list)?;
            }
            Ok(())
        };
        let (db, _) = match open_database(name, dbtype, &path, &restrict) {
            Ok(opened) => opened,
            Err(error) => return database_error(diag, error),
        };
        let dbtype_protein = dbtype == "prot";
        return retrieve::process(
            &db,
            args,
            &path,
            diag,
            stdin,
            &mut out,
            Some(&tax_ids),
            dbtype_protein,
        );
    }
    let (db, volumes) = match open_database(name, dbtype, &path, &|_| Ok(())) {
        Ok(opened) => opened,
        Err(error) => return database_error(diag, error),
    };
    if !args.given("info") {
        if args.given("metadata") {
            return metadata::metadata(&db, name, args, &path, diag, &mut out);
        }
        if args.given("tax_info") {
            return taxonomy::tax_info(&db, args, &path, diag, &mut out);
        }
        let protein = db.seqtype() == SeqType::Protein;
        return retrieve::process(&db, args, &path, diag, stdin, &mut out, None, protein);
    }
    // NCBI prints the date bytes of the index file; the reader keeps them as text
    // (`SeqDb::date`), with U+FFFD for bytes that are not UTF-8.
    if db.date().contains('\u{fffd}') {
        return diag.refuse(&refusal(&format!(
            "{name}: the date of the database is not UTF-8 text"
        )));
    }
    let mut text = Vec::new();
    write_info_head(&mut text, &db);
    let masks = db
        .mask_algorithms()
        .and_then(|algorithms| algorithms.descriptions());
    let tail = match &masks {
        Ok(masks) => {
            text.extend_from_slice(masks.as_bytes());
            write_info_volumes(&mut text, &volumes);
            None
        }
        Err(error) => Some(error),
    };
    // NCBI's stream ignores a failed write (no exception is enabled on it), so does this.
    // The lines before a bad mask column are written when the error is NCBI's own.
    if tail.is_none_or(|error| error.ncbi_message().is_some()) {
        let _ = out.write_all(&text);
        let _ = out.flush();
        // `cerr`'s tie flushed the output before the error: a closed pipe ends the program.
        if out.closed() {
            return 0;
        }
    }
    match masks {
        Ok(_) => 0,
        Err(error) => database_error(diag, error),
    }
}

/// LOSAT's refusal of a database (`algorithm::common::database::refusal`).
fn refusal(reason: &str) -> String {
    crate::algorithm::common::database::refusal(reason, PROGRAM).to_string()
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:110-153
// ```c++
// static CSeqDBImpl *
// s_SeqDBInit(const string       & dbname,
//             char                 prot_nucl,
//             ...
// {
//     CSeqDBImpl * impl = 0;
//
//     if (prot_nucl == '-') {
//         try {
//             prot_nucl = 'p';
//             impl = new CSeqDBImpl(dbname,
//                                   prot_nucl,
//                                   ...
//         }
//         catch(CSeqDBException &) {
//             prot_nucl = 'n';
//         }
//     }
//
//     if (! impl) {
//         impl = new CSeqDBImpl(dbname,
//                               prot_nucl,
// ```
// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:524-531
// ```c++
//     CSeqDB::ESeqType seqtype = ParseMoleculeTypeString(args[kArgDbType].AsString());
//     m_BlastDb.Reset(new CSeqDBExpert(args[kArgDb].AsString(), seqtype));
// ```
/// Opens the database `name` (the `-db` value as given, one name or several: blastdbcmd does not
/// trim it) of the type `dbtype` (`prot`, `nucl`, or `guess`: protein first, then nucleotide),
/// with the list of its volumes for `-info`, and `prepare`s it (the taxid list of `-target_only`, which NCBI's
/// constructor applies). With `guess`, an error that NCBI raises too (`ncbi_message`) moves
/// on to the nucleotide database, and the nucleotide error is the one reported; an error of
/// LOSAT's is reported at once, since NCBI may open that protein database.
fn open_database(
    name: &str,
    dbtype: &str,
    path: &SearchPath,
    prepare: &dyn Fn(&mut SeqDb) -> Result<(), SeqDbError>,
) -> Result<(SeqDb, Vec<String>), SeqDbError> {
    let types: &[SeqType] = match dbtype {
        "prot" => &[SeqType::Protein],
        "nucl" => &[SeqType::Nucleotide],
        _ => &[SeqType::Protein, SeqType::Nucleotide],
    };
    let mut types = types.iter().peekable();
    loop {
        let Some(&seqtype) = types.next() else {
            unreachable!("at least one sequence type is tried");
        };
        // The alias files' counts fail where NCBI's constructor reaches them, which the list of
        // `-target_only` moves (`SeqDb::check_alias_values`).
        let opened = SeqDb::open_search_in(name, seqtype, path).and_then(|mut db| {
            let volumes = SeqDb::info_volumes(name, seqtype, path)?;
            prepare(&mut db)?;
            db.check_alias_values()?;
            Ok((db, volumes))
        });
        match opened {
            Err(error) if types.peek().is_some() && error.ncbi_message().is_some() => {}
            opened => return opened,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203
// ```c++
//     catch (const CSeqDBException& e) {                                      \
//         ...
//         else {                                                              \
//         	LOG_POST(Error << "BLAST Database error: " << e.GetMsg());      \
//         	exit_code = BLAST_DATABASE_ERROR;                               \
//         }                                                                   \
//     }                                                                       \
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:247-250
// ```c++
//     catch (const CException& e) {                                           \
//         LOG_POST(Error << "Error: " << e.what());                           \
//         exit_code = BLAST_UNKNOWN_ERROR;                                    \
//     }                                                                       \
// ```
/// A reader error: NCBI's own (`SeqDbError::ncbi_message`) as `BLAST Database error: <text>`
/// with exit status 2, another exception of NCBI's (`SeqDbError::ncbi_exception`: an alias
/// file's value that does not convert) as `Error: <what()>` with exit status 255, else LOSAT's
/// refusal.
fn database_error(diag: &mut Diag, error: SeqDbError) -> i32 {
    if let Some(exception) = error.ncbi_exception() {
        diag.post(&format!("Error: {}", exception.what()));
        return 255;
    }
    match error.ncbi_message() {
        Some(text) => {
            diag.post(&format!("BLAST Database error: {text}"));
            2
        }
        None => diag.refuse(&refusal(&error.to_string())),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:533-573
// ```c++
//     static const NStr::TNumToStringFlags kFlags = NStr::fWithCommas;
//     const string kLetters = m_DbIsProtein ? "residues" : "bases";
//     const string kVersion = (m_BlastDb->GetBlastDbVersion() == EBlastDbVersion::eBDB_Version5) ? "5":"4";
//     ...
//     // Print basic database information
//     out << "Database: " << m_BlastDb->GetTitle() << endl
//         << "\t" << NStr::IntToString(m_BlastDb->GetNumSeqs(), kFlags)
//         << " sequences; ";
//         if(args["exact_length"])
//         	out << NStr::UInt8ToString(m_BlastDb->GetExactTotalLength(), kFlags);
//         else
//         	out << NStr::UInt8ToString(m_BlastDb->GetTotalLength(), kFlags);
//     out << " total " << kLetters << endl << endl
//         << "Date: " << m_BlastDb->GetDate()
//         << "\tLongest sequence: "
//         << NStr::IntToString(m_BlastDb->GetMaxLength(), kFlags) << " "
//         << kLetters << endl << endl;
//
//     out << "BLASTDB Version: " << kVersion << endl;
//     ...
//     out << m_BlastDb->GetAvailableMaskAlgorithmDescriptions();
//     ...
//     out << endl << "Volumes:" << endl;
//     ITERATE(vector<string>, file_name, volumes) {
//         out << "\t" << *file_name << endl;
//     }
// ```
/// The text of `-info` for `db` up to its filtering algorithms.
///
/// The counts are `GetNumSeqs` and `GetTotalLength` of the database opened without a list:
/// the walks of the alias tree (`SeqDb::printed_num_seqs`). `-exact_length` changes nothing
/// here: NCBI's exact total differs from the total only when its totals are counted over a
/// filtered list of sequences (an alias file with `GILIST`, `SEQIDLIST`, `OIDLIST`, ...,
/// seqdbimpl.cpp:104-131), which LOSAT's reader refuses; the OID ranges (`FIRST_OID`,
/// `LAST_OID`) do not make NCBI scan.
fn write_info_head(out: &mut Vec<u8>, db: &SeqDb) {
    let letters = match db.seqtype() {
        SeqType::Protein => "residues",
        SeqType::Nucleotide => "bases",
    };
    out.extend_from_slice(b"Database: ");
    out.extend_from_slice(db.title());
    out.push(b'\n');
    let _ = write!(
        out,
        "\t{} sequences; {} total {letters}\n\nDate: {}\tLongest sequence: {} {letters}\n\n\
         BLASTDB Version: {}\n",
        int_with_commas(db.printed_num_seqs()),
        with_commas(db.printed_total_length()),
        db.date(),
        // `GetMaxLength` is NCBI's `int` of the index files' longest length.
        int_with_commas(db.max_length() as i32),
        if db.is_version5() { "5" } else { "4" },
    );
}

/// The text of `-info` after the filtering algorithms: the `Volumes:` are `volumes`.
fn write_info_volumes(out: &mut Vec<u8>, volumes: &[String]) {
    out.extend_from_slice(b"\nVolumes:\n");
    for volume in volumes {
        let _ = writeln!(out, "\t{volume}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn commas_group_three_digits() {
        assert_eq!(with_commas(0), "0");
        assert_eq!(with_commas(999), "999");
        assert_eq!(with_commas(1_000), "1,000");
        assert_eq!(with_commas(186_063_409), "186,063,409");
    }

    // help.txt is NCBI's frozen `-help` text and is not edited or split: this keeps the version
    // that it names and the version that LOSAT writes elsewhere from drifting apart.
    #[test]
    fn the_frozen_help_names_the_version_that_losat_writes() {
        assert!(HELP.contains(&format!(
            "   BLAST database client, version {}\n",
            crate::NCBI_BLAST_VERSION
        )));
    }

    #[test]
    fn the_short_usage_is_the_head_of_the_help() {
        let usage = crate::ncbi_app::short_usage(HELP);
        assert_eq!(usage.len(), 780);
        assert!(usage.starts_with("USAGE\n  blastdbcmd [-h] [-help] [-db dbname]"));
        assert!(usage.ends_with(
            "   BLAST database client, version 2.17.0+\n\nUse '-help' to print detailed \
             descriptions of command line arguments\n"
        ));
        assert_eq!(HELP.len(), 9450);
    }

    /// An output whose every write fails (a full disk).
    struct FullDisk;

    impl Write for FullDisk {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("no space left"))
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Err(std::io::Error::other("no space left"))
        }
    }

    /// A pipe whose reader has closed.
    struct ClosedPipe;

    impl Write for ClosedPipe {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::from(std::io::ErrorKind::BrokenPipe))
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    /// The exit status and the diagnostics of `words` run with `input` as the standard input
    /// and `out` as the standard output.
    fn run_on(words: &[&str], input: &str, out: &mut dyn Write) -> (i32, String) {
        let words: Vec<OsString> = words.iter().map(OsString::from).collect();
        let mut err = Vec::new();
        let status = run_with(&words, &mut input.as_bytes(), out, &mut err);
        (status, String::from_utf8(err).unwrap())
    }

    fn fixture_db(name: &str) -> String {
        format!("{}/tests/{name}", env!("CARGO_MANIFEST_DIR"))
    }

    // NCBI 2.17.0 with `-out /dev/full` (review 2 of line RG item K1, F3): a failed write leaves
    // the retrieval going on, so an id that is not found after it is reported as without the
    // failure (`Skipped NOSUCHID`, exit status 1). The batch is the standard input (`-`).
    #[test]
    fn a_failed_write_leaves_the_later_entries_of_a_batch_and_their_errors() {
        let db = fixture_db("blastdbcmd/db/p5");
        let words = ["-db", db.as_str(), "-entry_batch", "-"];
        let batch = "XP_000001.1\nnosuchid\nP99999.2\n";
        let mut written = Vec::new();
        let expected = run_on(&words, batch, &mut written);
        assert_eq!(
            expected,
            (1, "Error: [blastdbcmd] Skipped NOSUCHID\n".to_string())
        );
        assert!(written.starts_with(b">"));
        assert_eq!(run_on(&words, batch, &mut FullDisk), expected);
        // A closed pipe ends NCBI's process by SIGPIPE (exit status 141) outside Windows before
        // the `Skipped` line, as `cerr`'s tie flushes the output first: the retrieval ends with
        // exit status 0 and no message (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 10). On
        // Windows it is a failed write like the others.
        let closed = if cfg!(windows) {
            expected
        } else {
            (0, String::new())
        };
        assert_eq!(run_on(&words, batch, &mut ClosedPipe), closed);
    }

    /// Standard output and standard error into one file (`2>&1`).
    #[derive(Clone)]
    struct Merged(std::rc::Rc<std::cell::RefCell<Vec<u8>>>);

    impl Write for Merged {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            self.0.borrow_mut().extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    // `cerr` is tied to `cout`: NCBI flushes the output before each diagnostic, so with `2>&1`
    // the sequences of the entries before a missing id come before its `Skipped` line, and an
    // unfinished line before an `-outfmt` error (review of line RG item K1b, R2-F1).
    #[test]
    fn diagnostics_come_after_the_output_written_before_them() {
        let db = fixture_db("blastdbcmd/db/p5");
        for (words, batch, error) in [
            (
                vec!["-db", db.as_str(), "-entry_batch", "-", "-outfmt", "%a"],
                "XP_000001.1\nnosuchid\n",
                "Error: [blastdbcmd] Skipped NOSUCHID\n",
            ),
            (
                vec![
                    "-db",
                    db.as_str(),
                    "-entry",
                    "XP_000001.1",
                    "-outfmt",
                    "%a %Z",
                ],
                "",
                "Error: [blastdbcmd] Unrecognized format specification: '%Z'\n",
            ),
        ] {
            let mut alone = Vec::new();
            let (status, err) = run_on(&words, batch, &mut alone);
            assert_eq!(err, error, "{words:?}");
            let merged = Merged(Default::default());
            let words: Vec<OsString> = words.iter().map(OsString::from).collect();
            let status_merged = run_with(
                &words,
                &mut batch.as_bytes(),
                &mut merged.clone(),
                &mut merged.clone(),
            );
            assert_eq!(status_merged, status);
            let mut expected = alone;
            expected.extend_from_slice(error.as_bytes());
            assert_eq!(*merged.0.borrow(), expected, "{words:?}");
        }
        // A bad mask algorithm id of a batch line, between the outputs of its neighbours
        // (NCBI 2.17.0 `-ipg_batch - -outfmt "%a %g" 2>&1`, review of K1b, R3-F1).
        let db = fixture_db("blastdbcmd/db/m/pmg");
        let words: Vec<OsString> = ["-db", db.as_str(), "-ipg_batch", "-", "-outfmt", "%a %g"]
            .iter()
            .map(OsString::from)
            .collect();
        let merged = Merged(Default::default());
        let status = run_with(
            &words,
            &mut &b"1001 21\n1002 99\n1003 21\n"[..],
            &mut merged.clone(),
            &mut merged.clone(),
        );
        assert_eq!(status, 1);
        assert_eq!(
            String::from_utf8(merged.0.borrow().clone()).unwrap(),
            "XP_1001.1 1001\nError: [blastdbcmd] Invalid filtering algorithm ID: 99\nError: \
             [blastdbcmd] Skipped IPG : 1002\nXP_1003.1 1003\n"
        );
    }

    // NCBI 2.17.0 with `-out /dev/full` (F2): the serial stream of `-metadata` reports the
    // failed write (five lines, exit status 6), which LOSAT refuses.
    #[test]
    fn a_failed_write_of_the_metadata_is_refused() {
        let db = fixture_db("dustmasker/db/basic");
        let words = ["-db", db.as_str(), "-metadata"];
        let mut written = Vec::new();
        assert_eq!(run_on(&words, "", &mut written), (0, String::new()));
        assert!(written.ends_with(b"}\n"));
        assert_eq!(
            run_on(&words, "", &mut FullDisk),
            (
                1,
                "Error: a failed write of the -metadata output (NCBI's serial stream reports it) \
                 is not supported by LOSAT's BLASTDBCMD\n"
                    .to_string()
            )
        );
        // A closed pipe ends the program silently outside Windows, where SIGPIPE ends NCBI's.
        if !cfg!(windows) {
            assert_eq!(run_on(&words, "", &mut ClosedPipe), (0, String::new()));
        }
    }
}
