//! `losat blastdbcmd`: NCBI BLAST+ 2.17.0 blastdbcmd, of which LOSAT serves `-info` (the
//! information of a database: `-db`, `-dbtype`, `-exact_length`, `-out`, `-logfile`), the
//! retrieval of sequences (`retrieve`: `-entry all`, `-entry`, `-entry_batch`, `-ipg`,
//! `-ipg_batch`, `-taxids`, `-taxidlist`, every `-outfmt` but the ASN.1 texts, `-range`,
//! `-strand`, `-line_length`, `-ctrl_a`, `-target_only`, `-get_dups`, `-long_seqids`), the
//! taxids of a database (`taxonomy`: `-tax_info`), its metadata (`metadata`: `-metadata`,
//! `-metadata_output_prefix`), the databases of a directory (`list`: `-list`,
//! `-list_outfmt`, `-recursive`, `-remove_redundant_dbs`), `-show_blastdb_search_path`, `-h`,
//! `-help` and `-version`. Its command line is NCBI's (`args`, NCBI's argument reader): `main`
//! hands it the words after `blastdbcmd`. The masks of `makeblastdb -mask_data`
//! are read (`-info`, `-mask_sequence_with`, `%m`). A `-db` of several names is read as NCBI's
//! reader reads it, as one database (`SeqDb::open_in`).

mod args;
mod bioseq;
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
use crate::ncbi_app::files::Cin;
use crate::ncbi_app::{database_error, Diag, Tool};
use crate::ncbi_str::{int_with_commas, with_commas};
use crate::report::ostream::CutBuffer;
use crate::seqdb::{SearchPath, SeqDb, SeqDbError, SeqType, TaxidList};

/// The name in LOSAT's refusals.
const PROGRAM: &str = "BLASTDBCMD";

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:907-909
/// The text of `blastdbcmd -help` (NCBI 2.17.0's, `CArgDescriptions::PrintUsage` of the
/// descriptions of `CBlastDBCmdApp::Init`). Its part before `OPTIONAL ARGUMENTS` is the usage
/// that `-h` and the argument errors print.
const HELP: &str = include_str!("help.txt");

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1175-1177
/// blastdbcmd as NCBI's application runs it (`ncbi_app`): its `-version` states LOSAT's version
/// where NCBI's states BLAST+'s (`CBlastDBCmdApp` sets `CBlastVersion`), its diagnostics carry the
/// prefix `[blastdbcmd]`.
static TOOL: Tool = Tool {
    program: "blastdbcmd",
    refusal_name: PROGRAM,
    help: HELP,
    post_prefix: Some("blastdbcmd"),
    components: &[],
    help_full: None,
    short_description: None,
    dry_run: None,
    closed_log_ends_run: true,
};

/// Runs `blastdbcmd` with `words` (the command line after `losat blastdbcmd`), after its
/// application start `start`, on the process's standard streams; the exit status.
pub fn run(
    words: &[OsString],
    start: &crate::blastinput::ncbi_environment::ApplicationStart,
) -> i32 {
    crate::ncbi_app::run(&TOOL, words, start, run_body)
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
    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1179-1180
    // `-metadata` (when no mode before it in NCBI's `Run` is given) writes through its own model
    // of NCBI's serializer and stream (`metadata::metadata`); every other mode through NCBI's
    // stream of `-out` (`files::OutputFile`).
    let metadata_mode = args.given("metadata")
        && !["show_blastdb_search_path", "list", "info"]
            .into_iter()
            .any(|mode| args.given(mode));
    let arg_error = |diag: &mut Diag, error: crate::ncbi_app::args::ArgException| {
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-180
        diag.post_arg_error(error.message());
        1
    };
    if metadata_mode {
        let device = match args.output_device("out", out) {
            Ok(device) => device,
            Err(error) => return arg_error(diag, error),
        };
        if let crate::ncbi_app::files::OutputDevice::Refused(reason) = &device {
            return diag.refuse(&output_refusal(reason));
        }
        let path = SearchPath::from_environment(settings.blastdb.as_deref());
        let name = args
            .bytes("db")
            .unwrap_or(std::borrow::Cow::Borrowed(b"nr"));
        let dbtype = args.value("dbtype").unwrap_or("guess");
        let (db, _) = match open_database(&name, dbtype, &path, settings, &|_| Ok(())) {
            Ok(opened) => opened,
            Err(error) => return database_error(diag, error),
        };
        return metadata::metadata(&db, &name, args, &path, diag, device);
    }
    let mut out_file = match args.output_file("out", out) {
        Ok(out) => out,
        Err(error) => return arg_error(diag, error),
    };
    // LOSAT Web's refusal of an output that it does not write (`OpenFiles::create`), as soon as
    // `-out` is opened, before the database (windowmasker's and makeprofiledb's place and words).
    if let Some(reason) = out_file.refusal() {
        return diag.refuse(&output_refusal(reason));
    }
    // NCBI's `ostream` of `-out` (`files::OutputFile`): a failed write drops the later ones, and
    // a closed pipe ends the program where SIGPIPE ends NCBI's; the texts on standard error come
    // after NCBI's writes of it (`Diag`).
    let out = &mut out_file;
    let path = SearchPath::from_environment(settings.blastdb.as_deref());
    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1181-1184
    if args.given("show_blastdb_search_path") {
        // The path's insertion, then `NcbiEndl`.
        let mut text = CutBuffer::default();
        text.insertion(path.as_bytes());
        text.endl();
        let _ = text.write_to(out);
        return 0;
    }
    if args.given("list") {
        return list::list(args, &path, diag, out);
    }
    // The database's name is the bytes of NCBI's `std::string`.
    let name = args
        .bytes("db")
        .unwrap_or(std::borrow::Cow::Borrowed(b"nr"));
    let dbtype = args.value("dbtype").unwrap_or("guess");
    let selected = ["info", "metadata", "tax_info"]
        .into_iter()
        .any(|selector| args.given(selector));
    if !selected && (args.given("taxids") || args.given("taxidlist")) {
        // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:358-369
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
        let (db, _) = match open_database(&name, dbtype, &path, settings, &restrict) {
            Ok(opened) => opened,
            Err(error) => return database_error(diag, error),
        };
        let dbtype_protein = dbtype == "prot";
        return retrieve::process(&db, args, diag, stdin, out, Some(&tax_ids), dbtype_protein);
    }
    let (db, volumes) = match open_database(&name, dbtype, &path, settings, &|_| Ok(())) {
        Ok(opened) => opened,
        Err(error) => return database_error(diag, error),
    };
    if !args.given("info") {
        if args.given("tax_info") {
            return taxonomy::tax_info(&db, args, diag, out);
        }
        let protein = db.seqtype() == SeqType::Protein;
        return retrieve::process(&db, args, diag, stdin, out, None, protein);
    }
    // NCBI prints the date bytes of the index file; the reader keeps them as text
    // (`SeqDb::date`), with U+FFFD for bytes that are not UTF-8.
    if db.date().contains('\u{fffd}') {
        return diag.refuse(&refusal(&format!(
            "{}: the date of the database is not UTF-8 text",
            String::from_utf8_lossy(&name)
        )));
    }
    let exact = if args.given("exact_length") {
        match db.exact_totals() {
            Ok(totals) => Some(totals),
            Err(error) => return database_error(diag, error),
        }
    } else {
        None
    };
    let mut text = CutBuffer::default();
    write_info_head(&mut text, &db, exact);
    let masks = db
        .mask_algorithms()
        .and_then(|algorithms| algorithms.descriptions());
    let tail = match &masks {
        Ok(masks) => {
            text.insertion(masks.as_bytes());
            write_info_volumes(&mut text, &volumes);
            None
        }
        Err(error) => Some(error),
    };
    // NCBI's stream ignores a failed write (no exception is enabled on it), so does this.
    // The lines before a bad mask column are written when the error is NCBI's own; its text
    // comes after NCBI's writes of them (`Diag`), and a closed pipe ends the program there.
    if tail.is_none_or(|error| error.ncbi_message().is_some()) && text.write_to(out).is_err() {
        return 0;
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

/// LOSAT's refusal of a `-out` that LOSAT Web does not write (`OutputFile::Refused`), in the words
/// of the other programs' refusal of it (`seqmasks_io::run::unsupported`), not a database's.
fn output_refusal(reason: &str) -> String {
    crate::seqmasks_io::run::unsupported(reason, PROGRAM)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:110-153
// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:524-531
/// Opens the database `name` (the `-db` value as given, one name or several: blastdbcmd does not
/// trim it) of the type `dbtype` (`prot`, `nucl`, or `guess`: protein first, then nucleotide),
/// with the list of its volumes for `-info`, and `prepare`s it (the taxid list of `-target_only`, which NCBI's
/// constructor applies), its reads of the definition lines under the checks of `settings`
/// (`SeqDb::refuse_unchecked_reads`). `guess` is NCBI's guess of the type (`SeqDb::open_guessing`).
fn open_database(
    name: &[u8],
    dbtype: &str,
    path: &SearchPath,
    settings: &ApplicationSettings,
    prepare: &dyn Fn(&mut SeqDb) -> Result<(), SeqDbError>,
) -> Result<(SeqDb, Vec<Vec<u8>>), SeqDbError> {
    let types: &[SeqType] = match dbtype {
        "prot" => &[SeqType::Protein],
        "nucl" => &[SeqType::Nucleotide],
        _ => &[SeqType::Protein, SeqType::Nucleotide],
    };
    SeqDb::open_guessing(types, |seqtype| {
        // The alias files' counts fail where NCBI's constructor reaches them, which the list of
        // `-target_only` moves (`SeqDb::check_alias_values`).
        let mut db = SeqDb::open_search_in(name, seqtype, path)?;
        let volumes = SeqDb::info_volumes(name, seqtype, path)?;
        db.refuse_unchecked_reads(&settings.read_verification);
        prepare(&mut db)?;
        db.check_alias_values()?;
        Ok((db, volumes))
    })
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:533-573
/// The text of `-info` for `db` up to its filtering algorithms.
///
/// The counts are `GetNumSeqs` and `GetTotalLength` of the database opened without a list:
/// the walks of the alias tree, or the scan of the included OIDs that an alias file's
/// `OIDLIST` asks for (`SeqDb::printed_num_seqs`). `exact` is `-exact_length`'s exact total and
/// the longest length after it (`SeqDb::exact_totals`), which differ from the others only where
/// NCBI scanned the totals (the OID ranges `FIRST_OID`, `LAST_OID` do not make it scan).
///
/// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:545-558
/// Each string and number is one of NCBI's insertions, each line end an `endl`.
fn write_info_head(out: &mut CutBuffer, db: &SeqDb, exact: Option<(u64, u32)>) {
    let (total, longest) = exact.unwrap_or((db.printed_total_length(), db.max_length()));
    let letters: &[u8] = match db.seqtype() {
        SeqType::Protein => b"residues",
        SeqType::Nucleotide => b"bases",
    };
    out.insertion(b"Database: ");
    out.insertion(db.title());
    out.endl();
    out.insertion(b"\t");
    out.insertion(int_with_commas(db.printed_num_seqs()).as_bytes());
    out.insertion(b" sequences; ");
    out.insertion(with_commas(total).as_bytes());
    out.insertion(b" total ");
    out.insertion(letters);
    out.endl();
    out.endl();
    out.insertion(b"Date: ");
    out.insertion(db.date().as_bytes());
    out.insertion(b"\tLongest sequence: ");
    // `GetMaxLength` is NCBI's `int` of the index files' longest length (of the scan's).
    out.insertion(int_with_commas(longest as i32).as_bytes());
    out.insertion(b" ");
    out.insertion(letters);
    out.endl();
    out.endl();
    out.insertion(b"BLASTDB Version: ");
    out.insertion(if db.is_version5() { b"5" } else { b"4" });
    out.endl();
}

/// The text of `-info` after the filtering algorithms: the `Volumes:` are `volumes`.
///
/// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:566-572
fn write_info_volumes(out: &mut CutBuffer, volumes: &[Vec<u8>]) {
    out.endl();
    out.insertion(b"Volumes:");
    out.endl();
    for volume in volumes {
        out.insertion(b"\t");
        out.insertion(volume);
        out.endl();
    }
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
