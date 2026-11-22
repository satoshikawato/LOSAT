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
mod bioseq;
mod fasta;
mod format;
mod letters;
mod list;
mod metadata;
mod range;
mod retrieve;
mod taxonomy;

use crate::wording::DATABASE_ERROR;
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
    let registry = match &settings.blastdb {
        Ok(value) => value.as_deref(),
        Err(reason) => return diag.refuse(&refusal(reason)),
    };
    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1179-1180
    let mut out = match args.output_file("out", out) {
        Ok(out) => out,
        Err(error) => {
            // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:177-180
            diag.post_bytes(
                &[b"Command line argument error: ".as_slice(), error.message()].concat(),
            );
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
            return taxonomy::tax_info(&db, args, diag, &mut out);
        }
        let protein = db.seqtype() == SeqType::Protein;
        return retrieve::process(&db, args, diag, stdin, &mut out, None, protein);
    }
    // NCBI prints the date bytes of the index file; the reader keeps them as text
    // (`SeqDb::date`), with U+FFFD for bytes that are not UTF-8.
    if db.date().contains('\u{fffd}') {
        return diag.refuse(&refusal(&format!(
            "{name}: the date of the database is not UTF-8 text"
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
    let mut text = Vec::new();
    write_info_head(&mut text, &db, exact);
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
// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:524-531
/// Opens the database `name` (the `-db` value as given, one name or several: blastdbcmd does not
/// trim it) of the type `dbtype` (`prot`, `nucl`, or `guess`: protein first, then nucleotide),
/// with the list of its volumes for `-info`, and `prepare`s it (the taxid list of `-target_only`, which NCBI's
/// constructor applies). With `guess`, an exception of NCBI's reader
/// (`SeqDbError::is_seqdb_exception`) moves on to the nucleotide database, and the nucleotide
/// error is the one reported; a refusal of LOSAT's is reported at once, since NCBI may open that
/// protein database.
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
            Err(error) if types.peek().is_some() && error.is_seqdb_exception() => {}
            opened => return opened,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:247-250
/// A reader error: NCBI's own (`SeqDbError::ncbi_message`) as `BLAST Database error: <text>`
/// with exit status 2, another exception of NCBI's (`SeqDbError::ncbi_exception`: an alias
/// file's value that does not convert) as `Error: <what()>` with exit status 255, else LOSAT's
/// refusal.
fn database_error(diag: &mut Diag, error: SeqDbError) -> i32 {
    if let Some(exception) = error.ncbi_exception() {
        diag.post_bytes(&[b"Error: ".as_slice(), &exception.what()].concat());
        return 255;
    }
    match error.ncbi_message() {
        Some(text) => {
            diag.post_bytes(&[DATABASE_ERROR.as_bytes(), b" ", &text].concat());
            2
        }
        None => diag.refuse(&crate::algorithm::common::database::losat_reader_error(
            &error, PROGRAM,
        )),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:533-573
/// The text of `-info` for `db` up to its filtering algorithms.
///
/// The counts are `GetNumSeqs` and `GetTotalLength` of the database opened without a list:
/// the walks of the alias tree, or the scan of the included OIDs that an alias file's
/// `OIDLIST` asks for (`SeqDb::printed_num_seqs`). `exact` is `-exact_length`'s exact total and
/// the longest length after it (`SeqDb::exact_totals`), which differ from the others only where
/// NCBI scanned the totals (the OID ranges `FIRST_OID`, `LAST_OID` do not make it scan).
fn write_info_head(out: &mut Vec<u8>, db: &SeqDb, exact: Option<(u64, u32)>) {
    let (total, longest) = exact.unwrap_or((db.printed_total_length(), db.max_length()));
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
        with_commas(total),
        db.date(),
        // `GetMaxLength` is NCBI's `int` of the index files' longest length (of the scan's).
        int_with_commas(longest as i32),
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
#[path = "tests.rs"]
mod tests;
