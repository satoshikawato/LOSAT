//! Registered inputs (docs/web/abi_v2.md §4 `register`, plan §4.6 R1): the records
//! read by the program's own FASTA reader (the engine's port of NCBI BLAST+'s
//! `CBlastFastaInputSource`), kept under a handle so that one subject serves several runs; and
//! the inputs of one run read together for its argv (`register_inputs`, a run with
//! `-parse_deflines`).

use std::collections::HashMap;
use std::fmt::Write as _;
use std::sync::Mutex;

use losat::api::FastaRecord;
use losat::web_adapter::{
    read_all, FastaInputSource, ParseErrorCode, ReadError, ReaderConfig, RecordSource,
};

use crate::json;
use crate::run::Program;
use crate::scan;

pub const ROLE_QUERY: u32 = 0;
pub const ROLE_SUBJECT: u32 = 1;

/// Deflines by record index, in the order of the records.
pub type KeptDeflines = Vec<(usize, Box<[u8]>)>;

pub struct Registered {
    pub role: u32,
    /// The program whose reader read the records; only its runs use them.
    pub program: Program,
    pub records: Vec<FastaRecord>,
    /// The deflines that the BLAST archive cannot make from the title (`read_deflines`); the
    /// records themselves keep no source.
    pub deflines: KeptDeflines,
    /// The argv of the one run whose inputs `register_inputs` read for it: the records are what the
    /// CLI's reader gives that argv (the ids of the deflines, the deflines that its archive shows),
    /// and only that run takes them. `None` for the records of `register`.
    pub run: Option<Vec<String>>,
}

#[derive(Default)]
struct Store {
    next: u32,
    entries: HashMap<u32, Registered>,
}

fn store() -> &'static Mutex<Store> {
    static STORE: std::sync::OnceLock<Mutex<Store>> = std::sync::OnceLock::new();
    STORE.get_or_init(Mutex::default)
}

/// The reader of a program's input of `role`: the program's name in LOSAT's messages and
/// whether the input is protein (scan kind 2) or nucleotide (scan kind 1).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:328-330
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:51,106
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2554
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:50,106
/// BLASTN and TBLASTX read nucleotides in both roles, BLASTP proteins in both, TBLASTN a
/// protein query and nucleotide subjects, BLASTX a nucleotide query and protein subjects.
fn reader_of(program: Program, role: u32) -> (&'static str, bool) {
    match (program, role) {
        (Program::Blastn, _) => ("BLASTN", false),
        (Program::Tblastx, _) => ("TBLASTX", false),
        (Program::Tblastn, ROLE_SUBJECT) => ("TBLASTN", false),
        (Program::Tblastn, _) => ("TBLASTN", true),
        (Program::Blastp, _) => ("BLASTP", true),
        (Program::Blastx, ROLE_SUBJECT) => ("BLASTX", true),
        (Program::Blastx, _) => ("BLASTX", false),
    }
}

/// The record ID of the *register* and *scan* responses: the title up to its first space
/// (empty without a title), with bytes that are not UTF-8 as U+FFFD.
pub fn record_id(title: &[u8]) -> String {
    let end = title
        .iter()
        .position(|&byte| byte == b' ')
        .unwrap_or(title.len());
    String::from_utf8_lossy(&title[..end]).into_owned()
}

/// An error of the engine's reader as the CLI reports it, without the CLI's line end.
fn reader_error(error: ReadError) -> String {
    format!("{:#}", error.into_app_error())
        .trim_end_matches('\n')
        .to_string()
}

/// The deflines that the BLAST archive (outfmt 11) cannot make from the title, by record index: the
/// archive's Bioseq shows the defline as read (`report::archive::defline_of`), which the CLI's
/// reader keeps when the run writes the archive (`ReaderConfig::with_defline`), and else `>` and
/// the title, which the title does not give for a defline with a control byte (the title ends
/// there), white space after the `>` or at the end. The index scan (`scanned`, which `check_scan`
/// matched with the records) tells the deflines whose line is not `>` and the title (a line never
/// holds a line end, so the line without its line end is the defline as read when it is `>` and the
/// title). The reader reads each of those lines alone, keeping the defline; the letters are not
/// read again, so the registration of any input holds what a reader without the archive holds. A
/// line that the reader alone reads as another title is read in the whole input instead.
fn read_deflines(
    bytes: &[u8],
    config: ReaderConfig,
    scanned: &[scan::ScanRecord],
    records: &[FastaRecord],
) -> Result<KeptDeflines, ReadError> {
    let made = |line: &[u8], title: &[u8]| line.first() == Some(&b'>') && &line[1..] == title;
    let mut quiet = |_: &[u8]| Ok(());
    let mut deflines = Vec::new();
    for (index, (scan, record)) in scanned.iter().zip(records).enumerate() {
        let range = scan.header_offset as usize..scan.sequence_offset as usize;
        let line = &bytes[range.clone()];
        let line = line.strip_suffix(b"\n").unwrap_or(line);
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.is_empty() || made(line, &record.title) {
            continue;
        }
        let mut alone = FastaInputSource::from_bytes(&bytes[range], config.with_defline(true));
        let read = alone.next_sequence(&mut quiet).ok();
        let Some(kept) = read.filter(|kept| kept.title == record.title) else {
            return read_all_deflines(bytes, config);
        };
        if let Some(defline) = kept
            .archive_defline()
            .filter(|defline| !made(defline, &record.title))
        {
            deflines.push((index, defline.into()));
        }
    }
    Ok(deflines)
}

/// `read_deflines` over the whole input: every record read again, one at a time, keeping the
/// deflines that `>` and the title do not give.
#[cold]
fn read_all_deflines(bytes: &[u8], config: ReaderConfig) -> Result<KeptDeflines, ReadError> {
    let made = |line: &[u8], title: &[u8]| line.first() == Some(&b'>') && &line[1..] == title;
    let mut source = FastaInputSource::from_bytes(bytes, config.with_defline(true));
    let mut quiet = |_: &[u8]| Ok(());
    let mut deflines = Vec::new();
    let mut index = 0;
    while !source.end() {
        match source.next_sequence(&mut quiet) {
            Ok(record) => {
                if let Some(defline) = record
                    .archive_defline()
                    .filter(|defline| !made(defline, &record.title))
                {
                    deflines.push((index, defline.into()));
                }
                index += 1;
            }
            Err(ReadError::Parse {
                code: ParseErrorCode::Eof,
                ..
            }) => break,
            Err(error) => return Err(error),
        }
    }
    Ok(deflines)
}

/// Reads `bytes` with the program's reader for `role`, checks the records against the
/// index scan of the same bytes (plan TD-8), keeps them, and returns the handle and the
/// *register* JSON.
///
/// The scan (kind 1 or 2, `scan/ncbi.rs`) comes first: it rejects what LOSAT Web cannot
/// index (a `>?` gap line, maintainer decision 3; lines that NCBI's line reader joins), a
/// first line that NCBI's data loaders would fetch as a Seq-id, and fails on the reader's
/// own error (`CheckDataLine`), which `register` reports with the CLI's text
/// (`BLAST query error: ...`). The engine's reader then reads the records with its data
/// loaders on (port plan Q3, as the CLI without `.ncbirc` or `DATA_LOADERS`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:78-80
pub fn register(program: &str, role: u32, bytes: &[u8]) -> Result<(u32, String), String> {
    let program = Program::parse(program)?;
    let role_name = match role {
        ROLE_QUERY => "query",
        ROLE_SUBJECT => "subject",
        other => return Err(format!("unknown input role {other}")),
    };
    let (name, protein) = reader_of(program, role);
    let config = if role == ROLE_QUERY {
        ReaderConfig::query(name, protein, true)
    } else {
        ReaderConfig::subject(name, protein, true)
    };
    let read = || {
        read_all(
            &mut FastaInputSource::from_bytes(bytes, config),
            &mut |_| Ok(()),
        )
    };
    let scanned = match scan::scan_ncbi(bytes, protein) {
        Ok(scanned) => scanned,
        // The scan stops at the reader's error with the reader's message: the CLI's text.
        Err(error) => {
            return Err(match read() {
                Err(parse @ ReadError::Parse { .. }) if parse.to_string() == error => {
                    reader_error(parse)
                }
                _ => error,
            })
        }
    };
    let records = read().map_err(reader_error)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:211-216
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:146-152
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
    // A query of white space only is NCBI's `Query is Empty!` at run time; a query with
    // blank and comment lines only is not empty, and its first batch (`eEOF`, no query)
    // fails every run, so `register` fails with NCBI's error (port plan Q2).
    if role == ROLE_QUERY {
        check_query_batch(bytes, &records)?;
    }
    check_scan(role_name, &scanned, &records)?;
    let deflines = read_deflines(bytes, config, &scanned, &records).map_err(reader_error)?;
    let mut response = String::new();
    let mut store = store().lock().expect("input store");
    let handle = store.next_handle()?;
    write_records(&mut response, handle, &scanned, &records, None);
    store.entries.insert(
        handle,
        Registered {
            role,
            program,
            records,
            deflines,
            run: None,
        },
    );
    Ok((handle, response))
}

/// NCBI's `BLAST engine error: Empty CBlastQueryVector`, without the CLI's line end: the error of
/// a query input whose first batch has no query (`check_query_batch`), and of a subject input
/// without records, which ends the CLI before it reads the query (`register_inputs`).
fn empty_query_vector() -> String {
    format!(
        "{:#}",
        losat::web_adapter::engine_error("Empty CBlastQueryVector")
    )
    .trim_end_matches('\n')
    .to_string()
}

/// NCBI's error of a query input whose first batch has no query although it is not white space
/// only (blank and comment lines; `register`).
fn check_query_batch(bytes: &[u8], records: &[FastaRecord]) -> Result<(), String> {
    if records.is_empty() && !losat::web_adapter::is_blank(bytes) {
        return Err(empty_query_vector());
    }
    Ok(())
}

impl Store {
    /// A handle that no input has had.
    fn next_handle(&mut self) -> Result<u32, String> {
        self.next = self
            .next
            .checked_add(1)
            .ok_or_else(|| "input handle space exhausted".to_string())?;
        Ok(self.next)
    }
}

/// The handle and the records of a *register* response (docs/web/abi_v2.md §9), `{"handle":…,
/// "records":[…]}`: each record's index, ID (the index scan's: the title up to its first space) and
/// length, and with `shown` (`register_inputs`) the ID that the tabular outputs show for it
/// (`tabular_id`).
fn write_records(
    out: &mut String,
    handle: u32,
    scanned: &[scan::ScanRecord],
    records: &[FastaRecord],
    shown: Option<(u32, bool)>,
) {
    let _ = write!(out, "{{\"handle\":{handle},\"records\":[");
    for (index, (scan, record)) in scanned.iter().zip(records).enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(out, "{{\"index\":{index},\"id\":");
        json::string(out, &scan.id);
        let _ = write!(out, ",\"length\":{}", record.seq().len());
        if let Some((role, protein)) = shown {
            out.push_str(",\"shown_id\":");
            json::string(out, &tabular_id(record, role, protein));
        }
        out.push('}');
    }
    out.push_str("]}");
}

/// The ID that the tabular outputs (outfmt 6 and 7) show for a record of `role` read by the reader
/// of `protein` letters: the `qaccver` of a query (`FastaRecord::shown_id`; with `-parse_deflines`
/// the accession and version of the ids that the reader read, `AB000001.1` of
/// `>gb|AB000001.1| ...`) and the `saccver` of a subject (`FastaRecord::subject_tabular_id`), the
/// first two columns of the CLI's default fields. Where the engine gives no `saccver` (a record
/// whose title NCBI cannot decode, or whose ids LOSAT refuses, which a report of it shows or refuses
/// as the CLI does), the subject's `shown_id`. Bytes that are not UTF-8 become U+FFFD.
fn tabular_id(record: &FastaRecord, role: u32, protein: bool) -> String {
    let saccver = match role {
        ROLE_SUBJECT => record.subject_tabular_id(protein).ok().flatten(),
        _ => None,
    };
    match saccver {
        Some(id) => String::from_utf8_lossy(&id).into_owned(),
        None => String::from_utf8_lossy(record.shown_id()).into_owned(),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:188-213
/// Reads the query and subject FASTA bytes of the run of `words` (program name first) together, as
/// the CLI reads its `-query` and `-subject` files for that argv: the library's reading of the
/// inputs of a search (`read_search_inputs`, which the library API's `read_inputs` reads with):
/// with `-parse_deflines` the ids of the deflines, read in one scope of ids, the subjects first,
/// and the query's `-strand`; and as the argv's reader, which keeps the deflines when the argv
/// writes the BLAST archive and then refuses a ranged id as it reads it. The records are checked
/// against the index scans of the same bytes, as `register` checks them, and kept under a query
/// handle and a subject handle that only the run of `words` takes (`with_inputs`). Returns the
/// handles and the *inputs* JSON (docs/web/abi_v2.md §9): the *register* response of each, with
/// each record's `shown_id` (`tabular_id`).
///
/// The errors are the CLI's texts. Their order is the CLI's for the subjects, which the CLI reads
/// whole before any query, and for a query that the CLI reads in one batch: every query record is
/// read here before the run, so a reader's error in a later query batch of the CLI (two records of
/// one Seq-id under `-parse_deflines`, for example) ends the run before the searches of the batches
/// before it, where the CLI ends with the first error of those (a refusal of outfmt 0 in batch 1;
/// review of line WQA's D7, L1, 2026-11-19, open). The scans come first: LOSAT Web's refusals and
/// the reader's own errors, as in `register`. A scan that fails gives the error of the reading when the
/// reading fails (the CLI fails with that error, whatever the scan refused), and its own error
/// otherwise. A subject without records ends the run with NCBI's `Empty CBlastQueryVector` before
/// the query is read, as it ends the CLI (`register` registers it, and the run gives the error).
pub fn register_inputs(
    words: &[&str],
    query: &[u8],
    subject: &[u8],
) -> Result<(u32, u32, String), String> {
    let (program, command, supplied) = crate::run::parse_imported(words)?;
    // An argv with `-import_search_strategy` takes the inputs that the strategy carries and the
    // argv does not override (`Supplied`), in place of the bytes that the host gives for them (it
    // gives none: the argv has no `-query` or `-subject` of its own).
    let (query, subject) = match &supplied {
        Some(supplied) => (
            if supplied.queries.is_empty() {
                query
            } else {
                &supplied.queries
            },
            if supplied.subjects.is_empty() {
                subject
            } else {
                &supplied.subjects
            },
        ),
        None => (query, subject),
    };
    let read = || {
        losat::web_adapter::read_search_inputs(
            program.api(),
            &command,
            Some(crate::run::writes_archive(words)),
            query,
            subject,
        )
        .map_err(|error| error.to_string())
    };
    let scan = |role: u32| {
        let bytes = if role == ROLE_QUERY { query } else { subject };
        scan::scan_ncbi(bytes, reader_of(program, role).1)
    };
    let subject_scan = scan(ROLE_SUBJECT);
    if subject_scan
        .as_ref()
        .is_ok_and(|records| records.is_empty())
    {
        return Err(empty_query_vector());
    }
    let (subject_scan, query_scan) = match (subject_scan, scan(ROLE_QUERY)) {
        (Ok(subject_scan), Ok(query_scan)) => (subject_scan, query_scan),
        (subject_scan, query_scan) => {
            let refusal = subject_scan.err().or(query_scan.err()).unwrap_or_default();
            return Err(read().err().unwrap_or(refusal));
        }
    };
    let (queries, subjects) = read()?;
    check_query_batch(query, &queries)?;
    check_scan("query", &query_scan, &queries)?;
    check_scan("subject", &subject_scan, &subjects)?;
    let run: Vec<String> = words.iter().map(|word| word.to_string()).collect();
    let mut store = store().lock().expect("input store");
    let (query_handle, subject_handle) = (store.next_handle()?, store.next_handle()?);
    let mut response = String::from("{\"query\":");
    let shown = |role: u32| Some((role, reader_of(program, role).1));
    write_records(
        &mut response,
        query_handle,
        &query_scan,
        &queries,
        shown(ROLE_QUERY),
    );
    response.push_str(",\"subject\":");
    write_records(
        &mut response,
        subject_handle,
        &subject_scan,
        &subjects,
        shown(ROLE_SUBJECT),
    );
    response.push('}');
    for (handle, role, records) in [
        (query_handle, ROLE_QUERY, queries),
        (subject_handle, ROLE_SUBJECT, subjects),
    ] {
        store.entries.insert(
            handle,
            Registered {
                role,
                program,
                records,
                deflines: Vec::new(),
                run: Some(run.clone()),
            },
        );
    }
    Ok((query_handle, subject_handle, response))
}

/// Fails unless the index scan found the reader's records: the same number, and for each
/// record the same ID, length and residue counts (the stored residues upper-cased; plan
/// TD-8).
fn check_scan(
    role: &str,
    scanned: &[scan::ScanRecord],
    records: &[FastaRecord],
) -> Result<(), String> {
    // A record whose ids the reader read from its defline (`-parse_deflines`) has the title after
    // them: the scan's ID, the defline's first word, is not the title's.
    let agrees = |scan: &scan::ScanRecord, record: &FastaRecord| {
        let mut counts = [0u64; 256];
        for &byte in record.seq() {
            counts[byte.to_ascii_uppercase() as usize] += 1;
        }
        (record.believed().is_some() || scan.id == record_id(record.title()))
            && scan.length == record.seq().len() as u64
            && *scan.residue_counts == counts
    };
    if scanned.len() != records.len()
        || scanned
            .iter()
            .zip(records)
            .any(|(scan, record)| !agrees(scan, record))
    {
        return Err(format!(
            "the index scan of the {role} FASTA disagrees with the reader"
        ));
    }
    Ok(())
}

pub fn release(handle: u32) -> Result<(), String> {
    store()
        .lock()
        .expect("input store")
        .entries
        .remove(&handle)
        .map(|_| ())
        .ok_or_else(|| format!("unknown input handle {handle}"))
}

/// Calls `work` with the records of a query handle and a subject handle registered for
/// `program`, for the run of the argv `words`. For a run that writes the BLAST archive
/// (`run::writes_archive`), every record of `register` has a source while
/// `work` runs, as the CLI's reader gives each record when the run writes the archive
/// (`ReaderConfig::with_defline`): its defline where `>` and the title do not give it
/// (`read_deflines`), and a record that a range cuts keeps its whole letters for the archive's
/// Bioseq (`FastaRecord::cut`, `archive_letters`). Afterwards the records have no source again, as
/// every other run sees them.
///
/// The handles of `register_inputs` are the inputs of the run of their argv alone, both of them,
/// and the records are as its reader gave them; a run that reads the ids of the deflines
/// (`run::inputs_together`) takes no other: its report would be that of a search without the option.
pub fn with_inputs<R>(
    program: Program,
    query: u32,
    subject: u32,
    words: &[&str],
    work: impl FnOnce(&[FastaRecord], &[FastaRecord]) -> R,
) -> Result<R, String> {
    let mut store = store().lock().expect("input store");
    for (handle, role, name) in [
        (query, ROLE_QUERY, "query"),
        (subject, ROLE_SUBJECT, "subject"),
    ] {
        let entry = store
            .entries
            .get(&handle)
            .filter(|entry| entry.role == role)
            .ok_or_else(|| format!("{handle} is not a registered {name} handle"))?;
        if entry.program != program {
            return Err(format!(
                "{name} handle {handle} was registered for {}, not {}",
                entry.program.name(),
                program.name()
            ));
        }
    }
    let read_for = |handle: u32| store.entries[&handle].run.as_deref();
    let read_together = match (read_for(query), read_for(subject)) {
        (None, None) if crate::run::inputs_together(words) => {
            return Err(
                "-parse_deflines needs the inputs of register_inputs, which reads the ids of the deflines"
                    .to_string(),
            )
        }
        (None, None) => false,
        (Some(query_run), Some(subject_run)) if query_run == words && subject_run == words => true,
        _ => {
            return Err(format!(
                "query handle {query} and subject handle {subject} are not the inputs that register_inputs read for this argv"
            ))
        }
    };
    let whole = !read_together && crate::run::writes_archive(words);
    if whole {
        set_sources(&mut store, [query, subject], true);
    }
    let result = work(
        &store.entries[&query].records,
        &store.entries[&subject].records,
    );
    if whole {
        set_sources(&mut store, [query, subject], false);
    }
    Ok(result)
}

/// Gives every record of `handles` a source with its kept defline (`give`), or takes the sources
/// back.
fn set_sources(store: &mut Store, handles: [u32; 2], give: bool) {
    for handle in handles {
        let entry = store.entries.get_mut(&handle).expect("a registered handle");
        let mut kept = entry.deflines.iter().peekable();
        for (index, record) in entry.records.iter_mut().enumerate() {
            record.source = give.then(|| {
                let defline = kept
                    .next_if(|(at, _)| *at == index)
                    .map(|(_, defline)| defline.to_vec())
                    .unwrap_or_default();
                Box::new(RecordSource {
                    defline,
                    ..RecordSource::default()
                })
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // `register` stops when the scan and the reader disagree on the records.
    #[test]
    fn a_scan_that_disagrees_with_the_reader_fails_register() {
        let bytes = b">a one\nACGT\n>b\nGG\n";
        let records = read_all(
            &mut FastaInputSource::from_bytes(bytes, ReaderConfig::query("BLASTN", false, true)),
            &mut |_| Ok(()),
        )
        .unwrap();
        let scanned = scan::scan_ncbi(bytes, false).unwrap();
        assert!(check_scan("query", &scanned, &records).is_ok());

        assert!(check_scan("query", &scanned[..1], &records).is_err());
        let mut renamed = scanned.clone();
        renamed[1].id = "c".to_string();
        assert!(check_scan("query", &renamed, &records).is_err());
        let mut longer = scanned.clone();
        longer[0].length += 1;
        assert!(check_scan("query", &longer, &records).is_err());
        let mut recounted = scanned.clone();
        recounted[0].residue_counts[b'A' as usize] += 1;
        assert!(check_scan("query", &recounted, &records).is_err());
    }

    /// The `records` of a *register* response.
    fn records(program: &str, role: u32, bytes: &[u8]) -> String {
        let (_, response) = register(program, role, bytes).unwrap();
        let at = response.find("\"records\":").unwrap();
        response[at..].to_string()
    }

    // The records are those of NCBI's reader: records without residues, a first record
    // without a defline and invalid residues (skipped with a warning) are read; the title
    // ends at a byte below 0x20; the ID is the title's first word, empty without a title,
    // lossy only in JSON.
    #[test]
    fn register_reads_what_ncbi_reads() {
        assert_eq!(
            records("blastn", ROLE_QUERY, b">q\tt\nACXGU\n"),
            r#""records":[{"index":0,"id":"q","length":4}]}"#
        );
        assert_eq!(
            records("tblastx", ROLE_SUBJECT, b">s0\n>s1 x\nACGT\n"),
            r#""records":[{"index":0,"id":"s0","length":0},{"index":1,"id":"s1","length":4}]}"#
        );
        assert_eq!(
            records("blastp", ROLE_QUERY, b"ACGT\n>x\nAC\n"),
            r#""records":[{"index":0,"id":"","length":4},{"index":1,"id":"x","length":2}]}"#
        );
        // BLASTX reads a nucleotide query (X is not a nucleotide letter) and protein subjects.
        assert_eq!(
            records("blastx", ROLE_QUERY, b">q\nACXGU\n"),
            r#""records":[{"index":0,"id":"q","length":4}]}"#
        );
        assert_eq!(
            records("blastx", ROLE_SUBJECT, b">s\nACXGU\n"),
            r#""records":[{"index":0,"id":"s","length":5}]}"#
        );
        assert_eq!(
            records("tblastn", ROLE_SUBJECT, b"\n>\xff\xfe a\nACGT\n"),
            "\"records\":[{\"index\":0,\"id\":\"\u{fffd}\u{fffd}\",\"length\":4}]}"
        );
        // White space only: no record (`Query is Empty!` at run time); comment lines only:
        // no record for a subject (the run's `Empty CBlastQueryVector`) and NCBI's error for
        // a query.
        assert_eq!(
            records("blastn", ROLE_QUERY, b" \n\t\n"),
            r#""records":[]}"#
        );
        assert_eq!(
            records("blastp", ROLE_SUBJECT, b";c\n\n"),
            r#""records":[]}"#
        );
        assert_eq!(
            register("blastp", ROLE_QUERY, b";c\n\n").unwrap_err(),
            "BLAST engine error: Empty CBlastQueryVector"
        );
    }

    // LOSAT Web's rejections (the scan's) and the reader's error with the CLI's text.
    #[test]
    fn register_rejects_what_ncbi_rejects_and_what_web_cannot_index() {
        let gap = register("blastn", ROLE_SUBJECT, b">s\nACGT\n>?10\nACGT\n").unwrap_err();
        assert!(gap.starts_with("line 3 is a gap line ('>?')"), "{gap}");
        assert!(gap.ends_with("not supported in the browser"), "{gap}");
        let seq_id = register("tblastx", ROLE_QUERY, b"AB123456\nACGT\n").unwrap_err();
        assert!(seq_id.contains("not supported in the browser"), "{seq_id}");
        let reader =
            register("blastp", ROLE_QUERY, b">q\n1234567890123456789012345\n").unwrap_err();
        assert_eq!(
            reader,
            "BLAST query error: CFastaReader: Near line 2, there's a line that doesn't look like plausible data, but it's not marked as defline or comment."
        );
    }

    // The handles of `register_inputs` are the inputs of the run of their argv alone, and a run with
    // `-parse_deflines` takes no others (line WQA's sweep, D7, 2026-11-19): handles of `register`
    // would bring records without the ids, a mix would pair inputs read apart.
    #[test]
    fn inputs_read_together_belong_to_their_argv() {
        let words = ["blastp", "-query", "q", "-subject", "s", "-parse_deflines"];
        // Test protein: the first residues of E. coli AsnA (RefSeq NP_418200.1).
        let (query_bytes, subject_bytes) =
            (b">gb|AB000001.1| q\nMKTAYIAKQR\n", b">s1 t\nMKTAYIAKQR\n");
        let (query, subject, response) =
            register_inputs(&words, query_bytes, subject_bytes).unwrap();
        assert!(
            response.contains(r#""id":"gb|AB000001.1|","length":10,"shown_id":"AB000001.1""#),
            "{response}"
        );
        assert!(
            response.contains(r#""id":"s1","length":10,"shown_id":"s1""#),
            "{response}"
        );
        let ids = |queries: &[FastaRecord], subjects: &[FastaRecord]| {
            (
                queries[0].believed().is_some(),
                subjects[0].believed().is_some(),
            )
        };
        assert_eq!(
            with_inputs(Program::Blastp, query, subject, &words, ids),
            Ok((true, true))
        );
        let other = [
            "blastp",
            "-query",
            "q",
            "-subject",
            "s",
            "-parse_deflines",
            "-seg",
            "yes",
        ];
        let error = with_inputs(Program::Blastp, query, subject, &other, |_, _| ()).unwrap_err();
        assert!(
            error.contains("are not the inputs that register_inputs read for this argv"),
            "{error}"
        );
        let (plain_query, _) = register("blastp", ROLE_QUERY, query_bytes).unwrap();
        let (plain_subject, _) = register("blastp", ROLE_SUBJECT, subject_bytes).unwrap();
        let mixed =
            with_inputs(Program::Blastp, plain_query, subject, &words, |_, _| ()).unwrap_err();
        assert!(
            mixed.contains("are not the inputs that register_inputs read for this argv"),
            "{mixed}"
        );
        assert_eq!(
            with_inputs(Program::Blastp, plain_query, plain_subject, &words, |_, _| ()).unwrap_err(),
            "-parse_deflines needs the inputs of register_inputs, which reads the ids of the deflines"
        );
        // Without the option the handles of `register` are the inputs, as before.
        let plain = &words[..5];
        assert!(with_inputs(
            Program::Blastp,
            plain_query,
            plain_subject,
            plain,
            |_, _| ()
        )
        .is_ok());
        for handle in [query, subject, plain_query, plain_subject] {
            release(handle).unwrap();
        }
    }

    // Another program's handle cannot bring a record to a run.
    #[test]
    fn a_handle_belongs_to_its_program() {
        let (query, _) = register("blastp", ROLE_QUERY, b">q t\nACXGT\n").unwrap();
        let (subject, _) = register("blastn", ROLE_SUBJECT, b">s\nACGT\n").unwrap();
        let words = ["blastn", "-query", "q", "-subject", "s"];
        let error = with_inputs(Program::Blastn, query, subject, &words, |_, _| ()).unwrap_err();
        assert!(
            error.contains("registered for blastp, not blastn"),
            "{error}"
        );
        let words = ["blastp", "-query", "q", "-subject", "s"];
        assert!(with_inputs(Program::Blastp, query, subject, &words, |_, _| ()).is_err());
    }
}
