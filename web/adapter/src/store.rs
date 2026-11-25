//! Registered inputs (docs/web/abi_v2.md §4 `register`, plan §4.6 R1): the records
//! read by the program's own FASTA reader (the engine's port of NCBI BLAST+'s
//! `CBlastFastaInputSource`), kept under a handle so that one subject serves several runs; and
//! the inputs of one run read together for its argv (`register_inputs`, a run with
//! `-parse_deflines`).

use std::collections::HashMap;
use std::fmt::Write as _;
use std::sync::Mutex;

use losat::api::{FastaRecord, QueryEnd};
use losat::web_adapter::{
    read_all, read_query_records, FastaInputSource, ParseErrorCode, ReadError, ReaderConfig,
    RecordSource,
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
    /// How the reading of a query ended after `records` (`QueryEnd::Input` for a subject): a run
    /// meets it where the CLI's query batches meet it, a reader's error in the batch that reads
    /// the record after `records` (line WQA's D10, 2026-11-19).
    pub end: QueryEnd,
    /// The deflines that the BLAST archive cannot make from the title (`read_deflines`); the
    /// records themselves keep no source.
    pub deflines: KeptDeflines,
    /// The argv of the one run whose inputs `register_inputs` read for it: the records are what the
    /// CLI's reader gives that argv (the ids of the deflines, the deflines that its archive shows),
    /// and only that run takes them. `None` for the records of `register`.
    pub run: Option<Vec<String>>,
    /// Whether the reader read the records with NCBI's data loaders on: then it refused a first
    /// line that NCBI would fetch as a Seq-id, and the records are those of a reader with them off
    /// too; read with them off, such a line is data, which a run with them on refuses
    /// (`with_inputs`).
    pub data_loaders: bool,
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

/// The reader of `program`'s input of `role` (`reader_of`), with NCBI's data loaders on or off
/// (`data_loaders`: the run's `ApplicationSettings::data_loaders`, which *validation* reports):
/// what `register` reads with and `scan_begin` scans with.
pub fn reader_config(program: Program, role: u32, data_loaders: bool) -> ReaderConfig {
    let (name, protein) = reader_of(program, role);
    if role == ROLE_QUERY {
        ReaderConfig::query(name, protein, data_loaders)
    } else {
        ReaderConfig::subject(name, protein, data_loaders)
    }
}

/// The name of an input role (`ROLE_QUERY`, `ROLE_SUBJECT`), or the error of another number.
pub fn role_name(role: u32) -> Result<&'static str, String> {
    match role {
        ROLE_QUERY => Ok("query"),
        ROLE_SUBJECT => Ok("subject"),
        other => Err(format!("unknown input role {other}")),
    }
}

/// An error of the engine's reader as the CLI reports it, without the CLI's line end.
pub(crate) fn reader_error(error: ReadError) -> String {
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
            return read_all_deflines(bytes, config, records.len());
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

/// `read_deflines` over the whole input: the `count` records read again, one at a time, keeping the
/// deflines that `>` and the title do not give.
#[cold]
fn read_all_deflines(
    bytes: &[u8],
    config: ReaderConfig,
    count: usize,
) -> Result<KeptDeflines, ReadError> {
    let made = |line: &[u8], title: &[u8]| line.first() == Some(&b'>') && &line[1..] == title;
    let mut source = FastaInputSource::from_bytes(bytes, config.with_defline(true));
    let mut quiet = |_: &[u8]| Ok(());
    let mut deflines = Vec::new();
    let mut index = 0;
    while index < count && !source.end() {
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

/// Reads `bytes` with the program's reader for `role`, with NCBI's data loaders on or off
/// (`data_loaders`, the run's setting: `[BLAST] DATA_LOADERS` of its registry), checks the records
/// against the index scan of the same bytes (plan TD-8), keeps them, and returns the handle and the
/// *register* JSON.
///
/// The scan (`scan/ncbi.rs`, with the same reader) comes first: it fails on the reader's refusals,
/// in the CLI's words (a first line that NCBI's data loaders would fetch as a Seq-id, a record
/// longer than LOSAT reads), before the reader holds such a record, and ends at the reader's own
/// error (`CheckDataLine`). The engine's reader then reads the records: a subject's every record,
/// or the reader's own error with the CLI's text (`LOSAT query error: ...`), as the CLI reads every
/// subject before its first query; a query's records as the CLI's query batches read them
/// (`read_query_records`): the records before the reader's error, and the error (`QueryEnd`),
/// which a run meets in the batch that reads the record, after the reports of the batches before,
/// and which the response gives (`reader_error`, line WQA's D10, 2026-11-19).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:78-80
pub fn register(
    program: &str,
    role: u32,
    data_loaders: bool,
    bytes: &[u8],
) -> Result<(u32, String), String> {
    let program = Program::parse(program)?;
    let role_name = role_name(role)?;
    let config = reader_config(program, role, data_loaders);
    // The reader's refusal, in its words: the reader refuses there too, after it held what it read
    // of the record (a record of 2^31 letters does not fit in the browser).
    let scanned = scan::scan_ncbi(bytes, config).map_err(|error| error.message().to_string())?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:211-216
    // A query of white space only is NCBI's `Query is Empty!` at run time; a query of blank and
    // comment lines only has a first batch without a query, NCBI's `Empty CBlastQueryVector`
    // at run time too (`QueryEnd::BlankLines`).
    let (records, end) = if role == ROLE_QUERY {
        read_query_records(bytes, config)
    } else {
        let records = read_all(
            &mut FastaInputSource::from_bytes(bytes, config),
            &mut |_| Ok(()),
        )
        .map_err(reader_error)?;
        (records, QueryEnd::Input)
    };
    check_scan(role_name, config.protein, &scanned, &records, &end)?;
    let deflines =
        read_deflines(bytes, config, &scanned.records, &records).map_err(reader_error)?;
    let mut response = String::new();
    let mut store = store().lock().expect("input store");
    let handle = store.next_handle()?;
    write_records(
        &mut response,
        handle,
        &scanned.records,
        &records,
        &end,
        None,
    );
    store.entries.insert(
        handle,
        Registered {
            role,
            program,
            records,
            end,
            deflines,
            run: None,
            data_loaders,
        },
    );
    Ok((handle, response))
}

/// NCBI's `BLAST engine error: Empty CBlastQueryVector`, without the CLI's line end: the error of a
/// subject input without records, which ends the CLI before it reads the query
/// (`register_inputs`).
fn empty_query_vector() -> String {
    format!(
        "{:#}",
        losat::web_adapter::engine_error("Empty CBlastQueryVector")
    )
    .trim_end_matches('\n')
    .to_string()
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
/// (`tabular_id`). The records are the scan's: past a reader's error that the scan does not meet
/// (an id under `-parse_deflines`), the scan's records that no run reads, with the scan's length and
/// no `shown_id` (`null`), as the host's record table lists them. A reader's error that ends the
/// query's records is `reader_error`, the CLI's text (line WQA's D10, 2026-11-19).
fn write_records(
    out: &mut String,
    handle: u32,
    scanned: &[scan::ScanRecord],
    records: &[FastaRecord],
    end: &QueryEnd,
    shown: Option<(u32, bool)>,
) {
    let _ = write!(out, "{{\"handle\":{handle},\"records\":[");
    for (index, scan) in scanned.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let record = records.get(index);
        let _ = write!(out, "{{\"index\":{index},\"id\":");
        json::string(out, &scan.id);
        let length = record.map_or(scan.length, |record| record.seq().len() as u64);
        let _ = write!(out, ",\"length\":{length}");
        if let Some((role, protein)) = shown {
            out.push_str(",\"shown_id\":");
            match record {
                Some(record) => json::string(out, &tabular_id(record, role, protein)),
                None => out.push_str("null"),
            }
        }
        out.push('}');
    }
    out.push(']');
    if let QueryEnd::Error(error) = end {
        out.push_str(",\"reader_error\":");
        json::string(out, &String::from_utf8_lossy(error.message()));
    }
    out.push('}');
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
/// The errors are the CLI's texts, in the CLI's order. A subject's error fails here, as the CLI
/// reads every subject before any query. The query's records are those that the CLI's query
/// batches read before a reader's error (a refused id under `-parse_deflines`, a
/// line that is not plausible data), and the error is kept with the query handle (`QueryEnd`): the
/// run meets it in the batch that reads the record, after the reports of the batches before, as
/// the CLI does (review of line WQA's D7, L1, 2026-11-19; line WQA's D10). The reader has NCBI's
/// data loaders on or off as the argv's registry sets them (`run::data_loaders`: `[BLAST]
/// DATA_LOADERS` of `-conffile`). The scans come first: the reader's refusals, as in `register`.
/// A scan that fails gives the error of the reading when the reading fails or ends with an error
/// (the CLI fails with that error, whatever the scan found), and its own error otherwise; where the
/// scan refused a record (a Seq-id first line, a record longer than LOSAT reads), the reading stops
/// before that record, as the CLI's reader does there (line WQA's D13, 2026-11-20). A subject
/// without records ends the run with NCBI's `Empty CBlastQueryVector` before the query is read, as
/// it ends the CLI (`register` registers it, and the run gives the error).
///
/// A run of `-db` (`run::searches_database`) has no subject input: its subject bytes must be empty,
/// the query is read alone (with the ids of its deflines and the argv's data loaders), and the subject handle is 0, with `null`
/// for the subject in the JSON (line WQA D12a, 2026-11-20).
pub fn register_inputs(
    words: &[&str],
    query: &[u8],
    subject: &[u8],
) -> Result<(u32, u32, String), String> {
    let (program, command) = crate::run::parse(words)?;
    let database = crate::run::searches_database(words);
    if database && !subject.is_empty() {
        return Err("a run of -db has no subject input: give no subject bytes".to_string());
    }
    let data_loaders = crate::run::data_loaders(words)?;
    let read = |query: &[u8], subject: &[u8]| {
        losat::web_adapter::read_search_inputs(
            program.api(),
            &command,
            Some(crate::run::keeps_whole_records(words)),
            data_loaders,
            query,
            subject,
        )
        .map_err(|error| error.to_string())
    };
    let parse_deflines = words.iter().any(|word| *word == "-parse_deflines");
    let scan = |role: u32| {
        let bytes = if role == ROLE_QUERY { query } else { subject };
        scan::scan_ncbi(
            bytes,
            reader_config(program, role, data_loaders).with_parse_deflines(parse_deflines),
        )
    };
    // A database is read by the run, not here: its subject scan is empty and is not checked.
    let subject_scan = if database {
        Ok(scan::Scanned {
            records: Vec::new(),
            stop: None,
        })
    } else {
        scan(ROLE_SUBJECT)
    };
    if !database
        && subject_scan
            .as_ref()
            .is_ok_and(|scanned| scanned.records.is_empty() && scanned.stop.is_none())
    {
        return Err(empty_query_vector());
    }
    let (subject_scan, query_scan) = match (subject_scan, scan(ROLE_QUERY)) {
        (Ok(subject_scan), Ok(query_scan)) => (subject_scan, query_scan),
        (Err(scan::ScanError::Refused { message, record }), _) => {
            // The CLI reads the subjects whole first: the query is not read.
            return Err(read(b"", &subject[..record as usize])
                .err()
                .unwrap_or(message));
        }
        (Ok(_), Err(scan::ScanError::Refused { message, record })) => {
            // The query records before the refused one: an error that ends them comes first.
            return Err(match read(&query[..record as usize], subject) {
                Err(error) => error,
                Ok(records) => match records.query_end {
                    QueryEnd::Error(error) => String::from_utf8_lossy(error.message()).into_owned(),
                    QueryEnd::Input | QueryEnd::BlankLines => message,
                },
            });
        }
    };
    let losat::api::InputRecords {
        query: queries,
        query_end,
        subject: subjects,
    } = read(query, subject)?;
    check_scan(
        "query",
        reader_of(program, ROLE_QUERY).1,
        &query_scan,
        &queries,
        &query_end,
    )?;
    check_scan(
        "subject",
        reader_of(program, ROLE_SUBJECT).1,
        &subject_scan,
        &subjects,
        &QueryEnd::Input,
    )?;
    let run: Vec<String> = words.iter().map(|word| word.to_string()).collect();
    let mut store = store().lock().expect("input store");
    let query_handle = store.next_handle()?;
    let subject_handle = if database { 0 } else { store.next_handle()? };
    let mut response = String::from("{\"query\":");
    let shown = |role: u32| Some((role, reader_of(program, role).1));
    write_records(
        &mut response,
        query_handle,
        &query_scan.records,
        &queries,
        &query_end,
        shown(ROLE_QUERY),
    );
    response.push_str(",\"subject\":");
    if database {
        response.push_str("null");
    } else {
        write_records(
            &mut response,
            subject_handle,
            &subject_scan.records,
            &subjects,
            &QueryEnd::Input,
            shown(ROLE_SUBJECT),
        );
    }
    response.push('}');
    let inputs = [
        (query_handle, ROLE_QUERY, queries, query_end),
        (subject_handle, ROLE_SUBJECT, subjects, QueryEnd::Input),
    ];
    for (handle, role, records, end) in inputs.into_iter().filter(|(handle, ..)| *handle != 0) {
        store.entries.insert(
            handle,
            Registered {
                role,
                program,
                records,
                end,
                deflines: Vec::new(),
                run: Some(run.clone()),
                data_loaders,
            },
        );
    }
    Ok((query_handle, subject_handle, response))
}

/// Fails unless the index scan found the reader's records: for each record the same ID, length and
/// residue counts (the stored residues upper-cased; plan TD-8) — or, when a `-parse_deflines`
/// Seq-id named the other molecule (`BelievedIds::named_molecule`), the exact alphabet relationship
/// between the scan's program alphabet and the reader's Seq-id alphabet — and the same number, or
/// where the reader's error ended the records (`end`), at least as many: the scan's records past
/// an error that it does not meet (an id under `-parse_deflines`), before its own stop
/// (`Scanned::stop`) or without one, and as many where the scan met the reader's own error. A
/// reader that stops where the scan reads on has no error (review M1, 2026-11-20).
fn check_scan(
    role: &str,
    protein: bool,
    scanned: &scan::Scanned,
    records: &[FastaRecord],
    end: &QueryEnd,
) -> Result<(), String> {
    let stopped = matches!(end, QueryEnd::Error(_));
    let count = match (stopped, &scanned.stop) {
        (true, _) => scanned.records.len() >= records.len(),
        (false, Some(_)) => false,
        (false, None) => scanned.records.len() == records.len(),
    };
    // A record whose ids the reader read from its defline (`-parse_deflines`) has the title after
    // them: the scan's ID, the defline's first word, is not the title's.
    let agrees = |scan: &scan::ScanRecord, record: &FastaRecord| {
        let mut counts = [0u64; 256];
        for &byte in record.seq() {
            counts[byte.to_ascii_uppercase() as usize] += 1;
        }
        let other_molecule = record
            .believed()
            .is_some_and(|ids| ids.named_molecule == Some(!protein));
        let letters = if !other_molecule {
            scan.length == record.seq().len() as u64 && *scan.residue_counts == counts
        } else if protein {
            // Protein program reading a nucleotide-named Seq-id: IUPACna letters stay as-is,
            // `>?` gap lines emit `X` in both while data-line `X` is dropped, and every other
            // protein-only letter is dropped.
            record.seq().len() as u64 <= scan.length
                && counts
                    .iter()
                    .zip(scan.residue_counts.iter())
                    .enumerate()
                    .all(|(byte, (&kept, &scanned))| match byte as u8 {
                        b'A' | b'B' | b'C' | b'D' | b'G' | b'H' | b'K' | b'M' | b'N' | b'R'
                        | b'S' | b'T' | b'U' | b'V' | b'W' | b'Y' => kept == scanned,
                        b'X' => kept <= scanned,
                        _ => kept == 0,
                    })
        } else {
            // Nucleotide program reading a protein-named Seq-id: non-IUPACna letters become `N`.
            record.seq().len() as u64 >= scan.length
                && counts
                    .iter()
                    .zip(scan.residue_counts.iter())
                    .enumerate()
                    .all(|(byte, (&kept, &scanned))| {
                        if byte == b'N' as usize {
                            kept >= scanned
                        } else {
                            kept == scanned
                        }
                    })
        };
        (record.believed().is_some() || scan.id == record_id(record.title())) && letters
    };
    if !count
        || scanned
            .records
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
/// `program`, for the run of the argv `words`. For a run that writes the BLAST archive or a search
/// strategy file (`run::keeps_whole_records`), every record of `register` has a source while
/// `work` runs, as the CLI's reader gives each record when the run writes either
/// (`ReaderConfig::with_defline`): its defline where `>` and the title do not give it
/// (`read_deflines`), and a record that a range cuts keeps its whole letters for the Bioseq of the
/// archive and of the strategy's request (`FastaRecord::cut`, `archive_letters`). Afterwards the records have no source again, as
/// every other run sees them.
///
/// The handles of `register_inputs` are the inputs of the run of their argv alone, both of them,
/// and the records are as its reader gave them; a run that reads the ids of the deflines
/// (`run::inputs_together`) takes no other: its report would be that of a search without the option.
/// A run with NCBI's data loaders on (`data_loaders`, `run::data_loaders` of the argv) takes no
/// records read with them off, which may hold a first line that its reader refuses as a Seq-id
/// (line WQA's D13 and D16, 2026-11-20); `None` where the argv's start fails, which fails the run.
/// `work` takes the query's end with its records (`Registered::end`).
///
/// A run of `-db` (`run::searches_database`) takes the subject handle 0: the engine reads the
/// database, and `work` gets no subject records (line WQA D12a, 2026-11-20).
pub fn with_inputs<R>(
    program: Program,
    query: u32,
    subject: u32,
    words: &[&str],
    data_loaders: Option<bool>,
    work: impl FnOnce(&[FastaRecord], &QueryEnd, &[FastaRecord]) -> R,
) -> Result<R, String> {
    let database = crate::run::searches_database(words);
    if database && subject != 0 {
        return Err(format!(
            "a run of -db takes no subject handle: give 0, not {subject}"
        ));
    }
    let mut store = store().lock().expect("input store");
    let roles = [
        (query, ROLE_QUERY, "query"),
        (subject, ROLE_SUBJECT, "subject"),
    ];
    // A run of -db has no subject (handle 0, checked above); its query is checked as any query.
    for &(handle, role, name) in roles
        .iter()
        .filter(|(_, role, _)| !(database && *role == ROLE_SUBJECT))
    {
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
        if data_loaders == Some(true) && !entry.data_loaders {
            return Err(format!(
                "{name} handle {handle} was read with the data loaders off (DATA_LOADERS=none), which this run has on"
            ));
        }
    }
    let read_for = |handle: u32| store.entries[&handle].run.as_deref();
    // A run of -db has the query's reading only: as its subject, the query's.
    let subject_read = if database {
        read_for(query)
    } else {
        read_for(subject)
    };
    let read_together = match (read_for(query), subject_read) {
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
    let whole = !read_together && crate::run::keeps_whole_records(words);
    let handles: &[u32] = if database {
        &[query]
    } else {
        &[query, subject]
    };
    if whole {
        set_sources(&mut store, handles, true);
    }
    let subjects: &[FastaRecord] = if database {
        &[]
    } else {
        &store.entries[&subject].records
    };
    let result = work(
        &store.entries[&query].records,
        &store.entries[&query].end,
        subjects,
    );
    if whole {
        set_sources(&mut store, handles, false);
    }
    Ok(result)
}

/// Gives every record of `handles` a source with its kept defline (`give`), or takes the sources
/// back.
fn set_sources(store: &mut Store, handles: &[u32], give: bool) {
    for &handle in handles {
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

    /// NCBI's text in LOSAT's words (`losat::web_adapter::losat_text`): the tests state what NCBI
    /// writes and compare it with what LOSAT writes.
    fn ncbi_in_losat_words(ncbi: &str) -> String {
        losat::web_adapter::losat_text(ncbi).into_owned()
    }

    // `register` stops when the scan and the reader disagree on the records.
    #[test]
    fn a_scan_that_disagrees_with_the_reader_fails_register() {
        let bytes = b">a one\nACGT\n>b\nGG\n";
        let records = read_all(
            &mut FastaInputSource::from_bytes(bytes, ReaderConfig::query("BLASTN", false, true)),
            &mut |_| Ok(()),
        )
        .unwrap();
        let scanned =
            scan::scan_ncbi(bytes, reader_config(Program::Blastn, ROLE_QUERY, true)).unwrap();
        let input = QueryEnd::Input;
        assert!(check_scan("query", false, &scanned, &records, &input).is_ok());

        let mut fewer = scanned.clone();
        fewer.records.truncate(1);
        assert!(check_scan("query", false, &fewer, &records, &input).is_err());
        let mut renamed = scanned.clone();
        renamed.records[1].id = "c".to_string();
        assert!(check_scan("query", false, &renamed, &records, &input).is_err());
        let mut longer = scanned.clone();
        longer.records[0].length += 1;
        assert!(check_scan("query", false, &longer, &records, &input).is_err());
        let mut recounted = scanned.clone();
        recounted.records[0].residue_counts[b'A' as usize] += 1;
        assert!(check_scan("query", false, &recounted, &records, &input).is_err());

        // Where the reader's error ends the records: the records before it, and the scan's stop
        // where the scan meets that error too (line WQA's D10, 2026-11-19).
        let bad = b">a one\nACGT\n>b\nGG\n>c\nACGT%%%%\n";
        let (read, end) = read_query_records(bad, ReaderConfig::query("BLASTN", false, true));
        assert!(matches!(end, QueryEnd::Error(_)));
        let stopped =
            scan::scan_ncbi(bad, reader_config(Program::Blastn, ROLE_QUERY, true)).unwrap();
        assert_eq!(stopped.records, scanned.records);
        assert_eq!(stopped.stop.as_ref().map(|stop| stop.offset), Some(18));
        assert!(check_scan("query", false, &stopped, &read, &end).is_ok());
        assert!(check_scan("query", false, &stopped, &read, &input).is_err());
        assert!(check_scan("query", false, &scanned, &read, &input).is_ok());
        // An error that the scan does not meet: the scan has more records, never fewer.
        assert!(check_scan("query", false, &scanned, &read[..1], &end).is_ok());
        assert!(check_scan("query", false, &fewer, &read, &end).is_err());
        let mut stopped_early = stopped.clone();
        stopped_early.records.truncate(1);
        assert!(check_scan("query", false, &stopped_early, &read[..1], &end).is_ok());
        assert!(check_scan("query", false, &stopped_early, &read, &end).is_err());
        // The reader's error before the scan's own stop (an id under `-parse_deflines`): the scan
        // has more records (review M1, 2026-11-20).
        assert!(check_scan("query", false, &stopped, &read[..1], &end).is_ok());
    }

    /// The `records` of a *register* response.
    fn records(program: &str, role: u32, bytes: &[u8]) -> String {
        let (_, response) = register(program, role, true, bytes).unwrap();
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
        // no record, and the run's `Empty CBlastQueryVector` (a subject without records; a
        // query whose first batch has none, `QueryEnd::BlankLines`, line WQA's D10, 2026-11-19).
        assert_eq!(
            records("blastn", ROLE_QUERY, b" \n\t\n"),
            r#""records":[]}"#
        );
        assert_eq!(
            records("blastp", ROLE_SUBJECT, b";c\n\n"),
            r#""records":[]}"#
        );
        assert_eq!(records("blastp", ROLE_QUERY, b";c\n\n"), r#""records":[]}"#);
    }

    // A `>?` gap line and the lines that NCBI's line reader joins are read as the CLI reads them,
    // where the scan refused them (line WQA's D13 (a) and (b), 2026-11-20): the gap's letters are
    // the record's (N for nucleotide input, X for protein input).
    #[test]
    fn register_reads_gap_lines_and_joined_lines() {
        assert_eq!(
            records("blastn", ROLE_SUBJECT, b">s\nACGT\n>?10\nACGT\n"),
            r#""records":[{"index":0,"id":"s","length":18}]}"#
        );
        assert_eq!(
            records("blastp", ROLE_QUERY, b">?5\nMKT\n>p\n>?unk3\nAC\n"),
            r#""records":[{"index":0,"id":"","length":8},{"index":1,"id":"p","length":5}]}"#
        );
        // The fixtures of NCBI's line reader whose lines it joins (`LOSAT/tests/fixtures/fasta_input`,
        // which the CLI reads as NCBI BLAST+ 2.17.0 does): `register` checks the scan's records
        // against the reader's.
        let fixtures = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../LOSAT/tests/fixtures/fasta_input"
        );
        for (name, program, role) in [
            ("eol_lfcr.q.nuc.fa", "blastn", ROLE_QUERY),
            ("eol_lfcr_defline.s.prot.faa", "blastp", ROLE_SUBJECT),
            ("eol_mixed_lf_cr.q.nuc.fa", "tblastx", ROLE_QUERY),
            ("eol_mixed_crlf_cr.s.prot.faa", "blastx", ROLE_SUBJECT),
            ("eol_mixed_lf_crlf_cr.q.prot.faa", "tblastn", ROLE_QUERY),
        ] {
            let bytes = std::fs::read(format!("{fixtures}/{name}")).unwrap();
            let (handle, _) = register(program, role, true, &bytes)
                .unwrap_or_else(|error| panic!("{name}: {error}"));
            release(handle).unwrap();
        }
        // A CR line end that NCBI's line reader joins to the `!` line after it, in a file of LF line ends,
        // within a record of six lines, which the scan refused ("NCBI BLAST+ reads two of its lines as one");
        // the CLI drops the joined line's residues (210 of 240 letters).
        assert_eq!(
            records(
                "blastn",
                ROLE_QUERY,
                b">q\nAAAAAAAAAA\nCCCCCCCCCC\r!x\nGGGGGGGGGG\nTTTTTTTTTT\nACGTA\n"
            ),
            records(
                "blastn",
                ROLE_QUERY,
                b">q\nAAAAAAAAAA\nCCCCCCCCCC\nTTTTTTTTTT\nACGTA\n"
            ),
        );
        let joined = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/d13/joined.nuc.fa"
        );
        assert_eq!(
            records("blastn", ROLE_SUBJECT, &std::fs::read(joined).unwrap()),
            r#""records":[{"index":0,"id":"q1","length":210}]}"#
        );
    }

    // The reader's error with the CLI's text, and its refusals in the CLI's words, not LOSAT Web's
    // (line WQA's D13 (c), 2026-11-20); with the data loaders off (`DATA_LOADERS=none`) a first line
    // that may be a Seq-id is data, as the CLI reads it (D16).
    #[test]
    fn register_refuses_what_the_cli_refuses_in_its_words() {
        let seq_id = register("tblastx", ROLE_QUERY, true, b"AB123456\nACGT\n").unwrap_err();
        assert_eq!(
            seq_id,
            "the first line of the query (\"AB123456\") is not a defline and may be a sequence identifier that NCBI BLAST+ fetches through a data loader (from GenBank or a BLAST database), which is not supported by LOSAT's TBLASTX (start the query with a '>' defline, or set DATA_LOADERS=none in the [BLAST] section of .ncbirc)"
        );
        let subject = register("blastp", ROLE_SUBJECT, true, b"2HHB\nMKTAYIAKQR\n").unwrap_err();
        assert!(
            subject.contains("the first line of the subject (\"2HHB\")")
                && subject.contains("LOSAT's BLASTP"),
            "{subject}"
        );
        let (handle, response) =
            register("blastp", ROLE_SUBJECT, false, b"2HHB\nMKTAYIAKQR\n").unwrap();
        assert!(
            response.ends_with(r#""records":[{"index":0,"id":"","length":13}]}"#),
            "{response}"
        );
        release(handle).unwrap();
        // A subject's reader error fails `register`; a query's ends its records, and a run meets it
        // in its batch (line WQA's D10, 2026-11-19): the response says it.
        let cli = ncbi_in_losat_words("BLAST query error: CFastaReader: Near line 4, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.");
        let bad = b">p\nMKV\n>q\n1234567890123456789012345\n>r\nMKV\n";
        assert_eq!(
            register("blastp", ROLE_SUBJECT, true, bad).unwrap_err(),
            cli
        );
        let (handle, response) = register("blastp", ROLE_QUERY, true, bad).unwrap();
        assert_eq!(
            &response[response.find("\"records\":").unwrap()..],
            format!(r#""records":[{{"index":0,"id":"p","length":3}}],"reader_error":"{cli}"}}"#)
        );
        release(handle).unwrap();
        // With the loaders off the reader's check of the data line ends the records, as the CLI's:
        // a subject's fails `register`, a query's is its first record's error.
        let data = ncbi_in_losat_words("BLAST query error: CFastaReader: Near line 1, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.");
        assert_eq!(
            register("blastn", ROLE_SUBJECT, false, b"NM_000518.5\n").unwrap_err(),
            data
        );
        let (handle, response) = register("blastn", ROLE_QUERY, false, b"NM_000518.5\n").unwrap();
        assert!(
            response.ends_with(&format!(r#""records":[],"reader_error":"{data}"}}"#)),
            "{response}"
        );
        release(handle).unwrap();
    }

    // A run reads its inputs with the data loaders that its registry sets (`-conffile` with
    // `[BLAST] DATA_LOADERS=none`; line WQA's D16, 2026-11-20): `register_inputs` reads such a first
    // line as data, and a run with them on takes no records read with them off.
    #[test]
    fn a_run_reads_with_the_data_loaders_of_its_registry() {
        let registry = b"[BLAST]\nDATA_LOADERS=none\n";
        let file =
            crate::files::begin(crate::files::KIND_REGISTRY, "d16-store.ini", registry.len())
                .unwrap();
        crate::files::chunk(file, registry).unwrap();
        crate::files::end(file).unwrap();
        let off = [
            "blastp",
            "-conffile",
            "d16-store.ini",
            "-query",
            "q",
            "-subject",
            "s",
        ];
        let on = [&off[..1], &off[3..]].concat();
        let (query, subject) = (&b"2HHB\nMKTAYIAKQR\n"[..], &b">s\nMKTAYIAKQR\n"[..]);
        assert_eq!(crate::run::data_loaders(&off), Ok(false));
        assert_eq!(crate::run::data_loaders(&on), Ok(true));
        let (query_handle, subject_handle, response) =
            register_inputs(&off, query, subject).unwrap();
        assert!(response.contains(r#""length":13"#), "{response}");
        let error = register_inputs(&on, query, subject).unwrap_err();
        assert!(
            error.contains("the first line of the query (\"2HHB\")"),
            "{error}"
        );
        release(query_handle).unwrap();
        release(subject_handle).unwrap();
        let (query_handle, _) = register("blastp", ROLE_QUERY, false, query).unwrap();
        let (subject_handle, _) = register("blastp", ROLE_SUBJECT, false, subject).unwrap();
        let plain = ["blastp", "-query", "q", "-subject", "s"];
        let refused = with_inputs(
            Program::Blastp,
            query_handle,
            subject_handle,
            &plain,
            Some(true),
            |_, _, _| (),
        )
        .unwrap_err();
        assert!(
            refused.contains("read with the data loaders off"),
            "{refused}"
        );
        assert!(with_inputs(
            Program::Blastp,
            query_handle,
            subject_handle,
            &plain,
            Some(false),
            |_, _, _| ()
        )
        .is_ok());
        release(query_handle).unwrap();
        release(subject_handle).unwrap();
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
        let ids = |queries: &[FastaRecord], _: &QueryEnd, subjects: &[FastaRecord]| {
            (
                queries[0].believed().is_some(),
                subjects[0].believed().is_some(),
            )
        };
        assert_eq!(
            with_inputs(Program::Blastp, query, subject, &words, Some(true), ids),
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
        let error = with_inputs(
            Program::Blastp,
            query,
            subject,
            &other,
            Some(true),
            |_, _, _| (),
        )
        .unwrap_err();
        assert!(
            error.contains("are not the inputs that register_inputs read for this argv"),
            "{error}"
        );
        let (plain_query, _) = register("blastp", ROLE_QUERY, true, query_bytes).unwrap();
        let (plain_subject, _) = register("blastp", ROLE_SUBJECT, true, subject_bytes).unwrap();
        let mixed = with_inputs(
            Program::Blastp,
            plain_query,
            subject,
            &words,
            Some(true),
            |_, _, _| (),
        )
        .unwrap_err();
        assert!(
            mixed.contains("are not the inputs that register_inputs read for this argv"),
            "{mixed}"
        );
        assert_eq!(
            with_inputs(Program::Blastp, plain_query, plain_subject, &words, Some(true), |_, _, _| ()).unwrap_err(),
            "-parse_deflines needs the inputs of register_inputs, which reads the ids of the deflines"
        );
        // Without the option the handles of `register` are the inputs, as before.
        let plain = &words[..5];
        assert!(with_inputs(
            Program::Blastp,
            plain_query,
            plain_subject,
            plain,
            Some(true),
            |_, _, _| ()
        )
        .is_ok());
        for handle in [query, subject, plain_query, plain_subject] {
            release(handle).unwrap();
        }
    }

    // Under `-parse_deflines` the reader's error can come before the scan's own stop (review M1,
    // 2026-11-20): a refused database-ordinal Seq-id (`blastp -parse_deflines`: "Error: the query
    // Seq-id ... names a database ordinal" on the command line), or a first record without a
    // defline, before a line that the scan stops at too. The inputs register with the reader's
    // error, as the command line ends there.
    #[test]
    fn the_readers_error_can_come_before_the_scans_stop() {
        let words = ["blastp", "-query", "q", "-subject", "s", "-parse_deflines"];
        let subject = b">s1\nMKVLAAGIVA\n";
        let cases: [(&[u8], &str, usize); 2] = [
            (
                b">p1\nMKVLAAGIVA\n>gnl|BL_ORD_ID|1\nMKVLAAGIVA\n>bad\nMKVL%%%%\n",
                "names a database ordinal",
                1,
            ),
            (
                b"MKVLAAGIVA\n>p2\nMKVLAAGIVA\n>bad\nMKVL%%%%\n",
                "doesn't start with a defline",
                0,
            ),
        ];
        for (query, message, read) in cases {
            let scanned =
                scan::scan_ncbi(query, reader_config(Program::Blastp, ROLE_QUERY, true)).unwrap();
            assert_eq!(scanned.records.len(), 2);
            assert!(scanned.stop.is_some());
            let (query_handle, subject_handle, response) =
                register_inputs(&words, query, subject).unwrap();
            assert!(response.contains(message), "{response}");
            let read_records = |queries: &[FastaRecord], end: &QueryEnd, _: &[FastaRecord]| {
                (queries.len(), matches!(end, QueryEnd::Error(_)))
            };
            assert_eq!(
                with_inputs(
                    Program::Blastp,
                    query_handle,
                    subject_handle,
                    &words,
                    Some(true),
                    read_records
                ),
                Ok((read, true))
            );
            for handle in [query_handle, subject_handle] {
                release(handle).unwrap();
            }
        }
    }

    // Under `-parse_deflines`, a Seq-id that names the other molecule is read with that molecule's
    // alphabet (line RG-XT item X7-5): `check_scan` accepts the resulting length and counts, and
    // the response reports the reader's length.
    #[test]
    fn other_molecule_parse_deflines_registers_with_readers_length() {
        let bp_words = ["blastp", "-query", "q", "-subject", "s", "-parse_deflines"];
        let (qh, sh, resp) = register_inputs(
            &bp_words,
            b">p1\nMKVLAAGIVA\n",
            b">gnl|TI|12345 trace\nMKVLEFAAGIVA\n",
        )
        .unwrap();
        // M, K, V, A, A, G, V, A kept (8); L, E, F, I dropped by nucleotide alphabet.
        assert!(resp.contains(r#""length":8"#), "blastp: {resp}");
        release(qh).unwrap();
        release(sh).unwrap();

        let bn_words = ["blastn", "-query", "q", "-subject", "s", "-parse_deflines"];
        let (qh, sh, resp) = register_inputs(
            &bn_words,
            b">n1\nACGTACGT\n",
            b">sp|P99995.1|NAME_W\nACGEFLPQ\n",
        )
        .unwrap();
        // All 8 letters kept by protein alphabet; E, F, L, P, Q become N.
        assert!(resp.contains(r#""length":8"#), "blastn: {resp}");
        release(qh).unwrap();
        release(sh).unwrap();

        // Protein-only letters on the first data line of a protein-named nucleotide record fill
        // `m_SeqData`, so a later non-plausible line skips `CheckDataLine` (NCBI exit 0).
        let (qh, sh, resp) = register_inputs(
            &bn_words,
            b">sp|P99995.1|NAME_W\nEEEEEE\n%%%%%%\n",
            b">sp|P99996.1|NAME_S\nEEEEEE\n%%%%%%\n",
        )
        .unwrap();
        assert!(resp.contains(r#""length":6"#), "{resp}");
        assert!(!resp.contains("reader_error"), "{resp}");
        release(qh).unwrap();
        release(sh).unwrap();

        // Conversely, a nucleotide-named protein record drops `EEEEEE` (`m_SeqData` stays empty),
        // so `%%%%%%` triggers `CheckDataLine` (subject fails registration; query records error).
        let bad_ti = b">gnl|TI|12345 trace\nEEEEEE\n%%%%%%\n";
        assert!(register_inputs(&bp_words, b">p1\nMKVLAAGIVA\n", bad_ti).is_err());
        let (qh, sh, resp) = register_inputs(&bp_words, bad_ti, b">p1\nMKVLAAGIVA\n").unwrap();
        assert!(resp.contains("doesn't look like plausible data"), "{resp}");
        release(qh).unwrap();
        release(sh).unwrap();
    }

    // Another program's handle cannot bring a record to a run.
    #[test]
    fn a_handle_belongs_to_its_program() {
        let (query, _) = register("blastp", ROLE_QUERY, true, b">q t\nACXGT\n").unwrap();
        let (subject, _) = register("blastn", ROLE_SUBJECT, true, b">s\nACGT\n").unwrap();
        let words = ["blastn", "-query", "q", "-subject", "s"];
        let error = with_inputs(
            Program::Blastn,
            query,
            subject,
            &words,
            Some(true),
            |_, _, _| (),
        )
        .unwrap_err();
        assert!(
            error.contains("registered for blastp, not blastn"),
            "{error}"
        );
        let words = ["blastp", "-query", "q", "-subject", "s"];
        assert!(with_inputs(
            Program::Blastp,
            query,
            subject,
            &words,
            Some(true),
            |_, _, _| ()
        )
        .is_err());
    }
}
