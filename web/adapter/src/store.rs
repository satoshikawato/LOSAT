//! Registered inputs (docs/web/abi_v2.md §4 `register`, plan §4.6 R1): the records
//! read by the program's own FASTA reader (the engine's port of NCBI BLAST+'s
//! `CBlastFastaInputSource`), kept under a handle so that one subject serves several runs.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::sync::Mutex;

use LOSAT::blastinput::fasta_reader::{
    read_all, FastaInputSource, FastaRecord, ReadError, ReaderConfig,
};

use crate::json;
use crate::run::Program;
use crate::scan;

pub const ROLE_QUERY: u32 = 0;
pub const ROLE_SUBJECT: u32 = 1;

pub struct Registered {
    pub role: u32,
    /// The program whose reader read the records; only its runs use them.
    pub program: Program,
    pub records: Vec<FastaRecord>,
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
/// ```c++
///     flags += (m_ReadProteins
///               ? CFastaReader::fAssumeProt
///               : CFastaReader::fAssumeNuc);
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:51,106
/// ```c++
///     const bool kQueryIsProtein = true;
/// ...
///     m_QueryOptsArgs.Reset(new CQueryOptionsArgs(kQueryIsProtein));
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2554
/// ```c++
///         m_Scope = ReadSequencesToBlast(*subj_input_stream, IsProtein(),
/// ```
/// BLASTN and TBLASTX read nucleotides in both roles, BLASTP proteins in both, TBLASTN a
/// protein query and nucleotide subjects.
fn reader_of(program: Program, role: u32) -> (&'static str, bool) {
    match (program, role) {
        (Program::Blastn, _) => ("BLASTN", false),
        (Program::Tblastx, _) => ("TBLASTX", false),
        (Program::Tblastn, ROLE_SUBJECT) => ("TBLASTN", false),
        (Program::Tblastn, _) => ("TBLASTN", true),
        (Program::Blastp, _) => ("BLASTP", true),
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
/// ```c++
///     static const string kDataLoadersConfig("DATA_LOADERS");
///
///     if (registry.HasEntry("BLAST", kDataLoadersConfig)) {
/// ```
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
    // ```c++
    //         if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())){
    //            	ERR_POST(Warning << "Query is Empty!");
    //            	return BLAST_EXIT_SUCCESS;
    //         }
    //         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
    //         CBlastInput input(&fasta, m_CmdLineArgs->GetQueryBatchSize());
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:146-152
    // ```c++
    //         try { q.Reset(m_Source->GetNextSequence(scope)); }
    //         catch (const CObjReaderParseException& e) {
    //             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
    //                 break;
    //             }
    //             throw;
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
    // ```c++
    // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
    //     : m_QueryVector(& queries)
    // {
    //     if (queries.Empty()) {
    //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
    //     }
    // ```
    // A query of white space only is NCBI's `Query is Empty!` at run time; a query with
    // blank and comment lines only is not empty, and its first batch (`eEOF`, no query)
    // fails every run, so `register` fails with NCBI's error (port plan Q2).
    if role == ROLE_QUERY && records.is_empty() && !LOSAT::blastinput::input_files::is_blank(bytes)
    {
        return Err(format!(
            "{:#}",
            LOSAT::blastinput::app::engine_error("Empty CBlastQueryVector")
        )
        .trim_end_matches('\n')
        .to_string());
    }
    check_scan(role_name, &scanned, &records)?;
    let mut response = String::from("{\"handle\":");
    let mut store = store().lock().expect("input store");
    store.next = store
        .next
        .checked_add(1)
        .ok_or_else(|| "input handle space exhausted".to_string())?;
    let handle = store.next;
    let _ = write!(response, "{handle},\"records\":[");
    for (index, record) in records.iter().enumerate() {
        if index > 0 {
            response.push(',');
        }
        let _ = write!(response, "{{\"index\":{index},\"id\":");
        json::string(&mut response, &record_id(record.title()));
        let _ = write!(response, ",\"length\":{}}}", record.seq().len());
    }
    response.push_str("]}");
    store.entries.insert(
        handle,
        Registered {
            role,
            program,
            records,
        },
    );
    Ok((handle, response))
}

/// Fails unless the index scan found the reader's records: the same number, and for each
/// record the same ID, length and residue counts (the stored residues upper-cased; plan
/// TD-8).
fn check_scan(
    role: &str,
    scanned: &[scan::ScanRecord],
    records: &[FastaRecord],
) -> Result<(), String> {
    let agrees = |scan: &scan::ScanRecord, record: &FastaRecord| {
        let mut counts = [0u64; 256];
        for &byte in record.seq() {
            counts[byte.to_ascii_uppercase() as usize] += 1;
        }
        scan.id == record_id(record.title())
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
/// `program`.
pub fn with_inputs<R>(
    program: Program,
    query: u32,
    subject: u32,
    work: impl FnOnce(&[FastaRecord], &[FastaRecord]) -> R,
) -> Result<R, String> {
    let store = store().lock().expect("input store");
    let get = |handle: u32, role: u32, name: &str| {
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
        Ok(entry.records.as_slice())
    };
    let queries = get(query, ROLE_QUERY, "query")?;
    let subjects = get(subject, ROLE_SUBJECT, "subject")?;
    Ok(work(queries, subjects))
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
        assert!(gap.ends_with("not supported by LOSAT Web"), "{gap}");
        let seq_id = register("tblastx", ROLE_QUERY, b"AB123456\nACGT\n").unwrap_err();
        assert!(seq_id.contains("not supported by LOSAT Web"), "{seq_id}");
        let reader =
            register("blastp", ROLE_QUERY, b">q\n1234567890123456789012345\n").unwrap_err();
        assert_eq!(
            reader,
            "BLAST query error: CFastaReader: Near line 2, there's a line that doesn't look like plausible data, but it's not marked as defline or comment."
        );
    }

    // Another program's handle cannot bring a record to a run.
    #[test]
    fn a_handle_belongs_to_its_program() {
        let (query, _) = register("blastp", ROLE_QUERY, b">q t\nACXGT\n").unwrap();
        let (subject, _) = register("blastn", ROLE_SUBJECT, b">s\nACGT\n").unwrap();
        let error = with_inputs(Program::Blastn, query, subject, |_, _| ()).unwrap_err();
        assert!(
            error.contains("registered for blastp, not blastn"),
            "{error}"
        );
        assert!(with_inputs(Program::Blastp, query, subject, |_, _| ()).is_err());
    }
}
