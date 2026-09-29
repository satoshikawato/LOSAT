//! Registered inputs (docs/web/abi_v2.md §4 `register`, plan §4.6 R1): the records
//! parsed by the program's own FASTA parser, kept under a handle so that one subject
//! serves several runs.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::sync::Mutex;

use bio::io::fasta;

use crate::json;
use crate::run::Program;
use crate::scan;

pub const ROLE_QUERY: u32 = 0;
pub const ROLE_SUBJECT: u32 = 1;

pub struct Registered {
    pub role: u32,
    /// The program whose parser and input checks read the records; only its runs use them.
    pub program: Program,
    pub records: Vec<fasta::Record>,
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

/// Parses `bytes` with the program's parser, checks the records against the index
/// scan of the same bytes (plan TD-8), keeps them, and returns the handle and the
/// *register* JSON.
pub fn register(program: &str, role: u32, bytes: &[u8]) -> Result<(u32, String), String> {
    let program = Program::parse(program)?;
    let role_name = match role {
        ROLE_QUERY => "query",
        ROLE_SUBJECT => "subject",
        other => return Err(format!("unknown input role {other}")),
    };
    use LOSAT::algorithm::blastn::input as blastn_input;
    let blastn = program == Program::Blastn;
    // BLASTN, as the CLI: a file of white space only has no record (NCBI's empty query, or
    // its error for no subject at run time).
    let records = if blastn && blastn_input::is_blank(bytes) {
        Vec::new()
    } else {
        // BLASTP, TBLASTN, BLASTN and TBLASTX read their inputs with bio::io::fasta.
        let records = fasta::Reader::new(bytes)
            .records()
            .collect::<Result<Vec<_>, _>>()
            .map_err(|error| {
                let unsupported = if blastn {
                    format!(" ({})", blastn_input::UNREADABLE_FASTA)
                } else {
                    String::new()
                };
                format!("failed to read {role_name} FASTA: {error}{unsupported}")
            })?;
        let scanned = scan::scan(bytes).map_err(|error| {
            format!("the index scan of the {role_name} FASTA disagrees with the parser: {error}")
        })?;
        check_scan(role_name, &scanned, &records)?;
        records
    };
    // BLASTN rejects the inputs that NCBI BLAST+ reads differently from bio.
    if blastn {
        blastn_input::check_deflines(bytes, role_name)
            .and_then(|()| blastn_input::check_residues(&records, role_name))
            .and_then(|()| blastn_input::check_records_have_residues(&records, role_name))
            .map_err(|error| format!("{error:#}"))?;
    }
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
        json::string(&mut response, record.id());
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

/// Fails unless the index scan found the parser's records: the same number, and for
/// each record the same ID, length and residue counts (plan TD-8).
fn check_scan(
    role: &str,
    scanned: &[scan::ScanRecord],
    records: &[fasta::Record],
) -> Result<(), String> {
    let agrees = |scan: &scan::ScanRecord, record: &fasta::Record| {
        let mut counts = [0u64; 256];
        for &byte in record.seq() {
            counts[byte as usize] += 1;
        }
        scan.id == record.id()
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
            "the index scan of the {role} FASTA disagrees with the parser"
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
    work: impl FnOnce(&[fasta::Record], &[fasta::Record]) -> R,
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

    // `register` stops when the scan and the parser disagree on the records.
    #[test]
    fn a_scan_that_disagrees_with_the_parser_fails_register() {
        let bytes = b">a one\nACGT\n>b\nGG\n";
        let records: Vec<fasta::Record> = fasta::Reader::new(&bytes[..])
            .records()
            .collect::<Result<_, _>>()
            .unwrap();
        let scanned = scan::scan(bytes).unwrap();
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

    // BLASTN refuses a record that NCBI BLAST+ reads differently; the other programs
    // keep their readers.
    #[test]
    fn blastn_register_rejects_records_that_ncbi_reads_differently() {
        let bytes = b">q\tt\nACXGT\n";
        let error = register("blastn", ROLE_QUERY, bytes).unwrap_err();
        assert!(error.contains("not supported by LOSAT's BLASTN"), "{error}");
        assert!(register("blastn", ROLE_QUERY, b">q t\nACGUT\n").is_ok());
        // Another program's handle cannot bring the record to a BLASTN run.
        let (query, _) = register("blastp", ROLE_QUERY, bytes).unwrap();
        let (subject, _) = register("blastn", ROLE_SUBJECT, b">s\nACGT\n").unwrap();
        let error = with_inputs(Program::Blastn, query, subject, |_, _| ()).unwrap_err();
        assert!(
            error.contains("registered for blastp, not blastn"),
            "{error}"
        );
        // FASTA that bio cannot read (text before the first defline, bytes that are not
        // UTF-8) says that LOSAT does not support it; white space only is a file without
        // records.
        let error = register("blastn", ROLE_SUBJECT, b">s0\n>s1\nACGT\n").unwrap_err();
        assert!(
            error.contains("subject record 1 (s0) has no residues"),
            "{error}"
        );
        for bytes in [&b"\n>q\nACGT\n"[..], b">q\nAC\xffGT\n"] {
            let error = register("blastn", ROLE_QUERY, bytes).unwrap_err();
            assert!(error.contains("not supported by LOSAT's BLASTN"), "{error}");
            assert!(error.contains("bytes that are not UTF-8"), "{error}");
        }
        let (_, response) = register("blastn", ROLE_QUERY, b" \n\t\n").unwrap();
        assert!(response.ends_with("\"records\":[]}"), "{response}");
    }
}
