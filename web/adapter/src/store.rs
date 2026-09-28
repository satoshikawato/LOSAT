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
    Program::parse(program)?;
    let role_name = match role {
        ROLE_QUERY => "query",
        ROLE_SUBJECT => "subject",
        other => return Err(format!("unknown input role {other}")),
    };
    // BLASTP, TBLASTN, BLASTN and TBLASTX read their inputs with bio::io::fasta.
    let records = fasta::Reader::new(bytes)
        .records()
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| format!("failed to read {role_name} FASTA: {error}"))?;
    let scanned = scan::scan(bytes).map_err(|error| {
        format!("the index scan of the {role_name} FASTA disagrees with the parser: {error}")
    })?;
    if scanned.len() != records.len()
        || scanned.iter().zip(&records).any(|(scan, record)| {
            scan.id != record.id() || scan.length != record.seq().len() as u64
        })
    {
        return Err(format!(
            "the index scan of the {role_name} FASTA disagrees with the parser"
        ));
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
    store.entries.insert(handle, Registered { role, records });
    Ok((handle, response))
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

/// Calls `work` with the records of a query handle and a subject handle.
pub fn with_inputs<R>(
    query: u32,
    subject: u32,
    work: impl FnOnce(&[fasta::Record], &[fasta::Record]) -> R,
) -> Result<R, String> {
    let store = store().lock().expect("input store");
    let get = |handle: u32, role: u32, name: &str| {
        store
            .entries
            .get(&handle)
            .filter(|entry| entry.role == role)
            .map(|entry| entry.records.as_slice())
            .ok_or_else(|| format!("{handle} is not a registered {name} handle"))
    };
    let queries = get(query, ROLE_QUERY, "query")?;
    let subjects = get(subject, ROLE_SUBJECT, "subject")?;
    Ok(work(queries, subjects))
}
