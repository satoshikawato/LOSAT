//! LOSAT Web ABI v2 (docs/web/abi_v2.md): the exports that the application's engine
//! worker calls. ABI v1 (`losat_web_*`, `LOSAT/src/web_api.rs`) is linked in unchanged
//! and shares no state with these exports (plan TD-1).
//!
//! Every export returns `0` or a handle on success and `-1` on failure, with the
//! message in `losat_web2_last_error_ptr/len`. Responses and outputs go to the host
//! through `losat_host.emit` (`emit.rs`).

pub mod defaults;
pub mod describe;
pub mod emit;
pub mod files;
pub mod inflate;
pub mod json;
pub mod run;
pub mod scan;
pub mod store;
pub mod tool;
pub mod zip;

use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

use emit::{emit, STREAM_RESPONSE};

fn last_error() -> &'static Mutex<Vec<u8>> {
    static LAST_ERROR: OnceLock<Mutex<Vec<u8>>> = OnceLock::new();
    LAST_ERROR.get_or_init(Mutex::default)
}

/// Runs an export body: records its error, or clears the last error on success.
fn guarded(body: impl FnOnce() -> Result<i32, String>) -> i32 {
    run::enter_run_directory();
    let result = body();
    let mut error = last_error().lock().expect("last error");
    match result {
        Ok(value) => {
            error.clear();
            value
        }
        Err(message) => {
            *error = message.into_bytes();
            -1
        }
    }
}

/// The bytes of a host buffer that stays valid for the duration of the call.
///
/// # Safety
/// `ptr` must point to `len` readable bytes, or `len` must be 0.
unsafe fn bytes<'a>(ptr: *const u8, len: usize) -> &'a [u8] {
    if len == 0 {
        &[]
    } else {
        std::slice::from_raw_parts(ptr, len)
    }
}

/// # Safety
/// As `bytes`.
unsafe fn text<'a>(ptr: *const u8, len: usize, name: &str) -> Result<&'a str, String> {
    std::str::from_utf8(bytes(ptr, len)).map_err(|error| format!("{name} is not UTF-8: {error}"))
}

/// The NUL-separated argv words (docs/web/abi_v2.md §7); a trailing NUL is allowed.
fn words(argv: &str) -> Vec<&str> {
    let mut words: Vec<&str> = argv.split('\0').collect();
    if words.last() == Some(&"") {
        words.pop();
    }
    words
}

/// The version of ABI v2 (docs/web/abi_v2.md §10): 14 since `losat_web2_scan_begin` and
/// `losat_web2_register` take `parse_deflines` and *validation* emits `parse_deflines` when the
/// argv parses so the index scan and the input check read records whose `-parse_deflines` Seq-id
/// names the other molecule as the engine's reader reads them (line RG-XT item X7-5, 2026-11-25);
/// 13 since `losat_web2_tool`, which runs `windowmasker` with the command line's words in the
/// run's working directory, streams its standard output on stream 11 and the files that it writes
/// on stream 12 as they are written, and gives its exit status on stream 2, and since
/// `losat_web2_validate` checks a `windowmasker` argv (line NP WW3a, 2026-11-22); 12 since a run
/// of `-db`, which takes no subject handle and reads the database from the run's working
/// directory, which the host preopens as `/work`, and the *start* field `database` (line WQA D12a,
/// 2026-11-20; 9 on its branch until 10 merged first, then 12 after D10's 11); 11 since a query is
/// read as the command line's query batches read it: the records before the reader's error, with
/// the error (`reader_error` of *register*, *inputs* and *scan*), which a run meets in the
/// record's batch (line WQA's D10, 2026-11-19); 10 since a run prints the text that the command
/// line prints in place of a run (its help, the XML of `-xmlhelp`, its versions) on stream 11,
/// and the *start* field `prints` (line WQA's D14, 2026-11-20); 8 since the index scan and
/// `register` read as the engine's reader with the run's data loaders (`scan_begin` and `register`
/// take the program, the role and the setting, which *validation* gives on stream 2; a `>?` gap
/// line and the lines that NCBI's line reader joins are indexed, and the reader's refusals come in
/// its words; line WQA's D13 and D16, 2026-11-20); 7 since a run that fails keeps the output of
/// the argv up to its error and drops the other outputs with the signal of stream 9 (line WQA's
/// sweep, D15, 2026-11-20; versions 5 and 6 are taken by other open branches); 4 since
/// `losat_web2_register_inputs`, which reads the inputs of a run with `-parse_deflines`, and the
/// signal of stream 9 for every output stream that failed (line WQA's sweep, D7, 2026-11-19); 3
/// since the signal of stream 9, which makes the host discard a stream that failed alone after a
/// chunk of it was sent (line RG item T7-W).
#[no_mangle]
pub extern "C" fn losat_web2_abi_version() -> i32 {
    14
}

fn layout(len: usize) -> Option<std::alloc::Layout> {
    std::alloc::Layout::array::<u8>(len.max(1)).ok()
}

/// Allocates `len` bytes for an input buffer, or returns null.
#[no_mangle]
pub extern "C" fn losat_web2_alloc(len: usize) -> *mut u8 {
    match layout(len) {
        // SAFETY: the layout has a non-zero size.
        Some(layout) => unsafe { std::alloc::alloc(layout) },
        None => std::ptr::null_mut(),
    }
}

/// # Safety
/// `ptr` and `len` must come from one `losat_web2_alloc` call.
#[no_mangle]
pub unsafe extern "C" fn losat_web2_dealloc(ptr: *mut u8, len: usize) {
    if let (false, Some(layout)) = (ptr.is_null(), layout(len)) {
        std::alloc::dealloc(ptr, layout);
    }
}

#[no_mangle]
pub extern "C" fn losat_web2_last_error_ptr() -> *const u8 {
    last_error().lock().expect("last error").as_ptr()
}

#[no_mangle]
pub extern "C" fn losat_web2_last_error_len() -> usize {
    last_error().lock().expect("last error").len()
}

/// # Safety
/// The pointer arguments must describe readable buffers of the given lengths.
#[no_mangle]
pub unsafe extern "C" fn losat_web2_describe(program_ptr: *const u8, program_len: usize) -> i32 {
    guarded(|| {
        let response = describe::describe(text(program_ptr, program_len, "program")?)?;
        emit(STREAM_RESPONSE, response.as_bytes());
        Ok(0)
    })
}

/// # Safety
/// The pointer arguments must describe readable buffers of the given lengths.
#[no_mangle]
pub unsafe extern "C" fn losat_web2_start(argv_ptr: *const u8, argv_len: usize) -> i32 {
    guarded(|| {
        run::emit_start(&words(text(argv_ptr, argv_len, "argv")?))?;
        Ok(0)
    })
}

/// Checks the argv as the CLI does before it reads its inputs (`run::validate`; for a tool's argv,
/// `windowmasker`, `tool::validate`, which emits nothing) and emits the *validation* JSON on stream
/// 2: whether the run reads its inputs with NCBI's data loaders on and, when the argv parses,
/// whether `-parse_deflines` is active (`{"data_loaders":true,"parse_deflines":false}`, or
/// `{"data_loaders":true}` when the argv does not parse; docs/web/abi_v2.md §9). It is emitted
/// once the start has read the argv's registry, before the option checks: a refused argv has it
/// too, unless its start or LOSAT's settings failed (line WQA's D13 round 3, L1, 2026-11-20).
///
/// # Safety
/// The pointer arguments must describe readable buffers of the given lengths.
#[no_mangle]
pub unsafe extern "C" fn losat_web2_validate(argv_ptr: *const u8, argv_len: usize) -> i32 {
    guarded(|| {
        let argv = words(text(argv_ptr, argv_len, "argv")?);
        // A tool's argv: the command line's checks before the tool opens a file, no *validation*
        // (line NP WW3a, D-WW3-7).
        if argv.first().is_some_and(|name| tool::TOOLS.contains(name)) {
            tool::validate(&argv)?;
            return Ok(0);
        }
        // The *validation* goes out once the start knows the settings, before the option checks, so
        // that a host has them also when they refuse the argv (line WQA's D13 round 3, L1).
        run::validate_with(&argv, &mut |validation| {
            let response = match validation.parse_deflines {
                Some(parse_deflines) => format!(
                    "{{\"data_loaders\":{},\"parse_deflines\":{parse_deflines}}}",
                    validation.data_loaders
                ),
                None => format!("{{\"data_loaders\":{}}}", validation.data_loaders),
            };
            emit(STREAM_RESPONSE, response.as_bytes());
        })?;
        Ok(0)
    })
}

/// A flag argument of the ABI: `0` or `1`.
fn flag(value: u32, name: &str) -> Result<bool, String> {
    match value {
        0 => Ok(false),
        1 => Ok(true),
        other => Err(format!("{name} must be 0 or 1, not {other}")),
    }
}

/// Reads the FASTA `bytes` of `program`'s input of `role` (`0` query, `1` subject) with NCBI's
/// data loaders on (`data_loaders` 1) or off (0) and `-parse_deflines` on (`parse_deflines` 1) or
/// off (0), as the run's *validation* gives them (`store::register`); emits the *register* JSON on
/// stream 2 and returns the handle.
///
/// # Safety
/// The pointer arguments must describe readable buffers of the given lengths.
#[no_mangle]
pub unsafe extern "C" fn losat_web2_register(
    program_ptr: *const u8,
    program_len: usize,
    role: u32,
    data_loaders: u32,
    parse_deflines: u32,
    bytes_ptr: *const u8,
    bytes_len: usize,
) -> i32 {
    guarded(|| {
        let program = text(program_ptr, program_len, "program")?;
        let data_loaders = flag(data_loaders, "data_loaders")?;
        let parse_deflines = flag(parse_deflines, "parse_deflines")?;
        let (handle, response) = store::register(
            program,
            role,
            data_loaders,
            parse_deflines,
            bytes(bytes_ptr, bytes_len),
        )?;
        emit(STREAM_RESPONSE, response.as_bytes());
        i32::try_from(handle).map_err(|_| "input handle space exhausted".to_string())
    })
}

/// Reads the query and the subject FASTA of the run of `argv` together, as the CLI reads its
/// `-query` and `-subject` files for that argv (`store::register_inputs`): the inputs of a run with
/// `-parse_deflines`, whose ids are read in one scope for the query and the subjects. Emits the
/// *inputs* JSON (both handles and their records) on stream 2; only the run of `argv` takes the
/// handles.
///
/// # Safety
/// The pointer arguments must describe readable buffers of the given lengths.
#[no_mangle]
pub unsafe extern "C" fn losat_web2_register_inputs(
    argv_ptr: *const u8,
    argv_len: usize,
    query_ptr: *const u8,
    query_len: usize,
    subject_ptr: *const u8,
    subject_len: usize,
) -> i32 {
    guarded(|| {
        let argv = words(text(argv_ptr, argv_len, "argv")?);
        let (_, _, response) = store::register_inputs(
            &argv,
            bytes(query_ptr, query_len),
            bytes(subject_ptr, subject_len),
        )?;
        emit(STREAM_RESPONSE, response.as_bytes());
        Ok(0)
    })
}

#[no_mangle]
pub extern "C" fn losat_web2_release(handle: u32) -> i32 {
    guarded(|| store::release(handle).map(|()| 0))
}

/// Starts a named input file of `kind` (`1`: a WindowMasker unit counts file, the value of
/// BLASTN's `-window_masker_db`) of `size` bytes; returns its handle.
///
/// # Safety
/// The pointer arguments must describe readable buffers of the given lengths.
#[no_mangle]
pub unsafe extern "C" fn losat_web2_file_begin(
    kind: u32,
    name_ptr: *const u8,
    name_len: usize,
    size: usize,
) -> i32 {
    guarded(|| {
        let handle = files::begin(kind, text(name_ptr, name_len, "file name")?, size)?;
        i32::try_from(handle).map_err(|_| "input file handle space exhausted".to_string())
    })
}

/// # Safety
/// The pointer arguments must describe readable buffers of the given lengths.
#[no_mangle]
pub unsafe extern "C" fn losat_web2_file_chunk(file: u32, ptr: *const u8, len: usize) -> i32 {
    guarded(|| files::chunk(file, bytes(ptr, len)).map(|()| 0))
}

/// Ends a named input file: emits a *file* JSON on stream 2.
#[no_mangle]
pub extern "C" fn losat_web2_file_end(file: u32) -> i32 {
    guarded(|| {
        let response = files::end(file)?;
        emit(STREAM_RESPONSE, response.as_bytes());
        Ok(0)
    })
}

#[no_mangle]
pub extern "C" fn losat_web2_file_release(file: u32) -> i32 {
    guarded(|| files::release(file).map(|()| 0))
}

fn scanners() -> &'static Mutex<(u32, HashMap<u32, scan::Scanner>)> {
    static SCANNERS: OnceLock<Mutex<(u32, HashMap<u32, scan::Scanner>)>> = OnceLock::new();
    SCANNERS.get_or_init(Mutex::default)
}

/// Starts an index scan of `program`'s input of `role` (`0` query, `1` subject) with NCBI's data
/// loaders on (`data_loaders` 1) or off (0) and `-parse_deflines` on (`parse_deflines` 1) or off
/// (0): the engine's reader of that input, as `register` / `register_inputs` read it (NCBI BLAST+'s
/// reader of nucleotide input for BLASTN, TBLASTX, the TBLASTN subject and the BLASTX query; of
/// protein input for BLASTP, the TBLASTN query and the BLASTX subject).
///
/// # Safety
/// The pointer arguments must describe readable buffers of the given lengths.
#[no_mangle]
pub unsafe extern "C" fn losat_web2_scan_begin(
    program_ptr: *const u8,
    program_len: usize,
    role: u32,
    data_loaders: u32,
    parse_deflines: u32,
) -> i32 {
    guarded(|| {
        let program = run::Program::parse(text(program_ptr, program_len, "program")?)?;
        store::role_name(role)?;
        let data_loaders = flag(data_loaders, "data_loaders")?;
        let parse_deflines = flag(parse_deflines, "parse_deflines")?;
        let scanner = scan::Scanner::new(
            store::reader_config(program, role, data_loaders).with_parse_deflines(parse_deflines),
        );
        let mut scanners = scanners().lock().expect("scanners");
        scanners.0 = scanners
            .0
            .checked_add(1)
            .ok_or_else(|| "scanner handle space exhausted".to_string())?;
        let id = scanners.0;
        scanners.1.insert(id, scanner);
        i32::try_from(id).map_err(|_| "scanner handle space exhausted".to_string())
    })
}

/// # Safety
/// The pointer arguments must describe readable buffers of the given lengths.
#[no_mangle]
pub unsafe extern "C" fn losat_web2_scan_chunk(scanner: u32, ptr: *const u8, len: usize) -> i32 {
    guarded(|| {
        let mut scanners = scanners().lock().expect("scanners");
        let state = scanners
            .1
            .get_mut(&scanner)
            .ok_or_else(|| format!("unknown scanner {scanner}"))?;
        state.feed(bytes(ptr, len));
        Ok(0)
    })
}

#[no_mangle]
pub extern "C" fn losat_web2_scan_end(scanner: u32) -> i32 {
    guarded(|| {
        let state = scanners()
            .lock()
            .expect("scanners")
            .1
            .remove(&scanner)
            .ok_or_else(|| format!("unknown scanner {scanner}"))?;
        let scanned = state
            .finish()
            .map_err(|error| error.message().to_string())?;
        emit(STREAM_RESPONSE, scan::to_json(&scanned).as_bytes());
        Ok(0)
    })
}

fn inflations() -> &'static Mutex<(u32, HashMap<u32, inflate::Inflation>)> {
    static INFLATIONS: OnceLock<Mutex<(u32, HashMap<u32, inflate::Inflation>)>> = OnceLock::new();
    INFLATIONS.get_or_init(Mutex::default)
}

/// Starts the inflation of a file named `name` (UTF-8), as the command line reads a file of that
/// name: `.gz`, `.bz2` or `.zst` in any case gives an inflation handle and an *inflate* JSON on
/// stream 2 (the compression and the name without the suffix); any other name returns `0`, no
/// handle: the file is read as it is.
///
/// # Safety
/// The pointer arguments must describe readable buffers of the given lengths.
#[no_mangle]
pub unsafe extern "C" fn losat_web2_inflate_begin(name_ptr: *const u8, name_len: usize) -> i32 {
    guarded(|| {
        let Some((inflation, response)) =
            inflate::Inflation::begin(text(name_ptr, name_len, "name")?)
        else {
            return Ok(0);
        };
        let handle = {
            let mut inflations = inflations().lock().expect("inflations");
            let id = inflations
                .0
                .checked_add(1)
                .filter(|id| i32::try_from(*id).is_ok())
                .ok_or_else(|| "inflation handle space exhausted".to_string())?;
            inflations.0 = id;
            inflations.1.insert(id, inflation);
            id as i32
        };
        emit(STREAM_RESPONSE, response.as_bytes());
        Ok(handle)
    })
}

/// Inflates the next bytes of the file; the inflated bytes go out on stream 5. A file that is not
/// what its name says fails with the command line's sentence, and its handle is dropped.
///
/// # Safety
/// The pointer arguments must describe readable buffers of the given lengths.
#[no_mangle]
pub unsafe extern "C" fn losat_web2_inflate_chunk(handle: u32, ptr: *const u8, len: usize) -> i32 {
    guarded(|| {
        // The inflation leaves the map while it runs, so that no lock is held while `emit` calls
        // the host (as `inflate_end`); it goes back only when the chunk inflated.
        let mut inflation = inflations()
            .lock()
            .expect("inflations")
            .1
            .remove(&handle)
            .ok_or_else(|| format!("unknown inflation {handle}"))?;
        inflation.chunk(bytes(ptr, len))?;
        inflations()
            .lock()
            .expect("inflations")
            .1
            .insert(handle, inflation);
        Ok(0)
    })
}

/// Ends the file: the rest of the inflated bytes go out on stream 5, or the call fails with the
/// command line's sentence of data cut short (or of an empty file). The handle is dropped either
/// way; a host that stops before the end of the file calls this and ignores the result.
#[no_mangle]
pub extern "C" fn losat_web2_inflate_end(handle: u32) -> i32 {
    guarded(|| {
        let inflation = inflations()
            .lock()
            .expect("inflations")
            .1
            .remove(&handle)
            .ok_or_else(|| format!("unknown inflation {handle}"))?;
        inflation.end()?;
        Ok(0)
    })
}

/// Runs the argv on the registered query and subject. A run of `-db` takes the subject handle 0:
/// the engine opens the database by name in the run's working directory, which the host preopens as
/// `/work` (docs/web/abi_v2.md §3, §4). The inputs are those read with the argv's data loaders
/// (`run::data_loaders`, §9).
///
/// # Safety
/// The pointer arguments must describe readable buffers of the given lengths.
#[no_mangle]
pub unsafe extern "C" fn losat_web2_run(
    argv_ptr: *const u8,
    argv_len: usize,
    query_handle: u32,
    subject_handle: u32,
) -> i32 {
    guarded(|| {
        let argv = words(text(argv_ptr, argv_len, "argv")?);
        let program = run::Program::parse(argv.first().copied().unwrap_or(""))?;
        // A run that prints a text in place of a search reads no input, as the command line reads
        // no file for it: its handles are not read (docs/web/abi_v2.md §4; line WQA's D14).
        if run::prints(&argv).is_some() {
            run::run(&argv, &[], &losat::api::QueryEnd::Input, &[])?;
            return Ok(0);
        }
        store::with_inputs(
            program,
            query_handle,
            subject_handle,
            &argv,
            run::data_loaders(&argv).ok(),
            |queries, end, subjects| run::run(&argv, queries, end, subjects),
        )??;
        Ok(0)
    })
}

/// Runs a tool of the command line with its words, `argv[0]` its name (`windowmasker`; the other
/// programs fail with `<name> is not available in the browser (ABI v2)`), in the run's working
/// directory, with no standard input: its standard output on stream 11 and its standard error on
/// stream 3 as they are written, each file that it writes on stream 12 as it is written (also when
/// the run fails), then its exit status on stream 2 as `{"status":N}` (docs/web/abi_v2.md §4, §5,
/// §9; `tool.rs`). Returns 0 for every exit status; -1 only for an argv that is not UTF-8 or a name
/// that is not a tool of LOSAT Web.
///
/// # Safety
/// `argv_ptr` must describe a readable buffer of `argv_len` bytes.
#[no_mangle]
pub unsafe extern "C" fn losat_web2_tool(argv_ptr: *const u8, argv_len: usize) -> i32 {
    guarded(|| {
        tool::run(&words(text(argv_ptr, argv_len, "argv")?))?;
        Ok(0)
    })
}
