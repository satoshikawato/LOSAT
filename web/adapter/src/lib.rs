//! LOSAT Web ABI v2 (docs/web/abi_v2.md): the exports that the application's engine
//! worker calls. ABI v1 (`losat_web_*`, `LOSAT/src/web_api.rs`) is linked in unchanged
//! and shares no state with these exports (plan TD-1).
//!
//! Every export returns `0` or a handle on success and `-1` on failure, with the
//! message in `losat_web2_last_error_ptr/len`. Responses and outputs go to the host
//! through `losat_host.emit` (`emit.rs`).

pub mod describe;
pub mod emit;
pub mod json;
pub mod run;
pub mod scan;
pub mod store;

use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

use emit::{emit, STREAM_RESPONSE};

fn last_error() -> &'static Mutex<Vec<u8>> {
    static LAST_ERROR: OnceLock<Mutex<Vec<u8>>> = OnceLock::new();
    LAST_ERROR.get_or_init(Mutex::default)
}

/// Runs an export body: records its error, or clears the last error on success.
fn guarded(body: impl FnOnce() -> Result<i32, String>) -> i32 {
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

#[no_mangle]
pub extern "C" fn losat_web2_abi_version() -> i32 {
    2
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
pub unsafe extern "C" fn losat_web2_validate(argv_ptr: *const u8, argv_len: usize) -> i32 {
    guarded(|| {
        run::validate(&words(text(argv_ptr, argv_len, "argv")?))?;
        Ok(0)
    })
}

/// # Safety
/// The pointer arguments must describe readable buffers of the given lengths.
#[no_mangle]
pub unsafe extern "C" fn losat_web2_register(
    program_ptr: *const u8,
    program_len: usize,
    role: u32,
    bytes_ptr: *const u8,
    bytes_len: usize,
) -> i32 {
    guarded(|| {
        let program = text(program_ptr, program_len, "program")?;
        let (handle, response) = store::register(program, role, bytes(bytes_ptr, bytes_len))?;
        emit(STREAM_RESPONSE, response.as_bytes());
        i32::try_from(handle).map_err(|_| "input handle space exhausted".to_string())
    })
}

#[no_mangle]
pub extern "C" fn losat_web2_release(handle: u32) -> i32 {
    guarded(|| store::release(handle).map(|()| 0))
}

fn scanners() -> &'static Mutex<(u32, HashMap<u32, scan::Scanner>)> {
    static SCANNERS: OnceLock<Mutex<(u32, HashMap<u32, scan::Scanner>)>> = OnceLock::new();
    SCANNERS.get_or_init(Mutex::default)
}

/// Starts an index scan with parser kind `0` (`bio::io::fasta`, the parser of BLASTP,
/// TBLASTN, BLASTN and TBLASTX). Kind `1` (BLASTX's NCBI-style reader) joins in SX.
#[no_mangle]
pub extern "C" fn losat_web2_scan_begin(parser: u32) -> i32 {
    guarded(|| {
        if parser != 0 {
            return Err(format!("unknown or unavailable FASTA parser kind {parser}"));
        }
        let mut scanners = scanners().lock().expect("scanners");
        scanners.0 = scanners
            .0
            .checked_add(1)
            .ok_or_else(|| "scanner handle space exhausted".to_string())?;
        let id = scanners.0;
        scanners.1.insert(id, scan::Scanner::new());
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
        let records = state.finish()?;
        emit(STREAM_RESPONSE, scan::to_json(&records).as_bytes());
        Ok(0)
    })
}

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
        store::with_inputs(
            program,
            query_handle,
            subject_handle,
            |queries, subjects| run::run(&argv, queries, subjects),
        )??;
        Ok(0)
    })
}
