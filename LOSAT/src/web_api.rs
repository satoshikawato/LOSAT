use std::path::PathBuf;
use std::slice;
use std::sync::{Mutex, OnceLock};

use crate::algorithm::{blastn, tblastx};

static LAST_RESULT: OnceLock<Mutex<Vec<u8>>> = OnceLock::new();
static LAST_ERROR: OnceLock<Mutex<Vec<u8>>> = OnceLock::new();

fn last_result() -> &'static Mutex<Vec<u8>> {
    LAST_RESULT.get_or_init(|| Mutex::new(Vec::new()))
}

fn last_error() -> &'static Mutex<Vec<u8>> {
    LAST_ERROR.get_or_init(|| Mutex::new(Vec::new()))
}

fn set_result(bytes: Vec<u8>) {
    *last_result().lock().expect("result mutex poisoned") = bytes;
    last_error().lock().expect("error mutex poisoned").clear();
}

fn set_error(message: impl AsRef<str>) -> i32 {
    last_result().lock().expect("result mutex poisoned").clear();
    *last_error().lock().expect("error mutex poisoned") = message.as_ref().as_bytes().to_vec();
    -1
}

unsafe fn read_str<'a>(ptr: *const u8, len: usize, name: &str) -> Result<&'a str, String> {
    if ptr.is_null() {
        if len == 0 {
            return Ok("");
        }
        return Err(format!("{name} pointer is null"));
    }
    // Safety: JS passes pointers allocated in this module's linear memory along
    // with their byte lengths; null is handled above and UTF-8 is validated below.
    let bytes = slice::from_raw_parts(ptr, len);
    std::str::from_utf8(bytes).map_err(|err| format!("{name} is not valid UTF-8: {err}"))
}

fn parse_extra_args(raw: &str) -> Vec<&str> {
    raw.split('\0').filter(|value| !value.is_empty()).collect()
}

fn next_arg<'a>(args: &'a [&str], index: &mut usize, flag: &str) -> Result<&'a str, String> {
    *index += 1;
    args.get(*index)
        .copied()
        .ok_or_else(|| format!("{flag} requires a value"))
}

fn parse_blastn_args(
    extra_args: &[&str],
    query: PathBuf,
    subject: PathBuf,
    out: PathBuf,
) -> Result<blastn::BlastnArgs, String> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
    // ```c
    // const string kArgQuery("query");
    // const string kArgOutput("out");
    // const string kArgSubject("subject");
    // const string kTask("task");
    // ```
    let mut args = blastn::BlastnArgs {
        query,
        subject,
        task: "megablast".to_string(),
        word_size: 28,
        num_threads: 1,
        evalue: 10.0,
        percent_identity: 0.0,
        min_hit_length: 0,
        max_target_seqs: None,
        hitlist_size: 500,
        limit_lookup: false,
        max_db_word_count: 30,
        max_hsps_per_subject: 0,
        min_diag_separation: 0,
        out: Some(out),
        reward: 1,
        penalty: -2,
        gap_open: 0,
        gap_extend: 0,
        dust: true,
        dust_level: 20,
        dust_window: 64,
        dust_linker: 1,
        lcase_masking: false,
        subject_besthit: false,
        verbose: false,
        scan_step: 0,
        outfmt: "6".to_string(),
    };

    let mut index = 0;
    while index < extra_args.len() {
        let flag = extra_args[index];
        if let Some(value) = flag.strip_prefix("--task=") {
            args.task = value.to_string();
        } else if flag == "--task" {
            args.task = next_arg(extra_args, &mut index, flag)?.to_string();
        } else if flag == "--outfmt" {
            args.outfmt = next_arg(extra_args, &mut index, flag)?.to_string();
        } else if let Some(value) = flag.strip_prefix("--outfmt=") {
            args.outfmt = value.to_string();
        } else if flag == "--word-size" || flag == "--word_size" {
            args.word_size = next_arg(extra_args, &mut index, flag)?
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if let Some(value) = flag.strip_prefix("--word-size=") {
            args.word_size = value
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if flag == "--evalue" {
            args.evalue = next_arg(extra_args, &mut index, flag)?
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if let Some(value) = flag.strip_prefix("--evalue=") {
            args.evalue = value
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else {
            return Err(format!("unsupported blastn argument for web API: {flag}"));
        }
        index += 1;
    }

    Ok(args)
}

fn parse_tblastx_args(
    extra_args: &[&str],
    query: PathBuf,
    subject: PathBuf,
    out: PathBuf,
) -> Result<tblastx::TblastxArgs, String> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/tblastx_args.cpp:66-118
    // ```c
    // arg.Reset(new CGenericSearchArgs( !kQueryIsProtein, false, false, true));
    // arg.Reset(new CGeneticCodeArgs(CGeneticCodeArgs::eQuery));
    // arg.Reset(new CGeneticCodeArgs(CGeneticCodeArgs::eDatabase));
    // arg.Reset(new CFormattingArgs);
    // ```
    let mut args = tblastx::TblastxArgs {
        query,
        subject,
        evalue: 10.0,
        threshold: 13,
        word_size: 3,
        num_threads: 1,
        out: Some(out),
        query_gencode: 1,
        db_gencode: 1,
        max_target_seqs: 500,
        seg: true,
        seg_window: 12,
        seg_locut: 2.2,
        seg_hicut: 2.5,
        window_size: 40,
        outfmt: "6".to_string(),
        culling_limit: 0,
    };

    let mut index = 0;
    while index < extra_args.len() {
        let flag = extra_args[index];
        if flag == "--query-gencode" {
            args.query_gencode = next_arg(extra_args, &mut index, flag)?
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if let Some(value) = flag.strip_prefix("--query-gencode=") {
            args.query_gencode = value
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if flag == "--db-gencode" {
            args.db_gencode = next_arg(extra_args, &mut index, flag)?
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if let Some(value) = flag.strip_prefix("--db-gencode=") {
            args.db_gencode = value
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if flag == "--outfmt" {
            args.outfmt = next_arg(extra_args, &mut index, flag)?.to_string();
        } else if let Some(value) = flag.strip_prefix("--outfmt=") {
            args.outfmt = value.to_string();
        } else if flag == "--evalue" {
            args.evalue = next_arg(extra_args, &mut index, flag)?
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if let Some(value) = flag.strip_prefix("--evalue=") {
            args.evalue = value
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else {
            return Err(format!("unsupported tblastx argument for web API: {flag}"));
        }
        index += 1;
    }

    Ok(args)
}

fn remove_file_if_present(path: &str) {
    let _ = std::fs::remove_file(path);
}

fn run_pair(
    program: &str,
    query_fasta: &str,
    subject_fasta: &str,
    outfmt: &str,
    extra_args: &str,
) -> Result<Vec<u8>, String> {
    let query_path = PathBuf::from(".losat_web_query.fa");
    let subject_path = PathBuf::from(".losat_web_subject.fa");
    let out_path = PathBuf::from(".losat_web_out.txt");

    remove_file_if_present(query_path.to_string_lossy().as_ref());
    remove_file_if_present(subject_path.to_string_lossy().as_ref());
    remove_file_if_present(out_path.to_string_lossy().as_ref());

    std::fs::write(&query_path, query_fasta)
        .map_err(|err| format!("failed to write query FASTA: {err}"))?;
    std::fs::write(&subject_path, subject_fasta)
        .map_err(|err| format!("failed to write subject FASTA: {err}"))?;

    let mut extra = parse_extra_args(extra_args);
    extra.push("--outfmt");
    extra.push(outfmt);

    match program {
        "blastn" => blastn::run(parse_blastn_args(
            &extra,
            query_path.clone(),
            subject_path.clone(),
            out_path.clone(),
        )?)
        .map_err(|err| err.to_string())?,
        "tblastx" => tblastx::run(parse_tblastx_args(
            &extra,
            query_path.clone(),
            subject_path.clone(),
            out_path.clone(),
        )?)
        .map_err(|err| err.to_string())?,
        other => return Err(format!("unsupported LOSAT program: {other}")),
    }

    let output =
        std::fs::read(&out_path).map_err(|err| format!("failed to read LOSAT output: {err}"))?;
    remove_file_if_present(query_path.to_string_lossy().as_ref());
    remove_file_if_present(subject_path.to_string_lossy().as_ref());
    remove_file_if_present(out_path.to_string_lossy().as_ref());
    Ok(output)
}

#[no_mangle]
pub extern "C" fn losat_web_alloc(len: usize) -> *mut u8 {
    let mut buffer = Vec::<u8>::with_capacity(len);
    let ptr = buffer.as_mut_ptr();
    std::mem::forget(buffer);
    ptr
}

#[no_mangle]
pub unsafe extern "C" fn losat_web_dealloc(ptr: *mut u8, len: usize) {
    if ptr.is_null() {
        return;
    }
    // Safety: `ptr` must come from `losat_web_alloc(len)`, which creates a Vec
    // with capacity `len` and transfers ownership to the caller.
    drop(Vec::from_raw_parts(ptr, 0, len));
}

#[no_mangle]
pub unsafe extern "C" fn losat_web_run_pair(
    program_ptr: *const u8,
    program_len: usize,
    query_ptr: *const u8,
    query_len: usize,
    subject_ptr: *const u8,
    subject_len: usize,
    outfmt_ptr: *const u8,
    outfmt_len: usize,
    extra_args_ptr: *const u8,
    extra_args_len: usize,
) -> i32 {
    let program = match read_str(program_ptr, program_len, "program") {
        Ok(value) => value,
        Err(err) => return set_error(err),
    };
    let query = match read_str(query_ptr, query_len, "query FASTA") {
        Ok(value) => value,
        Err(err) => return set_error(err),
    };
    let subject = match read_str(subject_ptr, subject_len, "subject FASTA") {
        Ok(value) => value,
        Err(err) => return set_error(err),
    };
    let outfmt = match read_str(outfmt_ptr, outfmt_len, "outfmt") {
        Ok(value) => value,
        Err(err) => return set_error(err),
    };
    let extra_args = match read_str(extra_args_ptr, extra_args_len, "extra args") {
        Ok(value) => value,
        Err(err) => return set_error(err),
    };

    match run_pair(program, query, subject, outfmt, extra_args) {
        Ok(output) => {
            set_result(output);
            0
        }
        Err(err) => set_error(err),
    }
}

#[no_mangle]
pub extern "C" fn losat_web_result_ptr() -> *const u8 {
    last_result()
        .lock()
        .expect("result mutex poisoned")
        .as_ptr()
}

#[no_mangle]
pub extern "C" fn losat_web_result_len() -> usize {
    last_result().lock().expect("result mutex poisoned").len()
}

#[no_mangle]
pub extern "C" fn losat_web_error_ptr() -> *const u8 {
    last_error().lock().expect("error mutex poisoned").as_ptr()
}

#[no_mangle]
pub extern "C" fn losat_web_error_len() -> usize {
    last_error().lock().expect("error mutex poisoned").len()
}

#[no_mangle]
pub extern "C" fn losat_web_clear() {
    last_result().lock().expect("result mutex poisoned").clear();
    last_error().lock().expect("error mutex poisoned").clear();
}

#[no_mangle]
pub extern "C" fn losat_web_has_direct_api() -> i32 {
    1
}
