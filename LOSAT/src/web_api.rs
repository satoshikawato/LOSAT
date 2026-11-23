use std::collections::HashMap;
use std::path::PathBuf;
use std::slice;
use std::sync::{Mutex, OnceLock};

// ABI v1's searches of BLASTN, BLASTP and TBLASTX (plan TD-1): the v1 layer makes its frozen
// checks of the request and reads its inputs with NCBI's reader (`v1_input`).
mod v1_blastn;
mod v1_blastp;
mod v1_input;
mod v1_tblastx;

use crate::algorithm::blastp::args::{BlastpCompBasedStats, BlastpSegSpec};
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
use crate::algorithm::{blastn, blastp, blastx, tblastn, tblastx};
use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};

static LAST_RESULT: OnceLock<Mutex<Vec<u8>>> = OnceLock::new();
static LAST_ERROR: OnceLock<Mutex<Vec<u8>>> = OnceLock::new();
static FASTA_STORE: OnceLock<Mutex<FastaStore>> = OnceLock::new();

fn last_result() -> &'static Mutex<Vec<u8>> {
    LAST_RESULT.get_or_init(|| Mutex::new(Vec::new()))
}

fn last_error() -> &'static Mutex<Vec<u8>> {
    LAST_ERROR.get_or_init(|| Mutex::new(Vec::new()))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:606-617
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1427
fn fasta_store() -> &'static Mutex<FastaStore> {
    FASTA_STORE.get_or_init(|| Mutex::new(FastaStore::default()))
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
struct StoredFasta {
    bytes: Vec<u8>,
    label: String,
}

#[derive(Default)]
struct FastaStore {
    next_handle: i32,
    entries: HashMap<i32, StoredFasta>,
}

impl FastaStore {
    fn insert(&mut self, records: usize, label: String, bytes: Vec<u8>) -> Result<i32, String> {
        if records == 0 {
            return Err("FASTA store input did not contain any records".to_string());
        }
        if self.next_handle <= 0 {
            self.next_handle = 1;
        }
        let handle = self.next_handle;
        self.next_handle = self
            .next_handle
            .checked_add(1)
            .ok_or_else(|| "FASTA store handle space exhausted".to_string())?;
        self.entries.insert(handle, StoredFasta { label, bytes });
        Ok(handle)
    }

    fn release(&mut self, handle: i32) -> Result<(), String> {
        if self.entries.remove(&handle).is_some() {
            Ok(())
        } else {
            Err(format!("unknown FASTA handle: {handle}"))
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1427
    fn clear(&mut self) {
        self.entries.clear();
        self.next_handle = 1;
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:260-266
// Catch standard exceptions and report their message through `LOG_POST`.
// Preserve the Rust engine context and its underlying cause across the ABI.
fn engine_error(err: anyhow::Error) -> String {
    format!("{err:#}")
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

unsafe fn read_bytes<'a>(ptr: *const u8, len: usize, name: &str) -> Result<&'a [u8], String> {
    if ptr.is_null() {
        if len == 0 {
            return Ok(&[]);
        }
        return Err(format!("{name} pointer is null"));
    }
    // Safety: JS passes pointers allocated in this module's linear memory along
    // with their byte lengths; null is handled above.
    Ok(slice::from_raw_parts(ptr, len))
}

/// The records of a stored input, which the store only counts: each run reads the stored bytes
/// with the reader of its program and role (`run_pair_handles`). The store does not know the
/// molecule, so it reads with protein flags, which take every nucleotide letter.
fn count_stored_records(bytes: &[u8]) -> Result<usize, String> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:485-651
    read_all(
        &mut FastaInputSource::from_bytes(bytes, ReaderConfig::subject("BLASTP", true, false)),
        &mut |_| Ok(()),
    )
    .map(|records| records.len())
    .map_err(|err| format!("failed to parse stored FASTA: {err}"))
}

fn parse_extra_args(raw: &str) -> Vec<&str> {
    raw.split('\0').filter(|value| !value.is_empty()).collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
//
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1100-1108
fn web_outfmt_or_default(outfmt: &str) -> &str {
    if outfmt.trim().is_empty() {
        "6"
    } else {
        outfmt
    }
}

fn next_arg<'a>(args: &'a [&str], index: &mut usize, flag: &str) -> Result<&'a str, String> {
    *index += 1;
    args.get(*index)
        .copied()
        .ok_or_else(|| format!("{flag} requires a value"))
}

fn parse_num_threads_arg(value: &str, flag: &str) -> Result<usize, String> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
    value
        .parse()
        .map_err(|err| format!("{flag} parse error: {err}"))
}

/// `-searchsp` and `-dbsize` of ABI v1's BLASTN, TBLASTX and BLASTP, as the CLI reads them (added
/// on gbdraw's request, 2026-11-10: a comparison's E-values take the scale of its genomes). The
/// value is the CLI's own parser's, in the two spellings of the other options of ABI v1 (the next
/// word, or joined with `=`), and goes to the same field of the args structure as on the CLI.
/// `Ok(true)` when the word at `index` was one of them (the value is consumed); `Ok(false)` for
/// any other word.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:234-238
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2570-2572
fn take_search_space_arg(
    args: &[&str],
    index: &mut usize,
    searchsp: &mut Option<i64>,
    dbsize: &mut Option<i64>,
) -> Result<bool, String> {
    use crate::blastinput::value_parsers::ncbi_int8;
    use crate::web_v1_argv::nonnegative_ncbi_int8;
    let flag = args[*index];
    let (name, joined) = match flag.split_once('=') {
        Some((name, value)) => (name, Some(value)),
        None => (flag, None),
    };
    let (target, parse): (&mut Option<i64>, fn(&str) -> Result<i64, String>) = match name {
        "-searchsp" => (searchsp, nonnegative_ncbi_int8),
        "-dbsize" => (dbsize, ncbi_int8),
        _ => return Ok(false),
    };
    let value = match joined {
        Some(value) => value,
        None => next_arg(args, index, flag)?,
    };
    *target = Some(parse(value).map_err(|err| format!("{name} parse error: {err}"))?);
    Ok(true)
}

fn parse_blastn_args(
    extra_args: &[&str],
    query: PathBuf,
    subject: PathBuf,
    out: PathBuf,
) -> Result<blastn::BlastnArgs, String> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
    let mut args = blastn::BlastnArgs {
        query,
        subject: Some(subject),
        // The frozen web ABI v1 (plan TD-1) has no database search.
        db: None,
        // The Web ABI searches a subject, not a database: no taxid list (line DI-8).
        db_filter: Default::default(),
        // ... and no database masks (series NP, M4a-2).
        db_mask: Default::default(),
        db_mask_resolved: None,
        // ... and no MegaBLAST database index (line RG item D5b).
        megablast_index_loaded: None,
        formatting: Default::default(),
        dbsize: None,
        show_gis: false,
        html: false,
        html_links: None,
        parse_deflines: false,
        id_scope: None,
        searchsp: None,
        task: "megablast".to_string(),
        template_type: None,
        template_length: None,
        word_size: None,
        num_threads: 1,
        mt_mode: None,
        use_index: None,
        index_name: None,
        evalue: Some(10.0),
        percent_identity: None,
        max_target_seqs: None,
        num_descriptions: None,
        num_alignments: None,
        report_counts: Default::default(),
        limit_lookup: false,
        max_db_word_count: 30,
        max_hsps_per_subject: None,
        qcov_hsp_perc: None,
        min_diag_separation: 0,
        out: Some(out),
        reward: None,
        penalty: None,
        gap_open: None,
        gap_extend: None,
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:46-48
        // Default DUST settings: `kDustLevel` 20, `kDustWindow` 64, `kDustLinker` 1.
        dust_filtering: None,
        dust: crate::blastinput::value_parsers::DustSpec::Yes,
        lcase_masking: false,
        // The frozen web ABI v1 has no `-soft_masking` (refused below as an unknown flag).
        soft_masking: None,
        // The frozen web ABI v1 (plan TD-1) has no WindowMasker.
        window_masker_taxid: None,
        window_masker_db: None,
        window_masker: None,
        window_masker_counts: None,
        window_masker_registry: None,
        filtering_db: None,
        repeat_filter: None,
        // ABI v1 writes no tabular field list.
        tabular_taxdb: None,
        subject_taxdb: None,
        taxdb_setting: None,
        strand: "both".to_string(),
        // The frozen web ABI v1 (plan TD-1) has no ranges.
        query_loc: None,
        subject_loc: None,
        subject_besthit: false,
        culling_limit: None,
        best_hit_overhang: None,
        best_hit_score_edge: None,
        // The frozen web ABI v1 has no `-ungapped`, X-drops, `-no_greedy`, `-window_size`,
        // `-off_diagonal_range` or `-min_raw_gapped_score` (they are refused below as unknown
        // flags).
        ungapped: false,
        xdrop_ungap: None,
        xdrop_gap: None,
        xdrop_gap_final: None,
        no_greedy: false,
        min_raw_gapped_score: None,
        window_size: None,
        off_diagonal_range: None,
        verbose: false,
        scan_step: 0,
        outfmt: "6".to_string(),
        archive: None,
        export_search_strategy: None,
        strategy: None,
        command_line: None,
        // The command line's progress bar; a web search writes nothing to standard error.
        progress_bar: false,
        one_query_batch: false,
        hits_only: false,
        first_query_ordinal: 0,
        scope_refusal: None,
    };

    let mut index = 0;
    while index < extra_args.len() {
        let flag = extra_args[index];
        if take_search_space_arg(extra_args, &mut index, &mut args.searchsp, &mut args.dbsize)? {
            index += 1;
            continue;
        }
        if let Some(value) = flag.strip_prefix("-task=") {
            args.task = value.to_string();
        } else if flag == "-task" {
            args.task = next_arg(extra_args, &mut index, flag)?.to_string();
        } else if flag == "-outfmt" {
            args.outfmt = next_arg(extra_args, &mut index, flag)?.to_string();
        } else if let Some(value) = flag.strip_prefix("-outfmt=") {
            args.outfmt = value.to_string();
        } else if flag == "-word_size" {
            args.word_size = v1_blastn_word_size(
                next_arg(extra_args, &mut index, flag)?
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("-word_size=") {
            args.word_size = v1_blastn_word_size(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if flag == "-num_threads" {
            args.num_threads =
                parse_num_threads_arg(next_arg(extra_args, &mut index, flag)?, flag)?;
        } else if let Some(value) = flag.strip_prefix("-num_threads=") {
            args.num_threads = parse_num_threads_arg(value, flag)?;
        } else if flag == "-evalue" {
            // The forms that NCBI's CArg_Double reads, as the CLI (`blastn_evalue`).
            args.evalue = Some(
                crate::blastinput::value_parsers::ncbi_double(next_arg(
                    extra_args, &mut index, flag,
                )?)
                .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("-evalue=") {
            args.evalue = Some(
                crate::blastinput::value_parsers::ncbi_double(value)
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else {
            return Err(format!("unsupported blastn argument for web API: {flag}"));
        }
        index += 1;
    }
    // Plan TD-1 (a fail-fast fix): the engine implements the tasks of the CLI only; v1
    // ran other names with LOSAT's own scoring, which NCBI does not use.
    if !matches!(args.task.as_str(), "megablast" | "blastn") {
        return Err(format!(
            "unsupported blastn task for web API: {}",
            args.task
        ));
    }

    Ok(args)
}

/// The v1 word size: plan TD-1 freezes v1, whose engine took a word size of 28 for the
/// megablast default and so gave `-task blastn -word_size 28` the blastn default (11).
fn v1_blastn_word_size(word_size: usize) -> Option<usize> {
    (word_size != 28).then_some(word_size)
}

fn parse_tblastx_args(
    extra_args: &[&str],
    query: PathBuf,
    subject: PathBuf,
    out: PathBuf,
) -> Result<tblastx::TblastxArgs, String> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastx_args.cpp:66-118
    let mut args = tblastx::TblastxArgs {
        query,
        subject: Some(subject),
        // A database search is not part of LOSAT Web yet (DB-7).
        db: None,
        // The Web ABI searches a subject, not a database: no taxid list (line DI-8).
        db_filter: Default::default(),
        // ... and no database masks (series NP, M4b).
        db_mask: Default::default(),
        formatting: Default::default(),
        dbsize: None,
        show_gis: false,
        html: false,
        parse_deflines: false,
        id_scope: None,
        searchsp: None,
        evalue: 10.0,
        // The frozen web ABI v1 (plan TD-1) has no -threshold, -matrix or -window_size: BLOSUM62's
        // suggestions, 13 and 40.
        threshold: None,
        matrix: None,
        word_size: None,
        num_threads: 1,
        out: Some(out),
        query_gencode: 1,
        db_gencode: 1,
        max_intron_length: None,
        max_target_seqs: None,
        num_descriptions: None,
        num_alignments: None,
        report_counts: Default::default(),
        html_links: None,
        max_hsps: None,
        qcov_hsp_perc: None,
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:83-85
        // Default SEG parameters from `kSegWindow`, `kSegLocut`, and `kSegHicut`.
        seg: "yes".to_string(),
        soft_masking: None,
        sum_stats: None,
        lcase_masking: false,
        // The frozen web ABI v1 (plan TD-1) searches both strands of the query.
        strand: "both".to_string(),
        xdrop_ungap: None,
        window_size: None,
        outfmt: "6".to_string(),
        culling_limit: None,
        best_hit_overhang: None,
        best_hit_score_edge: None,
        subject_besthit: false,
        // The frozen web ABI v1 (plan TD-1) has no ranges.
        query_loc: None,
        subject_loc: None,
        archive: None,
        export_search_strategy: None,
        strategy: None,
        tabular_taxdb: None,
        subject_taxdb: None,
        taxdb_setting: None,
        first_query_ordinal: 0,
        scope_refusal: None,
        matrix_file: None,
    };

    let mut index = 0;
    while index < extra_args.len() {
        let flag = extra_args[index];
        if take_search_space_arg(extra_args, &mut index, &mut args.searchsp, &mut args.dbsize)? {
            index += 1;
            continue;
        }
        if flag == "-query_gencode" {
            args.query_gencode = next_arg(extra_args, &mut index, flag)?
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if let Some(value) = flag.strip_prefix("-query_gencode=") {
            args.query_gencode = value
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if flag == "-db_gencode" {
            args.db_gencode = next_arg(extra_args, &mut index, flag)?
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if let Some(value) = flag.strip_prefix("-db_gencode=") {
            args.db_gencode = value
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if flag == "-num_threads" {
            args.num_threads =
                parse_num_threads_arg(next_arg(extra_args, &mut index, flag)?, flag)?;
        } else if let Some(value) = flag.strip_prefix("-num_threads=") {
            args.num_threads = parse_num_threads_arg(value, flag)?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2657-2660
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2800-2803
        // Plan TD-1 freezes ABI v1, whose TBLASTX wrote only outfmt 6: it keeps rejecting
        // any other format with the error that it gave before the engine implemented
        // outfmt 0 and 7 (session S08).
        } else if flag == "-outfmt" {
            args.outfmt = crate::blastinput::value_parsers::tblastx_v1_outfmt(next_arg(
                extra_args, &mut index, flag,
            )?)?;
        } else if let Some(value) = flag.strip_prefix("-outfmt=") {
            args.outfmt = crate::blastinput::value_parsers::tblastx_v1_outfmt(value)?;
        } else if flag == "-evalue" {
            args.evalue = next_arg(extra_args, &mut index, flag)?
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if let Some(value) = flag.strip_prefix("-evalue=") {
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:825-886
fn parse_blastp_comp_based_stats(value: &str) -> Result<BlastpCompBasedStats, String> {
    crate::algorithm::blastp::args::parse_comp_based_stats(value)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:386-407
fn parse_blastp_seg(value: &str) -> Result<BlastpSegSpec, String> {
    crate::blastinput::value_parsers::parse_seg_filtering(value)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastp_args.cpp:71-115
//
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:85-148
fn parse_blastp_args(
    extra_args: &[&str],
    query: PathBuf,
    subject: PathBuf,
    out: PathBuf,
) -> Result<blastp::BlastpArgs, String> {
    let mut args = blastp::BlastpArgs {
        query,
        subject: Some(subject),
        // The frozen web ABI v1 (plan TD-1) searches subject records, not a database.
        db: None,
        // The Web ABI searches a subject, not a database: no taxid list (line DI-8).
        db_filter: Default::default(),
        // ... and no database masks (series NP, M4b).
        db_mask: Default::default(),
        dbsize: None,
        show_gis: false,
        html: false,
        parse_deflines: false,
        id_scope: None,
        searchsp: None,
        task: "blastp".to_string(),
        evalue: None,
        threshold: None,
        word_size: None,
        num_threads: 1,
        mt_mode: None,
        out: Some(out),
        max_target_seqs: None,
        num_descriptions: None,
        num_alignments: None,
        max_hsps_per_subject: None,
        qcov_hsp_perc: None,
        culling_limit: None,
        best_hit_overhang: None,
        best_hit_score_edge: None,
        subject_besthit: false,
        ungapped: false,
        xdrop_ungap: None,
        xdrop_gap: None,
        xdrop_gap_final: None,
        window_size: None,
        matrix: None,
        gap_open: None,
        gap_extend: None,
        comp_based_stats: "2".to_string(),
        seg: None,
        soft_masking: None,
        lcase_masking: false,
        use_sw_tback: false,
        // The frozen web ABI v1 (plan TD-1) has no ranges.
        query_loc: None,
        subject_loc: None,
        outfmt: "6".to_string(),
        archive: None,
        export_search_strategy: None,
        strategy: None,
        request_inputs: None,
        app: blastp::args::BlastpApp::Blastp,
        taxdb_setting: None,
        matrix_search: None,
        inclusion_ethresh: None,
        gap_trigger: None,
        sum_stats: false,
        first_query_ordinal: 0,
        scope_refusal: None,
        psiblast_pssm: None,
        formatting: Default::default(),
        rps: None,
    };

    let mut index = 0;
    while index < extra_args.len() {
        let flag = extra_args[index];
        if take_search_space_arg(extra_args, &mut index, &mut args.searchsp, &mut args.dbsize)? {
            index += 1;
            continue;
        }
        if let Some(value) = flag.strip_prefix("-task=") {
            args.task = value.to_string();
        } else if flag == "-task" {
            args.task = next_arg(extra_args, &mut index, flag)?.to_string();
        } else if flag == "-outfmt" {
            args.outfmt = next_arg(extra_args, &mut index, flag)?.to_string();
        } else if let Some(value) = flag.strip_prefix("-outfmt=") {
            args.outfmt = value.to_string();
        } else if flag == "-evalue" {
            args.evalue = Some(
                next_arg(extra_args, &mut index, flag)?
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("-evalue=") {
            args.evalue = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if flag == "-threshold" {
            args.threshold = Some(
                next_arg(extra_args, &mut index, flag)?
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("-threshold=") {
            args.threshold = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if flag == "-word_size" {
            args.word_size = Some(
                next_arg(extra_args, &mut index, flag)?
                    .parse::<usize>()
                    .map_err(|err| format!("{flag} parse error: {err}"))?
                    .min(i32::MAX as usize) as i32,
            );
        } else if let Some(value) = flag.strip_prefix("-word_size=") {
            args.word_size = Some(
                value
                    .parse::<usize>()
                    .map_err(|err| format!("{flag} parse error: {err}"))?
                    .min(i32::MAX as usize) as i32,
            );
        } else if flag == "-num_threads" {
            args.num_threads = next_arg(extra_args, &mut index, flag)?
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if let Some(value) = flag.strip_prefix("-num_threads=") {
            args.num_threads = value
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if flag == "-max_target_seqs" {
            args.max_target_seqs = Some(
                next_arg(extra_args, &mut index, flag)?
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("-max_target_seqs=") {
            args.max_target_seqs = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if flag == "-max_hsps" {
            args.max_hsps_per_subject = Some(
                next_arg(extra_args, &mut index, flag)?
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("-max_hsps=") {
            args.max_hsps_per_subject = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if flag == "-window_size" {
            args.window_size = Some(
                next_arg(extra_args, &mut index, flag)?
                    .parse::<usize>()
                    .map_err(|err| format!("{flag} parse error: {err}"))?
                    .min(i32::MAX as usize) as i32,
            );
        } else if let Some(value) = flag.strip_prefix("-window_size=") {
            args.window_size = Some(
                value
                    .parse::<usize>()
                    .map_err(|err| format!("{flag} parse error: {err}"))?
                    .min(i32::MAX as usize) as i32,
            );
        } else if flag == "-matrix" {
            args.matrix = Some(next_arg(extra_args, &mut index, flag)?.to_string());
        } else if let Some(value) = flag.strip_prefix("-matrix=") {
            args.matrix = Some(value.to_string());
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:93-94:
        } else if flag == "-gapopen" {
            args.gap_open = Some(
                next_arg(extra_args, &mut index, flag)?
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("-gapopen=") {
            args.gap_open = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if flag == "-gapextend" {
            args.gap_extend = Some(
                next_arg(extra_args, &mut index, flag)?
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("-gapextend=") {
            args.gap_extend = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if flag == "-comp_based_stats" {
            let value = next_arg(extra_args, &mut index, flag)?;
            parse_blastp_comp_based_stats(value)?;
            args.comp_based_stats = value.to_string();
        } else if let Some(value) = flag.strip_prefix("-comp_based_stats=") {
            parse_blastp_comp_based_stats(value)?;
            args.comp_based_stats = value.to_string();
        } else if flag == "-seg" {
            let value = next_arg(extra_args, &mut index, flag)?;
            parse_blastp_seg(value)?;
            args.seg = Some(value.to_string());
        } else if let Some(value) = flag.strip_prefix("-seg=") {
            parse_blastp_seg(value)?;
            args.seg = Some(value.to_string());
        } else if flag == "-ungapped" {
            args.ungapped = true;
        } else if let Some(value) = flag.strip_prefix("-ungapped=") {
            args.ungapped = value.eq_ignore_ascii_case("true") || value == "1";
        } else if flag == "-use_sw_tback" {
            args.use_sw_tback = true;
        } else if let Some(value) = flag.strip_prefix("-use_sw_tback=") {
            args.use_sw_tback = value.eq_ignore_ascii_case("true") || value == "1";
        } else {
            return Err(format!("unsupported blastp argument for web API: {flag}"));
        }
        index += 1;
    }

    Ok(args)
}

fn run_pair(
    program: &str,
    query_fasta: &str,
    subject_fasta: &str,
    outfmt: &str,
    extra_args: &str,
) -> Result<Vec<u8>, String> {
    let mut extra = parse_extra_args(extra_args);
    let outfmt = web_outfmt_or_default(outfmt);
    extra.push("-outfmt");
    extra.push(outfmt);

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:485-651
    match program {
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
        "blastx" => blastx::web::run_pair(
            blastx::web::parse_args(&extra).map_err(engine_error)?,
            query_fasta.as_bytes(),
            subject_fasta.as_bytes(),
            "subject",
        )
        .map_err(engine_error),
        "blastp" => {
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastp_args.cpp:71-115
            v1_blastp::run_web_pair(
                parse_blastp_args(
                    &extra,
                    std::path::PathBuf::new(),
                    std::path::PathBuf::new(),
                    std::path::PathBuf::new(),
                )?,
                query_fasta,
                subject_fasta,
            )
            .map_err(engine_error)
        }
        "blastn" => v1_blastn::run_web_pair(
            parse_blastn_args(
                &extra,
                std::path::PathBuf::new(),
                std::path::PathBuf::new(),
                std::path::PathBuf::new(),
            )?,
            query_fasta,
            subject_fasta,
        )
        .map_err(engine_error),
        "tblastx" => v1_tblastx::run_web_pair(
            parse_tblastx_args(
                &extra,
                std::path::PathBuf::new(),
                std::path::PathBuf::new(),
                std::path::PathBuf::new(),
            )?,
            query_fasta,
            subject_fasta,
        )
        .map_err(engine_error),
        // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:213-219
        // Added to ABI v1 on the Owner's request (2026-11-04): the CLI's parser and run over
        // the call's bytes (`tblastn::web`).
        "tblastn" => tblastn::web::run_pair(
            tblastn::web::parse_args(&extra).map_err(engine_error)?,
            query_fasta.as_bytes(),
            subject_fasta.as_bytes(),
            "query",
            "subject",
        )
        .map_err(engine_error),
        other => Err(format!("unsupported LOSAT program: {other}")),
    }
}

fn run_pair_handles(
    program: &str,
    query_handle: i32,
    subject_handle: i32,
    outfmt: &str,
    extra_args: &str,
) -> Result<Vec<u8>, String> {
    let mut extra = parse_extra_args(extra_args);
    let outfmt = web_outfmt_or_default(outfmt);
    extra.push("-outfmt");
    extra.push(outfmt);

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:606-617
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1427
    match program {
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:277-279
        "blastx" => {
            let store = fasta_store()
                .lock()
                .map_err(|_| "FASTA store mutex poisoned".to_string())?;
            let query = store
                .entries
                .get(&query_handle)
                .ok_or_else(|| format!("unknown query FASTA handle: {query_handle}"))?;
            let subject = store
                .entries
                .get(&subject_handle)
                .ok_or_else(|| format!("unknown subject FASTA handle: {subject_handle}"))?;
            blastx::web::run_pair(
                blastx::web::parse_args(&extra).map_err(engine_error)?,
                &query.bytes,
                &subject.bytes,
                &subject.label,
            )
            .map_err(engine_error)
        }
        "blastp" => {
            let store = fasta_store()
                .lock()
                .map_err(|_| "FASTA store mutex poisoned".to_string())?;
            let query = store
                .entries
                .get(&query_handle)
                .ok_or_else(|| format!("unknown query FASTA handle: {query_handle}"))?;
            let subject = store
                .entries
                .get(&subject_handle)
                .ok_or_else(|| format!("unknown subject FASTA handle: {subject_handle}"))?;
            v1_blastp::run_web_pair_bytes(
                parse_blastp_args(
                    &extra,
                    std::path::PathBuf::new(),
                    std::path::PathBuf::new(),
                    std::path::PathBuf::new(),
                )?,
                &query.bytes,
                &subject.bytes,
                &query.label,
                &subject.label,
            )
            .map_err(engine_error)
        }
        // The stored bytes, read by the CLI's reader of each role (`tblastn::web`); the labels
        // are the names the report shows for the inputs.
        "tblastn" => {
            let store = fasta_store()
                .lock()
                .map_err(|_| "FASTA store mutex poisoned".to_string())?;
            let query = store
                .entries
                .get(&query_handle)
                .ok_or_else(|| format!("unknown query FASTA handle: {query_handle}"))?;
            let subject = store
                .entries
                .get(&subject_handle)
                .ok_or_else(|| format!("unknown subject FASTA handle: {subject_handle}"))?;
            tblastn::web::run_pair(
                tblastn::web::parse_args(&extra).map_err(engine_error)?,
                &query.bytes,
                &subject.bytes,
                &query.label,
                &subject.label,
            )
            .map_err(engine_error)
        }
        "blastn" | "tblastx" => Err(format!(
            "FASTA handle API is not yet ported for {program}; use losat_web_run_pair"
        )),
        other => Err(format!("unsupported LOSAT program: {other}")),
    }
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
pub unsafe extern "C" fn losat_web_store_fasta(
    fasta_ptr: *const u8,
    fasta_len: usize,
    label_ptr: *const u8,
    label_len: usize,
) -> i32 {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:485-651
    let fasta_bytes = match read_bytes(fasta_ptr, fasta_len, "FASTA") {
        Ok(value) => value,
        Err(err) => return set_error(err),
    };
    let label = match read_str(label_ptr, label_len, "FASTA label") {
        Ok(value) => value.to_string(),
        Err(err) => return set_error(err),
    };
    let records = match count_stored_records(fasta_bytes) {
        Ok(value) => value,
        Err(err) => return set_error(err),
    };
    match fasta_store()
        .lock()
        .map_err(|_| "FASTA store mutex poisoned".to_string())
        .and_then(|mut store| store.insert(records, label, fasta_bytes.to_vec()))
    {
        Ok(handle) => {
            last_error().lock().expect("error mutex poisoned").clear();
            handle
        }
        Err(err) => set_error(err),
    }
}

#[no_mangle]
pub extern "C" fn losat_web_release_fasta(handle: i32) -> i32 {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1427
    match fasta_store()
        .lock()
        .map_err(|_| "FASTA store mutex poisoned".to_string())
        .and_then(|mut store| store.release(handle))
    {
        Ok(()) => {
            last_error().lock().expect("error mutex poisoned").clear();
            0
        }
        Err(err) => set_error(err),
    }
}

#[no_mangle]
pub extern "C" fn losat_web_clear_fasta_store() -> i32 {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1427
    match fasta_store().lock() {
        Ok(mut store) => {
            store.clear();
            last_error().lock().expect("error mutex poisoned").clear();
            0
        }
        Err(_) => set_error("FASTA store mutex poisoned"),
    }
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
pub unsafe extern "C" fn losat_web_run_pair_handles(
    program_ptr: *const u8,
    program_len: usize,
    query_handle: i32,
    subject_handle: i32,
    outfmt_ptr: *const u8,
    outfmt_len: usize,
    extra_args_ptr: *const u8,
    extra_args_len: usize,
) -> i32 {
    let program = match read_str(program_ptr, program_len, "program") {
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

    match run_pair_handles(program, query_handle, subject_handle, outfmt, extra_args) {
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

#[cfg(test)]
#[path = "web_api_tests.rs"]
mod tests;
