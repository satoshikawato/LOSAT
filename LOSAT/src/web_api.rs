use std::collections::HashMap;
use std::path::PathBuf;
use std::slice;
use std::sync::{Mutex, OnceLock};

use bio::io::fasta;

use crate::algorithm::blastp::args::{BlastpCompBasedStats, BlastpCompositionMode, BlastpSegSpec};
use crate::algorithm::{blastn, blastp, tblastx};

static LAST_RESULT: OnceLock<Mutex<Vec<u8>>> = OnceLock::new();
static LAST_ERROR: OnceLock<Mutex<Vec<u8>>> = OnceLock::new();
static FASTA_STORE: OnceLock<Mutex<FastaStore>> = OnceLock::new();

fn last_result() -> &'static Mutex<Vec<u8>> {
    LAST_RESULT.get_or_init(|| Mutex::new(Vec::new()))
}

fn last_error() -> &'static Mutex<Vec<u8>> {
    LAST_ERROR.get_or_init(|| Mutex::new(Vec::new()))
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:606-617
// ```c
// sequence = queries.GetBlastSequence(index, encoding,
//                                     eNa_strand_unknown,
//                                     eSentinels,
//                                     &warnings);
// memcpy(&buf.get()[offset], sequence.data.get(), sequence.length);
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
// ```c
// db_length = BlastSeqSrcGetTotLen(seq_src);
// while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
//     if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
//         continue;
//     }
// }
// ```
fn fasta_store() -> &'static Mutex<FastaStore> {
    FASTA_STORE.get_or_init(|| Mutex::new(FastaStore::default()))
}

struct StoredFasta {
    records: Vec<fasta::Record>,
    label: String,
}

#[derive(Default)]
struct FastaStore {
    next_handle: i32,
    entries: HashMap<i32, StoredFasta>,
}

impl FastaStore {
    fn insert(&mut self, records: Vec<fasta::Record>, label: String) -> Result<i32, String> {
        if records.is_empty() {
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
        self.entries.insert(handle, StoredFasta { records, label });
        Ok(handle)
    }

    fn release(&mut self, handle: i32) -> Result<(), String> {
        if self.entries.remove(&handle).is_some() {
            Ok(())
        } else {
            Err(format!("unknown FASTA handle: {handle}"))
        }
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
    // ```c
    // db_length = BlastSeqSrcGetTotLen(seq_src);
    // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    // ```
    fn clear(&mut self) {
        self.entries.clear();
        self.next_handle = 1;
    }
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

fn parse_fasta_records(bytes: &[u8], name: &str) -> Result<Vec<fasta::Record>, String> {
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:486-651
    // ```c
    // SetupQueries_OMF(IBlastQuerySource& queries,
    //                  BlastQueryInfo* qinfo,
    //                  BLAST_SequenceBlk** seqblk,
    //                  EBlastProgramType prog,
    //                  ...)
    // ```
    fasta::Reader::new(bytes)
        .records()
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|err| format!("failed to parse {name}: {err}"))
}

fn parse_extra_args(raw: &str) -> Vec<&str> {
    raw.split('\0').filter(|value| !value.is_empty()).collect()
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
// ```c
// const string kArgOutputFormat("outfmt");
// const string kDfltArgOutputFormat("0");
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
// ```c
// ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
//     if (iter != m_FieldsToShow.begin())
//         m_Ostream << m_FieldDelimiter;
//     x_PrintField(*iter);
// }
// ```
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
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
    // ```c
    // const string kArgNumThreads("num_threads");
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
    // ```c
    // int num_threads = args[kArgNumThreads].AsInteger();
    // if (num_threads > kMaxValue) {
    //     m_NumThreads = kMaxValue;
    // } else {
    //     m_NumThreads = num_threads;
    // }
    // ```
    value
        .parse()
        .map_err(|err| format!("{flag} parse error: {err}"))
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
        } else if flag == "--num-threads"
            || flag == "--num_threads"
            || flag == "-num_threads"
            || flag == "-num-threads"
        {
            args.num_threads =
                parse_num_threads_arg(next_arg(extra_args, &mut index, flag)?, flag)?;
        } else if let Some(value) = flag.strip_prefix("--num-threads=") {
            args.num_threads = parse_num_threads_arg(value, flag)?;
        } else if let Some(value) = flag.strip_prefix("--num_threads=") {
            args.num_threads = parse_num_threads_arg(value, flag)?;
        } else if let Some(value) = flag.strip_prefix("-num_threads=") {
            args.num_threads = parse_num_threads_arg(value, flag)?;
        } else if let Some(value) = flag.strip_prefix("-num-threads=") {
            args.num_threads = parse_num_threads_arg(value, flag)?;
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
        } else if flag == "--num-threads"
            || flag == "--num_threads"
            || flag == "-num_threads"
            || flag == "-num-threads"
        {
            args.num_threads =
                parse_num_threads_arg(next_arg(extra_args, &mut index, flag)?, flag)?;
        } else if let Some(value) = flag.strip_prefix("--num-threads=") {
            args.num_threads = parse_num_threads_arg(value, flag)?;
        } else if let Some(value) = flag.strip_prefix("--num_threads=") {
            args.num_threads = parse_num_threads_arg(value, flag)?;
        } else if let Some(value) = flag.strip_prefix("-num_threads=") {
            args.num_threads = parse_num_threads_arg(value, flag)?;
        } else if let Some(value) = flag.strip_prefix("-num-threads=") {
            args.num_threads = parse_num_threads_arg(value, flag)?;
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

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:825-886
// ```c
// switch (comp_stat_string[0]) {
// case '0': case 'F': case 'f': compo_mode = eNoCompositionBasedStats; break;
// case '1': compo_mode = eCompositionBasedStats; break;
// case '2': compo_mode = eCompositionMatrixAdjust; break;
// case '3': compo_mode = eCompoForceFullMatrixAdjust; break;
// }
// if (program == eBlastp && compo_mode != eNoCompositionBasedStats &&
//     tolower(comp_stat_string[1]) == 'u') {
//     opt.SetUnifiedP(1);
// }
// ```
fn parse_blastp_comp_based_stats(value: &str) -> Result<BlastpCompBasedStats, String> {
    let mut chars = value.chars();
    let Some(first) = chars.next() else {
        return Err("composition-based statistics option cannot be empty".to_string());
    };

    let mode = match first {
        '0' | 'F' | 'f' => BlastpCompositionMode::NoCompositionBasedStats,
        '1' => BlastpCompositionMode::CompositionBasedStats,
        '2' | 'D' | 'd' | 'T' | 't' => BlastpCompositionMode::CompositionMatrixAdjust,
        '3' => BlastpCompositionMode::ForceFullMatrixAdjust,
        _ => {
            return Err(format!(
                "invalid composition-based statistics mode '{value}'"
            ))
        }
    };

    let unified_p = mode != BlastpCompositionMode::NoCompositionBasedStats
        && chars
            .next()
            .map(|c| c.eq_ignore_ascii_case(&'u'))
            .unwrap_or(false);

    Ok(BlastpCompBasedStats { mode, unified_p })
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:386-407
// ```c
// if (m_QueryIsProtein && args[kArgSegFiltering]) {
//     const string& seg_opts = args[kArgSegFiltering].AsString();
//     if (seg_opts == kDfltArgNoFiltering) {
//         opt.SetSegFiltering(false);
//     } else if (seg_opts == kDfltArgApplyFiltering) {
//         opt.SetSegFiltering(true);
//     } else {
//         x_TokenizeFilteringArgs(seg_opts, tokens);
//     }
// }
// ```
fn parse_blastp_seg(value: &str) -> Result<BlastpSegSpec, String> {
    if value.eq_ignore_ascii_case("no") || value.eq_ignore_ascii_case("false") || value == "0" {
        return Ok(BlastpSegSpec::No);
    }
    if value.eq_ignore_ascii_case("yes") || value.eq_ignore_ascii_case("true") || value == "1" {
        return Ok(BlastpSegSpec::Yes);
    }

    let tokens: Vec<&str> = value.split_whitespace().collect();
    if tokens.len() != 3 {
        return Err("invalid number of arguments to filtering option".to_string());
    }

    let window = tokens[0]
        .parse::<usize>()
        .map_err(|_| "invalid input for filtering parameters".to_string())?;
    let locut = tokens[1]
        .parse::<f64>()
        .map_err(|_| "invalid input for filtering parameters".to_string())?;
    let hicut = tokens[2]
        .parse::<f64>()
        .map_err(|_| "invalid input for filtering parameters".to_string())?;

    Ok(BlastpSegSpec::WindowLocutHicut {
        window,
        locut,
        hicut,
    })
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blastp_args.cpp:71-115
// ```c
// arg.Reset(new CGenericSearchArgs(kQueryIsProtein, false, false, false, false, true));
// arg.Reset(new CFilteringArgs(kQueryIsProtein, kFilterByDefault));
// arg.Reset(new CMatrixNameArg);
// arg.Reset(new CWordThresholdArg);
// arg.Reset(new CWindowSizeArg);
// arg.Reset(new CFormattingArgs);
// arg.Reset(new CGappedArgs);
// arg.Reset(new CCompositionBasedStatsArgs);
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_prot_options.cpp:85-148
// ```c
// SetWordSize(BLAST_WORDSIZE_PROT);
// SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP);
// SetMatrixName(BLAST_DEFAULT_MATRIX);
// SetGapOpeningCost(BLAST_GAP_OPEN_PROT);
// SetGapExtensionCost(BLAST_GAP_EXTN_PROT);
// SetHitlistSize(500);
// SetEvalueThreshold(BLAST_EXPECT_VALUE);
// ```
fn parse_blastp_args(
    extra_args: &[&str],
    query: PathBuf,
    subject: PathBuf,
    out: PathBuf,
) -> Result<blastp::BlastpArgs, String> {
    let mut args = blastp::BlastpArgs {
        query,
        subject,
        task: "blastp".to_string(),
        evalue: None,
        threshold: None,
        word_size: None,
        num_threads: 1,
        out: Some(out),
        max_target_seqs: 500,
        max_hsps_per_subject: 0,
        ungapped: false,
        window_size: None,
        matrix: None,
        gap_open: None,
        gap_extend: None,
        comp_based_stats: None,
        seg: None,
        use_sw_tback: false,
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
        } else if flag == "--evalue" {
            args.evalue = Some(
                next_arg(extra_args, &mut index, flag)?
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("--evalue=") {
            args.evalue = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if flag == "--threshold" {
            args.threshold = Some(
                next_arg(extra_args, &mut index, flag)?
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("--threshold=") {
            args.threshold = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if flag == "--word-size" || flag == "--word_size" {
            args.word_size = Some(
                next_arg(extra_args, &mut index, flag)?
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("--word-size=") {
            args.word_size = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("--word_size=") {
            args.word_size = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if flag == "--num-threads" || flag == "--num_threads" {
            args.num_threads = next_arg(extra_args, &mut index, flag)?
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if let Some(value) = flag.strip_prefix("--num-threads=") {
            args.num_threads = value
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if let Some(value) = flag.strip_prefix("--num_threads=") {
            args.num_threads = value
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if flag == "--max-target-seqs" || flag == "--max_target_seqs" {
            args.max_target_seqs = next_arg(extra_args, &mut index, flag)?
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if let Some(value) = flag.strip_prefix("--max-target-seqs=") {
            args.max_target_seqs = value
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if let Some(value) = flag.strip_prefix("--max_target_seqs=") {
            args.max_target_seqs = value
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if flag == "--max-hsps-per-subject" || flag == "--max_hsps_per_subject" {
            args.max_hsps_per_subject = next_arg(extra_args, &mut index, flag)?
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if let Some(value) = flag.strip_prefix("--max-hsps-per-subject=") {
            args.max_hsps_per_subject = value
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if let Some(value) = flag.strip_prefix("--max_hsps_per_subject=") {
            args.max_hsps_per_subject = value
                .parse()
                .map_err(|err| format!("{flag} parse error: {err}"))?;
        } else if flag == "--window-size" || flag == "--window_size" {
            args.window_size = Some(
                next_arg(extra_args, &mut index, flag)?
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("--window-size=") {
            args.window_size = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("--window_size=") {
            args.window_size = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if flag == "--matrix" {
            args.matrix = Some(next_arg(extra_args, &mut index, flag)?.to_string());
        } else if let Some(value) = flag.strip_prefix("--matrix=") {
            args.matrix = Some(value.to_string());
        } else if flag == "--gap-open" || flag == "--gapopen" || flag == "--gap_open" {
            args.gap_open = Some(
                next_arg(extra_args, &mut index, flag)?
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("--gap-open=") {
            args.gap_open = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("--gapopen=") {
            args.gap_open = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("--gap_open=") {
            args.gap_open = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if flag == "--gap-extend" || flag == "--gapextend" || flag == "--gap_extend" {
            args.gap_extend = Some(
                next_arg(extra_args, &mut index, flag)?
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("--gap-extend=") {
            args.gap_extend = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("--gapextend=") {
            args.gap_extend = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if let Some(value) = flag.strip_prefix("--gap_extend=") {
            args.gap_extend = Some(
                value
                    .parse()
                    .map_err(|err| format!("{flag} parse error: {err}"))?,
            );
        } else if flag == "--comp-based-stats" || flag == "--comp_based_stats" {
            args.comp_based_stats = Some(parse_blastp_comp_based_stats(next_arg(
                extra_args, &mut index, flag,
            )?)?);
        } else if let Some(value) = flag.strip_prefix("--comp-based-stats=") {
            args.comp_based_stats = Some(parse_blastp_comp_based_stats(value)?);
        } else if let Some(value) = flag.strip_prefix("--comp_based_stats=") {
            args.comp_based_stats = Some(parse_blastp_comp_based_stats(value)?);
        } else if flag == "--seg" {
            args.seg = Some(parse_blastp_seg(next_arg(extra_args, &mut index, flag)?)?);
        } else if let Some(value) = flag.strip_prefix("--seg=") {
            args.seg = Some(parse_blastp_seg(value)?);
        } else if flag == "--ungapped" {
            args.ungapped = true;
        } else if let Some(value) = flag.strip_prefix("--ungapped=") {
            args.ungapped = value.eq_ignore_ascii_case("true") || value == "1";
        } else if flag == "--use_sw_tback" || flag == "--use-sw-tback" {
            args.use_sw_tback = true;
        } else if let Some(value) = flag.strip_prefix("--use_sw_tback=") {
            args.use_sw_tback = value.eq_ignore_ascii_case("true") || value == "1";
        } else if let Some(value) = flag.strip_prefix("--use-sw-tback=") {
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
    extra.push("--outfmt");
    extra.push(outfmt);

    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:486-651
    // ```c
    // SetupQueries_OMF(IBlastQuerySource& queries,
    //                  BlastQueryInfo* qinfo,
    //                  BLAST_SequenceBlk** seqblk,
    //                  EBlastProgramType prog,
    //                  ...)
    // ```
    match program {
        "blastp" => {
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blastp_args.cpp:71-115
            // ```c
            // arg.Reset(new CGenericSearchArgs(kQueryIsProtein, false, false,
            //                                  false, false, true));
            // arg.Reset(new CFormattingArgs);
            // arg.Reset(new CGappedArgs);
            // arg.Reset(new CCompositionBasedStatsArgs);
            // ```
            blastp::run_web_pair(
                parse_blastp_args(
                    &extra,
                    std::path::PathBuf::new(),
                    std::path::PathBuf::new(),
                    std::path::PathBuf::new(),
                )?,
                query_fasta,
                subject_fasta,
            )
            .map_err(|err: anyhow::Error| err.to_string())
        }
        "blastn" => blastn::run_web_pair(
            parse_blastn_args(
                &extra,
                std::path::PathBuf::new(),
                std::path::PathBuf::new(),
                std::path::PathBuf::new(),
            )?,
            query_fasta,
            subject_fasta,
        )
        .map_err(|err: anyhow::Error| err.to_string()),
        "tblastx" => tblastx::run_web_pair(
            parse_tblastx_args(
                &extra,
                std::path::PathBuf::new(),
                std::path::PathBuf::new(),
                std::path::PathBuf::new(),
            )?,
            query_fasta,
            subject_fasta,
        )
        .map_err(|err: anyhow::Error| err.to_string()),
        other => return Err(format!("unsupported LOSAT program: {other}")),
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
    extra.push("--outfmt");
    extra.push(outfmt);

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:606-617
    // ```c
    // sequence = queries.GetBlastSequence(index, encoding,
    //                                     eNa_strand_unknown,
    //                                     eSentinels,
    //                                     &warnings);
    // memcpy(&buf.get()[offset], sequence.data.get(), sequence.length);
    // ```
    //
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
    // ```c
    // db_length = BlastSeqSrcGetTotLen(seq_src);
    // while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
    //     if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
    //         continue;
    //     }
    // }
    // ```
    match program {
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
            blastp::run_web_pair_records(
                parse_blastp_args(
                    &extra,
                    std::path::PathBuf::new(),
                    std::path::PathBuf::new(),
                    std::path::PathBuf::new(),
                )?,
                &query.records,
                &subject.records,
                &query.label,
                &subject.label,
            )
            .map_err(|err: anyhow::Error| err.to_string())
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
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:486-651
    // ```c
    // SetupQueries_OMF(IBlastQuerySource& queries,
    //                  BlastQueryInfo* qinfo,
    //                  BLAST_SequenceBlk** seqblk,
    //                  EBlastProgramType prog,
    //                  ...)
    // ```
    let fasta_bytes = match read_bytes(fasta_ptr, fasta_len, "FASTA") {
        Ok(value) => value,
        Err(err) => return set_error(err),
    };
    let label = match read_str(label_ptr, label_len, "FASTA label") {
        Ok(value) => value.to_string(),
        Err(err) => return set_error(err),
    };
    let records = match parse_fasta_records(fasta_bytes, "stored FASTA") {
        Ok(value) => value,
        Err(err) => return set_error(err),
    };
    match fasta_store()
        .lock()
        .map_err(|_| "FASTA store mutex poisoned".to_string())
        .and_then(|mut store| store.insert(records, label))
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
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
    // ```c
    // while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
    //     if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
    //         continue;
    //     }
    // }
    // ```
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
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
    // ```c
    // db_length = BlastSeqSrcGetTotLen(seq_src);
    // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    // ```
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
mod tests {
    use super::*;

    fn record(id: &str) -> fasta::Record {
        fasta::Record::with_attrs(id, None, b"ARND")
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:486-651
    // ```c
    // SetupQueries_OMF(IBlastQuerySource& queries,
    //                  BlastQueryInfo* qinfo,
    //                  BLAST_SequenceBlk** seqblk,
    //                  EBlastProgramType prog,
    //                  ...)
    // ```
    #[test]
    fn fasta_store_allocates_and_releases_handles() {
        let mut store = FastaStore::default();

        let first = store
            .insert(vec![record("q1")], "query".to_string())
            .expect("first handle");
        let second = store
            .insert(vec![record("s1")], "subject".to_string())
            .expect("second handle");

        assert_eq!(first, 1);
        assert_eq!(second, 2);
        assert!(store.entries.contains_key(&first));
        store.release(first).expect("release existing handle");
        assert!(!store.entries.contains_key(&first));
        assert!(store.release(first).is_err());
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
    // ```c
    // while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
    //     if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
    //         continue;
    //     }
    // }
    // ```
    #[test]
    fn fasta_store_rejects_empty_and_unknown_handles() {
        let mut store = FastaStore::default();

        assert!(store.insert(Vec::new(), "empty".to_string()).is_err());
        assert!(store.release(42).is_err());
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
    // ```c
    // db_length = BlastSeqSrcGetTotLen(seq_src);
    // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
    // ```
    #[test]
    fn fasta_store_clear_removes_entries_and_resets_handles() {
        let mut store = FastaStore::default();
        let handle = store
            .insert(vec![record("q1")], "query".to_string())
            .expect("stored handle");

        assert_eq!(handle, 1);
        store.clear();
        assert!(store.entries.is_empty());

        let next = store
            .insert(vec![record("q2")], "query2".to_string())
            .expect("stored handle after clear");
        assert_eq!(next, 1);
    }
}
