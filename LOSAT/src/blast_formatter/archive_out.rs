//! `-outfmt 11` of `blast_formatter` (archive re-emission): NCBI's `CBlastFormat::WriteArchive`
//! and `BlastBuildArchive` (c++/src/algo/blast/format/blast_format.cpp:1066-1133,
//! build_archive.cpp:84-203) over the archive object that `CRemoteBlast::LoadFromArchive` read.
//!
//! Each search object writes one `Blast4-archive` in the encoding named by `ARCHIVE_FORMAT`
//! (`ArchiveEncoding::from_env`: ASN.1 text by default, `xml` for XML, `bin*` for ASN.1 binary,
//! anything else writes nothing). Unlike formats 8, 9 and 12, `-outfmt 11` does not prune the
//! alignments (`PruneSeqalign`) and does not run `BlastFormatter_PreFetchSequenceData`; instead,
//! `-max_target_seqs` and `-num_alignments` update `HitlistSize` in `algorithm-options`, and
//! stored `messages` are never written back into the emitted archive.

use std::cell::RefCell;
use std::collections::BTreeSet;
use std::io::Write;

use super::archive::{invalid_archive, stop_error, Output, SkipWriter, Spill, SPILL_LAG};
use super::refusals::Refusal;
use crate::blastinput::app::engine_error;
use crate::ncbi_app::Diag;
use crate::report::archive::fields::{self, Field};
use crate::report::archive::options::{Param, ParamList, ParamValue};
use crate::report::archive::read::import::recorded_params;
use crate::report::archive::read::options::{BuiltOptions, EProgram, Setter};
use crate::report::archive::read::request::{Bioseq, Queries};
use crate::report::archive::read::results::MaskRange;
use crate::report::archive::read::Search;
use crate::report::archive::writer::{write_ka_block, write_param, KaBlock};
use crate::report::archive::ArchiveEncoding;
use crate::report::seqid::SeqId;
use crate::report::serial::asn_ber::AsnBinaryEncoder;
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::asn_types::view::V;
use crate::report::serial::asn_types::write::write_value;
use crate::report::serial::encoder::{Encoder, ReplacedChar, SerialError, Traced};
use crate::report::serial::xml::XmlEncoder;
use crate::seqid::ParsedSeqId;

/// Default repeat filtering database (`kDefaultRepeatFilterDb` in `blast_options.h`).
const DEFAULT_REPEAT_FILTER_DB: &str = "repeat/repeat_9606";

/// Default DUST parameters (`SDustOptionsNew`: 20, 64, 1).
const DEFAULT_DUST: (i32, i32, i32) = (20, 64, 1);

/// Default SEG parameters (`SSegOptionsNew`: 12, 2.2, 2.5).
const DEFAULT_SEG: (i32, f64, f64) = (12, 2.2, 2.5);

/// `BLASTOPTIONS_BUFFER_SIZE` in `blast_filter.c`.
const FILTER_BUF_SIZE: usize = 128;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:257-290,543-650
/// The remote-options state of NCBI's cached `CBlastOptionsHandle` (`m_CBOH`), created once on
/// the first search object and mutated across objects by `ExtractAlgorithmOptions` (hit list
/// size) and `GetResultSet` (RPS matrix name).
pub(super) struct AlgorithmOptions {
    program: EProgram,
    program_name: &'static str,
    service_name: &'static str,
    /// `m_DefaultsMode`: left `true` by `SetVecScreenDefaults`, which skips every `SetValue`
    /// except `EvalueThreshold` (blast_nucl_options.cpp:93-112, blast_options_cxx.cpp:545-555).
    defaults_mode: bool,
    params: ParamList,
}

impl AlgorithmOptions {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:124-200
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:70-134
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/disc_nucl_options.cpp:47-74
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_options.cpp:63-75
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/deltablast_options.cpp:50-67
    /// Builds the remote-options state of `CBlastOptionsFactory::Create` for `options.program`
    /// with locality `eBoth` and replays every setter that `CBlastOptionsBuilder` applied to it.
    pub(super) fn from_built(options: &BuiltOptions) -> Self {
        let (program_name, service_name, defaults_mode) =
            remote_program_and_service(options.program);
        let mut state = Self {
            program: options.program,
            program_name,
            service_name,
            defaults_mode,
            params: ParamList::default(),
        };
        state.init_factory_defaults();
        for setter in &options.setters {
            state.apply_setter(setter);
        }
        state
    }

    fn init_factory_defaults(&mut self) {
        use ParamValue::{Boolean, Evalue, Integer};
        let evalue = if self.program == EProgram::VecScreen {
            700.0
        } else {
            10.0
        };
        self.params.set(fields::EVALUE_THRESHOLD, Evalue(evalue));
        match self.program {
            EProgram::DiscMegablast => {
                self.params.set(fields::MB_TEMPLATE_TYPE, Integer(0));
                self.params.set(fields::MB_TEMPLATE_LENGTH, Integer(18));
                self.params.set(fields::WORD_SIZE, Integer(11));
                self.params.set(fields::WINDOW_SIZE, Integer(40));
            }
            EProgram::PsiBlast => {
                self.set_psi_defaults(0.002);
            }
            EProgram::PsiTblastn => {
                self.set_psi_defaults(0.002);
                self.params.set(fields::SMITH_WATERMAN_MODE, Boolean(false));
                self.params.set(fields::DB_GENETIC_CODE, Integer(1));
                self.params.set(fields::WORD_THRESHOLD, Integer(13));
            }
            EProgram::DeltaBlast => {
                self.set_psi_defaults(0.05);
            }
            _ => {}
        }
    }

    fn set_psi_defaults(&mut self, inclusion_threshold: f64) {
        self.params.set(
            fields::INCLUSION_THRESHOLD,
            ParamValue::Real(inclusion_threshold),
        );
        self.params
            .set(fields::PSEUDO_COUNT, ParamValue::Integer(0));
    }

    fn set_param(&mut self, field: Field, value: ParamValue) {
        if field == fields::EVALUE_THRESHOLD || !self.defaults_mode {
            self.params.set(field, value);
        }
    }

    fn apply_setter(&mut self, setter: &Setter) {
        match setter {
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:1042-1130
            Setter::FilterString(text) => self.apply_filter_string(text),
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:1352-1364
            Setter::WindowMaskerTaxId(0) => self.params.reset(fields::WINDOW_MASKER_TAX_ID),
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:2619-2634
            Setter::UseIndex { index_name } => {
                self.set_param(fields::FORCE_MB_INDEX, ParamValue::Boolean(true));
                if !index_name.is_empty() {
                    self.set_param(
                        fields::MB_INDEX_NAME,
                        ParamValue::Str(index_name.as_bytes().to_vec()),
                    );
                }
            }
            other => {
                for param in recorded_params(other) {
                    self.set_param(param.field(), param.value);
                }
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:1022-1130
    fn apply_filter_string(&mut self, text: &str) {
        // `ClearFilterOptions()`
        self.set_param(fields::DUST_FILTERING, ParamValue::Boolean(false));
        self.set_param(fields::SEG_FILTERING, ParamValue::Boolean(false));
        self.set_param(fields::REPEAT_FILTERING, ParamValue::Boolean(false));
        self.set_param(fields::MASK_AT_HASH, ParamValue::Boolean(false));
        self.params.reset(fields::WINDOW_MASKER_TAX_ID);
        self.params.reset(fields::WINDOW_MASKER_DATABASE);

        let pure_nucl = is_pure_nucleotide_query(self.program);
        let parsed =
            parse_filter_string(text, is_blastn_or_mapper(self.program)).unwrap_or_default();

        self.set_param(
            fields::MASK_AT_HASH,
            ParamValue::Boolean(parsed.mask_at_hash),
        );
        let do_seg = if !pure_nucl {
            let on = parsed.seg.is_some();
            self.set_param(fields::SEG_FILTERING, ParamValue::Boolean(on));
            on
        } else {
            self.params.reset(fields::SEG_FILTERING);
            self.params.reset(fields::SEG_FILTERING_WINDOW);
            self.params.reset(fields::SEG_FILTERING_LOCUT);
            self.params.reset(fields::SEG_FILTERING_HICUT);
            false
        };
        let (do_dust, do_rep) = if pure_nucl {
            let dust = parsed.dust.is_some();
            let rep = is_blastn_or_mapper(self.program) && parsed.repeat_db.is_some();
            self.set_param(fields::DUST_FILTERING, ParamValue::Boolean(dust));
            self.set_param(fields::REPEAT_FILTERING, ParamValue::Boolean(rep));
            (dust, rep)
        } else {
            self.params.reset(fields::DUST_FILTERING);
            self.params.reset(fields::DUST_FILTERING_LEVEL);
            self.params.reset(fields::DUST_FILTERING_WINDOW);
            self.params.reset(fields::DUST_FILTERING_LINKER);
            self.params.reset(fields::REPEAT_FILTERING);
            self.params.reset(fields::REPEAT_FILTERING_DB);
            (false, false)
        };
        if do_dust {
            if let Some((level, window, linker)) = parsed.dust {
                self.set_param(fields::DUST_FILTERING_LEVEL, ParamValue::Integer(level));
                self.set_param(fields::DUST_FILTERING_WINDOW, ParamValue::Integer(window));
                self.set_param(fields::DUST_FILTERING_LINKER, ParamValue::Integer(linker));
            }
        }
        if do_rep {
            if let Some(db) = parsed.repeat_db {
                self.set_param(fields::REPEAT_FILTERING_DB, ParamValue::Str(db));
            }
        }
        if do_seg {
            if let Some((window, locut, hicut)) = parsed.seg {
                self.set_param(fields::SEG_FILTERING_WINDOW, ParamValue::Integer(window));
                self.set_param(fields::SEG_FILTERING_LOCUT, ParamValue::Real(locut));
                self.set_param(fields::SEG_FILTERING_HICUT, ParamValue::Real(hicut));
            }
        }
    }

    /// Updates `HitlistSize` as `CFormattingArgs::ExtractAlgorithmOptions` calls `SetHitlistSize`.
    pub(super) fn set_hitlist_size(&mut self, size: i32) {
        self.set_param(fields::HITLIST_SIZE, ParamValue::Integer(size));
    }

    /// Updates `MatrixName` as `CRemoteBlast::GetResultSet` calls `SetMatrixName` for
    /// `rpsblast` and `rpstblastn` (remote_blast.cpp:2205-2215).
    fn set_matrix_name(&mut self, matrix: &[u8]) {
        self.set_param(fields::MATRIX_NAME, ParamValue::Str(matrix.to_vec()));
    }

    fn params(&self) -> &[Param] {
        self.params.params()
    }
}

/// `GetRemoteProgramAndService_Blast3` for the factory handle of `program`, and whether
/// `SetDefaults` leaves `m_DefaultsMode` set (`eVecScreen`).
fn remote_program_and_service(program: EProgram) -> (&'static str, &'static str, bool) {
    match program {
        EProgram::Blastn => ("blastn", "plain", false),
        EProgram::Megablast | EProgram::DiscMegablast | EProgram::PhiBlastn => {
            ("blastn", "megablast", false)
        }
        EProgram::VecScreen => ("blastn", "vecscreen", true),
        EProgram::Blastp | EProgram::PhiBlastp => ("blastp", "plain", false),
        EProgram::Blastx => ("blastx", "plain", false),
        EProgram::Tblastn => ("tblastn", "plain", false),
        EProgram::Tblastx => ("tblastx", "plain", false),
        EProgram::RpsBlast => ("blastp", "rpsblast", false),
        EProgram::RpsTblastn => ("tblastn", "rpsblast", false),
        EProgram::PsiBlast | EProgram::PsiTblastn => ("blastp", "psi", false),
        EProgram::DeltaBlast => ("blastp", "delta_blast", false),
        EProgram::Mapper => ("blastn", "plain", false),
    }
}

fn is_translated_query(program: EProgram) -> bool {
    matches!(
        program,
        EProgram::Blastx | EProgram::Tblastx | EProgram::RpsTblastn
    )
}

fn is_pure_nucleotide_query(program: EProgram) -> bool {
    program.query_is_nucleotide() && !is_translated_query(program)
}

fn is_blastn_or_mapper(program: EProgram) -> bool {
    matches!(
        program,
        EProgram::Blastn
            | EProgram::Megablast
            | EProgram::DiscMegablast
            | EProgram::VecScreen
            | EProgram::Mapper
    )
}

#[derive(Default)]
struct ParsedFilter {
    mask_at_hash: bool,
    dust: Option<(i32, i32, i32)>,
    seg: Option<(i32, f64, f64)>,
    repeat_db: Option<Vec<u8>>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:54-97
/// `s_LoadOptionsToBuffer`: copies bytes up to `;` or `FILTER_BUF_SIZE` chars, trimming leading
/// and trailing ASCII spaces, and returns `(buffer, rest)`.
fn load_filter_buffer(input: &[u8]) -> (Vec<u8>, &[u8]) {
    let mut not_started = true;
    let mut buf = Vec::new();
    let mut i = 0;
    let mut count = 0;
    while count < FILTER_BUF_SIZE && i < input.len() {
        let b = input[i];
        i += 1;
        count += 1;
        if b == b';' {
            break;
        }
        if not_started && b == b' ' {
            continue;
        }
        not_started = false;
        buf.push(b);
    }
    if !not_started {
        while buf.len() > 1 && buf.last() == Some(&b' ') {
            buf.pop();
        }
    }
    (buf, &input[i..])
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:435-605
/// `BlastFilteringOptionsFromString`: returns `None` on parse error (leaving cleared options).
fn parse_filter_string(instructions: &str, blastn_or_mapper: bool) -> Option<ParsedFilter> {
    if instructions.eq_ignore_ascii_case("F") {
        return Some(ParsedFilter::default());
    }
    let mut out = ParsedFilter::default();
    let mut rest = instructions.as_bytes();
    if rest.starts_with(b"m ") {
        out.mask_at_hash = true;
        rest = &rest[2..];
    }
    while let Some(&ch) = rest.first() {
        match ch {
            b'S' | b'D' | b'R' | b'W' => {
                let (buf, next) = load_filter_buffer(&rest[1..]);
                rest = next;
                match ch {
                    b'S' => {
                        let seg = if buf.is_empty() {
                            DEFAULT_SEG
                        } else {
                            parse_seg_buffer(&buf)?
                        };
                        out.seg = Some(seg);
                    }
                    b'D' => {
                        let dust = if buf.is_empty() {
                            DEFAULT_DUST
                        } else {
                            parse_dust_buffer(&buf)?
                        };
                        out.dust = Some(dust);
                    }
                    b'R' => {
                        let mut db = DEFAULT_REPEAT_FILTER_DB.as_bytes().to_vec();
                        if !buf.is_empty() {
                            if let Some(parsed_db) = parse_repeat_buffer(&buf) {
                                db = parsed_db;
                            }
                        }
                        out.repeat_db = Some(db);
                    }
                    _ => {}
                }
            }
            b'L' | b'T' => {
                if blastn_or_mapper {
                    out.dust = Some(DEFAULT_DUST);
                } else {
                    out.seg = Some(DEFAULT_SEG);
                }
                rest = &rest[1..];
            }
            b'm' => {
                out.mask_at_hash = true;
                rest = &rest[1..];
            }
            _ => rest = &rest[1..],
        }
    }
    Some(out)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:180-234
fn parse_dust_buffer(buf: &[u8]) -> Option<(i32, i32, i32)> {
    let tokens = split_filter_tokens(buf);
    if tokens.is_empty() {
        return Some((-1, -1, -1));
    }
    if tokens.len() != 3 {
        return None;
    }
    Some((
        scan_c_long(tokens[0]) as i32,
        scan_c_long(tokens[1]) as i32,
        scan_c_long(tokens[2]) as i32,
    ))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:243-294
fn parse_seg_buffer(buf: &[u8]) -> Option<(i32, f64, f64)> {
    let tokens = split_filter_tokens(buf);
    if tokens.is_empty() {
        return Some((0, 0.0, 0.0));
    }
    if tokens.len() != 3 {
        return None;
    }
    Some((
        scan_c_long(tokens[0]) as i32,
        scan_c_double(tokens[1]),
        scan_c_double(tokens[2]),
    ))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:104-122
fn parse_repeat_buffer(buf: &[u8]) -> Option<Vec<u8>> {
    let pos = buf.windows(2).position(|w| w == b"-d")?;
    let mut slice = &buf[pos + 2..];
    while matches!(slice.first(), Some(b' ' | b'\t')) {
        slice = &slice[1..];
    }
    Some(slice.to_vec())
}

fn split_filter_tokens(buf: &[u8]) -> Vec<&[u8]> {
    let mut tokens = Vec::new();
    let mut i = 0;
    let limit = buf.len().min(FILTER_BUF_SIZE);
    while i < limit {
        let start = i;
        while i < limit && buf[i] != b' ' {
            i += 1;
        }
        tokens.push(&buf[start..i]);
        while i < limit && buf[i] == b' ' {
            i += 1;
        }
    }
    tokens
}

fn trimmed_utf8(bytes: &[u8]) -> &str {
    std::str::from_utf8(bytes).unwrap_or("").trim_start()
}

/// Leading decimal integer prefix (matching `s_ParseDustOptions` / `s_ParseSegOptions`), returning 0 on no match.
fn scan_c_long(bytes: &[u8]) -> i64 {
    let s = trimmed_utf8(bytes);
    let b = s.as_bytes();
    let mut end = usize::from(matches!(b.first(), Some(b'+' | b'-')));
    let digits_start = end;
    while end < b.len() && b[end].is_ascii_digit() {
        end += 1;
    }
    if end == digits_start {
        return 0;
    }
    s[..end].parse::<i64>().unwrap_or(0)
}

/// Leading floating-point prefix (matching `s_ParseSegOptions`), returning 0.0 on no match.
fn scan_c_double(bytes: &[u8]) -> f64 {
    let s = trimmed_utf8(bytes);
    let mut end = s.len();
    while end > 0 {
        if let Ok(v) = s[..end].parse::<f64>() {
            if !v.is_nan() {
                return v;
            }
        }
        end -= 1;
    }
    0.0
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1094-1108
/// `Blast_SubjectIsNucleotide` / `IsLocalId`: local or general Seq-ids keep `bioseq-set` instead
/// of `seq-loc-list` when `CExportStrategy::x_Process_Query` builds `Blast4-queries`.
fn is_local_id(id: &ParsedSeqId) -> bool {
    use crate::seqid::accguide::{identify_seq_id, parse_flags, GENERAL, LOCAL};
    let info = identify_seq_id(
        id,
        parse_flags::ANY_RAW | parse_flags::ANY_LOCAL,
        &mut |_| {},
    );
    matches!(id, ParsedSeqId::Local(_)) || info == LOCAL || info == GENERAL
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:569-625
/// `CPsiBlastValidate::Pssm` (with `require_scores` false): validates a `results.pssm` before
/// `BlastBuildArchive` emits it as `queries pssm`.
fn validate_pssm(pssm_with_params: V<'_>) -> anyhow::Result<()> {
    let Some(pssm) = pssm_with_params.get("pssm") else {
        return Err(engine_error("Missing PSSM data"));
    };
    let missing_scores = pssm
        .get("finalData")
        .and_then(|final_data| final_data.get("scores"))
        .is_none_or(|scores| scores.len() == 0);
    let missing_freq_ratios = pssm
        .get("intermediateData")
        .and_then(|inter| inter.get("freqRatios"))
        .is_none_or(|ratios| ratios.len() == 0);
    if missing_freq_ratios && missing_scores {
        return Err(engine_error(
            "PSSM data must contain either scores or frequency ratios",
        ));
    }
    if !missing_scores {
        if let Some(scale) = pssm
            .get("finalData")
            .and_then(|final_data| final_data.get("scalingFactor"))
            .map(|value| value.int())
        {
            if scale != 1 {
                return Err(engine_error(format!(
                    "PSSM has a scaling factor of {scale}. PSI-BLAST does not accept scaled PSSMs"
                )));
            }
        }
    }
    let Some(query) = pssm.get("query") else {
        return Err(engine_error("Missing query sequence in PSSM"));
    };
    if query.choice().0 != "seq" {
        return Err(engine_error(
            "Query sequence in ASN.1 PSSM is not a single Bioseq",
        ));
    }
    let is_protein = pssm
        .get("isProtein")
        .is_none_or(|v| matches!(v.value, crate::report::serial::asn_types::Value::Bool(true)));
    if !is_protein {
        return Err(engine_error(
            "PSSM does not represent protein scoring matrix",
        ));
    }
    Ok(())
}

/// The query part of the emitted `Blast4-queue-search-request`.
enum EmitQueries<'a> {
    Pssm(V<'a>),
    Bioseqs(Vec<Bioseq<'a>>),
    SeqLocs(Vec<&'a ParsedSeqId>),
}

/// The subject part of the emitted `Blast4-queue-search-request`.
enum EmitSubject<'a> {
    Database(&'a [u8]),
    Sequences(Vec<Bioseq<'a>>),
}

/// One `Blast4-mask` of the emitted results (`CRemoteBlast::ConvertToRemoteMasks`).
struct EmitMask<'a> {
    ranges: Vec<MaskRange<'a>>,
    frame_name: &'static str,
    frame_value: i32,
}

/// Everything needed to serialize one `Blast4-archive` document.
struct ArchiveDoc<'a> {
    program: &'static str,
    service: &'static str,
    queries: EmitQueries<'a>,
    subject: EmitSubject<'a>,
    algorithm_options: &'a [Param],
    program_options: Vec<Param>,
    step_number: Option<i32>,
    alignments: Vec<V<'a>>,
    masks: Vec<EmitMask<'a>>,
    ka_blocks: [KaBlock; 2],
    search_space: i64,
    length_adjustment: i64,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1063-1138
/// `CRemoteBlast::ConvertToRemoteMasks` of the per-query mask vectors for `program`.
fn convert_to_remote_masks<'a>(
    per_query: Vec<Vec<MaskRange<'a>>>,
    program: EProgram,
) -> Vec<EmitMask<'a>> {
    let translated = is_translated_query(program);
    let default_frame = if program.query_is_nucleotide() || translated {
        ("plus1", 1)
    } else {
        ("notset", 0)
    };
    let mut out = Vec::new();
    for query_masks in per_query {
        if query_masks.is_empty() {
            continue;
        }
        let mut current_frame = query_masks[0].frame;
        let mut run = Vec::new();
        for range in query_masks {
            if translated && range.frame != current_frame {
                if !run.is_empty() {
                    let (frame_name, frame_value) = frame_number_to_network(current_frame);
                    out.push(EmitMask {
                        ranges: std::mem::take(&mut run),
                        frame_name,
                        frame_value,
                    });
                }
                current_frame = range.frame;
            }
            run.push(range);
        }
        if !run.is_empty() {
            let (frame_name, frame_value) = if translated {
                frame_number_to_network(current_frame)
            } else {
                default_frame
            };
            out.push(EmitMask {
                ranges: run,
                frame_name,
                frame_value,
            });
        }
    }
    out
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1965-1988
fn frame_number_to_network(frame: i32) -> (&'static str, i32) {
    match frame {
        1 => ("plus1", 1),
        2 => ("plus2", 2),
        3 => ("plus3", 3),
        -1 => ("minus1", 4),
        -2 => ("minus2", 5),
        -3 => ("minus3", 6),
        _ => ("notset", 0),
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:374-385
// NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:477-497
/// `program-options` for a database search: positive `TaxidList` (sorted, deduplicated) or an
/// empty `NegativeTaxidList` (matching `blast_formatter.cpp:381-385`, which constructs the
/// negative list from `taxid_list`).
fn build_program_options(search_options: &BuiltOptions, has_database: bool) -> Vec<Param> {
    if !has_database {
        return Vec::new();
    }
    let mut out = Vec::new();
    if let Some(list) = search_options
        .taxid_list
        .as_ref()
        .filter(|list| !list.is_empty())
    {
        let sorted: BTreeSet<i32> = list.iter().copied().collect();
        out.push(Param::new(
            fields::TAXID_LIST,
            ParamValue::IntegerList(sorted.into_iter().collect()),
        ));
    }
    if search_options
        .negative_taxid_list
        .as_ref()
        .is_some_and(|list| !list.is_empty())
    {
        out.push(Param::new(
            fields::NEGATIVE_TAXID_LIST,
            ParamValue::IntegerList(Vec::new()),
        ));
    }
    out
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:429-451
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1066-1133
// NCBI reference (598d8ae6): c++/src/algo/blast/format/build_archive.cpp:84-203
/// Formats one archive search object as `-outfmt 11` (`Blast4-archive`).
#[allow(clippy::too_many_arguments)]
pub(super) fn write_archive(
    search: &Search,
    search_options: &BuiltOptions,
    algo_opts: &mut AlgorithmOptions,
    database: Option<&[u8]>,
    psi_step: Option<Result<i32, ()>>,
    output: &mut Output<'_>,
    buf: &mut Vec<u8>,
    diag: &mut Diag,
) -> anyhow::Result<()> {
    if search_options.program == EProgram::Mapper {
        return Err(Refusal::Program(EProgram::Mapper.task_name()).into());
    }
    let sets = search
        .results
        .alignment_sets(&search.request.queries)
        .map_err(stop_error)?;
    let raw_masks = search
        .results
        .masks(search.request.queries.num_queries(), || {
            search.request.query_is_translated()
        })
        .map_err(stop_error)?;
    let stats = search.results.search_stats(&search.request.service);
    if let Some(matrix) = &stats.matrix {
        algo_opts.set_matrix_name(matrix);
    }

    let is_psiblast = search_options.program == EProgram::PsiBlast;
    let step_number = if is_psiblast {
        match psi_step {
            Some(Err(())) => return Err(invalid_archive()),
            Some(Ok(step)) if step != 0 && database.is_some() => Some(step),
            _ => None,
        }
    } else {
        None
    };

    let (queries, subject, service) =
        build_request_endpoints(search, is_psiblast, database, algo_opts.service_name)?;

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2243-2245
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/build_archive.cpp:95-100,124-130
    // `CBlastAncillaryData` fills the PSI Karlin blocks when `m_Task` is `"psiblast"` and the
    // regular ones otherwise; `s_BuildArchiveAll` reads the PSI blocks when the emitted request's
    // service is `"psi"` (so `psitblastn` reads null pointers and emits `-1, -1, -1`).
    let emitted_is_psi = service == "psi";
    let ka_blocks = if is_psiblast == emitted_is_psi {
        let kb = search.results.karlin_blocks();
        [
            KaBlock {
                lambda: kb.ungapped.0,
                k: kb.ungapped.1,
                h: kb.ungapped.2,
                gapped: false,
            },
            KaBlock {
                lambda: kb.gapped.0,
                k: kb.gapped.1,
                h: kb.gapped.2,
                gapped: true,
            },
        ]
    } else {
        [KaBlock::missing(false), KaBlock::missing(true)]
    };

    let all_alignments = search.results.alignments();
    let alignments: Vec<V<'_>> = sets
        .into_iter()
        .flat_map(|range| all_alignments[range].iter().copied())
        .collect();
    let masks = convert_to_remote_masks(raw_masks, search_options.program);
    let program_options = build_program_options(search_options, database.is_some());

    let doc = ArchiveDoc {
        program: algo_opts.program_name,
        service,
        queries,
        subject,
        algorithm_options: algo_opts.params(),
        program_options,
        step_number,
        alignments,
        masks,
        ka_blocks,
        search_space: stats.effective_search_space,
        length_adjustment: stats.length_adjustment,
    };

    let encoding = ArchiveEncoding::from_env();
    if encoding == ArchiveEncoding::Nothing {
        return Ok(());
    }
    output.write_object(diag, |sink, diagnostics| {
        write_archive_streamed(encoding, &doc, buf, sink, diagnostics)
    })
}

/// Builds the `queries`, `subject`, and `service` of the emitted request, validating `results.pssm`
/// and refusing `Queries::SeqLocList` where NCBI resolves locations through its data loaders.
fn build_request_endpoints<'a>(
    search: &'a Search,
    is_psiblast: bool,
    database: Option<&'a [u8]>,
    default_service: &'static str,
) -> anyhow::Result<(EmitQueries<'a>, EmitSubject<'a>, &'static str)> {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:434-446
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:587-598
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:408-420
    if is_psiblast {
        if let Some(pssm) = search.results.reply().get("pssm") {
            validate_pssm(pssm)?;
            let Some(name) = database else {
                return Err(engine_error("Empty reference for CSearchDatabase."));
            };
            if name.is_empty() {
                return Err(engine_error("Error: No database specified"));
            }
            return Ok((EmitQueries::Pssm(pssm), EmitSubject::Database(name), "psi"));
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:609-656
    let queries = match &search.request.queries {
        Queries::SeqLocList(_) => return Err(Refusal::SeqLocQuery.into()),
        Queries::BioseqSet(_) | Queries::Pssm(_) => {
            let bioseqs = search.request.queries.bioseqs();
            let any_local = bioseqs
                .iter()
                .any(|bioseq| bioseq.ids().next().is_some_and(is_local_id));
            if any_local {
                EmitQueries::Bioseqs(bioseqs)
            } else {
                let ids = bioseqs
                    .into_iter()
                    .map(|bioseq| bioseq.ids().next().expect("checked in check_queries"))
                    .collect();
                EmitQueries::SeqLocs(ids)
            }
        }
    };

    let subject = match database {
        Some(name) => {
            if name.is_empty() {
                return Err(engine_error("Error: No database specified"));
            }
            EmitSubject::Database(name)
        }
        None => EmitSubject::Sequences(search.request.subject.bioseqs()),
    };

    Ok((queries, subject, default_service))
}

/// Encodes `doc` into `enc`, spilling finalized bytes after each large element when `spill` is
/// present.
fn encode_archive<E: Encoder>(
    enc: &mut E,
    doc: &ArchiveDoc<'_>,
    spill: Option<&RefCell<Spill<'_>>>,
) {
    let flush = |enc: &mut E| {
        if let Some(spill) = spill {
            spill.borrow_mut().after_element(enc);
        }
    };
    enc.begin_document("Blast4-archive");
    enc.begin_sequence("Blast4-archive");

    encode_request(enc, doc, &flush);
    encode_results(enc, doc, &flush);

    enc.end_sequence();
    enc.end_document();
}

fn encode_request<E: Encoder>(enc: &mut E, doc: &ArchiveDoc<'_>, flush: &impl Fn(&mut E)) {
    enc.member("request", 0);
    enc.begin_sequence("Blast4-request");
    enc.member("ident", 0);
    enc.visible_str(crate::wording::LOSAT_VERSION.as_bytes());
    enc.end_member();

    enc.member("body", 1);
    enc.variant("Blast4-request-body", "queue-search", 8);
    enc.begin_sequence("Blast4-queue-search-request");
    enc.member("program", 0);
    enc.visible_str(doc.program.as_bytes());
    enc.end_member();
    enc.member("service", 1);
    enc.visible_str(doc.service.as_bytes());
    enc.end_member();

    enc.member("queries", 2);
    match &doc.queries {
        EmitQueries::Pssm(pssm) => {
            enc.variant("Blast4-queries", "pssm", 0);
            write_value(enc, pssm.ty, pssm.value);
            enc.end_variant();
            flush(enc);
        }
        EmitQueries::SeqLocs(ids) => {
            enc.variant("Blast4-queries", "seq-loc-list", 1);
            enc.begin_list(false);
            for id in ids {
                enc.element();
                enc.variant("Seq-loc", "whole", 2);
                SeqId::Db(id).write(enc);
                enc.end_variant();
                flush(enc);
            }
            enc.end_list();
            enc.end_variant();
        }
        EmitQueries::Bioseqs(bioseqs) => {
            enc.variant("Blast4-queries", "bioseq-set", 2);
            enc.begin_sequence("Bioseq-set");
            enc.member("seq-set", 7);
            enc.begin_list(false);
            for bioseq in bioseqs {
                enc.element();
                enc.variant("Seq-entry", "seq", 0);
                let v = bioseq.value();
                write_value(enc, v.ty, v.value);
                enc.end_variant();
                flush(enc);
            }
            enc.end_list();
            enc.end_member();
            enc.end_sequence();
            enc.end_variant();
        }
    }
    enc.end_member();

    enc.member("subject", 3);
    match &doc.subject {
        EmitSubject::Database(name) => {
            enc.variant("Blast4-subject", "database", 0);
            enc.visible_str(name);
            enc.end_variant();
        }
        EmitSubject::Sequences(bioseqs) => {
            enc.variant("Blast4-subject", "sequences", 1);
            enc.begin_list(false);
            for bioseq in bioseqs {
                enc.element();
                let v = bioseq.value();
                write_value(enc, v.ty, v.value);
                flush(enc);
            }
            enc.end_list();
            enc.end_variant();
        }
    }
    enc.end_member();

    encode_request_options(enc, doc);
    enc.end_sequence();
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
    enc.end_member();
}

fn encode_request_options<E: Encoder>(enc: &mut E, doc: &ArchiveDoc<'_>) {
    enc.member("algorithm-options", 5);
    enc.begin_named("Blast4-parameters");
    enc.begin_list(false);
    for param in doc.algorithm_options {
        write_param(enc, param);
    }
    enc.end_list();
    enc.end_named();
    enc.end_member();

    if !doc.program_options.is_empty() {
        enc.member("program-options", 6);
        enc.begin_named("Blast4-parameters");
        enc.begin_list(false);
        for param in &doc.program_options {
            write_param(enc, param);
        }
        enc.end_list();
        enc.end_named();
        enc.end_member();
    }

    if let Some(step) = doc.step_number {
        enc.member("format-options", 7);
        enc.begin_named("Blast4-parameters");
        enc.begin_list(false);
        write_param(
            enc,
            &Param::new(fields::WEB_STEP_NUMBER, ParamValue::Integer(step)),
        );
        enc.end_list();
        enc.end_named();
        enc.end_member();
    }
}

fn encode_results<E: Encoder>(enc: &mut E, doc: &ArchiveDoc<'_>, flush: &impl Fn(&mut E)) {
    enc.member("results", 1);
    enc.begin_sequence("Blast4-get-search-results-reply");

    enc.member("alignments", 0);
    enc.begin_named("Seq-align-set");
    enc.begin_list(true);
    for align in &doc.alignments {
        enc.element();
        write_value(enc, align.ty, align.value);
        flush(enc);
    }
    enc.end_list();
    enc.end_named();
    enc.end_member();

    enc.member("masks", 2);
    enc.begin_list(false);
    for mask in &doc.masks {
        enc.element();
        write_emit_mask(enc, mask);
        flush(enc);
    }
    enc.end_list();
    enc.end_member();

    enc.member("ka-blocks", 3);
    enc.begin_list(false);
    for block in &doc.ka_blocks {
        write_ka_block(enc, block);
    }
    enc.end_list();
    enc.end_member();

    enc.member("search-stats", 4);
    enc.begin_list(false);
    for text in [
        format!("Effective search space: {}", doc.search_space),
        format!("Effective search space used: {}", doc.search_space),
        format!("Length adjustment: {}", doc.length_adjustment),
    ] {
        enc.element();
        enc.visible_str(text.as_bytes());
    }
    enc.end_list();
    enc.end_member();

    enc.end_sequence();
    enc.end_member();
}

fn write_emit_mask<E: Encoder>(enc: &mut E, mask: &EmitMask<'_>) {
    enc.begin_sequence("Blast4-mask");
    enc.member("locations", 0);
    enc.begin_list(false);
    enc.element();
    enc.variant("Seq-loc", "packed-int", 4);
    enc.begin_named("Packed-seqint");
    enc.begin_list(false);
    for range in &mask.ranges {
        enc.element();
        enc.begin_sequence("Seq-interval");
        enc.member("from", 0);
        enc.big_int(range.from);
        enc.end_member();
        enc.member("to", 1);
        enc.big_int(range.to);
        enc.end_member();
        enc.member("id", 3);
        SeqId::Db(range.id).write(enc);
        enc.end_member();
        enc.end_sequence();
    }
    enc.end_list();
    enc.end_named();
    enc.end_variant();
    enc.end_list();
    enc.end_member();
    enc.member("frame", 1);
    enc.begin_named("Blast4-frame-type");
    enc.enumerated(mask.frame_name, mask.frame_value);
    enc.end_named();
    enc.end_member();
    enc.end_sequence();
}

/// Streams `doc` in `encoding` through the plain encoder with `Spill`; on
/// `SerialError::NonVisibleChar`, re-encodes through `Traced` (no prefix: `blast_formatter` sets
/// none) dropping the bytes already handed on.
fn write_archive_streamed(
    encoding: ArchiveEncoding,
    doc: &ArchiveDoc<'_>,
    buf: &mut Vec<u8>,
    out: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> anyhow::Result<()> {
    buf.clear();
    let (plain, held, written, io_error) = {
        let spill = RefCell::new(Spill::new(&mut *out, SPILL_LAG));
        let plain = run_plain(encoding, doc, buf, Some(&spill));
        let (held, written, io_error) = spill.into_inner().into_parts();
        (plain, held, written, io_error)
    };
    if let Some(error) = io_error {
        return Err(error.into());
    }
    match plain {
        Ok(()) => {
            out.write_all(&held)?;
            out.write_all(buf)?;
            Ok(())
        }
        Err(SerialError::NonVisibleChar(_)) => {
            let mut skip = SkipWriter {
                inner: out,
                skip: written,
            };
            write_archive_posting(encoding, doc, buf, &mut skip, diagnostics)
        }
        Err(error) => Err(unwritable_archive(error)),
    }
}

fn run_plain(
    encoding: ArchiveEncoding,
    doc: &ArchiveDoc<'_>,
    buf: &mut Vec<u8>,
    spill: Option<&RefCell<Spill<'_>>>,
) -> Result<(), SerialError> {
    match encoding {
        ArchiveEncoding::Text => {
            let mut enc = AsnTextEncoder::new(buf);
            encode_archive(&mut enc, doc, spill);
            enc.error().map_or(Ok(()), Err)
        }
        ArchiveEncoding::Binary => {
            let mut enc = AsnBinaryEncoder::new(buf);
            encode_archive(&mut enc, doc, spill);
            enc.error().map_or(Ok(()), Err)
        }
        ArchiveEncoding::Xml => {
            let mut enc = XmlEncoder::archive(buf);
            encode_archive(&mut enc, doc, spill);
            enc.error().map_or(Ok(()), Err)
        }
        ArchiveEncoding::Nothing => Ok(()),
    }
}

#[cold]
#[inline(never)]
fn write_archive_posting(
    encoding: ArchiveEncoding,
    doc: &ArchiveDoc<'_>,
    buf: &mut Vec<u8>,
    writer: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> anyhow::Result<()> {
    buf.clear();
    let post = |replaced: &ReplacedChar| {
        diagnostics.write_all(&replaced.post(None))?;
        diagnostics.flush()
    };
    match encoding {
        ArchiveEncoding::Text => write_traced_archive(
            Traced::new(AsnTextEncoder::new(buf), 0).posting(post),
            doc,
            writer,
        ),
        ArchiveEncoding::Binary => write_traced_archive(
            Traced::new(AsnBinaryEncoder::new(buf), 0).posting(post),
            doc,
            writer,
        ),
        ArchiveEncoding::Xml => write_traced_archive(
            Traced::xml(XmlEncoder::archive(buf), 0).posting(post),
            doc,
            writer,
        ),
        ArchiveEncoding::Nothing => Ok(()),
    }
}

fn write_traced_archive<E: Encoder>(
    mut enc: E,
    doc: &ArchiveDoc<'_>,
    writer: &mut dyn Write,
) -> anyhow::Result<()> {
    encode_archive(&mut enc, doc, None);
    if let Some(error) = enc.error() {
        return Err(unwritable_archive(error));
    }
    enc.flush_to(writer, true)?;
    Ok(())
}

fn unwritable_archive(error: SerialError) -> anyhow::Error {
    match error {
        SerialError::Unwritten(what) => Refusal::Unwritable(what, 11).into(),
        SerialError::OctetString => Refusal::Unwritable("an OCTET STRING", 11).into(),
        SerialError::NullElement => {
            Refusal::Undecided("an element stated as NULL in a list").into()
        }
        other => other.into(),
    }
}

#[cfg(test)]
#[path = "archive_out_tests.rs"]
mod archive_out_tests;
