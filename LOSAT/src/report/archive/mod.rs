//! The BLAST archive (outfmt 11, `Blast4-archive`), shared by the programs.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1051-1064
//! One `Blast4-archive` per query batch is written after the batch's search, one error
//! archive when the application stops with a diagnostic (`PrintErrorArchive`). Each archive
//! is written straight to the sink, element by element (`writer`), from data the program
//! holds for its other formats.
//!
//! - `options`: the request's `algorithm-options` (the fold of NCBI's `Set*` calls);
//! - `bioseq`: the Bioseq of a query or subject record;
//! - `masks`: the `Blast4-mask` lists (results and the `LCaseMask` options);
//! - `writer`: the document.
//!
//! What a run keeps from one archive to the next is `RunState`: the diagnostics saved since the
//! last archive (`DiagLog`) and what the searches have set in the options handle so far
//! (`HandleCarry`); what the run's database makes of every request is the run's `OptionsHandle`,
//! which the archive and the search strategy share.
//! - `requests`: the request inputs of a batch from its records.

pub mod bioseq;
pub mod masks;
pub mod options;
pub mod requests;
pub mod strategy;
pub mod writer;

#[cfg(test)]
mod option_walk_tests;

use crate::wording::{
    DATABASE_ERROR, ENGINE_ERROR, OPTIONS_ERROR, QUERY_ERROR, QUERY_OPTIONS_ERROR,
};
use std::borrow::Cow;
use std::ffi::OsStr;
use std::io::Write;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

use crate::blastinput::fasta_reader::FastaRecord;
use crate::report::diag_log::{DiagLog, Message, ERROR, WARNING};
use crate::report::serial::encoder::Stream;
use crate::stats::tables::KarlinParams;
use options::OptionInputs;

pub use writer::{
    write_archive, write_error_archive, AlignSource, ArchiveOut, BatchArchive, KaBlock,
    RequestInputs, SearchedContext0,
};

/// The Karlin-Altschul block of the results (`Blast4-ka-block`) that `params` makes.
pub(crate) fn ka_block(params: &KarlinParams, gapped: bool) -> KaBlock {
    KaBlock {
        lambda: params.lambda,
        k: params.k,
        h: params.h,
        gapped,
    }
}

/// The defline line of the Bioseq: the line as read, else `>` and the title.
pub(crate) fn defline_of(record: &FastaRecord) -> Cow<'_, [u8]> {
    match record.archive_defline() {
        Some(defline) => Cow::Borrowed(defline),
        None => {
            let mut line = vec![b'>'];
            line.extend_from_slice(&record.title);
            Cow::Owned(line)
        }
    }
}

/// How `ARCHIVE_FORMAT` makes NCBI write the archive.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ArchiveEncoding {
    /// ASN.1 text (the variable is unset or empty).
    Text,
    /// ASN.1 binary (a value that starts with `bin`, in any case).
    Binary,
    /// XML (the value `xml` in any case): NCBI's XML object stream (`XmlEncoder::archive`).
    Xml,
    /// Any other value: NCBI writes nothing.
    Nothing,
}

impl ArchiveEncoding {
    /// `PrintArchive`'s reading of `ARCHIVE_FORMAT`.
    pub fn of(value: Option<&OsStr>) -> Self {
        let Some(value) = value else {
            return ArchiveEncoding::Text;
        };
        let value = value.to_string_lossy();
        if value.is_empty() {
            ArchiveEncoding::Text
        } else if value.eq_ignore_ascii_case("xml") {
            ArchiveEncoding::Xml
        } else if value.len() >= 3 && value.as_bytes()[..3].eq_ignore_ascii_case(b"bin") {
            ArchiveEncoding::Binary
        } else {
            ArchiveEncoding::Nothing
        }
    }

    /// The encoding the environment asks for.
    pub fn from_env() -> Self {
        Self::of(std::env::var_os("ARCHIVE_FORMAT").as_deref())
    }

    /// Whether this encoding's stream writes `text` into a VisibleString of the archive without
    /// replacing a byte by `#` with NCBI's post 803.7 (`encoder::ReplacedChar`): the stream's own
    /// rule (`serial::encoder::Stream::replaces`). Text and binary replace any byte outside
    /// `0x20..=0x7E` (`GoodVisibleChar`). XML: the archive's stream has no encoding, so it writes
    /// every other byte as it is (TAB, LF and CR as character references) and replaces
    /// `0x01..=0x1F` without those three (`x_VerifyChar`, `BAD_CHAR`; objostrxml.cpp:548-565).
    /// `Nothing` writes no archive; it reads the text rule.
    pub fn writes_text(self, text: &[u8]) -> bool {
        let stream = match self {
            ArchiveEncoding::Xml => Stream::Xml,
            ArchiveEncoding::Binary => Stream::AsnBinary,
            ArchiveEncoding::Text | ArchiveEncoding::Nothing => Stream::AsnText,
        };
        !text.iter().any(|&byte| stream.replaces(byte))
    }

    /// Whether the string values of `params` are written as they are (`writes_text`).
    pub fn writes_params(self, params: &[options::Param]) -> bool {
        params.iter().all(|param| match &param.value {
            options::ParamValue::Str(value) => self.writes_text(value.as_bytes()),
            _ => true,
        })
    }
}

/// Whether NCBI's `PrintErrorArchive` writes the error archive for this `-outfmt` text.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:913-929
/// The conversion throws (and nothing is written) for text after the number, such as
/// `11 qseqid`; the batches' archives accept such a text (`-outfmt "11 <spec>"`), only the
/// error archive does not. `TruncateSpaces` (`eTrunc_Both`) and `StringToInt` (flags 0) are
/// `ncbi_str::truncate_spaces` and `ncbi_str::string_to_int`.
pub fn writes_error_archive(outfmt: &str) -> bool {
    use crate::ncbi_str::{string_to_int, truncate_spaces};
    string_to_int(truncate_spaces(outfmt.as_bytes()), 0)
        .is_ok_and(|value| value == crate::blastinput::app::output_format::ARCHIVE_FORMAT)
}

/// The `queue-search` request of a run, as the options handle and the arguments make it.
#[derive(Clone, Debug)]
pub struct ArchiveRequest {
    /// `blastn`, `blastp`, ...
    pub program: &'static str,
    /// The application that writes the request (`options::RequestKind::application`): the
    /// prefix of NCBI's posts of the bytes its stream replaces in it (`psiblast` for psiblast,
    /// whose `program` is `blastp`).
    pub application: &'static str,
    /// `megablast`, `psi`, `rpsblast` or `plain` (`options::service`).
    pub service: &'static str,
    /// The `algorithm-options` of the command line (what the search sets after them is the
    /// writer's `HandleCarry`).
    pub algorithm: Vec<options::Param>,
    /// `subject database "<-db as typed>"`; `None` for `-subject` (the sequences are written
    /// from the records).
    pub database: Option<Vec<u8>>,
    /// The `Web_StepNumber` format option: psiblast's `-num_iterations` with `-db`, when it is not
    /// 0.
    pub step_number: Option<i32>,
}

impl ArchiveRequest {
    /// The request of a run of `inputs.kind`: `database` is the `-db` argument as typed (`None`
    /// for `-subject`).
    pub(crate) fn for_inputs(inputs: &OptionInputs, database: Option<&str>) -> Self {
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:289-295,305-309
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:745-751
        // psiblast's `-num_iterations` (1 by default) with `-db` (`p9.ps.q.db.itdef.f6`), none
        // for 0 or with `-subject` (`p9.ps.q.db.it0.f6`, `p9.ps.q.sub.it2.f6`).
        let step_number = inputs
            .num_iterations
            .filter(|&iterations| iterations != 0 && database.is_some());
        Self {
            program: inputs.kind.program(),
            application: inputs.kind.application(),
            service: options::service(inputs.kind, &inputs.task),
            algorithm: options::algorithm_options(inputs),
            database: database.map(|name| name.as_bytes().to_vec()),
            step_number,
        }
    }

    /// Whether `encoding` writes the strings that the request copies from the command line as
    /// they are (`ArchiveEncoding::writes_text`): the `-db` name and the string options (the
    /// `-window_masker_db` path, the `-index_name`, the matrix name).
    pub(crate) fn written_as_is(&self, encoding: ArchiveEncoding) -> bool {
        self.database
            .as_deref()
            .is_none_or(|name| encoding.writes_text(name))
            && encoding.writes_params(&self.algorithm)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:164-173 (quoted
// at `OptionsHandle::set_database_filter`) and c++/src/algo/blast/api/effsearchspace_calc.cpp:111-121
// The request's `EffectiveSearchSpace` is `BLAST_CalcEffLengths` with the statistics' length and
// count before the setup's fallback to OID 0 (`BLAST_GapAlignSetUp`, blast_setup.c:929-937): the
// search's own value where the length is positive; nothing (0) for a length of 0 (blast_setup.c:
// 729-732); for a negative one (an alias file's `LENGTH` of 2^63 or more), no length adjustment
// (`BLAST_ComputeLengthAdjustment`'s `c < 0`, whatever the Karlin block) and an effective database
// length of 1, so the context's query length (NCBI 2.17.0: blastp 142, blastn 577 with
// `LENGTH 9223372036854775809`, probe run2 of line RG D1a-3, impl1).
/// How the `EffectiveSearchSpace` option of a searched batch's request is made.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum RequestLength {
    /// The search's own value of the batch's context 0: a positive statistics length, or
    /// `-searchsp`.
    #[default]
    Search,
    /// `BLAST_CalcEffLengths` of context 0 with the database numbers that
    /// `stats::setup::eff_lengths_database` gives for the statistics' length before the fallback
    /// (`None`: it returns at once, the option is 0), a length that is not positive
    /// (`note_archive_database_filter` refuses `-dbsize` there).
    Lengths(Option<(i64, i32)>),
}

impl RequestLength {
    /// The `EffectiveSearchSpace` option of a searched batch whose context 0 is `context0`.
    pub fn space(self, context0: writer::SearchedContext0) -> i64 {
        match self {
            Self::Search => context0.space,
            Self::Lengths(None) => 0,
            Self::Lengths(Some(_)) if context0.query_length <= 0 => 0,
            Self::Lengths(Some((db_length, num_seqs))) => {
                debug_assert!(db_length <= 0, "a request length that is not positive");
                // The Karlin block is not read for a database length that is not positive.
                crate::stats::setup::calc_eff_lengths(
                    context0.query_length,
                    db_length,
                    num_seqs,
                    &crate::stats::KarlinParams::default(),
                    0,
                    "",
                )
                .map_or(0, |lengths| lengths.eff_searchsp)
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:280-282,304-306
/// What NCBI's options handle (`m_OptsHndl`) holds for the requests of a run beyond the options
/// of the command line: the database list's parameter and how the database's statistics make the
/// `EffectiveSearchSpace` option of a searched batch (`RequestLength`), and the parameters that
/// the application and each batch's set-up set in it (`set_search_sets`). One
/// handle per run, which the archive (`ArchiveHandle`) and the search strategy
/// (`strategy::StrategyExport`) share, as NCBI's applications hand the one handle to
/// `SaveSearchStrategy` and `WriteArchive`; it is set once, where the database is opened
/// (`common::database::note_archive_database_filter`). What the searches set in the handle batch
/// by batch is carried by each writer (`HandleCarry`): the strategy's request of a batch is
/// written before the batch's search and its archive after it, and LOSAT writes some archives
/// after the searches of later batches, so each writer passes the batches at its own place.
#[derive(Clone, Debug, Default)]
pub struct OptionsHandle(Arc<HandleSets>);

/// What `OptionsHandle` holds: the database list's parameter, set where the database is opened;
/// the application's search parameters (`set_search_sets`); and the database masking's, set where
/// the formatter resolves it (after `Query is Empty!`).
#[derive(Debug, Default)]
struct HandleSets {
    filter: OnceLock<DatabaseFilter>,
    search: OnceLock<SearchSets>,
    /// `DbFilteringAlgorithmId` or `DbFilteringAlgorithmKey`, then `SubjectMaskingType`: the
    /// archive's (`blastinput::db_mask::ResolvedMask::archive_params`) and the search strategy's
    /// (`DbMaskOption::strategy_params`).
    masking: OnceLock<DatabaseMasking>,
}

/// The database masking's program options of the archive and of the search strategy, which
/// differ after the formatter's warning (`DbMaskOption::strategy_params`).
#[derive(Debug)]
struct DatabaseMasking {
    archive: Vec<options::Param>,
    strategy: Vec<options::Param>,
}

#[derive(Debug)]
struct DatabaseFilter {
    /// The `TaxidList` or `NegativeTaxidList` program option.
    list: Option<options::Param>,
    /// How the `EffectiveSearchSpace` option of a searched batch is made.
    request_length: RequestLength,
}

/// The parameters that a run sets in the handle besides the command line's and the batches'
/// `EffectiveSearchSpace` (`OptionsHandle::set_search_sets`).
#[derive(Debug)]
struct SearchSets {
    /// Set by the application once, before the first batch.
    before_batches: Vec<options::Param>,
    /// Set by each batch's set-up, before its search space.
    each_batch: Vec<options::Param>,
}

/// What the searches of the batches that one writer of a run's requests (the archive, the search
/// strategy) has passed set in the options handle (`OptionsHandle::pass`): the parameters in the
/// order NCBI's `x_AttachValue` keeps them (appended when first set, replaced in place), which the
/// writer applies by the same rule onto the request's options (`OptionsHandle::request_sets`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HandleCarry(options::ParamList);

impl OptionsHandle {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:164-173
    /// What the database's list and statistics make of the requests: its `program-options`
    /// parameter (`DbFilterArgs::archive_parameter`), and how the `EffectiveSearchSpace` option
    /// of a searched batch is made (`RequestLength`, `common::database::note_archive_database_filter`):
    /// for a statistics length that is not positive it is not the search's own value, though the
    /// statistics of the results read the length of OID 0 (`SeqDb::statistics`). Set once per
    /// run, where the database is opened; a second setting is LOSAT's internal error.
    pub fn set_database_filter(
        &self,
        list: Option<options::Param>,
        request_length: RequestLength,
    ) -> anyhow::Result<()> {
        let set = self.0.filter.set(DatabaseFilter {
            list,
            request_length,
        });
        anyhow::ensure!(
            set.is_ok(),
            "LOSAT internal error: the database's list is set twice in one run"
        );
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:499-504
    /// The database masking's program options of the archive (`ResolvedMask::archive_params`)
    /// and of the search strategy (`DbMaskOption::strategy_params`), which follow the list's: set
    /// once per run, where the formatter resolves the masking; a second setting is LOSAT's
    /// internal error.
    pub fn set_database_masking(
        &self,
        archive: Vec<options::Param>,
        strategy: Vec<options::Param>,
    ) -> anyhow::Result<()> {
        anyhow::ensure!(
            self.0
                .masking
                .set(DatabaseMasking { archive, strategy })
                .is_ok(),
            "LOSAT internal error: the database masking is set twice in one run"
        );
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:257-269
    /// The parameters that the run sets in the handle after the command line's options and
    /// besides the batches' `EffectiveSearchSpace`: `before_batches` once, before the first batch
    /// (rpsblast's `CompositionBasedStats` 0 for a database without frequency ratios), and
    /// `each_batch` at each batch's set-up (the RPS database's length and `.aux` costs,
    /// `rps_search::volumes::search_sets`). Set once per run, where the database is opened; a
    /// second setting is LOSAT's internal error.
    pub fn set_search_sets(
        &self,
        before_batches: Vec<options::Param>,
        each_batch: Vec<options::Param>,
    ) -> anyhow::Result<()> {
        let set = self.0.search.set(SearchSets {
            before_batches,
            each_batch,
        });
        anyhow::ensure!(
            set.is_ok(),
            "LOSAT internal error: the search's options are set twice in one run"
        );
        Ok(())
    }

    /// The archive's program options of the database: the `TaxidList` or `NegativeTaxidList`
    /// option (`set_database_filter`), then the masking's (`set_database_masking`).
    pub fn database_params(&self) -> Vec<options::Param> {
        self.params_with(|masking| &masking.archive)
    }

    /// `database_params` of the search strategy, whose masking options are its own.
    pub fn strategy_database_params(&self) -> Vec<options::Param> {
        self.params_with(|masking| &masking.strategy)
    }

    fn params_with(
        &self,
        masking: impl Fn(&DatabaseMasking) -> &Vec<options::Param>,
    ) -> Vec<options::Param> {
        self.0
            .filter
            .get()
            .and_then(|filter| filter.list.clone())
            .into_iter()
            .chain(
                self.0
                    .masking
                    .get()
                    .map(masking)
                    .into_iter()
                    .flatten()
                    .cloned(),
            )
            .collect()
    }

    /// What a request holds of the handle beyond the command line's options, after the batches
    /// that `carry` has passed: the application's parameters set before the batches, then the
    /// carried ones, folded by `x_AttachValue`. The writer applies them onto the request's options
    /// by the same rule (`writer::write_request_body`).
    pub fn request_sets(&self, carry: &HandleCarry) -> Vec<options::Param> {
        let mut sets = options::ParamList::default();
        if let Some(search) = self.0.search.get() {
            for param in &search.before_batches {
                sets.set(param.name, param.value.clone());
            }
        }
        for param in carry.0.params() {
            sets.set(param.name, param.value.clone());
        }
        sets.into_params()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:177-181
    /// Moves `carry` past a batch: what the batch's set-up sets in the handle (`set_search_sets`'s
    /// `each_batch`), then the batch's `EffectiveSearchSpace` option when the batch was searched
    /// (`SplitQuery_SetEffectiveSearchSpace` writes it into the handle): `searched` is the batch's
    /// context 0 (the option is the search's own value, or what `RequestLength` makes of it),
    /// `None` for a batch NCBI does not search, which leaves the option as the batch before left
    /// it.
    pub fn pass(&self, carry: &mut HandleCarry, searched: Option<SearchedContext0>) {
        if let Some(search) = self.0.search.get() {
            for param in &search.each_batch {
                carry.0.set(param.name, param.value.clone());
            }
        }
        if let Some(context0) = searched {
            let length = self
                .0
                .filter
                .get()
                .map_or(RequestLength::Search, |filter| filter.request_length);
            carry.0.set(
                "EffectiveSearchSpace",
                options::ParamValue::BigInteger(length.space(context0)),
            );
        }
    }

    /// Whether `self` and `other` are the same handle.
    pub fn same(&self, other: &OptionsHandle) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

#[derive(Clone, Debug, Default)]
struct RunState {
    log: DiagLog,
    /// What the searches set in the options handle, after the archives written so far.
    carry: HandleCarry,
    /// Whether the diagnostics that are written now are not the archive's (the warnings that
    /// another format of the same search writes, which a run of the archive alone does not).
    muted: bool,
}

#[derive(Default)]
struct Inner {
    real: RunState,
}

impl Inner {
    fn state(&mut self) -> &mut RunState {
        &mut self.real
    }
}

struct Shared {
    encoding: ArchiveEncoding,
    request: ArchiveRequest,
    /// The run's options handle (shared with the search strategy, `OptionsHandle`).
    options: OptionsHandle,
    /// Whether the error archives are written (`writes_error_archive` of the `-outfmt` text;
    /// true when the run has no such text).
    error_archive: AtomicBool,
    state: Mutex<Inner>,
}

/// What an archive-writing run shares between its diagnostics stream and its writer.
#[derive(Clone)]
pub struct ArchiveHandle(Arc<Shared>);

impl std::fmt::Debug for ArchiveHandle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ArchiveHandle")
            .field("encoding", &self.0.encoding)
            .finish_non_exhaustive()
    }
}

impl ArchiveHandle {
    /// The handle of a run of `inputs.kind` that writes the archive (`-outfmt 11`): the
    /// request from the recorded options, `ARCHIVE_FORMAT` read once. `database` is the
    /// `-db` argument as typed (`None` for `-subject`); `options` is the run's options handle
    /// (made once per run; its search strategy shares it).
    pub fn for_program(
        inputs: &OptionInputs,
        database: Option<&str>,
        options: &OptionsHandle,
    ) -> Self {
        Self::for_encoding(ArchiveEncoding::from_env(), inputs, database, options)
    }

    /// `for_program` for a given encoding (`ARCHIVE_FORMAT` read by the caller).
    pub fn for_encoding(
        encoding: ArchiveEncoding,
        inputs: &OptionInputs,
        database: Option<&str>,
        options: &OptionsHandle,
    ) -> Self {
        Self::new(
            encoding,
            ArchiveRequest::for_inputs(inputs, database),
            options.clone(),
        )
    }

    /// The handle of a run's archive, with the run's options handle `options`.
    pub fn new(encoding: ArchiveEncoding, request: ArchiveRequest, options: OptionsHandle) -> Self {
        Self(Arc::new(Shared {
            encoding,
            request,
            options,
            error_archive: AtomicBool::new(true),
            state: Mutex::new(Inner::default()),
        }))
    }

    /// The run's options handle, which the run's search strategy shares
    /// (`strategy::StrategyExport::for_run`).
    pub fn options(&self) -> &OptionsHandle {
        &self.0.options
    }

    pub fn encoding(&self) -> ArchiveEncoding {
        self.0.encoding
    }

    /// Ties the error archives to the `-outfmt` text of the command line
    /// (`writes_error_archive`): NCBI's `PrintErrorArchive` writes none for `11 <spec>`.
    /// The archives of the batches are not affected.
    #[must_use]
    pub fn with_outfmt_text(self, outfmt: &str) -> Self {
        self.0
            .error_archive
            .store(writes_error_archive(outfmt), Ordering::Relaxed);
        self
    }

    /// Whether `write_error_archive` writes (see `with_outfmt_text`).
    pub fn error_archive_enabled(&self) -> bool {
        self.0.error_archive.load(Ordering::Relaxed)
    }

    pub fn request(&self) -> &ArchiveRequest {
        &self.0.request
    }

    fn with_state<R>(&self, f: impl FnOnce(&mut RunState) -> R) -> R {
        let mut inner = self
            .0
            .state
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        f(inner.state())
    }

    /// Saves the lines that a write to the diagnostics stream completes.
    pub fn feed(&self, bytes: &[u8]) {
        self.with_state(|state| {
            if !state.muted {
                state.log.feed(bytes);
            }
        });
    }

    /// While muted, what is written to the diagnostics stream is not saved. The report of a
    /// format other than the archive is written with the handle muted: a search with several
    /// formats writes the warnings of the first one once, and the archive of the same search
    /// alone has none of them.
    pub fn set_muted(&self, muted: bool) {
        self.with_state(|state| state.muted = muted);
    }

    /// Saves a message (the error lines of an aborted run).
    pub fn post(&self, severity: i32, line: &[u8]) {
        self.with_state(|state| state.log.post(severity, line));
    }

    /// Writes `text`, one warning whose text has several lines (an `ERR_POST` of an
    /// exception's `what()`), to `diagnostics`, and saves it as one message, as NCBI's handler
    /// saves one `SDiagMessage` (`diag_log`); the diagnostics stream (`DiagTee`) would save
    /// each line. Nothing is saved while the handle is muted.
    pub fn write_warning(&self, diagnostics: &mut dyn Write, text: &[u8]) -> std::io::Result<()> {
        self.write_post(diagnostics, WARNING, text)
    }

    /// `write_warning` of a post of `severity` (`diag_log::Posts`).
    pub fn write_post(
        &self,
        diagnostics: &mut dyn Write,
        severity: i32,
        text: &[u8],
    ) -> std::io::Result<()> {
        let muted = self.with_state(|state| std::mem::replace(&mut state.muted, true));
        let written = diagnostics.write_all(text);
        self.with_state(|state| {
            state.muted = muted;
            if !muted {
                state.log.post(severity, text);
            }
        });
        written
    }

    /// Writes `text` to `diagnostics` without saving it: NCBI's posts while it writes an archive,
    /// which its `ResetMessages` after the archive drops (blastn_app.cpp:305-306).
    pub fn write_unsaved(&self, diagnostics: &mut dyn Write, text: &[u8]) -> std::io::Result<()> {
        let muted = self.with_state(|state| std::mem::replace(&mut state.muted, true));
        let written = diagnostics.write_all(text);
        self.with_state(|state| state.muted = muted);
        written
    }

    /// Whether messages are saved since the last archive.
    pub fn has_messages(&self) -> bool {
        self.with_state(|state| !state.log.is_empty())
    }

    /// The messages saved since the last archive.
    pub fn take_messages(&self) -> Vec<Message> {
        self.with_state(|state| state.log.take())
    }

    /// What the next archive's request holds of the options handle beyond the command line's
    /// options, after its batch (`OptionsHandle::request_sets`): its `EffectiveSearchSpace` option
    /// is the batch's own value when the batch was searched (`searched`, its context 0; kept for
    /// the batches after it), else the one of the last searched batch, which is absent before one
    /// (`OptionsHandle::pass`).
    pub fn request_sets(&self, searched: Option<writer::SearchedContext0>) -> Vec<options::Param> {
        self.with_state(|state| {
            self.0.options.pass(&mut state.carry, searched);
            self.0.options.request_sets(&state.carry)
        })
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1155-1172
/// Posts error 803.7 of each byte that NCBI's stream replaced in an archive of the run of the
/// application `application` (`encoder::ReplacedChar`, with its prefix: `SetDiagPostPrefix`
/// of the application name, such as "blastn"; `ArchiveRequest::application`) on `diagnostics`,
/// without saving them as messages of the run's archive `archive`:
/// NCBI's `ResetMessages` after the archive drops them (`ArchiveHandle::write_unsaved`). The
/// search strategy's posts are saved (`strategy::StrategyExport::take_posts`).
pub(crate) fn post_replaced(
    diagnostics: &mut dyn Write,
    archive: &ArchiveHandle,
    application: &str,
    replaced: &[crate::report::serial::encoder::ReplacedChar],
) -> std::io::Result<()> {
    for byte in replaced {
        archive.write_unsaved(diagnostics, &byte.post(Some(application)))?;
    }
    Ok(())
}

/// An encoder that writes nothing and finds whether `encoding` writes every VisibleString it is
/// given as it is (`ArchiveEncoding::writes_text`): the Seq-ids of a record are checked with the
/// writer that writes them (`report::seqid::SeqId::write`).
struct StringsAsIs {
    encoding: ArchiveEncoding,
    as_is: bool,
}

impl StringsAsIs {
    fn new(encoding: ArchiveEncoding) -> Self {
        Self {
            encoding,
            as_is: true,
        }
    }
}

impl crate::report::serial::encoder::Encoder for StringsAsIs {
    fn begin_document(&mut self, _: &str) {}
    fn end_document(&mut self) {}
    fn begin_sequence(&mut self, _: &str) {}
    fn end_sequence(&mut self) {}
    fn member(&mut self, _: &str, _: u8) {}
    fn end_member(&mut self) {}
    fn begin_list(&mut self, _: bool) {}
    fn element(&mut self) {}
    fn end_list(&mut self) {}
    fn variant(&mut self, _: &str, _: &str, _: u8) {}
    fn end_variant(&mut self) {}
    fn int(&mut self, _: i32) {}
    fn big_int(&mut self, _: i64) {}
    fn enumerated(&mut self, _: &str, _: i32) {}
    fn named_int(&mut self, _: &str, _: i32) {}
    fn boolean(&mut self, _: bool) {}
    fn visible_str(&mut self, value: &[u8]) {
        self.as_is &= self.encoding.writes_text(value);
    }
    fn real(&mut self, _: f64) {}
    fn octets(&mut self, _: &[u8]) {}
    fn null(&mut self) {}
    fn null_element(&mut self) {}
    fn utf8_str(&mut self, _: &[u8]) {}
    fn delayed(&mut self, _: &[u8]) {}
    fn error(&self) -> Option<crate::report::serial::encoder::SerialError> {
        None
    }
}

/// Whether `encoding` writes the strings of the Seq-id `id` as they are.
pub(crate) fn seq_id_written_as_is(
    encoding: ArchiveEncoding,
    id: crate::report::seqid::SeqId<'_>,
) -> bool {
    let mut check = StringsAsIs::new(encoding);
    id.write(&mut check);
    check.as_is
}

/// The error lines of a run that stopped with `error`, when NCBI writes an error archive for
/// it (`PrintErrorArchive` after its `CATCH_ALL`): the errors of the application's `Run`,
/// i.e. NCBI's own error texts, one message per diagnostic, severity `eDiag_Error`. LOSAT's own
/// refusals (and the argument errors raised before the application runs) write none.
/// `subject_file` is whether the error came before the output was opened (the taxid or seqid
/// list, opening the `-query` file, reading the `-subject` file); NCBI raises the "Command line
/// argument error" of those files inside `Run`.
pub fn error_lines(error: &anyhow::Error, subject_file: bool) -> Option<Vec<Message>> {
    let native = error.downcast_ref::<crate::cli::NativeError>()?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
    // `Error: NCBI C++ Exception:` is that report of an exception of NCBI's C++ Toolkit
    // (BLASTN's `-window_masker_db` file that cannot be opened).
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-220
    // `BLAST options error:` is that report of an options error that the search raises inside
    // `Run`: BLASTN's database mask that the database lacks (`DbMaskOption::validate_at_open`,
    // review 1 F1 of M4a-2a) and PSIBLAST's PSSM options (`IPssmInputData`; LOSAT's PSIBLAST
    // refuses `-outfmt 11`).
    // (An exception report, `Error: ` and NCBI's `what()` in LOSAT's words, has exit status 255:
    // `catch_all` below, `cli::NativeError::exception`.)
    const PREFIXES: [&[u8]; 5] = [
        QUERY_ERROR.as_bytes(),
        DATABASE_ERROR.as_bytes(),
        ENGINE_ERROR.as_bytes(),
        OPTIONS_ERROR.as_bytes(),
        QUERY_OPTIONS_ERROR.as_bytes(),
    ];
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:260-263
    // `CATCH_ALL`'s two reports of an exception, a `CException`'s and a `std::exception`'s
    // (`NativeError::exception`, `NativeError::std_exception`: exit status 255, `Error: `): the
    // latter is `Error: vector::reserve` of the `-in_pssm` reader's Dense-seg hook (NCBI BLAST+
    // 2.17.0 tblastn -outfmt 11, line RG D3-1 `e09_densereserve`). After NCBI's object manager
    // failed the main thread's diagnostic stream, that report is its line feed alone
    // (`cli::catch_all_text`).
    let catch_all = native.exit == 255
        && (native.message.starts_with(b"Error: ") || crate::ncbi_exception::main_diag_failbit());
    let known = catch_all
        || PREFIXES
            .iter()
            .any(|prefix| native.message.starts_with(prefix))
        || (subject_file && native.message.starts_with(b"Command line argument error:"));
    if !known {
        return None;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1155-1172 (cited in
    // `report::diag_log`: `Write` of `mess` for each post)
    // In a failed stream every post writes its header only (`ncbi_exception::post_line`) and saves
    // that text: `CATCH_ALL`'s posts, which have no header, save empty messages, one per line
    // feed that `cli::catch_all_text` writes (NCBI BLAST+ 2.17.0 tblastn -outfmt 11 of line RG
    // D3-1 `k02_failed_then_choice_exc`: `message ""`).
    if crate::ncbi_exception::main_diag_failbit() {
        let written = crate::cli::catch_all_text(&native.message);
        // One post per line feed (`split_terminator('\n')`: a final empty piece is not one).
        let pieces = if written.is_empty() {
            0
        } else {
            written
                .strip_suffix(b"\n")
                .unwrap_or(written)
                .split(|&byte| byte == b'\n')
                .count()
        };
        return Some((0..pieces).map(|_| (ERROR, Vec::new())).collect());
    }
    // One diagnostic per `ERR_POST`: the text of an error is one message even when it has
    // several lines (an exception's chain), and the macro of an invalid input posts the user
    // manual footer (`USER_MANUAL_FOOTER`, `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 12) as a
    // second one.
    let refer_line = crate::wording::USER_MANUAL_FOOTER.as_bytes();
    let text = &native.message[..native
        .message
        .iter()
        .rposition(|&byte| byte != b'\n')
        .map_or(0, |last| last + 1)];
    let (first, refer) = match text.strip_suffix(refer_line) {
        Some(head) if head.ends_with(b"\n") => (head, Some(refer_line)),
        _ => (text, None),
    };
    Some(
        std::iter::once(first)
            .chain(refer)
            .filter_map(|part| {
                let text = crate::report::diag_log::sanitize(part);
                (!text.is_empty()).then_some((ERROR, text))
            })
            .collect(),
    )
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:351-355
/// The error lines of a run that `CATCH_ALL` ended without a post (`cli::NativeError::silent_end`:
/// the failed write of the object manager's message after its diagnostic stream failed, a PSSM
/// file that ends after its header): none, and NCBI writes the error archive of the messages saved
/// before when there are some (`PrintErrorArchive` whatever the status). The command line's and
/// the library's ends of a run both read it after `error_lines`.
pub fn silent_end_lines(handle: &ArchiveHandle, error: &anyhow::Error) -> Option<Vec<Message>> {
    (crate::cli::is_silent_end(error) && handle.has_messages()).then(Vec::new)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:327-332
/// Writes the error archive of a command line run that stopped with `error` (`error_lines` of it:
/// `lines`) into `sink`, after the error's text. NCBI's `CATCH_ALL` writes the text to standard
/// error before `PrintErrorArchive` writes the archive (`cli::post_native_error`), so in one
/// stream of both (`2>&1`) the text comes between what the run wrote before into `sink` (the
/// archives of the batches before, flushed first) and the error archive. Without an error archive
/// (`11 <spec>`: `writes_error_archive`) the text is left to `main`, after the run's output; a
/// library caller of the command line entries (no `cli::post_errors_in_run`) gets the archive
/// alone and the text in the error, as before.
pub fn write_error_archive_after_post(
    handle: &ArchiveHandle,
    error: &anyhow::Error,
    lines: &[Message],
    sink: &mut dyn Write,
) -> anyhow::Result<()> {
    if handle.error_archive_enabled() {
        let _ = sink.flush();
        crate::cli::post_native_error(error);
    }
    write_error_archive(handle, lines, sink)
}

/// A diagnostics stream that saves what it writes as the messages of the next archive.
pub struct DiagTee<'a> {
    inner: &'a mut (dyn Write + Send),
    handle: ArchiveHandle,
}

impl<'a> DiagTee<'a> {
    pub fn new(inner: &'a mut (dyn Write + Send), handle: ArchiveHandle) -> Self {
        Self { inner, handle }
    }
}

impl Write for DiagTee<'_> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let written = self.inner.write(buf)?;
        self.handle.feed(&buf[..written]);
        Ok(written)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
