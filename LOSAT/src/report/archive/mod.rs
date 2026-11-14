//! The BLAST archive (outfmt 11, `Blast4-archive`), shared by the programs.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1051-1064
//! ```c++
//! void CBlastFormat::PrintArchive(CRef<objects::CBlast4_archive> archive, CNcbiOstream& out)
//! {
//!     if (archive.Empty()) {
//!         return;
//!     }
//!     string outfmt = CNcbiEnvironment().Get("ARCHIVE_FORMAT");
//!     if (outfmt.empty()) {
//!         out << MSerial_AsnText << *archive;
//!     } else if (!NStr::CompareNocase(outfmt, "xml")) {
//!         out << MSerial_Xml << *archive;
//!     } else if (NStr::StartsWith(outfmt, "bin", NStr::eNocase)) {
//!         out << MSerial_AsnBinary << *archive;
//!     }
//! }
//! ```
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
    RequestInputs,
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
}

/// Whether NCBI's `PrintErrorArchive` writes the error archive for this `-outfmt` text.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:913-929
/// ```c++
/// void PrintErrorArchive(const CArgs & a, const list<CRef<CBlast4_error> > & msg)
/// {
///     try {
///         if(NStr::StringToInt(NStr::TruncateSpaces(a[align_format::kArgOutputFormat].AsString())) == CFormattingArgs::eArchiveFormat) {
///             ...
///             CBlastFormat::PrintArchive(archive, a[kArgOutput].AsOutputFile());
///         }
///     } catch (...) {}
/// }
/// ```
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
        // ```c++
        //         else
        //         {
        //             if(num_iters != 0)
        //             	retval.Reset(new blast::CExportStrategy(queries, opts_hndl, search_db, client_id, num_iters));
        //             else
        //             	retval.Reset(new blast::CExportStrategy(queries, opts_hndl, search_db, client_id));
        //         }
        //     ...
        //         else
        //         {
        //             retval.Reset(new blast::CExportStrategy(queries, opts_hndl,
        //             								 db_args->GetSubjects(), client_id));
        //         }
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:745-751
        // ```c++
        // void CExportStrategy::x_AddPsiNumOfIterationsToFormatOptions(unsigned int num_iters)
        // {
        // 	CRef<CBlast4_parameter> p(new CBlast4_parameter);
        // 	p->SetName(CBlast4Field::GetName(eBlastOpt_Web_StepNumber));
        //
        // 	CRef<CBlast4_value> v(new CBlast4_value);
        // 	v->SetInteger(num_iters);
        // ```
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
            && self.algorithm.iter().all(|param| match &param.value {
                options::ParamValue::Str(value) => encoding.writes_text(value.as_bytes()),
                _ => true,
            })
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:280-282,304-306
// ```c++
//             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
//
//             SaveSearchStrategy(args, m_CmdLineArgs, queries, m_OptsHndl);
//     ...
//             if (isArchiveFormat) {
//                 formatter.WriteArchive(*queries, *m_OptsHndl, *results, 0, m_Bah.GetMessages());
// ```
/// What NCBI's options handle (`m_OptsHndl`) holds for the requests of a run beyond the options
/// of the command line: the database list's parameter and whether that list kept no sequence, and
/// the parameters that the application and each batch's set-up set in it (`set_search_sets`). One
/// handle per run, which the archive (`ArchiveHandle`) and the search strategy
/// (`strategy::StrategyExport`) share, as NCBI's applications hand the one handle to
/// `SaveSearchStrategy` and `WriteArchive`; it is set once, where the database is opened
/// (`common::database::note_archive_database_filter`). What the searches set in the handle batch
/// by batch is carried by each writer (`HandleCarry`): the strategy's request of a batch is
/// written before the batch's search and its archive after it, and LOSAT writes some archives
/// after the searches of later batches, so each writer passes the batches at its own place.
#[derive(Clone, Debug, Default)]
pub struct OptionsHandle(Arc<HandleSets>);

#[derive(Debug, Default)]
struct HandleSets {
    filter: OnceLock<DatabaseFilter>,
    search: OnceLock<SearchSets>,
}

#[derive(Debug)]
struct DatabaseFilter {
    /// The `TaxidList` or `NegativeTaxidList` program option.
    list: Option<options::Param>,
    /// A list kept no sequence and `-searchsp` is not given.
    empty: bool,
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
    // ```c++
    //     Int8 total_length = BlastSeqSrcGetTotLenStats(seqsrc);
    //     if (total_length <= 0)
    //         total_length = BlastSeqSrcGetTotLen(seqsrc);
    //     Int4 num_seqs = BlastSeqSrcGetNumSeqsStats(seqsrc);
    //     if (num_seqs <= 0)
    //         num_seqs = BlastSeqSrcGetNumSeqs(seqsrc);
    //
    //     CEffectiveSearchSpaceCalculator calc(full_query_fact, *options,
    //                                          num_seqs, total_length,
    //                                          full_data->m_ScoreBlk->GetPointer());
    // ```
    /// What the database's list makes of the requests: its `program-options` parameter
    /// (`DbFilterArgs::archive_parameter`), and `empty_database` when a list kept no sequence and
    /// `-searchsp` is not given. The `EffectiveSearchSpace` option is then 0
    /// (`BLAST_CalcEffLengths` returns at once for a database length of 0, blast_setup.c:729-732),
    /// though the statistics of the results read the length of OID 0 (`SeqDb::statistics_numbers`).
    /// Set once per run, where the database is opened; a second setting is LOSAT's internal
    /// error.
    pub fn set_database_filter(
        &self,
        list: Option<options::Param>,
        empty_database: bool,
    ) -> anyhow::Result<()> {
        let set = self.0.filter.set(DatabaseFilter {
            list,
            empty: empty_database,
        });
        anyhow::ensure!(
            set.is_ok(),
            "LOSAT internal error: the database's list is set twice in one run"
        );
        Ok(())
    }

    /// The `TaxidList` or `NegativeTaxidList` program option (`set_database_filter`).
    pub fn database_list(&self) -> Option<&options::Param> {
        self.0.filter.get().and_then(|filter| filter.list.as_ref())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:257-269
    // ```c++
    //     void x_AttachValue(CRef<objects::CBlast4_parameter> p)
    //     {
    //         typedef objects::CBlast4_parameter TParam;
    //
    //         NON_CONST_ITERATE(list< CRef<TParam> >, iter, m_ReqOpts->Set()) {
    //             if ((**iter).GetName() == p->GetName()) {
    //                 (*iter) = p;
    //                 return;
    //             }
    //         }
    //
    //         m_ReqOpts->Set().push_back(p);
    //     }
    // ```
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
    // ```c++
    //     vector<Int8> eff_searchsp;
    //     for (size_t index = 0; index <= (size_t)qinfo->last_context; index++) {
    //         eff_searchsp.push_back(calc.GetEffSearchSpaceForContext(index));
    //     }
    //     options->SetEffectiveSearchSpace(eff_searchsp);
    // ```
    /// Moves `carry` past a batch: what the batch's set-up sets in the handle (`set_search_sets`'s
    /// `each_batch`), then the batch's `EffectiveSearchSpace` option when the batch was searched
    /// (`SplitQuery_SetEffectiveSearchSpace` writes it into the handle): `searched`, `None` for a
    /// batch NCBI does not search, which leaves the option as the batch before left it; a list that
    /// kept no sequence makes it 0.
    pub fn pass(&self, carry: &mut HandleCarry, searched: Option<i64>) {
        if let Some(search) = self.0.search.get() {
            for param in &search.each_batch {
                carry.0.set(param.name, param.value.clone());
            }
        }
        if let Some(space) = searched {
            let empty = self.0.filter.get().is_some_and(|filter| filter.empty);
            carry.0.set(
                "EffectiveSearchSpace",
                options::ParamValue::BigInteger(if empty { 0 } else { space }),
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
    /// is the batch's own value when the batch was searched (kept for the batches after it), else
    /// the one of the last searched batch, which is absent before one (`OptionsHandle::pass`).
    pub fn request_sets(&self, searched: Option<i64>) -> Vec<options::Param> {
        self.with_state(|state| {
            self.0.options.pass(&mut state.carry, searched);
            self.0.options.request_sets(&state.carry)
        })
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1155-1172
// ```c++
// void CBlastAppDiagHandler::Post(const SDiagMessage & mess)
// {
// 	if(m_handler != NULL) {
// 		m_handler->Post(mess);
// 	}
// 	if(m_save) {
// 		CRef<CBlast4_error> d(new CBlast4_error);
// 		string m;
// 		mess.Write(m);
// 		d->SetMessage(NStr::Sanitize(m));
// 		d->SetCode((int)mess.m_Severity);
// 		...
// 			m_messages.push_back(d);
// ```
/// Posts error 803.7 of each byte that NCBI's stream replaced in an archive of the run of the
/// application `application` (`encoder::ReplacedChar`, with its prefix: `SetDiagPostPrefix(
/// "blastn")`, `ArchiveRequest::application`) on `diagnostics`, without saving them as messages of the run's archive `archive`:
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
    // ```c++
    //     catch (const CException& e) {                                           \
    //         LOG_POST(Error << "Error: " << e.what());                           \
    //         exit_code = BLAST_UNKNOWN_ERROR;                                    \
    //     }                                                                       \
    // ```
    // `Error: NCBI C++ Exception:` is that report of an exception of NCBI's C++ Toolkit
    // (BLASTN's `-window_masker_db` file that cannot be opened).
    const PREFIXES: [&str; 5] = [
        "BLAST query error:",
        "BLAST Database error:",
        "BLAST engine error:",
        "BLAST query/options error:",
        "Error: NCBI C++ Exception:",
    ];
    let known = PREFIXES
        .iter()
        .any(|prefix| native.message.starts_with(prefix))
        || (subject_file && native.message.starts_with("Command line argument error:"));
    if !known {
        return None;
    }
    // One diagnostic per `ERR_POST`: the text of an error is one message even when it has
    // several lines (an exception's chain), and the macro of an invalid input posts the line
    // "Please refer to the BLAST+ user manual." as a second one.
    const REFER: &str = "Please refer to the BLAST+ user manual.";
    let text = native.message.trim_end_matches('\n');
    let (first, refer) = match text.strip_suffix(REFER) {
        Some(head) if head.ends_with('\n') => (head, Some(REFER)),
        _ => (text, None),
    };
    Some(
        std::iter::once(first)
            .chain(refer)
            .filter_map(|part| {
                let text = crate::report::diag_log::sanitize(part.as_bytes());
                (!text.is_empty()).then_some((ERROR, text))
            })
            .collect(),
    )
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:327-332
// ```c++
//     } CATCH_ALL(status)
//
//     if(!m_Bah.GetMessages().empty()) {
//     	const CArgs & a = GetArgs();
//     	PrintErrorArchive(a, m_Bah.GetMessages());
//     }
// ```
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
mod tests {
    use super::*;
    use crate::blastinput::app::Program;

    /// `writes_text` and the encoders decide the same bytes: for every byte and each stream, a
    /// one-byte VisibleString that `writes_text` accepts is written by the default encoder
    /// without an error and by the replacing one with nothing replaced, and one that it refuses
    /// is the default encoder's `NonVisibleChar` and the replacing one's single replaced byte.
    #[test]
    fn writes_text_is_what_each_encoder_writes_as_it_is() {
        use crate::report::serial::asn_ber::AsnBinaryEncoder;
        use crate::report::serial::asn_text::AsnTextEncoder;
        use crate::report::serial::encoder::{Encoder, SerialError};
        use crate::report::serial::xml::XmlEncoder;
        for byte in 0..=u8::MAX {
            for encoding in [
                ArchiveEncoding::Text,
                ArchiveEncoding::Binary,
                ArchiveEncoding::Xml,
            ] {
                let as_is = encoding.writes_text(&[byte]);
                let (error, replaced) = match encoding {
                    ArchiveEncoding::Text => {
                        let mut buf = Vec::new();
                        let mut plain = AsnTextEncoder::new(&mut buf);
                        plain.visible_str(&[byte]);
                        let error = plain.error();
                        let mut buf = Vec::new();
                        let mut replacing = AsnTextEncoder::new(&mut buf);
                        replacing.visible_str_replacing(&[byte]);
                        (error, replacing.take_replaced())
                    }
                    ArchiveEncoding::Binary => {
                        let mut buf = Vec::new();
                        let mut plain = AsnBinaryEncoder::new(&mut buf);
                        plain.visible_str(&[byte]);
                        let error = plain.error();
                        let mut buf = Vec::new();
                        let mut replacing = AsnBinaryEncoder::new(&mut buf);
                        replacing.visible_str_replacing(&[byte]);
                        (error, replacing.take_replaced())
                    }
                    _ => {
                        let mut buf = Vec::new();
                        let mut plain = XmlEncoder::archive(&mut buf);
                        plain.visible_str(&[byte]);
                        let error = plain.error();
                        let mut buf = Vec::new();
                        let mut replacing = XmlEncoder::archive(&mut buf);
                        replacing.visible_str_replacing(&[byte]);
                        (error, replacing.take_replaced())
                    }
                };
                if as_is {
                    assert_eq!(error, None, "{encoding:?} 0x{byte:02X}");
                    assert!(replaced.is_empty(), "{encoding:?} 0x{byte:02X}");
                } else {
                    assert_eq!(
                        error,
                        Some(SerialError::NonVisibleChar(byte)),
                        "{encoding:?} 0x{byte:02X}"
                    );
                    assert_eq!(replaced.len(), 1, "{encoding:?} 0x{byte:02X}");
                    assert_eq!(replaced[0].byte, byte);
                }
            }
        }
    }

    #[test]
    fn archive_format_is_read_like_print_archive() {
        let of = |value: &str| ArchiveEncoding::of(Some(OsStr::new(value)));
        assert_eq!(ArchiveEncoding::of(None), ArchiveEncoding::Text);
        assert_eq!(of(""), ArchiveEncoding::Text);
        assert_eq!(of("xml"), ArchiveEncoding::Xml);
        assert_eq!(of("XML"), ArchiveEncoding::Xml);
        assert_eq!(of("Xml"), ArchiveEncoding::Xml);
        assert_eq!(of("bin"), ArchiveEncoding::Binary);
        assert_eq!(of("BIN"), ArchiveEncoding::Binary);
        assert_eq!(of("binary"), ArchiveEncoding::Binary);
        assert_eq!(of("bin foo"), ArchiveEncoding::Binary);
        assert_eq!(of(" bin"), ArchiveEncoding::Nothing);
        assert_eq!(of("asn"), ArchiveEncoding::Nothing);
        assert_eq!(of("text"), ArchiveEncoding::Nothing);
        assert_eq!(of("foo"), ArchiveEncoding::Nothing);
        assert_eq!(of("xmlx"), ArchiveEncoding::Nothing);
        assert_eq!(of("bi"), ArchiveEncoding::Nothing);
    }

    /// The strings the archive copies from the command line are written as they are where the
    /// encoding does: text and binary write `#` for a byte outside `0x20..=0x7E`, XML writes
    /// every byte but a control one as it is (`writes_text`); such a request is written through
    /// `Traced` (`writer::write_archive`).
    #[test]
    fn a_request_tells_whether_its_strings_are_written_as_they_are() {
        let as_is = |encoding, path: Option<&str>, db: Option<&str>| {
            let mut inputs = OptionInputs::new(Program::Blastn, "megablast");
            inputs.window_masker_db = path.map(str::to_string);
            ArchiveRequest::for_inputs(&inputs, db).written_as_is(encoding)
        };
        for encoding in [
            ArchiveEncoding::Text,
            ArchiveEncoding::Binary,
            ArchiveEncoding::Xml,
        ] {
            assert!(as_is(encoding, None, None));
            assert!(as_is(encoding, Some("wm&<\"'q/mg.obinary"), Some("d b&<")));
            // A control byte: NCBI writes `#` and a post, but TAB, LF and CR in XML.
            for bad in ["a\x01b", "a\x1fb", "a\x0bb", "a\x0eb"] {
                assert!(!as_is(encoding, Some(bad), None), "{encoding:?} {bad:?}");
                assert!(!as_is(encoding, None, Some(bad)), "{encoding:?} {bad:?}");
            }
        }
        for good in ["a\tb", "a\nb", "a\rb"] {
            assert!(as_is(ArchiveEncoding::Xml, Some(good), Some(good)));
            assert!(!as_is(ArchiveEncoding::Text, Some(good), None));
        }
        // Bytes from 0x7F up: replaced in text and binary, written as they are in XML.
        for high in ["wm\u{e9}/mg.obinary", "d\u{7f}b", "\u{1f600}"] {
            assert!(
                as_is(ArchiveEncoding::Xml, Some(high), Some(high)),
                "{high:?}"
            );
            for encoding in [ArchiveEncoding::Text, ArchiveEncoding::Binary] {
                assert!(!as_is(encoding, Some(high), None), "{encoding:?} {high:?}");
                assert!(!as_is(encoding, None, Some(high)), "{encoding:?} {high:?}");
            }
        }
    }

    /// A record's title and Seq-ids: in text and binary a byte outside `0x20..=0x7E`, in XML a
    /// control byte (which no FASTA title or id holds: NCBI's reader ends the title there).
    #[test]
    fn a_bioseq_tells_whether_its_strings_are_written_as_they_are() {
        use requests::RequestRecord;
        let caf = FastaRecord::new("Query_1", "caf\u{e9}".as_bytes(), b"ACGT");
        let del = FastaRecord::new("Query_1", b"ti\x7ftle", b"ACGT");
        let plain = FastaRecord::new("Query_1", b"title", b"ACGT");
        for encoding in [ArchiveEncoding::Text, ArchiveEncoding::Binary] {
            assert!(!caf.bioseq(false).written_as_is(encoding));
            assert!(!del.bioseq(false).written_as_is(encoding));
            assert!(plain.bioseq(false).written_as_is(encoding));
        }
        assert!(caf.bioseq(false).written_as_is(ArchiveEncoding::Xml));
        assert!(del.bioseq(false).written_as_is(ArchiveEncoding::Xml));
        // A parsed local id with DEL (NCBI's `CSeq_id::CheckLocalID` takes it).
        let ids = [crate::seqid::ParsedSeqId::Local(
            crate::seqid::model::ObjectId::Str(b"a\x7fb".to_vec()),
        )];
        let bioseq = bioseq::BioseqInput {
            ids: Some(&ids),
            ..plain.bioseq(false)
        };
        assert!(!bioseq.written_as_is(ArchiveEncoding::Text));
        assert!(bioseq.written_as_is(ArchiveEncoding::Xml));
    }

    #[test]
    fn the_search_space_option_carries_over_unsearched_batches() {
        let handle = ArchiveHandle::new(
            ArchiveEncoding::Text,
            ArchiveRequest {
                program: "blastp",
                application: "blastp",
                service: "plain",
                algorithm: Vec::new(),
                database: None,
                step_number: None,
            },
            OptionsHandle::default(),
        );
        let space = |searched| search_space(&handle.request_sets(searched));
        assert_eq!(space(None), None);
        assert_eq!(space(Some(8738)), Some(8738));
        assert_eq!(space(None), Some(8738));
        assert_eq!(space(Some(0)), Some(0));
        assert_eq!(space(None), Some(0));
    }

    /// The `EffectiveSearchSpace` option among a request's sets.
    pub(super) fn search_space(sets: &[options::Param]) -> Option<i64> {
        sets.iter()
            .find(|param| param.name == "EffectiveSearchSpace")
            .map(|param| match param.value {
                options::ParamValue::BigInteger(space) => space,
                ref other => panic!("EffectiveSearchSpace {other:?}"),
            })
    }

    /// The application's parameters come first and each batch's set-up's before its search
    /// space, by `x_AttachValue` (rpsblast's documents after the first, NCBI BLAST+ 2.17.0
    /// `p9.rps.3v.q40.n1.f6`, `p9.rps.err.nofreq`).
    #[test]
    fn the_handle_sets_fold_in_ncbis_order() {
        use options::{Param, ParamValue};
        let handle = OptionsHandle::default();
        handle
            .set_search_sets(
                vec![Param {
                    name: "CompositionBasedStats",
                    value: ParamValue::Integer(0),
                }],
                vec![
                    Param {
                        name: "DbLength",
                        value: ParamValue::BigInteger(560),
                    },
                    Param {
                        name: "MatrixName",
                        value: ParamValue::Str("BLOSUM62".to_string()),
                    },
                ],
            )
            .unwrap();
        assert!(handle.set_search_sets(Vec::new(), Vec::new()).is_err());
        let mut carry = HandleCarry::default();
        let names = |carry: &HandleCarry| -> Vec<&'static str> {
            handle
                .request_sets(carry)
                .iter()
                .map(|param| param.name)
                .collect()
        };
        assert_eq!(names(&carry), ["CompositionBasedStats"]);
        handle.pass(&mut carry, None);
        assert_eq!(
            names(&carry),
            ["CompositionBasedStats", "DbLength", "MatrixName"]
        );
        handle.pass(&mut carry, Some(85_680));
        handle.pass(&mut carry, Some(96_628));
        assert_eq!(
            names(&carry),
            [
                "CompositionBasedStats",
                "DbLength",
                "MatrixName",
                "EffectiveSearchSpace"
            ]
        );
        assert_eq!(search_space(&handle.request_sets(&carry)), Some(96_628));
    }

    #[test]
    fn the_tee_saves_lines_for_the_next_archive() {
        let handle = ArchiveHandle::new(
            ArchiveEncoding::Text,
            ArchiveRequest {
                program: "blastn",
                application: "blastn",
                service: "plain",
                algorithm: Vec::new(),
                database: None,
                step_number: None,
            },
            OptionsHandle::default(),
        );
        let mut seen: Vec<u8> = Vec::new();
        {
            let mut tee = DiagTee::new(&mut seen, handle.clone());
            tee.write_all(b"Warning: [blastn] a  b\n").unwrap();
        }
        assert_eq!(seen, b"Warning: [blastn] a  b\n");
        assert_eq!(
            handle.take_messages(),
            vec![(1, b"Warning: [blastn] a b".to_vec())]
        );
        assert!(handle.take_messages().is_empty());
    }

    /// NCBI 2.17.0 (`PrintErrorArchive`): each text below was run with a missing `-query` and
    /// wrote the error archive (`true`) or nothing (`false`).
    #[test]
    fn the_error_archive_follows_the_outfmt_text_like_print_error_archive() {
        for (text, writes) in [
            ("11", true),
            (" 11", true),
            ("11 ", true),
            ("\t11\n", true),
            ("\x0b11\x0c", true),
            ("011", true),
            ("+11", true),
            ("00000000000011", true),
            ("11 qseqid", false),
            (" 11 qseqid", false),
            ("11qseqid", false),
            ("11.0", false),
            ("1 1", false),
            ("-11", false),
            ("0x0B", false),
            ("1_1", false),
            ("11,0", false),
            ("4294967307", false),
            ("2147483648", false),
            ("", false),
            (" ", false),
            ("+", false),
            ("6", false),
        ] {
            assert_eq!(writes_error_archive(text), writes, "{text:?}");
        }
    }

    #[test]
    fn a_handle_tied_to_a_spec_writes_no_error_archive() {
        let make = || {
            ArchiveHandle::new(
                ArchiveEncoding::Text,
                ArchiveRequest {
                    program: "tblastx",
                    application: "tblastx",
                    service: "plain",
                    algorithm: Vec::new(),
                    database: None,
                    step_number: None,
                },
                OptionsHandle::default(),
            )
        };
        let lines = vec![(ERROR, b"Error: x".to_vec())];
        // No text (the API): the error archive is written.
        let mut sink = Vec::new();
        write_error_archive(&make(), &lines, &mut sink).unwrap();
        assert!(!sink.is_empty());
        // `11`: written; `11 qseqid`: nothing, and the messages stay unread.
        let mut sink = Vec::new();
        write_error_archive(&make().with_outfmt_text(" 11 "), &lines, &mut sink).unwrap();
        assert!(!sink.is_empty());
        let mut sink = Vec::new();
        let spec = make().with_outfmt_text("11 qseqid");
        assert!(!spec.error_archive_enabled());
        write_error_archive(&spec, &lines, &mut sink).unwrap();
        assert!(sink.is_empty());
    }
}
