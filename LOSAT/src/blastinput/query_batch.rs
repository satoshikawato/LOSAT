//! Query batches of a search (NCBI `CBlastInput::GetNextSeqBatch`).

use crate::blastinput::fasta_reader::{FastaRecord, ReadError};

/// `NStr::StringToInt` with the default flags (`ncbi_str::string_to_int`): an optional sign and
/// decimal digits, within the range of an `int`, and nothing else (no spaces; the digits end at
/// a NUL, as NCBI's `while (char ch = str[pos])` ends there); `None` where NCBI throws its
/// `CStringException`.
pub fn ncbi_string_to_int(value: &std::ffi::OsStr) -> Option<i32> {
    crate::ncbi_str::string_to_int(value.as_encoded_bytes(), 0).ok()
}

/// The query batch size of blastn from the environment variable `BATCH_SIZE` (0 when it
/// is not set: NCBI then uses `CBatchSizeMixer`). `Err` holds the value that NCBI cannot
/// convert (an empty value too).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input_aux.cpp:85-91
/// ```c
///     // used for experimentation purposes
///     char* batch_sz_str = getenv("BATCH_SIZE");
///     if (batch_sz_str) {
///         retval = NStr::StringToInt(batch_sz_str);
///         _TRACE("DEBUG: Using query batch size " << retval);
///         return retval;
///     }
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input_aux.cpp:99
/// ```c
///     if (! use_default) return 0;
/// ```
pub fn get_query_batch_size(value: Option<&std::ffi::OsStr>) -> Result<i32, String> {
    match value {
        None => Ok(0),
        Some(value) => {
            ncbi_string_to_int(value).ok_or_else(|| value.to_string_lossy().into_owned())
        }
    }
}

/// NCBI's adaptive query batch size of blastn: the first batch aims at 1/200 of the target
/// hits, and each later batch at the target over the (mixed) ratio of successful ungapped
/// extensions per residue of the batches before.
///
/// NCBI reference: c++/src/app/blast/blast_app_util.hpp:57-83
/// ```c++
/// class CBatchSizeMixer
/// {
/// private:
///     const double k_MixIn;        // mixing factor between batches
///     Int4 m_TargetHits;           // the target hits per batch
///     double m_Ratio;              // the hits to batch size ratio
///     Int4 m_BatchSize;            // the batch size for next run
///     const Int4 k_MinBatchSize;   // the minimum allowable batch size
///     const Int4 k_MaxBatchSize;   // the maximum allowable batch size
///     const Int4 k_MinTargetHits;  // the minimum target hits per batch
///
/// public:
///     CBatchSizeMixer(Int4 max_batch_size)
///           : k_MixIn        (0.3),
///             m_TargetHits   (1000000),
///             m_Ratio        (-1.0),
///             m_BatchSize    (5000),
///             k_MinBatchSize (100),
///             k_MaxBatchSize (max_batch_size),
/// 	    k_MinTargetHits(1000000) { }
///
///     void SetTargetHits(Int4 target) {
///         m_TargetHits = (target < k_MinTargetHits) ? k_MinTargetHits : target;
///         m_BatchSize = m_TargetHits / 200;
///     }
/// ```
#[derive(Clone, Debug)]
pub struct BatchSizeMixer {
    target_hits: i32,
    ratio: f64,
    batch_size: i32,
    max_batch_size: i32,
}

impl BatchSizeMixer {
    const MIX_IN: f64 = 0.3;
    const MIN_BATCH_SIZE: i32 = 100;
    const MIN_TARGET_HITS: i32 = 1_000_000;

    pub fn new(max_batch_size: i32) -> Self {
        Self {
            target_hits: Self::MIN_TARGET_HITS,
            ratio: -1.0,
            batch_size: 5000,
            max_batch_size,
        }
    }

    pub fn set_target_hits(&mut self, target: i32) {
        self.target_hits = target.max(Self::MIN_TARGET_HITS);
        self.batch_size = self.target_hits / 200;
    }

    /// The next batch size; `hits` is the number of successful ungapped extensions of the
    /// batch just searched (`None` before the first batch).
    ///
    /// NCBI reference: c++/src/app/blast/blast_app_util.cpp:65-81
    /// ```c++
    /// Int4 CBatchSizeMixer::GetBatchSize(Int4 hits)
    /// {
    ///      if (hits >= 0) {
    ///          double ratio = 1.0 * (hits+1) / m_BatchSize;
    ///          m_Ratio = (m_Ratio < 0) ? ratio
    ///                  : k_MixIn * ratio + (1.0 - k_MixIn) * m_Ratio;
    ///          m_BatchSize = (Int4) (1.0 * m_TargetHits / m_Ratio);
    ///      }
    ///      if (m_BatchSize > k_MaxBatchSize) {
    ///          m_BatchSize = k_MaxBatchSize;
    ///          m_Ratio = -1.0;
    ///      } else if (m_BatchSize < k_MinBatchSize) {
    ///          m_BatchSize = k_MinBatchSize;
    ///          m_Ratio = -1.0;
    ///      }
    ///      return m_BatchSize;
    /// }
    /// ```
    /// The `(Int4)` conversion is the x86-64 one of the oracle (`INT_MIN` beyond `Int4`), so a
    /// batch with fewer than about 2.3 extensions per million target hits and residue is
    /// followed by the smallest batch.
    pub fn batch_size(&mut self, hits: Option<i32>) -> i32 {
        if let Some(hits) = hits.filter(|&hits| hits >= 0) {
            let ratio = f64::from(hits.wrapping_add(1)) / f64::from(self.batch_size);
            self.ratio = if self.ratio < 0.0 {
                ratio
            } else {
                Self::MIX_IN * ratio + (1.0 - Self::MIX_IN) * self.ratio
            };
            self.batch_size = crate::core::blast_util::ncbi_int4_from_double(
                f64::from(self.target_hits) / self.ratio,
            );
        }
        if self.batch_size > self.max_batch_size {
            self.batch_size = self.max_batch_size;
            self.ratio = -1.0;
        } else if self.batch_size < Self::MIN_BATCH_SIZE {
            self.batch_size = Self::MIN_BATCH_SIZE;
            self.ratio = -1.0;
        }
        self.batch_size
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2295-2313
// ```c++
// void CBlastFormat::ResetScopeHistory()
// {
//     // Do not reset the scope for BLAST2Sequences or else we'll loose the
//     // sequence data! (see x_CreateSubjectBioseq)
//     if (m_IsBl2Seq){
//         return;
//     }
//     ...
//     if (m_FormatType != CFormattingArgs::eXml)
//     {
//         m_Scope->ResetDataAndHistory();
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-258 (each program's batch loop)
// ```c++
//         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
// ```
/// The query records that NCBI's scope holds while it reports a query batch under
/// `-parse_deflines` (`QueryScopes`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum QueryScopeSpan {
    /// The batch's records: the batch reader adds them to the scope, and the formatter empties
    /// it after the batch's reports.
    Batch,
    /// Every batch's: the formatter keeps them with `-subject` (bl2seq) and in outfmt 5, whose
    /// XML document it writes after the last batch (NCBI BLAST+ 2.17.0 then shows a reported
    /// database sequence with the id of a later batch's query as that query, or fails:
    /// `Q/SURVEY_Q1c.md` of line RG).
    Run,
}

impl QueryScopeSpan {
    /// The span of a search with `-subject` (`subject`) or a database, whose reports `xml`
    /// include outfmt 5 (where one call writes several formats, outfmt 5 among them gives the
    /// larger span).
    pub(crate) fn of(subject: bool, xml: bool) -> Self {
        if subject || xml {
            QueryScopeSpan::Run
        } else {
            QueryScopeSpan::Batch
        }
    }
}

/// Whether one of the `-outfmt` values `outfmts` is outfmt 5 (`QueryScopeSpan::of`).
pub(crate) fn outfmts_have_xml<'a>(outfmts: impl IntoIterator<Item = &'a str>) -> bool {
    use crate::blastinput::app::{parse_formatting_string, ReportFormat};
    outfmts.into_iter().any(|outfmt| {
        parse_formatting_string(outfmt)
            .is_ok_and(|choice| ReportFormat::from_number(choice.number) == Some(ReportFormat::Xml))
    })
}

/// The ids of the query records of one scope (`QueryScopes`): the keys of each id
/// (`seqid::defline_set::scope_keys`, and the own kind of `id_scope::version_keys` and
/// `id_scope::name_keys`).
#[derive(Debug, Default)]
pub(crate) struct QueryScope {
    keys: std::collections::HashSet<u64>,
}

impl QueryScope {
    /// Adds the ids `ids` of a query record, or LOSAT's refusal of a record whose id another
    /// record of the scope has (`believed::duplicate_id`).
    fn add(&mut self, ids: &[crate::seqid::ParsedSeqId], program: &str) -> anyhow::Result<()> {
        if let Some(refusal) = crate::blastinput::fasta_reader::believed::duplicate_id(
            ids,
            &mut self.keys,
            "query",
            program,
        ) {
            return Err(refusal);
        }
        // `duplicate_id` keeps the version's kind of two BLASTN queries only.
        self.keys.extend(
            ids.iter()
                .filter_map(|id| crate::seqid::id_scope::version_keys(id).map(|(own, _)| own)),
        );
        Ok(())
    }

    /// Adds the ids `ids` of query records (`BelievedIds`), whatever ids the scope has
    /// (`QueryScopes::with_whole_run`): the keys that `add` keeps of an id (`scope_keys` and the
    /// own kind of `name_keys` and `version_keys`), which `has_key` asks.
    pub(crate) fn add_ids_of<'r>(
        &mut self,
        ids: impl IntoIterator<Item = &'r crate::blastinput::fasta_reader::believed::BelievedIds>,
    ) {
        let mut keys = Vec::new();
        for id in ids.into_iter().flat_map(|believed| &believed.ids) {
            keys.clear();
            crate::seqid::defline_set::scope_keys(id, &mut keys);
            self.keys.extend(keys.iter().copied());
            self.keys
                .extend(crate::seqid::id_scope::name_keys(id).map(|(own, _)| own));
            self.keys
                .extend(crate::seqid::id_scope::version_keys(id).map(|(own, _)| own));
        }
    }

    /// Whether a query of the scope has an id with one of the keys `keys`.
    pub(crate) fn has_key(&self, keys: &[u64]) -> bool {
        keys.iter().any(|key| self.keys.contains(key))
    }
}

/// Whether the reports of a search show the ids of its query records `records`: under
/// `-parse_deflines` (`parse_deflines`), and with `-in_pssm` (`pssm_query`) when the PSSM's query is
/// shown with its ids (`pssm_search::PsiSearch::query_record`: an id that is not one local id).
/// The scopes of such queries (`QueryScopes::of_search`) refuse what NCBI's scope mixes up.
pub(crate) fn believed_queries(
    parse_deflines: bool,
    pssm_query: bool,
    records: &[FastaRecord],
) -> bool {
    parse_deflines || (pssm_query && records.iter().any(|record| record.believed().is_some()))
}

/// The scopes of the query records of a search under `-parse_deflines`, by searched query:
/// NCBI's scope finds a sequence by its id, so two queries of one scope with one id show as one
/// (the first: its letters are searched for both), and a reported database sequence with the id
/// of a query of its scope shows as that query (`DI7/DECISIONS.md` 8, 12 of the drop-in survey;
/// `common::database::refuse_id_of_a_query`). LOSAT refuses both.
#[derive(Debug)]
pub(crate) struct QueryScopes {
    span: QueryScopeSpan,
    /// The program, for the refusals.
    program: &'static str,
    scopes: Vec<QueryScope>,
    /// The first searched query of each scope.
    starts: Vec<usize>,
    /// The ids of every query of the run, for the reported database sequences where they come
    /// before the batches add theirs (`with_whole_run`).
    whole_run: Option<QueryScope>,
    /// Whether a read before the search added every batch (`read_ahead`).
    read_ahead: bool,
}

impl QueryScopes {
    /// No scope yet: the batches add theirs (`add_batch`).
    pub(crate) fn new(span: QueryScopeSpan, program: &'static str) -> Self {
        Self {
            span,
            program,
            scopes: Vec::new(),
            starts: Vec::new(),
            whole_run: None,
            read_ahead: false,
        }
    }

    /// The scopes that a read of the queries before the search made with the search's batches
    /// (BLASTX's, whose reader gives the batches one at a time): the batches of the search add
    /// nothing again (`add_batch_ids`), and a refusal came before the search.
    pub(crate) fn read_ahead(mut self) -> Self {
        self.read_ahead = true;
        self
    }

    /// With `QueryScopeSpan::Run`, `whole_run`: the ids of every query of the run (`QueryScope::
    /// add_ids_of`), which `of_query` gives. A reader that gives the batches one at a time
    /// (BLASTX's) collects them before the search, where NCBI writes outfmt 5 after the last batch:
    /// its scope then holds every batch's queries for the reported database sequences of the first.
    /// The batches still add theirs, for two queries with one id.
    pub(crate) fn with_whole_run(mut self, whole_run: QueryScope) -> Self {
        if self.span == QueryScopeSpan::Run {
            self.whole_run = Some(whole_run);
        }
        self
    }

    /// The span of the scopes.
    pub(crate) fn span(&self) -> QueryScopeSpan {
        self.span
    }

    /// `new` for a search whose reports show the ids of its queries (`believed`, `believed_queries`;
    /// `None` otherwise), with `QueryScopeSpan::of` of `subject` and `xml`.
    pub(crate) fn of_search(
        believed: bool,
        subject: bool,
        xml: bool,
        program: &'static str,
    ) -> Option<Self> {
        believed.then(|| Self::new(QueryScopeSpan::of(subject, xml), program))
    }

    /// The scopes of `records`, the input records of one call of a search whose queries
    /// (searched from the first, `-query_loc` cuts included) they all are: one batch.
    pub(crate) fn of_records(
        span: QueryScopeSpan,
        program: &'static str,
        records: &[FastaRecord],
    ) -> anyhow::Result<Self> {
        let mut scopes = Self::new(span, program);
        scopes.add_batch(records, 0)?;
        Ok(scopes)
    }

    /// Adds the batch of the input records `records`, whose first searched query is
    /// `first_searched` (a new scope, or with `Run` the one scope), or LOSAT's refusal of a
    /// record with the id of another of the scope.
    pub(crate) fn add_batch(
        &mut self,
        records: &[FastaRecord],
        first_searched: usize,
    ) -> anyhow::Result<()> {
        self.add_batch_ids(
            records.iter().filter_map(FastaRecord::believed),
            first_searched,
        )
    }

    /// `add_batch` of the ids `ids` of a batch's records (`BelievedIds`; BLASTX's reader has
    /// records of its own).
    pub(crate) fn add_batch_ids<'r>(
        &mut self,
        ids: impl IntoIterator<Item = &'r crate::blastinput::fasta_reader::believed::BelievedIds>,
        first_searched: usize,
    ) -> anyhow::Result<()> {
        if self.read_ahead {
            return Ok(());
        }
        if self.scopes.is_empty() || self.span == QueryScopeSpan::Batch {
            self.scopes.push(QueryScope::default());
            self.starts.push(first_searched);
        }
        let scope = self.scopes.last_mut().expect("pushed above");
        for believed in ids {
            scope.add(&believed.ids, self.program)?;
        }
        Ok(())
    }

    /// `with_whole_run` of the ids of the input records `records` (every query of the run), where
    /// NCBI's scope holds them when the reports read the database sequences (outfmt 5) and the
    /// batches are cut as the search goes (blastn's, tblastx's loops).
    pub(crate) fn with_whole_run_of(self, records: &[FastaRecord]) -> Self {
        if self.span != QueryScopeSpan::Run {
            return self;
        }
        let mut whole_run = QueryScope::default();
        whole_run.add_ids_of(records.iter().filter_map(FastaRecord::believed));
        self.with_whole_run(whole_run)
    }

    /// The scope of the searched query `query` (of the batches added; `whole_run` if set).
    pub(crate) fn of_query(&self, query: usize) -> &QueryScope {
        if let Some(whole_run) = &self.whole_run {
            return whole_run;
        }
        let after = self.starts.partition_point(|&start| start <= query);
        &self.scopes[after.saturating_sub(1)]
    }

    /// The scopes of the searched queries from `first`, as one call that searches them numbers
    /// them (`of_query` of its index).
    pub(crate) fn numbered_from(&self, first: usize) -> QueryScopesFrom<'_> {
        QueryScopesFrom {
            scopes: self,
            first,
        }
    }
}

/// `QueryScopes` of the searched queries from `first` (`QueryScopes::numbered_from`).
#[derive(Clone, Copy, Debug)]
pub(crate) struct QueryScopesFrom<'a> {
    scopes: &'a QueryScopes,
    first: usize,
}

impl QueryScopesFrom<'_> {
    /// The scope of the call's query `query`.
    pub(crate) fn of_query(&self, query: usize) -> &QueryScope {
        self.scopes.of_query(self.first + query)
    }
}

/// The query batches as NCBI's batch loop reads them, up to the batch that stops the run.
pub(crate) struct QueryBatches {
    /// The reader's messages of each batch's records (skipped ones too), before the
    /// report of the batch's first searched query, by searched query.
    pub(crate) warnings: Vec<Vec<u8>>,
    /// The input records and the searched queries of the batches that NCBI searches.
    pub(crate) input_done: usize,
    pub(crate) searched_done: usize,
    /// Where each of those batches starts: (first input record, first searched query).
    pub(crate) starts: Vec<(usize, usize)>,
    /// What stops the run at the batch after them: the messages that NCBI writes as it
    /// reads that batch, and its error.
    pub(crate) stop: Option<(Vec<u8>, anyhow::Error)>,
    /// The searched queries of that batch when NCBI's set-up of its search fails (`batch_error`:
    /// every query without letters): NCBI made the batch's query factory, so its search strategy
    /// request is written (`SaveSearchStrategy` before `CLocalBlast::Run`, blastp_app.cpp:259-276).
    /// `None` when the reading or an empty batch stops the run (no request).
    pub(crate) stop_batch: Option<std::ops::Range<usize>>,
    /// Under `-parse_deflines`, the scopes of those batches' queries (`QueryScopes`).
    pub(crate) query_scopes: Option<QueryScopes>,
    /// LOSAT's refusal of two queries of one scope with one id, in the first batch that has them
    /// before any batch that NCBI fails (`QueryScopes::add_batch`): the run returns it before any
    /// report, as a refusal follows no written byte; the other fields are then those of the
    /// batches before it.
    pub(crate) query_refusal: Option<anyhow::Error>,
}

impl QueryBatches {
    /// Where the batch of the searched query `query` starts: (first input record, first
    /// searched query). Every batch in `starts` holds at least one searched query.
    pub(crate) fn start_of(&self, query: usize) -> (usize, usize) {
        let after = self
            .starts
            .partition_point(|&(_, searched)| searched <= query);
        self.starts[after.saturating_sub(1)]
    }

    /// The searched queries of the batch of the searched query `query`.
    pub(crate) fn searched_span(&self, query: usize) -> std::ops::Range<usize> {
        let after = self
            .starts
            .partition_point(|&(_, searched)| searched <= query);
        let end = self
            .starts
            .get(after)
            .map_or(self.searched_done, |&(_, searched)| searched);
        self.start_of(query).1..end
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input.cpp:134-171
// ```c
// CRef<CBlastQueryVector>
// CBlastInput::GetNextSeqBatch(CScope& scope)
// {
//     CRef<CBlastQueryVector> retval(new CBlastQueryVector);
//     TSeqPos size_read = 0;
//
//     while (size_read < GetBatchSize()) {
//
//         if (End())
//             break;
//
//         CRef<CBlastSearchQuery> q;
//         try { q.Reset(m_Source->GetNextSequence(scope)); }
//         catch (const CObjReaderParseException& e) {
//             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
//                 break;
//             }
//             throw;
//         }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:378-380
// ```c
//     if (queries.Empty()) {
//         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
//     }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652
// ```c
//         } catch (const CException& e) {
//             ...
//             CRef<CSearchMessage> m
//                 (new CSearchMessage(eBlastSevWarning, index, e.GetMsg()));
//             messages[index].push_back(m);
//             s_InvalidateQueryContexts(qinfo, index);
//         }
//     ...
//     // Validate that at least one query context is valid
//     if (BlastSetup_Validate(qinfo, NULL) != 0 && messages.HasMessages()) {
//         NCBI_THROW(CBlastException, eSetup, messages.ToString());
//     }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_aux.cpp:1013-1025
// ```c
// TSearchMessages::ToString() const
// {
//     string retval;
//     ITERATE(vector<TQueryMessages>, qm, *this) {
//         if (qm->empty()) {
//             continue;
//         }
//         ITERATE(TQueryMessages, msg, *qm) {
//             retval += (*msg)->GetMessage() + " ";
//         }
//     }
//     return retval;
// }
// ```
/// Reads the query batches as NCBI's loop reads them (`input` and `searched` are the
/// records read and the queries searched, `seq_range::cut_queries`): a batch counts the
/// whole records up to the batch size, and the reader's messages of its records come
/// before the report of its first searched query. The run stops at the first batch that
/// NCBI fails, after the reports of the batches before (without the epilog) and the
/// messages of that batch's records:
/// - the reading ends there with a reader error (or LOSAT's rejection of a Seq-id line);
/// - the batch has no query (a batch size of 0, a batch of records that `-query_loc`
///   skips, or a batch that an `eEOF` ends before its first record): `Empty
///   CBlastQueryVector`;
/// - every query of the batch has no letters: the context of a protein query is valid at
///   the set-up exactly when it has letters, so NCBI's set-up fails with one `Sequence
///   contains no data` per query (`GetMessage` puts the severity first) and `CATCH_ALL`
///   writes `BLAST engine error: ` (exit 3).
///
/// A batch whose queries have letters but no Karlin-Altschul parameters is set up and not
/// searched (`run_resolved_in_pool`).
///
/// Under `-parse_deflines` (`query_scopes`, without a batch yet) each batch that NCBI searches
/// adds its records to the scopes (`QueryScopes::add_batch`); LOSAT's refusal of two of one scope
/// with one id ends the reading there (`QueryBatches::query_refusal`), after the batch's errors of
/// NCBI's own (a batch that NCBI does not search reports nothing of its scope).
pub(crate) fn read_query_batches(
    input_records: &[FastaRecord],
    searched: &[FastaRecord],
    input: &crate::blastinput::seq_range::QueryInput,
    reading: QueryReading,
    batch_size: u32,
    query_scopes: Option<QueryScopes>,
) -> QueryBatches {
    read_batches(
        input_records,
        searched,
        input,
        reading,
        batch_size,
        &no_data_batch_error,
        query_scopes,
    )
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:468
// ```c++
//             input.Reset(new CBlastInput(&*fasta, 1));
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_impl.cpp:88-92
// ```c++
//     // Either PSSM or a protein query must be provided
//     if (m_Pssm.NotEmpty()) {
//         CPsiBlastValidate::Pssm(*m_Pssm);
//     } else if (m_Query.NotEmpty()) {
//         CPsiBlastValidate::QueryFactory(m_Query, *m_OptsHandle);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:636-649
// ```c++
//     string excpt_msg("PSI-BLAST only accepts ");
//     if (qf_type == eQFT_Query) {
//         excpt_msg += "one protein sequence as query";
//     } else if (qf_type == eQFT_Subject) {
//         excpt_msg += "protein sequences as subjects";
//     } else {
//         abort();
//     }
//
//     if (qf_type == eQFT_Query) {
//         if (query_data->GetNumQueries() != 1) {
//             NCBI_THROW(CBlastException, eInvalidArgument, excpt_msg);
//         }
//     }
// ```
/// psiblast's query batches (`read_query_batches` with a batch size of 1): a batch holds one
/// query with letters and the queries without letters before it. Each batch is a `CPsiBlast` of
/// its own, which takes exactly one query: a batch of several (a query without letters and the
/// next) stops the run with `BLAST engine error: PSI-BLAST only accepts one protein sequence as
/// query` (exit 3). A batch of one query without letters (the last records, or a `-query_loc`
/// interval that starts just past its record's end) crashes NCBI 2.17.0's psiblast (SIGSEGV);
/// LOSAT rejects it there. A query whose Karlin-Altschul parameters cannot be computed stops the
/// run at its batch too (`blast_engine.rs` `search`).
pub(crate) fn read_psiblast_query_batches(
    input_records: &[FastaRecord],
    searched: &[FastaRecord],
    input: &crate::blastinput::seq_range::QueryInput,
    reading: QueryReading,
    query_scopes: Option<QueryScopes>,
) -> QueryBatches {
    read_batches(
        input_records,
        searched,
        input,
        reading,
        1,
        &psiblast_batch_error,
        query_scopes,
    )
}

/// `read_query_batches`, where the batch errors after those of NCBI's batch reader (a reader
/// error, `Empty CBlastQueryVector`) are `batch_error` of the batch's searched queries.
fn read_batches(
    input_records: &[FastaRecord],
    searched: &[FastaRecord],
    input: &crate::blastinput::seq_range::QueryInput,
    mut reading: QueryReading,
    batch_size: u32,
    batch_error: &dyn Fn(&[FastaRecord]) -> Option<anyhow::Error>,
    mut query_scopes: Option<QueryScopes>,
) -> QueryBatches {
    let mut warnings: Vec<Vec<u8>> = vec![Vec::new(); searched.len()];
    let mut starts = Vec::new();
    // Whether what ended the reading after the records (`reading`) has come.
    let mut reading_ended = false;
    let mut input_start = 0;
    // The reader is at the end of its input after the last record, unless an `eEOF` or an
    // error ends the reading there; that comes in the batch being read when the records
    // run out before the batch reaches its size, and otherwise in a batch of its own.
    while input_start < input.input_lengths.len()
        || (reading.reads_past_records() && !reading_ended)
    {
        let (input_end, reached_size) = crate::blastinput::seq_range::next_ranged_batch(
            &input.input_lengths,
            &input.skipped,
            input_start,
            batch_size,
        );
        let std::ops::Range { start, end } = input.searched_in(input_start..input_end);
        let reading_ends_here =
            input_end == input.input_lengths.len() && !reached_size && reading.reads_past_records();
        // `CFastaReader` writes its messages about a batch's records (their lines and
        // titles) when it reads them, after the report of the batch before.
        let mut batch_warnings: Vec<u8> = input_records[input_start..input_end]
            .iter()
            .flat_map(|record| record.warnings.iter().copied())
            .collect();
        let mut stop_error = None;
        if reading_ends_here {
            reading_ended = true;
            if let QueryReading::Error {
                error,
                warnings: error_warnings,
            } = std::mem::replace(&mut reading, QueryReading::End)
            {
                batch_warnings.extend_from_slice(&error_warnings);
                stop_error = Some(error.into_app_error());
            }
        }
        if stop_error.is_none() && start == end {
            stop_error = Some(
                crate::cli::NativeError {
                    exit: 3,
                    message: "BLAST engine error: Empty CBlastQueryVector\n".to_string(),
                }
                .into(),
            );
        }
        let mut stop_batch = None;
        if stop_error.is_none() {
            stop_error = batch_error(&searched[start..end]);
            stop_batch = stop_error.is_some().then_some(start..end);
        }
        if let Some(error) = stop_error {
            return QueryBatches {
                warnings,
                input_done: input_start,
                searched_done: start,
                starts,
                stop: Some((batch_warnings, error)),
                stop_batch,
                query_scopes,
                query_refusal: None,
            };
        }
        // LOSAT's refusal, which the records alone decide: the caller returns it before any
        // report (`QueryBatches::query_refusal`).
        if let Some(scopes) = query_scopes.as_mut() {
            if let Err(refusal) = scopes.add_batch(&input_records[input_start..input_end], start) {
                return QueryBatches {
                    warnings,
                    input_done: input_start,
                    searched_done: start,
                    starts,
                    stop: None,
                    stop_batch: None,
                    query_scopes,
                    query_refusal: Some(refusal),
                };
            }
        }
        warnings[start].extend_from_slice(&batch_warnings);
        starts.push((input_start, start));
        input_start = input_end;
    }
    QueryBatches {
        warnings,
        input_done: input_start,
        searched_done: searched.len(),
        starts,
        stop: None,
        stop_batch: None,
        query_scopes,
        query_refusal: None,
    }
}

/// The error of a blastp or tblastn batch whose queries all have no letters: one `Sequence
/// contains no data` per query (`read_query_batches`).
fn no_data_batch_error(batch: &[FastaRecord]) -> Option<anyhow::Error> {
    if !batch.iter().all(|record| record.seq().is_empty()) {
        return None;
    }
    let mut message = String::from("BLAST engine error: ");
    for _ in batch {
        message.push_str("Warning: ");
        message.push_str(NO_DATA_MESSAGE);
        message.push(' ');
    }
    message.push('\n');
    Some(crate::cli::NativeError { exit: 3, message }.into())
}

/// The error of a psiblast batch that `CPsiBlast` does not take (`read_psiblast_query_batches`).
fn psiblast_batch_error(batch: &[FastaRecord]) -> Option<anyhow::Error> {
    match batch {
        [query] if query.seq().is_empty() => Some(anyhow::anyhow!(
            "a query without letters alone in its batch (a last record without residues, or a -query_loc interval that starts just past the end of its record), with which NCBI BLAST+'s psiblast crashes, is not supported by LOSAT's PSIBLAST"
        )),
        [_] => None,
        _ => Some(
            crate::cli::NativeError {
                exit: 3,
                message: "BLAST engine error: PSI-BLAST only accepts one protein sequence as query\n"
                    .to_string(),
            }
            .into(),
        ),
    }
}

/// How the reading of the queries of a search ended, after its records.
#[derive(Debug)]
pub(crate) enum QueryReading {
    /// The records of `run_local`: an input without records is NCBI's `Query is Empty!`
    /// (the inputs of the ABI are read as regular files).
    Records,
    /// The CLI's reader reached the end of its input (`CBlastInput::End()`); without
    /// records, no query batch is read.
    End,
    /// The CLI's reader found only blank and comment lines where it expected a record
    /// (`eEOF`), which ends the batch being read.
    BlankLines,
    /// The CLI's reader failed on the record after the last one read (a reader error, or
    /// LOSAT's explicit rejection of a line that NCBI would fetch as a Seq-id), after
    /// writing `warnings`.
    Error { error: ReadError, warnings: Vec<u8> },
}

impl QueryReading {
    /// Whether NCBI's reader reads on after the last record (`End()` is false): an `eEOF`
    /// or an error then comes in the batch being read.
    pub(crate) fn reads_past_records(&self) -> bool {
        matches!(self, QueryReading::BlankLines | QueryReading::Error { .. })
    }
}

/// NCBI's message for a query or subject without letters.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup.hpp:190-198
/// ```c
///     TSeqPos size() const {
///         TSeqPos retval = x_Size();
///         if (retval == 0) {
///             NCBI_THROW(CBlastException, eInvalidArgument,
///                        "Sequence contains no data");
///         }
///         return retval;
///     }
/// ```
pub(crate) const NO_DATA_MESSAGE: &str = "Sequence contains no data";

#[cfg(test)]
mod tests {
    use super::{get_query_batch_size, ncbi_string_to_int, BatchSizeMixer};
    use crate::blastinput::seq_range::next_ranged_batch_end;

    /// The end of the batch that starts at `start` when no record is skipped. The batch size
    /// is a `TSeqPos` (32-bit unsigned): an `Int4` batch size converts to it.
    fn next_query_batch_end(lengths: &[usize], start: usize, batch_size: u32) -> usize {
        next_ranged_batch_end(lengths, &vec![false; lengths.len()], start, batch_size)
    }

    // blast_input.cpp:134-171: a batch takes records until their letters reach the batch
    // size; each searched query belongs to the batch that read it.
    #[test]
    fn every_searched_query_knows_where_its_batch_starts() {
        use super::QueryReading;
        use crate::blastinput::fasta_reader::FastaRecord;
        use crate::blastinput::seq_range::QueryInput;
        let records: Vec<FastaRecord> = [4usize, 7, 3, 12, 1]
            .iter()
            .enumerate()
            .map(|(i, &len)| FastaRecord::new(format!("q{i}"), b"", &vec![b'A'; len]))
            .collect();
        let input = QueryInput::whole(&records);
        let batches =
            super::read_query_batches(&records, &records, &input, QueryReading::Records, 10, None);
        assert!(batches.stop.is_none());
        assert_eq!(batches.starts, vec![(0, 0), (2, 2), (4, 4)]);
        let starts: Vec<(usize, usize)> = (0..5).map(|q| batches.start_of(q)).collect();
        assert_eq!(starts, vec![(0, 0), (0, 0), (2, 2), (2, 2), (4, 4)]);
    }

    /// The records of `fasta` read with `-parse_deflines` (BLASTP's reader).
    fn parsed_records(fasta: &[u8]) -> Vec<crate::blastinput::fasta_reader::FastaRecord> {
        use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};
        let config = ReaderConfig::query("BLASTP", true, false).with_parse_deflines(true);
        let mut source = FastaInputSource::from_bytes(fasta, config);
        read_all(&mut source, &mut |_| Ok(())).unwrap()
    }

    /// The scope keys of the local id `id`.
    fn local_keys(id: &str) -> Vec<u64> {
        use crate::seqid::model::{ObjectId, ParsedSeqId};
        let mut keys = Vec::new();
        crate::seqid::defline_set::scope_keys(
            &ParsedSeqId::Local(ObjectId::Str(id.as_bytes().to_vec())),
            &mut keys,
        );
        keys
    }

    // blast_format.cpp:2295-2313 (series Q of line RG, Q1c): a batch's scope holds its own queries;
    // a call that searches from a later query numbers them from there.
    #[test]
    fn query_scopes_number_the_searched_queries_of_each_batch() {
        use super::{QueryScopeSpan, QueryScopes};
        let records = parsed_records(b">a x\nAC\n>b x\nAC\n>a y\nAC\n>c x\nAC\n");
        let mut scopes = QueryScopes::new(QueryScopeSpan::Batch, "BLASTP");
        scopes.add_batch(&records[..2], 0).unwrap();
        scopes.add_batch(&records[2..], 2).unwrap();
        assert!(scopes.of_query(0).has_key(&local_keys("b")));
        assert!(scopes.of_query(1).has_key(&local_keys("b")));
        assert!(!scopes.of_query(2).has_key(&local_keys("b")));
        assert!(scopes.of_query(3).has_key(&local_keys("c")));
        assert!(!scopes.of_query(1).has_key(&local_keys("c")));
        let from = scopes.numbered_from(2);
        assert!(from.of_query(0).has_key(&local_keys("c")));
        assert!(!from.of_query(1).has_key(&local_keys("b")));
        // `-subject` and outfmt 5: one scope, in which the second `a` is refused.
        let mut run = QueryScopes::new(QueryScopeSpan::Run, "BLASTP");
        run.add_batch(&records[..2], 0).unwrap();
        let error = run.add_batch(&records[2..], 2).unwrap_err().to_string();
        assert!(
            error.contains("two query records have the Seq-id lcl|a"),
            "{error}"
        );
        // The ids of every query before the batches add theirs; a read ahead adds nothing again.
        let whole = QueryScopes::new(QueryScopeSpan::Run, "BLASTP").with_whole_run_of(&records);
        assert!(whole.of_query(0).has_key(&local_keys("c")));
        let mut ahead = QueryScopes::of_records(QueryScopeSpan::Run, "BLASTP", &records[1..])
            .unwrap()
            .read_ahead();
        ahead.add_batch(&records[1..], 0).unwrap();
        // A batch's scope has no ids of the whole run.
        let mut batch =
            QueryScopes::new(QueryScopeSpan::Batch, "BLASTP").with_whole_run_of(&records);
        batch.add_batch(&records[..2], 0).unwrap();
        assert!(!batch.of_query(0).has_key(&local_keys("c")));
    }

    // The refusal of two queries of one scope is a refusal of the run, not a stop of NCBI's batch
    // loop: the caller returns it before any report (`QueryBatches::query_refusal`).
    #[test]
    fn two_queries_of_one_scope_are_the_runs_refusal() {
        use super::{QueryReading, QueryScopeSpan, QueryScopes};
        use crate::blastinput::seq_range::QueryInput;
        let records = parsed_records(b">a x\nAC\n>a y\nAC\n");
        let input = QueryInput::whole(&records);
        for (span, refused) in [(QueryScopeSpan::Batch, false), (QueryScopeSpan::Run, true)] {
            let batches = super::read_query_batches(
                &records,
                &records,
                &input,
                QueryReading::Records,
                1,
                Some(QueryScopes::new(span, "BLASTP")),
            );
            assert!(batches.stop.is_none());
            assert_eq!(batches.query_refusal.is_some(), refused);
            assert_eq!(batches.starts.len(), if refused { 1 } else { 2 });
        }
        let batches = super::read_query_batches(
            &records,
            &records,
            &input,
            QueryReading::Records,
            10,
            Some(QueryScopes::new(QueryScopeSpan::Batch, "BLASTP")),
        );
        assert!(batches.query_refusal.is_some() && batches.starts.is_empty());
    }

    // psiblast_app.cpp:468 and psiblast_aux_priv.cpp:636-649: a batch of one letter takes one
    // query with letters and the queries without letters before it; `CPsiBlast` takes one query
    // only (NCBI 2.17.0: `tk, empty, pol` stops with this error after tk's report; a last record
    // without residues crashes it).
    #[test]
    fn a_psiblast_batch_holds_one_query() {
        use super::QueryReading;
        use crate::blastinput::fasta_reader::FastaRecord;
        use crate::blastinput::seq_range::QueryInput;
        let records = |lengths: &[usize]| -> Vec<FastaRecord> {
            lengths
                .iter()
                .enumerate()
                .map(|(i, &len)| FastaRecord::new(format!("q{i}"), b"", &vec![b'A'; len]))
                .collect()
        };
        let read = |records: &[FastaRecord]| {
            let input = QueryInput::whole(records);
            super::read_psiblast_query_batches(
                records,
                records,
                &input,
                QueryReading::Records,
                None,
            )
        };
        let whole = records(&[4, 7, 3]);
        let batches = read(&whole);
        assert!(batches.stop.is_none());
        assert_eq!(batches.starts, vec![(0, 0), (1, 1), (2, 2)]);
        let joined = records(&[4, 0, 3]);
        let batches = read(&joined);
        assert_eq!(batches.starts, vec![(0, 0)]);
        assert_eq!((batches.input_done, batches.searched_done), (1, 1));
        let error = batches.stop.expect("two queries in one batch").1;
        let native = error
            .downcast_ref::<crate::cli::NativeError>()
            .expect("an NCBI error");
        assert_eq!(native.exit, 3);
        assert_eq!(
            native.message,
            "BLAST engine error: PSI-BLAST only accepts one protein sequence as query\n"
        );
        for lengths in [&[4, 0][..], &[0][..]] {
            let alone = records(lengths);
            let batches = read(&alone);
            let error = batches.stop.expect("a query without letters alone").1;
            assert!(error.downcast_ref::<crate::cli::NativeError>().is_none());
            assert!(error
                .to_string()
                .contains("a query without letters alone in its batch"));
            assert_eq!(batches.starts.len(), lengths.len() - 1);
        }
        // blastp's batches of the same records: the record without letters is searched with the
        // next one, without an error.
        let input = QueryInput::whole(&joined);
        let blastp =
            super::read_query_batches(&joined, &joined, &input, QueryReading::Records, 1, None);
        assert!(blastp.stop.is_none());
        assert_eq!(blastp.starts, vec![(0, 0), (1, 1)]);
    }

    #[test]
    fn the_query_that_reaches_the_batch_size_stays_in_the_batch() {
        assert_eq!(next_query_batch_end(&[19_999, 2, 1], 0, 20_000), 2);
        assert_eq!(next_query_batch_end(&[19_999, 2, 1], 2, 20_000), 3);
    }

    #[test]
    fn batch_size_strings_are_read_as_ncbis_string_to_int() {
        use std::ffi::OsStr;
        for (text, value) in [
            ("100", Some(100)),
            ("+100", Some(100)),
            ("-1", Some(-1)),
            ("007", Some(7)),
            ("2147483647", Some(i32::MAX)),
            ("-2147483648", Some(i32::MIN)),
            ("2147483648", None),
            ("", None),
            (" 100", None),
            ("100 ", None),
            ("1e5", None),
            ("abc", None),
            ("+", None),
            ("1,000", None),
            ("0x10", None),
            // NCBI reads up to a NUL (`while (char ch = str[pos])`, ncbistr.cpp:841).
            ("12\u{0}x", Some(12)),
            ("\u{0}12", None),
        ] {
            assert_eq!(ncbi_string_to_int(OsStr::new(text)), value, "{text:?}");
        }
        assert_eq!(get_query_batch_size(None), Ok(0));
        assert_eq!(get_query_batch_size(Some(OsStr::new("5000"))), Ok(5000));
        assert_eq!(
            get_query_batch_size(Some(OsStr::new(""))),
            Err(String::new())
        );
    }

    #[test]
    fn a_negative_batch_size_is_a_huge_tseqpos() {
        assert_eq!(
            next_query_batch_end(&[3000, 1999, 2, 7], 0, -1i32 as u32),
            4
        );
        assert_eq!(next_query_batch_end(&[3000, 1999, 2, 7], 0, 1), 1);
    }

    #[test]
    fn the_next_batch_takes_queries_until_the_batch_size_is_reached() {
        assert_eq!(next_query_batch_end(&[3000, 1999, 2, 7], 0, 5000), 3);
        assert_eq!(next_query_batch_end(&[3000, 1999, 2, 7], 3, 5000), 4);
        assert_eq!(next_query_batch_end(&[3000, 2000, 2], 0, 5000), 2);
    }

    // The values follow NCBI's CBatchSizeMixer (blast_app_util.cpp:65-81).
    #[test]
    fn the_mixer_is_ncbis() {
        let mut mixer = BatchSizeMixer::new(1_000_000 - 1000);
        assert_eq!(mixer.batch_size(None), 5000);
        // 4999 hits: ratio 1.0, 1e6 residues, above the maximum.
        assert_eq!(mixer.batch_size(Some(4999)), 999_000);
        // 999 hits in 999000 residues: ratio 0.001 (reset after the clamp), 1e9, clamped.
        assert_eq!(mixer.batch_size(Some(998)), 999_000);
        // 9999 hits in 5000 residues: ratio 2.0, 500000 residues.
        let mut mixer = BatchSizeMixer::new(4_999_000);
        mixer.batch_size(None);
        assert_eq!(mixer.batch_size(Some(9999)), 500_000);
        // Mixed: 0.3 * (20000 / 500000) + 0.7 * 2.0 = 1.412 -> 708215.
        assert_eq!(mixer.batch_size(Some(19_999)), 708_215);
        // No hits: 1e6 / (1 / 5000) overflows Int4, INT_MIN on x86-64, so the smallest batch.
        let mut mixer = BatchSizeMixer::new(999_000);
        mixer.batch_size(None);
        assert_eq!(mixer.batch_size(Some(0)), 100);
        // A large subject set raises the target and the first batch, within the maximum.
        let mut mixer = BatchSizeMixer::new(4_999_000);
        mixer.set_target_hits(3_000_000_000_i64.min(i64::from(i32::MAX)) as i32);
        assert_eq!(mixer.batch_size(None), 4_999_000);
        let mut mixer = BatchSizeMixer::new(999_000);
        mixer.set_target_hits(60_000_000);
        assert_eq!(mixer.batch_size(None), 300_000);
    }
}
