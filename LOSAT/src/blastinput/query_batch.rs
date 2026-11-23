//! Query batches of a search (NCBI `CBlastInput::GetNextSeqBatch`).

use crate::blastinput::fasta_reader::believed::BelievedIds;
use crate::blastinput::fasta_reader::{FastaRecord, ReadError};
use crate::wording::ENGINE_ERROR;

/// The query batch size from the environment variable `BATCH_SIZE` (`value`; 0 when it is not
/// set: blastn then uses `CBatchSizeMixer`, the other programs their default). `Err` is the
/// `CStringException` of `NStr::StringToInt` (flags 0: an optional sign and decimal digits within
/// the range of an `int`, no spaces) for a value it cannot convert, an empty one too.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:85-91
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:99
pub(crate) fn get_query_batch_size(
    value: Option<&std::ffi::OsStr>,
) -> Result<i32, crate::ncbi_exception::NcbiException> {
    match value {
        None => Ok(0),
        Some(value) => crate::ncbi_str::string_to_int(value.as_encoded_bytes(), 0),
    }
}

/// NCBI's adaptive query batch size of blastn: the first batch aims at 1/200 of the target
/// hits, and each later batch at the target over the (mixed) ratio of successful ungapped
/// extensions per residue of the batches before.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:57-83
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
    /// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:65-81
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

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:261-273
// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:298-300
/// The query batch sizes of a run's batch loop: a fixed size, or those of BLASTN's application
/// without `BATCH_SIZE` (`mixed`), which its `CBatchSizeMixer` gives: the first from the
/// database's total length, each later one from the successful initial extensions of the batch
/// before (`searched`). Which one a run has is its running application's (`Program::
/// fixed_query_batch_size`), whatever program's options it runs.
#[derive(Clone, Debug)]
pub(crate) enum LoopBatchSizes {
    Fixed(u32),
    Mixed { size: u32, mixer: BatchSizeMixer },
}

impl LoopBatchSizes {
    /// The sizes of BLASTN's application: `fixed_batch_size` is `BATCH_SIZE` (0 when it is not
    /// set; the application takes a size of 0 as no fixed size), `max_batch_size` the mixer's
    /// maximum (the chunk size of the options' program less 1000) and `total_length` the database's
    /// (`GetDbTotalLength`, which sets the target hits only when it is positive).
    pub(crate) fn of_blastn_app(
        fixed_batch_size: i32,
        max_batch_size: i32,
        total_length: i64,
    ) -> Self {
        if fixed_batch_size != 0 {
            return Self::Fixed(fixed_batch_size as u32);
        }
        let mut mixer = BatchSizeMixer::new(max_batch_size);
        if total_length > 0 {
            mixer.set_target_hits((total_length / 3000) as i32);
        }
        // `SetBatchSize` of the input, given `GetBatchSize` of the mixer, converts the `Int4` to
        // `TSeqPos`.
        let size = mixer.batch_size(None) as u32;
        Self::Mixed { size, mixer }
    }

    /// The size of the next batch (`TSeqPos`).
    pub(crate) fn current(&self) -> u32 {
        match self {
            Self::Fixed(size) | Self::Mixed { size, .. } => *size,
        }
    }

    /// Whether the sizes follow the extensions of each batch.
    pub(crate) fn mixes(&self) -> bool {
        matches!(self, Self::Mixed { .. })
    }

    /// The batch just searched made `extensions` successful initial extensions (NCBI's `Int4`
    /// `good_init_extends` of the batch's diagnostics, `GetNumExtensions`).
    pub(crate) fn searched(&mut self, extensions: i32) {
        if let Self::Mixed { size, mixer } = self {
            *size = mixer.batch_size(Some(extensions)) as u32;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2295-2313
// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-259 (each program's batch loop)
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:157-163 (cited at
// `seq_range::next_ranged_batch_end`)
/// NCBI's reader's count of the query records into its batches under `-parse_deflines`
/// (`CBlastInput::GetNextSeqBatch`): each record by the length of its id in the scope, where a
/// repeated id (`Meets::Repeat`) finds its first record, so that the repeat counts the first
/// record's letters. The scope is NCBI's (`QueryScopeSpan`): the formatter empties it after each
/// batch (`ResetScopeHistory`) but with `-subject` and outfmt 5. Every cut of the query batches
/// counts so (`next_batch`; BLASTX's readers by `count` and `end_batch`): the batch of the first
/// repeat is NCBI's, with the records that its set-up reads and its reader's error, where its
/// reader meets one (line WQA D10 round 6, 2026-11-20; review 4 of D10, M1). The batches before
/// it hold no repeat, and none follows it.
///
/// NCBI BLAST+ 2.17.0 (`fix-reader-error/r6/probe` of line WQA, 2026-11-20): with `BATCH_SIZE=60`,
/// a `q1` of 60 letters, a `q1` of 55 and a line that the reader refuses, blastp with `-subject`
/// counts the second `q1` by 60, searches it in a batch of its own and then stops with the
/// reader's error; with `-db` it counts it by its own 55 letters and its batch meets the error. A
/// third `q1` after the scope was emptied counts the letters of the second (`chainB50`), not the
/// first's (`chainB5`). A record whose ids meet the scope otherwise (`Meets::Other`) counts its
/// own letters (`>gi|5|lcl|q1` after `>q1`, review 4's probe `oth1`). A record that a `-query_loc`
/// range skips adds nothing but joins the scope, as the reader added it before the range's check:
/// a later `q1` of 60 letters after a skipped `q1` of 10 counts by 10 (`r6/probe/skip2`; `skip2c`
/// with another id counts 60).
#[derive(Debug, Clone)]
pub(crate) struct BatchCount {
    span: QueryScopeSpan,
    program: &'static str,
    scope: QueryScope,
}

impl BatchCount {
    /// The count of a search whose scopes are `scopes` (`None` without `-parse_deflines`, whose
    /// records count their own letters).
    pub(crate) fn of(scopes: Option<&QueryScopes>) -> Option<Self> {
        scopes.map(|scopes| Self::new(scopes.span, scopes.program))
    }

    /// The count of a search as `QueryScopes::of_search` makes its scopes: `None` without
    /// `-parse_deflines` (`believed`).
    pub(crate) fn of_search(
        believed: bool,
        subject: bool,
        xml: bool,
        program: &'static str,
    ) -> Option<Self> {
        believed.then(|| Self::new(QueryScopeSpan::of(subject, xml), program))
    }

    /// The count of a search whose scopes span `span` (`QueryScopes::of_search`).
    pub(crate) fn new(span: QueryScopeSpan, program: &'static str) -> Self {
        Self {
            span,
            program,
            scope: QueryScope::default(),
        }
    }

    /// The letters by which NCBI's reader counts a record with the ids `believed` (`None`: an id
    /// of the reader's own) and `length` letters into its batch, 0 where a `-query_loc` range
    /// skips it (`skipped`); the record joins the scope either way.
    pub(crate) fn count(
        &mut self,
        believed: Option<&crate::blastinput::fasta_reader::believed::BelievedIds>,
        length: usize,
        skipped: bool,
    ) -> usize {
        let own = if skipped { 0 } else { length };
        let Some(believed) = believed else {
            return own;
        };
        let counted = match self.scope.meets(&believed.ids, self.program) {
            Meets::Repeat { first: Some(first) } if !skipped => first,
            _ => own,
        };
        let record = (self.scope.records, length);
        self.scope.records += 1;
        self.scope.add_after_refusal(&believed.ids, record);
        counted
    }

    /// The end of a batch: the formatter empties the scope after the batch's reports, but with
    /// `-subject` and outfmt 5 (`QueryScopeSpan::Batch`).
    pub(crate) fn end_batch(&mut self) {
        if self.span == QueryScopeSpan::Batch {
            self.scope = QueryScope::default();
        }
    }

    /// The query batch of the input that starts at input record `start`
    /// (`seq_range::next_ranged_batch`), each record counted by `count` (with `believed`, the ids
    /// of an input record), or by its own letters without one; whether it reached its size.
    pub(crate) fn next_batch(
        count: Option<&mut Self>,
        believed: &dyn Fn(
            usize,
        ) -> Option<
            std::sync::Arc<crate::blastinput::fasta_reader::believed::BelievedIds>,
        >,
        input_lengths: &[usize],
        skipped: &[bool],
        start: usize,
        batch_size: u32,
    ) -> (usize, bool) {
        let Some(count) = count else {
            return crate::blastinput::seq_range::next_ranged_batch(
                input_lengths,
                skipped,
                start,
                batch_size,
            );
        };
        let batch = crate::blastinput::seq_range::next_counted_batch(
            input_lengths.len(),
            start,
            batch_size,
            &mut |record| {
                count.count(
                    believed(record).as_deref(),
                    input_lengths[record],
                    skipped[record],
                )
            },
        );
        count.end_batch();
        batch
    }
}

/// The ids of the query records of one scope (`QueryScopes`): the keys of each id
/// (`seqid::defline_set::scope_keys`, and the own kind of `id_scope::version_keys` and
/// `id_scope::name_keys`).
#[derive(Clone, Debug, Default)]
pub(crate) struct QueryScope {
    keys: std::collections::HashSet<u64>,
    /// The record (its number in the scope, `records`) that is first in the scope with each key
    /// of `scope_keys`, and the letters of the whole record (`add`), by which NCBI's reader counts
    /// a later record with that id (`Meets::Repeat`).
    first_records: std::collections::HashMap<u64, (usize, usize)>,
    /// The records of the batches added to the scope before the one being added.
    records: usize,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:138-170
// (`CBlastInput::GetNextSeqBatch`, cited at `seq_range::next_ranged_batch_end`)
/// How the ids of a query record meet the ids of the records of its scope before it
/// (`QueryScope::meets`; line WQA D10 rounds 5 and 5r, 2026-11-20).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Meets {
    /// No id of an earlier record: the record is first in the scope with its ids, and NCBI's
    /// set-up reads its own letters.
    New,
    /// Every Seq-id of the record is one of an earlier record: NCBI's scope finds the earlier
    /// record by the id. Its reader counts the record into its batch by the earlier record's
    /// length (`GetLength` of the id in the scope, `BatchCount`), and its set-up reads the record's
    /// interval in the earlier record's letters. `first`: that length, where the ids find one
    /// earlier record; `None` where they find two, or ids that the scope holds without a record of
    /// their own (`QueryScopes::with_whole_run`). NCBI BLAST+ 2.17.0 (`r5r/probe` of line WQA D10,
    /// 2026-11-20): with `BATCH_SIZE=210`, a `q1` of 300 `X` (or `N`), then `q1` of 200 letters,
    /// `q4` of 200 `X`, `q2` of 400 letters and subjects too short for `-searchsp 1000`, every
    /// program writes `q1`'s Karlin-Altschul warning twice and the error of the subjects: the
    /// second `q1` fills a batch of its own, which is not searched, and `q4`'s warning goes with
    /// the failed batch of `q4` and `q2`.
    Repeat { first: Option<usize> },
    /// Any other meeting: some ids new and some not, a name with and without an accession, an
    /// accession with and without a version (BLASTN). NCBI's scope finds a record by its best
    /// choice of the ids (`FindBestChoice`): its reader counted `>gi|5|lcl|q1` after `>q1` by its
    /// own letters (review 4 of D10, probe `oth1`), as `BatchCount` counts such a record.
    Other,
}

impl QueryScope {
    /// Adds the ids `ids` of a query record (`record`: its number in the scope and its letters), or
    /// LOSAT's refusal of a record whose id another record of the scope has
    /// (`believed::duplicate_id`).
    fn add(
        &mut self,
        ids: &[crate::seqid::ParsedSeqId],
        record: (usize, usize),
        program: &str,
    ) -> anyhow::Result<()> {
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
        self.note_first(ids, record);
        Ok(())
    }

    /// Notes `record` (its number in the scope and its letters) for the keys of `ids` that no
    /// earlier record of the scope has.
    fn note_first(&mut self, ids: &[crate::seqid::ParsedSeqId], record: (usize, usize)) {
        let mut keys = Vec::new();
        for id in ids {
            keys.clear();
            crate::seqid::defline_set::scope_keys(id, &mut keys);
            for &key in &keys {
                self.first_records.entry(key).or_insert(record);
            }
        }
    }

    /// How the ids `ids` of a record meet those of the scope (`Meets`), by the kinds of
    /// `believed::duplicate_id`.
    fn meets(&self, ids: &[crate::seqid::ParsedSeqId], program: &str) -> Meets {
        let mut keys = Vec::new();
        let (mut known, mut new, mut other) = (false, false, false);
        // The earlier record of the scope that the ids find, while they find one only.
        let mut first: Option<Option<(usize, usize)>> = None;
        for id in ids {
            keys.clear();
            crate::seqid::defline_set::scope_keys(id, &mut keys);
            for key in &keys {
                let found = match self.first_records.get(key) {
                    Some(&record) => Some(record),
                    None if self.keys.contains(key) => None,
                    None => {
                        new = true;
                        continue;
                    }
                };
                known = true;
                first = match first {
                    None => Some(found),
                    Some(earlier) if earlier == found => Some(earlier),
                    Some(_) => Some(None),
                };
            }
            other |= crate::seqid::id_scope::name_keys(id)
                .is_some_and(|(_, name)| self.keys.contains(&name));
            other |= program == "BLASTN"
                && crate::seqid::id_scope::version_keys(id)
                    .is_some_and(|(_, version)| self.keys.contains(&version));
        }
        let first = first.flatten().map(|(_, letters)| letters);
        match (known, new, other) {
            (false, _, false) => Meets::New,
            (true, false, false) => Meets::Repeat { first },
            _ => Meets::Other,
        }
    }

    /// Adds the ids `ids` of a record (its number in the scope and its letters) whatever the scope
    /// holds (a record after LOSAT's refusal in its batch, `QueryScopes::add_batch_ids`).
    fn add_after_refusal(&mut self, ids: &[crate::seqid::ParsedSeqId], record: (usize, usize)) {
        self.add_ids(ids);
        self.note_first(ids, record);
    }

    /// Adds the ids `ids` of query records (`BelievedIds`), whatever ids the scope has
    /// (`QueryScopes::with_whole_run`): the keys that `add` keeps of an id (`scope_keys` and the
    /// own kind of `name_keys` and `version_keys`), which `has_key` asks.
    pub(crate) fn add_ids_of<'r>(
        &mut self,
        ids: impl IntoIterator<Item = &'r crate::blastinput::fasta_reader::believed::BelievedIds>,
    ) {
        for believed in ids {
            self.add_ids(&believed.ids);
        }
    }

    /// The scope of query records read without `-parse_deflines`, by their generated local ids
    /// `local_ids` (`Query_N`, without `lcl|`; `believed::generated_ids`): the keys that `add_ids_of`
    /// keeps of such an id (`scope_keys`: a local string id in any case, as NCBI's scope finds
    /// one).
    pub(crate) fn of_generated_ids<'r>(local_ids: impl IntoIterator<Item = &'r str>) -> Self {
        let mut scope = Self::default();
        let mut keys = Vec::new();
        for id in local_ids
            .into_iter()
            .flat_map(crate::blastinput::fasta_reader::believed::generated_ids)
        {
            keys.clear();
            crate::seqid::defline_set::scope_keys(&id, &mut keys);
            scope.keys.extend(keys.iter().copied());
        }
        scope
    }

    /// Whether a query of the scope has an id with one of the keys `keys`.
    pub(crate) fn has_key(&self, keys: &[u64]) -> bool {
        keys.iter().any(|key| self.keys.contains(key))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:364-367
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:395 (the
    // record joins the scope before `x_FastaToSeqLoc` checks its range, 449-452, whose exception
    // `CBlastInput::GetNextSeqBatch` drops the record for, blast_input.cpp:146-155)
    /// The scope of the query records `records` as NCBI's reader added them to its scope (each a
    /// record's local id and the ids that `-parse_deflines` read, `None` without the option): the
    /// generated local ids (`of_generated_ids`, the reader's `Query_N` of the whole input, which
    /// the records hold) where the first record has no believed ids, else the ids of the records
    /// (`add_ids_of`). One owner of the scope of query records: the scopes of the query-anchored
    /// formats (`AnchoredBatchScopes`) and the scope at a batch's prefetch
    /// (`common::database::PrefetchQueries`).
    pub(crate) fn of_shown<'r>(
        records: impl IntoIterator<
            Item = (
                &'r str,
                Option<&'r crate::blastinput::fasta_reader::believed::BelievedIds>,
            ),
        >,
    ) -> Self {
        let mut records = records.into_iter().peekable();
        let generated = records
            .peek()
            .is_none_or(|(_, believed)| believed.is_none());
        if generated {
            Self::of_generated_ids(records.map(|(local_id, _)| local_id))
        } else {
            let mut scope = Self::default();
            scope.add_ids_of(records.filter_map(|(_, believed)| believed));
            scope
        }
    }

    /// `of_shown` of the query records `records` (`FastaRecord::local_id`, `FastaRecord::believed`).
    pub(crate) fn of_shown_records(records: &[FastaRecord]) -> Self {
        Self::of_shown(
            records
                .iter()
                .map(|record| (record.local_id.as_str(), record.believed())),
        )
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:535-537 (c++/src/app/blast/tblastn_app.cpp:
    // 224-226 likewise)
    /// The scope of the query of a PSSM (`-in_pssm`), whose Bioseq NCBI adds to its scope with
    /// its Seq-ids `ids` (the keys that `add_ids_of` keeps of each): its local id too, which the
    /// record of the PSSM's query shows without `lcl|` (`pssm_search::PsiSearch::query_record`).
    pub(crate) fn of_bioseq_ids(ids: &[crate::seqid::ParsedSeqId]) -> Self {
        let mut scope = Self::default();
        scope.add_ids(ids);
        scope
    }

    /// Adds the keys that `add_ids_of` keeps of each of the Seq-ids `ids`.
    fn add_ids(&mut self, ids: &[crate::seqid::ParsedSeqId]) {
        let mut keys = Vec::new();
        for id in ids {
            keys.clear();
            crate::seqid::defline_set::scope_keys(id, &mut keys);
            self.keys.extend(keys.iter().copied());
            self.keys
                .extend(crate::seqid::id_scope::name_keys(id).map(|(own, _)| own));
            self.keys
                .extend(crate::seqid::id_scope::version_keys(id).map(|(own, _)| own));
        }
    }

    /// Adds the keys of the queries of `other`.
    pub(crate) fn extend(&mut self, other: &QueryScope) {
        self.keys.extend(other.keys.iter().copied());
    }
}

/// The scopes of the query batches of one report, for the refusals of the query-anchored formats
/// (outfmt 1 to 4; `common::database::refuse_query_ids_of_anchored_reports`), one owner for
/// BLASTN's and TBLASTN's reports: the scope of each batch (`spans`: the query records, by index,
/// that NCBI's formatter adds to its scope for the batch's reports) holds the generated local
/// ids (`QueryScope::of_generated_ids`) of queries read without `-parse_deflines` (`generated`:
/// the first query has no `believed` ids), or the ids of the queries that the reports show
/// (`QueryScope::add_ids_of`: `-parse_deflines`, a PSSM query shown with its ids).
pub(crate) struct AnchoredBatchScopes {
    generated: bool,
    spans: Vec<std::ops::Range<usize>>,
    scopes: Vec<QueryScope>,
    none: QueryScope,
}

impl AnchoredBatchScopes {
    /// The scopes of the batches `spans` of the query records `records` (`of_ids` of their ids).
    pub(crate) fn new(records: &[FastaRecord], spans: Vec<std::ops::Range<usize>>) -> Self {
        let ids: Vec<(&str, Option<&BelievedIds>)> = records
            .iter()
            .map(|record| (record.local_id.as_str(), record.believed()))
            .collect();
        Self::of_ids(&ids, spans)
    }

    /// `new` of the query records whose local ids and `-parse_deflines` Seq-ids are `ids` (BLASTX's
    /// records, which have their own type): a span past the records is cut to them;
    /// `QueryScope::of_shown` of each (every record is read with `-parse_deflines` or none is).
    pub(crate) fn of_ids(
        ids: &[(&str, Option<&BelievedIds>)],
        spans: Vec<std::ops::Range<usize>>,
    ) -> Self {
        let generated = ids.first().is_none_or(|(_, believed)| believed.is_none());
        let scopes = spans
            .iter()
            .map(|span| {
                QueryScope::of_shown(
                    ids[span.start.min(ids.len())..span.end.min(ids.len())]
                        .iter()
                        .copied(),
                )
            })
            .collect();
        Self {
            generated,
            spans,
            scopes,
            none: QueryScope::default(),
        }
    }

    /// The scopes of the batches that start at the query records `starts` (the first batch at 0
    /// whether or not `starts` lists it; each batch ends where the next starts, the last at the
    /// end of `records`).
    pub(crate) fn of_starts(records: &[FastaRecord], starts: &[usize]) -> Self {
        let mut starts: Vec<usize> = starts.to_vec();
        if starts.first() != Some(&0) {
            starts.insert(0, 0);
        }
        let spans = starts
            .iter()
            .enumerate()
            .map(|(index, &start)| start..starts.get(index + 1).copied().unwrap_or(records.len()))
            .collect();
        Self::new(records, spans)
    }

    /// Whether the queries were read without `-parse_deflines` (their generated ids).
    pub(crate) fn generated(&self) -> bool {
        self.generated
    }

    /// The scope of the batch that holds the query record `query` (empty past the batches).
    pub(crate) fn of_query(&self, query: usize) -> &QueryScope {
        self.spans
            .iter()
            .position(|span| span.contains(&query))
            .map_or(&self.none, |batch| &self.scopes[batch])
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
    /// The batches that the search has added, and LOSAT's refusal that a reading before the
    /// search met (`read_ahead`, `with_refusal`) with the number of batches added before the
    /// refusal's: the search's add of that batch gives it (`add_batch_ids`).
    added: usize,
    refusal_ahead: Option<(usize, anyhow::Error)>,
}

/// LOSAT's refusal of a query record with the id of another query record of its scope
/// (`QueryScope::add`), with the records of the refused record's batch (NCBI's batch,
/// `BatchCount`) that NCBI's set-up of that batch reads with their own letters. The run stops at
/// that batch in one form in every program, the command line's and the library's (line WQA D10
/// rounds 5 and 6, 2026-11-20; `refused_batch_error`).
#[derive(Clone, Debug)]
pub(crate) struct ScopeRefusal {
    message: String,
    /// By record of the batch, as `QueryScopes::add_batch_ids` took them: whether NCBI's set-up
    /// of the batch reads it with its own letters (a record first in the scope with its ids),
    /// where LOSAT can follow NCBI's set-up of every record of the batch; `None` otherwise.
    first_records: Option<Vec<bool>>,
}

impl std::fmt::Display for ScopeRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for ScopeRefusal {}

impl ScopeRefusal {
    /// The refusal that the error `error` is, if it is one.
    pub(crate) fn of(error: &anyhow::Error) -> Option<&ScopeRefusal> {
        error.downcast_ref::<ScopeRefusal>()
    }

    /// The records of the batch that NCBI's set-up reads with their own letters
    /// (`first_records`), where there is one; `None` otherwise.
    pub(crate) fn setup_records(&self) -> Option<&[bool]> {
        self.first_records
            .as_deref()
            .filter(|first| first.iter().any(|&first| first))
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:259-276 (each program's batch loop)
// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980 (`BLAST_GapAlignSetUp`)
/// The error of a query batch that LOSAT refuses for two queries of one scope with one id
/// (`refusal`), in NCBI's order of the batch (line WQA D10 round 5, 2026-11-20; review 3 of D10,
/// M1). The caller has written the reports of the batches before and the batch's reader messages,
/// and acted on the batch's checks before its set-up that read no query letters (NCBI's: the
/// environment of the query splitter, a batch of queries without letters, BLASTN's WindowMasker
/// file, TBLASTX's matrix name).
///
/// - With `-export_search_strategy` (`export`) the refusal comes before the batch's request: NCBI's
///   request holds the letters of the earlier record for the repeated id (`SaveSearchStrategy`),
///   so no request of the batch is written (with one batch, the file stays empty).
/// - Else NCBI sets the batch up (`BLAST_MainSetUp`, which reads the query letters, then
///   `BLAST_GapAlignSetUp`) and fails it where the set-up fails for its subjects: the average
///   subject length against `-searchsp` or `OLD_FSC`, the hit-saving parameters, the first read of
///   a database. Those errors are the run's, the same for every batch with a valid context, and a
///   batch without one is not searched (`CheckInternalData`), its warnings written: where the
///   records that NCBI's set-up of the batch reads with their own letters
///   (`ScopeRefusal::setup_records`: those first in the scope with their ids, where LOSAT can
///   follow NCBI's set-up of every record of the batch) have a valid context, NCBI's set-up of the
///   batch has one too and fails with that error, whatever the repeated records' letters, and the
///   failed batch's warnings are not written. `setup` is the program's set-up of those records (by record of the batch): NCBI's
///   error (`cli::NativeError`) where it fails, `None` where it passes or has no valid context.
/// - Else the refusal, before the batch's search.
///
/// NCBI BLAST+ 2.17.0 (`r5/probe` of line WQA D10): with a query of only X letters, a second query
/// of only X, and the first's id again with letters of its own, NCBI reads the first's letters
/// for the third, sets up no valid context and exits 0 with the Karlin-Altschul warnings, so a
/// batch whose records with their own letters have no valid context is refused (a repeat's letters
/// are those of its earlier record, of the batch or of a batch before, whose letters LOSAT does
/// not keep).
pub(crate) fn refused_batch_error(
    refusal: anyhow::Error,
    export: bool,
    setup: impl FnOnce(&[bool]) -> Option<anyhow::Error>,
) -> anyhow::Error {
    if export {
        return refusal;
    }
    let Some(first) = ScopeRefusal::of(&refusal).and_then(ScopeRefusal::setup_records) else {
        return refusal;
    };
    match setup(first) {
        Some(error) if error.is::<crate::cli::NativeError>() => error,
        _ => refusal,
    }
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
            added: 0,
            refusal_ahead: None,
        }
    }

    /// The scopes that a read of the queries before the search made with the search's batches
    /// (BLASTX's, whose reader gives the batches one at a time): the batches of the search add
    /// nothing again (`add_batch_ids`). Where the read met LOSAT's refusal of two queries with one
    /// id (`refusal`: the batches it added before the refusal's batch, and the refusal), the
    /// search's add of that batch gives the refusal: the run stops at that batch after the reports
    /// of the batches before, as where the batches add theirs, so that an error of NCBI's in the
    /// search of a batch before comes first (line WQA D10 round 4, 2026-11-20: blastx `-outfmt 5`
    /// with subjects too short for `-searchsp`, NCBI BLAST+ 2.17.0).
    pub(crate) fn read_ahead(mut self, refusal: Option<(usize, anyhow::Error)>) -> Self {
        self.read_ahead = true;
        self.added = 0;
        self.refusal_ahead = refusal;
        self
    }

    /// The scopes of a call of the library that searches one batch of a run (`api::run_batches`):
    /// `refusal` is LOSAT's refusal of the batch that the library's scope of the run's batches
    /// met (`ScopeRefusal`, with the batches before the call's), which the call's add of its batch
    /// gives, as the command line's add of that batch gives it (line WQA D10 round 5, 2026-11-20).
    pub(crate) fn with_refusal(mut self, refusal: Option<ScopeRefusal>) -> Self {
        self.refusal_ahead = refusal.map(|refusal| (0, refusal.into()));
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
    /// record with the id of another of the scope (`ScopeRefusal`). The batch is NCBI's
    /// (`BatchCount`).
    pub(crate) fn add_batch(
        &mut self,
        records: &[FastaRecord],
        first_searched: usize,
    ) -> anyhow::Result<()> {
        self.add_batch_ids(
            records
                .iter()
                .map(|record| (record.believed(), record.seq().len())),
            first_searched,
        )
    }

    /// `add_batch` of the records of a batch by their ids (`BelievedIds`, `None` without ids of
    /// their own) and the letters of each whole record (BLASTX's reader has records of its own).
    /// A refused batch's other records are added too, for what NCBI's set-up of the batch reads
    /// (`ScopeRefusal`).
    ///
    /// The batch is NCBI's: its reader counts a repeated record by its earlier record's letters
    /// (`BatchCount`, which every cut of the batches takes). NCBI's set-up of the batch reads each
    /// record that is first in the scope with its ids (`Meets::New`) with its own letters, and a
    /// repeat (`Meets::Repeat`) as its interval in the earlier record's letters. The refusal holds
    /// the first kind (`ScopeRefusal::setup_records`), where every repeat of the batch is one that
    /// LOSAT can follow there: no record meets the scope otherwise (`Meets::Other`, whose record
    /// NCBI finds by its best choice of the ids), every repeat finds one earlier record, and in
    /// BLASTN no repeat is longer than its earlier record (NCBI's set-up of blastn reads the
    /// repeat's interval past the earlier record's letters and stops there: `BLAST engine error:
    /// ... CSeqVector_CI: cannot locate segment`, `r5/probe`; the other programs set such a batch
    /// up as the others, `db.pd.q1c.bp.avglong.subj.6`). Else the batch is refused before its
    /// set-up (line WQA D10 round 6, 2026-11-20; review 4 of D10, M1).
    pub(crate) fn add_batch_ids<'r>(
        &mut self,
        records: impl IntoIterator<
            Item = (
                Option<&'r crate::blastinput::fasta_reader::believed::BelievedIds>,
                usize,
            ),
        >,
        first_searched: usize,
    ) -> anyhow::Result<()> {
        let batch = self.added;
        self.added += 1;
        match self.refusal_ahead.take() {
            Some((refused, refusal)) if refused == batch => return Err(refusal),
            kept => self.refusal_ahead = kept,
        }
        if self.read_ahead {
            return Ok(());
        }
        if self.scopes.is_empty() || self.span == QueryScopeSpan::Batch {
            // A scope of one batch takes no more records: its first records are not needed.
            if let Some(scope) = self.scopes.last_mut() {
                scope.first_records = Default::default();
            }
            self.scopes.push(QueryScope::default());
            self.starts.push(first_searched);
        }
        let scope = self.scopes.last_mut().expect("pushed above");
        let before = scope.records;
        let mut message: Option<String> = None;
        // By record of the batch, whether NCBI's set-up reads it with its own letters, while
        // every record is one that LOSAT can follow there.
        let mut own_letters = Some(Vec::new());
        let mut records_read = 0;
        for (place, (believed, length)) in records.into_iter().enumerate() {
            records_read = place + 1;
            let Some(believed) = believed else {
                if let Some(own) = own_letters.as_mut() {
                    own.push(true);
                }
                continue;
            };
            let follows = match scope.meets(&believed.ids, self.program) {
                Meets::New => Some(true),
                Meets::Repeat { first: Some(first) }
                    if self.program != "BLASTN" || length <= first =>
                {
                    Some(false)
                }
                _ => None,
            };
            own_letters = own_letters.zip(follows).map(|(mut own, follows)| {
                own.push(follows);
                own
            });
            let record = (before + place, length);
            if message.is_some() {
                scope.add_after_refusal(&believed.ids, record);
            } else if let Err(error) = scope.add(&believed.ids, record, self.program) {
                // `duplicate_id` stopped at the id it refuses: the record's other keys.
                scope.add_after_refusal(&believed.ids, record);
                message = Some(error.to_string());
            }
        }
        scope.records += records_read;
        let Some(message) = message else {
            return Ok(());
        };
        Err(ScopeRefusal {
            message,
            first_records: own_letters,
        }
        .into())
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

impl<'a> QueryScopesFrom<'a> {
    /// The scope of the call's query `query`.
    pub(crate) fn of_query(&self, query: usize) -> &'a QueryScope {
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
    /// Whether `stop` is LOSAT's refusal of two queries of one scope with one id, in the first
    /// batch that has them before any batch that NCBI fails at its reading (`QueryScopes::
    /// add_batch`). It stops the run at that batch as NCBI's errors do, after the reports of the
    /// batches before, so that an error of NCBI's in their search comes first (line WQA D10
    /// round 4, 2026-11-20: psiblast's query without Karlin-Altschul parameters, blastp's
    /// subjects too short for `-searchsp`); the caller's own stops of an earlier batch replace it.
    pub(crate) scope_refusal: bool,
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

    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:259-262
    /// Where the searched queries of the run's first batch with a search strategy request end:
    /// the first searched batch, else the batch whose set-up fails (`stop_batch`); `None` when
    /// the run stops before it makes a query factory.
    pub(crate) fn first_request_end(&self) -> Option<usize> {
        match self.starts.first() {
            Some(_) => Some(self.searched_span(0).end),
            None => self.stop_batch.as_ref().map(|span| span.end),
        }
    }

    /// Where the searched queries end of the batch whose formatting NCBI fails for a matrix name
    /// its BLAST XML report data does not know (`blast_xml::matrix::XmlMatrixFailure`): the
    /// search strategy's last request (`StrategyExport::stop_after`). The batch of the run's first
    /// searched query (as `first_request_end`), or of its second when the report is held
    /// (`XmlMatrixFailure::SecondQueryOrEnd`: NCBI's next query of another Seq-id); `None`
    /// when a held report fails at the run's end or the run stops before: every request that
    /// the run makes is written (probe K/k4f: `-subject` 13-16 with one query and an empty second
    /// batch wrote both requests and stopped with that batch's error).
    pub(crate) fn xml_matrix_request_end(
        &self,
        failure: crate::report::blast_xml::matrix::XmlMatrixFailure,
    ) -> Option<usize> {
        let query = failure.query();
        if query == 0 {
            return self.first_request_end();
        }
        (query < self.searched_done).then(|| self.searched_span(query).end)
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:134-171
// NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:378-380
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1013-1025
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
/// - the batch's query splitter, made before its setup, cannot convert `CHUNK_SIZE` or (where
///   it splits, `reads_overlap`: not RPS-BLAST) `OVERLAP_CHUNK_SIZE`
///   (`app::check_query_split_environment`): the first batch stops there;
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
/// with one id stops the run at that batch (`QueryBatches::scope_refusal`), after the batch's
/// errors of NCBI's own reading (a batch that NCBI does not search reports nothing of its scope).
pub(crate) fn read_query_batches(
    input_records: &[FastaRecord],
    searched: &[FastaRecord],
    input: &crate::blastinput::seq_range::QueryInput,
    reading: QueryReading,
    batch_size: u32,
    reads_overlap: bool,
    query_scopes: Option<QueryScopes>,
) -> QueryBatches {
    read_batches(
        input_records,
        searched,
        input,
        reading,
        batch_size,
        &|batch| {
            crate::blastinput::app::check_query_split_environment(false, reads_overlap)
                .err()
                .or_else(|| no_data_batch_error(batch))
        },
        query_scopes,
    )
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:468
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_impl.cpp:88-92
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:636-649
/// psiblast's query batches (`read_query_batches` with a batch size of 1): a batch holds one
/// query with letters and the queries without letters before it. Each batch is a `CPsiBlast` of
/// its own, which takes exactly one query: a batch of several (a query without letters and the
/// next) stops the run with `BLAST engine error: PSI-BLAST only accepts one protein sequence as
/// query` (exit 3). A batch of one query without letters (the last records, or a `-query_loc`
/// interval that starts just past its record's end) crashes NCBI 2.17.0's psiblast (SIGSEGV);
/// LOSAT rejects it there. A query whose Karlin-Altschul parameters cannot be computed stops the
/// run at its batch too (`blast_engine.rs` `search`). The query splitter of the batch's search
/// (`CPsiBlast::Run`, after the constructor's check) stops the first batch where it cannot convert
/// `CHUNK_SIZE` or `OVERLAP_CHUNK_SIZE` (`app::check_query_split_environment`).
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
        &|batch| {
            psiblast_batch_error(batch).or_else(|| {
                crate::blastinput::app::check_query_split_environment(false, true).err()
            })
        },
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
    // NCBI's reader's count of the records into the batches (`BatchCount`).
    let mut count = BatchCount::of(query_scopes.as_ref());
    // The reader is at the end of its input after the last record, unless an `eEOF` or an
    // error ends the reading there; that comes in the batch being read when the records
    // run out before the batch reaches its size, and otherwise in a batch of its own.
    while input_start < input.input_lengths.len()
        || (reading.reads_past_records() && !reading_ended)
    {
        let (input_end, reached_size) = BatchCount::next_batch(
            count.as_mut(),
            &|record| input_records[record].believed_shared(),
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
                    message: format!("{ENGINE_ERROR} Empty CBlastQueryVector\n").into(),
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
                scope_refusal: false,
            };
        }
        // LOSAT's refusal, which the records alone decide, stops the run at this batch before its
        // search, after the reports of the batches before and the messages of its records, as
        // NCBI's errors of a batch do (`QueryBatches::scope_refusal`; the form of the Owner's
        // decision Q-X2c of 2026-11-16 for TBLASTX's refusals). An error of NCBI's in the search of
        // a batch before comes first (line WQA D10 round 4, 2026-11-20).
        if let Some(scopes) = query_scopes.as_mut() {
            if let Err(refusal) = scopes.add_batch(&input_records[input_start..input_end], start) {
                return QueryBatches {
                    warnings,
                    input_done: input_start,
                    searched_done: start,
                    starts,
                    stop: Some((batch_warnings, refusal)),
                    stop_batch: None,
                    query_scopes,
                    scope_refusal: true,
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
        scope_refusal: false,
    }
}

/// The error of a blastp or tblastn batch whose queries all have no letters: one `Sequence
/// contains no data` per query (`read_query_batches`).
fn no_data_batch_error(batch: &[FastaRecord]) -> Option<anyhow::Error> {
    if !batch.iter().all(|record| record.seq().is_empty()) {
        return None;
    }
    let mut message = format!("{ENGINE_ERROR} ");
    for _ in batch {
        message.push_str("Warning: ");
        message.push_str(NO_DATA_MESSAGE);
        message.push(' ');
    }
    message.push('\n');
    Some(
        crate::cli::NativeError {
            exit: 3,
            message: message.into(),
        }
        .into(),
    )
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
                message: format!("{ENGINE_ERROR} PSI-BLAST only accepts one protein sequence as query\n")
                    .into(),
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup.hpp:190-198
pub(crate) const NO_DATA_MESSAGE: &str = "Sequence contains no data";

#[cfg(test)]
#[path = "query_batch_tests.rs"]
mod tests;
