//! blast_formatter's tabular formats (`-outfmt 6`, 7, 10 and 20, with any field list and
//! delimiter) of one archive object: NCBI's `CBlastFormat` with `PrintProlog`, the prefetch,
//! `PrintOneResultSet` (`x_PrintTabularReport`) of each query and `PrintEpilog`; the rows are
//! made by the one producer of tabular rows of Seq-aligns (`report::tabular_seqalign`) and
//! written by the shared writer (`report::tabular`), decision F3-13 of line RG.
//!
//! NCBI's scope finds each row's query among the request's queries, and its subject in the
//! database that the archive names (its BLAST database data loader: the first OID that the lookup
//! of the Seq-align's subject id gives, `SeqidLookup::seqid_to_oid`, decision F3-14) or among the
//! request's `-subject` Bioseqs. What LOSAT does not find as NCBI's scope would is refused by name
//! before the object's first byte: a subject that the database does not hold (NCBI asks its
//! servers, `Refusal::Database`), and the forms that this driver does not cover
//! (`Refusal::Tabular`): queries given as locations, and ids that NCBI's scope could find
//! elsewhere. The translated programs' Std-seg Seq-aligns are read as the others
//! (`report::tabular_seqalign`, PR 3b of line RG item F3), and `-subject` with `BL2SEQ_LEGACY` in
//! outfmt 7 names each query's subject (`Bl2seqSubjects`).

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use anyhow::Result;

use super::refusals::Refusal;
use super::PROGRAM;
use crate::algorithm::common::database::{
    database_error, parsed_refusal, prefetch_messages, prefetch_warning_of, read_header,
    subject_taxonomy, tabular_names, HeaderReader, ReportedHeader, StoredTitle, TaxdbRun,
    TaxdbSetting,
};
use crate::algorithm::common::tabular::{SubjectId, SubjectValues, TabularSubjects};
use crate::blastinput::app::FormatChoice;
use crate::blastinput::fasta_reader::believed::BelievedIds;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::config::ScoringMatrix;
use crate::record::FastaRecord;
use crate::report::archive::bioseq::SeqLocStrand;
use crate::report::archive::read::options::{BuiltOptions, EProgram, Setter};
use crate::report::archive::read::queries::{bioseq_letters, RenderError};
use crate::report::archive::read::request::{Bioseq, Queries, Request};
use crate::report::archive::read::{seq_align_id, Stop};
use crate::report::db_subject::{HitLoad, LoadsBefore, ParsedSubject, ScopeWalk};
use crate::report::serial::asn_types::view::V;
use crate::report::serial::seqalign_read::{read_set, ReadAlign, ReadStop};
use crate::report::tabular::{
    check_failed_lookups, failed_database_title, failed_record_title, taxdb_warning_line_of,
    HeaderTarget, LookupFailure, SubjectLookup, TabularKind, TabularNames, TabularQuery,
    TabularRows, TabularSpec, TabularTaxonomy, TabularWriter,
};
use crate::report::tabular_seqalign::{
    row_hsp, row_plan, row_scores, subject_scores, GeneticCodes, RowBuffers, RowContext,
    RowLetters, RowPlan, SubjectScores,
};
use crate::seqdb::lookup::SeqidLookup;
use crate::seqdb::{NucleotideAmbiguity, SeqDb, SeqType};
use crate::seqid::defline_set::{compare, match_keys, Sic};
use crate::seqid::model::ObjectId;
use crate::seqid::ParsedSeqId;

/// Where the object's report goes: the `-out` stream and the diagnostics (`archive::Output`,
/// `Diag`).
pub(super) trait ReportSink {
    /// Writes `bytes` to the stream.
    fn write(&mut self, bytes: &[u8]) -> Result<()>;
    /// `flush()` of the stream (the destructor of each `CBlastTabularInfo`).
    fn flush(&mut self) -> Result<()>;
    /// Posts the diagnostic lines `lines` (none when empty) after the stream's bytes before them.
    fn post(&mut self, lines: &[u8]);
}

/// What the object's report reads besides the archive's search.
pub(super) struct ObjectInputs<'a> {
    pub(super) settings: &'a ApplicationSettings,
    /// The options of the first search object (`m_CBOH`).
    pub(super) options: &'a BuiltOptions,
    /// blast_formatter's `-parse_deflines` (`m_BelieveQuery`).
    pub(super) believe: bool,
    pub(super) choice: &'a FormatChoice,
    /// The options' matrix name when the object's formatter is made (`opts.GetMatrixName()`),
    /// `None` for NULL.
    pub(super) matrix_name: Option<&'a [u8]>,
    /// The hit list size (`m_HitlistSize`) as `PruneSeqalign` takes it.
    pub(super) hitlist: u32,
    /// The database of a database object: its name and the restricted reader that opened it.
    pub(super) database: Option<(&'a [u8], &'a SeqDb)>,
    /// The run's record of the taxonomy database (NCBI's static `CTaxDBFileInfo`, the volumes'
    /// lookups of the run).
    pub(super) taxdb: &'a TaxdbRun,
    /// The psiblast round that `PrintOneResultSet` takes (`GetPsiNumberOfIterations` where the
    /// run's task is psiblast; `numeric_limits<unsigned int>::max()`, no line, otherwise), `Err`
    /// where its `Web_StepNumber` is no integer.
    pub(super) iteration: std::result::Result<Option<u32>, ()>,
    /// The run's index of the next subject that outfmt 7's header names where the formatter has
    /// no database name (`x_CreateSubjectBioseq`'s static `subj_index`: `-subject` archives with
    /// `BL2SEQ_LEGACY` set), kept across the objects.
    pub(super) subject_index: &'a std::cell::Cell<u32>,
}

/// The refusal of a form that the tabular formats do not write yet.
fn refused(what: &str) -> anyhow::Error {
    Refusal::Tabular(what.to_string()).into()
}

/// A stop of a Seq-align's read as the run's error.
fn read_stop(stop: ReadStop) -> anyhow::Error {
    match stop {
        ReadStop::Stop(stop) => super::archive::stop_error(stop),
        ReadStop::Unsupported(what) => refused(what),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:183
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3549-3562
// NCBI reference (598d8ae6): c++/src/util/tables/raw_scoremat.c:121-166
/// `GetAsciiProteinMatrix` of the options' matrix name: `NCBISM_GetStandardMatrix` (`blosum` and
/// 45, 50, 62, 80 or 90; `pam` and 250, 30 or 70; a name that begins with `identity`; the case of
/// the prefixes ignored), none for a NULL or blank name and for every other one (the positives
/// are then the identities).
pub(super) fn formatter_matrix(name: Option<&[u8]>) -> Option<ScoringMatrix> {
    let name = name?;
    let lower = name.to_ascii_lowercase();
    let rest = |prefix: &[u8]| lower.starts_with(prefix).then(|| &name[prefix.len()..]);
    match lower.first()? {
        b'b' => match rest(b"blosum")? {
            b"45" => Some(ScoringMatrix::Blosum45),
            b"50" => Some(ScoringMatrix::Blosum50),
            b"62" => Some(ScoringMatrix::Blosum62),
            b"80" => Some(ScoringMatrix::Blosum80),
            b"90" => Some(ScoringMatrix::Blosum90),
            _ => None,
        },
        b'p' => match rest(b"pam")? {
            b"250" => Some(ScoringMatrix::Pam250),
            b"30" => Some(ScoringMatrix::Pam30),
            b"70" => Some(ScoringMatrix::Pam70),
            _ => None,
        },
        b'i' => lower
            .starts_with(b"identity")
            .then_some(ScoringMatrix::Identity),
        _ => None,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:121
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:200
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.cpp:58
/// `GetMatrixName()` of the options of the first search object: the last `MatrixName` setter,
/// else the factory's: BLOSUM62 for the protein handles; none for the nucleotide handles
/// (`SetMatrixName` with a null name leaves the scoring options' NULL) and for RPS-BLAST's (whose
/// defaults set no matrix).
pub(super) fn options_matrix_name(options: &BuiltOptions) -> Option<Vec<u8>> {
    let set = options
        .setters
        .iter()
        .rev()
        .find_map(|setter| match setter {
            Setter::MatrixName(name) => Some(name.as_bytes().to_vec()),
            _ => None,
        });
    set.or_else(|| {
        matches!(
            options.program,
            EProgram::Blastp
                | EProgram::Blastx
                | EProgram::Tblastn
                | EProgram::Tblastx
                | EProgram::PsiBlast
                | EProgram::PsiTblastn
                | EProgram::PhiBlastp
                | EProgram::DeltaBlast
        )
        .then(|| b"BLOSUM62".to_vec())
    })
}

/// The programs whose Seq-aligns the tabular formats write: those of one Dense-seg or Dense-diag
/// per HSP (blastn's tasks, blastp, PSI-BLAST, DELTA-BLAST, PHI-BLAST of a protein query,
/// RPS-BLAST) and the translated ones of Std-segs (BLASTX, TBLASTN, TBLASTX, RPSTBLASTN,
/// PSI-TBLASTN).
fn covered(program: EProgram) -> bool {
    matches!(
        program,
        EProgram::Blastn
            | EProgram::Megablast
            | EProgram::DiscMegablast
            | EProgram::Blastp
            | EProgram::PsiBlast
            | EProgram::DeltaBlast
            | EProgram::PhiBlastp
            | EProgram::RpsBlast
            | EProgram::Blastx
            | EProgram::Tblastn
            | EProgram::Tblastx
            | EProgram::RpsTblastn
            | EProgram::PsiTblastn
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:237
// NCBI reference (598d8ae6): c++/src/algo/blast/api/tblastx_options.cpp:73
/// `m_IsUngappedSearch` (`!opts.GetGappedMode()`): the last `GappedMode` call of the options,
/// else the factory's: TBLASTX's handle is ungapped, every other one gapped.
fn ungapped_search(options: &BuiltOptions) -> bool {
    options
        .setters
        .iter()
        .rev()
        .find_map(|setter| match setter {
            Setter::GappedMode(gapped) => Some(!gapped),
            _ => None,
        })
        .unwrap_or(options.program == EProgram::Tblastx)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:233-234
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:526, 1621
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:556-559
/// The formatter's genetic codes (`GetQueryGeneticCode`, `GetDbGeneticCode`): the last setter of
/// each, else the factory's: `BLAST_GENETIC_CODE` (1) for the query; for the database 1 where the
/// program's subjects are translated (TBLASTN, TBLASTX, PSI-TBLASTN), else 0 (`SetProgram`
/// clears it).
fn genetic_codes(options: &BuiltOptions) -> GeneticCodes {
    let last = |query: bool| {
        options.setters.iter().rev().find_map(|setter| match setter {
            Setter::QueryGeneticCode(code) if query => Some(*code),
            Setter::DbGeneticCode(code) if !query => Some(*code),
            _ => None,
        })
    };
    let translated_subjects = matches!(
        options.program,
        EProgram::Tblastn | EProgram::Tblastx | EProgram::PsiTblastn
    );
    GeneticCodes::new(
        last(true).unwrap_or(1),
        last(false).unwrap_or(i32::from(translated_subjects)),
    )
}

/// Where a query's letters are: its Bioseq in the request, or for a query given as a location
/// the database sequence (its ordinal id) that the scope's loader of the database gives for its
/// Seq-id (the queries' own scope has no loader, blast_formatter.cpp:197).
#[derive(Clone, Copy)]
enum QuerySequence<'r> {
    Bioseq(Bioseq<'r>),
    Database(u32),
}

/// One query of the request as the rows show it: its Bioseq, or the database's Bioseq of a query
/// given as a location (`located_query`).
struct QueryEntry<'r> {
    sequence: QuerySequence<'r>,
    /// The ids of a database query's Bioseq (its definition line's, `GetId()`), by which the rows'
    /// query ids find it; empty for a query Bioseq of the request (its own ids).
    located_ids: Vec<ParsedSeqId>,
    /// `SetQueryId`: the ids after `s_ReplaceLocalId`, whose tabular ids and GI the rows print.
    ids: BelievedIds,
    length: i64,
    /// Whether the query is a nucleotide sequence (`bioseq_protein` false).
    nucleotide: bool,
    /// outfmt 7's `# Query:` text before its spaces are cut (`x_AcknowledgeBlastSequence`).
    title: Vec<u8>,
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:474-504
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4934-4942
/// `s_ReplaceLocalId` of each id of a Bioseq whose first title is `title`: a local id becomes the
/// first word of `GetTitle` (the titles, each followed by a space: the first title's first word,
/// empty where the first title is empty or begins with a space) unless that is empty or
/// `parse_local` is set, where it is its own text.
fn replaced_ids<'i>(
    ids: impl IntoIterator<Item = &'i ParsedSeqId>,
    title: &[u8],
    parse_local: bool,
) -> Vec<ParsedSeqId> {
    let word = &title[..title.iter().position(|&b| b == b' ').unwrap_or(title.len())];
    ids.into_iter()
        .map(|id| match id {
            ParsedSeqId::Local(local) => {
                let token = if word.is_empty() || parse_local {
                    match local {
                        ObjectId::Str(text) => text.clone(),
                        ObjectId::Id(n) => n.to_string().into_bytes(),
                    }
                } else {
                    word.to_vec()
                };
                ParsedSeqId::Local(ObjectId::Str(token))
            }
            other => other.clone(),
        })
        .collect()
}

impl<'r> QueryEntry<'r> {
    /// The query Bioseq `bioseq` (at least one id, `check_queries`) of a formatter whose
    /// `-parse_deflines` is `believe`, with its `# Query:` text (`acknowledged_title`).
    fn new(bioseq: Bioseq<'r>, believe: bool) -> Self {
        let shown = replaced_ids(bioseq.ids(), bioseq.title().unwrap_or_default(), believe);
        QueryEntry {
            sequence: QuerySequence::Bioseq(bioseq),
            located_ids: Vec::new(),
            ids: BelievedIds::of_bioseq(shown, SeqLocStrand::Unknown),
            length: bioseq.length().unwrap_or(0),
            nucleotide: !bioseq_protein(bioseq),
            title: acknowledged_title(bioseq, believe),
        }
    }

    /// The query's Bioseq in the request (`None` for a query given as a location).
    fn bioseq(&self) -> Option<Bioseq<'r>> {
        match self.sequence {
            QuerySequence::Bioseq(bioseq) => Some(bioseq),
            QuerySequence::Database(_) => None,
        }
    }

    fn query(&self) -> TabularQuery<'_> {
        TabularQuery {
            id: &self.ids.tabular[0],
            accession: &self.ids.tabular[1],
            accession_version: &self.ids.tabular[2],
            gi: self.ids.gi,
            length: usize::try_from(self.length).unwrap_or(0),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:629-664
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:669-683
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:713-740
/// outfmt 7's text of the Bioseq `bioseq` (`x_AcknowledgeBlastSequence`, `# Query:` and
/// `# Subject:`) before its spaces are cut, of a formatter whose `-parse_deflines` is `believe`:
/// `GetSeqIdString` (the id that `FindBestChoice` with `WorstRank` picks, by `GetBareId`, unless
/// it is local and the ids are not believed; `[BLAST] LONG_SEQID` is refused by every program), a
/// space and every title joined (`GetSeqDescrString`).
fn acknowledged_title(bioseq: Bioseq<'_>, believe: bool) -> Vec<u8> {
    acknowledged_text(bioseq.ids(), bioseq.titles(), believe)
}

/// `acknowledged_title` of a Bioseq of the ids `ids` and the titles `titles`.
fn acknowledged_text<'i>(
    ids: impl IntoIterator<Item = &'i ParsedSeqId>,
    titles: impl IntoIterator<Item = &'i [u8]>,
    believe: bool,
) -> Vec<u8> {
    let mut title = Vec::new();
    let best = crate::seqid::rank::find_best(ids, crate::seqid::rank::text_score);
    if let Some(id) = best {
        if !matches!(id, ParsedSeqId::Local(_)) || believe {
            crate::seqid::label::bare(id, &mut title);
        }
    }
    title.push(b' ');
    for each in titles {
        title.extend_from_slice(each);
    }
    title
}

/// The items of a list (the queries, the `-subject` Bioseqs, the Seq-align subject ids found) by
/// the keys of their ids (`defline_set::match_keys`): the candidates of a lookup by id, which
/// `compare` confirms, so that a lookup reads the few items that can hold the id, not every item.
#[derive(Default)]
struct IdIndex(HashMap<u64, Vec<usize>>);

impl IdIndex {
    /// The index of the ids of `bioseqs`, by position.
    fn of_bioseqs<'r>(bioseqs: impl IntoIterator<Item = Bioseq<'r>>) -> Self {
        let mut index = IdIndex::default();
        for (at, bioseq) in bioseqs.into_iter().enumerate() {
            for id in bioseq.ids() {
                index.add(at, id);
            }
        }
        index
    }

    /// Adds the item `item` (not below an item added before) with the id `id`.
    fn add(&mut self, item: usize, id: &ParsedSeqId) {
        let mut keys = Vec::new();
        match_keys(id, &mut keys);
        for key in keys {
            let items = self.0.entry(key).or_default();
            if items.last() != Some(&item) {
                items.push(item);
            }
        }
    }

    /// The items whose ids hold `id` (`holds` of an item confirms a candidate), in order, reusing
    /// caller-owned scratch buffers `keys` and `found`.
    fn holders_into(
        &self,
        id: &ParsedSeqId,
        holds: impl Fn(usize) -> bool,
        keys: &mut Vec<u64>,
        found: &mut Vec<usize>,
    ) {
        keys.clear();
        found.clear();
        match_keys(id, keys);
        for key in keys.iter() {
            if let Some(items) = self.0.get(key) {
                found.extend_from_slice(items);
            }
        }
        if found.len() > 1 {
            found.sort_unstable();
            found.dedup();
        }
        found.retain(|&item| holds(item));
    }

    /// The items whose ids hold `id` (`holds` of an item confirms a candidate), in order.
    fn holders(&self, id: &ParsedSeqId, holds: impl Fn(usize) -> bool) -> Vec<usize> {
        let mut keys = Vec::new();
        let mut found = Vec::new();
        self.holders_into(id, holds, &mut keys, &mut found);
        found
    }
}

/// The request's Bioseqs as NCBI's scope finds them by id: the queries' scope, which
/// blast_formatter adds before the database's loader and the `-subject` Bioseqs' scope (at
/// `kSubjectsDataLoaderPriority`), so that a query Bioseq is found first.
struct ScopeBioseqs<'b, 'r> {
    queries: &'b [QueryEntry<'r>],
    /// The ids of the query Bioseqs of the request (none where the queries are locations).
    query_ids: IdIndex,
    /// The ids of the database's Bioseqs of the queries given as locations.
    located_ids: IdIndex,
    subjects: &'b [Bioseq<'r>],
    subject_ids: IdIndex,
}

impl<'r> ScopeBioseqs<'_, 'r> {
    /// The scope of the queries `queries` and the `-subject` Bioseqs `subjects`.
    fn new<'b>(queries: &'b [QueryEntry<'r>], subjects: &'b [Bioseq<'r>]) -> ScopeBioseqs<'b, 'r> {
        let mut located_ids = IdIndex::default();
        for (at, query) in queries.iter().enumerate() {
            for id in &query.located_ids {
                located_ids.add(at, id);
            }
        }
        ScopeBioseqs {
            queries,
            query_ids: IdIndex::of_bioseqs(queries.iter().filter_map(QueryEntry::bioseq)),
            located_ids,
            subjects,
            subject_ids: IdIndex::of_bioseqs(subjects.iter().copied()),
        }
    }

    /// The query Bioseqs of the request that hold `id`, in order (none where the queries are
    /// locations, which only the database's loader gives).
    fn queries_holding(&self, id: &ParsedSeqId) -> Vec<usize> {
        self.query_ids.holders(id, |at| {
            holds(self.queries[at].bioseq().expect("a query Bioseq").ids(), id)
        })
    }

    /// The queries whose Bioseq the scope gives for a row's query id `id`: the request's query
    /// Bioseqs that hold it, or where the queries are locations those whose database Bioseq
    /// does.
    fn row_queries(&self, id: &ParsedSeqId) -> Vec<usize> {
        if self.queries.iter().any(|query| query.bioseq().is_none()) {
            self.located_ids
                .holders(id, |at| holds(&self.queries[at].located_ids, id))
        } else {
            self.queries_holding(id)
        }
    }

    /// The `-subject` Bioseqs that hold `id`, in order.
    fn subjects_holding(&self, id: &ParsedSeqId) -> Vec<usize> {
        self.subject_ids
            .holders(id, |at| holds(self.subjects[at].ids(), id))
    }
}

/// Where a subject of the object's rows comes from.
enum SubjectSource {
    /// A database sequence: its OID, ordinal id and header.
    Database {
        oid: u32,
        ordinal: u32,
        header: ReportedHeader,
    },
    /// The request's `-subject` Bioseq at this index.
    Bioseq(usize),
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:291, 486
    /// The query Bioseq at this index, which NCBI's scope finds before a database sequence or a
    /// `-subject` Bioseq of the same id (`ScopeBioseqs`): a self-hit of a search whose query ids
    /// are the subject's.
    Query(usize),
}

impl SubjectSource {
    /// The sequence of the source: one subject index serves every id of it, as NCBI's data
    /// loader gives the sequence it loaded for each id of its own (`RegisterIds`) and the scope
    /// one Bioseq for each of its ids.
    fn sequence(&self) -> (u8, usize) {
        match *self {
            SubjectSource::Database { oid, .. } => (0, oid as usize),
            SubjectSource::Bioseq(at) => (1, at),
            SubjectSource::Query(at) => (2, at),
        }
    }
}

/// The object's subjects in the order of their first Seq-align, with the Seq-align subject ids
/// found for each, and the database's lookup and header reader kept open for the object.
#[derive(Default)]
struct Subjects<'db> {
    sources: Vec<SubjectSource>,
    /// The Seq-align subject ids found, each with its subject index, by position.
    keys: Vec<(ParsedSeqId, usize)>,
    /// The index of `keys`.
    key_ids: IdIndex,
    /// The subject index of each sequence (`SubjectSource::sequence`).
    by_sequence: HashMap<(u8, usize), usize>,
    /// Index in `keys` of the most recently matched subject id (consecutive HSPs of one subject).
    last: Option<usize>,
    /// One `SeqidLookup`, `HeaderReader` and `ScopeWalk` opened lazily for the whole object.
    db_readers: Option<(SeqidLookup<'db>, HeaderReader<'db>, ScopeWalk)>,
    /// Whether each subject is a nucleotide sequence and its length (`GetSequenceType`,
    /// `GetBioseqLength`), by subject index.
    facts: Vec<(bool, i64)>,
    scratch_keys: Vec<u64>,
    scratch_found: Vec<usize>,
}

impl<'db> Subjects<'db> {
    /// The subject index of the Seq-align subject id `id`, found where no Seq-align before had
    /// it.
    fn index_of(
        &mut self,
        inputs: &ObjectInputs<'db>,
        scope: &ScopeBioseqs<'_, '_>,
        id: &ParsedSeqId,
    ) -> Result<usize> {
        if let Some(at) = self.last {
            if compare(&self.keys[at].0, id) == Sic::Yes {
                return Ok(self.keys[at].1);
            }
        }
        let keys = &self.keys;
        self.key_ids.holders_into(
            id,
            |at| compare(&keys[at].0, id) == Sic::Yes,
            &mut self.scratch_keys,
            &mut self.scratch_found,
        );
        if let Some(&at) = self.scratch_found.first() {
            self.last = Some(at);
            return Ok(self.keys[at].1);
        }
        let source = match scope.queries_holding(id).as_slice() {
            [query] => SubjectSource::Query(*query),
            [_, _, ..] => {
                return Err(refused(
                    "a subject id that names several query Bioseqs of the request",
                ))
            }
            [] => match inputs.database {
                None => match scope.subjects_holding(id).as_slice() {
                    [at] => SubjectSource::Bioseq(*at),
                    _ => {
                        return Err(refused(
                            "a subject id that names no -subject Bioseq of the request, or \
                             several",
                        ))
                    }
                },
                Some((name, db)) => {
                    if self.db_readers.is_none() {
                        let lookup = db
                            .seqid_lookup()
                            .map_err(|error| database_error(error, PROGRAM))?;
                        let headers = HeaderReader::new(db, false);
                        let walk = ScopeWalk::new(LoadsBefore::none(), Some(false));
                        self.db_readers = Some((lookup, headers, walk));
                    }
                    let (lookup, headers, walk) = self.db_readers.as_mut().expect("opened above");
                    database_subject(name, db, lookup, headers, walk, id)?
                }
            },
        };
        let index = match self.by_sequence.get(&source.sequence()) {
            Some(&index) => index,
            None => {
                let facts = source_facts(inputs, scope, &source)?;
                self.by_sequence
                    .insert(source.sequence(), self.sources.len());
                self.sources.push(source);
                self.facts.push(facts);
                self.sources.len() - 1
            }
        };
        let key_at = self.keys.len();
        self.key_ids.add(key_at, id);
        self.keys.push((id.clone(), index));
        self.last = Some(key_at);
        Ok(index)
    }
}

/// Whether the subject `source` is a nucleotide sequence and its length (`GetSequenceType`,
/// `GetBioseqLength`): the database's residue type and the sequence's length, or the Bioseq's.
fn source_facts(
    inputs: &ObjectInputs<'_>,
    scope: &ScopeBioseqs<'_, '_>,
    source: &SubjectSource,
) -> Result<(bool, i64)> {
    Ok(match *source {
        SubjectSource::Database { ordinal, .. } => {
            let (_, db) = inputs.database.expect("a database subject");
            let nucleotide = db.seqtype() != SeqType::Protein;
            let length = if nucleotide {
                db.nucleotide_length(ordinal).map(|length| length as i64)
            } else {
                db.protein_length(ordinal).map(i64::from)
            }
            .map_err(|error| database_error(error, PROGRAM))?;
            (nucleotide, length)
        }
        SubjectSource::Bioseq(at) => {
            let bioseq = scope.subjects[at];
            (!bioseq_protein(bioseq), bioseq.length().unwrap_or(0))
        }
        SubjectSource::Query(at) => {
            let query = &scope.queries[at];
            (query.nucleotide, query.length)
        }
    })
}

/// Whether one of `ids` is `id` for NCBI's scope (`Sic::Yes`).
fn holds<'i>(ids: impl IntoIterator<Item = &'i ParsedSeqId>, id: &ParsedSeqId) -> bool {
    ids.into_iter().any(|held| compare(held, id) == Sic::Yes)
}

// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:393
/// The database sequence that NCBI's BLAST database data loader serves for the Seq-align
/// subject id `id`: the first OID of the restricted database's lookup (`seqid_to_oid`), its
/// header read as the search programs read a reported sequence's (`read_header`). Not found is
/// NCBI's ask of its servers (`Refusal::Database`); a lookup that throws, a header with a
/// UTF8String, and a parsed header whose Seq-align id NCBI's lookup or scope may not give as it is
/// (`ParsedSubject::found`, `unfetchable`, `load_throws`, `scope_finds_none`, `ScopeWalk::load`,
/// or an id other than `id`) are refused.
fn database_subject(
    name: &[u8],
    db: &SeqDb,
    lookup: &SeqidLookup<'_>,
    headers: &HeaderReader<'_>,
    walk: &mut ScopeWalk,
    id: &ParsedSeqId,
) -> Result<SubjectSource> {
    let oid = match lookup.seqid_to_oid(id) {
        Ok(Some(oid)) => oid,
        Ok(None) => {
            return Err(Refusal::Database(String::from_utf8_lossy(name).into_owned()).into())
        }
        Err(_) => {
            return Err(refused(
                "a subject id whose lookup in the database throws, where NCBI BLAST+'s data \
                 loader may give a sequence loaded before",
            ))
        }
    };
    let ordinal = db.ordinal_of(oid);
    if db.utf8_string_path(ordinal).is_some() {
        return Err(refused(
            "a database sequence whose header holds a UTF8String",
        ));
    }
    let header = read_header(headers, ordinal, 0, PROGRAM)?;
    if let ReportedHeader::Parsed(subject) = &header {
        let plain = subject.found().is_none()
            && !subject.unfetchable()
            && !subject.load_throws()
            && !subject.scope_finds_none()
            && matches!(walk.load(subject), HitLoad::Loaded)
            && compare(id, subject.alignment_id()) == Sic::Yes;
        if !plain {
            return Err(refused(
                "a database sequence that NCBI BLAST+'s lookup of its Seq-align id may not give \
                 as it is",
            ));
        }
    }
    Ok(SubjectSource::Database {
        oid,
        ordinal,
        header,
    })
}

/// One subject of the object's rows, by subject index.
struct SubjectEntry<'d> {
    values: SubjectValues<'d>,
    nucleotide: bool,
}

/// The subjects that are Bioseqs of the request (a `-subject` Bioseq, or a query's) as the search
/// programs' records (`TabularSubjects::Records`), by subject index (a database subject's entry
/// is empty): a Bioseq of one local id has the id of `tabular_subject_id` of its text and title,
/// one without local ids its parsed subject.
struct BioseqRecords {
    ids: Vec<SubjectId>,
    lengths: Vec<usize>,
    parsed: Vec<Option<Arc<ParsedSubject>>>,
    records: Vec<FastaRecord>,
    /// The lengths as `TabularSubjects::Records` reads them (`lengths`, by subject index).
    length_of: Box<dyn Fn(usize) -> usize>,
    /// Whether each Bioseq is no nucleotide sequence (`bioseq_protein`), by subject index.
    protein: Vec<bool>,
}

/// Whether `bioseq` is no nucleotide sequence (`CSeq_inst::IsNa` of its `inst.mol` false): NCBI's
/// formatter takes the residue type of each subject from its scope (`GetSequenceType`), not of
/// the search's program, which is the first object's in an archive of several programs.
fn bioseq_protein(bioseq: Bioseq<'_>) -> bool {
    // `Seq-inst.mol`: 1 dna, 2 rna, 4 na.
    !matches!(bioseq.molecule(), 1 | 2 | 4)
}

/// Whether the descriptors `descr` (a Bioseq's or a Bioseq-set's, `None` for none) hold a
/// BioSource or an Org-ref, whose taxid `sequence::GetTaxId` reads.
fn taxon_descriptor(descr: Option<V<'_>>) -> bool {
    descr.is_some_and(|descr| {
        descr
            .items()
            .any(|desc| matches!(desc.choice().0, "source" | "org"))
    })
}

/// Whether the queries' Bioseq-set, or a set inside it, has a BioSource or an Org-ref
/// (`taxon_descriptor`), which `sequence::GetTaxId` of a query Bioseq inside reads.
fn query_sets_taxon(queries: &Queries) -> bool {
    let Some(set) = queries.bioseq_set() else {
        return false;
    };
    let mut todo = vec![set];
    while let Some(set) = todo.pop() {
        if taxon_descriptor(set.get("descr")) {
            return true;
        }
        todo.extend(
            set.field("seq-set")
                .items()
                .filter_map(|entry| match entry.choice() {
                    ("set", inner) => Some(inner),
                    _ => None,
                }),
        );
    }
    false
}

impl BioseqRecords {
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:787-790
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:582-603
    /// The records of the Bioseqs that `sources` name, by subject index, of an object whose
    /// formatter's `-parse_deflines` is `believe` and that searched a database (`database`) or
    /// `-subject` sequences. The local id of a subject id is kept where `s_ReplaceLocalId`'s flag
    /// (`m_ParseSubjectDefline` unset) is set: always for a database object, with `-parse_deflines`
    /// for a `-subject` object (bl2seq's). Where the rows print the taxonomy (`taxonomy`), a
    /// Bioseq with a BioSource or an Org-ref of its own or of a set around it (`sets_taxon`, the
    /// queries') is refused: NCBI shows its taxid (`x_SetTaxInfo`'s `GetTaxId`), which LOSAT does
    /// not read.
    fn of(
        sources: &[SubjectSource],
        scope: &ScopeBioseqs<'_, '_>,
        believe: bool,
        database: bool,
        taxonomy: bool,
        sets_taxon: bool,
    ) -> Result<Self> {
        let mut records = BioseqRecords {
            ids: Vec::new(),
            lengths: Vec::new(),
            parsed: Vec::new(),
            records: Vec::new(),
            length_of: Box::new(|_| 0),
            protein: Vec::new(),
        };
        for source in sources {
            let (bioseq, parse_local, in_taxon_set) = match *source {
                SubjectSource::Database { .. } => {
                    records.ids.push(None);
                    records.parsed.push(None);
                    records
                        .records
                        .push(FastaRecord::new(String::new(), b"", b""));
                    records.lengths.push(0);
                    records.protein.push(false);
                    continue;
                }
                SubjectSource::Bioseq(at) => (scope.subjects[at], believe, false),
                SubjectSource::Query(at) => (
                    scope.queries[at].bioseq().expect("a query Bioseq"),
                    database || believe,
                    sets_taxon,
                ),
            };
            if taxonomy && (in_taxon_set || taxon_descriptor(bioseq.value().get("descr"))) {
                return Err(refused(
                    "a subject Bioseq with a BioSource or an Org-ref (or in a set with one), \
                     whose taxid NCBI BLAST+ shows",
                ));
            }
            let protein = bioseq_protein(bioseq);
            records.protein.push(protein);
            let title = bioseq.title().unwrap_or_default();
            let ids: Vec<ParsedSeqId> = bioseq.ids().cloned().collect();
            let locals = ids
                .iter()
                .filter(|id| matches!(id, ParsedSeqId::Local(_)))
                .count();
            let length = usize::try_from(bioseq.length().unwrap_or(0)).unwrap_or(0);
            match (locals, ids.as_slice()) {
                (0, _) => {
                    let believed = BelievedIds::new(ids, SeqLocStrand::Unknown);
                    records
                        .parsed
                        .push(Some(Arc::new(ParsedSubject::of_record(&believed, title))));
                    records.ids.push(None);
                    records
                        .records
                        .push(FastaRecord::new(String::new(), title, b""));
                }
                (1, [ParsedSeqId::Local(_)]) if parse_local != believe => {
                    return Err(refused(
                        "a self-hit subject (a query Bioseq) of a local id in a database \
                         archive formatted without -parse_deflines, whose subject id and \
                         sallseqid NCBI BLAST+ replaces by different flags",
                    ))
                }
                (1, [ParsedSeqId::Local(local)]) => {
                    let text = match local {
                        ObjectId::Str(text) => String::from_utf8(text.clone())
                            .map_err(|_| refused("a subject Bioseq whose local id is not UTF-8"))?,
                        ObjectId::Id(n) => n.to_string(),
                    };
                    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:522-532
                    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:558-565
                    // `s_ReplaceLocalId` keeps the local id's text where its flag is set: for the
                    // subject id unless `m_ParseSubjectDefline` (a -subject object is bl2seq's,
                    // whose subject deflines are parsed unless `-parse_deflines`), for the `sall`
                    // ids where `m_ParseLocalIds` (`-parse_deflines`): with the option the record
                    // shows its local id as a record of `-parse_deflines` does (`believing`).
                    let id = if believe {
                        crate::report::defline::listed_subject_id(text.as_bytes(), title, protein)
                    } else {
                        crate::report::defline::tabular_subject_id(text.as_bytes(), title, protein)
                    }
                    .ok()
                    .map(Arc::from);
                    records.parsed.push(None);
                    records.ids.push(id);
                    let record = FastaRecord::new(text, title, b"");
                    records.records.push(if believe {
                        record.believing(vec![ParsedSeqId::Local(local.clone())])
                    } else {
                        record
                    });
                }
                _ => {
                    return Err(refused(
                        "a subject Bioseq of a local id and other ids, or of several local ids",
                    ))
                }
            }
            records.lengths.push(length);
        }
        let lengths = records.lengths.clone();
        records.length_of = Box::new(move |s_idx| lengths.get(s_idx).copied().unwrap_or(0));
        Ok(records)
    }
}

/// The values that the fields of `spec` show of each subject (`sources`, by subject index) that
/// a row prints (`printed`), `None` for the others.
fn subject_entries<'d>(
    db: Option<&SeqDb>,
    spec: &TabularSpec,
    sources: &'d [SubjectSource],
    records: &'d BioseqRecords,
    printed: &[bool],
) -> Result<Vec<Option<SubjectEntry<'d>>>> {
    let needs = spec.needs();
    let mut entries = Vec::with_capacity(sources.len());
    for (index, source) in sources.iter().enumerate() {
        if !printed[index] {
            entries.push(None);
            continue;
        }
        // The subject's residue type (`bioseq_protein`): the database's, or the Bioseq's.
        let protein = match (source, db) {
            (SubjectSource::Database { .. }, Some(db)) => db.seqtype() == SeqType::Protein,
            _ => records.protein[index],
        };
        let values = match (source, db) {
            (
                SubjectSource::Database {
                    ordinal, header, ..
                },
                Some(db),
            ) => {
                let length = if protein {
                    db.protein_length(*ordinal).map(|length| length as usize)
                } else {
                    db.nucleotide_length(*ordinal)
                }
                .map_err(|error| database_error(error, PROGRAM))?;
                let oid = db.oid(*ordinal);
                match header {
                    ReportedHeader::Parsed(subject) => {
                        let taxonomy = if needs.taxonomy {
                            TabularTaxonomy::of_deflines(subject.deflines())
                        } else {
                            TabularTaxonomy::default()
                        };
                        SubjectValues::of_parsed(
                            subject,
                            spec,
                            length,
                            std::borrow::Cow::Owned(taxonomy),
                            protein,
                            true,
                            index,
                            |problem| parsed_refusal(oid, &format!("has {problem}"), PROGRAM),
                        )?
                    }
                    ReportedHeader::Stored(StoredTitle { title, source, .. }) => {
                        let taxonomy = if needs.taxonomy {
                            let set = db
                                .defline_set(*ordinal)
                                .map_err(|error| database_error(error, PROGRAM))?;
                            TabularTaxonomy::of_deflines(&set.deflines)
                        } else {
                            TabularTaxonomy::default()
                        };
                        SubjectValues::of_stored(
                            title, source, oid, spec, length, taxonomy, protein, index, PROGRAM,
                        )?
                    }
                }
            }
            (SubjectSource::Bioseq(_) | SubjectSource::Query(_), _) => {
                match records.parsed[index].as_deref() {
                    Some(subject) => SubjectValues::of_parsed(
                        subject,
                        spec,
                        records.lengths[index],
                        std::borrow::Cow::Owned(TabularTaxonomy::default()),
                        protein,
                        false,
                        index,
                        |problem| {
                            refused(&format!(
                                "a subject Bioseq {} that has {problem}",
                                String::from_utf8_lossy(&subject.alignment_fasta())
                            ))
                        },
                    )?,
                    None => TabularSubjects::Records {
                        ids: &records.ids,
                        length: &*records.length_of,
                        parsed: &[],
                        records: &records.records,
                    }
                    .values(index, spec, PROGRAM)?,
                }
            }
            (SubjectSource::Database { .. }, None) => {
                unreachable!("a database subject of a database object")
            }
        };
        entries.push(Some(SubjectEntry {
            values,
            nucleotide: !protein,
        }));
    }
    Ok(entries)
}

/// The failure of the lookups of a subject (`check_failed_lookups`), `None` where they do not
/// fail.
fn lookup_failure(entry: &SubjectEntry<'_>) -> Option<LookupFailure> {
    match entry.values.subject(None, entry.nucleotide, -1, -1).lookup {
        SubjectLookup::Found => None,
        SubjectLookup::Failed(failed) => Some(failed.failure()),
        SubjectLookup::Throws(_) => Some(LookupFailure::thrown()),
    }
}

/// Whether the queries' Bioseq-set, or a set inside it, has a title: `GetTitle` of a query's
/// handle (`s_ReplaceLocalId`) reads the titles of the sets around the Bioseq after its own.
fn set_titles(queries: &Queries) -> bool {
    let Some(set) = queries.bioseq_set() else {
        return false;
    };
    let mut todo = vec![set];
    while let Some(set) = todo.pop() {
        let titled = set
            .get("descr")
            .is_some_and(|descr| descr.items().any(|desc| desc.choice().0 == "title"));
        if titled {
            return true;
        }
        todo.extend(
            set.field("seq-set")
                .items()
                .filter_map(|entry| match entry.choice() {
                    ("set", inner) => Some(inner),
                    _ => None,
                }),
        );
    }
    false
}

/// The object's report (`PrintFormattedOutput` of a tabular format after the results are
/// read): the prolog, the prefetch, each query's report, the epilog. Every refusal comes before
/// the object's first byte.
pub(super) fn print_object(
    inputs: &ObjectInputs<'_>,
    request: &Request,
    sets: &[std::ops::Range<usize>],
    alignments: &[V<'_>],
    sink: &mut dyn ReportSink,
) -> Result<()> {
    let options = inputs.options;
    if !covered(options.program) {
        return Err(refused(
            "a program other than blastn, blastp, blastx, tblastn, tblastx, psiblast, \
             deltablast, PHI-BLAST of a protein query, rpsblast and rpstblastn",
        ));
    }
    if set_titles(&request.queries) {
        return Err(refused("a title on the Bioseq-set of the queries"));
    }
    let spec = TabularSpec::of_formatter(inputs.choice).expect("a tabular format");
    let commented = spec.kind == TabularKind::Commented;
    let program = super::search_program::program_type(options.program);
    let blastn = program.blastn;
    let ungapped = ungapped_search(options);
    let codes = genetic_codes(options);
    let matrix = formatter_matrix(inputs.matrix_name);
    let version = format!(
        "{} {}",
        super::search_program::formatter_program(options.program).to_ascii_uppercase(),
        crate::wording::report_version!()
    );
    let writer = TabularWriter::new(&spec, &version, matrix);
    let database_text: &[u8] = match inputs.database {
        Some((name, _)) => name,
        None => b"User specified sequence set (Input: )",
    };
    let taxdb_found = match inputs.database {
        Some((_, db)) => inputs.taxdb.found(db),
        None => subject_taxonomy(Some(&TaxdbSetting::of(inputs.settings)), [&spec])
            .found
            .unwrap_or(false),
    };
    let warning = (spec.checks_taxdb() && !taxdb_found).then(|| taxdb_warning_line_of(None));
    let post_warning = |sink: &mut dyn ReportSink| {
        if let Some(line) = &warning {
            sink.post(line.as_bytes());
        }
    };
    // `PrintProlog`: outfmt 20's line of the field names, after its formatter's warning.
    let print_prolog = |sink: &mut dyn ReportSink| -> Result<()> {
        if spec.kind == TabularKind::CsvWithHeader {
            post_warning(sink);
            let mut line = Vec::new();
            writer.prolog(&mut line)?;
            sink.write(&line)?;
            sink.flush()?;
        }
        Ok(())
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:430, 455
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:728-734
    // The prefetch, after `PrintProlog` and before any Seq-align is read, converts
    // `PRE_FETCH_SEQS_LIMIT` first: a value that is no integer ends the run after the prolog.
    let prefetch_limit = match crate::blastinput::app::prefetch_seqs_limit() {
        Ok(limit) => limit,
        Err(error) => {
            print_prolog(sink)?;
            return Err(error);
        }
    };

    let queries: Vec<QueryEntry<'_>> = match &request.queries {
        Queries::SeqLocList(_) => located_queries(inputs, &request.queries)?,
        _ => request
            .queries
            .bioseqs()
            .into_iter()
            .map(|bioseq| QueryEntry::new(bioseq, inputs.believe))
            .collect(),
    };
    let subject_bioseqs = request.subject.bioseqs();
    if sets.len() > queries.len() {
        return Err(refused("a result set without its query Bioseq"));
    }
    let scope = ScopeBioseqs::new(&queries, &subject_bioseqs);
    let mut subjects = Subjects::default();

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:658-690, 725-786
    // `BlastFormatter_PreFetchSequenceData`: unless `PRE_FETCH_SEQS_LIMIT` is 0 or below the
    // count of the object's Seq-aligns (`INT_MAX` always loads), the ids and ranges of every
    // Seq-align as the archive holds it (`prefetch_reads`), then the loads of their subjects.
    // One of NCBI's exceptions there ends the reads: NCBI posts its warning and loads nothing
    // (the rows load their subjects).
    let seq_aligns: usize = sets.iter().map(|set| set.len()).sum();
    let prefetches = match prefetch_limit {
        None | Some(i32::MAX) => true,
        Some(0) => false,
        Some(limit) => i64::try_from(seq_aligns).unwrap_or(i64::MAX) <= i64::from(limit),
    };
    let prefetched = prefetches && prefetch_reads(sets, alignments)?;
    if prefetched {
        for set in sets {
            for &align in &alignments[set.clone()] {
                let id = seq_align_id(align, 1).map_err(super::archive::stop_error)?;
                subjects.index_of(inputs, &scope, id)?;
            }
        }
    }
    let prefetched_subjects = subjects.sources.len();

    // The queries' reports as NCBI makes them, read before any byte (every refusal comes
    // first): each query's Seq-aligns (the ungapped split where the result set has
    // alignments), `PruneSeqalign`, `InsertSubjectScores` and each row's reads, up to the first
    // of NCBI's exceptions, which ends the object after the rows before it. The round that
    // `PrintOneResultSet` takes is read for the first query: a `Web_StepNumber` that is no
    // integer throws there (`GetInteger`, a serial exception of the loop; probe E1 probe1 O8).
    let mut checks = PrePass {
        inputs,
        scope: &scope,
        spec: &spec,
        blastn,
        codes: &codes,
        checked: HashSet::new(),
    };
    let mut plans: Vec<QueryPlan<'_>> = Vec::with_capacity(sets.len());
    let mut stop: Option<NcbiStop> = None;
    if inputs.iteration.is_ok() {
        for (q, set) in sets.iter().enumerate() {
            let all = read_set(&alignments[set.clone()], ungapped && !set.is_empty());
            let mut plan = QueryPlan {
                count: all.len(),
                aligns: Vec::new(),
                scores: Vec::new(),
                rows: Vec::new(),
            };
            let read = super::archive::prune_seqalign(&all, inputs.hitlist)
                .map_err(ReadStop::Stop)
                .and_then(|pruned| {
                    let scores = subject_scores(
                        &pruned,
                        Some(queries[q].length),
                        &spec,
                        &inputs.choice.spec,
                        blastn,
                    )?;
                    Ok((pruned, scores))
                });
            match read.map_err(end_of) {
                Ok((pruned, scores)) => (plan.aligns, plan.scores) = (pruned, scores),
                Err(End::Refusal(error)) => return Err(error),
                Err(End::Ncbi(error)) => {
                    plans.push(plan);
                    stop = Some(NcbiStop {
                        error,
                        subject: None,
                    });
                    break;
                }
            }
            for align in &plan.aligns {
                match checks.row(align, &mut subjects) {
                    Ok(row) => plan.rows.push(row),
                    Err(RowEnd::Refusal(error)) => return Err(error),
                    Err(RowEnd::Ncbi(found)) => {
                        stop = Some(found);
                        break;
                    }
                }
            }
            plans.push(plan);
            if stop.is_some() {
                break;
            }
        }
    }

    let mut printed = vec![false; subjects.sources.len()];
    let mut row_subjects: Vec<Vec<usize>> = Vec::with_capacity(plans.len());
    for plan in &plans {
        let mut runs: Vec<usize> = plan.rows.iter().map(|&(_, s_idx)| s_idx as usize).collect();
        runs.dedup();
        row_subjects.push(runs);
    }
    // The row where NCBI's exception comes after its subject's values were read (its letters).
    if let (
        Some(NcbiStop {
            subject: Some(s_idx),
            ..
        }),
        Some(runs),
    ) = (&stop, row_subjects.last_mut())
    {
        if runs.last() != Some(s_idx) {
            runs.push(*s_idx);
        }
    }
    for &s_idx in row_subjects.iter().flatten() {
        printed[s_idx] = true;
    }
    let records = BioseqRecords::of(
        &subjects.sources,
        &scope,
        inputs.believe,
        inputs.database.is_some(),
        spec.needs().taxonomy,
        query_sets_taxon(&request.queries),
    )?;
    let db = inputs.database.map(|(_, db)| db);
    let entries = subject_entries(db, &spec, &subjects.sources, &records, &printed)?;
    let lookup_program = if blastn {
        crate::blastinput::app::Program::Blastn
    } else {
        crate::blastinput::app::Program::Blastp
    };
    check_failed_lookups(
        Some(&spec),
        row_subjects.iter().map(|runs| runs.iter().copied()),
        |s_idx| lookup_failure(entry_of(&entries, s_idx)),
        |s_idx, _failure| match &subjects.sources[s_idx] {
            SubjectSource::Database { oid, .. } => failed_database_title(oid),
            SubjectSource::Bioseq(at) | SubjectSource::Query(at) => {
                failed_record_title(*at, records.records[s_idx].shown_id())
            }
        },
        lookup_program,
        PROGRAM,
        taxdb_found,
    )?;
    let taxonomies: Vec<TabularTaxonomy> = entries
        .iter()
        .map(|entry| match entry {
            Some(entry) => entry
                .values
                .subject(None, entry.nucleotide, -1, -1)
                .taxonomy
                .clone(),
            None => TabularTaxonomy::default(),
        })
        .collect();
    let names: Option<TabularNames> = match inputs.database {
        Some((_, db)) => tabular_names(db, inputs.taxdb, &spec, &taxonomies, PROGRAM)?,
        None => None,
    };
    // The texts of the loads of the database subjects: at the prefetch where it loads them all,
    // else at the rows, which LOSAT does not write.
    let load_texts = match db {
        Some(db) => prefetch_messages(
            db,
            inputs.taxdb,
            subjects.sources.iter().filter_map(|source| match source {
                SubjectSource::Database {
                    ordinal, header, ..
                } => Some((
                    *ordinal,
                    match header {
                        ReportedHeader::Parsed(subject) => Some(&**subject),
                        ReportedHeader::Stored(_) => None,
                    },
                )),
                SubjectSource::Bioseq(_) | SubjectSource::Query(_) => None,
            }),
            PROGRAM,
        )?,
        None => Vec::new(),
    };
    if !load_texts.is_empty() && (!prefetched || subjects.sources.len() != prefetched_subjects) {
        return Err(refused(
            "texts of the taxonomy database at the loads of the subjects, which the rows make \
             where the prefetch does not load them (PRE_FETCH_SEQS_LIMIT, or Seq-aligns that \
             NCBI BLAST+'s prefetch cannot read)",
        ));
    }

    // outfmt 7's `# Subject:` line of each query where the formatter has no database name.
    let subject_headers = if commented
        && inputs.database.is_none()
        && std::env::var_os("BL2SEQ_LEGACY").is_some()
    {
        Some(subject_header_texts(
            inputs,
            &scope,
            &subject_bioseqs,
            plans.len(),
        )?)
    } else {
        None
    };

    print_prolog(sink)?;
    if prefetches && !prefetched {
        sink.post(&prefetch_warning_of(None));
    }
    sink.post(&load_texts);
    let Ok(iteration) = inputs.iteration else {
        return Err(super::archive::invalid_archive());
    };

    let mut letters = LetterCache::default();
    let mut buffers = RowBuffers::default();
    let last = plans.len().saturating_sub(1);
    for (q, plan) in plans.iter().enumerate() {
        post_warning(sink);
        let target = match &subject_headers {
            None => HeaderTarget::Database(database_text),
            Some(texts) => match texts.get(q) {
                Some(text) => HeaderTarget::Subject(text),
                // `CSeqVecSeqInfoSrc::GetId` of an index past the object's subjects.
                None => {
                    return Err(crate::blastinput::app::engine_error(
                        "Index out of range for id retrieval",
                    ))
                }
            },
        };
        let mut out = Vec::new();
        writer.header_of(
            &mut out,
            iteration,
            &queries[q].title,
            target,
            Some(plan.count),
        )?;
        let mut rows = TabularRows::default();
        for (n, &(q_idx, s_idx)) in plan.rows.iter().enumerate() {
            let (q_idx, s_idx) = (q_idx as usize, s_idx as usize);
            let (qcovs, qcovus) = plan.scores.get(n).copied().unwrap_or((-1, -1));
            let align = &plan.aligns[n];
            let query = &queries[q_idx];
            let ctx = row_context(&spec, blastn, &codes, query, subjects.facts[s_idx]);
            let planned = row_scores(align, &spec)
                .and_then(|scores| row_plan(align, scores, &ctx))
                .map_err(read_stop)?;
            let sequences = match letters.row(
                &planned,
                &ctx,
                &scope,
                inputs,
                &subjects.sources[s_idx],
                q_idx,
            ) {
                Ok(sequences) => sequences,
                Err(end) => return Err(end_at_row(sink, &out, end)),
            };
            let hsp = row_hsp(&planned, &ctx, &sequences, &mut buffers).map_err(read_stop)?;
            let entry = entry_of(&entries, s_idx);
            let mut subject = entry
                .values
                .subject(names.as_ref(), entry.nucleotide, qcovs, qcovus);
            let kept = rows.last_found().map(|kept| {
                let kept = entry_of(&entries, kept);
                kept.values
                    .subject(names.as_ref(), kept.nucleotide, qcovs, qcovus)
            });
            // A row that ends at the query block's `catch` (`RowPlan::leaves`) shows the
            // Seq-align's query id as it is (`SetQueryId` of `GetSeq_id(0)`), no subject strand
            // (`m_SubjectStrand` empty: `N/A`), and the subject's length only where the subject
            // block read it.
            let caught;
            let shown_query = if planned.leaves {
                subject.nucleotide = false;
                if !spec.reads_subject_bioseq() {
                    subject.length = 0;
                }
                let id = align.seq_id(0).map_err(super::archive::stop_error)?;
                caught = BelievedIds::of_bioseq(vec![id.clone()], SeqLocStrand::Unknown);
                TabularQuery {
                    id: &caught.tabular[0],
                    accession: &caught.tabular[1],
                    accession_version: &caught.tabular[2],
                    gi: caught.gi,
                    length: usize::try_from(query.length).unwrap_or(0),
                }
            } else {
                query.query()
            };
            // The texts of the row's lookups of the taxonomy database follow the rows before it.
            post_row_texts(sink, &mut out, names.as_ref(), &spec, &writer, &subject)?;
            writer.row(
                &mut out,
                &mut rows,
                &shown_query,
                &subject,
                kept.as_ref(),
                &hsp,
            );
        }
        if q == last {
            if let Some(NcbiStop { error, subject }) = stop.take() {
                // The row whose letters throw read its subject's values first.
                if let Some(s_idx) = subject {
                    let entry = entry_of(&entries, s_idx);
                    let shown = entry
                        .values
                        .subject(names.as_ref(), entry.nucleotide, -1, -1);
                    post_row_texts(sink, &mut out, names.as_ref(), &spec, &writer, &shown)?;
                }
                // `CBlastTabularInfo`'s destructor flushes the rows before the exception.
                sink.write(&out)?;
                sink.flush()?;
                return Err(error);
            }
        }
        sink.write(&out)?;
        sink.flush()?;
    }
    // `PrintEpilog`: outfmt 7's count of the queries formatted.
    if commented {
        post_warning(sink);
        let mut line = Vec::new();
        writer.epilog(&mut line, plans.len())?;
        sink.write(&line)?;
        sink.flush()?;
    }
    Ok(())
}

/// The refusal of queries given as locations that NCBI BLAST+ asks its servers for.
const LOCATED_ELSEWHERE: &str =
    "queries given as locations that the database does not hold, whose sequences NCBI BLAST+ \
     asks its servers for";

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:180-221
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1310-1470
/// The queries of a `seq-loc-list` (each a whole sequence or an interval, `check_queries`): the
/// queries' scope has no loader, so the scope's loader of the database gives each Seq-id's
/// sequence (`database_subject`, the subjects' lookup), whose Bioseq is the query: its
/// definition line's ids and title, the whole sequence's length and letters, the database's
/// residue type. A Seq-id that the database does not hold, or a `-subject` object, is NCBI's ask
/// of its servers (refused). A sequence of several definition lines, of taxids (the texts of
/// its load) or of a title or id that the reader repairs is refused.
fn located_queries<'r>(
    inputs: &ObjectInputs<'_>,
    queries: &'r Queries,
) -> Result<Vec<QueryEntry<'r>>> {
    let Some((name, db)) = inputs.database else {
        return Err(refused(LOCATED_ELSEWHERE));
    };
    let lookup = db
        .seqid_lookup()
        .map_err(|error| database_error(error, PROGRAM))?;
    let headers = HeaderReader::new(db, false);
    let mut walk = ScopeWalk::new(LoadsBefore::none(), Some(false));
    let protein = db.seqtype() == SeqType::Protein;
    let mut out = Vec::new();
    for loc in queries.seq_locs() {
        let id = match loc.choice() {
            ("whole", id) => id.seq_id(),
            ("int", int) => int.field("id").seq_id(),
            _ => return Err(refused("a query location other than a whole sequence or an interval")),
        };
        if !matches!(lookup.seqid_to_oid(id), Ok(Some(_))) {
            return Err(refused(LOCATED_ELSEWHERE));
        }
        let SubjectSource::Database {
            ordinal, header, ..
        } = database_subject(name, db, &lookup, &headers, &mut walk, id)?
        else {
            unreachable!("a database sequence");
        };
        let ReportedHeader::Parsed(subject) = header else {
            return Err(refused(
                "a query location whose database sequence has no parsed definition line",
            ));
        };
        let [defline] = subject.deflines() else {
            return Err(refused(
                "a query location whose database sequence has several definition lines",
            ));
        };
        if !subject.volume_lookups().is_empty()
            || defline.replaced_title_bytes > 0
            || defline.replaced_id_bytes > 0
        {
            return Err(refused(
                "a query location whose database sequence has taxids or a definition line that \
                 NCBI BLAST+ repairs",
            ));
        }
        let title = defline.title.as_deref().unwrap_or_default();
        let length = if protein {
            db.protein_length(ordinal).map(i64::from)
        } else {
            db.nucleotide_length(ordinal).map(|length| length as i64)
        }
        .map_err(|error| database_error(error, PROGRAM))?;
        let shown = replaced_ids(&defline.ids, title, inputs.believe);
        out.push(QueryEntry {
            sequence: QuerySequence::Database(ordinal),
            located_ids: defline.ids.clone(),
            ids: BelievedIds::of_bioseq(shown, SeqLocStrand::Unknown),
            length,
            nucleotide: !protein,
            title: acknowledged_text(&defline.ids, [title], inputs.believe),
        });
    }
    Ok(out)
}

/// The letters of the database sequence `ordinal` as a row reads them (IUPAC, or the protein's
/// FASTA letters).
fn database_letters(
    db: &SeqDb,
    ordinal: u32,
) -> std::result::Result<Vec<u8>, crate::seqdb::SeqDbError> {
    if db.seqtype() == SeqType::Protein {
        let block = db.read_sequence(ordinal)?;
        return Ok(crate::sequence::stdaa_letters::protein_fasta_letters(
            block.protein(0),
        ));
    }
    let ambiguity = db.nucleotide_ambiguity(ordinal)?;
    let length = ambiguity.length();
    db.nucleotide_iupac_range_with(ordinal, 0, length, &ambiguity)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:797-808, 1024-1047
// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:232-246, 276-279
// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqinfosrc_seqvec.cpp:61-77
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:700-746
/// `x_CreateSubjectBioseq` for the headers of `count` queries of a `-subject` object whose
/// formatter has no database name (`BL2SEQ_LEGACY` set: the subjects' adapter is not in its
/// scan mode): the text of `AcknowledgeBlastSubject` (`acknowledged_title`) of the Bioseq that
/// the scope finds for the id of the subject at the run's index (`FindBestChoice` of its ids by
/// `BestRank`; the query Bioseqs before the subjects), the index then moving on and starting
/// again at 0 past the object's last subject. The texts stop before the query whose index is
/// past the object's subjects (an index kept from an object of more subjects), where NCBI's
/// `GetId` throws (the run's engine error).
fn subject_header_texts(
    inputs: &ObjectInputs<'_>,
    scope: &ScopeBioseqs<'_, '_>,
    subjects: &[Bioseq<'_>],
    count: usize,
) -> Result<Vec<Vec<u8>>> {
    if subjects.is_empty() {
        return Err(refused("-subject archives without subject Bioseqs"));
    }
    let mut texts = Vec::with_capacity(count);
    let mut index = inputs.subject_index.get();
    for _ in 0..count {
        let Some(&subject) = subjects.get(index as usize) else {
            break;
        };
        index += 1;
        if index as usize >= subjects.len() {
            index = 0;
        }
        let Some(id) = subject.best_id() else {
            return Err(refused("a subject Bioseq without an id that FindBestChoice ranks"));
        };
        let shown = match scope.queries_holding(id).as_slice() {
            [at] => scope.queries[*at].bioseq().expect("a query Bioseq"),
            [] => match scope.subjects_holding(id).as_slice() {
                [at] => scope.subjects[*at],
                _ => return Err(refused("a subject id that names several -subject Bioseqs")),
            },
            _ => return Err(refused("a subject id that names several query Bioseqs")),
        };
        texts.push(acknowledged_title(shown, inputs.believe));
    }
    inputs.subject_index.set(index);
    Ok(texts)
}

/// The entry of the printed subject `s_idx`.
fn entry_of<'e, 'd>(entries: &'e [Option<SubjectEntry<'d>>], s_idx: usize) -> &'e SubjectEntry<'d> {
    entries[s_idx].as_ref().expect("a printed subject")
}

/// Posts the texts of the lookups of the taxonomy database of a row's subject `subject`, after
/// the rows before it (`out`, written first).
fn post_row_texts(
    sink: &mut dyn ReportSink,
    out: &mut Vec<u8>,
    names: Option<&TabularNames>,
    spec: &TabularSpec,
    writer: &TabularWriter,
    subject: &crate::report::tabular::TabularSubject<'_>,
) -> Result<()> {
    if let (Some(names), true) = (names, writer.sets_subject(subject)) {
        let text = names.row_messages(spec, subject.taxonomy);
        if !text.is_empty() {
            sink.write(out)?;
            out.clear();
            sink.post(&text);
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:658-690
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:153-160, 294-316
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:165-270
/// `s_ExtractSeqidsAndRanges` over the Seq-aligns of `sets` (`alignments`) as the archive holds
/// them: `GetSeq_id` and `GetSeqRange` of the subject's row and `GetSeqStrand` of the query's (a
/// Dense-seg's or Dense-diags', as the tabular formats read them, `ReadAlign::dense_seg`; a
/// Std-seg's own range, `ReadAlign::std_range`, after which its strands read).
/// `Ok(false)` where one of NCBI's exceptions ends them (the prefetch's warning); a stop that
/// NCBI has no result for, and a form that LOSAT does not read as NCBI, is refused.
fn prefetch_reads(sets: &[std::ops::Range<usize>], alignments: &[V<'_>]) -> Result<bool> {
    for set in sets {
        for align in read_set(&alignments[set.clone()], false) {
            let read = seq_align_id(align.align, 1)
                .map_err(ReadStop::Stop)
                .and_then(|_| match align.variant() {
                    "denseg" | "dendiag" => align.dense_seg(false).map(|_| ()),
                    "std" => align.std_range(1).map(|_| ()),
                    "disc" => Err(ReadStop::Unsupported(DISC_SEG)),
                    _ => Err(ReadStop::Unsupported(
                        "segments other than a Dense-seg, Dense-diags or Std-segs",
                    )),
                });
            match read.map_err(end_of) {
                Ok(()) => {}
                Err(End::Ncbi(_)) => return Ok(false),
                Err(End::Refusal(error)) => return Err(error),
            }
        }
    }
    Ok(true)
}

/// The refusal of a Seq-align of discontinuous segments, which no BLAST search writes.
const DISC_SEG: &str = "discontinuous Seq-aligns (Disc-seg)";

/// How a stop on a Seq-align ends the object (`archive::stop_error`): NCBI's exceptions
/// (`Ncbi`, a `CException` or a serial exception) where NCBI meets them, after the bytes before;
/// the stops that NCBI has no result for (a null pointer, an `abort`, memory the input does not
/// decide) and the forms LOSAT does not read as NCBI are refused, before any byte.
enum End {
    Ncbi(anyhow::Error),
    Refusal(anyhow::Error),
}

/// The end of `stop` (`End`).
fn end_of(stop: ReadStop) -> End {
    match stop {
        ReadStop::Stop(
            stop @ (Stop::Exception(_) | Stop::InvalidSelection { .. } | Stop::Unassigned { .. }),
        ) => End::Ncbi(super::archive::stop_error(stop)),
        other => End::Refusal(read_stop(other)),
    }
}

/// NCBI's exception that ends the object after the rows of the last query read
/// (`QueryPlan::rows`), and the subject of the row that throws where it read its values first.
struct NcbiStop {
    error: anyhow::Error,
    subject: Option<usize>,
}

/// How the reads of a row end the object (`End`, with the row's subject for `NcbiStop`).
enum RowEnd {
    Ncbi(NcbiStop),
    Refusal(anyhow::Error),
}

impl From<End> for RowEnd {
    fn from(end: End) -> Self {
        match end {
            End::Ncbi(error) => RowEnd::Ncbi(NcbiStop {
                error,
                subject: None,
            }),
            End::Refusal(error) => RowEnd::Refusal(error),
        }
    }
}

/// What the pre-pass found of one query's report: the count of its Seq-aligns as the formatter
/// reads them (outfmt 7's line), those `PruneSeqalign` keeps, their coverages, and the query and
/// subject index of each row that NCBI prints (to its exception).
struct QueryPlan<'a> {
    count: usize,
    aligns: Vec<ReadAlign<'a>>,
    scores: SubjectScores,
    rows: Vec<(u32, u32)>,
}

/// The reads of the rows before any byte (`PrePass::row`).
struct PrePass<'p, 'b, 'r> {
    inputs: &'p ObjectInputs<'p>,
    scope: &'p ScopeBioseqs<'b, 'r>,
    spec: &'p TabularSpec,
    blastn: bool,
    codes: &'p GeneticCodes,
    /// The Bioseqs whose letters read (`LetterCache::sequence` keys).
    checked: HashSet<(u8, usize)>,
}

impl<'p> PrePass<'p, '_, '_> {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:749-1096
    // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:893-904
    /// The reads of the row of `align` in `SetFields`' order: its scores, its query (`GetSeq_id`
    /// of row 0, among the request's queries), its subject (row 1, `Subjects::index_of`), the
    /// coverage and Dense-seg that it reads (`row_plan`; a Std-seg's Dense-seg is made after its
    /// subject's values are read), and where it reads letters (`GetWholeAlnSeqString` of row 0,
    /// then of row 1) those of each row's Bioseq (a query's, or a subject that is a Bioseq of the
    /// request; each Bioseq once; a database subject's at the row) and their translation
    /// (`RowPlan::translation_stop`). The row's query and subject index.
    fn row(
        &mut self,
        align: &ReadAlign<'_>,
        subjects: &mut Subjects<'p>,
    ) -> Result<(u32, u32), RowEnd> {
        if align.variant() == "disc" {
            return Err(RowEnd::Refusal(refused(DISC_SEG)));
        }
        let scores = row_scores(align, self.spec).map_err(end_of)?;
        let query_id = align
            .seq_id(0)
            .map_err(|stop| end_of(ReadStop::Stop(stop)))?;
        let q_idx = match self.scope.row_queries(query_id).as_slice() {
            [at] => *at,
            _ => {
                return Err(RowEnd::Refusal(refused(
                    "a query id that names no query Bioseq of the request, or several Bioseqs",
                )))
            }
        };
        let id = align
            .seq_id(1)
            .map_err(|stop| end_of(ReadStop::Stop(stop)))?;
        let s_idx = subjects
            .index_of(self.inputs, self.scope, id)
            .map_err(RowEnd::Refusal)?;
        let query = &self.scope.queries[q_idx];
        let ctx = row_context(self.spec, self.blastn, self.codes, query, subjects.facts[s_idx]);
        let at_row = |error| {
            RowEnd::Ncbi(NcbiStop {
                error,
                subject: Some(s_idx),
            })
        };
        let translated = align.variant() == "std";
        let plan = row_plan(align, scores, &ctx).map_err(|stop| match end_of(stop) {
            End::Ncbi(error) if translated => at_row(error),
            other => other.into(),
        })?;
        if plan.subject_range(&ctx).is_none() {
            return Ok((q_idx as u32, s_idx as u32));
        }
        let subject = match subjects.sources[s_idx] {
            SubjectSource::Database { .. } => None,
            SubjectSource::Bioseq(at) => Some((1, at, self.scope.subjects[at], "a -subject Bioseq")),
            SubjectSource::Query(at) => Some((
                2,
                at,
                self.scope.queries[at].bioseq().expect("a query Bioseq"),
                "a query Bioseq",
            )),
        };
        // A database query's letters do not throw (`LetterCache::row` reads them).
        let query_bioseq = query
            .bioseq()
            .map(|bioseq| (2, q_idx, bioseq, "a query Bioseq"));
        let rows = [query_bioseq, subject];
        for (row, bioseq) in rows.into_iter().enumerate() {
            if let Some((kind, at, bioseq, what)) = bioseq {
                if !self.checked.contains(&(kind, at)) {
                    match read_bioseq(bioseq, what) {
                        Ok(_) => {
                            self.checked.insert((kind, at));
                        }
                        Err(LetterEnd::Exception(exception)) => {
                            return Err(at_row(crate::cli::NativeError::exception(&exception).into()))
                        }
                        Err(LetterEnd::Refusal(error)) => return Err(RowEnd::Refusal(error)),
                    }
                }
            }
            match plan.translation_stop(&ctx, row) {
                Ok(None) => {}
                Ok(Some(exception)) => {
                    return Err(at_row(crate::cli::NativeError::exception(&exception).into()))
                }
                Err(stop) => return Err(end_of(stop).into()),
            }
        }
        Ok((q_idx as u32, s_idx as u32))
    }
}

/// What the row of the query `query` and a subject of `facts` (`Subjects::facts`) reads besides
/// its Seq-align (`RowContext`).
fn row_context<'s>(
    spec: &'s TabularSpec,
    blastn: bool,
    codes: &'s GeneticCodes,
    query: &QueryEntry<'_>,
    (subject_nucleotide, subject_length): (bool, i64),
) -> RowContext<'s> {
    RowContext {
        spec,
        no_fetch: blastn,
        query_length: query.length,
        subject_length,
        nucleotide: [query.nucleotide, subject_nucleotide],
        codes,
    }
}

/// How the read of a row's letters ends the run: NCBI's exception (`CSeqVector_CI`'s, as
/// `report::tabular::query_sequence_end`), or LOSAT's refusal.
enum LetterEnd {
    Exception(crate::ncbi_exception::NcbiException),
    Refusal(anyhow::Error),
}

/// The end of the run at a row whose letters do not read: the rows before it (`before`) are
/// written and flushed, then NCBI's exception (exit status 255); a refusal writes nothing.
fn end_at_row(sink: &mut dyn ReportSink, before: &[u8], end: LetterEnd) -> anyhow::Error {
    match end {
        LetterEnd::Exception(exception) => {
            if let Err(error) = sink.write(before).and_then(|()| sink.flush()) {
                return error;
            }
            crate::cli::NativeError::exception(&exception).into()
        }
        LetterEnd::Refusal(error) => error,
    }
}

/// The letters of a Bioseq as `CAlnVec` reads them (`bioseq_letters`); `what` names it in a
/// refusal.
fn read_bioseq(bioseq: Bioseq<'_>, what: &str) -> std::result::Result<Vec<u8>, LetterEnd> {
    match bioseq_letters(bioseq) {
        Ok(letters) => Ok(letters),
        Err(RenderError::Exception(exception)) => Err(LetterEnd::Exception(exception)),
        Err(RenderError::Unsupported(part)) => {
            Err(LetterEnd::Refusal(refused(&format!("{what} with {part}"))))
        }
        Err(RenderError::Write(error)) => Err(LetterEnd::Refusal(error.into())),
    }
}

/// The letters that the rows read, one sequence at a time: the current query's, the current
/// Bioseq subject's or database protein's whole, and of a database nucleotide only the row's
/// range (with the ambiguity runs of the current one), so that the object keeps no more than one
/// sequence of each kind.
#[derive(Default)]
struct LetterCache {
    query: Option<(usize, Vec<u8>)>,
    /// The letters of the current subject read whole: its key (`(kind, index)`: 0 a database
    /// protein's ordinal id, 1 a `-subject` Bioseq, 2 a query Bioseq) and its letters.
    whole: Option<((u8, usize), Vec<u8>)>,
    ambiguity: Option<(u32, NucleotideAmbiguity)>,
    /// A database nucleotide's letters of the row's range, and where they start.
    range: (Vec<u8>, u32),
}

impl LetterCache {
    /// The letters of the row whose plan is `plan` (`RowLetters`), of the query `q_idx` and the
    /// subject `source`; none where the row reads none.
    fn row<'c>(
        &'c mut self,
        plan: &RowPlan,
        ctx: &RowContext<'c>,
        scope: &ScopeBioseqs<'_, '_>,
        inputs: &ObjectInputs<'_>,
        source: &SubjectSource,
        q_idx: usize,
    ) -> std::result::Result<RowLetters<'c>, LetterEnd> {
        let Some((from, to)) = plan.subject_range(ctx) else {
            return Ok(RowLetters {
                query: &[],
                subject: &[],
                subject_from: 0,
                codes: ctx.codes,
            });
        };
        let db_failed = |error| LetterEnd::Refusal(database_error(error, PROGRAM));
        if self.query.as_ref().is_none_or(|(at, _)| *at != q_idx) {
            let letters = match scope.queries[q_idx].sequence {
                QuerySequence::Bioseq(bioseq) => read_bioseq(bioseq, "a query Bioseq")?,
                QuerySequence::Database(ordinal) => {
                    let (_, db) = inputs.database.expect("a database query");
                    database_letters(db, ordinal).map_err(db_failed)?
                }
            };
            self.query = Some((q_idx, letters));
        }
        let whole = match *source {
            SubjectSource::Database { ordinal, .. } => {
                let (_, db) = inputs.database.expect("a database subject");
                if db.seqtype() == SeqType::Protein {
                    Some(((0, ordinal as usize), None))
                } else {
                    if self.ambiguity.as_ref().is_none_or(|(at, _)| *at != ordinal) {
                        let ambiguity = db.nucleotide_ambiguity(ordinal).map_err(db_failed)?;
                        self.ambiguity = Some((ordinal, ambiguity));
                    }
                    let (_, ambiguity) = self.ambiguity.as_ref().expect("read above");
                    // The letters of the range that the subject holds (`RowLetters` reads none
                    // past its end).
                    let length = ambiguity.length();
                    let from = (from as usize).min(length);
                    let end = (to as usize).saturating_add(1).min(length);
                    let letters = if from < end {
                        db.nucleotide_iupac_range_with(ordinal, from, end, ambiguity)
                            .map_err(db_failed)?
                    } else {
                        Vec::new()
                    };
                    self.range = (letters, from as u32);
                    None
                }
            }
            SubjectSource::Bioseq(at) => Some(((1, at), Some(scope.subjects[at]))),
            SubjectSource::Query(at) => Some((
                (2, at),
                Some(scope.queries[at].bioseq().expect("a query Bioseq")),
            )),
        };
        let (subject, subject_from): (&[u8], u32) = match whole {
            Some((key, bioseq)) => {
                if self.whole.as_ref().is_none_or(|(at, _)| *at != key) {
                    let letters = match bioseq {
                        Some(bioseq) => read_bioseq(
                            bioseq,
                            if key.0 == 1 {
                                "a -subject Bioseq"
                            } else {
                                "a query Bioseq"
                            },
                        )?,
                        None => {
                            let (_, db) = inputs.database.expect("a database subject");
                            let block = db.read_sequence(key.1 as u32).map_err(db_failed)?;
                            crate::sequence::stdaa_letters::protein_fasta_letters(block.protein(0))
                        }
                    };
                    self.whole = Some((key, letters));
                }
                (&self.whole.as_ref().expect("read above").1, 0)
            }
            None => (&self.range.0, self.range.1),
        };
        Ok(RowLetters {
            query: &self.query.as_ref().expect("read above").1,
            subject,
            subject_from,
            codes: ctx.codes,
        })
    }
}

#[cfg(test)]
#[path = "tabular_tests.rs"]
mod tests;
