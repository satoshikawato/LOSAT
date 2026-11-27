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
//! (`Refusal::NotCovered`): the translated programs, queries given as locations, `-subject` with
//! `BL2SEQ_LEGACY` in outfmt 7, and ids that NCBI's scope could find elsewhere.

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
    LookupFailure, SubjectLookup, TabularKind, TabularNames, TabularQuery, TabularRows,
    TabularSpec, TabularTaxonomy, TabularWriter,
};
use crate::report::tabular_seqalign::{
    reads_letters, row_hsp, row_plan, row_scores, subject_scores, RowBuffers, RowContext,
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
    /// Sets the stream's `badbit` (an insertion of a null `const char*`, the pairwise epilog's
    /// matrix name, `report::ostream::ReportWrite::badbit`), which throws: the run's error.
    fn badbit(&mut self) -> anyhow::Error;
    /// The run's error of a write that failed in the alignments of outfmt 0-4
    /// (`report::ostream::Site::Alignments`; `flushing`: at a flush), after NCBI's posts there.
    fn alignments_failed(&mut self, error: anyhow::Error, flushing: bool) -> anyhow::Error;
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
}

/// The refusal of a form that the tabular formats do not write yet.
fn refused(what: &str) -> anyhow::Error {
    Refusal::NotCovered {
        format: super::refusals::TABULAR_FORMATS,
        what: what.to_string(),
    }
    .into()
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
/// RPS-BLAST); the translated programs come with their Std-segs (PR 3b of line RG item F3).
pub(super) fn covered(program: EProgram) -> bool {
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
    )
}

/// One query Bioseq of the request as the rows show it.
pub(super) struct QueryEntry<'r> {
    pub(super) bioseq: Bioseq<'r>,
    /// `SetQueryId`: the ids after `s_ReplaceLocalId`, whose tabular ids and GI the rows print.
    ids: BelievedIds,
    pub(super) length: i64,
    /// outfmt 7's `# Query:` text before its spaces are cut (`x_AcknowledgeBlastSequence`), the
    /// pairwise report's `Query=` text.
    pub(super) title: Vec<u8>,
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
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:629-664
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:669-683
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:713-740
    /// The query Bioseq `bioseq` (at least one id, `check_queries`) of a formatter whose
    /// `-parse_deflines` is `believe`. The `# Query:` text is `GetSeqIdString` (the id that
    /// `FindBestChoice` with `WorstRank` picks, by `GetBareId`, unless it is local and the ids are
    /// not believed; `[BLAST] LONG_SEQID` is refused by every program), a space and every title
    /// joined (`GetSeqDescrString`).
    pub(super) fn new(bioseq: Bioseq<'r>, believe: bool) -> Self {
        let titles: Vec<u8> = bioseq.titles().flatten().copied().collect();
        let shown = replaced_ids(bioseq.ids(), bioseq.title().unwrap_or_default(), believe);
        let mut title = Vec::new();
        let best = crate::seqid::rank::find_best(bioseq.ids(), crate::seqid::rank::text_score);
        if let Some(id) = best {
            if !matches!(id, ParsedSeqId::Local(_)) || believe {
                crate::seqid::label::bare(id, &mut title);
            }
        }
        title.push(b' ');
        title.extend_from_slice(&titles);
        QueryEntry {
            bioseq,
            ids: BelievedIds::of_bioseq(shown, SeqLocStrand::Unknown),
            length: bioseq.length().unwrap_or(0),
            title,
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

/// The items of a list (the queries, the `-subject` Bioseqs, the Seq-align subject ids found) by
/// the keys of their ids (`defline_set::match_keys`): the candidates of a lookup by id, which
/// `compare` confirms, so that a lookup reads the few items that can hold the id, not every item.
#[derive(Default)]
pub(super) struct IdIndex(HashMap<u64, Vec<usize>>);

impl IdIndex {
    /// The index of the ids of `bioseqs`, by position.
    pub(super) fn of_bioseqs<'r>(bioseqs: impl IntoIterator<Item = Bioseq<'r>>) -> Self {
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
pub(super) struct ScopeBioseqs<'b, 'r> {
    pub(super) queries: &'b [QueryEntry<'r>],
    pub(super) query_ids: IdIndex,
    pub(super) subjects: &'b [Bioseq<'r>],
    pub(super) subject_ids: IdIndex,
}

impl ScopeBioseqs<'_, '_> {
    /// The queries that hold `id`, in order.
    pub(super) fn queries_holding(&self, id: &ParsedSeqId) -> Vec<usize> {
        self.query_ids
            .holders(id, |at| holds(self.queries[at].bioseq.ids(), id))
    }

    /// The `-subject` Bioseqs that hold `id`, in order.
    fn subjects_holding(&self, id: &ParsedSeqId) -> Vec<usize> {
        self.subject_ids
            .holders(id, |at| holds(self.subjects[at].ids(), id))
    }
}

/// Where a subject of the object's rows comes from.
pub(super) enum SubjectSource {
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
pub(super) struct Subjects<'db> {
    pub(super) sources: Vec<SubjectSource>,
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
    scratch_keys: Vec<u64>,
    scratch_found: Vec<usize>,
}

impl<'db> Subjects<'db> {
    /// The subject index of the Seq-align subject id `id`, found where no Seq-align before had
    /// it.
    pub(super) fn index_of(
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
        let sources = &mut self.sources;
        let index = *self
            .by_sequence
            .entry(source.sequence())
            .or_insert_with(|| {
                sources.push(source);
                sources.len() - 1
            });
        let key_at = self.keys.len();
        self.key_ids.add(key_at, id);
        self.keys.push((id.clone(), index));
        self.last = Some(key_at);
        Ok(index)
    }
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
pub(super) struct BioseqRecords {
    ids: Vec<SubjectId>,
    pub(super) lengths: Vec<usize>,
    pub(super) parsed: Vec<Option<Arc<ParsedSubject>>>,
    pub(super) records: Vec<FastaRecord>,
    /// The lengths as `TabularSubjects::Records` reads them (`lengths`, by subject index).
    length_of: Box<dyn Fn(usize) -> usize>,
    /// Whether each Bioseq is no nucleotide sequence (`bioseq_protein`), by subject index.
    pub(super) protein: Vec<bool>,
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
pub(super) fn query_sets_taxon(queries: &Queries) -> bool {
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
    pub(super) fn of(
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
                SubjectSource::Query(at) => {
                    (scope.queries[at].bioseq, database || believe, sets_taxon)
                }
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
pub(super) fn set_titles(queries: &Queries) -> bool {
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
            "a translated search or a program other than blastn, blastp, psiblast, deltablast, \
             PHI-BLAST of a protein query and rpsblast",
        ));
    }
    if matches!(request.queries, Queries::SeqLocList(_)) {
        return Err(refused(
            "queries given as locations, whose sequences NCBI BLAST+ reads from the database",
        ));
    }
    if set_titles(&request.queries) {
        return Err(refused("a title on the Bioseq-set of the queries"));
    }
    let spec = TabularSpec::of_formatter(inputs.choice).expect("a tabular format");
    let commented = spec.kind == TabularKind::Commented;
    if inputs.database.is_none() && commented && std::env::var_os("BL2SEQ_LEGACY").is_some() {
        return Err(refused(
            "-subject sequences, BL2SEQ_LEGACY set and outfmt 7, whose lines NCBI BLAST+ writes \
             for each subject",
        ));
    }
    let program = super::search_program::program_type(options.program);
    let blastn = program.blastn;
    // `m_IsUngappedSearch`: the last `GappedMode` call of the options (gapped by default).
    let ungapped = options
        .setters
        .iter()
        .rev()
        .find_map(|setter| match setter {
            Setter::GappedMode(gapped) => Some(!gapped),
            _ => None,
        })
        .unwrap_or(false);
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

    let queries: Vec<QueryEntry<'_>> = request
        .queries
        .bioseqs()
        .into_iter()
        .map(|bioseq| QueryEntry::new(bioseq, inputs.believe))
        .collect();
    let subject_bioseqs = request.subject.bioseqs();
    if sets.len() > queries.len() {
        return Err(refused("a result set without its query Bioseq"));
    }
    let scope = ScopeBioseqs {
        query_ids: IdIndex::of_bioseqs(queries.iter().map(|query| query.bioseq)),
        queries: &queries,
        subject_ids: IdIndex::of_bioseqs(subject_bioseqs.iter().copied()),
        subjects: &subject_bioseqs,
    };
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
        let mut out = Vec::new();
        writer.round_query_header(
            &mut out,
            iteration,
            &queries[q].title,
            database_text,
            Some(plan.count),
        )?;
        let mut rows = TabularRows::default();
        for (n, &(q_idx, s_idx)) in plan.rows.iter().enumerate() {
            let (q_idx, s_idx) = (q_idx as usize, s_idx as usize);
            let (qcovs, qcovus) = plan.scores.get(n).copied().unwrap_or((-1, -1));
            let align = &plan.aligns[n];
            let query = &queries[q_idx];
            let ctx = RowContext {
                spec: &spec,
                no_fetch: blastn,
                query_length: query.length,
            };
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
            let subject = entry
                .values
                .subject(names.as_ref(), entry.nucleotide, qcovs, qcovus);
            let kept = rows.last_found().map(|kept| {
                let kept = entry_of(&entries, kept);
                kept.values
                    .subject(names.as_ref(), kept.nucleotide, qcovs, qcovus)
            });
            // The texts of the row's lookups of the taxonomy database follow the rows before it.
            post_row_texts(sink, &mut out, names.as_ref(), &spec, &writer, &subject)?;
            writer.row(
                &mut out,
                &mut rows,
                &query.query(),
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
/// Dense-seg's or Dense-diags', as the tabular formats read them, `ReadAlign::dense_seg`).
/// `Ok(false)` where one of NCBI's exceptions ends them (the prefetch's warning); a stop that
/// NCBI has no result for, and a form that LOSAT does not read as NCBI, is refused.
pub(super) fn prefetch_reads(
    sets: &[std::ops::Range<usize>],
    alignments: &[V<'_>],
) -> Result<bool> {
    for set in sets {
        for align in read_set(&alignments[set.clone()], false) {
            let read = seq_align_id(align.align, 1)
                .map_err(ReadStop::Stop)
                .and_then(|_| match align.variant() {
                    "denseg" | "dendiag" => align.dense_seg().map(|_| ()),
                    "std" | "disc" => Err(ReadStop::Unsupported(
                        "Std-seg or discontinuous Seq-aligns (those of a translated search)",
                    )),
                    _ => Err(ReadStop::Unsupported(
                        "segments other than a Dense-seg or Dense-diags",
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
    /// The Bioseqs whose letters read (`LetterCache::sequence` keys).
    checked: HashSet<(u8, usize)>,
}

impl<'p> PrePass<'p, '_, '_> {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:749-1096
    /// The reads of the row of `align` in `SetFields`' order: its scores, its query (`GetSeq_id`
    /// of row 0, among the request's queries), its subject (row 1, `Subjects::index_of`), the
    /// coverage and Dense-seg that it reads (`row_plan`), and where it reads letters
    /// (`reads_letters`) those of its query and of a subject that is a Bioseq of the request
    /// (each Bioseq once; a database subject's at the row). The row's query and subject index.
    fn row(
        &mut self,
        align: &ReadAlign<'_>,
        subjects: &mut Subjects<'p>,
    ) -> Result<(u32, u32), RowEnd> {
        if matches!(align.variant(), "std" | "disc") {
            return Err(RowEnd::Refusal(refused(
                "Std-seg or discontinuous Seq-aligns (those of a translated search)",
            )));
        }
        let scores = row_scores(align, self.spec).map_err(end_of)?;
        let query_id = align
            .seq_id(0)
            .map_err(|stop| end_of(ReadStop::Stop(stop)))?;
        let q_idx = match self.scope.queries_holding(query_id).as_slice() {
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
        let ctx = RowContext {
            spec: self.spec,
            no_fetch: self.blastn,
            query_length: self.scope.queries[q_idx].length,
        };
        row_plan(align, scores, &ctx).map_err(end_of)?;
        if reads_letters(self.spec, self.blastn) {
            let query = (2, q_idx, self.scope.queries[q_idx].bioseq, "a query Bioseq");
            let subject = match subjects.sources[s_idx] {
                SubjectSource::Database { .. } => None,
                SubjectSource::Bioseq(at) => {
                    Some((1, at, self.scope.subjects[at], "a -subject Bioseq"))
                }
                SubjectSource::Query(at) => {
                    Some((2, at, self.scope.queries[at].bioseq, "a query Bioseq"))
                }
            };
            for (kind, at, bioseq, what) in std::iter::once(query).chain(subject) {
                if self.checked.contains(&(kind, at)) {
                    continue;
                }
                match read_bioseq(bioseq, what) {
                    Ok(_) => {
                        self.checked.insert((kind, at));
                    }
                    Err(LetterEnd::Exception(exception)) => {
                        return Err(RowEnd::Ncbi(NcbiStop {
                            error: crate::cli::NativeError::exception(&exception).into(),
                            subject: Some(s_idx),
                        }))
                    }
                    Err(LetterEnd::Refusal(error)) => return Err(RowEnd::Refusal(error)),
                }
            }
        }
        Ok((q_idx as u32, s_idx as u32))
    }
}

/// How the read of a row's letters ends the run: NCBI's exception (`CSeqVector_CI`'s, as
/// `report::tabular::query_sequence_end`), or LOSAT's refusal.
pub(super) enum LetterEnd {
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
pub(super) fn read_bioseq(
    bioseq: Bioseq<'_>,
    what: &str,
) -> std::result::Result<Vec<u8>, LetterEnd> {
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
pub(super) struct LetterCache {
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
    fn row(
        &mut self,
        plan: &RowPlan,
        ctx: &RowContext<'_>,
        scope: &ScopeBioseqs<'_, '_>,
        inputs: &ObjectInputs<'_>,
        source: &SubjectSource,
        q_idx: usize,
    ) -> std::result::Result<RowLetters<'_>, LetterEnd> {
        self.letters(plan.subject_range(ctx), scope, inputs, source, q_idx)
    }

    /// The letters of a row of the query `q_idx` and the subject `source` that reads the
    /// subject's letters `range` (from the first to the last, plus strand; none where `None`).
    pub(super) fn letters(
        &mut self,
        range: Option<(u32, u32)>,
        scope: &ScopeBioseqs<'_, '_>,
        inputs: &ObjectInputs<'_>,
        source: &SubjectSource,
        q_idx: usize,
    ) -> std::result::Result<RowLetters<'_>, LetterEnd> {
        let Some((from, to)) = range else {
            return Ok(RowLetters {
                query: &[],
                subject: &[],
                subject_from: 0,
            });
        };
        if self.query.as_ref().is_none_or(|(at, _)| *at != q_idx) {
            let letters = read_bioseq(scope.queries[q_idx].bioseq, "a query Bioseq")?;
            self.query = Some((q_idx, letters));
        }
        let db_failed = |error| LetterEnd::Refusal(database_error(error, PROGRAM));
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
            SubjectSource::Query(at) => Some(((2, at), Some(scope.queries[at].bioseq))),
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
        })
    }
}

#[cfg(test)]
#[path = "tabular_tests.rs"]
mod tests;
