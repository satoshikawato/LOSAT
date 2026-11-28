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

use std::collections::{BTreeMap, HashMap};
use std::sync::Arc;

use anyhow::Result;

use super::formatter_options::FormatterOptions;
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
use crate::report::archive::read::options::{BuiltOptions, EProgram};
use crate::report::archive::read::queries::seq_vector::{RowVector, SeqMap};
use crate::report::archive::read::queries::RenderError;
use crate::report::archive::read::request::{Bioseq, Queries, Request};
use crate::report::archive::read::{seq_align_id, Stop};
use crate::report::db_subject::{HitLoad, LoadsBefore, ParsedSubject, ScopeWalk};
use crate::report::serial::asn_types::view::V;
use crate::report::serial::seqalign::Strand;
use crate::report::serial::seqalign_read::{read_set, PairDenseg, ReadAlign, ReadStop};
use crate::report::tabular::{
    failed_database_subject, failed_record_title, row_stop_refusal, taxdb_warning_line_of,
    MapLengthLookup, SubjectLookup, TabularField, TabularKind, TabularNames, TabularQuery,
    TabularRows, TabularSpec, TabularTaxonomy, TabularWriter,
};
use crate::report::tabular_seqalign::{
    caught_hsp, reads_letters, row_hsp, row_plan, row_reads, row_scores, subject_scores,
    RowBuffers, RowContext, RowLetters, RowPlan, RowSequence, RowSource, SubjectScores,
};
use crate::seqdb::lookup::SeqidLookup;
use crate::seqdb::{NucleotideAmbiguity, SeqDb, SeqType};
use crate::seqid::defline_set::{compare_ordered, tabular_id, TabularId};
use crate::seqid::id_scope::{scope_finds, ScopeIndex};
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
    /// The values of the options that the reports print (`FormatterOptions`, the run's: an
    /// `rpsblast` object's `GetResultSet` has set its matrix name).
    pub(super) formatter: &'a FormatterOptions,
    /// The options' matrix name when the object's formatter is made (`opts.GetMatrixName()`,
    /// before `GetResultSet`), `None` for NULL.
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
    refused_in(super::refusals::TABULAR_FORMATS, what)
}

/// The refusal of a form that the formats `format` (`refusals::TABULAR_FORMATS`,
/// `refusals::PAIRWISE_FORMATS`) do not write yet.
pub(super) fn refused_in(format: &'static str, what: &str) -> anyhow::Error {
    Refusal::NotCovered {
        format,
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

/// The scope index (`id_scope::ScopeIndex`) of the ids of `bioseqs`, by position.
pub(super) fn scope_index<'r>(bioseqs: impl IntoIterator<Item = Bioseq<'r>>) -> ScopeIndex {
    let mut index = ScopeIndex::default();
    for (at, bioseq) in bioseqs.into_iter().enumerate() {
        for id in bioseq.ids() {
            index.add(at, id);
        }
    }
    index
}

/// The request's Bioseqs as NCBI's scope finds them by id (`id_scope::ScopeIndex::finds`: a
/// Bioseq of the id itself, else one that the id matches weakly; a versioned text Seq-id finds no
/// Bioseq of that accession without its version): the queries' scope, which blast_formatter adds
/// before the database's loader and the `-subject` Bioseqs' scope (at
/// `kSubjectsDataLoaderPriority`), so that a query Bioseq is found first.
pub(super) struct ScopeBioseqs<'b, 'r> {
    pub(super) queries: &'b [QueryEntry<'r>],
    pub(super) query_ids: ScopeIndex,
    pub(super) subjects: &'b [Bioseq<'r>],
    pub(super) subject_ids: ScopeIndex,
}

impl ScopeBioseqs<'_, '_> {
    /// The queries that NCBI's scope finds for `id`, in order.
    pub(super) fn queries_holding(&self, id: &ParsedSeqId) -> Vec<usize> {
        self.query_ids.finds(id)
    }

    /// The `-subject` Bioseqs that NCBI's scope finds for `id`, in order.
    fn subjects_holding(&self, id: &ParsedSeqId) -> Vec<usize> {
        self.subject_ids.finds(id)
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
    /// The position in `keys` of each id. Ids are one key where they are equal: two ids that
    /// NCBI's scope may find as two sequences (`gb|X|` and `gb|X.1|`) are never one.
    key_at: BTreeMap<ParsedSeqId, usize>,
    /// The subject index of each sequence (`SubjectSource::sequence`).
    by_sequence: HashMap<(u8, usize), usize>,
    /// Index in `keys` of the most recently matched subject id (consecutive HSPs of one subject).
    last: Option<usize>,
    /// One `SeqidLookup`, `HeaderReader` and `ScopeWalk` opened lazily for the whole object.
    db_readers: Option<(SeqidLookup<'db>, HeaderReader<'db>, ScopeWalk)>,
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
            if self.keys[at].0 == *id {
                return Ok(self.keys[at].1);
            }
        }
        if let Some(&at) = self.key_at.get(id) {
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
        self.key_at.insert(id.clone(), key_at);
        self.keys.push((id.clone(), index));
        self.last = Some(key_at);
        Ok(index)
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:393
/// The database sequence that NCBI's BLAST database data loader serves for the Seq-align
/// subject id `id`: the first OID of the restricted database's lookup (`seqid_to_oid`), its
/// header read as the search programs read a reported sequence's (`read_header`). Not found is
/// NCBI's ask of its servers (`Refusal::Database`); a lookup that throws, a header with a
/// UTF8String, and a parsed header whose Seq-align id NCBI's lookup or scope may not give as it is
/// (`ParsedSubject::found`, `unfetchable`, `load_throws`, `scope_finds_none`, `ScopeWalk::load`,
/// or a Seq-align id that NCBI's scope does not find for `id`, `id_scope::scope_finds`) are
/// refused.
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
            && scope_finds(subject.alignment_id(), id);
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
                (1, [ParsedSeqId::Local(local)]) => {
                    let text = match local {
                        ObjectId::Str(text) => String::from_utf8(text.clone())
                            .map_err(|_| refused("a subject Bioseq whose local id is not UTF-8"))?,
                        ObjectId::Id(n) => n.to_string(),
                    };
                    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:522-532
                    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:558-565
                    // `s_ReplaceLocalId` keeps the local id's text where its flag is set: for the
                    // subject id unless `m_ParseSubjectDefline` (`parse_local`: always in a
                    // database object; a -subject object is bl2seq's, whose subject deflines are
                    // parsed unless `-parse_deflines`), for the `sall` ids where
                    // `m_ParseLocalIds` (`-parse_deflines`): with the option the record shows its
                    // local id as a record of `-parse_deflines` does (`believing`).
                    let id = if parse_local {
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
/// read): the prolog, the prefetch, each query's report, the epilog. Every refusal of the
/// archive by LOSAT comes before the object's first byte but that of a row where NCBI has no
/// result (`TabularWriter::row`'s stop), which follows NCBI's bytes before the stop as in the
/// search programs; NCBI's exceptions come where NCBI meets them, after the bytes before (the
/// pre-pass finds them first). An error of a read of the database (NCBI reads it at the row
/// too) or of a write of the output ends the run where it happens, after the queries written
/// before.
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
    let legacy = inputs.database.is_none() && std::env::var_os("BL2SEQ_LEGACY").is_some();
    if legacy && commented {
        return Err(refused(
            "-subject sequences, BL2SEQ_LEGACY set and outfmt 7, whose lines NCBI BLAST+ writes \
             for each subject",
        ));
    }
    let program = super::search_program::program_type(options.program);
    let blastn = program.blastn;
    // `m_IsUngappedSearch`: the last `GappedMode` call of the options (gapped by default).
    let ungapped = !inputs.formatter.gapped;
    let matrix = formatter_matrix(inputs.matrix_name);
    let version = format!(
        "{} {}",
        super::search_program::formatter_program(options.program).to_ascii_uppercase(),
        crate::wording::report_version!()
    );
    // The program whose rules decide where a row with an empty Bioseq handle stops
    // (`TabularWriter::row`): the formatter reads a nucleotide or a protein search's rows.
    let row_program = if blastn {
        crate::blastinput::app::Program::Blastn
    } else {
        crate::blastinput::app::Program::Blastp
    };
    let writer = TabularWriter::new(&spec, row_program, &version, matrix);
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
        query_ids: scope_index(queries.iter().map(|query| query.bioseq)),
        queries: &queries,
        subject_ids: scope_index(subject_bioseqs.iter().copied()),
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
    let prefetched =
        prefetches && prefetch_reads(sets, alignments, super::refusals::TABULAR_FORMATS)?;
    if prefetched {
        for set in sets {
            for &align in &alignments[set.clone()] {
                let id = seq_align_id(align, 1).map_err(super::archive::stop_error)?;
                subjects.index_of(inputs, &scope, id)?;
            }
        }
    }
    let prefetched_subjects = subjects.sources.len();
    // The sequence maps that the prefetch makes, and the subjects whose length then throws in
    // the rows' subject block (`MadeMaps::length_throws`).
    let mut made = MadeMaps::new(queries.len(), subject_bioseqs.len());
    if prefetched {
        prefetch_maps(sets, alignments, &mut subjects, inputs, &scope, &mut made)?;
    }
    let length_throws = made.length_throws(&subjects.sources, &scope);

    // The queries' reports as NCBI makes them, read before any byte (every refusal comes
    // first): each query's Seq-aligns (the ungapped split where the result set has
    // alignments), `PruneSeqalign`, `InsertSubjectScores` and each row's reads, up to the first
    // of NCBI's exceptions, which ends the object after the rows before it. The round that
    // `PrintOneResultSet` takes is read for the first query (`archive::psi_round`): a
    // `Web_StepNumber` that is no integer throws there (`GetInteger`, a serial exception of the
    // loop; probe E1 probe1 O8).
    let mut checks = PrePass {
        inputs,
        scope: &scope,
        spec: &spec,
        blastn,
        maps: BioseqMaps::default(),
        buffers: RowBuffers::default(),
        made,
        length_throws: &length_throws,
    };
    let mut plans: Vec<QueryPlan<'_>> = Vec::with_capacity(sets.len());
    let mut stop: Option<NcbiStop> = None;
    let round = super::archive::psi_round(inputs.iteration, sets.len());
    if round.is_ok() {
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
        let mut runs: Vec<usize> = plan.rows.iter().map(|row| row.subject as usize).collect();
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
    let iteration = round?;

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
        for (n, row) in plan.rows.iter().enumerate() {
            let (q_idx, s_idx) = (row.query as usize, row.subject as usize);
            let (qcovs, qcovus) = plan.scores.get(n).copied().unwrap_or((-1, -1));
            let align = &plan.aligns[n];
            let query = &queries[q_idx];
            let ctx = RowContext {
                spec: &spec,
                no_fetch: blastn,
                query_length: row.query_length.map_or(query.length, i64::from),
            };
            let planned = row_scores(align, &spec)
                .and_then(|scores| row_plan(align, scores, &ctx))
                .map_err(read_stop)?;
            let hsp = if row.length_throws {
                caught_hsp(&planned)
            } else {
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
                row_hsp(&planned, &ctx, &sequences, &mut buffers).map_err(read_stop)?
            };
            let entry = entry_of(&entries, s_idx);
            let mut subject = entry
                .values
                .subject(names.as_ref(), entry.nucleotide, qcovs, qcovus);
            if let Some(length) = row.subject_length {
                subject.length = length as usize;
            }
            // The catch pushes the row's Seq-align id (`MapLengthLookup`).
            let pushed: [Vec<u8>; 2];
            if row.length_throws {
                let id = align
                    .seq_id(1)
                    .map_err(|stop| read_stop(ReadStop::Stop(stop)))?;
                pushed = [TabularId::Full, TabularId::Accession].map(|kind| {
                    let mut out = Vec::new();
                    tabular_id(std::slice::from_ref(id), kind, &mut out);
                    out
                });
                subject.lookup = SubjectLookup::MapLength(MapLengthLookup {
                    full: &pushed[0],
                    gi: id.gi().unwrap_or(0),
                    accession: &pushed[1],
                    subject_id: matches!(subject.lookup, SubjectLookup::Failed(f) if f.subject_id),
                });
            }
            let mut query_row = query.query();
            if let Some(length) = row.query_length {
                query_row.length = length as usize;
            }
            let kept = rows.last_found().map(|kept| {
                let kept = entry_of(&entries, kept);
                kept.values
                    .subject(names.as_ref(), kept.nucleotide, qcovs, qcovus)
            });
            // The texts of the row's lookups of the taxonomy database follow the rows before it.
            post_row_texts(sink, &mut out, names.as_ref(), &spec, &writer, &subject)?;
            // A row where NCBI has no result ends the run after the bytes before its stop
            // (`report::tabular::row_stop_end`).
            if let Some(stop) = writer.row(
                &mut out,
                &mut rows,
                &query_row,
                &subject,
                kept.as_ref(),
                &hsp,
            ) {
                sink.write(&out)?;
                sink.flush()?;
                let describe = match &subjects.sources[stop.subject] {
                    SubjectSource::Database { oid, .. } => {
                        failed_database_subject(oid, None, None, stop.lineless)
                    }
                    SubjectSource::Bioseq(at) | SubjectSource::Query(at) => {
                        failed_record_title(*at, records.records[stop.subject].shown_id())
                    }
                };
                return Err(row_stop_refusal(stop, describe, PROGRAM));
            }
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
/// NCBI has no result for, and a form that LOSAT does not read as NCBI, is refused (in the formats
/// of the run, `format`).
pub(super) fn prefetch_reads(
    sets: &[std::ops::Range<usize>],
    alignments: &[V<'_>],
    format: &'static str,
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
            // A form LOSAT does not read is refused in the formats of the run (`format`).
            let read = read.map_err(|stop| match stop {
                ReadStop::Unsupported(what) => End::Refusal(refused_in(format, what)),
                stop => end_of(stop),
            });
            match read {
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
/// reads them (outfmt 7's line), those `PruneSeqalign` keeps, their coverages, and each row that
/// NCBI prints (to its exception).
struct QueryPlan<'a> {
    count: usize,
    aligns: Vec<ReadAlign<'a>>,
    scores: SubjectScores,
    rows: Vec<PlannedRow>,
}

/// A row that NCBI prints: its query and subject index, and what its `SetFields` reads of their
/// lengths where a sequence map makes them other than declared (`MadeMaps::map_length`).
#[derive(Clone, Copy, Debug)]
struct PlannedRow {
    query: u32,
    subject: u32,
    /// `m_QueryLength` (`qlen`, `qcovhsp`) and `m_SubjectLength` (`slen`) where a sequence map
    /// gives them; `None`: the declared lengths.
    query_length: Option<u32>,
    subject_length: Option<u32>,
    /// The subject block's catch: `GetBioseqLength` threw (`MapLengthLookup`), the row has the
    /// values of `x_ResetFields` but the scores, the query's and the subject's ids.
    length_throws: bool,
}

/// The reads of the rows before any byte (`PrePass::row`).
struct PrePass<'p, 'b, 'r> {
    inputs: &'p ObjectInputs<'p>,
    scope: &'p ScopeBioseqs<'b, 'r>,
    spec: &'p TabularSpec,
    blastn: bool,
    /// The sequence maps of the current query and subject Bioseq.
    maps: BioseqMaps<'r>,
    buffers: RowBuffers,
    /// Which Bioseqs have their sequence map made in NCBI's scope.
    made: MadeMaps,
    /// By subject index: whether its length throws in the subject block.
    length_throws: &'p [bool],
}

impl<'p> PrePass<'p, '_, '_> {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:749-1096
    /// The reads of the row of `align` in `SetFields`' order: its scores, its query (`GetSeq_id`
    /// of row 0, among the request's queries), its subject (row 1, `Subjects::index_of`), the
    /// coverage and Dense-seg that it reads (`row_plan`), and where it reads letters
    /// (`reads_letters`) the reads of its query and of a subject that is a Bioseq of the request
    /// (`row_reads`: each row's vector, `RowVector`, whose throw is NCBI's at the row; a
    /// database subject's letters are read at the row). The row's query and subject index.
    fn row(
        &mut self,
        align: &ReadAlign<'_>,
        subjects: &mut Subjects<'p>,
    ) -> Result<PlannedRow, RowEnd> {
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
        // `GetBioseqLength` of the query and of the subject where a sequence map gives it: the
        // query block's (a catch there is LOSAT's refusal), the subject block's (its catch ends
        // the row's fields), then `qlen`'s and `slen`'s own reads where those blocks did not set
        // them, outside any `try` (NCBI's exception at the row).
        let spec = self.spec;
        let query_key = BioseqKey::Query(q_idx);
        let subject_key = BioseqKey::of_subject(&subjects.sources[s_idx]);
        let query_block = spec.has_any(&[
            TabularField::QuerySeqId,
            TabularField::QueryGi,
            TabularField::QueryAccession,
            TabularField::QueryAccessionVersion,
            TabularField::QueryCovHsp,
        ]);
        let subject_block = spec.reads_subject_bioseq();
        let made = &self.made;
        let scope = self.scope;
        let read_length = |key: BioseqKey| match made.map_length(key, key.bioseq(scope)) {
            None => Ok(None),
            Some(Ok(length)) => Ok(Some(length)),
            Some(Err(error)) => Err(length_end(error, key, s_idx)),
        };
        if query_block
            && matches!(
                made.map_length(query_key, query_key.bioseq(scope)),
                Some(Err(_))
            )
        {
            return Err(RowEnd::Refusal(refused(
                "a query Bioseq declaring 4294967295 letters whose sequence map (made by the \
                 prefetch of a self-hit) throws at its length in a row's query block, where \
                 NCBI BLAST+'s catch shows the Seq-align's query id",
            )));
        }
        let length_throws = subject_block && self.length_throws.get(s_idx) == Some(&true);
        let mut row = PlannedRow {
            query: q_idx as u32,
            subject: s_idx as u32,
            query_length: None,
            subject_length: None,
            length_throws,
        };
        if query_block {
            row.query_length = read_length(query_key)?;
        }
        let ctx = RowContext {
            spec,
            no_fetch: self.blastn,
            query_length: row
                .query_length
                .map_or(self.scope.queries[q_idx].length, i64::from),
        };
        let plan = row_plan(align, scores, &ctx).map_err(end_of)?;
        if length_throws {
            return Ok(row);
        }
        if subject_block {
            if let Some(key) = subject_key {
                row.subject_length = read_length(key)?;
            }
        }
        let declared_zero = |key: BioseqKey| key.bioseq(self.scope).length() == Some(0);
        if spec.has(TabularField::QueryLength)
            && (!query_block
                || row
                    .query_length
                    .map_or(declared_zero(query_key), |l| l == 0))
        {
            row.query_length = read_length(query_key)?;
        }
        if let Some(key) = subject_key {
            if spec.has(TabularField::SubjectLength)
                && (!subject_block || row.subject_length.map_or(declared_zero(key), |l| l == 0))
            {
                row.subject_length = read_length(key)?;
            }
        }
        if reads_letters(spec, self.blastn) {
            let keys = [Some(query_key), subject_key];
            self.maps.set(self.scope, keys);
            let maps = &self.maps;
            let mut vectors: [Option<RowVector<'_, '_>>; 2] = [None, None];
            let read = row_reads(&plan, &ctx, &mut self.buffers, |row, plus, from, to| {
                let Some(key) = keys[row] else {
                    return Ok(());
                };
                let vector = match &mut vectors[row] {
                    Some(vector) => vector,
                    vector @ None => match maps.map(row) {
                        Some(map) => vector.insert(RowVector::new(map, plus)),
                        None => return Err(LetterError::Read(maps.error(row, self.scope), key)),
                    },
                };
                vector
                    .read(from, to)
                    .map_err(|error| LetterError::Read(error, key))
            });
            // `CAlnVec`'s vector of a row makes its Bioseq's sequence map.
            for (vector, key) in vectors.iter().zip(keys) {
                if let (Some(_), Some(key)) = (vector, key) {
                    self.made.make(key);
                }
            }
            match read {
                Ok(()) => {}
                Err(LetterError::Stop(stop)) => return Err(end_of(stop).into()),
                Err(LetterError::Read(error, key)) => return Err(length_end(error, key, s_idx)),
            }
        }
        Ok(row)
    }
}

/// How a read of the Bioseq `key` (its letters or its length) ends the object at a row of the
/// subject `s_idx`: NCBI's exception there, after the subject's values (`NcbiStop`), or LOSAT's
/// refusal.
fn length_end(error: RenderError, key: BioseqKey, s_idx: usize) -> RowEnd {
    match letter_end(error, key.what()) {
        LetterEnd::Exception(exception) => RowEnd::Ncbi(NcbiStop {
            error: crate::cli::NativeError::exception(&exception).into(),
            subject: Some(s_idx),
        }),
        LetterEnd::Refusal(error) => RowEnd::Refusal(error),
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

/// A Bioseq of the request that a row reads: a query's or a `-subject` Bioseq.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BioseqKey {
    Query(usize),
    Subject(usize),
}

impl BioseqKey {
    /// The Bioseq of the subject `source`, none for a database sequence.
    fn of_subject(source: &SubjectSource) -> Option<Self> {
        match *source {
            SubjectSource::Database { .. } => None,
            SubjectSource::Bioseq(at) => Some(BioseqKey::Subject(at)),
            SubjectSource::Query(at) => Some(BioseqKey::Query(at)),
        }
    }

    fn bioseq<'r>(self, scope: &ScopeBioseqs<'_, 'r>) -> Bioseq<'r> {
        match self {
            BioseqKey::Query(at) => scope.queries[at].bioseq,
            BioseqKey::Subject(at) => scope.subjects[at],
        }
    }

    /// The Bioseq's name in a refusal.
    fn what(self) -> &'static str {
        match self {
            BioseqKey::Query(_) => "a query Bioseq",
            BioseqKey::Subject(_) => "a -subject Bioseq",
        }
    }
}

/// The sequence maps of the Bioseqs of a row's query and subject (`SeqMap::of`, one Bioseq at a
/// time for each), or none where the map's constructor ends the run (at the row's first read of
/// it, `BioseqMaps::error`).
#[derive(Default)]
struct BioseqMaps<'r> {
    rows: [Option<(BioseqKey, Option<SeqMap<'r>>)>; 2],
}

impl<'r> BioseqMaps<'r> {
    /// The maps of the row's query and subject Bioseq (`keys`; none for a database subject), each
    /// made where the current one is another Bioseq's.
    fn set(&mut self, scope: &ScopeBioseqs<'_, 'r>, keys: [Option<BioseqKey>; 2]) {
        for (slot, key) in self.rows.iter_mut().zip(keys) {
            let Some(key) = key else {
                continue;
            };
            if slot.as_ref().is_none_or(|(at, _)| *at != key) {
                *slot = Some((key, SeqMap::of(key.bioseq(scope)).ok()));
            }
        }
    }

    /// The map of row `row` (0 the query, 1 the subject) after `set`.
    fn map(&self, row: usize) -> Option<&SeqMap<'r>> {
        self.rows[row].as_ref().and_then(|(_, map)| map.as_ref())
    }

    /// The end of the constructor of the map of row `row` (`map` is none).
    fn error(&self, row: usize, scope: &ScopeBioseqs<'_, 'r>) -> RenderError {
        let (key, _) = self.rows[row].as_ref().expect("a map set");
        match SeqMap::of(key.bioseq(scope)) {
            Err(error) => error,
            Ok(_) => unreachable!("a map whose constructor ended the run"),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objmgr/bioseq_info.cpp:676-685, 704-712
/// Which Bioseqs of the request (`BioseqKey`) have their sequence map made in NCBI's scope
/// (`CBioseq_Info::m_SeqMap`): the prefetch makes some (`prefetch_maps`), the vector of each row's
/// `CAlnVec` makes its Bioseq's. `GetBioseqLength` is the declared length, but where the map is
/// made the map's `GetLength` (`map_length`), which differs only for a Bioseq that declares
/// 4294967295 letters (`SeqMap::length`: it resolves the segments, and throws at one of
/// 4294967295 letters; probes L1 of agent E3L and t_dlen, t_huge* of agent E3, 2026-11-27).
struct MadeMaps {
    queries: Vec<bool>,
    subjects: Vec<bool>,
}

impl MadeMaps {
    /// No map made, of `queries` query and `subjects` -subject Bioseqs.
    fn new(queries: usize, subjects: usize) -> Self {
        MadeMaps {
            queries: vec![false; queries],
            subjects: vec![false; subjects],
        }
    }

    fn slot(&mut self, key: BioseqKey) -> &mut bool {
        match key {
            BioseqKey::Query(at) => &mut self.queries[at],
            BioseqKey::Subject(at) => &mut self.subjects[at],
        }
    }

    /// The map of `key` made.
    fn make(&mut self, key: BioseqKey) {
        *self.slot(key) = true;
    }

    /// `GetBioseqLength` of the Bioseq `key` (`bioseq`) where its made map gives it (a Bioseq
    /// that declares 4294967295 letters): the segments' end, or the map's throw; `None`: the
    /// declared length.
    fn map_length(
        &self,
        key: BioseqKey,
        bioseq: Bioseq<'_>,
    ) -> Option<std::result::Result<u32, RenderError>> {
        let made = match key {
            BioseqKey::Query(at) => self.queries[at],
            BioseqKey::Subject(at) => self.subjects[at],
        };
        (made && bioseq.length() == Some(i64::from(u32::MAX)))
            .then(|| SeqMap::of(bioseq).and_then(|map| map.length()))
    }

    /// By subject index of `sources`: whether a row's subject block throws at the subject's
    /// length (`GetBioseqLength` of a map that the prefetch made: a map made by a row's read
    /// throws at that read first, which ends the object).
    fn length_throws(&self, sources: &[SubjectSource], scope: &ScopeBioseqs<'_, '_>) -> Vec<bool> {
        sources
            .iter()
            .map(|source| {
                BioseqKey::of_subject(source).is_some_and(|key| {
                    matches!(self.map_length(key, key.bioseq(scope)), Some(Err(_)))
                })
            })
            .collect()
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:658-690
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1110-1153
// NCBI reference (598d8ae6): c++/src/objmgr/seq_map.cpp:986-1117
/// The sequence maps that the prefetch makes (`LoadSequencesToScope`, where it reads every
/// Seq-align, `prefetch_reads`): a virtual delta of one interval per Seq-align
/// (`s_ExtractSeqidsAndRanges`: the subject's id and `GetSeqRange` of its row, by id in
/// `CompareOrdered` order, each id's in the order of the Seq-aligns; the range of a minus query
/// row over a plus subject row taken back from the subject's length, `length - to_open` to
/// `length - from + 1`), over which `CanResolveRange` makes the map of each subject that it
/// reaches, interval by interval, and stops at the first throw (its catch, no warning): a map
/// whose constructor throws, or an iteration of the interval that resolves a segment of
/// 4294967295 letters (`SeqMap::plus_range_throws`). The maps change no length but of a Bioseq
/// that declares 4294967295 letters (`MadeMaps::map_length`); without one nothing is read. Where
/// one is, ranges that NCBI does not read as these intervals, and Seq-ids that NCBI orders by
/// their addresses, are refused.
fn prefetch_maps<'db, 'r>(
    sets: &[std::ops::Range<usize>],
    alignments: &[V<'_>],
    subjects: &mut Subjects<'db>,
    inputs: &ObjectInputs<'db>,
    scope: &ScopeBioseqs<'_, 'r>,
    made: &mut MadeMaps,
) -> Result<()> {
    let declares_huge = |bioseq: Bioseq<'_>| bioseq.length() == Some(i64::from(u32::MAX));
    if !scope
        .queries
        .iter()
        .any(|query| declares_huge(query.bioseq))
        && !scope.subjects.iter().copied().any(declares_huge)
    {
        return Ok(());
    }
    let unread = |what: &str| {
        refused(&format!(
            "a Bioseq declaring 4294967295 letters with {what}, which NCBI BLAST+'s prefetch \
             reads into its sequence maps"
        ))
    };
    let mut entries: Vec<(&ParsedSeqId, usize, u32, u32, bool)> = Vec::new();
    for set in sets {
        for align in read_set(&alignments[set.clone()], false) {
            let id = seq_align_id(align.align, 1).map_err(super::archive::stop_error)?;
            let ds = align.dense_seg().map_err(read_stop)?;
            let (from, to) = subject_range(&align, &ds).ok_or_else(|| {
                unread("a Seq-align whose subject row's range is not one interval")
            })?;
            let flipped = ds.strands == [Strand::Minus, Strand::Plus];
            let s_idx = subjects.index_of(inputs, scope, id)?;
            entries.push((id, s_idx, from, to, flipped));
        }
    }
    let mut unordered = false;
    entries.sort_by(|a, b| {
        compare_ordered(a.0, b.0).unwrap_or_else(|| {
            unordered = true;
            std::cmp::Ordering::Equal
        })
    });
    if unordered {
        return Err(unread(
            "subject ids that NCBI BLAST+ orders by their addresses",
        ));
    }
    let mut intervals = Vec::with_capacity(entries.len());
    let mut total = 0u64;
    for &(_, s_idx, from, to, flipped) in &entries {
        let key = BioseqKey::of_subject(&subjects.sources[s_idx]);
        let (from, to) = match (flipped, key) {
            (false, _) => (from, to),
            (true, Some(key)) => {
                let length = match made.map_length(key, key.bioseq(scope)) {
                    None => key
                        .bioseq(scope)
                        .length()
                        .and_then(|l| u32::try_from(l).ok()),
                    Some(length) => length.ok(),
                };
                length
                    .and_then(|length| {
                        Some((
                            length.checked_sub(to.checked_add(1)?)?,
                            (length - from).checked_add(1)?,
                        ))
                    })
                    .ok_or_else(|| {
                        unread("a minus query row over a plus subject row past its length")
                    })?
            }
            (true, None) => {
                return Err(unread(
                    "a minus query row over a plus row of a database subject",
                ))
            }
        };
        total += u64::from(to - from) + 1;
        intervals.push((key, from, to));
    }
    if total >= u64::from(u32::MAX) {
        return Err(unread("subject ranges of 2^32 - 1 letters or more"));
    }
    resolve_maps(&intervals, &|key| key.bioseq(scope), made)
}

/// `CanResolveRange` over the prefetch's intervals (`prefetch_maps`: the Bioseq of each, `None`
/// for a database sequence, and its first and last position): the map of each Bioseq reached
/// made, up to the first whose constructor throws or whose iteration of the interval resolves a
/// segment of 4294967295 letters. A map that LOSAT does not read is refused.
fn resolve_maps<'r>(
    intervals: &[(Option<BioseqKey>, u32, u32)],
    bioseq_of: &dyn Fn(BioseqKey) -> Bioseq<'r>,
    made: &mut MadeMaps,
) -> Result<()> {
    let mut current: Option<(BioseqKey, Option<SeqMap<'r>>)> = None;
    for &(key, from, to) in intervals {
        // A database sequence's map resolves its ranges.
        let Some(key) = key else {
            continue;
        };
        if current.as_ref().is_none_or(|(at, _)| *at != key) {
            let map = match SeqMap::of(bioseq_of(key)) {
                Ok(map) => Some(map),
                Err(RenderError::Exception(_)) => None,
                Err(error) => {
                    return Err(match letter_end(error, key.what()) {
                        LetterEnd::Refusal(error) => error,
                        LetterEnd::Exception(_) => unreachable!("an exception kept above"),
                    })
                }
            };
            current = Some((key, map));
        }
        let Some((_, Some(map))) = &current else {
            break;
        };
        made.make(key);
        if map.plus_range_throws(from, to + 1) {
            break;
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:153-169
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:165-237
/// `GetSeqRange` of the subject's row of `align` (Dense-seg `ds`) as `(from, to)`: a Dense-seg's
/// `GetSeqStart` to `GetSeqStop` (on a minus row from the last segment with letters back to the
/// first), the union of the Dense-diags' ranges. `None` where NCBI's range is no such interval: a
/// first start past the last stop, a Dense-diag start of -1, no letters.
fn subject_range(align: &ReadAlign<'_>, ds: &PairDenseg) -> Option<(u32, u32)> {
    let segs = || {
        ds.lens
            .iter()
            .enumerate()
            .map(|(seg, &len)| (ds.starts[2 * seg + 1], len))
    };
    let (from, to) = if align.variant() == "dendiag" {
        let mut range: Option<(u32, u32)> = None;
        for (start, len) in segs() {
            let start = u32::try_from(start).ok()?;
            if len > 0 {
                let stop = start.checked_add(len - 1)?;
                range = Some(range.map_or((start, stop), |(a, b)| (a.min(start), b.max(stop))));
            }
        }
        range?
    } else {
        let lettered = || segs().filter(|&(start, _)| start >= 0);
        let (first, last) = (lettered().next()?, lettered().next_back()?);
        let (at_start, at_stop) = if ds.strands[1] == Strand::Minus {
            (last, first)
        } else {
            (first, last)
        };
        let from = u32::try_from(at_start.0).ok()?;
        let to = u32::try_from(at_stop.0)
            .ok()?
            .checked_add(at_stop.1)?
            .checked_sub(1)?;
        (from, to)
    };
    (from <= to).then_some((from, to))
}

/// How the reads of a row's letters in the pre-pass end (`PrePass::row`): a stop of the row's
/// alignment map, or the end of a read of the Bioseq `BioseqKey`.
enum LetterError {
    Stop(ReadStop),
    Read(RenderError, BioseqKey),
}

impl From<ReadStop> for LetterError {
    fn from(stop: ReadStop) -> Self {
        LetterError::Stop(stop)
    }
}

/// How the read of the letters of a Bioseq that `what` names ends the run.
fn letter_end(error: RenderError, what: &str) -> LetterEnd {
    match error {
        RenderError::Exception(exception) => LetterEnd::Exception(exception),
        RenderError::Unsupported(part) => {
            LetterEnd::Refusal(refused(&format!("{what} with {part}")))
        }
        RenderError::Write(error) => LetterEnd::Refusal(error.into()),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:109-140
/// The sequence map of a Bioseq of the request that a report reads whole (the pairwise report's
/// rows and query-anchored sequences, `pairwise`), after a read of all its letters on the plus
/// strand (`RowVector::read`): the map's or the read's throw ends the run (`LetterEnd`), so that
/// `SeqMap::letters` gives every letter after it. `what` names the Bioseq in a refusal. A Bioseq
/// that declares 4294967295 letters, whose length NCBI takes from its map where a map is made
/// (`MadeMaps::map_length`), is refused.
pub(super) fn read_bioseq_map<'r>(
    bioseq: Bioseq<'r>,
    what: &str,
) -> std::result::Result<SeqMap<'r>, LetterEnd> {
    if bioseq.length() == Some(i64::from(u32::MAX)) {
        return Err(LetterEnd::Refusal(refused(&format!(
            "{what} declaring 4294967295 letters"
        ))));
    }
    let map = SeqMap::of(bioseq).map_err(|error| letter_end(error, what))?;
    let length = map.length().map_err(|error| letter_end(error, what))?;
    if length > 0 {
        RowVector::new(&map, true)
            .read(0, length - 1)
            .map_err(|error| letter_end(error, what))?;
    }
    Ok(map)
}

/// The letters of a Bioseq of the request whole, plus strand, as `CAlnVec`'s vector reads them
/// (`read_bioseq_map`); `what` names it in a refusal.
pub(super) fn read_bioseq(
    bioseq: Bioseq<'_>,
    what: &str,
) -> std::result::Result<Vec<u8>, LetterEnd> {
    let map = read_bioseq_map(bioseq, what)?;
    let mut letters = Vec::new();
    map.letters(0, map.length().unwrap_or(0), &mut letters);
    Ok(letters)
}

/// The letters that the rows read, one sequence at a time: the sequence maps of the current
/// query and subject Bioseq (`BioseqMaps`, whose reads the pre-pass ran: each row's letters are
/// those of its ranges, `SeqMap::letters`), a database protein's letters whole, and of a database
/// nucleotide only the row's range (with the ambiguity runs of the current one), so that the
/// object keeps no more than one sequence of each kind.
#[derive(Default)]
pub(super) struct LetterCache<'r> {
    maps: BioseqMaps<'r>,
    /// The letters of the current database protein (its ordinal id) whole.
    protein: Option<(u32, Vec<u8>)>,
    ambiguity: Option<(u32, NucleotideAmbiguity)>,
    /// A database nucleotide's letters of the row's range, and where they start.
    range: (Vec<u8>, u32),
}

impl<'r> LetterCache<'r> {
    /// The letters of the row whose plan is `plan` (`RowLetters`), of the query `q_idx` and the
    /// subject `source`; none where the row reads none.
    fn row(
        &mut self,
        plan: &RowPlan,
        ctx: &RowContext<'_>,
        scope: &ScopeBioseqs<'_, 'r>,
        inputs: &ObjectInputs<'_>,
        source: &SubjectSource,
        q_idx: usize,
    ) -> std::result::Result<RowLetters<'_>, LetterEnd> {
        self.letters(
            plan.subject_range(ctx),
            plan.rows_plus(),
            scope,
            inputs,
            source,
            q_idx,
        )
    }

    /// The letters of a row of the query `q_idx` and the subject `source` that reads the
    /// subject's letters `range` (from the first to the last, plus strand) with its rows on the
    /// strands `plus` (query, subject: plus where set); none where either is `None`.
    pub(super) fn letters(
        &mut self,
        range: Option<(u32, u32)>,
        plus: Option<[bool; 2]>,
        scope: &ScopeBioseqs<'_, 'r>,
        inputs: &ObjectInputs<'_>,
        source: &SubjectSource,
        q_idx: usize,
    ) -> std::result::Result<RowLetters<'_>, LetterEnd> {
        let (Some((from, to)), Some([query_plus, subject_plus])) = (range, plus) else {
            return Ok(RowLetters {
                query: RowSequence::NONE,
                subject: RowSequence::NONE,
            });
        };
        let keys = [Some(BioseqKey::Query(q_idx)), BioseqKey::of_subject(source)];
        self.maps.set(scope, keys);
        let db_failed = |error| LetterEnd::Refusal(database_error(error, PROGRAM));
        let held = match *source {
            SubjectSource::Database { ordinal, .. } => {
                let (_, db) = inputs.database.expect("a database subject");
                if db.seqtype() == SeqType::Protein {
                    if self.protein.as_ref().is_none_or(|(at, _)| *at != ordinal) {
                        let block = db.read_sequence(ordinal).map_err(db_failed)?;
                        let letters =
                            crate::sequence::stdaa_letters::protein_fasta_letters(block.protein(0));
                        self.protein = Some((ordinal, letters));
                    }
                    let (_, letters) = self.protein.as_ref().expect("read above");
                    Some((&letters[..], 0, letters.len()))
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
                    Some((&self.range.0[..], self.range.1, length))
                }
            }
            SubjectSource::Bioseq(_) | SubjectSource::Query(_) => None,
        };
        let subject = match held {
            Some((letters, offset, length)) => RowSequence {
                letters: RowSource::Held(letters, offset),
                size: u32::try_from(length).unwrap_or(u32::MAX),
                minus: !subject_plus,
            },
            None => map_sequence(self.maps.map(1), subject_plus),
        };
        Ok(RowLetters {
            query: map_sequence(self.maps.map(0), query_plus),
            subject,
        })
    }
}

/// The row's sequence of the Bioseq of `map` (`RowSequence`), on the plus strand where `plus`.
/// A row that reads none of its letters (no segment) leaves its vector unmade: its size is not
/// read.
fn map_sequence<'m>(map: Option<&'m SeqMap<'m>>, plus: bool) -> RowSequence<'m> {
    match map {
        Some(map) => RowSequence {
            letters: RowSource::Map(map),
            size: map.length().unwrap_or(0),
            minus: !plus,
        },
        None => RowSequence::NONE,
    }
}

#[cfg(test)]
#[path = "tabular_tests.rs"]
mod tests;
