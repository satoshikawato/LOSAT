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
//! (`Refusal::Tabular`): the translated programs, queries given as locations, `-subject` with
//! `BL2SEQ_LEGACY` in outfmt 7, and ids that NCBI's scope could find elsewhere.

use std::sync::Arc;

use anyhow::Result;

use super::refusals::Refusal;
use super::PROGRAM;
use crate::algorithm::common::database::{
    database_error, parsed_refusal, prefetch_messages, read_header, subject_taxonomy,
    tabular_names, HeaderReader, ReportedHeader, StoredTitle, TaxdbRun, TaxdbSetting,
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
use crate::report::db_subject::ParsedSubject;
use crate::report::serial::asn_types::view::V;
use crate::report::serial::seqalign_read::{read_set, ReadAlign, ReadStop};
use crate::report::tabular::{
    check_failed_lookups, failed_database_title, failed_record_title, taxdb_warning_line_of,
    LookupFailure, SubjectLookup, TabularKind, TabularNames, TabularQuery, TabularRows,
    TabularSpec, TabularTaxonomy, TabularWriter,
};
use crate::report::tabular_seqalign::{
    hsp, reads_letters, subject_scores, RowBuffers, RowContext, RowLetters,
};
use crate::seqdb::{SeqDb, SeqType};
use crate::seqid::defline_set::{compare, Sic};
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
/// RPS-BLAST); the translated programs come with their Std-segs (PR 3b of line RG item F3).
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
    )
}

/// One query Bioseq of the request as the rows show it.
struct QueryEntry<'r> {
    bioseq: Bioseq<'r>,
    /// `SetQueryId`: the ids after `s_ReplaceLocalId`, whose tabular ids and GI the rows print.
    ids: BelievedIds,
    length: i64,
    /// outfmt 7's `# Query:` text before its spaces are cut (`x_AcknowledgeBlastSequence`).
    title: Vec<u8>,
    letters: Option<Vec<u8>>,
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
    fn new(bioseq: Bioseq<'r>, believe: bool) -> Self {
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
            letters: None,
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
}

/// The object's subjects in the order of their first Seq-align, with the Seq-align subject ids
/// found for each.
#[derive(Default)]
struct Subjects {
    sources: Vec<SubjectSource>,
    keys: Vec<(ParsedSeqId, usize)>,
}

impl Subjects {
    /// The subject index of the Seq-align subject id `id`, found where no Seq-align before had
    /// it.
    fn index_of(
        &mut self,
        inputs: &ObjectInputs<'_>,
        queries: &[QueryEntry<'_>],
        subject_bioseqs: &[Bioseq<'_>],
        id: &ParsedSeqId,
    ) -> Result<usize> {
        if let Some(&(_, index)) = self
            .keys
            .iter()
            .find(|(key, _)| compare(key, id) == Sic::Yes)
        {
            return Ok(index);
        }
        if queries.iter().any(|query| holds(query.bioseq.ids(), id)) {
            return Err(refused(
                "a subject id that names a query Bioseq, which NCBI BLAST+'s scope finds first",
            ));
        }
        let source = match inputs.database {
            None => SubjectSource::Bioseq(only_holder(subject_bioseqs, id).ok_or_else(|| {
                refused("a subject id that names no -subject Bioseq of the request, or several")
            })?),
            Some((name, db)) => database_subject(name, db, id)?,
        };
        // One sequence by two ids: NCBI's data loader gives the sequence it loaded for every id
        // of its own (`RegisterIds`), which LOSAT does not follow.
        let same = |held: &SubjectSource| match (held, &source) {
            (SubjectSource::Database { oid: a, .. }, SubjectSource::Database { oid: b, .. }) => {
                a == b
            }
            (SubjectSource::Bioseq(a), SubjectSource::Bioseq(b)) => a == b,
            _ => false,
        };
        if self.sources.iter().any(same) {
            return Err(refused(
                "Seq-aligns that name one subject sequence by different ids",
            ));
        }
        let index = self.sources.len();
        self.sources.push(source);
        self.keys.push((id.clone(), index));
        Ok(index)
    }
}

/// Whether one of `ids` is `id` for NCBI's scope (`Sic::Yes`).
fn holds<'i>(ids: impl IntoIterator<Item = &'i ParsedSeqId>, id: &ParsedSeqId) -> bool {
    ids.into_iter().any(|held| compare(held, id) == Sic::Yes)
}

/// The index of the one item of `items` whose ids (`ids_of`) hold `id`, `None` for none or
/// several.
fn only_holder_by<'i, T>(
    items: &'i [T],
    ids_of: impl Fn(&'i T) -> Bioseq<'i>,
    id: &ParsedSeqId,
) -> Option<usize> {
    let mut found = items
        .iter()
        .enumerate()
        .filter(|&(_, item)| holds(ids_of(item).ids(), id))
        .map(|(at, _)| at);
    let first = found.next()?;
    found.next().is_none().then_some(first)
}

/// The index of the one Bioseq of `bioseqs` that holds `id`, `None` for none or several.
fn only_holder(bioseqs: &[Bioseq<'_>], id: &ParsedSeqId) -> Option<usize> {
    only_holder_by(bioseqs, |bioseq| *bioseq, id)
}

// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:393
/// The database sequence that NCBI's BLAST database data loader serves for the Seq-align
/// subject id `id`: the first OID of the restricted database's lookup (`seqid_to_oid`), its
/// header read as the search programs read a reported sequence's (`read_header`). Not found is
/// NCBI's ask of its servers (`Refusal::Database`); a lookup that throws, a header with a
/// UTF8String, and a parsed header whose Seq-align id NCBI's lookup may not give as it is
/// (`ParsedSubject::found`, `unfetchable`, `other_holder`, or an id other than `id`) are refused.
fn database_subject(name: &[u8], db: &SeqDb, id: &ParsedSeqId) -> Result<SubjectSource> {
    let lookup = db
        .seqid_lookup()
        .map_err(|error| database_error(error, PROGRAM))?;
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
    let headers = HeaderReader::new(db, false);
    let header = read_header(&headers, ordinal, 0, PROGRAM)?;
    if let ReportedHeader::Parsed(subject) = &header {
        let plain = subject.found().is_none()
            && !subject.unfetchable()
            && subject.other_holder().is_none()
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
    letters: Option<Vec<u8>>,
}

/// The `-subject` Bioseqs as the search programs' records (`TabularSubjects::Records`), by
/// subject index: a Bioseq of one local id has the id of `tabular_subject_id` of its text and
/// title, one without local ids its parsed subject.
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

impl BioseqRecords {
    /// The records of the `-subject` Bioseqs that `sources` name, by subject index.
    fn of(sources: &[SubjectSource], bioseqs: &[Bioseq<'_>], believe: bool) -> Result<Self> {
        let mut records = BioseqRecords {
            ids: Vec::new(),
            lengths: Vec::new(),
            parsed: Vec::new(),
            records: Vec::new(),
            length_of: Box::new(|_| 0),
            protein: Vec::new(),
        };
        for source in sources {
            let SubjectSource::Bioseq(at) = *source else {
                unreachable!("a -subject object's subjects are its Bioseqs")
            };
            let bioseq = bioseqs[at];
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
                        ObjectId::Str(text) => String::from_utf8(text.clone()).map_err(|_| {
                            refused("a -subject Bioseq whose local id is not UTF-8")
                        })?,
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
                        "a -subject Bioseq of a local id and other ids, or of several local ids",
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

/// The values of each subject (`sources`, by subject index) that the fields of `spec` show.
fn subject_entries<'d>(
    db: Option<&SeqDb>,
    spec: &TabularSpec,
    sources: &'d [SubjectSource],
    records: Option<&'d BioseqRecords>,
) -> Result<Vec<SubjectEntry<'d>>> {
    let needs = spec.needs();
    let mut entries = Vec::with_capacity(sources.len());
    for (index, source) in sources.iter().enumerate() {
        // The subject's residue type (`bioseq_protein`): the database's, or the Bioseq's.
        let protein = match (source, db) {
            (SubjectSource::Database { .. }, Some(db)) => db.seqtype() == SeqType::Protein,
            _ => records.is_some_and(|records| records.protein[index]),
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
                            taxonomy,
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
            (SubjectSource::Bioseq(_), _) => {
                let records = records.expect("the records of a -subject object");
                match records.parsed[index].as_deref() {
                    Some(subject) => SubjectValues::of_parsed(
                        subject,
                        spec,
                        records.lengths[index],
                        TabularTaxonomy::default(),
                        protein,
                        false,
                        index,
                        |problem| {
                            refused(&format!(
                                "a -subject Bioseq {} that has {problem}",
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
        entries.push(SubjectEntry {
            values,
            nucleotide: !protein,
            letters: None,
        });
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
        None => subject_taxonomy(Some(&TaxdbSetting::of(inputs.settings)), PROGRAM, [&spec])?
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

    let mut queries: Vec<QueryEntry<'_>> = request
        .queries
        .bioseqs()
        .into_iter()
        .map(|bioseq| QueryEntry::new(bioseq, inputs.believe))
        .collect();
    let subject_bioseqs = request.subject.bioseqs();

    // The Seq-aligns of each query as the formatter reads them (the ungapped split where the
    // result set has alignments), cut to the hit list size, and their subjects found before
    // any byte (NCBI's prefetch comes before the first query's report).
    let mut read: Vec<(Vec<ReadAlign<'_>>, usize)> = Vec::with_capacity(sets.len());
    for set in sets {
        let all = read_set(&alignments[set.clone()], ungapped && !set.is_empty());
        let count = all.len();
        let pruned = super::archive::prune_seqalign(&all, inputs.hitlist)
            .map_err(super::archive::stop_error)?;
        read.push((pruned, count));
    }
    if read.len() > queries.len() {
        return Err(refused("a result set without its query Bioseq"));
    }
    let mut subjects = Subjects::default();
    let mut row_subjects: Vec<Vec<usize>> = Vec::with_capacity(read.len());
    let mut row_queries: Vec<Vec<usize>> = Vec::with_capacity(read.len());
    for (aligns, _) in &read {
        let mut indices = Vec::with_capacity(aligns.len());
        let mut query_indices = Vec::with_capacity(aligns.len());
        for align in aligns {
            if matches!(align.variant(), "std" | "disc") {
                return Err(refused(
                    "Std-seg or discontinuous Seq-aligns (those of a translated search)",
                ));
            }
            let query_id = align.seq_id(0).map_err(super::archive::stop_error)?;
            let in_subjects = subject_bioseqs
                .iter()
                .any(|bioseq| holds(bioseq.ids(), query_id));
            let query = match only_holder_by(&queries, |query| query.bioseq, query_id) {
                Some(at) if !in_subjects => at,
                _ => {
                    return Err(refused(
                        "a query id that names no query Bioseq of the request, or several \
                         Bioseqs",
                    ))
                }
            };
            query_indices.push(query);
            let id = align.seq_id(1).map_err(super::archive::stop_error)?;
            indices.push(subjects.index_of(inputs, &queries, &subject_bioseqs, id)?);
        }
        row_subjects.push(indices);
        row_queries.push(query_indices);
    }
    let records = match inputs.database {
        None => Some(BioseqRecords::of(
            &subjects.sources,
            &subject_bioseqs,
            inputs.believe,
        )?),
        Some(_) => None,
    };
    let db = inputs.database.map(|(_, db)| db);
    let mut entries = subject_entries(db, &spec, &subjects.sources, records.as_ref())?;
    check_failed_lookups(
        Some(&spec),
        row_subjects.iter().map(|indices| {
            let mut runs = indices.clone();
            runs.dedup();
            runs
        }),
        |s_idx| lookup_failure(&entries[s_idx]),
        |s_idx| match &subjects.sources[s_idx] {
            SubjectSource::Database { oid, .. } => failed_database_title(oid),
            SubjectSource::Bioseq(at) => failed_record_title(
                *at,
                records
                    .as_ref()
                    .map_or(&[][..], |records| records.records[s_idx].shown_id()),
            ),
        },
        PROGRAM,
    )?;
    let taxonomies: Vec<TabularTaxonomy> = entries
        .iter()
        .map(|entry| {
            entry
                .values
                .subject(None, entry.nucleotide, -1, -1)
                .taxonomy
                .clone()
        })
        .collect();
    let names: Option<TabularNames> = match inputs.database {
        Some((_, db)) => tabular_names(db, inputs.taxdb, &spec, &taxonomies, PROGRAM)?,
        None => None,
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:725-750
    // `BlastFormatter_PreFetchSequenceData`: the loads of the subjects of every Seq-align of the
    // object, unless `PRE_FETCH_SEQS_LIMIT` is 0 or below their count (`INT_MAX` always loads).
    let seq_aligns: usize = sets.iter().map(|set| set.len()).sum();
    let prefetches = match prefetch_limit {
        None | Some(i32::MAX) => true,
        Some(0) => false,
        Some(limit) => i64::try_from(seq_aligns).unwrap_or(i64::MAX) <= i64::from(limit),
    };
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
                SubjectSource::Bioseq(_) => None,
            }),
            PROGRAM,
        )?,
        None => Vec::new(),
    };
    if !prefetches && !load_texts.is_empty() {
        return Err(refused(
            "texts of the taxonomy database at the loads of the subjects, which \
             PRE_FETCH_SEQS_LIMIT moves to the rows",
        ));
    }

    print_prolog(sink)?;
    sink.post(&load_texts);
    // The round that `PrintOneResultSet` takes is read for the first query: a `Web_StepNumber`
    // that is no integer throws there (`GetInteger`, a serial exception of the loop; probe
    // E1 probe1 O8).
    let Ok(iteration) = inputs.iteration else {
        return Err(super::archive::invalid_archive());
    };

    let letters = reads_letters(&spec, blastn);
    let mut buffers = RowBuffers::default();
    let mut queries_formatted = 0usize;
    for (q, (aligns, count)) in read.iter().enumerate() {
        queries_formatted += 1;
        post_warning(sink);
        let mut out = Vec::new();
        writer.round_query_header(
            &mut out,
            iteration,
            &queries[q].title,
            database_text,
            Some(*count),
        )?;
        let scores = subject_scores(aligns, Some(queries[q].length), &inputs.choice.spec, blastn)
            .map_err(read_stop)?;
        let mut rows = TabularRows::default();
        for (n, align) in aligns.iter().enumerate() {
            let (qcovs, qcovus) = scores.get(n).copied().unwrap_or((-1, -1));
            let s_idx = row_subjects[q][n];
            let q_idx = row_queries[q][n];
            if letters {
                let read = query_letters(&mut queries[q_idx]).and_then(|()| {
                    subject_letters(
                        &mut entries[s_idx],
                        inputs,
                        &subjects.sources[s_idx],
                        &subject_bioseqs,
                    )
                });
                if let Err(end) = read {
                    return Err(end_at_row(sink, &out, end));
                }
            }
            let query = &queries[q_idx];
            let entry = &entries[s_idx];
            let sequences = RowLetters {
                query: query.letters.as_deref().unwrap_or_default(),
                subject: entry.letters.as_deref().unwrap_or_default(),
            };
            let ctx = RowContext {
                spec: &spec,
                no_fetch: blastn,
                query_length: query.length,
            };
            let hsp = hsp(align, &ctx, &sequences, &mut buffers).map_err(read_stop)?;
            let subject = entry
                .values
                .subject(names.as_ref(), entry.nucleotide, qcovs, qcovus);
            let kept = rows.last_found().map(|kept| {
                entries[kept].values.subject(
                    names.as_ref(),
                    entries[kept].nucleotide,
                    qcovs,
                    qcovus,
                )
            });
            // The texts of the row's lookups of the taxonomy database follow the rows before it.
            if let (Some(names), true) = (&names, writer.sets_subject(&subject)) {
                let text = names.row_messages(&spec, subject.taxonomy);
                if !text.is_empty() {
                    sink.write(&out)?;
                    out.clear();
                    sink.post(&text);
                }
            }
            writer.row(
                &mut out,
                &mut rows,
                &query.query(),
                &subject,
                kept.as_ref(),
                &hsp,
            );
        }
        sink.write(&out)?;
        sink.flush()?;
    }
    // `PrintEpilog`: outfmt 7's count of the queries formatted.
    if commented {
        post_warning(sink);
        let mut line = Vec::new();
        writer.epilog(&mut line, queries_formatted)?;
        sink.write(&line)?;
        sink.flush()?;
    }
    Ok(())
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

/// The letters of the row's query, read once.
fn query_letters(query: &mut QueryEntry<'_>) -> std::result::Result<(), LetterEnd> {
    if query.letters.is_none() {
        query.letters = Some(read_bioseq(query.bioseq, "a query Bioseq")?);
    }
    Ok(())
}

/// The letters of the row's subject, read once: a database sequence's IUPAC or NCBIeaa letters,
/// or a `-subject` Bioseq's.
fn subject_letters(
    entry: &mut SubjectEntry<'_>,
    inputs: &ObjectInputs<'_>,
    source: &SubjectSource,
    subject_bioseqs: &[Bioseq<'_>],
) -> std::result::Result<(), LetterEnd> {
    if entry.letters.is_some() {
        return Ok(());
    }
    let letters = match *source {
        SubjectSource::Database { ordinal, .. } => {
            let (_, db) = inputs.database.expect("a database subject");
            let read = if entry.nucleotide {
                db.nucleotide_iupac(ordinal)
            } else {
                db.read_sequence(ordinal).map(|block| {
                    crate::sequence::stdaa_letters::protein_fasta_letters(block.protein(0))
                })
            };
            read.map_err(|error| LetterEnd::Refusal(database_error(error, PROGRAM)))?
        }
        SubjectSource::Bioseq(at) => read_bioseq(subject_bioseqs[at], "a -subject Bioseq")?,
    };
    entry.letters = Some(letters);
    Ok(())
}

#[cfg(test)]
#[path = "tabular_tests.rs"]
mod tests;
