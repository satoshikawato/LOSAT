//! The Bioseq of a query that NCBI's scope finds for a reported database sequence (line RG-KB
//! item 12, decision KB-12-1): one owner of "the query's Bioseq seen for this subject".
//!
//! NCBI's formatter adds the queries of a batch to its scope before the batch's reports and empties
//! it after them; under `-parse_deflines` a query's Bioseq holds the record's ids. A reported
//! database sequence whose Seq-align id is one of those ids (as the scope finds an id,
//! `seqid::defline_set::scope_keys`: a pdb id by its handle, so the query `pdb|1AB1|A` also holds
//! the id of a sequence stored as `pdb` 1AB1 with the chain 321 or 65 and no chain-id) is found as
//! that query in every report of the batch that reads the sequence's Bioseq: the tabular ids are the
//! query's, `stitle` and `salltitles` are `N/A` (the query's Bioseq has no Blast-def-line), `slen`
//! is the query record's length, and the taxonomy fields and the organism report have no taxid for
//! it. A batch of another scope shows the sequence as itself.
//!
//! [`QueryBioseqs`] holds, for each reported pair of a searched query and a database sequence, the
//! query whose Bioseq the scope finds ([`QueryBioseq`], the record as NCBI's scope holds it: a
//! Bioseq of those ids without a Blast-def-line, as a `-subject` record read under
//! `-parse_deflines`, `ParsedSubject::of_record`). The programs fill it where they check the
//! reported sequences against the batch's scope (`check`), and the
//! reports that read the sequence's Bioseq ask it by the query and the ordinal id
//! ([`SeenFrom::of`]); runs without such a sequence keep it empty and ask nothing else.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:836-866

use std::collections::HashMap;
use std::sync::Arc;

use crate::blastinput::fasta_reader::believed::BelievedIds;
use crate::report::db_subject::{AllTabularIds, ParsedSubject};
use crate::report::tabular::{DeflineIds, TabularTaxonomy};

/// A query record as NCBI's scope holds its Bioseq for a database sequence whose Seq-align id is
/// one of the record's ids ([`QueryBioseqs`]).
#[derive(Debug)]
pub struct QueryBioseq {
    /// `sseqid`, `sacc` and `saccver` of the record's ids (`ParsedSubject::tabular_ids`).
    tabular: [Vec<u8>; 3],
    /// `sgi`.
    gi: i64,
    /// `sallseqid`, `sallgi` and `sallacc`: the record's ids as one line.
    all_ids: Vec<DeflineIds>,
    /// `slen`: the letters of the whole record.
    length: usize,
    /// The record's one line in the organism report (outfmt 18): its best id, the title made
    /// from the record's title, no taxid (`ParsedSubject::tax_deflines`).
    tax_deflines: Vec<crate::report::tax_report::TaxDefline>,
}

/// A query record as the programs give it to [`QueryBioseqs::note`]: its ids (`-parse_deflines`),
/// title and the letters of the whole record.
#[derive(Clone, Copy, Debug)]
pub struct QueryRecordRef<'r> {
    pub believed: &'r BelievedIds,
    pub title: &'r [u8],
    pub length: usize,
}

impl<'r> QueryRecordRef<'r> {
    /// The query record `record` as the FASTA reader read it; `None` without `-parse_deflines`.
    pub fn of(record: &'r crate::record::FastaRecord) -> Option<Self> {
        Some(QueryRecordRef {
            believed: record.believed()?,
            title: &record.title,
            length: record.seq().len(),
        })
    }
}

impl QueryBioseq {
    /// The Bioseq of the record `record` (a protein query where `protein`), or
    /// `None` where LOSAT does not write the record's tabular ids as NCBI
    /// (`ParsedSubject::tabular_ids`), which the caller refuses.
    fn new(record: QueryRecordRef<'_>, protein: bool) -> Option<Self> {
        let subject = ParsedSubject::of_record(record.believed, record.title);
        let tabular = subject.tabular_ids(protein).ok()??;
        let AllTabularIds { ids, failed } = subject.all_tabular_ids(protein).ok()?;
        if failed {
            return None;
        }
        Some(QueryBioseq {
            gi: subject.tabular_gi(),
            tax_deflines: subject.tax_deflines(protein),
            tabular,
            all_ids: ids,
            length: record.length,
        })
    }

    /// `sseqid`, `sacc` and `saccver`.
    pub fn tabular_ids(&self) -> &[Vec<u8>; 3] {
        &self.tabular
    }

    /// `sgi`.
    pub fn gi(&self) -> i64 {
        self.gi
    }

    /// The ids of the record's one line (`sallseqid`, `sallgi`, `sallacc`).
    pub fn all_ids(&self) -> &[DeflineIds] {
        &self.all_ids
    }

    /// `slen`.
    pub fn length(&self) -> usize {
        self.length
    }

    /// The taxonomy of the tabular fields: none (`staxid` 0, `staxids` `N/A`), as a FASTA
    /// subject's.
    pub fn taxonomy(&self) -> TabularTaxonomy {
        TabularTaxonomy::default()
    }

    /// The record's line in the organism report (outfmt 18).
    pub fn tax_deflines(&self) -> &[crate::report::tax_report::TaxDefline] {
        &self.tax_deflines
    }
}

/// One owner of the query Bioseqs that NCBI's scope finds for a run's reported database
/// sequences (module docs): by the searched query of the run that reports a sequence and the
/// sequence's ordinal id, the query of its batch whose Bioseq the scope finds.
#[derive(Debug, Default)]
pub struct QueryBioseqs {
    seen: HashMap<(usize, u32), Arc<QueryBioseq>>,
    /// By the first searched query of a scope's span (`QueryScopes::span_of`): the searched
    /// queries of the span by the keys of their ids (`scope_keys`), made at the span's first
    /// sequence that a query holds, so that a batch's queries are walked once.
    holders: HashMap<usize, HashMap<u64, usize>>,
    /// The Bioseqs made, by searched query (`None`: one whose tabular ids LOSAT refuses).
    made: HashMap<usize, Option<Arc<QueryBioseq>>>,
}

/// Why [`QueryBioseqs::note`] does not note a sequence that the scope finds as a query: the
/// caller keeps LOSAT's refusal of it (`common::database::refuse_id_of_a_query`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NotNoted {
    /// No searched query of the span holds the keys (a record that a `-query_loc` range skips,
    /// which NCBI's scope holds too).
    NoHolder,
    /// LOSAT does not write the holder's tabular ids as NCBI (`QueryBioseq::new`).
    Refused,
}

impl QueryBioseqs {
    /// Whether no reported sequence is found as a query.
    pub fn is_empty(&self) -> bool {
        self.seen.is_empty()
    }

    /// The query Bioseq that the scope of the run's searched query `query` finds for the database
    /// sequence `ordinal`, where it finds one.
    pub fn of(&self, query: usize, ordinal: u32) -> Option<&QueryBioseq> {
        if self.seen.is_empty() {
            return None;
        }
        self.seen.get(&(query, ordinal)).map(Arc::as_ref)
    }

    /// The view of a call whose queries are the run's searched queries from `first`.
    pub fn numbered_from(&self, first: usize) -> SeenFrom<'_> {
        SeenFrom {
            seen: (!self.seen.is_empty()).then_some(self),
            first,
        }
    }

    /// Whether every reporting query in `queries` resolves `ordinal` to a noted query Bioseq:
    /// NCBI's scope then finds the query before calling the BLAST database data loader in every
    /// batch that reports `ordinal`.
    pub fn hides_ordinal(&self, queries: impl IntoIterator<Item = usize>, ordinal: u32) -> bool {
        if self.seen.is_empty() {
            return false;
        }
        let mut any = false;
        for query in queries {
            any = true;
            if self.of(query, ordinal).is_none() {
                return false;
            }
        }
        any
    }

    /// Notes that the scope of the run's searched query `query` (whose span of searched queries
    /// is `span`, `QueryScopes::span_of`) finds a query for the database sequence `ordinal`,
    /// whose Seq-align id has the scope keys `keys` (`QueryScopeKeys::own`, which the scope
    /// holds): the one searched query of the span whose ids have one of `keys`. `record(h)` is
    /// the run's searched query `h` for each `h` of `known` (the queries the caller has; an end
    /// of `usize::MAX`: to the run's last searched query, past which it is `None`), a protein
    /// query where `protein`. `Err` where it is not noted ([`NotNoted`]): also where `known`
    /// does not hold the whole span (a query that the caller does not have could hold the id).
    #[allow(clippy::too_many_arguments)]
    pub fn note<'r>(
        &mut self,
        query: usize,
        ordinal: u32,
        keys: &[u64],
        span: std::ops::Range<usize>,
        known: std::ops::Range<usize>,
        record: impl Fn(usize) -> Option<QueryRecordRef<'r>>,
        protein: bool,
    ) -> Result<(), NotNoted> {
        if span.start < known.start || span.end > known.end {
            return Err(NotNoted::NoHolder);
        }
        let holders = self
            .holders
            .entry(span.start)
            .or_insert_with(|| holders_of(span, &record));
        // A key that two queries of the span hold: which Bioseq the scope gives is not
        // determined (`AMBIGUOUS`), LOSAT keeps its refusal.
        let holder = keys
            .iter()
            .find_map(|key| holders.get(key).copied())
            .filter(|&holder| holder != AMBIGUOUS)
            .ok_or(NotNoted::NoHolder)?;
        let made = self
            .made
            .entry(holder)
            .or_insert_with(|| {
                record(holder)
                    .and_then(|found| QueryBioseq::new(found, protein))
                    .map(Arc::new)
            })
            .clone()
            .ok_or(NotNoted::Refused)?;
        self.seen.insert((query, ordinal), made);
        Ok(())
    }
}

/// The holder of a key that two queries of a span hold (`holders_of`).
const AMBIGUOUS: usize = usize::MAX;

/// The searched queries of `span` by the scope keys of their ids (`scope_keys`; `AMBIGUOUS`
/// for a key of two queries), `record(h)` the run's searched query `h` (`None` past the run's
/// last).
fn holders_of<'r>(
    span: std::ops::Range<usize>,
    record: &impl Fn(usize) -> Option<QueryRecordRef<'r>>,
) -> HashMap<u64, usize> {
    let mut holders: HashMap<u64, usize> = HashMap::new();
    let mut id_keys = Vec::new();
    for h in span {
        let Some(found) = record(h) else {
            break;
        };
        id_keys.clear();
        for id in &found.believed.ids {
            crate::seqid::defline_set::scope_keys(id, &mut id_keys);
        }
        id_keys.sort_unstable();
        id_keys.dedup();
        for &key in &id_keys {
            holders
                .entry(key)
                .and_modify(|holder| *holder = AMBIGUOUS)
                .or_insert(h);
        }
    }
    holders
}

/// [`QueryBioseqs`] as a call whose queries are the run's searched queries from `first` asks it
/// (`QueryBioseqs::numbered_from`); the default asks nothing.
#[derive(Clone, Copy, Debug, Default)]
pub struct SeenFrom<'a> {
    seen: Option<&'a QueryBioseqs>,
    first: usize,
}

impl<'a> SeenFrom<'a> {
    /// The query Bioseq that the scope of the call's query `query` finds for the database
    /// sequence `ordinal`.
    #[inline]
    pub fn of(&self, query: usize, ordinal: u32) -> Option<&'a QueryBioseq> {
        self.seen?.of(self.first + query, ordinal)
    }

    /// Whether some reported sequence is found as a query.
    #[inline]
    pub fn any(&self) -> bool {
        self.seen.is_some()
    }

    /// Whether every reporting query (in the call's numbering) in `queries` resolves `ordinal` to
    /// a noted query Bioseq (`QueryBioseqs::hides_ordinal`).
    pub fn hides_ordinal(&self, queries: impl IntoIterator<Item = usize>, ordinal: u32) -> bool {
        let Some(seen) = self.seen else {
            return false;
        };
        seen.hides_ordinal(queries.into_iter().map(|query| self.first + query), ordinal)
    }
}

/// The frame that NCBI's tabular formatter computes for a translated subject whose Bioseq handle
/// has length `bioseq_len`, from its 1-based alignment start `start1` (`s_start` in `tabular.cpp`,
/// which is the higher end on a minus strand) and strand `plus`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1079-1083
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1233-1245
pub fn ncbi_frame(plus: bool, start1: i32, bioseq_len: usize) -> i32 {
    let start0 = start1.wrapping_sub(1);
    if plus {
        (start0 % 3) + 1
    } else {
        let diff = (bioseq_len as i32).wrapping_sub(start0).wrapping_sub(1);
        -((diff % 3) + 1)
    }
}

// NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnvec.hpp:288-310
// NCBI reference (598d8ae6): c++/include/objmgr/seq_vector_ci.hpp:517-525
// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:631-644,720-728
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:893-901
/// Whether TBLASTX's `qseq` throws `CAlnException::eTranslateFailure` when `GetWholeAlnSeqString`
/// translates the subject row's alignment slice (`[start0, end)` in 0-based nucleotide coordinates,
/// strand `plus`) of a found query Bioseq of length `seen_len`: `CSeqVector` clips at `seen_len`
/// (empty when `start0 >= seen_len` on plus, or when a minus-strand slice passes `seen_len` and
/// its reversed start wraps below zero), so only a plus-strand HSP straddling `seen_len` with a
/// clipped length not divisible by 3 throws.
pub fn tblastx_qseq_throws(plus: bool, start0: usize, end: usize, seen_len: usize) -> bool {
    plus && start0 < seen_len && seen_len < end && !(seen_len - start0).is_multiple_of(3)
}

/// One report of a run as the table of line RG item X7 sorts it: what it shows of a reported
/// database sequence that NCBI's scope finds as a query (`FoundWriters`).
#[derive(Clone, Copy, Debug)]
pub(crate) enum FoundWriter<'s> {
    /// The pairwise report (outfmt 0, also with `-html`): the found query's title and letters.
    Pairwise,
    /// The query-anchored formats (outfmt 1 to 4): the found query's title and letters.
    Anchored,
    /// BLAST XML, XML2 and JSON (outfmt 5, 13 to 16): the found query's title and letters.
    Structured,
    /// A table (outfmt 6, 7, 10 and 20) or the organism report (outfmt 18, a `TabularSpec` of
    /// `TabularKind::TaxReport`) of these fields (`None`: the program's default fields of outfmt
    /// 6, 7 and 10).
    Table(Option<&'s crate::report::tabular::TabularSpec>),
    /// The organism report (outfmt 18) of a program that keeps no field list for it.
    TaxReport,
    /// SAM (outfmt 17): the target row's Bioseq, which is the subject's without `SR`, gives the
    /// CIGAR's trailing hard clip its length and, with `SQ`, SEQ its letters; with `SR` the
    /// subject's gives `@SQ LN:` its length (`target_letters`: `SQ` without `SR`).
    Sam { target_letters: bool },
    /// The Seq-annots and the archive (outfmt 8, 9, 11 and 12), which write the Seq-align id alone.
    Annots,
}

/// A run's reports as one owner of what they show for a reported database sequence that NCBI's
/// scope finds as a query (line RG item X7, decision XT-X7-1b): by the query's Seq-align id
/// (`query_batch::BatchQueries`) under `-parse_deflines` (`shows_shown`: KB-12's substitution,
/// `QueryBioseqs::note`) or by a generated id (`generated_effect`).
#[derive(Clone, Debug)]
pub(crate) struct FoundWriters<'s> {
    pub program: crate::blastinput::app::Program,
    pub writers: Vec<FoundWriter<'s>>,
    /// Whether the search itself reads the found query's letters: psiblast's calls that make a
    /// PSSM of a round's alignments, whose next round then differs in every format.
    pub search_reads: bool,
}

impl FoundWriters<'static> {
    /// No report: a run whose reports keep LOSAT's refusal of every such sequence (`SeenRun::REFUSED`).
    pub(crate) const NONE: Self = FoundWriters {
        program: crate::blastinput::app::Program::Blastp,
        writers: Vec::new(),
        search_reads: true,
    };
}

/// What NCBI does with a sequence found by a generated id in the pairwise report.
pub(crate) const PAIRWISE_EFFECT: &str =
    "shows that query in its place in the pairwise report (outfmt 0)";

/// What NCBI does with it in the query-anchored formats.
pub(crate) const QUERY_ANCHORED_EFFECT: &str =
    "shows that query in its place in the query-anchored formats (outfmt 1 to 4)";

/// What NCBI does with it in BLAST XML, XML2 and JSON (RPSBLAST gets the query's title for the
/// entry, BLASTX and RPSTBLASTN an exception from the Seq-align's conversion and half a document
/// or none; PS15a-1, PS15a-2); with outfmt 5 the formatter keeps its scope for the run, where a
/// later batch's query with that id is searched as the database sequence (`b1.bp.q2x.f5`).
pub(crate) const XML_EFFECT: &str = "takes that query for it in BLAST XML, XML2 and JSON (outfmt \
     5, 13 to 16), and with outfmt 5, whose scope it keeps for the run, searches a later batch's \
     query with that id as the database sequence";

/// What NCBI does with it in the tables and the organism report where they differ
/// (`tabular_differs_for_generated_id`).
pub(crate) const TABLE_EFFECT: &str = "takes that query for it in the tables (outfmt 6, 7, 10 and \
     20) whose fields show what it finds by that id and in the organism report (outfmt 18) of a \
     sequence with a taxid where a taxonomy database is found (BLASTX's and RPSTBLASTN's tables \
     can stop with an exception)";

/// What NCBI does with it in psiblast's rounds after the first and in a saved PSSM.
pub(crate) const SEARCH_EFFECT: &str = "makes the next round's PSSM, or the PSSM it saves, from \
     that query's letters in its place (psiblast)";

/// What NCBI does with it in SAM where the found query's length or letters are not the sequence's.
pub(crate) const SAM_EFFECT: &str = "takes that query's length (and with SQ its letters) for it \
     in SAM (outfmt 17)";

/// What NCBI does with it where the run's loads of a sequence with a taxid write texts of the
/// taxonomy database (`common::database::TaxdbRun::loads_write_texts`).
pub(crate) const LOAD_EFFECT: &str =
    "loads nothing of it, where LOSAT's load writes the taxonomy database's messages";

impl FoundWriters<'_> {
    /// Whether every report shows NCBI's bytes for a sequence that the scope finds by its own
    /// Seq-align id under `-parse_deflines` (KB-12, `QueryBioseqs::note`): a table that reads no
    /// letter of the found Bioseq and the organism report (`tabular_shows`), the Seq-annots and
    /// the archive. The callers add their programs' conditions.
    pub(crate) fn shows_shown(&self) -> bool {
        !self.search_reads
            && self.writers.iter().all(|writer| match *writer {
                FoundWriter::Annots | FoundWriter::TaxReport => true,
                FoundWriter::Table(spec) => tabular_shows(spec, self.program),
                _ => false,
            })
    }

    /// What NCBI does that LOSAT does not for the reported database sequence `subject` that the
    /// scope finds by a generated id (`query_batch::GeneratedQueries`), for an HSP with the
    /// Seq-align's identities `num_ident` (1 where the caller does not have them: `nident` then
    /// counts as read), with the taxonomy database found (`taxdb`) and `found()`, how the found
    /// query's record compares with the sequence over the HSPs ([`FoundMatch`], asked only where a
    /// report reads the found letters or length): the effect of the first report that differs,
    /// `None` where every report is NCBI's bytes (the Seq-annots and the archive; SAM where the
    /// found query has the sequence's length, and with `SQ` without `SR` its letters over the
    /// HSPs; tables of fields that read nothing of the found query, or only letters and a length
    /// that are the sequence's own there; the organism report of a sequence without a taxid, or
    /// without a taxonomy database; psiblast's PSSM of letters that are the sequence's over the
    /// HSPs; probes of line RG item X7, `x71/probe2`, `x71/fix2/ev2`).
    pub(crate) fn generated_effect(
        &self,
        subject: &ParsedSubject,
        num_ident: usize,
        taxdb: bool,
        found: &dyn Fn() -> FoundMatch,
    ) -> Option<&'static str> {
        if self.search_reads && !found().aligned {
            return Some(SEARCH_EFFECT);
        }
        let taxid = subject.first_source_taxid().is_some();
        self.writers.iter().find_map(|writer| match *writer {
            FoundWriter::Pairwise => Some(PAIRWISE_EFFECT),
            FoundWriter::Anchored => Some(QUERY_ANCHORED_EFFECT),
            FoundWriter::Structured => Some(XML_EFFECT),
            FoundWriter::Table(spec) => tabular_differs_for_generated_id(
                spec,
                self.program,
                num_ident,
                (taxid, taxdb),
                found,
            )
            .then_some(TABLE_EFFECT),
            FoundWriter::TaxReport => (taxid && taxdb).then_some(TABLE_EFFECT),
            FoundWriter::Sam { target_letters } => {
                let found = found();
                (!found.length || (target_letters && !found.aligned)).then_some(SAM_EFFECT)
            }
            FoundWriter::Annots => None,
        })
    }
}

/// What a table's fields read of a reported database sequence that NCBI's scope finds as a
/// query (one field rule for both kinds of ids, `tabular_shows` and
/// `tabular_differs_for_generated_id`).
struct FoundFields {
    /// Fields that show what the scope finds beyond the Seq-align id: its ids and titles
    /// (`sallseqid`, `sallacc`, `stitle`, `salltitles`) and `staxids` (`N/A` for a query).
    names: bool,
    /// Fields of the found Bioseq's length: `slen`, and in TBLASTN and TBLASTX the frame that it
    /// gives (`sframe`, `frames`).
    length: bool,
    /// Fields that read the found Bioseq's letters (`tabular_shows`).
    letters: bool,
    /// The sequence's taxid (`staxid`: `N/A` for a query).
    taxid: bool,
    /// Fields of the taxonomy database's names of the sequence's taxid (the names, the kingdoms:
    /// `N/A` for a query, and without a taxonomy database for the sequence too; probe
    /// `x71/fix2/rps18`).
    taxonomy: bool,
    /// In BLASTX (and RPSTBLASTN, whose driver is BLASTX's) a field that enters the alignment
    /// block, which throws (`TabularSpec::opens_counts`).
    throws: bool,
}

/// `FoundFields` of `spec` (`Some`) in `program`, with the Seq-align's identities `num_ident`.
fn found_fields(
    spec: &crate::report::tabular::TabularSpec,
    program: crate::blastinput::app::Program,
    num_ident: usize,
) -> FoundFields {
    use crate::blastinput::app::Program;
    use crate::report::tabular::TabularField::*;
    let translated = matches!(program, Program::Tblastn | Program::Tblastx);
    FoundFields {
        names: spec.has_any(&[
            SubjectAllSeqIds,
            SubjectAllAccessions,
            SubjectTitle,
            SubjectAllTitles,
            SubjectTaxIds,
        ]),
        length: spec.has(SubjectLength) || (translated && spec.has_any(&[SubjectFrame, Frames])),
        letters: spec.has_any(&[SubjectSeq, Positives, PercentPositives, Btop])
            || (program != Program::Blastn
                && (spec.has_any(&[Mismatches, PercentIdentical])
                    || (spec.has(NumIdentical) && num_ident > 0)))
            || (program == Program::Blastn
                && spec.has(QuerySeq)
                && spec.has_any(&[Mismatches, PercentIdentical, NumIdentical]))
            || (program == Program::Tblastn && spec.has(QuerySeq)),
        taxid: spec.has(SubjectTaxId),
        taxonomy: spec.has_any(&[
            SubjectSciName,
            SubjectCommonName,
            SubjectBlastName,
            SubjectSuperKingdom,
            SubjectSciNames,
            SubjectCommonNames,
            SubjectBlastNames,
            SubjectSuperKingdoms,
        ]),
        throws: program == Program::Blastx && spec.opens_counts(&num_ident),
    }
}

/// Whether `program`'s tabular report of `spec` (`None`: the default fields; the organism report
/// too) shows NCBI's bytes for a database sequence that the scope finds as a query by its own
/// Seq-align id under `-parse_deflines` (KB-12's substitution, `FoundWriters::shows_shown`):
/// - In BLASTX, every field that enters the alignment block (`TabularSpec::opens_counts`:
///   `qstart`, `qend`, `sstart`, `send`, `length`, `gapopen`, `gaps`, `qframe`, `sstrand`, and
///   the letter-reading fields) converts the alignment with `CreateTranslatedDensegFromNADenseg`,
///   which throws because the nucleotide query Bioseq makes `subject_is_na` true on a Dense-seg
///   that already has widths: LOSAT refuses those fields. `sframe` and `frames` alone stay
///   outside the block (`0`, `0/0`).
/// - Elsewhere, fields that read the found Bioseq's letters (`sseq`, `btop`, `positive`, `ppos`;
///   outside BLASTN also `mismatch`, `pident`, `nident`; in BLASTN `qseq` combined with
///   `mismatch`, `pident` or `nident`, which enters the sequence block and recomputes identities
///   from the found Bioseq; in TBLASTN `qseq`, which translates the protein query Bioseq and
///   throws) keep LOSAT's refusal. `qseq` alone in BLASTP, BLASTN and (except when a plus-strand
///   subject HSP straddles the found query Bioseq's length with a clipped length not divisible by
///   3, [`tblastx_qseq_throws`]) TBLASTX reads the real query's letters; `sframe` and `frames` in
///   TBLASTX and TBLASTN compute the frame from `QueryBioseq::length`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-1085
pub(crate) fn tabular_shows(
    spec: Option<&crate::report::tabular::TabularSpec>,
    program: crate::blastinput::app::Program,
) -> bool {
    use crate::blastinput::app::Program;
    let Some(spec) = spec else {
        // `6 qaccver saccver pident length mismatch ...`: BLASTN alone does not read letters.
        return program == Program::Blastn;
    };
    if spec.kind == crate::report::tabular::TabularKind::TaxReport {
        return true;
    }
    let fields = found_fields(spec, program, 1);
    if program == Program::Blastx {
        return !fields.throws;
    }
    !fields.letters
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:891-924
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:958,985-986
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:964-968
// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:996-1000
/// Whether `program`'s table of `spec` (`None`: the default fields of outfmt 6, 7 and 10) or
/// organism report differs for an HSP with the Seq-align's identities `num_ident` whose subject's
/// Seq-align id is a generated id (`lcl|Query_N`, or a PSSM query's local id) of a query of its
/// batch (`program` is the driver's: BLASTP for RPSBLAST and PSIBLAST, BLASTX for RPSTBLASTN).
/// NCBI's formatter finds that query by the id, so its row differs where it shows what it found
/// beyond the id (`FoundFields::names`), the found query's length (`length`) unless it is the
/// sequence's (`found().length`), the letters of the alignment's subject row (`letters`;
/// `nident` where the block opens for it, `num_ident > 0`) unless the found query has the
/// sequence's letters over the HSP (`found().aligned`: NCBI then reads the letters that
/// LOSAT reads of the sequence; BLASTP `q` against its copy `lcl|Query_1`,
/// `tests/non_utf8_names` R3_7 and R7; a found query shorter than the HSP's end reads gaps there,
/// TBLASTX `0.000`, probe `x71/fix2/ev2` `tx.tail.f6`), TBLASTX's `qseq` where an HSP straddles
/// the found query's end (`found().unstraddled`; it translates the found Bioseq and can throw,
/// `tblastx_qseq_throws`), and the sequence's taxid where it has
/// one (`subject_taxid`; a query has none: `N/A`), with a taxonomy database (`taxdb`) also its
/// names; in BLASTX and RPSTBLASTN, whose found query is nucleotide, every row that enters the
/// alignment block ends the run with an exception of the Seq-align's conversion. The organism
/// report differs for a subject with a taxid where a taxonomy database is found: without one,
/// NCBI lists neither (probe `tx.bp.qid.f18`; rows `db.x7.rps.of18coltax` and
/// `db.x7.rpstbn.of18coltax`). The subject's own ids (`sseqid`, `sacc`, `saccver`, `sgi`,
/// `sallgi`), the query's letters, the scores and the coordinates are NCBI's bytes (review S1
/// of PS15a-2; probes of line RG item X7, `x71/probe2`: BLASTN's default fields are, its
/// `sallseqid` and `stitle` are not).
pub(crate) fn tabular_differs_for_generated_id(
    spec: Option<&crate::report::tabular::TabularSpec>,
    program: crate::blastinput::app::Program,
    num_ident: usize,
    (subject_taxid, taxdb): (bool, bool),
    found: &dyn Fn() -> FoundMatch,
) -> bool {
    use crate::blastinput::app::Program;
    let Some(spec) = spec else {
        // `pident` and `mismatch` read the letters, but not BLASTN's; BLASTX's block throws.
        return match program {
            Program::Blastn => false,
            Program::Blastx => true,
            _ => !found().aligned,
        };
    };
    if spec.kind == crate::report::tabular::TabularKind::TaxReport {
        return subject_taxid && taxdb;
    }
    let fields = found_fields(spec, program, num_ident);
    let tblastx_qseq =
        program == Program::Tblastx && spec.has(crate::report::tabular::TabularField::QuerySeq);
    fields.names
        || (subject_taxid && (fields.taxid || (taxdb && fields.taxonomy)))
        || fields.throws
        || (fields.length && !found().length)
        || (fields.letters && !found().aligned)
        || (tblastx_qseq && !found().unstraddled)
}

/// How the query record that NCBI's scope finds by the generated id of a reported database
/// sequence compares with that sequence (`compare_found`; `FoundMatch::default()`: neither, where
/// the found query is not compared, as a query of the other molecule).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct FoundMatch {
    /// The record has the sequence's length (`slen`; the frames of TBLASTN and TBLASTX; SAM's
    /// trailing hard clip and `@SQ LN:`).
    pub length: bool,
    /// Over the subject range of each HSP the record has the sequence's letters (the fields that
    /// read the letters of the alignment's subject row; SAM's SEQ; psiblast's PSSM).
    pub aligned: bool,
    /// No HSP's subject range straddles the record's end (TBLASTX's `qseq`, which throws only
    /// where one does, `tblastx_qseq_throws`; a range wholly past the end translates nothing).
    pub unstraddled: bool,
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:946-947 (a masked residue is
// kept in upper case)
// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector.cpp:1082-1085 (a gap reads `-` in ncbieaa)
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:343,411-452 (the PSSM
// reads the subject over its alignment)
/// How the query record of length `length` whose whole letters are `record` (as the FASTA reader
/// keeps them: IUPAC, in lower case where the mask covers them; empty where the caller dropped
/// them, BLASTN's records of a run whose reports read no query letter, `x_query_drop`, which keeps
/// them wherever a report reads the found query's letters) compares with the database sequence
/// `ordinal` of `db` over the subject ranges `ranges` (0-based, end exclusive, on its plus
/// strand) of the reported HSPs, as NCBI's formatter reads the found query's Bioseq: its residues
/// in upper case, a protein's as ncbieaa letters (`tabular.cpp`'s `SetAaCoding`). A protein
/// record's `X` is never the sequence's: a `>?` gap line reads as `X` in the record and `-` in
/// NCBI's Bioseq. A range past the record's end, or past the letters kept, is not the sequence's
/// (NCBI reads gaps there). A sequence that the database cannot read compares as neither (LOSAT
/// refuses).
pub(crate) fn compare_found(
    record: &[u8],
    length: usize,
    db: &crate::seqdb::SeqDb,
    ordinal: u32,
    ranges: &[std::ops::Range<usize>],
) -> FoundMatch {
    let compare = |letters: &[u8], protein: bool| FoundMatch {
        length: letters.len() == length,
        unstraddled: ranges
            .iter()
            .all(|range| !(range.start < length && length < range.end)),
        aligned: ranges.iter().all(|range| {
            match (record.get(range.clone()), letters.get(range.clone())) {
                (Some(found), Some(letters)) => {
                    found.iter().zip(letters).all(|(&residue, &letter)| {
                        let residue = residue.to_ascii_uppercase();
                        residue == letter && !(protein && residue == b'X')
                    })
                }
                _ => false,
            }
        }),
    };
    match db.seqtype() {
        crate::seqdb::SeqType::Protein => db.read_sequence(ordinal).map_or_else(
            |_| FoundMatch::default(),
            |block| {
                compare(
                    &crate::sequence::stdaa_letters::protein_fasta_letters(block.protein(0)),
                    true,
                )
            },
        ),
        crate::seqdb::SeqType::Nucleotide => db.nucleotide_iupac(ordinal).map_or_else(
            |_| FoundMatch::default(),
            |letters| compare(&letters, false),
        ),
    }
}

/// The query record of `records` (a batch's input records, read without `-parse_deflines`) that
/// NCBI's scope finds by the Seq-align id of the reported sequence `subject`: the record whose
/// local id (`Query_N`) is that local string id, as the scope compares local ids, without case
/// (`query_batch::GeneratedQueries::finds`).
pub(crate) fn found_record<'r>(
    subject: &ParsedSubject,
    records: &'r [crate::record::FastaRecord],
) -> Option<&'r crate::record::FastaRecord> {
    records.get(found_index(subject, records)?)
}

/// The index in `records` of [`found_record`].
pub(crate) fn found_index(
    subject: &ParsedSubject,
    records: &[crate::record::FastaRecord],
) -> Option<usize> {
    let crate::seqid::ParsedSeqId::Local(crate::seqid::model::ObjectId::Str(text)) =
        subject.alignment_id()
    else {
        return None;
    };
    records
        .iter()
        .position(|record| record.local_id.as_bytes().eq_ignore_ascii_case(text))
}

/// What `check` reads of a run besides the sequence: whether its reports show NCBI's bytes for a
/// sequence that the scope finds as a query by its `-parse_deflines` id (`shows`,
/// `FoundWriters::shows_shown` of the run's reports, and no text that the taxonomy database
/// writes at a load, `common::database::TaxdbRun::loads_write_texts`), the queries the caller has
/// (`known`, ends `usize::MAX` where they reach the run's last searched query) and `record(h)` of
/// each, whether the queries are proteins (`protein`), the run's reports (`writers`), whether
/// its taxonomy database is found (`taxdb`), whether a load of a sequence with a taxid writes
/// the taxonomy database's texts (`load_texts`, `loads_write_texts`) and `found(subject,
/// ordinal, ranges)`, how the query record that the scope finds by the generated id of the
/// sequence `subject` (ordinal id `ordinal`) compares with it over the HSPs' subject ranges
/// `ranges` ([`compare_found`]; [`NOT_COMPARED`] where the found query's molecule is not the
/// sequence's), for a sequence found by a generated id.
pub(crate) struct SeenRun<'s, 'r> {
    pub shows: bool,
    pub known: std::ops::Range<usize>,
    pub record: &'s dyn Fn(usize) -> Option<QueryRecordRef<'r>>,
    pub protein: bool,
    pub writers: &'s FoundWriters<'s>,
    pub taxdb: bool,
    pub load_texts: bool,
    pub found: &'s FoundCompare<'s>,
}

/// `SeenRun::found`: how the query record found by the generated id of a reported sequence
/// compares with it (`compare_found`).
pub(crate) type FoundCompare<'s> =
    dyn Fn(&ParsedSubject, u32, &[std::ops::Range<usize>]) -> FoundMatch + 's;

/// `SeenRun::found` of a run whose found query is never compared: a query of the other molecule
/// (BLASTX, RPSTBLASTN, TBLASTN), or a run whose reports keep LOSAT's refusal.
pub(crate) fn not_compared(_: &ParsedSubject, _: u32, _: &[std::ops::Range<usize>]) -> FoundMatch {
    FoundMatch::default()
}

/// [`not_compared`] as `SeenRun::found`.
pub(crate) const NOT_COMPARED: &FoundCompare<'static> = &not_compared;

/// The HSPs of a reported sequence as `check` reads them: the Seq-align's identities
/// (`num_ident`, 1 where the caller does not have them: `nident` then counts as read) and their
/// subject ranges (`ranges`, 0-based on the plus strand, asked only where a report reads the found
/// query's letters).
pub(crate) struct ReportedHsps<'h> {
    pub num_ident: usize,
    pub ranges: &'h dyn Fn() -> Vec<std::ops::Range<usize>>,
}

impl SeenRun<'static, 'static> {
    /// A run whose reports keep LOSAT's refusal of every such sequence.
    pub(crate) const REFUSED: Self = SeenRun {
        shows: false,
        known: 0..0,
        record: &no_record,
        protein: false,
        writers: &FoundWriters::NONE,
        taxdb: false,
        load_texts: true,
        found: NOT_COMPARED,
    };
}

/// No query record (`SeenRun::REFUSED`).
fn no_record(_: usize) -> Option<QueryRecordRef<'static>> {
    None
}

/// The one check of line RG item X7 (decision XT-X7-1b) of the reported database sequence
/// `subject` (`oid` its OID, `ordinal` its ordinal id) of the run's searched query `query`, against
/// the queries that NCBI's scope holds for the query's batch (`queries`):
/// - `-parse_deflines` (`BatchQueries::Shown`): `common::database::refuse_id_of_a_query`, unless
///   the scope finds a query of its batch for it by its own Seq-align id and the run's reports
///   show NCBI's bytes for such a sequence (`run.shows`, `SeenRun`): that query is then noted in
///   `seen` ([`QueryBioseqs::note`]) and the reports take its Bioseq for the sequence. LOSAT
///   keeps its refusal where the sequence's load posts the warning of a UTF8String header
///   (`ParsedSubject::scope_utf8_read`: NCBI loads nothing of a sequence that its scope finds as
///   a query), where the id matches a query weakly, and where the holder is not one searched
///   query of the scope that the caller has, or LOSAT does not write its tabular ids
///   ([`NotNoted`]).
/// - A generated id (`BatchQueries::Generated`), for the HSPs `hsps` of the sequence: LOSAT's
///   refusal where a report differs for it (`FoundWriters::generated_effect`), where the sequence
///   has a taxid and the run's loads write the taxonomy database's texts (`SeenRun::load_texts`:
///   NCBI loads nothing of it; probe `bp.tax.badmagic` of `x71/fix2/ev2`), or where its load posts
///   the warning of a UTF8String header; else nothing (every report is NCBI's bytes).
#[allow(clippy::too_many_arguments)]
pub(crate) fn check(
    seen: &mut QueryBioseqs,
    subject: &ParsedSubject,
    oid: impl Fn() -> crate::seqdb::Oid,
    (query, ordinal): (usize, u32),
    hsps: &ReportedHsps<'_>,
    queries: &crate::blastinput::query_batch::BatchQueries<'_>,
    run: &SeenRun<'_, '_>,
    program: &str,
) -> anyhow::Result<()> {
    use crate::blastinput::query_batch::BatchQueries;
    let scopes = match queries {
        BatchQueries::Shown(scopes) => scopes,
        BatchQueries::Generated(generated) => {
            if !generated.finds(subject) {
                return Ok(());
            }
            let found = std::cell::OnceCell::new();
            let effect = run
                .writers
                .generated_effect(subject, hsps.num_ident, run.taxdb, &|| {
                    *found.get_or_init(|| (run.found)(subject, ordinal, &(hsps.ranges)()))
                })
                .or_else(|| {
                    (run.load_texts && subject.first_source_taxid().is_some())
                        .then_some(LOAD_EFFECT)
                })
                .or_else(|| {
                    subject.scope_utf8_read().is_some().then_some(
                        "loads nothing of it where LOSAT posts the warning of its header",
                    )
                });
            return match effect {
                Some(effect) => Err(generated_refusal(subject, oid(), program, effect)),
                None => Ok(()),
            };
        }
    };
    let scope = scopes.of_query(query);
    let keys = &subject.query_scope_keys().own;
    if run.shows
        && scope.has_key(keys)
        && subject.scope_utf8_read().is_none()
        && seen
            .note(
                query,
                ordinal,
                keys,
                scopes.span_of(query),
                run.known.clone(),
                run.record,
                run.protein,
            )
            .is_ok()
    {
        return Ok(());
    }
    crate::algorithm::common::database::refuse_id_of_a_query(subject, oid, scope, program)
}

/// What the refusals of a sequence found by a generated id call that id: the local id that LOSAT
/// generates for a FASTA query read without `-parse_deflines` (`Query_N`), or the id of the PSSM
/// query of `-in_pssm` that a report shows as a local id (`lcl|myprot`).
pub(crate) const GENERATED_QUERY_ID: &str = "the id of a query (generated, or the PSSM's)";

/// LOSAT's refusal of the reported database sequence `subject` (OID `oid`) whose Seq-align id
/// NCBI's scope finds as the generated id of a query of the batch that reports it, where NCBI
/// does `effect` with it (`FoundWriters::generated_effect`).
pub(crate) fn generated_refusal(
    subject: &ParsedSubject,
    oid: crate::seqdb::Oid,
    program: &str,
    effect: &str,
) -> anyhow::Error {
    crate::algorithm::common::database::parsed_refusal(
        oid,
        &format!(
            "has the Seq-id {}, {GENERATED_QUERY_ID} in the query batch that reports it: NCBI \
             BLAST+ {effect}",
            String::from_utf8_lossy(&subject.alignment_fasta())
        ),
        program,
    )
}

#[cfg(test)]
#[path = "query_bioseq_tests.rs"]
mod tests;
