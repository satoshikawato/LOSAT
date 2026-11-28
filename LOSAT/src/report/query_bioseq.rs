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
//! reported sequences against the batch's scope (`note_or_refuse`), and the
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

    /// Whether the scope of some query of the run finds a query for the database sequence
    /// `ordinal` (a noted pair of any query): NCBI's data loader may then not load it, so the
    /// lines of its record may not hold ids for others (`common::database::ShownHolders::shown`).
    pub fn finds(&self, ordinal: u32) -> bool {
        !self.seen.is_empty() && self.seen.keys().any(|&(_, seen)| seen == ordinal)
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

    /// Whether the scope of some query of the run finds a query for the database sequence
    /// `ordinal` (`QueryBioseqs::finds`).
    pub fn finds(&self, ordinal: u32) -> bool {
        self.seen.is_some_and(|seen| seen.finds(ordinal))
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

/// Whether `program`'s tabular report of `spec` (`None`: the default fields; the organism report
/// too) shows NCBI's bytes for a database sequence that the scope finds as a query:
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
    use crate::report::tabular::TabularField::*;
    let Some(spec) = spec else {
        // `6 qaccver saccver pident length mismatch ...`: BLASTN alone does not read letters.
        return program == Program::Blastn;
    };
    if spec.kind == crate::report::tabular::TabularKind::TaxReport {
        return true;
    }
    if program == Program::Blastx {
        return !spec.opens_counts(&1);
    }
    let reads_subject_letters = spec.has_any(&[SubjectSeq, Positives, PercentPositives, Btop])
        || (program != Program::Blastn
            && spec.has_any(&[Mismatches, PercentIdentical, NumIdentical]))
        || (program == Program::Blastn
            && spec.has(QuerySeq)
            && spec.has_any(&[Mismatches, PercentIdentical, NumIdentical]))
        || (program == Program::Tblastn && spec.has(QuerySeq));
    !reads_subject_letters
}

/// What `note_or_refuse` reads of a run besides the sequence: whether its reports show NCBI's
/// bytes for a sequence that the scope finds as a query (`shows`, `tabular_shows` of every
/// report, and no text that the taxonomy database writes at a load,
/// `common::database::TaxdbRun::loads_write_texts`), the queries the caller has (`known`, ends
/// `usize::MAX` where they reach the run's last searched query) and `record(h)` of each, and
/// whether the queries are proteins (`protein`).
pub(crate) struct SeenRun<'s, 'r> {
    pub shows: bool,
    pub known: std::ops::Range<usize>,
    pub record: &'s dyn Fn(usize) -> Option<QueryRecordRef<'r>>,
    pub protein: bool,
}

impl SeenRun<'static, 'static> {
    /// A run whose reports keep LOSAT's refusal of every such sequence.
    pub(crate) const REFUSED: Self = SeenRun {
        shows: false,
        known: 0..0,
        record: &no_record,
        protein: false,
    };
}

/// No query record (`SeenRun::REFUSED`).
fn no_record(_: usize) -> Option<QueryRecordRef<'static>> {
    None
}

/// `common::database::refuse_id_of_a_query` of the reported database sequence `subject` (`oid`
/// its OID, `ordinal` its ordinal id) of the run's searched query `query`, unless the scope of the
/// query (`scopes`) finds a query of its batch for it by its own Seq-align id and the run's
/// reports show NCBI's bytes for such a sequence (`run.shows`, `SeenRun`): that query is then
/// noted in `seen` ([`QueryBioseqs::note`]) and the reports take its Bioseq for the sequence.
/// LOSAT keeps its refusal where the sequence's load posts the warning of a UTF8String header
/// (`ParsedSubject::scope_utf8_read`: NCBI loads nothing of a sequence that its scope finds as a
/// query), where the id matches a query weakly, and where the holder is not one searched query
/// of the scope that the caller has, or LOSAT does not write its tabular ids ([`NotNoted`]).
#[allow(clippy::too_many_arguments)]
pub(crate) fn note_or_refuse(
    seen: &mut QueryBioseqs,
    subject: &ParsedSubject,
    oid: impl Fn() -> crate::seqdb::Oid,
    query: usize,
    ordinal: u32,
    scopes: &crate::blastinput::query_batch::QueryScopes,
    run: &SeenRun<'_, '_>,
    program: &str,
) -> anyhow::Result<()> {
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

#[cfg(test)]
#[path = "query_bioseq_tests.rs"]
mod tests;
