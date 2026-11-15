//! TBLASTX's sequences for the query-anchored reports (outfmt 1 to 4, `report::query_anchored`):
//! the letters, labels and masks of a query's section, read where outfmt 0's rows read them
//! (`report::PairwiseRows`).
//!
//! - Letters: the searched nucleotides of the queries and of the subjects (the records' letters,
//!   or a packed subject's, `report::SubjectLetters`), upper case, in their records' positions (the
//!   searched letters start at the `-query_loc` or `-subject_loc` interval); the section translates
//!   the query's rows with the query's genetic code (`-query_gencode`, NCBI's
//!   `m_MasterGeneticCode`) and the subjects' with the database's (`-db_gencode`,
//!   `m_SlaveGeneticCode`).
//! - Labels: `x_FillSeqid` of the Bioseq's Seq-ids (`query_anchored::record_label`,
//!   `query_anchored::database_label`).
//! - Masks: the query's SEG and `-lcase_masking` masks of every frame (`report::query_frame_masks`,
//!   as outfmt 0's rows show them), which a group's query row shows where they have its frame.
//!   The subjects' masks have no frame (`CSeqLocInfo::eFrameNotSet`), and a subject's row of
//!   width 3 has one, so NCBI shows none of them (`query_anchored::prepare_rows`): none are read.

use super::blast_engine::ReportDatabase;
use super::report::{DisplayCode, SubjectLetters};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::report::query_anchored::{
    database_label, query_record_letters, record_label, Mask, ReportSequences,
};

/// What the sections of a TBLASTX report read (module doc).
pub(crate) struct AnchoredSequences<'a> {
    pub(crate) queries: &'a [FastaRecord],
    pub(crate) query_placements: &'a Placements,
    pub(crate) subjects: &'a [FastaRecord],
    pub(crate) subject_letters: SubjectLetters<'a>,
    pub(crate) subject_placements: &'a Placements,
    /// `-db`: the reported sequences (their parsed Seq-ids and OIDs, by rank); `None` for
    /// `-subject` records.
    pub(crate) database: Option<ReportDatabase<'a>>,
    pub(crate) show_gis: bool,
    /// The masks of each query (`Mask`), by query index.
    pub(crate) query_masks: Vec<Vec<Mask>>,
    pub(crate) query_code: DisplayCode,
    pub(crate) db_code: DisplayCode,
}

impl ReportSequences for AnchoredSequences<'_> {
    /// A query's own record as a subject is refused before the reports
    /// (`report::refuse_query_records_in_sections`).
    fn subject_is_query(&self, _q_idx: usize, _s_idx: u32) -> bool {
        false
    }

    fn query_letters(&self, q_idx: usize, from: u32, to: u32, out: &mut Vec<u8>) {
        query_record_letters(
            &self.queries[q_idx],
            self.query_placements.offset(q_idx),
            from,
            to,
            out,
        );
    }

    fn subject_letters(&self, s_idx: u32, from: u32, to: u32, out: &mut Vec<u8>) {
        let s_idx = s_idx as usize;
        let shift = self.subject_placements.offset(s_idx);
        let length = self.subject_letters.length(s_idx);
        let end = (to as usize + 1 - shift).min(length);
        let start = (from as usize - shift).min(end);
        self.subject_letters.letters_of(s_idx, start, end, out);
    }

    fn query_label(&self, q_idx: usize, out: &mut Vec<u8>) {
        record_label(&self.queries[q_idx], self.show_gis, out);
    }

    fn subject_label(&self, s_idx: u32, out: &mut Vec<u8>) {
        let rank = s_idx as usize;
        match self.database {
            Some(database) => database_label(
                database
                    .reported
                    .parsed
                    .get(rank)
                    .and_then(Option::as_deref),
                database.oids.oid(database.reported.ordinals[rank]),
                out,
            ),
            None => record_label(&self.subjects[rank], self.show_gis, out),
        }
    }

    fn query_masks(&self, q_idx: usize, out: &mut Vec<Mask>) {
        if let Some(masks) = self.query_masks.get(q_idx) {
            out.extend_from_slice(masks);
        }
    }

    fn subject_masks(&self, _q_idx: usize, _s_idx: u32, _out: &mut Vec<Mask>) {}

    fn query_genetic_code(&self) -> Option<DisplayCode> {
        Some(self.query_code)
    }

    fn subject_genetic_code(&self) -> Option<DisplayCode> {
        Some(self.db_code)
    }
}

/// The refusals of the query-anchored formats (outfmt 1 to 4, `view`) of a TBLASTX report call,
/// before its first byte (`run`, its final hits `hits` in report order with record coordinates and
/// their Seq-aligns' identities `seqalign_identities`, `-sorthits` and `-sorthsps` `sort`):
///
/// - A query whose own record is one of its subjects (`-subject` read with `-parse_deflines`,
///   `query_anchored::subject_is_query_record`), with alignments to show (`-num_alignments`):
///   TBLASTX's sections add their Dense-segs with `fForceTranslation` (not `fPreserveRows`), so
///   NCBI's `CAlnMix` takes the query's row and that subject's for one sequence and merges them
///   (`aln_mix`: not ported). The refusal takes every query that reports its own record, also one
///   whose record lies past the first `-num_alignments` subjects.
/// - `-db`: a reported sequence that NCBI's scope finds as a query of its batch, by its Seq-align
///   id (a generated `lcl|Query_N`) or, in a description line, by its display id
///   (`common::database::refuse_query_ids_of_anchored_reports`, as BLASTN's, TBLASTN's and
///   BLASTX's): NCBI shows the query in its place.
#[allow(clippy::too_many_arguments)]
pub(crate) fn refuse_anchored_reports(
    run: &super::blast_engine::TblastxReportRun<'_>,
    hits: &[super::blast_engine::TblastxHsp],
    seqalign_identities: Option<&[u32]>,
    sort: crate::report::hit_sort::ReportSort,
    num_alignments: usize,
    num_descriptions: usize,
) -> anyhow::Result<()> {
    let Some(database) = run.database else {
        if num_alignments == 0 {
            return Ok(());
        }
        if let Some(hsp) = hits.iter().find(|hsp| {
            match (
                run.queries.get(hsp.hit.q_idx as usize),
                run.subjects.get(hsp.hit.s_idx as usize),
            ) {
                (Some(query), Some(subject)) => {
                    crate::report::query_anchored::subject_is_query_record(query, subject)
                }
                _ => false,
            }
        }) {
            let mut fasta = Vec::new();
            if let Some(ids) = run.queries[hsp.hit.q_idx as usize].believed() {
                crate::seqid::label::fasta(&ids.best, &mut fasta);
            }
            anyhow::bail!(
                "the record {} is both a query and a subject under -parse_deflines, and in the \
                 query-anchored formats (outfmt 1 to 4) NCBI BLAST+ merges its query row and its \
                 subject rows of TBLASTX into one sequence; this is not supported by LOSAT's \
                 TBLASTX",
                String::from_utf8_lossy(&fasta)
            );
        }
        return Ok(());
    };
    let reported = database.reported;
    let scopes = crate::blastinput::query_batch::AnchoredBatchScopes::new(
        run.queries,
        run.archive_spans.to_vec(),
    );
    let query_hits = |q_idx: usize| {
        let from = hits.partition_point(|hsp| (hsp.hit.q_idx as usize) < q_idx);
        let to = hits.partition_point(|hsp| (hsp.hit.q_idx as usize) <= q_idx);
        from..to.max(from)
    };
    crate::algorithm::common::database::refuse_query_ids_of_anchored_reports(
        run.queries.len(),
        |q_idx| {
            let mut ranks: Vec<u32> = Vec::new();
            for hsp in &hits[query_hits(q_idx)] {
                if !ranks.contains(&hsp.hit.s_idx) {
                    ranks.push(hsp.hit.s_idx);
                }
            }
            ranks
        },
        |q_idx| {
            // The rows of the query's hits as the report sorts them (`QueryGroup::sorted`, with
            // TBLASTX's alignment lengths in nucleotides), each with its Seq-align's identities.
            let range = query_hits(q_idx);
            let rows: Vec<crate::report::pairwise::PairwiseHit> = hits[range.clone()]
                .iter()
                .zip(range)
                .map(|(hsp, index)| {
                    let mut row = crate::report::pairwise::PairwiseHit::from(hsp.hit.clone());
                    if let Some(identities) = seqalign_identities {
                        row.seqalign_num_ident = identities[index] as usize;
                    }
                    row
                })
                .collect();
            let indexed: Vec<_> = rows.iter().enumerate().collect();
            Ok(crate::report::x_pairwise_par::QueryGroup::sorted(
                &indexed,
                &[],
                |previous, s_idx| {
                    crate::report::db_subject::rank_joins_as_copy(
                        &reported.parsed,
                        previous as usize,
                        s_idx as usize,
                    )
                },
                sort,
                true,
            )
            .map_or_else(Vec::new, |group| group.described(num_descriptions).to_vec()))
        },
        &reported.parsed,
        |rank| database.oids.oid(reported.ordinals[rank as usize]),
        |q_idx| scopes.of_query(q_idx),
        scopes.generated(),
        "TBLASTX",
    )
}
