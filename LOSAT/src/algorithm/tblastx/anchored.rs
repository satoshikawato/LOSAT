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
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1547-1548,1579-1584
//! (`CBlastFormat::PrintOneResultSet`, formats 0 to 4)
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2100-2101
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1244-1251
//! (`x_FillSeqid`, the labels: showalign.cpp:4545-4572, cited at `query_anchored::fill_seqid`.)

use super::blast_engine::ReportDatabase;
use super::report::{DisplayCode, SubjectLetters};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::report::query_anchored::{
    database_label, query_record_letters, record_label, subject_is_query_record, Mask,
    ReportSequences,
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
    /// A `-subject` record read with `-parse_deflines` that is the query's record
    /// (`subject_is_query_record`, as BLASTN's): its HSPs with the query are one `CAlnMixSeq` of
    /// both rows (`query_anchored::aln_mix`, `fForceTranslation`). A database sequence with the
    /// query's id is refused before the reports (`refuse_pairwise_reports`).
    fn subject_is_query(&self, q_idx: usize, s_idx: u32) -> bool {
        self.database.is_none()
            && match (self.queries.get(q_idx), self.subjects.get(s_idx as usize)) {
                (Some(query), Some(subject)) => subject_is_query_record(query, subject),
                _ => false,
            }
    }

    fn query_letters(&self, q_idx: usize, from: u32, to: u32, out: &mut Vec<u8>) {
        query_record_letters(
            crate::report::row_letters::QueryLetters::new(self.queries[q_idx].seq(), None),
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

/// The refusals of a pairwise format of a TBLASTX report call (outfmt 0 to 4, `view`), before
/// any byte of it that is not NCBI's own prefix (the prolog, and with the stream of
/// `report_stream` the earlier batches) (`run`, its final hits `hits` in report order with record coordinates and their
/// Seq-aligns' identities `seqalign_identities`, `-sorthits` and `-sorthsps` `sort`): a
/// description line that shows what LOSAT does not write as NCBI
/// (`common::database::refuse_described_lines`, line RG-KB item 12 review 1; of the `-subject`
/// records whose ids `-parse_deflines` read, `records`, `refuse_described_records`, review 2),
/// and with `-db` in the
/// query-anchored formats (outfmt 1 to 4) a reported sequence that NCBI's scope finds as a query
/// of its batch, by its Seq-align id (a generated `lcl|Query_N`) or, in a description line, by its
/// display id (`common::database::refuse_query_ids_of_anchored_reports`, as BLASTN's, TBLASTN's
/// and BLASTX's): NCBI shows the query in its place. With the stream of `report_stream` a call is
/// one query batch (NCBI's scope is the batch's), so the refusal of a generated id comes after the
/// reports of the batches before; a description line that LOSAT refuses turns that stream off
/// (`common::database::reports_may_refuse`), as does, under `-parse_deflines`, a database sequence
/// with the id of a query of the run (`TblastxDatabase::reports_may_refuse`), whose refusal
/// (`common::database::refuse_reported_ids_of_queries`, in the batch loop) comes before the first
/// batch's reports, as outfmt 0's.
///
/// A query whose own record is one of its `-subject` records (`-parse_deflines`) has no refusal
/// here: its section is NCBI's (`query_anchored::aln_mix`), also where NCBI's merge throws, which
/// ends the report after the query's description table with NCBI's exception
/// (`query_anchored::write_section`, `cli::NativeError::of_writer`).
pub(crate) fn refuse_pairwise_reports(
    run: &super::blast_engine::TblastxReportRun<'_>,
    hits: &[super::blast_engine::TblastxHsp],
    seqalign_identities: Option<&[u32]>,
    view: crate::report::pairwise::PairwiseView,
    sort: crate::report::hit_sort::ReportSort,
    num_descriptions: usize,
    records: &[Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
) -> anyhow::Result<()> {
    if run.database.is_none() && records.is_empty() {
        return Ok(());
    }
    // The report joins no hits of `-subject` records (those with one id are refused before it).
    let joined: &[Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>] =
        match run.database {
            Some(database) => &database.reported.parsed,
            None => &[],
        };
    let query_hits = |q_idx: usize| {
        let from = hits.partition_point(|hsp| (hsp.hit.q_idx as usize) < q_idx);
        let to = hits.partition_point(|hsp| (hsp.hit.q_idx as usize) <= q_idx);
        from..to.max(from)
    };
    let reported_ranks = |q_idx: usize| {
        let mut ranks: Vec<u32> = Vec::new();
        for hsp in &hits[query_hits(q_idx)] {
            if !ranks.contains(&hsp.hit.s_idx) {
                ranks.push(hsp.hit.s_idx);
            }
        }
        ranks
    };
    let described = |q_idx: usize| {
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
                    joined,
                    previous as usize,
                    s_idx as usize,
                )
            },
            sort,
            true,
        )
        .map_or_else(Vec::new, |group| group.described(num_descriptions).to_vec()))
    };
    // TBLASTX's subjects are nucleotide sequences; the refusals name the running application.
    let program = run.args.app().name();
    let Some(database) = run.database else {
        return crate::algorithm::common::database::refuse_described_records(
            run.queries.len(),
            reported_ranks,
            described,
            records,
            false,
            program,
        );
    };
    let reported = database.reported;
    let oid_of = |rank: u32| database.oids.oid(reported.ordinals[rank as usize]);
    crate::algorithm::common::database::refuse_described_lines(
        run.queries.len(),
        reported_ranks,
        described,
        &reported.parsed,
        |rank| reported.loads.load(rank as usize),
        oid_of,
        false,
        program,
    )?;
    let scopes = crate::blastinput::query_batch::AnchoredBatchScopes::new(
        run.queries,
        run.archive_spans.to_vec(),
    );
    crate::algorithm::common::database::refuse_query_ids_of_anchored_reports(
        run.queries.len(),
        reported_ranks,
        described,
        &reported.parsed,
        oid_of,
        |q_idx| scopes.of_query(q_idx),
        scopes.generated(),
        view != crate::report::pairwise::PairwiseView::Pairwise,
        program,
    )
}
