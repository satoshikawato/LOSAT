//! The BLASTN pairwise report (outfmt 0): the displayed alignment of every HSP of the final
//! hit list, and the per-query statistics that the report prints.

use std::iter::repeat_n;

use crate::blastinput::fasta_reader::FastaRecord;
use anyhow::{ensure, Context, Result};

use crate::common::GapEditOp;
#[cfg(test)]
use crate::core::blast_encoding::encode_iupac_to_blastna;
use crate::report::pairwise::PairwiseHit;
use crate::stats::{
    compute_karlin_params_ungapped, score_freq_profile_from_probabilities, KarlinParams,
};
use crate::utils::dust::MaskedInterval;

use super::alignment::build_blastna_matrix;
use super::hsp::{BlastnHitList, BlastnHsp, BlastnHspList};
use super::lookup::reverse_complement_in_place;
use super::x_subject_store::SubjectLetters;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1520-1534
// ```c++
// void RemapToQueryLoc(CRef<CSeq_align> sar, const CSeq_loc & query)
// {
//     _ASSERT(sar);
//     const int query_row = 0;
//
//     TSeqPos q_shift = 0;
//
//     if (query.IsInt()) {
//         q_shift = query.GetInt().GetFrom();
//     }
//
//     if (q_shift > 0) {
//         sar->OffsetRow(query_row, q_shift);
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:108-124
// ```c++
// CSeqVector& CAlnVec::x_GetSeqVector(TNumrow row) const
// {
//     TSeqVectorCache::iterator iter = m_SeqVectorCache.find(row);
//     CRef<CSeqVector> seq_vec;
//     if (iter != m_SeqVectorCache.end()) {
//         seq_vec = iter->second;
//     }
//     else {
//         CBioseq_Handle h = GetBioseqHandle(row);
//         CSeqVector vec = h.GetSeqVector
//             (CBioseq_Handle::eCoding_Iupac,
//              IsPositiveStrand(row) ?
//              CBioseq_Handle::eStrand_Plus :
//              CBioseq_Handle::eStrand_Minus);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:458-460
// ```c++
//     // set sequence range
//     retval->SetInt().SetFrom(from);
//     retval->SetInt().SetTo((to > 0 && to < seqlen) ? to : (seqlen-1));
// ```
/// The letters of a query's record as NCBI's reports read them (the `CSeqVector` of the query's
/// Bioseq: the whole record, also where `-query_loc` searched an interval of it) and the record
/// position of the searched letters' first (the interval's start, by which `RemapToQueryLoc`
/// moves the query row of every Seq-align; 0 without `-query_loc`).
#[derive(Clone, Copy)]
pub(crate) struct QueryRecord<'a> {
    pub letters: &'a [u8],
    pub offset: usize,
}

/// The query records of a report call as its rows read them (`QueryRecord`): the searched records,
/// and with `-query_loc` the input records, whole, whose intervals they are.
#[derive(Clone, Copy)]
pub(crate) struct QueryRecords<'a> {
    /// The searched records (with `-query_loc`, the intervals' letters).
    searched: &'a [FastaRecord],
    /// Where each searched query lies in its record, and its input position.
    query_input: &'a crate::blastinput::seq_range::QueryInput,
    /// The input records from input position `first` on (read with `-query_loc` only).
    input: &'a [FastaRecord],
    first: usize,
}

impl<'a> QueryRecords<'a> {
    /// The records of the searched queries `searched` placed by `query_input`, of which `input`
    /// holds the input records from input position `first` on.
    pub(crate) fn new(
        searched: &'a [FastaRecord],
        query_input: &'a crate::blastinput::seq_range::QueryInput,
        input: &'a [FastaRecord],
        first: usize,
    ) -> Self {
        Self {
            searched,
            query_input,
            input,
            first,
        }
    }

    /// The searched records.
    pub(crate) fn searched(&self) -> &'a [FastaRecord] {
        self.searched
    }

    /// The record of searched query `q_idx` (`None` for a query the call does not hold).
    pub(crate) fn record(&self, q_idx: usize) -> Option<QueryRecord<'a>> {
        let searched = self.searched.get(q_idx)?;
        if !self.query_input.placements.ranged() {
            return Some(QueryRecord {
                letters: searched.seq(),
                offset: 0,
            });
        }
        let record = self
            .input
            .get(self.query_input.ordinal(q_idx).checked_sub(self.first)?)?;
        Some(QueryRecord {
            letters: record.seq(),
            offset: self.query_input.placements.offset(q_idx),
        })
    }
}

/// The regions that the report prints in lowercase.
pub(crate) struct DisplayMasks<'a> {
    /// Per query, on its plus strand of the searched letters: DUST and, with
    /// `-lcase_masking`, the input lowercase.
    pub query: &'a [Vec<MaskedInterval>],
    /// Where each query's searched letters lie in its record (`-query_loc`).
    pub query_offsets: &'a crate::blastinput::seq_range::Placements,
    /// Per subject, with `-lcase_masking`: the input lowercase of the subject's record.
    pub subject: Option<&'a [Vec<MaskedInterval>]>,
    /// Per reported database sequence (by rank), with `-db_soft_mask`: its stored masks of the
    /// search's algorithm in the order they are stored (`BlastnDatabase::display_masks`).
    pub database: Option<&'a [Vec<MaskedInterval>]>,
}

impl<'a> Rows<'a> {
    /// The rows' inputs of a report call of `queries` and `subjects`, its shown masks made
    /// once for every window of its queries (`pairwise_hits`).
    // Kept out of line: this runs only for outfmt 0 or hit records, and inlining it into the
    // shared post-processing makes every BLASTN run compile it in Wasm hosts.
    #[inline(never)]
    pub(crate) fn new(
        queries: QueryRecords<'a>,
        subjects: SubjectLetters<'a>,
        masks: &DisplayMasks<'_>,
    ) -> Self {
        let query_masks: Vec<Vec<MaskedInterval>> = masks
            .query
            .iter()
            .zip(queries.searched())
            .enumerate()
            .map(|(q_idx, (query_masks, query))| {
                shown_query_masks(
                    query_masks,
                    query.seq().len(),
                    masks.query_offsets.offset(q_idx),
                )
            })
            .collect();
        // A subject set's lowercase is merged first; a database's masks are kept as stored
        // (`shown_subject_masks` merges what it keeps).
        let subject_masks: Option<Vec<Vec<MaskedInterval>>> = match (masks.database, masks.subject)
        {
            (Some(database), _) => Some(database.to_vec()),
            (None, Some(subject)) => Some(subject.iter().map(|masks| merged(masks)).collect()),
            (None, None) => None,
        };
        Rows {
            queries,
            subjects,
            query_masks,
            subject_masks,
        }
    }
}

/// Every HSP of the final hit lists `hit_lists` (those of all queries of `rows`, or of a window
/// of them), in the order of the tabular formats (queries, then the subjects and HSPs of each
/// query as ranked), with its displayed rows.
// Kept out of line: this runs only for outfmt 0 or hit records, and inlining it into the
// shared post-processing makes every BLASTN run compile it in Wasm hosts.
#[inline(never)]
pub(crate) fn pairwise_hits(
    hit_lists: &[Option<BlastnHitList>],
    rows: &Rows<'_>,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
) -> Result<Vec<PairwiseHit>> {
    let hsp_lists = hit_lists
        .iter()
        .flatten()
        .flat_map(|hit_list| &hit_list.hsplist_array);
    // EXPERIMENT (LOSAT_X_PAIRWISEPAR): each subject's list of each query on the pool, the rows
    // collected in the list order (`report::x_pairwise_par`); the same function of each list.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if let Some(pool) = x_pool {
        use rayon::prelude::*;
        let lists: Vec<&BlastnHspList> = hsp_lists.collect();
        let rows_of_lists: Vec<Result<Vec<PairwiseHit>>> = pool.install(|| {
            lists
                .par_iter()
                .map(|hsp_list| {
                    let mut hits = Vec::with_capacity(hsp_list.hsps.len());
                    rows.list(hsp_list, &mut hits)?;
                    Ok(hits)
                })
                .collect()
        });
        let mut hits = Vec::with_capacity(lists.iter().map(|list| list.hsps.len()).sum());
        for list_rows in rows_of_lists {
            hits.append(&mut list_rows?);
        }
        return Ok(hits);
    }
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    let _ = x_pool;
    let mut hits = Vec::new();
    for hsp_list in hsp_lists {
        rows.list(hsp_list, &mut hits)?;
    }
    Ok(hits)
}

/// The windows of the queries `0..num_queries` whose outfmt 0 rows are made and written
/// together (`x_pairwise_par::query_windows` of the estimates of their HSPs' lengths in a report
/// of `line_length` columns a line, the report's `PairwiseConfig::line_length`, before any row is
/// made, with the budget of `x_pool`, `x_pairwise_par::window_bytes`).
#[inline(never)]
pub(crate) fn pairwise_windows(
    hit_lists: &[Option<BlastnHitList>],
    num_queries: usize,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
    line_length: usize,
) -> Vec<std::ops::Range<usize>> {
    use crate::report::x_pairwise_par::{query_estimate_of, query_windows, window_bytes};
    let estimates: Vec<usize> = (0..num_queries)
        .map(|q_idx| {
            query_estimate_of(
                hit_lists
                    .get(q_idx)
                    .and_then(Option::as_ref)
                    .into_iter()
                    .flat_map(|hit_list| &hit_list.hsplist_array)
                    .flat_map(|hsp_list| &hsp_list.hsps)
                    .map(|hsp| hsp.length),
                line_length,
            )
        })
        .collect();
    query_windows(&estimates, window_bytes(x_pool))
}

/// The displayed rows of one HSP (the query on its plus strand, the subject reverse-complemented
/// with it for a minus-strand HSP, upper case, no masking) for BLAST XML, which recounts its
/// identities from them (`algorithm::blastn::blast_xml`). The same function of the HSP as the
/// outfmt 0 rows, without a `PairwiseHit`, written into two buffers that the caller reuses.
pub(crate) fn xml_rows_into(
    hsp: &BlastnHsp,
    queries: QueryRecords<'_>,
    subjects: SubjectLetters<'_>,
    query_row: &mut Vec<u8>,
    subject_row: &mut Vec<u8>,
) -> Result<()> {
    let query = queries
        .record(hsp.q_idx as usize)
        .context("BLASTN HSP of an unknown query")?;
    let s_idx = hsp.s_idx as usize;
    subjects
        .records()
        .get(s_idx)
        .context("BLASTN HSP of an unknown subject")?;
    displayed_rows_into(
        hsp,
        query,
        subjects.len(s_idx),
        |from, to| subjects.upper_letters(s_idx, from, to),
        query_row,
        subject_row,
    )
}

/// What the displayed rows of an HSP read besides the HSP (`pairwise_hits`).
pub(crate) struct Rows<'a> {
    queries: QueryRecords<'a>,
    subjects: SubjectLetters<'a>,
    /// Per query, the shown masks (`shown_query_masks`).
    query_masks: Vec<Vec<MaskedInterval>>,
    /// Per subject, the masks that `shown_subject_masks` shows with its HSPs: of a subject set's
    /// record with `-lcase_masking`, its merged lowercase (`DisplayMasks::subject`); of a reported
    /// database sequence with `-db_soft_mask`, its stored masks as stored
    /// (`DisplayMasks::database`).
    subject_masks: Option<Vec<Vec<MaskedInterval>>>,
}

impl Rows<'_> {
    /// The shown masks of query `q_idx` (`shown_query_masks`), the query's `TMaskedQueryRegions`
    /// in record positions (the query-anchored reports, `algorithm::blastn::anchored`).
    pub(crate) fn query_masks(&self, q_idx: usize) -> &[MaskedInterval] {
        self.query_masks.get(q_idx).map_or(&[], Vec::as_slice)
    }

    /// The masks of the subject of the HSPs `hsps` (one subject's list of a query) that NCBI's
    /// `GetSubjectMasks` gives with them (`subject_mask_regions`), in their order and unmerged:
    /// the query-anchored reports place each one on the first row that meets it
    /// (`query_anchored::fill_loc_list`). None without subject masks.
    pub(crate) fn subject_masks(&self, hsps: &[BlastnHsp]) -> Vec<MaskedInterval> {
        match (self.subject_masks.as_deref(), hsps.first()) {
            (Some(subject_masks), Some(hsp)) => {
                subject_mask_regions(&subject_masks[hsp.s_idx as usize], hsps).collect()
            }
            _ => Vec::new(),
        }
    }

    /// The rows of every HSP of one subject's list, appended to `hits` in the list's order.
    fn list(&self, hsp_list: &BlastnHspList, hits: &mut Vec<PairwiseHit>) -> Result<()> {
        let (queries, subjects) = (self.queries, self.subjects);
        let list_subject_masks = self.subject_masks.as_deref().and_then(|subject_masks| {
            let s_idx = hsp_list.hsps.first()?.s_idx as usize;
            Some(shown_subject_masks(&subject_masks[s_idx], &hsp_list.hsps))
        });
        for hsp in &hsp_list.hsps {
            let query = queries
                .record(hsp.q_idx as usize)
                .context("BLASTN HSP of an unknown query")?;
            let s_idx = hsp.s_idx as usize;
            let subject = subjects
                .records()
                .get(s_idx)
                .context("BLASTN HSP of an unknown subject")?;
            let subject_length = subjects.len(s_idx);
            let (mut query_row, mut subject_row) =
                displayed_rows(hsp, query, subject_length, |from, to| {
                    subjects.upper_letters(s_idx, from, to)
                })?;
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1796-1798
            // ```c++
            //     x_FillIdentityInfo(aln_vec_info->alnRowInfo->sequence[0],
            //                        aln_vec_info->alnRowInfo->sequence[1],
            //                        aln_vec_info->match,
            // ```
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2127-2140
            // ```c++
            //     int min_length=min<int>((int)sequence_standard.size(), (int)sequence.size());
            //     ...
            //     for(int i=0; i<min_length; i++){
            //         if(sequence_standard[i]==sequence[i]){
            //     ...
            //             match ++;
            // ```
            // The identities that outfmt 0 shows are those of the displayed rows (the stored
            // letters in upper case, before the masks' lower case), counted again (series NP,
            // M4a-2b decision H6): with `-db_hard_mask` the engine's `num_ident` counts a masked
            // letter ('B' in the traceback) as a mismatch. Elsewhere the two counts are equal.
            let displayed_ident = query_row
                .iter()
                .zip(&subject_row)
                .filter(|(query, subject)| query == subject)
                .count();
            lowercase_masked(
                &mut query_row,
                hsp.q_start,
                1,
                &self.query_masks[hsp.q_idx as usize],
            );
            if let Some(subject_masks) = &list_subject_masks {
                let step = if hsp.s_start > hsp.s_end { -1 } else { 1 };
                lowercase_masked(&mut subject_row, hsp.s_start, step, subject_masks);
            }
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2601-2611
            // ```c++
            //        for (int row=0; row<m_AV->GetNumRows(); row++) {
            //            CRef<CAlnMap::CAlnChunkVec> chunk_vec
            //                = m_AV->GetAlnChunks(row, m_AV->GetSeqAlnRange(0));
            //            for (int i=0; i<chunk_vec->size(); i++) {
            //                CConstRef<CAlnMap::CAlnChunk> chunk = (*chunk_vec)[i];
            //                if (chunk->IsGap()) {
            //                    gap += (chunk->GetAlnRange().GetTo()
            //                            - chunk->GetAlnRange().GetFrom() + 1);
            // ```
            // The gaps are the columns of the edit script's gap operations, not the gap
            // characters of the rows: a query letter past the query's end shows one too
            // (`context_letters`).
            let gaps = hsp.gap_info.as_deref().map_or(0, |ops| {
                ops.iter()
                    .filter(|op| !matches!(op, GapEditOp::Sub(_)))
                    .map(|op| op.num() as usize)
                    .sum()
            });
            let mut hit = hsp.clone().into_hit();
            hit.num_ident = displayed_ident;
            hits.push(PairwiseHit {
                hit,
                query_seq: Some(String::from_utf8(query_row).context("BLASTN query row")?),
                subject_seq: Some(String::from_utf8(subject_row).context("BLASTN subject row")?),
                query_frame: None,
                subject_frame: None,
                positives: None,
                gaps: Some(gaps),
                subject_length: Some(subject_length),
                // The title after its first word; the reports take the subject's title
                // bytes by subject index.
                subject_title: subject
                    .title
                    .iter()
                    .position(|&byte| byte == b' ')
                    .map(|space| String::from_utf8_lossy(&subject.title[space + 1..]).into_owned()),
                comp_adjust_method: None,
                sum_n: None,
                subject_stored_gaps: Vec::new(),
                // The Seq-align's `num_ident` (`seqalign::query_records` writes the same
                // `hsp.num_ident`): the engine's count, which "Max ident" and the sorts read.
                seqalign_num_ident: hsp.num_ident,
            });
        }
        Ok(())
    }
}

// The rows of one HSP as NCBI shows them: the query on its plus strand.
//
// NCBI reference: c++/src/objtools/align_format/showalign.cpp:1852-1858
// ```c++
//         if((ds->IsSetStrands()
//             && ds->GetStrands().front()==eNa_strand_minus)
//            && !(ds->IsSetWidths() && ds->GetWidths()[0] == 3)){
//             //show plus strand if master is minus for non-translated case
//             finalDenseg->Reverse();
// ```
// NCBI reference: c++/include/algo/blast/core/gapinfo.h:43-61
// ```c
// typedef enum EGapAlignOpType {
//    eGapAlignDel = 0, /**< Deletion: a gap in query */
//    ...
//    eGapAlignSub = 3, /**< Substitution */
//    ...
//    eGapAlignIns = 6, /**< Insertion: a gap in subject */
// ```
// The HSP offsets and its edit script are those of the query context (the reverse
// complement of the query for a minus-strand HSP) against the plus-strand subject.
// Reversing the alignment of a minus-strand HSP reverse-complements both rows. The
// residues are the input letters in uppercase (NCBI's IUPAC sequence data); masking is
// applied afterwards. Only the aligned segments are read: the segment of the minus-strand
// context `[qs, qe)` is the reverse complement of query `[len - qe, len - qs)`, read in the
// query's record (`context_letters`). The subject has `subject_length` letters, and
// `subject_letters(from, to)` gives its letters `[from, to)` in uppercase.
fn displayed_rows(
    hsp: &BlastnHsp,
    query: QueryRecord<'_>,
    subject_length: usize,
    subject_letters: impl FnOnce(usize, usize) -> Vec<u8>,
) -> Result<(Vec<u8>, Vec<u8>)> {
    let (mut query_row, mut subject_row) = (Vec::new(), Vec::new());
    displayed_rows_into(
        hsp,
        query,
        subject_length,
        subject_letters,
        &mut query_row,
        &mut subject_row,
    )?;
    Ok((query_row, subject_row))
}

/// `displayed_rows` into `query_row` and `subject_row` (cleared first; their capacity is kept),
/// reversed in place for a minus-strand HSP.
fn displayed_rows_into(
    hsp: &BlastnHsp,
    query: QueryRecord<'_>,
    subject_length: usize,
    subject_letters: impl FnOnce(usize, usize) -> Vec<u8>,
    query_row: &mut Vec<u8>,
    subject_row: &mut Vec<u8>,
) -> Result<()> {
    context_rows_into(
        hsp,
        query,
        subject_length,
        subject_letters,
        query_row,
        subject_row,
    )?;
    if hsp.query_frame < 0 {
        reverse_complement_in_place(query_row);
        reverse_complement_in_place(subject_row);
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:279-300
// ```c++
//             if (start >= 0) {
//                 stop = start + len - 1;
//
//                 // add regular sequence to buffer
//                 GetSeqString(buff, row, start, stop);
//                 TSeqPos buf_len = min<TSeqPos>(TSeqPos(buff.size()), seg_len);
//                 buffer += buff;
//                 if (buf_len < seg_len) {
//                     // Not enough chars in the sequence, add gap
//                     ...
//                     } else {
//                         fill_ch = GetGapChar(row);
//                     }
// ```
// NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnvec.hpp:304-305
// ```c++
//         if (IsPositiveStrand(row)) {
//             x_GetSeqVector(row).GetSeqData(seq_from, seq_to + 1, buffer);
// ```
// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:631-643, 720-728
// ```c++
// void CSeqVector_CI::x_SetPos(TSeqPos pos)
// {
//     TSeqPos size = x_GetSize();
//     if ( pos >= size ) {
//         ...
//         m_CachePos = size;
//         ...
//         return;
// ...
//     TSeqPos pos = GetPos();
//     _ASSERT(pos <= x_GetSize());
//     count = min(count, x_GetSize() - pos);
//     if ( !count ) {
//         return;
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1852-1857
// ```c++
//     if((ds->IsSetStrands()
//         && ds->GetStrands().front()==eNa_strand_minus)
//        && !(ds->IsSetWidths() && ds->GetWidths()[0] == 3)){
//         //show plus strand if master is minus for non-translated case
//         finalDenseg->Assign(*ds);
//         finalDenseg->Reverse();
// ```
/// The query letters of an HSP's context offsets `[query.offset, query.end)`, in the context's
/// direction (the complement for a minus-strand HSP, still to be reversed with the rows), as the
/// reports read them through the Seq-align (`CAlnVec::GetWholeAlnSeqString` on the plus strand of
/// the query's record, after `Reverse` for a minus-strand HSP): context offset `c` lies at record
/// position `offset + c` on the plus strand and `offset + length - 1 - c` on the minus strand
/// (`QueryRecord::offset`, `RemapToQueryLoc`; `length` the context's), and a position past the
/// record's end is not read (`CSeqVector` stops at its size): its column shows the gap character.
/// Only an HSP that NCBI keeps outside its query context (remaining-gaps X1b-o) reaches past the
/// searched letters: with `-query_loc` it reads the record's letters outside the interval, and
/// past the record's end every column is a gap character (X1b-p). A position below 0 never
/// reaches the rows (`run::left_of_query_refusal`); it would be a gap character too.
fn context_letters(hsp: &BlastnHsp, query: QueryRecord<'_>) -> Vec<u8> {
    let (qs, qe) = (
        i64::from(hsp.internal_q_offset_0),
        i64::from(hsp.internal_q_end_0),
    );
    let (shift, length) = (query.offset as i64, hsp.query_length as i64);
    let record_length = query.letters.len() as i64;
    let minus = hsp.query_frame < 0;
    // The record positions `[from, to)` of the context offsets `[qs, qe)` on the plus strand.
    let (from, to) = if minus {
        (shift + length - qe, shift + length - qs)
    } else {
        (shift + qs, shift + qe)
    };
    // The positions below 0 and those from the record's end on are not read; in the context's
    // direction the plus strand's low positions come first, the minus strand's high ones.
    let below = (to.min(0) - from).max(0) as usize;
    let past = (to - from.max(record_length)).max(0) as usize;
    let (first, last) = if minus { (past, below) } else { (below, past) };
    let (read_from, read_to) = (from.max(0), to.min(record_length));
    let mut context = Vec::with_capacity((qe - qs).max(0) as usize);
    context.resize(first, b'-');
    if read_from < read_to {
        let letters = &query.letters[read_from as usize..read_to as usize];
        if minus {
            let at = context.len();
            context.extend_from_slice(letters);
            reverse_complement_in_place(&mut context[at..]);
        } else {
            context.extend(letters.iter().map(u8::to_ascii_uppercase));
        }
    }
    context.resize(context.len() + last, b'-');
    context
}

/// The two rows in the direction of the HSP's query context (cleared first; their capacity is
/// kept): a minus-strand HSP's rows are still to be reverse-complemented.
fn context_rows_into(
    hsp: &BlastnHsp,
    query: QueryRecord<'_>,
    subject_length: usize,
    subject_letters: impl FnOnce(usize, usize) -> Vec<u8>,
    query_row: &mut Vec<u8>,
    subject_row: &mut Vec<u8>,
) -> Result<()> {
    query_row.clear();
    subject_row.clear();
    let (qs, qe) = (hsp.internal_q_offset_0, hsp.internal_q_end_0);
    let (ss, se) = (hsp.internal_s_offset_0, hsp.internal_s_end_0);
    // An HSP whose Seq-align has a query start below 0 in the record never reaches the rows
    // (`run::left_of_query_refusal`, `BlastnHsp::record_query_start`).
    ensure!(
        qs <= qe
            && query.offset as i64 + hsp.seqalign_query_start() >= 0
            && ss <= se
            && se <= subject_length,
        "BLASTN HSP outside its sequences"
    );
    let context = context_letters(hsp, query);
    let subject = subject_letters(ss, se);
    let ungapped = [GapEditOp::Sub((qe - qs) as u32)];
    let ops = hsp.gap_info.as_deref().unwrap_or(&ungapped);
    // The rows live in a `PairwiseHit` while the reports are written: their exact size, not
    // what growing leaves (a third more on the alignments of the BLASTN hit records). A
    // reused buffer keeps its capacity.
    let columns = ops.iter().map(|op| op.num() as usize).sum();
    query_row.reserve_exact(columns);
    subject_row.reserve_exact(columns);
    // Offsets into the two segments.
    let (mut q, mut s) = (0, 0);
    for op in ops {
        let n = op.num() as usize;
        let (take_query, take_subject) = match op {
            GapEditOp::Sub(_) => (true, true),
            GapEditOp::Del(_) => (false, true),
            GapEditOp::Ins(_) => (true, false),
        };
        if take_query {
            query_row.extend_from_slice(
                context
                    .get(q..q + n)
                    .context("BLASTN edit script past the query")?,
            );
            q += n;
        } else {
            query_row.extend(repeat_n(b'-', n));
        }
        if take_subject {
            subject_row.extend_from_slice(
                subject
                    .get(s..s + n)
                    .context("BLASTN edit script past the subject")?,
            );
            s += n;
        } else {
            subject_row.extend(repeat_n(b'-', n));
        }
    }
    ensure!(
        q == context.len() && s == subject.len(),
        "BLASTN edit script does not span its HSP"
    );
    Ok(())
}

// NCBI reference: c++/src/algo/blast/api/blast_aux.cpp:878-892
// ```c++
//         const TSeqRange kTarget((*query_interval)->GetFrom(),
//                                 (*query_interval)->GetTo());
//         ...
//             TSeqRange range(Map(kTarget, masked_range));
//             if (range.NotEmpty() && range != kTarget) {
//                 ...
//                 CRef<CSeqLocInfo> seqlocinfo
//                     (new CSeqLocInfo(seqint, CSeqLocInfo::eFrameNotSet));
// ```
// NCBI reference: c++/src/algo/blast/api/blast_aux.cpp:825-842
// ```c++
// template <class Position>
// CRange<Position> Map(const CRange<Position>& target,
//                      const CRange<Position>& range)
// {
//     if (target.Empty()) {
//         throw std::runtime_error("Target range is empty");
//     }
//
//     if (range.Empty() ||
//         (range.GetFrom() > target.GetTo()) ||
//         ((range.GetFrom() + target.GetFrom()) > target.GetTo())) {
//         return target;
//     }
//
//     CRange<Position> retval;
//     retval.SetFrom(max(target.GetFrom() + range.GetFrom(), target.GetFrom()));
//     retval.SetTo(min(target.GetFrom() + range.GetTo(), target.GetTo()));
//     return retval;
// }
// ```
// A mask that covers the whole query (its searched interval, `kTarget`) is not reported,
// so it is not shown. The masks of the search lie in the searched letters (`query_length`
// of them); `Map` moves them by the start of the query's interval (`offset`, 0 without
// `-query_loc`) into record coordinates, as the report's HSPs.
pub(crate) fn shown_query_masks(
    masks: &[MaskedInterval],
    query_length: usize,
    offset: usize,
) -> Vec<MaskedInterval> {
    merged(
        &masks
            .iter()
            .filter(|mask| !(mask.start == 0 && mask.end >= query_length))
            .map(|mask| MaskedInterval::new(mask.start + offset, mask.end + offset))
            .collect::<Vec<_>>(),
    )
}

// NCBI reference: c++/src/algo/blast/api/blast_seqalign.cpp:1588-1600
// ```c++
//         // Union subject sequence ranges
//         vector <TSeqRange> ranges;
//         for (int i=0; i<hsp_list->hspcnt; i++) {
//             const BlastHSP* hsp = hsp_list->hsp_array[i];
//             TSeqRange rg;
//             rg.SetFrom(hsp->subject.offset);
//             rg.SetTo(hsp->subject.end);
//             ranges.push_back(rg);
//         }
//
//         // Extract subject masks
//         TMaskedSubjRegions masks;
//         if (!ranges.empty() && seqinfo_src->GetMasks(kOid, ranges, masks)) {
// ```
// NCBI reference: c++/src/algo/blast/api/seqinfosrc_seqvec.cpp:109-125
// ```c++
// static void
// s_SeqIntervalToSeqLocInfo(CRef<CSeq_interval> interval,
//                           const vector <TSeqRange>& target_ranges,
//                           const CSeqLocInfo::ETranslationFrame frame,
//                           TMaskedSubjRegions& retval)
// {
//     TSeqRange loc(interval->GetFrom(), interval->GetTo());
//
//     for (size_t ir=0; ir< target_ranges.size(); ir++) {
//         if (target_ranges[ir] != TSeqRange::GetEmpty() &&
//            loc.IntersectingWith(target_ranges[ir])) {
//            CRef<CSeqLocInfo> sli(new CSeqLocInfo(interval, frame));
//            retval.push_back(sli);
//            return;
//         }
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqinfosrc_seqdb.cpp:121-133
// ```c++
//     CSeqDB::TSequenceRanges ranges;
//     m_iSeqDb->GetMaskData(index, m_FilteringAlgoId, ranges);
//     ITERATE(CSeqDB::TSequenceRanges, itr, ranges) {
//         for (size_t it=0; it<target.size(); it++) {
//             if (target[it] != TSeqRange::GetEmpty() &&
//                 target[it].IntersectingWith(*itr)) {
//                 CRef<CSeq_interval> si
//                     (new CSeq_interval(*id, itr->first, itr->second-1));
//                 CRef<CSeqLocInfo> sli(new CSeqLocInfo(si, kFrame));
//                 retval.push_back(sli);
//                 break;
//             }
//         }
//     }
// ```
// The subject masks shown with the HSPs of one subject (`GetMasks` of the subject set's
// `CSeqVecSeqInfoSrc` or of the database's `CSeqDbSeqInfoSrc`): each mask `[start, end)` that
// holds a letter (a stored database mask `TSeqRange(first, second-1)` of none is empty) and meets
// the closed range `[offset, end]` of an HSP in the searched letters (the end is one past the
// HSP), in the order of `masks`, merged for `lowercase_masked`, which marks each position once
// (the merged lowercase of a subject set's record is kept as it is). With `-subject_loc` the
// record's coordinates and the searched letters' differ by the interval's start, so a region of
// the record is kept by where it would lie in the interval, and shown where it lies in the
// record.
fn shown_subject_masks(masks: &[MaskedInterval], hsps: &[BlastnHsp]) -> Vec<MaskedInterval> {
    merged(&subject_mask_regions(masks, hsps).collect::<Vec<_>>())
}

/// The masks of `masks` that `GetSubjectMasks` gives with the HSPs `hsps` (`shown_subject_masks`
/// before the merge), in their order.
fn subject_mask_regions<'m>(
    masks: &'m [MaskedInterval],
    hsps: &'m [BlastnHsp],
) -> impl Iterator<Item = MaskedInterval> + 'm {
    masks
        .iter()
        .filter(|mask| {
            mask.start < mask.end
                && hsps.iter().any(|hsp| {
                    mask.start <= hsp.internal_s_end_0 && hsp.internal_s_offset_0 < mask.end
                })
        })
        .cloned()
}

/// The masked positions as disjoint intervals in order, for `lowercase_masked`.
fn merged(masks: &[MaskedInterval]) -> Vec<MaskedInterval> {
    let mut sorted = masks.to_vec();
    sorted.sort_by_key(|mask| mask.start);
    let mut merged: Vec<MaskedInterval> = Vec::with_capacity(sorted.len());
    for mask in sorted {
        match merged.last_mut() {
            Some(last) if mask.start <= last.end => last.end = last.end.max(mask.end),
            _ => merged.push(mask),
        }
    }
    merged
}

// NCBI reference: c++/src/objtools/align_format/showalign.cpp:2495-2521
// ```c++
//         if(id.Which() != CSeq_id::e_not_set){
//             /*only do this for sequence but not for others like middle line,
//               features*/
//             ...
//                 } else if (m_SeqLocChar==eLowerCase){
//                     actualSeq[i-start]=tolower((unsigned char) actualSeq[i-start]);
// ```
// The row's residue positions start at `first` (1-based) and move by `step`; a residue
// whose position lies in a mask is shown in lowercase. `masks` are disjoint and in order
// (`merged`), so the masks that overlap the row's span are found by binary search.
fn lowercase_masked(row: &mut [u8], first: usize, step: isize, masks: &[MaskedInterval]) {
    let residues = row.iter().filter(|&&residue| residue != b'-').count();
    if masks.is_empty() || residues == 0 {
        return;
    }
    let last = first as isize + step * (residues as isize - 1);
    let (low, high) = (
        (first as isize).min(last) - 1,
        (first as isize).max(last) - 1,
    );
    let (low, high) = (low as usize, high as usize);
    let candidates = &masks[masks.partition_point(|mask| mask.end <= low)
        ..masks.partition_point(|mask| mask.start <= high)];
    if candidates.is_empty() {
        return;
    }
    let mut position = first as isize;
    for residue in row.iter_mut().filter(|residue| **residue != b'-') {
        let zero_based = (position - 1) as usize;
        if candidates
            .iter()
            .any(|mask| mask.start <= zero_based && zero_based < mask.end)
        {
            residue.make_ascii_lowercase();
        }
        position += step;
    }
}

/// The ungapped Karlin-Altschul block of a blastn query, from its composition; `None`
/// when NCBI cannot compute it and marks the query invalid.
///
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:2778-2792
/// ```c
///       Blast_ResFreqString(sbp, rfp, (char*)buffer, query_length);
///       sbp->sfp[context] = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
///       BlastScoreFreqCalc(sbp, sbp->sfp[context], rfp, stdrfp);
///       sbp->kbp_std[context] = kbp = Blast_KarlinBlkNew();
///       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
///       if (loop_status) {
///           contexts[context].is_valid = FALSE;
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:2004-2013
/// ```c
///    for (lp = str, lpmax = lp+length; lp < lpmax; lp++)
///    {
///       ++rcp->comp[(int)(*lp & mask)];
///    }
///
///    /* Don't count ambig. residues. */
///    for (index=0; index<sbp->ambig_occupy; index++)
///    {
///       rcp->comp[sbp->ambiguous_res[index]] = 0;
///    }
/// ```
/// NCBI reference: c++/src/algo/blast/core/blast_setup.c:346-350
/// ```c
///     if (program_number == eBlastTypeBlastn ||
///         program_number == eBlastTypeMapping) {
///
///         BLAST_ScoreSetAmbigRes(sbp, 'N');
///         BLAST_ScoreSetAmbigRes(sbp, '-');
/// ```
/// The residues are counted in BLASTNA without `N` and the gap, normalized, and paired
/// with the standard composition (A, C, G, T at 0.25; `Blast_ResFreqStdComp`) through the
/// BLASTNA matrix (`BlastScoreFreqCalc`, blast_stat.c:2151-2205). The minus strand of a
/// query (its reverse complement) has the same residues in another order, so its block
/// can differ in the last bits; the search computes one per strand (`scoring.rs`). A
/// query of `N` only has no counted residue, so every score probability is 0 and the
/// calculation fails.
#[cfg(test)]
pub(crate) fn query_ungapped_karlin(
    query: &[u8],
    reward: i32,
    penalty: i32,
) -> Option<KarlinParams> {
    query_ungapped_karlin_blastna(
        &encode_iupac_to_blastna(&query.to_ascii_uppercase()),
        reward,
        penalty,
    )
}

/// `query_ungapped_karlin` of a context already in BLASTNA (the search holds its contexts
/// only in that encoding).
pub(crate) fn query_ungapped_karlin_blastna(
    blastna: &[u8],
    reward: i32,
    penalty: i32,
) -> Option<KarlinParams> {
    const BLASTNA_N: usize = 14;
    const BLASTNA_GAP: usize = 15;
    let mut counts = [0f64; 16];
    for &code in blastna {
        counts[(code & 0x0f) as usize] += 1.0;
    }
    counts[BLASTNA_N] = 0.0;
    counts[BLASTNA_GAP] = 0.0;
    let sum: f64 = counts.iter().sum();
    if sum == 0.0 {
        return None;
    }
    let matrix = build_blastna_matrix(reward, penalty);
    // At most 16 residues x 4 bases: a fixed array in place of a vector grown for each query.
    let mut probabilities = [(0i32, 0f64); 64];
    let mut num_probabilities = 0usize;
    for (residue, &count) in counts.iter().enumerate() {
        if count == 0.0 {
            continue;
        }
        for base in 0..4 {
            probabilities[num_probabilities] = (matrix[residue * 16 + base], count / sum * 0.25);
            num_probabilities += 1;
        }
    }
    let sfp = score_freq_profile_from_probabilities(
        penalty.min(reward),
        reward.max(penalty),
        &probabilities[..num_probabilities],
    );
    compute_karlin_params_ungapped(&sfp).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(seq: &[u8]) -> FastaRecord {
        FastaRecord::new("Query_1", b"q", seq)
    }

    /// The letters of a query searched whole (no `-query_loc`).
    fn whole(letters: &[u8]) -> QueryRecord<'_> {
        QueryRecord { letters, offset: 0 }
    }

    fn hsp(
        frame: i32,
        q: (usize, usize),
        s: (usize, usize),
        ops: Option<Vec<GapEditOp>>,
    ) -> BlastnHsp {
        BlastnHsp {
            identity: 0.0,
            length: 0,
            mismatch: 0,
            gapopen: 0,
            q_start: 0,
            q_end: 0,
            s_start: 0,
            s_end: 0,
            e_value: 0.0,
            bit_score: 0.0,
            num_ident: 0,
            query_frame: frame,
            // The context's length, as the engine sets it: the tests' queries have 8 letters.
            query_length: 8,
            q_idx: 0,
            s_idx: 0,
            raw_score: 0,
            internal_q_offset_0: q.0 as i32,
            internal_q_end_0: q.1 as i32,
            internal_s_offset_0: s.0,
            internal_s_end_0: s.1,
            internal_query_context_offset: 0,
            gap_info: ops,
            num_positives: 0,
        }
    }

    // seqinfosrc_seqdb.cpp:121-133 (quoted at `shown_subject_masks`): a stored mask is shown
    // when it holds a letter and meets the closed range `[offset, end]` of an HSP (the end one
    // past the HSP); overlapping and nested masks (hand-made masks of the M4a-2 edge databases)
    // are merged for the display. A subject set's merged lowercase is kept as it is.
    #[test]
    fn database_masks_meeting_an_hsp_are_shown() {
        let masks = [
            MaskedInterval::new(600, 700),
            MaskedInterval::new(100, 200),
            MaskedInterval::new(150, 250),
            MaskedInterval::new(250, 260),
            MaskedInterval::new(300, 300),
            MaskedInterval::new(900, 950),
        ];
        // The mask 600-700 holds 600..=699: it meets an HSP from 699 on, not one from 700 on.
        let hsps = [
            hsp(1, (0, 10), (199, 250), None),
            hsp(1, (0, 10), (700, 760), None),
        ];
        assert!(shown_subject_masks(&masks[..1], &hsps).is_empty());
        let hsps = [
            hsp(1, (0, 10), (199, 250), None),
            hsp(1, (0, 10), (699, 701), None),
        ];
        // 100-200, 150-250 and 250-260 meet the first HSP, 600-700 the second; the empty 300-300
        // and 900-950 none.
        assert_eq!(
            shown_subject_masks(&masks, &hsps),
            [MaskedInterval::new(100, 260), MaskedInterval::new(600, 700)]
        );
        // Merged lowercase runs are kept as they are.
        let runs = merged(&masks);
        assert_eq!(
            shown_subject_masks(&runs, &hsps),
            [MaskedInterval::new(100, 260), MaskedInterval::new(600, 700)]
        );
        let hsps = [hsp(1, (0, 10), (940, 960), None)];
        assert_eq!(
            shown_subject_masks(&runs, &hsps),
            [MaskedInterval::new(900, 950)]
        );
    }

    #[test]
    fn plus_strand_rows_follow_the_edit_script() {
        let query = record(b"ACGTTACG");
        let subject = record(b"ACGTACCG");
        let ops = vec![GapEditOp::Sub(4), GapEditOp::Ins(1), GapEditOp::Sub(3)];
        let (q, s) = displayed_rows(
            &hsp(1, (0, 8), (0, 7), Some(ops)),
            whole(query.seq()),
            subject.seq().len(),
            |from, to| subject.seq()[from..to].to_ascii_uppercase(),
        )
        .unwrap();
        assert_eq!(
            (q.as_slice(), s.as_slice()),
            (&b"ACGTTACG"[..], &b"ACGT-ACC"[..])
        );
    }

    // A minus-strand HSP is shown with the query on its plus strand: both rows of the
    // context alignment are reverse-complemented, gaps included.
    #[test]
    fn minus_strand_rows_are_reverse_complemented() {
        let query = record(b"AACCGGTT");
        let subject = record(b"CCGGT");
        // Context (reverse complement of the query) AACCGGTT[1..6] = ACCGG against CCGG-.
        let ops = vec![GapEditOp::Ins(1), GapEditOp::Sub(4)];
        let (q, s) = displayed_rows(
            &hsp(-1, (1, 6), (0, 4), Some(ops)),
            whole(query.seq()),
            subject.seq().len(),
            |from, to| subject.seq()[from..to].to_ascii_uppercase(),
        )
        .unwrap();
        assert_eq!((q.as_slice(), s.as_slice()), (&b"CCGGT"[..], &b"CCGG-"[..]));
    }

    // remaining-gaps X1b-p: the query letters of an HSP that NCBI keeps outside its query context
    // with its Seq-align's query start at or past 0 (`context_letters`, alnvec.cpp:279-300): a
    // column whose position on the query's plus strand lies past the query's end shows the gap
    // character, the others their letters. Plus strand: offsets past the end; minus strand:
    // offsets below 0 (positions past the end of the plus strand); a range that crosses the end
    // shows the letters before it. A start below 0 is refused (it never reaches the rows).
    #[test]
    fn query_letters_past_the_query_end_show_the_gap_character() {
        let query = record(b"AACCGGTT");
        let subject = record(b"ACGTACGTAC");
        let letters = |from: usize, to: usize| subject.seq()[from..to].to_ascii_uppercase();
        let outside = |frame: i32, offset: i32, end: i32, s: (usize, usize), ops| {
            let mut hsp = hsp(frame, (0, 0), s, Some(ops));
            (hsp.internal_q_offset_0, hsp.internal_q_end_0) = (offset, end);
            hsp
        };
        let gapped = || vec![GapEditOp::Sub(2), GapEditOp::Del(1), GapEditOp::Sub(1)];
        for (hsp, start, rows) in [
            (outside(1, 10, 13, (0, 4), gapped()), 10, ("----", "ACGT")),
            (outside(-1, -5, -2, (0, 4), gapped()), 10, ("----", "ACGT")),
            (
                outside(1, 6, 10, (0, 4), vec![GapEditOp::Sub(4)]),
                6,
                ("TT--", "ACGT"),
            ),
            // Context "--AAC" (record positions 9, 8, 7, 6, 5), shown reverse-complemented.
            (
                outside(-1, -2, 3, (0, 5), vec![GapEditOp::Sub(5)]),
                5,
                ("GTT--", "TACGT"),
            ),
        ] {
            assert_eq!(hsp.seqalign_query_start(), start);
            let (q, s) =
                displayed_rows(&hsp, whole(query.seq()), subject.seq().len(), letters).unwrap();
            assert_eq!(
                (q.as_slice(), s.as_slice()),
                (rows.0.as_bytes(), rows.1.as_bytes())
            );
        }
        let left = outside(1, -4, -1, (0, 3), vec![GapEditOp::Sub(3)]);
        assert_eq!(left.seqalign_query_start(), -4);
        assert!(displayed_rows(&left, whole(query.seq()), subject.seq().len(), letters).is_err());
        let right_minus = outside(-1, 9, 12, (0, 3), vec![GapEditOp::Sub(3)]);
        assert_eq!(right_minus.seqalign_query_start(), -4);
        assert!(displayed_rows(
            &right_minus,
            whole(query.seq()),
            subject.seq().len(),
            letters
        )
        .is_err());
    }

    // remaining-gaps X1b-p review 1 (M1): with `-query_loc` NCBI moves the query row of the
    // Seq-align by the interval's start (`RemapToQueryLoc`) and reads the letters of the query's
    // Bioseq, the whole record (alnvec.cpp:108-124): an HSP outside its context shows the record's
    // letters outside the interval, and the gap character only past the record's end. The record
    // GGAACCGGTTAC is searched from position 2 (the interval AACCGGTT, the context of 8 letters):
    // past the interval's end on the plus strand and left of the minus context (offsets below 0),
    // then a start below 0 that the interval's start moves to the record's first letters (an
    // ungapped search's Dense-diag, `BlastnHsp::record_query_start`), on both strands. A start
    // that stays below 0 in the record never reaches the rows.
    #[test]
    fn query_letters_outside_the_interval_are_the_records() {
        let query = record(b"GGAACCGGTTAC");
        let in_record = QueryRecord {
            letters: query.seq(),
            offset: 2,
        };
        let subject = record(b"ACGTACGTAC");
        let letters = |from: usize, to: usize| subject.seq()[from..to].to_ascii_uppercase();
        let outside = |frame: i32, offset: i32, end: i32| {
            let n = (end - offset) as usize;
            let mut hsp = hsp(frame, (0, 0), (0, n), Some(vec![GapEditOp::Sub(n as u32)]));
            (hsp.internal_q_offset_0, hsp.internal_q_end_0) = (offset, end);
            hsp
        };
        for (hsp, start, row) in [
            (outside(1, 6, 10), 8, "TTAC"),
            (outside(1, 8, 11), 10, "AC-"),
            (outside(1, 10, 14), 12, "----"),
            // Context offsets -2..3: record positions 11 down to 7.
            (outside(-1, -2, 3), 7, "GTTAC"),
            (outside(1, -2, 1), 0, "GGA"),
            // Context offsets 7..10: record positions 2 down to 0.
            (outside(-1, 7, 10), 0, "GGA"),
        ] {
            assert_eq!(hsp.record_query_start(in_record.offset, true), start);
            let (q, _) = displayed_rows(&hsp, in_record, subject.seq().len(), letters).unwrap();
            assert_eq!(q.as_slice(), row.as_bytes());
        }
        let left = outside(1, -4, -1);
        assert_eq!(left.record_query_start(2, true), -2);
        assert!(displayed_rows(&left, in_record, subject.seq().len(), letters).is_err());
        // A gapped search's Dense-seg keeps a start below 0 (`CDense_seg::OffsetRow`).
        assert_eq!(outside(1, -2, 1).record_query_start(2, false), -2);
    }

    // The records the rows read (`QueryRecords::record`): the searched record itself without
    // `-query_loc`; with it the input record of the query's input position (records that the
    // interval leaves out counted), from the call's first input record, and the interval's start.
    #[test]
    fn query_records_map_searched_queries_to_their_input_records() {
        use crate::blastinput::seq_range::{Placements, QueryInput, RecordPlacement};
        let input = [record(b"AAAA"), record(b"CC"), record(b"GGGGGG")];
        let searched = [record(b"AA"), record(b"GGG")];
        let ranged = QueryInput {
            placements: Placements(vec![
                RecordPlacement {
                    offset: 1,
                    length: 4,
                },
                RecordPlacement {
                    offset: 1,
                    length: 6,
                },
            ]),
            ordinals: vec![0, 2],
            input_lengths: vec![4, 2, 6],
            skipped: vec![false, true, false],
            range_length: Some(2),
        };
        let records = QueryRecords::new(&searched, &ranged, &input, 0);
        let second = records.record(1).unwrap();
        assert_eq!((second.letters, second.offset), (&b"GGGGGG"[..], 1));
        assert!(records.record(2).is_none());
        // A batch's records from input record 2 on.
        let batch = QueryRecords::new(&searched, &ranged, &input[2..], 2);
        assert_eq!(batch.record(1).unwrap().letters, b"GGGGGG");
        assert!(batch.record(0).is_none());
        let plain = QueryInput::whole(&searched);
        let unranged = QueryRecords::new(&searched, &plain, &[], 0)
            .record(1)
            .unwrap();
        assert_eq!((unranged.letters, unranged.offset), (&b"GGG"[..], 0));
    }

    // The XML rows are written into buffers that the caller reuses: stale contents are cleared
    // and the minus-strand reversal in place gives the rows of `displayed_rows`.
    #[test]
    fn rows_into_reused_buffers_equal_displayed_rows() {
        let query = record(b"AACCGGTT");
        let subject = record(b"CCGGT");
        let (mut q, mut s) = (b"stale row".to_vec(), b"x".to_vec());
        for (frame, window, ops) in [
            (-1, (1, 6), vec![GapEditOp::Ins(1), GapEditOp::Sub(4)]),
            (1, (0, 5), vec![GapEditOp::Sub(5)]),
        ] {
            let hsp = hsp(frame, window, (0, 4 + usize::from(frame > 0)), Some(ops));
            let letters = |from: usize, to: usize| subject.seq()[from..to].to_ascii_uppercase();
            let once =
                displayed_rows(&hsp, whole(query.seq()), subject.seq().len(), letters).unwrap();
            displayed_rows_into(
                &hsp,
                whole(query.seq()),
                subject.seq().len(),
                letters,
                &mut q,
                &mut s,
            )
            .unwrap();
            assert_eq!((q.clone(), s.clone()), once);
        }
    }

    #[test]
    fn masked_positions_are_lowercase_in_the_row_direction() {
        let mut row = b"AC-GT".to_vec();
        lowercase_masked(&mut row, 10, -1, &[MaskedInterval { start: 7, end: 9 }]);
        assert_eq!(row, b"Ac-gT");
    }

    // A context in BLASTNA (either letter case, no uppercase copy) gives the block of the
    // letters it was encoded from.
    #[test]
    fn karlin_of_blastna_equals_karlin_of_letters() {
        for seq in [
            &b"ACGTACGTTTGA"[..],
            b"acgtnACGTRYKMSWBDHVnn-",
            b"NNNN",
            b"",
            b"aaaaCCCCggggTTTTUu",
        ] {
            let blastna = encode_iupac_to_blastna(seq);
            let from_letters = query_ungapped_karlin(seq, 2, -3);
            let from_blastna = query_ungapped_karlin_blastna(&blastna, 2, -3);
            match (from_letters, from_blastna) {
                (None, None) => {}
                (Some(a), Some(b)) => {
                    assert_eq!(a.lambda.to_bits(), b.lambda.to_bits());
                    assert_eq!(a.k.to_bits(), b.k.to_bits());
                    assert_eq!(a.h.to_bits(), b.h.to_bits());
                }
                other => panic!("{seq:?}: {other:?}"),
            }
        }
    }

    #[test]
    fn an_all_n_query_has_no_ungapped_block() {
        assert!(query_ungapped_karlin(b"NNNN", 2, -3).is_none());
        let acgt = query_ungapped_karlin(b"ACGTACGTTTGA", 2, -3).unwrap();
        assert!((acgt.lambda - 0.634).abs() < 5e-4);
    }
}
