//! BLASTP's ungapped HSPs of length 0 (`-ungapped` with `-window_size 0` and a low `-threshold`,
//! which `Blast_HSPTest` keeps): what the reports do with one that the final hit lists show.
//!
//! The rows, blocks and documents of such an HSP are the ordinary writers' (their lengths and
//! coordinates already allow 0, `alignment::build_alignment_view_with_matrix` gives it empty
//! rows). Two cells differ:
//! - outfmt 1 to 4: the merge of the query-anchored alignments frees a segment that it still reads
//!   (NCBI crashes), which LOSAT refuses by name where `query_anchored::section_frees` finds it
//!   (`refuse_merge_frees`);
//! - `-subject`: an HSP at the first letter of the searched part of a subject makes the Seq-loc
//!   mapper throw an exception before any query's report (`subject_start_exception`), with the
//!   formats' prologs before it (`write_exception_prolog`).

use std::collections::BTreeSet;
use std::io::Write;
use std::ops::Range;
use std::sync::Arc;

use anyhow::Result;

use super::report_write::{write_pairwise, ReportView};
use super::BlastpReport;
use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::api::local_blast::{FormatOutput, FormatProbe};
use crate::blastinput::seq_range::Placements;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::pairwise::{PairwiseHit, PairwiseView};
use crate::report::query_warnings::QueryWarnings;

/// The refusal of outfmt 1 to 4 with an ungapped HSP of length 0 among the alignments shown where
/// NCBI BLAST+ crashes (`refuse_merge_frees`).
pub(super) const UNPORTED_ANCHORED_FREE: &str = "BLASTP -ungapped -window_size 0 with -outfmt 1 to 4: a one-hit ungapped alignment of length 0 among the alignments shown where NCBI BLAST+'s merge of the query-anchored alignments frees memory that it still reads (it crashed in our runs), which LOSAT does not support";

/// The HSPs that a final hit list shows (its first `hsplist_count` lists).
fn shown_hsps(
    hit_list: &BlastpHitList,
) -> impl Iterator<Item = &crate::algorithm::blastp::hsp::BlastpHsp> {
    hit_list
        .hsplist_array
        .iter()
        .take(hit_list.hsplist_count)
        .flat_map(|list| &list.hsps)
}

/// Whether a final hit list shows an HSP of length 0 (query `q_idx`'s list `hit_list`).
fn shows_empty(hit_list: &Option<BlastpHitList>) -> bool {
    hit_list
        .as_ref()
        .is_some_and(|list| shown_hsps(list).any(|hsp| hsp.length == 0))
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1285-1294
/// `CAnnotMapperException::eBadAlignment`, thrown by the Seq-loc mapper's Dense-diag of a mapped
/// Seq-align (`x_GetDstDendiag`) whose row starts at `kInvalidSeqPos`.
fn bad_alignment() -> NcbiException {
    NcbiException::new(
        ThrowSite {
            file: "src/objects/seq/seq_align_mapper_base.cpp",
            line: 1292,
            kind: "CAnnotMapperException::eBadAlignment",
            function: None,
        },
        "Mapped alignment contains gaps and can not be converted to dense-diag.",
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1541-1547
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:1403-1437
/// The exception that ends a `-subject` search (also with `-subject_loc`) whose final hit lists
/// show an HSP of length 0 at the first letter of the searched part of its subject: NCBI's
/// `s_RemapToSubjectLoc` maps each subject's Dense-diags into the subject's record, and the row
/// of a segment of length 0 at start 0 wraps. It is thrown for the first such HSP while the
/// results are made, before any query's report. `hit_lists` are in record coordinates
/// (`shift_hsp_to_records`), `placements` the subjects' searched parts; a database search maps
/// nothing.
pub(super) fn subject_start_exception(
    hit_lists: &[Option<BlastpHitList>],
    placements: &Placements,
    subject_set: bool,
) -> Option<NcbiException> {
    if !subject_set {
        return None;
    }
    hit_lists
        .iter()
        .flatten()
        .flat_map(shown_hsps)
        .any(|hsp| hsp.length == 0 && hsp.s_start - 1 == placements.offset(hsp.s_idx as usize))
        .then(bad_alignment)
}

/// Writes what NCBI has written when `subject_start_exception` ends the run: the prolog of the
/// format (`PrintProlog`: outfmt 0 to 4's banner and database line, outfmt 20's field names, the
/// opening of the XML2 and JSON documents) and nothing more. The BLAST archive's error document
/// is the run's own (`write_run_error_archive`). The pairwise report is written for no query
/// (`ReportView::prolog_only`).
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:244-248
#[allow(clippy::too_many_arguments)]
pub(super) fn write_exception_prolog(
    view: &ReportView<'_>,
    outfmt: BlastpReport,
    format: &mut FormatOutput<'_>,
    warnings: &mut Option<QueryWarnings<'_>>,
    probe: &mut Option<FormatProbe<'_>>,
    pairwise_shown: &Option<BTreeSet<usize>>,
    rows_of: &impl Fn(Range<usize>) -> Result<Vec<PairwiseHit>>,
) -> Result<()> {
    match outfmt {
        BlastpReport::Pairwise(pairwise) => write_pairwise(
            view,
            pairwise,
            format,
            warnings,
            probe,
            pairwise_shown,
            rows_of,
            None,
        ),
        BlastpReport::CustomTabular(spec) => {
            let mut writer = crate::cli::NcbiOstream::sink(&mut format.sink)?;
            crate::algorithm::blastp::tabular::write_tabular_prolog(spec, &mut writer)?;
            writer.flush()?;
            Ok(())
        }
        BlastpReport::Xml2(kind) => {
            let mut writer = format.sink.open()?;
            if view.reports.x_prolog {
                writer.write_all(&crate::report::blast_xml2::framing::prolog_bytes(kind))?;
            }
            writer.flush()?;
            Ok(())
        }
        BlastpReport::Xml
        | BlastpReport::Text(_)
        | BlastpReport::Seqalign(_)
        | BlastpReport::Archive => Ok(()),
    }
}

/// Refuses, before any byte is written, a query-anchored report (outfmt 1 to 4, `view`) of an
/// `-ungapped` search whose final hit lists show an HSP of length 0 where NCBI's merge of the
/// section frees a segment that it still reads (`pairwise::blastp_section_frees`), unless the
/// search ends with `subject_start_exception` first. `rows_of` makes the rows of a window of
/// queries (only those of a query with such an HSP are made); `joins` tells the copies of one
/// database record shown in one hit.
#[allow(clippy::too_many_arguments)]
pub(super) fn refuse_merge_frees(
    view: PairwiseView,
    sort: crate::report::hit_sort::ReportSort,
    num_alignments: usize,
    hit_lists: &[Option<BlastpHitList>],
    pairwise_hits: Option<&[PairwiseHit]>,
    rows_of: &impl Fn(Range<usize>) -> Result<Vec<PairwiseHit>>,
    subject_titles: &[Arc<[u8]>],
    joins: &impl Fn(u32, u32) -> bool,
    throws_first: bool,
) -> Result<()> {
    if throws_first || view == PairwiseView::Pairwise {
        return Ok(());
    }
    for (q_idx, _) in hit_lists
        .iter()
        .enumerate()
        .filter(|(_, hit_list)| shows_empty(hit_list))
    {
        let rows: std::borrow::Cow<'_, [PairwiseHit]> = match pairwise_hits {
            Some(hits) => std::borrow::Cow::Borrowed(hits),
            None => std::borrow::Cow::Owned(rows_of(q_idx..q_idx + 1)?),
        };
        if crate::report::pairwise::blastp_section_frees(
            view,
            sort,
            num_alignments,
            q_idx,
            &rows,
            subject_titles,
            joins,
        )? {
            anyhow::bail!(UNPORTED_ANCHORED_FREE);
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "empty_alignments_tests.rs"]
mod tests;
