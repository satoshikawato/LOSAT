//! psiblast's outfmt 0 report of PHI-BLAST's round (`-phi_pattern`, series PS, item PS8a):
//! `CBlastFormat::PrintPhiResult` of one query. The head is a round's (`Results from round 1`, the
//! query's header), then the pattern's line block (`PrintPhiInfo`), a description table for each
//! occurrence of the pattern in the query, and for each occurrence its alignments under
//! `Significant alignments for pattern occurrence i at position p`, whose rows show the
//! occurrence's query columns in a `pattern` row (`PatternRow`, `pairwise::RowDecoration`); then
//! the query's footer. When the first occurrence has no alignment the report says `No hits found`,
//! whatever the other occurrences found.
//!
//! The writer is a child of `psiblast_round` and writes through its writer type and its sections
//! (`RoundSections`), so that no default report's instance is shared.

use super::super::{
    cpp_default_double, write_description_table, write_no_hits_found, write_query_preamble,
    write_spaces, write_tblastn_unsearched_query_footer_spacing, BlastpPairwiseQuery,
    BlastpPairwiseReport, PairwiseConfig, PairwiseDatabase, PairwiseHit, RowDecoration,
};
use super::{ReportRound, RoundSections, RoundWriter};
use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::report::ostream::{DisplaySeqalign, ReportWrite};
use crate::report::query_warnings::ReportPlace;
use crate::report::x_pairwise_par::{QueryGroup, QuerySections};
use std::io::{self, Write};
use std::sync::Arc;

/// The pattern of a query's report (`SPHIQueryInfo`): the pattern's text as NCBI keeps it (upper
/// case), its probability and the query's occurrences (0-based offset and length, in the searched
/// range of the query).
pub(crate) struct PhiQueryInfo<'a> {
    pub(crate) text: &'a [u8],
    pub(crate) probability: f64,
    pub(crate) occurrences: Vec<(u32, u32)>,
}

/// The footer of a PHI-BLAST report (`x_PrintOneQueryFooter` of its query's context).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:445-478
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PhiFooter {
    /// The query is searched: the pattern's Karlin blocks and search space, which the caller has
    /// put in the query's footer values (`phi_hooks::phi_footers`).
    Searched,
    /// The report's query has no letters (its batch's query with letters is searched): no Karlin
    /// block, a search space of 0.
    NoLetters,
    /// NCBI's table has no Karlin block for the matrix and gap costs: -1 in both blocks, no Gumbel
    /// columns and a search space of 0, also with `-searchsp`, as for a query that is not
    /// searched.
    NoStatistics,
}

/// The `pattern` row of an occurrence's alignments: `*` over the alignment columns of the query's
/// positions `from..=to` (the occurrence's offsets, which NCBI's feature takes as positions of the
/// query's record also under `-query_loc`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2142-2150
pub(crate) struct PatternRow {
    pub(crate) from: u32,
    pub(crate) to: u32,
}

impl RowDecoration for PatternRow {
    const FEATURE_ID: &'static str = "pattern";
    const FEATURE_CHAR: u8 = b'*';

    fn feature(&self, rows: &crate::report::query_anchored::HspAlnVec) -> Option<(usize, usize)> {
        rows.query_columns(self.from, self.to)
    }
}

/// The outfmt 0 report of query `q_idx` of PHI-BLAST's round: `occurrence_hits[i]` are the rows of
/// the HSPs of occurrence `i` (its result set), in its result's order, after the warnings of the
/// queries `batch` (the batch's reading, from `q_idx` on) and the messages `result_warning` of the
/// results (`query_warnings::phi_result_warning`), and the footer `footer`. No prolog, no epilog.
/// `info` is `None` where NCBI's table had no Karlin block and nothing was searched: no pattern
/// block, no hits.
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2039-2047
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2068-2199
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_phi_round_report(
    occurrence_hits: &[Vec<PairwiseHit>],
    batch: std::ops::Range<usize>,
    result_warning: &[u8],
    info: Option<&PhiQueryInfo<'_>>,
    footer: PhiFooter,
    writer: &mut dyn ReportWrite,
    config: &PairwiseConfig,
    queries: &[BlastpPairwiseQuery],
    query_titles: &[Arc<[u8]>],
    subject_titles: &[Arc<[u8]>],
    parsed_subjects: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    subject_sources: &[crate::report::defline::SourceTaxname],
    loads: &crate::report::db_subject::HitLoads,
    report: &BlastpPairwiseReport,
    database: Option<&PairwiseDatabase>,
    mut probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
) -> io::Result<()> {
    let mut buffered = io::BufWriter::new(RoundWriter(writer));
    let writer = &mut buffered;
    let q_idx = batch.start;
    if let Some(warnings) = warnings.as_deref_mut() {
        for query in batch {
            warnings.before_query(query, writer)?;
        }
    }
    writer.write_all(result_warning)?;
    let sections = RoundSections {
        config,
        queries,
        query_titles,
        report,
        database,
        parsed_subjects,
        subject_sources,
        loads,
        round: ReportRound {
            number: 1,
            found_before: None,
        },
        anchored: None,
    };
    // Each occurrence's group, its HSPs numbered over the query's report.
    let mut first_index: HspIndex = 0;
    let indexed: Vec<Vec<(HspIndex, &PairwiseHit)>> = occurrence_hits
        .iter()
        .map(|hits| {
            let start = first_index;
            first_index += hits.len();
            hits.iter()
                .enumerate()
                .map(|(offset, hit)| (start + offset, hit))
                .collect()
        })
        .collect();
    let groups: Vec<Option<QueryGroup<'_>>> = indexed
        .iter()
        .map(|hits| {
            QueryGroup::new(hits, subject_titles, |previous, s_idx| {
                sections.joins(previous, s_idx)
            })
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2068-2080
    let query_title = query_titles
        .get(q_idx)
        .map_or(&b""[..], |title| title.as_ref());
    write_query_preamble(
        writer,
        Some(sections.round.number),
        query_title,
        queries[q_idx].query_length,
        false,
    )?;
    if let Some(info) = info {
        write_phi_info(writer, info)?;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2108-2115
    if groups.first().is_none_or(Option::is_none) {
        write_no_hits_found(writer)?;
        return write_phi_footer(writer, &sections, q_idx, footer).and_then(|()| writer.flush());
    }
    // NCBI loads each subject that the query's batch did not prefetch at its description line,
    // or one the tables do not show at its alignment's heading (`QueryWarnings::subject_read`).
    let mut place = |place: ReportPlace, writer: &mut io::BufWriter<RoundWriter<'_>>| {
        crate::report::x_pairwise_par::at_place(warnings.as_deref_mut(), None, q_idx, place, writer)
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2118-2128
    for group in &groups {
        match group {
            Some(group) => write_description_table(
                writer,
                group.described(report.num_descriptions),
                group.subject_order.len() <= report.num_descriptions,
                &group.subject_hits,
                false,
                // psiblast does not hand `-sorthits` to the formatter (psiblast_app.cpp:502).
                None,
                |s_idx| sections.label(group, s_idx),
                &mut place,
                // PHI-BLAST refuses `-html` (`phi_hooks::phi_limits`): the plain table.
                None,
                None,
            )?,
            None => write_empty_table(writer, report.num_descriptions)?,
        }
    }
    writeln!(writer)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2136-2196
    let occurrences = info.map_or(&[][..], |info| &info.occurrences[..]);
    for (index, (group, &(offset, length))) in groups.iter().zip(occurrences).enumerate() {
        write!(
            writer,
            "\nSignificant alignments for pattern occurrence {} at position {}\n\n",
            index + 1,
            1 + i64::from(offset)
        )?;
        if let Some(group) = group {
            let pattern = PatternRow {
                from: offset,
                to: (offset + length).saturating_sub(1),
            };
            let aligned = group.aligned(report.num_alignments);
            let end = crate::report::x_pairwise_par::aligned_end(
                config.after_throw(),
                aligned,
                |s_idx| sections.thrown(s_idx),
            );
            // Each occurrence's alignments are one `CDisplaySeqalign` (its region and its flush).
            let mut display = DisplaySeqalign::default();
            for &s_idx in &aligned[..end.map_or(aligned.len(), |(at, _)| at)] {
                place(ReportPlace::Load(s_idx), &mut *writer)?;
                display.enter(writer)?;
                sections.phi_subject(writer, probe.as_deref_mut(), group, s_idx, &pattern)?;
            }
            if let Some((_, exception)) = end {
                if let Some(warnings) = warnings.as_deref_mut() {
                    warnings.end_run();
                }
                writer.flush()?;
                return Err(exception.clone().into_io_error());
            }
            display.leave(writer)?;
        }
        writeln!(writer)?;
    }
    write_phi_footer(writer, &sections, q_idx, footer)?;
    writer.flush()
}

/// The footer `footer` of query `q_idx`'s report (`PhiFooter`).
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:445-478
fn write_phi_footer<W: ReportWrite>(
    writer: &mut W,
    sections: &RoundSections<'_>,
    q_idx: usize,
    footer: PhiFooter,
) -> io::Result<()> {
    match footer {
        PhiFooter::Searched => sections.foot(writer, q_idx),
        PhiFooter::NoLetters => RoundSections::invalid_foot(writer),
        PhiFooter::NoStatistics => write_tblastn_unsearched_query_footer_spacing(writer, false),
    }
}

/// `PrintPhiInfo`: the number of occurrences, the pattern, their 1-based positions and the
/// pattern's probability (a `double` on a C++ stream, `cpp_default_double`). The pattern is read
/// as a C string (the info's `pattern` member made a `std::string`): it ends at a NUL byte.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:757-779
fn write_phi_info<W: Write>(writer: &mut W, info: &PhiQueryInfo<'_>) -> io::Result<()> {
    writeln!(
        writer,
        "{} occurrence(s) of pattern: ",
        info.occurrences.len()
    )?;
    let text = info
        .text
        .split(|&byte| byte == 0)
        .next()
        .unwrap_or_default();
    writer.write_all(text)?;
    writer.write_all(b" at position(s) ")?;
    for (index, &(offset, _)) in info.occurrences.iter().enumerate() {
        if index > 0 {
            writer.write_all(b", ")?;
        }
        write!(writer, "{}", 1 + i64::from(offset))?;
    }
    writeln!(writer, " of query sequence")?;
    writeln!(
        writer,
        "pattern probability={}",
        cpp_default_double(info.probability)
    )
}

/// The description table of an occurrence without alignments (`CShowBlastDefline` of an empty
/// set): the header's two lines of the widths of `(Bits)` and `Value` and the blank line after
/// them (none with `-num_descriptions 0`, `eNoShowHeader`), then `x_DisplayDeflines`'s line feed.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:756-848
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:605-613
fn write_empty_table<W: Write>(writer: &mut W, num_descriptions: usize) -> io::Result<()> {
    if num_descriptions > 0 {
        let bits = super::super::ncbi_k_bits();
        let width = bits.len();
        write_spaces(writer, 70)?;
        writeln!(writer, "{:<width$}    E", "Score")?;
        write!(
            writer,
            "{:<69}",
            "Sequences producing significant alignments:"
        )?;
        writeln!(writer, "{bits:<width$}  Value")?;
        writeln!(writer)?;
    }
    writeln!(writer)
}

#[cfg(test)]
#[path = "phi_round_tests.rs"]
mod tests;
