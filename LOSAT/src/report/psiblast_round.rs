//! psiblast's outfmt 0 report of one round of a FASTA query (series PS, item PS5a): the report of
//! `write_blastp_pairwise_report` with the round's number in its head and, after the first round,
//! the description table in two parts, the subjects used in the PSSM's model and found again and
//! those found for the first time (`CBlastFormat::x_DisplayDeflines` with the previous round's
//! inclusion set). The query-anchored reports of outfmt 1 to 4 are the same report with the one
//! section of their view in place of each subject's alignments (item PS15b-2). The round driver
//! (`blastp::blast_engine::psi_round`) writes the epilog after the last round; this writer never
//! does.
//!
//! The writer is a child of `pairwise` and uses its parts. It writes through a type of its own
//! (`RoundWriter`), so that every generic part it calls is an instance of its own and the default
//! report's instances are not shared (the copy rule of DESIGN_ps4 §4.3), and formats the one query
//! of a round serially (no `LOSAT_X_PAIRWISEPAR` pieces).

use super::{
    alignment_id_local, blastp_family_program, database_description_label, shown_hsp, table_scores,
    write_blast_database_header, write_blastp_database_header_spacing,
    write_blastp_query_footer_spacing, write_blastp_query_header, write_blastp_subject_part,
    write_description_row_start, write_description_table, write_no_hits_found, write_prolog,
    write_spaces, write_tblastn_unsearched_query_footer_spacing, BlastpPairwiseQuery,
    BlastpPairwiseReport, HtmlRowMasks, HtmlRows, PairwiseConfig, PairwiseDatabase, PairwiseHit,
    PairwiseView,
};
use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::report::ostream::ReportWrite;
use crate::report::outfmt6::{format_bitscore_ncbi, format_evalue_ncbi};
use crate::report::query_warnings::{ReportHook, ReportPlace};
use crate::report::x_pairwise_par::{QueryGroup, QuerySections};
use std::io::{self, Write};
use std::sync::Arc;

/// The report's writer: the caller's stream behind a type that no default report writes through.
struct RoundWriter<'w>(&'w mut dyn ReportWrite);

impl ReportWrite for RoundWriter<'_> {
    fn site(&mut self, site: crate::report::ostream::Site) -> io::Result<()> {
        self.0.site(site)
    }

    fn ncbi_flush(&mut self) -> io::Result<()> {
        self.0.ncbi_flush()
    }
}

impl Write for RoundWriter<'_> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0.write(buf)
    }

    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        self.0.write_all(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.0.flush()
    }
}

/// The round of a report: its number (`itr.GetIterationNumber()`) and, when the previous round's
/// inclusion set (`itr.GetPreviouslyFoundSeqIds()`) is not empty, whether each subject (by subject
/// index, a database's by rank) is one of it.
#[derive(Clone, Copy)]
pub(crate) struct ReportRound<'a> {
    pub(crate) number: u32,
    pub(crate) found_before: Option<&'a [bool]>,
}

/// The outfmt 0 report of the queries `query_range` of one round of `queries` (the one query of a
/// round), whose `hits` are those queries' (the rows of a window of queries,
/// `pairwise::write_blastp_pairwise_window`): the prolog when `prolog` is set (the run's first
/// round) and the window starts the report, each query's warnings, then its section. No epilog.
/// `anchored`: the sequences of the query-anchored section of outfmt 1 to 4 (`config.view`;
/// `None` for outfmt 0).
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:278-286
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_psiblast_round_report(
    hits: &[PairwiseHit],
    query_range: std::ops::Range<usize>,
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
    probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut super::super::query_warnings::QueryWarnings<'_>>,
    prolog: bool,
    round: ReportRound<'_>,
    anchored: Option<&dyn crate::report::query_anchored::ReportSequences>,
) -> io::Result<()> {
    let mut buffered = io::BufWriter::new(RoundWriter(writer));
    let writer = &mut buffered;
    if prolog && query_range.start == 0 {
        let prolog_result = write_prolog(
            writer,
            blastp_family_program(&config.program),
            &report.version,
            false,
            false,
            false,
            super::html_with_diagnostics(config.html_links(), warnings.as_deref_mut())?,
        );
        if let Err(error) = prolog_result {
            let _ = writer.flush();
            return Err(error);
        }
        match database {
            Some(database) => write_blast_database_header(writer, database)?,
            None => write_blastp_database_header_spacing(
                writer,
                &report.database_name,
                report.database_num_sequences,
                report.database_total_letters,
                1,
            )?,
        }
    }
    let mut hits_by_query: Vec<Vec<(HspIndex, &PairwiseHit)>> = vec![Vec::new(); query_range.len()];
    for (hsp_index, hit) in hits.iter().enumerate() {
        if let Some(bucket) = (hit.hit.q_idx as usize)
            .checked_sub(query_range.start)
            .and_then(|offset| hits_by_query.get_mut(offset))
        {
            bucket.push((hsp_index, hit));
        }
    }
    let sections = RoundSections {
        config,
        queries,
        query_titles,
        report,
        database,
        parsed_subjects,
        subject_sources,
        loads,
        round,
        anchored,
    };
    let written = crate::report::x_pairwise_par::write(
        writer,
        &sections,
        query_range.start,
        &hits_by_query,
        subject_titles,
        probe,
        warnings,
        None,
    );
    if let Err(error) = written {
        let _ = writer.flush();
        return Err(error);
    }
    writer.flush()
}

/// The parts of one round's query report (`BlastpSections` with the round).
struct RoundSections<'a> {
    config: &'a PairwiseConfig,
    queries: &'a [BlastpPairwiseQuery],
    query_titles: &'a [Arc<[u8]>],
    report: &'a BlastpPairwiseReport,
    database: Option<&'a PairwiseDatabase>,
    parsed_subjects: &'a [Option<Arc<crate::report::db_subject::ParsedSubject>>],
    /// The sources of the others' Bioseqs (`pairwise::stored_source`).
    subject_sources: &'a [crate::report::defline::SourceTaxname],
    /// The loads of the report's subjects, by subject index (`BlastpSections::loads`).
    loads: &'a crate::report::db_subject::HitLoads,
    round: ReportRound<'a>,
    /// The sequences of the query-anchored section (outfmt 1 to 4); `None` for outfmt 0.
    anchored: Option<&'a dyn crate::report::query_anchored::ReportSequences>,
}

impl<'a> RoundSections<'a> {
    /// The `-html` inputs of the round's description tables (`PairwiseConfig::html_rows`).
    fn html_rows(&self) -> Option<HtmlRows<'a>> {
        self.config.html_rows().map(|(report, after)| HtmlRows {
            report,
            after,
            parsed: self.parsed_subjects,
            loads: self.loads,
        })
    }

    /// The subject of `s_idx` among the parsed subjects (`pairwise::parsed_subject`).
    fn parsed(&self, s_idx: u32) -> Option<&crate::report::db_subject::ParsedSubject> {
        self.parsed_subjects
            .get(s_idx as usize)
            .and_then(Option::as_deref)
            .filter(|subject| !subject.shown_as_record())
    }

    /// The description-table label of the subject `s_idx` of `group` (`BlastpSections::head`).
    fn label(&self, group: &QueryGroup<'_>, s_idx: u32) -> super::DescriptionLabel {
        let title = group
            .subject_titles
            .get(&s_idx)
            .map_or(&[][..], Vec::as_slice);
        database_description_label(
            title,
            super::stored_source(self.subject_sources, s_idx),
            self.parsed(s_idx),
            true,
            self.thrown(s_idx).is_some(),
        )
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:529-557
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:578-602
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:613
    /// The two-part table: the Seq-aligns (HSPs) in the report's order, each in the part of its
    /// subject, counted after it is put there, up to `-num_descriptions` (so the bound counts
    /// HSPs, and `-num_descriptions 0` keeps the first: `x_pairwise_par::two_part_lines`, which
    /// also gives the subjects whose description lines the table shows to the refusals and the
    /// loads decided before it, `x_pairwise_par::described_subjects`); each part is a table of its
    /// own over its HSPs (its own column widths), whose rows are its subjects; then the blank line
    /// of `x_DisplayDeflines`. `x_ConfigCShowBlastDefline` is called without the number of
    /// descriptions, so neither part leaves its header out.
    fn write_two_part_table<W: Write>(
        &self,
        writer: &mut W,
        group: &QueryGroup<'_>,
        found_before: &[bool],
        place: ReportHook<'_, W>,
    ) -> io::Result<()> {
        let [repeated, new] = crate::report::x_pairwise_par::two_part_lines(
            group.hsps(),
            self.report.num_descriptions,
            |hit: &&PairwiseHit| {
                found_before
                    .get(hit.hit.s_idx as usize)
                    .copied()
                    .unwrap_or(false)
            },
        );
        self.write_part(writer, group, &repeated, true, &mut *place)?;
        writeln!(writer)?;
        self.write_part(writer, group, &new, false, place)?;
        writeln!(writer)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:586-611
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:678-681
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:697-720
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:756-848
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1004-1011,1101,1134
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4261-4309
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:86-87
    /// One part: for an `-html` report, its own read of `./.ncbirc` (even when the part has no
    /// rows), then for the repeat pass the two lines of the Score/E header, then the part's line
    /// and a blank line, then its rows (`write_description_table`'s rows and widths, every row of
    /// the part counted: its `m_NumToShow` is the number of its Seq-aligns, and its `-html`
    /// `blast_rank` counts its rows from 1). A part without rows keeps the header's widths of
    /// "(Bits)" and "Value" (`x_InitDeflineTable`).
    fn write_part<W: Write>(
        &self,
        writer: &mut W,
        group: &QueryGroup<'_>,
        part: &[(u32, Vec<&PairwiseHit>)],
        repeat_pass: bool,
        place: ReportHook<'_, W>,
    ) -> io::Result<()> {
        let html = self.html_rows();
        if html.is_some() {
            place(
                ReportPlace::Html(crate::report::html_ids::HtmlPost::Ncbirc),
                writer,
            )?;
        }
        let rows: Vec<(u32, Option<&PairwiseHit>, f64)> = part
            .iter()
            .map(|(s_idx, hits)| {
                (
                    *s_idx,
                    shown_hsp(hits),
                    hits.iter().map(|hit| hit.hit.bit_score).sum(),
                )
            })
            .collect();
        // `write_description_table`'s widths (showdefline.cpp:1059-1142, cited there).
        let mut max_score = super::ncbi_k_bits().len();
        let mut max_evalue = "Value".len();
        let mut max_total = "Total".len();
        for (index, (_, best, total)) in rows.iter().enumerate() {
            let (bits, e_value) = table_scores(*best);
            max_score = max_score.max(format_bitscore_ncbi(bits).len());
            max_evalue = max_evalue.max(format_evalue_ncbi(e_value).len());
            let total_len = format_bitscore_ncbi(*total).len();
            if index + 1 == rows.len() {
                if max_total < total_len {
                    max_score = total_len;
                }
            } else {
                max_total = max_total.max(total_len);
            }
        }
        if repeat_pass {
            write_spaces(writer, 70)?;
            writeln!(writer, "{:<max_score$}    E", "Score")?;
            write!(
                writer,
                "{:<69}",
                "Sequences producing significant alignments:"
            )?;
            writeln!(writer, "{:<max_score$}  Value", super::ncbi_k_bits())?;
            writeln!(writer, "Sequences used in model and found again:")?;
        } else {
            writeln!(
                writer,
                "Sequences not found previously or not previously below threshold:"
            )?;
        }
        writeln!(writer)?;
        let end = html.and_then(|html| html.end(&rows.iter().map(|row| row.0).collect::<Vec<_>>()));
        for (index, (s_idx, best, _)) in rows.iter().enumerate() {
            // Each line loads its subject (`write_description_table`).
            place(ReportPlace::Load(*s_idx), writer)?;
            crate::report::x_pairwise_par::end_here(end, index)?;
            let (bits, e_value) = table_scores(*best);
            let bits = format_bitscore_ncbi(bits);
            write_description_row_start(
                writer,
                &self.label(group, *s_idx),
                (*s_idx, index + 1),
                &bits,
                max_score,
                html,
                &mut *place,
            )?;
            writeln!(writer, "  {:<max_evalue$}", format_evalue_ncbi(e_value))?;
        }
        Ok(())
    }
}

impl QuerySections for RoundSections<'_> {
    fn num_alignments(&self) -> usize {
        self.report.num_alignments
    }

    fn line_length(&self) -> usize {
        self.config.line_length
    }

    fn sort(&self) -> (crate::report::hit_sort::ReportSort, bool) {
        (self.config.sort, false)
    }

    /// `BlastpSections::joins`.
    fn joins(&self, previous: u32, s_idx: u32) -> bool {
        self.database.is_some()
            && crate::report::db_subject::rank_joins_as_copy(
                self.parsed_subjects,
                previous as usize,
                s_idx as usize,
            )
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1474-1476
    /// `Results from round N`, `\n\n`, the query's header, then "No hits found" and the footer,
    /// or the description table (in two parts after a round whose inclusion set is not empty)
    /// and the blank line after it.
    fn head<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: Option<&QueryGroup<'_>>,
        place: ReportHook<'_, W>,
    ) -> io::Result<()> {
        let (query, report) = (&self.queries[q_idx], self.report);
        writeln!(writer, "Results from round {}", self.round.number)?;
        writeln!(writer)?;
        writeln!(writer)?;
        let query_title = self
            .query_titles
            .get(q_idx)
            .map_or(&b""[..], |title| title.as_ref());
        write_blastp_query_header(
            writer,
            query_title,
            query.query_length,
            self.config.html.is_some(),
        )?;
        let Some(group) = group else {
            write_no_hits_found(writer)?;
            if query.batch_skipped {
                write_tblastn_unsearched_query_footer_spacing(writer, false)?;
            } else if !query.valid {
                // `x_PrintOneQueryFooter`'s first `NcbiEndl` flushes the stream.
                crate::report::ostream::ncbi_endl(writer)?;
                writeln!(writer)?;
                writeln!(writer)?;
                writeln!(writer, "Effective search space used: 0")?;
            } else {
                write_blastp_query_footer_spacing(
                    writer,
                    query.ungapped_karlin,
                    report.gapped_karlin,
                    report.gumbel,
                    report.gapped,
                    query.effective_search_space,
                    false,
                )?;
            }
            return Ok(());
        };
        match self.round.found_before {
            Some(found_before) => self.write_two_part_table(writer, group, found_before, place)?,
            None => write_description_table(
                writer,
                &group.subject_order[..group.subject_order.len().min(report.num_descriptions)],
                group.subject_order.len() <= report.num_descriptions,
                &group.subject_hits,
                false,
                // psiblast does not hand `-sorthits` to the formatter (psiblast_app.cpp:502).
                None,
                |s_idx| self.label(group, s_idx),
                place,
                self.html_rows(),
                None,
            )?,
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1540-1545
        writeln!(writer)
    }

    /// The alignments of the subject `s_idx` (`BlastpSections::subject`, without mask spans:
    /// NCBI's psiblast shows no `<font>` mask spans in any round).
    fn subject<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        probe: Option<&mut FormatProbe<'_>>,
        group: &QueryGroup<'_>,
        s_idx: u32,
        part: std::ops::Range<usize>,
        place: ReportHook<'_, W>,
    ) -> io::Result<()> {
        write_blastp_subject_part(
            writer,
            probe,
            self.config,
            self.database.is_some(),
            self.parsed_subjects,
            self.subject_sources,
            self.loads,
            group,
            s_idx,
            part,
            HtmlRowMasks::new(None, 0),
            place,
        )
    }

    fn html_links(&self) -> Option<&crate::report::html_ids::HtmlLinks> {
        self.config.html_links()
    }

    fn alignment_id_is_local(&self, s_idx: u32) -> bool {
        alignment_id_local(self.parsed_subjects, s_idx)
    }

    fn view(&self) -> PairwiseView {
        self.config.view
    }

    /// `BlastpSections::thrown`.
    fn thrown(&self, s_idx: u32) -> Option<(&crate::ncbi_exception::NcbiException, bool)> {
        super::thrown_load(self.loads, self.parsed_subjects, s_idx)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1550-1552
    /// The query-anchored section of the query's aligned subjects (`BlastpSections::alignments`);
    /// nothing for outfmt 0, whose subjects' alignments `subject` writes.
    fn alignments<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: &QueryGroup<'_>,
    ) -> io::Result<()> {
        if self.anchored.is_none() {
            return Ok(());
        }
        super::write_anchored_section(writer, self.config, self.anchored, q_idx, |sequences| {
            super::blastp_section_hsps(
                self.report.num_alignments,
                group,
                self.parsed_subjects,
                |s_idx| sequences.subject_is_query(q_idx, s_idx),
            )
        })
    }

    /// The footer of a query with hits (`BlastpSections::foot`).
    fn foot<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
    ) -> io::Result<()> {
        let query = &self.queries[q_idx];
        write_blastp_query_footer_spacing(
            writer,
            query.ungapped_karlin,
            self.report.gapped_karlin,
            self.report.gumbel,
            self.report.gapped,
            query.effective_search_space,
            false,
        )
    }
}

/// psiblast's outfmt 0 epilog after the last round (`PrintEpilog`): the blank lines and the
/// footer of `write_blastp_pairwise_report`'s epilog, with `</PRE>\n</BODY>\n</HTML>\n` when
/// `html` (`-html`) is set.
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:561-562
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2281-2292
pub(crate) fn write_psiblast_epilog(
    writer: &mut dyn ReportWrite,
    report: &BlastpPairwiseReport,
    database: Option<&PairwiseDatabase>,
    html: bool,
) -> io::Result<()> {
    let mut buffered = io::BufWriter::new(RoundWriter(writer));
    let writer = &mut buffered;
    // The blank lines before the epilog (blast_format.cpp:2249), each a flush.
    crate::report::ostream::ncbi_endl(writer)?;
    crate::report::ostream::ncbi_endl(writer)?;
    super::write_blastp_final_footer(writer, report, database, html)?;
    writer.flush()
}
