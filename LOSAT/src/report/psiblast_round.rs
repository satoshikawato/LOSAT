//! psiblast's outfmt 0 report of one round of a FASTA query (series PS, item PS5a): the report of
//! `write_blastp_pairwise_report` with the round's number in its head and, after the first round,
//! the description table in two parts, the subjects used in the PSSM's model and found again and
//! those found for the first time (`CBlastFormat::x_DisplayDeflines` with the previous round's
//! inclusion set). The round driver (`blastp::blast_engine::psi_round`) writes the epilog after the
//! last round; this writer never does.
//!
//! The writer is a child of `pairwise` and uses its parts. It writes through a type of its own
//! (`RoundWriter`), so that every generic part it calls is an instance of its own and the default
//! report's instances are not shared (the copy rule of DESIGN_ps4 §4.3), and formats the one query
//! of a round serially (no `LOSAT_X_PAIRWISEPAR` pieces).

use super::{
    database_description_label, write_alignment, write_blastp_pairwise_prolog,
    write_blastp_query_footer_spacing, write_blastp_query_header,
    write_database_subject_alignments_part, write_description_table, write_hsp_info,
    write_no_hits_found, write_spaces, write_subject_alignments_part,
    write_tblastn_unsearched_query_footer_spacing, BlastpPairwiseQuery, BlastpPairwiseReport,
    DatabaseHeading, PairwiseConfig, PairwiseDatabase, PairwiseHit,
};
use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::report::outfmt6::{format_bitscore_ncbi, format_evalue_ncbi};
use crate::report::query_warnings::{no_places, ReportHook, ReportPlace};
use crate::report::x_pairwise_par::QueryGroup;
use std::io::{self, Write};
use std::sync::Arc;

/// The report's writer: the caller's stream behind a type that no default report writes through.
struct RoundWriter<'w>(&'w mut dyn Write);

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
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:278-286
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_psiblast_round_report(
    hits: &[PairwiseHit],
    query_range: std::ops::Range<usize>,
    writer: &mut dyn Write,
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
    mut warnings: Option<&mut super::super::query_warnings::QueryWarnings<'_>>,
    prolog: bool,
    round: ReportRound<'_>,
) -> io::Result<()> {
    let mut buffered = io::BufWriter::new(RoundWriter(writer));
    let writer = &mut buffered;
    if prolog && query_range.start == 0 {
        write_blastp_pairwise_prolog(
            writer,
            &report.version,
            &config.program,
            &report.database_name,
            report.database_num_sequences,
            report.database_total_letters,
            database,
            super::html_with_diagnostics(config.html_links(), warnings.as_deref_mut())?,
        )?;
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
    };
    for (q_idx, query_hits) in query_range.zip(&hits_by_query) {
        // The query's warnings come before its preamble, after the report written so far
        // (`QueryWarnings::before_query`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, writer)?;
        }
        let group = QueryGroup::new(query_hits, subject_titles, |previous, s_idx| {
            sections.joins(previous, s_idx)
        });
        // NCBI loads each subject that the query's batch did not prefetch at its description
        // line, or one the table does not show at its alignment's heading
        // (`QueryWarnings::subject_read`).
        let mut place = |place: ReportPlace, writer: &mut io::BufWriter<RoundWriter<'_>>| {
            crate::report::x_pairwise_par::at_place(
                warnings.as_deref_mut(),
                None,
                q_idx,
                place,
                writer,
            )
        };
        sections.head(writer, q_idx, group.as_ref(), &mut place)?;
        if let Some(group) = &group {
            // An aligned subject's load may end the run where the view's answer says
            // (`PairwiseView::after_throw`, as `x_pairwise_par::alignment_end` asks it).
            let aligned = group.aligned(report.num_alignments);
            let end = crate::report::x_pairwise_par::aligned_end(
                config.view.after_throw(),
                aligned,
                |s_idx| sections.thrown(s_idx),
            );
            for &s_idx in &aligned[..end.map_or(aligned.len(), |(at, _)| at)] {
                place(ReportPlace::Load(s_idx), &mut *writer)?;
                sections.subject(writer, probe.as_deref_mut(), group, s_idx)?;
            }
            if let Some((_, exception)) = end {
                if let Some(warnings) = warnings.as_deref_mut() {
                    warnings.end_run();
                }
                writer.flush()?;
                return Err(exception.clone().into_io_error());
            }
            sections.foot(writer, q_idx)?;
        }
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
}

impl RoundSections<'_> {
    /// `BlastpSections::joins`.
    fn joins(&self, previous: u32, s_idx: u32) -> bool {
        self.database.is_some()
            && crate::report::db_subject::rank_joins_as_copy(
                self.parsed_subjects,
                previous as usize,
                s_idx as usize,
            )
    }

    /// `BlastpSections::thrown`.
    fn thrown(&self, s_idx: u32) -> Option<(&crate::ncbi_exception::NcbiException, bool)> {
        super::thrown_load(self.loads, self.parsed_subjects, s_idx)
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1474-1476
    /// `Results from round N`, `\n\n`, the query's header, then "No hits found" and the footer,
    /// or the description table (in two parts after a round whose inclusion set is not empty)
    /// and the blank line after it.
    fn head<W: Write>(
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
        // psiblast refuses `-html` (`cli::UNPORTED_PSIBLAST_OWN_ARGS`).
        write_blastp_query_header(writer, query_title, query.query_length, false)?;
        let Some(group) = group else {
            write_no_hits_found(writer)?;
            if query.batch_skipped {
                write_tblastn_unsearched_query_footer_spacing(writer, false)?;
            } else if !query.valid {
                writeln!(writer)?;
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
                // psiblast refuses `-html` (`cli::UNPORTED_PSIBLAST_OWN_ARGS`): the plain table.
                None,
                None,
            )?,
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1520-1528
        writeln!(writer)
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

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:756-848
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:85-87
    /// One part: for the repeat pass the two lines of the Score/E header, then the part's line
    /// and a blank line, then its rows (`write_description_table`'s rows and widths, every row of
    /// the part counted: its `m_NumToShow` is the number of its Seq-aligns). A part without rows
    /// keeps the header's widths of "(Bits)" and "Value" (`x_InitDeflineTable`).
    fn write_part<W: Write>(
        &self,
        writer: &mut W,
        group: &QueryGroup<'_>,
        part: &[(u32, Vec<&PairwiseHit>)],
        repeat_pass: bool,
        place: ReportHook<'_, W>,
    ) -> io::Result<()> {
        let rows: Vec<(u32, &PairwiseHit, f64)> = part
            .iter()
            .map(|(s_idx, hits)| {
                let mut best = hits[0];
                for hit in hits {
                    if hit.hit.bit_score > best.hit.bit_score {
                        best = hit;
                    }
                }
                (*s_idx, best, hits.iter().map(|hit| hit.hit.bit_score).sum())
            })
            .collect();
        // `write_description_table`'s widths (showdefline.cpp:1059-1142, cited there).
        let mut max_score = super::ncbi_k_bits().len();
        let mut max_evalue = "Value".len();
        let mut max_total = "Total".len();
        for (index, (_, best, total)) in rows.iter().enumerate() {
            max_score = max_score.max(format_bitscore_ncbi(best.hit.bit_score).len());
            max_evalue = max_evalue.max(format_evalue_ncbi(best.hit.e_value).len());
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
        for (s_idx, best, _) in &rows {
            // Each line loads its subject (`write_description_table`).
            place(ReportPlace::Load(*s_idx), writer)?;
            self.label(group, *s_idx).write_cut(writer)?;
            writeln!(
                writer,
                "  {:<max_score$}  {:<max_evalue$}",
                format_bitscore_ncbi(best.hit.bit_score),
                format_evalue_ncbi(best.hit.e_value)
            )?;
        }
        Ok(())
    }

    /// The alignments of the subject `s_idx` (`BlastpSections::subject`, all its HSPs).
    fn subject<W: Write>(
        &self,
        writer: &mut W,
        probe: Option<&mut FormatProbe<'_>>,
        group: &QueryGroup<'_>,
        s_idx: u32,
    ) -> io::Result<()> {
        let config = self.config;
        let hsps = group.subject_hits[&s_idx].len();
        if let Some(id) = config.unreadable_query_id(self.parsed_subjects, s_idx) {
            return write_database_subject_alignments_part(
                writer,
                probe,
                DatabaseHeading::Skipped(&id),
                true,
                &group.subject_hits[&s_idx],
                &group.subject_hit_indices[&s_idx],
                0..hsps,
                None,
                &mut no_places,
                |_, _, _| unreachable!("a skipped alignment is written by its line"),
            );
        }
        let parsed = self.parsed(s_idx);
        if self.database.is_some() || parsed.is_some() {
            let heading = match parsed {
                Some(subject) => DatabaseHeading::loaded(
                    subject,
                    true,
                    group.own_hsps(s_idx),
                    self.thrown(s_idx).is_some(),
                ),
                None => DatabaseHeading::Title(
                    &group.subject_titles[&s_idx],
                    super::stored_source(self.subject_sources, s_idx),
                ),
            };
            return write_database_subject_alignments_part(
                writer,
                probe,
                heading,
                true,
                &group.subject_hits[&s_idx],
                &group.subject_hit_indices[&s_idx],
                0..hsps,
                None,
                &mut no_places,
                |writer, hit, _| {
                    write_hsp_info(writer, hit, config)?;
                    // psiblast refuses `-html`: the rows have no spans.
                    write_alignment(writer, hit, None, config)
                },
            );
        }
        write_subject_alignments_part(
            writer,
            probe,
            s_idx,
            self.parsed_subjects,
            &group.subject_titles[&s_idx],
            true,
            &group.subject_hits[&s_idx],
            &group.subject_hit_indices[&s_idx],
            0..hsps,
            None,
            &mut no_places,
            |writer, hit, _| {
                write_hsp_info(writer, hit, config)?;
                write_alignment(writer, hit, None, config)
            },
        )
    }

    /// The footer of a query with hits (`BlastpSections::foot`).
    fn foot<W: Write>(&self, writer: &mut W, q_idx: usize) -> io::Result<()> {
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
/// footer of `write_blastp_pairwise_report`'s epilog.
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:561-562
pub(crate) fn write_psiblast_epilog(
    writer: &mut dyn Write,
    report: &BlastpPairwiseReport,
    database: Option<&PairwiseDatabase>,
) -> io::Result<()> {
    let mut buffered = io::BufWriter::new(RoundWriter(writer));
    let writer = &mut buffered;
    writeln!(writer)?;
    writeln!(writer)?;
    // psiblast refuses `-html` (`cli::UNPORTED_PSIBLAST_OWN_ARGS`): the plain epilog.
    super::write_blastp_final_footer(writer, report, database, false)?;
    writer.flush()
}
