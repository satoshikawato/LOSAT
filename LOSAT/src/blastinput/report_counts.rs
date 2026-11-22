//! The hit list size and the numbers of descriptions and alignments of a run
//! (`CFormattingArgs::ExtractAlgorithmOptions`), shared by the five programs.
//!
//! Every program calls [`report_counts`] once, before the search, and keeps the result.
//! The warnings and the error that NCBI raises with the numbers are [`ReportCounts::warnings`]
//! and [`ReportCounts::validate`].

use crate::blastinput::app::{options_error, output_format};

/// NCBI's default hit list size.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:160
pub const BLAST_HITLIST_SIZE: usize = 500;

/// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:219,221
pub const DEFAULT_NUM_DESCRIPTIONS: usize = 500;
/// See [`DEFAULT_NUM_DESCRIPTIONS`].
pub const DEFAULT_NUM_ALIGNMENTS: usize = 250;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2683-2722
/// The formatting handler's options besides the counts, as the command line gives them:
/// `-line_length`, `-sorthits` and `-sorthsps` (`CFormattingArgs`), one definition that every
/// program's parser flattens in (BLASTN, BLASTP, BLASTX, TBLASTN, TBLASTX, psiblast and
/// rpsblast); the default is none given. The five programs hand all three to the formatter
/// (`blastn_app.cpp:251-253` and the others' alike); psiblast and rpsblast only the line length
/// (`psiblast_app.cpp:502`, `rpsblast_app.cpp:187`): their `report_sort` is none.
#[derive(clap::Args, Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FormattingArgs {
    /// `-line_length`: the width of the pairwise alignment rows (60 when it is not given).
    #[arg(long = "line_length", value_name = "line_length", value_parser = formatting_line_length, help = "Line length for formatting alignments; Not applicable for outfmt > 4 [default: 60]")]
    pub line_length: Option<usize>,
    /// `-sorthits` and `-sorthsps` (0 to 4).
    #[arg(long = "sorthits", value_name = "sort_hits", value_parser = formatting_sort_option, help = "Sorting option for hits: alignment view options: 0 = Sort by evalue, 1 = Sort by bit score, 2 = Sort by total score, 3 = Sort by percent identity, 4 = Sort by query coverage; Not applicable for outfmt > 4")]
    pub sort_hits: Option<u8>,
    #[arg(long = "sorthsps", value_name = "sort_hsps", value_parser = formatting_sort_option, help = "Sorting option for hps: 0 = Sort by hsp evalue, 1 = Sort by hsp score, 2 = Sort by hsp query start, 3 = Sort by hsp percent identity, 4 = Sort by hsp subject start; Not applicable for outfmt != 0")]
    pub sort_hsps: Option<u8>,
}

impl FormattingArgs {
    /// The width of the pairwise alignment rows (`m_LineLength`): `-line_length`, or NCBI's
    /// `kDfltLineLength` (`report::pairwise::DEFAULT_LINE_LENGTH`).
    pub fn line_length(&self) -> usize {
        self.line_length
            .unwrap_or(crate::report::pairwise::DEFAULT_LINE_LENGTH)
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:251-253
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2942-2946
    /// `-sorthits` and `-sorthsps` as one of the five search programs' pairwise reports reads
    /// them (blastp_app.cpp:246-248, blastx_app.cpp:247-249, tblastn_app.cpp:263-265 and
    /// tblastx_app.cpp:166-168 hand them over as blastn_app.cpp does). `m_HitsSortOption` is read
    /// for the output formats up to 4 and `m_HspsSortOption` for format 0 only: formats 0 to 4 sort
    /// the hits, format 0 alone the HSPs (`PairwiseView::report_sort`; `ReportCounts::warnings_with`
    /// warns about the others).
    pub fn report_sort(&self) -> crate::report::hit_sort::ReportSort {
        crate::report::hit_sort::ReportSort {
            hits: self.sort_hits,
            hsps: self.sort_hsps,
        }
    }
}

/// A `-sorthits` or `-sorthsps` value: an integer from 0 to 4 (`CArgAllowValuesBetween(0, 4,
/// true)`, cited at `FormattingArgs`).
fn formatting_sort_option(value: &str) -> Result<u8, String> {
    let option = crate::blastinput::value_parsers::ncbi_integer_at_least(value, 0)?;
    u8::try_from(option)
        .ok()
        .filter(|option| *option <= 4)
        .ok_or_else(|| "Illegal value: expected an integer from 0 to 4".to_string())
}

/// A `-line_length` value: an integer of 1 or more (`CArgAllowValuesGreaterThanOrEqual(1)`,
/// cited at `FormattingArgs`).
fn formatting_line_length(value: &str) -> Result<usize, String> {
    let length = crate::blastinput::value_parsers::ncbi_integer_at_least(value, 1)?;
    usize::try_from(length).map_err(|_| "Illegal value: expected an integer >= 1".to_string())
}

/// The numbers that one output format of a run uses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReportCounts {
    /// The hit list size of the search (`opt.SetHitlistSize`).
    pub hitlist_size: usize,
    /// The number of one-line descriptions of the pairwise format (`m_NumDescriptions`).
    pub num_descriptions: usize,
    /// The number of alignments of the pairwise format (`m_NumAlignments`).
    pub num_alignments: usize,
}

impl ReportCounts {
    /// The counts of what the pairwise reports show of each query (`db_subject::Shown`):
    /// `num_descriptions` description lines, the alignments of the first `num_alignments` hits.
    pub fn shown(self) -> crate::report::db_subject::Shown {
        crate::report::db_subject::Shown {
            descriptions: self.num_descriptions,
            alignments: self.num_alignments,
        }
    }
}

impl Default for ReportCounts {
    /// Neither `-max_target_seqs`, `-num_descriptions` nor `-num_alignments`, outfmt 0.
    fn default() -> Self {
        Self {
            hitlist_size: BLAST_HITLIST_SIZE,
            num_descriptions: DEFAULT_NUM_DESCRIPTIONS,
            num_alignments: DEFAULT_NUM_ALIGNMENTS,
        }
    }
}

/// The numbers of a run with output format `outfmt`, given `-max_target_seqs` (`mts`),
/// `-num_descriptions` (`nd`) and `-num_alignments` (`na`) as the command line gives them.
/// `-max_target_seqs` excludes the other two (the argument parsers refuse both).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2894-2978
/// (`BLAST_HITLIST_SIZE` is 500, `m_DfltNumDescriptions` 500, `m_DfltNumAlignments` 250; the
/// import search strategy branch at 2900-2908 does not exist in LOSAT)
/// Note that `args[kArgNumDescriptions]` is false for an absent option, not for the value 0.
pub fn report_counts(
    outfmt: i32,
    max_target_seqs: Option<usize>,
    num_descriptions: Option<usize>,
    num_alignments: Option<usize>,
) -> ReportCounts {
    let mut hitlist_size = BLAST_HITLIST_SIZE;
    if outfmt <= output_format::FLAT_QUERY_ANCHORED_NO_IDENTITIES {
        let mut counts = ReportCounts::default();
        if let Some(n) = num_descriptions {
            counts.num_descriptions = n;
        }
        if let Some(n) = num_alignments {
            counts.num_alignments = n;
        }
        if let Some(n) = max_target_seqs {
            counts.num_descriptions = n;
            counts.num_alignments = n;
            hitlist_size = n;
        }
        if num_descriptions.is_some() || num_alignments.is_some() {
            hitlist_size = counts.num_descriptions.max(counts.num_alignments);
        }
        counts.hitlist_size = hitlist_size;
        counts
    } else {
        if let Some(n) = max_target_seqs {
            hitlist_size = n;
        } else if let Some(n) = num_alignments {
            hitlist_size = n;
        }
        ReportCounts {
            hitlist_size,
            num_descriptions: hitlist_size,
            num_alignments: hitlist_size,
        }
    }
}

impl ReportCounts {
    /// The warnings of the formatting handler for this run, in NCBI's order: the ignored
    /// `-num_descriptions` of an output format above 4, then the few matches. They go to
    /// standard error once per run, before any warning of a query, and not into
    /// `QueryWarnings` (which the XML format writes into its iterations). A run with
    /// `-line_length`, `-sorthits` or `-sorthsps` takes `warnings_with`.
    pub fn warnings(&self, program: &str, outfmt: i32, num_descriptions_given: bool) -> Vec<u8> {
        self.warnings_with(
            program,
            outfmt,
            num_descriptions_given,
            &FormattingArgs::default(),
        )
    }

    /// `warnings` with the formatting options `formatting`, in NCBI's order: for an
    /// output format above 4 the ignored `-num_descriptions`, `-line_length` and `-sorthits`,
    /// then the few matches, then the ignored `-sorthsps` of a format other than 0.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2947-2989
    pub fn warnings_with(
        &self,
        program: &str,
        outfmt: i32,
        num_descriptions_given: bool,
        formatting: &FormattingArgs,
    ) -> Vec<u8> {
        let mut text = Vec::new();
        let mut warn = |message: &str| {
            text.extend_from_slice(format!("Warning: [{program}] {message}\n").as_bytes());
        };
        if outfmt > output_format::FLAT_QUERY_ANCHORED_NO_IDENTITIES {
            if num_descriptions_given {
                warn("The parameter -num_descriptions is ignored for output formats > 4 . Use -max_target_seqs to control output");
            }
            if formatting.line_length.is_some() {
                warn("The parameter -line_length is not applicable for output formats > 4 .");
            }
            if formatting.sort_hits.is_some() {
                warn("The parameter -sorthits is ignored for output formats > 4.");
            }
        }
        if self.hitlist_size < 5 {
            text.extend_from_slice(&crate::report::query_warnings::few_matches_warning(program));
        }
        if formatting.sort_hsps.is_some() && outfmt != output_format::PAIRWISE {
            text.extend_from_slice(
                format!(
                    "Warning: [{program}] The parameter -sorthsps is ignored for output formats != 0.\n"
                )
                .as_bytes(),
            );
        }
        text
    }

    /// NCBI's error for a hit list size of 0, which only `-num_descriptions` and
    /// `-num_alignments` can give (`-max_target_seqs` is 1 or more). It is raised by the
    /// validation of the options, with the hit saving options' checks after those of the
    /// scoring, lookup table and word options and before the e-value's.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1511-1517
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.hitlist_size < 1 {
            return Err(options_error("No hits are being saved"));
        }
        Ok(())
    }
}

/// The numbers of a run that writes the output formats `outfmts` (one for the command
/// line), and the warnings of the formatting handler, which belong to the first format.
///
/// The search keeps one hit list per run, so output formats whose hit list sizes differ
/// (possible only with `-num_descriptions` / `-num_alignments`, for example the formats 0
/// and 6 with `-num_descriptions 3`: 250 and 500) cannot be written from one run and are
/// refused, never resolved for one of them. Formats of one hit list size may differ in
/// `num_descriptions` and `num_alignments` (format 0: 500 and 250; the tabular formats: the
/// hit list size), which only the pairwise format reads: the result carries those of the
/// first pairwise format (the first format when there is none).
pub fn run_report_counts(
    program: &str,
    outfmts: &[i32],
    max_target_seqs: Option<usize>,
    num_descriptions: Option<usize>,
    num_alignments: Option<usize>,
) -> anyhow::Result<(ReportCounts, Vec<u8>)> {
    run_report_counts_with(
        program,
        outfmts,
        max_target_seqs,
        num_descriptions,
        num_alignments,
        &FormattingArgs::default(),
    )
}

/// `run_report_counts` of a run with the formatting options `formatting` (`-line_length`,
/// `-sorthits`, `-sorthsps`), whose warnings come with the counts' in NCBI's order
/// (`ReportCounts::warnings_with`).
pub fn run_report_counts_with(
    program: &str,
    outfmts: &[i32],
    max_target_seqs: Option<usize>,
    num_descriptions: Option<usize>,
    num_alignments: Option<usize>,
    formatting: &FormattingArgs,
) -> anyhow::Result<(ReportCounts, Vec<u8>)> {
    let formats: Vec<ReportCounts> = outfmts
        .iter()
        .map(|&outfmt| report_counts(outfmt, max_target_seqs, num_descriptions, num_alignments))
        .collect();
    let counts = common_report_counts(outfmts, &formats)?;
    let warnings =
        formats[0].warnings_with(program, outfmts[0], num_descriptions.is_some(), formatting);
    Ok((counts, warnings))
}

/// The counts that the formats `outfmts` (with the counts `formats`, one each) share: see
/// `run_report_counts`.
pub fn common_report_counts(
    outfmts: &[i32],
    formats: &[ReportCounts],
) -> anyhow::Result<ReportCounts> {
    let first = *formats
        .first()
        .ok_or_else(|| anyhow::anyhow!("no output format was requested"))?;
    if formats
        .iter()
        .any(|counts| counts.hitlist_size != first.hitlist_size)
    {
        anyhow::bail!(
            "the requested output formats need different hit list sizes (-num_descriptions and -num_alignments are read per output format); LOSAT searches once per run, so request them in separate runs"
        );
    }
    Ok(outfmts
        .iter()
        .position(|&outfmt| outfmt <= output_format::FLAT_QUERY_ANCHORED_NO_IDENTITIES)
        .map_or(first, |index| formats[index]))
}

#[cfg(test)]
#[path = "report_counts_tests.rs"]
mod tests;
