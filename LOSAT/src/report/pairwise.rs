//! Pairwise alignment output (outfmt 0)
//!
//! This module implements the traditional BLAST pairwise alignment output format.
//!
//! Reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp

use super::outfmt6::{format_bitscore_ncbi, format_evalue_ncbi, ReportContext};
use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::common::Hit;
use crate::config::ScoringMatrix;
use crate::stats::{BlastGumbelBlk, KarlinParams};
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3563-3581
// ```c++
//     retval.Resize(k_NumAsciiChar, k_NumAsciiChar, -1000);
//
//     SNCBIFullScoreMatrix mtx;
//     NCBISM_Unpack(packed_mtx, &mtx);
//
//     for(int i = 0; i < ePMatrixSize; ++i){
//         for(int j = 0; j < ePMatrixSize; ++j){
//             retval((size_t)k_PSymbol[i], (size_t)k_PSymbol[j]) =
//                 mtx.s[(size_t)k_PSymbol[i]][(size_t)k_PSymbol[j]];
//         }
//     }
//     for(int i = 0; i < ePMatrixSize; ++i) {
//         retval((size_t)k_PSymbol[i], '*') = retval('*',(size_t)k_PSymbol[i]) = -4;
//     }
//     retval('*', '*') = 1;
//     // this is to count Selenocysteine to Cysteine matches as positive
//     retval('U', 'U') = retval('C', 'C');
//     retval('U', 'C') = retval('C', 'C');
//     retval('C', 'U') = retval('C', 'C');
// ```
use crate::utils::matrix::protein_display_score;
use std::io::{self, Write};
use std::sync::Arc;

// =============================================================================
// NCBI Pairwise Output Format (outfmt 0)
// =============================================================================

/// Line length for alignment display (NCBI default: 60)
pub const DEFAULT_LINE_LENGTH: usize = 60;

/// Configuration for pairwise output
#[derive(Debug, Clone)]
pub struct PairwiseConfig {
    /// Line length for sequence display
    pub line_length: usize,
    /// Show GI numbers if available
    pub show_gi: bool,
    /// Show frame information (for translated searches)
    pub show_frame: bool,
    /// Program type for proper formatting
    pub program: String,
    /// Protein scoring matrix for positives/midline rendering.
    pub protein_matrix: ScoringMatrix,
}

impl Default for PairwiseConfig {
    fn default() -> Self {
        Self {
            line_length: DEFAULT_LINE_LENGTH,
            show_gi: false,
            show_frame: true,
            program: "tblastx".to_string(),
            protein_matrix: ScoringMatrix::Blosum62,
        }
    }
}

/// Extended hit information for pairwise display
///
/// For full pairwise output, we need the actual aligned sequences.
/// This struct extends Hit with optional sequence data.
#[derive(Debug, Clone)]
pub struct PairwiseHit {
    /// Base hit information
    pub hit: Hit,
    /// Aligned query sequence (if available)
    pub query_seq: Option<String>,
    /// Aligned subject sequence (if available)
    pub subject_seq: Option<String>,
    /// Query frame (for translated searches)
    pub query_frame: Option<i8>,
    /// Subject frame (for translated searches)
    pub subject_frame: Option<i8>,
    /// Number of positive matches (for protein/translated)
    pub positives: Option<usize>,
    /// Number of gaps
    pub gaps: Option<usize>,
    /// Subject sequence length
    pub subject_length: Option<usize>,
    /// Subject description/title
    pub subject_title: Option<String>,
    // NCBI c++/src/algo/blast/core/blast_kappa.c:331-342;
    // c++/src/objtools/align_format/showalign.cpp:3599-3604:
    // eCompoScaleOldMatrix => Method 1; other adjusted matrices => Method 2.
    pub comp_adjust_method: Option<u8>,
    // NCBI c++/src/objtools/align_format/showalign.cpp:3595-3598:
    // if (aln_vec_info->sum_n > 0) out << "(" << sum_n << ")";
    pub sum_n: Option<i32>,
}

impl From<Hit> for PairwiseHit {
    fn from(hit: Hit) -> Self {
        let positives = hit.num_positives;
        let gaps = hit.gap_letters();
        Self {
            hit,
            query_seq: None,
            subject_seq: None,
            query_frame: None,
            subject_frame: None,
            positives: Some(positives),
            gaps: Some(gaps),
            subject_length: None,
            subject_title: None,
            comp_adjust_method: None,
            sum_n: None,
        }
    }
}

/// NCBI BLASTP pairwise per-query ancillary data.
///
/// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_results.cpp:72-115
/// ```c
/// CBlastAncillaryData::CBlastAncillaryData(EBlastProgramType program_type,
///                     int query_number,
///                     const BlastScoreBlk *sbp,
///                     const BlastQueryInfo *query_info)
/// {
///     ...
///     m_SearchSpace = ctx->eff_searchsp;
///     ...
///     s_InitializeKarlinBlk(sbp->kbp_std[ctx_index], &m_UngappedKarlinBlk);
/// }
/// ```
#[derive(Debug, Clone)]
pub struct BlastpPairwiseQuery {
    pub query_name: String,
    pub query_length: usize,
    pub ungapped_karlin: KarlinParams,
    pub effective_search_space: i64,
    /// Whether the query's Karlin-Altschul parameters could be computed (`is_valid`).
    pub valid: bool,
    /// Whether every query of the query's batch is invalid, so NCBI did not search the
    /// batch (local_blast.cpp:177-207: -1 parameters and no Seq-align set).
    pub batch_skipped: bool,
}

/// NCBI BLASTP pairwise run-level footer/header data.
///
/// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:2248-2285
/// ```c
/// if ( !m_IsBl2Seq || m_IsDbScan) {
///     CBlastFormatUtil::PrintDbReport(m_DbInfo, kFormatLineLength,
///                                     m_Outfile, false);
/// }
/// ...
/// m_Outfile << "\n\nMatrix: " << options.GetMatrixName() << "\n";
/// ...
/// m_Outfile << "Neighboring words threshold: "
///           << options.GetWordThreshold() << "\n";
/// m_Outfile << "Window for multiple hits: "
///           << options.GetWindowSize() << "\n";
/// ```
#[derive(Debug, Clone)]
pub struct BlastpPairwiseReport {
    pub version: String,
    pub database_name: String,
    pub database_num_sequences: usize,
    pub database_total_letters: usize,
    pub matrix_name: String,
    pub gap_open: i32,
    pub gap_extend: i32,
    /// `GetWordThreshold()`, a double.
    pub word_threshold: f64,
    pub window_size: i32,
    pub gapped_karlin: KarlinParams,
    pub gumbel: BlastGumbelBlk,
    /// Subjects in the description table and with alignments, per query (NCBI's
    /// `m_NumDescriptions` and `m_NumAlignments`).
    pub num_descriptions: usize,
    pub num_alignments: usize,
}

// =============================================================================
// Pairwise Output Writers
// =============================================================================

/// The defline of a subject as `bio` reads it: the ID, a space and the rest.
fn subject_defline(subject_id: &str, subject_title: Option<&str>) -> String {
    match subject_title.filter(|title| !title.is_empty()) {
        Some(title) => format!("{subject_id} {title}"),
        None => subject_id.to_string(),
    }
}

/// Write database/subject information header
///
/// Reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp
///
/// Format:
/// ```text
/// >subject_id subject_description
/// Length=XXXX
/// ```
pub fn write_subject_header<W: Write>(
    writer: &mut W,
    subject_id: &str,
    subject_title: Option<&str>,
    subject_length: Option<usize>,
) -> io::Result<()> {
    // NCBI reference: c++/src/objtools/align_format/showalign.cpp:2385-2388,2461-2471,340-359
    // out << ">"; if (out.tellp() > 1L) out << " "; s_WrapOutputLine(out, alnDispParams->title);
    let title = subject_title
        .filter(|s| !s.is_empty())
        .map_or_else(|| subject_id.to_string(), |t| format!("{subject_id} {t}"));
    writer.write_all(b"> ")?;
    let mut do_wrap = false;
    for (i, &c) in title.as_bytes().iter().enumerate() {
        if i > 0 && i % 60 == 0 {
            do_wrap = true;
        }
        writer.write_all(&[c])?;
        if do_wrap && (c.is_ascii_whitespace() || c == 11) {
            writer.write_all(b"\n")?;
            do_wrap = false;
        }
    }
    writeln!(writer)?;

    // Length line
    if let Some(len) = subject_length {
        writeln!(writer, "Length={}", len)?;
    }

    writeln!(writer)?;
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
// ```c++
// int CAlignFormatUtil::GetPercentMatch(int numerator, int denominator)
// {
//      if (numerator == denominator)
//         return 100;
//      else {
//        int retval =(int) (0.5 + 100.0*((double)numerator)/((double)denominator));
//        retval = min(99, retval);
//        return retval;
//      }
// }
// ```
fn ncbi_percent_match(numerator: usize, denominator: usize) -> usize {
    if numerator == denominator {
        100
    } else {
        ((0.5 + 100.0 * numerator as f64 / denominator as f64) as usize).min(99)
    }
}

/// Write HSP score information
///
/// Reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp x_DisplayAlignInfo()
///
/// Format:
/// ```text
///  Score = XXX bits (YYY),  Expect = Z.Ze-NN
///  Identities = AA/BB (CC%), Positives = DD/BB (EE%), Gaps = FF/BB (GG%)
///  Frame = +X/+Y
/// ```
pub fn write_hsp_info<W: Write>(
    writer: &mut W,
    hit: &PairwiseHit,
    config: &PairwiseConfig,
) -> io::Result<()> {
    let h = &hit.hit;

    // Score line
    // NCBI: " Score = XXX bits (YYY),  Expect = Z.Ze-NN"
    let bit_score_str = format_bitscore_ncbi(h.bit_score);
    let evalue_str = format_evalue_ncbi(h.e_value);
    if config.program == "blastp" {
        // NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:3600-3603
        // ```c
        //             if (aln_vec_info->comp_adj_method == 1)
        //             out << ", Method: Composition-based stats.";
        //             else if (aln_vec_info->comp_adj_method == 2)
        //             out << ", Method: Compositional matrix adjust.";
        // ```
        let method = match hit.comp_adjust_method {
            Some(1) => ", Method: Composition-based stats.",
            Some(2) => ", Method: Compositional matrix adjust.",
            _ => "",
        };
        writeln!(
            writer,
            " Score = {} bits ({}),  Expect = {}{}",
            bit_score_str, h.raw_score, evalue_str, method
        )?;
    } else {
        writeln!(
            writer,
            " Score = {} bits ({}),  Expect = {}",
            bit_score_str, h.raw_score, evalue_str
        )?;
    }

    // Identity/Positives/Gaps line
    let align_len = h.length;
    let num_ident = h.num_ident;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
    // ```c++
    // int CAlignFormatUtil::GetPercentMatch(int numerator, int denominator)
    // {
    //      if (numerator == denominator)
    //         return 100;
    //      else {
    //        int retval =(int) (0.5 + 100.0*((double)numerator)/((double)denominator));
    //        retval = min(99, retval);
    //        return retval;
    //      }
    // }
    // ```
    let ident_pct = ncbi_percent_match(num_ident, align_len);

    let positives = hit.positives.unwrap_or(h.num_positives.max(num_ident));
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
    // ```c++
    // int CAlignFormatUtil::GetPercentMatch(int numerator, int denominator)
    // {
    //      if (numerator == denominator)
    //         return 100;
    //      else {
    //        int retval =(int) (0.5 + 100.0*((double)numerator)/((double)denominator));
    //        retval = min(99, retval);
    //        return retval;
    //      }
    // }
    // ```
    let pos_pct = ncbi_percent_match(positives, align_len);

    let gaps = hit.gaps.unwrap_or_else(|| h.gap_letters());
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
    // ```c++
    // int CAlignFormatUtil::GetPercentMatch(int numerator, int denominator)
    // {
    //      if (numerator == denominator)
    //         return 100;
    //      else {
    //        int retval =(int) (0.5 + 100.0*((double)numerator)/((double)denominator));
    //        retval = min(99, retval);
    //        return retval;
    //      }
    // }
    // ```
    let gap_pct = ncbi_percent_match(gaps, align_len);

    // For protein/translated: show Identities, Positives, Gaps
    // For nucleotide: show Identities, Gaps (no Positives)
    if config.program.contains("blast") && !config.program.contains("blastn") {
        writeln!(
            writer,
            " Identities = {}/{} ({}%), Positives = {}/{} ({}%), Gaps = {}/{} ({}%)",
            num_ident,
            align_len,
            ident_pct,
            positives,
            align_len,
            pos_pct,
            gaps,
            align_len,
            gap_pct
        )?;
    } else if config.program == "blastn" {
        // NCBI reference: c++/src/objtools/align_format/showalign.cpp:1794-1805,2131-2141
        // ```c++
        //     x_FillIdentityInfo(aln_vec_info->alnRowInfo->sequence[0],
        //                        aln_vec_info->alnRowInfo->sequence[1],
        //                        aln_vec_info->match,
        //     ...
        //         aln_vec_info->gap = x_GetNumGaps();
        // ...
        //         if(sequence_standard[i]==sequence[i]){
        //             ...
        //             match ++;
        // ```
        // The displayed counts come from the two rows (before the lowercase masking of
        // the display), and the gaps are the gap columns of both rows.
        let (query_row, subject_row) = (
            hit.query_seq.as_deref().unwrap_or_default().as_bytes(),
            hit.subject_seq.as_deref().unwrap_or_default().as_bytes(),
        );
        let columns = query_row.len();
        let matches = query_row
            .iter()
            .zip(subject_row)
            .filter(|(q, s)| q.eq_ignore_ascii_case(s))
            .count();
        let gaps = query_row
            .iter()
            .chain(subject_row)
            .filter(|&&residue| residue == b'-')
            .count();
        // NCBI reference: c++/src/objtools/align_format/showalign.cpp:310-320
        // ```c++
        //     out<<" Identities = "<<match<<"/"<<(aln_stop+1)<<" ("<<identity<<"%"<<")";
        //     ...
        //     out<<", Gaps = "<<gap<<"/"<<(aln_stop+1)
        //        <<" ("<<CAlignFormatUtil::GetPercentMatch(gap, aln_stop+1)<<"%"<<")"<<"\n";
        //     if (!aln_is_prot){
        //         out<<" Strand="<<(master_strand==1 ? "Plus" : "Minus")
        //            <<"/"<<(slave_strand==1? "Plus" : "Minus")<<"\n";
        // ```
        // The query is always shown on its plus strand (showalign.cpp:1852-1858).
        writeln!(
            writer,
            " Identities = {}/{} ({}%), Gaps = {}/{} ({}%)",
            matches,
            columns,
            ncbi_percent_match(matches, columns),
            gaps,
            columns,
            ncbi_percent_match(gaps, columns)
        )?;
        // NCBI reference: c++/src/objtools/align_format/showalign.cpp:4011-4018
        // ```c++
        //         s_DisplayIdentityInfo(out,
        //                               ...
        //                               m_AV->StrandSign(0),
        //                               m_AV->StrandSign(1),
        // ```
        // The strand is the HSP's (its query frame; the shown query is on its plus
        // strand), not the order of the coordinates, which is the same for an HSP of one
        // letter.
        let subject_strand = if h.query_frame < 0 { "Minus" } else { "Plus" };
        writeln!(writer, " Strand=Plus/{subject_strand}")?;
    } else {
        writeln!(
            writer,
            " Identities = {}/{} ({:.0}%), Gaps = {}/{} ({:.0}%)",
            num_ident, align_len, ident_pct, gaps, align_len, gap_pct
        )?;
    }

    // Frame line (for translated searches)
    if config.show_frame {
        if let (Some(qf), Some(sf)) = (hit.query_frame, hit.subject_frame) {
            let qf_str = if qf > 0 {
                format!("+{}", qf)
            } else {
                format!("{}", qf)
            };
            let sf_str = if sf > 0 {
                format!("+{}", sf)
            } else {
                format!("{}", sf)
            };
            writeln!(writer, " Frame = {}/{}", qf_str, sf_str)?;
        }
    }

    writeln!(writer)?;
    Ok(())
}

/// Write alignment rows
///
/// Reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp x_DisplayRowData()
///
/// Format:
/// ```text
/// Query  1    MVLSPADKTN VKAAWGKVGA HAGEYGAEAL ERMFLSFPTT KTYFPHFDLS H  60
///             MVLSPADKTN VKAAWGKVGA HAGEYGAEAL ERMFLSFPTT KTYFPHFDLS H
/// Sbjct  1    MVLSPADKTN VKAAWGKVGA HAGEYGAEAL ERMFLSFPTT KTYFPHFDLS H  60
/// ```
pub fn write_alignment<W: Write>(
    writer: &mut W,
    hit: &PairwiseHit,
    config: &PairwiseConfig,
) -> io::Result<()> {
    let h = &hit.hit;

    // If we have actual sequences, display them
    if let (Some(ref qseq), Some(ref sseq)) = (&hit.query_seq, &hit.subject_seq) {
        write_alignment_with_sequences(writer, h, qseq, sseq, config)?;
    } else {
        // No sequences available - show placeholder
        writeln!(
            writer,
            "Query  {}  [... {} aa ...]  {}",
            h.q_start, h.length, h.q_end
        )?;
        writeln!(writer)?;
        writeln!(
            writer,
            "Sbjct  {}  [... {} aa ...]  {}",
            h.s_start, h.length, h.s_end
        )?;
    }

    writeln!(writer)?;
    Ok(())
}

/// Write alignment with actual sequences
fn write_alignment_with_sequences<W: Write>(
    writer: &mut W,
    hit: &Hit,
    query_seq: &str,
    subject_seq: &str,
    config: &PairwiseConfig,
) -> io::Result<()> {
    let line_len = config.line_length;
    let q_chars: Vec<char> = query_seq.chars().collect();
    let s_chars: Vec<char> = subject_seq.chars().collect();

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/showalign.cpp:96-98
    // ```c
    // static const int k_IdStartMargin = 2;
    // static const int k_SeqStopMargin = 2;
    // static const int k_StartSequenceMargin = 2;
    // ```
    //
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/showalign.cpp:1600-1625
    // ```c
    // CAlignFormatUtil::AddSpace(out, alnRoInfo->maxIdLen-alnRoInfo->seqidArray[row].size()
    //                            + k_IdStartMargin);
    // out << start;
    // CAlignFormatUtil::AddSpace(out, alnRoInfo->maxStartLen-startLen + k_StartSequenceMargin);
    // ...
    // CAlignFormatUtil::AddSpace(out, k_SeqStopMargin);
    // out << end;
    // ```
    let max_start_len = coordinate_width(hit);
    let max_id_len = "Query".len().max("Sbjct".len());
    // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1570-1573
    // ```c++
    //     if (m_Program == "blastn" || m_Program == "megablast") {
    //         display.SetMiddleLineStyle(CDisplaySeqalign::eBar);
    //         display.SetAlignType(CDisplaySeqalign::eNuc);
    // ```
    let nucleotide = config.program == "blastn";
    // A minus-strand subject is shown with decreasing coordinates (the query is always
    // on its plus strand, showalign.cpp:1852-1858).
    let subject_step: isize = if hit.s_start > hit.s_end { -1 } else { 1 };

    let mut q_next = hit.q_start as isize;
    let mut s_next = hit.s_start as isize;
    let mut offset = 0;

    while offset < q_chars.len() {
        let end = (offset + line_len).min(q_chars.len());
        let first_chunk = offset == 0;

        // NCBI reference: c++/src/objtools/align_format/showalign.cpp:2131-2155
        // ```c++
        //         if(sequence_standard[i]==sequence[i]){
        //             if(m_AlignOption & eShowMiddleLine){
        //                 if(m_MidLineStyle == eBar ) {
        //                     middle_line[i] = '|';
        //                 } else if (m_MidLineStyle == eChar){
        //                     middle_line[i] = sequence[i];
        //         ...
        //                 if (m_AlignOption & eShowMiddleLine){
        //                     middle_line[i] = ' ';
        // ```
        // The rows may carry lowercase masking, which NCBI applies only when it prints a
        // row (showalign.cpp:2495-2521), so they are compared without case, and the protein
        // middle line shows the uppercase residue.
        let middle: String = q_chars[offset..end]
            .iter()
            .zip(s_chars[offset..end].iter())
            .map(|(q, s)| {
                let q = &q.to_ascii_uppercase();
                if nucleotide {
                    if q.eq_ignore_ascii_case(s) {
                        '|'
                    } else {
                        ' '
                    }
                } else if q == s {
                    *q // Identity: show the character
                } else if is_positive_match(*q, *s, config.protein_matrix) {
                    '+' // Positive: show +
                } else {
                    ' ' // Mismatch: show space
                }
            })
            .collect();

        write_sequence_row(
            writer,
            "Query",
            &q_chars[offset..end],
            &mut q_next,
            1,
            max_start_len,
            first_chunk,
        )?;
        write_spaces(writer, max_id_len + 2 + max_start_len + 2)?;
        writeln!(writer, "{}", middle)?;
        write_sequence_row(
            writer,
            "Sbjct",
            &s_chars[offset..end],
            &mut s_next,
            subject_step,
            max_start_len,
            first_chunk,
        )?;
        writeln!(writer)?;

        offset = end;
    }

    Ok(())
}

/// Writes one row of one alignment chunk and advances `next`, the position of the row's
/// next residue, by `step` per residue.
///
/// NCBI reference: c++/src/objtools/align_format/showalign.cpp:1598-1626
/// ```c++
///     int start = alnRoInfo->seqStarts[row].front() + 1;  //+1 for 1 based
///     int end = alnRoInfo->seqStops[row].front() + 1;
///     ...
///     //not to display start and stop number for empty row
///     if ((j > 0 && end == prev_stop)
///         || (j == 0 && start == 1 && end == 1)) {
///         startLen = 0;
///     } else {
///         out << start;
///         startLen=NStr::IntToString(start).size();
///     }
///
///     CAlignFormatUtil::AddSpace(out, alnRoInfo->maxStartLen-startLen + k_StartSequenceMargin);
///     x_OutputSeq(alnRoInfo->sequence[row], m_AV->GetSeqId(row), j,
///     ...
///     CAlignFormatUtil::AddSpace(out, k_SeqStopMargin);
///
///      //not to display stop number for empty row in the middle
///     if (!(j > 0 && end == prev_stop)
///         && !(j == 0 && start == 1 && end == 1)) {
///         out << end;
///     }
///     out<<"\n";
/// ```
/// A later chunk without residues ends where the previous chunk ended
/// (`end == prev_stop`). `AddSpace` takes a `size_t`, so an omitted start gives
/// `maxStartLen + 2` spaces.
fn write_sequence_row<W: Write>(
    writer: &mut W,
    label: &str,
    residues: &[char],
    next: &mut isize,
    step: isize,
    width: usize,
    first_chunk: bool,
) -> io::Result<()> {
    let count = residues.iter().filter(|&&c| c != '-').count() as isize;
    let start = *next;
    let end = start + step * (count - 1);
    let shown = !(!first_chunk && count == 0 || first_chunk && start == 1 && end == 1);
    write!(writer, "{label}")?;
    write_spaces(writer, 2)?;
    if shown {
        write!(writer, "{start}")?;
        write_spaces(writer, width + 2 - digit_count(start.unsigned_abs()))?;
    } else {
        write_spaces(writer, width + 2)?;
    }
    for &residue in residues {
        write!(writer, "{residue}")?;
    }
    write_spaces(writer, 2)?;
    if shown {
        write!(writer, "{end}")?;
    }
    writeln!(writer)?;
    *next += step * count;
    Ok(())
}

/// Check if two amino acids are a positive match (similar)
///
/// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:2122-2149
/// ```c
/// if(sequence_standard[i]==sequence[i]){
///     match ++;
/// } else {
///     if ((m_AlignType&eProt)
///         && m_Matrix[(int)sequence_standard[i]][(int)sequence[i]] > 0){
///         positive ++;
///         if(m_AlignOption & eShowMiddleLine){
///             middle_line[i] = '+';
///         }
///     }
/// }
/// ```
fn is_positive_match(a: char, b: char, matrix: ScoringMatrix) -> bool {
    if a == b {
        return true;
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3563-3581
    // ```c++
    //     retval.Resize(k_NumAsciiChar, k_NumAsciiChar, -1000);
    //
    //     SNCBIFullScoreMatrix mtx;
    //     NCBISM_Unpack(packed_mtx, &mtx);
    //
    //     for(int i = 0; i < ePMatrixSize; ++i){
    //         for(int j = 0; j < ePMatrixSize; ++j){
    //             retval((size_t)k_PSymbol[i], (size_t)k_PSymbol[j]) =
    //                 mtx.s[(size_t)k_PSymbol[i]][(size_t)k_PSymbol[j]];
    //         }
    //     }
    //     for(int i = 0; i < ePMatrixSize; ++i) {
    //         retval((size_t)k_PSymbol[i], '*') = retval('*',(size_t)k_PSymbol[i]) = -4;
    //     }
    //     retval('*', '*') = 1;
    //     // this is to count Selenocysteine to Cysteine matches as positive
    //     retval('U', 'U') = retval('C', 'C');
    //     retval('U', 'C') = retval('C', 'C');
    //     retval('C', 'U') = retval('C', 'C');
    // ```
    protein_display_score(matrix, a as u8, b as u8) > 0
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:790-803
// ```c
// if (m_IsDbScan)
//     dbname = string("User specified sequence set (Input: ") + m_SubjectTag + string(")");
// else
//     dbname = m_DbName;
// ```
fn write_database_header<W: Write>(writer: &mut W, context: &ReportContext) -> io::Result<()> {
    if let Some(ref dbname) = context.database_name {
        writeln!(writer, "Database: {}", dbname)?;
        if let (Some(num_sequences), Some(total_letters)) = (
            context.database_num_sequences,
            context.database_total_letters,
        ) {
            writeln!(
                writer,
                "           {} sequences; {} total letters",
                num_sequences, total_letters
            )?;
        }
        writeln!(writer)?;
    } else if let Some(ref db) = context.subject_name {
        writeln!(writer, "Database: {}", db)?;
        writeln!(writer)?;
    }
    Ok(())
}

// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showdefline.cpp:75-83
// ```c
// static const char*  kHeader = "Sequences producing significant alignments:";
// ...
// static const char*  kBits = (getenv("CTOOLKIT_COMPATIBLE") ? "(bits)" : "(Bits)");
// static const size_t kBits_size = strlen(kBits);
// ...
// static const char*  kValue = "Value";
// ```
// `kBits` is a static initialized when the program starts, from whether the
// environment has CTOOLKIT_COMPATIBLE (any value, also an empty one).
fn ncbi_k_bits() -> &'static str {
    static K_BITS: std::sync::OnceLock<&'static str> = std::sync::OnceLock::new();
    K_BITS.get_or_init(|| ncbi_k_bits_for(std::env::var_os("CTOOLKIT_COMPATIBLE").is_some()))
}

fn ncbi_k_bits_for(ctoolkit_compatible: bool) -> &'static str {
    if ctoolkit_compatible {
        "(bits)"
    } else {
        "(Bits)"
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:830-837
// ```c++
//             if((m_Option & eShowSumN) || (m_Option & eShowPercentIdent)){
//                 CAlignFormatUtil::AddSpace(out, m_MaxEvalueLen - kValue_size);
//                 CAlignFormatUtil::AddSpace(out, kTwoSpaceMargin_size);
//             }
//             if(m_Option & eShowSumN){
//                 out << kN;
//             }
//             if (m_Option & eShowPercentIdent) {
// ```
fn write_subject_summary_table<W: Write>(
    writer: &mut W,
    subject_order: &[u32],
    subject_hits: &std::collections::HashMap<u32, Vec<&PairwiseHit>>,
    subject_ids: &[Arc<str>],
) -> io::Result<()> {
    write_subject_summary_table_with_sum_n(writer, subject_order, subject_hits, subject_ids, false)
}
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:830-837
// ```c++
//             if((m_Option & eShowSumN) || (m_Option & eShowPercentIdent)){
//                 CAlignFormatUtil::AddSpace(out, m_MaxEvalueLen - kValue_size);
//                 CAlignFormatUtil::AddSpace(out, kTwoSpaceMargin_size);
//             }
//             if(m_Option & eShowSumN){
//                 out << kN;
//             }
//             if (m_Option & eShowPercentIdent) {
// ```
fn write_subject_summary_table_with_sum_n<W: Write>(
    writer: &mut W,
    subject_order: &[u32],
    subject_hits: &std::collections::HashMap<u32, Vec<&PairwiseHit>>,
    subject_ids: &[Arc<str>],
    show_sum_n: bool,
) -> io::Result<()> {
    // NCBI reference: c++/src/objtools/align_format/showdefline.cpp:668-709,760-829,929-959
    // m_MaxScoreLen = kBits_size; m_MaxEvalueLen = kValue_size;
    // AddSpace(out, m_LineLen+2); out << kScore; AddSpace(out,m_MaxScoreLen-kScore_size);
    let max_score = subject_order
        .iter()
        .filter_map(|i| subject_hits.get(i)?.first())
        .map(|h| format_bitscore_ncbi(h.hit.bit_score).len())
        .max()
        .unwrap_or(6)
        .max(6);
    let max_evalue = subject_order
        .iter()
        .filter_map(|i| subject_hits.get(i)?.first())
        .map(|h| format_evalue_ncbi(h.hit.e_value).len())
        .max()
        .unwrap_or(5)
        .max(5);
    write_spaces(writer, 70)?;
    writeln!(writer, "{:<max_score$}    E", "Score")?;
    write!(
        writer,
        "{:<69}",
        "Sequences producing significant alignments:"
    )?;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:830-840
    // ```c++
    //             if((m_Option & eShowSumN) || (m_Option & eShowPercentIdent)){
    //                 CAlignFormatUtil::AddSpace(out, m_MaxEvalueLen - kValue_size);
    //                 CAlignFormatUtil::AddSpace(out, kTwoSpaceMargin_size);
    //             }
    //             if(m_Option & eShowSumN){
    //                 out << kN;
    //             }
    //             if (m_Option & eShowPercentIdent) {
    //                 out << kIdentLine2;//"ident"
    //             }
    //             out << "\n";
    // ```
    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/showdefline.cpp:811-813
    // ```c
    // out << kBits;
    // //in case m_MaxScoreLen > kBits.size()
    // CAlignFormatUtil::AddSpace(out, m_MaxScoreLen - kBits_size);
    // ```
    write!(writer, "{:<max_score$}  Value", ncbi_k_bits())?;
    let max_sum_n = subject_order
        .iter()
        .filter_map(|i| subject_hits.get(i)?.first())
        .map(|h| h.sum_n.unwrap_or(1).to_string().len())
        .max()
        .unwrap_or(1)
        .max(1);
    if show_sum_n {
        write_spaces(writer, max_evalue - 5 + 2)?;
        write!(writer, "N")?;
    }
    writeln!(writer)?;
    writeln!(writer)?;
    for s_idx in subject_order {
        let Some(shits) = subject_hits.get(s_idx) else {
            continue;
        };
        let Some(best_hit) = shits.first() else {
            continue;
        };
        let subject_id = subject_ids
            .get(*s_idx as usize)
            .map(|id| id.as_ref())
            .unwrap_or("unknown");
        let mut label = subject_id.to_string();
        if let Some(title) = best_hit.subject_title.as_deref() {
            label.push(' ');
            label.push_str(title);
        }
        // NCBI reference: c++/src/objtools/align_format/showdefline.cpp:914-931
        // actual_line_component = line_component.substr(0,m_LineLen-line_length-3); actual_line_component += kEllipsis;
        // String widths are byte counts, including non-ASCII FASTA titles.
        if label.len() > 68 {
            writer.write_all(&label.as_bytes()[..65])?;
            writer.write_all(b"...")?;
        } else {
            writer.write_all(label.as_bytes())?;
            write_spaces(writer, 68 - label.len())?;
        }
        write!(
            writer,
            "  {:<max_score$}  {:<max_evalue$}",
            format_bitscore_ncbi(best_hit.hit.bit_score),
            format_evalue_ncbi(best_hit.hit.e_value)
        )?;
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:961-965
        // ```c++
        //         if(m_Option & eShowSumN){
        //             out << kTwoSpaceMargin << (*iter)->sum_n;
        //             CAlignFormatUtil::AddSpace(out, m_MaxSumNLen -
        //                      NStr::IntToString((*iter)->sum_n).size());
        //         }
        // ```
        if show_sum_n {
            write!(writer, "  {:<max_sum_n$}", best_hit.sum_n.unwrap_or(1))?;
        }
        writeln!(writer)?;
    }
    writeln!(writer)?;
    Ok(())
}

#[inline]
fn digit_count(value: usize) -> usize {
    value.to_string().len()
}

/// Width of the start-coordinate column of an alignment block.
///
/// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:1366-1367
/// ```c++
///         size_t maxCood=max<size_t>(m_AV->GetSeqStart(row), m_AV->GetSeqStop(row));
///         maxStartLen = max<size_t>(NStr::SizetToString(maxCood).size(), maxStartLen);
/// ```
/// `GetSeqStart` and `GetSeqStop` are 0-based, so the width is that of the largest
/// 1-based coordinate minus one, and a printed 1-based start can have one digit
/// more. NCBI then pads with the `size_t` value `maxStartLen - startLen +
/// k_StartSequenceMargin`, which is one space in that case; callers write
/// `width + 2 - digits` (showalign.cpp:1614; align_format_util.cpp:932-936).
fn coordinate_width(hit: &Hit) -> usize {
    digit_count(hit.q_start.max(hit.q_end).max(hit.s_start).max(hit.s_end) - 1)
}

#[inline]
fn write_spaces<W: Write>(writer: &mut W, count: usize) -> io::Result<()> {
    for _ in 0..count {
        writer.write_all(b" ")?;
    }
    Ok(())
}

#[inline]
fn format_count_with_commas(value: i64) -> String {
    let digits = value.to_string();
    let mut formatted = String::with_capacity(digits.len() + digits.len() / 3);
    let len = digits.len();
    for (index, ch) in digits.chars().enumerate() {
        if index > 0 && (len - index) % 3 == 0 {
            formatted.push(',');
        }
        formatted.push(ch);
    }
    formatted
}

#[inline]
fn ensure_trailing_period(text: &str) -> String {
    if text.ends_with('.') {
        text.to_string()
    } else {
        format!("{text}.")
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:584-603
// ```c
// if (gapped) {
//     out << "Gapped" << "\n";
// }
// out << "Lambda      K        H";
// if (gbp) {
//     if (gapped) {
//         out << "        a         alpha    sigma";
//     } else {
//         out << "        a         alpha";
//     }
// }
// out << "\n";
// sprintf(buffer, "%#8.3g ", lambda);
// ```
// C's `%#.3g` (C99 7.21.6.1): with precision P = 3 and X the exponent of the
// `%e` conversion (after rounding to P digits), `%f` with precision P - 1 - X
// when P > X >= -4, otherwise `%e` with precision P - 1 and an exponent of a
// sign and at least two digits; `#` keeps the decimal point and the trailing
// zeros. Rust's precision formatting rounds the exact binary value half to
// even, as glibc does.
fn format_ncbi_ka_value(value: f64) -> String {
    const PRECISION: i32 = 3;
    if !value.is_finite() {
        let text = if value.is_nan() { "nan" } else { "inf" };
        return if value.is_sign_negative() {
            format!("-{text}")
        } else {
            text.to_string()
        };
    }
    let e_style = format!("{:.*e}", (PRECISION - 1) as usize, value);
    let (mantissa, exponent) = e_style
        .split_once('e')
        .expect("Rust exponent formatting has an 'e'");
    let exponent: i32 = exponent
        .parse()
        .expect("Rust exponent formatting has an integer exponent");
    if PRECISION > exponent && exponent >= -4 {
        let mut formatted = format!("{:.*}", (PRECISION - 1 - exponent) as usize, value);
        if !formatted.contains('.') {
            formatted.push('.');
        }
        formatted
    } else {
        let sign = if exponent < 0 { '-' } else { '+' };
        format!("{mantissa}e{sign}{:02}", exponent.unsigned_abs())
    }
}

#[inline]
fn write_ncbi_ka_field<W: Write>(writer: &mut W, value: f64) -> io::Result<()> {
    write!(writer, "{:>8} ", format_ncbi_ka_value(value))
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blastfmtutil.cpp:73-111
// ```c
// if (m_Program == "psiblast" || m_Program == "blastp") {
//     CBlastFormatUtil::BlastPrintReference(..., CReference::eCompBasedStats, ...);
// }
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/version.cpp:46-69
// ```c
// "Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schaffer, ..."
// "Alejandro A. Schaffer, L. Aravind, Thomas L. Madden, ..."
// ```
fn write_blastp_pairwise_intro<W: Write>(writer: &mut W, version: &str) -> io::Result<()> {
    writeln!(writer, "BLASTP {}", version)?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A."
    )?;
    writeln!(
        writer,
        "Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J."
    )?;
    writeln!(
        writer,
        "Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of"
    )?;
    writeln!(
        writer,
        "protein database search programs\", Nucleic Acids Res. 25:3389-3402."
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Reference for composition-based statistics: Alejandro A. Schaffer,"
    )?;
    writeln!(
        writer,
        "L. Aravind, Thomas L. Madden, Sergei Shavirin, John L. Spouge, Yuri"
    )?;
    writeln!(
        writer,
        "I. Wolf, Eugene V. Koonin, and Stephen F. Altschul (2001),"
    )?;
    writeln!(
        writer,
        "\"Improving the accuracy of PSI-BLAST protein database searches with"
    )?;
    writeln!(
        writer,
        "composition-based statistics and other refinements\", Nucleic Acids"
    )?;
    writeln!(writer, "Res. 29:2994-3005.")?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)?;
    Ok(())
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blastfmtutil.cpp:125-138
// ```c
// string dbString = (html) ? "<b>Database:</b> " : "Database: ";
// str << dbString << definition_line << endl;
// ...
// out << "           " << NStr::IntToString(nNumSeqs,NStr::fWithCommas)
//     << " sequences; " << NStr::UInt8ToString(nTotalLength,NStr::fWithCommas)
//     << " total letters" << endl;
// ```
// NCBI reference: c++/src/corelib/ncbistr.cpp:5088-5340 (WrapIt, fWrap_FlatFile)
// enum EScore { eForced, ePunct, eComma, eSpace, eNewline };
// if (score >= best_score && score_pos > pos0) { best_pos = score_pos; best_score = score; }
// NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:279-291
// NStr::Wrap(str, line_len, string_l, NStr::fWrap_FlatFile);
fn write_flatfile_wrapped(writer: &mut impl Write, text: &str, width: usize) -> io::Result<()> {
    let bytes = text.as_bytes();
    let mut pos = 0;
    // NCBI reference: c++/src/corelib/ncbistr.cpp:5100,5153-5157,5137
    // SIZE_TYPE pos=0,len=str.size(),nl_pos=0; if(nl_pos<=pos) nl_pos=str.find('\n',pos);
    // bool thisPartHasBackspace = false;
    let mut newline = 0;
    while pos < bytes.len() {
        let mut has_backspace = false;
        if newline <= pos {
            newline = bytes[pos..]
                .iter()
                .position(|&c| c == b'\n')
                .map_or(bytes.len(), |n| pos + n);
        }
        let pos0 = if newline - pos <= width { newline } else { pos };
        let mut scan = pos0;
        let mut column = 0;
        let mut best = pos;
        let mut best_score = 0;
        while scan < bytes.len() && column <= width {
            let c = bytes[scan];
            let mut score = 0;
            let mut score_pos = scan;
            if c == b'\n' {
                best = scan;
                best_score = 4;
                break;
            }
            if c.is_ascii_whitespace() || c == 11 {
                score = 3;
            } else if c == b',' && column < width && scan + 1 < bytes.len() {
                score = 2;
                score_pos += 1;
            } else if c == b'-' && column < width && scan + 1 < bytes.len() {
                score = 1;
                score_pos += 1;
            }
            if score >= best_score && score_pos > pos0 {
                best = score_pos;
                best_score = score;
            }
            while scan + 1 < bytes.len() && bytes[scan + 1] == 8 {
                scan += 1;
                column = column.saturating_sub(1);
                has_backspace = true;
            }
            scan += 1;
            column += 1;
        }
        // NCBI reference: c++/src/corelib/ncbistr.cpp:5260-5266,5277-5300
        // if (best_pos != len) { best_pos=len; thisPartHasBackspace=true; }
        // if (thisPartHasBackspace) { /* eat backspaces and preceding characters */ }
        if best_score != 4 && column <= width && best != bytes.len() {
            best = bytes.len();
            has_backspace = true;
        }
        let mut line = Vec::with_capacity(best - pos);
        for &c in &bytes[pos..best] {
            if has_backspace && c == 8 {
                line.pop();
            } else {
                line.push(c);
            }
        }
        writer.write_all(&line)?;
        writer.write_all(b"\n")?;
        pos = best;
        if best_score == 3 {
            while pos < bytes.len() && bytes[pos] == b' ' {
                pos += 1;
            }
            if pos < bytes.len() && bytes[pos] == b'\n' {
                pos += 1;
            }
        }
        if best_score == 4 {
            pos += 1;
        }
        while pos < bytes.len() && bytes[pos] == 8 {
            pos += 1;
        }
    }
    Ok(())
}

/// The outfmt 0 prolog of blastp (the version, the references and the database), which
/// NCBI writes before it reads the first query batch (blast_format.cpp, `PrintProlog`).
pub fn write_blastp_pairwise_prolog<W: Write>(
    writer: &mut W,
    version: &str,
    database_name: &str,
    database_num_sequences: usize,
    database_total_letters: usize,
) -> io::Result<()> {
    write_blastp_pairwise_intro(writer, version)?;
    // The blank lines before the first query are the query's (as TBLASTX's prolog).
    write_blastp_database_header_spacing(
        writer,
        database_name,
        database_num_sequences,
        database_total_letters,
        1,
    )
}

fn write_blastp_database_header<W: Write>(
    writer: &mut W,
    database_name: &str,
    database_num_sequences: usize,
    database_total_letters: usize,
) -> io::Result<()> {
    write_blastp_database_header_spacing(
        writer,
        database_name,
        database_num_sequences,
        database_total_letters,
        3,
    )
}
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:533-538
// ```c++
//         CAlignFormatUtil::AddSpace(out, 11);
//         out << NStr::Int8ToString(tot_num_seqs, NStr::fWithCommas) <<
//             " sequences; " <<
//             NStr::Int8ToString(tot_length, NStr::fWithCommas) <<
//             " total letters\n\n";
//         return;
// ```
fn write_blastp_database_header_spacing<W: Write>(
    writer: &mut W,
    database_name: &str,
    database_num_sequences: usize,
    database_total_letters: usize,
    trailing: usize,
) -> io::Result<()> {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:513-525
    // ```c++
    //         out << "Database: ";
    //
    //         string db_titles = dbinfo->definition;
    //         Int8 tot_num_seqs = static_cast<Int8>(dbinfo->number_seqs);
    //         Int8 tot_length = dbinfo->total_length;
    //
    //         for (size_t i = 1; i < dbinfo_list.size(); i++) {
    //             db_titles += "; " + dbinfo_list[i].definition;
    //             tot_num_seqs += static_cast<Int8>(dbinfo_list[i].number_seqs);
    //             tot_length += dbinfo_list[i].total_length;
    //         }
    //
    //         x_WrapOutputLine(db_titles, line_length, out);
    // ```
    write!(writer, "Database: ")?;
    write_flatfile_wrapped(writer, &ensure_trailing_period(database_name), 68)?;
    writeln!(
        writer,
        "           {} sequences; {} total letters",
        format_count_with_commas(i64::try_from(database_num_sequences).unwrap_or(0)),
        format_count_with_commas(i64::try_from(database_total_letters).unwrap_or(0))
    )?;
    for _ in 0..trailing {
        writeln!(writer)?;
    }
    Ok(())
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:712-746
// ```c
// out << label << "= ";
// ...
// out << "\nLength=";
// out << cbs.GetInst().GetLength() <<"\n";
// ```
fn write_blastp_query_header<W: Write>(
    writer: &mut W,
    query_name: &str,
    query_length: usize,
) -> io::Result<()> {
    // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:729-746
    // out << label << "= "; x_WrapOutputLine(all_id_str, line_len, out, html);
    // out << "\nLength=" << cbs.GetInst().GetLength() << "\n";
    write!(writer, "Query= ")?;
    write_flatfile_wrapped(writer, query_name, 68)?;
    writeln!(writer)?;
    writeln!(writer, "Length={}", query_length)?;
    Ok(())
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1510-1514
// ```c
// m_Outfile << "\n\n"
//           << "***** " << CBlastFormatUtil::kNoHitsFound << " *****" << "\n"
//           << "\n\n";
// ```
fn write_no_hits_found<W: Write>(writer: &mut W) -> io::Result<()> {
    // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1510-1514
    // m_Outfile << "\n\n" << "***** " << kNoHitsFound << " *****\n\n\n";
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer, "***** No hits found *****")?;
    writeln!(writer)?;
    writeln!(writer)?;
    Ok(())
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:445-477
// ```c
// m_Outfile << NcbiEndl;
// if (kbp_ungap) {
//     CBlastFormatUtil::PrintKAParameters(..., false, gbp);
// }
// ...
// if (kbp_gap) {
//     CBlastFormatUtil::PrintKAParameters(..., true, gbp);
// }
// m_Outfile << "\n";
// m_Outfile << "Effective search space used: "
//           << summary.GetSearchSpace() << "\n";
// ```
// NCBI c++/src/algo/blast/api/local_blast.cpp:177-180,204-208:
// if (status != 0) {
//     pair<double, double> tmp_pair(-1.0, -1.0);
//     CRef<CBlastAncillaryData> tmp_ancillary_data(
//         new CBlastAncillaryData(tmp_pair, tmp_pair, tmp_pair, 0));
// }
// NCBI c++/src/algo/blast/format/blast_format.cpp:451-478 and
// c++/src/objtools/align_format/align_format_util.cpp:581-613:
// PrintKAParameters emits both -1 Karlin blocks with no Gumbel columns.
fn write_tblastn_unsearched_query_footer<W: Write>(writer: &mut W) -> io::Result<()> {
    write_tblastn_unsearched_query_footer_spacing(writer, true)
}
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:477-478
// ```c++
//     m_Outfile << "Effective search space used: " <<
//                         summary.GetSearchSpace() << "\n";
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1491-1491
// ```c++
//     m_Outfile << "\n\n";
// ```
fn write_tblastn_unsearched_query_footer_spacing<W: Write>(
    writer: &mut W,
    trailing: bool,
) -> io::Result<()> {
    writeln!(writer)?;
    writeln!(writer, "Lambda      K        H")?;
    for _ in 0..3 {
        write_ncbi_ka_field(writer, -1.0)?;
    }
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer, "Gapped")?;
    writeln!(writer, "Lambda      K        H")?;
    for _ in 0..3 {
        write_ncbi_ka_field(writer, -1.0)?;
    }
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer, "Effective search space used: 0")?;
    if trailing {
        writeln!(writer)?;
        writeln!(writer)?;
    }
    Ok(())
}

// NCBI c++/src/algo/blast/format/blast_format.cpp:445-478:
// if (kbp_ungap) CBlastFormatUtil::PrintKAParameters(..., false, gbp);
// if (kbp_gap) CBlastFormatUtil::PrintKAParameters(..., true, gbp);
// m_Outfile << "Effective search space used: " << summary.GetSearchSpace() << "\n";
fn write_blastp_query_footer<W: Write>(
    writer: &mut W,
    ungapped_karlin: KarlinParams,
    gapped_karlin: KarlinParams,
    gumbel: BlastGumbelBlk,
    effective_search_space: i64,
) -> io::Result<()> {
    write_blastp_query_footer_spacing(
        writer,
        ungapped_karlin,
        gapped_karlin,
        gumbel,
        effective_search_space,
        true,
    )
}
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:477-478
// ```c++
//     m_Outfile << "Effective search space used: " <<
//                         summary.GetSearchSpace() << "\n";
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1491-1491
// ```c++
//     m_Outfile << "\n\n";
// ```
fn write_blastp_query_footer_spacing<W: Write>(
    writer: &mut W,
    ungapped_karlin: KarlinParams,
    gapped_karlin: KarlinParams,
    gumbel: BlastGumbelBlk,
    effective_search_space: i64,
    trailing: bool,
) -> io::Result<()> {
    writeln!(writer)?;
    writeln!(writer, "Lambda      K        H        a         alpha")?;
    write_ncbi_ka_field(writer, ungapped_karlin.lambda)?;
    write_ncbi_ka_field(writer, ungapped_karlin.k)?;
    write_ncbi_ka_field(writer, ungapped_karlin.h)?;
    write_ncbi_ka_field(writer, gumbel.a_un)?;
    write_ncbi_ka_field(writer, gumbel.alpha_un)?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer, "Gapped")?;
    writeln!(
        writer,
        "Lambda      K        H        a         alpha    sigma"
    )?;
    write_ncbi_ka_field(writer, gapped_karlin.lambda)?;
    write_ncbi_ka_field(writer, gapped_karlin.k)?;
    write_ncbi_ka_field(writer, gapped_karlin.h)?;
    write_ncbi_ka_field(writer, gumbel.a)?;
    write_ncbi_ka_field(writer, gumbel.alpha)?;
    write_ncbi_ka_field(writer, gumbel.sigma)?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Effective search space used: {}",
        effective_search_space
    )?;
    if trailing {
        writeln!(writer)?;
        writeln!(writer)?;
    }
    Ok(())
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:503-579
// ```c
// out << "  Database: ";
// ...
// out << "  Number of letters in database: ";
// out << NStr::Int8ToString(dbinfo->total_length, NStr::fWithCommas) << "\n";
// out << "  Number of sequences in database:  ";
// out << NStr::IntToString(dbinfo->number_seqs, NStr::fWithCommas) << "\n";
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:2258-2285
// ```c
// m_Outfile << "\n\nMatrix: " << options.GetMatrixName() << "\n";
// ...
// m_Outfile << "Gap Penalties: Existence: "
//         << options.GetGapOpeningCost() << ", Extension: "
//         << gap_extension << "\n";
// ...
// m_Outfile << "Neighboring words threshold: "
//         << options.GetWordThreshold() << "\n";
// m_Outfile << "Window for multiple hits: "
//         << options.GetWindowSize() << "\n";
// ```
fn write_blastp_final_footer<W: Write>(
    writer: &mut W,
    report: &BlastpPairwiseReport,
) -> io::Result<()> {
    write_final_database_report(
        writer,
        &report.database_name,
        report.database_num_sequences,
        report.database_total_letters,
    )?;
    writeln!(writer, "Matrix: {}", report.matrix_name)?;
    writeln!(
        writer,
        "Gap Penalties: Existence: {}, Extension: {}",
        report.gap_open, report.gap_extend
    )?;
    if report.word_threshold != 0.0 {
        // `GetWordThreshold()` is a double, written with the stream's default format.
        writeln!(
            writer,
            "Neighboring words threshold: {}",
            cpp_default_double(report.word_threshold)
        )?;
    }
    if report.window_size != 0 {
        writeln!(writer, "Window for multiple hits: {}", report.window_size)?;
    }
    Ok(())
}

/// The database block of the epilog and the blank lines before the matrix line
/// (align_format_util.cpp:503-579 and blast_format.cpp:2258-2262, quoted at
/// `write_blastp_final_footer`).
fn write_final_database_report<W: Write>(
    writer: &mut W,
    database_name: &str,
    database_num_sequences: usize,
    database_total_letters: usize,
) -> io::Result<()> {
    // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:543-544
    // out << "  Database: "; x_WrapOutputLine(dbinfo->definition, line_length, out);
    write!(writer, "  Database: ")?;
    write_flatfile_wrapped(writer, &ensure_trailing_period(database_name), 68)?;
    writeln!(writer, "    Posted date:  Unknown")?;
    writeln!(
        writer,
        "  Number of letters in database: {}",
        format_count_with_commas(i64::try_from(database_total_letters).unwrap_or(0))
    )?;
    writeln!(
        writer,
        "  Number of sequences in database:  {}",
        format_count_with_commas(i64::try_from(database_num_sequences).unwrap_or(0))
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:372-424
// ```c
// CBlastFormatUtil::BlastPrintVersionInfo(m_Program, m_IsHTML, m_Outfile);
// ...
// CBlastFormatUtil::BlastPrintReference(...);
// ...
// CBlastFormatUtil::BlastPrintReference(..., CReference::eCompBasedStats, ...);
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1520-1589
// ```c
// if ( (!m_IsBl2Seq || m_IsDbScan) && !(m_DisableKAStats || kIsGlobal) ) {
//     x_DisplayDeflines(aln_set, itr_num, prev_seqids);
// }
// ...
// display.DisplaySeqalign(m_Outfile);
// x_PrintOneQueryFooter(*results.GetAncillaryData());
// ```
pub fn write_blastp_pairwise_report<W: Write>(
    hits: &[PairwiseHit],
    writer: &mut W,
    config: &PairwiseConfig,
    queries: &[BlastpPairwiseQuery],
    subject_ids: &[Arc<str>],
    report: &BlastpPairwiseReport,
    mut probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut super::query_warnings::QueryWarnings<'_>>,
    epilog: bool,
) -> io::Result<()> {
    let mut buffered = io::BufWriter::new(writer);
    let writer = &mut buffered;

    write_blastp_pairwise_intro(writer, &report.version)?;
    // The prolog ends with one blank line; each query's preamble starts with two
    // (blast_format.cpp:1491), after the query's warnings.
    write_blastp_database_header_spacing(
        writer,
        &report.database_name,
        report.database_num_sequences,
        report.database_total_letters,
        1,
    )?;

    // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1411
    // ```c++
    // CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
    // ```
    // Each hit keeps its index in `hits` (the final result order), which identifies it
    // to the probe.
    let mut hits_by_query: Vec<Vec<(HspIndex, &PairwiseHit)>> = vec![Vec::new(); queries.len()];
    for (hsp_index, hit) in hits.iter().enumerate() {
        if let Some(bucket) = hits_by_query.get_mut(hit.hit.q_idx as usize) {
            bucket.push((hsp_index, hit));
        }
    }

    for (q_idx, query) in queries.iter().enumerate() {
        // The query's warnings come before its preamble (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, writer)?;
        }
        // NCBI blast_format.cpp:1491: m_Outfile << "\n\n";
        writeln!(writer)?;
        writeln!(writer)?;
        write_blastp_query_header(writer, &query.query_name, query.query_length)?;
        let query_hits = &hits_by_query[q_idx];
        if query_hits.is_empty() {
            write_no_hits_found(writer)?;
            // NCBI c++/src/algo/blast/api/local_blast.cpp:177-180,204-208:
            // if (m_PrelimSearch->CheckInternalData() != 0)
            //     new CBlastAncillaryData(tmp_pair, tmp_pair, tmp_pair, 0);
            // The all-invalid batch uses -1 sentinel blocks. An invalid query
            // within a searched batch has null blocks (blast_results.cpp:82-103),
            // so only the blank lines and the zero search space are written
            // (blast_format.cpp:445-477), as TBLASTN does.
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
                    query.effective_search_space,
                    false,
                )?;
            }
            continue;
        }

        use std::collections::HashMap;
        let mut subject_hits: HashMap<u32, Vec<&PairwiseHit>> = HashMap::new();
        let mut subject_hit_indices: HashMap<u32, Vec<HspIndex>> = HashMap::new();
        let mut subject_order: Vec<u32> = Vec::new();
        for &(hsp_index, hit) in query_hits {
            let s_idx = hit.hit.s_idx;
            if !subject_hits.contains_key(&s_idx) {
                subject_order.push(s_idx);
            }
            subject_hits.entry(s_idx).or_default().push(hit);
            subject_hit_indices
                .entry(s_idx)
                .or_default()
                .push(hsp_index);
        }

        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:606-611
        // ```c++
        //     CShowBlastDefline showdef(*aln_set, *m_Scope,
        //                               defline_length == -1 ? kFormatLineLength:defline_length,
        //                               m_NumSummary + additional);
        // ```
        // The table follows `x_InitDeflineTable` (`write_blastn_description_table`): the
        // subject's highest bit score and that HSP's E-value, and the protein title.
        let described = &subject_order[..subject_order.len().min(report.num_descriptions)];
        write_blastn_description_table(
            writer,
            described,
            subject_order.len() <= report.num_descriptions,
            &subject_hits,
            subject_ids,
            false,
            true,
        )?;
        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1520-1589
        // x_DisplayDeflines(aln_set, ...); ... display.DisplaySeqalign(m_Outfile);
        writeln!(writer)?;

        // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:1014-1040
        // (PruneSeqalign keeps the alignments of the first `m_NumAlignments` subjects).
        let aligned = &subject_order[..subject_order.len().min(report.num_alignments)];
        for &s_idx in aligned {
            let subject_id = subject_ids
                .get(s_idx as usize)
                .map(|id| id.as_ref())
                .unwrap_or("unknown");
            let shits = subject_hits
                .get(&s_idx)
                .expect("subject order must reference existing grouped hits");
            let first_hit = shits
                .first()
                .expect("subject group must contain at least one HSP");
            // NCBI reference: c++/src/objtools/align_format/showalign.cpp:3613-3632
            // ```c++
            //     if(show_defline) {
            // 		...
            // 				string deflines = x_PrintDefLine(bsp_handle, aln_vec_info);
            // 				out<< deflines;
            // 		...
            // 			out << "\n";
            // ```
            // The heading is written before the first HSP of the subject; the probe marks
            // it with that HSP's index without changing the written bytes.
            let first_index = subject_hit_indices[&s_idx][0];
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.subject_begin(first_index);
            }
            // NCBI reference: c++/src/objtools/align_format/showalign.cpp:2273
            // ```c++
            //             alnDispParams->title = CDeflineGenerator().GenerateDefline(bsp_handle);
            // ```
            // The subject (a protein sequence) is shown with its title (`defline.rs`).
            let heading = super::defline::ncbi_protein_title(
                &subject_defline(subject_id, first_hit.subject_title.as_deref()),
                false,
            );
            write_subject_header(writer, &heading, None, first_hit.subject_length)?;
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.subject_end(first_index);
            }
            // NCBI reference: c++/src/objtools/align_format/showalign.cpp:1970-1973
            // ```c++
            // subid=&(avRef->GetSeqId(1));
            // bool showDefLine = previousId.Empty() || !subid->Match(*previousId);
            // x_DisplayAlnvecInfo(out, alnvecInfo,showDefLine);
            // ```
            // The subject defline above is the part of the first x_DisplayAlnvecInfo call
            // that precedes the HSP; the probe marks the rest of each call (the score
            // block and the alignment) without changing the written bytes.
            for (hit, &hsp_index) in shits.iter().zip(&subject_hit_indices[&s_idx]) {
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.begin(hsp_index);
                }
                write_hsp_info(writer, hit, config)?;
                write_alignment(writer, hit, config)?;
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.end(hsp_index);
                }
            }
        }

        write_blastp_query_footer_spacing(
            writer,
            query.ungapped_karlin,
            report.gapped_karlin,
            report.gumbel,
            query.effective_search_space,
            false,
        )?;
    }

    // NCBI reference: c++/src/app/blast/blastp_app.cpp:294-295
    // ```c++
    //         BLAST_PROF_START( APP.POST );
    //         formatter.PrintEpilog(opt);
    // ```
    // A search that stops at a query batch (`Empty CBlastQueryVector` with `-query_loc`)
    // writes no epilog.
    if !epilog {
        return writer.flush();
    }

    // The blank lines before the epilog (blast_format.cpp:2249).

    writeln!(writer)?;

    writeln!(writer)?;

    write_blastp_final_footer(writer, report)?;
    writer.flush()
}

// =============================================================================
// BLASTN pairwise report
// =============================================================================

/// One BLASTN query of the pairwise report (NCBI `CBlastAncillaryData` and the query's
/// validity).
///
/// NCBI reference: c++/src/algo/blast/api/blast_results.cpp:82-104
/// ```c++
///     // find the first valid context corresponding to this query
///     ...
///     m_SearchSpace = ctx->eff_searchsp;
///     ...
///     s_InitializeKarlinBlk(sbp->kbp_std[ctx_index], &m_UngappedKarlinBlk);
/// ```
#[derive(Debug, Clone)]
pub struct BlastnPairwiseQuery {
    /// The FASTA defline without `>`.
    pub query_name: String,
    pub query_length: usize,
    /// The ungapped and the gapped block of the query's first valid context; `None` for
    /// an invalid query.
    pub karlin: Option<(KarlinParams, KarlinParams)>,
    pub effective_search_space: i64,
}

/// Run-level data of the BLASTN pairwise report.
#[derive(Debug, Clone)]
pub struct BlastnPairwiseReport {
    pub version: String,
    /// `-task megablast` (NCBI's `m_Megablast`; the program name stays `blastn`).
    pub megablast: bool,
    pub database_name: String,
    pub database_num_sequences: usize,
    pub database_total_letters: usize,
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    /// The task is megablast or blastn, for which a gap extension cost of 0 is printed as
    /// the PMID 10890397 value (NCBI compares the task name, `m_Program`).
    pub zero_gap_extension_formula: bool,
    /// The two-hit window (`options.GetWindowSize()`): 40 for dc-megablast.
    pub window_size: usize,
    /// Subjects in the description table and with alignments, per query.
    pub num_descriptions: usize,
    pub num_alignments: usize,
    /// Every query of the batch is invalid, so NCBI did not search it (all queries of such
    /// a batch get the `-1` footer, local_blast.cpp:177-207). Indexed like the queries.
    pub unsearched: Vec<bool>,
    /// Whether NCBI's epilog (the database and statistics footer, `PrintEpilog`) ends the
    /// report.
    pub epilog: bool,
    /// Whether the report starts with NCBI's prolog (`PrintProlog`); false when the caller
    /// wrote it before the search, as NCBI does.
    pub prolog: bool,
}

// The description table of the BLASTN report.
//
// NCBI reference: c++/src/objtools/align_format/showdefline.cpp:1004-1011
// ```c++
// void CShowBlastDefline::DisplayBlastDefline(CNcbiOstream & out)
// {
//     x_InitDeflineTable();
//     ...
//     x_DisplayDefline(out);
// ```
// NCBI reference: c++/src/objtools/align_format/showdefline.cpp:1100-1142
// ```c++
//         if(!is_first_aln && !(subid->Match(*previous_id))) {
//             SScoreInfo* sci = x_GetScoreInfoForTable(hit, num_align);
//             if(sci){
//                 m_ScoreList.push_back(sci);
//                 if(m_MaxScoreLen < sci->bit_string.size()){
//                     m_MaxScoreLen = sci->bit_string.size();
//                 }
//                 if(m_MaxTotalScoreLen < sci->total_bit_string.size()){
//                     m_MaxTotalScoreLen = sci->total_bit_string.size();
//                 }
//     ...
//         if (num_align < m_NumToShow) { //no adding if number to show already reached
//             hit.Set().push_back(*iter);
//         }
//     ...
//     //the last hit
//     SScoreInfo* sci = x_GetScoreInfoForTable(hit, num_align);
//     if(sci){
//          m_ScoreList.push_back(sci);
//         if(m_MaxScoreLen < sci->bit_string.size()){
//             m_MaxScoreLen = sci->bit_string.size();
//         }
//         if(m_MaxTotalScoreLen < sci->total_bit_string.size()){
//             m_MaxScoreLen = sci->total_bit_string.size();
//         }
// ```
// NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:4285-4310
// ```c++
//         total_bits += bits;
//     ...
//         if (bits > highest_bits) {
//             highest_bits = bits;
//             lowest_evalue = evalue;
//         }
//     ...
//     seqSetInfo->total_bit_score = total_bits;
//     seqSetInfo->bit_score = highest_bits;
//     seqSetInfo->evalue = lowest_evalue;
// ```
// A row shows the highest bit score of the subject's HSPs and that HSP's E-value. The
// widths start from "(Bits)" and "Value" (showdefline.cpp:1059-1062). The last row's
// total score sets the score width when it is wider than every earlier total (the
// width that NCBI assigns there is the score width). The last row exists only when the
// subjects fit in the table: otherwise the loop stops with an empty last group. A total
// is wider than "Total" only above 99999 bits (`%5.3le`), which `format_bitscore_ncbi`
// formats alike. The header and the rows are those of `write_subject_summary_table`.
//
// NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:4229-4230,4259
// ```c++
//     seqSetInfo->sum_n = sum_n == -1 ? 1:sum_n ;
//     ...
//     seqSetInfo = GetSeqAlignCalcParams(*(aln.Get().front()));
// ```
// NCBI reference: c++/src/objtools/align_format/showdefline.cpp:1061,1118-1120
// ```c++
//     m_MaxSumNLen =1;
//     ...
//                 if( m_MaxSumNLen < NStr::IntToString(sci->sum_n).size()){
//                     m_MaxSumNLen = NStr::IntToString(sci->sum_n).size();
//                 }
// ```
// With `show_sum_n` (NCBI's `eShowSumN`: ungapped sum statistics, blast_format.cpp:165-167,
// 503-504) a last column shows the linked-set size of the subject's first HSP; the HSP
// without a "sum_n" score (`num` of 1) shows 1.
fn write_blastn_description_table<W: Write>(
    writer: &mut W,
    described: &[u32],
    last_row_counted: bool,
    subject_hits: &std::collections::HashMap<u32, Vec<&PairwiseHit>>,
    subject_ids: &[Arc<str>],
    show_sum_n: bool,
    protein: bool,
) -> io::Result<()> {
    let rows: Vec<(u32, &PairwiseHit, f64, i32)> = described
        .iter()
        .filter_map(|s_idx| {
            let hits = subject_hits.get(s_idx)?;
            let first = *hits.first()?;
            let mut best = first;
            for hit in hits {
                if hit.hit.bit_score > best.hit.bit_score {
                    best = hit;
                }
            }
            let sum_n = first.sum_n.filter(|&num| num > 1).unwrap_or(1);
            Some((
                *s_idx,
                best,
                hits.iter().map(|hit| hit.hit.bit_score).sum(),
                sum_n,
            ))
        })
        .collect();
    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/showdefline.cpp:1059
    // ```c
    // m_MaxScoreLen = kBits_size;
    // ```
    let mut max_score = ncbi_k_bits().len();
    let mut max_evalue = "Value".len();
    let mut max_total = "Total".len();
    let mut max_sum_n = 1;
    for (index, (_, best, total, sum_n)) in rows.iter().enumerate() {
        max_score = max_score.max(format_bitscore_ncbi(best.hit.bit_score).len());
        max_evalue = max_evalue.max(format_evalue_ncbi(best.hit.e_value).len());
        max_sum_n = max_sum_n.max(sum_n.to_string().len());
        let total_len = format_bitscore_ncbi(*total).len();
        if index + 1 == rows.len() && last_row_counted {
            if max_total < total_len {
                max_score = total_len;
            }
        } else {
            max_total = max_total.max(total_len);
        }
    }
    write_spaces(writer, 70)?;
    writeln!(writer, "{:<max_score$}    E", "Score")?;
    write!(
        writer,
        "{:<69}",
        "Sequences producing significant alignments:"
    )?;
    write!(writer, "{:<max_score$}  Value", ncbi_k_bits())?;
    // NCBI reference: c++/src/objtools/align_format/showdefline.cpp:829-836
    // ```c++
    //             out << kValue;
    //             if((m_Option & eShowSumN) || (m_Option & eShowPercentIdent)){
    //                 CAlignFormatUtil::AddSpace(out, m_MaxEvalueLen - kValue_size);
    //                 CAlignFormatUtil::AddSpace(out, kTwoSpaceMargin_size);
    //             }
    //             if(m_Option & eShowSumN){
    //                 out << kN;
    //             }
    // ```
    if show_sum_n {
        write_spaces(writer, max_evalue - "Value".len() + 2)?;
        write!(writer, "N")?;
    }
    writeln!(writer)?;
    writeln!(writer)?;
    for (s_idx, best, _, sum_n) in &rows {
        let subject_id = subject_ids
            .get(*s_idx as usize)
            .map(|id| id.as_ref())
            .unwrap_or("unknown");
        // NCBI reference: c++/src/objtools/align_format/showdefline.cpp:498
        // The description keeps the prefixes (`fLeavePrefixSuffix`, `defline.rs`).
        let defline = subject_defline(subject_id, best.subject_title.as_deref());
        let label = if protein {
            super::defline::ncbi_protein_title(&defline, true)
        } else {
            super::defline::ncbi_nucleotide_title(&defline, true)
        };
        // NCBI reference: c++/src/objtools/align_format/showdefline.cpp:915-918,930
        // ```c++
        //         if(line_component.size()+line_length > m_LineLen){
        //             actual_line_component = line_component.substr(0, m_LineLen -
        //                                                           line_length - 3);
        //             actual_line_component += kEllipsis;
        //     ...
        //         CAlignFormatUtil::AddSpace(out, m_LineLen - line_length);
        // ```
        // String widths are byte counts, including non-ASCII FASTA titles.
        if label.len() > 68 {
            writer.write_all(&label.as_bytes()[..65])?;
            writer.write_all(b"...")?;
        } else {
            writer.write_all(label.as_bytes())?;
            write_spaces(writer, 68 - label.len())?;
        }
        write!(
            writer,
            "  {:<max_score$}  {:<max_evalue$}",
            format_bitscore_ncbi(best.hit.bit_score),
            format_evalue_ncbi(best.hit.e_value)
        )?;
        // NCBI reference: c++/src/objtools/align_format/showdefline.cpp:961-965
        // ```c++
        //         if(m_Option & eShowSumN){
        //             out << kTwoSpaceMargin << (*iter)->sum_n;
        //             CAlignFormatUtil::AddSpace(out, m_MaxSumNLen -
        //                      NStr::IntToString((*iter)->sum_n).size());
        //         }
        // ```
        if show_sum_n {
            write!(writer, "  {sum_n:<max_sum_n$}")?;
        }
        writeln!(writer)?;
    }
    writeln!(writer)
}

// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:387-407
// ```c++
//       CBlastFormatUtil::BlastPrintVersionInfo(m_Program, m_IsHTML,
//                                               m_Outfile);
//     ...
//     m_Outfile << NcbiEndl << NcbiEndl;
//     ...
//     if (m_Megablast)
//         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
//                                           m_Outfile, CReference::eMegaBlast);
// ```
// NCBI reference: c++/src/algo/blast/api/version.cpp:58-61
// ```c++
//     // eMegaBlast
//     "Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), \
// \"A greedy algorithm for aligning DNA sequences\", \
// J Comput Biol 2000; 7(1-2):203-14.",
// ```
// The reference is printed after "Reference: " and wrapped at 68 columns
// (blastfmtutil.cpp:116-121); the layout is that of `write_translated_pairwise_intro`.
fn write_megablast_pairwise_intro(writer: &mut impl Write, version: &str) -> io::Result<()> {
    writeln!(writer, "BLASTN {version}")?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Reference: Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb"
    )?;
    writeln!(
        writer,
        "Miller (2000), \"A greedy algorithm for aligning DNA sequences\", J"
    )?;
    writeln!(writer, "Comput Biol 2000; 7(1-2):203-14.")?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)
}

// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:451-478
// ```c++
//     m_Outfile << NcbiEndl;
//     if (kbp_ungap) {
//         CBlastFormatUtil::PrintKAParameters(kbp_ungap->Lambda, kbp_ungap->K,
//                                             kbp_ungap->H, kFormatLineLength,
//                                             m_Outfile, false, gbp);
//     }
//
//     m_Outfile << "\n";
//     if (kbp_gap) {
//         ...
//     }
//
//     m_Outfile << "Effective search space used: " <<
//         summary.GetSearchSpace() << "\n";
// ```
// NCBI reference: c++/src/algo/blast/core/blast_setup.c:472-480
// ```c
//     if (program_number == eBlastTypeBlastn ||
//     ...
//         /* disable new FSC rules for nucleotide case for now */
//         if (sbp && sbp->gbp) {
//             sfree(sbp->gbp);
// ```
// Without a Gumbel block, PrintKAParameters (align_format_util.cpp:585-619) prints only
// the Lambda, K and H columns. An invalid query has neither Karlin block (`karlin` is
// `None`) and a search space of 0.
fn write_nucleotide_query_footer<W: Write>(
    writer: &mut W,
    karlin: Option<(KarlinParams, KarlinParams)>,
    effective_search_space: i64,
) -> io::Result<()> {
    write_nucleotide_query_footer_spacing(writer, karlin, effective_search_space, true)
}

/// `write_nucleotide_query_footer`, with the two blank lines that follow it (the next
/// query's preamble or the epilog) when `trailing`.
fn write_nucleotide_query_footer_spacing<W: Write>(
    writer: &mut W,
    karlin: Option<(KarlinParams, KarlinParams)>,
    effective_search_space: i64,
    trailing: bool,
) -> io::Result<()> {
    writeln!(writer)?;
    if let Some((ungapped, _)) = karlin {
        writeln!(writer, "Lambda      K        H")?;
        write_ncbi_ka_field(writer, ungapped.lambda)?;
        write_ncbi_ka_field(writer, ungapped.k)?;
        write_ncbi_ka_field(writer, ungapped.h)?;
        writeln!(writer)?;
    }
    writeln!(writer)?;
    if let Some((_, gapped)) = karlin {
        writeln!(writer, "Gapped")?;
        writeln!(writer, "Lambda      K        H")?;
        write_ncbi_ka_field(writer, gapped.lambda)?;
        write_ncbi_ka_field(writer, gapped.k)?;
        write_ncbi_ka_field(writer, gapped.h)?;
        writeln!(writer)?;
    }
    writeln!(writer)?;
    writeln!(
        writer,
        "Effective search space used: {}",
        effective_search_space
    )?;
    // The two blank lines of the next query's preamble or of the epilog
    // (blast_format.cpp:1491, 2249).
    if trailing {
        writeln!(writer)?;
        writeln!(writer)?;
    }
    Ok(())
}

// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:2261-2279
// ```c++
//     if (m_Program == "blastn" || m_Program == "megablast") {
//         m_Outfile << "\n\nMatrix: " << "blastn matrix " <<
//                         options.GetMatchReward() << " " <<
//                         options.GetMismatchPenalty() << "\n";
//     }
//     ...
//     double gap_extension = (double) options.GetGapExtensionCost();
//     if ((m_Program == "megablast" || m_Program == "blastn") && options.GetGapExtensionCost() == 0)
//     {
//         // Formula from PMID 10890397 applies if both gap values are zero.
//         gap_extension = -2*options.GetMismatchPenalty() + options.GetMatchReward();
//         gap_extension /= 2.0;
//     }
//     m_Outfile << "Gap Penalties: Existence: "
//             << options.GetGapOpeningCost() << ", Extension: "
//             << gap_extension << "\n";
// ```
// The word threshold and the window size are 0 for blastn (blast_format.cpp:2281-2288),
// so those lines are not printed. The extension is a whole or half number, which the
// C++ stream and Rust print alike (`2`, `2.5`).
fn write_blastn_final_footer<W: Write>(
    writer: &mut W,
    report: &BlastnPairwiseReport,
) -> io::Result<()> {
    write_final_database_report(
        writer,
        &report.database_name,
        report.database_num_sequences,
        report.database_total_letters,
    )?;
    writeln!(
        writer,
        "Matrix: blastn matrix {} {}",
        report.reward, report.penalty
    )?;
    // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:2270-2288
    // ```c
    //     if (options.GetGappedMode() == true) {
    //         double gap_extension = (double) options.GetGapExtensionCost();
    //         if ((m_Program == "megablast" || m_Program == "blastn") && options.GetGapExtensionCost() == 0)
    //         { // Formula from PMID 10890397 applies if both gap values are zero.
    //                gap_extension = -2*options.GetMismatchPenalty() + options.GetMatchReward();
    //                gap_extension /= 2.0;
    //         }
    //         m_Outfile << "Gap Penalties: Existence: "
    //                 << options.GetGapOpeningCost() << ", Extension: "
    //                 << gap_extension << "\n";
    //     }
    //     if (options.GetWordThreshold()) {
    //         m_Outfile << "Neighboring words threshold: " <<
    //                         options.GetWordThreshold() << "\n";
    //     }
    //     if (options.GetWindowSize()) {
    //         m_Outfile << "Window for multiple hits: " <<
    //                         options.GetWindowSize() << "\n";
    //     }
    // ```
    // BLASTN's word threshold is 0 for every task (BLAST_WORD_THRESHOLD_BLASTN and
    // BLAST_WORD_THRESHOLD_MEGABLAST, blast_options.h).
    let gap_extension = if report.zero_gap_extension_formula && report.gap_extend == 0 {
        f64::from(-2 * report.penalty + report.reward) / 2.0
    } else {
        f64::from(report.gap_extend)
    };
    writeln!(
        writer,
        "Gap Penalties: Existence: {}, Extension: {}",
        report.gap_open, gap_extension
    )?;
    if report.window_size != 0 {
        writeln!(writer, "Window for multiple hits: {}", report.window_size)?;
    }
    Ok(())
}

/// The start of the BLASTN report, which NCBI writes before it searches: the program
/// and its reference, then the subjects.
///
/// NCBI reference: c++/src/app/blast/blastn_app.cpp:258-261
/// ```c
///         formatter.PrintProlog();
///
///         /*** Process the input ***/
///         CBatchSizeMixer mixer(SplitQuery_GetChunkSize(opt.GetProgram())-1000);
/// ```
pub fn write_blastn_pairwise_prolog<W: Write>(
    writer: &mut W,
    version: &str,
    megablast: bool,
    database_name: &str,
    database_num_sequences: usize,
    database_total_letters: usize,
) -> io::Result<()> {
    if megablast {
        write_megablast_pairwise_intro(writer, version)?;
    } else {
        write_translated_pairwise_intro(writer, "BLASTN", version)?;
    }
    write_blastp_database_header_spacing(
        writer,
        database_name,
        database_num_sequences,
        database_total_letters,
        1,
    )
}

/// Writes the BLASTN pairwise report (outfmt 0).
///
/// `hits` is the final hit list of the run in its order (queries in input order, then
/// subjects and HSPs in result order); each hit's index in it identifies it to the probe.
///
/// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1491-1589
/// ```c++
///     m_Outfile << "\n\n";
///     CBlastFormatUtil::AcknowledgeBlastQuery(*bioseq, kFormatLineLength,
///     ...
///         m_Outfile << "\n\n"
///                   << "***** " << CBlastFormatUtil::kNoHitsFound << " *****" << "\n"
///                   << "\n\n";
///         x_PrintOneQueryFooter(*results.GetAncillaryData());
///     ...
///         x_DisplayDeflines(aln_set, itr_num, prev_seqids);
///     ...
///     m_Outfile << "\n";
///     ...
///     CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_NumAlignments);
///     ...
///     display.DisplaySeqalign(m_Outfile);
///     ...
///     x_PrintOneQueryFooter(*results.GetAncillaryData());
/// ```
// Kept out of line: this runs only for outfmt 0 or hit records, and inlining it into the
// shared post-processing makes every BLASTN run compile it in Wasm hosts.
#[inline(never)]
pub fn write_blastn_pairwise_report<W: Write>(
    hits: &[PairwiseHit],
    writer: &mut W,
    config: &PairwiseConfig,
    queries: &[BlastnPairwiseQuery],
    subject_ids: &[Arc<str>],
    report: &BlastnPairwiseReport,
    mut probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut super::query_warnings::QueryWarnings<'_>>,
) -> io::Result<()> {
    let mut buffered = io::BufWriter::new(writer);
    let writer = &mut buffered;

    if report.prolog {
        write_blastn_pairwise_prolog(
            writer,
            &report.version,
            report.megablast,
            &report.database_name,
            report.database_num_sequences,
            report.database_total_letters,
        )?;
    }
    let mut hits_by_query: Vec<Vec<(HspIndex, &PairwiseHit)>> = vec![Vec::new(); queries.len()];
    for (hsp_index, hit) in hits.iter().enumerate() {
        if let Some(bucket) = hits_by_query.get_mut(hit.hit.q_idx as usize) {
            bucket.push((hsp_index, hit));
        }
    }

    for (q_idx, query) in queries.iter().enumerate() {
        // The query's warnings (and those of its batch's reading) come before its preamble
        // (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, writer)?;
        }
        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1489-1491
        // ```c++
        //     // print the preamble for this query
        //
        //     m_Outfile << "\n\n";
        // ```
        writer.write_all(b"\n\n")?;
        write_blastp_query_header(writer, &query.query_name, query.query_length)?;
        let query_hits = &hits_by_query[q_idx];
        if query_hits.is_empty() {
            write_no_hits_found(writer)?;
            if report.unsearched.get(q_idx).copied().unwrap_or(false) {
                write_tblastn_unsearched_query_footer_spacing(writer, false)?;
            } else {
                write_nucleotide_query_footer_spacing(
                    writer,
                    query.karlin,
                    query.effective_search_space,
                    false,
                )?;
            }
            continue;
        }

        use std::collections::HashMap;
        let mut subject_hits: HashMap<u32, Vec<&PairwiseHit>> = HashMap::new();
        let mut subject_hit_indices: HashMap<u32, Vec<HspIndex>> = HashMap::new();
        let mut subject_order: Vec<u32> = Vec::new();
        for &(hsp_index, hit) in query_hits {
            let s_idx = hit.hit.s_idx;
            if !subject_hits.contains_key(&s_idx) {
                subject_order.push(s_idx);
            }
            subject_hits.entry(s_idx).or_default().push(hit);
            subject_hit_indices
                .entry(s_idx)
                .or_default()
                .push(hsp_index);
        }

        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:606-611
        // ```c++
        //     CShowBlastDefline showdef(*aln_set, *m_Scope,
        //                               defline_length == -1 ? kFormatLineLength:defline_length,
        //                               m_NumSummary + additional);
        // ```
        let described = &subject_order[..subject_order.len().min(report.num_descriptions)];
        write_blastn_description_table(
            writer,
            described,
            subject_order.len() <= report.num_descriptions,
            &subject_hits,
            subject_ids,
            false,
            false,
        )?;
        writeln!(writer)?;

        // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:1014-1040
        // (PruneSeqalign keeps the alignments of the first `m_NumAlignments` subjects).
        let aligned = &subject_order[..subject_order.len().min(report.num_alignments)];
        for &s_idx in aligned {
            let subject_id = subject_ids
                .get(s_idx as usize)
                .map(|id| id.as_ref())
                .unwrap_or("unknown");
            let shits = &subject_hits[&s_idx];
            let first_hit = shits
                .first()
                .expect("subject group must contain at least one HSP");
            // NCBI reference: c++/src/objtools/align_format/showalign.cpp:3613-3632
            // ```c++
            //     if(show_defline) {
            // 		...
            // 				string deflines = x_PrintDefLine(bsp_handle, aln_vec_info);
            // 				out<< deflines;
            // 		...
            // 			out << "\n";
            // ```
            // The heading is written before the first HSP of the subject; the probe marks
            // it with that HSP's index without changing the written bytes.
            let first_index = subject_hit_indices[&s_idx][0];
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.subject_begin(first_index);
            }
            // NCBI reference: c++/src/objtools/align_format/showalign.cpp:2273
            // ```c++
            //             alnDispParams->title = CDeflineGenerator().GenerateDefline(bsp_handle);
            // ```
            // The subject's defline is its title (`defline.rs`).
            let defline = subject_defline(subject_id, first_hit.subject_title.as_deref());
            // The search rejects the titles that have none (`blastn` `search`).
            let heading = super::defline::ncbi_nucleotide_title(&defline, false);
            write_subject_header(writer, &heading, None, first_hit.subject_length)?;
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.subject_end(first_index);
            }
            // NCBI reference: c++/src/objtools/align_format/showalign.cpp:1970-1973
            // ```c++
            // subid=&(avRef->GetSeqId(1));
            // bool showDefLine = previousId.Empty() || !subid->Match(*previousId);
            // x_DisplayAlnvecInfo(out, alnvecInfo,showDefLine);
            // ```
            for (hit, &hsp_index) in shits.iter().zip(&subject_hit_indices[&s_idx]) {
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.begin(hsp_index);
                }
                write_hsp_info(writer, hit, config)?;
                write_alignment(writer, hit, config)?;
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.end(hsp_index);
                }
            }
        }

        write_nucleotide_query_footer_spacing(
            writer,
            query.karlin,
            query.effective_search_space,
            false,
        )?;
    }

    // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:317-318
    // ```c
    //         BLAST_PROF_START( APP.POST );
    //         formatter.PrintEpilog(opt);
    // ```
    // An error in a later query batch skips the epilog, and nothing follows the last
    // query's footer. The epilog starts with two blank lines (blast_format.cpp:2249).
    if report.epilog {
        writer.write_all(b"\n\n")?;
        write_blastn_final_footer(writer, report)?;
    }
    writer.flush()
}

// =============================================================================
// TBLASTX pairwise report (outfmt 0)
// =============================================================================

/// Per-query data of the TBLASTX pairwise report.
///
/// NCBI reference: c++/src/algo/blast/api/blast_results.cpp:82-100
/// ```c++
///     // find the first valid context corresponding to this query
///     for (i = 0; i < context_per_query; i++) {
///         BlastContextInfo *ctx = query_info->contexts +
///                                 query_number * context_per_query + i;
///         if (ctx->is_valid) {
///             m_SearchSpace = ctx->eff_searchsp;
///     ...
///     if (sbp->kbp_std) {
///         s_InitializeKarlinBlk(sbp->kbp_std[ctx_index], &m_UngappedKarlinBlk);
///     }
/// ```
#[derive(Debug, Clone)]
pub struct TblastxPairwiseQuery {
    /// The FASTA defline without `>`.
    pub query_name: String,
    pub query_length: usize,
    /// The ungapped block of the query's first valid context; `None` for an invalid query.
    pub karlin: Option<KarlinParams>,
    pub effective_search_space: i64,
}

/// Run-level data of the TBLASTX pairwise report.
#[derive(Debug, Clone)]
pub struct TblastxPairwiseReport {
    pub version: String,
    pub database_name: String,
    pub database_num_sequences: usize,
    pub database_total_letters: usize,
    /// `-threshold` and `-window_size` (the epilog prints them when they are not 0).
    pub word_threshold: f64,
    pub window_size: usize,
    /// Subjects in the description table and with alignments, per query.
    pub num_descriptions: usize,
    pub num_alignments: usize,
    /// The queries of the batches without a valid context, which NCBI does not search
    /// (local_blast.cpp:177-208). Indexed like the queries.
    pub unsearched: Vec<bool>,
    /// Whether NCBI's epilog (`PrintEpilog`) ends the report.
    pub epilog: bool,
    /// Whether the report starts with NCBI's prolog (`PrintProlog`); false when the caller
    /// wrote it before the search, as NCBI does.
    pub prolog: bool,
}

/// The start of the TBLASTX report, which NCBI writes before it searches: the program and
/// its reference, then the subjects.
///
/// NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:173-178
/// ```c
///         formatter.PrintProlog();
///
///         /*** Process the input ***/
///         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
///
///             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
/// ```
pub fn write_tblastx_pairwise_prolog<W: Write>(
    writer: &mut W,
    version: &str,
    database_name: &str,
    database_num_sequences: usize,
    database_total_letters: usize,
) -> io::Result<()> {
    write_translated_pairwise_intro(writer, "TBLASTX", version)?;
    write_blastp_database_header_spacing(
        writer,
        database_name,
        database_num_sequences,
        database_total_letters,
        1,
    )
}

// NCBI reference: c++/src/objtools/align_format/showalign.cpp:3593-3606
// ```c++
//             out<<" Score = "<<bit_score_buf<<" ";
//             out<<"bits ("<<aln_vec_info->score<<"),"<<"  ";
//             out<<"Expect";
//             if (aln_vec_info->sum_n > 0) {
//             out << "(" << aln_vec_info->sum_n << ")";
//             }
//             out << " = " << evalue_buf;
//     ...
//     out << "\n";
// ```
// NCBI reference: c++/src/objtools/align_format/showalign.cpp:310-332
// ```c++
//     out<<" Identities = "<<match<<"/"<<(aln_stop+1)<<" ("<<identity<<"%"<<")";
//     if(aln_is_prot){
//         out<<", Positives = "<<(positive + match)<<"/"<<(aln_stop+1)
// 			<<" ("<<CAlignFormatUtil::GetPercentMatch(positive + match, aln_stop+1)<<"%"<<")";
//     }
//     out<<", Gaps = "<<gap<<"/"<<(aln_stop+1)
//        <<" ("<<CAlignFormatUtil::GetPercentMatch(gap, aln_stop+1)<<"%"<<")"<<"\n";
//     ...
//     if(master_frame != 0 && slave_frame != 0) {
//         out <<" Frame = " << ((master_frame > 0) ? "+" : "")
//             << master_frame <<"/"<<((slave_frame > 0) ? "+" : "")
//             << slave_frame<<"\n";
//     ...
//     out<<"\n";
// ```
// NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:4229-4230
// ```c++
//     seqSetInfo->sum_n = sum_n == -1 ? 1:sum_n ;
// ```
// The Seq-align has a "sum_n" score only for an HSP of a linked set of two or more
// (blast_seqalign.cpp:1180-1183); without it the display's `sum_n` is -1, which prints no
// count. The counts are those of the displayed rows (`x_FillIdentityInfo`), which the
// caller puts in the record; a TBLASTX HSP has no gaps.
fn write_tblastx_hsp_info<W: Write>(writer: &mut W, hit: &PairwiseHit) -> io::Result<()> {
    let h = &hit.hit;
    let bits = format_bitscore_ncbi(h.bit_score);
    let evalue = format_evalue_ncbi(h.e_value);
    write!(writer, " Score = {bits} bits ({}),  Expect", h.raw_score)?;
    if let Some(sum_n) = hit.sum_n.filter(|&num| num > 1) {
        write!(writer, "({sum_n})")?;
    }
    writeln!(writer, " = {evalue}")?;
    let length = h.length;
    let positives = hit.positives.unwrap_or(h.num_positives);
    let gaps = hit.gaps.unwrap_or(0);
    writeln!(
        writer,
        " Identities = {}/{length} ({}%), Positives = {positives}/{length} ({}%), Gaps = {gaps}/{length} ({}%)",
        h.num_ident,
        ncbi_percent_match(h.num_ident, length),
        ncbi_percent_match(positives, length),
        ncbi_percent_match(gaps, length),
    )?;
    let frame = |frame: i8| {
        if frame > 0 {
            format!("+{frame}")
        } else {
            frame.to_string()
        }
    };
    writeln!(
        writer,
        " Frame = {}/{}",
        frame(hit.query_frame.unwrap_or(1)),
        frame(hit.subject_frame.unwrap_or(1))
    )?;
    writeln!(writer)
}

// NCBI reference: c++/src/objtools/align_format/showalign.cpp:1598-1626
// ```c++
//     int start = alnRoInfo->seqStarts[row].front() + 1;  //+1 for 1 based
//     int end = alnRoInfo->seqStops[row].front() + 1;
//     ...
//         out << start;
//         startLen=NStr::IntToString(start).size();
//     ...
//     CAlignFormatUtil::AddSpace(out, alnRoInfo->maxStartLen-startLen + k_StartSequenceMargin);
//     x_OutputSeq(alnRoInfo->sequence[row], m_AV->GetSeqId(row), j,
//     ...
//     CAlignFormatUtil::AddSpace(out, k_SeqStopMargin);
//     ...
//         out << end;
//     ...
//     out<<"\n";
// ```
// NCBI reference: c++/src/objtools/alnmgr/alnvec.cpp:230,305,358
// ```c++
//     scrn_width *= width;
//     ...
//                         scrn_lft_seq_pos = plus ? start : stop;
//     ...
//                     scrn_rgt_seq_pos = plus ? stop : start;
// ```
// NCBI reference: c++/src/objtools/align_format/showalign.cpp:1776-1784
// ```c++
//     CAlignFormatUtil:: AddSpace(out, alnRoInfo->maxIdLen + k_IdStartMargin + alnRoInfo->maxStartLen + k_StartSequenceMargin);
//     x_OutputSeq(alnRoInfo->middleLine, no_id, j, (int)actualLineLen, 0, row, false, alnRoInfo->masked_regions[row], out);
// ```
// NCBI reference: c++/src/objtools/align_format/showalign.cpp:2131-2150
// ```c++
//         if(sequence_standard[i]==sequence[i]){
//             if(m_AlignOption & eShowMiddleLine){
//                 if(m_MidLineStyle == eBar ) {
//                     middle_line[i] = '|';
//                 } else if (m_MidLineStyle == eChar){
//                     middle_line[i] = sequence[i];
//                 }
//             }
//             match ++;
//         } else {
//             if ((m_AlignType&eProt)
//                 && m_Matrix[(int)sequence_standard[i]][(int)sequence[i]] > 0){
//                 positive ++;
//                 if(m_AlignOption & eShowMiddleLine){
//                     if (m_MidLineStyle == eChar)
//                         middle_line[i] = '+';
//                 }
// ```
// Both rows are translated (widths 3): a 60-residue chunk spans 180 nucleotides, each row
// runs in the direction of its frame, and a minus-frame row is not reversed
// (showalign.cpp:1852-1858 applies only without widths). The coordinate width is that of
// the largest 0-based nucleotide coordinate (`coordinate_width`). The query row may carry
// lowercase masks, which NCBI applies only when it prints the row (showalign.cpp:2495-2521),
// so the middle line compares the uppercase residues.
fn write_tblastx_alignment<W: Write>(
    writer: &mut W,
    hit: &PairwiseHit,
    config: &PairwiseConfig,
) -> io::Result<()> {
    let (Some(qseq), Some(sseq)) = (&hit.query_seq, &hit.subject_seq) else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing TBLASTX alignment strings",
        ));
    };
    let q = qseq.as_bytes();
    let s = sseq.as_bytes();
    if q.len() != s.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "TBLASTX alignment string lengths differ",
        ));
    }
    let width = coordinate_width(&hit.hit);
    let direction = |frame: Option<i8>| -> isize {
        if frame.unwrap_or(1) < 0 {
            -1
        } else {
            1
        }
    };
    let q_direction = direction(hit.query_frame);
    let s_direction = direction(hit.subject_frame);
    let mut q_pos = hit.hit.q_start as isize;
    let mut s_pos = hit.hit.s_start as isize;
    for (qpart, spart) in q
        .chunks(config.line_length)
        .zip(s.chunks(config.line_length))
    {
        let q_end = q_pos + q_direction * (3 * qpart.len() as isize - 1);
        let s_end = s_pos + s_direction * (3 * spart.len() as isize - 1);
        write!(writer, "Query  {q_pos}")?;
        write_spaces(writer, width + 2 - digit_count(q_pos.unsigned_abs()))?;
        writer.write_all(qpart)?;
        writeln!(writer, "  {q_end}")?;
        write_spaces(writer, 7 + width + 2)?;
        for (&qc, &sc) in qpart.iter().zip(spart) {
            let query_residue = qc.to_ascii_uppercase();
            let mid = if query_residue == sc {
                query_residue
            } else if is_positive_match(query_residue as char, sc as char, config.protein_matrix) {
                b'+'
            } else {
                b' '
            };
            writer.write_all(&[mid])?;
        }
        writeln!(writer)?;
        write!(writer, "Sbjct  {s_pos}")?;
        write_spaces(writer, width + 2 - digit_count(s_pos.unsigned_abs()))?;
        writer.write_all(spart)?;
        writeln!(writer, "  {s_end}")?;
        writeln!(writer)?;
        q_pos = q_end + q_direction;
        s_pos = s_end + s_direction;
    }
    Ok(())
}

// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:451-478
// ```c++
//     m_Outfile << NcbiEndl;
//     if (kbp_ungap) {
//         CBlastFormatUtil::PrintKAParameters(kbp_ungap->Lambda,
//                                             kbp_ungap->K, kbp_ungap->H,
//                                             kFormatLineLength, m_Outfile,
//                                             false, gbp);
//     }
//     ...
//     m_Outfile << "\n";
//     if (kbp_gap) {
//     ...
//     m_Outfile << "\n";
//     m_Outfile << "Effective search space used: " <<
//                         summary.GetSearchSpace() << "\n";
// ```
// NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:583-603
// ```c++
//     char buffer[256];
//     if (gapped) {
//         out << "Gapped" << "\n";
//     }
//     out << "Lambda      K        H";
//     if (gbp) {
//     ...
//     out << "\n";
//     sprintf(buffer, "%#8.3g ", lambda);
//     out << buffer;
// ```
// TBLASTX is ungapped and has no Gumbel block: only the ungapped Karlin block of the
// query's first valid context is printed. A query without a valid context in a searched
// batch prints neither block and a search space of 0 (blast_results.cpp:92-94).
fn write_tblastx_query_footer<W: Write>(
    writer: &mut W,
    query: &TblastxPairwiseQuery,
) -> io::Result<()> {
    writeln!(writer)?;
    if let Some(karlin) = query.karlin {
        writeln!(writer, "Lambda      K        H")?;
        write_ncbi_ka_field(writer, karlin.lambda)?;
        write_ncbi_ka_field(writer, karlin.k)?;
        write_ncbi_ka_field(writer, karlin.h)?;
        writeln!(writer)?;
    }
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Effective search space used: {}",
        if query.karlin.is_some() {
            query.effective_search_space
        } else {
            0
        }
    )
}

// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:2249,2256-2259
// ```c++
//     m_Outfile << NcbiEndl << NcbiEndl;
//     ...
//     if ( !m_IsBl2Seq || m_IsDbScan) {
//         CBlastFormatUtil::PrintDbReport(m_DbInfo, kFormatLineLength,
//                                         m_Outfile, false);
//     }
// ```
// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:2266-2288
// ```c++
//         m_Outfile << "\n\nMatrix: " << options.GetMatrixName() << "\n";
//     }
//
//     if (options.GetGappedMode() == true) {
//     ...
//     }
//     if (options.GetWordThreshold()) {
//         m_Outfile << "Neighboring words threshold: " <<
//                         options.GetWordThreshold() << "\n";
//     }
//     if (options.GetWindowSize()) {
//         m_Outfile << "Window for multiple hits: " <<
//                         options.GetWindowSize() << "\n";
//     }
// ```
// TBLASTX is ungapped, so the gap penalties are not printed.
fn write_tblastx_epilog<W: Write>(
    writer: &mut W,
    report: &TblastxPairwiseReport,
) -> io::Result<()> {
    writer.write_all(b"\n\n")?;
    write_final_database_report(
        writer,
        &report.database_name,
        report.database_num_sequences,
        report.database_total_letters,
    )?;
    writeln!(writer, "Matrix: BLOSUM62")?;
    if report.word_threshold != 0.0 {
        // `GetWordThreshold()` is a double, written with the stream's default format.
        writeln!(
            writer,
            "Neighboring words threshold: {}",
            cpp_default_double(report.word_threshold)
        )?;
    }
    if report.window_size != 0 {
        writeln!(writer, "Window for multiple hits: {}", report.window_size)?;
    }
    Ok(())
}

/// A `double` written to a C++ stream with the default flags: `%g` with precision 6
/// (`1e+06` for 1000000, `1.23457e+06` for 1234567).
fn cpp_default_double(value: f64) -> String {
    fn trim(text: &str) -> &str {
        if text.contains('.') {
            text.trim_end_matches('0').trim_end_matches('.')
        } else {
            text
        }
    }
    if value == 0.0 {
        return "0".to_string();
    }
    // glibc's `%g` of an infinity or a NaN.
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.to_string();
    }
    if value.is_nan() {
        return if value.is_sign_negative() {
            "-nan"
        } else {
            "nan"
        }
        .to_string();
    }
    let scientific = format!("{value:.5e}");
    let (mantissa, exponent) = scientific.split_once('e').expect("exponent");
    let exponent: i32 = exponent.parse().expect("exponent digits");
    if !(-4..6).contains(&exponent) {
        let sign = if exponent < 0 { '-' } else { '+' };
        format!("{}e{sign}{:02}", trim(mantissa), exponent.abs())
    } else {
        trim(&format!("{value:.*}", (5 - exponent) as usize)).to_string()
    }
}

/// Writes the TBLASTX pairwise report (outfmt 0).
///
/// `hits` is the final hit list of the run in its order (queries in input order, then
/// subjects and HSPs in result order); each hit's index in it identifies it to the probe.
///
/// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1491-1589
/// ```c++
///     m_Outfile << "\n\n";
///     CBlastFormatUtil::AcknowledgeBlastQuery(*bioseq, kFormatLineLength,
///     ...
///         m_Outfile << "\n\n"
///                   << "***** " << CBlastFormatUtil::kNoHitsFound << " *****" << "\n"
///                   << "\n\n";
///         x_PrintOneQueryFooter(*results.GetAncillaryData());
///     ...
///         x_DisplayDeflines(aln_set, itr_num, prev_seqids);
///     ...
///     m_Outfile << "\n";
///     ...
///     CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_NumAlignments);
///     ...
///     display.DisplaySeqalign(m_Outfile);
///     ...
///     x_PrintOneQueryFooter(*results.GetAncillaryData());
/// ```
// Kept out of line: this runs only for outfmt 0.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub fn write_tblastx_pairwise_report<W: Write>(
    hits: &[PairwiseHit],
    writer: &mut W,
    config: &PairwiseConfig,
    queries: &[TblastxPairwiseQuery],
    subject_ids: &[Arc<str>],
    report: &TblastxPairwiseReport,
    mut probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut super::query_warnings::QueryWarnings<'_>>,
) -> io::Result<()> {
    let mut buffered = io::BufWriter::new(writer);
    let writer = &mut buffered;

    if report.prolog {
        write_tblastx_pairwise_prolog(
            writer,
            &report.version,
            &report.database_name,
            report.database_num_sequences,
            report.database_total_letters,
        )?;
    }
    let mut hits_by_query: Vec<Vec<(HspIndex, &PairwiseHit)>> = vec![Vec::new(); queries.len()];
    for (hsp_index, hit) in hits.iter().enumerate() {
        if let Some(bucket) = hits_by_query.get_mut(hit.hit.q_idx as usize) {
            bucket.push((hsp_index, hit));
        }
    }

    for (q_idx, query) in queries.iter().enumerate() {
        // The query's warnings come before its preamble (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, writer)?;
        }
        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1489-1491
        // ```c++
        //     // print the preamble for this query
        //
        //     m_Outfile << "\n\n";
        // ```
        writer.write_all(b"\n\n")?;
        write_blastp_query_header(writer, &query.query_name, query.query_length)?;
        let query_hits = &hits_by_query[q_idx];
        if query_hits.is_empty() {
            write_no_hits_found(writer)?;
            if report.unsearched.get(q_idx).copied().unwrap_or(false) {
                write_tblastn_unsearched_query_footer_spacing(writer, false)?;
            } else {
                write_tblastx_query_footer(writer, query)?;
            }
            continue;
        }

        use std::collections::HashMap;
        let mut subject_hits: HashMap<u32, Vec<&PairwiseHit>> = HashMap::new();
        let mut subject_hit_indices: HashMap<u32, Vec<HspIndex>> = HashMap::new();
        let mut subject_order: Vec<u32> = Vec::new();
        for &(hsp_index, hit) in query_hits {
            let s_idx = hit.hit.s_idx;
            if !subject_hits.contains_key(&s_idx) {
                subject_order.push(s_idx);
            }
            subject_hits.entry(s_idx).or_default().push(hit);
            subject_hit_indices
                .entry(s_idx)
                .or_default()
                .push(hsp_index);
        }

        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:606-613
        // ```c++
        //     CShowBlastDefline showdef(*aln_set, *m_Scope,
        //                               defline_length == -1 ? kFormatLineLength:defline_length,
        //                               m_NumSummary + additional);
        //     ...
        //     m_Outfile << "\n";
        // ```
        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:165-167
        // ```c++
        //     if (use_sum_statistics && m_IsUngappedSearch) {
        //         m_ShowLinkedSetSize = true;
        //     }
        // ```
        let described = &subject_order[..subject_order.len().min(report.num_descriptions)];
        write_blastn_description_table(
            writer,
            described,
            subject_order.len() <= report.num_descriptions,
            &subject_hits,
            subject_ids,
            true,
            false,
        )?;
        writeln!(writer)?;

        // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:1014-1040
        // (PruneSeqalign keeps the alignments of the first `m_NumAlignments` subjects).
        let aligned = &subject_order[..subject_order.len().min(report.num_alignments)];
        for &s_idx in aligned {
            let subject_id = subject_ids
                .get(s_idx as usize)
                .map(|id| id.as_ref())
                .unwrap_or("unknown");
            let shits = &subject_hits[&s_idx];
            let first_hit = shits
                .first()
                .expect("subject group must contain at least one HSP");
            // NCBI reference: c++/src/objtools/align_format/showalign.cpp:3613-3632
            // ```c++
            //     if(show_defline) {
            // 		...
            // 				string deflines = x_PrintDefLine(bsp_handle, aln_vec_info);
            // 				out<< deflines;
            // 		...
            // 			out << "\n";
            // ```
            // The heading is written before the first HSP of the subject; the probe marks
            // it with that HSP's index without changing the written bytes.
            let first_index = subject_hit_indices[&s_idx][0];
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.subject_begin(first_index);
            }
            // NCBI reference: c++/src/objtools/align_format/showalign.cpp:2273
            // ```c++
            //             alnDispParams->title = CDeflineGenerator().GenerateDefline(bsp_handle);
            // ```
            let defline = subject_defline(subject_id, first_hit.subject_title.as_deref());
            let heading = super::defline::ncbi_nucleotide_title(&defline, false);
            write_subject_header(writer, &heading, None, first_hit.subject_length)?;
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.subject_end(first_index);
            }
            // NCBI reference: c++/src/objtools/align_format/showalign.cpp:3971-3980
            // ```c++
            // 		x_ShowAlnvecInfo(out,aln_vec_info,show_defline);
            // 	...
            //     out<<"\n";
            // ```
            for (hit, &hsp_index) in shits.iter().zip(&subject_hit_indices[&s_idx]) {
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.begin(hsp_index);
                }
                write_tblastx_hsp_info(writer, hit)?;
                write_tblastx_alignment(writer, hit, config)?;
                writeln!(writer)?;
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.end(hsp_index);
                }
            }
        }

        write_tblastx_query_footer(writer, query)?;
    }

    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:209
    // ```c
    //         formatter.PrintEpilog(opt);
    // ```
    // An error in a later query batch skips the epilog, and nothing follows the last
    // query's footer.
    if report.epilog {
        write_tblastx_epilog(writer, report)?;
    }
    writer.flush()
}

// =============================================================================
// Main pairwise output function
// =============================================================================

/// Write hits in pairwise format (outfmt 0)
///
/// Reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp DisplaySeqalign()
pub fn write_pairwise<W: Write>(
    hits: &[PairwiseHit],
    writer: &mut W,
    config: &PairwiseConfig,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
    context: &ReportContext,
) -> io::Result<()> {
    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:1408-1469
    // ```c
    // string CDisplaySeqalign::x_DisplayRowData(SAlnRowInfo *alnRoInfo)
    // {
    //     ...
    //     CNcbiOstrstream out;
    //     ...
    // }
    //
    // void CDisplaySeqalign::x_DisplayRowData(SAlnRowInfo *alnRoInfo, CNcbiOstream& out)
    // {
    //     ...
    //     out << rowdata;
    // }
    // ```
    let mut buffered = io::BufWriter::new(writer);
    let writer = &mut buffered;

    // Write program header
    let version = context.version.as_deref().unwrap_or("0.1.0");
    writeln!(writer, "{} {}", context.program.to_uppercase(), version)?;
    writeln!(writer)?;

    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:790-803
    // ```c
    // dbname = string("User specified sequence set (Input: ") + m_SubjectTag + string(")");
    // ...
    // tabinfo.PrintHeader(..., dbname, ...);
    // ```
    write_database_header(writer, context)?;

    // Write query info
    if let Some(ref query) = context.query_name {
        writeln!(writer, "Query= {}", query)?;
        writeln!(writer)?;
        if let Some(query_length) = context.query_length {
            writeln!(writer, "Length={}", query_length)?;
        }
        writeln!(writer)?;
    }

    if hits.is_empty() {
        writeln!(writer, " ***** No hits found *****")?;
        writer.flush()?;
        return Ok(());
    }

    // Group hits by subject
    use std::collections::HashMap;
    let mut subject_hits: HashMap<u32, Vec<&PairwiseHit>> = HashMap::new();
    let mut subject_order: Vec<u32> = Vec::new();

    for hit in hits {
        // NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:2278-2315
        // ```c
        // string CDisplaySeqalign::x_PrintDefLine(...)
        // {
        //     ...
        //     out << ">";
        //     ...
        // }
        // ```
        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
        // ```c
        // typedef struct BlastHSPList {
        //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
        //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
        //                       Set to 0 if not applicable */
        // } BlastHSPList;
        // ```
        let s_idx = hit.hit.s_idx;
        if !subject_hits.contains_key(&s_idx) {
            subject_order.push(s_idx);
        }
        subject_hits.entry(s_idx).or_default().push(hit);
    }

    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/showdefline.cpp:75-83
    // ```c
    // static const char*  kHeader = "Sequences producing significant alignments:";
    // ...
    // static const char*  kBits = (getenv("CTOOLKIT_COMPATIBLE") ? "(bits)" : "(Bits)");
    // ...
    // static const char*  kValue = "Value";
    // ```
    write_subject_summary_table(writer, &subject_order, &subject_hits, subject_ids)?;

    // Write each subject's hits
    for s_idx in subject_order {
        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
        // ```c
        // typedef struct BlastHSPList {
        //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
        //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
        //                       Set to 0 if not applicable */
        // } BlastHSPList;
        // ```
        let subject_id = subject_ids
            .get(s_idx as usize)
            .map(|id| id.as_ref())
            .unwrap_or("unknown");
        let shits = subject_hits.get(&s_idx).unwrap();
        let first_hit = shits.first().unwrap();

        // Subject header
        write_subject_header(
            writer,
            subject_id,
            first_hit.subject_title.as_deref(),
            first_hit.subject_length,
        )?;

        // Each HSP
        for hit in shits {
            write_hsp_info(writer, hit, config)?;
            write_alignment(writer, hit, config)?;
        }
    }

    writer.flush()
}

/// Write hits in pairwise format (simplified version for Hit without sequences)
///
/// This is a convenience function when only Hit structures are available
/// without the extended sequence data.
pub fn write_pairwise_simple<W: Write>(
    hits: &[Hit],
    writer: &mut W,
    config: &PairwiseConfig,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
    context: &ReportContext,
) -> io::Result<()> {
    let pairwise_hits: Vec<PairwiseHit> = hits.iter().cloned().map(PairwiseHit::from).collect();
    write_pairwise(
        &pairwise_hits,
        writer,
        config,
        query_ids,
        subject_ids,
        context,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn k_bits_follows_ctoolkit_compatible() {
        assert_eq!(ncbi_k_bits_for(true), "(bits)");
        assert_eq!(ncbi_k_bits_for(false), "(Bits)");
        assert_eq!(ncbi_k_bits_for(true).len(), ncbi_k_bits_for(false).len());
    }

    #[test]
    fn ka_value_matches_c_percent_hash_8_3g() {
        // Expected strings from C's sprintf(buffer, "%#8.3g ", value) semantics
        // (Python's '%#.3g' % value gives the same strings).
        let cases: [(f64, &str); 32] = [
            (0.0, "0.00"),
            (-1.0, "-1.00"),
            (0.625, "0.625"),
            (0.41, "0.410"),
            (0.78, "0.780"),
            (1.37, "1.37"),
            (1.28, "1.28"),
            (0.46, "0.460"),
            (0.85, "0.850"),
            (0.99996, "1.00"),
            (9.9996, "10.0"),
            (999.5, "1.00e+03"),
            (1.23e-05, "1.23e-05"),
            (0.0001, "0.000100"),
            (9.999e-05, "0.000100"),
            (9.9994e-05, "0.000100"),
            (1000.0, "1.00e+03"),
            (99.95, "100."),
            (99.94, "99.9"),
            (0.125, "0.125"),
            (0.375, "0.375"),
            (12.5, "12.5"),
            (125.0, "125."),
            (1.5, "1.50"),
            (2.25e-07, "2.25e-07"),
            (123456.0, "1.23e+05"),
            (-0.00042, "-0.000420"),
            (5e-324, "4.94e-324"),
            (f64::MAX, "1.80e+308"),
            (0.0009995, "0.000999"),
            (0.00099949, "0.000999"),
            (-999.49, "-999."),
        ];
        for (value, expected) in cases {
            assert_eq!(format_ncbi_ka_value(value), expected, "value {value:e}");
            let mut field = Vec::new();
            write_ncbi_ka_field(&mut field, value).unwrap();
            assert_eq!(String::from_utf8(field).unwrap(), format!("{expected:>8} "));
        }
        assert_eq!(format_ncbi_ka_value(f64::NAN), "nan");
        assert_eq!(format_ncbi_ka_value(f64::INFINITY), "inf");
        assert_eq!(format_ncbi_ka_value(f64::NEG_INFINITY), "-inf");
    }

    fn make_hit() -> Hit {
        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
        // ```c
        // typedef struct BlastHSPList {
        //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
        //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
        //                       Set to 0 if not applicable */
        //    BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
        //    Int4 hspcnt; /**< Number of HSPs saved */
        //    ...
        // } BlastHSPList;
        // ```
        Hit {
            identity: 95.0,
            length: 100,
            mismatch: 5,
            gapopen: 0,
            q_start: 1,
            q_end: 100,
            s_start: 1,
            s_end: 100,
            e_value: 1e-50,
            bit_score: 185.5,
            num_ident: 95,
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
            // ```c
            // if (hsp->query.frame != hsp->subject.frame) {
            //    *q_end = query_length - hsp->query.offset;
            //    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
            // }
            // ```
            query_frame: 1,
            query_length: 0,
            q_idx: 0,
            s_idx: 0,
            raw_score: 200,
            // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-143
            // BlastHSP stores query/subject BlastSeg offsets used by HSP comparators.
            sort_query_offset: 0,
            sort_query_end: 0,
            sort_subject_offset: 0,
            sort_subject_end: 0,
            has_sort_offsets: false,
            gap_info: None,
            num_positives: 95,
        }
    }

    // NCBI reference: c++/src/corelib/ncbistr.cpp:5088-5340
    // enum EScore { eForced, ePunct, eComma, eSpace, eNewline };
    #[test]
    fn flatfile_wrapping_preserves_ncbi_break_boundaries() {
        for (input, width, expected) in [
            ("abcd", 4, "abcd\n"),
            ("abc-def", 4, "abc-\ndef\n"),
            ("a b c", 4, "a b\nc\n"),
            ("abcd ef", 4, "abcd\nef\n"),
            ("a\n\nb", 4, "a\n\nb\n"),
            ("abc\u{8}d", 4, "abd\n"),
            ("ab\u{8}cd\nz", 68, "ab\u{8}cd\nz\n"),
            ("abc,def", 4, "abc,\ndef\n"),
        ] {
            let mut output = Vec::new();
            write_flatfile_wrapped(&mut output, input, width).unwrap();
            assert_eq!(output, expected.as_bytes(), "{input:?}");
        }
    }

    #[test]
    fn test_write_subject_header() {
        let mut output = Vec::new();
        write_subject_header(&mut output, "seq1", Some("Test sequence"), Some(500)).unwrap();
        let output_str = String::from_utf8(output).unwrap();

        assert!(output_str.contains("> seq1 Test sequence"));
        assert!(output_str.contains("Length=500"));
    }

    #[test]
    fn test_write_hsp_info() {
        let hit = PairwiseHit::from(make_hit());
        let config = PairwiseConfig::default();

        let mut output = Vec::new();
        write_hsp_info(&mut output, &hit, &config).unwrap();
        let output_str = String::from_utf8(output).unwrap();

        assert!(output_str.contains("Score ="));
        assert!(output_str.contains("Expect ="));
        assert!(output_str.contains("Identities ="));
    }

    // NCBI reference: c++/src/objtools/align_format/showalign.cpp:4017-4018
    // ```c++
    //                               m_AV->StrandSign(0),
    //                               m_AV->StrandSign(1),
    // ```
    // A BLASTN HSP of one letter has the same start and end; its strand is its frame's.
    #[test]
    fn blastn_strand_comes_from_the_hsp_frame() {
        let config = PairwiseConfig {
            program: "blastn".to_string(),
            show_frame: false,
            ..PairwiseConfig::default()
        };
        for (query_frame, strand) in [(1, " Strand=Plus/Plus\n"), (-1, " Strand=Plus/Minus\n")] {
            let mut hit = PairwiseHit::from(Hit {
                length: 1,
                s_start: 10,
                s_end: 10,
                query_frame,
                ..make_hit()
            });
            hit.query_seq = Some("G".to_string());
            hit.subject_seq = Some("G".to_string());
            let mut output = Vec::new();
            write_hsp_info(&mut output, &hit, &config).unwrap();
            let output = String::from_utf8(output).unwrap();
            assert!(output.contains(strand), "{output}");
        }
    }

    #[test]
    fn test_write_pairwise_simple() {
        let hits = vec![make_hit()];
        let config = PairwiseConfig::default();
        let context = ReportContext {
            program: "tblastx".to_string(),
            query_name: Some("test_query".to_string()),
            query_length: Some(100),
            subject_name: Some("test_db".to_string()),
            database_name: Some("test_db".to_string()),
            database_num_sequences: Some(1),
            database_total_letters: Some(100),
            version: Some("0.1.0".to_string()),
        };
        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
        // ```c
        // typedef struct BlastHSPList {
        //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
        //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
        //                       Set to 0 if not applicable */
        // } BlastHSPList;
        // ```
        let query_ids = vec![Arc::<str>::from("query1")];
        let subject_ids = vec![Arc::<str>::from("subject1")];

        let mut output = Vec::new();
        write_pairwise_simple(
            &hits,
            &mut output,
            &config,
            &query_ids,
            &subject_ids,
            &context,
        )
        .unwrap();
        let output_str = String::from_utf8(output).unwrap();

        assert!(output_str.contains("TBLASTX"));
        assert!(output_str.contains("Query= test_query"));
        assert!(output_str.contains("> subject1"));
    }
}

// NCBI c++/src/algo/blast/format/blast_format.cpp:348-445,1490-1589:
// PrintProlog(); AcknowledgeBlastQuery(...); x_DisplayDeflines(...);
// display.DisplaySeqalign(...); x_PrintOneQueryFooter(...);
// This uses the local-subject database header and TBLASTN translated alignment.
/// The outfmt 0 prolog of tblastn (the version, the references and the database), which
/// NCBI writes before it reads the first query batch (as `write_tblastx_pairwise_prolog`).
pub fn write_tblastn_pairwise_prolog<W: Write>(
    writer: &mut W,
    version: &str,
    database_name: &str,
    database_num_sequences: usize,
    database_total_letters: usize,
) -> io::Result<()> {
    write_translated_pairwise_intro(writer, "TBLASTN", version)?;
    write_blastp_database_header_spacing(
        writer,
        database_name,
        database_num_sequences,
        database_total_letters,
        1,
    )
}

pub fn write_tblastn_pairwise_report<W: Write>(
    hits: &[PairwiseHit],
    writer: &mut W,
    config: &PairwiseConfig,
    queries: &[BlastpPairwiseQuery],
    query_validity: &[bool],
    query_batch_skipped: &[bool],
    subject_ids: &[Arc<str>],
    report: &BlastpPairwiseReport,
    mut probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut super::query_warnings::QueryWarnings<'_>>,
    epilog: bool,
) -> io::Result<()> {
    // NCBI c++/src/algo/blast/api/local_blast.cpp:177-224:
    // an all-invalid Run() batch carries -1 Karlin sentinel blocks only for
    // its own queries, even if another batch of the input has valid contexts.
    if query_batch_skipped.len() != queries.len() || query_validity.len() != queries.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "TBLASTN batch validity count mismatch",
        ));
    }
    let mut writer = io::BufWriter::new(writer);
    // NCBI c++/src/algo/blast/format/blast_format.cpp:372-424:
    // BlastPrintVersionInfo(m_Program,...); BlastPrintReference(...);
    write_translated_pairwise_intro(&mut writer, "TBLASTN", &report.version)?;
    // The prolog ends with one blank line; each query's preamble starts with two
    // (blast_format.cpp:1491), after the query's warnings.
    write_blastp_database_header_spacing(
        &mut writer,
        &report.database_name,
        report.database_num_sequences,
        report.database_total_letters,
        1,
    )?;

    // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1411
    // ```c++
    // CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
    // ```
    // Each hit keeps its index in `hits` (the final result order), which identifies it
    // to the probe.
    let mut hits_by_query: Vec<Vec<(HspIndex, &PairwiseHit)>> = vec![Vec::new(); queries.len()];
    for (hsp_index, hit) in hits.iter().enumerate() {
        if let Some(bucket) = hits_by_query.get_mut(hit.hit.q_idx as usize) {
            bucket.push((hsp_index, hit));
        }
    }
    for (q_idx, query) in queries.iter().enumerate() {
        // The query's warnings come before its preamble (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut writer)?;
        }
        // NCBI blast_format.cpp:1491: m_Outfile << "\n\n";
        writeln!(&mut writer)?;
        writeln!(&mut writer)?;
        write_blastp_query_header(&mut writer, &query.query_name, query.query_length)?;
        let query_hits = &hits_by_query[q_idx];
        if query_hits.is_empty() {
            write_no_hits_found(&mut writer)?;
            // NCBI blast_results.cpp:82-97 leaves both Karlin blocks null
            // when there is no valid context. blast_format.cpp:445-477
            // writes only the blank lines and zero search space in that case.
            // NCBI c++/src/algo/blast/api/local_blast.cpp:177-180,204-208:
            // if (m_PrelimSearch->CheckInternalData() != 0)
            //     new CBlastAncillaryData(tmp_pair, tmp_pair, tmp_pair, 0);
            // The all-invalid batch uses -1 sentinel blocks. An invalid query
            // within a searched batch has null blocks (blast_results.cpp:82-103).
            if query_batch_skipped[q_idx] {
                write_tblastn_unsearched_query_footer_spacing(&mut writer, false)?;
            } else if !query_validity[q_idx] {
                writeln!(writer)?;
                writeln!(writer)?;
                writeln!(writer)?;
                writeln!(writer, "Effective search space used: 0")?;
            } else {
                write_blastp_query_footer_spacing(
                    &mut writer,
                    query.ungapped_karlin,
                    report.gapped_karlin,
                    report.gumbel,
                    query.effective_search_space,
                    false,
                )?;
            }
            continue;
        }
        let mut subject_hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
            std::collections::HashMap::new();
        let mut subject_hit_indices: std::collections::HashMap<u32, Vec<HspIndex>> =
            std::collections::HashMap::new();
        let mut subject_order = Vec::new();
        for &(hsp_index, hit) in query_hits {
            if !subject_hits.contains_key(&hit.hit.s_idx) {
                subject_order.push(hit.hit.s_idx);
            }
            subject_hits.entry(hit.hit.s_idx).or_default().push(hit);
            subject_hit_indices
                .entry(hit.hit.s_idx)
                .or_default()
                .push(hsp_index);
        }
        // NCBI reference: c++/src/objtools/align_format/showdefline.cpp:1004-1011
        // ```c++
        // void CShowBlastDefline::DisplayBlastDefline(CNcbiOstream & out)
        // {
        //     x_InitDeflineTable();
        //     ...
        //     x_DisplayDefline(out);
        // ```
        // The table of every program follows `x_InitDeflineTable`
        // (`write_blastn_description_table`): the subject's highest bit score and that HSP's
        // E-value, and the title of the subject (`CDeflineGenerator`), for the first
        // `m_NumDescriptions` subjects.
        let described = &subject_order[..subject_order.len().min(report.num_descriptions)];
        write_blastn_description_table(
            &mut writer,
            described,
            subject_order.len() <= report.num_descriptions,
            &subject_hits,
            subject_ids,
            false,
            false,
        )?;
        writeln!(writer)?;
        // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:1014-1040
        // (PruneSeqalign keeps the alignments of the first `m_NumAlignments` subjects).
        let aligned = &subject_order[..subject_order.len().min(report.num_alignments)];
        for &s_idx in aligned {
            let shits = &subject_hits[&s_idx];
            let first = shits[0];
            // NCBI reference: c++/src/objtools/align_format/showalign.cpp:3613-3632
            // ```c++
            //     if(show_defline) {
            // 		...
            // 				string deflines = x_PrintDefLine(bsp_handle, aln_vec_info);
            // 				out<< deflines;
            // 		...
            // 			out << "\n";
            // ```
            // The heading is written before the first HSP of the subject; the probe marks
            // it with that HSP's index without changing the written bytes.
            let first_index = subject_hit_indices[&s_idx][0];
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.subject_begin(first_index);
            }
            // NCBI reference: c++/src/objtools/align_format/showalign.cpp:2273
            // ```c++
            //             alnDispParams->title = CDeflineGenerator().GenerateDefline(bsp_handle);
            // ```
            // The subject (a nucleotide sequence) is shown with its title (`defline.rs`).
            let defline =
                subject_defline(&subject_ids[s_idx as usize], first.subject_title.as_deref());
            let heading = super::defline::ncbi_nucleotide_title(&defline, false);
            write_subject_header(&mut writer, &heading, None, first.subject_length)?;
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.subject_end(first_index);
            }
            // NCBI reference: c++/src/objtools/align_format/showalign.cpp:1970-1973
            // ```c++
            // subid=&(avRef->GetSeqId(1));
            // bool showDefLine = previousId.Empty() || !subid->Match(*previousId);
            // x_DisplayAlnvecInfo(out, alnvecInfo,showDefLine);
            // ```
            // NCBI reference: c++/src/objtools/align_format/showalign.cpp:3971-3980
            // ```c++
            // x_ShowAlnvecInfo(out,aln_vec_info,show_defline);
            // ...
            // out<<"\n";
            // ```
            // The subject defline above is the part of the first x_DisplayAlnvecInfo call
            // that precedes the HSP; the probe marks the rest of each call (the score
            // block, the alignment and the closing blank line) without changing the
            // written bytes.
            for (hit, &hsp_index) in shits.iter().zip(&subject_hit_indices[&s_idx]) {
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.begin(hsp_index);
                }
                write_tblastn_hsp_info(&mut writer, hit)?;
                write_tblastn_alignment(&mut writer, hit, config)?;
                // NCBI c++/src/objtools/align_format/showalign.cpp:3650-3668:
                // display leaves an extra blank line after each translated HSP.
                writeln!(writer)?;
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.end(hsp_index);
                }
            }
        }
        write_blastp_query_footer_spacing(
            &mut writer,
            query.ungapped_karlin,
            report.gapped_karlin,
            report.gumbel,
            query.effective_search_space,
            false,
        )?;
    }
    // NCBI reference: c++/src/app/blast/tblastn_app.cpp:342
    // ```c++
    //         formatter.PrintEpilog(opt);
    // ```
    // A search that stops at a query batch (`Empty CBlastQueryVector` with `-query_loc`)
    // writes no epilog.
    if !epilog {
        return writer.flush();
    }
    // The blank lines before the epilog (blast_format.cpp:2249).
    writeln!(&mut writer)?;
    writeln!(&mut writer)?;
    write_blastp_final_footer(&mut writer, report)?;
    writer.flush()
}

// NCBI c++/src/objtools/align_format/showalign.cpp:3578-3604,2122-2149:
// out << " Score = " << bit_score << " bits (" << score << ")";
// if (comp_adjust) out << ", Method: Compositional matrix adjust.";
// out << " Identities = ..." << " Positives = ..." << " Gaps = ...";
// TBLASTN reports one translated subject frame, without a query frame.
fn write_tblastn_hsp_info<W: Write>(writer: &mut W, hit: &PairwiseHit) -> io::Result<()> {
    let h = &hit.hit;
    let bits = format_bitscore_ncbi(h.bit_score);
    let evalue = format_evalue_ncbi(h.e_value);
    // NCBI c++/src/objtools/align_format/showalign.cpp:3599-3604:
    // if (comp_adj_method == 1) out << ", Method: Composition-based stats.";
    // else if (comp_adj_method == 2) out << ", Method: Compositional matrix adjust.";
    let method = match hit.comp_adjust_method {
        Some(1) => ", Method: Composition-based stats.",
        Some(2) => ", Method: Compositional matrix adjust.",
        _ => "",
    };
    // NCBI c++/src/algo/blast/api/blast_seqalign.cpp:1180-1182:
    // if (hsp->num > 1) scores.push_back(s_MakeScore("sum_n", ..., hsp->num));
    // NCBI c++/src/objtools/align_format/showalign.cpp:3595-3598:
    // if (aln_vec_info->sum_n > 0) out << "(" << aln_vec_info->sum_n << ")";
    let expect = hit
        .sum_n
        .filter(|&count| count > 1)
        .map_or_else(|| "Expect".to_string(), |count| format!("Expect({count})"));
    writeln!(
        writer,
        " Score = {} bits ({}),  {} = {}{}",
        bits, h.raw_score, expect, evalue, method
    )?;
    let length = h.length;
    // NCBI align_format_util.cpp:3152-3161:
    // if (numerator == denominator) return 100;
    // int retval = (int)(0.5 + 100.0*numerator/denominator);
    // return min(99, retval);
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
    // ```c++
    // int CAlignFormatUtil::GetPercentMatch(int numerator, int denominator)
    // {
    //      if (numerator == denominator)
    //         return 100;
    //      else {
    //        int retval =(int) (0.5 + 100.0*((double)numerator)/((double)denominator));
    //        retval = min(99, retval);
    //        return retval;
    //      }
    // }
    // ```
    let percent = |n: usize| ncbi_percent_match(n, length);
    let positives = hit.positives.unwrap_or(h.num_positives);
    let gaps = hit.gaps.unwrap_or_else(|| h.gap_letters());
    writeln!(
        writer,
        " Identities = {}/{} ({:.0}%), Positives = {}/{} ({:.0}%), Gaps = {}/{} ({:.0}%)",
        h.num_ident,
        length,
        percent(h.num_ident),
        positives,
        length,
        percent(positives),
        gaps,
        length,
        percent(gaps)
    )?;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:321-332
    // ```c++
    //     if(master_frame != 0 && slave_frame != 0) {
    //         out <<" Frame = " << ((master_frame > 0) ? "+" : "")
    //             << master_frame <<"/"<<((slave_frame > 0) ? "+" : "")
    //             << slave_frame<<"\n";
    //     } else if (master_frame != 0){
    //         out <<" Frame = " << ((master_frame > 0) ? "+" : "")
    //             << master_frame << "\n";
    //     }  else if (slave_frame != 0){
    //         out <<" Frame = " << ((slave_frame > 0) ? "+" : "")
    //             << slave_frame <<"\n";
    //     }
    //     out<<"\n";
    // ```
    if let Some(frame) = hit.query_frame.or(hit.subject_frame) {
        if frame > 0 {
            writeln!(writer, " Frame = +{frame}")?;
        } else {
            writeln!(writer, " Frame = {frame}")?;
        }
    }
    writeln!(writer)?;
    Ok(())
}

// NCBI c++/src/objtools/align_format/showalign.cpp:1600-1625,2122-2149:
// AddSpace(out,maxStartLen-startLen+k_StartSequenceMargin); out << sequence;
// if (sequence_standard[i]==sequence[i]) middle_line[i]=sequence[i];
// else if (m_Matrix[...]>0) middle_line[i]='+';
// NCBI c++/src/algo/blast/core/blast_hits.c:1087-1105:
// translated subject offsets advance by CODON_LENGTH per aligned residue.
fn write_tblastn_alignment<W: Write>(
    writer: &mut W,
    hit: &PairwiseHit,
    config: &PairwiseConfig,
) -> io::Result<()> {
    let (Some(qseq), Some(sseq)) = (&hit.query_seq, &hit.subject_seq) else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing TBLASTN alignment strings",
        ));
    };
    let q = qseq.as_bytes();
    let s = sseq.as_bytes();
    if q.len() != s.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "TBLASTN alignment string lengths differ",
        ));
    }
    let width = coordinate_width(&hit.hit);
    let direction: isize = if hit.subject_frame.unwrap_or(1) < 0 {
        -1
    } else {
        1
    };
    let mut q_pos = hit.hit.q_start;
    let mut s_pos = hit.hit.s_start as isize;
    for (qpart, spart) in q
        .chunks(config.line_length)
        .zip(s.chunks(config.line_length))
    {
        let q_count = qpart.iter().filter(|&&c| c != b'-').count();
        let s_count = spart.iter().filter(|&&c| c != b'-').count();
        let q_end = q_pos + q_count.saturating_sub(1);
        let s_end = s_pos + direction * (3 * s_count as isize - 1);
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1605-1627
        // ```c++
        //     //not to display start and stop number for empty row
        //     if ((j > 0 && end == prev_stop)
        //         || (j == 0 && start == 1 && end == 1)) {
        //         startLen = 0;
        //     } else {
        //         out << start;
        // ...
        //      //not to display stop number for empty row in the middle
        //     if (!(j > 0 && end == prev_stop)
        // ```
        // A row of gaps only shows no coordinates (as `write_blastx_alignment`).
        write!(writer, "Query  ")?;
        if q_count > 0 {
            write!(writer, "{q_pos}")?;
            write_spaces(writer, width + 2 - digit_count(q_pos))?;
        } else {
            write_spaces(writer, width + 2)?;
        }
        writer.write_all(qpart)?;
        write!(writer, "  ")?;
        if q_count > 0 {
            write!(writer, "{q_end}")?;
        }
        writeln!(writer)?;
        write_spaces(writer, 7 + width + 2)?;
        for (&qc, &sc) in qpart.iter().zip(spart) {
            // NCBI c++/src/objtools/align_format/showalign.cpp:2122-2149:
            // the displayed query preserves lowercase masking, while the
            // protein matrix comparison uses its uppercase residue.
            let query_residue = qc.to_ascii_uppercase();
            let mid = if query_residue == sc {
                query_residue
            } else if is_positive_match(query_residue as char, sc as char, config.protein_matrix) {
                b'+'
            } else {
                b' '
            };
            writer.write_all(&[mid])?;
        }
        writeln!(writer)?;
        write!(writer, "Sbjct  ")?;
        if s_count > 0 {
            write!(writer, "{s_pos}")?;
            write_spaces(writer, width + 2 - digit_count(s_pos.unsigned_abs()))?;
        } else {
            write_spaces(writer, width + 2)?;
        }
        writer.write_all(spart)?;
        write!(writer, "  ")?;
        if s_count > 0 {
            write!(writer, "{s_end}")?;
        }
        writeln!(writer)?;
        writeln!(writer)?;
        q_pos += q_count;
        s_pos += direction * 3 * s_count as isize;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:387-407
// ```c++
//       CBlastFormatUtil::BlastPrintVersionInfo(m_Program, m_IsHTML,
//                                               m_Outfile);
//     }
//
//     if (m_IsBl2Seq && !m_IsDbScan) {
//         return;
//     }
//
//     m_Outfile << NcbiEndl << NcbiEndl;
//     if (m_Program == "deltablast") {
//         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
//                               m_Outfile, CReference::eDeltaBlast);
//         m_Outfile << "\n";
//     }
//
//     if (m_Megablast)
//         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
//                                           m_Outfile, CReference::eMegaBlast);
//     else
//         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
//                                           m_Outfile);
// ```
fn write_translated_pairwise_intro(
    writer: &mut impl Write,
    program: &str,
    version: &str,
) -> io::Result<()> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:66-73
    // ```c++
    // void CBlastFormatUtil::BlastPrintVersionInfo(const string program, bool html,
    //                                              CNcbiOstream& out)
    // {
    //     if (html)
    //         out << "<b>" << BlastGetVersion(program) << "</b>" << "\n";
    //     else
    //         out << BlastGetVersion(program) << "\n";
    // }
    // ```
    writeln!(writer, "{} {}", program, version)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:395-395
    // ```c++
    //     m_Outfile << NcbiEndl << NcbiEndl;
    // ```
    writeln!(writer)?;
    writer.flush()?;
    writeln!(writer)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:395-395
    // ```c++
    //     m_Outfile << NcbiEndl << NcbiEndl;
    // ```
    writer.flush()?;
    writeln!(
        writer,
        "Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A."
    )?;
    writeln!(
        writer,
        "Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J."
    )?;
    writeln!(
        writer,
        "Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of"
    )?;
    writeln!(
        writer,
        "protein database search programs\", Nucleic Acids Res. 25:3389-3402."
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)?;
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:275-282
// ```c++
//     if ((m_FormatType == CFormattingArgs::eXml2) || (m_FormatType == CFormattingArgs::eJson) ||
//         (m_FormatType == CFormattingArgs::eXml2_S) || (m_FormatType == CFormattingArgs::eJson_S)) {
//            m_AccumulatedQueries.Reset(new CBlastQueryVector());
//     }
//
//     if (opts.GetSumStatisticsMode() && m_IsUngappedSearch) {
//         m_ShowLinkedSetSize = true;
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1548-1565
// ```c++
//     results.GetMaskedQueryRegions(masklocs);
//
//     CSeq_align_set copy_aln_set;
//     CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_NumAlignments);
//
//     int flags = s_SetFlags(m_Program, m_FormatType, m_IsHTML, m_ShowGi,
//                            (m_IsBl2Seq && !m_IsDbScan), (m_DisableKAStats || kIsGlobal));
//
//     CDisplaySeqalign display(copy_aln_set, *m_Scope, &masklocs, NULL, m_MatrixName);
//     display.SetDbName(m_DbName);
//     display.SetDbType(!m_DbIsAA);
//     display.SetLineLen(m_LineLength);
//     int kAlignToShow=2000000000;  // Nice large number per SB-1817
//     display.SetNumAlignToShow(kAlignToShow);
//
//     // set the alignment flags
//     display.SetAlignOption(flags);
//
// ```
pub struct BlastxPairwiseOptions {
    pub word_threshold: f64,
    pub gapped: bool,
    pub sum_stats: bool,
    pub num_descriptions: usize,
    pub num_alignments: usize,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1490-1519
// ```c++
//
//     m_Outfile << "\n\n";
//     CBlastFormatUtil::AcknowledgeBlastQuery(*bioseq, kFormatLineLength,
//                                             m_Outfile, m_BelieveQuery,
//                                             m_IsHTML, kIsTabularOutput,
//                                             results.GetRID());
//
//     if (m_IsBl2Seq && !m_IsDbScan) {
//         m_Outfile << "\n";
//         // FIXME: this might be configurable in the future
//         const bool kBelieveSubject = false;
//         CConstRef<CBioseq> subject_bioseq = x_CreateSubjectBioseq();
//         CBlastFormatUtil::AcknowledgeBlastSubject(*subject_bioseq,
//                                                   kFormatLineLength,
//                                                   m_Outfile, kBelieveSubject,
//                                                   m_IsHTML, kIsTabularOutput);
//     }
//
//     // quit early if there are no hits
//     if ( !results.HasAlignments() ) {
//         m_Outfile << "\n\n"
//               << "***** " << CBlastFormatUtil::kNoHitsFound << " *****" << "\n"
//               << "\n\n";
//         x_PrintOneQueryFooter(*results.GetAncillaryData());
//         return;
//     }
//
//     CConstRef<CSeq_align_set> aln_set = results.GetSeqAlign();
//     _ASSERT(results.HasAlignments());
//     if (m_IsUngappedSearch) {
// ```
#[allow(clippy::too_many_arguments)]
pub fn write_blastx_pairwise_report<W: Write>(
    hits: &[PairwiseHit],
    writer: &mut W,
    stderr: &mut dyn Write,
    warnings: &[Vec<u8>],
    config: &PairwiseConfig,
    queries: &[BlastpPairwiseQuery],
    query_validity: &[bool],
    skipped: &[bool],
    subject_ids: &[Arc<str>],
    report: &BlastpPairwiseReport,
    options: &BlastxPairwiseOptions,
) -> io::Result<()> {
    let mut writer = io::BufWriter::new(writer);
    let mut by_query: Vec<Vec<&PairwiseHit>> = vec![Vec::new(); queries.len()];
    for hit in hits {
        by_query[hit.hit.q_idx as usize].push(hit);
    }
    for (qidx, query) in queries.iter().enumerate() {
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452
        // ```c++
        //     if (results.HasWarnings()) {
        //         ERR_POST(Warning << results.GetWarningStrings());
        //     }
        // ```
        stderr.write_all(&warnings[qidx])?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1491-1492
        // ```c++
        //     m_Outfile << "\n\n";
        //     CBlastFormatUtil::AcknowledgeBlastQuery(*bioseq, kFormatLineLength,
        // ```
        writeln!(writer)?;
        writeln!(writer)?;
        write_blastp_query_header(&mut writer, &query.query_name, query.query_length)?;
        let query_hits = &by_query[qidx];
        if query_hits.is_empty() {
            write_no_hits_found(&mut writer)?;
        } else {
            let mut grouped = std::collections::HashMap::new();
            let mut order = Vec::new();
            for hit in query_hits {
                if !grouped.contains_key(&hit.hit.s_idx) {
                    order.push(hit.hit.s_idx);
                }
                grouped
                    .entry(hit.hit.s_idx)
                    .or_insert_with(Vec::new)
                    .push(*hit);
            }
            let descriptions = &order[..order.len().min(options.num_descriptions)];
            write_subject_summary_table_with_sum_n(
                &mut writer,
                descriptions,
                &grouped,
                subject_ids,
                !options.gapped && options.sum_stats,
            )?;
            writeln!(writer)?;
            for oid in order.iter().take(options.num_alignments) {
                let shits = &grouped[oid];
                let first = shits[0];
                write_subject_header(
                    &mut writer,
                    &subject_ids[*oid as usize],
                    first.subject_title.as_deref(),
                    first.subject_length,
                )?;
                for hit in shits {
                    write_tblastn_hsp_info(&mut writer, hit)?;
                    write_blastx_alignment(&mut writer, hit, config)?;
                    writeln!(writer)?;
                }
            }
        }
        write_blastx_query_footer(
            &mut writer,
            query,
            report,
            query_validity[qidx],
            skipped[qidx],
            options.gapped,
        )?;
    }
    writer.flush()
}
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254-254
// ```c++
//         formatter.PrintProlog();
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:297-297
// ```c++
//         formatter.PrintEpilog(opt);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:387-388
// ```c++
//       CBlastFormatUtil::BlastPrintVersionInfo(m_Program, m_IsHTML,
//                                               m_Outfile);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:406-407
// ```c++
//         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
//                                           m_Outfile);
// ```
pub fn write_blastx_prolog(
    writer: &mut impl Write,
    report: &BlastpPairwiseReport,
) -> io::Result<()> {
    write_translated_pairwise_intro(writer, "BLASTX", &report.version)?;
    write_blastp_database_header_spacing(
        writer,
        &report.database_name,
        report.database_num_sequences,
        report.database_total_letters,
        1,
    )
}
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2264-2288
// ```c++
//                         options.GetMismatchPenalty() << "\n";
//     }
//     else {
//         m_Outfile << "\n\nMatrix: " << options.GetMatrixName() << "\n";
//     }
//
//     if (options.GetGappedMode() == true) {
//         double gap_extension = (double) options.GetGapExtensionCost();
//         if ((m_Program == "megablast" || m_Program == "blastn") && options.GetGapExtensionCost() == 0)
//         { // Formula from PMID 10890397 applies if both gap values are zero.
//                gap_extension = -2*options.GetMismatchPenalty() + options.GetMatchReward();
//                gap_extension /= 2.0;
//         }
//         m_Outfile << "Gap Penalties: Existence: "
//                 << options.GetGapOpeningCost() << ", Extension: "
//                 << gap_extension << "\n";
//     }
//     if (options.GetWordThreshold()) {
//         m_Outfile << "Neighboring words threshold: " <<
//                         options.GetWordThreshold() << "\n";
//     }
//     if (options.GetWindowSize()) {
//         m_Outfile << "Window for multiple hits: " <<
//                         options.GetWindowSize() << "\n";
//     }
// ```
pub fn write_blastx_epilog(
    writer: &mut impl Write,
    report: &BlastpPairwiseReport,
    options: &BlastxPairwiseOptions,
) -> io::Result<()> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2264-2288
    // ```c++
    //                         options.GetMismatchPenalty() << "\n";
    //     }
    //     else {
    //         m_Outfile << "\n\nMatrix: " << options.GetMatrixName() << "\n";
    //     }
    //
    //     if (options.GetGappedMode() == true) {
    //         double gap_extension = (double) options.GetGapExtensionCost();
    //         if ((m_Program == "megablast" || m_Program == "blastn") && options.GetGapExtensionCost() == 0)
    //         { // Formula from PMID 10890397 applies if both gap values are zero.
    //                gap_extension = -2*options.GetMismatchPenalty() + options.GetMatchReward();
    //                gap_extension /= 2.0;
    //         }
    //         m_Outfile << "Gap Penalties: Existence: "
    //                 << options.GetGapOpeningCost() << ", Extension: "
    //                 << gap_extension << "\n";
    //     }
    //     if (options.GetWordThreshold()) {
    //         m_Outfile << "Neighboring words threshold: " <<
    //                         options.GetWordThreshold() << "\n";
    //     }
    //     if (options.GetWindowSize()) {
    //         m_Outfile << "Window for multiple hits: " <<
    //                         options.GetWindowSize() << "\n";
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2249-2249
    // ```c++
    //     m_Outfile << NcbiEndl << NcbiEndl;
    // ```
    writeln!(writer)?;
    writeln!(writer)?;
    write!(writer, "  Database: ")?;
    write_flatfile_wrapped(writer, &ensure_trailing_period(&report.database_name), 68)?;
    writeln!(writer, "    Posted date:  Unknown")?;
    writeln!(
        writer,
        "  Number of letters in database: {}",
        format_count_with_commas(report.database_total_letters as i64)
    )?;
    writeln!(
        writer,
        "  Number of sequences in database:  {}",
        format_count_with_commas(report.database_num_sequences as i64)
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer, "Matrix: {}", report.matrix_name)?;
    if options.gapped {
        writeln!(
            writer,
            "Gap Penalties: Existence: {}, Extension: {}",
            report.gap_open, report.gap_extend
        )?;
    }
    if options.word_threshold != 0.0 {
        writeln!(
            writer,
            "Neighboring words threshold: {}",
            options.word_threshold
        )?;
    }
    if report.window_size != 0 {
        writeln!(writer, "Window for multiple hits: {}", report.window_size)?;
    }
    Ok(())
}
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:445-478
// ```c++
// CBlastFormat::x_PrintOneQueryFooter(const blast::CBlastAncillaryData& summary)
// {
//     /* Skip printing KA parameters if the program is rmblastn -RMH- */
//     if ( m_DisableKAStats )
//       return;
//
//     const Blast_KarlinBlk *kbp_ungap =
//         (m_Program == "psiblast" || m_Program == "deltablast")
//         ? summary.GetPsiUngappedKarlinBlk()
//         : summary.GetUngappedKarlinBlk();
//     const Blast_GumbelBlk *gbp = summary.GetGumbelBlk();
//     m_Outfile << NcbiEndl;
//     if (kbp_ungap) {
//         CBlastFormatUtil::PrintKAParameters(kbp_ungap->Lambda,
//                                             kbp_ungap->K, kbp_ungap->H,
//                                             kFormatLineLength, m_Outfile,
//                                             false, gbp);
//     }
//
//     const Blast_KarlinBlk *kbp_gap =
//         (m_Program == "psiblast" || m_Program == "deltablast")
//         ? summary.GetPsiGappedKarlinBlk()
//         : summary.GetGappedKarlinBlk();
//     m_Outfile << "\n";
//     if (kbp_gap) {
//         CBlastFormatUtil::PrintKAParameters(kbp_gap->Lambda,
//                                             kbp_gap->K, kbp_gap->H,
//                                             kFormatLineLength, m_Outfile,
//                                             true, gbp);
//     }
//
//     m_Outfile << "\n";
//     m_Outfile << "Effective search space used: " <<
//                         summary.GetSearchSpace() << "\n";
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:583-613
// ```c++
//
//     char buffer[256];
//     if (gapped) {
//         out << "Gapped" << "\n";
//     }
//     out << "Lambda      K        H";
//     if (gbp) {
//         if (gapped) {
//             out << "        a         alpha    sigma";
//         } else {
//             out << "        a         alpha";
//         }
//     }
//     out << "\n";
//     sprintf(buffer, "%#8.3g ", lambda);
//     out << buffer;
//     sprintf(buffer, "%#8.3g ", k);
//     out << buffer;
//     sprintf(buffer, "%#8.3g ", h);
//     out << buffer;
//     if (gbp) {
//         if (gapped) {
//             sprintf(buffer, "%#8.3g ", gbp->a);
//             out << buffer;
//             sprintf(buffer, "%#8.3g ", gbp->Alpha);
//             out << buffer;
//             sprintf(buffer, "%#8.3g ", gbp->Sigma);
//             out << buffer;
//         } else {
//             sprintf(buffer, "%#8.3g ", gbp->a_un);
//             out << buffer;
// ```
fn write_blastx_query_footer(
    writer: &mut impl Write,
    query: &BlastpPairwiseQuery,
    report: &BlastpPairwiseReport,
    valid: bool,
    skipped: bool,
    gapped: bool,
) -> io::Result<()> {
    if skipped {
        return write_tblastn_unsearched_query_footer_spacing(writer, false);
    }
    if valid && gapped {
        return write_blastp_query_footer_spacing(
            writer,
            query.ungapped_karlin,
            report.gapped_karlin,
            report.gumbel,
            query.effective_search_space,
            false,
        );
    }
    writeln!(writer)?;
    if valid {
        writeln!(writer, "Lambda      K        H")?;
        for v in [
            query.ungapped_karlin.lambda,
            query.ungapped_karlin.k,
            query.ungapped_karlin.h,
        ] {
            write_ncbi_ka_field(writer, v)?;
        }
        writeln!(writer)?;
    }
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Effective search space used: {}",
        if valid {
            query.effective_search_space
        } else {
            0
        }
    )?;
    Ok(())
}
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1596-1627
// ```c++
// {
//     size_t startLen=0;
//     int start = alnRoInfo->seqStarts[row].front() + 1;  //+1 for 1 based
//     int end = alnRoInfo->seqStops[row].front() + 1;
//     int j = alnRoInfo->currPrintSegment;
//     int actualLineLen = alnRoInfo->currActualLineLen;
//     //print out sequence line
//     //adjust space between id and start
//     CAlignFormatUtil::AddSpace(out, alnRoInfo->maxIdLen-alnRoInfo->seqidArray[row].size() + k_IdStartMargin);
//     //not to display start and stop number for empty row
//     if ((j > 0 && end == prev_stop)
//         || (j == 0 && start == 1 && end == 1)) {
//         startLen = 0;
//     } else {
//         out << start;
//         startLen=NStr::IntToString(start).size();
//     }
//
//     CAlignFormatUtil::AddSpace(out, alnRoInfo->maxStartLen-startLen + k_StartSequenceMargin);
//     x_OutputSeq(alnRoInfo->sequence[row], m_AV->GetSeqId(row), j,
//                 (int)actualLineLen, alnRoInfo->frame[row], row,
//                 (row > 0 && alnRoInfo->colorMismatch)?true:false,
//                 alnRoInfo->masked_regions[row], out);
//     CAlignFormatUtil::AddSpace(out, k_SeqStopMargin);
//
//      //not to display stop number for empty row in the middle
//     if (!(j > 0 && end == prev_stop)
//         && !(j == 0 && start == 1 && end == 1)) {
//         out << end;
//     }
//     out<<"\n";
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:120-134
// ```c++
//
//     if (strand == eNa_strand_minus) {
//
//         if (translate)
//             retval = original_length -
//                 CODON_LENGTH*(s_GetCurrPos(curr_pos, num) + num)
//                 + frame + 1;
//         else
//             retval = length - s_GetCurrPos(curr_pos, num) - num;
//
//     } else {
//
//         if (translate)
//             retval = frame - 1 + CODON_LENGTH*s_GetCurrPos(curr_pos, num);
//         else
// ```
fn write_blastx_alignment(
    writer: &mut impl Write,
    hit: &PairwiseHit,
    config: &PairwiseConfig,
) -> io::Result<()> {
    let (Some(qseq), Some(sseq)) = (&hit.query_seq, &hit.subject_seq) else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing BLASTX alignment strings",
        ));
    };
    let q = qseq.as_bytes();
    let s = sseq.as_bytes();
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1364-1367
    // ```c++
    //         x_FillSeqid(seqidArray[row], row);
    //         maxIdLen=max<size_t>(seqidArray[row].size(), maxIdLen);
    //         size_t maxCood=max<size_t>(m_AV->GetSeqStart(row), m_AV->GetSeqStop(row));
    //         maxStartLen = max<size_t>(NStr::SizetToString(maxCood).size(), maxStartLen);
    // ```
    let width = digit_count(
        hit.hit
            .q_start
            .max(hit.hit.q_end)
            .max(hit.hit.s_start)
            .max(hit.hit.s_end)
            - 1,
    );
    let direction: isize = if hit.query_frame.unwrap_or(1) < 0 {
        -1
    } else {
        1
    };
    let mut qp = hit.hit.q_start as isize;
    let mut sp = hit.hit.s_start;
    for (qpart, spart) in q
        .chunks(config.line_length)
        .zip(s.chunks(config.line_length))
    {
        let qn = qpart.iter().filter(|&&c| c != b'-').count();
        let sn = spart.iter().filter(|&&c| c != b'-').count();
        let qe = qp + direction * (3 * qn as isize - 1);
        let se = sp + sn.saturating_sub(1);
        write!(writer, "Query  ")?;
        if qn > 0 {
            write!(writer, "{qp}")?;
            write_spaces(writer, width + 2 - digit_count(qp.unsigned_abs()))?;
        } else {
            write_spaces(writer, width + 2)?;
        }
        writer.write_all(qpart)?;
        write!(writer, "  ")?;
        if qn > 0 {
            write!(writer, "{qe}")?;
        }
        writeln!(writer)?;
        write_spaces(writer, 7 + width + 2)?;
        for (&qc, &sc) in qpart.iter().zip(spart) {
            let qu = qc.to_ascii_uppercase();
            let middle = if qu == sc {
                qu
            } else if is_positive_match(qu as char, sc as char, config.protein_matrix) {
                b'+'
            } else {
                b' '
            };
            writer.write_all(&[middle])?;
        }
        writeln!(writer)?;
        write!(writer, "Sbjct  ")?;
        if sn > 0 {
            write!(writer, "{sp}")?;
            write_spaces(writer, width + 2 - digit_count(sp))?;
        } else {
            write_spaces(writer, width + 2)?;
        }
        writer.write_all(spart)?;
        write!(writer, "  ")?;
        if sn > 0 {
            write!(writer, "{se}")?;
        }
        writeln!(writer)?;
        writeln!(writer)?;
        qp += direction * 3 * qn as isize;
        sp += sn;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
// ```c++
// int CAlignFormatUtil::GetPercentMatch(int numerator, int denominator)
// {
//      if (numerator == denominator)
//         return 100;
//      else {
//        int retval =(int) (0.5 + 100.0*((double)numerator)/((double)denominator));
//        retval = min(99, retval);
//        return retval;
//      }
// }
// ```
#[cfg(test)]
#[test]
fn percent_match_matches_independent_pinned_cpp_edges() {
    let expected = include_str!("../../tests/unit/blastx_stage_e_percent_expected.tsv");
    for line in expected.lines() {
        let fields: Vec<usize> = line
            .split('\t')
            .map(|value| value.parse().unwrap())
            .collect();
        assert_eq!(fields.len(), 3);
        assert_eq!(
            ncbi_percent_match(fields[0], fields[1]),
            fields[2],
            "{line}"
        );
    }
}

#[cfg(test)]
mod tblastx_tests {
    use super::*;

    #[test]
    fn thresholds_are_written_as_a_cpp_stream_writes_a_double() {
        for (value, text) in [
            (11.0, "11"),
            (999999.0, "999999"),
            (1000000.0, "1e+06"),
            (1234567.0, "1.23457e+06"),
            (2147483647.0, "2.14748e+09"),
            (12.5, "12.5"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
        ] {
            assert_eq!(cpp_default_double(value), text, "{value}");
        }
    }

    fn tblastx_hit(sum_n: i32) -> PairwiseHit {
        let query: String = "ACDEFGHIKL".repeat(6) + "M";
        let subject: String = "ACDEFGHIKL".repeat(6) + "N";
        PairwiseHit {
            hit: Hit {
                identity: 0.0,
                length: 61,
                mismatch: 1,
                gapopen: 0,
                q_start: 200,
                q_end: 18,
                s_start: 3,
                s_end: 185,
                e_value: 2e-30,
                bit_score: 120.5,
                num_ident: 60,
                query_frame: -2,
                query_length: 0,
                q_idx: 0,
                s_idx: 0,
                raw_score: 255,
                sort_query_offset: 0,
                sort_query_end: 61,
                sort_subject_offset: 0,
                sort_subject_end: 61,
                has_sort_offsets: true,
                gap_info: None,
                num_positives: 60,
            },
            query_seq: Some(query),
            subject_seq: Some(subject),
            query_frame: Some(-2),
            subject_frame: Some(3),
            positives: Some(60),
            gaps: Some(0),
            subject_length: Some(400),
            subject_title: None,
            comp_adjust_method: None,
            sum_n: Some(sum_n),
        }
    }

    // NCBI c++/src/objtools/align_format/showalign.cpp:3593-3606 and 310-332: Expect(n)
    // only for a linked set of two or more, both frames, no Strand line.
    #[test]
    fn tblastx_hsp_info_prints_the_linked_set_and_both_frames() {
        let mut out = Vec::new();
        write_tblastx_hsp_info(&mut out, &tblastx_hit(4)).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            " Score = 120 bits (255),  Expect(4) = 2e-30\n Identities = 60/61 (98%), Positives = 60/61 (98%), Gaps = 0/61 (0%)\n Frame = -2/+3\n\n"
        );
        let mut out = Vec::new();
        write_tblastx_hsp_info(&mut out, &tblastx_hit(1)).unwrap();
        assert!(String::from_utf8(out)
            .unwrap()
            .starts_with(" Score = 120 bits (255),  Expect = 2e-30\n"));
    }

    // NCBI c++/src/objtools/alnmgr/alnvec.cpp:230,305,358 and showalign.cpp:1598-1626: 60
    // residues span 180 nucleotides, each row in the direction of its frame, and the width
    // is that of the largest 0-based coordinate (199: 3 digits).
    #[test]
    fn tblastx_rows_step_three_nucleotides_in_each_frame_direction() {
        let mut out = Vec::new();
        write_tblastx_alignment(&mut out, &tblastx_hit(1), &PairwiseConfig::default()).unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        let row = "ACDEFGHIKL".repeat(6);
        assert_eq!(lines[0], format!("Query  200  {row}  21"));
        assert_eq!(lines[1], format!("            {row}"));
        assert_eq!(lines[2], format!("Sbjct  3    {row}  182"));
        assert_eq!(lines[3], "");
        assert_eq!(lines[4], "Query  20   M  18");
        assert_eq!(lines[5], "             ");
        assert_eq!(lines[6], "Sbjct  183  N  185");
    }

    // NCBI c++/src/objtools/align_format/showdefline.cpp:829-836,961-965: the N column of
    // an ungapped sum-statistics search, from the first HSP of each subject.
    #[test]
    fn description_table_shows_the_linked_set_column() {
        let first = tblastx_hit(12);
        let mut second = tblastx_hit(1);
        second.hit.s_idx = 1;
        second.hit.bit_score = 30.2;
        second.hit.e_value = 0.5;
        let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
            [(0, vec![&first]), (1, vec![&second])]
                .into_iter()
                .collect();
        let ids: Vec<Arc<str>> = vec![Arc::from("s0"), Arc::from("s1")];
        let mut out = Vec::new();
        write_blastn_description_table(&mut out, &[0, 1], true, &hits, &ids, true, false).unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert!(lines[0].ends_with("Score     E"), "{text}");
        assert!(lines[1].ends_with("(Bits)  Value  N"), "{text}");
        assert!(lines[3].ends_with("  120     2e-30  12"), "{text}");
        assert!(lines[4].ends_with("  30.2    0.50   1 "), "{text}");
    }
}
