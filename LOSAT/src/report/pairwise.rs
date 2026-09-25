//! Pairwise alignment output (outfmt 0)
//!
//! This module implements the traditional BLAST pairwise alignment output format.
//!
//! Reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp

use super::outfmt6::{format_bitscore_ncbi, format_evalue_ncbi, ReportContext};
use crate::common::Hit;
use crate::config::ScoringMatrix;
use crate::stats::{BlastGumbelBlk, KarlinParams};
use crate::utils::matrix::{aa_char_to_ncbistdaa, protein_score};
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
    pub word_threshold: i32,
    pub window_size: i32,
    pub gapped_karlin: KarlinParams,
    pub gumbel: BlastGumbelBlk,
}

// =============================================================================
// Pairwise Output Writers
// =============================================================================

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
        // NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:3578-3604
        // ```c
        // out << " Score = " << bit_score
        //     << " bits (" << score << "),  Expect = " << evalue;
        // ...
        // out << ", Method: Compositional matrix adjust.";
        // ```
        writeln!(
            writer,
            " Score = {} bits ({}),  Expect = {}, Method: Compositional matrix adjust.",
            bit_score_str, h.raw_score, evalue_str
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
    let ident_pct = if align_len > 0 {
        (num_ident as f64 / align_len as f64) * 100.0
    } else {
        0.0
    };

    let positives = hit.positives.unwrap_or(h.num_positives.max(num_ident));
    let pos_pct = if align_len > 0 {
        (positives as f64 / align_len as f64) * 100.0
    } else {
        0.0
    };

    let gaps = hit.gaps.unwrap_or_else(|| h.gap_letters());
    let gap_pct = if align_len > 0 {
        (gaps as f64 / align_len as f64) * 100.0
    } else {
        0.0
    };

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
    let max_pos = hit.q_start.max(hit.q_end).max(hit.s_start).max(hit.s_end);
    let max_start_len = format!("{max_pos}").len();
    let max_id_len = "Query".len().max("Sbjct".len());

    let mut q_pos = hit.q_start;
    let mut s_pos = hit.s_start;
    let mut offset = 0;

    while offset < q_chars.len() {
        let end = (offset + line_len).min(q_chars.len());
        let chunk_q: String = q_chars[offset..end].iter().collect();
        let chunk_s: String = s_chars[offset..end].iter().collect();

        // Build middle line (identity markers)
        let middle: String = q_chars[offset..end]
            .iter()
            .zip(s_chars[offset..end].iter())
            .map(|(q, s)| {
                if q == s {
                    *q // Identity: show the character
                } else if is_positive_match(*q, *s, config.protein_matrix) {
                    '+' // Positive: show +
                } else {
                    ' ' // Mismatch: show space
                }
            })
            .collect();

        // Count non-gap characters in this chunk
        let q_non_gap = chunk_q.chars().filter(|&c| c != '-').count();
        let s_non_gap = chunk_s.chars().filter(|&c| c != '-').count();

        let q_end_pos = q_pos + q_non_gap.saturating_sub(1);
        let s_end_pos = s_pos + s_non_gap.saturating_sub(1);

        write!(writer, "Query")?;
        write_spaces(writer, max_id_len - "Query".len() + 2)?;
        write!(writer, "{}", q_pos)?;
        write_spaces(writer, max_start_len - digit_count(q_pos) + 2)?;
        write!(writer, "{}", chunk_q)?;
        write_spaces(writer, 2)?;
        writeln!(writer, "{}", q_end_pos)?;

        write_spaces(writer, max_id_len + 2 + max_start_len + 2)?;
        writeln!(writer, "{}", middle)?;

        write!(writer, "Sbjct")?;
        write_spaces(writer, max_id_len - "Sbjct".len() + 2)?;
        write!(writer, "{}", s_pos)?;
        write_spaces(writer, max_start_len - digit_count(s_pos) + 2)?;
        write!(writer, "{}", chunk_s)?;
        write_spaces(writer, 2)?;
        writeln!(writer, "{}", s_end_pos)?;

        writeln!(writer)?;

        q_pos = q_end_pos + 1;
        s_pos = s_end_pos + 1;
        offset = end;
    }

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

    let a_ncbi = aa_char_to_ncbistdaa(a as u8);
    let b_ncbi = aa_char_to_ncbistdaa(b as u8);
    protein_score(matrix, a_ncbi, b_ncbi) > 0
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

// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showdefline.cpp:69-79
// ```c
// static const char*  kHeader = "Sequences producing significant alignments:";
// static const char*  kBits = "(Bits)";
// static const char*  kValue = "Value";
// ```
fn write_subject_summary_table<W: Write>(
    writer: &mut W,
    subject_order: &[u32],
    subject_hits: &std::collections::HashMap<u32, Vec<&PairwiseHit>>,
    subject_ids: &[Arc<str>],
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
    writeln!(writer, "{:<max_score$}  Value", "(Bits)")?;
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
        writeln!(
            writer,
            "  {:<max_score$}  {:<max_evalue$}",
            format_bitscore_ncbi(best_hit.hit.bit_score),
            format_evalue_ncbi(best_hit.hit.e_value)
        )?;
    }
    writeln!(writer)?;
    Ok(())
}

#[inline]
fn digit_count(value: usize) -> usize {
    value.to_string().len()
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

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:581-613
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
// ...
// sprintf(buffer, "%#8.3g ", lambda);
// ```
fn format_ncbi_ka_value(value: f64) -> String {
    if value == 0.0 {
        return "0.00".to_string();
    }

    let abs = value.abs();
    let exponent = abs.log10().floor() as i32;
    if exponent <= -5 || exponent >= 3 {
        return format!("{value:.2e}");
    }

    let decimals = usize::try_from((2 - exponent).max(0)).unwrap_or(0);
    let mut formatted = format!("{value:.prec$}", prec = decimals);
    if !formatted.contains('.') {
        formatted.push('.');
    }
    formatted
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

fn write_blastp_database_header<W: Write>(
    writer: &mut W,
    database_name: &str,
    database_num_sequences: usize,
    database_total_letters: usize,
) -> io::Result<()> {
    // NCBI reference: c++/src/algo/blast/format/blastfmtutil.cpp:135-137
    // str << dbString << definition_line << endl; x_WrapOutputLine(str.str(),line_len,out);
    write_flatfile_wrapped(
        writer,
        &format!("Database: {}", ensure_trailing_period(database_name)),
        68,
    )?;
    writeln!(
        writer,
        "           {} sequences; {} total letters",
        format_count_with_commas(i64::try_from(database_num_sequences).unwrap_or(0)),
        format_count_with_commas(i64::try_from(database_total_letters).unwrap_or(0))
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)?;
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
    writeln!(writer)?;
    writeln!(writer)?;
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
    writeln!(writer)?;
    writeln!(writer)?;
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
    // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:543-544
    // out << "  Database: "; x_WrapOutputLine(dbinfo->definition, line_length, out);
    write!(writer, "  Database: ")?;
    write_flatfile_wrapped(writer, &ensure_trailing_period(&report.database_name), 68)?;
    writeln!(writer, "    Posted date:  Unknown")?;
    writeln!(
        writer,
        "  Number of letters in database: {}",
        format_count_with_commas(i64::try_from(report.database_total_letters).unwrap_or(0))
    )?;
    writeln!(
        writer,
        "  Number of sequences in database:  {}",
        format_count_with_commas(i64::try_from(report.database_num_sequences).unwrap_or(0))
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer, "Matrix: {}", report.matrix_name)?;
    writeln!(
        writer,
        "Gap Penalties: Existence: {}, Extension: {}",
        report.gap_open, report.gap_extend
    )?;
    if report.word_threshold != 0 {
        writeln!(
            writer,
            "Neighboring words threshold: {}",
            report.word_threshold
        )?;
    }
    if report.window_size != 0 {
        writeln!(writer, "Window for multiple hits: {}", report.window_size)?;
    }
    Ok(())
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
) -> io::Result<()> {
    let mut buffered = io::BufWriter::new(writer);
    let writer = &mut buffered;

    write_blastp_pairwise_intro(writer, &report.version)?;
    write_blastp_database_header(
        writer,
        &report.database_name,
        report.database_num_sequences,
        report.database_total_letters,
    )?;

    let mut hits_by_query: Vec<Vec<&PairwiseHit>> = vec![Vec::new(); queries.len()];
    for hit in hits {
        if let Some(bucket) = hits_by_query.get_mut(hit.hit.q_idx as usize) {
            bucket.push(hit);
        }
    }

    for (q_idx, query) in queries.iter().enumerate() {
        write_blastp_query_header(writer, &query.query_name, query.query_length)?;
        let query_hits = &hits_by_query[q_idx];
        if query_hits.is_empty() {
            write_no_hits_found(writer)?;
            write_blastp_query_footer(
                writer,
                query.ungapped_karlin,
                report.gapped_karlin,
                report.gumbel,
                query.effective_search_space,
            )?;
            continue;
        }

        use std::collections::HashMap;
        let mut subject_hits: HashMap<u32, Vec<&PairwiseHit>> = HashMap::new();
        let mut subject_order: Vec<u32> = Vec::new();
        for hit in query_hits {
            let s_idx = hit.hit.s_idx;
            if !subject_hits.contains_key(&s_idx) {
                subject_order.push(s_idx);
            }
            subject_hits.entry(s_idx).or_default().push(*hit);
        }

        write_subject_summary_table(writer, &subject_order, &subject_hits, subject_ids)?;
        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1520-1589
        // x_DisplayDeflines(aln_set, ...); ... display.DisplaySeqalign(m_Outfile);
        writeln!(writer)?;

        for s_idx in subject_order {
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
            write_subject_header(
                writer,
                subject_id,
                first_hit.subject_title.as_deref(),
                first_hit.subject_length,
            )?;
            for hit in shits {
                write_hsp_info(writer, hit, config)?;
                write_alignment(writer, hit, config)?;
            }
        }

        write_blastp_query_footer(
            writer,
            query.ungapped_karlin,
            report.gapped_karlin,
            report.gumbel,
            query.effective_search_space,
        )?;
    }

    write_blastp_final_footer(writer, report)?;
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

    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/showdefline.cpp:69-79
    // ```c
    // static const char*  kHeader = "Sequences producing significant alignments:";
    // static const char*  kBits = "(Bits)";
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
pub fn write_tblastn_pairwise_report<W: Write>(
    hits: &[PairwiseHit],
    writer: &mut W,
    config: &PairwiseConfig,
    queries: &[BlastpPairwiseQuery],
    query_validity: &[bool],
    query_batch_skipped: &[bool],
    subject_ids: &[Arc<str>],
    report: &BlastpPairwiseReport,
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
    writeln!(writer, "TBLASTN {}", report.version)?;
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
    writeln!(writer)?;
    write_blastp_database_header(
        &mut writer,
        &report.database_name,
        report.database_num_sequences,
        report.database_total_letters,
    )?;

    let mut hits_by_query: Vec<Vec<&PairwiseHit>> = vec![Vec::new(); queries.len()];
    for hit in hits {
        if let Some(bucket) = hits_by_query.get_mut(hit.hit.q_idx as usize) {
            bucket.push(hit);
        }
    }
    for (q_idx, query) in queries.iter().enumerate() {
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
                write_tblastn_unsearched_query_footer(&mut writer)?;
            } else if !query_validity[q_idx] {
                writeln!(writer)?;
                writeln!(writer)?;
                writeln!(writer)?;
                writeln!(writer, "Effective search space used: 0")?;
                writeln!(writer)?;
                writeln!(writer)?;
            } else {
                write_blastp_query_footer(
                    &mut writer,
                    query.ungapped_karlin,
                    report.gapped_karlin,
                    report.gumbel,
                    query.effective_search_space,
                )?;
            }
            continue;
        }
        let mut subject_hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
            std::collections::HashMap::new();
        let mut subject_order = Vec::new();
        for hit in query_hits {
            if !subject_hits.contains_key(&hit.hit.s_idx) {
                subject_order.push(hit.hit.s_idx);
            }
            subject_hits.entry(hit.hit.s_idx).or_default().push(*hit);
        }
        write_subject_summary_table(&mut writer, &subject_order, &subject_hits, subject_ids)?;
        writeln!(writer)?;
        for s_idx in subject_order {
            let shits = &subject_hits[&s_idx];
            let first = shits[0];
            write_subject_header(
                &mut writer,
                &subject_ids[s_idx as usize],
                first.subject_title.as_deref(),
                first.subject_length,
            )?;
            for hit in shits {
                write_tblastn_hsp_info(&mut writer, hit)?;
                write_tblastn_alignment(&mut writer, hit, config)?;
                // NCBI c++/src/objtools/align_format/showalign.cpp:3650-3668:
                // display leaves an extra blank line after each translated HSP.
                writeln!(writer)?;
            }
        }
        write_blastp_query_footer(
            &mut writer,
            query.ungapped_karlin,
            report.gapped_karlin,
            report.gumbel,
            query.effective_search_space,
        )?;
    }
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
    let percent = |n: usize| -> usize {
        if n == length {
            100
        } else if length == 0 {
            0
        } else {
            ((0.5 + 100.0 * n as f64 / length as f64) as usize).min(99)
        }
    };
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
    if let Some(frame) = hit.subject_frame {
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
    let max_pos = hit
        .hit
        .q_start
        .max(hit.hit.q_end)
        .max(hit.hit.s_start)
        .max(hit.hit.s_end);
    let width = digit_count(max_pos);
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
        write!(writer, "Query  {}", q_pos)?;
        write_spaces(writer, width - digit_count(q_pos) + 2)?;
        writer.write_all(qpart)?;
        write!(writer, "  {}\n", q_end)?;
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
        write!(writer, "Sbjct  {}", s_pos)?;
        write_spaces(writer, width - digit_count(s_pos.unsigned_abs()) + 2)?;
        writer.write_all(spart)?;
        write!(writer, "  {}\n", s_end)?;
        writeln!(writer)?;
        q_pos = q_end + 1;
        s_pos = s_end + direction;
    }
    Ok(())
}
