//! Pairwise alignment output (outfmt 0)
//!
//! This module implements the traditional BLAST pairwise alignment output format.
//!
//! Reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp

use crate::common::Hit;
use super::outfmt6::{format_evalue_ncbi, format_bitscore_ncbi, ReportContext};
use std::io::{self, Write};

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
}

impl Default for PairwiseConfig {
    fn default() -> Self {
        Self {
            line_length: DEFAULT_LINE_LENGTH,
            show_gi: false,
            show_frame: true,
            program: "tblastx".to_string(),
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
}

impl From<Hit> for PairwiseHit {
    fn from(hit: Hit) -> Self {
        Self {
            hit,
            query_seq: None,
            subject_seq: None,
            query_frame: None,
            subject_frame: None,
            positives: None,
            gaps: None,
            subject_length: None,
            subject_title: None,
        }
    }
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
    // Subject defline
    if let Some(title) = subject_title {
        writeln!(writer, ">{} {}", subject_id, title)?;
    } else {
        writeln!(writer, ">{}", subject_id)?;
    }
    
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
    writeln!(
        writer,
        " Score = {} bits ({}),  Expect = {}",
        bit_score_str, h.raw_score, evalue_str
    )?;
    
    // Identity/Positives/Gaps line
    let align_len = h.length;
    let num_ident = ((h.identity / 100.0) * align_len as f64).round() as usize;
    let ident_pct = h.identity;
    
    let positives = hit.positives.unwrap_or(num_ident);
    let pos_pct = if align_len > 0 {
        (positives as f64 / align_len as f64) * 100.0
    } else {
        0.0
    };
    
    let gaps = hit.gaps.unwrap_or(0);
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
            " Identities = {}/{} ({:.0}%), Positives = {}/{} ({:.0}%), Gaps = {}/{} ({:.0}%)",
            num_ident, align_len, ident_pct,
            positives, align_len, pos_pct,
            gaps, align_len, gap_pct
        )?;
    } else {
        writeln!(
            writer,
            " Identities = {}/{} ({:.0}%), Gaps = {}/{} ({:.0}%)",
            num_ident, align_len, ident_pct,
            gaps, align_len, gap_pct
        )?;
    }
    
    // Frame line (for translated searches)
    if config.show_frame {
        if let (Some(qf), Some(sf)) = (hit.query_frame, hit.subject_frame) {
            let qf_str = if qf > 0 { format!("+{}", qf) } else { format!("{}", qf) };
            let sf_str = if sf > 0 { format!("+{}", sf) } else { format!("{}", sf) };
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
        writeln!(writer, "Query  {}  [... {} aa ...]  {}", h.q_start, h.length, h.q_end)?;
        writeln!(writer)?;
        writeln!(writer, "Sbjct  {}  [... {} aa ...]  {}", h.s_start, h.length, h.s_end)?;
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
    
    // Calculate position width for formatting
    let max_pos = hit.q_end.max(hit.s_end);
    let pos_width = format!("{}", max_pos).len().max(4);
    
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
                    *q  // Identity: show the character
                } else if is_positive_match(*q, *s) {
                    '+'  // Positive: show +
                } else {
                    ' '  // Mismatch: show space
                }
            })
            .collect();
        
        // Count non-gap characters in this chunk
        let q_non_gap = chunk_q.chars().filter(|&c| c != '-').count();
        let s_non_gap = chunk_s.chars().filter(|&c| c != '-').count();
        
        let q_end_pos = q_pos + q_non_gap.saturating_sub(1);
        let s_end_pos = s_pos + s_non_gap.saturating_sub(1);
        
        // Query line
        writeln!(
            writer,
            "Query  {:>width$}  {}  {}",
            q_pos, chunk_q, q_end_pos,
            width = pos_width
        )?;
        
        // Middle line (identity markers)
        writeln!(
            writer,
            "       {:>width$}  {}",
            "", middle,
            width = pos_width
        )?;
        
        // Subject line
        writeln!(
            writer,
            "Sbjct  {:>width$}  {}  {}",
            s_pos, chunk_s, s_end_pos,
            width = pos_width
        )?;
        
        writeln!(writer)?;
        
        q_pos = q_end_pos + 1;
        s_pos = s_end_pos + 1;
        offset = end;
    }
    
    Ok(())
}

/// Check if two amino acids are a positive match (similar)
/// 
/// This is a simplified check - a full implementation would use the BLOSUM62 matrix
fn is_positive_match(a: char, b: char) -> bool {
    if a == b {
        return true;
    }
    
    // Simplified positive match groups based on BLOSUM62
    let groups = [
        "ILMV",      // Hydrophobic
        "FYW",       // Aromatic
        "KRH",       // Basic (positive)
        "DE",        // Acidic (negative)
        "STNQ",      // Polar
        "AG",        // Small
    ];
    
    for group in &groups {
        if group.contains(a) && group.contains(b) {
            return true;
        }
    }
    
    false
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
    context: &ReportContext,
) -> io::Result<()> {
    // Write program header
    let version = context.version.as_deref().unwrap_or("0.1.0");
    writeln!(writer, "{} {}", context.program.to_uppercase(), version)?;
    writeln!(writer)?;
    
    // Write query info
    if let Some(ref query) = context.query_name {
        writeln!(writer, "Query= {}", query)?;
        writeln!(writer)?;
    }
    
    // Write database info
    if let Some(ref db) = context.subject_name {
        writeln!(writer, "Database: {}", db)?;
        writeln!(writer)?;
    }
    
    if hits.is_empty() {
        writeln!(writer, " ***** No hits found *****")?;
        return Ok(());
    }
    
    // Group hits by subject
    use std::collections::HashMap;
    let mut subject_hits: HashMap<&str, Vec<&PairwiseHit>> = HashMap::new();
    let mut subject_order: Vec<&str> = Vec::new();
    
    for hit in hits {
        let sid = hit.hit.subject_id.as_str();
        if !subject_hits.contains_key(sid) {
            subject_order.push(sid);
        }
        subject_hits.entry(sid).or_default().push(hit);
    }
    
    // Write each subject's hits
    for subject_id in subject_order {
        let shits = subject_hits.get(subject_id).unwrap();
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
    
    Ok(())
}

/// Write hits in pairwise format (simplified version for Hit without sequences)
///
/// This is a convenience function when only Hit structures are available
/// without the extended sequence data.
pub fn write_pairwise_simple<W: Write>(
    hits: &[Hit],
    writer: &mut W,
    config: &PairwiseConfig,
    context: &ReportContext,
) -> io::Result<()> {
    let pairwise_hits: Vec<PairwiseHit> = hits.iter().cloned().map(PairwiseHit::from).collect();
    write_pairwise(&pairwise_hits, writer, config, context)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_hit() -> Hit {
        Hit {
            query_id: "query1".to_string(),
            subject_id: "subject1".to_string(),
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
            q_idx: 0,
            s_idx: 0,
            raw_score: 200,
        }
    }

    #[test]
    fn test_write_subject_header() {
        let mut output = Vec::new();
        write_subject_header(&mut output, "seq1", Some("Test sequence"), Some(500)).unwrap();
        let output_str = String::from_utf8(output).unwrap();
        
        assert!(output_str.contains(">seq1 Test sequence"));
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
            subject_name: Some("test_db".to_string()),
            version: Some("0.1.0".to_string()),
        };
        
        let mut output = Vec::new();
        write_pairwise_simple(&hits, &mut output, &config, &context).unwrap();
        let output_str = String::from_utf8(output).unwrap();
        
        assert!(output_str.contains("TBLASTX"));
        assert!(output_str.contains("Query= test_query"));
        assert!(output_str.contains(">subject1"));
    }
}

