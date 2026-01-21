use crate::common::Hit;
use std::io::{self, Write};

// =============================================================================
// NCBI Output Format Types
// =============================================================================

/// Output format enum matching NCBI BLAST -outfmt options
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    /// 0 = Pairwise alignment view (traditional BLAST output)
    Pairwise = 0,
    /// 6 = Tabular output (tab-separated values)
    Tabular = 6,
    /// 7 = Tabular output with comment lines (headers)
    TabularWithComments = 7,
}

impl Default for OutputFormat {
    fn default() -> Self {
        OutputFormat::Tabular
    }
}

impl OutputFormat {
    /// Parse output format from string (e.g., "0", "6", "7")
    /// Returns the format and any custom field specifications
    pub fn parse(s: &str) -> Result<(Self, Option<String>), String> {
        let s = s.trim();
        
        // Check if there are custom field specifications
        let (fmt_str, fields) = if let Some(idx) = s.find(char::is_whitespace) {
            let (f, rest) = s.split_at(idx);
            (f.trim(), Some(rest.trim().to_string()))
        } else {
            (s, None)
        };
        
        let format = match fmt_str {
            "0" => OutputFormat::Pairwise,
            "6" => OutputFormat::Tabular,
            "7" => OutputFormat::TabularWithComments,
            _ => return Err(format!("Unsupported output format: {}. Supported: 0, 6, 7", fmt_str)),
        };
        
        Ok((format, fields))
    }
}

/// Output format configuration
#[derive(Debug, Clone)]
pub struct OutputConfig {
    /// Whether to include header line
    pub include_header: bool,
    /// Whether to use NCBI-compatible E-value formatting
    pub ncbi_evalue_format: bool,
    /// Delimiter between fields (default: tab)
    pub delimiter: char,
    /// Number of decimal places for identity
    pub identity_decimals: usize,
    /// Number of decimal places for bit score
    pub bit_score_decimals: usize,
}

impl Default for OutputConfig {
    fn default() -> Self {
        Self {
            include_header: false,
            ncbi_evalue_format: false,
            delimiter: '\t',
            identity_decimals: 3,
            bit_score_decimals: 1,
        }
    }
}

impl OutputConfig {
    /// Configuration for NCBI BLAST compatible output
    pub fn ncbi_compat() -> Self {
        Self {
            include_header: false,
            ncbi_evalue_format: true,
            delimiter: '\t',
            identity_decimals: 3,
            bit_score_decimals: 1,
        }
    }

    /// Configuration with header line
    pub fn with_header() -> Self {
        Self {
            include_header: true,
            ..Default::default()
        }
    }
}

/// Context for report generation
#[derive(Debug, Clone)]
pub struct ReportContext {
    /// Query file name or description
    pub query_name: Option<String>,
    /// Subject/database file name or description
    pub subject_name: Option<String>,
    /// Program name (blastn, tblastx, etc.)
    pub program: String,
    /// Version string
    pub version: Option<String>,
}

impl Default for ReportContext {
    fn default() -> Self {
        Self {
            query_name: None,
            subject_name: None,
            program: "blast".to_string(),
            version: None,
        }
    }
}

/// Format E-value for output
///
/// NCBI format: uses scientific notation with specific formatting rules
/// from align_format_util.cpp:GetScoreString() and tabular.cpp:SetScores()
///
/// Reference: ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:940-1011
/// Reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1328-1338
pub fn format_evalue(e_value: f64, ncbi_compat: bool) -> String {
    if ncbi_compat {
        format_evalue_ncbi(e_value)
    } else {
        // LOSAT default formatting
        if e_value == 0.0 {
            "0.0".to_string()
        } else if e_value < 0.001 {
            format!("{:.2e}", e_value)
        } else {
            format!("{:.6}", e_value)
        }
    }
}

/// Format E-value exactly as NCBI BLAST does
///
/// This implements the exact formatting logic from:
/// 1. align_format_util.cpp:GetScoreString() - base formatting
/// 2. tabular.cpp:SetScores() - override for tabular output in [1e-180, 0.0009)
///
/// NCBI GetScoreString() (align_format_util.cpp:965-984):
/// ```c
/// if (evalue < 1.0e-180) {
///     snprintf(evalue_buf, sizeof(evalue_buf), "0.0");
/// } else if (evalue < 1.0e-99) {
///     snprintf(evalue_buf, sizeof(evalue_buf), "%2.0le", evalue);
/// } else if (evalue < 0.0009) {
///     snprintf(evalue_buf, sizeof(evalue_buf), "%3.0le", evalue);
/// } else if (evalue < 0.1) {
///     snprintf(evalue_buf, sizeof(evalue_buf), "%4.3lf", evalue);
/// } else if (evalue < 1.0) {
///     snprintf(evalue_buf, sizeof(evalue_buf), "%3.2lf", evalue);
/// } else if (evalue < 10.0) {
///     snprintf(evalue_buf, sizeof(evalue_buf), "%2.1lf", evalue);
/// } else {
///     snprintf(evalue_buf, sizeof(evalue_buf), "%2.0lf", evalue);
/// }
/// ```
///
/// tabular.cpp:SetScores() override (tabular.cpp:1335-1337):
/// ```c
/// if ((evalue >= 1.0e-180) && (evalue < 0.0009)){
///     m_Evalue = NStr::DoubleToString(evalue, 2, NStr::fDoubleScientific);
/// }
/// ```
pub fn format_evalue_ncbi(e_value: f64) -> String {
    if e_value == 0.0 || e_value < 1.0e-180 {
        // NCBI: snprintf(evalue_buf, sizeof(evalue_buf), "0.0");
        "0.0".to_string()
    } else if e_value < 1.0e-99 {
        // NCBI: "%2.0le" -> e.g., "1e-100"
        // But tabular.cpp overrides for [1e-180, 0.0009) to use %.2e
        // For < 1e-99, use "%2.0le" style (no decimal in mantissa)
        format!("{:.0e}", e_value)
    } else if e_value < 0.0009 {
        // tabular.cpp overrides: NStr::DoubleToString(evalue, 2, fDoubleScientific)
        // -> "%.2e" format, e.g., "1.23e-50"
        format!("{:.2e}", e_value)
    } else if e_value < 0.1 {
        // NCBI: "%4.3lf" -> e.g., "0.005"
        format!("{:.3}", e_value)
    } else if e_value < 1.0 {
        // NCBI: "%3.2lf" -> e.g., "0.50"
        format!("{:.2}", e_value)
    } else if e_value < 10.0 {
        // NCBI: "%2.1lf" -> e.g., "5.5"
        format!("{:.1}", e_value)
    } else {
        // NCBI: "%2.0lf" -> e.g., "100"
        format!("{:.0}", e_value)
    }
}

/// Format bit score exactly as NCBI BLAST does
///
/// Reference: ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:986-994
/// ```c
/// if (bit_score > 99999){
///     snprintf(bit_score_buf, sizeof(bit_score_buf), "%5.3le", bit_score);
/// } else if (bit_score > 99.9){
///     snprintf(bit_score_buf, sizeof(bit_score_buf), "%3.0ld", (long)bit_score);
/// } else {
///     snprintf(bit_score_buf, sizeof(bit_score_buf), kBitScoreFormat.c_str(), bit_score);
/// }
/// // where kBitScoreFormat = "%4.1lf"
/// ```
pub fn format_bitscore_ncbi(bit_score: f64) -> String {
    if bit_score > 99999.0 {
        // NCBI: "%5.3le" -> scientific notation
        format!("{:.3e}", bit_score)
    } else if bit_score > 99.9 {
        // NCBI: "%3.0ld" -> integer (no decimal)
        format!("{:.0}", bit_score)
    } else {
        // NCBI: "%4.1lf" -> one decimal place
        format!("{:.1}", bit_score)
    }
}

/// Format a single hit as outfmt 6 line
pub fn format_hit(hit: &Hit, config: &OutputConfig) -> String {
    let delim = config.delimiter;
    let identity_fmt = format!("{:.prec$}", hit.identity, prec = config.identity_decimals);
    
    // Use NCBI-compatible formatting when enabled
    let bit_score_fmt = if config.ncbi_evalue_format {
        format_bitscore_ncbi(hit.bit_score)
    } else {
        format!("{:.prec$}", hit.bit_score, prec = config.bit_score_decimals)
    };
    let evalue_fmt = format_evalue(hit.e_value, config.ncbi_evalue_format);

    format!(
        "{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}",
        hit.query_id,
        delim,
        hit.subject_id,
        delim,
        identity_fmt,
        delim,
        hit.length,
        delim,
        hit.mismatch,
        delim,
        hit.gapopen,
        delim,
        hit.q_start,
        delim,
        hit.q_end,
        delim,
        hit.s_start,
        delim,
        hit.s_end,
        delim,
        evalue_fmt,
        delim,
        bit_score_fmt
    )
}

/// Write a single hit as outfmt 6 line without building an intermediate String.
// NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
// ```c
// void CBlastTabularInfo::Print()
// {
//     ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
//         if (iter != m_FieldsToShow.begin())
//             m_Ostream << m_FieldDelimiter;
//         x_PrintField(*iter);
//     }
//     m_Ostream << "\n";
// }
// ```
pub fn write_hit_fields<W: Write>(
    writer: &mut W,
    query_id: &str,
    subject_id: &str,
    identity: f64,
    length: usize,
    mismatch: usize,
    gapopen: usize,
    q_start: usize,
    q_end: usize,
    s_start: usize,
    s_end: usize,
    e_value: f64,
    bit_score: f64,
    config: &OutputConfig,
) -> io::Result<()> {
    let delim = config.delimiter;
    let bit_score_fmt = if config.ncbi_evalue_format {
        format_bitscore_ncbi(bit_score)
    } else {
        format!("{:.prec$}", bit_score, prec = config.bit_score_decimals)
    };
    let evalue_fmt = format_evalue(e_value, config.ncbi_evalue_format);

    write!(writer, "{}{}{}", query_id, delim, subject_id)?;
    write!(writer, "{}{:.prec$}", delim, identity, prec = config.identity_decimals)?;
    write!(writer, "{}{}", delim, length)?;
    write!(writer, "{}{}", delim, mismatch)?;
    write!(writer, "{}{}", delim, gapopen)?;
    write!(writer, "{}{}", delim, q_start)?;
    write!(writer, "{}{}", delim, q_end)?;
    write!(writer, "{}{}", delim, s_start)?;
    write!(writer, "{}{}", delim, s_end)?;
    write!(writer, "{}{}", delim, evalue_fmt)?;
    write!(writer, "{}{}", delim, bit_score_fmt)?;
    writeln!(writer)
}

/// Write hits to a writer in outfmt 6 format
pub fn write_outfmt6<W: Write>(
    hits: &[Hit],
    writer: &mut W,
    config: &OutputConfig,
) -> io::Result<()> {
    // Write header if configured
    if config.include_header {
        let delim = config.delimiter;
        writeln!(
            writer,
            "qseqid{}sseqid{}pident{}length{}mismatch{}gapopen{}qstart{}qend{}sstart{}send{}evalue{}bitscore",
            delim, delim, delim, delim, delim, delim, delim, delim, delim, delim, delim
        )?;
    }

    // Write each hit
    for hit in hits {
        writeln!(writer, "{}", format_hit(hit, config))?;
    }

    Ok(())
}

// =============================================================================
// outfmt 7: Tabular with comment lines (NCBI BLAST compatible)
// Reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:PrintHeader()
// =============================================================================

/// Write outfmt 7 header (comment lines) to a writer
///
/// NCBI format (from tabular.cpp:1266-1284):
/// ```text
/// # TBLASTX 2.17.0+
/// # Query: query_name
/// # Database: database_name
/// # Fields: qaccver, saccver, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
/// # X hits found
/// ```
pub fn write_outfmt7_header<W: Write>(
    writer: &mut W,
    context: &ReportContext,
    num_hits: usize,
) -> io::Result<()> {
    // Program version line
    // NCBI: "# TBLASTX 2.17.0+"
    let version_str = context.version.as_deref().unwrap_or("0.1.0");
    writeln!(writer, "# {} {}", context.program.to_uppercase(), version_str)?;
    
    // Query line
    // NCBI: "# Query: query_name"
    if let Some(ref query) = context.query_name {
        writeln!(writer, "# Query: {}", query)?;
    }
    
    // Database line
    // NCBI: "# Database: database_name"
    if let Some(ref db) = context.subject_name {
        writeln!(writer, "# Database: {}", db)?;
    }
    
    // Fields line (only if there are hits)
    // NCBI: "# Fields: qaccver, saccver, pident, ..."
    if num_hits > 0 {
        writeln!(
            writer,
            "# Fields: qaccver, saccver, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore"
        )?;
    }
    
    // Hits count line
    // NCBI: "# X hits found"
    writeln!(writer, "# {} hits found", num_hits)?;
    
    Ok(())
}

/// Write hits to a writer in outfmt 7 format (tabular with comment headers)
///
/// Reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp
pub fn write_outfmt7<W: Write>(
    hits: &[Hit],
    writer: &mut W,
    config: &OutputConfig,
    context: &ReportContext,
) -> io::Result<()> {
    // Write header comments
    write_outfmt7_header(writer, context, hits.len())?;
    
    // Write each hit (same format as outfmt 6)
    for hit in hits {
        writeln!(writer, "{}", format_hit(hit, config))?;
    }
    
    Ok(())
}

/// Write hits grouped by query in outfmt 7 format
/// 
/// Each query gets its own header block followed by its hits.
/// This matches NCBI's output where headers are repeated per query.
pub fn write_outfmt7_grouped<W: Write>(
    hits: &[Hit],
    writer: &mut W,
    config: &OutputConfig,
    context: &ReportContext,
) -> io::Result<()> {
    use std::collections::HashMap;
    
    // Group hits by query
    let mut query_hits: HashMap<&str, Vec<&Hit>> = HashMap::new();
    let mut query_order: Vec<&str> = Vec::new();
    
    for hit in hits {
        // NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1264-1283
        // ```c
        // void CBlastTabularInfo::PrintHeader(...)
        // {
        //     x_PrintQueryAndDbNames(program_version, bioseq, dbname, rid, iteration, subj_bioseq);
        //     if (align_set) {
        //         int num_hits = align_set->Get().size();
        //         if (num_hits != 0) {
        //             PrintFieldNames(is_csv);
        //         }
        //         m_Ostream << "# " << num_hits << " hits found" << "\n";
        //     }
        // }
        // ```
        let qid = hit.query_id.as_ref();
        if !query_hits.contains_key(qid) {
            query_order.push(qid);
        }
        query_hits.entry(qid).or_default().push(hit);
    }
    
    // Write each query's header and hits
    for query_id in query_order {
        let query_context = ReportContext {
            query_name: Some(query_id.to_string()),
            ..context.clone()
        };
        
        let qhits = query_hits.get(query_id).unwrap();
        write_outfmt7_header(writer, &query_context, qhits.len())?;
        
        for hit in qhits {
            writeln!(writer, "{}", format_hit(hit, config))?;
        }
    }
    
    Ok(())
}

/// Write hits to a file in outfmt 6 format
pub fn write_to_file(hits: &[Hit], path: &str, config: &OutputConfig) -> io::Result<()> {
    let file = std::fs::File::create(path)?;
    let mut writer = std::io::BufWriter::new(file);
    write_outfmt6(hits, &mut writer, config)
}

/// Write hits to stdout in outfmt 6 format
pub fn write_to_stdout(hits: &[Hit], config: &OutputConfig) -> io::Result<()> {
    let stdout = io::stdout();
    let mut writer = stdout.lock();
    write_outfmt6(hits, &mut writer, config)
}

/// Generate a summary report of the search results
pub fn generate_summary(hits: &[Hit], context: &ReportContext) -> String {
    let mut summary = String::new();

    // Header
    summary.push_str(&format!(
        "# {} Search Results\n",
        context.program.to_uppercase()
    ));

    if let Some(ref version) = context.version {
        summary.push_str(&format!("# Version: {}\n", version));
    }

    if let Some(ref query) = context.query_name {
        summary.push_str(&format!("# Query: {}\n", query));
    }

    if let Some(ref subject) = context.subject_name {
        summary.push_str(&format!("# Database: {}\n", subject));
    }

    summary.push_str(&format!("# Total hits: {}\n", hits.len()));

    if !hits.is_empty() {
        // Statistics
        let total_bit_score: f64 = hits.iter().map(|h| h.bit_score).sum();
        let avg_bit_score = total_bit_score / hits.len() as f64;
        let avg_identity: f64 = hits.iter().map(|h| h.identity).sum::<f64>() / hits.len() as f64;
        let min_evalue = hits.iter().map(|h| h.e_value).fold(f64::INFINITY, f64::min);

        summary.push_str(&format!("# Average bit score: {:.1}\n", avg_bit_score));
        summary.push_str(&format!("# Average identity: {:.1}%\n", avg_identity));
        summary.push_str(&format!("# Best E-value: {:.2e}\n", min_evalue));

        // Unique query-subject pairs
        let mut pairs = std::collections::HashSet::new();
        for hit in hits {
            pairs.insert((&hit.query_id, &hit.subject_id));
        }
        summary.push_str(&format!("# Unique query-subject pairs: {}\n", pairs.len()));
    }

    summary
}

/// Standard BLAST outfmt 6 column names
pub const OUTFMT6_COLUMNS: &[&str] = &[
    "qseqid",   // Query sequence ID
    "sseqid",   // Subject sequence ID
    "pident",   // Percentage of identical matches
    "length",   // Alignment length
    "mismatch", // Number of mismatches
    "gapopen",  // Number of gap openings
    "qstart",   // Start of alignment in query
    "qend",     // End of alignment in query
    "sstart",   // Start of alignment in subject
    "send",     // End of alignment in subject
    "evalue",   // Expect value
    "bitscore", // Bit score
];

/// Get column description
pub fn get_column_description(column: &str) -> Option<&'static str> {
    match column {
        "qseqid" => Some("Query sequence ID"),
        "sseqid" => Some("Subject sequence ID"),
        "pident" => Some("Percentage of identical matches"),
        "length" => Some("Alignment length"),
        "mismatch" => Some("Number of mismatches"),
        "gapopen" => Some("Number of gap openings"),
        "qstart" => Some("Start of alignment in query"),
        "qend" => Some("End of alignment in query"),
        "sstart" => Some("Start of alignment in subject"),
        "send" => Some("End of alignment in subject"),
        "evalue" => Some("Expect value"),
        "bitscore" => Some("Bit score"),
        _ => None,
    }
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
            query_id: "query1".into(),
            subject_id: "subject1".into(),
            identity: 95.123,
            length: 100,
            mismatch: 5,
            gapopen: 0,
            q_start: 1,
            q_end: 100,
            s_start: 1,
            s_end: 100,
            e_value: 1e-50,
            bit_score: 185.5,
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
            gap_info: None,
        }
    }

    // ==========================================================================
    // E-value formatting tests (NCBI parity)
    // Reference: align_format_util.cpp:GetScoreString() + tabular.cpp:SetScores()
    // ==========================================================================
    
    #[test]
    fn test_format_evalue_ncbi_zero() {
        // Values at or below 0.0 should return "0.0"
        assert_eq!(format_evalue_ncbi(0.0), "0.0");
    }

    #[test]
    fn test_format_evalue_ncbi_very_small() {
        // < 1e-180: "0.0"
        assert_eq!(format_evalue_ncbi(1e-200), "0.0");
        assert_eq!(format_evalue_ncbi(1e-181), "0.0");
    }

    #[test]
    fn test_format_evalue_ncbi_small_scientific() {
        // < 1e-99: "%2.0le" style (no decimal in mantissa)
        assert_eq!(format_evalue_ncbi(1e-100), "1e-100");
        assert_eq!(format_evalue_ncbi(1e-150), "1e-150");
        // Note: Rust may format 5e-100 differently
    }

    #[test]
    fn test_format_evalue_ncbi_medium_scientific() {
        // [1e-99, 0.0009): "%.2e" format (tabular.cpp override)
        assert_eq!(format_evalue_ncbi(1e-50), "1.00e-50");
        assert_eq!(format_evalue_ncbi(1.23e-20), "1.23e-20");
        assert_eq!(format_evalue_ncbi(5e-5), "5.00e-5");
    }

    #[test]
    fn test_format_evalue_ncbi_decimal() {
        // [0.0009, 0.1): "%4.3lf" - 3 decimal places
        assert_eq!(format_evalue_ncbi(0.001), "0.001");
        assert_eq!(format_evalue_ncbi(0.005), "0.005");
        assert_eq!(format_evalue_ncbi(0.099), "0.099");
        
        // [0.1, 1.0): "%3.2lf" - 2 decimal places
        assert_eq!(format_evalue_ncbi(0.1), "0.10");
        assert_eq!(format_evalue_ncbi(0.5), "0.50");
        assert_eq!(format_evalue_ncbi(0.99), "0.99");
        
        // [1.0, 10.0): "%2.1lf" - 1 decimal place
        assert_eq!(format_evalue_ncbi(1.0), "1.0");
        assert_eq!(format_evalue_ncbi(5.5), "5.5");
        assert_eq!(format_evalue_ncbi(9.9), "9.9");
    }

    #[test]
    fn test_format_evalue_ncbi_integer() {
        // >= 10.0: "%2.0lf" - no decimal
        assert_eq!(format_evalue_ncbi(10.0), "10");
        assert_eq!(format_evalue_ncbi(100.0), "100");
        assert_eq!(format_evalue_ncbi(1000.0), "1000");
    }

    // ==========================================================================
    // Bit score formatting tests (NCBI parity)
    // Reference: align_format_util.cpp:GetScoreString()
    // ==========================================================================

    #[test]
    fn test_format_bitscore_ncbi_small() {
        // <= 99.9: "%4.1lf" - one decimal place
        assert_eq!(format_bitscore_ncbi(0.0), "0.0");
        assert_eq!(format_bitscore_ncbi(50.5), "50.5");
        assert_eq!(format_bitscore_ncbi(99.9), "99.9");
    }

    #[test]
    fn test_format_bitscore_ncbi_medium() {
        // > 99.9 and <= 99999: "%3.0ld" - integer (no decimal)
        assert_eq!(format_bitscore_ncbi(100.0), "100");
        assert_eq!(format_bitscore_ncbi(185.5), "186"); // Rounds
        assert_eq!(format_bitscore_ncbi(692.0), "692");
        assert_eq!(format_bitscore_ncbi(99999.0), "99999");
    }

    #[test]
    fn test_format_bitscore_ncbi_large() {
        // > 99999: "%5.3le" - scientific notation
        assert_eq!(format_bitscore_ncbi(100000.0), "1.000e5");
        assert_eq!(format_bitscore_ncbi(123456.789), "1.235e5");
    }

    // ==========================================================================
    // Legacy tests
    // ==========================================================================

    #[test]
    fn test_format_evalue_compat() {
        // Test the wrapper function with ncbi_compat flag
        assert_eq!(format_evalue(0.0, true), "0.0");
        assert_eq!(format_evalue(1e-100, true), "1e-100");
        assert_eq!(format_evalue(1e-50, true), "1.00e-50");
        assert_eq!(format_evalue(0.5, true), "0.50");
        assert_eq!(format_evalue(5.5, true), "5.5");
        assert_eq!(format_evalue(100.0, true), "100");
    }

    #[test]
    fn test_format_hit() {
        let hit = make_hit();
        let config = OutputConfig::default();
        let line = format_hit(&hit, &config);

        assert!(line.contains("query1"));
        assert!(line.contains("subject1"));
        assert!(line.contains("95.123"));
        assert!(line.contains("185.5"));
    }

    #[test]
    fn test_format_hit_ncbi_compat() {
        let hit = make_hit();
        let config = OutputConfig::ncbi_compat();
        let line = format_hit(&hit, &config);

        assert!(line.contains("query1"));
        assert!(line.contains("subject1"));
        assert!(line.contains("95.123"));
        // Bit score 185.5 > 99.9 should be formatted as integer "186"
        assert!(line.contains("186"));
        // E-value 1e-50 should be "1.00e-50"
        assert!(line.contains("1.00e-50"));
    }

    #[test]
    fn test_write_with_header() {
        let hits = vec![make_hit()];
        let config = OutputConfig::with_header();

        let mut output = Vec::new();
        write_outfmt6(&hits, &mut output, &config).unwrap();

        let output_str = String::from_utf8(output).unwrap();
        assert!(output_str.starts_with("qseqid"));
        assert!(output_str.contains("query1"));
    }

    #[test]
    fn test_generate_summary() {
        let hits = vec![make_hit()];
        let context = ReportContext {
            program: "blastn".to_string(),
            query_name: Some("test.fasta".to_string()),
            ..Default::default()
        };

        let summary = generate_summary(&hits, &context);
        assert!(summary.contains("BLASTN"));
        assert!(summary.contains("Total hits: 1"));
    }

    // ==========================================================================
    // OutputFormat tests
    // ==========================================================================

    #[test]
    fn test_output_format_parse() {
        assert_eq!(OutputFormat::parse("0").unwrap(), (OutputFormat::Pairwise, None));
        assert_eq!(OutputFormat::parse("6").unwrap(), (OutputFormat::Tabular, None));
        assert_eq!(OutputFormat::parse("7").unwrap(), (OutputFormat::TabularWithComments, None));
        
        // With custom fields
        let (fmt, fields) = OutputFormat::parse("6 qaccver saccver").unwrap();
        assert_eq!(fmt, OutputFormat::Tabular);
        assert_eq!(fields, Some("qaccver saccver".to_string()));
        
        // Invalid format
        assert!(OutputFormat::parse("99").is_err());
    }
}
