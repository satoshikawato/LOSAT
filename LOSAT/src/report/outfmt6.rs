use crate::common::Hit;
use std::io::{self, Write};

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
/// NCBI format: uses scientific notation with specific formatting
/// BLEMIR format: uses Rust's default scientific notation
pub fn format_evalue(e_value: f64, ncbi_compat: bool) -> String {
    if ncbi_compat {
        // NCBI-style formatting
        if e_value == 0.0 {
            "0.0".to_string()
        } else if e_value < 1e-180 {
            "0.0".to_string()
        } else if e_value < 1e-99 {
            format!("{:.0e}", e_value)
        } else if e_value < 0.01 {
            format!("{:.2e}", e_value)
        } else if e_value < 1.0 {
            format!("{:.2}", e_value)
        } else if e_value < 10.0 {
            format!("{:.1}", e_value)
        } else {
            format!("{:.0}", e_value)
        }
    } else {
        // BLEMIR default formatting
        if e_value == 0.0 {
            "0.0".to_string()
        } else if e_value < 0.001 {
            format!("{:.2e}", e_value)
        } else {
            format!("{:.6}", e_value)
        }
    }
}

/// Format a single hit as outfmt 6 line
pub fn format_hit(hit: &Hit, config: &OutputConfig) -> String {
    let delim = config.delimiter;
    let identity_fmt = format!("{:.prec$}", hit.identity, prec = config.identity_decimals);
    let bit_score_fmt = format!("{:.prec$}", hit.bit_score, prec = config.bit_score_decimals);
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
        Hit {
            query_id: "query1".to_string(),
            subject_id: "subject1".to_string(),
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
        }
    }

    #[test]
    fn test_format_evalue_ncbi() {
        assert_eq!(format_evalue(0.0, true), "0.0");
        assert_eq!(format_evalue(1e-100, true), "1e-100");
        assert_eq!(format_evalue(1e-50, true), "1.00e-50");
        assert_eq!(format_evalue(0.005, true), "5.00e-3");
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
}
