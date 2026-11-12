use std::io::{self, Write};

// =============================================================================
// NCBI Output Format Types
// =============================================================================

/// Output format enum matching NCBI BLAST -outfmt options
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutputFormat {
    /// 0 = Pairwise alignment view (traditional BLAST output)
    Pairwise = 0,
    /// 6 = Tabular output (tab-separated values)
    #[default]
    Tabular = 6,
    /// 7 = Tabular output with comment lines (headers)
    TabularWithComments = 7,
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
            _ => {
                return Err(format!(
                    "Unsupported output format: {}. Supported: 0, 6, 7",
                    fmt_str
                ))
            }
        };

        Ok((format, fields))
    }
}

/// Output format configuration
#[derive(Debug, Clone)]
pub struct OutputConfig {
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
            ncbi_evalue_format: true,
            delimiter: '\t',
            identity_decimals: 3,
            bit_score_decimals: 1,
        }
    }
}

/// Context for report generation
#[derive(Debug, Clone)]
pub struct ReportContext {
    /// Subject/database file name or description
    pub subject_name: Option<String>,
    /// Fully formatted database label for pairwise formatting
    pub database_name: Option<String>,
    /// Program name (blastn, tblastx, etc.)
    pub program: String,
}

impl Default for ReportContext {
    fn default() -> Self {
        Self {
            subject_name: None,
            database_name: None,
            program: "blast".to_string(),
        }
    }
}

/// Format scientific notation with NCBI-style exponent padding.
///
/// NCBI reference: ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:965-988
/// ```c
/// if (evalue < 1.0e-99) {
///     snprintf(evalue_buf, sizeof(evalue_buf), "%2.0le", evalue);
/// } else if (evalue < 0.0009) {
///     snprintf(evalue_buf, sizeof(evalue_buf), "%3.0le", evalue);
/// }
/// if (bit_score > 99999){
///     snprintf(bit_score_buf, sizeof(bit_score_buf), "%5.3le", bit_score);
/// }
/// ```
/// NCBI reference: ncbi-blast/c++/src/corelib/ncbistr.cpp:2389-2402
/// ```c
/// switch (flags & fDoubleGeneral) {
///     case fDoubleScientific:
///         format = "%.*e";
///         break;
/// }
/// n = ::snprintf(buffer, kMaxDoubleStringSize, format, (int)precision, value);
/// ```
fn format_scientific_ncbi(value: f64, precision: usize) -> String {
    let mut raw = format!("{:.*e}", precision, value);
    let exp_pos = match raw.find('e').or_else(|| raw.find('E')) {
        Some(pos) => pos,
        None => return raw,
    };

    // Force lowercase exponent marker to match NCBI output.
    if raw.as_bytes()[exp_pos] == b'E' {
        raw.replace_range(exp_pos..=exp_pos, "e");
    }

    let exp_start = exp_pos + 1;
    if exp_start >= raw.len() {
        return raw;
    }

    let has_sign = matches!(raw.as_bytes()[exp_start], b'+' | b'-');
    if !has_sign {
        raw.insert(exp_start, '+');
    }
    // The exponent's digits start after its sign, whether inserted or already there.
    let digits_start = exp_start + 1;

    let digits_len = raw.len().saturating_sub(digits_start);
    if digits_len < 2 {
        raw.insert(digits_start, '0');
    }

    raw
}

/// Write scientific notation with NCBI-style exponent padding.
///
/// NCBI reference: ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:965-988
/// ```c
/// if (evalue < 1.0e-99) {
///     snprintf(evalue_buf, sizeof(evalue_buf), "%2.0le", evalue);
/// } else if (evalue < 0.0009) {
///     snprintf(evalue_buf, sizeof(evalue_buf), "%3.0le", evalue);
/// }
/// if (bit_score > 99999){
///     snprintf(bit_score_buf, sizeof(bit_score_buf), "%5.3le", bit_score);
/// }
/// ```
fn write_scientific_ncbi<W: Write>(writer: &mut W, value: f64, precision: usize) -> io::Result<()> {
    let formatted = format_scientific_ncbi(value, precision);
    writer.write_all(formatted.as_bytes())
}

/// Format E-value exactly as NCBI BLAST does (base GetScoreString logic).
///
/// NCBI reference: ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:965-983
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
pub fn format_evalue_ncbi(e_value: f64) -> String {
    if e_value == 0.0 || e_value < 1.0e-180 {
        // NCBI: snprintf(evalue_buf, sizeof(evalue_buf), "0.0");
        "0.0".to_string()
    } else if e_value < 1.0e-99 {
        // NCBI: "%2.0le" -> e.g., "1e-100"
        format_scientific_ncbi(e_value, 0)
    } else if e_value < 0.0009 {
        // NCBI: "%3.0le" -> e.g., "1e-50"
        format_scientific_ncbi(e_value, 0)
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

/// Format E-value exactly as NCBI BLAST does for tabular output.
///
/// NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1332-1337
/// ```c
/// CAlignFormatUtil::GetScoreString(evalue, bit_score, 0, score, m_Evalue, ...);
/// if ((evalue >= 1.0e-180) && (evalue < 0.0009)){
///     m_Evalue = NStr::DoubleToString(evalue, 2, NStr::fDoubleScientific);
/// }
/// ```
pub fn format_evalue_ncbi_tabular(e_value: f64) -> String {
    if e_value == 0.0 || e_value < 1.0e-180 {
        "0.0".to_string()
    } else if e_value < 0.0009 {
        format_scientific_ncbi(e_value, 2)
    } else {
        format_evalue_ncbi(e_value)
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
        format_scientific_ncbi(bit_score, 3)
    } else if bit_score > 99.9 {
        // NCBI: "%3.0ld" with `(long)bit_score` truncates toward zero before formatting,
        // in a field of 3 (a score in (99.9, 100) is " 99").
        format!("{:3}", bit_score as i64)
    } else {
        // NCBI: "%4.1lf" -> one decimal place
        format!("{:4.1}", bit_score)
    }
}

/// The total bit score of a description table row (`total_bit_string`), which `GetScoreString`
/// formats as the bit score but with `%2.1lf` below 99.9.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:994-1003
/// ```c++
///     if (total_bit_score > 99999){
///         snprintf(total_bit_score_buf, sizeof(total_bit_score_buf), "%5.3le",
///             total_bit_score);
///     } else if (total_bit_score > 99.9){
///         snprintf(total_bit_score_buf, sizeof(total_bit_score_buf), "%3.0ld",
///             (long)total_bit_score);
///     } else {
///         snprintf(total_bit_score_buf, sizeof(total_bit_score_buf), "%2.1lf",
///             total_bit_score);
///     }
/// ```
pub fn format_total_bitscore_ncbi(total_bit_score: f64) -> String {
    if total_bit_score > 99.9 {
        format_bitscore_ncbi(total_bit_score)
    } else {
        format!("{:2.1}", total_bit_score)
    }
}

/// Format percent identity exactly as NCBI tabular output does.
///
/// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/objtools/align_format/tabular.hpp:437-441
/// ```c
/// inline void CBlastTabularInfo::x_PrintPercentIdentical(void)
/// {
///     double perc_ident =
///         (m_AlignLength > 0 ? ((double)m_NumIdent)/m_AlignLength * 100 : 0);
///     m_Ostream << NStr::DoubleToString(perc_ident, 3);
/// }
/// ```
/// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/corelib/ncbistr.cpp:2388-2401
/// ```c
/// case fDoubleFixed: // default
/// default:
///     format = "%.*f";
///     break;
/// ...
/// n = ::snprintf(buffer, kMaxDoubleStringSize, format, (int)precision, value);
/// ```
pub fn format_percent_identity_ncbi(
    num_ident: usize,
    alignment_length: usize,
    precision: usize,
) -> String {
    let perc_ident = percent_identity_ncbi(num_ident, alignment_length);
    format!("{:.prec$}", perc_ident, prec = precision)
}

/// The value that `format_percent_identity_ncbi` prints, not rounded: the expression of
/// `x_PrintPercentIdentical` (c++/include/objtools/align_format/tabular.hpp:437-442, above).
pub fn percent_identity_ncbi(num_ident: usize, alignment_length: usize) -> f64 {
    if alignment_length > 0 {
        (num_ident as f64 / alignment_length as f64) * 100.0
    } else {
        0.0
    }
}

/// Write E-value exactly as NCBI BLAST does (base GetScoreString logic).
///
/// NCBI reference: ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:965-983
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
fn write_evalue_ncbi_base<W: Write>(writer: &mut W, e_value: f64) -> io::Result<()> {
    if e_value == 0.0 || e_value < 1.0e-180 {
        write!(writer, "0.0")
    } else if e_value < 0.0009 {
        // NCBI's `%2.0le` (below 1.0e-99) and `%3.0le` (below 0.0009) print the same bytes: the
        // `%.0le` text (`de-XX`, at least 5 characters) is wider than either width.
        write_scientific_ncbi(writer, e_value, 0)
    } else if e_value < 0.1 {
        write!(writer, "{:.3}", e_value)
    } else if e_value < 1.0 {
        write!(writer, "{:.2}", e_value)
    } else if e_value < 10.0 {
        write!(writer, "{:.1}", e_value)
    } else {
        write!(writer, "{:.0}", e_value)
    }
}

/// Write E-value exactly as NCBI BLAST does for tabular output.
///
/// NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1332-1337
/// ```c
/// CAlignFormatUtil::GetScoreString(evalue, bit_score, 0, score, m_Evalue, ...);
/// if ((evalue >= 1.0e-180) && (evalue < 0.0009)){
///     m_Evalue = NStr::DoubleToString(evalue, 2, NStr::fDoubleScientific);
/// }
/// ```
fn write_evalue_ncbi_tabular<W: Write>(writer: &mut W, e_value: f64) -> io::Result<()> {
    if e_value == 0.0 || e_value < 1.0e-180 {
        write!(writer, "0.0")
    } else if e_value < 0.0009 {
        write_scientific_ncbi(writer, e_value, 2)
    } else {
        write_evalue_ncbi_base(writer, e_value)
    }
}

/// Write bit score exactly as NCBI BLAST does.
///
/// NCBI reference: ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:986-993
/// ```c
/// if (bit_score > 99999){
///     snprintf(bit_score_buf, sizeof(bit_score_buf), "%5.3le", bit_score);
/// } else if (bit_score > 99.9){
///     snprintf(bit_score_buf, sizeof(bit_score_buf), "%3.0ld",
///         (long)bit_score);
/// } else {
///     snprintf(bit_score_buf, sizeof(bit_score_buf), kBitScoreFormat.c_str(),
///         bit_score);
/// }
/// ```
fn write_bitscore_ncbi<W: Write>(writer: &mut W, bit_score: f64) -> io::Result<()> {
    if bit_score > 99999.0 {
        write_scientific_ncbi(writer, bit_score, 3)
    } else if bit_score > 99.9 {
        write!(writer, "{:3}", bit_score as i64)
    } else {
        write!(writer, "{:4.1}", bit_score)
    }
}

/// The twelve fields of a default tabular row (`qaccver saccver pident length mismatch gapopen
/// qstart qend sstart send evalue bitscore`) with the query and subject IDs as bytes, written as
/// they are (NCBI writes the IDs that it takes from the titles without converting them).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:214-217
/// ```c++
/// void CBlastTabularInfo::x_PrintQueryAccessionVersion()
/// {
///     m_Ostream << s_GetSeqIdListString(m_QueryId, eAccVersion);
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:252-255
/// ```c++
/// void CBlastTabularInfo::x_PrintSubjectAccessionVersion(void)
/// {
///     m_Ostream << s_GetSeqIdListString(m_SubjectId, eAccVersion);
/// }
/// ```
#[allow(clippy::too_many_arguments)]
pub fn write_hit_fields_bytes<W: Write>(
    writer: &mut W,
    query_id: &[u8],
    subject_id: &[u8],
    identity: f64,
    num_ident: usize,
    length: usize,
    mismatch: i64,
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
    writer.write_all(query_id)?;
    write!(writer, "{}", delim)?;
    writer.write_all(subject_id)?;
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/objtools/align_format/tabular.hpp:437-441
    // ```c
    // double perc_ident =
    //     (m_AlignLength > 0 ? ((double)m_NumIdent)/m_AlignLength * 100 : 0);
    // m_Ostream << NStr::DoubleToString(perc_ident, 3);
    // ```
    let identity_fmt = if config.ncbi_evalue_format {
        format_percent_identity_ncbi(num_ident, length, config.identity_decimals)
    } else {
        format!("{:.prec$}", identity, prec = config.identity_decimals)
    };
    write!(writer, "{}{}", delim, identity_fmt)?;
    write!(writer, "{}{}", delim, length)?;
    write!(writer, "{}{}", delim, mismatch)?;
    write!(writer, "{}{}", delim, gapopen)?;
    write!(writer, "{}{}", delim, q_start)?;
    write!(writer, "{}{}", delim, q_end)?;
    write!(writer, "{}{}", delim, s_start)?;
    write!(writer, "{}{}", delim, s_end)?;
    write!(writer, "{}", delim)?;
    if config.ncbi_evalue_format {
        write_evalue_ncbi_tabular(writer, e_value)?;
    } else if e_value == 0.0 {
        write!(writer, "0.0")?;
    } else if e_value < 0.001 {
        write!(writer, "{:.2e}", e_value)?;
    } else {
        write!(writer, "{:.6}", e_value)?;
    }
    write!(writer, "{}", delim)?;
    if config.ncbi_evalue_format {
        write_bitscore_ncbi(writer, bit_score)?;
    } else {
        write!(
            writer,
            "{:.prec$}",
            bit_score,
            prec = config.bit_score_decimals
        )?;
    }
    writeln!(writer)
}

/// `NStr::TruncateSpaces(text)` (`eTrunc_Both`): `text` without the white space at its
/// ends (`isspace` of an `unsigned char` in the C locale: 0x09-0x0D and the space).
///
/// NCBI reference: c++/src/corelib/ncbistr.cpp:3145-3175
/// ```c
/// template <class TStr>
/// TStr s_TruncateSpaces(const TStr& str, NStr::ETrunc where,
///                       const TStr& empty_str)
/// {
///     SIZE_TYPE length = str.length();
///     if (length == 0) {
///         return empty_str;
///     }
///     SIZE_TYPE beg = 0;
///     if (where == NStr::eTrunc_Begin  ||  where == NStr::eTrunc_Both) {
///         _ASSERT(beg < length);
///         while ( isspace((unsigned char) str[beg]) ) {
///             if (++beg == length) {
///                 return empty_str;
///             }
///         }
///     }
///     SIZE_TYPE end = length;
///     if ( where == NStr::eTrunc_End  ||  where == NStr::eTrunc_Both ) {
///         _ASSERT(beg < end);
///         while (isspace((unsigned char) str[--end])) {
///             if (beg == end) {
///                 return empty_str;
///             }
///         }
///         _ASSERT(beg <= end  &&  !isspace((unsigned char) str[end]));
///         ++end;
///     }
/// ```
pub fn truncate_spaces(text: &[u8]) -> &[u8] {
    let is_space = |byte: &u8| matches!(byte, 0x09..=0x0D | b' ');
    let begin = text
        .iter()
        .position(|byte| !is_space(byte))
        .unwrap_or(text.len());
    let end = text
        .iter()
        .rposition(|byte| !is_space(byte))
        .map_or(begin, |last| last + 1);
    &text[begin..end]
}

/// The `# Query: ` line of the outfmt 7 header of a query with the title `query_title`
/// (the bytes as they are; the ID of a local query is not shown).
///
/// NCBI reference: c++/src/objtools/align_format/tabular.cpp:1305-1308
/// ```c
///     CAlignFormatUtil::AcknowledgeBlastQuery(bioseq, kLineLength, m_Ostream,
///                                             m_ParseLocalIds, kHtmlFormat,
///                                             kTabularFormat, rid);
/// ```
/// NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:725-743
/// ```c
///     } else if (tabular) {
///         out << "# " << label << ": ";
///     } else {
///         out << label << "= ";
///     }
///
///     string all_id_str = GetSeqIdString(cbs, believe_query);
///     all_id_str += " ";
///     all_id_str = NStr::TruncateSpaces(all_id_str + GetSeqDescrString(cbs));
///
///     // For tabular output, there is no limit on the line length.
///     // There is also no extra line with the sequence length.
///     if (tabular) {
///         out << all_id_str;
/// ```
/// The line ends with the `\n` that `x_PrintQueryAndDbNames` writes before `# Database:`
/// (tabular.cpp:1310-1311).
pub fn write_outfmt7_query_line<W: Write + ?Sized>(
    writer: &mut W,
    query_title: &[u8],
) -> io::Result<()> {
    writer.write_all(b"# Query: ")?;
    writer.write_all(truncate_spaces(query_title))?;
    writer.write_all(b"\n")
}

// =============================================================================
// outfmt 7: Tabular with comment lines (NCBI BLAST compatible)
// Reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:PrintHeader()
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ==========================================================================
    // E-value formatting tests (NCBI parity)
    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:965-988
    // ```c
    // if (evalue < 1.0e-99) {
    //     snprintf(evalue_buf, sizeof(evalue_buf), "%2.0le", evalue);
    // } else if (evalue < 0.0009) {
    //     snprintf(evalue_buf, sizeof(evalue_buf), "%3.0le", evalue);
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1328-1337
    // ```c
    // CAlignFormatUtil::GetScoreString(evalue, bit_score, 0, score, m_Evalue, ...);
    // if ((evalue >= 1.0e-180) && (evalue < 0.0009)){
    //     m_Evalue = NStr::DoubleToString(evalue, 2, NStr::fDoubleScientific);
    // }
    // ```
    // ==========================================================================
    // ==========================================================================

    #[test]
    fn test_format_evalue_ncbi_zero() {
        // Values at or below 0.0 should return "0.0"
        assert_eq!(format_evalue_ncbi_tabular(0.0), "0.0");
    }

    #[test]
    fn test_format_evalue_ncbi_very_small() {
        // < 1e-180: "0.0"
        assert_eq!(format_evalue_ncbi_tabular(1e-200), "0.0");
        assert_eq!(format_evalue_ncbi_tabular(1e-181), "0.0");
    }

    #[test]
    fn test_format_evalue_ncbi_small_scientific() {
        // [1e-180, 0.0009): tabular override uses "%.2e"
        assert_eq!(format_evalue_ncbi_tabular(1e-100), "1.00e-100");
        assert_eq!(format_evalue_ncbi_tabular(1e-150), "1.00e-150");
    }

    #[test]
    fn test_format_evalue_ncbi_medium_scientific() {
        // [1e-180, 0.0009): "%.2e" format (tabular.cpp override)
        assert_eq!(format_evalue_ncbi_tabular(1e-50), "1.00e-50");
        assert_eq!(format_evalue_ncbi_tabular(1.23e-20), "1.23e-20");
        assert_eq!(format_evalue_ncbi_tabular(5e-5), "5.00e-05");
    }

    #[test]
    fn test_format_evalue_ncbi_decimal() {
        // [0.0009, 0.1): "%4.3lf" - 3 decimal places
        assert_eq!(format_evalue_ncbi_tabular(0.001), "0.001");
        assert_eq!(format_evalue_ncbi_tabular(0.005), "0.005");
        assert_eq!(format_evalue_ncbi_tabular(0.099), "0.099");

        // [0.1, 1.0): "%3.2lf" - 2 decimal places
        assert_eq!(format_evalue_ncbi_tabular(0.1), "0.10");
        assert_eq!(format_evalue_ncbi_tabular(0.5), "0.50");
        assert_eq!(format_evalue_ncbi_tabular(0.99), "0.99");

        // [1.0, 10.0): "%2.1lf" - 1 decimal place
        assert_eq!(format_evalue_ncbi_tabular(1.0), "1.0");
        assert_eq!(format_evalue_ncbi_tabular(5.5), "5.5");
        assert_eq!(format_evalue_ncbi_tabular(9.9), "9.9");
    }

    #[test]
    fn test_format_evalue_ncbi_integer() {
        // >= 10.0: "%2.0lf" - no decimal
        assert_eq!(format_evalue_ncbi_tabular(10.0), "10");
        assert_eq!(format_evalue_ncbi_tabular(100.0), "100");
        assert_eq!(format_evalue_ncbi_tabular(1000.0), "1000");
    }

    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:965-975
    // ```c
    // if (evalue < 1.0e-99) {
    //     snprintf(evalue_buf, sizeof(evalue_buf), "%2.0le", evalue);
    // } else if (evalue < 0.0009) {
    //     snprintf(evalue_buf, sizeof(evalue_buf), "%3.0le", evalue);
    // }
    // ```
    #[test]
    fn test_format_evalue_ncbi_base_scientific() {
        // Base GetScoreString formatting uses no decimals for < 0.0009
        assert_eq!(format_evalue_ncbi(1e-100), "1e-100");
        assert_eq!(format_evalue_ncbi(1e-50), "1e-50");
    }

    // ==========================================================================
    // Bit score formatting tests (NCBI parity)
    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:986-993
    // ```c
    // if (bit_score > 99999){
    //     snprintf(bit_score_buf, sizeof(bit_score_buf), "%5.3le", bit_score);
    // } else if (bit_score > 99.9){
    //     snprintf(bit_score_buf, sizeof(bit_score_buf), "%3.0ld",
    //         (long)bit_score);
    // } else {
    //     snprintf(bit_score_buf, sizeof(bit_score_buf), kBitScoreFormat.c_str(),
    //         bit_score);
    // }
    // ```
    // ==========================================================================

    #[test]
    fn test_format_bitscore_ncbi_small() {
        // <= 99.9: "%4.1lf" - one decimal place
        assert_eq!(format_bitscore_ncbi(0.0), " 0.0");
        assert_eq!(format_bitscore_ncbi(9.8), " 9.8");
        let mut buffer = Vec::new();
        write_bitscore_ncbi(&mut buffer, 9.8).unwrap();
        assert_eq!(buffer, b" 9.8");
        assert_eq!(format_bitscore_ncbi(50.5), "50.5");
        assert_eq!(format_bitscore_ncbi(99.9), "99.9");
    }

    #[test]
    fn test_format_bitscore_ncbi_medium() {
        // > 99.9 and <= 99999: "%3.0ld" - integer (no decimal)
        assert_eq!(format_bitscore_ncbi(100.0), "100");
        // "%3.0ld" pads (long)99.95 = 99 to three characters.
        assert_eq!(format_bitscore_ncbi(99.95), " 99");
        let mut buffer = Vec::new();
        write_bitscore_ncbi(&mut buffer, 99.95).unwrap();
        assert_eq!(buffer, b" 99");
        assert_eq!(format_bitscore_ncbi(185.5), "185");
        assert_eq!(format_bitscore_ncbi(598.9), "598");
        assert_eq!(format_bitscore_ncbi(692.0), "692");
        assert_eq!(format_bitscore_ncbi(99999.0), "99999");
    }

    #[test]
    fn test_format_bitscore_ncbi_large() {
        // > 99999: "%5.3le" - scientific notation
        assert_eq!(format_bitscore_ncbi(100000.0), "1.000e+05");
        assert_eq!(format_bitscore_ncbi(123456.789), "1.235e+05");
    }

    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:965-988
    // ```c
    // if (evalue < 1.0e-99) {
    //     snprintf(evalue_buf, sizeof(evalue_buf), "%2.0le", evalue);
    // } else if (evalue < 0.0009) {
    //     snprintf(evalue_buf, sizeof(evalue_buf), "%3.0le", evalue);
    // }
    // if (bit_score > 99999){
    //     snprintf(bit_score_buf, sizeof(bit_score_buf), "%5.3le", bit_score);
    // }
    // ```
    #[test]
    fn test_exponent_padding_without_generated_fixture() {
        assert_eq!(format_evalue_ncbi_tabular(5e-5), "5.00e-05");
        assert_eq!(format_bitscore_ncbi(100000.0), "1.000e+05");
    }

    // ==========================================================================
    // Legacy tests
    // ==========================================================================

    #[test]
    fn test_format_percent_identity_ncbi_divides_before_scaling() {
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/objtools/align_format/tabular.hpp:437-441
        // ```c
        // double perc_ident =
        //     (m_AlignLength > 0 ? ((double)m_NumIdent)/m_AlignLength * 100 : 0);
        // m_Ostream << NStr::DoubleToString(perc_ident, 3);
        // ```
        assert_eq!(format_percent_identity_ncbi(179, 320, 3), "55.937");
        assert_eq!(format_percent_identity_ncbi(0, 0, 3), "0.000");
    }

    // ==========================================================================
    // OutputFormat tests
    // ==========================================================================

    #[test]
    fn test_output_format_parse() {
        assert_eq!(
            OutputFormat::parse("0").unwrap(),
            (OutputFormat::Pairwise, None)
        );
        assert_eq!(
            OutputFormat::parse("6").unwrap(),
            (OutputFormat::Tabular, None)
        );
        assert_eq!(
            OutputFormat::parse("7").unwrap(),
            (OutputFormat::TabularWithComments, None)
        );

        // With custom fields
        let (fmt, fields) = OutputFormat::parse("6 qaccver saccver").unwrap();
        assert_eq!(fmt, OutputFormat::Tabular);
        assert_eq!(fields, Some("qaccver saccver".to_string()));

        // Invalid format
        assert!(OutputFormat::parse("99").is_err());
    }

    // NCBI BLAST+ 2.17.0, `blastn -outfmt 7` (inventory range RP,
    // `scratch_RP/out/<case>/blastn_f7.out`): the query's title as the reader keeps it,
    // raw, after `# Query: `.
    #[test]
    fn outfmt7_query_lines_follow_the_ncbi_oracle() {
        use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};
        for (defline, line) in [
            (&b"id desc"[..], &b"# Query: id desc\n"[..]),
            (b"id\tdesc", b"# Query: id\n"),
            (b"", b"# Query: \n"),
            (b"id\xe9x desc", b"# Query: id\xe9x desc\n"),
            (
                b"id\xffx desc \xff end",
                b"# Query: id\xffx desc \xff end\n",
            ),
            (
                b"id &amp; &lt;b&gt; \"q\" desc",
                b"# Query: id &amp; &lt;b&gt; \"q\" desc\n",
            ),
            (b"\x01id desc", b"# Query: \x01id desc\n"),
            (b"id desc   ", b"# Query: id desc\n"),
            (b"TPA: id desc", b"# Query: TPA: id desc\n"),
            (b"id desc\x01more", b"# Query: id desc\n"),
        ] {
            let bytes = [&b">"[..], defline, b"\nACGT\n"].concat();
            let mut source =
                FastaInputSource::from_bytes(&bytes, ReaderConfig::query("BLASTN", false, false));
            let title = read_all(&mut source, &mut |_| Ok(())).unwrap()[0]
                .title
                .clone();
            let mut out = Vec::new();
            write_outfmt7_query_line(&mut out, &title).unwrap();
            assert_eq!(out, line, "{defline:?}");
        }
        // `NStr::TruncateSpaces`: 0x09-0x0D and the space at both ends.
        assert_eq!(truncate_spaces(b" \x0b\x0cid x\t\r "), b"id x");
        assert_eq!(truncate_spaces(b" \n "), b"");
        assert_eq!(truncate_spaces(b"\xa0id\xa0"), b"\xa0id\xa0");
    }
}
