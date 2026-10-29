//! The tabular reports with a field list (outfmt 6, 7, 10 and 20 with `-outfmt "<n> <fields>"`
//! or `delim=`): NCBI's `CBlastTabularInfo`, shared by the programs. A program's adapter hands
//! over, for each HSP of a query's final hit list in report order, the values of the HSP
//! (`TabularHsp`), its subject (`TabularSubject`) and its query (`TabularQuery`); this module
//! reads the field list, writes the comment lines of outfmt 7, the header line of outfmt 20 and
//! the rows.
//!
//! The default reports (outfmt 6 and 7 with the `std` fields and a tab) keep each program's own
//! writer: `TabularSpec::of` returns `None` for them.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:770-786
//! ```c++
//!     if (m_FormatType == CFormattingArgs::eTabular ||
//!         m_FormatType == CFormattingArgs::eTabularWithComments ||
//!         m_FormatType == CFormattingArgs::eCommaSeparatedValues ||
//!         m_FormatType == CFormattingArgs::eCommaSeparatedValuesWithHeader) {
//!
//!       const CBlastTabularInfo::EFieldDelimiter kDelim =
//!             ((m_FormatType == CFormattingArgs::eCommaSeparatedValues ||
//!               m_FormatType == CFormattingArgs::eCommaSeparatedValuesWithHeader)
//!              ? CBlastTabularInfo::eComma : CBlastTabularInfo::eTab);
//!
//!         CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec, kDelim);
//!         if(!m_CustomDelim.empty()) {
//!             tabinfo.SetCustomDelim(m_CustomDelim);
//!         }
//! ```

use std::io::{self, Write};

use crate::blastinput::app::{FormatChoice, ReportFormat};
use crate::report::blast_xml::matrix::FormattingMatrix;
use crate::report::outfmt6::{
    format_bitscore_ncbi, format_evalue_ncbi_tabular, format_percent_identity_ncbi,
    write_outfmt7_query_line,
};

/// NCBI's `ETabularField`: the fields of a tabular report.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum TabularField {
    QuerySeqId,
    QueryGi,
    QueryAccession,
    QueryAccessionVersion,
    QueryLength,
    SubjectSeqId,
    SubjectAllSeqIds,
    SubjectGi,
    SubjectAllGis,
    SubjectAccession,
    SubjectAccessionVersion,
    SubjectAllAccessions,
    SubjectLength,
    QueryStart,
    QueryEnd,
    SubjectStart,
    SubjectEnd,
    QuerySeq,
    SubjectSeq,
    Evalue,
    BitScore,
    Score,
    AlignmentLength,
    PercentIdentical,
    NumIdentical,
    Mismatches,
    Positives,
    GapOpenings,
    Gaps,
    PercentPositives,
    Frames,
    QueryFrame,
    SubjectFrame,
    Btop,
    SubjectTaxId,
    SubjectSciName,
    SubjectCommonName,
    SubjectBlastName,
    SubjectSuperKingdom,
    SubjectTaxIds,
    SubjectSciNames,
    SubjectCommonNames,
    SubjectBlastNames,
    SubjectSuperKingdoms,
    SubjectTitle,
    SubjectAllTitles,
    SubjectStrand,
    QueryCovSubject,
    QueryCovHsp,
    QueryCovUniqSubject,
}

use TabularField::*;

/// The field names of `-outfmt` (`sc_FormatSpecifiers`) and the names of outfmt 7's
/// `# Fields:` line (`PrintFieldNames`), by field.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:47-198
/// ```c++
/// const SFormatSpec sc_FormatSpecifiers[kNumTabularOutputFormatSpecifiers] = {
///     SFormatSpec("qseqid",
///                 "Query Seq-id",
///                 eQuerySeqId),
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1111-1236
/// ```c++
///         switch (*iter) {
///         case eQuerySeqId:
///             m_Ostream << "query id"; break;
/// ```
const FIELDS: [(&str, TabularField, &str); 50] = [
    ("qseqid", QuerySeqId, "query id"),
    ("qgi", QueryGi, "query gi"),
    ("qacc", QueryAccession, "query acc."),
    ("qaccver", QueryAccessionVersion, "query acc.ver"),
    ("qlen", QueryLength, "query length"),
    ("sseqid", SubjectSeqId, "subject id"),
    ("sallseqid", SubjectAllSeqIds, "subject ids"),
    ("sgi", SubjectGi, "subject gi"),
    ("sallgi", SubjectAllGis, "subject gis"),
    ("sacc", SubjectAccession, "subject acc."),
    ("saccver", SubjectAccessionVersion, "subject acc.ver"),
    ("sallacc", SubjectAllAccessions, "subject accs."),
    ("slen", SubjectLength, "subject length"),
    ("qstart", QueryStart, "q. start"),
    ("qend", QueryEnd, "q. end"),
    ("sstart", SubjectStart, "s. start"),
    ("send", SubjectEnd, "s. end"),
    ("qseq", QuerySeq, "query seq"),
    ("sseq", SubjectSeq, "subject seq"),
    ("evalue", Evalue, "evalue"),
    ("bitscore", BitScore, "bit score"),
    ("score", Score, "score"),
    ("length", AlignmentLength, "alignment length"),
    ("pident", PercentIdentical, "% identity"),
    ("nident", NumIdentical, "identical"),
    ("mismatch", Mismatches, "mismatches"),
    ("positive", Positives, "positives"),
    ("gapopen", GapOpenings, "gap opens"),
    ("gaps", Gaps, "gaps"),
    ("ppos", PercentPositives, "% positives"),
    ("frames", Frames, "query/sbjct frames"),
    ("qframe", QueryFrame, "query frame"),
    ("sframe", SubjectFrame, "sbjct frame"),
    ("btop", Btop, "BTOP"),
    ("staxid", SubjectTaxId, "subject tax id"),
    ("ssciname", SubjectSciName, "subject sci name"),
    ("scomname", SubjectCommonName, "subject com names"),
    ("sblastname", SubjectBlastName, "subject blast name"),
    ("sskingdom", SubjectSuperKingdom, "subject super kingdom"),
    ("staxids", SubjectTaxIds, "subject tax ids"),
    ("sscinames", SubjectSciNames, "subject sci names"),
    ("scomnames", SubjectCommonNames, "subject com names"),
    ("sblastnames", SubjectBlastNames, "subject blast names"),
    ("sskingdoms", SubjectSuperKingdoms, "subject super kingdoms"),
    ("stitle", SubjectTitle, "subject title"),
    ("salltitles", SubjectAllTitles, "subject titles"),
    ("sstrand", SubjectStrand, "subject strand"),
    ("qcovs", QueryCovSubject, "% query coverage per subject"),
    ("qcovhsp", QueryCovHsp, "% query coverage per hsp"),
    (
        "qcovus",
        QueryCovUniqSubject,
        "% query coverage per uniq subject",
    ),
];

/// The fields of `std` and of a list without a known field.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:38-43
/// ```c++
/// const char* kDfltArgTabularOutputFmt =
///     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
///     "evalue bitscore";
/// const char* kDfltArgTabularOutputFmtTag("std");
/// ```
pub const STD_FIELDS: [TabularField; 12] = [
    QueryAccessionVersion,
    SubjectAccessionVersion,
    PercentIdentical,
    AlignmentLength,
    Mismatches,
    GapOpenings,
    QueryStart,
    QueryEnd,
    SubjectStart,
    SubjectEnd,
    Evalue,
    BitScore,
];

impl TabularField {
    /// The field of an `-outfmt` field name.
    pub fn from_name(name: &str) -> Option<Self> {
        FIELDS
            .iter()
            .find(|(spec, _, _)| *spec == name)
            .map(|&(_, field, _)| field)
    }

    fn entry(self) -> &'static (&'static str, TabularField, &'static str) {
        &FIELDS[self as usize]
    }

    /// The `-outfmt` name (outfmt 20's header line).
    pub fn name(self) -> &'static str {
        self.entry().0
    }

    /// The name on outfmt 7's `# Fields:` line.
    pub fn header_name(self) -> &'static str {
        self.entry().2
    }
}

/// NCBI's field list of a custom specification (`x_SetFieldsToShow`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:70-99
/// ```c++
///     vector<string> format_tokens;
///     NStr::Split(format, " ", format_tokens);
///
///     if (format_tokens.empty())
///         x_AddDefaultFieldsToShow();
///
///     ITERATE (vector<string>, iter, format_tokens) {
///         if (*iter == kDfltArgTabularOutputFmtTag)
///             x_AddDefaultFieldsToShow();
///         else if ((*iter)[0] == '-') {
///             string field = (*iter).substr(1);
///             if (m_FieldMap.count(field) > 0)
///                 x_DeleteFieldToShow(m_FieldMap[field]);
///         } else {
///             if (m_FieldMap.count(*iter) > 0)
///                 x_AddFieldToShow(m_FieldMap[*iter]);
///         }
///     }
///
///     if (m_FieldsToShow.empty()) {
///         x_AddDefaultFieldsToShow();
///     }
/// ```
pub fn parse_fields(spec: &str) -> Vec<TabularField> {
    let mut fields: Vec<TabularField> = Vec::new();
    let mut add = |fields: &mut Vec<TabularField>, field| {
        if !fields.contains(&field) {
            fields.push(field);
        }
    };
    for token in spec.split(' ') {
        if token == "std" {
            for field in STD_FIELDS {
                add(&mut fields, field);
            }
        } else if let Some(name) = token.strip_prefix('-') {
            if let Some(field) = TabularField::from_name(name) {
                fields.retain(|shown| *shown != field);
            }
        } else if let Some(field) = TabularField::from_name(token) {
            add(&mut fields, field);
        }
    }
    if fields.is_empty() {
        fields = STD_FIELDS.to_vec();
    }
    fields
}

/// NCBI's check that a custom delimiter is not the separator of a field's entries.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2753-2792
/// ```c++
/// static void s_ValidateCustomDelim(const string& customFmtSpec,const string& customDelim)
/// {
///     bool error = false;
///     string checkfield;
///     auto custom_fmt_spec = NStr::TruncateSpaces(customFmtSpec);
///     if(custom_fmt_spec.empty()) return;
///
///     //Check if delim is already used
///     const string kFieldsWithSemicolSeparator = "sallseqid staxids sscinames scomnames sblastnames sskingdoms";//sep = ";"
///     const string kFramesField = "frames"; //sep = "/"
///     const string kAllTitlesField ="salltitles"; //sep = "<>""
///
///     if(customDelim == ";") {
///         ...
///             if(NStr::Find(custom_fmt_spec,tokens[i]) != NPOS) {
///     ...
///     if(error) {
///         string msg("Your custom record separator (" + customDelim + ") is also used by the format specifier (" + checkfield +
///                                                                 ") to separate multiple entries. Please use a different record separator (delim keyword).");
///         NCBI_THROW(CInputException, eInvalidInput, msg);
/// ```
/// Returns NCBI's message; the field list (`spec`) is searched as text, as NCBI does.
pub fn custom_delimiter_conflict(spec: &str, delimiter: &str) -> Option<String> {
    let spec = super::outfmt6::truncate_spaces(spec.as_bytes());
    if spec.is_empty() {
        return None;
    }
    let contains = |needle: &str| spec.windows(needle.len()).any(|w| w == needle.as_bytes());
    let field = match delimiter {
        ";" => [
            "sallseqid",
            "staxids",
            "sscinames",
            "scomnames",
            "sblastnames",
            "sskingdoms",
        ]
        .into_iter()
        .find(|field| contains(field)),
        "/" => Some("frames").filter(|field| contains(field)),
        "<>" => Some("salltitles").filter(|field| contains(field)),
        _ => None,
    }?;
    Some(format!(
        "Your custom record separator ({delimiter}) is also used by the format specifier ({field}) to separate multiple entries. Please use a different record separator (delim keyword)."
    ))
}

/// Which of the tabular reports.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TabularKind {
    /// outfmt 6.
    Tabular,
    /// outfmt 7: comment lines before each query's rows and after the last query.
    Commented,
    /// outfmt 10.
    Csv,
    /// outfmt 20: outfmt 10 after a line of the field names.
    CsvWithHeader,
}

/// A tabular report with its field list and delimiter.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TabularSpec {
    pub kind: TabularKind,
    pub fields: Vec<TabularField>,
    /// The field delimiter of the rows (and of outfmt 20's header line).
    pub delimiter: String,
}

/// What a program's adapter must give the writer besides the values every row has.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct TabularNeeds {
    /// The alignment rows (`TabularHsp::rows`): `qseq`, `sseq`, `btop`, `positive`, `ppos`.
    pub rows: bool,
    /// `qcovs` (`TabularSubject::qcovs`).
    pub qcovs: bool,
    /// `qcovus` (`TabularSubject::qcovus`).
    pub qcovus: bool,
    /// Any of the taxonomy fields (`TabularSubject::taxonomy`).
    pub taxonomy: bool,
    /// The taxonomy names (singular or plural), which need the taxonomy database.
    pub taxonomy_names: bool,
    /// The fields of every definition line (`sallseqid`, `sallgi`, `sallacc`, `salltitles`).
    pub all_deflines: bool,
}

impl TabularSpec {
    /// The tabular report of `format` with the custom specification of `choice`, or `None` for
    /// a format that is not tabular and for the default outfmt 6 and 7 (the `std` fields and a
    /// tab), which the programs write with their own writers.
    pub fn of(format: ReportFormat, choice: &FormatChoice) -> Option<Self> {
        let (kind, default_delimiter) = match format {
            ReportFormat::Tabular => (TabularKind::Tabular, "\t"),
            ReportFormat::TabularWithComments => (TabularKind::Commented, "\t"),
            ReportFormat::CommaSeparatedValues => (TabularKind::Csv, ","),
            ReportFormat::CommaSeparatedValuesWithHeader => (TabularKind::CsvWithHeader, ","),
            _ => return None,
        };
        let fields = parse_fields(&choice.spec);
        let delimiter = if choice.delimiter.is_empty() {
            default_delimiter.to_string()
        } else {
            choice.delimiter.clone()
        };
        let default_report = matches!(kind, TabularKind::Tabular | TabularKind::Commented)
            && fields == STD_FIELDS
            && delimiter == "\t";
        (!default_report).then_some(TabularSpec {
            kind,
            fields,
            delimiter,
        })
    }

    /// The spec, kept for the rest of the process (one per distinct spec), so that the
    /// programs' `Copy` format values can hold it.
    pub fn intern(self) -> &'static TabularSpec {
        static SPECS: std::sync::Mutex<Vec<&'static TabularSpec>> =
            std::sync::Mutex::new(Vec::new());
        let mut specs = SPECS
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        if let Some(spec) = specs.iter().find(|spec| ***spec == self) {
            return spec;
        }
        let spec: &'static TabularSpec = Box::leak(Box::new(self));
        specs.push(spec);
        spec
    }

    pub fn has(&self, field: TabularField) -> bool {
        self.fields.contains(&field)
    }

    fn has_any(&self, fields: &[TabularField]) -> bool {
        fields.iter().any(|field| self.has(*field))
    }

    pub fn needs(&self) -> TabularNeeds {
        const NAMES: [TabularField; 8] = [
            SubjectSciName,
            SubjectCommonName,
            SubjectBlastName,
            SubjectSuperKingdom,
            SubjectSciNames,
            SubjectCommonNames,
            SubjectBlastNames,
            SubjectSuperKingdoms,
        ];
        TabularNeeds {
            rows: self.has_any(&[QuerySeq, SubjectSeq, Btop, Positives, PercentPositives]),
            qcovs: self.has(QueryCovSubject),
            qcovus: self.has(QueryCovUniqSubject),
            taxonomy: self.has_any(&NAMES) || self.has_any(&[SubjectTaxId, SubjectTaxIds]),
            taxonomy_names: self.has_any(&NAMES),
            all_deflines: self.has_any(&[
                SubjectAllSeqIds,
                SubjectAllGis,
                SubjectAllAccessions,
                SubjectAllTitles,
            ]),
        }
    }

    /// Whether NCBI writes the taxonomy database's warning each time it builds the report's
    /// formatter (`x_CheckTaxDB`: only for the plural names) when `taxdb.bti` is not found.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1093-1104
    /// ```c++
    /// void CBlastTabularInfo::x_CheckTaxDB()
    /// {
    ///    	if( x_IsFieldRequested(eSubjectSciNames) ||
    ///     	x_IsFieldRequested(eSubjectCommonNames) ||
    ///     	x_IsFieldRequested(eSubjectBlastNames) ||
    ///     	x_IsFieldRequested(eSubjectSuperKingdoms))
    ///    	{
    ///    		string resolved = SeqDB_ResolveDbPath("taxdb.bti");
    ///    		if(resolved.empty())
    ///    			ERR_POST(Warning << "Taxonomy name lookup from taxid requires installation of taxdb database with ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz");
    ///    	}
    /// }
    /// ```
    pub fn checks_taxdb(&self) -> bool {
        self.has_any(&[
            SubjectSciNames,
            SubjectCommonNames,
            SubjectBlastNames,
            SubjectSuperKingdoms,
        ])
    }

    /// `x_IsCsv`: the names are quoted when the delimiter is a comma.
    fn is_csv(&self) -> bool {
        self.delimiter == ","
    }
}

/// The warning of `x_CheckTaxDB` (`TabularSpec::checks_taxdb`), as NCBI's `ERR_POST` writes it
/// for `program` (lower case).
pub fn taxdb_warning_line(program: &str) -> String {
    format!(
        "Warning: [{program}] Taxonomy name lookup from taxid requires installation of taxdb database with ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz\n"
    )
}

/// A query of the report.
#[derive(Clone, Copy, Debug)]
pub struct TabularQuery<'a> {
    /// `qseqid`, `qacc` and `qaccver`: LOSAT's queries have local ids, shown as the first word of
    /// the title (`FastaRecord::shown_id`); `qgi` is 0.
    pub id: &'a [u8],
    /// `qlen`: the searched length (of the `-query_loc` interval).
    pub length: usize,
}

/// The ids of one definition line of a subject (`sallseqid`, `sallgi`, `sallacc`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DeflineIds {
    pub id: Vec<u8>,
    pub gi: i64,
    pub accession: Vec<u8>,
}

/// The taxonomy of a subject (`x_SetTaxInfo`, `x_SetTaxInfoAll`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct TabularTaxonomy {
    /// `staxid`: the first non-zero taxid of the definition lines (0 for none: `N/A`).
    pub taxid: i32,
    /// `staxids`: the taxids of the definition lines, sorted, without repeats (a database
    /// sequence without a taxid has `{0}`; a FASTA subject none: `N/A`).
    pub taxids: Vec<i32>,
}

impl TabularTaxonomy {
    /// The taxonomy of a database sequence with the definition lines `deflines`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:590-603
    /// ```c++
    ///     if (bdlRef.NotEmpty() && bdlRef->CanGet() && bdlRef->IsSet() && !bdlRef->Get().empty()){
    ///     	ITERATE(CBlast_def_line_set::Tdata, itr, bdlRef->Get()) {
    ///             if((*itr)->IsSetTaxid()) {
    ///                     if((*itr)->GetTaxid() != ZERO_TAX_ID) {
    ///                     	m_SubjectTaxId = (*itr)->GetTaxid();
    ///                     	break;
    ///                     }
    ///             }
    ///     	}
    ///     }
    ///
    ///     if(m_SubjectTaxId == ZERO_TAX_ID) {
    ///           m_SubjectTaxId = sequence::GetTaxId(handle);
    ///     }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:647-655
    /// ```c++
    ///     if (bdlRef.NotEmpty() && bdlRef->CanGet() && bdlRef->IsSet() && !bdlRef->Get().empty()){
    ///     	ITERATE(CBlast_def_line_set::Tdata, itr, bdlRef->Get()) {
    ///     		CBlast_def_line::TTaxIds t = (*itr)->GetTaxIds();
    ///     		m_SubjectTaxIds.insert(t.begin(), t.end());
    ///     	}
    ///     }
    ///
    ///     if(m_SubjectTaxIds.empty()) {
    ///             CSeqdesc_CI desc_s(handle, CSeqdesc::e_Source);
    /// ```
    /// The fallbacks read the Bioseq's sources, which the BLAST database loader makes of the
    /// definition lines' positive `taxid` members only (`CSeqDBVol::x_GetTaxonomy`), so they add
    /// nothing that the definition lines did not give.
    pub fn of_deflines(deflines: &[crate::seqid::BlastDefline]) -> Self {
        let taxid = deflines
            .iter()
            .filter_map(|defline| defline.taxid)
            .find(|&taxid| taxid != 0)
            .unwrap_or(0);
        let mut taxids = std::collections::BTreeSet::new();
        for defline in deflines {
            taxids.extend(defline.tax_ids());
        }
        TabularTaxonomy {
            taxid,
            taxids: taxids.into_iter().collect(),
        }
    }
}

/// A subject of the report.
#[derive(Clone, Copy, Debug)]
pub struct TabularSubject<'a> {
    /// `sseqid`, `sacc`, `saccver` and `sgi`.
    pub id: &'a [u8],
    pub accession: &'a [u8],
    pub accession_version: &'a [u8],
    pub gi: i64,
    /// The ids of each definition line; empty for a subject without definition lines (a FASTA
    /// subject or a database sequence without parsed ids), whose `sall*` are its own ids.
    pub all_ids: &'a [DeflineIds],
    /// `slen`.
    pub length: usize,
    /// The titles of the definition lines (`stitle` = the first, `salltitles` joined by `<>`,
    /// an empty one `N/A`); `None` for a subject without definition lines (a FASTA subject):
    /// `N/A`.
    pub titles: Option<&'a [&'a [u8]]>,
    pub taxonomy: &'a TabularTaxonomy,
    /// Whether the subject is a nucleotide sequence (`sstrand` is `plus` / `minus`, else `N/A`).
    pub nucleotide: bool,
    /// `qcovs` and `qcovus` of the subject's HSPs; -1 for `N/A`.
    pub qcovs: i32,
    pub qcovus: i32,
}

/// An HSP of the report, with the values as NCBI's tabular report shows them.
#[derive(Clone, Copy, Debug)]
pub struct TabularHsp<'a> {
    /// `qstart`, `qend`, `sstart`, `send` (1-based; a minus strand or frame has start > end).
    pub q_start: i64,
    pub q_end: i64,
    pub s_start: i64,
    pub s_end: i64,
    pub evalue: f64,
    pub bit_score: f64,
    pub score: i32,
    /// `length`, `gapopen`, `gaps`: the alignment's columns, gap openings and gap letters.
    pub length: usize,
    pub gap_opens: usize,
    pub gaps: usize,
    /// The engine's identities (counted again from `rows` when they are given).
    pub num_ident: usize,
    /// `qframe` and `sframe` as NCBI shows them.
    pub q_frame: i32,
    pub s_frame: i32,
    /// The query and subject rows of the alignment (`-` for a gap; translated rows in
    /// residues), when `TabularNeeds::rows`.
    pub rows: Option<(&'a [u8], &'a [u8])>,
    /// `qcovhsp` (the integer of NCBI's `hsp_percent_coverage`).
    pub qcovhsp: i32,
}

/// `NStr::Quote` with its default quote (`"`) and escape (`\`) characters.
fn write_quoted(out: &mut Vec<u8>, text: &[u8]) {
    if text.is_empty() {
        return;
    }
    out.push(b'"');
    for &byte in text {
        if byte == b'"' || byte == b'\\' {
            out.push(b'\\');
        }
        out.push(byte);
    }
    out.push(b'"');
}

const NA: &[u8] = b"N/A";

/// The writer of a tabular report.
pub struct TabularWriter<'a> {
    pub spec: &'a TabularSpec,
    /// The `# BLASTN 2.17.0+` line's program and version.
    pub program_version: &'a str,
    /// The positive pairs of a protein search (`None` for BLASTN: positives are identities).
    pub matrix: Option<&'a FormattingMatrix>,
}

impl TabularWriter<'_> {
    /// outfmt 20's line of the field names (`PrintProlog` → `PrintFieldSpecs`).
    pub fn prolog(&self, out: &mut dyn Write) -> io::Result<()> {
        if self.spec.kind != TabularKind::CsvWithHeader {
            return Ok(());
        }
        let mut line = Vec::new();
        for (i, field) in self.spec.fields.iter().enumerate() {
            if i > 0 {
                line.extend_from_slice(self.spec.delimiter.as_bytes());
            }
            line.extend_from_slice(field.name().as_bytes());
        }
        line.push(b'\n');
        out.write_all(&line)
    }

    /// outfmt 7's lines before a query's rows (`PrintHeader`): the program, the query's title,
    /// the database, and when `num_hits` is known the fields (if there are hits) and the count.
    pub fn query_header(
        &self,
        out: &mut dyn Write,
        query_title: &[u8],
        database: &[u8],
        num_hits: Option<usize>,
    ) -> io::Result<()> {
        if self.spec.kind != TabularKind::Commented {
            return Ok(());
        }
        writeln!(out, "# {}", self.program_version)?;
        write_outfmt7_query_line(out, query_title)?;
        out.write_all(b"# Database: ")?;
        out.write_all(database)?;
        out.write_all(b"\n")?;
        let Some(num_hits) = num_hits else {
            return Ok(());
        };
        if num_hits > 0 {
            let mut line = b"# Fields: ".to_vec();
            for (i, field) in self.spec.fields.iter().enumerate() {
                if i > 0 {
                    line.extend_from_slice(b", ");
                }
                line.extend_from_slice(field.header_name().as_bytes());
            }
            line.push(b'\n');
            out.write_all(&line)?;
        }
        writeln!(out, "# {num_hits} hits found")
    }

    /// outfmt 7's last line (`PrintNumProcessed`).
    pub fn epilog(&self, out: &mut dyn Write, num_queries: usize) -> io::Result<()> {
        if self.spec.kind != TabularKind::Commented {
            return Ok(());
        }
        writeln!(out, "# BLAST processed {num_queries} queries")
    }

    /// One row (`SetFields` then `Print`), appended to `out`.
    pub fn row(
        &self,
        out: &mut Vec<u8>,
        query: &TabularQuery<'_>,
        subject: &TabularSubject<'_>,
        hsp: &TabularHsp<'_>,
    ) {
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:973-1027
        // (the identities and positives counted again from the rows, and the BTOP string).
        let (num_ident, positives) = match hsp.rows {
            Some((qrow, srow)) => {
                let (mut ident, mut positive) = (0usize, 0usize);
                for (&q, &s) in qrow.iter().zip(srow) {
                    if q == s {
                        ident += 1;
                        positive += 1;
                    } else if self.matrix.is_some_and(|matrix| matrix.positive(q, s)) {
                        positive += 1;
                    }
                }
                (ident, positive)
            }
            None => (hsp.num_ident, hsp.num_ident),
        };
        let spec = self.spec;
        let csv = spec.is_csv();
        for (i, field) in spec.fields.iter().enumerate() {
            if i > 0 {
                out.extend_from_slice(spec.delimiter.as_bytes());
            }
            match field {
                QuerySeqId | QueryAccession | QueryAccessionVersion => {
                    out.extend_from_slice(query.id)
                }
                QueryGi => out.push(b'0'),
                QueryLength => push_num(out, query.length),
                SubjectSeqId => out.extend_from_slice(subject.id),
                SubjectAccession => out.extend_from_slice(subject.accession),
                SubjectAccessionVersion => out.extend_from_slice(subject.accession_version),
                SubjectGi => push_num(out, subject.gi),
                SubjectAllSeqIds | SubjectAllGis | SubjectAllAccessions => {
                    if subject.all_ids.is_empty() {
                        match field {
                            SubjectAllSeqIds => out.extend_from_slice(subject.id),
                            SubjectAllGis => push_num(out, subject.gi),
                            _ => out.extend_from_slice(subject.accession),
                        }
                    }
                    for (n, ids) in subject.all_ids.iter().enumerate() {
                        if n > 0 {
                            out.push(b';');
                        }
                        match field {
                            SubjectAllSeqIds => out.extend_from_slice(&ids.id),
                            SubjectAllGis => push_num(out, ids.gi),
                            _ => out.extend_from_slice(&ids.accession),
                        }
                    }
                }
                SubjectLength => push_num(out, subject.length),
                QueryStart => push_num(out, hsp.q_start),
                QueryEnd => push_num(out, hsp.q_end),
                SubjectStart => push_num(out, hsp.s_start),
                SubjectEnd => push_num(out, hsp.s_end),
                QuerySeq => out.extend_from_slice(hsp.rows.map_or(&[][..], |rows| rows.0)),
                SubjectSeq => out.extend_from_slice(hsp.rows.map_or(&[][..], |rows| rows.1)),
                Evalue => out.extend_from_slice(format_evalue_ncbi_tabular(hsp.evalue).as_bytes()),
                BitScore => out.extend_from_slice(format_bitscore_ncbi(hsp.bit_score).as_bytes()),
                Score => push_num(out, hsp.score),
                AlignmentLength => push_num(out, hsp.length),
                PercentIdentical => out.extend_from_slice(
                    format_percent_identity_ncbi(num_ident, hsp.length, 3).as_bytes(),
                ),
                NumIdentical => push_num(out, num_ident),
                // NCBI reference (598d8ae6): c++/include/objtools/align_format/tabular.hpp:543-547
                // ```c++
                //     int num_mismatches = m_AlignLength - m_NumIdent - m_NumGaps;
                // ```
                Mismatches => push_num(out, hsp.length as i64 - num_ident as i64 - hsp.gaps as i64),
                Positives => push_num(out, positives),
                GapOpenings => push_num(out, hsp.gap_opens),
                Gaps => push_num(out, hsp.gaps),
                PercentPositives => out.extend_from_slice(
                    format_percent_identity_ncbi(positives, hsp.length, 2).as_bytes(),
                ),
                Frames => {
                    push_num(out, hsp.q_frame);
                    out.push(b'/');
                    push_num(out, hsp.s_frame);
                }
                QueryFrame => push_num(out, hsp.q_frame),
                SubjectFrame => push_num(out, hsp.s_frame),
                Btop => {
                    if let Some((qrow, srow)) = hsp.rows {
                        write_btop(out, qrow, srow);
                    }
                }
                SubjectTaxId => match subject.taxonomy.taxid {
                    0 => out.extend_from_slice(NA),
                    taxid => push_num(out, taxid),
                },
                SubjectTaxIds => {
                    if subject.taxonomy.taxids.is_empty() {
                        out.extend_from_slice(NA);
                    }
                    for (n, taxid) in subject.taxonomy.taxids.iter().enumerate() {
                        if n > 0 {
                            out.push(b';');
                        }
                        push_num(out, *taxid);
                    }
                }
                // No name of the taxonomy database is known (a program refuses a report that
                // would read one, DECISIONS 3): every lookup fails. With `staxids` among the
                // fields NCBI then fills the scientific and common names with `N/A`: the singular
                // ones for a non-zero `staxid`, the plural ones once per taxid of `staxids`
                // (`x_SetTaxInfo`, `x_SetTaxInfoAll`); a filled singular name and the plural
                // names (blast names too) are quoted in a CSV, an empty singular name and the
                // super kingdoms are not.
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:626-633
                // ```c++
                //         } catch (const CException&) {
                //         	//only put fillers in if we are going to show tax id
                //             // the fillers are put in so that the name list would
                //             // match the taxid list
                //             if(x_IsFieldRequested(eSubjectTaxIds) ) {
                //             	m_SubjectSciName = NA;
                //                 m_SubjectCommonName = NA;
                //             }
                // ```
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:337-359
                // ```c++
                // void CBlastTabularInfo::x_PrintSubjectSciName()
                // {
                // 	if(m_SubjectSciName == kEmptyStr) {
                // 		m_Ostream << NA;
                // 		return;
                // 	}
                // 	m_Ostream << (x_IsCsv() ? NStr::Quote(m_SubjectSciName) : m_SubjectSciName);
                // }
                //
                // void CBlastTabularInfo::x_PrintSubjectSciNames()
                // {
                //     CNcbiOstrstream oss;
                // 	if(m_SubjectSciNames.empty()) {
                // 		oss << NA;
                // 	} else {
                //         ITERATE(vector<string>, iter, m_SubjectSciNames) {
                //              if (iter != m_SubjectSciNames.begin())
                //                  oss << ";";
                //              oss << *iter;
                //         }
                //     }
                //     const string output = CNcbiOstrstreamToString(oss);
                //     m_Ostream << (x_IsCsv() ? NStr::Quote(output) : output);
                // ```
                SubjectSciName | SubjectCommonName => {
                    if subject.taxonomy.taxid != 0 && spec.has(SubjectTaxIds) && csv {
                        write_quoted(out, NA);
                    } else {
                        out.extend_from_slice(NA);
                    }
                }
                SubjectSciNames | SubjectCommonNames | SubjectBlastNames => {
                    let fillers = match field {
                        SubjectBlastNames => 0,
                        _ if spec.has(SubjectTaxIds) => subject.taxonomy.taxids.len(),
                        _ => 0,
                    };
                    let mut names = NA.to_vec();
                    for _ in 1..fillers {
                        names.push(b';');
                        names.extend_from_slice(NA);
                    }
                    if csv {
                        write_quoted(out, &names);
                    } else {
                        out.extend_from_slice(&names);
                    }
                }
                SubjectBlastName | SubjectSuperKingdom | SubjectSuperKingdoms => {
                    out.extend_from_slice(NA)
                }
                SubjectTitle => match subject.titles.and_then(|titles| titles.first()) {
                    Some(title) if !title.is_empty() => out.extend_from_slice(title),
                    _ => out.extend_from_slice(NA),
                },
                SubjectAllTitles => match subject.titles.filter(|titles| !titles.is_empty()) {
                    Some(titles) => {
                        for (n, title) in titles.iter().enumerate() {
                            if n > 0 {
                                out.extend_from_slice(b"<>");
                            }
                            out.extend_from_slice(if title.is_empty() { NA } else { title });
                        }
                    }
                    None => out.extend_from_slice(NA),
                },
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1064-1070
                // ```c++
                //         	if(!subject_is_na)
                //         		m_SubjectStrand = NA;
                //         	else
                //         		m_SubjectStrand = ((s_start - s_end) > 0 )? "minus":"plus";
                // ```
                SubjectStrand => out.extend_from_slice(if !subject.nucleotide {
                    NA
                } else if hsp.s_start - hsp.s_end > 0 {
                    b"minus"
                } else {
                    b"plus"
                }),
                QueryCovSubject => push_coverage(out, subject.qcovs),
                QueryCovUniqSubject => push_coverage(out, subject.qcovus),
                QueryCovHsp => push_coverage(out, hsp.qcovhsp),
            }
        }
        out.push(b'\n');
    }
}

fn push_num(out: &mut Vec<u8>, value: impl std::fmt::Display) {
    let _ = write!(out, "{value}");
}

fn push_coverage(out: &mut Vec<u8>, value: i32) {
    if value < 0 {
        out.extend_from_slice(NA);
    } else {
        push_num(out, value);
    }
}

/// The BTOP string of two rows.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1001-1027
/// ```c++
///             for (unsigned int i = 0;
///                  i < min(m_QuerySeq.size(), m_SubjectSeq.size());
///                  ++i) {
///                 if (m_QuerySeq[i] == m_SubjectSeq[i]) {
///                     ++num_ident;
///                     ++num_positives;
///                     ++num_matches;
///                 } else {
///                     if(num_matches > 0) {
///                         btop_string +=  NStr::Int8ToString(num_matches);
///                         num_matches=0;
///                     }
///                     btop_string += m_QuerySeq[i];
///                     btop_string += m_SubjectSeq[i];
/// ```
pub fn write_btop(out: &mut Vec<u8>, qrow: &[u8], srow: &[u8]) {
    let mut matches = 0usize;
    for (&q, &s) in qrow.iter().zip(srow) {
        if q == s {
            matches += 1;
        } else {
            if matches > 0 {
                push_num(out, matches);
                matches = 0;
            }
            out.push(q);
            out.push(s);
        }
    }
    if matches > 0 {
        push_num(out, matches);
    }
}

/// `seq_percent_coverage` of one subject (`qcovs`) from the query intervals of its HSPs
/// (inclusive, 0-based, any order): the length of their union in percent of `query_length`,
/// 0.5 added below 99, cut to an integer; -1 (`N/A`) when nothing is covered.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:629-640
/// ```c++
///        	if((check_type & eQueryCovPerSubj)) {
///        		int master_coverage = align_format::CAlignFormatUtil::GetMasterCoverage(tmp_align_set);
///        		if (master_coverage)
///         	{
///            		double subj_coverage = 100.0 * (double) master_coverage/ (double) query_len;
///            		if(subj_coverage < 99)
///            			subj_coverage +=0.5;
///            		(*left_it)->SetNamedScore ("seq_percent_coverage", (int) subj_coverage);
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:2940-2985
/// (`s_MergeRangeList`: sorted by start, a range that intersects the last is merged).
pub fn coverage_percent(covered: usize, query_length: usize) -> i32 {
    if covered == 0 || query_length == 0 {
        return -1;
    }
    let mut pct = 100.0 * covered as f64 / query_length as f64;
    if pct < 99.0 {
        pct += 0.5;
    }
    pct as i32
}

/// The length of the union of inclusive intervals (`GetMasterCoverage` after
/// `s_MergeRangeList`).
pub fn union_length(ranges: &mut [(u64, u64)]) -> usize {
    ranges.sort_unstable();
    let mut total = 0u64;
    let mut current: Option<(u64, u64)> = None;
    for &(from, to) in ranges.iter() {
        match current {
            Some((start, end)) if from <= end => current = Some((start, end.max(to))),
            Some((start, end)) => {
                total += end - start + 1;
                current = Some((from, to));
            }
            None => current = Some((from, to)),
        }
    }
    if let Some((start, end)) = current {
        total += end - start + 1;
    }
    total as usize
}

/// The HSP's query coverage in percent cut to an integer (`qcovhsp`): `x_SetQueryCovSeqalign`
/// reads the Seq-align's `hsp_percent_coverage` (`Blast_HSPGetQueryCoverage` of the HSP's query
/// context: `query_span` letters of a context of `query_length`).
pub fn hsp_coverage(query_span: i32, query_length: i32) -> i32 {
    crate::report::serial::seqalign::query_coverage(0, query_span, query_length) as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fields_follow_ncbis_rules() {
        assert_eq!(parse_fields(""), STD_FIELDS.to_vec());
        assert_eq!(parse_fields("bogus"), STD_FIELDS.to_vec());
        assert_eq!(
            parse_fields("qseqid  sseqid qseqid -qseqid evalue"),
            vec![SubjectSeqId, Evalue]
        );
        let mut std_less = STD_FIELDS.to_vec();
        std_less.retain(|field| *field != Evalue);
        std_less.push(QueryCovHsp);
        assert_eq!(parse_fields("std -evalue qcovhsp"), std_less);
        for (index, (name, field, _)) in FIELDS.iter().enumerate() {
            assert_eq!(*field as usize, index);
            assert_eq!(TabularField::from_name(name), Some(*field));
        }
    }

    #[test]
    fn delimiter_conflicts_are_ncbis() {
        assert!(custom_delimiter_conflict("qseqid staxids", ";")
            .unwrap()
            .contains("(staxids)"));
        assert!(custom_delimiter_conflict("std -frames", "/").is_some());
        assert!(custom_delimiter_conflict("salltitles", "<>").is_some());
        assert!(custom_delimiter_conflict("qseqid", ";").is_none());
        assert!(custom_delimiter_conflict("  ", ";").is_none());
    }

    #[test]
    fn default_reports_have_no_spec() {
        let choice = |spec: &str, delimiter: &str| FormatChoice {
            number: 6,
            spec: spec.to_string(),
            delimiter: delimiter.to_string(),
        };
        assert!(TabularSpec::of(ReportFormat::Tabular, &choice("", "")).is_none());
        assert!(TabularSpec::of(ReportFormat::Tabular, &choice("std", "\t")).is_none());
        assert!(TabularSpec::of(ReportFormat::TabularWithComments, &choice("", "")).is_none());
        let csv = TabularSpec::of(ReportFormat::CommaSeparatedValues, &choice("", "")).unwrap();
        assert_eq!((csv.delimiter.as_str(), csv.fields.len()), (",", 12));
        let piped = TabularSpec::of(ReportFormat::Tabular, &choice("std", "|")).unwrap();
        assert_eq!(piped.delimiter, "|");
    }

    #[test]
    fn coverage_merges_touching_ranges() {
        assert_eq!(union_length(&mut [(10, 19), (0, 9), (30, 39)]), 30);
        assert_eq!(union_length(&mut [(0, 9), (5, 14)]), 15);
        assert_eq!(coverage_percent(0, 100), -1);
        assert_eq!(coverage_percent(50, 100), 50);
        assert_eq!(coverage_percent(995, 1000), 99);
        assert_eq!(coverage_percent(1000, 1000), 100);
    }

    // NCBI tabular.cpp:626-633 and 693-700: without the taxonomy database's names, `staxids`
    // brings `N/A` fillers, quoted in a CSV (the singular names only for a non-zero taxid, the
    // plural ones once per taxid).
    #[test]
    fn taxonomy_names_get_ncbis_fillers() {
        let row = |format: ReportFormat, fields: &str, taxonomy: &TabularTaxonomy| {
            let choice = FormatChoice {
                number: 10,
                spec: fields.to_string(),
                delimiter: String::new(),
            };
            let spec = TabularSpec::of(format, &choice).unwrap();
            let writer = TabularWriter {
                spec: &spec,
                program_version: "BLASTN 2.17.0+",
                matrix: None,
            };
            let subject = TabularSubject {
                id: b"s",
                accession: b"s",
                accession_version: b"s",
                gi: 0,
                all_ids: &[],
                length: 10,
                titles: None,
                taxonomy,
                nucleotide: true,
                qcovs: -1,
                qcovus: -1,
            };
            let hsp = TabularHsp {
                q_start: 1,
                q_end: 10,
                s_start: 1,
                s_end: 10,
                evalue: 1e-5,
                bit_score: 20.0,
                score: 10,
                length: 10,
                gap_opens: 0,
                gaps: 0,
                num_ident: 10,
                q_frame: 1,
                s_frame: 1,
                rows: None,
                qcovhsp: 100,
            };
            let mut out = Vec::new();
            writer.row(
                &mut out,
                &TabularQuery {
                    id: b"q",
                    length: 10,
                },
                &subject,
                &hsp,
            );
            String::from_utf8(out).unwrap()
        };
        let two = TabularTaxonomy {
            taxid: 9606,
            taxids: vec![9606, 10090],
        };
        let fields =
            "ssciname scomname sblastname sskingdom sscinames scomnames sblastnames sskingdoms";
        assert_eq!(
            row(
                ReportFormat::CommaSeparatedValues,
                &format!("{fields} staxids"),
                &two
            ),
            "\"N/A\",\"N/A\",N/A,N/A,\"N/A;N/A\",\"N/A;N/A\",\"N/A\",N/A,9606;10090\n"
        );
        assert_eq!(
            row(ReportFormat::Tabular, &format!("{fields} staxids"), &two),
            "N/A\tN/A\tN/A\tN/A\tN/A;N/A\tN/A;N/A\tN/A\tN/A\t9606;10090\n"
        );
        assert_eq!(
            row(ReportFormat::CommaSeparatedValues, fields, &two),
            "N/A,N/A,N/A,N/A,\"N/A\",\"N/A\",\"N/A\",N/A\n"
        );
        let none = TabularTaxonomy {
            taxid: 0,
            taxids: vec![0],
        };
        assert_eq!(
            row(
                ReportFormat::CommaSeparatedValues,
                "ssciname sscinames staxid staxids",
                &none
            ),
            "N/A,\"N/A\",N/A,0\n"
        );
    }

    #[test]
    fn btop_counts_runs() {
        let mut out = Vec::new();
        write_btop(&mut out, b"ACGT-A", b"ACCTTA");
        assert_eq!(out, b"2GC1-T1");
    }
}
