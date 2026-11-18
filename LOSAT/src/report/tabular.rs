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
use crate::config::ScoringMatrix;
use crate::report::outfmt6::{
    format_bitscore_ncbi, format_evalue_ncbi_tabular, format_percent_identity_ncbi,
    write_outfmt7_query_line,
};
use crate::report::pairwise::end_gap_span;

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

    /// The field's bit in a set of fields (`TabularSpec::has_any`): every field is one of the
    /// `FIELDS`, by its index (`fields_follow_ncbis_rules`).
    #[inline]
    const fn bit(self) -> u64 {
        const _: () = assert!(FIELDS.len() <= u64::BITS as usize);
        1 << self as u32
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
    let add = |fields: &mut Vec<TabularField>, field| {
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
    let spec = crate::ncbi_str::truncate_spaces(spec.as_bytes());
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
    /// outfmt 18, the organism report (`report::tax_report`): no fields; it reads the final hit
    /// lists as the tabular reports do (NCBI's `PrintOneResultSet` writes it after the same
    /// checks of the results, from the same Seq-align set), and each subject's Bioseq,
    /// definition lines and the names of their taxids.
    TaxReport,
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

/// Whether `program`'s tabular report of the field list `spec` reads the sequence that NCBI's
/// formatter finds by a subject's Seq-align id ([`TabularSpec::reads_subject`]); `None` is the
/// default field list, whose `saccver` reads it.
pub fn reads_subject(spec: Option<&TabularSpec>, program: crate::blastinput::app::Program) -> bool {
    spec.is_none_or(|spec| spec.reads_subject(program))
}

/// [`TabularSpec::writes_rows`] of a tabular report of `spec`, `None` for the default fields:
/// the one rule by which the programs make (and refuse) the tabular subject ids of a report.
pub fn writes_rows(spec: Option<&TabularSpec>) -> bool {
    spec.is_none_or(TabularSpec::writes_rows)
}

/// Which rows of a tabular report read the subject's Bioseq in NCBI's formatter
/// (`TabularSpec::subject_fetch`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SubjectFetch {
    /// No row.
    Never,
    /// Every row.
    Every,
    /// The rows whose HSP has identities (the engine's, the Seq-align's `num_ident`): `nident`
    /// is the only field that enters the alignment block.
    WithIdentities,
}

impl SubjectFetch {
    /// Whether the row of an HSP with the engine's identities `num_ident` reads the subject.
    pub fn reads(self, num_ident: usize) -> bool {
        match self {
            Self::Never => false,
            Self::Every => true,
            Self::WithIdentities => num_ident > 0,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:612-660
// (`CBlastFormatUtil::InsertSubjectScores`: after a group's first Seq-align it takes the following
// ones whose subject id `Match`es that first id, and stores the coverage of the whole group on the
// first; `CBlastTabularInfo::x_SetQueryCovSubject`, tabular.cpp:707-719, writes it on the later
// rows of the same FASTA id and `N/A` on the others).
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:808-811
// ```c++
//     // Extract the full list of subject ids
//     bool setSubjectIds = (x_IsFieldRequested(eSubjectAllSeqIds) ||
//     					  x_IsFieldRequested(eSubjectAllGis) ||
//                           x_IsFieldRequested(eSubjectAllAccessions));
// ```
/// The subject ids that a tabular report prints of those that NCBI's `GetSeqIdList` makes, in
/// which the first word of the subject Bioseq's generated title can replace an id
/// (`ParsedSubject::check_printed_word`): `sseqid`, `sacc` and `saccver` of the subject (the
/// default fields print `saccver`), `sallseqid` and `sallacc` of every definition line
/// (`sallgi` prints the GIs, `salltitles` the stored titles).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PrintedIds {
    /// `sseqid`, `sacc`, `saccver`.
    pub subject: bool,
    /// `sallseqid`, `sallacc`.
    pub all: bool,
}

impl PrintedIds {
    /// The ids of the report of `spec`, `None` for the default fields; none of the organism
    /// report (`TabularSpec::writes_rows`).
    pub fn of(spec: Option<&TabularSpec>) -> Self {
        let Some(spec) = spec else {
            return Self {
                subject: true,
                all: false,
            };
        };
        if !spec.writes_rows() {
            return Self::default();
        }
        Self {
            subject: spec.has_any(&[
                TabularField::SubjectSeqId,
                TabularField::SubjectAccession,
                TabularField::SubjectAccessionVersion,
            ]),
            all: spec.has_any(&[
                TabularField::SubjectAllSeqIds,
                TabularField::SubjectAllAccessions,
            ]),
        }
    }

    /// The ids that either report prints.
    pub fn or(self, other: Self) -> Self {
        Self {
            subject: self.subject || other.subject,
            all: self.all || other.all,
        }
    }
}

impl TabularNeeds {
    /// Whether the rows count following subjects whose Seq-align ids match as one subject:
    /// `qcovs` and `qcovus` (`ParsedShown::joined`).
    pub fn joins_subjects(&self) -> bool {
        self.qcovs || self.qcovus
    }
}

impl TabularSpec {
    /// The tabular report of `format` with the custom specification of `choice`, or `None` for
    /// a format that is not tabular and for the default outfmt 6 and 7 (the `std` fields and a
    /// tab), which the programs write with their own writers. The delimiter is the custom one:
    /// the report of an application whose formatter gets it (`TabularSpec::of_app`).
    pub fn of(format: ReportFormat, choice: &FormatChoice) -> Option<Self> {
        Self::of_fields(format, &choice.spec, &choice.delimiter)
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:245-249
    // ```c++
    //         formatter.SetQueryRange(query_opts->GetRange());
    //         formatter.SetLineLength(fmt_args->GetLineLength());
    //         formatter.SetHitsSortOption(fmt_args->GetHitsSortOption());
    //         formatter.SetHspsSortOption(fmt_args->GetHspsSortOption());
    //         formatter.SetCustomDelimiter(fmt_args->GetCustomDelimiter());
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:501-505
    // ```c++
    //         formatter.SetQueryRange(query_opts->GetRange());
    //         formatter.SetLineLength(fmt_args->GetLineLength());
    //         if(UseXInclude(*fmt_args, args[kArgOutput].AsString())) {
    //         	formatter.SetBaseFile(args[kArgOutput].AsString());
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/rpsblast_app.cpp:186-189
    // ```c++
    //         formatter.SetQueryRange(query_opts->GetRange());
    //         formatter.SetLineLength(fmt_args->GetLineLength());
    //         if(UseXInclude(*fmt_args, args[kArgOutput].AsString())) {
    //         	formatter.SetBaseFile(args[kArgOutput].AsString());
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/rpstblastn_app.cpp:183-186
    // ```c++
    //         formatter.SetQueryRange(query_opts->GetRange());
    //         formatter.SetLineLength(fmt_args->GetLineLength());
    //         if(UseXInclude(*fmt_args, args[kArgOutput].AsString())) {
    //         	formatter.SetBaseFile(args[kArgOutput].AsString());
    // ```
    /// The tabular report of `format` as the formatter of `app`'s application writes it
    /// ([`TabularSpec::of`]). The applications of blastn, blastp, blastx, tblastn (PSI-TBLASTN
    /// too) and tblastx hand `-outfmt`'s custom delimiter to their formatter
    /// (`SetCustomDelimiter`, blastn_app.cpp:254, blastx_app.cpp:250, tblastn_app.cpp:266,
    /// tblastx_app.cpp:169); psiblast's, rpsblast's and rpstblastn's never do, so their reports
    /// keep the format's own delimiter, the tab or the comma (blast_format.cpp:783-785 above), and
    /// a `delim=` with the `std` fields leaves the default outfmt 6 and 7 (`None`). Their `delim=`
    /// is still checked against the fields (`custom_delimiter_conflict`). The one rule of the
    /// command line and the library, for every program.
    pub fn of_app(
        format: ReportFormat,
        choice: &FormatChoice,
        app: crate::report::archive::options::RequestKind,
    ) -> Option<Self> {
        use crate::report::archive::options::RequestKind;
        let delimiter = match app {
            RequestKind::Psiblast | RequestKind::Rpsblast | RequestKind::Rpstblastn => "",
            RequestKind::Blastn
            | RequestKind::Blastp
            | RequestKind::Blastx
            | RequestKind::Tblastn
            | RequestKind::PsiTblastn
            | RequestKind::Tblastx => &choice.delimiter,
        };
        Self::of_fields(format, &choice.spec, delimiter)
    }

    /// [`TabularSpec::of`] with the delimiter `delimiter` (empty: the format's own).
    fn of_fields(format: ReportFormat, spec: &str, delimiter: &str) -> Option<Self> {
        let (kind, default_delimiter) = match format {
            ReportFormat::Tabular => (TabularKind::Tabular, "\t"),
            ReportFormat::TabularWithComments => (TabularKind::Commented, "\t"),
            ReportFormat::CommaSeparatedValues => (TabularKind::Csv, ","),
            ReportFormat::CommaSeparatedValuesWithHeader => (TabularKind::CsvWithHeader, ","),
            // NCBI's `ParseFormattingString` keeps no field list for outfmt 18.
            ReportFormat::TaxFormat => {
                return Some(TabularSpec {
                    kind: TabularKind::TaxReport,
                    fields: Vec::new(),
                    delimiter: String::new(),
                })
            }
            _ => return None,
        };
        let fields = parse_fields(spec);
        let delimiter = if delimiter.is_empty() {
            default_delimiter
        } else {
            delimiter
        }
        .to_string();
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

    /// Whether any of `fields` is among the spec's fields: one pass over the spec's fields, each
    /// tested against the set of `fields` (`TabularField::bit`), instead of one pass for each of
    /// `fields` (no NCBI counterpart: code placement only).
    #[inline]
    pub(crate) fn has_any(&self, fields: &[TabularField]) -> bool {
        let set = fields.iter().fold(0, |set, field| set | field.bit());
        self.fields.iter().any(|field| set & field.bit() != 0)
    }

    /// Whether NCBI's formatter reads the query's Bioseq for a row, and with it the HSP's query
    /// coverage (`x_SetQueryCovSeqalign`, which reads the Seq-align's query range where the
    /// Seq-align has no `hsp_percent_coverage`): one of the query's id fields or `qcovhsp`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:778-800
    /// ```c++
    ///     if (x_IsFieldRequested(eQuerySeqId) || x_IsFieldRequested(eQueryGi) ||
    ///         x_IsFieldRequested(eQueryAccession) ||
    ///         x_IsFieldRequested(eQueryAccessionVersion) ||
    ///         x_IsFieldRequested(eQueryCovSeqalign)) {
    ///         try {
    ///             ...
    ///             x_SetQueryCovSeqalign(align, m_QueryLength);
    ///         } catch (const CException&) {
    ///             list<CRef<CSeq_id> > query_ids;
    ///             CRef<CSeq_id> id(new CSeq_id());
    ///             id->Assign(align.GetSeq_id(0));
    ///             query_ids.push_back(id);
    ///             SetQueryId(query_ids);
    ///             bioseqs_found = false;
    ///         }
    ///     }
    /// ```
    pub fn reads_query(&self) -> bool {
        self.has_any(&[
            QuerySeqId,
            QueryGi,
            QueryAccession,
            QueryAccessionVersion,
            QueryCovHsp,
        ])
    }

    /// Whether NCBI's formatter reads the subject's Bioseq for a row before its alignment block,
    /// which sets the subject's length (`slen`) there: one of the subject's id, taxonomy, title or
    /// strand fields. A row that leaves before the alignment block (`bioseqs_found` false) shows
    /// an `slen` of 0 without one of them.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:833-842
    /// ```c++
    ///     if(setSubjectIds || setSubjectTaxInfo || setSubjectTaxInfoAll || setSubjectTitle ||
    ///        x_IsFieldRequested(eSubjectStrand) || setSubjectId)
    ///     {
    ///         try {
    ///        		const CBioseq_Handle& subject_bh =
    ///                 scope.GetBioseqHandle(align.GetSeq_id(1));
    ///             if(setSubjectId) {
    ///             	SetSubjectId(subject_bh);
    ///             }
    ///             m_SubjectLength = subject_bh.GetBioseqLength();
    /// ```
    pub fn reads_subject_bioseq(&self) -> bool {
        self.has_any(&[
            SubjectAllSeqIds,
            SubjectAllGis,
            SubjectAllAccessions,
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
            SubjectSeqId,
            SubjectGi,
            SubjectAccession,
            SubjectAccessionVersion,
        ])
    }

    /// Whether NCBI's formatter opens its alignment block for a row whose Seq-align counts
    /// `num_ident` identities (the engine's count): one of its fields, or `nident` with
    /// identities. Only the block counts the identities on the aligned sequences and sets the
    /// frames (`SetCounts`); outside it `nident`, `frames` and `sframe` keep the 0 of
    /// `x_ResetFields`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-908
    /// ```c++
    ///     if (x_IsFieldRequested(eQueryStart) || x_IsFieldRequested(eQueryEnd) ||
    ///         x_IsFieldRequested(eSubjectStart) || x_IsFieldRequested(eSubjectEnd) ||
    ///         x_IsFieldRequested(eAlignmentLength) || x_IsFieldRequested(eGaps) ||
    ///         x_IsFieldRequested(eGapOpenings) || x_IsFieldRequested(eQuerySeq) ||
    ///         x_IsFieldRequested(eSubjectSeq) || (x_IsFieldRequested(eNumIdentical) && num_ident > 0 ) ||
    ///         x_IsFieldRequested(ePositives) || x_IsFieldRequested(eMismatches) ||
    ///         x_IsFieldRequested(ePercentPositives) || x_IsFieldRequested(ePercentIdentical) ||
    ///         x_IsFieldRequested(eQueryFrame) || x_IsFieldRequested(eBTOP) ||
    ///         x_IsFieldRequested(eSubjectStrand)) {
    /// ```
    /// Inline, and the count by reference so that it is read only where the condition reaches it:
    /// `TabularWriter::row` evaluates it for every row as the expression it replaced (no NCBI
    /// counterpart: code placement only).
    #[inline]
    pub fn opens_counts(&self, num_ident: &usize) -> bool {
        self.has_any(&[
            QueryStart,
            QueryEnd,
            SubjectStart,
            SubjectEnd,
            AlignmentLength,
            Gaps,
            GapOpenings,
            QuerySeq,
            SubjectSeq,
            Positives,
            Mismatches,
            PercentPositives,
            PercentIdentical,
            QueryFrame,
            Btop,
            SubjectStrand,
        ]) || (self.has(NumIdentical) && *num_ident > 0)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:833-837 and 973-980
    // ```c++
    //     if(setSubjectIds || setSubjectTaxInfo || setSubjectTaxInfoAll || setSubjectTitle ||
    //        x_IsFieldRequested(eSubjectStrand) || setSubjectId)
    //     {
    //         try {
    //        		const CBioseq_Handle& subject_bh =
    //     ...
    //     if (x_IsFieldRequested(eQuerySeq) ||
    //         x_IsFieldRequested(eSubjectSeq) ||
    //         x_IsFieldRequested(ePositives) ||
    //         x_IsFieldRequested(ePercentPositives) ||
    //         x_IsFieldRequested(eBTOP) ||
    //         (x_IsFieldRequested(eNumIdentical) && !kNoFetchSequence) ||
    //         (x_IsFieldRequested(eMismatches) && !kNoFetchSequence) ||
    //         (x_IsFieldRequested(ePercentIdentical) && !kNoFetchSequence)) {
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-908,915,1078-1082
    // ```c++
    //     if (x_IsFieldRequested(eQueryStart) || x_IsFieldRequested(eQueryEnd) ||
    //         x_IsFieldRequested(eSubjectStart) || x_IsFieldRequested(eSubjectEnd) ||
    //         x_IsFieldRequested(eAlignmentLength) || x_IsFieldRequested(eGaps) ||
    //         x_IsFieldRequested(eGapOpenings) || x_IsFieldRequested(eQuerySeq) ||
    //         x_IsFieldRequested(eSubjectSeq) || (x_IsFieldRequested(eNumIdentical) && num_ident > 0 ) ||
    //     ...
    //         x_IsFieldRequested(eQueryFrame) || x_IsFieldRequested(eBTOP) ||
    //         x_IsFieldRequested(eSubjectStrand)) {
    //     ...
    //     const bool kTranslated = align.GetSegs().IsStd();
    //     ...
    //         if (x_IsFieldRequested(eSubjFrame) || x_IsFieldRequested(eFrames)) {
    //             subject_frame = CAlignFormatUtil::
    //                 GetFrame (s_start - 1, ds.GetSeqStrand(kSubjectRow),
    //                           scope.GetBioseqHandle(align.GetSeq_id(1)));
    // ```
    /// Which of `program`'s rows read the subject's Bioseq in NCBI's formatter, which loads its
    /// taxonomy where no prefetch loaded it before (`PRE_FETCH_SEQS_LIMIT`): every row for the
    /// subject's ids, taxonomy, titles and strand, and the alignment's letters, also for the
    /// identities unless the formatter does not fetch sequences (BLASTN's,
    /// blast_format.cpp:791-792); in BLASTX, TBLASTN and TBLASTX, whose alignments are Std-segs,
    /// also the subject's frame (`sframe`, `frames`) where a field that enters that block is asked
    /// for (a coordinate, `length`, `gaps`, `gapopen` or `qframe`; the block's other fields fetch
    /// on their own). `slen` is read without the Bioseq (`sequence::GetLength`). NCBI enters the
    /// block for `nident` alone only on a row whose HSP has identities (the Seq-align's
    /// `num_ident` above 0): where no other field enters it, only those rows read the subject,
    /// for the letters of `nident` and the frames ([`SubjectFetch::WithIdentities`]).
    pub fn subject_fetch(&self, program: crate::blastinput::app::Program) -> SubjectFetch {
        use crate::blastinput::app::Program;
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:996-1000
        // ```c++
        //         try {
        //             const CBioseq_Handle& handle = m_Scope.GetBioseqHandle(*subid);
        //             if( !handle ) continue;
        //
        //             const CRef<CBlast_def_line_set> bdlRef =  CSeqDB::ExtractBlastDefline(handle);
        // ```
        // The organism report reads the Bioseq of each HSP that it does not skip.
        if self.kind == TabularKind::TaxReport {
            return SubjectFetch::Every;
        }
        let no_fetch = program == Program::Blastn;
        let translated = matches!(
            program,
            Program::Blastx | Program::Tblastn | Program::Tblastx
        );
        // The fields that enter the block on every row, and whether the block reads the subject
        // through `nident` or the frames.
        let block = self.has_any(&[
            QueryStart,
            QueryEnd,
            SubjectStart,
            SubjectEnd,
            AlignmentLength,
            Gaps,
            GapOpenings,
            QuerySeq,
            SubjectSeq,
            Positives,
            Mismatches,
            PercentPositives,
            PercentIdentical,
            QueryFrame,
            Btop,
            SubjectStrand,
        ]);
        let block_reads = (!no_fetch && self.has(NumIdentical))
            || (translated && self.has_any(&[SubjectFrame, Frames]));
        let every = self.has_any(&[
            SubjectSeqId,
            SubjectGi,
            SubjectAccession,
            SubjectAccessionVersion,
            SubjectAllSeqIds,
            SubjectAllGis,
            SubjectAllAccessions,
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
            QuerySeq,
            SubjectSeq,
            Positives,
            PercentPositives,
            Btop,
        ]) || (!no_fetch && self.has_any(&[Mismatches, PercentIdentical]))
            || (block && block_reads);
        if every {
            SubjectFetch::Every
        } else if block_reads && self.has(NumIdentical) {
            SubjectFetch::WithIdentities
        } else {
            SubjectFetch::Never
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:891-898
    // ```c++
    //     if (x_IsFieldRequested(eSubjectLength) && m_SubjectLength == 0) {
    //         ...
    //         m_SubjectLength = sequence::GetLength(align.GetSeq_id(1), &scope);
    //     }
    // ```
    /// Whether `program`'s rows read the sequence that NCBI's formatter finds by the subject's
    /// Seq-align id: its Bioseq or letters on some row ([`TabularSpec::subject_fetch`]), or its
    /// length (`slen`, which the scope looks up by that id without loading the Bioseq). A table
    /// of other fields (`6 qseqid evalue`, the coordinates, `length`, `gaps`, `qcovs`; in
    /// BLASTN also `pident`, `nident` and `mismatch`) writes nothing of that sequence, so it is
    /// the same whatever the lookup finds; only the batches' prefetch reads it
    /// (`common::database::ParsedShown::prefetched`).
    pub fn reads_subject(&self, program: crate::blastinput::app::Program) -> bool {
        self.subject_fetch(program) != SubjectFetch::Never || self.has(SubjectLength)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:833-838, 862-875
    // ```c++
    //     if(setSubjectIds || setSubjectTaxInfo || setSubjectTaxInfoAll || setSubjectTitle ||
    //        x_IsFieldRequested(eSubjectStrand) || setSubjectId)
    //     {
    //         try {
    //        		const CBioseq_Handle& subject_bh =
    //                 scope.GetBioseqHandle(align.GetSeq_id(1));
    //     ...
    //         } catch (const CException&) {
    //             list<CRef<CSeq_id> > subject_ids;
    //             CRef<CSeq_id> id(new CSeq_id());
    //             id->Assign(align.GetSeq_id(1));
    //             subject_ids.push_back(id);
    //             SetSubjectId(subject_ids);
    //             bioseqs_found = false;
    //         }
    //
    //     }
    //     ...
    //     if (!bioseqs_found)
    //         return -1;
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203 (the `CATCH_ALL` of
    // every application, quoted at `cli::NativeError::exception`)
    // Where a batch's prefetch is off, the formatter loads each row's subject by its Seq-align id
    // and a lookup that throws on a file throws at each such load (probe of line DBID pfoff-2,
    // `R/pfoff/p2`; DECISIONS 62 (12)).
    /// What `program`'s row of an HSP with the engine's identities `num_ident` (the Seq-align's
    /// `num_ident`) does where its subject's lookup throws (`report::db_subject::AfterThrow`): a
    /// field of the subject block (`reads_subject_bioseq`) loads the subject inside its `try`,
    /// whose `catch` ends the row's fields, and the row is written (`WritesRow`,
    /// `TabularWriter::row`); a row that reads the subject only after that block (`slen` by
    /// `sequence::GetLength`, the letters and frames of `subject_fetch`) throws outside any
    /// `catch`, and the application's `CATCH_ALL` ends the run there (`EndsRun`); `None`: the row
    /// reads nothing of the subject. The one place that decides it per row.
    pub fn row_throw(
        &self,
        program: crate::blastinput::app::Program,
        num_ident: usize,
    ) -> Option<crate::report::db_subject::AfterThrow> {
        use crate::report::db_subject::AfterThrow;
        if self.reads_subject_bioseq() {
            Some(AfterThrow::WritesRow)
        } else if self.throw_rows(program).reads(num_ident) {
            Some(AfterThrow::EndsRun)
        } else {
            None
        }
    }

    /// Which of `program`'s rows load the subject where its batch's prefetch is off, so that a
    /// lookup that throws throws there (`row_throw`): every row of a field of the subject block,
    /// of `slen` (`sequence::GetLength`, on every row) or of a field that reads the subject on
    /// every row (`subject_fetch`); for `nident` alone among the fields that read the subject, the
    /// rows of HSPs with identities (`SubjectFetch::WithIdentities`): a row of no identities reads
    /// nothing of it and prints the 0 of `x_ResetFields`.
    pub fn throw_rows(&self, program: crate::blastinput::app::Program) -> SubjectFetch {
        if self.reads_subject_bioseq() || self.has(SubjectLength) {
            SubjectFetch::Every
        } else {
            self.subject_fetch(program)
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:995-1023 (quoted at
    // `report::tax_report::organism_report`: the `catch` of each HSP's load, `continue`)
    /// What the report does with a subject whose lookup throws where its batch's prefetch is off
    /// (`report::db_subject::FormatThrow`, `row_throw` of its rows, at its `throw_rows`): a table
    /// (outfmt 6, 7, 10 and 20) as its fields say, a table that reads nothing of the subject never
    /// looks it up; the organism report (outfmt 18) leaves the HSPs of such a subject out and
    /// goes on (`WritesRow`: each HSP looks its own Seq-align id up).
    pub fn format_throw(
        &self,
        program: crate::blastinput::app::Program,
    ) -> crate::report::db_subject::FormatThrow {
        use crate::report::db_subject::{AfterThrow, FormatThrow, RunThrow};
        // The organism report's HSPs, and a table's rows with a field of the subject block, catch
        // the throw and go on.
        if self.kind == TabularKind::TaxReport || self.reads_subject_bioseq() {
            FormatThrow::Writes(RunThrow::every_row(AfterThrow::WritesRow))
        } else if self.reads_subject(program) {
            FormatThrow::Writes(RunThrow {
                after: AfterThrow::EndsRun,
                rows: self.throw_rows(program),
            })
        } else {
            FormatThrow::Unread
        }
    }

    pub fn needs(&self) -> TabularNeeds {
        // The organism report shows the names of each definition line's taxid.
        if self.kind == TabularKind::TaxReport {
            return TabularNeeds {
                taxonomy: true,
                taxonomy_names: true,
                ..TabularNeeds::default()
            };
        }
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

    /// Whether the report writes rows of fields: every tabular report but the organism report
    /// (outfmt 18), whose table shows the labels and titles of the subjects' definition lines
    /// (`report::tax_report`, taxFormat.cpp:935-976) and never a row's subject ids (`sseqid`,
    /// `saccver`: `GetSeqIdList` and `s_ReplaceLocalId`, tabular.cpp:522-532), so that the
    /// refusals of those ids (`PrintedIds`, the programs' tabular subject ids) do not apply to it.
    /// The Bioseq reads (`subject_fetch`) and the taxonomy (`needs`) apply to both.
    pub fn writes_rows(&self) -> bool {
        self.kind != TabularKind::TaxReport
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:808-811,828-851
    // ```c++
    //     bool setSubjectIds = (x_IsFieldRequested(eSubjectAllSeqIds) ||
    //     					  x_IsFieldRequested(eSubjectAllGis) ||
    //                           x_IsFieldRequested(eSubjectAllAccessions));
    //     ...
    //     bool setSubjectId = (x_IsFieldRequested(eSubjectSeqId) ||
    //     		 	 	 	 x_IsFieldRequested(eSubjectGi) ||
    //     		 	 	 	 x_IsFieldRequested(eSubjectAccession) ||
    //     		 	 	 	 x_IsFieldRequested(eSubjAccessionVersion));
    //     ...
    //             if(setSubjectId) {
    //             	SetSubjectId(subject_bh);
    //             }
    //             m_SubjectLength = subject_bh.GetBioseqLength();
    //
    //             if(setSubjectIds || setSubjectTaxInfo || setSubjectTitle || setSubjectTaxInfoAll) {
    //             	CRef<CBlast_def_line_set> bdlRef =
    //             			CSeqDB::ExtractBlastDefline(subject_bh);
    //             	if(setSubjectIds) {
    //             		x_SetSubjectIds(subject_bh, bdlRef);
    //             	}
    // ```
    /// Where `SetFields` stops for a row of a subject with the lookup `lookup`: `None` where the
    /// lookup does not fail or no field of the list makes it (the organism report writes no rows).
    pub fn failed_stage(&self, lookup: SubjectLookup<'_>) -> Option<FailedStage> {
        match lookup {
            SubjectLookup::Found => None,
            SubjectLookup::Failed(failed) => self.failure_stage(failed.failure()),
            SubjectLookup::Throws(_) => self.failure_stage(LookupFailure::thrown()),
        }
    }

    /// `failed_stage` of a subject whose lookups fail as `failure` says. The load of the subject
    /// (`failure.bioseq`) comes before the other lookups, in the subject block.
    pub fn failure_stage(&self, failure: LookupFailure) -> Option<FailedStage> {
        if !self.writes_rows() {
            None
        } else if failure.bioseq {
            self.reads_subject_bioseq().then_some(FailedStage::Bioseq)
        } else if failure.subject_id
            && self.has_any(&[
                SubjectSeqId,
                SubjectGi,
                SubjectAccession,
                SubjectAccessionVersion,
            ])
        {
            Some(FailedStage::SubjectId)
        } else if failure.all_ids
            && self.has_any(&[SubjectAllSeqIds, SubjectAllGis, SubjectAllAccessions])
        {
            Some(FailedStage::SubjectIds)
        } else {
            None
        }
    }

    /// Whether a row of this report can stop where NCBI has no result (`failed_lookup_stop`):
    /// the report asks `sacc` or `saccver`, or `staxid` (a failed first row: `staxid` loads the
    /// subject, whose load may throw, `FailedStage::Bioseq`). A report that cannot stop needs no
    /// check of its rows before it is written.
    pub fn may_stop_at_failed_lookup(&self) -> bool {
        self.writes_rows()
            && self.has_any(&[SubjectAccession, SubjectAccessionVersion, SubjectTaxId])
    }

    /// Whether the report looks every taxid of its subjects up in the taxonomy database, 0 too:
    /// the plural names (`x_SetTaxInfoAll`, `checks_taxdb`) and the organism report (outfmt 18,
    /// `CTaxFormat::x_InitBlastDBTaxInfo`); its rows write the texts of those lookups.
    pub fn looks_up_every_taxid(&self) -> bool {
        self.checks_taxdb() || self.kind == TabularKind::TaxReport
    }

    /// Whether a singular name is among the fields (`x_SetTaxInfo` looks `staxid` up).
    pub fn has_singular_names(&self) -> bool {
        self.has_any(&[
            SubjectSciName,
            SubjectCommonName,
            SubjectBlastName,
            SubjectSuperKingdom,
        ])
    }

    /// `x_IsCsv`: the names are quoted when the delimiter is a comma.
    fn is_csv(&self) -> bool {
        self.delimiter == ","
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:174-185
// ```c++
//     case CBlastTabularInfo::eAccession:
//     {
//         CConstRef<CSeq_id> accid = FindBestChoice(id, CSeq_id::Score);
//         accid->GetLabel(&id_str, CSeq_id::eContent, 0);
//         break;
//     }
//     case CBlastTabularInfo::eAccVersion:
//     {
//         CConstRef<CSeq_id> accid = FindBestChoice(id, CSeq_id::Score);
//         accid->GetLabel(&id_str, CSeq_id::eContent, CSeq_id::fLabel_Version);
//         break;
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:262-268
// ```c++
// void CBlastTabularInfo::x_PrintSubjectTaxId()
// {
// 	if(m_SubjectTaxId == ZERO_TAX_ID) {
// 		m_Ostream << NA;
// 		return;
// 	}
// 	m_Ostream << m_SubjectTaxId;
// ```
/// Why NCBI has no result for a tabular report of `spec` (`None`: the default fields, whose
/// `saccver` is the first case) with a row whose subject lookups fail as `failure` says, the row
/// being the first of its result set when `first`: `sacc` and `saccver` of the empty
/// `m_SubjectId` (`SetSubjectId` threw) take the label of a null `FindBestChoice` (NCBI BLAST+
/// 2.17.0 stops with `CCoreException::eNullPtr`, exit 255, after the row's earlier fields), and
/// `staxid` of a first row prints `m_SubjectTaxId` before any row set it (an uninitialised `int`:
/// a different number in each run). `None`: the report has a result (`TabularWriter::row`).
pub fn failed_lookup_stop(
    spec: Option<&TabularSpec>,
    failure: LookupFailure,
    first: bool,
) -> Option<&'static str> {
    let Some(spec) = spec else {
        return (failure.subject_id || (failure.bioseq && first))
            .then_some(FAILED_LOOKUP_NULL_POINTER);
    };
    let stage = spec.failure_stage(failure)?;
    // A load that throws leaves `m_SubjectId` as the last row that found its subject set it: empty
    // only before such a row (the result set's first row; every row before a later one that
    // throws has thrown too, which this check meets first).
    let empty_id = match stage {
        FailedStage::SubjectId => true,
        FailedStage::Bioseq => first,
        FailedStage::SubjectIds => false,
    };
    if empty_id && spec.has_any(&[SubjectAccession, SubjectAccessionVersion]) {
        Some(FAILED_LOOKUP_NULL_POINTER)
    } else if first && spec.has(SubjectTaxId) {
        Some("NCBI BLAST+ prints an uninitialised staxid in the row when it is the first row of a query's report (a different number in each run)")
    } else {
        None
    }
}

const FAILED_LOOKUP_NULL_POINTER: &str = "NCBI BLAST+ stops with a null pointer exception (exit 255) at the sacc or saccver of its row (the default fields of outfmt 6, 7, 10 and 20 print saccver)";

/// Whether a tabular report of `spec` (`None`: the default fields) can stop where NCBI has no
/// result (`TabularSpec::may_stop_at_failed_lookup`), so that its rows are checked
/// (`check_failed_lookups`) before it is written.
pub fn may_stop_at_failed_lookup(spec: Option<&TabularSpec>) -> bool {
    spec.is_none_or(TabularSpec::may_stop_at_failed_lookup)
}

/// The check, before a tabular report of `spec` (`None`: the default fields) writes the rows of
/// `result_sets` (each the indices of its subjects in row order, a subject's run of rows once,
/// one set per query or psiblast round), that NCBI has a result for each row whose subject
/// lookups fail (`failure`, `None` where they do not): `failed_lookup_stop`'s cases are LOSAT's
/// refusals, which name the subject (`describe`: what its title is) and the program.
pub fn check_failed_lookups<S>(
    spec: Option<&TabularSpec>,
    result_sets: impl IntoIterator<Item = S>,
    failure: impl Fn(usize) -> Option<LookupFailure>,
    describe: impl Fn(usize) -> String,
    program: &str,
) -> anyhow::Result<()>
where
    S: IntoIterator<Item = usize>,
{
    if !may_stop_at_failed_lookup(spec) {
        return Ok(());
    }
    for set in result_sets {
        for (row, s_idx) in set.into_iter().enumerate() {
            let Some(failure) = failure(s_idx) else {
                continue;
            };
            if let Some(reason) = failed_lookup_stop(spec, failure, row == 0) {
                return Err(FailedLookupRefusal(format!(
                    "{}; {reason}, which is not supported by LOSAT's {program}",
                    describe(s_idx)
                ))
                .into());
            }
        }
    }
    Ok(())
}

/// LOSAT's refusal of `check_failed_lookups` (exit status 1 with its message in every program; a
/// type of its own, so that a program that reports its engine's errors otherwise, BLASTX's
/// `native::report_error`, passes it on as a refusal).
#[derive(Debug)]
pub struct FailedLookupRefusal(String);

impl std::fmt::Display for FailedLookupRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for FailedLookupRefusal {}

/// What `check_failed_lookups` says of the `-subject` record `s_idx` (0-based) shown as
/// `shown_id`, whose lookup fails (a `Subject_` first word of its title, or a parsed
/// `lcl|Subject_N` id).
pub fn failed_record_title(s_idx: usize, shown_id: &[u8]) -> String {
    format!(
        "subject record {} ({}) has a title whose non-UTF-8 bytes are in no encoding that NCBI BLAST+ recognises, from which NCBI BLAST+ makes the id of a 'Subject_' record",
        s_idx + 1,
        String::from_utf8_lossy(shown_id)
    )
}

/// What `check_failed_lookups` says of the database sequence `oid` whose lookup fails: where its
/// load throws `thrown` (a batch whose prefetch is off, `SubjectLookup::Throws`), that load;
/// otherwise its title (`failed_database_title`).
pub fn failed_database_subject(
    oid: impl std::fmt::Display,
    thrown: Option<&crate::ncbi_exception::NcbiException>,
) -> String {
    match thrown {
        Some(exception) => format!(
            "database sequence {oid} has a Seq-align id whose lookup throws '{}' in a batch whose prefetch PRE_FETCH_SEQS_LIMIT turns off, where NCBI BLAST+ 2.17.0 reads the reported sequences one by one",
            exception.get_msg()
        ),
        None => failed_database_title(oid),
    }
}

/// What `check_failed_lookups` says of the database sequence `oid` whose lookup fails.
pub fn failed_database_title(oid: impl std::fmt::Display) -> String {
    format!(
        "database sequence {oid} has a title (or a name in the taxonomy database from which NCBI BLAST+ makes an empty one or the end of a protein title ending in ']') in no encoding that NCBI BLAST+ recognises"
    )
}

/// The warning of `x_CheckTaxDB` (`TabularSpec::checks_taxdb`), as NCBI's `ERR_POST` writes it
/// for `program` (lower case).
/// The post's header alone after NCBI's main-thread diagnostic stream failed
/// (`ncbi_exception::post_line`; line RG D3-1 `om/mute/` e09).
pub fn taxdb_warning_line(program: &str) -> String {
    String::from_utf8(crate::ncbi_exception::post_line(
        crate::ncbi_exception::Severity::Warning,
        None,
        Some(program),
        b"Taxonomy name lookup from taxid requires installation of taxdb database with ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz",
    ))
    .expect("ASCII")
}

/// A query of the report.
#[derive(Clone, Copy, Debug)]
pub struct TabularQuery<'a> {
    /// `qseqid`: a query with a local id shows the first word of the title
    /// (`FastaRecord::shown_id`) in `id`, `accession` and `accession_version`, and `gi` is 0; with
    /// `-parse_deflines` these are the record's ids (`TabularQuery::of_record`).
    pub id: &'a [u8],
    /// `qacc`.
    pub accession: &'a [u8],
    /// `qaccver`.
    pub accession_version: &'a [u8],
    /// `qgi`.
    pub gi: i64,
    /// `qlen`: the searched length (of the `-query_loc` interval).
    pub length: usize,
}

impl<'a> TabularQuery<'a> {
    /// The query with the shown id `id` (every id field) and the searched length `length`.
    pub fn local(id: &'a [u8], length: usize) -> Self {
        Self {
            id,
            accession: id,
            accession_version: id,
            gi: 0,
            length,
        }
    }

    /// The query `query`, shown as `shown_id` when `-parse_deflines` did not read its ids.
    pub fn of_ids(
        believed: Option<&'a crate::blastinput::fasta_reader::believed::BelievedIds>,
        shown_id: &'a [u8],
        length: usize,
    ) -> Self {
        match believed {
            Some(ids) => Self {
                id: &ids.tabular[0],
                accession: &ids.tabular[1],
                accession_version: &ids.tabular[2],
                gi: ids.gi,
                length,
            },
            None => Self::local(shown_id, length),
        }
    }
}

/// The ids of one definition line of a subject (`sallseqid`, `sallgi`, `sallacc`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DeflineIds {
    pub id: Vec<u8>,
    pub gi: i64,
    pub accession: Vec<u8>,
}

impl DeflineIds {
    /// The one list of `x_SetSubjectIds` of a `-subject` record (no definition lines): its local
    /// id after `s_ReplaceLocalId`, the first word of its title (an empty title: the local id;
    /// `FastaRecord::shown_id`), as `word`. `GetSeqIdList`, which makes `sseqid`, does not run
    /// there, so a title that `GenerateDefline` changes (`Subject_4.`) keeps its first word.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:557-564
    /// ```c++
    ///     } else {
    ///         // Blast-def-line is not filled, hence retrieve all Seq-ids directly
    ///         // from the Bioseq handle's Seq-id.
    ///         list<CRef<objects::CSeq_id> > subject_id_list;
    ///         ITERATE(CBioseq_Handle::TId, itr, bh.GetId()) {
    ///             CRef<CSeq_id> next_id = s_ReplaceLocalId(bh, itr->GetSeqId(), m_ParseLocalIds);
    ///             subject_id_list.push_back(next_id);
    ///         }
    ///         m_SubjectIds.push_back(subject_id_list);
    /// ```
    pub fn of_record(word: &[u8]) -> Vec<DeflineIds> {
        vec![DeflineIds {
            id: word.to_vec(),
            gi: 0,
            accession: word.to_vec(),
        }]
    }
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

/// The names of the taxonomy database that a report's rows look up (`x_SetTaxInfoAll`, then
/// `x_SetTaxInfo`: `CSeqDB::GetTaxInfo` of each taxid of `staxids` for the plural names, of
/// `staxid` for the singular ones), asked once per taxid before the report's first byte, so that
/// a lookup the reader refuses comes before it, with the texts that NCBI writes to standard error
/// at each lookup. `None` in a report (`TabularSubject::names`) is a run without `taxdb.bti` on
/// the search path, or one that looks nothing up: every lookup fails, as NCBI's does without the
/// database.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:669-695
/// ```c++
///             for(; itr !=  m_SubjectTaxIds.end(); ++itr) {
///             	try {
///                     SSeqDBTaxInfo taxinfo;
///                     CSeqDB::GetTaxInfo(*itr, taxinfo);
///                     m_SubjectSciNames.push_back(taxinfo.scientific_name);
///                     m_SubjectCommonNames.push_back(taxinfo.common_name);
///                     if(s_IsValidName(taxinfo.blast_name))
///                             m_SubjectBlastNames.insert(taxinfo.blast_name);
///
///                     if(s_IsValidName(taxinfo.s_kingdom))
///                             m_SubjectSuperKingdoms.insert(taxinfo.s_kingdom);
///
///             	} catch (const CException&) {
/// ```
#[derive(Debug)]
pub struct TabularNames {
    /// The program's name in lower case, as NCBI's messages show it.
    program: String,
    /// The texts of opening the database, when the lookups of the report opened it: NCBI writes
    /// them at the first lookup of the process.
    open: Vec<String>,
    opened: std::cell::Cell<bool>,
    /// What `GetTaxInfo` returned for each taxid (`None` when it threw), and the texts of its
    /// lookup, which every lookup of the taxid writes.
    answers: std::collections::HashMap<i32, (Option<crate::seqdb::taxdb::TaxNames>, Vec<String>)>,
    /// The organism report's (outfmt 18) database, which it asks itself (`for_report`).
    report: Option<crate::report::tax_report::TaxNamesDb>,
}

impl TabularNames {
    /// The lookups of the rows of the subjects `taxonomies` (in any order, repeats allowed) for
    /// the fields of `spec`, through `taxdb`; `program` names the messages' program. An open or a
    /// lookup that the reader refuses (`SeqDbError::Unsupported`) is the error, and so is a super
    /// kingdom that a field of `spec` writes and that NCBI reads past its mapping of `taxdb.btd`
    /// (`TaxNames::super_kingdom`).
    pub fn look_up<'t>(
        spec: &TabularSpec,
        taxdb: &crate::seqdb::taxdb::LazyTaxDb,
        taxonomies: impl IntoIterator<Item = &'t TabularTaxonomy>,
        program: &str,
    ) -> Result<Self, crate::seqdb::SeqDbError> {
        let mut names = TabularNames {
            program: program.to_ascii_lowercase(),
            open: Vec::new(),
            opened: std::cell::Cell::new(false),
            answers: Default::default(),
            report: None,
        };
        let (plural, singular) = (spec.checks_taxdb(), spec.has_singular_names());
        // Each lookup, and whether a field writes its super kingdom.
        let (plural_kingdom, singular_kingdom) = (
            spec.has(SubjectSuperKingdoms),
            spec.has(SubjectSuperKingdom),
        );
        let mut asked: Vec<(i32, bool)> = Vec::new();
        for taxonomy in taxonomies {
            if plural {
                asked.extend(taxonomy.taxids.iter().map(|&taxid| (taxid, plural_kingdom)));
            }
            if singular && taxonomy.taxid != 0 {
                asked.push((taxonomy.taxid, singular_kingdom));
            }
        }
        if asked.is_empty() {
            return Ok(names);
        }
        names.open = taxdb.open()?;
        for (taxid, kingdom) in asked {
            let (found, _) = match names.answers.entry(taxid) {
                std::collections::hash_map::Entry::Occupied(entry) => entry.into_mut(),
                std::collections::hash_map::Entry::Vacant(entry) => {
                    let mut messages = Vec::new();
                    let found = taxdb.names(taxid, &mut messages)?;
                    entry.insert((found, messages))
                }
            };
            if let (true, Some(Err(past))) = (kingdom, found.as_ref().map(|f| &f.super_kingdom)) {
                return Err(past.error("super kingdom"));
            }
        }
        Ok(names)
    }

    /// Whether a row of the report looks a name up (a taxid was asked for): only then are the
    /// texts of opening written (`row_messages`).
    pub fn looks_up(&self) -> bool {
        !self.answers.is_empty()
    }

    /// These lookups with `texts`, the texts of an open of the database by the rows before them
    /// (an earlier query batch's), written before the texts of their own open.
    pub fn opened_before(mut self, mut texts: Vec<String>) -> Self {
        texts.append(&mut self.open);
        self.open = texts;
        self
    }

    /// What the lookup of `taxid` found (`None`: not found, or not looked up).
    fn get(&self, taxid: i32) -> Option<&crate::seqdb::taxdb::TaxNames> {
        self.answers
            .get(&taxid)
            .and_then(|(found, _)| found.as_ref())
    }

    /// The names of the organism report (outfmt 18), which looks each taxid up itself where NCBI's
    /// report does (`report::tax_report::TaxLookups`) in the run's database `db`.
    pub fn for_report(db: crate::report::tax_report::TaxNamesDb, program: &str) -> Self {
        TabularNames {
            program: program.to_ascii_lowercase(),
            open: Vec::new(),
            opened: std::cell::Cell::new(false),
            answers: Default::default(),
            report: Some(db),
        }
    }

    /// The organism report's database (`for_report`).
    pub fn report_db(&self) -> Option<&crate::report::tax_report::TaxNamesDb> {
        self.report.as_ref()
    }

    /// The standard-error lines of the lookups of one row of a subject of `taxonomy`, in NCBI's
    /// order (the plural names' lookups, then the singular names'), the texts of opening before
    /// the first lookup of the report that has them.
    pub fn row_messages(&self, spec: &TabularSpec, taxonomy: &TabularTaxonomy) -> Vec<u8> {
        let mut lookups: Vec<i32> = Vec::new();
        if spec.checks_taxdb() {
            lookups.extend(&taxonomy.taxids);
        }
        if spec.has_singular_names() && taxonomy.taxid != 0 {
            lookups.push(taxonomy.taxid);
        }
        let mut text = Vec::new();
        let mut push = |line: &str| {
            text.extend_from_slice(
                crate::seqdb::taxdb::err_post_line(&self.program, line).as_bytes(),
            );
            text.push(b'\n');
        };
        if !lookups.is_empty() && !self.opened.replace(true) {
            self.open.iter().for_each(|line| push(line));
        }
        for taxid in lookups {
            if let Some((_, messages)) = self.answers.get(&taxid) {
                messages.iter().for_each(|line| push(line));
            }
        }
        text
    }
}

/// `s_IsValidName`: a blast name or super kingdom of `-` or `unclassified` is not shown.
fn valid_name(name: &[u8]) -> bool {
    name != b"-" && name != b"unclassified"
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
    /// The names of the taxonomy database (`None`: every lookup fails).
    pub names: Option<&'a TabularNames>,
    /// Whether the subject is a nucleotide sequence (`sstrand` is `plus` / `minus`, else `N/A`).
    pub nucleotide: bool,
    /// `qcovs` and `qcovus` of the subject's HSPs; -1 for `N/A`.
    pub qcovs: i32,
    pub qcovus: i32,
    /// Whether NCBI's `GetSeqIdList` finds the subject's ids (`SubjectLookup`).
    pub lookup: SubjectLookup<'a>,
    /// The caller's index of the subject (its subject index or rank), which `TabularRows` keeps
    /// for the last row that found its subject (`TabularRows::last_found`).
    pub index: usize,
}

impl TabularSubject<'_> {
    /// Whether the subject's load throws (`SubjectLookup::Throws`): NCBI loads nothing for its
    /// rows, so no text of a load or of a taxonomy lookup comes with them.
    pub fn throws(&self) -> bool {
        matches!(self.lookup, SubjectLookup::Throws(_))
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:222-230
// ```c++
//         if (((*itr)->IsGeneral() &&
//             (*itr)->AsFastaString().find("gnl|BL_ORD_ID")
//             != string::npos) ||
// 		(*itr)->AsFastaString().find("lcl|Subject_") != string::npos) {
//             vector<string> title_tokens;
//             string defline = sequence::CDeflineGenerator().GenerateDefline(bh);
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:862-875
// ```c++
//         } catch (const CException&) {
//             list<CRef<CSeq_id> > subject_ids;
//             CRef<CSeq_id> id(new CSeq_id());
//             id->Assign(align.GetSeq_id(1));
//             subject_ids.push_back(id);
//             SetSubjectId(subject_ids);
//             bioseqs_found = false;
//         }
//
//     }
//
//     // If Bioseq has not been found for one or both of the sequences, all
//     // subsequent computations cannot proceed. Hence don't set any of the other
//     // fields.
//     if (!bioseqs_found)
//         return -1;
// ```
/// Whether NCBI's `GetSeqIdList` makes a subject's tabular ids, or throws: it puts the first word
/// of the subject Bioseq's generated title in the place of a `gnl|BL_ORD_ID` id or an id with
/// `lcl|Subject_`, and `GenerateDefline` throws for a title in no encoding that
/// `CUtf8::GuessEncoding` recognises (`report::defline::html_decode`). The catch of `SetFields`
/// then ends the row's fields (`FailedLookup`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SubjectLookup<'a> {
    #[default]
    Found,
    Failed(FailedLookup<'a>),
    /// The load of the subject throws (`ThrownLookup`), before `GetSeqIdList` runs.
    Throws(ThrownLookup<'a>),
}

/// A subject whose ids `GetSeqIdList` cannot make (`SubjectLookup::Failed`): the subject's
/// Seq-align id, which the catch pushes onto `m_SubjectIds` (`sallseqid`, `sallgi`, `sallacc`),
/// and where the lookup throws.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FailedLookup<'a> {
    /// The Seq-align id as `sallseqid` prints it (`eFullId`: `Subject_1`, `gnl|BL_ORD_ID|0`)
    /// and as `sallacc` prints it (`eAccession`: `Subject_1`, `BL_ORD_ID:0`); `sallgi` prints 0.
    pub full: &'a [u8],
    pub accession: &'a [u8],
    /// Whether `SetSubjectId` throws: `GetSeqIdList` of the subject Bioseq's own ids replaces
    /// one (`gnl|BL_ORD_ID`, `lcl|Subject_`) with the first word of a title that
    /// `GenerateDefline` cannot make. A database sequence whose Seq-align id is on a definition
    /// line without such an id has its ids (`gb|X.1|` before `lcl|Subject_9`).
    pub subject_id: bool,
    /// Where `x_SetSubjectIds` throws: `None` where it does not (a `-subject` record has no
    /// definition lines, and its ids are only `s_ReplaceLocalId`'s, tabular.cpp:557-564), else
    /// after the ids of the definition lines before the first whose `GetSeqIdList` throws
    /// (tabular.cpp:546-556: each line's list is pushed as it is made).
    pub all_ids: Option<&'a [DeflineIds]>,
}

impl FailedLookup<'_> {
    /// Which of NCBI's two lookups throw.
    pub fn failure(&self) -> LookupFailure {
        LookupFailure {
            subject_id: self.subject_id,
            all_ids: self.all_ids.is_some(),
            bioseq: false,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:833-838, 862-875 (quoted at
// `TabularSpec::row_throw`)
/// A subject whose load by its Seq-align id throws where the query's batch has its prefetch off
/// (`SubjectLookup::Throws`; `report::db_subject::HitLoad::Throws`): the Seq-align id that the
/// `catch` pushes onto `m_SubjectIds`, and the exception, which ends the run at a row that reads
/// the subject outside the subject block (`TabularWriter::ends_run`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThrownLookup<'a> {
    /// The Seq-align id as `sallseqid`, `sallgi` and `sallacc` print it (`eFullId`, `eGi`,
    /// `eAccession`; `ParsedSubject::alignment_tabular_ids`, `alignment_gi`).
    pub full: &'a [u8],
    pub gi: i64,
    pub accession: &'a [u8],
    pub exception: &'a crate::ncbi_exception::NcbiException,
}

/// Which of NCBI's lookups of a subject throw (`FailedLookup` without its ids): `SetSubjectId`
/// (`subject_id`) and `x_SetSubjectIds` (`all_ids`), the checks of a report's field list
/// (`failed_lookup_stop`) read only these.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct LookupFailure {
    pub subject_id: bool,
    pub all_ids: bool,
    /// Whether the load of the subject's Bioseq by its Seq-align id throws, before both
    /// (`SubjectLookup::Throws`: a lookup that fails on a file where the batch's prefetch is off).
    pub bioseq: bool,
}

impl LookupFailure {
    /// The failure of a subject without parsed Seq-ids whose generated title throws: a `-subject`
    /// record's (`SetSubjectId` only), or a database sequence's (`database`: its one definition
    /// line holds the `gnl|BL_ORD_ID` id, so `x_SetSubjectIds` throws at it too).
    pub fn of_plain(database: bool) -> Self {
        Self {
            subject_id: true,
            all_ids: database,
            bioseq: false,
        }
    }

    /// The failure of a subject whose load throws (`LookupFailure::bioseq`), whatever its other
    /// lookups would do (NCBI makes none of them).
    pub fn thrown() -> Self {
        Self {
            bioseq: true,
            ..Self::default()
        }
    }
}

/// Where `SetFields` stops for a row whose subject lookup failed (`TabularSpec::failed_stage`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FailedStage {
    /// The load of the subject (`GetBioseqHandle`, a field of the subject block asked): every
    /// member that the block sets stays as the rows before left it (`m_SubjectId`, the taxonomy,
    /// the titles), `m_SubjectLength` 0, and the Seq-align id is added to `m_SubjectIds`.
    Bioseq,
    /// `SetSubjectId` (`sseqid`, `sgi`, `sacc`, `saccver` asked): `m_SubjectId` stays empty,
    /// `m_SubjectLength` 0, and the Seq-align id is added to the `m_SubjectIds` of the rows before.
    SubjectId,
    /// `x_SetSubjectIds` (`sallseqid`, `sallgi`, `sallacc` asked, none of the above): after
    /// `m_SubjectIds.clear()`, so the Seq-align id is its only list; the length is set.
    SubjectIds,
}

/// The fields of a row that a row whose subject lookup failed leaves as the last row of the result
/// set that found its subject set them (`x_SetTaxInfoAll`, `x_SetTaxInfo` and `m_SubjectDefline`
/// come after the lookups; `x_ResetFields` does not reset them).
fn kept_by_failed_row(field: TabularField) -> bool {
    matches!(
        field,
        SubjectTaxId
            | SubjectSciName
            | SubjectCommonName
            | SubjectBlastName
            | SubjectSuperKingdom
            | SubjectTaxIds
            | SubjectSciNames
            | SubjectCommonNames
            | SubjectBlastNames
            | SubjectSuperKingdoms
            | SubjectTitle
            | SubjectAllTitles
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:780-786
// ```c++
//         CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec, kDelim);
//         if(!m_CustomDelim.empty()) {
//             tabinfo.SetCustomDelim(m_CustomDelim);
//         }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:101-114
// ```c++
// void CBlastTabularInfo::x_ResetFields()
// {
//     m_QueryLength = m_SubjectLength = 0U;
//     m_Score = m_AlignLength = m_NumGaps = m_NumGapOpens = m_NumIdent =
//     m_NumPositives = m_QueryStart = m_QueryEnd = m_SubjectStart =
//     m_SubjectEnd = m_QueryFrame = m_SubjectFrame = 0;
//     m_BitScore = NcbiEmptyString;
//     m_Evalue = NcbiEmptyString;
//     m_QuerySeq = NcbiEmptyString;
//     m_SubjectSeq = NcbiEmptyString;
//     m_BTOP = NcbiEmptyString;
//     m_SubjectStrand = NcbiEmptyString;
//     m_QueryCovSeqalign = -1;
// }
// ```
/// What NCBI's formatter of one result set (a query's report, or a psiblast round's: one
/// `CBlastTabularInfo` each) keeps from row to row: the subject ids of `sall*` (`m_SubjectIds`,
/// to which a row whose lookup fails adds its Seq-align id) and the fields that such a row leaves
/// as the last row that found its subject set them (`kept_by_failed_row`). Made empty for each
/// result set. A row that finds its subject stores only its subject's index (`last_found`); a
/// row whose lookup fails is given that subject again (`TabularWriter::row`'s `kept`) and writes
/// the kept fields from it, so that the rows that find their subjects copy nothing.
#[derive(Debug, Default)]
pub struct TabularRows {
    /// `TabularSubject::index` of the subject of the last row that found its subject.
    last_found: Option<usize>,
    /// `sallseqid`, `sallgi`, `sallacc` as the rows whose lookup failed since the last row that
    /// found its subject printed them; `None`: that row's own (its subject's ids).
    failed_ids: Option<[Vec<u8>; 3]>,
    /// Whether a row since that one emptied `m_SubjectId` (`FailedStage::SubjectId`), which a
    /// row whose load throws (`FailedStage::Bioseq`) then prints as it is.
    subject_id_emptied: bool,
}

impl TabularRows {
    /// The index (`TabularSubject::index`) of the subject of the result set's last row that found
    /// its subject, which a row whose lookup fails is given as `kept` (`TabularWriter::row`).
    pub fn last_found(&self) -> Option<usize> {
        self.last_found
    }

    /// A row that finds its subject, of index `index`: it sets `m_SubjectId`, the titles, the
    /// taxonomy and `m_SubjectIds` anew, which the rows after it keep where their lookups fail.
    /// The one record of it for the custom fields (`TabularWriter::row`) and the default fields
    /// (the programs' writers of `write_thrown_default` rows).
    pub fn found(&mut self, index: usize) {
        self.last_found = Some(index);
        self.failed_ids = None;
        self.subject_id_emptied = false;
    }

    /// The subject whose `m_SubjectId` a row whose load throws prints (`FailedStage::Bioseq`): the
    /// result set's last row that found its subject, `None` (empty: `sseqid` `Unknown`, `sgi` 0)
    /// before such a row or after a row that emptied it (`FailedStage::SubjectId`).
    pub fn kept_id(&self) -> Option<usize> {
        self.last_found.filter(|_| !self.subject_id_emptied)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:833-838, 862-875
    // (quoted at `TabularSpec::row_throw`: the subject block's `catch`, `return -1`)
    /// Writes the default fields' row (outfmt 6 and 7 without a field list,
    /// `outfmt6::write_hit_fields_bytes`) of the query `query_id` whose subject's load throws where
    /// its batch's prefetch is off (`report::db_subject::HitLoad::Throws`): the catch of the
    /// subject block ends the fields after the scores (`e_value`, `bit_score`), so the counts and
    /// coordinates keep the 0 of `x_ResetFields` (`pident` 0.000 of an alignment length 0), and
    /// `saccver` prints `m_SubjectId` as the kept subject set it (`kept_id`, whose `saccver` the
    /// caller gives, `saccver(index)`); before such a row `m_SubjectId` is NCBI's null pointer,
    /// refused before the report (`failed_lookup_stop`). The one writer of such rows of every
    /// program's default fields.
    pub fn write_thrown_default<W: Write + ?Sized, S: AsRef<[u8]>>(
        &self,
        writer: &mut W,
        query_id: &[u8],
        saccver: impl FnOnce(usize) -> S,
        e_value: f64,
        bit_score: f64,
        config: &crate::report::outfmt6::OutputConfig,
    ) -> io::Result<()> {
        let kept = self.kept_id().map(saccver);
        crate::report::outfmt6::write_hit_fields_bytes(
            writer,
            query_id,
            kept.as_ref().map_or(&[][..], AsRef::as_ref),
            0.0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            e_value,
            bit_score,
            config,
        )
    }
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
    /// Whether a `-` that begins or ends a row of `rows` is an end gap (`end_gap_span`), which NCBI
    /// shows as its end character, a space (`push_alignment_row`), as in the Smith-Waterman
    /// alignments of rpsblast's `-use_sw_tback`. NCBI decides by where the `-` comes from
    /// (alnvec.cpp:279-297): a gap segment outside the row's aligned segments is the end
    /// character, a fill of a segment that has a start is the gap character. `true` for every
    /// alignment of the engines; `false` only for `-subject`'s mapped row of BLASTX
    /// (`blastx::empty_rows::MappedAlignment`), whose `-` past the query's end fills its segment.
    pub end_gaps: bool,
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

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203 (the `CATCH_ALL` of
// every application, quoted at `cli::NativeError::exception`)
/// The end of the run at a row whose subject's load throws outside the subject block
/// (`TabularWriter::ends_run`): writes the report's rows before it, `buffered`, to `writer`,
/// ends the reads (`end_reads`: `QueryWarnings::end_run`, `ReadCursor::end_run`; NCBI loads
/// nothing after it) and gives NCBI's exception as the run's error
/// (`blast_xml::writer::run_ending_exception`: `BLAST Database error: <message>`, exit status 2).
pub fn run_end(
    writer: &mut dyn Write,
    buffered: &[u8],
    exception: &crate::ncbi_exception::NcbiException,
    end_reads: impl FnOnce(),
) -> anyhow::Error {
    if let Err(error) = writer.write_all(buffered).and_then(|()| writer.flush()) {
        return error.into();
    }
    end_reads();
    crate::report::blast_xml::writer::run_ending_exception(exception)
}

/// The writer of a tabular report.
pub struct TabularWriter<'a> {
    pub spec: &'a TabularSpec,
    /// The `# BLASTN 2.17.0+` line's program and version.
    pub program_version: &'a str,
    /// The matrix of a protein search's positives (`None` for BLASTN: positives are
    /// identities): the formatter's `m_ScoringMatrix`, which `CBlastFormat` makes with
    /// `GetAsciiProteinMatrix` (outfmt 0's matrix, `utils::matrix::protein_display_score`), not
    /// the XML's `CBlastFormattingMatrix`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:183
    /// ```c++
    ///     CAlignFormatUtil::GetAsciiProteinMatrix(m_MatrixName, m_ScoringMatrix);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:828-830
    /// ```c++
    ///                     tabinfo.SetFields(s, *m_Scope, &m_ScoringMatrix);
    /// ```
    pub matrix: Option<ScoringMatrix>,
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
        self.round_query_header(out, None, query_title, database, num_hits)
    }

    /// `query_header` of a report of psiblast's round `iteration` (`Some`), whose lines carry the
    /// round after the program's.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1294-1298
    /// ```c++
    ///     m_Ostream << "# ";
    ///     m_Ostream << program_version << "\n";
    ///
    ///     if (iteration != numeric_limits<unsigned int>::max())
    ///         m_Ostream << "# Iteration: " << iteration << "\n";
    /// ```
    pub fn round_query_header(
        &self,
        out: &mut dyn Write,
        iteration: Option<u32>,
        query_title: &[u8],
        database: &[u8],
        num_hits: Option<usize>,
    ) -> io::Result<()> {
        if self.spec.kind != TabularKind::Commented {
            return Ok(());
        }
        writeln!(out, "# {}", self.program_version)?;
        if let Some(iteration) = iteration {
            writeln!(out, "# Iteration: {iteration}")?;
        }
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

    /// Whether the row of `subject` sets the subject's values, as a row whose lookup does not fail
    /// (`TabularSpec::failed_stage`): only then does NCBI look its taxids up
    /// (`TabularNames::row_messages`).
    pub fn sets_subject(&self, subject: &TabularSubject<'_>) -> bool {
        self.spec.failed_stage(subject.lookup).is_none()
    }

    /// The exception that ends the run at `program`'s row of `subject` with the engine's
    /// identities `num_ident` (`TabularSpec::row_throw` is `EndsRun` for a subject whose load
    /// throws, `SubjectLookup::Throws`: its first row that reads it, the first row of HSPs with
    /// identities where `nident` alone reads it, as `report::db_subject::HitLoads::of` takes the
    /// run to end there): the row is not written, and the caller writes the rows before it and
    /// returns the exception (`run_end`); `None` for every other row (a row of such a subject that
    /// reads nothing of it is written as any other).
    pub fn ends_run<'s>(
        &self,
        subject: &TabularSubject<'s>,
        program: crate::blastinput::app::Program,
        num_ident: usize,
    ) -> Option<&'s crate::ncbi_exception::NcbiException> {
        let SubjectLookup::Throws(thrown) = subject.lookup else {
            return None;
        };
        (self.spec.row_throw(program, num_ident)
            == Some(crate::report::db_subject::AfterThrow::EndsRun))
        .then_some(thrown.exception)
    }

    /// Whether a row of `subject` fails its lookups (`TabularSpec::failed_stage`), so that it is
    /// written with the subject of the result set's last row that found its subject
    /// (`TabularRows::last_found`, the `kept` of `row`).
    pub fn fails(&self, subject: &TabularSubject<'_>) -> bool {
        !self.sets_subject(subject)
    }

    /// One row (`SetFields` then `Print`) of the result set whose formatter state is `rows`,
    /// appended to `out`. A row whose subject lookup fails (`fails`) is given `kept`, the subject
    /// of `rows.last_found()` (`None` before a row that found its subject), whose titles,
    /// taxonomy and `sall*` it shows as NCBI's members keep them; other rows ignore it.
    pub fn row(
        &self,
        out: &mut Vec<u8>,
        rows: &mut TabularRows,
        query: &TabularQuery<'_>,
        subject: &TabularSubject<'_>,
        kept: Option<&TabularSubject<'_>>,
        hsp: &TabularHsp<'_>,
    ) {
        if let Some(stage) = self.spec.failed_stage(subject.lookup) {
            self.failed_row(out, rows, query, subject, kept, hsp, stage);
            return;
        }
        rows.found(subject.index);
        self.found_row(out, query, subject, hsp);
    }

    /// The row of a subject whose lookup fails at `stage`: the fields before the catch (the
    /// scores, the query's ids and `qcovhsp`, `qcovs`, `qcovus`), the values of `x_ResetFields`,
    /// and the kept members of the result set (`TabularRows`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:771-797
    /// ```c++
    ///     if (x_IsFieldRequested(eQuerySeqId) || x_IsFieldRequested(eQueryGi) ||
    ///         x_IsFieldRequested(eQueryAccession) ||
    ///         x_IsFieldRequested(eQueryAccessionVersion) ||
    ///         x_IsFieldRequested(eQueryCovSeqalign)) {
    ///         try {
    ///             // FIXME: do this only if the query has changed
    ///             const CBioseq_Handle& query_bh =
    ///                 scope.GetBioseqHandle(align.GetSeq_id(0));
    ///             SetQueryId(query_bh);
    ///             if(m_QueryRange.NotEmpty())
    ///             	m_QueryLength = m_QueryRange.GetLength();
    ///             else
    ///             	m_QueryLength = query_bh.GetBioseqLength();
    ///             x_SetQueryCovSeqalign(align, m_QueryLength);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:169-194
    /// ```c++
    ///     case CBlastTabularInfo::eAccession:
    ///     {
    ///         CConstRef<CSeq_id> accid = FindBestChoice(id, CSeq_id::Score);
    ///         accid->GetLabel(&id_str, CSeq_id::eContent, 0);
    ///         break;
    ///     }
    ///     ...
    ///     case CBlastTabularInfo::eGi:
    ///         id_str = NStr::NumericToString(FindGi(id));
    ///         break;
    ///     default: break;
    ///     }
    ///
    ///     if (id_str == NcbiEmptyString)
    ///         id_str = "Unknown";
    /// ```
    /// `sseqid` of the empty `m_SubjectId` (`SetSubjectId` threw) is `Unknown` and `sgi` 0;
    /// `sacc` and `saccver` read a null `FindBestChoice` (NCBI stops: `failed_lookup_stop`, before
    /// the report), so no row of this writer prints them. Where only `x_SetSubjectIds` threw, the
    /// subject's own ids and length were set. Where the load threw (`FailedStage::Bioseq`), the
    /// subject ids are the kept subject's (the empty ones before a row that found its subject,
    /// or after a row that emptied them).
    #[allow(clippy::too_many_arguments)]
    #[cold]
    #[inline(never)]
    fn failed_row(
        &self,
        out: &mut Vec<u8>,
        rows: &mut TabularRows,
        query: &TabularQuery<'_>,
        subject: &TabularSubject<'_>,
        kept: Option<&TabularSubject<'_>>,
        hsp: &TabularHsp<'_>,
        stage: FailedStage,
    ) {
        let spec = self.spec;
        // The Seq-align id that the catch pushes, and where `x_SetSubjectIds` threw, the lists
        // of the definition lines before the one that throws.
        let (pushed, pushed_gi, lines): ([&[u8]; 2], i64, _) = match subject.lookup {
            SubjectLookup::Failed(failed) => ([failed.full, failed.accession], 0, failed.all_ids),
            SubjectLookup::Throws(thrown) => ([thrown.full, thrown.accession], thrown.gi, None),
            SubjectLookup::Found => ([&[][..], &[][..]], 0, None),
        };
        if stage == FailedStage::SubjectId {
            rows.subject_id_emptied = true;
        }
        // `m_SubjectId` as the kept subject set it (a load that throws leaves it), or empty.
        let kept_id = kept.filter(|_| stage == FailedStage::Bioseq && !rows.subject_id_emptied);
        let all_ids = spec.has_any(&[SubjectAllSeqIds, SubjectAllGis, SubjectAllAccessions]);
        if all_ids {
            // `m_SubjectIds` as the last row left it: the kept subject's lists, or those of the
            // failed rows since.
            let lists = rows.failed_ids.get_or_insert_with(|| {
                let mut lists: [Vec<u8>; 3] = Default::default();
                if let Some(kept) = kept {
                    for (list, field) in lists.iter_mut().zip([
                        SubjectAllSeqIds,
                        SubjectAllGis,
                        SubjectAllAccessions,
                    ]) {
                        push_all_ids(list, field, kept);
                    }
                }
                lists
            });
            // `x_SetSubjectIds` clears the lists and pushes the lists of the definition lines
            // before the one that throws; `SetSubjectId(list)` of the catch pushes the Seq-align id
            // (tabular.cpp:1374-1378).
            if stage == FailedStage::SubjectIds {
                lists.iter_mut().for_each(Vec::clear);
                for line in lines.unwrap_or_default() {
                    for (list, field) in lists.iter_mut().zip([
                        SubjectAllSeqIds,
                        SubjectAllGis,
                        SubjectAllAccessions,
                    ]) {
                        if !list.is_empty() {
                            list.push(b';');
                        }
                        push_defline_id(list, field, line);
                    }
                }
            }
            for (n, list) in lists.iter_mut().enumerate() {
                if !list.is_empty() {
                    list.push(b';');
                }
                match n {
                    0 => list.extend_from_slice(pushed[0]),
                    1 => push_num(list, pushed_gi),
                    _ => list.extend_from_slice(pushed[1]),
                }
            }
        }
        // The query's block of `SetFields` sets `m_QueryLength` and `qcovhsp`; `qlen`'s own
        // block comes after the catch's `return -1`.
        let query_block = spec.has_any(&[
            QuerySeqId,
            QueryGi,
            QueryAccession,
            QueryAccessionVersion,
            QueryCovHsp,
        ]);
        let names = kept.map(|kept| RowNames::of(spec, kept));
        let subject_set = stage == FailedStage::SubjectIds;
        for (i, field) in spec.fields.iter().enumerate() {
            if i > 0 {
                out.extend_from_slice(spec.delimiter.as_bytes());
            }
            if kept_by_failed_row(*field) {
                match (kept, &names) {
                    (Some(kept), Some(names)) => self.subject_field(out, *field, kept, names),
                    // Before a row that found its subject: the members as constructed (empty).
                    // `staxid` reads an uninitialised `m_SubjectTaxId` there
                    // (`failed_lookup_stop`, before the report).
                    _ => {
                        let quoted = spec.is_csv()
                            && matches!(
                                field,
                                SubjectSciNames | SubjectCommonNames | SubjectBlastNames
                            );
                        if quoted {
                            write_quoted(out, NA);
                        } else {
                            out.extend_from_slice(NA);
                        }
                    }
                }
                continue;
            }
            if let Some(kept) = kept_id {
                match field {
                    SubjectSeqId => out.extend_from_slice(kept.id),
                    SubjectAccession => out.extend_from_slice(kept.accession),
                    SubjectAccessionVersion => out.extend_from_slice(kept.accession_version),
                    SubjectGi => push_num(out, kept.gi),
                    _ => {}
                }
                if matches!(
                    field,
                    SubjectSeqId | SubjectAccession | SubjectAccessionVersion | SubjectGi
                ) {
                    continue;
                }
            }
            let lists = rows.failed_ids.as_ref();
            let list = |n: usize| lists.map_or(&[][..], |lists| lists[n].as_slice());
            match field {
                QuerySeqId => out.extend_from_slice(query.id),
                QueryAccession => out.extend_from_slice(query.accession),
                QueryAccessionVersion => out.extend_from_slice(query.accession_version),
                QueryGi => push_num(out, query.gi as usize),
                QueryLength => push_num(out, if query_block { query.length } else { 0 }),
                SubjectSeqId if subject_set => out.extend_from_slice(subject.id),
                SubjectAccession if subject_set => out.extend_from_slice(subject.accession),
                SubjectAccessionVersion if subject_set => {
                    out.extend_from_slice(subject.accession_version)
                }
                SubjectGi if subject_set => push_num(out, subject.gi),
                SubjectSeqId => out.extend_from_slice(b"Unknown"),
                SubjectGi => push_num(out, 0),
                SubjectAccession | SubjectAccessionVersion => {}
                SubjectAllSeqIds => out.extend_from_slice(list(0)),
                SubjectAllGis => out.extend_from_slice(list(1)),
                SubjectAllAccessions => out.extend_from_slice(list(2)),
                SubjectLength => push_num(out, if subject_set { subject.length } else { 0 }),
                Evalue => out.extend_from_slice(format_evalue_ncbi_tabular(hsp.evalue).as_bytes()),
                BitScore => out.extend_from_slice(format_bitscore_ncbi(hsp.bit_score).as_bytes()),
                Score => push_num(out, hsp.score),
                PercentIdentical => out.extend_from_slice(b"0.000"),
                PercentPositives => out.extend_from_slice(b"0.00"),
                Frames => out.extend_from_slice(b"0/0"),
                QuerySeq | SubjectSeq | Btop => {}
                SubjectStrand => out.extend_from_slice(NA),
                QueryCovSubject => push_coverage(out, subject.qcovs),
                QueryCovUniqSubject => push_coverage(out, subject.qcovus),
                QueryCovHsp => push_coverage(out, hsp.qcovhsp),
                // The other counts and coordinates keep `x_ResetFields`' 0.
                _ => push_num(out, 0),
            }
        }
        out.push(b'\n');
    }

    /// The row of a subject whose lookup does not fail. Inline in its one caller (`row`), so that
    /// such a row costs no call more than before the rows whose lookups fail (no NCBI
    /// counterpart: code placement only).
    #[inline(always)]
    fn found_row(
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
                    } else if self.matrix.is_some_and(|matrix| {
                        crate::utils::matrix::protein_display_score(matrix, q, s) > 0
                    }) {
                        positive += 1;
                    }
                }
                (ident, positive)
            }
            None => (hsp.num_ident, hsp.num_ident),
        };
        let spec = self.spec;
        // Only NCBI's alignment block (`TabularSpec::opens_counts`) sets the frames (`SetCounts`);
        // without one of its fields `frames` and `sframe` keep the 0 of `x_ResetFields`.
        let counts_set = spec.opens_counts(&hsp.num_ident);
        let frames = |hsp: &TabularHsp<'_>| {
            if counts_set {
                (hsp.q_frame, hsp.s_frame)
            } else {
                (0, 0)
            }
        };
        let names = RowNames::of(spec, subject);
        for (i, field) in spec.fields.iter().enumerate() {
            if i > 0 {
                out.extend_from_slice(spec.delimiter.as_bytes());
            }
            match field {
                QuerySeqId => out.extend_from_slice(query.id),
                QueryAccession => out.extend_from_slice(query.accession),
                QueryAccessionVersion => out.extend_from_slice(query.accession_version),
                QueryGi => push_num(out, query.gi as usize),
                QueryLength => push_num(out, query.length),
                SubjectSeqId => out.extend_from_slice(subject.id),
                SubjectAccession => out.extend_from_slice(subject.accession),
                SubjectAccessionVersion => out.extend_from_slice(subject.accession_version),
                SubjectGi => push_num(out, subject.gi),
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:224-231
                // ```c++
                // void CBlastTabularInfo::x_PrintSubjectAllSeqIds(void)
                // {
                //     ITERATE(vector<list<CRef<CSeq_id> > >, iter, m_SubjectIds) {
                //         if (iter != m_SubjectIds.begin())
                //             m_Ostream << ";";
                //         m_Ostream << s_GetSeqIdListString(*iter, eFullId);
                //     }
                // }
                // ```
                // `x_PrintSubjectAllGis` and `x_PrintSubjectAllAccessions` (tabular.cpp:238-264)
                // are the same loop with `eGi` and `eAccession`.
                SubjectAllSeqIds | SubjectAllGis | SubjectAllAccessions => {
                    push_all_ids(out, *field, subject)
                }
                SubjectLength => push_num(out, subject.length),
                QueryStart => push_num(out, hsp.q_start),
                QueryEnd => push_num(out, hsp.q_end),
                SubjectStart => push_num(out, hsp.s_start),
                SubjectEnd => push_num(out, hsp.s_end),
                QuerySeq => {
                    push_alignment_row(out, hsp.rows.map_or(&[][..], |rows| rows.0), hsp.end_gaps)
                }
                SubjectSeq => {
                    push_alignment_row(out, hsp.rows.map_or(&[][..], |rows| rows.1), hsp.end_gaps)
                }
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
                    let (q_frame, s_frame) = frames(hsp);
                    push_num(out, q_frame);
                    out.push(b'/');
                    push_num(out, s_frame);
                }
                QueryFrame => push_num(out, hsp.q_frame),
                SubjectFrame => push_num(out, frames(hsp).1),
                Btop => {
                    if let Some((qrow, srow)) = hsp.rows {
                        write_btop(out, qrow, srow, hsp.end_gaps);
                    }
                }
                SubjectTaxId | SubjectTaxIds | SubjectSciName | SubjectCommonName
                | SubjectBlastName | SubjectSuperKingdom | SubjectSciNames | SubjectCommonNames
                | SubjectBlastNames | SubjectSuperKingdoms | SubjectTitle | SubjectAllTitles => {
                    self.subject_field(out, *field, subject, &names)
                }
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

    /// A field of `subject` that a row whose lookup fails keeps from the last row that found its
    /// subject (`kept_by_failed_row`: the taxonomy and the titles), with the row's names `names`.
    #[inline(always)]
    fn subject_field(
        &self,
        out: &mut Vec<u8>,
        field: TabularField,
        subject: &TabularSubject<'_>,
        names: &RowNames<'_>,
    ) {
        let csv = self.spec.is_csv();
        match field {
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
            // The names of `RowNames`: an empty singular name is `N/A`, an empty list of
            // plural names `N/A`; the scientific, common and blast names are quoted in a CSV
            // (`NStr::Quote`, so `N/A` too and an empty name not), the super kingdoms never.
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
            SubjectSciName | SubjectCommonName | SubjectBlastName | SubjectSuperKingdom => {
                let name = match field {
                    SubjectSciName => names.sci,
                    SubjectCommonName => names.common,
                    SubjectBlastName => names.blast,
                    _ => names.kingdom,
                };
                if name.is_empty() {
                    out.extend_from_slice(NA);
                } else if csv && field != SubjectSuperKingdom {
                    write_quoted(out, name);
                } else {
                    out.extend_from_slice(name);
                }
            }
            SubjectSciNames | SubjectCommonNames | SubjectBlastNames | SubjectSuperKingdoms => {
                let list = match field {
                    SubjectSciNames => &names.scis,
                    SubjectCommonNames => &names.commons,
                    SubjectBlastNames => &names.blasts,
                    _ => &names.kingdoms,
                };
                let joined = if list.is_empty() {
                    NA.to_vec()
                } else {
                    list.join(&b';')
                };
                if csv && field != SubjectSuperKingdoms {
                    write_quoted(out, &joined);
                } else {
                    out.extend_from_slice(&joined);
                }
            }
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:387-440
            // ```c++
            // void CBlastTabularInfo::x_PrintSubjectAllTitles ()
            // {
            // 	if(m_SubjectDefline.NotEmpty() && m_SubjectDefline->CanGet() &&
            // 	   m_SubjectDefline->IsSet() && !m_SubjectDefline->Get().empty())
            // 	{
            // 		const list<CRef<CBlast_def_line> > & defline = m_SubjectDefline->Get();
            // 		list<CRef<CBlast_def_line> >::const_iterator iter = defline.begin();
            // 		for(; iter != defline.end(); ++iter)
            // 		{
            // 			if (iter != defline.begin())
            // 						 m_Ostream << "<>";
            //
            // 			if((*iter)->IsSetTitle())
            // 			{
            // 				if((*iter)->GetTitle().empty())
            // 					m_Ostream << NA;
            // 				else
            // 					m_Ostream << (*iter)->GetTitle();
            // 			}
            // 			else
            // 				m_Ostream << NA;
            // 		}
            // 	}
            // 	else
            // 		m_Ostream << NA;
            //
            // }
            //
            // void CBlastTabularInfo::x_PrintSubjectTitle()
            // ...
            // 			if(defline.front()->IsSetTitle())
            // 			{
            // 				if(defline.front()->GetTitle().empty())
            // 					m_Ostream << NA;
            // 				else
            // 					m_Ostream << defline.front()->GetTitle();
            // 			}
            // 			else
            // 				m_Ostream << NA;
            // ```
            // `x_PrintSubjectTitle` is the `stitle` arm below, `x_PrintSubjectAllTitles` the
            // `salltitles` arm after it. An absent title is empty in `titles`, so both read as
            // `N/A`.
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
            _ => unreachable!("not a field that a failed row keeps"),
        }
    }
}

/// `sallseqid`, `sallgi` or `sallacc` (`field`) of `subject` as its row prints them: the ids of
/// its definition lines joined by `;`, or a subject's own ids where it has none.
#[inline(always)]
fn push_all_ids(out: &mut Vec<u8>, field: TabularField, subject: &TabularSubject<'_>) {
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
        push_defline_id(out, field, ids);
    }
}

/// The id of one definition line that `sallseqid`, `sallgi` or `sallacc` (`field`) prints.
#[inline(always)]
fn push_defline_id(out: &mut Vec<u8>, field: TabularField, ids: &DeflineIds) {
    match field {
        SubjectAllSeqIds => out.extend_from_slice(&ids.id),
        SubjectAllGis => push_num(out, ids.gi),
        _ => out.extend_from_slice(&ids.accession),
    }
}

/// The taxonomy names of one row (`x_SetTaxInfoAll` and `x_SetTaxInfo`, after `x_ResetFields`):
/// empty unless a field shows them. A lookup that fails (no taxonomy database, a taxid it does not
/// know) leaves the names empty, except that with `staxids` among the fields NCBI fills the
/// scientific and common names with `N/A` so that their list matches the taxids.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:612-633
/// ```c++
///         try {
///         	SSeqDBTaxInfo taxinfo;
///             CSeqDB::GetTaxInfo(m_SubjectTaxId, taxinfo);
///             m_SubjectSciName = taxinfo.scientific_name;
///             m_SubjectCommonName = taxinfo.common_name;
///             if(s_IsValidName(taxinfo.blast_name)) {
///             	m_SubjectBlastName = taxinfo.blast_name;
///             }
///
///             if(s_IsValidName(taxinfo.s_kingdom)) {
///             	m_SubjectSuperKingdom = taxinfo.s_kingdom;
///             }
///
///         } catch (const CException&) {
///         	//only put fillers in if we are going to show tax id
///             // the fillers are put in so that the name list would
///             // match the taxid list
///             if(x_IsFieldRequested(eSubjectTaxIds) ) {
///             	m_SubjectSciName = NA;
///                 m_SubjectCommonName = NA;
///             }
///         }
/// ```
#[derive(Default)]
struct RowNames<'a> {
    sci: &'a [u8],
    common: &'a [u8],
    blast: &'a [u8],
    kingdom: &'a [u8],
    /// The plural names: by taxid of `staxids` (scientific, common), and the sorted sets of the
    /// blast names and super kingdoms (`set<string>`).
    scis: Vec<&'a [u8]>,
    commons: Vec<&'a [u8]>,
    blasts: Vec<&'a [u8]>,
    kingdoms: Vec<&'a [u8]>,
}

impl<'a> RowNames<'a> {
    #[inline(always)]
    fn of(spec: &TabularSpec, subject: &TabularSubject<'a>) -> Self {
        let mut names = RowNames::default();
        let lookup = |taxid: i32| subject.names.and_then(|names| names.get(taxid));
        let fillers = spec.has(SubjectTaxIds);
        if spec.checks_taxdb() {
            for &taxid in &subject.taxonomy.taxids {
                match lookup(taxid) {
                    Some(found) => {
                        names.scis.push(&found.scientific_name);
                        names.commons.push(&found.common_name);
                        if valid_name(&found.blast_name) {
                            names.blasts.push(&found.blast_name);
                        }
                        // A kingdom past NCBI's mapping of taxdb.btd is not written
                        // (`TabularNames::look_up` refuses a field that writes it).
                        if let Ok(kingdom) = found.super_kingdom.as_deref() {
                            if valid_name(kingdom) {
                                names.kingdoms.push(kingdom);
                            }
                        }
                    }
                    None if fillers => {
                        names.scis.push(NA);
                        names.commons.push(NA);
                    }
                    None => {}
                }
            }
            for set in [&mut names.blasts, &mut names.kingdoms] {
                set.sort_unstable();
                set.dedup();
            }
        }
        let taxid = subject.taxonomy.taxid;
        if spec.has_singular_names() && taxid != 0 {
            match lookup(taxid) {
                Some(found) => {
                    names.sci = &found.scientific_name;
                    names.common = &found.common_name;
                    if valid_name(&found.blast_name) {
                        names.blast = &found.blast_name;
                    }
                    if let Ok(kingdom) = found.super_kingdom.as_deref() {
                        if valid_name(kingdom) {
                            names.kingdom = kingdom;
                        }
                    }
                }
                None if fillers => {
                    names.sci = NA;
                    names.common = NA;
                }
                None => {}
            }
        }
        names
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
/// The rows' end gaps are spaces (`push_alignment_row`), so an end-gap column is the pair of a
/// letter and a space. Rows without end gaps (every alignment but some Smith-Waterman redos', and
/// rows whose `-` fill a segment, `end_gaps` false: `TabularHsp::end_gaps`) take the columns as
/// they are; `end_gap_span` stops at each row's first letter from each end.
pub fn write_btop(out: &mut Vec<u8>, qrow: &[u8], srow: &[u8], end_gaps: bool) {
    let columns = qrow.iter().copied().zip(srow.iter().copied());
    if !end_gaps {
        write_btop_columns(out, columns);
        return;
    }
    let (q_span, s_span) = (end_gap_span(qrow, |_| false), end_gap_span(srow, |_| false));
    if q_span == (0, qrow.len()) && s_span == (0, srow.len()) {
        write_btop_columns(out, columns);
    } else {
        write_btop_columns(
            out,
            columns
                .enumerate()
                .map(|(column, (q, s))| (end_char(q, column, q_span), end_char(s, column, s_span))),
        );
    }
}

/// `write_btop` of the column pairs `columns` (query letter, subject letter).
#[inline(always)]
fn write_btop_columns(out: &mut Vec<u8>, columns: impl Iterator<Item = (u8, u8)>) {
    let mut matches = 0usize;
    for (q, s) in columns {
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

/// The letter of column `column` of a row whose columns between its end gaps are `span`
/// (`end_gap_span`): NCBI's end character, a space, outside them.
fn end_char(letter: u8, column: usize, span: (usize, usize)) -> u8 {
    if column < span.0 || column >= span.1 {
        b' '
    } else {
        letter
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:982-986
// ```c++
//         alnVec->SetGapChar('-');
//         alnVec->SetGenCode(m_QueryGeneticCode, 0);
//         alnVec->SetGenCode(m_DbGeneticCode, 1);
//         alnVec->GetWholeAlnSeqString(0, m_QuerySeq);
//         alnVec->GetWholeAlnSeqString(1, m_SubjectSeq);
// ```
/// `qseq` or `sseq`: the row with `-` for its inner gaps and a space (NCBI's end character, which
/// the tabular writer does not set) for its end gaps (`end_gap_span`), where it has some
/// (`end_gaps`, `TabularHsp::end_gaps`).
fn push_alignment_row(out: &mut Vec<u8>, row: &[u8], end_gaps: bool) {
    let span = if end_gaps {
        end_gap_span(row, |_| false)
    } else {
        (0, row.len())
    };
    if span == (0, row.len()) {
        out.extend_from_slice(row);
        return;
    }
    out.extend(
        row.iter()
            .enumerate()
            .map(|(column, &letter)| end_char(letter, column, span)),
    );
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

    // Line DBID P1: with NCBI's prefetch off, the rows read the subject's Bioseq (and look its
    // taxids up) only for these fields (review 1 of P1, field by field against NCBI 2.17.0); P2
    // (review 1): in the translated programs also for the subject's frame with a coordinate; L1
    // (review 2 of P2): `nident` without another field of the block only on rows with identities.
    #[test]
    fn rows_read_the_subject_for_its_ids_taxonomy_titles_strand_and_letters() {
        use crate::blastinput::app::Program;
        let spec = |fields: &str| TabularSpec {
            kind: TabularKind::Tabular,
            fields: parse_fields(fields),
            delimiter: "\t".into(),
        };
        let fetch = |fields: &str, no_fetch: bool| {
            spec(fields).subject_fetch(if no_fetch {
                Program::Blastn
            } else {
                Program::Blastp
            })
        };
        for (fields, translated) in [
            ("qseqid sframe", SubjectFetch::Never),
            ("qseqid frames qcovs", SubjectFetch::Never),
            ("qseqid sstart sframe", SubjectFetch::Every),
            ("qseqid qframe frames", SubjectFetch::Every),
            ("qseqid length sframe", SubjectFetch::Every),
            ("qseqid qstart qframe slen", SubjectFetch::Never),
            ("qseqid nident sframe", SubjectFetch::WithIdentities),
            ("qseqid nident frames sstart", SubjectFetch::Every),
        ] {
            for program in [Program::Blastx, Program::Tblastn, Program::Tblastx] {
                assert_eq!(spec(fields).subject_fetch(program), translated, "{fields}");
            }
            let protein = match translated {
                SubjectFetch::Every if !fields.contains("nident") => SubjectFetch::Never,
                fetch => fetch,
            };
            assert_eq!(
                spec(fields).subject_fetch(Program::Blastp),
                protein,
                "{fields}"
            );
            assert_eq!(
                spec(fields).subject_fetch(Program::Blastn),
                SubjectFetch::Never,
                "{fields}"
            );
        }
        for fields in [
            "qseqid evalue",
            "qseqid qacc qlen slen qstart qend sstart send evalue bitscore score length gaps gapopen",
            "qseqid qframe sframe frames qcovs qcovhsp qcovus",
        ] {
            assert_eq!(fetch(fields, false), SubjectFetch::Never, "{fields}");
        }
        for fields in [
            "std",
            "qseqid sseqid",
            "qseqid sallacc",
            "qseqid staxid",
            "qseqid sskingdoms",
            "qseqid salltitles",
            "qseqid sstrand",
            "qseqid sseq",
            "qseqid btop",
            "qseqid ppos",
            "qseqid nident sseqid",
        ] {
            assert_eq!(fetch(fields, false), SubjectFetch::Every, "{fields}");
            assert_eq!(fetch(fields, true), SubjectFetch::Every, "{fields}");
        }
        for fields in ["qseqid mismatch", "qseqid pident", "qseqid nident length"] {
            assert_eq!(fetch(fields, false), SubjectFetch::Every, "{fields}");
            assert_eq!(fetch(fields, true), SubjectFetch::Never, "{fields}");
        }
        for fields in ["qseqid nident", "nident evalue qcovs slen"] {
            assert_eq!(
                fetch(fields, false),
                SubjectFetch::WithIdentities,
                "{fields}"
            );
            assert_eq!(fetch(fields, true), SubjectFetch::Never, "{fields}");
        }
        assert!(!SubjectFetch::WithIdentities.reads(0) && SubjectFetch::WithIdentities.reads(1));
        assert!(SubjectFetch::Every.reads(0) && !SubjectFetch::Never.reads(1));
    }

    // Line DBID tabfields: a table reads the sequence that NCBI finds by the subject's Seq-align
    // id where its rows read the Bioseq or letters, or its length (`slen`); NCBI 2.17.0 writes
    // every other list the same for a subject it finds as another sequence or not at all (probe
    // `blastdb-ids-20261021/R/tabfields`).
    #[test]
    fn tables_read_the_found_sequence_for_the_bioseq_its_letters_and_slen() {
        use crate::blastinput::app::Program;
        let spec = |fields: &str| TabularSpec {
            kind: TabularKind::Tabular,
            fields: parse_fields(fields),
            delimiter: "\t".into(),
        };
        let all = [
            Program::Blastp,
            Program::Blastn,
            Program::Blastx,
            Program::Tblastn,
            Program::Tblastx,
        ];
        assert!(all.iter().all(|&program| reads_subject(None, program)));
        for fields in [
            "qseqid evalue",
            "qseqid qlen qstart qend sstart send evalue bitscore score length gaps gapopen qcovhsp",
            "qseqid qcovs qcovus",
        ] {
            for program in all {
                assert!(!spec(fields).reads_subject(program), "{fields} {program:?}");
                assert!(!reads_subject(Some(&spec(fields)), program), "{fields}");
            }
        }
        for fields in [
            "std",
            "qseqid sseqid",
            "qseqid slen",
            "evalue stitle",
            "qseqid sseq",
        ] {
            for program in all {
                assert!(spec(fields).reads_subject(program), "{fields} {program:?}");
            }
        }
        // BLASTN's formatter counts the identities without the letters; the translated
        // programs read the subject for its frame.
        for fields in ["qseqid pident", "qseqid nident", "qseqid mismatch length"] {
            assert!(!spec(fields).reads_subject(Program::Blastn), "{fields}");
            assert!(spec(fields).reads_subject(Program::Blastp), "{fields}");
        }
        assert!(spec("qseqid sstart sframe").reads_subject(Program::Tblastn));
        assert!(!spec("qseqid sstart sframe").reads_subject(Program::Blastp));
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

    // NCBI 2.17.0 (line RG F/d36 probes ps.f10at, ps.f20at, ps.f6at, ps.pssm.f10at, bp.f20at,
    // bp.f6at, rps.f10at, rps.f6at, rpst.f10at): `-outfmt "10 delim=@ qseqid sseqid evalue"`
    // writes `,` in psiblast, rpsblast and rpstblastn, `@` in blastp; `6 delim=@` is the default
    // outfmt 6 of the first three.
    #[test]
    fn only_the_apps_that_set_the_delimiter_take_it() {
        use crate::report::archive::options::RequestKind;
        let choice = |number: i32, spec: &str| FormatChoice {
            number,
            spec: spec.to_string(),
            delimiter: "@".to_string(),
        };
        let csv = choice(10, "qseqid sseqid evalue");
        let std6 = choice(6, "");
        let header = choice(20, "qseqid");
        for app in [
            RequestKind::Blastn,
            RequestKind::Blastp,
            RequestKind::Blastx,
            RequestKind::Tblastn,
            RequestKind::PsiTblastn,
            RequestKind::Tblastx,
        ] {
            let spec = TabularSpec::of_app(ReportFormat::CommaSeparatedValues, &csv, app);
            assert_eq!(
                spec,
                TabularSpec::of(ReportFormat::CommaSeparatedValues, &csv)
            );
            assert_eq!(spec.unwrap().delimiter, "@", "{app:?}");
            let spec = TabularSpec::of_app(ReportFormat::Tabular, &std6, app).unwrap();
            assert_eq!(spec.delimiter, "@", "{app:?}");
        }
        for app in [
            RequestKind::Psiblast,
            RequestKind::Rpsblast,
            RequestKind::Rpstblastn,
        ] {
            let spec = TabularSpec::of_app(ReportFormat::CommaSeparatedValues, &csv, app).unwrap();
            assert_eq!(
                (spec.kind, spec.delimiter.as_str(), spec.fields.len()),
                (TabularKind::Csv, ",", 3),
                "{app:?}"
            );
            let spec =
                TabularSpec::of_app(ReportFormat::CommaSeparatedValuesWithHeader, &header, app);
            assert_eq!(spec.unwrap().delimiter, ",", "{app:?}");
            assert!(TabularSpec::of_app(ReportFormat::Tabular, &std6, app).is_none());
            let commented = choice(7, "qseqid");
            let spec = TabularSpec::of_app(ReportFormat::TabularWithComments, &commented, app);
            assert_eq!(spec.unwrap().delimiter, "\t", "{app:?}");
        }
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
        let row = |format: ReportFormat,
                   fields: &str,
                   taxonomy: &TabularTaxonomy,
                   names: Option<&TabularNames>| {
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
                names,
                nucleotide: true,
                qcovs: -1,
                qcovus: -1,
                lookup: SubjectLookup::Found,
                index: 0,
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
                end_gaps: true,
                qcovhsp: 100,
            };
            let mut out = Vec::new();
            writer.row(
                &mut out,
                &mut TabularRows::default(),
                &TabularQuery::local(b"q", 10),
                &subject,
                None,
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
                &two,
                None
            ),
            "\"N/A\",\"N/A\",N/A,N/A,\"N/A;N/A\",\"N/A;N/A\",\"N/A\",N/A,9606;10090\n"
        );
        assert_eq!(
            row(
                ReportFormat::Tabular,
                &format!("{fields} staxids"),
                &two,
                None
            ),
            "N/A\tN/A\tN/A\tN/A\tN/A;N/A\tN/A;N/A\tN/A\tN/A\t9606;10090\n"
        );
        assert_eq!(
            row(ReportFormat::CommaSeparatedValues, fields, &two, None),
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
                &none,
                None
            ),
            "N/A,\"N/A\",N/A,0\n"
        );
        // NCBI tabular.cpp:612-624 and 669-680 with the names found (RESULTS.md §3 of line OF
        // taxdb): the plural blast names and super kingdoms are sorted sets, `-` and
        // `unclassified` are dropped, an empty common name stays empty in the plural list and is
        // `N/A` alone; a taxid not found gets the fillers.
        let entry = |taxid, sci: &str, common: &str, blast: &str, kingdom: &str| {
            let names = crate::seqdb::taxdb::TaxNames {
                taxid,
                scientific_name: sci.into(),
                common_name: common.into(),
                blast_name: blast.into(),
                super_kingdom: Ok(kingdom.into()),
            };
            (taxid, (Some(names), Vec::new()))
        };
        let found = TabularNames {
            program: "blastp".into(),
            open: Vec::new(),
            opened: std::cell::Cell::new(false),
            answers: [
                entry(9606, "Homo sapiens", "human", "primates", "Eukaryota"),
                entry(7777, "Quote \"x\", y; z", "", "-", "Eukaryota"),
                entry(32630, "synthetic construct", "", "unclassified", "-"),
            ]
            .into_iter()
            .collect(),
            report: None,
        };
        let three = TabularTaxonomy {
            taxid: 9606,
            taxids: vec![7777, 9606, 10090, 32630],
        };
        assert_eq!(
            row(
                ReportFormat::CommaSeparatedValues,
                &format!("{fields} staxids"),
                &three,
                Some(&found)
            ),
            "\"Homo sapiens\",\"human\",\"primates\",Eukaryota,\"Quote \\\"x\\\", y; z;Homo sapiens;N/A;synthetic construct\",\";human;N/A;\",\"primates\",Eukaryota,7777;9606;10090;32630\n"
        );
        let quote = TabularTaxonomy {
            taxid: 7777,
            taxids: vec![7777],
        };
        assert_eq!(
            row(ReportFormat::Tabular, fields, &quote, Some(&found)),
            "Quote \"x\", y; z\tN/A\tN/A\tEukaryota\tQuote \"x\", y; z\t\tN/A\tEukaryota\n"
        );
    }

    // NCBI writes the texts of opening at the first lookup of the report, and a lookup's own
    // texts at each lookup: the plural names' lookups before the singular one's.
    #[test]
    fn row_messages_follow_the_lookups() {
        let spec = TabularSpec::of(
            ReportFormat::Tabular,
            &FormatChoice {
                number: 6,
                spec: "ssciname sscinames".to_string(),
                delimiter: String::new(),
            },
        )
        .unwrap();
        let names = TabularNames {
            program: "blastn".into(),
            open: vec!["open".into()],
            opened: std::cell::Cell::new(false),
            answers: [
                (99990, (None, vec!["offset".into()])),
                (9606, (None, Vec::new())),
            ]
            .into_iter()
            .collect(),
            report: None,
        };
        let none = TabularTaxonomy::default();
        assert!(names.row_messages(&spec, &none).is_empty());
        let last = TabularTaxonomy {
            taxid: 99990,
            taxids: vec![9606, 99990],
        };
        assert_eq!(
            String::from_utf8(names.row_messages(&spec, &last)).unwrap(),
            "Error: [blastn] open\nError: [blastn] offset\nError: [blastn] offset\n"
        );
        assert_eq!(
            String::from_utf8(names.row_messages(&spec, &last)).unwrap(),
            "Error: [blastn] offset\nError: [blastn] offset\n"
        );
    }

    // NCBI tabular.cpp:900-908 (the frames only with a field of `SetCounts`' block, `num_ident >
    // 0` for `nident`), tabular.cpp:1003-1021 (the positives of the rows by the formatter's
    // matrix, `GetAsciiProteinMatrix`: U and C count as C-C, a J pair not, `*` with `*` is 1)
    // and tabular.cpp:1293-1298 (psiblast's round).
    #[test]
    fn protein_rows_frames_and_rounds_are_ncbis() {
        let spec = |format: ReportFormat, fields: &str| {
            TabularSpec::of(
                format,
                &FormatChoice {
                    number: 6,
                    spec: fields.to_string(),
                    delimiter: String::new(),
                },
            )
            .unwrap()
        };
        let taxonomy = TabularTaxonomy::default();
        let subject = TabularSubject {
            id: b"s",
            accession: b"s",
            accession_version: b"s",
            gi: 0,
            all_ids: &[],
            length: 10,
            titles: None,
            taxonomy: &taxonomy,
            names: None,
            nucleotide: false,
            qcovs: -1,
            qcovus: -1,
            lookup: SubjectLookup::Found,
            index: 0,
        };
        let row = |spec: &TabularSpec, num_ident: usize, rows: Option<(&[u8], &[u8])>| {
            let writer = TabularWriter {
                spec,
                program_version: "BLASTP 2.17.0+",
                matrix: Some(ScoringMatrix::Blosum62),
            };
            let hsp = TabularHsp {
                q_start: 1,
                q_end: 6,
                s_start: 1,
                s_end: 6,
                evalue: 1e-5,
                bit_score: 20.0,
                score: 10,
                length: 6,
                gap_opens: 0,
                gaps: 0,
                num_ident,
                q_frame: 1,
                s_frame: 1,
                rows,
                end_gaps: true,
                qcovhsp: 100,
            };
            let mut out = Vec::new();
            writer.row(
                &mut out,
                &mut TabularRows::default(),
                &TabularQuery::local(b"q", 6),
                &subject,
                None,
                &hsp,
            );
            String::from_utf8(out).unwrap()
        };
        let frames = spec(ReportFormat::Tabular, "qseqid frames sframe");
        assert_eq!(row(&frames, 3, None), "q\t0/0\t0\n");
        let counted = spec(ReportFormat::Tabular, "nident frames sframe");
        assert_eq!(row(&counted, 0, None), "0\t0/0\t0\n");
        assert_eq!(row(&counted, 3, None), "3\t1/1\t1\n");
        let positives = spec(ReportFormat::Tabular, "nident positive mismatch sstrand");
        assert_eq!(
            row(&positives, 1, Some((b"UCJ*KA", b"CUJ*IA"))),
            "3\t5\t3\tN/A\n"
        );
        assert_eq!(row(&positives, 1, Some((b"UJ", b"UI"))), "1\t1\t5\tN/A\n");
        let commented = spec(ReportFormat::TabularWithComments, "qseqid");
        let writer = TabularWriter {
            spec: &commented,
            program_version: "PSIBLAST 2.17.0+",
            matrix: Some(ScoringMatrix::Blosum62),
        };
        let mut out = Vec::new();
        writer
            .round_query_header(&mut out, Some(1), b"q1 title", b"db", Some(0))
            .unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "# PSIBLAST 2.17.0+\n# Iteration: 1\n# Query: q1 title\n# Database: db\n# 0 hits found\n"
        );
    }

    // NCBI tabular.cpp:219-255 (`sseqid`, `sacc`, `saccver`, `sgi`: the subject's own ids) and
    // :387-439 (`salltitles` joined by `<>`, `stitle` the first title, an empty one `N/A`); a local
    // query's id is its `qseqid`, `qacc` and `qaccver` alike.
    #[test]
    fn id_columns_are_the_subject_s_own() {
        let spec = TabularSpec::of(
            ReportFormat::Tabular,
            &FormatChoice {
                number: 6,
                spec: "qseqid qacc qaccver sseqid sacc saccver sgi stitle salltitles".to_string(),
                delimiter: String::new(),
            },
        )
        .unwrap();
        let taxonomy = TabularTaxonomy::default();
        let titles: [&[u8]; 2] = [b"first title", b""];
        let subject = TabularSubject {
            id: b"gb|AE000111.1|",
            accession: b"AE000111",
            accession_version: b"AE000111.1",
            gi: 42,
            all_ids: &[],
            length: 10,
            titles: Some(&titles),
            taxonomy: &taxonomy,
            names: None,
            nucleotide: false,
            qcovs: -1,
            qcovus: -1,
            lookup: SubjectLookup::Found,
            index: 0,
        };
        let hsp = TabularHsp {
            q_start: 1,
            q_end: 6,
            s_start: 1,
            s_end: 6,
            evalue: 1e-5,
            bit_score: 20.0,
            score: 10,
            length: 6,
            gap_opens: 0,
            gaps: 0,
            num_ident: 6,
            q_frame: 1,
            s_frame: 1,
            rows: None,
            end_gaps: true,
            qcovhsp: 100,
        };
        let mut out = Vec::new();
        TabularWriter {
            spec: &spec,
            program_version: "BLASTP 2.17.0+",
            matrix: None,
        }
        .row(
            &mut out,
            &mut TabularRows::default(),
            &TabularQuery::local(b"q1", 6),
            &subject,
            None,
            &hsp,
        );
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "q1\tq1\tq1\tgb|AE000111.1|\tAE000111\tAE000111.1\t42\tfirst title\tfirst title<>N/A\n"
        );
    }

    /// The rows of `subjects` (failed or found, in order) of one query under `fields`, with a
    /// nucleotide HSP of 450 letters (qcovs 99, qcovhsp 50) and query length 900.
    fn failed_rows(fields: &str, subjects: &[TabularSubject<'_>]) -> String {
        let spec = TabularSpec::of(
            ReportFormat::Tabular,
            &FormatChoice {
                number: 6,
                spec: fields.to_string(),
                delimiter: String::new(),
            },
        )
        .unwrap();
        let hsp = TabularHsp {
            q_start: 1,
            q_end: 450,
            s_start: 1,
            s_end: 450,
            evalue: 0.0,
            bit_score: 812.0,
            score: 900,
            length: 450,
            gap_opens: 0,
            gaps: 0,
            num_ident: 450,
            q_frame: 1,
            s_frame: 1,
            rows: None,
            end_gaps: true,
            qcovhsp: 50,
        };
        let writer = TabularWriter {
            spec: &spec,
            program_version: "BLASTN 2.17.0+",
            matrix: None,
        };
        let mut rows = TabularRows::default();
        let mut out = Vec::new();
        for subject in subjects {
            // The kept subject of a failed row: the last found one, by its index.
            let kept = rows
                .last_found()
                .and_then(|index| subjects.iter().find(|subject| subject.index == index));
            writer.row(
                &mut out,
                &mut rows,
                &TabularQuery::local(b"nq1", 900),
                subject,
                kept.filter(|_| writer.fails(subject)),
                &hsp,
            );
        }
        String::from_utf8(out).unwrap()
    }

    // NCBI tabular.cpp:828-875 and 1374-1378 (the catch of `SetFields` and `SetSubjectId(list)`):
    // NCBI BLAST+ 2.17.0 `blastn -subject` with two titles `Subject_9 bad x\x81y title`,
    // `Subject_7 other \x8dz` (probe `remaining-gaps-20261112/T/t7`, `s.bb.blastn.f6wide`), and a
    // found subject before a failed one (`s.gb`): `sseqid` Unknown, `sgi` 0, the Seq-align ids added
    // to `sall*` row by row, the counts of `x_ResetFields`, `slen` 0, the scores, the query's ids
    // and coverages real, the titles and taxonomy of the last found row.
    #[test]
    fn failed_lookups_write_ncbis_catch_row() {
        let taxonomy = TabularTaxonomy::default();
        let subject = |id: &'static [u8], failed: Option<&'static [u8]>, index| TabularSubject {
            id,
            accession: id,
            accession_version: id,
            gi: 0,
            all_ids: &[],
            length: 1350,
            titles: None,
            taxonomy: &taxonomy,
            names: None,
            nucleotide: true,
            qcovs: 99,
            qcovus: -1,
            lookup: match failed {
                Some(align_id) => SubjectLookup::Failed(FailedLookup {
                    full: align_id,
                    accession: align_id,
                    subject_id: true,
                    all_ids: None,
                }),
                None => SubjectLookup::Found,
            },
            index,
        };
        let wide = "qseqid sseqid evalue bitscore pident length mismatch qstart sstart slen qlen \
                    nident sgi sallseqid sallgi sallacc stitle salltitles staxids sscinames \
                    sstrand qcovs qcovhsp";
        let bad = |n: &'static [u8]| subject(b"", Some(n), 9);
        let fallback = "\t0.0\t812\t0.000\t0\t0\t0\t0\t0\t900\t0\t0\t";
        let tail = "\tN/A\tN/A\tN/A\tN/A\tN/A\t99\t50\n";
        assert_eq!(
            failed_rows(
                wide,
                &[bad(b"Subject_1"), bad(b"Subject_1"), bad(b"Subject_2")]
            ),
            format!(
                "nq1\tUnknown{fallback}Subject_1\t0\tSubject_1{tail}\
                 nq1\tUnknown{fallback}Subject_1;Subject_1\t0;0\tSubject_1;Subject_1{tail}\
                 nq1\tUnknown{fallback}Subject_1;Subject_1;Subject_2\t0;0;0\t\
                 Subject_1;Subject_1;Subject_2{tail}"
            )
        );
        let rows = failed_rows(wide, &[subject(b"s2", None, 1), bad(b"Subject_2")]);
        assert_eq!(
            rows.lines().nth(1).unwrap(),
            format!(
                "nq1\tUnknown{fallback}s2;Subject_2\t0;0\ts2;Subject_2{}",
                tail.trim_end()
            )
        );
        // Without `sseqid`, `sgi`, `sacc` and `saccver` a `-subject` record's lookup does not run
        // (`x_SetSubjectIds` of a record has no `GetSeqIdList`): the row is the found one's.
        assert_eq!(
            failed_rows(
                "qseqid sallseqid evalue",
                &[subject(b"Subject_9", Some(b"Subject_1"), 0)]
            ),
            "nq1\tSubject_9\t0.0\n"
        );
    }

    // NCBI tabular.cpp:828-875 (the catch of the subject block, `SetSubjectId(list)`, `return -1`)
    // where the load of the subject throws (`DB contains no accession info.`, a batch whose
    // prefetch is off; probe of line DBID pfoff-2, `R/pfoff/p2`, `bp.U2.*`, DECISIONS 62 (12)):
    // the row keeps the last loaded row's ids, titles and taxonomy, `sall*` add the Seq-align id
    // (its GI in `sallgi`), `slen` and the counts are `x_ResetFields`' 0, `qlen` comes only with
    // a field of the query block; before any loaded row `sseqid` is `Unknown`, the titles and
    // names `N/A` (quoted in CSV).
    #[test]
    fn a_load_that_throws_writes_the_kept_row() {
        let exception = crate::report::db_subject::LookupFailure::ncbi(
            crate::seqdb::tax_files::NO_ACCESSION_INFO,
        )
        .exception()
        .unwrap();
        let taxonomy = TabularTaxonomy {
            taxid: 9606,
            taxids: vec![9606],
        };
        let none = TabularTaxonomy::default();
        let titles: [&[u8]; 1] = [b"found title"];
        let all_ids = [DeflineIds {
            id: b"gb|AAA00001.1|".to_vec(),
            gi: 7,
            accession: b"AAA00001".to_vec(),
        }];
        let found = TabularSubject {
            id: b"gb|AAA00001.1|",
            accession: b"AAA00001",
            accession_version: b"AAA00001.1",
            gi: 7,
            all_ids: &all_ids,
            length: 1350,
            titles: Some(&titles),
            taxonomy: &taxonomy,
            names: None,
            nucleotide: true,
            qcovs: 99,
            qcovus: -1,
            lookup: SubjectLookup::Found,
            index: 0,
        };
        let thrown = TabularSubject {
            id: b"",
            accession: b"",
            accession_version: b"",
            gi: 0,
            all_ids: &[],
            length: 0,
            titles: None,
            taxonomy: &none,
            names: None,
            nucleotide: true,
            qcovs: 98,
            qcovus: -1,
            lookup: SubjectLookup::Throws(ThrownLookup {
                full: b"gi|5|gb|BBB00002.1|",
                gi: 5,
                accession: b"BBB00002",
                exception: &exception,
            }),
            index: 1,
        };
        assert!(thrown.throws() && !found.throws());
        assert_eq!(
            failed_rows(
                "qseqid sseqid sacc saccver sgi slen qlen length pident evalue stitle staxids \
                 sstrand qcovs",
                &[found, thrown]
            )
            .lines()
            .nth(1)
            .unwrap(),
            "nq1\tgb|AAA00001.1|\tAAA00001\tAAA00001.1\t7\t0\t900\t0\t0.000\t0.0\tfound title\t9606\tN/A\t98"
        );
        // `qlen` without a field of the query block keeps `x_ResetFields`' 0.
        assert_eq!(
            failed_rows("sseqid qlen slen evalue", &[found, thrown])
                .lines()
                .nth(1)
                .unwrap(),
            "gb|AAA00001.1|\t0\t0\t0.0"
        );
        // `sall*`: the kept lists and the Seq-align id, row by row.
        assert_eq!(
            failed_rows("qseqid sallseqid sallgi sallacc", &[found, thrown, thrown]),
            "nq1\tgb|AAA00001.1|\t7\tAAA00001\n\
             nq1\tgb|AAA00001.1|;gi|5|gb|BBB00002.1|\t7;5\tAAA00001;BBB00002\n\
             nq1\tgb|AAA00001.1|;gi|5|gb|BBB00002.1|;gi|5|gb|BBB00002.1|\t7;5;5\t\
             AAA00001;BBB00002;BBB00002\n"
        );
        // Before any loaded row: `Unknown`, `sgi` 0, the titles `N/A`.
        assert_eq!(
            failed_rows("qseqid sseqid sgi stitle sallseqid evalue", &[thrown]),
            "nq1\tUnknown\t0\tN/A\tgi|5|gb|BBB00002.1|\t0.0\n"
        );
        // The refusals before the report: a first row with `sacc`, `saccver` (the default fields)
        // or `staxid`; a later one has a result.
        let spec = |fields: &str| {
            TabularSpec::of(
                ReportFormat::Tabular,
                &FormatChoice {
                    number: 6,
                    spec: fields.to_string(),
                    delimiter: String::new(),
                },
            )
            .unwrap()
        };
        let thrown = LookupFailure::thrown();
        assert!(failed_lookup_stop(None, thrown, true).is_some());
        assert!(failed_lookup_stop(None, thrown, false).is_none());
        assert!(failed_lookup_stop(Some(&spec("qseqid saccver")), thrown, true).is_some());
        assert!(failed_lookup_stop(Some(&spec("qseqid saccver")), thrown, false).is_none());
        assert!(failed_lookup_stop(Some(&spec("staxid")), thrown, true).is_some());
        assert!(failed_lookup_stop(Some(&spec("qseqid sseqid")), thrown, true).is_none());
        assert!(may_stop_at_failed_lookup(Some(&spec("staxid ssciname"))));
    }

    // NCBI tabular.cpp:828-875 (the catch of the subject block) for the default fields (probe of
    // line DBID pfoff-2, `R/pfoff/p2`, `*.U2.dflt.6`; review 2, M3): the one writer of the
    // default tables' rows of a subject whose load throws (BLASTP's, BLASTN's, BLASTX's, TBLASTN's
    // and TBLASTX's, `TabularRows::write_thrown_default`) prints `saccver` of the query's last
    // loaded row (`TabularRows::found`), the counts and coordinates 0 and the scores; a new
    // query's formatter keeps none. The end of the run (`run_end`) writes the rows before it,
    // ends the reads and gives NCBI's exception (exit status 2).
    #[test]
    fn the_default_row_of_a_throw_keeps_the_last_loaded_subject() {
        let config = crate::report::outfmt6::OutputConfig::ncbi_compat();
        let saccver = |index: usize| ["AAA00001.1", "BBB00002.1"][index];
        let mut rows = TabularRows::default();
        assert!(rows.kept_id().is_none());
        rows.found(1);
        assert_eq!(rows.kept_id(), Some(1));
        let mut out = Vec::new();
        rows.write_thrown_default(&mut out, b"q1", saccver, 1e-30, 55.4, &config)
            .unwrap();
        rows.found(0);
        rows.write_thrown_default(&mut out, b"q1", saccver, 0.5, 120.7, &config)
            .unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "q1\tBBB00002.1\t0.000\t0\t0\t0\t0\t0\t0\t0\t1.00e-30\t55.4\n\
             q1\tAAA00001.1\t0.000\t0\t0\t0\t0\t0\t0\t0\t0.50\t120\n"
        );
        let exception = crate::report::db_subject::LookupFailure::ncbi(
            crate::seqdb::tax_files::NO_ACCESSION_INFO,
        )
        .exception()
        .unwrap();
        let (mut written, mut ended) = (Vec::new(), false);
        let error = run_end(&mut written, b"rows before\n", &exception, || ended = true);
        assert_eq!(written, b"rows before\n");
        assert!(ended);
        let native = error.downcast::<crate::cli::NativeError>().unwrap();
        assert_eq!(native.exit, 2);
        assert!(
            native.message.contains("DB contains no accession info."),
            "{}",
            native.message
        );
    }

    // NCBI tabular.cpp:833-898: which rows end the run where the load throws (`slen`, the
    // letters after the subject block), which are written (a field of that block), and which
    // read nothing; `nident` alone reads only on rows with identities.
    #[test]
    fn a_rows_fields_say_what_its_throw_does() {
        use crate::blastinput::app::Program;
        use crate::report::db_subject::{AfterThrow, FormatThrow, RunThrow};
        let spec = |fields: &str| {
            TabularSpec::of(
                ReportFormat::CommaSeparatedValues,
                &FormatChoice {
                    number: 10,
                    spec: fields.to_string(),
                    delimiter: String::new(),
                },
            )
            .unwrap()
        };
        let row = |fields: &str, program, ident| spec(fields).row_throw(program, ident);
        assert_eq!(
            row("qseqid stitle slen", Program::Blastp, 1),
            Some(AfterThrow::WritesRow)
        );
        assert_eq!(
            row("qseqid slen", Program::Blastn, 1),
            Some(AfterThrow::EndsRun)
        );
        assert_eq!(
            row("qseqid sseq", Program::Tblastx, 1),
            Some(AfterThrow::EndsRun)
        );
        assert_eq!(
            row("qseqid nident", Program::Blastx, 1),
            Some(AfterThrow::EndsRun)
        );
        assert_eq!(row("qseqid nident", Program::Blastx, 0), None);
        assert_eq!(row("qseqid nident", Program::Blastn, 1), None);
        assert_eq!(row("qseqid evalue qcovs", Program::Tblastn, 1), None);
        assert_eq!(
            spec("qseqid nident").format_throw(Program::Blastx),
            FormatThrow::Writes(RunThrow {
                after: AfterThrow::EndsRun,
                rows: SubjectFetch::WithIdentities,
            })
        );
        // `slen` and a field of the alignment block read the subject at every row.
        for every in ["qseqid nident slen", "qseqid nident qstart"] {
            assert_eq!(spec(every).throw_rows(Program::Blastx), SubjectFetch::Every);
            assert_eq!(
                spec(every).format_throw(Program::Blastx),
                FormatThrow::Writes(RunThrow::every_row(AfterThrow::EndsRun)),
                "{every}"
            );
            assert_eq!(
                spec(every).row_throw(Program::Blastx, 0),
                Some(AfterThrow::EndsRun)
            );
        }
        // The organism report (outfmt 18) leaves out the HSPs of a throwing subject and goes on.
        let report = TabularSpec::of(
            ReportFormat::TaxFormat,
            &FormatChoice {
                number: 18,
                spec: String::new(),
                delimiter: String::new(),
            },
        )
        .unwrap();
        assert_eq!(
            report.format_throw(Program::Blastp),
            FormatThrow::Writes(RunThrow::every_row(AfterThrow::WritesRow))
        );
        assert_eq!(
            spec("qseqid evalue").format_throw(Program::Blastx),
            FormatThrow::Unread
        );
    }

    // NCBI BLAST+ 2.17.0 `blastn -db` of a sequence with an empty title and a taxid whose
    // name has 0x81 after one of Homo sapiens (probe `t.tx_bh.blastn.f6all|f6wide`): without
    // `sseqid` the lookup fails in `x_SetSubjectIds`, which clears `m_SubjectIds` first, so the
    // Seq-align id is the list alone and `slen` is set; the taxonomy is the last found row's.
    #[test]
    fn a_database_lookup_failing_in_x_set_subject_ids_keeps_the_length() {
        let homo = TabularTaxonomy {
            taxid: 9606,
            taxids: vec![9606],
        };
        let none = TabularTaxonomy::default();
        let titles: [&[u8]; 1] = [b""];
        let subject = |id: &'static [u8], taxonomy, failed: bool, index| TabularSubject {
            id,
            accession: id,
            accession_version: id,
            gi: 0,
            all_ids: &[],
            length: 660,
            titles: Some(&titles),
            taxonomy,
            names: None,
            nucleotide: true,
            qcovs: 100,
            qcovus: -1,
            lookup: if failed {
                SubjectLookup::Failed(FailedLookup {
                    full: b"gnl|BL_ORD_ID|0",
                    accession: b"BL_ORD_ID:0",
                    subject_id: true,
                    all_ids: Some(&[]),
                })
            } else {
                SubjectLookup::Found
            },
            index,
        };
        let found = subject(b"Homo", &homo, false, 1);
        let failed = subject(b"", &none, true, 0);
        assert_eq!(
            failed_rows("qseqid sallseqid sallacc slen staxids", &[found, failed]),
            "nq1\tHomo\tHomo\t660\t9606\nnq1\tgnl|BL_ORD_ID|0\tBL_ORD_ID:0\t660\t9606\n"
        );
        assert_eq!(
            failed_rows("qseqid sseqid slen staxids", &[found, failed]),
            "nq1\tHomo\t660\t9606\nnq1\tUnknown\t0\t9606\n"
        );
    }

    // NCBI BLAST+ 2.17.0 `blastn -db` of a `-parse_seqids` sequence with the definition lines
    // `gb|AB000001.1|` and `lcl|Subject_9`, an empty title and a taxid whose name has 0x81 (probe
    // `T/t7/review` `B.gbs9.blastn.f6wide`): `SetSubjectId` finds the Bioseq's `gb` id, and
    // `x_SetSubjectIds` throws at the second line after it pushed the first one's ids, so `sall*`
    // are the first line's and the Seq-align id; `slen` is set, the counts are `x_ResetFields`'.
    #[test]
    fn a_lookup_failing_after_the_first_definition_line_keeps_its_ids() {
        let none = TabularTaxonomy::default();
        let titles: [&[u8]; 2] = [b"", b""];
        let before = [DeflineIds {
            id: b"gb|AB000001.1|".to_vec(),
            gi: 0,
            accession: b"AB000001".to_vec(),
        }];
        let subject = TabularSubject {
            id: b"gb|AB000001.1|",
            accession: b"AB000001",
            accession_version: b"AB000001.1",
            gi: 0,
            all_ids: &before,
            length: 1350,
            titles: Some(&titles),
            taxonomy: &none,
            names: None,
            nucleotide: true,
            qcovs: 99,
            qcovus: -1,
            lookup: SubjectLookup::Failed(FailedLookup {
                full: b"gb|AB000001.1|",
                accession: b"AB000001",
                subject_id: false,
                all_ids: Some(&before),
            }),
            index: 0,
        };
        assert_eq!(
            failed_rows("qseqid sseqid slen sallseqid sallgi sallacc evalue", &[subject; 2]),
            "nq1\tgb|AB000001.1|\t1350\tgb|AB000001.1|;gb|AB000001.1|\t0;0\tAB000001;AB000001\t0.0\n"
                .repeat(2)
        );
        // `saccver` is the Bioseq's: no stop where only `x_SetSubjectIds` throws.
        assert_eq!(
            failed_rows("qseqid saccver sallseqid", &[subject]),
            "nq1\tAB000001.1\tgb|AB000001.1|;gb|AB000001.1|\n"
        );
    }

    // NCBI tabular.cpp:174-185 (`FindBestChoice` of the empty `m_SubjectId`, a null pointer:
    // NCBI BLAST+ 2.17.0 exits 255) and :262-268 (`staxid` of an `m_SubjectTaxId` that no row set:
    // a different number in each run, probe of line RG T7 `6 qseqid sseqid staxid`).
    #[test]
    fn failed_lookups_stop_where_ncbi_has_no_result() {
        let spec = |fields: &str| {
            TabularSpec::of(
                ReportFormat::CommaSeparatedValues,
                &FormatChoice {
                    number: 10,
                    spec: fields.to_string(),
                    delimiter: String::new(),
                },
            )
            .unwrap()
        };
        let stop = |fields: Option<&str>, all_ids: bool, first: bool| {
            let spec = fields.map(spec);
            let failure = LookupFailure {
                subject_id: true,
                all_ids,
                bioseq: false,
            };
            failed_lookup_stop(spec.as_ref(), failure, first).is_some()
        };
        assert!(stop(None, false, false));
        assert!(stop(Some(""), false, false));
        assert!(stop(Some("qseqid sacc"), false, false));
        assert!(stop(Some("qseqid saccver"), true, false));
        assert!(!stop(Some("qseqid sseqid evalue"), false, true));
        assert!(stop(Some("qseqid sseqid staxid"), false, true));
        assert!(!stop(Some("qseqid sseqid staxid"), false, false));
        assert!(!stop(Some("qseqid sallseqid staxid"), false, true));
        assert!(stop(Some("qseqid sallseqid staxid"), true, true));
        // `SetSubjectId` finds the ids where only `x_SetSubjectIds` throws: no null pointer.
        let only_all = LookupFailure {
            subject_id: false,
            all_ids: true,
            bioseq: false,
        };
        assert!(
            failed_lookup_stop(Some(&spec("qseqid saccver sallseqid")), only_all, false).is_none()
        );
        assert!(failed_lookup_stop(None, only_all, true).is_none());
        assert!(!may_stop_at_failed_lookup(Some(&spec(
            "qseqid sseqid sallseqid evalue"
        ))));
        assert!(may_stop_at_failed_lookup(None));
        let record = |s_idx| (s_idx == 1).then(|| LookupFailure::of_plain(false));
        let result = check_failed_lookups(
            Some(&spec("qseqid sseqid staxid")),
            [vec![0, 1], vec![1, 0]],
            record,
            |s_idx| failed_record_title(s_idx, b"Subject_9"),
            "BLASTN",
        );
        assert!(result
            .unwrap_err()
            .to_string()
            .starts_with("subject record 2 (Subject_9) has a title"));
        assert!(check_failed_lookups(
            Some(&spec("qseqid sseqid staxid")),
            [vec![0, 1], vec![0, 1]],
            record,
            |s_idx| failed_record_title(s_idx, b"x"),
            "BLASTN",
        )
        .is_ok());
    }

    #[test]
    fn btop_counts_runs() {
        let mut out = Vec::new();
        write_btop(&mut out, b"ACGT-A", b"ACCTTA", true);
        assert_eq!(out, b"2GC1-T1");
    }

    // tabular.cpp:982-1027 (quoted at `push_alignment_row` and `write_btop`): a row's end gaps
    // are spaces, its inner gaps `-`; an end-gap column of `btop` is a letter and a space.
    #[test]
    fn end_gaps_are_spaces_in_the_rows_and_btop() {
        let mut out = Vec::new();
        push_alignment_row(&mut out, b"---AC-DE--", true);
        assert_eq!(out, b"   AC-DE  ");
        out.clear();
        push_alignment_row(&mut out, b"MKV-L", true);
        assert_eq!(out, b"MKV-L");
        out.clear();
        write_btop(&mut out, b"MQLAC-D", b"---ACKD", true);
        assert_eq!(out, b"M Q L 2-K1");
        out.clear();
        write_btop(&mut out, b"AC-DKK", b"ACGD--", true);
        assert_eq!(out, b"2-G1K K ");
        // Rows whose `-` fill a segment (BLASTX's mapped row past a query's end, db.x2g.bx.s1o.subj.6w
        // and .10: NCBI's `TCGGCATGGCTG-` and BTOP `...-X`).
        out.clear();
        push_alignment_row(&mut out, b"TCG-", false);
        assert_eq!(out, b"TCG-");
        out.clear();
        write_btop(&mut out, b"TCG-", b"TANX", false);
        assert_eq!(out, b"1CAGN-X");
    }
}
