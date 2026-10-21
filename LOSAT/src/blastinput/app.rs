//! NCBI's application layer shared by the BLASTP, TBLASTN and TBLASTX command lines
//! (`CBlastAppArgs::SetOptions` and the `CATCH_ALL` of the applications): the error
//! frames, the `-outfmt` string (read so by all five programs), the `-seg` and
//! `-comp_based_stats` strings that NCBI reads in its option handlers, after the input and
//! output files are opened.

use std::path::Path;

use crate::blastinput::value_parsers::SegSpec;
use crate::cli::NativeError;

/// An error that NCBI raises as `CInputException` (an option handler or `Validate`).
///
/// NCBI reference: c++/src/app/blast/blast_app_util.hpp:172-176
/// ```c
///     catch (const blast::CInputException& e) {                               \
///         LOG_POST(Error << "BLAST query/options error: " << e.GetMsg());     \
///         LOG_POST(Error << "Please refer to the BLAST+ user manual.");       \
///         exit_code = BLAST_INPUT_ERROR;                                      \
///     }                                                                       \
/// ```
pub fn options_error(message: &str) -> anyhow::Error {
    NativeError {
        exit: 1,
        message: format!(
            "BLAST query/options error: {message}\nPlease refer to the BLAST+ user manual.\n"
        ),
    }
    .into()
}

/// An error that NCBI raises as a `CBlastException` that is neither an option nor a
/// memory error (exit 3).
///
/// NCBI reference: c++/src/app/blast/blast_app_util.hpp:225-227
/// ```c
///         } else {                                                            \
///             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
///             exit_code = BLAST_ENGINE_ERROR;                                 \
/// ```
pub fn engine_error(message: &str) -> anyhow::Error {
    NativeError {
        exit: 3,
        message: format!("BLAST engine error: {message}\n"),
    }
    .into()
}

/// NCBI's error for a search without `-subject` (and without `-db`), which the handler of
/// the database arguments raises before the query and the output are opened.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2558-2562
/// ```c
///     } else if (!m_IsIgBlast){
///         // IgBlast permits use of germline database
///         NCBI_THROW(CInputException, eInvalidInput,
///            "Either a BLAST database or subject sequence(s) must be specified");
///     }
/// ```
pub fn missing_subject_error() -> anyhow::Error {
    options_error("Either a BLAST database or subject sequence(s) must be specified")
}

/// NCBI's error for a subject file without records, raised when the subjects are read.
///
/// NCBI reference: c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
/// ```c
/// CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
///     : m_QueryVector(& queries)
/// {
///     if (queries.Empty()) {
///         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
///     }
/// ```
pub fn empty_subjects_error() -> anyhow::Error {
    engine_error("Empty CBlastQueryVector")
}

/// NCBI's output format numbers (`CFormattingArgs::EOutputFormat`).
///
/// NCBI reference: c++/include/algo/blast/blastinput/blast_args.hpp:1023-1074
/// ```c
///     /// Defines the output formats supported by our command line formatter
///     enum EOutputFormat {
///         /// Standard pairwise alignments
///         ePairwise = 0,
///         ///< Query anchored showing identities
///         eQueryAnchoredIdentities,
///         ///< Query anchored no identities
///         eQueryAnchoredNoIdentities,
///         ///< Flat query anchored showing identities
///         eFlatQueryAnchoredIdentities,
///         ///< Flat query anchored no identities
///         eFlatQueryAnchoredNoIdentities,
///         /// XML output
///         eXml,
///         /// Tabular output
///         eTabular,
///         /// Tabular output with comments
///         eTabularWithComments,
///         /// ASN.1 text output
///         eAsnText,
///         /// ASN.1 binary output
///         eAsnBinary,
///         /// Comma-separated values
///         eCommaSeparatedValues,
///         /// BLAST archive format
///         eArchiveFormat,
///         /// JSON seq-align
///         eJsonSeqalign,
///         /// JSON XInclude
///         eJson,
///         /// XML2 XInclude
///         eXml2,
///         /// JSON2 single file
///         eJson_S,
///         /// XML2 single file
///         eXml2_S,
///         /// SAM format
///         eSAM,
///
///         eTaxFormat,
///
///         ///igblast AIRR rearrangement, 19
///         eAirrRearrangement,
///         /// Comma-separated values with a header
///         eCommaSeparatedValuesWithHeader,
///
///         /// unaligned reads in magicblast
///         eFasta,
///         /// Sentinel value for error checking
///         eEndValue
///
///     };
/// ```
pub mod output_format {
    pub const PAIRWISE: i32 = 0;
    pub const QUERY_ANCHORED_IDENTITIES: i32 = 1;
    pub const QUERY_ANCHORED_NO_IDENTITIES: i32 = 2;
    pub const FLAT_QUERY_ANCHORED_IDENTITIES: i32 = 3;
    pub const FLAT_QUERY_ANCHORED_NO_IDENTITIES: i32 = 4;
    pub const XML: i32 = 5;
    pub const TABULAR: i32 = 6;
    pub const TABULAR_WITH_COMMENTS: i32 = 7;
    pub const ASN_TEXT: i32 = 8;
    pub const ASN_BINARY: i32 = 9;
    pub const COMMA_SEPARATED_VALUES: i32 = 10;
    pub const ARCHIVE_FORMAT: i32 = 11;
    pub const JSON_SEQALIGN: i32 = 12;
    pub const JSON: i32 = 13;
    pub const XML2: i32 = 14;
    pub const JSON_S: i32 = 15;
    pub const XML2_S: i32 = 16;
    pub const SAM: i32 = 17;
    pub const TAX_FORMAT: i32 = 18;
    pub const AIRR_REARRANGEMENT: i32 = 19;
    pub const COMMA_SEPARATED_VALUES_WITH_HEADER: i32 = 20;
    pub const FASTA: i32 = 21;
    pub const END_VALUE: i32 = 22;
}

use output_format::*;

/// A `-outfmt` value as NCBI's `ParseFormattingString` reads it: the format number and the
/// custom specification (fields and a delimiter) that remains for it (only the tabular,
/// comma-separated and SAM formats keep one).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FormatChoice {
    pub number: i32,
    /// The custom fields (`custom_fmt_spec`), without a `delim=` token.
    pub spec: String,
    /// The custom delimiter (`custom_delim`), empty when none is given.
    pub delimiter: String,
}

impl FormatChoice {
    /// Whether a custom specification (fields or a delimiter) remains.
    pub fn custom(&self) -> bool {
        !self.spec.is_empty() || !self.delimiter.is_empty()
    }

    /// The format as LOSAT's reports read it: the number and the custom fields.
    pub fn normalized(&self) -> String {
        if self.spec.is_empty() {
            self.number.to_string()
        } else {
            format!("{} {}", self.number, self.spec)
        }
    }
}

/// The reports of NCBI's search programs: every `EOutputFormat` that blastn, blastp, blastx,
/// tblastn or tblastx writes (blast_args.hpp:1023-1074 above). The AIRR rearrangement (19)
/// and FASTA (21) formats are igblastn's and magicblast's (`formatting_handler_check`).
/// Which of them LOSAT writes is `Program::writes`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum ReportFormat {
    Pairwise = PAIRWISE,
    QueryAnchoredIdentities = QUERY_ANCHORED_IDENTITIES,
    QueryAnchoredNoIdentities = QUERY_ANCHORED_NO_IDENTITIES,
    FlatQueryAnchoredIdentities = FLAT_QUERY_ANCHORED_IDENTITIES,
    FlatQueryAnchoredNoIdentities = FLAT_QUERY_ANCHORED_NO_IDENTITIES,
    Xml = XML,
    Tabular = TABULAR,
    TabularWithComments = TABULAR_WITH_COMMENTS,
    AsnText = ASN_TEXT,
    AsnBinary = ASN_BINARY,
    CommaSeparatedValues = COMMA_SEPARATED_VALUES,
    ArchiveFormat = ARCHIVE_FORMAT,
    JsonSeqalign = JSON_SEQALIGN,
    Json = JSON,
    Xml2 = XML2,
    JsonS = JSON_S,
    Xml2S = XML2_S,
    Sam = SAM,
    TaxFormat = TAX_FORMAT,
    CommaSeparatedValuesWithHeader = COMMA_SEPARATED_VALUES_WITH_HEADER,
}

impl ReportFormat {
    /// The report of NCBI's format number; `None` for 19 and 21 (the formats of other
    /// programs) and for a number outside `EOutputFormat`.
    pub fn from_number(number: i32) -> Option<Self> {
        use ReportFormat::*;
        Some(match number {
            PAIRWISE => Pairwise,
            QUERY_ANCHORED_IDENTITIES => QueryAnchoredIdentities,
            QUERY_ANCHORED_NO_IDENTITIES => QueryAnchoredNoIdentities,
            FLAT_QUERY_ANCHORED_IDENTITIES => FlatQueryAnchoredIdentities,
            FLAT_QUERY_ANCHORED_NO_IDENTITIES => FlatQueryAnchoredNoIdentities,
            XML => Xml,
            TABULAR => Tabular,
            TABULAR_WITH_COMMENTS => TabularWithComments,
            ASN_TEXT => AsnText,
            ASN_BINARY => AsnBinary,
            COMMA_SEPARATED_VALUES => CommaSeparatedValues,
            ARCHIVE_FORMAT => ArchiveFormat,
            JSON_SEQALIGN => JsonSeqalign,
            JSON => Json,
            XML2 => Xml2,
            JSON_S => JsonS,
            XML2_S => Xml2S,
            SAM => Sam,
            TAX_FORMAT => TaxFormat,
            COMMA_SEPARATED_VALUES_WITH_HEADER => CommaSeparatedValuesWithHeader,
            _ => return None,
        })
    }

    /// NCBI's number of the format (`-outfmt`).
    pub fn number(self) -> i32 {
        self as i32
    }

    /// Whether NCBI writes a query's result warnings (an invalid query, O replaced by X,
    /// ...) with the query's report: not for the structured formats, whose report returns
    /// before them. (The reader's warnings of a query batch are written when it is read,
    /// whatever the format.)
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1431-1452
    /// ```c++
    ///     if (m_FormatType == CFormattingArgs::eAsnText
    ///       || m_FormatType == CFormattingArgs::eAsnBinary
    ///       || m_FormatType == CFormattingArgs::eXml
    ///       || m_FormatType == CFormattingArgs::eXml2
    ///       || m_FormatType == CFormattingArgs::eJson
    ///       || m_FormatType == CFormattingArgs::eXml2_S
    ///       || m_FormatType == CFormattingArgs::eJson_S
    ///       || m_FormatType == CFormattingArgs::eJsonSeqalign
    ///       || m_FormatType == CFormattingArgs::eSAM)
    ///     {
    ///         x_PrintStructuredReport(results, queries);
    ///         return;
    ///     }
    ///
    ///     if (results.HasErrors()) {
    ///         ERR_POST(Error << results.GetErrorStrings());
    ///         return; // errors are deemed fatal
    ///     }
    ///     if (results.HasWarnings()) {
    ///         ERR_POST(Warning << results.GetWarningStrings());
    ///     }
    /// ```
    pub fn writes_query_warnings(self) -> bool {
        use ReportFormat::*;
        !matches!(
            self,
            AsnText | AsnBinary | Xml | Xml2 | Json | Xml2S | JsonS | JsonSeqalign | Sam
        )
    }
}

/// LOSAT's search programs, for what depends on the program in the reading of `-outfmt`.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blastn_args.cpp:104
/// ```c
///     m_FormattingArgs.Reset(new CFormattingArgs(false, CFormattingArgs::eIsSAM));
/// ```
/// NCBI reference: c++/src/algo/blast/blastinput/blastp_args.cpp:95 (blastx_args.cpp:113,
/// tblastn_args.cpp:110 and tblastx_args.cpp:104 are the same)
/// ```c
///     m_FormattingArgs.Reset(new CFormattingArgs);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Program {
    Blastn,
    Blastp,
    Blastx,
    Tblastn,
    Tblastx,
}

impl Program {
    /// The name in LOSAT's rejections ("... is not supported by LOSAT's BLASTN").
    pub fn name(self) -> &'static str {
        match self {
            Program::Blastn => "BLASTN",
            Program::Blastp => "BLASTP",
            Program::Blastx => "BLASTX",
            Program::Tblastn => "TBLASTN",
            Program::Tblastx => "TBLASTX",
        }
    }

    /// Whether NCBI's formatting arguments of the program have `eIsSAM` (blastn only).
    pub fn is_sam(self) -> bool {
        self == Program::Blastn
    }

    /// Whether LOSAT's program writes custom tabular field lists (BLASTP and BLASTX check
    /// their fields themselves).
    pub fn custom_fields(self) -> bool {
        matches!(self, Program::Blastp | Program::Blastx)
    }

    /// Whether LOSAT's program writes the format: the one table of the formats of every
    /// program. One arm per format, grouped by the line OF track that adds it, so that a
    /// change for one format edits its own line. No NCBI counterpart: NCBI's programs write
    /// every format of `ReportFormat` (blastn alone writes SAM, `report_format`).
    pub fn writes(self, format: ReportFormat) -> bool {
        use ReportFormat::*;
        match format {
            // Written today by every program.
            Pairwise => true,
            Tabular => true,
            TabularWithComments => true,
            // Line OF track T (text formats): none yet.
            QueryAnchoredIdentities => false,
            QueryAnchoredNoIdentities => false,
            FlatQueryAnchoredIdentities => false,
            FlatQueryAnchoredNoIdentities => false,
            CommaSeparatedValues => false,
            TaxFormat => false,
            CommaSeparatedValuesWithHeader => false,
            // Line OF track S (structured formats): 5 for BLASTP and BLASTN (OF4a), BLASTX
            // (OF4b part 1), TBLASTN (OF4b part 3) and TBLASTX (OF4b part 4): every program;
            // 15 and 16 (JSON_S, XML2_S) for BLASTP (OF5a part 1), BLASTN (OF5a part 2),
            // TBLASTN (OF5b part 2) and TBLASTX (OF5b part 3).
            Xml => true,
            AsnText => true,
            AsnBinary => true,
            ArchiveFormat => false,
            JsonSeqalign => true,
            Json => false,
            Xml2 => false,
            JsonS => matches!(
                self,
                Program::Blastp | Program::Blastn | Program::Tblastn | Program::Tblastx
            ),
            Xml2S => matches!(
                self,
                Program::Blastp | Program::Blastn | Program::Tblastn | Program::Tblastx
            ),
            Sam => false,
        }
    }
}

/// NCBI's reading of a `-outfmt` value, which `SetOptions` does before every other handler
/// (`ArchiveFormatRequested`), with NCBI's errors.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2801-2851
/// ```c
///     if (args[kArgOutputFormat]) {
///         string fmt_choice =
///             NStr::TruncateSpaces(args[kArgOutputFormat].AsString());
///         string::size_type pos;
///         if ( (pos = fmt_choice.find_first_of(' ')) != string::npos) {
///             custom_fmt_spec.assign(fmt_choice, pos+1,
///                                    fmt_choice.size()-(pos+1));
///             fmt_choice.erase(pos);
///         }
///         if(!custom_fmt_spec.empty()) {
///             if(NStr::StartsWith(custom_fmt_spec, "delim")) {
///                 vector <string> tokens;
///                 NStr::Split(custom_fmt_spec," ",tokens);
///                 if(tokens.size() > 0) {
///                     string tag;
///                     bool isValid = NStr::SplitInTwo(tokens[0],"=",tag,custom_delim);
///                     if(!isValid) {
///                         string msg("Delimiter format is invalid. Valid format is delim=<delimiter value>");
///                         NCBI_THROW(CInputException, eInvalidInput, msg);
///                     }
///                     else {
///                         custom_fmt_spec = NStr::Replace(custom_fmt_spec,tokens[0],"");
///                         custom_fmt_spec = NStr::TruncateSpaces(custom_fmt_spec);
///                     }
///                 }
///             }
///         }
///         int val = 0;
///         try { val = NStr::StringToInt(fmt_choice); }
///         catch (const CStringException&) {   // probably a conversion error
///             CNcbiOstrstream os;
///             os << "'" << fmt_choice << "' is not a valid output format";
///             string msg = CNcbiOstrstreamToString(os);
///             NCBI_THROW(CInputException, eInvalidInput, msg);
///         }
///         if (val < 0 || val >= static_cast<int>(eEndValue)) {
///             string msg("Formatting choice is out of range");
///             throw std::out_of_range(msg);
///         }
///         fmt_type = static_cast<EOutputFormat>(val);
///         if ( !(fmt_type == eTabular ||
///                fmt_type == eTabularWithComments ||
///                fmt_type == eCommaSeparatedValues ||
///                fmt_type == eCommaSeparatedValuesWithHeader ||
///                fmt_type == eSAM) ) {
///                custom_fmt_spec.clear();
///         }
///     }
/// ```
/// NCBI reference: c++/src/app/blast/blast_app_util.hpp:260-263
/// ```c
///     catch (const std::exception& e) {                                       \
///         LOG_POST(Error << "Error: " << e.what());                           \
///         exit_code = BLAST_UNKNOWN_ERROR;                                    \
///     }                                                                       \
/// ```
/// `NStr::TruncateSpaces` removes the characters of C's `isspace`, and `NStr::StringToInt`
/// reads what `i32::from_str` reads.
pub fn parse_formatting_string(spec: &str) -> anyhow::Result<FormatChoice> {
    let is_space = |c: char| matches!(c, ' ' | '\t' | '\n' | '\x0b' | '\x0c' | '\r');
    let choice = spec.trim_matches(is_space);
    let (choice, custom) = choice.split_once(' ').unwrap_or((choice, ""));
    let mut custom = custom.to_string();
    let mut delimiter = String::new();
    if custom.starts_with("delim") {
        let token = custom.split(' ').next().unwrap_or_default().to_string();
        let Some((_, value)) = token.split_once('=') else {
            return Err(options_error(
                "Delimiter format is invalid. Valid format is delim=<delimiter value>",
            ));
        };
        delimiter = value.to_string();
        custom = custom
            .replace(&token, "")
            .trim_matches(is_space)
            .to_string();
    }
    let Ok(number) = choice.parse::<i32>() else {
        return Err(options_error(&format!(
            "'{choice}' is not a valid output format"
        )));
    };
    if !(0..END_VALUE).contains(&number) {
        return Err(NativeError {
            exit: 255,
            message: "Error: Formatting choice is out of range\n".to_string(),
        }
        .into());
    }
    let keeps_spec = matches!(
        number,
        TABULAR
            | TABULAR_WITH_COMMENTS
            | COMMA_SEPARATED_VALUES
            | COMMA_SEPARATED_VALUES_WITH_HEADER
            | SAM
    );
    if !keeps_spec {
        custom.clear();
        delimiter.clear();
    }
    Ok(FormatChoice {
        number,
        spec: custom,
        delimiter,
    })
}

/// The report of a format choice. A format that NCBI runs and LOSAT's program does not
/// write (`Program::writes`) is rejected here, before the inputs are read; `None` is a
/// format that NCBI itself rejects later (`formatting_handler_check`, `xinclude_check`),
/// where the run then stops: SAM (17) in a program without `eIsSAM`, AIRR rearrangement (19),
/// FASTA (21), and the XInclude formats (13, 14) written to standard output. A custom
/// specification is rejected unless the program writes custom fields
/// (`Program::custom_fields`), and a custom delimiter always.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2874-2886
/// ```c
///     ParseFormattingString(args, m_OutputFormat, m_CustomOutputFormatSpec,m_CustomDelim);
///     if((m_OutputFormat == eSAM) && !(m_FormatFlags & eIsSAM) ){
///     		NCBI_THROW(CInputException, eInvalidInput,
///     		                        "SAM format is only applicable to blastn" );
///     }
/// ```
/// NCBI reference: c++/src/app/blast/blast_app_util.cpp:891-893
/// ```c
/// 	if((fmt ==  CFormattingArgs::eXml2) || (fmt ==  CFormattingArgs::eJson)) {
/// 	   if (s == "-"){
/// 		   string f_str = (fmt == CFormattingArgs::eXml2) ? "14.": "13.";
/// ```
pub fn report_format(
    choice: &FormatChoice,
    program: Program,
    out: Option<&Path>,
) -> anyhow::Result<Option<ReportFormat>> {
    let to_standard_output = out.is_none_or(|path| path.as_os_str() == "-");
    let name = program.name();
    let format = match ReportFormat::from_number(choice.number) {
        Some(ReportFormat::Sam) if !program.is_sam() => return Ok(None),
        Some(ReportFormat::Json | ReportFormat::Xml2) if to_standard_output => return Ok(None),
        Some(format) if program.writes(format) => format,
        None if matches!(choice.number, AIRR_REARRANGEMENT | FASTA) => return Ok(None),
        _ => anyhow::bail!(
            "output format {} is not supported by LOSAT's {name}",
            choice.number
        ),
    };
    if choice.custom() && !program.custom_fields() {
        anyhow::bail!(
            "a custom output format specification (fields or a delimiter) is not supported by LOSAT's {name}"
        );
    }
    if !choice.delimiter.is_empty() {
        anyhow::bail!(
            "a custom delimiter (delim=) in the output format is not supported by LOSAT's {name}"
        );
    }
    Ok(Some(format))
}

/// NCBI's errors for the formats of other programs, which the formatting handler raises
/// after the input and output files are opened.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2874-2886
/// ```c
///     ParseFormattingString(args, m_OutputFormat, m_CustomOutputFormatSpec,m_CustomDelim);
///     if((m_OutputFormat == eSAM) && !(m_FormatFlags & eIsSAM) ){
///     		NCBI_THROW(CInputException, eInvalidInput,
///     		                        "SAM format is only applicable to blastn" );
///     }
///     if((m_OutputFormat == eAirrRearrangement) && !(m_FormatFlags & eIsAirrRearrangement) ){
///         NCBI_THROW(CInputException, eInvalidInput,
///                    "AIRR rearrangement format is only applicable to igblastn" );
///     }
///     if (m_OutputFormat == eFasta) {
///         NCBI_THROW(CInputException, eInvalidInput,
///                    "FASTA output format is only applicable to magicblast");
///     }
/// ```
pub fn formatting_handler_check(choice: &FormatChoice, program: Program) -> anyhow::Result<()> {
    match choice.number {
        SAM if !program.is_sam() => Err(options_error("SAM format is only applicable to blastn")),
        AIRR_REARRANGEMENT => Err(options_error(
            "AIRR rearrangement format is only applicable to igblastn",
        )),
        FASTA => Err(options_error(
            "FASTA output format is only applicable to magicblast",
        )),
        _ => Ok(()),
    }
}

/// NCBI's error for the XInclude formats (13 and 14) written to standard output, after
/// `Query is Empty!`.
///
/// NCBI reference: c++/src/app/blast/blast_app_util.cpp:888-899
/// ```c
/// UseXInclude(const CFormattingArgs & f, const string & s)
/// {
/// 	CFormattingArgs::EOutputFormat fmt = f.GetFormattedOutputChoice();
/// 	if((fmt ==  CFormattingArgs::eXml2) || (fmt ==  CFormattingArgs::eJson)) {
/// 	   if (s == "-"){
/// 		   string f_str = (fmt == CFormattingArgs::eXml2) ? "14.": "13.";
/// 		   NCBI_THROW(CInputException, eEmptyUserInput,
/// 		              "Please provide a file name for outfmt " + f_str);
/// 	   }
/// 	   return true;
/// 	}
/// 	return false;
/// }
/// ```
pub fn xinclude_check(choice: &FormatChoice) -> anyhow::Result<()> {
    match choice.number {
        JSON => Err(options_error("Please provide a file name for outfmt 13.")),
        XML2 => Err(options_error("Please provide a file name for outfmt 14.")),
        _ => Ok(()),
    }
}

/// A `-seg` value as NCBI's filtering handler reads it, with NCBI's errors. A locut or
/// hicut that `NStr::StringToDouble` reads but that is not a finite decimal number (an
/// infinity, NaN or a hexadecimal number) is rejected with LOSAT's message.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:396-406
/// ```c
///         if (m_QueryIsProtein && args[kArgSegFiltering]) {
///             const string& seg_opts = args[kArgSegFiltering].AsString();
///             if (seg_opts == kDfltArgNoFiltering) {
///                 opt.SetSegFiltering(false);
///             } else if (seg_opts == kDfltArgApplyFiltering) {
///                 opt.SetSegFiltering(true);
///             } else {
///                 x_TokenizeFilteringArgs(seg_opts, tokens);
///                 opt.SetSegFilteringWindow(NStr::StringToInt(tokens[0]));
///                 opt.SetSegFilteringLocut(NStr::StringToDouble(tokens[1]));
///                 opt.SetSegFilteringHicut(NStr::StringToDouble(tokens[2]));
///             }
///         }
/// ```
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:375-384
/// ```c
/// CFilteringArgs::x_TokenizeFilteringArgs(const string& filtering_args,
///                                         vector<string>& output) const
/// {
///     output.clear();
///     NStr::Split(filtering_args, " ", output);
///     if (output.size() != 3) {
///         NCBI_THROW(CInputException, eInvalidInput,
///                    "Invalid number of arguments to filtering option");
///     }
/// }
/// ```
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:423-428
/// ```c
///     } catch (const CStringException& e) {
///         if (e.GetErrCode() == CStringException::eConvert) {
///             NCBI_THROW(CInputException, eInvalidInput,
///                        "Invalid input for filtering parameters");
///         }
///     }
/// ```
pub fn parse_seg_option(value: &str, program: &str) -> anyhow::Result<SegSpec> {
    match value {
        "no" => return Ok(SegSpec::No),
        "yes" => return Ok(SegSpec::Yes),
        _ => {}
    }
    let tokens: Vec<&str> = value.split(' ').collect();
    if tokens.len() != 3 {
        return Err(options_error(
            "Invalid number of arguments to filtering option",
        ));
    }
    let invalid = || options_error("Invalid input for filtering parameters");
    let window = tokens[0].parse::<i32>().map_err(|_| invalid())?;
    let mut cuts = [0.0; 2];
    for (cut, token) in cuts.iter_mut().zip(&tokens[1..]) {
        use crate::blastinput::value_parsers::{
            ncbi_string_to_double, strtod_reads_whole, NcbiDoubleError,
        };
        // NCBI reads the values with `strtod` and fails on one that it does not read to the
        // end (`strtod_reads_whole`); LOSAT reads the finite decimal ones.
        *cut = match ncbi_string_to_double(token) {
            Ok(number) if number.is_finite() => number,
            Ok(_) | Err(NcbiDoubleError::Unsupported) if strtod_reads_whole(token) => {
                anyhow::bail!(
                    "the SEG locut or hicut {token:?} (not a finite decimal number) is not supported by LOSAT's {program}"
                )
            }
            Ok(_) | Err(_) => return Err(invalid()),
        };
    }
    Ok(SegSpec::WindowLocutHicut {
        window,
        locut: cuts[0],
        hicut: cuts[1],
    })
}

/// NCBI's composition-based statistics modes (`ECompoAdjustModes`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompositionMode {
    NoCompositionBasedStats = 0,
    CompositionBasedStats = 1,
    CompositionMatrixAdjust = 2,
    CompoForceFullMatrixAdjust = 3,
}

/// A `-comp_based_stats` value as NCBI reads it for BLASTP or TBLASTN (not RPS-BLAST or
/// DELTA-BLAST): the first character chooses the mode, any other first character (or an
/// empty value, whose first character is the terminating NUL) is mode 0, and for BLASTP a
/// `u` or `U` second character asks for unified P-values with an adjusting mode.
/// Returns the mode and whether unified P-values are asked for.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:834-892
/// ```c
///         ECompoAdjustModes compo_mode = eNoCompositionBasedStats;
///
///         switch (comp_stat_string[0]) {
///             case '0': case 'F': case 'f':
///                 compo_mode = eNoCompositionBasedStats;
///                 break;
///             case '1':
///                 compo_mode = eCompositionBasedStats;
///                 break;
///             case 'D': case 'd':
/// ...
///                 else {
///                     compo_mode = eCompositionMatrixAdjust;
///                 }
///                 break;
///             case '2':
///                 compo_mode = eCompositionMatrixAdjust;
///                 break;
///             case '3':
///                 compo_mode = eCompoForceFullMatrixAdjust;
///                 break;
///             case 'T': case 't':
///                 compo_mode = (program == eRPSBlast || program == eRPSTblastn || program == eDeltaBlast) ?
///                     eCompositionBasedStats : eCompositionMatrixAdjust;
///                 break;
///         }
/// ...
///         if (ungapped && *ungapped && compo_mode != eNoCompositionBasedStats) {
///             NCBI_THROW(CInputException, eInvalidInput,
///                        "Composition-adjusted searched are not supported with "
///                        "an ungapped search, please add -comp_based_stats F or "
///                        "do a gapped search");
///         }
///
///         opt.SetCompositionBasedStats(compo_mode);
///         if (program == eBlastp &&
///             compo_mode != eNoCompositionBasedStats &&
///             tolower(comp_stat_string[1]) == 'u') {
///             opt.SetUnifiedP(1);
///         }
/// ```
pub fn parse_comp_based_stats(
    value: &str,
    blastp: bool,
    ungapped: bool,
) -> anyhow::Result<(CompositionMode, bool)> {
    let bytes = value.as_bytes();
    let mode = match bytes.first() {
        Some(b'1') => CompositionMode::CompositionBasedStats,
        Some(b'D' | b'd' | b'2' | b'T' | b't') => CompositionMode::CompositionMatrixAdjust,
        Some(b'3') => CompositionMode::CompoForceFullMatrixAdjust,
        _ => CompositionMode::NoCompositionBasedStats,
    };
    if ungapped && mode != CompositionMode::NoCompositionBasedStats {
        return Err(options_error(
            "Composition-adjusted searched are not supported with an ungapped search, please add -comp_based_stats F or do a gapped search",
        ));
    }
    let unified_p = blastp
        && mode != CompositionMode::NoCompositionBasedStats
        && bytes.get(1).is_some_and(|c| c.eq_ignore_ascii_case(&b'u'));
    Ok((mode, unified_p))
}

/// The query splitter that NCBI sets up for every query batch (`CBlastPrelimSearch`, in
/// `BlastSetupPreliminarySearchEx`), which reads `CHUNK_SIZE` and then `OVERLAP_CHUNK_SIZE`
/// but does not split an ungapped search. The outer `Err` is LOSAT's rejection of a value
/// that NCBI cannot convert (its `CStringException` names its build's source files, exit
/// 255); the inner one is NCBI's error for a chunk size that is not divisible by 3.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:55-61
/// ```c
///     m_ChunkSize = SplitQuery_GetChunkSize(m_Options->GetProgram());
///     m_LocalQueryData = m_QueryFactory->MakeLocalQueryData(m_Options);
///     m_TotalQueryLength = m_LocalQueryData->GetSumOfSequenceLengths();
///     m_NumChunks = SplitQuery_CalculateNumChunks(m_Options->GetProgramType(),
///         &m_ChunkSize, m_TotalQueryLength, m_LocalQueryData->GetNumQueries());
///     /* No split for ungapped mode JIRA SB-1082 */
///     if (!options->GetGappedMode()) m_NumChunks = 1;
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp:98-103
/// ```c
///     const EBlastProgramType prog_type(EProgramToEBlastProgramType(program));
///     if (Blast_QueryIsTranslated(prog_type) && !Blast_SubjectIsPssm(prog_type) &&
///         (retval % CODON_LENGTH) != 0) {
///         NCBI_THROW(CBlastException, eInvalidArgument,
///                    "Split query chunk size must be divisible by 3");
///     }
/// ```
/// `SplitQuery_ShouldSplit` is true for blastp, tblastn and tblastx
/// (split_query_aux_priv.cpp:73-97), so `SplitQuery_CalculateNumChunks` always reads the
/// overlap (`blastn/query_split.rs` `SplitSizes`; the `int` becomes a 64-bit `size_t`). The
/// chunk size must be divisible by 3 only for a translated query (`translated_query`,
/// tblastx). Whether a blastp or tblastn query batch is split is
/// `common/protein_query_split.rs`.
pub fn check_query_split_environment(
    program: &str,
    translated_query: bool,
) -> anyhow::Result<anyhow::Result<()>> {
    let read = |variable: &str| read_query_split_variable(program, variable);
    if let Some(chunk_size) = read("CHUNK_SIZE")? {
        if translated_query && (i64::from(chunk_size) as u64) % 3 != 0 {
            return Ok(Err(NativeError {
                exit: 3,
                message: "BLAST engine error: Split query chunk size must be divisible by 3\n"
                    .to_string(),
            }
            .into()));
        }
    }
    read("OVERLAP_CHUNK_SIZE")?;
    Ok(Ok(()))
}

/// `CHUNK_SIZE` or `OVERLAP_CHUNK_SIZE` as NCBI reads it (`None` when unset or blank).
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp:58-61
/// ```c
///     char* chunk_sz_str = getenv("CHUNK_SIZE");
///     if (chunk_sz_str && !NStr::IsBlank(chunk_sz_str)) {
///         retval = NStr::StringToInt(chunk_sz_str);
/// ```
fn read_query_split_variable(program: &str, variable: &str) -> anyhow::Result<Option<i32>> {
    match std::env::var_os(variable) {
        Some(value) if !crate::blastinput::input_files::is_blank(value.as_encoded_bytes()) => {
            match crate::blastinput::query_batch::ncbi_string_to_int(&value) {
                Some(number) => Ok(Some(number)),
                None => anyhow::bail!(
                    "the environment variable {variable} has the value {:?}, which NCBI BLAST+ cannot convert to an int (it stops with a CStringException that names its build's source files); this is not supported by LOSAT's {program}",
                    value.to_string_lossy()
                ),
            }
        }
        _ => Ok(None),
    }
}

/// NCBI's query chunk size and overlap of a blastp or tblastn search
/// (`common/protein_query_split.rs` `protein_split_sizes`), with the environment variables
/// `CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE` read as NCBI reads them.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/split_query_cxx.cpp:55-58
/// ```c
///     m_ChunkSize = SplitQuery_GetChunkSize(m_Options->GetProgram());
///     m_LocalQueryData = m_QueryFactory->MakeLocalQueryData(m_Options);
///     m_TotalQueryLength = m_LocalQueryData->GetSumOfSequenceLengths();
///     m_NumChunks = SplitQuery_CalculateNumChunks(m_Options->GetProgramType(),
/// ```
pub fn protein_query_split_sizes(
    program: &str,
    default_chunk_size: u32,
) -> anyhow::Result<crate::algorithm::blastn::query_split::SplitSizes> {
    Ok(
        crate::algorithm::common::protein_query_split::protein_split_sizes(
            default_chunk_size,
            read_query_split_variable(program, "CHUNK_SIZE")?,
            read_query_split_variable(program, "OVERLAP_CHUNK_SIZE")?,
        ),
    )
}

/// The environment that changes NCBI's blastp, tblastn or tblastx in a way that LOSAT does
/// not reproduce.
///
/// NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.cpp:206-210
/// ```c
/// 	char* bl2seq_legacy = getenv("BL2SEQ_LEGACY");
/// 	if (bl2seq_legacy)
///         	db_adapter.Reset(new CLocalDbAdapter(subjects, opts_hndl, false));
/// 	else
///         	db_adapter.Reset(new CLocalDbAdapter(subjects, opts_hndl, true));
/// ```
///
/// NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.cpp:732-737
/// ```c
/// 		char * pre_fetch_limit_str = getenv("PRE_FETCH_SEQS_LIMIT");
/// 		if (pre_fetch_limit_str) {
/// 			int pre_fetch_limit = NStr::StringToInt(pre_fetch_limit_str);
/// 			if(pre_fetch_limit == 0) {
/// 				return false;
/// 			}
/// ```
/// An integer only decides whether the report's sequences are fetched ahead; a value
/// that `NStr::StringToInt` cannot convert (an empty one too) stops NCBI before the report
/// of each query batch with a `CStringException` that names its build's source files (as
/// BLASTN's `check_unsupported_environment`).
pub fn check_unsupported_environment(program: &str) -> anyhow::Result<()> {
    if std::env::var_os("BL2SEQ_LEGACY").is_some() {
        anyhow::bail!(
            "the environment variable BL2SEQ_LEGACY, which makes NCBI BLAST+ search each subject on its own and write its legacy bl2seq report, is not supported by LOSAT's {program}"
        );
    }
    if let Some(value) = std::env::var_os("PRE_FETCH_SEQS_LIMIT") {
        if crate::blastinput::query_batch::ncbi_string_to_int(&value).is_none() {
            anyhow::bail!(
                "the environment variable PRE_FETCH_SEQS_LIMIT has the value {:?}, which NCBI BLAST+ cannot convert to an int (it stops with a CStringException that names its build's source files); this is not supported by LOSAT's {program}",
                value.to_string_lossy()
            );
        }
    }
    Ok(())
}

/// The query batch size of blastp, tblastn or tblastx (`GetQueryBatchSize`, with the
/// program's `default`: tblastx 10002 nucleotides, tblastn 20000 and blastp 10000
/// residues), as the `TSeqPos` that `CBlastInput` compares with.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_input_aux.cpp:117-119
/// ```c
///     case eTblastn:
///         retval = 20000;
///         break;
/// ```
pub fn query_batch_size(program: &str, default: u32) -> anyhow::Result<u32> {
    match crate::blastinput::query_batch::get_query_batch_size(
        std::env::var_os("BATCH_SIZE").as_deref(),
    ) {
        Ok(0) if std::env::var_os("BATCH_SIZE").is_none() => Ok(default),
        // `CBlastInput` stores the `int` as a `TSeqPos` (blast_input.hpp:313,364).
        Ok(size) => Ok(size as u32),
        Err(value) => Err(anyhow::anyhow!(
            "the BATCH_SIZE value '{value}' is not an integer; NCBI stops with a \
             CStringException for it, which is not supported by LOSAT's {program}"
        )),
    }
}

/// The environment variable `OLD_FSC`, which makes NCBI's protein programs compute without
/// the Gumbel block (finite-size correction) that LOSAT's blastp and tblastn use.
///
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:922-923
/// ```c
///     use_old_fsc = getenv("OLD_FSC");
///     if (!use_old_fsc) sbp->gbp = s_BlastGumbelBlkNew();
/// ```
pub fn check_old_fsc(program: &str) -> anyhow::Result<()> {
    if std::env::var_os("OLD_FSC").is_some() {
        anyhow::bail!(
            "the environment variable OLD_FSC, which makes NCBI BLAST+ compute without its finite-size correction (the Gumbel block), is not supported by LOSAT's {program}"
        );
    }
    Ok(())
}

/// The window of the two-hit word finder where NCBI's diagonal table never gets a length.
///
/// NCBI doubles the length of the table from 1 while it is less than the query length plus
/// the window, in `Int4`. When the sum is over 2^30 the length reaches 2^31, wraps to
/// `INT_MIN` and then 0, and the loop never ends; when the sum itself passes 2^31 − 1 it
/// wraps negative, the table gets one cell, and the later `Int4` offsets wrap (NCBI reports
/// no hits). Neither has a use, and LOSAT rejects both (decisions D8 and D11 of
/// docs/evidence/losat_web_e2e/AUTHORITY.md). `query_length` is the length of the
/// concatenated query of the batch (`BLAST_SequenceBlk->length`).
///
/// NCBI reference: c++/src/algo/blast/core/blast_extend.c:52-57
/// ```c
///                 diag_array_length = 1;
///                 /* What power of 2 is just longer than the query? */
///                 while (diag_array_length < (qlen+window_size))
///                 {
///                         diag_array_length = diag_array_length << 1;
///                 }
/// ```
pub fn check_diag_table_window(
    query_length: i32,
    window: i32,
    program: &str,
) -> anyhow::Result<()> {
    if i64::from(query_length) + i64::from(window) > 1 << 30 {
        anyhow::bail!(
            "a -window_size of {window} with a query of {query_length} letters (their sum is over 2^30, where NCBI BLAST+ never ends or wraps a 32-bit integer) is not supported by LOSAT's {program}"
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diag_table_windows_over_2_to_the_30_are_rejected() {
        // NCBI blast_extend.c:52-57: the doubling loop ends only when qlen + window <= 2^30.
        assert!(check_diag_table_window(1000, (1 << 30) - 1000, "BLASTP").is_ok());
        for window in [(1 << 30) - 999, i32::MAX] {
            let error = format!(
                "{:#}",
                check_diag_table_window(1000, window, "TBLASTX").unwrap_err()
            );
            assert!(
                error.contains("not supported by LOSAT's TBLASTX"),
                "{window}: {error}"
            );
        }
    }

    #[test]
    fn seg_values_that_strtod_does_not_read_are_ncbi_errors() {
        // NCBI blast_args.cpp:405-406,423-428 (StringToDouble with default flags).
        for value in [
            "12 2.5x 2.5",
            "12 2,2 2.5",
            "12 1e 2.5",
            "12 0x 2.5",
            "12 . 2.5",
        ] {
            let error = format!("{:#}", parse_seg_option(value, "TBLASTX").unwrap_err());
            assert!(
                error.contains("Invalid input for filtering parameters"),
                "{value}: {error}"
            );
        }
        for value in ["12 0x10 2.5", "12 +inf 2.5", "12 1e400 2.5", "12 2.2 -nan"] {
            let error = format!("{:#}", parse_seg_option(value, "TBLASTX").unwrap_err());
            assert!(
                error.contains("not supported by LOSAT's TBLASTX"),
                "{value}: {error}"
            );
        }
        assert!(parse_seg_option("12 2.2 2.5", "TBLASTX").is_ok());
    }

    fn message(error: anyhow::Error) -> String {
        match error.downcast::<NativeError>() {
            Ok(native) => format!("{} {}", native.exit, native.message),
            Err(other) => format!("losat {other}"),
        }
    }

    // NCBI BLAST+ 2.17.0 -outfmt spellings (blast_args.cpp:2801-2851).
    #[test]
    fn formatting_strings_are_read_as_ncbi_reads_them() {
        assert_eq!(
            parse_formatting_string(" 6 ").unwrap(),
            FormatChoice {
                number: 6,
                spec: String::new(),
                delimiter: String::new(),
            }
        );
        assert_eq!(parse_formatting_string("+6").unwrap().number, 6);
        assert_eq!(parse_formatting_string("06").unwrap().normalized(), "6");
        assert!(parse_formatting_string("6 std").unwrap().custom());
        assert!(!parse_formatting_string("0 std").unwrap().custom());
        assert_eq!(parse_formatting_string("0 std").unwrap().normalized(), "0");
        let delimited = parse_formatting_string("7 delim=, qseqid  sseqid ").unwrap();
        assert_eq!(
            (delimited.spec.as_str(), delimited.delimiter.as_str()),
            ("qseqid  sseqid", ",")
        );
        assert_eq!(delimited.normalized(), "7 qseqid  sseqid");
        assert!(!parse_formatting_string("6 delim=").unwrap().custom());
        assert!(message(parse_formatting_string("6 delim").unwrap_err())
            .contains("Delimiter format is invalid"));
        assert_eq!(
            message(parse_formatting_string("abc").unwrap_err()),
            "1 BLAST query/options error: 'abc' is not a valid output format\nPlease refer to the BLAST+ user manual.\n"
        );
        assert_eq!(
            message(parse_formatting_string("22").unwrap_err()),
            "255 Error: Formatting choice is out of range\n"
        );
        assert_eq!(
            message(parse_formatting_string("-1").unwrap_err()),
            "255 Error: Formatting choice is out of range\n"
        );
    }

    const PROGRAMS: [Program; 5] = [
        Program::Blastn,
        Program::Blastp,
        Program::Blastx,
        Program::Tblastn,
        Program::Tblastx,
    ];

    // NCBI blast_args.hpp:1023-1074: the numbers of `EOutputFormat`; 19 (AIRR) and 21 (FASTA)
    // are igblastn's and magicblast's, 22 is `eEndValue`.
    #[test]
    fn report_formats_are_ncbis_format_numbers() {
        for number in (0..=18).chain([20]) {
            let format = ReportFormat::from_number(number).unwrap();
            assert_eq!(format.number(), number);
        }
        for number in [-1, AIRR_REARRANGEMENT, FASTA, END_VALUE, 99] {
            assert_eq!(ReportFormat::from_number(number), None, "{number}");
        }
        assert_eq!(ReportFormat::from_number(17), Some(ReportFormat::Sam));
        assert_eq!(
            ReportFormat::from_number(20),
            Some(ReportFormat::CommaSeparatedValuesWithHeader)
        );
    }

    // The one table of the formats LOSAT writes (`Program::writes`): 0, 6, 7, 8, 9 and 12 in
    // every program (8 and 12: line OF track S, OF3b; 9: OF3c), and 5 in every program
    // (BLASTN, BLASTP: OF4a; BLASTX, TBLASTN, TBLASTX: OF4b); 15 and 16 in BLASTP and BLASTN
    // (OF5a), TBLASTN (OF5b part 2) and TBLASTX (OF5b part 3).
    #[test]
    fn every_program_writes_its_formats_only() {
        for program in PROGRAMS {
            let written: Vec<i32> = (0..END_VALUE)
                .filter_map(ReportFormat::from_number)
                .filter(|&format| program.writes(format))
                .map(ReportFormat::number)
                .collect();
            let expected: &[i32] = if matches!(
                program,
                Program::Blastp | Program::Blastn | Program::Tblastn | Program::Tblastx
            ) {
                &[0, 5, 6, 7, 8, 9, 12, 15, 16]
            } else {
                &[0, 5, 6, 7, 8, 9, 12]
            };
            assert_eq!(written, expected, "{program:?}");
        }
        assert_eq!(
            PROGRAMS.map(Program::name),
            ["BLASTN", "BLASTP", "BLASTX", "TBLASTN", "TBLASTX"]
        );
        assert_eq!(
            PROGRAMS.map(Program::is_sam),
            [true, false, false, false, false]
        );
        assert_eq!(
            PROGRAMS.map(Program::custom_fields),
            [false, true, true, false, false]
        );
    }

    // NCBI blast_format.cpp:1431-1452: the structured reports return before the result
    // warnings.
    #[test]
    fn structured_formats_write_no_query_warnings() {
        let silent: Vec<i32> = (0..END_VALUE)
            .filter_map(ReportFormat::from_number)
            .filter(|format| !format.writes_query_warnings())
            .map(ReportFormat::number)
            .collect();
        assert_eq!(silent, [5, 8, 9, 12, 13, 14, 15, 16, 17]);
    }

    #[test]
    fn report_formats_defer_ncbis_own_errors() {
        let choice = |number| FormatChoice {
            number,
            spec: String::new(),
            delimiter: String::new(),
        };
        let file = Some(Path::new("o.json"));
        for program in PROGRAMS {
            let name = program.name();
            assert_eq!(
                report_format(&choice(0), program, None).unwrap(),
                Some(ReportFormat::Pairwise)
            );
            // NCBI rejects these later (`formatting_handler_check`, `xinclude_check`).
            for number in [AIRR_REARRANGEMENT, FASTA, JSON, XML2] {
                assert_eq!(
                    report_format(&choice(number), program, None).unwrap(),
                    None,
                    "{name} {number}"
                );
            }
            assert_eq!(
                report_format(&choice(JSON), program, Some(Path::new("-"))).unwrap(),
                None
            );
            // Formats that NCBI writes and LOSAT does not, before the inputs are read.
            for number in (1..=16)
                .filter(|&number| !matches!(number, PAIRWISE | TABULAR | TABULAR_WITH_COMMENTS))
                .filter(|&number| {
                    let format = ReportFormat::from_number(number).unwrap();
                    !program.writes(format)
                })
                .chain([18, 20])
            {
                assert_eq!(
                    report_format(&choice(number), program, file)
                        .unwrap_err()
                        .to_string(),
                    format!("output format {number} is not supported by LOSAT's {name}")
                );
            }
            let sam = report_format(&choice(SAM), program, None);
            if program.is_sam() {
                assert_eq!(
                    sam.unwrap_err().to_string(),
                    "output format 17 is not supported by LOSAT's BLASTN"
                );
            } else {
                assert_eq!(sam.unwrap(), None);
            }
            let delimited = parse_formatting_string("6 delim=, std").unwrap();
            assert!(report_format(&delimited, program, None)
                .unwrap_err()
                .to_string()
                .ends_with(&format!("is not supported by LOSAT's {name}")));
        }
        // A number outside NCBI's enum (no `-outfmt` reads so) is LOSAT's rejection.
        assert!(report_format(&choice(END_VALUE), Program::Blastp, None).is_err());
    }

    // NCBI blast_args.cpp:2874-2886 and blast_app_util.cpp:888-899.
    #[test]
    fn ncbis_later_format_errors_are_its_messages() {
        let choice = |number| FormatChoice {
            number,
            spec: String::new(),
            delimiter: String::new(),
        };
        let options = |text: &str| {
            format!(
                "1 BLAST query/options error: {text}\nPlease refer to the BLAST+ user manual.\n"
            )
        };
        for program in PROGRAMS {
            let sam = formatting_handler_check(&choice(SAM), program);
            if program.is_sam() {
                assert!(sam.is_ok());
            } else {
                assert_eq!(
                    message(sam.unwrap_err()),
                    options("SAM format is only applicable to blastn")
                );
            }
            assert_eq!(
                message(
                    formatting_handler_check(&choice(AIRR_REARRANGEMENT), program).unwrap_err()
                ),
                options("AIRR rearrangement format is only applicable to igblastn")
            );
            assert_eq!(
                message(formatting_handler_check(&choice(FASTA), program).unwrap_err()),
                options("FASTA output format is only applicable to magicblast")
            );
            for number in (0..=16).chain([18, 20]) {
                assert!(formatting_handler_check(&choice(number), program).is_ok());
            }
        }
        assert_eq!(
            message(xinclude_check(&choice(JSON)).unwrap_err()),
            options("Please provide a file name for outfmt 13.")
        );
        assert_eq!(
            message(xinclude_check(&choice(XML2)).unwrap_err()),
            options("Please provide a file name for outfmt 14.")
        );
        for number in (0..=12).chain(15..=21) {
            assert!(xinclude_check(&choice(number)).is_ok());
        }
    }

    #[test]
    fn seg_values_follow_ncbis_filtering_handler() {
        assert_eq!(parse_seg_option("no", "BLASTP").unwrap(), SegSpec::No);
        assert!(message(parse_seg_option("12 2.2", "BLASTP").unwrap_err())
            .contains("Invalid number of arguments to filtering option"));
        assert!(
            message(parse_seg_option("12  2.2 2.5", "BLASTP").unwrap_err())
                .contains("Invalid number of arguments to filtering option")
        );
        assert!(
            message(parse_seg_option("x 2.2 2.5", "BLASTP").unwrap_err())
                .contains("Invalid input for filtering parameters")
        );
        assert!(
            message(parse_seg_option("12 inf 2.5", "BLASTP").unwrap_err())
                .contains("Invalid input for filtering parameters")
        );
        assert!(
            message(parse_seg_option("12 +inf 2.5", "BLASTP").unwrap_err()).starts_with("losat ")
        );
        assert_eq!(
            parse_seg_option("+12 2.2 .5", "BLASTP").unwrap(),
            SegSpec::WindowLocutHicut {
                window: 12,
                locut: 2.2,
                hicut: 0.5
            }
        );
    }

    #[test]
    fn composition_modes_read_the_first_character() {
        for (value, mode) in [
            ("0", CompositionMode::NoCompositionBasedStats),
            ("F", CompositionMode::NoCompositionBasedStats),
            ("", CompositionMode::NoCompositionBasedStats),
            ("x", CompositionMode::NoCompositionBasedStats),
            ("4", CompositionMode::NoCompositionBasedStats),
            ("1", CompositionMode::CompositionBasedStats),
            ("D", CompositionMode::CompositionMatrixAdjust),
            ("t", CompositionMode::CompositionMatrixAdjust),
            ("2xyz", CompositionMode::CompositionMatrixAdjust),
            ("3", CompositionMode::CompoForceFullMatrixAdjust),
        ] {
            assert_eq!(
                parse_comp_based_stats(value, true, false).unwrap().0,
                mode,
                "{value}"
            );
        }
        assert!(parse_comp_based_stats("2U", true, false).unwrap().1);
        assert!(!parse_comp_based_stats("2U", false, false).unwrap().1);
        assert!(!parse_comp_based_stats("0u", true, false).unwrap().1);
        assert!(parse_comp_based_stats("2", true, true).is_err());
        assert!(parse_comp_based_stats("x", true, true).is_ok());
    }
}
