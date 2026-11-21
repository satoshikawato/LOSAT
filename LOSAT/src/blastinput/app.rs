//! NCBI's application layer shared by the BLASTP, TBLASTN and TBLASTX command lines
//! (`CBlastAppArgs::SetOptions` and the `CATCH_ALL` of the applications): the error
//! frames, the `-outfmt` string (read so by all five programs), the `-seg` and
//! `-comp_based_stats` strings that NCBI reads in its option handlers, after the input and
//! output files are opened.

use std::path::Path;

use crate::blastinput::value_parsers::SegSpec;
use crate::cli::NativeError;

use crate::wording::USER_MANUAL_FOOTER;

/// An error that NCBI raises as `CInputException` (an option handler or `Validate`).
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-176
/// Exception (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`, exception 12): directs users to
/// `USER_MANUAL_FOOTER` ("Please refer to the LOSAT user manual.").
pub fn options_error(message: &str) -> anyhow::Error {
    options_error_bytes(message.as_bytes())
}

/// [`options_error`] of a message's bytes (a matrix name that NCBI's buffer cuts inside a
/// character: `stats::protein_options::print_matrix_message`).
pub fn options_error_bytes(message: &[u8]) -> anyhow::Error {
    NativeError {
        exit: 1,
        message: [
            b"BLAST query/options error: ".as_slice(),
            message,
            format!("\n{USER_MANUAL_FOOTER}\n").as_bytes(),
        ]
        .concat(),
    }
    .into()
}

/// An error that NCBI raises as a `CBlastException` that is neither an option nor a
/// memory error (exit 3).
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
pub fn engine_error(message: &str) -> anyhow::Error {
    NativeError {
        exit: 3,
        message: format!("BLAST engine error: {message}\n").into(),
    }
    .into()
}

/// An error of a reader of NCBI's objtools (`CObjReaderParseException`: the FASTA reader's, the
/// alignment reader's of psiblast's `-in_msa`), which the applications report with the reader's
/// message and exit 1. The one producer of that line.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:181-184
pub fn query_error(message: &str) -> anyhow::Error {
    NativeError {
        exit: 1,
        message: format!("BLAST query error: {message}\n").into(),
    }
    .into()
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275,289-291,300-302 (TBLASTX's
// loop is the same, tblastx_app.cpp:176-207)
/// An error of a query batch's search that comes after the reports of the batches before, for
/// every program: NCBI searches a batch and prints its reports before it reads the next one, so
/// when the search of a batch stops NCBI has printed the reports of the batches before, without
/// the epilog. The error held is the run's error as it is: NCBI's own (`engine_error`, exit 3) or
/// LOSAT's refusal of a step where NCBI crashes, does not finish or stops (exit 1). The batch loop
/// that meets it writes those reports (for the first batch, the prologs), then the failing batch's
/// reader messages and search strategy posts, and returns the error it holds, its message
/// unchanged: TBLASTX's (`tblastx::blast_engine::run_impl`, the X2c form, the Owner's decision of
/// 2026-11-16 on line RG's Q-X2c: "the results of the query batches before are written") and
/// TBLASTN's (`tblastn::tblastn_app::batch_failure`). Any other refusal of a batch's search comes
/// before the run's first byte.
#[derive(Debug)]
pub(crate) struct ErrorAfterBatchReports(pub(crate) anyhow::Error);

impl std::fmt::Display for ErrorAfterBatchReports {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for ErrorAfterBatchReports {}

/// NCBI's error for a search without `-subject` (and without `-db`), which the handler of
/// the database arguments raises before the query and the output are opened.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2558-2562
pub fn missing_subject_error() -> anyhow::Error {
    options_error("Either a BLAST database or subject sequence(s) must be specified")
}

/// NCBI's error for a subject file without records, raised when the subjects are read.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
pub fn empty_subjects_error() -> anyhow::Error {
    engine_error("Empty CBlastQueryVector")
}

/// NCBI's output format numbers (`CFormattingArgs::EOutputFormat`).
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_args.hpp:1023-1074
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

/// A `-outfmt` value after the early reading that every program's `run` makes before its option
/// handlers (`early_format_choice`): the choice it read, or the value left to the formatting
/// handler. A deferred value has no format number and no custom specification, so the type has
/// nothing to read for it: a place that needs the choice takes it from `parsed`, `number` or
/// `at_formatting_handler`, and says what it does for a value not read yet.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EarlyFormat {
    /// A value that NCBI's `ParseFormattingString` read.
    Parsed(FormatChoice),
    /// A value that does not parse, of a run that writes a search strategy file: NCBI reads no
    /// `-outfmt` before its option handlers, and the program's formatting handler reads the
    /// value again and stops with its error (`at_formatting_handler`).
    Deferred,
}

impl EarlyFormat {
    /// The choice that was read, none for a deferred value.
    pub fn parsed(&self) -> Option<&FormatChoice> {
        match self {
            EarlyFormat::Parsed(choice) => Some(choice),
            EarlyFormat::Deferred => None,
        }
    }

    /// The format number that was read, none for a deferred value.
    pub fn number(&self) -> Option<i32> {
        self.parsed().map(|choice| choice.number)
    }

    /// The choice at the program's formatting handler, the last of its option handlers: a
    /// deferred value reads the value `outfmt` again, with NCBI's error (the one owner of where
    /// the value's error comes), any other is as read.
    pub fn at_formatting_handler(&self, outfmt: &str) -> anyhow::Result<FormatChoice> {
        match self {
            EarlyFormat::Parsed(choice) => Ok(choice.clone()),
            EarlyFormat::Deferred => parse_formatting_string(outfmt),
        }
    }

    /// `report_format` of the choice; a deferred value has no report (its error comes at the
    /// formatting handler).
    pub fn report_format(
        &self,
        program: Program,
        out: Option<&Path>,
    ) -> anyhow::Result<Option<ReportFormat>> {
        match self {
            EarlyFormat::Parsed(choice) => report_format(choice, program, out),
            EarlyFormat::Deferred => Ok(None),
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
    pub fn writes_query_warnings(self) -> bool {
        use ReportFormat::*;
        !matches!(
            self,
            AsnText
                | AsnBinary
                | Xml
                | Xml2
                | Json
                | Xml2S
                | JsonS
                | JsonSeqalign
                | Sam
                | ArchiveFormat
        )
    }

    /// Whether the report's writer, as the first format, writes the texts of the subjects that a
    /// query's report reads one by one (where the batch's prefetch is off) where it reads each
    /// (`QueryWarnings::subject_read`): outfmt 0 (description lines and alignment headings), 1 to 4
    /// (description lines, and the subjects with alignments before the section, which NCBI's
    /// merge loads first, showalign.cpp:2003), 5 (`Hit`s), the tabular rows (6, 7, 10 and 20) and
    /// the organism report (18, at each Seq-align it does not skip, `report::tax_report`); the
    /// others write them before the query's report.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-838
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:997
    pub fn places_subject_reads(self) -> bool {
        use ReportFormat::*;
        matches!(
            self,
            Pairwise
                | QueryAnchoredIdentities
                | QueryAnchoredNoIdentities
                | FlatQueryAnchoredIdentities
                | FlatQueryAnchoredNoIdentities
                | Xml
                | Tabular
                | TabularWithComments
                | CommaSeparatedValues
                | TaxFormat
                | CommaSeparatedValuesWithHeader
        )
    }
}

/// LOSAT's search programs, for what depends on the program in the reading of `-outfmt`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastn_args.cpp:104
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastp_args.cpp:95 (c++/src/algo/blast/blastinput/blastx_args.cpp:113,
/// tblastn_args.cpp:110 and tblastx_args.cpp:104 are the same)
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

    /// Whether LOSAT's program writes custom tabular field lists (BLASTN, BLASTP and BLASTX
    /// write them with `report::tabular`).
    pub fn custom_fields(self) -> bool {
        self.shared_tabular()
    }

    /// Whether the program writes its tabular reports with a field list or a delimiter, and
    /// outfmt 10 and 20, through the shared `report::tabular` (every field, `delim=`): BLASTN
    /// (line CLI C2, PR 1), TBLASTN and TBLASTX (PR 2), BLASTP, with psiblast, which runs
    /// BLASTP's engine, and BLASTX (PR 3).
    pub fn shared_tabular(self) -> bool {
        matches!(
            self,
            Program::Blastn
                | Program::Blastp
                | Program::Blastx
                | Program::Tblastn
                | Program::Tblastx
        )
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
            // Line RG T3 (the query-anchored formats 1 to 4, `report::query_anchored`): BLASTP
            // (T3a), BLASTN and TBLASTN (T3b), BLASTX and TBLASTX (T3c): every program.
            QueryAnchoredIdentities
            | QueryAnchoredNoIdentities
            | FlatQueryAnchoredIdentities
            | FlatQueryAnchoredNoIdentities => true,
            // 10 and 20 with the shared tabular writer (line CLI C2).
            CommaSeparatedValues => self.shared_tabular(),
            // Line RG T2: the organism report, which reads the hit lists as the shared tabular
            // writer does (`report::tax_report`).
            TaxFormat => self.shared_tabular(),
            CommaSeparatedValuesWithHeader => self.shared_tabular(),
            // Line OF track S (structured formats): 5 for BLASTP and BLASTN (OF4a), BLASTX
            // (OF4b part 1), TBLASTN (OF4b part 3) and TBLASTX (OF4b part 4): every program;
            // 15 and 16 (JSON_S, XML2_S) for every program: BLASTP (OF5a part 1), BLASTN (OF5a part 2),
            // BLASTX (OF5b part 1), TBLASTN (OF5b part 2) and TBLASTX (OF5b part 3).
            Xml => true,
            AsnText => true,
            AsnBinary => true,
            // OF8a: the BLAST archive for BLASTN (part 1) and BLASTP (part 2); OF8b: BLASTX
            // (part 1), TBLASTN (part 2) and TBLASTX (part 3): every program.
            ArchiveFormat => true,
            JsonSeqalign => true,
            // 13 and 14 (BLAST JSON and XML2 as side files and a master file): every program
            // (BLASTP and BLASTN: OF5c part 1; BLASTX, TBLASTN, TBLASTX: part 2).
            Json => true,
            Xml2 => true,
            JsonS => true,
            Xml2S => true,
            // Line OF track S: SAM for BLASTN alone (OF6).
            Sam => self == Program::Blastn,
        }
    }
}

/// NCBI's reading of a `-outfmt` value, which `SetOptions` does before every other handler
/// (`ArchiveFormatRequested`), with NCBI's errors.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2801-2851
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:260-263
/// `NStr::TruncateSpaces` (`eTrunc_Both`) is `ncbi_str::truncate_spaces_str`; `NStr::StringToInt`
/// is `ncbi_str::string_to_int` (its digits end at a NUL, as NCBI's do).
pub fn parse_formatting_string(spec: &str) -> anyhow::Result<FormatChoice> {
    use crate::ncbi_str::truncate_spaces_str;
    let choice = truncate_spaces_str(spec);
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
        custom = truncate_spaces_str(&custom.replace(&token, "")).to_string();
    }
    let Ok(number) = crate::ncbi_str::string_to_int(choice.as_bytes(), 0) else {
        return Err(options_error(&format!(
            "'{choice}' is not a valid output format"
        )));
    };
    if !(0..END_VALUE).contains(&number) {
        return Err(NativeError {
            exit: 255,
            message: "Error: Formatting choice is out of range\n".into(),
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
/// (`Program::custom_fields`), and a custom delimiter unless it writes them with the shared
/// tabular writer (`Program::shared_tabular`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2874-2886
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:891-893
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
    // SAM reads its own words and ignores a delimiter (`CBlast_SAM_Formatter::x_ProcessCustomSpec`).
    let is_sam = format == ReportFormat::Sam;
    if choice.custom() && !program.custom_fields() && !is_sam {
        anyhow::bail!(
            "a custom output format specification (fields or a delimiter) is not supported by LOSAT's {name}"
        );
    }
    if !choice.delimiter.is_empty() && !is_sam && !program.shared_tabular() {
        anyhow::bail!(
            "a custom delimiter (delim=) in the output format is not supported by LOSAT's {name}"
        );
    }
    Ok(Some(format))
}

/// NCBI's errors for the formats of other programs, which the formatting handler raises
/// after the input and output files are opened.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2874-2886
///
/// Then NCBI's check of a custom delimiter (`report::tabular::custom_delimiter_conflict`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2887
pub fn formatting_handler_check(choice: &FormatChoice, program: Program) -> anyhow::Result<()> {
    match choice.number {
        SAM if !program.is_sam() => Err(options_error("SAM format is only applicable to blastn")),
        AIRR_REARRANGEMENT => Err(options_error(
            "AIRR rearrangement format is only applicable to igblastn",
        )),
        FASTA => Err(options_error(
            "FASTA output format is only applicable to magicblast",
        )),
        _ => {
            match crate::report::tabular::custom_delimiter_conflict(&choice.spec, &choice.delimiter)
            {
                Some(message) => Err(options_error(&message)),
                None => Ok(()),
            }
        }
    }
}

/// NCBI's error for the XInclude formats (13 and 14) written to standard output, after
/// `Query is Empty!`.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:888-899
pub fn xinclude_check(choice: &FormatChoice) -> anyhow::Result<()> {
    xinclude_check_for(choice, None)
}

/// `xinclude_check` for a run with `-out` `out`: the XInclude formats written to a file are
/// not an error (`UseXInclude` returns true and the formatter takes the file name as the base of
/// its side files).
pub fn xinclude_check_for(choice: &FormatChoice, out: Option<&Path>) -> anyhow::Result<()> {
    if out.is_some_and(|path| path.as_os_str() != "-") {
        return Ok(());
    }
    match choice.number {
        JSON => Err(options_error("Please provide a file name for outfmt 13.")),
        XML2 => Err(options_error("Please provide a file name for outfmt 14.")),
        _ => Ok(()),
    }
}

/// Whether NCBI refuses the format choice `choice` of a `program` run with `-out` `out` itself,
/// later than its parse: SAM, AIRR and FASTA outside their programs and a custom delimiter that
/// the field list uses (`formatting_handler_check`), and an XInclude format (13, 14) without a
/// file name (`xinclude_check_for`). LOSAT's refusals of the formats that psiblast, rpsblast and
/// rpstblastn do not write yet, which come at the parse, leave these to NCBI's own error at its
/// place (line RG F/d36 probes rps.conf20, rps.o14, ps.o13, ps.o14).
pub fn ncbi_refuses_format(choice: &FormatChoice, program: Program, out: Option<&Path>) -> bool {
    formatting_handler_check(choice, program).is_err() || xinclude_check_for(choice, out).is_err()
}

/// The early reading of the `-outfmt` value `outfmt` that every program's `run` makes before its
/// option handlers (NCBI's `ArchiveFormatRequested`), with NCBI's errors: the choice, or
/// `EarlyFormat::Deferred` where the value does not parse and the run writes a search strategy
/// file (`strategy`). NCBI then reads no `-outfmt` before its option handlers: the strategy file
/// is created first and the early read is skipped (the `||` of that reference), so the value's error comes
/// from the formatting handler, the last of the handlers, after the arguments' own errors in
/// their order (`-subject` opened and `-subject_loc` read, `-query` opened, `-out` created,
/// `-query_loc` read; line RG F/d36 review 2c B1). The program's own formatting handler reads the
/// value again and raises it there; nothing before it can use a deferred choice (`EarlyFormat`
/// has no number for it).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
pub fn early_format_choice(outfmt: &str, strategy: bool) -> anyhow::Result<EarlyFormat> {
    match parse_formatting_string(outfmt) {
        Ok(choice) => Ok(EarlyFormat::Parsed(choice)),
        Err(_) if strategy => Ok(EarlyFormat::Deferred),
        Err(error) => Err(error),
    }
}

/// A `-seg` value as NCBI's filtering handler reads it, with NCBI's errors. Every locut and
/// hicut that `NStr::StringToDouble` reads is taken, an infinity or NaN too: the query filter
/// keeps the default for one that is not above 0 (NaN included, `SegSpec::params`), SEG
/// compares an infinite one as NCBI does (`+inf` masks every window), and the filter string of
/// the reports prints it as glibc does (`report::blast_xml::filter`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:396-406
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:375-384
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:423-428
pub fn parse_seg_option(value: &str) -> anyhow::Result<SegSpec> {
    read_seg_option(value).map_err(options_error)
}

/// `parse_seg_option` with NCBI's error as its message (a `CInputException`), for a program that
/// reports it with its other deferred errors (BLASTX: `algorithm::blastx::args`).
pub fn read_seg_option(value: &str) -> Result<SegSpec, &'static str> {
    match value {
        "no" => return Ok(SegSpec::No),
        "yes" => return Ok(SegSpec::Yes),
        _ => {}
    }
    let tokens: Vec<&str> = value.split(' ').collect();
    if tokens.len() != 3 {
        return Err("Invalid number of arguments to filtering option");
    }
    let invalid = || "Invalid input for filtering parameters";
    // `NStr::StringToInt(tokens[0])`, flags 0.
    let window = crate::ncbi_str::string_to_int(tokens[0].as_bytes(), 0).map_err(|_| invalid())?;
    // NCBI reads the values with `NStr::StringToDouble(token)` (default flags,
    // `ncbi_str::string_to_double`: glibc's `strtod`, decimal or hexadecimal, which must end at
    // the end of the token).
    let mut cuts = [0.0; 2];
    for (cut, token) in cuts.iter_mut().zip(&tokens[1..]) {
        *cut = crate::ncbi_str::string_to_double(token.as_bytes()).map_err(|_| invalid())?;
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:834-892
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

/// The error of the query splitter that NCBI sets up for every query batch (`CQuerySplitter`, made
/// with the setup data at the start of `BlastSetupPreliminarySearchEx`), which reads `CHUNK_SIZE`
/// and then (where it can split) `OVERLAP_CHUNK_SIZE` but does not split an ungapped search. The
/// first query batch meets it after the reading of its records (their reader errors and `Empty
/// CBlastQueryVector` come first) and before its setup (a query without letters comes after): a
/// value that `NStr::StringToInt` cannot convert is its `CStringException`
/// ([`environment_int`]), and a chunk size that is not divisible by 3 is NCBI's engine error. Each
/// program returns it through the error frame of its first batch, after the outfmt 0 prolog.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:55-61
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.hpp:104-106
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:98-103
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:78-111
/// `SplitQuery_CalculateNumChunks` reads the overlap (`getenv("OVERLAP_CHUNK_SIZE")` in
/// `SplitQuery_GetOverlapChunkSize`) only after `SplitQuery_ShouldSplit` is true, and that is
/// false for RPS-BLAST (`Blast_SubjectIsPssm`, "do not split RPS-BLAST"), so `rpsblast` reads
/// `CHUNK_SIZE` (the constructor above) but never `OVERLAP_CHUNK_SIZE`: `reads_overlap` is false
/// for it and an overlap that is not an integer does not stop the run. It is true for blastp,
/// psiblast, tblastn and tblastx (`blastn/query_split.rs` `SplitSizes`; the `int` becomes a
/// 64-bit `size_t`). The chunk size must be divisible by 3 only for a translated query
/// (`translated_query`, tblastx). Whether a blastp or tblastn query batch is split is
/// `blastp/protein_query_split.rs`.
pub fn check_query_split_environment(
    translated_query: bool,
    reads_overlap: bool,
) -> anyhow::Result<()> {
    if let Some(chunk_size) = read_query_split_variable("CHUNK_SIZE")? {
        if translated_query && !(i64::from(chunk_size) as u64).is_multiple_of(3) {
            return Err(engine_error(
                "Split query chunk size must be divisible by 3",
            ));
        }
    }
    if reads_overlap {
        read_query_split_variable("OVERLAP_CHUNK_SIZE")?;
    }
    Ok(())
}

/// `CHUNK_SIZE` or `OVERLAP_CHUNK_SIZE` as NCBI reads it (`None` when unset or blank); `Err` is
/// the `CStringException` of a value that it cannot convert ([`environment_int`]).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:58-61
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:54-57
pub(crate) fn read_query_split_variable(variable: &str) -> anyhow::Result<Option<i32>> {
    match std::env::var_os(variable) {
        Some(value) if !crate::ncbi_str::is_blank(value.as_encoded_bytes()) => {
            environment_int(&value).map(Some)
        }
        _ => Ok(None),
    }
}

/// An environment variable's value as NCBI's search converts it, `NStr::StringToInt(str)` with
/// flags 0 (`ncbi_str::string_to_int`: an optional sign and decimal digits within the range of an
/// `int`, no spaces). A value that it cannot convert (an empty one too) throws its
/// `CStringException`, which ends the program through `CATCH_ALL` (`NativeError::exception`: `Error:
/// ` and the exception's report, exit status 255). The bytes are those of the value as the
/// process got it.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:635-643 (cited at
/// `ncbi_str::string_to_int`)
pub(crate) fn environment_int(value: &std::ffi::OsStr) -> anyhow::Result<i32> {
    crate::ncbi_str::string_to_int(value.as_encoded_bytes(), 0)
        .map_err(|e| NativeError::exception(&e).into())
}

/// The environment that changes NCBI's search or report in a way that LOSAT does not
/// reproduce: `BL2SEQ_LEGACY` in a search of `-subject` sequences (`fasta_subjects`; NCBI
/// searches each subject on its own and writes another report, a `-db` search ignores it).
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:206-210
/// `BL2SEQ_LEGACY` (any value, also an empty one) turns off the scan of the subjects as a
/// database: `CLocalBlast` then searches each subject on its own (local_blast.cpp:189,289)
/// and `CBlastFormat` writes the legacy bl2seq report (`m_IsBl2Seq && !m_IsDbScan`).
///
/// LOSAT refuses before the search, so the refusal is wider than NCBI's: a `-subject` search
/// with `BL2SEQ_LEGACY` whose report would not change (no hits, or one subject record in outfmt
/// 5, 6 or 10: NCBI 2.17.0 writes the same bytes as without the variable there) is refused too.
pub fn check_search_environment(program: &str, fasta_subjects: bool) -> anyhow::Result<()> {
    if fasta_subjects && std::env::var_os("BL2SEQ_LEGACY").is_some() {
        anyhow::bail!(
            "the environment variable BL2SEQ_LEGACY, which makes NCBI BLAST+ search each subject on its own and write its legacy bl2seq report, is not supported by LOSAT's {program}"
        );
    }
    Ok(())
}

/// The error of the prefetch of the report's sequences before the reports of a query batch
/// (`BlastFormatter_PreFetchSequenceData`): `PRE_FETCH_SEQS_LIMIT` (any set value, an empty one
/// too) is converted by `NStr::StringToInt`, and a value that it cannot convert stops NCBI with its
/// `CStringException` ([`environment_int`]). The first batch meets it after its search (whose
/// errors come first) and before anything of its report (the query warnings included); each
/// program returns it through the error frame of that batch, after the outfmt 0 prolog and the
/// reader's messages of the batch's records. An integer only decides whether the sequences are
/// fetched ahead (`common::database::prefetch_off`).
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:732-737
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:756-765
/// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:283-288
/// The archive (`archive`: every requested format is outfmt 11, `ReportOutputs::only_archive`) is
/// written without the prefetch, so NCBI never converts the value there; blastp, tblastn,
/// tblastx, rpsblast and blastn write it so too (blastp_app.cpp:280-285, tblastn_app.cpp:294-299
/// and 330-335, tblastx_app.cpp:197-202, rpsblast_app.cpp:217-222, blastn_app.cpp:304-309), and
/// so do the nodes of a search split by queries (blastp_node.cpp:178-182, blastx_node.cpp:182-186,
/// tblastn_node.cpp:206-210 and 243-247, blastn_node.cpp:204-208). psiblast's iterations fetch
/// ahead before they choose the archive, so PSIBLAST passes `archive` false.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:262-266
/// When NCBI splits a `-db` search by queries (blastp, blastx, tblastn, blastn: `-mt_mode 1`, or
/// the automatic choice at `-num_threads` 2 or more for a query file over the program's size and a
/// database under its length, `s_GetMT1Cutoffs`), its nodes catch the exception: NCBI writes
/// nothing to standard output and exits 0; LOSAT gives NCBI's `-num_threads 1` output at any thread
/// count (env-parity DECISIONS 11).
pub fn check_prefetch_environment(archive: bool) -> anyhow::Result<()> {
    if archive {
        return Ok(());
    }
    prefetch_seqs_limit().map(drop)
}

/// `PRE_FETCH_SEQS_LIMIT` as NCBI's prefetch converts it (`s_PreFetchSeqs`, cited at
/// `check_prefetch_environment`): `None` when it is not set; `Err` is the `CStringException` of a
/// value that `NStr::StringToInt` cannot convert ([`environment_int`]). The one reader of the
/// variable: the check of the first batch's prefetch, and the limit that decides whether a batch's
/// prefetch is on (`common::database::prefetch_off`).
pub(crate) fn prefetch_seqs_limit() -> anyhow::Result<Option<i32>> {
    std::env::var_os("PRE_FETCH_SEQS_LIMIT")
        .map(|value| environment_int(&value))
        .transpose()
}

/// The query batch size of blastp, tblastn or tblastx (`GetQueryBatchSize`, with the
/// program's `default`: tblastx 10002 nucleotides, tblastn 20000 and blastp 10000
/// residues), as the `TSeqPos` that `CBlastInput` compares with.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:117-119
/// A value that `NStr::StringToInt` cannot convert (an empty one too) stops NCBI there, before the
/// formatter and anything else is written, with its `CStringException` (`NativeError::exception`,
/// exit status 255; with outfmt 11 NCBI's error archive follows).
pub fn query_batch_size(default: u32) -> anyhow::Result<u32> {
    let value = std::env::var_os("BATCH_SIZE");
    match crate::blastinput::query_batch::get_query_batch_size(value.as_deref()) {
        Ok(_) if value.is_none() => Ok(default),
        // `CBlastInput` stores the `int` as a `TSeqPos` (blast_input.hpp:313,364).
        Ok(size) => Ok(size as u32),
        Err(e) => Err(NativeError::exception(&e).into()),
    }
}

/// Whether the environment variable `OLD_FSC` is set (to any value, the empty one too), with
/// which NCBI's protein programs compute without the Gumbel block (finite-size correction).
/// NCBI sets it itself for `-searchsp` (c++/src/algo/blast/blastinput/blast_args.cpp:303-307), so BLASTP, PSIBLAST, RPSBLAST, BLASTX
/// and TBLASTN read it as they read `-searchsp`: their `old_fsc`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:922-923
pub fn old_fsc_environment() -> bool {
    std::env::var_os("OLD_FSC").is_some()
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-57
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
#[path = "app_tests.rs"]
mod tests;
