//! The tabular reports with a field list (outfmt 6, 7, 10 and 20 with `-outfmt "<n> <fields>"`
//! or `delim=`): NCBI's `CBlastTabularInfo`, shared by the programs. A program's adapter hands
//! over, for each HSP of a query's final hit list in report order, the values of the HSP
//! (`TabularHsp`), its subject (`TabularSubject`) and its query (`TabularQuery`); this module
//! reads the field list, writes the comment lines of outfmt 7, the header line of outfmt 20 and
//! the rows.
//!
//! `TabularSpec::of_app` returns `None` for the default reports (outfmt 6 and 7 with the `std`
//! fields and a tab): BLASTX writes them here too (`TabularSpec::std`), the other programs with
//! their own writers.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:770-786

use std::io::{self, Write};

use crate::blastinput::app::{FormatChoice, ReportFormat};
use crate::common::GapEditOp;
use crate::config::ScoringMatrix;
use crate::report::outfmt6::{
    format_bitscore_ncbi, format_evalue_ncbi_tabular, format_percent_identity_ncbi,
    write_outfmt7_query_line,
};
use crate::report::serial::seqalign::{row_end_gaps, AlignRow};

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
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:44-198
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1111-1236
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
pub fn parse_fields(spec: impl AsRef<[u8]>) -> Vec<TabularField> {
    let mut fields: Vec<TabularField> = Vec::new();
    let add = |fields: &mut Vec<TabularField>, field| {
        if !fields.contains(&field) {
            fields.push(field);
        }
    };
    for token in spec.as_ref().split(|&b| b == b' ') {
        let Ok(token) = std::str::from_utf8(token) else {
            continue;
        };
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
/// Returns NCBI's message; the field list (`spec`) is searched as text, as NCBI does.
pub fn custom_delimiter_conflict(
    spec: impl AsRef<[u8]>,
    delimiter: impl AsRef<[u8]>,
) -> Option<String> {
    let spec = crate::ncbi_str::truncate_spaces(spec.as_ref());
    if spec.is_empty() {
        return None;
    }
    let delimiter = delimiter.as_ref();
    let contains = |needle: &str| spec.windows(needle.len()).any(|w| w == needle.as_bytes());
    let (shown, field) = match delimiter {
        b";" => (
            ";",
            [
                "sallseqid",
                "staxids",
                "sscinames",
                "scomnames",
                "sblastnames",
                "sskingdoms",
            ]
            .into_iter()
            .find(|field| contains(field)),
        ),
        b"/" => ("/", Some("frames").filter(|field| contains(field))),
        b"<>" => ("<>", Some("salltitles").filter(|field| contains(field))),
        _ => return None,
    };
    let field = field?;
    Some(format!(
        "Your custom record separator ({shown}) is also used by the format specifier ({field}) to separate multiple entries. Please use a different record separator (delim keyword)."
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
    pub delimiter: Vec<u8>,
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
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:245-249
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:501-505
    // NCBI reference (598d8ae6): c++/src/app/blast/rpsblast_app.cpp:186-189
    // NCBI reference (598d8ae6): c++/src/app/blast/rpstblastn_app.cpp:183-186
    /// The tabular report of `format` with the custom specification of `choice` as the formatter
    /// of `app`'s application writes it, or `None` for a format that is not tabular and for the
    /// default outfmt 6 and 7 (the `std` fields and a tab), which BLASTX writes with `std`'s spec
    /// and the other programs with their own writers. The applications of blastn, blastp, blastx, tblastn (PSI-TBLASTN too) and
    /// tblastx hand `-outfmt`'s custom delimiter to their formatter (`SetCustomDelimiter`,
    /// blastn_app.cpp:254, blastx_app.cpp:250, tblastn_app.cpp:266, tblastx_app.cpp:169);
    /// psiblast's, rpsblast's and rpstblastn's never do, so their reports keep the format's own
    /// delimiter, the tab or the comma (blast_format.cpp:783-785 above), and a `delim=` with the
    /// `std` fields leaves the default outfmt 6 and 7 (`None`). Their `delim=` is still checked
    /// against the fields (`custom_delimiter_conflict`). The one reading of a tabular `-outfmt`
    /// value, of the command line and the library, for every program: no other function makes a
    /// spec of a format choice.
    pub fn of_app(
        format: ReportFormat,
        choice: &FormatChoice,
        app: crate::report::archive::options::RequestKind,
    ) -> Option<Self> {
        use crate::report::archive::options::RequestKind;
        let delimiter: &[u8] = match app {
            RequestKind::Psiblast
            | RequestKind::Phiblastp
            | RequestKind::Rpsblast
            | RequestKind::Rpstblastn => b"",
            RequestKind::Blastn
            | RequestKind::Blastp
            | RequestKind::Blastx
            | RequestKind::Tblastn
            | RequestKind::PsiTblastn
            | RequestKind::Tblastx => &choice.delimiter,
        };
        Self::of_fields(format, &choice.spec, delimiter)
    }

    /// [`TabularSpec::of_app`] with the delimiter `delimiter` (empty: the format's own).
    fn of_fields(format: ReportFormat, spec: &[u8], delimiter: &[u8]) -> Option<Self> {
        let (kind, default_delimiter): (_, &[u8]) = match format {
            ReportFormat::Tabular => (TabularKind::Tabular, b"\t"),
            ReportFormat::TabularWithComments => (TabularKind::Commented, b"\t"),
            ReportFormat::CommaSeparatedValues => (TabularKind::Csv, b","),
            ReportFormat::CommaSeparatedValuesWithHeader => (TabularKind::CsvWithHeader, b","),
            // NCBI's `ParseFormattingString` keeps no field list for outfmt 18.
            ReportFormat::TaxFormat => {
                return Some(TabularSpec {
                    kind: TabularKind::TaxReport,
                    fields: Vec::new(),
                    delimiter: Vec::new(),
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
        .to_vec();
        let default_report = matches!(kind, TabularKind::Tabular | TabularKind::Commented)
            && fields == STD_FIELDS
            && delimiter == b"\t";
        (!default_report).then_some(TabularSpec {
            kind,
            fields,
            delimiter,
        })
    }

    /// The default report of outfmt 6 or 7 (`commented`): the `std` fields and a tab, for which
    /// `of_app` gives `None` (BLASTX writes it with this module as any field list,
    /// `algorithm::blastx::args::ResolvedOptions::tabular_report`). Made once per process.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:78-81
    pub fn std(commented: bool) -> &'static TabularSpec {
        static TABULAR: std::sync::OnceLock<&'static TabularSpec> = std::sync::OnceLock::new();
        static COMMENTED: std::sync::OnceLock<&'static TabularSpec> = std::sync::OnceLock::new();
        let (cell, kind) = if commented {
            (&COMMENTED, TabularKind::Commented)
        } else {
            (&TABULAR, TabularKind::Tabular)
        };
        cell.get_or_init(|| {
            TabularSpec {
                kind,
                fields: STD_FIELDS.to_vec(),
                delimiter: b"\t".to_vec(),
            }
            .intern()
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

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:973-980
    /// Whether NCBI's formatter reads the query's and the subject's letters for a row of this field
    /// list (`CAlnVec::GetWholeAlnSeqString`) in psiblast, whose formatter fetches sequences
    /// (`kNoFetchSequence` is false): `qseq`, `sseq`, `positive`, `ppos`, `btop`, `nident`,
    /// `mismatch` and `pident`, and so the default fields.
    pub fn reads_sequences(&self) -> bool {
        self.has_any(&[
            QuerySeq,
            SubjectSeq,
            Positives,
            PercentPositives,
            Btop,
            NumIdentical,
            Mismatches,
            PercentIdentical,
        ])
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:833-837,973-980 and 973-980
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-908,915,1078-1082
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

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:891-924
    // NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:964-968
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:996-1000
    /// Whether an RPS program refuses its organism report (outfmt 18) or comma-separated values
    /// with a header line (outfmt 20) of these fields for an HSP with the Seq-align's identities
    /// `num_ident` whose profile's Seq-align id is the generated id (`lcl|Query_N`) of a query of
    /// its batch (`common::database::refuse_generated_query_id`; `program` is the driver's:
    /// BLASTP for RPSBLAST, BLASTX for RPSTBLASTN). NCBI's formatter finds that query by the id,
    /// so its row differs where it shows what it found beyond the id: the found sequence's ids,
    /// titles, taxids and length (`sallseqid`, `sallacc`, `stitle`, `salltitles`, `staxids`:
    /// `N/A` for a query, `slen`), the subject's letters of the alignment (`sseq`, `btop`,
    /// `positive`, `ppos`, `mismatch`, `pident`, and `nident` where the block opens for it), and
    /// the taxonomy (`staxid`, the names and the organism report) where the profile has a taxid
    /// (`subject_taxid`; a query has none, as a profile without one). In RPSTBLASTN, whose found
    /// query is nucleotide, every row that enters the alignment block (`opens_counts`) ends the
    /// run with an exception of the Seq-align's conversion after the header line. The profile's
    /// own ids (`sseqid`, `sacc`, `saccver`, `sgi`, `sallgi`), the query's letters and the
    /// strand of RPSBLAST, the scores and coordinates are NCBI's bytes (review S1 of PS15a-2:
    /// every field alone, `tmp/rev15a2/col`). Outfmt 6, 7 and 10 of such an HSP are line RG
    /// item X7's.
    pub fn refuses_generated_id(
        &self,
        program: crate::blastinput::app::Program,
        num_ident: usize,
        subject_taxid: bool,
    ) -> bool {
        match self.kind {
            TabularKind::TaxReport => subject_taxid,
            TabularKind::CsvWithHeader => {
                let found = self.has_any(&[
                    SubjectAllSeqIds,
                    SubjectAllAccessions,
                    SubjectTitle,
                    SubjectAllTitles,
                    SubjectTaxIds,
                    SubjectLength,
                ]);
                let letters = self.has_any(&[
                    SubjectSeq,
                    Btop,
                    Positives,
                    PercentPositives,
                    Mismatches,
                    PercentIdentical,
                ]) || (self.has(NumIdentical) && num_ident > 0);
                let taxonomy = subject_taxid
                    && self.has_any(&[
                        SubjectTaxId,
                        SubjectSciName,
                        SubjectCommonName,
                        SubjectBlastName,
                        SubjectSuperKingdom,
                        SubjectSciNames,
                        SubjectCommonNames,
                        SubjectBlastNames,
                        SubjectSuperKingdoms,
                    ]);
                let throws = program == crate::blastinput::app::Program::Blastx
                    && self.opens_counts(&num_ident);
                found || letters || taxonomy || throws
            }
            TabularKind::Tabular | TabularKind::Commented | TabularKind::Csv => false,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:833-878
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203 (the `CATCH_ALL` of
    // every application, cited at `cli::NativeError::exception`)
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

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:979-1034
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:1050-1055
    /// Whether `program`'s report of this field list writes what it writes with a readable
    /// taxonomy database where every load of a database sequence with a positive taxid and every
    /// taxid lookup throw (a file of that database does not open, `SeqDbError::TaxdbThrows`): a
    /// table whose rows load no subject's Bioseq ([`TabularSpec::subject_fetch`]; `slen` is read
    /// without it), and the organism report, whose `catch` leaves out each HSP whose load or lookup
    /// throws, as it leaves out one whose taxid is not found, as LOSAT does
    /// (`tax_report::TaxLookups`; such a lookup finds nothing, and every load that throws is of a
    /// subject with a positive taxid, whose lookup throws too).
    pub fn unchanged_by_load_throws(&self, program: crate::blastinput::app::Program) -> bool {
        self.kind == TabularKind::TaxReport || self.subject_fetch(program) == SubjectFetch::Never
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:995-1023 (cited at
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
            FormatThrow::Writes(RunThrow {
                first_row_ends: first_thrown_row_ends_run(Some(self)),
                ..RunThrow::every_row(AfterThrow::WritesRow)
            })
        } else if self.reads_subject(program) {
            FormatThrow::Writes(RunThrow {
                rows: self.throw_rows(program),
                ..RunThrow::every_row(AfterThrow::EndsRun)
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
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:126-135
    pub fn checks_taxdb(&self) -> bool {
        self.has_any(&[
            SubjectSciNames,
            SubjectCommonNames,
            SubjectBlastNames,
            SubjectSuperKingdoms,
        ])
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:139-158
    /// The warning that the formatter NCBI makes for each query's report of this field list
    /// writes when it is made (`x_CheckTaxDB`: `checks_taxdb` and no `taxdb.bti` found,
    /// `taxdb_found`), as `program` (lower case) posts it: the one owner of that rule for the
    /// report writers and for what a refused report writes before its refusal
    /// (`FailedLookupRefusal::formatter_warning`).
    pub fn formatter_warning(&self, taxdb_found: bool, program: &str) -> Option<String> {
        (self.checks_taxdb() && !taxdb_found).then(|| taxdb_warning_line(program))
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
    /// subject, whose load may throw, `FailedStage::Bioseq`), or `sallacc` (a subject with a
    /// definition line without a Seq-id, `LookupFailure::lineless`). A report that cannot stop
    /// needs no check of its rows before it is written.
    pub fn may_stop_at_failed_lookup(&self) -> bool {
        self.writes_rows()
            && self.has_any(&[
                SubjectAccession,
                SubjectAccessionVersion,
                SubjectTaxId,
                SubjectAllAccessions,
            ])
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
        self.delimiter == b","
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:174-185
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:262-272
/// Why NCBI has no result for a tabular report of `spec` (`None`: the default fields, whose
/// `saccver` is the first case) with a row whose subject lookups fail as `failure` says, the row
/// being the first of its result set when `first`: `sacc` and `saccver` of the empty
/// `m_SubjectId` (`SetSubjectId` threw) take the label of a null `FindBestChoice` (NCBI BLAST+
/// 2.17.0 stops with `CCoreException::eNullPtr`, exit 255, after the row's earlier fields), and
/// `staxid` of a first row prints `m_SubjectTaxId` before any row set it (an uninitialised `int`:
/// a different number in each run). `None`: the report has a result (`TabularWriter::row`).
///
/// A subject with a definition line without a Seq-id (`failure.lineless`) and no other failure
/// has rows whose `sallacc` takes the label of that line's best id, null: NCBI's
/// `x_PrintSubjectAllAccessions` has no `catch` and stops (exit 255) after the row's earlier fields,
/// at any row (`LINELESS_ACCESSION_NULL_POINTER`).
pub fn failed_lookup_stop(
    spec: Option<&TabularSpec>,
    failure: LookupFailure,
    first: bool,
) -> Option<&'static str> {
    let Some(spec) = spec else {
        return (failure.subject_id || (failure.bioseq && first))
            .then_some(FAILED_LOOKUP_NULL_POINTER);
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:174-185
    if failure.lineless
        && !failure.bioseq
        && spec.writes_rows()
        && spec.has(SubjectAllAccessions)
        && spec.failure_stage(failure).is_none()
    {
        return Some(LINELESS_ACCESSION_NULL_POINTER);
    }
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

const FAILED_LOOKUP_NULL_POINTER: &str = "NCBI BLAST+ stops with a null pointer exception (exit 255) at the sacc or saccver of its row, which reads the accession of the subject's Bioseq that its load did not give (the default fields of outfmt 6, 7, 10 and 20 print saccver), and reports the addresses of its own binary";

/// The reason of `failed_lookup_stop` for the `sallacc` of a subject with a definition line
/// without a Seq-id (`LookupFailure::lineless`).
pub const LINELESS_ACCESSION_NULL_POINTER: &str = "NCBI BLAST+ stops with a null pointer exception (exit 255) at the sallacc of each row of such a sequence";

/// Whether NCBI's run of a tabular report of `spec` (`None`: the default fields) ends where the
/// first row of a query's report loads a subject whose load throws (`failed_lookup_stop` of
/// `LookupFailure::thrown`: the null pointer of `sacc` or `saccver`, exit 255), for
/// `report::db_subject::RunThrow::first_row_ends`.
pub fn first_thrown_row_ends_run(spec: Option<&TabularSpec>) -> bool {
    failed_lookup_stop(spec, LookupFailure::thrown(), true) == Some(FAILED_LOOKUP_NULL_POINTER)
}

/// Whether a tabular report of `spec` (`None`: the default fields) can stop where NCBI has no
/// result (`TabularSpec::may_stop_at_failed_lookup`), so that its rows are checked
/// (`check_failed_lookups`) before it is written.
pub fn may_stop_at_failed_lookup(spec: Option<&TabularSpec>) -> bool {
    spec.is_none_or(TabularSpec::may_stop_at_failed_lookup)
}

/// Whether a tabular report of `spec` (`None`: the default fields) can stop at a row whose
/// subject lookup throws (`sacc`, `saccver`, default fields, or `staxid`), apart from a
/// lineless-only definition line (`stops_at_lineless_accession`).
pub fn stops_at_lookup_failure(spec: Option<&TabularSpec>) -> bool {
    spec.is_none_or(|spec| {
        spec.writes_rows()
            && spec.has_any(&[SubjectAccession, SubjectAccessionVersion, SubjectTaxId])
    })
}

/// Whether a tabular report of `spec` (`None`: the default fields) stops at a subject that has a
/// definition line without a Seq-id (`LookupFailure::lineless`, `failed_lookup_stop`: `sallacc`).
pub fn stops_at_lineless_accession(spec: Option<&TabularSpec>) -> bool {
    spec.is_some_and(|spec| spec.writes_rows() && spec.has(SubjectAllAccessions))
}

/// The check, before a tabular report of `spec` (`None`: the default fields) writes the rows of
/// `result_sets` (each the indices of its subjects in row order, a subject's run of rows once,
/// one set per query or psiblast round), that NCBI has a result for each row whose subject
/// lookups fail (`failure`, `None` where they do not): `failed_lookup_stop`'s cases are LOSAT's
/// refusals, which name the subject (`describe` of its index and failure: what its title or
/// definition line is) and the program (upper case; its formatter's warning, `taxdb_found` as
/// `TabularSpec::formatter_warning` reads it, names it in lower case).
pub fn check_failed_lookups<S>(
    spec: Option<&TabularSpec>,
    result_sets: impl IntoIterator<Item = S>,
    failure: impl Fn(usize) -> Option<LookupFailure>,
    describe: impl Fn(usize, LookupFailure) -> String,
    program: &str,
    taxdb_found: bool,
) -> anyhow::Result<()>
where
    S: IntoIterator<Item = usize>,
{
    if !may_stop_at_failed_lookup(spec) {
        return Ok(());
    }
    let mut rows = Vec::new();
    for (index, set) in result_sets.into_iter().enumerate() {
        rows.clear();
        for (row, s_idx) in set.into_iter().enumerate() {
            rows.push(s_idx);
            let Some(failure) = failure(s_idx) else {
                continue;
            };
            if let Some(reason) = failed_lookup_stop(spec, failure, row == 0) {
                return Err(FailedLookupRefusal {
                    message: format!(
                        "{}; {reason}, which is not supported by LOSAT's {program}",
                        describe(s_idx, failure)
                    ),
                    set: index,
                    rows,
                    formatter_warning: spec.and_then(|spec| {
                        spec.formatter_warning(taxdb_found, &program.to_ascii_lowercase())
                    }),
                }
                .into());
            }
        }
    }
    Ok(())
}

/// LOSAT's refusal of `check_failed_lookups` (exit status 1 with its message in every program; a
/// type of its own, so that a program that reports its engine's errors otherwise, BLASTX's
/// `native::report_error`, passes it on as a refusal), with the result set at which NCBI stops
/// (`set`, 0-based in the order given): a program that checks the result sets of several query
/// batches at once fails the batch of that set, after the reports of the batches before, as
/// NCBI's run does, and after what NCBI writes to the standard error of that batch's reports up
/// to the row (`report::query_warnings::posts_before_stop`: the subjects of the set's rows up to
/// it, `rows`, and its formatter's warning, `formatter_warning`).
#[derive(Debug)]
pub struct FailedLookupRefusal {
    message: String,
    set: usize,
    /// The subjects of the set's rows up to that row, in order (the row's last), as the result
    /// sets name them: the rows that NCBI's report of the set makes before it stops, whose reads
    /// of their subjects come before the stop.
    rows: Vec<usize>,
    formatter_warning: Option<String>,
}

impl FailedLookupRefusal {
    /// The result set (0-based) whose row NCBI has no result for.
    pub fn set(&self) -> usize {
        self.set
    }

    /// The warning of the formatter that NCBI makes for each query's report of the refused
    /// field list (`TabularSpec::formatter_warning`), written after each query's posts.
    pub fn formatter_warning(&self) -> Option<&str> {
        self.formatter_warning.as_deref()
    }

    /// What NCBI writes to the standard error of the refused batch's reports before it stops at
    /// the row (`report::query_warnings::posts_before_stop`), of `posts` and `reads`, the posts
    /// and subject reads of the batch's queries from its first up to the row's, whose reads name a
    /// subject as `read_subject` names a subject of the result sets (`rows`).
    pub fn posts_before_stop(
        &self,
        posts: &[crate::report::diag_log::Posts],
        reads: &[crate::report::query_warnings::QueryReads],
        read_subject: impl Fn(usize) -> u32,
    ) -> crate::report::diag_log::Posts {
        let rows: Vec<u32> = self.rows.iter().map(|&s_idx| read_subject(s_idx)).collect();
        crate::report::query_warnings::posts_before_stop(
            posts,
            reads,
            &rows,
            self.formatter_warning(),
        )
    }

    /// The first query of the batch whose report `error` stops (a `FailedLookupRefusal` whose
    /// result sets are the queries, in order; `starts`: the first query of each batch, in order):
    /// the reports of the queries before it are written before the refusal (TBLASTN, TBLASTX,
    /// which check every batch's rows at once). 0 for the first batch, one batch, or another
    /// error.
    pub fn batch_start(error: &anyhow::Error, starts: impl IntoIterator<Item = usize>) -> usize {
        error.downcast_ref::<Self>().map_or(0, |refusal| {
            starts
                .into_iter()
                .take_while(|&start| start <= refusal.set)
                .last()
                .unwrap_or(0)
        })
    }
}

impl std::fmt::Display for FailedLookupRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
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

/// What `check_failed_lookups` says of the database sequence `oid` whose lookup fails as `failure`
/// says: where its load throws `thrown` (a batch whose prefetch is off, or fails at a lookup that
/// throws on the id: `SubjectLookup::Throws`), that load; where no lookup fails and a definition
/// line has no Seq-id (`LookupFailure::lineless`), that line; otherwise its title
/// (`failed_database_title`).
pub fn failed_database_subject(
    oid: impl std::fmt::Display,
    thrown: Option<&crate::ncbi_exception::NcbiException>,
    failure: LookupFailure,
) -> String {
    match thrown {
        Some(exception) => format!(
            "database sequence {oid} has a Seq-align id whose lookup throws '{}' where NCBI BLAST+ 2.17.0 reads the reported sequences one by one (a batch whose prefetch PRE_FETCH_SEQS_LIMIT or the format turns off, or whose prefetch fails at a Seq-align id whose lookup throws on the id itself, with its warning Error pre-fetching sequence data)",
            String::from_utf8_lossy(exception.get_msg())
        ),
        None if failure.lineless && !failure.subject_id && !failure.all_ids => format!(
            "database sequence {oid} has a definition line without a Seq-id (a hand-edited database)"
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
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:557-565
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
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:647-655
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
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:862-878
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
            lineless: false,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:833-838, 862-875 (cited at
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
    /// Whether a definition line of the subject has no Seq-id
    /// (`db_subject::ParsedSubject::has_lineless_line`): no lookup fails, but `sallacc` prints the
    /// label of that line's best id, which is null (`failed_lookup_stop`).
    pub lineless: bool,
}

impl LookupFailure {
    /// The failure of a subject without parsed Seq-ids whose generated title throws: a `-subject`
    /// record's (`SetSubjectId` only), or a database sequence's (`database`: its one definition
    /// line holds the `gnl|BL_ORD_ID` id, so `x_SetSubjectIds` throws at it too).
    pub fn of_plain(database: bool) -> Self {
        Self {
            subject_id: true,
            all_ids: database,
            ..Self::default()
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
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:101-114
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
    // (cited at `TabularSpec::row_throw`: the subject block's `catch`, `return -1`)
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
    // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:279-297
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1100-1122
    /// The HSP's edit script (`gap_info`; the rows of `rows` are its columns), whose rows' aligned
    /// ranges NCBI's table reads (`serial::seqalign::row_end_gaps`): `length`, `gapopen` and
    /// `gaps` count the columns inside the query row's (`TabularHsp::within_query_range`, which
    /// the writer applies), and `qseq`, `sseq` and `btop` show a row's end gaps as NCBI's end
    /// character, a space (`push_alignment_row`), as in the Smith-Waterman alignments of
    /// `-use_sw_tback`. NCBI decides by where a `-` comes from: a gap segment outside the row's
    /// aligned segments is the end character; a `-` that is a residue (a PSSM query's gap residue,
    /// a database's stored gap character) or fills a segment that has a start is the row's
    /// letter. Empty for an ungapped HSP and for the rows whose `-` at an end fills a segment past
    /// the query's end: BLASTN's (whose X-drop and greedy alignments start and end with a
    /// substitution: `algorithm::blastn::tabular`, remaining-gaps X1b-p) and `-subject`'s mapped
    /// row of BLASTX (`blastx::empty_rows::MappedAlignment`).
    pub ops: &'a [GapEditOp],
    /// `qcovhsp` (the integer of NCBI's `hsp_percent_coverage`).
    pub qcovhsp: i32,
}

impl TabularHsp<'_> {
    /// The row with `length`, `gaps` and `gap_opens` as NCBI's tables count them
    /// (`serial::seqalign::AlignLengths::within_query_range` of `ops`, from the whole alignment's;
    /// `pident` and `ppos` then divide by the shorter length, `mismatch` stays). The writer reads
    /// every row so (`TabularWriter::row`).
    pub(crate) fn within_query_range(self) -> Self {
        let crate::report::serial::seqalign::AlignLengths {
            length,
            gaps,
            gap_opens,
        } = crate::report::serial::seqalign::AlignLengths {
            length: self.length,
            gaps: self.gaps,
            gap_opens: self.gap_opens,
        }
        .within_query_range(self.ops);
        Self {
            length,
            gaps,
            gap_opens,
            ..self
        }
    }
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

// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:548
/// The exception that ends a tabular report at its first row when the query's Bioseq has the
/// molecule `other` (psiblast's `-in_msa` with `[molecule=other]`, or `-in_pssm` of such a PSSM):
/// the object manager reads the query's letters of the coding number `coding` (iupacaa 2, ncbieaa 8,
/// ncbistdaa 10) as the coding of a nucleotide sequence
/// (`Incompatible sequence codings: <coding> -> 1`, `CSeqVector_CI::x_FillCache`), where the rows' letters
/// are read ([`TabularSpec::reads_sequences`]). The report's text before the row stays; the run
/// ends with NCBI's exception (exit status 255).
pub fn query_sequence_exception(coding: i32) -> crate::ncbi_exception::NcbiException {
    crate::ncbi_exception::NcbiException::new(
        crate::ncbi_exception::ThrowSite {
            file: "src/objmgr/seq_vector_ci.cpp",
            line: 548,
            kind: "CSeqVectorException::eCodingError",
            function: Some("ncbi::objects::CSeqVector_CI::x_FillCache"),
        },
        format!("Incompatible sequence codings: {coding} -> 1"),
    )
}

/// The end of a report at a row that reads the letters of a query whose molecule is `other`
/// ([`query_sequence_exception`]): what the report holds so far is written (`buffered`), the
/// exception's unwinding flushes the stream (NCBI's flush, as at `run_end`), and the exception is
/// the run's error.
pub fn query_sequence_end(
    writer: &mut dyn crate::report::ostream::ReportWrite,
    buffered: &[u8],
    coding: i32,
) -> anyhow::Error {
    if let Err(error) = writer
        .write_all(buffered)
        .and_then(|()| writer.ncbi_flush())
    {
        return error.into();
    }
    crate::cli::NativeError::exception(&query_sequence_exception(coding)).into()
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203 (the `CATCH_ALL` of
// every application, cited at `cli::NativeError::exception`)
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:160-163 (cited at
// `blastn::tabular::write_blastn_tabular`)
/// The end of the run at a row whose subject's load throws outside the subject block
/// (`TabularWriter::ends_run`): writes the report's rows before it, `buffered`, to `writer`, and
/// the exception's unwinding runs `~CBlastTabularInfo`, which flushes the stream (NCBI's flush,
/// the one of the query's end: a device that fails there ends the run as any tabular failure,
/// not at the exception), ends the reads (`end_reads`: `QueryWarnings::end_run`,
/// `ReadCursor::end_run`; NCBI loads nothing after it) and gives NCBI's exception as the run's
/// error (`blast_xml::writer::run_ending_exception`: `BLAST Database error: <message>`, exit
/// status 2).
pub fn run_end(
    writer: &mut dyn crate::report::ostream::ReportWrite,
    buffered: &[u8],
    exception: &crate::ncbi_exception::NcbiException,
    end_reads: impl FnOnce(),
) -> anyhow::Error {
    if let Err(error) = writer
        .write_all(buffered)
        .and_then(|()| writer.ncbi_flush())
    {
        return error.into();
    }
    end_reads();
    crate::report::blast_xml::writer::run_ending_exception(exception)
}

/// The writer of a tabular report (`TabularWriter::new`).
pub struct TabularWriter<'a> {
    spec: &'a TabularSpec,
    /// The `# BLASTN 2.17.0+` line's program and version.
    program_version: &'a str,
    /// The matrix that scores the rows' pairs of letters for the positives, where the report
    /// prints them (`new`).
    matrix: Option<ScoringMatrix>,
}

impl<'a> TabularWriter<'a> {
    /// The writer of the report `spec`, whose `# BLASTX 2.17.0+` line names `program_version`;
    /// `matrix` is the matrix of a protein search's positives (`None` for BLASTN and for a
    /// formatter without one: positives are identities): the formatter's `m_ScoringMatrix`,
    /// which `CBlastFormat` makes with `GetAsciiProteinMatrix` (outfmt 0's matrix,
    /// `utils::matrix::protein_display_score`), not the XML's `CBlastFormattingMatrix`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:183
    /// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:828-830
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1017-1020
    /// NCBI reference (598d8ae6): c++/include/objtools/align_format/tabular.hpp:444-449
    /// The count of the positives is read only by `positive` and `ppos` (`x_PrintNumPositives`,
    /// `x_PrintPercentPositives`): the writer of a report without them keeps no matrix, so its
    /// rows score no pair of letters (the same bytes; decided once for the writer, not per row).
    pub fn new(
        spec: &'a TabularSpec,
        program_version: &'a str,
        matrix: Option<ScoringMatrix>,
    ) -> Self {
        Self {
            spec,
            program_version,
            matrix: matrix.filter(|_| spec.has_any(&[Positives, PercentPositives])),
        }
    }

    /// outfmt 20's line of the field names (`PrintProlog` → `PrintFieldSpecs`).
    pub fn prolog(&self, out: &mut dyn Write) -> io::Result<()> {
        if self.spec.kind != TabularKind::CsvWithHeader {
            return Ok(());
        }
        let mut line = Vec::new();
        for (i, field) in self.spec.fields.iter().enumerate() {
            if i > 0 {
                line.extend_from_slice(&self.spec.delimiter);
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
        // The counts inside the query row's aligned range (`TabularHsp::ops`).
        let hsp = &hsp.within_query_range();
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
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:169-194
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
            // before the one that throws; `SetSubjectId` of `list` in the catch pushes the Seq-align id
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
                out.extend_from_slice(&spec.delimiter);
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
        // (the identities and positives counted again from the rows, and the BTOP string). A pair
        // of letters is scored only by the matrix of a writer whose report prints the positives
        // (`new`).
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
                out.extend_from_slice(&spec.delimiter);
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
                QuerySeq => push_alignment_row(
                    out,
                    hsp.rows.map_or(&[][..], |rows| rows.0),
                    hsp.ops,
                    AlignRow::Query,
                ),
                SubjectSeq => push_alignment_row(
                    out,
                    hsp.rows.map_or(&[][..], |rows| rows.1),
                    hsp.ops,
                    AlignRow::Subject,
                ),
                Evalue => out.extend_from_slice(format_evalue_ncbi_tabular(hsp.evalue).as_bytes()),
                BitScore => out.extend_from_slice(format_bitscore_ncbi(hsp.bit_score).as_bytes()),
                Score => push_num(out, hsp.score),
                AlignmentLength => push_num(out, hsp.length),
                PercentIdentical => out.extend_from_slice(
                    format_percent_identity_ncbi(num_ident, hsp.length, 3).as_bytes(),
                ),
                NumIdentical => push_num(out, num_ident),
                // NCBI reference (598d8ae6): c++/include/objtools/align_format/tabular.hpp:478
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
                        write_btop(out, qrow, srow, hsp.ops);
                    }
                }
                SubjectTaxId | SubjectTaxIds | SubjectSciName | SubjectCommonName
                | SubjectBlastName | SubjectSuperKingdom | SubjectSciNames | SubjectCommonNames
                | SubjectBlastNames | SubjectSuperKingdoms | SubjectTitle | SubjectAllTitles => {
                    self.subject_field(out, *field, subject, &names)
                }
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1063-1070
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

/// The BTOP string of two rows, the columns of the edit script `ops` (`TabularHsp::ops`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1001-1027
/// The rows' end gaps are spaces (`push_alignment_row`), so an end-gap column is the pair of a
/// letter and a space. Rows without end gaps (every alignment but some Smith-Waterman redos', and
/// the rows of an empty script) take the columns as they are.
pub fn write_btop(out: &mut Vec<u8>, qrow: &[u8], srow: &[u8], ops: &[GapEditOp]) {
    let columns = qrow.iter().copied().zip(srow.iter().copied());
    let q_span = row_end_gaps(ops, AlignRow::Query).span(qrow.len());
    let s_span = row_end_gaps(ops, AlignRow::Subject).span(srow.len());
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
/// (`serial::seqalign::RowEndGaps::span`): NCBI's end character, a space, outside them.
fn end_char(letter: u8, column: usize, span: (usize, usize)) -> u8 {
    if column < span.0 || column >= span.1 {
        b' '
    } else {
        letter
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:982-986
/// `qseq` or `sseq`: the row `row` (`which`) of the edit script `ops` with `-` for its inner gaps
/// and a space (NCBI's end character, which the tabular writer does not set) for its end gaps
/// (`serial::seqalign::row_end_gaps`, `TabularHsp::ops`), where it has some.
fn push_alignment_row(out: &mut Vec<u8>, row: &[u8], ops: &[GapEditOp], which: AlignRow) {
    let span = row_end_gaps(ops, which).span(row.len());
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:629-643
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

/// The query coverage of the HSPs of one subject in percent cut to an integer (`qcovhsp`), the
/// one rule of every program: `x_SetQueryCovSeqalign` reads the `hsp_percent_coverage` of the
/// Seq-align that the formatter gets, which `Blast_HSPGetQueryCoverage` makes of the HSP's span,
/// its letters in its query context, and the context's length, the length that the program's
/// Seq-align builder gives it (the searched letters; a translated query's span is residues of its
/// frame). Where that Seq-align has no such score, it counts the Seq-align's own span of the
/// query (letters of the record: a translated query's nucleotides) against the tabular query
/// length (`m_QueryLength`: the `-query_loc` range as given, tabular.cpp:786-790).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1397-1447
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1308-1309,1382-1383
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3162-3205
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:518-533
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:737-747
/// A gapped search makes one Seq-align per HSP, with its scores. An ungapped search makes one
/// Seq-align of a subject's HSPs, each a Dense-diag (BLASTN, BLASTP) or a Std-seg (the
/// translated programs) with its scores. The formatter splits a Seq-align of several into one per
/// HSP that takes those scores, and keeps one of a single HSP as it is: its scores are on the
/// Dense-diag or Std-seg, and `x_GetNamedScore` reads only the Seq-align's own, so the span
/// counts against the tabular query length. In BLASTN, BLASTP and TBLASTN both spans are the
/// HSP's query letters, so the counts differ only where the range as given passes the record's
/// end.
///
/// Everything the choice reads belongs to the subject (`of_subject`), so a report makes it once
/// per subject and its rows only divide (`of_hsp`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct QueryCoveragePerHsp {
    /// The Seq-align's own span counts (a subject's single HSP of an ungapped search), not the
    /// HSP's span in its context.
    seqalign: bool,
    /// The length that the counted span divides.
    length: i32,
}

impl QueryCoveragePerHsp {
    /// The rule of a subject's `hsps` HSPs of a search that is `ungapped` or not: their spans
    /// in a query context of `context_length`, or the Seq-align's span against `query_length`.
    pub fn of_subject(ungapped: bool, hsps: usize, context_length: i32, query_length: i32) -> Self {
        if ungapped && hsps == 1 {
            Self {
                seqalign: true,
                length: query_length,
            }
        } else {
            Self {
                seqalign: false,
                length: context_length,
            }
        }
    }

    /// The `qcovhsp` of one of the subject's HSPs: `span` letters in its query context,
    /// `seqalign_span` in its Seq-align.
    #[inline]
    pub fn of_hsp(self, span: i32, seqalign_span: i32) -> i32 {
        let covered = if self.seqalign { seqalign_span } else { span };
        crate::report::serial::seqalign::query_coverage(0, covered, self.length) as i32
    }
}

#[cfg(test)]
#[path = "tabular_tests.rs"]
mod tests;
