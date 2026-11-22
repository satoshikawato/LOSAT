//! The organism report (outfmt 18, NCBI's "Tax BLAST report", `CTaxFormat` in text mode without
//! the taxonomy server): for each query a caption, the outfmt 0 `Query=` block, and a table of
//! the query's subjects grouped by the taxid of each definition line, under the names that the
//! taxonomy database gives for it. One writer for every program; a program's adapter hands over,
//! for each query, its title and length and its HSPs in report order (the final hit list as the
//! tabular reports read it, `report::tabular`), each with its subject, and the definition lines
//! of the subjects as the report shows them (`TaxDefline`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1461-1467
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:873-902
//! LOSAT's HSPs are NCBI's Seq-aligns one by one (`PrepareBlastUngappedSeqalign` splits an ungapped
//! Seq-align of several HSPs into one per HSP, as a gapped search has them).
//!
//! The column widths are NCBI's `unsigned int`s: the description's width wraps when the other
//! columns take more than 96 bytes, and NCBI then pads with about 4 GiB of spaces. The report is
//! therefore built of segments (`Text`: bytes, or a run of spaces), on which the template calls of
//! NCBI (`MapTemplate`, `NStr::Replace`) work as on NCBI's strings (a token holds no space, so it
//! never reaches into a run), and written by streaming the runs.

use std::io::{self, Write};

use anyhow::Result;

use crate::report::outfmt6::{format_bitscore_ncbi, format_evalue_ncbi};
use crate::report::query_warnings::ReadCursor;
use crate::report::tabular::TabularNames;
use crate::seqdb::taxdb::TaxNames;

/// The greater of `kMinTaxFormatLineLength` and `kFormatLineLength`: the report's line length (`m_LineLength`).
///
/// NCBI reference (598d8ae6): c++/include/objtools/align_format/taxFormat.hpp:58
const LINE_LENGTH: u32 = 100;

/// The longest run of spaces kept as bytes; a longer one is a `Seg::Spaces` (only the widths
/// that wrap make one).
const INLINE_SPACES: u64 = 4096;

/// One segment of report text.
#[derive(Clone, Debug, PartialEq, Eq)]
enum Seg {
    Bytes(Vec<u8>),
    Spaces(u64),
}

/// Report text: bytes, with the runs of more than `INLINE_SPACES` spaces kept as counts. Two
/// `Seg::Bytes` never follow each other, so a token (which holds no space) lies within one.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct Text(Vec<Seg>);

impl Text {
    fn of(bytes: &[u8]) -> Self {
        let mut text = Text::default();
        text.push_bytes(bytes);
        text
    }

    fn push_bytes(&mut self, bytes: &[u8]) {
        if bytes.is_empty() {
            return;
        }
        match self.0.last_mut() {
            Some(Seg::Bytes(last)) => last.extend_from_slice(bytes),
            _ => self.0.push(Seg::Bytes(bytes.to_vec())),
        }
    }

    fn push_spaces(&mut self, count: u64) {
        if count == 0 {
            return;
        }
        if count <= INLINE_SPACES {
            match self.0.last_mut() {
                Some(Seg::Bytes(last)) => last.resize(last.len() + count as usize, b' '),
                _ => self.0.push(Seg::Bytes(vec![b' '; count as usize])),
            }
            return;
        }
        match self.0.last_mut() {
            Some(Seg::Spaces(last)) => *last += count,
            _ => self.0.push(Seg::Spaces(count)),
        }
    }

    fn push_text(&mut self, text: &Text) {
        for seg in &text.0 {
            match seg {
                Seg::Bytes(bytes) => self.push_bytes(bytes),
                Seg::Spaces(count) => self.push_spaces(*count),
            }
        }
    }

    /// Whether the text holds `token` (which holds no space).
    fn contains(&self, token: &[u8]) -> bool {
        self.0.iter().any(|seg| match seg {
            Seg::Bytes(bytes) => find(bytes, token, 0).is_some(),
            Seg::Spaces(_) => false,
        })
    }

    /// `CAlignFormatUtil::MapTemplate` of `text`, `name` and `value`: every `<@name@>` replaced by `value`,
    /// left to right; the text put in is not searched again by this call.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3593-3599
    fn map(&self, name: &str, value: &Text) -> Text {
        let token = [b"<@", name.as_bytes(), b"@>"].concat();
        let mut out = Text::default();
        for seg in &self.0 {
            match seg {
                Seg::Spaces(count) => out.push_spaces(*count),
                Seg::Bytes(bytes) => {
                    let mut start = 0;
                    while let Some(at) = find(bytes, &token, start) {
                        out.push_bytes(&bytes[start..at]);
                        out.push_text(value);
                        start = at + token.len();
                    }
                    out.push_bytes(&bytes[start..]);
                }
            }
        }
        out
    }

    fn map_bytes(&self, name: &str, value: &[u8]) -> Text {
        self.map(name, &Text::of(value))
    }

    /// The bytes of a text without runs (a short one: the organism line before its centring).
    fn flat(&self) -> Vec<u8> {
        let mut out = Vec::new();
        for seg in &self.0 {
            match seg {
                Seg::Bytes(bytes) => out.extend_from_slice(bytes),
                Seg::Spaces(count) => out.resize(out.len() + *count as usize, b' '),
            }
        }
        out
    }

    /// The number of bytes of the text.
    fn len(&self) -> u64 {
        self.0
            .iter()
            .map(|seg| match seg {
                Seg::Bytes(bytes) => bytes.len() as u64,
                Seg::Spaces(count) => *count,
            })
            .sum()
    }

    /// The insertion of `text` into `out`: one insertion of NCBI's, streamed (`cli::insert_streamed`), so that a run
    /// of spaces is never held whole.
    fn insert_into(&self, out: &mut crate::cli::NcbiOstream<'_>) -> io::Result<()> {
        crate::cli::insert_streamed(out, self.len(), |insertion| {
            for seg in &self.0 {
                match seg {
                    Seg::Bytes(bytes) => insertion.bytes(bytes)?,
                    Seg::Spaces(count) => insertion.spaces(*count)?,
                }
            }
            Ok(())
        })
    }
}

/// The first `needle` in `haystack` at or after `from`.
fn find(haystack: &[u8], needle: &[u8], from: usize) -> Option<usize> {
    if needle.is_empty() || haystack.len() < needle.len() {
        return None;
    }
    (from..=haystack.len() - needle.len()).find(|&at| haystack[at..at + needle.len()] == *needle)
}

/// Where `AddSpaces` puts the spaces (`CAlignFormatUtil::ESpacesFormat`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Place {
    /// `eSpacePosAtLineEnd` (`MapSpaceTemplate`'s default).
    End,
    /// `eSpacePosToCenter`.
    Center,
}

/// `CAlignFormatUtil::AddSpaces` of `value`, `width` and `place` with `eAddEOLAtLineStart` (`eol_start`)
/// and `eAddEOLAtLineEnd` (`eol_end`). `width` is NCBI's `size_t`: a width of 0 to 2 cuts
/// nothing (`width - 3` wraps) and still adds `...`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3610-3640
fn add_spaces(value: &[u8], width: u64, place: Place, eol_start: bool, eol_end: bool) -> Text {
    let (value, spaces): (std::borrow::Cow<'_, [u8]>, u64) = if width >= value.len() as u64 {
        let spaces = width - value.len() as u64 + 1;
        (
            value.into(),
            if place == Place::Center {
                spaces / 2
            } else {
                spaces
            },
        )
    } else {
        // `substr(0, width - 3)`: `width - 3` wraps below 3 (`substr` then keeps the whole value).
        let cut = width.wrapping_sub(3).min(value.len() as u64) as usize;
        ([&value[..cut], b"..."].concat().into(), 1)
    };
    let mut text = Text::default();
    if eol_start {
        text.push_bytes(b"\n");
    }
    if place == Place::Center {
        text.push_spaces(spaces);
    }
    text.push_bytes(&value);
    text.push_spaces(spaces);
    if eol_end {
        text.push_bytes(b"\n");
    }
    text
}

/// The title of a definition line as the report makes it (`x_FillTaxDispParams`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TaxTitle {
    /// The stored title, or the generated one where the line has none.
    Text(Vec<u8>),
    /// `GenerateDefline` throws: `NStr::HtmlDecode` (create_defline.cpp:4066) cannot guess the
    /// title's encoding (ncbistr.cpp:4528-4532). The report's `catch` (taxFormat.cpp:1033-1035)
    /// then drops the rest of the HSP's subject before this line widens the columns
    /// (`x_InitTextFormatInfo` comes after the title, taxFormat.cpp:926-927), and the subject is
    /// not the previous id: its next HSP is read again.
    Throws,
    /// A title that LOSAT does not write as NCBI: refused where NCBI makes it (a line that the
    /// report reaches), or with `shown_only` where NCBI shows it (a line whose taxid the taxonomy
    /// database names; NCBI drops the row of any other line before it shows the title).
    Refused { reason: String, shown_only: bool },
}

/// One definition line of a subject as the report shows it (`CTaxFormat::SSeqInfo`): the label
/// (`CAlignFormatUtil::GetLabel` of its `WorstRank`-best Seq-id, without version), the title (the
/// generated title of the subject's Bioseq where it has none) and the taxid (0 for none).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TaxDefline {
    pub label: Vec<u8>,
    pub title: TaxTitle,
    pub taxid: i32,
}

impl TaxDefline {
    /// The definition line of a `-subject` record read without `-parse_deflines`: its Bioseq has
    /// the local id `Subject_N` (`local_id`) and no `Blast-def-line`, so its title is the
    /// generated one of the record's title `title`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:911-929
    pub fn of_record(local_id: &[u8], title: &[u8], protein: bool) -> Self {
        TaxDefline {
            label: local_id.to_vec(),
            title: generated_title(
                title,
                protein,
                &crate::report::defline::BioseqFacts::default(),
            ),
            taxid: 0,
        }
    }

    /// The definition line of the database sequence `oid` without parsed Seq-ids: the label
    /// of its `gnl|BL_ORD_ID|<oid>` id is the OID, its title the stored one (`title`)
    /// or, where that is empty, the Bioseq's generated title, which NCBI makes from the name of
    /// its source (`source`, the taxonomy database's name of the taxid:
    /// `report::defline::BioseqFacts::of_source`; refused where the report shows it for a source
    /// that LOSAT refuses).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:961-971
    pub fn of_stored(
        oid: crate::seqdb::Oid,
        title: &[u8],
        source: &crate::report::defline::SourceTaxname,
        taxid: i32,
        protein: bool,
    ) -> Self {
        let title = if !title.is_empty() {
            TaxTitle::Text(title.to_vec())
        } else {
            match source.name() {
                Ok(name) => generated_title(
                    title,
                    protein,
                    &crate::report::defline::BioseqFacts::of_source(name),
                ),
                Err(reason) => TaxTitle::Refused {
                    reason: format!("has {reason}"),
                    shown_only: true,
                },
            }
        };
        TaxDefline {
            label: oid.get().to_string().into_bytes(),
            title,
            taxid,
        }
    }
}

/// `CDeflineGenerator().GenerateDefline` (flags 0) of a Bioseq with the title `title` and the
/// Seq-ids `ids`: the text; `Throws` where NCBI's `HtmlDecode` cannot guess the encoding (the
/// report skips the subject there); refused where NCBI's `x_CleanAndCompress` reads past the end
/// of a protein title of punctuation (NCBI crashes; approved exception 2 of
/// PD-LOSAT-NCBI-DEFECTS has a nucleotide title stop at the end of the string).
pub(crate) fn generated_title(
    title: &[u8],
    protein: bool,
    ids: &crate::report::defline::BioseqFacts,
) -> TaxTitle {
    match crate::report::defline::generate_bioseq_defline(title, protein, false, ids) {
        Ok(defline) if protein && defline.reads_past_end => TaxTitle::Refused {
            reason: concat!("has a title ", crate::report::defline::punctuation_title!()).into(),
            shown_only: false,
        },
        Ok(defline) => TaxTitle::Text(defline.text),
        Err(_) => TaxTitle::Throws,
    }
}

/// One HSP of a query's final hit list, in report order: its subject (an index that the adapter's
/// `TaxSubjects` reads) and its scores.
#[derive(Clone, Copy, Debug)]
pub struct TaxHsp {
    pub subject: usize,
    pub evalue: f64,
    pub bit_score: f64,
}

/// The subjects of a report, as a program's adapter gives them (`SubjectTable`).
pub trait TaxSubjects {
    /// The definition lines of `subject` that the report shows (those a positive list keeps,
    /// `CAlignFormatUtil::MatchSeqInSeqList`; one for a subject without `Blast-def-line`s).
    fn deflines(&self, subject: usize) -> &[TaxDefline];

    /// Whether the Seq-align id of `subject` matches that of `previous` (`CSeq_id::Match`): the
    /// report skips such an HSP after one of `previous` that it took.
    fn same(&self, previous: usize, subject: usize) -> bool {
        previous == subject
    }

    /// LOSAT's refusal of a title of `subject` (`TaxTitle::Refused`), which names it.
    fn refused(&self, subject: usize, reason: &str) -> anyhow::Error;

    /// Whether the load of `subject` throws where its query's batch has its prefetch off
    /// (`db_subject::HitLoads::thrown`), or at every load where a file of the taxonomy database
    /// does not open and the load looks a positive taxid up
    /// (`SubjectTable::with_loads_that_throw`): the report leaves its HSPs out (`organism_report`).
    fn throws(&self, _subject: usize) -> bool {
        false
    }
}

/// LOSAT's refusal of an organism report that it does not write as NCBI (a title, a taxonomy
/// name, a `./.ncbirc`), made before the first byte of the report's query batch (each batch's
/// report is made before it is written), or, for a query that NCBI does not search, after the
/// bytes that NCBI writes before it fails (`TaxReport::write_query`): no error of NCBI's engine.
#[derive(Debug)]
pub struct TaxReportRefusal(String);

impl std::fmt::Display for TaxReportRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for TaxReportRefusal {}

/// A subject without parsed Seq-ids as the adapter has it (`TaxSource`).
#[derive(Clone, Copy, Debug)]
pub enum TaxPlain<'a> {
    /// A `-subject` record: its local id `Subject_N` and title.
    Record { local_id: &'a [u8], title: &'a [u8] },
    /// The database sequence `oid`: its stored title, the name of its source (from which
    /// NCBI makes an empty title, `algorithm::common::database::source_name`) and its taxid.
    Stored {
        oid: crate::seqdb::Oid,
        title: &'a [u8],
        source: &'a crate::report::defline::SourceTaxname,
        taxid: i32,
    },
}

/// Where the definition lines of a subject come from: its parsed Seq-ids where it has them (a
/// database sequence's definition lines, or a `-subject` record read under `-parse_deflines`),
/// else the record or stored sequence `plain`, which also names the subject in LOSAT's refusals.
#[derive(Clone, Copy, Debug)]
pub struct TaxSource<'a> {
    pub parsed: Option<&'a crate::report::db_subject::ParsedSubject>,
    pub plain: TaxPlain<'a>,
}

impl TaxSource<'_> {
    /// The definition lines of the subject as the report shows them (`TaxDefline`), of a protein
    /// (`protein`) or nucleotide subject.
    fn deflines(&self, protein: bool) -> Vec<TaxDefline> {
        match (self.parsed, self.plain) {
            (Some(subject), _) => subject.tax_deflines(protein),
            (None, TaxPlain::Record { local_id, title }) => {
                vec![TaxDefline::of_record(local_id, title, protein)]
            }
            (
                None,
                TaxPlain::Stored {
                    oid,
                    title,
                    source,
                    taxid,
                },
            ) => vec![TaxDefline::of_stored(oid, title, source, taxid, protein)],
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1130-1155
    /// Whether NCBI's load of the subject Bioseq looks a positive taxid up
    /// (`CSeqDBVol::x_GetTaxonomy`, at each definition line with one): a database sequence with
    /// such a line (`ParsedSubject::first_source_taxid`; the stored taxid without parsed Seq-ids).
    /// A `-subject` record is not loaded from a database.
    fn looks_taxid_up(&self) -> bool {
        match (self.parsed, self.plain) {
            (Some(subject), _) => subject.first_source_taxid().is_some(),
            (None, TaxPlain::Stored { taxid, .. }) => taxid > 0,
            (None, TaxPlain::Record { .. }) => false,
        }
    }

    /// How LOSAT's refusals name the subject.
    fn name(&self) -> SubjectName {
        match (self.parsed, self.plain) {
            (_, TaxPlain::Stored { oid, .. }) => SubjectName::Database(oid),
            (Some(subject), TaxPlain::Record { .. }) => SubjectName::Record(
                String::from_utf8_lossy(&subject.alignment_fasta()).into_owned(),
            ),
            (None, TaxPlain::Record { title, .. }) => {
                SubjectName::Record(String::from_utf8_lossy(title).into_owned())
            }
        }
    }
}

/// How LOSAT's refusals name a subject (`TaxSource::name`).
#[derive(Debug)]
enum SubjectName {
    Database(crate::seqdb::Oid),
    Record(String),
}

/// The subjects of a report by subject index (a subject set's index, a database's rank): their
/// definition lines and names, and the parsed subjects, whose Seq-align ids `Match` where they
/// hold equal ids (`db_subject::joined_in_reports`).
pub struct SubjectTable<'p> {
    subjects: std::collections::HashMap<usize, (Vec<TaxDefline>, SubjectName)>,
    /// The subjects whose load looks a positive taxid up (`TaxSource::looks_taxid_up`).
    taxid_loads: std::collections::HashSet<usize>,
    /// Whether every load of a database sequence that looks a positive taxid up throws
    /// (`with_loads_that_throw`).
    taxid_loads_throw: bool,
    parsed: &'p [Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
    /// What the reports load of each subject (`db_subject::HitLoads`; `None`: every subject).
    loads: Option<&'p crate::report::db_subject::HitLoads>,
    program: String,
}

impl<'p> SubjectTable<'p> {
    /// The table of the subjects of `queries`' HSPs: `parsed` are the parsed subjects by subject
    /// index (empty when none has parsed Seq-ids), `plain(index)` a subject as the adapter has it
    /// (`TaxPlain`); the subjects are proteins where `protein`; `loads` what the reports load of
    /// each subject (`None`: every subject). `program` names the program in LOSAT's refusals.
    pub fn new<'a>(
        queries: &[TaxQuery<'_>],
        parsed: &'p [Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
        plain: impl Fn(usize) -> TaxPlain<'a>,
        protein: bool,
        loads: Option<&'p crate::report::db_subject::HitLoads>,
        program: &str,
    ) -> Self {
        let mut subjects = std::collections::HashMap::new();
        let mut taxid_loads = std::collections::HashSet::new();
        for hsp in queries.iter().flat_map(|query| query.hsps.iter().flatten()) {
            subjects.entry(hsp.subject).or_insert_with(|| {
                let source = TaxSource {
                    parsed: parsed.get(hsp.subject).and_then(Option::as_deref),
                    plain: plain(hsp.subject),
                };
                if source.looks_taxid_up() {
                    taxid_loads.insert(hsp.subject);
                }
                (source.deflines(protein), source.name())
            });
        }
        SubjectTable {
            subjects,
            taxid_loads,
            taxid_loads_throw: false,
            parsed,
            loads,
            program: program.to_owned(),
        }
    }

    /// The table where every load of a database sequence that looks a positive taxid up throws
    /// (`throws`: a file of the run's taxonomy database does not open,
    /// `common::database::TaxdbRun::load_throws`, `TaxNamesDb::loads_throw`).
    pub fn with_loads_that_throw(self, throws: bool) -> Self {
        SubjectTable {
            taxid_loads_throw: throws,
            ..self
        }
    }

    /// The parsed subject of `subject`, if it has one.
    fn parsed(&self, subject: usize) -> Option<&'p crate::report::db_subject::ParsedSubject> {
        self.parsed.get(subject).and_then(Option::as_deref)
    }
}

impl TaxSubjects for SubjectTable<'_> {
    fn deflines(&self, subject: usize) -> &[TaxDefline] {
        self.subjects
            .get(&subject)
            .map_or(&[], |(lines, _)| lines.as_slice())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:989-993
    fn same(&self, previous: usize, subject: usize) -> bool {
        previous == subject
            || matches!(
                (self.parsed(previous), self.parsed(subject)),
                (Some(a), Some(b)) if crate::report::db_subject::joined_in_reports(a, b)
            )
    }

    fn refused(&self, subject: usize, reason: &str) -> anyhow::Error {
        let program = &self.program;
        anyhow::Error::new(TaxReportRefusal(match self.subjects.get(&subject) {
            Some((_, SubjectName::Database(oid))) => format!(
                "database sequence {} {reason}; this is not supported by LOSAT's {program}",
                oid.get()
            ),
            Some((_, SubjectName::Record(name))) => format!(
                "the subject record {name} {reason}, which is not supported by LOSAT's {program}"
            ),
            None => format!("a subject {reason}; this is not supported by LOSAT's {program}"),
        }))
    }

    fn throws(&self, subject: usize) -> bool {
        (self.taxid_loads_throw && self.taxid_loads.contains(&subject))
            || self
                .loads
                .is_some_and(|loads| loads.thrown(subject).is_some())
    }
}

/// The taxonomy database that a run's organism reports look names up in (`CSeqDB::GetTaxInfo`,
/// the static `CTaxDBFileInfo` of `GetTaxNames`, opened at the first lookup: one per run in
/// LOSAT, a `-db` run's `TaxdbRun` and a `-subject` run's `SubjectTaxdb`), handed to the report as
/// `TabularNames::for_report`. The report asks it for the taxids that NCBI's report looks up,
/// in NCBI's order, and nothing else (`TaxLookups`).
#[derive(Clone)]
pub struct TaxNamesDb {
    taxdb: std::sync::Arc<crate::seqdb::taxdb::LazyTaxDb>,
    /// The texts of an open of the database by the rows of an earlier query batch
    /// (`TaxdbRun::open_for_rows`), which the run's first lookup of a report writes.
    opened_before: Option<std::sync::Arc<std::sync::Mutex<Vec<String>>>>,
    /// LOSAT's error of an open or a lookup that the reader refuses or cannot read, naming the
    /// program (`common::database::taxdb_error`).
    error: fn(crate::seqdb::SeqDbError, &str) -> anyhow::Error,
    /// Whether every load of a database sequence that looks a positive taxid up throws
    /// (`where_loads_throw`).
    loads_throw: bool,
}

impl TaxNamesDb {
    /// The lookups through `taxdb`; `opened_before` and `error` as the fields say.
    pub(crate) fn new(
        taxdb: std::sync::Arc<crate::seqdb::taxdb::LazyTaxDb>,
        opened_before: Option<std::sync::Arc<std::sync::Mutex<Vec<String>>>>,
        error: fn(crate::seqdb::SeqDbError, &str) -> anyhow::Error,
    ) -> Self {
        TaxNamesDb {
            taxdb,
            opened_before,
            error,
            loads_throw: false,
        }
    }

    /// The lookups of a `-db` run where every load of a database sequence that looks a positive
    /// taxid up throws (`throws`: a file of the database does not open,
    /// `common::database::TaxdbRun::load_throws`): the report leaves such a subject's HSPs out at
    /// its load (`SubjectTable::with_loads_that_throw`).
    pub(crate) fn where_loads_throw(self, throws: bool) -> Self {
        TaxNamesDb {
            loads_throw: throws,
            ..self
        }
    }
}

impl std::fmt::Debug for TaxNamesDb {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TaxNamesDb")
            .field("taxdb", &self.taxdb)
            .finish_non_exhaustive()
    }
}

/// The lookups of one query batch's reports (`CSeqDB::GetTaxInfo`), made as NCBI's report makes
/// them (`organism_report`): what each found, the text each writes to standard error each time
/// it is called, the texts of opening the database, which the batch's first lookup writes where
/// the run opened it then, and what the batch's reports read of `./.ncbirc`.
pub struct TaxLookups {
    /// The run's database; `None` without `taxdb.bti` on the search path: every lookup fails,
    /// silently.
    db: Option<TaxNamesDb>,
    program: String,
    /// The program's name in lower case, as NCBI's messages show it.
    lower: String,
    /// What each lookup of a taxid found (`None`: it threw) and the texts it writes.
    answers: std::collections::HashMap<i32, (Option<TaxNames>, Vec<String>)>,
    /// The texts of opening that the batch's first lookup writes.
    open: Vec<String>,
    /// Whether the batch has made a lookup.
    looked_up: bool,
    /// `<@taxBrowserURL@>` (`tax_browser_url`).
    browser_url: Vec<u8>,
    /// What one read of `./.ncbirc` by the report posts (`read_ncbirc`); the report of each query
    /// reads it twice.
    registry_posts: Vec<u8>,
}

impl TaxLookups {
    /// The lookups of a batch through the run's database that `names` carries
    /// (`TabularNames::for_report`; `None`: no `taxdb.bti` on the search path), with what the
    /// reports read of `./.ncbirc`. LOSAT refuses a `./.ncbirc` that NCBI's report cannot read
    /// (it reads the file for every query) before the batch's first byte.
    pub fn new(names: Option<&TabularNames>, program: &str) -> Result<Self> {
        Ok(Self::with_registry(
            names,
            program,
            read_ncbirc(program)?.as_ref(),
        ))
    }

    /// `new` with the read of `./.ncbirc` `registry` (`None`: no file).
    fn with_registry(
        names: Option<&TabularNames>,
        program: &str,
        registry: Option<&crate::blastinput::ncbi_registry::RegistryRead>,
    ) -> Self {
        let lower = program.to_ascii_lowercase();
        TaxLookups {
            db: names.and_then(TabularNames::report_db).cloned(),
            program: program.to_owned(),
            registry_posts: registry.map_or_else(Vec::new, |read| {
                super::formatter_registry::registry_posts(read, &lower)
            }),
            browser_url: tax_browser_url(registry),
            lower,
            answers: Default::default(),
            open: Vec::new(),
            looked_up: false,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:979-1034
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:1050-1055
    // A file of the database that does not open (`SeqDbError::TaxdbThrows`) throws `CMemoryFile`'s
    // exception at every lookup (and, with `-db`, at every load of a subject with a positive taxid
    // before it), which the report's `catch` takes as it takes a taxid that is not found: it leaves
    // the HSP out and writes nothing (review 1 of line WQA D18b, L1).
    /// `GetTaxInfo` of `taxid`: what it finds (`None`: it throws, also where a file of the database
    /// does not open). The database is opened at the
    /// run's first lookup (the texts of opening, with those of an open by earlier rows, go before
    /// the batch's first lookup's), each taxid is read once per batch. An open or a lookup that
    /// LOSAT's reader refuses is LOSAT's error (`TaxNamesDb::error`).
    fn look_up(&mut self, taxid: i32) -> Result<Option<&TaxNames>> {
        if !self.answers.contains_key(&taxid) {
            let answer = match &self.db {
                None => (None, Vec::new()),
                Some(db) => {
                    let opened = if self.looked_up {
                        Ok(())
                    } else {
                        if let Some(before) = &db.opened_before {
                            self.open = std::mem::take(&mut *before.lock().expect("open texts"));
                        }
                        db.taxdb.open().map(|opened| self.open.extend(opened))
                    };
                    let mut messages = Vec::new();
                    match opened.and_then(|()| db.taxdb.names(taxid, &mut messages)) {
                        Ok(found) => (found, messages),
                        // A lookup that throws where a file does not open (above).
                        Err(crate::seqdb::SeqDbError::TaxdbThrows { .. }) => (None, Vec::new()),
                        Err(error) => return Err((db.error)(error, &self.program)),
                    }
                }
            };
            self.answers.insert(taxid, answer);
        }
        self.looked_up = true;
        Ok(self.answers[&taxid].0.as_ref())
    }

    /// The standard-error lines of one lookup of `taxid`, after the texts of opening where it is
    /// the batch's first (`first`).
    fn messages(&self, taxid: i32, first: bool) -> Vec<u8> {
        let mut text = Vec::new();
        let mut push = |line: &str| {
            text.extend_from_slice(
                crate::seqdb::taxdb::err_post_line(&self.lower, line).as_bytes(),
            );
            text.push(b'\n');
        };
        if first {
            self.open.iter().for_each(|line| push(line));
        }
        if let Some((_, messages)) = self.answers.get(&taxid) {
            messages.iter().for_each(|line| push(line));
        }
        text
    }
}

/// `./.ncbirc` as NCBI's organism report reads it (`CTaxFormat::x_InitTaxFormat` and
/// `CAlignFormatUtil::GetProtocol` read the file of the current directory each time a report is
/// made: `formatter_registry::read_cwd_ncbirc`): `None` without the file. LOSAT refuses a file
/// whose read throws (NCBI's report then fails) and one that reads other files; the messages of a
/// read that does not throw are NCBI's posts of each report (`formatter_registry::registry_posts`).
fn read_ncbirc(program: &str) -> Result<Option<crate::blastinput::ncbi_registry::RegistryRead>> {
    super::formatter_registry::read_cwd_ncbirc().map_err(|what| {
        anyhow::Error::new(TaxReportRefusal(format!(
            "the registry file ./.ncbirc {what}, which NCBI BLAST+'s organism report (outfmt 18) reads for each query; this is not supported by LOSAT's {program}"
        )))
    })
}

/// `m_TaxBrowserURL`: the browser URL of `./.ncbirc` (`registry`,
/// `formatter_registry::tax_browser_template`) with `<@protocol@>` replaced by NCBI's scheme
/// (`formatter_registry::protocol`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:176-178
fn tax_browser_url(registry: Option<&crate::blastinput::ncbi_registry::RegistryRead>) -> Vec<u8> {
    use super::formatter_registry::{protocol, tax_browser_template};
    Text::of(tax_browser_template(registry))
        .map_bytes("protocol", protocol(registry))
        .flat()
}

/// The organism line of `names` before its centring and the replacements after it
/// (`x_MapTaxInfoTemplate` up to `AddSpaces`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:1160-1170
/// The template (`ORG_REPORT_TXT_ORG_HEADER_NO_TAX_CON`, taxFormat.cpp:98) has no
/// `<@blast_name_link@>`: its replacement changes nothing.
fn organism_line(names: &TaxNames) -> Text {
    let common: Vec<u8> = if names.common_name != names.scientific_name {
        [b"(", names.common_name.as_slice(), b")"].concat()
    } else {
        Vec::new()
    };
    let line = Text::of(b"<@scientific_name@> <@common_name@> [<@blast_name@>]")
        .map_bytes("scientific_name", &names.scientific_name)
        .map_bytes("common_name", &common)
        .map_bytes("blast_name", &names.blast_name)
        .flat();
    add_spaces(&line, u64::from(LINE_LENGTH), Place::Center, true, true)
}

/// A row of the table (`SSeqInfo` as the text report shows it).
struct Row<'s> {
    label: &'s [u8],
    title: &'s [u8],
    bit_score: std::rc::Rc<[u8]>,
    evalue: std::rc::Rc<[u8]>,
}

/// The rows of one taxid (`STaxInfo`), in the order of the HSPs, under its organism line
/// (`organism_line`).
struct Group<'s> {
    taxid: i32,
    line: Text,
    rows: Vec<Row<'s>>,
}

/// The column widths (`m_MaxAccLength`, `m_MaxScoreLength`, `m_MaxEvalLength`,
/// `m_MaxDescrLength`), NCBI's `unsigned int`s, 0 until a definition line is shown.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Widths {
    acc: u32,
    score: u32,
    evalue: u32,
    descr: u32,
}

impl Widths {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:898-909
    /// The widths after a definition line of the label `label`, with the scores `bit_score` and
    /// `evalue` (the width of the description, set between, is overwritten).
    fn take(&mut self, label: &[u8], bit_score: &[u8], evalue: &[u8]) {
        let len = |bytes: &[u8]| bytes.len() as u32;
        self.acc = self.acc.max(len(label)).max(len(b"Accession"));
        self.score = self.score.max(len(bit_score)).max(len(b"Score"));
        self.evalue = self.evalue.max(len(evalue)).max(len(b"E-value"));
        self.descr = LINE_LENGTH
            .wrapping_sub(self.acc)
            .wrapping_sub(self.score)
            .wrapping_sub(self.evalue)
            .wrapping_sub(4);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-181
// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:204-205
/// One query of a report: its title (outfmt 0's `Query=` text), the length of its record, and
/// its HSPs in report order; `None` where NCBI did not search the query's batch (no query of it
/// has a valid context), whose result NCBI gives a null Seq-align-set (a searched query without
/// hits has an empty one).
pub struct TaxQuery<'a> {
    pub title: &'a [u8],
    pub length: usize,
    pub hsps: Option<Vec<TaxHsp>>,
}

/// What NCBI's report of a query does that can write to standard error, in its order
/// (`organism_report`): the read of a subject's Bioseq at each Seq-align it does not skip
/// (`GetBioseqHandle`: where the batch's prefetch is off, the first load of a blob writes the texts
/// of its taxonomy lookups, `QueryWarnings::subject_read`) and the lookup of a taxid
/// (`GetTaxInfo`, `TaxLookups::messages`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Step {
    /// The read of the subject (`TaxHsp::subject`).
    Read(usize),
    /// The lookup of the taxid.
    LookUp(i32),
}

/// The organism reports of one query batch, made before its first byte (`TaxReport::new`), so
/// that LOSAT's refusals come before it: per query the reads and lookups of its report, in order
/// (`Step`), and its text (`None`: a query that NCBI did not search, `TaxQuery::hsps`).
pub struct TaxReport {
    lookups: TaxLookups,
    queries: Vec<Option<(Vec<Step>, Text)>>,
    /// Whether the batch's first lookup has written its texts (the texts of opening with it).
    wrote_first: std::cell::Cell<bool>,
}

impl TaxReport {
    /// The reports of the queries `queries` of a batch, of the subjects `subjects`, whose taxids
    /// the lookups `lookups` look up (`organism_report`).
    pub fn new(
        queries: &[TaxQuery<'_>],
        subjects: &dyn TaxSubjects,
        mut lookups: TaxLookups,
    ) -> Result<Self> {
        let queries = queries
            .iter()
            .map(|query| {
                query
                    .hsps
                    .as_ref()
                    .map(|hsps| organism_report(hsps, subjects, &mut lookups))
                    .transpose()
            })
            .collect::<Result<_>>()?;
        Ok(TaxReport {
            lookups,
            queries,
            wrote_first: std::cell::Cell::new(false),
        })
    }

    /// Writes the report of query `q_idx` (of the `queries` given to `new`): the caption, the
    /// `Query=` block of its title `title` (outfmt 0's, `pairwise::write_blastp_query_header`)
    /// with its record's length `length`, and its organism report; the texts of the report's reads
    /// of `./.ncbirc`, of its subjects' reads and of its lookups go to `diagnostics` (after `out`
    /// is flushed), in the order of NCBI's report (`Step`). `reads`: where the batch's prefetch is
    /// off, the texts of the subjects that the report reads one by one and the query's index among
    /// them (`ReadCursor::subject_read`). A query that NCBI did not search (`TaxQuery::hsps`) is
    /// refused after its `Query=` block (`TaxReportRefusal`), where NCBI's report stops.
    pub fn write_query(
        &self,
        q_idx: usize,
        out: &mut crate::cli::NcbiOstream<'_>,
        diagnostics: &mut dyn Write,
        mut reads: Option<(&mut ReadCursor<'_>, usize)>,
        title: &[u8],
        length: usize,
    ) -> Result<()> {
        // The report's insertions are NCBI's (`cli::insert`, `Text::insert_into`): `reportCaption`
        // and a line feed into `m_Outfile` (module docs), the `Query=` block's and `orgReportData`
        // into `out` (`organism_report`).
        add_spaces(
            b"Tax BLAST report",
            u64::from(LINE_LENGTH),
            Place::Center,
            true,
            true,
        )
        .insert_into(out)?;
        crate::cli::insert(out, b"\n")?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:885-888
        // outfmt 0's `Query=` block (`kFormatLineLength`, 68).
        crate::report::pairwise::write_blastp_query_header(out, title, length)?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:893-899
        // The Seq-align-set of a query whose batch NCBI did not search is null: `*new_aln_set`
        // throws (`CCoreException::eNullPtr`, "Attempt to access NULL pointer", exit status 255,
        // with a stack trace of NCBI's build; line RG T2 review 5, fix 5 `T/t2/fix5/probe_d2.sh`),
        // after the `Query=` block, before `CTaxFormat` reads `./.ncbirc`. No result to check: LOSAT
        // refuses there (PD-LOSAT-NCBI-DEFECTS rule 3).
        let Some((steps, text)) = &self.queries[q_idx] else {
            return Err(TaxReportRefusal(format!(
                "the organism report (outfmt 18) of a query that is not searched (no query of its batch has a valid context: \"{}\"), where NCBI BLAST+ stops with a null pointer exception; this is not supported by LOSAT's {}",
                crate::report::query_warnings::INVALID_QUERY_MESSAGE,
                self.lookups.program
            ))
            .into());
        };
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:165-173
        // The report's two reads of `./.ncbirc` (`x_InitTaxFormat`) post before its reads and
        // lookups (`x_InitTaxInfoMap`, taxFormat.cpp:140-141).
        let lookups = &self.lookups;
        if !lookups.registry_posts.is_empty() {
            out.flush()?;
            for _ in 0..2 {
                diagnostics.write_all(&lookups.registry_posts)?;
            }
            diagnostics.flush()?;
        }
        for step in steps {
            match *step {
                Step::Read(subject) => {
                    if let Some((cursor, query)) = reads.as_mut() {
                        cursor.subject_read(*query, subject as u32, out, diagnostics)?;
                    }
                }
                Step::LookUp(taxid) => {
                    let text = lookups.messages(taxid, !self.wrote_first.replace(true));
                    if !text.is_empty() {
                        out.flush()?;
                        diagnostics.write_all(&text)?;
                        diagnostics.flush()?;
                    }
                }
            }
        }
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:515
        Ok(text.insert_into(out)?)
    }
}

/// Writes the reports of `queries` (`TaxReport::write_query` of `report`), each query's warnings
/// before its report and the warnings after the last query at the end (`QueryWarnings`); the
/// texts of the subjects' reads (`QueryWarnings::with_subject_reads`) and of the lookups go to the
/// warnings' sink (none without it).
pub fn write_report(
    writer: &mut crate::cli::NcbiOstream<'_>,
    mut warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
    queries: &[TaxQuery<'_>],
    report: &TaxReport,
) -> Result<()> {
    for (q_idx, query) in queries.iter().enumerate() {
        match warnings.as_deref_mut() {
            Some(warnings) => {
                warnings.before_query(q_idx, &mut *writer)?;
                let (reads, sink) = warnings.reads_and_sink();
                report.write_query(
                    q_idx,
                    &mut *writer,
                    sink,
                    reads.map(|cursor| (cursor, q_idx)),
                    query.title,
                    query.length,
                )?;
            }
            None => report.write_query(
                q_idx,
                &mut *writer,
                &mut io::sink(),
                None,
                query.title,
                query.length,
            )?,
        }
    }
    if let Some(warnings) = warnings {
        warnings.after_last_query(queries.len(), &mut *writer)?;
    }
    Ok(())
}

/// The organism report of one query batch of a program's adapter (`algorithm::common::tabular`,
/// BLASTP's, BLASTX's): the table of the subjects of `queries`' HSPs (`SubjectTable::new` of
/// `parsed`, `plain` and `protein`), the lookups of the run's database that
/// `names` carries (`TaxLookups::new`) and the reports (`TaxReport::new`), all made before the
/// batch's first byte; `loads` what the reports load of each subject (`None`: every subject).
/// `program` names the program in LOSAT's refusals and NCBI's messages.
#[allow(clippy::too_many_arguments)]
pub fn report_batch<'a>(
    queries: &[TaxQuery<'_>],
    parsed: &[Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
    plain: impl Fn(usize) -> TaxPlain<'a>,
    protein: bool,
    loads: Option<&crate::report::db_subject::HitLoads>,
    names: Option<&TabularNames>,
    program: &str,
) -> Result<TaxReport> {
    let table = SubjectTable::new(queries, parsed, plain, protein, loads, program)
        .with_loads_that_throw(
            names
                .and_then(TabularNames::report_db)
                .is_some_and(|db| db.loads_throw),
        );
    TaxReport::new(queries, &table, TaxLookups::new(names, program)?)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:979-1034
// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:1040-1055
/// The organism report of `hsps` (`CTaxFormat::x_InitTaxInfoMap`, then `DisplayOrgReport`) and
/// its reads and lookups, in order (`Step`): the read of the subject of each HSP that it does not
/// skip, then the lookups of its definition lines (`lookups`: a taxid not yet in the query's map,
/// at each definition line that reaches it). A title that throws drops the rest of its subject's
/// HSP before it widens the columns; a lookup that throws drops it after; either way the subject
/// is not the previous id, so its next HSP is read again. An HSP whose subject's load throws
/// (`TaxSubjects::throws`, the load of `GetBioseqHandle` where its batch's prefetch is off) is
/// left out at the load, and is not the previous id either. LOSAT's refusals: a title that it does
/// not write as NCBI where NCBI makes it (or, `shown_only`, shows it), and an organism line that
/// holds `<@bl_taxid@>` (NCBI writes an uninitialised `int` there).
fn organism_report(
    hsps: &[TaxHsp],
    subjects: &dyn TaxSubjects,
    lookups: &mut TaxLookups,
) -> Result<(Vec<Step>, Text)> {
    let mut widths = Widths::default();
    let mut groups: Vec<Group<'_>> = Vec::new();
    let mut group_of: std::collections::HashMap<i32, usize> = Default::default();
    let mut steps = Vec::new();
    let mut previous: Option<usize> = None;
    for hsp in hsps {
        if previous.is_some_and(|previous| subjects.same(previous, hsp.subject)) {
            continue;
        }
        // `m_Scope.GetBioseqHandle` of the id `subid` (taxFormat.cpp:997), before the definition lines; a
        // load that throws goes to the `catch`, `continue`, and reads nothing.
        if subjects.throws(hsp.subject) {
            continue;
        }
        steps.push(Step::Read(hsp.subject));
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:950-951
        let evalue: std::rc::Rc<[u8]> = format_evalue_ncbi(hsp.evalue).into_bytes().into();
        let bit_score: std::rc::Rc<[u8]> = format_bitscore_ncbi(hsp.bit_score).into_bytes().into();
        let mut threw = false;
        for defline in subjects.deflines(hsp.subject) {
            // The title, or the refusal of a title that the report shows.
            let (title, refused_if_shown) = match &defline.title {
                TaxTitle::Text(title) => (title.as_slice(), None),
                TaxTitle::Throws => {
                    threw = true;
                    break;
                }
                TaxTitle::Refused {
                    reason,
                    shown_only: false,
                } => return Err(subjects.refused(hsp.subject, reason)),
                TaxTitle::Refused { reason, .. } => (&[][..], Some(reason)),
            };
            widths.take(&defline.label, &bit_score, &evalue);
            let index = match group_of.get(&defline.taxid) {
                Some(&index) => index,
                None => {
                    steps.push(Step::LookUp(defline.taxid));
                    let Some(names) = lookups.look_up(defline.taxid)? else {
                        threw = true;
                        break;
                    };
                    let line = organism_line(names);
                    if line.contains(b"<@bl_taxid@>") {
                        return Err(TaxReportRefusal(format!(
                            "the taxonomy database names taxid {} with the template token <@bl_taxid@>, which NCBI BLAST+'s organism report (outfmt 18) replaces with an uninitialised number; this is not supported by LOSAT's {}",
                            defline.taxid, lookups.program
                        ))
                        .into());
                    }
                    group_of.insert(defline.taxid, groups.len());
                    groups.push(Group {
                        taxid: defline.taxid,
                        line,
                        rows: Vec::new(),
                    });
                    groups.len() - 1
                }
            };
            if let Some(reason) = refused_if_shown {
                return Err(subjects.refused(hsp.subject, reason));
            }
            groups[index].rows.push(Row {
                label: &defline.label,
                title,
                bit_score: bit_score.clone(),
                evalue: evalue.clone(),
            });
        }
        if !threw {
            previous = Some(hsp.subject);
        }
    }
    let text = org_report_text(&groups, widths, &lookups.browser_url);
    Ok((steps, text))
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:450-517 (`DisplayOrgReport`)
// `x_InitTaxIdToSeqsMap` fills a map that only the HTML report shows.
/// The text of the organism report of `groups` with the widths `widths`; `browser_url` replaces
/// `<@taxBrowserURL@>`.
fn org_report_text(groups: &[Group<'_>], widths: Widths, browser_url: &[u8]) -> Text {
    const DISABLED: &[u8] = b"disabled=\"disabled\"";
    let mut data = Text::default();
    for (i, group) in groups.iter().enumerate() {
        let last = i + 1 == groups.len();
        let number = |index: usize| groups[index].taxid.to_string().into_bytes();
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:1172-1187
        // `<@bl_taxid@>` is refused where the group is made (`organism_report`); the line holds none.
        let header = group
            .line
            .map_bytes("taxid", group.taxid.to_string().as_bytes())
            .map_bytes("taxBrowserURL", browser_url)
            .map_bytes("rid", b"0")
            .map_bytes("numhits", group.rows.len().to_string().as_bytes())
            .map_bytes("depth", b"")
            .map_bytes("next_taxid", &if last { Vec::new() } else { number(i + 1) })
            .map_bytes("disable_nexttaxid", if last { DISABLED } else { b"" })
            .map_bytes(
                "prev_taxid",
                &if i == 0 { Vec::new() } else { number(i - 1) },
            )
            .map_bytes("disable_prevtaxid", if i == 0 { DISABLED } else { b"" })
            .map_bytes("disable_top", if i == 0 { DISABLED } else { b"" });
        data.push_text(&header);
        for row in &group.rows {
            data.push_text(&row_text(row, widths));
        }
    }
    Text::of(
        b"<@org_report_caption@>\n<@acc_hd@><@descr_hd@><@score_hd@><@evalue_hd@>\n<@table_rows@>",
    )
    .map("table_rows", &data)
    .map(
        "org_report_caption",
        &add_spaces(
            b"Organism Report",
            u64::from(LINE_LENGTH),
            Place::Center,
            true,
            false,
        ),
    )
    .map(
        "acc_hd",
        &add_spaces(
            b"Accession",
            u64::from(widths.acc),
            Place::End,
            false,
            false,
        ),
    )
    .map(
        "descr_hd",
        &add_spaces(
            b"Description",
            u64::from(widths.descr),
            Place::Center,
            false,
            false,
        ),
    )
    .map(
        "score_hd",
        &add_spaces(b"Score", u64::from(widths.score), Place::End, false, false),
    )
    .map(
        "evalue_hd",
        &add_spaces(
            b"E-value",
            u64::from(widths.evalue),
            Place::End,
            false,
            false,
        ),
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:1130-1150
/// One row (`ORG_REPORT_TXT_TABLE_ROW`, ` <@acc@><@descr_text@><@score@><@evalue@>\n`). The
/// replacements of `gi`, `disp_gi`, `descr_abbr` and `rid` come before any value is put in and
/// the template holds none of them: they change nothing.
fn row_text(row: &Row<'_>, widths: Widths) -> Text {
    let column =
        |value: &[u8], width: u32| add_spaces(value, u64::from(width), Place::End, false, false);
    Text::of(b" <@acc@><@descr_text@><@score@><@evalue@>\n")
        .map("acc", &column(row.label, widths.acc))
        .map("descr_text", &column(row.title, widths.descr))
        .map("score", &column(&row.bit_score, widths.score))
        .map("evalue", &column(&row.evalue, widths.evalue))
}

#[cfg(test)]
#[path = "tax_report_tests.rs"]
mod tests;
