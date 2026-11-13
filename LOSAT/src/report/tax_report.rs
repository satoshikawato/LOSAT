//! The organism report (outfmt 18, NCBI's "Tax BLAST report", `CTaxFormat` in text mode without
//! the taxonomy server): for each query a caption, the outfmt 0 `Query=` block, and a table of
//! the query's subjects grouped by the taxid of each definition line, under the names that the
//! taxonomy database gives for it. One writer for every program; a program's adapter hands over,
//! for each query, its title and length and its HSPs in report order (the final hit list as the
//! tabular reports read it, `report::tabular`), each with its subject, and the definition lines
//! of the subjects as the report shows them (`TaxDefline`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1461-1467
//! ```c++
//!     if (m_FormatType == CFormattingArgs::eTaxFormat) {
//!         string reportCaption = "Tax BLAST report";
//!         reportCaption = m_IsHTML ? "<h1>" + reportCaption + "</h1>" : CAlignFormatUtil::AddSpaces(reportCaption,max(kMinTaxFormatLineLength,kFormatLineLength),CAlignFormatUtil::eSpacePosToCenter | CAlignFormatUtil::eAddEOLAtLineStart | CAlignFormatUtil::eAddEOLAtLineEnd);
//!         m_Outfile << reportCaption;
//!         x_PrintTaxReport(results);
//!         return;
//!     }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:873-902
//! ```c++
//! CBlastFormat::x_PrintTaxReport(const blast::CSearchResults& results)
//! {
//!     ...
//!         m_Outfile <<  "\n";
//!     ...
//!     CBlastFormatUtil::AcknowledgeBlastQuery(*bioseq, kFormatLineLength,
//!                                             m_Outfile, m_BelieveQuery,
//!                                             m_IsHTML, false,
//!                                             results.GetRID());
//!     ...
//!     CConstRef<CSeq_align_set> aln_set = results.GetSeqAlign();
//!     if (m_IsUngappedSearch && results.HasAlignments()) {
//!         aln_set.Reset(CDisplaySeqalign::PrepareBlastUngappedSeqalign(*aln_set));
//!     }
//!     ...
//!     CTaxFormat *taxFormatRes = new CTaxFormat (*new_aln_set,*m_Scope,m_IsHTML ? CTaxFormat::eHtml:CTaxFormat::eText, false,max(kMinTaxFormatLineLength,kFormatLineLength));
//!     taxFormatRes->DisplayOrgReport(m_Outfile);
//! }
//! ```
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
use crate::report::tabular::TabularNames;
use crate::seqdb::taxdb::TaxNames;

/// `max(kMinTaxFormatLineLength, kFormatLineLength)`: the report's line length (`m_LineLength`).
///
/// NCBI reference (598d8ae6): c++/include/objtools/align_format/taxFormat.hpp:58
/// ```c++
/// const unsigned int kMinTaxFormatLineLength = 100;
/// ```
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

    /// `CAlignFormatUtil::MapTemplate(text, name, value)`: every `<@name@>` replaced by `value`,
    /// left to right; the text put in is not searched again by this call.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3593-3599
    /// ```c++
    /// string CAlignFormatUtil::MapTemplate(string inpString,string tmplParamName,string templParamVal)
    /// {
    ///     string outString;
    ///     string tmplParam = "<@" + tmplParamName + "@>";
    ///     NStr::Replace(inpString,tmplParam,templParamVal,outString);
    ///     return outString;
    /// }
    /// ```
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

    /// The length in bytes.
    #[cfg(test)]
    fn len(&self) -> u64 {
        self.0
            .iter()
            .map(|seg| match seg {
                Seg::Bytes(bytes) => bytes.len() as u64,
                Seg::Spaces(count) => *count,
            })
            .sum()
    }

    fn write_to(&self, out: &mut dyn Write) -> io::Result<()> {
        const BLANK: [u8; 1 << 16] = [b' '; 1 << 16];
        for seg in &self.0 {
            match seg {
                Seg::Bytes(bytes) => out.write_all(bytes)?,
                Seg::Spaces(count) => {
                    let mut left = *count;
                    while left > 0 {
                        let n = left.min(BLANK.len() as u64) as usize;
                        out.write_all(&BLANK[..n])?;
                        left -= n as u64;
                    }
                }
            }
        }
        Ok(())
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

/// `CAlignFormatUtil::AddSpaces(value, width, place)` with `eAddEOLAtLineStart` (`eol_start`)
/// and `eAddEOLAtLineEnd` (`eol_end`). `width` is NCBI's `size_t`: a width of 0 to 2 cuts
/// nothing (`width - 3` wraps) and still adds `...`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3610-3640
/// ```c++
/// string CAlignFormatUtil::AddSpaces(string paramVal, size_t maxParamValLength, int spacesFormatFlag)
/// {
///     //if(!spacePos.empty()) {
///         string spaceString;
///         if(maxParamValLength >= paramVal.size()) {
///             size_t numSpaces = maxParamValLength - paramVal.size() + 1;
///             if(spacesFormatFlag & eSpacePosToCenter) {
///                 numSpaces = numSpaces/2;
///             }
/// 	    spaceString.assign(numSpaces,' ');
///         }
///         else {
///             paramVal = paramVal.substr(0, maxParamValLength - 3) + "...";
///             spaceString += " ";
///         }
///         if(spacesFormatFlag & eSpacePosAtLineEnd) {
///             paramVal = paramVal + spaceString;
///         }
///         else if(spacesFormatFlag & eSpacePosToCenter) {
///             paramVal = spaceString + paramVal + spaceString;
///         }
///         else {
///             paramVal = spaceString + paramVal;
///         }
///         if(spacesFormatFlag & eAddEOLAtLineStart) paramVal = "\n" + paramVal;
///         if(spacesFormatFlag & eAddEOLAtLineEnd) paramVal = paramVal + "\n";
///     //}
///
///     return paramVal;
/// }
/// ```
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

/// One definition line of a subject as the report shows it (`CTaxFormat::SSeqInfo`): the label
/// (`CAlignFormatUtil::GetLabel` of its `WorstRank`-best Seq-id, without version), the title (the
/// generated title of the subject's Bioseq where it has none) and the taxid (0 for none).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TaxDefline {
    pub label: Vec<u8>,
    pub title: Vec<u8>,
    pub taxid: i32,
}

impl TaxDefline {
    /// The definition line of a `-subject` record read without `-parse_deflines`: its Bioseq has
    /// the local id `Subject_N` (`local_id`) and no `Blast-def-line`, so its title is the
    /// generated one of the record's title `title`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:911-929
    /// ```c++
    /// CTaxFormat::SSeqInfo *CTaxFormat::x_FillTaxDispParams(const CBioseq_Handle& bsp_handle,
    ///                                     double bits,
    ///                                     double evalue)
    /// {
    ///     SSeqInfo *seqInfo = new SSeqInfo;
    /// 	seqInfo->gi = FindGi(bsp_handle.GetBioseqCore()->GetId());
    /// 	seqInfo->seqID = FindBestChoice(bsp_handle.GetBioseqCore()->GetId(),CSeq_id::WorstRank);
    /// 	seqInfo->label =  CAlignFormatUtil::GetLabel(seqInfo->seqID);//Just accession without db part like ref| or pdbd|
    ///     ...
    ///     seqInfo->taxid = ZERO_TAX_ID;
    /// 	seqInfo->title = CDeflineGenerator().GenerateDefline(bsp_handle);
    /// ```
    pub fn of_record(local_id: &[u8], title: &[u8], protein: bool) -> Result<Self, String> {
        Ok(TaxDefline {
            label: local_id.to_vec(),
            title: generated_title(
                title,
                protein,
                &crate::report::defline::BioseqIds::default(),
            )?,
            taxid: 0,
        })
    }

    /// The definition line of the database sequence `ordinal` without parsed Seq-ids: the label
    /// of its `gnl|BL_ORD_ID|<ordinal>` id is the ordinal id, its title the stored one (`title`)
    /// or, where that is empty, the Bioseq's generated title. `Err` where NCBI makes that title
    /// from the taxid and the taxonomy database (`taxdb_found`), which LOSAT does not.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:962-971
    /// ```c++
    /// 		if(bdl->IsSetTaxid() &&  bdl->CanGetTaxid()){
    /// 		    taxid = bdl->GetTaxid();
    /// 		}
    ///         seqInfo->taxid = taxid;
    ///
    /// 		if(bdl->IsSetTitle()){
    /// 			seqInfo->title = bdl->GetTitle();
    /// 		}
    ///         if(seqInfo->title.empty()) {
    ///             seqInfo->title = CDeflineGenerator().GenerateDefline(bsp_handle);
    ///         }
    /// ```
    pub fn of_stored(
        ordinal: u32,
        title: &[u8],
        taxid: i32,
        protein: bool,
        taxdb_found: bool,
    ) -> Result<Self, String> {
        let ids = crate::report::defline::BioseqIds::default();
        let title = if !title.is_empty() {
            title.to_vec()
        } else if taxid > 0
            && taxdb_found
            && crate::report::defline::bioseq_title_from_source(title, &ids)
        {
            return Err(TITLE_FROM_TAXONOMY.into());
        } else {
            generated_title(title, protein, &ids)?
        };
        Ok(TaxDefline {
            label: ordinal.to_string().into_bytes(),
            title,
            taxid,
        })
    }
}

/// LOSAT's refusal of a generated title that NCBI makes from the taxonomy database's names.
pub(crate) const TITLE_FROM_TAXONOMY: &str = "has an empty title that NCBI BLAST+'s organism report (outfmt 18) makes from its taxids and the taxonomy database";

/// `CDeflineGenerator().GenerateDefline` (flags 0) of a Bioseq with the title `title` and the
/// Seq-ids `ids`. `Err`: a protein title of punctuation that NCBI's `x_CleanAndCompress` reads
/// past the end of (NCBI crashes; approved exception 2 of PD-LOSAT-NCBI-DEFECTS has a nucleotide
/// title stop at the end of the string), and a title whose encoding NCBI's `HtmlDecode` cannot
/// guess (NCBI's report then skips the subject).
pub(crate) fn generated_title(
    title: &[u8],
    protein: bool,
    ids: &crate::report::defline::BioseqIds,
) -> Result<Vec<u8>, String> {
    match crate::report::defline::generate_bioseq_defline(title, protein, false, ids) {
        Ok(defline) if protein && defline.reads_past_end => Err("has a title of punctuation that NCBI BLAST+'s x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes)".into()),
        Ok(defline) => Ok(defline.text),
        Err(_) => Err("has a title whose encoding NCBI BLAST+'s HtmlDecode cannot guess (its organism report then skips the subject)".into()),
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

/// The subjects of a report, as a program's adapter gives them.
pub trait TaxSubjects {
    /// The definition lines of `subject` that the report shows (those a positive list keeps,
    /// `CAlignFormatUtil::MatchSeqInSeqList`; one for a subject without `Blast-def-line`s).
    fn deflines(&self, subject: usize) -> &[TaxDefline];

    /// Whether the Seq-align id of `subject` matches that of `previous` (`CSeq_id::Match`): the
    /// report skips such an HSP after one of `previous` that it took.
    fn same(&self, previous: usize, subject: usize) -> bool {
        previous == subject
    }
}

/// LOSAT's refusal of an organism report that it does not write as NCBI (a title, a taxonomy
/// name, a `./.ncbirc`), made before the report's first byte: no error of NCBI's engine.
#[derive(Debug)]
pub struct TaxReportRefusal(String);

impl std::fmt::Display for TaxReportRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for TaxReportRefusal {}

/// Where the definition lines of a subject come from (`subject_deflines`).
#[derive(Clone, Copy, Debug)]
pub enum TaxSource<'a> {
    /// A subject with parsed Seq-ids: the database sequence `ordinal`, or a `-subject` record read
    /// under `-parse_deflines` (`ordinal` `None`).
    Parsed {
        subject: &'a crate::report::db_subject::ParsedSubject,
        ordinal: Option<u32>,
    },
    /// A `-subject` record read without `-parse_deflines`: its local id `Subject_N` and title.
    Record { local_id: &'a [u8], title: &'a [u8] },
    /// The database sequence `ordinal` without parsed Seq-ids: its stored title and taxid.
    Stored {
        ordinal: u32,
        title: &'a [u8],
        taxid: i32,
    },
}

/// The definition lines of the subject `source` as the report shows them (`TaxDefline`), of a
/// protein (`protein`) or nucleotide subject; `taxdb_found`: whether `taxdb.bti` is on the search
/// path. A title that LOSAT does not write as NCBI is LOSAT's refusal (`program` names the
/// program).
pub fn subject_deflines(
    source: TaxSource<'_>,
    protein: bool,
    taxdb_found: bool,
    program: &str,
) -> Result<Vec<TaxDefline>> {
    let refused = |reason: String| {
        anyhow::Error::new(TaxReportRefusal(match source {
            TaxSource::Parsed {
                ordinal: Some(ordinal),
                ..
            }
            | TaxSource::Stored { ordinal, .. } => format!(
                "database sequence {ordinal} {reason}; this is not supported by LOSAT's {program}"
            ),
            TaxSource::Parsed { subject, .. } => format!(
                "the subject record {} {reason}, which is not supported by LOSAT's {program}",
                String::from_utf8_lossy(&subject.alignment_fasta())
            ),
            TaxSource::Record { title, .. } => format!(
                "the subject record {} {reason}, which is not supported by LOSAT's {program}",
                String::from_utf8_lossy(title)
            ),
        }))
    };
    match source {
        TaxSource::Parsed { subject, .. } => subject.tax_deflines(protein, taxdb_found),
        TaxSource::Record { local_id, title } => {
            TaxDefline::of_record(local_id, title, protein).map(|line| vec![line])
        }
        TaxSource::Stored {
            ordinal,
            title,
            taxid,
        } => TaxDefline::of_stored(ordinal, title, taxid, protein, taxdb_found)
            .map(|line| vec![line]),
    }
    .map_err(refused)
}

/// The subjects of a report by subject index (a subject set's index, a database's rank): their
/// definition lines, made before the report's first byte, and the parsed subjects, whose Seq-align
/// ids `Match` where they hold equal ids (`db_subject::joined_in_reports`).
pub struct SubjectTable<'p> {
    lines: std::collections::HashMap<usize, Vec<TaxDefline>>,
    parsed: &'p [Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
}

impl<'p> SubjectTable<'p> {
    /// An empty table; `parsed` are the parsed subjects by subject index (empty when none).
    pub fn new(
        parsed: &'p [Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
    ) -> Self {
        SubjectTable {
            lines: Default::default(),
            parsed,
        }
    }

    /// Adds the definition lines of `subject` by `make` unless it has them.
    pub fn add(
        &mut self,
        subject: usize,
        make: impl FnOnce() -> Result<Vec<TaxDefline>>,
    ) -> Result<()> {
        if let std::collections::hash_map::Entry::Vacant(entry) = self.lines.entry(subject) {
            entry.insert(make()?);
        }
        Ok(())
    }

    /// The parsed subject of `subject`, if it has one.
    pub fn parsed(&self, subject: usize) -> Option<&'p crate::report::db_subject::ParsedSubject> {
        self.parsed.get(subject).and_then(Option::as_deref)
    }
}

impl TaxSubjects for SubjectTable<'_> {
    fn deflines(&self, subject: usize) -> &[TaxDefline] {
        self.lines.get(&subject).map_or(&[], Vec::as_slice)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:991-993
    // ```c++
    //         subid = &((*iter)->GetSeq_id(1));
    //
    //         if(!previousId.Empty() && subid->Match(*previousId)) continue;
    // ```
    fn same(&self, previous: usize, subject: usize) -> bool {
        previous == subject
            || matches!(
                (self.parsed(previous), self.parsed(subject)),
                (Some(a), Some(b)) if crate::report::db_subject::joined_in_reports(a, b)
            )
    }
}

/// The lookups of the report's taxids (`CSeqDB::GetTaxInfo`): what each finds, and the text it
/// writes to standard error each time it is called.
pub struct TaxLookups<'a> {
    names: Option<&'a TabularNames>,
    /// The names of `[TAX_BROWSER] BLASTFMTUTIL`'s URL (`<@taxBrowserURL@>`), read where a name
    /// holds the token.
    browser_url: Vec<u8>,
}

impl<'a> TaxLookups<'a> {
    /// The lookups of `names` (`None`: every lookup fails silently, no `taxdb.bti` on the search
    /// path). LOSAT's refusals of the report come here, before its first byte: a name that holds
    /// `<@bl_taxid@>` (NCBI writes an uninitialised `int` there) or `<@taxBrowserURL@>` with a
    /// `./.ncbirc` that LOSAT does not read as NCBI's report does (`browser_url`), and a
    /// `./.ncbirc` that NCBI's report cannot read (it reads the file for every query).
    pub fn new(names: Option<&'a TabularNames>, program: &str) -> Result<Self> {
        let registry = read_ncbirc(program)?;
        let mut browser_url = Vec::new();
        if let Some(names) = names {
            for found in names.found() {
                let line = organism_line(found);
                if line.contains(b"<@bl_taxid@>") {
                    return Err(TaxReportRefusal(format!(
                        "the taxonomy database names taxid {} with the template token <@bl_taxid@>, which NCBI BLAST+'s organism report (outfmt 18) replaces with an uninitialised number; this is not supported by LOSAT's {program}",
                        found.taxid
                    ))
                    .into());
                }
                if line.contains(b"<@taxBrowserURL@>") && browser_url.is_empty() {
                    browser_url = tax_browser_url(registry.as_ref());
                }
            }
        }
        Ok(TaxLookups { names, browser_url })
    }

    /// What `GetTaxInfo(taxid)` finds (`None`: it throws), and the text it writes to standard
    /// error.
    fn look_up(&self, taxid: i32) -> (Option<&'a TaxNames>, Vec<u8>) {
        match self.names {
            Some(names) => (names.get(taxid), names.lookup_messages(taxid)),
            None => (None, Vec::new()),
        }
    }
}

/// `./.ncbirc` as NCBI's organism report reads it (`CTaxFormat::x_InitTaxFormat` and
/// `CAlignFormatUtil::GetProtocol` read the file of the current directory each time a report is
/// made, with a registry of their own): `None` without the file. LOSAT refuses a file whose read
/// throws (NCBI's report then fails), writes a message (a stream that fails, an entry set again,
/// which NCBI writes at each report) or reads other files (`[NCBI] .Inherits`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:171-177
/// ```c++
///     //set config file
///     m_ConfigFile = new CNcbiIfstream(".ncbirc");
///     m_Reg = new CNcbiRegistry(*m_ConfigFile);
///     if(m_Reg) {
///         m_TaxBrowserURL = m_Reg->Get("TAX_BROWSER","BLASTFMTUTIL");
///     }
/// ```
fn read_ncbirc(program: &str) -> Result<Option<crate::blastinput::ncbi_registry::RegistryRead>> {
    use crate::blastinput::ncbi_registry::{read_registry, Layer};
    let path = std::path::Path::new(".ncbirc");
    let Ok(bytes) = std::fs::read(path) else {
        return Ok(None);
    };
    let refused = |what: &str| {
        anyhow::Error::new(TaxReportRefusal(format!(
            "the registry file ./.ncbirc {what}, which NCBI BLAST+'s organism report (outfmt 18) reads for each query; this is not supported by LOSAT's {program}"
        )))
    };
    let read = read_registry(&bytes, Layer::Application).map_err(|reason| refused(&reason))?;
    if read.error.is_some() || read.stream_failed {
        return Err(refused("cannot be read as a registry"));
    }
    if !read.repeated.is_empty() {
        return Err(refused("sets an entry more than once"));
    }
    if read.last_value("NCBI", ".Inherits").is_some() {
        return Err(refused("sets [NCBI] .Inherits"));
    }
    Ok(Some(read))
}

/// `m_TaxBrowserURL`: `[TAX_BROWSER] BLASTFMTUTIL` of `./.ncbirc` (`registry`), else NCBI's URL,
/// with `<@protocol@>` replaced by `[BLASTFMTUTIL] PROTOCOL` (`https:` where it is not set or
/// empty: `HasEntry` counts no empty value).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:176-178
/// ```c++
///     m_TaxBrowserURL = (m_TaxBrowserURL.empty()) ? s_TagToConstString("TAX_BROWSER_URL") : m_TaxBrowserURL;
///     m_TaxBrowserURL = CAlignFormatUtil::MapTemplate(m_TaxBrowserURL,"protocol",m_Protocol);
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3643-3655
/// ```c++
/// string CAlignFormatUtil::GetProtocol()
/// {
///     CNcbiIfstream  config_file(".ncbirc");
///     CNcbiRegistry config_reg(config_file);
///     string httpProt = "https:";
///     if(!config_reg.Empty()) {
///         if(config_reg.HasEntry("BLASTFMTUTIL","PROTOCOL")) {
///             httpProt = config_reg.Get("BLASTFMTUTIL","PROTOCOL");
///         }
///     }
///     return httpProt;
/// }
/// ```
fn tax_browser_url(registry: Option<&crate::blastinput::ncbi_registry::RegistryRead>) -> Vec<u8> {
    const DEFAULT: &[u8] = b"<@protocol@>//www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi";
    let value = |section: &str, name: &str| {
        registry
            .and_then(|read| read.last_value(section, name))
            .filter(|value| !value.is_empty())
    };
    let url = value("TAX_BROWSER", "BLASTFMTUTIL").unwrap_or(DEFAULT);
    let protocol = value("BLASTFMTUTIL", "PROTOCOL").unwrap_or(b"https:");
    Text::of(url).map_bytes("protocol", protocol).flat()
}

/// The organism line of `names` before its centring and the replacements after it
/// (`x_MapTaxInfoTemplate` up to `AddSpaces`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:1160-1170
/// ```c++
/// string CTaxFormat::x_MapTaxInfoTemplate(string tableRowTemplate,STaxInfo &taxInfo,unsigned int depth)
/// {
///     string reportTableRow = CAlignFormatUtil::MapTemplate(tableRowTemplate,"blast_name_link",m_TaxFormatTemplates->blastNameLink);
///     reportTableRow = CAlignFormatUtil::MapTemplate(reportTableRow,"scientific_name",taxInfo.scientificName);
///     string commonName = (taxInfo.commonName != taxInfo.scientificName) ? "(" + taxInfo.commonName + ")" : "";
///     reportTableRow = CAlignFormatUtil::MapTemplate(reportTableRow,"common_name",commonName);
///     reportTableRow = CAlignFormatUtil::MapTemplate(reportTableRow,"blast_name",taxInfo.blastName);
///     if(m_DisplayOption == eText) {
///         string txtOrgValue = reportTableRow;//Contains <@scientific_name@>[<@blast_name@>] mapped
///         reportTableRow = CAlignFormatUtil::AddSpaces(txtOrgValue,m_LineLength,CAlignFormatUtil::eSpacePosToCenter | CAlignFormatUtil::eAddEOLAtLineStart | CAlignFormatUtil::eAddEOLAtLineEnd);
///     }
/// ```
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

/// The rows of one taxid (`STaxInfo`), in the order of the HSPs.
struct Group<'s, 'n> {
    taxid: i32,
    names: &'n TaxNames,
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
    // ```c++
    // void CTaxFormat::x_InitTextFormatInfo(CTaxFormat::SSeqInfo *seqInfo)
    // {
    //     string orgAccTxtTableHeader = s_TagToConstString("ORG_ACC_TXT_TABLE_HEADER");
    //     m_MaxAccLength = max(max(m_MaxAccLength,(unsigned int)seqInfo->label.size()),(unsigned int)orgAccTxtTableHeader.size());
    //     ...
    //     m_MaxDescrLength =  m_LineLength - m_MaxAccLength - m_MaxScoreLength - m_MaxEvalLength - 4;//4 for spaces in between
    //  }
    // ```
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

/// One query of a report: its title (outfmt 0's `Query=` text), the length of its record, and
/// its HSPs in report order.
pub struct TaxQuery<'a> {
    pub title: &'a [u8],
    pub length: usize,
    pub hsps: Vec<TaxHsp>,
}

/// Writes the reports of `queries` (`write_query`), each query's warnings before its report and
/// the warnings after the last query at the end (`QueryWarnings`); the texts of the lookups go
/// to the warnings' sink (none without it).
pub fn write_report(
    writer: &mut dyn Write,
    mut warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
    queries: &[TaxQuery<'_>],
    subjects: &dyn TaxSubjects,
    lookups: &TaxLookups<'_>,
) -> io::Result<()> {
    let mut sink = io::sink();
    for (q_idx, query) in queries.iter().enumerate() {
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        let diagnostics: &mut dyn Write = match warnings.as_deref_mut() {
            Some(warnings) => &mut *warnings.sink,
            None => &mut sink,
        };
        write_query(
            &mut *writer,
            diagnostics,
            query.title,
            query.length,
            &query.hsps,
            subjects,
            lookups,
        )?;
    }
    if let Some(warnings) = warnings {
        warnings.after_last_query(queries.len(), &mut *writer)?;
    }
    Ok(())
}

/// Writes the report of one query: the caption, the `Query=` block of the title `title` (outfmt
/// 0's, `pairwise::write_blastp_query_header`) with the length `length`, and the organism report of the
/// HSPs `hsps` (report order) of the subjects `subjects`, whose taxids `lookups` looks up; the
/// texts of the lookups go to `diagnostics` (after `out` is flushed).
pub fn write_query(
    out: &mut dyn Write,
    diagnostics: &mut dyn Write,
    title: &[u8],
    length: usize,
    hsps: &[TaxHsp],
    subjects: &dyn TaxSubjects,
    lookups: &TaxLookups<'_>,
) -> io::Result<()> {
    add_spaces(
        b"Tax BLAST report",
        u64::from(LINE_LENGTH),
        Place::Center,
        true,
        true,
    )
    .write_to(out)?;
    out.write_all(b"\n")?;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:729-746
    // (`AcknowledgeBlastQuery` with `kFormatLineLength`, 68: outfmt 0's `Query=` block).
    let mut out = out;
    crate::report::pairwise::write_blastp_query_header(&mut out, title, length)?;
    organism_report(diagnostics, out, hsps, subjects, lookups)?.write_to(out)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:979-1023
// ```c++
//     ITERATE(CSeq_align_set::Tdata, iter, m_SeqalignSetRef->Get()){
//
//         subid = &((*iter)->GetSeq_id(1));
//
//         if(!previousId.Empty() && subid->Match(*previousId)) continue;
//         ...
//         try {
//             const CBioseq_Handle& handle = m_Scope.GetBioseqHandle(*subid);
//             if( !handle ) continue;
//             ...
//             if(!bdl.empty()){ //no blast defline struct, should be no such case now
//                 ...
//                 for(list< CRef< CBlast_def_line > >::const_iterator iter = bdl.begin(); iter != bdl.end(); iter++){
//                     SSeqInfo *seqInfo = x_FillTaxDispParams(*iter,handle,bits,evalue,use_this_seqid);
//                     if(seqInfo) {
//                         ...
//                         x_InitBlastDBTaxInfo(seqInfo);
//                     }
//                 }
//             }
//             else
//             {
//                 SSeqInfo *seqInfo = x_FillTaxDispParams(handle,bits,evalue);
//                 if(seqInfo) {
//                     x_InitBlastDBTaxInfo(seqInfo);
//                 }
//             }
//             previousId = subid;
//         }
//         catch (const CException&){
//             continue;
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:1040-1055
// ```c++
// void CTaxFormat::x_InitBlastDBTaxInfo(CTaxFormat::SSeqInfo *seqInfo)
// {
//     TTaxId taxid = seqInfo->taxid;
//     if(m_BlastResTaxInfo->seqTaxInfoMap.count(taxid) == 0){
//         SSeqDBTaxInfo taxInfo;
//         CSeqDB::GetTaxInfo(taxid, taxInfo);
// ```
/// The organism report of `hsps` (`CTaxFormat::x_InitTaxInfoMap`, then `DisplayOrgReport`): a
/// subject whose taxid lookup throws keeps its rows before that definition line, widens the
/// columns, and is not the previous id, so its next HSP is read again.
fn organism_report(
    diagnostics: &mut dyn Write,
    out: &mut dyn Write,
    hsps: &[TaxHsp],
    subjects: &dyn TaxSubjects,
    lookups: &TaxLookups<'_>,
) -> io::Result<Text> {
    let mut widths = Widths::default();
    let mut groups: Vec<Group<'_, '_>> = Vec::new();
    let mut previous: Option<usize> = None;
    for hsp in hsps {
        if previous.is_some_and(|previous| subjects.same(previous, hsp.subject)) {
            continue;
        }
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:950-951
        // ```c++
        //         string total_bit_score_buf, raw_score_buf;
        //         CAlignFormatUtil::GetScoreString(evalue, bits, 0, 0, seqInfo->evalue, seqInfo->bit_score, total_bit_score_buf,raw_score_buf);
        // ```
        let evalue: std::rc::Rc<[u8]> = format_evalue_ncbi(hsp.evalue).into_bytes().into();
        let bit_score: std::rc::Rc<[u8]> = format_bitscore_ncbi(hsp.bit_score).into_bytes().into();
        let mut threw = false;
        for defline in subjects.deflines(hsp.subject) {
            widths.take(&defline.label, &bit_score, &evalue);
            let row = Row {
                label: &defline.label,
                title: &defline.title,
                bit_score: bit_score.clone(),
                evalue: evalue.clone(),
            };
            match groups.iter_mut().find(|group| group.taxid == defline.taxid) {
                Some(group) => group.rows.push(row),
                None => {
                    let (found, text) = lookups.look_up(defline.taxid);
                    if !text.is_empty() {
                        out.flush()?;
                        diagnostics.write_all(&text)?;
                        diagnostics.flush()?;
                    }
                    let Some(names) = found else {
                        threw = true;
                        break;
                    };
                    groups.push(Group {
                        taxid: defline.taxid,
                        names,
                        rows: vec![row],
                    });
                }
            }
        }
        if !threw {
            previous = Some(hsp.subject);
        }
    }
    Ok(org_report_text(&groups, widths, &lookups.browser_url))
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:450-517 (`DisplayOrgReport`)
// ```c++
//     for(size_t i = 0; i < m_BlastResTaxInfo->orderedTaxids.size(); i++) {
//         ...
//         orgHeader = x_MapTaxInfoTemplate(orgHeader,seqsForTaxID);
//         ...
//         orgHeader = CAlignFormatUtil::MapTemplate(orgHeader,"next_taxid",nextTaxid);
//         orgHeader = CAlignFormatUtil::MapTemplate(orgHeader,"disable_nexttaxid",hideNextTaxid);
//         orgHeader = CAlignFormatUtil::MapTemplate(orgHeader,"prev_taxid",prevTaxid);
//         orgHeader = CAlignFormatUtil::MapTemplate(orgHeader,"disable_prevtaxid",hidePrevTaxid);
//         orgHeader = CAlignFormatUtil::MapTemplate(orgHeader,"disable_top",hideTop);
//
//         string orgReportSeqInfo;
//         for(size_t j = 0; j < seqsForTaxID.seqInfoList.size(); j++) {
//             SSeqInfo *seqInfo = seqsForTaxID.seqInfoList[j];
//             string orgReportTableRow = x_MapSeqTemplate(m_TaxFormatTemplates->orgReportTableRow,seqInfo);
//             orgReportSeqInfo += orgReportTableRow;
//         }
//
//         string oneOrgInfo = orgHeader + orgReportSeqInfo;
//         orgReportData += oneOrgInfo;
//     }
//     orgReportData = CAlignFormatUtil::MapTemplate(m_TaxFormatTemplates->orgReportTable,"table_rows",orgReportData);
//     if(m_DisplayOption == eText) {
//         string orgReportTxtTableCaption = s_TagToConstString("ORG_REPORT_TXT_TABLE_CAPTION");
//         string reportCaption = CAlignFormatUtil::AddSpaces(orgReportTxtTableCaption,m_LineLength,CAlignFormatUtil::eSpacePosToCenter | CAlignFormatUtil::eAddEOLAtLineStart);
//         orgReportData = CAlignFormatUtil::MapTemplate(orgReportData,"org_report_caption",reportCaption);
//         string orgAccTxtTableHeader = s_TagToConstString("ORG_ACC_TXT_TABLE_HEADER");
//         orgReportData = CAlignFormatUtil::MapSpaceTemplate(orgReportData,"acc_hd",orgAccTxtTableHeader,m_MaxAccLength);
//         string orgDescrTxtTableHeader = s_TagToConstString("ORG_DESCR_TXT_TABLE_HEADER");
//         orgReportData = CAlignFormatUtil::MapSpaceTemplate(orgReportData,"descr_hd",orgDescrTxtTableHeader,m_MaxDescrLength,CAlignFormatUtil::eSpacePosToCenter);
//         string orgScoreTxtTableHeader = s_TagToConstString("ORG_SCORE_TXT_TABLE_HEADER");
//         orgReportData = CAlignFormatUtil::MapSpaceTemplate(orgReportData,"score_hd",orgScoreTxtTableHeader,m_MaxScoreLength);
//         string orgEValueTxtTableHeader = s_TagToConstString("ORG_EVALUE_TXT_TABLE_HEADER");
//         orgReportData = CAlignFormatUtil::MapSpaceTemplate(orgReportData,"evalue_hd",orgEValueTxtTableHeader,m_MaxEvalLength);
//     }
// ```
// `x_InitTaxIdToSeqsMap` fills a map that only the HTML report shows.
/// The text of the organism report of `groups` with the widths `widths`; `browser_url` replaces
/// `<@taxBrowserURL@>`.
fn org_report_text(groups: &[Group<'_, '_>], widths: Widths, browser_url: &[u8]) -> Text {
    const DISABLED: &[u8] = b"disabled=\"disabled\"";
    let mut data = Text::default();
    for (i, group) in groups.iter().enumerate() {
        let last = i + 1 == groups.len();
        let number = |index: usize| groups[index].taxid.to_string().into_bytes();
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:1172-1187
        // ```c++
        //     reportTableRow = CAlignFormatUtil::MapTemplate(reportTableRow,"bl_taxid",taxInfo.blNameTaxid);
        //     reportTableRow = CAlignFormatUtil::MapTemplate(reportTableRow,"taxid",TAX_ID_TO(int, taxInfo.taxid));
        //     reportTableRow = CAlignFormatUtil::MapTemplate(reportTableRow,"taxBrowserURL",m_TaxBrowserURL);
        //     reportTableRow = CAlignFormatUtil::MapTemplate(reportTableRow,"rid",m_Rid);
        //     int numHits = taxInfo.seqInfoList.size();
        //     numHits = (numHits > 0) ? numHits : taxInfo.numHits;
        //     reportTableRow = CAlignFormatUtil::MapTemplate(reportTableRow,"numhits",numHits);
        //
        //     string shift;
        //     for(size_t i = 0; i < depth; i++) {
        //         shift += ".";
        //     }
        //     reportTableRow = CAlignFormatUtil::MapTemplate(reportTableRow,"depth",shift);
        // ```
        // `<@bl_taxid@>` is refused before the report (`TaxLookups::new`); the line holds none.
        let header = organism_line(group.names)
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
// ```c++
// string CTaxFormat::x_MapSeqTemplate(string seqTemplate, SSeqInfo *seqInfo)
// {
//     string reportTableRow = CAlignFormatUtil::MapTemplate(seqTemplate,"gi",NStr::NumericToString(seqInfo->gi));
//     ...
//     reportTableRow = CAlignFormatUtil::MapTemplate(reportTableRow,"descr_abbr",seqInfo->title.substr(0,60));//for standalone output
//     reportTableRow = CAlignFormatUtil::MapTemplate(reportTableRow,"rid",m_Rid);
//
//     if(m_DisplayOption == eText) {
//         reportTableRow = CAlignFormatUtil::MapSpaceTemplate(reportTableRow,"acc",seqInfo->label,m_MaxAccLength);
//         reportTableRow = CAlignFormatUtil::MapSpaceTemplate(reportTableRow,"descr_text",seqInfo->title,m_MaxDescrLength);
//         reportTableRow = CAlignFormatUtil::MapSpaceTemplate(reportTableRow,"score",seqInfo->bit_score,m_MaxScoreLength);
//         reportTableRow = CAlignFormatUtil::MapSpaceTemplate(reportTableRow,"evalue",seqInfo->evalue,m_MaxEvalLength);
//     }
// ```
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
