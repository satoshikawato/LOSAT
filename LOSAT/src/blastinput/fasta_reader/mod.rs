//! NCBI BLAST+'s reading of FASTA input for BLASTN, TBLASTX, TBLASTN and BLASTP: the
//! input source (`CBlastFastaInputSource`, `CBlastInput`), the readers
//! (`CBlastInputReader`, `CCustomizedFastaReader`, `CFastaReader`) and the line reader
//! (`CStreamLineReader`). Nucleotide input uses `fAssumeNuc`, protein input `fAssumeProt`;
//! the other flags are the same.
//!
//! A record keeps what the search and the reports use: the generated local ID
//! (`Query_N`, `Subject_N`), the title bytes, the residues, and the reader's messages. The
//! lowercase mask (`x_OpenMask`/`x_CloseMask`) is kept as the case of the residues: a
//! residue in a masked interval is in lower case, every other residue in upper case.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:316-361
//! `-parse_deflines` (`fParseRawID`) is `believed`; `fSkipCheck` is not set by the BLAST programs.

pub mod believed;
pub(crate) mod reader;
pub(crate) mod seq_id;
pub(crate) mod stream;
#[cfg(test)]
mod tests;

use std::borrow::Cow;
use std::io::{Read, Seek};

pub(crate) use stream::FastaStream;
use stream::LineReader;

// The record and its public methods: `crate::record` (documented as `api::FastaRecord`).
pub use crate::record::FastaRecord;

/// How an input is read: the molecule, the role's local ID prefix and whether the data
/// loaders are configured (the Seq-id path of `CBlastInputReader`).
#[derive(Clone, Copy, Debug)]
pub struct ReaderConfig {
    /// `fAssumeProt` (protein input) or `fAssumeNuc`.
    pub protein: bool,
    /// `SetQueryLocalIdMode` (`Query_`) or `SetSubjectLocalIdMode` (`Subject_`).
    pub id_prefix: &'static str,
    /// `SDataLoaderConfig::UseDataLoaders()`: `CBlastInputReader` tries the first line of
    /// a record as a Seq-id.
    pub data_loaders: bool,
    /// The program and role named in LOSAT's rejections.
    pub program: &'static str,
    pub role: &'static str,
    /// Keep the defline line of each record (`FastaRecord::source`), for the BLAST archive.
    pub keep_defline: bool,
    /// `-parse_deflines` (`believed`).
    pub parse_deflines: bool,
    /// The strand of a query's Seq-loc (`CBlastInputSourceConfig::GetStrand`), read with
    /// `-parse_deflines`.
    pub strand: crate::report::archive::bioseq::SeqLocStrand,
    /// Whether the reader knows that the run writes the archive when it keeps the deflines
    /// (`keep_defline`), and refuses what the archive cannot write as it reads
    /// (`refuse_range_in_archive`). The library API reads before it knows the formats.
    pub archive_known: bool,
}

impl ReaderConfig {
    /// The same reader, keeping the defline lines (`-outfmt 11`).
    pub fn with_defline(mut self, keep: bool) -> Self {
        self.keep_defline = keep;
        self
    }

    // No NCBI counterpart: NCBI reads for one format; a caller of the library reads the records
    // before it asks for its formats (`api::read_fasta`), and the refusal comes with the search.
    /// The reader of a caller that keeps every defline, whatever the formats of its searches
    /// (`archive_known`).
    pub fn with_archive_unknown(mut self) -> Self {
        self.archive_known = false;
        self
    }

    /// The reader with `-parse_deflines` or without it (`believed`).
    pub fn with_parse_deflines(mut self, parse_deflines: bool) -> Self {
        self.parse_deflines = parse_deflines;
        self
    }

    /// The reader of a nucleotide query with `-strand value`.
    pub fn with_strand(mut self, value: &str) -> Self {
        self.strand = crate::report::archive::bioseq::SeqLocStrand::of_option(value);
        self
    }

    /// The query reader of `program`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:66-74
    pub fn query(program: &'static str, protein: bool, data_loaders: bool) -> Self {
        Self {
            protein,
            id_prefix: "Query_",
            data_loaders,
            program,
            role: "query",
            keep_defline: false,
            archive_known: true,
            parse_deflines: false,
            strand: if protein {
                crate::report::archive::bioseq::SeqLocStrand::Unknown
            } else {
                crate::report::archive::bioseq::SeqLocStrand::Both
            },
        }
    }

    /// The subject reader of `program`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:230-237
    pub fn subject(program: &'static str, protein: bool, data_loaders: bool) -> Self {
        Self {
            protein,
            id_prefix: "Subject_",
            data_loaders,
            program,
            role: "subject",
            keep_defline: false,
            archive_known: true,
            parse_deflines: false,
            strand: if protein {
                crate::report::archive::bioseq::SeqLocStrand::Unknown
            } else {
                crate::report::archive::bioseq::SeqLocStrand::Both
            },
        }
    }
}

/// The longest record that LOSAT reads, gaps included: NCBI keeps positions in `TSeqPos` (32
/// bits) and BLAST's sequence lengths in `Int4`.
pub const MAX_RECORD_LETTERS: u64 = i32::MAX as u64;

/// LOSAT's refusal of a record longer than [`MAX_RECORD_LETTERS`] letters, which the reader of
/// `config` gives as it reads one (`FastaReader::check_length`) and LOSAT Web's index scan as it
/// scans one (`web/adapter/src/scan/ncbi.rs`, which cannot read such a record into the reader).
pub fn record_too_long(config: &ReaderConfig) -> ReadError {
    ReadError::Unsupported(anyhow::anyhow!(
        "a {} record longer than {MAX_RECORD_LETTERS} letters is not supported by LOSAT's {}",
        config.role,
        config.program
    ))
}

/// The letters of the gap line `line` (`>?` and what follows on its line, after NCBI's
/// `TruncateSpaces` and its change of `>?_` into `>`): the run of `N` or `X` that the reader puts
/// in the record for it (`FastaReader::parse_gap_line`, `reader::gap_line_size`), without the
/// line's warnings; `None` for a line that does not start with `>?`. LOSAT Web's index scan
/// counts a gap line's letters with it, so that its lengths are the reader's
/// (`web/adapter/src/scan/ncbi.rs`).
pub fn gap_line_letters(line: &[u8]) -> Option<u32> {
    if !line.starts_with(b">?") {
        return None;
    }
    // Only the warnings can fail, and these are dropped.
    reader::gap_line_size(line, 0, &mut |_| Ok(())).ok()
}

/// The parts of a record that only the BLAST archive writes.
///
/// NCBI keeps the defline line as read in the Bioseq's `CFastaReader` user object, and the
/// whole sequence in the Bioseq of a query that `-query_loc` cuts.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct RecordSource {
    /// The defline line as read, with its `>`, without its line end (empty when not kept).
    pub defline: Vec<u8>,
    /// The record's letters before a cut to a range (`FastaRecord::cut`).
    pub uncut: Option<Vec<u8>>,
    /// The Seq-ids `-parse_deflines` read (`believed`).
    pub believed: Option<std::sync::Arc<believed::BelievedIds>>,
    /// The start of the letters after a cut to a range, in the record.
    pub cut_from: usize,
    /// The length of the record before a cut to a range (0 for a record that no range cut);
    /// BLASTX's reports read where the letters lie in the record from it and `cut_from`. A plain
    /// `usize` keeps the source in glibc's 80-byte chunks.
    pub uncut_length: usize,
}

impl RecordSource {
    /// The source of a part of the record (a chunk of a split query): all of it but the whole
    /// letters, which only the archive's Bioseq of the record reads.
    pub fn without_uncut(&self) -> Box<Self> {
        Box::new(Self {
            defline: self.defline.clone(),
            uncut: None,
            believed: self.believed.clone(),
            cut_from: self.cut_from,
            uncut_length: self.uncut_length,
        })
    }
}

// What the engine reads of a record: outside the library API (`api::FastaRecord`), whose types
// these would bring along.
#[doc(hidden)]
impl FastaRecord {
    /// The Seq-ids `-parse_deflines` read, `None` without the option.
    pub fn believed(&self) -> Option<&believed::BelievedIds> {
        self.source.as_ref()?.believed.as_deref()
    }

    /// The record of the Bioseq whose Seq-ids are `ids` (at least one), shown as `-parse_deflines`
    /// shows the ids it read, whose Seq-aligns take the first (`BelievedIds::of_bioseq`): the query
    /// of psiblast's and tblastn's `-in_pssm`, the PSSM's Bioseq (`pssm_search::PsiSearch::
    /// query_record`; NCBI BLAST+ 2.17.0: `Query= gi|4504347 <title>`, `qseqid` `gi|4504347`, and
    /// with `-parse_deflines` `Query= Query_1 <title>` of a local id).
    pub fn believing(mut self, ids: Vec<crate::seqid::ParsedSeqId>) -> Self {
        let believed = believed::BelievedIds::of_bioseq(
            ids,
            crate::report::archive::bioseq::SeqLocStrand::Unknown,
        );
        self.local_id = believed.bare_id();
        self.source.get_or_insert_with(Default::default).believed =
            Some(std::sync::Arc::new(believed));
        self
    }

    /// The record's names (`RecordNames`): its local ID, its title and the Seq-ids
    /// `-parse_deflines` read.
    pub fn names(&self) -> RecordNames<'_> {
        RecordNames {
            local_id: &self.local_id,
            title: &self.title,
            believed: self.believed(),
        }
    }

    /// What `Query=` and `# Query:` show (`RecordNames::query_line`).
    pub fn query_line(&self) -> Cow<'_, [u8]> {
        self.names().query_line()
    }

    /// The id of the record in NCBI's warnings (`RecordNames::warning_id`).
    pub fn warning_id(&self) -> &[u8] {
        self.names().warning_id()
    }

    /// The Seq-id of the record in the Seq-aligns and the archive (`RecordNames::seq_id`).
    pub fn seq_id(&self) -> crate::report::seqid::SeqId<'_> {
        self.names().seq_id()
    }

    /// The record as the Seq-loc of the archive's `seq-loc-list` (`-parse_deflines`, an id that is
    /// not local): the searched interval (`x_FastaToSeqLoc`; the `to` of an empty record is -1).
    pub fn archive_seq_loc(&self) -> Option<crate::report::archive::bioseq::QuerySeqLoc<'_>> {
        let source = self.source.as_ref()?;
        let from = source.cut_from as i32;
        source
            .believed
            .as_ref()?
            .seq_loc(from, from + self.sequence.len() as i32 - 1)
    }

    /// The ID that the tabular formats and the reports show (`RecordNames::shown_id`).
    pub fn shown_id(&self) -> &[u8] {
        self.names().shown_id()
    }

    /// The subject ID of the tabular formats of a `-subject` record that the default report
    /// prints (`RecordNames::subject_tabular_id`).
    pub fn subject_tabular_id(&self, protein: bool) -> Result<Option<Vec<u8>>, &'static str> {
        self.names().subject_tabular_id(protein)
    }

    /// The record cut to the residues `from..to_exclusive` of a range (`-query_loc`,
    /// `-subject_loc`): the same local ID and title, the letters (their case too) of the
    /// interval, and no messages (the reader wrote them when it read the whole record).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:458-466
    /// NCBI keeps the whole record and searches the interval of its location; LOSAT keeps
    /// the interval's letters and where it lies (`seq_range::Placements`).
    pub fn cut(&self, from: usize, to_exclusive: usize) -> Self {
        Self {
            local_id: self.local_id.clone(),
            title: self.title.clone(),
            sequence: self.sequence[from..to_exclusive].to_vec(),
            warnings: Vec::new(),
            source: self.source.as_ref().map(|source| {
                Box::new(RecordSource {
                    defline: source.defline.clone(),
                    uncut: Some(
                        source
                            .uncut
                            .clone()
                            .unwrap_or_else(|| self.sequence.clone()),
                    ),
                    believed: source.believed.clone(),
                    cut_from: source.cut_from + from,
                    uncut_length: source.uncut_length.max(self.sequence.len()),
                })
            }),
        }
    }

    /// The defline line as read (with its `>`), when the reader kept it and the record has one.
    pub fn archive_defline(&self) -> Option<&[u8]> {
        self.source
            .as_ref()
            .map(|source| source.defline.as_slice())
            .filter(|defline| !defline.is_empty())
    }

    /// The letters of the record as the archive's Bioseq holds them: the whole record
    /// even when a range cut it.
    pub fn archive_letters(&self) -> &[u8] {
        self.source
            .as_ref()
            .and_then(|source| source.uncut.as_deref())
            .unwrap_or(&self.sequence)
    }

    /// The record `cut` makes, in place: the letters before `from` and from `to_exclusive` on are
    /// dropped (the room of the whole letters stays until the sequence is replaced). A record that
    /// keeps its source for outfmt 11 keeps the whole letters there first, as `cut` does.
    pub fn cut_in_place(&mut self, from: usize, to_exclusive: usize) {
        if let Some(source) = self.source.as_mut() {
            source.uncut.get_or_insert_with(|| self.sequence.clone());
            source.uncut_length = source.uncut_length.max(self.sequence.len());
            source.cut_from += from;
        }
        self.sequence.truncate(to_exclusive);
        self.sequence.drain(..from);
        self.warnings = Vec::new();
    }
}

/// The names of a record that the reports show, made from its local ID, its title bytes and the
/// Seq-ids `-parse_deflines` read: the one owner of these rules, for the shared reader's records
/// (`FastaRecord::names`) and BLASTX's (`algorithm::blastx::input::FastaRecord::names`).
#[derive(Clone, Copy, Debug)]
pub struct RecordNames<'a> {
    /// The generated local ID (`Query_N`, `Subject_N`), or with `-parse_deflines` the bare id
    /// (`BelievedIds::bare_id`).
    pub local_id: &'a str,
    /// The title as read, bytes (also bytes that are not UTF-8).
    pub title: &'a [u8],
    /// The Seq-ids `-parse_deflines` read, `None` without the option.
    pub believed: Option<&'a believed::BelievedIds>,
}

impl<'a> RecordNames<'a> {
    /// What `Query=` and `# Query:` show: the title, and with `-parse_deflines` the id before it
    /// (`local_id`, `BelievedIds::bare_id`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:732-734
    pub fn query_line(self) -> Cow<'a, [u8]> {
        if self.believed.is_none() {
            return Cow::Borrowed(self.title);
        }
        let mut line = self.local_id.as_bytes().to_vec();
        if !self.title.is_empty() {
            line.push(b' ');
            line.extend_from_slice(self.title);
        }
        Cow::Owned(line)
    }

    /// The id of the record in NCBI's warnings (`report::query_warnings::warning_query_id`).
    pub fn warning_id(self) -> &'a [u8] {
        match self.believed {
            Some(ids) => &ids.warning_id,
            None => self.local_id.as_bytes(),
        }
    }

    /// The query ID of NCBI's warnings of a query (`report::query_warnings::warning_query_id`):
    /// the id (`warning_id`), then a space and the title's bytes, cut to 25 bytes and `.. ` when
    /// longer than 35.
    pub fn warning_query_id(self) -> Vec<u8> {
        crate::report::query_warnings::warning_query_id(self.warning_id(), self.title)
    }

    /// The Seq-id of the record in the Seq-aligns and the archive.
    pub fn seq_id(self) -> crate::report::seqid::SeqId<'a> {
        match self.believed {
            Some(ids) => crate::report::seqid::SeqId::Db(&ids.best),
            None => crate::report::seqid::SeqId::Local(self.local_id),
        }
    }

    /// The ID that the tabular formats and the reports show for a record with a local
    /// ID: the title up to its first space, or the local ID when there is no title. With
    /// `-parse_deflines` the `qaccver` of its ids.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:474-504
    /// (`s_ReplaceLocalId`; see `report::` for its callers)
    pub fn shown_id(self) -> &'a [u8] {
        if let Some(ids) = self.believed {
            return &ids.tabular[2];
        }
        if self.title.is_empty() {
            self.local_id.as_bytes()
        } else {
            let end = self
                .title
                .iter()
                .position(|&byte| byte == b' ')
                .unwrap_or(self.title.len());
            &self.title[..end]
        }
    }

    /// The subject ID of the tabular formats of a `-subject` record that the default report
    /// prints (`saccver`; also `sseqid` and `sacc` of a record read without `-parse_deflines`):
    /// `report::defline::tabular_subject_id`, or with `-parse_deflines` the `saccver` of the
    /// record's ids as `report::db_subject::ParsedSubject` makes them (NCBI keeps a parsed local
    /// id, `s_ReplaceLocalId` with `parse_local`, and `GetSeqIdList` replaces an `lcl|Subject_`
    /// one with the first word of the title generated from the rest of the defline).
    /// `Ok(None)` where `GetSeqIdList` throws (`GenerateDefline` cannot decode the title,
    /// `report::tabular::SubjectLookup::Failed`); `Err` where LOSAT refuses the ids of a record
    /// read with `-parse_deflines` (`ParsedSubject::tabular_ids`; never a nucleotide record's).
    pub fn subject_tabular_id(self, protein: bool) -> Result<Option<Vec<u8>>, &'static str> {
        match self.believed {
            Some(ids) => Ok(
                crate::report::db_subject::ParsedSubject::of_record(ids, self.title)
                    .tabular_ids(protein)?
                    .map(|[_, _, saccver]| saccver),
            ),
            None => Ok(crate::report::defline::tabular_subject_id(
                self.local_id.as_bytes(),
                self.title,
                protein,
            )
            .ok()),
        }
    }
}

/// A record as the programs use it: its residues, the record cut to a range, and its title
/// bytes (the ranges in `seq_range.rs`, the query warnings, the masks and the lookup
/// tables).
pub trait InputRecord: Sized {
    /// The residues (the search's letters, lower case where the lowercase mask covers
    /// them).
    fn seq(&self) -> &[u8];

    /// The record cut to the residues `from..to_exclusive` (`FastaRecord::cut`).
    fn cut(&self, from: usize, to_exclusive: usize) -> Self;

    /// The title (`CBlastQuerySourceOM::GetTitle`: the record's title descriptor; empty
    /// without one).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:366-372
    fn title_bytes(&self) -> Cow<'_, [u8]>;

    /// The id of the record in NCBI's warnings when `-parse_deflines` read it.
    fn believed_warning_id(&self) -> Option<&[u8]>;
}

impl InputRecord for FastaRecord {
    fn seq(&self) -> &[u8] {
        &self.sequence
    }

    fn cut(&self, from: usize, to_exclusive: usize) -> Self {
        FastaRecord::cut(self, from, to_exclusive)
    }

    fn title_bytes(&self) -> Cow<'_, [u8]> {
        Cow::Borrowed(&self.title)
    }

    fn believed_warning_id(&self) -> Option<&[u8]> {
        self.believed().map(|ids| ids.warning_id.as_slice())
    }
}

/// `CObjReaderParseException`'s codes that the reader throws with these flags.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ParseErrorCode {
    /// `eEOF`: the batch ends (`GetNextSeqBatch`, `GetAllSeqs`).
    Eof,
    /// `eFormat` (`CheckDataLine`).
    Format,
    /// `eNoDefline`: the first line that is neither blank nor a comment is not a defline, under
    /// `-parse_deflines` (a reader without `fDLOptional`). The queries end with it as with any
    /// reader error; the subjects end with `GetAllSeqs`'s own error (`ReadError::into_all_seqs_error`).
    NoDefline,
}

/// Why a record was not read.
#[derive(Debug)]
pub enum ReadError {
    /// A `CObjReaderParseException` with its message (`GetMsg()`) and line.
    Parse {
        code: ParseErrorCode,
        message: String,
        line: u64,
    },
    /// Input that LOSAT rejects explicitly (`not supported by LOSAT's <PROGRAM>`).
    Unsupported(anyhow::Error),
    /// The warnings could not be written.
    Write(std::io::Error),
}

impl std::fmt::Display for ReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReadError::Parse { message, .. } => f.write_str(message),
            ReadError::Unsupported(error) => write!(f, "{error:#}"),
            ReadError::Write(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for ReadError {}

impl ReadError {
    /// The error as the programs end with it: a reader exception is NCBI's `BLAST query
    /// error: <message>` with exit 1 (`blastinput::app::query_error`, the one producer of that
    /// line), a rejection is LOSAT's, and a failed write of the messages is the write error.
    pub fn into_app_error(self) -> anyhow::Error {
        match self {
            ReadError::Parse { message, .. } => crate::blastinput::app::query_error(&message),
            ReadError::Unsupported(error) => error,
            ReadError::Write(error) => error.into(),
        }
    }

    /// The error as the programs end with it when it stops a read of all records
    /// (`CBlastInput::GetAllSeqs`, the subjects): `eNoDefline` is `GetAllSeqs`'s own input
    /// error (`no_defline_error`), every other error is `into_app_error`'s.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:205-214
    pub fn into_all_seqs_error(self) -> anyhow::Error {
        match self {
            ReadError::Parse {
                code: ParseErrorCode::NoDefline,
                line,
                ..
            } => no_defline_error(line),
            error => error.into_app_error(),
        }
    }
}

/// The error of `CBlastInput::GetAllSeqs` (the subjects) whose reader stopped at line `line`
/// with `eNoDefline`: a `CInputException` (`BLAST query/options error: ...`, exit 1). The text
/// ends with the NUL that `ends` writes into the stream; the error archive drops it
/// (`report::diag_log::sanitize`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:209-214
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-176
pub fn no_defline_error(line: u64) -> anyhow::Error {
    crate::blastinput::app::options_error(&format!(
        "Query input doesn't start with a defline or comment, line {line}\0"
    ))
}

/// NCBI's `CBlastFastaInputSource` over one opened input.
pub struct FastaInputSource<R: Read> {
    reader: reader::FastaReader<R>,
}

impl FastaInputSource<crate::blastinput::input_files::InputFile> {
    /// An input source over an opened file or standard input (`CStreamLineReader` over
    /// the argument's stream, which reads as `stream::InputBackend::of` says).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:291-300
    /// `GetConvertGapsToNs` is set only by the mapper (magicblast).
    pub fn from_file(
        file: impl Into<crate::blastinput::input_files::InputFile>,
        config: ReaderConfig,
    ) -> Self {
        Self::from_stream(FastaStream::from_file(file), config)
    }
}

impl<'a> FastaInputSource<&'a [u8]> {
    /// An input source over bytes in memory (the ABI's registered inputs), read as a
    /// regular file with the same bytes (`FastaStream::from_bytes`; test
    /// `bytes_read_as_a_file_with_the_same_bytes`).
    pub fn from_bytes(bytes: &'a [u8], config: ReaderConfig) -> Self {
        Self::from_stream(FastaStream::from_bytes(bytes), config)
    }
}

#[cfg(test)]
impl<'a> FastaInputSource<&'a [u8]> {
    /// `from_bytes` with LOSAT's bulk paths off: one byte at a time, the reference that
    /// the tests compare the bulk paths with.
    pub(crate) fn from_bytes_byte_at_a_time(bytes: &'a [u8], config: ReaderConfig) -> Self {
        let mut stream = FastaStream::from_bytes(bytes);
        stream.bulk = false;
        Self::from_stream(stream, config)
    }
}

impl<R: Read + Seek> FastaInputSource<R> {
    /// `IsIStreamEmpty` on the input, which the programs call on the query's stream before
    /// they build the input source: true when the input has only white space (`Query is
    /// Empty!`), false for a stream without a position, such as a pipe, whatever it holds.
    /// Call it before the first record is read.
    ///
    /// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:845-875
    /// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:209-214
    /// LOSAT takes the Linux branch on every platform: the oracle is NCBI BLAST+ on Linux
    /// (`AUTHORITY.md` §C1).
    pub fn stream_is_empty(&mut self) -> bool {
        self.reader.lines.stream_is_empty()
    }
}

impl<R: Read> FastaInputSource<R> {
    fn from_stream(stream: FastaStream<R>, config: ReaderConfig) -> Self {
        Self {
            reader: reader::FastaReader::new(LineReader::new(stream), config),
        }
    }

    /// LOSAT_X_QSTREAM's first reading in parts: a reader over `input`, a part of a query file that
    /// starts at a record start (`x_query_scan_par.rs`); its stream is a file's without the
    /// readable-byte hint, which only the pushback of a `\r` reads.
    #[cfg_attr(not(all(unix, feature = "parallel")), allow(dead_code))]
    pub(crate) fn x_from_reader(input: R, config: ReaderConfig) -> Self {
        Self::from_stream(FastaStream::new(input), config)
    }

    /// LOSAT_X_QSTREAM's second reading in parts: `x_from_reader` for a part whose first record the
    /// file's reader reads after `lines` lines, with the local ID `next_id`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:48-56
    /// The BLAST programs leave `local_id_counter` at 1 (`FastaReader::new`); a part's reader starts
    /// it where the file's one reader is (`FastaReader::x_start_at`).
    #[cfg_attr(not(all(unix, feature = "parallel")), allow(dead_code))]
    pub(crate) fn x_from_reader_at(
        input: R,
        config: ReaderConfig,
        lines: u64,
        next_id: i64,
    ) -> Self {
        let mut source = Self::x_from_reader(input, config);
        source.reader.x_start_at(lines, next_id);
        source
    }

    /// LOSAT_X_QSTREAM's reading in parts: the lines read so far (`GetLineNumber`,
    /// line_reader.cpp:208-216, `LineReader::line_number`).
    #[cfg_attr(not(all(unix, feature = "parallel")), allow(dead_code))]
    pub(crate) fn x_line_number(&self) -> u64 {
        self.reader.lines.line_number()
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:370-374
    pub fn end(&mut self) -> bool {
        self.reader.lines.at_eof()
    }

    /// The search whose records this source reads under `-parse_deflines`: a query and a subject
    /// with one id are compared (`crate::seqid::id_scope::IdScope`).
    pub fn set_id_scope(
        &mut self,
        id_scope: Option<std::sync::Arc<crate::seqid::id_scope::IdScope>>,
    ) {
        self.reader.set_id_scope(id_scope);
    }

    /// The next record (`GetNextSequence`, without the range and strand of
    /// `x_FastaToSeqLoc`, which the programs apply: `seq_range.rs`). The reader's messages
    /// go to `warn` as NCBI writes them, and are also kept in the record.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:376-387
    /// The molecule checks after it (`Nucleotide FASTA provided for protein sequence` and
    /// the reverse) cannot fail: `AssignMolType` sets the molecule of the input.
    ///
    /// A line that NCBI's data loaders would fetch as a Seq-id gives
    /// `ReadError::Unsupported`; the source is then past that one line, with the local-ID
    /// counter unchanged, so a further call reads from the next line, as NCBI's next
    /// `ReadOneSeq` does. The programs stop at the rejection (`read_all`, `read_queries`),
    /// where NCBI skips a query it cannot fetch and fails on such a subject.
    /// EXPERIMENT (LOSAT_X_BNQPACK): pack the letters of the nucleotide records read from now on;
    /// a record whose letters packed has an empty `sequence`, and `x_take_packed` gives them.
    pub fn x_pack_letters(&mut self) {
        self.reader.x_pack_letters();
    }

    /// EXPERIMENT (LOSAT_X_BNQPACK): the packed letters of the record read last, if it was packed.
    pub fn x_take_packed(&mut self) -> Option<crate::sequence::packed_letters::PackedLetters> {
        self.reader.x_take_packed()
    }

    pub fn next_sequence(
        &mut self,
        warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
    ) -> Result<FastaRecord, ReadError> {
        let mut warnings = Vec::new();
        let mut record = self.reader.read_one_seq(&mut |message: &[u8]| {
            warnings.extend_from_slice(message);
            warn(message)
        })?;
        record.warnings = warnings;
        Ok(record)
    }
}

/// All records of an input (`CBlastInput::GetAllSeqs`, the subjects): the messages go to
/// `warn` as they are written; a reader error stops the reading.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:198-219
/// `eNoDefline` needs a reader without `fDLOptional`; the error is returned as the reader gave
/// it, and a caller that reads subjects ends with `ReadError::into_all_seqs_error`.
pub fn read_all<R: Read>(
    source: &mut FastaInputSource<R>,
    warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
) -> Result<Vec<FastaRecord>, ReadError> {
    let mut records = Vec::new();
    while !source.end() {
        match source.next_sequence(warn) {
            Ok(record) => records.push(record),
            Err(ReadError::Parse {
                code: ParseErrorCode::Eof,
                ..
            }) => break,
            Err(error) => return Err(error),
        }
    }
    Ok(records)
}

/// The subjects of `-subject` as NCBI reads them while it processes the arguments
/// (`ReadSequencesToBlast`): every record, with the reader's messages written to `warn` as
/// each record is read; a record whose interval of `range` (`-subject_loc`) starts more than
/// one letter past its end stops the reading after that record, with NCBI's range error. The
/// records are returned whole (`seq_range::cut_subjects` cuts them to `range`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:221-247
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:198-219
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:450-453
/// `eNoDefline` needs a reader without `fDLOptional` (`-parse_deflines`): it is `GetAllSeqs`'s
/// `CInputException` (`no_defline_error`). The range's `CInputException` is also
/// `BLAST query/options error: ...` (`seq_range::record_interval`, `app::options_error`); any other
/// reader exception is `BLAST query error: ...` (`ReadError::into_all_seqs_error`). The source's
/// `ReaderConfig::subject` carries `read_proteins` and the data loaders of `dlconfig`;
/// `gaps_to_Ns` is the mapper's; LOSAT keeps the lowercase mask in every record
/// (`FastaRecord::sequence`).
pub fn read_subjects<R: Read>(
    source: &mut FastaInputSource<R>,
    range: Option<&crate::blastinput::seq_range::SequenceRange>,
    warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
) -> anyhow::Result<Vec<FastaRecord>> {
    read_subjects_map(source, range, warn, Ok)
}

/// `read_subjects`, with each record given to `convert` as soon as it has been read and its
/// range checked; the results are kept in record order. A caller that keeps a smaller form of a
/// record (`sequence::packed_subject`) never holds the letters of the whole set at once. An
/// error of `convert` ends the reading with that error.
pub fn read_subjects_map<R: Read, T>(
    source: &mut FastaInputSource<R>,
    range: Option<&crate::blastinput::seq_range::SequenceRange>,
    warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
    mut convert: impl FnMut(FastaRecord) -> anyhow::Result<T>,
) -> anyhow::Result<Vec<T>> {
    use crate::blastinput::seq_range::{record_interval, RecordInterval};
    let mut records = Vec::new();
    while !source.end() {
        match source.next_sequence(warn) {
            Ok(record) => {
                if let Some(range) = range {
                    if record_interval(range, record.sequence.len()) == RecordInterval::PastEnd {
                        return Err(crate::blastinput::app::options_error(
                            "Invalid from coordinate (greater than sequence length)",
                        ));
                    }
                }
                records.push(convert(record)?);
            }
            Err(ReadError::Parse {
                code: ParseErrorCode::Eof,
                ..
            }) => break,
            Err(error) => return Err(error.into_all_seqs_error()),
        }
    }
    Ok(records)
}

/// The queries of an input read to its end, as NCBI's batches read them one record at a
/// time (`CBlastInput::GetNextSeqBatch`): the records, and how the reading ended. The
/// caller cuts the batches (`query_batch.rs`) and writes each record's messages when its
/// batch is read.
#[derive(Debug)]
pub struct QueryRecords {
    pub records: Vec<FastaRecord>,
    pub end: QueryEnd,
}

/// How the reading of the queries ended.
#[derive(Debug)]
pub enum QueryEnd {
    /// `End()`: the line reader is at the end of the input.
    Input,
    /// `eEOF` from the reader (only blank and comment lines were left; no record follows):
    /// the batch being read ends there.
    EofError,
    /// Another reader error, thrown while reading the record after `records`, with the
    /// messages written before it.
    Error { error: ReadError, warnings: Vec<u8> },
}

/// Reads all queries (`GetNextSeqBatch` until `End()`): the messages are kept with the
/// records (written later, batch by batch, by the caller).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:134-171
/// The `catch (const exception&)` skips a query whose Seq-id was not found or whose
/// range does not fit (`seq_range.rs`); LOSAT rejects Seq-id lines.
pub fn read_queries<R: Read>(source: &mut FastaInputSource<R>) -> QueryRecords {
    read_queries_with(source, |_, _| {})
}

/// EXPERIMENT (LOSAT_X_BNQPACK): `read_queries` with the letters of each record packed as it is
/// read (`FastaInputSource::x_pack_letters`), so the letters of the whole input are never held:
/// the records (a packed one without letters) and the packed letters of the packed ones, with
/// their positions among the records, in order.
pub fn x_read_queries_packed<R: Read>(
    source: &mut FastaInputSource<R>,
) -> (
    QueryRecords,
    Vec<(usize, crate::sequence::packed_letters::PackedLetters)>,
) {
    source.x_pack_letters();
    let mut packed = Vec::new();
    let records = read_queries_with(source, |source, index| {
        if let Some(letters) = source.x_take_packed() {
            packed.push((index, letters));
        }
    });
    packed.shrink_to_fit();
    (records, packed)
}

/// `read_queries`, with `read` called after each record is read (with its position), before the
/// next one.
fn read_queries_with<R: Read>(
    source: &mut FastaInputSource<R>,
    mut read: impl FnMut(&mut FastaInputSource<R>, usize),
) -> QueryRecords {
    let mut records = Vec::new();
    loop {
        if source.end() {
            return QueryRecords {
                records,
                end: QueryEnd::Input,
            };
        }
        let mut warnings = Vec::new();
        match source.next_sequence(&mut |message: &[u8]| {
            warnings.extend_from_slice(message);
            Ok(())
        }) {
            Ok(record) => {
                records.push(record);
                read(source, records.len() - 1);
            }
            Err(ReadError::Parse {
                code: ParseErrorCode::Eof,
                ..
            }) => {
                return QueryRecords {
                    records,
                    end: QueryEnd::EofError,
                }
            }
            Err(error) => {
                return QueryRecords {
                    records,
                    end: QueryEnd::Error { error, warnings },
                }
            }
        }
    }
}
