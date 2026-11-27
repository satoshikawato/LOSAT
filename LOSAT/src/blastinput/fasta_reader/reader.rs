//! NCBI's `CFastaReader` with the flags of `CBlastFastaInputSource`
//! (`fNoParseID` and `fDLOptional` or, with `-parse_deflines`, `fParseRawID`, `fAssumeNuc` or
//! `fAssumeProt`, `fNoSplit`,
//! `fHyphensIgnoreAndWarn`, `fDisableNoResidues`, `fQuickIDCheck`), as overridden by
//! `CCustomizedFastaReader` and `CBlastInputReader`.

use std::io::Read;

use super::stream::LineReader;
use super::{FastaRecord, ReadError, ReaderConfig, UpperRecord};
use crate::ncbi_str::{
    is_blank, is_space, truncate_spaces, truncate_spaces_begin, truncate_spaces_end,
};

/// The kinds of reader problems (`ILineError::EProblem`) that LOSAT meets, for
/// `PostWarning`'s list of ignored problems.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:358-360
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Problem {
    ModifierFoundButNoneExpected,
    TooLong,
    TooManyAmbiguousResidues,
    UnexpectedNucResidues,
    UnexpectedAminoAcids,
    IgnoredResidue,
    InvalidResidue,
    NonPositiveLength,
    ParsingModifiers,
    UnrecognizedQualifierName,
    ContradictoryModifiers,
    ExtraModifierFound,
    ExpectedModifierMissing,
}

/// The problems that `CBlastFastaInputSource::x_InitInputReader` tells the reader to ignore.
const IGNORED_PROBLEMS: [Problem; 3] = [
    Problem::ModifierFoundButNoneExpected,
    Problem::TooLong,
    Problem::TooManyAmbiguousResidues,
];

/// `CSeq_gap::EType` of the gap types that `ParseGapLine` reads (`Seq_gap.cpp:180-191`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum GapType {
    Unknown,
    ShortArm,
    Heterochromatin,
    Centromere,
    Telomere,
    Repeat,
    Contig,
    Scaffold,
    Contamination,
}

/// `CSeq_gap::ELinkEvid`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum LinkEvid {
    UnspecifiedOnly,
    Forbidden,
    Required,
}

/// NCBI reference (598d8ae6): c++/src/objects/seq/Seq_gap.cpp:176-194
/// `NameToGapTypeInfo` looks the name up after `CanonicalizeString` as a C string: the map
/// compares `const char*` keys with `strcmp` (`PCase_CStr`), so a value is looked up up to
/// its first NUL byte and the bytes after it are ignored (audit finding A-1 of session SFc).
/// The warning for an unknown value still writes the whole value (`gap_mod_warnings`).
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/Seq_gap.cpp:157-172
/// NCBI reference (598d8ae6): c++/include/objects/seq/Seq_gap.hpp:90-93
/// `Seq_gap.cpp`'s own `CanonicalizeString` (:51-70) maps the bytes as `fasta.cpp`'s
/// (`canonicalize_string`). The other lookups of a gap line compare whole byte strings:
/// the keys with `std::string::operator==` and the linkage evidence in a
/// `map<CTempString, ..., PQuickStringLess>` (serial/enumvalues.hpp:55).
pub(crate) fn gap_type_info(name: &[u8]) -> Option<(GapType, LinkEvid)> {
    let canonical = canonicalize_string(name);
    let c_string = canonical
        .iter()
        .position(|&byte| byte == 0)
        .map_or(&canonical[..], |nul| &canonical[..nul]);
    Some(match c_string {
        b"between-scaffolds" => (GapType::Contig, LinkEvid::Required),
        b"centromere" => (GapType::Centromere, LinkEvid::Forbidden),
        b"contamination" => (GapType::Contamination, LinkEvid::Required),
        b"heterochromatin" => (GapType::Heterochromatin, LinkEvid::Forbidden),
        b"repeat-between-scaffolds" => (GapType::Repeat, LinkEvid::Forbidden),
        b"repeat-within-scaffold" => (GapType::Repeat, LinkEvid::Required),
        b"short-arm" => (GapType::ShortArm, LinkEvid::Forbidden),
        b"telomere" => (GapType::Telomere, LinkEvid::Forbidden),
        b"unknown" => (GapType::Unknown, LinkEvid::UnspecifiedOnly),
        b"within-scaffold" => (GapType::Scaffold, LinkEvid::Forbidden),
        _ => return None,
    })
}

/// The values of `CLinkage_evidence::EType` by their ASN.1 names
/// (`NameToValue` of the enum's info, from `ENUM_METHOD_NAME`).
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:383-397
pub(crate) fn linkage_evidence_value(name: &[u8]) -> Option<u8> {
    Some(match name {
        b"paired-ends" => 0,
        b"align-genus" => 1,
        b"align-xgenus" => 2,
        b"align-trnscpt" => 3,
        b"within-clone" => 4,
        b"clone-contig" => 5,
        b"map" => 6,
        b"strobe" => 7,
        b"unspecified" => 8,
        b"pcr" => 9,
        b"proximity-ligation" => 10,
        b"other" => 255,
        _ => return None,
    })
}

/// `CLinkage_evidence::eType_unspecified`.
pub(crate) const LINKAGE_UNSPECIFIED: u8 = 8;

/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2128-2146
pub(crate) fn canonicalize_string(value: &[u8]) -> Vec<u8> {
    value
        .iter()
        .map(|&ch| match ch {
            b'A'..=b'Z' => ch.to_ascii_lowercase(),
            b' ' | b'_' => b'-',
            _ => ch,
        })
        .collect()
}

/// One `>?` line of a record (`CFastaReader::SGap`): where it is (raw residues before it)
/// and its length.
struct Gap {
    position: u32,
    length: u32,
}

/// The title of a defline (`SLineTextAndLoc`; its line number is not read).
struct LineText {
    text: Vec<u8>,
}

/// What a residue byte of a data line is (`ParseDataLine`'s `ECharType`, with the
/// residues and the white space that the first switch handles itself).
#[derive(Clone, Copy, PartialEq, Eq)]
enum CharType {
    Residue,
    MaskedResidue,
    HyphenToIgnoreAndWarn,
    Comment,
    Space,
    Bad,
}

/// The character types of a data line for nucleotide (`bIsNuc`) or protein
/// parsing, without `fParseGaps` (no gap characters).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:809-814
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:854-935
/// `bIsNuc` comes from the molecule set by `xSetSeqMol` under `-parse_deflines`, else
/// `fAssumeNuc`; `bAllowLetterGaps` and `bHyphensAreGaps` are false (no `fParseGaps`),
/// `bHyphensIgnoreAndWarn` is true. LOSAT looks the type up in a table (speed only).
const fn char_types(protein: bool) -> [CharType; 256] {
    let mut table = [CharType::Bad; 256];
    let nucleotide = b"ABCDGHKMRSTUVWYN";
    let mut index = 0;
    while index < nucleotide.len() {
        table[nucleotide[index] as usize] = CharType::Residue;
        table[nucleotide[index].to_ascii_lowercase() as usize] = CharType::MaskedResidue;
        index += 1;
    }
    if protein {
        let protein_only = b"EFIJLOPQZX";
        let mut index = 0;
        while index < protein_only.len() {
            table[protein_only[index] as usize] = CharType::Residue;
            table[protein_only[index].to_ascii_lowercase() as usize] = CharType::MaskedResidue;
            index += 1;
        }
        table[b'*' as usize] = CharType::Residue;
    }
    table[b'-' as usize] = CharType::HyphenToIgnoreAndWarn;
    table[b';' as usize] = CharType::Comment;
    let spaces = b"\t\n\x0b\x0c\r ";
    let mut index = 0;
    while index < spaces.len() {
        table[spaces[index] as usize] = CharType::Space;
        index += 1;
    }
    table
}

const NUCLEOTIDE_CHAR_TYPES: [CharType; 256] = char_types(false);
const PROTEIN_CHAR_TYPES: [CharType; 256] = char_types(true);

/// 1 for the bytes of `types` that are `kind`, 0 for the others.
const fn kind_table(types: &[CharType; 256], kind: CharType) -> [u8; 256] {
    let mut table = [0; 256];
    let mut index = 0;
    while index < 256 {
        table[index] = (types[index] as u8 == kind as u8) as u8;
        index += 1;
    }
    table
}

const NUCLEOTIDE_RESIDUES: [u8; 256] = kind_table(&NUCLEOTIDE_CHAR_TYPES, CharType::Residue);
const NUCLEOTIDE_MASKED: [u8; 256] = kind_table(&NUCLEOTIDE_CHAR_TYPES, CharType::MaskedResidue);
const PROTEIN_RESIDUES: [u8; 256] = kind_table(&PROTEIN_CHAR_TYPES, CharType::Residue);
const PROTEIN_MASKED: [u8; 256] = kind_table(&PROTEIN_CHAR_TYPES, CharType::MaskedResidue);

/// The number of bytes at the start of `s` that `table` marks: first 32 at a time while
/// `sure` (a subset of the table's bytes, tested without it) holds for them, then eight at
/// a time with the table (LOSAT's speed only).
fn run_length(s: &[u8], table: &[u8; 256], sure: impl Fn(u8) -> bool) -> usize {
    let mut length = 0;
    for block in s.chunks_exact(32) {
        if !block.iter().fold(true, |all, &byte| all & sure(byte)) {
            break;
        }
        length += 32;
    }
    for chunk in s[length..].chunks_exact(8) {
        let all = table[chunk[0] as usize]
            & table[chunk[1] as usize]
            & table[chunk[2] as usize]
            & table[chunk[3] as usize]
            & table[chunk[4] as usize]
            & table[chunk[5] as usize]
            & table[chunk[6] as usize]
            & table[chunk[7] as usize];
        if all == 0 {
            break;
        }
        length += 8;
    }
    length
        + s[length..]
            .iter()
            .take_while(|&&byte| table[byte as usize] != 0)
            .count()
}

/// Nucleotide residues and masked residues that `run_length` tests without the table.
fn sure_nucleotide(byte: u8) -> bool {
    (byte == b'A') | (byte == b'C') | (byte == b'G') | (byte == b'T') | (byte == b'N')
}

fn sure_nucleotide_masked(byte: u8) -> bool {
    (byte == b'a') | (byte == b'c') | (byte == b'g') | (byte == b't') | (byte == b'n')
}

/// Protein residues (every upper-case letter and `*`) and masked residues (every
/// lower-case letter).
fn sure_protein(byte: u8) -> bool {
    (byte.wrapping_sub(b'A') < 26) | (byte == b'*')
}

fn sure_protein_masked(byte: u8) -> bool {
    byte.wrapping_sub(b'a') < 26
}

/// How `CheckDataLine` counts a byte.
#[derive(Clone, Copy)]
enum CheckClass {
    Good,
    Neutral,
    Comment,
    Bad,
}

/// `CheckDataLine`'s classes of the bytes (`check_data_line`), looked up in a table
/// (LOSAT's speed only): letters and `*` are good; hyphens (`fHyphensIgnoreAndWarn`),
/// white space and digits count for nothing; `;` ends the check; the others are bad.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:730-750
const CHECK_CLASSES: [CheckClass; 256] = {
    let mut table = [CheckClass::Bad; 256];
    let mut index = 0;
    while index < 256 {
        let byte = index as u8;
        table[index] = if byte.is_ascii_alphabetic() || byte == b'*' {
            CheckClass::Good
        } else if byte == b'-' || is_space(byte) || byte.is_ascii_digit() {
            CheckClass::Neutral
        } else if byte == b';' {
            CheckClass::Comment
        } else {
            CheckClass::Bad
        };
        index += 1;
    }
    table
};

/// Records shorter than this give a copy of the residue buffer, which is kept for the
/// next record; longer ones take the buffer.
const SEQ_DATA_KEPT: usize = 1 << 20;

/// NCBI's `CFastaReader` state of the record being read (`ReadOneSeq` resets it).
pub(crate) struct FastaReader<R: Read> {
    pub(crate) lines: LineReader<R>,
    config: ReaderConfig,
    /// `CSeqIdGenerator::m_Counter`.
    id_counter: i64,
    /// `m_SeqData` up to `m_CurrentPos` (the raw residues, upper case).
    seq_data: Vec<u8>,
    /// `m_Gaps`.
    gaps: Vec<Gap>,
    /// `m_TotalGapLength`.
    total_gap_length: u32,
    /// `m_MaskRangeStart` (`kInvalidSeqPos` is `None`).
    mask_range_start: Option<u32>,
    /// `m_CurrentMask`'s intervals: the closed ranges of `x_CloseMask`, in positions with
    /// gaps (`ePosWithGapsAndSegs`, `m_SegmentBase` 0).
    mask: Vec<(u32, u32)>,
    /// `m_CurrentSeqTitles`.
    titles: Vec<LineText>,
    /// The title of the record (`x_ApplyMods`' `processed_title`).
    title: Vec<u8>,
    /// The local ID of the record (`GenerateID`), once its defline is read; with
    /// `-parse_deflines` the id that `Query=` shows.
    local_id: Option<String>,
    /// The Seq-ids of the record under `-parse_deflines`, once its defline is read.
    believed: Option<std::sync::Arc<super::believed::BelievedIds>>,
    /// Whether the ids of the record name a protein (`believed::named_molecule`).
    believed_molecule: Option<bool>,
    /// The keys of the ids of the subjects read so far (`believed::duplicate_id`).
    believed_keys: std::collections::HashSet<u64>,
    /// The ids of the other role's records of the search (`set_id_scope`).
    id_scope: Option<std::sync::Arc<crate::seqid::id_scope::IdScope>>,
    /// The defline of the record, as written, when its ids are too (`IdScope`'s records without a
    /// mark have generated ids).
    believed_written: Vec<u8>,
    /// The defline line as read, kept when the configuration asks for it
    /// (`ReaderConfig::keep_defline`).
    defline: Vec<u8>,
    /// Bad positions of a data line, reused.
    bad_positions: Vec<u32>,
    /// EXPERIMENT (LOSAT_X_BNQPACK): a nucleotide record's letters are packed as it is assembled
    /// (`x_packed`), its `sequence` left empty (`FastaInputSource::x_pack_letters`).
    x_pack: bool,
    /// The packed letters of the record assembled last, until the caller takes them.
    x_packed: Option<crate::sequence::packed_letters::PackedLetters>,
}

/// Where a reader message goes: a warning posted with `LOG_POST_X` (subcode 1) with no listener.
pub(crate) type Warn<'a> = &'a mut dyn FnMut(&[u8]) -> std::io::Result<()>;

impl<R: Read> FastaReader<R> {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:364-367
    pub(crate) fn new(lines: LineReader<R>, config: ReaderConfig) -> Self {
        Self {
            lines,
            config,
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:48-56
            id_counter: 1,
            seq_data: Vec::new(),
            gaps: Vec::new(),
            total_gap_length: 0,
            mask_range_start: None,
            mask: Vec::new(),
            titles: Vec::new(),
            title: Vec::new(),
            local_id: None,
            believed: None,
            believed_molecule: None,
            believed_keys: std::collections::HashSet::new(),
            id_scope: None,
            believed_written: Vec::new(),
            defline: Vec::new(),
            bad_positions: Vec::new(),
            x_pack: false,
            x_packed: None,
        }
    }

    /// LOSAT_X_QSTREAM's second reading in parts (`FastaInputSource::x_from_reader_at`): a reader
    /// that has read nothing counts its lines from `lines` and its local IDs from `next_id`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:364-367
    /// No NCBI counterpart for the start itself: NCBI reads the query file once, with one reader
    /// whose counter starts at the configuration's value (1) and whose line count starts at 0
    /// (`LineReader::new`); a reader of a part starts both where that one reader has them at the
    /// part's first byte. It reads the same records (scheduling only).
    #[cfg_attr(not(all(unix, feature = "parallel")), allow(dead_code))]
    pub(crate) fn x_start_at(&mut self, lines: u64, next_id: i64) {
        self.lines.x_start_at(lines);
        self.id_counter = next_id;
    }

    /// The search whose records this reader reads under `-parse_deflines`.
    pub(crate) fn set_id_scope(
        &mut self,
        id_scope: Option<std::sync::Arc<crate::seqid::id_scope::IdScope>>,
    ) {
        self.id_scope = id_scope;
    }

    /// EXPERIMENT (LOSAT_X_BNQPACK): pack the letters of the nucleotide records read from now on.
    pub(crate) fn x_pack_letters(&mut self) {
        self.x_pack = true;
    }

    /// EXPERIMENT (LOSAT_X_BNQPACK): the packed letters of the record read last, if it was packed
    /// (its `sequence` is then empty).
    pub(crate) fn x_take_packed(
        &mut self,
    ) -> Option<crate::sequence::packed_letters::PackedLetters> {
        self.x_packed.take()
    }

    /// `LineNumber()`.
    fn line_number(&self) -> u64 {
        self.lines.line_number()
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2196-2214
    /// The BLAST input source gives the reader no listener, so a warning is the message
    /// on a line of its own (`LOG_POST` writes no severity).
    fn post_warning(
        &self,
        problem: Problem,
        message: &[u8],
        warn: Warn<'_>,
    ) -> Result<(), ReadError> {
        post_line_warning(problem, message, warn)
    }

    /// `FASTA_ERROR`: `PostWarning` with `eDiag_Error` and no listener throws
    /// `CObjReaderParseException`, which gets the error code, the message and the line number.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2215-2218
    fn fasta_error(&self, code: super::ParseErrorCode, message: String) -> ReadError {
        ReadError::Parse {
            code,
            message,
            line: self.line_number(),
        }
    }

    /// One record (`ReadOneSeq` of `CBlastInputReader`, then of `CFastaReader`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:130-162
    /// NCBI reads a line that is a Seq-id as a sequence to fetch with its data loaders
    /// (GenBank over the network, or a BLAST database), which LOSAT does not do: such a
    /// line is an explicit rejection (`seq_id::reject_seq_id_line`). The rejection leaves
    /// the reader where NCBI's Seq-id record leaves it: the line is consumed (no
    /// `UngetLine` on that path) and the local-ID counter has not moved, so a further call
    /// reads from the next line (oracle BI net2, net5).
    pub(crate) fn read_one_seq(&mut self, warn: Warn<'_>) -> Result<FastaRecord, ReadError> {
        self.read_one_seq_inner(None, warn)
            .map(|(record, _)| record)
    }

    /// Reads one record (`read_one_seq`) keeping `sequence` in upper case and returning the
    /// lowercase mask intervals directly (`self.mask`) when `keep_masks` is set (empty otherwise),
    /// avoiding a lowercase-and-rescan round trip for callers whose record stores upper-case
    /// letters with interval masks (`algorithm::blastx::input_source`).
    pub(crate) fn read_one_seq_upper(
        &mut self,
        keep_masks: bool,
        warn: Warn<'_>,
    ) -> Result<UpperRecord, ReadError> {
        self.read_one_seq_inner(Some(keep_masks), warn)
    }

    fn read_one_seq_inner(
        &mut self,
        upper: Option<bool>,
        warn: Warn<'_>,
    ) -> Result<UpperRecord, ReadError> {
        if self.config.data_loaders {
            let (program, role) = (self.config.program, self.config.role);
            let line = self.lines.next_line();
            if let Some(rejection) = super::seq_id::reject_seq_id_line(line, program, role) {
                return Err(ReadError::Unsupported(rejection));
            }
            self.lines.unget_line();
        }
        self.read_fasta_seq(upper, warn)
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:312-440
    /// `FASTA_PROGRESS` posts to a listener only, and there is none. `CBadResiduesException`
    /// comes only with `fValidate`, which is not set.
    fn read_fasta_seq(
        &mut self,
        upper: Option<bool>,
        warn: Warn<'_>,
    ) -> Result<UpperRecord, ReadError> {
        self.seq_data.clear();
        self.gaps.clear();
        self.mask_range_start = None;
        self.mask.clear();
        self.total_gap_length = 0;
        self.titles.clear();
        self.title.clear();
        self.local_id = None;
        self.believed = None;
        self.believed_molecule = None;

        let mut need_defline = true;
        while !self.lines.at_eof() {
            let c = self.lines.peek_char();
            if self.lines.at_eof() {
                return Err(self.fasta_error(
                    super::ParseErrorCode::Eof,
                    format!(
                        "CFastaReader: Unexpected end-of-file around line {}",
                        self.line_number()
                    ),
                ));
            }
            if c == Some(b'>') {
                self.lines.next_line();
                let raw = self.lines.take_line();
                let modified;
                let next_line = match raw.strip_prefix(b">?_") {
                    Some(rest) => {
                        modified = [b">".as_slice(), rest].concat();
                        &modified
                    }
                    None => &raw[..],
                };
                if next_line.starts_with(b">?") {
                    self.lines.put_line(raw);
                    self.lines.unget_line();
                } else if need_defline {
                    let parsed = if self.config.parse_deflines {
                        self.parse_believed_def_line(next_line, warn)
                    } else {
                        self.parse_def_line(next_line);
                        Ok(())
                    };
                    self.lines.put_line(raw);
                    parsed?;
                    need_defline = false;
                    continue;
                } else {
                    self.lines.put_line(raw);
                    self.lines.unget_line();
                    break;
                }
            }
            self.lines.next_line();
            let raw = self.lines.take_line();
            // `NStr::TruncateSpaces_Unsafe` of the line (`eTrunc_Both`, fasta.cpp:376).
            let line = truncate_spaces(&raw);
            if line.is_empty() || matches!(line[0], b'!' | b'#' | b';') {
                self.lines.put_line(raw);
                continue;
            }
            if need_defline {
                if self.config.parse_deflines {
                    let line_number = self.line_number();
                    self.lines.put_line(raw);
                    return Err(ReadError::Parse {
                        code: super::ParseErrorCode::NoDefline,
                        message: format!(
                            "CFastaReader: Input doesn't start with a defline or comment around \
                             line {line_number}"
                        ),
                        line: line_number,
                    });
                }
                self.parse_def_line(b">");
                need_defline = false;
            }
            let parsed = match line.strip_prefix(b">?_") {
                Some(rest) => self.parse_data_line(&[b">".as_slice(), rest].concat(), warn),
                None => self.parse_data_line(line, warn),
            };
            self.lines.put_line(raw);
            parsed?;
        }
        if need_defline && self.lines.at_eof() {
            return Err(self.fasta_error(
                super::ParseErrorCode::Eof,
                format!(
                    "CFastaReader: Expected defline around line {}",
                    self.line_number()
                ),
            ));
        }
        self.assemble_seq(upper, warn)
    }

    /// The title of a defline (`ParseDefLine` with `fNoParseID`), and the generated local
    /// ID (`PostProcessIDs` without IDs).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:146-226
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:566-574
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:619-631
    /// The `CFastaReader` user object with the raw defline is not output.
    fn parse_def_line(&mut self, defline: &[u8]) {
        self.titles.clear();
        if self.config.keep_defline {
            self.defline.clear();
            self.defline.extend_from_slice(defline);
        }
        let len = defline.len();
        if len > 1 && !is_blank(&defline[1..]) {
            let start = 1 + defline[1..]
                .iter()
                .position(|&byte| !is_space(byte))
                .unwrap_or(len - 1);
            let title_start = start
                + defline[start..]
                    .iter()
                    .position(|&byte| !is_space(byte))
                    .unwrap_or(len - start);
            if title_start < len {
                let end = title_start
                    + 1
                    + defline[title_start + 1..]
                        .iter()
                        .position(|&byte| byte < b' ')
                        .unwrap_or(len - title_start - 1);
                self.titles.push(LineText {
                    text: defline[title_start..end].to_vec(),
                });
            }
        }
        self.local_id = Some(self.generate_id());
    }

    /// `ParseDefLine` with `fParseRawID` (`-parse_deflines`, `believed`): the ids of the first word
    /// (a generated local id for none, and for a word with a range), and the title after it.
    fn parse_believed_def_line(&mut self, defline: &[u8], warn: Warn<'_>) -> Result<(), ReadError> {
        use super::believed::{
            duplicate_id, generated_ids, named_molecule, parse_believed_defline,
            refuse_range_in_archive, refused_id, BelievedIds,
        };
        self.titles.clear();
        if self.config.keep_defline {
            self.defline.clear();
            self.defline.extend_from_slice(defline);
        }
        let line_number = self.line_number();
        let mut written = Ok(());
        let parsed = parse_believed_defline(
            defline,
            self.config.protein,
            self.config.program,
            line_number,
            &mut |message: &[u8]| {
                let mut line = message.to_vec();
                line.push(b'\n');
                if written.is_ok() {
                    written = warn(&line);
                }
            },
        );
        written.map_err(ReadError::Write)?;
        let parsed = parsed.map_err(|message| ReadError::Parse {
            code: super::ParseErrorCode::Format,
            message,
            line: line_number,
        })?;
        let len = defline.len();
        let mut title_start = parsed.title_search;
        while title_start < len && is_space(defline[title_start]) {
            title_start += 1;
        }
        if title_start < len {
            let end = title_start
                + 1
                + defline[title_start + 1..]
                    .iter()
                    .position(|&byte| byte < b' ')
                    .unwrap_or(len - title_start - 1);
            self.titles.push(LineText {
                text: defline[title_start..end].to_vec(),
            });
        }
        if parsed.range && self.config.keep_defline && self.config.archive_known {
            return Err(ReadError::Unsupported(refuse_range_in_archive(
                self.config.role,
                self.config.program,
            )));
        }
        self.believed_molecule = named_molecule(&parsed.ids);
        self.believed_written.clear();
        if self.id_scope.is_some() && !parsed.ids.is_empty() && !parsed.range {
            self.believed_written.extend_from_slice(defline);
        }
        let ids = if parsed.ids.is_empty() || parsed.range {
            generated_ids(&self.generate_id())
        } else {
            parsed.ids
        };
        // Two queries with one id are the scopes' of the query batches (`QueryScopes`), which
        // NCBI's formatter empties after each batch against a database.
        let refusal = refused_id(
            &ids,
            self.config.role,
            "-parse_deflines",
            self.config.program,
        )
        .or_else(|| {
            (self.config.role != "query")
                .then(|| {
                    duplicate_id(
                        &ids,
                        &mut self.believed_keys,
                        self.config.role,
                        self.config.program,
                    )
                })
                .flatten()
        });
        if let Some(refusal) = refusal {
            return Err(ReadError::Unsupported(refusal));
        }
        let believed =
            BelievedIds::new(ids, self.config.strand).with_named_molecule(self.believed_molecule);
        self.local_id = Some(believed.bare_id());
        self.believed = Some(std::sync::Arc::new(believed));
        Ok(())
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:466-483
    fn generate_id(&mut self) -> String {
        let n = self.id_counter;
        self.id_counter += 1;
        format!("{}{n}", self.config.id_prefix)
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:713-772
    /// `bIsNuc` only feeds `eProblem_TooManyAmbiguousResidues`, which `CBlastFastaInputSource`
    /// ignores (`blast_fasta_input.cpp:360`), and the other warning needs `good == 0`, which
    /// the format error before has already thrown.
    fn check_data_line(&self, s: &[u8]) -> Result<(), ReadError> {
        if !self.seq_data.is_empty() {
            return Ok(());
        }
        let (mut good, mut bad) = (0usize, 0usize);
        let len_to_check = s.len().min(70);
        for &c in &s[..len_to_check] {
            match CHECK_CLASSES[c as usize] {
                CheckClass::Good => good += 1,
                CheckClass::Neutral => {}
                CheckClass::Comment => break,
                CheckClass::Bad => bad += 1,
            }
        }
        if bad >= good / 3 && (len_to_check > 3 || good == 0 || bad > good) {
            let line = self.line_number();
            return Err(self.fasta_error(
                super::ParseErrorCode::Format,
                format!(
                    "CFastaReader: Near line {line}, there's a line that doesn't look like plausible data, but it's not marked as defline or comment."
                ),
            ));
        }
        Ok(())
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:774-1016 (the
    /// character switch is cited at `char_types`)
    /// `CloseGap` does nothing (no gap characters without `fParseGaps`, and
    /// `CCustomizedFastaReader::x_CloseGap` is empty); `x_NucOrProt()` names the molecule set
    /// by `xSetSeqMol` under `-parse_deflines` and is empty otherwise. The fast copy of
    /// `fSkipCheck` is not taken (no `fSkipCheck`).
    fn parse_data_line(&mut self, s: &[u8], warn: Warn<'_>) -> Result<(), ReadError> {
        if s.starts_with(b">?") {
            return self.parse_gap_line(s, warn);
        }
        self.check_data_line(s)?;
        let parse_protein =
            DataLineAlphabet::new(self.believed_molecule, self.config.protein).is_protein();
        let (types, residues, masked) = if parse_protein {
            (&PROTEIN_CHAR_TYPES, &PROTEIN_RESIDUES, &PROTEIN_MASKED)
        } else {
            (
                &NUCLEOTIDE_CHAR_TYPES,
                &NUCLEOTIDE_RESIDUES,
                &NUCLEOTIDE_MASKED,
            )
        };
        let mut hyphen_seen = false;
        self.bad_positions.clear();
        let bulk = self.lines.bulk();
        let mut pos = 0;
        while pos < s.len() {
            let c = s[pos];
            match types[c as usize] {
                // After the first residue of a run, `CloseMask` (`OpenMask`) has nothing
                // left to do for the others, which are copied at once (LOSAT's speed).
                CharType::Residue => {
                    self.seq_data.push(c);
                    self.close_mask();
                    let rest = &s[pos + 1..];
                    let run = if !bulk {
                        0
                    } else if parse_protein {
                        run_length(rest, residues, sure_protein)
                    } else {
                        run_length(rest, residues, sure_nucleotide)
                    };
                    self.seq_data.extend_from_slice(&s[pos + 1..pos + 1 + run]);
                    pos += run;
                }
                CharType::MaskedResidue => {
                    self.seq_data.push(c.to_ascii_uppercase());
                    self.open_mask();
                    let rest = &s[pos + 1..];
                    let run = if !bulk {
                        0
                    } else if parse_protein {
                        run_length(rest, masked, sure_protein_masked)
                    } else {
                        run_length(rest, masked, sure_nucleotide_masked)
                    };
                    self.seq_data
                        .extend(s[pos + 1..pos + 1 + run].iter().map(u8::to_ascii_uppercase));
                    pos += run;
                }
                CharType::HyphenToIgnoreAndWarn => hyphen_seen = true,
                CharType::Comment => break,
                CharType::Space => {}
                CharType::Bad => self.bad_positions.push(pos as u32),
            }
            pos += 1;
        }
        self.check_length()?;
        if hyphen_seen {
            let line = self.line_number();
            self.post_warning(
                Problem::IgnoredResidue,
                format!("CFastaReader: Hyphens are invalid and will be ignored around line {line}")
                    .as_bytes(),
                warn,
            )?;
        }
        if !self.bad_positions.is_empty() {
            let mut message = format!(
                "FASTA-Reader: Ignoring invalid {}residues at position(s): ",
                super::believed::nuc_or_prot(self.believed_molecule)
            )
            .into_bytes();
            convert_bad_indexes_to_string(&mut message, self.line_number(), &self.bad_positions);
            self.post_warning(Problem::InvalidResidue, &message, warn)?;
        }
        Ok(())
    }

    /// LOSAT's limit: NCBI keeps positions in `TSeqPos` (32 bits) and BLAST's sequence
    /// lengths in `Int4`; a record longer than that is an explicit rejection
    /// (`super::record_too_long`).
    fn check_length(&self) -> Result<(), ReadError> {
        if self.seq_data.len() as u64 + self.total_gap_length as u64 > i32::MAX as u64 {
            return Err(super::record_too_long(&self.config));
        }
        Ok(())
    }

    /// `GetCurrentPos` with `ePosWithGapsAndSegs`: the raw position plus the gaps before it
    /// (`m_SegmentBase` is 0 without `fInSegSet`).
    ///
    /// NCBI reference (598d8ae6): c++/include/objtools/readers/fasta.hpp:502-516
    fn current_pos_with_gaps(&self) -> u32 {
        self.seq_data.len() as u32 + self.total_gap_length
    }

    /// The mask is always recorded (LOSAT marks the masked residues by case, and the
    /// programs apply the mask only with `-lcase_masking`, where NCBI's
    /// `x_FastaToSeqLoc` calls `SaveMask` before the record is read).
    ///
    /// NCBI reference (598d8ae6): c++/include/objtools/readers/fasta.hpp:494-499
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1079-1083
    /// `x_OpenMask` runs before `m_CurrentPos` is incremented, so the start is the masked residue.
    fn open_mask(&mut self) {
        if self.mask_range_start.is_none() {
            self.mask_range_start = Some(self.current_pos_with_gaps() - 1);
        }
    }

    /// NCBI reference (598d8ae6): c++/include/objtools/readers/fasta.hpp:291-292
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1085-1092
    /// In `ParseDataLine` `CloseMask` runs after the residue is stored and before
    /// `m_CurrentPos` is incremented: the interval ends before the residue that closes it (LOSAT has
    /// already pushed that residue, hence the `- 2`). `AssembleSeq` closes it at the end.
    fn close_mask(&mut self) {
        if let Some(start) = self.mask_range_start.take() {
            self.mask.push((start, self.current_pos_with_gaps() - 2));
        }
    }

    /// `ParseGapLine` (`gap_line_size`): the gap of the line, as a run of `N` or `X` that
    /// `AssembleSeq` puts at this position.
    fn parse_gap_line(&mut self, line: &[u8], warn: Warn<'_>) -> Result<(), ReadError> {
        let position = self.seq_data.len() as u32;
        let gap_size = gap_line_size(line, self.line_number(), warn)?;
        self.gaps.push(Gap {
            position,
            length: gap_size,
        });
        self.total_gap_length = self.total_gap_length.saturating_add(gap_size);
        self.check_length()?;
        Ok(())
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1351-1448
    /// The molecule is the input's (`CCustomizedFastaReader::AssignMolType` below
    /// `m_SeqLenThreshold`, which the BLAST input source sets to `UINT_MAX`), so the gap
    /// residue is `N` for nucleotides and `X` for proteins, and `x_NucOrProt` is not used
    /// after this point. Without `fNoSplit` (`BLASTINPUT_GEN_DELTA_SEQ` set) the gaps have
    /// already become residues, so the record is the same sequence in delta pieces.
    /// `CSeqportUtil::Pack` packs nucleotides into `ncbi2na`/`ncbi4na`, where `U` and `T`
    /// are the same code and any non-IUPACna byte read under a protein Seq-id maps to `N`
    /// (`CIupacnaTo4na::scm_Table`, `sequtil_tables.cpp:392-457`), so LOSAT stores `U` as `T`.
    /// The residues in the lowercase mask are returned in lower case (LOSAT's representation of
    /// the mask), or as mask intervals when `upper_masks` is `Some`.
    fn assemble_seq(
        &mut self,
        upper_masks: Option<bool>,
        warn: Warn<'_>,
    ) -> Result<UpperRecord, ReadError> {
        // CloseMask at the end of the record: the interval ends at the last position.
        if let Some(start) = self.mask_range_start.take() {
            self.mask.push((start, self.current_pos_with_gaps() - 1));
        }
        let mut titles = std::mem::take(&mut self.titles);
        for title in titles.drain(..) {
            self.parse_title(title, warn)?;
        }
        self.titles = titles;
        if !self.config.protein && self.believed_molecule == Some(true) {
            pack_protein_read_nucleotide(&mut self.seq_data);
        }
        self.add_to_id_scope()?;
        let gap_char = if self.config.protein { b'X' } else { b'N' };
        // EXPERIMENT (LOSAT_X_BNQPACK): the letters the branches below would give, packed from the
        // buffer (`U` as `T` first, the mask as lower-case intervals) instead of copied out of it.
        self.x_packed = None;
        if self.x_pack && self.total_gap_length == 0 && !self.config.protein {
            for residue in self.seq_data.iter_mut() {
                *residue = if *residue == b'U' { b'T' } else { *residue };
            }
            self.x_packed = crate::sequence::packed_letters::PackedLetters::pack_masked(
                &self.seq_data,
                &self.mask,
            );
        }
        let mut sequence = if self.x_packed.is_some() {
            if self.seq_data.len() < SEQ_DATA_KEPT {
                self.seq_data.clear();
            } else {
                self.seq_data = Vec::new();
            }
            Vec::new()
        } else if self.total_gap_length > 0 {
            let mut new_data =
                Vec::with_capacity(self.seq_data.len() + self.total_gap_length as usize);
            let mut pos = 0usize;
            for gap in &self.gaps {
                let gap_pos = gap.position as usize;
                if gap_pos > pos {
                    new_data.extend_from_slice(&self.seq_data[pos..gap_pos]);
                    pos = gap_pos;
                }
                new_data.resize(new_data.len() + gap.length as usize, gap_char);
            }
            if self.seq_data.len() > pos {
                new_data.extend_from_slice(&self.seq_data[pos..]);
            }
            new_data
        } else if self.seq_data.len() < SEQ_DATA_KEPT {
            // The buffer stays for the next record (LOSAT's speed: no regrowth per record).
            let sequence = self.seq_data.clone();
            self.seq_data.clear();
            sequence
        } else {
            // LOSAT's memory: the record keeps no growth slack of the buffer.
            let mut sequence = std::mem::take(&mut self.seq_data);
            sequence.shrink_to_fit();
            sequence
        };
        if !self.config.protein {
            // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1444-1447
            for residue in sequence.iter_mut() {
                // A select rather than a conditional store (LOSAT's speed: vectorized).
                *residue = if *residue == b'U' { b'T' } else { *residue };
            }
        }
        // A packed record's mask is in its intervals; callers asking for upper-case letters with
        // mask intervals (`upper_masks`) take them straight from `self.mask`.
        let masks = match upper_masks {
            Some(keep_masks) => {
                if keep_masks && !self.mask.is_empty() {
                    // `check_length` enforces `<= i32::MAX`.
                    self.mask
                        .iter()
                        .map(|&(from, to)| (from as i32, to as i32))
                        .collect()
                } else {
                    Vec::new()
                }
            }
            None => {
                if self.x_packed.is_none() {
                    for &(from, to) in &self.mask {
                        for residue in &mut sequence[from as usize..=to as usize] {
                            *residue = residue.to_ascii_lowercase();
                        }
                    }
                }
                Vec::new()
            }
        };
        Ok((
            FastaRecord {
                local_id: self.local_id.take().unwrap_or_default(),
                title: std::mem::take(&mut self.title),
                sequence,
                warnings: Vec::new(),
                source: (self.config.keep_defline || self.believed.is_some()).then(|| {
                    Box::new(super::RecordSource {
                        defline: std::mem::take(&mut self.defline),
                        uncut: None,
                        believed: self.believed.take(),
                        cut_from: 0,
                        uncut_length: 0,
                    })
                }),
            },
            masks,
        ))
    }

    /// `-parse_deflines` (`IdScope::add_record`): the record's letters are complete and its letters
    /// in `seq_data` are in upper case, whether or not they are packed next. Cold: no reader without
    /// the option has an id scope.
    fn add_to_id_scope(&mut self) -> Result<(), ReadError> {
        use crate::seqid::id_scope::RecordMark;
        let (Some(scope), Some(believed)) = (&self.id_scope, &mut self.believed) else {
            return Ok(());
        };
        let mark = (!self.believed_written.is_empty() && !believed.subject_local()).then(|| {
            RecordMark::new(
                &self.believed_written,
                self.config.protein,
                &self.seq_data,
                self.gaps.iter().map(|gap| (gap.position, gap.length)),
                !self.mask.is_empty(),
            )
        });
        let shared = scope
            .add_record(&believed.ids, mark, self.config.role, self.config.program)
            .map_err(ReadError::Unsupported)?;
        // The defline made the only handle of the ids.
        if let Some(ids) = std::sync::Arc::get_mut(believed) {
            ids.shared = shared;
        }
        Ok(())
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:670-686
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2018-2043
    /// The long-title and modifier warnings are ignored problems.
    fn parse_title(&mut self, title: LineText, warn: Warn<'_>) -> Result<(), ReadError> {
        self.create_warnings_for_seq_data_in_title(&title.text, warn)?;
        // `TruncateSpacesInPlace` (`eTrunc_Both`) on the title's own bytes.
        let mut text = title.text;
        text.truncate(truncate_spaces_end(&text).len());
        let leading = text.len() - truncate_spaces_begin(&text).len();
        text.drain(..leading);
        self.title = text;
        Ok(())
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1623-1680
    /// Exactly one of `fAssumeNuc` and `fAssumeProt` is set, so exactly one test runs, from
    /// the end of the title.
    fn create_warnings_for_seq_data_in_title(
        &self,
        text: &[u8],
        warn: Warn<'_>,
    ) -> Result<(), ReadError> {
        if let Some((message, protein)) = seq_data_in_title_warning(text, Some(self.config.protein))
        {
            let problem = if protein {
                Problem::UnexpectedAminoAcids
            } else {
                Problem::UnexpectedNucResidues
            };
            self.post_warning(problem, message, warn)?;
        }
        Ok(())
    }
}

/// The alphabet molecule `CFastaReader::ParseDataLine` checks a line against: the Seq-id's
/// named molecule when set (`xSetSeqMol`), otherwise the program's.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:809-814
#[derive(Clone, Copy)]
pub struct DataLineAlphabet {
    parse_protein: bool,
}

impl DataLineAlphabet {
    #[inline]
    pub fn new(believed_molecule: Option<bool>, program_protein: bool) -> Self {
        Self {
            parse_protein: believed_molecule.unwrap_or(program_protein),
        }
    }

    #[inline]
    pub fn is_protein(self) -> bool {
        self.parse_protein
    }
}

/// How `CSeqportUtil::Pack` maps upper-case residues read under a protein Seq-id into a
/// nucleotide sequence (`CIupacnaTo4na::scm_Table`): IUPACna letters stay, `U` becomes `T`,
/// and every other byte (`E`, `F`, `I`, `J`, `L`, `O`, `P`, `Q`, `X`, `Z`, `*`) becomes `N`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1444-1448
/// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:392-457
#[inline]
pub fn pack_protein_read_nucleotide(seq_data: &mut [u8]) {
    for residue in seq_data.iter_mut() {
        *residue = match *residue {
            b'A' | b'B' | b'C' | b'D' | b'G' | b'H' | b'K' | b'M' | b'N' | b'R' | b'S' | b'T'
            | b'V' | b'W' | b'Y' => *residue,
            b'U' => b'T',
            _ => b'N',
        };
    }
}

/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1094-1349
/// The rest of the function (the modifiers, cited at `gap_mod_warnings`) only warns:
/// without `fParseGaps` `AssembleSeq` turns each gap into a run of `N` or `X`.
///
/// The size of the gap line `line` (it starts with `>?`), line `line_number` of the input, after
/// its warnings: the reader of the four programs (`FastaReader::parse_gap_line`) and BLASTX's
/// (`blastx::input_source`) read gap lines with it.
pub(crate) fn gap_line_size(
    line: &[u8],
    line_number: u64,
    warn: Warn<'_>,
) -> Result<u32, ReadError> {
    let mut remaining = truncate_spaces(&line[2..]);
    if let Some(rest) = remaining.strip_prefix(b"unk") {
        remaining = truncate_spaces_begin(rest);
    }
    let digits = remaining
        .iter()
        .position(|byte| !byte.is_ascii_digit())
        .unwrap_or(remaining.len());
    // `NStr::StringToUInt` of `sDigits` with `NStr::fConvErr_NoThrow`: 0 for an empty string or a value
    // over UINT_MAX.
    let mut gap_size = crate::ncbi_str::string_to_uint(&remaining[..digits]).unwrap_or(0);
    if gap_size == 0 {
        post_line_warning(
            Problem::NonPositiveLength,
            format!("CFastaReader: Bad gap size at line {line_number}").as_bytes(),
            warn,
        )?;
        gap_size = 1;
    }
    remaining = truncate_spaces_begin(&remaining[digits..]);
    let mut mods: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    while !remaining.is_empty() {
        let close = if remaining.starts_with(b"[") {
            remaining[1..]
                .iter()
                .position(|&byte| byte == b'=')
                .map(|equal| equal + 1)
                .and_then(|equal| {
                    remaining[equal + 1..]
                        .iter()
                        .position(|&byte| byte == b']')
                        .map(|close| (equal, equal + 1 + close))
                })
        } else {
            None
        };
        let Some((equal, close)) = close else {
            post_line_warning(
                Problem::ParsingModifiers,
                format!("CFastaReader: Problem parsing gap mods at line {line_number}").as_bytes(),
                warn,
            )?;
            break;
        };
        let key = truncate_spaces(&remaining[1..equal]).to_vec();
        let value = truncate_spaces(&remaining[equal + 1..close]).to_vec();
        mods.push((key, value));
        remaining = truncate_spaces_begin(&remaining[close + 1..]);
    }
    // `multimap` iterates in key order (stable for equal keys).
    mods.sort_by(|left, right| left.0.cmp(&right.0));
    gap_mod_warnings(line_number, &mods, warn)?;
    Ok(gap_size)
}

/// The warnings of `ParseGapLine` about the gap modifiers. Without `fParseGaps` the
/// gap type and linkage evidence are not kept (`AssembleSeq` warns that the
/// modifiers are ignored, an ignored problem).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1190-1340
/// The caller's default gap type (`m_gap_type`) is not set by the BLAST input source.
/// `NStr::Split` with `fSplit_MergeDelimiters` and `fSplit_Truncate` drops empty pieces.
fn gap_mod_warnings(
    line: u64,
    mods: &[(Vec<u8>, Vec<u8>)],
    warn: Warn<'_>,
) -> Result<(), ReadError> {
    let mut conflicting = false;
    let mut gap_type: Option<GapType> = None;
    let mut link_evid = LinkEvid::UnspecifiedOnly;
    let mut evidence: std::collections::BTreeSet<u8> = std::collections::BTreeSet::new();
    for (key, value) in mods {
        let canonical_key = canonicalize_string(key);
        if canonical_key == b"gap-type" {
            match gap_type_info(value) {
                Some((value_type, value_evid)) => match gap_type {
                    None => {
                        gap_type = Some(value_type);
                        link_evid = value_evid;
                    }
                    Some(current) if current != value_type => conflicting = true,
                    Some(_) => {}
                },
                None => {
                    let mut message = b"Unknown gap-type: ".to_vec();
                    message.extend_from_slice(value);
                    post_line_warning(Problem::ParsingModifiers, &message, warn)?;
                }
            }
            continue;
        }
        if canonical_key == b"linkage-evidence" {
            for piece in value
                .split(|&byte| byte == b';')
                .filter(|piece| !piece.is_empty())
            {
                match linkage_evidence_value(&canonicalize_string(piece)) {
                    Some(evidence_value) => {
                        evidence.insert(evidence_value);
                    }
                    None => {
                        let mut message = b"Unknown linkage-evidence: ".to_vec();
                        message.extend_from_slice(value);
                        post_line_warning(Problem::ParsingModifiers, &message, warn)?;
                    }
                }
            }
            continue;
        }
        let mut message = b"Unknown gap modifier name(s): ".to_vec();
        message.extend_from_slice(key);
        post_line_warning(Problem::UnrecognizedQualifierName, &message, warn)?;
    }
    if conflicting {
        post_line_warning(
            Problem::ContradictoryModifiers,
            format!("There were conflicting gap-types around line {line}").as_bytes(),
            warn,
        )?;
    }
    if gap_type.is_none() && !evidence.is_empty() {
        gap_type = Some(GapType::Unknown);
    }
    match link_evid {
        LinkEvid::UnspecifiedOnly => {
            if !evidence.is_empty()
                && (evidence.len() > 1 || !evidence.contains(&LINKAGE_UNSPECIFIED))
            {
                post_line_warning(
                    Problem::ExtraModifierFound,
                    b"FASTA-Reader: Unknown gap-type can have linkage-evidence of type 'unspecified' only.",
                    warn,
                )?;
            }
        }
        LinkEvid::Forbidden => {
            if !evidence.is_empty() {
                post_line_warning(
                    Problem::ModifierFoundButNoneExpected,
                    b"FASTA-Reader: This gap-type cannot have any linkage-evidence specified, so any will be ignored.",
                    warn,
                )?;
            }
        }
        LinkEvid::Required => {
            if evidence.is_empty() {
                evidence.insert(LINKAGE_UNSPECIFIED);
            }
            if evidence.len() == 1 && evidence.contains(&LINKAGE_UNSPECIFIED) {
                post_line_warning(
                    Problem::ExpectedModifierMissing,
                    b"CFastaReader: This gap-type should have at least one specified linkage-evidence.",
                    warn,
                )?;
            }
        }
    }
    let _ = gap_type;
    Ok(())
}

/// `PostWarning` without a listener (`FastaReader::post_warning`): the message on a line of its own,
/// unless the reader ignores the problem.
fn post_line_warning(problem: Problem, message: &[u8], warn: Warn<'_>) -> Result<(), ReadError> {
    if IGNORED_PROBLEMS.contains(&problem) {
        return Ok(());
    }
    let mut line = Vec::with_capacity(message.len() + 1);
    line.extend_from_slice(message);
    line.push(b'\n');
    warn(&line).map_err(ReadError::Write)
}

/// The message of `CreateWarningsForSeqDataInTitle` for the title text `text` (see
/// `FastaReader::create_warnings_for_seq_data_in_title`, which posts it), or `None`;
/// `protein` is `fAssumeProt` (else `fAssumeNuc`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1624-1643
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1651-1673
/// Exactly one of `fAssumeNuc` and `fAssumeProt` is set, so exactly one test runs.
///
/// `assumed_protein` is `Some(true)` for `fAssumeProt`, `Some(false)` for `fAssumeNuc` and
/// `None` for neither (convert2blastmask's reader), which makes both checks: the amino acid
/// check goes on from where the nucleotide check stopped, at a byte that is no unambiguous
/// nucleotide, and the bytes after it are letters, so it is the check of the last 50 bytes.
/// The warning and whether it is the amino acid one.
pub(crate) fn seq_data_in_title_warning(
    text: &[u8],
    assumed_protein: Option<bool>,
) -> Option<(&'static [u8], bool)> {
    let length = text.len();
    if assumed_protein != Some(true)
        && length > 20
        && text[length - 20..]
            .iter()
            .all(|byte| matches!(byte, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
    {
        return Some((b"FASTA-Reader: Title ends with at least 20 valid nucleotide characters.  Was the sequence accidentally put in the title line?", false));
    }
    if assumed_protein != Some(false)
        && length > 50
        && text[length - 50..].iter().all(u8::is_ascii_alphabetic)
    {
        return Some((b"FASTA-Reader: Title ends with at least 50 valid amino acid characters.  Was the sequence accidentally put in the title line?", true));
    }
    None
}

/// The positions of the bad residues of one line, as ranges of 1-based positions.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_exception.cpp:81-150
/// `maxRanges` is the default 1000 (fasta_exception.hpp:90-92); `iRangesFound` never
/// exceeds it, so ", and more" is not written.
pub(crate) fn convert_bad_indexes_to_string(out: &mut Vec<u8>, line: u64, positions: &[u32]) {
    const MAX_RANGES: usize = 1000;
    let mut ranges: Vec<(u32, u32)> = Vec::new();
    for &idx in positions {
        if let Some(last) = ranges.last_mut() {
            if idx == last.1 + 1 {
                last.1 += 1;
                continue;
            }
            if ranges.len() >= MAX_RANGES {
                break;
            }
        }
        ranges.push((idx, idx));
    }
    out.extend_from_slice(format!("On line {line}: ").as_bytes());
    for (index, &(first, second)) in ranges.iter().enumerate() {
        if index > 0 {
            out.extend_from_slice(b", ");
        }
        out.extend_from_slice((first + 1).to_string().as_bytes());
        if first != second {
            out.extend_from_slice(format!("-{}", second + 1).as_bytes());
        }
    }
}

// The tables and sets of LOSAT's fast paths against the definitions they replace.
#[cfg(test)]
#[test]
fn fast_path_tables_match_their_definitions() {
    for byte in 0..=255u8 {
        for (types, residues, masked, sure, sure_masked) in [
            (
                &NUCLEOTIDE_CHAR_TYPES,
                &NUCLEOTIDE_RESIDUES,
                &NUCLEOTIDE_MASKED,
                sure_nucleotide as fn(u8) -> bool,
                sure_nucleotide_masked as fn(u8) -> bool,
            ),
            (
                &PROTEIN_CHAR_TYPES,
                &PROTEIN_RESIDUES,
                &PROTEIN_MASKED,
                sure_protein,
                sure_protein_masked,
            ),
        ] {
            let kind = types[byte as usize];
            assert_eq!(residues[byte as usize] == 1, kind == CharType::Residue);
            assert_eq!(masked[byte as usize] == 1, kind == CharType::MaskedResidue);
            assert!(!sure(byte) || kind == CharType::Residue, "{byte}");
            assert!(
                !sure_masked(byte) || kind == CharType::MaskedResidue,
                "{byte}"
            );
        }
        // CheckDataLine (fasta.cpp:730-750).
        let class = if byte.is_ascii_alphabetic() || byte == b'*' {
            "good"
        } else if byte == b'-' || is_space(byte) || byte.is_ascii_digit() {
            "neutral"
        } else if byte == b';' {
            "comment"
        } else {
            "bad"
        };
        let table = match CHECK_CLASSES[byte as usize] {
            CheckClass::Good => "good",
            CheckClass::Neutral => "neutral",
            CheckClass::Comment => "comment",
            CheckClass::Bad => "bad",
        };
        assert_eq!(table, class, "{byte}");
    }
    // Protein: the sure sets are the whole kinds.
    for byte in 0..=255u8 {
        assert_eq!(
            sure_protein(byte),
            PROTEIN_CHAR_TYPES[byte as usize] == CharType::Residue
        );
        assert_eq!(
            sure_protein_masked(byte),
            PROTEIN_CHAR_TYPES[byte as usize] == CharType::MaskedResidue
        );
    }
}
