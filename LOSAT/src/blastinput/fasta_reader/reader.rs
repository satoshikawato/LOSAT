//! NCBI's `CFastaReader` with the flags of `CBlastFastaInputSource`
//! (`fNoParseID | fDLOptional`, `fAssumeNuc` or `fAssumeProt`, `fNoSplit`,
//! `fHyphensIgnoreAndWarn`, `fDisableNoResidues`, `fQuickIDCheck`), as overridden by
//! `CCustomizedFastaReader` and `CBlastInputReader`.

use std::io::Read;

use super::stream::{input_space, trim_input_end, trim_input_space, trim_input_start, LineReader};
use super::{FastaRecord, ReadError, ReaderConfig};

/// The kinds of reader problems (`ILineError::EProblem`) that LOSAT meets, for
/// `PostWarning`'s list of ignored problems.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:358-360
/// ```c++
///     m_InputReader->IgnoreProblem(ILineError::eProblem_ModifierFoundButNoneExpected);
///     m_InputReader->IgnoreProblem(ILineError::eProblem_TooLong);
///     m_InputReader->IgnoreProblem(ILineError::eProblem_TooManyAmbiguousResidues);
/// ```
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
enum GapType {
    Unknown,
    Fragment,
    Clone,
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
enum LinkEvid {
    UnspecifiedOnly,
    Forbidden,
    Required,
}

/// NCBI reference (598d8ae6): c++/src/objects/seq/Seq_gap.cpp:176-194
/// ```c++
///     static const TGapTypeElem sc_gap_type_map[] = {
///         { "between-scaffolds", { CSeq_gap::eType_contig, eLinkEvid_Required } },
///         { "centromere", { CSeq_gap::eType_centromere, eLinkEvid_Forbidden } },
///         { "contamination", { eType_contamination, eLinkEvid_Required } },
///         { "heterochromatin", { CSeq_gap::eType_heterochromatin, eLinkEvid_Forbidden } },
///         { "repeat-between-scaffolds", { CSeq_gap::eType_repeat, eLinkEvid_Forbidden } },
///         { "repeat-within-scaffold", { CSeq_gap::eType_repeat, eLinkEvid_Required } },
///         { "short-arm", { CSeq_gap::eType_short_arm, eLinkEvid_Forbidden } },
///         { "telomere", { CSeq_gap::eType_telomere, eLinkEvid_Forbidden } },
///         { "unknown", { CSeq_gap::eType_unknown, eLinkEvid_UnspecifiedOnly } },
///         { "within-scaffold", { CSeq_gap::eType_scaffold, eLinkEvid_Forbidden } },
///     };
/// ```
/// `NameToGapTypeInfo` looks the name up after `CanonicalizeString` (Seq_gap.cpp:158-172).
fn gap_type_info(name: &[u8]) -> Option<(GapType, LinkEvid)> {
    Some(match canonicalize_string(name).as_slice() {
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
/// (`ENUM_METHOD_NAME(EType)()->NameToValue()`).
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:383-397
/// ```asn
/// Linkage-evidence ::= SEQUENCE {
///     type INTEGER {
///         paired-ends(0),
///         align-genus(1),
///         align-xgenus(2),
///         align-trnscpt(3),
///         within-clone(4),
///         clone-contig(5),
///         map(6),
///         strobe(7),
///         unspecified(8),
///         pcr(9),
///         proximity-ligation(10),
///         other(255)
///     }
/// }
/// ```
fn linkage_evidence_value(name: &[u8]) -> Option<u8> {
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
const LINKAGE_UNSPECIFIED: u8 = 8;

/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2128-2146
/// ```c++
/// string CFastaReader::CanonicalizeString(const TStr & sValue)
/// {
///     string newString;
///     newString.reserve(sValue.length());
///
///     ITERATE_0_IDX(ii, sValue.length()) {
///         const char ch = sValue[ii];
///         if( isupper(ch) ) {
///             newString.push_back(tolower(ch));
///         } else if( ch == ' ' || ch == '_' ) {
///             newString.push_back('-');
///         } else {
///             newString.push_back(ch);
///         }
///     }
///
///     return newString;
/// }
/// ```
fn canonicalize_string(value: &[u8]) -> Vec<u8> {
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

/// The title of a defline and its line number (`SLineTextAndLoc`).
struct LineText {
    text: Vec<u8>,
    line: u64,
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

/// The character types of a data line for nucleotide (`fAssumeNuc`) or protein
/// (`fAssumeProt`) input, without `fParseGaps` (no gap characters).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:854-935
/// ```c++
///         case 'A': case 'B': case 'C': case 'D':
///         case 'G': case 'H':
///         case 'K':
///         case 'M':
///         case 'R': case 'S': case 'T': case 'U': case 'V': case 'W':
///         case 'Y':
///             CloseGap(pos == 0);
///             m_SeqData[m_CurrentPos] = c;
///             CloseMask();
///             ++m_CurrentPos;
///             continue;
///         case 'a': case 'b': case 'c': case 'd':
///         case 'g': case 'h':
///         case 'k':
///         case 'm':
///         case 'r': case 's': case 't': case 'u': case 'v': case 'w':
///         case 'y':
///             char_type = eCharType_MaskedNonGap;
///             break;
///
///         case 'E': case 'F':
///         case 'I': case 'J':
///         case 'L':
///         case 'O': case 'P': case 'Q':
///         case 'Z':
///         case '*':
///             if( bIsNuc ) {
///                 char_type = eCharType_Bad;
///             } else {
///                 CloseGap(pos == 0);
///                 m_SeqData[m_CurrentPos] = c;
///                 CloseMask();
///                 ++m_CurrentPos;
///                 continue;
///             }
///             break;
///         case 'e': case 'f':
///         case 'i': case 'j':
///         case 'l':
///         case 'o': case 'p': case 'q':
///         case 'z':
///             char_type = (bIsNuc ? eCharType_Bad : eCharType_MaskedNonGap );
///             break;
///
///         case 'N':
///             char_type = ( bIsNuc && bAllowLetterGaps ?
///                      eCharType_Gap : eCharType_NormalNonGap );
///             break;
///         case 'n':
///             char_type = ( bIsNuc && bAllowLetterGaps ?
///                      eCharType_Gap : eCharType_MaskedNonGap );
///             break;
///
///         case 'X':
///             char_type = ( bIsNuc ? eCharType_Bad :
///                      bAllowLetterGaps ? eCharType_Gap :
///                      eCharType_NormalNonGap);
///             break;
///         case 'x':
///             char_type = ( bIsNuc ? eCharType_Bad :
///                      bAllowLetterGaps ? eCharType_Gap :
///                      eCharType_MaskedNonGap);
///             break;
///
///         case '-':
///             char_type = (
///                 bHyphensAreGaps ? eCharType_Gap :
///                 bHyphensIgnoreAndWarn ? eCharType_HyphenToIgnoreAndWarn :
///                 eCharType_NormalNonGap );
///             break;
///         case ';':
///             char_type = eCharType_Comment;
///             break;
///
///         case '\t': case '\n': case '\v': case '\f': case '\r': case ' ':
///             continue;
///
///         default:
///             char_type = eCharType_Bad;
///             break;
///         }
/// ```
/// `bIsNuc` is `fAssumeNuc` (the molecule is not set before `AssembleSeq` with
/// `fNoParseID`), `bAllowLetterGaps` and `bHyphensAreGaps` are false (no `fParseGaps`),
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
/// ```c++
///     for (size_t pos = 0;  pos < len_to_check;  ++pos) {
///         unsigned char c = s[pos];
///         if (s_ASCII_IsAlpha(c) ||  c == '*') {
///             ++good;
///             if( bIsNuc && s_ASCII_IsAmbigNuc(c) ) {
///                 ++ambig_nuc;
///             }
///         } else if( c == '-' ) {
///             if( ! bIgnoreHyphens ) {
///                 ++good;
///             }
///             // if bIgnoreHyphens == true, the "hyphens are ignored" warning
///             // will be triggered elsewhere
///         } else if (isspace(c)  ||  (c >= '0' && c <= '9')) {
///             // treat whitespace and digits as neutral
///         } else if (c == ';') {
///             break; // comment -- ignore rest of line
///         } else {
///             ++bad;
///         }
///     }
/// ```
const CHECK_CLASSES: [CheckClass; 256] = {
    let mut table = [CheckClass::Bad; 256];
    let mut index = 0;
    while index < 256 {
        let byte = index as u8;
        table[index] = if byte.is_ascii_alphabetic() || byte == b'*' {
            CheckClass::Good
        } else if byte == b'-' || input_space(byte) || byte.is_ascii_digit() {
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
    /// The local ID of the record (`GenerateID`), once its defline is read.
    local_id: Option<String>,
    /// Bad positions of a data line, reused.
    bad_positions: Vec<u32>,
}

/// Where a reader message goes: `LOG_POST_X(1, Warning << message)` with no listener.
pub(crate) type Warn<'a> = &'a mut dyn FnMut(&[u8]) -> std::io::Result<()>;

impl<R: Read> FastaReader<R> {
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:364-367
    /// ```c++
    ///     CRef<CSeqIdGenerator> idgen
    ///         (new CSeqIdGenerator(m_Config.GetLocalIdCounterInitValue(),
    ///                              m_Config.GetLocalIdPrefix()));
    ///     m_InputReader->SetIDGenerator(*idgen);
    /// ```
    pub(crate) fn new(lines: LineReader<R>, config: ReaderConfig) -> Self {
        Self {
            lines,
            config,
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:48-56
            // ```c++
            // CBlastInputSourceConfig::CBlastInputSourceConfig
            //     (const SDataLoaderConfig& dlconfig,
            //      objects::ENa_strand strand         /* = objects::eNa_strand_other */,
            //      bool lowercase                     /* = false */,
            //      bool believe_defline               /* = false */,
            //      TSeqRange range                    /* = TSeqRange() */,
            //      bool retrieve_seq_data             /* = true */,
            //      int local_id_counter               /* = 1 */,
            // ```
            id_counter: 1,
            seq_data: Vec::new(),
            gaps: Vec::new(),
            total_gap_length: 0,
            mask_range_start: None,
            mask: Vec::new(),
            titles: Vec::new(),
            title: Vec::new(),
            local_id: None,
            bad_positions: Vec::new(),
        }
    }

    /// `LineNumber()`.
    fn line_number(&self) -> u64 {
        self.lines.line_number()
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2196-2214
    /// ```c++
    /// void CFastaReader::PostWarning(
    ///             ILineErrorListener * pMessageListener,
    ///             EDiagSev _eSeverity, size_t _uLineNum, CTempString _MessageStrmOps, CObjReaderParseException::EErrCode _eErrCode, ILineError::EProblem _eProblem, CTempString _sFeature, CTempString _sQualName, CTempString _sQualValue) const
    /// {
    ///     if (find(m_ignorable.begin(), m_ignorable.end(), _eProblem) != m_ignorable.end())
    ///         // this is a problem that should be ignored
    ///         return;
    ///
    ///     string sSeqId = ( m_BestID ? m_BestID->AsFastaString() : kEmptyStr);
    ///     AutoPtr<CObjReaderLineException> pLineExpt(
    ///         CObjReaderLineException::Create(
    ///         (_eSeverity), static_cast<unsigned>(_uLineNum),
    ///         _MessageStrmOps,
    ///         (_eProblem),
    ///         sSeqId, (_sFeature),
    ///         (_sQualName), (_sQualValue),
    ///         _eErrCode) );
    ///     if ( ! pMessageListener && (_eSeverity) <= eDiag_Warning ) {
    ///         LOG_POST_X(1, Warning << pLineExpt->Message());
    /// ```
    /// The BLAST input source gives the reader no listener, so a warning is the message
    /// on a line of its own (`LOG_POST` writes no severity).
    fn post_warning(
        &self,
        problem: Problem,
        message: &[u8],
        warn: Warn<'_>,
    ) -> Result<(), ReadError> {
        if IGNORED_PROBLEMS.contains(&problem) {
            return Ok(());
        }
        let mut line = Vec::with_capacity(message.len() + 1);
        line.extend_from_slice(message);
        line.push(b'\n');
        warn(&line).map_err(ReadError::Write)
    }

    /// `FASTA_ERROR`: `PostWarning` with `eDiag_Error` and no listener throws
    /// `CObjReaderParseException(..., _eErrCode, _MessageStrmOps, _uLineNum, ...)`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2215-2218
    /// ```c++
    ///     } else if ( ! pMessageListener || ! pMessageListener->PutError( *pLineExpt ) )
    ///     {
    ///         throw CObjReaderParseException(DIAG_COMPILE_INFO, 0, _eErrCode, _MessageStrmOps, _uLineNum, _eSeverity);
    ///     }
    /// ```
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
    /// ```c++
    ///     virtual CRef<CSeq_entry> ReadOneSeq(ILineErrorListener * pMessageListener) {
    ///
    ///         const string line = NStr::TruncateSpaces_Unsafe(*++GetLineReader());
    ///         if ( !line.empty() && isalnum(line.data()[0]&0xff) ) {
    ///             try {
    ///                 CRef<CSeq_id> id(new CSeq_id(line, (CSeq_id::fParse_AnyRaw |
    /// 							CSeq_id::fParse_ValidLocal)));
    /// 		if (id->IsLocal()  &&  !NStr::StartsWith(line, "lcl|") ) {
    ///                     // Expected to throw an exception.
    ///                     id.Reset(new CSeq_id(line));
    /// 		}
    ///                 CRef<CBioseq> bioseq(x_CreateBioseq(id));
    ///                 CRef<CSeq_entry> retval(new CSeq_entry());
    ///                 retval->SetSeq(*bioseq);
    ///                 return retval;
    ///             } catch (const CSeqIdException& e) {
    ///                 if (NStr::Find(e.GetMsg(), "Malformatted ID") != NPOS) {
    ///                     // This is probably just plain fasta, so just
    ///                     // defer to CFastaReader
    ///                 } else {
    ///                     throw;
    ///                 }
    ///             } catch (const exception&) {
    ///                 throw;
    ///             } catch (...) {
    ///                 // in case of other exceptions, just defer to CFastaReader
    ///             }
    ///         } // end if ( !line.empty() ...
    ///
    ///         // If all fails, fall back to parent's implementation
    ///         GetLineReader().UngetLine();
    ///         return CFastaReader::ReadOneSeq(pMessageListener);
    ///     }
    /// ```
    /// NCBI reads a line that is a Seq-id as a sequence to fetch with its data loaders
    /// (GenBank over the network, or a BLAST database), which LOSAT does not do: such a
    /// line is an explicit rejection (`seq_id::reject_seq_id_line`). The rejection leaves
    /// the reader where NCBI's Seq-id record leaves it: the line is consumed (no
    /// `UngetLine` on that path) and the local-ID counter has not moved, so a further call
    /// reads from the next line (oracle BI net2, net5).
    pub(crate) fn read_one_seq(&mut self, warn: Warn<'_>) -> Result<FastaRecord, ReadError> {
        if self.config.data_loaders {
            let line = trim_input_space(self.lines.next_line()).to_vec();
            if !line.is_empty() && line[0].is_ascii_alphanumeric() {
                if let Some(rejection) = super::seq_id::reject_seq_id_line(&line, &self.config) {
                    return Err(ReadError::Unsupported(rejection));
                }
            }
            self.lines.unget_line();
        }
        self.read_fasta_seq(warn)
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:312-440
    /// ```c++
    /// CRef<CSeq_entry> CFastaReader::ReadOneSeq(ILineErrorListener * pMessageListener)
    /// {
    ///     m_CurrentSeq.Reset(new CBioseq);
    ///     // m_CurrentMask.Reset();
    ///     m_SeqData.erase();
    ///     m_Gaps.clear();
    ///     m_CurrentPos = 0;
    ///     m_BestID.Reset();
    ///     m_MaskRangeStart = kInvalidSeqPos;
    ///     if ( !TestFlag(fInSegSet) ) {
    ///         if (m_MaskVec  &&  m_NextMask.IsNull()) {
    ///             m_MaskVec->push_back(SaveMask());
    ///         }
    ///         m_CurrentMask.Reset(m_NextMask);
    ///         if (m_CurrentMask) {
    ///             m_CurrentMask->SetNull();
    ///         }
    ///         m_NextMask.Reset();
    ///         m_SegmentBase = 0;
    ///         m_Offset = 0;
    ///     }
    ///     m_CurrentGapLength = m_TotalGapLength = 0;
    ///     m_CurrentGapChar = '\0';
    ///     m_CurrentSeqTitles.clear();
    ///
    ///     bool need_defline = true;
    ///     CBadResiduesException::SBadResiduePositions bad_residue_positions;
    ///     while ( !GetLineReader().AtEOF() ) {
    ///         char c = GetLineReader().PeekChar();
    ///         if( LineNumber() % 10000 == 0 && LineNumber() != 0 ) {
    ///             FASTA_PROGRESS("Processing line " << LineNumber());
    ///         }
    ///         if (GetLineReader().AtEOF()) {
    ///             FASTA_ERROR(LineNumber(),
    ///                         "CFastaReader: Unexpected end-of-file around line " << LineNumber(),
    ///                         CObjReaderParseException::eEOF );
    ///             break;
    ///         }
    ///         if (c == '>' ) {
    ///             CTempString next_line = *++GetLineReader();
    ///             string strmodified;
    ///             if( NStr::StartsWith(next_line, ">?_") ) {
    ///                 CTempString tmp = next_line.substr(3);
    ///                 strmodified = ">";
    ///                 strmodified.append(tmp.data(), tmp.length());
    ///                 next_line = strmodified;
    ///             }
    ///             if( NStr::StartsWith(next_line, ">?") ) {
    ///                 // This is actually a data line. an assembly gap, in particular, which
    ///                 // we handle farther below
    ///                 GetLineReader().UngetLine();
    ///             } else {
    ///                 if(need_defline) {
    ///                     ParseDefLine(next_line, pMessageListener);
    ///                     need_defline = false;
    ///                     continue;
    ///                 } else {
    ///                     GetLineReader().UngetLine();
    ///                     // start of the next sequence
    ///                     break;
    ///                 }
    ///             }
    ///         }
    ///
    ///         CTempString line = NStr::TruncateSpaces_Unsafe(*++GetLineReader());
    ///
    ///         if (line.empty()) {
    ///             continue; // ignore lines containing only whitespace
    ///         }
    ///         c = line[0];
    ///
    ///         if (c == '!'  ||  c == '#' || c == ';') {
    ///             // no content, just a comment or blank line
    ///             continue;
    ///         } else if (need_defline) {
    ///             if (TestFlag(fDLOptional)) {
    ///                 ParseDefLine(">", pMessageListener);
    ///                 need_defline = false;
    ///             } else {
    ///                 ...
    ///             }
    ///         }
    ///
    ///         if ( !TestFlag(fNoSeqData) ) {
    ///             try {
    ///                 string strmodified;
    ///                 if( NStr::StartsWith(line, ">?_") ) {
    ///                     CTempString tmp = line.substr(3);
    ///                     strmodified = ">";
    ///                     strmodified.append(tmp.data(), tmp.length());
    ///                     line = strmodified;
    ///                 }
    ///                 ParseDataLine(line, pMessageListener);
    ///             } catch(const CBadResiduesException & e) {
    ///                 ...
    ///             }
    ///         }
    ///     }
    ///     ...
    ///     if (need_defline  &&  GetLineReader().AtEOF()) {
    ///         FASTA_ERROR(LineNumber(),
    ///             "CFastaReader: Expected defline around line " << LineNumber(),
    ///             CObjReaderParseException::eEOF);
    ///     }
    ///
    ///     AssembleSeq(pMessageListener);
    /// ```
    /// `FASTA_PROGRESS` posts to a listener only, and there is none. `CBadResiduesException`
    /// comes only with `fValidate`, which is not set.
    fn read_fasta_seq(&mut self, warn: Warn<'_>) -> Result<FastaRecord, ReadError> {
        self.seq_data.clear();
        self.gaps.clear();
        self.mask_range_start = None;
        self.mask.clear();
        self.total_gap_length = 0;
        self.titles.clear();
        self.title.clear();
        self.local_id = None;

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
                    self.parse_def_line(next_line);
                    self.lines.put_line(raw);
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
            let line = trim_input_space(&raw);
            if line.is_empty() || matches!(line[0], b'!' | b'#' | b';') {
                self.lines.put_line(raw);
                continue;
            }
            if need_defline {
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
        self.assemble_seq(warn)
    }

    /// The title of a defline (`ParseDefLine` with `fNoParseID`), and the generated local
    /// ID (`PostProcessIDs` without IDs).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:146-226
    /// ```c++
    ///     const size_t len = defline.length();
    ///     if (len <= 1 ||
    ///         NStr::IsBlank(defline.substr(1))) {
    ///         return;
    ///     }
    ///
    ///     if (defline[0] != '>') {
    ///         NCBI_THROW2(CObjReaderParseException, eFormat,
    ///             "Invalid defline. First character is not '>'", 0);
    ///     }
    ///
    ///     // ignore spaces between '>' and the sequence ID
    ///     size_t start;
    ///     for(start = 1 ; start < len; ++start ) {
    ///         if( ! isspace(defline[start]) ) {
    ///             break;
    ///         }
    ///     }
    ///
    ///     size_t pos;
    ///     size_t title_start = NPOS;
    ///     if ((fFastaFlags & CFastaReader::fNoParseID)) {
    ///         title_start = start;
    ///     }
    ///     ...
    ///     // trim leading whitespace from title (is this appropriate?)
    ///     while (title_start < len
    ///         &&  isspace((unsigned char)defline[title_start])) {
    ///         ++title_start;
    ///     }
    ///
    ///     if (title_start < len) {
    ///         for (pos = title_start + 1;  pos < len;  ++pos) {
    ///             if ((unsigned char)defline[pos] < ' ') {
    ///             break;
    ///             }
    ///         }
    ///         // Parse the title elsewhere - after the molecule has been deduced
    ///         data.titles.push_back(
    ///             SLineTextAndLoc(
    ///                 defline.substr(title_start, pos - title_start), lineNumber));
    ///     }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:566-574
    /// ```c++
    /// void CFastaReader::ParseDefLine(const TStr& s, ILineErrorListener * pMessageListener)
    /// {
    ///     SDefLineParseInfo parseInfo;
    ///     x_SetDeflineParseInfo(parseInfo);
    ///
    ///     CFastaDeflineReader::SDeflineData data;
    ///     CFastaDeflineReader::ParseDefline(s, parseInfo, data, pMessageListener, m_fIdCheck);
    ///
    ///     m_CurrentSeqTitles = std::move(data.titles);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:619-631
    /// ```c++
    /// void CFastaReader::PostProcessIDs(
    ///     const CBioseq::TId& defline_ids,
    ///     const string& /*defline*/,
    ///     const bool has_range,
    ///     const TSeqPos range_start,
    ///     const TSeqPos range_end)
    /// {
    ///     if (defline_ids.empty()) {
    ///         GenerateID();
    ///     }
    /// ```
    /// The `CFastaReader` user object with the raw defline is not output.
    fn parse_def_line(&mut self, defline: &[u8]) {
        self.titles.clear();
        let len = defline.len();
        if len > 1 && !defline[1..].iter().all(|&byte| input_space(byte)) {
            let start = 1 + defline[1..]
                .iter()
                .position(|&byte| !input_space(byte))
                .unwrap_or(len - 1);
            let title_start = start
                + defline[start..]
                    .iter()
                    .position(|&byte| !input_space(byte))
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
                    line: self.line_number(),
                });
            }
        }
        self.local_id = Some(self.generate_id());
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:466-483
    /// ```c++
    /// CRef<CSeq_id> CSeqIdGenerator::GenerateID(const string& defline, const bool advance)
    /// {
    ///     CRef<CSeq_id> seq_id(new CSeq_id);
    ///     auto n = m_Counter.load();
    ///     if (advance)
    ///         m_Counter++;
    ///
    ///     if (m_Prefix.empty()  &&  m_Suffix.empty()) {
    ///         seq_id->SetLocal().SetId(n);
    ///     } else {
    ///         string& id = seq_id->SetLocal().SetStr();
    ///         id.reserve(128);
    ///         id += m_Prefix;
    ///         id += NStr::IntToString(n);
    ///         id += m_Suffix;
    ///     }
    ///     return seq_id;
    /// }
    /// ```
    fn generate_id(&mut self) -> String {
        let n = self.id_counter;
        self.id_counter += 1;
        format!("{}{n}", self.config.id_prefix)
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:713-772
    /// ```c++
    ///     // make sure the first data line has at least SOME resemblance to
    ///     // actual sequence data.
    ///     if (TestFlag(fSkipCheck)  ||  ! m_SeqData.empty() ) {
    ///         return;
    ///     }
    ///     const bool bIgnoreHyphens = TestFlag(fHyphensIgnoreAndWarn);
    ///     size_t good = 0, bad = 0;
    ///     // in case the data has huge sequences all on the first line we do need
    ///     // a cutoff and "70" seems reasonable since it's the default width of
    ///     // CFastaOstream (as of 2017-03-09)
    ///     size_t len_to_check = min(s.length(),
    ///                               static_cast<size_t>(70));
    ///     const bool bIsNuc = (
    ///         ( TestFlag(fAssumeNuc) && TestFlag(fForceType) ) ||
    ///         ( m_CurrentSeq && m_CurrentSeq->IsSetInst() &&
    ///           m_CurrentSeq->GetInst().IsSetMol() &&  m_CurrentSeq->IsNa() ) );
    ///     size_t ambig_nuc = 0;
    ///     for (size_t pos = 0;  pos < len_to_check;  ++pos) {
    ///         unsigned char c = s[pos];
    ///         if (s_ASCII_IsAlpha(c) ||  c == '*') {
    ///             ++good;
    ///             if( bIsNuc && s_ASCII_IsAmbigNuc(c) ) {
    ///                 ++ambig_nuc;
    ///             }
    ///         } else if( c == '-' ) {
    ///             if( ! bIgnoreHyphens ) {
    ///                 ++good;
    ///             }
    ///             // if bIgnoreHyphens == true, the "hyphens are ignored" warning
    ///             // will be triggered elsewhere
    ///         } else if (isspace(c)  ||  (c >= '0' && c <= '9')) {
    ///             // treat whitespace and digits as neutral
    ///         } else if (c == ';') {
    ///             break; // comment -- ignore rest of line
    ///         } else {
    ///             ++bad;
    ///         }
    ///     }
    ///     if (bad >= good / 3  &&
    ///         (len_to_check > 3  ||  good == 0  ||  bad > good))
    ///     {
    ///         FASTA_ERROR( LineNumber(),
    ///             "CFastaReader: Near line " << LineNumber()
    ///             << ", there's a line that doesn't look like plausible data, "
    ///             "but it's not marked as defline or comment.",
    ///             CObjReaderParseException::eFormat);
    ///     }
    ///     // warn if more than a certain percentage is ambiguous nucleotides
    ///     const static size_t kWarnPercentAmbiguous = 40; // e.g. "40" means "40%"
    ///     const size_t percent_ambig = (good == 0)?100:((ambig_nuc * 100) / good);
    ///     if( len_to_check > 3 && percent_ambig > kWarnPercentAmbiguous ) {
    ///         FASTA_WARNING(LineNumber(),
    ///             "FASTA-Reader: Start of first data line in seq is about "
    ///             << percent_ambig << "% ambiguous nucleotides (shouldn't be over "
    ///             << kWarnPercentAmbiguous << "%)",
    ///             ILineError::eProblem_TooManyAmbiguousResidues,
    ///             "first data line");
    ///     }
    /// ```
    /// `bIsNuc` is false (no `fForceType`; the molecule is set in `AssembleSeq`), so no
    /// letter counts as ambiguous, and the warning (an ignored problem) needs `good == 0`,
    /// which the error before has already thrown.
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
    /// character switch is quoted at `char_types`)
    /// ```c++
    /// void CFastaReader::ParseDataLine(
    ///     const TStr& s, ILineErrorListener * pMessageListener)
    /// {
    ///     if( NStr::StartsWith(s, ">?") ) {
    ///         ParseGapLine(s, pMessageListener);
    ///         return;
    ///     }
    ///
    ///     CheckDataLine(s, pMessageListener);
    ///     ...
    ///     switch(char_type) {
    ///     case eCharType_NormalNonGap:
    ///         CloseGap(pos == 0);
    ///         m_SeqData[m_CurrentPos] = c;
    ///         CloseMask();
    ///         ++m_CurrentPos;
    ///         break;
    ///     case eCharType_MaskedNonGap:
    ///         CloseGap(pos == 0);
    ///         m_SeqData[m_CurrentPos] = s_ASCII_MustBeLowerToUpper(c);
    ///         OpenMask();
    ///         ++m_CurrentPos;
    ///         break;
    ///     ...
    ///     case eCharType_HyphenToIgnoreAndWarn:
    ///         bIgnorableHyphenSeen = true;
    ///         break;
    ///     case eCharType_Comment:
    ///         // artificially advance pos to the end to break the pos loop
    ///         pos = s_len;
    ///         break;
    ///     case eCharType_Bad:
    ///         if( bad_pos_line_num < 0 ) {
    ///             bad_pos_line_num = LineNumber();
    ///         }
    ///         bad_pos_vec.push_back(pos);
    ///         break;
    ///     default:
    ///         _TROUBLE;
    ///     }
    /// }
    ///
    /// m_SeqData.resize(m_CurrentPos);
    ///
    /// if( bIgnorableHyphenSeen ) {
    ///     _ASSERT( bHyphensIgnoreAndWarn );
    ///     FASTA_WARNING_EX(LineNumber(),
    ///         "CFastaReader: Hyphens are invalid and will be ignored around line " << LineNumber(),
    ///         ILineError::eProblem_IgnoredResidue,
    ///         kEmptyStr, kEmptyStr, "-" );
    /// }
    ///
    /// // before throwing, be sure that we're in a valid state so that callers can
    /// // parse multiple lines and get the invalid residues in all of them.
    ///
    /// if( ! bad_pos_vec.empty() ) {
    ///     if (TestFlag(fValidate)) {
    ///         ...
    ///     } else {
    ///         stringstream warn_strm;
    ///         warn_strm << "FASTA-Reader: Ignoring invalid " << x_NucOrProt()
    ///             << "residues at position(s): ";
    ///         CBadResiduesException::SBadResiduePositions(
    ///             m_BestID, bad_pos_vec, bad_pos_line_num ).ConvertBadIndexesToString(warn_strm);
    ///
    ///         FASTA_WARNING(0,
    ///             warn_strm.str(),
    ///             ILineError::eProblem_InvalidResidue,
    ///             kEmptyStr );
    ///     }
    /// }
    /// ```
    /// `CloseGap` does nothing (no gap characters without `fParseGaps`, and
    /// `CCustomizedFastaReader::x_CloseGap` is empty); `x_NucOrProt()` is empty, because
    /// the molecule is set only in `AssembleSeq`. The fast copy of `fSkipCheck` is not
    /// taken (no `fSkipCheck`).
    fn parse_data_line(&mut self, s: &[u8], warn: Warn<'_>) -> Result<(), ReadError> {
        if s.starts_with(b">?") {
            return self.parse_gap_line(s, warn);
        }
        self.check_data_line(s)?;
        let (types, residues, masked) = if self.config.protein {
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
                    } else if self.config.protein {
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
                    } else if self.config.protein {
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
            let mut message = b"FASTA-Reader: Ignoring invalid residues at position(s): ".to_vec();
            convert_bad_indexes_to_string(&mut message, self.line_number(), &self.bad_positions);
            self.post_warning(Problem::InvalidResidue, &message, warn)?;
        }
        Ok(())
    }

    /// LOSAT's limit: NCBI keeps positions in `TSeqPos` (32 bits) and BLAST's sequence
    /// lengths in `Int4`; a record longer than that is an explicit rejection.
    fn check_length(&self) -> Result<(), ReadError> {
        if self.seq_data.len() as u64 + self.total_gap_length as u64 > i32::MAX as u64 {
            return Err(ReadError::Unsupported(anyhow::anyhow!(
                "a {} record longer than 2147483647 letters is not supported by LOSAT's {}",
                self.config.role,
                self.config.program
            )));
        }
        Ok(())
    }

    /// `GetCurrentPos(ePosWithGapsAndSegs)`: the raw position plus the gaps before it
    /// (`m_SegmentBase` is 0 without `fInSegSet`).
    ///
    /// NCBI reference (598d8ae6): c++/include/objtools/readers/fasta.hpp:502-516
    /// ```c++
    /// inline
    /// TSeqPos CFastaReader::GetCurrentPos(EPosType pos_type)
    /// {
    ///     TSeqPos pos = m_CurrentPos;
    ///     switch (pos_type) {
    ///     case ePosWithGapsAndSegs:
    ///         return pos + m_SegmentBase + m_TotalGapLength;
    ///     case ePosWithGaps:
    ///         return pos + m_TotalGapLength;
    ///     case eRawPos:
    ///         return pos;
    /// ```
    fn current_pos_with_gaps(&self) -> u32 {
        self.seq_data.len() as u32 + self.total_gap_length
    }

    /// The mask is always recorded (LOSAT marks the masked residues by case, and the
    /// programs apply the mask only with `-lcase_masking`, where NCBI's
    /// `x_FastaToSeqLoc` calls `SaveMask` before the record is read).
    ///
    /// NCBI reference (598d8ae6): c++/include/objtools/readers/fasta.hpp:494-499
    /// ```c++
    /// inline
    /// void CFastaReader::OpenMask()
    /// {
    ///     if (m_MaskRangeStart == kInvalidSeqPos  &&  m_CurrentMask.NotEmpty()) {
    ///         x_OpenMask();
    ///     }
    /// }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1079-1083
    /// ```c++
    /// void CFastaReader::x_OpenMask(void)
    /// {
    ///     _ASSERT(m_MaskRangeStart == kInvalidSeqPos);
    ///     m_MaskRangeStart = GetCurrentPos(ePosWithGapsAndSegs);
    /// }
    /// ```
    /// `x_OpenMask` runs before `++m_CurrentPos`, so the start is the masked residue.
    fn open_mask(&mut self) {
        if self.mask_range_start.is_none() {
            self.mask_range_start = Some(self.current_pos_with_gaps() - 1);
        }
    }

    /// NCBI reference (598d8ae6): c++/include/objtools/readers/fasta.hpp:291-292
    /// ```c++
    ///     void CloseMask(void)
    ///         { if (m_MaskRangeStart != kInvalidSeqPos) { x_CloseMask(); } }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1085-1092
    /// ```c++
    /// void CFastaReader::x_CloseMask(void)
    /// {
    ///     _ASSERT(m_MaskRangeStart != kInvalidSeqPos);
    ///     m_CurrentMask->SetPacked_int().AddInterval
    ///         (GetBestID(), m_MaskRangeStart, GetCurrentPos(ePosWithGapsAndSegs) - 1,
    ///          eNa_strand_plus);
    ///     m_MaskRangeStart = kInvalidSeqPos;
    /// }
    /// ```
    /// In `ParseDataLine` `CloseMask` runs after the residue is stored and before
    /// `++m_CurrentPos`: the interval ends before the residue that closes it (LOSAT has
    /// already pushed that residue, hence the `- 2`). `AssembleSeq` closes it at the end.
    fn close_mask(&mut self) {
        if let Some(start) = self.mask_range_start.take() {
            self.mask.push((start, self.current_pos_with_gaps() - 2));
        }
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1094-1349
    /// ```c++
    /// bool CFastaReader::ParseGapLine(
    ///     const TStr& line, ILineErrorListener * pMessageListener)
    /// {
    ///     _ASSERT( NStr::StartsWith(line, ">?") );
    ///
    ///     // just in case there's a gap before this one,
    ///     // even though somewhat unusual
    ///     CloseGap();
    ///
    ///     // sRemainingLine will hold the part of the line left to parse
    ///     TStr sRemainingLine = line.substr(2);
    ///     NStr::TruncateSpacesInPlace(sRemainingLine);
    ///
    ///     const TSeqPos uPos = GetCurrentPos(eRawPos);
    ///
    ///     // check if size is unknown
    ///     SGap::EKnownSize eIsKnown = SGap::eKnownSize_Yes;
    ///     if( NStr::StartsWith(sRemainingLine, "unk") ) {
    ///         eIsKnown = SGap::eKnownSize_No;
    ///         sRemainingLine = sRemainingLine.substr(3);
    ///         NStr::TruncateSpacesInPlace(sRemainingLine, NStr::eTrunc_Begin);
    ///     }
    ///
    ///     // extract the gap size
    ///     TSeqPos uGapSize = 0;
    ///     {
    ///         // find how many digits in number
    ///         TStr::size_type uNumDigits = 0;
    ///         while( uNumDigits != sRemainingLine.size() &&
    ///                ::isdigit(sRemainingLine[uNumDigits]) )
    ///         {
    ///             ++uNumDigits;
    ///         }
    ///         TStr sDigits = sRemainingLine.substr(
    ///             0, uNumDigits);
    ///         uGapSize = NStr::StringToUInt(sDigits, NStr::fConvErr_NoThrow);
    ///         if( uGapSize <= 0 ) {
    ///             FASTA_WARNING(LineNumber(),
    ///                 "CFastaReader: Bad gap size at line " << LineNumber(),
    ///                 ILineError::eProblem_NonPositiveLength,
    ///                 "gapline" );
    ///             // try to continue the best we can
    ///             uGapSize = 1;
    ///         }
    ///         sRemainingLine = sRemainingLine.substr(sDigits.length());
    ///         NStr::TruncateSpacesInPlace(sRemainingLine, NStr::eTrunc_Begin);
    ///     }
    ///
    ///     // extract the raw key-value pairs for the gap
    ///     typedef multimap<TStr, TStr> TModKeyValueMultiMap;
    ///     TModKeyValueMultiMap modKeyValueMultiMap;
    ///     while( ! sRemainingLine.empty() ) {
    ///         TStr::size_type uOpenBracketPos = TStr::npos;
    ///         if ( NStr::StartsWith(sRemainingLine, "[") ) {
    ///             uOpenBracketPos = 0;
    ///         }
    ///         TStr::size_type uPosOfEqualSign = TStr::npos;
    ///         if( uOpenBracketPos != TStr::npos ) {
    ///             // uses "1" to skip the '['
    ///             uPosOfEqualSign = sRemainingLine.find('=', uOpenBracketPos + 1);
    ///         }
    ///         TStr::size_type uCloseBracketPos = TStr::npos;
    ///         if( uPosOfEqualSign != TStr::npos ) {
    ///             uCloseBracketPos = sRemainingLine.find(']', uPosOfEqualSign + 1);
    ///         }
    ///         if( uCloseBracketPos == TStr::npos )
    ///         {
    ///             FASTA_WARNING(LineNumber(),
    ///                 "CFastaReader: Problem parsing gap mods at line "
    ///                 << LineNumber(),
    ///                 ILineError::eProblem_ParsingModifiers,
    ///                 "gapline" );
    ///             break; // give up on mod-parsing
    ///         }
    ///
    ///         // extract the key and the value
    ///         TStr sKey = NStr::TruncateSpaces_Unsafe(
    ///             sRemainingLine.substr(uOpenBracketPos + 1,
    ///                 (uPosOfEqualSign - uOpenBracketPos - 1) ) );
    ///         TStr sValue = NStr::TruncateSpaces_Unsafe(
    ///             sRemainingLine.substr(uPosOfEqualSign + 1,
    ///                 uCloseBracketPos - uPosOfEqualSign - 1) );
    ///
    ///         // remember what we saw
    ///         modKeyValueMultiMap.insert(
    ///             TModKeyValueMultiMap::value_type(sKey, sValue) );
    ///
    ///         // prepare for the next loop around
    ///         sRemainingLine = sRemainingLine.substr(uCloseBracketPos + 1);
    ///         NStr::TruncateSpacesInPlace(sRemainingLine, NStr::eTrunc_Begin);
    ///     }
    /// ```
    /// The rest of the function (the modifiers, quoted at `gap_mod_warnings`) only warns:
    /// without `fParseGaps` `AssembleSeq` turns each gap into a run of `N` or `X`.
    fn parse_gap_line(&mut self, line: &[u8], warn: Warn<'_>) -> Result<(), ReadError> {
        let mut remaining = trim_input_space(&line[2..]);
        let position = self.seq_data.len() as u32;
        if let Some(rest) = remaining.strip_prefix(b"unk") {
            remaining = trim_input_start(rest);
        }
        let digits = remaining
            .iter()
            .position(|byte| !byte.is_ascii_digit())
            .unwrap_or(remaining.len());
        // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp (NStr::StringToUInt with
        // fConvErr_NoThrow): an empty string or a value over UINT_MAX converts to 0.
        let mut gap_size = std::str::from_utf8(&remaining[..digits])
            .ok()
            .and_then(|text| text.parse::<u32>().ok())
            .unwrap_or(0);
        if gap_size == 0 {
            let line = self.line_number();
            self.post_warning(
                Problem::NonPositiveLength,
                format!("CFastaReader: Bad gap size at line {line}").as_bytes(),
                warn,
            )?;
            gap_size = 1;
        }
        remaining = trim_input_start(&remaining[digits..]);
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
                let line = self.line_number();
                self.post_warning(
                    Problem::ParsingModifiers,
                    format!("CFastaReader: Problem parsing gap mods at line {line}").as_bytes(),
                    warn,
                )?;
                break;
            };
            let key = trim_input_space(&remaining[1..equal]).to_vec();
            let value = trim_input_space(&remaining[equal + 1..close]).to_vec();
            mods.push((key, value));
            remaining = trim_input_start(&remaining[close + 1..]);
        }
        // `multimap` iterates in key order (stable for equal keys).
        mods.sort_by(|left, right| left.0.cmp(&right.0));
        self.gap_mod_warnings(&mods, warn)?;
        self.gaps.push(Gap {
            position,
            length: gap_size,
        });
        self.total_gap_length = self.total_gap_length.saturating_add(gap_size);
        self.check_length()?;
        Ok(())
    }

    /// The warnings of `ParseGapLine` about the gap modifiers. Without `fParseGaps` the
    /// gap type and linkage evidence are not kept (`AssembleSeq` warns that the
    /// modifiers are ignored, an ignored problem).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1190-1340
    /// ```c++
    ///     // remember if there is a gap-type conflict
    ///     bool bConflictingGapTypes = false;
    ///     // extract the mods, if any
    ///     SGap::TNullableGapType pGapType;
    ///     CSeq_gap::ELinkEvid eLinkEvid = CSeq_gap::eLinkEvid_UnspecifiedOnly;
    ///     set<CLinkage_evidence::EType> setOfLinkageEvidence;
    ///
    ///     if (m_gap_type && modKeyValueMultiMap.empty()) // fall back to default values coming from caller
    ///     {
    ///         ...
    ///     }
    ///
    ///     ITERATE(TModKeyValueMultiMap, modKeyValue_it, modKeyValueMultiMap) {
    ///         const TStr & sKey   = modKeyValue_it->first;
    ///         const TStr & sValue = modKeyValue_it->second;
    ///
    ///         string canonicalKey = CanonicalizeString(sKey);
    ///         if(  canonicalKey == "gap-type") {
    ///
    ///             const CSeq_gap::SGapTypeInfo *pGapTypeInfo = CSeq_gap::NameToGapTypeInfo(sValue);
    ///             if( pGapTypeInfo ) {
    ///                  CSeq_gap::EType eGapType = pGapTypeInfo->m_eType;
    ///
    ///                 if( ! pGapType ) {
    ///                     pGapType.Reset( new SGap::TGapTypeObj(eGapType) );
    ///                     eLinkEvid = pGapTypeInfo->m_eLinkEvid;
    ///                 } else if( eGapType != *pGapType ) {
    ///                     // check if pGapType already set and different
    ///                     bConflictingGapTypes = true;
    ///                 }
    ///             } else {
    ///                 FASTA_WARNING_EX(
    ///                     LineNumber(),
    ///                     "Unknown gap-type: " << sValue,
    ///                     ILineError::eProblem_ParsingModifiers,
    ///                     "gapline",
    ///                     "gap-type",
    ///                     sValue );
    ///             }
    ///             continue;
    ///
    ///         }
    ///
    ///         if( canonicalKey == "linkage-evidence") {
    ///             // could be semi-colon separated
    ///             vector<CTempString> arrLinkageEvidences;
    ///             NStr::Split(sValue, ";", arrLinkageEvidences,
    ///                         NStr::fSplit_MergeDelimiters | NStr::fSplit_Truncate);
    ///
    ///             ITERATE(vector<CTempString>, link_evid_it, arrLinkageEvidences) {
    ///                 CTempString sLinkEvid = *link_evid_it;
    ///                 CEnumeratedTypeValues::TNameToValue::const_iterator find_iter =
    ///                     linkage_evidence_to_value_map.find(CanonicalizeString(sLinkEvid));
    ///                 if( find_iter != linkage_evidence_to_value_map.end() ) {
    ///                     setOfLinkageEvidence.insert(
    ///                         static_cast<CLinkage_evidence::EType>(
    ///                         find_iter->second));
    ///                 } else {
    ///                     FASTA_WARNING_EX(
    ///                         LineNumber(),
    ///                         "Unknown linkage-evidence: " << sValue,
    ///                         ILineError::eProblem_ParsingModifiers,
    ///                         "gapline",
    ///                         "linkage-evidence",
    ///                         sValue );
    ///                 }
    ///             }
    ///             continue;
    ///         }
    ///
    ///         // unknown mod.
    ///         FASTA_WARNING_EX(
    ///             LineNumber(),
    ///             "Unknown gap modifier name(s): " << sKey,
    ///             ILineError::eProblem_UnrecognizedQualifierName,
    ///             "gapline", sKey, kEmptyStr );
    ///     }
    ///
    ///     if( bConflictingGapTypes ) {
    ///         FASTA_WARNING_EX(LineNumber(),
    ///             "There were conflicting gap-types around line " << LineNumber(),
    ///             ILineError::eProblem_ContradictoryModifiers,
    ///             "gapline", "gap-type", kEmptyStr );
    ///     }
    ///
    ///     // check validation beyond basic parsing problems
    ///
    ///     // if no gap-type set (but linkage-evidence explicitly set, use "unknown")
    ///     if( ! pGapType && ! setOfLinkageEvidence.empty() ) {
    ///         pGapType.Reset( new SGap::TGapTypeObj(CSeq_gap::eType_unknown) );
    ///     }
    ///
    ///     // check if linkage-evidence(s) compatible with gap-type
    ///     switch( eLinkEvid ) {
    ///     case CSeq_gap::eLinkEvid_UnspecifiedOnly:
    ///         if( setOfLinkageEvidence.empty() ) {
    ///             if( pGapType ) {
    ///                 // silently add the required "unspecified"
    ///                 setOfLinkageEvidence.insert(CLinkage_evidence::eType_unspecified);
    ///             }
    ///         } else if( setOfLinkageEvidence.size() > 1 ||
    ///             *setOfLinkageEvidence.begin() != CLinkage_evidence::eType_unspecified )
    ///         {
    ///             // only "unspecified" is allowed
    ///             FASTA_WARNING(
    ///                 LineNumber(),
    ///                 "FASTA-Reader: Unknown gap-type can have linkage-evidence "
    ///                     "of type 'unspecified' only.",
    ///                 ILineError::eProblem_ExtraModifierFound,
    ///                 "gapline");
    ///             setOfLinkageEvidence.clear();
    ///             setOfLinkageEvidence.insert(CLinkage_evidence::eType_unspecified);
    ///         }
    ///         break;
    ///     case CSeq_gap::eLinkEvid_Forbidden:
    ///         if( ! setOfLinkageEvidence.empty() ) {
    ///             FASTA_WARNING(LineNumber(),
    ///                 "FASTA-Reader: This gap-type cannot have any "
    ///                 "linkage-evidence specified, so any will be ignored.",
    ///                 ILineError::eProblem_ModifierFoundButNoneExpected,
    ///                 "gapline" );
    ///             setOfLinkageEvidence.clear();
    ///         }
    ///         break;
    ///     case CSeq_gap::eLinkEvid_Required:
    ///         if( setOfLinkageEvidence.empty() ) {
    ///             setOfLinkageEvidence.insert(CLinkage_evidence::eType_unspecified);
    ///         }
    ///         if( setOfLinkageEvidence.size() == 1 &&
    ///             *setOfLinkageEvidence.begin() == CLinkage_evidence::eType_unspecified)
    ///         {
    ///             FASTA_WARNING(LineNumber(),
    ///                 "CFastaReader: This gap-type should have at least one "
    ///                 "specified linkage-evidence.",
    ///                 ILineError::eProblem_ExpectedModifierMissing,
    ///                 "gapline" );
    ///         }
    ///         break;
    ///         // intentionally omitted "default:" so a compiler warning will
    ///         // hopefully let us know if we've forgotten a case
    ///     }
    /// ```
    /// The caller's default gap type (`m_gap_type`) is not set by the BLAST input source.
    /// `NStr::Split` with `fSplit_MergeDelimiters | fSplit_Truncate` drops empty pieces.
    fn gap_mod_warnings(
        &self,
        mods: &[(Vec<u8>, Vec<u8>)],
        warn: Warn<'_>,
    ) -> Result<(), ReadError> {
        let line = self.line_number();
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
                        self.post_warning(Problem::ParsingModifiers, &message, warn)?;
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
                            self.post_warning(Problem::ParsingModifiers, &message, warn)?;
                        }
                    }
                }
                continue;
            }
            let mut message = b"Unknown gap modifier name(s): ".to_vec();
            message.extend_from_slice(key);
            self.post_warning(Problem::UnrecognizedQualifierName, &message, warn)?;
        }
        if conflicting {
            self.post_warning(
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
                    self.post_warning(
                        Problem::ExtraModifierFound,
                        b"FASTA-Reader: Unknown gap-type can have linkage-evidence of type 'unspecified' only.",
                        warn,
                    )?;
                }
            }
            LinkEvid::Forbidden => {
                if !evidence.is_empty() {
                    self.post_warning(
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
                    self.post_warning(
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

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1351-1448
    /// ```c++
    /// void CFastaReader::AssembleSeq(ILineErrorListener * pMessageListener)
    /// {
    ///     CSeq_inst& inst = m_CurrentSeq->SetInst();
    ///
    ///     CloseGap();
    ///     CloseMask();
    ///     if (TestFlag(fInSegSet)) {
    ///         m_SegmentBase += GetCurrentPos(ePosWithGaps);
    ///     }
    ///     AssignMolType(pMessageListener);
    ///
    ///     // apply source mods *after* figuring out mol type
    ///     ITERATE(vector<SLineTextAndLoc>, title_ci, m_CurrentSeqTitles) {
    ///         ParseTitle(*title_ci, pMessageListener);
    ///     }
    ///     m_CurrentSeqTitles.clear();
    ///     ...
    ///     if ( !TestFlag(fParseGaps)  &&  m_TotalGapLength > 0 ) {
    ///         // Encountered >? lines; substitute runs of Ns or Xs as appropriate.
    ///         string    new_data;
    ///         char      gap_char(inst.IsAa() ? 'X' : 'N');
    ///         SIZE_TYPE pos = 0;
    ///         new_data.reserve(GetCurrentPos(ePosWithGaps));
    ///         ITERATE (TGaps, it, m_Gaps) {
    ///             ...
    ///             if ((*it)->m_uPos > pos) {
    ///                 new_data.append(m_SeqData, pos, (*it)->m_uPos - pos);
    ///                 pos = (*it)->m_uPos;
    ///             }
    ///             new_data.append((*it)->m_uLen, gap_char);
    ///         }
    ///         if (m_CurrentPos > pos) {
    ///             new_data.append(m_SeqData, pos, m_CurrentPos - pos);
    ///         }
    ///         swap(m_SeqData, new_data);
    ///         m_Gaps.clear();
    ///         m_CurrentPos += m_TotalGapLength;
    ///         m_TotalGapLength = 0;
    ///         m_CurrentGapChar = '\0';
    ///     }
    ///
    ///     if (m_Gaps.empty() && m_SeqData.empty()) {
    ///
    ///         _ASSERT(m_TotalGapLength == 0);
    ///             inst.SetLength(0);
    ///             inst.SetRepr(CSeq_inst::eRepr_virtual);
    ///             // empty sequence triggers warning if seq data was expected
    ///             if( ! TestFlag(fDisableNoResidues) &&
    ///                 ! TestFlag(fNoSeqData) ) {
    ///                 FASTA_ERROR(LineNumber(),
    ///                     "FASTA-Reader: No residues given",
    ///                     CObjReaderParseException::eNoResidues);
    ///             }
    ///     }
    ///     else
    ///     if (m_Gaps.empty() && TestFlag(fNoSplit)) {
    ///         inst.SetLength(GetCurrentPos(eRawPos));
    ///         inst.SetRepr(CSeq_inst::eRepr_raw);
    ///         CRef<CSeq_data> data(new CSeq_data(m_SeqData, format));
    ///         if ( !TestFlag(fLeaveAsText) ) {
    ///             CSeqportUtil::Pack(data, inst.GetLength());
    ///         }
    ///         inst.SetSeq_data(*data);
    /// ```
    /// The molecule is the input's (`CCustomizedFastaReader::AssignMolType` below
    /// `m_SeqLenThreshold`, which the BLAST input source sets to `UINT_MAX`), so the gap
    /// residue is `N` for nucleotides and `X` for proteins, and `x_NucOrProt` is not used
    /// after this point. Without `fNoSplit` (`BLASTINPUT_GEN_DELTA_SEQ` set) the gaps have
    /// already become residues, so the record is the same sequence in delta pieces.
    /// `CSeqportUtil::Pack` packs nucleotides into `ncbi2na`/`ncbi4na`, where `U` and `T`
    /// are the same code, so LOSAT stores `U` as `T`. The residues in the lowercase mask
    /// are returned in lower case (LOSAT's representation of the mask).
    fn assemble_seq(&mut self, warn: Warn<'_>) -> Result<FastaRecord, ReadError> {
        // CloseMask at the end of the record: the interval ends at the last position.
        if let Some(start) = self.mask_range_start.take() {
            self.mask.push((start, self.current_pos_with_gaps() - 1));
        }
        let mut titles = std::mem::take(&mut self.titles);
        for title in titles.drain(..) {
            self.parse_title(title, warn)?;
        }
        self.titles = titles;
        let gap_char = if self.config.protein { b'X' } else { b'N' };
        let mut sequence = if self.total_gap_length > 0 {
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
            std::mem::take(&mut self.seq_data)
        };
        if !self.config.protein {
            // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1444-1447
            // ```c++
            //         CRef<CSeq_data> data(new CSeq_data(m_SeqData, format));
            //         if ( !TestFlag(fLeaveAsText) ) {
            //             CSeqportUtil::Pack(data, inst.GetLength());
            //         }
            // ```
            for residue in sequence.iter_mut() {
                // A select rather than a conditional store (LOSAT's speed: vectorized).
                *residue = if *residue == b'U' { b'T' } else { *residue };
            }
        }
        for &(from, to) in &self.mask {
            for residue in &mut sequence[from as usize..=to as usize] {
                *residue = residue.to_ascii_lowercase();
            }
        }
        Ok(FastaRecord {
            local_id: self.local_id.take().unwrap_or_default(),
            title: std::mem::take(&mut self.title),
            sequence,
            warnings: Vec::new(),
        })
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:670-686
    /// ```c++
    /// void CFastaReader::ParseTitle(
    ///     const SLineTextAndLoc & lineInfo, ILineErrorListener * pMessageListener)
    /// {
    ///     const static size_t kWarnTitleLength = 1000;
    ///     if( lineInfo.m_sLineText.length() > kWarnTitleLength ) {
    ///         FASTA_WARNING(lineInfo.m_iLineNum,
    ///             "FASTA-Reader: Title is very long: " << lineInfo.m_sLineText.length()
    ///             << " characters (max is " << kWarnTitleLength << ")",
    ///             ILineError::eProblem_TooLong, "defline");
    ///     }
    ///
    ///     CreateWarningsForSeqDataInTitle(lineInfo.m_sLineText,lineInfo.m_iLineNum, pMessageListener);
    ///
    ///     CTempString title(lineInfo.m_sLineText.data(), lineInfo.m_sLineText.length());
    ///     x_ApplyMods(title, lineInfo.m_iLineNum, *m_CurrentSeq, pMessageListener);
    /// }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2018-2043
    /// ```c++
    /// void CFastaReader::x_ApplyMods(
    ///      const string& title,
    ///      TSeqPos line_number,
    ///      CBioseq& bioseq,
    ///      ILineErrorListener* pMessageListener )
    /// {
    ///     string processed_title = title;
    ///     if (TestFlag(fAddMods)) {
    ///         x_AddMods(line_number, bioseq, processed_title, pMessageListener);
    ///     }
    ///     else
    ///     if (!TestFlag(fIgnoreMods) &&
    ///         CTitleParser::HasMods(title)) {
    ///         FASTA_WARNING(line_number,
    ///         "FASTA-Reader: Ignoring FASTA modifier(s) found because "
    ///         "the input was not expected to have any.",
    ///         ILineError::eProblem_ModifierFoundButNoneExpected,
    ///         "defline");
    ///     }
    ///
    ///     NStr::TruncateSpacesInPlace(processed_title);
    ///     if (!processed_title.empty()) {
    ///         auto pDesc = Ref(new CSeqdesc());
    ///         pDesc->SetTitle() = processed_title;
    ///         bioseq.SetDescr().Set().push_back(std::move(pDesc));
    ///     }
    /// }
    /// ```
    /// The long-title and modifier warnings are ignored problems.
    fn parse_title(&mut self, title: LineText, warn: Warn<'_>) -> Result<(), ReadError> {
        self.create_warnings_for_seq_data_in_title(&title.text, warn)?;
        // `TruncateSpacesInPlace` on the title's own bytes.
        let mut text = title.text;
        text.truncate(trim_input_end(&text).len());
        let leading = text.len() - trim_input_start(&text).len();
        text.drain(..leading);
        self.title = text;
        Ok(())
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1623-1680
    /// ```c++
    ///     // check for nuc or aa sequences at the end of the title
    ///     const static size_t kWarnNumNucCharsAtEnd = 20;
    ///     const static size_t kWarnAminoAcidCharsAtEnd = 50;
    ///
    ///     const size_t length = sLineText.length();
    ///     SIZE_TYPE pos_to_check = length-1;
    ///
    ///     if((length > kWarnNumNucCharsAtEnd) && !TestFlag(fAssumeProt)) {
    ///         // find last non-nuc character, within the last kWarnNumNucCharsAtEnd characters
    ///         const SIZE_TYPE last_pos_to_check_for_nuc = (sLineText.length() - kWarnNumNucCharsAtEnd);
    ///         for( ; pos_to_check >= last_pos_to_check_for_nuc; --pos_to_check ) {
    ///             if( ! s_ASCII_IsUnAmbigNuc(sLineText[pos_to_check]) ) {
    ///                 // found a character which is not an unambiguous nucleotide
    ///                 break;
    ///             }
    ///         }
    ///         if( pos_to_check < last_pos_to_check_for_nuc ) {
    ///             FASTA_WARNING(iLineNum,
    ///                 "FASTA-Reader: Title ends with at least " << kWarnNumNucCharsAtEnd
    ///                 << " valid nucleotide characters.  Was the sequence "
    ///                 << "accidentally put in the title line?",
    ///                 ILineError::eProblem_UnexpectedNucResidues,
    ///                 "defline"
    ///                 );
    ///             return true; // found problem
    ///         }
    ///     }
    ///
    ///     if((length > kWarnAminoAcidCharsAtEnd) && !TestFlag(fAssumeNuc)) {
    ///         // check for aa's at the end of the title
    ///         // for efficiency, continue where the nuc search left off, since
    ///         // we know that nucs can be amino acids, also
    ///         const SIZE_TYPE last_pos_to_check_for_amino_acid =
    ///             ( sLineText.length() - kWarnAminoAcidCharsAtEnd );
    ///         for( ; pos_to_check >= last_pos_to_check_for_amino_acid; --pos_to_check ) {
    ///             // can't just use "isalpha" in case it includes characters
    ///             // with diacritics (an accent, tilde, umlaut, etc.)
    ///             const char ch = sLineText[pos_to_check];
    ///             if( ( ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') ) {
    ///                 // potential amino acid, so keep going
    ///             } else {
    ///                 // non-amino-acid found
    ///                 break;
    ///             }
    ///         }
    ///
    ///         if( pos_to_check < last_pos_to_check_for_amino_acid ) {
    ///             FASTA_WARNING(iLineNum,
    ///                 "FASTA-Reader: Title ends with at least " << kWarnAminoAcidCharsAtEnd
    ///                 << " valid amino acid characters.  Was the sequence "
    ///                 << "accidentally put in the title line?",
    ///                 ILineError::eProblem_UnexpectedAminoAcids,
    ///                 "defline");
    ///             return true; // found problem
    ///         }
    ///     }
    ///
    ///     return false;
    /// ```
    /// Exactly one of `fAssumeNuc` and `fAssumeProt` is set, so exactly one test runs, from
    /// the end of the title.
    fn create_warnings_for_seq_data_in_title(
        &self,
        text: &[u8],
        warn: Warn<'_>,
    ) -> Result<(), ReadError> {
        if let Some(message) = seq_data_in_title_warning(text, self.config.protein) {
            let problem = if self.config.protein {
                Problem::UnexpectedAminoAcids
            } else {
                Problem::UnexpectedNucResidues
            };
            self.post_warning(problem, message, warn)?;
        }
        Ok(())
    }
}

/// The message of `CreateWarningsForSeqDataInTitle` for the title text `text` (see
/// `FastaReader::create_warnings_for_seq_data_in_title`, which posts it), or `None`;
/// `protein` is `fAssumeProt` (else `fAssumeNuc`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1624-1643
/// ```c++
///     const static size_t kWarnNumNucCharsAtEnd = 20;
///     const static size_t kWarnAminoAcidCharsAtEnd = 50;
///
///     const size_t length = sLineText.length();
///     SIZE_TYPE pos_to_check = length-1;
///
///     if((length > kWarnNumNucCharsAtEnd) && !TestFlag(fAssumeProt)) {
///         // find last non-nuc character, within the last kWarnNumNucCharsAtEnd characters
///         const SIZE_TYPE last_pos_to_check_for_nuc = (sLineText.length() - kWarnNumNucCharsAtEnd);
///         for( ; pos_to_check >= last_pos_to_check_for_nuc; --pos_to_check ) {
///             if( ! s_ASCII_IsUnAmbigNuc(sLineText[pos_to_check]) ) {
///                 // found a character which is not an unambiguous nucleotide
///                 break;
///             }
///         }
///         if( pos_to_check < last_pos_to_check_for_nuc ) {
///             FASTA_WARNING(iLineNum,
///                 "FASTA-Reader: Title ends with at least " << kWarnNumNucCharsAtEnd
///                 << " valid nucleotide characters.  Was the sequence "
///                 << "accidentally put in the title line?",
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1651-1673
/// ```c++
///     if((length > kWarnAminoAcidCharsAtEnd) && !TestFlag(fAssumeNuc)) {
///         ...
///         for( ; pos_to_check >= last_pos_to_check_for_amino_acid; --pos_to_check ) {
///             // can't just use "isalpha" in case it includes characters
///             // with diacritics (an accent, tilde, umlaut, etc.)
///             const char ch = sLineText[pos_to_check];
///             if( ( ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') ) {
///                 // potential amino acid, so keep going
///             } else {
///                 // non-amino-acid found
///                 break;
///             }
///         }
///
///         if( pos_to_check < last_pos_to_check_for_amino_acid ) {
///             FASTA_WARNING(iLineNum,
///                 "FASTA-Reader: Title ends with at least " << kWarnAminoAcidCharsAtEnd
///                 << " valid amino acid characters.  Was the sequence "
///                 << "accidentally put in the title line?",
/// ```
/// Exactly one of `fAssumeNuc` and `fAssumeProt` is set, so exactly one test runs.
pub(crate) fn seq_data_in_title_warning(text: &[u8], protein: bool) -> Option<&'static [u8]> {
    let length = text.len();
    if !protein {
        if length > 20
            && text[length - 20..]
                .iter()
                .all(|byte| matches!(byte, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
        {
            return Some(b"FASTA-Reader: Title ends with at least 20 valid nucleotide characters.  Was the sequence accidentally put in the title line?");
        }
    } else if length > 50 && text[length - 50..].iter().all(u8::is_ascii_alphabetic) {
        return Some(b"FASTA-Reader: Title ends with at least 50 valid amino acid characters.  Was the sequence accidentally put in the title line?");
    }
    None
}

/// The positions of the bad residues of one line, as ranges of 1-based positions.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_exception.cpp:81-150
/// ```c++
/// void CBadResiduesException::SBadResiduePositions::ConvertBadIndexesToString(
///         CNcbiOstream & out,
///         unsigned int maxRanges ) const
/// {
///     const char *line_prefix = "";
///     unsigned int iRangesFound = 0;
///     ITERATE( SBadResiduePositions::TBadIndexMap, index_map_iter, m_BadIndexMap ) {
///         const int lineNum = index_map_iter->first;
///         const vector<TSeqPos> & badIndexesOnLine = index_map_iter->second;
///         ...
///         TRangeVec rangesFound;
///
///         ITERATE( vector<TSeqPos>, idx_iter, badIndexesOnLine ) {
///             const TSeqPos idx = *idx_iter;
///
///             // first one
///             if( rangesFound.empty() ) {
///                 rangesFound.push_back(TRange(idx, idx));
///                 ++iRangesFound;
///                 continue;
///             }
///
///             const TSeqPos last_idx = rangesFound.back().second;
///             if( idx == (last_idx+1) ) {
///                 // extend previous range
///                 ++rangesFound.back().second;
///                 continue;
///             }
///
///             if( iRangesFound >= maxRanges ) {
///                 break;
///             }
///
///             // create new range
///             rangesFound.push_back(TRange(idx, idx));
///             ++iRangesFound;
///         }
///
///         // turn the ranges found on this line into a string
///         out << line_prefix << "On line " << lineNum << ": ";
///         line_prefix = ", ";
///
///         const char *pos_prefix = "";
///         for( unsigned int rng_idx = 0;
///             ( rng_idx < rangesFound.size() );
///             ++rng_idx )
///         {
///             out << pos_prefix;
///             const TRange &range = rangesFound[rng_idx];
///             out << (range.first + 1); // "+1" because 1-based for user
///             if( range.first != range.second ) {
///                 out << "-" << (range.second + 1); // "+1" because 1-based for user
///             }
///
///             pos_prefix = ", ";
///         }
///         if (iRangesFound > maxRanges) {
///             out << ", and more";
///             return;
///         }
///     }
/// }
/// ```
/// `maxRanges` is the default 1000 (fasta_exception.hpp:90-92); `iRangesFound` never
/// exceeds it, so ", and more" is not written.
fn convert_bad_indexes_to_string(out: &mut Vec<u8>, line: u64, positions: &[u32]) {
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
        } else if byte == b'-' || input_space(byte) || byte.is_ascii_digit() {
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
