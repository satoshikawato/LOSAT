//! NCBI's `CFastaReader` with the flags of `CFastaBioseqSource` without `-parse_seqids`:
//! `fForceType`, `fAssumeProt` (protein) or `fAssumeNuc | fParseGaps` (nucleotide),
//! `fNoParseID`, `fQuickIDCheck`, `fDisableNoResidues`.
//!
//! What differs from the reader of the BLAST programs (`blastinput::fasta_reader`, with
//! `fDLOptional | fNoSplit | fHyphensIgnoreAndWarn`): a definition line is required; a
//! hyphen is a residue in a protein and a gap in a nucleotide sequence; the residues are
//! split and packed into a delta sequence (`pack`); the definition line is kept in a user
//! object for the database's title (`title`); no lowercase mask is recorded
//! (`m_CurrentMask` is empty: `SaveMask` is never called).

use std::io::Read;

use super::line_reader::BufferedLineReader;
use super::pack::{self, Coding, Segment};
use super::title::printable_string;
use super::Warn;
use crate::blastinput::fasta_reader::reader::{
    canonicalize_string, convert_bad_indexes_to_string, gap_type_info, linkage_evidence_value,
    seq_data_in_title_warning, GapType, LinkEvid, Problem, LINKAGE_UNSPECIFIED,
};
use crate::blastinput::fasta_reader::stream::{
    input_space, trim_input_end, trim_input_space, trim_input_start,
};
use crate::writedb::bioseq::SeqData;

/// The problems that `CFastaBioseqSource` tells the reader to ignore.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:735-737
/// ```c++
///     m_FastaReader->IgnoreProblem(ILineError::eProblem_ModifierFoundButNoneExpected);
///     m_FastaReader->IgnoreProblem(ILineError::eProblem_TooManyAmbiguousResidues);
///     m_FastaReader->IgnoreProblem(ILineError::eProblem_TooLong);
/// ```
const IGNORED_PROBLEMS: [Problem; 3] = [
    Problem::ModifierFoundButNoneExpected,
    Problem::TooManyAmbiguousResidues,
    Problem::TooLong,
];

/// A `CObjReaderParseException` of the reader: its message (`GetMsg()`).
#[derive(Debug)]
pub(super) struct ParseError {
    pub(super) message: String,
}

/// One gap of a record (`CFastaReader::SGap`): the raw residues before it and its
/// length; length 0 is a gap of unknown length (a lone hyphen at the end of a line).
/// The gap type and the linkage evidence of a `>?` line only decide whether the delta
/// literal of the gap gets a `Seq-data.gap`, which the database is built from like a
/// literal without data (`SeqInst`), so they are not kept; nor is the size's fuzz.
#[derive(Clone, Copy)]
struct Gap {
    position: usize,
    length: u64,
}

/// One piece of a delta sequence (`CDelta_seq`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum DeltaSeq {
    /// `Seq-loc null`: a gap of unknown length.
    LocNull,
    /// A literal without residues (no `Seq-data`, or `Seq-data.gap`): a gap.
    Gap { length: u64 },
    /// A literal with NCBI2na or NCBI4na data: residues `start..start + length`.
    Literal(Segment),
}

/// `Seq-inst` of a record as `CFastaReader::AssembleSeq` leaves it.
#[derive(PartialEq, Eq, Debug)]
pub(super) enum SeqInst {
    /// `repr virtual`: no residues and no gaps; length 0, no `seq-data`.
    Virtual,
    /// `repr raw`: a delta of one piece, whose `Seq-data` became the sequence's.
    Raw { length: u64, seq_data: SeqData },
    /// `repr delta` of more than one piece over the IUPACna `residues`.
    Delta {
        length: u64,
        residues: Vec<u8>,
        pieces: Vec<DeltaSeq>,
    },
    /// A record longer than 2,147,483,647 residues, which `CBuildDatabase::AddSequences`
    /// rejects by its length before it reads the sequence data: LOSAT does not build it.
    TooLong { length: u64 },
}

/// A record as `CFastaReader::ReadOneSeq` returns it.
#[derive(Debug)]
pub(super) struct Record {
    /// `Seq-id local id`, the number that `CSeqIdGenerator` gave the record.
    pub(super) local_id: u64,
    /// The `DefLine` field of the `CFastaReader` user object.
    pub(super) defline_field: Vec<u8>,
    /// The title descriptor (`x_ApplyMods`), if the record has one.
    pub(super) title_descr: Option<Vec<u8>>,
    pub(super) inst: SeqInst,
}

/// The title of a definition line (`SLineTextAndLoc`; its line number is used for the
/// reader's listener only).
struct LineText {
    text: Vec<u8>,
}

/// What `ParseDataLine` does with a byte: a value from `RESIDUE_MIN` up is the residue
/// to store (the letter in upper case, `*`, or a protein's `-`).
const BAD: u8 = 0;
const SPACE: u8 = 1;
const COMMENT: u8 = 2;
const GAP: u8 = 3;
const RESIDUE_MIN: u8 = b'*';

/// The character switch of `ParseDataLine` for protein (`fAssumeProt`) or nucleotide
/// (`fAssumeNuc | fParseGaps`) input.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:807-814,822-826,856-936
/// ```c++
///     // we're stricter with nucs, so try to determine if we should
///     // assume this is a nuc
///     const bool bIsNuc = (
///         (! TestFlag(fForceType) &&
///             FIELD_CHAIN_OF_2_IS_SET(*m_CurrentSeq, Inst, Mol))
///         ? m_CurrentSeq->IsNa()
///         : TestFlag(fAssumeNuc)
///     );
/// ...
///     const bool bHyphensIgnoreAndWarn = TestFlag(fHyphensIgnoreAndWarn);
///     const bool bHyphensAreGaps =
///         ( TestFlag(fParseGaps) && ! bHyphensIgnoreAndWarn );
///     const bool bAllowLetterGaps =
///         ( TestFlag(fParseGaps) && TestFlag(fLetterGaps) );
/// ...
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
/// `bIsNuc` is `fAssumeNuc` (`fForceType`). Neither input has `fLetterGaps` or
/// `fHyphensIgnoreAndWarn`: `bAllowLetterGaps` is false, and `bHyphensAreGaps` is
/// `fParseGaps`, so a hyphen is a gap in nucleotide input and a residue in protein input.
/// A lower-case letter is stored in upper case (`eCharType_MaskedNonGap`; `OpenMask` does
/// nothing without a mask). LOSAT looks the action up in a table (speed only).
const fn char_actions(protein: bool) -> [u8; 256] {
    let mut table = [BAD; 256];
    let letters: &[u8] = if protein {
        b"ABCDGHKMRSTUVWYEFIJLOPQZNX"
    } else {
        b"ABCDGHKMRSTUVWYN"
    };
    let mut index = 0;
    while index < letters.len() {
        table[letters[index] as usize] = letters[index];
        table[letters[index].to_ascii_lowercase() as usize] = letters[index];
        index += 1;
    }
    if protein {
        table[b'*' as usize] = b'*';
        table[b'-' as usize] = b'-';
    } else {
        table[b'-' as usize] = GAP;
    }
    table[b';' as usize] = COMMENT;
    let spaces = b"\t\n\x0b\x0c\r ";
    let mut index = 0;
    while index < spaces.len() {
        table[spaces[index] as usize] = SPACE;
        index += 1;
    }
    table
}

const NUCLEOTIDE_ACTIONS: [u8; 256] = char_actions(false);
const PROTEIN_ACTIONS: [u8; 256] = char_actions(true);

/// The state of the record being read (`ReadOneSeq` resets it).
struct RecordState {
    /// `fAssumeProt` (else `fAssumeNuc | fParseGaps`).
    protein: bool,
    /// `m_SeqData` up to `m_CurrentPos`: the raw residues, in upper case.
    seq_data: Vec<u8>,
    /// `m_Gaps`.
    gaps: Vec<Gap>,
    /// `m_CurrentGapLength`: the hyphens read since the last residue (their character,
    /// `m_CurrentGapChar`, is always `-`).
    current_gap_length: u64,
    /// `m_TotalGapLength`.
    total_gap_length: u64,
    /// `m_CurrentSeqTitles`.
    titles: Vec<LineText>,
    /// The `DefLine` field of the `CFastaReader` user object.
    defline_field: Vec<u8>,
    /// The title descriptor.
    title_descr: Option<Vec<u8>>,
    /// The local ID of the record.
    local_id: u64,
    /// The bad positions of a data line, reused.
    bad_positions: Vec<u32>,
}

/// NCBI's `CFastaReader` over a `CBufferedLineReader`.
pub(super) struct FastaReader<R: Read> {
    pub(super) lines: BufferedLineReader<R>,
    /// `CSeqIdGenerator::m_Counter`.
    id_counter: u64,
    state: RecordState,
}

/// `PostWarning` for a warning: nothing for an ignored problem, else the message.
///
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
/// `CFastaBioseqSource::GetNext` gives the reader no listener, so a warning is the
/// message on a line of its own (`LOG_POST` writes no severity).
fn post_warning(problem: Problem, message: &[u8], warn: Warn<'_>) {
    if !IGNORED_PROBLEMS.contains(&problem) {
        warn(message);
    }
}

impl<R: Read> FastaReader<R> {
    /// `CFastaReader(*m_LineReader, flags)`; the IDs come from the reader's own
    /// `CSeqIdGenerator`, which counts from 1 and has no prefix.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:701-737
    /// ```c++
    /// CFastaBioseqSource::CFastaBioseqSource(CNcbiIstream & fasta_file,
    ///                                        bool is_protein,
    ///                                        bool parse_ids,
    ///                                        bool long_ids)
    ///     : m_FastaReader(NULL)
    /// {
    ///     m_LineReader.Reset(new CBufferedLineReader(fasta_file));
    ///     typedef CFastaReader::EFlags TFlags;
    ///
    ///     int iflags = CFastaReader::fForceType;
    ///
    ///     if (is_protein) {
    ///         iflags |= CFastaReader::fAssumeProt;
    ///     } else {
    ///         iflags |= CFastaReader::fAssumeNuc;
    ///         iflags |=  CFastaReader::fParseGaps;
    ///     }
    ///
    ///     if (parse_ids) {
    ///         iflags |= CFastaReader::fRequireID;
    ///         // parse bare accessions
    ///         if (!long_ids) {
    ///             iflags |= CFastaReader::fParseRawID;
    ///         }
    ///     } else {
    ///         iflags |= CFastaReader::fNoParseID;
    ///     }
    ///
    ///     iflags |= CFastaReader::fQuickIDCheck;
    ///     iflags |= CFastaReader::fDisableNoResidues;
    ///
    ///     TFlags flags = (TFlags) iflags;
    ///
    ///     m_FastaReader = new CFastaReader(*m_LineReader, flags);
    ///     m_FastaReader->IgnoreProblem(ILineError::eProblem_ModifierFoundButNoneExpected);
    ///     m_FastaReader->IgnoreProblem(ILineError::eProblem_TooManyAmbiguousResidues);
    ///     m_FastaReader->IgnoreProblem(ILineError::eProblem_TooLong);
    /// ```
    /// NCBI reference (598d8ae6): c++/include/objtools/readers/fasta_reader_utils.hpp:160-165
    /// ```c++
    ///     CSeqIdGenerator(TCount count = 1,
    ///                     const string& prefix = kEmptyStr,
    ///                     const string& suffix = kEmptyStr) :
    ///         m_Prefix(prefix),
    ///         m_Suffix(suffix),
    ///         m_Counter(count)
    /// ```
    /// The registry value `BLAST/MAX_SEQID_LENGTH` (build_db.cpp:738-747) is not read:
    /// LOSAT has no `.ncbirc`, and without `-parse_seqids` no ID is parsed.
    pub(super) fn new(input: R, protein: bool) -> Self {
        Self {
            lines: BufferedLineReader::new(input),
            id_counter: 1,
            state: RecordState {
                protein,
                seq_data: Vec::new(),
                gaps: Vec::new(),
                current_gap_length: 0,
                total_gap_length: 0,
                titles: Vec::new(),
                defline_field: Vec::new(),
                title_descr: None,
                local_id: 0,
                bad_positions: Vec::new(),
            },
        }
    }

    /// `CBioseq::IsAa()` of every record: the molecule is the input's (`assemble_seq`).
    pub(super) fn is_protein(&self) -> bool {
        self.state.protein
    }

    /// One record.
    ///
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
    ///                 const auto lineNum = LineNumber();
    ///                 GetLineReader().UngetLine();
    ///                 NCBI_THROW2(CObjReaderParseException, eNoDefline,
    ///                             "CFastaReader: Input doesn't start with"
    ///                             " a defline or comment around line " + NStr::NumericToString(lineNum),
    ///                              lineNum);
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
    ///                 // we have to catch this exception so we can build up
    ///                 // information on all lines, not just the first line
    ///                 // with a bad residue
    ///                 bad_residue_positions.m_SeqId = e.GetBadResiduePositions().m_SeqId;
    ///                 bad_residue_positions.AddBadIndexMap(e.GetBadResiduePositions().m_BadIndexMap);
    ///             }
    ///         }
    ///     }
    ///
    ///     if( ! bad_residue_positions.m_BadIndexMap.empty() ) {
    ///         // bad residues unconditionally throws, for now.
    ///         // (not worth the refactoring at this time)
    ///         NCBI_THROW2(CBadResiduesException, eBadResidues,
    ///             "CFastaReader: There are invalid " + x_NucOrProt() + "residue(s) in input sequence",
    ///             bad_residue_positions );
    ///     }
    ///
    ///     if (need_defline  &&  GetLineReader().AtEOF()) {
    ///         FASTA_ERROR(LineNumber(),
    ///             "CFastaReader: Expected defline around line " << LineNumber(),
    ///             CObjReaderParseException::eEOF);
    ///     }
    ///
    ///     AssembleSeq(pMessageListener);
    ///     CRef<CSeq_entry> entry(new CSeq_entry);
    ///     entry->SetSeq(*m_CurrentSeq);
    ///
    ///     entry->Parentize();
    ///     return entry;
    /// }
    /// ```
    /// `PeekChar` of `CBufferedLineReader` does not read, so the reader is not at the end
    /// of the input after it. `FASTA_PROGRESS` posts to a listener only, and there is
    /// none. Without `fDLOptional` data before a definition line is the `eNoDefline`
    /// error. `CBadResiduesException` comes only with `fValidate`, which is not set.
    pub(super) fn read_one_seq(&mut self, warn: Warn<'_>) -> Result<Record, ParseError> {
        let state = &mut self.state;
        state.seq_data.clear();
        state.gaps.clear();
        state.current_gap_length = 0;
        state.total_gap_length = 0;
        state.titles.clear();
        state.defline_field.clear();
        state.title_descr = None;

        let mut need_defline = true;
        while !self.lines.at_eof() {
            let c = self.lines.peek_char();
            if c == b'>' {
                self.lines.next_line();
                let raw = self.lines.line();
                let modified;
                let next_line = match raw.strip_prefix(b">?_") {
                    Some(rest) => {
                        modified = [b">".as_slice(), rest].concat();
                        &modified[..]
                    }
                    None => raw,
                };
                if next_line.starts_with(b">?") {
                    // A gap line: a data line, read again below.
                    self.lines.unget_line();
                } else if need_defline {
                    state.parse_def_line(next_line, &mut self.id_counter);
                    need_defline = false;
                    continue;
                } else {
                    // The start of the next sequence.
                    self.lines.unget_line();
                    break;
                }
            }

            self.lines.next_line();
            let line_number = self.lines.line_number();
            let raw = self.lines.line();
            let line = trim_input_space(raw);
            if line.is_empty() || matches!(line[0], b'!' | b'#' | b';') {
                continue;
            }
            if need_defline {
                self.lines.unget_line();
                return Err(ParseError {
                    message: format!(
                        "CFastaReader: Input doesn't start with a defline or comment around line {line_number}"
                    ),
                });
            }
            let raw_length = raw.len();
            match line.strip_prefix(b">?_") {
                Some(rest) => state.parse_data_line(
                    &[b">".as_slice(), rest].concat(),
                    raw_length,
                    line_number,
                    warn,
                )?,
                None => state.parse_data_line(line, raw_length, line_number, warn)?,
            }
        }

        if need_defline && self.lines.at_eof() {
            return Err(ParseError {
                message: format!(
                    "CFastaReader: Expected defline around line {}",
                    self.lines.line_number()
                ),
            });
        }

        let inst = state.assemble_seq(warn);
        Ok(Record {
            local_id: state.local_id,
            defline_field: std::mem::take(&mut state.defline_field),
            title_descr: state.title_descr.take(),
            inst,
        })
    }
}

impl RecordState {
    /// `ParseDefLine` with `fNoParseID`: the title for the reader's checks
    /// (`CFastaDeflineReader::ParseDefline`), the generated local ID (`PostProcessIDs`
    /// without IDs), and the user object with the definition line.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:157-180,210-226
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
    /// ...
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
    /// }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:566-574,591-603
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
    /// ...
    ///     PostProcessIDs(data.ids, s, data.has_range, data.range_start, data.range_end);
    ///
    ///     m_BestID = FindBestChoice(GetIDs(), CSeq_id::BestRank);
    ///     FASTA_PROGRESS("Processing Seq-id: " <<
    ///         ( m_BestID ? m_BestID->AsFastaString() : string("UNKNOWN") ) );
    ///
    ///     if ( !TestFlag(fNoUserObjs) ) {
    ///         // store the raw defline in a User-object for reference
    ///         CRef<CSeqdesc> desc(new CSeqdesc);
    ///         desc->SetUser().SetType().SetStr("CFastaReader");
    ///         desc->SetUser().AddField("DefLine", NStr::PrintableString(s));
    ///         m_CurrentSeq->SetDescr().Set().push_back(desc);
    ///     }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:626-628
    /// ```c++
    ///     if (defline_ids.empty()) {
    ///         GenerateID();
    ///     }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:466-474
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
    /// ```
    /// `defline` starts with `>` except for an empty line (`BufferedLineReader`'s module
    /// comment), which returns before the test of its first character. The molecule is
    /// not set here (`fForceType`), and no range is parsed (`fNoParseID`).
    fn parse_def_line(&mut self, defline: &[u8], id_counter: &mut u64) {
        self.titles.clear();
        let len = defline.len();
        if len > 1 && !defline[1..].iter().all(|&byte| input_space(byte)) {
            let start = 1 + defline[1..]
                .iter()
                .position(|&byte| !input_space(byte))
                .unwrap_or(len - 1);
            let title_start = start;
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
        self.local_id = *id_counter;
        *id_counter += 1;
        self.defline_field = printable_string(defline).into_owned();
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:710-771
    /// ```c++
    /// void CFastaReader::CheckDataLine(
    ///     const TStr& s, ILineErrorListener * pMessageListener)
    /// {
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
    /// }
    /// ```
    /// `bIgnoreHyphens` is false: a hyphen counts as good. The warning about ambiguous
    /// nucleotides is an ignored problem.
    fn check_data_line(&self, s: &[u8], line_number: u64) -> Result<(), ParseError> {
        if !self.seq_data.is_empty() {
            return Ok(());
        }
        let (mut good, mut bad) = (0usize, 0usize);
        let len_to_check = s.len().min(70);
        for &c in &s[..len_to_check] {
            if c.is_ascii_alphabetic() || c == b'*' || c == b'-' {
                good += 1;
            } else if input_space(c) || c.is_ascii_digit() {
                // White space and digits are neutral.
            } else if c == b';' {
                break;
            } else {
                bad += 1;
            }
        }
        if bad >= good / 3 && (len_to_check > 3 || good == 0 || bad > good) {
            return Err(ParseError {
                message: format!(
                    "CFastaReader: Near line {line_number}, there's a line that doesn't look like plausible data, but it's not marked as defline or comment."
                ),
            });
        }
        Ok(())
    }

    /// `s` is the line without white space at its ends, `raw_length` the length of the
    /// line as read (`(*GetLineReader()).length()`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:774-782,816-816,938-1016
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
    /// ...
    ///     m_SeqData.resize(m_CurrentPos + s_len);
    /// ...
    ///         switch(char_type) {
    ///         case eCharType_NormalNonGap:
    ///             CloseGap(pos == 0);
    ///             m_SeqData[m_CurrentPos] = c;
    ///             CloseMask();
    ///             ++m_CurrentPos;
    ///             break;
    ///         case eCharType_MaskedNonGap:
    ///             CloseGap(pos == 0);
    ///             m_SeqData[m_CurrentPos] = s_ASCII_MustBeLowerToUpper(c);
    ///             OpenMask();
    ///             ++m_CurrentPos;
    ///             break;
    ///         case eCharType_Gap: {
    ///             CloseMask();
    ///             // open a gap
    ///
    ///             size_t pos2 = pos + 1;
    ///             while( pos2 < s_len && s[pos2] == c ) {
    ///                 ++pos2;
    ///             }
    ///             _ASSERT(pos2 <= s_len);
    ///             m_CurrentGapLength += pos2 - pos;
    ///             m_CurrentGapChar = toupper(c);
    ///             pos = pos2 - 1; // `- 1` compensates for the `++pos` in the `for`
    ///             break;
    ///         }
    ///         case eCharType_JustIgnore:
    ///             break;
    ///         case eCharType_HyphenToIgnoreAndWarn:
    ///             bIgnorableHyphenSeen = true;
    ///             break;
    ///         case eCharType_Comment:
    ///             // artificially advance pos to the end to break the pos loop
    ///             pos = s_len;
    ///             break;
    ///         case eCharType_Bad:
    ///             if( bad_pos_line_num < 0 ) {
    ///                 bad_pos_line_num = LineNumber();
    ///             }
    ///             bad_pos_vec.push_back(pos);
    ///             break;
    ///         default:
    ///             _TROUBLE;
    ///         }
    ///     }
    ///
    ///     m_SeqData.resize(m_CurrentPos);
    ///
    ///     if( bIgnorableHyphenSeen ) {
    ///         _ASSERT( bHyphensIgnoreAndWarn );
    ///         FASTA_WARNING_EX(LineNumber(),
    ///             "CFastaReader: Hyphens are invalid and will be ignored around line " << LineNumber(),
    ///             ILineError::eProblem_IgnoredResidue,
    ///             kEmptyStr, kEmptyStr, "-" );
    ///     }
    ///
    ///     // before throwing, be sure that we're in a valid state so that callers can
    ///     // parse multiple lines and get the invalid residues in all of them.
    ///
    ///     if( ! bad_pos_vec.empty() ) {
    ///         if (TestFlag(fValidate)) {
    ///                         NCBI_THROW2(CBadResiduesException, eBadResidues,
    ///                 "CFastaReader: There are invalid " + x_NucOrProt() + "residue(s) in input sequence",
    ///                 CBadResiduesException::SBadResiduePositions( m_BestID, bad_pos_vec, bad_pos_line_num ) );
    ///         } else {
    ///             stringstream warn_strm;
    ///             warn_strm << "FASTA-Reader: Ignoring invalid " << x_NucOrProt()
    ///                 << "residues at position(s): ";
    ///             CBadResiduesException::SBadResiduePositions(
    ///                 m_BestID, bad_pos_vec, bad_pos_line_num ).ConvertBadIndexesToString(warn_strm);
    ///
    ///             FASTA_WARNING(0,
    ///                 warn_strm.str(),
    ///                 ILineError::eProblem_InvalidResidue,
    ///                 kEmptyStr );
    ///         }
    ///     }
    /// }
    /// ```
    /// The character switch is quoted at `char_actions`. `m_SeqData.resize(m_CurrentPos +
    /// s_len)` makes the string's length the residues before the line plus the line's
    /// length while the line is read; `x_CloseGap` reads that length (`close_gap`). The
    /// fast copy of `fSkipCheck` is not taken (no `fSkipCheck`). `x_NucOrProt()` is empty:
    /// the molecule is set only in `AssembleSeq`.
    fn parse_data_line(
        &mut self,
        s: &[u8],
        raw_length: usize,
        line_number: u64,
        warn: Warn<'_>,
    ) -> Result<(), ParseError> {
        if s.starts_with(b">?") {
            self.parse_gap_line(s, line_number, warn);
            return Ok(());
        }
        self.check_data_line(s, line_number)?;

        let actions = if self.protein {
            &PROTEIN_ACTIONS
        } else {
            &NUCLEOTIDE_ACTIONS
        };
        let seq_data_length = self.seq_data.len() + s.len();
        self.bad_positions.clear();
        let mut pos = 0;
        while pos < s.len() {
            let c = s[pos];
            match actions[c as usize] {
                SPACE => {}
                COMMENT => break,
                GAP => {
                    let run = 1 + s[pos + 1..]
                        .iter()
                        .position(|&byte| byte != c)
                        .unwrap_or(s.len() - pos - 1);
                    self.current_gap_length += run as u64;
                    pos += run - 1;
                }
                BAD => self.bad_positions.push(pos as u32),
                _ => {
                    self.close_gap(pos == 0, seq_data_length, raw_length);
                    // The residues that follow are stored at once (LOSAT's speed: no gap
                    // is open after the first one).
                    let run = 1 + s[pos + 1..]
                        .iter()
                        .position(|&byte| actions[byte as usize] < RESIDUE_MIN)
                        .unwrap_or(s.len() - pos - 1);
                    self.seq_data
                        .extend(s[pos..pos + run].iter().map(|&byte| actions[byte as usize]));
                    pos += run - 1;
                }
            }
            pos += 1;
        }

        if !self.bad_positions.is_empty() {
            let mut message = b"FASTA-Reader: Ignoring invalid residues at position(s): ".to_vec();
            convert_bad_indexes_to_string(&mut message, line_number, &self.bad_positions);
            post_warning(Problem::InvalidResidue, &message, warn);
        }
        Ok(())
    }

    /// `CloseGap(atStartOfLine)`: ends the run of hyphens, if one is open.
    /// `seq_data_length` is `m_SeqData.length()` (see `parse_data_line`; outside a data
    /// line it is the number of residues), `line_length` the length of the reader's
    /// current line.
    ///
    /// NCBI reference (598d8ae6): c++/include/objtools/readers/fasta.hpp:284-289
    /// ```c++
    ///     void CloseGap(bool atStartOfLine=false, ILineErrorListener* pMessageListener = nullptr) {
    ///         if (m_CurrentGapLength > 0) {
    ///             x_CloseGap(m_CurrentGapLength, atStartOfLine, pMessageListener);
    ///             m_CurrentGapLength = 0;
    ///         }
    ///     }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1018-1077
    /// ```c++
    /// void CFastaReader::x_CloseGap(
    ///     TSeqPos len, bool atStartOfLine, ILineErrorListener * /*pMessageListener*/)
    /// {
    ///     _ASSERT(len > 0  &&  TestFlag(fParseGaps));
    ///
    ///     if (m_CurrentGapLength < m_gapNmin)
    ///     {
    ///         // the run of N is too short to be assumed as a gap, put N as is.
    ///         m_SeqData.resize(m_SeqData.size() + m_CurrentGapLength, 'X');
    ///         memset(&m_SeqData.at(m_CurrentPos), m_CurrentGapChar, m_CurrentGapLength);
    ///         m_CurrentPos += m_CurrentGapLength;
    ///         return;
    ///     }
    ///
    ///     if (TestFlag(fAligning)) {
    ///         TSeqPos pos = GetCurrentPos(ePosWithGapsAndSegs);
    ///         m_Starts[pos + m_Offset][m_Row] = CFastaAlignmentBuilder::kNoPos;
    ///         m_Offset += len;
    ///         m_Starts[pos + m_Offset][m_Row] = pos;
    ///     } else {
    ///         TSeqPos pos = GetCurrentPos(eRawPos);
    ///         SGap::EKnownSize eKnownSize = SGap::eKnownSize_Yes;
    ///         if (len == m_gap_Unknown_length)
    ///         {
    ///             eKnownSize = SGap::eKnownSize_No;
    ///             //len = 0;
    ///         }
    ///         else
    ///         {
    ///             // Special case -- treat a lone hyphen at the end of a line as
    ///             // a gap of unknown length.
    ///             // (do NOT treat a lone 'N' or 'X' as unknown length)
    ///             if (len == 1 && m_CurrentGapChar == '-' ) {
    ///                 TSeqPos l = m_SeqData.length();
    ///                 if ((l == pos) || (l == pos + (*GetLineReader()).length() && atStartOfLine)) {
    ///                     //and it's not the first col of the line
    ///                     len = 0;
    ///                     eKnownSize = SGap::eKnownSize_No;
    ///                 }
    ///             }
    ///         }
    ///
    ///         const auto cit = m_GapsizeToLinkageEvidence.find(len);
    ///         const auto& gap_linkage_evidence =
    ///             (cit != m_GapsizeToLinkageEvidence.end()) ?
    ///             cit->second :
    ///             m_DefaultLinkageEvidence;
    ///
    ///         TGapRef pGap( new SGap(
    ///             pos, len,
    ///             eKnownSize,
    ///             LineNumber(),
    ///             m_gap_type,
    ///             gap_linkage_evidence));
    ///
    ///         m_Gaps.push_back(pGap);
    ///         m_TotalGapLength += len;
    ///         m_CurrentGapLength = 0;
    ///     }
    /// }
    /// ```
    /// `m_gapNmin` and `m_gap_Unknown_length` are 0 (`SetMinGaps` is not called), the
    /// reader is not aligning, and the gap character is always `-`. So a lone hyphen is a
    /// gap of unknown length (length 0) when it is closed by the end of the record or by
    /// a `>?` line (`l == pos`), or by a residue that is the first byte of a line without
    /// white space at its ends (the line's residues are not stored yet, so `l` is `pos`
    /// plus the line's length); closed by any other residue it is a gap of one.
    fn close_gap(&mut self, at_start_of_line: bool, seq_data_length: usize, line_length: usize) {
        if self.current_gap_length == 0 {
            return;
        }
        let mut len = self.current_gap_length;
        let pos = self.seq_data.len();
        if len == 1
            && (seq_data_length == pos
                || (seq_data_length == pos + line_length && at_start_of_line))
        {
            len = 0;
        }
        self.gaps.push(Gap {
            position: pos,
            length: len,
        });
        self.total_gap_length += len;
        self.current_gap_length = 0;
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1094-1184,1342-1349
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
    /// ...
    ///     TGapRef pGap( new SGap(
    ///         uPos, uGapSize, eIsKnown, LineNumber(), pGapType,
    ///         setOfLinkageEvidence ) );
    ///
    ///     m_Gaps.push_back(pGap);
    ///     m_TotalGapLength += pGap->m_uLen;
    ///     return true;
    /// }
    /// ```
    /// The modifiers (quoted at `gap_mod_warnings`) only warn here (`Gap`). The caller's
    /// default gap type (`m_gap_type`) is not set by `CFastaBioseqSource`.
    fn parse_gap_line(&mut self, line: &[u8], line_number: u64, warn: Warn<'_>) {
        self.close_gap(false, self.seq_data.len(), 0);

        let mut remaining = trim_input_space(&line[2..]);
        let position = self.seq_data.len();
        if let Some(rest) = remaining.strip_prefix(b"unk") {
            remaining = trim_input_start(rest);
        }
        let digits = remaining
            .iter()
            .position(|byte| !byte.is_ascii_digit())
            .unwrap_or(remaining.len());
        // `NStr::StringToUInt` with `fConvErr_NoThrow`: 0 for an empty string and for a
        // value above 4,294,967,295.
        let mut gap_size = std::str::from_utf8(&remaining[..digits])
            .ok()
            .and_then(|text| text.parse::<u32>().ok())
            .unwrap_or(0);
        if gap_size == 0 {
            post_warning(
                Problem::NonPositiveLength,
                format!("CFastaReader: Bad gap size at line {line_number}").as_bytes(),
                warn,
            );
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
                post_warning(
                    Problem::ParsingModifiers,
                    format!("CFastaReader: Problem parsing gap mods at line {line_number}")
                        .as_bytes(),
                    warn,
                );
                break;
            };
            let key = trim_input_space(&remaining[1..equal]).to_vec();
            let value = trim_input_space(&remaining[equal + 1..close]).to_vec();
            mods.push((key, value));
            remaining = trim_input_start(&remaining[close + 1..]);
        }
        // `multimap` iterates in key order (the order of insertion for equal keys).
        mods.sort_by(|left, right| left.0.cmp(&right.0));
        gap_mod_warnings(&mods, line_number, warn);

        self.gaps.push(Gap {
            position,
            length: u64::from(gap_size),
        });
        self.total_gap_length += u64::from(gap_size);
    }

    /// `AssembleSeq`: the sequence of the record.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1351-1373,1386-1391,1411-1461
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
    ///
    ///     CSeq_data::E_Choice format
    ///         = ( inst.IsAa() ?
    ///             ( TestFlag(fUseIupacaa) ?
    ///               CSeq_data::e_Iupacaa :
    ///               CSeq_data::e_Ncbieaa ) :
    ///             CSeq_data::e_Iupacna);
    /// ...
    ///     if ( !TestFlag(fParseGaps)  &&  m_TotalGapLength > 0 ) {
    ///         // Encountered >? lines; substitute runs of Ns or Xs as appropriate.
    ///         string    new_data;
    ///         char      gap_char(inst.IsAa() ? 'X' : 'N');
    ///         SIZE_TYPE pos = 0;
    ///         new_data.reserve(GetCurrentPos(ePosWithGaps));
    /// ...
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
    ///     } else {
    ///         CDelta_ext& delta_ext = inst.SetExt().SetDelta();
    ///         inst.SetRepr(CSeq_inst::eRepr_delta);
    ///         inst.SetLength(GetCurrentPos(ePosWithGaps));
    ///         SIZE_TYPE n = m_Gaps.size();
    ///
    ///         if (n==0 || m_Gaps[0]->m_uPos > 0)
    ///         {
    ///             TStr chunk(m_SeqData, 0, (n>0 && m_Gaps[0]->m_uPos > 0)?m_Gaps[0]->m_uPos : inst.GetLength());
    ///             delta_ext.AddAndSplit(chunk, format, chunk.length(), false, !TestFlag(fLeaveAsText));
    ///         }
    ///
    ///         for (SIZE_TYPE i = 0;  i < n;  ++i) {
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1463-1471,1497-1516
    /// ```c++
    ///             // add delta-seq
    ///             CRef<CDelta_seq> gap_ds(new CDelta_seq);
    ///             if (m_Gaps[i]->m_uLen == 0) { // totally unknown
    ///                 gap_ds->SetLoc().SetNull();
    ///             } else { // has a nominal length (normally 100)
    ///                 gap_ds->SetLiteral().SetLength(m_Gaps[i]->m_uLen);
    ///                 if( m_Gaps[i]->m_eKnownSize == SGap::eKnownSize_No ) {
    ///                     gap_ds->SetLiteral().SetFuzz().SetLim(CInt_fuzz::eLim_unk);
    ///                 }
    /// ...
    ///             }
    ///             delta_ext.Set().push_back(std::move(gap_ds));
    ///
    ///             TSeqPos next_start = (i == n-1) ? m_CurrentPos : m_Gaps[i+1]->m_uPos;
    ///             if (next_start != m_Gaps[i]->m_uPos) {
    ///                 TSeqPos seq_len = next_start - m_Gaps[i]->m_uPos;
    ///                 TStr chunk(m_SeqData, m_Gaps[i]->m_uPos, seq_len);
    ///                 delta_ext.AddAndSplit(chunk, format, chunk.length(), false, !TestFlag(fLeaveAsText));
    ///             }
    ///         }
    ///         if (delta_ext.Get().size() == 1) {
    ///             // simplify -- just one piece
    ///             inst.SetRepr(CSeq_inst::eRepr_raw);
    ///             inst.SetSeq_data(delta_ext.Set().front()
    ///                                 ->SetLiteral().SetSeq_data());
    ///             inst.ResetExt();
    ///         }
    ///
    ///     }
    /// }
    /// ```
    /// The molecule is the input's (`AssignMolType` with `fForceType`):
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1569-1573
    /// ```c++
    ///     if (TestFlag(fForceType)) {
    ///         _ASSERT(default_mol != CSeq_inst::eMol_not_set);
    ///         inst.SetMol(default_mol);
    ///         return;
    ///     }
    /// ```
    /// Protein input (no `fParseGaps`) turns its `>?` gaps into runs of `X`, so it has no
    /// gaps left; its residues are one NCBIeaa literal, because `CSeqConvert::Pack` copies
    /// NCBIeaa when the target takes no gaps, and that one piece becomes the sequence:
    ///
    /// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1087-1092,1108-1112
    /// ```c++
    ///     case CSeqUtil::e_Iupacaa:
    ///     case CSeqUtil::e_Ncbieaa:
    ///         if (gaps_ok) { // no point otherwise
    ///             best_coding = best_codings.ncbieaa;
    ///         }
    ///         break;
    /// ...
    ///     } else {
    ///         memcpy(dst.NewSegment(src_coding, length), src,
    ///                GetBytesNeeded(src_coding, length));
    ///         return length;
    ///     }
    /// ```
    /// Nucleotide input keeps its gaps as delta pieces (a gap of length 0 as `Seq-loc
    /// null`), and each run of residues between them is packed into NCBI2na and NCBI4na
    /// literals (`pack::pack_iupacna`). A delta of one piece becomes a raw sequence with
    /// that piece's `Seq-data`: `SetLiteral().SetSeq_data()` of a gap piece creates a
    /// `Seq-data` without a choice when it has none, so a record that is one gap has a
    /// `Seq-data` that the database writer cannot read (`SeqData::Other`).
    fn assemble_seq(&mut self, warn: Warn<'_>) -> SeqInst {
        self.close_gap(false, self.seq_data.len(), 0);

        let titles = std::mem::take(&mut self.titles);
        for title in &titles {
            self.parse_title(&title.text, warn);
        }

        let length = self.seq_data.len() as u64 + self.total_gap_length;
        if length > i32::MAX as u64 {
            return SeqInst::TooLong { length };
        }

        if self.protein && self.total_gap_length > 0 {
            let mut new_data = Vec::with_capacity(length as usize);
            let mut pos = 0;
            for gap in &self.gaps {
                if gap.position > pos {
                    new_data.extend_from_slice(&self.seq_data[pos..gap.position]);
                    pos = gap.position;
                }
                new_data.resize(new_data.len() + gap.length as usize, b'X');
            }
            if self.seq_data.len() > pos {
                new_data.extend_from_slice(&self.seq_data[pos..]);
            }
            self.seq_data = new_data;
            self.gaps.clear();
            self.total_gap_length = 0;
        }

        if self.gaps.is_empty() && self.seq_data.is_empty() {
            return SeqInst::Virtual;
        }
        if self.protein {
            return SeqInst::Raw {
                length,
                seq_data: SeqData::Ncbieaa(std::mem::take(&mut self.seq_data)),
            };
        }

        let mut pieces = Vec::new();
        let n = self.gaps.len();
        if n == 0 || self.gaps[0].position > 0 {
            let end = if n > 0 {
                self.gaps[0].position
            } else {
                self.seq_data.len()
            };
            add_and_split(&mut pieces, &self.seq_data, 0, end);
        }
        for i in 0..n {
            let gap = self.gaps[i];
            pieces.push(if gap.length == 0 {
                DeltaSeq::LocNull
            } else {
                DeltaSeq::Gap { length: gap.length }
            });
            let next_start = if i == n - 1 {
                self.seq_data.len()
            } else {
                self.gaps[i + 1].position
            };
            if next_start != gap.position {
                add_and_split(&mut pieces, &self.seq_data, gap.position, next_start);
            }
        }

        if let [piece] = pieces[..] {
            let seq_data = match piece {
                DeltaSeq::Literal(segment) => match segment.coding {
                    Coding::Ncbi2na => SeqData::Ncbi2na(pack::iupacna_to_ncbi2na(&self.seq_data)),
                    Coding::Ncbi4na => SeqData::Ncbi4na(pack::iupacna_to_ncbi4na(&self.seq_data)),
                },
                DeltaSeq::LocNull | DeltaSeq::Gap { .. } => SeqData::Other,
            };
            return SeqInst::Raw { length, seq_data };
        }
        SeqInst::Delta {
            length,
            residues: std::mem::take(&mut self.seq_data),
            pieces,
        }
    }

    /// `ParseTitle`: the reader's warning about residues at the end of a title, and the
    /// title descriptor.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:670-685
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
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2018-2044
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
    /// The long-title and modifier warnings are ignored problems; the test for residues
    /// at the end of the title is `fasta_reader`'s (`seq_data_in_title_warning`, one of
    /// its two tests runs because exactly one of `fAssumeNuc` and `fAssumeProt` is set).
    fn parse_title(&mut self, text: &[u8], warn: Warn<'_>) {
        if let Some(message) = seq_data_in_title_warning(text, self.protein) {
            let problem = if self.protein {
                Problem::UnexpectedAminoAcids
            } else {
                Problem::UnexpectedNucResidues
            };
            post_warning(problem, message, warn);
        }
        let processed_title = trim_input_end(trim_input_start(text));
        if !processed_title.is_empty() {
            self.title_descr = Some(processed_title.to_vec());
        }
    }
}

/// `CDelta_ext::AddAndSplit` of the IUPACna residues `residues[start..end]`: one literal
/// for each segment of `CSeqConvert::Pack`.
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/Delta_ext.cpp:141-147,163-164,183-185,201-205
/// ```c++
/// char* CDelta_ext_PackTarget::NewSegment(CSeqUtil::TCoding coding,
///                                         TSeqPos length)
/// {
///     CRef<CDelta_seq> ds(new CDelta_seq);
///     CSeq_literal&    lit = ds->SetLiteral();
///     lit.SetLength(length);
///     m_Obj.Set().push_back(ds);
/// ...
///     CODING_CASE_EX(e_Ncbi2na, CNCBI2na, SetNcbi2na, (length + 3) / 4)
///     CODING_CASE_EX(e_Ncbi4na, CNCBI4na, SetNcbi4na, (length + 1) / 2)
/// ...
/// void CDelta_ext::AddAndSplit(const CTempString& src, CSeq_data::E_Choice in_format,
///                              TSeqPos length /* in residues */,
///                              bool gaps_ok, bool allow_packing)
/// ...
///     if (allow_packing)
///     {
///         CDelta_ext_PackTarget dst(*this, gaps_ok);
///         CSeqConvert::Pack(src.data(), length, coding, dst);
///     }
/// ```
fn add_and_split(pieces: &mut Vec<DeltaSeq>, residues: &[u8], start: usize, end: usize) {
    pieces.extend(
        pack::pack_iupacna(&residues[start..end])
            .into_iter()
            .map(|segment| {
                DeltaSeq::Literal(Segment {
                    start: segment.start + start as u32,
                    ..segment
                })
            }),
    );
}

/// The warnings of `ParseGapLine` about the modifiers of a `>?` line.
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
///         pGapType = m_gap_type;
///         if (!m_GapsizeToLinkageEvidence.empty()) {
///             setOfLinkageEvidence = m_GapsizeToLinkageEvidence.begin()->second;
///         }
///         for (const auto& rec : CSeq_gap::GetNameToGapTypeInfoMap())
///         {
///             if (rec.second.m_eType == *pGapType)
///             {
///                 eLinkEvid = rec.second.m_eLinkEvid;
///             }
///         }
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
/// `NStr::Split` with `fSplit_MergeDelimiters | fSplit_Truncate` drops empty pieces. The
/// names are looked up by `fasta_reader`'s tables (`gap_type_info`,
/// `linkage_evidence_value`).
fn gap_mod_warnings(mods: &[(Vec<u8>, Vec<u8>)], line_number: u64, warn: Warn<'_>) {
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
                    post_warning(Problem::ParsingModifiers, &message, warn);
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
                        post_warning(Problem::ParsingModifiers, &message, warn);
                    }
                }
            }
            continue;
        }
        let mut message = b"Unknown gap modifier name(s): ".to_vec();
        message.extend_from_slice(key);
        post_warning(Problem::UnrecognizedQualifierName, &message, warn);
    }
    if conflicting {
        post_warning(
            Problem::ContradictoryModifiers,
            format!("There were conflicting gap-types around line {line_number}").as_bytes(),
            warn,
        );
    }
    match link_evid {
        LinkEvid::UnspecifiedOnly => {
            if !evidence.is_empty()
                && (evidence.len() > 1 || !evidence.contains(&LINKAGE_UNSPECIFIED))
            {
                post_warning(
                    Problem::ExtraModifierFound,
                    b"FASTA-Reader: Unknown gap-type can have linkage-evidence of type 'unspecified' only.",
                    warn,
                );
            }
        }
        LinkEvid::Forbidden => {
            if !evidence.is_empty() {
                post_warning(
                    Problem::ModifierFoundButNoneExpected,
                    b"FASTA-Reader: This gap-type cannot have any linkage-evidence specified, so any will be ignored.",
                    warn,
                );
            }
        }
        LinkEvid::Required => {
            if evidence.is_empty() {
                evidence.insert(LINKAGE_UNSPECIFIED);
            }
            if evidence.len() == 1 && evidence.contains(&LINKAGE_UNSPECIFIED) {
                post_warning(
                    Problem::ExpectedModifierMissing,
                    b"CFastaReader: This gap-type should have at least one specified linkage-evidence.",
                    warn,
                );
            }
        }
    }
}
