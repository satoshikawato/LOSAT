//! NCBI's `CFastaReader` with the flags of `CFastaBioseqSource` (makeblastdb), of
//! `CMaskFastaReader` (dustmasker, segmasker, windowmasker), of convert2blastmask's
//! `CMaskFromFasta` or of makembindex's `CSequenceIStreamFasta`: [`ReaderFlags`].
//!
//! What differs from the reader of the BLAST programs (`blastinput::fasta_reader`, with
//! `fDLOptional`, `fNoSplit` and `fHyphensIgnoreAndWarn`): a definition line is required; a
//! hyphen is a residue, or a gap in a nucleotide sequence with `fParseGaps`; the residues
//! are split and packed into a delta sequence (`pack`); the definition line is kept in a
//! user object for the database's title (`title`); lowercase masks are recorded only for
//! convert2blastmask and makembindex (`ReaderFlags::lower_case_masks`: the others never call
//! `SaveMask`, so `m_CurrentMask` stays empty).

use std::io::Read;

use super::defline_ids;
use super::line_reader::BufferedLineReader;
use super::pack::{self, Coding, Segment};
use crate::blastinput::fasta_reader::reader::{
    canonicalize_string, convert_bad_indexes_to_string, gap_type_info, linkage_evidence_value,
    seq_data_in_title_warning, GapType, LinkEvid, Problem, LINKAGE_UNSPECIFIED,
};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::ncbi_str::{
    is_blank, is_space, printable_string, string_to_uint, truncate_spaces, truncate_spaces_begin,
};
use crate::seqid::accguide::{identify_seq_id, parse_flags, SeqIdWarning, F_NUC, F_PROT};
use crate::seqid::model::{ObjectId, ParsedSeqId};
use crate::seqid::rank::{best_rank_score, find_best};
use crate::utils::byte_search::find_first;
use crate::writedb::bioseq::SeqData;

/// A warning of the reader, as NCBI posts it for a reader without a listener (neither
/// `CFastaBioseqSource` nor `CMaskFastaReader` gives the reader one). Whether and how a line
/// is written is the application's: its diagnostic post level and prefix.
#[derive(Debug)]
pub(crate) enum ReaderWarning<'a> {
    /// The warning logged with `LOG_POST_X` (subcode 1) from `pLineExpt->Message()` in
    /// `CFastaReader::PostWarning` (cited at `post_warning`) and in `s_PostWarning`
    /// (`defline_ids::process_ids`): the message, which NCBI writes on a line of its own.
    Line(&'a [u8]),
    /// The warning posted with `ERR_POST_X` (subcode 7) from `e.what()` and the other warnings of
    /// `CSeq_id::ParseIDs` (Seq_id.cpp, `NCBI_ERRCODE_X` 1306), written with the application's
    /// prefix (`defline_ids::post_seq_id_warning`).
    SeqId(SeqIdWarning),
}

/// Where the reader posts its warnings.
pub(crate) type ReaderWarn<'a> = &'a mut dyn FnMut(ReaderWarning<'_>);

/// The problems that `CFastaBioseqSource` tells the reader to ignore.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:735-737
const FASTA_BIOSEQ_SOURCE_IGNORED_PROBLEMS: [Problem; 3] = [
    Problem::ModifierFoundButNoneExpected,
    Problem::TooManyAmbiguousResidues,
    Problem::TooLong,
];

/// The molecule flags of a caller: `fAssumeNuc` or `fAssumeProt` with `fForceType`
/// (makeblastdb, the mask reader: every record is of the input's molecule), or none of the
/// three (convert2blastmask: the ids of a definition line, else the residues, tell it,
/// `AssignMolType`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Molecule {
    Nucleotide,
    Protein,
    Guessed,
}

/// The `CFastaReader::EFlags` of a caller, the problems it tells the reader to ignore, and
/// whether the reader stops at makeblastdb's length limit. Where ids are parsed, makeblastdb
/// and the mask reader set `fParseRawID` (`-parse_seqids`; makeblastdb's `long_ids` is false),
/// and `fParseGaps` is set only for nucleotide input.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ReaderFlags {
    /// `fAssumeNuc` or `fAssumeProt` with `fForceType`, or none of them.
    pub(crate) molecule: Molecule,
    /// `fNoParseID` is not set: the ids of each definition line are parsed.
    pub(crate) parse_ids: bool,
    /// `fParseRawID`: an id without a `|` is parsed as an accession (`fParse_RawText`), not as
    /// a local id.
    pub(crate) parse_raw_id: bool,
    /// `SaveMask` before each record (`m_CurrentMask` is set, so lower case letters open
    /// intervals): the reader keeps the masked intervals of each record
    /// (`FastaReader::take_masks`). convert2blastmask's reader and makembindex's FASTA input
    /// (`dbindex::stream`).
    pub(crate) lower_case_masks: bool,
    /// The reader of convert2blastmask (`CMaskFromFasta`, a subclass, with
    /// `lower_case_masks`): it also reads a data line that starts with a digit as an interval.
    pub(crate) mask_from_fasta: bool,
    /// `fParseGaps`: a hyphen of a nucleotide sequence is a gap, and `>?` lines stay gaps.
    pub(crate) parse_gaps: bool,
    /// `fRequireID`: a definition line without ids is an error (with `parse_ids`).
    pub(crate) require_id: bool,
    /// `fQuickIDCheck`: only the first byte of a local id is checked.
    pub(crate) quick_id_check: bool,
    /// `fDisableNoResidues`: a record without residues is no error.
    pub(crate) disable_no_residues: bool,
    /// `IgnoreProblem`: the warnings that are not posted.
    pub(crate) ignored_problems: &'static [Problem],
    /// A record above 2,147,483,647 residues is returned without its residues
    /// (`SeqInst::TooLong`): `CBuildDatabase::AddSequences` rejects it by its length.
    pub(crate) build_db_length_limit: bool,
}

impl ReaderFlags {
    /// The flags of `CFastaBioseqSource` of the file, `is_protein` and `parse_ids`, with `long_ids`
    /// false.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:701-737
    pub(crate) const fn fasta_bioseq_source(protein: bool, parse_ids: bool) -> Self {
        Self {
            molecule: if protein {
                Molecule::Protein
            } else {
                Molecule::Nucleotide
            },
            parse_ids,
            parse_raw_id: parse_ids,
            lower_case_masks: false,
            mask_from_fasta: false,
            parse_gaps: !protein,
            require_id: parse_ids,
            quick_id_check: true,
            disable_no_residues: true,
            ignored_problems: &FASTA_BIOSEQ_SOURCE_IGNORED_PROBLEMS,
            build_db_length_limit: true,
        }
    }

    /// The flags of `CMaskFastaReader` of the input, `is_nucl` and `parse_seqids`: no `fParseGaps`,
    /// no `fRequireID`, no `fQuickIDCheck`, no `fDisableNoResidues`, and no ignored problem.
    /// `fOneSeq` makes `ReadSet` read one record, which is the caller's loop
    /// (`seqmasks_io::fasta_reader`).
    ///
    /// NCBI reference (598d8ae6): c++/include/objtools/seqmasks_io/mask_fasta_reader.hpp:57-94
    /// No length limit applies: the masking programs read a record of any length.
    pub(crate) const fn mask_fasta_reader(protein: bool, parse_seqids: bool) -> Self {
        Self {
            molecule: if protein {
                Molecule::Protein
            } else {
                Molecule::Nucleotide
            },
            parse_ids: parse_seqids,
            parse_raw_id: parse_seqids,
            lower_case_masks: false,
            mask_from_fasta: false,
            parse_gaps: false,
            require_id: false,
            quick_id_check: false,
            disable_no_residues: false,
            ignored_problems: &[],
            build_db_length_limit: false,
        }
    }

    /// The flags of `CMaskFromFasta` of the input and `parse_seqids` (convert2blastmask):
    /// `fNoParseID` without `-parse_seqids`, else none: no `fForceType`, no `fAssumeNuc` or
    /// `fAssumeProt` and no `fParseRawID`.
    ///
    /// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:56-57
    pub(crate) const fn mask_from_fasta(parse_seqids: bool) -> Self {
        Self {
            molecule: Molecule::Guessed,
            parse_ids: parse_seqids,
            parse_raw_id: false,
            lower_case_masks: true,
            mask_from_fasta: true,
            parse_gaps: false,
            require_id: false,
            quick_id_check: false,
            disable_no_residues: false,
            ignored_problems: &[],
            build_db_length_limit: false,
        }
    }

    /// The flags of makembindex's FASTA input (`CSequenceIStreamFasta`): `fForceType`,
    /// `fAssumeNuc`, `fNoParseID` and `fAllSeqIds` (no longer used by `CFastaReader`), with
    /// `SaveMask` before each record (`lower_case_masks`). They read as dustmasker's
    /// nucleotide reader without `-parse_seqids` (`mask_fasta_reader`: its `fOneSeq` acts on
    /// `ReadSet` only, and `fParseRawID` on ids that are not parsed).
    // NCBI reference (598d8ae6): c++/include/objtools/readers/fasta.hpp:90-93
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/sequence_istream_fasta.cpp:47-52
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/sequence_istream_fasta.cpp:107-113
    pub(crate) const fn dbindex_fasta() -> Self {
        Self {
            lower_case_masks: true,
            ..Self::mask_fasta_reader(false, false)
        }
    }
}

/// `CObjReaderParseException::EErrCode` of the exceptions that the reader throws.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ErrCode {
    Format,
    Eof,
    NoDefline,
    NoIds,
    Ambiguous,
    IdTooLong,
    NoResidues,
    InvalidId,
}

impl ErrCode {
    /// `GetType()::GetErrCodeString()`.
    ///
    /// NCBI reference (598d8ae6): c++/include/objtools/readers/reader_exception.hpp:64-80
    const fn kind(self) -> &'static str {
        match self {
            ErrCode::Format => "CObjReaderParseException::eFormat",
            ErrCode::Eof => "CObjReaderParseException::eEOF",
            ErrCode::NoDefline => "CObjReaderParseException::eNoDefline",
            ErrCode::NoIds => "CObjReaderParseException::eNoIDs",
            ErrCode::Ambiguous => "CObjReaderParseException::eAmbiguous",
            ErrCode::IdTooLong => "CObjReaderParseException::eIDTooLong",
            ErrCode::NoResidues => "CObjReaderParseException::eNoResidues",
            ErrCode::InvalidId => "CObjReaderParseException::eInvalidID",
        }
    }
}

/// Where NCBI throws a `CObjReaderParseException` of the reader. The places, and the
/// functions that the reports name, are those of segmasker 2.17.0's reports (`what()`;
/// probes of 2026-11-12, `tests/dustmasker/reader/ncbi/seg_*.err`); dustmasker's report
/// prints neither.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Thrower {
    /// `FASTA_ERROR` without a listener: the throw of `CFastaReader::PostWarning`, whose
    /// report names no function.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2215-2218
    PostWarning,
    /// `NCBI_THROW2` of `CFastaReader::ReadOneSeq` (fasta.cpp:393, cited at
    /// `FastaReader::read_one_seq`).
    ReadOneSeq,
    /// `NCBI_THROW2` of `CFastaDeflineReader::ParseDefline` (fasta_reader_utils.cpp:195,
    /// cited at `RecordState::parse_def_line`).
    ParseDefline,
    /// `s_PostError` without a listener (`defline_ids`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:77-82
    PostError,
    /// `NCBI_THROW2` of `CFastaReader::AssignMolType` (fasta.cpp:1606, cited at
    /// `RecordState::assign_mol_type`).
    AssignMolType,
}

/// A `CObjReaderParseException` of the reader.
#[derive(Debug)]
pub(crate) struct ParseError {
    /// `GetMsg()`.
    pub(crate) message: String,
    pub(crate) code: ErrCode,
    /// `m_Pos`: the line number that NCBI passes, 0 for `ParseDefline`'s throws.
    pub(crate) pos: u64,
    pub(crate) thrower: Thrower,
}

impl ParseError {
    pub(crate) fn new(thrower: Thrower, code: ErrCode, message: String, pos: u64) -> Self {
        Self {
            message,
            code,
            pos,
            thrower,
        }
    }

    /// `FASTA_ERROR` with the line, the message and the code, of a reader without a listener:
    /// `PostWarning` throws, with the line as the position.
    pub(crate) fn fasta_error(code: ErrCode, message: String, line: u64) -> Self {
        Self::new(Thrower::PostWarning, code, message, line)
    }

    /// The exception as NCBI reports it: `CObjReaderParseException` derives from
    /// `CParseTemplException`, whose report adds the position (`m_Pos = <pos>`).
    pub(crate) fn exception(&self) -> NcbiException {
        let (file, line, function) = match self.thrower {
            Thrower::PostWarning => ("src/objtools/readers/fasta.cpp", 2217, None),
            Thrower::ReadOneSeq => (
                "src/objtools/readers/fasta.cpp",
                393,
                Some("ncbi::objects::CFastaReader::ReadOneSeq"),
            ),
            Thrower::ParseDefline => (
                "src/objtools/readers/fasta_reader_utils.cpp",
                195,
                Some("ncbi::objects::CFastaDeflineReader::ParseDefline"),
            ),
            Thrower::PostError => (
                "src/objtools/readers/fasta_reader_utils.cpp",
                77,
                Some("ncbi::objects::s_PostError"),
            ),
            Thrower::AssignMolType => (
                "src/objtools/readers/fasta.cpp",
                1606,
                Some("ncbi::objects::CFastaReader::AssignMolType"),
            ),
        };
        let site = ThrowSite {
            file,
            line,
            kind: self.code.kind(),
            function,
        };
        NcbiException::with_pos(site, self.message.as_str(), self.pos as usize)
    }
}

/// What ends `ReadOneSeq` without a record.
#[derive(Debug)]
pub(crate) enum ReadError {
    /// The reader's own exception.
    Parse(ParseError),
    /// Another exception of NCBI's thrown while a record is read: the `CStringException` of
    /// an interval line of `CMaskFromFasta` (`ReaderFlags::mask_from_fasta`).
    Ncbi(Box<NcbiException>),
    /// An allocation whose size the input states could not be made
    /// (`std::bad_alloc` in NCBI, which no reader catches).
    OutOfMemory,
    /// The `std::out_of_range` of `basic_string::append` from position `pos` of a string of
    /// `size` (`RecordState::substituted`), which no reader catches either.
    OutOfRange { pos: u64, size: u64 },
    /// A nucleotide piece that NCBI's reader packs into a segment whose byte count wraps its
    /// `TSeqPos` to 0: `CDelta_ext_PackTarget::NewSegment` sizes an NCBI2na segment of `length`
    /// bases as `(length + 3) / 4` bytes and an NCBI4na one as `(length + 1) / 2`, so one of
    /// 4294967293 to 4294967295 bases (NCBI2na) or of 4294967295 (NCBI4na) gets 0 bytes, and
    /// the conversion writes past them (SIGSEGV in every program; survey 11d dm4, probe 2 s3n
    /// and s7n, probe 3 b4n and b5n). LOSAT's programs refuse it (`pack_overflow_what`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1456-1458,1504
    /// NCBI reference (598d8ae6): c++/src/objects/seq/Delta_ext.cpp:141-164,201-205
    PackOverflow,
}

/// What LOSAT refuses for `ReadError::PackOverflow` (before "is not supported by LOSAT's ...").
pub(crate) fn pack_overflow_what() -> String {
    "a nucleotide FASTA record with a stretch of 4,294,967,293 letters or more that NCBI packs \
     as one piece (NCBI counts the piece's bytes in 32 bits, gets 0 and crashes)"
        .to_string()
}

impl From<ParseError> for ReadError {
    fn from(error: ParseError) -> Self {
        ReadError::Parse(error)
    }
}

/// One gap of a record (`CFastaReader::SGap`): the raw residues before it and its
/// length; length 0 is a gap of unknown length (a lone hyphen at the end of a line).
/// The gap type and the linkage evidence of a `>?` line only decide whether the delta
/// literal of the gap gets a `Seq-data.gap`, which the database is built from like a
/// literal without data (`SeqInst`), and whether a gap that becomes N's or X's warns
/// (`mods`), so only that is kept; nor is the size's fuzz.
#[derive(Clone, Copy)]
struct Gap {
    position: usize,
    length: u64,
    /// A gap type other than `unknown`, or a linkage evidence other than `unspecified`
    /// alone (`gap_mod_warnings`): `AssembleSeq` warns of it when the gap becomes N's or X's.
    mods: bool,
}

/// One piece of a delta sequence (`CDelta_seq`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum DeltaSeq {
    /// `Seq-loc null`: a gap of unknown length.
    LocNull,
    /// A literal without residues (no `Seq-data`, or `Seq-data.gap`): a gap.
    Gap { length: u64 },
    /// A literal with NCBI2na or NCBI4na data: residues `start..start + length`.
    Literal(Segment),
}

/// `Seq-inst` of a record as `CFastaReader::AssembleSeq` leaves it.
#[derive(PartialEq, Eq, Debug)]
pub(crate) enum SeqInst {
    /// `repr virtual`: no residues and no gaps; length 0, no `seq-data`.
    Virtual,
    /// `repr raw`: a delta of one piece, whose `Seq-data` became the sequence's.
    Raw { length: u64, seq_data: SeqData },
    /// `repr delta` of more than one piece over the IUPACna `residues`. Past 2^32 letters
    /// (`assemble_seq`) gaps that stay pieces can make the pieces longer than `length`, and
    /// their residues repeat or skip residues where the gaps' wrapped positions do.
    Delta {
        length: u64,
        residues: Vec<u8>,
        pieces: Vec<DeltaSeq>,
    },
    /// A record longer than 2,147,483,647 residues, which `CBuildDatabase::AddSequences`
    /// rejects by its length before it reads the sequence data: LOSAT does not build it.
    /// `length` is NCBI's, the record's count modulo 2^32 (`assemble_seq`).
    TooLong { length: u64 },
}

// The text is observed, not taken from a library's source: NCBI's 2.17.0 binaries write this
// `what()` for the `std::out_of_range` of `basic_string::append` (probe of line RG item KC-11d,
// remaining-gaps-20261112/K/kc/p11d/p2).
/// The `what()` of `ReadError::OutOfRange`.
pub(crate) fn out_of_range_what(pos: u64, size: u64) -> String {
    format!("basic_string::append: __pos (which is {pos}) > this->size() (which is {size})")
}

/// A record as `CFastaReader::ReadOneSeq` returns it.
#[derive(Debug)]
pub(crate) struct Record {
    /// `Seq-id local id`, the number that `CSeqIdGenerator` gave the record (without
    /// `-parse_seqids`).
    pub(crate) local_id: u64,
    /// The ids of the record (with `-parse_seqids`).
    pub(crate) ids: Vec<ParsedSeqId>,
    /// The other id of the record's Seq-hist alignment (an id with a range).
    pub(crate) hist_old_id: Option<ParsedSeqId>,
    /// The `DefLine` field of the `CFastaReader` user object.
    pub(crate) defline_field: Vec<u8>,
    /// The title descriptor (`x_ApplyMods`), if the record has one.
    pub(crate) title_descr: Option<Vec<u8>>,
    /// `CBioseq::IsAa()`: the molecule of the input (`fForceType`), or the one
    /// `AssignMolType` found.
    pub(crate) is_aa: bool,
    pub(crate) inst: SeqInst,
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
/// The least residue value of the switch; short-run scans check the action table against it
/// directly and longer runs use `Alphabet::is_residue` (tested against it for every byte).
const RESIDUE_MIN: u8 = b'*';

/// The character switch of `ParseDataLine` for protein (`fAssumeProt`) or nucleotide
/// (`fAssumeNuc` with `fParseGaps`) input.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:807-814,822-826,856-936
/// `bIsNuc` is `fAssumeNuc` (`fForceType`), or without `fForceType` the nucleotide that the
/// ids named (`RecordState::is_nuc`). No caller sets `fLetterGaps` or
/// `fHyphensIgnoreAndWarn`: `bAllowLetterGaps` is false, and `bHyphensAreGaps` is
/// `fParseGaps`, so a hyphen is a gap with `fParseGaps` (makeblastdb's nucleotide input)
/// and a residue without it (makeblastdb's protein input, the mask reader's input).
/// A lower-case letter is stored in upper case (`eCharType_MaskedNonGap`; `OpenMask` does
/// nothing without a mask). LOSAT looks the action up in a table (speed only).
const fn char_actions(protein: bool, parse_gaps: bool) -> [u8; 256] {
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
    }
    table[b'-' as usize] = if parse_gaps { GAP } else { b'-' };
    table[b';' as usize] = COMMENT;
    let spaces = b"\t\n\x0b\x0c\r ";
    let mut index = 0;
    while index < spaces.len() {
        table[spaces[index] as usize] = SPACE;
        index += 1;
    }
    table
}

const NUCLEOTIDE_GAP_ACTIONS: [u8; 256] = char_actions(false, true);
const NUCLEOTIDE_ACTIONS: [u8; 256] = char_actions(false, false);
const PROTEIN_ACTIONS: [u8; 256] = char_actions(true, false);

/// The residues of a data line, following the protein, nucleotide-with-gaps, and
/// nucleotide-without-gaps cases of `CFastaReader::ParseDataLine`'s character switch.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:807-826
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Alphabet {
    Protein,
    NucleotideGaps,
    Nucleotide,
}

/// Whether `byte` is a letter accepted by protein input in either case, following the letter
/// cases of `CFastaReader::ParseDataLine`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:856-918
const fn is_protein_letter(byte: u8) -> bool {
    (byte & 0xdf).wrapping_sub(b'A') < 26
}

/// Whether `byte` is a nucleotide letter in either case (`ABCDGHKMNRSTUVWY`), following the
/// nucleotide branch of `CFastaReader::ParseDataLine` by the ranges of their offsets from `A`:
/// A–D, G–H, K, M–N, R–W and Y. Clearing bit 5 folds lower case to upper case without mapping
/// any non-letter into `A`–`Z`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:856-874,900-907
const fn is_nucleotide_letter(byte: u8) -> bool {
    let offset = (byte & 0xdf).wrapping_sub(b'A');
    (offset <= 3)
        | (offset.wrapping_sub(6) <= 1)
        | (offset == 10)
        | (offset.wrapping_sub(12) <= 1)
        | (offset.wrapping_sub(17) <= 5)
        | (offset == 24)
}

impl Alphabet {
    /// The character switch of the alphabet.
    const fn actions(self) -> &'static [u8; 256] {
        match self {
            Alphabet::Protein => &PROTEIN_ACTIONS,
            Alphabet::NucleotideGaps => &NUCLEOTIDE_GAP_ACTIONS,
            Alphabet::Nucleotide => &NUCLEOTIDE_ACTIONS,
        }
    }

    /// Whether `CFastaReader::ParseDataLine`'s character switch stores `byte` as a residue
    /// (`actions()[byte] >= RESIDUE_MIN`), computed by arithmetic so `find_first` can test a
    /// block of bytes at once (`tests::residue_tests_match_the_tables` checks every byte).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:856-925
    const fn is_residue(self, byte: u8) -> bool {
        match self {
            Alphabet::Protein => is_protein_letter(byte) | (byte == b'*') | (byte == b'-'),
            Alphabet::NucleotideGaps => is_nucleotide_letter(byte),
            Alphabet::Nucleotide => is_nucleotide_letter(byte) | (byte == b'-'),
        }
    }

    /// The number of consecutive residues that `s` starts with, matching the residue cases of
    /// `CFastaReader::ParseDataLine`. Short slices and the first 11 bytes use the action table
    /// directly so spaced or digit-broken lines exit without extra arithmetic per byte; longer
    /// runs use [`find_first`].
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:841-936
    fn residue_run(self, s: &[u8]) -> usize {
        let actions = self.actions();
        if s.len() < 28 {
            return s
                .iter()
                .position(|&byte| actions[byte as usize] < RESIDUE_MIN)
                .unwrap_or(s.len());
        }
        for (i, &byte) in s[..11].iter().enumerate() {
            if actions[byte as usize] < RESIDUE_MIN {
                return i;
            }
        }
        let rest = &s[11..];
        let end = match self {
            Alphabet::Protein => find_first(rest, |byte| !Alphabet::Protein.is_residue(byte)),
            Alphabet::NucleotideGaps => {
                find_first(rest, |byte| !Alphabet::NucleotideGaps.is_residue(byte))
            }
            Alphabet::Nucleotide => find_first(rest, |byte| !Alphabet::Nucleotide.is_residue(byte)),
        };
        11 + end.unwrap_or(rest.len())
    }
}

/// The bytes that `extend_upper` converts at once.
const UPPER_BLOCK: usize = 32;

/// Appends `run` to `out` in upper case (`to_ascii_uppercase` of every byte), matching the
/// upper-case storage of accepted residues in `CFastaReader::ParseDataLine`. Short runs convert
/// in a single pass; runs of at least `UPPER_BLOCK` bytes are copied whole and converted in
/// blocks of `UPPER_BLOCK`, with the final block overlapping the previous one (converting a byte
/// twice leaves it upper case).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:939-950
fn extend_upper(out: &mut Vec<u8>, run: &[u8]) {
    let len = run.len();
    if len < UPPER_BLOCK {
        out.extend(run.iter().map(u8::to_ascii_uppercase));
        return;
    }
    let start = out.len();
    out.extend_from_slice(run);
    let added = &mut out[start..];
    for block in added.as_chunks_mut::<UPPER_BLOCK>().0 {
        block.make_ascii_uppercase();
    }
    added[len - UPPER_BLOCK..].as_chunks_mut::<UPPER_BLOCK>().0[0].make_ascii_uppercase();
}

/// The state of the record being read (`ReadOneSeq` resets it).
struct RecordState {
    flags: ReaderFlags,
    /// `m_SeqData` up to `m_CurrentPos`: the raw residues, in upper case.
    seq_data: Vec<u8>,
    /// `m_Gaps`.
    gaps: Vec<Gap>,
    /// `m_CurrentGapLength`: the hyphens read since the last residue (their character,
    /// `m_CurrentGapChar`, is always `-`), a `TSeqPos` (a run of 2^32 hyphens is none).
    current_gap_length: u32,
    /// `m_TotalGapLength`, counted on (NCBI's `TSeqPos` wraps at 2^32: `assemble_seq`).
    total_gap_length: u64,
    /// `m_CurrentSeqTitles`.
    titles: Vec<LineText>,
    /// The `DefLine` field of the `CFastaReader` user object.
    defline_field: Vec<u8>,
    /// The title descriptor.
    title_descr: Option<Vec<u8>>,
    /// The local ID of the record (without `-parse_seqids`).
    local_id: u64,
    /// The Bioseq's ids (with `-parse_seqids`).
    ids: Vec<ParsedSeqId>,
    /// The id of the definition line that a range replaced (the other row of the
    /// Seq-hist alignment).
    hist_old_id: Option<ParsedSeqId>,
    /// The bad positions of a data line, reused.
    bad_positions: Vec<u32>,
    /// `Seq-inst.mol` of the record before `AssignMolType` (`ReaderFlags::molecule` is
    /// `Guessed`): `Some(true)` for a protein, `Some(false)` for a nucleotide, set by
    /// `ParseDefLine` from the ids (`xSetSeqMol`); `None` while it is not set.
    mol: Option<bool>,
    /// `m_MaskRangeStart` (`kInvalidSeqPos` is `None`).
    mask_start: Option<u32>,
    /// `CMaskFromFasta::m_mask`.
    masks: Vec<(u32, u32)>,
    /// 2^32, the range of NCBI's `TSeqPos`, in which `assemble_seq` and `close_gap` count the
    /// positions of a record (a smaller power of two in the tests of the wrap:
    /// `FastaReader::with_seq_pos_modulus`).
    seq_pos_modulus: u64,
}

/// The range of NCBI's `TSeqPos` (32 bits).
const SEQ_POS_MODULUS: u64 = 1 << 32;

/// NCBI's `CFastaReader` over a `CBufferedLineReader`.
pub(crate) struct FastaReader<R: Read> {
    lines: BufferedLineReader<R>,
    /// `CSeqIdGenerator::m_Counter`.
    id_counter: u64,
    state: RecordState,
}

/// `PostWarning` for a warning: nothing for an ignored problem, else the message.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2196-2214
/// Neither `CFastaBioseqSource::GetNext` nor `CMaskFastaReader::GetNextSequence` gives the
/// reader a listener, so a warning is the message on a line of its own (`LOG_POST` writes
/// no severity) where the application's post level lets it through.
fn post_warning(ignored: &[Problem], problem: Problem, message: &[u8], warn: ReaderWarn<'_>) {
    if !ignored.contains(&problem) {
        warn(ReaderWarning::Line(message));
    }
}

impl<R: Read> FastaReader<R> {
    /// `CFastaReader` of `*m_LineReader` and the flags (`CFastaBioseqSource`) or of a stream and
    /// the flags (`CMaskFastaReader`; `ILineReader::New` makes the same line reader, cited at
    /// `seqmasks_io::fasta_reader`) with the caller's `flags`; the IDs come from the reader's
    /// own `CSeqIdGenerator`, which counts from 1 and has no prefix.
    ///
    /// NCBI reference (598d8ae6): c++/include/objtools/readers/fasta_reader_utils.hpp:160-165
    /// The registry value `BLAST/MAX_SEQID_LENGTH` (build_db.cpp:738-747, makeblastdb's
    /// `SetMaxIDLength`) is not read: LOSAT has no `.ncbirc`, and makeblastdb refuses the
    /// setting (`check_registry_settings`). The mask reader never sets a maximum.
    pub(crate) fn new(input: R, flags: ReaderFlags) -> Self {
        Self {
            lines: BufferedLineReader::new(input),
            id_counter: 1,
            state: RecordState {
                flags,
                ids: Vec::new(),
                hist_old_id: None,
                seq_data: Vec::new(),
                gaps: Vec::new(),
                current_gap_length: 0,
                total_gap_length: 0,
                titles: Vec::new(),
                defline_field: Vec::new(),
                title_descr: None,
                local_id: 0,
                bad_positions: Vec::new(),
                mol: None,
                mask_start: None,
                masks: Vec::new(),
                seq_pos_modulus: SEQ_POS_MODULUS,
            },
        }
    }

    /// The reader with the positions of a record counted modulo `modulus`, a power of two
    /// below 2^32, so that a test shows the wrap at 2^32 on a short record.
    #[cfg(test)]
    pub(crate) fn with_seq_pos_modulus(mut self, modulus: u64) -> Self {
        assert!(modulus.is_power_of_two() && modulus <= SEQ_POS_MODULUS);
        self.state.seq_pos_modulus = modulus;
        self
    }

    /// The reader's molecule is protein (`fAssumeProt` with `fForceType`).
    pub(crate) fn is_protein(&self) -> bool {
        self.state.flags.molecule == Molecule::Protein
    }

    /// The ids of each definition line are parsed (`-parse_seqids`).
    pub(crate) fn parses_ids(&self) -> bool {
        self.state.flags.parse_ids
    }

    /// With `ReaderFlags::lower_case_masks`, the masks of the record that `read_one_seq` read
    /// last (after an error, what it had read of that record): the intervals of its lower case
    /// letters and, with `ReaderFlags::mask_from_fasta`, of its interval lines (0-based,
    /// inclusive), in the order in which they were closed or read (convert2blastmask's
    /// `CMaskFromFasta::GetMask()`, makembindex's `SaveMask` location); empty otherwise. The
    /// intervals stay with the reader (not in `Record`) so that the records of the other
    /// programs, which never have any, are not made larger.
    pub(crate) fn take_masks(&mut self) -> Vec<(u32, u32)> {
        std::mem::take(&mut self.state.masks)
    }

    /// The residues (IUPACna, upper case) that the record `read_one_seq` returned last was
    /// packed from when it is a `SeqInst::Raw` of NCBI2na or NCBI4na data, which the reader
    /// keeps until it reads the next record (a delta takes them along as its `residues`).
    /// They are handed over with their buffer, so that a masking program holds its letters
    /// without a second copy of the sequence; the next record is read into a buffer of the
    /// same capacity, as before.
    pub(crate) fn take_raw_residues(&mut self) -> Vec<u8> {
        let capacity = self.state.seq_data.len();
        std::mem::replace(&mut self.state.seq_data, Vec::with_capacity(capacity))
    }

    /// `GetLineReader().AtEOF()`.
    pub(crate) fn at_eof(&self) -> bool {
        self.lines.at_eof()
    }

    /// The error of a failed read of the input, once (the line reader then reads as at
    /// the end of the input).
    pub(crate) fn take_read_error(&mut self) -> Option<std::io::Error> {
        self.lines.take_read_error()
    }

    /// One record.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:312-440
    /// `PeekChar` of `CBufferedLineReader` does not read, so the reader is not at the end
    /// of the input after it. `FASTA_PROGRESS` posts to a listener only, and there is
    /// none. Without `fDLOptional` data before a definition line is the `eNoDefline`
    /// error. `CBadResiduesException` comes only with `fValidate`, which is not set.
    pub(crate) fn read_one_seq(&mut self, warn: ReaderWarn<'_>) -> Result<Record, ReadError> {
        debug_assert!(
            self.state.flags.lower_case_masks || !self.state.flags.mask_from_fasta,
            "convert2blastmask's reader keeps the lower case masks"
        );
        match (
            self.state.flags.lower_case_masks,
            self.state.flags.mask_from_fasta,
        ) {
            (false, _) => self.read_one_seq_of::<false, false>(warn),
            (true, false) => self.read_one_seq_of::<true, false>(warn),
            (true, true) => self.read_one_seq_of::<true, true>(warn),
        }
    }

    /// `read_one_seq` with `ReaderFlags::lower_case_masks` as `MASKS` and
    /// `ReaderFlags::mask_from_fasta` as `MASK_FROM_FASTA`: the records of the other readers
    /// take no step of the masks but `close_mask`, which has nothing open to close for them (one
    /// choice per record).
    fn read_one_seq_of<const MASKS: bool, const MASK_FROM_FASTA: bool>(
        &mut self,
        warn: ReaderWarn<'_>,
    ) -> Result<Record, ReadError> {
        let state = &mut self.state;
        state.seq_data.clear();
        state.gaps.clear();
        state.current_gap_length = 0;
        state.total_gap_length = 0;
        state.titles.clear();
        state.defline_field.clear();
        state.title_descr = None;
        state.mol = None;
        // `m_MaskRangeStart` is set to `kInvalidSeqPos`; `CMaskFromFasta::GetNextSequence` clears
        // its intervals (`m_mask.clear()`) before `ReadOneSeq`.
        state.mask_start = None;
        state.masks.clear();

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
                    let line_number = self.lines.line_number();
                    state.parse_def_line(next_line, &mut self.id_counter, line_number, warn)?;
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
            // `NStr::TruncateSpaces_Unsafe` of the line (`eTrunc_Both`, fasta.cpp:376).
            let line = truncate_spaces(raw);
            if line.is_empty() || matches!(line[0], b'!' | b'#' | b';') {
                continue;
            }
            if need_defline {
                self.lines.unget_line();
                return Err(ReadError::Parse(ParseError::new(
                    Thrower::ReadOneSeq,
                    ErrCode::NoDefline,
                    format!(
                        "CFastaReader: Input doesn't start with a defline or comment around line {line_number}"
                    ),
                    line_number,
                )));
            }
            let raw_length = raw.len();
            match line.strip_prefix(b">?_") {
                Some(rest) => state.parse_data_line::<MASKS, MASK_FROM_FASTA>(
                    &[b">".as_slice(), rest].concat(),
                    raw_length,
                    line_number,
                    warn,
                )?,
                None => state.parse_data_line::<MASKS, MASK_FROM_FASTA>(
                    line,
                    raw_length,
                    line_number,
                    warn,
                )?,
            }
        }

        if need_defline && self.lines.at_eof() {
            let line_number = self.lines.line_number();
            return Err(ReadError::Parse(ParseError::fasta_error(
                ErrCode::Eof,
                format!("CFastaReader: Expected defline around line {line_number}"),
                line_number,
            )));
        }

        let (inst, is_aa) = state.assemble_seq(self.lines.line_number(), warn)?;
        Ok(Record {
            local_id: state.local_id,
            ids: std::mem::take(&mut state.ids),
            hist_old_id: state.hist_old_id.take(),
            defline_field: std::mem::take(&mut state.defline_field),
            title_descr: state.title_descr.take(),
            is_aa,
            inst,
        })
    }
}

/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:189-207
fn is_ambig_nuc(c: u8) -> bool {
    b"URYSWKMBDHVN".contains(&c.to_ascii_uppercase())
}

/// General counting fallback for [`data_line_counts`] when `checked` contains a non-`ACGT`
/// nucleotide byte or a non-residue protein byte, following `CFastaReader::CheckDataLine`'s
/// classification loop and `;` comment stop.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:730-750
#[cold]
#[inline(never)]
fn data_line_counts_slow(checked: &[u8], is_nuc: bool) -> (usize, usize, usize) {
    if checked.len() < 32 {
        let (mut good, mut bad, mut ambig_nuc) = (0usize, 0usize, 0usize);
        for &c in checked {
            if c.is_ascii_alphabetic() || c == b'*' {
                good += 1;
                if is_nuc && is_ambig_nuc(c) {
                    ambig_nuc += 1;
                }
            } else if c == b'-' {
                good += 1;
            } else if is_space(c) || c.is_ascii_digit() {
            } else if c == b';' {
                break;
            } else {
                bad += 1;
            }
        }
        return (good, bad, ambig_nuc);
    }
    fn count(bytes: &[u8], test: impl Fn(u8) -> bool) -> usize {
        // At most 70 bytes: the count fits a byte.
        usize::from(
            bytes
                .iter()
                .fold(0u8, |count, &c| count + u8::from(test(c))),
        )
    }
    let checked = &checked[..find_first(checked, |c| c == b';').unwrap_or(checked.len())];
    let good = |c: u8| c.is_ascii_alphabetic() | (c == b'*') | (c == b'-');
    let goods = count(checked, good);
    let bads = count(checked, |c| !(good(c) | is_space(c) | c.is_ascii_digit()));
    let ambiguous = if is_nuc {
        count(checked, |c| c.is_ascii_alphabetic() && is_ambig_nuc(c))
    } else {
        0
    };
    (goods, bads, ambiguous)
}

/// Counts good bytes (letters, `*`, `-`), bad bytes (neither good, white space, nor digit), and
/// for nucleotide input ambiguous nucleotide letters over `checked` (at most 70 bytes) up to the
/// first `;`, following `CFastaReader::CheckDataLine`. A line consisting solely of `A`/`C`/`G`/`T`
/// (nucleotide, either case) or accepted protein residues returns `(checked.len(), 0, 0)` after
/// one [`find_first`] scan; any other line uses [`data_line_counts_slow`]
/// (`tests::data_line_counts_match_the_loop` compares against the byte-by-byte loop).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:719-750
fn data_line_counts(checked: &[u8], is_nuc: bool) -> (usize, usize, usize) {
    debug_assert!(checked.len() <= 70);
    if is_nuc {
        let is_acgt = |c: u8| {
            let u = c & 0xdf;
            (u == b'A') | (u == b'C') | (u == b'G') | (u == b'T')
        };
        if find_first(checked, |c| !is_acgt(c)).is_none() {
            return (checked.len(), 0, 0);
        }
    } else if find_first(checked, |c| !Alphabet::Protein.is_residue(c)).is_none() {
        return (checked.len(), 0, 0);
    }
    data_line_counts_slow(checked, is_nuc)
}

impl RecordState {
    /// `ParseDefLine`: the title for the reader's checks
    /// (`CFastaDeflineReader::ParseDefline`), the ids (with `-parse_seqids`) or the
    /// generated local ID (`PostProcessIDs` without IDs), and the user object with the
    /// definition line.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:157-180,210-226
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:566-574,591-603
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:626-628
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:466-474
    /// `defline` starts with `>` except for an empty line (`BufferedLineReader`'s module
    /// comment), which returns before the test of its first character. With `fForceType`
    /// the molecule is not set here; without it (`Molecule::Guessed`) the ids set it
    /// (`xSetSeqMol`, before a ranged id is replaced). Without `-parse_seqids` no range is
    /// parsed (`fNoParseID`).
    fn parse_def_line(
        &mut self,
        defline: &[u8],
        id_counter: &mut u64,
        line_number: u64,
        warn: ReaderWarn<'_>,
    ) -> Result<(), ParseError> {
        self.titles.clear();
        self.ids.clear();
        self.hist_old_id = None;
        let len = defline.len();
        let mut has_range = false;
        if len > 1 && !is_blank(&defline[1..]) {
            let start = 1 + defline[1..]
                .iter()
                .position(|&byte| !is_space(byte))
                .unwrap_or(len - 1);
            let mut title_start = start;
            if self.flags.parse_ids {
                // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:181-207
                // `defline[pos]` is a `char`, signed on the platforms makeblastdb is
                // compared with (x86-64): a byte from 0x80 up is below ' ' and ends the id.
                let mut pos = start;
                while pos < len && (defline[pos] as i8) > (b' ' as i8) {
                    pos += 1;
                }
                let (range_len, _, _) = defline_ids::parse_range(&defline[start..pos]);
                let id_string = &defline[start..pos - range_len];
                if is_blank(id_string) {
                    return Err(ParseError::new(
                        Thrower::ParseDefline,
                        ErrCode::Format,
                        "Unable to locate sequence id in definition line".to_string(),
                        0,
                    ));
                }
                title_start = pos;
                self.ids = defline_ids::process_ids(
                    id_string,
                    self.flags.molecule == Molecule::Protein,
                    self.flags.quick_id_check,
                    self.flags.parse_raw_id,
                    line_number,
                    warn,
                )?;
                has_range = range_len > 0;
            }
            while title_start < len && is_space(defline[title_start]) {
                title_start += 1;
            }
            if title_start < len {
                let end = title_start
                    + 1
                    + find_first(&defline[title_start + 1..], |byte| byte < b' ')
                        .unwrap_or(len - title_start - 1);
                self.titles.push(LineText {
                    text: defline[title_start..end].to_vec(),
                });
            }
        }
        if self.flags.parse_ids {
            // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:576-583
            // The ids are empty only for a definition line that is blank after its `>`
            // (`ParseDefline` returns before it reads ids); without `fRequireID` (the mask
            // reader) `PostProcessIDs` generates the record's id.
            //
            // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:584-589
            if !self.ids.is_empty() && self.flags.molecule == Molecule::Guessed {
                self.mol = mol_of_ids(&self.ids, warn);
            }
            if self.ids.is_empty() {
                if self.flags.require_id {
                    return Err(ParseError::fasta_error(
                        ErrCode::NoIds,
                        format!(
                            "CFastaReader: Defline lacks a proper ID around line {line_number}"
                        ),
                        line_number,
                    ));
                }
                self.ids
                    .push(ParsedSeqId::Local(ObjectId::Id(*id_counter as i32)));
                *id_counter += 1;
            }
            // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:619-644
            // The alignment is kept as its other id (`CBioseq::GetNonLocalId` reads it).
            if has_range {
                let old_id = find_best(&self.ids, best_rank_score).cloned();
                self.ids.clear();
                self.ids
                    .push(ParsedSeqId::Local(ObjectId::Id(*id_counter as i32)));
                *id_counter += 1;
                self.hist_old_id = old_id;
            }
        } else {
            self.local_id = *id_counter;
            *id_counter += 1;
        }
        self.defline_field = printable_string(defline).into_owned();
        Ok(())
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:710-771
    /// `bIgnoreHyphens` is false: a hyphen counts as good. `bIsNuc` is true for nucleotide
    /// input (`fAssumeNuc` with `fForceType`). The warning about ambiguous nucleotides is an
    /// ignored problem of makeblastdb.
    #[inline(always)]
    fn check_data_line(
        &self,
        s: &[u8],
        line_number: u64,
        warn: ReaderWarn<'_>,
    ) -> Result<(), ParseError> {
        if !self.seq_data.is_empty() {
            return Ok(());
        }
        let len_to_check = s.len().min(70);
        let (good, bad, ambig_nuc) = data_line_counts(&s[..len_to_check], self.is_nuc());
        if bad >= good / 3 && (len_to_check > 3 || good == 0 || bad > good) {
            return Err(ParseError::fasta_error(
                ErrCode::Format,
                format!(
                    "CFastaReader: Near line {line_number}, there's a line that doesn't look like plausible data, but it's not marked as defline or comment."
                ),
                line_number,
            ));
        }
        const WARN_PERCENT_AMBIGUOUS: usize = 40;
        let percent_ambig = if good == 0 {
            100
        } else {
            ambig_nuc * 100 / good
        };
        if len_to_check > 3 && percent_ambig > WARN_PERCENT_AMBIGUOUS {
            post_warning(
                self.flags.ignored_problems,
                Problem::TooManyAmbiguousResidues,
                format!(
                    "FASTA-Reader: Start of first data line in seq is about {percent_ambig}% ambiguous nucleotides (shouldn't be over {WARN_PERCENT_AMBIGUOUS}%)"
                )
                .as_bytes(),
                warn,
            );
        }
        Ok(())
    }

    /// `s` is the line without white space at its ends, `raw_length` the length of the
    /// line as read (`(*GetLineReader()).length()`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:774-782,816-816,938-1016
    /// The character switch is cited at `char_actions`. `m_SeqData` is resized to
    /// `m_CurrentPos` plus `s_len`, which makes the string's length the residues before the line
    /// plus the line's length while the line is read; `x_CloseGap` reads that length (`close_gap`).
    /// The fast copy of `fSkipCheck` is not taken (no `fSkipCheck`). `x_NucOrProt()` is empty: the
    /// molecule is set only in `AssembleSeq`.
    fn parse_data_line<const MASKS: bool, const MASK_FROM_FASTA: bool>(
        &mut self,
        s: &[u8],
        raw_length: usize,
        line_number: u64,
        warn: ReaderWarn<'_>,
    ) -> Result<(), ReadError> {
        // `CMaskFromFasta::ParseDataLine`: an interval line is read as the residue `A`.
        let s = if MASK_FROM_FASTA && s[0].is_ascii_digit() {
            self.parse_interval(s)?;
            &b"A"[..]
        } else {
            s
        };
        if s.starts_with(b">?") {
            self.parse_gap_line(s, line_number, warn);
            return Ok(());
        }
        self.check_data_line(s, line_number, warn)?;

        // `fParseGaps` is set only for nucleotide input (`ReaderFlags`).
        let alphabet = if !self.is_nuc() {
            Alphabet::Protein
        } else if self.flags.parse_gaps {
            Alphabet::NucleotideGaps
        } else {
            Alphabet::Nucleotide
        };
        let actions = alphabet.actions();
        let seq_data_length = self.seq_data.len() + s.len();
        self.bad_positions.clear();
        let mut pos = 0;
        while pos < s.len() {
            let c = s[pos];
            match actions[c as usize] {
                SPACE => {}
                COMMENT => break,
                GAP => {
                    if MASKS {
                        self.close_mask();
                    }
                    let run = 1 + s[pos + 1..]
                        .iter()
                        .position(|&byte| byte != c)
                        .unwrap_or(s.len() - pos - 1);
                    // `m_CurrentGapLength` grows by `pos2` minus `pos`, a `TSeqPos`.
                    self.current_gap_length = self.current_gap_length.wrapping_add(run as u32);
                    pos += run - 1;
                }
                BAD => self.bad_positions.push(pos as u32),
                residue => {
                    self.close_gap(pos == 0, seq_data_length, raw_length);
                    // The residues that follow are stored at once (LOSAT's speed: no gap
                    // is open after the first one). The switch stores a residue in upper case
                    // (`tests::residue_tests_match_the_tables`).
                    let run = 1 + alphabet.residue_run(&s[pos + 1..]);
                    debug_assert!(run >= 1);
                    if MASKS {
                        self.mask_run(&s[pos..pos + run]);
                    }
                    if run == 1 {
                        self.seq_data.push(residue);
                    } else if run < UPPER_BLOCK {
                        self.seq_data
                            .extend(s[pos..pos + run].iter().map(|&byte| actions[byte as usize]));
                    } else {
                        extend_upper(&mut self.seq_data, &s[pos..pos + run]);
                    }
                    pos += run - 1;
                }
            }
            pos += 1;
        }
        // NCBI resizes `m_SeqData` to the `TSeqPos` position `m_CurrentPos` (fasta.cpp:985): past
        // 2^32 residues the string keeps those written since the position last wrapped.
        if self.seq_data.len() as u64 >= self.seq_pos_modulus {
            self.wrap_seq_data();
        }

        if !self.bad_positions.is_empty() {
            let mut message = b"FASTA-Reader: Ignoring invalid residues at position(s): ".to_vec();
            convert_bad_indexes_to_string(&mut message, line_number, &self.bad_positions);
            post_warning(
                self.flags.ignored_problems,
                Problem::InvalidResidue,
                &message,
                warn,
            );
        }
        Ok(())
    }

    /// The residues of a record past `seq_pos_modulus` (2^32) of them: those since the count
    /// last wrapped, the last (count modulo 2^32), as NCBI's `m_SeqData` after a data line.
    /// `CheckDataLine` then tests the next line again where none are left, as NCBI's does
    /// for an empty `m_SeqData`.
    #[cold]
    fn wrap_seq_data(&mut self) {
        let kept = (self.seq_data.len() as u64 % self.seq_pos_modulus) as usize;
        self.seq_data.drain(..self.seq_data.len() - kept);
    }

    /// `bIsNuc` of `ParseDataLine` and of `CheckDataLine`: `fAssumeNuc` (with `fForceType`),
    /// or, without `fForceType`, a `Seq-inst.mol` of `na` that the ids set (`ParseDefLine`).
    /// The two tests agree for every caller's flags (cited at `check_data_line` and at the
    /// action tables).
    fn is_nuc(&self) -> bool {
        match self.flags.molecule {
            Molecule::Nucleotide => true,
            Molecule::Protein => false,
            Molecule::Guessed => self.mol == Some(false),
        }
    }

    /// `GetCurrentPos` with `ePosWithGapsAndSegs` (no `fInSegSet`, so `m_SegmentBase` is 0), a
    /// `TSeqPos`.
    fn current_pos_with_gaps(&self) -> u32 {
        (self.seq_data.len() as u64 + self.total_gap_length) as u32
    }

    /// `OpenMask` and `CloseMask` of the residues of `run`, which are about to be stored: a
    /// lower case letter is `eCharType_MaskedNonGap` (`OpenMask`), any other residue
    /// `eCharType_NormalNonGap` (`CloseMask`, before the position moves past it).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:939-955
    /// NCBI reference (598d8ae6): c++/include/objtools/readers/fasta.hpp:494-500
    /// `CMaskFromFasta` keeps its own copy of each interval (`x_OpenMask`, `x_CloseMask`,
    /// cited at `close_mask`). Only `parse_data_line` with `MASKS` (convert2blastmask,
    /// makembindex) calls it, once for each run of residues, so it is not cold: the other
    /// readers' instance has no call.
    fn mask_run(&mut self, run: &[u8]) {
        let base = self.current_pos_with_gaps();
        for (i, byte) in run.iter().enumerate() {
            let pos = base.wrapping_add(i as u32);
            if byte.is_ascii_lowercase() {
                // `m_MaskRangeStart` is `kInvalidSeqPos` (4294967295) while no interval is
                // open: one opened there stays closed (`x_OpenMask` stores that value).
                if self.mask_start.is_none() && pos != u32::MAX {
                    self.mask_start = Some(pos);
                }
            } else if let Some(start) = self.mask_start.take() {
                self.masks.push((start, pos.wrapping_sub(1)));
            }
        }
    }

    /// `CloseMask()`: the open interval ends before the current position.
    ///
    /// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:94-105
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1079-1092
    /// (convert2blastmask keeps its own copy of the interval; `CFastaReader::x_CloseMask` adds
    /// the same interval to `m_CurrentMask`, the location of `SaveMask` that makembindex
    /// reads.)
    fn close_mask(&mut self) {
        if let Some(start) = self.mask_start.take() {
            let end = self.current_pos_with_gaps().wrapping_sub(1);
            self.masks.push((start, end));
        }
    }

    /// `CMaskFromFasta::ParseDataLine` of a line that starts with a digit: the interval
    /// `from - to` (`parse_data_line` then reads a residue `A` in its place).
    ///
    /// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:76-92
    /// The `find` of a hyphen in `s` is never 0 after a digit (`npos` without one), so every such
    /// line is an interval: without `-` the second number is the empty string, an error. GCC
    /// evaluates the two conversions from the right: the second number's error is the one reported
    /// (convert2blastmask 2.17.0: `1abc` reports `''`, `0x10 - 5` reports `'0x10'`). `s` is the
    /// line without its spaces at both ends (`ReadOneSeq`, cited at `FastaReader::read_one_seq`).
    #[cold]
    fn parse_interval(&mut self, s: &[u8]) -> Result<(), ReadError> {
        let (s1, s2) = match s.iter().position(|&byte| byte == b'-') {
            Some(dash) => (&s[..dash], &s[dash + 1..]),
            None => (s, &[][..]),
        };
        let to = string_to_uint(truncate_spaces(s2))
            .map_err(|error| ReadError::Ncbi(Box::new(error)))?;
        let from = string_to_uint(truncate_spaces(s1))
            .map_err(|error| ReadError::Ncbi(Box::new(error)))?;
        self.masks.push((from, to));
        Ok(())
    }

    /// `AssignMolType` without `fForceType`, `fAssumeNuc` or `fAssumeProt`: whether the record
    /// is a protein (`IsAa()`; a `mol` of `not-set`, for a record without residues, is not).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1575-1610
    /// (`sRefineNaMol` makes an `na` a `dna` or an `rna`, which is not a protein either; no
    /// caller sets `fStrictGuess` or `fLaxGuess`.)
    fn assign_mol_type(&self, line_number: u64) -> Result<bool, ParseError> {
        if let Some(is_aa) = self.mol {
            return Ok(is_aa);
        }
        if self.seq_data.is_empty() {
            return Ok(false);
        }
        let length = self.seq_data.len().min(4096);
        match sequence_type(&self.seq_data[..length]) {
            Some(is_aa) => Ok(is_aa),
            None => Err(ParseError::new(
                Thrower::AssignMolType,
                ErrCode::Ambiguous,
                format!(
                    "CFastaReader: Unable to determine sequence type (is it nucleotide? protein?) around line {line_number}"
                ),
                line_number,
            )),
        }
    }

    /// `CloseGap` with `atStartOfLine`: ends the run of hyphens, if one is open.
    /// `seq_data_length` is `m_SeqData.length()` (see `parse_data_line`; outside a data
    /// line it is the number of residues), `line_length` the length of the reader's
    /// current line.
    ///
    /// NCBI reference (598d8ae6): c++/include/objtools/readers/fasta.hpp:284-289
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1018-1077
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
        let mut len = u64::from(self.current_gap_length);
        let pos = self.seq_data.len();
        // `l` is the `TSeqPos` of `m_SeqData.length()` and `pos` the wrapped `m_CurrentPos`;
        // `pos + length()` is a `size_t` sum.
        let l = seq_data_length as u64 % self.seq_pos_modulus;
        let pos_wrapped = pos as u64 % self.seq_pos_modulus;
        if len == 1
            && (l == pos_wrapped || (l == pos_wrapped + line_length as u64 && at_start_of_line))
        {
            len = 0;
        }
        self.gaps.push(Gap {
            position: pos,
            length: len,
            mods: false,
        });
        self.total_gap_length += len;
        self.current_gap_length = 0;
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1094-1184,1342-1349
    /// The modifiers (cited at `gap_mod_warnings`) only warn here (`Gap`). The caller's
    /// default gap type (`m_gap_type`) is set neither by `CFastaBioseqSource` nor by
    /// `CMaskFastaReader`.
    fn parse_gap_line(&mut self, line: &[u8], line_number: u64, warn: ReaderWarn<'_>) {
        self.close_gap(false, self.seq_data.len(), 0);

        let mut remaining = truncate_spaces(&line[2..]);
        let position = self.seq_data.len();
        if let Some(rest) = remaining.strip_prefix(b"unk") {
            remaining = truncate_spaces_begin(rest);
        }
        let digits = remaining
            .iter()
            .position(|byte| !byte.is_ascii_digit())
            .unwrap_or(remaining.len());
        // `NStr::StringToUInt` of `sDigits` with `NStr::fConvErr_NoThrow`: 0 for an empty string
        // and for a value above 4,294,967,295.
        let mut gap_size = string_to_uint(&remaining[..digits]).unwrap_or(0);
        if gap_size == 0 {
            post_warning(
                self.flags.ignored_problems,
                Problem::NonPositiveLength,
                format!("CFastaReader: Bad gap size at line {line_number}").as_bytes(),
                warn,
            );
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
                post_warning(
                    self.flags.ignored_problems,
                    Problem::ParsingModifiers,
                    format!("CFastaReader: Problem parsing gap mods at line {line_number}")
                        .as_bytes(),
                    warn,
                );
                break;
            };
            let key = truncate_spaces(&remaining[1..equal]).to_vec();
            let value = truncate_spaces(&remaining[equal + 1..close]).to_vec();
            mods.push((key, value));
            remaining = truncate_spaces_begin(&remaining[close + 1..]);
        }
        // `multimap` iterates in key order (the order of insertion for equal keys).
        mods.sort_by(|left, right| left.0.cmp(&right.0));
        let specified = gap_mod_warnings(&mods, line_number, self.flags.ignored_problems, warn);

        self.gaps.push(Gap {
            position,
            length: u64::from(gap_size),
            mods: specified,
        });
        self.total_gap_length += u64::from(gap_size);
    }

    /// `AssembleSeq`: the sequence of the record.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1351-1373,1386-1391,1411-1461
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1463-1471,1497-1516
    /// The molecule is the input's (`AssignMolType` with `fForceType`), or without it the
    /// one that the ids named or the letters show (`assign_mol_type`):
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1569-1573
    /// Input without `fParseGaps` (protein input, and every input of the mask reader) turns
    /// its `>?` gaps into runs of `X` (protein) or `N` (nucleotide), so it has no gaps left.
    /// The residues of a protein are one NCBIeaa literal, because `CSeqConvert::Pack` copies
    /// NCBIeaa when the target takes no gaps, and that one piece becomes the sequence:
    ///
    /// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:1087-1092,1108-1112
    /// Nucleotide input keeps its gaps as delta pieces (a gap of length 0 as `Seq-loc
    /// null`), and each run of residues between them is packed into NCBI2na and NCBI4na
    /// literals (`pack::pack_iupacna`). A delta of one piece becomes a raw sequence with
    /// that piece's `Seq-data`: `SetLiteral().SetSeq_data()` of a gap piece creates a
    /// `Seq-data` without a choice when it has none, so a record that is one gap has a
    /// `Seq-data` that the database writer cannot read (`SeqData::Other`).
    ///
    /// The run of `X` or `N` is as long as the `>?` lines say, up to 2,147,483,647 residues
    /// for a few bytes of input (makeblastdb's limit; the mask reader's is 4,294,967,295, the
    /// most a wrapped length can be: `substituted`). NCBI
    /// reserves the whole string at once and a failure is a `std::bad_alloc` (makeblastdb:
    /// "BLAST ran out of memory"); here it is the error of this function
    /// (`ReadError::OutOfMemory`), where an allocation that aborts the process would lose
    /// that outcome.
    ///
    /// A record without residues and gaps is the `eNoResidues` error at `line_number` (the
    /// reader's `LineNumber()`) unless the caller sets `fDisableNoResidues`; with it, the
    /// record is virtual (makeblastdb).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1388-1391
    fn assemble_seq(
        &mut self,
        line_number: u64,
        warn: ReaderWarn<'_>,
    ) -> Result<(SeqInst, bool), ReadError> {
        self.close_gap(false, self.seq_data.len(), 0);
        self.close_mask();
        // NCBI keeps the positions of a record in `TSeqPos` (`m_CurrentPos`, `SGap::m_uPos`,
        // `m_TotalGapLength`), which wrap at 2^32 (`RecordState::seq_pos_modulus`); past 2^32
        // residues `m_SeqData` holds the last (count modulo 2^32), which all that follows
        // reads (`wrap_seq_data`, after each data line).
        //
        // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:816,863-865,985
        let modulus = self.seq_pos_modulus;
        let is_aa = match self.flags.molecule {
            Molecule::Protein => true,
            Molecule::Nucleotide => false,
            Molecule::Guessed => self.assign_mol_type(line_number)?,
        };

        let titles = std::mem::take(&mut self.titles);
        for title in &titles {
            self.parse_title(&title.text, warn);
        }

        // `GetCurrentPos` with `ePosWithGaps`, `m_CurrentPos` plus `m_TotalGapLength`: the record's
        // length modulo 2^32, which makeblastdb tests and the `Seq-inst` takes.
        let current = self.seq_data.len() as u64 % modulus;
        let total_gap_length = self.total_gap_length % modulus;
        let length = (current + total_gap_length) % modulus;

        // The substituted string holds the residues and the runs, never empty; its pieces are
        // those of its first `length` letters, and none for a length that wrapped to 0. NCBI
        // substitutes before it packs the data and before makeblastdb tests the length, so the
        // substitution's warnings and `std::out_of_range` come first, then the packing's crash
        // (`ReadError::PackOverflow`), which only a nucleotide piece of `modulus - 3` letters
        // or more can meet. A record that ends at the length keeps no letters here unless its
        // packing is to be checked.
        let substituted = !self.flags.parse_gaps && total_gap_length > 0;
        let packed = if substituted { length } else { current };
        let check_packing = !is_aa && packed + 3 >= modulus;
        let too_long = self.flags.build_db_length_limit && length > i32::MAX as u64;
        if substituted {
            let keep = if too_long && !check_packing {
                0
            } else {
                length
            };
            self.seq_data = self.substituted(keep, is_aa, warn)?;
            self.gaps.clear();
            self.total_gap_length = 0;
        }
        if too_long && !check_packing {
            return Ok((SeqInst::TooLong { length }, is_aa));
        }
        if !substituted && self.gaps.is_empty() && self.seq_data.is_empty() {
            if !self.flags.disable_no_residues {
                return Err(ReadError::Parse(ParseError::fasta_error(
                    ErrCode::NoResidues,
                    "FASTA-Reader: No residues given".to_string(),
                    line_number,
                )));
            }
            return Ok((SeqInst::Virtual, is_aa));
        }
        if is_aa && self.gaps.is_empty() && !self.seq_data.is_empty() {
            let seq_data = SeqData::Ncbieaa(std::mem::take(&mut self.seq_data));
            return Ok((SeqInst::Raw { length, seq_data }, is_aa));
        }

        // The pieces between the gaps, at the gaps' wrapped positions; a piece is a `TStr` of
        // `m_SeqData` from its position for its length, which `CTempString` cuts to the residues (empty from a position at or past
        // their end), and the length after the last gap is a `TSeqPos` difference. The residues
        // of a protein are one NCBIeaa literal each (gaps stay a protein's only when their sum
        // wraps to 0).
        //
        // NCBI reference (598d8ae6): c++/include/corelib/tempstr.hpp:377-387
        let size = self.seq_data.len() as u64;
        let span = |pos: u64, len: u64| -> (usize, usize) {
            if pos >= size {
                (0, 0)
            } else {
                (pos as usize, (pos + len.min(size - pos)) as usize)
            }
        };
        let mut pieces = Vec::new();
        let add = |pieces: &mut Vec<DeltaSeq>,
                   residues: &[u8],
                   (start, end): (usize, usize)|
         -> Result<(), ReadError> {
            if is_aa {
                if end > start {
                    pieces.push(DeltaSeq::Literal(Segment {
                        coding: Coding::Ncbieaa,
                        start: start as u32,
                        length: (end - start) as u32,
                    }));
                }
                Ok(())
            } else {
                add_and_split(pieces, residues, start, end, modulus)
            }
        };
        let gap_pos = |gap: &Gap| gap.position as u64 % modulus;
        let n = self.gaps.len();
        if n == 0 || gap_pos(&self.gaps[0]) > 0 {
            let len = if n > 0 {
                gap_pos(&self.gaps[0])
            } else {
                length
            };
            add(&mut pieces, &self.seq_data, span(0, len))?;
        }
        for i in 0..n {
            let gap = self.gaps[i];
            pieces.push(if gap.length == 0 {
                DeltaSeq::LocNull
            } else {
                DeltaSeq::Gap { length: gap.length }
            });
            let next_start = if i == n - 1 {
                current
            } else {
                gap_pos(&self.gaps[i + 1])
            };
            if next_start != gap_pos(&gap) {
                let seq_len = (next_start + modulus - gap_pos(&gap)) % modulus;
                add(&mut pieces, &self.seq_data, span(gap_pos(&gap), seq_len))?;
            }
        }

        if too_long {
            return Ok((SeqInst::TooLong { length }, is_aa));
        }
        if let [piece] = pieces[..] {
            let seq_data = match piece {
                DeltaSeq::Literal(segment) => {
                    debug_assert_eq!(segment.length as usize, self.seq_data.len());
                    match segment.coding {
                        Coding::Ncbi2na => {
                            SeqData::Ncbi2na(pack::iupacna_to_ncbi2na(&self.seq_data))
                        }
                        Coding::Ncbi4na => {
                            SeqData::Ncbi4na(pack::iupacna_to_ncbi4na(&self.seq_data))
                        }
                        Coding::Ncbieaa => {
                            unreachable!("a protein of one literal has no gaps")
                        }
                    }
                }
                DeltaSeq::LocNull | DeltaSeq::Gap { .. } => SeqData::Other,
            };
            return Ok((SeqInst::Raw { length, seq_data }, is_aa));
        }
        let residues = std::mem::take(&mut self.seq_data);
        Ok((
            SeqInst::Delta {
                length,
                residues,
                pieces,
            },
            is_aa,
        ))
    }

    /// `AssembleSeq` without `fParseGaps` (protein input, and every input of the mask readers):
    /// the `>?` gaps become runs of `X` (protein) or `N` (nucleotide), and the record's residues
    /// are the first `length` letters (`GetCurrentPos` with `ePosWithGaps`) of the string that this
    /// builds, of which `keep` are stored (`length`, or 0 for a record that ends before its data
    /// are used). Its letters before a gap are appended where the gap's (wrapped) position is past
    /// the last one, cut to the residues; an append from a position past them throws
    /// `std::out_of_range` (`ReadError::OutOfRange`; only past 2^32 residues, two gaps before
    /// the wrap beyond the residues kept). Only the first `length` letters are stored: NCBI
    /// reserves `length` and grows the string to all of them, here they are not made (a record
    /// whose gaps wrap past 2^32 is shorter than its runs).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1386-1424
    fn substituted(
        &self,
        keep: u64,
        is_aa: bool,
        warn: ReaderWarn<'_>,
    ) -> Result<Vec<u8>, ReadError> {
        let gap_char = if is_aa { b'X' } else { b'N' };
        let residues = &self.seq_data[..];
        let size = residues.len() as u64;
        let mut new_data = Vec::new();
        new_data
            .try_reserve_exact(keep as usize)
            .map_err(|_| ReadError::OutOfMemory)?;
        let room = |new_data: &Vec<u8>| keep as usize - new_data.len();
        let mut pos = 0u64;
        for gap in &self.gaps {
            // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1393-1410
            if gap.mods {
                post_warning(
                    self.flags.ignored_problems,
                    Problem::ModifierFoundButNoneExpected,
                    b"CFastaReader: Gap mods are ignored because gaps are becoming N's or X's in this case.",
                    warn,
                );
            }
            let gap_pos = gap.position as u64 % self.seq_pos_modulus;
            if gap_pos > pos {
                if pos > size {
                    return Err(ReadError::OutOfRange { pos, size });
                }
                let end = gap_pos.min(size) as usize;
                let take = (end - pos as usize).min(room(&new_data));
                new_data.extend_from_slice(&residues[pos as usize..pos as usize + take]);
                pos = gap_pos;
            }
            let run = (gap.length as usize).min(room(&new_data));
            new_data.resize(new_data.len() + run, gap_char);
        }
        if size > pos {
            let take = ((size - pos) as usize).min(room(&new_data));
            new_data.extend_from_slice(&residues[pos as usize..pos as usize + take]);
        }
        Ok(new_data)
    }

    /// `ParseTitle`: the reader's warning about residues at the end of a title, and the
    /// title descriptor.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:670-685
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2018-2044
    /// The long-title and modifier warnings are ignored problems of makeblastdb; the test
    /// for residues at the end of the title is `fasta_reader`'s (`seq_data_in_title_warning`,
    /// one of its two tests runs because exactly one of `fAssumeNuc` and `fAssumeProt` is
    /// set). No caller sets `fAddMods` or `fIgnoreMods`.
    fn parse_title(&mut self, text: &[u8], warn: ReaderWarn<'_>) {
        const WARN_TITLE_LENGTH: usize = 1000;
        let ignored = self.flags.ignored_problems;
        if text.len() > WARN_TITLE_LENGTH {
            post_warning(
                ignored,
                Problem::TooLong,
                format!(
                    "FASTA-Reader: Title is very long: {} characters (max is {WARN_TITLE_LENGTH})",
                    text.len()
                )
                .as_bytes(),
                warn,
            );
        }
        let assumed = match self.flags.molecule {
            Molecule::Nucleotide => Some(false),
            Molecule::Protein => Some(true),
            Molecule::Guessed => None,
        };
        if let Some((message, protein)) = seq_data_in_title_warning(text, assumed) {
            let problem = if protein {
                Problem::UnexpectedAminoAcids
            } else {
                Problem::UnexpectedNucResidues
            };
            post_warning(ignored, problem, message, warn);
        }
        if crate::blastinput::title_parser::has_mods(text) {
            post_warning(
                ignored,
                Problem::ModifierFoundButNoneExpected,
                b"FASTA-Reader: Ignoring FASTA modifier(s) found because the input was not expected to have any.",
                warn,
            );
        }
        // `NStr::TruncateSpacesInPlace` of `processed_title` (`eTrunc_Both`, fasta.cpp:2038).
        let processed_title = truncate_spaces(text);
        if !processed_title.is_empty() {
            self.title_descr = Some(processed_title.to_vec());
        }
    }
}

/// `xSetSeqMol`: the molecule that the first id with one names (`IdentifyAccession`):
/// `Some(true)` for a protein, `Some(false)` for a nucleotide.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:548-563
fn mol_of_ids(ids: &[ParsedSeqId], warn: ReaderWarn<'_>) -> Option<bool> {
    ids.iter().find_map(|id| {
        let info = identify_seq_id(
            id,
            parse_flags::ANY_RAW | parse_flags::ANY_LOCAL,
            &mut |w| warn(ReaderWarning::SeqId(w)),
        );
        if info & F_NUC != 0 {
            Some(false)
        } else if info & F_PROT != 0 {
            Some(true)
        } else {
            None
        }
    })
}

/// `CFormatGuess::SequenceType` of the data and its length with `eST_Default`, for the stored
/// residues (upper case letters, `*` and `-`): `Some(false)` for `eNucleotide`, `Some(true)` for
/// `eProtein`, `None` for `eUndefined`.
///
/// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:250-374
/// The stored residues have no space or digit.
pub(crate) fn sequence_type(data: &[u8]) -> Option<bool> {
    let (mut main_nuc, mut ambig, mut bad_nuc) = (0usize, 0usize, 0usize);
    let (mut amino_acid, mut exotic_aa, mut bad_aa) = (0usize, 0usize, 0usize);
    for &c in data {
        let c = c.to_ascii_uppercase();
        if b"ACGNTU-".contains(&c) {
            main_nuc += 1;
        } else if b"BDHKMRSVWY".contains(&c) {
            ambig += 1;
        } else if !(is_space(c) || c.is_ascii_digit()) {
            bad_nuc += 1;
        }
        if b"ACDEFGHIKLMNPQRSTVWYBZX*-".contains(&c) {
            amino_acid += 1;
        } else if c.is_ascii_alphabetic() {
            exotic_aa += 1;
        } else if !(is_space(c) || c.is_ascii_digit()) {
            bad_aa += 1;
        }
    }
    if bad_nuc + ambig <= main_nuc / 9
        || (bad_nuc + ambig <= main_nuc / 3 && bad_nuc <= (main_nuc + ambig) / 19)
    {
        Some(false)
    } else if bad_aa + exotic_aa <= amino_acid / 9 {
        Some(true)
    } else if bad_nuc == 0 && ambig <= main_nuc / 3 {
        Some(false)
    } else if bad_aa == 0 && exotic_aa <= amino_acid / 9 {
        Some(true)
    } else {
        None
    }
}

/// `CDelta_ext::AddAndSplit` of the IUPACna residues `residues[start..end]`: one literal
/// for each segment of `CSeqConvert::Pack`, or `ReadError::PackOverflow` at a segment whose
/// byte count (`(length + 3) / 4` for NCBI2na, `(length + 1) / 2` for NCBI4na) wraps the
/// `TSeqPos` (2^32, `modulus`) to 0.
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/Delta_ext.cpp:141-147,163-164,183-185,201-205
fn add_and_split(
    pieces: &mut Vec<DeltaSeq>,
    residues: &[u8],
    start: usize,
    end: usize,
    modulus: u64,
) -> Result<(), ReadError> {
    for segment in pack::pack_iupacna(&residues[start..end]) {
        let rounding = match segment.coding {
            Coding::Ncbi2na => 3,
            Coding::Ncbi4na => 1,
            Coding::Ncbieaa => 0,
        };
        if u64::from(segment.length) + rounding >= modulus {
            return Err(ReadError::PackOverflow);
        }
        pieces.push(DeltaSeq::Literal(Segment {
            start: segment.start + start as u32,
            ..segment
        }));
    }
    Ok(())
}

/// The warnings of `ParseGapLine` about the modifiers of a `>?` line.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1190-1340
/// `NStr::Split` with the flags `fSplit_MergeDelimiters` and `fSplit_Truncate` drops empty pieces.
/// The names are looked up by `fasta_reader`'s tables (`gap_type_info`, `linkage_evidence_value`).
///
/// Whether the gap has a specified gap type or linkage evidence, as `AssembleSeq` asks of
/// the `SGap` (quoted at its use): the first known gap type (`unknown` when only evidence is
/// given), and the evidence the switch leaves (cleared to `unspecified` for a type that allows
/// only that, cleared for one that forbids any).
fn gap_mod_warnings(
    mods: &[(Vec<u8>, Vec<u8>)],
    line_number: u64,
    ignored: &[Problem],
    warn: ReaderWarn<'_>,
) -> bool {
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
                    post_warning(ignored, Problem::ParsingModifiers, &message, warn);
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
                        post_warning(ignored, Problem::ParsingModifiers, &message, warn);
                    }
                }
            }
            continue;
        }
        let mut message = b"Unknown gap modifier name(s): ".to_vec();
        message.extend_from_slice(key);
        post_warning(ignored, Problem::UnrecognizedQualifierName, &message, warn);
    }
    if conflicting {
        post_warning(
            ignored,
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
                    ignored,
                    Problem::ExtraModifierFound,
                    b"FASTA-Reader: Unknown gap-type can have linkage-evidence of type 'unspecified' only.",
                    warn,
                );
            }
        }
        LinkEvid::Forbidden => {
            if !evidence.is_empty() {
                post_warning(
                    ignored,
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
                    ignored,
                    Problem::ExpectedModifierMissing,
                    b"CFastaReader: This gap-type should have at least one specified linkage-evidence.",
                    warn,
                );
            }
        }
    }
    let specified_type = gap_type.is_some_and(|value_type| value_type != GapType::Unknown);
    let specified_evidence = matches!(link_evid, LinkEvid::Required)
        && !(evidence.len() == 1 && evidence.contains(&LINKAGE_UNSPECIFIED));
    specified_type || specified_evidence
}

#[cfg(test)]
#[path = "reader_tests.rs"]
mod tests;
