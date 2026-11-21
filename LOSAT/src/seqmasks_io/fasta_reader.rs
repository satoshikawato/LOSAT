//! The FASTA input of dustmasker, segmasker and windowmasker: `CMaskFastaReader`
//! (`c++/src/objtools/seqmasks_io/mask_fasta_reader.cpp`) over NCBI's `CFastaReader`, whose
//! one port is `writedb::fasta_source::reader` (here with
//! [`ReaderFlags::mask_fasta_reader`]), and what a masking program takes from each `Bioseq`:
//! its Seq-ids, its title descriptor and its letters (`CSeqVector`, `eCoding_Iupac`).
//!
//! # API
//!
//! - [`MaskFastaReader::new`]`(input: R, protein: bool, parse_seqids: bool) ->
//!   Result<MaskFastaReader<R>, std::io::Error>`: `CMaskFastaReader(stream, is_nucl =
//!   !protein, parse_seqids)`. NCBI's line reader reads its first buffer here; a failed read
//!   is the error (NCBI throws the stream's `std::ios_base::failure` from the constructor,
//!   after the program has opened its output).
//! - [`MaskFastaReader::next`]`(&mut self, warn: ReaderWarn<'_>) ->
//!   Result<Option<MaskRecord>, MaskReadError>`: `GetNextSequence`. `Ok(None)` at the end of
//!   the input, also for input of blank and comment lines only (`eEOF`, which `ReadSet`
//!   catches). After `None` or an error the reader returns `Ok(None)`.
//! - [`MaskRecord`]: `ids: Vec<ParsedSeqId>` (`bsh.GetId()`; without `-parse_seqids` one
//!   `ParsedSeqId::Local(ObjectId::Id(n))`, `n` counting every record from 1; the writers name
//!   the sequence by `writer::access_seq_id` of them, `bsh.GetSeqId()`), `title:
//!   Option<Vec<u8>>` (the Seq-descr title for `CDeflineGenerator`, as the reader stores it)
//!   and `letters: Vec<u8>` (upper-case IUPACna or IUPACaa; a hyphen or a `>?` gap of a
//!   nucleotide sequence is `N`, of a protein `X`). A record has at least one letter: an empty
//!   one is the `eNoResidues` error.
//! - [`MaskReadError`]: `Ncbi(NcbiException)`, a `CObjReaderParseException` (dustmasker
//!   writes `e.report_default_ex(106, 16, "Application's execution failed")` and exits with
//!   3; segmasker writes `e.what()` and a line feed and exits with 1); `Read(std::io::Error)`,
//!   a failed read of the input; `OutOfMemory`, the `std::bad_alloc` of a gap line's length;
//!   `Refused(what)`, a record past NCBI's 32-bit length (`SeqInst::PastSeqPos`), which LOSAT
//!   refuses.
//! - Warnings (`writedb::fasta_source::reader::ReaderWarning`, passed to a
//!   `FnMut(ReaderWarning<'_>)`): `Line(message)`, which NCBI writes as the line `message`
//!   (`LOG_POST`), and `SeqId(warning)`, which it writes as `Warning: (1306.<subcode>)
//!   <message>` with the application's prefix. Both are of severity Warning: dustmasker and
//!   segmasker (NCBI's release post level, Error) write none, windowmasker
//!   (`SetDiagPostLevel(eDiag_Warning)`) writes them. The texts are those of dustmasker
//!   2.17.0 with `DIAG_POST_LEVEL=Warning` (`tests/dustmasker/reader`, `*_warn*` cases).

use std::io::Read;

use crate::ncbi_exception::NcbiException;
use crate::seqid::model::{ObjectId, ParsedSeqId};
use crate::writedb::bioseq::SeqData;
use crate::writedb::fasta_source::pack::{self, Coding};
use crate::writedb::fasta_source::reader::{
    past_seq_pos_refusal, DeltaSeq, ErrCode, FastaReader, ReadError, ReaderFlags, ReaderWarn,
    Record, SeqInst,
};

// libstdc++ (GCC 14, the compiler of NCBI's 2.17.0 build): bits/fstream.tcc
// (`basic_filebuf::xsgetn`), which NCBI's line reader calls through `istream::read`.
// (`what()` is the text, `: ` and the `strerror` of `errno`: `-in` naming a directory gives
// `basic_filebuf::xsgetn error reading the file: Is a directory`, dustmasker and segmasker
// fixtures `in_dir`.)
/// The `what()` of the `std::ios_base::failure` that a failed read of the input
/// (`MaskReadError::Read`) is in NCBI's stream; `None` for an error without an OS error code,
/// whose text LOSAT does not know.
pub(crate) fn ios_failure_what(error: &std::io::Error) -> Option<String> {
    let code = error.raw_os_error()?;
    let text = error.to_string();
    let strerror = text
        .strip_suffix(&format!(" (os error {code})"))
        .unwrap_or(&text);
    Some(format!(
        "basic_filebuf::xsgetn error reading the file: {strerror}"
    ))
}

/// What a masking program takes from one `Bioseq` of the reader.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct MaskRecord {
    /// `CBioseq::GetId()`: the ids of the definition line (`-parse_seqids`), or the id that
    /// the reader generated.
    pub(crate) ids: Vec<ParsedSeqId>,
    /// The title descriptor (`CSeqdesc::e_Title`): the title of the definition line without
    /// white space at its ends, `None` for none.
    pub(crate) title: Option<Vec<u8>>,
    /// `bsh.GetSeqVector(CBioseq_Handle::eCoding_Iupac)`.
    pub(crate) letters: Vec<u8>,
}

/// Why `GetNextSequence` ended without a record.
#[derive(Debug)]
pub(crate) enum MaskReadError {
    /// The reader's `CObjReaderParseException`.
    Ncbi(NcbiException),
    /// A failed read of the input (`std::ios_base::failure` of the stream in NCBI).
    Read(std::io::Error),
    /// The `std::bad_alloc` of a record whose `>?` lines state more residues than can be
    /// allocated.
    OutOfMemory,
    /// What LOSAT refuses (the text before `is not supported by LOSAT's <program>`): a record
    /// past NCBI's `TSeqPos` (`SeqInst::PastSeqPos`).
    Refused(String),
}

/// `CMaskFastaReader`.
pub(crate) struct MaskFastaReader<R: Read> {
    /// `fasta_reader_`; `None` once `GetNextSequence` has returned no record or thrown.
    reader: Option<FastaReader<R>>,
}

impl<R: Read> MaskFastaReader<R> {
    /// `CMaskFastaReader(newInputStream, is_nucl, parse_seqids)` (its flags are cited at
    /// `ReaderFlags::mask_fasta_reader`). The `CFastaReader` constructor over a stream makes
    /// the line reader of `ILineReader::New`, the `CBufferedLineReader` of makeblastdb:
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:269-270
    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:69-74
    /// `CBufferedLineReader` reads its first buffer in its constructor. The stream's state
    /// test after it (`eBadStream`) does not throw for a stream that NCBI's argument
    /// parser opened: a failed read throws first (dustmasker 2.17.0 with `-in` a directory:
    /// `Error: (106.7) Application's execution failed: basic_filebuf::xsgetn error reading
    /// the file: Is a directory`, fixture `in_dir`).
    pub(crate) fn new(input: R, protein: bool, parse_seqids: bool) -> Result<Self, std::io::Error> {
        let mut reader =
            FastaReader::new(input, ReaderFlags::mask_fasta_reader(protein, parse_seqids));
        match reader.take_read_error() {
            Some(error) => Err(error),
            None => Ok(Self {
                reader: Some(reader),
            }),
        }
    }

    /// `GetNextSequence`: the next record, or `None` at the end of the input.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_fasta_reader.cpp:43-62
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:442-472
    /// After `eEOF` `ReadSet` returns an entry that is neither a Seq nor a set, which ends
    /// the reading. With `fForceType` every record is of the reader's molecule, so the
    /// molecule test always passes. A failed read throws from within the reader in NCBI;
    /// LOSAT's line reader keeps the error and reads as at the end of the input, and the
    /// error is returned after the record that was being read (that record is lost in
    /// both).
    pub(crate) fn next(
        &mut self,
        warn: ReaderWarn<'_>,
    ) -> Result<Option<MaskRecord>, MaskReadError> {
        let Some(reader) = &mut self.reader else {
            return Ok(None);
        };
        if reader.at_eof() {
            let read_error = reader.take_read_error();
            self.reader = None;
            return match read_error {
                Some(error) => Err(MaskReadError::Read(error)),
                None => Ok(None),
            };
        }
        let result = reader.read_one_seq(warn);
        if let Some(error) = reader.take_read_error() {
            self.reader = None;
            return Err(MaskReadError::Read(error));
        }
        match result {
            Ok(record) => {
                if let SeqInst::PastSeqPos { length, .. } = record.inst {
                    let record = mask_record(
                        Record {
                            inst: SeqInst::Virtual,
                            ..record
                        },
                        reader,
                    );
                    self.reader = None;
                    let mut id = Vec::new();
                    if let Some(first) = record.ids.first() {
                        crate::seqid::label::fasta(first, &mut id);
                    }
                    return Err(MaskReadError::Refused(past_seq_pos_refusal(
                        &String::from_utf8_lossy(&id),
                        length,
                    )));
                }
                Ok(Some(mask_record(record, reader)))
            }
            Err(ReadError::Parse(error)) if error.code == ErrCode::Eof => {
                self.reader = None;
                Ok(None)
            }
            Err(ReadError::Parse(error)) => {
                self.reader = None;
                Err(MaskReadError::Ncbi(error.exception()))
            }
            Err(ReadError::Ncbi(exception)) => {
                self.reader = None;
                Err(MaskReadError::Ncbi(*exception))
            }
            Err(ReadError::OutOfMemory) => {
                self.reader = None;
                Err(MaskReadError::OutOfMemory)
            }
        }
    }
}

/// The ids, the title and the letters of a record that `reader` read.
///
/// Without `-parse_seqids` the reader generated the record's id (`CSeqIdGenerator`, cited at
/// `reader::RecordState::parse_def_line`): a local id, counting the records from 1.
fn mask_record<R: Read>(record: Record, reader: &mut FastaReader<R>) -> MaskRecord {
    let ids = if reader.parses_ids() {
        record.ids
    } else {
        vec![ParsedSeqId::Local(ObjectId::Id(record.local_id as i32))]
    };
    MaskRecord {
        ids,
        title: record.title_descr,
        letters: iupac_letters(record.inst, reader),
    }
}

/// The letters of `CSeqVector(eCoding_Iupac)` over the sequence that `AssembleSeq` made: a
/// raw sequence of one literal's data, or a delta of literals (the reader without
/// `fParseGaps` keeps no gaps: `>?` gaps became `N` or `X`, and a hyphen is a residue). Each
/// literal of a nucleotide sequence reads back from the NCBI2na or NCBI4na data that
/// `AddAndSplit` packs it into (`pack`), so a `U` comes out as `T` and a hyphen (NCBI4na 15)
/// as `N`; a protein is NCBIeaa, whose `*` and `-` are `X` in IUPACaa.
///
/// The letters are made in the buffer of the residues the reader stored (a delta's
/// `residues`, a protein's data, or for a raw nucleotide sequence the residues it packed,
/// `FastaReader::take_raw_residues`), each residue replaced by the letter it reads back as
/// ([`IUPACNA_FROM`]): the sequence is held once.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector.cpp:1281-1284
///
/// A virtual sequence (`fDisableNoResidues`) and a record above makeblastdb's length limit do
/// not come from this reader's flags, and a record past NCBI's `TSeqPos` is refused before
/// (`MaskFastaReader::next`).
pub(crate) fn iupac_letters<R: Read>(inst: SeqInst, reader: &mut FastaReader<R>) -> Vec<u8> {
    match inst {
        SeqInst::Raw {
            seq_data: SeqData::Ncbieaa(mut residues),
            ..
        } => {
            iupacaa_letters(&mut residues);
            residues
        }
        SeqInst::Raw {
            length,
            seq_data: SeqData::Ncbi2na(packed),
        } => {
            drop(packed);
            raw_letters(reader, length, Coding::Ncbi2na)
        }
        SeqInst::Raw {
            length,
            seq_data: SeqData::Ncbi4na(packed),
        } => {
            drop(packed);
            raw_letters(reader, length, Coding::Ncbi4na)
        }
        SeqInst::Raw {
            seq_data: SeqData::Other,
            ..
        } => Vec::new(),
        SeqInst::Delta {
            residues, pieces, ..
        } => {
            // The literals are in order and do not overlap, so each one's letters move to
            // the end of the letters before it (they are already there when the literals
            // cover the residues, as without `fParseGaps`).
            let mut letters = residues;
            let mut end = 0;
            for piece in pieces {
                if let DeltaSeq::Literal(segment) = piece {
                    let start = segment.start as usize;
                    let length = segment.length as usize;
                    debug_assert!(start >= end);
                    if start != end {
                        letters.copy_within(start..start + length, end);
                    }
                    read_back(&mut letters[end..end + length], segment.coding);
                    end += length;
                }
            }
            letters.truncate(end);
            letters
        }
        SeqInst::Virtual | SeqInst::TooLong { .. } | SeqInst::PastSeqPos { .. } => Vec::new(),
    }
}

/// The letters of a raw nucleotide sequence of `length` residues packed in `coding`: the
/// residues that the reader packed, read back in place.
fn raw_letters<R: Read>(reader: &mut FastaReader<R>, length: u64, coding: Coding) -> Vec<u8> {
    let mut letters = reader.take_raw_residues();
    assert_eq!(
        letters.len() as u64,
        length,
        "the reader's residues are the record's"
    );
    read_back(&mut letters, coding);
    letters
}

/// Replaces each IUPACna residue of `residues` by the letter it reads back as from `coding`.
fn read_back(residues: &mut [u8], coding: Coding) {
    let table = match coding {
        Coding::Ncbi2na => &IUPACNA_FROM[0],
        Coding::Ncbi4na => &IUPACNA_FROM[1],
    };
    for residue in residues {
        *residue = table[usize::from(*residue)];
    }
}

/// For each byte, the IUPACna letter that a residue reads back as after `pack` stores it in
/// NCBI2na (row 0) or NCBI4na (row 1): the letter of its code in NCBI's code tables. The
/// code of one residue fills its own bits of a packed byte (at most 3 in NCBI2na, at most
/// 15 in NCBI4na), so the letters of a literal are this map of its residues.
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/seqport_util.cpp:6468-6477 (NCBI2na to
/// IUPACna: `A`, `C`, `G`, `T`), :6489-6508 (NCBI4na to IUPACna: a gap, code 0, is `N`)
const IUPACNA_FROM: [[u8; 256]; 2] = {
    const NCBI2NA: [u8; 4] = *b"ACGT";
    const NCBI4NA: [u8; 16] = *b"NACMGRSVTWYHKDBN";
    let mut table = [[0; 256]; 2];
    let mut byte = 0;
    while byte < 256 {
        table[0][byte] = NCBI2NA[pack::ncbi2na_code(byte as u8) as usize & 3];
        table[1][byte] = NCBI4NA[pack::ncbi4na_code(byte as u8) as usize & 15];
        byte += 1;
    }
    table
};

/// Replaces each NCBIeaa residue of a protein by its IUPACaa letter. Out of line: inlined into
/// the reader, LTO vectorized the loop at half the width of its own function (segmasker on
/// proteins, 0.05 % more instructions than the copy it replaces).
#[inline(never)]
fn iupacaa_letters(residues: &mut [u8]) {
    for residue in residues {
        *residue = ncbieaa_to_iupacaa(*residue);
    }
}

/// The IUPACaa letter of an NCBIeaa residue of the reader (a letter, `*` or `-`): the
/// `ncbieaa` to `iupacaa` map of NCBI's code tables, from `*` (42) on.
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/seqport_util.cpp:6566-6575,6594-6595
/// The letters `A` to `Z` map to themselves (`J`, `O` and `U` too); the reader stores no
/// other byte.
fn ncbieaa_to_iupacaa(residue: u8) -> u8 {
    match residue {
        b'*' | b'-' => b'X',
        _ => residue,
    }
}

#[cfg(test)]
mod tests;
