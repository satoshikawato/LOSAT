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
//!   `ParsedSeqId::Local(ObjectId::Id(n))`, `n` counting every record from 1; the first id
//!   is `bsh.GetSeqId()`), `title: Option<Vec<u8>>` (the Seq-descr title for
//!   `CDeflineGenerator`, as the reader stores it) and `letters: Vec<u8>` (upper-case
//!   IUPACna or IUPACaa; a hyphen or a `>?` gap of a nucleotide sequence is `N`, of a
//!   protein `X`). A record has at least one letter: an empty one is the `eNoResidues`
//!   error.
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
use crate::writedb::fasta_source::pack::{iupacna_to_ncbi2na, iupacna_to_ncbi4na, Coding};
use crate::writedb::fasta_source::reader::{
    past_seq_pos_refusal, DeltaSeq, ErrCode, FastaReader, ReadError, ReaderFlags, ReaderWarn,
    Record, SeqInst,
};

// NCBI reference: libstdc++'s `basic_filebuf::xsgetn` (GCC 14, the compiler of NCBI's 2.17.0
// build), which NCBI's line reader calls through `istream::read`:
// ```c++
//           __throw_ios_failure(__N("basic_filebuf::xsgetn "
//                                   "error reading the file"), errno);
// ```
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
    /// What LOSAT refuses (the text before "is not supported by LOSAT's <program>"): a record
    /// past NCBI's `TSeqPos` (`SeqInst::PastSeqPos`).
    Refused(String),
}

/// `CMaskFastaReader`.
pub(crate) struct MaskFastaReader<R: Read> {
    /// `fasta_reader_`; `None` once `GetNextSequence` has returned no record or thrown.
    reader: Option<FastaReader<R>>,
}

impl<R: Read> MaskFastaReader<R> {
    /// `CMaskFastaReader(newInputStream, is_nucl, parse_seqids)` (its flags are quoted at
    /// `ReaderFlags::mask_fasta_reader`). The `CFastaReader` constructor over a stream makes
    /// the line reader of `ILineReader::New`, the `CBufferedLineReader` of makeblastdb:
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:269-270
    /// ```c++
    /// CFastaReader::CFastaReader(CNcbiIstream& in, TFlags flags, FIdCheck f_idcheck)
    ///     : CFastaReader(*(ILineReader::New(in)), flags, f_idcheck) {}
    /// ```
    /// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:69-74
    /// ```c++
    /// CRef<ILineReader> ILineReader::New(CNcbiIstream& is, EOwnership ownership)
    /// {
    ///     CRef<ILineReader> lr;
    ///     lr.Reset(new CBufferedLineReader(is, ownership));
    ///     return lr;
    /// }
    /// ```
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
    /// ```c++
    /// CRef< CSeq_entry > CMaskFastaReader::GetNextSequence()
    /// {
    ///     while( !fasta_reader_.AtEOF() )
    ///     {
    ///         CRef< CSeq_entry > aSeqEntry( null );
    ///         aSeqEntry = fasta_reader_.ReadSet( 1 );
    ///
    ///         if( !input_stream && !input_stream.eof() ) {
    ///             NCBI_THROW( Exception, eBadStream,
    ///                     "error reading input stream" );
    ///         }
    ///
    ///         if( aSeqEntry != 0 && aSeqEntry->IsSeq() &&
    ///             aSeqEntry->GetSeq().IsNa() == is_nucleotide_)
    ///             return aSeqEntry;
    ///         else break;
    ///     }
    ///
    ///     return CRef< CSeq_entry >( null );
    /// }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:442-472
    /// ```c++
    /// CRef<CSeq_entry> CFastaReader::ReadSet(int max_seqs, ILineErrorListener * pMessageListener)
    /// {
    ///     CRef<CSeq_entry> entry(new CSeq_entry);
    ///     if (TestFlag(fOneSeq)) {
    ///         max_seqs = 1;
    ///     }
    ///     for (int i = 0;  i < max_seqs  &&  !GetLineReader().AtEOF();  ++i) {
    ///         try {
    ///             CRef<CSeq_entry> entry2(ReadOneSeq(pMessageListener));
    ///             if (max_seqs == 1) {
    ///                 return entry2;
    ///             }
    ///             if (entry2.NotEmpty())
    ///               entry->SetSet().SetSeq_set().push_back(entry2);
    ///         } catch (const CObjReaderParseException& e) {
    ///             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
    ///                 break;
    ///             } else {
    ///                 throw;
    ///             }
    ///         }
    ///     }
    ///
    ///
    ///     if (entry->IsSet()  &&  entry->GetSet().GetSeq_set().size() == 1) {
    ///         return entry->SetSet().SetSeq_set().front();
    ///     } else {
    ///         entry->Parentize();
    ///         return entry;
    ///     }
    /// }
    /// ```
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
                if let SeqInst::PastSeqPos { length } = record.inst {
                    self.reader = None;
                    return Err(MaskReadError::Refused(past_seq_pos_refusal(length)));
                }
                let parsed_ids = reader.parses_ids();
                Ok(Some(mask_record(record, parsed_ids)))
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
                Err(MaskReadError::Ncbi(exception))
            }
            Err(ReadError::OutOfMemory) => {
                self.reader = None;
                Err(MaskReadError::OutOfMemory)
            }
        }
    }
}

/// The ids, the title and the letters of a record.
///
/// Without `-parse_seqids` the reader generated the record's id (`CSeqIdGenerator`, quoted at
/// `reader::RecordState::parse_def_line`): a local id, counting the records from 1.
fn mask_record(record: Record, parsed_ids: bool) -> MaskRecord {
    let ids = if parsed_ids {
        record.ids
    } else {
        vec![ParsedSeqId::Local(ObjectId::Id(record.local_id as i32))]
    };
    MaskRecord {
        ids,
        title: record.title_descr,
        letters: iupac_letters(record.inst),
    }
}

/// The letters of `CSeqVector(eCoding_Iupac)` over the sequence that `AssembleSeq` made: a
/// raw sequence of one literal's data, or a delta of literals (the reader without
/// `fParseGaps` keeps no gaps: `>?` gaps became `N` or `X`, and a hyphen is a residue). Each
/// literal of a nucleotide sequence is decoded from the NCBI2na or NCBI4na data that
/// `AddAndSplit` packs it into (`pack`), so a `U` comes out as `T` and a hyphen (NCBI4na 15)
/// as `N`; a protein is NCBIeaa, whose `*` and `-` are `X` in IUPACaa.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector.cpp:1281-1284
/// ```c++
/// void CSeqVector::SetIupacCoding(void)
/// {
///     SetCoding(IsProtein()? CSeq_data::e_Iupacaa: CSeq_data::e_Iupacna);
/// }
/// ```
/// A virtual sequence (`fDisableNoResidues`) and a record above makeblastdb's length limit do
/// not come from this reader's flags, and a record past NCBI's `TSeqPos` is refused before
/// (`MaskFastaReader::next`).
fn iupac_letters(inst: SeqInst) -> Vec<u8> {
    match inst {
        SeqInst::Raw { length, seq_data } => match seq_data {
            SeqData::Ncbieaa(residues) => residues.iter().map(|&c| ncbieaa_to_iupacaa(c)).collect(),
            SeqData::Ncbi2na(data) => ncbi2na_letters(&data, length as usize),
            SeqData::Ncbi4na(data) => ncbi4na_letters(&data, length as usize),
            SeqData::Other => Vec::new(),
        },
        SeqInst::Delta {
            length,
            residues,
            pieces,
        } => {
            let mut letters = Vec::with_capacity(length as usize);
            for piece in pieces {
                if let DeltaSeq::Literal(segment) = piece {
                    let start = segment.start as usize;
                    let literal = &residues[start..start + segment.length as usize];
                    letters.extend(match segment.coding {
                        Coding::Ncbi2na => {
                            ncbi2na_letters(&iupacna_to_ncbi2na(literal), literal.len())
                        }
                        Coding::Ncbi4na => {
                            ncbi4na_letters(&iupacna_to_ncbi4na(literal), literal.len())
                        }
                    });
                }
            }
            letters
        }
        SeqInst::Virtual | SeqInst::TooLong { .. } | SeqInst::PastSeqPos { .. } => Vec::new(),
    }
}

/// The IUPACna letters of `length` residues of NCBI2na data, four to a byte from its high
/// bits (the `ncbi2na` to `iupacna` map of NCBI's code tables).
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/seqport_util.cpp:6468-6477
/// ```c++
///     " {\n",
///     " from ncbi2na ,\n",
///     " to iupacna ,\n",
///     " num 4 ,\n",
///     " table {\n",
///     " 65,     -- A\n",
///     " 67,     -- C\n",
///     " 71,     -- G\n",
///     " 84 }    -- T\n",
///     " } ,\n",
/// ```
fn ncbi2na_letters(data: &[u8], length: usize) -> Vec<u8> {
    const LETTERS: [u8; 4] = *b"ACGT";
    (0..length)
        .map(|i| LETTERS[usize::from(data[i / 4] >> (6 - 2 * (i % 4)) & 3)])
        .collect()
}

/// The IUPACna letters of `length` residues of NCBI4na data, two to a byte from its high
/// half (the `ncbi4na` to `iupacna` map of NCBI's code tables: a gap, code 0, is `N`).
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/seqport_util.cpp:6489-6508
/// ```c++
///     " from ncbi4na ,\n",
///     " to iupacna ,\n",
///     " num 16 ,\n",
///     " table {\n",
///     " 78,    -- gap -> N\n",
///     " 65,    -- A\n",
///     " 67,    -- C\n",
///     " 77,    -- M\n",
///     " 71,    -- G\n",
///     " 82,    -- R\n",
///     " 83,    -- S\n",
///     " 86,    -- V\n",
///     " 84,    -- T\n",
///     " 87,    -- W\n",
///     " 89,    -- Y\n",
///     " 72,    -- H\n",
///     " 75,    -- K\n",
///     " 68,    -- D\n",
///     " 66,    -- B\n",
///     " 78 }   -- N\n",
/// ```
fn ncbi4na_letters(data: &[u8], length: usize) -> Vec<u8> {
    const LETTERS: [u8; 16] = *b"NACMGRSVTWYHKDBN";
    (0..length)
        .map(|i| LETTERS[usize::from(data[i / 2] >> (4 - 4 * (i % 2)) & 0xf)])
        .collect()
}

/// The IUPACaa letter of an NCBIeaa residue of the reader (a letter, `*` or `-`): the
/// `ncbieaa` to `iupacaa` map of NCBI's code tables, from `*` (42) on.
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/seqport_util.cpp:6566-6575,6594-6595
/// ```c++
///     " from ncbieaa ,\n",
///     " to iupacaa ,\n",
///     " num 49 ,\n",
///     " start-at 42 ,\n",
///     " table {\n",
///     " 88 ,   -- termination -> X\n",
///     " 255,\n",
///     " 255,\n",
///     " 88,    -- Gap -> X\n",
///     " 255,\n",
/// ...
///     " 65 ,    -- from here they map directly\n",
///     " 66 ,\n",
/// ```
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
