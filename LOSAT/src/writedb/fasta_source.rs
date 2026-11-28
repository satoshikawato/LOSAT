//! makeblastdb's reading of a FASTA file (`CFastaBioseqSource` of
//! `c++/src/objtools/blast/seqdb_writer/build_db.cpp` over `CFastaReader`), without
//! `-parse_seqids`.
//!
//! A record becomes the `Bioseq` that `CBuildDatabase::AddSequences` adds: read by
//! `CFastaReader::ReadOneSeq` with makeblastdb's flags (`reader`, over NCBI's
//! `CBufferedLineReader`, `line_reader`), with the title that
//! `CWriteDB::ExtractBioseqDeflines` gives its Blast-def-line (`title`) and the sequence
//! data that `s_FixBioseqDeltas` leaves (`fix_bioseq_deltas`; the reader packs nucleotide
//! residues with `CSeqConvert::Pack`, `pack`). The nucleotide reader keeps the residues as
//! NCBI4na codes (`reader::Ncbi4naResidues`), from which the packing and the delta's data are
//! made; NCBI's keeps their letters and converts them.
//!
//! What `AddSequences` does with the sequence afterwards is the caller's: it skips a
//! sequence of the other molecule, stops at a length above 2,147,483,647, and ignores a
//! sequence of length 0 or without sequence data.

pub(crate) mod defline_ids;
mod line_reader;
pub(crate) mod pack;
pub(crate) mod reader;
#[cfg(test)]
mod tests;
mod title;

pub(crate) use defline_ids::post_seq_id_warning;
pub(crate) use title::fasta_reader_fasta;

use std::io::Read;

use super::bioseq::{Bioseq, DeltaFix, ParsedIds, SeqData};
use pack::Ncbi4naBuilder;
use reader::{
    DeltaSeq, FastaReader, Ncbi4naResidues, ReadError, ReaderFlags, ReaderWarning, ResidueStore,
    SeqInst,
};

/// An error that ends the reading: its message and how makeblastdb reports it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceError {
    pub message: Vec<u8>,
    pub kind: SourceErrorKind,
}

/// How makeblastdb reports a `SourceError`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourceErrorKind {
    /// NCBI's error: the message after "BLAST Database creation error: ".
    Creation,
    /// A failed allocation for a length that the input states (`std::bad_alloc` in NCBI):
    /// makeblastdb prints the message alone, "BLAST ran out of memory", and exits with its
    /// status for that.
    OutOfMemory,
    /// Another standard exception (`message` is its `what()`), which makeblastdb's general
    /// handler prints after "Error: " (`WriteDbError::Other`): the `std::out_of_range` of a
    /// record past NCBI's 32-bit positions (`reader::ReadError::OutOfRange`).
    Std,
    /// LOSAT's refusal (`message` ends with "is not supported by LOSAT's makeblastdb";
    /// `WriteDbError::Unsupported`): a record on which NCBI's reader crashes
    /// (`reader::ReadError::PackOverflow`).
    Refused,
}

impl SourceError {
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:256-259
    // `CFastaBioseqSource::GetNext` catches the reader's parse exception only and
    // `s_FixBioseqDeltas` only `CInvalidChoiceSelection`, so a `std::bad_alloc` of
    // either reaches this clause (after `AddFasta` has ended the build).
    fn out_of_memory() -> Self {
        Self {
            message: crate::wording::OUT_OF_MEMORY.as_bytes().to_vec(),
            kind: SourceErrorKind::OutOfMemory,
        }
    }
}

/// A warning line of the reader (without the final newline).
pub type Warn<'a> = &'a mut dyn FnMut(&[u8]);

/// `CFastaBioseqSource`: the sequences of one FASTA input, one at a time.
pub struct FastaBioseqSource<R: Read> {
    /// `m_FastaReader` with its `m_LineReader`; `None` once `GetNext` has reset the line
    /// reader (the end of the input) or thrown.
    reader: Option<SourceReader<R>>,
}

/// The reader of a source, with the store of residues that its sequences' data is made from.
enum SourceReader<R: Read> {
    /// Nucleotide input: NCBI4na codes (`reader::Ncbi4naResidues`).
    Nucleotide(FastaReader<R, Ncbi4naResidues>),
    /// Protein input: NCBIeaa letters, which are the sequence data.
    Protein(FastaReader<R>),
    /// Nucleotide input read into letters and fixed from them
    /// (`tests::fix_bioseq_deltas_of_letters`): the tests' reference for `Nucleotide`.
    #[cfg(test)]
    NucleotideLetters(FastaReader<R>),
}

/// The length (`Seq-inst.length`), the sequence data and what the fix made of a delta
/// (`fix_bioseq_deltas`).
type Fixed = (u64, Option<SeqData>, Option<DeltaFix>);

impl<R: Read> FastaBioseqSource<R> {
    /// `CFastaBioseqSource` of the file, `is_protein` and `parse_ids`, with `long_ids` false: with
    /// `parse_ids` (`-parse_seqids`) the reader parses the ids of each definition line
    /// (flags `fRequireID` and `fParseRawID`) and the Bioseq carries them (`Bioseq::parsed`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:701-734
    /// The flags are cited at `reader::ReaderFlags::fasta_bioseq_source`.
    pub fn with_parse_ids(input: R, is_protein: bool, parse_ids: bool) -> Self {
        let reader = if is_protein {
            SourceReader::Protein(FastaReader::new(
                input,
                ReaderFlags::fasta_bioseq_source(true, parse_ids),
            ))
        } else {
            SourceReader::Nucleotide(FastaReader::nucleotide_db(input, parse_ids))
        };
        Self {
            reader: Some(reader),
        }
    }

    /// The nucleotide source of `with_parse_ids` whose reader keeps the residues as letters,
    /// with the fix from letters (`tests::fix_bioseq_deltas_of_letters`): the tests' reference.
    #[cfg(test)]
    pub(crate) fn nucleotide_letters(input: R, parse_ids: bool) -> Self {
        Self {
            reader: Some(SourceReader::NucleotideLetters(FastaReader::new(
                input,
                ReaderFlags::fasta_bioseq_source(false, parse_ids),
            ))),
        }
    }

    /// The source with its reader's positions counted modulo `modulus`
    /// (`FastaReader::with_seq_pos_modulus`).
    #[cfg(test)]
    pub(crate) fn with_seq_pos_modulus(mut self, modulus: u64) -> Self {
        self.reader = self.reader.map(|reader| match reader {
            SourceReader::Nucleotide(reader) => {
                SourceReader::Nucleotide(reader.with_seq_pos_modulus(modulus))
            }
            SourceReader::Protein(reader) => {
                SourceReader::Protein(reader.with_seq_pos_modulus(modulus))
            }
            SourceReader::NucleotideLetters(reader) => {
                SourceReader::NucleotideLetters(reader.with_seq_pos_modulus(modulus))
            }
        });
        self
    }

    /// `CFastaBioseqSource::GetNext`, then the title of
    /// `CWriteDB::ExtractBioseqDeflines` and `s_FixBioseqDeltas`: the next
    /// sequence, or `None` at the end of the input.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:755-792
    /// `x_EditAndAddBioseq` then takes the deflines and fixes the deltas:
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:473-488
    /// `AddSequences` tests the length of the sequence (`Seq-inst.length`, `Bioseq::length`)
    /// before that call, and `s_FixBioseqDeltas` gives the raw sequence it makes the length
    /// of its data (`Bioseq::delta_fix`). A record above 2,147,483,647 residues is returned
    /// with its length and no sequence data: NCBI stops at that length, before
    /// `x_EditAndAddBioseq`. NCBI counts in 32 bits: a record of 4,294,967,296 residues or
    /// more (gaps included) has the length of its count modulo 2^32, and the pieces of its
    /// delta those of its wrapped positions (`reader::RecordState::assemble_seq`).
    ///
    /// A failed read of the input ends the reading with "Read error" and the system's
    /// message (LOSAT's text: NCBI's line reader throws a `CIOException` "Read error").
    /// After an error the source returns `None`.
    pub fn next_bioseq(&mut self, warn: Warn<'_>) -> Result<Option<Bioseq>, SourceError> {
        let next = match &mut self.reader {
            None => return Ok(None),
            Some(SourceReader::Nucleotide(reader)) => next_of(reader, warn, fix_bioseq_deltas),
            Some(SourceReader::Protein(reader)) => {
                next_of(reader, warn, |inst| Ok(fix_protein(inst)))
            }
            #[cfg(test)]
            Some(SourceReader::NucleotideLetters(reader)) => {
                next_of(reader, warn, tests::fix_bioseq_deltas_of_letters)
            }
        };
        // The reader ends with the input and at an error.
        if !matches!(next, Ok(Some(_))) {
            self.reader = None;
        }
        next
    }
}

/// `next_bioseq` of `reader`, whose deltas `fix` fixes (`s_FixBioseqDeltas`): the next sequence,
/// `None` at the end of the input, or the error that ends the reading.
fn next_of<R: Read, S: ResidueStore>(
    reader: &mut FastaReader<R, S>,
    warn: Warn<'_>,
    fix: impl FnOnce(SeqInst<S::Delta>) -> Result<Fixed, std::collections::TryReserveError>,
) -> Result<Option<Bioseq>, SourceError> {
    if reader.at_eof() {
        return match reader.take_read_error() {
            Some(error) => Err(read_error_message(&error)),
            None => Ok(None),
        };
    }
    // makeblastdb writes a warning of the reader as its line, and a warning of the
    // Seq-id parser with its prefix.
    let record = reader.read_one_seq(&mut |warning| match warning {
        ReaderWarning::Line(line) => warn(line),
        ReaderWarning::SeqId(warning) => post_seq_id_warning(warning, &mut *warn),
    });
    if let Some(error) = reader.take_read_error() {
        return Err(read_error_message(&error));
    }
    let record = match record {
        Ok(record) => record,
        Err(ReadError::OutOfMemory) => return Err(SourceError::out_of_memory()),
        Err(ReadError::OutOfRange { pos, size }) => {
            return Err(SourceError {
                message: reader::out_of_range_what(pos, size).into_bytes(),
                kind: SourceErrorKind::Std,
            });
        }
        Err(ReadError::PackOverflow) => {
            return Err(SourceError {
                message: format!(
                    "{} is not supported by LOSAT's makeblastdb",
                    reader::pack_overflow_what()
                )
                .into_bytes(),
                kind: SourceErrorKind::Refused,
            });
        }
        Err(ReadError::Ncbi(exception)) => {
            // Only convert2blastmask's reader throws other exceptions.
            return Err(SourceError {
                message: exception.get_msg().to_vec(),
                kind: SourceErrorKind::Creation,
            });
        }
        Err(ReadError::Parse(error)) => {
            let message = error.message;
            let message = message
                .strip_prefix("CFastaReader: ")
                .map(str::to_string)
                .unwrap_or(message);
            return Err(SourceError {
                message: message.into_bytes(),
                kind: SourceErrorKind::Creation,
            });
        }
    };
    let is_aa = record.is_aa;
    // With `-parse_seqids` the definition lines depend on the ids that
    // `AddSequences` edits first: the writer builds them (`writedb::deflines`).
    let (title, parsed) = if reader.parses_ids() {
        let parsed = ParsedIds {
            ids: record.ids,
            hist_old_id: record.hist_old_id,
            defline_field: record.defline_field,
            title_descr: record.title_descr,
        };
        (Vec::new(), Some(parsed))
    } else {
        match extract_title(&record.defline_field, record.title_descr) {
            Ok(title) => (title, None),
            Err(message) => {
                return Err(SourceError {
                    message: message.to_string().into_bytes(),
                    kind: SourceErrorKind::Creation,
                });
            }
        }
    };
    let Ok((length, seq_data, delta_fix)) = fix(record.inst) else {
        return Err(SourceError::out_of_memory());
    };
    Ok(Some(Bioseq {
        local_id: record.local_id.to_string(),
        title,
        is_aa,
        length,
        seq_data,
        parsed,
        delta_fix,
    }))
}

fn read_error_message(error: &std::io::Error) -> SourceError {
    SourceError {
        message: format!("Read error: {error}").into_bytes(),
        kind: SourceErrorKind::Creation,
    }
}

/// The title of the Blast-def-line that `CWriteDB::ExtractBioseqDeflines` of the bioseq `bs`
/// (`parse_ids` false, `long_ids` false, `scan` true) builds for a record of `CFastaReader`: the
/// one of the reader's user object (`title::fasta_reader_title`), or, when that object has an empty
/// definition line, the title descriptor (`x_BuildDeflinesFromBioseq`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:802-829
/// `s_UseFastaReaderDeflines` keeps the deflines of the user object, because the record
/// has only a local ID:
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:736-750
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:610-613,694
/// The definition line is empty only for the empty line of `line_reader`'s module
/// comment read as a definition line; such a record has no residues.
fn extract_title(
    defline_field: &[u8],
    title_descr: Option<Vec<u8>>,
) -> Result<Vec<u8>, &'static str> {
    Ok(match title::fasta_reader_title(defline_field)? {
        Some(title) => title,
        None => title_descr.unwrap_or_default(),
    })
}

/// `s_FixBioseqDeltas` of a protein: a sequence with `seq-data` (raw) or without a delta
/// (virtual) is returned as it is, and a delta (gaps that wrapped to 0 past 2^32) is the
/// caller's (`DeltaFix::ProteinDelta`: NCBI throws for it).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:290-295,306-311
fn fix_protein(inst: SeqInst) -> Fixed {
    match inst {
        SeqInst::Virtual => (0, None, None),
        SeqInst::TooLong { length } => (length, None, None),
        SeqInst::Raw { length, seq_data } => (length, Some(seq_data), None),
        SeqInst::Delta { length, .. } => (length, None, Some(DeltaFix::ProteinDelta)),
    }
}

/// `s_FixBioseqDeltas` of a nucleotide sequence: the length (`Seq-inst.length`), the sequence
/// data and what the fix made of a delta of a record, whose residues are NCBI4na codes
/// (`reader::Ncbi4naResidues`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:290-295,306-311,318-334,343-382,395-398
/// A sequence with `seq-data` (raw) or without a delta (virtual) is returned as it is.
/// A literal without a length or a `Seq-loc` other than null does not occur. A delta whose
/// first piece is not a literal (a gap of unknown length) is returned as it is, without
/// `seq-data`. The NCBI8na code of a literal's residue is the NCBI4na code of its letter
/// (`C2naTo8na` of an NCBI2na literal's base, `C4naTo8na` of an NCBI4na one; cited at
/// `pack::Ncbi4naBuilder`), which is what the reader stored.
///
/// A delta of literals only (no gap: `pack`'s segments of the residues, in order) is the
/// residues: its data are the stored codes, taken as they are. NCBI appends each literal to
/// `seq8na` and converts it, the same bytes.
///
/// A delta with gap pieces: NCBI appends every piece to `seq8na`, a gap as that many NCBI8na 15
/// and a gap of unknown length as one, and takes its size, an `int` (`length`), as the
/// `TSeqPos` of the letters it converts and of `Seq-inst.length`: the sum of the pieces modulo
/// 2^32 (an `int` that is negative there is the same `TSeqPos`). Past 2^32 letters (gaps that
/// wrapped NCBI's `TSeqPos`) that is fewer than the pieces hold, and it can differ from
/// `Seq-inst.length` (one letter for each gap of unknown length, and the pieces of a record
/// whose wrapped gap positions repeat or skip residues). Only those letters are made here
/// (NCBI makes all of the pieces' letters first), each literal's copied from the stored codes
/// (`Ncbi4naBuilder::push_ncbi4na`).
///
/// The length of a delta is the one its gap lines state, up to 2,147,483,647 residues
/// for a few bytes of input, and NCBI reserves `seq8na` for it at once, then its NCBI4na
/// copy. Where an allocation fails, NCBI throws `std::bad_alloc` (makeblastdb: "BLAST ran
/// out of memory"); LOSAT reserves only the NCBI4na bytes (half as many), and where that
/// fails it is the error of this function, where an allocation that aborts the process
/// would lose that outcome. The space is reserved for all letters, so that no later growth
/// of the vector can fail.
fn fix_bioseq_deltas(
    inst: SeqInst<Ncbi4naBuilder>,
) -> Result<Fixed, std::collections::TryReserveError> {
    Ok(match inst {
        SeqInst::Virtual => (0, None, None),
        SeqInst::TooLong { length } => (length, None, None),
        SeqInst::Raw { length, seq_data } => (length, Some(seq_data), None),
        SeqInst::Delta {
            length,
            residues,
            pieces,
        } => {
            if pieces.first() == Some(&DeltaSeq::LocNull) {
                return Ok((length, None, None));
            }
            if pieces
                .iter()
                .all(|piece| matches!(piece, DeltaSeq::Literal(_)))
            {
                // Without gaps the record's length is its residues' count (below 2^32), and
                // its literals are the residues in order.
                debug_assert_eq!(residues.letters() as u64, length);
                return Ok((length, Some(SeqData::Ncbi4na(residues.into_bytes())), None));
            }
            let size = pieces.iter().fold(0u64, |size, piece| {
                size + match piece {
                    DeltaSeq::LocNull => 1,
                    DeltaSeq::Gap { length } => *length,
                    DeltaSeq::Literal(segment) => u64::from(segment.length),
                }
            });
            // Below 2^32, so a `usize` on every target.
            let letters = (size % (1 << 32)) as usize;
            let mut seq4na = Ncbi4naBuilder::with_letters(letters)?;
            // The letters still to push (`letters` less those pushed), counted down by each piece.
            let mut room = letters;
            for piece in pieces {
                if room == 0 {
                    break;
                }
                room -= match piece {
                    DeltaSeq::LocNull => {
                        seq4na.push_run(0x0f, 1);
                        1
                    }
                    DeltaSeq::Gap { length } => {
                        let run = (length as usize).min(room);
                        seq4na.push_run(0x0f, run);
                        run
                    }
                    DeltaSeq::Literal(segment) => {
                        let run = (segment.length as usize).min(room);
                        // NCBI converts each literal from its coding (`C2naTo8na`, `C4naTo8na`).
                        seq4na.push_ncbi4na(residues.bytes(), segment.start as usize, run);
                        run
                    }
                };
            }
            debug_assert_eq!(seq4na.letters(), letters - room);
            let seq4na = seq4na.into_bytes();
            let fix = (letters as u64 != length).then_some(DeltaFix::Length(letters as u64));
            (length, Some(SeqData::Ncbi4na(seq4na)), fix)
        }
    })
}
