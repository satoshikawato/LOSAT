//! makeblastdb's reading of a FASTA file (`CFastaBioseqSource` of
//! `c++/src/objtools/blast/seqdb_writer/build_db.cpp` over `CFastaReader`), without
//! `-parse_seqids`.
//!
//! A record becomes the `Bioseq` that `CBuildDatabase::AddSequences` adds: read by
//! `CFastaReader::ReadOneSeq` with makeblastdb's flags (`reader`, over NCBI's
//! `CBufferedLineReader`, `line_reader`), with the title that
//! `CWriteDB::ExtractBioseqDeflines` gives its Blast-def-line (`title`) and the sequence
//! data that `s_FixBioseqDeltas` leaves (`fix_bioseq_deltas`; the reader packs nucleotide
//! residues with `CSeqConvert::Pack`, `pack`).
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

use super::bioseq::{Bioseq, ParsedIds, PastSeqPos, SeqData};
use reader::{
    DeltaSeq, FastaReader, PastSeqPosForm, ReadError, ReaderFlags, ReaderWarning, SeqInst,
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
    reader: Option<FastaReader<R>>,
}

impl<R: Read> FastaBioseqSource<R> {
    /// `CFastaBioseqSource` of the file, `is_protein` and `parse_ids`, with `long_ids` false: with
    /// `parse_ids` (`-parse_seqids`) the reader parses the ids of each definition line
    /// (flags `fRequireID` and `fParseRawID`) and the Bioseq carries them (`Bioseq::parsed`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:701-734
    /// The flags are cited at `reader::ReaderFlags::fasta_bioseq_source`.
    pub fn with_parse_ids(input: R, is_protein: bool, parse_ids: bool) -> Self {
        Self {
            reader: Some(FastaReader::new(
                input,
                ReaderFlags::fasta_bioseq_source(is_protein, parse_ids),
            )),
        }
    }

    /// `CFastaBioseqSource::GetNext`, then the title of
    /// `CWriteDB::ExtractBioseqDeflines` and `s_FixBioseqDeltas`: the next
    /// sequence, or `None` at the end of the input.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:755-792
    /// `x_EditAndAddBioseq` then takes the deflines and fixes the deltas:
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:473-488
    /// `AddSequences` tests the length of the sequence before that call, and
    /// `s_FixBioseqDeltas` can make a sequence longer (a gap of unknown length becomes one
    /// `N`) but not empty or the reverse, so `Bioseq::length` is the length after it. A
    /// record above 2,147,483,647 residues is returned with its length and no sequence
    /// data: NCBI stops at that length, before `x_EditAndAddBioseq`. NCBI counts in 32
    /// bits: a record of 4,294,967,296 residues or more (gaps included) has the length of
    /// its count modulo 2^32 (`reader::SeqInst::PastSeqPos`, `Bioseq::past_seq_pos`).
    ///
    /// A failed read of the input ends the reading with "Read error" and the system's
    /// message (LOSAT's text: NCBI's line reader throws a `CIOException` "Read error").
    /// After an error the source returns `None`.
    pub fn next_bioseq(&mut self, warn: Warn<'_>) -> Result<Option<Bioseq>, SourceError> {
        let Some(reader) = &mut self.reader else {
            return Ok(None);
        };
        if reader.at_eof() {
            let read_error = reader.take_read_error();
            self.reader = None;
            return match read_error {
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
            self.reader = None;
            return Err(read_error_message(&error));
        }
        let record = match record {
            Ok(record) => record,
            Err(ReadError::OutOfMemory) => {
                self.reader = None;
                return Err(SourceError::out_of_memory());
            }
            Err(ReadError::Ncbi(exception)) => {
                // Only convert2blastmask's reader throws other exceptions.
                self.reader = None;
                return Err(SourceError {
                    message: exception.get_msg().to_vec(),
                    kind: SourceErrorKind::Creation,
                });
            }
            Err(ReadError::Parse(error)) => {
                self.reader = None;
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
                    self.reader = None;
                    return Err(SourceError {
                        message: message.to_string().into_bytes(),
                        kind: SourceErrorKind::Creation,
                    });
                }
            }
        };
        // A record past NCBI's 32-bit length comes with NCBI's (wrapped) length and no data:
        // `AddSequences` skips it or stops at its length first, as any record
        // (`CBuildDatabase::add_sequences`), and `s_FixBioseqDeltas` then throws for a protein
        // delta, returns a delta that starts with `Seq-loc null` unchanged (no data), or
        // builds the data, which LOSAT refuses (`PastSeqPos::Refused`).
        let (length, seq_data, past_seq_pos) = match record.inst {
            SeqInst::PastSeqPos {
                length,
                seq_pos_length,
                form,
            } => {
                let past = match form {
                    PastSeqPosForm::Delta | PastSeqPosForm::DeltaFromNull if is_aa => {
                        Some(PastSeqPos::ProteinDelta)
                    }
                    PastSeqPosForm::DeltaFromNull => None,
                    PastSeqPosForm::Delta | PastSeqPosForm::Residues => {
                        Some(PastSeqPos::Refused { length })
                    }
                };
                (seq_pos_length, None, past)
            }
            inst => {
                let Ok((length, seq_data)) = fix_bioseq_deltas(inst) else {
                    self.reader = None;
                    return Err(SourceError::out_of_memory());
                };
                (length, seq_data, None)
            }
        };
        Ok(Some(Bioseq {
            local_id: record.local_id.to_string(),
            title,
            is_aa,
            length,
            seq_data,
            parsed,
            past_seq_pos,
        }))
    }
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

/// `s_FixBioseqDeltas`: the length and the sequence data of a record.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:290-295,306-311,318-334,343-382,395-398
/// A sequence with `seq-data` (raw) or without a delta (virtual) is returned as it is.
/// The delta of a protein does not occur (`reader::SeqInst`), nor a literal without a
/// length or a `Seq-loc` other than null. A delta whose first piece is not a literal (a
/// gap of unknown length) is returned as it is, without `seq-data`. The NCBI8na codes of
/// a literal are those of its letters (`pack::append_ncbi8na`).
///
/// The length of a delta is the one its gap lines state, up to 2,147,483,647 residues
/// for a few bytes of input, and NCBI reserves `seq8na` for it at once. Where that
/// allocation or the one of the NCBI4na copy fails, NCBI throws `std::bad_alloc`
/// (makeblastdb: "BLAST ran out of memory"); here it is the error of this function,
/// where an allocation that aborts the process would lose that outcome. The space is
/// reserved for all pieces, so that no later growth of the vector can fail.
fn fix_bioseq_deltas(
    inst: SeqInst,
) -> Result<(u64, Option<SeqData>), std::collections::TryReserveError> {
    Ok(match inst {
        SeqInst::Virtual => (0, None),
        SeqInst::TooLong { length } => (length, None),
        // `next_bioseq` hands a record past NCBI's `TSeqPos` over without calling this.
        SeqInst::PastSeqPos { seq_pos_length, .. } => (seq_pos_length, None),
        SeqInst::Raw { length, seq_data } => (length, Some(seq_data)),
        SeqInst::Delta {
            length,
            residues,
            pieces,
        } => {
            if pieces.first() == Some(&DeltaSeq::LocNull) {
                return Ok((length, None));
            }
            let size = pieces.iter().fold(0usize, |size, piece| {
                size.saturating_add(match piece {
                    DeltaSeq::LocNull => 1,
                    DeltaSeq::Gap { length } => *length as usize,
                    DeltaSeq::Literal(segment) => segment.length as usize,
                })
            });
            let mut seq8na = Vec::new();
            seq8na.try_reserve_exact(size)?;
            for piece in pieces {
                match piece {
                    DeltaSeq::LocNull => seq8na.push(0x0f),
                    DeltaSeq::Gap { length } => {
                        seq8na.resize(seq8na.len() + length as usize, 0x0f);
                    }
                    DeltaSeq::Literal(segment) => {
                        let start = segment.start as usize;
                        pack::append_ncbi8na(
                            &mut seq8na,
                            &residues[start..start + segment.length as usize],
                        );
                    }
                }
            }
            let seq4na = pack::ncbi8na_to_ncbi4na(&seq8na)?;
            (seq8na.len() as u64, Some(SeqData::Ncbi4na(seq4na)))
        }
    })
}
