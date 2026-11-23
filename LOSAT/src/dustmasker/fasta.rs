//! dustmasker's FASTA input: a `CMaskFastaReader` over the input stream, with `true` and the
//! `parse_seqids` argument (`seqmasks_io::fasta_reader`).

use std::io::Read;

use crate::seqmasks_io::fasta_reader::{
    seq_map_exception, MaskFastaReader, MaskReadError, SeqMapError,
};

use super::work::{Failure, Record, SequenceSource};

/// The reader of `input`, or the failure of its first read (NCBI's line reader reads its first
/// buffer in the constructor).
pub(super) fn source<'a, R: Read + 'a>(
    input: R,
    parse_seqids: bool,
) -> Result<Box<dyn SequenceSource + 'a>, Failure> {
    let reader = MaskFastaReader::new(input, false, parse_seqids).map_err(read_failure)?;
    Ok(Box::new(FastaSource { reader }))
}

struct FastaSource<R: Read> {
    reader: MaskFastaReader<R>,
}

impl<R: Read> SequenceSource for FastaSource<R> {
    fn next_record(&mut self) -> Result<Option<Record>, Failure> {
        // dustmasker posts at NCBI's release level (Error): the reader's warnings are not
        // written.
        let record = self
            .reader
            .next(&mut |_warning| {})
            .map_err(|error| match error {
                MaskReadError::Ncbi(exception) => Failure::Ncbi(exception),
                MaskReadError::Read(error) => read_failure(error),
                MaskReadError::OutOfMemory => Failure::Std("std::bad_alloc".to_string()),
                MaskReadError::Std(what) => Failure::Std(what),
                // NCBI's dustmasker crashes on a raw sequence of 4294967295 letters (SIGSEGV,
                // survey 11d dm4).
                MaskReadError::SeqMap {
                    error: SeqMapError::InvalidRawLength,
                    id,
                    ..
                } => Failure::Refused(super::work::unsupported(&format!(
                    "a FASTA record of 4,294,967,295 letters, counting the letters of its gap \
                     lines ({id}; NCBI's object manager takes that length for an invalid \
                     position)"
                ))),
                MaskReadError::SeqMap { error, .. } => Failure::Ncbi(seq_map_exception(error)),
            })?;
        Ok(record.map(|record| Record {
            ids: record.ids,
            title: record.title,
            taxname: None,
            has_source: false,
            iupac: record.letters,
            texts: Vec::new(),
        }))
    }
}

/// A failed read of `-in`: the `std::ios_base::failure` of NCBI's stream.
fn read_failure(error: std::io::Error) -> Failure {
    crate::seqmasks_io::run::read_failure(error, super::PROGRAM)
}
