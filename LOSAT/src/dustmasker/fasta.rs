//! dustmasker's FASTA input: `CMaskFastaReader(input_stream, true, args["parse_seqids"])`
//! (`seqmasks_io::fasta_reader`).

use std::io::Read;

use crate::seqmasks_io::fasta_reader::{MaskFastaReader, MaskReadError};

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
                MaskReadError::Refused(what) => Failure::Refused(super::work::unsupported(&what)),
            })?;
        Ok(record.map(|record| Record {
            ids: record.ids,
            title: record.title,
            iupac: record.letters,
        }))
    }
}

/// A failed read of `-in`: the `std::ios_base::failure` of NCBI's stream.
fn read_failure(error: std::io::Error) -> Failure {
    crate::seqmasks_io::run::read_failure(error, super::PROGRAM)
}
