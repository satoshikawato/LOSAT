//! makeblastdb's reading of a FASTA file (`CFastaBioseqSource` of
//! `c++/src/objtools/blast/seqdb_writer/build_db.cpp` over `CFastaReader`).
//!
//! Not implemented yet (stage DB-2, step "FASTA source").

use std::io::Read;

use super::bioseq::Bioseq;

/// An error that ends the reading: the message makeblastdb prints after
/// "BLAST Database creation error: ".
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceError {
    pub message: String,
}

/// A warning line of the reader (without the final newline).
pub type Warn<'a> = &'a mut dyn FnMut(&[u8]);

/// `CFastaBioseqSource`: the sequences of one FASTA input, one at a time.
pub struct FastaBioseqSource<R: Read> {
    input: R,
    is_protein: bool,
}

impl<R: Read> FastaBioseqSource<R> {
    /// `CFastaBioseqSource(fasta_file, is_protein, parse_ids = false, long_ids = false)`.
    pub fn new(input: R, is_protein: bool) -> Self {
        Self { input, is_protein }
    }

    /// `CFastaBioseqSource::GetNext`, then the title of
    /// `CWriteDB::ExtractBioseqDeflines` and `s_FixBioseqDeltas`: the next
    /// sequence, or `None` at the end of the input.
    pub fn next_bioseq(&mut self, _warn: Warn<'_>) -> Result<Option<Bioseq>, SourceError> {
        let _ = (&mut self.input, self.is_protein);
        Err(SourceError {
            message: "the FASTA source is not implemented".to_string(),
        })
    }
}
