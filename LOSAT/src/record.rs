//! The record of a FASTA input that a search reads (`api::FastaRecord`).
//!
//! No NCBI counterpart of the module: it only gives the type a place that rustdoc documents.
//! The type is NCBI's record as `blastinput::fasta_reader` reads it (`CBlastFastaInputSource`,
//! see that module's NCBI reference), and its reader and the engine's methods stay there. rustdoc
//! drops every item of a `#[doc(hidden)]` module such as `blastinput`, and inlines a re-export of
//! a private module's item: `api` and `blastinput::fasta_reader` re-export the type from here
//! (`docs/losat_library_surface_plan.md` §14).

use crate::blastinput::fasta_reader::RecordSource;

/// One record as NCBI's FASTA input source reads it: [`read_fasta`](crate::api::read_fasta) reads
/// them, [`FastaRecord::new`] makes one from its parts.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
#[non_exhaustive]
pub struct FastaRecord {
    /// The generated local ID (`lcl|Query_N`, `lcl|Subject_N`, without `lcl|`).
    pub local_id: String,
    /// The title: the defline after `>` and the white space that follows it, up to the
    /// first byte below 0x20, without white space at its ends (`x_ApplyMods`); empty when
    /// the record has none. Any other bytes, including bytes that are not UTF-8.
    pub title: Vec<u8>,
    /// The residues (IUPAC letters; nucleotide `U` as `T`; `>?` gaps as `N` or `X`), in
    /// lower case where the lowercase mask covers them.
    pub sequence: Vec<u8>,
    /// The reader's messages written while the record was read, in order (each one a
    /// line; NCBI writes them to standard error as it reads).
    pub warnings: Vec<u8>,
    /// What the BLAST archive (outfmt 11) writes besides the letters; `None` (8 bytes, no
    /// allocation) unless the reader was configured to keep it (`ReaderConfig::keep_defline`).
    #[doc(hidden)]
    pub source: Option<Box<RecordSource>>,
}

impl FastaRecord {
    /// A record from its parts: the local ID, the title bytes and the residues, with no reader
    /// messages.
    pub fn new(local_id: impl Into<String>, title: &[u8], sequence: &[u8]) -> Self {
        Self {
            local_id: local_id.into(),
            title: title.to_vec(),
            sequence: sequence.to_vec(),
            warnings: Vec::new(),
            source: None,
        }
    }

    /// The residues.
    pub fn seq(&self) -> &[u8] {
        &self.sequence
    }

    /// The title bytes (empty without one).
    pub fn title(&self) -> &[u8] {
        &self.title
    }
}
