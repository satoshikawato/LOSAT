//! The parts of a `CBioseq` that makeblastdb's FASTA path hands to CWriteDB.
//!
//! `fasta_source` builds this value from a FASTA record as `CFastaReader`,
//! `CWriteDB::ExtractBioseqDeflines` and `s_FixBioseqDeltas` build NCBI's
//! `CBioseq`; the writer reads it as `CWriteDB_Impl` reads the `CBioseq`.

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1007-1036
/// `CSeq_data` of a sequence: the choice and its bytes exactly as NCBI holds
/// them (one letter per residue for the text codings, 4 bases per byte for
/// ncbi2na, 2 bases per byte for ncbi4na; the FASTA source pads the last byte with zero
/// bits, a PSSM's query has the bytes and the byte count of its file).
///
/// iupacna, iupacaa and ncbistdaa are codings of a PSSM's query only (makeprofiledb), which
/// also gives the other three.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SeqData {
    Ncbi2na(Vec<u8>),
    Ncbi4na(Vec<u8>),
    Iupacna(Vec<u8>),
    Ncbieaa(Vec<u8>),
    Iupacaa(Vec<u8>),
    Ncbistdaa(Vec<u8>),
    /// Any other choice (ncbi8na, ncbipna, ncbi8aa, ncbipaa, gap): CWriteDB
    /// stops with "Unable to process sequence for entry \[id\].".
    Other,
}

use crate::seqid::ParsedSeqId;

/// One sequence of the input as `CBuildDatabase::AddSequences` sees it after
/// `s_FixBioseqDeltas`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bioseq {
    /// The content of the first Seq-id, a local id: `lcl|<local_id>` is
    /// `AsFastaString()` (the "Ignoring sequence 'lcl|2' ..." log line) and
    /// `<local_id>` is `GetSeqIdString` without the version ("Unable to process sequence
    /// for entry \[6\]."). Without `-parse_seqids` CFastaReader numbers the
    /// records of the input from 1.
    pub local_id: String,
    /// The title that `CWriteDB::ExtractBioseqDeflines` puts into the
    /// Blast-def-line, before ASN.1 serialization (the serializer replaces
    /// bytes a VisibleString cannot hold; the writer does that).
    pub title: Vec<u8>,
    /// `CBioseq::IsAa()`.
    pub is_aa: bool,
    /// `CBioseq::GetLength()` (`Seq-inst.length`).
    pub length: u64,
    /// `Seq-inst.seq-data` after `s_FixBioseqDeltas`; `None` when the
    /// Seq-inst has no `seq-data` (`CanGetSeq_data` is false).
    pub seq_data: Option<SeqData>,
    /// With `-parse_seqids`: the parsed ids and what `ExtractBioseqDeflines` reads; the
    /// definition lines are built from them (`deflines`), and `local_id` and `title` are
    /// not used.
    pub parsed: Option<ParsedIds>,
    /// What `s_FixBioseqDeltas` makes of a delta sequence that `length` (its `Seq-inst.length`,
    /// which `AddSequences` tests) does not describe; `None` for every other record (and for a
    /// delta that starts with a gap of unknown length, which NCBI leaves without data:
    /// `seq_data` is `None`).
    pub delta_fix: Option<DeltaFix>,
}

/// See `Bioseq::delta_fix`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeltaFix {
    /// A protein whose gaps stayed a delta (their sum wrapped to 0 past 2^32): NCBI throws
    /// "Protein delta sequences are not supported." (`CMultisourceException`).
    ProteinDelta,
    /// The length of the raw sequence made of the delta, `int length = seq8na.size()`, the sum
    /// of its pieces (a gap of unknown length is one letter) modulo 2^32, which the database
    /// takes; `seq_data` holds that many letters.
    Length(u64),
    /// That sum modulo 2^32 is 2^31 or more, a negative `int` that NCBI converts with as a
    /// `TSeqPos`: LOSAT refuses it (`size` is the sum).
    PastInt { size: u64 },
}

/// The parts of a Bioseq of `CFastaReader` with the flags `fRequireID` and `fParseRawID` that the
/// writer reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParsedIds {
    /// `CBioseq::GetId()`: the ids of the first parse (`CFastaDeflineReader`), or the
    /// generated `lcl|N` of an id with a range; `CBuildDatabase::AddSequences` edits them.
    pub ids: Vec<ParsedSeqId>,
    /// The id that a range replaced: the other row of the Seq-hist alignment
    /// (`CFastaReader::xCreateAlignment`), which `CBioseq::GetNonLocalId` reads.
    pub hist_old_id: Option<ParsedSeqId>,
    /// The `DefLine` field of the `CFastaReader` user object.
    pub defline_field: Vec<u8>,
    /// The title descriptor (`x_ApplyMods`).
    pub title_descr: Option<Vec<u8>>,
}

impl Bioseq {
    /// `GetSeqIdString` without the version, of the first of `GetId()`.
    pub fn first_id_string(&self) -> String {
        match self.parsed.as_ref().and_then(|parsed| parsed.ids.first()) {
            Some(id) => {
                let mut out = Vec::new();
                crate::seqid::label::content(id, false, &mut out);
                String::from_utf8_lossy(&out).into_owned()
            }
            None => self.local_id.clone(),
        }
    }

    /// `GetId().front()->AsFastaString()` ("Unknown" for none), as
    /// `CBuildDatabase::AddSequences` names the sequence in its messages.
    pub fn first_id_fasta(&self) -> String {
        match &self.parsed {
            Some(parsed) => match parsed.ids.first() {
                Some(id) => {
                    let mut out = Vec::new();
                    crate::seqid::label::fasta(id, &mut out);
                    String::from_utf8_lossy(&out).into_owned()
                }
                None => "Unknown".to_string(),
            },
            None => format!("lcl|{}", self.local_id),
        }
    }
}
