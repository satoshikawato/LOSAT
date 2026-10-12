//! The parts of a `CBioseq` that makeblastdb's FASTA path hands to CWriteDB.
//!
//! `fasta_source` builds this value from a FASTA record as `CFastaReader`,
//! `CWriteDB::ExtractBioseqDeflines` and `s_FixBioseqDeltas` build NCBI's
//! `CBioseq`; the writer reads it as `CWriteDB_Impl` reads the `CBioseq`.

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1005-1038
// ```c++
//         switch(sd.Which()) {
//         case CSeq_data::e_Ncbistdaa:
//             WriteDB_StdaaToBinary(si, m_Sequence);
//             break;
//
//         case CSeq_data::e_Ncbieaa:
//             WriteDB_EaaToBinary(si, m_Sequence);
//             break;
//
//         case CSeq_data::e_Iupacaa:
//             WriteDB_IupacaaToBinary(si, m_Sequence);
//             break;
//
//         case CSeq_data::e_Ncbi2na:
//             WriteDB_Ncbi2naToBinary(si, m_Sequence);
//             break;
//
//         case CSeq_data::e_Ncbi4na:
//             WriteDB_Ncbi4naToBinary(si, m_Sequence, m_Ambig);
//             break;
//
//         case CSeq_data::e_Iupacna:
//              WriteDB_IupacnaToBinary(si, m_Sequence, m_Ambig);
//              break;
//
//         default:
//             msg = "Unable to process sequence for entry [";
//             msg += (m_Bioseq->GetId().front())->GetSeqIdString(false);
//             msg += "].";
//         }
// ```
/// `CSeq_data` of a sequence: the choice and its bytes exactly as NCBI holds
/// them (one letter per residue for the text codings, 4 bases per byte for
/// ncbi2na, 2 bases per byte for ncbi4na; the last byte is padded with zero
/// bits).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SeqData {
    Iupacna(Vec<u8>),
    Iupacaa(Vec<u8>),
    Ncbi2na(Vec<u8>),
    Ncbi4na(Vec<u8>),
    Ncbieaa(Vec<u8>),
    Ncbistdaa(Vec<u8>),
    /// Any other choice (ncbi8na, ncbipna, ncbi8aa, ncbipaa, gap): CWriteDB
    /// stops with "Unable to process sequence for entry [id].".
    Other,
}

/// One sequence of the input as `CBuildDatabase::AddSequences` sees it after
/// `s_FixBioseqDeltas`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bioseq {
    /// The content of the first Seq-id, a local id: `lcl|<local_id>` is
    /// `AsFastaString()` (the "Ignoring sequence 'lcl|2' ..." log line) and
    /// `<local_id>` is `GetSeqIdString(false)` ("Unable to process sequence
    /// for entry [6]."). Without `-parse_seqids` CFastaReader numbers the
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
    /// Seq-inst has no `seq-data` (`!CanGetSeq_data()`).
    pub seq_data: Option<SeqData>,
}
