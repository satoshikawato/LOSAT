//! The first line of a record that `CBlastInputReader` tries as a Seq-id.
//!
//! NCBI reads such a line as the ID of a sequence to fetch with its data loaders (GenBank
//! over the network, or a BLAST database), which LOSAT does not do; it is an explicit
//! rejection. A line that NCBI's `CSeq_id` parser rejects as a malformatted ID is read as
//! FASTA data, as NCBI does.

use super::ReaderConfig;

/// LOSAT's rejection of a first line that NCBI may read as a Seq-id, or `None` when NCBI
/// reads it as FASTA (`CSeq_id` throws "Malformatted ID").
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:131-151
/// ```c++
///         const string line = NStr::TruncateSpaces_Unsafe(*++GetLineReader());
///         if ( !line.empty() && isalnum(line.data()[0]&0xff) ) {
///             try {
///                 CRef<CSeq_id> id(new CSeq_id(line, (CSeq_id::fParse_AnyRaw |
/// 							CSeq_id::fParse_ValidLocal)));
/// 		if (id->IsLocal()  &&  !NStr::StartsWith(line, "lcl|") ) {
///                     // Expected to throw an exception.
///                     id.Reset(new CSeq_id(line));
/// 		}
///                 CRef<CBioseq> bioseq(x_CreateBioseq(id));
///                 CRef<CSeq_entry> retval(new CSeq_entry());
///                 retval->SetSeq(*bioseq);
///                 return retval;
///             } catch (const CSeqIdException& e) {
///                 if (NStr::Find(e.GetMsg(), "Malformatted ID") != NPOS) {
///                     // This is probably just plain fasta, so just
///                     // defer to CFastaReader
///                 } else {
///                     throw;
///                 }
/// ```
/// LOSAT reads the line as FASTA only where `CSeq_id` certainly fails: a line of letters
/// only (no accession or GI is made of letters only, and a local ID made of letters is
/// retried without `fParse_ValidLocal`, which fails). Any other line is rejected.
pub(super) fn reject_seq_id_line(line: &[u8], config: &ReaderConfig) -> Option<anyhow::Error> {
    if line.iter().all(u8::is_ascii_alphabetic) {
        return None;
    }
    Some(anyhow::anyhow!(
        "the first line of the {role} ({:?}) is not a defline, and NCBI BLAST+ may read it as a sequence identifier to fetch from GenBank or a BLAST database, which is not supported by LOSAT's {program} (start the {role} with a '>' defline, or set DATA_LOADERS=none in the [BLAST] section of .ncbirc)",
        String::from_utf8_lossy(line),
        role = config.role,
        program = config.program,
    ))
}
