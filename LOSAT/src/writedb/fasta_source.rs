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

mod line_reader;
mod pack;
mod reader;
#[cfg(test)]
mod tests;
mod title;

use std::io::Read;

use super::bioseq::{Bioseq, SeqData};
use reader::{DeltaSeq, FastaReader, SeqInst};

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
    /// `m_FastaReader` with its `m_LineReader`; `None` once `GetNext` has reset the line
    /// reader (the end of the input) or thrown.
    reader: Option<FastaReader<R>>,
}

impl<R: Read> FastaBioseqSource<R> {
    /// `CFastaBioseqSource(fasta_file, is_protein, parse_ids = false, long_ids = false)`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:701-734
    /// ```c++
    /// CFastaBioseqSource::CFastaBioseqSource(CNcbiIstream & fasta_file,
    ///                                        bool is_protein,
    ///                                        bool parse_ids,
    ///                                        bool long_ids)
    ///     : m_FastaReader(NULL)
    /// {
    ///     m_LineReader.Reset(new CBufferedLineReader(fasta_file));
    ///     typedef CFastaReader::EFlags TFlags;
    ///
    ///     int iflags = CFastaReader::fForceType;
    ///
    ///     if (is_protein) {
    ///         iflags |= CFastaReader::fAssumeProt;
    ///     } else {
    ///         iflags |= CFastaReader::fAssumeNuc;
    ///         iflags |=  CFastaReader::fParseGaps;
    ///     }
    ///
    ///     if (parse_ids) {
    ///         iflags |= CFastaReader::fRequireID;
    ///         // parse bare accessions
    ///         if (!long_ids) {
    ///             iflags |= CFastaReader::fParseRawID;
    ///         }
    ///     } else {
    ///         iflags |= CFastaReader::fNoParseID;
    ///     }
    ///
    ///     iflags |= CFastaReader::fQuickIDCheck;
    ///     iflags |= CFastaReader::fDisableNoResidues;
    ///
    ///     TFlags flags = (TFlags) iflags;
    ///
    ///     m_FastaReader = new CFastaReader(*m_LineReader, flags);
    /// ```
    /// The flags are those of `reader`.
    pub fn new(input: R, is_protein: bool) -> Self {
        Self {
            reader: Some(FastaReader::new(input, is_protein)),
        }
    }

    /// `CFastaBioseqSource::GetNext`, then the title of
    /// `CWriteDB::ExtractBioseqDeflines` and `s_FixBioseqDeltas`: the next
    /// sequence, or `None` at the end of the input.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:755-792
    /// ```c++
    /// CConstRef<CBioseq> CFastaBioseqSource::GetNext()
    /// {
    ///     CConstRef<CBioseq> rv;
    ///
    ///     if (m_LineReader.NotEmpty() && ! m_LineReader->AtEOF()) {
    ///         CRef<CSeq_entry> entry;
    ///         try { entry = m_FastaReader->ReadOneSeq(); }
    ///         catch (const CObjReaderParseException& e) {
    ///             static const string kKeyword("m_Pos = ");
    ///             SIZE_TYPE start = NStr::Find(e.what(), kKeyword);
    ///             SIZE_TYPE end = NStr::Find(e.what(), ")", start);
    ///             string pos("unknown");
    ///             if (start != NPOS && end != NPOS) {
    ///                 start += kKeyword.size();
    ///                 pos = string(e.what()).substr(start, end-start);
    ///             }
    ///             string msg = e.GetMsg();
    ///             const string extra_string("CFastaReader: ");
    ///             if (NStr::StartsWith(msg, extra_string)) {
    ///                 msg.erase(0, extra_string.size());
    ///             }
    ///             NCBI_THROW(CWriteDBException, eFileErr, msg);
    ///         }
    ///
    ///         if (entry.NotEmpty()) {
    ///             _ASSERT(entry->IsSeq());
    ///             rv.Reset(& entry->GetSeq());
    ///         }
    ///     }
    ///
    ///     // Any failure to read a Bioseq is considered an EOF.
    ///
    ///     if (rv.Empty()) {
    ///         m_LineReader.Reset();
    ///     }
    ///
    ///     return rv;
    /// }
    /// ```
    /// `x_EditAndAddBioseq` then takes the deflines and fixes the deltas:
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:473-488
    /// ```c++
    ///     CRef<CBlast_def_line_set> headers =
    ///         CWriteDB::ExtractBioseqDeflines(*bs, m_ParseIDs, m_LongIDs,
    ///                 m_ScanBioseq4CFastaReaderUsrObjct);
    ///
    ///     x_EditHeaders(headers);
    ///
    ///     // Add the sequence
    ///     if (sv) {
    ///         m_OutputDb->AddSequence(*bs, *sv);
    ///     } else {
    ///         bs = s_FixBioseqDeltas(bs);
    ///         if(bs->GetInst().CanGetSeq_data())
    ///         	m_OutputDb->AddSequence(*bs);
    ///         else
    ///         	return false;
    ///     }
    /// ```
    /// `AddSequences` tests the length of the sequence before that call, and
    /// `s_FixBioseqDeltas` can make a sequence longer (a gap of unknown length becomes one
    /// `N`) but not empty or the reverse, so `Bioseq::length` is the length after it. A
    /// record above 2,147,483,647 residues is returned with its length and no sequence
    /// data: NCBI stops at that length, before `x_EditAndAddBioseq`. (NCBI counts in 32
    /// bits and a record of 4,294,967,296 residues or more wraps there; LOSAT counts on.)
    ///
    /// A failed read of the input ends the reading with "Read error" and the system's
    /// message (LOSAT's text: NCBI's line reader throws a `CIOException` "Read error").
    /// After an error the source returns `None`.
    pub fn next_bioseq(&mut self, warn: Warn<'_>) -> Result<Option<Bioseq>, SourceError> {
        let Some(reader) = &mut self.reader else {
            return Ok(None);
        };
        if reader.lines.at_eof() {
            let read_error = reader.lines.take_read_error();
            self.reader = None;
            return match read_error {
                Some(error) => Err(read_error_message(&error)),
                None => Ok(None),
            };
        }
        let record = reader.read_one_seq(warn);
        if let Some(error) = reader.lines.take_read_error() {
            self.reader = None;
            return Err(read_error_message(&error));
        }
        let record = match record {
            Ok(record) => record,
            Err(error) => {
                self.reader = None;
                let message = error.message;
                let message = message
                    .strip_prefix("CFastaReader: ")
                    .map(str::to_string)
                    .unwrap_or(message);
                return Err(SourceError { message });
            }
        };
        let title = match extract_title(&record.defline_field, record.title_descr) {
            Ok(title) => title,
            Err(message) => {
                self.reader = None;
                return Err(SourceError {
                    message: message.to_string(),
                });
            }
        };
        let is_aa = reader.is_protein();
        let (length, seq_data) = fix_bioseq_deltas(record.inst);
        Ok(Some(Bioseq {
            local_id: record.local_id.to_string(),
            title,
            is_aa,
            length,
            seq_data,
        }))
    }
}

fn read_error_message(error: &std::io::Error) -> SourceError {
    SourceError {
        message: format!("Read error: {error}"),
    }
}

/// The title of the Blast-def-line that `CWriteDB::ExtractBioseqDeflines(bs, parse_ids =
/// false, long_ids = false, scan = true)` builds for a record of `CFastaReader`: the one
/// of the reader's user object (`title::fasta_reader_title`), or, when that object has an
/// empty definition line, the title descriptor (`x_BuildDeflinesFromBioseq`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:802-829
/// ```c++
///         if (bin_hdr.empty()) {
/// 	    try {
///             x_GetFastaReaderDeflines(*bioseq,
///                                      deflines,
///                                      membbits,
///                                      linkouts,
///                                      pig,
///                                      false,
///                                      parse_ids,
///                                      long_ids,
///                                      scan_bioseq_4_cfastareader_usrobj);
/// 	      } catch ( const CSeqIdException&e ) {
/// 		//LOG_POST(Info << "x_GetFastaReaderDeflines " << e.GetMsg() );
/// 	    }
///         }
///
///         if(!s_UseFastaReaderDeflines(bioseq, deflines, long_ids)) {
///         	deflines.Reset();
///         }
///
///         if (bin_hdr.empty() && deflines.Empty()) {
///             x_BuildDeflinesFromBioseq(*bioseq,
///                                       deflines,
///                                       membbits,
///                                       linkouts,
///                                       pig);
///         }
///     }
/// ```
/// `s_UseFastaReaderDeflines` keeps the deflines of the user object, because the record
/// has only a local ID:
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:736-750
/// ```c++
/// static bool s_UseFastaReaderDeflines(CConstRef<CBioseq> & bioseq, CConstRef<CBlast_def_line_set> & deflines, bool long_seqid)
/// {
/// 	if(deflines.Empty())
/// 		return false;
///
/// 	const CSeq_id * bioseq_id = bioseq->GetNonLocalId();
///
/// 	if(bioseq_id == NULL ||
///        // For bare pir and prf ids go with the one from defline.
///        // This is to parse bare ids as local ones. The bare pdb ids are pdb in
///        // bioseq (parsed by CFastaReader), but local in deflines (parsed by
///        // CSeq_id).
///        (!long_seqid && (bioseq_id->IsPrf() || bioseq_id->IsPir()))) {
///         return true;
///        }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:610-613,694
/// ```c++
///         if (desc.IsTitle()) {
///             //defline->SetTitle((**iter)->GetTitle());
///             titles = (**iter).GetTitle();
///         }
/// ...
///         defline->SetTitle(titles);
/// ```
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
/// ```c++
/// static CConstRef<CBioseq> s_FixBioseqDeltas(CConstRef<objects::CBioseq> bs)
/// {
///     if ((! bs->CanGetInst()) || bs->GetInst().CanGetSeq_data() ||
///     	 ! bs->GetInst().CanGetExt() || ! bs->GetInst().GetExt().IsDelta()) {
///         return bs;
///     }
/// ...
///     try {
///         const CDelta_ext & dext = bs->GetInst().GetExt().GetDelta();
///
///         if(dext.Get().front()->Which() != CDelta_seq::e_Literal)
///         	return bs;
///
/// ...
///         string seq8na;
///         if (bs->GetInst().CanGetLength()) {
///             seq8na.reserve(bs->GetInst().GetLength());
///         }
///
///         string na8;
///
///         ITERATE(TItems, item, dext.Get()) {
///         	if(((**item).IsLoc()) && ((**item).GetLoc().IsNull())) {
///                 seq8na.append(1, 0x0f);
///         		continue;
///         	}
///
///             const CSeq_literal & L = (**item).GetLiteral();
///
///             if (!L.CanGetSeq_data()) {
///                 if (L.CanGetLength()){
/// ...
///             if (L.GetSeq_data().IsNcbi2na()) {
///                 CSeqConvert::Convert(L.GetSeq_data().GetNcbi2na(),
///                                      CSeqUtil::e_Ncbi2na,
///                                      0,
///                                      L.GetLength(),
///                                      na8,
///                                      CSeqUtil::e_Ncbi8na);
///             } else if (L.GetSeq_data().IsNcbi4na()) {
///                 CSeqConvert::Convert(L.GetSeq_data().GetNcbi4na(),
///                                      CSeqUtil::e_Ncbi4na,
///                                      0,
///                                      L.GetLength(),
///                                      na8,
///                                      CSeqUtil::e_Ncbi8na);
/// 	    } else if (L.GetSeq_data().IsGap()) {
///                 seq8na.append(L.GetLength(), 0x0f);
///             } else {
///                 NCBI_THROW(CMultisourceException, eArg,
///                            "Unhandled type of sequence data encountered.");
///             }
///
///             seq8na += na8;
///             na8.resize(0);
///         }
///
///         // Now convert back to 4na, since WriteDB does not yet handle
///         // 8na sequences.
///
///         int length = seq8na.size();
///         vector<char> seq4na;
///         CSeqConvert::Convert(seq8na,
///                              CSeqUtil::e_Ncbi8na,
///                              0,
///                              length,
///                              seq4na,
///                              CSeqUtil::e_Ncbi4na);
///
///         // Copy the needed fields of the CBioseq (but remove the delta
///         // sequence) and add a Seq-data.
///
/// ...
///         inst->SetSeq_data().SetNcbi4na().Set().swap(seq4na);
///         inst->SetMol(CSeq_inst::eMol_na);
///         inst->SetLength(length);
///         inst->SetRepr(CSeq_inst::eRepr_raw);
/// ```
/// A sequence with `seq-data` (raw) or without a delta (virtual) is returned as it is.
/// The delta of a protein does not occur (`reader::SeqInst`), nor a literal without a
/// length or a `Seq-loc` other than null. A delta whose first piece is not a literal (a
/// gap of unknown length) is returned as it is, without `seq-data`. The NCBI8na codes of
/// a literal are those of its letters (`pack::append_ncbi8na`).
fn fix_bioseq_deltas(inst: SeqInst) -> (u64, Option<SeqData>) {
    match inst {
        SeqInst::Virtual => (0, None),
        SeqInst::TooLong { length } => (length, None),
        SeqInst::Raw { length, seq_data } => (length, Some(seq_data)),
        SeqInst::Delta {
            length,
            residues,
            pieces,
        } => {
            if pieces.first() == Some(&DeltaSeq::LocNull) {
                return (length, None);
            }
            let mut seq8na = Vec::with_capacity(length as usize);
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
            let seq4na = pack::ncbi8na_to_ncbi4na(&seq8na);
            (seq8na.len() as u64, Some(SeqData::Ncbi4na(seq4na)))
        }
    }
}
