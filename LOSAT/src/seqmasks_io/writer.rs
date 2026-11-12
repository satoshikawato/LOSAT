//! The mask writers of NCBI's masking tools: `CMaskWriterInt` (`interval`),
//! `CMaskWriterTabular` (`acclist`), `CMaskWriterFasta` (`fasta`), `CMaskWriterSeqLoc`
//! (`seqloc_asn1_bin`, `seqloc_asn1_text`, `seqloc_xml`) and `CMaskWriterBlastDbMaskInfo`
//! (`maskinfo_asn1_bin`, `maskinfo_asn1_text`, `maskinfo_xml`).
//!
//! A mask is the list of masked intervals of one sequence, 0-based and inclusive, in order
//! (`CSymDustMasker::TMaskList`).

use std::io::{self, Write};

use crate::report::archive::masks::write_ranges_seq_loc;
use crate::report::db_subject::bioseq_ids;
use crate::report::defline::generate_bioseq_defline;
use crate::report::seqid::SeqId;
use crate::report::serial::asn_ber::AsnBinaryEncoder;
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::encoder::{Encoder, SerialError};
use crate::report::serial::xml::XmlEncoder;
use crate::seqid::label::bioseq_fasta;
use crate::seqid::ParsedSeqId;

/// One sequence as the writers see it (`CBioseq_Handle`): its Seq-ids in order (the first is
/// `GetSeqId()`), its title descriptor and its letters in `eCoding_Iupac`.
pub(crate) struct MaskedBioseq<'a> {
    pub ids: &'a [ParsedSeqId],
    pub title: Option<&'a [u8]>,
    pub iupac: &'a [u8],
}

/// The serial formats of the `seqloc_*` and `maskinfo_*` writers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Serial {
    AsnBinary,
    AsnText,
    Xml,
}

/// `Blast-filter-program` (blastdb.asn), an INTEGER with named values.
///
/// NCBI reference (598d8ae6): c++/src/objects/blastdb/blastdb.asn:51-59
/// ```text
/// Blast-filter-program ::= INTEGER {
///     not-set             (0),
///     dust                (10),
///     seg                 (20),
///     windowmasker        (30),
///     repeat              (40),
///     other               (100),
///     max                 (255)
/// }
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum FilterProgram {
    Dust,
}

impl FilterProgram {
    fn name_and_value(self) -> (&'static str, i32) {
        match self {
            FilterProgram::Dust => ("dust", 10),
        }
    }
}

/// Why a sequence could not be written.
#[derive(Debug)]
pub(crate) enum WriteError {
    Io(io::Error),
    /// A case LOSAT does not reproduce: what it is, as a noun phrase that the tool completes
    /// with "is not supported by LOSAT's <program>".
    Refused(String),
}

impl From<io::Error> for WriteError {
    fn from(error: io::Error) -> Self {
        WriteError::Io(error)
    }
}

/// The writer of one `-outfmt`.
pub(crate) enum MaskWriter {
    /// `CMaskWriterInt`.
    Interval,
    /// `CMaskWriterTabular`.
    Tabular,
    /// `CMaskWriterFasta` (`-hard_masking`: `N` instead of lower case).
    Fasta { hard_masking: bool },
    /// `CMaskWriterSeqLoc`.
    SeqLoc(Serial),
    /// `CMaskWriterBlastDbMaskInfo`.
    MaskInfo(MaskInfoWriter),
}

/// `CMaskWriterBlastDbMaskInfo`: the masks of every sequence, written as one
/// `Blast-db-mask-info` when the writer is finished.
pub(crate) struct MaskInfoWriter {
    serial: Serial,
    algo_id: i32,
    program: FilterProgram,
    options: Vec<u8>,
    /// `m_ListOfMasks`: the first Seq-id and the mask of each sequence with a mask.
    masks: Vec<(ParsedSeqId, Vec<(u32, u32)>)>,
}

impl MaskWriter {
    /// `CMaskWriterBlastDbMaskInfo(os, format, algo_id, filt_program, algo_options)`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_blastdb_maskinfo.cpp:45-54
    /// ```c++
    ///     m_BlastDbMaskInfo.Reset(new CBlast_db_mask_info());
    ///     m_BlastDbMaskInfo->SetAlgo_id(algo_id);
    ///     m_BlastDbMaskInfo->SetAlgo_program((int)filt_program);
    ///     m_BlastDbMaskInfo->SetAlgo_options(algo_options);
    /// ```
    pub(crate) fn mask_info(
        serial: Serial,
        algo_id: i32,
        program: FilterProgram,
        options: Vec<u8>,
    ) -> Self {
        MaskWriter::MaskInfo(MaskInfoWriter {
            serial,
            algo_id,
            program,
            options,
            masks: Vec::new(),
        })
    }

    /// `Print(bsh, mask, parsed_id)` of the writer.
    pub(crate) fn print(
        &mut self,
        seq: &MaskedBioseq<'_>,
        mask: &[(u32, u32)],
        parsed_id: bool,
        out: &mut dyn Write,
    ) -> Result<(), WriteError> {
        match self {
            MaskWriter::Interval => {
                // NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_int.cpp:41-56
                // ```c++
                //     PrintId( bsh, parsed_id );
                //     os << endl;
                //     PrintMasks(os, mask);
                // ...
                //     ITERATE(TMaskList, i, mask) {
                //         os << i->first << " - " << i->second << "\n";
                //     }
                // ```
                let mut text = id_to_string(seq, parsed_id)?;
                text.push(b'\n');
                for &(from, to) in mask {
                    text.extend_from_slice(format!("{from} - {to}\n").as_bytes());
                }
                out.write_all(&text)?;
            }
            MaskWriter::Tabular => {
                // NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_tab.cpp:40-49
                // ```c++
                //     const string id = IdToString(bsh, parsed_id);
                //     ITERATE(TMaskList, i, mask) {
                //         os << id << "\t" << i->first << "\t" << i->second << "\n";
                //     }
                // ```
                let id = id_to_string(seq, parsed_id)?;
                let mut text = Vec::new();
                for &(from, to) in mask {
                    text.extend_from_slice(&id);
                    text.extend_from_slice(format!("\t{from}\t{to}\n").as_bytes());
                }
                out.write_all(&text)?;
            }
            MaskWriter::Fasta { hard_masking } => {
                let mut text = id_to_string(seq, parsed_id)?;
                text.push(b'\n');
                fasta_letters(seq.iupac, mask, *hard_masking, &mut text);
                out.write_all(&text)?;
            }
            MaskWriter::SeqLoc(serial) => {
                // NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_seqloc.cpp:63-90
                // ```c++
                //     if (mask.empty()) {
                //         return;
                //     }
                //     ...
                //     CConstRef<CSeq_id> id = bsh.GetSeqId();
                //     CSeq_loc seqloc(const_cast<CSeq_id&>(*id), masked_ranges);
                //     switch (m_OutputFormat) {
                //     case eSerial_AsnBinary:
                //         os << MSerial_AsnBinary << seqloc;
                //     ...
                // ```
                if mask.is_empty() {
                    return Ok(());
                }
                let mut bytes = Vec::new();
                write_document(*serial, &Document::SeqLoc(first_id(seq), mask), &mut bytes)?;
                out.write_all(&bytes)?;
            }
            MaskWriter::MaskInfo(writer) => {
                // NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_blastdb_maskinfo.cpp:138-176
                // ```c++
                // void CMaskWriterBlastDbMaskInfo::Print( const objects::CSeq_id& id,
                //                                         const TMaskList & mask )
                // {
                //     if (mask.empty()) {
                //         return;
                //     }
                //     ...
                //     CRef<CSeq_loc> seqloc(new CSeq_loc(const_cast<CSeq_id&>(id),
                //                                        masked_ranges));
                //     CRef<CBlast_mask_list> mask_list(new CBlast_mask_list);
                //     mask_list->SetMasks().push_back(seqloc);
                //     mask_list->SetMore(true);
                //     m_ListOfMasks.push_back(mask_list);
                // ...
                //     Print(*bsh.GetSeqId(), mask);
                // ```
                if !mask.is_empty() {
                    writer.masks.push((first_id(seq).clone(), mask.to_vec()));
                }
            }
        }
        Ok(())
    }

    /// The destructor of the writer: `CMaskWriterBlastDbMaskInfo` writes its object (one
    /// `Blast-mask-list` with every Seq-loc and `more FALSE`, or an empty one); the others
    /// write nothing.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_blastdb_maskinfo.cpp:89-106
    /// ```c++
    ///     if (m_ListOfMasks.empty()) {
    ///         CRef<CBlast_mask_list> empty_list(new CBlast_mask_list);
    ///         empty_list->SetMasks();
    ///         empty_list->SetMore(false);
    ///         m_ListOfMasks.push_back(empty_list);
    ///     }
    ///
    ///     x_ConsolidateListOfMasks();
    ///     m_BlastDbMaskInfo->SetMasks(*m_ListOfMasks.front());
    ///     s_WriteObject(m_BlastDbMaskInfo, os, m_OutputFormat);
    /// ```
    /// (`x_ConsolidateListOfMasks`, :108-136, moves the Seq-loc of each list into one list and
    /// sets its `more` to false.)
    pub(crate) fn finish(self, out: &mut dyn Write) -> Result<(), WriteError> {
        let MaskWriter::MaskInfo(writer) = self else {
            return Ok(());
        };
        let mut bytes = Vec::new();
        write_document(writer.serial, &Document::MaskInfo(&writer), &mut bytes)?;
        out.write_all(&bytes)?;
        Ok(())
    }
}

/// `bsh.GetSeqId()`: the first Seq-id of the sequence.
fn first_id<'a>(seq: &MaskedBioseq<'a>) -> &'a ParsedSeqId {
    seq.ids.first().expect("a sequence has a Seq-id")
}

/// The intervals as the `Seq-interval` integers (`TSeqPos` stored in an INTEGER).
fn intervals(mask: &[(u32, u32)]) -> Vec<(i32, i32)> {
    mask.iter()
        .map(|&(from, to)| (from as i32, to as i32))
        .collect()
}

/// The object of one document.
enum Document<'a> {
    /// A `Seq-loc` `packed-int` of one sequence.
    SeqLoc(&'a ParsedSeqId, &'a [(u32, u32)]),
    /// The `Blast-db-mask-info` of every sequence.
    MaskInfo(&'a MaskInfoWriter),
}

/// `os << MSerial_<format> << object`: one document in `serial`.
fn write_document(
    serial: Serial,
    document: &Document<'_>,
    out: &mut Vec<u8>,
) -> Result<(), WriteError> {
    let error = match serial {
        Serial::AsnText => encode(&mut AsnTextEncoder::new(out), document),
        Serial::AsnBinary => encode(&mut AsnBinaryEncoder::new(out), document),
        Serial::Xml => {
            let module = match document {
                Document::SeqLoc(..) => "NCBI-Seqloc",
                Document::MaskInfo(_) => "NCBI-BlastDL",
            };
            encode(&mut XmlEncoder::ncbi_object(out, module), document)
        }
    };
    match error {
        None => Ok(()),
        Some(error) => Err(WriteError::Refused(format!(
            "an object on which LOSAT's serial writer stops ({error})"
        ))),
    }
}

/// Writes `document` through `enc` and returns the encoder's error.
fn encode<E: Encoder>(enc: &mut E, document: &Document<'_>) -> Option<SerialError> {
    match document {
        Document::SeqLoc(id, mask) => {
            enc.begin_document("Seq-loc");
            write_ranges_seq_loc(enc, &SeqId::Db(id), &intervals(mask));
        }
        Document::MaskInfo(writer) => {
            enc.begin_document("Blast-db-mask-info");
            enc.begin_sequence("Blast-db-mask-info");
            enc.member("algo-id", 0);
            enc.int(writer.algo_id);
            enc.end_member();
            enc.member("algo-program", 1);
            enc.begin_named("Blast-filter-program");
            let (name, value) = writer.program.name_and_value();
            enc.named_int(name, value);
            enc.end_named();
            enc.end_member();
            enc.member("algo-options", 2);
            enc.visible_str(&writer.options);
            enc.end_member();
            enc.member("masks", 3);
            enc.begin_sequence("Blast-mask-list");
            enc.member("masks", 0);
            enc.begin_list(false);
            for (id, mask) in &writer.masks {
                enc.element();
                write_ranges_seq_loc(enc, &SeqId::Db(id), &intervals(mask));
            }
            enc.end_list();
            enc.end_member();
            enc.member("more", 1);
            enc.boolean(false);
            enc.end_member();
            enc.end_sequence();
            enc.end_member();
            enc.end_sequence();
        }
    }
    enc.end_document();
    enc.error()
}

/// `CMaskWriter::IdToString`: `>`, with `-parse_seqids` the Seq-ids in FASTA form
/// (`GetStringDescr(eFormat_FastA)`) and a space, then `CDeflineGenerator::GenerateDefline`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer.cpp:46-55
/// ```c++
/// string CMaskWriter::IdToString( objects::CBioseq_Handle& bsh, bool parsed_id )
/// {
///     CNcbiOstrstream oss;
///     oss << ">";
///     if( parsed_id ) {
///         oss << CSeq_id::GetStringDescr(*bsh.GetCompleteBioseq(),
///                                        CSeq_id::eFormat_FastA) + " ";
///     }
///     oss << sequence::CDeflineGenerator().GenerateDefline(bsh);
///     return CNcbiOstrstreamToString(oss);
/// }
/// ```
fn id_to_string(seq: &MaskedBioseq<'_>, parsed_id: bool) -> Result<Vec<u8>, WriteError> {
    let mut text = vec![b'>'];
    if parsed_id {
        bioseq_fasta(seq.ids, false, &mut text);
        text.push(b' ');
    }
    let defline = generate_bioseq_defline(
        seq.title.unwrap_or_default(),
        false,
        false,
        &bioseq_ids(seq.ids),
    )
    .map_err(|_| WriteError::Refused("a title whose encoding NCBI cannot tell".to_string()))?;
    if defline.reads_past_end {
        return Err(WriteError::Refused(
            "a title on which NCBI's defline generator reads past its end".to_string(),
        ));
    }
    text.extend_from_slice(&defline.text);
    Ok(text)
}

/// The letters of `CMaskWriterFasta::Print`: lines of 60, the masked letters in lower case
/// (`N` with `-hard_masking`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_fasta.cpp:57-100
/// ```c++
///         string accumulator;
///         TMaskList::const_iterator imask = mask.begin();
///
///         for( TSeqPos i = 0; i < data.size(); ++i )
///         {
///             char letter = data[i];
///
///             if( imask != mask.end() && i >= imask->first ) {
///                 if( i <= imask->second ) {
///                 	if(m_HardMasking) {
///                 		letter = 'N';
///                 	}
///                 	else {
///                 		letter = tolower((unsigned char) letter);
///                 	}
/// 				}
///                 else
///                 {
///                     ++imask;
///
///                     if(imask != mask.end()
///                         && i >= imask->first && i <= imask->second ) {
///                     	... (the same)
/// 					}
///                 }
///             }
///
///             accumulator.append( 1, letter );
///
///             if( !((i + 1)%60) )
///             {
///                 os << accumulator << "\n";
///                 accumulator = "";
///             }
///         }
///
///         if( accumulator.length() ) os << accumulator << "\n";
/// ```
fn fasta_letters(iupac: &[u8], mask: &[(u32, u32)], hard_masking: bool, out: &mut Vec<u8>) {
    let masked = |letter: u8| {
        if hard_masking {
            b'N'
        } else {
            letter.to_ascii_lowercase()
        }
    };
    let mut imask = 0;
    let mut line = 0;
    for (i, &letter) in iupac.iter().enumerate() {
        let i_pos = i as u64;
        let mut letter = letter;
        if let Some(&(first, second)) = mask.get(imask) {
            if i_pos >= u64::from(first) {
                if i_pos <= u64::from(second) {
                    letter = masked(letter);
                } else {
                    imask += 1;
                    if let Some(&(first, second)) = mask.get(imask) {
                        if i_pos >= u64::from(first) && i_pos <= u64::from(second) {
                            letter = masked(letter);
                        }
                    }
                }
            }
        }
        out.push(letter);
        line += 1;
        if line == 60 {
            out.push(b'\n');
            line = 0;
        }
    }
    if line > 0 {
        out.push(b'\n');
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fasta_letters_follow_ncbis_one_step_mask_iterator() {
        // A mask list where one position passes two intervals: NCBI advances the iterator
        // once per position, so the second interval [3, 3] after [0, 1] is reached at 2 and
        // its letter 3 is masked; a third interval right after is reached only at its own
        // positions.
        let mut out = Vec::new();
        fasta_letters(b"ACGTACGT", &[(0, 1), (3, 3), (4, 5)], false, &mut out);
        assert_eq!(out, b"acGtacGT\n");
        let mut out = Vec::new();
        fasta_letters(&[b'A'; 61], &[(59, 60)], true, &mut out);
        assert_eq!(&out[..], [&[b'A'; 59][..], b"N\nN\n"].concat());
    }
}
