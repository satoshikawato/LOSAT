//! The mask writers of NCBI's masking tools: `CMaskWriterInt` (`interval`),
//! `CMaskWriterTabular` (`acclist`), `CMaskWriterFasta` (`fasta`), `CMaskWriterSeqLoc`
//! (`seqloc_asn1_bin`, `seqloc_asn1_text`, `seqloc_xml`) and `CMaskWriterBlastDbMaskInfo`
//! (`maskinfo_asn1_bin`, `maskinfo_asn1_text`, `maskinfo_xml`).
//!
//! A mask is the list of masked intervals of one sequence, 0-based and inclusive, in order
//! (`CSymDustMasker::TMaskList`).

use std::io::{self, Write};

use crate::ncbi_app::Diag;
use crate::report::archive::masks::write_ranges_seq_loc;
use crate::report::db_subject::bioseq_facts;
use crate::report::defline::generate_bioseq_defline;
use crate::report::seqid::SeqId;
use crate::report::serial::asn_ber::AsnBinaryEncoder;
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::encoder::{Encoder, Posted, SerialError, Traced};
use crate::report::serial::xml::XmlEncoder;
use crate::seqdb::mask::FilterProgram;
use crate::seqid::label::bioseq_fasta;
use crate::seqid::ParsedSeqId;

/// One sequence as the writers see it (`CBioseq_Handle`): its Seq-ids in order (the first is
/// `GetSeqId()`), its title descriptor and its letters in `eCoding_Iupac`.
pub(crate) struct MaskedBioseq<'a> {
    pub ids: &'a [ParsedSeqId],
    pub title: Option<&'a [u8]>,
    /// The name of the organism of the Bioseq's first source (a database sequence's), from which
    /// `GenerateDefline` makes an empty title (`report::defline::BioseqFacts::taxname`).
    pub taxname: Option<&'a std::sync::Arc<[u8]>>,
    /// Whether the Bioseq has a source (`report::defline::BioseqFacts::has_source`): a database
    /// protein's title ending in `]` gets its end from it.
    pub has_source: bool,
    pub iupac: &'a [u8],
    /// Whether the sequence is a protein (`CDeflineGenerator` names a protein without a title
    /// `unnamed protein product`).
    pub protein: bool,
}

/// The serial formats of the `seqloc_*` and `maskinfo_*` writers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Serial {
    AsnBinary,
    AsnText,
    Xml,
}

/// Whether the writer of `outfmt` writes through NCBI's serial stream (`CObjectOStream`: the
/// `seqloc_*` and `maskinfo_*` formats), which reports a failed write; the `interval`,
/// `acclist` and `fasta` writers write to the `ostream` itself and go on silently.
pub(crate) fn writes_serial_objects(outfmt: &str) -> bool {
    outfmt.starts_with("seqloc_") || outfmt.starts_with("maskinfo_")
}

/// Why a sequence could not be written.
#[derive(Debug)]
pub(crate) enum WriteError {
    Io(io::Error),
    /// A case LOSAT does not reproduce: what it is, as a noun phrase that the tool completes
    /// with "is not supported by LOSAT's `<program>`".
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
/// `Blast-db-mask-info` when the writer is finished; or, made with the format `interval`
/// (convert2blastmask's), the intervals of each sequence as they are printed.
pub(crate) struct MaskInfoWriter {
    /// `m_OutputFormat`: `None` is `eSerial_None`, the signal of the `interval` format.
    serial: Option<Serial>,
    algo_id: i32,
    program: FilterProgram,
    options: Vec<u8>,
    /// `m_ListOfMasks`: the first Seq-id and the mask of each sequence with a mask.
    masks: Vec<(ParsedSeqId, Vec<(u32, u32)>)>,
}

impl MaskWriter {
    /// `CMaskWriterBlastDbMaskInfo(os, format, algo_id, filt_program, algo_options)`;
    /// `serial` is `None` for the format `interval`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_blastdb_maskinfo.cpp:45-67
    pub(crate) fn mask_info(
        serial: Option<Serial>,
        algo_id: i32,
        program: FilterProgram,
        options: Vec<u8>,
    ) -> Self {
        MaskWriter::MaskInfo(MaskInfoWriter::new(serial, algo_id, program, options))
    }

    /// `Print(bsh, mask, parsed_id)` of the writer; `diag` takes the posts of the serial
    /// stream (`write_document`).
    pub(crate) fn print(
        &mut self,
        seq: &MaskedBioseq<'_>,
        mask: &[(u32, u32)],
        parsed_id: bool,
        out: &mut dyn Write,
        diag: &mut Diag,
    ) -> Result<(), WriteError> {
        match self {
            MaskWriter::Interval => {
                // NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_int.cpp:41-56
                let mut text = id_to_string(seq, parsed_id)?;
                text.push(b'\n');
                interval_lines(mask, &mut text);
                out.write_all(&text)?;
            }
            MaskWriter::Tabular => {
                // NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_tab.cpp:40-49
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
                // NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_seqloc.cpp:61-90
                if mask.is_empty() {
                    return Ok(());
                }
                let mut bytes = Vec::new();
                let document = Document::SeqLoc(first_id(seq), mask);
                write_document(*serial, &document, &mut bytes, out, diag)?;
                out.write_all(&bytes)?;
            }
            MaskWriter::MaskInfo(writer) => {
                // NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_blastdb_maskinfo.cpp:137-176
                writer.print_id(first_id(seq), mask, out)?;
            }
        }
        Ok(())
    }

    /// The destructor of the writer: `CMaskWriterBlastDbMaskInfo` writes its object (one
    /// `Blast-mask-list` with every Seq-loc and `more FALSE`, or an empty one); the others
    /// write nothing.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_blastdb_maskinfo.cpp:89-106
    /// (`x_ConsolidateListOfMasks`, :108-136, moves the Seq-loc of each list into one list and
    /// sets its `more` to false.)
    /// (`m_OutputFormat` `eSerial_None`, the format `interval`, writes nothing either.)
    pub(crate) fn finish(self, out: &mut dyn Write, diag: &mut Diag) -> Result<(), WriteError> {
        let MaskWriter::MaskInfo(writer) = self else {
            return Ok(());
        };
        let Some(serial) = writer.serial else {
            return Ok(());
        };
        let mut bytes = Vec::new();
        write_document(serial, &Document::MaskInfo(&writer), &mut bytes, out, diag)?;
        out.write_all(&bytes)?;
        Ok(())
    }
}

impl MaskInfoWriter {
    /// `CMaskWriterBlastDbMaskInfo(os, format, algo_id, filt_program, algo_options)`, cited at
    /// `MaskWriter::mask_info`.
    pub(crate) fn new(
        serial: Option<Serial>,
        algo_id: i32,
        program: FilterProgram,
        options: Vec<u8>,
    ) -> Self {
        MaskInfoWriter {
            serial,
            algo_id,
            program,
            options,
            masks: Vec::new(),
        }
    }

    /// `CMaskWriterBlastDbMaskInfo::Print(id, mask)`: nothing for an empty mask; with the
    /// format `interval`, the intervals at once (`CMaskWriterInt::PrintMasks`, no id); else
    /// the mask is kept for the object that `finish` writes.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_blastdb_maskinfo.cpp:137-150
    pub(crate) fn print_id(
        &mut self,
        id: &ParsedSeqId,
        mask: &[(u32, u32)],
        out: &mut dyn Write,
    ) -> Result<(), WriteError> {
        if mask.is_empty() {
            return Ok(());
        }
        if self.serial.is_none() {
            let mut text = Vec::new();
            interval_lines(mask, &mut text);
            out.write_all(&text)?;
            return Ok(());
        }
        self.masks.push((id.clone(), mask.to_vec()));
        Ok(())
    }
}

/// `CMaskWriterInt::PrintMasks`: one line `from - to` for each interval.
///
/// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_int.cpp:49-56
/// The lines are added to `text`.
fn interval_lines(mask: &[(u32, u32)], text: &mut Vec<u8>) {
    for &(from, to) in mask {
        // Writing to a `Vec<u8>` does not fail.
        let _ = writeln!(text, "{from} - {to}");
    }
}

/// `bsh.GetSeqId()`: the first Seq-id of the sequence.
fn first_id<'a>(seq: &MaskedBioseq<'a>) -> &'a ParsedSeqId {
    seq.ids.first().expect("a sequence has a Seq-id")
}

/// The object of one document.
enum Document<'a> {
    /// A `Seq-loc` `packed-int` of one sequence.
    SeqLoc(&'a ParsedSeqId, &'a [(u32, u32)]),
    /// The `Blast-db-mask-info` of every sequence.
    MaskInfo(&'a MaskInfoWriter),
}

/// `os << MSerial_<format> << object`: one document in `serial` appended to `out`. A byte that
/// NCBI's stream writes as `#` with an error message (`ReplacedChar`) makes the encoder report
/// `NonVisibleChar`; such a document is written again through `Traced` and its messages are
/// posted (`rewrite_and_post`, which writes the start of the document to `stream` before them),
/// before the caller writes the rest of `out` to `stream`.
fn write_document(
    serial: Serial,
    document: &Document<'_>,
    out: &mut Vec<u8>,
    stream: &mut dyn Write,
    diag: &mut Diag,
) -> Result<(), WriteError> {
    let start = out.len();
    let mut error = match serial {
        Serial::AsnText => encode(&mut AsnTextEncoder::new(out), document),
        Serial::AsnBinary => encode(&mut AsnBinaryEncoder::new(out), document),
        Serial::Xml => encode(
            &mut XmlEncoder::ncbi_object(out, module(document)),
            document,
        ),
    };
    if let Some(SerialError::NonVisibleChar(_)) = error {
        error = rewrite_and_post(serial, document, out, start, stream, diag)?;
    }
    match error {
        None => Ok(()),
        Some(error) => Err(WriteError::Refused(format!(
            "an object on which LOSAT's serial writer stops ({error})"
        ))),
    }
}

/// The DTD module of the XML form of `document`.
fn module(document: &Document<'_>) -> &'static str {
    match document {
        Document::SeqLoc(..) => "NCBI-Seqloc",
        Document::MaskInfo(_) => "NCBI-BlastDL",
    }
}

/// `document` written again from `start` of `out` through `Traced`, and NCBI's posts of the
/// bytes that its stream replaced (`Diag::post_replaced`), which come while the document is
/// written, NCBI's `cerr` tied to `cout`: each after the output written before and the bytes of
/// the document that NCBI's serial buffer (`COStreamBuffer`, `Traced`) had written to the output
/// by then (`Posted::written`), which go to `stream` first and leave `out`. The error of the
/// second writing, if any.
#[cold]
#[inline(never)]
fn rewrite_and_post(
    serial: Serial,
    document: &Document<'_>,
    out: &mut Vec<u8>,
    start: usize,
    stream: &mut dyn Write,
    diag: &mut Diag,
) -> io::Result<Option<SerialError>> {
    out.truncate(start);
    let (error, posted) = match serial {
        Serial::AsnText => encode_traced(Traced::new(AsnTextEncoder::new(out), start), document),
        Serial::AsnBinary => {
            encode_traced(Traced::new(AsnBinaryEncoder::new(out), start), document)
        }
        Serial::Xml => encode_traced(
            Traced::xml(XmlEncoder::ncbi_object(out, module(document)), start),
            document,
        ),
    };
    if error.is_none() {
        let mut handed = start;
        for posted in &posted {
            let written = start + posted.written as usize;
            if written > handed {
                stream.write_all(&out[handed..written])?;
                handed = written;
            }
            let _ = stream.flush();
            diag.post_replaced(std::slice::from_ref(&posted.replaced));
        }
        out.drain(start..handed);
    }
    Ok(error)
}

/// `encode` through `enc` and the replaced bytes.
fn encode_traced<E: Encoder>(
    mut enc: Traced<E>,
    document: &Document<'_>,
) -> (Option<SerialError>, Vec<Posted>) {
    let error = encode(&mut enc, document);
    (error, enc.take_posted())
}

/// Writes `document` through `enc` and returns the encoder's error.
fn encode<E: Encoder>(enc: &mut E, document: &Document<'_>) -> Option<SerialError> {
    match document {
        Document::SeqLoc(id, mask) => {
            enc.begin_document("Seq-loc");
            write_ranges_seq_loc(enc, &SeqId::Db(id), mask);
        }
        Document::MaskInfo(writer) => {
            enc.begin_document("Blast-db-mask-info");
            enc.begin_sequence("Blast-db-mask-info");
            enc.member("algo-id", 0);
            enc.int(writer.algo_id);
            enc.end_member();
            enc.member("algo-program", 1);
            enc.begin_named("Blast-filter-program");
            enc.named_int(writer.program.name(), writer.program.value());
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
                write_ranges_seq_loc(enc, &SeqId::Db(id), mask);
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
/// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer.cpp:46-65
fn id_to_string(seq: &MaskedBioseq<'_>, parsed_id: bool) -> Result<Vec<u8>, WriteError> {
    let mut text = vec![b'>'];
    if parsed_id {
        bioseq_fasta(seq.ids, false, &mut text);
        text.push(b' ');
    }
    let mut facts = bioseq_facts(seq.ids);
    facts.taxname = seq.taxname.cloned();
    facts.has_source = seq.has_source;
    let defline =
        generate_bioseq_defline(seq.title.unwrap_or_default(), seq.protein, false, &facts)
            .map_err(|_| {
                WriteError::Refused("a title whose encoding NCBI cannot tell".to_string())
            })?;
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
/// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_fasta.cpp:57-101
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
