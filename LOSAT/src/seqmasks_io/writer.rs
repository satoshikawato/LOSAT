//! The mask writers of NCBI's masking tools: `CMaskWriterInt` (`interval`),
//! `CMaskWriterTabular` (`acclist`), `CMaskWriterFasta` (`fasta`), `CMaskWriterSeqLoc`
//! (`seqloc_asn1_bin`, `seqloc_asn1_text`, `seqloc_xml`) and `CMaskWriterBlastDbMaskInfo`
//! (`maskinfo_asn1_bin`, `maskinfo_asn1_text`, `maskinfo_xml`).
//!
//! A mask is the list of masked intervals of one sequence, 0-based and inclusive, in order
//! (`CSymDustMasker::TMaskList`).

use std::io;

use crate::ncbi_app::Diag;
use crate::report::archive::masks::write_ranges_seq_loc;
use crate::report::defline::bioseq_facts;
use crate::report::defline::generate_bioseq_defline;
use crate::report::ostream::{CutBuffer, ReportWrite};
use crate::report::seqid::SeqId;
use crate::report::serial::asn_ber::AsnBinaryEncoder;
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::encoder::{Encoder, Posted, SerialError, Traced};
use crate::report::serial::xml::XmlEncoder;
use crate::seqdb::mask::FilterProgram;
use crate::seqid::label::bioseq_fasta;
use crate::seqid::ParsedSeqId;
use crate::utils::decimal::push_decimal;

/// One sequence as the writers see it (`CBioseq_Handle`): its Seq-ids in order, its title
/// descriptor and its letters in `eCoding_Iupac`.
pub(crate) struct MaskedBioseq<'a> {
    pub ids: &'a [ParsedSeqId],
    /// `GetSeqId()` where a handle that outlives the sequence made the info of the sequence's
    /// access handle ([`access_seq_id`]): windowmasker's id lists matched as Seq-ids
    /// (`windowmasker::ids::IdLists::handle_seq_id`). `None`: the access Seq-id as it is.
    pub handle_seq_id: Option<ParsedSeqId>,
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
    Fasta {
        hard_masking: bool,
        /// Reusable buffer across records.
        buffer: CutBuffer,
    },
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
    /// `m_ListOfMasks`: the Seq-id that names each sequence with a mask (`seq_id`) and its mask.
    masks: Vec<(ParsedSeqId, Vec<(u32, u32)>)>,
}

impl MaskWriter {
    /// `CMaskWriterFasta` (`-hard_masking`: `N` instead of lower case).
    pub(crate) fn fasta(hard_masking: bool) -> Self {
        MaskWriter::Fasta {
            hard_masking,
            buffer: CutBuffer::default(),
        }
    }

    /// The `CMaskWriterBlastDbMaskInfo` constructor (stream, format, algorithm id, filtering
    /// program and algorithm options);
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

    /// `Print` of the writer, given `bsh`, the mask and `parsed_id`, into `out` in NCBI's
    /// insertions and flushes (`CutBuffer`: a text on standard error between two records comes
    /// after NCBI's last write, `ncbi_app::Diag`); `diag` takes the posts of the serial stream
    /// (`write_document`).
    ///
    /// The definition line of the `interval` and `fasta` writers is one insertion followed by
    /// `endl` (a `put` of the line feed and a flush). Each interval line is four insertions (the
    /// two numbers, ` - ` and the line feed); each `acclist` line six (the id, a tab, a number, a
    /// tab, a number, the line feed); each FASTA line of letters two (the letters, the line feed).
    pub(crate) fn print(
        &mut self,
        seq: &MaskedBioseq<'_>,
        mask: &[(u32, u32)],
        parsed_id: bool,
        out: &mut dyn ReportWrite,
        diag: &mut Diag,
    ) -> Result<(), WriteError> {
        match self {
            MaskWriter::Interval => {
                // NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_int.cpp:41-56
                let mut text = CutBuffer::default();
                text.insertion(&id_to_string(seq, parsed_id)?);
                text.endl();
                interval_lines(mask, &mut text);
                text.write_to(out)?;
            }
            MaskWriter::Tabular => {
                // NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_tab.cpp:40-49
                let id = id_to_string(seq, parsed_id)?;
                let mut text = CutBuffer::default();
                let mut number = Vec::with_capacity(10);
                for &(from, to) in mask {
                    text.insertion(&id);
                    text.insertion(b"\t");
                    number.clear();
                    push_decimal(&mut number, i64::from(from));
                    text.insertion(&number);
                    text.insertion(b"\t");
                    number.clear();
                    push_decimal(&mut number, i64::from(to));
                    text.insertion(&number);
                    text.insertion(b"\n");
                }
                text.write_to(out)?;
            }
            MaskWriter::Fasta {
                hard_masking,
                buffer,
            } => {
                let id = id_to_string(seq, parsed_id)?;
                fasta_letters(&id, seq.iupac, mask, *hard_masking, buffer, out)?;
            }
            MaskWriter::SeqLoc(serial) => {
                // NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_seqloc.cpp:61-90
                if mask.is_empty() {
                    return Ok(());
                }
                let mut bytes = Vec::new();
                let document = Document::SeqLoc(seq_id(seq), mask);
                write_document(*serial, &document, &mut bytes, out, diag)?;
                // The object stream's buffer, then its `Close`'s flush of the stream.
                out.write_all(&bytes)?;
                out.ncbi_flush()?;
            }
            MaskWriter::MaskInfo(writer) => {
                // NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_blastdb_maskinfo.cpp:137-176
                writer.print_id(seq_id(seq), mask, out)?;
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
    pub(crate) fn finish(
        self,
        out: &mut dyn ReportWrite,
        diag: &mut Diag,
    ) -> Result<(), WriteError> {
        let MaskWriter::MaskInfo(writer) = self else {
            return Ok(());
        };
        let Some(serial) = writer.serial else {
            return Ok(());
        };
        let mut bytes = Vec::new();
        write_document(serial, &Document::MaskInfo(&writer), &mut bytes, out, diag)?;
        // The object stream's buffer, then its `Close`'s flush of the stream.
        out.write_all(&bytes)?;
        out.ncbi_flush()?;
        Ok(())
    }
}

impl MaskInfoWriter {
    /// The `CMaskWriterBlastDbMaskInfo` constructor (stream, format, algorithm id, filtering
    /// program and algorithm options), cited at
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

    /// `CMaskWriterBlastDbMaskInfo::Print` of an id and a mask: nothing for an empty mask; with the
    /// format `interval`, the intervals at once (`CMaskWriterInt::PrintMasks`, no id); else
    /// the mask is kept for the object that `finish` writes.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_blastdb_maskinfo.cpp:137-150
    pub(crate) fn print_id(
        &mut self,
        id: &ParsedSeqId,
        mask: &[(u32, u32)],
        out: &mut dyn ReportWrite,
    ) -> Result<(), WriteError> {
        if mask.is_empty() {
            return Ok(());
        }
        if self.serial.is_none() {
            let mut text = CutBuffer::default();
            interval_lines(mask, &mut text);
            text.write_to(out)?;
            return Ok(());
        }
        self.masks.push((id.clone(), mask.to_vec()));
        Ok(())
    }
}

/// `CMaskWriterInt::PrintMasks`: one line `from - to` for each interval, four insertions (the
/// numbers, ` - `, the line feed).
///
/// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_int.cpp:49-56
/// The lines are added to `text`.
fn interval_lines(mask: &[(u32, u32)], text: &mut CutBuffer) {
    let mut number = Vec::with_capacity(10);
    for &(from, to) in mask {
        number.clear();
        push_decimal(&mut number, i64::from(from));
        text.insertion(&number);
        text.insertion(b" - ");
        number.clear();
        push_decimal(&mut number, i64::from(to));
        text.insertion(&number);
        text.insertion(b"\n");
    }
}

/// `bsh.GetSeqId()`: the Seq-id of the sequence's access handle, as the handles that outlive the
/// sequence give it back (`MaskedBioseq::handle_seq_id`).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/bioseq_handle.cpp:128-131
fn seq_id<'s>(seq: &'s MaskedBioseq<'_>) -> &'s ParsedSeqId {
    seq.handle_seq_id
        .as_ref()
        .unwrap_or_else(|| access_seq_id(seq.ids))
}

/// The Seq-id of `CBioseq_Handle::GetAccessSeq_id_Handle()` of a sequence whose Seq-ids are `ids`
/// (a handle of `CBioseq_CI` has no Seq-id of its own): its first gi, else its first Textseq-id,
/// else its first Seq-id. A gi of 0 is not a gi there (`IsGi()` asks for a non-zero packed
/// value, and gi 0's handle packs 0): `>gi|0|gb|AB000001.1|` is named `genbank`, `>lcl|x|gi|0`
/// `local`, and only a sequence whose first id is gi 0 and that has no Textseq-id `gi 0`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/bioseq_handle.cpp:939-985
/// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:2102-2106
/// NCBI reference (598d8ae6): c++/include/objects/seq/seq_id_handle.hpp:279-282
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:609-627
pub(crate) fn access_seq_id(ids: &[ParsedSeqId]) -> &ParsedSeqId {
    ids.iter()
        .find(|id| matches!(id, ParsedSeqId::Gi(gi) if *gi != 0))
        .or_else(|| ids.iter().find(|id| id.textseq().is_some()))
        .or_else(|| ids.first())
        .expect("a sequence has a Seq-id")
}

/// The object of one document.
enum Document<'a> {
    /// A `Seq-loc` `packed-int` of one sequence.
    SeqLoc(&'a ParsedSeqId, &'a [(u32, u32)]),
    /// The `Blast-db-mask-info` of every sequence.
    MaskInfo(&'a MaskInfoWriter),
}

/// The stream `os` given `MSerial_<format>` and then the object: one document in `serial` appended
/// to `out`. A byte that NCBI's stream writes as `#` with an error message (`ReplacedChar`) makes
/// the encoder report `NonVisibleChar`; such a document is written again through `Traced` and its
/// messages are posted (`rewrite_and_post`, which writes the start of the document to `stream`
/// before them), before the caller writes the rest of `out` to `stream`.
fn write_document(
    serial: Serial,
    document: &Document<'_>,
    out: &mut Vec<u8>,
    stream: &mut dyn ReportWrite,
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
/// written: each after the bytes of the document that NCBI's serial buffer (`COStreamBuffer`,
/// `Traced`) had written to the output by then (`Posted::written`, one insertion), which go to
/// `stream` first and leave `out`; the post makes NCBI's writes of the output first (`cerr`'s tie,
/// `Diag`). The error of the second writing, if any.
#[cold]
#[inline(never)]
fn rewrite_and_post(
    serial: Serial,
    document: &Document<'_>,
    out: &mut Vec<u8>,
    start: usize,
    stream: &mut dyn ReportWrite,
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
            diag.post_replaced(std::slice::from_ref(&posted.replaced));
            if diag.closed() {
                return Err(io::ErrorKind::BrokenPipe.into());
            }
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
/// (`GetStringDescr` with `eFormat_FastA`, which ranks the ids for a protein when the Bioseq is one:
/// segmasker's `lcl|a|sp|P12345|` is `sp|P12345|`) and a space, then
/// `CDeflineGenerator::GenerateDefline`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer.cpp:46-65
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2347-2375
fn id_to_string(seq: &MaskedBioseq<'_>, parsed_id: bool) -> Result<Vec<u8>, WriteError> {
    let mut text = vec![b'>'];
    if parsed_id {
        bioseq_fasta(seq.ids, seq.protein, &mut text);
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

/// `CMaskWriterFasta::Print` of the definition line `id` and the letters `iupac`: the line and
/// `endl`, then lines of 60 letters, the masked letters in lower case (`N` with `-hard_masking`),
/// each line two insertions (the letters, then a line feed of one byte).
///
/// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_fasta.cpp:46-101
/// The lines go to `out` a block at a time (each line is one of NCBI's insertions into its
/// buffered stream, `CutBuffer`), so the text of a record is never held whole. Each line is copied
/// whole and then masked in place ([`mask_line`]): the loop over the letters touches only the
/// line, not the growing text.
fn fasta_letters(
    id: &[u8],
    iupac: &[u8],
    mask: &[(u32, u32)],
    hard_masking: bool,
    text: &mut CutBuffer,
    out: &mut dyn ReportWrite,
) -> io::Result<()> {
    const LINE: usize = 60;
    const BLOCK: usize = 1 << 16;
    text.clear();
    text.insertion(id);
    text.endl();
    text.write_to(out)?;
    if iupac.is_empty() {
        return Ok(());
    }
    let mut imask = 0;
    let mut pos = 0u64;
    let seq_bytes = iupac.len() + iupac.len().div_ceil(LINE);
    // After the id line's `endl()`, NCBI's buffer is empty; if the entire sequence body fits
    // within `ncbi_room()`, no line insertion can trigger a write, so per-line cuts are skipped.
    if seq_bytes < out.ncbi_room() {
        text.bytes.reserve(seq_bytes);
        for letters in iupac.chunks(LINE) {
            let at = text.bytes.len();
            text.bytes.extend_from_slice(letters);
            mask_line(&mut text.bytes[at..], pos, mask, &mut imask, hard_masking);
            text.bytes.push(b'\n');
            pos += letters.len() as u64;
        }
        return out.write_all(&text.bytes);
    }
    for letters in iupac.chunks(LINE) {
        let at = text.bytes.len();
        text.bytes.extend_from_slice(letters);
        mask_line(&mut text.bytes[at..], pos, mask, &mut imask, hard_masking);
        text.end_insertion();
        text.insertion(b"\n");
        pos += letters.len() as u64;
        if text.bytes.len() >= BLOCK {
            text.write_to(out)?;
        }
    }
    text.write_to(out)
}

/// NCBI's walk over the mask for the letters of one line, the first at `pos`: a letter at or
/// past the start of interval `imask` is masked while within it; past its end, the walk takes
/// the next interval (one step per letter) and masks the letter if it is within that one
/// (NCBI's loop cited at [`fasta_letters`]).
fn mask_line(
    line: &mut [u8],
    pos: u64,
    mask: &[(u32, u32)],
    imask: &mut usize,
    hard_masking: bool,
) {
    // A letter before the start of interval `imask` changes neither itself nor `imask`, so a
    // line that ends before that start is left as it is.
    match mask.get(*imask) {
        None => return,
        Some(&(first, _)) if pos + line.len() as u64 <= u64::from(first) => return,
        Some(_) => {}
    }
    let masked = |letter: u8| {
        if hard_masking {
            b'N'
        } else {
            letter.to_ascii_lowercase()
        }
    };
    let mut i_pos = pos;
    for letter in line {
        if let Some(&(first, second)) = mask.get(*imask) {
            if i_pos >= u64::from(first) {
                if i_pos <= u64::from(second) {
                    *letter = masked(*letter);
                } else {
                    *imask += 1;
                    if let Some(&(first, second)) = mask.get(*imask) {
                        if i_pos >= u64::from(first) && i_pos <= u64::from(second) {
                            *letter = masked(*letter);
                        }
                    }
                }
            }
        }
        i_pos += 1;
    }
}

#[cfg(test)]
#[path = "writer_tests.rs"]
mod tests;
