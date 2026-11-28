//! The Bioseq of a database sequence as NCBI's formatters read it (`CSeqDB::GetBioseq`: the
//! FASTA formatter and `CBlastDB_BioseqFormatter` call it twice for each sequence), `-outfmt %b`
//! (the Bioseq as an ASN.1 text document) and the ASN.1 text documents of the retrieval (`%b`,
//! and `%d` of `format::FormatSpec`) with NCBI's posts of the bytes its text stream replaces.

use std::sync::Arc;

use crate::report::seqid::SeqId;
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::bioseq::{
    write_bioseq, OrgSource, SeqData, SeqInst, Seqdesc, UserData, UserObject,
};
use crate::report::serial::encoder::{Encoder, ReplacedChar, SerialError, Traced};
use crate::seqdb::taxdb::SourceNames;
use crate::seqid::accguide::parse_flags::ANY_RAW;
use crate::seqid::{DeflineSet, ParsedSeqId};

use super::retrieve::{parse_target_with, HeaderSet, Output, Source, Stop};

/// What NCBI's two `GetBioseq` calls of one `Write` of `oid`, `config` and `target` read of the
/// sequence `oid` (`with`): the first `GetBioseq` of `oid`, then `GetBioseq` of `oid`, `gi` and
/// `seqid` for a `-target_only` target or `GetBioseq` of `oid` again.
pub(super) struct DbBioseq<'s> {
    /// `GetHdr` of `oid`: every definition line the database keeps.
    pub set: HeaderSet<'s>,
    /// The user object's header (`x_GetFilteredBinaryHeader`), read for `%b` alone.
    pub header: Vec<u8>,
    /// The lines of the Bioseq's ids and title for a target: the first that holds it.
    kept: Option<DeflineSet>,
    /// The target's Seq-id.
    pub target: Option<ParsedSeqId>,
    /// `GetSeqLength` of `oid`, at least 1 (NCBI makes no Bioseq of no letters).
    pub length: u32,
}

impl DbBioseq<'_> {
    /// The definition lines of the Bioseq's ids (the first line's) and title: the target's line,
    /// or every line.
    pub fn shown(&self) -> &DeflineSet {
        self.kept.as_ref().unwrap_or(&self.set)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1306-1473
    /// NCBI's two `GetBioseq` calls for the sequence `oid` and the target text `target`, whose
    /// Seq-id `target_id` reads: the texts of each call's taxonomy lookups posted to `out` for a
    /// Bioseq that NCBI makes (one with definition lines and letters), before any byte of the
    /// sequence, then `write` of that Bioseq; nothing where NCBI gives no Bioseq (no definition
    /// line, or no letters). With `with_header` the user object's header is read too (`%b`).
    /// The Bioseq is lent to `write`, not returned, so each formatter's copy of this function
    /// builds it in place: a returned one cost a copy of it and a call for each sequence of a
    /// FASTA run (callgrind against dev, line RG-Dbc DECISIONS 10).
    #[inline]
    pub fn with<'a, 'b>(
        source: &Source,
        oid: u32,
        with_header: bool,
        target: Option<&[u8]>,
        target_id: fn(&[u8]) -> Result<ParsedSeqId, Stop>,
        out: &mut Output<'a, 'b>,
        write: impl FnOnce(&DbBioseq<'_>, &mut Output<'a, 'b>) -> Result<(), Stop>,
    ) -> Result<(), Stop> {
        let read = if with_header {
            source.deflines_with_header(oid)?
        } else {
            source.deflines(oid)?.map(|set| (set, Vec::new()))
        };
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1328-1348
        // NCBI's null set (a stored header of no bytes): no Bioseq, and the second call's read
        // of a target's line throws (probes za, zb, zd of line RG-Dbc).
        let Some((set, header)) = read else {
            if let Some(target) = target {
                target_id(target)?;
                return Err(super::retrieve::null_pointer());
            }
            return Ok(());
        };
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1362-1386
        // `GetBioseq` looks the taxids up for a Bioseq it gives, one with definition lines and
        // letters, the first call for `oid` and the second for the target or `oid` again.
        let length = if set.deflines.is_empty() {
            0i64
        } else {
            source.length(oid)?
        };
        let made = length >= 1;
        if !made && !set.deflines.is_empty() {
            source.check_non_positive_protein_range(oid)?;
        }
        if made {
            source.post_bioseq_texts(oid, &set, None, out)?;
        }
        // `GetBioseq` of `oid`, `gi` and `seqid`: the target's line for the ids and the title,
        // every line the filter keeps for the sources (`x_GetTaxonomy` of `x_GetTaxDefline`, the
        // target's line first).
        let target = match target {
            Some(target) => Some(target_id(target)?),
            None => None,
        };
        let kept = match &target {
            Some(target) => Some(target_set(&set, target)?),
            None => None,
        };
        if made {
            source.post_bioseq_texts(oid, &set, target.as_ref(), out)?;
        }
        // A front line without Seq-ids (makeblastdb writes none; an edited header) makes a Bioseq
        // of no ids, as NCBI's: each formatter reads it as NCBI's does (`fasta::short_header`
        // writes no header line, a range or a mask of it throws, `%b` writes the empty list;
        // line RG D1c-1c-b).
        if !made || kept.as_ref().unwrap_or(&set).deflines.is_empty() {
            return Ok(());
        }
        let bioseq = DbBioseq {
            set,
            header,
            kept,
            target,
            length: length as u32,
        };
        write(&bioseq, out)
    }
}

/// The definition lines NCBI's `GetBioseq` keeps for `-target_only` (the Bioseq's ids and
/// title): the first that holds the target id `target_id` (`s_SeqDB_SeqIdIn`), or the local id
/// of its label.
fn target_set(set: &DeflineSet, target_id: &ParsedSeqId) -> Result<DeflineSet, Stop> {
    let Some(index) = set.bioseq_defline(target_id) else {
        return Err(Stop::Ncbi(
            "Error: oid headers do not contain target gi/seq_id.".to_string(),
        ));
    };
    Ok(DeflineSet {
        deflines: vec![set.deflines[index].clone()],
        utf8: set.utf8,
    })
}

/// A document of ASN.1 text that the retrieval writes (`Output::asn_text`).
pub(super) trait AsnTextDocument {
    fn write<E: Encoder>(&self, enc: &mut E);
}

/// `CBlastDB_BioseqFormatter` (`-outfmt %b`).
pub(super) struct BioseqFormatter;

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:553-580
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2264-2283
/// The Seq-id of a `-target_only` target as `CBlastDB_BioseqFormatter` gets it: `CSeq_id` of `text`
/// (`fParse_AnyRaw`: a number is a GI), whose GI `SeqDB_SimplifySeqid` gives as one
/// (`GetBioseq` of `oid` and `gi`).
fn bioseq_target_id(text: &[u8]) -> Result<ParsedSeqId, Stop> {
    parse_target_with(text, ANY_RAW)
}

impl BioseqFormatter {
    /// Writes the Bioseq of sequence `oid` as an ASN.1 text document to `out` (nothing for a
    /// sequence NCBI's `GetBioseq` gives no Bioseq for), after the texts of the taxonomy lookups
    /// of NCBI's two `GetBioseq` calls.
    pub fn write(
        &self,
        source: &Source,
        oid: u32,
        target: Option<&[u8]>,
        out: &mut Output<'_, '_>,
    ) -> Result<(), Stop> {
        DbBioseq::with(
            source,
            oid,
            true,
            target,
            bioseq_target_id,
            out,
            |bioseq, out| write_db_bioseq(source, oid, bioseq, out),
        )
    }
}

/// The ASN.1 text document of the Bioseq `bioseq` of sequence `oid` (`BioseqFormatter::write`).
fn write_db_bioseq(
    source: &Source,
    oid: u32,
    bioseq: &DbBioseq,
    out: &mut Output<'_, '_>,
) -> Result<(), Stop> {
    let shown = bioseq.shown();
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:955-1018
    let title = shown.bioseq_title();
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1099-1189
    // `x_GetTaxonomy`: a source for each positive taxid of the lines, the target's first,
    // with the names of the taxonomy database where its lookup finds the taxid.
    let sources = bioseq
        .set
        .taxonomy_lookups(bioseq.target.as_ref())
        .into_iter()
        .map(|taxid| Ok((taxid, source.tax_names(taxid)?)))
        .collect::<Result<Vec<(i32, Arc<SourceNames>)>, Stop>>()?;
    let data = source.seq_data(oid)?;
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1269-1303
    // `x_GetAsnDefline`: the header is never empty here (a Bioseq is made of a stored header
    // with bytes, `DbBioseq::with`), so the user object is always written.
    let mut descr = vec![
        Seqdesc::Title(&title),
        Seqdesc::User(UserObject {
            type_str: ASN_DEFLINE_LABEL,
            label: ASN_DEFLINE_LABEL,
            num: Some(1),
            data: UserData::Oss(&bioseq.header),
        }),
    ];
    for (taxid, names) in &sources {
        let names = names.names();
        descr.push(Seqdesc::Source(OrgSource {
            taxname: names.map(|names| names.scientific_name.as_slice()),
            common: names.map(|names| names.common_name.as_slice()),
            taxid: *taxid,
        }));
    }
    let ids: Vec<SeqId<'_>> = shown.deflines[0].ids.iter().map(SeqId::Db).collect();
    let inst = SeqInst {
        protein: source.protein,
        length: bioseq.length,
        data: Some(match &data {
            DbSeqData::Ncbi2na(bytes) => SeqData::Ncbi2na(bytes),
            DbSeqData::Ncbi4na(bytes) => SeqData::Ncbi4na(bytes),
            DbSeqData::Ncbistdaa(bytes) => SeqData::Ncbistdaa(bytes),
        }),
    };
    out.asn_text(&BioseqDocument {
        ids: &ids,
        descr: &descr,
        inst: &inst,
    })
}

/// `kAsnDeflineObjLabel`, the type and the field label of the user object that holds the
/// definition lines' header.
const ASN_DEFLINE_LABEL: &[u8] = b"ASN1_BlastDefLine";

/// The `Seq-data` of `GetBioseq` (`Source::seq_data`).
pub(super) enum DbSeqData {
    Ncbi2na(Vec<u8>),
    Ncbi4na(Vec<u8>),
    Ncbistdaa(Vec<u8>),
}

/// The document that `MSerial_AsnText` writes for the Bioseq.
struct BioseqDocument<'a> {
    ids: &'a [SeqId<'a>],
    descr: &'a [Seqdesc<'a>],
    inst: &'a SeqInst<'a>,
}

impl AsnTextDocument for BioseqDocument<'_> {
    fn write<E: Encoder>(&self, enc: &mut E) {
        enc.begin_document("Bioseq");
        write_bioseq(enc, self.ids, self.descr, self.inst);
        enc.end_document();
    }
}

// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:413-444
// NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:843-908 (the object stream's buffer
// writes to the output when it is full and at `Close`; the object stream's `Close` flushes the
// output first, `Output::asn_text`)
/// The ASN.1 text of `document` as NCBI's stream writes it to the output: the writes of its
/// buffer (`Traced::pieces`, each one insertion; a document that fits the buffer is one), or,
/// where a VisibleString holds a byte that NCBI's stream replaces by `#` (a name of the taxonomy
/// database outside ' ' ..= '~'), written again through `Traced`, which gives `post` NCBI's post
/// of each replaced byte after the bytes NCBI's buffer had written to `out` by then (`cerr` is
/// tied to `cout`; `Err` from `post` where that met a closed pipe). Err: the closed pipe of a
/// write, or LOSAT's refusal of another error of the stream.
pub(super) fn write_asn_text<D: AsnTextDocument>(
    document: &D,
    out: &mut dyn std::io::Write,
    mut post: impl FnMut(&ReplacedChar) -> std::io::Result<()>,
) -> Result<(), Stop> {
    let mut buf = Vec::new();
    let mut enc = AsnTextEncoder::new(&mut buf);
    document.write(&mut enc);
    match enc.error() {
        None if buf.len() <= OBJECT_BUFFER => return out.write_all(&buf).map_err(Stop::write),
        None => return write_pieces(document, &buf, out).map_err(Stop::write),
        Some(SerialError::NonVisibleChar(_)) => {}
        Some(error) => return Err(serial_refusal(error)),
    }
    buf.clear();
    let mut enc = Traced::new(AsnTextEncoder::new(&mut buf), 0)
        .posting(|replaced: &ReplacedChar| post(replaced));
    document.write(&mut enc);
    if let Some(error) = enc.error() {
        return Err(serial_refusal(error));
    }
    enc.flush_to(out, true).map_err(Stop::write)
}

/// The size of NCBI's object stream's buffer (`COStreamBuffer`): a document of at most this many
/// bytes reaches the output in one write, at `Close`.
const OBJECT_BUFFER: usize = 4096;

/// The document `document`, whose text is `text`, larger than the object stream's buffer: the
/// writes of NCBI's buffer to `out` (`Traced::pieces`, from a second encoding), then the rest at
/// `Close`.
#[cold]
#[inline(never)]
fn write_pieces<D: AsnTextDocument>(
    document: &D,
    text: &[u8],
    out: &mut dyn std::io::Write,
) -> std::io::Result<()> {
    let mut again = Vec::new();
    let mut enc = Traced::new(AsnTextEncoder::new(&mut again), 0);
    document.write(&mut enc);
    let mut at = 0;
    for piece in enc.pieces() {
        let end = (at + piece.size).min(text.len());
        out.write_all(&text[at..end])?;
        at = end;
    }
    if at < text.len() {
        out.write_all(&text[at..])?;
    }
    Ok(())
}

/// LOSAT's refusal of an error of the text stream that NCBI's stream does not raise.
fn serial_refusal(error: SerialError) -> Stop {
    Stop::refused(format!(
        "an ASN.1 text output that LOSAT does not write ({error})"
    ))
}
