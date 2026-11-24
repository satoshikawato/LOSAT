//! Header records of a database volume (`.phr` / `.nhr`).
//!
//! Each sequence has one BER-encoded `Blast-def-line-set`. Two decoders read it:
//! `decode_full` gives every member and every Seq-id (`seqid::ber`), and the
//! decoders below give the part a database built without `-parse_seqids` uses,
//! the title and the `gnl|BL_ORD_ID|n` id that makeblastdb gives the sequence
//! (`n` its number in the volume), without building the Seq-ids. For any other
//! Seq-id, a `gnl|BL_ORD_ID|n` with another number included, they return
//! `SeqDbError::Unsupported` (`PARSED_SEQIDS`): their callers read the ordinal
//! ids only.
//!
//! ASN.1 definitions (NCBI 598d8ae6):
//! `c++/src/objects/blastdb/blastdb.asn:31-47` (`Blast-def-line-set`,
//! `Blast-def-line`), `c++/src/objects/seqloc/seqloc.asn:24-45` (`Seq-id`),
//! `c++/src/objects/general/general.asn:71-80` (`Dbtag`, `Object-id`).

use std::path::Path;

use super::SeqDbError;
use crate::seqid::ber::{
    fix_visible_chars, int_list, BerError, Reader, DEFLINE_LAST_MEMBER, DEFLINE_MEMBERS,
    DEFLINE_SEQID, DEFLINE_TAXID, DEFLINE_TITLE, TAG_CONTEXT_CONSTRUCTED, TAG_SEQUENCE,
};
use crate::seqid::DeflineSet;

/// One definition line of a sequence in a database without parsed Seq-ids.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Defline {
    /// The title: the whole FASTA definition line after `>`.
    pub title: Vec<u8>,
    /// `n` of `gnl|BL_ORD_ID|n` as stored, that is, relative to the volume: the number of the
    /// sequence in its volume (another number is a parsed Seq-id).
    pub volume_ordinal: u32,
    /// Number of stored title bytes outside ' ' ..= '~' that were read as '#'
    /// (`fix_visible_chars`). NCBI posts an error message for each of them every
    /// time it reads the definition line; makeblastdb stores none.
    pub replaced_title_bytes: usize,
    /// `taxid` of the definition line (NCBI's `TTaxId`, an `int`), 0 when the line has none or
    /// the decode did not read it (`decode_deflines_as` without `want_taxid`).
    pub taxid: i32,
    /// Whether the record (of every definition line) holds a UTF8String where a VisibleString
    /// is expected, which NCBI reads with its warning once per process
    /// (`seqid::ber::utf8_string_path` names the member).
    pub utf8: bool,
}

/// `general` is the eleventh alternative of `Seq-id` (seqloc.asn:25-35).
const SEQID_GENERAL: u8 = 10;
/// Members of `Dbtag` and alternatives of `Object-id` (general.asn:71-80).
const DBTAG_DB: u8 = 0;
const DBTAG_TAG: u8 = 1;
const OBJECT_ID_ID: u8 = 0;

/// Why the ordinal-id decoders stop.
enum HeaderError {
    /// The record is not a valid encoding.
    Ber(BerError),
    /// The record holds a Seq-id other than the `gnl|BL_ORD_ID|n` of the sequence's own number
    /// in the volume (`PARSED_SEQIDS`).
    ParsedSeqIds,
}

impl From<BerError> for HeaderError {
    fn from(error: BerError) -> Self {
        HeaderError::Ber(error)
    }
}

/// The reader's error of a record of the header file `path`: NCBI's exception where NCBI's
/// read of the record fails and its report is ported (`SeqDbError::Ncbi`, which NCBI's
/// applications report as they report any `CException`), else LOSAT's refusal.
fn format_error(error: BerError, path: &Path) -> SeqDbError {
    match error.exception() {
        Some(exception) => SeqDbError::Ncbi(exception),
        None => SeqDbError::Format {
            path: path.to_path_buf(),
            message: error.to_string(),
        },
    }
}

impl HeaderError {
    /// The error of the ordinal-id decoders for `record`. Where the record is no valid
    /// encoding, the error is that of NCBI's read of the whole record (`decode_full_set`, the
    /// one port of it: NCBI reads every member, `x_GetHdrAsn1`), which fails wherever these
    /// decoders' checks do; only a check of their own that NCBI's read passes (an ordinal id
    /// that is negative) keeps their text.
    #[cold]
    #[inline(never)]
    fn into_seqdb(self, record: &[u8], path: &Path) -> SeqDbError {
        match self {
            HeaderError::Ber(error) => match crate::seqid::ber::decode_full_set(record) {
                Err(full) => format_error(full, path),
                Ok(_) => format_error(error, path),
            },
            HeaderError::ParsedSeqIds => SeqDbError::Unsupported {
                path: path.to_path_buf(),
                message: PARSED_SEQIDS.to_string(),
            },
        }
    }
}

pub(super) const PARSED_SEQIDS: &str =
    "the database stores parsed Seq-ids (makeblastdb -parse_seqids, or a \
     database distributed by NCBI), which the decoder of the ordinal ids does not read";

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2208-2218
/// Decodes the `Blast-def-line-set` of the sequence `oid` (volume-relative). An empty record
/// gives no definition line, as in NCBI.
#[cfg(test)]
pub fn decode_defline_set(
    record: &[u8],
    path: &Path,
    oid: u32,
) -> Result<Vec<Defline>, SeqDbError> {
    decode_deflines_as(record, path, oid, false)
}

/// The `taxid` of each definition line of one sequence's record (0 for none): the decode of the
/// definition lines with the taxid (`decode_deflines_as`, as `Volume::deflines` reads them), of
/// which it keeps the taxids, for the reads that need nothing else (the XML2 and JSON `Hit`
/// description, `Volume::defline_taxids`, and the taxonomy pre-pass of their reports).
pub fn decode_defline_taxids(record: &[u8], path: &Path, oid: u32) -> Result<Vec<i32>, SeqDbError> {
    Ok(decode_deflines_as(record, path, oid, true)?
        .into_iter()
        .map(|defline| defline.taxid)
        .collect())
}

/// The definition lines of the record of the sequence `oid` (volume-relative), with their
/// `taxid` when `want_taxid`. A `gnl|BL_ORD_ID|n` with n other than `oid` is a parsed Seq-id
/// (see `decode_defline`).
pub fn decode_deflines_as(
    record: &[u8],
    path: &Path,
    oid: u32,
    want_taxid: bool,
) -> Result<Vec<Defline>, SeqDbError> {
    decode_ordinal_deflines(record, oid, want_taxid).map_err(|error| error.into_seqdb(record, path))
}

/// Decodes the `Blast-def-line-set` of one sequence with every member and every Seq-id
/// (`seqid::ber::decode_full_set`); the `general` ids `BL_ORD_ID` keep their
/// volume-relative tag.
pub fn decode_full(record: &[u8], path: &Path) -> Result<DeflineSet, SeqDbError> {
    crate::seqid::ber::decode_full_set(record).map_err(|error| format_error(error, path))
}

/// [`decode_full`] written over `set`, whose buffers are kept
/// (`seqid::ber::decode_full_set_in`).
pub fn decode_full_in(record: &[u8], path: &Path, set: &mut DeflineSet) -> Result<(), SeqDbError> {
    crate::seqid::ber::decode_full_set_in(record, set).map_err(|error| format_error(error, path))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqinfosrc_aux.cpp:51-69
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1579-1586
// The Seq-align builder takes the best of the ids of every definition line (`GetSeqIDs`) and
// makes no Seq-align of a sequence whose lines hold none. makeblastdb writes no line without a
// Seq-id; an edited header can hold one (line RG D1c-1c-b, `SURVEY_D1c1cb.md` 1.1).
/// Whether the `Blast-def-line-set` `record` holds a Seq-id in any of its definition lines: read
/// up to the first `seqid` list that is not empty, every other member skipped by its length (no
/// Seq-id or title is decoded). An empty record counts as holding one (no read of NCBI's was seen
/// on it; LOSAT reads it as before). The reader's error of a record that is no valid encoding up
/// to there.
pub fn holds_a_seqid(record: &[u8], path: &Path) -> Result<bool, SeqDbError> {
    if record.is_empty() {
        return Ok(true);
    }
    seqid_in_set(record).map_err(|error| format_error(error, path))
}

fn seqid_in_set(record: &[u8]) -> Result<bool, BerError> {
    let mut ber = Reader::new(record);
    ber.expect_tag(TAG_SEQUENCE, "Blast-def-line-set")?;
    let set_end = ber.open()?;
    while !ber.at_end(set_end)? {
        ber.expect_tag(TAG_SEQUENCE, "Blast-def-line")?;
        let defline_end = ber.open()?;
        while !ber.at_end(defline_end)? {
            let tag = ber.tag()?;
            if tag.class != TAG_CONTEXT_CONSTRUCTED || tag.number > u32::from(DEFLINE_LAST_MEMBER) {
                return Err(ber.error(format!(
                    "tag {:#04x}/{} is not a member of Blast-def-line",
                    tag.class, tag.number
                )));
            }
            let member_end = ber.open()?;
            if tag.number == u32::from(DEFLINE_SEQID) {
                ber.expect_tag(TAG_SEQUENCE, "the Seq-id list")?;
                let list_end = ber.open()?;
                if !ber.at_end(list_end)? {
                    return Ok(true);
                }
            }
            while !ber.at_end(member_end)? {
                ber.skip_value(0)?;
            }
        }
    }
    Ok(false)
}

fn decode_ordinal_deflines(
    record: &[u8],
    oid: u32,
    want_taxid: bool,
) -> Result<Vec<Defline>, HeaderError> {
    let mut deflines = Vec::new();
    if record.is_empty() {
        return Ok(deflines);
    }
    let mut ber = Reader::new(record);
    ber.expect_tag(TAG_SEQUENCE, "Blast-def-line-set")?;
    let set_end = ber.open()?;
    while !ber.at_end(set_end)? {
        deflines.push(decode_defline(&mut ber, oid, want_taxid)?);
    }
    // Bytes after the set are not read, as in NCBI (`CObjectIStream::EndOfRead`,
    // `seqid::ber::decode_full_set`).
    if ber.read_utf8() {
        mark_utf8(&mut deflines);
    }
    Ok(deflines)
}

/// Marks the definition lines of a record that holds a UTF8String (`Defline::utf8`).
#[cold]
#[inline(never)]
fn mark_utf8(deflines: &mut [Defline]) {
    for defline in deflines {
        defline.utf8 = true;
    }
}

fn decode_defline(
    ber: &mut Reader<'_>,
    oid: u32,
    want_taxid: bool,
) -> Result<Defline, HeaderError> {
    ber.expect_tag(TAG_SEQUENCE, "Blast-def-line")?;
    let defline_end = ber.open()?;
    // Made after the ids: a record with parsed Seq-ids gets no title here.
    let mut stored_title: &[u8] = &[];
    let mut volume_ordinal = None;
    let mut taxid = 0;
    let mut next_member = DEFLINE_TITLE;
    while !ber.at_end(defline_end)? {
        let tag = ber.tag()?;
        if tag.class != TAG_CONTEXT_CONSTRUCTED
            || tag.number < u32::from(next_member)
            || tag.number > u32::from(DEFLINE_LAST_MEMBER)
        {
            return Err(ber
                .error(format!(
                    "tag {:#04x}/{} is not a member of Blast-def-line",
                    tag.class, tag.number
                ))
                .into());
        }
        // Below `DEFLINE_LAST_MEMBER`.
        let member = tag.number as u8;
        next_member = member + 1;
        let member_end = ber.open()?;
        match member {
            DEFLINE_TITLE => {
                stored_title = ber.visible_string()?;
                ber.close(member_end, "the title")?;
            }
            DEFLINE_SEQID => {
                ber.expect_tag(TAG_SEQUENCE, "the Seq-id list")?;
                let list_end = ber.open()?;
                while !ber.at_end(list_end)? {
                    if volume_ordinal.is_some() {
                        return Err(HeaderError::ParsedSeqIds);
                    }
                    volume_ordinal = Some(decode_ordinal_id(ber)?);
                }
                ber.close(member_end, "the Seq-id list")?;
            }
            // With `want_taxid` the record is read as NCBI's `x_GetHdrAsn1` reads it (the whole
            // `Blast-def-line-set`, seqdbvol.cpp:2208-2218): the taxid as a `TTaxId` (an `int`),
            // memberships, links and other-info as lists of `int` (not kept; `Blast_def_line_.hpp`,
            // generated from blastdb.asn:32-47), so that a value NCBI cannot read fails here.
            DEFLINE_TAXID if want_taxid => {
                taxid = ber.int()?;
                ber.close(member_end, "the taxid")?;
            }
            _ if want_taxid => {
                let name = DEFLINE_MEMBERS[usize::from(member)];
                int_list(ber, name)?;
                ber.close(member_end, name)?;
            }
            // Without it (the check of a volume's first header, `Volume::stored_deflines`, of a
            // sequence that NCBI may never read): taxid, memberships, links, other-info skipped.
            _ => {
                while !ber.at_end(member_end)? {
                    ber.skip_value(0)?;
                }
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:903-908
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:866-872
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2712-2715
    // Without `-parse_seqids` makeblastdb gives each sequence the id `gnl|BL_ORD_ID|<its number in
    // the volume>`, cooked again when the sequence goes to a new volume (writedb_impl.cpp:1214-1215).
    // With `-parse_seqids` it stores a FASTA id `gnl|BL_ORD_ID|n` as written. NCBI reads n plus the
    // first ordinal id of the volume (`SeqDb::defline_set`, seqdbvol.cpp:2222-2241) and its
    // formatters look the sequence up by that number (`seqdb::lookup::SeqidLookup::seqid_to_oid`):
    // another n than the sequence's own finds no sequence or another one (NCBI 2.17.0 writes
    // `Unknown` or another sequence's names, or stops; probes of line DBID N2a,
    // `blastdb-ids-20261021/R/ordid/p1`). Such an id is read as a parsed Seq-id, whose reports LOSAT
    // checks (`report::db_subject`).
    let Some(volume_ordinal) = volume_ordinal.filter(|&ordinal| ordinal == oid) else {
        return Err(HeaderError::ParsedSeqIds);
    };
    let mut title = Vec::with_capacity(stored_title.len());
    let replaced_title_bytes = fix_visible_chars(stored_title, &mut title);
    Ok(Defline {
        title,
        volume_ordinal,
        replaced_title_bytes,
        taxid,
        utf8: false,
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2231-2236
/// Decodes one `Seq-id`, which must be `general { db "BL_ORD_ID", tag id n }`.
fn decode_ordinal_id(ber: &mut Reader<'_>) -> Result<u32, HeaderError> {
    if !ber.take_tag(TAG_CONTEXT_CONSTRUCTED | SEQID_GENERAL)? {
        return Err(HeaderError::ParsedSeqIds);
    }
    let choice_end = ber.open()?;
    ber.expect_tag(TAG_SEQUENCE, "Dbtag")?;
    let dbtag_end = ber.open()?;

    let db_end = ber.open_context(DBTAG_DB, "Dbtag.db")?;
    let db = ber.visible_string()?;
    ber.close(db_end, "Dbtag.db")?;
    if db != b"BL_ORD_ID" {
        return Err(HeaderError::ParsedSeqIds);
    }

    let tag_end = ber.open_context(DBTAG_TAG, "Dbtag.tag")?;
    if !ber.take_tag(TAG_CONTEXT_CONSTRUCTED | OBJECT_ID_ID)? {
        return Err(HeaderError::ParsedSeqIds);
    }
    let id_end = ber.open()?;
    // `Object-id.id` is an `int` (`CObject_id::TId`).
    let id = ber.int()?;
    ber.close(id_end, "Object-id.id")?;
    ber.close(tag_end, "Dbtag.tag")?;
    ber.close(dbtag_end, "Dbtag")?;
    ber.close(choice_end, "Seq-id")?;

    u32::try_from(id).map_err(|_| {
        ber.error(format!("BL_ORD_ID {id} is not a valid ordinal id"))
            .into()
    })
}
