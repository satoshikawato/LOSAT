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
    fix_visible_chars, replaced_bytes, taxid_value, BerError, Reader, DEFLINE_LAST_MEMBER,
    DEFLINE_SEQID, DEFLINE_TAXID, DEFLINE_TITLE, TAG_CONTEXT_CONSTRUCTED, TAG_NUMBER_MASK,
    TAG_SEQUENCE,
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
    /// `taxid` of the definition line, 0 unless `decode_defline_taxids` read it (a value
    /// above `i32::MAX`, which makeblastdb does not write, reads as `i32::MAX`).
    pub taxid: i32,
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

/// The reader's error of a record of the header file `path`.
fn format_error(error: BerError, path: &Path) -> SeqDbError {
    SeqDbError::Format {
        path: path.to_path_buf(),
        message: error.to_string(),
    }
}

impl HeaderError {
    fn into_seqdb(self, path: &Path) -> SeqDbError {
        match self {
            HeaderError::Ber(error) => format_error(error, path),
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
// ```c++
//     CTempString raw_data = x_GetHdrAsn1Binary(oid);
//
//     if (! raw_data.size()) {
//         return bdls;
//     }
//
//     // Now create an ASN.1 object from the memory chunk provided here.
//
//     CObjectIStreamAsnBinary inpstr(raw_data.data(), raw_data.size());
//
//     bdls.Reset(new objects::CBlast_def_line_set);
// ```
/// Decodes the `Blast-def-line-set` of the sequence `oid` (volume-relative). An empty record
/// gives no definition line, as in NCBI.
pub fn decode_defline_set(
    record: &[u8],
    path: &Path,
    oid: u32,
) -> Result<Vec<Defline>, SeqDbError> {
    decode_deflines_as(record, path, oid, false)
}

/// The `taxid` of each definition line of one sequence's record (0 for none), which the XML2
/// and JSON formats show (`Hit` description). It is decoded apart from `decode_defline_set`,
/// which skips it, so that no other report reads it.
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
    decode_ordinal_deflines(record, oid, want_taxid).map_err(|error| error.into_seqdb(path))
}

/// Decodes the `Blast-def-line-set` of one sequence with every member and every Seq-id
/// (`seqid::ber::decode_full_set`); the `general` ids `BL_ORD_ID` keep their
/// volume-relative tag.
pub fn decode_full(record: &[u8], path: &Path) -> Result<DeflineSet, SeqDbError> {
    crate::seqid::ber::decode_full_set(record).map_err(|error| format_error(error, path))
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
    if !ber.at_record_end() {
        return Err(ber.error("unexpected data after Blast-def-line-set").into());
    }
    Ok(deflines)
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
        let member = tag & TAG_NUMBER_MASK;
        if tag & !TAG_NUMBER_MASK != TAG_CONTEXT_CONSTRUCTED
            || member < next_member
            || member > DEFLINE_LAST_MEMBER
        {
            return Err(ber
                .error(format!("tag 0x{tag:02x} is not a member of Blast-def-line"))
                .into());
        }
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
            // The taxid is read only when asked for (`decode_defline_taxids`).
            DEFLINE_TAXID if want_taxid => {
                taxid = taxid_value(ber.integer()?);
                ber.close(member_end, "the taxid")?;
            }
            // taxid, memberships, links, other-info: not used by LOSAT.
            _ => {
                while !ber.at_end(member_end)? {
                    ber.skip_value(0)?;
                }
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:903-908
    // ```c++
    // void CWriteDB_Impl::x_CookHeader()
    // {
    //     int OID = -1;
    //     if (! m_ParseIDs) {
    //         OID = (m_Volume ) ? m_Volume->GetOID() : 0;
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:866-872
    // ```c++
    //     if (OID>=0) {
    //         // Re-inject the BL_ORD_ID
    //         CRef<CSeq_id> gnl_id(new CSeq_id);
    //         gnl_id->SetGeneral().SetDb("BL_ORD_ID");
    //         gnl_id->SetGeneral().SetTag().SetId(OID);
    //         CRef<CBlast_def_line_set> bdls = s_EditDeflineSet(deflines);
    //         bdls->Set().front()->SetSeqid().front() = gnl_id;
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2712-2715
    // ```c++
    //     case eOID:
    //         // Converted to OID directly.
    //         oids.push_back((int) ident);
    //         break;
    // ```
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
    Ok(Defline {
        title: fix_visible_chars(stored_title),
        volume_ordinal,
        replaced_title_bytes: replaced_bytes(stored_title),
        taxid,
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2231-2236
// ```c++
//                 if (seqid.Which() == CSeq_id::e_General) {
//                     CDbtag & dbt = seqid.SetGeneral();
//
//                     if (dbt.GetDb() == "BL_ORD_ID") {
//                         int vol_oid = dbt.GetTag().GetId();
//                         dbt.SetTag().SetId(m_VolStart + vol_oid);
// ```
/// Decodes one `Seq-id`, which must be `general { db "BL_ORD_ID", tag id n }`.
fn decode_ordinal_id(ber: &mut Reader<'_>) -> Result<u32, HeaderError> {
    let tag = ber.tag()?;
    if tag != TAG_CONTEXT_CONSTRUCTED | SEQID_GENERAL {
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
    let object_id = ber.tag()?;
    if object_id != TAG_CONTEXT_CONSTRUCTED | OBJECT_ID_ID {
        return Err(HeaderError::ParsedSeqIds);
    }
    let id_end = ber.open()?;
    let id = ber.integer()?;
    ber.close(id_end, "Object-id.id")?;
    ber.close(tag_end, "Dbtag.tag")?;
    ber.close(dbtag_end, "Dbtag")?;
    ber.close(choice_end, "Seq-id")?;

    match u32::try_from(id) {
        Ok(ordinal) if ordinal <= i32::MAX as u32 => Ok(ordinal),
        _ => Err(ber
            .error(format!("BL_ORD_ID {id} is not a valid ordinal id"))
            .into()),
    }
}
