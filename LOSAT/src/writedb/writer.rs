//! The database being written (`CWriteDB_Impl`,
//! `c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp`): sequences are
//! published one behind the caller, volumes are cut when a file would pass
//! the size limit, and closing renames a single volume or writes the alias
//! file of several.

use std::collections::BTreeSet;

use super::bioseq::Bioseq;
use super::column::{unsupported_late_column, BlobPair, ColumnDates, ColumnMeta};
use super::convert::{
    find_sequence_length, raw_sequence_hash, seq_data_to_binary, sequence_hash, BinarySequence,
};
use super::cout::Cout;
pub use super::files::DbVersion;
use super::files::OffsetArrays;
use super::header::{
    defline_set, encode_defline_set, encode_defline_set_checked, BadChar, BadTitleChar,
    NullSeqIdWriteError,
};
use super::lmdb::{
    lmdb_file_name, lmdb_map_size_from_env, LmdbEnvSettings, WriteDbLmdb, WriteDbTaxId,
};
use super::mask_data::{MaskWriter, MaskedRanges};
use super::ofstream::{OfStream, StreamBuf};
use super::platform::Platform;
pub use super::volume::Indices;
use super::volume::{IndexScratch, Volume};
use super::WriteDbError;
use crate::seqdb::column::MASK_DATA_TITLE;
use crate::seqid::model::{Dbtag, ObjectId};
use crate::seqid::{DeflineSet, ParsedSeqId};

/// `CWriteDB_Impl` of format version 4 or 5, with the indices of `-parse_seqids` and
/// `-hash_index`.
pub struct WriteDb {
    dbname: Vec<u8>,
    protein: bool,
    title: Vec<u8>,
    date: Vec<u8>,
    max_file_size: u64,
    /// `m_DbVersion`.
    version: DbVersion,
    /// `m_VolumeList`; the last element is `m_Volume`.
    volumes: Vec<Volume>,
    closed: bool,
    /// `m_Bioseq`, or `m_Sequence` and `m_Ambig` of `AddSequence` of `seq` and `ambig`, while
    /// `m_HaveSequence` is set.
    pending: Option<Pending>,
    /// What NCBI's serializer posts to the diagnostic stream, in order.
    diagnostics: Vec<u8>,
    /// Application prefix for diagnostic posts (`SetDiagPostPrefix`: `Some("makeblastdb")` in
    /// `makeblastdb`, `None` in `makeprofiledb`).
    diag_prefix: Option<&'static str>,
    /// What NCBI's writer writes to the standard output (`cout`), in order: the lines of the
    /// conversion of NCBI4na data whose size is not the one its length needs
    /// (`convert::ncbi4na_to_binary`).
    stdout: Vec<u8>,
    /// `m_Lmdbdb` (version 5, from the first published sequence until `Close`).
    lmdb: Option<WriteDbLmdb>,
    /// `m_Taxdb` (as `m_Lmdbdb`).
    taxdb: Option<WriteDbTaxId>,
    /// `m_LmdbOid`: the OID of the next sequence over all volumes.
    lmdb_oid: u32,
    /// The platform whose makeblastdb the files follow.
    platform: Platform,
    /// `m_Indices`; `indices.full` is `m_ParseIDs` (`-parse_seqids`).
    indices: Indices,
    /// `m_OidMasks` (`-oid_masks`): the OID mask files each volume writes.
    oid_masks: u64,
    /// `m_Deflines` of the kept sequence (`SetDeflines`: makeblastdb's with `-parse_seqids`,
    /// makeprofiledb's with or without its `-index`).
    pending_deflines: Option<DeflineSet>,
    /// Positions `(defline_index, slot_index)` of `NULL` `Seq-id` elements inside
    /// `pending_deflines` (from binary definition lines `ASN1_BlastDefLine`).
    pending_null_seqids: Vec<(usize, usize)>,
    /// The taxid of the kept sequence's definition line without `-parse_seqids`
    /// (`FixTaxId` of `CBuildDatabase::x_EditHeaders`, `set_taxid`).
    pending_taxid: i32,
    /// `m_Hash` of the kept sequence (`-hash_index`).
    pending_hash: u32,
    /// `m_Pig` of the kept sequence (`SetPig`).
    pending_pig: i32,
    /// `m_Ids` of the sequence published last: the ids of all its definition lines
    /// (`x_CookIds`), which the volume checks for duplicates and indexes, and which the
    /// version 5 LMDB maps to its OID (makedb-v5 S3c).
    ids: Vec<ParsedSeqId>,
    /// `m_TaxIds` of the sequence published last: the taxids of its definition lines,
    /// ascending, each once.
    tax_ids: Vec<i32>,
    /// Every taxid given to `m_Taxdb` (version 5): `GetDBTaxIds` of the new database.
    db_tax_ids: BTreeSet<i32>,
    /// Version 5 with `-parse_seqids`: whether a sequence gave `m_Lmdbdb` keys, and the first
    /// sequence after such a one whose ids were all GIs (its OID and the FASTA form of its
    /// first id), until a sequence with keys follows.
    lmdb_had_keys: bool,
    gi_only_after_keys: Option<(u32, Vec<u8>)>,
    /// The first sequence of GIs only between two sequences with keys: makeblastdb 2.17.0
    /// aborts on it when it writes the `.pos` file (`WriteDb::close`).
    gi_only_between_keys: Option<(u32, Vec<u8>)>,
    /// The mask data: the registered algorithms, the mask column and the blobs of the
    /// sequence being written. Empty for a database without mask data.
    mask: MaskWriter,
    /// `m_ColumnTitles`, `m_ColumnMetas`, `m_HaveBlob` and `m_Blobs`: the columns in the order
    /// of their ids.
    columns: Vec<ColumnSlot>,
    /// The creation dates of the columns' index files.
    column_dates: ColumnDates,
}

/// The sequence kept until the next one is added (`x_Publish`).
#[derive(Debug)]
#[allow(
    clippy::large_enum_variant,
    reason = "one value per writer; a box would add an allocation for each FASTA sequence"
)]
enum Pending {
    /// `AddSequence` of a `CBioseq`.
    Bioseq(Bioseq),
    /// `AddSequence` of `seq` and `ambig`: a sequence in the form of the sequence file.
    Raw(BinarySequence),
}

/// One column of the database.
#[derive(Debug)]
enum ColumnSlot {
    /// The mask column (`m_MaskDataColumn`): its metadata and blobs are the mask writer's.
    Mask,
    /// A column of `CreateUserColumn`: its title, its metadata, and the blobs of the sequence
    /// being written (`SetBlobData`; `have` is `m_HaveBlob`), of which a volume writes the first.
    User {
        title: Vec<u8>,
        meta: ColumnMeta,
        blobs: [Vec<u8>; 2],
        have: u8,
    },
}

/// "Error: provided column ID is not valid" (`AddColumnMetaData`, `SetBlobData`).
fn invalid_column_id() -> WriteDbError {
    WriteDbError::Creation("Error: provided column ID is not valid".to_string())
}

impl std::fmt::Debug for WriteDb {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("WriteDb")
            .field("dbname", &String::from_utf8_lossy(&self.dbname))
            .field("protein", &self.protein)
            .field("version", &self.version)
            .field("volumes", &self.volumes.len())
            .field("closed", &self.closed)
            .field("lmdb_oid", &self.lmdb_oid)
            .field("platform", &self.platform)
            .field("indices", &self.indices)
            .finish_non_exhaustive()
    }
}

// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1864-1885
// makeblastdb prints, for each byte of a title that a VisibleString cannot
// hold and each time the definition line is serialized (`post_bad_chars`: of any
// string of a parsed definition line, under the path of its member):
// "Error: (803.7) [makeblastdb] Blast-def-line-set.E.title", "Bad char [0xC3]
// in string at byte 21" and the title.
fn post_bad_title_chars(
    diagnostics: &mut Vec<u8>,
    diag_prefix: Option<&str>,
    title: &[u8],
    bad: &[BadTitleChar],
) {
    let header = bad_char_header(diag_prefix);
    for bad_char in bad {
        let message = [
            format!(
                "Blast-def-line-set.E.title\nBad char [0x{:X}] in string at byte {}\n",
                bad_char.byte, bad_char.stream_offset
            )
            .as_bytes(),
            title,
        ]
        .concat();
        diagnostics.extend_from_slice(&crate::ncbi_exception::post_text(&header, &message));
        diagnostics.push(b'\n');
    }
}

/// `post_bad_title_chars` for the strings of definition lines with parsed Seq-ids:
/// `Error: (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.local.str`, `Bad char
/// [0xE9] in string at byte 38` and the string (makeblastdb 2.17.0, survey B fixture
/// p23_nonascii; the member paths of the other alternatives from the S3a probes).
fn post_bad_chars(diagnostics: &mut Vec<u8>, diag_prefix: Option<&str>, bad: &[BadChar]) {
    let header = bad_char_header(diag_prefix);
    for bad_char in bad {
        let message = [
            format!(
                "{}\nBad char [0x{:X}] in string at byte {}\n",
                bad_char.path, bad_char.byte, bad_char.stream_offset
            )
            .as_bytes(),
            &bad_char.string,
        ]
        .concat();
        diagnostics.extend_from_slice(&crate::ncbi_exception::post_text(&header, &message));
        diagnostics.push(b'\n');
    }
}

/// The header of the 803.7 posts of `post_bad_title_chars` and `post_bad_chars`.
fn bad_char_header(diag_prefix: Option<&str>) -> Vec<u8> {
    let tag = diag_prefix.map_or_else(String::new, |prefix| format!("[{prefix}] "));
    format!("Error: (803.7) {tag}").into_bytes()
}

// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:1232-1256
// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:373-384,522-557
/// Posts NCBI's `802.5` `CObjectOStream` error for a `NULL` `Seq-id` inside
/// `Blast-def-line-set.E.seqid.E` at `byte` and returns the matching `CUnassignedMember::eWrite`
/// exception.
fn null_seqid_write_error(
    diagnostics: &mut Vec<u8>,
    diag_prefix: Option<&str>,
    byte: usize,
) -> WriteDbError {
    let null = crate::report::serial::encoder::NullElement {
        stream: crate::report::serial::encoder::Stream::AsnBinary,
        element_frame: true,
        trace: "Blast-def-line-set.E.seqid.E".to_string(),
        position: crate::report::serial::encoder::StreamPosition::Byte(byte as u64),
        frames: Vec::new(),
    };
    diagnostics.extend_from_slice(&null.post(diag_prefix));
    WriteDbError::Exception(null.exception())
}

impl WriteDb {
    /// `CWriteDB` of the database name, the sequence type, the title and the database version
    /// (among other arguments), followed by
    /// `SetMaxFileSize`. `date` is the creation date of the index files, in
    /// the form of `current_date`'s argument; `platform` the platform whose
    /// makeblastdb the files follow ([`Platform::host`] for a build);
    /// `indices` is `indexing` (`eFullIndex` for `-parse_seqids`, whose definition lines
    /// are given by `set_deflines` and written as they are; `eAddHash` for
    /// `-hash_index`).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        dbname: &[u8],
        protein: bool,
        title: &[u8],
        date: &[u8],
        max_file_size: u64,
        version: DbVersion,
        platform: Platform,
        indices: Indices,
    ) -> Self {
        Self {
            dbname: dbname.to_vec(),
            protein,
            title: title.to_vec(),
            date: date.to_vec(),
            max_file_size,
            version,
            volumes: Vec::new(),
            closed: false,
            pending: None,
            diagnostics: Vec::new(),
            diag_prefix: Some("makeblastdb"),
            stdout: Vec::new(),
            lmdb: None,
            taxdb: None,
            lmdb_oid: 0,
            platform,
            indices,
            oid_masks: 0,
            pending_deflines: None,
            pending_null_seqids: Vec::new(),
            pending_taxid: 0,
            pending_hash: 0,
            pending_pig: 0,
            ids: Vec::new(),
            tax_ids: Vec::new(),
            db_tax_ids: BTreeSet::new(),
            lmdb_had_keys: false,
            gi_only_after_keys: None,
            gi_only_between_keys: None,
            mask: MaskWriter::default(),
            columns: Vec::new(),
            column_dates: ColumnDates::default(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1563-1566
    /// `SetMaxFileSize`: the limit of the files of the volumes created from
    /// now on.
    pub fn set_max_file_size(&mut self, max_file_size: u64) {
        self.max_file_size = max_file_size;
    }

    /// Sets the application prefix of serializer diagnostics (`None` for `makeprofiledb`, which
    /// calls no `SetDiagPostPrefix`).
    pub(crate) fn set_diag_prefix(&mut self, prefix: Option<&'static str>) {
        self.diag_prefix = prefix;
    }

    /// Takes the diagnostics posted since the last call (text for the
    /// standard error stream).
    pub fn take_diagnostics(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.diagnostics)
    }

    /// Takes the text written to the standard output since the last call. Only a PSSM's query
    /// gives such text (makeprofiledb): makeblastdb's FASTA source gives NCBI4na data of the size
    /// its length needs, and a copied sequence is not converted.
    pub fn take_stdout(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.stdout)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:473-475
    // The call serializes the definition line once (writedb_impl.cpp:887-893,
    // with the sequence's own id), which posts the diagnostics of its title;
    // the bytes are not kept.
    /// `CWriteDB::ExtractBioseqDeflines` as far as it is visible: the
    /// diagnostics of the title.
    pub fn extract_bioseq_deflines(&mut self, bioseq: &Bioseq) {
        let mut bad = Vec::new();
        defline_set(&bioseq.title, 0, 0, &mut bad);
        post_bad_title_chars(&mut self.diagnostics, self.diag_prefix, &bioseq.title, &bad);
    }

    /// The serialization of `ExtractBioseqDeflines` with `-parse_seqids`: the
    /// diagnostics of the definition lines it built (`writedb::deflines`), before
    /// `FixTaxId` gives them a taxid.
    pub fn post_extracted_deflines(&mut self, deflines: &DeflineSet) {
        let mut bad = Vec::new();
        encode_defline_set(deflines, &mut bad);
        post_bad_chars(&mut self.diagnostics, self.diag_prefix, &bad);
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1242-1249
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:560-589
    // `s_CheckEmptyLists` removes empty membership and link lists, which the definition lines
    // of FASTA input do not have (a copied database's may).
    /// `SetDeflines`: the definition lines of the sequence added last.
    pub fn set_deflines(&mut self, mut deflines: DeflineSet) {
        super::deflines::strip_empty_lists(&mut deflines);
        self.pending_deflines = Some(deflines);
        self.pending_null_seqids.clear();
    }

    /// [`WriteDb::set_deflines`], keeping the `(defline_index, slot_index)` positions of any
    /// `NULL` `Seq-id` elements from binary definition lines (`ASN1_BlastDefLine`).
    pub(crate) fn set_deflines_with_null_seqids(
        &mut self,
        deflines: DeflineSet,
        null_seqids: Vec<(usize, usize)>,
    ) {
        self.set_deflines(deflines);
        self.pending_null_seqids = null_seqids;
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1558-1561
    /// `SetPig`: the PIG of the sequence added last (indexed in a protein volume with indices
    /// when it is not 0).
    pub fn set_pig(&mut self, pig: i32) {
        self.pending_pig = pig;
    }

    /// `SetDeflines` without `-parse_seqids`: the definition line is built at the
    /// publication (`cook_header`) with this taxid (`FixTaxId`).
    pub fn set_taxid(&mut self, taxid: i32) {
        self.pending_taxid = taxid;
    }

    /// The taxids of the database (version 5; `GetDBTaxIds` of the new database), each once.
    pub fn db_tax_ids(&self) -> &BTreeSet<i32> {
        &self.db_tax_ids
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:158-180
    /// `AddSequence` of a `CBioseq`: publishes the previous sequence and
    /// keeps this one.
    pub fn add_sequence(&mut self, bioseq: Bioseq) -> Result<(), WriteDbError> {
        self.add_sequence_inner(bioseq, None)
    }

    /// `AddSequence` of a `CBioseq`, flushing any `WriteDB_Ncbi4naToBinary` stdout line to `cout`
    /// inside `x_Publish` (before creating or writing the volume) and stopping with
    /// [`WriteDbError::ClosedPipe`] when `cout` hits a closed pipe where SIGPIPE ends NCBI.
    pub(crate) fn add_sequence_with_cout(
        &mut self,
        bioseq: Bioseq,
        cout: &mut Cout<'_>,
    ) -> Result<(), WriteDbError> {
        self.add_sequence_inner(bioseq, Some(cout))
    }

    fn add_sequence_inner(
        &mut self,
        bioseq: Bioseq,
        cout: Option<&mut Cout<'_>>,
    ) -> Result<(), WriteDbError> {
        self.publish_inner(cout)?;
        self.reset_sequence_data();
        if bioseq.is_aa != self.protein {
            let name = |protein: bool| if protein { "protein" } else { "nucleotide" };
            return Err(WriteDbError::Creation(format!(
                "Invalid molecule type of sequence added ({}); expected {}",
                name(bioseq.is_aa),
                name(self.protein)
            )));
        }
        if self.indices.hash {
            if let Some(seq_data) = &bioseq.seq_data {
                // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1709-1712
                self.pending_hash = match sequence_hash(seq_data, bioseq.length) {
                    Ok(hash) => hash.unwrap_or(0),
                    Err(super::convert::SequenceHashError::OutOfRange(exception)) => {
                        return Err(WriteDbError::Exception(exception));
                    }
                    Err(super::convert::SequenceHashError::PastEnd(coding)) => {
                        return Err(WriteDbError::Unsupported(format!(
                            "the sequence hash (-hash_index) of {coding} data shorter than its \
                             length (NCBI makeblastdb reads past the end of the buffer)"
                        )));
                    }
                };
            }
        }
        self.pending = Some(Pending::Bioseq(bioseq));
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:139-156
    // makeblastdb adds such sequences with `-parse_seqids` indices (`eFullIndex`), and their
    // definition lines come from `SetDeflines`: the header of a sequence without them (built
    // from a `CBioseq` with the injected `gnl|BL_ORD_ID|<oid>`) has no counterpart here.
    /// `AddSequence` of `seq` and `ambig`: publishes the previous sequence and keeps this one,
    /// given in the form of the sequence file (a protein's NCBIstdaa bytes; a nucleotide's packed
    /// bases and its ambiguity data), as a database stores it.
    pub fn add_raw_sequence(&mut self, seq: &[u8], ambig: &[u8]) -> Result<(), WriteDbError> {
        self.publish()?;
        self.reset_sequence_data();
        if !self.indices.full {
            return Err(WriteDbError::Unsupported(
                "a sequence in the form of the sequence file is added with -parse_seqids indices \
                 only (as makeblastdb adds it)"
                    .to_string(),
            ));
        }
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1690-1700
        // NCBI reads the ambiguity data of the hash past the end of its vector where they are
        // malformed, which LOSAT refuses.
        if self.indices.hash {
            self.pending_hash = raw_sequence_hash(self.protein, seq, ambig).map_err(|message| {
                WriteDbError::Unsupported(format!(
                    "the ambiguity data of a copied sequence cannot be read for its hash \
                         (-hash_index): {message}; NCBI makeblastdb reads past them there"
                ))
            })?;
        }
        self.pending = Some(Pending::Raw(BinarySequence {
            seq: seq.to_vec(),
            ambig: ambig.to_vec(),
        }));
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:139-149
    /// `AddSequence` of a raw protein whose signed 32-bit `GetRawSeqAndAmbig` length wrapped
    /// negative (`>= 2^31` residues): publishes the previous sequence, resets the pending sequence
    /// state, and returns the first error (`publish` failure if any, otherwise NCBI's
    /// `std::length_error` `basic_string::_M_create`) without allocating, so `~CWriteDB_Impl`
    /// closes the partial database containing prior sequences.
    #[must_use]
    pub fn add_raw_negative_length(&mut self) -> WriteDbError {
        if let Err(error) = self.publish() {
            return error;
        }
        self.reset_sequence_data();
        WriteDbError::LengthError("basic_string::_M_create")
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:110-137
    /// `x_ResetSequenceData`: the next sequence starts without the data of the previous one
    /// (its definition lines, taxid, hash, PIG, ids and column blobs).
    fn reset_sequence_data(&mut self) {
        self.pending_deflines = None;
        self.pending_null_seqids.clear();
        self.pending_taxid = 0;
        self.pending_hash = 0;
        self.pending_pig = 0;
        self.ids.clear();
        self.tax_ids.clear();
        self.mask.reset_sequence();
        for slot in &mut self.columns {
            if let ColumnSlot::User { blobs, have, .. } = slot {
                blobs[0].clear();
                blobs[1].clear();
                *have = 0;
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:903-908
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:779-781,887-893,925-952
    // With `-parse_seqids` (OID -1) the definition lines given by `SetDeflines` are
    // serialized once; a second `x_CookHeader` for a new volume returns at once.
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:895-900
    // `GetTaxIds` is the definition line's taxid (set on every one by `FixTaxId`) and its
    // `links` (the leaf taxids of a copied database; FASTA input has none).
    /// `x_CookHeader` and `x_CookIds` with `-parse_seqids`.
    fn cook_parsed_header(&mut self) -> Result<Vec<u8>, WriteDbError> {
        let deflines = self.pending_deflines.take().unwrap_or_default();
        let null_seqids = std::mem::take(&mut self.pending_null_seqids);
        let mut bad = Vec::new();
        let encoded = encode_defline_set_checked(&deflines, &null_seqids, &mut bad);
        post_bad_chars(&mut self.diagnostics, self.diag_prefix, &bad);
        let binhdr = match encoded {
            Ok(binhdr) => binhdr,
            Err(NullSeqIdWriteError { byte }) => {
                return Err(null_seqid_write_error(
                    &mut self.diagnostics,
                    self.diag_prefix,
                    byte,
                ));
            }
        };
        self.tax_ids = union_of_tax_ids(&deflines.deflines);
        self.ids = deflines
            .deflines
            .into_iter()
            .flat_map(|defline| defline.ids)
            .collect();
        Ok(binhdr)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:925-929
    // Without `-parse_seqids` the definition line's first id is the injected
    // `gnl|BL_ORD_ID|<oid>` (writedb_impl.cpp:866-874), so `m_Ids` is that id. When the
    // sequence goes to a new volume, the header is cooked again with OID 0 but `m_Ids`
    // keeps the id of the first cooking, the OID the sequence would have had in the full
    // volume: makeblastdb 2.17.0 indexes `gnl|bl_ord_id|<that oid>` for OID 0 of the new
    // volume and checks it for duplicates there (s3b/probe: `-hash_index -max_file_sz 4000`
    // on 3 long and 20 short proteins stops with "Duplicate seq_ids are found:
    // GNL|BL_ORD_ID|6"). The ids are needed by a volume with indices and by the OID mask of
    // `-oid_masks` (`writedb::oid_list::is_model`).
    /// `x_CookHeader`: the header record with the number the sequence gets
    /// in the current volume; the first cooking of a sequence also sets `m_Ids` (with
    /// indices or OID masks) and `m_TaxIds`.
    fn cook_header(&mut self, title: &[u8]) -> Result<Vec<u8>, WriteDbError> {
        let oid = self.volumes.last().map_or(0, Volume::oid);
        if self.pending_deflines.is_some() {
            return self.cook_given_header(oid);
        }
        let mut bad = Vec::new();
        let binhdr = defline_set(title, oid, self.pending_taxid, &mut bad);
        post_bad_title_chars(&mut self.diagnostics, self.diag_prefix, title, &bad);
        if self.tax_ids.is_empty() {
            self.tax_ids.push(self.pending_taxid);
        }
        // The OID mask reads the ids too (`m_Ids` of every definition line, makeblastdb 2.17.0
        // probe `p_nops1`: a sequence without `-parse_seqids` is no model).
        if (self.indices.any() || self.oid_masks != 0) && self.ids.is_empty() {
            self.ids.push(ParsedSeqId::General(Dbtag {
                db: b"BL_ORD_ID".to_vec(),
                tag: ObjectId::Id(oid as i32),
            }));
        }
        Ok(binhdr)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:779-781,866-893,925-952
    // Definition lines given by `SetDeflines` without `eFullIndex` (makeprofiledb `-index
    // false`): every cooking takes the given lines with the first id of the first line
    // replaced by `gnl|BL_ORD_ID|<oid>` and serializes them; `m_Ids` (every line's ids after
    // the replacement) and `m_TaxIds` are set by the first cooking only. When the first slot
    // of the first line was `NULL`, replacing it with `gnl|BL_ORD_ID|<oid>` removes that `NULL`
    // element before serialization.
    /// `x_CookHeader` and `x_CookIds` of given definition lines without `-parse_seqids`.
    fn cook_given_header(&mut self, oid: u32) -> Result<Vec<u8>, WriteDbError> {
        let mut deflines = self
            .pending_deflines
            .clone()
            .expect("the caller checked the given definition lines");
        let mut null_seqids = self.pending_null_seqids.clone();
        let gnl = ParsedSeqId::General(Dbtag {
            db: b"BL_ORD_ID".to_vec(),
            tag: ObjectId::Id(oid as i32),
        });
        if let Some(pos) = null_seqids.iter().position(|&pair| pair == (0, 0)) {
            null_seqids.remove(pos);
            if let Some(defline) = deflines.deflines.first_mut() {
                defline.ids.insert(0, gnl);
            }
        } else if let Some(first) = deflines
            .deflines
            .first_mut()
            .and_then(|defline| defline.ids.first_mut())
        {
            *first = gnl;
        }
        let mut bad = Vec::new();
        let encoded = encode_defline_set_checked(&deflines, &null_seqids, &mut bad);
        post_bad_chars(&mut self.diagnostics, self.diag_prefix, &bad);
        let binhdr = match encoded {
            Ok(binhdr) => binhdr,
            Err(NullSeqIdWriteError { byte }) => {
                return Err(null_seqid_write_error(
                    &mut self.diagnostics,
                    self.diag_prefix,
                    byte,
                ));
            }
        };
        if self.tax_ids.is_empty() {
            self.tax_ids = union_of_tax_ids(&deflines.deflines);
        }
        if self.ids.is_empty() {
            self.ids = deflines
                .deflines
                .into_iter()
                .flat_map(|defline| defline.ids)
                .collect();
        }
        Ok(binhdr)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1127-1133
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1096-1098
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1156-1239
    // The size check of the old volume uses the header cooked with the old
    // volume's number; the record written to a new volume is cooked again
    // with number 0.
    //
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1169-1175
    // The same follows the write into a new volume (writedb_impl.cpp:1226-1232).
    // Without `-parse_seqids` no Seq-id goes to `m_Lmdbdb`. `m_TaxIds` is the
    // set of the tax ids of the sequence's definition lines: makeblastdb sets
    // tax id 0 on each of them without `-taxid` and `-taxid_map`
    // (`CTaxIdSet::FixTaxId`, taxid_set.cpp:84-136), and `InsertEntries` records
    // tax id 0 for the set {0} as for an empty set (writedb_lmdb.cpp:446-464).
    /// `x_Publish`: writes the kept sequence, into a new volume when the
    /// current one is full.
    fn publish(&mut self) -> Result<(), WriteDbError> {
        self.publish_inner(None)
    }

    fn publish_inner(&mut self, cout: Option<&mut Cout<'_>>) -> Result<(), WriteDbError> {
        let Some(pending) = self.pending.take() else {
            return Ok(());
        };
        if self.version == DbVersion::V5 && self.lmdb.is_none() {
            self.open_lmdb_files()?;
        }
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:831-837
        // Definition lines given as an empty set (makeprofiledb: binary definition lines that
        // hold none) stop the publication before the header is made.
        if self
            .pending_deflines
            .as_ref()
            .is_some_and(|set| set.deflines.is_empty())
        {
            return Err(WriteDbError::Creation(
                "Error: No deflines provided.".to_string(),
            ));
        }
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:865-874
        let first_has_null = self.pending_null_seqids.iter().any(|&(d, _)| d == 0);
        if !self.indices.full
            && !first_has_null
            && self
                .pending_deflines
                .as_ref()
                .and_then(|set| set.deflines.first())
                .is_some_and(|defline| defline.ids.is_empty())
        {
            return Err(WriteDbError::OrdinalIdWithoutIds);
        }
        let title: &[u8] = match &pending {
            Pending::Bioseq(bioseq) => &bioseq.title,
            Pending::Raw(_) => &[],
        };
        let mut binhdr = if self.indices.full {
            self.cook_parsed_header()?
        } else {
            self.cook_header(title)?
        };
        let cooked;
        let binary = match &pending {
            Pending::Bioseq(bioseq) => {
                cooked = cook_sequence(bioseq, &mut self.stdout)?;
                // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_convert.cpp:464-466
                // `WriteDB_Ncbi4naToBinary` writes `ib=...,n4sz=...` to `cout` with `endl` before
                // the volume is opened or written: on a closed pipe, SIGPIPE ends NCBI right here.
                if let Some(cout) = cout {
                    let text = std::mem::take(&mut self.stdout);
                    if !text.is_empty() {
                        cout.insert(&text);
                        cout.flush();
                        if cout.rdbuf().closed_pipe() {
                            return Err(WriteDbError::ClosedPipe);
                        }
                    }
                }
                &cooked
            }
            Pending::Raw(binary) => binary,
        };

        let mut done = false;
        if !self.volumes.is_empty() {
            done = self.write_to_volume(binary, &binhdr)?;
        }
        if !done {
            let index = self.volumes.len() as u32;
            // The new volume's index file takes the arrays of the closed one
            // (`OffsetArrays`), its indices their memory (`IndexScratch`).
            let mut arrays = OffsetArrays::default();
            let mut scratch = IndexScratch::default();
            if let Some(volume) = self.volumes.last_mut() {
                volume.close()?;
                arrays = volume.take_offset_arrays();
                scratch = volume.take_index_scratch();
            }
            // The dates of the columns are read before the volume's files are created, so
            // that a clock that cannot be read leaves none.
            let column_dates = (0..self.columns.len())
                .map(|_| self.column_dates.date(index as usize))
                .collect::<Result<Vec<_>, _>>()?;
            self.volumes.push(Volume::create(
                &self.dbname,
                self.protein,
                &self.title,
                &self.date,
                index,
                self.max_file_size,
                self.version,
                self.platform.sizeof_long,
                arrays,
                self.indices,
                scratch,
                self.oid_masks,
            )?);
            let volume = self.volumes.last_mut().expect("a volume was pushed");
            for (slot, date) in self.columns.iter().zip(&column_dates) {
                match slot {
                    ColumnSlot::Mask => {
                        let meta = self.mask.column_meta().cloned().unwrap_or_default();
                        volume.create_column(MASK_DATA_TITLE, &meta, date, self.max_file_size)?;
                    }
                    ColumnSlot::User { title, meta, .. } => {
                        volume.create_column(title, meta, date, self.max_file_size)?;
                    }
                }
            }
            if !self.indices.full {
                binhdr = self.cook_header(title)?;
            }
            done = self.write_to_volume(binary, &binhdr)?;
            if !done {
                return Err(WriteDbError::Creation(
                    "Cannot write sequence to volume.".to_string(),
                ));
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1158-1175
    /// `WriteSequence` of `m_Volume` with the kept sequence's data, then its LMDB entries when
    /// it was written; `false` when it does not fit the current volume.
    fn write_to_volume(
        &mut self,
        binary: &BinarySequence,
        binhdr: &[u8],
    ) -> Result<bool, WriteDbError> {
        let blobs: Vec<BlobPair<'_>> = self
            .columns
            .iter()
            .map(|slot| match slot {
                ColumnSlot::Mask => self.mask.blobs(),
                ColumnSlot::User { blobs, .. } => BlobPair {
                    blob: &blobs[0],
                    blob2: &blobs[1],
                },
            })
            .collect();
        let mask_column = self
            .columns
            .iter()
            .position(|slot| matches!(slot, ColumnSlot::Mask));
        let volume = self.volumes.last_mut().expect("a volume is open");
        let done = volume.write_sequence(
            &binary.seq,
            &binary.ambig,
            binhdr,
            &self.ids,
            self.pending_pig,
            self.pending_hash,
            &blobs,
            mask_column,
            &mut self.diagnostics,
        )?;
        drop(blobs);
        if done {
            self.insert_lmdb_entries()?;
        }
        Ok(done)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.hpp:655
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdbcommon.hpp:58
    // An OID is a 4-byte signed integer in NCBI's writer and in its LMDB files, so a database
    // has at most 2^31 sequences; NCBI's `m_LmdbOid` overflows when it is incremented past that.
    // LOSAT stops with an error at the sequence that would get OID 2^31 (no NCBI overflow
    // reproduced).
    /// The keys of the ids (`-parse_seqids`) and the tax ids of the sequence just written,
    /// under the next OID of the database (version 5).
    fn insert_lmdb_entries(&mut self) -> Result<(), WriteDbError> {
        if let Some(taxdb) = self.taxdb.as_mut() {
            if self.lmdb_oid > i32::MAX as u32 {
                return Err(WriteDbError::Other(format!(
                    "a version 5 database holds at most {} sequences (an OID is an Int4)",
                    i32::MAX as u64 + 1
                )));
            }
            // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1169-1173
            if self.indices.full {
                if let Some(lmdb) = self.lmdb.as_mut() {
                    let keys = lmdb.insert_entries(&self.ids, self.lmdb_oid)?;
                    // A sequence of GIs only between two sequences with keys stops
                    // makeblastdb's `.pos` file (`close`): the first one is kept to name it.
                    if keys > 0 {
                        if self.gi_only_between_keys.is_none() {
                            self.gi_only_between_keys = self.gi_only_after_keys.take();
                        }
                        self.gi_only_after_keys = None;
                        self.lmdb_had_keys = true;
                    } else if self.lmdb_had_keys && self.gi_only_after_keys.is_none() {
                        let mut name = Vec::new();
                        if let Some(id) = self.ids.first() {
                            crate::seqid::label::fasta(id, &mut name);
                        }
                        self.gi_only_after_keys = Some((self.lmdb_oid, name));
                    }
                }
            }
            taxdb.insert_entries(&self.tax_ids, self.lmdb_oid);
            if self.tax_ids.is_empty() {
                self.db_tax_ids.insert(0);
            } else {
                self.db_tax_ids.extend(self.tax_ids.iter().copied());
            }
            self.lmdb_oid += 1;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1136-1154
    // The page size of LMDB, its default map sizes and the length of its files on Windows are
    // the platform's (DECISIONS 8 of the makedb-v5 line, `Platform`).
    /// Creates the LMDB files of a version 5 database (`<db>.pdb` and
    /// `<db>.ptf`, or `.ndb` and `.ntf`) at its first published sequence.
    fn open_lmdb_files(&mut self) -> Result<(), WriteDbError> {
        let lmdb_fname_w_path = build_lmdb_file_name(&self.dbname, self.protein)?;
        let map_size = lmdb_map_size_from_env()?;
        let settings = LmdbEnvSettings::from_env();
        let platform = self.platform.lmdb()?;
        let (lmdb_map, tax_map) = if map_size > 0 {
            (map_size, map_size)
        } else {
            (
                self.platform.default_lmdb_map_size,
                self.platform.default_taxid_map_size,
            )
        };
        let lmdb = WriteDbLmdb::new(
            &crate::utils::os_bytes::path_of_bytes(&lmdb_fname_w_path),
            lmdb_map,
            &platform,
            &settings,
        )?;
        self.lmdb = Some(lmdb);
        let taxdb = WriteDbTaxId::new(
            &lmdb_file_name(
                &crate::utils::os_bytes::path_of_bytes(&lmdb_fname_w_path),
                "tf",
            ),
            tax_map,
            &platform,
            &settings,
        )?;
        self.taxdb = Some(taxdb);
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:326-354
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:378-394
    // The volume names are the final ones (after `RenameSingle` or
    // `RenameFileIndex`). `Reset` runs the destructors of the two LMDB
    // writers, which write their lookup files and their LMDB transactions
    // (`WriteDbLmdb::close`, `WriteDbTaxId::close`). An exception that leaves `Close` before
    // them (the duplicate id of the last sequence, found when `x_Publish` writes it) leaves
    // the last volume open and both writers alive: their destructors run when the writer is
    // destroyed, after `x_EndBuild` has removed the files of the build
    // (`destroy_after_failed_close`). An exception that leaves one of these destructors ends
    // makeblastdb (`std::terminate`): LOSAT ends with an error there, and no destructor runs
    // after it (`lmdb` and `taxdb` are dropped without their work).
    //
    // The GI masks close between the last volume and the renames, and the alias file follows
    // them, before the LMDB step: a version 5 build whose LMDB step fails (duplicate Seq-ids)
    // leaves its GI-mask files and its alias file, which `list_files` does not name. When
    // `x_Publish` fails, the GI masks are not closed: their data files hold what was written
    // (the streams write it when they are dropped, as NCBI's destructors do), and no offset or
    // index file is made.
    /// `Close`: publishes the last sequence, closes the last volume and the GI masks, gives
    /// the files their final names, writes the alias file of a database of several volumes or
    /// with GI masks, and the LMDB files of version 5.
    pub fn close(&mut self) -> Result<(), WriteDbError> {
        self.close_inner(None)
    }

    /// `Close`, flushing any `WriteDB_Ncbi4naToBinary` stdout line emitted while publishing the
    /// final pending sequence directly to `cout` so a closed pipe stops before volume close,
    /// rename, or LMDB/taxdb close (as SIGPIPE ends NCBI inside `x_Publish`).
    pub(crate) fn close_with_cout(&mut self, cout: &mut Cout<'_>) -> Result<(), WriteDbError> {
        self.close_inner(Some(cout))
    }

    fn close_inner(&mut self, cout: Option<&mut Cout<'_>>) -> Result<(), WriteDbError> {
        if self.closed {
            return Ok(());
        }
        self.closed = true;
        self.publish_inner(cout)?;
        let count = self.volumes.len();
        if count == 0 {
            return Ok(());
        }
        self.volumes[count - 1].close()?;
        drop(self.volumes[count - 1].take_offset_arrays());
        drop(self.volumes[count - 1].take_index_scratch());
        self.mask.close_gi_masks();
        if count == 1 {
            self.volumes[0].rename_single();
        } else if count > 100 {
            let num_digits = count.ilog10() + 1;
            for volume in &mut self.volumes {
                volume.rename_file_index(num_digits);
            }
        }
        if count > 1 || self.mask.uses_gi_masks() {
            self.make_alias();
        }
        if self.version == DbVersion::V5 {
            if let Some(lmdb) = self.lmdb.as_mut() {
                let vol_names: Vec<&[u8]> = self
                    .volumes
                    .iter()
                    .map(|volume| file_name(volume.volume_name()))
                    .collect();
                let vol_num_oids: Vec<u32> = self.volumes.iter().map(Volume::oid).collect();
                lmdb.insert_volumes_info(&vol_names, &vol_num_oids)?;
                // `m_Lmdbdb.Reset()`: makeblastdb aborts in this destructor (the `.pos` file).
                if let Some((oid, name)) = self.gi_only_between_keys.take() {
                    self.lmdb = None;
                    self.taxdb = None;
                    return Err(gi_only_between_keys_error(oid, &name));
                }
                if let Some(lmdb) = self.lmdb.take() {
                    if let Err(error) = lmdb.close() {
                        self.taxdb = None;
                        return Err(error);
                    }
                }
                if let Some(taxdb) = self.taxdb.take() {
                    taxdb.close()?;
                }
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:100-108
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.hpp:629-648
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:146-151
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:65-72,437-444
    // After a `Close` that threw, the destructor's `Close` returns at once (`m_Closed`), and the
    // members are destroyed in the reverse order of their declaration: `m_Taxdb` (the `.pot`
    // and `.pto` files, the commit into `.ptf`), `m_Lmdbdb` (the `.pos` file, the commit into the
    // `.pdb` file, which `x_EndBuild` removed: it goes to the removed file), then
    // `m_VolumeList` and `m_Volume`, whose last volume closes: its index file is written to the
    // removed one, and the files it makes at the close (the ISAM files, the OID to GI file) are
    // made under the volume's numbered name, since no rename ran. makeblastdb destroys the
    // writer when it ends, after its message; the files are those of these steps either way.
    // A destructor that throws ends makeblastdb (`std::terminate`): the steps after it do not
    // run (the tax id writer of a database whose first sequence failed: "No tax info for any
    // oid"). An error of the last volume's close is LOSAT's 4 GiB stop, after the build's error.
    /// The destructors of the writer's members after a `close` that failed (`x_EndBuild` has
    /// removed the files of the build).
    /// `Err` is the error of the destructor that throws, after which nothing runs.
    pub fn destroy_after_failed_close(&mut self) -> Result<(), WriteDbError> {
        self.destroy_lmdb_writers()?;
        if let Some(volume) = self.volumes.last_mut() {
            let _ = volume.close();
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:324-397
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:437-444,533-536
    // A `Close` that found no volume leaves the LMDB writers that the publication of a first
    // sequence created before it failed; their destructors run with the writer's, and the tax
    // id writer's, with no OID, throws "No tax info for any oid" (makeprofiledb 2.17.0 then ends
    // with `std::terminate`). After any other `Close` they are gone and nothing runs.
    /// The destructors of the writer's members after a `close` that succeeded; `Err` is the error
    /// of the destructor that throws.
    pub fn destroy_after_close(&mut self) -> Result<(), WriteDbError> {
        self.destroy_lmdb_writers()
    }

    /// The destructors of the tax id writer, then of the LMDB writer, when they are alive; `Err`
    /// is the error of the one that throws (the other's work is not done).
    fn destroy_lmdb_writers(&mut self) -> Result<(), WriteDbError> {
        if let Some(taxdb) = self.taxdb.take() {
            if let Err(error) = taxdb.close() {
                self.lmdb = None;
                return Err(error);
            }
        }
        if let Some(lmdb) = self.lmdb.take() {
            lmdb.close()?;
        }
        Ok(())
    }

    /// `RegisterMaskAlgorithm` of a program, its options and a name: registers a masking algorithm
    /// of the filtering program `program` (a `Blast-filter-program` number) and returns its id.
    /// `name` is the name of the mask file for makeblastdb (with `-gi_mask`, the name of the
    /// algorithm's GI-mask files). The first registration creates the mask column; with
    /// `-gi_mask` every registration creates the GI masks of the algorithm instead.
    pub fn register_mask_algorithm(
        &mut self,
        program: i32,
        options: &[u8],
        name: &[u8],
    ) -> Result<i32, WriteDbError> {
        let id = self.mask.register_program(
            program,
            options,
            name,
            !self.volumes.is_empty(),
            self.max_file_size,
            &self.column_dates,
        )?;
        self.note_mask_column();
        Ok(id)
    }

    /// The `use_gi_mask` argument of `CWriteDB` (`makeblastdb -gi_mask`): the masking algorithms
    /// write GI-indexed files (`super::gi_mask`) instead of the mask column, and the database
    /// gets an alias file whatever its number of volumes. Called before any registration.
    pub fn use_gi_masks(&mut self) {
        self.mask.use_gi_masks();
    }

    /// `RegisterMaskAlgorithm` of an id, a description and options: the form of
    /// `makeblastdb -mask_id`.
    pub fn register_mask_algorithm_id(
        &mut self,
        id: &[u8],
        description: &[u8],
        options: &[u8],
    ) -> Result<i32, WriteDbError> {
        let id = self
            .mask
            .register_id(id, description, options, !self.volumes.is_empty())?;
        self.note_mask_column();
        Ok(id)
    }

    /// `SetMaskData`: the masked ranges of the sequence added last (its mask blobs
    /// are written with it when it is published; with `-gi_mask` they go to the GI masks under
    /// the sequence's GIs `gis`, and nowhere for a sequence without one). The ranges are those
    /// `MaskedRangeSet::get_ranges` gives.
    pub fn set_mask_data(
        &mut self,
        ranges: &[MaskedRanges],
        gis: &[i64],
    ) -> Result<(), WriteDbError> {
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:964-987
        // `x_ComputeSeqLength` of a sequence from a `CBioseq` is the length of its Seq-inst, of
        // a sequence in the form of the sequence file `WriteDB_FindSequenceLength`.
        let seq_length = match &self.pending {
            Some(Pending::Bioseq(bioseq)) => u32::try_from(bioseq.length).ok(),
            Some(Pending::Raw(binary)) if !binary.seq.is_empty() => {
                Some(find_sequence_length(self.protein, &binary.seq))
            }
            _ => None,
        };
        self.mask
            .set_mask_data(ranges, gis, seq_length, !self.volumes.is_empty())?;
        self.note_mask_column();
        Ok(())
    }

    /// `m_UseGiMask`: whether the masks go to GI-indexed files (`use_gi_masks`).
    pub fn uses_gi_masks(&self) -> bool {
        self.mask.uses_gi_masks()
    }

    /// Puts the mask column in the list of columns when the mask writer has just made it (the
    /// first registration or `SetMaskData`): it takes the next column id.
    fn note_mask_column(&mut self) {
        if self.mask.column_meta().is_some()
            && !self
                .columns
                .iter()
                .any(|slot| matches!(slot, ColumnSlot::Mask))
        {
            self.columns.push(ColumnSlot::Mask);
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1482-1491
    /// `FindColumn`: the id of the column titled `title`, the mask column's included.
    pub fn find_column(&self, title: &[u8]) -> Option<usize> {
        self.columns.iter().position(|slot| match slot {
            ColumnSlot::Mask => title == MASK_DATA_TITLE,
            ColumnSlot::User { title: own, .. } => own == title,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1493-1519
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb.cpp:180-183
    // A volume has at most 36 columns: the error comes when the volume makes them
    // (`Volume::create_column`).
    /// `CreateUserColumn`: a new column titled `title` (not one that `find_column` finds), and
    /// its id.
    pub fn create_user_column(&mut self, title: &[u8]) -> Result<usize, WriteDbError> {
        if !self.volumes.is_empty() {
            return Err(unsupported_late_column());
        }
        self.columns.push(ColumnSlot::User {
            title: title.to_vec(),
            meta: ColumnMeta::new(),
            blobs: [Vec::new(), Vec::new()],
            have: 0,
        });
        Ok(self.columns.len() - 1)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1521-1535
    /// `AddColumnMetaData`: sets `key` to `value` in the metadata of column `col_id`.
    pub fn add_column_meta_data(
        &mut self,
        col_id: usize,
        key: &[u8],
        value: &[u8],
    ) -> Result<(), WriteDbError> {
        let volumes_exist = !self.volumes.is_empty();
        let meta = match self.columns.get_mut(col_id) {
            None => return Err(invalid_column_id()),
            Some(_) if volumes_exist => return Err(unsupported_late_column()),
            Some(ColumnSlot::User { meta, .. }) => meta,
            Some(ColumnSlot::Mask) => self
                .mask
                .existing_column_meta_mut()
                .ok_or_else(invalid_column_id)?,
        };
        meta.insert(key.to_vec(), value.to_vec());
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1537-1555
    // The mask column's blobs are made by `SetMaskData`.
    /// `SetBlobData` of `col_id` filled with `data` (`Clear`, `WriteRaw`): the next blob of column
    /// `col_id` for the sequence added last.
    pub fn set_blob_data(&mut self, col_id: usize, data: &[u8]) -> Result<(), WriteDbError> {
        match self.columns.get_mut(col_id) {
            None => Err(invalid_column_id()),
            Some(ColumnSlot::Mask) => Err(WriteDbError::Unsupported(
                "the blobs of the mask column are made from the masked ranges (SetMaskData)"
                    .to_string(),
            )),
            Some(ColumnSlot::User { blobs, have, .. }) => {
                if *have > 1 {
                    return Err(WriteDbError::Creation(
                        "Error: Already have blob for this sequence and column".to_string(),
                    ));
                }
                *have += 1;
                let blob = &mut blobs[usize::from(*have) - 1];
                blob.clear();
                blob.extend_from_slice(data);
                Ok(())
            }
        }
    }

    /// Test hook: the date that the columns of volume `index` store instead of the time of
    /// their creation.
    #[cfg(test)]
    pub fn set_column_date(&mut self, index: usize, date: &[u8]) {
        self.column_dates.set(index, date);
    }

    /// `oid_masks` of the `CWriteDB` constructor (`-oid_masks`), before the first sequence.
    pub fn set_oid_masks(&mut self, oid_masks: u64) {
        self.oid_masks = oid_masks;
    }

    /// The date for a column when the clock cannot be read (the time of the build).
    pub fn set_column_date_fallback(&mut self, date: &[u8]) {
        self.column_dates.set_fallback(date);
    }

    fn alias_name(&self) -> Vec<u8> {
        let extension: &[u8] = if self.protein { b".pal" } else { b".nal" };
        [self.dbname.as_slice(), extension].concat()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:405-438
    // The `ofstream` is not tested (`super::ofstream`): an alias file that cannot be opened or
    // written leaves the build going on, and the new database cannot be opened afterwards. The
    // stream writes the file when it is destroyed.
    //
    // With `-gi_mask` a database of one volume gets an alias file too. NCBI's `DBLIST` there is
    // the database's own absolute name (`m_Dbname`): the alias file names itself through a
    // directory, which NCBI's reader takes as another alias file, and every read of the database
    // ends with "DB alias files are mutually recursive" (NCBI's own makeblastdb too, at its `.njs`
    // step: exit status 2). LOSAT writes the database's file name without a directory there,
    // as the list of several volumes names each volume: NCBI's reader takes a `DBLIST` entry
    // without a directory and with the alias file's own name as the volume
    // (`CSeqDBAliasNode::x_ExpandAliases`, seqdbalias.cpp:764-797), so the database opens and
    // the build goes on as a build of several volumes does (approved exception 25 of
    // PD-LOSAT-NCBI-DEFECTS, the Owner's decision of 2026-11-20; NCBI BLAST+ 2.17.0 reads such a
    // database: `docs/evidence/gi_mask_single_volume/`). `MASKLIST` names the GI masks that got
    // GIs, each followed by a space.
    /// `x_MakeAlias`: the alias file of a database of several volumes, or of a database with GI
    /// masks.
    fn make_alias(&self) {
        let mut dblist = Vec::new();
        if self.volumes.len() > 1 {
            for (i, volume) in self.volumes.iter().enumerate() {
                if i > 0 {
                    dblist.push(b' ');
                }
                dblist.extend_from_slice(file_name(volume.volume_name()));
            }
        } else {
            dblist.extend_from_slice(file_name(&self.dbname));
        }
        let mut masklist = Vec::new();
        for name in self.mask.gi_mask_names().filter(|name| !name.is_empty()) {
            masklist.extend_from_slice(name);
            masklist.push(b' ');
        }
        let mut alias =
            OfStream::create(&crate::utils::os_bytes::path_of_bytes(&self.alias_name()));
        alias.insert(b"#\n# Alias file created: ");
        alias.insert(&self.date);
        alias.insert(b"\n#\n");
        alias.insert(b"TITLE ");
        alias.insert(&self.title);
        alias.insert(b"\n");
        alias.insert(b"DBLIST ");
        alias.insert(&dblist);
        alias.insert(b"\n");
        if !masklist.is_empty() {
            alias.insert(b"MASKLIST ");
            alias.insert(&masklist);
            alias.insert(b"\n");
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1655-1674
    /// `ListVolumes`: the volumes written so far (paths without extension).
    pub fn list_volumes(&self) -> Vec<Vec<u8>> {
        self.volumes
            .iter()
            .map(|volume| volume.volume_name().to_vec())
            .collect()
    }

    // The LMDB file is named for version 5 whether it was created or not; the
    // other files of `CWriteDB_LMDB` and `CWriteDB_TaxID` (`.pos`, `.ptf`,
    // `.pot`, `.pto`) are not, so a failed build of version 5 leaves them
    // (makeblastdb 2.17.0, an error in the FASTA input after the first
    // sequence: `db.pot`, `db.ptf` and `db.pto` stay). Nor are the files of the GI masks, or the
    // alias file of a database of one volume with GI masks: a failed build with `-gi_mask` leaves
    // them (survey b1/w30, b2/x03_oob_gi_gimask).
    /// `ListFiles`: every file of the database written so far.
    pub fn list_files(&self) -> Vec<Vec<u8>> {
        let mut files = Vec::new();
        for volume in &self.volumes {
            volume.list_files(&mut files);
        }
        if self.volumes.len() > 1 {
            files.push(self.alias_name());
        }
        // makeblastdb's names are never empty (`build_lmdb_file_name`'s error).
        if self.version == DbVersion::V5 {
            files.extend(build_lmdb_file_name(&self.dbname, self.protein).ok());
        }
        files
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1032-1041
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1042-1048
/// `x_CookSequence`, with the conversion's text on the standard output appended to `cout`.
fn cook_sequence(bioseq: &Bioseq, cout: &mut Vec<u8>) -> Result<BinarySequence, WriteDbError> {
    let Some(seq_data) = &bioseq.seq_data else {
        return Err(WriteDbError::Creation(
            "No sequence data in Bioseq, and no Bioseq_Handle available.".to_string(),
        ));
    };
    seq_data_to_binary(seq_data, bioseq.length as u32, cout).map_err(|_| {
        WriteDbError::Creation(format!(
            "Unable to process sequence for entry [{}].",
            bioseq.first_id_string()
        ))
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_lmdb.cpp:383-392
// A sequence whose ids are all GIs gives no element of `m_list` (`x_InsertEntry` returns
// at a GI), so one between two sequences with other ids leaves a gap in the OIDs of the list.
// `x_CreateOidToSeqidsLookupFile` runs in `~CWriteDB_LMDB`, where the exception ends
// makeblastdb 2.17.0 with SIGABRT after the volumes are written (survey B p11_gionly_mid).
// LOSAT refuses the database there with an error that names the sequence (makedb-v5
// DECISIONS 10: no NCBI abort reproduced). GI-only sequences before the first sequence with
// keys or after the last one leave no gap and are written as makeblastdb writes them.
/// The refusal of a version 5 database with `-parse_seqids` in which the sequence of OID `oid`
/// (first id `name`) has GIs only and lies between sequences with other ids.
fn gi_only_between_keys_error(oid: u32, name: &[u8]) -> WriteDbError {
    WriteDbError::GiOnlyBetweenKeys {
        oid,
        name: name.to_vec(),
    }
}

/// The text of [`WriteDbError::GiOnlyBetweenKeys`] (after makeblastdb's creation error prefix).
pub(crate) fn gi_only_between_keys_message(oid: u32, name: &[u8]) -> Vec<u8> {
    let mut message = b"Input id list not in ascending oid order: the sequence ".to_vec();
    message.extend_from_slice(name);
    message.extend_from_slice(
        format!(
            " (OID {oid}) has only GI ids and lies between sequences with other ids, which \
             makeblastdb 2.17.0 cannot write into the .pos file of a version 5 database (it \
             aborts); LOSAT refuses the database (-blastdb_version 4 writes it)"
        )
        .as_bytes(),
    );
    message
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:749-761
// The writer calls it with the defaults of `use_index` (false) and `index` (0)
// (include/objtools/blast/seqdb_reader/impl/seqdb_lmdb.hpp:118): one LMDB file for the
// database, beside its volumes.
// An empty name throws `std::invalid_argument` (makeprofiledb `-out ''`: its version 5 writer
// stops at the first sequence it publishes).
/// `BuildLMDBFileName` of `basename` and `is_protein`: `<db>.pdb` or `<db>.ndb`.
pub(crate) fn build_lmdb_file_name(
    basename: &[u8],
    is_protein: bool,
) -> Result<Vec<u8>, WriteDbError> {
    if basename.is_empty() {
        return Err(WriteDbError::InvalidArgument("Basename is empty"));
    }
    let extension: &[u8] = if is_protein { b".pdb" } else { b".ndb" };
    Ok([basename, extension].concat())
}

// NCBI reference (598d8ae6): c++/src/objects/blastdb/Blast_def_line.cpp:97-125
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:895-900
// Each line's values are those of `BlastDefline::tax_ids`, gathered in one vector for the whole
// set (no set per line or per sequence): a line drops its 0s only when it has another value.
/// `m_TaxIds` of a sequence with `-parse_seqids`: the union of `GetTaxIds` of each definition
/// line, ascending.
fn union_of_tax_ids(deflines: &[crate::seqid::BlastDefline]) -> Vec<i32> {
    let mut tax_ids = Vec::new();
    for defline in deflines {
        let start = tax_ids.len();
        tax_ids.extend(defline.taxid);
        if let Some(links) = &defline.links {
            tax_ids.extend_from_slice(links);
        }
        let line = &tax_ids[start..];
        if line.iter().any(|&id| id != line[0]) {
            let mut kept = start;
            for at in start..tax_ids.len() {
                if tax_ids[at] != 0 {
                    tax_ids[kept] = tax_ids[at];
                    kept += 1;
                }
            }
            tax_ids.truncate(kept);
        }
    }
    tax_ids.sort_unstable();
    tax_ids.dedup();
    tax_ids
}

/// `GetName` of a `CDirEntry` of `path`: the part after the last path separator.
fn file_name(path: &[u8]) -> &[u8] {
    let is_separator = |c: &u8| *c == b'/' || (cfg!(windows) && (*c == b'\\' || *c == b':'));
    match path.iter().rposition(is_separator) {
        Some(at) => &path[at + 1..],
        None => path,
    }
}

#[cfg(test)]
#[path = "writer_tests.rs"]
mod tests;
