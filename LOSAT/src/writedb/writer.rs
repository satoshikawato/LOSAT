//! The database being written (`CWriteDB_Impl`,
//! `c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp`): sequences are
//! published one behind the caller, volumes are cut when a file would pass
//! the size limit, and closing renames a single volume or writes the alias
//! file of several.

use std::path::PathBuf;

use super::bioseq::Bioseq;
use super::convert::{seq_data_to_binary, BinarySequence};
pub use super::files::DbVersion;
use super::files::OffsetArrays;
use super::header::{defline_set, BadTitleChar};
use super::lmdb::{
    lmdb_file_name, lmdb_map_size_from_env, LmdbEnvSettings, WriteDbLmdb, WriteDbTaxId,
};
use super::platform::Platform;
use super::volume::Volume;
use super::WriteDbError;

/// `CWriteDB_Impl` for a database without parsed Seq-ids, of format version 4
/// or 5.
pub struct WriteDb {
    dbname: String,
    protein: bool,
    title: Vec<u8>,
    date: Vec<u8>,
    max_file_size: u64,
    /// `m_DbVersion`.
    version: DbVersion,
    /// `m_VolumeList`; the last element is `m_Volume`.
    volumes: Vec<Volume>,
    closed: bool,
    /// `m_Bioseq` while `m_HaveSequence` is set.
    pending: Option<Bioseq>,
    /// What NCBI's serializer posts to the diagnostic stream, in order.
    diagnostics: Vec<u8>,
    /// `m_Lmdbdb` (version 5, from the first published sequence until `Close`).
    lmdb: Option<WriteDbLmdb>,
    /// `m_Taxdb` (as `m_Lmdbdb`).
    taxdb: Option<WriteDbTaxId>,
    /// `m_LmdbOid`: the OID of the next sequence over all volumes.
    lmdb_oid: u32,
    /// The platform whose makeblastdb the files follow.
    platform: Platform,
}

impl std::fmt::Debug for WriteDb {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("WriteDb")
            .field("dbname", &self.dbname)
            .field("protein", &self.protein)
            .field("version", &self.version)
            .field("volumes", &self.volumes.len())
            .field("closed", &self.closed)
            .field("lmdb_oid", &self.lmdb_oid)
            .field("platform", &self.platform)
            .finish_non_exhaustive()
    }
}

// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1864-1885
// ```c++
//     if ( fix_method != eFNP_Replace ) {
//         string message;
//         if (io != NULL) {
//             message += io->GetStackTrace() + "\n";
//         }
//         message += "Bad char [0x" +
//                    NStr::NumericToString((unsigned char)c,0,16)+
//                    "] in string";
//         if (io != NULL) {
//             message += " at " + io->GetPosition();
//         }
//         if (!str.empty()) {
//             message += "\n" + str;
//         }
//         switch (fix_method) {
//         case eFNP_ReplaceAndWarn:
// #if 0
//             CNcbiDiag(eDiag_Error, eDPF_Default)
//                 << ErrCode(NCBI_ERRCODE_X, 7) << message << Endm;
// #else
//             ERR_POST_X(7, message);
// #endif
// ```
// makeblastdb prints, for each byte of a title that a VisibleString cannot
// hold and each time the definition line is serialized:
// "Error: (803.7) [makeblastdb] Blast-def-line-set.E.title", "Bad char [0xC3]
// in string at byte 21" and the title.
fn post_bad_title_chars(diagnostics: &mut Vec<u8>, title: &[u8], bad: &[BadTitleChar]) {
    for bad_char in bad {
        diagnostics.extend_from_slice(
            format!(
                "Error: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0x{:X}] in \
                 string at byte {}\n",
                bad_char.byte, bad_char.stream_offset
            )
            .as_bytes(),
        );
        diagnostics.extend_from_slice(title);
        diagnostics.push(b'\n');
    }
}

impl WriteDb {
    /// `CWriteDB(dbname, seqtype, title, ..., dbver, ...)` followed by
    /// `SetMaxFileSize`. `date` is the creation date of the index files, in
    /// the form of `current_date`'s argument; `platform` the platform whose
    /// makeblastdb the files follow ([`Platform::host`] for a build).
    pub fn new(
        dbname: &str,
        protein: bool,
        title: &[u8],
        date: &[u8],
        max_file_size: u64,
        version: DbVersion,
        platform: Platform,
    ) -> Self {
        Self {
            dbname: dbname.to_string(),
            protein,
            title: title.to_vec(),
            date: date.to_vec(),
            max_file_size,
            version,
            volumes: Vec::new(),
            closed: false,
            pending: None,
            diagnostics: Vec::new(),
            lmdb: None,
            taxdb: None,
            lmdb_oid: 0,
            platform,
        }
    }

    /// `m_DbVersion`.
    pub fn version(&self) -> DbVersion {
        self.version
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1563-1566
    // ```c++
    // void CWriteDB_Impl::SetMaxFileSize(Uint8 sz)
    // {
    //     m_MaxFileSize = sz;
    // }
    // ```
    /// `SetMaxFileSize`: the limit of the files of the volumes created from
    /// now on.
    pub fn set_max_file_size(&mut self, max_file_size: u64) {
        self.max_file_size = max_file_size;
    }

    /// Takes the diagnostics posted since the last call (text for the
    /// standard error stream).
    pub fn take_diagnostics(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.diagnostics)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:473-475
    // ```c++
    //     CRef<CBlast_def_line_set> headers =
    //         CWriteDB::ExtractBioseqDeflines(*bs, m_ParseIDs, m_LongIDs,
    //                 m_ScanBioseq4CFastaReaderUsrObjct);
    // ```
    // The call serializes the definition line once (writedb_impl.cpp:887-893,
    // with the sequence's own id), which posts the diagnostics of its title;
    // the bytes are not kept.
    /// `CWriteDB::ExtractBioseqDeflines` as far as it is visible: the
    /// diagnostics of the title.
    pub fn extract_bioseq_deflines(&mut self, bioseq: &Bioseq) {
        let mut bad = Vec::new();
        defline_set(&bioseq.title, 0, &mut bad);
        post_bad_title_chars(&mut self.diagnostics, &bioseq.title, &bad);
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:158-180
    // ```c++
    // void CWriteDB_Impl::AddSequence(const CBioseq & bs)
    // {
    //     // Publish previous sequence
    //     x_Publish();
    //
    //     // Blank slate for new sequence.
    //     x_ResetSequenceData();
    //
    //     m_Bioseq.Reset(& bs);
    //     if (m_Bioseq->GetInst().CanGetMol() && (m_Bioseq->IsAa() != m_Protein)) {
    //         CNcbiOstrstream msg;
    //         msg << "Invalid molecule type of sequence added ("
    //             << (m_Bioseq->IsAa() ? "protein" : "nucleotide")
    //             << "); expected " << (m_Protein ? "protein" : "nucleotide");
    //         NCBI_THROW(CWriteDBException, eArgErr, CNcbiOstrstreamToString(msg));
    //     }
    //
    //     if (m_Indices & CWriteDB::eAddHash) {
    //         x_ComputeHash(bs);
    //     }
    //
    //     x_SetHaveSequence();
    // }
    // ```
    /// `AddSequence(const CBioseq &)`: publishes the previous sequence and
    /// keeps this one.
    pub fn add_sequence(&mut self, bioseq: Bioseq) -> Result<(), WriteDbError> {
        self.publish()?;
        if bioseq.is_aa != self.protein {
            let name = |protein: bool| if protein { "protein" } else { "nucleotide" };
            return Err(WriteDbError::Creation(format!(
                "Invalid molecule type of sequence added ({}); expected {}",
                name(bioseq.is_aa),
                name(self.protein)
            )));
        }
        self.pending = Some(bioseq);
        Ok(())
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
    /// `x_CookHeader`: the header record with the number the sequence gets
    /// in the current volume.
    fn cook_header(&mut self, title: &[u8]) -> Vec<u8> {
        let oid = self.volumes.last().map_or(0, Volume::oid);
        let mut bad = Vec::new();
        let binhdr = defline_set(title, oid, &mut bad);
        post_bad_title_chars(&mut self.diagnostics, title, &bad);
        binhdr
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1127-1133
    // ```c++
    //     if (x_HaveSequence()) {
    //         _ASSERT(! (m_Bioseq.Empty() && m_Sequence.empty()));
    //
    //         x_ClearHaveSequence();
    //     } else {
    //         return;
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1096-1098
    // ```c++
    //     x_CookHeader();
    //     x_CookSequence();
    //     x_CookColumns();
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1156-1239
    // ```c++
    //     x_CookData();
    //
    //     bool done = false;
    //
    //     if (! m_Volume.Empty()) {
    //         done = m_Volume->WriteSequence(m_Sequence,
    //                                        m_Ambig,
    //                                        m_BinHdr,
    //                                        m_Ids,
    //                                        m_Pig,
    //                                        m_Hash,
    //                                        m_Blobs,
    //                                        m_MaskDataColumn);
    //     ...
    //     if (! done) {
    //         int index = (int) m_VolumeList.size();
    //
    //         if (m_Volume.NotEmpty()) {
    //             m_Volume->Close();
    //         }
    //
    //         {
    //
    //             m_Volume.Reset(new CWriteDB_Volume(m_Dbname,
    //                                                m_Protein,
    //                                                m_Title,
    //                                                m_Date,
    //                                                index,
    //                                                m_MaxFileSize,
    //                                                m_MaxVolumeLetters,
    //                                                m_Indices,
    //                                                m_DbVersion,
    //                                                m_OidMasks));
    //
    //             m_VolumeList.push_back(m_Volume);
    //     ...
    //         // need to reset OID,  hense recalculate the header and id
    //         x_CookHeader();
    //
    //         done = m_Volume->WriteSequence(m_Sequence,
    //                                        m_Ambig,
    //                                        m_BinHdr,
    //                                        m_Ids,
    //                                        m_Pig,
    //                                        m_Hash,
    //                                        m_Blobs,
    //                                        m_MaskDataColumn);
    //     ...
    //         if (! done) {
    //             NCBI_THROW(CWriteDBException,
    //                        eArgErr,
    //                        "Cannot write sequence to volume.");
    //         }
    //     }
    // ```
    // The size check of the old volume uses the header cooked with the old
    // volume's number; the record written to a new volume is cooked again
    // with number 0.
    //
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1169-1175
    // ```c++
    //         if (done  &&  (m_DbVersion == eBDB_Version5)  &&  m_Lmdbdb) {
    //         	if (m_ParseIDs) {
    //         		m_Lmdbdb->InsertEntries(m_Ids,m_LmdbOid);
    //         	}
    //             m_Taxdb->InsertEntries(m_TaxIds, m_LmdbOid);
    //             m_LmdbOid++;
    //         }
    // ```
    // The same follows the write into a new volume (writedb_impl.cpp:1226-1232).
    // Without `-parse_seqids` no Seq-id goes to `m_Lmdbdb`. `m_TaxIds` is the
    // set of the tax ids of the sequence's definition lines: makeblastdb sets
    // tax id 0 on each of them without `-taxid` and `-taxid_map`
    // (`CTaxIdSet::FixTaxId`, taxid_set.cpp:84-136), and `InsertEntries` records
    // tax id 0 for the set {0} as for an empty set (writedb_lmdb.cpp:446-464).
    /// `x_Publish`: writes the kept sequence, into a new volume when the
    /// current one is full.
    fn publish(&mut self) -> Result<(), WriteDbError> {
        let Some(bioseq) = self.pending.take() else {
            return Ok(());
        };
        if self.version == DbVersion::V5 && self.lmdb.is_none() {
            self.open_lmdb_files()?;
        }
        let mut binhdr = self.cook_header(&bioseq.title);
        let binary = self.cook_sequence(&bioseq)?;

        let mut done = false;
        if let Some(volume) = self.volumes.last_mut() {
            done = volume.write_sequence(&binary.seq, &binary.ambig, &binhdr)?;
            if done {
                self.insert_lmdb_entries()?;
            }
        }
        if !done {
            let index = self.volumes.len() as u32;
            // The new volume's index file takes the arrays of the closed one
            // (`OffsetArrays`).
            let mut arrays = OffsetArrays::default();
            if let Some(volume) = self.volumes.last_mut() {
                volume.close()?;
                arrays = volume.take_offset_arrays();
            }
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
            )?);
            binhdr = self.cook_header(&bioseq.title);
            let volume = self.volumes.last_mut().expect("a volume was pushed");
            done = volume.write_sequence(&binary.seq, &binary.ambig, &binhdr)?;
            if done {
                self.insert_lmdb_entries()?;
            }
            if !done {
                return Err(WriteDbError::Creation(
                    "Cannot write sequence to volume.".to_string(),
                ));
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.hpp:655
    // ```c++
    //     int m_LmdbOid;
    // ```
    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdbcommon.hpp:58
    // ```c++
    // typedef Int4 TOid;
    // ```
    // An OID is a 4-byte signed integer in NCBI's writer and in its LMDB files, so a database
    // has at most 2^31 sequences; NCBI's `m_LmdbOid++` overflows past that. LOSAT stops with
    // an error at the sequence that would get OID 2^31 (no NCBI overflow reproduced).
    /// The tax ids of the sequence just written, under the next OID of the
    /// database (version 5).
    fn insert_lmdb_entries(&mut self) -> Result<(), WriteDbError> {
        if let Some(taxdb) = self.taxdb.as_mut() {
            if self.lmdb_oid > i32::MAX as u32 {
                return Err(WriteDbError::Other(format!(
                    "a version 5 database holds at most {} sequences (an OID is an Int4)",
                    i32::MAX as u64 + 1
                )));
            }
            taxdb.insert_entries(&[], self.lmdb_oid);
            self.lmdb_oid += 1;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1136-1154
    // ```c++
    //     if(m_DbVersion == eBDB_Version5 && m_Lmdbdb.Empty()) {
    //         const string lmdb_fname_w_path = BuildLMDBFileName(m_Dbname, m_Protein);
    //         Uint8 map_size = 0;
    //         char* map_sz_str = getenv("BLASTDB_LMDB_MAP_SIZE");
    //         if (map_sz_str) {
    //         	map_size = NStr::StringToUInt8(map_sz_str);
    //         }
    //         if(map_size > 0){
    //         	m_Lmdbdb.Reset(new CWriteDB_LMDB(lmdb_fname_w_path, map_size));
    //         	m_Taxdb.Reset(new CWriteDB_TaxID(
    //         		          GetFileNameFromExistingLMDBFile(lmdb_fname_w_path, ELMDBFileType::eTaxId2Offsets),
    //         		          map_size));
    //         }
    //         else {
    //         	m_Lmdbdb.Reset(new CWriteDB_LMDB(lmdb_fname_w_path));
    //         	m_Taxdb.Reset(new CWriteDB_TaxID(
    //         		          GetFileNameFromExistingLMDBFile(lmdb_fname_w_path, ELMDBFileType::eTaxId2Offsets)));
    //         }
    //     }
    // ```
    // The page size of LMDB, its default map sizes and the length of its files on Windows are
    // the platform's (DECISIONS 8 of the makedb-v5 line, `Platform`).
    /// Creates the LMDB files of a version 5 database (`<db>.pdb` and
    /// `<db>.ptf`, or `.ndb` and `.ntf`) at its first published sequence.
    fn open_lmdb_files(&mut self) -> Result<(), WriteDbError> {
        let lmdb_fname_w_path = build_lmdb_file_name(&self.dbname, self.protein);
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
            std::path::Path::new(&lmdb_fname_w_path),
            lmdb_map,
            &platform,
            &settings,
        )?;
        self.lmdb = Some(lmdb);
        let taxdb = WriteDbTaxId::new(
            &lmdb_file_name(std::path::Path::new(&lmdb_fname_w_path), "tf"),
            tax_map,
            &platform,
            &settings,
        )?;
        self.taxdb = Some(taxdb);
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1032-1041
    // ```c++
    //         default:
    //             msg = "Unable to process sequence for entry [";
    //             msg += (m_Bioseq->GetId().front())->GetSeqIdString(false);
    //             msg += "].";
    //         }
    //
    //         if (! msg.empty()) {
    //             NCBI_THROW(CWriteDBException, eArgErr, msg);
    //         }
    //     } else {
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1042-1048
    // ```c++
    //         int sz = m_SeqVector.size();
    //
    //         if (sz == 0) {
    //             NCBI_THROW(CWriteDBException,
    //                        eArgErr,
    //                        "No sequence data in Bioseq, "
    //                        "and no Bioseq_Handle available.");
    // ```
    /// `x_CookSequence`.
    fn cook_sequence(&self, bioseq: &Bioseq) -> Result<BinarySequence, WriteDbError> {
        let Some(seq_data) = &bioseq.seq_data else {
            return Err(WriteDbError::Creation(
                "No sequence data in Bioseq, and no Bioseq_Handle available.".to_string(),
            ));
        };
        seq_data_to_binary(seq_data, bioseq.length as u32).map_err(|_| {
            WriteDbError::Creation(format!(
                "Unable to process sequence for entry [{}].",
                bioseq.local_id
            ))
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:326-354
    // ```c++
    //     if (m_Closed)
    //         return;
    //
    //     m_Closed = true;
    //
    //     x_Publish();
    //     m_Sequence.erase();
    //     m_Ambig.erase();
    //
    //     if (! m_Volume.Empty()) {
    //         m_Volume->Close();
    //
    //         if (m_UseGiMask) {
    //             for (unsigned int i=0; i<m_GiMasks.size(); ++i) {
    //                 m_GiMasks[i]->Close();
    //             }
    //         }
    //
    //         if (m_VolumeList.size() == 1) {
    //             m_Volume->RenameSingle();
    //         }
    //         else if(m_VolumeList.size() > 100){
    //         	unsigned int num_digits = log10(m_VolumeList.size()) +1;
    //         	for(unsigned i=0; i < m_VolumeList.size(); i++) {
    //         		CRef<CWriteDB_Volume> & v = m_VolumeList[i];
    //         		v->RenameFileIndex(num_digits);
    //         	}
    //         	LOG_POST(Info << "Rename files index to " << num_digits << " digits");
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:378-394
    // ```c++
    //         if (m_VolumeList.size() > 1 || m_UseGiMask) {
    //             x_MakeAlias();
    //         }
    //         if ((m_DbVersion == eBDB_Version5)  &&  m_Lmdbdb) {
    //         	vector<string> vol_names(m_VolumeList.size());
    //         	vector<blastdb::TOid> vol_num_oids(m_VolumeList.size());
    //         	for(unsigned i=0; i < m_VolumeList.size(); i++) {
    //         		CRef<CWriteDB_Volume> & v = m_VolumeList[i];
    //         		vol_names[i] = CDirEntry(v->GetVolumeName()).GetName();
    //         		vol_num_oids[i] = v->GetOID();
    //         	}
    //             m_Lmdbdb->InsertVolumesInfo(vol_names, vol_num_oids);
    //             m_Lmdbdb.Reset();
    //             m_Taxdb.Reset();
    //         }
    //
    //         m_Volume.Reset();
    //     }
    // ```
    // The volume names are the final ones (after `RenameSingle` or
    // `RenameFileIndex`). `Reset` runs the destructors of the two LMDB
    // writers, which write their lookup files and their LMDB transactions
    // (`WriteDbLmdb::close`, `WriteDbTaxId::close`).
    /// `Close`: publishes the last sequence, closes the last volume, gives
    /// the files their final names, and writes the LMDB files of version 5.
    pub fn close(&mut self) -> Result<(), WriteDbError> {
        if self.closed {
            return Ok(());
        }
        self.closed = true;
        self.publish()?;
        let count = self.volumes.len();
        if count == 0 {
            return Ok(());
        }
        self.volumes[count - 1].close()?;
        drop(self.volumes[count - 1].take_offset_arrays());
        if count == 1 {
            self.volumes[0].rename_single()?;
        } else if count > 100 {
            let num_digits = count.ilog10() + 1;
            for volume in &mut self.volumes {
                volume.rename_file_index(num_digits)?;
            }
        }
        if count > 1 {
            self.make_alias()?;
        }
        if self.version == DbVersion::V5 {
            if let Some(mut lmdb) = self.lmdb.take() {
                let vol_names: Vec<String> = self
                    .volumes
                    .iter()
                    .map(|volume| file_name(volume.volume_name()).to_string())
                    .collect();
                let vol_num_oids: Vec<u32> = self.volumes.iter().map(Volume::oid).collect();
                lmdb.insert_volumes_info(&vol_names, &vol_num_oids)?;
                lmdb.close()?;
                if let Some(taxdb) = self.taxdb.take() {
                    taxdb.close()?;
                }
            }
        }
        Ok(())
    }

    /// The volumes written: the path of each without extension, its number
    /// of sequences and its letters (what NCBI's reader takes from the index
    /// files for the totals of the database).
    pub fn volume_totals(&self) -> Vec<(String, u32, u64)> {
        self.volumes
            .iter()
            .map(|volume| {
                (
                    volume.volume_name().to_string(),
                    volume.oid(),
                    volume.letters(),
                )
            })
            .collect()
    }

    fn alias_name(&self) -> String {
        format!(
            "{}.{}",
            self.dbname,
            if self.protein { "pal" } else { "nal" }
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:405-434
    // ```c++
    //     string dblist;
    //     if (m_VolumeList.size() > 1) {
    //         for(unsigned i = 0; i < m_VolumeList.size(); i++) {
    //             if (dblist.size())
    //                 dblist += " ";
    //
    //         		CRef<CWriteDB_Volume> & v = m_VolumeList[i];
    //         		dblist += CDirEntry(v->GetVolumeName()).GetName();
    //         }
    //     } else {
    //         dblist = m_Dbname;
    //     }
    //     ...
    //     string nm = x_MakeAliasName();
    //
    //     ofstream alias(nm.c_str());
    //
    //     alias << "#\n# Alias file created: " << m_Date  << "\n#\n"
    //           << "TITLE "        << m_Title << "\n"
    //           << "DBLIST "       << dblist  << "\n";
    // ```
    /// `x_MakeAlias`: the alias file of a database of several volumes.
    fn make_alias(&self) -> Result<(), WriteDbError> {
        let mut text = Vec::new();
        text.extend_from_slice(b"#\n# Alias file created: ");
        text.extend_from_slice(&self.date);
        text.extend_from_slice(b"\n#\nTITLE ");
        text.extend_from_slice(&self.title);
        text.extend_from_slice(b"\nDBLIST ");
        for (i, volume) in self.volumes.iter().enumerate() {
            if i > 0 {
                text.push(b' ');
            }
            text.extend_from_slice(file_name(volume.volume_name()).as_bytes());
        }
        text.push(b'\n');
        let name = self.alias_name();
        std::fs::write(&name, text).map_err(|source| WriteDbError::Io {
            path: PathBuf::from(name),
            source,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1655-1674
    // ```c++
    // void CWriteDB_Impl::ListVolumes(vector<string> & vols)
    // {
    //     vols.clear();
    //
    //     ITERATE(vector< CRef<CWriteDB_Volume> >, iter, m_VolumeList) {
    //         vols.push_back((**iter).GetVolumeName());
    //     }
    // }
    //
    // void CWriteDB_Impl::ListFiles(vector<string> & files)
    // {
    //     files.clear();
    //
    //     ITERATE(vector< CRef<CWriteDB_Volume> >, iter, m_VolumeList) {
    //         (**iter).ListFiles(files);
    //     }
    //
    //     if (m_VolumeList.size() > 1) {
    //         files.push_back(x_MakeAliasName());
    //     }
    // ```
    /// `ListVolumes`: the volumes written so far (paths without extension).
    pub fn list_volumes(&self) -> Vec<String> {
        self.volumes
            .iter()
            .map(|volume| volume.volume_name().to_string())
            .collect()
    }

    // The LMDB file is named for version 5 whether it was created or not; the
    // other files of `CWriteDB_LMDB` and `CWriteDB_TaxID` (`.pos`, `.ptf`,
    // `.pot`, `.pto`) are not, so a failed build of version 5 leaves them
    // (makeblastdb 2.17.0, an error in the FASTA input after the first
    // sequence: `db.pot`, `db.ptf` and `db.pto` stay).
    /// `ListFiles`: every file of the database written so far.
    pub fn list_files(&self) -> Vec<String> {
        let mut files = Vec::new();
        for volume in &self.volumes {
            volume.list_files(&mut files);
        }
        if self.volumes.len() > 1 {
            files.push(self.alias_name());
        }
        if self.version == DbVersion::V5 {
            files.push(build_lmdb_file_name(&self.dbname, self.protein));
        }
        files
    }

    /// Drops every open file without flushing (before a failed build removes
    /// its files).
    pub fn abandon(&mut self) {
        self.closed = true;
        self.pending = None;
        for volume in &mut self.volumes {
            volume.abandon();
        }
        self.lmdb = None;
        self.taxdb = None;
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:749-761
// ```c++
// string BuildLMDBFileName(const string& basename, bool is_protein, bool use_index, unsigned int index)
// {
//     if (basename.empty()) {
//         throw invalid_argument("Basename is empty");
//     }
//
//     string vol_str=kEmptyStr;
//     if(use_index) {
//     	vol_str = (index > 9) ?".": ".0";
//     	vol_str += NStr::UIntToString(index);
//     }
//     return basename + vol_str + (is_protein ? ".pdb" : ".ndb");
// }
// ```
// The writer calls it with the defaults `use_index = false`, `index = 0`
// (include/objtools/blast/seqdb_reader/impl/seqdb_lmdb.hpp:118): one LMDB file for the
// database, beside its volumes.
/// `BuildLMDBFileName(basename, is_protein)`: `<db>.pdb` or `<db>.ndb`.
pub(crate) fn build_lmdb_file_name(basename: &str, is_protein: bool) -> String {
    format!("{basename}{}", if is_protein { ".pdb" } else { ".ndb" })
}

/// `CDirEntry(path).GetName()`: the part after the last path separator.
fn file_name(path: &str) -> &str {
    let is_separator = |c: char| c == '/' || (cfg!(windows) && (c == '\\' || c == ':'));
    match path.rfind(is_separator) {
        Some(at) => &path[at + 1..],
        None => path,
    }
}
