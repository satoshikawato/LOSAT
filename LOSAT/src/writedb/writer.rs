//! The database being written (`CWriteDB_Impl`,
//! `c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp`): sequences are
//! published one behind the caller, volumes are cut when a file would pass
//! the size limit, and closing renames a single volume or writes the alias
//! file of several.

use std::path::PathBuf;

use super::bioseq::Bioseq;
use super::convert::{seq_data_to_binary, BinarySequence};
use super::header::{defline_set, BadTitleChar};
use super::volume::Volume;
use super::WriteDbError;

/// `CWriteDB_Impl` for a version 4 database without parsed Seq-ids.
#[derive(Debug)]
pub struct WriteDb {
    dbname: String,
    protein: bool,
    title: Vec<u8>,
    date: Vec<u8>,
    max_file_size: u64,
    /// `m_VolumeList`; the last element is `m_Volume`.
    volumes: Vec<Volume>,
    closed: bool,
    /// `m_Bioseq` while `m_HaveSequence` is set.
    pending: Option<Bioseq>,
    /// What NCBI's serializer posts to the diagnostic stream, in order.
    diagnostics: Vec<u8>,
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
    /// `CWriteDB(dbname, seqtype, title, ...)` followed by `SetMaxFileSize`.
    /// `date` is the creation date of the index files, in the form of
    /// `current_date`'s argument.
    pub fn new(dbname: &str, protein: bool, title: &[u8], date: &[u8], max_file_size: u64) -> Self {
        Self {
            dbname: dbname.to_string(),
            protein,
            title: title.to_vec(),
            date: date.to_vec(),
            max_file_size,
            volumes: Vec::new(),
            closed: false,
            pending: None,
            diagnostics: Vec::new(),
        }
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
    /// `x_Publish`: writes the kept sequence, into a new volume when the
    /// current one is full.
    fn publish(&mut self) -> Result<(), WriteDbError> {
        let Some(bioseq) = self.pending.take() else {
            return Ok(());
        };
        let mut binhdr = self.cook_header(&bioseq.title);
        let binary = self.cook_sequence(&bioseq)?;

        let mut done = false;
        if let Some(volume) = self.volumes.last_mut() {
            done = volume.write_sequence(&binary.seq, &binary.ambig, &binhdr)?;
        }
        if !done {
            let index = self.volumes.len() as u32;
            if let Some(volume) = self.volumes.last_mut() {
                volume.close()?;
            }
            self.volumes.push(Volume::create(
                &self.dbname,
                self.protein,
                &self.title,
                &self.date,
                index,
                self.max_file_size,
            )?);
            binhdr = self.cook_header(&bioseq.title);
            let volume = self.volumes.last_mut().expect("a volume was pushed");
            done = volume.write_sequence(&binary.seq, &binary.ambig, &binhdr)?;
            if !done {
                return Err(WriteDbError::Creation(
                    "Cannot write sequence to volume.".to_string(),
                ));
            }
        }
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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:378-380
    // ```c++
    //         if (m_VolumeList.size() > 1 || m_UseGiMask) {
    //             x_MakeAlias();
    //         }
    // ```
    /// `Close`: publishes the last sequence, closes the last volume, and
    /// gives the files their final names.
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
        Ok(())
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

    /// `ListFiles`: every file of the database written so far.
    pub fn list_files(&self) -> Vec<String> {
        let mut files = Vec::new();
        for volume in &self.volumes {
            volume.list_files(&mut files);
        }
        if self.volumes.len() > 1 {
            files.push(self.alias_name());
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
    }
}

/// `CDirEntry(path).GetName()`: the part after the last path separator.
fn file_name(path: &str) -> &str {
    let is_separator = |c: char| c == '/' || (cfg!(windows) && (c == '\\' || c == ':'));
    match path.rfind(is_separator) {
        Some(at) => &path[at + 1..],
        None => path,
    }
}
