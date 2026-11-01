//! The lookup of an id in the tables of a database for `losat blastdbcmd`: NCBI's
//! `CSeqDBImpl` calls `AccessionToOids`, `AccessionsToOids` and `GiToOidwFilterCheck`, over the
//! LMDB files of a version 5 database (`CSeqDBLMDBSet`) or the ISAM files of the volumes of a
//! version 4 database (`CSeqDBVol::AccessionToOids`). The search programs look up the Seq-align
//! id of a reported sequence as NCBI's formatters find it again (`SeqidLookup`, `SeqidToOid`).

use std::sync::OnceLock;

use super::file::DbFile;
use super::filter::{entries_accessions_to_oids, LmdbEntry};
use super::isam::{seqdb_simplify_seqid, SeqDbIdType, SimplifiedAccession, StringIsam};
use super::isam_num::NumericIsam;
use super::lmdb_read::Env;
use super::volume::with_extension;
use super::{SeqDb, SeqDbError, SeqType};
use crate::seqid::{label, ParsedSeqId};

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2724-2728
// ```c++
//     if ((! fits_in_four) && needs_four) {
//         NCBI_THROW(CSeqDBException,
//                    eArgErr,
//                    "ID overflows range of specified type.");
//     }
// ```
/// The exception of a version 4 lookup whose number does not fit 32 bits.
pub const ID_OVERFLOWS: &str = "ID overflows range of specified type.";

/// The lookup tables of a database. A database restricted by a list (`SeqDb::restrict`) finds
/// only the OIDs the list includes, as NCBI's `x_CheckOrFindOID`; `losat blastdbcmd` restricts
/// none.
pub struct IdLookup<'a> {
    db: &'a SeqDb,
    /// The LMDB files of a version 5 database; none for version 4.
    entries: Vec<LmdbEntry>,
    isams: Isams<'a>,
}

/// The ISAMs of each volume (GIs; for version 4 also the strings and the PIGs), each opened at
/// the first lookup that reaches the volume.
pub(super) struct VolumeIsams {
    gi: Vec<OnceLock<Option<NumericIsam>>>,
    strings: Vec<OnceLock<Option<StringIsam>>>,
    pig: Vec<OnceLock<Option<NumericIsam>>>,
    /// Whether a string ISAM keeps its index file in memory (`StringIsam::keep_index`), as NCBI
    /// maps it: the search programs' lookups of the reported sequences (`SeqDb::seqid_lookup`).
    index_in_memory: bool,
}

impl VolumeIsams {
    fn new(volumes: usize, index_in_memory: bool) -> Self {
        Self {
            gi: cells(volumes),
            strings: cells(volumes),
            pig: cells(volumes),
            index_in_memory,
        }
    }
}

impl std::fmt::Debug for VolumeIsams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("VolumeIsams")
            .field("volumes", &self.gi.len())
            .finish()
    }
}

/// The ISAMs of an `IdLookup`: its own (`losat blastdbcmd`) or the database's, kept for the
/// search (`SeqDb::seqid_lookup`).
enum Isams<'a> {
    Own(VolumeIsams),
    Shared(&'a VolumeIsams),
}

impl std::ops::Deref for Isams<'_> {
    type Target = VolumeIsams;

    fn deref(&self) -> &VolumeIsams {
        match self {
            Isams::Own(isams) => isams,
            Isams::Shared(isams) => isams,
        }
    }
}

impl SeqDb {
    /// The lookup tables of this database (for version 5, its LMDB files are read here).
    pub fn id_lookup(&self) -> Result<IdLookup<'_>, SeqDbError> {
        let entries = if self.is_version5() {
            self.lmdb_entries()?
        } else {
            Vec::new()
        };
        Ok(IdLookup {
            db: self,
            entries,
            isams: Isams::Own(VolumeIsams::new(self.volumes.len(), false)),
        })
    }
}

fn cells<T>(count: usize) -> Vec<OnceLock<T>> {
    (0..count).map(|_| OnceLock::new()).collect()
}

/// The value of `cell`, opened by `open` when the cell is empty.
fn opened<T>(
    cell: &OnceLock<Option<T>>,
    open: impl FnOnce() -> Result<Option<T>, SeqDbError>,
) -> Result<Option<&T>, SeqDbError> {
    if cell.get().is_none() {
        let _ = cell.set(open()?);
    }
    Ok(cell.get().and_then(Option::as_ref))
}

impl IdLookup<'_> {
    /// `x_CheckOrFindOID(oid) && oid == oid0`: an OID of the database that its list, if any,
    /// includes.
    fn visible(&self, oid: u32) -> bool {
        match &self.db.filter {
            Some(filter) => filter.contains(oid),
            None => oid < self.db.volume_starts[self.db.volumes.len()],
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1237-1246
    // ```c++
    //     if (m_LMDBSet.IsBlastDBVersion5()) {
    //     	vector<int>  tmp;
    //     	m_LMDBSet.AccessionToOids(acc, tmp);
    //    	    for(unsigned int i=0; i < tmp.size(); i++) {
    //    	    	int oid2 = tmp[i];
    //             if (x_CheckOrFindOID(oid2, locked) && (tmp[i] == oid2)) {
    //             	oids.push_back(tmp[i]);
    //             }
    //         }
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:275-286
    // ```c++
    // void CSeqDBLMDBSet::AccessionToOids(const string & acc, vector<TOid>  & oids)const
    // {
    // 	m_LMDBEntrySet[0]->AccessionToOids(acc, oids);
    // 	vector <TOid> tmp;
    // 	for(unsigned int i=1; i < m_LMDBEntrySet.size(); i++) {
    // 		m_LMDBEntrySet[i]->AccessionToOids(acc, tmp);
    // 		if(!tmp.empty()){
    // 			oids.insert(oids.end(), tmp.begin(), tmp.end());
    // 			tmp.clear();
    // 		}
    // 	}
    // }
    // ```
    /// Every OID of the key `acc` in a version 5 database: the LMDB files in the database's
    /// order, inside a file the order of LMDB's duplicates (`memcmp` of the 4 little-endian
    /// bytes).
    pub fn accession_to_oids(&self, acc: &[u8]) -> Result<Vec<u32>, SeqDbError> {
        let mut oids = Vec::new();
        for entry in &self.entries {
            for oid in entry.oids(acc)? {
                if let Some(oid) = entry.database_oid(oid).filter(|&oid| self.visible(oid)) {
                    oids.push(oid);
                }
            }
        }
        Ok(oids)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1347-1361
    // ```c++
    // void CSeqDBImpl::AccessionsToOids(const vector<string>& accs, vector<blastdb::TOid>& oids)
    // {
    //     CHECK_MARKER();
    //     oids.clear();
    //     oids.resize(accs.size());
    //     if (m_LMDBSet.IsBlastDBVersion5()) {
    //     	m_LMDBSet.AccessionsToOids(accs,oids);
    //     	CSeqDBLockHold locked(m_Atlas);
    //     	for(unsigned int i=0; i < oids.size(); i++) {
    //     		if(oids[i] == kSeqDBEntryNotFound) {
    //     			continue;
    //     		}
    //     		blastdb::TOid oid2 = oids[i];
    //     		if (!(x_CheckOrFindOID(oid2, locked) && (oids[i] == oid2))) {
    //     			oids[i] = kSeqDBEntryNotFound;
    //     		}
    //     	}
    //     }
    // ```
    /// One OID for each of `accs` in a version 5 database (`None`: not found), the first
    /// duplicate of the key in the last LMDB file that has it. On an error `oids` holds what
    /// NCBI's vector holds at its exception: the answers of the LMDB files before the one that
    /// failed, unfiltered.
    pub fn accessions_to_oids(
        &self,
        accs: &[Vec<u8>],
        oids: &mut Vec<Option<u32>>,
    ) -> Result<(), SeqDbError> {
        entries_accessions_to_oids(&self.entries, accs, oids)?;
        for slot in oids.iter_mut() {
            if slot.is_some_and(|oid| !self.visible(oid)) {
                *slot = None;
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1247-1275
    // ```c++
    //     else {
    //         vector<int> vol_oids;
    //
    //         for(int vol_idx = 0; vol_idx < m_VolSet.GetNumVols(); vol_idx++) {
    //             // Append any additional OIDs from this volume's indices.
    //             m_VolSet.GetVol(vol_idx)->AccessionToOids(acc,
    //                                                       vol_oids,
    //                                                       locked);
    //
    //             if (vol_oids.empty()) {
    //                 continue;
    //             }
    //
    //             int vol_start = m_VolSet.GetVolOIDStart(vol_idx);
    //
    //             ITERATE(vector<int>, iter, vol_oids) {
    //                 int oid1 = ((*iter) + vol_start);
    //                 int oid2 = oid1;
    //
    //                 // Remove OIDs already found in OIDs.
    //
    //                 if (find(oids.begin(), oids.end(), oid1) != oids.end()) {
    //                     continue;
    //                 }
    //
    //                 // Filter out any oids not in the virtual oid bitmaps.
    //
    //                 if (x_CheckOrFindOID(oid2, locked) && (oid1 == oid2)) {
    //                     oids.push_back(oid1);
    //                 }
    //             }
    //
    //             vol_oids.clear();
    //         }
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2791-2802
    // ```c++
    // void CSeqDBVol::AccessionToOids(const string         & acc,
    //                                 vector<int>          & oids,
    //                                 CSeqDBLockHold       & locked) const
    // {
    //     Int8   ident   (-1);
    //     string str_id;
    //     bool   simpler (false);
    //
    //     ESeqDBIdType id_type = SeqDB_SimplifyAccession(acc, ident, str_id, simpler);
    //
    //     x_StringToOids(acc, id_type, ident, str_id, simpler, oids);
    //
    // }
    // ```
    /// Every OID of the user's id `acc` in a version 4 database, `simplified` its
    /// `SeqDB_SimplifyAccession` (which every volume computes alike): the volumes in order, in a
    /// volume the order of its ISAM (the lines of the string ISAM sort the OIDs of one key as
    /// decimal texts, `10` before `9`), an OID found before left out. `None` is NCBI's exception
    /// [`ID_OVERFLOWS`], after the lookup in the first volume.
    ///
    /// `gnl|BL_ORD_ID|n` is the one sequence `n` of the whole database (line BDC DECISIONS 12):
    /// NCBI adds `n` to the start of every volume.
    pub fn accession_to_oids_v4(
        &self,
        acc: &[u8],
        simplified: &SimplifiedAccession,
    ) -> Result<Option<Vec<u32>>, SeqDbError> {
        if simplified.id_type == SeqDbIdType::Oid {
            let oid = u32::try_from(simplified.num_id).ok();
            return Ok(Some(
                oid.filter(|&oid| self.visible(oid)).into_iter().collect(),
            ));
        }
        let mut oids = Vec::new();
        for at in 0..self.db.volumes.len() {
            let vol_oids = self.volume_string_to_oids(at, acc, simplified)?;
            let ident = simplified.num_id;
            if !(ident == -1 || ident >> 32 == 0) {
                return Ok(None);
            }
            let vol_start = i64::from(self.db.volume_starts[at]);
            for oid in vol_oids {
                let Some(oid) = u32::try_from(i64::from(oid) + vol_start)
                    .ok()
                    .filter(|&oid| self.visible(oid))
                else {
                    continue;
                };
                if !oids.contains(&oid) {
                    oids.push(oid);
                }
            }
        }
        Ok(Some(oids))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2638-2722
    // ```c++
    //     switch(id_type) {
    //     case eStringId:
    //         x_OpenStrFile();
    //         if (! m_IsamStr.Empty()) {
    //             // Not simplified
    //             vcheck = true;
    //             m_IsamStr->StringToOids(str_id, oids, simpler, vcheck);
    //             x_UnleaseStrFile();
    //         }
    //         break;
    //
    //     case ePigId:
    //         // Converted to PIG type.
    //         x_OpenPigFile();
    //         if (! m_IsamPig.Empty()) {
    //             int oid(-1);
    //
    //             if (m_IsamPig->PigToOid((int) ident, oid)) {
    //                 oids.push_back(oid);
    //             }
    //             x_UnleasePigFile();
    //         }
    //         break;
    //
    //     case eGiId:
    //         // Converted to GI type.
    //         x_OpenGiFile();
    //         if (! m_IsamGi.Empty()) {
    //             int oid(-1);
    //
    //             if (m_IsamGi->IdToOid(ident, oid)) {
    //                 oids.push_back(oid);
    //             }
    //             x_UnleaseGiFile();
    //         }
    //         break;
    //
    //     case eTiId:
    //         // Converted to TI type.
    //          x_OpenTiFile();
    //         if (! m_IsamTi.Empty()) {
    //             ...
    //         }
    //         else {
    //         	x_OpenStrFile();
    //         	if (m_IsamStr.NotEmpty()) {
    //             	// Not every database with TIs has a TI index, so fall
    //             	// back to a string comparison if the first attempt fails.
    //             	//
    //             	// 1. TI's don't have versions.
    //             	// 2. Specify "adjusted" as true, because lookup of
    //             	//    "gb|.." and similar tricks are not needed for TIs.
    //         		m_IsamStr->StringToOids(acc, oids, true, vcheck);
    //         		x_UnleaseStrFile();
    //         	}
    //         }
    //         break;
    // ```
    /// `CSeqDBVol::x_StringToOids` of the volume `at` without its final checks (the size of the
    /// number, then `x_CheckVersions`, which [`StringIsam::string_id_to_oids`] runs): the OIDs
    /// local to the volume. A volume with a trace-id ISAM (`.nti`, which makeblastdb never
    /// writes) is refused when a `ti` id reaches it.
    fn volume_string_to_oids(
        &self,
        at: usize,
        acc: &[u8],
        simplified: &SimplifiedAccession,
    ) -> Result<Vec<i32>, SeqDbError> {
        let volume = &self.db.volumes[at];
        if volume.num_oids() == 0 {
            return Ok(Vec::new());
        }
        let seqtype = self.db.seqtype;
        let protein = seqtype == SeqType::Protein;
        // NCBI reads the index of a volume whose files fail `x_InitSearch` all the same
        // (`x_FindIndexBounds` in the constructor) and crashes (probe I of line BDC: an empty,
        // short or garbage index, another version or type, a data file of another length).
        let string_isam = || {
            opened(&self.isams.strings[at], || {
                let name = volume.base().to_string_lossy();
                if !StringIsam::index_exists(&name, protein) {
                    return Ok(None);
                }
                let isam = StringIsam::open_volume(&name, protein)?;
                if !isam.is_initialized() {
                    return Err(SeqDbError::Unsupported {
                        path: volume.base().to_path_buf(),
                        message: "a string ISAM whose header does not fit its files (NCBI \
                                  crashes on it)"
                            .to_string(),
                    });
                }
                if self.isams.index_in_memory {
                    return Ok(Some(isam.keep_index()?));
                }
                Ok(Some(isam))
            })
        };
        let local = |oids: Vec<u32>| oids.into_iter().map(|oid| oid as i32).collect();
        let found = |oid: Option<i32>| oid.into_iter().collect();
        Ok(match simplified.id_type {
            SeqDbIdType::String => match string_isam()? {
                Some(isam) => local(isam.string_id_to_oids(acc, simplified, &mut |oid| {
                    let set = volume.defline_set(oid)?;
                    Ok(set
                        .deflines
                        .into_iter()
                        .flat_map(|defline| defline.ids)
                        .collect())
                })?),
                None => Vec::new(),
            },
            SeqDbIdType::Pig => {
                let isam = opened(&self.isams.pig[at], || {
                    NumericIsam::open_pig_volume(volume.base(), seqtype)
                })?;
                match isam {
                    Some(isam) => found(isam.id_to_oid(i64::from(simplified.num_id as i32))?),
                    None => Vec::new(),
                }
            }
            SeqDbIdType::Gi => match self.gi_isam(at)? {
                Some(isam) => found(isam.id_to_oid(simplified.num_id)?),
                None => Vec::new(),
            },
            SeqDbIdType::Ti => {
                let index = with_extension(volume.base(), seqtype, "ti");
                if index.exists() && with_extension(volume.base(), seqtype, "td").exists() {
                    return Err(SeqDbError::Unsupported {
                        path: index,
                        message: "an ISAM of trace ids".to_string(),
                    });
                }
                match string_isam()? {
                    Some(isam) => local(isam.string_to_oids(acc, true, false)?.oids),
                    None => Vec::new(),
                }
            }
            SeqDbIdType::Oid => vec![simplified.num_id as i32],
        })
    }

    /// The GI ISAM of the volume `at` (`x_OpenGiFile`).
    fn gi_isam(&self, at: usize) -> Result<Option<&NumericIsam>, SeqDbError> {
        let volume = &self.db.volumes[at];
        if volume.num_oids() == 0 {
            return Ok(None);
        }
        opened(&self.isams.gi[at], || {
            NumericIsam::open_volume(volume.base(), self.db.seqtype)
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1190-1205
    // ```c++
    // bool CSeqDBImpl::GiToOidwFilterCheck(TGi gi, int & oid)
    // {
    //     CHECK_MARKER();
    //     CSeqDBLockHold locked(m_Atlas);
    //
    //     for(int i = 0; i < m_VolSet.GetNumVols(); i++) {
    //     	oid =-1;
    //         if (m_VolSet.GetVol(i)->GiToOid(gi, oid, locked)) {
    //            	oid += m_VolSet.GetVolOIDStart(i);
    //             int oid0 = oid;
    //             if (CheckOrFindOID(oid) && (oid==oid0)) {
    //             	return true;
    //             }
    //         }
    //     }
    //     return false;
    // }
    // ```
    /// The OID of the GI `gi`: the first volume whose file has it.
    pub fn gi_to_oid(&self, gi: i64) -> Result<Option<u32>, SeqDbError> {
        for at in 0..self.db.volumes.len() {
            let Some(isam) = self.gi_isam(at)? else {
                continue;
            };
            let Some(oid) = isam.id_to_oid(gi)? else {
                continue;
            };
            let oid = i64::from(oid) + i64::from(self.db.volume_starts[at]);
            if let Some(oid) = u32::try_from(oid).ok().filter(|&oid| self.visible(oid)) {
                return Ok(Some(oid));
            }
        }
        Ok(None)
    }
}

/// The LMDB files of a version 5 database, open for every [`SeqidLookup`] of a search
/// (`SeqDb::seqid_lookup`).
pub(super) struct LookupFiles {
    entries: Vec<LmdbEntry>,
    files: Vec<DbFile>,
}

impl std::fmt::Debug for LookupFiles {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LookupFiles")
            .field("files", &self.files.len())
            .finish()
    }
}

/// NCBI's lookup of a reported sequence's Seq-align id (`seqid_to_oid`), with the files of the
/// database kept open for the search (`SeqDb::seqid_lookup`): the LMDB files of version 5, here
/// as environments that keep the pages they read, and the ISAMs of the volumes, whose string
/// index stays in memory.
pub struct SeqidLookup<'a> {
    /// The ISAMs of the volumes: the strings of version 4, the GIs and the PIGs.
    volumes: IdLookup<'a>,
    /// Each LMDB file of the database with its environment (version 5), in the database's order.
    envs: Vec<(&'a LmdbEntry, Env<'a, DbFile>)>,
}

impl SeqDb {
    /// The lookup of the Seq-align ids of the reported sequences of a search. The database's
    /// LMDB files (version 5) and ISAMs are opened at the first lookup of the search and stay
    /// open; each `SeqidLookup` reads the LMDB files through environments of its own.
    pub fn seqid_lookup(&self) -> Result<SeqidLookup<'_>, SeqDbError> {
        let mut envs = Vec::new();
        if self.is_version5() {
            let files = match self.lookup_files.get() {
                Some(files) => files,
                None => {
                    let entries = self.lmdb_entries()?;
                    let files = entries
                        .iter()
                        .map(LmdbEntry::open_file)
                        .collect::<Result<Vec<_>, _>>()?;
                    let _ = self.lookup_files.set(LookupFiles { entries, files });
                    self.lookup_files.get().expect("set above")
                }
            };
            for (entry, file) in files.entries.iter().zip(&files.files) {
                envs.push((entry, entry.open_env(file)?));
            }
        }
        let isams = self
            .lookup_isams
            .get_or_init(|| VolumeIsams::new(self.volumes.len(), true));
        Ok(SeqidLookup {
            volumes: IdLookup {
                db: self,
                entries: Vec::new(),
                isams: Isams::Shared(isams),
            },
            envs,
        })
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2679-2696
// ```c++
// bool IsStringId(const CSeq_id & id)
// {
// 	 switch(id.Which()) {
// 	    case CSeq_id::e_Gi:
// 	    	return false;
// 	    break;
// 	    case CSeq_id::e_General:
// 	    {
// 	    	const CDbtag & dbt = id.GetGeneral();
//             if (dbt.CanGetDb() && (dbt.GetDb() == "PIG")) {
//             	return false;
//             }
// 	    }
// 	    default:
// 	    	return true;
// 	    break;
// 	 };
// }
// ```
/// `IsStringId`: whether a version 5 database looks `id` up in its LMDB files (a `general` id
/// other than `PIG` falls through to the default).
fn is_string_id(id: &ParsedSeqId) -> bool {
    match id {
        ParsedSeqId::Gi(_) => false,
        ParsedSeqId::General(tag) => tag.db != b"PIG",
        _ => true,
    }
}

impl SeqidLookup<'_> {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:903-916
    // ```c++
    // bool CSeqDB::SeqidToOid(const CSeq_id & seqid, int & oid) const
    // {
    //     ////m_Impl->Verify();
    //     bool rv = false;
    //
    //     oid = -1;
    //
    //     vector<int> oids;
    //     m_Impl->SeqidToOids(seqid, oids, false);
    //
    //     if (! oids.empty()) {
    //         rv = true;
    //         oid = oids[0];
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1381-1459
    // ```c++
    // void CSeqDBImpl::SeqidToOids(const CSeq_id & seqid_in,
    //                              vector<int>   & oids,
    //                              bool            multi)
    // {
    //     ...
    //     if (m_LMDBSet.IsBlastDBVersion5() && (!is_BL_ORD_ID)) {
    //     	if(IsStringId(seqid_in)) {
    //     		vector<int>  tmp;
    //     		if(seqid_in.IsPir() || seqid_in.IsPrf()) {
    //     			m_LMDBSet.AccessionToOids(seqid_in.AsFastaString(), tmp);
    //     		}
    //     		else {
    //     			m_LMDBSet.AccessionToOids(seqid_in.GetSeqIdString(true), tmp);
    //     		}
    //     		for(unsigned int i=0; i < tmp.size(); i++) {
    //     			int oid2 = tmp[i];
    //     			if (x_CheckOrFindOID(oid2, locked) && (tmp[i] == oid2)) {
    //     				oids.push_back(tmp[i]);
    //     			}
    //     		}
    //     		return;
    //     	}
    //     }
    //     ...
    //     for(int vol_idx = 0; vol_idx < m_VolSet.GetNumVols(); vol_idx++) {
    //        	// Append any additional OIDs from this volume's indices.
    //        	m_VolSet.GetVol(vol_idx)->SeqidToOids(seqid, vol_oids, locked);
    //
    //        	if (vol_oids.empty()) {
    //            	continue;
    //        	}
    //
    //        	int vol_start = m_VolSet.GetVolOIDStart(vol_idx);
    //
    //        	ITERATE(vector<int>, iter, vol_oids) {
    //            	int oid1 = ((*iter) + vol_start);
    //            	int oid2 = oid1;
    //
    //            	// Filter out any oids not in the virtual oid bitmaps.
    //
    //            	if (x_CheckOrFindOID(oid2, locked) && (oid1 == oid2)) {
    //                	oids.push_back(oid1);
    //
    //                	if (! multi) {
    //                    	return;
    //                	}
    //            	}
    //        	}
    // ```
    /// `CSeqDB::SeqidToOid`: the OID that NCBI's lookup of the Seq-id `seqid` gives first, the
    /// one whose Bioseq NCBI's BLAST database data loader serves for it
    /// (`CBlastDbDataLoader::x_GetOid`; the search's own database, so its list's OIDs only).
    /// Version 5 looks a text id up by its key in every LMDB file (the files in order, in a file
    /// LMDB's order of the duplicates), version 4 and the GIs and PIGs in the ISAMs of the
    /// volumes in order, and `gnl|BL_ORD_ID|n` is `n` from the start of each volume in turn.
    /// `None`: not found, also where NCBI's lookup throws on a number that does not fit 32 bits.
    /// An error of a file is NCBI's exception of the lookup.
    pub fn seqid_to_oid(&self, seqid: &ParsedSeqId) -> Result<Option<u32>, SeqDbError> {
        let db = self.volumes.db;
        let bl_ord_id = matches!(seqid, ParsedSeqId::General(tag) if tag.db == b"BL_ORD_ID");
        if db.is_version5() && !bl_ord_id && is_string_id(seqid) {
            let mut key = Vec::new();
            if matches!(seqid, ParsedSeqId::Pir(_) | ParsedSeqId::Prf(_)) {
                label::fasta(seqid, &mut key);
            } else {
                label::content(seqid, true, &mut key);
            }
            let mut oids = Vec::new();
            for (entry, env) in &self.envs {
                let found = entry.oids_in(env, &key)?;
                oids.extend(found.into_iter().filter_map(|oid| entry.database_oid(oid)));
            }
            return Ok(oids.into_iter().find(|&oid| self.volumes.visible(oid)));
        }
        let simplified =
            seqdb_simplify_seqid(seqid).map_err(|error| SeqDbError::NotFound(error.to_string()))?;
        let mut fasta = Vec::new();
        label::fasta(seqid, &mut fasta);
        for at in 0..db.volumes.len() {
            let vol_oids = self
                .volumes
                .volume_string_to_oids(at, &fasta, &simplified)?;
            let ident = simplified.num_id;
            if !(ident == -1 || ident >> 32 == 0) {
                return Ok(None);
            }
            let vol_start = i64::from(db.volume_starts[at]);
            for oid in vol_oids {
                if let Some(oid) = u32::try_from(i64::from(oid) + vol_start)
                    .ok()
                    .filter(|&oid| self.volumes.visible(oid))
                {
                    return Ok(Some(oid));
                }
            }
        }
        Ok(None)
    }
}
