//! The lookup of an id in the tables of a version 5 database: NCBI's `CSeqDBLMDBSet` and the
//! `CSeqDBImpl` calls over it (`AccessionToOids`, `AccessionsToOids`,
//! `GiToOidwFilterCheck`), for `losat blastdbcmd`.

use std::cell::OnceCell;

use super::filter::{entries_accessions_to_oids, LmdbEntry};
use super::isam_num::NumericIsam;
use super::{SeqDb, SeqDbError};

/// The lookup tables of a database that no list restricts.
pub struct IdLookup<'a> {
    db: &'a SeqDb,
    entries: Vec<LmdbEntry>,
    /// The GI ISAM of each volume, opened at the first lookup that reaches the volume.
    gi_isams: Vec<OnceCell<Option<NumericIsam>>>,
}

impl SeqDb {
    /// The lookup tables of this version 5 database.
    pub fn id_lookup(&self) -> Result<IdLookup<'_>, SeqDbError> {
        Ok(IdLookup {
            db: self,
            entries: self.lmdb_entries()?,
            gi_isams: self.volumes.iter().map(|_| OnceCell::new()).collect(),
        })
    }
}

impl IdLookup<'_> {
    /// `x_CheckOrFindOID(oid) && oid == oid0` of a database without a list.
    fn visible(&self, oid: u32) -> bool {
        oid < self.db.volume_starts[self.db.volumes.len()]
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
    /// Every OID of the key `acc`: the LMDB files in the database's order, inside a file the
    /// order of LMDB's duplicates (`memcmp` of the 4 little-endian bytes).
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
    /// One OID for each of `accs` (`None`: not found), the first duplicate of the key in the
    /// last LMDB file that has it. On an error `oids` holds what NCBI's vector holds at its
    /// exception: the answers of the LMDB files before the one that failed, unfiltered.
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
        for (at, volume) in self.db.volumes.iter().enumerate() {
            if volume.num_oids() == 0 {
                continue;
            }
            let cell = &self.gi_isams[at];
            if cell.get().is_none() {
                let opened = NumericIsam::open_volume(volume.base(), self.db.seqtype)?;
                let _ = cell.set(opened);
            }
            let Some(Some(isam)) = cell.get() else {
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
