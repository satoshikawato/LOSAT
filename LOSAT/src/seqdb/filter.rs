//! A database seen through a list of taxids (line DI-8: `-taxids`, `-taxidlist`,
//! `-negative_taxids`, `-negative_taxidlist`): the OIDs that NCBI's `CSeqDB` visits when it
//! is opened with a user list or a negative list (`CSeqDBOIDList`), the totals it scans for
//! them (`CSeqDBImpl::x_ScanTotals`) and the definition lines its filtered header read keeps
//! (`CSeqDBVol::x_GetFilteredHeader`).
//!
//! A restricted [`SeqDb`] shows its included OIDs as the ordinal ids `0..k`, in OID order
//! (DI-8 DECISIONS 7): every table of a program that an ordinal id indexes stays dense, and
//! [`SeqDb::oid`] gives the OID that NCBI prints as `gnl|BL_ORD_ID|n`. The OIDs are held as a
//! bitmap with a rank table, as NCBI holds them (`CSeqDBBitVector`), not as a list.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use super::file::DbFile;
use super::lmdb_read::Env;
use super::{tax_files, SeqDb, SeqDbError, SeqType};
use crate::seqid::{BlastDefline, DeflineSet};

/// The taxids of a search's list after NCBI's expansion (`s_GetTaxIDList`): `-taxids` and
/// `-taxidlist` give a positive list, `-negative_taxids` and `-negative_taxidlist` a negative
/// one. An empty set is no list at all (`CSearchDatabase::x_InitializeDb` opens the database
/// without one); [`SeqDb::restrict_taxids`] is not called for it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TaxidList {
    /// Search only the OIDs of these taxids (`CSeqDBGiList` with taxids).
    Positive(BTreeSet<i32>),
    /// Search every OID except those whose taxids all are among these
    /// (`CSeqDBNegativeList` with taxids).
    Negative(BTreeSet<i32>),
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:242-244
// ```c++
//     	if((user_list->GetNumTaxIds() > 0) && !(lmdb_set.IsBlastDBVersion5())) {
//     		NCBI_THROW(CSeqDBException, eArgErr, "Taxonomy filtering is not supported in v4 BLAST dbs");
//     	}
// ```
const V4_TAXID_ERROR: &str = "Taxonomy filtering is not supported in v4 BLAST dbs";

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:47-49
// ```c++
// 	if(vol_names.size()  > tmp_name.size()) {
// 		NCBI_THROW(CSeqDBException, eArgErr, "Input db vol size does not match lmdb vol size");
// 	}
// ```
const VOLUME_COUNT_ERROR: &str = "Input db vol size does not match lmdb vol size";

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:65-67
// ```c++
// 	if(m_OIDEnd == 0) {
// 		NCBI_THROW(CSeqDBException, eArgErr, "Input db vol does not match lmdb vol");
// 	}
// ```
const VOLUME_NAME_ERROR: &str = "Input db vol does not match lmdb vol";

/// The included OIDs of a restricted database: bit `oid % 64` of word `oid / 64`, with the
/// number of included OIDs before every 512th OID.
#[derive(Debug)]
pub(crate) struct OidFilter {
    words: Vec<u64>,
    /// `supers[s]`: the included OIDs below OID `512 * s`.
    supers: Vec<u32>,
    count: u32,
}

/// Words of one entry of `OidFilter::supers`.
const SUPER_WORDS: usize = 8;

impl OidFilter {
    /// The filter of the OIDs whose bits are set in `words` (bits at `num_oids` and above
    /// are cleared).
    fn new(mut words: Vec<u64>, num_oids: u32) -> Self {
        words.resize((num_oids as usize).div_ceil(64), 0);
        if num_oids % 64 != 0 {
            if let Some(last) = words.last_mut() {
                *last &= (1u64 << (num_oids % 64)) - 1;
            }
        }
        let mut supers = Vec::with_capacity(words.len().div_ceil(SUPER_WORDS));
        let mut count = 0u32;
        for (at, word) in words.iter().enumerate() {
            if at % SUPER_WORDS == 0 {
                supers.push(count);
            }
            count += word.count_ones();
        }
        Self {
            words,
            supers,
            count,
        }
    }

    /// The number of included OIDs: the ordinal ids are `0..count`.
    pub(crate) fn count(&self) -> u32 {
        self.count
    }

    /// Whether `oid` is included.
    pub(crate) fn contains(&self, oid: u32) -> bool {
        self.words
            .get((oid / 64) as usize)
            .is_some_and(|word| word >> (oid % 64) & 1 == 1)
    }

    /// The number of included OIDs below `oid`: the ordinal id of `oid` when it is included.
    pub(crate) fn rank(&self, oid: u32) -> u32 {
        let word = (oid / 64) as usize;
        if word >= self.words.len() {
            return self.count;
        }
        let first = word - word % SUPER_WORDS;
        let mut rank = self.supers[first / SUPER_WORDS]
            + self.words[first..word]
                .iter()
                .map(|w| w.count_ones())
                .sum::<u32>();
        let bits = oid % 64;
        if bits > 0 {
            rank += (self.words[word] & ((1u64 << bits) - 1)).count_ones();
        }
        rank
    }

    /// The OID of ordinal id `ordinal` (below `count`).
    pub(crate) fn select(&self, ordinal: u32) -> u32 {
        let block = self.supers.partition_point(|&rank| rank <= ordinal) - 1;
        let mut rest = ordinal - self.supers[block];
        let first = block * SUPER_WORDS;
        for (at, &word) in self.words[first..].iter().enumerate() {
            let ones = word.count_ones();
            if rest < ones {
                let mut word = word;
                for _ in 0..rest {
                    word &= word - 1;
                }
                return ((first + at) * 64) as u32 + word.trailing_zeros();
            }
            rest -= ones;
        }
        unreachable!("ordinal id {ordinal} of a filter of {} OIDs", self.count)
    }

    /// The first OID from `oid` on that is not included, `limit` when every OID of
    /// `oid..limit` is.
    pub(crate) fn run_end(&self, oid: u32, limit: u32) -> u32 {
        let mut at = oid;
        while at < limit {
            let word = (at / 64) as usize;
            let Some(&bits) = self.words.get(word) else {
                return at;
            };
            let excluded = !bits >> (at % 64);
            if excluded != 0 {
                return (at + excluded.trailing_zeros()).min(limit);
            }
            at = (word as u32 + 1) * 64;
        }
        limit
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1910-1941
// ```c++
// bool s_IncludeDefline_Taxid(const CBlast_def_line & def, const set<TTaxId> & user_tax_ids)
// {
// 	CBlast_def_line::TTaxIds tax_ids;
// 	if (def.IsSetTaxid()) {
// 	    tax_ids.insert(def.GetTaxid());
// 	}
// 	if(def.IsSetLinks()) {
// 		CBlast_def_line::TLinks leaf_ids = def.GetLinks();
// #ifdef NCBI_STRICT_TAX_ID
//         ITERATE(CBlast_def_line::TLinks, it, leaf_ids) tax_ids.insert(TAX_ID_FROM(int, *it));
// #else
//         tax_ids.insert(leaf_ids.begin(), leaf_ids.end());
// #endif
//     }
//
// 	if(user_tax_ids.size() > tax_ids.size()) {
// 		ITERATE(CBlast_def_line::TTaxIds, itr, tax_ids) {
// 			if(user_tax_ids.find(*itr) != user_tax_ids.end()) {
// 				return true;
// 			}
// 		}
//
// 	}
// 	else {
// 		ITERATE(set<TTaxId>, itr, user_tax_ids) {
// 			if(tax_ids.find(*itr) != tax_ids.end()) {
// 				return true;
// 			}
// 		}
// 	}
// 	return false;
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1956-1968
// ```c++
// bool s_IncludeDefline_NegativeTaxid(const CBlast_def_line & def, const set<TTaxId> & user_tax_ids)
// {
// 	CBlast_def_line::TTaxIds taxid_set = def.GetTaxIds();
// 	if(taxid_set.size() > user_tax_ids.size()) {
// 		return true;
// 	}
// 	ITERATE(CBlast_def_line::TTaxIds, itr, taxid_set) {
// 		if(user_tax_ids.find(*itr) == user_tax_ids.end()) {
// 			return true;
// 		}
// 	}
// 	return false;
// }
// ```
// The sets are the taxids that the database holds of those requested: NCBI swaps the list's
// set for them (seqdblmdbset.cpp:315-333, 335-354).
/// Which definition lines of an included OID a restricted database keeps.
#[derive(Debug)]
pub(crate) enum DeflineRule {
    /// A positive list: the lines whose `taxid` or `links` hold one of these taxids.
    Taxids(BTreeSet<i32>),
    /// A negative list: the lines with a taxid (`GetTaxIds`) outside these.
    NegativeTaxids(BTreeSet<i32>),
}

impl DeflineRule {
    /// Whether NCBI's filtered header read keeps `defline`.
    fn keeps(&self, defline: &BlastDefline) -> bool {
        match self {
            DeflineRule::Taxids(found) => defline
                .taxid
                .into_iter()
                .chain(defline.links.iter().flatten().copied())
                .any(|taxid| found.contains(&taxid)),
            DeflineRule::NegativeTaxids(found) => {
                defline.tax_ids().iter().any(|taxid| !found.contains(taxid))
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2059-2066
    // ```c++
    //             if (! have_memb) {
    //                 TBDLLIter eraseme = iter++;
    //                 dl.erase(eraseme);
    //                 asn_changed = true;
    //             } else {
    //                 iter++;
    //             }
    // ```
    // NCBI reads no Seq-id and no title for an OID whose every line is erased (`GetBioseq`
    // returns null, seqdbvol.cpp:1370-1376). makeblastdb writes the taxid lookup files from
    // the same lines, so an OID of the list has a line of the list; LOSAT keeps every line of
    // an OID that would have none.
    /// Which of the definition lines of `set` the read keeps, in order: every one when the
    /// rule would keep none.
    pub(crate) fn keep_mask(&self, set: &DeflineSet) -> Vec<bool> {
        let mask: Vec<bool> = set.deflines.iter().map(|d| self.keeps(d)).collect();
        if mask.iter().any(|&keep| keep) {
            mask
        } else {
            vec![true; mask.len()]
        }
    }
}

/// The definition lines of one sequence in their stored order, as a reader decodes them
/// (`Volume::deflines`, `defline_taxids`, `defline_set`), shortened by a `DeflineRule`.
pub(crate) trait DeflineList {
    /// Keeps the lines `i` with `keep[i]`; a list of another length is left as it is.
    fn retain_lines(&mut self, keep: &[bool]);
}

impl<T> DeflineList for Vec<T> {
    fn retain_lines(&mut self, keep: &[bool]) {
        if keep.len() == self.len() {
            let mut at = 0;
            self.retain(|_| {
                at += 1;
                keep[at - 1]
            });
        }
    }
}

impl DeflineList for DeflineSet {
    fn retain_lines(&mut self, keep: &[bool]) {
        self.deflines.retain_lines(keep);
    }
}

/// The volumes of a database that share one LMDB file (`CSeqDBLMDBEntry`).
struct LmdbEntry {
    /// The LMDB file (`.pdb`, `.ndb`).
    lmdb: PathBuf,
    /// The OID of the entry's first volume (`m_OIDStart`).
    oid_start: u32,
    /// The OIDs of the LMDB file's volumes, in its order: the end of each volume's OIDs
    /// (`max_oid`) and their count when the database does not hold the volume
    /// (`skipped_oids`, 0 for a held volume).
    volumes: Vec<(u32, u32)>,
    /// Whether the database holds only some of the LMDB file's volumes (`m_isPartial`).
    partial: bool,
    /// The number of OIDs of all the LMDB file's volumes (`m_NumOids`).
    num_oids: u32,
}

impl LmdbEntry {
    /// The file beside the LMDB file with the last two letters of its extension replaced
    /// (`GetFileNameFromExistingLMDBFile`, seqdb_lmdb.cpp:763-790): `tf`, `to`, `ot`.
    fn lookup_file(&self, letters: &str) -> PathBuf {
        let mut name = self.lmdb.as_os_str().to_string_lossy().into_owned();
        name.truncate(name.len().saturating_sub(2));
        name.push_str(letters);
        PathBuf::from(name)
    }

    /// The name NCBI's messages give the entry: the LMDB file's name without directory and
    /// extension (`GetDBFileNames`, `CSeqDB_Path::FindBaseName`).
    fn base_name(path: &Path) -> String {
        path.file_stem()
            .map(|stem| stem.to_string_lossy().into_owned())
            .unwrap_or_default()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:118-143
    // ```c++
    // void CSeqDBLMDBEntry::x_AdjustOidsOffset_TaxList(vector<TOid> & oids) const
    // {
    // 	if ((m_OIDStart > 0) && (!m_isPartial)) {
    // 		for(unsigned int i=0; i < oids.size(); i++) {
    // 				oids[i] +=m_OIDStart;
    // 		}
    // 	}
    //
    // 	if(m_isPartial) {
    // 		vector<TOid> tmp;
    // 		for(unsigned int i=0; i < oids.size(); i++) {
    // 			TOid skipped_oids = 0;
    // 			for(unsigned int j=0;j < m_VolInfo.size(); j++) {
    // 				if(oids[i] >=  m_VolInfo[j].max_oid) {
    // 					skipped_oids += m_VolInfo[j].skipped_oids;
    // 				}
    // 				else {
    // 					if(m_VolInfo[j].skipped_oids > 0) {
    // 					}
    // 					else {
    // 						tmp.push_back(oids[i] + m_OIDStart - skipped_oids) ;
    // 					}
    // 					break;
    // 				}
    // 			}
    // 		}
    // 		oids.swap(tmp);
    // 	}
    //
    // }
    // ```
    /// The database OIDs of the LMDB file's OIDs `oids`: the OIDs of volumes that the
    /// database does not hold are dropped.
    fn database_oids(&self, oids: Vec<u32>) -> Vec<u32> {
        if !self.partial {
            return oids.into_iter().map(|oid| oid + self.oid_start).collect();
        }
        let mut out = Vec::with_capacity(oids.len());
        for oid in oids {
            let mut skipped = 0;
            for &(max_oid, skipped_oids) in &self.volumes {
                if oid >= max_oid {
                    skipped += skipped_oids;
                } else {
                    if skipped_oids == 0 {
                        out.push(oid + self.oid_start - skipped);
                    }
                    break;
                }
            }
        }
        out
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:165-175
    // ```c++
    // void CSeqDBLMDBEntry::TaxIdsToOids(const set<TTaxId>& tax_ids, vector<blastdb::TOid>& rv, vector<TTaxId> & tax_ids_found) const
    // {
    // 	m_LMDB->GetOidsForTaxIds(tax_ids, rv, tax_ids_found);
    // 	x_AdjustOidsOffset_TaxList(rv);
    // }
    //
    // void CSeqDBLMDBEntry::NegativeTaxIdsToOids(const set<TTaxId>& tax_ids, vector<blastdb::TOid>& rv, vector<TTaxId> & tax_ids_found) const
    // {
    // 	m_LMDB->NegativeTaxIdsToOids(tax_ids, rv, tax_ids_found);
    // 	x_AdjustOidsOffset_TaxList(rv);
    // }
    // ```
    /// The database OIDs of `tax_ids` in this entry (`negative`: the OIDs the list
    /// excludes) and the taxids the entry holds of `tax_ids`.
    fn tax_oids(
        &self,
        negative: bool,
        tax_ids: &BTreeSet<i32>,
    ) -> Result<tax_files::TaxOids, SeqDbError> {
        let ntf_path = self.lookup_file("tf");
        if !ntf_path.is_file() {
            return Err(missing_file(&ntf_path));
        }
        let ntf_file = DbFile::open(&ntf_path)?;
        let ntf = Env::open(&ntf_file).map_err(|message| SeqDbError::Format {
            path: ntf_path.clone(),
            message,
        })?;
        let nto_path = self.lookup_file("to");
        let nto = DbFile::open(&nto_path)?;
        let tax_error = |error: tax_files::TaxError, path: &Path| match error {
            tax_files::TaxError::Data(message) => SeqDbError::Format {
                path: path.to_path_buf(),
                message,
            },
            other => SeqDbError::NotFound(other.message(&Self::base_name(&ntf_path))),
        };
        let mut found = if negative {
            let not_path = self.lookup_file("ot");
            let not = DbFile::open(&not_path)?;
            tax_files::negative_tax_ids_to_oids(&ntf, &nto, &not, self.num_oids, tax_ids)
                .map_err(|error| tax_error(error, &not_path))?
        } else {
            tax_files::get_oids_for_tax_ids(&ntf, &nto, self.num_oids, tax_ids)
                .map_err(|error| tax_error(error, &nto_path))?
        };
        found.oids = self.database_oids(std::mem::take(&mut found.oids));
        Ok(found)
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:85-89
// ```c++
//         CFile tf(fname);
// 		Uint8 readMapSize = (tf.GetLength()/10000 + 1) *10000;
//         if (!tf.Exists()) {
//             NCBI_THROW( CSeqDBException, eFileErr, "File " + fname + " not found. If you renamed any BLAST database files, please use original file names, and makeblastdb to rename the database. If you deleted any BLAST database files, you need to recreate the database.");
//         }
// ```
fn missing_file(path: &Path) -> SeqDbError {
    SeqDbError::NotFound(format!(
        "File {} not found. If you renamed any BLAST database files, please use original file \
         names, and makeblastdb to rename the database. If you deleted any BLAST database files, \
         you need to recreate the database.",
        path.display()
    ))
}

/// The OIDs of the sequences `first..end` of a volume read per window of the index.
const SCAN_WINDOW: u32 = 1 << 16;

impl SeqDb {
    /// The OID of ordinal id `ordinal`: `ordinal` itself unless the database is restricted
    /// ([`SeqDb::restrict_taxids`]). NCBI prints it as `n` of `gnl|BL_ORD_ID|n`.
    pub fn oid(&self, ordinal: u32) -> u32 {
        match &self.filter {
            Some(filter) => filter.select(ordinal),
            None => ordinal,
        }
    }

    /// Whether the database is restricted by a taxid list ([`SeqDb::restrict_taxids`]).
    pub fn is_restricted(&self) -> bool {
        self.filter.is_some()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:216-262
    // ```c++
    // CSeqDBLMDBSet::CSeqDBLMDBSet( const CSeqDBVolSet & volSet)
    // {
    // 	string lmdb_full_name = kEmptyStr;
    // 	vector<string> lmdb_vol_names;
    // 	...
    // 		CSeqDB_Path p(vol->GetVolName());
    // 		string vol_name;
    // 		p.FindBaseName().GetString(vol_name);
    // 		string vol_lmdb_full_name;
    // 		CSeqDB_Substring fn(vol_lmdb_name);
    // 		SeqDB_CombinePath(p.FindDirName(),fn, NULL, vol_lmdb_full_name);
    // 		if( 0 == i ) {
    // 			lmdb_full_name = vol_lmdb_full_name;
    // 		}
    //
    // 		if(vol_lmdb_full_name == lmdb_full_name) {
    // 			lmdb_vol_names.push_back(vol_name);
    // 			continue;
    // 		}
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:44-71
    // ```c++
    // 	m_LMDB->GetVolumesInfo(tmp_name, tmp_num_oids);
    // 	m_VolInfo.resize(tmp_name.size());
    // 	vector<string>::const_iterator itr = vol_names.begin();
    // 	if(vol_names.size()  > tmp_name.size()) {
    // 		NCBI_THROW(CSeqDBException, eArgErr, "Input db vol size does not match lmdb vol size");
    // 	}
    // 	// The input vol_names should be sorted in the same order as
    // 	for(unsigned int i=0; i < tmp_name.size(); i++) {
    // 		SVolumeInfo & v = m_VolInfo[i];
    // 		v.vol_name = tmp_name[i];
    // 		total_num_oids_lmdb += tmp_num_oids[i];
    // 		v.max_oid = total_num_oids_lmdb;
    // 		if((itr != vol_names.end()) && (*itr == v.vol_name)) {
    // 			v.skipped_oids = 0;
    // 			m_OIDEnd  += tmp_num_oids[i];
    // 			itr++;
    // 		}
    // 		else {
    // 			v.skipped_oids = tmp_num_oids[i];
    // 		}
    // 	}
    //
    // 	if(m_OIDEnd == 0) {
    // 		NCBI_THROW(CSeqDBException, eArgErr, "Input db vol does not match lmdb vol");
    // 	}
    // 	if(m_OIDEnd != total_num_oids_lmdb) {
    // 		m_isPartial = true;
    // 	}
    // ```
    // NCBI opens these entries with the database, LOSAT only for a taxid list. A volume of
    // the database that the LMDB file does not list, or lists with another number of OIDs,
    // leaves NCBI's OIDs off by the difference; LOSAT refuses such a database.
    /// The LMDB files of a version 5 database with the volumes each one lists.
    fn lmdb_entries(&self) -> Result<Vec<LmdbEntry>, SeqDbError> {
        let mut groups: Vec<(PathBuf, std::ops::Range<usize>)> = Vec::new();
        for (at, volume) in self.volumes.iter().enumerate() {
            let name = volume.index().lmdb_file.as_deref().unwrap_or_default();
            let lmdb = volume
                .base()
                .with_file_name(String::from_utf8_lossy(name).as_ref());
            match groups.last_mut() {
                Some((path, range)) if *path == lmdb => range.end = at + 1,
                _ => groups.push((lmdb, at..at + 1)),
            }
        }
        let mut entries = Vec::with_capacity(groups.len());
        for (lmdb, range) in groups {
            let file = DbFile::open(&lmdb)?;
            let env = Env::open(&file).map_err(|message| SeqDbError::Format {
                path: lmdb.clone(),
                message,
            })?;
            let listed = tax_files::volumes(&env).map_err(|message| SeqDbError::Format {
                path: lmdb.clone(),
                message,
            })?;
            let held = &self.volumes[range.clone()];
            if held.len() > listed.len() {
                return Err(SeqDbError::NotFound(VOLUME_COUNT_ERROR.to_string()));
            }
            let mut next = 0;
            let mut volumes = Vec::with_capacity(listed.len());
            let (mut total, mut included) = (0u64, 0u64);
            for (name, count) in &listed {
                total += u64::from(*count);
                let holds = held.get(next).is_some_and(|volume| {
                    volume
                        .base()
                        .file_name()
                        .is_some_and(|base| base.as_encoded_bytes() == name.as_slice())
                });
                if holds {
                    if held[next].num_oids() != *count {
                        return Err(SeqDbError::Unsupported {
                            path: lmdb.clone(),
                            message: format!(
                                "the LMDB file gives volume {} {count} sequences, its index {}",
                                String::from_utf8_lossy(name),
                                held[next].num_oids()
                            ),
                        });
                    }
                    included += u64::from(*count);
                    next += 1;
                }
                let max_oid = u32::try_from(total).map_err(|_| SeqDbError::Unsupported {
                    path: lmdb.clone(),
                    message: "the LMDB file lists more sequences than a signed 32-bit count"
                        .to_string(),
                })?;
                volumes.push((max_oid, if holds { 0 } else { *count }));
            }
            if included == 0 {
                return Err(SeqDbError::NotFound(VOLUME_NAME_ERROR.to_string()));
            }
            if next < held.len() {
                return Err(SeqDbError::Unsupported {
                    path: lmdb.clone(),
                    message: format!(
                        "the LMDB file does not list volume {}",
                        held[next].base().display()
                    ),
                });
            }
            entries.push(LmdbEntry {
                lmdb,
                oid_start: self.volume_starts[range.start],
                partial: included != total,
                num_oids: total as u32,
                volumes,
            });
        }
        Ok(entries)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:315-333
    // ```c++
    // void CSeqDBLMDBSet::TaxIdsToOids(set<TTaxId>& tax_ids, vector<blastdb::TOid>& rv) const
    // {
    // 	vector<TTaxId> tax_ids_found;
    // 	set<TTaxId> rv_tax_ids;
    // 	m_LMDBEntrySet[0]->TaxIdsToOids(tax_ids, rv, tax_ids_found);
    // 	rv_tax_ids.insert(tax_ids_found.begin(), tax_ids_found.end());
    // 	for(unsigned int i=1; i < m_LMDBEntrySet.size(); i++) {
    // 		vector <TOid> tmp;
    // 		m_LMDBEntrySet[i]->TaxIdsToOids(tax_ids, tmp, tax_ids_found);
    // 		rv.insert(rv.end(), tmp.begin(), tmp.end());
    // 		if(rv_tax_ids.size() < tax_ids.size()) {
    // 			rv_tax_ids.insert(tax_ids_found.begin(), tax_ids_found.end());
    // 		}
    // 	}
    // 	if(rv.size() == 0) {
    // 		NCBI_THROW(CSeqDBException, eTaxidErr, "Taxonomy ID(s) not found in the " + GetDBFileNames() + " database" + (m_LMDBEntrySet.size()>1?"s":"")+ ".");
    // 	}
    // 	tax_ids.swap(rv_tax_ids);
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:349-351
    // ```c++
    // 	if(rv.size() == 0) {
    // 		LOG_POST(Warning << "Taxonomy ID(s) not found in the " + GetDBFileNames() + " database" + (m_LMDBEntrySet.size()>1?"s":"")+ ".");
    // 	}
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:109-123
    // ```c++
    //     try {
    //         if (gi_list || neg_list || m_Aliases.NeedTotalsScan(m_VolSet)) {
    //             m_NeedTotalsScan = true;
    //             x_InitIdSet();
    //         }
    //
    //         if ((! m_OidListSetup) && (oid_begin || oid_end)) {
    //             m_NeedTotalsScan = true;
    //         }
    //
    //         if (m_NeedTotalsScan) {
    //             CSeqDBLockHold locked(m_Atlas);
    //
    //             // This is a whole-database scan; it's always done in
    //             // approximate length mode.
    //             x_ScanTotals(true, & m_NumSeqs, & m_TotalLength,
    //                                & m_MaxLength, & m_MinLength, locked);
    // ```
    // A union of the found taxids is what the loop keeps: they are taxids of `tax_ids`.
    /// Restricts the database to the OIDs of a taxid list, as NCBI's `CSeqDB` opened with the
    /// list: from here on the ordinal ids are those of the included OIDs, the number of
    /// sequences and the total length are NCBI's scanned totals of them (`num_seqs`,
    /// `total_length`, `max_length`; the reports print `printed_num_seqs` and
    /// `printed_total_length`, the whole database's), and the definition lines are those of
    /// NCBI's filtered header read. Returns the warning NCBI writes (a negative list of
    /// taxids the database does not hold, which excludes nothing), to be written as it is.
    ///
    /// Errors: NCBI's (`SeqDbError::NotFound`): a version 4 database, a positive list that
    /// selects no OID, a missing lookup file; the reader's for damaged lookup files.
    pub fn restrict_taxids(&mut self, list: &TaxidList) -> Result<Option<String>, SeqDbError> {
        let (negative, tax_ids) = match list {
            TaxidList::Positive(tax_ids) => (false, tax_ids),
            TaxidList::Negative(tax_ids) => (true, tax_ids),
        };
        if !self.is_version5() {
            return Err(SeqDbError::NotFound(V4_TAXID_ERROR.to_string()));
        }
        let entries = self.lmdb_entries()?;
        let mut oids = Vec::new();
        let mut found = BTreeSet::new();
        for entry in &entries {
            let part = entry.tax_oids(negative, tax_ids)?;
            oids.extend(part.oids);
            found.extend(part.tax_ids_found);
        }
        let mut warning = None;
        if oids.is_empty() {
            let names: Vec<String> = entries
                .iter()
                .map(|entry| LmdbEntry::base_name(&entry.lmdb))
                .collect();
            let text = format!(
                "Taxonomy ID(s) not found in the {} database{}.",
                names.join(" "),
                if entries.len() > 1 { "s" } else { "" }
            );
            if !negative {
                return Err(SeqDbError::NotFound(text));
            }
            warning = Some(text);
        }
        let num_oids = self.volume_starts[self.volumes.len()];
        let fill = if negative { u64::MAX } else { 0 };
        let mut words = vec![fill; (num_oids as usize).div_ceil(64)];
        for oid in oids {
            if oid >= num_oids {
                return Err(SeqDbError::Format {
                    path: entries[0].lmdb.clone(),
                    message: format!(
                        "the taxid lookup gives OID {oid} of a database of {num_oids} sequences"
                    ),
                });
            }
            let (word, bit) = ((oid / 64) as usize, oid % 64);
            if negative {
                words[word] &= !(1u64 << bit);
            } else {
                words[word] |= 1u64 << bit;
            }
        }
        let filter = OidFilter::new(words, num_oids);
        let (total_length, max_length) = self.scan_totals(&filter)?;
        let rule = if negative {
            (!found.is_empty()).then(|| DeflineRule::NegativeTaxids(found))
        } else {
            Some(DeflineRule::Taxids(found))
        };
        if let Some(rule) = rule.map(Arc::new) {
            for volume in &mut self.volumes {
                volume.set_defline_rule(Arc::clone(&rule));
            }
        }
        self.num_seqs = filter.count();
        self.total_length = total_length;
        self.max_length = max_length;
        self.filter = Some(Box::new(filter));
        Ok(warning)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1569-1590
    // ```c++
    //     for(int oid = 0; x_CheckOrFindOID(oid, locked); oid++) {
    //         ++ oid_count;
    //
    //         int vol_oid = 0;
    //
    //         volp = m_VolSet.FindVol(oid, vol_oid);
    //
    //         _ASSERT(volp);
    //
    //         if (totlen || maxlen || minlen) {
    //             int len;
    //             if ('p' == m_SeqType) {
    //                 len = volp->GetSeqLengthProt(vol_oid);
    //             } else {
    //                 if (approx) {
    //                     len = volp->GetSeqLengthApprox(vol_oid);
    //                 } else {
    //                     len = volp->GetSeqLengthExact(vol_oid);
    //                 }
    //             }
    //             max_count = max(len, max_count);
    //             min_count = min(len, min_count);
    //             base_count += len;
    //         }
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:325-344
    // ```c++
    // int CSeqDBVol::GetSeqLengthApprox(int oid) const
    // {
    //     TIndx start_offset = 0;
    //     TIndx end_offset   = 0;
    //
    //     //m_Atlas.Lock(locked);
    //     m_Idx->GetSeqStartEnd(oid, start_offset, end_offset);
    //
    //     _ASSERT(m_Idx->GetSeqType() == 'n');
    //
    //     int whole_bytes = int(end_offset - start_offset - 1);
    //
    //     // Same principle as below - but use lower bits of oid
    //     // instead of fetching the actual last byte.  this should
    //     // correct for bias, unless sequence length modulo 4 has a
    //     // significant statistical bias, which seems unlikely to
    //     // me.
    //
    //     return (whole_bytes * 4) + (oid & 0x03);
    // }
    // ```
    // A protein's length is exact (`GetSeqLengthProt`: the bytes up to the next sequence
    // less the NUL between them); a nucleotide's is approximate, from its bytes and the
    // volume-relative OID, not from the count in its last byte.
    /// NCBI's totals of the included OIDs: the total length and the longest length.
    fn scan_totals(&self, filter: &OidFilter) -> Result<(u64, u32), SeqDbError> {
        let (mut total, mut longest) = (0u64, 0u32);
        for (at, volume) in self.volumes.iter().enumerate() {
            let start = self.volume_starts[at];
            let mut first = 0;
            while first < volume.num_oids() {
                let end = volume.num_oids().min(first.saturating_add(SCAN_WINDOW));
                if filter.rank(start + end) == filter.rank(start + first) {
                    first = end;
                    continue;
                }
                let index = volume.index();
                let seq = index.sequence_offsets(first, end)?;
                let amb = index.ambiguity_offsets(first, end)?;
                for i in 0..(end - first) as usize {
                    let oid = first + i as u32;
                    if !filter.contains(start + oid) {
                        continue;
                    }
                    let length = match self.seqtype {
                        SeqType::Protein => i64::from(seq[i + 1]) - i64::from(seq[i]) - 1,
                        SeqType::Nucleotide => {
                            4 * (i64::from(amb[i]) - i64::from(seq[i]) - 1) + i64::from(oid & 3)
                        }
                    };
                    let length = u32::try_from(length.max(0)).unwrap_or(u32::MAX);
                    total += u64::from(length);
                    longest = longest.max(length);
                }
                first = end;
            }
        }
        Ok((total, longest))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn filter_of(num_oids: u32, included: &[u32]) -> OidFilter {
        let mut words = vec![0u64; (num_oids as usize).div_ceil(64)];
        for &oid in included {
            words[(oid / 64) as usize] |= 1 << (oid % 64);
        }
        OidFilter::new(words, num_oids)
    }

    /// A deterministic pseudo-random sequence (xorshift).
    fn next(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    #[test]
    fn rank_select_and_runs_agree_with_a_list() {
        let mut state = 0x2026_1030_u64;
        for &(num_oids, per_mille) in &[
            (1u32, 1000u64),
            (63, 500),
            (64, 1000),
            (65, 0),
            (511, 900),
            (512, 10),
            (513, 500),
            (5000, 2),
            (5000, 998),
            (20_000, 300),
        ] {
            let included: Vec<u32> = (0..num_oids)
                .filter(|_| next(&mut state) % 1000 < per_mille)
                .collect();
            let filter = filter_of(num_oids, &included);
            assert_eq!(filter.count() as usize, included.len(), "{num_oids}");
            for (ordinal, &oid) in included.iter().enumerate() {
                assert_eq!(
                    filter.select(ordinal as u32),
                    oid,
                    "{num_oids} select {ordinal}"
                );
            }
            for oid in 0..=num_oids {
                let below = included.partition_point(|&o| o < oid) as u32;
                assert_eq!(filter.rank(oid), below, "{num_oids} rank {oid}");
                assert_eq!(filter.contains(oid), included.binary_search(&oid).is_ok());
                let limit = num_oids.min(oid + 130);
                let mut end = oid;
                while end < limit && included.binary_search(&end).is_ok() {
                    end += 1;
                }
                assert_eq!(
                    filter.run_end(oid, limit),
                    end.max(oid.min(limit)),
                    "{num_oids} run {oid}"
                );
            }
        }
    }

    #[test]
    fn bits_past_the_last_oid_are_not_included() {
        let filter = OidFilter::new(vec![u64::MAX; 2], 70);
        assert_eq!(filter.count(), 70);
        assert!(!filter.contains(70));
        assert_eq!(filter.run_end(60, 128), 70);
        assert_eq!(filter.select(69), 69);
    }

    fn tx_dir() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/lmdb")
    }

    /// NCBI's OIDs of `tx` for each taxid set (`blastdbcmd -taxids S -no_taxid_expansion` for a
    /// positive set, the sequences `blastp -negative_taxids S -no_taxid_expansion` searched for a
    /// negative one: `tests/seqdb/lmdb/expected/tx_oracle.txt`).
    #[test]
    fn restricted_oids_are_ncbis_on_the_tx_fixture() {
        let oracle = std::fs::read_to_string(tx_dir().join("expected/tx_oracle.txt")).unwrap();
        let name = tx_dir().join("tx");
        let mut checked = 0;
        for line in oracle.lines().filter(|line| !line.trim().is_empty()) {
            let mut fields = line.split(' ').filter(|field| !field.is_empty());
            let kind = fields.next().unwrap();
            let set: BTreeSet<i32> = fields
                .next()
                .unwrap()
                .split(',')
                .map(|taxid| taxid.parse().unwrap())
                .collect();
            let rest: Vec<&str> = fields.collect();
            let list = match kind {
                "P" => TaxidList::Positive(set),
                "N" => TaxidList::Negative(set),
                other => panic!("{other}"),
            };
            let mut db = SeqDb::open(name.to_str().unwrap(), SeqType::Protein).unwrap();
            let whole = (db.num_seqs(), db.total_length());
            let result = db.restrict_taxids(&list);
            if rest.first() == Some(&"Error:") {
                let Err(SeqDbError::NotFound(text)) = result else {
                    panic!("{line}: {result:?}");
                };
                assert_eq!(
                    text, "Taxonomy ID(s) not found in the tx database.",
                    "{line}"
                );
                continue;
            }
            let warning = result.unwrap_or_else(|error| panic!("{line}: {error}"));
            let want: Vec<u32> = rest.iter().map(|oid| oid.parse().unwrap()).collect();
            let got: Vec<u32> = (0..db.num_seqs()).map(|ordinal| db.oid(ordinal)).collect();
            assert_eq!(got, want, "{line}");
            assert_eq!(
                (db.printed_num_seqs(), db.printed_total_length()),
                whole,
                "{line}"
            );
            // A negative list of taxids that `tx` does not hold excludes nothing and warns.
            assert_eq!(
                warning.is_some(),
                matches!(&list, TaxidList::Negative(set) if set.iter().all(|t| *t == 999)),
                "{line}"
            );
            // The sequences read through the ordinal ids are those of the OIDs.
            let mut whole_db = SeqDb::open(name.to_str().unwrap(), SeqType::Protein).unwrap();
            whole_db.filter = None;
            let mut lengths = 0u64;
            for (ordinal, &oid) in want.iter().enumerate() {
                let filtered = db.read_sequence(ordinal as u32).unwrap();
                let plain = whole_db.read_sequence(oid).unwrap();
                assert_eq!(filtered.protein(0), plain.protein(0), "{line} {oid}");
                lengths += filtered.protein(0).len() as u64;
            }
            assert_eq!(db.total_length(), lengths, "{line}");
            let mut first = 0;
            let mut seen = Vec::new();
            while let Some(range) = db.next_block(first, 1 << 20) {
                for at in range.first..range.end {
                    seen.push(db.volume_starts[range.volume] + at);
                }
                first = range.end_ordinal();
            }
            assert_eq!(seen, want, "{line}: blocks");
            checked += 1;
        }
        assert!(checked >= 30, "{checked}");
    }

    #[test]
    fn a_version_4_database_refuses_a_taxid_list() {
        let name = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/db/prot_v4");
        let mut db = SeqDb::open(name.to_str().unwrap(), SeqType::Protein).unwrap();
        let error = db
            .restrict_taxids(&TaxidList::Positive([9606].into_iter().collect()))
            .unwrap_err();
        assert_eq!(
            error.ncbi_message().as_deref(),
            Some("Taxonomy filtering is not supported in v4 BLAST dbs")
        );
    }

    fn defline(taxid: Option<i32>, links: Option<Vec<i32>>) -> BlastDefline {
        BlastDefline {
            taxid,
            links,
            ..BlastDefline::default()
        }
    }

    #[test]
    fn definition_lines_follow_ncbis_rules() {
        let set = |ids: &[i32]| ids.iter().copied().collect::<BTreeSet<i32>>();
        let positive = DeflineRule::Taxids(set(&[9606, 10090]));
        assert!(positive.keeps(&defline(Some(9606), None)));
        assert!(positive.keeps(&defline(Some(1), Some(vec![10090]))));
        assert!(!positive.keeps(&defline(Some(1), Some(vec![2]))));
        // `IsSetTaxid`: an absent taxid is not 0.
        assert!(!DeflineRule::Taxids(set(&[0])).keeps(&defline(None, None)));
        assert!(DeflineRule::Taxids(set(&[0])).keeps(&defline(Some(0), None)));

        let negative = DeflineRule::NegativeTaxids(set(&[9606]));
        assert!(!negative.keeps(&defline(Some(9606), None)));
        assert!(negative.keeps(&defline(Some(9606), Some(vec![10090]))));
        // `GetTaxIds` drops 0 beside another taxid.
        assert!(!negative.keeps(&defline(Some(0), Some(vec![9606]))));
        assert!(negative.keeps(&defline(Some(0), None)));
        // No taxid at all: an empty set is a subset of every list.
        assert!(!negative.keeps(&defline(None, None)));
    }

    #[test]
    fn a_rule_that_keeps_no_line_keeps_every_line() {
        let rule = DeflineRule::Taxids([5].into_iter().collect());
        let set = DeflineSet {
            deflines: vec![defline(Some(1), None), defline(Some(5), None)],
        };
        assert_eq!(rule.keep_mask(&set), vec![false, true]);
        let none = DeflineSet {
            deflines: vec![defline(Some(1), None), defline(Some(2), None)],
        };
        assert_eq!(rule.keep_mask(&none), vec![true, true]);
        let mut lines = vec![1, 2, 3];
        lines.retain_lines(&[true, false, true]);
        assert_eq!(lines, vec![1, 3]);
        lines.retain_lines(&[false]);
        assert_eq!(lines, vec![1, 3]);
    }
}
