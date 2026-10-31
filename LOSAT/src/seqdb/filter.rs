//! A database seen through a list of taxids (line DI-8a: `-taxids`, `-taxidlist`,
//! `-negative_taxids`, `-negative_taxidlist`) or of seqids (DI-8b: `-seqidlist`,
//! `-negative_seqidlist`): the OIDs that NCBI's `CSeqDB` visits when it is opened with a
//! user list or a negative list (`CSeqDBOIDList`), the totals it scans for them
//! (`CSeqDBImpl::x_ScanTotals`) and the definition lines its filtered header read keeps
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
use super::isam;
use super::lmdb_read::Env;
use super::seqidlist::{self, SeqidListFile};
use super::{tax_files, SeqDb, SeqDbError, SeqType};
use crate::seqid::model::ObjectId;
use crate::seqid::{label, BlastDefline, DeflineSet, ParsedSeqId};

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

/// The seqid list of a search as read from its file (`CSeqDBFileGiList(fname, eSiList)`,
/// `seqidlist::read_si_list`): `-seqidlist` gives a positive list, `-negative_seqidlist`
/// a negative one. A file without an id is no list at all (`CSearchDatabase::x_InitializeDb`
/// opens the database without one); [`SeqDb::restrict`] is not called for it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeqidList {
    pub negative: bool,
    pub file: SeqidListFile,
}

/// The list a search restricts its database to: NCBI reads one of them
/// (`CBlastDatabaseArgs::ExtractAlgorithmOptions`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DbList {
    Taxids(TaxidList),
    Seqids(SeqidList),
}

/// A warning NCBI writes while it applies the list, to be written as it is.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ListWarning {
    /// `LOG_POST(Warning << text)` or a line written to `cerr`: the text alone.
    Log(String),
    /// `ERR_POST(Warning << text)`: with the diagnostic prefix of the program.
    Posted(String),
    /// `ERR_POST(text)` at the default severity (an error that is reported before NCBI
    /// throws): with the diagnostic prefix of the program.
    Error(String),
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:206-226
// ```c++
// bool s_VerifySeqidlist(const SBlastSeqIdListInfo & list_info, const CSeqDBVolSet & volset, const CSeqDBLMDBSet & lmdb_set)
// {
// 	if(list_info.is_v4 && lmdb_set.IsBlastDBVersion5()) {
// 		ERR_POST(Warning << "To obtain better run time performance, " \
// 				            "please run blastdb_aliastool -seqid_file_in <INPUT_FILE_NAME> " \
// 				            "-seqid_file_out <OUT_FILE_NAME> and use <OUT_FILE_NAME> as the argument to -seqidlist");
// 		return true;
// 	}
//
// 	if((!list_info.is_v4) && (!lmdb_set.IsBlastDBVersion5())) {
// 		NCBI_THROW(CSeqDBException, eArgErr,
// 		           "Seqidlist is not in BLAST db v4 format");
// 	}
//
// 	if((list_info.db_vol_length > 0) &&
// 	   (list_info.db_vol_length != volset.GetVolumeSetLength())) {
// 		ERR_POST(Warning << "Seqidlist file db info does not match input db");
// 	}
//
// 	 return false;
// }
// ```
const TEXT_LIST_WARNING: &str = "To obtain better run time performance, please run \
     blastdb_aliastool -seqid_file_in <INPUT_FILE_NAME> -seqid_file_out <OUT_FILE_NAME> and use \
     <OUT_FILE_NAME> as the argument to -seqidlist";
const BINARY_LIST_V4_ERROR: &str = "Seqidlist is not in BLAST db v4 format";
const DB_INFO_WARNING: &str = "Seqidlist file db info does not match input db";
const EMPTY_LIST_WARNING: &str = "Empty seqid list";

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
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbvol.hpp:943-962
// ```c++
//     void x_FilterHasId(const CSeq_id & id,
//                        bool          & have_user,
//                        bool          & have_vol) const
//     {
//         if (! have_user) {
//             if (m_UserGiList.NotEmpty() && m_UserGiList->HasIdFilters()) {
//                 have_user |= x_ListIncludesId(*m_UserGiList, id);
//             } else if (m_NegativeList.NotEmpty() && m_NegativeList->HasIdFilters() ) {
//                 have_user |= x_ListIncludesId(*m_NegativeList, id);
//             } else {
//                 have_user = true;
//             }
//         }
//
//         if (! have_vol) {
//             if (m_VolumeGiLists.empty()) {
//                 have_vol = true;
//             } else {
// ```
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbvol.hpp:997-1008
// ```c++
//     bool x_ListIncludesId(CSeqDBNegativeList & L, const CSeq_id & id) const
//     {
//         // A defline is included IFF either a GI or TI is found, and
//         // that ID is not on the list.
//
//         // I use the terms 'included' and 'mentioned' to describe the
//         // negative list processing as follows: "A negative list
//         // INCLUDES a TI or GI if that ID is not MENTIONED in the
//         // negative list."
//
//         bool match_type = false;
//         bool found = L.FindId(id, match_type);
//
//         return (! found) && match_type;
// ```
// A list of seqids holds no gi and no ti (`HasIdFilters` is its `GetNumSis() > 0`), and the
// volumes have no lists of their own (alias `GILIST` keys are refused, DI-8 DECISIONS 6).
/// Which definition lines of an included OID a restricted database keeps.
#[derive(Debug)]
pub(crate) enum DeflineRule {
    /// A positive list: the lines whose `taxid` or `links` hold one of these taxids.
    Taxids(BTreeSet<i32>),
    /// A negative list: the lines with a taxid (`GetTaxIds`) outside these.
    NegativeTaxids(BTreeSet<i32>),
    /// A positive seqid list (sorted, unique, as looked up): the lines with a Seq-id in it
    /// (`CSeqDBGiList::FindId`).
    Seqids(Vec<Vec<u8>>),
    /// A negative seqid list (sorted, unique, as looked up): the lines with a Seq-id that is
    /// not a gi or a ti and not in it (`CSeqDBNegativeList::FindId`).
    NegativeSeqids(Vec<Vec<u8>>),
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
            DeflineRule::Seqids(list) => defline.ids.iter().any(|id| find_id(list, id)),
            DeflineRule::NegativeSeqids(list) => defline.ids.iter().any(|id| {
                let (found, match_type) = negative_find_id(list, id);
                !found && match_type
            }),
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
    // returns null, seqdbvol.cpp:1370-1376): its reports leave the sequence's hits out
    // without a word, in every format (dropin-survey-20261029/DI8 probe R01-R07, R14, R21:
    // a seqid list naming a PDB molecule without its chain, a local id in other letters, or
    // a Swiss-Prot name on a version 4 database selects the OID through the lookup table
    // but matches none of its lines). makeblastdb writes the taxid lookup files from the
    // same lines, so an OID of a taxid list has a line of the list; LOSAT keeps every line
    // of such an OID that would have none, and refuses a seqid list that leaves an OID
    // without a line (DI-8b DECISIONS 17).
    /// Which of the definition lines of `set` the read keeps, in order: every one when a
    /// taxid rule would keep none, `None` when a seqid rule would keep none.
    pub(crate) fn keep_mask(&self, set: &DeflineSet) -> Option<Vec<bool>> {
        let mask: Vec<bool> = set.deflines.iter().map(|d| self.keeps(d)).collect();
        if mask.iter().any(|&keep| keep) {
            Some(mask)
        } else if matches!(
            self,
            DeflineRule::Seqids(_) | DeflineRule::NegativeSeqids(_)
        ) {
            None
        } else {
            Some(vec![true; mask.len()])
        }
    }
}

/// LOSAT's refusal of an included OID whose every definition line a seqid list removes.
pub(crate) fn no_line_kept(oid: u32) -> String {
    format!(
        "the seqid list keeps sequence {oid} of the volume but none of its definition lines, \
         and NCBI leaves its hits out of the report without notice"
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:733-737
// ```c++
// bool CSeqDBGiList::FindSi(const string &si) const
// {
//     int oid(0), index(0);
//     return (const_cast<CSeqDBGiList *>(this))->SiToOid(si, oid, index);
// }
// ```
// `SiToOid` (:745-767) and `CSeqDBNegativeList::FindSi` (:1542-1560) are binary searches of
// the sorted ids by `std::string` order (bytes).
/// Whether the sorted list holds `si`.
fn find_si(list: &[Vec<u8>], si: &[u8]) -> bool {
    list.binary_search_by(|entry| entry.as_slice().cmp(si))
        .is_ok()
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2698-2705
// ```c++
// string GetBlastSeqIdString(const CSeq_id & seqid, bool version)
// {
// 	if(seqid.IsPir() || seqid.IsPrf()) {
// 		return seqid.AsFastaString();
// 	}
//
// 	return seqid.GetSeqIdString(version);
// }
// ```
fn blast_seq_id_string(id: &ParsedSeqId, version: bool) -> Vec<u8> {
    let mut out = Vec::new();
    match id {
        ParsedSeqId::Pir(_) | ParsedSeqId::Prf(_) => label::fasta(id, &mut out),
        _ => label::content(id, version, &mut out),
    }
    out
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2264-2400
// ```c++
// ESeqDBIdType SeqDB_SimplifySeqid(CSeq_id       & bestid,
//                                  const string  * acc,
//                                  Int8          & num_id,
//                                  string        & str_id,
//                                  bool          & simpler)
// {
//     ESeqDBIdType result = eStringId;
//
//     const CTextseq_id * tsip = 0;
//
//     bool matched = true;
//
//     switch(bestid.Which()) {
//     case CSeq_id::e_Gi:
//         simpler = true;
//         num_id = GI_TO(Int8, bestid.GetGi());
//         result = eGiId;
//         break;
//
//     case CSeq_id::e_Gibbsq:    /* gibbseq */
//         simpler = true;
//         result = eStringId;
//         str_id = NStr::UIntToString(bestid.GetGibbsq());
//         break;
//
//     case CSeq_id::e_General:
//         {
//             const CDbtag & dbt = bestid.GetGeneral();
//
//             if (dbt.CanGetDb()) {
//                 if (dbt.GetDb() == "BL_ORD_ID") {
//                     simpler = true;
//                     num_id = dbt.GetTag().GetId();
//                     result = eOID;
//                     break;
//                 }
//
//                 if (dbt.GetDb() == "PIG") {
//                     simpler = true;
//                     num_id = dbt.GetTag().GetId();
//                     result = ePigId;
//                     break;
//                 }
//
//                 if (dbt.GetDb() == "ti") {
//                     simpler = true;
//                     num_id = (dbt.GetTag().IsStr()
//                               ? NStr::StringToInt8(dbt.GetTag().GetStr())
//                               : dbt.GetTag().GetId());
//
//                     result = eTiId;
//                     break;
//                 }
//
//
//                 if (NStr::CompareNocase(dbt.GetDb(), "GNOMON") == 0) {
//                     str_id = bestid.AsFastaString();
//                     str_id = NStr::ToLower(str_id);
//                     result = eStringId;
//                     break;
//                 }
//             }
//
//             if (dbt.CanGetTag() && dbt.GetTag().IsStr()) {
//                 result = eStringId;
//                 str_id = dbt.GetTag().GetStr();
//                 str_id = NStr::ToLower(str_id);
//             } else {
//                 // Use the default logic.
//                 matched = false;
//             }
//         }
//         break;
//
//     case CSeq_id::e_Local:     /* local */
//         simpler = true;
//         result = eStringId;
//         {
//             const CObject_id & objid = bestid.GetLocal();
//
//             if (objid.IsStr()) {
//                 // sparse version will leave "lcl|" off.
//                 str_id = objid.GetStr();
//                 str_id = NStr::ToLower(str_id);
//             } else {
//                 // Local numeric ids are stored as strings.
//                 str_id = "lcl|" + NStr::IntToString(objid.GetId());
//             }
//         }
//         break;
//
//         // tsip types
//
//     case CSeq_id::e_Embl:      /* embl */
//     case CSeq_id::e_Ddbj:      /* ddbj */
//     case CSeq_id::e_Genbank:   /* genbank */
//     case CSeq_id::e_Tpg:       /* Third Party Annot/Seq Genbank */
//     case CSeq_id::e_Tpe:       /* Third Party Annot/Seq EMBL */
//     case CSeq_id::e_Tpd:       /* Third Party Annot/Seq DDBJ */
//     case CSeq_id::e_Other:     /* other */
//     case CSeq_id::e_Swissprot: /* swissprot (now with versions) */
//     case CSeq_id::e_Gpipe:     /* internal NCBI genome pipeline */
//         tsip = bestid.GetTextseq_Id();
//         break;
//
//     case CSeq_id::e_Pir:       /* pir   */
//     case CSeq_id::e_Prf:       /* prf   */
//         tsip = bestid.GetTextseq_Id();
//         break;
//
//     default:
//         matched = false;
//     }
//
//     // Default: if we have a string, use it; if we only have seqid,
//     // create a string.  This should not happen if the seqid matches
//     // one of the cases above, which currently correspond to all the
//     // supported seqid types.
//
//     CSeq_id::ELabelFlags label_flags = (CSeq_id::ELabelFlags)
//         (CSeq_id::fLabel_GeneralDbIsContent | CSeq_id::fLabel_Version);
//
//     if (! matched) {
//         // (should not happen normally)
//
//         simpler = false;
//         result  = eStringId;
//
//         if (acc) {
//             str_id = *acc;
//             str_id = NStr::ToLower(str_id);
//         } else {
//             bestid.GetLabel(& str_id, CSeq_id::eFasta, label_flags);
//             str_id = NStr::ToLower(str_id);
//         }
//     }
//
//     if (tsip) {
//         bool found = false;
//
//         if (tsip->CanGetAccession()) {
//             str_id = tsip->GetAccession();
//             str_id = NStr::ToLower(str_id);
//             found = true;
//
//             if (tsip->CanGetVersion()) {
//                 str_id += ".";
//                 str_id += NStr::UIntToString(tsip->GetVersion());
//             }
//         } else if (tsip->CanGetName()) {
//             str_id = tsip->GetName();
//             str_id = NStr::ToLower(str_id);
//             found = true;
//         }
//
//         if (found) {
//             simpler = true;
// ```
// The list lookups call it with no accession string (`acc` null) and read only `str_id`:
// empty when the id is a number (a gi, a `BL_ORD_ID`, `PIG` or `ti` general id). A `ti`
// general id with a non-numeric tag throws in NCBI; it gives the empty string here.
/// `str_id` of `SeqDB_SimplifySeqid(id, NULL, ...)`, lower-cased as NCBI leaves it.
fn simplify_seqid_string(id: &ParsedSeqId) -> Vec<u8> {
    let fasta_lower = |id: &ParsedSeqId| {
        let mut out = Vec::new();
        label::fasta(id, &mut out);
        out.make_ascii_lowercase();
        out
    };
    match id {
        ParsedSeqId::Gi(_) => Vec::new(),
        ParsedSeqId::Gibbsq(value) => value.to_string().into_bytes(),
        ParsedSeqId::General(dbtag) => {
            if dbtag.db == b"BL_ORD_ID" || dbtag.db == b"PIG" || dbtag.db == b"ti" {
                return Vec::new();
            }
            if dbtag.db.eq_ignore_ascii_case(b"GNOMON") {
                return fasta_lower(id);
            }
            match &dbtag.tag {
                ObjectId::Str(tag) => tag.to_ascii_lowercase(),
                ObjectId::Id(_) => fasta_lower(id),
            }
        }
        ParsedSeqId::Local(ObjectId::Str(tag)) => tag.to_ascii_lowercase(),
        ParsedSeqId::Local(ObjectId::Id(value)) => format!("lcl|{value}").into_bytes(),
        ParsedSeqId::Embl(text)
        | ParsedSeqId::Ddbj(text)
        | ParsedSeqId::Genbank(text)
        | ParsedSeqId::Tpg(text)
        | ParsedSeqId::Tpe(text)
        | ParsedSeqId::Tpd(text)
        | ParsedSeqId::Other(text)
        | ParsedSeqId::Swissprot(text)
        | ParsedSeqId::Gpipe(text)
        | ParsedSeqId::Pir(text)
        | ParsedSeqId::Prf(text) => {
            if let Some(accession) = &text.accession {
                let mut str_id = accession.to_ascii_lowercase();
                if let Some(version) = text.version {
                    str_id.push(b'.');
                    str_id.extend_from_slice(version.to_string().as_bytes());
                }
                str_id
            } else if let Some(name) = &text.name {
                name.to_ascii_lowercase()
            } else {
                Vec::new()
            }
        }
        _ => fasta_lower(id),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1635-1665
// ```c++
// bool CSeqDBGiList::FindId(const CSeq_id & id)
// {
//     if (id.IsGi()) {
//         return FindGi(id.GetGi());
//     } else if (id.IsGeneral() && id.GetGeneral().GetDb() == "ti") {
//         const CObject_id & obj = id.GetGeneral().GetTag();
//
//         TTi ti = (obj.IsId()
//                    ? (TTi) obj.GetId()
//                    : (TTi) NStr::StringToInt8(obj.GetStr()));
//
//         return FindTi(ti);
//     } else {
//     	if(FindSi(GetBlastSeqIdString(id, true))) return true;
//     	if(FindSi(GetBlastSeqIdString(id, false))) return true;
//
//     	/// For isam lookup
//         Int8 num_id;
//         string str_id;
//         bool simpler;
//         SeqDB_SimplifySeqid(*(const_cast<CSeq_id *>(&id)), 0, num_id, str_id, simpler);
//         if (FindSi(str_id)) return true;
//
//         // We may have to strip the version to find it...
//         size_t pos = str_id.find(".");
//         if (pos != str_id.npos) {
//             string nover(str_id, 0, pos);
//             return FindSi(nover);
//         }
//     }
//     return false;
// }
// ```
// A seqid list holds no gi and no ti: `FindGi` and `FindTi` find nothing.
/// `CSeqDBGiList::FindId` of a list of seqids.
fn find_id(list: &[Vec<u8>], id: &ParsedSeqId) -> bool {
    match id {
        ParsedSeqId::Gi(_) => false,
        ParsedSeqId::General(dbtag) if dbtag.db == b"ti" => false,
        _ => find_si_forms(list, id),
    }
}

/// The lookups of a Seq-id that is not a gi or a ti: its string with and without version,
/// then its simplified string with and without version.
fn find_si_forms(list: &[Vec<u8>], id: &ParsedSeqId) -> bool {
    if find_si(list, &blast_seq_id_string(id, true))
        || find_si(list, &blast_seq_id_string(id, false))
    {
        return true;
    }
    let str_id = simplify_seqid_string(id);
    if find_si(list, &str_id) {
        return true;
    }
    match str_id.iter().position(|&c| c == b'.') {
        Some(pos) => find_si(list, &str_id[..pos]),
        None => false,
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1565-1610
// ```c++
// bool CSeqDBNegativeList::FindId(const CSeq_id & id, bool & match_type)
// {
//     if (id.IsGi()) {
//         match_type = (GetNumGis() > 0) ? true : false;
//         if(match_type) {
//             return(FindGi(id.GetGi()));
//         }
//     } else if (id.IsGeneral() && id.GetGeneral().GetDb() == "ti") {
//         match_type = (GetNumTis() > 0) ? true : false;
//
//         if(match_type) {
//             const CObject_id & obj = id.GetGeneral().GetTag();
//
//             Int8 ti = (obj.IsId()
//                    ? obj.GetId()
//                    : NStr::StringToInt8(obj.GetStr()));
//
//             return FindTi(ti);
//         }
//     } else {
//         match_type = (GetNumSis() > 0) ? true : false;
//
//         if(match_type) {
//         	if(FindSi(GetBlastSeqIdString(id, true))) return true;
//         	if(FindSi(GetBlastSeqIdString(id, false))) return true;
//
//         	// For isam lookup
//             Int8 num_id;
//             string str_id;
//             bool simpler;
//
//             SeqDB_SimplifySeqid(*(const_cast<CSeq_id *>(&id)), 0, num_id, str_id, simpler);
//
//             if (FindSi(str_id)) {
//                 return true;
//             }
//
//             // We may have to strip the version to find it...
//             size_t pos = str_id.find(".");
//             if (pos != str_id.npos) {
//                 string nover(str_id, 0, pos);
//                 return FindSi(nover);
//             }
//         }
//     }
//     return false;
// }
// ```
/// `CSeqDBNegativeList::FindId` of a list of seqids: `(found, match_type)`.
fn negative_find_id(list: &[Vec<u8>], id: &ParsedSeqId) -> (bool, bool) {
    match id {
        ParsedSeqId::Gi(_) => (false, false),
        ParsedSeqId::General(dbtag) if dbtag.db == b"ti" => (false, false),
        _ => {
            let match_type = !list.is_empty();
            (match_type && find_si_forms(list, id), match_type)
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

    /// The LMDB file opened for a lookup, handed to `f`.
    fn with_env<T>(
        &self,
        f: impl FnOnce(&Env<DbFile>) -> Result<T, SeqDbError>,
    ) -> Result<T, SeqDbError> {
        let file = DbFile::open(&self.lmdb)?;
        let env = Env::open(&file).map_err(|message| SeqDbError::Format {
            path: self.lmdb.clone(),
            message,
        })?;
        f(&env)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:81-113
    // ```c++
    // void CSeqDBLMDBEntry::x_AdjustOidsOffset(vector<TOid> & oids) const
    // {
    // 	if ((m_OIDStart > 0) && (!m_isPartial)) {
    // 		for(unsigned int i=0; i < oids.size(); i++) {
    // 			if (oids[i] != kSeqDBEntryNotFound) {
    // 				oids[i] +=m_OIDStart;
    // 			}
    // 		}
    // 	}
    //
    // 	if(m_isPartial) {
    // 		for(unsigned int i=0; i < oids.size(); i++) {
    // 			if (oids[i] != kSeqDBEntryNotFound) {
    // 				TOid skipped_oids = 0;
    // 				for(unsigned int j=0;j < m_VolInfo.size(); j++) {
    // 					if(oids[i] >=  m_VolInfo[j].max_oid) {
    // 						skipped_oids += m_VolInfo[j].skipped_oids;
    // 					}
    // 					else {
    // 						if(m_VolInfo[j].skipped_oids > 0) {
    // 							// ids found but in skipped vol
    // 							oids[i] = kSeqDBEntryNotFound;
    // 						}
    // 						else {
    // 							oids[i] = oids[i] + m_OIDStart - skipped_oids ;
    // 						}
    // 						break;
    // 					}
    // 				}
    // 			}
    // 		}
    // 	}
    // }
    // ```
    /// The database OID of the LMDB file's OID `oid` of a seqid lookup: `None` when its
    /// volume is one the database does not hold.
    fn database_oid(&self, oid: u32) -> Option<u32> {
        if !self.partial {
            return Some(oid + self.oid_start);
        }
        let mut skipped = 0;
        for &(max_oid, skipped_oids) in &self.volumes {
            if oid >= max_oid {
                skipped += skipped_oids;
            } else {
                return (skipped_oids == 0).then(|| oid + self.oid_start - skipped);
            }
        }
        None
    }

    /// `CSeqDBLMDB::GetOids` of one id: the first OID of the key (`tax_files::acc_first_oid`),
    /// with NCBI's messages for a database without the table (`DB contains no accession
    /// info.`, the text of NCBI's LMDB manager for a `.pdb`/`.ndb` without `acc2oid`: a
    /// database made without `-parse_seqids`, probe S48) and for a damaged one.
    fn first_oid(&self, env: &Env<DbFile>, acc: &[u8]) -> Result<Option<u32>, SeqDbError> {
        tax_files::acc_first_oid(env, acc).map_err(|message| {
            let name = Self::base_name(&self.lmdb);
            if message == tax_files::NO_ACCESSION_INFO {
                SeqDbError::NotFound(message)
            } else {
                SeqDbError::Format {
                    path: self.lmdb.clone(),
                    message: format!("{} ({message})", tax_files::accessions_lookup_error(&name)),
                }
            }
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:153-157
    // ```c++
    // void CSeqDBLMDBEntry::AccessionsToOids(const vector<string>& accs, vector<TOid>& oids) const
    // {
    // 	m_LMDB->GetOids(accs, oids);
    // 	x_AdjustOidsOffset(oids);
    // }
    // ```
    // `GetOids` (seqdb_lmdb.cpp:360-400) reads the first duplicate of each key with
    // `MDB_SET` and `MDB_GET_CURRENT`: an accession of more than about 500 OIDs, whose
    // duplicates LMDB keeps in a sub-database, makes NCBI's read fail (`Accessions to Oids
    // lookup error`); LOSAT reads its first duplicate like any other key (DI-8 DECISIONS 3,
    // PD-LOSAT-NCBI-DEFECTS).
    /// The database OID of each of `accs` in this entry (`None`: not found).
    fn accessions_to_oids(&self, accs: &[Vec<u8>]) -> Result<Vec<Option<u32>>, SeqDbError> {
        self.with_env(|env| {
            accs.iter()
                .map(|acc| {
                    Ok(self
                        .first_oid(env, acc)?
                        .and_then(|oid| self.database_oid(oid)))
                })
                .collect()
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:402-413
    // ```c++
    // struct SOidSeqIdPair
    // {
    // 	SOidSeqIdPair(blastdb::TOid o, const string & i) : oid(o), id(i) {}
    // 	blastdb::TOid oid;
    // 	string id;
    // 	static bool cmp_oid(const SOidSeqIdPair & v, const SOidSeqIdPair & k) {
    // 		if(v.oid == k.oid) {
    // 			return (v.id < k.id);
    // 		}
    // 		return (v.oid < k.oid );
    // 	}
    // };
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:515-556
    // ```c++
    // void
    // CSeqDBLMDB::NegativeSeqIdsToOids(const vector<string>& ids, vector<blastdb::TOid>& rv) const
    // {
    // 	rv.clear();
    // 	vector<blastdb::TOid> oids;
    // 	GetOids(ids, oids);
    // 	vector<SOidSeqIdPair> pairs;
    // 	for (unsigned int i=0; i < ids.size(); i++) {
    // 		if(oids[i] == kSeqDBEntryNotFound) {
    // 			continue;
    // 		}
    // 		else {
    // 			SOidSeqIdPair p(oids[i], ids[i]);
    // 			pairs.push_back(p);
    // 		}
    // 	}
    //
    // 	if(pairs.size() == 0) {
    // 		return;
    // 	}
    //
    // 	sort (pairs.begin(), pairs.end(), SOidSeqIdPair::cmp_oid);
    //
    // 	CMemoryFile oid_file(m_Oid2SeqIdsFile);
    // 	CLookupSeqIds lookup(oid_file);
    // 	blastdb::TOid current_oid = 0;
    // 	unsigned int i=0;
    // 	while (i < pairs.size()) {
    // 		vector<string> file_idlist;
    // 		vector<string> input_idlist;
    // 		current_oid = pairs[i].oid;
    // 		lookup.GetSeqIdListForOid(current_oid, file_idlist);
    // 		while ((i < pairs.size()) && (current_oid == pairs[i].oid)) {
    // 			input_idlist.push_back(pairs[i].id);
    // 			i++;
    // 		}
    // 		if(s_CompareIdList(file_idlist, input_idlist)) {
    // 			rv.push_back(current_oid);
    // 		}
    //
    // 	}
    //
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:159-163
    // ```c++
    // void CSeqDBLMDBEntry::NegativeSeqIdsToOids(const vector<string>& ids, vector<blastdb::TOid>& rv) const
    // {
    // 	m_LMDB->NegativeSeqIdsToOids(ids, rv);
    // 	x_AdjustOidsOffset(rv);
    // }
    // ```
    // An excluded OID of a volume the database does not hold becomes `kSeqDBEntryNotFound`,
    // which NCBI then clears as a bit (`ClearBit(-1)`, undefined); LOSAT drops it.
    /// The database OIDs that the negative list `ids` excludes in this entry: those whose
    /// every id of the `.nos`/`.pos` file is covered by the listed ids that have them as
    /// first OID.
    fn negative_seqids_to_oids(&self, ids: &[Vec<u8>]) -> Result<Vec<u32>, SeqDbError> {
        let oids = self.with_env(|env| {
            ids.iter()
                .map(|id| self.first_oid(env, id))
                .collect::<Result<Vec<_>, _>>()
        })?;
        let mut pairs: Vec<(u32, &[u8])> = oids
            .iter()
            .zip(ids)
            .filter_map(|(oid, id)| oid.map(|oid| (oid, id.as_slice())))
            .collect();
        if pairs.is_empty() {
            return Ok(Vec::new());
        }
        pairs.sort();
        let nos_path = self.lookup_file("os");
        let nos_file = DbFile::open(&nos_path)?;
        let lookup =
            tax_files::OidSeqIds::open(&nos_file).map_err(|message| SeqDbError::Format {
                path: nos_path.clone(),
                message,
            })?;
        let mut excluded = Vec::new();
        let mut i = 0;
        while i < pairs.len() {
            let current_oid = pairs[i].0;
            let file_ids = lookup
                .seqids(current_oid)
                .map_err(|message| SeqDbError::Format {
                    path: nos_path.clone(),
                    message,
                })?;
            let mut input_ids = Vec::new();
            while i < pairs.len() && pairs[i].0 == current_oid {
                input_ids.push(pairs[i].1);
                i += 1;
            }
            if compare_id_list(&file_ids, &input_ids).map_err(|message| SeqDbError::Format {
                path: nos_path.clone(),
                message,
            })? {
                excluded.push(current_oid);
            }
        }
        Ok(excluded
            .into_iter()
            .filter_map(|oid| self.database_oid(oid))
            .collect())
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:461-513
// ```c++
// bool s_CompareIdList(vector<string> & file_idlist, vector<string> &input_idlist)
// {
// 	bool rv = false;
// 	vector<string>::iterator f_itr = file_idlist.begin();
// 	vector<string>::iterator i_itr = input_idlist.begin();
// 	while(f_itr != file_idlist.end() && i_itr != input_idlist.end()) {
// 		if(*i_itr == *f_itr) {
// 			i_itr++;
// 			f_itr++;
// 			continue;
// 		}
// 		else {
// 			CSeq_id seq_id(*i_itr, CSeq_id::fParse_RawText | CSeq_id::fParse_AnyLocal | CSeq_id::fParse_PartialOK);
// 			// Input id is PDB with just mol id
// 			if(seq_id.IsPdb() &&  !seq_id.GetPdb().IsSetChain_id()) {
// 				CSeq_id file_seqid(*f_itr, CSeq_id::fParse_RawText | CSeq_id::fParse_AnyLocal | CSeq_id::fParse_PartialOK);
// 				if (file_seqid.IsPdb() && file_seqid.GetPdb().GetMol().Get() == *i_itr) {
// 					f_itr++;
// 					string tmp_pdb = *i_itr;
// 					while ((f_itr != file_idlist.end()) && ((*f_itr).find_first_of(tmp_pdb) == 0)){
// 						f_itr++;
// 					}
// 					// Skip pdb id in input list but with chain id
// 					while ((i_itr != input_idlist.end()) && ((*i_itr).find_first_of(tmp_pdb) == 0)){
// 						i_itr++;
// 					}
// 					continue;
// 				}
// 			}
// 			else {
// 				CSeq_id file_seq_id(*f_itr, CSeq_id::fParse_RawText | CSeq_id::fParse_AnyLocal | CSeq_id::fParse_PartialOK);
// 				if( file_seq_id.GetSeqIdString(false) == *i_itr) {
// 					i_itr++;
// 					// Skip identical id in input list but with version
// 					if((i_itr != input_idlist.end()) && (file_seq_id.GetSeqIdString(true) == *i_itr)){
// 						i_itr++;
// 					}
// 					f_itr++;
// 					continue;
// 				}
// 			}
// 			break;
// 		}
// 	}
// 	if(f_itr == file_idlist.end()){
// 		rv=true;
// 	}
//
// 	file_idlist.clear();
// 	input_idlist.clear();
// 	return rv;
// }
// ```
// `find_first_of(tmp_pdb) == 0` asks whether the first character of the id is one of the
// characters of the molecule name, not whether the id starts with the name: kept as it is.
// A PDB id of the list with a molecule name that is not the file id's falls through to
// `break` (the `else` is not taken), as in NCBI.
/// Whether the sorted ids of the file (`file_ids`, an OID's `.nos` entries) are all covered
/// by the sorted listed ids `input_ids` that resolved to the OID. `Err`: an id that does not
/// parse (NCBI throws `CSeqIdException`).
fn compare_id_list(file_ids: &[Vec<u8>], input_ids: &[&[u8]]) -> Result<bool, String> {
    let parse = |text: &[u8]| {
        seqidlist::parse_list_id(text).map_err(|error| {
            format!(
                "{}: {}",
                String::from_utf8_lossy(text),
                error.what().trim_end()
            )
        })
    };
    let starts_with_one_of = |s: &[u8], set: &[u8]| s.first().is_some_and(|c| set.contains(c));
    let (mut f, mut i) = (0usize, 0usize);
    while f < file_ids.len() && i < input_ids.len() {
        if input_ids[i] == file_ids[f].as_slice() {
            i += 1;
            f += 1;
            continue;
        }
        let seq_id = parse(input_ids[i])?;
        if let ParsedSeqId::Pdb(pdb) = &seq_id {
            if pdb.chain_id.is_none() {
                let file_seqid = parse(&file_ids[f])?;
                if matches!(&file_seqid, ParsedSeqId::Pdb(file_pdb) if file_pdb.mol == input_ids[i])
                {
                    f += 1;
                    let tmp_pdb = input_ids[i];
                    while f < file_ids.len() && starts_with_one_of(&file_ids[f], tmp_pdb) {
                        f += 1;
                    }
                    while i < input_ids.len() && starts_with_one_of(input_ids[i], tmp_pdb) {
                        i += 1;
                    }
                    continue;
                }
            }
        } else {
            let file_seq_id = parse(&file_ids[f])?;
            if blast_seq_id_string_plain(&file_seq_id, false) == input_ids[i] {
                i += 1;
                if i < input_ids.len()
                    && blast_seq_id_string_plain(&file_seq_id, true) == input_ids[i]
                {
                    i += 1;
                }
                f += 1;
                continue;
            }
        }
        break;
    }
    Ok(f == file_ids.len())
}

/// `CSeq_id::GetSeqIdString(with_version)`.
fn blast_seq_id_string_plain(id: &ParsedSeqId, with_version: bool) -> Vec<u8> {
    let mut out = Vec::new();
    label::content(id, with_version, &mut out);
    out
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

/// A lookup file names an OID the database does not have.
fn bad_lookup_oid(lmdb: &Path, oid: u32, num_oids: u32) -> SeqDbError {
    SeqDbError::Format {
        path: lmdb.to_path_buf(),
        message: format!("the seqid lookup gives OID {oid} of a database of {num_oids} sequences"),
    }
}

/// The OIDs of the sequences `first..end` of a volume read per window of the index.
const SCAN_WINDOW: u32 = 1 << 16;

impl SeqDb {
    /// The OID of ordinal id `ordinal`: `ordinal` itself unless the database is restricted
    /// ([`SeqDb::restrict`]). NCBI prints it as `n` of `gnl|BL_ORD_ID|n`.
    pub fn oid(&self, ordinal: u32) -> u32 {
        match &self.filter {
            Some(filter) => filter.select(ordinal),
            None => ordinal,
        }
    }

    /// Whether the database is restricted by a list ([`SeqDb::restrict`]).
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
        let rule = if negative {
            (!found.is_empty()).then(|| DeflineRule::NegativeTaxids(found))
        } else {
            Some(DeflineRule::Taxids(found))
        };
        self.apply_filter(words, num_oids, rule)?;
        Ok(warning)
    }

    /// Restricts the database to `list` ([`SeqDb::restrict_taxids`], [`SeqDb::restrict_seqids`]);
    /// the warnings NCBI writes go to `warnings`, before an error too.
    pub fn restrict(
        &mut self,
        list: &DbList,
        warnings: &mut Vec<ListWarning>,
    ) -> Result<(), SeqDbError> {
        match list {
            DbList::Taxids(list) => {
                warnings.extend(self.restrict_taxids(list)?.map(ListWarning::Log));
                Ok(())
            }
            DbList::Seqids(list) => self.restrict_seqids(list, warnings),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:236-240
    // ```c++
    //     if (m_UserList.NotEmpty() && m_UserList->NotEmpty()) {
    //     	if(user_list->GetNumSis() > 0) {
    //     		if(s_VerifySeqidlist(user_list->GetListInfo(), volset, lmdb_set)) {
    //     			s_ProcessPositiveSeqIDsForV5(user_list);
    //     		}
    //     	}
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:245-252
    // ```c++
    //     	if (lmdb_set.IsBlastDBVersion5()) {
    //     		if(user_list->GetNumSis() > 0) {
    //     			vector<string> accs;
    //     			vector<blastdb::TOid> oids;
    //     			user_list->GetSiList(accs);
    //     			lmdb_set.AccessionsToOids(accs, oids);
    //     			for(unsigned int i=0; i < accs.size(); i++) {
    //     				user_list->SetSiTranslation(i, oids[i]);
    //     			}
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:294-302
    // ```c++
    // void CSeqDBOIDList::x_ApplyUserGiList(CSeqDBGiList   & gis)
    //
    // {
    //     //m_Atlas.Lock(locked);
    //
    //     if (gis.Empty()) {
    //         x_ClearBitRange(0, m_NumOIDs);
    //         m_NumOIDs = 0;
    //         return;
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:323-329
    // ```c++
    //     if(gis.GetNumSis()) {
    //         for(j = 0; j < gis.GetNumSis(); j++) {
    //             int oid = gis.GetSiOid(j).oid;
    //             if ((oid != -1) && (oid < m_NumOIDs)) {
    //                 gilist_oids->SetBit(oid);
    //             }
    //         }
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:328-331
    // ```c++
    //     		if(m_NegativeList->GetNumSis() > 0) {
    //     			vector<blastdb::TOid> & oids = m_NegativeList->SetExcludedOids();
    //     			lmdb_set.NegativeSeqIdsToOids(m_NegativeList->GetSiList(), oids);
    //
    //     		}
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:365-374
    // ```c++
    // void CSeqDBOIDList::x_ApplyNegativeList(CSeqDBNegativeList & nlist, bool is_v5)
    //
    // {
    //     // We require a normalized list in order to turn bits off.
    // 	m_AllBits->Normalize();
    //     const vector<blastdb::TOid> & excluded_oids = nlist.GetExcludedOids();
    // 	for(unsigned int i=0; i < excluded_oids.size(); i++) {
    // 	    m_AllBits->ClearBit(excluded_oids[i]);
    // 	}
    // ```
    // A text list that lost every id (all `gi|N`): NCBI warns `Empty seqid list` and searches
    // the whole database, positive or negative (probe S17, S18, R08, R09 of
    // dropin-survey-20261029/DI8: 500 hits of 514 sequences, the statistics of the whole
    // database; a positive list with no id would clear every OID by the code above, which
    // the run does not show). LOSAT does what the run shows: the warnings, every OID kept,
    // and, as NCBI's `CSeqDBImpl` still scans the totals for a list (`x_ScanTotals`,
    // approximate nucleotide lengths, the shortest length), the scanned totals.
    // The ids of a binary list are used as written (the alias tool wrote them as the
    // processing writes them). The lists the volumes keep for the filtered header read are
    // the processed ids.
    /// Restricts the database to the OIDs of a seqid list, as NCBI's `CSeqDB` opened with the
    /// list (`CSeqDBGiListSet::x_ResolvePositiveList` / `x_ResolveNegativeList`, then
    /// `CSeqDBOIDList`): the ordinal ids, totals and definition lines as
    /// [`SeqDb::restrict_taxids`] leaves them. The warnings NCBI writes (a text list on a
    /// version 5 database, a list with no id, a binary list of another database) go to
    /// `warnings`, before an error too.
    ///
    /// Errors: NCBI's (`SeqDbError::NotFound`): a binary list on a version 4 database, an
    /// id that does not parse, a database without the accession table; the reader's for
    /// damaged lookup files.
    pub fn restrict_seqids(
        &mut self,
        list: &SeqidList,
        warnings: &mut Vec<ListWarning>,
    ) -> Result<(), SeqDbError> {
        if !self.is_version5() {
            return self.restrict_seqids_v4(list, warnings);
        }
        let ids: Vec<Vec<u8>> = if list.file.info.is_v4 {
            warnings.push(ListWarning::Posted(TEXT_LIST_WARNING.to_string()));
            let ids = seqidlist::process_for_v5(&list.file.ids).map_err(|invalid| {
                warnings.push(ListWarning::Error(invalid.error.message.clone()));
                SeqDbError::NotFound(invalid.message())
            })?;
            if ids.is_empty() {
                warnings.push(ListWarning::Posted(EMPTY_LIST_WARNING.to_string()));
                let num_oids = self.volume_starts[self.volumes.len()];
                let words = vec![u64::MAX; (num_oids as usize).div_ceil(64)];
                return self.apply_filter(words, num_oids, None);
            }
            ids
        } else {
            if list.file.info.db_vol_length > 0
                && list.file.info.db_vol_length != self.whole_total_length
            {
                warnings.push(ListWarning::Posted(DB_INFO_WARNING.to_string()));
            }
            let mut ids = list.file.ids.clone();
            ids.sort();
            ids
        };
        let num_oids = self.volume_starts[self.volumes.len()];
        let words = self.seqid_words_v5(&ids, list.negative, num_oids)?;
        let rule = if list.negative {
            DeflineRule::NegativeSeqids(ids)
        } else {
            DeflineRule::Seqids(ids)
        };
        self.apply_filter(words, num_oids, Some(rule))?;
        if !list.negative {
            self.check_lines_kept()?;
        }
        Ok(())
    }

    /// Reads the definition lines of every included OID of a positive seqid list, so that
    /// an OID the list keeps without a line is refused before the search (`no_line_kept`);
    /// a negative list cannot leave a kept OID without a line short of NCBI's corner cases
    /// of `s_CompareIdList`, which the read of the line refuses in the same words.
    fn check_lines_kept(&self) -> Result<(), SeqDbError> {
        for ordinal in 0..self.num_seqs {
            let oid = self.oid(ordinal);
            let (at, vol_oid) = self.volume_of_oid(oid);
            self.volumes[at].defline_set(vol_oid)?;
        }
        Ok(())
    }

    /// The volume holding OID `oid` and the OID's number in it.
    fn volume_of_oid(&self, oid: u32) -> (usize, u32) {
        let at =
            self.volume_starts[..self.volumes.len()].partition_point(|&start| start <= oid) - 1;
        (at, oid - self.volume_starts[at])
    }

    /// The bitmap of the OIDs a seqid list keeps in a version 5 database (`ids`: sorted and
    /// unique, as the lookup wants them).
    fn seqid_words_v5(
        &self,
        ids: &[Vec<u8>],
        negative: bool,
        num_oids: u32,
    ) -> Result<Vec<u64>, SeqDbError> {
        let entries = self.lmdb_entries()?;
        let word_count = (num_oids as usize).div_ceil(64);
        if negative {
            let mut words = vec![u64::MAX; word_count];
            for entry in &entries {
                for oid in entry.negative_seqids_to_oids(ids)? {
                    if oid >= num_oids {
                        return Err(bad_lookup_oid(&entry.lmdb, oid, num_oids));
                    }
                    words[(oid / 64) as usize] &= !(1u64 << (oid % 64));
                }
            }
            return Ok(words);
        }
        let mut words = vec![0u64; word_count];
        let mut oids: Vec<Option<u32>> = vec![None; ids.len()];
        for (at, entry) in entries.iter().enumerate() {
            let found = entry.accessions_to_oids(ids)?;
            for (slot, oid) in oids.iter_mut().zip(found) {
                if at == 0 || oid.is_some() {
                    *slot = oid;
                }
            }
        }
        for oid in oids.into_iter().flatten() {
            if oid < num_oids {
                words[(oid / 64) as usize] |= 1u64 << (oid % 64);
            }
        }
        Ok(words)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:266-301
    // ```c++
    //     	if((user_list->GetNumSis() > 0) && !(lmdb_set.IsBlastDBVersion5())) {
    //
    //     		user_list->PreprocessIdsForISAMSiLookup();
    //     	}
    //
    //         typedef SSeqDB_IndexCountPair TIndexCount;
    //         vector<TIndexCount> OidsPerVolume;
    //
    //         // Build a list of volumes sorted by OID count.
    //
    //         for(int i = 0; i < volset.GetNumVols(); i++) {
    //             const CSeqDBVolEntry * vol = volset.GetVolEntry(i);
    //
    //             TIndexCount vol_oids;
    //             vol_oids.m_Index = i;
    //             vol_oids.m_Count = vol->OIDEnd() - vol->OIDStart();
    //
    //             OidsPerVolume.push_back(vol_oids);
    //         }
    //
    //         // The largest volumes should be used first, to minimize the
    //         // number of failed GI->OID conversion attempts.  Searching input
    //         // GIs against larger volumes first should eliminate most of the
    //         // GIs by the time smaller volumes are searched, thus reducing the
    //         // total number of lookups.
    //
    //         std::sort(OidsPerVolume.begin(), OidsPerVolume.end());
    //
    //         for(int i = 0; i < (int)OidsPerVolume.size(); i++) {
    //             int vol_idx = OidsPerVolume[i].m_Index;
    //
    //             const CSeqDBVolEntry * vol = volset.GetVolEntry(vol_idx);
    //
    //             // Note: The implied ISAM lookups will sort by GI/TI.
    //
    //             vol->Vol()->IdsToOids(*m_UserList, locked);
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:345-347
    // ```c++
    //     	if((m_NegativeList->GetNumSis() > 0) && !(lmdb_set.IsBlastDBVersion5())) {
    //     		m_NegativeList->PreprocessIdsForISAMSiLookup();
    //     	}
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:369-375
    // ```c++
    //         for(int i = 0; i < volset.GetNumVols(); i++) {
    //             const CSeqDBVolEntry * vol = volset.GetVolEntry(i);
    //
    //             // Note: The implied ISAM lookups will sort by GI/TI.
    //
    //             vol->Vol()->IdsToOids(*m_NegativeList, locked);
    //         }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2540-2550
    // ```c++
    //     if (ids.GetNumSis() && (GetLMDBFileName() == kEmptyStr)) {
    //         x_OpenStrFile();
    //         if (m_IsamStr.NotEmpty()) {
    //             m_IsamStr->IdsToOids(m_VolStart, m_VolEnd, ids);
    //             x_UnleaseStrFile();
    //         } else {
    //             NCBI_THROW(CSeqDBException,
    //                        eArgErr,
    //                        "SI list specified but no ISAM file found for SI in " + m_VolName);
    //         }
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:375-402
    // ```c++
    // 	if((!is_v5 && nlist.GetNumSis() > 0) || nlist.GetNumGis() > 0 || nlist.GetNumTis() >  0) {
    //
    //     // Intersect the user GI list with the OID bit map.
    //
    //     // Iterate over the bitmap, clearing bits we find there but not in
    //     // the bool vector.  For very dense OID bit maps, it might be
    //     // faster to use two similarly implemented bitmaps and AND them
    //     // together word-by-word.
    //
    //     int max = nlist.GetNumOids();
    //
    //     // Clear any OIDs after the included range.
    //
    //     if (max < m_NumOIDs) {
    //         CSeqDB_BitSet new_range(0, max, CSeqDB_BitSet::eAllSet);
    //         m_AllBits->IntersectWith(new_range, true);
    //     }
    //
    //     // If a 'get next included oid' method was added to the negative
    //     // list, the following loop could be made a bit faster.
    //
    //     for(int oid = 0; oid < max; oid++) {
    //         if (! nlist.GetOidStatus(oid)) {
    //             m_AllBits->ClearBit(oid);
    //         }
    //     }
    // 	}
    // ```
    // A text list is not rewritten for a version 4 database (`s_VerifySeqidlist` returns
    // false): every id, simplified and lower-cased, is a key of the string index; an id that
    // devolves to a number (a gi) becomes the empty key, which no line has (probe R12, R13).
    // An exception of `SeqDB_SimplifyAccession` (an id such as `sp||`) is not a database
    // error in NCBI (`Error: <what>`, exit status 255); LOSAT refuses it with the text.
    /// [`SeqDb::restrict_seqids`] on a version 4 database: the string ISAM lookups of each
    /// volume (`CSeqDBVol::IdsToOids`), after `PreprocessIdsForISAMSiLookup`.
    fn restrict_seqids_v4(
        &mut self,
        list: &SeqidList,
        _warnings: &mut [ListWarning],
    ) -> Result<(), SeqDbError> {
        if !list.file.info.is_v4 {
            return Err(SeqDbError::NotFound(BINARY_LIST_V4_ERROR.to_string()));
        }
        let protein = matches!(self.seqtype, SeqType::Protein);
        let mut ids: Vec<Vec<u8>> = Vec::with_capacity(list.file.ids.len());
        for id in &list.file.ids {
            ids.push(isam::simplify_accession_for_isam(id).map_err(|error| {
                SeqDbError::Unsupported {
                    path: self.volumes[0].base().to_path_buf(),
                    message: format!(
                        "the id {} of the seqid list throws in NCBI's SeqDB_SimplifyAccession: {}",
                        String::from_utf8_lossy(id),
                        error.to_string().trim_end()
                    ),
                }
            })?);
        }
        let string_isam = |at: usize| -> Result<isam::StringIsam, SeqDbError> {
            let base = self.volumes[at].base();
            let name = base.to_string_lossy();
            if !isam::StringIsam::index_exists(&name, protein) {
                return Err(SeqDbError::NotFound(format!(
                    "SI list specified but no ISAM file found for SI in {name}"
                )));
            }
            isam::StringIsam::open_volume(&name, protein)
        };
        let num_oids = self.volume_starts[self.volumes.len()];
        let word_count = (num_oids as usize).div_ceil(64);
        let words = if list.negative {
            let mut out = isam::NegativeOids::new();
            for at in 0..self.volumes.len() {
                string_isam(at)?.negative_scan(
                    self.volume_starts[at],
                    self.volume_starts[at + 1],
                    &mut ids,
                    &mut out,
                )?;
            }
            let mut words = vec![u64::MAX; word_count];
            let max = out.num_oids().min(num_oids);
            for oid in 0..num_oids {
                if oid >= max || !out.get_oid_status(oid) {
                    words[(oid / 64) as usize] &= !(1u64 << (oid % 64));
                }
            }
            words
        } else {
            let mut entries: Vec<isam::SiOid> =
                ids.iter().map(|id| isam::SiOid::new(id.clone())).collect();
            let mut order: Vec<usize> = (0..self.volumes.len()).collect();
            order.sort_by_key(|&at| std::cmp::Reverse(self.volumes[at].num_oids()));
            for at in order {
                string_isam(at)?.ids_to_oids(self.volume_starts[at], &mut entries)?;
            }
            let mut words = vec![0u64; word_count];
            for oid in entries.iter().filter_map(|entry| entry.oid) {
                if oid < num_oids {
                    words[(oid / 64) as usize] |= 1u64 << (oid % 64);
                }
            }
            words
        };
        ids.sort();
        let rule = if list.negative {
            DeflineRule::NegativeSeqids(ids)
        } else {
            DeflineRule::Seqids(ids)
        };
        self.apply_filter(words, num_oids, Some(rule))?;
        if !list.negative {
            self.check_lines_kept()?;
        }
        Ok(())
    }

    /// Makes `words` (bit `oid` set for an included OID) the database's OID set: the
    /// ordinal ids, NCBI's scanned totals and the definition lines `rule` keeps.
    fn apply_filter(
        &mut self,
        words: Vec<u64>,
        num_oids: u32,
        rule: Option<DeflineRule>,
    ) -> Result<(), SeqDbError> {
        let filter = OidFilter::new(words, num_oids);
        let (total_length, max_length, min_length) = self.scan_totals(&filter)?;
        if let Some(rule) = rule.map(Arc::new) {
            for volume in &mut self.volumes {
                volume.set_defline_rule(Arc::clone(&rule));
            }
        }
        // OID 0 is the first OID of the first volume that has any.
        let first_volume = self.volumes.iter().find(|volume| volume.num_oids() > 0);
        self.empty_list_first_length = match first_volume {
            Some(volume) if filter.count() == 0 => Some(match self.seqtype {
                SeqType::Protein => u64::from(volume.protein_length(0)?),
                SeqType::Nucleotide => volume.nucleotide_length(0, 0)? as u64,
            }),
            _ => None,
        };
        self.num_seqs = filter.count();
        self.total_length = total_length;
        self.max_length = max_length;
        self.min_length = min_length;
        self.filter = Some(Box::new(filter));
        Ok(())
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
    /// NCBI's totals of the included OIDs: the total length, the longest length and the
    /// shortest (`INT4_MAX` when no OID is included, as `min_count` starts).
    fn scan_totals(&self, filter: &OidFilter) -> Result<(u64, u32, u32), SeqDbError> {
        let (mut total, mut longest, mut shortest) = (0u64, 0u32, i32::MAX as u32);
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
                    shortest = shortest.min(length);
                }
                first = end;
            }
        }
        Ok((total, longest, shortest))
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

    fn ids_db(name: &str, seqtype: SeqType) -> SeqDb {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/db_search/ids")
            .join(name);
        SeqDb::open(path.to_str().unwrap(), seqtype).unwrap()
    }

    fn list_of(text: &[u8], negative: bool) -> SeqidList {
        let mut warnings = Vec::new();
        SeqidList {
            negative,
            file: seqidlist::read_si_list(text, &mut warnings).unwrap(),
        }
    }

    /// The Seq-id strings (`GetSeqIdString(true)`) of the definition lines the restricted
    /// database keeps, per ordinal id.
    fn kept_ids(db: &SeqDb) -> Vec<Vec<String>> {
        (0..db.num_seqs())
            .map(|ordinal| {
                let (at, vol_oid) = db.volume_of_oid(db.oid(ordinal));
                let set = db.volumes[at].defline_set(vol_oid).unwrap();
                set.deflines
                    .iter()
                    .flat_map(|d| d.ids.iter())
                    .map(|id| {
                        let mut out = Vec::new();
                        label::content(id, true, &mut out);
                        String::from_utf8(out).unwrap()
                    })
                    .collect()
            })
            .collect()
    }

    #[test]
    fn a_positive_seqid_list_keeps_the_named_oids_and_their_lines() {
        let whole = ids_db("prot_v5", SeqType::Protein);
        let mut db = ids_db("prot_v5", SeqType::Protein);
        let mut warnings = Vec::new();
        db.restrict(
            &DbList::Seqids(list_of(b">DAA10001.1\n", false)),
            &mut warnings,
        )
        .unwrap();
        assert_eq!(
            warnings,
            [ListWarning::Posted(TEXT_LIST_WARNING.to_string())]
        );
        assert!(db.is_restricted());
        assert_eq!(db.num_seqs(), 1);
        assert_eq!(kept_ids(&db), [vec!["DAA10001.1".to_string()]]);
        assert!(db.total_length() > 0 && db.total_length() < whole.total_length());
        assert_eq!(db.printed_num_seqs(), whole.num_seqs());

        // Versionless accessions, a pdb chain, a general id as db:tag, an unknown id, a gi.
        let mut db = ids_db("prot_v5", SeqType::Protein);
        db.restrict(
            &DbList::Seqids(list_of(
                b"NP_100006 XP_000011\n2XYZ_B mydb:prot10\nUNKNOWN1\ngi|1414\n",
                false,
            )),
            &mut Vec::new(),
        )
        .unwrap();
        assert_eq!(db.num_seqs(), 4);
        for ids in kept_ids(&db) {
            assert!(
                ids.iter().any(|id| {
                    ["NP_100006.1", "XP_000011.2", "2XYZ_B", "mydb:prot10"].contains(&id.as_str())
                }),
                "{ids:?}"
            );
        }
    }

    #[test]
    fn a_negative_seqid_list_excludes_an_oid_only_when_every_id_is_covered() {
        let whole = ids_db("prot_v5", SeqType::Protein).num_seqs();
        // `lcl|Subject_28|tpg|DAA10008.1|`: one id listed keeps the OID, both exclude it (the
        // versionless `DAA10008` covers the stored `DAA10008.1`).
        let mut db = ids_db("prot_v5", SeqType::Protein);
        db.restrict(
            &DbList::Seqids(list_of(b"Subject_28\n", true)),
            &mut Vec::new(),
        )
        .unwrap();
        assert_eq!(db.num_seqs(), whole);
        assert!(db.is_restricted());
        let mut db = ids_db("prot_v5", SeqType::Protein);
        db.restrict(
            &DbList::Seqids(list_of(b"Subject_28\nDAA10008\n", true)),
            &mut Vec::new(),
        )
        .unwrap();
        assert_eq!(db.num_seqs(), whole - 1);
        // A definition line whose only ids are gis or `ti` general ids is dropped by a
        // negative seqid list (`CSeqDBNegativeList::FindId` matches no type): an OID of such
        // lines alone keeps none, which the read refuses (DECISIONS 17); every other kept OID
        // reads.
        let unrestricted = ids_db("prot_v5", SeqType::Protein);
        let gi_only: Vec<u32> = (0..unrestricted.num_seqs())
            .filter(|&oid| {
                let (at, vol_oid) = unrestricted.volume_of_oid(oid);
                let set = unrestricted.volumes[at].defline_set(vol_oid).unwrap();
                set.deflines.iter().all(|d| {
                    d.ids.iter().all(|id| {
                        matches!(id, ParsedSeqId::Gi(_))
                            || matches!(id, ParsedSeqId::General(tag) if tag.db == b"ti")
                    })
                })
            })
            .collect();
        assert!(!gi_only.is_empty());
        for ordinal in 0..db.num_seqs() {
            let oid = db.oid(ordinal);
            let (at, vol_oid) = db.volume_of_oid(oid);
            match db.volumes[at].defline_set(vol_oid) {
                Ok(set) => {
                    assert!(!gi_only.contains(&oid), "{oid}");
                    assert!(set.deflines.iter().all(|d| {
                        d.ids.iter().all(|id| {
                            let mut out = Vec::new();
                            label::content(id, true, &mut out);
                            out != b"DAA10008.1"
                        })
                    }));
                }
                Err(error) => {
                    assert!(gi_only.contains(&oid), "{oid}: {error}");
                    assert!(error.to_string().contains("none of its definition lines"));
                }
            }
        }
        // A nucleotide local id of its own OID.
        let whole = ids_db("nucl_v5", SeqType::Nucleotide).num_seqs();
        let mut db = ids_db("nucl_v5", SeqType::Nucleotide);
        db.restrict(
            &DbList::Seqids(list_of(b"mynuc6\nUNKNOWN1\n", true)),
            &mut Vec::new(),
        )
        .unwrap();
        assert_eq!(db.num_seqs(), whole - 1);
    }

    #[test]
    fn a_list_that_loses_every_id_leaves_a_version_5_database_whole() {
        let whole = ids_db("prot_v5", SeqType::Protein).num_seqs();
        for negative in [false, true] {
            let mut db = ids_db("prot_v5", SeqType::Protein);
            let mut warnings = Vec::new();
            db.restrict(
                &DbList::Seqids(list_of(b"gi|1414\ngi|2\n", negative)),
                &mut warnings,
            )
            .unwrap();
            assert_eq!(
                warnings,
                [
                    ListWarning::Posted(TEXT_LIST_WARNING.to_string()),
                    ListWarning::Posted(EMPTY_LIST_WARNING.to_string())
                ]
            );
            assert!(db.is_restricted());
            assert_eq!(db.num_seqs(), whole);
            assert_eq!(kept_ids(&db).len(), whole as usize);
            assert!(db.min_length().is_some_and(|min| min > 0));
        }
    }

    // Probes T08-T12, C01-C97 of line DI-8: 0 OIDs, NCBI's statistics read OID 0 and one sequence.
    #[test]
    fn a_list_without_a_kept_sequence_has_the_statistics_of_oid_zero() {
        for (name, seqtype) in [
            ("prot_v5", SeqType::Protein),
            ("nucl_v5", SeqType::Nucleotide),
        ] {
            let mut db = ids_db(name, seqtype);
            let first = match seqtype {
                SeqType::Protein => u64::from(db.protein_length(0).unwrap()),
                SeqType::Nucleotide => db.nucleotide_length(0).unwrap() as u64,
            };
            let whole = (db.total_length(), db.num_seqs());
            assert!(!db.filtered_empty());
            assert_eq!(db.statistics_numbers(), whole, "{name}");
            let mut warnings = Vec::new();
            db.restrict(
                &DbList::Seqids(list_of(b"UNKNOWN1\n", false)),
                &mut warnings,
            )
            .unwrap();
            assert!(db.filtered_empty(), "{name}");
            assert_eq!((db.num_seqs(), db.total_length()), (0, 0), "{name}");
            assert_eq!(db.statistics_numbers(), (first, 1), "{name}");
            assert_eq!(db.printed_num_seqs(), whole.1, "{name}");
            assert!(db.next_block(0, 1 << 20).is_none(), "{name}");
        }
    }

    #[test]
    fn a_binary_list_selects_what_its_text_list_selects() {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/seqidlist");
        let text = std::fs::read(dir.join("sl_nucl.txt")).unwrap();
        let binary = std::fs::read(dir.join("sl_nucl.bsl")).unwrap();
        let mut from_text = ids_db("nucl_v5", SeqType::Nucleotide);
        let mut warnings = Vec::new();
        from_text
            .restrict(&DbList::Seqids(list_of(&text, false)), &mut warnings)
            .unwrap();
        assert_eq!(warnings.len(), 1);
        let mut from_binary = ids_db("nucl_v5", SeqType::Nucleotide);
        let mut warnings = Vec::new();
        from_binary
            .restrict(&DbList::Seqids(list_of(&binary, false)), &mut warnings)
            .unwrap();
        assert!(warnings.is_empty());
        assert!(from_text.num_seqs() >= 6);
        assert_eq!(from_text.num_seqs(), from_binary.num_seqs());
        assert_eq!(kept_ids(&from_text), kept_ids(&from_binary));
        assert_eq!(from_text.total_length(), from_binary.total_length());
        // A recorded database length that is not this database's: NCBI's warning, the list used.
        let badlen = std::fs::read(dir.join("sl_tfp_badlen.bsl")).unwrap();
        let mut db = ids_db("prot_v5", SeqType::Protein);
        let mut warnings = Vec::new();
        let result = db.restrict(&DbList::Seqids(list_of(&badlen, false)), &mut warnings);
        assert_eq!(warnings, [ListWarning::Posted(DB_INFO_WARNING.to_string())]);
        // tfp's ids are not in prot_v5: no OID is kept (the search programs refuse that).
        assert!(result.is_ok() && db.num_seqs() == 0);
    }

    #[test]
    fn a_sequence_kept_without_a_definition_line_is_refused() {
        // `lcl|Subject_29|pdb|3ABC|C`: the molecule name is a lookup key, no line's id.
        let mut db = ids_db("prot_v5", SeqType::Protein);
        let error = db
            .restrict(
                &DbList::Seqids(list_of(b"3ABC\nDAA10001.1\n", false)),
                &mut Vec::new(),
            )
            .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("but none of its definition lines, and NCBI leaves its hits out"),
            "{error}"
        );
        assert!(error.ncbi_message().is_none());
    }

    #[test]
    fn a_version_4_database_refuses_a_binary_seqid_list() {
        let binary = std::fs::read(
            Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/seqidlist/sl_tfp.bsl"),
        )
        .unwrap();
        let mut db = ids_db("prot_v4", SeqType::Protein);
        let error = db
            .restrict(&DbList::Seqids(list_of(&binary, false)), &mut Vec::new())
            .unwrap_err();
        assert_eq!(error.ncbi_message().as_deref(), Some(BINARY_LIST_V4_ERROR));
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
        assert_eq!(rule.keep_mask(&set), Some(vec![false, true]));
        let none = DeflineSet {
            deflines: vec![defline(Some(1), None), defline(Some(2), None)],
        };
        assert_eq!(rule.keep_mask(&none), Some(vec![true, true]));
        let mut lines = vec![1, 2, 3];
        lines.retain_lines(&[true, false, true]);
        assert_eq!(lines, vec![1, 3]);
        lines.retain_lines(&[false]);
        assert_eq!(lines, vec![1, 3]);
    }
}
