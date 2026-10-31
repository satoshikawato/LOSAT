//! A seqid list file (line DI-8b: `-seqidlist`, `-negative_seqidlist`), read as NCBI's
//! `CSeqDBFileGiList(fname, eSiList)` reads it (`SeqDB_ReadSiList`): a binary list made by
//! `blastdb_aliastool -seqid_file_in` (its first byte is NUL; `CSeqidlistRead`), else text
//! (`SeqDB_ReadMemorySiList`: white-space separated tokens, a leading `>` skipped, `#` to
//! the end of the line a comment). The ids of a text list are then parsed and rewritten
//! for a version 5 database as `s_ProcessSeqIDsForV5` rewrites them: `gi|N` dropped, `pir`
//! and `prf` ids in FASTA form, every other id as `CSeq_id::GetSeqIdString(true)` (the
//! upper-cased accession with its version when given), sorted and unique.

use crate::seqid::accguide::parse_flags::{ANY_LOCAL, PARTIAL_OK, RAW_TEXT};
use crate::seqid::label;
use crate::seqid::parse::{seq_id_from_text, SeqIdError};
use crate::seqid::ParsedSeqId;

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdbcommon.hpp:134-145
// ```c++
// struct NCBI_XOBJREAD_EXPORT SBlastSeqIdListInfo {
// 		SBlastSeqIdListInfo() : is_v4(true), file_size(0), num_ids(0), create_date(kEmptyStr),
// 				      db_vol_length(0), db_create_date(kEmptyStr), db_vol_names(kEmptyStr) {}
// 		bool is_v4;
// 		Uint8 file_size;
// 		Uint8 num_ids;
// 		string title;
// 		string create_date;
// 		Uint8 db_vol_length;
// 		string db_create_date;
// 		string db_vol_names;
// };
// ```
/// What the list file says of itself (`SBlastSeqIdListInfo`): a text list is `is_v4`, a binary
/// list carries its header.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeqidListInfo {
    /// A text list (`is_v4 = true`, the constructor's default); `false` for a binary list.
    pub is_v4: bool,
    pub file_size: u64,
    pub num_ids: u64,
    pub title: Vec<u8>,
    pub create_date: Vec<u8>,
    /// The total length of the database the binary list was made for (0: none recorded).
    pub db_vol_length: u64,
    pub db_create_date: Vec<u8>,
    pub db_vol_names: Vec<u8>,
}

impl Default for SeqidListInfo {
    fn default() -> Self {
        SeqidListInfo {
            is_v4: true,
            file_size: 0,
            num_ids: 0,
            title: Vec::new(),
            create_date: Vec::new(),
            db_vol_length: 0,
            db_create_date: Vec::new(),
            db_vol_names: Vec::new(),
        }
    }
}

/// A seqid list file as read: its ids as the file gives them and its header.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeqidListFile {
    /// The ids: the tokens of a text list (as written), the ids of a binary list.
    pub ids: Vec<Vec<u8>>,
    pub info: SeqidListInfo,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqidlist_reader.cpp:58-85
// ```c++
// CSeqidlistRead::CSeqidlistRead (CMemoryFile & file) : m_Ptr((char*) file.GetPtr()), m_EndPtr((char*) file.GetPtr()) {
// 	if(m_Ptr == NULL) {
// 		NCBI_THROW(CSeqDBException, eArgErr, "Failed to map seqidlist file ");
// 	}
//
// 	char null_byte = x_GetChar();
// 	if (null_byte == 0) {
// 		m_info.is_v4 = false;
// 		Uint8 file_size = file.GetFileSize();
// 		m_info.file_size = x_GetUint8();
// 		if (m_info.file_size != file_size) {
// 			NCBI_THROW(CSeqDBException, eArgErr, "Invalid seqidlist file");
// 		}
// 		m_EndPtr += file_size;
// 		m_info.num_ids = x_GetUint8();
// 		Uint4 title_length = x_GetUint4();
// 		x_GetString(m_info.title, title_length);
// 		char file_create_date_length = x_GetChar();
// 		x_GetString(m_info.create_date, file_create_date_length);
// 		m_info.db_vol_length = x_GetUint8();
// 		if(m_info.db_vol_length != 0) {
// 			char file_db_create_date_length = x_GetChar();
// 			x_GetString(m_info.db_create_date, file_db_create_date_length);
// 			Uint4 file_vol_names_length = x_GetUint4();
// 			x_GetString(m_info.db_vol_names, file_vol_names_length);
// 		}
// 	}
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqidlist_reader.cpp:87-108
// ```c++
// int CSeqidlistRead::GetIds(vector<CSeqDBGiList::SSiOid>  & idlist)
// {
// 	const unsigned char byte_max = 0xFF;
// 	unsigned int i = 0;
// 	idlist.clear();
// 	idlist.resize(m_info.num_ids);
// 	for(; (m_Ptr < m_EndPtr) && (i < m_info.num_ids); i++) {
// 		unsigned char id_len = (unsigned char) x_GetChar();
// 		if(id_len == byte_max) {
// 			Uint4 long_id_len = x_GetUint4();
// 			x_GetString(idlist[i].si, long_id_len);
// 		}
// 		else {
// 			x_GetString(idlist[i].si, id_len);
// 		}
// 	}
// 	if(i != m_info.num_ids) {
// 		NCBI_THROW(CSeqDBException, eArgErr, "Invalid total num of ids in seqidlist file");
// 	}
//
// 	return i;
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqidlist_reader.cpp:111-120
// ```c++
// int CBlastSeqidlistFile::GetSeqidlist(CMemoryFile & file, vector<CSeqDBGiList::SSiOid>  & idlist,
// 		                              SBlastSeqIdListInfo & list_info)
// {
//
// 	CSeqidlistRead list(file);
// 	list.GetListInfo(list_info);
// 	list.GetIds(idlist);
//
// 	return list_info.num_ids;
// }
// ```
// A 0-byte file has no memory map (`CMemoryFile` leaves the pointer null), hence the first
// error. NCBI reads the header and the ids from the map without bounds; a header or an id
// cut off by the end of the file reads memory past it. LOSAT gives the second error for such
// a file instead.
/// The header of a binary list and its ids: `Ok(None)` for a file that is not binary (its
/// first byte is not NUL), `Err` with NCBI's message for a damaged one.
fn read_binary(bytes: &[u8]) -> Result<Option<(SeqidListInfo, Vec<Vec<u8>>)>, String> {
    const INVALID: &str = "Invalid seqidlist file";
    if bytes.first() != Some(&0) {
        return Ok(None);
    }
    let mut cursor = Cursor { bytes, at: 1 };
    let mut info = SeqidListInfo {
        is_v4: false,
        ..SeqidListInfo::default()
    };
    info.file_size = cursor.u64().ok_or(INVALID)?;
    if info.file_size != bytes.len() as u64 {
        return Err(INVALID.to_string());
    }
    info.num_ids = cursor.u64().ok_or(INVALID)?;
    let title_length = cursor.u32().ok_or(INVALID)? as usize;
    info.title = cursor.take(title_length).ok_or(INVALID)?.to_vec();
    let create_date_length = usize::from(cursor.u8().ok_or(INVALID)?);
    info.create_date = cursor.take(create_date_length).ok_or(INVALID)?.to_vec();
    info.db_vol_length = cursor.u64().ok_or(INVALID)?;
    if info.db_vol_length != 0 {
        let db_create_date_length = usize::from(cursor.u8().ok_or(INVALID)?);
        info.db_create_date = cursor.take(db_create_date_length).ok_or(INVALID)?.to_vec();
        let vol_names_length = cursor.u32().ok_or(INVALID)? as usize;
        info.db_vol_names = cursor.take(vol_names_length).ok_or(INVALID)?.to_vec();
    }
    let mut ids = Vec::with_capacity(usize::try_from(info.num_ids).unwrap_or(0).min(1 << 20));
    while cursor.at < bytes.len() && (ids.len() as u64) < info.num_ids {
        let id_len = cursor.u8().ok_or(INVALID)?;
        let len = if id_len == 0xFF {
            cursor.u32().ok_or(INVALID)? as usize
        } else {
            usize::from(id_len)
        };
        ids.push(cursor.take(len).ok_or(INVALID)?.to_vec());
    }
    if ids.len() as u64 != info.num_ids {
        return Err("Invalid total num of ids in seqidlist file".to_string());
    }
    Ok(Some((info, ids)))
}

/// A read position in the binary list (`m_Ptr`); every read is bounded by the file's end.
struct Cursor<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, len: usize) -> Option<&'a [u8]> {
        let end = self
            .at
            .checked_add(len)
            .filter(|&end| end <= self.bytes.len())?;
        let piece = &self.bytes[self.at..end];
        self.at = end;
        Some(piece)
    }

    fn u8(&mut self) -> Option<u8> {
        self.take(1).map(|b| b[0])
    }

    fn u32(&mut self) -> Option<u32> {
        self.take(4)
            .map(|b| u32::from_le_bytes(b.try_into().expect("4 bytes")))
    }

    fn u64(&mut self) -> Option<u64> {
        self.take(8)
            .map(|b| u64::from_le_bytes(b.try_into().expect("8 bytes")))
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1284-1322
// ```c++
// void SeqDB_ReadMemorySiList(const char * fbeginp,
//                             const char * fendp,
//                             vector<CSeqDBGiList::SSiOid> & sis,
//                             bool * in_order)
// {
//     Int8 file_size = fendp - fbeginp;
//
//     // We would prefer to do only one allocation, so assume
//     // average seqid is 6 digits plus newline.  A few extra will be
//     // allocated, but this is preferable to letting the vector
//     // double itself (which it still will do if needed).
//
//     sis.reserve(sis.size() + int(file_size / 7));
//
//     const char * p = fbeginp;
//     const char * head;
//     while ( p < fendp) {
//         // find the head of the seqid
//         while (p< fendp && (*p=='>' || *p==' ' || *p=='\t' || *p=='\n' || *p=='\r')) ++p;
//         if (p< fendp && *p == '#') {
//             // anything beyond this point in the line is a comment
//             while (p< fendp && *p!='\n') ++p;
//             continue;
//         }
//         head = p;
//         while (p< fendp && *p!=' ' && *p!='\t' && *p!='\n' && *p!='\r') ++p;
//         if (p > head) {
//             string acc(head, p);
//             string str_id = NStr::TruncateSpaces(acc, NStr::eTrunc_Both);
//             if (str_id != "") {
//                 sis.push_back(str_id);
//             } else {
//                 cerr << "WARNING:  " << acc
//                      << " is not a valid seqid string." << endl;
//             }
//         }
//     }
//     if (in_order) *in_order = false;
// }
// ```
// A token holds no space, tab, line feed or carriage return, so `TruncateSpaces` (which
// trims `isspace`: also form feed and vertical tab) can empty it only when it is made of
// those two characters; such a token is reported and dropped.
/// The tokens of a text list.
fn read_memory_si_list(bytes: &[u8], warnings: &mut Vec<String>) -> Vec<Vec<u8>> {
    let mut sis = Vec::with_capacity(bytes.len() / 7);
    let mut p = 0;
    while p < bytes.len() {
        while p < bytes.len() && matches!(bytes[p], b'>' | b' ' | b'\t' | b'\n' | b'\r') {
            p += 1;
        }
        if p < bytes.len() && bytes[p] == b'#' {
            while p < bytes.len() && bytes[p] != b'\n' {
                p += 1;
            }
            continue;
        }
        let head = p;
        while p < bytes.len() && !matches!(bytes[p], b' ' | b'\t' | b'\n' | b'\r') {
            p += 1;
        }
        if p > head {
            let acc = &bytes[head..p];
            let str_id = crate::seqid::parse::truncate_spaces(acc);
            if !str_id.is_empty() {
                sis.push(str_id.to_vec());
            } else {
                warnings.push(format!(
                    "WARNING:  {} is not a valid seqid string.",
                    String::from_utf8_lossy(acc)
                ));
            }
        }
    }
    sis
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1476-1491
// ```c++
// void SeqDB_ReadSiList(const string & fname, vector<CSeqDBGiList::SSiOid> & sis, bool * in_order, SBlastSeqIdListInfo & db_info)
// {
// 	CMemoryFile mfile(SeqDB_MakeOSPath(fname));
//     if (CBlastSeqidlistFile::GetSeqidlist(mfile, sis, db_info)) {
//     	*in_order = true;
//     	return;
//     }
//     else {
//     	Int8 file_size = mfile.GetSize();
//     	const char *fbeginp = (char*) mfile.GetPtr();
//     	const char *fendp   = fbeginp + file_size;
//     	SeqDB_ReadMemorySiList(fbeginp, fendp, sis, in_order);
//     }
// }
// ```
// A binary list of no ids (`GetSeqidlist` returns 0) is read again as text: its header
// bytes become tokens. The header of a text file is the default (`is_v4`).
/// Reads a seqid list file's bytes: a binary list, else a text list. `Err`: NCBI's
/// `CSeqDBException` text (`BLAST Database error: ...`, exit status 2). The lines NCBI
/// writes to `cerr` for tokens it drops are pushed to `warnings`.
pub fn read_si_list(bytes: &[u8], warnings: &mut Vec<String>) -> Result<SeqidListFile, String> {
    if bytes.is_empty() {
        return Err("Failed to map seqidlist file ".to_string());
    }
    if let Some((info, ids)) = read_binary(bytes)? {
        if info.num_ids > 0 {
            return Ok(SeqidListFile { ids, info });
        }
    }
    Ok(SeqidListFile {
        ids: read_memory_si_list(bytes, warnings),
        info: SeqidListInfo::default(),
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:148-181
// ```c++
// void s_ProcessSeqIDsForV5 (vector<string> & idlist)
// {
// 	 vector<string> tmplist;
// 	 tmplist.reserve(idlist.size());
//     for(unsigned int i=0; i < idlist.size(); i++) {
//             try {
//             CSeq_id seqid(idlist[i], CSeq_id::fParse_RawText | CSeq_id::fParse_AnyLocal | CSeq_id::fParse_PartialOK);
//             if(seqid.IsGi()) {
//                     continue;
//             }
//
//             if(seqid.IsPir() || seqid.IsPrf()) {
//                     string id = seqid.AsFastaString();
//                     tmplist.push_back(id);
//                     continue;
//             }
//             tmplist.push_back(seqid.GetSeqIdString(true));
//             } catch (CException & e) {
//                     ERR_POST(e.GetMsg());
//                     NCBI_THROW(CSeqDBException, eArgErr, "Invalid seq id: " + idlist[i]);
//
//             }
//     }
//
//     if (tmplist.size() == 0) {
//         ERR_POST(Warning << "Empty seqid list");
//     }
//     else {
//     	sort(tmplist.begin(), tmplist.end());
//     	vector<string>::iterator it = unique (tmplist.begin(), tmplist.end());
//     	tmplist.resize(distance(tmplist.begin(),it));
//     }
//    	idlist.swap(tmplist);
// }
// ```
/// The parse of a token failed: NCBI posts the exception's message, then throws
/// `Invalid seq id: <token>`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InvalidSeqId {
    pub token: Vec<u8>,
    pub error: SeqIdError,
}

impl InvalidSeqId {
    /// The `CSeqDBException` text (`BLAST Database error: Invalid seq id: <token>`).
    pub fn message(&self) -> String {
        format!("Invalid seq id: {}", String::from_utf8_lossy(&self.token))
    }
}

/// The ids of a text list as a version 5 database looks them up (`s_ProcessSeqIDsForV5`):
/// sorted and unique. An empty result is NCBI's `Empty seqid list` warning (the caller
/// writes it).
pub fn process_for_v5(ids: &[Vec<u8>]) -> Result<Vec<Vec<u8>>, InvalidSeqId> {
    let mut out: Vec<Vec<u8>> = Vec::with_capacity(ids.len());
    for token in ids {
        let id = parse_list_id(token).map_err(|error| InvalidSeqId {
            token: token.clone(),
            error,
        })?;
        let mut text = Vec::new();
        match &id {
            ParsedSeqId::Gi(_) => continue,
            ParsedSeqId::Pir(_) | ParsedSeqId::Prf(_) => label::fasta(&id, &mut text),
            _ => label::content(&id, true, &mut text),
        }
        out.push(text);
    }
    out.sort();
    out.dedup();
    Ok(out)
}

/// `CSeq_id(text, fParse_RawText | fParse_AnyLocal | fParse_PartialOK)`, the parse of a
/// list id (the list's and the lookup files' ids alike). The accession guide's warnings are
/// not posted here (NCBI posts them through the diagnostic stream, which the search
/// programs leave at its default).
pub fn parse_list_id(text: &[u8]) -> Result<ParsedSeqId, SeqIdError> {
    seq_id_from_text(text, RAW_TEXT | ANY_LOCAL | PARTIAL_OK, &mut |_| {})
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tokens(text: &str) -> Vec<String> {
        let mut warnings = Vec::new();
        read_memory_si_list(text.as_bytes(), &mut warnings)
            .into_iter()
            .map(|t| String::from_utf8(t).unwrap())
            .collect()
    }

    #[test]
    fn text_tokens_as_ncbi_reads_them() {
        assert_eq!(
            tokens(">AAA00001.1\n# comment line\naaa00001.2 ref|NP_000001.1| # trailing\r\nlcl|x\tgnl|db|t\n"),
            ["AAA00001.1", "aaa00001.2", "ref|NP_000001.1|", "lcl|x", "gnl|db|t"]
        );
        assert_eq!(tokens("a#1 b"), ["a#1", "b"]);
        assert_eq!(tokens(">>  >id\r\n"), ["id"]);
        assert_eq!(tokens("P69905"), ["P69905"]);
        assert!(tokens("# only\n\n").is_empty());
        assert!(tokens("").is_empty());
    }

    #[test]
    fn a_token_of_vertical_tabs_is_reported_and_dropped() {
        let mut warnings = Vec::new();
        let ids = read_memory_si_list(b"\x0b\x0c ok\n", &mut warnings);
        assert_eq!(ids, [b"ok".to_vec()]);
        assert_eq!(
            warnings,
            ["WARNING:  \x0b\x0c is not a valid seqid string."]
        );
    }

    fn binary(ids: &[&[u8]], db_vol_length: u64) -> Vec<u8> {
        let mut body = Vec::new();
        body.extend_from_slice(&(ids.len() as u64).to_le_bytes());
        body.extend_from_slice(&5u32.to_le_bytes());
        body.extend_from_slice(b"title");
        body.push(3);
        body.extend_from_slice(b"Oct");
        body.extend_from_slice(&db_vol_length.to_le_bytes());
        if db_vol_length != 0 {
            body.push(2);
            body.extend_from_slice(b"Dt");
            body.extend_from_slice(&4u32.to_le_bytes());
            body.extend_from_slice(b"vols");
        }
        for id in ids {
            if id.len() >= 0xFF {
                body.push(0xFF);
                body.extend_from_slice(&(id.len() as u32).to_le_bytes());
            } else {
                body.push(id.len() as u8);
            }
            body.extend_from_slice(id);
        }
        let size = 1 + 8 + body.len();
        let mut file = vec![0u8];
        file.extend_from_slice(&(size as u64).to_le_bytes());
        file.extend_from_slice(&body);
        file
    }

    #[test]
    fn a_binary_list_gives_its_ids_and_header() {
        let long = vec![b'A'; 300];
        let file = binary(&[b"AAA00001.2", b"Q00002", &long], 12345);
        let mut warnings = Vec::new();
        let list = read_si_list(&file, &mut warnings).unwrap();
        assert_eq!(
            list.ids,
            [b"AAA00001.2".to_vec(), b"Q00002".to_vec(), long.clone()]
        );
        assert_eq!(
            list.info,
            SeqidListInfo {
                is_v4: false,
                file_size: file.len() as u64,
                num_ids: 3,
                title: b"title".to_vec(),
                create_date: b"Oct".to_vec(),
                db_vol_length: 12345,
                db_create_date: b"Dt".to_vec(),
                db_vol_names: b"vols".to_vec(),
            }
        );
        assert!(warnings.is_empty());
        let plain = binary(&[b"X"], 0);
        let list = read_si_list(&plain, &mut warnings).unwrap();
        assert_eq!(list.ids, [b"X".to_vec()]);
        assert_eq!(list.info.db_vol_length, 0);
        assert!(list.info.db_vol_names.is_empty());
    }

    #[test]
    fn damaged_binary_lists_give_ncbis_messages() {
        let mut warnings = Vec::new();
        assert_eq!(
            read_si_list(b"", &mut warnings).unwrap_err(),
            "Failed to map seqidlist file "
        );
        let mut file = binary(&[b"X"], 0);
        file.push(0);
        assert_eq!(
            read_si_list(&file, &mut warnings).unwrap_err(),
            "Invalid seqidlist file"
        );
        let mut short = binary(&[b"X", b"Y"], 0);
        let cut = short.len() - 2;
        short.truncate(cut);
        short[1..9].copy_from_slice(&(cut as u64).to_le_bytes());
        assert_eq!(
            read_si_list(&short, &mut warnings).unwrap_err(),
            "Invalid total num of ids in seqidlist file"
        );
        // A binary list of no ids is read as text: its header bytes are the tokens.
        let empty = binary(&[], 0);
        let list = read_si_list(&empty, &mut warnings).unwrap();
        assert!(list.info.is_v4);
        assert!(!list.ids.is_empty());
    }

    fn processed(tokens: &[&str]) -> Vec<String> {
        let ids: Vec<Vec<u8>> = tokens.iter().map(|t| t.as_bytes().to_vec()).collect();
        process_for_v5(&ids)
            .unwrap()
            .into_iter()
            .map(|t| String::from_utf8(t).unwrap())
            .collect()
    }

    #[test]
    fn version5_ids_are_rewritten_as_ncbi_rewrites_them() {
        assert_eq!(
            processed(&[
                "aaa00001.2",
                "ref|NP_000001.1|",
                "np_000002",
                "lcl|myprot9",
                "gnl|mydb|tag1",
                "mydb:tag1",
                "12345",
                "gi|999",
                "sp|Q00001.1|N1_HUMAN",
                "MixedCase1",
                "pdb|1ABC|A",
                "1ABC_A",
                "1ABC",
                "pir||S00007",
                "AAA00001.2",
            ]),
            [
                "12345",
                "1ABC",
                "1ABC_A",
                "AAA00001.2",
                "MixedCase1",
                "NP_000001.1",
                "NP_000002",
                "Q00001.1",
                "mydb:tag1",
                "myprot9",
                "pir||S00007",
            ]
        );
        assert!(processed(&["gi|999", "gi|998"]).is_empty());
        // An unknown FASTA tag is a local id under `fParse_AnyLocal` (probe S23 of
        // dropin-survey-20261029/DI8: no error, the id is looked up as it is).
        assert_eq!(processed(&["xyz|ABC|DEF"]), ["xyz|ABC|DEF"]);
        let bad = process_for_v5(&[b"gi|abc".to_vec()]).unwrap_err();
        assert_eq!(bad.message(), "Invalid seq id: gi|abc");
    }
}
