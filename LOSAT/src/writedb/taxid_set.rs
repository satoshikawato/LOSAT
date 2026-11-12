//! `-taxid` and `-taxid_map`: the tax ids makeblastdb puts on the definition lines
//! (`CTaxIdSet`, `c++/src/objtools/blast/seqdb_writer/taxid_set.cpp`) and the keys of the map
//! (`AccessionToKey`, `GetSeqIdKey`, `multisource_util.cpp`).

use std::io::BufRead;

use crate::ncbi_exception::NcbiException;
use crate::ncbi_str::{string_to_int, string_to_int8, ALLOW_LEADING_SPACES};
use crate::seqid::accguide::{parse_flags, SeqIdWarning};
use crate::seqid::parse::{seq_id_from_text, truncate_spaces};
use crate::seqid::{label, BlastDefline, DeflineSet, ParsedSeqId};

use super::id_set::KeyMap;

/// `CTaxIdSet`.
#[derive(Debug, Clone)]
pub struct TaxIdSet {
    /// `m_GlobalTaxId` (`-taxid`); 0 is `kTaxIdNotSet`.
    global_taxid: i32,
    /// `m_TaxIdMap` (`-taxid_map`): the keys in an arena (`KeyMap`, about the key's bytes +
    /// 16 to 26 bytes per entry; a `map<string, TTaxId>` node and string in NCBI).
    map: KeyMap,
    /// `m_Matched`.
    matched: bool,
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/taxid_set.hpp:51-53
// ```c++
//     CTaxIdSet(TTaxId global_taxid = kTaxIdNotSet)
//         : m_GlobalTaxId(global_taxid),
//           m_Matched(true) {}
// ```
impl Default for TaxIdSet {
    fn default() -> Self {
        Self::with_global(0)
    }
}

/// Why a `-taxid_map` file could not be read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TaxIdMapError {
    /// An exception of NCBI's `NStr` that makeblastdb reports with "Error: " and its report
    /// (`CException::what()`), exit status 255.
    Exception(NcbiException),
    /// A line of white space only (1-based number): makeblastdb 2.17.0 reads the first token
    /// of an empty token list and crashes (segmentation fault); LOSAT refuses the file.
    BlankLine(usize),
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/multisource_util.cpp:157-160
// ```c++
// void GetSeqIdKey(const objects::CSeq_id & id, string & key)
// {
//     id.GetLabel(& key, CSeq_id::eContent, CSeq_id::fLabel_GeneralDbIsContent);
// }
// ```
/// `GetSeqIdKey`.
pub fn seq_id_key(id: &ParsedSeqId) -> Vec<u8> {
    let mut key = Vec::new();
    label::content_general_db_is_content(id, &mut key);
    key
}

/// What `CheckAccession` found: a GI, or another Seq-id.
enum Checked {
    Gi(i64),
    SeqId(ParsedSeqId),
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/multisource_util.cpp:108-155
// ```c++
// bool CheckAccession(const string  & acc,
//                     TGi           & gi,
//                     CRef<objects::CSeq_id> & seqid,
//                     bool          & specific)
// {
//     specific = true;
//     gi = ZERO_GI;
//     seqid.Reset();
//
//     // Case 1: Numeric GI
//
//     bool digits = !! acc.size();
//
//     for(unsigned i = 0; i < acc.size(); i++) {
//         if (! isdigit(acc[i])) {
//             digits = false;
//             break;
//         }
//     }
//
//     if (digits) {
//         gi = NStr::StringToNumeric<TGi>(acc);
//         return true;
//     }
//
//     try {
//         seqid.Reset(new CSeq_id(acc));
//     }
//     catch(CException &) {
//         return false;
//     }
//
//     if (seqid->IsGi()) {
//         gi = seqid->GetGi();
//         seqid.Reset();
//         return true;
//     }
//     ...
//     return true;
// }
// ```
// `CSeq_id(acc)` parses with its default flags, `fParse_AnyRaw` (Seq_id.hpp:135-136); a
// `TGi` is an `Int8`, whose conversion of too many digits throws a `CStringException` that
// nothing catches (makeblastdb: exit status 255).
/// `CheckAccession`: `Ok(None)` when `acc` is not a Seq-id.
fn check_accession(
    acc: &[u8],
    warn: &mut dyn FnMut(SeqIdWarning),
) -> Result<Option<Checked>, NcbiException> {
    if !acc.is_empty() && acc.iter().all(u8::is_ascii_digit) {
        return string_to_int8(acc, 0).map(|gi| Some(Checked::Gi(gi)));
    }
    Ok(match seq_id_from_text(acc, parse_flags::ANY_RAW, warn) {
        Err(_) => None,
        Ok(ParsedSeqId::Gi(gi)) => Some(Checked::Gi(gi)),
        Ok(id) => Some(Checked::SeqId(id)),
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/multisource_util.cpp:82-106
// ```c++
// string AccessionToKey(const string & acc)
// {
//     TGi gi(ZERO_GI);
//     CRef<CSeq_id> seqid;
//     bool specific(false);
//
//     string str;
//
//     if (CheckAccession(acc, gi, seqid, specific)) {
//         if (seqid.Empty()) {
//             if (gi != ZERO_GI) {
//                 str = NStr::NumericToString(gi);
//             }
//         } else {
//             GetSeqIdKey(*seqid, str);
//         }
//     } else {
//         string tempID = "lcl|" + acc;
//         if(CheckAccession(tempID, gi, seqid, specific)) {
//             GetSeqIdKey(*seqid, str);
//         }
//     }
//
//     return str;
// }
// ```
/// `AccessionToKey`: the key of a map line's Seq-id, empty for none.
pub fn accession_to_key(
    acc: &[u8],
    warn: &mut dyn FnMut(SeqIdWarning),
) -> Result<Vec<u8>, NcbiException> {
    let gi_key = |gi: i64| {
        if gi != 0 {
            gi.to_string().into_bytes()
        } else {
            Vec::new()
        }
    };
    Ok(match check_accession(acc, warn)? {
        Some(Checked::Gi(gi)) => gi_key(gi),
        Some(Checked::SeqId(id)) => seq_id_key(&id),
        None => {
            let mut temp_id = b"lcl|".to_vec();
            temp_id.extend_from_slice(acc);
            match check_accession(&temp_id, warn)? {
                Some(Checked::SeqId(id)) => seq_id_key(&id),
                // `GetSeqIdKey(*seqid)` of a GI: the retried id is a local id.
                Some(Checked::Gi(gi)) => gi.to_string().into_bytes(),
                None => Vec::new(),
            }
        }
    })
}

impl TaxIdSet {
    /// `CTaxIdSet(global_taxid)`.
    pub fn with_global(global_taxid: i32) -> Self {
        Self {
            global_taxid,
            map: KeyMap::default(),
            matched: true,
        }
    }

    /// `m_GlobalTaxId`.
    pub fn global_taxid(&self) -> i32 {
        self.global_taxid
    }

    /// `HasEverFixedId`.
    pub fn has_ever_fixed_id(&self) -> bool {
        self.matched
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/taxid_set.cpp:45-74
    // ```c++
    // void CTaxIdSet::SetMappingFromFile(CNcbiIstream & f)
    // {
    //     while(f && (! f.eof())) {
    //         string s;
    //         NcbiGetlineEOL(f, s);
    //
    //         if (s.empty())
    //             continue;
    //
    //         // Remove leading/trailing spaces.
    //         s = NStr::TruncateSpaces(s);
    //
    //         vector<string> tokens;
    //         NStr::Split(s, " \t", tokens);
    //
    //         string gi_str = tokens.front();
    //         string tx_str;
    //         if (tokens.size() == 2) {
    //             tx_str = tokens.back();
    //         }
    //
    //         if (gi_str.size() && tx_str.size()) {
    //             TTaxId taxid = NStr::StringToNumeric<TTaxId>(tx_str, NStr::fAllowLeadingSpaces);
    //             string key = AccessionToKey(gi_str);
    //
    //             m_TaxIdMap[key] = taxid;
    //         }
    //     }
    //     m_Matched = (m_GlobalTaxId != kTaxIdNotSet) || m_TaxIdMap.empty();
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:196-208
    // ```c++
    // CNcbiIstream& NcbiGetlineEOL(CNcbiIstream& is, string& str, SIZE_TYPE* count)
    // {
    // #if   defined(NCBI_OS_MSWIN)
    //     NcbiGetline(is, str, '\n', count);
    //     if (!str.empty()  &&  str[str.length() - 1] == '\r')
    //         str.resize(str.length() - 1);
    // #elif defined(NCBI_OS_DARWIN)
    //     NcbiGetline(is, str, "\r\n", count);
    // #else /* assume UNIX-like EOLs */
    //     NcbiGetline(is, str, '\n', count);
    // #endif //NCBI_OS_...
    //     return is;
    // }
    // ```
    // Lines end at '\n' (the Linux reading; a '\r' before it is white space that
    // `TruncateSpaces` removes, as on Windows). `NStr::Split` without flags merges no
    // delimiters: two blanks make an empty token, so the line has three tokens and is
    // ignored, as is a line of one token or of more than two. A line of white space only
    // leaves no token, and `tokens.front()` reads past an empty vector (makeblastdb 2.17.0
    // crashes); LOSAT refuses the file there.
    // A read error (a directory opened as the file, EIO) sets the stream's badbit: the line
    // read up to it is still taken (the loop's body does not test the stream), then the loop
    // ends and the map keeps what it has, with no message (makeblastdb 2.17.0 with a
    // directory as `-taxid_map`: an empty map, exit status 0). `read_until` keeps the bytes
    // read before its error in `line` as `NcbiGetline` keeps them in `s`.
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:153-176
    // ```c++
    //     do {
    //         CT_INT_TYPE nextc = is.get();
    //         if (CT_EQ_INT_TYPE(nextc, CT_EOF)
    //             ||  CT_EQ_INT_TYPE(nextc, CT_TO_INT_TYPE(delim))) {
    //             ++size;
    //             break;
    //         }
    //         if ( !is.unget() )
    //             break;
    //         ...
    //         is.get(buf, n < sizeof(buf) ? n : sizeof(buf), delim);
    //         n = (size_t) is.gcount();
    //         str.append(buf, n);
    //         size += n;
    //         _ASSERT(size == str.length());
    //     } while ( is.good() );
    // #endif
    //
    //     if (is.rdstate() == NcbiEofbit  &&  str.empty())
    //         is.setstate(NcbiFailbit);
    // ```
    // The stream's state at the end matters when it is the standard input that `-in -` reads
    // next (`AddFasta` reads only a stream that has not failed): a read that finds nothing
    // (the end after a line feed) or an error fails it; a last line without a line feed ends
    // at the end of the stream with only `eofbit` set, and the loop stops there.
    /// `SetMappingFromFile`: the map of a `-taxid_map` file. `warn` takes the warnings of
    /// the Seq-id parser. `Ok(true)` when the stream has failed at the end (`!f`).
    pub fn set_mapping_from_file(
        &mut self,
        file: &mut dyn BufRead,
        warn: &mut dyn FnMut(SeqIdWarning),
    ) -> Result<bool, TaxIdMapError> {
        let mut line = Vec::new();
        let mut number = 0usize;
        let mut failed = false;
        let mut at_end = false;
        while !failed && !at_end {
            line.clear();
            let read = match file.read_until(b'\n', &mut line) {
                Ok(read) => read,
                Err(_) => {
                    failed = true;
                    line.len()
                }
            };
            if read == 0 {
                failed = true;
                break;
            }
            number += 1;
            if line.last() == Some(&b'\n') {
                line.pop();
            } else {
                at_end = true;
            }
            if line.is_empty() {
                continue;
            }
            let s = truncate_spaces(&line);
            if s.is_empty() {
                return Err(TaxIdMapError::BlankLine(number));
            }
            let tokens: Vec<&[u8]> = s.split(|&c| c == b' ' || c == b'\t').collect();
            let gi_str = tokens[0];
            let tx_str: &[u8] = if tokens.len() == 2 { tokens[1] } else { b"" };
            if !gi_str.is_empty() && !tx_str.is_empty() {
                // NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:3582-3587
                // ```c++
                //     template <typename TNumeric>
                //     static typename enable_if< is_integral<TNumeric>::value && is_signed<TNumeric>::value && (sizeof(TNumeric) == sizeof(int) && !is_same<TNumeric, long>::value), TNumeric>::type
                //     x_StringToNumeric(const CTempString str, TStringToNumFlags flags, int base)
                //     {
                //         return StringToInt(str, flags, base);
                //     }
                // ```
                // `TTaxId` is an `int`: `StringToNumeric<TTaxId>` is `NStr::StringToInt`.
                let taxid = string_to_int(tx_str, ALLOW_LEADING_SPACES)
                    .map_err(TaxIdMapError::Exception)?;
                let key = accession_to_key(gi_str, warn).map_err(TaxIdMapError::Exception)?;
                self.map.set(&key, taxid);
            }
        }
        self.matched = self.global_taxid != 0 || self.map.is_empty();
        Ok(failed)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/taxid_set.cpp:84-128
    // ```c++
    // TTaxId CTaxIdSet::x_SelectBestTaxid(const objects::CBlast_def_line & defline)
    // {
    //     TTaxId retval = m_GlobalTaxId;
    //
    //     if (retval != kTaxIdNotSet) {
    //         return retval;
    //     }
    //
    //     if ( !m_TaxIdMap.empty() ) {
    //         vector<string> keys;
    //         GetDeflineKeys(defline, keys);
    //
    //         ITERATE(vector<string>, key, keys) {
    //             if (key->empty())
    //                 continue;
    //
    //             map<string, TTaxId>::const_iterator item = m_TaxIdMap.find(*key);
    //
    //             if (item != m_TaxIdMap.end()) {
    //                 retval = item->second;
    //                 m_Matched = true;
    //                 break;
    //             }
    //
    // 	    // try removing version to see if strings match.
    // 	    // this is most likely to help if the Query ID was parsed as local ID but is really accession.
    // 	    string keyNoVersion;
    // 	    string version;
    // 	    if (NStr::SplitInTwo(*key, ".", keyNoVersion, version) == true)
    // 	    {
    //             	item = m_TaxIdMap.find(keyNoVersion);
    //
    //             	if (item != m_TaxIdMap.end()) {
    //                	 retval = item->second;
    //                	 m_Matched = true;
    //                	 break;
    // 		}
    //             }
    //         }
    //     } else if (defline.IsSetTaxid()) {
    //         retval = defline.GetTaxid();
    //     }
    //
    //     return retval;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/multisource_util.cpp:162-172
    // ```c++
    // void GetDeflineKeys(const objects::CBlast_def_line & defline,
    //                     vector<string>        & keys)
    // {
    //     keys.clear();
    //
    //     ITERATE(CBlast_def_line::TSeqid, iter, defline.GetSeqid()) {
    //         string key;
    //         GetSeqIdKey(**iter, key);
    //         keys.push_back(key);
    //     }
    // }
    // ```
    /// `x_SelectBestTaxid`.
    fn select_best_taxid(&mut self, defline: &BlastDefline) -> i32 {
        if self.global_taxid != 0 {
            return self.global_taxid;
        }
        if self.map.is_empty() {
            return defline.taxid.unwrap_or(0);
        }
        for id in &defline.ids {
            let key = seq_id_key(id);
            if key.is_empty() {
                continue;
            }
            if let Some(taxid) = self.map.get(&key) {
                self.matched = true;
                return taxid;
            }
            if let Some(dot) = key.iter().position(|&c| c == b'.') {
                if let Some(taxid) = self.map.get(&key[..dot]) {
                    self.matched = true;
                    return taxid;
                }
            }
        }
        0
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/taxid_set.cpp:130-136
    // ```c++
    // void
    // CTaxIdSet::FixTaxId(CRef<objects::CBlast_def_line_set> deflines)
    // {
    //     NON_CONST_ITERATE(CBlast_def_line_set::Tdata, itr, deflines->Set()) {
    //         (*itr)->SetTaxid(x_SelectBestTaxid(**itr));
    //     }
    // }
    // ```
    /// `FixTaxId`: the taxid of every definition line.
    pub fn fix_tax_id(&mut self, deflines: &mut DeflineSet) {
        for defline in &mut deflines.deflines {
            defline.taxid = Some(self.select_best_taxid(defline));
        }
    }
}

#[cfg(test)]
#[path = "taxid_set_tests.rs"]
mod tests;
