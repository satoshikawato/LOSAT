//! `-taxid` and `-taxid_map`: the tax ids makeblastdb puts on the definition lines
//! (`CTaxIdSet`, `c++/src/objtools/blast/seqdb_writer/taxid_set.cpp`) and the keys of the map
//! (`AccessionToKey`, `GetSeqIdKey`, `multisource_util.cpp`).

use std::io::BufRead;

use rustc_hash::FxHashMap;

use crate::seqid::accguide::{parse_flags, SeqIdWarning};
use crate::seqid::parse::{is_space, seq_id_from_text, truncate_spaces};
use crate::seqid::{label, BlastDefline, DeflineSet, ParsedSeqId};
use crate::utils::winmask::NCBI_BUILD_ROOT;

use super::data_size::printable_string;

/// `CTaxIdSet`.
#[derive(Debug, Clone)]
pub struct TaxIdSet {
    /// `m_GlobalTaxId` (`-taxid`); 0 is `kTaxIdNotSet`.
    global_taxid: i32,
    /// `m_TaxIdMap` (`-taxid_map`).
    map: FxHashMap<Vec<u8>, i32>,
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
    /// An `NCBI C++ Exception` makeblastdb reports with "Error: " and exit status 255: the
    /// report (`CException::what()`, ending with a line feed).
    Exception(String),
    /// A line of white space only (1-based number): makeblastdb 2.17.0 reads the first token
    /// of an empty token list and crashes (segmentation fault); LOSAT refuses the file.
    BlankLine(usize),
}

/// A `CStringException` of `NStr`'s number conversions, as `CException::what()` prints it.
fn string_exception(line: u32, function: &str, message: String, pos: usize) -> String {
    format!(
        "NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/corelib/ncbistr.cpp\", line {line}: \
         Error: (CStringException::eConvert) ncbi::NStr::{function}() - {message} (m_Pos = {pos})\n"
    )
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:798-871
// ```c++
// Int8 NStr::StringToInt8(const CTempString str, TStringToNumFlags flags, int base)
// {
//     S2N_CONVERT_GUARD(flags);
//
//     // Current position in the string
//     SIZE_TYPE pos = 0;
//
//     // Skip allowed leading symbols
//     if (flags & fAllowLeadingSymbols) {
//         bool spaces = ((flags & fAllowLeadingSymbols) == fAllowLeadingSpaces);
//         s_SkipAllowedSymbols(str, pos, spaces ? eSkipSpacesOnly : eSkipAllAllowed, flags);
//     }
//     // Determine sign
//     bool sign = false;
//     switch (str[pos]) {
//     case '-':
//         sign = true;
//         /*FALLTHRU*/
//     case '+':
//         pos++;
//         break;
//     ...
//     }
//     SIZE_TYPE pos0 = pos;
//     ...
//     Int8 n = 0;
//     Int8 limdiv = base==10? kMax_I8 / 10: kMax_I8 / base;
//     Int8 limoff = (base==10? kMax_I8 % 10: kMax_I8 % base) + (sign ? 1 : 0);
//     ...
//     while (char ch = str[pos]) {
//         int  delta;   // corresponding numeric value of 'ch'
//         ...
//         if ( !s_IsGoodCharForRadix(ch, base, &delta) ) {
//             break;
//         }
//         // Overflow check
//         if ( n >= limdiv  &&  (n > limdiv  ||  delta > limoff) ) {
//             S2N_CONVERT_ERROR_OVERFLOW(Int8);
//         }
//         n *= base;
//         n += delta;
//         pos++;
//     }
//
//     // Last checks
//     if ( pos == pos0  || ((comma >= 0)  &&  (comma != 3)) ) {
//         S2N_CONVERT_ERROR_INVAL(Int8);
//     }
//     ...
//     // Assign sign before the end pointer check
//     n = sign ? -n : n;
//     CHECK_ENDPTR(Int8);
//
//     return n;
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:740-759
// ```c++
// void s_SkipAllowedSymbols(const CTempString       str,
//                           SIZE_TYPE&              pos,
//                           ESkipMode               skip_mode,
//                           NStr::TStringToNumFlags flags)
// {
//     ...
//     for ( SIZE_TYPE len = str.length();  pos < len;  ++pos ) {
//         unsigned char ch = str[pos];
//         if ( isdigit(ch)  ||  ch == '+' ||  ch == '-'  ||  s_IsDecimalPoint(ch,flags) ) {
//             break;
//         }
//         if ( (skip_mode == eSkipSpacesOnly)  &&  !isspace(ch) ) {
//             break;
//         }
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:546-557
// ```c++
//     s += "Cannot convert string '";
//     s += NStr::PrintableString(str);
//     s += "' to ";
//     s += to_type;
//     if ( !msg.empty() ) {
//         s += ", ";
//         s += msg;
//     }
// ```
// The string is a `std::string`'s: the loop ends at its terminating NUL, and at an
// embedded one (whose end check then passes). The errors are thrown at lines 852
// (overflow), 861 (no digit) and 870 (`CHECK_ENDPTR`).
/// `NStr::StringToInt8(str, flags)` in base 10, `leading_spaces` for
/// `fAllowLeadingSpaces`: the number, or the report of NCBI's exception.
fn string_to_int8(s: &[u8], leading_spaces: bool) -> Result<i64, String> {
    let at = |pos: usize| s.get(pos).copied().unwrap_or(0);
    let error = |line: u32, extra: &str, pos: usize| {
        string_exception(
            line,
            "StringToInt8",
            format!(
                "Cannot convert string '{}' to Int8{extra}",
                printable_string(s)
            ),
            pos,
        )
    };
    let mut pos = 0;
    if leading_spaces {
        while pos < s.len() {
            let ch = s[pos];
            if ch.is_ascii_digit() || ch == b'+' || ch == b'-' || !is_space(ch) {
                break;
            }
            pos += 1;
        }
    }
    let sign = match at(pos) {
        b'-' => {
            pos += 1;
            true
        }
        b'+' => {
            pos += 1;
            false
        }
        _ => false,
    };
    let pos0 = pos;
    let limdiv = i64::MAX / 10;
    let limoff = i64::MAX % 10 + i64::from(sign);
    let mut n: i64 = 0;
    loop {
        let ch = at(pos);
        if ch == 0 || !ch.is_ascii_digit() {
            break;
        }
        let delta = i64::from(ch - b'0');
        if n >= limdiv && (n > limdiv || delta > limoff) {
            return Err(error(852, ", overflow", pos));
        }
        // `i64::MIN` passes the check as its magnitude, one past `i64::MAX`.
        n = n.wrapping_mul(10).wrapping_add(delta);
        pos += 1;
    }
    if pos == pos0 {
        return Err(error(861, "", pos));
    }
    let n = if sign { n.wrapping_neg() } else { n };
    if at(pos) != 0 {
        return Err(error(870, "", pos));
    }
    Ok(n)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:635-643
// ```c++
// int NStr::StringToInt(const CTempString str, TStringToNumFlags flags, int base)
// {
//     S2N_CONVERT_GUARD_EX(flags);
//     Int8 value = StringToInt8(str, flags, base);
//     if ( value < kMin_Int  ||  value > kMax_Int ) {
//         S2N_CONVERT_ERROR(int, "overflow", ERANGE, 0);
//     }
//     return (int) value;
// }
// ```
/// `NStr::StringToNumeric<TTaxId>(str, fAllowLeadingSpaces)` (`TTaxId` is `int`:
/// `StringToInt`).
fn string_to_taxid(s: &[u8]) -> Result<i32, String> {
    let value = string_to_int8(s, true)?;
    i32::try_from(value).map_err(|_| {
        string_exception(
            640,
            "StringToInt",
            format!(
                "Cannot convert string '{}' to int, overflow",
                printable_string(s)
            ),
            0,
        )
    })
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
) -> Result<Option<Checked>, String> {
    if !acc.is_empty() && acc.iter().all(u8::is_ascii_digit) {
        return string_to_int8(acc, false).map(|gi| Some(Checked::Gi(gi)));
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
pub fn accession_to_key(acc: &[u8], warn: &mut dyn FnMut(SeqIdWarning)) -> Result<Vec<u8>, String> {
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
            map: FxHashMap::default(),
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
    /// `SetMappingFromFile`: the map of a `-taxid_map` file. `warn` takes the warnings of
    /// the Seq-id parser.
    pub fn set_mapping_from_file(
        &mut self,
        file: &mut dyn BufRead,
        warn: &mut dyn FnMut(SeqIdWarning),
    ) -> Result<(), TaxIdMapError> {
        let mut line = Vec::new();
        let mut number = 0usize;
        loop {
            line.clear();
            let read = file
                .read_until(b'\n', &mut line)
                .map_err(|error| TaxIdMapError::Exception(error.to_string()))?;
            if read == 0 {
                break;
            }
            number += 1;
            if line.last() == Some(&b'\n') {
                line.pop();
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
                let taxid = string_to_taxid(tx_str).map_err(TaxIdMapError::Exception)?;
                let key = accession_to_key(gi_str, warn).map_err(TaxIdMapError::Exception)?;
                self.map.insert(key, taxid);
            }
        }
        self.matched = self.global_taxid != 0 || self.map.is_empty();
        Ok(())
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
            if let Some(&taxid) = self.map.get(&key) {
                self.matched = true;
                return taxid;
            }
            if let Some(dot) = key.iter().position(|&c| c == b'.') {
                if let Some(&taxid) = self.map.get(&key[..dot]) {
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
