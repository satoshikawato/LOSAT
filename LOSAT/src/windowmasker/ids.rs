//! `-ids` and `-exclude_ids`: the id lists of `CWinMaskConfig::FillIdList`, matched as text
//! (`CWinMaskUtil::CIdSet_TextMatch`) against the first word of each sequence's title
//! (`sequence::GetTitle`), and `CWinMaskUtil::consider`.

use std::collections::HashSet;

use crate::seqid::ParsedSeqId;

/// `CWinMaskUtil::CIdSet_TextMatch`: the ids of a list, by their number of `|`-words
/// (`nword_sets_`).
#[derive(Debug, Default)]
pub(super) struct TextMatchIds {
    nword_sets: Vec<HashSet<Vec<u8>>>,
}

/// An id of a list that NCBI cannot store: `CIdSet_TextMatch::insert` posts `bad id` and then
/// indexes its sets at -1 (undefined behaviour in C++; NCBI BLAST+ 2.17.0 dies with SIGSEGV).
#[derive(Debug, PartialEq, Eq)]
pub(super) struct EmptyId {
    pub(super) id: Vec<u8>,
}

/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_util.cpp:75-100
/// ```c++
/// const vector< Uint4 >
/// CWinMaskUtil::CIdSet_TextMatch::split( const string & id_str )
/// {
///     vector< Uint4 > result;
///     string tmp = id_str;
///
///     if( !tmp.empty() && tmp[tmp.length() - 1] == '|' )
///         tmp = tmp.substr( 0, tmp.length() - 1 );
///
///     if( !tmp.empty() ) {
///         string::size_type pos = ( tmp[0] == '>' ) ? 1 : 0;
///         string::size_type len = tmp.length();
///
///         while( pos != string::npos && pos < len ) {
///             result.push_back( pos );
///
///             if( (pos = tmp.find_first_of( "|", pos )) != string::npos ) {
///                 ++pos;
///             }
///         }
///     }
///
///     result.push_back( tmp.length() + 1 );
///     return result;
/// }
/// ```
fn split(id: &[u8]) -> Vec<usize> {
    let tmp = id.strip_suffix(b"|").unwrap_or(id);
    let mut result = Vec::new();
    if !tmp.is_empty() {
        let mut pos = usize::from(tmp[0] == b'>');
        while pos < tmp.len() {
            result.push(pos);
            match tmp[pos..].iter().position(|&byte| byte == b'|') {
                Some(bar) => pos += bar + 1,
                None => break,
            }
        }
    }
    result.push(tmp.len() + 1);
    result
}

impl TextMatchIds {
    /// `CIdSet_TextMatch::insert( id_str )`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_util.cpp:103-123
    /// ```c++
    ///     Uint4 nwords = split( id_str ).size() - 1;
    ///
    ///     if( nwords == 0 ) {
    ///         ERR_POST( Error
    ///             << "CWinMaskConfig::CIdSet_TextMatch::insert(): bad id: "
    ///             << id_str << ": ignoring");
    ///     }
    ///
    ///     if( nword_sets_.size() < nwords ) {
    ///         nword_sets_.resize( nwords );
    ///     }
    ///
    ///     if( id_str[id_str.length() - 1] != '|' ) {
    ///         nword_sets_[nwords - 1].insert( id_str );
    ///     }else {
    ///         nword_sets_[nwords - 1].insert(
    ///                 id_str.substr( 0, id_str.length() - 1 ) );
    ///     }
    /// ```
    /// An id of no words (empty, `>`, `|`, `>|`) goes on to `nword_sets_[-1]`: `EmptyId`.
    pub(super) fn insert(&mut self, id: &[u8]) -> Result<(), EmptyId> {
        let nwords = split(id).len() - 1;
        if nwords == 0 {
            return Err(EmptyId { id: id.to_vec() });
        }
        if self.nword_sets.len() < nwords {
            self.nword_sets.resize_with(nwords, HashSet::new);
        }
        let stored = id.strip_suffix(b"|").unwrap_or(id);
        self.nword_sets[nwords - 1].insert(stored.to_vec());
        Ok(())
    }

    /// `empty()`: no id was stored.
    fn is_empty(&self) -> bool {
        self.nword_sets.is_empty()
    }

    /// `CIdSet_TextMatch::find( id_str )`: whether a run of consecutive words of `id` is an id
    /// of the list with as many words.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_util.cpp:154-175
    /// ```c++
    ///     vector< Uint4 > word_starts = split( id_str );
    ///
    ///     for( Uint4 i = 0;
    ///             i < nword_sets_.size() && i < word_starts.size() - 1; ++i ) {
    ///         if( !nword_sets_[i].empty() ) {
    ///             for( Uint4 j = 0; j < word_starts.size() - i - 1; ++j ) {
    ///                 string pattern = id_str.substr(
    ///                         word_starts[j],
    ///                         word_starts[j + i + 1] - word_starts[j] - 1 );
    ///
    ///                 if( find( pattern, i ) ) {
    ///                     return true;
    ///                 }
    ///             }
    ///         }
    ///     }
    /// ```
    fn find_text(&self, id: &[u8]) -> bool {
        let starts = split(id);
        let words = starts.len() - 1;
        for i in 0..self.nword_sets.len().min(words) {
            let set = &self.nword_sets[i];
            if set.is_empty() {
                continue;
            }
            for j in 0..words - i {
                let from = starts[j].min(id.len());
                let len = starts[j + i + 1] - starts[j] - 1;
                let pattern = &id[from..(from + len).min(id.len())];
                if set.contains(pattern) {
                    return true;
                }
            }
        }
        false
    }

    /// `CIdSet_TextMatch::find( bsh )` of a sequence whose `GetTitle` is `title`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_util.cpp:126-145
    /// ```c++
    ///     CConstRef< CBioseq > seq = bsh.GetCompleteBioseq();
    ///     string id_str = sequence::GetTitle( bsh );
    ///
    ///     if( !id_str.empty() ) {
    ///         string::size_type pos = id_str.find_first_of( " \t" );
    ///         id_str = id_str.substr( 0, pos );
    ///     }
    ///
    ///     if( find( id_str ) ) return true;
    ///     else if( id_str.substr( 0, 4 ) == "lcl|" ) {
    ///         id_str = id_str.substr( 4, string::npos );
    ///         return find( id_str );
    ///     }
    ///     else return false;
    /// ```
    fn find(&self, title: &[u8]) -> bool {
        let end = title
            .iter()
            .position(|&byte| byte == b' ' || byte == b'\t')
            .unwrap_or(title.len());
        let id = &title[..end];
        if self.find_text(id) {
            return true;
        }
        match id.strip_prefix(b"lcl|") {
            Some(rest) => self.find_text(rest),
            None => false,
        }
    }
}

/// `CWinMaskConfig::FillIdList( file_name, id_list )`: the ids of the file `name`, one per
/// line; a file that does not open gives none.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:471-484
/// ```c++
///     CNcbiIfstream file( file_name.c_str() );
///     string line;
///
///     while( NcbiGetlineEOL( file, line ) ) {
///         if( !line.empty() )
///         {
///             string::size_type stop( line.find_first_of( " \t" ) );
///             string::size_type start( line[0] == '>' ? 1 : 0 );
///             string id_str = line.substr( start, stop - start );
///             id_list.insert( id_str );
///         }
///     }
/// ```
/// `NcbiGetlineEOL` ends a line at `\n`, `\r` or `\r\n` (a `\r` alone ends one too), and the
/// last line needs no end. A line that starts with a space or a tab gives an empty id; one
/// whose first space comes at 0 after `>` (`> x`) too.
pub(super) fn fill_id_list(name: &str) -> Result<TextMatchIds, EmptyId> {
    let mut ids = TextMatchIds::default();
    let Ok(bytes) = std::fs::read(name) else {
        return Ok(ids);
    };
    for line in getline_eol(&bytes) {
        if line.is_empty() {
            continue;
        }
        let stop = line
            .iter()
            .position(|&byte| byte == b' ' || byte == b'\t')
            .unwrap_or(line.len());
        let start = usize::from(line[0] == b'>');
        let id = if stop >= start {
            &line[start..stop]
        } else {
            &[][..]
        };
        ids.insert(id)?;
    }
    Ok(ids)
}

/// The lines of `NcbiGetlineEOL`: ended by `\n`, `\r\n` or a `\r` alone.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:191-224
/// ```c++
/// CNcbiIstream& NcbiGetlineEOL(CNcbiIstream& is, string& str, string::size_type* count)
/// {
///     return s_NcbiGetline(is, str, "\r\n", count);
/// }
/// ```
/// (`s_NcbiGetline` stops at either delimiter and, after a `\r`, takes a `\n` that follows as
/// part of the same end.)
fn getline_eol(bytes: &[u8]) -> Vec<&[u8]> {
    let mut lines = Vec::new();
    let mut rest = bytes;
    while !rest.is_empty() {
        match rest.iter().position(|&byte| byte == b'\n' || byte == b'\r') {
            Some(end) => {
                lines.push(&rest[..end]);
                let mut next = end + 1;
                if rest[end] == b'\r' && rest.get(next) == Some(&b'\n') {
                    next += 1;
                }
                rest = &rest[next..];
            }
            None => {
                lines.push(rest);
                rest = &[];
            }
        }
    }
    lines
}

/// The `-ids` and `-exclude_ids` lists of the run (`CWinMaskConfig::Ids`, `ExcludeIds`).
#[derive(Debug, Default)]
pub(super) struct IdLists {
    pub(super) ids: Option<TextMatchIds>,
    pub(super) exclude_ids: Option<TextMatchIds>,
}

impl IdLists {
    /// Whether a list asks for the titles (`consider` returns at once without one).
    pub(super) fn active(&self) -> bool {
        let used = |list: &Option<TextMatchIds>| list.as_ref().is_some_and(|ids| !ids.is_empty());
        used(&self.ids) || used(&self.exclude_ids)
    }

    /// `CWinMaskUtil::consider( bsh, ids, exclude_ids )` of a sequence whose `GetTitle` is
    /// `title` (asked only when a list is not empty).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_util.cpp:218-245
    /// ```c++
    ///     if( (ids == 0 || ids->empty()) &&
    ///         (exclude_ids == 0 || exclude_ids->empty()) ) {
    ///         return true;
    ///     }
    ///
    ///     bool result = true;
    ///
    ///     if( ids != 0 && !ids->empty() )
    ///     {
    ///         result = false;
    ///
    ///         if( ids->find( bsh ) ) {
    ///             result = true;
    ///         }
    ///     }
    ///
    ///     if( exclude_ids != 0 && !exclude_ids->empty() )
    ///     {
    ///         if( exclude_ids->find( bsh ) ) {
    ///             result = false;
    ///         }
    ///     }
    ///
    ///     return result;
    /// ```
    pub(super) fn consider(&self, title: &[u8]) -> bool {
        let mut result = true;
        if let Some(ids) = self.ids.as_ref().filter(|ids| !ids.is_empty()) {
            result = ids.find(title);
        }
        if let Some(exclude) = self.exclude_ids.as_ref().filter(|ids| !ids.is_empty()) {
            if exclude.find(title) {
                result = false;
            }
        }
        result
    }
}

/// What `sequence::GetTitle( bsh )` gives for a sequence of the FASTA or database reader, or
/// the case LOSAT does not reproduce.
#[derive(Debug, PartialEq, Eq)]
pub(super) enum Title {
    Text(Vec<u8>),
    /// An empty title of a sequence with a source (a database sequence with a taxid), which NCBI
    /// makes from the source's names (`s_TitleFromBioSource`, or "No definition line found").
    FromSource,
}

/// `sequence::GetTitle( bsh )` (flags 0) of a sequence with the Seq-ids `ids`, the title
/// descriptor `title` and a source when `has_source`: the readers' Bioseqs have no MolInfo,
/// no GenBank or EMBL block and no PDB descriptor, and are nucleotides of the raw or delta
/// representation, so of the function only these steps apply.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/seqtitle.cpp:131-178,286-296,371-412
/// ```c++
///     ITERATE (CBioseq_Handle::TId, idh, hnd.GetId()) {
///         ...
///         case CSeq_id::e_Tpg:
///         case CSeq_id::e_Tpe:
///         case CSeq_id::e_Tpd:
///             third_party = true;
///             break;
///         ...
///         case CSeq_id::e_Patent:
///             pat_id = &id->GetPatent();
///             break;
/// ...
///         CSeqdesc_CI it(hnd, CSeqdesc::e_Title, search_depth);
///         if (it) {
///             title = it->GetTitle();
///         }
/// ...
///     // originally further down, but moved up to match the C version
///     while (NStr::EndsWith(title, ".")  ||  NStr::EndsWith(title, " ")) {
///         title.erase(title.end() - 1);
///     }
/// ...
///     if (title.empty()  &&  pat_id.NotEmpty()) {
///         title = "Sequence " + NStr::IntToString(pat_id->GetSeqid())
///             + " from Patent " + pat_id->GetCit().GetCountry()
///             + ' ' + pat_id->GetCit().GetSomeNumber();
///     }
/// ...
///     if (title.empty()  &&  !htg_tech  &&  source.NotEmpty()) {
///         title = s_TitleFromBioSource(*source, tech);
///         if (title.empty()) {
///             title = "No definition line found";
///         }
///     }
///
///     if (is_tsa  &&  !title.empty() ) {
///         prefix = "TSA: ";
///     } else if (third_party  &&  !title.empty() ) {
///         bool tpa_start = NStr::StartsWith(title, "TPA: ", NStr::eNocase);
///         if (tpa_exp) {
///         ...
///         } else if ( !tpa_start ) {
///             prefix = "TPA: ";
///         }
///     }
/// ...
///     return prefix + title + suffix;
/// ```
/// (The steps for an empty title from the source of an NC, NM or NR accession, from a PDB
/// descriptor or from a protein do not apply: they need a MolInfo, a feature or a descriptor
/// the readers do not make, or a protein. The flags add no organism.)
pub(super) fn get_title(ids: &[ParsedSeqId], title: Option<&[u8]>, has_source: bool) -> Title {
    let mut third_party = false;
    let mut patent = None;
    for id in ids {
        match id {
            ParsedSeqId::Tpg(_) | ParsedSeqId::Tpe(_) | ParsedSeqId::Tpd(_) => third_party = true,
            ParsedSeqId::Patent(pat) => patent = Some(pat),
            _ => {}
        }
    }
    let mut text = title.unwrap_or_default().to_vec();
    while matches!(text.last(), Some(b'.' | b' ')) {
        text.pop();
    }
    if text.is_empty() {
        if let Some(pat) = patent {
            let number = match &pat.cit.id {
                crate::seqid::model::PatentNumber::Number(number)
                | crate::seqid::model::PatentNumber::AppNumber(number) => number,
            };
            text = format!("Sequence {} from Patent ", pat.seqid).into_bytes();
            text.extend_from_slice(&pat.cit.country);
            text.push(b' ');
            text.extend_from_slice(number);
        }
    }
    if text.is_empty() && has_source {
        return Title::FromSource;
    }
    if third_party && !text.is_empty() && !starts_with_nocase(&text, b"TPA: ") {
        let mut prefixed = b"TPA: ".to_vec();
        prefixed.extend_from_slice(&text);
        text = prefixed;
    }
    Title::Text(text)
}

/// `NStr::StartsWith(s, prefix, NStr::eNocase)`.
fn starts_with_nocase(text: &[u8], prefix: &[u8]) -> bool {
    text.len() >= prefix.len() && text[..prefix.len()].eq_ignore_ascii_case(prefix)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn list(ids: &[&str]) -> TextMatchIds {
        let mut set = TextMatchIds::default();
        for id in ids {
            set.insert(id.as_bytes()).unwrap();
        }
        set
    }

    #[test]
    fn split_gives_the_word_starts_and_the_end() {
        assert_eq!(split(b"gb|AB000001.1|"), vec![0, 3, 14]);
        assert_eq!(split(b">rec"), vec![1, 5]);
        assert_eq!(split(b"a||b"), vec![0, 2, 3, 5]);
        assert_eq!(split(b""), vec![1]);
        assert_eq!(split(b"|"), vec![1]);
        assert_eq!(split(b">"), vec![2]);
    }

    #[test]
    fn ids_of_no_words_are_refused() {
        let mut set = TextMatchIds::default();
        for id in ["", ">", "|", ">|"] {
            assert_eq!(
                set.insert(id.as_bytes()),
                Err(EmptyId {
                    id: id.as_bytes().to_vec()
                })
            );
        }
    }

    #[test]
    fn runs_of_words_of_the_first_word_match() {
        let set = list(&["AB000001.1", "P12345|ABC_HUMAN", "NC_000099.1|", "rec2"]);
        assert!(set.find(b"gb|AB000001.1| GenBank style"));
        assert!(set.find(b"sp|P12345|ABC_HUMAN protein-like id"));
        assert!(set.find(b"ref|NC_000099.1|"));
        assert!(set.find(b"lcl|rec2 local"));
        assert!(!set.find(b"rec22"));
        assert!(!set.find(b"x AB000001.1"));
        let two = list(&["lcl|rec2"]);
        assert!(two.find(b"lcl|rec2 local"));
        assert!(!two.find(b"rec2"));
    }

    #[test]
    fn the_lines_end_at_any_line_end() {
        assert_eq!(
            getline_eol(b"a\r\nb\rc\n\nd"),
            vec![&b"a"[..], b"b", b"c", b"", b"d"]
        );
    }

    #[test]
    fn titles_lose_trailing_dots_and_spaces_and_third_party_ids_get_tpa() {
        let tpg = ParsedSeqId::Tpg(Default::default());
        assert_eq!(
            get_title(&[], Some(b"rec4. ."), false),
            Title::Text(b"rec4".to_vec())
        );
        assert_eq!(
            get_title(std::slice::from_ref(&tpg), Some(b"third party"), false),
            Title::Text(b"TPA: third party".to_vec())
        );
        assert_eq!(
            get_title(std::slice::from_ref(&tpg), Some(b"tpa: x"), false),
            Title::Text(b"tpa: x".to_vec())
        );
        assert_eq!(get_title(&[tpg], None, false), Title::Text(Vec::new()));
        assert_eq!(get_title(&[], Some(b"."), true), Title::FromSource);
    }
}
