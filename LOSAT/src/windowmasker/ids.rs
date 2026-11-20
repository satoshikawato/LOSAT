//! `-ids` and `-exclude_ids`: the id lists of `CWinMaskConfig::FillIdList`, matched as text
//! (`CWinMaskUtil::CIdSet_TextMatch`, the default) against the first word of each sequence's
//! title (`sequence::GetTitle`), or as Seq-ids (`CWinMaskUtil::CIdSet_SeqId`, `-text_match false`
//! with `-infmt blastdb`) against the sequence's Seq-ids, and `CWinMaskUtil::consider`.

use std::collections::HashSet;

use crate::ncbi_app::files::{input_lines, LineEnds};
use crate::ncbi_app::Diag;
use crate::seqid::accguide::parse_flags;
use crate::seqid::mask_handle::{mask_handle_key, GiimTree, TreeHandle};
use crate::seqid::parse::{seq_id_from_text, SeqIdError};
use crate::seqid::ParsedSeqId;

use super::work::{post_seq_id_warning, Record};

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

/// `CWinMaskConfig::FillIdList( file_name, id_list )`: `insert` of each id of the file `name`,
/// one per line, in order; a file that does not open gives none.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:471-484
/// The lines are those of `NcbiGetlineEOL` on the host (`ncbi_app::files::input_lines`: on Linux
/// a `\r` stays in the line, so a CRLF list matches nothing as text, as NCBI's Linux build), and
/// the empty ones are skipped. A line that starts with a space or a tab gives an empty id; one
/// whose first space comes at 0 after `>` (`> x`) too.
pub(super) fn fill_id_list<E>(
    name: &str,
    mut insert: impl FnMut(&[u8]) -> Result<(), E>,
) -> Result<(), E> {
    let Ok(bytes) = std::fs::read(name) else {
        return Ok(());
    };
    for line in input_lines(&bytes, LineEnds::HOST) {
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
        insert(id)?;
    }
    Ok(())
}

/// `CWinMaskUtil::CIdSet_SeqId`: the `CSeq_id_Handle`s of a list's ids (`idset`), as NCBI's id
/// mapper gives them (`seqid::mask_handle`).
#[derive(Debug, Default)]
pub(super) struct SeqIdSet {
    /// The handles that a pure function knows (`mask_handle_key`).
    keys: HashSet<Vec<u8>>,
    /// The `giim` tree of NCBI's mapper while the set holds its handles (`giim_handles`): the
    /// handle of a `giim` id depends on the handles made before it (`TreeHandle::Giim`).
    giim: GiimTree,
    giim_handles: HashSet<usize>,
    /// The handles that equal no other (a Textseq-id with neither accession nor name): they
    /// match nothing, and the set is not empty.
    fresh: usize,
}

impl SeqIdSet {
    /// `CIdSet_SeqId::insert( id_str )`: the handle of `CSeq_id( id_str )`; an id that it cannot
    /// read (`CSeqIdException`) is posted as an error and skipped. The parse's own warnings (its
    /// `ERR_POST_X` of error code 1306) are posted as they come.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_util.cpp:47-57
    /// NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_id.hpp:135-136
    /// The constructor parses with `fParse_AnyRaw`. The post carries the exception's `what()`
    /// (NCBI's report, a line feed at its end) in LOSAT's words (`crate::wording`).
    pub(super) fn insert(&mut self, id: &[u8], diag: &mut Diag) {
        let parsed = seq_id_from_text(id, parse_flags::ANY_RAW, &mut |warning| {
            post_seq_id_warning(diag, &warning)
        });
        match parsed {
            Ok(seq_id) => self.add(&seq_id),
            Err(error) => diag.err_post_bytes(&cannot_read(id, &error)),
        }
    }

    /// `idset.insert( CSeq_id_Handle::GetHandle( *id ) )`.
    fn add(&mut self, id: &ParsedSeqId) {
        match mask_handle_key(id) {
            Ok(Some(key)) => {
                self.keys.insert(key);
            }
            Ok(None) => self.fresh += 1,
            Err(TreeHandle::Giim { .. }) => {
                let ParsedSeqId::Giim(giim) = id else {
                    unreachable!("only a giim Seq-id has a tree handle")
                };
                let (handle, _) = self.giim.find_or_create(giim);
                self.giim_handles.insert(handle);
            }
        }
    }

    /// `empty()`.
    fn is_empty(&self) -> bool {
        self.keys.is_empty() && self.giim_handles.is_empty() && self.fresh == 0
    }

    /// `CIdSet_SeqId::find( bsh )` of a sequence whose Seq-ids are `ids`: whether one of their
    /// handles is in the set.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_util.cpp:60-73
    /// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp:244-257
    /// The handles are made when the sequence's scope adds its Bioseq, in the order of its
    /// Seq-ids, and dropped with that scope after the sequence: a `giim` id's handle is found
    /// in the tree of the list's handles and those made before it for the same sequence.
    fn find(&self, ids: &[ParsedSeqId]) -> bool {
        let mut giim = None;
        ids.iter().any(|id| match mask_handle_key(id) {
            Ok(Some(key)) => self.keys.contains(&key),
            Ok(None) => false,
            Err(TreeHandle::Giim { .. }) => {
                let ParsedSeqId::Giim(id) = id else {
                    unreachable!("only a giim Seq-id has a tree handle")
                };
                let tree = giim.get_or_insert_with(|| self.giim.clone());
                let (handle, _) = tree.find_or_create(id);
                self.giim_handles.contains(&handle)
            }
        })
    }
}

/// The error that `CIdSet_SeqId::insert` posts for the id `id` that `CSeq_id` cannot read
/// (`error`), in LOSAT's words: NCBI's message with the exception's `what()`, mapped as a whole.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_util.cpp:52-55
fn cannot_read(id: &[u8], error: &SeqIdError) -> Vec<u8> {
    let message = [
        b"CWinMaskConfig::FillIdList(): can't understand id: ".as_slice(),
        id,
        b": ",
        &error.ncbi_what_bytes(),
        b": ignoring",
    ]
    .concat();
    crate::wording::losat_bytes(&message).into_owned()
}

/// An id list of `CWinMaskConfig` (`CIdSet`): matched as text, or as Seq-ids.
#[derive(Debug)]
pub(super) enum IdSet {
    Text(TextMatchIds),
    SeqId(SeqIdSet),
}

/// An id list matched as text that asks for the title of a sequence of which LOSAT does not
/// make it (`Title::FromSource`).
pub(super) struct TitleFromSource;

impl IdSet {
    /// `empty()`.
    fn is_empty(&self) -> bool {
        match self {
            IdSet::Text(ids) => ids.is_empty(),
            IdSet::SeqId(ids) => ids.is_empty(),
        }
    }

    /// `find( bsh )` of `record`.
    fn find(&self, record: &Record) -> Result<bool, TitleFromSource> {
        match self {
            IdSet::Text(ids) => {
                match get_title(&record.ids, record.title.as_deref(), record.any_source) {
                    Title::Text(title) => Ok(ids.find(&title)),
                    Title::FromSource => Err(TitleFromSource),
                }
            }
            IdSet::SeqId(ids) => Ok(ids.find(&record.ids)),
        }
    }
}

/// The `-ids` and `-exclude_ids` lists of the run (`CWinMaskConfig::Ids`, `ExcludeIds`).
#[derive(Debug, Default)]
pub(super) struct IdLists {
    pub(super) ids: Option<IdSet>,
    pub(super) exclude_ids: Option<IdSet>,
}

impl IdLists {
    /// `CWinMaskUtil::consider( bsh, ids, exclude_ids )` of `record`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_util.cpp:218-245
    pub(super) fn consider(&self, record: &Record) -> Result<bool, TitleFromSource> {
        let mut result = true;
        if let Some(ids) = self.ids.as_ref().filter(|ids| !ids.is_empty()) {
            result = ids.find(record)?;
        }
        if let Some(exclude) = self.exclude_ids.as_ref().filter(|ids| !ids.is_empty()) {
            if exclude.find(record)? {
                result = false;
            }
        }
        Ok(result)
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

    fn seq_id(text: &str) -> ParsedSeqId {
        seq_id_from_text(text.as_bytes(), parse_flags::ANY_RAW, &mut |_| {}).unwrap()
    }

    fn seq_id_list(ids: &[&str]) -> SeqIdSet {
        let mut set = SeqIdSet::default();
        for id in ids {
            set.add(&seq_id(id));
        }
        set
    }

    fn record(ids: Vec<ParsedSeqId>) -> Record {
        Record {
            ids,
            title: Some(b"title".to_vec()),
            taxname: None,
            has_source: false,
            any_source: false,
            iupac: b"ACGT".to_vec(),
            texts: Vec::new(),
        }
    }

    // NCBI BLAST+ 2.17.0 windowmasker -text_match F on a database made with -parse_seqids (probe
    // s1 of line NP M5, 2026-11-21: tmf_wm_ps_ids_case, _unver, _forms, _lcl).
    #[test]
    fn seq_ids_match_as_ncbi_s_handles() {
        let set = seq_id_list(&["ab000001.1", "nc_000099.1", "lcl|REC2", "lcl|rec4."]);
        assert!(set.find(&[seq_id("gb|AB000001.1|")]));
        assert!(set.find(&[seq_id("ref|NC_000099.1|")]));
        assert!(set.find(&[seq_id("lcl|rec2")]));
        assert!(set.find(&[seq_id("lcl|rec4.")]));
        assert!(!set.find(&[seq_id("lcl|rec1")]));
        let bare = seq_id_list(&["BK000001.1", "P12345"]);
        assert!(bare.find(&[seq_id("tpg|BK000001.1|")]));
        let unversioned = seq_id_list(&["AB000001", "NC_000099", "AB000001.2"]);
        assert!(!unversioned.find(&[seq_id("gb|AB000001.1|")]));
        assert!(!unversioned.find(&[seq_id("ref|NC_000099.1|")]));
    }

    // tmf_wm_ids_ord of the same probe: `gnl|BL_ORD_ID|0` names OID 0 by a string tag "0"
    // (`CObject_id::SetStrOrId` takes 1-9 first), which is not the integer tag of the database.
    #[test]
    fn an_ordinal_id_of_zero_is_a_string_tag() {
        let set = seq_id_list(&["gnl|BL_ORD_ID|0", "gnl|BL_ORD_ID|3"]);
        let ord = |oid| {
            ParsedSeqId::General(crate::seqid::model::Dbtag {
                db: b"BL_ORD_ID".to_vec(),
                tag: crate::seqid::model::ObjectId::Id(oid),
            })
        };
        assert!(!set.find(&[ord(0)]));
        assert!(set.find(&[ord(3)]));
        assert!(set.find(&[ord(1), ord(3)]));
    }

    // The giim tree finds the first handle whose Giimport-id is not `Equals` (seq_id_tree.cpp
    // 2926-2941): an equal giim id of a sequence gets a new handle, another one of the same
    // number the list's.
    #[test]
    fn a_giim_id_finds_the_handle_of_another_giimport_id() {
        let set = seq_id_list(&["gim|5", "gim|5"]);
        assert!(!set.is_empty());
        assert!(!set.find(&[seq_id("gim|5")]));
        assert!(!set.find(&[seq_id("gim|6")]));
        let other = ParsedSeqId::Giim(crate::seqid::model::GiimportId {
            id: 5,
            db: Some(b"x".to_vec()),
            release: None,
        });
        assert!(set.find(&[other]));
    }

    // A list of which no id was read is empty: `consider` keeps every sequence.
    #[test]
    fn an_empty_seq_id_list_keeps_every_sequence() {
        let lists = IdLists {
            ids: Some(IdSet::SeqId(SeqIdSet::default())),
            exclude_ids: None,
        };
        assert!(matches!(
            lists.consider(&record(vec![seq_id("lcl|x")])),
            Ok(true)
        ));
        let lists = IdLists {
            ids: Some(IdSet::SeqId(seq_id_list(&["lcl|y"]))),
            exclude_ids: None,
        };
        assert!(matches!(
            lists.consider(&record(vec![seq_id("lcl|x")])),
            Ok(false)
        ));
        assert!(matches!(
            lists.consider(&record(vec![seq_id("lcl|Y")])),
            Ok(true)
        ));
        let lists = IdLists {
            ids: None,
            exclude_ids: Some(IdSet::SeqId(seq_id_list(&["lcl|y"]))),
        };
        assert!(matches!(
            lists.consider(&record(vec![seq_id("lcl|y")])),
            Ok(false)
        ));
    }

    // tmf_wm_ids_bad of probe s1: NCBI's post (its build tree aside), here in LOSAT's words; the
    // id keeps its bytes.
    #[test]
    fn an_id_that_cseq_id_cannot_read_is_posted_with_ncbi_s_report() {
        let parse = |id: &[u8]| seq_id_from_text(id, parse_flags::ANY_RAW, &mut |_| {});
        let error = parse(b"gb|").unwrap_err();
        assert_eq!(
            String::from_utf8(error.ncbi_what_bytes()).unwrap(),
            format!(
                "NCBI C++ Exception:\n    T0 \"{}/src/objects/seq/../seqloc/Seq_id.cpp\", line \
                 2501: Error: (CSeqIdException::eFormat) ncbi::objects::CSeq_id::Set() - \
                 Malformatted ID gb|\n",
                crate::ncbi_exception::NCBI_BUILD_ROOT
            )
        );
        assert_eq!(
            cannot_read(b"gb|", &error),
            b"CWinMaskConfig::FillIdList(): can't understand id: gb|: Malformatted ID gb|\n: ignoring"
        );
        let error = parse(b"caf\xe9\r").unwrap_err();
        assert_eq!(
            cannot_read(b"caf\xe9\r", &error),
            b"CWinMaskConfig::FillIdList(): can't understand id: caf\xe9\r: Malformatted ID caf\xe9\n: ignoring"
        );
        let error = parse(b"").unwrap_err();
        assert_eq!(
            cannot_read(b"", &error),
            b"CWinMaskConfig::FillIdList(): can't understand id: : Empty bare accession supplied\n: ignoring"
        );
    }
}
