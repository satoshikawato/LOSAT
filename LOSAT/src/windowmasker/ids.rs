//! `-ids` and `-exclude_ids`: the id lists of `CWinMaskConfig::FillIdList`, matched as text
//! (`CWinMaskUtil::CIdSet_TextMatch`, the default) against the first word of each sequence's
//! title (`sequence::GetTitle`), or as Seq-ids (`CWinMaskUtil::CIdSet_SeqId`, `-text_match false`
//! with `-infmt blastdb`) against the sequence's Seq-ids, and `CWinMaskUtil::consider`.

use std::collections::{HashMap, HashSet};

use crate::ncbi_app::files::{input_lines, LineEnds};
use crate::ncbi_app::Diag;
use crate::seqid::accguide::parse_flags;
use crate::seqid::mask_handle::{
    handle_seq_id, mask_handle_key, mask_info_key, GiimTree, TreeHandle,
};
use crate::seqid::parse::{seq_id_from_text, SeqIdError};
use crate::seqid::ParsedSeqId;

use crate::seqmasks_io::writer::access_seq_id;

use super::work::{post_seq_id_warning, Record};

/// `CWinMaskUtil::CIdSet_TextMatch`: the ids of a list, by their number of `|`-words
/// (`nword_sets_`).
#[derive(Debug, Default, PartialEq, Eq)]
pub(super) struct TextMatchIds {
    nword_sets: Vec<HashSet<Vec<u8>>>,
}

/// An id of a list that NCBI cannot store, one of no words: `CIdSet_TextMatch::insert` posts
/// `bad id` and then stores the id at its `Uint4` count of words less 1, 4294967295, far past the
/// end of its sets (undefined behaviour in C++; NCBI BLAST+ 2.17.0 dies with SIGSEGV). LOSAT posts the same error and skips the id, as the error says
/// (`PD-LOSAT-NCBI-DEFECTS` approved exception 22).
#[derive(Debug, PartialEq, Eq)]
pub(super) struct EmptyId {
    pub(super) id: Vec<u8>,
}

impl EmptyId {
    /// The error that `CIdSet_TextMatch::insert` posts for the id, in LOSAT's words (NCBI's
    /// message through `crate::wording`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_util.cpp:107-111
    pub(super) fn post(&self) -> Vec<u8> {
        let message = [
            b"CWinMaskConfig::CIdSet_TextMatch::insert(): bad id: ".as_slice(),
            &self.id,
            b": ignoring",
        ]
        .concat();
        crate::wording::losat_bytes(&message).into_owned()
    }
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
    /// `CIdSet_TextMatch::insert` of `id`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_util.cpp:103-123
    /// An id of no words (empty, `>`, `|`, `>|`) goes on to `nword_sets_[nwords - 1]` with `nwords`
    /// 0, past the end of the sets: `EmptyId`, the set
    /// unchanged, as NCBI's of the list without that line (`resize` to 0 words changes nothing).
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

    /// `CIdSet_TextMatch::find` of an id: whether a run of consecutive words of `id` is an id
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

    /// `CIdSet_TextMatch::find` of a sequence whose `GetTitle` is `title`.
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

/// `CWinMaskConfig::FillIdList`: `insert` of each id of the file `name`,
/// one per line, in order; a file that does not open gives none.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:467-482
/// The lines are those of `NcbiGetlineEOL` on the host (`ncbi_app::files::input_lines`: on Linux
/// a `\r` stays in the line, so a CRLF list matches nothing as text, as NCBI's Linux build), and
/// the empty ones are skipped. A line that starts with a space or a tab gives an empty id; one
/// whose first space comes at 0 after `>` (`> x`) too.
pub(super) fn fill_id_list(name: &str, mut insert: impl FnMut(&[u8])) {
    let Ok(bytes) = std::fs::read(name) else {
        return;
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
        insert(id);
    }
}

/// `CWinMaskUtil::CIdSet_SeqId`: the `CSeq_id_Handle`s of a list's ids (`idset`), as NCBI's id
/// mapper gives them (`seqid::mask_handle`). They live for the whole run (the set is a member of
/// `CWinMaskConfig`), so the info of each is the one its list id made, also for the handles of
/// the sequences that share it.
#[derive(Debug, Default)]
pub(super) struct SeqIdSet {
    /// The handles that a pure function knows (`mask_handle_key`).
    keys: HashSet<Vec<u8>>,
    /// The list id that made each info of those handles (`mask_info_key`), the first to reach it.
    made_from: HashMap<Vec<u8>, ParsedSeqId>,
    /// The `giim` tree of NCBI's mapper while the set holds its handles (`giim_handles`): the
    /// handle of a `giim` id depends on the handles made before it (`TreeHandle::Giim`).
    giim: GiimTree,
    giim_handles: HashSet<usize>,
    /// The handles that equal no other (a Textseq-id with neither accession nor name): they
    /// match nothing, and the set is not empty. No list id reaches it today: NCBI's parse throws
    /// for such a Textseq-id (Textseq_id.cpp:128-132), as LOSAT's does.
    fresh: usize,
}

impl SeqIdSet {
    /// `CIdSet_SeqId::insert`: the handle of the Seq-id that `CSeq_id`'s constructor reads from
    /// `id`; an id that it cannot read (`CSeqIdException`) is posted as an error and skipped. The
    /// parse's own warnings (error code 1306) are posted as they come.
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

    /// The handle of `id` (`CSeq_id_Handle::GetHandle`), added to the set.
    fn add(&mut self, id: &ParsedSeqId) {
        match mask_handle_key(id) {
            Ok(Some(key)) => {
                self.keys.insert(key);
                if let Ok(Some(info)) = mask_info_key(id) {
                    self.made_from.entry(info).or_insert_with(|| id.clone());
                }
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

    /// `CIdSet_SeqId::find` of a sequence whose Seq-ids are `ids`: whether one of their handles is
    /// in the set.
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

    /// `GetSeqId()` of the access handle of a sequence whose Seq-ids are `ids`
    /// (`seqmasks_io::writer::access_seq_id`) where a list id made its info; `None` where the
    /// sequence's own Seq-id did. The access Seq-id is the sequence's first gi, Textseq-id or
    /// other Seq-id, so no Seq-id of the sequence before it can share its info: only the list's.
    ///
    /// NCBI reference (598d8ae6): c++/include/objects/seq/seq_id_handle.hpp:300-312
    fn handle_seq_id(&self, ids: &[ParsedSeqId]) -> Option<ParsedSeqId> {
        let id = access_seq_id(ids);
        match mask_info_key(id) {
            Ok(Some(info)) => self
                .made_from
                .get(&info)
                .map(|first| handle_seq_id(first, id)),
            Ok(None) => None,
            Err(TreeHandle::Giim { id: number }) => {
                let ParsedSeqId::Giim(value) = id else {
                    unreachable!("only a giim Seq-id has a tree handle")
                };
                let handle = self.giim.find(value)?;
                Some(ParsedSeqId::Giim(
                    self.giim.made_from(number, handle).clone(),
                ))
            }
        }
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

impl IdSet {
    /// `empty()`.
    fn is_empty(&self) -> bool {
        match self {
            IdSet::Text(ids) => ids.is_empty(),
            IdSet::SeqId(ids) => ids.is_empty(),
        }
    }

    /// `find` of `record`: a list matched as text against its title (`Record::get_title`).
    fn find(&self, record: &Record) -> bool {
        match self {
            IdSet::Text(ids) => ids.find(
                record
                    .get_title
                    .as_deref()
                    .expect("the sources give the titles of a run with a list matched as text"),
            ),
            IdSet::SeqId(ids) => ids.find(&record.ids),
        }
    }

    /// Whether the list is matched as text, against the titles (`CIdSet_TextMatch`).
    fn matches_titles(&self) -> bool {
        matches!(self, IdSet::Text(_))
    }
}

/// The `-ids` and `-exclude_ids` lists of the run (`CWinMaskConfig::Ids`, `ExcludeIds`).
#[derive(Debug, Default)]
pub(super) struct IdLists {
    pub(super) ids: Option<IdSet>,
    pub(super) exclude_ids: Option<IdSet>,
}

impl IdLists {
    /// `CWinMaskUtil::consider` of `record` and the two lists.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_util.cpp:218-245
    pub(super) fn consider(&self, record: &Record) -> bool {
        let mut result = true;
        if let Some(ids) = self.ids.as_ref().filter(|ids| !ids.is_empty()) {
            result = ids.find(record);
        }
        if let Some(exclude) = self.exclude_ids.as_ref().filter(|ids| !ids.is_empty()) {
            if exclude.find(record) {
                result = false;
            }
        }
        result
    }

    /// Whether `consider` reads the titles of the sequences (`Record::get_title`): a list that is
    /// not empty is matched as text.
    pub(super) fn match_titles(&self) -> bool {
        [&self.ids, &self.exclude_ids]
            .into_iter()
            .flatten()
            .any(|list| list.matches_titles() && !list.is_empty())
    }

    /// `GetSeqId()` of the access handle of a sequence whose Seq-ids are `ids` where the handle of
    /// an id of a list matched as Seq-ids made its info (`SeqIdSet::handle_seq_id`), as the
    /// `seqloc_*` and `maskinfo_*` writers print it; `None` otherwise.
    pub(super) fn handle_seq_id(&self, ids: &[ParsedSeqId]) -> Option<ParsedSeqId> {
        [&self.ids, &self.exclude_ids]
            .into_iter()
            .flatten()
            .find_map(|list| match list {
                IdSet::SeqId(set) => set.handle_seq_id(ids),
                IdSet::Text(_) => None,
            })
    }
}

#[cfg(test)]
#[path = "ids_tests.rs"]
mod tests;
