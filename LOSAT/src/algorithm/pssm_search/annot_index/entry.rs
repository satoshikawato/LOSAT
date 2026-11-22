//! The PSSM's query Seq-entry as NCBI's object manager adds it to the scope (`AddTopLevelSeqEntry`):
//! a Bioseq, or a Bioseq-set with its members at any depth. The TSE is made with every object's
//! information (phase C: a Bioseq-set makes its members, then its own annotations), then attached
//! (each Bioseq-set's integer id and each Bioseq's Seq-ids, which NCBI keeps once per TSE), then
//! its annotations are indexed (phase I: a Bioseq-set's own annotations, then its members').
//!
//! NCBI reference (598d8ae6): c++/src/objmgr/bioseq_set_info.cpp:246-264 (phase C of a set)
//! NCBI reference (598d8ae6): c++/src/objmgr/tse_info.cpp:150-167 (the TSE is made, then attached)
//! NCBI reference (598d8ae6): c++/src/objmgr/bioseq_set_info.cpp:139-154 (attach of a set)
//! NCBI reference (598d8ae6): c++/src/objmgr/bioseq_set_info.cpp:564-570 (phase I of a set)

use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};

use super::range::{Handle, LiveHandles};
use super::view::V;
use super::{Refusal, Stop};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::serial::asn_types::table::named;
use crate::report::serial::pssm::SeqEntry;
use crate::seqid::ParsedSeqId;

/// The TSE of a query Seq-entry: its annotations in phase I's order, phase C's order of what it
/// makes, and what attaches, in NCBI's order.
pub(super) struct Tse<'a> {
    pub(super) annots: Vec<V<'a>>,
    /// What phase C makes, in its order.
    pub(super) construction: Vec<Made>,
    attach: Vec<Attach<'a>>,
}

/// What phase C makes: a Bioseq's Seq-id handles (`CBioseq_Info::m_Id`, bioseq_info.cpp:300-307),
/// at its place in `Tse::attach`, or an annotation's information, at its place in `Tse::annots`.
#[derive(Clone, Copy, Debug)]
pub(super) enum Made {
    Bioseq(usize),
    Annot(usize),
}

/// One object that attaches to the TSE.
enum Attach<'a> {
    /// A Bioseq-set's `id`, when it is an integer.
    SetId(i64),
    /// A Bioseq's Seq-ids.
    Bioseq(Vec<&'a ParsedSeqId>),
}

/// A Bioseq-set being walked (`Tse::add_set`): its members and the next one, and its own
/// annotations (phase C makes them after the members).
struct SetWalk<'a> {
    members: Option<V<'a>>,
    next: usize,
    own: Vec<Made>,
}

impl<'a> Tse<'a> {
    /// The TSE of `entry`.
    pub(super) fn of(entry: &'a SeqEntry) -> Self {
        let mut tse = Tse {
            annots: Vec::new(),
            construction: Vec::new(),
            attach: Vec::new(),
        };
        match entry {
            SeqEntry::Seq(query) => {
                tse.construction.push(Made::Bioseq(tse.attach.len()));
                tse.attach.push(Attach::Bioseq(query.ids.iter().collect()));
                for annot in query.annot.iter().flatten() {
                    tse.construction.push(Made::Annot(tse.annots.len()));
                    tse.annots.push(V::new(named::SEQ_ANNOT, annot));
                }
            }
            SeqEntry::Set(set) => {
                let mut construction = Vec::new();
                tse.add_set(V::new(named::BIOSEQ_SET, set), &mut construction);
                tse.construction = construction;
            }
        }
        tse
    }

    /// Adds the set `set` (phase I's order of `annots` and the attach order are the order of the
    /// walk: the set's own annotations and id, then its members); its phase C order (members, each
    /// Bioseq's Seq-ids before its annotations, then the set's own annotations) goes to
    /// `construction`. The sets inside it are walked from a list on the heap (`SetWalk`), to any
    /// depth the reader read.
    fn add_set(&mut self, set: V<'a>, construction: &mut Vec<Made>) {
        let mut walks = vec![self.enter_set(set)];
        while let Some(walk) = walks.last_mut() {
            let Some(member) = walk.members.and_then(|members| members.item(walk.next)) else {
                let done = walks.pop().expect("a set");
                construction.extend(done.own);
                continue;
            };
            walk.next += 1;
            match member.choice() {
                ("seq", bioseq) => {
                    construction.push(Made::Bioseq(self.attach.len()));
                    self.attach.push(Attach::Bioseq(
                        bioseq.field("id").items().map(V::seq_id).collect(),
                    ));
                    for annot in bioseq.get("annot").into_iter().flat_map(V::items) {
                        construction.push(Made::Annot(self.annots.len()));
                        self.annots.push(annot);
                    }
                }
                (_, set) => {
                    let walk = self.enter_set(set);
                    walks.push(walk);
                }
            }
        }
    }

    /// The set `set` entered: its id and its own annotations added, its members to walk.
    fn enter_set(&mut self, set: V<'a>) -> SetWalk<'a> {
        if let Some(id) = set.get("id") {
            if let ("id", number) = id.choice() {
                self.attach.push(Attach::SetId(number.int()));
            }
        }
        let own: Vec<Made> = set
            .get("annot")
            .into_iter()
            .flat_map(V::items)
            .map(|annot| {
                self.annots.push(annot);
                Made::Annot(self.annots.len() - 1)
            })
            .collect();
        SetWalk {
            members: set.get("seq-set"),
            next: 0,
            own,
        }
    }

    /// The number of objects that attach (Bioseq-set ids and Bioseqs).
    pub(super) fn objects(&self) -> usize {
        self.attach.len()
    }

    /// The Seq-ids of the Bioseq at `index` of the attach (`Made::Bioseq`).
    pub(super) fn bioseq_ids(&self, index: usize) -> &[&'a ParsedSeqId] {
        match &self.attach[index] {
            Attach::Bioseq(ids) => ids,
            Attach::SetId(_) => unreachable!("phase C makes a Bioseq's handles"),
        }
    }

    /// The attach of the TSE's objects (`CTSE_Info::x_SetBioseq_setId`, `x_SetBioseqIds`): NCBI's
    /// exception of a Bioseq-set id or a Bioseq's Seq-id handle that the TSE has already, in one
    /// pass over a map of the handles. `ids` holds each Bioseq's handles (phase C, `Made::Bioseq`),
    /// `handles` the infos they were made from. `ids_known`: no other sequence is in NCBI's scope (a
    /// search of `-db`), so the first Seq-id of the entry with a handle made its info; else (the
    /// subjects of `-subject` came first) LOSAT refuses a message in which a subject's Seq-id may
    /// show through (`seqid::mask_handle::handle_fasta_is_own`).
    pub(super) fn attach(
        &self,
        ids: &[Vec<Handle>],
        handles: &LiveHandles,
        ids_known: bool,
    ) -> Result<(), Stop> {
        let mut set_ids = HashSet::new();
        let mut held: HashMap<&Handle, usize> = HashMap::new();
        for (index, object) in self.attach.iter().enumerate() {
            match object {
                Attach::SetId(id) => {
                    if *id >= 0 && !set_ids.insert(*id) {
                        return Err(Stop::Exception(NcbiException::new(
                            SET_ID_DUPLICATE,
                            format!(" duplicate Bioseq_set id '{id}"),
                        )));
                    }
                }
                Attach::Bioseq(seq_ids) => {
                    for (&id, handle) in seq_ids.iter().zip(&ids[index]) {
                        let first = match held.entry(handle) {
                            Entry::Vacant(slot) => {
                                slot.insert(index);
                                continue;
                            }
                            Entry::Occupied(slot) => *slot.get(),
                        };
                        if !ids_known
                            && [first, index].iter().any(|&object| {
                                self.bioseq_ids(object)
                                    .iter()
                                    .any(|id| !crate::seqid::mask_handle::handle_fasta_is_own(id))
                            })
                        {
                            return Err(Stop::Refuse(Refusal::DuplicateIdInScope));
                        }
                        return Err(Stop::Exception(NcbiException::new(
                            BIOSEQ_ID_DUPLICATE,
                            format!(
                                "duplicate Bioseq id {} present in\n  seq1: {}\n  seq2: {}",
                                handles.as_string(handle, id),
                                self.id_string(ids, handles, first),
                                self.id_string(ids, handles, index)
                            ),
                        )));
                    }
                }
            }
        }
        Ok(())
    }

    /// `CBioseq_Info::IdString()` of the Bioseq attached at `index`: its handles' `AsString()`.
    // NCBI reference (598d8ae6): c++/src/objmgr/bioseq_info.cpp:1221-1230
    fn id_string(&self, ids: &[Vec<Handle>], handles: &LiveHandles, index: usize) -> String {
        self.bioseq_ids(index)
            .iter()
            .zip(&ids[index])
            .map(|(id, handle)| handles.as_string(handle, id))
            .collect::<Vec<_>>()
            .join(" | ")
    }
}

// NCBI reference (598d8ae6): c++/src/objmgr/tse_info.cpp:891-905
const BIOSEQ_ID_DUPLICATE: ThrowSite = ThrowSite {
    file: "src/objmgr/tse_info.cpp",
    line: 901,
    kind: "CObjMgrException::eAddDataError",
    function: Some("ncbi::objects::CTSE_Info::x_SetBioseqIds"),
};

// NCBI reference (598d8ae6): c++/src/objmgr/tse_info.cpp:941-953
const SET_ID_DUPLICATE: ThrowSite = ThrowSite {
    file: "src/objmgr/tse_info.cpp",
    line: 951,
    kind: "CObjMgrException::eAddDataError",
    function: Some("ncbi::objects::CTSE_Info::x_SetBioseq_setId"),
};
