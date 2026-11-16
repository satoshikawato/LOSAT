//! The PSSM's query Seq-entry as NCBI's object manager adds it to the scope (`AddTopLevelSeqEntry`):
//! a Bioseq, or a Bioseq-set with its members at any depth. The TSE is made with every object's
//! information (phase C: a Bioseq-set makes its members, then its own annotations), then attached
//! (each Bioseq-set's integer id and each Bioseq's Seq-ids, which NCBI keeps once per TSE), then
//! its annotations are indexed (phase I: a Bioseq-set's own annotations, then its members').
//!
//! NCBI reference (598d8ae6): c++/src/objmgr/bioseq_set_info.cpp:246-264 (phase C of a set)
//! ```c++
//!     if ( obj.IsSetSeq_set() ) {
//!         NON_CONST_ITERATE ( TObject::TSeq_set, it, obj.SetSeq_set() ) {
//!             CRef<CSeq_entry_Info> info(new CSeq_entry_Info(**it));
//!             m_Seq_set.push_back(info);
//!             x_AttachEntry(info);
//!         }
//!     }
//!     if ( obj.IsSetAnnot() ) {
//!         x_SetAnnot();
//!     }
//! ```
//! NCBI reference (598d8ae6): c++/src/objmgr/tse_info.cpp:150-167 (the TSE is made, then attached)
//! ```c++
//!     SetSeq_entry(entry);
//!     m_LoadState = eLoaded;
//!
//!     x_TSEAttach(*this);
//! ```
//! NCBI reference (598d8ae6): c++/src/objmgr/bioseq_set_info.cpp:139-154 (attach of a set)
//! ```c++
//!     if ( IsSetId() ) {
//!         m_Bioseq_set_Id = x_GetBioseq_set_Id(GetId());
//!         if ( m_Bioseq_set_Id >= 0 ) {
//!             tse.x_SetBioseq_setId(m_Bioseq_set_Id, this);
//!         }
//!     }
//!     SetBioObjectId(tse.x_IndexBioseq_set(this));
//!     // members
//!     NON_CONST_ITERATE ( TSeq_set, it, m_Seq_set ) {
//!         (*it)->x_TSEAttach(tse);
//!     }
//! ```
//! NCBI reference (598d8ae6): c++/src/objmgr/bioseq_set_info.cpp:564-570 (phase I of a set)
//! ```c++
//! void CBioseq_set_Info::x_UpdateAnnotIndexContents(CTSE_Info& tse)
//! {
//!     TParent::x_UpdateAnnotIndexContents(tse);
//!     for ( size_t i = 0; i < m_Seq_set.size(); ++i ) {
//!         m_Seq_set[i]->x_UpdateAnnotIndex(tse);
//!     }
//! }
//! ```

use std::collections::HashSet;

use super::view::V;
use super::{Refusal, Stop};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::serial::asn_types::table::named;
use crate::report::serial::pssm::SeqEntry;
use crate::seqid::ParsedSeqId;

/// The TSE of a query Seq-entry: its annotations in phase I's order, the order of phase C over them,
/// and what attaches, in NCBI's order.
pub(super) struct Tse<'a> {
    pub(super) annots: Vec<V<'a>>,
    /// The positions in `annots` in phase C's order.
    pub(super) construction: Vec<usize>,
    attach: Vec<Attach<'a>>,
}

/// One object that attaches to the TSE.
enum Attach<'a> {
    /// A Bioseq-set's `id`, when it is an integer.
    SetId(i64),
    /// A Bioseq's Seq-ids.
    Bioseq(Vec<&'a ParsedSeqId>),
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
                tse.attach.push(Attach::Bioseq(query.ids.iter().collect()));
                for annot in query.annot.iter().flatten() {
                    tse.construction.push(tse.annots.len());
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
    /// walk: the set's own annotations and id, then its members); its phase C order (members, then
    /// its own annotations) goes to `construction`.
    fn add_set(&mut self, set: V<'a>, construction: &mut Vec<usize>) {
        if let Some(id) = set.get("id") {
            if let ("id", number) = id.choice() {
                self.attach.push(Attach::SetId(number.int()));
            }
        }
        let own: Vec<usize> = set
            .get("annot")
            .into_iter()
            .flat_map(V::items)
            .map(|annot| {
                self.annots.push(annot);
                self.annots.len() - 1
            })
            .collect();
        for member in set.get("seq-set").into_iter().flat_map(V::items) {
            match member.choice() {
                ("seq", bioseq) => {
                    self.attach.push(Attach::Bioseq(
                        bioseq.field("id").items().map(V::seq_id).collect(),
                    ));
                    for annot in bioseq.get("annot").into_iter().flat_map(V::items) {
                        construction.push(self.annots.len());
                        self.annots.push(annot);
                    }
                }
                (_, set) => self.add_set(set, construction),
            }
        }
        construction.extend(own);
    }

    /// The Seq-ids of the TSE's Bioseqs, in order.
    pub(super) fn ids(&self) -> Vec<ParsedSeqId> {
        self.attach
            .iter()
            .flat_map(|object| match object {
                Attach::Bioseq(ids) => ids.as_slice(),
                Attach::SetId(_) => &[],
            })
            .map(|&id| id.clone())
            .collect()
    }

    /// The attach of the TSE's objects: NCBI's exception of a Bioseq-set id or a Bioseq id that the
    /// TSE has already (`CTSE_Info::x_SetBioseq_setId`, `x_SetBioseqIds`), when LOSAT can make its
    /// text (`ids_known`: no other Seq-id is in NCBI's scope, as with `-db`; LOSAT tells equal
    /// Seq-ids, and Seq-ids that NCBI's scope never takes as one); else LOSAT's refusal.
    pub(super) fn attach(&self, ids_known: bool) -> Result<(), Stop> {
        let ids: Vec<&ParsedSeqId> = self
            .attach
            .iter()
            .flat_map(|object| match object {
                Attach::Bioseq(ids) => ids.as_slice(),
                Attach::SetId(_) => &[],
            })
            .copied()
            .collect();
        for (n, a) in ids.iter().enumerate() {
            if ids[..n].iter().any(|b| a != b && may_be_one(a, b)) {
                return Err(Stop::Refuse(Refusal::IdsMaybeOne));
            }
        }
        let mut set_ids = HashSet::new();
        let mut held: Vec<(&ParsedSeqId, usize)> = Vec::new();
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
                Attach::Bioseq(ids) => {
                    for &id in ids {
                        if let Some(&(_, first)) = held.iter().find(|(other, _)| *other == id) {
                            if !ids_known {
                                return Err(Stop::Refuse(Refusal::DuplicateIdInScope));
                            }
                            return Err(Stop::Exception(NcbiException::new(
                                BIOSEQ_ID_DUPLICATE,
                                format!(
                                    "duplicate Bioseq id {} present in\n  seq1: {}\n  seq2: {}",
                                    fasta(id),
                                    self.id_string(first),
                                    self.id_string(index)
                                ),
                            )));
                        }
                        held.push((id, index));
                    }
                }
            }
        }
        Ok(())
    }

    /// `CBioseq_Info::IdString()` of the Bioseq attached at `index`.
    fn id_string(&self, index: usize) -> String {
        match &self.attach[index] {
            Attach::Bioseq(ids) => ids
                .iter()
                .map(|id| fasta(id))
                .collect::<Vec<_>>()
                .join(" | "),
            Attach::SetId(_) => unreachable!("a Bioseq holds Seq-ids"),
        }
    }
}

/// `CSeq_id_Handle::AsString()` of the handle of `id`, whose Seq-id is `id` (the first of the TSE
/// with it: `attach` refuses Seq-ids that NCBI may take as one and that differ).
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_handle.cpp:220-233
// ```c++
//     if ( IsGi() ) {
//         os << "gi|" << m_Packed;
//     }
//     else if ( m_Info ) {
//         GetSeqId()->WriteAsFasta(os);
// ```
fn fasta(id: &ParsedSeqId) -> String {
    let mut out = Vec::new();
    crate::seqid::label::fasta(id, &mut out);
    String::from_utf8_lossy(&out).into_owned()
}

/// Whether NCBI's scope may take `a` and `b` as one sequence (the keys of `seqid::id_scope`).
fn may_be_one(a: &ParsedSeqId, b: &ParsedSeqId) -> bool {
    let keys = |id: &ParsedSeqId| {
        let mut own = Vec::new();
        crate::seqid::defline_set::scope_keys(id, &mut own);
        let mut other = Vec::new();
        for pair in [
            crate::seqid::id_scope::version_keys(id),
            crate::seqid::id_scope::name_keys(id),
        ]
        .into_iter()
        .flatten()
        {
            own.push(pair.0);
            other.push(pair.1);
        }
        (own, other)
    };
    let (own_a, other_a) = keys(a);
    let (own_b, other_b) = keys(b);
    own_a
        .iter()
        .any(|key| own_b.contains(key) || other_b.contains(key))
        || other_a.iter().any(|key| own_b.contains(key))
}

// NCBI reference (598d8ae6): c++/src/objmgr/tse_info.cpp:891-905
// ```c++
//         ITERATE ( CBioseq_Info::TId, it, info->GetId() ) {
//             pair<TBioseqs::iterator, bool> ins =
//                 m_Bioseqs.insert(TBioseqs::value_type(*it, info));
//             if ( !ins.second ) {
//                 // No duplicate bioseqs in the same TSE
//                 NCBI_THROW(CObjMgrException, eAddDataError,
//                            "duplicate Bioseq id "+it->AsString()+" present in"+
//                            "\n  seq1: " + ins.first->second->IdString()+
//                            "\n  seq2: " + info->IdString());
// ```
const BIOSEQ_ID_DUPLICATE: ThrowSite = ThrowSite {
    file: "src/objmgr/tse_info.cpp",
    line: 901,
    kind: "CObjMgrException::eAddDataError",
    function: Some("ncbi::objects::CTSE_Info::x_SetBioseqIds"),
};

// NCBI reference (598d8ae6): c++/src/objmgr/tse_info.cpp:941-953
// ```c++
//     pair<TBioseq_sets::iterator, bool> ins =
//         m_Bioseq_sets.insert(TBioseq_sets::value_type(key, info));
//     ...
//         NCBI_THROW(CObjMgrException, eAddDataError,
//                    " duplicate Bioseq_set id '"+NStr::IntToString(key));
// ```
const SET_ID_DUPLICATE: ThrowSite = ThrowSite {
    file: "src/objmgr/tse_info.cpp",
    line: 951,
    kind: "CObjMgrException::eAddDataError",
    function: Some("ncbi::objects::CTSE_Info::x_SetBioseq_setId"),
};
