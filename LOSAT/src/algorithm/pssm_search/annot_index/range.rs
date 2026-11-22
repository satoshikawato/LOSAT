//! The ranges of NCBI's object manager: `COpenRange<TSeqPos>` (`Range`), a Seq-id's ranges of a
//! location (`CHandleRange`, `HandleRange`) and a location's ranges by Seq-id
//! (`CHandleRangeMap`, `HandleRangeMap`), which the annotation index asks for each Seq-id's
//! overlapping range.

use super::view::V;
use std::collections::HashMap;

use crate::seqid::mask_handle::{
    handle_seq_id, mask_handle_key, mask_handle_order, mask_handle_packed, mask_info_key, GiimTree,
    TreeHandle,
};
use crate::seqid::model::GiimportId;
use crate::seqid::ParsedSeqId;

/// `GetPositionMax()` of `TSeqPos`.
pub(super) const MAX: u32 = u32::MAX;

// NCBI reference (598d8ae6): c++/include/util/range.hpp:148-151, 255-308, 345-357, 398-401
/// `COpenRange<TSeqPos>`: `from` and the open end `to_open`, empty when `to_open <= from`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Range {
    pub(super) from: u32,
    pub(super) to_open: u32,
}

impl Range {
    /// `GetEmpty()`.
    pub(super) const EMPTY: Range = Range {
        from: MAX,
        to_open: MAX,
    };
    /// `GetWhole()`.
    pub(super) const WHOLE: Range = Range {
        from: 0,
        to_open: MAX,
    };

    /// `CRange<TSeqPos>` of `from` and `to`: the open end is `to` plus 1, wrapping (`to` 4294967295
    /// makes an empty range).
    pub(super) fn new(from: u32, to: u32) -> Range {
        Range {
            from,
            to_open: to.wrapping_add(1),
        }
    }

    pub(super) fn is_empty(self) -> bool {
        self.to_open <= self.from
    }

    /// `operator+=` (`CombineWith`).
    pub(super) fn combine(&mut self, other: Range) {
        if !other.is_empty() {
            if !self.is_empty() {
                self.from = self.from.min(other.from);
                self.to_open = self.to_open.max(other.to_open);
            } else {
                *self = other;
            }
        }
    }
}

/// `ENa_strand` (`Na-strand`): unknown 0, plus 1, minus 2, both 3, both-rev 4, other 255; a
/// Seq-table column's value of any `int`.
pub(super) type Strand = i64;
/// `eNa_strand_unknown`.
pub(super) const UNKNOWN: Strand = 0;
const MINUS: Strand = 2;
const BOTH: Strand = 3;
const BOTH_REV: Strand = 4;

// NCBI reference (598d8ae6): c++/include/objmgr/impl/handle_range.hpp:185-199
fn includes_plus(strand: Strand) -> bool {
    strand != MINUS
}

fn includes_minus(strand: Strand) -> bool {
    matches!(strand, UNKNOWN | MINUS | BOTH | BOTH_REV)
}

// NCBI reference (598d8ae6): c++/include/objects/seqloc/Na_strand.hpp:75-79
fn is_reverse(strand: Strand) -> bool {
    matches!(strand, MINUS | BOTH_REV)
}

/// `CHandleRange`: one Seq-id's ranges of a location, with their totals by strand.
#[derive(Clone, Debug)]
pub(super) struct HandleRange {
    ranges: Vec<(Range, Strand)>,
    plus: Range,
    minus: Range,
    pub(super) circular: bool,
    single_strand: bool,
    more_before: bool,
    more_after: bool,
}

impl Default for HandleRange {
    fn default() -> Self {
        HandleRange {
            ranges: Vec::new(),
            plus: Range::EMPTY,
            minus: Range::EMPTY,
            circular: false,
            single_strand: true,
            more_before: false,
            more_after: false,
        }
    }
}

impl HandleRange {
    // NCBI reference (598d8ae6): c++/src/objmgr/handle_range.cpp:115-206
    // The map's `AddRange` (handle_range_map.cpp:277) passes `more_before` and `more_after`
    // false: it sets the flags itself (`HandleRangeMap::add`).
    /// `AddRange` of `range`, `strand` and `circular_rna`, with `more_before` and `more_after`
    /// false.
    pub(super) fn add(&mut self, range: Range, strand: Strand, circular_rna: bool) {
        if let Some(&(_, first)) = self.ranges.first() {
            if self.single_strand {
                if strand != first {
                    if self.circular {
                        let mut total = self.plus;
                        total.combine(self.minus);
                        self.plus = if includes_plus(first) {
                            total
                        } else {
                            Range::EMPTY
                        };
                        self.minus = if includes_minus(first) {
                            total
                        } else {
                            Range::EMPTY
                        };
                        self.circular = false;
                    }
                    self.single_strand = false;
                } else if !self.circular && !circular_rna {
                    if let Some(&(last, _)) = self
                        .ranges
                        .iter()
                        .rev()
                        .find(|(range, _)| !range.is_empty())
                    {
                        self.circular = if includes_plus(strand) {
                            range.from < last.from
                        } else {
                            range.from > last.from
                        };
                    }
                    if self.circular {
                        let minus = self.minus;
                        self.plus.combine(minus);
                        self.minus = Range::EMPTY;
                    }
                }
            }
        }
        self.ranges.push((range, strand));
        if !self.circular {
            if includes_plus(strand) {
                self.plus.combine(range);
            }
            if includes_minus(strand) {
                self.minus.combine(range);
            }
        } else {
            self.minus.combine(range);
        }
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/handle_range.cpp:331-376
    /// `GetOverlappingRange()` (`flags` `eStrandAny`).
    pub(super) fn overlapping(&self) -> Range {
        if self.circular {
            return Range::WHOLE;
        }
        let mut ret = self.plus;
        ret.combine(self.minus);
        if self.single_strand && (self.more_before || self.more_after) {
            let first = self.first_strand();
            if includes_plus(first) {
                if self.more_before {
                    ret.from = 0;
                }
                if self.more_after {
                    ret.to_open = MAX;
                }
            } else {
                if self.more_after {
                    ret.from = 0;
                }
                if self.more_before {
                    ret.to_open = MAX;
                }
            }
        }
        ret
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/handle_range.cpp:442-447
    pub(super) fn has_gaps(&self) -> bool {
        self.ranges.len() > 1 || self.more_before || self.more_after
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/handle_range.cpp:378-417
    /// The strand of the first range (`m_Ranges.front().second`; the callers' ranges are never
    /// empty).
    fn first_strand(&self) -> Strand {
        self.ranges.first().map_or(UNKNOWN, |&(_, strand)| strand)
    }

    /// `GetCircularRangeStart()` (a circular range).
    pub(super) fn circular_start(&self) -> Range {
        let mut ret = self.plus;
        if !is_reverse(self.first_strand()) {
            ret.to_open = MAX;
        } else {
            ret.from = 0;
        }
        ret
    }

    /// `GetCircularRangeEnd()` (a circular range).
    pub(super) fn circular_end(&self) -> Range {
        let mut ret = self.minus;
        if !is_reverse(self.first_strand()) {
            ret.from = 0;
        } else {
            ret.to_open = MAX;
        }
        ret
    }
}

/// `CSeq_id_Handle`: the key under which the object manager keeps a Seq-id's ranges
/// (`CSeq_id_Handle::GetHandle`, the tree of the Seq-id's alternative). Two Seq-ids are one handle
/// exactly when `seqid::mask_handle::mask_handle_key` gives them one key (NCBI's identity of every
/// alternative: accessions, names, local and `general` strings, patent numbers without case,
/// versions and releases exactly, GenBank/EMBL/DDBJ accessions alone across the three; `om/` cases
/// `s_i*` of line RG D3-1); a text Seq-id without accession and name, which the tree finds nowhere,
/// is a new handle at each lookup (`Unique`); a `giim` Seq-id is the handle that NCBI's inverted
/// `giim` tree gives it (`seqid::mask_handle::GiimTree`, kept in `LiveHandles`).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(super) enum Handle {
    /// The handle of `mask_handle_key`, with its `m_Packed` (`mask_handle_packed`).
    Key { key: Vec<u8>, packed: i64 },
    /// A handle that no other lookup finds.
    Unique(u64),
    /// A `giim` handle of the number `id`: the info (numbered) that NCBI's tree gives it.
    Giim { id: i32, info: usize },
}

/// The lookups and infos of Seq-ids that no other lookup finds, numbered.
static LOOKUPS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

fn next_lookup() -> u64 {
    LOOKUPS.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:1186-1203 (a text Seq-id that
// does not pack: `x_FindStrInfo` looks in `m_ByAcc` when the accession is set, else in `m_ByName`,
// else nowhere; the new info is put in the maps of what is set, so one without accession and name
// is never found again: probe p4 r0261, `om/` case `s_i26_pir_empty`)
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:2926-2941, 2953-2964 (a `giim`
// Seq-id finds the first held handle of its number whose Seq-id is NOT `Equals` to it, else makes
// one: two equal `giim` Seq-ids are two handles, two different ones of a number one handle, unless
// other `giim` handles of the number are alive; `om/` cases `s_i72_giim_same`, `s_i73_giim_db_diff`)
impl Handle {
    /// `CSeq_id_Handle::GetHandle` of `id`, a `giim` Seq-id found in `giim`.
    pub(super) fn of(id: &ParsedSeqId, handles: &mut LiveHandles) -> Handle {
        match mask_handle_key(id) {
            Ok(Some(key)) => Handle::Key {
                key,
                packed: mask_handle_packed(id),
            },
            Ok(None) => Handle::Unique(next_lookup()),
            Err(TreeHandle::Giim { id: number }) => {
                let ParsedSeqId::Giim(value) = id else {
                    unreachable!("only a giim Seq-id has a tree handle")
                };
                Handle::Giim {
                    id: number,
                    info: handles.find_or_create(value, Holder::Object),
                }
            }
        }
    }

    /// `CSeq_id_Handle::GetHandle` of a handle that lives as long as the run (a Bioseq's Seq-id,
    /// `CBioseq_Info::m_Id`; a Seq-table column's default, `CSeqTableLocColumns::m_DefaultIdHandle`):
    /// the first such lookup of a keyed handle makes its info (`LiveHandles::made_from`: the info
    /// of an accession is shared with every other number of its prefix, digit count and version).
    pub(super) fn of_run(id: &ParsedSeqId, handles: &mut LiveHandles) -> Handle {
        let handle = Handle::of(id, handles);
        if let (Handle::Key { .. }, Ok(Some(info))) = (&handle, mask_info_key(id)) {
            handles.made_from.entry(info).or_insert_with(|| id.clone());
        }
        handles.mapped(&handle);
        handles.end_object();
        handle
    }

    /// The lookup of `id` by a `CSeq_loc_CI` (its handle lives until `LiveHandles::end_iterator`).
    pub(super) fn look_up_in_iterator(id: &ParsedSeqId, handles: &mut LiveHandles) {
        if let ParsedSeqId::Giim(value) = id {
            handles.find_or_create(value, Holder::Iterator);
        }
    }

    /// `CSeq_id_Handle::GetGiHandle` of `gi`.
    pub(super) fn gi(gi: i64) -> Handle {
        match mask_handle_key(&ParsedSeqId::Gi(gi)) {
            Ok(Some(key)) => Handle::Key {
                key,
                packed: mask_handle_packed(&ParsedSeqId::Gi(gi)),
            },
            _ => unreachable!("a gi is a key"),
        }
    }

    /// The handle's place in `CSeq_id_Handle::operator<` among handles of another `m_Packed`
    /// (`mask_handle_order`); handles of one `m_Packed` follow the addresses of their infos.
    pub(super) fn order(&self) -> u64 {
        mask_handle_order(match self {
            Handle::Key { packed, .. } => *packed,
            Handle::Unique(_) | Handle::Giim { .. } => 0,
        })
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_handle.cpp:89-92 and
// seq_id_tree.cpp:220-236 (an info leaves its tree when its last handle goes:
// `CSeq_id_Info::x_RemoveLastLock` -> `DropInfo` -> `x_Unindex`)
// NCBI reference (598d8ae6): c++/src/objmgr/bioseq_info.cpp:300-307 (the query's Seq-ids are
// looked up when its Bioseq joins the scope, before its annotations; their handles live as long
// as the scope)
/// The Seq-id handles of NCBI's object manager during the query entry's step (`index_entry`): the
/// `giim` tree (`seqid::mask_handle::GiimTree`, cited at `Handle`) with what holds each of its live
/// handles, and the Seq-id that made the info of each handle that lives as long as the run
/// (`made_from`). A lookup in the `giim` tree gets the first live handle of its number whose
/// Seq-id is not `Equals` to it, else a new handle; a handle leaves the tree when the last thing
/// holding it goes. What holds a handle, as NCBI's object manager keeps it:
/// - the run: the query entry's Seq-ids (`CBioseq_Info::m_Id`, made in phase C), a Seq-table
///   column's default (`m_DefaultIdHandle`, phase C), and every handle mapped into the index
///   (`x_Map`: the TSE's annotation maps and the object index keep its key);
/// - the annotation object being indexed: the maps of a feature (location and product, built
///   before they are mapped), a graph, an alignment (all its rows; a `disc` all its parts), a
///   location set, a table location, a table row's location or product map;
/// - a Seq-table's `key` (`x_InitFeatTableKeys`, seq_annot_info.cpp:1087, 1164): a simple row's handle
///   (`SetTableKeyAndIndex`) until the next row puts another one there;
/// - a `CSeq_loc_CI` of a Std-seg location (Seq_loc.cpp:1420-1425), until it goes.
///
/// The other `key` and `hrmaps` that NCBI keeps from one object to the next
/// (seq_annot_info.cpp:726-729, annot_object.cpp:429-444) only hold handles of non-empty ranges,
/// which are mapped (an empty one is LOSAT's refusal, `Indexer::map_location`). The subjects of a
/// `-subject` run hold no `giim` Seq-id: LOSAT refuses a record with one
/// (`blastinput::fasta_reader::believed::refused_id`), as it refuses a PSSM query with one.
#[derive(Debug, Default)]
pub(super) struct LiveHandles {
    tree: GiimTree,
    /// What holds each live `giim` handle: its number and holders, by its info.
    holds: HashMap<usize, (i32, Holds)>,
    /// The Seq-id that made each info of the keyed handles of the run (by `mask_info_key`: the
    /// info of a packing Textseq-id is shared by the accessions of one prefix, digit count and
    /// version), the first of the run to have it (the phase C handles: no other handle is made
    /// before the attach).
    made_from: HashMap<Vec<u8>, ParsedSeqId>,
}

/// What holds a live `giim` handle.
#[derive(Clone, Copy, Debug, Default)]
struct Holds {
    run: bool,
    object: bool,
    key: bool,
    iterator: bool,
}

impl Holds {
    fn hold(&mut self, holder: Holder) {
        match holder {
            Holder::Object => self.object = true,
            Holder::Iterator => self.iterator = true,
        }
    }

    fn any(self) -> bool {
        self.run || self.object || self.key || self.iterator
    }
}

/// What a lookup's handle is held by.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Holder {
    /// The annotation object's maps.
    Object,
    /// A `CSeq_loc_CI`.
    Iterator,
}

impl LiveHandles {
    /// The handles after a query Bioseq with the Seq-ids `query_ids` joined (phase C).
    #[cfg(test)]
    pub(super) fn new(query_ids: &[ParsedSeqId]) -> LiveHandles {
        let mut handles = LiveHandles::default();
        for id in query_ids {
            Handle::of_run(id, &mut handles);
        }
        handles
    }

    /// `CSeq_id_Giim_Tree::FindOrCreate` of the value, the handle held by `holder`: its info.
    fn find_or_create(&mut self, value: &GiimportId, holder: Holder) -> usize {
        let (info, _) = self.tree.find_or_create(value);
        self.holds
            .entry(info)
            .or_insert((value.id, Holds::default()))
            .1
            .hold(holder);
        info
    }

    fn holds(&mut self, handle: &Handle) -> Option<&mut Holds> {
        let Handle::Giim { info, .. } = handle else {
            return None;
        };
        self.holds.get_mut(info).map(|(_, holds)| holds)
    }

    /// `handle` is mapped into the index (`x_Map` of a non-empty range): it lives with the run.
    pub(super) fn mapped(&mut self, handle: &Handle) {
        if let Some(holds) = self.holds(handle) {
            holds.run = true;
        }
    }

    /// Every handle of the maps `maps` is mapped (an object whose ranges are all non-empty).
    pub(super) fn mapped_all<'m>(&mut self, maps: impl IntoIterator<Item = &'m HandleRangeMap>) {
        for map in maps {
            for (handle, _) in &map.map {
                self.mapped(handle);
            }
        }
    }

    /// The Seq-table's `key` takes `handle` (`None`: a handle that lives anyway, or none).
    pub(super) fn set_key(&mut self, handle: Option<&Handle>) {
        for (_, holds) in self.holds.values_mut() {
            holds.key = false;
        }
        if let Some(holds) = handle.and_then(|handle| self.holds(handle)) {
            holds.key = true;
        }
        self.prune();
    }

    /// The end of an annotation object: its maps go.
    pub(super) fn end_object(&mut self) {
        for (_, holds) in self.holds.values_mut() {
            holds.object = false;
        }
        self.prune();
    }

    /// The end of a Seq-table's index: its `key` goes.
    pub(super) fn end_table(&mut self) {
        self.set_key(None);
    }

    /// The end of a `CSeq_loc_CI`.
    pub(super) fn end_iterator(&mut self) {
        for (_, holds) in self.holds.values_mut() {
            holds.iterator = false;
        }
        self.prune();
    }

    fn prune(&mut self) {
        let tree = &mut self.tree;
        self.holds.retain(|&info, &mut (number, holds)| {
            if !holds.any() {
                tree.drop_handle(number, info);
            }
            holds.any()
        });
    }

    // NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_handle.cpp:220-233
    // NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.hpp:359-388 (the info of a packed
    // Textseq-id is shared by its packed key, not by its accession: `seqid::mask_handle::mask_info_key`)
    /// `CSeq_id_Handle::AsString()` of `handle`, the handle of `id` that lives as long as the run:
    /// `GetSeqId()` (`seqid::mask_handle::handle_seq_id` of the Seq-id that made its info) in
    /// FASTA form (`gi|` and the number for a gi).
    pub(super) fn as_string(&self, handle: &Handle, id: &ParsedSeqId) -> String {
        let seq_id = match handle {
            Handle::Key { .. } => match mask_info_key(id) {
                Ok(Some(info)) => self
                    .made_from
                    .get(&info)
                    .map_or_else(|| id.clone(), |first| handle_seq_id(first, id)),
                _ => id.clone(),
            },
            Handle::Giim { id: number, info } => {
                ParsedSeqId::Giim(self.tree.made_from(*number, *info).clone())
            }
            Handle::Unique(_) => id.clone(),
        };
        let mut out = Vec::new();
        crate::seqid::label::fasta(&seq_id, &mut out);
        String::from_utf8_lossy(&out).into_owned()
    }
}

/// `CHandleRangeMap::ETransSplicing`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Mode {
    /// `eNoTransSplicing`.
    None,
    /// `eTransSplicing`.
    TransSplicing,
    /// `eCircularRNA`.
    CircularRna,
}

/// `CHandleRangeMap`: a location's ranges by Seq-id, in the order the Seq-ids came; NCBI's map
/// iterates them by `CSeq_id_Handle::operator<` (`HandleRangeMap::ncbi_order`).
#[derive(Default, Debug)]
pub(super) struct HandleRangeMap {
    pub(super) map: Vec<(Handle, HandleRange)>,
}

/// `SAddState`: the Seq-id of the range added before, and the mode; with the `giim` tree of the
/// lookups.
struct AddState<'g> {
    prev: Option<usize>,
    mode: Mode,
    handles: &'g mut LiveHandles,
}

impl HandleRangeMap {
    /// `AddLocation` of `loc` in `mode` on a fresh map.
    pub(super) fn of_location(loc: V<'_>, mode: Mode, handles: &mut LiveHandles) -> HandleRangeMap {
        let mut map = HandleRangeMap::default();
        map.add_location(
            loc,
            &mut AddState {
                prev: None,
                mode,
                handles,
            },
        );
        map
    }

    fn index_of(&mut self, handle: Handle) -> usize {
        match self.map.iter().position(|(h, _)| *h == handle) {
            Some(at) => at,
            None => {
                self.map.push((handle, HandleRange::default()));
                self.map.len() - 1
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/handle_range_map.cpp:246-281
    // The query is a raw Bioseq: its TSE has no master segments (`m_MasterSeq` null).
    fn add(&mut self, handle: Handle, range: Range, strand: Strand, state: &mut AddState<'_>) {
        let at = self.index_of(handle);
        if state.mode == Mode::None {
            if let Some(prev) = state.prev.filter(|&prev| prev != at) {
                if let Some((_, hr)) = self.map.get_mut(prev) {
                    hr.more_after = true;
                }
                if let Some((_, hr)) = self.map.get_mut(at) {
                    hr.more_before = true;
                }
            }
        }
        if let Some((_, hr)) = self.map.get_mut(at) {
            hr.add(range, strand, state.mode == Mode::CircularRna);
        }
        state.prev = Some(at);
    }

    /// `AddRange` of an id, `range` and `strand` with a state of its own (no Seq-id before it): the
    /// rows of an alignment (annot_object.cpp:591-842).
    pub(super) fn add_range(
        &mut self,
        id: &ParsedSeqId,
        range: Range,
        strand: Strand,
        handles: &mut LiveHandles,
    ) {
        let handle = Handle::of(id, handles);
        self.add_handle(handle, range, strand, handles);
    }

    /// `AddRange` of `handle`, `range` and `strand` with a state of its own.
    pub(super) fn add_handle(
        &mut self,
        handle: Handle,
        range: Range,
        strand: Strand,
        handles: &mut LiveHandles,
    ) {
        self.add(
            handle,
            range,
            strand,
            &mut AddState {
                prev: None,
                mode: Mode::None,
                handles,
            },
        );
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/handle_range_map.cpp:91-205
    fn add_location<'v>(&mut self, loc: V<'v>, state: &mut AddState<'_>) {
        // The `mix` and `equiv` locations being walked (to any depth the reader read): the
        // location's parts, the next one, and whether it is an `equiv`.
        let mut sets: Vec<(V<'v>, usize, bool)> = Vec::new();
        let mut next = Some(loc);
        loop {
            let loc = match next.take() {
                Some(loc) => loc,
                None => {
                    let Some((parts, index, equiv)) = sets.last_mut() else {
                        return;
                    };
                    let Some(part) = parts.item(*index) else {
                        sets.pop();
                        continue;
                    };
                    if *equiv && *index > 0 {
                        state.prev = None;
                    }
                    *index += 1;
                    part
                }
            };
            let (variant, value) = loc.choice();
            match variant {
                "empty" => {
                    let handle = Handle::of(value.seq_id(), state.handles);
                    self.add(handle, Range::EMPTY, UNKNOWN, state)
                }
                "whole" => {
                    let handle = Handle::of(value.seq_id(), state.handles);
                    self.add(handle, Range::WHOLE, UNKNOWN, state)
                }
                "int" => self.add_interval(value, state),
                "packed-int" => {
                    for interval in value.items() {
                        self.add_interval(interval, state);
                    }
                }
                "pnt" => self.add_point(value, state),
                "packed-pnt" => {
                    // One handle for all the points (`idh`).
                    let handle = Handle::of(value.field("id").seq_id(), state.handles);
                    let strand = value.get("strand").map_or(UNKNOWN, V::int);
                    for point in value.field("points").ints() {
                        let point = point as u32;
                        self.add(handle.clone(), Range::new(point, point), strand, state);
                    }
                }
                "mix" => sets.push((value, 0, false)),
                "equiv" => sets.push((value, 0, true)),
                "bond" => {
                    self.add_point(value.field("a"), state);
                    if let Some(b) = value.get("b") {
                        self.add_point(b, state);
                    }
                }
                // `null`, `feat`.
                _ => {}
            }
        }
    }

    fn add_interval(&mut self, interval: V<'_>, state: &mut AddState<'_>) {
        let (id, range, strand) = interval_parts(interval);
        let handle = Handle::of(id, state.handles);
        self.add(handle, range, strand, state);
    }

    fn add_point(&mut self, point: V<'_>, state: &mut AddState<'_>) {
        let (id, range, strand) = point_parts(point);
        let handle = Handle::of(id, state.handles);
        self.add(handle, range, strand, state);
    }

    // NCBI reference (598d8ae6): c++/include/objmgr/impl/handle_range_map.hpp:91 (the map's order
    // is `CSeq_id_Handle`'s `operator<`)
    // NCBI reference (598d8ae6): c++/include/objects/seq/seq_id_handle.hpp:201-207
    /// The map's entries in NCBI's order, in groups of one `m_Packed`: the
    /// groups in NCBI's order, the handles of a group in the order of their infos' addresses,
    /// which LOSAT does not know (it keeps the order they came in). A gi and the packed accessions
    /// and `general` tags are ordered by their numbers (`om/` cases `s_o*`); the handles of one
    /// number, and every unpacked handle (local, PDB, patent, a text Seq-id with a name or a
    /// release or whose accession does not pack, ...), are one group.
    pub(super) fn ncbi_order(&self) -> Vec<Vec<&(Handle, HandleRange)>> {
        let mut entries: Vec<(u64, &(Handle, HandleRange))> = self
            .map
            .iter()
            .map(|entry| (entry.0.order(), entry))
            .collect();
        entries.sort_by_key(|&(order, _)| order);
        entries
            .chunk_by(|a, b| a.0 == b.0)
            .map(|group| group.iter().map(|&(_, entry)| entry).collect())
            .collect()
    }
}

/// A Seq-interval's Seq-id, range (a `CRange` from `from` to `to`) and strand (`unknown` when not
/// set).
pub(super) fn interval_parts(interval: V<'_>) -> (&ParsedSeqId, Range, Strand) {
    let from = interval.field("from").int() as u32;
    let to = interval.field("to").int() as u32;
    (
        interval.field("id").seq_id(),
        Range::new(from, to),
        interval.get("strand").map_or(UNKNOWN, V::int),
    )
}

/// A Seq-point's Seq-id, range (a `CRange` of `point` alone) and strand.
pub(super) fn point_parts(point: V<'_>) -> (&ParsedSeqId, Range, Strand) {
    let at = point.field("point").int() as u32;
    (
        point.field("id").seq_id(),
        Range::new(at, at),
        point.get("strand").map_or(UNKNOWN, V::int),
    )
}

#[cfg(test)]
#[path = "range_tests.rs"]
mod range_tests;
