//! How two Seq-ids relate, as NCBI's formatters and its scope ask it: the one owner of these
//! questions (line RG-KB item 12d; decisions KB-12c-3 and XT-X7-split-1 of
//! `remaining-gaps-20261112/DECISIONS.md`).
//!
//! - [`matches`]: `CSeq_id::Match`, `defline_set::compare` (the port of `CSeq_id::Compare`) giving
//!   `Sic::Yes`: a text Seq-id's accession, else its name, in any letter case, the versions
//!   compared only where both ids have one, across the GenBank-like alternatives; every other
//!   alternative by its own `Match`.
//! - [`row_matches`] and [`RowPairs`]: `matches` of the Seq-ids of two reported sequences
//!   (`report::seqid::SeqId`: a generated local id, a database's ordinal id, a parsed id), what
//!   `CDisplaySeqalign::x_FillLocList` asks of a mask's id and an alignment row's: a mask goes to
//!   the first row whose id matches its own and whose range meets it
//!   (`report::query_anchored::fill_row_loc_list`, `fill_loc_list`).
//! - The scope's finding: [`scope_finds`], [`ScopeIndex`] and the keys they index
//!   ([`weak_keys`], [`kind_keys`], [`version_keys`], [`name_keys`], [`match_form`]; with
//!   `defline_set::scope_keys`, an id's own key): the Bioseq that NCBI's scope finds for an id, by
//!   the id itself or weakly.
//!
//! The questions of a reported subject against the queries of its scope (whether the scope finds
//! it as a query, by its own or a weak key, and by its display id) come here too (lines XT X7-2,
//! RG-KB X7-3). `seqid::id_scope` re-exports the scope's keys for its callers.

use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::hash::{DefaultHasher, Hash, Hasher};

use super::defline_set::{compare, scope_keys, Sic};
use super::model::{Dbtag, ObjectId, ParsedSeqId};
use crate::report::seqid::SeqId;

// NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_id.hpp:1065-1070
/// `CSeq_id::Match`: whether `CSeq_id::Compare` (`defline_set::compare`) finds `a` and `b` the
/// same.
pub(crate) fn matches(a: &ParsedSeqId, b: &ParsedSeqId) -> bool {
    compare(a, b) == Sic::Yes
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4305-4306
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:364-367
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2221-2238
/// Whether the Seq-ids `a` and `b` of two reported sequences match (`CSeq_id::Match`, [`matches`]):
/// a generated id (`Query_N`, `Subject_N`) as the local id of that string, a sequence of a
/// database without parsed Seq-ids as the `general` id of the database `BL_ORD_ID` and its OID.
pub(crate) fn row_matches(a: SeqId<'_>, b: SeqId<'_>) -> bool {
    match (a, b) {
        // Two `general` ids of one database: their tags, two numbers (`CDbtag::Match`).
        (SeqId::BlOrdId(a), SeqId::BlOrdId(b)) => a == b,
        // A `local` and a `general` id: `Compare` gives `e_DIFF`.
        (SeqId::Local(_), SeqId::BlOrdId(_)) | (SeqId::BlOrdId(_), SeqId::Local(_)) => false,
        _ => matches(&as_parsed(a), &as_parsed(b)),
    }
}

/// The Seq-id `id` as a [`ParsedSeqId`] (`row_matches`).
fn as_parsed(id: SeqId<'_>) -> Cow<'_, ParsedSeqId> {
    match id {
        SeqId::Db(id) => Cow::Borrowed(id),
        SeqId::Local(text) => {
            Cow::Owned(ParsedSeqId::Local(ObjectId::Str(text.as_bytes().to_vec())))
        }
        SeqId::BlOrdId(oid) => Cow::Owned(ParsedSeqId::General(Dbtag {
            db: b"BL_ORD_ID".to_vec(),
            // An OID is below 2^31 (an `Int4` in NCBI, `report::seqid::SeqId::write`).
            tag: ObjectId::Id(oid.get() as i32),
        })),
    }
}

// No NCBI counterpart: the pairs whose rows match, found once for a report, so that the rows of
// each HSP ask a set (empty in nearly every search) and not the ids.
/// The pairs (query, subject) of a report whose Seq-ids match ([`row_matches`]): the pairs of
/// whose HSPs `x_FillLocList` may place a query's mask on the subject's row, and a subject's on the
/// query's.
#[derive(Clone, Debug, Default)]
pub(crate) struct RowPairs {
    pairs: HashSet<(u32, u32)>,
}

impl RowPairs {
    /// The pairs of `pairs` (query and subject index, an HSP's; a pair may come again) whose ids
    /// `query_id` and `subject_id` match.
    pub(crate) fn of<'i>(
        pairs: impl IntoIterator<Item = (u32, u32)>,
        query_id: impl Fn(u32) -> SeqId<'i>,
        subject_id: impl Fn(u32) -> SeqId<'i>,
    ) -> Self {
        let mut found = HashSet::new();
        // The HSPs of a pair come together: each pair is asked once per run of them.
        let mut last = None;
        for pair in pairs {
            if last == Some(pair) {
                continue;
            }
            last = Some(pair);
            if row_matches(query_id(pair.0), subject_id(pair.1)) {
                found.insert(pair);
            }
        }
        Self { pairs: found }
    }

    /// Whether the rows of query `q_idx` and subject `s_idx` match.
    #[inline(always)]
    pub(crate) fn matches(&self, q_idx: u32, s_idx: u32) -> bool {
        !self.pairs.is_empty() && self.pairs.contains(&(q_idx, s_idx))
    }
}

/// For a text Seq-id with a name, the key of its name with its kind (with or without an accession)
/// and the key of the other kind. NCBI's scope finds `gb|AAA00007.1|NAMEA` by `gb||NAMEA`, though
/// `scope_keys` tells them apart and `gb|AAA00008.1|NAMEA` is another sequence (probe of the review
/// of DECISIONS 15, 2026-11-02).
pub(crate) fn name_keys(id: &ParsedSeqId) -> Option<(u64, u64)> {
    let text = id.textseq()?;
    let name = text.name.as_ref()?.to_ascii_uppercase();
    let accession = text.accession.is_some();
    Some((
        name_kind_key(&name, accession),
        name_kind_key(&name, !accession),
    ))
}

/// The key of the name `name` (in upper case) of a text Seq-id with or without an accession
/// (`name_keys`).
fn name_kind_key(name: &[u8], accession: bool) -> u64 {
    let mut hasher = DefaultHasher::new();
    (b"name kind", accession, name).hash(&mut hasher);
    hasher.finish()
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:865-882
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:1494-1507
// `x_FindMatchByAcc` takes every version of an accession without one (and the packed path at
// 1469-1479 alike), `x_FindMatchByName` with an accession only the name-only ids (1367-1372),
// and without one any id of that name.
/// Adds to `keys` the keys (`version_keys`, `name_keys`) of the text Seq-ids that NCBI's scope
/// finds by `id` where no sequence it holds has `id` itself (`GetMatchingHandles` with
/// `eAllowWeakMatch`): for an accession without a version, that accession with one; for a name
/// with an accession, the name alone; for a name alone, the name with an accession; for a `pir` or
/// `prf` accession alone, that text as a name alone; in any letter case and across the alternatives
/// (`CSeq_id_Mapper::GetMatchingHandles` asks every text tree), without the release that NCBI also
/// compares for a name (a superset). A sequence that holds one of them has it among
/// its `kind_keys`. NCBI 2.17.0 shows a query read with `-parse_deflines` that has one for a
/// reported database sequence, and, where it reads the reported sequences one by one
/// (`PRE_FETCH_SEQS_LIMIT`), another reported sequence read before (probe of line DBID N2c,
/// `blastdb-ids-20261021/R/weak/p1.out`). Not the other way: `gb|X|` for `gb|X.1|`, `pir|ABCD|`
/// for `pir||ABCD`.
pub(crate) fn weak_keys(id: &ParsedSeqId, keys: &mut Vec<u64>) {
    let Some(text) = id.textseq() else {
        return;
    };
    if text.version.is_none() {
        keys.extend(version_keys(id).map(|(_, versioned)| versioned));
    }
    keys.extend(name_keys(id).map(|(_, other)| other));
    if let (ParsedSeqId::Pir(_) | ParsedSeqId::Prf(_), Some(accession), None, None, None) =
        (id, &text.accession, text.version, &text.name, &text.release)
    {
        keys.push(name_kind_key(&accession.to_ascii_uppercase(), false));
    }
}

/// Adds to `keys` the keys of `id` that `weak_keys` gives for the ids that find it.
pub(crate) fn kind_keys(id: &ParsedSeqId, keys: &mut Vec<u64>) {
    keys.extend(version_keys(id).map(|(own, _)| own));
    keys.extend(name_keys(id).map(|(own, _)| own));
}

/// For a text Seq-id with an accession, the key of its accession with its kind of version (with or
/// without one) and the key of the other kind. NCBI's BLASTN writes `CScope::GetSynonyms: ... is
/// resolved to another Bioseq` for the pair (probe `db.pd.blastn.ver`, 2026-10-31).
pub(crate) fn version_keys(id: &ParsedSeqId) -> Option<(u64, u64)> {
    let text = id.textseq()?;
    let accession = text.accession.as_ref()?.to_ascii_uppercase();
    let versioned = text.version.is_some();
    Some((
        version_kind_key(&accession, versioned),
        version_kind_key(&accession, !versioned),
    ))
}

/// One key of `version_keys` of a text Seq-id with an accession: of its own kind of version
/// (`own`), or of the other kind.
pub(crate) fn version_key(id: &ParsedSeqId, own: bool) -> Option<u64> {
    let text = id.textseq()?;
    let accession = text.accession.as_ref()?.to_ascii_uppercase();
    Some(version_kind_key(&accession, text.version.is_some() == own))
}

/// The key of the accession `accession` (in upper case) of a text Seq-id with or without a version
/// (`version_keys`).
fn version_kind_key(accession: &[u8], versioned: bool) -> u64 {
    let mut hasher = DefaultHasher::new();
    (b"version kind", versioned, accession).hash(&mut hasher);
    hasher.finish()
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Textseq_id.cpp:141-156
// NCBI reference (598d8ae6): c++/src/objects/seqloc/PDB_seq_id.cpp:53-67
/// The bit of the form of a Seq-id among the Seq-ids that `CSeq_id::Match` takes for one though
/// `defline_set::scope_keys` tells them apart (0 for any other Seq-id): for a text Seq-id with an
/// accession, its accession with or without a version (`MATCH_FORMS_ACCESSION`); for a PDB Seq-id,
/// its molecule with or without a `chain` and a `chain_id` (`MATCH_FORMS_PDB`). Two Seq-ids that
/// `defline_set::compare` finds the same (`Sic::Yes`) with different `scope_keys` share the
/// accession (`version_keys`) or the molecule and differ in the form: ids of one form match only
/// where they are equal there (one accession with one version, one molecule with the same chains).
pub(crate) fn match_form(id: &ParsedSeqId) -> u8 {
    match id {
        ParsedSeqId::Pdb(pdb) => {
            1 << (2 + u8::from(pdb.chain.is_some()) + 2 * u8::from(pdb.chain_id.is_some()))
        }
        _ => id
            .textseq()
            .filter(|text| text.accession.is_some())
            .map_or(0, |text| 1 << u8::from(text.version.is_some())),
    }
}

/// The forms of `match_form` of a text Seq-id's accession.
pub(crate) const MATCH_FORMS_ACCESSION: u8 = 0b11;

/// The forms of `match_form` of a PDB Seq-id.
pub(crate) const MATCH_FORMS_PDB: u8 = 0b11_1100;

// NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:865-882
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:1253-1290
/// The Bioseqs of one data source of NCBI's scope (an archive's queries, its `-subject` Bioseqs)
/// by the keys of their ids, for the lookups of `CDataSource_ScopeInfo::x_GetSeqMatch`
/// ([`ScopeIndex::finds`]): `defline_set::scope_keys` of each id (the Bioseqs that hold the id
/// itself: a text Seq-id's accession with its version, so `gb|X.1|` is not `gb|X|`) and
/// `kind_keys` of each id (the Bioseqs that an id matches weakly, `weak_keys`).
#[derive(Debug, Default)]
pub(crate) struct ScopeIndex {
    own: HashMap<u64, Vec<usize>>,
    weak: HashMap<u64, Vec<usize>>,
}

impl ScopeIndex {
    /// Adds the Bioseq `item` (not below a Bioseq added before) with the id `id`.
    pub(crate) fn add(&mut self, item: usize, id: &ParsedSeqId) {
        let mut keys = Vec::new();
        scope_keys(id, &mut keys);
        Self::insert(&mut self.own, &keys, item);
        keys.clear();
        kind_keys(id, &mut keys);
        Self::insert(&mut self.weak, &keys, item);
    }

    fn insert(map: &mut HashMap<u64, Vec<usize>>, keys: &[u64], item: usize) {
        for &key in keys {
            let items = map.entry(key).or_default();
            if items.last() != Some(&item) {
                items.push(item);
            }
        }
    }

    /// The Bioseqs that NCBI's scope finds for `id` in this data source, in order: those with an
    /// id that the scope takes for `id` itself (`x_FindBestTSE`), else those with an id that `id`
    /// matches weakly (`GetMatchingHandles` with `eAllowWeakMatch`: an accession without a version
    /// finds every version of it, one with a version only that version, `x_FindMatchByAcc`).
    /// Of several, NCBI takes the better weak match (`IsBetter`) or throws its conflict of TSEs;
    /// the callers refuse them.
    pub(crate) fn finds(&self, id: &ParsedSeqId) -> Vec<usize> {
        let mut keys = Vec::new();
        scope_keys(id, &mut keys);
        let found = Self::found(&self.own, &keys);
        if !found.is_empty() {
            return found;
        }
        keys.clear();
        weak_keys(id, &mut keys);
        Self::found(&self.weak, &keys)
    }

    fn found(map: &HashMap<u64, Vec<usize>>, keys: &[u64]) -> Vec<usize> {
        let mut found: Vec<usize> = keys
            .iter()
            .filter_map(|key| map.get(key))
            .flatten()
            .copied()
            .collect();
        found.sort_unstable();
        found.dedup();
        found
    }
}

/// Whether NCBI's scope finds a Bioseq that has the id `held` for the id `id` (the lookup of
/// [`ScopeIndex::finds`] in a data source of that one Bioseq).
pub(crate) fn scope_finds(held: &ParsedSeqId, id: &ParsedSeqId) -> bool {
    let (mut held_keys, mut id_keys) = (Vec::new(), Vec::new());
    scope_keys(held, &mut held_keys);
    scope_keys(id, &mut id_keys);
    if held_keys.iter().any(|key| id_keys.contains(key)) {
        return true;
    }
    held_keys.clear();
    id_keys.clear();
    kind_keys(held, &mut held_keys);
    weak_keys(id, &mut id_keys);
    held_keys.iter().any(|key| id_keys.contains(key))
}

#[cfg(test)]
#[path = "relation_tests.rs"]
mod tests;
