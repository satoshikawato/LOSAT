//! The Seq-ids of the queries and subjects of one search under `-parse_deflines`, which NCBI keeps
//! in one scope: it finds a sequence again by its id, so a query and a subject with one id show as
//! one sequence, and a database sequence with the id of a query shows as that query
//! (`dropin-survey-20261029/DI7/DECISIONS.md` 8, 12).
//!
//! A record given as both query and subject (its defline and letters: `RecordMark`) is one
//! sequence to the scope too, and the reports are those of two records, except where outfmt 0
//! draws its masks (`DECISIONS.md` 15, 19). A record of a nucleotide and of a protein reader never
//! is one record: NCBI stops with an exception when it reports such a hit (`eBadAlignment`).

use std::collections::HashMap;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::sync::Mutex;

use super::defline_set::scope_keys;
use super::label;
use super::model::ParsedSeqId;

/// What tells a record from another of the same ids, taken while it is read.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RecordMark {
    /// Of its molecule, its defline as written (the archive shows the query's for both: `>a` and
    /// `>lcl|a` are one id, not one record) and its letters in upper case (`U` as `T` in a
    /// nucleotide record: NCBI's packed letters do not tell them apart).
    digest: u64,
    /// Whether it has lower-case letters.
    lowercase: bool,
}

impl RecordMark {
    /// The mark of the record with the defline `defline`, whose letters are `letters` in upper case
    /// and `gaps` (position, length) between them.
    pub fn new(
        defline: &[u8],
        protein: bool,
        letters: &[u8],
        gaps: impl Iterator<Item = (u32, u32)>,
        lowercase: bool,
    ) -> Self {
        let mut hasher = DefaultHasher::new();
        protein.hash(&mut hasher);
        defline.hash(&mut hasher);
        for gap in gaps {
            gap.hash(&mut hasher);
        }
        // `write` takes the bytes as one stream, in any pieces.
        let mut piece = [0u8; 4096];
        for part in letters.chunks(piece.len()) {
            for (out, &letter) in piece.iter_mut().zip(part) {
                *out = if letter == b'U' && !protein {
                    b'T'
                } else {
                    letter
                };
            }
            hasher.write(&piece[..part.len()]);
        }
        Self {
            digest: hasher.finish(),
            lowercase,
        }
    }
}

/// A record that a reader of the other role has read too.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SharedRecord {
    /// Whether one of the two copies has lower-case letters.
    pub lowercase: bool,
}

/// The ids of the records that the readers of one search have read (`IdScope::add_record`).
#[derive(Debug, Default)]
pub struct IdScope {
    /// By `scope_keys` and `version_keys`' own key: the queries, the subjects. A record without a
    /// mark is never one record with another: BLASTX's, whose reader takes none, and a record whose
    /// id is generated or local `Subject_N`, which the reports show another way.
    records: Mutex<[HashMap<u64, Option<RecordMark>>; 2]>,
}

impl IdScope {
    /// LOSAT's refusal of a record with the ids `ids` and the mark `mark` when a record of the other
    /// role has one of the ids and is not the same record; else that record, if there is one. Adds
    /// the ids.
    pub fn add_record(
        &self,
        ids: &[ParsedSeqId],
        mark: Option<RecordMark>,
        role: &str,
        program: &str,
    ) -> anyhow::Result<Option<SharedRecord>> {
        let mut records = self
            .records
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        let [queries, subjects] = &mut *records;
        let (own, other) = if role == "query" {
            (queries, subjects)
        } else {
            (subjects, queries)
        };
        let mut shared = None;
        let mut keys = Vec::new();
        for id in ids {
            keys.clear();
            scope_keys(id, &mut keys);
            let versions = version_keys(id);
            let names = name_keys(id);
            let found = keys.iter().find_map(|key| other.get(key));
            // The same accession with and without a version, or name with and without an
            // accession, is another id, whatever the record.
            let other_kind = versions.is_some_and(|(_, kind)| other.contains_key(&kind))
                || names.is_some_and(|(_, kind)| other.contains_key(&kind));
            if found.is_some() || other_kind {
                let copy = found
                    .copied()
                    .flatten()
                    .zip(mark)
                    .filter(|(theirs, mine)| !other_kind && theirs.digest == mine.digest);
                let Some((theirs, mine)) = copy else {
                    let mut fasta = Vec::new();
                    label::fasta(id, &mut fasta);
                    anyhow::bail!(
                        "a query and a subject record have the Seq-id {} under -parse_deflines \
                         (NCBI BLAST+ shows them as one sequence); this is not supported by \
                         LOSAT's {program}",
                        String::from_utf8_lossy(&fasta)
                    );
                };
                shared = Some(SharedRecord {
                    lowercase: theirs.lowercase || mine.lowercase,
                });
            }
            own.extend(keys.iter().map(|&key| (key, mark)));
            own.extend(versions.map(|(kind, _)| (kind, mark)));
            own.extend(names.map(|(kind, _)| (kind, mark)));
        }
        Ok(shared)
    }

    /// Whether a query read so far has an id with one of the keys `keys`. NCBI would show a database
    /// sequence with such an id as the query (`report::db_subject::ParsedSubject`). A query of
    /// another batch counts too, though NCBI's scope holds one batch: a refusal, not another report.
    pub fn query_has_key(&self, keys: &[u64]) -> bool {
        let records = self
            .records
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        keys.iter().any(|key| records[0].contains_key(key))
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
// ```c++
// SSeqMatch_Scope CDataSource_ScopeInfo::x_GetSeqMatch(const CSeq_id_Handle& idh)
// {
//     SSeqMatch_Scope ret = x_FindBestTSE(idh);
//     if ( !ret && idh.HaveMatchingHandles() ) {
//         CSeq_id_Handle::TMatches ids;
//         idh.GetMatchingHandles(ids, eAllowWeakMatch);
//         ITERATE ( CSeq_id_Handle::TMatches, it, ids ) {
//             ...
//             if ( SSeqMatch_Scope match = x_FindBestTSE(*it) ) {
//                 ret = match;
//             }
// ```
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:1494-1507
// ```c++
//         if ( tid->IsSetAccession() ) {
//             x_FindMatchByAcc(id_list, tid->GetAccession(), tid);
//         }
//         if ( tid->IsSetName() ) {
//             x_FindMatchByName(id_list, tid->GetName(), tid);
//         }
//         // special case for PIR and PRF ids - matching accesstion to name
//         // the id should have only accession and no other fields (version, name, release)
//         if ( (type == CSeq_id::e_Pir || type == CSeq_id::e_Prf) &&
//              tid->IsSetAccession() &&
//              !tid->IsSetVersion() && !tid->IsSetName() && !tid->IsSetRelease() &&
//              !m_ByName.empty() ) {
//             x_FindMatchByName(id_list, tid->GetAccession(), nullptr, eOnlyName);
//         }
// ```
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
    let key = |versioned: bool| {
        let mut hasher = DefaultHasher::new();
        (b"version kind", versioned, &accession).hash(&mut hasher);
        hasher.finish()
    };
    let versioned = text.version.is_some();
    Some((key(versioned), key(!versioned)))
}
