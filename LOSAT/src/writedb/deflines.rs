//! The definition lines of a sequence read with `-parse_seqids`, as CWriteDB builds them
//! (`CWriteDB_Impl::ExtractBioseqDeflines`, `c++/src/objtools/blast/seqdb_writer/
//! writedb_impl.cpp`): the second parse of the raw definition line, split at ^A
//! (`x_GetFastaReaderDeflines`), the choice between it and the Bioseq's ids
//! (`s_UseFastaReaderDeflines`, `x_BuildDeflinesFromBioseq`), and the edits of
//! `CBuildDatabase` around it (the ids of the other molecule, `FixTaxId`).

use super::bioseq::{Bioseq, ParsedIds};
use super::fasta_source::{fasta_reader_fasta, post_seq_id_warning, Warn};
use crate::seqid::accguide::{identify_seq_id, parse_flags, F_NUC, F_PROT};
use crate::seqid::label::content;
use crate::seqid::model::ParsedSeqId;
use crate::seqid::parse::{parse_ids, seq_id_from_text, set_typed, SeqIdError};
use crate::seqid::rank::{best_rank_score, find_best};
use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::{member_index, members, Value};
use crate::report::serial::pssm::Seqdesc;
use crate::seqid::{BlastDefline, DeflineSet};

/// `GetSeqIdString(true)`.
fn seq_id_string_with_version(id: &ParsedSeqId) -> Vec<u8> {
    let mut out = Vec::new();
    content(id, true, &mut out);
    out
}

/// `CSeq_id(CSeq_id::e_Local, label)`: the typed `Set` of a local id (`SetStrOrId` of the
/// trimmed label), which does not fail.
fn local_from_label(label: &[u8]) -> ParsedSeqId {
    set_typed(crate::seqid::accguide::SEQID_LOCAL, label, b"", 0, b"")
        .ok()
        .flatten()
        .expect("a local id is always made")
}

/// Whether `IdentifyAccession()` of `id` names the other molecule than `is_protein`:
/// `(info & (fAcc_prot | fAcc_nuc)) && is_protein == !!(info & fAcc_nuc)`.
fn of_other_molecule(id: &ParsedSeqId, is_protein: bool, warn: Warn<'_>) -> bool {
    let info = identify_seq_id(
        id,
        parse_flags::ANY_RAW | parse_flags::ANY_LOCAL,
        &mut |w| post_seq_id_warning(w, &mut *warn),
    );
    info & (F_PROT | F_NUC) != 0 && is_protein == (info & F_NUC != 0)
}

/// The ids edit of `CBuildDatabase::AddSequences` (without `[BLAST] LONG_SEQID`): every
/// id that is neither local nor a GI and whose accession is of the other molecule
/// becomes the local id of its label.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:840-856
/// `IdentifyAccession()` runs for every id, local and GI ones included (its diagnostic,
/// if any, is posted for them too).
pub fn localize_other_molecule(ids: &mut [ParsedSeqId], is_protein: bool, warn: Warn<'_>) {
    for id in ids {
        let other = of_other_molecule(id, is_protein, warn);
        if !matches!(id, ParsedSeqId::Local(_) | ParsedSeqId::Gi(_)) && other {
            *id = local_from_label(&seq_id_string_with_version(id));
        }
    }
}

/// The GI check of `CBuildDatabase::AddSequences`: when an id is a GI above 4294967295,
/// the sequence is skipped and the GI of the last GI id is logged.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:812-833
/// `IdentifyAccession()` of a GI is `eAcc_gi` (its `E_Choice`), of no other id. The
/// diagnostics of `IdentifyAccession` for the other ids are posted to `warn`.
pub fn large_gi_to_skip(ids: &[ParsedSeqId], warn: Warn<'_>) -> Option<i64> {
    let mut skip_this = false;
    let mut check_gi = 0;
    for id in ids {
        identify_seq_id(
            id,
            parse_flags::ANY_RAW | parse_flags::ANY_LOCAL,
            &mut |w| post_seq_id_warning(w, &mut *warn),
        );
        if let ParsedSeqId::Gi(gi) = *id {
            check_gi = gi;
            if gi > 0xFFFF_FFFF {
                skip_this = true;
            }
        }
    }
    skip_this.then_some(check_gi)
}

/// The second parse of `x_GetFastaReaderDeflines(bioseq, ..., accept_gt = false,
/// parse_ids = true, long_seqids = false, scan = true)`: one definition line per part of
/// the raw definition line between ^A; `None` where NCBI's parse throws a
/// `CSeqIdException` (which `x_ExtractDeflines` catches) or finds no definition line.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1766-1767,1807-1880,1896-1900
/// `ParseFastaIds(seqids, ids)` is `ParseIDs` with `fParse_RawText | fParse_AnyLocal`
/// (no `fParse_PartialOK`, Seq_id.cpp:2603-2611): any exception goes up. No membership,
/// link or PIG is given (`ExtractBioseqDeflines` passes none). The `fasta` string is the
/// one `fasta_reader_title` describes; its escapes do not fail.
pub fn second_parse(parsed: &ParsedIds, is_aa: bool, warn: Warn<'_>) -> Option<DeflineSet> {
    let mut fasta = fasta_reader_fasta(&parsed.defline_field).ok()?;
    if fasta.is_empty() {
        return None;
    }
    let find = |text: &[u8], byte: u8, from: usize| {
        text.get(from..)
            .and_then(|rest| rest.iter().position(|&c| c == byte))
            .map(|offset| from + offset)
    };
    let mut deflines = Vec::new();
    let mut skip = 1usize;
    while !fasta.is_empty() {
        let id_start = skip;
        let mut pos_title = find(&fasta, b' ', skip);
        let mut pos_next = find(&fasta, 0x01, skip);
        skip = 1;
        if pos_next.is_none() {
            pos_next = Some(fasta.len());
            skip = 0;
        }
        let pos_next = pos_next.expect("set above");
        if pos_title.is_none_or(|pos| pos >= pos_next) {
            pos_title = Some(pos_next);
        }
        let mut pos_title = pos_title.expect("set above");
        let ids = fasta[id_start.min(pos_title)..pos_title].to_vec();
        if pos_title == pos_next {
            // `pos_title--` of 0 wraps; the title is then from 0 (`pos_title + 1`).
            pos_title = pos_title.wrapping_sub(1);
        }
        let title_start = pos_title.wrapping_add(1);
        let title = fasta[title_start.min(pos_next)..pos_next].to_vec();
        fasta.drain(..pos_next);
        let seqids = match parse_second_ids(&ids, is_aa, warn) {
            Ok(seqids) => seqids,
            Err(_) => return None,
        };
        deflines.push(BlastDefline {
            title: Some(title),
            ids: seqids,
            ..BlastDefline::default()
        });
    }
    Some(DeflineSet { deflines })
}

/// The ids of one part of the second parse (`x_GetFastaReaderDeflines`, cited there).
fn parse_second_ids(
    ids: &[u8],
    is_aa: bool,
    warn: Warn<'_>,
) -> Result<Vec<ParsedSeqId>, SeqIdError> {
    let mut on_warning = |w| post_seq_id_warning(w, &mut *warn);
    if ids.contains(&b'|') {
        let parsed = parse_ids(
            ids,
            parse_flags::RAW_TEXT | parse_flags::ANY_LOCAL,
            &mut on_warning,
        )?;
        // Without fParse_PartialOK no id stays unset (an exception goes up first).
        return Ok(parsed.into_iter().flatten().collect());
    }
    let mut id = seq_id_from_text(
        ids,
        parse_flags::RAW_TEXT | parse_flags::VALID_LOCAL,
        &mut on_warning,
    )?;
    if matches!(id, ParsedSeqId::Prf(_) | ParsedSeqId::Pir(_)) {
        id = local_from_label(&seq_id_string_with_version(&id));
    }
    if of_other_molecule(&id, is_aa, warn) {
        id = local_from_label(&seq_id_string_with_version(&id));
    }
    Ok(vec![id])
}

/// `CBioseq::GetNonLocalId()` of a FASTA Bioseq: its best id by `s_BestNonLocalRank` if
/// that is not local, else the non-local row of its Seq-hist alignment (an id with a
/// range).
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/Bioseq.cpp:280-323
/// The alignment of `CFastaReader::xCreateAlignment` has the generated `lcl|N` as its
/// first row and the replaced id as its second; the generated id is the Bioseq's only id
/// (the edits of `AddSequences` keep local ids), so the second row is returned when it is
/// not local.
fn non_local_id(parsed: &ParsedIds) -> Option<&ParsedSeqId> {
    let score = |id: &ParsedSeqId| {
        if matches!(id, ParsedSeqId::Local(_)) {
            i32::MAX - 1
        } else {
            best_rank_score(id)
        }
    };
    let best = find_best(&parsed.ids, score)?;
    if !matches!(best, ParsedSeqId::Local(_)) {
        return Some(best);
    }
    let old = parsed.hist_old_id.as_ref()?;
    let generated = parsed.ids.first()?;
    (generated == best && !matches!(old, ParsedSeqId::Local(_))).then_some(old)
}

/// `s_UseFastaReaderDeflines(bioseq, deflines, long_seqid = false)`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:736-762
fn use_fasta_reader_deflines(parsed: &ParsedIds, deflines: &DeflineSet) -> bool {
    match non_local_id(parsed) {
        None | Some(ParsedSeqId::Prf(_)) | Some(ParsedSeqId::Pir(_)) => true,
        Some(_) => deflines.deflines.iter().any(|defline| {
            find_best(&defline.ids, best_rank_score)
                .is_some_and(|id| !matches!(id, ParsedSeqId::Local(_)))
        }),
    }
}

/// `CWriteDB::ExtractBioseqDeflines(bioseq, parse_ids = true, long_ids = false, scan =
/// true)` for a Bioseq of the FASTA path (after the edits of `AddSequences`), without the
/// serialization (the writer posts its diagnostics).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:802-836
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:606-614,648,658-662,694,715-719
/// (`x_BuildDeflinesFromBioseq`: the one definition line of the Bioseq's ids and title
/// descriptor; a FASTA Bioseq has no taxon.) The Bioseq always has an id, so the "No
/// deflines" error does not occur.
pub fn extract_bioseq_deflines(bioseq: &Bioseq, warn: Warn<'_>) -> DeflineSet {
    let parsed = bioseq
        .parsed
        .as_ref()
        .expect("a Bioseq read with -parse_seqids");
    if let Some(deflines) = second_parse(parsed, bioseq.is_aa, warn) {
        if use_fasta_reader_deflines(parsed, &deflines) {
            return deflines;
        }
    }
    deflines_of_bioseq(
        parsed.ids.clone(),
        parsed.title_descr.clone().unwrap_or_default(),
        &[],
    )
    .expect("a FASTA Bioseq has an id")
}

/// `x_BuildDeflinesFromBioseq` of a Bioseq of the ids `ids` whose descriptors give the title
/// `title` (the last title descriptor; empty without one) and the taxids `taxids`
/// ([`descr_taxids`]): one definition line of all the ids, the title and the first taxid;
/// `None` for a Bioseq without ids, which leaves NCBI without definition lines ("Error: No
/// deflines provided.", `x_ExtractDeflines`). The FASTA path ([`extract_bioseq_deflines`]) and
/// the Bioseqs of PSSMs (makeprofiledb, [`bioseq_deflines`]) build theirs here.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:592-703
pub(crate) fn deflines_of_bioseq(
    ids: Vec<ParsedSeqId>,
    title: Vec<u8>,
    taxids: &[i32],
) -> Option<DeflineSet> {
    if ids.is_empty() {
        return None;
    }
    Some(DeflineSet {
        deflines: vec![BlastDefline {
            title: Some(title),
            ids,
            taxid: taxids.first().copied(),
            ..BlastDefline::default()
        }],
    })
}

/// Why `CWriteDB::ExtractBioseqDeflines` makes no definition lines of a Bioseq.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum BioseqDeflinesError {
    /// The Bioseq has no ids: NCBI's `CWriteDBException` "Error: No deflines provided.".
    NoDeflines,
    /// A user object `ASN1_BlastDefLine` holds the binary definition lines
    /// (`x_GetBioseqBinaryHeader`), which NCBI decodes from ASN.1 binary (BER): LOSAT does not
    /// read them yet.
    BinaryHeader,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1574-1596
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:777-836
// With `scan_bioseq_4_cfastareader_usrobj` false `x_GetFastaReaderDeflines` finds no definition
// line, so the definition lines are the binary ones of the descriptors or
// `x_BuildDeflinesFromBioseq`'s.
/// `CWriteDB::ExtractBioseqDeflines(bioseq, parse_ids, long_ids = false, scan = false)` of a
/// Bioseq with descriptors (`descr`) and the ids `ids`, without the serialization (the writer
/// posts its diagnostics: `WriteDb::post_extracted_deflines`).
pub(crate) fn bioseq_deflines(
    ids: &[ParsedSeqId],
    descr: &[Seqdesc],
) -> Result<DeflineSet, BioseqDeflinesError> {
    if has_binary_header(descr) {
        return Err(BioseqDeflinesError::BinaryHeader);
    }
    let title = descr
        .iter()
        .rev()
        .find_map(|descriptor| match descriptor {
            Seqdesc::Title(title) => Some(title.clone()),
            Seqdesc::Other(_) => None,
        })
        .unwrap_or_default();
    deflines_of_bioseq(ids.to_vec(), title, &descr_taxids(descr))
        .ok_or(BioseqDeflinesError::NoDeflines)
}

/// The variant `name` of the `Seqdesc` CHOICE `descriptor`, when it is that variant.
fn seqdesc_variant<'a>(descriptor: &'a Seqdesc, name: &str) -> Option<&'a Value> {
    match descriptor {
        Seqdesc::Other(Value::Choice(index, value))
            if usize::from(*index) == member_index(named::SEQDESC, name) =>
        {
            Some(value)
        }
        _ => None,
    }
}

/// The `Object-id` `value`'s `id`, when it is that variant.
fn object_id_number(value: &Value) -> Option<i64> {
    match value {
        Value::Choice(index, id) if usize::from(*index) == member_index(named::OBJECT_ID, "id") => {
            match id.as_ref() {
                Value::Int(number) => Some(*number),
                _ => None,
            }
        }
        _ => None,
    }
}

/// The `Object-id` `value`'s `str`, when it is that variant.
fn object_id_str(value: &Value) -> Option<&[u8]> {
    match value {
        Value::Choice(index, text)
            if usize::from(*index) == member_index(named::OBJECT_ID, "str") =>
        {
            match text.as_ref() {
                Value::Str(text) => Some(text),
                _ => None,
            }
        }
        _ => None,
    }
}

/// The taxids of the `taxon` Dbtags of the `Org-ref` value `org` whose tag is an `id`, in order.
fn org_taxids(org: &Value) -> impl Iterator<Item = i32> + '_ {
    let tags = match org.member(member_index(named::ORG_REF, "db")) {
        Some(Value::List(tags)) => tags.as_slice(),
        _ => &[],
    };
    tags.iter().filter_map(|tag| {
        let db = tag.member(member_index(named::DBTAG, "db"));
        if !matches!(db, Some(Value::Str(db)) if db.as_slice() == b"taxon") {
            return None;
        }
        // `CObject_id::TId` is an `int`; the reader keeps the INTEGER of Int4 storage.
        object_id_number(tag.member(member_index(named::DBTAG, "tag"))?).map(|id| id as i32)
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:608-640
/// The taxids that `x_BuildDeflinesFromBioseq` collects from `descr`: those of the `taxon` tags
/// of each `source` descriptor's organism and each `org` descriptor, in order.
pub(crate) fn descr_taxids(descr: &[Seqdesc]) -> Vec<i32> {
    let mut taxids = Vec::new();
    for descriptor in descr {
        let org = seqdesc_variant(descriptor, "org").or_else(|| {
            seqdesc_variant(descriptor, "source")?.member(member_index(named::BIOSOURCE, "org"))
        });
        if let Some(org) = org {
            taxids.extend(org_taxids(org));
        }
    }
    taxids
}

// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Org_ref.cpp:71-90
/// `GetOrg().GetTaxId()` of the first `org` descriptor of `descr` (the first `taxon` tag whose tag
/// is an `id`, 0 without one); `None` without an `org` descriptor.
pub(crate) fn first_org_taxid(descr: &[Seqdesc]) -> Option<i32> {
    let org = descr
        .iter()
        .find_map(|descriptor| seqdesc_variant(descriptor, "org"))?;
    Some(org_taxids(org).next().unwrap_or(0))
}

/// `kAsnDeflineObjLabel`.
const ASN_DEFLINE_OBJ_LABEL: &[u8] = b"ASN1_BlastDefLine";

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:441-489
/// Whether `x_GetBioseqBinaryHeader` finds binary definition lines in `descr`: the first user
/// object of type `ASN1_BlastDefLine` whose first field has that label and octet strings
/// (`oss`), the first of which is not empty.
fn has_binary_header(descr: &[Seqdesc]) -> bool {
    let field_data = members(named::USER_FIELD)[member_index(named::USER_FIELD, "data")].ty;
    let oss = member_index(field_data, "oss");
    for descriptor in descr {
        let Some(object) = seqdesc_variant(descriptor, "user") else {
            continue;
        };
        let is_label = |value: Option<&Value>| {
            value.and_then(object_id_str) == Some(ASN_DEFLINE_OBJ_LABEL)
        };
        if !is_label(object.member(member_index(named::USER_OBJECT, "type"))) {
            continue;
        }
        let Some(Value::List(fields)) = object.member(member_index(named::USER_OBJECT, "data"))
        else {
            continue;
        };
        let Some(field) = fields.first() else {
            continue;
        };
        if !is_label(field.member(member_index(named::USER_FIELD, "label"))) {
            continue;
        }
        if let Some(Value::Choice(index, data)) =
            field.member(member_index(named::USER_FIELD, "data"))
        {
            if usize::from(*index) == oss {
                return match data.as_ref() {
                    Value::List(strings) => {
                        matches!(strings.first(), Some(Value::Octets(bytes)) if !bytes.is_empty())
                    }
                    _ => false,
                };
            }
        }
    }
    false
}

#[cfg(test)]
#[path = "deflines_tests.rs"]
mod tests;
