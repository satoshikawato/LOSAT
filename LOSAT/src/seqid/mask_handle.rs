//! When two Seq-ids are one `CSeq_id_Handle` in NCBI BLAST+ 2.17.0: the identity that
//! `makeblastdb`'s `CMaskedRangeSet` matches a mask's Seq-id and a database sequence's Seq-ids
//! by (`c++/src/app/blastdb/masked_range_set.hpp`, `.cpp`: a `map<CSeq_id_Handle, CConstRef<CSeq_loc>>`
//! per algorithm, filled with `CSeq_id_Handle::GetHandle` of the id and looked up the same way), and
//! that the object manager keys the ranges of a PSSM query's annotations by
//! (`algorithm::pssm_search::annot_index`, `CHandleRangeMap`; psiblast and tblastn `-in_pssm`
//! probes of line RG D3-1, `om/` cases `s_i*` and `s_o*`, agree with the rules below and with
//! [`mask_handle_order`]).
//!
//! `CSeq_id_Handle::GetHandle` is `CSeq_id_Mapper::GetHandle` = the `FindOrCreate` of the
//! alternative's tree (seq_id_handle.cpp:113-116, seq_id_mapper.cpp:121-125); two handles are
//! equal when `m_Packed` and `m_Info` are (seq_id_handle.hpp:193-200), so two Seq-ids are one
//! handle exactly when their tree gives them the same info and packed number. The rules of every
//! alternative, with the probes of NCBI 2.17.0 that settled them (98 cases, ASCII case only):
//!
//! - Case never matters in accessions, names, local string ids, `general` databases and string
//!   tags, patent country and number. It does in a Textseq-id's `release`, a PDB `chain-id`, and
//!   in PDB molecule names except their 2nd to 4th characters.
//! - What is present matters: a version, name or release set on one side only gives two handles.
//! - GenBank, EMBL and DDBJ are one tree. A Textseq-id with an accession only (no name, no
//!   release) whose accession packs (`CSeq_id_Textseq_Info::ParseAcc`, seq_id_tree.cpp:824-880)
//!   is one handle across the three; any other Textseq-id needs the same alternative too.
//! - A string never equals an integer: `local str "5"` is not `local id 5`, `general` tag `str
//!   "123"` is not tag `id 123` (the matching of `FindMatch` is not used).
//! - `swissprot` compares its `release` exactly (`sp|` reads as `"reviewed"`, `tr|` as
//!   `"unreviewed"`), `patent` country, number kind, number and sequence number, `pdb` the ids
//!   after `CSeq_id_PDB_Info::Normalize`, `gi`, `gibbsq` and `gibbmt` their number.
//! - A Textseq-id with neither accession nor name gets a new handle on every call: it matches
//!   nothing, itself included.
//!
//! Two things a pure function cannot give: the handle of a `giim` id ([`TreeHandle`]: the Giim
//! tree's `x_FindInfo` takes the first stored handle of the number whose Giimport-id is *not*
//! `Equals`, seq_id_tree.cpp:2926-2941, so the handle depends on the ids looked up before; the
//! caller keeps that tree, `writedb::mask_data::range_set`), and the order of the ranges that come
//! from two different handles with the same packed number ([`mask_handle_packed`]: their order is
//! the info's address; the caller refuses it where the two orders write different ranges).
//!
//! NCBI's `OBJECTS_PACK_TEXTID` and `OBJECTS_PACK_GENERAL` parameters (environment, or
//! `[OBJECTS] PACK_TEXTID` in the registry; seq_id_tree.cpp:692-707) turn packing off, which would
//! change the identity of EMBL/DDBJ ids and the order of ranges; LOSAT refuses them
//! (`blastinput::ncbi_environment`), so packing is always on here.
//!
//! The probes are in `tests/seqdb/mask_ids/` (`README.md` says how they were made).

use std::collections::HashMap;

use super::model::{
    Date, DateStd, GiimportId, ObjectId, ParsedSeqId, PatentNumber, PdbSeqId, TextseqId,
};

/// A Seq-id whose handle a pure function cannot know: it is found in a tree whose state the
/// caller keeps.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum TreeHandle {
    /// A `giim` id: `CSeq_id_Giim_Tree::x_FindInfo` returns the first stored handle of the same
    /// `id` whose Giimport-id is NOT `Equals` (seq_id_tree.cpp:2926-2941), so the handle depends
    /// on the giim ids of that number made before it and their order. Ids of other numbers never
    /// share its handle.
    Giim { id: i32 },
}

// NCBI reference (598d8ae6): c++/include/objects/seq/seq_id_handle.hpp:193-207
/// The key of the `CSeq_id_Handle` that `CSeq_id_Handle::GetHandle` of the id gives in NCBI 2.17.0
/// (`CSeq_id_Mapper::GetHandle` -> `FindOrCreate`): two Seq-ids are one handle exactly when
/// their keys are equal. `Ok(None)`: a Textseq-id with neither accession nor name, which gets a
/// new handle on every call (seq_id_tree.cpp:1187-1198) and so is one handle with no Seq-id,
/// itself included.
pub(crate) fn mask_handle_key(id: &ParsedSeqId) -> Result<Option<Vec<u8>>, TreeHandle> {
    let mut key = vec![handle_tree(id)];
    match id {
        // CSeq_id_Local_Tree: strings without case, integers exactly, a string never an integer
        // (seq_id_tree.cpp:1998-2011, seq_id_tree.hpp:767-768).
        ParsedSeqId::Local(tag) => push_object_id(&mut key, tag),
        // CSeq_id_int_Tree: the integer (seq_id_tree.cpp:412-424, 502-527).
        ParsedSeqId::Gibbsq(number) | ParsedSeqId::Gibbmt(number) => {
            key.extend(number.to_le_bytes());
        }
        // CSeq_id_Gi_Tree: the GI; GI 0 has an info of its own (seq_id_tree.cpp:609-627).
        ParsedSeqId::Gi(gi) => key.extend(gi.to_le_bytes()),
        ParsedSeqId::Giim(giim) => return Err(TreeHandle::Giim { id: giim.id }),
        // CSeq_id_Patent_Tree: country and number (or application number) without case, the
        // kind of number, the sequence number; doc-type is not indexed (seq_id_tree.cpp:3074-3146).
        ParsedSeqId::Patent(patent) => {
            push_bytes(&mut key, &patent.cit.country.to_ascii_lowercase());
            let (kind, number) = match &patent.cit.id {
                PatentNumber::Number(number) => (b'n', number),
                PatentNumber::AppNumber(number) => (b'a', number),
            };
            key.push(kind);
            push_bytes(&mut key, &number.to_ascii_lowercase());
            key.extend(patent.seqid.to_le_bytes());
        }
        // CSeq_id_General_Tree: db without case; a string tag without case, an integer tag exactly,
        // a string never an integer. Packed or plain, the classes are the same
        // (seq_id_tree.cpp:2578-2655, seq_id_tree.hpp:813-818, 918-925).
        ParsedSeqId::General(dbtag) => {
            push_bytes(&mut key, &dbtag.db.to_ascii_lowercase());
            push_object_id(&mut key, &dbtag.tag);
        }
        // CSeq_id_PDB_Tree: `Equals` of the normalized ids (seq_id_tree.cpp:3408-3464, 3502-3523).
        ParsedSeqId::Pdb(pdb) => push_pdb_normalized(&mut key, pdb),
        _ => {
            let Some(text) = id.textseq() else {
                unreachable!("the other alternatives are Textseq-ids")
            };
            return Ok(textseq_key(id, text, key));
        }
    }
    Ok(Some(key))
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.hpp:359-388
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:1161-1185
/// The key of the info (`CSeq_id_Info`) of the handle that `CSeq_id_Handle::GetHandle` of the id gives:
/// two Seq-ids whose keys are equal share one info, and `GetSeqId()` of either handle is made from
/// the Seq-id that made it ([`handle_seq_id`]). It is [`mask_handle_key`] but for a Textseq-id that
/// packs: its info is shared by every accession with the same prefix (without case), number of
/// digits and version (the packed tree's `TKey`), whatever the number and, in the GenBank tree,
/// whatever the alternative (GenBank, EMBL or DDBJ).
pub(crate) fn mask_info_key(id: &ParsedSeqId) -> Result<Option<Vec<u8>>, TreeHandle> {
    if let Some(text) = id.textseq().filter(|text| textseq_packs(text)) {
        let accession = text
            .accession
            .as_deref()
            .expect("a packed Textseq-id has an accession");
        let (prefix_len, digits) =
            textseq_packed_split(accession).expect("a packed Textseq-id splits");
        let mut key = vec![handle_tree(id), b'I'];
        push_bytes(&mut key, &accession[..prefix_len].to_ascii_lowercase());
        key.extend((digits as u32).to_le_bytes());
        push_opt_i32(&mut key, text.version);
        return Ok(Some(key));
    }
    mask_handle_key(id)
}

/// The key of a Textseq-id (`key` holds its tree).
fn textseq_key(id: &ParsedSeqId, text: &TextseqId, mut key: Vec<u8>) -> Option<Vec<u8>> {
    if text.accession.is_none() && text.name.is_none() {
        return None;
    }
    // Packed (seq_id_tree.cpp:1165-1185): accession only (name and release unset) and `ParseAcc`
    // gives a key: the accession without case and the version; the alternative is not part of the
    // key, so GenBank, EMBL and DDBJ share it (one tree, 253-257).
    if text.name.is_none() && text.release.is_none() {
        if let Some(accession) = &text.accession {
            if textseq_packed_split(accession).is_some() {
                key.push(b'P');
                push_bytes(&mut key, &accession.to_ascii_lowercase());
                push_opt_i32(&mut key, text.version);
                return Some(key);
            }
        }
    }
    // Plain (`x_FindStrInfo` + `x_Equals`, seq_id_tree.cpp:1064-1131): the alternative, each
    // member's presence, accession and name without case, release and version exactly.
    key.push(b'L');
    key.push(id.alternative());
    push_opt_bytes(
        &mut key,
        text.accession
            .as_ref()
            .map(|a| a.to_ascii_lowercase())
            .as_deref(),
    );
    push_opt_bytes(
        &mut key,
        text.name
            .as_ref()
            .map(|n| n.to_ascii_lowercase())
            .as_deref(),
    );
    push_opt_bytes(&mut key, text.release.as_deref());
    push_opt_i32(&mut key, text.version);
    Some(key)
}

/// `m_Packed` of the handle (0: the id has an info of its own). `CSeq_id_Handle::operator<` orders
/// handles by one less than `m_Packed`, as an unsigned number ([`mask_handle_order`]), then by the
/// info's address; `CSeq_loc::Add` with `fMerge_All` and `fSort` puts the ranges of each handle
/// together in that order (Seq_loc.cpp:4278, 4547-4596, 4641-4692, 5011-5021). Two different
/// handles with the same `m_Packed` (both 0 included) are ordered by address, which no pure
/// function knows.
pub(crate) fn mask_handle_packed(id: &ParsedSeqId) -> i64 {
    match id {
        ParsedSeqId::Gi(gi) => *gi,
        ParsedSeqId::General(dbtag) => match &dbtag.tag {
            // CSeq_id_General_Id_Info::Pack (seq_id_tree.cpp:2156-2169), TPacked = Int8
            ObjectId::Id(tag) => {
                let tag = i64::from(*tag);
                if tag <= 0 {
                    tag - 1
                } else {
                    tag
                }
            }
            // CSeq_id_General_Str_Info::Parse / Pack (seq_id_tree.cpp:2262-2335), at least 3 digits
            ObjectId::Str(text) => match general_str_packed_split(text) {
                Some((prefix_len, digits)) => {
                    let number = parse_digits(&text[prefix_len..prefix_len + digits]);
                    if number <= 0 {
                        number - 1
                    } else {
                        number
                    }
                }
                None => 0,
            },
        },
        _ => match id.textseq() {
            Some(text) if text.name.is_none() && text.release.is_none() => {
                match text
                    .accession
                    .as_ref()
                    .and_then(|acc| textseq_packed_split(acc).map(|split| (acc, split)))
                {
                    Some((acc, (prefix_len, digits))) => {
                        parse_digits(&acc[prefix_len..prefix_len + digits])
                    }
                    None => 0,
                }
            }
            _ => 0,
        },
    }
}

/// The position of a handle in `CSeq_id_Handle::operator<` among handles of different `m_Packed`:
/// one less than `m_Packed`, as an unsigned number (`TUintId`, 64 bits), so an unpacked handle (0)
/// comes last.
pub(crate) fn mask_handle_order(packed: i64) -> u64 {
    (packed as u64).wrapping_sub(1)
}

// NCBI reference (598d8ae6): c++/include/objects/seq/seq_id_handle.hpp:300-313 (the Seq-id of a
// handle: its info's own, or the one its packed number and case variant give back)
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:736-786 (a case variant has
// one bit per letter of the info's Seq-id, at most 64: `TVariant` is `Uint8`)
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:938-966, 1022-1029 (a packed
// Textseq-id is made again with the info's alternative and the accession's case; a plain one is
// the info's Seq-id with the accession's case only), 1934-1950 (a local string id), 2172-2181,
// 2338-2357, 2466-2481 (a general id: the database, then a string tag), 3408-3463 (a PDB id: the
// variant of `CSeq_id_PDB_Info::Normalize` gives back the molecule's case, `chain` and `chain-id`)
/// `CSeq_id_Handle::GetSeqId()` of the handle that `CSeq_id_Handle::GetHandle` of the id gives when the
/// handle's info was made from `first` (the two share the info: [`mask_info_key`], or the `giim`
/// tree's [`GiimTree::made_from`]): the Seq-id that `AsString()` writes and the masking tools'
/// `seqloc_*` and `maskinfo_*` writers print. The info keeps `first`;
/// the handle's case variant gives back the case of `id`'s first 64 letters where the tree keeps
/// one:
/// - a local string id; a general id's database, then its string tag;
/// - a Textseq-id's accession only (a name keeps `first`'s case); a packed one (accession only)
///   takes `first`'s alternative (GenBank, EMBL and DDBJ are one tree);
/// - a PDB id as a whole (molecule case, `chain-id`, release date) but for its `chain`, the
///   normalized one ([`pdb_handle_seq_id`]);
/// - nothing of a patent or a `giim` id (`first` as it is); gi, gibbsq, gibbmt and local integer
///   ids are their number.
pub(crate) fn handle_seq_id(first: &ParsedSeqId, id: &ParsedSeqId) -> ParsedSeqId {
    match (first, id) {
        (ParsedSeqId::Local(ObjectId::Str(text)), ParsedSeqId::Local(ObjectId::Str(own))) => {
            let mut text = text.clone();
            restore_case(&mut [&mut text], &[own]);
            ParsedSeqId::Local(ObjectId::Str(text))
        }
        (ParsedSeqId::General(dbtag), ParsedSeqId::General(own)) => {
            let mut dbtag = dbtag.clone();
            match (&mut dbtag.tag, &own.tag) {
                (ObjectId::Str(text), ObjectId::Str(own_text)) => {
                    restore_case(&mut [&mut dbtag.db, text], &[&own.db, own_text]);
                }
                _ => restore_case(&mut [&mut dbtag.db], &[&own.db]),
            }
            ParsedSeqId::General(dbtag)
        }
        (ParsedSeqId::Patent(_) | ParsedSeqId::Giim(_), _) => first.clone(),
        (_, ParsedSeqId::Pdb(own)) => ParsedSeqId::Pdb(pdb_handle_seq_id(own)),
        _ => match (first.textseq(), id.textseq()) {
            (Some(text), Some(own)) => {
                let mut out = first.clone();
                let target = out.textseq_mut().expect("a Textseq-id");
                if textseq_packs(own) {
                    *target = own.clone();
                } else if let (Some(accession), Some(own)) =
                    (target.accession.as_mut(), own.accession.as_ref())
                {
                    debug_assert_eq!(text.accession.as_ref().map(Vec::len), Some(own.len()));
                    restore_case(&mut [accession], &[own]);
                }
                out
            }
            _ => id.clone(),
        },
    }
}

/// Whether `AsString()` of `id`'s handle writes `id`'s FASTA text whatever Seq-id of the handle
/// made its info (`handle_seq_id`): false where the info's Seq-id shows through, namely a GenBank,
/// EMBL or DDBJ accession that packs (the alternative), a Textseq-id name with a letter, a patent
/// country or number with a letter, and more than 64 letters in a local string id, a general id's
/// database and string tag, or a Textseq-id's accession (the case variant's 64 bits). A `giim` id's
/// text is its number.
pub(crate) fn handle_fasta_is_own(id: &ParsedSeqId) -> bool {
    let letters = |text: &[u8]| text.iter().filter(|c| c.is_ascii_alphabetic()).count();
    match id {
        ParsedSeqId::Local(ObjectId::Str(text)) => letters(text) <= VARIANT_LETTERS,
        ParsedSeqId::General(dbtag) => {
            let tag = match &dbtag.tag {
                ObjectId::Str(text) => letters(text),
                ObjectId::Id(_) => 0,
            };
            letters(&dbtag.db) + tag <= VARIANT_LETTERS
        }
        ParsedSeqId::Patent(patent) => {
            let number = match &patent.cit.id {
                PatentNumber::Number(number) | PatentNumber::AppNumber(number) => number,
            };
            letters(&patent.cit.country) + letters(number) == 0
        }
        _ => match id.textseq() {
            Some(text) if textseq_packs(text) => !matches!(
                id,
                ParsedSeqId::Genbank(_) | ParsedSeqId::Embl(_) | ParsedSeqId::Ddbj(_)
            ),
            Some(text) => {
                text.name.as_deref().is_none_or(|name| letters(name) == 0)
                    && text
                        .accession
                        .as_deref()
                        .is_none_or(|acc| letters(acc) <= VARIANT_LETTERS)
            }
            None => true,
        },
    }
}

/// The letters a case variant covers (`CSeq_id_Handle::TVariant` = `Uint8`, one bit each).
const VARIANT_LETTERS: usize = 64;

/// `s_RestoreCaseVariant` of `s_ParseCaseVariant`: the first 64 letters of `to` (its strings in
/// order) take the case of the same letters of `from` (equal to `to` but for case).
fn restore_case(to: &mut [&mut Vec<u8>], from: &[&Vec<u8>]) {
    let mut left = VARIANT_LETTERS;
    for (to, from) in to.iter_mut().zip(from) {
        for (c, &own) in to.iter_mut().zip(from.iter()) {
            if left == 0 {
                return;
            }
            if c.is_ascii_alphabetic() {
                *c = own;
                left -= 1;
            }
        }
    }
}

/// Whether a Textseq-id packs (`FindOrCreate`, seq_id_tree.cpp:1165-1185): an accession only that
/// `ParseAcc` takes.
fn textseq_packs(text: &TextseqId) -> bool {
    text.name.is_none()
        && text.release.is_none()
        && text
            .accession
            .as_deref()
            .is_some_and(|acc| textseq_packed_split(acc).is_some())
}

/// The handle's tree (`CSeq_id_Which_Tree::Initialize`, seq_id_tree.cpp:244-275): the alternative,
/// with GenBank, EMBL and DDBJ in the GenBank tree.
fn handle_tree(id: &ParsedSeqId) -> u8 {
    match id {
        ParsedSeqId::Embl(_) | ParsedSeqId::Ddbj(_) => {
            ParsedSeqId::Genbank(TextseqId::default()).alternative()
        }
        _ => id.alternative(),
    }
}

/// `CSeq_id_Textseq_Info::ParseAcc` (seq_id_tree.cpp:824-880): `Some((prefix length, digits))`
/// when the accession packs: the trailing run of digits must hold a non-zero digit, have 2-12
/// digits, at most 9 from its first non-zero digit, and the prefix before it must be at most 4
/// characters and at most twice the digit count (the `prefix_len == 3` branch at 851 sits after
/// `prefix_len <= 4` and never runs). With more than 6 digits and leading zeros, the zeros
/// beyond 6 move into the prefix, which must then be at most 7 characters.
fn textseq_packed_split(acc: &[u8]) -> Option<(usize, usize)> {
    let len = acc.len();
    let mut prefix_len = len;
    let mut most_significant = None;
    while prefix_len > 0 {
        prefix_len -= 1;
        let c = acc[prefix_len];
        if (b'1'..=b'9').contains(&c) {
            most_significant = Some(prefix_len);
        } else if c != b'0' {
            prefix_len += 1;
            break;
        }
    }
    let most_significant = most_significant?;
    let mut acc_digits = len - prefix_len;
    let real_digits = len - most_significant;
    if !(2..=12).contains(&acc_digits)
        || real_digits > 9
        || acc_digits * 2 < prefix_len
        || prefix_len > 4
    {
        return None;
    }
    if acc_digits > 6 && real_digits < acc_digits {
        acc_digits = real_digits.max(6);
        prefix_len = len - acc_digits;
    }
    (prefix_len <= 7).then_some((prefix_len, acc_digits))
}

/// `CSeq_id_General_Str_Info::Parse` (seq_id_tree.cpp:2262-2320) with `kMinGeneralStrDigits`
/// (2530): `Some((prefix length, digits))` when the string tag packs.
fn general_str_packed_split(text: &[u8]) -> Option<(usize, usize)> {
    let len = text.len();
    let (mut prefix_len, mut str_digits) = (len, 0usize);
    let (mut cur_digits, mut total_digits) = (0usize, 0usize);
    let mut i = len as isize;
    while i >= 0 {
        i -= 1;
        let c = if i < 0 { 0 } else { text[i as usize] };
        if c.is_ascii_digit() {
            total_digits += 1;
            cur_digits += 1;
        } else {
            if str_digits == 0 || cur_digits > str_digits + 2 {
                str_digits = cur_digits;
                prefix_len = (i + 1) as usize;
            }
            cur_digits = 0;
        }
    }
    if str_digits > 9 {
        prefix_len += str_digits - 9;
        total_digits += str_digits - 9;
        str_digits = 9;
    }
    if str_digits * 3 < total_digits * 2 || str_digits < 3 {
        return None;
    }
    Some((prefix_len, str_digits))
}

fn parse_digits(digits: &[u8]) -> i64 {
    digits
        .iter()
        .fold(0i64, |n, &c| n * 10 + i64::from(c - b'0'))
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:3408-3464
/// The PDB id that `CSeq_id_PDB_Info::Normalize` keys the tree with:
/// `None` for an id with a release date that does not normalize, which is kept as it is; any other
/// id loses its date, gets the molecule's 2nd-4th characters in upper case when the molecule has a
/// lower-case letter, `chain_id` from `chain` (as a `char`) and `chain` from `chain_id[0]` (a
/// signed char) when one of them is set (`chain` only for a 1-character `chain_id`). The tree
/// compares the results with `Equals` (values and set flags; 3502-3523).
fn pdb_normalized(pdb: &PdbSeqId) -> Option<PdbSeqId> {
    let date_packs = match &pdb.rel {
        None => true,
        Some(Date::Std(date)) => pdb_date_normalizes(date),
        Some(Date::Str(_)) => false,
    };
    if !date_packs {
        return None;
    }
    let has_chain = pdb.chain.is_some() || pdb.chain_id.as_ref().is_some_and(|c| c.len() == 1);
    let has_chain_id = pdb.chain_id.is_some() || pdb.chain.is_some();
    let need_upcase = pdb.mol.iter().any(u8::is_ascii_lowercase);
    if pdb.rel.is_none()
        && !need_upcase
        && pdb.chain.is_some() == has_chain
        && pdb.chain_id.is_some() == has_chain_id
    {
        return Some(pdb.clone());
    }
    let mut mol = pdb.mol.clone();
    if need_upcase {
        mol.iter_mut()
            .skip(1)
            .take(3)
            .for_each(u8::make_ascii_uppercase);
    }
    let chain_id = has_chain_id.then(|| {
        pdb.chain_id
            .clone()
            .unwrap_or_else(|| vec![pdb.chain.unwrap_or(0) as u8])
    });
    let chain = has_chain.then(|| {
        chain_id
            .as_ref()
            .and_then(|c| c.first())
            .map_or(0, |&c| i32::from(c as i8))
    });
    Some(PdbSeqId {
        mol,
        chain,
        rel: None,
        chain_id,
    })
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:3293-3346
/// The Seq-id of the handle of the PDB id `pdb` (`CSeq_id_Handle::GetSeqId()`, which
/// `GetPackedSeqId` makes from the normalized id, [`pdb_normalized`]): the members that `pdb`
/// lacks are unset again and its molecule's case and release date come back, so the handle's
/// Seq-id is `pdb` but for `chain`, the normalized one where `pdb` has a `chain`: a `chain` outside
/// -128..=127 without a `chain-id` (the chain's low byte read back as a signed char), or a
/// `chain` that is not `chain-id[0]` of an id that `Normalize` builds anew (a lower-case molecule
/// letter, or a release date). The object manager looks the handle's Seq-id up, so a database whose stored id is
/// `pdb` does not serve it where that differs (`ParsedSubject::scope_finds_none`).
pub(crate) fn pdb_handle_seq_id(pdb: &PdbSeqId) -> PdbSeqId {
    let mut handle = pdb.clone();
    if let (Some(chain), Some(normalized)) = (handle.chain.as_mut(), pdb_normalized(pdb)) {
        *chain = normalized.chain.unwrap_or(*chain);
    }
    handle
}

/// The key of a PDB id: the normalized id ([`pdb_normalized`]), or the id as it is (with its
/// release date) where the date does not normalize.
fn push_pdb_normalized(key: &mut Vec<u8>, pdb: &PdbSeqId) {
    let Some(normalized) = pdb_normalized(pdb) else {
        key.push(b'R');
        push_bytes(key, &pdb.mol);
        push_opt_i32(key, pdb.chain);
        if let Some(rel) = &pdb.rel {
            push_date(key, rel);
        }
        push_opt_bytes(key, pdb.chain_id.as_deref());
        return;
    };
    key.push(b'N');
    push_bytes(key, &normalized.mol);
    push_opt_i32(key, normalized.chain);
    push_opt_bytes(key, normalized.chain_id.as_deref());
}

/// `CSeq_id_PDB_Info::x_NormalizeDate` (seq_id_tree.cpp:3381-3405): year 1-4095, month 1-15,
/// day 1-31, hour 0-30, minute and second 0-62, no season.
fn pdb_date_normalizes(date: &DateStd) -> bool {
    let in_range =
        |value: Option<i32>, low: i32, high: i32| value.is_none_or(|v| (low..=high).contains(&v));
    (1..=4095).contains(&date.year)
        && in_range(date.month, 1, 15)
        && in_range(date.day, 1, 31)
        && date.season.is_none()
        && in_range(date.hour, 0, 30)
        && in_range(date.minute, 0, 62)
        && in_range(date.second, 0, 62)
}

fn push_object_id(key: &mut Vec<u8>, tag: &ObjectId) {
    match tag {
        ObjectId::Id(id) => {
            key.push(b'i');
            key.extend(id.to_le_bytes());
        }
        ObjectId::Str(text) => {
            key.push(b's');
            push_bytes(key, &text.to_ascii_lowercase());
        }
    }
}

fn push_date(key: &mut Vec<u8>, date: &Date) {
    match date {
        Date::Str(text) => {
            key.push(b's');
            push_bytes(key, text);
        }
        Date::Std(d) => {
            key.push(b'd');
            key.extend(d.year.to_le_bytes());
            push_opt_i32(key, d.month);
            push_opt_i32(key, d.day);
            push_opt_bytes(key, d.season.as_deref());
            push_opt_i32(key, d.hour);
            push_opt_i32(key, d.minute);
            push_opt_i32(key, d.second);
        }
    }
}

fn push_bytes(key: &mut Vec<u8>, bytes: &[u8]) {
    key.extend((bytes.len() as u32).to_le_bytes());
    key.extend_from_slice(bytes);
}

fn push_opt_bytes(key: &mut Vec<u8>, bytes: Option<&[u8]>) {
    match bytes {
        None => key.push(0),
        Some(bytes) => {
            key.push(1);
            push_bytes(key, bytes);
        }
    }
}

fn push_opt_i32(key: &mut Vec<u8>, value: Option<i32>) {
    match value {
        None => key.push(0),
        Some(value) => {
            key.push(1);
            key.extend(value.to_le_bytes());
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:2926-2964
/// NCBI's `CSeq_id_Giim_Tree` ([`TreeHandle::Giim`]): the live handles of each `giim` number, in
/// the order they were made, each with the Giimport-id it was made from (`GetSeqId()` of the
/// info). The one model of the tree: `makeblastdb`'s mask ranges (`writedb::mask_data::range_set`)
/// and the object manager's index and attach of a PSSM query (`algorithm::pssm_search::annot_index`)
/// keep their handles in it and drop each one when the last thing holding it goes.
#[derive(Clone, Debug, Default)]
pub(crate) struct GiimTree {
    lists: HashMap<i32, Vec<(usize, GiimportId)>>,
    made: usize,
}

impl GiimTree {
    /// `x_FindInfo`: the first live handle whose Giimport-id is not `Equals` to `id`.
    pub(crate) fn find(&self, id: &GiimportId) -> Option<usize> {
        self.lists
            .get(&id.id)?
            .iter()
            .find(|(_, made_from)| made_from != id)
            .map(|(handle, _)| *handle)
    }

    /// `FindOrCreate`: the handle, and whether it is new.
    pub(crate) fn find_or_create(&mut self, id: &GiimportId) -> (usize, bool) {
        if let Some(handle) = self.find(id) {
            return (handle, false);
        }
        self.made += 1;
        self.lists
            .entry(id.id)
            .or_default()
            .push((self.made, id.clone()));
        (self.made, true)
    }

    /// The Giimport-id the handle was made from.
    pub(crate) fn made_from(&self, number: i32, handle: usize) -> &GiimportId {
        let list = &self.lists[&number];
        &list
            .iter()
            .find(|(h, _)| *h == handle)
            .expect("a live handle")
            .1
    }

    /// `x_Unindex` of a handle that nothing holds any more.
    pub(crate) fn drop_handle(&mut self, number: i32, handle: usize) {
        if let Some(list) = self.lists.get_mut(&number) {
            list.retain(|(h, _)| *h != handle);
            if list.is_empty() {
                self.lists.remove(&number);
            }
        }
    }
}

#[cfg(test)]
#[path = "mask_handle_tests.rs"]
mod tests;
