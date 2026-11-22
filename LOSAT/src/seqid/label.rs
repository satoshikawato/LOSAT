//! The printings of a Seq-id that NCBI's formatters use: the FASTA form (`AsFastaString`,
//! `ref|XP_000001.1|`), the content label (`GetSeqIdString`, `XP_000001.1`), the formatter's
//! label (`CAlignFormatUtil::GetLabel`: a `general` id as its tag), the bare id
//! (`CAlignFormatUtil::GetBareId`), the id list of a definition line
//! (`CShowBlastDefline::GetSeqIdListString`) and the first GI (`FindGi`).
//!
//! Each function appends to `out`; the strings are the bytes the reader gives (`seqid::ber`).

use std::io::Write;

use super::model::{ObjectId, ParsedSeqId, PatentNumber, PatentSeqId, PdbSeqId, TextseqId};
use super::rank::{find_best, text_score};

/// Appends `value` in decimal (`NStr::NumericToString`, `operator<<`).
pub(super) fn number(out: &mut Vec<u8>, value: impl std::fmt::Display) {
    let _ = write!(out, "{value}");
}

// NCBI reference (598d8ae6): c++/src/objects/general/Object_id.cpp:202-215
fn object_id(tag: &ObjectId, out: &mut Vec<u8>) {
    match tag {
        ObjectId::Id(id) => number(out, id),
        ObjectId::Str(text) => out.extend_from_slice(text),
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:572-596
/// The FASTA type tag of each alternative, in the order of `model::SEQ_ID_ALTERNATIVES`.
const FASTA_TAGS: [&[u8]; 20] = [
    b"lcl", b"bbs", b"bbm", b"gim", b"gb", b"emb", b"pir", b"sp", b"pat", b"ref", b"gnl", b"gi",
    b"dbj", b"prf", b"pdb", b"tpg", b"tpe", b"tpd", b"gpp", b"nat",
];

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2164-2186
// The tag of a patent id with an application number is `pat|` too (the constant named
// `pgp` holds "pat|").
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2347-2375
/// `CSeq_id::WriteAsFasta` of the Bioseq (`GetStringDescr` with `eFormat_FastA`): the first
/// GI of the Bioseq's ids `ids`, then `|` and its best id for FASTA (`FastaAARank` for a
/// protein, else `FastaNARank`) when that is not a GI: `gi|123|gb|ABC12345.1|`.
pub fn bioseq_fasta(ids: &[ParsedSeqId], protein: bool, out: &mut Vec<u8>) {
    let best = if protein {
        find_best(ids, super::rank::fasta_aa_score)
    } else {
        find_best(ids, super::rank::fasta_na_score)
    };
    let gi = ids.iter().find(|id| id.gi().is_some());
    if let Some(gi) = gi {
        fasta(gi, out);
    }
    if let Some(best) = best.filter(|best| best.gi().is_none()) {
        if gi.is_some() {
            out.push(b'|');
        }
        fasta(best, out);
    }
}

/// `CSeq_id::AsFastaString` (`WriteAsFasta`): `ref|XP_000001.1|`, `tr|A0A000|A0A000_YEAST`,
/// `pdb|1ABC|A`, `gnl|mydb|prot10`, `gi|1414`, ...
pub fn fasta(id: &ParsedSeqId, out: &mut Vec<u8>) {
    match id {
        ParsedSeqId::Swissprot(text) if text.release.as_deref() == Some(b"unreviewed") => {
            out.extend_from_slice(b"tr|");
        }
        _ => {
            out.extend_from_slice(FASTA_TAGS[usize::from(id.alternative())]);
            out.push(b'|');
        }
    }
    fasta_content(id, out);
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2188-2262
/// `CSeq_id::x_WriteContentAsFasta`: the FASTA form without the type tag.
fn fasta_content(id: &ParsedSeqId, out: &mut Vec<u8>) {
    match id {
        ParsedSeqId::Local(tag) => object_id(tag, out),
        ParsedSeqId::Gibbsq(value) | ParsedSeqId::Gibbmt(value) => number(out, value),
        ParsedSeqId::Giim(giim) => number(out, giim.id),
        ParsedSeqId::Patent(patent) => patent_fasta(patent, out),
        ParsedSeqId::General(dbtag) => {
            out.extend_from_slice(&dbtag.db);
            out.push(b'|');
            object_id(&dbtag.tag, out);
        }
        ParsedSeqId::Gi(gi) => number(out, gi),
        ParsedSeqId::Pdb(pdb) => pdb_fasta(pdb, out),
        _ => textseq_fasta(id.textseq().expect("a Textseq-id alternative"), out),
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Textseq_id.cpp:228-247
fn textseq_fasta(id: &TextseqId, out: &mut Vec<u8>) {
    if let Some(accession) = &id.accession {
        out.extend_from_slice(accession);
        if let Some(version) = id.version.filter(|&version| version != 0) {
            out.push(b'.');
            number(out, version);
        }
    }
    out.push(b'|');
    if let Some(name) = &id.name {
        out.extend_from_slice(name);
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/PDB_seq_id.cpp:94-108
fn pdb_fasta(id: &PdbSeqId, out: &mut Vec<u8>) {
    out.extend_from_slice(&id.mol);
    out.push(b'|');
    if let Some(chain_id) = &id.chain_id {
        out.extend_from_slice(chain_id);
        return;
    }
    match id.chain_or_default() as u8 {
        b'|' => {}
        0 => out.push(b' '),
        chain => out.push(chain),
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Patent_seq_id.cpp:71-88
fn patent_fasta(id: &PatentSeqId, out: &mut Vec<u8>) {
    out.extend_from_slice(&id.cit.country);
    out.push(b'|');
    match &id.cit.id {
        PatentNumber::Number(number) => out.extend_from_slice(number),
        PatentNumber::AppNumber(number) => {
            out.extend_from_slice(number);
            if let Some(doc_type) = &id.cit.doc_type {
                out.extend_from_slice(doc_type);
            }
        }
    }
    out.push(b'|');
    number(out, id.seqid);
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1926-2038
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2145-2154
/// `CSeq_id::GetSeqIdString` with `with_version` (`GetLabel` of `eContent` with the flag
/// `fLabel_Version` or none): the upper-cased accession (else the name) of a Textseq-id with
/// `.version` when asked and stored, the number of a `gi`, `db:tag` of a `general` id,
/// `mol_chain` of a `pdb` id, `country number_seqid` of a patent.
pub fn content(id: &ParsedSeqId, with_version: bool, out: &mut Vec<u8>) {
    label_content(id, with_version, false, out);
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/multisource_util.cpp:157-160
/// `GetLabel` of `eContent` with `fLabel_GeneralDbIsContent` (makeblastdb's `GetSeqIdKey`, the keys
/// of `-taxid_map`): `content` without the version, and only the tag of a `general` id (its
/// database is part of the type, not of the content).
pub fn content_general_db_is_content(id: &ParsedSeqId, out: &mut Vec<u8>) {
    label_content(id, false, true, out);
}

/// `x_GetLabel_Content` with the flags `fLabel_Version` (`with_version`) and
/// `fLabel_GeneralDbIsContent` (`general_db_is_content`).
fn label_content(
    id: &ParsedSeqId,
    with_version: bool,
    general_db_is_content: bool,
    out: &mut Vec<u8>,
) {
    if let Some(text) = id.textseq() {
        let start = out.len();
        if let Some(accession) = &text.accession {
            out.extend(accession.iter().map(u8::to_ascii_uppercase));
        } else if let Some(name) = &text.name {
            out.extend_from_slice(name);
        }
        if out.len() > start && with_version {
            if let Some(version) = text.version {
                out.push(b'.');
                number(out, version);
            }
        }
        return;
    }
    match id {
        ParsedSeqId::Local(tag) => object_id(tag, out),
        ParsedSeqId::Gibbsq(value) | ParsedSeqId::Gibbmt(value) => number(out, value),
        ParsedSeqId::Giim(giim) => number(out, giim.id),
        ParsedSeqId::General(dbtag) => {
            if !general_db_is_content {
                out.extend_from_slice(&dbtag.db);
                out.push(b':');
            }
            object_id(&dbtag.tag, out);
        }
        ParsedSeqId::Patent(patent) => {
            out.extend_from_slice(&patent.cit.country);
            match &patent.cit.id {
                PatentNumber::Number(number) | PatentNumber::AppNumber(number) => {
                    out.extend_from_slice(number)
                }
            }
            out.push(b'_');
            number(out, patent.seqid);
        }
        ParsedSeqId::Gi(gi) => number(out, gi),
        ParsedSeqId::Pdb(pdb) => {
            out.extend_from_slice(&pdb.mol);
            if let Some(chain_id) = &pdb.chain_id {
                out.push(b'_');
                out.extend_from_slice(chain_id);
            } else if let Some(chain) = pdb.chain {
                let chain = chain as u8;
                if chain > b' ' {
                    out.push(b'_');
                    out.push(chain);
                }
            }
        }
        _ => unreachable!("Textseq-id alternatives are handled above"),
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1890-1920
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2059-2084
// Without `fLabel_GeneralDbIsContent` the type of a `general` id is "gnl" and of a patent
// "pat", which are the FASTA tags; `s_TextId` names a swissprot id "sp" whatever its release.
/// `GetLabel` of `eDefault` with `fLabel_Default` and `fLabel_UpperCase` (`eDefault` is `eBoth`,
/// `fLabel_Default` is `fLabel_Version`): the label of CWriteDB_Volume's duplicate check,
/// `REF|NP_000549.1`, `SP|P69905`, `LCL|ABC`, `GNL|DB|TAG`, `PDB|1ABC|a`, `GI|123`,
/// `PAT|US|1234567|1`.
pub fn both_upper_version(id: &ParsedSeqId, out: &mut Vec<u8>) {
    let start = out.len();
    out.extend_from_slice(FASTA_TAGS[usize::from(id.alternative())]);
    out.push(b'|');
    out[start..].make_ascii_uppercase();
    compose_oslt_label(id, out);
    if let Some(version) = id.textseq().and_then(|text| text.version) {
        out.push(b'.');
        number(out, version);
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2059-2086 (cited at
// `both_upper_version`): without `fLabel_UpperCase`, `eBoth` is the type, `|` and
// `x_GetLabel_Content` (`label_content`).
/// `GetLabel` of `eDefault` with `fLabel_Default` (`eBoth` with `fLabel_Version`, not upper-cased):
/// `lcl|Query_1`, `ref|NP_000549.1`, `gnl|CDD:1234`, `gi|123` (rpsbproc's query label,
/// `CSeq_id::GetStringDescr` with `eFormat_BestWithVersion`).
pub fn both_version(id: &ParsedSeqId, out: &mut Vec<u8>) {
    out.extend_from_slice(FASTA_TAGS[usize::from(id.alternative())]);
    out.push(b'|');
    label_content(id, true, false, out);
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:3416-3550
// `GetLabel` asks with a secondary list and `fAllowLocalId` (no `fGpipeAddSecondary`), and
// takes the primary id, or the secondary one when the primary is empty: so the label part is
// the primary id when there is one, else the (only) secondary id, both upper-cased (ASCII,
// `NStr::ToUpper` in the C locale) except the chain of a pdb id.
/// The part of `both_upper_version` that `ComposeOSLT` gives.
fn compose_oslt_label(id: &ParsedSeqId, out: &mut Vec<u8>) {
    let start = out.len();
    let mut upper_from = start;
    match id {
        ParsedSeqId::Giim(giim) => number(out, giim.id),
        ParsedSeqId::Gibbsq(value) | ParsedSeqId::Gibbmt(value) => number(out, value),
        ParsedSeqId::Pir(text) | ParsedSeqId::Prf(text) => {
            if let Some(name) = text.name.as_ref().or(text.accession.as_ref()) {
                out.extend_from_slice(name);
            }
        }
        ParsedSeqId::Patent(patent) => {
            out.extend_from_slice(&patent.cit.country);
            out.push(b'|');
            match &patent.cit.id {
                PatentNumber::Number(number) | PatentNumber::AppNumber(number) => {
                    out.extend_from_slice(number)
                }
            }
            out.push(b'|');
            number(out, patent.seqid);
        }
        ParsedSeqId::Pdb(pdb) => {
            out.extend_from_slice(&pdb.mol);
            out[start..].make_ascii_uppercase();
            if let Some(chain_id) = &pdb.chain_id {
                out.push(b'|');
                out.extend_from_slice(chain_id);
            } else if let Some(chain) = pdb.chain.filter(|&chain| chain != i32::from(b' ')) {
                out.push(b'|');
                out.push(chain as u8);
            }
            // `mixed_case`: the chain keeps its case.
            upper_from = out.len();
        }
        ParsedSeqId::General(dbtag) => {
            let suffix_start = out.len();
            out.extend_from_slice(&dbtag.db);
            out.push(b'|');
            let tag_start = out.len();
            object_id(&dbtag.tag, out);
            if out.len() == tag_start {
                // An empty suffix gives no secondary id.
                out.truncate(suffix_start);
            }
        }
        ParsedSeqId::Gi(gi) => number(out, gi),
        ParsedSeqId::Local(tag) => object_id(tag, out),
        _ => {
            let text = id.textseq().expect("a Textseq-id alternative");
            match &text.accession {
                Some(accession) if !accession.is_empty() => out.extend_from_slice(accession),
                // The primary id is empty: the name, when there is one, is the secondary id.
                _ => {
                    if let Some(name) = &text.name {
                        out.extend_from_slice(name);
                    }
                }
            }
        }
    }
    out[upper_from..].make_ascii_uppercase();
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:907-930
/// `CAlignFormatUtil::GetLabel` of `id` with `with_version`: a `general` id as its tag alone, every
/// other id (and a `general` id with an empty tag) as `content`.
pub fn format_label(id: &ParsedSeqId, with_version: bool, out: &mut Vec<u8>) {
    if let ParsedSeqId::General(dbtag) = id {
        let start = out.len();
        object_id(&dbtag.tag, out);
        if out.len() > start {
            return;
        }
    }
    content(id, with_version, out);
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4944-4956
/// `CAlignFormatUtil::GetBareId`: the FASTA form of a `gi`, `prf` or `pir` id, the content
/// with its version for every other id.
pub fn bare(id: &ParsedSeqId, out: &mut Vec<u8>) {
    match id {
        ParsedSeqId::Gi(_) | ParsedSeqId::Prf(_) | ParsedSeqId::Pir(_) => fasta(id, out),
        _ => content(id, true, out),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:155-189
/// `CShowBlastDefline::GetSeqIdListString` of `ids` with `show_gi`: the first GI in FASTA form
/// (with `show_gi`), then `|` and the FASTA form of the `Score`-best id unless that is the GI (a
/// local id as its content, without `lcl|`): `gi|123456|ref|NP_000002.1|`, `myprot9`.
pub fn seqid_list_string(ids: &[ParsedSeqId], show_gi: bool, out: &mut Vec<u8>) {
    let best = find_best(ids, text_score);
    let mut found_gi = false;
    if show_gi {
        if let Some(gi) = ids.iter().find(|id| id.gi().is_some()) {
            fasta(gi, out);
            found_gi = true;
        }
    }
    if let Some(best) = best.filter(|best| best.gi().is_none()) {
        if found_gi {
            out.push(b'|');
        }
        if let ParsedSeqId::Local(_) = best {
            content(best, false, out);
        } else {
            fasta(best, out);
        }
    }
}

// NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_id.hpp:1041-1047
/// `FindGi`: the first GI of `ids`, 0 for none.
pub fn first_gi<'a>(ids: impl IntoIterator<Item = &'a ParsedSeqId>) -> i64 {
    ids.into_iter().find_map(ParsedSeqId::gi).unwrap_or(0)
}
