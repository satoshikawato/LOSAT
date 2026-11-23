//! `sequence::GetTitle` (flags 0): the title of a nucleotide Bioseq of the FASTA or BLAST
//! database readers, as windowmasker matches its `-ids` against it (`CWinMaskUtil::consider`)
//! and makembindex writes its first word to the `-idmap` file (`extractSeqVector`). Unlike
//! `report::defline::generate_bioseq_defline` (`CDeflineGenerator`), it decodes no HTML and
//! trims only the dots and spaces at the end.

use std::sync::Arc;

use crate::seqid::model::PatentNumber;
use crate::seqid::ParsedSeqId;

/// The title NCBI gives a Bioseq whose source makes none (`GetTitle`).
const NO_DEFINITION_LINE: &[u8] = b"No definition line found";

/// `sequence::GetTitle` (flags 0) of a nucleotide Bioseq with the Seq-ids `ids` and the title
/// descriptor `title`, which has a source when `has_source`; `taxname` gives the scientific name
/// of its first source (`None` without one), asked only when the title is made from the source.
/// The readers' Bioseqs have no MolInfo, no GenBank or EMBL block and no PDB descriptor, are
/// raw or delta, and a database's source has only its organism's names and taxid
/// (seqdbvol.cpp:1165-1174), so of the function only these steps apply.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/seqtitle.cpp:131-178,286-296,371-436
/// (The steps for an empty title from the source of an NC, NM or NR accession, from a PDB
/// descriptor or from a protein do not apply: they need a MolInfo, a feature or a descriptor
/// the readers do not make, or a protein. The flags add no organism.)
pub(crate) fn get_title<E>(
    ids: &[ParsedSeqId],
    title: &[u8],
    has_source: bool,
    taxname: impl FnOnce() -> Result<Option<Arc<[u8]>>, E>,
) -> Result<Vec<u8>, E> {
    let mut third_party = false;
    let mut patent = None;
    for id in ids {
        match id {
            ParsedSeqId::Tpg(_) | ParsedSeqId::Tpe(_) | ParsedSeqId::Tpd(_) => third_party = true,
            ParsedSeqId::Patent(pat) => patent = Some(pat),
            _ => {}
        }
    }
    let mut text = title.to_vec();
    while matches!(text.last(), Some(b'.' | b' ')) {
        text.pop();
    }
    if text.is_empty() {
        if let Some(pat) = patent {
            let number = match &pat.cit.id {
                PatentNumber::Number(number) | PatentNumber::AppNumber(number) => number,
            };
            text = format!("Sequence {} from Patent ", pat.seqid).into_bytes();
            text.extend_from_slice(&pat.cit.country);
            text.push(b' ');
            text.extend_from_slice(number);
        }
    }
    if text.is_empty() && has_source {
        text = title_from_bio_source(taxname()?.as_deref());
        if text.is_empty() {
            text = NO_DEFINITION_LINE.to_vec();
        }
    }
    if third_party && !text.is_empty() && !starts_with_nocase(&text, b"TPA: ") {
        text.splice(0..0, b"TPA: ".iter().copied());
    }
    Ok(text)
}

/// `s_TitleFromBioSource` of a source whose organism is named `taxname` and has no subtype, no
/// modifier and no MolInfo technique (a database's): the name without the white space around it
/// (`NStr::TruncateSpaces`, C locale `isspace`), its first letter upper-cased (C locale).
// NCBI reference (598d8ae6): c++/src/objmgr/util/seqtitle.cpp:919-979
fn title_from_bio_source(taxname: Option<&[u8]>) -> Vec<u8> {
    let name = taxname.unwrap_or_default();
    let is_space = |byte: &u8| matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r');
    let start = name.iter().position(|byte| !is_space(byte));
    let mut title = match start {
        Some(start) => {
            let end = name
                .iter()
                .rposition(|byte| !is_space(byte))
                .unwrap_or(start);
            name[start..=end].to_vec()
        }
        None => Vec::new(),
    };
    if let Some(first) = title.first_mut() {
        first.make_ascii_uppercase();
    }
    title
}

/// `NStr::StartsWith` with `NStr::eNocase`.
fn starts_with_nocase(text: &[u8], prefix: &[u8]) -> bool {
    text.len() >= prefix.len() && text[..prefix.len()].eq_ignore_ascii_case(prefix)
}

/// `get_title` of a Bioseq without a source (the FASTA readers').
pub(crate) fn get_title_without_source(ids: &[ParsedSeqId], title: &[u8]) -> Vec<u8> {
    let Ok(text) = get_title(ids, title, false, || {
        Ok::<_, std::convert::Infallible>(None)
    });
    text
}

#[cfg(test)]
#[path = "title_tests.rs"]
mod tests;
