//! The subject title of the pairwise report (outfmt 0) and the subject ID of the tabular
//! formats as NCBI makes them for a sequence read from FASTA:
//! `CDeflineGenerator::GenerateDefline` over the title of the record (the defline after
//! `>`, its first word included, as bytes), and `CShowBlastDefline::GetSeqIdList` over it.
//!
//! NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3952-3962
//! NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:4050-4073
//! NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:4092-4101
//! A local FASTA record has no MolInfo, source or metagenome data, so `mag`, `prefix` and
//! `suffix` are empty, and only an empty title is generated (`capitalize` stays true only
//! then): a protein's is `unnamed protein product` (`x_SetTitleFromProtein`), a nucleotide's
//! is empty, so nothing is capitalized. The alignment heading removes the prefixes; the
//! description table keeps them (`fLeavePrefixSuffix`, showdefline.cpp:498).
//! `NStr::HtmlDecode` guesses the encoding of the whole title (`CUtf8::GuessEncoding`) and
//! throws when it cannot; NCBI then shows `Unknown` and skips the alignments of the subject
//! (`UnknownEncoding`, see `report::pairwise`). For some titles of commas, semicolons,
//! tildes and spaces, NCBI's `x_CleanAndCompress` lets its count of the remaining letters
//! wrap and reads past the end of the string (and crashes); `clean_and_compress` stops at
//! the end of the string there (approved exception 2 of PD-LOSAT-NCBI-DEFECTS).

use crate::ncbi_str::replace_in_place;

/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3431-3446
const TPA_PREFIXES: [&str; 14] = [
    "MAG ",
    "MAG:",
    "MULTISPECIES:",
    "TLS:",
    "TPA:",
    "TPA_exp:",
    "TPA_inf:",
    "TPA_reasm:",
    "TPA_asm:",
    "TSA:",
    "UNVERIFIED_ORG:",
    "UNVERIFIED_ASMBLY:",
    "UNVERIFIED_CONTAM:",
    "UNVERIFIED:",
];

/// The title that NCBI generates for a protein record without one.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2508-2509
/// (`x_SetTitleFromProtein`, called for an empty title of a protein at
/// create_defline.cpp:4020-4021.)
const UNNAMED_PROTEIN_PRODUCT: &[u8] = b"unnamed protein product";

/// `GenerateDefline` threw: `NStr::HtmlDecode` could not guess the encoding of the title
/// (`CUtf8::GuessEncoding` gave `eEncoding_Unknown`: a title that is not UTF-8 with one of
/// the bytes 0x81, 0x8D, 0x8F, 0x90 and 0x9D).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4528-4534
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct UnknownEncoding;

/// A title made by `generate_defline`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Defline {
    /// The bytes NCBI writes.
    pub text: Vec<u8>,
    /// Whether NCBI's `x_CleanAndCompress` reads past the end of the title (NCBI crashes
    /// there; `text` stops at the end of the string).
    pub reads_past_end: bool,
}

/// `CDeflineGenerator::GenerateDefline` of a record from FASTA with the title `title` (see
/// the module): for the alignment heading (`leave_prefix_suffix` false, flags 0) or for
/// the description table (true, `fLeavePrefixSuffix`), of a protein (`m_IsAA`) or a
/// nucleotide record.
pub fn generate_defline(
    title: &[u8],
    protein: bool,
    leave_prefix_suffix: bool,
) -> Result<Defline, UnknownEncoding> {
    generate_bioseq_defline(title, protein, leave_prefix_suffix, &BioseqFacts::default())
}

/// What `CDeflineGenerator::x_SetFlags` takes from the Seq-ids of a Bioseq (none for a record
/// read from FASTA, whose id is local): the Bioseq of a database sequence with parsed Seq-ids
/// has the ids of one definition line (`report::db_subject`).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:745-829
/// The other flags of `x_SetFlags` come from descriptors that the Bioseq of a database
/// sequence does not have (MolInfo, user objects, GenBank and EMBL blocks), or change no title
/// that a search writes.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1004-1011
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1264-1274
/// The Bioseq of a database sequence has a source for each definition line with a positive
/// taxid (`CSeqDBVol::x_GetTaxonomy`, seqdbvol.cpp:1136-1180), whose Org-ref has only the taxid's
/// dbtag and, where the taxonomy database names the taxid, its names: no genome, subsource or
/// orgname (so no organelle, strain or metagenome source), and the Bioseq has no MolInfo (so it
/// is complete) and no features. `taxname` is the name of the first source
/// (`common::database::HeaderReader::source_taxname`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BioseqFacts {
    /// `m_ThirdParty`: a `tpg`, `tpe` or `tpd` id (the prefix `TPA: `).
    pub third_party: bool,
    /// `m_IsNC` or `m_IsNR`: an accession of a RefSeq chromosome or ncRNA, whose empty title
    /// NCBI makes with `x_SetTitleFromNC` or `x_SetTitleFromNR` and then from the source (a
    /// protein's is then not `unnamed protein product`; `source_title`).
    pub refseq_chromosome_or_ncrna: bool,
    /// `m_IsNM`: an accession of a RefSeq mRNA (`NM_`, `XM_`), whose empty title NCBI makes with
    /// `x_SetTitleFromNM` and then from the source, unless `local_annots_only` (`source_title`).
    pub refseq_mrna: bool,
    /// `m_IsPDB`, `m_PDBChainID` and `m_PDBChain`: whether a `pdb` id was seen and its chain.
    pub pdb: Option<(Vec<u8>, i32)>,
    /// `m_IsPatent`, `m_PatentSequence`, `m_PatentCountry` and `m_PatentNumber`.
    pub patent: Option<(i32, Vec<u8>, Vec<u8>)>,
    /// `m_Taxname`: the scientific name of the Bioseq's first source, `None` without a source
    /// or a name (never empty).
    pub taxname: Option<std::sync::Arc<[u8]>>,
    /// `m_Source` is set: some definition line has a positive taxid
    /// (`DeflineSet::first_source_taxid`), with or without a name. `x_AdjustProteinTitleSuffix`
    /// runs only then (`adjust_protein_title_suffix`); a name implies it.
    pub has_source: bool,
    /// `m_LocalAnnotsOnly` (`fLocalAnnotsOnly`, set by blastdbcmd's `CFastaOstream` with
    /// `fNoExpensiveOps`): an empty protein title is not `unnamed protein product`.
    pub local_annots_only: bool,
}

/// The name of the organism of a database sequence's first source (`BioseqFacts::taxname`), as
/// the reader of the reported headers finds it (`common::database::HeaderReader::source_taxname`,
/// asked only for a title that `GenerateDefline` makes from the source, where a taxonomy database
/// opens).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum SourceTaxname {
    /// No name: the title is not made from the source, the Bioseq has no source, no taxonomy
    /// database opens, or it gives no scientific name for the taxid. The title is made as of a
    /// Bioseq without a source.
    #[default]
    Unnamed,
    /// The scientific name of the taxid (not empty).
    Named(std::sync::Arc<[u8]>),
    /// The name `name` of a database sequence without parsed Seq-ids whose title NCBI cannot
    /// make: `NStr::HtmlDecode` cannot guess its encoding and throws (`UnknownEncoding`), so the
    /// outfmt 0 alignment heading names the sequence by `seq_id`, the `CSeq_id::GetSeqIdString`
    /// of its `gnl|BL_ORD_ID|n` (`BL_ORD_ID:n`), in place of its alignments
    /// (`report::pairwise::DatabaseHeading::skipped_id`). The other reports read it as `Named`.
    Undecodable {
        name: std::sync::Arc<[u8]>,
        seq_id: std::sync::Arc<[u8]>,
    },
    /// A source without a name: the definition lines have a positive taxid that no taxonomy
    /// database names. Asked only for a database sequence without parsed Seq-ids whose protein
    /// title's first word can change where `x_AdjustProteinTitleSuffix` rewrites its end
    /// (`protein_title_word_needs_source`; `common::database::HeaderReader::stored_source`), which
    /// runs for any source: it then drops a ` [ ]` end. Elsewhere it reads as `Unnamed`.
    Nameless,
    /// LOSAT's refusal where a report shows the title, with what the sequence has (after
    /// "database sequence N has"): a lookup of the taxid whose scientific name NCBI BLAST+ reads
    /// past its mapping of `taxdb.btd` (an empty file: NCBI faults; otherwise its bytes, or
    /// whether NCBI faults, depend on its memory map), or a definition line whose taxid does not
    /// read.
    Refused(&'static str),
}

/// `SourceTaxname::Unnamed`, for a subject without a source.
pub static UNNAMED_SOURCE: SourceTaxname = SourceTaxname::Unnamed;

impl SourceTaxname {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1165-1169
    // NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1271-1273
    /// The source of a Bioseq whose first source's taxid the taxonomy database gives `names`
    /// (`seqdb::taxdb::TaxCache::names`): its scientific name, none where the lookup found none or
    /// it is empty (an empty `m_Taxname` makes the title of a Bioseq without one,
    /// `seqdb::taxdb::SourceNames::taxname`).
    pub fn of_names(names: &crate::seqdb::taxdb::SourceNames) -> SourceTaxname {
        match names.taxname() {
            Some(name) => SourceTaxname::Named(std::sync::Arc::clone(name)),
            None => SourceTaxname::Unnamed,
        }
    }

    /// The name for `BioseqFacts::taxname`; `Err` with the reason of a refusal.
    pub fn name(&self) -> Result<Option<std::sync::Arc<[u8]>>, &'static str> {
        match self {
            SourceTaxname::Unnamed | SourceTaxname::Nameless => Ok(None),
            SourceTaxname::Named(name) | SourceTaxname::Undecodable { name, .. } => {
                Ok(Some(std::sync::Arc::clone(name)))
            }
            SourceTaxname::Refused(reason) => Err(reason),
        }
    }

    /// Whether the source is known to exist (`BioseqFacts::has_source`): a name, or `Nameless`.
    /// `Unnamed` is no source, or one that no report asked about.
    pub fn has_source(&self) -> bool {
        matches!(
            self,
            SourceTaxname::Named(_) | SourceTaxname::Undecodable { .. } | SourceTaxname::Nameless
        )
    }

    /// The source of the database sequence without parsed Seq-ids whose Seq-id text is
    /// `seq_id` (`SourceTaxname::Undecodable`, made where the title made from the name throws).
    pub fn of_stored(self, seq_id: impl FnOnce() -> Vec<u8>) -> SourceTaxname {
        match self {
            SourceTaxname::Named(name)
                if generate_bioseq_defline(
                    b"",
                    false,
                    true,
                    &BioseqFacts::of_source(Some(std::sync::Arc::clone(&name))),
                )
                .is_err() =>
            {
                SourceTaxname::Undecodable {
                    name,
                    seq_id: std::sync::Arc::from(seq_id()),
                }
            }
            source => source,
        }
    }

    /// The Seq-id text that the outfmt 0 alignment heading of a sequence whose title NCBI cannot
    /// make names (`SourceTaxname::Undecodable`).
    pub fn undecodable_id(&self) -> Option<&[u8]> {
        match self {
            SourceTaxname::Undecodable { seq_id, .. } => Some(seq_id),
            _ => None,
        }
    }
}

impl BioseqFacts {
    /// The facts of a database sequence without parsed Seq-ids (its one id is
    /// `gnl|BL_ORD_ID|n`, which sets no flag) whose first source names the organism `taxname`
    /// (a name implies a source; a source without a name is `SourceTaxname::has_source`).
    pub fn of_source(taxname: Option<std::sync::Arc<[u8]>>) -> Self {
        Self {
            has_source: taxname.is_some(),
            taxname,
            ..Self::default()
        }
    }
}

/// `GenerateDefline` of a Bioseq with the title `title` and the facts `ids` (`BioseqFacts`):
/// `generate_defline` with the parts of the title that come from the Seq-ids and the source, the
/// generated title of a `pdb` or `patent` record without a title, of a record without a title
/// from its source's name (`title_before_decoding_of`), the prefix `TPA: ` of a third-party
/// record (alignment heading only, as the other prefixes), and the first letter of a nucleotide
/// title made from nothing upper-cased.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:4075-4093
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2815-2840
pub fn generate_bioseq_defline(
    title: &[u8],
    protein: bool,
    leave_prefix_suffix: bool,
    ids: &BioseqFacts,
) -> Result<Defline, UnknownEncoding> {
    let (mut decoded, _) = html_decode(&title_before_decoding_of(
        title,
        protein,
        leave_prefix_suffix,
        ids,
    ))?;
    trim_end_of(&mut decoded, b",;~ ");
    if !leave_prefix_suffix && ids.third_party {
        decoded.splice(0..0, b"TPA: ".iter().copied());
    }
    let (mut text, reads_past_end) = clean_and_compress(&decoded, protein);
    // NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3934-3962
    // NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:4095-4099
    // A title of punctuation is not erased, so `capitalize` stays true only for an empty one;
    // the Bioseq of a database sequence is raw (`m_IsSeg` false). `islower` and `toupper` are
    // the C locale's: ASCII only.
    if title.is_empty() && !protein && ids.pdb.is_none() && ids.patent.is_none() {
        if let Some(first) = text.first_mut() {
            first.make_ascii_uppercase();
        }
    }
    Ok(Defline {
        text,
        reads_past_end,
    })
}

/// Whether `generate_bioseq_defline` makes the title of `title` and `ids` from the source of
/// the Bioseq (`x_SetBioSrc`, `x_SetTitleFromBioSrc`, ...: a title that is empty once trimmed,
/// of a record that is neither a `pdb` record with a chain nor a `patent` one). A database
/// sequence's Bioseq has a source for each definition line with a taxid other than 0, with the
/// names of the taxonomy database when NCBI finds one (seqdbvol.cpp:1055-1160): the title is
/// then the organism's (`unnamed protein product [Homo sapiens]`). `ids` gives the Bioseq's
/// facts, asked only for a title that is empty once trimmed (a titled sequence, the readers'
/// usual one, needs none).
pub fn bioseq_title_from_source(title: &[u8], ids: impl FnOnce() -> BioseqFacts) -> bool {
    // `trim_end_of` without a copy: the stored titles of every reported sequence are asked.
    let main_title = match title.iter().rposition(|byte| !b".,;~ ".contains(byte)) {
        Some(pos) => &title[..=pos],
        None => title,
    };
    main_title.is_empty() && pdb_or_patent_title(&ids()).is_empty()
}

// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1815-1858
// The Bioseq of a database sequence has no PDB block and no comment: `m_PDBCompound` and
// `m_Comment` are empty.
/// The title that `x_SetTitleFromPDB` or `x_SetTitleFromPatent` makes for an empty one (empty
/// for a `pdb` record without a chain, and for a record of neither kind).
fn pdb_or_patent_title(ids: &BioseqFacts) -> Vec<u8> {
    if let Some((chain_id, chain)) = &ids.pdb {
        let mut title = b"Chain ".to_vec();
        if !chain_id.is_empty() {
            title.extend_from_slice(chain_id);
        } else if (0x20..0x7f).contains(&(*chain as u8)) {
            // `isprint` of `m_PDBChain` cast to `unsigned char`, and the one-character string of it cast
            // to `char`: the low byte of the integer (`chain 321` is `A`).
            title.push(*chain as u8);
        } else {
            return Vec::new();
        }
        title.extend_from_slice(b", ");
        return title;
    }
    if let Some((sequence, country, number)) = &ids.patent {
        let mut title = format!("Sequence {sequence} from Patent ").into_bytes();
        title.extend_from_slice(country);
        title.push(b' ');
        title.extend_from_slice(number);
        return title;
    }
    Vec::new()
}

/// The subject ID of the tabular formats (`sseqid`, `sacc`, `saccver`) for a subject with
/// the local ID `local_id` (`Subject_N`) and the title `title`: the title's first word
/// (`s_ReplaceLocalId`; the local ID for an empty title), and when that names a local
/// `Subject_` ID, the first word of the title that `GenerateDefline` makes (an empty
/// protein title gives `unnamed`). `Err` where `GenerateDefline` throws (NCBI stops with
/// an exception; RP-20).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:787-789
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:522-532
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:482-496
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:218-242
/// A local ID is written `lcl|` and its text (`CSeq_id::WriteAsFasta`, Seq_id.cpp:2164-2196,
/// `CObject_id::AsString`, Object_id.cpp:202-210). Of `lcl|` and the text, `lcl|Subject_` is
/// found where the text starts with `Subject_` or holds it all: the `lcl|` before the text takes
/// part in no other match.
#[inline]
pub fn tabular_subject_id(
    local_id: &[u8],
    title: &[u8],
    protein: bool,
) -> Result<Vec<u8>, UnknownEncoding> {
    let replaced = first_word(title);
    let replaced = if replaced.is_empty() {
        local_id
    } else {
        replaced
    };
    if replaced.starts_with(b"Subject_") || contains(replaced, b"lcl|Subject_") {
        let defline = generate_defline(title, protein, false)?.text;
        let token = first_word(&defline);
        if !token.is_empty() {
            return Ok(token.to_vec());
        }
    }
    Ok(replaced.to_vec())
}

/// The words of LOSAT's refusal of a protein title of punctuation that NCBI BLAST+'s
/// `x_CleanAndCompress` reads past the end of (NCBI BLAST+ 2.17.0 crashes; approved exception 2
/// of PD-LOSAT-NCBI-DEFECTS has a nucleotide title stop at the end of the string), after
/// "a title ": one text for the tabular ids (`database_subject_id`,
/// `ParsedSubject::generated_title_word`) and the organism report (`tax_report::generated_title`).
macro_rules! punctuation_title {
    () => {
        "of punctuation that NCBI BLAST+'s x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes)"
    };
}
pub(crate) use punctuation_title;

/// The subject ID of the tabular formats (`sseqid`, `sacc`, `saccver`) of a database
/// sequence without parsed Seq-ids (its only id is `gnl|BL_ORD_ID|n`) with the stored
/// title `title`: the first word of the title that `GenerateDefline` makes of it
/// (flags 0), which `GetSeqIdList` puts in the place of the id.
///
/// `SubjectWord::General` for a title whose generated first word is empty (a title of spaces;
/// `MAG:` alone; a nucleotide's empty title whose source has no name; `MAG: (plastid) [ ]` of a
/// protein with a source, see below): NCBI keeps the general id and writes it in two forms,
/// `gnl|BL_ORD_ID|n` for `sseqid` and `BL_ORD_ID:n` for `sacc` and `saccver`, n the sequence's
/// OID (`algorithm::common::database::TabularSubjectId`). Where the taxonomy database names the
/// taxid of an empty title, NCBI makes the title from the name (`bioseq_title_from_source`,
/// `ids`), whose first word is the id (`Homo` for taxid 9606).
/// `SubjectWord::Unknown` for a title (or name) whose encoding NCBI cannot guess:
/// `GenerateDefline` throws and `GetSeqIdList` makes no id (`report::tabular::SubjectLookup`).
///
/// `Err` names what LOSAT does not write as NCBI: a title of punctuation that NCBI's
/// `x_CleanAndCompress` reads past the end of (NCBI crashes).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:218-242
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:522-532
/// The Bioseq of a database sequence carries its stored title as its title descriptor
/// (seqdbvol.cpp:1445-1463) and the source of its taxid (`ids`, `BioseqFacts::of_source`), so
/// `GenerateDefline` reads the title as it reads the one of a FASTA record, and makes an empty
/// one from the source's name (`Homo` of `Homo sapiens`; `unnamed` of a protein). The
/// empty word comes after the title's prefix goes (create_defline.cpp:4050-4062, `MAG: (plastid)
/// [ ]` adjusted to `MAG:`, then emptied: probe of the review of line RG item T5c,
/// `remaining-gaps-20261112/T/t5c/review/out/ncbi/sall.*`) or of a title that stays empty.
pub fn database_subject_id(
    title: &[u8],
    ids: &BioseqFacts,
    protein: bool,
) -> Result<SubjectWord, &'static str> {
    let Ok(defline) = generate_bioseq_defline(title, protein, false, ids) else {
        return Ok(SubjectWord::Unknown);
    };
    if defline.reads_past_end {
        return Err(punctuation_title!());
    }
    let token = first_word(&defline.text);
    Ok(if token.is_empty() {
        SubjectWord::General
    } else {
        SubjectWord::Word(token.to_vec())
    })
}

/// What `GetSeqIdList` puts in the place of a database sequence's general id
/// (`database_subject_id`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SubjectWord {
    /// The first word of the generated title.
    Word(Vec<u8>),
    /// The general id itself: the generated title's first word is empty.
    General,
    /// No id: `GenerateDefline` throws (a title in no encoding that NCBI guesses).
    Unknown,
}

/// The first token that `NStr::Split` makes of `text` at a space: the bytes before the first space.
fn first_word(text: &[u8]) -> &[u8] {
    let end = text
        .iter()
        .position(|&byte| byte == b' ')
        .unwrap_or(text.len());
    &text[..end]
}

fn contains(text: &[u8], pattern: &[u8]) -> bool {
    text.windows(pattern.len()).any(|window| window == pattern)
}

/// The title that `GenerateDefline` gives `NStr::HtmlDecode` (`m_MainTitle`, see the
/// module), for the alignment heading (`leave_prefix_suffix` false) or the description
/// table, of a Bioseq with the facts `ids`: an empty title of a `pdb` or `patent` record is
/// generated from the id, else from the source (`source_title`) (`generate_bioseq_defline`).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3994-4022
fn title_before_decoding_of(
    title: &[u8],
    protein: bool,
    leave_prefix_suffix: bool,
    ids: &BioseqFacts,
) -> Vec<u8> {
    let mut main_title = title.to_vec();
    trim_end_of(&mut main_title, b".,;~ ");
    // NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3964-3972
    // `CDeflineGenerator()` has no feature index (`m_Idx` null). The adjustment comes after the
    // trim and before a title is made for an empty one (which is never adjusted), the prefixes
    // are removed and the title is decoded; it returns at once without a source.
    if protein && !leave_prefix_suffix && ids.has_source && !main_title.is_empty() {
        adjust_protein_title_suffix(&mut main_title, ids.taxname.as_deref().unwrap_or_default());
    }
    if main_title.is_empty() {
        main_title = pdb_or_patent_title(ids);
    }
    if main_title.is_empty() {
        main_title = source_title(protein, ids);
    }
    if !leave_prefix_suffix {
        for prefix in TPA_PREFIXES {
            if main_title.len() >= prefix.len()
                && main_title[..prefix.len()].eq_ignore_ascii_case(prefix.as_bytes())
            {
                main_title.drain(..prefix.len());
                trim_leading_spaces(&mut main_title);
            }
        }
    }
    trim_leading_spaces(&mut main_title);
    main_title
}

// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3090-3110
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2933-2977 (`NStr::Find`, `eNocase`)
// The reverse search gives the last place at or before the title's length less the name's (`slen`
// less `plen`) where the pattern matches without case (ASCII: the C locale): an empty pattern has
// no first letter, so `find_last_of` finds nothing (an empty name never matches), and a stop at 0
// without a match is never `idx + 2` (at least 2). So the place is `idx + 2` (one less than `slen`
// less `plen`) exactly where the name is there and not at `slen` less `plen`: only those two places
// are compared.
/// `s_TitleEndsInOrganism`: the position of `[` where `title` ends in ` [<taxname>]` (the name in
/// any ASCII case).
fn title_ends_in_organism(title: &[u8], taxname: &[u8]) -> Option<usize> {
    let (len1, len2) = (title.len(), taxname.len());
    if len2 == 0 || len1 <= len2 + 4 {
        return None;
    }
    let idx = len1 - len2 - 3;
    if title[idx] != b' ' || title[idx + 1] != b'[' || title[len1 - 1] != b']' {
        return None;
    }
    let found_at = |at: usize| title[at..at + len2].eq_ignore_ascii_case(taxname);
    (!found_at(len1 - len2) && found_at(idx + 2)).then_some(idx + 1)
}

// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2123-2147
// The commented-out names are not used: `(virus)` stays (probe of line RG T5c, the title "nu (x)(virus)
// [Bacillus sp.]" unchanged).
/// `s_proteinOrganellePrefix` from `chloroplast` (2) to `chromatophore` (22), the non-empty ones.
const PROTEIN_ORGANELLE_PREFIXES: [&[u8]; 15] = [
    b"chloroplast",
    b"chromoplast",
    b"kinetoplast",
    b"mitochondrion",
    b"plastid",
    b"macronuclear",
    b"plasmid",
    b"cyanelle",
    b"nucleomorph",
    b"apicoplast",
    b"leucoplast",
    b"protoplast",
    b"endogenous virus",
    b"hydrogenosome",
    b"chromatophore",
];

// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3263-3414
// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3049-3055 (`s_TrimMainTitle`,
// `trim_end_of(.., b".,;~ ")`)
// The source of a database Bioseq (`CSeqDBVol::x_GetTaxonomy`, seqdbvol.cpp:1161-1176) has the
// taxid's dbtag and, where the taxonomy database names the taxid, its names: no genome and no
// orgname, so no genus or species (the binomial is `" "`, which matches a title ending in
// ` [ ]`), no organelle and not cross-kingdom; the Bioseq has no MolInfo (`x_IsComplete()`,
// :3416-3428, holds: no `, partial` back) and no CDS feature for `x_GetSourceFeatViaCDS`
// (:2044-2056); `m_OmitTaxonomicName` comes from a flag no caller sets. NCBI BLAST+ 2.17.0 writes
// these bytes with and without a taxonomy database (probe of line RG T5c, `T/t5c/out/ncbi`:
// `alpha, partial [Bacillus sp.]` -> `alpha [Bacillus sp.]`, `zeta protein, partial [ ]` ->
// `zeta protein [Bacillus sp.]` or, without a name, `zeta protein`).
/// `x_AdjustProteinTitleSuffix` of the title `title` of a database protein Bioseq with a source
/// whose scientific name is `taxname` (empty without one): a title that ends in ` [taxname]` (any
/// ASCII case) or ` [ ]` loses that, then a `(organelle)` and `, partial` before it, and gets
/// ` [taxname]` back where the name is not empty.
pub(crate) fn adjust_protein_title_suffix(title: &mut Vec<u8>, taxname: &[u8]) {
    trim_end_of(title, b".,;~ ");
    if title.len() <= taxname.len() + 4 {
        return;
    }
    let Some(tpos) =
        title_ends_in_organism(title, taxname).or_else(|| title_ends_in_organism(title, b" "))
    else {
        return;
    };
    title.truncate(tpos);
    trim_end_of(title, b".,;~ ");
    let mut opos = None;
    if title.len() > 2 && title[title.len() - 1] == b')' {
        if let Some(pos) = title.iter().rposition(|&c| c == b'(') {
            // `NStr::EndsWith` of `m_MainTitle` and `str` in parentheses; the title ends in `)`.
            let len = title.len();
            let ends_in_paren = |prefix: &[u8]| {
                len >= prefix.len() + 2
                    && title[len - prefix.len() - 2] == b'('
                    && &title[len - prefix.len() - 1..len - 1] == prefix
            };
            if PROTEIN_ORGANELLE_PREFIXES
                .iter()
                .any(|prefix| ends_in_paren(prefix))
            {
                opos = Some(pos);
            }
        }
        trim_end_of(title, b".,;~ ");
    }
    if let Some(opos) = opos {
        title.truncate(opos);
        trim_end_of(title, b".,;~ ");
    }
    if title.ends_with(b", partial") {
        title.truncate(title.len() - 9);
        trim_end_of(title, b".,;~ ");
    }
    if !taxname.is_empty() {
        title.extend_from_slice(b" [");
        title.extend_from_slice(taxname);
        title.push(b']');
    }
}

/// Whether `adjust_protein_title_suffix` can change `title` with some source: its main title
/// (after `s_TrimMainTitle`) ends in `]` (`s_TitleEndsInOrganism` asks for it). The readers that
/// show the whole title ask for the name of a database sequence's source only then (and for an
/// empty title); the search's tabular ids, which show its first word, only where
/// `protein_title_word_needs_source`.
pub fn protein_title_suffix_adjustable(title: &[u8]) -> bool {
    title
        .iter()
        .rposition(|byte| !b".,;~ ".contains(byte))
        .is_some_and(|end| title[end] == b']')
}

// No NCBI counterpart: which titles `x_AdjustProteinTitleSuffix` (cited at
// `adjust_protein_title_suffix`) cannot change the first word of; it changes no value NCBI
// computes. Let M be the main title after `s_TrimMainTitle`, ending in `]`, w its bytes before its
// first space (at s) and c the byte after that space. The adjustment changes M only in two ways:
// (a) M ends in ` [X]` with X the source's name in some ASCII case (`s_TitleEndsInOrganism`,
// a name that is not empty): M becomes P + ` [name]`; (b) M ends in ` [ ]` (the binomial
// `" "` of a source without genus and species): M becomes P, or P + ` [name]` for any name. P is
// M cut at the `[` of a ` [`, then stripped from its end of `.,;~ `, of one `(organelle)` and of
// one `, partial`, all ASCII. Where w is not empty and holds only the "plain" bytes below (no
// `.,;~` that the strip takes, no `(`, `)` or `,` that its organelle and `, partial` end in), and c
// is plain or `[` (no `(` that would begin an organelle at s + 1), each stage of the strip stops
// at w (the space at s goes, and c with it only where c is the cut `[`) or keeps w, the space and
// c: P is w, or begins with w, a space and c. In (a) the bytes of M' and M at and above 0x80 are
// the same, in runs that ASCII bytes end in both (P's and the name's), so
// `CUtf8::GuessEncoding` gives both an encoding that fails or both one that keeps ASCII bytes as
// they are (UTF-8 and CESU-8 keep every byte). Then the passes after the adjustment read w, the
// space and c alike in both: no TPA prefix matches (all but `MAG ` end in `:`, which w lacks, and
// w is not `MAG`), `NStr::HtmlDecode` keeps them (no `&`), `x_CleanAndCompress` writes w and then
// the space (no pair of its patterns: `,,`, two spaces, `( `, ` )`, ` ,`, ` ;`, `, `, `; `), which
// its `. [` and `, [` replacements do not reach (w ends in neither `.` nor `,`), and it does not
// read past the end of a title ending in `]`. So the first word is w with any source, or the
// title is undecodable with any source, as without one. (b) can put any name in M, whose bytes
// can make the encoding fail: such a title needs the source.
/// Whether the first word of the title that `GenerateDefline` makes of the stored protein title
/// `title` (`database_subject_id`, the tabular ids of a database sequence without parsed
/// Seq-ids) can depend on the Bioseq's source (its name, or whether it has one), which
/// `x_AdjustProteinTitleSuffix` reads (`adjust_protein_title_suffix`): a title ending in `]`
/// unless its first word is made of bytes from `!` to `~` other than `& , ; . ( ) : ~`, is not
/// `MAG` (any case), and is followed by a space and such a byte or `[`, and the title does not end
/// in ` [ ]` (see above).
pub fn protein_title_word_needs_source(title: &[u8]) -> bool {
    /// The plain bytes: `!` to `~` but `& , ; . ( ) : ~`.
    const PLAIN: [bool; 256] = {
        let mut plain = [false; 256];
        let mut byte = b'!';
        while byte <= b'~' {
            plain[byte as usize] =
                !matches!(byte, b'&' | b',' | b';' | b'.' | b'(' | b')' | b':' | b'~');
            byte += 1;
        }
        plain
    };
    let Some(end) = title.iter().rposition(|byte| !b".,;~ ".contains(byte)) else {
        return false;
    };
    let main = &title[..=end];
    if main[end] != b']' {
        return false;
    }
    if main.ends_with(b" [ ]") {
        return true;
    }
    // The first word, up to the first space; the title ends in `]`, so a byte follows that space.
    let Some(space) = main.iter().position(|&byte| !PLAIN[byte as usize]) else {
        return true;
    };
    let next = main[space + 1..].first().copied().unwrap_or(b' ');
    !(space > 0
        && main[space] == b' '
        && !main[..space].eq_ignore_ascii_case(b"MAG")
        && (PLAIN[next as usize] || next == b'['))
}

// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1545-1550
// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1598-1600
// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:91-99 (`CDefLineJoiner::Add`)
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3197-3234 (`NStr::TruncateSpacesInPlace`
// with `eTrunc_Both`: the `isspace` bytes at both ends, of the C locale)
// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2508-2510, 2532-2541, 2572-2574
// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1603-1610, 1715-1716, 1763-1764
// The Bioseq of a database sequence has no features, no MolInfo (so `x_IsComplete()` holds:
// no `, partial`), no genome (no organelle), no CDS (so no source feature for `synthetic
// construct` and the like) and no orgname (not cross-kingdom); of `x_SetTitleFromBioSrc`'s
// parts only the organism remains. A RefSeq chromosome or ncRNA record, and an mRNA record
// without `m_LocalAnnotsOnly`, skips `x_SetTitleFromProtein` for `x_SetTitleFromNC`,
// `x_SetTitleFromNR` or `x_SetTitleFromNM`, which make nothing of a database Bioseq
// (`x_SetTitleFromNC` returns at its MolInfo test, :1606-1607; `x_SetTitleFromNM` needs one gene
// and one CDS feature, :1741; `x_SetTitleFromNR` sets its title in its gene feature loop only,
// :1766-1768), and then gets `x_SetTitleFromBioSrc` (create_defline.cpp:4033-4036): the name,
// also of a protein (probe of line RG T5, `T/t5/rs`: `ref|NC_000001.1|` and `ref|NM_000001.1|`
// proteins `Homo sapiens`, blastdbcmd `-long_seqids` of the `NM_` one `[Homo sapiens]`).
/// The title that `GenerateDefline` makes of the source of a Bioseq without a title (`ids`):
/// a protein's is `unnamed protein product` (none with `local_annots_only`) and ` [<taxname>]`
/// where the source names the organism, a nucleotide's (and a RefSeq chromosome's or ncRNA's,
/// and an mRNA's but with `local_annots_only`) the name without white space at its ends.
fn source_title(protein: bool, ids: &BioseqFacts) -> Vec<u8> {
    let taxname = ids.taxname.as_deref().unwrap_or_default();
    let refseq = ids.refseq_chromosome_or_ncrna || (ids.refseq_mrna && !ids.local_annots_only);
    if protein && !refseq {
        let mut title = if ids.local_annots_only {
            Vec::new()
        } else {
            UNNAMED_PROTEIN_PRODUCT.to_vec()
        };
        if !taxname.is_empty() {
            title.extend_from_slice(b" [");
            title.extend_from_slice(taxname);
            title.push(b']');
        }
        return title;
    }
    // C `isspace` of the C locale: space, and tab to carriage return.
    let space = |byte: &u8| matches!(byte, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r');
    let start = taxname.iter().position(|byte| !space(byte));
    let end = taxname.iter().rposition(|byte| !space(byte));
    match (start, end) {
        (Some(start), Some(end)) => taxname[start..=end].to_vec(),
        _ => Vec::new(),
    }
}

/// The entities of `NStr::HtmlDecode`, in NCBI's order.
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4221-4232
/// (and the rows after them, to `{ 9830, "diams" }` at ncbistr.cpp:4508.)
#[rustfmt::skip]
const NCBI_HTML_ENTITIES: [(u32, &str); 280] = [
    (9, "Tab"), (10, "NewLine"), (33, "excl"), (34, "quot"), (35, "num"), (36, "dollar"),
    (37, "percnt"), (38, "amp"), (39, "apos"), (40, "lpar"), (41, "rpar"), (42, "ast"),
    (43, "plus"), (44, "comma"), (46, "period"), (47, "sol"), (58, "colon"), (59, "semi"),
    (60, "lt"), (61, "equals"), (62, "gt"), (63, "quest"), (64, "commat"), (91, "lsqb"),
    (92, "bsol"), (93, "rsqb"), (94, "Hat"), (95, "lowbar"), (96, "grave"), (123, "lcub"),
    (124, "verbar"), (125, "rcub"), (160, "nbsp"), (161, "iexcl"), (162, "cent"),
    (163, "pound"), (164, "curren"), (165, "yen"), (166, "brvbar"), (167, "sect"),
    (168, "uml"), (169, "copy"), (170, "ordf"), (171, "laquo"), (172, "not"), (173, "shy"),
    (174, "reg"), (175, "macr"), (176, "deg"), (177, "plusmn"), (178, "sup2"), (179, "sup3"),
    (180, "acute"), (181, "micro"), (182, "para"), (183, "middot"), (184, "cedil"),
    (185, "sup1"), (186, "ordm"), (187, "raquo"), (188, "frac14"), (189, "frac12"),
    (190, "frac34"), (191, "iquest"), (192, "Agrave"), (193, "Aacute"), (194, "Acirc"),
    (195, "Atilde"), (196, "Auml"), (197, "Aring"), (198, "AElig"), (199, "Ccedil"),
    (200, "Egrave"), (201, "Eacute"), (202, "Ecirc"), (203, "Euml"), (204, "Igrave"),
    (205, "Iacute"), (206, "Icirc"), (207, "Iuml"), (208, "ETH"), (209, "Ntilde"),
    (210, "Ograve"), (211, "Oacute"), (212, "Ocirc"), (213, "Otilde"), (214, "Ouml"),
    (215, "times"), (216, "Oslash"), (217, "Ugrave"), (218, "Uacute"), (219, "Ucirc"),
    (220, "Uuml"), (221, "Yacute"), (222, "THORN"), (223, "szlig"), (224, "agrave"),
    (225, "aacute"), (226, "acirc"), (227, "atilde"), (228, "auml"), (229, "aring"),
    (230, "aelig"), (231, "ccedil"), (232, "egrave"), (233, "eacute"), (234, "ecirc"),
    (235, "euml"), (236, "igrave"), (237, "iacute"), (238, "icirc"), (239, "iuml"),
    (240, "eth"), (241, "ntilde"), (242, "ograve"), (243, "oacute"), (244, "ocirc"),
    (245, "otilde"), (246, "ouml"), (247, "divide"), (248, "oslash"), (249, "ugrave"),
    (250, "uacute"), (251, "ucirc"), (252, "uuml"), (253, "yacute"), (254, "thorn"),
    (255, "yuml"), (338, "OElig"), (339, "oelig"), (352, "Scaron"), (353, "scaron"),
    (376, "Yuml"), (402, "fnof"), (710, "circ"), (732, "tilde"), (913, "Alpha"),
    (914, "Beta"), (915, "Gamma"), (916, "Delta"), (917, "Epsilon"), (918, "Zeta"),
    (919, "Eta"), (920, "Theta"), (921, "Iota"), (922, "Kappa"), (923, "Lambda"), (924, "Mu"),
    (925, "Nu"), (926, "Xi"), (927, "Omicron"), (928, "Pi"), (929, "Rho"), (931, "Sigma"),
    (932, "Tau"), (933, "Upsilon"), (934, "Phi"), (935, "Chi"), (936, "Psi"), (937, "Omega"),
    (945, "alpha"), (946, "beta"), (947, "gamma"), (948, "delta"), (949, "epsilon"),
    (950, "zeta"), (951, "eta"), (952, "theta"), (953, "iota"), (954, "kappa"),
    (955, "lambda"), (956, "mu"), (957, "nu"), (958, "xi"), (959, "omicron"), (960, "pi"),
    (961, "rho"), (962, "sigmaf"), (963, "sigma"), (964, "tau"), (965, "upsilon"),
    (966, "phi"), (967, "chi"), (968, "psi"), (969, "omega"), (977, "thetasym"),
    (978, "upsih"), (982, "piv"), (8194, "ensp"), (8195, "emsp"), (8201, "thinsp"),
    (8204, "zwnj"), (8205, "zwj"), (8206, "lrm"), (8207, "rlm"), (8211, "ndash"),
    (8212, "mdash"), (8216, "lsquo"), (8217, "rsquo"), (8218, "sbquo"), (8220, "ldquo"),
    (8221, "rdquo"), (8222, "bdquo"), (8224, "dagger"), (8225, "Dagger"), (8226, "bull"),
    (8230, "hellip"), (8240, "permil"), (8242, "prime"), (8243, "Prime"), (8249, "lsaquo"),
    (8250, "rsaquo"), (8254, "oline"), (8260, "frasl"), (8364, "euro"), (8472, "weierp"),
    (8465, "image"), (8476, "real"), (8482, "trade"), (8501, "alefsym"), (8592, "larr"),
    (8593, "uarr"), (8594, "rarr"), (8595, "darr"), (8596, "harr"), (8629, "crarr"),
    (8656, "lArr"), (8657, "uArr"), (8658, "rArr"), (8659, "dArr"), (8660, "hArr"),
    (8704, "forall"), (8706, "part"), (8707, "exist"), (8709, "empty"), (8711, "nabla"),
    (8712, "isin"), (8713, "notin"), (8715, "ni"), (8719, "prod"), (8721, "sum"),
    (8722, "minus"), (8727, "lowast"), (8730, "radic"), (8733, "prop"), (8734, "infin"),
    (8736, "ang"), (8743, "and"), (8744, "or"), (8745, "cap"), (8746, "cup"), (8747, "int"),
    (8756, "there4"), (8764, "sim"), (8773, "cong"), (8776, "asymp"), (8800, "ne"),
    (8801, "equiv"), (8804, "le"), (8805, "ge"), (8834, "sub"), (8835, "sup"), (8836, "nsub"),
    (8838, "sube"), (8839, "supe"), (8853, "oplus"), (8855, "otimes"), (8869, "perp"),
    (8901, "sdot"), (8968, "lceil"), (8969, "rceil"), (8970, "lfloor"), (8971, "rfloor"),
    (9001, "lang"), (9002, "rang"), (9674, "loz"), (9824, "spades"), (9827, "clubs"),
    (9829, "hearts"), (9830, "diams"),
];

/// `NStr::HtmlDecode` of `text` with its default `eEncoding_Unknown`: the decoded bytes and
/// whether a character reference was decoded (`fHtmlDec_CharRef_Entity` or
/// `fHtmlDec_CharRef_Numeric`). The text is read as `char`s under the C locale
/// (`isalpha` and the others are ASCII only). A numeric reference wraps modulo 2^32 and is
/// written by `x_AppendChar` without any check (`&#0;` gives a NUL byte, `&#xD800;` a
/// surrogate's bytes).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4523-4618
/// NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:5729-5759 (`CUtf8::x_Append` of a
/// `TUnicodeSymbol` buffer: `x_TCharToUnicodeSymbol` takes the symbol as it is for a
/// 4-byte character, then `x_AppendChar` of `ch` to `u8str`).
fn html_decode(text: &[u8]) -> Result<(Vec<u8>, bool), UnknownEncoding> {
    let encoding = guess_encoding(text);
    if encoding == Encoding::Unknown {
        return Err(UnknownEncoding);
    }
    let mut decoded = Vec::with_capacity(text.len());
    let mut char_refs = false;
    let e = text.len();
    let mut i = 0;
    while i < e {
        let ch = text[i];
        i += 1;
        if i != e && ch == b'&' {
            let mut itmp = i;
            let mut end_of_entity = i;
            let mut ent = text[itmp].is_ascii_alphabetic();
            let mut dec = false;
            if !ent && text[itmp] == b'#' {
                itmp += 1;
                dec = itmp != e && text[itmp].is_ascii_digit();
            }
            let mut hex = false;
            if !dec && itmp != e && matches!(text[itmp], b'x' | b'X') {
                itmp += 1;
                hex = itmp != e && text[itmp].is_ascii_hexdigit();
            }
            let start_of_entity = itmp;
            if itmp != e && (ent || dec || hex) {
                let mut len = 0;
                while len < 16 && itmp != e {
                    let c = text[itmp];
                    if c == b'&' || c == b'#' {
                        break;
                    }
                    if c == b';' {
                        end_of_entity = itmp;
                        break;
                    }
                    ent = ent && c.is_ascii_alphanumeric();
                    dec = dec && c.is_ascii_digit();
                    hex = hex && c.is_ascii_hexdigit();
                    len += 1;
                    itmp += 1;
                }
                if end_of_entity != i && (ent || dec || hex) {
                    let mut uch: u32 = 0;
                    let mut parsed = false;
                    if ent {
                        let entity = &text[start_of_entity..end_of_entity];
                        if let Some(&(symbol, _)) = NCBI_HTML_ENTITIES
                            .iter()
                            .find(|(_, name)| name.as_bytes() == entity)
                        {
                            uch = symbol;
                            parsed = true;
                        }
                    } else {
                        parsed = true;
                        for &digit in &text[start_of_entity..end_of_entity] {
                            let mut ud = u32::from(digit);
                            if dec {
                                uch = uch.wrapping_mul(10).wrapping_add(ud.wrapping_sub(0x30));
                            } else if hex {
                                if digit.is_ascii_digit() {
                                    ud -= u32::from(b'0');
                                } else if (b'a'..=b'f').contains(&digit) {
                                    ud = ud - u32::from(b'a') + 10;
                                } else if (b'A'..=b'F').contains(&digit) {
                                    ud = ud - u32::from(b'A') + 10;
                                }
                                uch = uch.wrapping_mul(16).wrapping_add(ud);
                            }
                        }
                    }
                    if parsed {
                        char_refs = true;
                        append_char(&mut decoded, uch);
                        i = end_of_entity + 1;
                        continue;
                    }
                }
            }
        }
        match encoding {
            Encoding::Utf8 | Encoding::Ascii => decoded.push(ch),
            _ => append_as_utf8(&mut decoded, ch, encoding),
        }
    }
    Ok((decoded, char_refs))
}

/// NCBI's `EEncoding` values that `CUtf8::GuessEncoding` returns.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Encoding {
    Unknown,
    Utf8,
    Ascii,
    Iso8859_1,
    Windows1252,
    Cesu8,
}

/// `CUtf8::GuessEncoding`. The UTF-8 test checks the lead bytes and the number of
/// continuation bytes only (not the ranges of the second byte after E0, ED, F0 and F4).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:6711-6769
fn guess_encoding(src: &[u8]) -> Encoding {
    let mut more = 0;
    let (mut cp1252, mut iso1, mut ascii, mut utf8, mut cesu8) = (true, true, true, true, false);
    for (index, &ch) in src.iter().enumerate() {
        let mut skip = false;
        if more != 0 {
            if eval_next(ch) {
                more -= 1;
                if more == 0 {
                    ascii = false;
                }
                skip = true;
            } else {
                more = 0;
                utf8 = false;
            }
        }
        if ch > 0x7F {
            ascii = false;
            if ch < 0xA0 {
                iso1 = false;
                if matches!(ch, 0x81 | 0x8D | 0x8F | 0x90 | 0x9D) {
                    cp1252 = false;
                }
            }
            if !skip && utf8 && !eval_first(ch, &mut more) {
                utf8 = false;
            }
            if utf8 && !cesu8 && ch == 0xED && src.len() - index > 5 {
                let c1 = src[index + 1];
                let c3 = src[index + 3];
                let c4 = src[index + 4];
                if (c1 & 0xA0) == 0xA0 && c3 == 0xED && (c4 & 0xB0) == 0xB0 {
                    cesu8 = true;
                }
            }
        }
    }
    if more != 0 {
        utf8 = false;
    }
    if ascii {
        Encoding::Ascii
    } else if utf8 {
        if cesu8 {
            Encoding::Cesu8
        } else {
            Encoding::Utf8
        }
    } else if cp1252 {
        if iso1 {
            Encoding::Iso8859_1
        } else {
            Encoding::Windows1252
        }
    } else {
        Encoding::Unknown
    }
}

/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7143-7166
fn eval_first(ch: u8, more: &mut usize) -> bool {
    *more = 0;
    if ch & 0x80 != 0 {
        if ch & 0xE0 == 0xC0 {
            if ch & 0xFE == 0xC0 {
                return false;
            }
            *more = 1;
        } else if ch & 0xF0 == 0xE0 {
            *more = 2;
        } else if ch & 0xF8 == 0xF0 {
            if ch > 0xF4 {
                return false;
            }
            *more = 3;
        } else {
            return false;
        }
    }
    true
}

/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7169-7172
fn eval_next(ch: u8) -> bool {
    ch & 0xC0 == 0x80
}

/// `CUtf8::AsUTF8` of a one-byte `CTempString` and `encoding`, for a byte of a title that
/// `GuessEncoding` did not find to be ASCII or UTF-8: a CESU-8 title keeps the byte (a
/// one-byte string has no six-byte pair), an ISO-8859-1 or Windows-1252 title converts it
/// (`CharToSymbol`, then `x_AppendChar`).
///
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7062-7115
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:6851-6881
fn append_as_utf8(out: &mut Vec<u8>, ch: u8, encoding: Encoding) {
    const CP1252_TABLE: [u32; 32] = [
        0x20AC, 0x003F, 0x201A, 0x0192, 0x201E, 0x2026, 0x2020, 0x2021, 0x02C6, 0x2030, 0x0160,
        0x2039, 0x0152, 0x003F, 0x017D, 0x003F, 0x003F, 0x2018, 0x2019, 0x201C, 0x201D, 0x2022,
        0x2013, 0x2014, 0x02DC, 0x2122, 0x0161, 0x203A, 0x0153, 0x003F, 0x017E, 0x0178,
    ];
    match encoding {
        Encoding::Cesu8 => out.push(ch),
        Encoding::Windows1252 if (0x80..0xA0).contains(&ch) => {
            append_char(out, CP1252_TABLE[usize::from(ch - 0x80)]);
        }
        _ => append_char(out, u32::from(ch)),
    }
}

/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:7040-7060
fn append_char(out: &mut Vec<u8>, ch: u32) {
    if ch < 0x80 {
        out.push(ch as u8);
    } else if ch < 0x800 {
        out.push(((ch >> 6) | 0xC0) as u8);
        out.push(((ch & 0x3F) | 0x80) as u8);
    } else if ch < 0x10000 {
        out.push(((ch >> 12) | 0xE0) as u8);
        out.push((((ch >> 6) & 0x3F) | 0x80) as u8);
        out.push(((ch & 0x3F) | 0x80) as u8);
    } else {
        // A cast to `Uint1` keeps the low byte of a symbol above 0x1FFFFF, as `as u8` does.
        out.push(((ch >> 18) | 0xF0) as u8);
        out.push((((ch >> 12) & 0x3F) | 0x80) as u8);
        out.push((((ch >> 6) & 0x3F) | 0x80) as u8);
        out.push(((ch & 0x3F) | 0x80) as u8);
    }
}

/// `find_last_not_of(chars)` then `erase(pos + 1)`; a text made only of `chars` stays.
fn trim_end_of(text: &mut Vec<u8>, chars: &[u8]) {
    if let Some(pos) = text.iter().rposition(|byte| !chars.contains(byte)) {
        text.truncate(pos + 1);
    }
}

/// `erase(0, find_first_not_of(' '))`: a text of spaces only becomes empty.
fn trim_leading_spaces(text: &mut Vec<u8>) {
    let first = text
        .iter()
        .position(|&byte| byte != b' ')
        .unwrap_or(text.len());
    text.drain(..first);
}

/// NCBI's `x_CleanAndCompress`, byte for byte, except that it stops at the end of the
/// string where NCBI's count of the remaining letters wraps (see below); the flag is
/// whether it wrapped. The bytes are C++ `char`s, which are signed: the last byte is
/// written only when it is between 1 and 0x7F (and not a space), so a title that ends with
/// a byte of 0x80 or more loses that byte. In the middle of the title the sign does not
/// matter, because every pattern ends with an ASCII byte.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:219-312
fn clean_and_compress(input: &[u8], is_protein: bool) -> (Vec<u8>, bool) {
    let mut start = 0;
    let mut end = input.len();
    while start < end && input[start] == b' ' {
        start += 1;
    }
    while end > start && input[end - 1] == b' ' {
        end -= 1;
    }
    let text = &input[start..end];
    let mut out = Vec::with_capacity(text.len());
    if text.is_empty() {
        return (out, false);
    }
    // Past the end reads the terminating NUL of the C++ string.
    let at = |index: usize| text.get(index).copied().unwrap_or(0);
    let mut index = 0;
    let mut left = text.len();
    let mut wrapped = false;
    let mut curr = at(index);
    index += 1;
    left -= 1;
    let mut two_chars = u16::from(curr);
    while left > 0 {
        let mut next = at(index);
        index += 1;
        two_chars = (two_chars << 8) | u16::from(next);
        match two_chars.to_be_bytes() {
            [b',', b','] => {
                out.push(curr);
                next = b' ';
            }
            [b' ', b' '] | [b' ', b')'] => {}
            [b'(', b' '] => {
                next = curr;
                two_chars = u16::from(curr);
            }
            [b' ', b','] | [b' ', b';'] => {
                out.push(next);
                next = curr;
                two_chars = u16::from(curr);
            }
            [separator @ (b',' | b';'), b' '] => {
                out.push(curr);
                out.push(b' ');
                while next == b' ' || next == separator {
                    next = at(index);
                    index += 1;
                    // Where NCBI's `left` (a size_t) wraps, for a run of spaces and
                    // separators that reaches the end of the string, its loop runs past the
                    // string (it reads the terminating NUL and on, and crashes). LOSAT stops
                    // at the end of the string (approved exception 2 of
                    // PD-LOSAT-NCBI-DEFECTS): `next` is the NUL and nothing more is written.
                    wrapped |= left == 0;
                    left = left.saturating_sub(1);
                }
                two_chars = u16::from(next);
            }
            _ => out.push(curr),
        }
        curr = next;
        wrapped |= left == 0;
        left = left.saturating_sub(1);
    }
    // `curr > 0` of a signed `char`.
    if (1..0x80).contains(&curr) && curr != b' ' {
        out.push(curr);
    }
    if is_protein {
        replace_in_place(&mut out, b". [", b" [");
        replace_in_place(&mut out, b", [", b" [");
    }
    (out, wrapped)
}

#[cfg(test)]
#[path = "defline_tests.rs"]
mod tests;
