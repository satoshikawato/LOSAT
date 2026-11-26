//! NCBI `CDeflineGenerator::GenerateDefline` for FASTA, BLAST database, and full ASN.1 Bioseq
//! records, and `CShowBlastDefline::GetSeqIdList` tabular subject IDs.
//!
//! NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3928-4107
//!
//! All callers (search reports, `blastdbcmd`, and `makeblastdb -input_type asn1_txt|asn1_bin`)
//! route through [`generate_defline_core`], the single owner of `CDeflineGenerator`.
//! For some titles of commas, semicolons, tildes and spaces, NCBI's `x_CleanAndCompress` lets
//! its count of the remaining letters wrap and reads past the end of the string (and crashes);
//! `clean_and_compress` stops at the end of the string there (approved exception 2 of
//! PD-LOSAT-NCBI-DEFECTS).

pub(crate) mod flags;
mod html;
mod prefix_suffix;
pub(crate) mod titles;

use crate::seqid::accguide::{
    identify_accession, parse_flags, SeqIdWarning, CHROMOSOME, CON, DIVISION_MASK, F_MASTER, F_NUC,
    F_PREDICTED, F_PROT, MRNA, NCRNA, SEQID_REFSEQ, TSA, UNIQUE, WGS,
};
use crate::seqid::model::{ObjectId, ParsedSeqId, PatentNumber};

use flags::{mol_info, DeflineFeatures, DeflineInput};
use html::{clean_and_compress, contains_bytes, html_decode, trim_end_of, trim_leading_spaces};
#[cfg(test)]
use prefix_suffix::adjust_protein_title_suffix;
use prefix_suffix::{adjust_protein_title_suffix_full, set_prefix, set_suffix, TPA_PREFIXES};
use titles::{
    organelle_name, set_title_from_biosrc, set_title_from_gpipe, set_title_from_map,
    set_title_from_nc, set_title_from_nm, set_title_from_nr, set_title_from_patent,
    set_title_from_pdb, set_title_from_protein, set_title_from_seg_seq, set_title_from_wgs,
};

/// The title that NCBI generates for a protein record without one.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2508-2509
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

/// `CDeflineGenerator::GenerateDefline` of a record from FASTA with the title `title`: for the
/// alignment heading (`leave_prefix_suffix` false, flags 0) or for the description table (true,
/// `fLeavePrefixSuffix`), of a protein or a nucleotide record.
pub fn generate_defline(
    title: &[u8],
    protein: bool,
    leave_prefix_suffix: bool,
) -> Result<Defline, UnknownEncoding> {
    generate_bioseq_defline(title, protein, leave_prefix_suffix, &BioseqFacts::default())
}

/// What `CDeflineGenerator::x_SetFlags` takes from the Seq-ids and source of a database Bioseq
/// (none for a record read from FASTA, whose id is local).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:745-829,1004-1011,1264-1274
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BioseqFacts {
    /// `m_ThirdParty`: a `tpg`, `tpe` or `tpd` id (the prefix `TPA: `).
    pub third_party: bool,
    /// `m_IsNC` or `m_IsNR`: an accession of a RefSeq chromosome or ncRNA.
    pub refseq_chromosome_or_ncrna: bool,
    /// `m_IsNM`: an accession of a RefSeq mRNA (`NM_`, `XM_`).
    pub refseq_mrna: bool,
    /// `m_IsNZ`: an accession of a RefSeq contig (`NT_`, `NZ_...`).
    pub refseq_contig: bool,
    /// `m_IsWP`: an accession of a RefSeq unique protein (`WP_...`).
    pub refseq_unique_prot: bool,
    /// `m_IsPDB`, `m_PDBChainID` and `m_PDBChain`: whether a `pdb` id was seen and its chain.
    pub pdb: Option<(Vec<u8>, i32)>,
    /// `m_IsPatent`, `m_PatentSequence`, `m_PatentCountry` and `m_PatentNumber`.
    pub patent: Option<(i32, Vec<u8>, Vec<u8>)>,
    /// `m_Taxname`: the scientific name of the Bioseq's first source, `None` without a source
    /// or a name (never empty).
    pub taxname: Option<std::sync::Arc<[u8]>>,
    /// `m_Source` is set: some definition line has a positive taxid, with or without a name.
    pub has_source: bool,
    /// `m_LocalAnnotsOnly` (`fLocalAnnotsOnly`, set by blastdbcmd's `CFastaOstream` with
    /// `fNoExpensiveOps`): an empty protein title is not `unnamed protein product`.
    pub local_annots_only: bool,
}

/// The name of the organism of a database sequence's first source (`BioseqFacts::taxname`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum SourceTaxname {
    #[default]
    Unnamed,
    Named(std::sync::Arc<[u8]>),
    Undecodable {
        name: std::sync::Arc<[u8]>,
        seq_id: std::sync::Arc<[u8]>,
    },
    Nameless,
    Refused(&'static str),
}

/// `SourceTaxname::Unnamed`, for a subject without a source.
pub static UNNAMED_SOURCE: SourceTaxname = SourceTaxname::Unnamed;

impl SourceTaxname {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1165-1169
    // NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1271-1273
    pub fn of_names(names: &crate::seqdb::taxdb::SourceNames) -> SourceTaxname {
        match names.taxname() {
            Some(name) => SourceTaxname::Named(std::sync::Arc::clone(name)),
            None => SourceTaxname::Unnamed,
        }
    }

    pub fn name(&self) -> Result<Option<std::sync::Arc<[u8]>>, &'static str> {
        match self {
            SourceTaxname::Unnamed | SourceTaxname::Nameless => Ok(None),
            SourceTaxname::Named(name) | SourceTaxname::Undecodable { name, .. } => {
                Ok(Some(std::sync::Arc::clone(name)))
            }
            SourceTaxname::Refused(reason) => Err(reason),
        }
    }

    pub fn has_source(&self) -> bool {
        matches!(
            self,
            SourceTaxname::Named(_) | SourceTaxname::Undecodable { .. } | SourceTaxname::Nameless
        )
    }

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

    pub fn undecodable_id(&self) -> Option<&[u8]> {
        match self {
            SourceTaxname::Undecodable { seq_id, .. } => Some(seq_id),
            _ => None,
        }
    }
}

/// Seq-id classification shared between lightweight `BioseqFacts` and full ASN.1 `DeflineInput`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:745-838
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct FullSeqIdFlags {
    pub(crate) is_nc: bool,
    pub(crate) is_nm: bool,
    pub(crate) is_nr: bool,
    pub(crate) is_nz: bool,
    pub(crate) is_wp: bool,
    pub(crate) third_party: bool,
    pub(crate) wgs_master: bool,
    pub(crate) tsa_master: bool,
    pub(crate) general_str: Vec<u8>,
    pub(crate) general_id: i32,
    pub(crate) is_pdb: bool,
    pub(crate) pdb_chain_id: Vec<u8>,
    pub(crate) pdb_chain: i32,
    pub(crate) is_patent: bool,
    pub(crate) patent_sequence: i32,
    pub(crate) patent_country: Vec<u8>,
    pub(crate) patent_number: Vec<u8>,
}

/// Classifies all Seq-ids of a Bioseq in NCBI's exact order, posting `IdentifyAccession` warnings.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:745-838
pub(crate) fn classify_seq_ids_full(
    ids: &[ParsedSeqId],
    is_virtual: bool,
    on_warning: &mut dyn FnMut(SeqIdWarning),
) -> FullSeqIdFlags {
    let mut out = FullSeqIdFlags::default();
    for id in ids {
        match id {
            ParsedSeqId::Tpg(text)
            | ParsedSeqId::Tpe(text)
            | ParsedSeqId::Tpd(text)
            | ParsedSeqId::Other(text)
            | ParsedSeqId::Genbank(text)
            | ParsedSeqId::Embl(text)
            | ParsedSeqId::Ddbj(text) => {
                if matches!(
                    id,
                    ParsedSeqId::Tpg(_) | ParsedSeqId::Tpe(_) | ParsedSeqId::Tpd(_)
                ) {
                    out.third_party = true;
                }
                if let Some(acc) = text.accession.as_deref() {
                    let info = identify_accession(acc, parse_flags::DEFAULT, on_warning);
                    let div = info & DIVISION_MASK;
                    if div == WGS {
                        if (info & F_MASTER) != 0 {
                            out.wgs_master = true;
                        }
                    } else if div == TSA {
                        if (info & F_MASTER) != 0 && is_virtual {
                            out.tsa_master = true;
                        }
                    } else if info == (SEQID_REFSEQ | CHROMOSOME | F_NUC) {
                        out.is_nc = true;
                    } else if info == (SEQID_REFSEQ | MRNA | F_NUC)
                        || info == (SEQID_REFSEQ | MRNA | F_NUC | F_PREDICTED)
                    {
                        out.is_nm = true;
                    } else if info == (SEQID_REFSEQ | NCRNA | F_NUC) {
                        out.is_nr = true;
                    } else if info == (SEQID_REFSEQ | CON | F_NUC) {
                        out.is_nz = true;
                    } else if info == (SEQID_REFSEQ | UNIQUE | F_PROT) {
                        out.is_wp = true;
                    }
                }
            }
            ParsedSeqId::General(dbtag) => {
                if !flags::is_skippable_db(&dbtag.db) {
                    match &dbtag.tag {
                        ObjectId::Str(s) => out.general_str = s.clone(),
                        ObjectId::Id(n) => out.general_id = *n,
                    }
                }
            }
            ParsedSeqId::Pdb(pdb) => {
                out.is_pdb = true;
                if let Some(chain_id) = &pdb.chain_id {
                    out.pdb_chain_id = chain_id.clone();
                } else if let Some(chain) = pdb.chain {
                    out.pdb_chain = chain;
                }
            }
            ParsedSeqId::Patent(patent) => {
                out.is_patent = true;
                out.patent_sequence = patent.seqid;
                out.patent_country = patent.cit.country.clone();
                out.patent_number = match &patent.cit.id {
                    PatentNumber::Number(n) | PatentNumber::AppNumber(n) => n.clone(),
                };
            }
            _ => {}
        }
    }
    out
}

impl BioseqFacts {
    pub fn of_source(taxname: Option<std::sync::Arc<[u8]>>) -> Self {
        Self {
            has_source: taxname.is_some(),
            taxname,
            ..Self::default()
        }
    }

    fn to_defline_input<'a>(
        &'a self,
        title: &'a [u8],
        protein: bool,
        leave_prefix_suffix: bool,
    ) -> DeflineInput<'a> {
        let (is_pdb, pdb_chain_id, pdb_chain) = match &self.pdb {
            Some((chain_id, chain)) => (true, chain_id.as_slice(), *chain),
            None => (false, &[][..], 0),
        };
        let (is_patent, patent_sequence, patent_country, patent_number) = match &self.patent {
            Some((seq, country, num)) => (true, *seq, country.as_slice(), num.as_slice()),
            None => (false, 0, &[][..], &[][..]),
        };
        DeflineInput {
            leave_prefix_suffix,
            local_annots_only: self.local_annots_only,
            is_aa: protein,
            is_na: !protein,
            is_nc: self.refseq_chromosome_or_ncrna,
            is_nm: self.refseq_mrna,
            is_nz: self.refseq_contig,
            is_wp: self.refseq_unique_prot,
            third_party: self.third_party,
            main_title: title,
            is_pdb,
            pdb_chain_id,
            pdb_chain,
            use_pdb_compound_for_defline: is_pdb,
            is_patent,
            patent_sequence,
            patent_country,
            patent_number,
            has_source: self.has_source,
            taxname: self.taxname.as_deref().unwrap_or_default(),
            ..DeflineInput::default()
        }
    }
}

/// `BioseqFacts` of the Seq-ids `ids` of a Bioseq, posting any `CSeq_id::IdentifyAccession`
/// warnings through `on_warning`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:745-832
pub(crate) fn bioseq_facts_with_warn(
    ids: &[ParsedSeqId],
    on_warning: &mut dyn FnMut(SeqIdWarning),
) -> BioseqFacts {
    let full = classify_seq_ids_full(ids, false, on_warning);
    BioseqFacts {
        third_party: full.third_party,
        refseq_chromosome_or_ncrna: full.is_nc || full.is_nr,
        refseq_mrna: full.is_nm,
        refseq_contig: full.is_nz,
        refseq_unique_prot: full.is_wp,
        pdb: full.is_pdb.then_some((full.pdb_chain_id, full.pdb_chain)),
        patent: full.is_patent.then_some((
            full.patent_sequence,
            full.patent_country,
            full.patent_number,
        )),
        taxname: None,
        has_source: false,
        local_annots_only: false,
    }
}

/// `BioseqFacts` of the Seq-ids `ids` of a Bioseq, discarding accession-guide warnings.
pub(crate) fn bioseq_facts(ids: &[ParsedSeqId]) -> BioseqFacts {
    bioseq_facts_with_warn(ids, &mut |_| {})
}

/// Builds the empty-title fallback via `x_SetBioSrc` and the specialized RefSeq/protein/seg/WGS/
/// map/GPipe/BioSrc title generators, returning `(main_title, organelle, multispecies, append_complete)`.
fn build_empty_source_title(
    inp: &DeflineInput<'_>,
    feats: Option<&dyn DeflineFeatures>,
) -> (Vec<u8>, &'static [u8], bool, bool) {
    let organelle = organelle_name(inp);
    let multispecies = inp.multispecies_from_source;
    let has_clone = inp.has_clone || feats.is_some_and(DeflineFeatures::has_biosrc_clone_feature);

    let mut main_title = if inp.is_nc {
        set_title_from_nc(inp, organelle)
    } else if inp.is_nm && !inp.local_annots_only {
        set_title_from_nm(inp, feats)
    } else if inp.is_nr {
        set_title_from_nr(inp, feats)
    } else if inp.is_aa {
        set_title_from_protein(inp, feats)
    } else if inp.is_seg && !inp.is_est_sts_gss {
        set_title_from_seg_seq(inp, feats)
    } else if inp.is_tsa || (inp.is_wgs && !inp.wgs_master) || (inp.is_tls && !inp.tls_master) {
        set_title_from_wgs(inp, has_clone)
    } else if inp.is_map {
        set_title_from_map(inp)
    } else {
        Vec::new()
    };

    if main_title.is_empty() && inp.gpipe_mode {
        main_title = set_title_from_gpipe(inp, organelle, has_clone);
    }

    let mut append_complete = false;
    if main_title.is_empty() {
        main_title = set_title_from_biosrc(inp, organelle, has_clone);
        if inp.mi_completeness == mol_info::COMPLETENESS_COMPLETE && !main_title.is_empty() {
            append_complete = true;
        }
    }
    (main_title, organelle, multispecies, append_complete)
}

/// Builds the main title prior to `NStr::HtmlDecode`, returning
/// `(main_title, organelle, multispecies, capitalize, append_complete)`.
fn build_main_title_before_decoding(
    inp: &DeflineInput<'_>,
    feats: Option<&dyn DeflineFeatures>,
) -> (Vec<u8>, &'static [u8], bool, bool, bool) {
    let mut capitalize = true;
    let mut append_complete = false;
    let mut main_title = if inp.reconstruct {
        Vec::new()
    } else {
        let mut t = inp.main_title.to_vec();
        trim_end_of(&mut t, b".,;~ ");
        if !t.is_empty() {
            capitalize = false;
        }
        t
    };

    if inp.is_aa && !main_title.is_empty() && !inp.leave_prefix_suffix {
        adjust_protein_title_suffix_full(&mut main_title, inp, feats);
    }

    let mut organelle: &'static [u8] = b"";
    let mut multispecies = false;

    if main_title.is_empty() {
        if inp.is_pdb {
            main_title = set_title_from_pdb(inp);
        } else if inp.is_patent {
            main_title = set_title_from_patent(inp);
        }

        if main_title.is_empty() {
            (main_title, organelle, multispecies, append_complete) =
                build_empty_source_title(inp, feats);
        }
    }

    if !inp.leave_prefix_suffix {
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
    (
        main_title,
        organelle,
        multispecies,
        capitalize,
        append_complete,
    )
}

/// Single-owner `CDeflineGenerator::GenerateDefline` implementation used by every caller
/// (FASTA/DB search reports, `blastdbcmd`, and `makeblastdb -input_type asn1_txt|asn1_bin`).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3928-4107
pub(crate) fn generate_defline_core(
    inp: &DeflineInput<'_>,
    feats: Option<&dyn DeflineFeatures>,
) -> Result<Defline, UnknownEncoding> {
    let (main_title, organelle, multispecies, capitalize, append_complete) =
        build_main_title_before_decoding(inp, feats);

    let (mut decoded, _) = html_decode(&main_title)?;
    trim_end_of(&mut decoded, b",;~ ");

    // NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:4066-4093
    // `x_SetPrefix` and `x_SetSuffix` inspect pre-decoded `m_MainTitle`, while `decoded` is used
    // only to assemble `penult`. Avoid cloning `decoded` when no prefix/suffix/MAG applies.
    let penult = if !inp.leave_prefix_suffix {
        let prefix = set_prefix(inp, &main_title, multispecies);
        let suffix = set_suffix(inp, &main_title, organelle, append_complete);
        let mag: &[u8] = if !inp.metagenome_source.is_empty() {
            if prefix.is_empty() {
                b"MAG: "
            } else {
                b"MAG "
            }
        } else {
            b""
        };
        if mag.is_empty() && prefix.is_empty() && suffix.is_empty() {
            decoded
        } else {
            let mut out =
                Vec::with_capacity(mag.len() + prefix.len() + decoded.len() + suffix.len());
            out.extend_from_slice(mag);
            out.extend_from_slice(prefix);
            out.extend_from_slice(&decoded);
            out.extend_from_slice(&suffix);
            out
        }
    } else {
        decoded
    };

    let (mut text, reads_past_end) = clean_and_compress(&penult, inp.is_aa);

    if !inp.is_pdb && !inp.is_patent && !inp.is_aa && !inp.is_seg && capitalize {
        if let Some(first) = text.first_mut() {
            if first.is_ascii_lowercase() {
                first.make_ascii_uppercase();
            }
        }
    }

    Ok(Defline {
        text,
        reads_past_end,
    })
}

/// `GenerateDefline` of a Bioseq with the title `title` and the facts `ids` (`BioseqFacts`).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3928-4107
pub fn generate_bioseq_defline(
    title: &[u8],
    protein: bool,
    leave_prefix_suffix: bool,
    ids: &BioseqFacts,
) -> Result<Defline, UnknownEncoding> {
    let inp = ids.to_defline_input(title, protein, leave_prefix_suffix);
    generate_defline_core(&inp, None)
}

/// Whether `generate_bioseq_defline` makes the title of `title` and `ids` from the source of
/// the Bioseq.
pub fn bioseq_title_from_source(title: &[u8], ids: impl FnOnce() -> BioseqFacts) -> bool {
    let main_title = match title.iter().rposition(|byte| !b".,;~ ".contains(byte)) {
        Some(pos) => &title[..=pos],
        None => title,
    };
    if !main_title.is_empty() {
        return false;
    }
    let facts = ids();
    let inp = facts.to_defline_input(b"", false, false);
    if inp.is_pdb {
        set_title_from_pdb(&inp).is_empty()
    } else if inp.is_patent {
        set_title_from_patent(&inp).is_empty()
    } else {
        true
    }
}

/// The subject ID of the tabular formats for a subject with local ID `local_id` and title `title`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:787-789
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:482-496,522-532
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:218-242
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
    if replaced.starts_with(b"Subject_") || contains_bytes(replaced, b"lcl|Subject_") {
        let defline = generate_defline(title, protein, false)?.text;
        let token = first_word(&defline);
        if !token.is_empty() {
            return Ok(token.to_vec());
        }
    }
    Ok(replaced.to_vec())
}

/// The words of LOSAT's refusal of a protein title of punctuation that NCBI BLAST+'s
/// `x_CleanAndCompress` reads past the end of.
macro_rules! punctuation_title {
    () => {
        "of punctuation that NCBI BLAST+'s x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes)"
    };
}
pub(crate) use punctuation_title;

/// The subject ID of the tabular formats of a database sequence without parsed Seq-ids.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:218-242
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:522-532
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

/// What `GetSeqIdList` puts in the place of a database sequence's general id.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SubjectWord {
    Word(Vec<u8>),
    General,
    Unknown,
}

fn first_word(text: &[u8]) -> &[u8] {
    let end = text
        .iter()
        .position(|&byte| byte == b' ')
        .unwrap_or(text.len());
    &text[..end]
}

/// Whether `adjust_protein_title_suffix` can change `title` with some source.
pub fn protein_title_suffix_adjustable(title: &[u8]) -> bool {
    title
        .iter()
        .rposition(|byte| !b".,;~ ".contains(byte))
        .is_some_and(|end| title[end] == b']')
}

/// Whether the first word of the generated protein title can depend on the Bioseq's source.
pub fn protein_title_word_needs_source(title: &[u8]) -> bool {
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
    let Some(space) = main.iter().position(|&byte| !PLAIN[byte as usize]) else {
        return true;
    };
    let next = main[space + 1..].first().copied().unwrap_or(b' ');
    !(space > 0
        && main[space] == b' '
        && !main[..space].eq_ignore_ascii_case(b"MAG")
        && (PLAIN[next as usize] || next == b'['))
}

#[cfg(test)]
#[path = "defline_tests.rs"]
mod tests;
