//! Full `CDeflineGenerator` flag state, NCBI ASN.1 enum constants, and lazy feature-lookup trait.
//!
//! NCBI reference (598d8ae6): c++/include/objmgr/util/create_defline.hpp:234-373
//! NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:612-1434

use super::html::contains_bytes;

/// `CSeq_inst::ERepr` values from `seq.asn`.
pub(crate) mod seq_repr {
    pub(crate) const VIRTUAL: i32 = 1;
    pub(crate) const SEG: i32 = 3;
    pub(crate) const MAP: i32 = 7;
    pub(crate) const DELTA: i32 = 8;
}

/// `CSeq_inst::ETopology` values from `seq.asn`.
pub(crate) mod topology {
    pub(crate) const LINEAR: i32 = 1;
    pub(crate) const CIRCULAR: i32 = 2;
}

/// `CMolInfo` enum values from `seq.asn`.
pub(crate) mod mol_info {
    pub(crate) const BIOMOL_GENOMIC: i32 = 1;
    pub(crate) const BIOMOL_PRE_RNA: i32 = 2;
    pub(crate) const BIOMOL_MRNA: i32 = 3;
    pub(crate) const BIOMOL_RRNA: i32 = 4;
    pub(crate) const BIOMOL_TRNA: i32 = 5;
    pub(crate) const BIOMOL_SNRNA: i32 = 6;
    pub(crate) const BIOMOL_SCRNA: i32 = 7;
    pub(crate) const BIOMOL_OTHER_GENETIC: i32 = 9;
    pub(crate) const BIOMOL_CRNA: i32 = 11;
    pub(crate) const BIOMOL_SNORNA: i32 = 12;
    pub(crate) const BIOMOL_TRANSCRIBED_RNA: i32 = 13;
    pub(crate) const BIOMOL_NCRNA: i32 = 14;
    pub(crate) const BIOMOL_TMRNA: i32 = 15;

    pub(crate) const TECH_EST: i32 = 2;
    pub(crate) const TECH_STS: i32 = 3;
    pub(crate) const TECH_SURVEY: i32 = 4;
    pub(crate) const TECH_HTGS_1: i32 = 14;
    pub(crate) const TECH_HTGS_2: i32 = 15;
    pub(crate) const TECH_HTGS_3: i32 = 16;
    pub(crate) const TECH_HTGS_0: i32 = 18;
    pub(crate) const TECH_WGS: i32 = 20;
    pub(crate) const TECH_TSA: i32 = 23;
    pub(crate) const TECH_TARGETED: i32 = 24;

    pub(crate) const COMPLETENESS_COMPLETE: i32 = 1;
    pub(crate) const COMPLETENESS_PARTIAL: i32 = 2;
    pub(crate) const COMPLETENESS_NO_LEFT: i32 = 3;
    pub(crate) const COMPLETENESS_NO_RIGHT: i32 = 4;
    pub(crate) const COMPLETENESS_NO_ENDS: i32 = 5;
}

/// `CBioSource::EGenome` values from `seqfeat.asn`.
pub(crate) mod genome {
    pub(crate) const CHLOROPLAST: i32 = 2;
    pub(crate) const CHROMOPLAST: i32 = 3;
    pub(crate) const KINETOPLAST: i32 = 4;
    pub(crate) const MITOCHONDRION: i32 = 5;
    pub(crate) const PLASTID: i32 = 6;
    pub(crate) const MACRONUCLEAR: i32 = 7;
    pub(crate) const EXTRACHROM: i32 = 8;
    pub(crate) const PLASMID: i32 = 9;
    pub(crate) const CYANELLE: i32 = 12;
    pub(crate) const PROVIRAL: i32 = 13;
    pub(crate) const VIRION: i32 = 14;
    pub(crate) const NUCLEOMORPH: i32 = 15;
    pub(crate) const APICOPLAST: i32 = 16;
    pub(crate) const LEUCOPLAST: i32 = 17;
    pub(crate) const PROPLASTID: i32 = 18;
    pub(crate) const ENDOGENOUS_VIRUS: i32 = 19;
    pub(crate) const HYDROGENOSOME: i32 = 20;
    pub(crate) const CHROMOSOME: i32 = 21;
    pub(crate) const CHROMATOPHORE: i32 = 22;
}

/// `CSubSource::ESubtype` values used by `CDeflineGenerator`.
pub(crate) mod subsource {
    pub(crate) const CHROMOSOME: i32 = 1;
    pub(crate) const MAP: i32 = 2;
    pub(crate) const CLONE: i32 = 3;
    pub(crate) const PLASMID_NAME: i32 = 19;
    pub(crate) const SEGMENT: i32 = 24;
    pub(crate) const LINKAGE_GROUP: i32 = 39;
}

/// `COrgMod::ESubtype` values used by `CDeflineGenerator`.
pub(crate) mod orgmod {
    pub(crate) const STRAIN: i32 = 2;
    pub(crate) const SUBSTRAIN: i32 = 3;
    pub(crate) const CULTIVAR: i32 = 10;
    pub(crate) const ISOLATE: i32 = 17;
    pub(crate) const SPECIMEN_VOUCHER: i32 = 23;
    pub(crate) const BREED: i32 = 31;
    pub(crate) const METAGENOME_SOURCE: i32 = 37;
}

/// Protein feature data returned by `DeflineFeatures::longest_protein`.
#[derive(Clone, Debug, Default)]
pub(crate) struct ProteinFeatInfo {
    pub(crate) names: Vec<Vec<u8>>,
    pub(crate) desc: Option<Vec<u8>>,
    pub(crate) first_activity: Option<Vec<u8>>,
}

/// CDS data returned by `DeflineFeatures::cds_for_product`.
#[derive(Clone, Debug, Default)]
pub(crate) struct CdsForProductInfo {
    pub(crate) has_low_quality_exception: bool,
}

/// Gene data returned by `DeflineFeatures::gene_via_cds`.
#[derive(Clone, Debug, Default)]
pub(crate) struct GeneViaCdsInfo {
    pub(crate) locus: Option<Vec<u8>>,
    pub(crate) first_syn: Option<Vec<u8>>,
    pub(crate) desc: Option<Vec<u8>>,
    pub(crate) locus_tag: Option<Vec<u8>>,
}

/// Segmented-sequence CDS info returned by `DeflineFeatures::segseq_info_via_cds`.
#[derive(Clone, Debug, Default)]
pub(crate) struct SegSeqCdsInfo {
    pub(crate) product: Vec<u8>,
    pub(crate) locus: Vec<u8>,
    pub(crate) is_partial: bool,
}

/// Lazy feature and annotation callbacks used by `CDeflineGenerator` when generating titles from
/// full ASN.1 `Seq-entry` objects. Database/FASTA callers pass `None`.
pub(crate) trait DeflineFeatures {
    /// Whether any `Biosrc` feature on the Bioseq carries a `SubSource` of subtype `clone` with a
    /// set name (`x_SetBioSrc`, create_defline.cpp:1409-1433).
    fn has_biosrc_clone_feature(&self) -> bool {
        false
    }

    /// For `x_SetTitleFromNM`: when the Bioseq has exactly one `Gene` and one `Cdregion` feature,
    /// returns their `(cds_label, gene_label)` computed with `feature::GetLabel` (`fFGL_Content`).
    fn nm_single_gene_and_cds(&self) -> Option<(Vec<u8>, Vec<u8>)> {
        None
    }

    /// For `x_SetTitleFromNR`: label of the first `Gene` feature on the Bioseq (`fFGL_Content`).
    fn nr_first_gene_label(&self) -> Option<Vec<u8>> {
        None
    }

    /// For `x_SetTitleFromProtein`: `x_GetLongestProtein`.
    fn longest_protein(&self) -> Option<ProteinFeatInfo> {
        None
    }

    /// For `x_SetTitleFromProtein`: `GetMappedCDSForProduct` and `x_CDShasLowQualityException`.
    fn cds_for_product(&self) -> Option<CdsForProductInfo> {
        None
    }

    /// For `x_SetTitleFromProtein`: `x_GetGeneRefViaCDS` of the CDS that `cds_for_product`
    /// finds, read only where the title needs it, as NCBI does.
    /// `x_SetSuffix`'s walk of the sequence map of a `delta` sequence (HTGS-1/2 or `UNORDERED`),
    /// which looks up the length of each `whole` segment and may throw.
    fn suffix_length_walk(&self) {}

    fn gene_via_cds(&self) -> Option<GeneViaCdsInfo> {
        None
    }

    /// For `x_SetTitleFromProtein` and `x_AdjustProteinTitleSuffix`: `x_GetSourceFeatViaCDS`.
    /// Returns `Some(taxname)` when an overlapping `Biosrc` feature with a set `taxname` exists.
    fn source_taxname_via_cds(&self) -> Option<Vec<u8>> {
        None
    }

    /// For `x_SetTitleFromSegSeq`: `x_GetSegSeqInfoViaCDS`.
    fn segseq_info_via_cds(&self) -> Option<SegSeqCdsInfo> {
        None
    }
}

/// Complete borrowed inputs for `CDeflineGenerator::GenerateDefline`, populated either from a
/// lightweight `BioseqFacts` (zero allocations) or from a typed ASN.1 Bioseq and its ancestors.
#[derive(Clone, Debug, Default)]
pub(crate) struct DeflineInput<'a> {
    // User flags
    pub(crate) reconstruct: bool,
    pub(crate) all_prot_names: bool,
    pub(crate) local_annots_only: bool,
    pub(crate) gpipe_mode: bool,
    pub(crate) omit_taxonomic_name: bool,
    pub(crate) leave_prefix_suffix: bool,
    pub(crate) fasta_format: bool,

    // Seq-inst fields
    pub(crate) is_na: bool,
    pub(crate) is_aa: bool,
    pub(crate) topology: i32,
    pub(crate) length: u64,
    pub(crate) is_seg: bool,
    pub(crate) is_delta: bool,
    pub(crate) is_virtual: bool,
    pub(crate) is_map: bool,

    // Seq-id fields
    pub(crate) is_nc: bool,
    pub(crate) is_nm: bool,
    pub(crate) is_nr: bool,
    pub(crate) is_nz: bool,
    pub(crate) is_patent: bool,
    pub(crate) is_pdb: bool,
    pub(crate) is_wp: bool,
    pub(crate) third_party: bool,
    pub(crate) wgs_master: bool,
    pub(crate) tsa_master: bool,
    pub(crate) tls_master: bool,

    pub(crate) main_title: &'a [u8],
    pub(crate) general_str: &'a [u8],
    pub(crate) general_id: i32,
    pub(crate) patent_country: &'a [u8],
    pub(crate) patent_number: &'a [u8],
    pub(crate) patent_sequence: i32,
    pub(crate) pdb_chain: i32,
    pub(crate) pdb_chain_id: &'a [u8],

    // MolInfo fields
    pub(crate) mi_biomol: i32,
    pub(crate) mi_tech: i32,
    pub(crate) mi_completeness: i32,
    pub(crate) htg_tech: bool,
    pub(crate) htgs_unfinished: bool,
    pub(crate) is_tls: bool,
    pub(crate) is_tsa: bool,
    pub(crate) is_wgs: bool,
    pub(crate) is_est_sts_gss: bool,

    // GenBank / EMBL keyword flags
    pub(crate) htgs_cancelled: bool,
    pub(crate) htgs_draft: bool,
    pub(crate) htgs_pooled: bool,
    pub(crate) tpa_exp: bool,
    pub(crate) tpa_inf: bool,
    pub(crate) tpa_reasm: bool,
    pub(crate) unordered: bool,

    // PDB block fields
    pub(crate) pdb_compound: &'a [u8],
    pub(crate) use_pdb_compound_for_defline: bool,

    // BioSource fields
    pub(crate) has_source: bool,
    pub(crate) taxname: &'a [u8],
    pub(crate) genus: &'a [u8],
    pub(crate) species: &'a [u8],
    pub(crate) multispecies_from_source: bool,
    pub(crate) genome: i32,
    pub(crate) is_plasmid: bool,
    pub(crate) is_chromosome: bool,
    pub(crate) first_super_kingdom: &'a [u8],
    pub(crate) second_super_kingdom: &'a [u8],
    pub(crate) is_cross_kingdom: bool,

    // SubSource fields
    pub(crate) chromosome: &'a [u8],
    pub(crate) linkage_group: &'a [u8],
    pub(crate) clone: &'a [u8],
    pub(crate) has_clone: bool,
    pub(crate) map: &'a [u8],
    pub(crate) plasmid: &'a [u8],
    pub(crate) segment: &'a [u8],

    // OrgMod fields
    pub(crate) breed: &'a [u8],
    pub(crate) cultivar: &'a [u8],
    pub(crate) specimen_voucher: &'a [u8],
    pub(crate) isolate: &'a [u8],
    pub(crate) strain: &'a [u8],
    pub(crate) substrain: &'a [u8],
    pub(crate) metagenome_source: &'a [u8],

    // User object fields
    pub(crate) is_unverified: bool,
    pub(crate) unverified_prefix: &'static [u8],
    pub(crate) is_unreviewed: bool,
    pub(crate) unreviewed_prefix: &'static [u8],
    pub(crate) targeted_locus: &'a [u8],

    // Comment fields
    pub(crate) comment: &'a [u8],
    pub(crate) is_pseudogene: bool,

    // Map extension restriction enzyme
    pub(crate) r_enzyme: &'a [u8],

    // Number of gaps reported by CSeqMap_CI with fFindGap for delta sequences
    pub(crate) delta_gap_count: u32,
}

/// `CDbtag::IsSkippable`: case-insensitive match against BankIt, NCBIFILE, TMSMART.
///
/// NCBI reference (598d8ae6): c++/src/objects/general/Dbtag.cpp:75-80,281-285
pub(crate) fn is_skippable_db(db: &[u8]) -> bool {
    db.eq_ignore_ascii_case(b"BankIt")
        || db.eq_ignore_ascii_case(b"NCBIFILE")
        || db.eq_ignore_ascii_case(b"TMSMART")
}

/// Determines `m_UsePDBCompoundForDefline` from PDB/molecule/length/comment state.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1170-1192
pub(crate) fn compute_use_pdb_compound(
    is_pdb: bool,
    is_na: bool,
    length: u64,
    comment: &[u8],
) -> bool {
    if !is_pdb {
        return false;
    }
    if comment.is_empty() {
        return true;
    }
    if is_na {
        length < 25
            || contains_bytes(comment, b"COMPLETE GENOME")
            || contains_bytes(comment, b"CHROMOSOME XII")
            || contains_bytes(comment, b"Dna (5'")
            || contains_bytes(comment, b"SEQRES")
    } else {
        contains_bytes(comment, b"hypothetical protein")
            || contains_bytes(comment, b"uncharacterized protein")
            || contains_bytes(comment, b"putative uncharacterized protein")
            || contains_bytes(comment, b"putative protein")
            || contains_bytes(comment, b"SEQRES")
    }
}

/// Computes `m_UnverifiedPrefix` from the `Type` string fields of an `Unverified` user object.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:946-969
/// NCBI reference (598d8ae6): c++/src/objects/general/User_object.cpp:53-56,659-689,791-848
pub(crate) fn unverified_prefix_from_types(
    has_feature: bool,
    has_misassembled: bool,
    has_contaminant: bool,
    has_organism: bool,
) -> &'static [u8] {
    let mut prefix: &'static [u8] = b"UNVERIFIED: ";
    if has_feature {
        prefix = b"UNVERIFIED: ";
    }
    if has_misassembled {
        prefix = b"UNVERIFIED_ASMBLY: ";
    }
    if has_contaminant {
        prefix = b"UNVERIFIED_CONTAM: ";
    }
    if has_organism {
        prefix = b"UNVERIFIED_ORG: ";
    }
    prefix
}

/// Computes `m_UnreviewedPrefix` from whether an `Unreviewed` user object has an `Unannotated`
/// `Type` field.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:970-975
/// NCBI reference (598d8ae6): c++/src/objects/general/User_object.cpp:58,725-755,863-866
pub(crate) fn unreviewed_prefix_from_unannotated(has_unannotated: bool) -> &'static [u8] {
    if has_unannotated {
        b"UNREVIEWED_UNANNOT: "
    } else {
        b"UNREVIEWED: "
    }
}

/// Applies `MolInfo.tech` flags to `input`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:890-930
pub(crate) fn apply_molinfo_tech(input: &mut DeflineInput<'_>) {
    match input.mi_tech {
        mol_info::TECH_HTGS_0 | mol_info::TECH_HTGS_1 | mol_info::TECH_HTGS_2 => {
            input.htgs_unfinished = true;
            input.reconstruct = true;
            input.main_title = b"";
            input.htg_tech = true;
        }
        mol_info::TECH_HTGS_3 => {
            input.htg_tech = true;
        }
        mol_info::TECH_EST | mol_info::TECH_STS | mol_info::TECH_SURVEY => {
            input.is_est_sts_gss = true;
        }
        mol_info::TECH_WGS => {
            input.is_wgs = true;
        }
        mol_info::TECH_TSA => {
            input.is_tsa = true;
            if input.is_virtual {
                input.tsa_master = true;
            }
        }
        mol_info::TECH_TARGETED => {
            input.is_tls = true;
            if input.is_virtual {
                input.tls_master = true;
            }
        }
        _ => {}
    }
}

/// Applies GenBank/EMBL keyword clauses (`;` tokenized) to `input`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1122-1150
pub(crate) fn apply_keyword_clauses<'a>(
    input: &mut DeflineInput<'_>,
    clauses: impl IntoIterator<Item = &'a [u8]>,
) {
    for clause in clauses {
        for part in clause.split(|&b| b == b';').filter(|p| !p.is_empty()) {
            if part.eq_ignore_ascii_case(b"UNORDERED") {
                input.unordered = true;
            }
            if !input.htg_tech && !input.third_party {
                continue;
            }
            if part.eq_ignore_ascii_case(b"HTGS_DRAFT") {
                input.htgs_draft = true;
            } else if part.eq_ignore_ascii_case(b"HTGS_CANCELLED") {
                input.htgs_cancelled = true;
            } else if part.eq_ignore_ascii_case(b"HTGS_POOLED_MULTICLONE") {
                input.htgs_pooled = true;
            } else if part.eq_ignore_ascii_case(b"TPA:experimental") {
                input.tpa_exp = true;
            } else if part.eq_ignore_ascii_case(b"TPA:inferential") {
                input.tpa_inf = true;
            } else if part.eq_ignore_ascii_case(b"TPA:reassembly")
                || part.eq_ignore_ascii_case(b"TPA:assembly")
            {
                input.tpa_reasm = true;
            }
        }
    }
}
