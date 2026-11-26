//! `CDeflineGenerator` prefix, suffix, completeness, and protein-title suffix adjustment.
//!
//! NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2123-2147,2798-3446

use super::flags::{genome, mol_info, DeflineFeatures, DeflineInput};
use super::html::{contains_bytes, find_nocase_rev, trim_end_of};

/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3431-3446
pub(super) const TPA_PREFIXES: [&str; 14] = [
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

/// `s_proteinOrganellePrefix` from `chloroplast` (2) to `chromatophore` (22), the non-empty ones.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2123-2147
pub(super) const PROTEIN_ORGANELLE_PREFIXES: [&[u8]; 15] = [
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

/// `s_proteinOrganellePrefix` indexed by `CBioSource::EGenome` (0..=22).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2123-2147
pub(super) fn protein_organelle_prefix(g: i32) -> &'static [u8] {
    match g {
        genome::CHLOROPLAST => b"chloroplast",
        genome::CHROMOPLAST => b"chromoplast",
        genome::KINETOPLAST => b"kinetoplast",
        genome::MITOCHONDRION => b"mitochondrion",
        genome::PLASTID => b"plastid",
        genome::MACRONUCLEAR => b"macronuclear",
        genome::PLASMID => b"plasmid",
        genome::CYANELLE => b"cyanelle",
        genome::NUCLEOMORPH => b"nucleomorph",
        genome::APICOPLAST => b"apicoplast",
        genome::LEUCOPLAST => b"leucoplast",
        genome::PROPLASTID => b"protoplast",
        genome::ENDOGENOUS_VIRUS => b"endogenous virus",
        genome::HYDROGENOSOME => b"hydrogenosome",
        genome::CHROMATOPHORE => b"chromatophore",
        _ => b"",
    }
}

/// `CDeflineGenerator::x_IsComplete`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3416-3428
pub(super) fn is_complete(completeness: i32) -> bool {
    !matches!(
        completeness,
        mol_info::COMPLETENESS_PARTIAL
            | mol_info::COMPLETENESS_NO_LEFT
            | mol_info::COMPLETENESS_NO_RIGHT
            | mol_info::COMPLETENESS_NO_ENDS
    )
}

/// `CDeflineGenerator::x_SetPrefix` (default constructor path where `m_Idx` is null).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2799-2866
pub(super) fn set_prefix(
    inp: &DeflineInput<'_>,
    main_title: &[u8],
    multispecies: bool,
) -> &'static [u8] {
    if inp.is_unverified {
        if !contains_bytes(main_title, b"UNVERIFIED") {
            return inp.unverified_prefix;
        }
        return b"";
    }
    if inp.is_unreviewed {
        if !contains_bytes(main_title, b"UNREVIEWED") {
            return inp.unreviewed_prefix;
        }
        return b"";
    }
    if inp.third_party {
        if inp.tpa_exp {
            return b"TPA_exp: ";
        }
        if inp.tpa_inf {
            return b"TPA_inf: ";
        }
        if inp.tpa_reasm {
            return b"TPA_asm: ";
        }
        return b"TPA: ";
    }
    if inp.is_tsa {
        return b"TSA: ";
    }
    if inp.is_tls {
        return b"TLS: ";
    }
    if multispecies && inp.is_wp {
        return b"MULTISPECIES: ";
    }
    if inp.is_pseudogene && !contains_bytes(main_title, b"PUTATIVE PSEUDOGENE") {
        return b"PUTATIVE PSEUDOGENE: ";
    }
    b""
}

/// `CDeflineGenerator::x_SetSuffix`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2869-3047
pub(super) fn set_suffix(
    inp: &DeflineInput<'_>,
    main_title: &[u8],
    organelle: &[u8],
    append_complete: bool,
) -> Vec<u8> {
    let mut kind = Vec::new();
    let mut study: &[u8] = b"";
    let mut comp: &[u8] = b"";

    match inp.mi_tech {
        mol_info::TECH_HTGS_0 => {
            if !contains_bytes(main_title, b"LOW-PASS") {
                kind.extend_from_slice(b", LOW-PASS SEQUENCE SAMPLING");
            }
        }
        mol_info::TECH_HTGS_1 | mol_info::TECH_HTGS_2 => {
            if inp.htgs_draft {
                if !contains_bytes(main_title, b"WORKING DRAFT") {
                    kind.extend_from_slice(b", WORKING DRAFT SEQUENCE");
                }
            } else if !inp.htgs_cancelled && !contains_bytes(main_title, b"SEQUENCING IN") {
                kind.extend_from_slice(b", *** SEQUENCING IN PROGRESS ***");
            }
            let un = if inp.mi_tech == mol_info::TECH_HTGS_1 {
                "un"
            } else {
                ""
            };
            if inp.is_delta {
                let pieces = inp.delta_gap_count.saturating_add(1);
                if pieces != 1 {
                    kind.extend_from_slice(format!(", {pieces} {un}ordered pieces").as_bytes());
                }
            }
        }
        mol_info::TECH_HTGS_3 => {
            if !contains_bytes(main_title, b"complete sequence") {
                kind.extend_from_slice(b", complete sequence");
            }
        }
        mol_info::TECH_EST => {
            if !contains_bytes(main_title, b"mRNA sequence") {
                kind.extend_from_slice(b", mRNA sequence");
            }
        }
        mol_info::TECH_STS => {
            if !contains_bytes(main_title, b"sequence tagged site") {
                kind.extend_from_slice(b", sequence tagged site");
            }
        }
        mol_info::TECH_SURVEY => {
            if !contains_bytes(main_title, b"genomic survey sequence") {
                kind.extend_from_slice(b", genomic survey sequence");
            }
        }
        mol_info::TECH_WGS => {
            if inp.wgs_master {
                if !contains_bytes(main_title, b"whole genome shotgun sequencing project") {
                    kind.extend_from_slice(b", whole genome shotgun sequencing project");
                }
            } else if !contains_bytes(main_title, b"whole genome shotgun sequence") {
                if !organelle.is_empty() && !contains_bytes(main_title, organelle) {
                    let is_mito = organelle.eq_ignore_ascii_case(b"mitochondrial")
                        || organelle.eq_ignore_ascii_case(b"mitochondrion");
                    let has_mito = contains_bytes(main_title, b"mitochondrial")
                        || contains_bytes(main_title, b"mitochondrion");
                    let is_chr = organelle.eq_ignore_ascii_case(b"chromosome");
                    let has_chr = contains_bytes(main_title, b"linkage group")
                        || contains_bytes(main_title, b"chromosome");
                    if !(is_mito && has_mito || is_chr && has_chr) {
                        kind.push(b' ');
                        kind.extend_from_slice(organelle);
                    }
                }
                kind.extend_from_slice(b", whole genome shotgun sequence");
            }
        }
        mol_info::TECH_TSA => {
            if inp.tsa_master {
                if !contains_bytes(main_title, b"transcriptome shotgun assembly") {
                    kind.extend_from_slice(b", transcriptome shotgun assembly");
                }
            } else if !contains_bytes(main_title, b"RNA sequence") {
                match inp.mi_biomol {
                    mol_info::BIOMOL_MRNA => kind.extend_from_slice(b", mRNA sequence"),
                    mol_info::BIOMOL_RRNA => kind.extend_from_slice(b", rRNA sequence"),
                    mol_info::BIOMOL_NCRNA => kind.extend_from_slice(b", ncRNA sequence"),
                    mol_info::BIOMOL_PRE_RNA
                    | mol_info::BIOMOL_SNRNA
                    | mol_info::BIOMOL_SCRNA
                    | mol_info::BIOMOL_CRNA
                    | mol_info::BIOMOL_SNORNA
                    | mol_info::BIOMOL_TRANSCRIBED_RNA => {
                        kind.extend_from_slice(b", transcribed RNA sequence");
                    }
                    _ => {}
                }
            }
        }
        mol_info::TECH_TARGETED => {
            if inp.tls_master {
                if !contains_bytes(main_title, b"targeted locus study") {
                    kind.extend_from_slice(b", targeted locus study");
                }
            } else if !contains_bytes(main_title, b"sequence") {
                kind.extend_from_slice(b", sequence");
            }
            if !inp.targeted_locus.is_empty() && !contains_bytes(main_title, inp.targeted_locus) {
                study = inp.targeted_locus;
            }
        }
        _ => {}
    }

    if append_complete
        && !contains_bytes(main_title, b"complete")
        && !contains_bytes(main_title, b"partial")
        && inp.mi_completeness == mol_info::COMPLETENESS_COMPLETE
    {
        if inp.is_plasmid {
            comp = b", complete sequence";
        } else if matches!(
            inp.genome,
            genome::MITOCHONDRION
                | genome::CHLOROPLAST
                | genome::KINETOPLAST
                | genome::PLASTID
                | genome::APICOPLAST
        ) {
            comp = b", complete genome";
        } else if inp.is_chromosome {
            if !inp.chromosome.is_empty() || !inp.linkage_group.is_empty() {
                comp = b", complete sequence";
            } else {
                comp = b", complete genome";
            }
        }
    }

    if inp.unordered && inp.is_delta && inp.delta_gap_count > 0 {
        let pieces = inp.delta_gap_count.saturating_add(1);
        kind.extend_from_slice(format!(", {pieces} unordered pieces").as_bytes());
    }

    let mut suffix = Vec::with_capacity(study.len() + kind.len() + comp.len() + 2);
    if !study.is_empty() {
        suffix.push(b' ');
        suffix.extend_from_slice(study);
        suffix.push(b' ');
    }
    suffix.extend_from_slice(&kind);
    suffix.extend_from_slice(comp);
    suffix
}

/// `s_TitleEndsInOrganism`: the position of `[` where `title` ends in ` [<taxname>]` (the name in
/// any ASCII case).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3090-3110
/// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2933-2977
pub(super) fn title_ends_in_organism(title: &[u8], taxname: &[u8]) -> Option<usize> {
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

/// Whether `taxname` allows replacing the protein's organism suffix from an overlapping `Biosrc`
/// feature via CDS (`x_GetSourceFeatViaCDS`).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2558-2563,3383-3388
pub(super) fn may_use_cds_source_taxname(taxname: &[u8], local_annots_only: bool) -> bool {
    !local_annots_only
        && (taxname.is_empty()
            || (!taxname.eq_ignore_ascii_case(b"synthetic construct")
                && !taxname.eq_ignore_ascii_case(b"artificial sequence")
                && !contains_bytes(taxname, b"vector")
                && !contains_bytes(taxname, b"Vector")))
}

/// Appends the protein `(organelle)` and `[taxname]` / `[sk1][sk2]` suffix shared by
/// `x_SetTitleFromProtein` and `x_AdjustProteinTitleSuffix`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2543-2574,3368-3413
pub(super) fn append_protein_org_suffix(
    title: &mut Vec<u8>,
    inp: &DeflineInput<'_>,
    feats: Option<&dyn DeflineFeatures>,
) {
    if inp.omit_taxonomic_name {
        return;
    }
    let organelle = protein_organelle_prefix(inp.genome);
    if !organelle.is_empty() && !inp.taxname.is_empty() {
        title.extend_from_slice(b" (");
        title.extend_from_slice(organelle);
        title.push(b')');
    }

    let cds_taxname = if may_use_cds_source_taxname(inp.taxname, inp.local_annots_only) {
        feats.and_then(|f| f.source_taxname_via_cds())
    } else {
        None
    };
    let taxname = cds_taxname.as_deref().unwrap_or(inp.taxname);

    if inp.is_cross_kingdom
        && !inp.first_super_kingdom.is_empty()
        && !inp.second_super_kingdom.is_empty()
    {
        title.extend_from_slice(b" [");
        title.extend_from_slice(inp.first_super_kingdom);
        title.extend_from_slice(b"][");
        title.extend_from_slice(inp.second_super_kingdom);
        title.push(b']');
    } else if !taxname.is_empty() {
        title.extend_from_slice(b" [");
        title.extend_from_slice(taxname);
        title.push(b']');
    }
}

/// Full `CDeflineGenerator::x_AdjustProteinTitleSuffix` (used by both database/FASTA callers and
/// full ASN.1 Bioseqs).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3263-3414
pub(super) fn adjust_protein_title_suffix_full(
    title: &mut Vec<u8>,
    inp: &DeflineInput<'_>,
    feats: Option<&dyn DeflineFeatures>,
) {
    if !inp.has_source {
        return;
    }
    trim_end_of(title, b".,;~ ");
    if title.len() <= inp.taxname.len() + 4 {
        return;
    }

    let mut tpos = title_ends_in_organism(title, inp.taxname);
    if tpos.is_none() {
        let mut binomial = Vec::with_capacity(inp.genus.len() + 1 + inp.species.len());
        binomial.extend_from_slice(inp.genus);
        binomial.push(b' ');
        binomial.extend_from_slice(inp.species);
        tpos = title_ends_in_organism(title, &binomial);
        if tpos.is_none() && inp.is_cross_kingdom {
            if let Some(pos) = find_nocase_rev(title, b"][") {
                title.truncate(pos + 1);
                trim_end_of(title, b".,;~ ");
                tpos = title_ends_in_organism(title, inp.taxname);
            }
        }
    }
    let Some(tpos) = tpos else {
        return;
    };

    title.truncate(tpos);
    trim_end_of(title, b".,;~ ");

    let mut opos = None;
    if title.len() > 2 && title[title.len() - 1] == b')' {
        if let Some(pos) = title.iter().rposition(|&c| c == b'(') {
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
    if !is_complete(inp.mi_completeness) {
        title.extend_from_slice(b", partial");
    }
    append_protein_org_suffix(title, inp, feats);
}

/// `x_AdjustProteinTitleSuffix` for a database protein Bioseq with only `taxname`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3263-3414
#[cfg(test)]
pub(super) fn adjust_protein_title_suffix(title: &mut Vec<u8>, taxname: &[u8]) {
    let inp = DeflineInput {
        has_source: true,
        taxname,
        ..DeflineInput::default()
    };
    adjust_protein_title_suffix_full(title, &inp, None);
}
