//! All `CDeflineGenerator::x_SetTitleFrom*` builders for BioSource, NC, NM, NR, Patent, PDB,
//! Protein, SegSeq, WGS/TSA/TLS, and Map records.
//!
//! NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:64-120,316-426,1436-2796

use crate::ncbi_str::replace_in_place;

use super::flags::{genome, mol_info, topology, DeflineFeatures, DeflineInput};
use super::html::{
    clean_and_compress, contains_bytes, contains_nocase, find_bytes_from, find_nocase_rev,
    is_c_punct, is_c_space, trim_end_of, truncate_spaces_both,
};
use super::prefix_suffix::{append_protein_org_suffix, is_complete};
use super::UNNAMED_PROTEIN_PRODUCT;

/// `CDefLineJoiner` with modifier-prefix display disabled (`m_ShowMods` false).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:64-120
struct DeflineJoiner {
    buf: Vec<u8>,
}

impl DeflineJoiner {
    fn new() -> Self {
        Self { buf: Vec::new() }
    }

    fn add(&mut self, name: &[u8], value: &[u8], hide_type: bool) {
        if !hide_type && !name.is_empty() {
            self.buf.push(b' ');
            self.buf.extend_from_slice(name);
        }
        if !value.is_empty() {
            self.buf.push(b' ');
            self.buf.extend_from_slice(value);
        }
    }

    fn into_bytes(self) -> Vec<u8> {
        self.buf
    }
}

/// `s_IsVirusOrPhage`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:316-320
fn is_virus_or_phage(taxname: &[u8]) -> bool {
    contains_nocase(taxname, b"virus") || contains_nocase(taxname, b"phage")
}

/// `CDeflineGenerator::x_OrganelleName`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:323-426
pub(super) fn organelle_name(inp: &DeflineInput<'_>) -> &'static [u8] {
    let has_plasmid = !inp.plasmid.is_empty();
    match inp.genome {
        genome::CHLOROPLAST => b"chloroplast",
        genome::CHROMOPLAST => b"chromoplast",
        genome::KINETOPLAST => b"kinetoplast",
        genome::MITOCHONDRION => {
            if !inp.fasta_format && (has_plasmid || inp.is_wgs) {
                b"mitochondrial"
            } else {
                b"mitochondrion"
            }
        }
        genome::PLASTID => b"plastid",
        genome::MACRONUCLEAR => b"macronuclear",
        genome::EXTRACHROM => {
            if !inp.is_wgs {
                b"extrachromosomal"
            } else {
                b""
            }
        }
        genome::PLASMID => {
            if !inp.is_wgs {
                b"plasmid"
            } else {
                b""
            }
        }
        genome::CYANELLE => b"cyanelle",
        genome::PROVIRAL => {
            if !is_virus_or_phage(inp.taxname) {
                if has_plasmid || inp.is_wgs {
                    b"proviral"
                } else {
                    b"provirus"
                }
            } else {
                b""
            }
        }
        genome::VIRION => {
            if !is_virus_or_phage(inp.taxname) {
                b"virus"
            } else {
                b""
            }
        }
        genome::NUCLEOMORPH => {
            if !inp.is_wgs {
                b"nucleomorph"
            } else {
                b""
            }
        }
        genome::APICOPLAST => b"apicoplast",
        genome::LEUCOPLAST => b"leucoplast",
        genome::PROPLASTID => b"proplastid",
        genome::ENDOGENOUS_VIRUS => b"endogenous virus",
        genome::HYDROGENOSOME => b"hydrogenosome",
        genome::CHROMOSOME => b"chromosome",
        genome::CHROMATOPHORE => b"chromatophore",
        _ => b"",
    }
}

/// `CDeflineGenerator::x_DescribeClones`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1437-1468
fn append_clones(inp: &DeflineInput<'_>, has_clone: bool, joiner: &mut DeflineJoiner) {
    if inp.htgs_unfinished && inp.htgs_pooled && has_clone {
        joiner.add(b"clone", b", pooled multiple clones", true);
        return;
    }
    if inp.clone.is_empty() {
        return;
    }
    let count = 1 + inp.clone.iter().filter(|&&b| b == b';').count();
    if count > 3 {
        let num = count.to_string();
        joiner.add(b"clone", b", ", true);
        joiner.add(b"clone", num.as_bytes(), true);
        joiner.add(b"clone", b" clones", true);
    } else {
        joiner.add(b"clone", b" clone ", true);
        joiner.add(b"clone", inp.clone, true);
    }
}

/// `s_EndsWithStrain`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1470-1502
fn ends_with_strain(taxname: &[u8], strain: &[u8]) -> bool {
    if strain.len() >= taxname.len() {
        return false;
    }
    let Some(first_space) = taxname.iter().position(|&b| b == b' ') else {
        return false;
    };
    if !taxname[first_space + 1..].contains(&b' ') {
        return false;
    }
    let Some(pos) = find_nocase_rev(taxname, strain) else {
        return false;
    };
    if pos + strain.len() == taxname.len() {
        if pos > 0 {
            let ch = taxname[pos - 1];
            if is_c_punct(ch) || is_c_space(ch) {
                return true;
            }
        }
    } else if pos + strain.len() + 1 == taxname.len()
        && pos > 0
        && taxname[pos - 1] == b'\''
        && taxname[taxname.len() - 1] == b'\''
    {
        return true;
    }
    false
}

fn remove_colons_and_ws(bytes: &[u8]) -> Vec<u8> {
    bytes
        .iter()
        .copied()
        .filter(|&c| !matches!(c, b':' | b' ' | b'\t'))
        .collect()
}

fn remove_ws(bytes: &[u8]) -> Vec<u8> {
    bytes
        .iter()
        .copied()
        .filter(|&c| !matches!(c, b' ' | b'\t'))
        .collect()
}

fn before_semicolon(bytes: &[u8]) -> &[u8] {
    match bytes.iter().position(|&b| b == b';') {
        Some(pos) => &bytes[..pos],
        None => bytes,
    }
}

/// `s_AddVoucherAndIsolate`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1523-1542
fn add_voucher_and_isolate(inp: &DeflineInput<'_>, joiner: &mut DeflineJoiner) {
    if !inp.specimen_voucher.is_empty()
        && (inp.strain.is_empty()
            || remove_colons_and_ws(inp.strain) != remove_colons_and_ws(inp.specimen_voucher))
    {
        joiner.add(b"voucher", inp.specimen_voucher, false);
    }
    if !inp.isolate.is_empty()
        && inp.isolate != inp.specimen_voucher
        && !ends_with_strain(inp.taxname, inp.isolate)
        && remove_colons_and_ws(inp.specimen_voucher) != remove_ws(inp.isolate)
    {
        joiner.add(b"isolate", inp.isolate, false);
    }
}

/// `CDeflineGenerator::x_SetTitleFromBioSrc`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1545-1600
pub(super) fn set_title_from_biosrc(
    inp: &DeflineInput<'_>,
    organelle: &[u8],
    has_clone: bool,
) -> Vec<u8> {
    let mut joiner = DeflineJoiner::new();
    joiner.add(b"organism", inp.taxname, true);

    if !inp.strain.is_empty() {
        let add = before_semicolon(inp.strain);
        if !ends_with_strain(inp.taxname, add) {
            joiner.add(b"strain", add, false);
        }
    }
    if !inp.substrain.is_empty() {
        let add = before_semicolon(inp.substrain);
        if !ends_with_strain(inp.taxname, add) {
            joiner.add(b"substr.", add, false);
        }
    }
    if !inp.breed.is_empty() {
        joiner.add(b"breed", before_semicolon(inp.breed), false);
    }
    if !inp.cultivar.is_empty() {
        joiner.add(b"cultivar", before_semicolon(inp.cultivar), false);
    }

    add_voucher_and_isolate(inp, &mut joiner);

    if !inp.chromosome.is_empty() {
        joiner.add(b"location", b"chromosome", true);
        joiner.add(b"chromosome", inp.chromosome, true);
    } else if !inp.linkage_group.is_empty() {
        joiner.add(b"location", b"linkage group", true);
        joiner.add(b"linkage group", inp.linkage_group, true);
    } else if !inp.plasmid.is_empty() {
        joiner.add(b"location", organelle, true);
        joiner.add(b"plasmid name", inp.plasmid, true);
    } else if !organelle.is_empty() {
        joiner.add(b"location", organelle, true);
    }

    if has_clone {
        append_clones(inp, has_clone, &mut joiner);
    }
    if !inp.map.is_empty() {
        joiner.add(b"map", inp.map, false);
    }

    let mut out = joiner.into_bytes();
    truncate_spaces_both(&mut out);
    out
}

/// `CDeflineGenerator::x_SetTitleFromNC`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1603-1679
pub(super) fn set_title_from_nc(inp: &DeflineInput<'_>, organelle: &[u8]) -> Vec<u8> {
    if inp.mi_biomol != mol_info::BIOMOL_GENOMIC && inp.mi_biomol != mol_info::BIOMOL_OTHER_GENETIC
    {
        return Vec::new();
    }
    if inp.taxname.is_empty() {
        return Vec::new();
    }

    let mut joiner = DeflineJoiner::new();
    joiner.add(b"organism", inp.taxname, true);

    let mut add_gen_tag = false;
    if contains_nocase(inp.taxname, b"plasmid") {
        // Taxname already mentions plasmid; no location/genome tag.
    } else if inp.is_plasmid || !inp.plasmid.is_empty() {
        if inp.plasmid.is_empty() {
            joiner.add(b"", b"unnamed plasmid", true);
        } else {
            if !inp.is_plasmid {
                joiner.add(b"location", organelle, true);
            }
            if !contains_nocase(inp.plasmid, b"plasmid")
                && !contains_nocase(inp.plasmid, b"element")
            {
                joiner.add(b"plasmid", inp.plasmid, false);
            } else {
                joiner.add(b"", inp.plasmid, true);
            }
        }
    } else if !organelle.is_empty() {
        if inp.chromosome.is_empty() {
            if matches!(
                inp.genome,
                genome::MITOCHONDRION
                    | genome::CHLOROPLAST
                    | genome::KINETOPLAST
                    | genome::PLASTID
                    | genome::APICOPLAST
            ) {
                joiner.add(b"location", organelle, true);
            }
        } else {
            if !inp.is_chromosome {
                joiner.add(b"location", organelle, true);
            }
            joiner.add(b"chromosome", inp.chromosome, false);
        }
    } else if !inp.segment.is_empty() {
        if !contains_bytes(inp.segment, b"DNA")
            && !contains_bytes(inp.segment, b"RNA")
            && !contains_bytes(inp.segment, b"segment")
            && !contains_bytes(inp.segment, b"Segment")
        {
            joiner.add(b"segment", inp.segment, false);
        } else {
            joiner.add(b"", inp.segment, true);
        }
    } else if !inp.chromosome.is_empty() {
        joiner.add(b"chromosome", inp.chromosome, false);
    } else {
        add_gen_tag = true;
    }

    let complete = is_complete(inp.mi_completeness);
    let tag: &[u8] = match (add_gen_tag, complete) {
        (true, true) => b", complete genome",
        (true, false) => b", partial genome",
        (false, true) => b", complete sequence",
        (false, false) => b", partial sequence",
    };
    joiner.add(b"completeness", tag, true);

    let mut out = joiner.into_bytes();
    replace_in_place(&mut out, b"Plasmid", b"plasmid");
    replace_in_place(&mut out, b"Element", b"element");
    out
}

/// `x_FlyCG_PtoR`: converts `CG<digits>-P<alpha>` tokens to `CG<digits>-R<alpha>`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1682-1704
fn fly_cg_p_to_r(s: &mut [u8]) {
    let len = s.len();
    let mut pos = 0;
    while pos + 3 < len {
        let Some(found) = find_bytes_from(s, b"CG", pos) else {
            break;
        };
        pos = found;
        if pos > 0 && !is_c_space(s[pos - 1]) {
            pos += 2;
            continue;
        }
        pos += 2;
        while pos + 3 < len && s[pos].is_ascii_digit() {
            pos += 1;
        }
        if pos + 2 < len
            && s[pos] == b'-'
            && s[pos + 1] == b'P'
            && s[pos + 2].is_ascii_alphabetic()
            && (pos + 3 == len || matches!(s[pos + 3], b' ' | b',' | b';'))
        {
            s[pos + 1] = b'R';
        }
    }
}

/// `CDeflineGenerator::x_SetTitleFromNM`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1706-1755
pub(super) fn set_title_from_nm(
    inp: &DeflineInput<'_>,
    feats: Option<&dyn DeflineFeatures>,
) -> Vec<u8> {
    if inp.taxname.is_empty() {
        return Vec::new();
    }
    let Some((mut cds_label, gene_label)) = feats.and_then(|f| f.nm_single_gene_and_cds()) else {
        return Vec::new();
    };
    if inp.taxname.eq_ignore_ascii_case(b"Drosophila melanogaster") {
        fly_cg_p_to_r(&mut cds_label);
    }
    replace_in_place(&mut cds_label, b"isoform ", b"transcript variant ");
    let mut out = Vec::with_capacity(inp.taxname.len() + cds_label.len() + gene_label.len() + 12);
    out.extend_from_slice(inp.taxname);
    out.push(b' ');
    out.extend_from_slice(&cds_label);
    out.extend_from_slice(b" (");
    out.extend_from_slice(&gene_label);
    out.extend_from_slice(b"), mRNA");
    out
}

/// `CDeflineGenerator::x_SetTitleFromNR`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1758-1812
pub(super) fn set_title_from_nr(
    inp: &DeflineInput<'_>,
    feats: Option<&dyn DeflineFeatures>,
) -> Vec<u8> {
    if inp.taxname.is_empty() {
        return Vec::new();
    }
    let Some(gene_label) = feats.and_then(|f| f.nr_first_gene_label()) else {
        return Vec::new();
    };
    let mut out = Vec::with_capacity(inp.taxname.len() + gene_label.len() + 20);
    out.extend_from_slice(inp.taxname);
    out.push(b' ');
    out.extend_from_slice(&gene_label);
    out.extend_from_slice(b", ");
    let rna: &[u8] = match inp.mi_biomol {
        mol_info::BIOMOL_PRE_RNA => b"precursorRNA",
        mol_info::BIOMOL_MRNA => b"mRNA",
        mol_info::BIOMOL_RRNA => b"rRNA",
        mol_info::BIOMOL_TRNA => b"tRNA",
        mol_info::BIOMOL_SNRNA => b"snRNA",
        mol_info::BIOMOL_SCRNA => b"scRNA",
        mol_info::BIOMOL_CRNA => b"cRNA",
        mol_info::BIOMOL_SNORNA => b"snoRNA",
        mol_info::BIOMOL_TRANSCRIBED_RNA => b"miscRNA",
        mol_info::BIOMOL_NCRNA => b"ncRNA",
        mol_info::BIOMOL_TMRNA => b"tmRNA",
        _ => b"",
    };
    out.extend_from_slice(rna);
    out
}

/// `CDeflineGenerator::x_SetTitleFromPatent`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1815-1823
pub(super) fn set_title_from_patent(inp: &DeflineInput<'_>) -> Vec<u8> {
    let mut out = format!("Sequence {} from Patent ", inp.patent_sequence).into_bytes();
    out.extend_from_slice(inp.patent_country);
    out.push(b' ');
    out.extend_from_slice(inp.patent_number);
    out
}

/// `CDeflineGenerator::x_SetTitleFromPDB`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1825-1859
pub(super) fn set_title_from_pdb(inp: &DeflineInput<'_>) -> Vec<u8> {
    let chain_byte = [inp.pdb_chain as u8];
    let chain: Option<&[u8]> = if !inp.pdb_chain_id.is_empty() {
        Some(inp.pdb_chain_id)
    } else if (0x20..0x7f).contains(&chain_byte[0]) {
        Some(&chain_byte)
    } else {
        None
    };

    let Some(chain) = chain else {
        return inp.pdb_compound.to_vec();
    };

    let mut out =
        Vec::with_capacity(8 + chain.len() + inp.pdb_compound.len().max(inp.comment.len()));
    out.extend_from_slice(b"Chain ");
    out.extend_from_slice(chain);
    out.extend_from_slice(b", ");
    if inp.use_pdb_compound_for_defline {
        out.extend_from_slice(inp.pdb_compound);
    } else {
        let tail = match inp.comment.iter().position(|b| !b.is_ascii_digit()) {
            Some(pos) if inp.comment[pos] == b' ' => &inp.comment[pos..],
            _ => inp.comment,
        };
        out.extend_from_slice(tail);
    }
    out
}

/// `CDeflineGenerator::x_SetTitleFromGPipe`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1861-1913
pub(super) fn set_title_from_gpipe(
    inp: &DeflineInput<'_>,
    organelle: &[u8],
    has_clone: bool,
) -> Vec<u8> {
    let mut joiner = DeflineJoiner::new();
    joiner.add(b"organism", inp.taxname, true);
    if !organelle.is_empty() && contains_nocase(organelle, b"plasmid") {
        joiner.add(b"location", organelle, true);
    }
    if !inp.strain.is_empty() {
        let add = before_semicolon(inp.strain);
        if !ends_with_strain(inp.taxname, add) {
            joiner.add(b"strain", add, false);
        }
        let sub = before_semicolon(inp.substrain);
        if !ends_with_strain(inp.taxname, sub) {
            joiner.add(b"substr.", sub, false);
        }
    }
    if !inp.chromosome.is_empty() {
        joiner.add(b"chromosome", inp.chromosome, false);
    }
    if has_clone {
        append_clones(inp, has_clone, &mut joiner);
    }
    if !inp.map.is_empty() {
        joiner.add(b"map", inp.map, false);
    }
    if !inp.plasmid.is_empty() {
        if !contains_nocase(inp.plasmid, b"plasmid") && !contains_nocase(inp.plasmid, b"element") {
            joiner.add(b"plasmid", inp.plasmid, false);
        } else {
            joiner.add(b"", inp.plasmid, false);
        }
    }
    if is_complete(inp.mi_completeness) {
        joiner.add(b"completeness", b", complete sequence", true);
    }
    let mut out = joiner.into_bytes();
    truncate_spaces_both(&mut out);
    out
}

/// `CDeflineGenerator::x_CDShasLowQualityException` (`CLowQualityTextFsm`, case-insensitive).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2098-2121,4176-4183
pub(crate) fn cds_except_is_low_quality(except_text: &[u8]) -> bool {
    contains_nocase(except_text, b"heterogeneous population sequenced")
        || contains_nocase(except_text, b"low-quality sequence region")
        || contains_nocase(except_text, b"unextendable partial coding region")
}

/// `s_RemoveBracketedOrgFromEnd`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2149-2171
fn remove_bracketed_org_from_end(input: &[u8], taxname: &[u8]) -> Vec<u8> {
    if input.len() < 5 || taxname.is_empty() || input.last() != Some(&b']') {
        return input.to_vec();
    }
    let Some(cp) = input.iter().rposition(|&b| b == b'[') else {
        return input.to_vec();
    };
    let suffix = &input[cp + 1..];
    if suffix.starts_with(b"NAD") {
        return input.to_vec();
    }
    if suffix.len() == taxname.len() + 1 && suffix.starts_with(taxname) {
        let (cleaned, _) = clean_and_compress(&input[..cp], true);
        return cleaned;
    }
    input.to_vec()
}

/// `CDeflineGenerator::x_SetTitleFromProtein`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2386-2575
pub(super) fn set_title_from_protein(
    inp: &DeflineInput<'_>,
    feats: Option<&dyn DeflineFeatures>,
) -> Vec<u8> {
    let prot = feats.and_then(|f| f.longest_protein());
    let cds = feats.and_then(|f| f.cds_for_product());
    // NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2455-2464,2483-2502,2507-2517
    let gene = || feats.and_then(|f| f.gene_via_cds());
    let mut main_title = Vec::new();

    if let Some(prot) = &prot {
        let mut prefix: &[u8] = b"";
        for name in &prot.names {
            let trimmed = remove_bracketed_org_from_end(name, inp.taxname);
            main_title.extend_from_slice(prefix);
            main_title.extend_from_slice(&trimmed);
            if !inp.all_prot_names {
                break;
            }
            prefix = b"; ";
        }

        if !main_title.is_empty() {
            trim_end_of(&mut main_title, b".,;~ ");
            let mut offset = 0usize;
            let mut delta = 0usize;
            let mut comma: &[u8] = b"";
            let mut isoform: Vec<u8> = Vec::new();
            if main_title.starts_with(b"hypothetical protein") {
                offset = 20;
            } else if main_title.starts_with(b"uncharacterized protein") {
                offset = 23;
            }
            if offset > 0 && offset < main_title.len() {
                if main_title[offset] == b',' && main_title.get(offset + 1).copied() == Some(b' ') {
                    comma = b", isoform ";
                    delta = 2;
                }
                if main_title[offset] == b' ' {
                    comma = b" isoform ";
                    delta = 1;
                }
                if main_title
                    .get(offset + delta..)
                    .is_some_and(|rest| rest.starts_with(b"isoform "))
                {
                    isoform = main_title[offset + delta + 8..].to_vec();
                    main_title.truncate(offset);
                }
            }
            if main_title.eq_ignore_ascii_case(b"hypothetical protein")
                || main_title.eq_ignore_ascii_case(b"uncharacterized protein")
            {
                if let Some(lt) = gene().and_then(|g| g.locus_tag) {
                    if !lt.is_empty() {
                        main_title.push(b' ');
                        main_title.extend_from_slice(&lt);
                        main_title.extend_from_slice(comma);
                        main_title.extend_from_slice(&isoform);
                    }
                }
            }
        }

        if main_title.is_empty() && !inp.local_annots_only {
            if let Some(desc) = &prot.desc {
                main_title.extend_from_slice(desc);
            }
        }
        if main_title.is_empty() && !inp.local_annots_only {
            if let Some(act) = &prot.first_activity {
                main_title.extend_from_slice(act);
            }
        }
    }

    if main_title.is_empty() && !inp.local_annots_only {
        if let Some(gene) = gene() {
            if let Some(locus) = &gene.locus {
                main_title.extend_from_slice(locus);
            }
            if main_title.is_empty() {
                if let Some(syn) = &gene.first_syn {
                    main_title.extend_from_slice(syn);
                }
            }
            if main_title.is_empty() {
                if let Some(desc) = &gene.desc {
                    main_title.extend_from_slice(desc);
                }
            }
        }
        if !main_title.is_empty() {
            main_title.extend_from_slice(b" gene product");
        }
    }

    if main_title.is_empty() && !inp.local_annots_only {
        main_title.extend_from_slice(UNNAMED_PROTEIN_PRODUCT);
        if let Some(lt) = gene().and_then(|g| g.locus_tag) {
            if !lt.is_empty() {
                main_title.push(b' ');
                main_title.extend_from_slice(&lt);
            }
        }
    }

    if cds.as_ref().is_some_and(|c| c.has_low_quality_exception) {
        const LOW_QUAL: &[u8] = b"LOW QUALITY PROTEIN: ";
        if !contains_nocase(&main_title, LOW_QUAL) {
            main_title.splice(0..0, LOW_QUAL.iter().copied());
        }
    }

    trim_end_of(&mut main_title, b".,;~ ");

    if !is_complete(inp.mi_completeness) {
        main_title.extend_from_slice(b", partial");
    }

    append_protein_org_suffix(&mut main_title, inp, feats);
    main_title
}

/// `CDeflineGenerator::x_SetTitleFromSegSeq`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2646-2691
pub(super) fn set_title_from_seg_seq(
    inp: &DeflineInput<'_>,
    feats: Option<&dyn DeflineFeatures>,
) -> Vec<u8> {
    let mut joiner = DeflineJoiner::new();
    let taxname: &[u8] = if inp.taxname.is_empty() {
        b"Unknown"
    } else {
        inp.taxname
    };
    joiner.add(b"organism", taxname, true);

    let cds_info = if !inp.local_annots_only {
        feats.and_then(|f| f.segseq_info_via_cds())
    } else {
        None
    };

    let (product, locus, completeness): (&[u8], &[u8], &[u8]) = match &cds_info {
        Some(info) => (
            &info.product,
            &info.locus,
            if info.is_partial {
                b"partial"
            } else {
                b"complete"
            },
        ),
        None => {
            if !inp.strain.is_empty() && !ends_with_strain(taxname, inp.strain) {
                joiner.add(b"strain", inp.strain, false);
            } else if !inp.clone.is_empty() {
                append_clones(inp, true, &mut joiner);
            } else if !inp.isolate.is_empty() {
                joiner.add(b"isolate", inp.isolate, false);
            }
            (b"", b"", b"complete")
        }
    };

    if !product.is_empty() {
        joiner.add(b"product", product, true);
    }
    let mut out = joiner.into_bytes();
    if !locus.is_empty() {
        out.extend_from_slice(b" (");
        out.extend_from_slice(locus);
        out.push(b')');
    }
    if !product.is_empty() || !locus.is_empty() {
        out.extend_from_slice(b" gene, ");
        out.extend_from_slice(completeness);
        out.extend_from_slice(b" cds");
    }
    truncate_spaces_both(&mut out);
    out
}

/// `CDeflineGenerator::x_SetTitleFromWGS`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2694-2756
pub(super) fn set_title_from_wgs(inp: &DeflineInput<'_>, has_clone: bool) -> Vec<u8> {
    let mut joiner = DeflineJoiner::new();
    joiner.add(b"organism", inp.taxname, true);

    if !inp.strain.is_empty() {
        if !ends_with_strain(inp.taxname, inp.strain) {
            joiner.add(b"strain", before_semicolon(inp.strain), false);
        }
        if !inp.substrain.is_empty() && !ends_with_strain(inp.taxname, inp.substrain) {
            joiner.add(b"substr.", before_semicolon(inp.substrain), false);
        }
    } else if !inp.breed.is_empty() {
        joiner.add(b"breed", before_semicolon(inp.breed), false);
    } else if !inp.cultivar.is_empty() {
        joiner.add(b"cultivar", before_semicolon(inp.cultivar), false);
    }

    add_voucher_and_isolate(inp, &mut joiner);

    if !inp.chromosome.is_empty() {
        joiner.add(b"chromosome", inp.chromosome, false);
    } else if !inp.linkage_group.is_empty() {
        joiner.add(b"linkage group", inp.linkage_group, false);
    }
    if !inp.clone.is_empty() {
        append_clones(inp, has_clone, &mut joiner);
    }
    if !inp.map.is_empty() {
        joiner.add(b"map", inp.map, false);
    }
    if !inp.plasmid.is_empty() && inp.is_wgs {
        joiner.add(b"plasmid", inp.plasmid, false);
    }

    if inp.genome == genome::PLASMID && inp.topology == topology::CIRCULAR {
        // Omitted for circular plasmids.
    } else if inp.genome == genome::CHROMOSOME {
        // Omitted for chromosomes.
    } else if !inp.general_str.is_empty() {
        if inp.general_str != inp.chromosome && (!inp.is_wgs || inp.general_str != inp.plasmid) {
            joiner.add(b"", inp.general_str, true);
        }
    } else if inp.general_id > 0 {
        let tmp = inp.general_id.to_string();
        let tmp_b = tmp.as_bytes();
        if tmp_b != inp.chromosome && (!inp.is_wgs || tmp_b != inp.plasmid) {
            joiner.add(b"", tmp_b, true);
        }
    }

    let mut out = joiner.into_bytes();
    truncate_spaces_both(&mut out);
    out
}

/// `CDeflineGenerator::x_SetTitleFromMap`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2759-2796
pub(super) fn set_title_from_map(inp: &DeflineInput<'_>) -> Vec<u8> {
    let mut joiner = DeflineJoiner::new();
    joiner.add(b"organism", inp.taxname, true);

    if !inp.strain.is_empty() && !ends_with_strain(inp.taxname, inp.strain) {
        joiner.add(b"strain", before_semicolon(inp.strain), false);
    }
    if !inp.substrain.is_empty() && !ends_with_strain(inp.taxname, inp.substrain) {
        joiner.add(b"substr.", before_semicolon(inp.substrain), false);
    }
    if !inp.chromosome.is_empty() {
        joiner.add(b"chromosome", inp.chromosome, false);
    } else if inp.is_chromosome {
        joiner.add(b"location", b"chromosome", true);
    }
    if !inp.plasmid.is_empty() {
        joiner.add(b"plasmid", inp.plasmid, false);
    } else if inp.is_plasmid {
        joiner.add(b"location", b"plasmid", true);
    }
    if !inp.isolate.is_empty() {
        joiner.add(b"isolate", inp.isolate, false);
    }

    let mut out = joiner.into_bytes();
    if !inp.r_enzyme.is_empty() {
        out.extend_from_slice(b", ");
        out.extend_from_slice(inp.r_enzyme);
        out.extend_from_slice(b" whole genome map");
    }
    truncate_spaces_both(&mut out);
    out
}
