//! Typed ASN.1 Bioseq descriptor/extension extraction and compact TSE feature index for
//! `CDeflineGenerator` in `losat makeblastdb -input_type asn1_txt|asn1_bin`.
//!
//! NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:612-1434,1916-2096,2578-2645,3928-4107
//! NCBI reference (598d8ae6): c++/src/objmgr/util/feature.cpp:103-248,586-776,2731-3340
//! NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:1060-1275,2548-2589
//! NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:88-115,481-513,1653-1850
//! NCBI reference (598d8ae6): c++/src/objmgr/seq_map.cpp:219-292,1359-1435

use std::cell::RefCell;
use std::collections::HashMap;
use std::sync::OnceLock;

use super::fasta_source::{post_seq_id_warning, Warn};
use super::WriteDbError;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::defline::flags::{
    apply_keyword_clauses, apply_molinfo_tech, compute_use_pdb_compound, genome, orgmod, seq_repr,
    subsource, topology, unreviewed_prefix_from_unannotated, unverified_prefix_from_types,
    CdsForProductInfo, DeflineFeatures, DeflineInput, ProteinFeatInfo, SegSeqCdsInfo,
};
use crate::report::defline::titles::cds_except_is_low_quality;
use crate::report::defline::{classify_seq_ids_full, generate_defline_core};
use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::{member_index, members, Value};
use crate::seqid::label;
use crate::seqid::mask_handle::{mask_handle_key, GiimTree, TreeHandle};
use crate::seqid::model::ParsedSeqId;

/// Cached ASN.1 table member indices for `CDeflineGenerator` descriptor/feature extraction.
struct DeflineSchema {
    seqdesc_title: usize,
    seqdesc_source: usize,
    seqdesc_molinfo: usize,
    seqdesc_user: usize,
    seqdesc_comment: usize,
    seqdesc_genbank: usize,
    seqdesc_embl: usize,
    seqdesc_pdb: usize,
    molinfo_biomol: usize,
    molinfo_tech: usize,
    molinfo_completeness: usize,
    user_obj_type: usize,
    user_obj_data: usize,
    user_field_label: usize,
    user_field_data: usize,
    user_field_str: usize,
    object_id_id: usize,
    object_id_str: usize,
    dbtag_db: usize,
    dbtag_tag: usize,
    biosrc_genome: usize,
    biosrc_org: usize,
    biosrc_subtype: usize,
    subsrc_subtype: usize,
    subsrc_name: usize,
    org_ref_taxname: usize,
    org_ref_orgname: usize,
    orgname_name: usize,
    orgname_binomial: usize,
    orgname_partial: usize,
    orgname_mod: usize,
    binomial_genus: usize,
    binomial_species: usize,
    taxelem_fixed_level: usize,
    taxelem_level: usize,
    taxelem_name: usize,
    orgmod_subtype: usize,
    orgmod_subname: usize,
    gb_block_keywords: usize,
    embl_block_keywords: usize,
    pdb_block_compound: usize,
    seq_inst_topology: usize,
    seq_inst_seq_data: usize,
    seq_inst_ext: usize,
    seq_data_gap: usize,
    seq_ext_seg: usize,
    seq_ext_map: usize,
    seq_ext_delta: usize,
    delta_loc: usize,
    delta_literal: usize,
    seq_lit_data: usize,
    seq_annot_data: usize,
    seq_annot_ftable: usize,
    seq_feat_id: usize,
    seq_feat_data: usize,
    seq_feat_partial: usize,
    seq_feat_except: usize,
    seq_feat_comment: usize,
    seq_feat_product: usize,
    seq_feat_location: usize,
    seq_feat_qual: usize,
    seq_feat_xref: usize,
    seq_feat_except_text: usize,
    feat_id_local: usize,
    gb_qual_qual: usize,
    gb_qual_val: usize,
    xref_id: usize,
    xref_data: usize,
    feat_gene: usize,
    feat_cdregion: usize,
    feat_prot: usize,
    feat_rsite: usize,
    feat_biosrc: usize,
    rsite_str: usize,
    prot_name: usize,
    prot_desc: usize,
    prot_activity: usize,
    prot_db: usize,
    prot_processed: usize,
    gene_locus: usize,
    gene_desc: usize,
    gene_maploc: usize,
    gene_db: usize,
    gene_syn: usize,
    gene_locus_tag: usize,
    cdregion_orf: usize,
    cdregion_frame: usize,
    seqloc_null: usize,
    seqloc_empty: usize,
    seqloc_whole: usize,
    seqloc_int: usize,
    seqloc_packed_int: usize,
    seqloc_pnt: usize,
    seqloc_packed_pnt: usize,
    seqloc_mix: usize,
    seqloc_equiv: usize,
    seqloc_bond: usize,
    seqint_from: usize,
    seqint_to: usize,
    seqint_strand: usize,
    seqint_id: usize,
    seqpnt_point: usize,
    seqpnt_strand: usize,
    seqpnt_id: usize,
    packed_pnt_strand: usize,
    packed_pnt_id: usize,
    packed_pnt_points: usize,
    seqbond_a: usize,
    seqbond_b: usize,
}

fn dschema() -> &'static DeflineSchema {
    static CACHE: OnceLock<DeflineSchema> = OnceLock::new();
    CACHE.get_or_init(|| {
        let user_field_data = member_index(named::USER_FIELD, "data");
        let udata_ty = members(named::USER_FIELD)[user_field_data].ty;
        let orgname_name = member_index(named::ORGNAME, "name");
        let on_name_ty = members(named::ORGNAME)[orgname_name].ty;
        let seq_annot_data = member_index(named::SEQ_ANNOT, "data");
        let annot_data_ty = members(named::SEQ_ANNOT)[seq_annot_data].ty;
        DeflineSchema {
            seqdesc_title: member_index(named::SEQDESC, "title"),
            seqdesc_source: member_index(named::SEQDESC, "source"),
            seqdesc_molinfo: member_index(named::SEQDESC, "molinfo"),
            seqdesc_user: member_index(named::SEQDESC, "user"),
            seqdesc_comment: member_index(named::SEQDESC, "comment"),
            seqdesc_genbank: member_index(named::SEQDESC, "genbank"),
            seqdesc_embl: member_index(named::SEQDESC, "embl"),
            seqdesc_pdb: member_index(named::SEQDESC, "pdb"),
            molinfo_biomol: member_index(named::MOLINFO, "biomol"),
            molinfo_tech: member_index(named::MOLINFO, "tech"),
            molinfo_completeness: member_index(named::MOLINFO, "completeness"),
            user_obj_type: member_index(named::USER_OBJECT, "type"),
            user_obj_data: member_index(named::USER_OBJECT, "data"),
            user_field_label: member_index(named::USER_FIELD, "label"),
            user_field_data,
            user_field_str: member_index(udata_ty, "str"),
            object_id_id: member_index(named::OBJECT_ID, "id"),
            object_id_str: member_index(named::OBJECT_ID, "str"),
            dbtag_db: member_index(named::DBTAG, "db"),
            dbtag_tag: member_index(named::DBTAG, "tag"),
            biosrc_genome: member_index(named::BIOSOURCE, "genome"),
            biosrc_org: member_index(named::BIOSOURCE, "org"),
            biosrc_subtype: member_index(named::BIOSOURCE, "subtype"),
            subsrc_subtype: member_index(named::SUBSOURCE, "subtype"),
            subsrc_name: member_index(named::SUBSOURCE, "name"),
            org_ref_taxname: member_index(named::ORG_REF, "taxname"),
            org_ref_orgname: member_index(named::ORG_REF, "orgname"),
            orgname_name,
            orgname_binomial: member_index(on_name_ty, "binomial"),
            orgname_partial: member_index(on_name_ty, "partial"),
            orgname_mod: member_index(named::ORGNAME, "mod"),
            binomial_genus: member_index(named::BINOMIALORGNAME, "genus"),
            binomial_species: member_index(named::BINOMIALORGNAME, "species"),
            taxelem_fixed_level: member_index(named::TAXELEMENT, "fixed-level"),
            taxelem_level: member_index(named::TAXELEMENT, "level"),
            taxelem_name: member_index(named::TAXELEMENT, "name"),
            orgmod_subtype: member_index(named::ORGMOD, "subtype"),
            orgmod_subname: member_index(named::ORGMOD, "subname"),
            gb_block_keywords: member_index(named::GB_BLOCK, "keywords"),
            embl_block_keywords: member_index(named::EMBL_BLOCK, "keywords"),
            pdb_block_compound: member_index(named::PDB_BLOCK, "compound"),
            seq_inst_topology: member_index(named::SEQ_INST, "topology"),
            seq_inst_seq_data: member_index(named::SEQ_INST, "seq-data"),
            seq_inst_ext: member_index(named::SEQ_INST, "ext"),
            seq_data_gap: member_index(named::SEQ_DATA, "gap"),
            seq_ext_seg: member_index(named::SEQ_EXT, "seg"),
            seq_ext_map: member_index(named::SEQ_EXT, "map"),
            seq_ext_delta: member_index(named::SEQ_EXT, "delta"),
            delta_loc: member_index(named::DELTA_SEQ, "loc"),
            delta_literal: member_index(named::DELTA_SEQ, "literal"),
            seq_lit_data: member_index(named::SEQ_LITERAL, "seq-data"),
            seq_annot_data,
            seq_annot_ftable: member_index(annot_data_ty, "ftable"),
            seq_feat_id: member_index(named::SEQ_FEAT, "id"),
            seq_feat_data: member_index(named::SEQ_FEAT, "data"),
            seq_feat_partial: member_index(named::SEQ_FEAT, "partial"),
            seq_feat_except: member_index(named::SEQ_FEAT, "except"),
            seq_feat_comment: member_index(named::SEQ_FEAT, "comment"),
            seq_feat_product: member_index(named::SEQ_FEAT, "product"),
            seq_feat_location: member_index(named::SEQ_FEAT, "location"),
            seq_feat_qual: member_index(named::SEQ_FEAT, "qual"),
            seq_feat_xref: member_index(named::SEQ_FEAT, "xref"),
            seq_feat_except_text: member_index(named::SEQ_FEAT, "except-text"),
            feat_id_local: member_index(named::FEAT_ID, "local"),
            gb_qual_qual: member_index(named::GB_QUAL, "qual"),
            gb_qual_val: member_index(named::GB_QUAL, "val"),
            xref_id: member_index(named::SEQFEATXREF, "id"),
            xref_data: member_index(named::SEQFEATXREF, "data"),
            feat_gene: member_index(named::SEQFEATDATA, "gene"),
            feat_cdregion: member_index(named::SEQFEATDATA, "cdregion"),
            feat_prot: member_index(named::SEQFEATDATA, "prot"),
            feat_rsite: member_index(named::SEQFEATDATA, "rsite"),
            feat_biosrc: member_index(named::SEQFEATDATA, "biosrc"),
            rsite_str: member_index(named::RSITE_REF, "str"),
            prot_name: member_index(named::PROT_REF, "name"),
            prot_desc: member_index(named::PROT_REF, "desc"),
            prot_activity: member_index(named::PROT_REF, "activity"),
            prot_db: member_index(named::PROT_REF, "db"),
            prot_processed: member_index(named::PROT_REF, "processed"),
            gene_locus: member_index(named::GENE_REF, "locus"),
            gene_desc: member_index(named::GENE_REF, "desc"),
            gene_maploc: member_index(named::GENE_REF, "maploc"),
            gene_db: member_index(named::GENE_REF, "db"),
            gene_syn: member_index(named::GENE_REF, "syn"),
            gene_locus_tag: member_index(named::GENE_REF, "locus-tag"),
            cdregion_orf: member_index(named::CDREGION, "orf"),
            cdregion_frame: member_index(named::CDREGION, "frame"),
            seqloc_null: member_index(named::SEQ_LOC, "null"),
            seqloc_empty: member_index(named::SEQ_LOC, "empty"),
            seqloc_whole: member_index(named::SEQ_LOC, "whole"),
            seqloc_int: member_index(named::SEQ_LOC, "int"),
            seqloc_packed_int: member_index(named::SEQ_LOC, "packed-int"),
            seqloc_pnt: member_index(named::SEQ_LOC, "pnt"),
            seqloc_packed_pnt: member_index(named::SEQ_LOC, "packed-pnt"),
            seqloc_mix: member_index(named::SEQ_LOC, "mix"),
            seqloc_equiv: member_index(named::SEQ_LOC, "equiv"),
            seqloc_bond: member_index(named::SEQ_LOC, "bond"),
            seqint_from: member_index(named::SEQ_INTERVAL, "from"),
            seqint_to: member_index(named::SEQ_INTERVAL, "to"),
            seqint_strand: member_index(named::SEQ_INTERVAL, "strand"),
            seqint_id: member_index(named::SEQ_INTERVAL, "id"),
            seqpnt_point: member_index(named::SEQ_POINT, "point"),
            seqpnt_strand: member_index(named::SEQ_POINT, "strand"),
            seqpnt_id: member_index(named::SEQ_POINT, "id"),
            packed_pnt_strand: member_index(named::PACKED_SEQPNT, "strand"),
            packed_pnt_id: member_index(named::PACKED_SEQPNT, "id"),
            packed_pnt_points: member_index(named::PACKED_SEQPNT, "points"),
            seqbond_a: member_index(named::SEQ_BOND, "a"),
            seqbond_b: member_index(named::SEQ_BOND, "b"),
        }
    })
}

fn as_choice(val: &Value) -> Option<(usize, &Value)> {
    match val {
        Value::Choice(c, inner) => Some((usize::from(*c), inner.as_ref())),
        _ => None,
    }
}

fn as_list(val: &Value) -> Option<&[Value]> {
    match val {
        Value::List(items) => Some(items.as_slice()),
        _ => None,
    }
}

fn obj_id_str(val: &Value) -> Option<&[u8]> {
    let s = dschema();
    match as_choice(val) {
        Some((c, Value::Str(b))) if c == s.object_id_str => Some(b),
        _ => None,
    }
}

fn obj_id_int(val: &Value) -> Option<i32> {
    let s = dschema();
    match as_choice(val) {
        Some((c, Value::Int(n))) if c == s.object_id_id => Some(*n as i32),
        _ => None,
    }
}

/// `CSeq_id_Handle` identity key shared between duplicate checks and TSE feature resolution.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum HandleKey {
    Bytes(Vec<u8>),
    Giim(i32, usize),
}

/// Interner mapping `ParsedSeqId` to `HandleKey` (including `Giim` handles) and formatted FASTA.
#[derive(Default)]
pub(crate) struct SeqIdInterner {
    giim: GiimTree,
}

impl SeqIdInterner {
    pub(crate) fn key_and_fasta(&mut self, id: &ParsedSeqId) -> (Option<HandleKey>, Vec<u8>) {
        let hkey = self.key_only(id);
        let mut fasta = Vec::new();
        label::fasta(id, &mut fasta);
        (hkey, fasta)
    }

    pub(crate) fn key_only(&mut self, id: &ParsedSeqId) -> Option<HandleKey> {
        match mask_handle_key(id) {
            Ok(Some(key)) => Some(HandleKey::Bytes(key)),
            Ok(None) => None,
            Err(TreeHandle::Giim { id: num }) => {
                let ParsedSeqId::Giim(g) = id else {
                    unreachable!("TreeHandle::Giim only comes from ParsedSeqId::Giim");
                };
                let (h, _) = self.giim.find_or_create(g);
                Some(HandleKey::Giim(num, h))
            }
        }
    }

    pub(crate) fn giim_find_or_create(
        &mut self,
        g: &crate::seqid::model::GiimportId,
        num: i32,
    ) -> (usize, ParsedSeqId) {
        let (h, _) = self.giim.find_or_create(g);
        let first = ParsedSeqId::Giim(self.giim.made_from(num, h).clone());
        (h, first)
    }
}

/// One interval/whole/point piece of a `Seq-loc`.
#[derive(Clone, Debug)]
struct LocPiece {
    id_key: Option<HandleKey>,
    from: u64,
    /// Exclusive upper bound (`to + 1`) when `is_whole` is false.
    to_open: u64,
    is_whole: bool,
    /// `ENa_strand` (0 unknown, 1 plus, 2 minus, 3 both, 4 both-rev, 255 other).
    strand: i32,
}

/// Compact representation of a `Seq-loc` extracted during the TSE walk.
#[derive(Clone, Debug, Default)]
pub(crate) struct CompactLoc {
    pieces: Vec<LocPiece>,
    /// `sequence::GetStrand` value for `s_GetCdregionLabel` ORF formatting and overlap checks.
    strand: i32,
    /// Whether the top-level `Seq-loc` is `mix` or `packed-int` (`CAnnotObjectType_Less`).
    is_complex: bool,
    /// FASTA string of the first Seq-id (for `cannot find sequence` error on missing product).
    first_id_fasta: Option<Vec<u8>>,
}

#[derive(Clone, Debug, Default)]
struct CompactGeneRef {
    locus: Option<Vec<u8>>,
    desc: Option<Vec<u8>>,
    locus_tag: Option<Vec<u8>>,
    first_syn: Option<Vec<u8>>,
    first_db_label: Option<Vec<u8>>,
    maploc: Option<Vec<u8>>,
}

impl CompactGeneRef {
    /// `CGene_ref::GetLabel`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Gene_ref.cpp:57-72
    fn get_label(&self) -> Vec<u8> {
        if let Some(v) = &self.locus {
            return v.clone();
        }
        if let Some(v) = &self.desc {
            return v.clone();
        }
        if let Some(v) = &self.locus_tag {
            return v.clone();
        }
        if let Some(v) = &self.first_syn {
            return v.clone();
        }
        if let Some(v) = &self.first_db_label {
            return v.clone();
        }
        if let Some(v) = &self.maploc {
            return v.clone();
        }
        Vec::new()
    }
}

#[derive(Clone, Debug)]
struct ProtFeature {
    seq_order: usize,
    loc: CompactLoc,
    names: Vec<Vec<u8>>,
    desc: Option<Vec<u8>>,
    first_activity: Option<Vec<u8>>,
    first_db_label: Option<Vec<u8>>,
    processed: i32,
}

impl ProtFeature {
    /// `CProt_ref::GetLabel`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Prot_ref.cpp:62-71
    fn get_label(&self) -> Vec<u8> {
        if let Some(first) = self.names.first() {
            return first.clone();
        }
        if let Some(desc) = &self.desc {
            return desc.clone();
        }
        if let Some(db) = &self.first_db_label {
            return db.clone();
        }
        Vec::new()
    }
}

#[derive(Clone, Debug)]
struct GeneFeature {
    seq_order: usize,
    feat_id: Option<i32>,
    loc: CompactLoc,
    gene: CompactGeneRef,
    qual_comment_fallback: Vec<u8>,
}

impl GeneFeature {
    /// `feature::GetLabel` with `fFGL_Content` for a `Gene` feature.
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/feature.cpp:598-600,710-776
    fn content_label(&self) -> Vec<u8> {
        let lbl = self.gene.get_label();
        if !lbl.is_empty() {
            return lbl;
        }
        if !self.qual_comment_fallback.is_empty() {
            return self.qual_comment_fallback.clone();
        }
        b"Gene".to_vec()
    }
}

#[derive(Clone, Debug)]
struct CdsFeature {
    seq_order: usize,
    loc: CompactLoc,
    product: Option<CompactLoc>,
    partial_is_set: bool,
    low_quality_except: bool,
    trans_splicing: bool,
    orf: bool,
    frame: i32,
    xref_prot_label: Option<Vec<u8>>,
    xref_gene: Option<CompactGeneRef>,
    feat_id_xrefs: Vec<i32>,
    qual_comment_fallback: Vec<u8>,
}

#[derive(Clone, Debug)]
struct BiosrcFeature {
    seq_order: usize,
    loc: CompactLoc,
    taxname: Option<Vec<u8>>,
    has_clone: bool,
    qual_comment_fallback: Vec<u8>,
}

/// Sorted interval index per Bioseq enabling `O(log N + k)` covering-feature lookup
/// (matching NCBI `CFeatTree` / `CFeat_CI` interval-tree complexity rather than `O(N^2)`).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/feature.cpp:2731-2920
#[derive(Default)]
struct SpanIndex {
    entries: Vec<SpanEntry>,
}

struct SpanEntry {
    min_from: u64,
    max_to: u64,
    prefix_max_to: u64,
    idx: usize,
}

impl SpanIndex {
    fn from_spans(mut entries: Vec<SpanEntry>) -> Self {
        entries.sort_unstable_by_key(|e| (e.min_from, e.idx));
        let mut running_max = 0u64;
        for e in &mut entries {
            if e.max_to > running_max {
                running_max = e.max_to;
            }
            e.prefix_max_to = running_max;
        }
        Self { entries }
    }

    fn iter_ids(&self) -> impl Iterator<Item = usize> + '_ {
        self.entries.iter().map(|e| e.idx)
    }

    /// Calls `f(feat_idx)` for every indexed feature whose total span on this Bioseq
    /// contains `[from, to_open)`.
    fn for_each_covering(&self, from: u64, to_open: u64, mut f: impl FnMut(usize)) {
        let end = self.entries.partition_point(|e| e.min_from <= from);
        for e in self.entries[..end].iter().rev() {
            if e.prefix_max_to < to_open {
                break;
            }
            if e.max_to >= to_open {
                f(e.idx);
            }
        }
    }
}

/// Compact per-TSE index of Bioseqs and defline-relevant features (`Prot`, `Cdregion`, `Gene`,
/// `Biosrc`), dropping all raw `Seq-annot` `Value` trees immediately during `collect_tree`.
#[derive(Default)]
pub(crate) struct TseFeatureIndex {
    pub(crate) interner: SeqIdInterner,
    pub(crate) id_to_bioseq: HashMap<HandleKey, usize>,
    pub(crate) bioseq_lengths: Vec<Option<u64>>,
    pub(crate) bioseq_seg_locs: HashMap<usize, CompactLoc>,
    next_order: usize,
    prots: Vec<ProtFeature>,
    genes: Vec<GeneFeature>,
    cds: Vec<CdsFeature>,
    biosrcs: Vec<BiosrcFeature>,
    prots_by_bioseq: HashMap<usize, Vec<usize>>,
    genes_by_bioseq: HashMap<usize, SpanIndex>,
    cds_by_bioseq: HashMap<usize, Vec<usize>>,
    cds_by_product: HashMap<usize, Vec<usize>>,
    biosrcs_by_bioseq: HashMap<usize, SpanIndex>,
    gene_by_feat_id: HashMap<i32, Vec<usize>>,
}

impl TseFeatureIndex {
    /// Registers a Bioseq's length and optional `Seg-ext` locations during `collect_tree`.
    pub(crate) fn record_bioseq(&mut self, bio_idx: usize, length: Option<u64>, inst: &Value) {
        if self.bioseq_lengths.len() <= bio_idx {
            self.bioseq_lengths.resize(bio_idx + 1, None);
        }
        self.bioseq_lengths[bio_idx] = length;
        let s = dschema();
        if let Some((ext_c, ext_v)) = inst.member(s.seq_inst_ext).and_then(as_choice) {
            if ext_c == s.seq_ext_seg {
                if let Some(items) = as_list(ext_v) {
                    let mut seg_loc = CompactLoc::default();
                    for item in items {
                        let sub = extract_loc(item, &mut self.interner);
                        seg_loc.pieces.extend(sub.pieces);
                    }
                    if !seg_loc.pieces.is_empty() {
                        self.bioseq_seg_locs.insert(bio_idx, seg_loc);
                    }
                }
            }
        }
    }

    /// Scans one `Seq-annot` list and extracts only the four feature types read by
    /// `CDeflineGenerator`, checking for `NULL` list elements as NCBI does.
    pub(crate) fn scan_annots(&mut self, val: &Value) -> Result<(), WriteDbError> {
        let Some(annots) = as_list(val) else {
            return Ok(());
        };
        let s = dschema();
        for annot in annots {
            if matches!(annot, Value::Null) {
                return Err(WriteDbError::Unsupported(
                    "a NULL element in Seq-annot list (NCBI makeblastdb dereferences a NULL pointer there)"
                        .to_string(),
                ));
            }
            let Some((choice, inner)) = annot.member(s.seq_annot_data).and_then(as_choice) else {
                continue;
            };
            if choice != s.seq_annot_ftable {
                continue;
            }
            let Some(feats) = as_list(inner) else {
                continue;
            };
            for feat in feats {
                if matches!(feat, Value::Null) {
                    return Err(WriteDbError::Unsupported(
                        "a NULL element in Seq-annot.data.ftable (NCBI makeblastdb dereferences a NULL pointer there)"
                            .to_string(),
                    ));
                }
                self.extract_feature(feat);
            }
        }
        Ok(())
    }

    fn extract_feature(&mut self, feat: &Value) {
        let s = dschema();
        let Some((fchoice, fval)) = feat.member(s.seq_feat_data).and_then(as_choice) else {
            return;
        };
        if fchoice != s.feat_prot
            && fchoice != s.feat_gene
            && fchoice != s.feat_cdregion
            && fchoice != s.feat_biosrc
        {
            return;
        }
        let seq_order = self.next_order;
        self.next_order += 1;
        let loc = feat
            .member(s.seq_feat_location)
            .map(|l| extract_loc(l, &mut self.interner))
            .unwrap_or_default();

        if fchoice == s.feat_prot {
            let names = fval
                .member(s.prot_name)
                .and_then(as_list)
                .map(|list| {
                    list.iter()
                        .filter_map(|v| match v {
                            Value::Str(b) => Some(b.clone()),
                            _ => None,
                        })
                        .collect()
                })
                .unwrap_or_default();
            let desc = match fval.member(s.prot_desc) {
                Some(Value::Str(b)) => Some(b.clone()),
                _ => None,
            };
            let first_activity = fval
                .member(s.prot_activity)
                .and_then(as_list)
                .and_then(|list| list.first())
                .and_then(|v| match v {
                    Value::Str(b) => Some(b.clone()),
                    _ => None,
                });
            let first_db_label = fval
                .member(s.prot_db)
                .and_then(as_list)
                .and_then(|list| list.first())
                .and_then(dbtag_label);
            let processed = match fval.member(s.prot_processed) {
                Some(Value::Enum(v)) => *v,
                _ => 0,
            };
            self.prots.push(ProtFeature {
                seq_order,
                loc,
                names,
                desc,
                first_activity,
                first_db_label,
                processed,
            });
        } else if fchoice == s.feat_gene {
            let feat_id = feat
                .member(s.seq_feat_id)
                .and_then(as_choice)
                .filter(|(c, _)| *c == s.feat_id_local)
                .and_then(|(_, v)| obj_id_int(v));
            let gene = extract_gene_ref(fval);
            let qual_comment_fallback = extract_qual_comment_fallback(feat);
            self.genes.push(GeneFeature {
                seq_order,
                feat_id,
                loc,
                gene,
                qual_comment_fallback,
            });
        } else if fchoice == s.feat_cdregion {
            let product = feat
                .member(s.seq_feat_product)
                .map(|p| extract_loc(p, &mut self.interner));
            let partial_is_set = feat.member(s.seq_feat_partial).is_some();
            let is_except = matches!(feat.member(s.seq_feat_except), Some(Value::Bool(true)));
            let except_text = match feat.member(s.seq_feat_except_text) {
                Some(Value::Str(t)) => Some(t.as_slice()),
                _ => None,
            };
            let low_quality_except =
                is_except && except_text.is_some_and(cds_except_is_low_quality);
            let trans_splicing = except_text.is_some_and(|t| contains_bytes(t, b"trans-splicing"));
            let orf = matches!(fval.member(s.cdregion_orf), Some(Value::Bool(true)));
            let frame = match fval.member(s.cdregion_frame) {
                Some(Value::Enum(v)) => *v,
                _ => 0,
            };
            let mut xref_prot_label = None;
            let mut xref_gene = None;
            let mut feat_id_xrefs = Vec::new();
            if let Some(xrefs) = feat.member(s.seq_feat_xref).and_then(as_list) {
                for xr in xrefs {
                    if let Some((ic, iv)) = xr.member(s.xref_id).and_then(as_choice) {
                        if ic == s.feat_id_local {
                            if let Some(id_num) = obj_id_int(iv) {
                                feat_id_xrefs.push(id_num);
                            }
                        }
                    }
                    if let Some((dc, dv)) = xr.member(s.xref_data).and_then(as_choice) {
                        if dc == s.feat_prot {
                            xref_prot_label = Some(prot_ref_label(dv));
                        } else if dc == s.feat_gene {
                            xref_gene = Some(extract_gene_ref(dv));
                        }
                    }
                }
            }
            let qual_comment_fallback = extract_qual_comment_fallback(feat);
            self.cds.push(CdsFeature {
                seq_order,
                loc,
                product,
                partial_is_set,
                low_quality_except,
                trans_splicing,
                orf,
                frame,
                xref_prot_label,
                xref_gene,
                feat_id_xrefs,
                qual_comment_fallback,
            });
        } else if fchoice == s.feat_biosrc {
            let taxname = fval
                .member(s.biosrc_org)
                .and_then(|o| o.member(s.org_ref_taxname))
                .and_then(|v| match v {
                    Value::Str(b) => Some(b.clone()),
                    _ => None,
                });
            let has_clone = fval
                .member(s.biosrc_subtype)
                .and_then(as_list)
                .is_some_and(|subs| {
                    subs.iter().any(|sbs| {
                        matches!(
                            sbs.member(s.subsrc_subtype),
                            Some(Value::Enum(subsource::CLONE))
                        ) && sbs.member(s.subsrc_name).is_some()
                    })
                });
            let qual_comment_fallback = extract_qual_comment_fallback(feat);
            self.biosrcs.push(BiosrcFeature {
                seq_order,
                loc,
                taxname,
                has_clone,
                qual_comment_fallback,
            });
        }
    }

    /// Finalizes inverted indices once all Bioseq IDs are mapped in `id_to_bioseq`.
    pub(crate) fn build_indices(&mut self) {
        if self.prots.is_empty()
            && self.genes.is_empty()
            && self.cds.is_empty()
            && self.biosrcs.is_empty()
        {
            return;
        }
        for (idx, feat) in self.prots.iter().enumerate() {
            index_loc_bioseqs(
                &feat.loc,
                &self.id_to_bioseq,
                &mut self.prots_by_bioseq,
                idx,
            );
        }
        let mut raw_genes: HashMap<usize, Vec<usize>> = HashMap::new();
        for (idx, feat) in self.genes.iter().enumerate() {
            if let Some(fid) = feat.feat_id {
                self.gene_by_feat_id.entry(fid).or_default().push(idx);
            }
            index_loc_bioseqs(&feat.loc, &self.id_to_bioseq, &mut raw_genes, idx);
        }
        for (bi, ids) in raw_genes {
            let entries = ids
                .into_iter()
                .map(|idx| {
                    let (min_from, max_to) = self
                        .loc_total_range_on(&self.genes[idx].loc, bi)
                        .unwrap_or((0, 0));
                    SpanEntry {
                        min_from,
                        max_to,
                        prefix_max_to: 0,
                        idx,
                    }
                })
                .collect();
            self.genes_by_bioseq
                .insert(bi, SpanIndex::from_spans(entries));
        }
        for (idx, feat) in self.cds.iter().enumerate() {
            index_loc_bioseqs(&feat.loc, &self.id_to_bioseq, &mut self.cds_by_bioseq, idx);
            if let Some(prod) = &feat.product {
                index_loc_bioseqs(prod, &self.id_to_bioseq, &mut self.cds_by_product, idx);
            }
        }
        let mut raw_biosrcs: HashMap<usize, Vec<usize>> = HashMap::new();
        for (idx, feat) in self.biosrcs.iter().enumerate() {
            index_loc_bioseqs(&feat.loc, &self.id_to_bioseq, &mut raw_biosrcs, idx);
        }
        for (bi, ids) in raw_biosrcs {
            let entries = ids
                .into_iter()
                .map(|idx| {
                    let (min_from, max_to) = self
                        .loc_total_range_on(&self.biosrcs[idx].loc, bi)
                        .unwrap_or((0, 0));
                    SpanEntry {
                        min_from,
                        max_to,
                        prefix_max_to: 0,
                        idx,
                    }
                })
                .collect();
            self.biosrcs_by_bioseq
                .insert(bi, SpanIndex::from_spans(entries));
        }
    }

    /// Resolves the ordered list of Bioseq indices targeted by `bio_idx` (itself, plus any
    /// segment Bioseqs when `repr == seg`).
    fn target_bioseqs(&self, bio_idx: usize) -> Vec<usize> {
        let mut out = vec![bio_idx];
        if let Some(seg_loc) = self.bioseq_seg_locs.get(&bio_idx) {
            for p in &seg_loc.pieces {
                if let Some(key) = &p.id_key {
                    if let Some(&seg_bio) = self.id_to_bioseq.get(key) {
                        if !out.contains(&seg_bio) {
                            out.push(seg_bio);
                        }
                    }
                }
            }
        }
        out
    }

    fn loc_length(&self, loc: &CompactLoc) -> u64 {
        let mut total = 0u64;
        for p in &loc.pieces {
            if p.is_whole {
                if let Some(len) = p
                    .id_key
                    .as_ref()
                    .and_then(|k| self.id_to_bioseq.get(k))
                    .and_then(|&bi| self.bioseq_lengths.get(bi).copied().flatten())
                {
                    total = total.saturating_add(len);
                }
            } else {
                total = total.saturating_add(p.to_open.saturating_sub(p.from));
            }
        }
        total
    }

    fn loc_total_range_on(&self, loc: &CompactLoc, bio_idx: usize) -> Option<(u64, u64)> {
        let mut min_from = u64::MAX;
        let mut max_to = 0u64;
        let mut any = false;
        for p in &loc.pieces {
            let Some(key) = &p.id_key else { continue };
            if self.id_to_bioseq.get(key).copied() != Some(bio_idx) {
                continue;
            }
            let (from, to_open) = if p.is_whole {
                let len = self
                    .bioseq_lengths
                    .get(bio_idx)
                    .copied()
                    .flatten()
                    .unwrap_or(0);
                (0, len)
            } else {
                (p.from, p.to_open)
            };
            min_from = min_from.min(from);
            max_to = max_to.max(to_open);
            any = true;
        }
        any.then_some((min_from, max_to))
    }

    fn loc_single_bioseq(&self, loc: &CompactLoc) -> Option<usize> {
        let mut bio = None;
        for p in &loc.pieces {
            let key = p.id_key.as_ref()?;
            let &bi = self.id_to_bioseq.get(key)?;
            match bio {
                None => bio = Some(bi),
                Some(prev) if prev != bi => return None,
                _ => {}
            }
        }
        bio
    }

    fn loc_ranges_by_rev(&self, loc: &CompactLoc) -> Vec<(usize, bool, u64, u64)> {
        let mut out = Vec::with_capacity(loc.pieces.len());
        for p in &loc.pieces {
            let Some(key) = &p.id_key else { continue };
            let Some(&bi) = self.id_to_bioseq.get(key) else {
                continue;
            };
            let (from, to_open) = if p.is_whole {
                let len = self.bioseq_lengths.get(bi).copied().flatten().unwrap_or(0);
                (0, len)
            } else {
                (p.from, p.to_open)
            };
            out.push((bi, is_reverse_strand(p.strand), from, to_open));
        }
        out
    }

    /// `s_SeqLocToTotalRangesInfoMapByStrand`: collapses consecutive same-Bioseq, same-orientation,
    /// properly-ordered pieces into total ranges per `(bio_idx, is_reverse)`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:1385-1476
    fn total_ranges_by_rev(&self, loc: &CompactLoc) -> Vec<(usize, bool, u64, u64)> {
        let pieces = self.loc_ranges_by_rev(loc);
        if pieces.is_empty() {
            return Vec::new();
        }
        let mut out = Vec::new();
        let (mut last_bi, mut last_rev, mut cur_from, mut cur_to) = pieces[0];
        let mut last_from = cur_from;
        for &(bi, rev, from, to_open) in &pieces[1..] {
            let break_range = bi != last_bi || rev != last_rev;
            let bad_order = !break_range
                && if rev {
                    from > last_from
                } else {
                    from < last_from
                };
            if break_range || bad_order {
                out.push((last_bi, last_rev, cur_from, cur_to));
                last_bi = bi;
                last_rev = rev;
                cur_from = from;
                cur_to = to_open;
            } else {
                cur_from = cur_from.min(from);
                cur_to = cur_to.max(to_open);
            }
            last_from = from;
        }
        out.push((last_bi, last_rev, cur_from, cur_to));
        out.sort_unstable();
        out
    }

    /// Single-Bioseq strand shortcut of `s_TestForOverlapEx`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:1874-1891
    fn single_seq_strand_incompatible(&self, outer: &CompactLoc, inner: &CompactLoc) -> bool {
        let single_seq =
            self.loc_single_bioseq(outer).is_some() && self.loc_single_bioseq(inner).is_some();
        single_seq
            && !test_strands(outer.strand, inner.strand)
            && outer.strand != 255
            && inner.strand != 255
    }

    /// `TestForOverlap64` under `eOverlap_Subset` (also used by `eOverlap_SubsetRev` after
    /// swapping arguments).
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:1535-1589,1874-1909
    fn test_for_overlap_subset(&self, outer: &CompactLoc, inner: &CompactLoc) -> Option<u64> {
        if self.single_seq_strand_incompatible(outer, inner) {
            return None;
        }
        let o_ranges = self.loc_ranges_by_rev(outer);
        let i_ranges = self.loc_ranges_by_rev(inner);
        if !ranges_contain_all(&o_ranges, &i_ranges) {
            return None;
        }
        Some(coverage_len(&o_ranges).saturating_sub(coverage_len(&i_ranges)))
    }

    /// `TestForOverlap64` under `eOverlap_Contained` plus the `Compare` (`fCompareOverlapping`)
    /// guard in `GetOverlappingFeatures`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:1036-1150,1655-1755,1874-1902
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:1144-1189
    fn test_for_overlap_contained(&self, outer: &CompactLoc, inner: &CompactLoc) -> Option<u64> {
        if self.single_seq_strand_incompatible(outer, inner) {
            return None;
        }
        let o_tot = self.total_ranges_by_rev(outer);
        let i_tot = self.total_ranges_by_rev(inner);
        if !ranges_contain_all(&o_tot, &i_tot) {
            return None;
        }
        // `Compare` with `fCompareOverlapping` requires `eContains` or `eSame`: every individual
        // interval of `inner` must be contained in some interval of `outer` on the same Bioseq.
        let o_int = self.loc_ranges_by_rev(outer);
        let i_int = self.loc_ranges_by_rev(inner);
        for &(bi, _, from, to_open) in &i_int {
            let contained = o_int
                .iter()
                .any(|&(obi, _, ofrom, oto)| obi == bi && ofrom <= from && oto >= to_open);
            if !contained {
                return None;
            }
        }
        Some(uncovered_len_by_rev(&o_tot, &i_tot))
    }

    fn sorted_bioseq_intervals(&self, loc: &CompactLoc) -> Vec<(usize, u64, u64)> {
        let mut out: Vec<(usize, u64, u64)> = self
            .loc_ranges_by_rev(loc)
            .into_iter()
            .map(|(bi, _, from, to)| (bi, from, to))
            .collect();
        out.sort_unstable();
        out
    }

    /// Whether `sequence::Compare` with `fCompareOverlapping` returns `eSame`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:1045-1072
    fn locs_same_intervals(&self, a: &CompactLoc, b: &CompactLoc) -> bool {
        self.sorted_bioseq_intervals(a) == self.sorted_bioseq_intervals(b)
    }

    /// Whether `CFeat_CI` with `SetOverlapTotalRange` collects `cand_loc` into
    /// `CFeatTree`'s `m_InfoMap` (shared Bioseq with intersecting total range and intersecting
    /// `CHandleRange::GetStrandsFlag`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/annot_collector.cpp:2737-2740
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/feature.cpp:2124-2130,3171-3207
    fn total_range_and_strand_intersects(
        &self,
        cds_loc: &CompactLoc,
        cand_loc: &CompactLoc,
    ) -> bool {
        if (handle_range_strands_flag(cds_loc) & handle_range_strands_flag(cand_loc)) == 0 {
            return false;
        }
        for p in &cds_loc.pieces {
            let Some(key) = &p.id_key else { continue };
            let Some(&bi) = self.id_to_bioseq.get(key) else {
                continue;
            };
            if let (Some((cf, ct)), Some((gf, gt))) = (
                self.loc_total_range_on(cds_loc, bi),
                self.loc_total_range_on(cand_loc, bi),
            ) {
                if gf < ct && cf < gt {
                    return true;
                }
            }
        }
        false
    }

    /// `s_GetCdregionLabel` + `feature::GetLabel` (`fFGL_Content`) for a `Cdregion` feature.
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/feature.cpp:139-247,710-776
    fn cdregion_content_label(&self, cds: &CdsFeature, warn: Warn<'_>) -> Vec<u8> {
        let mut tlabel = Vec::new();
        if let Some(prot_lbl) = &cds.xref_prot_label {
            tlabel.extend_from_slice(prot_lbl);
        } else {
            let mut found_from_product = false;
            if let Some(prod) = &cds.product {
                if let Some(key) = prod.pieces.first().and_then(|p| p.id_key.as_ref()) {
                    if let Some(&prod_bio) = self.id_to_bioseq.get(key) {
                        if let Some(prot_ids) = self.prots_by_bioseq.get(&prod_bio) {
                            let mut sorted = prot_ids.clone();
                            sorted.sort_by_key(|&pi| {
                                let f = &self.prots[pi];
                                let (from, to_open) =
                                    self.loc_total_range_on(&f.loc, prod_bio).unwrap_or((0, 0));
                                (from, std::cmp::Reverse(to_open), f.seq_order)
                            });
                            if let Some(&first_pi) = sorted.first() {
                                tlabel.extend_from_slice(&self.prots[first_pi].get_label());
                                found_from_product = true;
                            }
                        }
                    } else if let Some(fasta) = &prod.first_id_fasta {
                        let mut msg = b"Error: [makeblastdb] cannot find sequence: ".to_vec();
                        msg.extend_from_slice(fasta);
                        warn(&msg);
                    }
                }
            }
            if !found_from_product {
                if let Some(gref) = &cds.xref_gene {
                    tlabel.extend_from_slice(&gref.get_label());
                }
                if cds.orf {
                    tlabel.extend_from_slice(b"open reading frame: ");
                    let frame_str: &[u8] = match cds.frame {
                        1 => b"frame 1; ",
                        2 => b"frame 2; ",
                        3 => b"frame 3; ",
                        _ => b"frame not set; ",
                    };
                    tlabel.extend_from_slice(frame_str);
                    let strand_str: &[u8] = match cds.loc.strand {
                        1 => b"positive strand",
                        2 => b"negative strand",
                        3 => b"both strands",
                        4 => b"both strands (reverse)",
                        _ => b"strand unknown",
                    };
                    tlabel.extend_from_slice(strand_str);
                }
            }
        }
        if !tlabel.is_empty() {
            return tlabel;
        }
        if !cds.qual_comment_fallback.is_empty() {
            return cds.qual_comment_fallback.clone();
        }
        b"CDS".to_vec()
    }

    /// `CFeatTree::GetBestGene` for a `Cdregion` feature (`AddGenesForCds` + `x_LookupParentByRef`
    /// + `s_CollectBestOverlaps` with `PLessByEnd` and `CAnnotObject_Less` tie-breaking).
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/feature.cpp:1749-1759,2106-2142,2234-2250,2601-2700,3127-3143,3250-3256,3321-3337
    /// NCBI reference (598d8ae6): c++/src/objmgr/annot_collector.cpp:1054-1183,1375-1386,2737-2740
    fn best_gene_for_cds(&self, cds: &CdsFeature) -> Option<&GeneFeature> {
        // Explicit Feat-id xref only matches a gene that `AddGenesForCds` collected into `m_InfoMap`.
        for &fid in &cds.feat_id_xrefs {
            if let Some(gis) = self.gene_by_feat_id.get(&fid) {
                for &gi in gis {
                    let g = &self.genes[gi];
                    if self.total_range_and_strand_intersects(&cds.loc, &g.loc) {
                        return Some(g);
                    }
                }
            }
        }
        let cds_ranges = self.loc_ranges_by_rev(&cds.loc);
        let cds_bio = cds_ranges.first().map(|r| r.0)?;
        let gene_index = self.genes_by_bioseq.get(&cds_bio)?;
        let (cds_from, cds_to) = ranges_span_on(&cds_ranges, cds_bio)?;
        let cds_sflag = handle_range_strands_flag(&cds.loc);

        let mut best: Option<(i8, u64, u64, u64, GeneAnnotTieKey, usize)> = None;
        gene_index.for_each_covering(cds_from, cds_to, |gi| {
            let g = &self.genes[gi];
            if (cds_sflag & handle_range_strands_flag(&g.loc)) == 0 {
                return;
            }
            let (neg_quality, diff) =
                if let Some(d) = self.test_for_overlap_subset(&g.loc, &cds.loc) {
                    (0i8, d)
                } else if cds.trans_splicing && cds.loc.strand == 255 && g.loc.strand != 255 {
                    let pstrand = if g.loc.strand == 0 { 1 } else { g.loc.strand };
                    let forced = force_loc_strand(&cds.loc, pstrand);
                    let Some(d) = self.test_for_overlap_subset(&g.loc, &forced) else {
                        return;
                    };
                    (1i8, d)
                } else {
                    return;
                };
            let (g_from, g_to) = self.loc_total_range_on(&g.loc, cds_bio).unwrap_or((0, 0));
            let tie = gene_annot_tie_key(g);
            let is_better = match &best {
                None => true,
                Some((bq, bd, bto, bfrom, btie, _)) => {
                    (neg_quality, diff, g_to, g_from, &tie) < (*bq, *bd, *bto, *bfrom, btie)
                }
            };
            if is_better {
                best = Some((neg_quality, diff, g_to, g_from, tie, gi));
            }
        });
        best.map(|(_, _, _, _, _, gi)| &self.genes[gi])
    }

    /// Best overlapping `Biosrc` feature under `eOverlap_SubsetRev` with reverse-strand fallback
    /// (`GetBestOverlappingFeat`, `x_GetSourceFeatViaCDS`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2044-2096
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:908-943,964-974,1092,1158-1205
    /// NCBI reference (598d8ae6): c++/src/objmgr/annot_collector.cpp:1054-1183,1375-1386,2700-2735
    fn best_biosrc_for_cds(&self, cds: &CdsFeature) -> Option<&BiosrcFeature> {
        let cds_ranges = self.loc_ranges_by_rev(&cds.loc);
        let cds_bio = cds_ranges.first().map(|r| r.0)?;
        let src_index = self.biosrcs_by_bioseq.get(&cds_bio)?;
        let (cds_from, cds_to) = ranges_span_on(&cds_ranges, cds_bio)?;

        let find_for = |target: &CompactLoc| -> Option<&BiosrcFeature> {
            let mut best: Option<(u64, u64, std::cmp::Reverse<u64>, BiosrcAnnotTieKey, usize)> =
                None;
            src_index.for_each_covering(cds_from, cds_to, |si| {
                let src = &self.biosrcs[si];
                if !cfeat_ci_intervals_strand_match(target, &src.loc) {
                    return;
                }
                let Some(diff) = self.test_for_overlap_subset(&src.loc, target) else {
                    return;
                };
                let (s_from, s_to) = self.loc_total_range_on(&src.loc, cds_bio).unwrap_or((0, 0));
                let tie = biosrc_annot_tie_key(src);
                let is_better = match &best {
                    None => true,
                    Some((bd, bfrom, bto, btie, _)) => {
                        (diff, s_from, std::cmp::Reverse(s_to), &tie) < (*bd, *bfrom, *bto, btie)
                    }
                };
                if is_better {
                    best = Some((diff, s_from, std::cmp::Reverse(s_to), tie, si));
                }
            });
            best.map(|(_, _, _, _, si)| &self.biosrcs[si])
        };

        if let Some(found) = find_for(&cds.loc) {
            return Some(found);
        }
        let rev_loc = reverse_loc(&cds.loc);
        find_for(&rev_loc)
    }

    /// Best overlapping `Gene` feature under `eOverlap_Contained` (`GetBestOverlappingFeat`,
    /// `x_GetSegSeqInfoViaCDS`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:2621-2638
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:908-943,958-962,1144-1205
    /// NCBI reference (598d8ae6): c++/src/objmgr/annot_collector.cpp:1054-1183,1375-1386,2737-2740
    fn best_contained_gene(&self, cds: &CdsFeature) -> Option<&GeneFeature> {
        let cds_ranges = self.loc_ranges_by_rev(&cds.loc);
        let cds_bio = cds_ranges.first().map(|r| r.0)?;
        let gene_index = self.genes_by_bioseq.get(&cds_bio)?;
        let (cds_from, cds_to) = ranges_span_on(&cds_ranges, cds_bio)?;
        let cds_sflag = handle_range_strands_flag(&cds.loc);

        let mut cands: Vec<(u64, u64, std::cmp::Reverse<u64>, GeneAnnotTieKey, usize)> = Vec::new();
        gene_index.for_each_covering(cds_from, cds_to, |gi| {
            let g = &self.genes[gi];
            if (cds_sflag & handle_range_strands_flag(&g.loc)) == 0 {
                return;
            }
            let Some(diff) = self.test_for_overlap_contained(&g.loc, &cds.loc) else {
                return;
            };
            let (g_from, g_to) = self.loc_total_range_on(&g.loc, cds_bio).unwrap_or((0, 0));
            cands.push((
                diff,
                g_from,
                std::cmp::Reverse(g_to),
                gene_annot_tie_key(g),
                gi,
            ));
        });
        if cands.is_empty() {
            return None;
        }
        // Order into `CFeat_CI` (`CAnnotObject_Less`) order first, then stable-sort by
        // `COverlapPairLess` (diff ASC, and `CGene_ref::GetLabel` ASC on `eSame` locations).
        cands.sort_by(|a, b| (&a.1, &a.2, &a.3).cmp(&(&b.1, &b.2, &b.3)));
        cands.sort_by(|a, b| {
            if a.0 != b.0 {
                return a.0.cmp(&b.0);
            }
            let ga = &self.genes[a.4];
            let gb = &self.genes[b.4];
            if self.locs_same_intervals(&ga.loc, &gb.loc) {
                return ga.gene.get_label().cmp(&gb.gene.get_label());
            }
            std::cmp::Ordering::Equal
        });
        cands.first().map(|c| &self.genes[c.4])
    }
}

/// `CAnnotObjectType_Less` tie-break key for `Gene` features with equal `(from, to_open)`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/annot_collector.cpp:1054-1183
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct GeneAnnotTieKey {
    is_rev: bool,
    complex_tie: bool,
    locus_lower: Vec<u8>,
    desc_lower: Vec<u8>,
    has_qc: bool,
    qc_lower: Vec<u8>,
    seq_order: usize,
}

fn gene_annot_tie_key(g: &GeneFeature) -> GeneAnnotTieKey {
    let is_rev = is_reverse_strand(g.loc.strand);
    GeneAnnotTieKey {
        is_rev,
        complex_tie: g.loc.is_complex ^ is_rev,
        locus_lower: g
            .gene
            .locus
            .as_deref()
            .unwrap_or_default()
            .to_ascii_lowercase(),
        desc_lower: g
            .gene
            .desc
            .as_deref()
            .unwrap_or_default()
            .to_ascii_lowercase(),
        has_qc: !g.qual_comment_fallback.is_empty(),
        qc_lower: g.qual_comment_fallback.to_ascii_lowercase(),
        seq_order: g.seq_order,
    }
}

/// `CAnnotObjectType_Less` tie-break key for `Biosrc` features with equal `(from, to_open)`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/annot_collector.cpp:1054-1091,1153-1183
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct BiosrcAnnotTieKey {
    is_rev: bool,
    complex_tie: bool,
    has_qc: bool,
    qc_lower: Vec<u8>,
    seq_order: usize,
}

fn biosrc_annot_tie_key(s: &BiosrcFeature) -> BiosrcAnnotTieKey {
    let is_rev = is_reverse_strand(s.loc.strand);
    BiosrcAnnotTieKey {
        is_rev,
        complex_tie: s.loc.is_complex ^ is_rev,
        has_qc: !s.qual_comment_fallback.is_empty(),
        qc_lower: s.qual_comment_fallback.to_ascii_lowercase(),
        seq_order: s.seq_order,
    }
}

/// `IsReverse` for `ENa_strand`: true for `minus` (2) and `both-rev` (4).
///
/// NCBI reference (598d8ae6): c++/include/objects/seqloc/Na_strand.hpp:75-79
fn is_reverse_strand(s: i32) -> bool {
    matches!(s, 2 | 4)
}

/// `Reverse` for `ENa_strand` used by `SeqLocRevCmpl`.
///
/// NCBI reference (598d8ae6): c++/include/objects/seqloc/Na_strand.hpp:90-105
fn reverse_strand(s: i32) -> i32 {
    match s {
        0 | 1 => 2,
        2 => 1,
        3 => 4,
        4 => 3,
        other => other,
    }
}

/// `s_Test_Strands` in `TestForOverlap64`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:1246-1261
fn test_strands(s1: i32, s2: i32) -> bool {
    if s1 == 255 || s2 == 255 {
        return false;
    }
    s1 == s2 || s1 == 3 || s2 == 3 || (s1 == 0 && s2 != 2) || (s2 == 0 && s1 != 2)
}

/// `CHandleRange::GetStrandsFlag`: bit 1 (`eStrandPlus`) when any piece is not `minus`,
/// bit 2 (`eStrandMinus`) when any piece is `unknown`/`minus`/`both`/`both-rev`.
///
/// NCBI reference (598d8ae6): c++/include/objmgr/impl/handle_range.hpp:185-199
/// NCBI reference (598d8ae6): c++/src/objmgr/handle_range.cpp:81-106
fn handle_range_strands_flag(loc: &CompactLoc) -> u8 {
    let mut flags = 0u8;
    for p in &loc.pieces {
        if p.strand != 2 {
            flags |= 1;
        }
        if matches!(p.strand, 0 | 2 | 3 | 4) {
            flags |= 2;
        }
    }
    flags
}

/// `CAnnot_Collector::x_MatchRange` strand check for `eOverlap_Intervals` (`GetBestOverlappingFeat`
/// with `eOverlap_SubsetRev`).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/annot_collector.cpp:2700-2735
/// NCBI reference (598d8ae6): c++/src/objmgr/handle_range.cpp:219-228
fn cfeat_ci_intervals_strand_match(query: &CompactLoc, cand: &CompactLoc) -> bool {
    let q_strand = query.strand;
    if cand.pieces.len() <= 1 {
        let cand_s = match handle_range_strands_flag(cand) {
            1 => 1,
            2 => 2,
            _ => 0,
        };
        return q_strand == 0 || cand_s == 0 || q_strand == cand_s;
    }
    cand.pieces
        .iter()
        .any(|p| q_strand == 0 || p.strand == 0 || q_strand == p.strand)
}

/// `SeqLocRevCmpl` on a `CompactLoc`: flips every piece's strand via `Reverse` and reverses order.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:970-974
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_interval.cpp:244-251
fn reverse_loc(loc: &CompactLoc) -> CompactLoc {
    let mut pieces: Vec<LocPiece> = loc
        .pieces
        .iter()
        .rev()
        .map(|p| {
            let mut q = p.clone();
            if !q.is_whole {
                q.strand = reverse_strand(q.strand);
            }
            q
        })
        .collect();
    // Preserve `whole` strand (`both` -> `both`, `FlipStrand` does not alter `whole`).
    let strand = if loc.pieces.len() == 1 && loc.pieces[0].is_whole {
        loc.strand
    } else {
        reverse_strand(loc.strand)
    };
    if pieces.is_empty() {
        pieces = Vec::new();
    }
    CompactLoc {
        pieces,
        strand,
        is_complex: loc.is_complex,
        first_id_fasta: loc.first_id_fasta.clone(),
    }
}

/// Forces all pieces and top strand of `loc` to `strand` (matching `SetStrand` on the cloned
/// child location for trans-spliced features).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/feature.cpp:2676-2684
fn force_loc_strand(loc: &CompactLoc, strand: i32) -> CompactLoc {
    let pieces = loc
        .pieces
        .iter()
        .map(|p| {
            let mut q = p.clone();
            q.strand = strand;
            q
        })
        .collect();
    CompactLoc {
        pieces,
        strand,
        is_complex: loc.is_complex,
        first_id_fasta: loc.first_id_fasta.clone(),
    }
}

/// `sequence::GetCoverage`: total union length of intervals per Bioseq (ignoring strand).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:1908-1909
fn coverage_len(ranges: &[(usize, bool, u64, u64)]) -> u64 {
    let mut sorted: Vec<(usize, u64, u64)> = ranges
        .iter()
        .map(|&(bi, _, from, to_open)| (bi, from, to_open))
        .collect();
    sorted.sort_unstable();
    let mut total = 0u64;
    let mut cur: Option<(usize, u64, u64)> = None;
    for (bi, from, to_open) in sorted {
        match cur {
            Some((cbi, cf, ct)) if cbi == bi && from <= ct => {
                cur = Some((cbi, cf, ct.max(to_open)));
            }
            Some((_, cf, ct)) => {
                total = total.saturating_add(ct.saturating_sub(cf));
                cur = Some((bi, from, to_open));
            }
            None => cur = Some((bi, from, to_open)),
        }
    }
    if let Some((_, cf, ct)) = cur {
        total = total.saturating_add(ct.saturating_sub(cf));
    }
    total
}

/// `s_GetUncoveredLength`: sum of subranges of `outer` not covered by `inner` on the same
/// `(bio_idx, is_reverse)`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:1479-1532
fn uncovered_len_by_rev(
    outer: &[(usize, bool, u64, u64)],
    inner: &[(usize, bool, u64, u64)],
) -> u64 {
    let mut diff = 0u64;
    for &(bi, rev, mut from, to_open) in outer {
        for &(ibi, irev, ifrom, ito) in inner {
            if ibi != bi || irev != rev || ito <= from {
                continue;
            }
            if ifrom >= to_open {
                break;
            }
            if ifrom > from {
                diff = diff.saturating_add(ifrom - from);
            }
            if ito < to_open {
                from = ito;
            } else {
                from = to_open;
                break;
            }
        }
        diff = diff.saturating_add(to_open.saturating_sub(from));
    }
    diff
}

/// Checks that `outer` and `inner` are non-empty and every `(bio, rev, from, to_open)` in `inner`
/// lies inside some range in `outer` on the same Bioseq and strand orientation.
fn ranges_contain_all(
    outer: &[(usize, bool, u64, u64)],
    inner: &[(usize, bool, u64, u64)],
) -> bool {
    !outer.is_empty()
        && !inner.is_empty()
        && inner.iter().all(|&(bi, rev, from, to_open)| {
            outer.iter().any(|&(obi, orev, ofrom, oto)| {
                obi == bi && orev == rev && ofrom <= from && oto >= to_open
            })
        })
}

fn ranges_span_on(ranges: &[(usize, bool, u64, u64)], bio_idx: usize) -> Option<(u64, u64)> {
    let mut min_from = u64::MAX;
    let mut max_to = 0u64;
    let mut any = false;
    for &(bi, _, from, to_open) in ranges {
        if bi == bio_idx {
            min_from = min_from.min(from);
            max_to = max_to.max(to_open);
            any = true;
        }
    }
    any.then_some((min_from, max_to))
}

fn index_loc_bioseqs(
    loc: &CompactLoc,
    id_to_bioseq: &HashMap<HandleKey, usize>,
    map: &mut HashMap<usize, Vec<usize>>,
    feat_idx: usize,
) {
    let mut added = Vec::new();
    for p in &loc.pieces {
        let Some(key) = &p.id_key else { continue };
        if let Some(&bi) = id_to_bioseq.get(key) {
            if !added.contains(&bi) {
                added.push(bi);
                map.entry(bi).or_default().push(feat_idx);
            }
        }
    }
}

/// Extracts `Dbtag::GetLabel` from a `Dbtag` `Value`.
///
/// NCBI reference (598d8ae6): c++/src/objects/general/Dbtag.cpp:187-209
fn dbtag_label(val: &Value) -> Option<Vec<u8>> {
    let s = dschema();
    let Value::Str(db) = val.member(s.dbtag_db)? else {
        return None;
    };
    let tag = val.member(s.dbtag_tag)?;
    if let Some(str_val) = obj_id_str(tag) {
        if str_val.len() > db.len()
            && str_val[db.len()] == b':'
            && str_val[..db.len()].eq_ignore_ascii_case(db)
        {
            return Some(str_val.to_vec());
        }
        let mut out = Vec::with_capacity(db.len() + 2 + str_val.len());
        out.extend_from_slice(db);
        out.extend_from_slice(b": ");
        out.extend_from_slice(str_val);
        return Some(out);
    }
    if let Some(id_val) = obj_id_int(tag) {
        let mut out = db.clone();
        out.extend_from_slice(b": ");
        out.extend_from_slice(id_val.to_string().as_bytes());
        return Some(out);
    }
    Some(db.clone())
}

fn extract_gene_ref(val: &Value) -> CompactGeneRef {
    let s = dschema();
    let get_str = |idx: usize| match val.member(idx) {
        Some(Value::Str(b)) => Some(b.clone()),
        _ => None,
    };
    let first_syn = val
        .member(s.gene_syn)
        .and_then(as_list)
        .and_then(|list| list.first())
        .and_then(|v| match v {
            Value::Str(b) => Some(b.clone()),
            _ => None,
        });
    let first_db_label = val
        .member(s.gene_db)
        .and_then(as_list)
        .and_then(|list| list.first())
        .and_then(dbtag_label);
    CompactGeneRef {
        locus: get_str(s.gene_locus),
        desc: get_str(s.gene_desc),
        locus_tag: get_str(s.gene_locus_tag),
        first_syn,
        first_db_label,
        maploc: get_str(s.gene_maploc),
    }
}

fn prot_ref_label(val: &Value) -> Vec<u8> {
    let s = dschema();
    if let Some(Value::Str(first)) = val
        .member(s.prot_name)
        .and_then(as_list)
        .and_then(|l| l.first())
    {
        return first.clone();
    }
    if let Some(Value::Str(desc)) = val.member(s.prot_desc) {
        return desc.clone();
    }
    if let Some(lbl) = val
        .member(s.prot_db)
        .and_then(as_list)
        .and_then(|l| l.first())
        .and_then(dbtag_label)
    {
        return lbl;
    }
    Vec::new()
}

/// Extracts `Seq-feat.qual` + `Seq-feat.comment` fallback label (`feature.cpp:717-738`).
fn extract_qual_comment_fallback(feat: &Value) -> Vec<u8> {
    let s = dschema();
    let mut out = Vec::new();
    if let Some(quals) = feat.member(s.seq_feat_qual).and_then(as_list) {
        let mut prefix: &[u8] = b"/";
        for q in quals {
            if let Some(Value::Str(qname)) = q.member(s.gb_qual_qual) {
                out.extend_from_slice(prefix);
                out.extend_from_slice(qname);
                prefix = b" /";
                if let Some(Value::Str(qval)) = q.member(s.gb_qual_val) {
                    if !qval.is_empty() {
                        out.push(b'=');
                        out.extend_from_slice(qval);
                    }
                }
            }
        }
    }
    if let Some(Value::Str(comment)) = feat.member(s.seq_feat_comment) {
        if out.is_empty() {
            out.extend_from_slice(comment);
        } else {
            out.extend_from_slice(b"; ");
            out.extend_from_slice(comment);
        }
    }
    out
}

/// Extracts a `CompactLoc` iteratively from a `Seq-loc` `Value`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:88-115,481-513
fn extract_loc(root: &Value, interner: &mut SeqIdInterner) -> CompactLoc {
    let s = dschema();
    let mut out = CompactLoc {
        is_complex: matches!(as_choice(root), Some((c, _)) if c == s.seqloc_mix || c == s.seqloc_packed_int),
        ..CompactLoc::default()
    };

    let mut stack = vec![root];
    while let Some(loc) = stack.pop() {
        let Some((choice, inner)) = as_choice(loc) else {
            continue;
        };
        if choice == s.seqloc_whole {
            if let Value::SeqId(sid) = inner {
                let (key, fasta) = interner.key_and_fasta(sid);
                if out.first_id_fasta.is_none() {
                    out.first_id_fasta = Some(fasta);
                }
                out.pieces.push(LocPiece {
                    id_key: key,
                    from: 0,
                    to_open: 0,
                    is_whole: true,
                    strand: 0,
                });
            }
        } else if choice == s.seqloc_int {
            push_seq_int(inner, interner, &mut out);
        } else if choice == s.seqloc_packed_int {
            if let Some(ints) = as_list(inner) {
                for iv in ints {
                    push_seq_int(iv, interner, &mut out);
                }
            }
        } else if choice == s.seqloc_pnt {
            push_seq_pnt(inner, interner, &mut out);
        } else if choice == s.seqloc_packed_pnt {
            let strand = match inner.member(s.packed_pnt_strand) {
                Some(Value::Enum(v)) => *v,
                _ => 0,
            };
            let key = match inner.member(s.packed_pnt_id) {
                Some(Value::SeqId(sid)) => {
                    let (k, fasta) = interner.key_and_fasta(sid);
                    if out.first_id_fasta.is_none() {
                        out.first_id_fasta = Some(fasta);
                    }
                    k
                }
                _ => None,
            };
            if let Some(pts) = inner.member(s.packed_pnt_points) {
                match pts {
                    Value::List(items) => {
                        for p in items {
                            if let Value::Int(pt) = p {
                                let pos = *pt as u64;
                                out.pieces.push(LocPiece {
                                    id_key: key.clone(),
                                    from: pos,
                                    to_open: pos.saturating_add(1),
                                    is_whole: false,
                                    strand,
                                });
                            }
                        }
                    }
                    Value::Ints(items) => {
                        for &pt in items {
                            let pos = pt as u64;
                            out.pieces.push(LocPiece {
                                id_key: key.clone(),
                                from: pos,
                                to_open: pos.saturating_add(1),
                                is_whole: false,
                                strand,
                            });
                        }
                    }
                    Value::Int4s(items) => {
                        for &pt in items {
                            let pos = pt as u64;
                            out.pieces.push(LocPiece {
                                id_key: key.clone(),
                                from: pos,
                                to_open: pos.saturating_add(1),
                                is_whole: false,
                                strand,
                            });
                        }
                    }
                    _ => {}
                }
            }
        } else if choice == s.seqloc_mix || choice == s.seqloc_equiv {
            if let Some(items) = as_list(inner) {
                for item in items.iter().rev() {
                    stack.push(item);
                }
            }
        } else if choice == s.seqloc_bond {
            if let Some(a) = inner.member(s.seqbond_a) {
                push_seq_pnt(a, interner, &mut out);
            }
            if let Some(b) = inner.member(s.seqbond_b) {
                push_seq_pnt(b, interner, &mut out);
            }
        }
    }
    out.strand = compute_top_strand(root, &out.pieces);
    out
}

fn push_seq_int(iv: &Value, interner: &mut SeqIdInterner, out: &mut CompactLoc) {
    let s = dschema();
    let from = match iv.member(s.seqint_from) {
        Some(Value::Int(v)) => *v as u64,
        _ => 0,
    };
    let to = match iv.member(s.seqint_to) {
        Some(Value::Int(v)) => *v as u64,
        _ => from,
    };
    let strand = match iv.member(s.seqint_strand) {
        Some(Value::Enum(v)) => *v,
        _ => 0,
    };
    let key = match iv.member(s.seqint_id) {
        Some(Value::SeqId(sid)) => {
            let (k, fasta) = interner.key_and_fasta(sid);
            if out.first_id_fasta.is_none() {
                out.first_id_fasta = Some(fasta);
            }
            k
        }
        _ => None,
    };
    out.pieces.push(LocPiece {
        id_key: key,
        from,
        to_open: to.saturating_add(1),
        is_whole: false,
        strand,
    });
}

fn push_seq_pnt(pt: &Value, interner: &mut SeqIdInterner, out: &mut CompactLoc) {
    let s = dschema();
    let pos = match pt.member(s.seqpnt_point) {
        Some(Value::Int(v)) => *v as u64,
        _ => 0,
    };
    let strand = match pt.member(s.seqpnt_strand) {
        Some(Value::Enum(v)) => *v,
        _ => 0,
    };
    let key = match pt.member(s.seqpnt_id) {
        Some(Value::SeqId(sid)) => {
            let (k, fasta) = interner.key_and_fasta(sid);
            if out.first_id_fasta.is_none() {
                out.first_id_fasta = Some(fasta);
            }
            k
        }
        _ => None,
    };
    out.pieces.push(LocPiece {
        id_key: key,
        from: pos,
        to_open: pos.saturating_add(1),
        is_whole: false,
        strand,
    });
}

/// `sequence::GetStrand` + `s_GetStrand` for ORF CDS labels and `TestForOverlap64`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:393-513
fn compute_top_strand(loc: &Value, pieces: &[LocPiece]) -> i32 {
    let s = dschema();
    let Some((choice, inner)) = as_choice(loc) else {
        return 0;
    };
    if choice == s.seqloc_int {
        return match inner.member(s.seqint_strand) {
            Some(Value::Enum(v)) => *v,
            _ => 0,
        };
    }
    if choice == s.seqloc_whole {
        return 3; // eNa_strand_both
    }
    if choice == s.seqloc_pnt {
        return match inner.member(s.seqpnt_strand) {
            Some(Value::Enum(v)) => *v,
            _ => 0,
        };
    }
    if choice == s.seqloc_packed_pnt {
        return match inner.member(s.packed_pnt_strand) {
            Some(Value::Enum(v)) => *v,
            _ => 0,
        };
    }
    // `GetStrand` checks `IsOneBioseq` before calling `s_GetStrand`.
    if !pieces_are_one_id(pieces) {
        return 0;
    }
    eval_s_get_strand(loc)
}

fn pieces_are_one_id(pieces: &[LocPiece]) -> bool {
    let mut first: Option<&HandleKey> = None;
    for p in pieces {
        let Some(k) = p.id_key.as_ref() else {
            return false;
        };
        match first {
            None => first = Some(k),
            Some(prev) if prev != k => return false,
            _ => {}
        }
    }
    first.is_some()
}

/// Iterative evaluation of `s_GetStrand` (`seq_loc_util.cpp:393-477`).
fn eval_s_get_strand(root: &Value) -> i32 {
    let s = dschema();
    let mut stack: Vec<(&[Value], usize, i32, bool)> = Vec::new();
    let mut cur = root;
    loop {
        let result = match as_choice(cur) {
            Some((c, inner)) if c == s.seqloc_int => match inner.member(s.seqint_strand) {
                Some(Value::Enum(v)) => *v,
                _ => 0,
            },
            Some((c, _)) if c == s.seqloc_whole => 3,
            Some((c, inner)) if c == s.seqloc_pnt => match inner.member(s.seqpnt_strand) {
                Some(Value::Enum(v)) => *v,
                _ => 0,
            },
            Some((c, inner)) if c == s.seqloc_packed_pnt => {
                match inner.member(s.packed_pnt_strand) {
                    Some(Value::Enum(v)) => *v,
                    _ => 0,
                }
            }
            Some((c, inner)) if c == s.seqloc_bond => {
                let mut a_strand = inner
                    .member(s.seqbond_a)
                    .and_then(|a| a.member(s.seqpnt_strand))
                    .and_then(|v| match v {
                        Value::Enum(x) => Some(*x),
                        _ => None,
                    })
                    .unwrap_or(0);
                let mut b_strand = inner
                    .member(s.seqbond_b)
                    .and_then(|b| b.member(s.seqpnt_strand))
                    .and_then(|v| match v {
                        Value::Enum(x) => Some(*x),
                        _ => None,
                    })
                    .unwrap_or(0);
                if a_strand == 0 && b_strand != 0 {
                    a_strand = b_strand;
                } else if a_strand != 0 && b_strand == 0 {
                    b_strand = a_strand;
                }
                if a_strand != b_strand {
                    255
                } else {
                    a_strand
                }
            }
            Some((c, inner)) if c == s.seqloc_packed_int => {
                let mut strand = 0i32;
                let mut strand_set = false;
                if let Some(ints) = as_list(inner) {
                    for iv in ints {
                        let istrand = match iv.member(s.seqint_strand) {
                            Some(Value::Enum(x)) => *x,
                            _ => 0,
                        };
                        if !advance_strand_state(&mut strand, &mut strand_set, istrand) {
                            strand = 255;
                            break;
                        }
                    }
                }
                strand
            }
            Some((c, inner)) if c == s.seqloc_mix => {
                if let Some(items) = as_list(inner) {
                    if let Some((next_idx, next_val)) = next_non_empty_mix_item(items, 0, s) {
                        stack.push((items, next_idx + 1, 0, false));
                        cur = next_val;
                        continue;
                    }
                }
                0
            }
            _ => 0,
        };

        let mut val = result;
        loop {
            let Some((items, next_idx, mut strand, mut strand_set)) = stack.pop() else {
                return val;
            };
            if !advance_strand_state(&mut strand, &mut strand_set, val) {
                val = 255;
                continue;
            }
            if let Some((idx, next_val)) = next_non_empty_mix_item(items, next_idx, s) {
                stack.push((items, idx + 1, strand, strand_set));
                cur = next_val;
                break;
            }
            val = strand;
        }
    }
}

fn next_non_empty_mix_item<'a>(
    items: &'a [Value],
    start: usize,
    s: &DeflineSchema,
) -> Option<(usize, &'a Value)> {
    for (offset, item) in items[start..].iter().enumerate() {
        if matches!(as_choice(item), Some((c, _)) if c == s.seqloc_null || c == s.seqloc_empty) {
            continue;
        }
        return Some((start + offset, item));
    }
    None
}

/// Exact `s_GetStrand` state transition (`seq_loc_util.cpp:435-446, 459-470`).
/// Returns `false` when the strand becomes `eNa_strand_other` (255).
fn advance_strand_state(strand: &mut i32, strand_set: &mut bool, istrand: i32) -> bool {
    if *strand == 0 && istrand == 1 {
        *strand = 1;
        *strand_set = true;
    } else if *strand == 1 && istrand == 0 {
        *strand_set = true;
    } else if !*strand_set {
        *strand = istrand;
        *strand_set = true;
    } else if istrand != *strand {
        return false;
    }
    true
}

/// Counts `CSeqMap::eSeqGap` segments in `Seq-inst` (`CSeqMap_CI` with `fFindGap`).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/seq_map.cpp:219-292,1359-1435
fn count_seq_map_gaps(inst: &Value, repr: i32, length: u64) -> u32 {
    let s = dschema();
    if let Some(sd) = inst.member(s.seq_inst_seq_data) {
        return u32::from(matches!(as_choice(sd), Some((c, _)) if c == s.seq_data_gap));
    }
    if let Some((ec, ev)) = inst.member(s.seq_inst_ext).and_then(as_choice) {
        if ec == s.seq_ext_delta {
            let mut gaps = 0u32;
            if let Some(items) = as_list(ev) {
                for item in items {
                    let Some((dc, dv)) = as_choice(item) else {
                        continue;
                    };
                    if dc == s.delta_literal {
                        let is_gap = match dv.member(s.seq_lit_data) {
                            None => true,
                            Some(sd) => {
                                matches!(as_choice(sd), Some((c, _)) if c == s.seq_data_gap)
                            }
                        };
                        if is_gap {
                            gaps = gaps.saturating_add(1);
                        }
                    } else if dc == s.delta_loc {
                        gaps = gaps.saturating_add(count_loc_gaps(dv));
                    }
                }
            }
            return gaps;
        }
        if ec == s.seq_ext_seg {
            let mut gaps = 0u32;
            if let Some(items) = as_list(ev) {
                for loc in items {
                    gaps = gaps.saturating_add(count_loc_gaps(loc));
                }
            }
            return gaps;
        }
        if ec == s.seq_ext_map {
            return 1;
        }
    }
    if repr == seq_repr::VIRTUAL || repr == seq_repr::MAP || (repr != 0 && length != 0) {
        return 1;
    }
    0
}

fn count_loc_gaps(loc: &Value) -> u32 {
    let s = dschema();
    let mut gaps = 0u32;
    let mut stack = vec![loc];
    while let Some(cur) = stack.pop() {
        let Some((c, v)) = as_choice(cur) else {
            gaps = gaps.saturating_add(1);
            continue;
        };
        if c == s.seqloc_null || c == s.seqloc_empty {
            gaps = gaps.saturating_add(1);
        } else if c == s.seqloc_mix || c == s.seqloc_equiv {
            if let Some(items) = as_list(v) {
                for item in items {
                    stack.push(item);
                }
            }
        }
    }
    gaps
}

/// Extracts the last restriction enzyme string from `Seq-inst.ext.map` (`x_SetFlags` lines 1152-1168).
fn extract_map_r_enzyme(inst: &Value) -> Option<&[u8]> {
    let s = dschema();
    let (ec, ev) = inst.member(s.seq_inst_ext).and_then(as_choice)?;
    if ec != s.seq_ext_map {
        return None;
    }
    let mut enzyme = None;
    for feat in as_list(ev)? {
        if let Some((fc, fv)) = feat.member(s.seq_feat_data).and_then(as_choice) {
            if fc == s.feat_rsite {
                if let Some((rc, Value::Str(text))) = as_choice(fv) {
                    if rc == s.rsite_str {
                        enzyme = Some(text.as_slice());
                    }
                }
            }
        }
    }
    enzyme
}

/// Feature query adapter bridging `TseFeatureIndex` to `report::defline::flags::DeflineFeatures`.
struct BioseqFeatureView<'a, 'w> {
    index: &'a TseFeatureIndex,
    bio_idx: usize,
    warn: RefCell<Warn<'w>>,
}

impl DeflineFeatures for BioseqFeatureView<'_, '_> {
    fn has_biosrc_clone_feature(&self) -> bool {
        for bi in self.index.target_bioseqs(self.bio_idx) {
            if let Some(ids) = self.index.biosrcs_by_bioseq.get(&bi) {
                if ids.iter_ids().any(|si| self.index.biosrcs[si].has_clone) {
                    return true;
                }
            }
        }
        false
    }

    fn nm_single_gene_and_cds(&self) -> Option<(Vec<u8>, Vec<u8>)> {
        let targets = self.index.target_bioseqs(self.bio_idx);
        let mut gene_idx = None;
        let mut gene_count = 0usize;
        let mut seen_genes = Vec::new();
        for &bi in &targets {
            if let Some(ids) = self.index.genes_by_bioseq.get(&bi) {
                for gi in ids.iter_ids() {
                    if !seen_genes.contains(&gi) {
                        seen_genes.push(gi);
                        gene_count += 1;
                        gene_idx = Some(gi);
                    }
                }
            }
        }
        let mut cds_idx = None;
        let mut cds_count = 0usize;
        let mut seen_cds = Vec::new();
        for &bi in &targets {
            if let Some(ids) = self.index.cds_by_bioseq.get(&bi) {
                for &ci in ids {
                    if !seen_cds.contains(&ci) {
                        seen_cds.push(ci);
                        cds_count += 1;
                        cds_idx = Some(ci);
                    }
                }
            }
        }
        if gene_count != 1 || cds_count != 1 {
            return None;
        }
        let cds = &self.index.cds[cds_idx?];
        let gene = &self.index.genes[gene_idx?];
        let mut warn = self.warn.borrow_mut();
        let cds_label = self.index.cdregion_content_label(cds, *warn);
        let gene_label = gene.content_label();
        Some((cds_label, gene_label))
    }

    fn nr_first_gene_label(&self) -> Option<Vec<u8>> {
        let targets = self.index.target_bioseqs(self.bio_idx);
        let mut candidates: Vec<(u64, std::cmp::Reverse<u64>, usize, usize)> = Vec::new();
        let mut seen = Vec::new();
        for &bi in &targets {
            if let Some(ids) = self.index.genes_by_bioseq.get(&bi) {
                for gi in ids.iter_ids() {
                    if seen.contains(&gi) {
                        continue;
                    }
                    seen.push(gi);
                    let f = &self.index.genes[gi];
                    let (from, to_open) =
                        self.index.loc_total_range_on(&f.loc, bi).unwrap_or((0, 0));
                    candidates.push((from, std::cmp::Reverse(to_open), f.seq_order, gi));
                }
            }
        }
        candidates.sort();
        let &(_, _, _, gi) = candidates.first()?;
        Some(self.index.genes[gi].content_label())
    }

    fn longest_protein(&self) -> Option<ProteinFeatInfo> {
        let targets = self.index.target_bioseqs(self.bio_idx);
        let mut candidates: Vec<(u64, std::cmp::Reverse<u64>, usize, usize)> = Vec::new();
        let mut seen = Vec::new();
        for &bi in &targets {
            if let Some(ids) = self.index.prots_by_bioseq.get(&bi) {
                for &pi in ids {
                    if seen.contains(&pi) {
                        continue;
                    }
                    seen.push(pi);
                    let f = &self.index.prots[pi];
                    let (from, to_open) =
                        self.index.loc_total_range_on(&f.loc, bi).unwrap_or((0, 0));
                    candidates.push((from, std::cmp::Reverse(to_open), f.seq_order, pi));
                }
            }
        }
        candidates.sort();

        let mut longest = 0u64;
        let mut best_processed = 0i32;
        let mut best_pi: Option<usize> = None;
        for (_, _, _, pi) in candidates {
            let f = &self.index.prots[pi];
            let plen = self.index.loc_length(&f.loc);
            if plen > longest
                || (plen == longest && best_pi.is_some() && f.processed < best_processed)
            {
                best_pi = Some(pi);
                longest = plen;
                best_processed = f.processed;
            }
        }
        let pi = best_pi?;
        let f = &self.index.prots[pi];
        Some(ProteinFeatInfo {
            names: f.names.clone(),
            desc: f.desc.clone(),
            first_activity: f.first_activity.clone(),
        })
    }

    fn cds_for_product(&self) -> Option<CdsForProductInfo> {
        let &ci = self.index.cds_by_product.get(&self.bio_idx)?.first()?;
        let cds = &self.index.cds[ci];
        let gene = cds
            .xref_gene
            .as_ref()
            .or_else(|| self.index.best_gene_for_cds(cds).map(|g| &g.gene));
        Some(CdsForProductInfo {
            has_low_quality_exception: cds.low_quality_except,
            gene_locus: gene.and_then(|g| g.locus.clone()),
            gene_first_syn: gene.and_then(|g| g.first_syn.clone()),
            gene_desc: gene.and_then(|g| g.desc.clone()),
            gene_locus_tag: gene.and_then(|g| g.locus_tag.clone()),
        })
    }

    fn source_taxname_via_cds(&self) -> Option<Vec<u8>> {
        let &ci = self.index.cds_by_product.get(&self.bio_idx)?.first()?;
        let cds = &self.index.cds[ci];
        let src = self.index.best_biosrc_for_cds(cds)?;
        src.taxname.clone()
    }

    fn segseq_info_via_cds(&self) -> Option<SegSeqCdsInfo> {
        let seg_loc = self.index.bioseq_seg_locs.get(&self.bio_idx)?;
        let mut candidates: Vec<(u64, std::cmp::Reverse<u64>, usize, usize)> = Vec::new();
        let mut seen = Vec::new();
        for p in &seg_loc.pieces {
            let Some(key) = &p.id_key else { continue };
            let Some(&bi) = self.index.id_to_bioseq.get(key) else {
                continue;
            };
            if let Some(ids) = self.index.cds_by_bioseq.get(&bi) {
                for &ci in ids {
                    if seen.contains(&ci) {
                        continue;
                    }
                    seen.push(ci);
                    let f = &self.index.cds[ci];
                    let (from, to_open) =
                        self.index.loc_total_range_on(&f.loc, bi).unwrap_or((0, 0));
                    candidates.push((from, std::cmp::Reverse(to_open), f.seq_order, ci));
                }
            }
        }
        candidates.sort();
        let &(_, _, _, ci) = candidates.first()?;
        let cds = &self.index.cds[ci];
        let mut warn = self.warn.borrow_mut();
        let product = self.index.cdregion_content_label(cds, *warn);
        let mut locus = if let Some(g) = &cds.xref_gene {
            g.locus
                .clone()
                .or_else(|| g.first_syn.clone())
                .unwrap_or_default()
        } else {
            Vec::new()
        };
        if locus.is_empty() {
            if let Some(gene) = self.index.best_contained_gene(cds) {
                locus = gene.content_label();
            }
        }
        Some(SegSeqCdsInfo {
            product,
            locus,
            is_partial: cds.partial_is_set,
        })
    }
}

/// Runs full `CDeflineGenerator::GenerateDefline` over a typed ASN.1 Bioseq, its ancestor
/// descriptors, and the TSE feature index, refusing only `AutoDef` (`D3c-6`).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:612-1434,3928-4107
pub(crate) fn generate_asn_bioseq_title(
    bio_idx: usize,
    ids: &[ParsedSeqId],
    bio_descr: &[Value],
    ancestors: &[&[Value]],
    inst: &Value,
    is_aa: bool,
    is_na: bool,
    repr: i32,
    length: u64,
    index: &TseFeatureIndex,
    warn: Warn<'_>,
) -> Result<Vec<u8>, WriteDbError> {
    let s = dschema();
    let is_seg = repr == seq_repr::SEG;
    let is_delta = repr == seq_repr::DELTA;
    let is_virtual = repr == seq_repr::VIRTUAL;
    let is_map = repr == seq_repr::MAP;
    let topo = match inst.member(s.seq_inst_topology) {
        Some(Value::Enum(v)) => *v,
        _ => topology::LINEAR,
    };

    let sid = classify_seq_ids_full(ids, is_virtual, &mut |w| {
        post_seq_id_warning(w, warn);
    });

    let mut inp = DeflineInput {
        is_na,
        is_aa,
        topology: topo,
        length,
        is_seg,
        is_delta,
        is_virtual,
        is_map,
        is_nc: sid.is_nc,
        is_nm: sid.is_nm,
        is_nr: sid.is_nr,
        is_nz: sid.is_nz,
        is_wp: sid.is_wp,
        third_party: sid.third_party,
        wgs_master: sid.wgs_master,
        tsa_master: sid.tsa_master,
        general_str: &sid.general_str,
        general_id: sid.general_id,
        is_pdb: sid.is_pdb,
        pdb_chain_id: &sid.pdb_chain_id,
        pdb_chain: sid.pdb_chain,
        is_patent: sid.is_patent,
        patent_sequence: sid.patent_sequence,
        patent_country: &sid.patent_country,
        patent_number: &sid.patent_number,
        ..DeflineInput::default()
    };

    let has_autodef = scan_descriptors_into(bio_descr, ancestors, &mut inp);

    if is_map {
        if let Some(enz) = extract_map_r_enzyme(inst) {
            inp.r_enzyme = enz;
        }
    }
    inp.use_pdb_compound_for_defline =
        compute_use_pdb_compound(inp.is_pdb, inp.is_na, inp.length, inp.comment);
    if is_delta {
        inp.delta_gap_count = count_seq_map_gaps(inst, repr, length);
    }

    // D3c-6 boundary: refuse only when NCBI's `GenerateDefline` actually enters `CAutoDef`
    // (create_defline.cpp:3974-3991).
    let main_empty_after_trim =
        inp.reconstruct || trim_end_slice(inp.main_title, b".,;~ ").is_empty();
    if main_empty_after_trim && inp.is_na && has_autodef && !inp.is_tls && !inp.is_nz {
        return Err(WriteDbError::Unsupported(
            "AutoDef title generation is not supported by LOSAT's makeblastdb".to_string(),
        ));
    }

    let view = BioseqFeatureView {
        index,
        bio_idx,
        warn: RefCell::new(warn),
    };

    let defline = generate_defline_core(&inp, Some(&view)).map_err(|_| {
        // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4530-4533
        WriteDbError::Exception(NcbiException::new(
            ThrowSite {
                file: "src/corelib/ncbistr.cpp",
                line: 4531,
                kind: "CStringException::eBadArgs",
                function: Some("ncbi::NStr::HtmlDecode"),
            },
            "Unable to guess the source string encoding",
        ))
    })?;

    if defline.reads_past_end {
        return Err(WriteDbError::Unsupported(
            "a title of punctuation that NCBI's x_CleanAndCompress reads past the end of \
             (NCBI makeblastdb crashes)"
                .to_string(),
        ));
    }

    Ok(defline.text)
}

/// Populates `inp` from Bioseq and ancestor descriptors in `CSeqdesc_CI` order, returning whether
/// any `AutodefOptions` user object (case-insensitive) is present.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:840-1150,1264-1406
fn scan_descriptors_into<'a>(
    bio_descr: &'a [Value],
    ancestors: &[&'a [Value]],
    inp: &mut DeflineInput<'a>,
) -> bool {
    const F_MOLINFO: u8 = 1 << 0;
    const F_USER: u8 = 1 << 1;
    const F_SOURCE: u8 = 1 << 2;
    const F_GENBANK: u8 = 1 << 3;
    const F_EMBL: u8 = 1 << 4;
    const F_TITLE: u8 = 1 << 5;
    const F_PDB: u8 = 1 << 6;

    let mut needed = F_MOLINFO | F_USER | F_SOURCE | F_GENBANK | F_EMBL | F_TITLE;
    if inp.is_pdb {
        needed |= F_PDB;
    }

    let s = dschema();
    let mut keywords: Option<&'a [Value]> = None;
    let mut num_sk = 0usize;
    let mut sk_diff = false;
    let mut has_autodef = false;

    for (level_idx, level) in std::iter::once(bio_descr)
        .chain(ancestors.iter().copied())
        .enumerate()
    {
        let on_bioseq = level_idx == 0;
        for (choice, inner) in level.iter().filter_map(as_choice) {
            if choice == s.seqdesc_molinfo {
                if needed & F_MOLINFO == 0 {
                    continue;
                }
                inp.mi_biomol = match inner.member(s.molinfo_biomol) {
                    Some(Value::Enum(v)) => *v,
                    _ => 0,
                };
                inp.mi_tech = match inner.member(s.molinfo_tech) {
                    Some(Value::Enum(v)) => *v,
                    _ => 0,
                };
                inp.mi_completeness = match inner.member(s.molinfo_completeness) {
                    Some(Value::Enum(v)) => *v,
                    _ => 0,
                };
                apply_molinfo_tech(inp);
                if inp.htgs_unfinished {
                    needed &= !F_TITLE;
                }
                needed &= !F_MOLINFO;
            } else if choice == s.seqdesc_user {
                if let Some(ty) = inner.member(s.user_obj_type).and_then(obj_id_str) {
                    if ty.eq_ignore_ascii_case(b"AutodefOptions") {
                        has_autodef = true;
                    }
                    if needed & F_USER != 0 {
                        if ty.eq_ignore_ascii_case(b"Unverified") {
                            inp.is_unverified = true;
                            needed &= !F_USER;
                            let (mut feat, mut mis, mut cont, mut org) =
                                (false, false, false, false);
                            if let Some(fields) = inner.member(s.user_obj_data).and_then(as_list) {
                                for fld in fields {
                                    if fld.member(s.user_field_label).and_then(obj_id_str)
                                        != Some(b"Type")
                                    {
                                        continue;
                                    }
                                    if let Some((dc, Value::Str(val))) =
                                        fld.member(s.user_field_data).and_then(as_choice)
                                    {
                                        if dc == s.user_field_str {
                                            match val.as_slice() {
                                                b"Features" => feat = true,
                                                b"Misassembled" => mis = true,
                                                b"Contaminant" => cont = true,
                                                b"Organism" => org = true,
                                                _ => {}
                                            }
                                        }
                                    }
                                }
                            }
                            inp.unverified_prefix =
                                unverified_prefix_from_types(feat, mis, cont, org);
                        } else if ty.eq_ignore_ascii_case(b"Unreviewed") {
                            inp.is_unreviewed = true;
                            let mut unannot = false;
                            if let Some(fields) = inner.member(s.user_obj_data).and_then(as_list) {
                                for fld in fields {
                                    if fld.member(s.user_field_label).and_then(obj_id_str)
                                        != Some(b"Type")
                                    {
                                        continue;
                                    }
                                    if let Some((dc, Value::Str(val))) =
                                        fld.member(s.user_field_data).and_then(as_choice)
                                    {
                                        if dc == s.user_field_str
                                            && val.as_slice() == b"Unannotated"
                                        {
                                            unannot = true;
                                        }
                                    }
                                }
                            }
                            inp.unreviewed_prefix = unreviewed_prefix_from_unannotated(unannot);
                        } else if ty == b"AutodefOptions" {
                            if let Some(fields) = inner.member(s.user_obj_data).and_then(as_list) {
                                for fld in fields {
                                    let Some(lbl) =
                                        fld.member(s.user_field_label).and_then(obj_id_str)
                                    else {
                                        continue;
                                    };
                                    if !lbl.eq_ignore_ascii_case(b"Targeted Locus Name") {
                                        continue;
                                    }
                                    if let Some((dc, Value::Str(val))) =
                                        fld.member(s.user_field_data).and_then(as_choice)
                                    {
                                        if dc == s.user_field_str {
                                            inp.targeted_locus = val.as_slice();
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else if choice == s.seqdesc_comment {
                if let Value::Str(c) = inner {
                    inp.comment = c.as_slice();
                    if contains_bytes(c, b"[CAUTION] Could be the product of a pseudogene") {
                        inp.is_pseudogene = true;
                    }
                }
            } else if choice == s.seqdesc_source {
                if needed & F_SOURCE != 0 {
                    needed &= !F_SOURCE;
                    inp.has_source = true;
                    populate_primary_biosource(inner, inp);
                }
                scan_biosource_for_mag_and_kingdom(inner, inp, &mut num_sk, &mut sk_diff);
            } else if choice == s.seqdesc_title {
                if needed & F_TITLE != 0 {
                    if inp.is_pdb || on_bioseq {
                        if let Value::Str(t) = inner {
                            inp.main_title = t.as_slice();
                        }
                    }
                    needed &= !F_TITLE;
                }
            } else if choice == s.seqdesc_genbank {
                if needed & F_GENBANK != 0 {
                    if let Some(kw) = inner.member(s.gb_block_keywords).and_then(as_list) {
                        keywords = Some(kw);
                    }
                    needed &= !(F_GENBANK | F_EMBL);
                }
            } else if choice == s.seqdesc_embl {
                if needed & F_EMBL != 0 {
                    if let Some(kw) = inner.member(s.embl_block_keywords).and_then(as_list) {
                        keywords = Some(kw);
                    }
                    needed &= !F_EMBL;
                }
            } else if choice == s.seqdesc_pdb && (needed & F_PDB != 0) {
                if let Some(compounds) = inner.member(s.pdb_block_compound).and_then(as_list) {
                    for cp in compounds {
                        if inp.pdb_compound.is_empty() {
                            if let Value::Str(text) = cp {
                                inp.pdb_compound = text.as_slice();
                                needed &= !F_PDB;
                            }
                        }
                    }
                }
            }
        }
    }

    if let Some(kw_list) = keywords {
        apply_keyword_clauses(
            inp,
            kw_list.iter().filter_map(|v| match v {
                Value::Str(b) => Some(b.as_slice()),
                _ => None,
            }),
        );
    }

    has_autodef
}

fn populate_primary_biosource<'a>(bsrc: &'a Value, inp: &mut DeflineInput<'a>) {
    let s = dschema();
    if let Some(Value::Enum(g)) = bsrc.member(s.biosrc_genome) {
        inp.genome = *g;
        inp.is_plasmid = *g == genome::PLASMID;
        inp.is_chromosome = *g == genome::CHROMOSOME;
    }
    if let Some(subs) = bsrc.member(s.biosrc_subtype).and_then(as_list) {
        for sbs in subs {
            let Some(Value::Enum(sub)) = sbs.member(s.subsrc_subtype) else {
                continue;
            };
            let Some(Value::Str(name)) = sbs.member(s.subsrc_name) else {
                continue;
            };
            match *sub {
                subsource::CHROMOSOME => inp.chromosome = name.as_slice(),
                subsource::CLONE => {
                    inp.clone = name.as_slice();
                    inp.has_clone = true;
                }
                subsource::MAP => inp.map = name.as_slice(),
                subsource::PLASMID_NAME => inp.plasmid = name.as_slice(),
                subsource::SEGMENT => inp.segment = name.as_slice(),
                subsource::LINKAGE_GROUP => inp.linkage_group = name.as_slice(),
                _ => {}
            }
        }
    }
    let Some(org) = bsrc.member(s.biosrc_org) else {
        return;
    };
    if let Some(Value::Str(taxname)) = org.member(s.org_ref_taxname) {
        inp.taxname = taxname.as_slice();
    }
    let Some(orgname) = org.member(s.org_ref_orgname) else {
        return;
    };
    if let Some((choice, inner)) = orgname.member(s.orgname_name).and_then(as_choice) {
        if choice == s.orgname_binomial {
            if let Some(Value::Str(g)) = inner.member(s.binomial_genus) {
                inp.genus = g.as_slice();
            }
            if let Some(Value::Str(sp)) = inner.member(s.binomial_species) {
                inp.species = sp.as_slice();
            }
        } else if choice == s.orgname_partial {
            if let Some(elems) = as_list(inner) {
                for te in elems {
                    let Some(Value::Enum(fl)) = te.member(s.taxelem_fixed_level) else {
                        continue;
                    };
                    if *fl > 0 {
                        inp.multispecies_from_source = true;
                    } else if let Some(Value::Str(lvl)) = te.member(s.taxelem_level) {
                        if !lvl.eq_ignore_ascii_case(b"species") {
                            inp.multispecies_from_source = true;
                        }
                    }
                }
            }
        }
    }
    if let Some(mods) = orgname.member(s.orgname_mod).and_then(as_list) {
        for m in mods {
            let Some(Value::Enum(sub)) = m.member(s.orgmod_subtype) else {
                continue;
            };
            let Some(Value::Str(val)) = m.member(s.orgmod_subname) else {
                continue;
            };
            match *sub {
                orgmod::STRAIN if inp.strain.is_empty() => inp.strain = val.as_slice(),
                orgmod::SUBSTRAIN if inp.substrain.is_empty() => inp.substrain = val.as_slice(),
                orgmod::CULTIVAR if inp.cultivar.is_empty() => inp.cultivar = val.as_slice(),
                orgmod::SPECIMEN_VOUCHER if inp.specimen_voucher.is_empty() => {
                    inp.specimen_voucher = val.as_slice();
                }
                orgmod::ISOLATE if inp.isolate.is_empty() => inp.isolate = val.as_slice(),
                orgmod::BREED if inp.breed.is_empty() => inp.breed = val.as_slice(),
                orgmod::METAGENOME_SOURCE if inp.metagenome_source.is_empty() => {
                    inp.metagenome_source = val.as_slice();
                }
                _ => {}
            }
        }
    }
}

fn scan_biosource_for_mag_and_kingdom<'a>(
    bsrc: &'a Value,
    inp: &mut DeflineInput<'a>,
    num_sk: &mut usize,
    sk_diff: &mut bool,
) {
    let s = dschema();
    let Some(orgname) = bsrc
        .member(s.biosrc_org)
        .and_then(|o| o.member(s.org_ref_orgname))
    else {
        return;
    };
    if let Some(om_list) = orgname.member(s.orgname_mod).and_then(as_list) {
        for m in om_list {
            if matches!(
                m.member(s.orgmod_subtype),
                Some(Value::Enum(orgmod::METAGENOME_SOURCE))
            ) {
                if let Some(Value::Str(val)) = m.member(s.orgmod_subname) {
                    if inp.metagenome_source.is_empty() {
                        inp.metagenome_source = val.as_slice();
                    }
                }
            }
        }
    }
    if inp.is_wp {
        if let Some((choice, partial_val)) = orgname.member(s.orgname_name).and_then(as_choice) {
            if choice == s.orgname_partial {
                if let Some(tax_elems) = as_list(partial_val) {
                    for te in tax_elems {
                        if !matches!(te.member(s.taxelem_fixed_level), Some(Value::Enum(0))) {
                            continue;
                        }
                        let Some(Value::Str(lvl)) = te.member(s.taxelem_level) else {
                            continue;
                        };
                        if !lvl.eq_ignore_ascii_case(b"superkingdom") {
                            continue;
                        }
                        *num_sk += 1;
                        if let Some(Value::Str(name)) = te.member(s.taxelem_name) {
                            if inp.first_super_kingdom.is_empty() {
                                inp.first_super_kingdom = name.as_slice();
                            } else if !inp.first_super_kingdom.eq_ignore_ascii_case(name)
                                && inp.second_super_kingdom.is_empty()
                            {
                                *sk_diff = true;
                                inp.second_super_kingdom = name.as_slice();
                            }
                        }
                        if *num_sk > 1 && *sk_diff {
                            inp.is_cross_kingdom = true;
                        }
                    }
                }
            }
        }
    }
}

fn contains_bytes(haystack: &[u8], needle: &[u8]) -> bool {
    haystack.windows(needle.len()).any(|w| w == needle)
}

fn trim_end_slice<'a>(bytes: &'a [u8], chars: &[u8]) -> &'a [u8] {
    match bytes.iter().rposition(|b| !chars.contains(b)) {
        Some(pos) => &bytes[..=pos],
        None => bytes,
    }
}
