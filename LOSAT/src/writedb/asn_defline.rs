//! Typed ASN.1 Bioseq descriptor/extension extraction and compact TSE feature index for
//! `CDeflineGenerator` in `losat makeblastdb -input_type asn1_txt|asn1_bin`.
//!
//! NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:612-1434,1916-2096,2578-2645,3928-4107
//! NCBI reference (598d8ae6): c++/src/objmgr/util/feature.cpp:103-248,586-776,2106-2142,2731-3340
//! NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:1060-1275,2548-2589
//! NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:88-115,300-513,1653-1850
//! NCBI reference (598d8ae6): c++/include/objmgr/impl/tse_info.hpp:335-346
//! NCBI reference (598d8ae6): c++/src/objmgr/tse_info.cpp:1514-1531,1644-1654
//! NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:1177-1253
//! NCBI reference (598d8ae6): c++/src/objmgr/seq_map.cpp:219-292,1359-1435

use std::cell::RefCell;
use std::collections::HashMap;
use std::hash::Hasher;
use std::sync::OnceLock;

use rustc_hash::{FxHashMap, FxHasher};

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
use crate::seqid::defline_set::{compare, object_id_number, Sic};
use crate::seqid::label;
use crate::seqid::mask_handle::{
    handle_tree, mask_handle_key, pdb_normalized, textseq_pack_prefix, textseq_packs, GiimTree,
    TreeHandle,
};
use crate::seqid::model::{Date, Dbtag, ObjectId, ParsedSeqId, PdbSeqId, TextseqId};

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
    seq_feat_ids: usize,
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
            seq_feat_ids: member_index(named::SEQ_FEAT, "ids"),
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

/// `CObject_id` key (`id` or `str`) for `CTSE_Info` local Feat-id lookup.
///
/// NCBI reference (598d8ae6): c++/include/objmgr/impl/tse_info.hpp:335-346
/// NCBI reference (598d8ae6): c++/src/objmgr/tse_info.cpp:1514-1531,1644-1654
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum FeatObjId {
    Id(i32),
    Str(Vec<u8>),
}

fn feat_obj_id(val: &Value) -> Option<FeatObjId> {
    let s = dschema();
    match as_choice(val) {
        Some((c, Value::Int(n))) if c == s.object_id_id => Some(FeatObjId::Id(*n as i32)),
        Some((c, Value::Str(b))) if c == s.object_id_str => Some(FeatObjId::Str(b.clone())),
        _ => None,
    }
}

fn local_feat_id(val: &Value) -> Option<FeatObjId> {
    let s = dschema();
    as_choice(val)
        .filter(|(c, _)| *c == s.feat_id_local)
        .and_then(|(_, v)| feat_obj_id(v))
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

/// Index of a Seq-id of the TSE's feature and segment locations in `SeqIdInterner`.
type IdRef = u32;

/// The Seq-ids of the TSE's feature and segment locations with their `CSeq_id_Handle` keys
/// (`None` for a Textseq-id with neither accession nor name). One entry per distinct Seq-id; a
/// `giim` id or one without a key gets an entry per occurrence, since every
/// `CSeq_id_Handle::GetHandle` call may give it another handle (`mask_handle_key`, `GiimTree`).
#[derive(Default)]
pub(crate) struct SeqIdInterner {
    giim: GiimTree,
    ids: Vec<ParsedSeqId>,
    keys: Vec<Option<HandleKey>>,
    by_id: FxHashMap<ParsedSeqId, IdRef>,
    last: Option<IdRef>,
    /// The Seq-ids of the location being extracted, in `CSeq_loc_CI` order (kept for reuse).
    ci: Vec<IdRef>,
}

impl SeqIdInterner {
    fn key_only(&mut self, id: &ParsedSeqId) -> Option<HandleKey> {
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

    /// The entry of `id`, made on its first occurrence (its handle key then taken in the order
    /// NCBI's object manager takes it while it indexes the TSE).
    fn intern(&mut self, id: &ParsedSeqId) -> IdRef {
        if let Some(last) = self.last {
            let at = last as usize;
            if matches!(self.keys[at], Some(HandleKey::Bytes(_))) && self.ids[at] == *id {
                return last;
            }
        }
        if !matches!(id, ParsedSeqId::Giim(_)) {
            if let Some(&found) = self.by_id.get(id) {
                self.last = Some(found);
                return found;
            }
        }
        let key = self.key_only(id);
        let new = IdRef::try_from(self.ids.len()).expect("fewer than 2^32 Seq-ids in one TSE");
        if matches!(key, Some(HandleKey::Bytes(_))) {
            self.by_id.insert(id.clone(), new);
        }
        self.ids.push(id.clone());
        self.keys.push(key);
        self.last = Some(new);
        new
    }

    fn id(&self, r: IdRef) -> &ParsedSeqId {
        &self.ids[r as usize]
    }

    fn key(&self, r: IdRef) -> Option<&HandleKey> {
        self.keys[r as usize].as_ref()
    }

    /// `CSeq_id::Equals` (`CSeq_loc::GetId` compares every Seq-id with the first so).
    fn equals(&self, a: IdRef, b: IdRef) -> bool {
        a == b || self.id(a) == self.id(b)
    }

    /// `IsSameBioseq` without a scope: `CSeq_id::Match`, or one `CSeq_id_Handle`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:300-321
    fn same_free(&self, a: IdRef, b: IdRef) -> bool {
        let key = self.key(a);
        (a == b && key.is_some())
            || compare(self.id(a), self.id(b)) == Sic::Yes
            || (key.is_some() && key == self.key(b))
    }

    /// Whether `s_GetId` (or `CSeq_loc::GetId`) without a scope finds one Bioseq for the Seq-ids
    /// `ids` of a location in `CSeq_loc_CI` order: `IsOneBioseq` with no scope.
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:325-357,385-389
    fn one_free(&self, ids: &[IdRef]) -> bool {
        let Some((&first, rest)) = ids.split_first() else {
            return false;
        };
        rest.iter().all(|&x| self.equals(first, x))
            || rest.iter().all(|&x| self.same_free(first, x))
    }

    /// `CSeq_loc::GetId` and `s_GetId` of a location whose Seq-ids in `CSeq_loc_CI` order are
    /// `ids`: the first Seq-id, whether every Seq-id `Equals` it, and the Seq-ids that only the
    /// scope can find on its Bioseq.
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:325-357
    /// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:3541-3630
    fn loc_ids(&self, ids: &[IdRef]) -> (Option<IdRef>, bool, LocIds) {
        let Some((&first, rest)) = ids.split_first() else {
            return (None, false, LocIds::None);
        };
        let mut exact = true;
        let mut scoped: Vec<IdRef> = Vec::new();
        for &x in rest {
            exact = exact && self.equals(first, x);
            if !self.same_free(first, x) {
                scoped.push(x);
            }
        }
        let form = if scoped.is_empty() {
            LocIds::Free
        } else {
            LocIds::Scoped(Box::new(ScopedIds {
                strand: 0,
                ids: scoped.into_boxed_slice(),
            }))
        };
        (Some(first), exact, form)
    }
}

/// One interval/whole/point piece of a `Seq-loc`.
#[derive(Clone, Debug)]
struct LocPiece {
    id: Option<IdRef>,
    from: u64,
    /// Exclusive upper bound (`to + 1`) when `is_whole` is false.
    to_open: u64,
    is_whole: bool,
    /// `ENa_strand` (0 unknown, 1 plus, 2 minus, 3 both, 4 both-rev, 255 other).
    strand: i32,
}

/// How `s_GetId` finds the Seq-ids of a location on one Bioseq.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:300-357
#[derive(Clone, Debug, Default)]
enum LocIds {
    /// The location has no Seq-id (`s_GetId` finds none).
    #[default]
    None,
    /// Every Seq-id `Match`es the first or has its handle: one Bioseq without a scope.
    Free,
    /// Some Seq-ids only the scope can find on the first one's Bioseq (`build_indices`).
    Scoped(Box<ScopedIds>),
}

/// The Seq-ids of a location that neither `Match` its first Seq-id nor have its handle, and the
/// `s_GetStrand` that `GetStrand` gives when the scope finds them on its Bioseq.
#[derive(Clone, Debug)]
struct ScopedIds {
    strand: i32,
    ids: Box<[IdRef]>,
}

/// Compact representation of a `Seq-loc` extracted during the TSE walk.
#[derive(Clone, Debug, Default)]
pub(crate) struct CompactLoc {
    pieces: Vec<LocPiece>,
    /// `sequence::GetStrand` with the TSE scope (`s_GetCdregionLabel` ORF strand); set in
    /// `build_indices` for `LocIds::Scoped`.
    strand: i32,
    /// Whether the top-level `Seq-loc` is `mix` or `packed-int` (`CAnnotObjectType_Less`).
    is_complex: bool,
    /// Whether every Seq-id `Equals` the first (`CSeq_loc::GetId` is set).
    exact: bool,
    /// The first Seq-id in `CSeq_loc_CI` order, `empty` sub-locations included (`s_GetId`).
    first_id: Option<IdRef>,
    ids: LocIds,
    /// Whether `strand` waits for the picks of ties (`TseFeatureIndex::loc_strand`).
    deferred: Deferred,
}

impl CompactLoc {
    /// `sequence::GetStrand` without a scope (`s_TestForOverlapEx`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:481-513,1874-1891
    fn strand_free(&self) -> i32 {
        if !self.exact && matches!(self.ids, LocIds::Scoped(_)) {
            0
        } else {
            self.strand
        }
    }
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
    feat_ids: Vec<FeatObjId>,
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
    feat_id_xrefs: Vec<FeatObjId>,
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
    gene_by_feat_id: HashMap<FeatObjId, Vec<usize>>,
    /// The number of Bioseqs of the TSE (keys of unresolved Seq-ids in `seq_key` follow them).
    n_bioseqs: usize,
    /// Per interned Seq-id: the Bioseq the scope finds for it (`resolve_ids`).
    bio_of: Vec<Option<u32>>,
    /// Per interned Seq-id: the Bioseqs whose features include a location on it, as a range of
    /// `attach_bios` (`attach_start[r]..attach_start[r + 1]`).
    attach_start: Vec<u32>,
    attach_bios: Vec<u32>,
    /// Per interned Seq-id: the first Seq-id with its handle (its key when the scope finds none).
    unresolved_key: Vec<u32>,
    /// The interned Seq-ids for which NCBI's scope may find one of several Bioseqs, or none, by
    /// Seq-id, with what it may find (`resolve_ids`).
    ties: Vec<Tie>,
    /// The pick of every tie that a title is made with, and the ties it read (`Scope::bio`).
    tie_state: RefCell<TieState>,
    /// Whether the TSE has Seq-ids that NCBI's object manager holds handles for and LOSAT does
    /// not intern: locations and products of other features and other `Seq-annot`s, products of
    /// genes, proteins and sources, and segments (`resolve_ids`).
    unindexed_ids: bool,
}

/// What NCBI's scope finds for a Seq-id that no Bioseq of the TSE has.
enum Found {
    None,
    One(u32),
    /// Any of these (`None`: no Bioseq), by an order or a state LOSAT does not know.
    Picks(Vec<Option<u32>>),
}

/// A Seq-id for which NCBI's scope may find any of `picks` (`None`: no Bioseq).
struct Tie {
    id: IdRef,
    picks: Box<[Option<u32>]>,
}

/// The pick of every tie (an index into its `picks`) and whether the title being made read it.
#[derive(Default)]
struct TieState {
    pick: Vec<usize>,
    read: Vec<bool>,
}

/// Whether the `GetStrand` with the scope of a location waits for the picks of its ties
/// (`Scope::apply_strand`), and whether it is then reversed (`reverse_loc`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum Deferred {
    #[default]
    No,
    Strand,
    Reversed,
}

/// A Bioseq Seq-id that NCBI's scope finds for another Seq-id (`GetMatchingHandles`): its Bioseq,
/// its tree and its version (`IsBetterVersion` ranks only within one tree).
#[derive(Clone, Copy)]
struct Candidate {
    bio: u32,
    tree: u8,
    version: Option<i32>,
}

/// Bioseq Seq-ids by the upper-case hash of a text field (accession or name), for the lookups of
/// `resolve_ids`.
struct FieldIndex<'a> {
    entries: Vec<(u64, u32, &'a ParsedSeqId)>,
}

impl<'a> FieldIndex<'a> {
    fn new(bioseq_ids: &[&'a [ParsedSeqId]], field: fn(&TextseqId) -> Option<&[u8]>) -> Self {
        let mut entries = Vec::new();
        for (bio, ids) in bioseq_ids.iter().enumerate() {
            let bio = u32::try_from(bio).expect("fewer than 2^32 Bioseqs in one TSE");
            for id in ids.iter() {
                if let Some(text) = id.textseq().and_then(field) {
                    entries.push((fold_hash(text), bio, id));
                }
            }
        }
        entries.sort_unstable_by_key(|e| (e.0, e.1));
        Self { entries }
    }

    /// The Bioseq Seq-ids whose field is `text` in any letter case.
    fn find<'s>(
        &'s self,
        text: &'s [u8],
        field: fn(&TextseqId) -> Option<&[u8]>,
    ) -> impl Iterator<Item = (u32, &'a ParsedSeqId)> + 's {
        let hash = fold_hash(text);
        let start = self.entries.partition_point(|e| e.0 < hash);
        self.entries[start..]
            .iter()
            .take_while(move |e| e.0 == hash)
            .filter(move |e| {
                e.2.textseq()
                    .and_then(field)
                    .is_some_and(|t| t.eq_ignore_ascii_case(text))
            })
            .map(|e| (e.1, e.2))
    }
}

fn text_accession(text: &TextseqId) -> Option<&[u8]> {
    text.accession.as_deref()
}

fn text_name(text: &TextseqId) -> Option<&[u8]> {
    text.name.as_deref()
}

/// A hash of `bytes` in upper case.
fn fold_hash(bytes: &[u8]) -> u64 {
    let mut hasher = FxHasher::default();
    for &c in bytes {
        hasher.write_u8(c.to_ascii_uppercase());
    }
    hasher.finish()
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
            if ext_c == s.seq_ext_seg || ext_c == s.seq_ext_delta {
                self.unindexed_ids = true;
            }
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
                self.unindexed_ids = true;
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
            self.unindexed_ids = true;
            return;
        };
        if fchoice != s.feat_prot
            && fchoice != s.feat_gene
            && fchoice != s.feat_cdregion
            && fchoice != s.feat_biosrc
        {
            self.unindexed_ids = true;
            return;
        }
        if fchoice != s.feat_cdregion && feat.member(s.seq_feat_product).is_some() {
            self.unindexed_ids = true;
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
            // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:1234-1241
            let mut feat_ids = Vec::new();
            if let Some(fid) = feat.member(s.seq_feat_id).and_then(local_feat_id) {
                feat_ids.push(fid);
            }
            if let Some(ids) = feat.member(s.seq_feat_ids).and_then(as_list) {
                for item in ids {
                    if let Some(fid) = local_feat_id(item) {
                        feat_ids.push(fid);
                    }
                }
            }
            let gene = extract_gene_ref(fval);
            let qual_comment_fallback = extract_qual_comment_fallback(feat);
            self.genes.push(GeneFeature {
                seq_order,
                feat_ids,
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
                    if let Some(fid) = xr.member(s.xref_id).and_then(local_feat_id) {
                        feat_id_xrefs.push(fid);
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

    /// Finds for every interned Seq-id the Bioseq that NCBI's scope finds (`x_GetSeqMatch`: the
    /// Bioseq with its handle, else the best of `GetMatchingHandles` with weak matching) and the
    /// Bioseqs whose features include a location on it (`GetTSESetWithAnnots`: the reverse
    /// matches of each Bioseq Seq-id in its own tree).
    ///
    /// The scope finds for a text Seq-id X that no Bioseq has, in every text tree: the Bioseq
    /// Seq-ids Y with the accession of X in any letter case and, when X has a version, that
    /// version; when X has a name, the Y with that name in any letter case (not when both have an
    /// accession) and, when X has a release, that release (a `swissprot` Y without one also for
    /// `reviewed` and `unreviewed`); when X is a `pir` or `prf` accession alone, the Y that are
    /// that text as a name alone, in a tree with no handle of an accession that is not packed and
    /// starts with the prefix of X when X is packed. A PDB X finds the Y of its molecule and
    /// chain, those with its release date when X has one. Of several Bioseqs, `IsBetter` keeps a
    /// Y with a version over one without and a higher version, within one tree; the order of the
    /// handles in memory decides the rest, so each of those Bioseqs is a pick of a tie (`ties`), as
    /// are both outcomes of a packed `pir` or `prf` X where Seq-ids LOSAT does not intern could
    /// hold such an accession (`unindexed_ids`). A `local` Seq-id finds the integer of its number
    /// string and the reverse, a `general` one the same within its database.
    ///
    /// A feature location on a text Seq-id X is under a Bioseq Seq-id Y of the tree of X with its
    /// accession in any letter case when X has no version or the version of Y; on a PDB Seq-id
    /// without a release date, under a Y of its molecule and chain with one; on another Seq-id,
    /// under the Bioseq with its handle.
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:865-882
    /// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_mapper.cpp:186-211,301-311
    /// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:1161-1200,1245-1295
    /// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:1355-1394,1438-1510
    /// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:1564-1625,1845-1890,2034-2058
    /// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:2696-2760,3467-3481,3552-3615
    /// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.hpp:440-442
    /// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:1019-1056
    fn resolve_ids(&mut self, bioseq_ids: &[&[ParsedSeqId]]) {
        self.n_bioseqs = bioseq_ids.len();
        let n = self.interner.ids.len();
        if n == 0 {
            return;
        }
        let accessions = FieldIndex::new(bioseq_ids, text_accession);
        let mut names: Option<FieldIndex<'_>> = None;
        let mut pdbs: Option<PdbIndex> = None;
        let mut bio_of = Vec::with_capacity(n);
        let mut attach_start = Vec::with_capacity(n + 1);
        let mut attach_bios = Vec::with_capacity(n);
        let mut first_with_key: FxHashMap<&HandleKey, u32> = FxHashMap::default();
        let mut unresolved_key = Vec::with_capacity(n);
        let mut ties = Vec::new();
        let mut cands: Vec<Candidate> = Vec::new();
        let mut name_cands: Vec<Candidate> = Vec::new();
        let bio_with_id = |id: &ParsedSeqId| match mask_handle_key(id) {
            Ok(Some(key)) => self
                .id_to_bioseq
                .get(&HandleKey::Bytes(key))
                .map(|&b| b as u32),
            _ => None,
        };
        for (r, id) in self.interner.ids.iter().enumerate() {
            let r32 = r as u32;
            let key = self.interner.keys[r].as_ref();
            let own = key
                .and_then(|k| self.id_to_bioseq.get(k))
                .map(|&b| b as u32);
            unresolved_key.push(match key {
                Some(k) => *first_with_key.entry(k).or_insert(r32),
                None => r32,
            });
            let start = attach_bios.len();
            attach_start.push(start as u32);
            attach_bios.extend(own);
            cands.clear();
            let found = match id {
                ParsedSeqId::Local(tag) => Found::of(own.or_else(|| {
                    let other = match tag {
                        ObjectId::Id(number) => ObjectId::Str(number.to_string().into_bytes()),
                        ObjectId::Str(_) => {
                            ObjectId::Id(object_id_number(tag).and_then(|v| i32::try_from(v).ok())?)
                        }
                    };
                    bio_with_id(&ParsedSeqId::Local(other))
                })),
                ParsedSeqId::General(dbtag) => Found::of(own.or_else(|| {
                    let other = match &dbtag.tag {
                        ObjectId::Id(number) if *number >= 0 => {
                            ObjectId::Str(number.to_string().into_bytes())
                        }
                        ObjectId::Id(_) => return None,
                        ObjectId::Str(text) => ObjectId::Id(canonical_non_negative(text)?),
                    };
                    bio_with_id(&ParsedSeqId::General(Dbtag {
                        db: dbtag.db.clone(),
                        tag: other,
                    }))
                })),
                ParsedSeqId::Pdb(pdb) => {
                    let pdbs = pdbs.get_or_insert_with(|| PdbIndex::new(bioseq_ids));
                    let stored = pdb_stored(pdb);
                    if stored.rel.is_none() {
                        for &bio in pdbs.dated.get(&pdb_str_key(&stored)).into_iter().flatten() {
                            if !attach_bios[start..].contains(&bio) {
                                attach_bios.push(bio);
                            }
                        }
                    }
                    match own {
                        Some(bio) => Found::One(bio),
                        None => {
                            for (bio, rel) in
                                pdbs.stored.get(&pdb_str_key(pdb)).into_iter().flatten()
                            {
                                if pdb.rel.is_none() || (rel.is_some() && *rel == pdb.rel) {
                                    cands.push(Candidate {
                                        bio: *bio,
                                        tree: 0,
                                        version: None,
                                    });
                                }
                            }
                            pick_candidate(&cands)
                        }
                    }
                }
                _ => match id.textseq() {
                    None => Found::of(own),
                    Some(text) => {
                        let tree = handle_tree(id);
                        if let Some(acc) = text.accession.as_deref() {
                            for (bio, y) in accessions.find(acc, text_accession) {
                                let y_text = y.textseq().expect("a text Seq-id");
                                if text.version.is_some() && text.version != y_text.version {
                                    continue;
                                }
                                let y_tree = handle_tree(y);
                                if y_tree == tree && !attach_bios[start..].contains(&bio) {
                                    attach_bios.push(bio);
                                }
                                cands.push(Candidate {
                                    bio,
                                    tree: y_tree,
                                    version: y_text.version,
                                });
                            }
                        }
                        if let Some(bio) = own {
                            Found::One(bio)
                        } else {
                            if let Some(name) = text.name.as_deref() {
                                let names = names
                                    .get_or_insert_with(|| FieldIndex::new(bioseq_ids, text_name));
                                for (bio, y) in names.find(name, text_name) {
                                    let y_text = y.textseq().expect("a text Seq-id");
                                    if text.accession.is_some() && y_text.accession.is_some() {
                                        continue;
                                    }
                                    if let Some(rel) = &text.release {
                                        let default_swissprot = y_text.release.is_none()
                                            && matches!(y, ParsedSeqId::Swissprot(_))
                                            && matches!(
                                                rel.as_slice(),
                                                b"reviewed" | b"unreviewed"
                                            );
                                        if !default_swissprot
                                            && y_text.release.as_ref() != Some(rel)
                                        {
                                            continue;
                                        }
                                    }
                                    cands.push(Candidate {
                                        bio,
                                        tree: handle_tree(y),
                                        version: y_text.version,
                                    });
                                }
                            }
                            name_cands.clear();
                            if let (
                                ParsedSeqId::Pir(_) | ParsedSeqId::Prf(_),
                                Some(acc),
                                None,
                                None,
                                None,
                            ) = (id, &text.accession, text.version, &text.name, &text.release)
                            {
                                let names = names
                                    .get_or_insert_with(|| FieldIndex::new(bioseq_ids, text_name));
                                for (bio, y) in names.find(acc, text_name) {
                                    let y_text = y.textseq().expect("a text Seq-id");
                                    if y_text.accession.is_none()
                                        && y_text.version.is_none()
                                        && y_text.release.is_none()
                                    {
                                        name_cands.push(Candidate {
                                            bio,
                                            tree: handle_tree(y),
                                            version: None,
                                        });
                                    }
                                }
                            }
                            match textseq_pack_prefix(text) {
                                Some(prefix) if !name_cands.is_empty() => {
                                    let unsure = self.accession_name_match(
                                        bioseq_ids,
                                        prefix,
                                        &mut name_cands,
                                    );
                                    cands.append(&mut name_cands);
                                    pick_candidate_unsure(&cands, &unsure)
                                }
                                _ => {
                                    cands.append(&mut name_cands);
                                    pick_candidate(&cands)
                                }
                            }
                        }
                    }
                },
            };
            match found {
                Found::None => bio_of.push(None),
                Found::One(bio) => bio_of.push(Some(bio)),
                Found::Picks(picks) => {
                    bio_of.push(picks[0]);
                    ties.push(Tie {
                        id: r32,
                        picks: picks.into_boxed_slice(),
                    });
                }
            }
        }
        attach_start.push(attach_bios.len() as u32);
        self.bio_of = bio_of;
        self.attach_start = attach_start;
        self.attach_bios = attach_bios;
        self.unresolved_key = unresolved_key;
        let n_ties = ties.len();
        self.ties = ties;
        self.tie_state = RefCell::new(TieState {
            pick: vec![0; n_ties],
            read: vec![false; n_ties],
        });
    }

    /// Keeps of the name matches `name_cands` of a packed `pir` or `prf` accession with the
    /// prefix `prefix` those that NCBI's `FindMatch` takes: in a tree where no handle has an
    /// accession that is not packed and starts with `prefix` in any letter case (the first such
    /// accession at or after `prefix` sends the search to the accessions). Moves out and gives
    /// the name matches that stay unsure: where Seq-ids that LOSAT does not intern may hold such
    /// an accession (`unindexed_ids`), and all of them where segments are, whose handles NCBI
    /// makes when it reads the segments.
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:1161-1200,1438-1490
    fn accession_name_match(
        &self,
        bioseq_ids: &[&[ParsedSeqId]],
        prefix: &[u8],
        name_cands: &mut Vec<Candidate>,
    ) -> Vec<Candidate> {
        if !self.bioseq_seg_locs.is_empty() {
            return std::mem::take(name_cands);
        }
        let blocked = |tree: u8| {
            bioseq_ids
                .iter()
                .flat_map(|ids| ids.iter())
                .chain(self.interner.ids.iter())
                .any(|y| {
                    handle_tree(y) == tree
                        && y.textseq().is_some_and(|t| {
                            !textseq_packs(t)
                                && t.accession.as_deref().is_some_and(|a| {
                                    a.len() >= prefix.len()
                                        && a[..prefix.len()].eq_ignore_ascii_case(prefix)
                                })
                        })
                })
        };
        let mut open_trees: Vec<u8> = name_cands.iter().map(|c| c.tree).collect();
        open_trees.sort_unstable();
        open_trees.dedup();
        open_trees.retain(|&tree| !blocked(tree));
        name_cands.retain(|c| open_trees.contains(&c.tree));
        if self.unindexed_ids {
            std::mem::take(name_cands)
        } else {
            Vec::new()
        }
    }

    /// The Bioseqs whose features include a location on the interned Seq-id `r`.
    fn attached(&self, r: IdRef) -> &[u32] {
        let r = r as usize;
        &self.attach_bios[self.attach_start[r] as usize..self.attach_start[r + 1] as usize]
    }

    /// The Bioseq the scope finds for the interned Seq-id `r` (`Scope::bio`).
    fn bio_of(&self, r: IdRef) -> Option<usize> {
        self.scope().bio(r).map(|b| b as usize)
    }

    /// `GetStrand` of `loc` with the scope (`Scope::strand`).
    fn loc_strand(&self, loc: &CompactLoc) -> i32 {
        self.scope().strand(loc)
    }

    /// Makes a title with `make` at the first pick of every tie and, when it read ties, again at
    /// every combination of the picks of the ties it read; gives the title when all are the
    /// same (NCBI makes one of them), `None` when they differ.
    fn same_at_every_pick<T: PartialEq>(&self, mut make: impl FnMut() -> T) -> Option<T> {
        if self.ties.is_empty() {
            return Some(make());
        }
        let set_picks = |read: &[usize], digits: &[usize]| {
            let mut state = self.tie_state.borrow_mut();
            state.pick.iter_mut().for_each(|p| *p = 0);
            state.read.iter_mut().for_each(|r| *r = false);
            for (&t, &d) in read.iter().zip(digits) {
                state.pick[t] = d;
            }
        };
        let newly_read = |read: &[usize]| -> Vec<usize> {
            let state = self.tie_state.borrow();
            (0..state.read.len())
                .filter(|t| state.read[*t] && !read.contains(t))
                .collect()
        };
        set_picks(&[], &[]);
        let first = make();
        let mut read = newly_read(&[]);
        let mut same = true;
        'restart: while !read.is_empty() {
            let mut digits = vec![0usize; read.len()];
            loop {
                let mut k = 0;
                loop {
                    if k == digits.len() {
                        break 'restart;
                    }
                    digits[k] += 1;
                    if digits[k] < self.ties[read[k]].picks.len() {
                        break;
                    }
                    digits[k] = 0;
                    k += 1;
                }
                set_picks(&read, &digits);
                if make() != first {
                    same = false;
                    break 'restart;
                }
                let new = newly_read(&read);
                if !new.is_empty() {
                    read.extend(new);
                    continue 'restart;
                }
            }
        }
        set_picks(&[], &[]);
        same.then_some(first)
    }

    /// The key of the sequence of the interned Seq-id `r` in `TestForOverlap64`'s ranges
    /// (`s_GetSynHandle`): its Bioseq, or after the Bioseqs its handle where the scope finds none.
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:985-1003
    fn seq_key(&self, r: IdRef) -> usize {
        match self.bio_of(r) {
            Some(bio) => bio,
            None => self.n_bioseqs + self.unresolved_key[r as usize] as usize,
        }
    }

    /// `IsSameBioseq` and `GetId` with the TSE scope.
    fn scope(&self) -> Scope<'_> {
        Scope {
            interner: &self.interner,
            bio_of: &self.bio_of,
            ties: &self.ties,
            tie_state: &self.tie_state,
        }
    }

    /// Finalizes inverted indices once all Bioseq IDs are mapped in `id_to_bioseq`; `bioseq_ids`
    /// are the Seq-ids of each Bioseq of the TSE in order.
    pub(crate) fn build_indices(&mut self, bioseq_ids: &[&[ParsedSeqId]]) {
        self.resolve_ids(bioseq_ids);
        if self.prots.is_empty()
            && self.genes.is_empty()
            && self.cds.is_empty()
            && self.biosrcs.is_empty()
        {
            return;
        }
        let scope = Scope {
            interner: &self.interner,
            bio_of: &self.bio_of,
            ties: &self.ties,
            tie_state: &self.tie_state,
        };
        let locs = self.genes.iter_mut().map(|g| &mut g.loc);
        let locs = locs.chain(self.cds.iter_mut().map(|c| &mut c.loc));
        for loc in locs.chain(self.biosrcs.iter_mut().map(|b| &mut b.loc)) {
            scope.apply_strand(loc);
        }
        for (idx, feat) in self.prots.iter().enumerate() {
            index_loc_bioseqs(
                &feat.loc,
                &self.attach_start,
                &self.attach_bios,
                &mut self.prots_by_bioseq,
                idx,
            );
        }
        let mut raw_genes: HashMap<usize, Vec<usize>> = HashMap::new();
        for (idx, feat) in self.genes.iter().enumerate() {
            for fid in &feat.feat_ids {
                self.gene_by_feat_id
                    .entry(fid.clone())
                    .or_default()
                    .push(idx);
            }
            index_loc_bioseqs(
                &feat.loc,
                &self.attach_start,
                &self.attach_bios,
                &mut raw_genes,
                idx,
            );
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
            index_loc_bioseqs(
                &feat.loc,
                &self.attach_start,
                &self.attach_bios,
                &mut self.cds_by_bioseq,
                idx,
            );
            if let Some(prod) = &feat.product {
                index_loc_bioseqs(
                    prod,
                    &self.attach_start,
                    &self.attach_bios,
                    &mut self.cds_by_product,
                    idx,
                );
            }
        }
        let mut raw_biosrcs: HashMap<usize, Vec<usize>> = HashMap::new();
        for (idx, feat) in self.biosrcs.iter().enumerate() {
            index_loc_bioseqs(
                &feat.loc,
                &self.attach_start,
                &self.attach_bios,
                &mut raw_biosrcs,
                idx,
            );
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
    /// segment Bioseqs of a segmented Bioseq, found through the scope).
    fn target_bioseqs(&self, bio_idx: usize) -> Vec<usize> {
        let mut out = vec![bio_idx];
        if let Some(seg_loc) = self.bioseq_seg_locs.get(&bio_idx) {
            for p in &seg_loc.pieces {
                if let Some(seg_bio) = p.id.and_then(|r| self.bio_of(r)) {
                    if !out.contains(&seg_bio) {
                        out.push(seg_bio);
                    }
                }
            }
        }
        out
    }

    /// `sequence::GetLength`: a `whole` piece has the length of the Bioseq the scope finds.
    fn loc_length(&self, loc: &CompactLoc) -> u64 {
        let mut total = 0u64;
        for p in &loc.pieces {
            if p.is_whole {
                if let Some(len) =
                    p.id.and_then(|r| self.bio_of(r))
                        .and_then(|bi| self.bioseq_lengths.get(bi).copied().flatten())
                {
                    total = total.saturating_add(len);
                }
            } else {
                total = total.saturating_add(p.to_open.saturating_sub(p.from));
            }
        }
        total
    }

    /// The total range on Bioseq `bio_idx` of the pieces of `loc` whose Seq-id `on` takes.
    fn loc_total_range_by(
        &self,
        loc: &CompactLoc,
        bio_idx: usize,
        on: impl Fn(IdRef) -> bool,
    ) -> Option<(u64, u64)> {
        let mut min_from = u64::MAX;
        let mut max_to = 0u64;
        let mut any = false;
        for p in &loc.pieces {
            let Some(r) = p.id else { continue };
            if !on(r) {
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

    /// The total range of a feature location `loc` among the features of Bioseq `bio_idx`
    /// (`CFeat_CI`): its pieces under that Bioseq's Seq-ids.
    fn loc_total_range_on(&self, loc: &CompactLoc, bio_idx: usize) -> Option<(u64, u64)> {
        self.loc_total_range_by(loc, bio_idx, |r| {
            self.attached(r).contains(&(bio_idx as u32))
        })
    }

    /// The total range of `loc` on Bioseq `bio_idx` as the scope finds its pieces' Seq-ids.
    fn loc_total_range_found(&self, loc: &CompactLoc, bio_idx: usize) -> Option<(u64, u64)> {
        self.loc_total_range_by(loc, bio_idx, |r| self.bio_of(r) == Some(bio_idx))
    }

    /// The ranges of `loc` by sequence (`seq_key`) and orientation for `TestForOverlap64`; all on
    /// one sequence when `single` (`s_TestForOverlapEx` keys both locations with the first
    /// location's Seq-id when both have one).
    fn loc_ranges_by_rev(&self, loc: &CompactLoc, single: bool) -> Vec<(usize, bool, u64, u64)> {
        let mut out = Vec::with_capacity(loc.pieces.len());
        for p in &loc.pieces {
            let Some(r) = p.id else { continue };
            let (from, to_open) = if p.is_whole {
                let len = self
                    .bio_of(r)
                    .and_then(|bi| self.bioseq_lengths.get(bi).copied().flatten())
                    .unwrap_or(0);
                (0, len)
            } else {
                (p.from, p.to_open)
            };
            let key = if single { 0 } else { self.seq_key(r) };
            out.push((key, is_reverse_strand(p.strand), from, to_open));
        }
        out
    }

    /// The first range of `loc` on a Bioseq the scope finds, and that Bioseq.
    fn first_found_bio(&self, ranges: &[(usize, bool, u64, u64)]) -> Option<usize> {
        ranges.iter().map(|r| r.0).find(|&key| key < self.n_bioseqs)
    }

    /// `s_SeqLocToTotalRangesInfoMapByStrand`: collapses consecutive same-Bioseq, same-orientation,
    /// properly-ordered pieces into total ranges per `(bio_idx, is_reverse)`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:1385-1476
    fn total_ranges_by_rev(&self, loc: &CompactLoc, single: bool) -> Vec<(usize, bool, u64, u64)> {
        let pieces = self.loc_ranges_by_rev(loc, single);
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

    /// The start of `s_TestForOverlapEx`: `None` (no overlap) when both locations have one Seq-id
    /// (`CSeq_loc::GetId`) that `IsSameBioseq` with the scope tells apart, or when they are on one
    /// sequence and their strands (`GetStrand` without a scope) cannot overlap; else whether both
    /// have one Seq-id (one sequence).
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:1831-1891
    fn overlap_start(&self, outer: &CompactLoc, inner: &CompactLoc) -> Option<bool> {
        let single = outer.exact && inner.exact;
        if single {
            if let (Some(a), Some(b)) = (outer.first_id, inner.first_id) {
                if !self.scope().same(a, b) {
                    return None;
                }
            }
        }
        let (s1, s2) = (outer.strand_free(), inner.strand_free());
        if single && !test_strands(s1, s2) && s1 != 255 && s2 != 255 {
            return None;
        }
        Some(single)
    }

    /// `TestForOverlap64` under `eOverlap_Subset` (also used by `eOverlap_SubsetRev` after
    /// swapping arguments).
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:1535-1589,1874-1909
    fn test_for_overlap_subset(&self, outer: &CompactLoc, inner: &CompactLoc) -> Option<u64> {
        let single = self.overlap_start(outer, inner)?;
        let o_ranges = self.loc_ranges_by_rev(outer, single);
        let i_ranges = self.loc_ranges_by_rev(inner, single);
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
        let single = self.overlap_start(outer, inner)?;
        let o_tot = self.total_ranges_by_rev(outer, single);
        let i_tot = self.total_ranges_by_rev(inner, single);
        if !ranges_contain_all(&o_tot, &i_tot) {
            return None;
        }
        // `Compare` with `fCompareOverlapping` requires `eContains` or `eSame`: every individual
        // interval of `inner` must be contained in some interval of `outer` on the same Bioseq.
        let o_int = self.loc_ranges_by_rev(outer, false);
        let i_int = self.loc_ranges_by_rev(inner, false);
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
            .loc_ranges_by_rev(loc, false)
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
            let Some(bi) = p.id.and_then(|r| self.bio_of(r)) else {
                continue;
            };
            if let (Some((cf, ct)), Some((gf, gt))) = (
                self.loc_total_range_found(cds_loc, bi),
                self.loc_total_range_on(cand_loc, bi),
            ) {
                if gf < ct && cf < gt {
                    return true;
                }
            }
        }
        false
    }

    /// `s_GetCdregionLabel` + `feature::GetLabel` (`fFGL_Content`) for a `Cdregion` feature: the
    /// product's protein name comes from the Bioseq the scope finds for the product's `GetId`
    /// with the scope.
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/feature.cpp:139-247,710-776
    fn cdregion_content_label(&self, cds: &CdsFeature, warn: Warn<'_>) -> Vec<u8> {
        let mut tlabel = Vec::new();
        if let Some(prot_lbl) = &cds.xref_prot_label {
            tlabel.extend_from_slice(prot_lbl);
        } else {
            let mut found_from_product = false;
            if let Some(id) = cds.product.as_ref().and_then(|p| self.scope().loc_id(p)) {
                if let Some(prod_bio) = self.bio_of(id) {
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
                } else {
                    let mut msg = b"Error: [makeblastdb] cannot find sequence: ".to_vec();
                    label::fasta(self.interner.id(id), &mut msg);
                    warn(&msg);
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
                    let strand_str: &[u8] = match self.loc_strand(&cds.loc) {
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
        for fid in &cds.feat_id_xrefs {
            if let Some(gis) = self.gene_by_feat_id.get(fid) {
                for &gi in gis {
                    let g = &self.genes[gi];
                    if self.total_range_and_strand_intersects(&cds.loc, &g.loc) {
                        return Some(g);
                    }
                }
            }
        }
        let cds_ranges = self.loc_ranges_by_rev(&cds.loc, false);
        let cds_bio = self.first_found_bio(&cds_ranges)?;
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
                } else if cds.trans_splicing
                    && cds.loc.strand_free() == 255
                    && self.loc_strand(&g.loc) != 255
                {
                    let g_strand = self.loc_strand(&g.loc);
                    let pstrand = if g_strand == 0 { 1 } else { g_strand };
                    let forced = force_loc_strand(&cds.loc, pstrand);
                    let Some(d) = self.test_for_overlap_subset(&g.loc, &forced) else {
                        return;
                    };
                    (1i8, d)
                } else {
                    return;
                };
            let (g_from, g_to) = self.loc_total_range_on(&g.loc, cds_bio).unwrap_or((0, 0));
            let tie = gene_annot_tie_key(g, self.loc_strand(&g.loc));
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
        let cds_ranges = self.loc_ranges_by_rev(&cds.loc, false);
        let cds_bio = self.first_found_bio(&cds_ranges)?;
        let src_index = self.biosrcs_by_bioseq.get(&cds_bio)?;
        let (cds_from, cds_to) = ranges_span_on(&cds_ranges, cds_bio)?;

        let find_for = |target: &CompactLoc| -> Option<&BiosrcFeature> {
            let mut best: Option<(u64, u64, std::cmp::Reverse<u64>, BiosrcAnnotTieKey, usize)> =
                None;
            src_index.for_each_covering(cds_from, cds_to, |si| {
                let src = &self.biosrcs[si];
                if !cfeat_ci_intervals_strand_match(self.loc_strand(target), &src.loc) {
                    return;
                }
                let Some(diff) = self.test_for_overlap_subset(&src.loc, target) else {
                    return;
                };
                let (s_from, s_to) = self.loc_total_range_on(&src.loc, cds_bio).unwrap_or((0, 0));
                let tie = biosrc_annot_tie_key(src, self.loc_strand(&src.loc));
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
        let cds_ranges = self.loc_ranges_by_rev(&cds.loc, false);
        let cds_bio = self.first_found_bio(&cds_ranges)?;
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
                gene_annot_tie_key(g, self.loc_strand(&g.loc)),
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

fn gene_annot_tie_key(g: &GeneFeature, strand: i32) -> GeneAnnotTieKey {
    let is_rev = is_reverse_strand(strand);
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

fn biosrc_annot_tie_key(s: &BiosrcFeature, strand: i32) -> BiosrcAnnotTieKey {
    let is_rev = is_reverse_strand(strand);
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
fn cfeat_ci_intervals_strand_match(q_strand: i32, cand: &CompactLoc) -> bool {
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
    let pieces: Vec<LocPiece> = loc
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
    CompactLoc {
        pieces,
        strand,
        is_complex: loc.is_complex,
        exact: loc.exact,
        first_id: loc.first_id,
        ids: loc.ids.clone(),
        deferred: match loc.deferred {
            Deferred::No => Deferred::No,
            Deferred::Strand => Deferred::Reversed,
            Deferred::Reversed => Deferred::Strand,
        },
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
        exact: loc.exact,
        first_id: loc.first_id,
        ids: loc.ids.clone(),
        deferred: Deferred::No,
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

/// Adds feature `feat_idx` to `map` under every Bioseq whose features include a piece of `loc`
/// (`TseFeatureIndex::attached`; `empty` sub-locations add none).
fn index_loc_bioseqs(
    loc: &CompactLoc,
    attach_start: &[u32],
    attach_bios: &[u32],
    map: &mut HashMap<usize, Vec<usize>>,
    feat_idx: usize,
) {
    for p in &loc.pieces {
        let Some(r) = p.id else { continue };
        let r = r as usize;
        for &bio in &attach_bios[attach_start[r] as usize..attach_start[r + 1] as usize] {
            let list = map.entry(bio as usize).or_default();
            if list.last() != Some(&feat_idx) {
                list.push(feat_idx);
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

/// `IsSameBioseq` and `GetId` with the TSE scope (`TseFeatureIndex::resolve_ids`).
struct Scope<'a> {
    interner: &'a SeqIdInterner,
    bio_of: &'a [Option<u32>],
    ties: &'a [Tie],
    tie_state: &'a RefCell<TieState>,
}

impl Scope<'_> {
    /// The Bioseq the scope finds for the interned Seq-id `r`; for a tie, the pick that the title
    /// is made with, noting that the title read it.
    fn bio(&self, r: IdRef) -> Option<u32> {
        match self.tie(r) {
            None => self.bio_of[r as usize],
            Some(t) => {
                let mut state = self.tie_state.borrow_mut();
                state.read[t] = true;
                self.ties[t].picks[state.pick[t]]
            }
        }
    }

    fn tie(&self, r: IdRef) -> Option<usize> {
        if self.ties.is_empty() {
            return None;
        }
        self.ties.binary_search_by_key(&r, |t| t.id).ok()
    }

    /// `IsSameBioseq` with a scope: `CSeq_id::Match`, one handle, or one Bioseq that the scope
    /// finds for both.
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:300-321
    /// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:1880-1893
    fn same(&self, a: IdRef, b: IdRef) -> bool {
        self.interner.same_free(a, b) || self.bio(a).is_some_and(|bio| self.bio(b) == Some(bio))
    }

    /// `GetId` with the scope: the first Seq-id of `loc` when `s_GetId` finds every Seq-id on its
    /// Bioseq.
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:325-370
    fn loc_id(&self, loc: &CompactLoc) -> Option<IdRef> {
        let first = loc.first_id?;
        match &loc.ids {
            LocIds::None => None,
            LocIds::Free => Some(first),
            LocIds::Scoped(scoped) => scoped
                .ids
                .iter()
                .all(|&x| self.same(first, x))
                .then_some(first),
        }
    }

    /// `GetStrand` with the scope for a location that is on one Bioseq only through the scope:
    /// `s_GetStrand` when `IsOneBioseq` holds. A location on a tie waits for the picks
    /// (`TseFeatureIndex::loc_strand`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:385-389,481-513
    fn apply_strand(&self, loc: &mut CompactLoc) {
        if let LocIds::Scoped(scoped) = &loc.ids {
            if loc.exact {
                return;
            }
            if loc.first_id.is_some_and(|f| self.tie(f).is_some())
                || scoped.ids.iter().any(|&x| self.tie(x).is_some())
            {
                loc.deferred = Deferred::Strand;
            } else if self.loc_id(loc).is_some() {
                loc.strand = scoped.strand;
            }
        }
    }

    /// `GetStrand` with the scope of `loc` (`apply_strand`, at the picks of its ties).
    fn strand(&self, loc: &CompactLoc) -> i32 {
        if loc.deferred == Deferred::No {
            return loc.strand;
        }
        let strand = match &loc.ids {
            LocIds::Scoped(scoped) if self.loc_id(loc).is_some() => scoped.strand,
            _ => 0,
        };
        if loc.deferred == Deferred::Reversed {
            reverse_strand(strand)
        } else {
            strand
        }
    }
}

/// The Bioseqs that `x_GetSeqMatch` may keep among the Bioseq Seq-ids `cands` that the scope
/// finds for a Seq-id, going through them in the order of their handles (`None` when there are
/// none): it takes each one unless the one it holds `IsBetter`, which holds only within one tree
/// for a version over none or a higher version. So within one tree the Bioseqs of the best
/// version remain; across trees, every one.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:865-882
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_mapper.cpp:301-311
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:99-103,1564-1572
fn pick_set(cands: &[Candidate], picks: &mut Vec<Option<u32>>) {
    let Some(first) = cands.first() else {
        if !picks.contains(&None) {
            picks.push(None);
        }
        return;
    };
    let across = cands.iter().any(|c| c.tree != first.tree);
    let best = cands.iter().map(|c| c.version).max().flatten();
    for c in cands {
        if (across || c.version == best) && !picks.contains(&Some(c.bio)) {
            picks.push(Some(c.bio));
        }
    }
}

/// What the scope finds among `cands` (`pick_set`).
fn pick_candidate(cands: &[Candidate]) -> Found {
    match cands {
        [] => Found::None,
        [one] => Found::One(one.bio),
        _ if cands.iter().all(|c| c.bio == cands[0].bio) => Found::One(cands[0].bio),
        _ => {
            let mut picks = Vec::new();
            pick_set(cands, &mut picks);
            Found::from_picks(picks)
        }
    }
}

/// What the scope finds among `cands` and any of the trees of `unsure` (whose name matches NCBI
/// may or may not take, `accession_name_match`).
fn pick_candidate_unsure(cands: &[Candidate], unsure: &[Candidate]) -> Found {
    if unsure.is_empty() {
        return pick_candidate(cands);
    }
    let mut trees: Vec<u8> = unsure.iter().map(|c| c.tree).collect();
    trees.sort_unstable();
    trees.dedup();
    let mut picks = Vec::new();
    let mut all: Vec<Candidate> = Vec::with_capacity(cands.len() + unsure.len());
    for mask in 0u64..(1u64 << trees.len()) {
        all.clear();
        all.extend_from_slice(cands);
        all.extend(unsure.iter().copied().filter(|c| {
            let at = trees
                .iter()
                .position(|&t| t == c.tree)
                .expect("a tree of unsure");
            (mask >> at) & 1 == 1
        }));
        pick_set(&all, &mut picks);
    }
    Found::from_picks(picks)
}

impl Found {
    fn of(bio: Option<u32>) -> Self {
        bio.map_or(Found::None, Found::One)
    }

    fn from_picks(picks: Vec<Option<u32>>) -> Self {
        match picks.as_slice() {
            [one] => Found::of(*one),
            _ => Found::Picks(picks),
        }
    }
}

/// `NStr::StringToNonNegativeInt` of `text` when `NStr::IntToString` gives `text` back.
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:2750-2758
fn canonical_non_negative(text: &[u8]) -> Option<i32> {
    if text.is_empty()
        || !text.iter().all(u8::is_ascii_digit)
        || (text.len() > 1 && text[0] == b'0')
    {
        return None;
    }
    std::str::from_utf8(text).ok()?.parse::<i32>().ok()
}

/// The Seq-id that NCBI's PDB tree keeps for the handle of `pdb` (`Normalize`): the normalized
/// one, without its release date, or `pdb` itself when its date does not normalize.
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:3408-3463,3502-3523
fn pdb_stored(pdb: &PdbSeqId) -> PdbSeqId {
    pdb_normalized(pdb).unwrap_or_else(|| pdb.clone())
}

/// `x_IdToStrKey`: the molecule, then `_` and the `chain-id` or else the `chain` as a byte.
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:3467-3481
fn pdb_str_key(pdb: &PdbSeqId) -> Vec<u8> {
    let mut key = pdb.mol.clone();
    if let Some(chain_id) = &pdb.chain_id {
        key.push(b'_');
        key.extend_from_slice(chain_id);
    } else if let Some(chain) = pdb.chain {
        key.push(b'_');
        key.push(chain as u8);
    }
    key
}

/// The Bioseq PDB Seq-ids of a TSE for `resolve_ids`: by the key of the Seq-id NCBI keeps
/// (`stored`, with that Seq-id's release date), and those with a release date by their own key
/// (`dated`, whose reverse matches add the Seq-ids kept without a date, `FindReverseMatch`).
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:3552-3615
struct PdbIndex {
    stored: FxHashMap<Vec<u8>, Vec<(u32, Option<Date>)>>,
    dated: FxHashMap<Vec<u8>, Vec<u32>>,
}

impl PdbIndex {
    fn new(bioseq_ids: &[&[ParsedSeqId]]) -> Self {
        let mut index = PdbIndex {
            stored: FxHashMap::default(),
            dated: FxHashMap::default(),
        };
        for (bio, ids) in bioseq_ids.iter().enumerate() {
            for id in ids.iter() {
                if let ParsedSeqId::Pdb(pdb) = id {
                    let stored = pdb_stored(pdb);
                    index
                        .stored
                        .entry(pdb_str_key(&stored))
                        .or_default()
                        .push((bio as u32, stored.rel));
                    if pdb.rel.is_some() {
                        let list = index.dated.entry(pdb_str_key(pdb)).or_default();
                        if !list.contains(&(bio as u32)) {
                            list.push(bio as u32);
                        }
                    }
                }
            }
        }
        index
    }
}

/// What `GetStrand` without a scope gives a sub-location, as `s_GetStrand` of a `mix` folds it.
enum NodeStrand {
    /// `null` and `empty` (skipped by the fold).
    Skip,
    /// The strand, whatever the Seq-ids (`int`, `pnt`, `packed-pnt`, `whole`, `equiv`, ...).
    Direct(i32),
    /// `s_GetStrand` of a `packed-int` or `bond`, given only when it is on one Bioseq.
    OneBioseq(i32),
}

/// The `s_GetStrand` fold of a `mix` (or of a `packed-int`): `other` once two strands disagree.
#[derive(Clone, Copy, Default)]
struct StrandFold {
    strand: i32,
    set: bool,
    other: bool,
    /// An `equiv`, whose `s_GetStrand` is unknown whatever its parts.
    sink: bool,
}

impl StrandFold {
    fn add(&mut self, istrand: i32) {
        if !self.sink
            && !self.other
            && !advance_strand_state(&mut self.strand, &mut self.set, istrand)
        {
            self.other = true;
        }
    }

    fn value(&self) -> i32 {
        if self.sink {
            0
        } else if self.other {
            255
        } else {
            self.strand
        }
    }
}

/// Extracts a `CompactLoc` from a `Seq-loc` `Value`: its pieces, its Seq-ids in `CSeq_loc_CI`
/// order (`empty` sub-locations included, `null` skipped; `CSeq_loc::GetId`, `s_GetId`) and
/// `sequence::GetStrand`, whose `s_GetStrand` gives a nested `packed-int`, `mix` or `bond` its
/// strand only when it is on one Bioseq without a scope (`GetStrand` with no scope); the top
/// level asks the scope in `build_indices` when only the scope can tell.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_loc_util.cpp:300-513
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_loc.cpp:3541-3630
fn extract_loc(root: &Value, interner: &mut SeqIdInterner) -> CompactLoc {
    let s = dschema();
    let mut out = CompactLoc::default();
    let Some((choice, inner)) = as_choice(root) else {
        return out;
    };
    out.is_complex = choice == s.seqloc_mix || choice == s.seqloc_packed_int;
    let mut ci = std::mem::take(&mut interner.ci);
    ci.clear();
    let node = if choice == s.seqloc_mix || choice == s.seqloc_equiv {
        NodeStrand::OneBioseq(extract_nested(root, interner, &mut out.pieces, &mut ci))
    } else {
        visit_flat(choice, inner, interner, &mut out.pieces, &mut ci)
    };
    let (first_id, exact, ids) = interner.loc_ids(&ci);
    interner.ci = ci;
    out.first_id = first_id;
    out.exact = exact;
    out.strand = match node {
        NodeStrand::Skip => 0,
        NodeStrand::Direct(strand) => strand,
        NodeStrand::OneBioseq(strand) => match &ids {
            LocIds::None => 0,
            LocIds::Free => strand,
            LocIds::Scoped(_) if exact => strand,
            LocIds::Scoped(_) => 0,
        },
    };
    out.ids = match ids {
        LocIds::Scoped(mut scoped) => {
            if let NodeStrand::OneBioseq(strand) = node {
                scoped.strand = strand;
            }
            LocIds::Scoped(scoped)
        }
        other => other,
    };
    out
}

/// Walks a `mix` or `equiv` location depth first in `CSeq_loc_CI` order, adding its pieces and
/// Seq-ids, and returns its `s_GetStrand`.
fn extract_nested(
    root: &Value,
    interner: &mut SeqIdInterner,
    pieces: &mut Vec<LocPiece>,
    ci: &mut Vec<IdRef>,
) -> i32 {
    enum Frame<'a> {
        Visit(&'a Value),
        /// The end of a `mix` or `equiv` whose Seq-ids start at `ci[start]`.
        Close {
            start: usize,
        },
    }
    let s = dschema();
    let mut frames = vec![Frame::Visit(root)];
    let mut folds: Vec<StrandFold> = Vec::new();
    let mut top = 0;
    while let Some(frame) = frames.pop() {
        let value = match frame {
            Frame::Visit(loc) => match as_choice(loc) {
                Some((c, inner)) if c == s.seqloc_mix || c == s.seqloc_equiv => {
                    frames.push(Frame::Close { start: ci.len() });
                    folds.push(StrandFold {
                        sink: c == s.seqloc_equiv,
                        ..StrandFold::default()
                    });
                    if let Some(items) = as_list(inner) {
                        frames.extend(items.iter().rev().map(Frame::Visit));
                    }
                    continue;
                }
                Some((c, inner)) => {
                    let start = ci.len();
                    match visit_flat(c, inner, interner, pieces, ci) {
                        NodeStrand::Skip => continue,
                        NodeStrand::Direct(strand) => strand,
                        NodeStrand::OneBioseq(strand) => {
                            if interner.one_free(&ci[start..]) {
                                strand
                            } else {
                                0
                            }
                        }
                    }
                }
                None => 0,
            },
            Frame::Close { start } => {
                let fold = folds.pop().expect("an open mix or equiv");
                if folds.is_empty() {
                    top = fold.value();
                    continue;
                }
                if fold.sink || !interner.one_free(&ci[start..]) {
                    0
                } else {
                    fold.value()
                }
            }
        };
        if let Some(parent) = folds.last_mut() {
            parent.add(value);
        }
    }
    top
}

/// Adds the pieces and Seq-ids of a location with no nested locations and returns what
/// `GetStrand` gives it.
fn visit_flat(
    choice: usize,
    inner: &Value,
    interner: &mut SeqIdInterner,
    pieces: &mut Vec<LocPiece>,
    ci: &mut Vec<IdRef>,
) -> NodeStrand {
    let s = dschema();
    if choice == s.seqloc_empty {
        if let Value::SeqId(sid) = inner {
            ci.push(interner.intern(sid));
        }
        NodeStrand::Skip
    } else if choice == s.seqloc_null {
        NodeStrand::Skip
    } else if choice == s.seqloc_whole {
        if let Value::SeqId(sid) = inner {
            let r = interner.intern(sid);
            ci.push(r);
            pieces.push(LocPiece {
                id: Some(r),
                from: 0,
                to_open: 0,
                is_whole: true,
                strand: 0,
            });
        }
        NodeStrand::Direct(3) // eNa_strand_both
    } else if choice == s.seqloc_int {
        NodeStrand::Direct(push_seq_int(inner, interner, pieces, ci))
    } else if choice == s.seqloc_pnt {
        NodeStrand::Direct(push_seq_pnt(inner, interner, pieces, ci))
    } else if choice == s.seqloc_packed_int {
        let mut fold = StrandFold::default();
        if let Some(ints) = as_list(inner) {
            for iv in ints {
                fold.add(push_seq_int(iv, interner, pieces, ci));
            }
        }
        NodeStrand::OneBioseq(fold.value())
    } else if choice == s.seqloc_packed_pnt {
        let strand = match inner.member(s.packed_pnt_strand) {
            Some(Value::Enum(v)) => *v,
            _ => 0,
        };
        let id = match inner.member(s.packed_pnt_id) {
            Some(Value::SeqId(sid)) => Some(interner.intern(sid)),
            _ => None,
        };
        let mut push = |pos: u64| {
            if let Some(r) = id {
                ci.push(r);
            }
            pieces.push(LocPiece {
                id,
                from: pos,
                to_open: pos.saturating_add(1),
                is_whole: false,
                strand,
            });
        };
        match inner.member(s.packed_pnt_points) {
            Some(Value::List(items)) => {
                for p in items {
                    if let Value::Int(pt) = p {
                        push(*pt as u64);
                    }
                }
            }
            Some(Value::Ints(items)) => items.iter().for_each(|&pt| push(pt as u64)),
            Some(Value::Int4s(items)) => items.iter().for_each(|&pt| push(pt as u64)),
            _ => {}
        }
        NodeStrand::Direct(strand)
    } else if choice == s.seqloc_bond {
        let mut a_strand = inner
            .member(s.seqbond_a)
            .map_or(0, |a| push_seq_pnt(a, interner, pieces, ci));
        let mut b_strand = inner
            .member(s.seqbond_b)
            .map_or(0, |b| push_seq_pnt(b, interner, pieces, ci));
        if a_strand == 0 && b_strand != 0 {
            a_strand = b_strand;
        } else if a_strand != 0 && b_strand == 0 {
            b_strand = a_strand;
        }
        NodeStrand::OneBioseq(if a_strand != b_strand { 255 } else { a_strand })
    } else {
        NodeStrand::Direct(0)
    }
}

/// Adds a `Seq-interval` piece and returns its strand.
fn push_seq_int(
    iv: &Value,
    interner: &mut SeqIdInterner,
    pieces: &mut Vec<LocPiece>,
    ci: &mut Vec<IdRef>,
) -> i32 {
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
    let id = match iv.member(s.seqint_id) {
        Some(Value::SeqId(sid)) => Some(interner.intern(sid)),
        _ => None,
    };
    ci.extend(id);
    pieces.push(LocPiece {
        id,
        from,
        to_open: to.saturating_add(1),
        is_whole: false,
        strand,
    });
    strand
}

/// Adds a `Seq-point` piece and returns its strand.
fn push_seq_pnt(
    pt: &Value,
    interner: &mut SeqIdInterner,
    pieces: &mut Vec<LocPiece>,
    ci: &mut Vec<IdRef>,
) -> i32 {
    let s = dschema();
    let pos = match pt.member(s.seqpnt_point) {
        Some(Value::Int(v)) => *v as u64,
        _ => 0,
    };
    let strand = match pt.member(s.seqpnt_strand) {
        Some(Value::Enum(v)) => *v,
        _ => 0,
    };
    let id = match pt.member(s.seqpnt_id) {
        Some(Value::SeqId(sid)) => Some(interner.intern(sid)),
        _ => None,
    };
    ci.extend(id);
    pieces.push(LocPiece {
        id,
        from: pos,
        to_open: pos.saturating_add(1),
        is_whole: false,
        strand,
    });
    strand
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
struct BioseqFeatureView<'a> {
    index: &'a TseFeatureIndex,
    bio_idx: usize,
    /// The posts of the title, written once the title is (`generate_asn_bioseq_title`).
    posts: RefCell<Vec<Vec<u8>>>,
}

impl BioseqFeatureView<'_> {
    fn content_label(&self, cds: &CdsFeature) -> Vec<u8> {
        self.index.cdregion_content_label(cds, &mut |msg: &[u8]| {
            self.posts.borrow_mut().push(msg.to_vec())
        })
    }
}

impl DeflineFeatures for BioseqFeatureView<'_> {
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
        let cds_label = self.content_label(cds);
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
            let Some(bi) = p.id.and_then(|r| self.index.bio_of(r)) else {
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
        let product = self.content_label(cds);
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

    // The posts are written once the title is; a refusal writes none.
    let mut posts: Vec<Vec<u8>> = Vec::new();
    let sid = classify_seq_ids_full(ids, is_virtual, &mut |w| {
        post_seq_id_warning(w, &mut |msg: &[u8]| posts.push(msg.to_vec()));
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

    // D3c-6 boundary: refuse only when NCBI's `GenerateDefline` enters `CAutoDef`: with no
    // title (or none kept for unfinished HTGS); trimming the trailing punctuation never empties
    // a title, as it keeps a title made only of it.
    // NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3952-3991
    let main_title_empty = inp.reconstruct || inp.main_title.is_empty();
    if main_title_empty && inp.is_na && has_autodef && !inp.is_tls && !inp.is_nz {
        return Err(WriteDbError::Unsupported(
            "AutoDef title generation is not supported by LOSAT's makeblastdb".to_string(),
        ));
    }

    let make = || {
        let view = BioseqFeatureView {
            index,
            bio_idx,
            posts: RefCell::new(Vec::new()),
        };
        let defline = generate_defline_core(&inp, Some(&view))
            .map(|d| (d.text, d.reads_past_end))
            .map_err(|_| ());
        (defline, view.posts.into_inner())
    };
    let Some((defline, title_posts)) = index.same_at_every_pick(make) else {
        return Err(WriteDbError::Unsupported(
            "a title made from the features that changes with the Bioseq NCBI's object manager \
             finds for a location Seq-id, which two Bioseqs of the entry match alike (NCBI picks \
             one by the order of Seq-id handles in memory) or, for a pir or prf accession alone, \
             by a name where Seq-ids that LOSAT does not read (other annotations, segments) may \
             hold another accession of its prefix, is not supported by LOSAT's makeblastdb"
                .to_string(),
        ));
    };
    match defline {
        Ok((_, true)) => Err(WriteDbError::Unsupported(
            "a title of punctuation that NCBI's x_CleanAndCompress reads past the end of \
             (NCBI makeblastdb crashes)"
                .to_string(),
        )),
        Ok((text, false)) => {
            for msg in posts.iter().chain(&title_posts) {
                warn(&msg[..]);
            }
            Ok(text)
        }
        Err(()) => {
            for msg in posts.iter().chain(&title_posts) {
                warn(&msg[..]);
            }
            // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4530-4533
            Err(WriteDbError::Exception(NcbiException::new(
                ThrowSite {
                    file: "src/corelib/ncbistr.cpp",
                    line: 4531,
                    kind: "CStringException::eBadArgs",
                    function: Some("ncbi::NStr::HtmlDecode"),
                },
                "Unable to guess the source string encoding",
            )))
        }
    }
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
