//! ASN.1 `Seq-entry` / `Bioseq-set` input source for `losat makeblastdb`
//! (`-input_type asn1_txt` and `-input_type asn1_bin`).
//!
//! NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:367-578
//! NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:290-415
//! NCBI reference (598d8ae6): c++/src/objmgr/tse_info.cpp:891-906
//! NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:725-1375,2800-3045,3928-4107

use std::collections::HashMap;
use std::io::Read;
use std::sync::{Arc, OnceLock};

use super::bioseq::SeqData;
use super::deflines::best_non_local_rank;
use super::fasta_source::pack::Ncbi4naBuilder;
use super::fasta_source::{post_seq_id_warning, Warn};
use super::WriteDbError;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::defline::{
    bioseq_facts_with_warn, generate_bioseq_defline, protein_title_suffix_adjustable,
};
use crate::report::serial::asn_ber_reader::BerReader;
use crate::report::serial::asn_input::{Names, TextInput, Utf8TagOnce};
use crate::report::serial::asn_text_reader::StreamSource;
use crate::report::serial::asn_types::read::{ObjectStream, Post, ReadError};
use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::{member_index, members, Value};
use crate::report::serial::object_frame::UNASSIGNED_MEMBER;
use crate::seqid::label;
use crate::seqid::mask_handle::{
    handle_seq_id, mask_handle_key, mask_info_key, GiimTree, TreeHandle,
};
use crate::seqid::model::ParsedSeqId;
use crate::seqid::rank::{best_rank_score, find_best};

/// Cached ASN.1 table member indices used while walking `Seq-entry` / `Bioseq-set` trees.
struct Schema {
    seq_inst_repr: usize,
    seq_inst_mol: usize,
    seq_inst_length: usize,
    seq_inst_seq_data: usize,
    seq_inst_ext: usize,
    seq_inst_hist: usize,
    seq_entry_seq: usize,
    seq_entry_set: usize,
    bioseq_set_class: usize,
    bioseq_set_descr: usize,
    bioseq_set_seq_set: usize,
    bioseq_set_annot: usize,
    bioseq_id: usize,
    bioseq_descr: usize,
    bioseq_inst: usize,
    bioseq_annot: usize,
    seq_annot_data: usize,
    seq_annot_ftable: usize,
    seq_feat_data: usize,
    feat_gene: usize,
    feat_cdregion: usize,
    feat_prot: usize,
    feat_biosrc: usize,
    seqdesc_title: usize,
    seqdesc_org: usize,
    seqdesc_source: usize,
    seqdesc_molinfo: usize,
    seqdesc_user: usize,
    seqdesc_comment: usize,
    seqdesc_genbank: usize,
    seqdesc_embl: usize,
    user_obj_type: usize,
    user_obj_data: usize,
    user_field_label: usize,
    user_field_data: usize,
    user_field_oss: usize,
    user_field_str: usize,
    object_id_id: usize,
    object_id_str: usize,
    org_ref_taxname: usize,
    org_ref_db: usize,
    org_ref_orgname: usize,
    dbtag_db: usize,
    dbtag_tag: usize,
    biosrc_genome: usize,
    biosrc_org: usize,
    biosrc_subtype: usize,
    subsrc_subtype: usize,
    subsrc_name: usize,
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
    molinfo_tech: usize,
    molinfo_completeness: usize,
    gb_block_keywords: usize,
    embl_block_keywords: usize,
    seq_hist_assembly: usize,
    seq_align_segs: usize,
    seq_align_denseg: usize,
    dense_seg_dim: usize,
    dense_seg_ids: usize,
    seq_data_iupacna: usize,
    seq_data_iupacaa: usize,
    seq_data_ncbi2na: usize,
    seq_data_ncbi4na: usize,
    seq_data_ncbieaa: usize,
    seq_data_ncbistdaa: usize,
    seq_data_gap: usize,
    seq_ext_delta: usize,
    delta_loc: usize,
    delta_literal: usize,
    seq_loc_null: usize,
    seq_lit_length: usize,
    seq_lit_data: usize,
}

fn schema() -> &'static Schema {
    static CACHE: OnceLock<Schema> = OnceLock::new();
    CACHE.get_or_init(|| {
        let seq_annot_data = member_index(named::SEQ_ANNOT, "data");
        let annot_data_ty = members(named::SEQ_ANNOT)[seq_annot_data].ty;
        let user_field_data = member_index(named::USER_FIELD, "data");
        let udata_ty = members(named::USER_FIELD)[user_field_data].ty;
        let orgname_name = member_index(named::ORGNAME, "name");
        let on_name_ty = members(named::ORGNAME)[orgname_name].ty;
        let seq_align_segs = member_index(named::SEQ_ALIGN, "segs");
        let segs_ty = members(named::SEQ_ALIGN)[seq_align_segs].ty;
        Schema {
            seq_inst_repr: member_index(named::SEQ_INST, "repr"),
            seq_inst_mol: member_index(named::SEQ_INST, "mol"),
            seq_inst_length: member_index(named::SEQ_INST, "length"),
            seq_inst_seq_data: member_index(named::SEQ_INST, "seq-data"),
            seq_inst_ext: member_index(named::SEQ_INST, "ext"),
            seq_inst_hist: member_index(named::SEQ_INST, "hist"),
            seq_entry_seq: member_index(named::SEQ_ENTRY, "seq"),
            seq_entry_set: member_index(named::SEQ_ENTRY, "set"),
            bioseq_set_class: member_index(named::BIOSEQ_SET, "class"),
            bioseq_set_descr: member_index(named::BIOSEQ_SET, "descr"),
            bioseq_set_seq_set: member_index(named::BIOSEQ_SET, "seq-set"),
            bioseq_set_annot: member_index(named::BIOSEQ_SET, "annot"),
            bioseq_id: member_index(named::BIOSEQ, "id"),
            bioseq_descr: member_index(named::BIOSEQ, "descr"),
            bioseq_inst: member_index(named::BIOSEQ, "inst"),
            bioseq_annot: member_index(named::BIOSEQ, "annot"),
            seq_annot_data,
            seq_annot_ftable: member_index(annot_data_ty, "ftable"),
            seq_feat_data: member_index(named::SEQ_FEAT, "data"),
            feat_gene: member_index(named::SEQFEATDATA, "gene"),
            feat_cdregion: member_index(named::SEQFEATDATA, "cdregion"),
            feat_prot: member_index(named::SEQFEATDATA, "prot"),
            feat_biosrc: member_index(named::SEQFEATDATA, "biosrc"),
            seqdesc_title: member_index(named::SEQDESC, "title"),
            seqdesc_org: member_index(named::SEQDESC, "org"),
            seqdesc_source: member_index(named::SEQDESC, "source"),
            seqdesc_molinfo: member_index(named::SEQDESC, "molinfo"),
            seqdesc_user: member_index(named::SEQDESC, "user"),
            seqdesc_comment: member_index(named::SEQDESC, "comment"),
            seqdesc_genbank: member_index(named::SEQDESC, "genbank"),
            seqdesc_embl: member_index(named::SEQDESC, "embl"),
            user_obj_type: member_index(named::USER_OBJECT, "type"),
            user_obj_data: member_index(named::USER_OBJECT, "data"),
            user_field_label: member_index(named::USER_FIELD, "label"),
            user_field_data,
            user_field_oss: member_index(udata_ty, "oss"),
            user_field_str: member_index(udata_ty, "str"),
            object_id_id: member_index(named::OBJECT_ID, "id"),
            object_id_str: member_index(named::OBJECT_ID, "str"),
            org_ref_taxname: member_index(named::ORG_REF, "taxname"),
            org_ref_db: member_index(named::ORG_REF, "db"),
            org_ref_orgname: member_index(named::ORG_REF, "orgname"),
            dbtag_db: member_index(named::DBTAG, "db"),
            dbtag_tag: member_index(named::DBTAG, "tag"),
            biosrc_genome: member_index(named::BIOSOURCE, "genome"),
            biosrc_org: member_index(named::BIOSOURCE, "org"),
            biosrc_subtype: member_index(named::BIOSOURCE, "subtype"),
            subsrc_subtype: member_index(named::SUBSOURCE, "subtype"),
            subsrc_name: member_index(named::SUBSOURCE, "name"),
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
            molinfo_tech: member_index(named::MOLINFO, "tech"),
            molinfo_completeness: member_index(named::MOLINFO, "completeness"),
            gb_block_keywords: member_index(named::GB_BLOCK, "keywords"),
            embl_block_keywords: member_index(named::EMBL_BLOCK, "keywords"),
            seq_hist_assembly: member_index(named::SEQ_HIST, "assembly"),
            seq_align_segs,
            seq_align_denseg: member_index(segs_ty, "denseg"),
            dense_seg_dim: member_index(named::DENSE_SEG, "dim"),
            dense_seg_ids: member_index(named::DENSE_SEG, "ids"),
            seq_data_iupacna: member_index(named::SEQ_DATA, "iupacna"),
            seq_data_iupacaa: member_index(named::SEQ_DATA, "iupacaa"),
            seq_data_ncbi2na: member_index(named::SEQ_DATA, "ncbi2na"),
            seq_data_ncbi4na: member_index(named::SEQ_DATA, "ncbi4na"),
            seq_data_ncbieaa: member_index(named::SEQ_DATA, "ncbieaa"),
            seq_data_ncbistdaa: member_index(named::SEQ_DATA, "ncbistdaa"),
            seq_data_gap: member_index(named::SEQ_DATA, "gap"),
            seq_ext_delta: member_index(named::SEQ_EXT, "delta"),
            delta_loc: member_index(named::DELTA_SEQ, "loc"),
            delta_literal: member_index(named::DELTA_SEQ, "literal"),
            seq_loc_null: member_index(named::SEQ_LOC, "null"),
            seq_lit_length: member_index(named::SEQ_LITERAL, "length"),
            seq_lit_data: member_index(named::SEQ_LITERAL, "seq-data"),
        }
    })
}

fn as_choice(val: &Value) -> Option<(usize, &Value)> {
    match val {
        Value::Choice(choice, inner) => Some((usize::from(*choice), inner.as_ref())),
        _ => None,
    }
}

fn as_choice_mut(val: &mut Value) -> Option<(usize, &mut Value)> {
    match val {
        Value::Choice(choice, inner) => Some((usize::from(*choice), inner.as_mut())),
        _ => None,
    }
}

fn as_list(val: &Value) -> Option<&[Value]> {
    match val {
        Value::List(items) => Some(items.as_slice()),
        _ => None,
    }
}

fn take_list(mut val: Option<Value>) -> Option<Vec<Value>> {
    match val.as_mut() {
        Some(Value::List(items)) => Some(std::mem::take(items)),
        _ => None,
    }
}

fn take_str(val: &mut Value) -> Option<Vec<u8>> {
    match val {
        Value::Str(bytes) => Some(std::mem::take(bytes)),
        _ => None,
    }
}

fn take_octets(val: &mut Value) -> Option<Vec<u8>> {
    match val {
        Value::Octets(bytes) => Some(std::mem::take(bytes)),
        _ => None,
    }
}

fn ncbi_err(site: ThrowSite, message: impl Into<Vec<u8>>) -> WriteDbError {
    WriteDbError::Exception(NcbiException::new(site, message))
}

fn require_literal_len(length: Option<u32>) -> Result<usize, WriteDbError> {
    length.map(|n| n as usize).ok_or_else(|| {
        WriteDbError::Exception(unassigned_member(
            "CSeq_literal::GetLength(): Attempt to get unassigned member NCBI-Sequence::Seq-literal.length",
        ))
    })
}

/// ASN.1 serialization format of `-input_type asn1_txt` or `asn1_bin`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AsnFormat {
    Text,
    Binary,
}

/// Streaming reader of consecutive `Seq-entry` / `Bioseq-set` objects from one input stream.
pub(crate) enum AsnStream<'a> {
    Text(ObjectStream<TextInput<StreamSource<&'a mut dyn Read>>>),
    Binary(ObjectStream<BerReader<StreamSource<&'a mut dyn Read>>>),
}

impl<'a> AsnStream<'a> {
    /// Opens an ASN.1 object stream over `input` in `format`.
    pub(crate) fn new(input: &'a mut dyn Read, format: AsnFormat) -> Self {
        let src = StreamSource::new(input, &[]);
        match format {
            AsnFormat::Text => AsnStream::Text(ObjectStream::text(src, Names::BIOSEQ)),
            AsnFormat::Binary => AsnStream::Binary(ObjectStream::binary(src, Names::BIOSEQ)),
        }
    }

    /// Reads the next top-level object (`Bioseq-set` when the first non-whitespace byte is `'B'`
    /// or `'0'`, else `Seq-entry`), drains reader posts to `warn`, and constructs a
    /// [`SeqEntrySource`] (including `Unverified` marking and intra-object duplicate Seq-id
    /// checks). Returns `Ok(None)` on clean EOF or top-level `CEofException::eEof`.
    ///
    /// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:437-497,559-578
    pub(crate) fn next_source(
        &mut self,
        skip_unver: bool,
        warn: Warn<'_>,
        utf8_once: &mut Utf8TagOnce,
    ) -> Result<Option<SeqEntrySource>, WriteDbError> {
        let first = match self {
            AsnStream::Text(stream) => stream.skip_raw_spaces_and_peek(),
            AsnStream::Binary(stream) => stream.skip_raw_spaces_and_peek(),
        };
        let Some(ch) = first else {
            return Ok(None);
        };
        let is_set = matches!(ch, b'B' | b'0');
        let root = if is_set {
            named::BIOSEQ_SET
        } else {
            named::SEQ_ENTRY
        };
        let (result, posts) = match self {
            AsnStream::Text(stream) => stream.read_object(root),
            AsnStream::Binary(stream) => stream.read_object(root),
        };
        write_serial_posts(posts, warn, utf8_once);
        let value = match result {
            Ok(value) => value,
            Err(ReadError::Eof) => return Ok(None),
            Err(ReadError::Serial(e) | ReadError::Exception(e)) => {
                return Err(WriteDbError::Exception(e));
            }
            Err(ReadError::OutOfMemory) => return Err(WriteDbError::OutOfMemory),
            Err(ReadError::Std(what)) => return Err(WriteDbError::Other(what.to_string())),
            Err(ReadError::Unsupported(what)) => return Err(WriteDbError::Unsupported(what)),
        };
        SeqEntrySource::from_root(value, is_set, skip_unver).map(Some)
    }
}

/// Posts serial reader diagnostics in order, deduplicating warning 803.10 across the run.
pub(crate) fn write_serial_posts(posts: Vec<Post>, warn: Warn<'_>, utf8_once: &mut Utf8TagOnce) {
    for post in posts {
        if matches!(post, Post::Utf8Tag(_)) && !utf8_once.first() {
            continue;
        }
        let line = post.line(Some("makeblastdb"));
        warn(line.strip_suffix(b"\n").unwrap_or(&line));
    }
}

/// One element of an ASN.1 `Delta-ext` sequence, extracted with sequence data moved out of `Value`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum AsnDeltaItem {
    NullElement,
    LocNull,
    LocOther,
    LiteralNoData(Option<u32>),
    LiteralNcbi2na { length: Option<u32>, bytes: Vec<u8> },
    LiteralNcbi4na { length: Option<u32>, bytes: Vec<u8> },
    LiteralGap { length: Option<u32> },
    LiteralUnhandled,
}

/// A `CBioseq` yielded by [`SeqEntrySource::next_bioseq`] after `CDeflineGenerator` and
/// `s_GetTaxId`, ready for `CBuildDatabase::AddSequences`.
#[derive(Debug)]
pub(crate) struct AsnBioseq {
    pub(crate) ids: Vec<ParsedSeqId>,
    pub(crate) hist_old_id: Option<ParsedSeqId>,
    pub(crate) is_aa: bool,
    pub(crate) is_na: bool,
    pub(crate) length: u64,
    pub(crate) seq_data: Option<SeqData>,
    pub(crate) delta: Option<Vec<AsnDeltaItem>>,
    pub(crate) binary_header: Option<Vec<u8>>,
    pub(crate) cfastareader_defline: Option<Vec<u8>>,
    pub(crate) title: Vec<u8>,
    pub(crate) descr_taxids: Vec<i32>,
}

impl AsnBioseq {
    /// `GetId().front()->AsFastaString()` (or `"Unknown"` when the Bioseq has no Seq-ids).
    pub(crate) fn first_id_fasta(&self) -> String {
        match self.ids.first() {
            Some(id) => {
                let mut out = Vec::new();
                label::fasta(id, &mut out);
                String::from_utf8_lossy(&out).into_owned()
            }
            None => "Unknown".to_string(),
        }
    }
}

/// Summary of `Seq-annot` feature table kinds present anywhere inside the top-level `Seq-entry`
/// (used to detect feature-driven `CDeflineGenerator` paths that belong to D3c-5/D3c-6).
#[derive(Debug, Default, Clone, Copy)]
struct TseAnnotSummary {
    has_prot: bool,
    has_cdregion: bool,
    has_gene: bool,
    has_biosrc: bool,
    has_biosrc_clone: bool,
}

/// Descriptors of an ancestor `Bioseq-set` along the parent chain of a `Bioseq`.
#[derive(Debug)]
struct SetNode {
    parent: Option<usize>,
    descr: Vec<Value>,
}

/// One `Bioseq` discovered during the depth-first pre-order walk of a top-level `Seq-entry`.
#[derive(Debug)]
struct BioseqSlot {
    parent_set: Option<usize>,
    unverified: bool,
    ids: Vec<ParsedSeqId>,
    descr: Vec<Value>,
    inst: Value,
}

/// `CSeqEntrySource`: holds one parsed top-level `Seq-entry` / `Bioseq-set` and yields its
/// Bioseqs one at a time while moving large sequence buffers out without cloning.
pub(crate) struct SeqEntrySource {
    sets: Vec<SetNode>,
    bioseqs: Vec<Option<BioseqSlot>>,
    cursor: usize,
    skip_unver: bool,
    annots: TseAnnotSummary,
}

impl SeqEntrySource {
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:466-497
    // NCBI reference (598d8ae6): c++/src/objmgr/tse_info.cpp:891-906
    fn from_root(root: Value, is_set: bool, skip_unver: bool) -> Result<Self, WriteDbError> {
        let mut sets = Vec::new();
        let mut bioseqs = Vec::new();
        let mut annots = TseAnnotSummary::default();
        collect_tree(root, is_set, &mut sets, &mut bioseqs, &mut annots)?;
        check_duplicate_bioseq_ids(&bioseqs)?;
        Ok(Self {
            sets,
            bioseqs: bioseqs.into_iter().map(Some).collect(),
            cursor: 0,
            skip_unver,
            annots,
        })
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:499-548
    /// `CSeqEntrySource::GetNext`: skips marked unverified Bioseqs when `skip_unver` is true,
    /// runs `CDeflineGenerator::GenerateDefline` and `s_GetTaxId`, and returns the next Bioseq.
    pub(crate) fn next_bioseq(
        &mut self,
        warn: Warn<'_>,
    ) -> Result<Option<AsnBioseq>, WriteDbError> {
        while self.cursor < self.bioseqs.len() {
            let slot = self.bioseqs[self.cursor].take();
            self.cursor += 1;
            let Some(slot) = slot else {
                continue;
            };
            if self.skip_unver && slot.unverified {
                continue;
            }
            return Ok(Some(self.materialize_bioseq(slot, warn)?));
        }
        Ok(None)
    }

    fn materialize_bioseq(
        &self,
        mut slot: BioseqSlot,
        warn: Warn<'_>,
    ) -> Result<AsnBioseq, WriteDbError> {
        let s = schema();

        // NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:729-733
        // `x_SetFlags` calls `bsh.GetInst_Length()` before returning to `AddSequences`, throwing
        // `CUnassignedMember::eGet` whenever `Seq-inst.length` is omitted.
        let Some(Value::Int(length)) = slot.inst.member(s.seq_inst_length) else {
            return Err(WriteDbError::Exception(unassigned_member(
                "CSeq_inst::GetLength(): Attempt to get unassigned member NCBI-Sequence::Seq-inst.length",
            )));
        };
        let length = *length as u64;
        let repr = match slot.inst.member(s.seq_inst_repr) {
            Some(Value::Enum(v)) => *v,
            _ => 0,
        };
        let mol = match slot.inst.member(s.seq_inst_mol) {
            Some(Value::Enum(v)) => *v,
            _ => 0,
        };
        let is_aa = mol == 3;
        let is_na = matches!(mol, 1 | 2 | 4);

        let generated_title = generate_asn_title(
            &slot.ids,
            &slot.descr,
            self.ancestor_descrs(slot.parent_set),
            is_aa,
            is_na,
            repr,
            self.annots,
            warn,
        )?;

        // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:520-532
        // `s_HasTitle`: if Bioseq has any title descriptor, overwrite every title descriptor on
        // the Bioseq with `generated_title`; otherwise append a new title descriptor.
        let mut had_title = false;
        for desc in &mut slot.descr {
            if let Some((idx, inner)) = as_choice_mut(desc) {
                if idx == s.seqdesc_title {
                    *inner = Value::Str(generated_title.clone());
                    had_title = true;
                }
            }
        }
        if !had_title {
            slot.descr.push(Value::Choice(
                s.seqdesc_title as u16,
                Box::new(Value::Str(generated_title)),
            ));
        }

        // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:534-541
        // When the Bioseq's own taxid is zero, walk `s_GetTaxId` and append an `org` descriptor
        // when non-zero.
        let mut extra_taxid = None;
        if bioseq_own_taxid(&slot.descr) == 0 {
            let inherited = s_get_taxid(&slot.descr, self.ancestor_descrs(slot.parent_set));
            if inherited != 0 {
                extra_taxid = Some(inherited);
            }
        }

        let binary_header = extract_binary_header(&slot.descr);
        let cfastareader_defline = extract_cfastareader_defline(&slot.descr);
        let mut descr_taxids = collect_descr_taxids(&slot.descr);
        if let Some(taxid) = extra_taxid {
            descr_taxids.push(taxid);
        }
        // `x_BuildDeflinesFromBioseq` takes the last title descriptor on the Bioseq.
        let title = slot
            .descr
            .iter()
            .rev()
            .find_map(|desc| match as_choice(desc) {
                Some((idx, Value::Str(text))) if idx == s.seqdesc_title => Some(text.clone()),
                _ => None,
            })
            .unwrap_or_default();

        let hist_old_id = slot
            .inst
            .member(s.seq_inst_hist)
            .and_then(|hist| extract_hist_old_id(hist, &slot.ids));

        // Move `Seq-data` / `Delta-ext` buffers out of `Value` without cloning.
        let seq_data = slot
            .inst
            .take_member(s.seq_inst_seq_data)
            .map(value_to_seq_data);
        let delta = if seq_data.is_none() {
            slot.inst
                .take_member(s.seq_inst_ext)
                .and_then(extract_delta_ext)
        } else {
            None
        };

        Ok(AsnBioseq {
            ids: slot.ids,
            hist_old_id,
            is_aa,
            is_na,
            length,
            seq_data,
            delta,
            binary_header,
            cfastareader_defline,
            title,
            descr_taxids,
        })
    }

    fn ancestor_descrs(&self, mut parent: Option<usize>) -> Vec<&[Value]> {
        let mut out = Vec::new();
        while let Some(idx) = parent {
            out.push(self.sets[idx].descr.as_slice());
            parent = self.sets[idx].parent;
        }
        out
    }
}

/// Walks a parsed `Seq-entry` or `Bioseq-set` iteratively (no recursion) to collect all Bioseqs,
/// ancestor descriptor chains, `Unverified` flags, and annotation summaries.
///
/// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:466-496
fn collect_tree(
    root: Value,
    is_set: bool,
    sets: &mut Vec<SetNode>,
    bioseqs: &mut Vec<BioseqSlot>,
    annots: &mut TseAnnotSummary,
) -> Result<(), WriteDbError> {
    enum Work {
        Entry {
            val: Value,
            parent_set: Option<usize>,
            inherited_unver: bool,
        },
        Set {
            val: Value,
            parent_set: Option<usize>,
            inherited_unver: bool,
        },
    }

    let s = schema();

    let mut stack = vec![if is_set {
        Work::Set {
            val: root,
            parent_set: None,
            inherited_unver: false,
        }
    } else {
        Work::Entry {
            val: root,
            parent_set: None,
            inherited_unver: false,
        }
    }];

    while let Some(item) = stack.pop() {
        match item {
            Work::Entry {
                mut val,
                parent_set,
                inherited_unver,
            } => {
                if matches!(val, Value::Null) {
                    return Err(null_element_refusal("Bioseq-set.seq-set"));
                }
                let Some((choice, inner)) = as_choice_mut(&mut val) else {
                    continue;
                };
                let inner = inner.take();
                if choice == s.seq_entry_seq {
                    let mut bio = inner;
                    let ids = extract_bioseq_ids(bio.take_member(s.bioseq_id))?;
                    let descr = extract_descr_list(bio.take_member(s.bioseq_descr))?;
                    if let Some(annot_val) = bio.member(s.bioseq_annot) {
                        scan_annots(annot_val, annots)?;
                    }
                    let own_unver = is_unverified_descr(&descr);
                    let inst = bio.take_member(s.bioseq_inst).unwrap_or(Value::Null);
                    bioseqs.push(BioseqSlot {
                        parent_set,
                        unverified: inherited_unver || own_unver,
                        ids,
                        descr,
                        inst,
                    });
                } else if choice == s.seq_entry_set {
                    stack.push(Work::Set {
                        val: inner,
                        parent_set,
                        inherited_unver,
                    });
                }
            }
            Work::Set {
                mut val,
                parent_set,
                inherited_unver,
            } => {
                let class = match val.member(s.bioseq_set_class) {
                    Some(Value::Enum(c)) => *c,
                    _ => 0,
                };
                let descr = extract_descr_list(val.take_member(s.bioseq_set_descr))?;
                if let Some(annot_val) = val.member(s.bioseq_set_annot) {
                    scan_annots(annot_val, annots)?;
                }
                // `CBioseq_set::eClass_nuc_prot` (1): Unverified on a nuc-prot set marks all
                // Bioseqs beneath it (makeblastdb.cpp:474-480).
                let nuc_prot_unver = class == 1 && is_unverified_descr(&descr);
                let set_idx = sets.len();
                sets.push(SetNode {
                    parent: parent_set,
                    descr,
                });
                if let Some(children) = take_list(val.take_member(s.bioseq_set_seq_set)) {
                    for child in children.into_iter().rev() {
                        stack.push(Work::Entry {
                            val: child,
                            parent_set: Some(set_idx),
                            inherited_unver: inherited_unver || nuc_prot_unver,
                        });
                    }
                }
            }
        }
    }
    Ok(())
}

fn null_element_refusal(where_in: &str) -> WriteDbError {
    WriteDbError::Unsupported(format!(
        "a NULL element in {where_in} (NCBI makeblastdb dereferences a NULL pointer there)"
    ))
}

fn extract_bioseq_ids(val: Option<Value>) -> Result<Vec<ParsedSeqId>, WriteDbError> {
    let Some(items) = take_list(val) else {
        return Ok(Vec::new());
    };
    let mut ids = Vec::with_capacity(items.len());
    for item in items {
        if matches!(item, Value::Null) {
            return Err(null_element_refusal("Bioseq.id"));
        }
        if let Some(id) = item.into_seq_id() {
            ids.push(id);
        }
    }
    Ok(ids)
}

fn extract_descr_list(val: Option<Value>) -> Result<Vec<Value>, WriteDbError> {
    let Some(items) = take_list(val) else {
        return Ok(Vec::new());
    };
    if items.iter().any(|v| matches!(v, Value::Null)) {
        return Err(null_element_refusal("Seq-descr"));
    }
    Ok(items)
}

fn scan_annots(val: &Value, summary: &mut TseAnnotSummary) -> Result<(), WriteDbError> {
    let Some(annots) = as_list(val) else {
        return Ok(());
    };
    let s = schema();

    for annot in annots {
        if matches!(annot, Value::Null) {
            return Err(null_element_refusal("Seq-annot list"));
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
                return Err(null_element_refusal("Seq-annot.data.ftable"));
            }
            let Some((fchoice, fval)) = feat.member(s.seq_feat_data).and_then(as_choice) else {
                continue;
            };
            if fchoice == s.feat_prot {
                summary.has_prot = true;
            } else if fchoice == s.feat_cdregion {
                summary.has_cdregion = true;
            } else if fchoice == s.feat_gene {
                summary.has_gene = true;
            } else if fchoice == s.feat_biosrc {
                summary.has_biosrc = true;
                if let Some(subs) = fval.member(s.biosrc_subtype).and_then(as_list) {
                    for sbs in subs {
                        let is_clone = matches!(sbs.member(s.subsrc_subtype), Some(Value::Enum(3)));
                        if is_clone && sbs.member(s.subsrc_name).is_some() {
                            summary.has_biosrc_clone = true;
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

/// `CSeqEntrySource::IsUnverified`: whether `descr` holds a `User-object` whose `type` is a
/// string matching `"Unverified"` without ASCII case.
///
/// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:413-435
fn is_unverified_descr(descr: &[Value]) -> bool {
    let s = schema();
    descr.iter().any(|desc| {
        as_choice(desc)
            .filter(|(c, _)| *c == s.seqdesc_user)
            .and_then(|(_, inner)| inner.member(s.user_obj_type))
            .and_then(object_id_str)
            .is_some_and(|name| name.eq_ignore_ascii_case(b"Unverified"))
    })
}

/// `CTSE_Info::x_SetBioseqIds`: checks intra-object duplicate Seq-ids under `CSeq_id_Handle`
/// identity (`seqid::mask_handle`), throwing NCBI's exact `CObjMgrException::eAddDataError`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/tse_info.cpp:891-906
fn check_duplicate_bioseq_ids(bioseqs: &[BioseqSlot]) -> Result<(), WriteDbError> {
    const DUP_SITE: ThrowSite = ThrowSite {
        file: "src/objmgr/tse_info.cpp",
        line: 901,
        kind: "CObjMgrException::eAddDataError",
        function: Some("ncbi::objects::CTSE_Info::x_SetBioseqIds"),
    };

    #[derive(Clone, PartialEq, Eq, Hash)]
    enum HandleKey {
        Bytes(Vec<u8>),
        Giim(i32, usize),
    }

    let mut info_first: HashMap<Vec<u8>, ParsedSeqId> = HashMap::new();
    let mut giim = GiimTree::default();
    let mut seen: HashMap<HandleKey, usize> = HashMap::new();
    let mut id_strings: Vec<Vec<u8>> = Vec::with_capacity(bioseqs.len());

    for (b_idx, slot) in bioseqs.iter().enumerate() {
        let mut handles: Vec<(Option<HandleKey>, Vec<u8>)> = Vec::with_capacity(slot.ids.len());
        let mut joined = Vec::new();
        for (i, id) in slot.ids.iter().enumerate() {
            let (hkey, shown_id) = match mask_handle_key(id) {
                Ok(Some(key)) => {
                    let info_key = mask_info_key(id)
                        .expect("non-giim info key")
                        .expect("non-empty info key");
                    let first = info_first.entry(info_key).or_insert_with(|| id.clone());
                    (Some(HandleKey::Bytes(key)), handle_seq_id(first, id))
                }
                Ok(None) => (None, id.clone()),
                Err(TreeHandle::Giim { id: num }) => {
                    let ParsedSeqId::Giim(g) = id else {
                        unreachable!("TreeHandle::Giim only comes from ParsedSeqId::Giim");
                    };
                    let (h, _) = giim.find_or_create(g);
                    let first = ParsedSeqId::Giim(giim.made_from(num, h).clone());
                    (Some(HandleKey::Giim(num, h)), handle_seq_id(&first, id))
                }
            };
            let mut fasta = Vec::new();
            label::fasta(&shown_id, &mut fasta);
            if i > 0 {
                joined.extend_from_slice(b" | ");
            }
            joined.extend_from_slice(&fasta);
            handles.push((hkey, fasta));
        }
        id_strings.push(joined);

        for (hkey, fasta) in handles {
            let Some(hkey) = hkey else {
                continue;
            };
            if let Some(&prev_idx) = seen.get(&hkey) {
                let mut msg = Vec::new();
                msg.extend_from_slice(b"duplicate Bioseq id ");
                msg.extend_from_slice(&fasta);
                msg.extend_from_slice(b" present in\n  seq1: ");
                msg.extend_from_slice(&id_strings[prev_idx]);
                msg.extend_from_slice(b"\n  seq2: ");
                msg.extend_from_slice(&id_strings[b_idx]);
                return Err(WriteDbError::Exception(NcbiException::new(DUP_SITE, msg)));
            }
            seen.insert(hkey, b_idx);
        }
    }
    Ok(())
}

/// `CUnassignedMember::eGet` at `src/serial/serialobject.cpp:306` (`object_frame::UNASSIGNED_MEMBER`).
fn unassigned_member(message: &str) -> NcbiException {
    NcbiException::new(UNASSIGNED_MEMBER, message)
}

fn object_id_str(val: &Value) -> Option<&[u8]> {
    let s = schema();
    match as_choice(val) {
        Some((choice, Value::Str(bytes))) if choice == s.object_id_str => Some(bytes),
        _ => None,
    }
}

fn object_id_int(val: &Value) -> Option<i64> {
    let s = schema();
    match as_choice(val) {
        Some((choice, Value::Int(n))) if choice == s.object_id_id => Some(*n),
        _ => None,
    }
}

/// `COrg_ref::GetTaxId`: the first `taxon` dbtag whose tag is an integer (`id`), else 0.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Org_ref.cpp:71-90
fn org_taxids(org: &Value) -> impl Iterator<Item = i32> + '_ {
    let s = schema();
    let tags = org.member(s.org_ref_db).and_then(as_list).unwrap_or(&[]);
    tags.iter().filter_map(move |tag| {
        if !matches!(tag.member(s.dbtag_db), Some(Value::Str(db)) if db.as_slice() == b"taxon") {
            return None;
        }
        object_id_int(tag.member(s.dbtag_tag)?).map(|id| id as i32)
    })
}

/// `CBioseq::GetTaxId`: scans Bioseq's own descriptors; `source` takes precedence over `org`,
/// stopping at the first non-zero `source` taxid, else returning the last `org` taxid (or 0).
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/Bioseq.cpp:177-199
fn bioseq_own_taxid(descr: &[Value]) -> i32 {
    let s = schema();
    let mut from_source = 0i32;
    let mut from_org = 0i32;
    for desc in descr {
        let Some((choice, inner)) = as_choice(desc) else {
            continue;
        };
        if choice == s.seqdesc_org {
            from_org = org_taxids(inner).next().unwrap_or(0);
        } else if choice == s.seqdesc_source {
            if let Some(org) = inner.member(s.biosrc_org) {
                from_source = org_taxids(org).next().unwrap_or(0);
            }
        }
        if from_source != 0 {
            break;
        }
    }
    if from_source != 0 {
        from_source
    } else {
        from_org
    }
}

/// `s_GetTaxId`: starts at the Bioseq's own descriptors (`bio.GetParentEntry()`), then walks
/// ancestor `Bioseq-set` descriptors upward, returning the `GetTaxId()` of the very first
/// `source` or `org` descriptor found (even if 0).
///
/// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:367-389
fn s_get_taxid(bio_descr: &[Value], ancestors: Vec<&[Value]>) -> i32 {
    let s = schema();
    for level in std::iter::once(bio_descr).chain(ancestors) {
        for (tag, inner) in level.iter().filter_map(as_choice) {
            if tag == s.seqdesc_source {
                return inner
                    .member(s.biosrc_org)
                    .and_then(|org| org_taxids(org).next())
                    .unwrap_or(0);
            }
            if tag == s.seqdesc_org {
                return org_taxids(inner).next().unwrap_or(0);
            }
        }
    }
    0
}

/// Collects all `taxon` dbtag integer IDs from `org` and `source` descriptors on the Bioseq, in
/// order (`x_BuildDeflinesFromBioseq`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:608-640
fn collect_descr_taxids(descr: &[Value]) -> Vec<i32> {
    let s = schema();
    let mut taxids = Vec::new();
    for (choice, inner) in descr.iter().filter_map(as_choice) {
        let org = if choice == s.seqdesc_org {
            Some(inner)
        } else if choice == s.seqdesc_source {
            inner.member(s.biosrc_org)
        } else {
            None
        };
        if let Some(org) = org {
            taxids.extend(org_taxids(org));
        }
    }
    taxids
}

/// `x_GetBioseqBinaryHeader`: first non-empty octet string of the first `ASN1_BlastDefLine` user
/// object whose first field has label `ASN1_BlastDefLine` and `data.oss`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:441-489
fn extract_binary_header(descr: &[Value]) -> Option<Vec<u8>> {
    const LABEL: &[u8] = b"ASN1_BlastDefLine";
    let s = schema();
    for (choice, obj) in descr.iter().filter_map(as_choice) {
        if choice != s.seqdesc_user {
            continue;
        }
        if obj.member(s.user_obj_type).and_then(object_id_str) != Some(LABEL) {
            continue;
        }
        let Some(first) = obj
            .member(s.user_obj_data)
            .and_then(as_list)
            .and_then(|f| f.first())
        else {
            continue;
        };
        if first.member(s.user_field_label).and_then(object_id_str) != Some(LABEL) {
            continue;
        }
        if let Some((fchoice, data)) = first.member(s.user_field_data).and_then(as_choice) {
            if fchoice == s.user_field_oss {
                if let Some(Value::Octets(bytes)) = as_list(data).and_then(|oss| oss.first()) {
                    if !bytes.is_empty() {
                        return Some(bytes.clone());
                    }
                }
                return None;
            }
        }
    }
    None
}

/// Extracts the raw `DefLine` string from the last `CFastaReader` user object on the Bioseq
/// (first `DefLine` string field of each `CFastaReader` user object overwrites earlier ones).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1737-1765
fn extract_cfastareader_defline(descr: &[Value]) -> Option<Vec<u8>> {
    let s = schema();
    let mut found = None;
    for (tag, obj) in descr.iter().filter_map(as_choice) {
        if tag != s.seqdesc_user {
            continue;
        }
        if obj.member(s.user_obj_type).and_then(object_id_str) != Some(b"CFastaReader") {
            continue;
        }
        let Some(fields) = obj.member(s.user_obj_data).and_then(as_list) else {
            continue;
        };
        for field in fields {
            if field.member(s.user_field_label).and_then(object_id_str) != Some(b"DefLine") {
                continue;
            }
            if let Some((fchoice, Value::Str(text))) =
                field.member(s.user_field_data).and_then(as_choice)
            {
                if fchoice == s.user_field_str {
                    found = Some(text.clone());
                    break;
                }
            }
        }
    }
    found
}

/// Extracts the non-local Seq-id paired with a local Bioseq ID in `Seq-hist.assembly`
/// (`CBioseq::GetNonLocalId` fallback for `s_UseFastaReaderDeflines`).
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/Bioseq.cpp:292-325
fn extract_hist_old_id(hist: &Value, bio_ids: &[ParsedSeqId]) -> Option<ParsedSeqId> {
    let best_local = find_best(bio_ids, best_non_local_rank)?;
    if !matches!(best_local, ParsedSeqId::Local(_)) {
        return None;
    }
    let s = schema();
    let aligns = hist.member(s.seq_hist_assembly).and_then(as_list)?;
    for align in aligns {
        let Some((choice, seg)) = align.member(s.seq_align_segs).and_then(as_choice) else {
            continue;
        };
        if choice != s.seq_align_denseg {
            continue;
        }
        let dim = match seg.member(s.dense_seg_dim) {
            Some(Value::Int(d)) => *d,
            None => 2,
            _ => continue,
        };
        if dim != 2 {
            continue;
        }
        let Some(ids) = seg.member(s.dense_seg_ids).and_then(as_list) else {
            continue;
        };
        if ids.len() != 2 {
            continue;
        }
        let (Value::SeqId(id1), Value::SeqId(id2)) = (&ids[0], &ids[1]) else {
            continue;
        };
        let eq_local = |a: &ParsedSeqId| {
            matches!(a, ParsedSeqId::Local(_))
                && mask_handle_key(a).ok().flatten() == mask_handle_key(best_local).ok().flatten()
        };
        if eq_local(id1) && !matches!(id2.as_ref(), ParsedSeqId::Local(_)) {
            return Some((**id2).clone());
        }
        if eq_local(id2) && !matches!(id1.as_ref(), ParsedSeqId::Local(_)) {
            return Some((**id1).clone());
        }
    }
    None
}

/// Moves `Seq-data` bytes out of a parsed `Value` without copying.
fn value_to_seq_data(mut val: Value) -> SeqData {
    let s = schema();
    let Some((idx, inner)) = as_choice_mut(&mut val) else {
        return SeqData::Other;
    };
    if idx == s.seq_data_iupacna {
        return take_str(inner).map_or(SeqData::Other, SeqData::Iupacna);
    }
    if idx == s.seq_data_iupacaa {
        return take_str(inner).map_or(SeqData::Other, SeqData::Iupacaa);
    }
    if idx == s.seq_data_ncbi2na {
        return take_octets(inner).map_or(SeqData::Other, SeqData::Ncbi2na);
    }
    if idx == s.seq_data_ncbi4na {
        return take_octets(inner).map_or(SeqData::Other, SeqData::Ncbi4na);
    }
    if idx == s.seq_data_ncbieaa {
        return take_str(inner).map_or(SeqData::Other, SeqData::Ncbieaa);
    }
    if idx == s.seq_data_ncbistdaa {
        return take_octets(inner).map_or(SeqData::Other, SeqData::Ncbistdaa);
    }
    SeqData::Other
}

/// Extracts `Delta-ext` items from `Seq-inst.ext`, moving literal byte buffers out of `Value`.
fn extract_delta_ext(mut ext: Value) -> Option<Vec<AsnDeltaItem>> {
    let s = schema();
    let (choice, inner) = as_choice_mut(&mut ext)?;
    if choice != s.seq_ext_delta {
        return None;
    }
    let items = take_list(Some(inner.take()))?;
    let mut out = Vec::with_capacity(items.len());
    for mut item in items {
        if matches!(item, Value::Null) {
            out.push(AsnDeltaItem::NullElement);
            continue;
        }
        let Some((dchoice, dval)) = as_choice_mut(&mut item) else {
            continue;
        };
        if dchoice == s.delta_loc {
            let is_null = matches!(as_choice(dval), Some((lc, _)) if lc == s.seq_loc_null);
            out.push(if is_null {
                AsnDeltaItem::LocNull
            } else {
                AsnDeltaItem::LocOther
            });
        } else if dchoice == s.delta_literal {
            let mut lit = dval.take();
            let length = match lit.member(s.seq_lit_length) {
                Some(Value::Int(n)) => Some(*n as u32),
                _ => None,
            };
            let entry = match lit.take_member(s.seq_lit_data) {
                None => AsnDeltaItem::LiteralNoData(length),
                Some(mut data) => match as_choice_mut(&mut data) {
                    Some((c, v)) if c == s.seq_data_ncbi2na => match take_octets(v) {
                        Some(bytes) => AsnDeltaItem::LiteralNcbi2na { length, bytes },
                        None => AsnDeltaItem::LiteralUnhandled,
                    },
                    Some((c, v)) if c == s.seq_data_ncbi4na => match take_octets(v) {
                        Some(bytes) => AsnDeltaItem::LiteralNcbi4na { length, bytes },
                        None => AsnDeltaItem::LiteralUnhandled,
                    },
                    Some((c, _)) if c == s.seq_data_gap => AsnDeltaItem::LiteralGap { length },
                    _ => AsnDeltaItem::LiteralUnhandled,
                },
            };
            out.push(entry);
        }
    }
    Some(out)
}

/// `s_FixBioseqDeltas` for an ASN.1 `Bioseq`: expands a nucleotide literal-only delta sequence
/// directly into packed NCBI4na bytes (`Ncbi4naBuilder`), matching every NCBI throw site and
/// wrapping the resulting length to 32 bits.
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:290-415
pub(crate) fn fix_asn_delta(bioseq: &mut AsnBioseq) -> Result<(), WriteDbError> {
    if bioseq.seq_data.is_some() {
        return Ok(());
    }
    let Some(items) = bioseq.delta.take() else {
        return Ok(());
    };
    const fn delta_site(line: u32) -> ThrowSite {
        ThrowSite {
            file: "src/objtools/blast/seqdb_writer/build_db.cpp",
            line,
            kind: "CMultisourceException::eArgErr",
            function: Some("BLASTDB::ncbi::s_FixBioseqDeltas"),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:297-304
    if !bioseq.is_na {
        let mut best = Vec::new();
        if let Some(id) = find_best(&bioseq.ids, best_rank_score) {
            label::fasta(id, &mut best);
        }
        return Err(ncbi_err(
            delta_site(303),
            [
                best.as_slice(),
                b": Protein delta sequences are not supported.",
            ]
            .concat(),
        ));
    }

    let Some(first) = items.first() else {
        return Err(WriteDbError::Unsupported(
            "an empty Delta-ext list (NCBI makeblastdb dereferences an empty list there)"
                .to_string(),
        ));
    };
    if matches!(first, AsnDeltaItem::LocNull | AsnDeltaItem::LocOther) {
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:309-310
        return Ok(());
    }

    // Validate items in NCBI's exact iteration order and compute total letter count before
    // reserving the single NCBI4na output buffer (avoiding reallocation and handling OOM cleanly).
    let mut total_letters = 0usize;
    for item in &items {
        let add = match item {
            AsnDeltaItem::NullElement => return Err(null_element_refusal("Delta-ext")),
            AsnDeltaItem::LocNull => 1usize,
            AsnDeltaItem::LocOther => {
                return Err(ncbi_err(
                    delta_site(409),
                    "Bioseq must have Seq-data or Delta containing only literals.",
                ));
            }
            AsnDeltaItem::LiteralNoData(Some(len)) => *len as usize,
            AsnDeltaItem::LiteralNoData(None) => {
                return Err(ncbi_err(
                    delta_site(338),
                    "Part of the delta sequence, including its length, is un-available.",
                ));
            }
            AsnDeltaItem::LiteralNcbi2na { length, bytes } => {
                require_literal_len(*length)?.min(bytes.len().saturating_mul(4))
            }
            AsnDeltaItem::LiteralNcbi4na { length, bytes } => {
                require_literal_len(*length)?.min(bytes.len().saturating_mul(2))
            }
            AsnDeltaItem::LiteralGap { length } => require_literal_len(*length)?,
            AsnDeltaItem::LiteralUnhandled => {
                return Err(ncbi_err(
                    delta_site(360),
                    "Unhandled type of sequence data encountered.",
                ));
            }
        };
        total_letters = total_letters
            .checked_add(add)
            .ok_or(WriteDbError::OutOfMemory)?;
    }

    let mut builder =
        Ncbi4naBuilder::with_letters(total_letters).map_err(|_| WriteDbError::OutOfMemory)?;
    for item in items {
        match item {
            AsnDeltaItem::LocNull => builder.push_run(0x0f, 1),
            AsnDeltaItem::LiteralNoData(Some(len))
            | AsnDeltaItem::LiteralGap { length: Some(len) } => {
                builder.push_run(0x0f, len as usize);
            }
            AsnDeltaItem::LiteralNcbi2na {
                length: Some(len),
                bytes,
            } => builder.push_ncbi2na(&bytes, len as usize),
            AsnDeltaItem::LiteralNcbi4na {
                length: Some(len),
                bytes,
            } => builder.push_ncbi4na(&bytes, 0, len as usize),
            _ => {}
        }
    }
    let letters = builder.letters() as u32;
    bioseq.seq_data = Some(SeqData::Ncbi4na(builder.into_bytes()));
    bioseq.length = u64::from(letters);
    Ok(())
}

/// Generates a Bioseq title using `report::defline::generate_bioseq_defline` for all cases within
/// D3c-4 scope, and explicitly refuses D3c-5/D3c-6 title generation paths (`MolInfo`
/// prefix/suffix, `Unverified`/`Unreviewed` prefix, feature-derived protein names, `BioSource`
/// modifiers, specialized RefSeq/WGS/TSA/TLS/seg/map/PDB-compound/AutoDef paths).
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:725-1375,2800-3045,3928-4107
#[allow(clippy::too_many_arguments)]
fn generate_asn_title(
    ids: &[ParsedSeqId],
    bio_descr: &[Value],
    ancestors: Vec<&[Value]>,
    is_aa: bool,
    is_na: bool,
    repr: i32,
    annots: TseAnnotSummary,
    warn: Warn<'_>,
) -> Result<Vec<u8>, WriteDbError> {
    let is_seg = repr == 3;
    let is_map = repr == 7;
    let is_delta = repr == 8;

    let mut facts = bioseq_facts_with_warn(ids, &mut |w| {
        post_seq_id_warning(w, warn);
    });

    let flags = scan_defline_descriptors(
        bio_descr,
        &ancestors,
        facts.pdb.is_some(),
        facts.refseq_unique_prot,
    );
    facts.has_source = flags.has_source;
    facts.taxname = flags
        .taxname
        .as_deref()
        .filter(|t| !t.is_empty())
        .map(Arc::from);

    let main_trimmed = trim_end_slice(&flags.raw_title, b".,;~ ");
    let has_explicit = !main_trimmed.is_empty();

    // Check prefix / suffix triggers that belong to D3c-5/D3c-6.
    if flags.molinfo_tech_unsupported {
        return Err(title_refusal("MolInfo tech prefix/suffix title generation"));
    }
    if flags.is_unverified {
        return Err(title_refusal("Unverified title prefix generation"));
    }
    if flags.is_unreviewed {
        return Err(title_refusal("Unreviewed title prefix generation"));
    }
    if facts.third_party && flags.has_tpa_keyword {
        return Err(title_refusal("TPA keyword title prefix generation"));
    }
    if facts.refseq_unique_prot && flags.is_multispecies {
        return Err(title_refusal("MULTISPECIES WP title prefix generation"));
    }
    if flags.is_pseudogene {
        return Err(title_refusal("pseudogene comment title prefix generation"));
    }
    if flags.has_metagenome_source {
        return Err(title_refusal(
            "MAG metagenome-source title prefix generation",
        ));
    }
    if flags.has_unordered_kw && is_delta {
        return Err(title_refusal("UNORDERED delta title suffix generation"));
    }

    if has_explicit {
        if is_aa
            && protein_title_suffix_adjustable(main_trimmed)
            && (!flags.is_complete
                || is_protein_organelle(flags.genome)
                || flags.has_binomial
                || flags.is_cross_kingdom
                || annots.has_biosrc)
        {
            return Err(title_refusal(
                "protein title suffix adjustment with partial/organelle/binomial/feature modifiers",
            ));
        }
    } else {
        // Empty title path: check D3c-5/D3c-6 generators before delegating to report::defline.
        if is_na && flags.has_autodef && !facts.refseq_contig {
            return Err(title_refusal("AutoDef title generation"));
        }
        if facts.pdb.is_some() {
            return Err(title_refusal("empty-title PDB title generation"));
        }
        if facts.patent.is_some() {
            return Err(title_refusal("empty-title patent title generation"));
        }
        if facts.refseq_chromosome_or_ncrna || facts.refseq_mrna {
            return Err(title_refusal(
                "empty-title RefSeq NC/NM/NR title generation",
            ));
        }
        if is_seg || is_map {
            return Err(title_refusal("empty-title segmented/map title generation"));
        }
        if is_aa {
            if !flags.is_complete
                || is_protein_organelle(flags.genome)
                || flags.is_cross_kingdom
                || annots.has_prot
                || annots.has_cdregion
                || annots.has_gene
                || annots.has_biosrc
            {
                return Err(title_refusal(
                    "empty-title protein title generation with features/partial/organelle modifiers",
                ));
            }
        } else if flags.has_biosrc_modifiers
            || !matches!(flags.genome, 0 | 1)
            || annots.has_biosrc_clone
        {
            return Err(title_refusal(
                "empty-title BioSource modifier/organelle title generation",
            ));
        }
    }

    let defline =
        generate_bioseq_defline(&flags.raw_title, is_aa, false, &facts).map_err(|_| {
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

fn title_refusal(feature: &str) -> WriteDbError {
    WriteDbError::Unsupported(format!("{feature} is not supported by LOSAT's makeblastdb"))
}

#[derive(Default)]
struct DeflineFlags {
    raw_title: Vec<u8>,
    molinfo_tech_unsupported: bool,
    is_complete: bool,
    is_unverified: bool,
    is_unreviewed: bool,
    has_autodef: bool,
    is_pseudogene: bool,
    has_tpa_keyword: bool,
    has_unordered_kw: bool,
    has_source: bool,
    taxname: Option<Vec<u8>>,
    genome: i32,
    has_binomial: bool,
    has_metagenome_source: bool,
    is_multispecies: bool,
    is_cross_kingdom: bool,
    has_biosrc_modifiers: bool,
}

/// Scans Bioseq + ancestor descriptors in `CSeqdesc_CI` order, matching
/// `CDeflineGenerator::x_SetFlags` and `x_SetBioSrc`.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:840-1150,1264-1375
fn scan_defline_descriptors(
    bio_descr: &[Value],
    ancestors: &[&[Value]],
    is_pdb: bool,
    is_wp: bool,
) -> DeflineFlags {
    const F_MOLINFO: u8 = 1 << 0;
    const F_USER: u8 = 1 << 1;
    const F_SOURCE: u8 = 1 << 2;
    const F_GENBANK: u8 = 1 << 3;
    const F_EMBL: u8 = 1 << 4;
    const F_TITLE: u8 = 1 << 5;
    const F_COMMENT: u8 = 1 << 6;

    let mut needed = F_MOLINFO | F_USER | F_SOURCE | F_GENBANK | F_EMBL | F_TITLE | F_COMMENT;
    let mut flags = DeflineFlags {
        is_complete: true,
        ..DeflineFlags::default()
    };

    let s = schema();
    let mut keywords: Option<&[Value]> = None;
    let mut first_sk: Option<&[u8]> = None;
    let mut num_sk = 0usize;
    let mut sk_diff = false;

    for (level_idx, level) in std::iter::once(bio_descr)
        .chain(ancestors.iter().copied())
        .enumerate()
    {
        let on_bioseq = level_idx == 0;
        for (choice, inner) in level.iter().filter_map(as_choice) {
            if choice == s.seqdesc_user {
                if let Some(ty) = inner.member(s.user_obj_type).and_then(object_id_str) {
                    if ty == b"AutodefOptions" {
                        flags.has_autodef = true;
                    }
                    if needed & F_USER != 0 {
                        if ty.eq_ignore_ascii_case(b"Unverified") {
                            flags.is_unverified = true;
                            needed &= !F_USER;
                        } else if ty.eq_ignore_ascii_case(b"Unreviewed") {
                            flags.is_unreviewed = true;
                        }
                    }
                }
            } else if needed == 0 {
                continue;
            } else if choice == s.seqdesc_molinfo && (needed & F_MOLINFO != 0) {
                let tech = match inner.member(s.molinfo_tech) {
                    Some(Value::Enum(v)) => *v,
                    _ => 0,
                };
                let comp = match inner.member(s.molinfo_completeness) {
                    Some(Value::Enum(v)) => *v,
                    _ => 0,
                };
                if matches!(tech, 2 | 3 | 4 | 14 | 15 | 16 | 18 | 20 | 23 | 24) {
                    flags.molinfo_tech_unsupported = true;
                }
                if matches!(comp, 2..=5) {
                    flags.is_complete = false;
                }
                needed &= !F_MOLINFO;
            } else if choice == s.seqdesc_comment && (needed & F_COMMENT != 0) {
                if let Value::Str(comment) = inner {
                    if contains_bytes(comment, b"[CAUTION] Could be the product of a pseudogene") {
                        flags.is_pseudogene = true;
                    }
                }
            } else if choice == s.seqdesc_source {
                if needed & F_SOURCE != 0 {
                    needed &= !F_SOURCE;
                    flags.has_source = true;
                    inspect_primary_biosource(inner, &mut flags);
                }
                inspect_any_biosource(
                    inner,
                    is_wp,
                    &mut flags,
                    &mut first_sk,
                    &mut num_sk,
                    &mut sk_diff,
                );
            } else if choice == s.seqdesc_title && (needed & F_TITLE != 0) {
                if is_pdb || on_bioseq {
                    if let Value::Str(t) = inner {
                        flags.raw_title = t.clone();
                    }
                }
                needed &= !F_TITLE;
            } else if choice == s.seqdesc_genbank && (needed & F_GENBANK != 0) {
                if let Some(kw) = inner.member(s.gb_block_keywords).and_then(as_list) {
                    keywords = Some(kw);
                }
                needed &= !(F_GENBANK | F_EMBL);
            } else if choice == s.seqdesc_embl && (needed & F_EMBL != 0) {
                if let Some(kw) = inner.member(s.embl_block_keywords).and_then(as_list) {
                    keywords = Some(kw);
                }
                needed &= !F_EMBL;
            }
        }
    }

    if let Some(kw_list) = keywords {
        for kw in kw_list {
            let Value::Str(clause) = kw else {
                continue;
            };
            for part in clause.split(|&b| b == b';').filter(|p| !p.is_empty()) {
                if part.eq_ignore_ascii_case(b"UNORDERED") {
                    flags.has_unordered_kw = true;
                } else if part.eq_ignore_ascii_case(b"TPA:experimental")
                    || part.eq_ignore_ascii_case(b"TPA:inferential")
                    || part.eq_ignore_ascii_case(b"TPA:reassembly")
                    || part.eq_ignore_ascii_case(b"TPA:assembly")
                {
                    flags.has_tpa_keyword = true;
                }
            }
        }
    }

    flags
}

fn inspect_primary_biosource(bsrc: &Value, flags: &mut DeflineFlags) {
    let s = schema();
    if let Some(Value::Enum(g)) = bsrc.member(s.biosrc_genome) {
        flags.genome = *g;
    }
    if let Some(subs) = bsrc.member(s.biosrc_subtype).and_then(as_list) {
        for sbs in subs {
            let Some(Value::Enum(sub)) = sbs.member(s.subsrc_subtype) else {
                continue;
            };
            let has_name =
                matches!(sbs.member(s.subsrc_name), Some(Value::Str(v)) if !v.is_empty());
            if has_name && matches!(*sub, 1 | 2 | 3 | 19 | 24 | 39) {
                flags.has_biosrc_modifiers = true;
            }
        }
    }
    let Some(org) = bsrc.member(s.biosrc_org) else {
        return;
    };
    if let Some(Value::Str(taxname)) = org.member(s.org_ref_taxname) {
        flags.taxname = Some(taxname.clone());
    }
    let Some(orgname) = org.member(s.org_ref_orgname) else {
        return;
    };
    if let Some((choice, inner)) = orgname.member(s.orgname_name).and_then(as_choice) {
        if choice == s.orgname_binomial {
            let has_g =
                matches!(inner.member(s.binomial_genus), Some(Value::Str(v)) if !v.is_empty());
            let has_sp =
                matches!(inner.member(s.binomial_species), Some(Value::Str(v)) if !v.is_empty());
            if has_g || has_sp {
                flags.has_binomial = true;
            }
        } else if choice == s.orgname_partial {
            if let Some(elems) = as_list(inner) {
                for te in elems {
                    let Some(Value::Enum(fl)) = te.member(s.taxelem_fixed_level) else {
                        continue;
                    };
                    if *fl > 0 {
                        flags.is_multispecies = true;
                    } else if let Some(Value::Str(lvl)) = te.member(s.taxelem_level) {
                        if !lvl.eq_ignore_ascii_case(b"species") {
                            flags.is_multispecies = true;
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
            let has_val =
                matches!(m.member(s.orgmod_subname), Some(Value::Str(v)) if !v.is_empty());
            if has_val && matches!(*sub, 2 | 3 | 10 | 17 | 23 | 31) {
                flags.has_biosrc_modifiers = true;
            }
        }
    }
}

fn inspect_any_biosource<'a>(
    bsrc: &'a Value,
    is_wp: bool,
    flags: &mut DeflineFlags,
    first_sk: &mut Option<&'a [u8]>,
    num_sk: &mut usize,
    sk_diff: &mut bool,
) {
    let s = schema();
    let Some(orgname) = bsrc
        .member(s.biosrc_org)
        .and_then(|o| o.member(s.org_ref_orgname))
    else {
        return;
    };
    if let Some(om_list) = orgname.member(s.orgname_mod).and_then(as_list) {
        for m in om_list {
            let is_mg = matches!(m.member(s.orgmod_subtype), Some(Value::Enum(37)));
            if is_mg && matches!(m.member(s.orgmod_subname), Some(Value::Str(v)) if !v.is_empty()) {
                flags.has_metagenome_source = true;
            }
        }
    }
    if is_wp {
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
                            match *first_sk {
                                None => *first_sk = Some(name.as_slice()),
                                Some(first) if !first.eq_ignore_ascii_case(name) => *sk_diff = true,
                                _ => {}
                            }
                        }
                        if *num_sk > 1 && *sk_diff {
                            flags.is_cross_kingdom = true;
                        }
                    }
                }
            }
        }
    }
}

/// Whether `genome` has a non-empty entry in `s_proteinOrganellePrefix` (create_defline.cpp:2123).
fn is_protein_organelle(genome: i32) -> bool {
    matches!(
        genome,
        2 | 3 | 4 | 5 | 6 | 7 | 9 | 12 | 15 | 16 | 17 | 18 | 19 | 20 | 22
    )
}

fn trim_end_slice<'a>(bytes: &'a [u8], chars: &[u8]) -> &'a [u8] {
    match bytes.iter().rposition(|b| !chars.contains(b)) {
        Some(pos) => &bytes[..=pos],
        None => bytes,
    }
}

fn contains_bytes(haystack: &[u8], needle: &[u8]) -> bool {
    haystack.windows(needle.len()).any(|w| w == needle)
}

#[cfg(test)]
#[path = "asn_source_tests.rs"]
mod tests;
