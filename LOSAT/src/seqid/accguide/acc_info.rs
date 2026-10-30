//! `CSeq_id::EAccessionInfo`: the values and the names that the accession guide uses
//! (`sc_AccInfoArray`).
//!
//! NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_id.hpp:246-260,262-303,305-318
//! ```c++
//!         // Mask for Seq_id type; allow 8 bits to be safe
//!         eAcc_type_mask = 0xff,
//!
//!         // Useful general flags (not inherent in choice of division,
//!         // albeit not necessarily applicable to all divisions).
//!         fAcc_nuc       = 0x80000000,
//!         fAcc_prot      = 0x40000000,
//!         fAcc_seq       = 0,
//!         fAcc_predicted = 0x20000000, // only for refseq
//!         fAcc_specials  = 0x10000000, // has special cases; only used internally
//!         fAcc_master    = 0x08000000,
//!         fAcc_ncbo      = 0x04000000, // for refseq pathogen detection pipeline
//!         fAcc_fallback  = 0x02000000, // is a fallback; only used internally
//!         fAcc_vdb_only  = 0x01000000,
//!         eAcc_flag_mask = 0xff000000,
//!
//!         // Divisions and categories (multiples of 1 << 8; always
//!         // globally unique nowadays, no matter how specialized)
//!         eAcc_other         =   0 << 8, // no further classification
//!         eAcc_est           =   1 << 8, // expressed sequence tag (mRNA)
//!         ...
//!         eAcc_targeted      = 139 << 8, // targeted genomic project
//!         eAcc_division_mask = 0xff00,
//!
//!         // Internal values combinining a division and a flag for convenience.
//!         eAcc_wgs_master          = eAcc_wgs          | fAcc_master,
//!         ...
//!         eAcc_targeted_vdb_master = eAcc_targeted | fAcc_master | fAcc_vdb_only,
//!
//!         // Internal macro, left defined only temporarily
//! #define NCBI_ACC(type, div, mol) eSeqId_##type | eAcc_##div | fAcc_##mol
//! ```
//! The values of the table below are the enumerators of Seq_id.hpp:323-541 written with
//! that macro (`acc`), each under its name in `sc_AccInfoArray` (Seq_id.cpp:646-860).

// NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_id_.hpp (generated from
// seqloc.asn:24-45): `E_Choice` = 0 for `e_not_set`, then 1 + the position of the
// alternative (`e_Local` 1 ... `e_Named_annot_track` 20); `eSeqId_<type>` = `e_<Type>`
// (Seq_id.hpp:221-244, `eSeqId_refseq` = `eSeqId_other`).
pub const SEQID_NOT_SET: u32 = 0;
pub const SEQID_LOCAL: u32 = 1;
pub const SEQID_GIBBSQ: u32 = 2;
pub const SEQID_GIBBMT: u32 = 3;
pub const SEQID_GIIM: u32 = 4;
pub const SEQID_GENBANK: u32 = 5;
pub const SEQID_EMBL: u32 = 6;
pub const SEQID_PIR: u32 = 7;
pub const SEQID_SWISSPROT: u32 = 8;
pub const SEQID_PATENT: u32 = 9;
pub const SEQID_OTHER: u32 = 10;
pub const SEQID_REFSEQ: u32 = SEQID_OTHER;
pub const SEQID_GENERAL: u32 = 11;
pub const SEQID_GI: u32 = 12;
pub const SEQID_DDBJ: u32 = 13;
pub const SEQID_PRF: u32 = 14;
pub const SEQID_PDB: u32 = 15;
pub const SEQID_TPG: u32 = 16;
pub const SEQID_TPE: u32 = 17;
pub const SEQID_TPD: u32 = 18;
pub const SEQID_GPIPE: u32 = 19;
pub const SEQID_NAMED_ANNOT_TRACK: u32 = 20;

pub const TYPE_MASK: u32 = 0xff;
pub const F_NUC: u32 = 0x8000_0000;
pub const F_PROT: u32 = 0x4000_0000;
pub const F_SEQ: u32 = 0;
pub const F_PREDICTED: u32 = 0x2000_0000;
pub const F_SPECIALS: u32 = 0x1000_0000;
pub const F_MASTER: u32 = 0x0800_0000;
pub const F_NCBO: u32 = 0x0400_0000;
pub const F_FALLBACK: u32 = 0x0200_0000;
pub const F_VDB_ONLY: u32 = 0x0100_0000;
pub const FLAG_MASK: u32 = 0xff00_0000;

pub const OTHER: u32 = 0 << 8;
pub const EST: u32 = 1 << 8;
pub const DIRSUB: u32 = 2 << 8;
pub const DIV_PATENT: u32 = 3 << 8;
pub const MRNA: u32 = 4 << 8;
pub const NCRNA: u32 = 5 << 8;
pub const GSDB_DS: u32 = 6 << 8;
pub const GSDB: u32 = 7 << 8;
pub const BACKBONE: u32 = 8 << 8;
pub const TSA: u32 = 9 << 8;
pub const SEGSET: u32 = 10 << 8;
pub const GSS: u32 = 11 << 8;
pub const UNIQUE: u32 = 12 << 8;
pub const AMBIG: u32 = 112 << 8;
pub const AMBIG_G: u32 = 113 << 8;
pub const AMBIG_E: u32 = 114 << 8;
pub const AMBIG_GE: u32 = 115 << 8;
pub const AMBIG_D: u32 = 116 << 8;
pub const AMBIG_GD: u32 = 117 << 8;
pub const AMBIG_ED: u32 = 118 << 8;
pub const AMBIG_GED: u32 = 119 << 8;
pub const UNRESERVED: u32 = 127 << 8;
pub const GENOME: u32 = 128 << 8;
pub const HTGS: u32 = 129 << 8;
pub const CON: u32 = 130 << 8;
pub const WGS: u32 = 131 << 8;
pub const CHROMOSOME: u32 = 133 << 8;
pub const GENOMIC_RGN: u32 = 134 << 8;
pub const WGS_INTERMED: u32 = 135 << 8;
pub const STS: u32 = 136 << 8;
pub const MGA: u32 = 137 << 8;
pub const OPTICAL_MAP: u32 = 138 << 8;
pub const TARGETED: u32 = 139 << 8;
pub const DIVISION_MASK: u32 = 0xff00;

pub const WGS_MASTER: u32 = WGS | F_MASTER;
pub const WGS_INTERMED_MASTER: u32 = WGS_INTERMED | F_MASTER;
pub const TSA_MASTER: u32 = TSA | F_MASTER;
pub const TARGETED_MASTER: u32 = TARGETED | F_MASTER;
pub const WGS_VDB_ONLY: u32 = WGS | F_VDB_ONLY;
pub const WGS_INTERMED_VDB_ONLY: u32 = WGS_INTERMED | F_VDB_ONLY;
pub const TSA_VDB_ONLY: u32 = TSA | F_VDB_ONLY;
pub const TARGETED_VDB_ONLY: u32 = TARGETED | F_VDB_ONLY;
pub const WGS_VDB_MASTER: u32 = WGS | F_MASTER | F_VDB_ONLY;
pub const WGS_INTERMED_VDB_MASTER: u32 = WGS_INTERMED | F_MASTER | F_VDB_ONLY;
pub const TSA_VDB_MASTER: u32 = TSA | F_MASTER | F_VDB_ONLY;
pub const TARGETED_VDB_MASTER: u32 = TARGETED | F_MASTER | F_VDB_ONLY;

pub(super) const NUC: u32 = F_NUC;
pub(super) const PROT: u32 = F_PROT;
pub(super) const SEQ: u32 = F_SEQ;

/// `NCBI_ACC(type, div, mol)`.
pub(super) const fn acc(seq_id: u32, division: u32, mol: u32) -> u32 {
    seq_id | division | mol
}

/// `eAcc_unknown`.
pub const UNKNOWN: u32 = acc(SEQID_NOT_SET, OTHER, SEQ);
/// `eAcc_general`, `eAcc_general_nuc`.
pub const GENERAL: u32 = acc(SEQID_GENERAL, OTHER, SEQ);
pub const GENERAL_NUC: u32 = acc(SEQID_GENERAL, OTHER, NUC);
/// `eAcc_gi`.
pub const GI: u32 = acc(SEQID_GI, OTHER, SEQ);

/// `sc_AccInfoArray`: the names of the values, in its order (case-insensitive, for a
/// binary search).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:645-651,856-864
/// ```c++
/// typedef SStaticPair<const char*, CSeq_id::EAccessionInfo> TAccInfoMapEntry;
/// // used for binary searching; must be in order.
/// static const TAccInfoMapEntry sc_AccInfoArray[] = {
///     { "ambiguous_nuc",           CSeq_id::eAcc_ambiguous_nuc },
///     { "ddbj_con",                CSeq_id::eAcc_ddbj_con },
///     ...
///     { "unreserved_prot",         CSeq_id::eAcc_unreserved_prot }
/// };
/// typedef CStaticPairArrayMap<CTempString, CSeq_id::EAccessionInfo,
///                             PNocase_Generic<CTempString> > TAccInfoMap;
/// DEFINE_STATIC_ARRAY_MAP_WITH_COPY(TAccInfoMap, sc_AccInfoMap, sc_AccInfoArray);
/// ```
pub(super) const ACC_INFO_NAMES: [(&str, u32); 213] = [
    ("ambiguous_nuc", acc(SEQID_NOT_SET, AMBIG, NUC)),
    ("ddbj_con", acc(SEQID_DDBJ, CON, NUC)),
    ("ddbj_dirsub", acc(SEQID_DDBJ, DIRSUB, NUC)),
    ("ddbj_est", acc(SEQID_DDBJ, EST, NUC)),
    ("ddbj_genome", acc(SEQID_DDBJ, GENOME, NUC)),
    ("ddbj_gss", acc(SEQID_DDBJ, GSS, NUC)),
    ("ddbj_htgs", acc(SEQID_DDBJ, HTGS, NUC)),
    ("ddbj_mga", acc(SEQID_DDBJ, MGA, NUC)),
    ("ddbj_mrna", acc(SEQID_DDBJ, MRNA, NUC)),
    ("ddbj_other", acc(SEQID_DDBJ, OTHER, SEQ)),
    ("ddbj_other_nuc", acc(SEQID_DDBJ, OTHER, NUC)),
    ("ddbj_patent", acc(SEQID_DDBJ, DIV_PATENT, NUC)),
    ("ddbj_prot", acc(SEQID_DDBJ, OTHER, PROT)),
    ("ddbj_targeted_nuc", acc(SEQID_DDBJ, TARGETED, NUC)),
    ("ddbj_targetedm_nuc", acc(SEQID_DDBJ, TARGETED_MASTER, NUC)),
    (
        "ddbj_targetedv_nuc",
        acc(SEQID_DDBJ, TARGETED_VDB_ONLY, NUC),
    ),
    (
        "ddbj_targetedvm_nuc",
        acc(SEQID_DDBJ, TARGETED_VDB_MASTER, NUC),
    ),
    ("ddbj_tpa_chromosome", acc(SEQID_TPD, CHROMOSOME, NUC)),
    ("ddbj_tpa_con", acc(SEQID_TPD, CON, NUC)),
    ("ddbj_tpa_nuc", acc(SEQID_TPD, OTHER, NUC)),
    ("ddbj_tpa_other", acc(SEQID_TPD, OTHER, SEQ)),
    ("ddbj_tpa_prot", acc(SEQID_TPD, OTHER, PROT)),
    ("ddbj_tpa_targeted_nuc", acc(SEQID_TPD, TARGETED, NUC)),
    (
        "ddbj_tpa_targetedm_nuc",
        acc(SEQID_TPD, TARGETED_MASTER, NUC),
    ),
    (
        "ddbj_tpa_targetedv_nuc",
        acc(SEQID_TPD, TARGETED_VDB_ONLY, NUC),
    ),
    (
        "ddbj_tpa_targetedvm_nuc",
        acc(SEQID_TPD, TARGETED_VDB_MASTER, NUC),
    ),
    ("ddbj_tpa_tsa_nuc", acc(SEQID_TPD, TSA, NUC)),
    ("ddbj_tpa_tsa_prot", acc(SEQID_TPD, TSA, PROT)),
    ("ddbj_tpa_tsam_nuc", acc(SEQID_TPD, TSA_MASTER, NUC)),
    ("ddbj_tpa_tsam_prot", acc(SEQID_TPD, TSA_MASTER, PROT)),
    ("ddbj_tpa_tsav_nuc", acc(SEQID_TPD, TSA_VDB_ONLY, NUC)),
    ("ddbj_tpa_tsav_prot", acc(SEQID_TPD, TSA_VDB_ONLY, PROT)),
    ("ddbj_tpa_tsavm_nuc", acc(SEQID_TPD, TSA_VDB_MASTER, NUC)),
    ("ddbj_tpa_tsavm_prot", acc(SEQID_TPD, TSA_VDB_MASTER, PROT)),
    ("ddbj_tpa_wgs_nuc", acc(SEQID_TPD, WGS, NUC)),
    ("ddbj_tpa_wgs_prot", acc(SEQID_TPD, WGS, PROT)),
    ("ddbj_tpa_wgsm_nuc", acc(SEQID_TPD, WGS_MASTER, NUC)),
    ("ddbj_tpa_wgsm_prot", acc(SEQID_TPD, WGS_MASTER, PROT)),
    ("ddbj_tpa_wgsv_nuc", acc(SEQID_TPD, WGS_VDB_ONLY, NUC)),
    ("ddbj_tpa_wgsv_prot", acc(SEQID_TPD, WGS_VDB_ONLY, PROT)),
    ("ddbj_tpa_wgsvm_nuc", acc(SEQID_TPD, WGS_VDB_MASTER, NUC)),
    ("ddbj_tpa_wgsvm_prot", acc(SEQID_TPD, WGS_VDB_MASTER, PROT)),
    ("ddbj_tsa_nuc", acc(SEQID_DDBJ, TSA, NUC)),
    ("ddbj_tsa_prot", acc(SEQID_DDBJ, TSA, PROT)),
    ("ddbj_tsam_nuc", acc(SEQID_DDBJ, TSA_MASTER, NUC)),
    ("ddbj_tsam_prot", acc(SEQID_DDBJ, TSA_MASTER, PROT)),
    ("ddbj_tsav_nuc", acc(SEQID_DDBJ, TSA_VDB_ONLY, NUC)),
    ("ddbj_tsav_prot", acc(SEQID_DDBJ, TSA_VDB_ONLY, PROT)),
    ("ddbj_tsavm_nuc", acc(SEQID_DDBJ, TSA_VDB_MASTER, NUC)),
    ("ddbj_tsavm_prot", acc(SEQID_DDBJ, TSA_VDB_MASTER, PROT)),
    ("ddbj_wgs_nuc", acc(SEQID_DDBJ, WGS, NUC)),
    ("ddbj_wgs_prot", acc(SEQID_DDBJ, WGS, PROT)),
    ("ddbj_wgsm_nuc", acc(SEQID_DDBJ, WGS_MASTER, NUC)),
    ("ddbj_wgsm_prot", acc(SEQID_DDBJ, WGS_MASTER, PROT)),
    ("ddbj_wgsv_nuc", acc(SEQID_DDBJ, WGS_VDB_ONLY, NUC)),
    ("ddbj_wgsv_prot", acc(SEQID_DDBJ, WGS_VDB_ONLY, PROT)),
    ("ddbj_wgsvm_nuc", acc(SEQID_DDBJ, WGS_VDB_MASTER, NUC)),
    ("ddbj_wgsvm_prot", acc(SEQID_DDBJ, WGS_VDB_MASTER, PROT)),
    ("embl_con", acc(SEQID_EMBL, CON, NUC)),
    ("embl_ddbj", acc(SEQID_NOT_SET, AMBIG_ED, NUC)),
    ("embl_dirsub", acc(SEQID_EMBL, DIRSUB, NUC)),
    ("embl_est", acc(SEQID_EMBL, EST, NUC)),
    ("embl_genome", acc(SEQID_EMBL, GENOME, NUC)),
    ("embl_gss", acc(SEQID_EMBL, GSS, NUC)),
    ("embl_htgs", acc(SEQID_EMBL, HTGS, NUC)),
    ("embl_mga", acc(SEQID_EMBL, MGA, NUC)),
    ("embl_other", acc(SEQID_EMBL, OTHER, SEQ)),
    ("embl_other_nuc", acc(SEQID_EMBL, OTHER, NUC)),
    ("embl_patent", acc(SEQID_EMBL, DIV_PATENT, NUC)),
    ("embl_prot", acc(SEQID_EMBL, OTHER, PROT)),
    ("embl_tpa_nuc", acc(SEQID_TPE, OTHER, NUC)),
    ("embl_tpa_other", acc(SEQID_TPE, OTHER, SEQ)),
    ("embl_tpa_prot", acc(SEQID_TPE, OTHER, PROT)),
    ("embl_tpa_tsa_nuc", acc(SEQID_TPE, TSA, NUC)),
    ("embl_tpa_tsa_prot", acc(SEQID_TPE, TSA, PROT)),
    ("embl_tpa_tsam_nuc", acc(SEQID_TPE, TSA_MASTER, NUC)),
    ("embl_tpa_tsam_prot", acc(SEQID_TPE, TSA_MASTER, PROT)),
    ("embl_tpa_tsav_nuc", acc(SEQID_TPE, TSA_VDB_ONLY, NUC)),
    ("embl_tpa_tsav_prot", acc(SEQID_TPE, TSA_VDB_ONLY, PROT)),
    ("embl_tpa_tsavm_nuc", acc(SEQID_TPE, TSA_VDB_MASTER, NUC)),
    ("embl_tpa_tsavm_prot", acc(SEQID_TPE, TSA_VDB_MASTER, PROT)),
    ("embl_tpa_wgs_nuc", acc(SEQID_TPE, WGS, NUC)),
    ("embl_tpa_wgs_prot", acc(SEQID_TPE, WGS, PROT)),
    ("embl_tpa_wgsm_nuc", acc(SEQID_TPE, WGS_MASTER, NUC)),
    ("embl_tpa_wgsm_prot", acc(SEQID_TPE, WGS_MASTER, PROT)),
    ("embl_tpa_wgsv_nuc", acc(SEQID_TPE, WGS_VDB_ONLY, NUC)),
    ("embl_tpa_wgsv_prot", acc(SEQID_TPE, WGS_VDB_ONLY, PROT)),
    ("embl_tpa_wgsvm_nuc", acc(SEQID_TPE, WGS_VDB_MASTER, NUC)),
    ("embl_tpa_wgsvm_prot", acc(SEQID_TPE, WGS_VDB_MASTER, PROT)),
    ("embl_tsa_nuc", acc(SEQID_EMBL, TSA, NUC)),
    ("embl_tsa_prot", acc(SEQID_EMBL, TSA, PROT)),
    ("embl_tsam_nuc", acc(SEQID_EMBL, TSA_MASTER, NUC)),
    ("embl_tsam_prot", acc(SEQID_EMBL, TSA_MASTER, PROT)),
    ("embl_tsav_nuc", acc(SEQID_EMBL, TSA_VDB_ONLY, NUC)),
    ("embl_tsav_prot", acc(SEQID_EMBL, TSA_VDB_ONLY, PROT)),
    ("embl_tsavm_nuc", acc(SEQID_EMBL, TSA_VDB_MASTER, NUC)),
    ("embl_tsavm_prot", acc(SEQID_EMBL, TSA_VDB_MASTER, PROT)),
    ("embl_wgs_nuc", acc(SEQID_EMBL, WGS, NUC)),
    ("embl_wgs_prot", acc(SEQID_EMBL, WGS, PROT)),
    ("embl_wgsm_nuc", acc(SEQID_EMBL, WGS_MASTER, NUC)),
    ("embl_wgsm_prot", acc(SEQID_EMBL, WGS_MASTER, PROT)),
    ("embl_wgsv_nuc", acc(SEQID_EMBL, WGS_VDB_ONLY, NUC)),
    ("embl_wgsv_prot", acc(SEQID_EMBL, WGS_VDB_ONLY, PROT)),
    ("embl_wgsvm_nuc", acc(SEQID_EMBL, WGS_VDB_MASTER, NUC)),
    ("embl_wgsvm_prot", acc(SEQID_EMBL, WGS_VDB_MASTER, PROT)),
    ("gb_backbone", acc(SEQID_GENBANK, BACKBONE, NUC)),
    ("gb_cdna", acc(SEQID_GENBANK, MRNA, NUC)),
    ("gb_chromosome", acc(SEQID_GENBANK, CHROMOSOME, NUC)),
    ("gb_con", acc(SEQID_GENBANK, CON, NUC)),
    ("gb_ddbj", acc(SEQID_NOT_SET, AMBIG_GD, NUC)),
    ("gb_dirsub", acc(SEQID_GENBANK, DIRSUB, NUC)),
    ("gb_embl", acc(SEQID_NOT_SET, AMBIG_GE, NUC)),
    ("gb_embl_ddbj", acc(SEQID_NOT_SET, AMBIG_GED, NUC)),
    ("gb_est", acc(SEQID_GENBANK, EST, NUC)),
    ("gb_genome", acc(SEQID_GENBANK, GENOME, NUC)),
    ("gb_gsdb", acc(SEQID_GENBANK, GSDB, NUC)),
    ("gb_gss", acc(SEQID_GENBANK, GSS, NUC)),
    ("gb_htgs", acc(SEQID_GENBANK, HTGS, NUC)),
    ("gb_mga", acc(SEQID_GENBANK, MGA, NUC)),
    ("gb_optical_map", acc(SEQID_GENBANK, OPTICAL_MAP, NUC)),
    ("gb_other", acc(SEQID_GENBANK, OTHER, SEQ)),
    ("gb_other_nuc", acc(SEQID_GENBANK, OTHER, NUC)),
    ("gb_patent", acc(SEQID_GENBANK, DIV_PATENT, SEQ)),
    ("gb_patent_prot", acc(SEQID_GENBANK, DIV_PATENT, PROT)),
    ("gb_prot", acc(SEQID_GENBANK, OTHER, PROT)),
    ("gb_segset", acc(SEQID_GENBANK, SEGSET, SEQ)),
    ("gb_sts", acc(SEQID_GENBANK, STS, NUC)),
    ("gb_targeted_nuc", acc(SEQID_GENBANK, TARGETED, NUC)),
    ("gb_tpa_chromosome", acc(SEQID_TPG, CHROMOSOME, NUC)),
    ("gb_tpa_con", acc(SEQID_TPG, CON, NUC)),
    ("gb_tpa_nuc", acc(SEQID_TPG, OTHER, NUC)),
    ("gb_tpa_other", acc(SEQID_TPG, OTHER, SEQ)),
    ("gb_tpa_prot", acc(SEQID_TPG, OTHER, PROT)),
    ("gb_tpa_segset", acc(SEQID_TPG, SEGSET, NUC)),
    ("gb_tpa_wgs_nuc", acc(SEQID_TPG, WGS, NUC)),
    ("gb_tpa_wgs_prot", acc(SEQID_TPG, WGS, PROT)),
    ("gb_tpa_wgsm_nuc", acc(SEQID_TPG, WGS_MASTER, NUC)),
    ("gb_tpa_wgsm_prot", acc(SEQID_TPG, WGS_MASTER, PROT)),
    ("gb_tpa_wgsv_nuc", acc(SEQID_TPG, WGS_VDB_ONLY, NUC)),
    ("gb_tpa_wgsv_prot", acc(SEQID_TPG, WGS_VDB_ONLY, PROT)),
    ("gb_tpa_wgsvm_nuc", acc(SEQID_TPG, WGS_VDB_MASTER, NUC)),
    ("gb_tpa_wgsvm_prot", acc(SEQID_TPG, WGS_VDB_MASTER, PROT)),
    ("gb_tsa_nuc", acc(SEQID_GENBANK, TSA, NUC)),
    ("gb_tsa_prot", acc(SEQID_GENBANK, TSA, PROT)),
    ("gb_tsam_nuc", acc(SEQID_GENBANK, TSA_MASTER, NUC)),
    ("gb_tsam_prot", acc(SEQID_GENBANK, TSA_MASTER, PROT)),
    ("gb_tsav_nuc", acc(SEQID_GENBANK, TSA_VDB_ONLY, NUC)),
    ("gb_tsav_prot", acc(SEQID_GENBANK, TSA_VDB_ONLY, PROT)),
    ("gb_tsavm_nuc", acc(SEQID_GENBANK, TSA_VDB_MASTER, NUC)),
    ("gb_tsavm_prot", acc(SEQID_GENBANK, TSA_VDB_MASTER, PROT)),
    ("gb_wgs_nuc", acc(SEQID_GENBANK, WGS, NUC)),
    ("gb_wgs_prot", acc(SEQID_GENBANK, WGS, PROT)),
    ("gb_wgsm_nuc", acc(SEQID_GENBANK, WGS_MASTER, NUC)),
    ("gb_wgsm_prot", acc(SEQID_GENBANK, WGS_MASTER, PROT)),
    ("gb_wgsv_nuc", acc(SEQID_GENBANK, WGS_VDB_ONLY, NUC)),
    ("gb_wgsv_prot", acc(SEQID_GENBANK, WGS_VDB_ONLY, PROT)),
    ("gb_wgsvm_nuc", acc(SEQID_GENBANK, WGS_VDB_MASTER, NUC)),
    ("gb_wgsvm_prot", acc(SEQID_GENBANK, WGS_VDB_MASTER, PROT)),
    ("general", acc(SEQID_GENERAL, OTHER, SEQ)),
    ("general_nuc", acc(SEQID_GENERAL, OTHER, NUC)),
    ("general_prot", acc(SEQID_GENERAL, OTHER, PROT)),
    ("gi", acc(SEQID_GI, OTHER, SEQ)),
    ("gibbmt", acc(SEQID_GIBBMT, OTHER, SEQ)),
    ("gibbsq", acc(SEQID_GIBBSQ, OTHER, SEQ)),
    ("giim", acc(SEQID_GIIM, OTHER, SEQ)),
    ("gpipe_chromosome", acc(SEQID_GPIPE, CHROMOSOME, NUC)),
    ("gpipe_genomic", acc(SEQID_GPIPE, GENOMIC_RGN, NUC)),
    ("gpipe_mrna", acc(SEQID_GPIPE, MRNA, NUC)),
    ("gpipe_ncrna", acc(SEQID_GPIPE, NCRNA, NUC)),
    ("gpipe_other_nuc", acc(SEQID_GPIPE, OTHER, NUC)),
    ("gpipe_prot", acc(SEQID_GPIPE, OTHER, PROT)),
    ("gpipe_scaffold", acc(SEQID_GPIPE, CON, NUC)),
    ("gpipe_unreserved", acc(SEQID_GPIPE, UNRESERVED, SEQ)),
    ("gsdb_dirsub", acc(SEQID_GENBANK, GSDB_DS, NUC)),
    ("local", acc(SEQID_LOCAL, OTHER, SEQ)),
    ("maybe_ddbj", acc(SEQID_NOT_SET, AMBIG_D, NUC)),
    ("maybe_embl", acc(SEQID_NOT_SET, AMBIG_E, NUC)),
    ("maybe_gb", acc(SEQID_NOT_SET, AMBIG_G, NUC)),
    (
        "named_annot_track",
        acc(SEQID_NAMED_ANNOT_TRACK, OTHER, SEQ),
    ),
    ("patent", acc(SEQID_PATENT, DIV_PATENT, SEQ)),
    ("pdb", acc(SEQID_PDB, OTHER, SEQ)),
    ("pir", acc(SEQID_PIR, OTHER, PROT)),
    ("prf", acc(SEQID_PRF, OTHER, PROT)),
    ("refseq_chromosome", acc(SEQID_REFSEQ, CHROMOSOME, NUC)),
    (
        "refseq_chromosome_ncbo",
        acc(SEQID_REFSEQ, CHROMOSOME, NUC) | F_NCBO,
    ),
    ("refseq_contig", acc(SEQID_REFSEQ, CON, NUC)),
    ("refseq_contig_ncbo", acc(SEQID_REFSEQ, CON, NUC) | F_NCBO),
    ("refseq_genome", acc(SEQID_REFSEQ, GENOME, NUC)),
    ("refseq_genomic", acc(SEQID_REFSEQ, GENOMIC_RGN, NUC)),
    ("refseq_mrna", acc(SEQID_REFSEQ, MRNA, NUC)),
    (
        "refseq_mrna_predicted",
        acc(SEQID_REFSEQ, MRNA, NUC) | F_PREDICTED,
    ),
    ("refseq_ncrna", acc(SEQID_REFSEQ, NCRNA, NUC)),
    (
        "refseq_ncrna_predicted",
        acc(SEQID_REFSEQ, NCRNA, NUC) | F_PREDICTED,
    ),
    ("refseq_prot", acc(SEQID_REFSEQ, OTHER, PROT)),
    (
        "refseq_prot_predicted",
        acc(SEQID_REFSEQ, OTHER, PROT) | F_PREDICTED,
    ),
    ("refseq_unique_prot", acc(SEQID_REFSEQ, UNIQUE, PROT)),
    ("refseq_unreserved", acc(SEQID_REFSEQ, UNRESERVED, SEQ)),
    ("refseq_wgs_intermed", acc(SEQID_REFSEQ, WGS_INTERMED, NUC)),
    ("refseq_wgs_nuc", acc(SEQID_REFSEQ, WGS, NUC)),
    ("refseq_wgs_prot", acc(SEQID_REFSEQ, WGS, PROT)),
    (
        "refseq_wgsm_intermed",
        acc(SEQID_REFSEQ, WGS_INTERMED_MASTER, NUC),
    ),
    ("refseq_wgsm_nuc", acc(SEQID_REFSEQ, WGS_MASTER, NUC)),
    ("refseq_wgsm_prot", acc(SEQID_REFSEQ, WGS_MASTER, PROT)),
    (
        "refseq_wgsv_intermed",
        acc(SEQID_REFSEQ, WGS_INTERMED_VDB_ONLY, NUC),
    ),
    ("refseq_wgsv_nuc", acc(SEQID_REFSEQ, WGS_VDB_ONLY, NUC)),
    ("refseq_wgsv_prot", acc(SEQID_REFSEQ, WGS_VDB_ONLY, PROT)),
    (
        "refseq_wgsvm_intermed",
        acc(SEQID_REFSEQ, WGS_INTERMED_VDB_MASTER, NUC),
    ),
    ("refseq_wgsvm_nuc", acc(SEQID_REFSEQ, WGS_VDB_MASTER, NUC)),
    ("refseq_wgsvm_prot", acc(SEQID_REFSEQ, WGS_VDB_MASTER, PROT)),
    ("swissprot", acc(SEQID_SWISSPROT, OTHER, PROT)),
    ("unknown", acc(SEQID_NOT_SET, OTHER, SEQ)),
    ("unreserved_nuc", acc(SEQID_NOT_SET, UNRESERVED, NUC)),
    ("unreserved_prot", acc(SEQID_NOT_SET, UNRESERVED, PROT)),
];

/// `sc_AccInfoMap.find(name)`: the value of a type name, compared without case.
pub(super) fn find_acc_info(name: &[u8]) -> Option<u32> {
    let lower = name.to_ascii_lowercase();
    ACC_INFO_NAMES
        .binary_search_by(|(known, _)| known.as_bytes().cmp(lower.as_slice()))
        .ok()
        .map(|index| ACC_INFO_NAMES[index].1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_names_are_in_order_and_the_values_are_ncbi_s() {
        assert!(ACC_INFO_NAMES.windows(2).all(|pair| pair[0].0 < pair[1].0));
        // Values of the enumerators (NCBI_ACC of Seq_id.hpp).
        assert_eq!(find_acc_info(b"gb_prot"), Some(0x4000_0005));
        assert_eq!(find_acc_info(b"REFSEQ_PROT_PREDICTED"), Some(0x6000_000a));
        assert_eq!(find_acc_info(b"refseq_chromosome"), Some(0x8000_850a));
        assert_eq!(find_acc_info(b"gb_wgsvm_nuc"), Some(0x8900_8305));
        assert_eq!(find_acc_info(b"pdb"), Some(15));
        assert_eq!(find_acc_info(b"unknown"), Some(UNKNOWN));
        assert_eq!(find_acc_info(b"no_such_type"), None);
    }
}
