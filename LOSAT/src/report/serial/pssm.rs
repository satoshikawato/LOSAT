//! NCBI-ScoreMat's `PssmWithParameters` (objects/scoremat/scoremat.asn): the value that psiblast
//! writes with `-out_pssm` (`MSerial_AsnText`), and its writer over an `Encoder`, in the member
//! order of the ASN.1 module, so that the ASN.1 binary and the archive can write the same value.
//! The reader of the ASN.1 text is `blastinput::pssm_text`.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:385-396
//! ```c++
//! CPsiBlastApp::SavePssmToFile(CRef<CPssmWithParameters> pssm,
//!                              const CPsiBlastIterationState* itr)
//! {
//!     if (pssm.Empty()) {
//!         return;
//!     }
//!
//!     if (m_CmdLineArgs->SaveCheckpoint() &&
//!         (itr == NULL || // this is true in the case of remote PSI-BLAST
//!          itr->GetIterationNumber() >= 1)) {
//!         *m_CmdLineArgs->GetCheckpointStream() << MSerial_AsnText << *pssm;
//!     }
//! ```
//!
//! A member is written when it is set: an OPTIONAL member that is `None` is not, and neither is a
//! member with a DEFAULT that was not set (`is_protein` and `by_row` hold whether and how they were
//! set). The members that LOSAT neither writes nor reads (`constraints` and `annotatedSites` of
//! the parameters, the descriptors of the query other than its titles, the members of the query's
//! `Seq-inst` after `seq-data` and before it other than `repr`, `mol` and `length`) are not kept:
//! the reader skips them.

use super::encoder::Encoder;

/// `PssmWithParameters`.
// NCBI reference (598d8ae6): c++/src/objects/scoremat/scoremat.asn:345-366
// ```asn
// PssmWithParameters ::= SEQUENCE {
//     ...
//     pssm        Pssm,
//     ...
//     params      PssmParameters OPTIONAL
// }
// ```
#[derive(Debug, Clone, PartialEq)]
pub struct PssmWithParameters {
    pub pssm: Pssm,
    pub params: Option<PssmParameters>,
}

/// `Pssm`.
// NCBI reference (598d8ae6): c++/src/objects/scoremat/scoremat.asn:259-297
// ```asn
// Pssm ::= SEQUENCE {
//     ...
//     isProtein       BOOLEAN DEFAULT TRUE,
//     ...
//     identifier      Object-id OPTIONAL,
//     ...
//     numRows         INTEGER,    -- number of rows
//     numColumns      INTEGER,    -- number of columns
//     ...
//     rowLabels       SEQUENCE OF VisibleString OPTIONAL,
//     ...
//     byRow           BOOLEAN DEFAULT FALSE,
//     ...
//     query           Seq-entry OPTIONAL,
//     ...
//     intermediateData    PssmIntermediateData OPTIONAL,
//     ...
//     finalData           PssmFinalData OPTIONAL
// }
// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Pssm {
    /// `isProtein`: `Some` when set (written), `None` for the default TRUE (not written).
    pub is_protein: Option<bool>,
    pub identifier: Option<ObjectId>,
    pub num_rows: i32,
    pub num_columns: i32,
    pub row_labels: Option<Vec<Vec<u8>>>,
    /// `byRow`: `Some` when set (written), `None` for the default FALSE (not written).
    pub by_row: Option<bool>,
    pub query: Option<PssmQuery>,
    pub intermediate_data: Option<PssmIntermediateData>,
    pub final_data: Option<PssmFinalData>,
}

/// `Object-id`.
// NCBI reference (598d8ae6): c++/src/objects/general/general.asn:78-80
// ```asn
// Object-id ::= CHOICE {
//     id INTEGER ,
//     str VisibleString }
// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ObjectId {
    Id(i32),
    Str(Vec<u8>),
}

/// The PSSM's query: the `Seq-entry` `seq` variant, a Bioseq with its ids, the titles among its
/// descriptors (`None`: no `descr` member, or one without a title) and its `Seq-inst`.
// NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:36-40
// ```asn
// Bioseq ::= SEQUENCE {
//     id SET OF Seq-id ,            -- equivalent identifiers
//     descr Seq-descr OPTIONAL , -- descriptors
//     inst Seq-inst ,            -- the sequence data
//     annot SET OF Seq-annot OPTIONAL }
// ```
#[derive(Debug, Clone, PartialEq)]
pub struct PssmQuery {
    pub ids: Vec<PssmSeqId>,
    pub titles: Option<Vec<Vec<u8>>>,
    pub inst: SeqInst,
}

/// The `Seq-id` variants that the reader types (`local`, `gi`, the `Textseq-id` ones and
/// `general`); the reader refuses the others.
// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:24-45
// ```asn
// Seq-id ::= CHOICE {
//     local Object-id ,            -- local use
//     gibbsq INTEGER ,             -- Geninfo backbone seqid
//     gibbmt INTEGER ,             -- Geninfo backbone moltype
//     giim Giimport-id ,           -- Geninfo import id
//     genbank Textseq-id ,
//     embl Textseq-id ,
//     pir Textseq-id ,
//     swissprot Textseq-id ,
//     patent Patent-seq-id ,
//     other Textseq-id ,           -- for historical reasons, 'other' = 'refseq'
//     general Dbtag ,              -- for other databases
//     gi INTEGER ,                 -- GenInfo Integrated Database
//     ddbj Textseq-id ,            -- DDBJ
//     prf Textseq-id ,             -- PRF SEQDB
//     pdb PDB-seq-id ,             -- PDB sequence
//     tpg Textseq-id ,             -- Third Party Annot/Seq Genbank
//     tpe Textseq-id ,             -- Third Party Annot/Seq EMBL
//     tpd Textseq-id ,             -- Third Party Annot/Seq DDBJ
//     gpipe Textseq-id ,           -- Internal NCBI genome pipeline processing ID
//     named-annot-track Textseq-id -- Internal named annotation tracking ID
// }
// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PssmSeqId {
    Local(ObjectId),
    Gi(i64),
    /// A `Textseq-id` variant: its name and position in the CHOICE, and its members.
    Textseq {
        variant: &'static str,
        tag: u8,
        id: TextseqId,
    },
    /// `general`: a `Dbtag`.
    General {
        db: Vec<u8>,
        tag: ObjectId,
    },
}

/// `Textseq-id`.
// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:54-58
// ```asn
// Textseq-id ::= SEQUENCE {
//     name VisibleString OPTIONAL ,
//     accession VisibleString OPTIONAL ,
//     release VisibleString OPTIONAL ,
//     version INTEGER OPTIONAL }
// ```
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TextseqId {
    pub name: Option<Vec<u8>>,
    pub accession: Option<Vec<u8>>,
    pub release: Option<Vec<u8>>,
    pub version: Option<i32>,
}

/// The `Seq-id` variants of `Textseq-id` type, with their positions in the CHOICE.
pub const TEXTSEQ_VARIANTS: [(&str, u8); 12] = [
    ("genbank", 4),
    ("embl", 5),
    ("pir", 6),
    ("swissprot", 7),
    ("other", 9),
    ("ddbj", 12),
    ("prf", 13),
    ("tpg", 15),
    ("tpe", 16),
    ("tpd", 17),
    ("gpipe", 18),
    ("named-annot-track", 19),
];

/// `Seq-inst` without the members LOSAT does not keep.
// NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:262-297
// ```asn
// Seq-inst ::= SEQUENCE {            -- the sequence data itself
//     repr ENUMERATED {              -- representation class
//         ...
//     mol ENUMERATED {               -- molecule class in living organism
//         ...
//     length INTEGER OPTIONAL ,      -- length of sequence in residues
//     fuzz Int-fuzz OPTIONAL ,       -- length uncertainty
//     topology ENUMERATED {          -- topology of molecule
//         ...
//         other (255) } DEFAULT linear ,
//     strand ENUMERATED {            -- strandedness in living organism
//         ...
//         other (255) } OPTIONAL ,   -- default ds for DNA, ss for RNA, pept
//     seq-data Seq-data OPTIONAL ,   -- the sequence
//     ext Seq-ext OPTIONAL ,         -- extensions for special types
//     hist Seq-hist OPTIONAL }       -- sequence history
// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeqInst {
    /// `repr`: the value's name and number (`raw`, 2).
    pub repr: (&'static str, i32),
    /// `mol`: the value's name and number (`aa`, 3).
    pub mol: (&'static str, i32),
    pub length: Option<i32>,
    pub seq_data: Option<SeqData>,
}

/// The values of `Seq-inst.repr`.
pub const SEQ_INST_REPR: [(&str, i32); 10] = [
    ("not-set", 0),
    ("virtual", 1),
    ("raw", 2),
    ("seg", 3),
    ("const", 4),
    ("ref", 5),
    ("consen", 6),
    ("map", 7),
    ("delta", 8),
    ("other", 255),
];

/// The values of `Seq-inst.mol`.
pub const SEQ_INST_MOL: [(&str, i32); 6] = [
    ("not-set", 0),
    ("dna", 1),
    ("rna", 2),
    ("aa", 3),
    ("na", 4),
    ("other", 255),
];

/// The protein `Seq-data` variants that the reader types; it refuses the others.
// NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:346-358
// ```asn
// Seq-data ::= CHOICE {              -- sequence representations
//     iupacna IUPACna ,              -- IUPAC 1 letter nuc acid code
//     iupacaa IUPACaa ,              -- IUPAC 1 letter amino acid code
//     ...
//     ncbieaa NCBIeaa ,              -- extended ASCII 1 letter aa codes
//     ncbipaa NCBIpaa ,              -- amino acid probabilities
//     ncbistdaa NCBIstdaa,           -- consecutive codes for std aas
//     gap Seq-gap                    -- gap types
// }
// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SeqData {
    /// `iupacaa` (`IUPACaa ::= StringStore`).
    Iupacaa(Vec<u8>),
    /// `ncbieaa` (`NCBIeaa ::= StringStore`).
    Ncbieaa(Vec<u8>),
    /// `ncbistdaa` (`NCBIstdaa ::= OCTET STRING`).
    Ncbistdaa(Vec<u8>),
}

/// `PssmIntermediateData`.
// NCBI reference (598d8ae6): c++/src/objects/scoremat/scoremat.asn:197-243
// ```asn
// PssmIntermediateData ::= SEQUENCE {
//     ...
//     resFreqsPerPos              SEQUENCE OF INTEGER OPTIONAL,
//     ...
//     weightedResFreqsPerPos      SEQUENCE OF REAL OPTIONAL,
//     ...
//     freqRatios                  SEQUENCE OF REAL,
//     ...
//     informationContent          SEQUENCE OF REAL OPTIONAL,
//     ...
//     gaplessColumnWeights        SEQUENCE OF REAL OPTIONAL,
//     ...
//     sigma                       SEQUENCE OF REAL OPTIONAL,
//     ...
//     intervalSizes               SEQUENCE OF INTEGER OPTIONAL,
//     ...
//     numMatchingSeqs             SEQUENCE OF INTEGER OPTIONAL,
//     ...
//     numIndeptObsr               SEQUENCE OF REAL OPTIONAL
// }
// ```
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PssmIntermediateData {
    pub res_freqs_per_pos: Option<Vec<i32>>,
    pub weighted_res_freqs_per_pos: Option<Vec<f64>>,
    pub freq_ratios: Vec<f64>,
    pub information_content: Option<Vec<f64>>,
    pub gapless_column_weights: Option<Vec<f64>>,
    pub sigma: Option<Vec<f64>>,
    pub interval_sizes: Option<Vec<i32>>,
    pub num_matching_seqs: Option<Vec<i32>>,
    pub num_indept_obsr: Option<Vec<f64>>,
}

/// `PssmFinalData`.
// NCBI reference (598d8ae6): c++/src/objects/scoremat/scoremat.asn:154-191
// ```asn
// PssmFinalData ::= SEQUENCE {
//     ...
//     scores              SEQUENCE OF INTEGER,
//     ...
//     lambda              REAL,
//     ...
//     kappa               REAL,
//     ...
//     h                   REAL,
//     ...
//     scalingFactor       INTEGER DEFAULT 1,
//     ...
//     lambdaUngapped      REAL OPTIONAL,
//     ...
//     kappaUngapped       REAL OPTIONAL,
//     ...
//     hUngapped           REAL OPTIONAL,
//     ...
//     wordScoreThreshold    REAL OPTIONAL
// }
// ```
#[derive(Debug, Clone, PartialEq)]
pub struct PssmFinalData {
    pub scores: Vec<i32>,
    pub lambda: f64,
    pub kappa: f64,
    pub h: f64,
    /// `scalingFactor`: `Some` when set, `None` for the default 1 (not written).
    pub scaling_factor: Option<i32>,
    pub lambda_ungapped: Option<f64>,
    pub kappa_ungapped: Option<f64>,
    pub h_ungapped: Option<f64>,
    pub word_score_threshold: Option<f64>,
}

/// `PssmParameters`, without `constraints` and `annotatedSites`.
// NCBI reference (598d8ae6): c++/src/objects/scoremat/scoremat.asn:319-341
// ```asn
// PssmParameters ::= SEQUENCE {
//     ...
//     pseudocount INTEGER OPTIONAL,
//     ...
//     rpsdbparams     FormatRpsDbParameters OPTIONAL,
//     ...
//     constraints     CoreDef OPTIONAL,
//     ...
//     bitScoreThresh  REAL OPTIONAL,
//     ...
//     bitScoreReportingThresh  REAL OPTIONAL,
//     ...
//     annotatedSites  Site-annot-set OPTIONAL
// }
// ```
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PssmParameters {
    pub pseudocount: Option<i32>,
    pub rpsdbparams: Option<FormatRpsDbParameters>,
    pub bit_score_thresh: Option<f64>,
    pub bit_score_reporting_thresh: Option<f64>,
}

/// `FormatRpsDbParameters`.
// NCBI reference (598d8ae6): c++/src/objects/scoremat/scoremat.asn:303-315
// ```asn
// FormatRpsDbParameters ::= SEQUENCE {
//     ...
//     matrixName   VisibleString,
//     ...
//     gapOpen      INTEGER OPTIONAL,
//     ...
//     gapExtend    INTEGER OPTIONAL
//     ...
// }
// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FormatRpsDbParameters {
    pub matrix_name: Vec<u8>,
    pub gap_open: Option<i32>,
    pub gap_extend: Option<i32>,
}

/// The six Karlin-Altschul values of a `finalData`, written in place of the object's: psiblast's
/// K/A write-back into the PSSM it searched (`psiblast_impl.cpp:150-168`, quoted at
/// `algorithm::psiblast::in_pssm_files::InPssmRounds::written_back`), which sets all six.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FinalKarlin {
    pub lambda: f64,
    pub kappa: f64,
    pub h: f64,
    pub lambda_ungapped: f64,
    pub kappa_ungapped: f64,
    pub h_ungapped: f64,
}

/// The document of a PSSM (`WriteFileHeader`, the object, `EndOfWrite` and the separator), with
/// `karlin` in place of the six values of `finalData`.
fn write_document<E: Encoder>(
    enc: &mut E,
    value: &PssmWithParameters,
    karlin: Option<FinalKarlin>,
) {
    enc.begin_document("PssmWithParameters");
    enc.begin_sequence("PssmWithParameters");
    enc.member("pssm", 0);
    write_pssm(enc, &value.pssm, karlin);
    enc.end_member();
    if let Some(params) = &value.params {
        enc.member("params", 1);
        write_params(enc, params);
        enc.end_member();
    }
    enc.end_sequence();
    enc.end_document();
}

fn write_object_id<E: Encoder>(enc: &mut E, id: &ObjectId) {
    match id {
        ObjectId::Id(value) => {
            enc.variant("Object-id", "id", 0);
            enc.int(*value);
        }
        ObjectId::Str(value) => {
            enc.variant("Object-id", "str", 1);
            enc.visible_str(value);
        }
    }
    enc.end_variant();
}

fn write_ints<E: Encoder>(enc: &mut E, values: &[i32]) {
    enc.begin_list(false);
    for &value in values {
        enc.element();
        enc.int(value);
    }
    enc.end_list();
}

fn write_reals<E: Encoder>(enc: &mut E, values: &[f64]) {
    enc.begin_list(false);
    for &value in values {
        enc.element();
        enc.real(value);
    }
    enc.end_list();
}

fn write_pssm<E: Encoder>(enc: &mut E, pssm: &Pssm, karlin: Option<FinalKarlin>) {
    enc.begin_sequence("Pssm");
    if let Some(is_protein) = pssm.is_protein {
        enc.member("isProtein", 0);
        enc.boolean(is_protein);
        enc.end_member();
    }
    if let Some(identifier) = &pssm.identifier {
        enc.member("identifier", 1);
        write_object_id(enc, identifier);
        enc.end_member();
    }
    enc.member("numRows", 2);
    enc.int(pssm.num_rows);
    enc.end_member();
    enc.member("numColumns", 3);
    enc.int(pssm.num_columns);
    enc.end_member();
    if let Some(labels) = &pssm.row_labels {
        enc.member("rowLabels", 4);
        enc.begin_list(false);
        for label in labels {
            enc.element();
            enc.visible_str(label);
        }
        enc.end_list();
        enc.end_member();
    }
    if let Some(by_row) = pssm.by_row {
        enc.member("byRow", 5);
        enc.boolean(by_row);
        enc.end_member();
    }
    if let Some(query) = &pssm.query {
        enc.member("query", 6);
        enc.variant("Seq-entry", "seq", 0);
        write_query(enc, query);
        enc.end_variant();
        enc.end_member();
    }
    if let Some(data) = &pssm.intermediate_data {
        enc.member("intermediateData", 7);
        write_intermediate_data(enc, data);
        enc.end_member();
    }
    if let Some(data) = &pssm.final_data {
        enc.member("finalData", 8);
        write_final_data(enc, data, karlin);
        enc.end_member();
    }
    enc.end_sequence();
}

fn write_seq_id<E: Encoder>(enc: &mut E, id: &PssmSeqId) {
    match id {
        PssmSeqId::Local(id) => {
            enc.variant("Seq-id", "local", 0);
            write_object_id(enc, id);
        }
        PssmSeqId::Gi(gi) => {
            enc.variant("Seq-id", "gi", 11);
            enc.big_int(*gi);
        }
        PssmSeqId::Textseq { variant, tag, id } => {
            enc.variant("Seq-id", variant, *tag);
            enc.begin_sequence("Textseq-id");
            for (name, tag, value) in [
                ("name", 0, &id.name),
                ("accession", 1, &id.accession),
                ("release", 2, &id.release),
            ] {
                if let Some(value) = value {
                    enc.member(name, tag);
                    enc.visible_str(value);
                    enc.end_member();
                }
            }
            if let Some(version) = id.version {
                enc.member("version", 3);
                enc.int(version);
                enc.end_member();
            }
            enc.end_sequence();
        }
        PssmSeqId::General { db, tag } => {
            enc.variant("Seq-id", "general", 10);
            enc.begin_sequence("Dbtag");
            enc.member("db", 0);
            enc.visible_str(db);
            enc.end_member();
            enc.member("tag", 1);
            write_object_id(enc, tag);
            enc.end_member();
            enc.end_sequence();
        }
    }
    enc.end_variant();
}

/// The query's Bioseq (`CPsiBlastInputData::x_ExtractQueryForPssm` makes psiblast's).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:144-175
// ```c++
//     m_QueryBioseq.Reset(new CBioseq);
//     // set the sequence id
//     CRef<CSeq_align> aln =
//         const_cast<CSeq_align_set*>(&*m_SeqAlignSet)->Set().front();
//     CRef<CSeq_id> query_id(const_cast<CSeq_id*>(&aln->GetSeq_id(0)));
//     m_QueryBioseq->SetId().push_back(query_id);
//
//     CRef<CSeqdesc> desc(new CSeqdesc);
//     desc->SetTitle(m_QueryTitle);
//     m_QueryBioseq->SetDescr().Set().push_back(desc);
//
//     // set required Seq-inst fields
//     m_QueryBioseq->SetInst().SetRepr(CSeq_inst::eRepr_raw);
//     m_QueryBioseq->SetInst().SetMol(CSeq_inst::eMol_aa);
//     m_QueryBioseq->SetInst().SetLength(GetQueryLength());
//
//     // set the sequence data in ncbistdaa format
//     CNCBIstdaa& seq = m_QueryBioseq->SetInst().SetSeq_data().SetNcbistdaa();
// ```
fn write_query<E: Encoder>(enc: &mut E, query: &PssmQuery) {
    enc.begin_sequence("Bioseq");
    enc.member("id", 0);
    enc.begin_list(true);
    for id in &query.ids {
        enc.element();
        write_seq_id(enc, id);
    }
    enc.end_list();
    enc.end_member();
    if let Some(titles) = &query.titles {
        enc.member("descr", 1);
        enc.begin_named("Seq-descr");
        enc.begin_list(true);
        for title in titles {
            enc.element();
            enc.variant("Seqdesc", "title", 4);
            enc.visible_str(title);
            enc.end_variant();
        }
        enc.end_list();
        enc.end_named();
        enc.end_member();
    }
    enc.member("inst", 2);
    enc.begin_sequence("Seq-inst");
    enc.member("repr", 0);
    enc.enumerated(query.inst.repr.0, query.inst.repr.1);
    enc.end_member();
    enc.member("mol", 1);
    enc.enumerated(query.inst.mol.0, query.inst.mol.1);
    enc.end_member();
    if let Some(length) = query.inst.length {
        enc.member("length", 2);
        enc.int(length);
        enc.end_member();
    }
    if let Some(data) = &query.inst.seq_data {
        enc.member("seq-data", 6);
        match data {
            SeqData::Iupacaa(letters) => {
                enc.variant("Seq-data", "iupacaa", 1);
                enc.begin_named("IUPACaa");
                enc.string_store(letters);
                enc.end_named();
            }
            SeqData::Ncbieaa(letters) => {
                enc.variant("Seq-data", "ncbieaa", 7);
                enc.begin_named("NCBIeaa");
                enc.string_store(letters);
                enc.end_named();
            }
            SeqData::Ncbistdaa(codes) => {
                enc.variant("Seq-data", "ncbistdaa", 9);
                enc.begin_named("NCBIstdaa");
                enc.octets(codes);
                enc.end_named();
            }
        }
        enc.end_variant();
        enc.end_member();
    }
    enc.end_sequence();
    enc.end_member();
    enc.end_sequence();
}

fn write_intermediate_data<E: Encoder>(enc: &mut E, data: &PssmIntermediateData) {
    enc.begin_sequence("PssmIntermediateData");
    if let Some(values) = &data.res_freqs_per_pos {
        enc.member("resFreqsPerPos", 0);
        write_ints(enc, values);
        enc.end_member();
    }
    if let Some(values) = &data.weighted_res_freqs_per_pos {
        enc.member("weightedResFreqsPerPos", 1);
        write_reals(enc, values);
        enc.end_member();
    }
    enc.member("freqRatios", 2);
    write_reals(enc, &data.freq_ratios);
    enc.end_member();
    for (name, tag, values) in [
        ("informationContent", 3, &data.information_content),
        ("gaplessColumnWeights", 4, &data.gapless_column_weights),
        ("sigma", 5, &data.sigma),
    ] {
        if let Some(values) = values {
            enc.member(name, tag);
            write_reals(enc, values);
            enc.end_member();
        }
    }
    for (name, tag, values) in [
        ("intervalSizes", 6, &data.interval_sizes),
        ("numMatchingSeqs", 7, &data.num_matching_seqs),
    ] {
        if let Some(values) = values {
            enc.member(name, tag);
            write_ints(enc, values);
            enc.end_member();
        }
    }
    if let Some(values) = &data.num_indept_obsr {
        enc.member("numIndeptObsr", 8);
        write_reals(enc, values);
        enc.end_member();
    }
    enc.end_sequence();
}

fn write_final_data<E: Encoder>(enc: &mut E, data: &PssmFinalData, karlin: Option<FinalKarlin>) {
    let (lambda, kappa, h, lambda_ungapped, kappa_ungapped, h_ungapped) = match karlin {
        Some(karlin) => (
            karlin.lambda,
            karlin.kappa,
            karlin.h,
            Some(karlin.lambda_ungapped),
            Some(karlin.kappa_ungapped),
            Some(karlin.h_ungapped),
        ),
        None => (
            data.lambda,
            data.kappa,
            data.h,
            data.lambda_ungapped,
            data.kappa_ungapped,
            data.h_ungapped,
        ),
    };
    enc.begin_sequence("PssmFinalData");
    enc.member("scores", 0);
    write_ints(enc, &data.scores);
    enc.end_member();
    for (name, tag, value) in [("lambda", 1, lambda), ("kappa", 2, kappa), ("h", 3, h)] {
        enc.member(name, tag);
        enc.real(value);
        enc.end_member();
    }
    if let Some(factor) = data.scaling_factor {
        enc.member("scalingFactor", 4);
        enc.int(factor);
        enc.end_member();
    }
    for (name, tag, value) in [
        ("lambdaUngapped", 5, lambda_ungapped),
        ("kappaUngapped", 6, kappa_ungapped),
        ("hUngapped", 7, h_ungapped),
        ("wordScoreThreshold", 8, data.word_score_threshold),
    ] {
        if let Some(value) = value {
            enc.member(name, tag);
            enc.real(value);
            enc.end_member();
        }
    }
    enc.end_sequence();
}

fn write_params<E: Encoder>(enc: &mut E, params: &PssmParameters) {
    enc.begin_sequence("PssmParameters");
    if let Some(pseudocount) = params.pseudocount {
        enc.member("pseudocount", 0);
        enc.int(pseudocount);
        enc.end_member();
    }
    if let Some(rps) = &params.rpsdbparams {
        enc.member("rpsdbparams", 1);
        enc.begin_sequence("FormatRpsDbParameters");
        enc.member("matrixName", 0);
        enc.visible_str(&rps.matrix_name);
        enc.end_member();
        if let Some(gap_open) = rps.gap_open {
            enc.member("gapOpen", 1);
            enc.int(gap_open);
            enc.end_member();
        }
        if let Some(gap_extend) = rps.gap_extend {
            enc.member("gapExtend", 2);
            enc.int(gap_extend);
            enc.end_member();
        }
        enc.end_sequence();
        enc.end_member();
    }
    for (name, tag, value) in [
        ("bitScoreThresh", 3, params.bit_score_thresh),
        (
            "bitScoreReportingThresh",
            4,
            params.bit_score_reporting_thresh,
        ),
    ] {
        if let Some(value) = value {
            enc.member(name, tag);
            enc.real(value);
            enc.end_member();
        }
    }
    enc.end_sequence();
}

/// `value` in ASN.1 text (`MSerial_AsnText`), as one document; `Err` when the encoder refuses a
/// value (a byte outside the visible ASCII range in a VisibleString, which the callers refuse
/// first).
pub fn asn_text(value: &PssmWithParameters) -> Result<Vec<u8>, super::encoder::SerialError> {
    asn_text_with_karlin(value, None)
}

/// `asn_text` of `value` with `karlin`, when given, in place of the six values of its `finalData`
/// (the object is not copied).
pub fn asn_text_with_karlin(
    value: &PssmWithParameters,
    karlin: Option<FinalKarlin>,
) -> Result<Vec<u8>, super::encoder::SerialError> {
    let mut out = Vec::new();
    let mut enc = super::asn_text::AsnTextEncoder::new(&mut out);
    write_document(&mut enc, value, karlin);
    if let Some(error) = enc.error() {
        return Err(error);
    }
    Ok(out)
}
