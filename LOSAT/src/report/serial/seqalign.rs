//! The Seq-annot of Seq-aligns that NCBI writes for each query in outfmt 8 (ASN.1 text) and
//! 12 (JSON): NCBI's conversion of a query's final hit list to Seq-aligns
//! (algo/blast/api/blast_seqalign.cpp) and its wrapping in a Seq-annot
//! (algo/blast/format/blastfmtutil.cpp), written through an `Encoder` in the member order of
//! the ASN.1 modules.
//!
//! The documents are made by `write_query`, one per query with alignments, from the report
//! records of `report::record`. The segment computation is that of `s_CollectSeqAlignData`
//! for every program; the forms are NCBI's:
//! - a gapped search whose rows are not translated (BLASTN, BLASTP): one Seq-align per HSP
//!   with a Dense-seg (`write_dense_seg`, `s_CreateDenseg`);
//! - a gapped search with a translated row (BLASTX, TBLASTN): one Seq-align per HSP with one
//!   Std-seg per segment (`write_std_segs`, `s_CreateStdSegs`): Seq-intervals in nucleotides
//!   on the translated row, `empty` for a gap row;
//! - an ungapped translated search (TBLASTX, BLASTX `-ungapped`): one `diags` Seq-align per
//!   subject, one Std-seg with its scores per HSP (`write_seq_align_diags`,
//!   `BLASTUngappedHspListToSeqAlign` / `x_UngappedHSPToStdSeg`).
//! The Dense-diag of an ungapped BLASTN or BLASTP search is not made: neither program runs an
//! ungapped search (`writes`). A `-subject` search's Seq-aligns go through NCBI's mapper of
//! `s_RemapToSubjectLoc` (`remap_strands` for a Dense-seg, `std_seg_kept` for Std-segs).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:677-707
//! ```c++
//! void
//! CBlastFormat::x_PrintStructuredReport(const blast::CSearchResults& results,
//!               CConstRef<blast::CBlastQueryVector> queries)
//! {
//!     string db_title;
//!     if (!m_DbInfo.empty()) {
//!         db_title = m_DbInfo.front().definition;
//!         for (size_t i=1;i < m_DbInfo.size();i++) {
//!             db_title += "; ";
//!             db_title += m_DbInfo[i].definition;
//!         }
//!     }
//!
//!     // ASN.1 formatting is straightforward
//!     if (m_FormatType == CFormattingArgs::eAsnText || m_FormatType == CFormattingArgs::eJsonSeqalign) {
//!         if (results.HasAlignments()) {
//!     	    CRef<CSeq_align_set> aln_set (new CSeq_align_set);
//!             CBlastFormatUtil::PruneSeqalign(*(results.GetSeqAlign()), *aln_set, m_HitlistSize);
//!             if(m_FormatType == CFormattingArgs::eAsnText)
//!             	m_Outfile << MSerial_AsnText << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
//!             else
//!             	m_Outfile << MSerial_Json << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
//!         }
//!         return;
//! ```
//! `PruneSeqalign` keeps the Seq-aligns of the first `m_HitlistSize` subjects, which is the
//! search's own hit list size for these formats, so it keeps every Seq-align of the search's
//! hit lists (outfmt-20261019 S/DESIGN_OF3.md §2.1).

use crate::blastinput::app::Program;
use crate::common::GapEditOp;
use crate::report::record::{HspRecord, QueryRecord, RunRecord, SubjectRecord};
use crate::report::seqid::SeqId;

use super::asn_ber::AsnBinaryEncoder;
use super::asn_text::AsnTextEncoder;
use super::encoder::{Encoder, SerialError};
use super::json::JsonEncoder;

/// The serialized formats of a Seq-annot of Seq-aligns.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SerialFormat {
    /// outfmt 8: `MSerial_AsnText`.
    AsnText,
    /// outfmt 9: `MSerial_AsnBinary`.
    AsnBinary,
    /// outfmt 12: `MSerial_Json`.
    Json,
}

impl SerialFormat {
    /// The Seq-annot format of a report format (`x_PrintStructuredReport` above: `eAsnText`
    /// and `eJsonSeqalign`), `None` for the others.
    pub fn of(format: crate::blastinput::app::ReportFormat) -> Option<Self> {
        use crate::blastinput::app::ReportFormat;
        match format {
            ReportFormat::AsnText => Some(SerialFormat::AsnText),
            ReportFormat::AsnBinary => Some(SerialFormat::AsnBinary),
            ReportFormat::JsonSeqalign => Some(SerialFormat::Json),
            _ => None,
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:61-71
/// ```c++
/// #ifndef SMALLEST_EVALUE
/// /// Threshold below which e-values are saved as 0
/// #define SMALLEST_EVALUE 1.0e-180
/// #endif
/// #ifndef GAP_VALUE
/// /// Value in the Dense-seg indicating a gap
/// #define GAP_VALUE -1
/// #endif
///
/// /// BLAST alignments have always 2 dimensions (i.e.: Pairwise alignment is 2
/// /// dimensional)
/// static const TSeqPos kBlastAlignmentDim = 2;
/// ```
const SMALLEST_EVALUE: f64 = 1.0e-180;

/// The e-value of a Seq-align score list: 0 below `SMALLEST_EVALUE` (shared with BLAST XML,
/// which reads the same Seq-align, `s_SeqAlignSetToXMLHsps`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1186
/// ```c++
///     double evalue = (hsp->evalue < SMALLEST_EVALUE) ? 0.0 : hsp->evalue;
/// ```
pub(crate) fn seqalign_evalue(evalue: f64) -> f64 {
    if evalue < SMALLEST_EVALUE {
        0.0
    } else {
        evalue
    }
}
const GAP_VALUE: i32 = -1;
const BLAST_ALIGNMENT_DIM: i32 = 2;

/// `CODON_LENGTH` (include/algo/blast/core/blast_def.h:62-63: `#define CODON_LENGTH 3`).
const CODON_LENGTH: i32 = 3;

/// The values of `Na-strand` that BLAST writes.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:115-121
/// ```text
/// Na-strand ::= ENUMERATED {          -- strand of nucleic acid
///     unknown (0) ,
///     plus (1) ,
///     minus (2) ,
///     both (3) ,                -- in forward orientation
///     both-rev (4) ,            -- in reverse orientation
///     other (255) }
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    Unknown,
    Plus,
    Minus,
}

impl Strand {
    fn name(self) -> &'static str {
        match self {
            Strand::Unknown => "unknown",
            Strand::Plus => "plus",
            Strand::Minus => "minus",
        }
    }

    fn value(self) -> i32 {
        match self {
            Strand::Unknown => 0,
            Strand::Plus => 1,
            Strand::Minus => 2,
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:73-83
    /// ```c++
    /// /// Converts a frame into the appropriate strand
    /// static ENa_strand
    /// s_Frame2Strand(short frame)
    /// {
    ///     if (frame > 0)
    ///         return eNa_strand_plus;
    ///     else if (frame < 0)
    ///         return eNa_strand_minus;
    ///     else
    ///         return eNa_strand_unknown;
    /// }
    /// ```
    fn of_frame(frame: i8) -> Self {
        if frame > 0 {
            Strand::Plus
        } else if frame < 0 {
            Strand::Minus
        } else {
            Strand::Unknown
        }
    }
}

/// The segments of one HSP: `starts` (2 per segment, query then subject, `GAP_VALUE` for a
/// gap), `lens` and `strands` (2 per segment), as `s_CollectSeqAlignData` fills them, in the
/// searched letters (before the remap of `RemapToQueryLoc` / `s_RemapToSubjectLoc`, which
/// `write_dense_seg` adds). Kept from HSP to HSP to reuse the vectors.
#[derive(Default)]
pub struct Segments {
    starts: Vec<i32>,
    lens: Vec<i32>,
    strands: Vec<Strand>,
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:89-100
/// ```c++
/// static int
/// s_GetCurrPos(int& pos, int pos2advance)
/// {
///     int retval;
///
///     if (pos < 0) /// @todo FIXME: is this condition possible?
///         retval = -(pos + pos2advance - 1);
///     else
///         retval = pos;
///     pos += pos2advance;
///     return retval;
/// }
/// ```
fn get_curr_pos(pos: &mut i32, pos2advance: i32) -> i32 {
    let retval = if *pos < 0 {
        -(*pos + pos2advance - 1)
    } else {
        *pos
    };
    *pos += pos2advance;
    retval
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:114-145
/// ```c++
/// static TSeqPos
/// s_GetAlignmentStart(int& curr_pos, int num,
///         ENa_strand strand, bool translate, int length, int original_length,
///         short frame)
/// {
///     TSeqPos retval;
///
///     if (strand == eNa_strand_minus) {
///
///         if (translate)
///             retval = original_length -
///                 CODON_LENGTH*(s_GetCurrPos(curr_pos, num) + num)
///                 + frame + 1;
///         else
///             retval = length - s_GetCurrPos(curr_pos, num) - num;
///
///     } else {
///
///         if (translate)
///             retval = frame - 1 + CODON_LENGTH*s_GetCurrPos(curr_pos, num);
///         else
///             retval = s_GetCurrPos(curr_pos, num);
///
///     }
///
///     return retval;
/// }
/// ```
/// (The result, a `TSeqPos`, is stored in a `TSignedSeqPos`; it is not negative.)
#[allow(clippy::too_many_arguments)]
fn get_alignment_start(
    curr_pos: &mut i32,
    num: i32,
    strand: Strand,
    translate: bool,
    length: i32,
    original_length: i32,
    frame: i8,
) -> i32 {
    let frame = i32::from(frame);
    if strand == Strand::Minus {
        if translate {
            original_length - CODON_LENGTH * (get_curr_pos(curr_pos, num) + num) + frame + 1
        } else {
            length - get_curr_pos(curr_pos, num) - num
        }
    } else if translate {
        frame - 1 + CODON_LENGTH * get_curr_pos(curr_pos, num)
    } else {
        get_curr_pos(curr_pos, num)
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:147-151
/// ```c++
/// static Int4
/// s_GetProteinFrameLength(Int4 nuc_length, Int2 frame)
/// {
///     return (nuc_length - (ABS(frame)-1)%CODON_LENGTH) / CODON_LENGTH;
/// }
/// ```
fn get_protein_frame_length(nuc_length: i32, frame: i8) -> i32 {
    (nuc_length - (i32::from(frame).abs() - 1) % CODON_LENGTH) / CODON_LENGTH
}

/// The segments of `hsp` (`s_CollectSeqAlignData` over the whole edit script): one per
/// edit-script operation, in script order.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:168-268
/// ```c++
/// static void
/// s_CollectSeqAlignData(const BlastHSP* hsp, const GapEditScript* esp,
///                       unsigned int first, unsigned int nsegs,
///                       CDense_seg::TStarts& starts,
///                       CDense_seg::TLens& lengths,
///                       CDense_seg::TStrands& strands,
///                       Int4 query_length, Int4 subject_length,
///                       bool translate1, bool translate2)
/// {
///     ...
///     int start1 = hsp->query.offset;   // start of alignment on master sequence
///     int start2 = hsp->subject.offset; // start of alignment on slave sequence
///     int length1 = query_length;
///     int length2 = subject_length;
///     ...
///     if (translate1)
///         length1 = s_GetProteinFrameLength(length1, hsp->query.frame);
///     if (translate2)
///         length2 = s_GetProteinFrameLength(length2, hsp->subject.frame);
///
///     m_strand = s_Frame2Strand(hsp->query.frame);
///     s_strand = s_Frame2Strand(hsp->subject.frame);
///
///     for (unsigned int esp_index = first; esp_index< (unsigned int)esp->size && esp_index < (unsigned int)(first+nsegs); esp_index++) {
///         switch (esp->op_type[esp_index]) {
///         case eGapAlignDecline:
///         case eGapAlignSub:
///             m_start =
///                 s_GetAlignmentStart(start1, esp->num[esp_index], m_strand, translate1, length1,
///                                     query_length, hsp->query.frame);
///
///             s_start =
///                 s_GetAlignmentStart(start2, esp->num[esp_index], s_strand, translate2, length2,
///                                     subject_length, hsp->subject.frame);
///
///             strands.push_back(m_strand);
///             strands.push_back(s_strand);
///             starts.push_back(m_start);
///             starts.push_back(s_start);
///             break;
///
///         // Insertion on the master sequence (gap on slave)
///         case eGapAlignIns:
///             m_start =
///                 s_GetAlignmentStart(start1, esp->num[esp_index], m_strand, translate1, length1,
///                                     query_length, hsp->query.frame);
///
///             s_start = GAP_VALUE;
///
///             strands.push_back(m_strand);
///             strands.push_back(esp_index == 0 ? eNa_strand_unknown : s_strand);
///             starts.push_back(m_start);
///             starts.push_back(s_start);
///             break;
///
///         // Deletion on master sequence (gap; insertion on slave)
///         case eGapAlignDel:
///             m_start = GAP_VALUE;
///
///             s_start =
///                 s_GetAlignmentStart(start2, esp->num[esp_index], s_strand, translate2, length2,
///                                     subject_length, hsp->subject.frame);
///
///             strands.push_back(esp_index == 0 ? eNa_strand_unknown : m_strand);
///             strands.push_back(s_strand);
///             starts.push_back(m_start);
///             starts.push_back(s_start);
///             break;
///
///         default:
///             break;
///         }
///
///         lengths.push_back(esp->num[esp_index]);
///     }
/// ```
/// An empty script (no NCBI counterpart: every gapped `BlastHSP` has one) stands for one
/// substitution over the query's extent (`HspRecord::ops`).
fn collect_seq_align_data(
    hsp: &HspRecord<'_>,
    query_length: i32,
    subject_length: i32,
    translate1: bool,
    translate2: bool,
    segments: &mut Segments,
) {
    segments.starts.clear();
    segments.lens.clear();
    segments.strands.clear();
    let mut start1 = hsp.q_offset;
    let mut start2 = hsp.s_offset;
    let mut length1 = query_length;
    let mut length2 = subject_length;
    if translate1 {
        length1 = get_protein_frame_length(length1, hsp.q_frame);
    }
    if translate2 {
        length2 = get_protein_frame_length(length2, hsp.s_frame);
    }
    let m_strand = Strand::of_frame(hsp.q_frame);
    let s_strand = Strand::of_frame(hsp.s_frame);
    let whole = [GapEditOp::Sub((hsp.q_end - hsp.q_offset) as u32)];
    let ops = if hsp.ops.is_empty() {
        &whole[..]
    } else {
        hsp.ops
    };
    for (esp_index, op) in ops.iter().enumerate() {
        let num = op.num() as i32;
        match *op {
            GapEditOp::Sub(_) => {
                let m_start = get_alignment_start(
                    &mut start1,
                    num,
                    m_strand,
                    translate1,
                    length1,
                    query_length,
                    hsp.q_frame,
                );
                let s_start = get_alignment_start(
                    &mut start2,
                    num,
                    s_strand,
                    translate2,
                    length2,
                    subject_length,
                    hsp.s_frame,
                );
                segments.strands.extend([m_strand, s_strand]);
                segments.starts.extend([m_start, s_start]);
            }
            GapEditOp::Ins(_) => {
                let m_start = get_alignment_start(
                    &mut start1,
                    num,
                    m_strand,
                    translate1,
                    length1,
                    query_length,
                    hsp.q_frame,
                );
                segments.strands.extend([
                    m_strand,
                    if esp_index == 0 {
                        Strand::Unknown
                    } else {
                        s_strand
                    },
                ]);
                segments.starts.extend([m_start, GAP_VALUE]);
            }
            GapEditOp::Del(_) => {
                let s_start = get_alignment_start(
                    &mut start2,
                    num,
                    s_strand,
                    translate2,
                    length2,
                    subject_length,
                    hsp.s_frame,
                );
                segments.strands.extend([
                    if esp_index == 0 {
                        Strand::Unknown
                    } else {
                        m_strand
                    },
                    s_strand,
                ]);
                segments.starts.extend([GAP_VALUE, s_start]);
            }
        }
        segments.lens.push(num);
    }
}

/// `query_length` of the HSP coverage score (`BLASTHspListToSeqAlign` passes a translated
/// query's length in protein letters, Int4 division; the ungapped conversion passes the
/// length itself).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1497-1503
/// ```c++
///         if  (Blast_QueryIsTranslated(program)) {
///         	s_AddScoresToSeqAlign(seqalign, hsp, emptyList, query_length/3);
///         }
///         else {
///         	// Pass in empty list until removed.
///         	s_AddScoresToSeqAlign(seqalign, hsp, emptyList, query_length);
///         }
/// ```
pub fn coverage_query_length(run: &RunRecord<'_>, searched_length: i32) -> i32 {
    if run.gapped && query_is_translated(run.program) {
        searched_length / 3
    } else {
        searched_length
    }
}

/// `Blast_QueryIsTranslated` of the program (blastx, tblastx).
fn query_is_translated(program: Program) -> bool {
    matches!(program, Program::Blastx | Program::Tblastx)
}

/// `Blast_SubjectIsTranslated` of the program (tblastn, tblastx).
fn subject_is_translated(program: Program) -> bool {
    matches!(program, Program::Tblastn | Program::Tblastx)
}

/// The HSP's coverage of the query, in percent (`hsp_percent_coverage`; NCBI's tabular
/// `qcovhsp` is this value cut to an integer).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1034-1043
/// ```c
/// double Blast_HSPGetQueryCoverage(const BlastHSP* hsp, Int4 query_length)
/// {
/// 	double pct = 0;
///     if(query_length > 0) {
///             pct = 100.0 * (double) (hsp->query.end - hsp->query.offset)/ (double) query_length;
///             if(pct < 99)
///             	pct +=0.5;
///     }
///     return pct;
/// }
/// ```
pub fn query_coverage(q_offset: i32, q_end: i32, query_length: i32) -> f64 {
    let mut pct = 0.0;
    if query_length > 0 {
        pct = 100.0 * f64::from(q_end - q_offset) / f64::from(query_length);
        if pct < 99.0 {
            pct += 0.5;
        }
    }
    pct
}

/// The "Blast Type" of the program: `EProgramToTaskName(ProgramNameToEnum(program name))`
/// and the `EProgram` value. The program name of every BLASTN task is "blastn"
/// (`Blast_ProgramNameFromType`, blast_format.cpp:94), and of every BLASTP task "blastp".
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_types.hpp:56-62
/// ```c++
/// enum EProgram {
///     eBlastNotSet = 0,   ///< Not yet set.
///     eBlastn,            ///< Nucl-Nucl (traditional blastn)
///     eBlastp,            ///< Protein-Protein
///     eBlastx,            ///< Translated nucl-Protein
///     eTblastn,           ///< Protein-Translated nucl
///     eTblastx,           ///< Translated nucl-Translated nucl
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:676-686
/// ```c++
/// string EProgramToTaskName(EProgram p)
/// {
///     string retval;
///     switch (p) {
///     case eBlastn:           retval.assign("blastn"); break;
///     case eMegablast:        retval.assign("megablast"); break;
///     case eDiscMegablast:    retval.assign("dc-megablast"); break;
///     case eBlastp:           retval.assign("blastp"); break;
///     case eBlastx:           retval.assign("blastx"); break;
///     case eTblastn:          retval.assign("tblastn"); break;
///     case eTblastx:          retval.assign("tblastx"); break;
/// ```
fn blast_type(program: Program) -> (&'static str, i32) {
    match program {
        Program::Blastn => ("blastn", 1),
        Program::Blastp => ("blastp", 2),
        Program::Blastx => ("blastx", 3),
        Program::Tblastn => ("tblastn", 4),
        Program::Tblastx => ("tblastx", 5),
    }
}

/// `Blast_SubjectIsNucleotide` of the program.
fn subject_is_nucleotide(program: Program) -> bool {
    matches!(
        program,
        Program::Blastn | Program::Tblastn | Program::Tblastx
    )
}

/// The value of a `User-field`'s `data` that BLAST sets.
enum UserData {
    Bool(bool),
    Int(i32),
}

/// A `user` Annotdesc with one field (`CUser_object::AddField`).
///
/// NCBI reference (598d8ae6): c++/src/objects/general/general.asn:128-136
/// ```text
/// User-object ::= SEQUENCE {
///     class VisibleString OPTIONAL ,   -- endeavor which designed this object
///     type Object-id ,                 -- type of object within class
///     data SEQUENCE OF User-field }    -- the object itself
///
/// User-field ::= SEQUENCE {
///     label Object-id ,                -- field label
///     num INTEGER OPTIONAL ,           -- required for strs, ints, reals, oss
///     data CHOICE {                    -- field contents
/// ```
/// (`data`: `str` 0, `int` 1, `real` 2, `bool` 3, general.asn:136-140.)
fn write_user_object<E: Encoder>(enc: &mut E, type_str: &[u8], label: &[u8], data: UserData) {
    // Annotdesc ::= CHOICE { name, title, comment, pub, user User-object, ... }
    // (objects/seq/seq.asn:440-445).
    enc.variant("Annotdesc", "user", 4);
    enc.begin_sequence("User-object");
    enc.member("type", 1);
    enc.variant("Object-id", "str", 1);
    enc.visible_str(type_str);
    enc.end_variant();
    enc.end_member();
    enc.member("data", 2);
    enc.begin_list(false);
    enc.element();
    enc.begin_sequence("User-field");
    enc.member("label", 0);
    enc.variant("Object-id", "str", 1);
    enc.visible_str(label);
    enc.end_variant();
    enc.end_member();
    enc.member("data", 2);
    match data {
        UserData::Int(value) => {
            enc.variant("User-field_data", "int", 1);
            enc.int(value);
        }
        UserData::Bool(value) => {
            enc.variant("User-field_data", "bool", 3);
            enc.boolean(value);
        }
    }
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
    enc.end_list();
    enc.end_member();
    enc.end_sequence();
    enc.end_variant();
}

/// The Seq-annot's `desc` (an Annot-descr, SET OF Annotdesc): the three user objects.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:331-378
/// ```c++
/// CRef<objects::CSeq_annot>
/// CBlastFormatUtil::CreateSeqAnnotFromSeqAlignSet(const objects::CSeq_align_set & alnset,
/// 												blast::EProgram program,
/// 												const string & db_name,
///                                                 const string & db_title,
/// 												bool vdb_search)
/// {
///     CRef<CSeq_annot> retval(new CSeq_annot);
///
///     //Fill in Hist Seqalign
///     CRef<CUser_object> hist_align_obj(new CUser_object);
///     static const string kHistSeqalign("Hist Seqalign");
///     hist_align_obj->SetType().SetStr(kHistSeqalign);
///     hist_align_obj->AddField(kHistSeqalign, true);
///     retval->AddUserObject(*hist_align_obj);
///
///     //Fill in Blast Type
///     CRef<CUser_object> blast_type(new CUser_object);
///     static const string kBlastType("Blast Type");
///     blast_type->SetType().SetStr(kBlastType);
///     blast_type->AddField(blast::EProgramToTaskName(program), program);
///     retval->AddUserObject(*blast_type);
///
///     //Fill in DB Title
///     CRef<CUser_object> blast_db_info(new CUser_object);
///     if(vdb_search)
///     {
///     	...
///     }
///     else
///     {
///     	static const string kBlastDBTitle("Blast Database Title");
///     	blast_db_info->SetType().SetStr(kBlastDBTitle);
///     	if(0 == db_name.size() || 0 == NStr::CompareNocase(db_name, "n/a"))
///     	{
///     		blast_db_info->AddField( "n/a", false );
///     	}
///     	else if(0 == NStr::CompareNocase(db_name, "SRA"))
///     	{
///     		blast_db_info->AddField( db_name, true );
///     	}
///     	else
///     	{
///     		bool is_nucl = Blast_SubjectIsNucleotide(EProgramToEBlastProgramType(program));
///     		blast_db_info->AddField( db_title, is_nucl );
///     	}
///     }
///
///    	retval->AddUserObject(*blast_db_info);
/// ```
/// A `-subject` search has no database name (`n/a`). A `-db` named `n/a` or `SRA` is
/// refused by the programs (`check_database`), so its title is written.
fn write_desc<E: Encoder>(enc: &mut E, run: &RunRecord<'_>) {
    enc.member("desc", 3);
    enc.begin_list(true);
    enc.element();
    write_user_object(
        enc,
        b"Hist Seqalign",
        b"Hist Seqalign",
        UserData::Bool(true),
    );
    let (task, program) = blast_type(run.program);
    enc.element();
    write_user_object(enc, b"Blast Type", task.as_bytes(), UserData::Int(program));
    enc.element();
    match run.db_title {
        None => write_user_object(enc, b"Blast Database Title", b"n/a", UserData::Bool(false)),
        Some(title) => write_user_object(
            enc,
            b"Blast Database Title",
            title,
            UserData::Bool(subject_is_nucleotide(run.program)),
        ),
    }
    enc.end_list();
    enc.end_member();
}

/// LOSAT's refusals of a database search (`-db name`, title `title`: NCBI's
/// `m_DbInfo.front().definition`, `PairwiseDatabase::definition`) for a Seq-annot format,
/// made before the search:
/// - a database named `n/a` or `SRA` (compared without case), for which NCBI writes another
///   third user object (`CreateSeqAnnotFromSeqAlignSet` above);
/// - for the ASN.1 text, a title with a byte outside 0x20..0x7E, which NCBI's stream
///   replaces by `#` with an error message on standard error (the default
///   `eFNP_ReplaceAndWarn`).
///
/// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:268-279
/// ```c++
/// NCBI_PARAM_ENUM_ARRAY(EFixNonPrint, SERIAL, WRONG_CHARS_WRITE)
/// {
///     {"SKIP",               eFNP_Skip},
///     ...
/// };
/// NCBI_PARAM_ENUM_DECL(EFixNonPrint, SERIAL, WRONG_CHARS_WRITE);
/// NCBI_PARAM_ENUM_DEF(EFixNonPrint, SERIAL, WRONG_CHARS_WRITE, eFNP_ReplaceAndWarn);
/// ```
/// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1855-1871
/// ```c++
/// char ReplaceVisibleChar(char c, EFixNonPrint fix_method,
///     const CObjectStack* io, const CTempString& str, char subst)
/// {
///     ...
///     if ( fix_method != eFNP_Replace ) {
///         string message;
///         if (io != NULL) {
///             message += io->GetStackTrace() + "\n";
///         }
///         message += "Bad char [0x" +
///                    NStr::NumericToString((unsigned char)c,0,16)+
///                    "] in string";
/// ```
pub fn check_database(
    format: SerialFormat,
    name: &str,
    title: &[u8],
    program: Program,
) -> anyhow::Result<()> {
    let number = match format {
        SerialFormat::AsnText => 8,
        SerialFormat::AsnBinary => 9,
        SerialFormat::Json => 12,
    };
    if name.eq_ignore_ascii_case("n/a") || name.eq_ignore_ascii_case("SRA") {
        anyhow::bail!(
            "a database named \"{name}\" (which NCBI BLAST+ describes as no database or as SRA in output format {number}) is not supported by LOSAT's {}",
            program.name()
        );
    }
    if matches!(format, SerialFormat::AsnText | SerialFormat::AsnBinary)
        && title.iter().any(|byte| !(0x20..=0x7e).contains(byte))
    {
        anyhow::bail!(
            "a database title with a character outside the visible ASCII range (which NCBI BLAST+ replaces by '#' with an error message in output format {number}) is not supported by LOSAT's {}",
            program.name()
        );
    }
    Ok(())
}

/// A Score with an integer value.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1098-1110
/// ```c++
/// static CRef<CScore>
/// s_MakeScore(const string& ident_string, double d, int i, bool is_integer)
/// {
///     CRef<CScore> retval(new CScore());
///     retval->SetId().SetStr(ident_string);
///
///     if (is_integer)
///         retval->SetValue().SetInt(i);
///     else
///         retval->SetValue().SetReal(d);
///
///     return retval;
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:283-289
/// ```text
/// Score ::= SEQUENCE {
///     id Object-id OPTIONAL ,
///     value CHOICE {
///         real REAL ,
///         int INTEGER
///     }
/// }
/// ```
fn write_score_int<E: Encoder>(enc: &mut E, name: &str, value: i32) {
    write_score_id(enc, name);
    enc.member("value", 1);
    enc.variant("Score_value", "int", 1);
    enc.int(value);
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
}

/// A Score with a real value (`s_MakeScore(..., false)`, see `write_score_int`).
fn write_score_real<E: Encoder>(enc: &mut E, name: &str, value: f64) {
    write_score_id(enc, name);
    enc.member("value", 1);
    enc.variant("Score_value", "real", 0);
    enc.real(value);
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
}

/// The start of a Score: its `id str`.
fn write_score_id<E: Encoder>(enc: &mut E, name: &str) {
    enc.element();
    enc.begin_sequence("Score");
    enc.member("id", 0);
    enc.variant("Object-id", "str", 1);
    enc.visible_str(name.as_bytes());
    enc.end_variant();
    enc.end_member();
}

/// The scores of an HSP, in NCBI's order and under its conditions, as the elements of a
/// SET OF Score that the caller has opened (`score` of a gapped Seq-align, `scores` of an
/// ungapped Std-seg).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1162-1231
/// ```c++
/// static void
/// s_BuildScoreList(const BlastHSP     * hsp,
///                  CSeq_align::TScore & scores,
///                  const vector<string>  & seqid_list,
///                  Int4 query_length)
/// {
///     ...
///     if (hsp->score) {
///         static const string kScore("score");
///         scores.push_back(s_MakeScore(kScore, 0.0, hsp->score, true));
///         static const string kBlastScore("blast_score");
///         scores.push_back(s_MakeScore(kBlastScore, 0.0, hsp->score, true));
///     }
///
///     if (hsp->num > 1) {
///         static const string kSumN("sum_n");
///         scores.push_back(s_MakeScore(kSumN, 0.0, hsp->num, true));
///     }
///
///     // Set the E-Value
///     double evalue = (hsp->evalue < SMALLEST_EVALUE) ? 0.0 : hsp->evalue;
///     if (evalue >= 0.0) {
///         string score_type = (hsp->num <= 1) ? "e_value" : "sum_e";
///         scores.push_back(s_MakeScore(score_type, evalue, 0, false));
///     }
///
///     // Calculate the bit score from the raw score
///
///     if (hsp->bit_score >= 0.0) {
///         static const string kBitScore("bit_score");
///         scores.push_back(s_MakeScore(kBitScore, hsp->bit_score, 0, false));
///     }
///
///     // Set the identity score
///     if (hsp->num_ident >= 0) {
///         static const string kNumIdent("num_ident");
///         scores.push_back(s_MakeScore(kNumIdent, 0.0, hsp->num_ident, true));
///     }
///
///     if (hsp->comp_adjustment_method > 0) {
///         static const string kCompAdj("comp_adjustment_method");
///         scores.push_back(s_MakeScore(kCompAdj, 0.0,
///                                      hsp->comp_adjustment_method, true));
///     }
///
///     if ( !seqid_list.empty() ) {
///         ...
///     }
///
///     if (hsp->num_positives > 0) {
///             static const string kNumPositives("num_positives");
///             scores.push_back(s_MakeScore(kNumPositives, 0.0, hsp->num_positives, true));
///     }
///
///     if(query_length > 0) {
///             static const string kQueryCovHsp("hsp_percent_coverage");
///             double hsp_coverage = Blast_HSPGetQueryCoverage( hsp, query_length);
///             scores.push_back(s_MakeScore(kQueryCovHsp, hsp_coverage, 0, false));
///     }
///     return;
/// }
/// ```
/// The `seqid_list` entries exist only with a GI or Seq-id list, which LOSAT does not have.
fn write_score_list<E: Encoder>(enc: &mut E, hsp: &HspRecord<'_>, query_length: i32) {
    if hsp.score != 0 {
        write_score_int(enc, "score", hsp.score);
        write_score_int(enc, "blast_score", hsp.score);
    }
    if hsp.num > 1 {
        write_score_int(enc, "sum_n", hsp.num);
    }
    let evalue = seqalign_evalue(hsp.evalue);
    if evalue >= 0.0 {
        let score_type = if hsp.num <= 1 { "e_value" } else { "sum_e" };
        write_score_real(enc, score_type, evalue);
    }
    if hsp.bit_score >= 0.0 {
        write_score_real(enc, "bit_score", hsp.bit_score);
    }
    if hsp.num_ident >= 0 {
        write_score_int(enc, "num_ident", hsp.num_ident);
    }
    if hsp.comp_adjust > 0 {
        write_score_int(enc, "comp_adjustment_method", i32::from(hsp.comp_adjust));
    }
    if hsp.num_positives > 0 {
        write_score_int(enc, "num_positives", hsp.num_positives);
    }
    if query_length > 0 {
        write_score_real(
            enc,
            "hsp_percent_coverage",
            query_coverage(hsp.q_offset, hsp.q_end, query_length),
        );
    }
}

/// The Dense-seg of an HSP of a non-translated search, with the remap of the query row
/// (`RemapToQueryLoc`: the `-query_loc` start added to every start that is not a gap) and
/// of the subject row (`s_RemapToSubjectLoc` of a `-subject` search: the `-subject_loc`
/// start, on the plus strand that LOSAT's ranges have).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:271-301
/// ```c++
/// static void
/// s_CreateDenseg(CDense_seg& dense_seg,
///                CRef<CSeq_id> master, CRef<CSeq_id> slave,
///                CDense_seg::TStarts& starts,
///                CDense_seg::TLens& lengths,
///                CDense_seg::TStrands& strands)
/// {
///     ...
///     dense_seg.SetDim(kBlastAlignmentDim);
///
///     // Set the sequence ids
///     CDense_seg::TIds & ids = dense_seg.SetIds();
///     ids.reserve(kBlastAlignmentDim);
///
///     ids.push_back(master);
///     ids.push_back(slave);
///
///     dense_seg.SetNumseg((int) lengths.size());
///     dense_seg.SetLens().swap(lengths);
///     dense_seg.SetStrands().swap(strands);
///     dense_seg.SetStarts().swap(starts);
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:890-918
/// ```c++
/// void CDense_seg::OffsetRow(TDim row,
///                            TSignedSeqPos offset)
/// {
///     if (offset == 0) return;
///     ...
///     // Modify positions
///     for (TNumseg seg = 0, pos = row;
///          seg < GetNumseg();
///          ++seg, pos += GetDim()) {
///         if (GetStarts()[pos] >= 0) {
///             SetStarts()[pos] += offset;
///         }
///     }
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:74-81
/// ```text
/// Dense-seg ::= SEQUENCE {          -- for (multiway) global or partial alignments
///     dim INTEGER DEFAULT 2 ,       -- dimensionality
///     numseg INTEGER ,              -- number of segments here
///     ids SEQUENCE OF Seq-id ,      -- sequences in order
///     starts SEQUENCE OF INTEGER ,  -- start OFFSETS in ids order within segs
///     lens SEQUENCE OF INTEGER ,    -- lengths in ids order within segs
///     strands SEQUENCE OF Na-strand OPTIONAL ,
///     scores SEQUENCE OF Score OPTIONAL }  -- score for each seg
/// ```
fn write_dense_seg<E: Encoder>(
    enc: &mut E,
    query_id: &SeqId<'_>,
    subject_id: &SeqId<'_>,
    segments: &Segments,
    offsets: [i32; 2],
) {
    enc.begin_sequence("Dense-seg");
    enc.member("dim", 0);
    enc.int(BLAST_ALIGNMENT_DIM);
    enc.end_member();
    enc.member("numseg", 1);
    enc.int(segments.lens.len() as i32);
    enc.end_member();
    enc.member("ids", 2);
    enc.begin_list(false);
    enc.element();
    query_id.write(enc);
    enc.element();
    subject_id.write(enc);
    enc.end_list();
    enc.end_member();
    enc.member("starts", 3);
    enc.begin_list(false);
    for (index, &start) in segments.starts.iter().enumerate() {
        enc.element();
        let offset = offsets[index % 2];
        enc.int(if start >= 0 { start + offset } else { start });
    }
    enc.end_list();
    enc.end_member();
    enc.member("lens", 4);
    enc.begin_list(false);
    for &len in &segments.lens {
        enc.element();
        enc.int(len);
    }
    enc.end_list();
    enc.end_member();
    enc.member("strands", 5);
    enc.begin_list(false);
    for &strand in &segments.strands {
        enc.element();
        enc.enumerated(strand.name(), strand.value());
    }
    enc.end_list();
    enc.end_member();
    enc.end_sequence();
}

/// The strands of a Dense-seg after `s_RemapToSubjectLoc` of a `-subject` search: NCBI maps
/// every Seq-align of a subject whose source can return partial sequences through a
/// `CSeq_loc_Mapper` (the subject's location: the record or its `-subject_loc` interval),
/// which rebuilds the Dense-seg. A row that is not a gap keeps its strand, `unknown` becoming
/// `plus`; a gap takes the strand of the row's first segment that is not a gap (`unknown`
/// becoming `plus`). (The starts move by the interval's start, `write_dense_seg`.)
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1541-1547
/// ```c++
/// static void
/// s_RemapToSubjectLoc(CRef<CSeq_align> & subj_aligns, const CSeq_loc& subj_loc)
/// {
///     const int kSubjDimension = 1;
///     _ASSERT(subj_loc.IsInt() || subj_loc.IsWhole());
///     subj_aligns.Reset(RemapAlignToLoc(*subj_aligns, kSubjDimension, subj_loc));
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1242-1269
/// ```c++
/// void CSeq_align_Mapper_Base::x_FillKnownStrands(TStrands& strands) const
/// {
///     ...
///     for (size_t r_idx = 0; r_idx < max_rows; r_idx++) {
///         ENa_strand strand = eNa_strand_unknown;
///         // Skip gaps, try find a row with mapped strand
///         ITERATE(TSegments, seg_it, m_Segs) {
///             // Make sure the row exists in the current segment.
///             if (seg_it->m_Rows.size() <= r_idx) continue;
///             if (seg_it->m_Rows[r_idx].GetSegStart() != -1) {
///                 strand = seg_it->m_Rows[r_idx].m_Strand;
///                 break;
///             }
///         }
///         // Store the strand.
///         strands.push_back(strand == eNa_strand_unknown ?
///                           eNa_strand_plus : strand);
///     }
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1397-1402
/// ```c++
///             if (m_HaveStrands) {
///                 // For gaps use the strand of the first mapped row
///                 dseg.SetStrands().
///                     push_back((TSeqPos)row->GetSegStart() != kInvalidSeqPos ?
///                     (row->m_Strand != eNa_strand_unknown ?
///                     row->m_Strand : eNa_strand_plus): strands[str_idx]);
///             }
/// ```
/// A `-subject` record's location is an interval with the reader's strand
/// (blast_fasta_input.cpp:418-430: `unknown` for a protein, `both` for a nucleotide), so the
/// mapping keeps the strand of a mapped plus or minus row.
fn remap_strands(segments: &mut Segments) {
    let numseg = segments.lens.len();
    let mut known = [Strand::Plus; 2];
    for (row, known) in known.iter_mut().enumerate() {
        if let Some(seg) = (0..numseg).find(|&seg| segments.starts[2 * seg + row] != GAP_VALUE) {
            let strand = segments.strands[2 * seg + row];
            *known = if strand == Strand::Unknown {
                Strand::Plus
            } else {
                strand
            };
        }
    }
    for (index, strand) in segments.strands.iter_mut().enumerate() {
        *strand = if segments.starts[index] != GAP_VALUE {
            if *strand == Strand::Unknown {
                Strand::Plus
            } else {
                *strand
            }
        } else {
            known[index % 2]
        };
    }
}

/// A Seq-loc `int`: the Seq-interval `from`..`to` on `strand` of `id` (BLAST sets the strand
/// of every interval it makes, `unknown` included, so it is always written).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:79-99
/// ```text
/// Seq-loc ::= CHOICE {
///     null NULL ,           -- not placed
///     empty Seq-id ,        -- to NULL one Seq-id in a collection
///     whole Seq-id ,        -- whole sequence
///     int Seq-interval ,    -- from to
///     ...
/// Seq-interval ::= SEQUENCE {
///     from INTEGER ,
///     to INTEGER ,
///     strand Na-strand OPTIONAL ,
///     id Seq-id ,    -- WARNING: this used to be optional
///     fuzz-from Int-fuzz OPTIONAL ,
///     fuzz-to Int-fuzz OPTIONAL }
/// ```
fn write_seq_interval<E: Encoder>(enc: &mut E, id: &SeqId<'_>, from: i32, to: i32, strand: Strand) {
    enc.variant("Seq-loc", "int", 3);
    enc.begin_sequence("Seq-interval");
    enc.member("from", 0);
    enc.int(from);
    enc.end_member();
    enc.member("to", 1);
    enc.int(to);
    enc.end_member();
    enc.member("strand", 2);
    enc.enumerated(strand.name(), strand.value());
    enc.end_member();
    enc.member("id", 3);
    id.write(enc);
    enc.end_member();
    enc.end_sequence();
    enc.end_variant();
}

/// A Seq-loc `empty`: the gap row of a Std-seg (`master_loc->SetEmpty(*master)`, see
/// `write_std_segs`).
fn write_seq_loc_empty<E: Encoder>(enc: &mut E, id: &SeqId<'_>) {
    enc.variant("Seq-loc", "empty", 1);
    id.write(enc);
    enc.end_variant();
}

/// The start of a Std-seg (an element of the caller's SEQUENCE OF Std-seg): `dim` and `ids`,
/// then `loc` opened for its two Seq-locs (closed by `end_std_seg_locs`).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:94-98
/// ```text
/// Std-seg ::= SEQUENCE {
///     dim INTEGER DEFAULT 2 ,       -- dimensionality
///     ids SEQUENCE OF Seq-id OPTIONAL ,
///     loc SEQUENCE OF Seq-loc ,
///     scores SET OF Score OPTIONAL }
/// ```
/// (BLAST sets `dim` to `kBlastAlignmentDim`, so it is written.)
fn begin_std_seg<E: Encoder>(enc: &mut E, query_id: &SeqId<'_>, subject_id: &SeqId<'_>) {
    enc.element();
    enc.begin_sequence("Std-seg");
    enc.member("dim", 0);
    enc.int(BLAST_ALIGNMENT_DIM);
    enc.end_member();
    enc.member("ids", 1);
    enc.begin_list(false);
    enc.element();
    query_id.write(enc);
    enc.element();
    subject_id.write(enc);
    enc.end_list();
    enc.end_member();
    enc.member("loc", 2);
    enc.begin_list(false);
}

/// The end of a Std-seg's `loc` (`begin_std_seg`).
fn end_std_seg_locs<E: Encoder>(enc: &mut E) {
    enc.end_list();
    enc.end_member();
}

/// Whether the Std-seg of segment `seg` (of `numseg`, `gap` telling whether a segment has a
/// gap row) is kept by the mapper of `s_RemapToSubjectLoc` (a `-subject` search; see
/// `remap_strands`), which rebuilds the Std-segs and leaves out a first and a last segment
/// with fewer than two aligned rows. (The segments of a BLAST edit script start and end with
/// a substitution, so this keeps every segment of LOSAT's HSPs.)
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1568-1580
/// ```c++
///         seg_n++;
///         // Ignore starting segments with no actually aligned sequences.
///         if (seg_it == m_Segs.begin()  &&  non_gap_count < 2) {
///             continue;
///         }
///         // Ignore empty segments
///         if ( std_seg->GetLoc().empty() ) continue;
///         std_segs.push_back(std_seg);
///     }
///     if (non_gap_count < 2  &&  !std_segs.empty()) {
///         // Remove the last segment if there are no aligned sequences in it.
///         std_segs.pop_back();
///     }
/// ```
fn std_seg_kept(seg: usize, numseg: usize, gap: impl Fn(usize) -> bool) -> bool {
    !(seg == 0 && gap(0)) && !(seg + 1 == numseg && gap(seg))
}

/// The Std-segs of an HSP of a gapped translated search (`segs std`), one per segment, with
/// the remap of the query row (`RemapToQueryLoc`, `CStd_seg::OffsetRow`: the `-query_loc`
/// start added to an interval's ends) and of the subject row (`s_RemapToSubjectLoc` of a
/// `-subject` search, `remapped`: the `-subject_loc` start, on the plus strand that LOSAT's
/// ranges have; the mapper keeps each interval's strand, `CMappingRange::Map_Strand` with
/// the strand set, and `std_seg_kept`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:306-368
/// ```c++
/// static CSeq_align::C_Segs::TStd
/// s_CreateStdSegs(CRef<CSeq_id> master, CRef<CSeq_id> slave,
///                 CDense_seg::TStarts& starts,
///                 CDense_seg::TLens& lengths,
///                 CDense_seg::TStrands& strands,
///                 bool translate_master,
///                 bool translate_slave)
/// {
///     ...
///     for (int i = 0; i < nsegs; i++) {
///         CRef<CStd_seg> std_seg(new CStd_seg());
///         CRef<CSeq_loc> master_loc(new CSeq_loc());
///         CRef<CSeq_loc> slave_loc(new CSeq_loc());
///
///         std_seg->SetDim(kBlastAlignmentDim);
///
///         // Set master seqloc
///         if ( (m_start = starts[2*i]) != GAP_VALUE) {
///             master_loc->SetInt().SetId(*master);
///             master_loc->SetInt().SetFrom(m_start);
///             if (translate_master)
///                 m_stop = m_start + CODON_LENGTH*lengths[i] - 1;
///             else
///                 m_stop = m_start + lengths[i] - 1;
///             master_loc->SetInt().SetTo(m_stop);
///             master_loc->SetInt().SetStrand(strands[2*i]);
///         } else {
///             master_loc->SetEmpty(*master);
///         }
///
///         // Set slave seqloc
///         if ( (s_start = starts[2*i+1]) != GAP_VALUE) {
///             slave_loc->SetInt().SetId(*slave);
///             slave_loc->SetInt().SetFrom(s_start);
///             if (translate_slave)
///                 s_stop = s_start + CODON_LENGTH*lengths[i] - 1;
///             else
///                 s_stop = s_start + lengths[i] - 1;
///             slave_loc->SetInt().SetTo(s_stop);
///             slave_loc->SetInt().SetStrand(strands[2*i+1]);
///         } else {
///             slave_loc->SetEmpty(*slave);
///         }
///
///         std_seg->SetIds().reserve(kBlastAlignmentDim);
///         std_seg->SetLoc().reserve(kBlastAlignmentDim);
///         std_seg->SetIds().push_back(master);
///         std_seg->SetIds().push_back(slave);
///         std_seg->SetLoc().push_back(master_loc);
///         std_seg->SetLoc().push_back(slave_loc);
///
///         retval.push_back(std_seg);
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/Std_seg.cpp:170-188
/// ```c++
/// void CStd_seg::OffsetRow(TDim row,
///                          TSignedSeqPos offset)
///
/// {
///     if (offset == 0) return;
///
///     CSeq_loc& src_loc = *SetLoc()[row];
///
///     switch (src_loc.Which()) {
///     case CSeq_loc::e_Int:
///         ...
///         src_loc.SetInt().SetFrom() += offset;
///         src_loc.SetInt().SetTo() += offset;
///         break;
/// ```
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq_loc_mapper_base.cpp:312-327
/// ```c++
/// bool CMappingRange::Map_Strand(bool        is_set_strand,
///                                ENa_strand  src,
///                                ENa_strand* dst) const
/// {
///     _ASSERT(dst);
///     if ( m_Reverse ) {
///         // Always convert to reverse strand, even if the source
///         // strand is unknown.
///         *dst = Reverse(src);
///         return true;
///     }
///     if (is_set_strand) {
///         // Use original strand if set
///         *dst = src;
///         return true;
///     }
/// ```
fn write_std_segs<E: Encoder>(
    enc: &mut E,
    query_id: &SeqId<'_>,
    subject_id: &SeqId<'_>,
    segments: &Segments,
    translate: [bool; 2],
    offsets: [i32; 2],
    remapped: bool,
) {
    let numseg = segments.lens.len();
    let gap = |seg: usize| {
        segments.starts[2 * seg] == GAP_VALUE || segments.starts[2 * seg + 1] == GAP_VALUE
    };
    let ids = [query_id, subject_id];
    enc.begin_list(false);
    for seg in 0..numseg {
        if remapped && !std_seg_kept(seg, numseg, gap) {
            continue;
        }
        begin_std_seg(enc, query_id, subject_id);
        for row in 0..2 {
            enc.element();
            let start = segments.starts[2 * seg + row];
            if start == GAP_VALUE {
                write_seq_loc_empty(enc, ids[row]);
            } else {
                let len = segments.lens[seg];
                let stop = if translate[row] {
                    start + CODON_LENGTH * len - 1
                } else {
                    start + len - 1
                };
                write_seq_interval(
                    enc,
                    ids[row],
                    start + offsets[row],
                    stop + offsets[row],
                    segments.strands[2 * seg + row],
                );
            }
        }
        end_std_seg_locs(enc);
        enc.end_sequence();
    }
    enc.end_list();
}

/// The interval `(from, to)` of one row of an ungapped HSP (`offset`, `end` exclusive, in
/// the frame's letters; `length` the row's searched nucleotide length for a frame).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1348-1361
/// ```c++
///     if (hsp->query.frame == 0) {
///        query_loc->SetInt().SetFrom(hsp->query.offset);
///        query_loc->SetInt().SetTo(hsp->query.end - 1);
///     } else if (hsp->query.frame > 0) {
///        query_loc->SetInt().SetFrom(CODON_LENGTH*(hsp->query.offset) +
///                                    hsp->query.frame - 1);
///        query_loc->SetInt().SetTo(CODON_LENGTH*(hsp->query.end) +
///                                  hsp->query.frame - 2);
///     } else {
///        query_loc->SetInt().SetFrom(query_length -
///            CODON_LENGTH*(hsp->query.end) + hsp->query.frame + 1);
///        query_loc->SetInt().SetTo(query_length - CODON_LENGTH*hsp->query.offset
///                                  + hsp->query.frame);
///     }
/// ```
/// (The subject's interval, :1363-1377, is the same with `subject_length`.)
fn ungapped_interval(offset: i32, end: i32, frame: i8, length: i32) -> (i32, i32) {
    let frame = i32::from(frame);
    if frame == 0 {
        (offset, end - 1)
    } else if frame > 0 {
        (
            CODON_LENGTH * offset + frame - 1,
            CODON_LENGTH * end + frame - 2,
        )
    } else {
        (
            length - CODON_LENGTH * end + frame + 1,
            length - CODON_LENGTH * offset + frame,
        )
    }
}

/// The `diags` Seq-align of one subject's HSPs in an ungapped translated search: no `dim`, no
/// `score`, one Std-seg per HSP with the HSP's scores (the coverage over the query length
/// of the construction, a translated query's in nucleotides), and the remaps of
/// `write_std_segs` (an ungapped Std-seg has no gap row, so the mapper keeps them all).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1323-1346,1379-1386
/// ```c++
/// CRef<CStd_seg>
/// x_UngappedHSPToStdSeg(BlastHSP* hsp, CRef<CSeq_id> query_id,
///                       CRef<CSeq_id> subject_id,
///                       Int4 query_length, Int4 subject_length,
///                       const vector<string> & seqid_list)
/// {
///     CRef<CStd_seg> retval(new CStd_seg());
///
///     retval->SetDim(kBlastAlignmentDim);
///     retval->SetLoc().reserve(kBlastAlignmentDim);
///
///     CRef<CSeq_loc> query_loc(new CSeq_loc());
///     CRef<CSeq_loc> subject_loc(new CSeq_loc());
///     query_loc->SetInt().SetId(*query_id);
///     subject_loc->SetInt().SetId(*subject_id);
///
///     // Set the sequence ids
///     CStd_seg::TIds& ids = retval->SetIds();
///     ids.reserve(kBlastAlignmentDim);
///     ids.push_back(query_id);
///     ids.push_back(subject_id);
///
///     query_loc->SetInt().SetStrand(s_Frame2Strand(hsp->query.frame));
///     subject_loc->SetInt().SetStrand(s_Frame2Strand(hsp->subject.frame));
///     ...
///     retval->SetLoc().push_back(query_loc);
///     retval->SetLoc().push_back(subject_loc);
///
///     CSeq_align::TScore& score_list = retval->SetScores();
///     s_BuildScoreList(hsp, score_list, seqid_list, query_length);
///
///     return retval;
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1406-1447
/// ```c++
///     CRef<CSeq_align> seqalign(new CSeq_align());
///     ...
///     seqalign->SetType(CSeq_align::eType_diags);
///     ...
///     } else { // Translated search
///         for (index=0; index<hsp_list->hspcnt; index++) {
///             BlastHSP* hsp = hsp_array[index];
///             seqalign->SetSegs().SetStd().push_back(
///                 x_UngappedHSPToStdSeg(hsp,
///                                       query_id,
///                                       subject_id,
/// 		                              query_length,
///                                       subject_length,
///                                       emptyList));
///         }
///     }
///     s_AddUserObjectToSeqAlign(seqalign, seqid_list);
///     sa_vector.push_back(seqalign);
/// ```
/// (`s_AddUserObjectToSeqAlign` adds nothing without a Seq-id list.)
fn write_seq_align_diags<E: Encoder>(
    enc: &mut E,
    query: &QueryRecord<'_>,
    subject: &SubjectRecord<'_>,
    hsps: &[HspRecord<'_>],
) {
    enc.element();
    enc.begin_sequence("Seq-align");
    enc.member("type", 0);
    enc.enumerated("diags", 2);
    enc.end_member();
    enc.member("segs", 3);
    enc.variant("Seq-align_segs", "std", 2);
    enc.begin_list(false);
    for hsp in hsps {
        let (q_from, q_to) =
            ungapped_interval(hsp.q_offset, hsp.q_end, hsp.q_frame, query.searched_length);
        let (s_from, s_to) = ungapped_interval(
            hsp.s_offset,
            hsp.s_end,
            hsp.s_frame,
            subject.searched_length,
        );
        begin_std_seg(enc, &query.id, &subject.id);
        enc.element();
        write_seq_interval(
            enc,
            &query.id,
            q_from + query.offset,
            q_to + query.offset,
            Strand::of_frame(hsp.q_frame),
        );
        enc.element();
        write_seq_interval(
            enc,
            &subject.id,
            s_from + subject.offset,
            s_to + subject.offset,
            Strand::of_frame(hsp.s_frame),
        );
        end_std_seg_locs(enc);
        enc.member("scores", 3);
        enc.begin_list(true);
        write_score_list(enc, hsp, query.searched_length);
        enc.end_list();
        enc.end_member();
        enc.end_sequence();
    }
    enc.end_list();
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
}

/// The Seq-align of one HSP of a gapped search, or nothing for an HSP of score 0.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:627-656
/// ```c++
/// static CRef<CSeq_align>
/// s_CreateSeqAlign(CRef<CSeq_id> master, CRef<CSeq_id> slave,
///                  CDense_seg::TStarts starts,
///                  CDense_seg::TLens lengths,
///                  CDense_seg::TStrands strands,
///                  bool translate_master, bool translate_slave)
/// {
///     CRef<CSeq_align> sar(new CSeq_align());
///     sar->SetType(CSeq_align::eType_partial);
///     sar->SetDim(kBlastAlignmentDim);
///
///     if (translate_master || translate_slave) {
///         sar->SetSegs().SetStd() =
///             s_CreateStdSegs(master, slave, starts, lengths, strands,
///                             translate_master, translate_slave);
///     } else {
///         s_CreateDenseg(sar->SetSegs().SetDenseg(), master, slave, starts,
///                        lengths, strands);
///     }
///
///     return sar;
/// }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:659-690
/// ```c++
/// static CRef<CSeq_align>
/// s_BlastHSP2SeqAlign(EBlastProgramType program, BlastHSP* hsp,
///                     CRef<CSeq_id> id1, CRef<CSeq_id> id2,
///                     Int4 query_length, Int4 subject_length)
/// {
///     ...
///     if (hsp->score == 0) {
///         return CRef<CSeq_align>();
///     }
///     ...
///     translate1 = (program == eBlastTypeBlastx || program == eBlastTypeTblastx ||
///                   program == eBlastTypeRpsTblastn);
///     translate2 = (program == eBlastTypeTblastn || program == eBlastTypePsiTblastn ||
///                   program == eBlastTypeTblastx);
/// ```
/// (No LOSAT edit script has `eGapAlignDecline`, so the `disc` branch does not occur.)
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:22-40
/// ```text
/// Seq-align ::= SEQUENCE {
///     type ENUMERATED {
///         not-set (0) ,
///         global (1) ,
///         diags (2) ,     -- unbroken, but not ordered, diagonals
///         partial (3) ,   -- mapping pieces together
///         disc (4) ,      -- discontinuous alignment
///         other (255) } ,
///     dim INTEGER OPTIONAL ,     -- dimensionality
///     score SET OF Score OPTIONAL ,   -- for whole alignment
///     segs CHOICE {                   -- alignment data
///         dendiag SEQUENCE OF Dense-diag ,
///         denseg              Dense-seg ,
///         std     SEQUENCE OF Std-seg ,
/// ```
fn write_seq_align_gapped<E: Encoder>(
    enc: &mut E,
    run: &RunRecord<'_>,
    query: &QueryRecord<'_>,
    subject: &SubjectRecord<'_>,
    hsp: &HspRecord<'_>,
    segments: &mut Segments,
) {
    let translate1 = query_is_translated(run.program);
    let translate2 = subject_is_translated(run.program);
    collect_seq_align_data(
        hsp,
        query.searched_length,
        subject.searched_length,
        translate1,
        translate2,
        segments,
    );
    let translated = translate1 || translate2;
    // The mapper of `s_RemapToSubjectLoc` rebuilds a Dense-seg's strands (`remap_strands`);
    // it keeps a Std-seg's (`write_std_segs`).
    if run.remaps_subjects() && !translated {
        remap_strands(segments);
    }
    enc.element();
    enc.begin_sequence("Seq-align");
    enc.member("type", 0);
    enc.enumerated("partial", 3);
    enc.end_member();
    enc.member("dim", 1);
    enc.int(BLAST_ALIGNMENT_DIM);
    enc.end_member();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1252-1260
    // ```c++
    // static void
    // s_AddScoresToSeqAlign(CRef<CSeq_align>  & seqalign,
    //                       const BlastHSP    * hsp,
    //                       const vector<string> & seqid_list,
    //                       Int4 query_length)
    // {
    //     // Add the scores for this HSP
    //     CSeq_align::TScore& score_list = seqalign->SetScore();
    //     s_BuildScoreList(hsp, score_list, seqid_list, query_length);
    // }
    // ```
    enc.member("score", 2);
    enc.begin_list(true);
    write_score_list(enc, hsp, coverage_query_length(run, query.searched_length));
    enc.end_list();
    enc.end_member();
    enc.member("segs", 3);
    if translated {
        enc.variant("Seq-align_segs", "std", 2);
        write_std_segs(
            enc,
            &query.id,
            &subject.id,
            segments,
            [translate1, translate2],
            [query.offset, subject.offset],
            run.remaps_subjects(),
        );
    } else {
        enc.variant("Seq-align_segs", "denseg", 1);
        write_dense_seg(
            enc,
            &query.id,
            &subject.id,
            segments,
            [query.offset, subject.offset],
        );
    }
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
}

/// Whether `write_query` writes the Seq-aligns of the run's program and mode: every gapped
/// search, and the ungapped translated searches (`BLASTUngappedHspListToSeqAlign` makes the
/// Std-segs of `write_seq_align_diags` for them, a Dense-diag for an ungapped BLASTN or
/// BLASTP search, which LOSAT does not run).
pub fn writes(run: &RunRecord<'_>) -> bool {
    run.gapped || query_is_translated(run.program) || subject_is_translated(run.program)
}

/// Whether the HSPs of a query make a Seq-align (`results.HasAlignments()`): in a gapped
/// search an HSP of score 0 makes none (`s_BlastHSP2SeqAlign`); an ungapped search makes one
/// Seq-align of every subject's HSPs.
fn has_alignments(run: &RunRecord<'_>, hsps: &[HspRecord<'_>]) -> bool {
    if run.gapped {
        hsps.iter().any(|hsp| hsp.score != 0)
    } else {
        !hsps.is_empty()
    }
}

/// The Seq-annot of one query's Seq-aligns: the `desc` user objects and `data align` with,
/// subject by subject in the order of `hsps`, one Seq-align per HSP of score other than 0
/// (a gapped search) or one `diags` Seq-align per subject (an ungapped search).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1606-1626
/// ```c++
///         GetFilteredRedundantSeqids(*seqinfo_src, hsp_list->oid, seqid_list, subject_id->IsGi());
///         // stores a CSeq_align for each matching sequence
///         vector<CRef<CSeq_align > > hit_align;
///         if (is_gapped) {
///                 BLASTHspListToSeqAlign(prog,
///                                        hsp_list,
///                                        query_id,
///                                        subject_id,
///                                        query_length,
///                                        subj_length,
///                                        is_ooframe,
///                                        seqid_list,
///                                        hit_align);
///         } else {
///                 BLASTUngappedHspListToSeqAlign(prog,
///                                                hsp_list,
///                                                query_id,
///                                                subject_id,
///                                                query_length,
///                                                subj_length,
///                                                seqid_list,
/// ```
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:460-476
/// ```text
/// Seq-annot ::= SEQUENCE {
///     id SET OF Annot-id OPTIONAL ,
///     db INTEGER {                 -- source of annotation
///         ...
///         other (255) } OPTIONAL ,
///     name VisibleString OPTIONAL ,-- source if "other" above
///     desc Annot-descr OPTIONAL ,  -- used only for stand alone Seq-annots
///     data CHOICE {
///         ftable SET OF Seq-feat ,
///         align SET OF Seq-align ,
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:381-385
/// ```c++
///    	//Fill in data -- Seq align
/// 	retval->SetData().SetAlign();
///     ITERATE(CSeq_align_set::Tdata, itr, alnset.Get()) {
///         retval->SetData().SetAlign().push_back(*itr);
///     }
/// ```
fn write_seq_annot<E: Encoder>(
    enc: &mut E,
    run: &RunRecord<'_>,
    query: &QueryRecord<'_>,
    subjects: &[SubjectRecord<'_>],
    hsps: &[HspRecord<'_>],
) {
    let mut segments = Segments::default();
    enc.begin_document("Seq-annot");
    enc.begin_sequence("Seq-annot");
    write_desc(enc, run);
    enc.member("data", 4);
    enc.variant("Seq-annot_data", "align", 1);
    enc.begin_list(true);
    if run.gapped {
        for hsp in hsps.iter().filter(|hsp| hsp.score != 0) {
            let subject = &subjects[hsp.subject as usize];
            write_seq_align_gapped(enc, run, query, subject, hsp, &mut segments);
        }
    } else {
        for list in hsps.chunk_by(|a, b| a.subject == b.subject) {
            write_seq_align_diags(enc, query, &subjects[list[0].subject as usize], list);
        }
    }
    enc.end_list();
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
    enc.end_document();
}

/// An error of `write_query`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeqAlignError {
    /// NCBI's object stream raised an error (the document is not written).
    Serial(SerialError),
    /// A program or mode whose Seq-aligns this writer does not make yet (`writes`).
    Unsupported(Program),
}

impl std::fmt::Display for SeqAlignError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SeqAlignError::Serial(error) => error.fmt(f),
            SeqAlignError::Unsupported(program) => write!(
                f,
                "the Seq-align output (outfmt 8, 12) is not supported by LOSAT's {}",
                program.name()
            ),
        }
    }
}

impl std::error::Error for SeqAlignError {}

/// Appends to `buf` the document of one query (`format`): the Seq-annot of its Seq-aligns,
/// from `hsps` (in NCBI's Seq-align order: subject by subject in hit-list order, each
/// subject's HSPs sorted by e-value) and the `subjects` they index. Returns false and
/// appends nothing when the query has no Seq-align (`results.HasAlignments()` is false:
/// NCBI writes nothing for it). On an error, `buf` is left as it was.
pub fn write_query(
    format: SerialFormat,
    run: &RunRecord<'_>,
    query: &QueryRecord<'_>,
    subjects: &[SubjectRecord<'_>],
    hsps: &[HspRecord<'_>],
    buf: &mut Vec<u8>,
) -> Result<bool, SeqAlignError> {
    if !writes(run) {
        return Err(SeqAlignError::Unsupported(run.program));
    }
    if !has_alignments(run, hsps) {
        return Ok(false);
    }
    let start = buf.len();
    let error = match format {
        SerialFormat::AsnText => {
            let mut enc = AsnTextEncoder::new(buf);
            write_seq_annot(&mut enc, run, query, subjects, hsps);
            enc.error()
        }
        SerialFormat::AsnBinary => {
            let mut enc = AsnBinaryEncoder::new(buf);
            write_seq_annot(&mut enc, run, query, subjects, hsps);
            enc.error()
        }
        SerialFormat::Json => {
            let mut enc = JsonEncoder::new(buf);
            write_seq_annot(&mut enc, run, query, subjects, hsps);
            enc.error()
        }
    };
    if let Some(error) = error {
        buf.truncate(start);
        return Err(SeqAlignError::Serial(error));
    }
    Ok(true)
}
