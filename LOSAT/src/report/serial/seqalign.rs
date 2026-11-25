//! The Seq-annot of Seq-aligns that NCBI writes for each query in outfmt 8 (ASN.1 text) and
//! 12 (JSON): NCBI's conversion of a query's final hit list to Seq-aligns
//! (algo/blast/api/blast_seqalign.cpp) and its wrapping in a Seq-annot
//! (algo/blast/format/blastfmtutil.cpp), written through an `Encoder` in the member order of
//! the ASN.1 modules.
//!
//! The documents are made by `write_query`, one per query with alignments, from the report
//! records of `report::record`; `write_query_to` and `write_hit_list_seqaligns` are the
//! per-query write step and hit-list loop the programs' writers share. The segment computation is that of `s_CollectSeqAlignData`
//! for every program; the forms are NCBI's:
//! - a gapped search whose rows are not translated (BLASTN, BLASTP): one Seq-align per HSP
//!   with a Dense-seg (`write_dense_seg`, `s_CreateDenseg`);
//! - a gapped search with a translated row (BLASTX, TBLASTN): one Seq-align per HSP with one
//!   Std-seg per segment (`write_std_segs`, `s_CreateStdSegs`): Seq-intervals in nucleotides
//!   on the translated row, `empty` for a gap row;
//! - an ungapped translated search (TBLASTX, BLASTX `-ungapped`): one `diags` Seq-align per
//!   subject, one Std-seg with its scores per HSP (`write_seq_align_diags`,
//!   `BLASTUngappedHspListToSeqAlign` / `x_UngappedHSPToStdSeg`).
//!
//! - an ungapped search whose rows are not translated (BLASTN `-ungapped`, BLASTP `-ungapped`): one
//!   `diags` Seq-align per subject, one Dense-diag with its scores per HSP
//!   (`write_seq_align_dendiag`, `BLASTUngappedHspListToSeqAlign` / `x_UngappedHSPToDenseDiag`).
//!
//! A `-subject` search's Seq-aligns go through NCBI's mapper of `s_RemapToSubjectLoc`
//! (`remap_strands` for a Dense-seg; for Std-segs `mapped_script`, which BLASTX and TBLASTN apply
//! to their report records once, so that every format reads the mapped HSP).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:676-707
//! `PruneSeqalign` keeps the Seq-aligns of the first `m_HitlistSize` subjects, which is the
//! search's own hit list size for these formats, so it keeps every Seq-align of the search's
//! hit lists (outfmt-20261019 S/DESIGN_OF3.md §2.1).

use std::io::Write;

use crate::blastinput::app::Program;
use crate::common::GapEditOp;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::ostream::{ReportWrite, Site};
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::{HspRecord, QueryNames, QueryRecord, RunRecord, SubjectRecord};
use crate::report::seqid::SeqId;

use super::asn_ber::AsnBinaryEncoder;
use super::asn_text::AsnTextEncoder;
use super::encoder::{Encoder, ReplacedChar, SerialError, Traced};
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

    /// The `-outfmt` number of the format: 8, 9 or 12.
    pub fn outfmt(self) -> i32 {
        match self {
            SerialFormat::AsnText => 8,
            SerialFormat::AsnBinary => 9,
            SerialFormat::Json => 12,
        }
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:60-71
const SMALLEST_EVALUE: f64 = 1.0e-180;

/// The e-value of a Seq-align score list: 0 below `SMALLEST_EVALUE` (shared with BLAST XML,
/// which reads the same Seq-align, `s_SeqAlignSetToXMLHsps`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1186
pub(crate) fn seqalign_evalue(evalue: f64) -> f64 {
    if evalue < SMALLEST_EVALUE {
        0.0
    } else {
        evalue
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1199-1203, 1490-1496
/// The `num_ident` score of an HSP's Seq-align (`s_AddScoresToSeqAlign`), from the engine's
/// identities `num_ident` of a search of `program_type` (`RunRecord::program_type`): the count,
/// none where it is below 0, and none for psiblast's (`ProgramType::PsiBlast`, every round) count
/// of 0, which NCBI sets to -1 before it adds the scores of a gapped HSP (psiblast's searches are
/// gapped; PSI-TBLASTN's keep their 0). The one owner of that rule: the Seq-annot formats and the
/// archive write the score (`write_seq_align_gapped`), BLAST XML's `catch` reads it back
/// (`read_num_ident`).
pub(crate) fn seq_align_num_ident(
    num_ident: i64,
    program_type: Option<crate::report::record::ProgramType>,
) -> Option<i32> {
    let psiblast = program_type == Some(crate::report::record::ProgramType::PsiBlast);
    if num_ident < 0 || (psiblast && num_ident == 0) {
        None
    } else {
        Some(i32::try_from(num_ident).unwrap_or(i32::MAX))
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:810-849
/// The `num_ident` of an HSP's Seq-align as `GetAlnScores` reads it back (BLAST XML's `catch` of
/// an HSP whose subject's lookup throws, `blast_xml::XmlHsp::alignment_identities`): the score
/// that the Seq-align has (`seq_align_num_ident`), or -1 where it has none.
pub(crate) fn read_num_ident(
    num_ident: i64,
    program_type: Option<crate::report::record::ProgramType>,
) -> i32 {
    seq_align_num_ident(num_ident, program_type).unwrap_or(-1)
}

pub(crate) const GAP_VALUE: i32 = -1;
const BLAST_ALIGNMENT_DIM: i32 = 2;

/// `CODON_LENGTH` (include/algo/blast/core/blast_def.h:62-63: `#define CODON_LENGTH 3`).
const CODON_LENGTH: i32 = 3;

/// The values of `Na-strand` that BLAST writes.
///
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:115-121
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
    fn of_frame(frame: i8) -> Self {
        if frame > 0 {
            Strand::Plus
        } else if frame < 0 {
            Strand::Minus
        } else {
            Strand::Unknown
        }
    }

    /// The strand after `CDense_seg::Reverse` (outfmt 0's minus-strand query row,
    /// `report::query_anchored::HspAlnVec::set`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:647-654
    pub(crate) fn reversed(self) -> Self {
        match self {
            Strand::Plus => Strand::Minus,
            Strand::Minus => Strand::Plus,
            Strand::Unknown => Strand::Unknown,
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

impl Segments {
    /// The numbers of the segments' rows: `starts` (2 per segment, query then subject,
    /// `GAP_VALUE` for a gap; in the searched letters, before the remaps of the query's and
    /// the subject's ranges), `lens` and `strands` (2 per segment).
    pub(crate) fn starts(&self) -> &[i32] {
        &self.starts
    }

    pub(crate) fn lens(&self) -> &[i32] {
        &self.lens
    }

    pub(crate) fn strands(&self) -> &[Strand] {
        &self.strands
    }
}

/// The Dense-seg segments of a gapped HSP of a non-translated search as the Seq-align holds
/// them after `RemapToQueryLoc` and `s_RemapToSubjectLoc` (for a `-subject` search), minus the
/// range starts that the caller adds (`QueryRecord::offset`, `SubjectRecord::offset`): what
/// `write_seq_align_gapped` writes into the Dense-seg. The SAM formatter reads them
/// (`report::sam`).
pub(crate) fn gapped_segments(
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
    if run.remaps_subjects() && !(translate1 || translate2) {
        remap_strands(segments);
    }
}

/// The Dense-seg segments of a hit of a protein search without translation (BLASTP) whose
/// coordinates are in its records (`Hit::q_start` and `s_start`, 1-based, after the shift to the
/// records): `collect_seq_align_data` of the HSP at those coordinates. NCBI makes the segments in
/// the searched letters and `RemapToQueryLoc` and `s_RemapToSubjectLoc` add the starts of the
/// ranges to the rows: the same starts, as a protein row has no strand to remap (the query-anchored
/// reports, `report::query_anchored`).
pub(crate) fn protein_row_segments(hit: &crate::common::Hit, segments: &mut Segments) {
    let record = HspRecord {
        subject: 0,
        q_offset: hit.q_start as i32 - 1,
        q_end: hit.q_end as i32,
        q_left_out: 0,
        s_offset: hit.s_start as i32 - 1,
        s_end: hit.s_end as i32,
        q_frame: 0,
        s_frame: 0,
        score: hit.raw_score,
        num: 0,
        evalue: hit.e_value,
        bit_score: hit.bit_score,
        num_ident: hit.num_ident as i32,
        num_positives: hit.num_positives as i32,
        comp_adjust: 0,
        ops: hit.gap_info.as_deref().unwrap_or(&[]),
    };
    // A protein row's strand is unknown: the lengths are not read.
    collect_seq_align_data(&record, 0, 0, false, false, segments);
}

/// The Dense-seg segments of a BLASTN hit whose coordinates are in its records (`Hit::q_start` to
/// `q_end` on the query's plus strand; `s_start` above `s_end` for a minus-strand HSP), and whether
/// its rows' strands differ (a minus-strand HSP: the query row on the minus strand, the subject's
/// on the plus strand). NCBI counts a minus-strand query row's starts from the end of the searched
/// query (`s_GetAlignmentStart`: `length - curr_pos - num`, the HSP's offset in its minus-strand
/// context) and `RemapToQueryLoc` moves them by the start of the `-query_loc` interval: the same
/// starts as counted from the HSP's end in the record, which `collect_seq_align_data` gives with
/// that end as the length and an offset of 0. The plus-strand subject row's starts are the
/// record's (`s_RemapToSubjectLoc` adds the `-subject_loc` start). The query-anchored reports read
/// them (`report::query_anchored`).
pub(crate) fn nucleotide_row_segments(
    hit: &crate::common::Hit,
    segments: &mut Segments,
) -> crate::report::query_anchored::HspStrands {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:114-140
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1520-1534
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1541-1547
    let minus = hit.query_frame < 0;
    let (q_offset, q_end, query_length) = if minus {
        (0, (hit.q_end + 1 - hit.q_start) as i32, hit.q_end as i32)
    } else {
        (hit.q_start as i32 - 1, hit.q_end as i32, 0)
    };
    let (s_low, s_high) = (hit.s_start.min(hit.s_end), hit.s_start.max(hit.s_end));
    let record = HspRecord {
        subject: 0,
        q_offset,
        q_end,
        q_left_out: 0,
        s_offset: s_low as i32 - 1,
        s_end: s_high as i32,
        q_frame: if minus { -1 } else { 1 },
        s_frame: 1,
        score: hit.raw_score,
        num: 0,
        evalue: hit.e_value,
        bit_score: hit.bit_score,
        num_ident: hit.num_ident as i32,
        num_positives: hit.num_positives as i32,
        comp_adjust: 0,
        ops: hit.gap_info.as_deref().unwrap_or(&[]),
    };
    collect_seq_align_data(&record, query_length, 0, false, false, segments);
    crate::report::query_anchored::HspStrands {
        query_minus: minus,
        ..Default::default()
    }
}

/// The Dense-seg segments of a hit of a translated search (`CreateDensegFromStdseg` of its
/// Std-segs: the starts of `s_CollectSeqAlignData`, a translated row's in nucleotides, and the
/// lengths in residues) whose coordinates are in its records, and its rows' strands and query
/// frame. A translated row is given by its frame in its record and the record's length in
/// nucleotides (`seq_range::record_frame`, the frame the report shows): the query's of BLASTX and
/// TBLASTX (`query`), the subject's of TBLASTN and TBLASTX (`subject`); an untranslated row is a
/// protein's (`None`), whose positions are the hit's. NCBI makes them in the searched letters and
/// the frame of the search and `RemapToQueryLoc` / `s_RemapToSubjectLoc` add the starts of the
/// `-query_loc` / `-subject_loc` intervals: the same nucleotides as `s_GetAlignmentStart` counts
/// in the record's frame (on the minus strand from the record's end, `original_length - 3 *
/// (curr_pos + num) + frame + 1`), where the HSP's offset in that frame's residues is exact (the
/// frame is that of the row's first nucleotide). A row of a minus frame is on the minus strand
/// (`s_Frame2Strand`); the query frame is the query's group (`query_anchored::Hsps`).
pub(crate) fn translated_row_segments(
    hit: &crate::common::Hit,
    query: Option<(i8, usize)>,
    subject: Option<(i8, usize)>,
    segments: &mut Segments,
) -> crate::report::query_anchored::HspStrands {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:114-140
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1520-1534
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1541-1547
    // A row's offset and end in its frame's residues (`HspRecord::q_offset`, `q_end`): the
    // residue of the HSP's first codon is `(start - (frame - 1)) / 3` on the plus strand, where
    // `frame - 1 + 3 * offset` is its first nucleotide, and `(length + frame + 1 - high) / 3` on
    // the minus strand, where `length - 3 * offset + frame + 1` is its end; an untranslated row's
    // positions are its residues.
    let row = |start: usize, end: usize, translated: Option<(i8, usize)>| {
        let (low, high) = (start.min(end) as i32, start.max(end) as i32);
        match translated {
            Some((frame, length)) => {
                let (frame, length) = (i32::from(frame), length as i32);
                let offset = if frame > 0 {
                    (low - 1 - (frame - 1)) / CODON_LENGTH
                } else {
                    (length + frame + 1 - high) / CODON_LENGTH
                };
                (
                    offset,
                    offset + (high - low + 1) / CODON_LENGTH,
                    frame as i8,
                    length,
                )
            }
            None => (low - 1, high, 0, 0),
        }
    };
    let (q_offset, q_end, q_frame, query_length) = row(hit.q_start, hit.q_end, query);
    let (s_offset, s_end, s_frame, subject_length) = row(hit.s_start, hit.s_end, subject);
    let record = HspRecord {
        subject: 0,
        q_offset,
        q_end,
        q_left_out: 0,
        s_offset,
        s_end,
        q_frame,
        s_frame,
        score: hit.raw_score,
        num: 0,
        evalue: hit.e_value,
        bit_score: hit.bit_score,
        num_ident: hit.num_ident as i32,
        num_positives: hit.num_positives as i32,
        comp_adjust: 0,
        ops: hit.gap_info.as_deref().unwrap_or(&[]),
    };
    collect_seq_align_data(
        &record,
        query_length,
        subject_length,
        query.is_some(),
        subject.is_some(),
        segments,
    );
    dense_seg_strands(segments);
    // A query row that starts with a gap (`-db`'s; `-subject`'s mapping drops such a first
    // segment, `mapped_script`) gives its first Std-seg an empty query interval, whose frame is
    // the record's length's (`HspStrands::master_frame`).
    let leading_gap = query.is_some() && row_end_gaps(record.ops, AlignRow::Query).leading > 0;
    crate::report::query_anchored::HspStrands {
        query_minus: q_frame < 0,
        subject_minus: s_frame < 0,
        query_frame: crate::report::query_anchored::HspStrands::master_frame(
            q_frame,
            leading_gap,
            query_length as usize,
        ),
        no_widths: false,
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:820-842
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:937-942
/// The strands of the Dense-seg that `CreateDensegFromStdseg` makes of the Std-segs `segments`
/// (`collect_seq_align_data`): each row takes, in every segment, the strand of its first
/// segment with one; a first segment's gap row has none (`Strand::Unknown`, the Std-seg's empty
/// location), so a Smith-Waterman alignment of `-use_sw_tback` that starts with a gap in a minus
/// row is shown on the minus strand. Each X-drop row has one strand already.
fn dense_seg_strands(segments: &mut Segments) {
    for row in 0..2 {
        let mut row_strands = segments.strands.iter().skip(row).step_by(2);
        let Some(&strand) = row_strands.find(|&&strand| strand != Strand::Unknown) else {
            continue;
        };
        for known in segments.strands.iter_mut().skip(row).step_by(2) {
            *known = strand;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:779-791
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:909-935
/// The Dense-seg that `CreateDensegFromStdseg` makes of the empty rows of a BLASTX HSP of length 0
/// (`ungapped_rows_empty`) as the engine made them (`-db`, and `-subject`'s Seq-aligns that the
/// Seq-loc mapper leaves, `mapper_rewrites_empty_rows`): each row `(from, from - 1)` has the length
/// `to - from + 1` of 0 in the unsigned `TSeqPos`, so its one segment has length 0 and starts at
/// the rows' `from`, and no width is determined (both are 0), which resets the widths
/// (`HspStrands::no_widths`). The query row keeps the strand of its frame `query_frame`, the protein
/// row has none (`write_empty_rows`). The rows' `from` in the records are the lower of each row's
/// shown ends (the query's `(from + 1, from)` on a plus frame, `(from, from + 1)` on a minus frame,
/// `blastx::report::query_endpoints`; the subject's `(from + 1, from)`), which hold the
/// `-query_loc` and `-subject_loc` starts that `RemapToQueryLoc` and the mapper add.
pub(crate) fn empty_row_segments(
    hit: &crate::common::Hit,
    query_frame: i8,
    segments: &mut Segments,
) -> crate::report::query_anchored::HspStrands {
    let q_from = hit.q_start.min(hit.q_end);
    let s_from = hit.s_start.min(hit.s_end);
    segments.starts.clear();
    segments.lens.clear();
    segments.strands.clear();
    segments.starts.extend([q_from as i32, s_from as i32]);
    segments.lens.push(0);
    segments
        .strands
        .extend([Strand::of_frame(query_frame), Strand::Unknown]);
    crate::report::query_anchored::HspStrands {
        query_minus: query_frame < 0,
        subject_minus: false,
        query_frame,
        no_widths: true,
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:89-100
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
fn get_protein_frame_length(nuc_length: i32, frame: i8) -> i32 {
    (nuc_length - (i32::from(frame).abs() - 1) % CODON_LENGTH) / CODON_LENGTH
}

/// The segments of `hsp` (`s_CollectSeqAlignData` over the whole edit script): one per
/// edit-script operation, in script order.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:168-268
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

/// The "Blast Type" of the run's program: `EProgramToTaskName` of `ProgramNameToEnum` of the
/// program name and the `EProgram` value (`RunRecord::program_type` where the program does not
/// name it: psiblast, PSI-TBLASTN, rpsblast, rpstblastn). The program name of every BLASTN task
/// is "blastn" (`Blast_ProgramNameFromType`, blast_format.cpp:94), and of every BLASTP task
/// "blastp".
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_types.hpp:56-68
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:676-694
fn blast_type(run: &RunRecord<'_>) -> (&'static str, i32) {
    use crate::report::record::ProgramType;
    match run.program_type {
        Some(program_type) => (
            program_type.name(),
            match program_type {
                ProgramType::RpsBlast => 6,
                ProgramType::RpsTblastn => 7,
                ProgramType::PsiBlast => 10,
                ProgramType::PsiTblastn => 11,
            },
        ),
        None => match run.program {
            Program::Blastn => ("blastn", 1),
            Program::Blastp => ("blastp", 2),
            Program::Blastx => ("blastx", 3),
            Program::Tblastn => ("tblastn", 4),
            Program::Tblastx => ("tblastx", 5),
        },
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
/// (`data`: `str` 0, `int` 1, `real` 2, `bool` 3, general.asn:136-140.)
fn write_user_object<E: Encoder>(enc: &mut E, type_str: &[u8], label: &[u8], data: UserData) {
    // `Annotdesc` choice variant `user` (`User-object`, index 4;
    // objects/seq/seq.asn:440-445).
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:330-381
/// A `-subject` search has no database name (`n/a`). A `-db` named `n/a` or `SRA` is
/// refused by the programs (`check_database`), so its title is written (a byte that NCBI's
/// ASN.1 text or binary stream replaces is written as `#` with NCBI's post, `write_query_to`).
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
    let (task, program) = blast_type(run);
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

/// LOSAT's refusal of a database search (`-db name`) for a Seq-annot format, made before the
/// search: a database named `n/a` or `SRA` (compared without case), for which NCBI writes another
/// third user object (`CreateSeqAnnotFromSeqAlignSet` above). The title of any other database is
/// written as NCBI writes it, also with a byte that NCBI's stream replaces (`write_query_to`).
pub fn check_database(format: SerialFormat, name: &str, program: Program) -> anyhow::Result<()> {
    if name.eq_ignore_ascii_case("n/a") || name.eq_ignore_ascii_case("SRA") {
        anyhow::bail!(
            "a database named \"{name}\" (which NCBI BLAST+ describes as no database or as SRA in output format {}) is not supported by LOSAT's {}",
            format.outfmt(),
            program.name()
        );
    }
    Ok(())
}

/// A Score with an integer value.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1098-1110
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:283-289
fn write_score_int<E: Encoder>(enc: &mut E, name: &str, value: i32) {
    write_score_id(enc, name);
    enc.member("value", 1);
    enc.variant("Score_value", "int", 1);
    enc.int(value);
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
}

/// A Score with a real value (`s_MakeScore` with `false`, see `write_score_int`).
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
/// Every caller passes an empty `seqid_list` (`emptyList` in the callers below; the list goes to
/// `write_use_this_seqid`).
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
            query_coverage(hsp.q_offset, hsp.q_end + hsp.q_left_out, query_length),
        );
    }
}

/// The Dense-seg of an HSP of a non-translated search, with the remap of the query row
/// (`RemapToQueryLoc`: the `-query_loc` start added to every start that is not a gap) and
/// of the subject row (`s_RemapToSubjectLoc` of a `-subject` search: the `-subject_loc`
/// start, on the plus strand that LOSAT's ranges have).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:271-301
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:890-918
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:74-81
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
        enc.begin_named("Na-strand");
        enc.enumerated(strand.name(), strand.value());
        enc.end_named();
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
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1242-1269
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1397-1402
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
    enc.begin_named("Na-strand");
    enc.enumerated(strand.name(), strand.value());
    enc.end_named();
    enc.end_member();
    enc.member("id", 3);
    id.write(enc);
    enc.end_member();
    enc.end_sequence();
    enc.end_variant();
}

/// A Seq-loc `empty`: the gap row of a Std-seg (`master_loc`'s `SetEmpty` of `master`, see
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
/// `remap_strands`), which rebuilds the Std-segs and leaves out the first and the last segment
/// when it has fewer than two aligned rows (once each: a segment that becomes first or last
/// stays). An X-drop edit script starts and ends with a substitution, so this keeps its every
/// segment; a Smith-Waterman alignment of `-use_sw_tback` can start with a gap (line PS item
/// PS12d, line RG item X4b-gap), whose first segment the mapper leaves out.
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1568-1580
fn std_seg_kept(seg: usize, numseg: usize, gap: impl Fn(usize) -> bool) -> bool {
    (seg != 0 || !gap(0)) && (seg + 1 != numseg || !gap(seg))
}

/// What `-subject`'s Seq-loc mapper keeps of a gapped translated HSP whose edit script is `ops`
/// (one Std-seg per operation, `collect_seq_align_data`; `std_seg_kept`): the operations `kept`,
/// and the letters of the operations it leaves out before them (`start`) and after them (`end`),
/// each `[query, subject]` in the rows' letters (a `Del` has the subject's, an `Ins` the
/// query's). The mapped HSP starts `start` letters later in each row and ends `end` letters
/// earlier; its scores are the HSP's, its coverage too (`MappedScript::query_left_out`). BLASTX
/// and TBLASTN apply it once to the records of a `-subject` report
/// (`algorithm::blastx::runtime::BatchResults::put_in_report_order`,
/// `algorithm::tblastn::stage_d_results::KappaResultHitList::put_in_report_order`, which refuse
/// a second mapping), so that every format, the Seq-aligns (`write_std_segs`) among them, reads
/// the mapped HSP. An empty script (an ungapped HSP) keeps everything.
pub(crate) fn mapped_script(ops: &[GapEditOp]) -> MappedScript {
    let numseg = ops.len();
    let gap = |seg: usize| !ops[seg].is_sub();
    let first = usize::from(numseg > 0 && !std_seg_kept(0, numseg, gap));
    let last = if numseg > first && !std_seg_kept(numseg - 1, numseg, gap) {
        numseg - 1
    } else {
        numseg
    };
    let letters = |dropped: &[GapEditOp]| {
        dropped.iter().fold([0u32; 2], |[q, s], op| match *op {
            GapEditOp::Sub(n) => [q + n, s + n],
            GapEditOp::Ins(n) => [q + n, s],
            GapEditOp::Del(n) => [q, s + n],
        })
    };
    MappedScript {
        kept: first..last,
        start: letters(&ops[..first]),
        end: letters(&ops[last..]),
    }
}

/// The part of an edit script that `-subject`'s Seq-loc mapper keeps (`mapped_script`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct MappedScript {
    /// The kept operations' indices.
    pub(crate) kept: std::ops::Range<usize>,
    /// The `[query, subject]` letters of the operations left out before `kept`.
    pub(crate) start: [u32; 2],
    /// The `[query, subject]` letters of the operations left out after `kept`.
    pub(crate) end: [u32; 2],
}

impl MappedScript {
    /// Whether the mapper leaves out an operation.
    pub(crate) fn drops(&self, ops: &[GapEditOp]) -> bool {
        self.kept != (0..ops.len())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1221-1225
    // NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:184-189
    /// The query letters that the mapper leaves out (of an end `Ins`), which the HSP's
    /// `hsp_percent_coverage` (`query_coverage`, the tables' `qcovhsp`) still counts: NCBI
    /// computes that score on the BlastHSP before `s_RemapToSubjectLoc`, whose mapper keeps the
    /// Seq-align's scores (`HspRecord::q_left_out`).
    pub(crate) fn query_left_out(&self) -> u32 {
        self.start[0] + self.end[0]
    }
}

/// The Std-segs of an HSP of a gapped translated search (`segs std`), one per segment, with
/// the remap of the query row (`RemapToQueryLoc`, `CStd_seg::OffsetRow`: the `-query_loc`
/// start added to an interval's ends) and of the subject row (`s_RemapToSubjectLoc` of a
/// `-subject` search: the `-subject_loc` start of `offsets`, on the plus strand that LOSAT's
/// ranges have; the mapper keeps each interval's strand, `CMappingRange::Map_Strand` with
/// the strand set). The records of a `-subject` report come mapped (`mapped_script`, which
/// leaves out the end segments that the mapper drops): every segment is written.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:306-368
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/Std_seg.cpp:170-189
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq_loc_mapper_base.cpp:312-327
fn write_std_segs<E: Encoder>(
    enc: &mut E,
    query_id: &SeqId<'_>,
    subject_id: &SeqId<'_>,
    segments: &Segments,
    translate: [bool; 2],
    offsets: [i32; 2],
) {
    let numseg = segments.lens.len();
    let ids = [query_id, subject_id];
    enc.begin_list(false);
    for seg in 0..numseg {
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
/// (The subject's interval, :1363-1377, is the same with `subject_length`.)
pub(crate) fn ungapped_interval(offset: i32, end: i32, frame: i8, length: i32) -> (i32, i32) {
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

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:226-230
/// Whether the Std-seg rows of an ungapped HSP from `offset` to `end` (exclusive, in its frame's
/// letters) are empty: an HSP of length 0, whose rows `ungapped_interval` ends one before they
/// start (a `CRange` from `from` to one less than `from` is empty, so `CSeq_align::GetSeqRange`
/// throws on them). BLASTX's one-hit word finder with `-window_size 0` reports such HSPs
/// (`algorithm::blastx::empty_rows`); every other reported HSP has rows of its letters. The one
/// predicate of every writer (`HspRecord::rows_empty`, `PairwiseHit::rows_empty`,
/// `empty_rows::rows_empty`).
#[inline]
pub(crate) fn ungapped_rows_empty(offset: i64, end: i64) -> bool {
    end == offset
}

/// A row of a pairwise Seq-align (`CSeq_align::TDim`): the query's, 0, and the subject's, 1.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum AlignRow {
    Query,
    Subject,
}

/// Why `CSeq_align::GetSeqRange` throws on a row of a translated search's Seq-align (`segs std`,
/// one Std-seg per segment of an HSP's edit script, `write_std_segs`): the exceptions of
/// `Seq_align.cpp`, which `row_range_exception` makes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RowRangeError {
    /// `Row is empty` (Seq_align.cpp:227): no segment has a location of letters in the row, the
    /// rows of an ungapped HSP of length 0 (`ungapped_rows_empty`).
    Empty,
    /// `Row seqids not consistent. Cannot determine range.` (Seq_align.cpp:186): the row's first
    /// Std-seg has no location of letters (`empty`, the gap row of a segment), so the row takes no
    /// Seq-id there, and the first segment with an interval compares its Seq-id with that unset
    /// one. A gapped HSP whose edit script starts with a gap in the row: an alignment of the
    /// Smith-Waterman redo (`-use_sw_tback`), which `ALIGN_EX` runs from the Smith-Waterman start;
    /// an X-drop alignment starts with a substitution. A gap at the row's end is read as any.
    SeqidsNotConsistent,
}

/// The `__FILE__` of `Seq_align.cpp` in BLAST+'s build (the `seq` library compiles it from
/// `src/objects/seq`), as NCBI's reports name it.
const SEQ_ALIGN_CPP: &str = "src/objects/seq/../seqalign/Seq_align.cpp";

/// `CSeq_align::GetSeqRange`'s exception `error` (`RowRangeError`) where a formatter asks for the
/// range of such a row: `Row is empty` (Seq_align.cpp:227; the deflines of outfmt 0, `qcovs`, BLAST
/// XML with `-db`), `Row seqids not consistent. Cannot determine range.` (Seq_align.cpp:186; the
/// same readers of a row that starts with a gap, and the prefetch).
pub(crate) fn row_range_exception(error: RowRangeError) -> NcbiException {
    let (line, kind, message) = match error {
        RowRangeError::Empty => (
            227,
            "CSeqalignException::eInvalidAlignment",
            "CSeq_align::GetSeqRange(): Row is empty",
        ),
        RowRangeError::SeqidsNotConsistent => (
            186,
            "CSeqalignException::eInvalidRowNumber",
            "CSeq_align::GetSeqRange(): Row seqids not consistent. Cannot determine range.",
        ),
    };
    NcbiException::new(
        ThrowSite {
            file: SEQ_ALIGN_CPP,
            line,
            kind,
            function: None,
        },
        message,
    )
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:788-803
/// `CreateDensegFromStdseg(): Std-seg segment lengths not accurate!` (Seq_align.cpp:801): the
/// Std-seg of `-subject` that the Seq-loc mapper rewrote, a query row of one
/// nucleotide and an empty subject row (0 is neither 3 nor a third of 1, nor 1), where a formatter
/// makes a Dense-seg of it (the tabular rows' alignment block, BLAST XML).
pub(crate) fn std_seg_lengths_not_accurate() -> NcbiException {
    NcbiException::new(
        ThrowSite {
            file: SEQ_ALIGN_CPP,
            line: 801,
            kind: "CSeqalignException::eInvalidInputAlignment",
            function: None,
        },
        "CreateDensegFromStdseg(): Std-seg segment lengths not accurate!",
    )
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:171-191
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:226-230
/// What `CSeq_align::GetSeqRange` of the row finds on the Std-segs of a translated search's HSP from the
/// query offset `offset` to `end` (exclusive, in its frame's letters) whose edit script is `ops`
/// (empty for an ungapped HSP): the exception it throws (`RowRangeError`), `None` where it gives
/// the range. An `Ins` keeps the query's letters against a gap of the subject's row, a `Del` the
/// subject's against a gap of the query's (`collect_seq_align_data`). The one predicate of the
/// writers that read a row's range (BLASTX's and RPSTBLASTN's reports through
/// `algorithm::blastx::empty_rows::row_range`, the pairwise report's `PairwiseHit::row_range`).
pub(crate) fn std_seg_row_range(
    offset: i64,
    end: i64,
    ops: &[GapEditOp],
    row: AlignRow,
) -> Option<RowRangeError> {
    if ungapped_rows_empty(offset, end) {
        return Some(RowRangeError::Empty);
    }
    gapped_row_range(ops, row)
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:171-191
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:226-230
/// What `CSeq_align::GetSeqRange` of `row` finds on a reported (already `-subject`-mapped,
/// `mapped_script`) Std-seg from `offset` to `end` with edit script `ops`: an empty row
/// (`RowRangeError::Empty`) unless `empty_query_rewritten` holds for `AlignRow::Query`
/// (`-subject`'s Seq-loc mapper rewrote the empty query row into a non-empty nucleotide,
/// `mapper_rewrites_empty_rows`), and any leading gap left in the reported script
/// (`RowRangeError::SeqidsNotConsistent`, for `-db` and `-subject` alike). The one owner of
/// what the formatters read on a reported Std-seg row (`report::pairwise::seen_row_range`,
/// `blastx::empty_rows::{reported_row_range, reported_rows_range}`, `tblastn::blast_xml`,
/// `tblastn::tabular`).
#[inline]
pub(crate) fn reported_std_seg_row_range(
    offset: i64,
    end: i64,
    ops: &[GapEditOp],
    row: AlignRow,
    empty_query_rewritten: bool,
) -> Option<RowRangeError> {
    let error = std_seg_row_range(offset, end, ops, row)?;
    match (error, row) {
        (RowRangeError::Empty, AlignRow::Query) if empty_query_rewritten => None,
        _ => Some(error),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:2907-2930
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:243-245
/// What a formatter that reads the query row's range of a reported Std-seg and then the
/// subject's finds (`reported_std_seg_row_range`), query first: the description table of
/// outfmt 0 to 4 (`s_ProcessAlignSet`) and BLAST XML (`s_SeqAlignToXMLHsp`).
#[inline]
pub(crate) fn reported_std_seg_rows_range(
    offset: i64,
    end: i64,
    ops: &[GapEditOp],
    empty_query_rewritten: bool,
) -> Option<RowRangeError> {
    reported_std_seg_row_range(offset, end, ops, AlignRow::Query, empty_query_rewritten).or_else(
        || reported_std_seg_row_range(offset, end, ops, AlignRow::Subject, empty_query_rewritten),
    )
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:171-191
/// `std_seg_row_range` of a gapped HSP whose edit script is `ops` (never of length 0): the row that
/// the script starts with a gap in (`RowRangeError::SeqidsNotConsistent`, `row_end_gaps`).
#[inline]
pub(crate) fn gapped_row_range(ops: &[GapEditOp], row: AlignRow) -> Option<RowRangeError> {
    (row_end_gaps(ops, row).leading > 0).then_some(RowRangeError::SeqidsNotConsistent)
}

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:969-995
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1100-1122
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2576-2613
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:279-297
/// The end gaps of the row `row` of an HSP whose edit script is `ops`: the gap operations of the
/// row (a `Del` in the query's, an `Ins` in the subject's) before the script's first operation
/// with a letter of the row and after its last, the columns outside the row's aligned range
/// (`CAlnMap::GetSeqAlnRange`). NCBI's readers of an alignment see them so:
/// - the tables' `length`, `gaps` and `gapopen` (`GetAlignLengths`) and outfmt 0's `Gaps`
///   (`x_GetNumGaps`) count the columns inside the query row's range only (`RowEndGaps::columns`
///   and `RowEndGaps::opens` of `AlignRow::Query` come off the alignment's counts);
/// - the rows' strings (`GetWholeAlnSeqString`) show them as the end character, a space (outfmt 0,
///   the tables' `qseq`/`sseq`/`btop`, BLAST XML and the XML2 and JSON documents);
/// - a Std-seg row that starts with one has no Seq-id there (`gapped_row_range`).
///
/// Only a Smith-Waterman alignment of `-use_sw_tback` has them (`ALIGN_EX` from the
/// Smith-Waterman start may start with a gap); an X-drop alignment starts and ends with a
/// substitution. A script without a letter of the row (never made) has none.
pub(crate) fn row_end_gaps(ops: &[GapEditOp], row: AlignRow) -> RowEndGaps {
    let has_letter = |op: &GapEditOp| {
        !matches!(
            (op, row),
            (GapEditOp::Del(_), AlignRow::Query) | (GapEditOp::Ins(_), AlignRow::Subject)
        )
    };
    let (Some(first), Some(last)) = (
        ops.iter().position(has_letter),
        ops.iter().rposition(has_letter),
    ) else {
        return RowEndGaps::default();
    };
    let columns = |gaps: &[GapEditOp]| gaps.iter().map(|op| op.num() as usize).sum();
    RowEndGaps {
        leading: columns(&ops[..first]),
        trailing: columns(&ops[last + 1..]),
        opens: first + (ops.len() - 1 - last),
    }
}

/// The columns, gap letters and gap openings of an alignment (`length`, `gaps` and `gapopen` of
/// NCBI's tables).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct AlignLengths {
    pub(crate) length: usize,
    pub(crate) gaps: usize,
    pub(crate) gap_opens: usize,
}

impl AlignLengths {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1100-1122
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:960-968
    /// The counts of the whole alignment of an HSP whose edit script is `ops` as NCBI's tables
    /// count them (`GetAlignLengths` over the query row's aligned range): the query row's end
    /// gaps (`row_end_gaps`, a Smith-Waterman alignment of `-use_sw_tback`) come off the three.
    /// The one count of every table writer (`report::tabular::TabularHsp::within_query_range`, the
    /// default fields of BLASTP and TBLASTN), of outfmt 0's `Gaps` and of the library's
    /// `HspRecord` (`report::pairwise::PairwiseHit::query_range_lengths`); an X-drop alignment
    /// has no end gaps.
    pub(crate) fn within_query_range(self, ops: &[GapEditOp]) -> Self {
        let end_gaps = row_end_gaps(ops, AlignRow::Query);
        Self {
            length: self.length.saturating_sub(end_gaps.columns()),
            gaps: self.gaps.saturating_sub(end_gaps.columns()),
            gap_opens: self.gap_opens.saturating_sub(end_gaps.opens),
        }
    }
}

/// The end gaps of a row of an alignment (`row_end_gaps`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct RowEndGaps {
    /// The columns before the row's first letter.
    pub(crate) leading: usize,
    /// The columns after the row's last letter.
    pub(crate) trailing: usize,
    /// The gap operations among them (each one gap opening).
    pub(crate) opens: usize,
}

impl RowEndGaps {
    /// The columns of the end gaps.
    pub(crate) fn columns(self) -> usize {
        self.leading + self.trailing
    }

    /// The columns of the row of `len` columns between its end gaps (`lo..hi`).
    pub(crate) fn span(self, len: usize) -> (usize, usize) {
        let lo = self.leading.min(len);
        (lo, len.saturating_sub(self.trailing).max(lo))
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:453-460
/// Whether `-subject`'s Seq-loc mapper (`remaps`, `RunRecord::remaps_subjects`) rewrites the
/// empty rows of a Seq-align whose Std-segs' rows are empty or not (`rows_empty`, one per Std-seg)
/// into the query nucleotide of `write_empty_rows` (`mapped_empty_query_nucleotide`): it learns
/// the sequence types from a Std-seg with rows of letters only (`x_Init`'s first pass, above;
/// an empty row's length 0 is skipped). A Seq-align of empty rows only gets no type and keeps its
/// rows where their stop does not wrap (`algorithm::blastx::empty_rows::SubjectLoc::of_list`;
/// where it does, BLASTX refused the report before). The one rule of the pairwise report
/// (`pairwise::mapper_rewrites`) and the Seq-aligns (`write_seq_align_diags`).
pub(crate) fn mapper_rewrites_empty_rows(
    remaps: bool,
    mut rows_empty: impl Iterator<Item = bool>,
) -> bool {
    remaps && rows_empty.any(|empty| !empty)
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1538-1559
/// The query row of the Std-seg that `-subject`'s Seq-loc mapper makes of an HSP of length 0
/// (`ungapped_rows_empty`) at a subject offset above 0, whose raw query row starts at `from`
/// (`ungapped_interval`, in the searched letters): the one nucleotide `from` on a plus frame,
/// `from - 1` on a minus frame (the mapper takes the frames from the empty protein row and pads
/// the nucleotide row to codons, `x_GetDstStd`); the subject row stays empty. Probes: c0054 and
/// P1's p0017, p0077, p0095, p0106 (`X/SURVEY_X2g_formats.md` 1). The one owner of that row
/// (`write_empty_rows`, BLASTX's `qcovs`).
pub(crate) fn mapped_empty_query_nucleotide(from: i32, frame: i8) -> i32 {
    if frame > 0 {
        from
    } else {
        from - 1
    }
}

/// The `diags` Seq-align of one subject's HSPs in an ungapped translated search: no `dim`, no
/// `score`, one Std-seg per HSP with the HSP's scores (the coverage over the query length
/// of the construction, a translated query's in nucleotides), and the remaps of
/// `write_std_segs` (an ungapped Std-seg has no gap row, so the mapper keeps them all).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1323-1346,1379-1386
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1406-1447
/// (`write_use_this_seqid` is `s_AddUserObjectToSeqAlign`.)
fn write_seq_align_diags<E: Encoder>(
    enc: &mut E,
    query: &QueryRecord<'_>,
    subject: &SubjectRecord<'_>,
    hsps: &[HspRecord<'_>],
    remaps: bool,
    after_each: &mut impl FnMut(&mut E),
) {
    enc.element();
    enc.begin_sequence("Seq-align");
    enc.member("type", 0);
    enc.enumerated("diags", 2);
    enc.end_member();
    enc.member("segs", 3);
    enc.variant("Seq-align_segs", "std", 2);
    enc.begin_list(false);
    let remaps = mapper_rewrites_empty_rows(remaps, hsps.iter().map(HspRecord::rows_empty));
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
        if hsp.rows_empty() {
            write_empty_rows(
                enc,
                query,
                subject,
                hsp,
                remaps,
                (q_from + query.offset, s_from + subject.offset),
            );
        } else {
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
        }
        end_std_seg_locs(enc);
        enc.member("scores", 3);
        enc.begin_list(true);
        write_score_list(enc, hsp, query.searched_length);
        enc.end_list();
        enc.end_member();
        enc.end_sequence();
        // A subject's HSPs are one Seq-align: hand the final bytes over after each Std-seg.
        after_each(enc);
    }
    enc.end_list();
    enc.end_variant();
    enc.end_member();
    write_use_this_seqid(enc, subject.use_this_seqid);
    enc.end_sequence();
}

/// The `diags` Seq-align of one subject's HSPs in an ungapped search whose rows are not
/// translated (BLASTN `-ungapped`, BLASTP `-ungapped`): no `dim`, no `score`, one Dense-diag per
/// HSP with the HSP's scores (the coverage over the query length of the construction), and the
/// starts that `RemapToQueryLoc` and `s_RemapToSubjectLoc` leave (see below).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1273-1311
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1406-1431
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:51-59
/// A BLASTN HSP's subject row is on the plus strand (frame +1), so its start is the offset;
/// the query row of a minus-strand HSP starts at `query_length - end` (the searched length:
/// `query.searched_length`, with `-query_loc` the interval's). A BLASTP HSP has frame 0 on both
/// rows, so both strands are `unknown` and both starts are the offsets.
///
/// The starts after `RemapToQueryLoc` (`CDense_diag::OffsetRow`: the `-query_loc` start added to
/// the query row) and, where the subject source can return partial sequences,
/// `s_RemapToSubjectLoc` (the `-subject_loc` start added to the subject row, `remaps`).
///
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_diag.cpp:72-96
fn write_seq_align_dendiag<E: Encoder>(
    enc: &mut E,
    query: &QueryRecord<'_>,
    subject: &SubjectRecord<'_>,
    hsps: &[HspRecord<'_>],
    after_each: &mut impl FnMut(&mut E),
) {
    enc.element();
    enc.begin_sequence("Seq-align");
    enc.member("type", 0);
    enc.enumerated("diags", 2);
    enc.end_member();
    enc.member("segs", 3);
    enc.variant("Seq-align_segs", "dendiag", 0);
    enc.begin_list(false);
    for hsp in hsps {
        let DendiagSegment {
            starts: [q_start, s_start],
            len,
            strands,
        } = dendiag_segment(query, subject, hsp);
        enc.element();
        enc.begin_sequence("Dense-diag");
        enc.member("dim", 0);
        enc.int(BLAST_ALIGNMENT_DIM);
        enc.end_member();
        enc.member("ids", 1);
        enc.begin_list(false);
        enc.element();
        query.id.write(enc);
        enc.element();
        subject.id.write(enc);
        enc.end_list();
        enc.end_member();
        enc.member("starts", 2);
        enc.begin_list(false);
        enc.element();
        enc.int(q_start + query.offset);
        enc.element();
        enc.int(s_start + subject.offset);
        enc.end_list();
        enc.end_member();
        enc.member("len", 3);
        enc.int(len);
        enc.end_member();
        enc.member("strands", 4);
        enc.begin_list(false);
        for strand in strands {
            enc.element();
            enc.begin_named("Na-strand");
            enc.enumerated(strand.name(), strand.value());
            enc.end_named();
        }
        enc.end_list();
        enc.end_member();
        enc.member("scores", 5);
        enc.begin_list(true);
        write_score_list(enc, hsp, query.searched_length);
        enc.end_list();
        enc.end_member();
        enc.end_sequence();
        // A subject's HSPs are one Seq-align: hand the final bytes over after each Dense-diag.
        after_each(enc);
    }
    enc.end_list();
    enc.end_variant();
    enc.end_member();
    write_use_this_seqid(enc, subject.use_this_seqid);
    enc.end_sequence();
}

/// The Dense-diag of an ungapped HSP of an untranslated search, without its ids and scores.
struct DendiagSegment {
    /// Query, subject: in the searched letters, minus the range starts that the caller adds
    /// (`QueryRecord::offset`, `SubjectRecord::offset`).
    starts: [i32; 2],
    len: i32,
    /// Query, subject.
    strands: [Strand; 2],
}

/// The Dense-diag of an ungapped HSP (`x_UngappedHSPToDenseDiag`): a row of a plus frame (or of
/// no frame) starts at its offset, a row of a minus frame at `searched_length - end`; the length
/// is the query row's; each row's strand is its frame's. The one source of the Dense-diags that
/// outfmt 8, 9, 11 and 12 write (`write_seq_align_dendiag`) and of the segments that SAM joins
/// (`ungapped_segments`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1289-1306
fn dendiag_segment(
    query: &QueryRecord<'_>,
    subject: &SubjectRecord<'_>,
    hsp: &HspRecord<'_>,
) -> DendiagSegment {
    let q_start = if hsp.q_frame >= 0 {
        hsp.q_offset
    } else {
        query.searched_length - hsp.q_end
    };
    let s_start = if hsp.s_frame >= 0 {
        hsp.s_offset
    } else {
        subject.searched_length - hsp.s_end
    };
    DendiagSegment {
        starts: [q_start, s_start],
        len: hsp.q_end - hsp.q_offset,
        strands: [Strand::of_frame(hsp.q_frame), Strand::of_frame(hsp.s_frame)],
    }
}

/// The Dense-seg that NCBI's SAM formatter makes of the Seq-align of an ungapped search's
/// subject (`CreateDensegFromDendiag`): one segment per HSP in the Seq-align's order, each the
/// HSP's Dense-diag (`dendiag_segment`), minus the range starts that the caller adds. The SAM
/// formatter reads them as it reads a gapped HSP's (`gapped_segments`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/sam.cpp:67-82
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1149-1208
pub(crate) fn ungapped_segments(
    query: &QueryRecord<'_>,
    subject: &SubjectRecord<'_>,
    hsps: &[HspRecord<'_>],
    segments: &mut Segments,
) {
    segments.starts.clear();
    segments.lens.clear();
    segments.strands.clear();
    for hsp in hsps {
        let segment = dendiag_segment(query, subject, hsp);
        segments.starts.extend(segment.starts);
        segments.lens.push(segment.len);
        segments.strands.extend(segment.strands);
    }
}

/// The two Seq-intervals of the Std-seg of an ungapped HSP of length 0 (BLASTX's one-hit word
/// finder with `-window_size 0`, `algorithm::blastx::empty_rows`), whose rows start at `from`
/// (query, subject) and end one before (`to = from - 1`, in NCBI's unsigned `TSeqPos`: a row
/// starting at 0 ends at 4294967295). With `-subject` (`remaps`, a Seq-align with rows of letters
/// too) NCBI's Seq-loc mapper rewrites the query row to one nucleotide with a fuzz
/// (`mapped_empty_query_nucleotide`); the mapper's form of a subject row starting at 0 is refused
/// before (`algorithm::blastx::empty_rows`).
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1538-1559
/// Probes (`X/SURVEY_X2g_formats.md` 1): a plus frame gives the query row `(Q, Q)` with
/// `fuzz-to alt {Q-2}`, a minus frame `(Q-1, Q-1)` with `fuzz-from alt {Q}`
/// (`mapped_empty_query_nucleotide`); the subject row stays.
#[cold]
#[inline(never)]
fn write_empty_rows<E: Encoder>(
    enc: &mut E,
    query: &QueryRecord<'_>,
    subject: &SubjectRecord<'_>,
    hsp: &HspRecord<'_>,
    remaps: bool,
    (q_from, s_from): (i32, i32),
) {
    let q_strand = Strand::of_frame(hsp.q_frame);
    let q_from = q_from as u32;
    enc.element();
    if remaps {
        let at = mapped_empty_query_nucleotide(q_from as i32, hsp.q_frame) as u32;
        let fuzz = if hsp.q_frame > 0 {
            (None, Some(q_from.wrapping_sub(2)))
        } else {
            (Some(q_from), None)
        };
        write_seq_interval_pos(enc, &query.id, (at, at), q_strand, fuzz);
    } else {
        write_seq_interval_pos(
            enc,
            &query.id,
            (q_from, q_from.wrapping_sub(1)),
            q_strand,
            (None, None),
        );
    }
    let s_from = s_from as u32;
    enc.element();
    write_seq_interval_pos(
        enc,
        &subject.id,
        (s_from, s_from.wrapping_sub(1)),
        Strand::of_frame(hsp.s_frame),
        (None, None),
    );
}

/// A Seq-interval whose `from` and `to` are NCBI's unsigned `TSeqPos`, with the `alt` fuzzes
/// `fuzz` (`fuzz-from`, `fuzz-to`) of the Seq-loc mapper (`write_empty_rows`).
///
/// NCBI reference (598d8ae6): c++/src/objects/general/general.asn:106-120
/// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:93-99
fn write_seq_interval_pos<E: Encoder>(
    enc: &mut E,
    id: &SeqId<'_>,
    (from, to): (u32, u32),
    strand: Strand,
    fuzz: (Option<u32>, Option<u32>),
) {
    enc.variant("Seq-loc", "int", 3);
    enc.begin_sequence("Seq-interval");
    enc.member("from", 0);
    enc.big_int(i64::from(from));
    enc.end_member();
    enc.member("to", 1);
    enc.big_int(i64::from(to));
    enc.end_member();
    enc.member("strand", 2);
    enc.begin_named("Na-strand");
    enc.enumerated(strand.name(), strand.value());
    enc.end_named();
    enc.end_member();
    enc.member("id", 3);
    id.write(enc);
    enc.end_member();
    for (name, tag, alt) in [("fuzz-from", 4, fuzz.0), ("fuzz-to", 5, fuzz.1)] {
        let Some(alt) = alt else { continue };
        enc.member(name, tag);
        enc.variant("Int-fuzz", "alt", 4);
        enc.begin_list(true);
        enc.element();
        enc.big_int(i64::from(alt));
        enc.end_list();
        enc.end_variant();
        enc.end_member();
    }
    enc.end_sequence();
    enc.end_variant();
}

/// The `ext` of a Seq-align under a positive list (`seqids`, `SubjectRecord::use_this_seqid`):
/// nothing when `seqids` is empty.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1233-1245
/// (`num` is the number of the strings, `CUser_field::SetValue`; `data`: `strs` 6,
/// general.asn:136-147.)
fn write_use_this_seqid<E: Encoder>(enc: &mut E, seqids: &[Vec<u8>]) {
    if seqids.is_empty() {
        return;
    }
    enc.member("ext", 6);
    enc.begin_list(false);
    enc.element();
    enc.begin_sequence("User-object");
    enc.member("type", 1);
    enc.variant("Object-id", "str", 1);
    enc.visible_str(b"use_this_seqid");
    enc.end_variant();
    enc.end_member();
    enc.member("data", 2);
    enc.begin_list(false);
    enc.element();
    enc.begin_sequence("User-field");
    enc.member("label", 0);
    enc.variant("Object-id", "str", 1);
    enc.visible_str(b"SEQIDS");
    enc.end_variant();
    enc.end_member();
    enc.member("num", 1);
    enc.int(seqids.len() as i32);
    enc.end_member();
    enc.member("data", 2);
    enc.variant_string_list("User-field_data", "strs", 6, seqids);
    enc.end_member();
    enc.end_sequence();
    enc.end_list();
    enc.end_member();
    enc.end_sequence();
    enc.end_list();
    enc.end_member();
}

/// The Seq-align of one HSP of a gapped search, or nothing for an HSP of score 0.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:627-656
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:659-690
/// (No LOSAT edit script has `eGapAlignDecline`, so the `disc` branch does not occur.)
/// NCBI reference (598d8ae6): c++/src/objects/seqalign/seqalign.asn:22-40
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
    enc.member("score", 2);
    enc.begin_list(true);
    // psiblast's Seq-aligns leave out a `num_ident` of 0 (`seq_align_num_ident`).
    let hsp = HspRecord {
        num_ident: seq_align_num_ident(i64::from(hsp.num_ident), run.program_type).unwrap_or(-1),
        ..*hsp
    };
    write_score_list(enc, &hsp, coverage_query_length(run, query.searched_length));
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
    write_use_this_seqid(enc, subject.use_this_seqid);
    enc.end_sequence();
}

/// Whether neither row of the program's searches is translated (BLASTN, BLASTP): the programs
/// whose ungapped Seq-aligns are Dense-diags (`Blast_QueryIsTranslated`, `Blast_SubjectIsTranslated`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1417-1431
/// (Every other program, `else // Translated search`, makes Std-segs. RPS-BLAST, which LOSAT's
/// reports of this writer do not include, makes Dense-diags too.)
fn untranslated(program: Program) -> bool {
    !query_is_translated(program) && !subject_is_translated(program)
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

/// A query's Seq-aligns as runs of its HSPs (`hsps` grouped by subject in the Seq-align order):
/// in a gapped search each HSP of score other than 0 alone (`s_BlastHSP2SeqAlign` makes no
/// Seq-align of an HSP of score 0), in an ungapped search a subject's HSPs together (one `diags`
/// Seq-align per HSP list). The one grouping of the Seq-aligns that outfmt 8, 9, 11 and 12 write
/// (`write_seq_aligns`) and that SAM formats (`report::sam`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:672-674
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1397-1446
pub(crate) fn seq_align_hsps<'h, 'a>(
    run: &RunRecord<'_>,
    hsps: &'h [HspRecord<'a>],
) -> impl Iterator<Item = &'h [HspRecord<'a>]> {
    let gapped = run.gapped;
    hsps.chunk_by(move |a, b| !gapped && a.subject == b.subject)
        .filter(move |align| !gapped || align[0].score != 0)
}

/// The Seq-aligns of one query's HSPs, as the elements of the list that the caller opened
/// (`Seq-annot` `data align` of outfmt 8, `alignments` of the archive of outfmt 11).
/// `after_each` is called after each Seq-align, and after each Std-seg of an ungapped search's
/// Seq-align (it holds all of a subject's HSPs); the archive hands the encoder's final bytes to
/// its sink there.
pub(crate) fn write_seq_aligns<E: Encoder>(
    enc: &mut E,
    run: &RunRecord<'_>,
    query: &QueryRecord<'_>,
    subjects: &[SubjectRecord<'_>],
    hsps: &[HspRecord<'_>],
    segments: &mut Segments,
    mut after_each: impl FnMut(&mut E),
) {
    if run.gapped {
        for align in seq_align_hsps(run, hsps) {
            let hsp = &align[0];
            let subject = &subjects[hsp.subject as usize];
            write_seq_align_gapped(enc, run, query, subject, hsp, segments);
            after_each(enc);
        }
    } else if untranslated(run.program) {
        for list in seq_align_hsps(run, hsps) {
            write_seq_align_dendiag(
                enc,
                query,
                &subjects[list[0].subject as usize],
                list,
                &mut after_each,
            );
            after_each(enc);
        }
    } else {
        for list in seq_align_hsps(run, hsps) {
            write_seq_align_diags(
                enc,
                query,
                &subjects[list[0].subject as usize],
                list,
                run.remaps_subjects(),
                &mut after_each,
            );
            after_each(enc);
        }
    }
}

/// The Seq-annot of one query's Seq-aligns: the `desc` user objects and `data align` with,
/// subject by subject in the order of `hsps`, one Seq-align per HSP of score other than 0
/// (a gapped search) or one `diags` Seq-align per subject (an ungapped search).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1606-1626
///
/// NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:460-476
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:381-386
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
    write_seq_aligns(enc, run, query, subjects, hsps, &mut segments, |_| {});
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
}

impl std::fmt::Display for SeqAlignError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SeqAlignError::Serial(error) => error.fmt(f),
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

/// Writes the document of one query to `writer` through `buf`, which is cleared first: what
/// `write_query` makes, or nothing when the query has no Seq-align. A document with a string byte
/// that NCBI's ASN.1 text or binary stream writes as `#` with an error message (the default
/// `eFNP_ReplaceAndWarn`: a database title, a `-parse_deflines` local id with DEL, ...) is
/// written again through `Traced` (`write_query_posting`), with NCBI's posts of those bytes on
/// `diagnostics`; a document without one is written as before.
///
/// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:268-279
///
/// Each document is one object of NCBI's object stream on the report's stream
/// (`ostream::Site::Serial`): its bytes, then `~CObjectOStream`'s flush. The plain document is
/// inserted whole; where its device write fails, NCBI's pieces of it (`document_pieces`) decide
/// NCBI's end.
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_query_to(
    format: SerialFormat,
    run: &RunRecord<'_>,
    query: &QueryRecord<'_>,
    subjects: &[SubjectRecord<'_>],
    hsps: &[HspRecord<'_>],
    buf: &mut Vec<u8>,
    writer: &mut dyn ReportWrite,
    diagnostics: &mut dyn Write,
) -> anyhow::Result<()> {
    buf.clear();
    let document = QueryDocument {
        run,
        query,
        subjects,
        hsps,
    };
    match write_query(format, run, query, subjects, hsps, buf) {
        Ok(true) => {
            writer.site(Site::Serial { whole: true })?;
            match writer.write_all(buf).and_then(|()| writer.ncbi_flush()) {
                Err(error) if crate::report::ostream::is_undecided_report_write_failure(&error) => {
                    writer.resolve_object(&document_pieces(format, &document))?;
                }
                written => written?,
            }
            writer.site(Site::Report)?;
        }
        Ok(false) => {}
        Err(SeqAlignError::Serial(SerialError::NonVisibleChar(byte))) => {
            writer.site(Site::Serial { whole: false })?;
            write_query_posting(format, byte, &document, buf, writer, diagnostics)?;
            writer.ncbi_flush()?;
            writer.site(Site::Report)?;
        }
        Err(error) => return Err(error.into()),
    }
    Ok(())
}

/// `document_pieces` of one query's document (for the tests).
#[cfg(test)]
pub(super) fn query_pieces(
    format: SerialFormat,
    run: &RunRecord<'_>,
    query: &QueryRecord<'_>,
    subjects: &[SubjectRecord<'_>],
    hsps: &[HspRecord<'_>],
) -> Vec<usize> {
    document_pieces(
        format,
        &QueryDocument {
            run,
            query,
            subjects,
            hsps,
        },
    )
}

/// NCBI's insertions of `document` into the report's stream (`Traced::pieces`): the writes of
/// its object stream's buffer, from a second encoding of it; none for a stream without a model of
/// its buffer.
#[cold]
#[inline(never)]
fn document_pieces(format: SerialFormat, document: &QueryDocument<'_, '_>) -> Vec<usize> {
    fn pieces<E: Encoder>(mut enc: Traced<E>, document: &QueryDocument<'_, '_>) -> Vec<usize> {
        write_seq_annot(
            &mut enc,
            document.run,
            document.query,
            document.subjects,
            document.hsps,
        );
        enc.piece_sizes()
    }
    let mut bytes = Vec::new();
    match format {
        SerialFormat::AsnText => pieces(Traced::new(AsnTextEncoder::new(&mut bytes), 0), document),
        SerialFormat::AsnBinary => {
            pieces(Traced::new(AsnBinaryEncoder::new(&mut bytes), 0), document)
        }
        SerialFormat::Json => pieces(Traced::new(JsonEncoder::new(&mut bytes), 0), document),
    }
}

/// The records of one query's Seq-annot.
struct QueryDocument<'r, 'a> {
    run: &'r RunRecord<'a>,
    query: &'r QueryRecord<'a>,
    subjects: &'r [SubjectRecord<'a>],
    hsps: &'r [HspRecord<'a>],
}

/// `write_query_to` of a document that the plain encoder of `format` stopped on (the byte `byte`,
/// `NonVisibleChar`): written again through `Traced` into `buf`, and handed to `writer` with
/// NCBI's post of each replaced byte (`ReplacedChar::post` with the application's prefix,
/// `options::RequestKind::application`, the one table of the applications' names) on
/// `diagnostics` where NCBI posts it: after the bytes of the document
/// that NCBI's buffer (`COStreamBuffer`, `out_buffer`) had written to its output by then, the
/// output flushed (`cerr` is tied to `cout`; `Traced::flush_to`). NCBI writes each query's
/// Seq-annot with an object stream of its own (`m_Outfile` takes `MSerial_AsnText` and then
/// `annot`, cited at the module), so the positions (`line N`, `byte N`) and the buffer start again
/// with each query's document. The JSON stream escapes such a byte (objostrjson.cpp:731-760), so
/// its encoder never stops on one.
///
/// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:413-444
/// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:843-861
#[cold]
#[inline(never)]
fn write_query_posting(
    format: SerialFormat,
    byte: u8,
    document: &QueryDocument<'_, '_>,
    buf: &mut Vec<u8>,
    writer: &mut dyn ReportWrite,
    diagnostics: &mut dyn Write,
) -> anyhow::Result<()> {
    buf.clear();
    let prefix =
        crate::report::archive::options::RequestKind::from(document.run.program).application();
    let post = |replaced: &ReplacedChar| {
        diagnostics.write_all(&replaced.post(Some(prefix)))?;
        diagnostics.flush()
    };
    match format {
        SerialFormat::AsnText => write_traced(
            Traced::new(AsnTextEncoder::new(buf), 0).posting(post),
            document,
            writer,
        ),
        SerialFormat::AsnBinary => write_traced(
            Traced::new(AsnBinaryEncoder::new(buf), 0).posting(post),
            document,
            writer,
        ),
        SerialFormat::Json => Err(SeqAlignError::Serial(SerialError::NonVisibleChar(byte)).into()),
    }
}

/// The document written through `enc` (`Traced` with its post function) and handed to `writer`
/// with the posts; nothing on an error of the encoder.
fn write_traced<E: Encoder>(
    mut enc: E,
    document: &QueryDocument<'_, '_>,
    writer: &mut dyn ReportWrite,
) -> anyhow::Result<()> {
    write_seq_annot(
        &mut enc,
        document.run,
        document.query,
        document.subjects,
        document.hsps,
    );
    if let Some(error) = enc.error() {
        return Err(SeqAlignError::Serial(error).into());
    }
    enc.flush_to(writer, true)?;
    Ok(())
}

/// Writes the documents of every query of `hit_lists` (indexed by query) to `writer`, each
/// query's warnings before it and NCBI's posts of the bytes that its stream replaces in the query's
/// document after them (`write_query_to`, on the warnings' sink). `records` fills the subject and
/// HSP records of the hit list of query `q_idx` (replacing their contents) for `write_query_to`; a
/// query without a hit list writes nothing but its warnings.
pub(crate) fn write_hit_list_seqaligns<'a, L>(
    format: SerialFormat,
    run: &RunRecord<'_>,
    hit_lists: &'a [Option<L>],
    queries: &QueryNames<'_>,
    writer: &mut dyn ReportWrite,
    warnings: &mut QueryWarnings<'_>,
    mut records: impl FnMut(usize, &'a L, &mut Vec<SubjectRecord<'a>>, &mut Vec<HspRecord<'a>>),
) -> anyhow::Result<()> {
    let mut subject_records: Vec<SubjectRecord<'a>> = Vec::new();
    let mut hsp_records: Vec<HspRecord<'a>> = Vec::new();
    let mut buf: Vec<u8> = Vec::new();
    for (q_idx, hit_list) in hit_lists.iter().enumerate() {
        // The query's warnings come before its document (`QueryWarnings`).
        warnings.before_query(q_idx, &mut *writer)?;
        let Some(hit_list) = hit_list else {
            continue;
        };
        records(q_idx, hit_list, &mut subject_records, &mut hsp_records);
        write_query_to(
            format,
            run,
            &queries.record(q_idx),
            &subject_records,
            &hsp_records,
            &mut buf,
            &mut *writer,
            &mut *warnings.sink,
        )?;
    }
    Ok(())
}
