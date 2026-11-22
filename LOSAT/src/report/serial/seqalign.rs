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
//! (`remap_strands` for a Dense-seg, `std_seg_kept` for Std-segs).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:676-707
//! `PruneSeqalign` keeps the Seq-aligns of the first `m_HitlistSize` subjects, which is the
//! search's own hit list size for these formats, so it keeps every Seq-align of the search's
//! hit lists (outfmt-20261019 S/DESIGN_OF3.md §2.1).

use std::io::Write;

use crate::blastinput::app::Program;
use crate::common::GapEditOp;
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
    crate::report::query_anchored::HspStrands {
        query_minus: q_frame < 0,
        subject_minus: s_frame < 0,
        query_frame: q_frame,
        no_widths: false,
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

/// The "Blast Type" of the program: `EProgramToTaskName(ProgramNameToEnum(program name))`
/// and the `EProgram` value. The program name of every BLASTN task is "blastn"
/// (`Blast_ProgramNameFromType`, blast_format.cpp:94), and of every BLASTP task "blastp".
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_types.hpp:56-62
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:676-686
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:330-381
/// A `-subject` search has no database name (`n/a`). A `-db` named `n/a` or `SRA` is
/// refused by the search programs (`check_database`), so their title is written (a byte that
/// NCBI's ASN.1 text or binary stream replaces is written as `#` with NCBI's post,
/// `write_query_to`); blast_formatter's database name is the archive's (`DatabaseDesc::of_name`).
fn write_desc<E: Encoder>(enc: &mut E, desc: &AnnotDesc<'_>) {
    enc.member("desc", 3);
    enc.begin_list(true);
    enc.element();
    write_user_object(
        enc,
        b"Hist Seqalign",
        b"Hist Seqalign",
        UserData::Bool(true),
    );
    enc.element();
    write_user_object(
        enc,
        b"Blast Type",
        desc.task.as_bytes(),
        UserData::Int(desc.program),
    );
    enc.element();
    match desc.database {
        DatabaseDesc::None => {
            write_user_object(enc, b"Blast Database Title", b"n/a", UserData::Bool(false))
        }
        DatabaseDesc::Sra(name) => {
            write_user_object(enc, b"Blast Database Title", name, UserData::Bool(true))
        }
        DatabaseDesc::Title { title, nucleotide } => write_user_object(
            enc,
            b"Blast Database Title",
            title,
            UserData::Bool(nucleotide),
        ),
    }
    enc.end_list();
    enc.end_member();
}

/// What the Seq-annot's `desc` says of a search: `CreateSeqAnnotFromSeqAlignSet`'s `program`
/// and its database (quoted at `write_desc`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct AnnotDesc<'a> {
    /// `EProgramToTaskName(program)` and the number of `program` (`EProgram`): the label and
    /// value of the "Blast Type" user object.
    pub(crate) task: &'static str,
    pub(crate) program: i32,
    /// The "Blast Database Title" user object.
    pub(crate) database: DatabaseDesc<'a>,
}

impl<'a> AnnotDesc<'a> {
    /// The `desc` of a search program's run `run`: its program's type, and its database's title
    /// (`RunRecord::db_title`), or `n/a` for `-subject`.
    pub(crate) fn of_run(run: &RunRecord<'a>) -> Self {
        let (task, program) = blast_type(run.program);
        AnnotDesc {
            task,
            program,
            database: match run.db_title {
                None => DatabaseDesc::None,
                Some(title) => DatabaseDesc::Title {
                    title,
                    nucleotide: subject_is_nucleotide(run.program),
                },
            },
        }
    }
}

/// The value of the "Blast Database Title" user object (`CreateSeqAnnotFromSeqAlignSet`, quoted
/// at `write_desc`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum DatabaseDesc<'a> {
    /// No database name (a search of `-subject` sequences), or the name `n/a`: `n/a`, false.
    None,
    /// The database name `SRA` (any case): the name, true.
    Sra(&'a [u8]),
    /// Another database: its title (`db_title`: the titles of the names of its list, each its
    /// name when it has none, joined by `; `) and whether the program's subjects are nucleotide
    /// (`Blast_SubjectIsNucleotide`).
    Title { title: &'a [u8], nucleotide: bool },
}

impl<'a> DatabaseDesc<'a> {
    /// The value for the database name `name` (`m_DbName`, empty for `-subject` sequences) whose
    /// title is `title`, with nucleotide subjects or not: NCBI's comparisons of the name.
    pub(crate) fn of_name(name: &'a [u8], title: &'a [u8], nucleotide: bool) -> Self {
        if name.is_empty() || name.eq_ignore_ascii_case(b"n/a") {
            DatabaseDesc::None
        } else if name.eq_ignore_ascii_case(b"SRA") {
            DatabaseDesc::Sra(name)
        } else {
            DatabaseDesc::Title { title, nucleotide }
        }
    }
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
fn std_seg_kept(seg: usize, numseg: usize, gap: impl Fn(usize) -> bool) -> bool {
    (seg != 0 || !gap(0)) && (seg + 1 != numseg || !gap(seg))
}

/// The Std-segs of an HSP of a gapped translated search (`segs std`), one per segment, with
/// the remap of the query row (`RemapToQueryLoc`, `CStd_seg::OffsetRow`: the `-query_loc`
/// start added to an interval's ends) and of the subject row (`s_RemapToSubjectLoc` of a
/// `-subject` search, `remapped`: the `-subject_loc` start, on the plus strand that LOSAT's
/// ranges have; the mapper keeps each interval's strand, `CMappingRange::Map_Strand` with
/// the strand set, and `std_seg_kept`).
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
/// start (`CRange(from, from - 1)` is empty, so `CSeq_align::GetSeqRange` throws on them).
/// BLASTX's one-hit word finder with `-window_size 0` reports such HSPs
/// (`algorithm::blastx::empty_rows`); every other reported HSP has rows of its letters. The one
/// predicate of every writer (`HspRecord::rows_empty`, `PairwiseHit::rows_empty`,
/// `empty_rows::rows_empty`).
#[inline]
pub(crate) fn ungapped_rows_empty(offset: i64, end: i64) -> bool {
    end == offset
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
fn write_seq_annot<E: Encoder, D: SeqAnnotDocument>(enc: &mut E, document: &D) {
    enc.begin_document("Seq-annot");
    enc.begin_sequence("Seq-annot");
    write_desc(enc, &document.desc());
    enc.member("data", 4);
    enc.variant("Seq-annot_data", "align", 1);
    enc.begin_list(true);
    document.write_aligns(enc);
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
    let document = QueryDocument {
        run,
        query,
        subjects,
        hsps,
    };
    write_document(format, &document, buf).map_err(SeqAlignError::Serial)?;
    Ok(true)
}

/// Appends to `buf` the Seq-annot of `document` in `format` (the one wrapper: `write_seq_annot`).
/// On an error of NCBI's object stream, `buf` is left as it was.
pub(crate) fn write_document<D: SeqAnnotDocument>(
    format: SerialFormat,
    document: &D,
    buf: &mut Vec<u8>,
) -> Result<(), SerialError> {
    let start = buf.len();
    let error = match format {
        SerialFormat::AsnText => {
            let mut enc = AsnTextEncoder::new(buf);
            write_seq_annot(&mut enc, document);
            enc.error()
        }
        SerialFormat::AsnBinary => {
            let mut enc = AsnBinaryEncoder::new(buf);
            write_seq_annot(&mut enc, document);
            enc.error()
        }
        SerialFormat::Json => {
            let mut enc = JsonEncoder::new(buf);
            write_seq_annot(&mut enc, document);
            enc.error()
        }
    };
    match error {
        Some(error) => {
            buf.truncate(start);
            Err(error)
        }
        None => Ok(()),
    }
}

/// Writes the document of one query to `writer` through `buf`, which is cleared first: what
/// `write_query` makes, or nothing when the query has no Seq-align. A document with a string byte
/// that NCBI's ASN.1 text or binary stream writes as `#` with an error message (the default
/// `eFNP_ReplaceAndWarn`: a database title, a `-parse_deflines` local id with DEL, ...) is
/// written again through `Traced` (`write_document_to`), with NCBI's posts of those bytes (the
/// search program's prefix) on `diagnostics`; a document without one is written as before.
///
/// NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:268-279
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_query_to(
    format: SerialFormat,
    run: &RunRecord<'_>,
    query: &QueryRecord<'_>,
    subjects: &[SubjectRecord<'_>],
    hsps: &[HspRecord<'_>],
    buf: &mut Vec<u8>,
    writer: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> anyhow::Result<()> {
    buf.clear();
    match write_query(format, run, query, subjects, hsps, buf) {
        Ok(true) => writer.write_all(buf)?,
        Ok(false) => {}
        Err(SeqAlignError::Serial(SerialError::NonVisibleChar(byte))) => {
            let document = QueryDocument {
                run,
                query,
                subjects,
                hsps,
            };
            // The prefix of the posts: `options::RequestKind::application`, the one table of
            // the applications' names.
            let prefix =
                crate::report::archive::options::RequestKind::from(run.program).application();
            write_document_posting(
                format,
                byte,
                &document,
                Some(prefix),
                buf,
                writer,
                diagnostics,
            )?;
        }
        Err(error) => return Err(error.into()),
    }
    Ok(())
}

/// Writes the Seq-annot of `document` (`write_document`) to `writer` through `buf`, which is
/// cleared first. A document with a string byte that NCBI's ASN.1 text or binary stream writes
/// as `#` with an error message is written again through `Traced` (`write_document_posting`),
/// with NCBI's posts of those bytes (the diagnostic prefix `prefix`: the search program's,
/// none for blast_formatter) on `diagnostics`; a document without one is written as it is.
pub(crate) fn write_document_to<D: SeqAnnotDocument>(
    format: SerialFormat,
    document: &D,
    prefix: Option<&str>,
    buf: &mut Vec<u8>,
    writer: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> anyhow::Result<()> {
    buf.clear();
    match write_document(format, document, buf) {
        Ok(()) => writer.write_all(buf)?,
        Err(SerialError::NonVisibleChar(byte)) => {
            write_document_posting(format, byte, document, prefix, buf, writer, diagnostics)?;
        }
        Err(error) => return Err(SeqAlignError::Serial(error).into()),
    }
    Ok(())
}

/// What one Seq-annot document holds (`CreateSeqAnnotFromSeqAlignSet`'s arguments): its `desc`
/// and its Seq-aligns. The one wrapper (`write_seq_annot`) writes the search programs' records
/// (`QueryDocument`) and blast_formatter's Seq-aligns as the archive holds them
/// (`blast_formatter::archive`).
pub(crate) trait SeqAnnotDocument {
    /// The Seq-annot's `desc`.
    fn desc(&self) -> AnnotDesc<'_>;
    /// The Seq-aligns, as the elements of the list `data align` that the wrapper opened.
    fn write_aligns<E: Encoder>(&self, enc: &mut E);
}

/// The records of one query's Seq-annot.
struct QueryDocument<'r, 'a> {
    run: &'r RunRecord<'a>,
    query: &'r QueryRecord<'a>,
    subjects: &'r [SubjectRecord<'a>],
    hsps: &'r [HspRecord<'a>],
}

impl SeqAnnotDocument for QueryDocument<'_, '_> {
    fn desc(&self) -> AnnotDesc<'_> {
        AnnotDesc::of_run(self.run)
    }

    /// `write_seq_aligns`: subject by subject in the order of `hsps`, one Seq-align per HSP of
    /// score other than 0 (a gapped search) or one `diags` Seq-align per subject (an ungapped
    /// search).
    fn write_aligns<E: Encoder>(&self, enc: &mut E) {
        let mut segments = Segments::default();
        write_seq_aligns(
            enc,
            self.run,
            self.query,
            self.subjects,
            self.hsps,
            &mut segments,
            |_| {},
        );
    }
}

/// `write_document_to` of a document that the plain encoder of `format` stopped on (the byte
/// `byte`, `NonVisibleChar`): written again through `Traced` into `buf`, and handed to `writer`
/// with NCBI's post of each replaced byte (`ReplacedChar::post` with the application's prefix
/// `prefix`) on `diagnostics` where NCBI posts it: after the bytes of the document
/// that NCBI's buffer (`COStreamBuffer`, `out_buffer`) had written to its output by then, the
/// output flushed (`cerr` is tied to `cout`; `Traced::flush_to`). NCBI writes each query's
/// Seq-annot with an object stream of its own (`m_Outfile << MSerial_AsnText << *annot`, cited at
/// the module), so the positions (`line N`, `byte N`) and the buffer start again with each query's
/// document. The JSON stream escapes such a byte (objostrjson.cpp:731-760), so its encoder never
/// stops on one.
///
/// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:413-444
/// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:843-861
#[cold]
#[inline(never)]
fn write_document_posting<D: SeqAnnotDocument>(
    format: SerialFormat,
    byte: u8,
    document: &D,
    prefix: Option<&str>,
    buf: &mut Vec<u8>,
    writer: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> anyhow::Result<()> {
    buf.clear();
    let post = |replaced: &ReplacedChar| {
        diagnostics.write_all(&replaced.post(prefix))?;
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
fn write_traced<E: Encoder, D: SeqAnnotDocument>(
    mut enc: E,
    document: &D,
    writer: &mut dyn Write,
) -> anyhow::Result<()> {
    write_seq_annot(&mut enc, document);
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
    writer: &mut dyn Write,
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
