//! The reports of BLASTX's ungapped HSPs of length 0: the one-hit word finder of `-ungapped
//! -window_size 0` with a low `-threshold` saves HSPs whose end lies at or before their start
//! (aa_ungapped.c:1054,1083); `Blast_HSPTest` deletes those of negative length and keeps those
//! of length 0 (`runtime.rs`, `core::blast_hits::hsp_test`), which are reported. Their Std-seg
//! rows end one before they start (`to = from - 1`, `x_UngappedHSPToStdSeg`): NCBI's formatters
//! find them empty, and each writes them in a way of its own (the tabular rows, `tabular.rs`;
//! outfmt 0, `report.rs`; the XML and JSON documents, `blast_xml.rs`; the Seq-aligns,
//! `report::serial::seqalign`). This module is the one owner of what those writers share: the
//! predicate of an empty row, the Seq-loc mapper's rows of `-subject`, and NCBI's exceptions.

use crate::ncbi_exception::{NcbiException, ThrowSite};

use super::preliminary::PreliminaryHsp;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1348-1377
// ```c++
//     if (hsp->query.frame == 0) {
//        query_loc->SetInt().SetFrom(hsp->query.offset);
//        query_loc->SetInt().SetTo(hsp->query.end - 1);
//     } else if (hsp->query.frame > 0) {
//        query_loc->SetInt().SetFrom(CODON_LENGTH*(hsp->query.offset) +
//                                    hsp->query.frame - 1);
//        query_loc->SetInt().SetTo(CODON_LENGTH*(hsp->query.end) +
//                                  hsp->query.frame - 2);
//     } else {
//        query_loc->SetInt().SetFrom(query_length -
//            CODON_LENGTH*(hsp->query.end) + hsp->query.frame + 1);
//        query_loc->SetInt().SetTo(query_length - CODON_LENGTH*hsp->query.offset
//                                  + hsp->query.frame);
//     }
//     ...
//     if (hsp->subject.frame == 0) {
//        subject_loc->SetInt().SetFrom(hsp->subject.offset);
//        subject_loc->SetInt().SetTo(hsp->subject.end - 1);
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:224-228
// ```c++
//             if (rng.Empty()) {
//                 NCBI_THROW(CSeqalignException, eInvalidAlignment,
//                            "CSeq_align::GetSeqRange(): Row is empty");
//             }
//             return rng;
// ```
/// Whether the Std-seg rows of the reported ungapped HSP `hsp` are empty: an HSP of length 0,
/// whose rows end one before they start (`CRange(from, from - 1)` is empty, so
/// `CSeq_align::GetSeqRange` throws on them). Every reported HSP of another length has rows of its
/// letters (gapped HSPs, and `Blast_HSPTest` deleted the ungapped ones of negative length).
#[inline]
pub(crate) fn rows_empty(hsp: &PreliminaryHsp) -> bool {
    hsp.q_end == hsp.q_start
}

/// The `__FILE__` of `Seq_align.cpp` in BLAST+'s build (the `seq` library compiles it from
/// `src/objects/seq`), as NCBI's reports name it.
const SEQ_ALIGN_CPP: &str = "src/objects/seq/../seqalign/Seq_align.cpp";

/// `CSeq_align::GetSeqRange(): Row is empty` (Seq_align.cpp:227): where a formatter asks for the
/// range of an empty row (the deflines of outfmt 0, `qcovs`, BLAST XML with `-db`).
pub(crate) fn row_is_empty() -> NcbiException {
    NcbiException::new(
        ThrowSite {
            file: SEQ_ALIGN_CPP,
            line: 227,
            kind: "CSeqalignException::eInvalidAlignment",
            function: None,
        },
        "CSeq_align::GetSeqRange(): Row is empty",
    )
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:788-803
// ```c++
//                 int width = 0;
//                 // try to determine/check the seg_len and width
//                 if (!seg_len) {
//                     width = 0;
//                     seg_len = row_len;
//                 } else {
//                     if (row_len * 3 == seg_len) {
//                         seg_len /= 3;
//                         width = 1;
//                     } else if (row_len / 3 == seg_len) {
//                         width = 3;
//                     } else if (row_len != seg_len) {
//                         NCBI_THROW(CSeqalignException, eInvalidInputAlignment,
//                                    "CreateDensegFromStdseg(): "
//                                    "Std-seg segment lengths not accurate!");
// ```
/// `CreateDensegFromStdseg(): Std-seg segment lengths not accurate!` (Seq_align.cpp:801): the
/// Std-seg of `-subject` that the Seq-loc mapper rewrote (`MappedRows`), a query row of one
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

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1541-1547
// ```c++
// static void
// s_RemapToSubjectLoc(CRef<CSeq_align> & subj_aligns, const CSeq_loc& subj_loc)
// {
//     const int kSubjDimension = 1;
//     _ASSERT(subj_loc.IsInt() || subj_loc.IsWhole());
//     subj_aligns.Reset(RemapAlignToLoc(*subj_aligns, kSubjDimension, subj_loc));
// }
// ```
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:424-427
// ```c++
//             int len = loc.GetTotalRange().GetLength();
//             if (len == 0  ||  loc.IsWhole()) {
//                 continue; // ignore unknown lengths
//             }
// ```
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1535-1557
// ```c++
//                         if ( start_frame ) {
//                             if (start >= start_frame) {
//                                 // Add 1 or 2 bases to match complete codon.
//                                 loc->SetInt().SetFuzz_from().SetAlt().push_back(start);
//                                 start -= start_frame;
//                             }
//                             ...
//                         if ( stop_frame ) {
//                             loc->SetInt().SetFuzz_to().SetAlt().push_back(stop - 1);
//                             stop += stop_frame;
//                         }
//                     }
//                 }
//
//                 loc->SetInt().SetFrom(start);
//                 // len may be 0 after dividing by width, check it before
//                 // decrementing stop.
//                 loc->SetInt().SetTo(stop > start ? stop - 1 : stop);
// ```
/// Which Seq-loc NCBI's Seq-aligns have for the subject: `-db`'s whole sequence, whose
/// Seq-aligns stay as the engine made them, or `-subject`'s interval
/// (blast_fasta_input.cpp:459-460), to which `s_RemapToSubjectLoc` maps every Seq-align with
/// `CSeq_loc_Mapper_Base`. The mapper keeps the rows of an HSP with letters; it rewrites the
/// empty rows of an HSP of length 0 (`MappedRows`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SubjectLoc {
    Whole,
    Interval,
}

impl SubjectLoc {
    /// The subjects of the report `shown`.
    pub(crate) fn of(shown: &super::report::ShownSubjects<'_>) -> Self {
        match shown {
            super::report::ShownSubjects::Fasta { .. } => SubjectLoc::Interval,
            super::report::ShownSubjects::Database { .. } => SubjectLoc::Whole,
        }
    }
}

/// The query row of the Std-seg that the Seq-loc mapper of `-subject` makes of an HSP of length
/// 0 at a subject offset above 0, in nucleotides of the record (0-based): one nucleotide, at the
/// engine's `from` on a plus frame and one before it on a minus frame (the mapper takes the
/// frames from the empty protein row and pads the nucleotide row to codons, `x_GetDstStd`); the
/// subject row stays empty. Probes: c0054 and P1's p0017, p0077, p0095, p0106 (`X/SURVEY_X2g.md`
/// 1.3).
pub(crate) fn mapped_query_nucleotide(raw_from: i32, frame: i8) -> i32 {
    if frame > 0 {
        raw_from
    } else {
        raw_from - 1
    }
}

/// LOSAT's refusal of a report of an ungapped HSP of length 0 whose form LOSAT does not
/// reproduce: not an error of NCBI's engine or formatters.
#[derive(Debug)]
pub(crate) struct EmptyRowsRefusal(pub(crate) String);

impl std::fmt::Display for EmptyRowsRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for EmptyRowsRefusal {}

/// The query's fields of a row that NCBI's formatter makes from the Seq-align's query id
/// (`align.GetSeq_id(0)`, after the query's Bioseq failed it): the record's local id (`Query_N`)
/// without `-parse_deflines`, else the id of its Seq-aligns (`BelievedIds::best`), as
/// `s_GetSeqIdListString` shows a list of that one id.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:794-798
/// ```c++
///             list<CRef<CSeq_id> > query_ids;
///             CRef<CSeq_id> id(new CSeq_id());
///             id->Assign(align.GetSeq_id(0));
///             query_ids.push_back(id);
///             SetQueryId(query_ids);
/// ```
pub(crate) struct AlignQueryIds {
    /// `qseqid`, `qacc` and `qaccver`.
    pub(crate) ids: [Vec<u8>; 3],
    /// `qgi`.
    pub(crate) gi: i64,
}

impl AlignQueryIds {
    pub(crate) fn of(query: &super::input::FastaRecord) -> Self {
        match query.believed() {
            Some(believed) => {
                let best = std::slice::from_ref(&believed.best);
                let field = |kind| {
                    let mut out = Vec::new();
                    crate::seqid::defline_set::tabular_id(best, kind, &mut out);
                    out
                };
                use crate::seqid::defline_set::TabularId;
                AlignQueryIds {
                    ids: [
                        field(TabularId::Full),
                        field(TabularId::Accession),
                        field(TabularId::AccessionVersion),
                    ],
                    gi: crate::seqid::label::first_gi(best),
                }
            }
            None => {
                let id = query.internal_id.as_bytes().to_vec();
                AlignQueryIds {
                    ids: [id.clone(), id.clone(), id],
                    gi: 0,
                }
            }
        }
    }
}

/// The refusal of `-subject` where a reported HSP of length 0 lies at the subject's first letter.
pub(crate) const UNPORTED_SUBJECT_START: &str = "BLASTX -ungapped -window_size 0 with -subject: a one-hit ungapped alignment of length 0 at the first letter of its subject among the hits to report, which NCBI BLAST+'s Seq-loc mapper cuts into two alignments without scores (tabular rows of the query's Seq-align id with an E-value of 0.0 and a bit score of -1.0, or an exception, depending on the frame and the fields), which LOSAT does not support";

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:1430-1437
// ```c++
//     CSeq_loc src_loc(*id, 0, len - 1);
//     ENa_strand strand = loc.GetStrand();
//     if (strand != eNa_strand_unknown) {
//         src_loc.SetStrand(strand);
//     }
//     CSeq_loc_Mapper_Base mapper(src_loc, loc);
//     return mapper.Map(align, row);
// ```
// No NCBI counterpart: LOSAT's refusal (decision X2g-4 of line RG). With `-subject` (also with
// `-subject_loc`) NCBI's Seq-loc mapper maps the empty subject row `(0, 4294967295)` of an HSP of
// length 0 at the first letter of the searched range through the interval of that range: one
// Std-seg of the range's length and one of the rest, both without scores
// (`X/SURVEY_X2g_formats.md` 1, 4). NCBI's formatters then write rows of the query's Seq-align
// id, throw where they read the rest's empty or wrapped rows, or run out of memory on its
// segment of about 2^32 letters (`BLAST ran out of memory`, exit status 4, which depends on the
// memory limit). The organism report (outfmt 18) does not read the rows.
/// Refuses the reports of `batches` with `-subject` (other than outfmt 18) where an HSP of length 0
/// that they show starts at the first letter of its subject's searched range.
pub(crate) fn refuse_mapped_at_subject_start(
    batches: &[super::runtime::BatchResults],
    gapped: bool,
) -> anyhow::Result<()> {
    let found = batches
        .iter()
        .flat_map(|batch| batch.queries.iter().flatten())
        .flat_map(|list| &list.hsps)
        .any(|hsp| {
            super::report::has_seqalign(hsp, gapped) && rows_empty(&hsp.hsp) && hsp.hsp.s_start == 0
        });
    if found {
        return Err(EmptyRowsRefusal(UNPORTED_SUBJECT_START.to_string()).into());
    }
    Ok(())
}
