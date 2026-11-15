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

/// Whether the Std-seg rows of the reported ungapped HSP `hsp` are empty: an HSP of length 0
/// (`serial::seqalign::ungapped_rows_empty`, the one predicate of the writers). Every reported HSP
/// of another length has rows of its letters (gapped HSPs, and `Blast_HSPTest` deleted the
/// ungapped ones of negative length).
#[inline]
pub(crate) fn rows_empty(hsp: &PreliminaryHsp) -> bool {
    crate::report::serial::seqalign::ungapped_rows_empty(
        i64::from(hsp.q_start),
        i64::from(hsp.q_end),
    )
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
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:425-428
// ```c++
//             int len = loc.GetTotalRange().GetLength();
//             if (len == 0  ||  loc.IsWhole()) {
//                 continue; // ignore unknown lengths
//             }
// ```
/// Which Seq-loc NCBI's Seq-aligns have for the subject: `-db`'s whole sequence, whose
/// Seq-aligns stay as the engine made them, or `-subject`'s interval
/// (blast_fasta_input.cpp:459-460), to which `s_RemapToSubjectLoc` maps every Seq-align with
/// `CSeq_loc_Mapper_Base`. The mapper keeps the rows of an HSP with letters; it rewrites the
/// empty rows of an HSP of length 0 (`serial::seqalign::mapped_empty_query_nucleotide`).
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

/// How NCBI's tabular formatter (`CBlastTabularInfo::SetFields`) writes the row of an HSP with
/// empty rows (`rows_empty`): the one decision of BLASTX's two tabular writers (the field lists,
/// `tabular::write_blastx_tabular`, and the default fields of outfmt 6 and 7, `report.rs`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum TabularEmptyRow {
    /// `-subject` with a field of the alignment block, a subject with two or more shown HSPs (or
    /// the mapped rows of a subject's only HSP whose query row wraps, `MappedRows::QueryRowWraps`,
    /// without a query field): `CreateDensegFromStdseg` on the mapped Std-seg
    /// (`std_seg_lengths_not_accurate`), after the rows before it.
    Throws,
    /// `-db`, a subject whose only shown HSP it is, with a query field (or `-subject`'s mapped
    /// rows whose query row wraps): the row that leaves at the query's fields (the query's
    /// Seq-align id, `AlignQueryIds`; the scores, `empty_row_scores`, and the subject's fields;
    /// `qcovhsp` `N/A`; every field of the alignment block at `x_ResetFields`' 0).
    Leaves,
    /// The row as of any HSP (`-db`: the empty ranges; `-subject` without a field of the
    /// alignment block: its scores and ids, the mapped nucleotide in `qcovs` and `qcovhsp`).
    Row,
    /// `-subject`, a subject whose only shown HSP it is: the row of the Seq-loc mapper's
    /// alignment without scores (`MappedAlignment`, decision X2g-7 of line RG).
    Mapped,
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:737-746
// ```c++
// void CBlastTabularInfo::x_SetQueryCovSeqalign(const CSeq_align & align, int query_len)
// {
// 	double tmp = 0;
// 	if(!align.GetNamedScore("hsp_percent_coverage", tmp)) {
// 		int len = abs((int) (align.GetSeqStop(0) - align.GetSeqStart(0))) + 1;
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3170-3184
// ```c++
//             if(seg.GetStd().size() > 1){
//                 //has more than one stdseg. Need to seperate as each
//                 //is a distinct HSP
//                 ITERATE (CSeq_align::C_Segs::TStd, iterStdseg, seg.GetStd()){
//                     CRef<CSeq_align> aln(new CSeq_align);
//                     if((*iterStdseg)->IsSetScores()){
//                         aln->SetScore() = (*iterStdseg)->GetScores();
//                     }
//                     aln->SetSegs().SetStd().push_back(*iterStdseg);
//                     alnSetRef->Set().push_back(aln);
//                 }
//
//             } else {
//                 alnSetRef->Set().push_back(*iter);
//             }
// ```
/// The row form of an HSP with empty rows of a subject with `subject_hsps` shown HSPs, where the
/// report's fields enter NCBI's alignment block (`opens_block`, `TabularSpec::opens_counts`,
/// tabular.cpp:900-908) and read the query (`reads_query`, `TabularSpec::reads_query`,
/// tabular.cpp:778-800). A subject's only HSP keeps its Seq-align unsplit, with the scores on its
/// Std-seg and none named on the Seq-align, so that `x_SetQueryCovSeqalign` reads the query
/// row's range, which throws on `-db`'s empty row, caught by `SetFields`. With `-subject` a
/// subject's only HSP of length 0 has a Seq-align of no row of letters, which the mapper makes
/// into an alignment of letters without scores (`MappedRows`), unless its query row wraps
/// (`query_row_wraps`): then the query row is empty as `-db`'s, and its length not the subject
/// row's.
pub(crate) fn tabular_empty_row(
    loc: SubjectLoc,
    opens_block: bool,
    subject_hsps: usize,
    reads_query: bool,
    query_row_wraps: bool,
) -> TabularEmptyRow {
    match loc {
        SubjectLoc::Interval if subject_hsps == 1 && !query_row_wraps => TabularEmptyRow::Mapped,
        SubjectLoc::Interval if subject_hsps == 1 && reads_query => TabularEmptyRow::Leaves,
        SubjectLoc::Interval if opens_block => TabularEmptyRow::Throws,
        SubjectLoc::Whole if subject_hsps == 1 && reads_query => TabularEmptyRow::Leaves,
        _ => TabularEmptyRow::Row,
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:820-833
// ```c++
//     bool hasScore = false;
//     score = -1;
//     bits = -1;
//     evalue = -1;
//     sum_n = -1;
//     num_ident = -1;
//     comp_adj_method = 0;
//
//     //look for scores at seqalign level first
//     hasScore = s_GetBlastScore(aln.GetScore(), score, bits, evalue,
//                                sum_n, num_ident, use_this_gi, comp_adj_method);
//
//     //look at the seg level
//     if(!hasScore){
// ```
/// The scores (`evalue`, `bitscore`, `score`) of the tabular row of an HSP of length 0 of a
/// subject with `subject_hsps` shown HSPs: the HSP's, but `-subject`'s mapped rows of a subject's
/// only HSP have none (`MappedRows`), which `GetAlnScores` reads as -1 each (an E-value of -1
/// shows as `0.0`).
pub(crate) fn empty_row_scores(
    loc: SubjectLoc,
    subject_hsps: usize,
    hsp: &super::results::Hsp,
) -> (f64, f64, i32) {
    if loc == SubjectLoc::Interval && subject_hsps == 1 {
        (-1.0, -1.0, -1)
    } else {
        (hsp.evalue, hsp.bit_score, hsp.hsp.score)
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:279-297
// ```c++
//             if (start >= 0) {
//                 stop = start + len - 1;
//
//                 // add regular sequence to buffer
//                 GetSeqString(buff, row, start, stop);
//                 TSeqPos buf_len = min<TSeqPos>(TSeqPos(buff.size()), seg_len);
//                 buffer += buff;
//                 if (buf_len < seg_len) {
//                     // Not enough chars in the sequence, add gap
//                     buf_len = seg_len - buf_len;
//                     char fill_ch;
//
//                     if (seg < left_seg  ||  seg > right_seg) {
//                         fill_ch = GetEndChar();
//                     } else {
//                         fill_ch = GetGapChar(row);
//                     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1034-1043
// ```c++
//         if (kTranslated && ds.GetSeqStrand(kQueryRow) == eNa_strand_minus) {
//             q_start = alnVec->GetSeqStop(kQueryRow) + 1;
//             q_end = alnVec->GetSeqStart(kQueryRow) + 1;
//         } else {
//             q_start = alnVec->GetSeqStart(kQueryRow) + 1;
//             q_end = alnVec->GetSeqStop(kQueryRow) + 1;
//         }
// ```
/// The alignment that `-subject`'s Seq-loc mapper makes of the only shown HSP of a subject when
/// it has length 0 (decision X2g-5 of line RG; its tabular rows ported by X2g-7). The Std-seg has
/// no row of letters, so the mapper takes both rows as nucleotides of one segment of unknown
/// length and maps the subject row from its start to the end of the subject's interval
/// (`refuse_mapped_without_scores`, `seq_align_mapper_base.cpp:422-503`, `:1002-1232`): `n`
/// residues of the subject against `n` nucleotides of the query on the frame's strand
/// (`MappedRows`), without scores.
/// NCBI's tabular formatter reads it as any Seq-align of letters (`CreateDensegFromStdseg` finds
/// equal lengths, width 1): the rows are the query's nucleotides (the complement on a minus
/// strand) and the subject's residues, which it compares as letters, and past the record's last
/// nucleotide the query row has the gap character (`alnvec.cpp:279-297`). Probes (ws, NCBI BLAST+
/// 2.17.0): `X/x2g/fix3` (`fix3_grid.py`), the base `c01241_q c01241_s 25.000 4 3 0 4 7 2 5 0.0
/// -1.0`, the minus twin `0.000 4 4 0 11 8 2 5`.
pub(crate) struct MappedAlignment {
    /// `qstart` and `qend` in the record (1-based; above `qend` on a minus strand).
    pub(crate) q_start: i64,
    pub(crate) q_end: i64,
    /// `sstart` and `send` in the record.
    pub(crate) s_start: i64,
    pub(crate) s_end: i64,
    /// `qframe`: `GetFrame` of `qstart` in the record (`seq_range::record_frame`).
    pub(crate) q_frame: i32,
    /// The rows (`GetWholeAlnSeqString`), of `n` letters each; no query row where it reads
    /// letters of a query record cut to `-query_loc` that its cut did not keep
    /// (`input::FastaRecord::cut_in_place`): letters LOSAT does not have.
    pub(crate) qrow: Option<Vec<u8>>,
    pub(crate) srow: Vec<u8>,
}

impl MappedAlignment {
    /// The alignment's length, also the query row's span in nucleotides, `qcovs`' and
    /// `qcovhsp`'s count (`x_SetQueryCovSeqalign` reads `GetSeqStop(0) - GetSeqStart(0)`), past
    /// the record's end too.
    pub(crate) fn length(&self) -> usize {
        self.srow.len()
    }

    /// The rows, where LOSAT has the query row's letters.
    pub(crate) fn rows(&self) -> Option<(&[u8], &[u8])> {
        Some((self.qrow.as_deref()?, &self.srow))
    }
}

/// What `-subject`'s Seq-loc mapper makes of the only shown HSP of a subject when it has length 0
/// (`MappedAlignment`): `n` residues of the subject from the HSP's start to the end of its
/// interval against `n` nucleotides of the query (in the record: `RemapToQueryLoc` moves the row
/// there first, blast_seqalign.cpp:1634-1636), from the start of its empty row up on a plus
/// strand, and ending two below it on a minus strand.
pub(crate) enum MappedRows {
    /// An alignment of letters.
    Letters(MappedAlignment),
    /// A minus strand's row that would start before the record's first nucleotide: its start, an
    /// unsigned `TSeqPos`, wraps, so NCBI's formatters find the query row empty as `-db`'s
    /// (`row_is_empty` where they read its range: `qcovs`, the query's fields), and its length
    /// not the subject row's (`std_seg_lengths_not_accurate` where they make a Dense-seg of it).
    /// Probes: `fix3_oob.py`'s minus-strand pairs.
    QueryRowWraps,
}

impl MappedRows {
    /// The mapped rows of the HSP `hsp` of length 0 of the query `query` and the subject
    /// `subject` (cut to their `-query_loc` / `-subject_loc` intervals); `None` for an HSP past
    /// its subject's end, which the engine does not report.
    pub(crate) fn of(
        hsp: &PreliminaryHsp,
        query: &super::input::FastaRecord,
        subject: &super::input::FastaRecord,
    ) -> Option<Self> {
        let s_from = usize::try_from(hsp.s_start).ok()?;
        let srow = subject.sequence.get(s_from..)?.to_vec();
        let n = srow.len();
        let place = query.placement();
        let searched = query.sequence.len();
        // The query row's start in the record: `RemapToQueryLoc` moves the row first.
        let start = usize::try_from(
            crate::report::serial::seqalign::ungapped_interval(
                hsp.q_start,
                hsp.q_end,
                hsp.frame,
                searched as i32,
            )
            .0,
        )
        .ok()?
            + place.offset;
        let plus = hsp.frame > 0;
        // NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1176-1184
        // ```c++
        //             else {
        //                 // Adjust starts of all other rows.
        //                 if (mrow.m_Start != kInvalidSeqPos) {
        //                     if (mrow.SameStrand(aln_row)) {
        //                         mrow.m_Start += left_shift;
        //                     }
        //                     else {
        //                         mrow.m_Start +=
        //                             seg.m_Len - mseg.m_Len - left_shift;
        // ```
        // The segment's length is the unknown -1 (`seglens`, :499) as an unsigned `TSeqPos`, the
        // mapped one `n`: a plus row (the strand of the subject's unset one) keeps its start, a
        // minus row's moves `n + 1` down, so that it ends two before its start (probes: qstart
        // one below dev's on the minus frames of every `fix3_grid.py` case).
        let (from, to) = if plus {
            (start, start + n.checked_sub(1)?)
        } else {
            match start.checked_sub(n + 1) {
                Some(from) => (from, start - 2),
                None => return Some(MappedRows::QueryRowWraps),
            }
        };
        // The record's letters: the searched ones, or the whole record that a cut kept (lower case
        // where it was masked); the gap character past its end.
        let whole = query
            .source
            .as_deref()
            .and_then(|source| source.uncut.as_deref());
        let letter = |at: usize| -> Option<u8> {
            if at >= place.length {
                return Some(b'-');
            }
            match whole {
                Some(whole) => whole.get(at).map(u8::to_ascii_uppercase),
                None => at
                    .checked_sub(place.offset)
                    .and_then(|at| query.sequence.get(at).copied()),
            }
        };
        let qrow = if plus {
            (from..=to).map(letter).collect::<Option<Vec<u8>>>()
        } else {
            (from..=to)
                .rev()
                .map(|at| letter(at).map(complement_iupacna))
                .collect::<Option<Vec<u8>>>()
        };
        let (q_start, q_end) = if plus { (from, to) } else { (to, from) };
        let s_offset = subject.placement().offset;
        Some(MappedRows::Letters(MappedAlignment {
            q_start: q_start as i64 + 1,
            q_end: q_end as i64 + 1,
            s_start: (s_offset + s_from) as i64 + 1,
            s_end: (s_offset + s_from + n) as i64,
            q_frame: crate::blastinput::seq_range::record_frame(plus, q_start, place.length),
            qrow,
            srow,
        }))
    }
}

/// The mapped rows (`MappedRows`) of the hit list `list` of the query `query` against the subject
/// `subject` of a report of `loc`'s subjects: a `-subject` subject whose only shown HSP has
/// length 0; `None` for every other list, and where LOSAT does not have the letters
/// (`MappedRows::of`, which `refuse_mapped_without_scores` refused before the report).
pub(crate) fn subject_mapped_rows(
    loc: SubjectLoc,
    list: &super::results::HspList,
    gapped: bool,
    query: &super::input::FastaRecord,
    subject: &super::input::FastaRecord,
) -> Option<MappedRows> {
    if loc != SubjectLoc::Interval || !subject_rows_all_empty(list, gapped) {
        return None;
    }
    let mut shown = list
        .hsps
        .iter()
        .filter(|hsp| super::report::has_seqalign(hsp, gapped));
    let only = shown.next()?;
    if shown.next().is_some() {
        return None;
    }
    MappedRows::of(&only.hsp, query, subject)
}

/// The complement of an IUPAC nucleotide letter (or `-`), as NCBI reads a minus strand: the
/// reversed bits of its `ncbi4na` value (`blastx::report::display_query_residue`).
fn complement_iupacna(letter: u8) -> u8 {
    const IUPACNA: &[u8; 16] = b"-ACMGRSVTWYHKDBN";
    match IUPACNA.iter().position(|&c| c == letter) {
        Some(mask) => {
            let mask = mask as u8;
            IUPACNA[(((mask & 1) << 3) | ((mask & 2) << 1) | ((mask & 4) >> 1) | ((mask & 8) >> 3))
                as usize]
        }
        None => letter,
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
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:793-797
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

/// The refusal of `-subject` where every reported HSP of a subject has length 0.
pub(crate) const UNPORTED_SUBJECT_ALL_EMPTY: &str = "BLASTX -ungapped -window_size 0 with -subject: a subject among the hits to report whose every one-hit ungapped alignment has length 0, of which NCBI BLAST+'s Seq-loc mapper, finding no sequence type in its rows, makes alignments without scores of nucleotides against the rest of the subject's letters (rows with an E-value of 0.0 and a bit score of -1.0, or, for two or more such alignments, more rows, an exception or an out-of-memory stop), which LOSAT does not support where the report reads the alignment";

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:422-427,499-503
// ```c++
//         ITERATE( CStd_seg::TLoc, it_loc, (*it)->GetLoc()) {
//             const CSeq_loc& loc = **it_loc;
//             const CSeq_id* id = loc.GetId();
//             int len = loc.GetTotalRange().GetLength();
//             if (len == 0  ||  loc.IsWhole()) {
//                 continue; // ignore unknown lengths
//             }
//     ...
//         seglens.push_back(maxlen == 0 ? -1 : maxlen);
//     }
//     // By this point all possible sequence types should be detected and
//     // stored in the loc-mapper.
//     // All unknown types are treated as nucs.
// ```
/// Whether every Std-seg of the Seq-align of the hit list `list` (one for each HSP with a Seq-align,
/// `report::has_seqalign`, of an ungapped search) has empty rows (`rows_empty`), and there is one:
/// a subject whose every reported HSP has length 0. Its Seq-align has no row of letters, so
/// `CSeq_align::GetSeqRange` finds its subject row empty (`row_is_empty`; a row of length 0 at
/// subject offset 0 ends at `TSeqPos` 4294967295, whose open end wraps to 0: empty too), and
/// `-subject`'s Seq-loc mapper finds neither sequence type in it.
pub(crate) fn subject_rows_all_empty(list: &super::results::HspList, gapped: bool) -> bool {
    let mut aligned = list
        .hsps
        .iter()
        .filter(|hsp| super::report::has_seqalign(hsp, gapped))
        .peekable();
    !gapped && aligned.peek().is_some() && aligned.all(|hsp| rows_empty(&hsp.hsp))
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:677-679
// ```c++
//         ITERATE(CSeq_align_set::Tdata, aln, (*result)->GetSeqAlign()->Get()) {
//             CConstRef<CSeq_id> subj(&(*aln)->GetSeq_id(kSubjRow));
//             TSeqRange subj_range((*aln)->GetSeqRange(kSubjRow));
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:226-230
// ```c++
//             if (rng.Empty()) {
//                 NCBI_THROW(CSeqalignException, eInvalidAlignment,
//                            "CSeq_align::GetSeqRange(): Row is empty");
//             }
//             return rng;
// ```
/// Whether NCBI's prefetch of the query batch `batches` (where it reads the Seq-aligns,
/// `common::database::prefetch_reads_seq_aligns`) throws before it loads anything: a Seq-align
/// whose subject row is empty (`subject_rows_all_empty`), where `s_ExtractSeqidsAndRanges` asks
/// for its range. The prefetch catches it with its warning (`common::database::prefetch_warning`)
/// and has read no subject, as a batch whose prefetch is off: the reports then read the subjects
/// one by one.
pub(crate) fn prefetch_reads_empty_row(
    batches: &[super::runtime::BatchResults],
    gapped: bool,
) -> bool {
    batches
        .iter()
        .flat_map(|batch| batch.queries.iter().flatten())
        .any(|list| subject_rows_all_empty(list, gapped))
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:1431-1437
// ```c++
//     CSeq_loc src_loc(*id, 0, len - 1);
//     ENa_strand strand = loc.GetStrand();
//     if (strand != eNa_strand_unknown) {
//         src_loc.SetStrand(strand);
//     }
//     CSeq_loc_Mapper_Base mapper(src_loc, loc);
//     return mapper.Map(align, row);
// ```
// No NCBI counterpart: LOSAT's refusals (decisions X2g-4 and X2g-5 of line RG). With `-subject`
// (also with `-subject_loc`) NCBI's Seq-loc mapper maps two kinds of Seq-aligns of HSPs of length
// 0 into alignments without scores (`X/SURVEY_X2g_formats.md` 1, 4; `X/x2g/REPORT_impl2.md`, ws
// gate):
// - X2g-4: the empty subject row `(0, 4294967295)` of an HSP of length 0 at the first letter of the
//   searched range, mapped through the interval of that range: one Std-seg of the range's length
//   and one of the rest. NCBI's formatters then write rows of the query's Seq-align id, throw
//   where they read the rest's empty or wrapped rows, or run out of memory on its segment of
//   about 2^32 letters (`BLAST ran out of memory`, exit status 4, which depends on the memory
//   limit).
// - X2g-5: a subject whose every HSP has length 0 (`subject_rows_all_empty`). No row of its
//   Std-segs has letters, so the mapper finds no sequence type (both rows are taken as
//   nucleotides, `seq_align_mapper_base.cpp:422-503`) and each segment's length stays -1: the
//   subject row maps from its start to the end of the interval (`x_ConvertSegment`,
//   `seq_align_mapper_base.cpp:1002-1232`, where `dr` trims the wrapped stop), so one HSP becomes
//   an alignment of that many query nucleotides against the subject's letters, without scores
//   (rows with an E-value of 0.0 and a bit score of -1.0, identities where a nucleotide letter
//   equals an amino acid letter), and two or more give more rows, exceptions and, on the 4 GiB
//   limit of the line's probes, `BLAST ran out of memory` (exit status 4).
// The organism report (outfmt 18) does not read the rows; the tabular reports write X2g-5's one
// alignment as NCBI's formatter reads it (`MappedAlignment`, decision X2g-7).
/// The tabular reports of a `-subject` run (`refuse_mapped_without_scores`): the queries and
/// subjects of their rows, and whether the rows read the alignment's letters
/// (`TabularSpec::reads_row_letters`; the default fields' `pident` and `mismatch` do).
#[derive(Clone, Copy)]
pub(crate) struct MappedTables<'a> {
    pub(crate) queries: &'a [super::input::FastaRecord],
    pub(crate) subjects: &'a [super::input::FastaRecord],
    pub(crate) reads_letters: bool,
}

/// Refuses the reports of `batches` with `-subject` (other than outfmt 18) where every HSP that they
/// show of a subject has length 0 (X2g-5), unless there is one such HSP, the report is a table
/// (`tabular`) and LOSAT has the letters of its mapped rows that the table reads
/// (`MappedRows::of`), or where, among the shown HSPs of a subject with letters in its rows,
/// one of length 0 starts at the first letter of the subject's searched range (X2g-4). The first
/// refusal in the order of the hit lists.
pub(crate) fn refuse_mapped_without_scores(
    batches: &[super::runtime::BatchResults],
    gapped: bool,
    tabular: Option<MappedTables<'_>>,
) -> anyhow::Result<()> {
    for batch in batches {
        for (qidx, list) in batch
            .queries
            .iter()
            .enumerate()
            .flat_map(|(qidx, lists)| lists.iter().map(move |list| (qidx, list)))
        {
            let mut shown = list
                .hsps
                .iter()
                .filter(|hsp| super::report::has_seqalign(hsp, gapped));
            // X2g-5 first: a subject of empty rows only has no sequence type for the mapper, also
            // where one of them lies at the first letter.
            let refusal = if subject_rows_all_empty(list, gapped) {
                let first = shown.next().expect("a shown HSP");
                let written = shown.next().is_none()
                    && tabular.is_some_and(|tables| {
                        match MappedRows::of(
                            &first.hsp,
                            &tables.queries[batch.query_ordinal + qidx],
                            &tables.subjects[list.oid],
                        ) {
                            Some(MappedRows::Letters(mapped)) => {
                                mapped.qrow.is_some() || !tables.reads_letters
                            }
                            Some(MappedRows::QueryRowWraps) => true,
                            None => false,
                        }
                    });
                (!written).then_some(UNPORTED_SUBJECT_ALL_EMPTY)
            } else {
                shown
                    .any(|hsp| rows_empty(&hsp.hsp) && hsp.hsp.s_start == 0)
                    .then_some(UNPORTED_SUBJECT_START)
            };
            if let Some(refusal) = refusal {
                return Err(EmptyRowsRefusal(refusal.to_string()).into());
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ncbi_exception::NCBI_BUILD_ROOT;

    // NCBI BLAST+ 2.17.0 (the official Linux x64 release), rows db.x2g.bx.f0616.0 and
    // db.x2g.bx.f0616.subj.6: `CATCH_ALL`'s report of the two exceptions.
    #[test]
    fn the_exceptions_are_ncbis_reports() {
        assert_eq!(
            crate::cli::NativeError::exception(&row_is_empty()).message,
            format!(
                "Error: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/objects/seq/../seqalign/Seq_align.cpp\", line 227: Error: (CSeqalignException::eInvalidAlignment) CSeq_align::GetSeqRange(): Row is empty\n\n"
            )
        );
        assert_eq!(
            crate::cli::NativeError::exception(&std_seg_lengths_not_accurate()).message,
            format!(
                "Error: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/objects/seq/../seqalign/Seq_align.cpp\", line 801: Error: (CSeqalignException::eInvalidInputAlignment) CreateDensegFromStdseg(): Std-seg segment lengths not accurate!\n\n"
            )
        );
    }

    fn ungapped(q_start: i32, q_end: i32) -> super::super::results::Hsp {
        super::super::results::Hsp {
            hsp: PreliminaryHsp {
                context: 0,
                frame: 1,
                score: 5,
                q_start,
                q_end,
                q_gapped_start: q_start,
                s_start: 1,
                s_end: 1 + q_end - q_start,
                s_gapped_start: 1,
            },
            num: 1,
            evalue: 1.4,
            bit_score: 5.2,
            identity: 0,
            positive: 0,
            edit_script: Vec::new(),
            composition_method: 0,
        }
    }

    // The ws gate's S1 rows (db.x2g.bx.s1.*): a subject whose every HSP has length 0 has a
    // Seq-align of no row of letters (NCBI's prefetch writes its warning, `-subject`'s mapper
    // finds no sequence type); one HSP of another length gives it a row of letters.
    #[test]
    fn a_subject_of_empty_rows_only_needs_every_hsp_of_length_0() {
        let list = |hsps| super::super::results::HspList {
            query_index: 0,
            oid: 0,
            hsps,
            best_evalue: 1.4,
        };
        assert!(subject_rows_all_empty(&list(vec![ungapped(1, 1)]), false));
        assert!(subject_rows_all_empty(
            &list(vec![ungapped(1, 1), ungapped(4, 4)]),
            false
        ));
        assert!(!subject_rows_all_empty(
            &list(vec![ungapped(1, 1), ungapped(2, 5)]),
            false
        ));
        assert!(!subject_rows_all_empty(&list(Vec::new()), false));
        // A gapped search has no ungapped HSP of length 0.
        assert!(!subject_rows_all_empty(&list(vec![ungapped(1, 1)]), true));
    }

    // X/SURVEY_X2g_formats.md 3 (set C): -subject throws only with a field of the alignment
    // block; -db leaves at the query's fields only for a subject's only HSP with a query field.
    // A -subject subject's only HSP of length 0 (decisions X2g-5 and X2g-7): the mapper's
    // alignment of letters, no exception.
    #[test]
    fn the_tabular_row_of_empty_rows_follows_the_mode_and_the_fields() {
        use SubjectLoc::{Interval, Whole};
        use TabularEmptyRow::*;
        assert_eq!(tabular_empty_row(Interval, true, 2, true, false), Throws);
        assert_eq!(tabular_empty_row(Interval, true, 3, false, false), Throws);
        assert_eq!(tabular_empty_row(Interval, false, 2, true, false), Row);
        assert_eq!(tabular_empty_row(Interval, true, 1, true, false), Mapped);
        assert_eq!(tabular_empty_row(Interval, false, 1, false, false), Mapped);
        // The mapped query row wraps (fix3_oob.py's minus-strand pairs): as `-db`'s empty row
        // with a query field, else the Dense-seg's 801 with a field of the block.
        assert_eq!(tabular_empty_row(Interval, true, 1, true, true), Leaves);
        assert_eq!(tabular_empty_row(Interval, true, 1, false, true), Throws);
        assert_eq!(tabular_empty_row(Interval, false, 1, false, true), Row);
        assert_eq!(tabular_empty_row(Whole, true, 1, true, false), Leaves);
        assert_eq!(tabular_empty_row(Whole, false, 1, true, false), Leaves);
        assert_eq!(tabular_empty_row(Whole, true, 2, true, false), Row);
        assert_eq!(tabular_empty_row(Whole, true, 1, false, false), Row);
    }

    fn record(sequence: &str) -> super::super::input::FastaRecord {
        super::super::input::FastaRecord {
            internal_id: "Query_1".to_string(),
            title: String::new(),
            sequence: sequence.as_bytes().to_vec(),
            lowercase_masks: Vec::new(),
            warnings: Vec::new(),
            source: None,
        }
    }

    fn mapped(query: &super::super::input::FastaRecord, subject: &str, frame: i8) -> MappedRows {
        let mut hsp = ungapped(1, 1).hsp;
        hsp.frame = frame;
        MappedRows::of(&hsp, query, &record(subject)).expect("an HSP inside its subject")
    }

    // Decision X2g-7, NCBI BLAST+ 2.17.0 on ws (`X/x2g/fix3`, fix3_grid.py): the rows of the
    // mapped alignment of a subject's only HSP of length 0 at query offset 1 and subject offset 1.
    #[test]
    fn the_mapped_rows_are_ncbis() {
        // x2g_s1 (c01241, frame 1): `c01241_q c01241_s 4 7 2 5 4 TCCC TLWI ... 1`.
        let MappedRows::Letters(plus) = mapped(&record("TTCTCCCTGCCCATA"), "GTLWI", 1) else {
            panic!("letters");
        };
        assert_eq!(
            (plus.q_start, plus.q_end, plus.s_start, plus.s_end),
            (4, 7, 2, 5)
        );
        assert_eq!(plus.rows(), Some((&b"TCCC"[..], &b"TLWI"[..])));
        assert_eq!((plus.q_frame, plus.length()), (1, 4));
        // Its minus-frame twin (an HSP of frame -1, whose empty row starts at 12):
        // `11 8 2 5 4 CCCT TLWI ... -2`, the row two below that start.
        let MappedRows::Letters(minus) = mapped(&record("TATGGGCAGGGAGAA"), "GTLWI", -1) else {
            panic!("letters");
        };
        assert_eq!((minus.q_start, minus.q_end, minus.q_frame), (11, 8, -2));
        assert_eq!(minus.rows(), Some((&b"CCCT"[..], &b"TLWI"[..])));
        // x2g_s1o (o00016): one nucleotide past the query's end, the gap character.
        let MappedRows::Letters(past) = mapped(&record("TGCTCGGCATGGCTG"), "ETANLXXXXXXXXX", 1)
        else {
            panic!("letters");
        };
        assert_eq!((past.q_start, past.q_end, past.length()), (4, 16, 13));
        assert_eq!(past.rows().unwrap().0, b"TCGGCATGGCTG-");
        // x2g_s1w: the minus twin would start before the query's first nucleotide.
        assert!(matches!(
            mapped(&record("CAGCCATGCCGAGCA"), "ETANLXXXXXXXXX", -1),
            MappedRows::QueryRowWraps
        ));
        // x2g_s1q: -query_loc 1-15 of 43 nucleotides; the row's 16th is the record's, not kept.
        let mut cut = record("TGCTCGGCATGGCTGACGTACGTACGTACGTACGTACGTACGT");
        cut.cut_in_place(0, 15, false);
        let MappedRows::Letters(unread) = mapped(&cut, "ETANLXXXXXXXXX", 1) else {
            panic!("letters");
        };
        assert_eq!((unread.q_start, unread.q_end), (4, 16));
        assert_eq!(unread.rows(), None);
        let mut kept = record("TGCTCGGCATGGCTGACGTACGTACGTACGTACGTACGTACGT");
        kept.cut_in_place(0, 15, true);
        let MappedRows::Letters(read) = mapped(&kept, "ETANLXXXXXXXXX", 1) else {
            panic!("letters");
        };
        assert_eq!(read.rows().unwrap().0, b"TCGGCATGGCTGA");
        assert_eq!(complement_iupacna(b'R'), b'Y');
        assert_eq!(complement_iupacna(b'N'), b'N');
    }

    // X/SURVEY_X2g_formats.md 1: p0095 (frame -1, raw query row (119, 118)) maps to 118..118;
    // a plus frame keeps its start.
    #[test]
    fn the_mapped_nucleotide_follows_the_strand() {
        use crate::report::serial::seqalign::mapped_empty_query_nucleotide;
        assert_eq!(mapped_empty_query_nucleotide(119, -1), 118);
        assert_eq!(mapped_empty_query_nucleotide(81, 1), 81);
        assert_eq!(mapped_empty_query_nucleotide(0, 1), 0);
    }
}
