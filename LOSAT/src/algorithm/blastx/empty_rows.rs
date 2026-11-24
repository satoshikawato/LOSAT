//! The reports of BLASTX's ungapped HSPs of length 0, and of the rows that the reports cannot
//! take the range of. The one-hit word finder of `-ungapped -window_size 0` with a low
//! `-threshold` saves HSPs whose end lies at or before their start (aa_ungapped.c:1054,1083);
//! `Blast_HSPTest` deletes those of negative length and keeps those of length 0 (`runtime.rs`,
//! `core::blast_hits::hsp_test`), which are reported. Their Std-seg rows end one before they start
//! (`to = from - 1`, `x_UngappedHSPToStdSeg`): NCBI's formatters find them empty, and each writes
//! them in a way of its own (the tabular rows, `tabular.rs`; outfmt 0, `report.rs`; the XML and
//! JSON documents, `blast_xml.rs`; the Seq-aligns, `report::serial::seqalign`). An alignment of
//! RPSTBLASTN's Smith-Waterman redo (`-use_sw_tback`, series PR item PS12d) can start with a gap in
//! the profile's row, whose range they cannot take either (`row_range`). This module is the one
//! owner of what those writers share: the predicates of an empty row and of a row's range, the
//! Seq-loc mapper's rows of `-subject`, and NCBI's exceptions (`Seq_align.cpp`'s,
//! `report::serial::seqalign::row_range_exception` and `std_seg_lengths_not_accurate`, and the
//! pairwise writer's stops, `report::pairwise::empty_rows_stop`, which TBLASTN's writers share).

use crate::report::pairwise::EmptyRowsRefusal;
use crate::report::serial::seqalign::{AlignRow, RowRangeError};

use super::preliminary::PreliminaryHsp;

/// Whether the Std-seg rows of the reported ungapped HSP `hsp` are empty: an HSP of length 0
/// (`serial::seqalign::ungapped_rows_empty`, the one predicate of the writers). Every reported HSP
/// of another length has rows of its letters (gapped HSPs, by their Seq-align's ends
/// `Hsp::seqalign_ends`, and `Blast_HSPTest` deleted the ungapped ones of negative length).
#[inline]
pub(crate) fn rows_empty(hsp: &super::results::Hsp) -> bool {
    crate::report::serial::seqalign::ungapped_rows_empty(
        i64::from(hsp.hsp.q_start),
        i64::from(hsp.seqalign_ends().0),
    )
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:171-191
/// What `CSeq_align::GetSeqRange` of the row finds on the Std-segs of the reported HSP `hsp`
/// (`serial::seqalign::std_seg_row_range`, the one predicate of the writers): the empty rows of an
/// HSP of length 0 (`rows_empty`), or the row that a gapped HSP's edit script starts with a gap in,
/// as the engine made them (`-db`; `-subject`'s Seq-loc mapper rewrites some, `SubjectLoc`).
#[inline]
pub(crate) fn row_range(hsp: &super::results::Hsp, row: AlignRow) -> Option<RowRangeError> {
    crate::report::serial::seqalign::std_seg_row_range(
        i64::from(hsp.hsp.q_start),
        i64::from(hsp.seqalign_ends().0),
        &hsp.edit_script,
        row,
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:2907-2930
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:243-245
/// What a formatter that reads the query's range of `hsp`'s Seq-align and then the subject's finds
/// (`row_range`): the first exception, `None` where both are taken. The description table of
/// outfmt 0 to 4 (`s_ProcessAlignSet`) and BLAST XML (`s_SeqAlignToXMLHsp`) read them so; an HSP of
/// length 0 has both rows empty, a gapped one at most one row that starts with a gap.
#[inline]
pub(crate) fn rows_range(hsp: &super::results::Hsp) -> Option<RowRangeError> {
    row_range(hsp, AlignRow::Query).or_else(|| row_range(hsp, AlignRow::Subject))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1541-1547
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:425-428
/// Which Seq-loc NCBI's Seq-aligns have for the subject: `-db`'s whole sequence, whose
/// Seq-aligns stay as the engine made them, or `-subject`'s interval
/// (blast_fasta_input.cpp:459-460), to which `s_RemapToSubjectLoc` maps every Seq-align with
/// `CSeq_loc_Mapper_Base`. The mapper keeps the rows of an HSP with letters; it rewrites the
/// empty rows of an HSP of length 0 (`serial::seqalign::mapped_empty_query_nucleotide`), except
/// those of a subject whose every HSP has length 0 that its arithmetic leaves as they are
/// (`SubjectLoc::of_list`): such a subject's Seq-align is `Whole`'s, in the record's coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SubjectLoc {
    Whole,
    Interval,
}

impl SubjectLoc {
    /// The subjects of the report `shown`.
    pub(crate) fn of(shown: &super::report::ShownSubjects<'_>) -> Self {
        Self::of_run(matches!(
            shown,
            super::report::ShownSubjects::Database { .. }
        ))
    }

    /// The subjects of a run that searches a database (`-db`), or `-subject`'s records: the one
    /// test of whether the run's Seq-aligns are mapped (`BatchResults::put_in_report_order`).
    pub(crate) fn of_run(database: bool) -> Self {
        if database {
            SubjectLoc::Whole
        } else {
            SubjectLoc::Interval
        }
    }

    /// The Seq-loc that NCBI's formatters find in the Seq-align of the subject of the hit list
    /// `list` of a report of this Seq-loc's subjects: `-subject`'s mapper leaves the Seq-align of
    /// a subject whose every HSP has length 0 (`subject_rows_all_empty`) as the engine made it,
    /// moved to the record's coordinates, unless it rewrites one of its empty rows
    /// (`mapper_wraps_empty_row`, an HSP at the first or second letter of the searched range):
    /// `Whole`, which `-db`'s ported writers write (rows of empty ranges with the HSP's scores,
    /// NCBI's prefetch warning, `row_range_exception`). Probes (ws, NCBI BLAST+ 2.17.0): `X/x2g/rev4`
    /// (review 4 M1: 248 pairs with the HSP at offset 2 or more, NCBI `-subject` = NCBI `-db` in
    /// every table and outfmt 18), `rev4c` (the offset counts from the `-subject_loc` interval's
    /// start).
    pub(crate) fn of_list(self, list: &super::results::HspList, gapped: bool) -> Self {
        let unmapped = || {
            subject_rows_all_empty(list, gapped)
                && !list
                    .hsps
                    .iter()
                    .filter(|hsp| super::report::has_seqalign(hsp, gapped))
                    .any(|hsp| mapper_wraps_empty_row(hsp.hsp.s_start, false))
        };
        match self {
            SubjectLoc::Interval if unmapped() => SubjectLoc::Whole,
            loc => loc,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1067-1068
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1134-1146
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:539-545
/// Whether `-subject`'s Seq-loc mapper (`x_ConvertSegment`) rewrites the empty subject row of an
/// HSP of length 0 that starts at `s_start` of the subject's searched range. The row's segment
/// has the unknown length -1 (`seglens`, seq_align_mapper_base.cpp:499), so its stop is two below
/// its start in the unsigned `TSeqPos`: only where that wraps does `dr` trim the stop to the
/// interval's end and the mapped segment get letters and lose its scores. Elsewhere the mapped
/// segment keeps the length -1 and the scores, and `x_GetDstStd` writes the empty rows again
/// (`to = from - 1`, :1556-1559). The start is in the mapper's units: the subject's letters,
/// times 3 where the mapper knows that the subject is a protein (`protein`: another Std-seg of
/// the Seq-align has rows of letters, `x_Init`'s first pass, :454-496); a Seq-align without a row
/// of letters gives it no type (both rows are taken as nucleotides). So a subject with letters
/// wraps at its first letter only (decision X2g-4 of line RG), a subject of empty rows only at
/// its first or second (decision X2g-5).
fn mapper_wraps_empty_row(s_start: i32, protein: bool) -> bool {
    /// `seglens`' -1 of a Std-seg whose rows have no letters.
    const UNKNOWN_LENGTH: i32 = -1;
    let width: u32 = if protein { 3 } else { 1 };
    let start = (s_start as u32).wrapping_mul(width);
    let stop = start.wrapping_add(UNKNOWN_LENGTH as u32).wrapping_sub(1);
    stop >= start
}

/// How NCBI's tabular formatter (`CBlastTabularInfo::SetFields`) writes the row of an HSP with
/// empty rows (`rows_empty`), for BLASTX's tabular writer (`tabular::write_blastx_tabular`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum TabularEmptyRow {
    /// `-subject` with a field of the alignment block, a subject with two or more shown HSPs (or
    /// the mapped rows of a subject's only HSP whose query row wraps, `MappedRows::QueryRowWraps`,
    /// without a query field): `CreateDensegFromStdseg` on the mapped Std-seg
    /// (`serial::seqalign::std_seg_lengths_not_accurate`), after the rows before it.
    Throws,
    /// `-db` (and `-subject`'s Seq-aligns as made, `SubjectLoc::of_list`), a subject whose only
    /// shown HSP it is, with a query field (or `-subject`'s mapped
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
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3170-3184
/// The row form of an HSP with empty rows of a subject of `loc` (the hit list's,
/// `SubjectLoc::of_list`) with `subject_hsps` shown HSPs, where the
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
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1034-1043
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
    /// The rows (`GetWholeAlnSeqString`), of `n` letters each: the query row reads the record's
    /// letters, also those of a record cut to `-query_loc` outside its interval
    /// (`input::FastaRecord::record_letter`, which the cut keeps for it, `query_letters_around`).
    pub(crate) qrow: Vec<u8>,
    pub(crate) srow: Vec<u8>,
}

impl MappedAlignment {
    /// The alignment's length, also the query row's span in nucleotides, `qcovs`' and
    /// `qcovhsp`'s count (`x_SetQueryCovSeqalign` reads `GetSeqStop` of row 0 less `GetSeqStart` of row 0), past
    /// the record's end too.
    pub(crate) fn length(&self) -> usize {
        self.srow.len()
    }

    /// The rows.
    pub(crate) fn rows(&self) -> (&[u8], &[u8]) {
        (&self.qrow, &self.srow)
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
    /// (`row_range_exception` where they read its range: `qcovs`, the query's fields), and its length
    /// not the subject row's (`serial::seqalign::std_seg_lengths_not_accurate` where they make a Dense-seg of it).
    /// Probes: `fix3_oob.py`'s minus-strand pairs.
    QueryRowWraps,
}

impl MappedRows {
    /// The mapped rows of the HSP `hsp` of length 0 of the query `query` and the subject
    /// `subject` (cut to their `-query_loc` / `-subject_loc` intervals), whose empty subject row
    /// the mapper rewrites (`mapper_wraps_empty_row`); `None` for an HSP past its subject's end,
    /// which the engine does not report, and for a query letter that the record does not hold
    /// (`input::FastaRecord::record_letter`; its cut keeps every letter the rows can read,
    /// `query_letters_around`).
    pub(crate) fn of(
        hsp: &PreliminaryHsp,
        query: &super::input::FastaRecord,
        subject: &super::input::FastaRecord,
    ) -> Option<Self> {
        let s_from = usize::try_from(hsp.s_start).ok()?;
        let srow = subject.sequence.get(s_from..)?.to_vec();
        let n = srow.len();
        let last = n.checked_sub(1)?;
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
        // NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1175-1183
        // The segment's length is the unknown -1 (`seglens`, :499) as an unsigned `TSeqPos`, the
        // mapped one `n`: a plus row (the strand of the subject's unset one) keeps its start, a
        // minus row's moves `n + 1` down, so that it ends two before its start (probes: qstart
        // one below dev's on the minus frames of every `fix3_grid.py` case).
        let (from, to) = if plus {
            (start, start + last)
        } else {
            match start.checked_sub(n + 1) {
                Some(from) => (from, start - 2),
                None => return Some(MappedRows::QueryRowWraps),
            }
        };
        // The record's letters; the gap character past its end.
        let letter = |at: usize| -> Option<u8> {
            if at >= place.length {
                return Some(b'-');
            }
            query.record_letter(at)
        };
        let qrow = if plus {
            (from..=to).map(letter).collect::<Option<Vec<u8>>>()?
        } else {
            (from..=to)
                .rev()
                .map(|at| letter(at).map(crate::utils::genetic_code::complement))
                .collect::<Option<Vec<u8>>>()?
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
/// `subject`, of `loc` (the hit list's, `SubjectLoc::of_list`): a `-subject` subject whose only
/// shown HSP has length 0 and whose empty row the mapper rewrites; `None` for every other list
/// (and where `MappedRows::of` finds none, which `refuse_mapped_without_scores` refused before
/// the report).
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

/// The query's fields of a row that NCBI's formatter makes from the Seq-align's query id
/// (`GetSeq_id` of row 0 of the alignment, after the query's Bioseq failed it): the record's local id (`Query_N`)
/// without `-parse_deflines`, else the id of its Seq-aligns (`BelievedIds::best`), as
/// `s_GetSeqIdListString` shows a list of that one id.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:793-797
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

/// The refusal of `-subject` where a reported HSP of length 0 lies at the first letter of the
/// searched part of a subject with HSPs of letters (decision X2g-4, `mapper_wraps_empty_row`).
pub(crate) const UNPORTED_SUBJECT_START: &str = "BLASTX -ungapped -window_size 0 with -subject: a one-hit ungapped alignment of length 0 at the first letter of the searched part of a subject with alignments of letters among the hits to report, which NCBI BLAST+'s Seq-loc mapper cuts into two alignments without scores (tabular rows of the query's Seq-align id with an E-value of 0.0 and a bit score of -1.0, or an exception, depending on the frame and the fields), which LOSAT does not support";

/// The refusal of `-subject` where every reported HSP of a subject has length 0 and the mapper
/// rewrites one (`mapper_wraps_empty_row`).
pub(crate) const UNPORTED_SUBJECT_ALL_EMPTY: &str = "BLASTX -ungapped -window_size 0 with -subject: a subject among the hits to report whose every one-hit ungapped alignment has length 0, one at the first or second letter of its subject's searched part, of which NCBI BLAST+'s Seq-loc mapper, finding no sequence type in its rows, makes alignments without scores of nucleotides against the rest of the subject's letters (rows with an E-value of 0.0 and a bit score of -1.0, or, for two or more such alignments, more rows, an exception or an out-of-memory stop), which LOSAT does not support where the report reads the alignment";

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:422-427,499-503
/// Whether every Std-seg of the Seq-align of the hit list `list` (one for each HSP with a Seq-align,
/// `report::has_seqalign`, of an ungapped search) has empty rows (`rows_empty`), and there is one:
/// a subject whose every reported HSP has length 0. Its Seq-align has no row of letters, so
/// `CSeq_align::GetSeqRange` finds its subject row empty (`RowRangeError::Empty`; a row of length 0 at
/// subject offset 0 ends at `TSeqPos` 4294967295, whose open end wraps to 0: empty too), and
/// `-subject`'s Seq-loc mapper finds neither sequence type in it.
pub(crate) fn subject_rows_all_empty(list: &super::results::HspList, gapped: bool) -> bool {
    let mut aligned = list
        .hsps
        .iter()
        .filter(|hsp| super::report::has_seqalign(hsp, gapped))
        .peekable();
    !gapped && aligned.peek().is_some() && aligned.all(rows_empty)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:677-679
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:226-230
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:171-191
/// Whether NCBI's prefetch of the query batch `batches` of a report of `loc`'s subjects (where it
/// reads the Seq-aligns, `common::database::prefetch_reads_seq_aligns`) throws before it loads
/// anything: `s_ExtractSeqidsAndRanges` asks for the subject row's range of every Seq-align
/// (`row_range`), which throws on an empty one, of a subject whose every HSP of an ungapped search
/// has length 0 (`subject_rows_all_empty`, and with `-subject` one that the mapper leaves so,
/// `SubjectLoc::of_list`; the prefetch runs with `-subject` too, blastx_app.cpp:286-287), and on
/// one that starts with a gap, of a gapped HSP as the engine made it (`-db`: RPSTBLASTN's
/// Smith-Waterman redo, also among the HSPs that left its lists past the formatter's cut,
/// `BatchResults::unlisted_subject_row_throws`; `-subject`'s mapper leaves out a first segment
/// with a gap row, `serial::seqalign::mapped_script`). The prefetch catches it with its warning
/// (`common::database::prefetch_warning`) and has read no subject, as a batch whose prefetch is
/// off: the reports then read the subjects one by one (`-subject`'s are in the scope already).
pub(crate) fn prefetch_reads_unranged_row(
    batches: &[super::runtime::BatchResults],
    gapped: bool,
    loc: SubjectLoc,
) -> bool {
    batches.iter().any(|batch| {
        (loc == SubjectLoc::Whole && batch.unlisted_subject_row_throws)
            || batch.queries.iter().flatten().any(|list| {
                if gapped {
                    loc == SubjectLoc::Whole
                        && list.hsps.iter().any(|hsp| {
                            super::report::has_seqalign(hsp, gapped)
                                && row_range(hsp, AlignRow::Subject).is_some()
                        })
                } else {
                    subject_rows_all_empty(list, gapped)
                        && loc.of_list(list, gapped) == SubjectLoc::Whole
                }
            })
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1352-1358
/// The letters on each side of a `-query_loc` interval that a query record cut to it keeps
/// (`input::FastaRecord::keep_around`) for `-subject`'s Seq-loc mapper (`MappedRows::of`), in a
/// search of `options` against the subjects `subjects` (each cut to its `-subject_loc` interval):
/// a mapped query row has as many nucleotides, `n`, as its subject has letters from the HSP on
/// (at most the longest subject's, `L`). Its empty row starts in the interval, or just past its
/// last nucleotide for an HSP at the end of a plus frame (`from` above). A plus row runs up from
/// that start: at most `L` letters past the interval. A minus row ends two below it, as the mapper moves its start
/// `n + 1` down (`seg.m_Len` less `mseg.m_Len`, the unknown length -1 as `TSeqPos`:
/// seq_align_mapper_base.cpp:1175-1183, in `MappedRows::of`): at most `L + 1` letters before the
/// interval, for an HSP at the end of its minus frame (its `from` is the interval's first
/// nucleotide) and the first letter of the longest subject. NCBI's scope holds the whole record.
/// 0 where no row is mapped: a query searched whole, and a search that reports no HSP of length
/// 0 (gapped, or the two-hit word finder of a `-window_size` above 0: `s_BlastAaExtendOneHit`
/// alone makes them, aa_ungapped.c:1054).
pub(crate) fn query_letters_around(
    options: &super::args::ResolvedOptions,
    subjects: &[super::input::FastaRecord],
) -> usize {
    if options.query_range.is_none() || options.gapped || options.window_size != 0 {
        return 0;
    }
    subjects
        .iter()
        .map(|subject| subject.sequence.len())
        .max()
        .map_or(0, |longest| longest + 1)
}

// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:1431-1437
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
// - X2g-5: a subject whose every HSP has length 0 (`subject_rows_all_empty`), one of them at the
//   first or second letter of the searched range. No row of its Std-segs has letters, so the
//   mapper finds no sequence type (both rows are taken as nucleotides,
//   `seq_align_mapper_base.cpp:422-503`) and each segment's length stays -1: the subject row of
//   such an HSP maps from its start to the end of the interval (`x_ConvertSegment`, where `dr`
//   trims the wrapped stop, `mapper_wraps_empty_row`), so one HSP becomes an alignment of that
//   many query nucleotides against the subject's letters, without scores (rows with an E-value of
//   0.0 and a bit score of -1.0, identities where a nucleotide letter equals an amino acid
//   letter), and two or more give more rows, exceptions and, on the 4 GiB limit of the line's
//   probes, `BLAST ran out of memory` (exit status 4). Where every HSP of such a subject lies
//   further in, the mapper leaves the Seq-align as made (`SubjectLoc::of_list`): `-db`'s output,
//   not refused.
// The organism report (outfmt 18) does not read the rows; the tabular reports write X2g-5's one
// alignment as NCBI's formatter reads it (`MappedAlignment`, decision X2g-7).
/// The tabular reports of a `-subject` run (`refuse_mapped_without_scores`): the queries and
/// subjects of their rows.
#[derive(Clone, Copy)]
pub(crate) struct MappedTables<'a> {
    pub(crate) queries: &'a [super::input::FastaRecord],
    pub(crate) subjects: &'a [super::input::FastaRecord],
}

/// Refuses the reports of `batches` with `-subject` (other than outfmt 18) where every HSP that they
/// show of a subject has length 0 and the mapper rewrites one (X2g-5, `SubjectLoc::of_list`),
/// unless there is one such HSP and the report is a table (`tabular`, `MappedRows::of`), or where,
/// among the shown HSPs of a subject with letters in its rows, the mapper rewrites the empty row
/// of one of length 0 (X2g-4: at the first letter of the subject's searched range,
/// `mapper_wraps_empty_row`). The first refusal in the order of the hit lists.
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
                // Left as made where the mapper rewrites none of them: `-db`'s output.
                if SubjectLoc::Interval.of_list(list, gapped) == SubjectLoc::Whole {
                    continue;
                }
                let first = shown.next().expect("a shown HSP");
                let written = shown.next().is_none()
                    && tabular.is_some_and(|tables| {
                        MappedRows::of(
                            &first.hsp,
                            &tables.queries[batch.query_ordinal + qidx],
                            &tables.subjects[list.oid],
                        )
                        .is_some()
                    });
                (!written).then_some(UNPORTED_SUBJECT_ALL_EMPTY)
            } else {
                shown
                    .any(|hsp| rows_empty(hsp) && mapper_wraps_empty_row(hsp.hsp.s_start, true))
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
#[path = "empty_rows_tests.rs"]
mod tests;
