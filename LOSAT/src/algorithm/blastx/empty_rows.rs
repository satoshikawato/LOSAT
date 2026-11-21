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
/// of another length has rows of its letters (gapped HSPs, by their Seq-align's ends
/// `Hsp::seqalign_ends`, and `Blast_HSPTest` deleted the ungapped ones of negative length).
#[inline]
pub(crate) fn rows_empty(hsp: &super::results::Hsp) -> bool {
    crate::report::serial::seqalign::ungapped_rows_empty(
        i64::from(hsp.hsp.q_start),
        i64::from(hsp.seqalign_ends().0),
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
        match shown {
            super::report::ShownSubjects::Fasta { .. } => SubjectLoc::Interval,
            super::report::ShownSubjects::Database { .. } => SubjectLoc::Whole,
        }
    }

    /// The Seq-loc that NCBI's formatters find in the Seq-align of the subject of the hit list
    /// `list` of a report of this Seq-loc's subjects: `-subject`'s mapper leaves the Seq-align of
    /// a subject whose every HSP has length 0 (`subject_rows_all_empty`) as the engine made it,
    /// moved to the record's coordinates, unless it rewrites one of its empty rows
    /// (`mapper_wraps_empty_row`, an HSP at the first or second letter of the searched range):
    /// `Whole`, which `-db`'s ported writers write (rows of empty ranges with the HSP's scores,
    /// NCBI's prefetch warning, `row_is_empty`). Probes (ws, NCBI BLAST+ 2.17.0): `X/x2g/rev4`
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
    /// (`std_seg_lengths_not_accurate`), after the rows before it.
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
    /// `qcovhsp`'s count (`x_SetQueryCovSeqalign` reads `GetSeqStop(0) - GetSeqStart(0)`), past
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
    /// (`row_is_empty` where they read its range: `qcovs`, the query's fields), and its length
    /// not the subject row's (`std_seg_lengths_not_accurate` where they make a Dense-seg of it).
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

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:598-606 (cited at
// `report::query_anchored::aln_mix::Mix::place`), c++/src/objtools/alnmgr/alnsegments.cpp:
// 310-327 (cited at `report::query_anchored::aln_mix::Mix::fill_unaligned_regions`)
/// The refusal of a query-anchored report (outfmt 1 to 4) whose section NCBI's merge frees a
/// segment in (`report::pairwise::EmptyRowsStop::MergeFreesSegment`,
/// `report::query_anchored::aln_mix::MixError::NcbiFreesSegment`): the Dense-seg of an HSP of
/// length 0 has no widths, so `DisplaySeqalign` merges it with the HSPs of the query's frame +1
/// (`report::query_anchored::HspStrands::no_widths`), as a segment of length 0 of a query row.
/// One whose query row starts at a codon of one of those HSPs splits that HSP's segment at a
/// length of 0, which replaces every reference to the segment while the merge still walks it
/// (`Mix::place`). One that its query row's next segment does not follow at once is where
/// `FillUnalignedRegions` starts to fill the row, which replaces the row's reference to its
/// segment: the last one unless the segment is the last that the merge made (`x_Merge`'s local
/// `CRef`), so the segment is freed while `m_Segments` still lists it and `x_CreateDenseg`
/// reads it (`Mix::fill_unaligned_regions`). PD-LOSAT-NCBI-DEFECTS rule 3.
pub(crate) const UNPORTED_ANCHORED_FREE: &str = "BLASTX -ungapped -window_size 0 with -outfmt 1 to 4: a one-hit ungapped alignment of length 0 among the alignments shown where NCBI BLAST+'s merge of the query-anchored alignments frees memory that it still reads (it crashed in our runs): one whose query row starts at the first nucleotide of a codon of an alignment in the query's reading frame +1, whose segment the merge splits at a length of 0, or one followed on its query row by unaligned letters, whose segment of length 0 the merge's FillUnalignedRegions replaces, unless it is the last segment the merge made, which LOSAT does not support";

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:422-427,499-503
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
    !gapped && aligned.peek().is_some() && aligned.all(rows_empty)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:677-679
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:226-230
/// Whether NCBI's prefetch of the query batch `batches` of a report of `loc`'s subjects (where it
/// reads the Seq-aligns, `common::database::prefetch_reads_seq_aligns`) throws before it loads
/// anything: a Seq-align whose subject row is empty (`subject_rows_all_empty`, and with
/// `-subject` one that the mapper leaves so, `SubjectLoc::of_list`; the prefetch runs with
/// `-subject` too, blastx_app.cpp:286-287), where `s_ExtractSeqidsAndRanges` asks for its range.
/// The prefetch catches it with its warning (`common::database::prefetch_warning`) and has read
/// no subject, as a batch whose prefetch is off: the reports then read the subjects one by one
/// (`-subject`'s are in the scope already).
pub(crate) fn prefetch_reads_empty_row(
    batches: &[super::runtime::BatchResults],
    gapped: bool,
    loc: SubjectLoc,
) -> bool {
    batches
        .iter()
        .flat_map(|batch| batch.queries.iter().flatten())
        .any(|list| {
            subject_rows_all_empty(list, gapped) && loc.of_list(list, gapped) == SubjectLoc::Whole
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
/// `n + 1` down (`seg.m_Len - mseg.m_Len`, the unknown length -1 as `TSeqPos`:
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
mod tests {
    use super::*;
    use crate::ncbi_exception::NCBI_BUILD_ROOT;

    // NCBI BLAST+ 2.17.0 (the official Linux x64 release), rows db.x2g.bx.f0616.0 and
    // db.x2g.bx.f0616.subj.6: `CATCH_ALL`'s report of the two exceptions.
    #[test]
    fn the_exceptions_are_ncbis_reports() {
        // NCBI's reports in LOSAT's words (`crate::wording::losat_text`).
        assert_eq!(
            crate::cli::NativeError::exception(&row_is_empty()).message,
            crate::wording::losat_text(&format!(
                "Error: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/objects/seq/../seqalign/Seq_align.cpp\", line 227: Error: (CSeqalignException::eInvalidAlignment) CSeq_align::GetSeqRange(): Row is empty\n\n"
            ))
            .as_bytes()
        );
        assert_eq!(
            crate::cli::NativeError::exception(&std_seg_lengths_not_accurate()).message,
            crate::wording::losat_text(&format!(
                "Error: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/objects/seq/../seqalign/Seq_align.cpp\", line 801: Error: (CSeqalignException::eInvalidInputAlignment) CreateDensegFromStdseg(): Std-seg segment lengths not accurate!\n\n"
            ))
            .as_bytes()
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
            title: Vec::new(),
            sequence: sequence.as_bytes().to_vec(),
            lowercase_masks: Vec::new(),
            warnings: Vec::new(),
            source: None,
            around: None,
            skipped: None,
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
        assert_eq!(plus.rows(), (&b"TCCC"[..], &b"TLWI"[..]));
        assert_eq!((plus.q_frame, plus.length()), (1, 4));
        // Its minus-frame twin (an HSP of frame -1, whose empty row starts at 12):
        // `11 8 2 5 4 CCCT TLWI ... -2`, the row two below that start.
        let MappedRows::Letters(minus) = mapped(&record("TATGGGCAGGGAGAA"), "GTLWI", -1) else {
            panic!("letters");
        };
        assert_eq!((minus.q_start, minus.q_end, minus.q_frame), (11, 8, -2));
        assert_eq!(minus.rows(), (&b"CCCT"[..], &b"TLWI"[..]));
        // x2g_s1o (o00016): one nucleotide past the query's end, the gap character.
        let MappedRows::Letters(past) = mapped(&record("TGCTCGGCATGGCTG"), "ETANLXXXXXXXXX", 1)
        else {
            panic!("letters");
        };
        assert_eq!((past.q_start, past.q_end, past.length()), (4, 16, 13));
        assert_eq!(past.rows().0, b"TCGGCATGGCTG-");
        // x2g_s1w: the minus twin would start before the query's first nucleotide.
        assert!(matches!(
            mapped(&record("CAGCCATGCCGAGCA"), "ETANLXXXXXXXXX", -1),
            MappedRows::QueryRowWraps
        ));
        // x2g_s1q: -query_loc 1-15 of 43 nucleotides; the row's 16th is the record's, which the
        // cut keeps around the interval (`query_letters_around`: the subject's 14 letters and one).
        let mut cut = record("TGCTCGGCATGGCTGACGTACGTACGTACGTACGTACGTACGT");
        cut.keep_around(0, 15, 15);
        cut.cut_in_place(0, 15, false);
        let MappedRows::Letters(read) = mapped(&cut, "ETANLXXXXXXXXX", 1) else {
            panic!("letters");
        };
        assert_eq!((read.q_start, read.q_end), (4, 16));
        assert_eq!(read.rows().0, b"TCGGCATGGCTGA");
        // Without the kept letters the record does not hold the 16th: no rows (refused before).
        let mut bare = record("TGCTCGGCATGGCTGACGTACGTACGTACGTACGTACGTACGT");
        bare.cut_in_place(0, 15, false);
        assert!(MappedRows::of(&ungapped(1, 1).hsp, &bare, &record("ETANLXXXXXXXXX")).is_none());
        // An interval starting past the record's first letter (review 4 S2): `-query_loc 5-11` of
        // 12 nucleotides, a minus-frame empty row at 8 of the record; the row ends two below it
        // and starts one before the interval (the record's 4th letter, kept around it).
        let mut minus_cut = record("TTTAGGGTACGT");
        minus_cut.keep_around(4, 11, 6);
        minus_cut.cut_in_place(4, 11, false);
        let MappedRows::Letters(shifted) = mapped(&minus_cut, "GTLWI", -1) else {
            panic!("letters");
        };
        assert_eq!((shifted.q_start, shifted.q_end), (7, 4));
        assert_eq!(shifted.rows().0, b"CCCT");
        assert_eq!(
            (0..12)
                .map(|at| minus_cut.record_letter(at))
                .collect::<Option<Vec<u8>>>(),
            Some(b"TTTAGGGTACGT".to_vec())
        );
        assert_eq!(minus_cut.record_letter(12), None);
        // The farthest a minus row reads before the interval (review 5 N3): an HSP at the end of
        // frame -1 of `-query_loc 7-9` (its empty row starts at the interval's first nucleotide)
        // and at the first letter of a subject of 5 letters reads from 6 (5 + 1) letters before
        // the interval: `query_letters_around`'s one more than the subject's letters.
        let farthest = |around: usize| {
            let mut cut = record("ACGTACGGGTT");
            cut.keep_around(6, 9, around);
            cut.cut_in_place(6, 9, false);
            let mut hsp = ungapped(1, 1).hsp;
            (hsp.frame, hsp.s_start, hsp.s_end) = (-1, 0, 0);
            MappedRows::of(&hsp, &cut, &record("GTLWI"))
        };
        let Some(MappedRows::Letters(far)) = farthest(6) else {
            panic!("letters");
        };
        assert_eq!((far.q_start, far.q_end, far.length()), (5, 1, 5));
        assert_eq!(far.rows().0, b"TACGT");
        assert!(farthest(5).is_none());
    }

    // Review 4 M1 (ws, NCBI BLAST+ 2.17.0, `X/x2g/rev4`): `x_ConvertSegment`'s stop wraps, and the
    // mapper rewrites an empty row, only at the first or second letter of a subject without a row
    // of letters, and only at the first of a protein subject (X2g-4).
    #[test]
    fn the_mapper_rewrites_an_empty_row_only_where_its_stop_wraps() {
        assert!(mapper_wraps_empty_row(0, false));
        assert!(mapper_wraps_empty_row(1, false));
        assert!(!mapper_wraps_empty_row(2, false));
        assert!(!mapper_wraps_empty_row(7, false));
        assert!(mapper_wraps_empty_row(0, true));
        assert!(!mapper_wraps_empty_row(1, true));
        let at = |s_start: i32| {
            let mut hsp = ungapped(1, 1);
            hsp.hsp.s_start = s_start;
            hsp.hsp.s_end = s_start;
            hsp
        };
        let list = |hsps| super::super::results::HspList {
            query_index: 0,
            oid: 0,
            hsps,
            best_evalue: 1.4,
        };
        use SubjectLoc::{Interval, Whole};
        assert_eq!(Interval.of_list(&list(vec![at(1)]), false), Interval);
        assert_eq!(Interval.of_list(&list(vec![at(3)]), false), Whole);
        assert_eq!(Interval.of_list(&list(vec![at(3), at(5)]), false), Whole);
        assert_eq!(Interval.of_list(&list(vec![at(3), at(0)]), false), Interval);
        // A subject with letters keeps the run's Seq-loc.
        assert_eq!(
            Interval.of_list(&list(vec![at(3), ungapped(2, 5)]), false),
            Interval
        );
        assert_eq!(Whole.of_list(&list(vec![at(1)]), false), Whole);
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
