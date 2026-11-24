//! The one producer of the tabular rows of Seq-aligns (`-outfmt 6`, 7, 10 and 20 of
//! blast_formatter, decision F3-13 of line RG): NCBI's `CBlastTabularInfo::SetFields` over a
//! Seq-align as `report::serial::seqalign_read` reads it, giving the shared writer's
//! `TabularHsp` (`report::tabular`), and `CBlastFormatUtil::InsertSubjectScores` (`qcovs`,
//! `qcovus`) over a result set's Seq-aligns. The Seq-aligns of the programs that are not
//! translated (a Dense-seg, or Dense-diags) only: a Std-seg is refused before any row
//! (`blast_formatter::refusals`).
//!
//! What a row reads of a Seq-align, in NCBI's order: its scores (`GetAlnScores`), the query's
//! coverage (`hsp_percent_coverage`, else the span over the query's length), and where a field
//! needs the alignment block its Dense-seg (a Dense-diag Seq-align's made by
//! `CreateDensegFromDendiag`), reversed where its first strand is minus, as `CAlnVec` reads it
//! (`query_anchored::aln_vec::AlnMap`: the lengths, the rows' letters and the coordinates).

use crate::report::query_anchored::aln_vec::{AlnMap, MapBuffers};
use crate::report::query_anchored::{Denseg, SeqKey, Sequences};
use crate::report::serial::seqalign::Strand;
use crate::report::serial::seqalign_read::{PairDenseg, ReadAlign, ReadStop};
use crate::report::tabular::{TabularField, TabularHsp, TabularSpec};
use crate::seqid::defline_set::{compare, Sic};
use TabularField::*;

/// What a row reads besides its Seq-align.
pub(crate) struct RowContext<'s> {
    pub(crate) spec: &'s TabularSpec,
    /// `SetNoFetch` set: the program is BLASTN (blast_format.cpp:791-792), whose rows are not
    /// read for `nident`, `mismatch` and `pident`.
    pub(crate) no_fetch: bool,
    /// The query Bioseq's length (`GetBioseqLength`), for `qcovhsp`.
    pub(crate) query_length: i64,
}

/// The vectors that the rows of a report reuse.
#[derive(Default)]
pub(crate) struct RowBuffers {
    map: MapBuffers,
    query: Vec<u8>,
    subject: Vec<u8>,
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:972-980
/// Whether the rows of a report of `spec` read the alignment's letters (`GetWholeAlnSeqString`
/// of both rows): `qseq`, `sseq`, `positive`, `ppos`, `btop`, and where the program is not
/// BLASTN (`no_fetch`) `nident`, `mismatch` and `pident`.
pub(crate) fn reads_letters(spec: &TabularSpec, no_fetch: bool) -> bool {
    spec.has_any(&[QuerySeq, SubjectSeq, Positives, PercentPositives, Btop])
        || (!no_fetch && spec.has_any(&[NumIdentical, Mismatches, PercentIdentical]))
}

/// The fields whose block reads the scores (tabular.cpp:764-773).
const SCORE_FIELDS: [TabularField; 6] = [
    Evalue,
    BitScore,
    Score,
    NumIdentical,
    Mismatches,
    PercentIdentical,
];

/// The fields of the query block (tabular.cpp:776-780).
const QUERY_FIELDS: [TabularField; 5] = [
    QuerySeqId,
    QueryGi,
    QueryAccession,
    QueryAccessionVersion,
    QueryCovHsp,
];

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:749-1096
/// The row of `align` (`SetFields`), its letters from `sequences` (the query row's and the
/// subject row's, `SeqKey::Query` and `SeqKey::Subject(0)`), which are read only where
/// `reads_letters`. The scores that a printed field needs are refused where the Seq-align has
/// none (NCBI prints its -1, which BLAST's archives never hold); so is a count of identities
/// that BLASTN prints from the Seq-align without one.
pub(crate) fn hsp<'r>(
    align: &ReadAlign<'_>,
    ctx: &RowContext<'_>,
    sequences: &dyn Sequences,
    buffers: &'r mut RowBuffers,
) -> Result<TabularHsp<'r>, ReadStop> {
    let spec = ctx.spec;
    let mut num_ident: i64 = -1;
    let (mut score, mut bit_score, mut evalue) = (0i32, 0.0f64, 0.0f64);
    if spec.has_any(&SCORE_FIELDS) {
        let scores = align.scores()?;
        let missing = (spec.has(Score) && scores.score == -1)
            || (spec.has(BitScore) && scores.bit_score == -1.0)
            || (spec.has(Evalue) && scores.evalue == -1.0);
        if missing {
            return Err(ReadStop::Unsupported(
                "a Seq-align without the score, bit score or e-value that a field prints",
            ));
        }
        score = i32::try_from(scores.score)
            .map_err(|_| ReadStop::Unsupported("a score past the int range"))?;
        bit_score = scores.bit_score;
        evalue = scores.evalue;
        num_ident = scores.num_ident;
    }
    let counts = spec.opens_counts(&usize::try_from(num_ident).unwrap_or(0));
    let query_block = spec.has_any(&QUERY_FIELDS);
    if !counts && !query_block {
        return Ok(TabularHsp {
            evalue,
            bit_score,
            score,
            num_ident: usize::try_from(num_ident).unwrap_or(0),
            ..empty_hsp()
        });
    }
    // The coverage of the query block (`x_SetQueryCovSeqalign`) reads the Seq-align's range
    // of the query, which needs its Dense-seg.
    let ds = align.dense_seg()?;
    let qcovhsp = if query_block {
        query_coverage_per_hsp(align, &ds, ctx.query_length)?
    } else {
        0
    };
    let mut hsp = TabularHsp {
        evalue,
        bit_score,
        score,
        num_ident: usize::try_from(num_ident).unwrap_or(0),
        qcovhsp,
        ..empty_hsp()
    };
    if !counts {
        return Ok(hsp);
    }
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:943-952
    let denseg = Denseg::pair(
        [SeqKey::Query, SeqKey::Subject(0)],
        &ds.starts,
        &ds.lens,
        ds.strands,
        ds.first_minus,
    );
    let RowBuffers {
        map,
        query,
        subject,
    } = buffers;
    let aln = AlnMap::new_in(&denseg, &mut *map).map_err(ReadStop::Unsupported)?;
    if spec.has_any(&[
        AlignmentLength,
        Gaps,
        Mismatches,
        PercentPositives,
        PercentIdentical,
        GapOpenings,
    ]) {
        let (length, gaps, gap_opens) = aln.align_lengths();
        hsp.length = length as usize;
        hsp.gaps = gaps as usize;
        hsp.gap_opens = gap_opens as usize;
    }
    let letters = reads_letters(spec, ctx.no_fetch);
    if letters {
        let row = |row: usize, out: &mut Vec<u8>| -> Result<(), ReadStop> {
            let string = aln
                .whole_aln_seq_string(row, 0, sequences)
                .map_err(|_| ReadStop::Unsupported("a translated alignment row"))?;
            *out = string.buffer;
            Ok(())
        };
        row(0, query)?;
        row(1, subject)?;
    } else if spec.has(NumIdentical) && num_ident < 0 {
        return Err(ReadStop::Unsupported(
            "a BLASTN Seq-align without the count of identities that nident prints",
        ));
    }
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1029-1060
    // (not translated: the query as it is, the subject swapped where its strand, or the
    // query's, is minus.)
    hsp.q_start = i64::from(aln.seq_start(0)) + 1;
    hsp.q_end = i64::from(aln.seq_stop(0)) + 1;
    let (s_first, s_last) = (
        i64::from(aln.seq_start(1)) + 1,
        i64::from(aln.seq_stop(1)) + 1,
    );
    if ds.strands[1] == Strand::Minus || ds.strands[0] == Strand::Minus {
        (hsp.s_start, hsp.s_end) = (s_last, s_first);
    } else {
        (hsp.s_start, hsp.s_end) = (s_first, s_last);
    }
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1070-1093
    hsp.q_frame = 1;
    hsp.s_frame = if spec.has_any(&[SubjectFrame, Frames]) && hsp.s_start > hsp.s_end {
        -1
    } else {
        1
    };
    if letters {
        let (query, subject): (&'r Vec<u8>, &'r Vec<u8>) = (query, subject);
        hsp.rows = Some((query, subject));
    }
    Ok(hsp)
}

/// The values of `x_ResetFields` (0, no rows).
fn empty_hsp<'r>() -> TabularHsp<'r> {
    TabularHsp {
        q_start: 0,
        q_end: 0,
        s_start: 0,
        s_end: 0,
        evalue: 0.0,
        bit_score: 0.0,
        score: 0,
        length: 0,
        gap_opens: 0,
        gaps: 0,
        num_ident: 0,
        q_frame: 0,
        s_frame: 0,
        rows: None,
        end_gaps: false,
        stored_gaps: [&[], &[]],
        qcovhsp: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:737-747
/// `x_SetQueryCovSeqalign`: the Seq-align's `hsp_percent_coverage` cut to an `int`, else its
/// span of the query (the absolute difference of row 0's `GetSeqStop` and `GetSeqStart`, plus 1) in percent of the query's
/// length, 0.5 added below 99 (`serial::seqalign::query_coverage`). A value whose cast is
/// undefined (past the `int` range, a query of length 0) is refused.
fn query_coverage_per_hsp(
    align: &ReadAlign<'_>,
    ds: &PairDenseg,
    query_length: i64,
) -> Result<i32, ReadStop> {
    let value = match align.named_real("hsp_percent_coverage") {
        Some(value) => value,
        None => {
            let (from, to) = align.seq_range(ds, 0);
            let span = (to.wrapping_sub(from) as i32).unsigned_abs() as i64 + 1;
            let query_length = i32::try_from(query_length)
                .ok()
                .filter(|&length| length > 0)
                .ok_or(ReadStop::Unsupported(
                    "a query whose length is 0 or past 2^31",
                ))?;
            let span = i32::try_from(span)
                .map_err(|_| ReadStop::Unsupported("an alignment span past 2^31"))?;
            crate::report::serial::seqalign::query_coverage(0, span, query_length)
        }
    };
    int_of(value)
}

/// `(int)` of a `double` whose value fits (the cast of every other is undefined: refused).
fn int_of(value: f64) -> Result<i32, ReadStop> {
    if value.is_nan() || value <= f64::from(i32::MIN) - 1.0 || value >= 2147483648.0 {
        return Err(ReadStop::Unsupported("a coverage past the int range"));
    }
    Ok(value.trunc() as i32)
}

/// The coverages of each row (`qcovs`, `qcovus`; -1 for `N/A`), by Seq-align.
pub(crate) type SubjectScores = Vec<(i32, i32)>;

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:812-821
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:572-660
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:707-735
/// `InsertSubjectScores` over the pruned Seq-aligns `aligns` of a query of length
/// `query_length`, then `x_SetQueryCovSubject` and `x_SetQueryCovUniqSubject` of each row: the
/// coverages that the rows print. `qcovs` is computed where `spec_text` (the format's field
/// list as typed) holds `qcovs`, `qcovus` where it holds `qcovus` and the program is BLASTN
/// (`blastn`), each unless the first Seq-align has its score; per group of adjacent Seq-aligns
/// whose subject ids match, set on the group's first where it is not 0. A row's coverage is its
/// Seq-align's score, else `N/A` where its subject's FASTA id is not that of the last row with a
/// score, else that row's.
pub(crate) fn subject_scores(
    aligns: &[ReadAlign<'_>],
    query_length: Option<i64>,
    spec_text: &str,
    blastn: bool,
) -> Result<SubjectScores, ReadStop> {
    let names = ["seq_percent_coverage", "uniq_seq_percent_coverage"];
    let mut set: Vec<[Option<i32>; 2]> = aligns
        .iter()
        .map(|align| -> Result<_, ReadStop> {
            let mut scores = [None, None];
            for (score, name) in scores.iter_mut().zip(names) {
                *score = match align.named_int(name)? {
                    Some(value) => Some(
                        i32::try_from(value)
                            .map_err(|_| ReadStop::Unsupported("a coverage past the int range"))?,
                    ),
                    None => None,
                };
            }
            Ok(scores)
        })
        .collect::<Result<_, _>>()?;
    let mut compute = [
        spec_text.contains("qcovs"),
        spec_text.contains("qcovus") && blastn,
    ];
    if let Some(first) = set.first() {
        for (kind, score) in compute.iter_mut().zip(first) {
            *kind &= score.is_none();
        }
    }
    let query_length = query_length.unwrap_or(0);
    if compute.iter().any(|&kind| kind) && query_length > 0 {
        let query_length = i32::try_from(query_length)
            .map_err(|_| ReadStop::Unsupported("a query whose length is past 2^31"))?;
        let mut left = 0;
        while left < aligns.len() {
            let id = aligns[left].seq_id(1)?;
            let mut right = left + 1;
            while right < aligns.len() && compare(aligns[right].seq_id(1)?, id) == Sic::Yes {
                right += 1;
            }
            let group = &aligns[left..right];
            let covered = [
                if compute[0] {
                    master_coverage(group)?
                } else {
                    0
                },
                if compute[1] {
                    u64::from(uniq_seq_coverage(group)?)
                } else {
                    0
                },
            ];
            for kind in 0..2 {
                if compute[kind] && covered[kind] != 0 {
                    let covered = i32::try_from(covered[kind])
                        .map_err(|_| ReadStop::Unsupported("a coverage past 2^31 letters"))?;
                    set[left][kind] = Some(int_of(
                        crate::report::serial::seqalign::query_coverage(0, covered, query_length),
                    )?);
                }
            }
            left = right;
        }
    }
    // The rows' state (`m_QueryCovSubject`, `m_QueryCovUniqSubject`): `N/A` and -1 at first.
    let mut kept: [(Vec<u8>, i32); 2] = [(b"N/A".to_vec(), -1), (b"N/A".to_vec(), -1)];
    let mut out = Vec::with_capacity(aligns.len());
    for (align, scores) in aligns.iter().zip(&set) {
        let mut fasta = Vec::new();
        crate::seqid::label::fasta(align.seq_id(1)?, &mut fasta);
        for (state, score) in kept.iter_mut().zip(scores) {
            match score {
                Some(value) => *state = (fasta.clone(), *value),
                None if state.0 != fasta => *state = (b"N/A".to_vec(), -1),
                None => {}
            }
        }
        out.push((kept[0].1, kept[1].1));
    }
    Ok(out)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:2969-2993
/// `GetMasterCoverage`: the length of the union of the Seq-aligns' query ranges.
fn master_coverage(group: &[ReadAlign<'_>]) -> Result<u64, ReadStop> {
    let mut ranges = Vec::with_capacity(group.len());
    for align in group {
        let ds = align.dense_seg()?;
        let (from, to) = align.seq_range(&ds, 0);
        let (from, to) = if from > to { (to, from) } else { (from, to) };
        ranges.push((u64::from(from), u64::from(to)));
    }
    Ok(crate::report::tabular::union_length(&mut ranges) as u64)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4838-4896
/// `GetUniqSeqCoverage`: the query letters of the Seq-aligns of one subject, where a part of a
/// Seq-align whose subject letters a Seq-align before it covered does not count (its query
/// positions mapped through the alignment's columns: the Dense-seg as it is, a Dense-diag
/// Seq-align's made by `CreateDensegFromDendiag`, which the first Seq-align's form decides).
/// The ranges are BLASTN's (`algorithm::blastn::tabular`).
pub(crate) fn uniq_seq_coverage(group: &[ReadAlign<'_>]) -> Result<u32, ReadStop> {
    use crate::algorithm::blastn::tabular::{RangeCollection, SeqRange};
    let Some((first, rest)) = group.split_first() else {
        return Ok(0);
    };
    let variant = first.variant();
    let range = |align: &ReadAlign<'_>, row: usize| -> Result<SeqRange, ReadStop> {
        let ds = align.dense_seg()?;
        let (from, to) = align.seq_range(&ds, row);
        Ok(SeqRange::new(from, to))
    };
    let mut subject_ranges = RangeCollection(vec![range(first, 1)?]);
    let mut query_ranges = RangeCollection(vec![range(first, 0)?.fix_minus_strand()]);
    let mut buffers = MapBuffers::default();
    for align in rest {
        let align_subject = range(align, 1)?;
        let mut parts = RangeCollection(vec![align_subject]);
        for &covered in &subject_ranges.0 {
            parts.subtract(covered);
        }
        if parts.0.is_empty() {
            continue;
        }
        if parts.0[0] == align_subject {
            query_ranges.combine_with(range(align, 0)?.fix_minus_strand());
            subject_ranges.combine_with(align_subject);
            continue;
        }
        if align.variant() != variant {
            return Err(ReadStop::Unsupported(
                "the Seq-aligns of a subject in both Dense-seg and Dense-diag forms",
            ));
        }
        let ds = align.dense_seg()?;
        let denseg = Denseg::pair(
            [SeqKey::Query, SeqKey::Subject(0)],
            &ds.starts,
            &ds.lens,
            ds.strands,
            false,
        );
        let map = AlnMap::new_in(&denseg, &mut buffers).map_err(ReadStop::Unsupported)?;
        for part in parts.0 {
            // `SetFrom` and `SetTo` of a `TSignedSeqPos` into the `TSeqPos` range.
            let position = |subject: u32| {
                map.seq_pos_from_aln_pos(0, map.aln_pos_from_seq_pos_none(1, subject) as u32) as u32
            };
            let query_part = SeqRange::new(position(part.from), position(part.to()));
            query_ranges.combine_with(query_part.fix_minus_strand());
            subject_ranges.combine_with(part);
        }
    }
    Ok(query_ranges.covered_length())
}

/// The letters of a row's query and subject, plus strand, as `CSeqVector` gives them in
/// `ncbieaa` or IUPAC (`Sequences` of one row).
pub(crate) struct RowLetters<'a> {
    pub(crate) query: &'a [u8],
    pub(crate) subject: &'a [u8],
}

impl Sequences for RowLetters<'_> {
    fn letters(&self, seq: SeqKey, from: u32, to: u32, out: &mut Vec<u8>) {
        let letters = match seq {
            SeqKey::Query => self.query,
            SeqKey::Subject(_) => self.subject,
        };
        let end = (to as usize + 1).min(letters.len());
        if let Some(range) = letters.get(from as usize..end) {
            out.extend_from_slice(range);
        }
    }

    fn label(&self, _seq: SeqKey, _out: &mut Vec<u8>) {}

    fn masks(&self, _seq: SeqKey, _out: &mut Vec<crate::report::query_anchored::Mask>) {}
}

#[cfg(test)]
#[path = "tabular_seqalign_tests.rs"]
mod tests;
