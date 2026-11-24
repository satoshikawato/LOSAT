//! NCBI's `CAlnMix` (`objtools/alnmgr`) for the flags of the query-anchored reports: `Add` with
//! `fPreserveRows` or `fForceTranslation` (`super::AddFlags`), `Merge` with `fMinGap`,
//! `fQuerySeqMergeOnly` and `fFillUnalignedRegions` (showalign.cpp:2058-2083, cited at
//! `super::write_section`), on Dense-segs of two rows whose first row is the query. The result is
//! the merged Dense-seg of `x_CreateDenseg`.
//!
//! Ported: `CAlnMixSequences::Add` (alnseq.cpp:105-196), `CAlnMixMatches::Add` (alnmatch.cpp:
//! 97-275), `CAlnMixMerger::x_Merge` (alnmerger.cpp:101-830), `x_SecondRowFits` up to its
//! `fQuerySeqMergeOnly` return (:835-860), `BuildRows` (alnseq.cpp:266-283), `Build`
//! (alnsegments.cpp:58-280), `x_MinimizeGaps` (:438-523), `FillUnalignedRegions` (:286-338) and
//! `x_CreateDenseg` (alnmerger.cpp:1156-1266).
//!
//! Not ported, because these flags never reach them: `fTruncateOverlaps` (the planes and the rest
//! of `x_SecondRowFits`), `fAllowTranslocation`, `fGapJoin` (`x_ConsolidateGaps`), `fCalcScore`,
//! `fSortSeqsByScore` and `fSortInputByScore` (the chain scores), `fRemoveLeadTrailGaps`, and the
//! scope-less identification by Seq-id.
//!
//! One Bioseq on both rows (line RG T3c fix 1). A subject that is the query's own Bioseq
//! (`Hsps::set_query_bioseq`: one record read as the query and as a subject under
//! `-parse_deflines`) is, with `fPreserveRows` (BLASTP, BLASTN), the query's `m_AnotherRow`: a
//! sequence of its own. With `fForceTranslation` (TBLASTX) `Add` makes no other row, so the HSP's
//! two rows are the query's one `CAlnMixSeq` (`Mix::add`): no sequence is then in every Dense-seg
//! as often as the mix has Dense-segs (no single refseq, `Mix::merge`), the record's subject
//! positions are placed on the query's own row or, from the second such HSP on, on its extra rows
//! (`Mix::merge_match`), and the rows after the first of that sequence are shown as the subject
//! (`Mix::create_denseg`). A merge that NCBI's code stops with an exception (`CAlnException`,
//! which `DisplaySeqalign` does not catch) is `MixError::NcbiThrows`, with NCBI's exception
//! (`X_MERGE_START_ITERATOR`, `BUILD_INFINITE_LOOP`).
//!
//! A query whose Bioseq is not a protein (line PS, PS15b-2 fix 1): the PSSM query of the molecule
//! `other` of psiblast, whose Dense-segs are those of a protein query and protein subjects, is a
//! nucleotide sequence to the mix (`m_IsAA` is `IsProtein`), so its rows have the width 3 and the
//! subjects' 1 (`super::Hsps::mix_widths`): each match goes to the query's row of its start
//! modulo 3 (`x_SetSeqFrame`), where it covers three positions for each of its letters. Two
//! matches of one HSP whose spans meet on such a row put the HSP's subject row twice on a
//! segment (`X_MERGE_START_ITERATOR`), and segments shared by HSPs on two such rows can follow
//! each other in two orders (`BUILD_INFINITE_LOOP`): `MixError::NcbiThrows`. A match of the
//! query's letters alone (a gap in the subject) whose start modulo 3 is a frame where no match of
//! two rows opened a row of the query opens a row whose first match has no second row, and NCBI
//! writes through that null row (`MixError::NcbiNullSecondRow`). The matches of two rows come
//! first (`SortByScore`), so where one of them throws, that exception is NCBI's end; `Build`
//! comes after `x_Merge`. `super::write_unreadable_section` refuses each of them.
//!
//! Strands and widths (lines RG T3b and T3c). A row's strand is set where it gets its first start.
//! The query's (the refseq of the first pass) is plus with `fPreserveRows`
//! (`fNegativeStrand` is not in `m_MergeFlags`: no merge flag has it) and, with `fForceTranslation`, the
//! sign of its strand score (the sum of its matches' lengths, negative on the minus strand), so a
//! BLASTX or TBLASTX query of a minus frame has a minus-strand row. A subject's row takes the
//! query's strand, or the other where the HSP's strands differ (`m_StrandsDiffer`), so the merge's
//! branches for rows of the other strand run. A nucleotide sequence translated against proteins
//! (BLASTX's query, TBLASTN's subjects, TBLASTX's both: `m_ContainsAA` and `m_ContainsNA` are set, or
//! `fForceTranslation`) has rows of width 3 in their frames (`x_SetSeqFrame`). With
//! `fForceTranslation` an HSP whose strand on a sequence differs from the sequence's strand score
//! so far is left out (`CAlnMixMatches::Add` throws, `Mix::add`).

use std::collections::{BTreeMap, BTreeSet};

use super::{AddFlags, Denseg, Hsps, SeqKey};
use crate::ncbi_exception::{NcbiException, ThrowSite};

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:787-794 (cited at
// `Mix::merge_match`)
// NCBI BLAST+ 2.17.0 reports it as `Error: NCBI C++ Exception:\n    T0 "<build
// tree>/src/objtools/alnmgr/alnmerger.cpp", line 790: Error: (CAlnException::eMergeFailure)
// ncbi::objects::CAlnMixMerger::x_Merge() - CAlnMixMerger::x_Merge(): Internal error: Start
// iterator already exists for seq2.` (exit 255; the review of line RG T3c fix 1, `T/t3/review2_t3c/
// probes/buf/file.err`).
/// Where `x_Merge` throws when the walk of a match's second row creates a start whose segment
/// already has that row's start iterator (`Mix::merge_match`).
const X_MERGE_START_ITERATOR: ThrowSite = ThrowSite {
    file: "src/objtools/alnmgr/alnmerger.cpp",
    line: 790,
    kind: "CAlnException::eMergeFailure",
    function: Some("ncbi::objects::CAlnMixMerger::x_Merge"),
};

/// The message of `X_MERGE_START_ITERATOR`.
const X_MERGE_START_ITERATOR_MESSAGE: &str =
    "CAlnMixMerger::x_Merge(): Internal error: Start iterator already exists for seq2.";

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnsegments.cpp:162-169 (cited at
// `Mix::build`)
// NCBI BLAST+ 2.17.0 reports it at line 168, in `ncbi::objects::CAlnMixSegments::Build()` (exit
// 255; the review of line RG T3c fix 1, case `mg0.sx.f1` of `T/t3/review2_t3c/cases_mg.tsv`).
/// Where `Build` throws when the segments on the left of a refseq segment lead back to a segment
/// already on its stack (`Mix::build`).
const BUILD_INFINITE_LOOP: ThrowSite = ThrowSite {
    file: "src/objtools/alnmgr/alnsegments.cpp",
    line: 168,
    kind: "CAlnException::eMergeFailure",
    function: Some("ncbi::objects::CAlnMixSegments::Build"),
};

/// The message of `BUILD_INFINITE_LOOP`.
const BUILD_INFINITE_LOOP_MESSAGE: &str =
    "CAlnMixSegments::Build(): Internal error: Infinite loop detected.";

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:349-376 (cited at
// `Mix::merge_match`)
// NCBI BLAST+ 2.17.0 reports it at line 372, in `ncbi::objects::CAlnMixMerger::x_Merge()` (exit
// 255; line PS, PS15b-2 fix 1, `PS/tmp/fix15b2_t0/t2a_s1.1.err`: psiblast `-in_pssm` of the
// molecule `other`).
/// Where `x_Merge` throws when no row of a subject's chain is the row of the Dense-seg of a match
/// of a later pass (`Mix::merge_match`).
const X_MERGE_PROPER_ROW: ThrowSite = ThrowSite {
    file: "src/objtools/alnmgr/alnmerger.cpp",
    line: 372,
    kind: "CAlnException::eMergeFailure",
    function: Some("ncbi::objects::CAlnMixMerger::x_Merge"),
};

/// The message of `X_MERGE_PROPER_ROW`.
const X_MERGE_PROPER_ROW_MESSAGE: &str = "CAlnMixMerger::x_Merge(): Proper row not found for the \
     match. Cannot use fQuerySeqMergeOnly?";

/// An index into `Mix::seqs` (a `CAlnMixSeq`).
type SeqIx = usize;
/// An index into `Mix::segs` (a `CAlnMixSegment`).
type SegIx = usize;
/// A `CAlnMixStarts` iterator as the key it points to; `None` is `end()`.
type It = Option<u32>;

/// Why the merge stops.
#[derive(Debug)]
pub(crate) enum MixError {
    /// An internal error: a state that NCBI's code reaches only with input that a BLAST search
    /// does not make (each named where it is raised).
    Internal(&'static str),
    /// NCBI aborts: `BuildRows` left out a row that a segment holds, and `FillUnalignedRegions`
    /// indexes its vectors with that row's `m_RowIdx` of -1 (`Mix::fill_unaligned_regions`).
    /// Only the query's own Bioseq gets there (`Mix::place`); the caller refuses the section.
    NcbiAborts,
    /// NCBI's merge throws `CAlnException` (`eMergeFailure`: this exception), which
    /// `DisplaySeqalign` does not catch (only `Add` is in its `try`, showalign.cpp:2060-2083, quoted
    /// at `super::write_section`): the exception leaves the formatter before the section's first
    /// byte, and the application's `CATCH_ALL` ends the run with `Error: NCBI C++ Exception` (exit
    /// 255). Only the query's own Bioseq added with `fForceTranslation` (TBLASTX) gets there
    /// (`Mix::merge_match`, `Mix::build`), and a query read at the width 3 whose letters are a
    /// protein's (module doc); `super::write_section` fails with it
    /// (`NcbiException::into_io_error`), `super::write_unreadable_section` refuses it.
    NcbiThrows(NcbiException),
    /// NCBI frees a segment that it still reads, where a write into a row's map replaces the last
    /// `CRef` of a segment (a `CAlnMixStarts` holds its segments by `CRef`, `m_Segments` by raw
    /// pointer, alnsegments.hpp:67,80,89), in one of two places:
    /// - The merge's split (`Mix::place`): a match of length 0 (BLASTX's HSP of length 0,
    ///   `serial::seqalign::empty_row_segments`) that starts at the start of a longer segment of
    ///   its row, or inside one (whose split gives it that start), "splits" that segment at its own
    ///   length of 0. The new segment's start on each row of the old one is the old one's:
    ///   `seg` stored in `SetStarts()` at `start` replaces the `CRef` in that row's map, and on a row of the
    ///   first row's strand (every row of such a segment of a BLASTX section) nothing else holds
    ///   the old segment, so the replacement on its last row frees it while `ITERATE` still walks
    ///   its `m_StartIts` (alnmerger.cpp:598-614).
    /// - `FillUnalignedRegions` (`Mix::fill_unaligned_regions`): the fill after a segment of
    ///   length 0 of a row starts at that segment's own start, so `seg` stored in
    ///   `row->SetStarts()` at `new_start` (alnsegments.cpp:323) replaces the row's `CRef` of it;
    ///   unless it is the last segment
    ///   that `x_Merge` created (`Mix::merge_seg`), the only other `CRef` that such a segment has
    ///   there (`Mix::held`), it is freed while `m_Segments` still lists it, and `x_CreateDenseg`
    ///   reads it (alnmerger.cpp:1196-1210).
    ///
    /// A read of freed memory (SIGSEGV in our runs). The program refuses the report before its
    /// first byte (PD-LOSAT-NCBI-DEFECTS rule 3, `super::section_frees`).
    NcbiFreesSegment,
    /// NCBI writes through a null pointer: the first match placed on a row has no second row, and
    /// `x_Merge` stores the match's Dense-seg into that row's `m_DsIdx` (`fQuerySeqMergeOnly`,
    /// `Mix::place`), so official NCBI BLAST+ 2.17.0 ends with SIGSEGV. Only a query read at the
    /// width 3 whose letters are a protein's gets there (a PSSM query of the molecule `other`,
    /// `super::Hsps::mix_widths`): a match of the query's letters alone (a gap in the subject)
    /// whose start modulo 3 is a frame where no match of two rows opened a row of the query, so
    /// `x_SetSeqFrame` gives it a new row (module doc). `super::write_unreadable_section` refuses
    /// the report there.
    NcbiNullSecondRow,
}

/// NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnseq.hpp:139-185
/// `m_Score` and `m_ChainScore` sort and check nothing under these flags.
struct MixSeq {
    /// `m_BioseqHandle`: the sequence (`x_IdentifyAlnMixSeq`).
    key: SeqKey,
    /// `m_PositiveStrand`.
    positive_strand: bool,
    /// `m_Width`: 3 for a nucleotide sequence translated against proteins, else 1.
    width: u32,
    /// `m_StrandScore`: the lengths of the sequence's matches added so far, negative on the minus
    /// strand (`Mix::add`).
    strand_score: i64,
    /// `m_Frame`: the frame (`start % 3`) of a row of width 3 (`x_SetSeqFrame`), -1 until set.
    frame: i32,
    /// `m_DsCnt`: the Dense-segs that hold the sequence.
    ds_cnt: u32,
    seq_idx: i32,
    child_idx: i32,
    row_idx: i32,
    /// `m_DsIdx`: the Dense-seg whose subject row this row is (`fQuerySeqMergeOnly`).
    ds_idx: u32,
    extra_row: Option<SeqIx>,
    ref_by: Option<SeqIx>,
    /// `m_MatchList`: the matches of other refseqs that hold this sequence.
    match_list: Vec<usize>,
    /// `m_Starts`: the sequence's positions, each the start of a segment.
    starts: BTreeMap<u32, SegIx>,
    /// `m_Starts.current` (`Build`).
    current: It,
}

impl MixSeq {
    fn new(key: SeqKey) -> Self {
        // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnseq.cpp:353-375 (`CAlnMixSeq()`:
        // `m_DsCnt`, `m_StrandScore`, `m_ExtraRow`, `m_DsIdx` and `m_ChildIdx` start at 0,
        // `m_SeqIdx`, `m_RowIdx` and `m_Frame` at -1, `m_PositiveStrand` as true, `m_Width` at 1)
        Self {
            key,
            positive_strand: true,
            width: 1,
            strand_score: 0,
            frame: -1,
            ds_cnt: 0,
            seq_idx: -1,
            child_idx: 0,
            row_idx: -1,
            ds_idx: 0,
            extra_row: None,
            ref_by: None,
            match_list: Vec::new(),
            starts: BTreeMap::new(),
            current: None,
        }
    }
}

/// `CAlnMixMatch`: one segment of an input Dense-seg (`CAlnMixMatches::Add`).
#[derive(Clone, Copy)]
struct MixMatch {
    seq1: SeqIx,
    start1: u32,
    /// `m_AlnSeq2`, null for a segment of one row.
    seq2: Option<SeqIx>,
    start2: u32,
    len: u32,
    ds_idx: u32,
    score: i32,
    /// `m_StrandsDiffer`: the HSP's rows are on different strands (a match of two rows).
    strands_differ: bool,
}

/// `CAlnMixSegment`: `m_StartIts` keyed as `SSeqComp` orders it, by the row's
/// `(m_SeqIdx, m_ChildIdx)` (alnsegments.hpp:99-106), each with its row and the position its
/// iterator points to.
struct MixSegment {
    len: u32,
    start_its: BTreeMap<(i32, i32), (SeqIx, u32)>,
}

/// The state of one `CAlnMix` (its sequences, matches and merger).
struct Mix {
    seqs: Vec<MixSeq>,
    /// `CAlnMixSequences::m_Seqs`, in order.
    order: Vec<SeqIx>,
    /// `m_BioseqHandles`.
    by_key: BTreeMap<SeqKey, SeqIx>,
    matches: Vec<MixMatch>,
    segs: Vec<MixSegment>,
    /// `m_DsCnt`.
    ds_cnt: u32,
    /// `m_ExtraRows`: the rows that `x_Merge` made for the HSPs of a subject after its first.
    extra_rows: Vec<SeqIx>,
    /// `m_Rows` (`BuildRows`).
    rows: Vec<SeqIx>,
    /// `m_Segments` (`Build`).
    segments: Vec<SegIx>,
    /// The widths of the query's and of the subjects' rows (`super::Hsps::mix_widths`).
    query_width: u32,
    subject_width: u32,
    /// `m_AlnMixMatches->m_AddFlags`: the flags of the last `Add`, which the merge reads (all the
    /// Dense-segs of a frame group are added with one flag).
    add_flags: AddFlags,
    /// The subject that is the query's own Bioseq (`super::Hsps::set_query_bioseq`).
    query_bioseq: Option<u32>,
    /// `m_SingleRefseq` (`Mix::merge`; `m_IndependentDSs` is false for every section, there).
    single_refseq: bool,
    /// `x_Merge`'s `seg`, a `CRef` of `CAlnMixSegment` (alnmerger.cpp:161): the segment that the merge
    /// created last (`Mix::merge_segment`), which that `CRef` holds until `x_CreateDenseg` has
    /// returned (alnmerger.cpp:822-830, cited at `Mix::held`).
    merge_seg: Option<SegIx>,
}

/// The first key of `map` (`begin()`), `None` (`end()`) when it is empty.
fn begin(map: &BTreeMap<u32, SegIx>) -> It {
    map.keys().next().copied()
}

/// Moves `it` on by one (`it` is not `end()`).
fn next(map: &BTreeMap<u32, SegIx>, it: It) -> It {
    let key = it?;
    map.range(key + 1..).next().map(|(&key, _)| key)
}

/// The last key of `map` (`--end()`), `None` when it is empty.
fn last(map: &BTreeMap<u32, SegIx>) -> It {
    map.keys().next_back().copied()
}

/// `--it` (`it` is not `begin()`): from `end()` the last key.
fn prev(map: &BTreeMap<u32, SegIx>, it: It) -> It {
    match it {
        None => map.keys().next_back().copied(),
        Some(key) => map.range(..key).next_back().map(|(&key, _)| key),
    }
}

/// `lower_bound(key)`.
fn lower_bound(map: &BTreeMap<u32, SegIx>, key: u32) -> It {
    map.range(key..).next().map(|(&key, _)| key)
}

impl Mix {
    /// A new `CAlnMixSegment` of length `len`. The arena keeps every segment; NCBI frees one when
    /// its last `CRef` goes (`Mix::held`).
    fn new_segment(&mut self, len: u32) -> SegIx {
        self.segs.push(MixSegment {
            len,
            start_its: BTreeMap::new(),
        });
        self.segs.len() - 1
    }

    /// The new `CAlnMixSegment` that `x_Merge` gives to `seg` (alnmerger.cpp:492, 528, 589, 654, cited at
    /// `Mix::place`): a new segment, which `seg` now holds (`merge_seg`).
    fn merge_segment(&mut self, len: u32) -> SegIx {
        let seg = self.new_segment(len);
        self.merge_seg = Some(seg);
        seg
    }

    /// Whether a `CRef` still holds `seg` after `FillUnalignedRegions` replaced one in a row's
    /// map: the maps of the rows (`CAlnMixStarts`, a `map<TSeqPos, CRef<CAlnMixSegment> >`) of
    /// every sequence (`m_Seqs` and `m_ExtraRows` keep them all to the end of the merge), and
    /// `x_Merge`'s `seg`, which is alive until `x_CreateDenseg` has returned. Nothing else holds a
    /// segment then: `m_Segments` and `m_StartIts` hold raw pointers and map iterators
    /// (alnsegments.hpp:67,80,106,119), and the `CRef`s of `Build` (`gapped_segs`, `seg_stack`,
    /// `seg_set`) and `x_MinimizeGaps` (`seg`, `len_map`, `new_segs`) are locals that went with
    /// their functions (alnsegments.cpp:69,111-113,442-448).
    ///
    /// The replaced row's own map can hold it too, at another start: where a run of
    /// `x_MinimizeGaps` ends at a segment that shares a row with an earlier one of the run, each
    /// new segment is bound at the starts of both on that row (alnsegments.cpp:478-480,526, quoted
    /// at `minimize_gaps`), and a later fill can write at the earlier one, which frees nothing (a
    /// segment of letters: in `tests_x2g`, a subject row whose inserted letters lie on either side
    /// of a query letter that the subject lacks, an insertion next to a deletion). A segment of
    /// length 0 has no such second start where the fill replaces it: on the first row, whose
    /// segments `x_MinimizeGaps` never takes (alnsegments.cpp:254-256), it has one start (the
    /// merge's walk of a match's second row does not run for a length of 0, alnmerger.cpp:752),
    /// and the other rows of the first frame group hold segments of length 0 only, which
    /// `x_MinimizeGaps` makes into one segment per run that nothing of their rows follows, so no
    /// fill starts at one there. Only `seg` can keep a segment of length 0.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:161,822-830
    fn held(&self, seg: SegIx) -> bool {
        self.merge_seg == Some(seg)
            || self
                .seqs
                .iter()
                .any(|seq| seq.starts.values().any(|&held| held == seg))
    }

    // NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnsegments.hpp:108-115
    fn set_start_iterator(&mut self, seg: SegIx, seq: SeqIx, at: u32) {
        let key = (self.seqs[seq].seq_idx, self.seqs[seq].child_idx);
        self.segs[seg].start_its.insert(key, (seq, at));
    }

    /// `x_IdentifyAlnMixSeq`: the sequence of `key`, added to `m_Seqs` when new.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnseq.cpp:230-262
    fn identify(&mut self, key: SeqKey) -> SeqIx {
        if let Some(&seq) = self.by_key.get(&key) {
            return seq;
        }
        self.seqs.push(MixSeq::new(key));
        let seq = self.seqs.len() - 1;
        self.by_key.insert(key, seq);
        self.order.push(seq);
        seq
    }

    /// `CAlnMix::Add` of `ds` with `fPreserveRows`: `CAlnMixSequences::Add` and `CAlnMixMatches::Add` of
    /// one HSP (`subject`, `segments`; the query is row 0, the subject row 1).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnseq.cpp:105-208
    /// `strand1` and `strand2` are the rows' strands in the Dense-seg (`strands_exist` holds for
    /// BLAST's Dense-segs and `CreateDensegFromStdseg`'s): `Hsp::strands_differ` (a match of one row
    /// has `m_StrandsDiffer` false).
    ///
    /// `SeqKey::Query` is only row 0 and `SeqKey::Subject` only row 1. A subject that is the
    /// query's own Bioseq (`query_bioseq`) is the query's `CAlnMixSeq` in NCBI, so its row 1 takes
    /// the branch for `row_index` differing from `row`: an `m_AnotherRow` that is a sequence of its own in
    /// `m_Seqs` (as the subject's here) except for its `m_DsIdx`, which is the Dense-seg that made
    /// it, where another subject's is 0 until its first match (`x_SecondRowFits`): the first HSP
    /// added keeps the subject's first row, not the HSP of the longest segment.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnseq.cpp:165-202
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmatch.cpp:116-256
    ///
    /// With `fForceTranslation` (`flags`) the rows are not preserved (no `m_AnotherRow`), and each
    /// match of two rows first compares its strands with the strand scores of its sequences: a
    /// strand that differs from a sequence's score so far throws, which `DisplaySeqalign` catches
    /// and skips the Seq-align (showalign.cpp:2073-2076). The throw comes after the increment of `m_DsCnt` of the
    /// sequences and of the matches and after the matches of the Dense-seg's segments before it (a
    /// BLAST HSP's first segment holds both rows, so there are none): the Dense-seg still counts
    /// (`m_DsCnt`, the `m_DsIdx` of the later ones), its matches before the throw stay, and none of
    /// the throwing segment's scores are added. Returns whether the Dense-seg was added whole
    /// (false: left out).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmatch.cpp:145-210
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2073-2076
    /// So a TBLASTX subject's HSP of the other strand than the subject's first HSP in the same query
    /// frame is not shown in formats 1 to 4, where the description table and outfmt 0 show it (a
    /// deterministic result of NCBI BLAST+ 2.17.0, written as it is: PD-LOSAT-NCBI-DEFECTS rule 2,
    /// DECISIONS T3-6). A BLASTX HSP of a minus query frame (`fNegativeStrand`) never throws: its
    /// query is minus in every HSP of the group and a protein row has no strand.
    ///
    /// A TBLASTX HSP of the query with its own Bioseq (`query_bioseq`, added with
    /// `fForceTranslation`): without `fPreserveRows` the loop that runs while `row_index` differs
    /// from `row` (in the code cited above) is not entered, so both rows are the query's
    /// `CAlnMixSeq` (`x_IdentifyAlnMixSeq` finds the one Bioseq handle): its `m_DsCnt` counts both
    /// rows, and the strand check reads its one strand score for both rows before either is added
    /// to it.
    fn add(&mut self, hsp: super::Hsp<'_>, flags: AddFlags) -> Result<bool, MixError> {
        let (subject, segments) = (hsp.subject, hsp.segments);
        let preserve_rows = flags == AddFlags::PreserveRows;
        let query_bioseq = self.query_bioseq;
        let subject_key = if !preserve_rows && query_bioseq == Some(subject) {
            SeqKey::Query
        } else {
            SeqKey::Subject(subject)
        };
        self.ds_cnt += 1;
        self.add_flags = flags;
        let mut ds_seq = [0; 2];
        for (row, key) in [SeqKey::Query, subject_key].into_iter().enumerate() {
            let seq = self.identify(key);
            // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnseq.cpp:168-170
            // (The widths of `CreateDensegFromStdseg`'s Dense-seg of a TBLASTN or BLASTX HSP: 1 for
            // the protein row, 3 for the nucleotide row; of a TBLASTX HSP none, which
            // `x_ExtendDSWithWidths` sets to 3 for `fForceTranslation`, its lengths in residues;
            // `x_Merge` sets the same, alnmerger.cpp:138-145: 1 where `m_IsAA` is set, else 3. A Dense-seg of
            // two untranslated rows has none: 1.)
            // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmix.cpp:132-143
            // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmix.cpp:198-220
            self.seqs[seq].width = if row == 0 {
                self.query_width
            } else {
                self.subject_width
            };
            if preserve_rows && self.seqs[seq].row_idx == -1 {
                self.seqs[seq].row_idx = row as i32;
                if row == 1 && query_bioseq == Some(subject) {
                    // The query's `m_AnotherRow`.
                    self.seqs[seq].child_idx = 1;
                    self.seqs[seq].ds_idx = self.ds_cnt;
                }
            }
            self.seqs[seq].ds_cnt += 1;
            ds_seq[row] = seq;
        }
        let (minus1, minus2) = (hsp.strands.query_minus, hsp.strands.subject_minus);
        for segment in segments {
            let len = segment.len;
            match (segment.query, segment.subject) {
                (Some(start1), Some(start2)) => {
                    if !preserve_rows {
                        let mixed =
                            |score: i64, minus: bool| score > 0 && minus || score < 0 && !minus;
                        if mixed(self.seqs[ds_seq[0]].strand_score, minus1)
                            || mixed(self.seqs[ds_seq[1]].strand_score, minus2)
                        {
                            return Ok(false);
                        }
                    }
                    let signed = |minus: bool| {
                        if minus {
                            -i64::from(len)
                        } else {
                            i64::from(len)
                        }
                    };
                    self.seqs[ds_seq[0]].strand_score += signed(minus1);
                    self.seqs[ds_seq[1]].strand_score += signed(minus2);
                    self.matches.push(MixMatch {
                        seq1: ds_seq[0],
                        start1,
                        seq2: Some(ds_seq[1]),
                        start2,
                        len,
                        ds_idx: self.ds_cnt,
                        score: len as i32,
                        strands_differ: hsp.strands.differ(),
                    })
                }
                (Some(start1), None) => self.matches.push(MixMatch {
                    seq1: ds_seq[0],
                    start1,
                    seq2: None,
                    start2: 0,
                    len,
                    ds_idx: self.ds_cnt,
                    score: 0,
                    strands_differ: false,
                }),
                (None, Some(start1)) => self.matches.push(MixMatch {
                    seq1: ds_seq[1],
                    start1,
                    seq2: None,
                    start2: 0,
                    len,
                    ds_idx: self.ds_cnt,
                    score: 0,
                    strands_differ: false,
                }),
                (None, None) => {}
            }
        }
        Ok(true)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:835-860
    /// `x_SecondRowFits` under `fQuerySeqMergeOnly` without `fTruncateOverlaps`: whether row
    /// `seq2` is the row of Dense-seg `ds_idx` (`eSecondRowFitsOk`; a row of none takes it),
    /// else `eForceSeparateRow`. The rest of the function (strand, frame and overlap checks)
    /// is not reached.
    fn second_row_fits(&mut self, seq2: SeqIx, ds_idx: u32) -> bool {
        let seq = &mut self.seqs[seq2];
        if seq.ds_idx != 0 {
            seq.ds_idx == ds_idx
        } else {
            seq.ds_idx = ds_idx;
            true
        }
    }

    /// The part of `x_Merge` that places a match of `len` at `start1` on the starts of `seq1`:
    /// its segment map is split where the match starts and ends, and a segment made for each
    /// stretch the map does not hold yet. Returns `lo_start_i`, the start of the first segment
    /// of the match (where the second row's walk starts).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:466-702
    /// The `else` branches of the two splits (a row of the other strand than `seq1`, whose letters
    /// run the other way: the first part of the segment is its end):
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:550-557 (the look back)
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:614-621 (the overlap loop)
    /// Returns `(lo_start_i, hi_start_i)`: the first and the last segment of the match on `seq1`,
    /// where the second row's walk starts (`hi_start_i` when the strands differ). The `m_DsIdx` of
    /// a segment is read only by the truncation (`fTruncateOverlaps`), which these flags leave out.
    ///
    /// The first match placed on `seq1` (the query's first, the longest segment of all, ties in
    /// the order of `Add`) sets the `m_DsIdx` of its second row `seq2` to its Dense-seg. For an
    /// ordinary subject that is the value `x_SecondRowFits` has just set. For a subject that is the
    /// query's own Bioseq, whose `m_AnotherRow` `Add` preset to the HSP added first, it is the
    /// HSP of that match when another HSP of the subject holds it: the match itself still goes to
    /// an extra row (`x_SecondRowFits` returned `eForceSeparateRow` before), the HSP's other
    /// segments then fit the `m_AnotherRow`, and the HSP added first goes to the next extra row.
    /// When that HSP has no other segment the `m_AnotherRow` keeps no start, which NCBI does not
    /// survive (`MixError::NcbiAborts`, raised by `fill_unaligned_regions`). A first match
    /// without a second row has NCBI write through a null `seq2` (`MixError::NcbiNullSecondRow`):
    /// the matches of one row come after those of two (score 0, `SortByScore`) and every row of a
    /// later pass has the starts of its Dense-seg's matches with the query, so only a query read
    /// at the width 3 whose letters are a protein's makes one: its match of the query's letters
    /// alone opens a row of its own where no match of two rows has the frame of its start (each
    /// segment of a translated row has its HSP's frame).
    fn place(
        &mut self,
        seq1: SeqIx,
        seq2: Option<SeqIx>,
        matched: &MixMatch,
    ) -> Result<(It, It), MixError> {
        let (start1, len) = (matched.start1, matched.len);
        let width1 = self.seqs[seq1].width;
        let mut curr_len = len;
        let mut start_i: It = None;
        let mut lo_start_i: It;
        let mut hi_start_i: It = None;
        if self.seqs[seq1].starts.is_empty() {
            // alnmerger.cpp:483-489, cited above: the sign of the strand score with
            // `fForceTranslation`, else whether `m_MergeFlags` lacks `fNegativeStrand`, which no merge flag
            // here has (module doc).
            self.seqs[seq1].positive_strand = match self.add_flags {
                AddFlags::ForceTranslation => self.seqs[seq1].strand_score >= 0,
                AddFlags::PreserveRows => true,
            };
            let seg = self.merge_segment(len);
            self.seqs[seq1].starts.insert(start1, seg);
            self.set_start_iterator(seg, seq1, start1);
            lo_start_i = Some(start1);
            hi_start_i = Some(start1);
            // With `fQuerySeqMergeOnly` in `m_MergeFlags`, `seq2->m_DsIdx` is set from `match->m_DsIdx`
            // (alnmerger.cpp:499-501, cited above), also where `seq2` is null.
            let seq2 = seq2.ok_or(MixError::NcbiNullSecondRow)?;
            self.seqs[seq2].ds_idx = matched.ds_idx;
        } else {
            let starts = &self.seqs[seq1].starts;
            start_i = lower_bound(starts, start1);
            lo_start_i = start_i;
            if start_i.is_none_or(|key| start1 < key) && lo_start_i != begin(starts) {
                lo_start_i = prev(starts, lo_start_i);
            }
            if start_i != lo_start_i {
                let lo_key = lo_start_i.expect("a start before the match's");
                let prev_seg = self.seqs[seq1].starts[&lo_key];
                if lo_key + self.segs[prev_seg].len * width1 > start1 {
                    let len1 = (start1 - lo_key) / width1;
                    let len2 = self.segs[prev_seg].len - len1;
                    let seg = self.merge_segment(len2);
                    self.seqs[seq1].starts.insert(start1, seg);
                    self.split_rows(seq1, prev_seg, seg, len1, len2);
                    self.segs[prev_seg].len = len1;
                    if start_i != begin(&self.seqs[seq1].starts) {
                        start_i = prev(&self.seqs[seq1].starts, start_i);
                    }
                }
                if lo_start_i.is_some() {
                    lo_start_i = next(&self.seqs[seq1].starts, lo_start_i);
                }
            }
        }
        let mut start = start1;
        while hi_start_i.is_none() {
            match start_i {
                Some(key) if key == start => {
                    let prev_seg = self.seqs[seq1].starts[&key];
                    let prev_len = self.segs[prev_seg].len;
                    if prev_len > curr_len && curr_len == 0 {
                        // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:598-606
                        // A match of length 0 at a longer segment's start (`NcbiFreesSegment`):
                        // `curr_len` times `seq->m_Width` is 0, so each row's `CRef` of `prev_seg` is
                        // replaced. Only rows of the first row's strand reach this (a row of the
                        // other strand would keep a reference).
                        let one_strand = self.segs[prev_seg].start_its.values().all(|&(seq, _)| {
                            self.seqs[seq].positive_strand == self.seqs[seq1].positive_strand
                        });
                        return Err(if one_strand {
                            MixError::NcbiFreesSegment
                        } else {
                            MixError::Internal("a match of length 0 in a segment of two strands")
                        });
                    }
                    if prev_len > curr_len {
                        let len1 = prev_len - curr_len;
                        let seg = self.merge_segment(len1);
                        start += curr_len * width1;
                        self.segs[prev_seg].len = curr_len;
                        self.split_rows(seq1, prev_seg, seg, curr_len, len1);
                        hi_start_i = start_i;
                    } else if curr_len == prev_len {
                        hi_start_i = start_i;
                    } else {
                        start += prev_len * width1;
                        curr_len -= prev_len;
                        start_i = next(&self.seqs[seq1].starts, start_i);
                    }
                }
                _ => {
                    let seg = self.merge_segment(0);
                    self.seqs[seq1].starts.insert(start, seg);
                    let starts = &self.seqs[seq1].starts;
                    let mut tmp_start_i = start_i;
                    if tmp_start_i != begin(starts) {
                        tmp_start_i = prev(starts, tmp_start_i);
                    }
                    let tmp_key = tmp_start_i.expect("the new start");
                    let seg_len = match start_i {
                        Some(key) if start + curr_len * width1 > key => (key - start) / width1,
                        _ => {
                            // `hi_start_i` takes `start_i` and steps back by one unless it is
                            // the first start: the new start, which is before `start_i`.
                            hi_start_i = Some(tmp_key);
                            curr_len
                        }
                    };
                    self.set_start_iterator(seg, seq1, tmp_key);
                    self.segs[seg].len = seg_len;
                    start += seg_len * width1;
                    curr_len -= seg_len;
                    let starts = &self.seqs[seq1].starts;
                    if lo_start_i == start_i && lo_start_i != begin(starts) {
                        lo_start_i = prev(starts, lo_start_i);
                    }
                }
            }
        }
        Ok((lo_start_i, hi_start_i))
    }

    /// The rows of `prev_seg` cut after `len1` of its letters (its first `len1` stay in
    /// `prev_seg`, the `len2` after them go to `seg`): the `create rows info` loops of `place`. A
    /// row of `seq1`'s strand starts `seg` after its first `len1` letters; a row of the other
    /// strand, whose letters run the other way, starts `prev_seg` after its first `len2` and `seg`
    /// at the old start.
    fn split_rows(&mut self, seq1: SeqIx, prev_seg: SegIx, seg: SegIx, len1: u32, len2: u32) {
        let plus1 = self.seqs[seq1].positive_strand;
        let rows: Vec<(SeqIx, u32)> = self.segs[prev_seg].start_its.values().copied().collect();
        for (seq, tmp_start) in rows {
            let width = self.seqs[seq].width;
            if self.seqs[seq].positive_strand == plus1 {
                let at = tmp_start + len1 * width;
                self.seqs[seq].starts.insert(at, seg);
                let after = next(&self.seqs[seq].starts, Some(tmp_start)).expect("the new start");
                self.set_start_iterator(seg, seq, after);
            } else {
                let at = tmp_start + len2 * width;
                self.seqs[seq].starts.insert(at, prev_seg);
                self.seqs[seq].starts.insert(tmp_start, seg);
                self.set_start_iterator(seg, seq, tmp_start);
                let after = next(&self.seqs[seq].starts, Some(tmp_start)).expect("the new start");
                self.set_start_iterator(prev_seg, seq, after);
            }
        }
    }

    /// `x_Merge` for these flags, then `BuildRows`, `Build` with false, true and false,
    /// `FillUnalignedRegions` and `x_CreateDenseg` (alnmerger.cpp:822-830).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:101-133
    /// The query is the first sequence and in every Dense-seg: it is the single refseq, unless a
    /// TBLASTX HSP of the query's own Bioseq counted it twice (`Mix::add`); then no sequence is
    /// the single refseq and the query, first in `m_Seqs`, stays the refseq of the first pass. Its
    /// `m_DsCnt` is above 1 whenever the mix has more than one Dense-seg (and also with one
    /// Dense-seg of its own Bioseq), so the Dense-segs are never independent.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:160-277
    /// The query's pass takes every match in score order; a match of a subject's letters alone
    /// (an insertion) goes to that subject's `m_MatchList`, which the next passes take, one
    /// subject at a time in the order of `m_Seqs`. `m_AlnSeq2` is never the refseq: the query is
    /// always `m_AlnSeq1`, and a subject's matches in its list have no second row. A match of the
    /// query with its own Bioseq (`m_AlnSeq1` is `m_AlnSeq2`) is the query's in the first pass, where
    /// the branch for `m_AlnSeq1` being `refseq` takes it; its second increment of `match_list_i` belongs
    /// to the later passes, which a mix without a single refseq never reaches (TBLASTX's HSPs
    /// have no gaps, so no match of one row).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:470-481
    /// (In the `!starts.size()` branch of a later pass: not reached, as above; an internal error
    /// if it were.)
    fn merge(&mut self) -> Result<(), MixError> {
        // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmix.cpp:279-286
        // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmatch.cpp:62-87
        // (A stable sort of indices: `sort_by_key` keeps equal scores in input order.)
        let mut order: Vec<usize> = (0..self.matches.len()).collect();
        order.sort_by_key(|&index| std::cmp::Reverse(self.matches[index].score));
        // `Find the refseq (if such exists)`, cited above.
        // (`m_IndependentDSs`, false here, is left out.)
        self.single_refseq = false;
        let mut first = true;
        for index in 0..self.order.len() {
            if self.seqs[self.order[index]].ds_cnt == self.ds_cnt {
                self.single_refseq = true;
                if !first {
                    let refseq = self.order.remove(index);
                    self.order.insert(0, refseq);
                }
                break;
            }
            first = false;
        }
        if self.order.first().map(|&seq| self.seqs[seq].key) != Some(SeqKey::Query) {
            // The query is in every Dense-seg, so no other sequence is in all of them more
            // often than it.
            return Err(MixError::Internal(
                "the query is not the refseq of the merge",
            ));
        }
        for (seq_idx, &seq) in self.order.iter().enumerate() {
            self.seqs[seq].seq_idx = seq_idx as i32;
        }
        let mut refseq = self.order[0];
        let mut first_refseq = true;
        let mut match_i = 0usize;
        let mut match_list: Vec<usize> = Vec::new();
        let mut match_list_i = 0usize;
        loop {
            let at_end = if first_refseq {
                match_i == order.len() && !order.is_empty()
            } else {
                match_list.is_empty() || match_list_i == match_list.len()
            };
            if at_end {
                if !first_refseq {
                    // The matches of the list were erased from it as they were used
                    // (`erase` of `match_list_iter1` from `seq1->m_MatchList`).
                    self.seqs[refseq].match_list.clear();
                }
                self.seqs[refseq].ref_by = None;
                if let Some(&seq) = self.order.iter().find(|&&seq| {
                    !self.seqs[seq].match_list.is_empty() && self.seqs[seq].ref_by == Some(refseq)
                }) {
                    refseq = seq;
                }
                if self.seqs[refseq].ref_by.is_none() {
                    if let Some(&seq) = self
                        .order
                        .iter()
                        .find(|&&seq| !self.seqs[seq].match_list.is_empty())
                    {
                        refseq = seq;
                    }
                }
                if self.seqs[refseq].match_list.is_empty() {
                    break;
                }
                if !self.single_refseq {
                    // A match of one row in a mix of the query's own Bioseq added with
                    // `fForceTranslation`: TBLASTX's HSPs have none (the doc above).
                    return Err(MixError::Internal(
                        "a later pass of a merge without a single refseq",
                    ));
                }
                first_refseq = false;
                match_list = self.seqs[refseq].match_list.clone();
                match_list_i = 0;
                continue;
            }
            let index = if first_refseq {
                match_i += 1;
                order[match_i - 1]
            } else {
                match_list_i += 1;
                match_list[match_list_i - 1]
            };
            let matched = self.matches[index];
            if matched.seq1 != refseq {
                // A subject's letters alone: a candidate for that subject's pass.
                self.seqs[matched.seq1].match_list.push(index);
                continue;
            }
            self.merge_match(matched, first_refseq)?;
        }
        Ok(())
    }

    /// One match of the refseq (`seq1`), from `// save the match info into the segments map`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:349-380
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:379-416
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:687-818
    /// A width-3 refseq (a TBLASTN subject's row in a later pass, its insertions) first goes to the
    /// row of its frame (alnmerger.cpp:335-349, `x_SetSeqFrame` of `match` and `match->m_AlnSeq1`), from
    /// which the proper row is looked for. Every row of a TBLASTN HSP has its frame (`start % 3`
    /// of each of its segments is the same), so that row is at or before the HSP's own.
    ///
    /// The row of an ordinary subject's HSP holds that HSP's aligned segments alone, which never
    /// share a subject position: there the walk never meets a position of the row.
    ///
    /// A query read at the width 3 whose letters are a protein's (module doc) has a row per frame,
    /// and the first match placed on each row of the query sets the `m_DsIdx` of its second row
    /// (`place`), which is the subject's first row whatever the match's Dense-seg (`x_SecondRowFits`
    /// has not moved it yet). So a later frame's first match gives the subject's first row to the
    /// HSP of that match: the HSPs of the subject after it can then fit that row too (their
    /// segments meet the earlier HSP's there, `X_MERGE_START_ITERATOR`, or join them), and a match
    /// of a later pass that looks for the row of the first HSP finds none
    /// (`X_MERGE_PROPER_ROW`, `MixError::NcbiThrows`).
    ///
    /// The query's own Bioseq added with `fForceTranslation` (TBLASTX, `Mix::add`): a match of
    /// such an HSP has the query's sequence on both sides. `x_SetSeqFrame` of the second row
    /// compares the row with the match's first row (whether `seq` is `match->m_AlnSeq1`, both the query's
    /// row of the group's frame) and takes the frame of the query's start, so the second row is
    /// the query's own row; `x_SecondRowFits` gives it the first such HSP (its `m_DsIdx`) and
    /// every later one a new extra row of the query. The first one's walk runs on the row of the
    /// query's starts: at a start that exists it joins the two segments (`this position already
    /// exists`, which is how the record's HSP with itself in the same frame folds into the query's
    /// row), and a start that does not exist throws, because the query's segment there already has
    /// the row's start iterator (`X_MERGE_START_ITERATOR`, `MixError::NcbiThrows`). A join of two
    /// segments of the query's row (a first such HSP whose subject positions are starts of the
    /// query's other matches, as in a tandem repeat) can leave the row's segments in a loop, which
    /// `build` then meets (`BUILD_INFINITE_LOOP`).
    fn merge_match(&mut self, matched: MixMatch, first_refseq: bool) -> Result<(), MixError> {
        let mut seq1 = matched.seq1;
        if self.seqs[seq1].width == 3 {
            seq1 = self.set_seq_frame(seq1, matched.start1, matched.ds_idx);
        }
        if !first_refseq {
            loop {
                if self.seqs[seq1].ds_idx == matched.ds_idx {
                    break;
                }
                match self.seqs[seq1].extra_row {
                    Some(row) => seq1 = row,
                    None => {
                        return Err(MixError::NcbiThrows(NcbiException::new(
                            X_MERGE_PROPER_ROW,
                            X_MERGE_PROPER_ROW_MESSAGE,
                        )))
                    }
                }
            }
        }
        let mut seq2 = matched.seq2;
        let mut fits = true;
        if let Some(mut row) = seq2 {
            self.seqs[row].ref_by = Some(self.order[0]);
            if self.seqs[row].width == 3 {
                // `frame` is `m_Start1` where `seq` is `match->m_AlnSeq1`, else `m_Start2`, modulo 3 (alnmerger.cpp:
                // 1272-1276, cited at `set_seq_frame`): the query's start where the second row
                // is the first row's sequence (its own Bioseq, `Mix::add`).
                let start = if row == seq1 {
                    matched.start1
                } else {
                    matched.start2
                };
                row = self.set_seq_frame(row, start, matched.ds_idx);
            }
            fits = self.second_row_fits(row, matched.ds_idx);
            seq2 = Some(row);
        }
        let (lo_start_i, hi_start_i) = self.place(seq1, seq2, &matched)?;
        let Some(mut row) = seq2 else {
            return Ok(());
        };
        while !fits {
            match self.seqs[row].extra_row {
                None => {
                    let mut extra = MixSeq::new(self.seqs[row].key);
                    extra.width = self.seqs[row].width;
                    extra.frame = (matched.start2 % 3) as i32;
                    extra.seq_idx = self.seqs[row].seq_idx;
                    extra.child_idx = self.seqs[row].child_idx + 1;
                    extra.ds_idx = matched.ds_idx;
                    self.seqs.push(extra);
                    let extra = self.seqs.len() - 1;
                    self.extra_rows.push(extra);
                    self.seqs[row].extra_row = Some(extra);
                    row = extra;
                    break;
                }
                Some(extra) => {
                    row = extra;
                    fits = self.second_row_fits(row, matched.ds_idx);
                }
            }
        }
        if self.seqs[row].starts.is_empty() {
            self.seqs[row].positive_strand = if self.seqs[seq1].positive_strand {
                !matched.strands_differ
            } else {
                matched.strands_differ
            };
        }
        let width2 = self.seqs[row].width;
        let mut start = matched.start2;
        let mut start2_i = lower_bound(&self.seqs[row].starts, matched.start2);
        let mut start_i = if matched.strands_differ {
            hi_start_i
        } else {
            lo_start_i
        };
        while start < matched.start2 + matched.len * width2 {
            // `start_i` stays on the first row's segments of the match (their lengths add up to
            // the match's); `end()` here would be NCBI's dereference of `end()`.
            let seg = self.seqs[seq1].starts[&start_i.ok_or(MixError::Internal(
                "the second row's walk passed the first row's segments",
            ))?];
            if start2_i == Some(start) {
                // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:753-771
                // (Only a row of the query's own Bioseq has the start: the query's row, whose
                // walk meets the starts of the query's matches. The segment left without rows
                // is in no list.)
                let existing = self.seqs[row].starts[&start];
                if existing != seg {
                    let rows: Vec<(SeqIx, u32)> =
                        self.segs[existing].start_its.values().copied().collect();
                    for (tmp_seq, at) in rows {
                        self.seqs[tmp_seq].starts.insert(at, seg);
                        self.set_start_iterator(seg, tmp_seq, at);
                    }
                }
            } else {
                self.seqs[row].starts.insert(start, seg);
                start2_i = Some(start);
                let key = (self.seqs[row].seq_idx, self.seqs[row].child_idx);
                if self.segs[seg].start_its.contains_key(&key) {
                    return Err(MixError::NcbiThrows(NcbiException::new(
                        X_MERGE_START_ITERATOR,
                        X_MERGE_START_ITERATOR_MESSAGE,
                    )));
                }
                self.set_start_iterator(seg, row, start);
            }
            start += self.segs[seg].len * width2;
            start2_i = next(&self.seqs[row].starts, start2_i);
            let starts = &self.seqs[seq1].starts;
            if matched.strands_differ {
                if start_i != begin(starts) {
                    start_i = prev(starts, start_i);
                }
            } else {
                start_i = next(starts, start_i);
            }
        }
        Ok(())
    }

    /// `x_SetSeqFrame`: the row of `seq`'s chain (`m_ExtraRow`) in the frame `start % 3` of a
    /// match at `start` of Dense-seg `ds_idx`: `seq` itself while it has no start (it takes that
    /// frame), else the first row of the chain in that frame, else a new row at the chain's end.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:1268-1302
    /// (`match` is the merge's copy of the match, alnmerger.cpp:417-422: the row it moves to is
    /// that copy's `m_AlnSeq2`.)
    fn set_seq_frame(&mut self, mut seq: SeqIx, start: u32, ds_idx: u32) -> SeqIx {
        let frame = (start % 3) as i32;
        if self.seqs[seq].starts.is_empty() {
            self.seqs[seq].frame = frame;
            return seq;
        }
        while self.seqs[seq].frame != frame {
            match self.seqs[seq].extra_row {
                None => {
                    let mut new_seq = MixSeq::new(self.seqs[seq].key);
                    new_seq.positive_strand = self.seqs[seq].positive_strand;
                    new_seq.width = self.seqs[seq].width;
                    new_seq.frame = frame;
                    new_seq.seq_idx = self.seqs[seq].seq_idx;
                    new_seq.child_idx = self.seqs[seq].child_idx + 1;
                    new_seq.ds_idx = ds_idx;
                    self.seqs.push(new_seq);
                    let new_seq = self.seqs.len() - 1;
                    self.extra_rows.push(new_seq);
                    self.seqs[seq].extra_row = Some(new_seq);
                    return new_seq;
                }
                Some(extra) => seq = extra,
            }
        }
        seq
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnseq.cpp:265-283
    fn build_rows(&mut self) {
        self.rows.clear();
        let mut count = 0;
        for index in 0..self.order.len() {
            let mut seq = self.order[index];
            if !self.seqs[seq].starts.is_empty() {
                self.rows.push(seq);
                self.seqs[seq].row_idx = count;
                count += 1;
                while let Some(extra) = self.seqs[seq].extra_row {
                    seq = extra;
                    self.seqs[seq].row_idx = count;
                    count += 1;
                    self.rows.push(seq);
                }
            }
        }
    }

    /// `Build` with `gap_join` false, `min_gap` true and `remove_leading_and_trailing_gaps` false.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnsegments.cpp:58-271
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnseq.cpp:286-294
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnseq.cpp:286-296
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnsegments.cpp:216-224
    /// A minus-strand row is read from its last start down. `InitExtraRowsStartIts` sets the same
    /// of every extra row (`m_ExtraRows`), also of one that `BuildRows` left out with the sequence
    /// of its first row (`MixError::NcbiAborts`).
    fn build(&mut self) -> Result<(), MixError> {
        for index in 0..self.rows.len() {
            let row = self.rows[index];
            self.seqs[row].current = self.first_current(row);
        }
        for index in 0..self.extra_rows.len() {
            let row = self.extra_rows[index];
            self.seqs[row].current = self.first_current(row);
        }
        let mut gapped: Vec<SegIx> = Vec::new();
        let mut refseq_it = 0;
        let mut orig_refseq = true;
        loop {
            let mut refseq = None;
            while refseq_it < self.rows.len() {
                let row = self.rows[refseq_it];
                refseq_it += 1;
                if self.seqs[row].current.is_some() {
                    refseq = Some(row);
                    break;
                }
            }
            let Some(refseq) = refseq else {
                if !gapped.is_empty() {
                    self.minimize_gaps(&mut gapped);
                    self.segments.append(&mut gapped);
                }
                break;
            };
            while let Some(current) = self.seqs[refseq].current {
                let first = self.seqs[refseq].starts[&current];
                let mut stack = vec![first];
                let mut on_stack = BTreeSet::from([first]);
                while let Some(&top) = stack.last() {
                    let mut pop_seg = true;
                    for &(row, at) in self.segs[top].start_its.values() {
                        let current = self.seqs[row].current;
                        if current != Some(at) {
                            let left = self.seqs[row].starts[&current
                                .ok_or(MixError::Internal("a row's segments ended early"))?];
                            stack.push(left);
                            if !on_stack.insert(left) {
                                // NCBI's throw (cited above), which only the query's own Bioseq
                                // with `fForceTranslation` (`merge_match`) and a query read at
                                // the width 3 whose letters are a protein's (module doc) reach.
                                return Err(MixError::NcbiThrows(NcbiException::new(
                                    BUILD_INFINITE_LOOP,
                                    BUILD_INFINITE_LOOP_MESSAGE,
                                )));
                            }
                            pop_seg = false;
                            break;
                        }
                    }
                    if !pop_seg {
                        continue;
                    }
                    let rows: Vec<SeqIx> = self.segs[top]
                        .start_its
                        .values()
                        .map(|&(row, _)| row)
                        .collect();
                    for row in rows {
                        let (current, starts) = (self.seqs[row].current, &self.seqs[row].starts);
                        self.seqs[row].current = if self.seqs[row].positive_strand {
                            next(starts, current)
                        } else if current == begin(starts) {
                            None
                        } else {
                            prev(starts, current)
                        };
                    }
                    stack.pop();
                    if !stack.is_empty() {
                        gapped.push(top);
                    } else {
                        if !gapped.is_empty() {
                            self.minimize_gaps(&mut gapped);
                            if orig_refseq {
                                self.segments.append(&mut gapped);
                            }
                        }
                        if orig_refseq {
                            self.segments.push(top);
                        } else {
                            gapped.push(top);
                        }
                    }
                }
            }
            orig_refseq = false;
        }
        Ok(())
    }

    /// The first start of `row` that `Build` reads (`InitRowsStartIts`): its first, or the last
    /// of a minus-strand row; `None` (`end()`) for a row without starts.
    fn first_current(&self, row: SeqIx) -> It {
        let starts = &self.seqs[row].starts;
        if self.seqs[row].positive_strand {
            begin(starts)
        } else {
            last(starts)
        }
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnsegments.cpp:438-543
    /// The run of segments that ends at a segment sharing a row with an earlier one includes that
    /// segment (`seg_i_end` is advanced before the run is cut), as NCBI's does.
    fn minimize_gaps(&mut self, gapped: &mut Vec<SegIx>) {
        let mut new_segs = Vec::new();
        let mut seg_i_begin = 0;
        let mut seg_i_end = 0;
        let mut len_map: BTreeMap<u32, SegIx> = BTreeMap::new();
        let mut lens: BTreeSet<u32> = BTreeSet::new();
        while seg_i_end < gapped.len() {
            lens.insert(self.segs[gapped[seg_i_end]].len);
            let possible = gapped[seg_i_begin..seg_i_end].iter().all(|&seg1| {
                self.segs[gapped[seg_i_end]]
                    .start_its
                    .keys()
                    .all(|key| !self.segs[seg1].start_its.contains_key(key))
            });
            seg_i_end += 1;
            if !possible || seg_i_end == gapped.len() {
                let mut len_so_far = 0;
                for &len in &lens {
                    let seg = self.new_segment(len - len_so_far);
                    len_so_far += self.segs[seg].len;
                    len_map.insert(len, seg);
                }
                for &orig in &gapped[seg_i_begin..seg_i_end] {
                    let orig_len = self.segs[orig].len;
                    let rows: Vec<(SeqIx, u32)> =
                        self.segs[orig].start_its.values().copied().collect();
                    for (seq, orig_start) in rows {
                        let (plus, width) = (self.seqs[seq].positive_strand, self.seqs[seq].width);
                        let mut len_so_far = 0;
                        for (_, &seg) in len_map.range(..=orig_len) {
                            let offset = if plus {
                                len_so_far
                            } else {
                                orig_len - len_so_far - self.segs[seg].len
                            };
                            let this_start = orig_start + offset * width;
                            self.seqs[seq].starts.insert(this_start, seg);
                            self.set_start_iterator(seg, seq, this_start);
                            len_so_far += self.segs[seg].len;
                        }
                    }
                }
                new_segs.extend(len_map.values().copied());
                len_map.clear();
                lens.clear();
                seg_i_begin = seg_i_end;
            }
        }
        *gapped = new_segs;
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnsegments.cpp:286-338
    /// The new segment's start iterator (`start_i` moved back by one, or on by one on the minus strand) is its own
    /// start: no start of the row lies between it and the segment being read. The fill of a row's
    /// stretch goes before the segment being read, whose own rows are then
    /// still read from it (`seg_i` points to it again after the insertion).
    ///
    /// After a segment of length 0 of a plus-strand row the fill starts at that segment's own
    /// start (`prev_start + prev_len * width`), so the write into the row's map replaces the row's
    /// `CRef` of it. Where that was its last (`Mix::held`: it is not `x_Merge`'s `seg`, the only
    /// other `CRef` that a segment of length 0 has there), NCBI frees it there, and
    /// `x_CreateDenseg` reads it through `m_Segments` (alnmerger.cpp:1196-1210, cited at
    /// `create_denseg`): a read of freed memory, `MixError::NcbiFreesSegment`, at the first such
    /// replacement (official NCBI BLAST+ 2.17.0 crashed with SIGSEGV in our runs; valgrind names
    /// the block "free'd by `FillUnalignedRegions`", the review 3 of line RG T3c, case r11_001). A
    /// fill that writes at a start where the row's map still holds the replaced segment at another
    /// one frees nothing (`Mix::held`). A segment that
    /// `m_Segments` does not list is read by nothing after its free. Only BLASTX's HSPs of length
    /// 0 make segments of length 0, on the query's rows of the first frame group, which are plus
    /// (`fPreserveRows`, `place`): no minus-strand row holds one (there the fill starts at the end
    /// of the segment being read, that segment's own start only for a length of 0). The only
    /// other stop of this function, `MixError::NcbiAborts`, needs the query's own Bioseq, which
    /// BLASTX's protein subjects never are.
    ///
    /// A row that `BuildRows` left out keeps the `m_RowIdx` of -1 of `CAlnMixSeq()`: the extra
    /// rows of a sequence without starts, the query's own Bioseq whose `m_AnotherRow` kept none
    /// (`place`). NCBI's `starts[rowidx]` and `lens[rowidx]` then write before their vectors'
    /// buffers (the stores of the previous start and the previous length, a length of at least
    /// 1), and glibc stops NCBI when it frees them (`double free or corruption (out)`, exit
    /// 134): no section to write, `MixError::NcbiAborts` (line RG T3a, PD-LOSAT-NCBI-DEFECTS
    /// rule 3).
    fn fill_unaligned_regions(&mut self) -> Result<(), MixError> {
        let mut starts = vec![-1i64; self.rows.len()];
        let mut lens = vec![0u32; self.rows.len()];
        let mut seg_i = 0;
        while seg_i < self.segments.len() {
            let current = self.segments[seg_i];
            let len = self.segs[current].len;
            let rows: Vec<(SeqIx, u32)> = self.segs[current].start_its.values().copied().collect();
            for (row, start) in rows {
                let Ok(rowidx) = usize::try_from(self.seqs[row].row_idx) else {
                    return Err(MixError::NcbiAborts);
                };
                let (prev_start, prev_len) = (starts[rowidx], lens[rowidx]);
                let (plus, width) = (self.seqs[row].positive_strand, self.seqs[row].width);
                let prev_start_plus_len = (prev_start as u32).wrapping_add(prev_len * width);
                let start_plus_len = start + len * width;
                if prev_start >= 0
                    && (plus && prev_start_plus_len < start
                        || !plus && start_plus_len < prev_start as u32)
                {
                    let (new_start, seg_len) = if plus {
                        (prev_start_plus_len, (start - prev_start_plus_len) / width)
                    } else {
                        (start_plus_len, (prev_start as u32 - start_plus_len) / width)
                    };
                    let seg = self.new_segment(seg_len);
                    if let Some(old) = self.seqs[row].starts.insert(new_start, seg) {
                        // The map's `CRef` of `old` goes (`MixError::NcbiFreesSegment`): its last
                        // one frees it, which `x_CreateDenseg` reads where `m_Segments` lists it.
                        if !self.held(old) && self.segments.contains(&old) {
                            return Err(MixError::NcbiFreesSegment);
                        }
                    }
                    self.set_start_iterator(seg, row, new_start);
                    self.segments.insert(seg_i, seg);
                    seg_i += 1;
                }
                starts[rowidx] = i64::from(start);
                lens[rowidx] = len;
            }
            seg_i += 1;
        }
        Ok(())
    }

    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:1156-1252
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:1215-1240
    /// A row's strand is the same in every segment (kept once per row). Every row of a segment
    /// has a row index here (`fill_unaligned_regions` stops at one without).
    fn create_denseg(&self) -> Result<Denseg, MixError> {
        let numrows = self.rows.len();
        let mut starts = vec![-1i32; numrows * self.segments.len()];
        let mut lens = Vec::with_capacity(self.segments.len());
        let mut row_empty = vec![true; numrows];
        for (numseg, &seg) in self.segments.iter().enumerate() {
            lens.push(self.segs[seg].len);
            for &(row, start) in self.segs[seg].start_its.values() {
                let row_idx = usize::try_from(self.seqs[row].row_idx)
                    .map_err(|_| MixError::Internal("a segment's row has no row index"))?;
                starts[numseg * numrows + row_idx] = start as i32;
                row_empty[row_idx] = false;
            }
        }
        let keep: Vec<usize> = (0..numrows).filter(|&row| !row_empty[row]).collect();
        // The rows of the query's own Bioseq after the first (`fForceTranslation`, `Mix::add`)
        // hold the record's subject positions: shown as that subject (its letters, label and
        // genetic code: `CAlnVec::SetGenCode` gives `m_SlaveGeneticCode` to every row but row 0,
        // showalign.cpp:2100-2101, cited at `super::write_section`). The query's first row is
        // row 0: it is the first sequence and has starts.
        let ids: Vec<SeqKey> = keep
            .iter()
            .enumerate()
            .map(
                |(index, &row)| match (self.seqs[self.rows[row]].key, self.query_bioseq) {
                    (SeqKey::Query, Some(subject)) if index > 0 => SeqKey::Subject(subject),
                    (key, _) => key,
                },
            )
            .collect();
        let plus: Vec<bool> = keep
            .iter()
            .map(|&row| self.seqs[self.rows[row]].positive_strand)
            .collect();
        // (A query whose Bioseq is not a protein with nucleotide subjects, TBLASTN's of the
        // molecule `other`, has no widths in NCBI's Dense-seg, which holds nucleotides only: its
        // section reads nothing of the merged Dense-seg, `super::write_unreadable_section`.)
        let widths: Vec<u32> = if self.query_width == 3 || self.subject_width == 3 {
            keep.iter()
                .map(|&row| self.seqs[self.rows[row]].width)
                .collect()
        } else {
            Vec::new()
        };
        let starts = if keep.len() == numrows {
            starts
        } else {
            (0..lens.len())
                .flat_map(|seg| keep.iter().map(move |&row| (seg, row)))
                .map(|(seg, row)| starts[seg * numrows + row])
                .collect()
        };
        Ok(Denseg {
            ids,
            lens,
            starts,
            plus,
            widths,
        })
    }
}

/// The merged Dense-seg of the HSPs `group` of a query's section (one frame group, in order):
/// `CAlnMix::Add` of each with the group's flag (`Hsps::add_flags`), then `Merge` with `fMinGap`,
/// `fQuerySeqMergeOnly` and `fFillUnalignedRegions`, and `GetDenseg`; `None` when every `Add` threw
/// (`hasAln` is false: no merge, and `GetDenseg` throws for the group).
pub(crate) fn merge(hsps: &Hsps, group: &[usize]) -> Result<Option<Denseg>, MixError> {
    // The rows' widths as `x_Merge` sets them (`Hsps::mix_widths`; the same as the Dense-segs'
    // but for a query whose Bioseq is not a protein).
    let (query_width, subject_width) = hsps.mix_widths();
    let mut mix = Mix {
        seqs: Vec::new(),
        order: Vec::new(),
        by_key: BTreeMap::new(),
        matches: Vec::new(),
        segs: Vec::new(),
        ds_cnt: 0,
        extra_rows: Vec::new(),
        rows: Vec::new(),
        segments: Vec::new(),
        query_width,
        subject_width,
        add_flags: AddFlags::PreserveRows,
        query_bioseq: hsps.query_bioseq,
        single_refseq: true,
        merge_seg: None,
    };
    let mut has_aln = false;
    for &index in group {
        let hsp = hsps.hsp(index);
        has_aln |= mix.add(hsp, hsps.add_flags(hsp.strands))?;
    }
    if !has_aln {
        return Ok(None);
    }
    mix.merge()?;
    mix.build_rows();
    mix.build()?;
    mix.fill_unaligned_regions()?;
    mix.create_denseg().map(Some)
}
