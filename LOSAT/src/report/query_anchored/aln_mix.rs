//! NCBI's `CAlnMix` (`objtools/alnmgr`) for the flags of the query-anchored reports: `Add` with
//! `fPreserveRows` or `fForceTranslation` (`super::AddFlags`), `Merge` with `fMinGap |
//! fQuerySeqMergeOnly | fFillUnalignedRegions` (showalign.cpp:2058-2083, quoted at
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
//! which `DisplaySeqalign` does not catch) is `MixError::NcbiThrows`.
//!
//! Strands and widths (lines RG T3b and T3c). A row's strand is set where it gets its first start.
//! The query's (the refseq of the first pass) is plus with `fPreserveRows`
//! (`!(m_MergeFlags & fNegativeStrand)`: no merge flag has it) and, with `fForceTranslation`, the
//! sign of its strand score (the sum of its matches' lengths, negative on the minus strand), so a
//! BLASTX or TBLASTX query of a minus frame has a minus-strand row. A subject's row takes the
//! query's strand, or the other where the HSP's strands differ (`m_StrandsDiffer`), so the merge's
//! branches for rows of the other strand run. A nucleotide sequence translated against proteins
//! (BLASTX's query, TBLASTN's subjects, TBLASTX's both: `m_ContainsAA && m_ContainsNA`, or
//! `fForceTranslation`) has rows of width 3 in their frames (`x_SetSeqFrame`). With
//! `fForceTranslation` an HSP whose strand on a sequence differs from the sequence's strand score
//! so far is left out (`CAlnMixMatches::Add` throws, `Mix::add`).

use std::collections::{BTreeMap, BTreeSet};

use super::{AddFlags, Denseg, Hsps, SeqKey};

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
    /// NCBI's merge throws `CAlnException` (`eMergeFailure`, its message here), which
    /// `DisplaySeqalign` does not catch (only `Add` is in its `try`, showalign.cpp:2060-2083, quoted
    /// at `super::write_section`): the run stops with `Error: NCBI C++ Exception` (exit 255; probe
    /// `tx.t3c.d.qid.f1` of line RG T3c). Only the query's own Bioseq added with
    /// `fForceTranslation` gets there (`Mix::merge_match`, `Mix::build`); the caller refuses the
    /// section, as for `NcbiAborts`.
    NcbiThrows(&'static str),
}

/// NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnseq.hpp:139-185
/// ```c++
/// class NCBI_XALNMGR_EXPORT CAlnMixSeq : public CObject
/// {
///     ...
///     int                   m_DsCnt;
///     const CBioseq_Handle* m_BioseqHandle;
///     ...
///     bool                  m_PositiveStrand;
///     CAlnMixSeq *          m_RefBy;
///     CAlnMixSeq *          m_ExtraRow;
///     int                   m_ExtraRowIdx;
///     CAlnMixSeq *          m_AnotherRow;
///     int                   m_DsIdx; // used by the truncate algorithm
///     int                   m_SeqIdx;
///     int                   m_ChildIdx;
///     int                   m_RowIdx;
///     TMatchList            m_MatchList;
///     ...
/// ```
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
        // `m_DsCnt(0)`, `m_StrandScore(0)`, `m_ExtraRow(0)`, `m_DsIdx(0)`, `m_SeqIdx(-1)`,
        // `m_ChildIdx(0)`, `m_RowIdx(-1)`, `m_PositiveStrand(true)`, `m_Width(1)`, `m_Frame(-1)`)
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
    /// The widths of the query's and of the subjects' rows (`super::Hsps`).
    query_width: u32,
    subject_width: u32,
    /// `m_AlnMixMatches->m_AddFlags`: the flags of the last `Add`, which the merge reads (all the
    /// Dense-segs of a frame group are added with one flag).
    add_flags: AddFlags,
    /// The subject that is the query's own Bioseq (`super::Hsps::set_query_bioseq`).
    query_bioseq: Option<u32>,
    /// `m_SingleRefseq` (`Mix::merge`; `m_IndependentDSs` is false for every section, there).
    single_refseq: bool,
}

/// The first key of `map` (`begin()`), `None` (`end()`) when it is empty.
fn begin(map: &BTreeMap<u32, SegIx>) -> It {
    map.keys().next().copied()
}

/// `++it` (`it` is not `end()`).
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
    fn new_segment(&mut self, len: u32) -> SegIx {
        self.segs.push(MixSegment {
            len,
            start_its: BTreeMap::new(),
        });
        self.segs.len() - 1
    }

    // NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnsegments.hpp:108-112
    // ```c++
    //     void SetStartIterator(CAlnMixSeq* seq, CAlnMixStarts::iterator iter)
    //     {
    //         m_StartIts.insert(TStartIterators::value_type(seq, iter))
    //             .first->second = iter;
    //     }
    // ```
    fn set_start_iterator(&mut self, seg: SegIx, seq: SeqIx, at: u32) {
        let key = (self.seqs[seq].seq_idx, self.seqs[seq].child_idx);
        self.segs[seg].start_its.insert(key, (seq, at));
    }

    /// `x_IdentifyAlnMixSeq`: the sequence of `key`, added to `m_Seqs` when new.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnseq.cpp:230-262
    /// ```c++
    ///     TBioseqHandleMap::iterator it = m_BioseqHandles.find(bioseq_handle);
    ///     if (it == m_BioseqHandles.end()) {
    ///         // add this bioseq handle
    ///         aln_seq = new CAlnMixSeq();
    ///         m_BioseqHandles[bioseq_handle] = aln_seq;
    ///         ...
    ///         aln_seq->m_DsCnt = 0;
    ///
    ///         // add this sequence
    ///         m_Seqs.push_back(aln_seq);
    ///         ...
    ///     } else {
    ///         aln_seq = it->second;
    ///     }
    /// ```
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

    /// `CAlnMix::Add(ds, fPreserveRows)`: `CAlnMixSequences::Add` and `CAlnMixMatches::Add` of
    /// one HSP (`subject`, `segments`; the query is row 0, the subject row 1).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnseq.cpp:105-112,158-195
    /// ```c++
    /// CAlnMixSequences::Add(const CDense_seg& ds, TAddFlags flags)
    /// {
    ///     m_DsCnt++;
    ///
    ///     vector<CRef<CAlnMixSeq> >& ds_seq = m_DsSeq[&ds];
    ///
    ///     for (CAlnMap::TNumrow row = 0;  row < ds.GetDim();  row++) {
    ///         ...
    ///             x_IdentifyAlnMixSeq(aln_seq, *(ds.GetIds())[row]);
    ///         ...
    ///         if (flags & fPreserveRows) {
    ///             int row_index = aln_seq->m_RowIdx;
    ///             if (row_index == -1) {
    ///                 // initialization
    ///                 aln_seq->m_RowIdx = row;
    ///             } else while (row_index != row) {
    ///                 ...
    ///             }
    ///         }
    ///
    ///         aln_seq->m_DsCnt++;
    ///         ds_seq.push_back(aln_seq);
    ///     }
    /// }
    /// ```
    /// `strand1` and `strand2` are the rows' strands in the Dense-seg (`strands_exist` holds for
    /// BLAST's Dense-segs and `CreateDensegFromStdseg`'s): `Hsp::strands_differ` (a match of one row
    /// has `m_StrandsDiffer` false).
    ///
    /// `SeqKey::Query` is only row 0 and `SeqKey::Subject` only row 1. A subject that is the
    /// query's own Bioseq (`query_bioseq`) is the query's `CAlnMixSeq` in NCBI, so its row 1 takes
    /// the `while (row_index != row)` branch: an `m_AnotherRow` that is a sequence of its own in
    /// `m_Seqs` (as the subject's here) except for its `m_DsIdx`, which is the Dense-seg that made
    /// it, where another subject's is 0 until its first match (`x_SecondRowFits`): the first HSP
    /// added keeps the subject's first row, not the HSP of the longest segment.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnseq.cpp:165-186
    /// ```c++
    ///             } else while (row_index != row) {
    ///                 if (aln_seq->m_AnotherRow) {
    ///                     aln_seq   = aln_seq->m_AnotherRow;
    ///                     row_index = aln_seq->m_RowIdx;
    ///                 } else {
    ///                     CRef<CAlnMixSeq> another_row (new CAlnMixSeq);
    ///
    ///                     another_row->m_BioseqHandle = aln_seq->m_BioseqHandle;
    ///                     another_row->m_SeqId        = aln_seq->m_SeqId;
    ///                     another_row->m_Width        = aln_seq->m_Width;
    ///                     another_row->m_SeqIdx       = aln_seq->m_SeqIdx;
    ///                     another_row->m_ChildIdx     = aln_seq->m_ChildIdx + 1;
    ///                     another_row->m_DsIdx        = int(m_DsCnt);
    ///                     another_row->m_RowIdx       = row;
    ///
    ///                     m_Seqs.push_back(another_row);
    ///
    ///                     aln_seq = aln_seq->m_AnotherRow = another_row;
    ///
    ///                     break;
    ///                 }
    ///             }
    /// ```
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmatch.cpp:120-252
    /// ```c++
    ///     for (CAlnMap::TNumseg seg =0;  seg < ds.GetNumseg();  seg++) {
    ///         len = ds.GetLens()[seg];
    ///         single_chunk = true;
    ///
    ///         for (CAlnMap::TNumrow row1 = 0;  row1 < ds.GetDim();  row1++) {
    ///             if ((start1 = ds.GetStarts()[seg_off + row1]) >= 0) {
    ///                 ...
    ///                 for (CAlnMap::TNumrow row2 = row1+1;
    ///                      row2 < ds.GetDim();  row2++) {
    ///                     if ((start2 = ds.GetStarts()[seg_off + row2]) >= 0) {
    ///                         //match found
    ///                         if (single_chunk) {
    ///                             single_chunk = false;
    ///                             first_non_gapped_row_found = row1;
    ///                         }
    ///                         ...
    ///                         } else {
    ///                             score = len;
    ///                         }
    ///                         ...
    ///                         if (row1 == first_non_gapped_row_found) {
    ///                             CRef<CAlnMixMatch> match(new CAlnMixMatch);
    ///                             match->m_AlnSeq1 = ds_seq[row1];
    ///                             ...
    ///                             match->m_Start1 = start1;
    ///                             match->m_AlnSeq2 = ds_seq[row2];
    ///                             ...
    ///                             match->m_Start2 = start2;
    ///                             match->m_Len = len;
    ///                             match->m_DsIdx = int(m_DsCnt);
    ///                             match->m_StrandsDiffer = false;
    ///                             if (strands_exist) {
    ///                                 if ((strand1 == eNa_strand_minus  &&
    ///                                      strand2 != eNa_strand_minus)  ||
    ///                                     (strand1 != eNa_strand_minus  &&
    ///                                      strand2 == eNa_strand_minus)) {
    ///
    ///                                     match->m_StrandsDiffer = true;
    ///                                 }
    ///                             }
    ///                             match->m_Score = score;
    ///                             ...
    ///                             m_Matches.push_back(match);
    ///                             ...
    ///                         }
    ///                     }
    ///                 }
    ///                 if (single_chunk) {
    ///                     //record it
    ///                     CRef<CAlnMixMatch> match(new CAlnMixMatch);
    ///                     match->m_Score = 0;
    ///                     match->m_AlnSeq1 = ds_seq[row1];
    ///                     ...
    ///                     match->m_Start1 = start1;
    ///                     match->m_AlnSeq2 = 0;
    ///                     match->m_Start2 = 0;
    ///                     match->m_Len = len;
    ///                     ...
    ///                     match->m_DsIdx = int(m_DsCnt);
    ///                     ...
    ///                     m_Matches.push_back(match);
    ///                 }
    ///             }
    ///         }
    ///         seg_off += ds.GetDim();
    ///     }
    /// ```
    ///
    /// With `fForceTranslation` (`flags`) the rows are not preserved (no `m_AnotherRow`), and each
    /// match of two rows first compares its strands with the strand scores of its sequences: a
    /// strand that differs from a sequence's score so far throws, which `DisplaySeqalign` catches
    /// and skips the Seq-align (showalign.cpp:2073-2076). The throw comes after `m_DsCnt++` of the
    /// sequences and of the matches and after the matches of the Dense-seg's segments before it (a
    /// BLAST HSP's first segment holds both rows, so there are none): the Dense-seg still counts
    /// (`m_DsCnt`, the `m_DsIdx` of the later ones), its matches before the throw stay, and none of
    /// the throwing segment's scores are added. Returns whether the Dense-seg was added whole
    /// (false: left out).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmatch.cpp:150-210
    /// ```c++
    ///                         // determine the strand
    ///                         ENa_strand strand1 = eNa_strand_plus;
    ///                         ENa_strand strand2 = eNa_strand_plus;
    ///                         if (strands_exist) {
    ///                             if (ds.GetStrands()[seg_off + row1]
    ///                                 == eNa_strand_minus) {
    ///                                 strand1 = eNa_strand_minus;
    ///                             }
    ///                             if (ds.GetStrands()[seg_off + row2]
    ///                                 == eNa_strand_minus) {
    ///                                 strand2 = eNa_strand_minus;
    ///                             }
    ///                         }
    ///                         ...
    ///                         // in case of fForceTranslation,
    ///                         // check if strands are not mixed by
    ///                         // comparing current strand to the prevailing one
    ///                         if (flags & fForceTranslation  &&
    ///                             (aln_seq1->m_StrandScore > 0  &&
    ///                              strand1 == eNa_strand_minus ||
    ///                              aln_seq1->m_StrandScore < 0  &&
    ///                              strand1 != eNa_strand_minus ||
    ///                              aln_seq2->m_StrandScore > 0  &&
    ///                              strand2 == eNa_strand_minus ||
    ///                              aln_seq2->m_StrandScore < 0  &&
    ///                              strand2 != eNa_strand_minus)) {
    ///                             NCBI_THROW(CAlnException, eMergeFailure,
    ///                                        "CAlnMixMatches::Add(): "
    ///                                        "Unable to mix strands when "
    ///                                        "forcing translation!");
    ///                         }
    ///
    ///                         // add to the prevailing strand
    ///                         aln_seq1->m_StrandScore += (strand1 == eNa_strand_minus ?
    ///                                                     - score : score);
    ///                         aln_seq2->m_StrandScore += (strand2 == eNa_strand_minus ?
    ///                                                     - score : score);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2073-2076
    /// ```c++
    ///                 } catch (const CException& e){
    ///                     _TRACE("Warning: " << e.what());
    ///                     continue;
    ///                 }
    /// ```
    /// So a TBLASTX subject's HSP of the other strand than the subject's first HSP in the same query
    /// frame is not shown in formats 1 to 4, where the description table and outfmt 0 show it (a
    /// deterministic result of NCBI BLAST+ 2.17.0, written as it is: PD-LOSAT-NCBI-DEFECTS rule 2,
    /// DECISIONS T3-6). A BLASTX HSP of a minus query frame (`fNegativeStrand`) never throws: its
    /// query is minus in every HSP of the group and a protein row has no strand.
    ///
    /// A TBLASTX HSP of the query with its own Bioseq (`query_bioseq`, added with
    /// `fForceTranslation`): without `fPreserveRows` the `while (row_index != row)` branch quoted
    /// above is not taken, so both rows are the query's `CAlnMixSeq` (`x_IdentifyAlnMixSeq` finds the
    /// one Bioseq handle): its `m_DsCnt` counts both rows, and the strand check reads its one strand
    /// score for both rows before either is added to it.
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
            // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnseq.cpp:157-160
            // ```c++
            //         // Set the width
            //         if (ds.IsSetWidths()) {
            //             aln_seq->m_Width = ds.GetWidths()[row];
            //         }
            // ```
            // (The widths of `CreateDensegFromStdseg`'s Dense-seg of a TBLASTN or BLASTX HSP: 1 for
            // the protein row, 3 for the nucleotide row; of a TBLASTX HSP none, which
            // `x_ExtendDSWithWidths` sets to 3 for `fForceTranslation`, its lengths in residues;
            // `x_Merge` sets the same, alnmerger.cpp:138-145, `m_IsAA ? 1 : 3`. A Dense-seg of
            // two untranslated rows has none: 1.)
            // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmix.cpp:133-143
            // ```c++
            //     // translate (extend with widths) the dense-seg if necessary
            //     if (flags & fForceTranslation  && !dsp->IsSetWidths()) {
            //         if ( !m_Scope ) {
            //             ...
            //         } else {
            //             m_InputDSs.push_back(x_ExtendDSWithWidths(*dsp));
            //             dsp = m_InputDSs.back();
            //         }
            // ```
            // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmix.cpp:198-220
            // ```c++
            //     if (contains_NA) {
            //         // fix the lengths
            //         ...
            //                 new_lens[numseg] = lens[numseg] / 3;
            //         ...
            //     // add the widths
            //     CDense_seg::TWidths&  new_widths  = new_ds->SetWidths();
            //     new_widths.resize(ds.GetDim(), contains_NA ? 3 : 1);
            // ```
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
    // ```c++
    //     // subject sequences go on separate rows if requested
    //     if (m_MergeFlags & fQuerySeqMergeOnly) {
    //         if (seq2->m_DsIdx) {
    //             if ( !(m_MergeFlags & fTruncateOverlaps) ) {
    //                 if (seq2->m_DsIdx == match->m_DsIdx) {
    //                     return result;
    //                 } else {
    //                     return eForceSeparateRow;
    //                 }
    //             }
    //         } else {
    //             seq2->m_DsIdx = match->m_DsIdx;
    //             if ( !(m_MergeFlags & fTruncateOverlaps) ) {
    //                 return result;
    //             }
    //         }
    //     }
    // ```
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
    /// ```c++
    ///             if (!starts.size()) {
    ///                 // no starts exist yet
    ///                 ...
    ///                 // this seq has not yet been used, set the strand
    ///                 if (m_AlnMixMatches->m_AddFlags & CAlnMixMatches::fForceTranslation) {
    ///                     seq1->m_PositiveStrand = (seq1->m_StrandScore >= 0);
    ///                 } else {
    ///                     seq1->m_PositiveStrand = ! (m_MergeFlags & fNegativeStrand);
    ///                 }
    ///
    ///                 //create the first one
    ///                 seg = new CAlnMixSegment();
    ///                 seg->m_Len = len;
    ///                 seg->m_DsIdx = match->m_DsIdx;
    ///                 starts[start1] = seg;
    ///                 seg->SetStartIterator
    ///                     (seq1, lo_start_i = hi_start_i = starts.begin());
    ///
    ///                 if (m_MergeFlags & fQuerySeqMergeOnly) {
    ///                     seq2->m_DsIdx = match->m_DsIdx;
    ///                 }
    ///                 // DONE!
    ///             } else {
    ///                 // some starts exist already
    ///
    ///                 // look ahead
    ///                 if ((lo_start_i = start_i = starts.lower_bound(start1))
    ///                     == starts.end()  ||
    ///                     start1 < start_i->first) {
    ///                     // the start position does not exist
    ///                     if (lo_start_i != starts.begin()) {
    ///                         --lo_start_i;
    ///                     }
    ///                 }
    ///
    ///                 // look back
    ///                 if (hi_start_i == starts.end()  &&  start_i != lo_start_i) {
    ///                     CAlnMixSegment * prev_seg = lo_start_i->second;
    ///                     if (lo_start_i->first + prev_seg->m_Len * width1 >
    ///                         start1) {
    ///                         // x----..   becomes  x-x--..
    ///                         //   x--..
    ///
    ///                         TSeqPos len1 = (start1 - lo_start_i->first) / width1;
    ///                         TSeqPos len2 = prev_seg->m_Len - len1;
    ///
    ///                         // create the second seg
    ///                         seg = new CAlnMixSegment();
    ///                         seg->m_Len = len2;
    ///                         seg->m_DsIdx = match->m_DsIdx;
    ///                         starts[start1] = seg;
    ///
    ///                         // create rows info
    ///                         ITERATE (CAlnMixSegment::TStartIterators, it,
    ///                                  prev_seg->m_StartIts) {
    ///                             CAlnMixSeq * seq = it->first;
    ///                             tmp_start_i = it->second;
    ///                             if (seq->m_PositiveStrand ==
    ///                                 seq1->m_PositiveStrand) {
    ///                                 seq->SetStarts()
    ///                                     [tmp_start_i->first + len1 * seq->m_Width]
    ///                                     = seg;
    ///                                 if (tmp_start_i != seq->SetStarts().end()) {
    ///                                     seg->SetStartIterator(seq, ++tmp_start_i);
    ///                                 } else {
    ///                                     ...
    ///                                 }
    ///                             } else {
    ///                                 ...
    ///                             }
    ///                         }
    ///
    ///                         // truncate the first seg
    ///                         prev_seg->m_Len = len1;
    ///
    ///                         if (start_i != starts.begin()) {
    ///                             start_i--; // point to the newly created start
    ///                         }
    ///                     }
    ///                     if (lo_start_i != starts.end()) {
    ///                         lo_start_i++;
    ///                     }
    ///                 }
    ///             }
    ///
    ///             // loop through overlapping segments
    ///             start = start1;
    ///             while (hi_start_i == starts.end()) {
    ///                 if (start_i != starts.end()  &&  start_i->first == start) {
    ///                     CAlnMixSegment * prev_seg = start_i->second;
    ///                     if (prev_seg->m_Len > curr_len) {
    ///                         // x-------)  becomes  x----)x--)
    ///                         // x----)
    ///
    ///                         // create the second seg
    ///                         seg = new CAlnMixSegment();
    ///                         TSeqPos len1 =
    ///                             seg->m_Len = prev_seg->m_Len - curr_len;
    ///                         start += curr_len * width1;
    ///
    ///                         // truncate the first seg
    ///                         prev_seg->m_Len = curr_len;
    ///
    ///                         // create rows info
    ///                         ITERATE (CAlnMixSegment::TStartIterators, it,
    ///                                 prev_seg->m_StartIts) {
    ///                             CAlnMixSeq * seq = it->first;
    ///                             tmp_start_i = it->second;
    ///                             if (seq->m_PositiveStrand ==
    ///                                 seq1->m_PositiveStrand) {
    ///                                 seq->SetStarts()[tmp_start_i->first +
    ///                                              curr_len * seq->m_Width]
    ///                                     = seg;
    ///                                 if (tmp_start_i != seq->SetStarts().end()) {
    ///                                     seg->SetStartIterator(seq, ++tmp_start_i);
    ///                                 } else {
    ///                                     ...
    ///                                 }
    ///                             } else{
    ///                                 ...
    ///                             }
    ///                         }
    ///                         ...
    ///                         hi_start_i = start_i; // DONE!
    ///                     } else if (curr_len == prev_seg->m_Len) {
    ///                         // x----)
    ///                         // x----)
    ///                         hi_start_i = start_i; // DONE!
    ///                     } else {
    ///                         // x----)     becomes  x----)x--)
    ///                         // x-------)
    ///                         start += prev_seg->m_Len * width1;
    ///                         curr_len -= prev_seg->m_Len;
    ///                         if (start_i != starts.end()) {
    ///                             start_i++;
    ///                         }
    ///                     }
    ///                 } else {
    ///                     seg = new CAlnMixSegment();
    ///                     starts[start] = seg;
    ///                     tmp_start_i = start_i;
    ///                     if (tmp_start_i != starts.begin()) {
    ///                         tmp_start_i--;
    ///                     }
    ///                     seg->SetStartIterator(seq1, tmp_start_i);
    ///                     if (start_i != starts.end()  &&
    ///                         start + curr_len * width1 > start_i->first) {
    ///                         //       x--..
    ///                         // x--------..
    ///                         seg->m_Len = (start_i->first - start) / width1;
    ///                         seg->m_DsIdx = match->m_DsIdx;
    ///                     } else {
    ///                         //       x-----)
    ///                         // x---)
    ///                         seg->m_Len = curr_len;
    ///                         seg->m_DsIdx = match->m_DsIdx;
    ///                         hi_start_i = start_i;
    ///                         if (hi_start_i != starts.begin()) {
    ///                             hi_start_i--; // DONE!
    ///                         }
    ///                     }
    ///                     start += seg->m_Len * width1;
    ///                     curr_len -= seg->m_Len;
    ///                     if (lo_start_i == start_i) {
    ///                         if (lo_start_i != starts.begin()) {
    ///                             lo_start_i--;
    ///                         }
    ///                     }
    ///                 }
    ///             }
    /// ```
    /// The `else` branches of the two splits (a row of the other strand than `seq1`, whose letters
    /// run the other way: the first part of the segment is its end):
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:550-557 (the look back)
    /// ```c++
    ///                             } else {
    ///                                 seq->SetStarts()
    ///                                     [tmp_start_i->first + len2 * seq->m_Width]
    ///                                     = prev_seg;
    ///                                 seq->SetStarts()[tmp_start_i->first] = seg;
    ///                                 seg->SetStartIterator(seq, tmp_start_i);
    ///                                 if (tmp_start_i != seq->SetStarts().end()) {
    ///                                     prev_seg->SetStartIterator(seq, ++tmp_start_i);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:614-621 (the overlap loop)
    /// ```c++
    ///                             } else{
    ///                                 seq->SetStarts()[tmp_start_i->first +
    ///                                              len1 * seq->m_Width]
    ///                                     = prev_seg;
    ///                                 seq->SetStarts()[tmp_start_i->first] = seg;
    ///                                 seg->SetStartIterator(seq, tmp_start_i);
    ///                                 if (tmp_start_i != seq->SetStarts().end()) {
    ///                                     prev_seg->SetStartIterator(seq, ++tmp_start_i);
    /// ```
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
    /// without a second row would have NCBI write through a null `seq2`: the matches of one row
    /// come after those of two (score 0, `SortByScore`) and every row of a later pass has the
    /// starts of its Dense-seg's matches with the query, so BLAST's HSPs never make one.
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
            // alnmerger.cpp:483-489, quoted above: the sign of the strand score with
            // `fForceTranslation`, else `! (m_MergeFlags & fNegativeStrand)`, which no merge flag
            // here has (module doc).
            self.seqs[seq1].positive_strand = match self.add_flags {
                AddFlags::ForceTranslation => self.seqs[seq1].strand_score >= 0,
                AddFlags::PreserveRows => true,
            };
            let seg = self.new_segment(len);
            self.seqs[seq1].starts.insert(start1, seg);
            self.set_start_iterator(seg, seq1, start1);
            lo_start_i = Some(start1);
            hi_start_i = Some(start1);
            // `if (m_MergeFlags & fQuerySeqMergeOnly) { seq2->m_DsIdx = match->m_DsIdx; }`
            // (alnmerger.cpp:499-501, quoted above).
            let seq2 = seq2.ok_or(MixError::Internal(
                "the first match of a row has no second row",
            ))?;
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
                    let seg = self.new_segment(len2);
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
                    if prev_len > curr_len {
                        let len1 = prev_len - curr_len;
                        let seg = self.new_segment(len1);
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
                    let seg = self.new_segment(0);
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
                            // `hi_start_i = start_i; if (hi_start_i != starts.begin())
                            // hi_start_i--;`: the new start, which is before `start_i`.
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

    /// `x_Merge` for these flags, then `BuildRows`, `Build(false, true, false)`,
    /// `FillUnalignedRegions` and `x_CreateDenseg` (alnmerger.cpp:822-830).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:101-130
    /// ```c++
    ///     // Find the refseq (if such exists)
    ///     {{
    ///         m_SingleRefseq = false;
    ///         m_IndependentDSs = m_DsCnt > 1;
    ///
    ///         unsigned int ds_cnt;
    ///         NON_CONST_ITERATE (TSeqs, it, m_Seqs){
    ///             ds_cnt = (*it)->m_DsCnt;
    ///             if (ds_cnt > 1) {
    ///                 m_IndependentDSs = false;
    ///             }
    ///             if (ds_cnt == m_DsCnt) {
    ///                 m_SingleRefseq = true;
    ///                 if ( !first_refseq ) {
    ///                     CRef<CAlnMixSeq> refseq = *it;
    ///                     m_Seqs.erase(it);
    ///                     m_Seqs.insert(m_Seqs.begin(), refseq);
    ///                 }
    ///                 break;
    ///             }
    ///             first_refseq = false;
    ///         }
    ///     }}
    ///
    ///     // Index the sequences
    ///     {{
    ///         int seq_idx=0;
    ///         NON_CONST_ITERATE (TSeqs, seq_i, m_Seqs) {
    ///             (*seq_i)->m_SeqIdx = seq_idx++;
    ///         }
    ///     }}
    /// ```
    /// The query is the first sequence and in every Dense-seg: it is the single refseq, unless a
    /// TBLASTX HSP of the query's own Bioseq counted it twice (`Mix::add`); then no sequence is
    /// the single refseq and the query, first in `m_Seqs`, stays the refseq of the first pass. Its
    /// `m_DsCnt` is above 1 whenever the mix has more than one Dense-seg (and also with one
    /// Dense-seg of its own Bioseq), so the Dense-segs are never independent.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:160-252
    /// ```c++
    ///     while (true) {
    ///         ...
    ///         // reached the end?
    ///         if (first_refseq ?
    ///             match_i == m_Matches.end()  &&  m_Matches.size() :
    ///             refseq->m_MatchList.empty()  ||
    ///             match_list_i == refseq->m_MatchList.end()) {
    ///
    ///             // move on to the next refseq
    ///             refseq->m_RefBy = 0;
    ///
    ///             // try to find the best scoring 'connected' candidate
    ///             NON_CONST_ITERATE (TSeqs, it, m_Seqs){
    ///                 if ( !((*it)->m_MatchList.empty())  &&
    ///                      (*it)->m_RefBy == refseq) {
    ///                     ...
    ///                     refseq = *it;
    ///                     break;
    ///                 }
    ///             }
    ///             if (refseq->m_RefBy == 0) {
    ///                 // no candidate was found 'connected' to the refseq
    ///                 // continue with the highest scoring candidate
    ///                 NON_CONST_ITERATE (TSeqs, it, m_Seqs){
    ///                     if ( !((*it)->m_MatchList.empty()) ) {
    ///                         ...
    ///                         refseq = *it;
    ///                         break;
    ///                     }
    ///                 }
    ///             }
    ///
    ///             if (refseq->m_MatchList.empty()) {
    ///                 break; // done
    ///             } else {
    ///                 first_refseq = false;
    ///                 match_list_i = refseq->m_MatchList.begin();
    ///             }
    ///             continue;
    ///         } else {
    ///             // iterate
    ///             if (first_refseq) {
    ///                 match = *match_i;
    ///                 ++match_i;
    ///             } else {
    ///                 match = *match_list_i;
    ///                 ++match_list_i;
    ///                 ...
    ///             }
    ///         }
    ///
    ///         curr_len = len = match->m_Len;
    ///
    ///         // is it a match with this refseq?
    ///         if (match->m_AlnSeq1 == refseq) {
    ///             seq1 = match->m_AlnSeq1;
    ///             start1 = match->m_Start1;
    ///             ...
    ///             seq2 = match->m_AlnSeq2;
    ///             start2 = match->m_Start2;
    ///             ...
    ///         } else if (match->m_AlnSeq2 == refseq) {
    ///             ...
    ///         } else {
    ///             seq1 = match->m_AlnSeq1;
    ///             seq2 = match->m_AlnSeq2;
    ///
    ///             // mark the two refseqs, they are candidates for next refseq(s)
    ///             _ASSERT(seq1);
    ///             match->m_MatchIter1 =
    ///                 seq1->m_MatchList.insert(seq1->m_MatchList.end(), match);
    ///             if (seq2) {
    ///                 ...
    ///             }
    ///             ...
    ///             // mark that there's no match with this refseq
    ///             seq1 = 0;
    ///         }
    /// ```
    /// The query's pass takes every match in score order; a match of a subject's letters alone
    /// (an insertion) goes to that subject's `m_MatchList`, which the next passes take, one
    /// subject at a time in the order of `m_Seqs`. `m_AlnSeq2` is never the refseq: the query is
    /// always `m_AlnSeq1`, and a subject's matches in its list have no second row. A match of the
    /// query with its own Bioseq (`m_AlnSeq1 == m_AlnSeq2`) is the query's in the first pass, where
    /// the `m_AlnSeq1 == refseq` branch takes it; its second increment of `match_list_i` belongs
    /// to the later passes, which a mix without a single refseq never reaches (TBLASTX's HSPs
    /// have no gaps, so no match of one row).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:470-481
    /// ```c++
    ///                 if ( !m_IndependentDSs ) {
    ///                     // Stop the algorithm after merging with respect to the first refseq.
    ///                     // Clear all MatchLists and exit
    ///                     if ( !(m_SingleRefseq  ||  first_refseq) ) {
    ///                         NON_CONST_ITERATE (TSeqs, it, m_Seqs){
    ///                             if ( !((*it)->m_MatchList.empty()) ) {
    ///                                 (*it)->m_MatchList.clear();
    ///                             }
    ///                         }
    ///                         break;
    ///                     }
    ///                 }
    /// ```
    /// (In the `!starts.size()` branch of a later pass: not reached, as above; an internal error
    /// if it were.)
    fn merge(&mut self) -> Result<(), MixError> {
        // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmix.cpp:279-284
        // ```c++
        //     if (flags & fSortInputByScore) {
        //         m_AlnMixMatches->SortByChainScore();
        //     } else {
        //         m_AlnMixMatches->SortByScore();
        //     }
        // ```
        // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmatch.cpp:62-87
        // ```c++
        //     return match1->m_Score > match2->m_Score;
        //     ...
        //     stable_sort(m_Matches.begin(), m_Matches.end(), x_CompareScores);
        // ```
        // (A stable sort of indices: `sort_by_key` keeps equal scores in input order.)
        let mut order: Vec<usize> = (0..self.matches.len()).collect();
        order.sort_by_key(|&index| std::cmp::Reverse(self.matches[index].score));
        // `Find the refseq (if such exists)`, quoted above.
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
                    // (`seq1->m_MatchList.erase(match_list_iter1)`).
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
            if matched.len == 0 {
                // A Dense-seg segment of no letters: BLAST's Seq-aligns have none.
                return Err(MixError::Internal("a segment of length 0"));
            }
            self.merge_match(matched, first_refseq)?;
        }
        Ok(())
    }

    /// One match of the refseq (`seq1`), from `// save the match info into the segments map`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:349-380
    /// ```c++
    ///             // if a subject sequence place it in the proper row
    ///             if ( !first_refseq  &&  m_MergeFlags & fQuerySeqMergeOnly) {
    ///                 bool proper_row_found = false;
    ///                 while (true) {
    ///                     if (seq1->m_DsIdx == match->m_DsIdx) {
    ///                         proper_row_found = true;
    ///                         break;
    ///                     } else {
    ///                         if (seq1->m_ExtraRow) {
    ///                             seq1 = match->m_AlnSeq1 = seq1->m_ExtraRow;
    ///                         } else {
    ///                             break;
    ///                         }
    ///                     }
    ///                 }
    ///                 if ( !proper_row_found   &&
    ///                      !(m_MergeFlags & fTruncateOverlaps)) {
    ///                     NCBI_THROW(CAlnException, eMergeFailure,
    ///                                "CAlnMixMerger::x_Merge(): "
    ///                                "Proper row not found for the match. "
    ///                                "Cannot use fQuerySeqMergeOnly?");
    ///                 }
    ///             }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:381-416
    /// ```c++
    ///             CAlnMixStarts& starts = seq1->SetStarts();
    ///             if (seq2) {
    ///                 // mark it, it is linked to the refseq
    ///                 seq2->m_RefBy = refseq;
    ///                 ...
    ///             }
    ///             ...
    ///             if (seq2) {
    ///                 if (width2 == 3) {
    ///                     // Set the frame if necessary
    ///                     x_SetSeqFrame(match, match->m_AlnSeq2);
    ///                 }
    ///                 // check if the second row fits
    ///                 // this will truncate the match if
    ///                 // there's an inconsistent overlap
    ///                 // and truncation was requested
    ///                 second_row_fits = x_SecondRowFits(match);
    ///                 ...
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:705-818
    /// ```c++
    ///             // try to resolve the second row
    ///             if (seq2) {
    ///
    ///                 while (second_row_fits != eSecondRowFitsOk  &&
    ///                        second_row_fits != eIgnoreMatch) {
    ///                     if (!seq2->m_ExtraRow) {
    ///                         // create an extra row
    ///                         CRef<CAlnMixSeq> row (new CAlnMixSeq);
    ///                         row->m_BioseqHandle = seq2->m_BioseqHandle;
    ///                         row->m_SeqId = seq2->m_SeqId;
    ///                         row->m_Width = seq2->m_Width;
    ///                         row->m_Frame = start2 % 3;
    ///                         row->m_SeqIdx = seq2->m_SeqIdx;
    ///                         row->m_ChildIdx = seq2->m_ChildIdx + 1;
    ///                         if (m_MergeFlags & fQuerySeqMergeOnly) {
    ///                             row->m_DsIdx = match->m_DsIdx;
    ///                         }
    ///                         m_ExtraRows.push_back(row);
    ///                         row->m_ExtraRowIdx = seq2->m_ExtraRowIdx + 1;
    ///                         seq2 = match->m_AlnSeq2 = seq2->m_ExtraRow = row;
    ///                         break;
    ///                     }
    ///                     seq2 = match->m_AlnSeq2 = seq2->m_ExtraRow;
    ///
    ///                     second_row_fits = x_SecondRowFits(match);
    ///                     ...
    ///                 }
    ///                 ...
    ///                 // set the strand if first time
    ///                 if (seq2->GetStarts().empty()) {
    ///                     seq2->m_PositiveStrand =
    ///                         (seq1->m_PositiveStrand ?
    ///                          !match->m_StrandsDiffer :
    ///                          match->m_StrandsDiffer);
    ///                 }
    ///
    ///                 // create row info
    ///                 CAlnMixStarts& starts2 = match->m_AlnSeq2->SetStarts();
    ///                 start = start2;
    ///                 CAlnMixStarts::iterator start2_i
    ///                     = starts2.lower_bound(start2);
    ///                 start_i = match->m_StrandsDiffer ? hi_start_i : lo_start_i;
    ///
    ///                 while(start < start2 + len * width2) {
    ///                     if (start2_i != starts2.end() &&
    ///                         start2_i->first == start) {
    ///                         // this position already exists
    ///                         ...
    ///                     } else {
    ///                         // this position does not exist, create it
    ///                         seq2->SetStarts()[start] = start_i->second;
    ///
    ///                         // start2_i != starts.begin() because we just
    ///                         // made an insertion, so decrement is ok
    ///                         start2_i--;
    ///
    ///                         // point this segment's row start iterator
    ///                         if (start_i->second->m_StartIts.find(seq2) !=
    ///                             start_i->second->m_StartIts.end()) {
    ///                             // consistency check fails
    ///                             NCBI_THROW(CAlnException, eMergeFailure,
    ///                                        "CAlnMixMerger::x_Merge(): "
    ///                                        "Internal error: "
    ///                                        "Start iterator already exists for seq2.");
    ///                         }
    ///                         start_i->second->SetStartIterator(seq2, start2_i);
    ///                         ...
    ///                     }
    ///
    ///                     // increment values
    ///                     start += start_i->second->m_Len * width2;
    ///                     if (start2_i != starts2.end()) {
    ///                         start2_i++;
    ///                     }
    ///                     if (match->m_StrandsDiffer) {
    ///                         if (start_i != starts.begin()) {
    ///                             start_i--;
    ///                         }
    ///                     } else {
    ///                         if (start_i != starts.end()) {
    ///                             start_i++;
    ///                         }
    ///                     }
    ///                 }
    ///             }
    /// ```
    /// A width-3 refseq (a TBLASTN subject's row in a later pass, its insertions) first goes to the
    /// row of its frame (alnmerger.cpp:335-349, `x_SetSeqFrame(match, match->m_AlnSeq1)`), from
    /// which the proper row is looked for. Every row of a TBLASTN HSP has its frame (`start % 3`
    /// of each of its segments is the same), so that row is at or before the HSP's own.
    ///
    /// The row of an ordinary subject's HSP holds that HSP's aligned segments alone, which never
    /// share a subject position: there the walk never meets a position of the row.
    ///
    /// The query's own Bioseq added with `fForceTranslation` (TBLASTX, `Mix::add`): a match of
    /// such an HSP has the query's sequence on both sides. `x_SetSeqFrame` of the second row
    /// compares the row with the match's first row (`seq == match->m_AlnSeq1`, both the query's
    /// row of the group's frame) and takes the frame of the query's start, so the second row is
    /// the query's own row; `x_SecondRowFits` gives it the first such HSP (its `m_DsIdx`) and
    /// every later one a new extra row of the query. The first one's walk runs on the row of the
    /// query's starts: at a start that exists it joins the two segments (`this position already
    /// exists`, which is how the record's HSP with itself in the same frame folds into the query's
    /// row), and a start that does not exist throws, because the query's segment there already has
    /// the row's start iterator (`MixError::NcbiThrows`).
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
                    None => return Err(MixError::Internal("proper row not found for the match")),
                }
            }
        }
        let mut seq2 = matched.seq2;
        let mut fits = true;
        if let Some(mut row) = seq2 {
            self.seqs[row].ref_by = Some(self.order[0]);
            if self.seqs[row].width == 3 {
                // `frame = (seq == match->m_AlnSeq1 ? m_Start1 : m_Start2) % 3` (alnmerger.cpp:
                // 1272-1276, quoted at `set_seq_frame`): the query's start where the second row
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
                // ```c++
                //                     if (start2_i != starts2.end() &&
                //                         start2_i->first == start) {
                //                         // this position already exists
                //
                //                         if (start2_i->second != start_i->second) {
                //                             // merge the two segments
                //
                //                             // store the seg in a CRef to delay its deletion
                //                             // until after the iteration on it is finished
                //                             CRef<CAlnMixSegment> tmp_seg(start2_i->second);
                //
                //                             ITERATE (CAlnMixSegment::TStartIterators,
                //                                      it,
                //                                      tmp_seg->m_StartIts) {
                //                                 CAlnMixSeq * tmp_seq = it->first;
                //                                 tmp_start_i = it->second;
                //                                 tmp_start_i->second = start_i->second;
                //                                 start_i->second->SetStartIterator(tmp_seq, tmp_start_i);
                //                             }
                // ```
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
                    return Err(MixError::NcbiThrows(
                        "CAlnMixMerger::x_Merge(): Internal error: Start iterator already exists \
                         for seq2.",
                    ));
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
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:1269-1302
    /// ```c++
    /// void
    /// CAlnMixMerger::x_SetSeqFrame(CAlnMixMatch* match, CAlnMixSeq*& seq)
    /// {
    ///     unsigned frame;
    ///     if (seq == match->m_AlnSeq1) {
    ///         frame = match->m_Start1 % 3;
    ///     } else {
    ///         frame = match->m_Start2 % 3;
    ///     }
    ///     if (seq->GetStarts().empty()) {
    ///         seq->m_Frame = frame;
    ///     } else {
    ///         while ((unsigned)seq->m_Frame != frame) {
    ///             if (!seq->m_ExtraRow) {
    ///                 // create an extra frame
    ///                 CRef<CAlnMixSeq> new_seq (new CAlnMixSeq);
    ///                 new_seq->m_BioseqHandle = seq->m_BioseqHandle;
    ///                 new_seq->m_SeqId = seq->m_SeqId;
    ///                 new_seq->m_PositiveStrand = seq->m_PositiveStrand;
    ///                 new_seq->m_Width = seq->m_Width;
    ///                 new_seq->m_Frame = frame;
    ///                 new_seq->m_SeqIdx = seq->m_SeqIdx;
    ///                 new_seq->m_ChildIdx = seq->m_ChildIdx + 1;
    ///                 if (m_MergeFlags & fQuerySeqMergeOnly) {
    ///                     new_seq->m_DsIdx = match->m_DsIdx;
    ///                 }
    ///                 m_ExtraRows.push_back(new_seq);
    ///                 new_seq->m_ExtraRowIdx = seq->m_ExtraRowIdx + 1;
    ///                 seq = seq->m_ExtraRow = new_seq;
    ///                 break;
    ///             }
    ///             seq = seq->m_ExtraRow;
    ///         }
    ///     }
    /// }
    /// ```
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

    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnseq.cpp:266-283
    /// ```c++
    /// void
    /// CAlnMixSequences::BuildRows()
    /// {
    ///     m_Rows.clear();
    ///
    ///     int count = 0;
    ///     NON_CONST_ITERATE (TSeqs, i, m_Seqs) {
    ///         CRef<CAlnMixSeq>& seq = *i;
    ///
    ///         if ( !seq->GetStarts().empty() ) {
    ///             m_Rows.push_back(seq);
    ///             seq->m_RowIdx = count++;
    ///             while ((seq = seq->m_ExtraRow) != NULL ) {
    ///                 seq->m_RowIdx = count++;
    ///                 m_Rows.push_back(seq);
    ///             }
    ///         }
    ///     }
    /// }
    /// ```
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

    /// `Build(gap_join = false, min_gap = true, remove_leading_and_trailing_gaps = false)`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnsegments.cpp:58-115,171-271
    /// ```c++
    ///     m_AlnMixSequences->InitRowsStartIts();
    ///     m_AlnMixSequences->InitExtraRowsStartIts();
    ///     ...
    ///     TSegmentsContainer gapped_segs;
    ///
    ///     CAlnMixSequences::TSeqs::iterator refseq_it = m_Rows.begin();
    ///     bool orig_refseq = true;
    ///     while (true) {
    ///         CAlnMixSeq * refseq = 0;
    ///         while (refseq_it != m_Rows.end()) {
    ///             refseq = *(refseq_it++);
    ///             if (refseq->GetStarts().current != refseq->GetStarts().end()) {
    ///                 break;
    ///             } else {
    ///                 refseq = 0;
    ///             }
    ///         }
    ///         if ( !refseq ) {
    ///             // Done
    ///
    ///             // add the gapped segments if any
    ///             if (gapped_segs.size()) {
    ///                 if (gap_join) {
    ///                     ...
    ///                 } else if (min_gap) {
    ///                     // request to try to align
    ///                     // all gapped segments
    ///                     x_MinimizeGaps(gapped_segs);
    ///                 }
    ///                 NON_CONST_ITERATE (TSegmentsContainer,
    ///                                    seg_i, gapped_segs) {
    ///                     m_Segments.push_back(&**seg_i);
    ///                 }
    ///                 gapped_segs.clear();
    ///             }
    ///             break; // from the main loop
    ///         }
    ///         ...
    ///         // for each refseq segment
    ///         while (refseq->GetStarts().current != refseq->GetStarts().end()) {
    ///             stack< CRef<CAlnMixSegment> > seg_stack;
    ///             // To prevent infinite loop, remember all segments already on the stack.
    ///             set< CRef<CAlnMixSegment> > seg_set;
    ///             seg_stack.push(refseq->GetStarts().current->second);
    ///             seg_set.insert(refseq->GetStarts().current->second);
    ///             ...
    ///             while ( !seg_stack.empty() ) {
    ///
    ///                 bool pop_seg = true;
    ///                 // check the gapped segments on the left
    ///                 ITERATE (CAlnMixSegment::TStartIterators, start_its_i,
    ///                          seg_stack.top()->m_StartIts) {
    ///
    ///                     CAlnMixSeq * row = start_its_i->first;
    ///
    ///                     if (row->GetStarts().current != start_its_i->second) {
    ///                         ...
    ///                         seg_stack.push(row->GetStarts().current->second);
    ///                         // If the segment is already in the set, we have infinite loop
    ///                         // and must terminate it.
    ///                         if (!seg_set.insert(row->GetStarts().current->second).second) {
    ///                             string errstr =
    ///                                 string("CAlnMixSegments::Build():")
    ///                                 + " Internal error: Infinite loop detected.";
    ///                             NCBI_THROW(CAlnException, eMergeFailure, errstr);
    ///                         }
    ///                         ...
    ///                         pop_seg = false;
    ///                         break;
    ///                     }
    ///                 }
    ///
    ///                 if (pop_seg) {
    ///
    ///                     // inc/dec iterators for each row of the seg
    ///                     ITERATE (CAlnMixSegment::TStartIterators, start_its_i,
    ///                              seg_stack.top()->m_StartIts) {
    ///                         CAlnMixSeq * row = start_its_i->first;
    ///                         ...
    ///                         if (row->m_PositiveStrand) {
    ///                             row->SetStarts().current++;
    ///                         } else {
    ///                             ...
    ///                         }
    ///                     }
    ///
    ///                     if (seg_stack.size() > 1) {
    ///                         // add to the gapped segments
    ///                         gapped_segs.push_back(seg_stack.top());
    ///                         seg_stack.pop();
    ///                         ...
    ///                     } else {
    ///                         // add the gapped segments if any
    ///                         if (gapped_segs.size()) {
    ///                             if (gap_join) {
    ///                                 ...
    ///                             } else if (min_gap) {
    ///                                 // request to try to align
    ///                                 // all gapped segments
    ///                                 x_MinimizeGaps(gapped_segs);
    ///                             }
    ///                             if (orig_refseq) {
    ///                                 NON_CONST_ITERATE (TSegmentsContainer,
    ///                                                    seg_i, gapped_segs) {
    ///                                     m_Segments.push_back(&**seg_i);
    ///                                 }
    ///                                 gapped_segs.clear();
    ///                             }
    ///                         }
    ///                         // add the refseq segment
    ///                         if (orig_refseq) {
    ///                             m_Segments.push_back(seg_stack.top());
    ///                         } else {
    ///                             gapped_segs.push_back(seg_stack.top());
    ///                         }
    ///                         seg_stack.pop();
    ///                         ...
    ///                     } // if (seg_stack.size() > 1)
    ///                 } // if (popseg)
    ///             } // while ( !seg_stack.empty() )
    ///         } // while (refseq->GetStarts().current != refseq->GetStarts().end())
    ///         orig_refseq = false;
    ///     } // while (true)
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnseq.cpp:286-293
    /// ```c++
    ///     NON_CONST_ITERATE (TSeqs, row_i, m_Rows) {
    ///         CAlnMixSeq * row = *row_i;
    ///         if ( !row->GetStarts().empty() ) {
    ///             if (row->m_PositiveStrand) {
    ///                 row->SetStarts().current = row->GetStarts().begin();
    ///             } else {
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnseq.cpp:286-296
    /// ```c++
    ///             if (row->m_PositiveStrand) {
    ///                 row->SetStarts().current = row->GetStarts().begin();
    ///             } else {
    ///                 row->SetStarts().current = row->GetStarts().end();
    ///                 row->SetStarts().current--;
    ///             }
    ///         } else {
    ///             row->SetStarts().current = row->GetStarts().end();
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnsegments.cpp:216-224
    /// ```c++
    ///                         if (row->m_PositiveStrand) {
    ///                             row->SetStarts().current++;
    ///                         } else {
    ///                             if (row->SetStarts().current == row->GetStarts().begin()) {
    ///                                 row->SetStarts().current = row->GetStarts().end();
    ///                             } else {
    ///                                 row->SetStarts().current--;
    ///                             }
    ///                         }
    /// ```
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
                                // NCBI's throw (quoted above), which only the query's own Bioseq
                                // with `fForceTranslation` can reach (`merge_match`).
                                return Err(MixError::NcbiThrows(
                                    "CAlnMixSegments::Build(): Internal error: Infinite loop \
                                     detected.",
                                ));
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

    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnsegments.cpp:438-523
    /// ```c++
    ///     seg_i_begin = seg_i_end = seg_i = gapped_segs.begin();
    ///
    ///     typedef map<TSeqPos, CRef<CAlnMixSegment> > TLenMap;
    ///     TLenMap len_map;
    ///
    ///     while (seg_i_end != gapped_segs.end()) {
    ///
    ///         len_map[(*seg_i_end)->m_Len];
    ///
    ///         // check if this seg can possibly be minimized
    ///         bool possible = true;
    ///
    ///         seg_i = seg_i_begin;
    ///         while (seg_i != seg_i_end) {
    ///             seg1 = *seg_i;
    ///             seg2 = *seg_i_end;
    ///
    ///             ITERATE (CAlnMixSegment::TStartIterators,
    ///                      st_it,
    ///                      seg2->m_StartIts) {
    ///                 seq = st_it->first;
    ///                 // check if this seq was already used
    ///                 if (seg1->m_StartIts.find(seq) != seg1->m_StartIts.end()) {
    ///                     possible = false;
    ///                     break;
    ///                 }
    ///             }
    ///             if ( !possible ) {
    ///                 break;
    ///             }
    ///             seg_i++;
    ///         }
    ///         seg_i_end++;
    ///
    ///         if ( !possible  ||  seg_i_end == gapped_segs.end()) {
    ///             // use the accumulated len info to create the new segs
    ///
    ///             // create new segs with appropriate len
    ///             TSeqPos len_so_far = 0;
    ///             TLenMap::iterator len_i = len_map.begin();
    ///             while (len_i != len_map.end()) {
    ///                 len_i->second = new CAlnMixSegment();
    ///                 len_i->second->m_Len = len_i->first - len_so_far;
    ///                 len_so_far += len_i->second->m_Len;
    ///                 len_i++;
    ///             }
    ///
    ///             // loop through the accumulated orig segs.
    ///             // copy info from them into the new segs
    ///             TLenMap::iterator len_i_end;
    ///             seg_i = seg_i_begin;
    ///             while (seg_i != seg_i_end) {
    ///                 TSeqPos orig_len = (*seg_i)->m_Len;
    ///
    ///                 // determine the span of the current seg
    ///                 len_i_end = len_map.find(orig_len);
    ///                 len_i_end++;
    ///
    ///                 // loop through its sequences
    ///                 NON_CONST_ITERATE (CAlnMixSegment::TStartIterators,
    ///                                    st_it,
    ///                                    (*seg_i)->m_StartIts) {
    ///
    ///                     seq = st_it->first;
    ///                     TSeqPos orig_start = st_it->second->first;
    ///
    ///                     len_i = len_map.begin();
    ///                     len_so_far = 0;
    ///                     // loop through the new segs
    ///                     while (len_i != len_i_end) {
    ///                         seg = len_i->second;
    ///
    ///                         // calc the start
    ///                         TSeqPos this_start = orig_start +
    ///                             (seq->m_PositiveStrand ?
    ///                              len_so_far :
    ///                              orig_len - len_so_far - seg->m_Len) *
    ///                             seq->m_Width;
    ///
    ///                         // create the bindings:
    ///                         seq->SetStarts()[this_start] = seg;
    ///                         seg->SetStartIterator(seq, seq->SetStarts().find(this_start));
    ///                         len_i++;
    ///                         len_so_far += seg->m_Len;
    ///                     }
    ///                 }
    ///                 seg_i++;
    ///             }
    ///             NON_CONST_ITERATE (TLenMap, len_it, len_map) {
    ///                 new_segs.push_back(len_it->second);
    ///             }
    ///             len_map.clear();
    ///             seg_i_begin = seg_i_end;
    ///         }
    ///     }
    ///     gapped_segs.clear();
    ///     ITERATE (TSegmentsContainer, new_seg_i, new_segs) {
    ///         gapped_segs.push_back(*new_seg_i);
    ///     }
    /// ```
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
    /// ```c++
    ///     vector<TSignedSeqPos> starts;
    ///     vector<TSeqPos> lens;
    ///     starts.resize(m_Rows.size(), -1);
    ///     lens.resize(m_Rows.size(), 0);
    ///
    ///     TSeqPos len = 0;
    ///     CAlnMap::TNumrow rowidx;
    ///
    ///     TSegments::iterator seg_i = m_Segments.begin();
    ///     while (seg_i != m_Segments.end()) {
    ///         len = (*seg_i)->m_Len;
    ///         ITERATE (CAlnMixSegment::TStartIterators, start_its_i,
    ///                  (*seg_i)->m_StartIts) {
    ///             CAlnMixSeq * row = start_its_i->first;
    ///             rowidx = row->m_RowIdx;
    ///             TSignedSeqPos& prev_start = starts[rowidx];
    ///             TSeqPos& prev_len = lens[rowidx];
    ///             TSeqPos start = start_its_i->second->first;
    ///             const bool plus = row->m_PositiveStrand;
    ///             const int& width = row->m_Width;
    ///             TSeqPos prev_start_plus_len = prev_start + prev_len * width;
    ///             TSeqPos start_plus_len = start + len * width;
    ///             if (prev_start >= 0) {
    ///                 if (plus  &&  prev_start_plus_len < start  ||
    ///                     !plus  &&  start_plus_len < (TSeqPos) prev_start) {
    ///                     // create a new seg
    ///                     CRef<CAlnMixSegment> seg (new CAlnMixSegment);
    ///                     TSeqPos new_start;
    ///                     if (row->m_PositiveStrand) {
    ///                         new_start = prev_start + prev_len * width;
    ///                         seg->m_Len = (start - new_start) / width;
    ///                     } else {
    ///                         new_start = start_plus_len;
    ///                         seg->m_Len = (prev_start - new_start) / width;
    ///                     }
    ///                     row->SetStarts()[new_start] = seg;
    ///                     CAlnMixStarts::iterator start_i =
    ///                         start_its_i->second;
    ///                     seg->SetStartIterator(row,
    ///                                           row->m_PositiveStrand ?
    ///                                           --start_i :
    ///                                           ++start_i);
    ///
    ///                     seg_i = m_Segments.insert(seg_i, seg);
    ///                     seg_i++;
    ///                 }
    ///             }
    ///             prev_start = start;
    ///             prev_len = len;
    ///         }
    ///         seg_i++;
    ///     }
    /// ```
    /// The new segment's start iterator (`--start_i`, or `++start_i` on the minus strand) is its own
    /// start: no start of the row lies between it and the segment being read. The fill of a row's
    /// stretch goes before the segment being read, whose own rows are then
    /// still read from it (`seg_i` points to it again after the insertion).
    ///
    /// A row that `BuildRows` left out keeps the `m_RowIdx` of -1 of `CAlnMixSeq()`: the extra
    /// rows of a sequence without starts, the query's own Bioseq whose `m_AnotherRow` kept none
    /// (`place`). NCBI's `starts[rowidx]` and `lens[rowidx]` then write before their vectors'
    /// buffers (`prev_start = start; prev_len = len;`, a length of at least 1), and glibc stops
    /// NCBI when it frees them (`double free or corruption (out)`, exit 134): no section to
    /// write, `MixError::NcbiAborts` (line RG T3a, PD-LOSAT-NCBI-DEFECTS rule 3).
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
                    self.seqs[row].starts.insert(new_start, seg);
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
    /// ```c++
    ///     ids.resize(numrows);
    ///     lens.resize(numsegs);
    ///     starts.resize(num, -1);
    ///     strands.resize(num, eNa_strand_minus);
    ///
    ///     vector<bool> row_empty(numrows, true);
    ///
    ///     // ids
    ///     numrow = 0;
    ///     ITERATE(TSeqs, row_i, m_Rows) {
    ///         ids[numrow++] = (*row_i)->m_SeqId;
    ///     }
    ///
    ///     int offset = 0;
    ///     numseg = 0;
    ///
    ///     ITERATE(CAlnMixSegments::TSegments,
    ///             seg_i,
    ///             m_AlnMixSegments->m_Segments) {
    ///
    ///         // lens
    ///         lens[numseg] = (*seg_i)->m_Len;
    ///
    ///         // starts
    ///         ITERATE (CAlnMixSegment::TStartIterators, start_its_i,
    ///                 (*seg_i)->m_StartIts) {
    ///             starts[offset + start_its_i->first->m_RowIdx] =
    ///                 start_its_i->second->first;
    ///
    ///             if (start_its_i->second->first != kInvalidSeqPos) {
    ///                 row_empty[start_its_i->first->m_RowIdx] = false;
    ///             }
    ///         }
    ///         ...
    ///         // next segment
    ///         numseg++; offset += numrows;
    ///         ...
    ///     }
    ///     ...
    ///     //remove empty rows
    ///     for(int row = numrows-1; row >=0; --row) {
    ///         if (row_empty[row]) {
    ///             ids.erase(ids.begin()+row);
    ///             ...
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:1215-1240
    /// ```c++
    ///         // strands
    ///         numrow = 0;
    ///         ITERATE(TSeqs, row_i, m_Rows) {
    ///             if ((*row_i)->m_PositiveStrand) {
    ///                 strands[offset + numrow] = eNa_strand_plus;
    ///             }
    ///             numrow++;
    ///         }
    ///         ...
    ///     // widths
    ///     CDense_seg::TWidths*  widths = NULL;
    ///     if ((m_AlnMixMatches->m_ContainsNA  &&  m_AlnMixMatches->m_ContainsAA)  ||
    ///         (m_AlnMixMatches->m_AddFlags & CAlnMixMatches::fForceTranslation)) {
    ///         widths = &m_DS->SetWidths();
    ///         widths->resize(numrows);
    ///         numrow = 0;
    ///         ITERATE (TSeqs, row_i, m_Rows) {
    ///             (*widths)[numrow++] = (*row_i)->m_Width;
    ///         }
    ///     }
    /// ```
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
        // showalign.cpp:2100-2101, quoted at `super::write_section`). The query's first row is
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
/// `CAlnMix::Add` of each with the group's flag (`Hsps::add_flags`), then `Merge(fMinGap |
/// fQuerySeqMergeOnly | fFillUnalignedRegions)` and `GetDenseg`; `None` when every `Add` threw
/// (`hasAln` is false: no merge, and `GetDenseg` throws for the group).
pub(crate) fn merge(hsps: &Hsps, group: &[usize]) -> Result<Option<Denseg>, MixError> {
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
        query_width: hsps.query_width,
        subject_width: hsps.subject_width,
        add_flags: AddFlags::PreserveRows,
        query_bioseq: hsps.query_bioseq,
        single_refseq: true,
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
