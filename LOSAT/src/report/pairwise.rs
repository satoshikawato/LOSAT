//! Pairwise alignment output (outfmt 0)
//!
//! This module implements the traditional BLAST pairwise alignment output format.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp

use super::html_ids::HtmlHeading;
use super::outfmt6::{format_bitscore_ncbi, format_evalue_ncbi, format_total_bitscore_ncbi};
use super::serial::asn_text::itoa_i32;
use super::serial::seqalign::{AlignRow, RowRangeError};
use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::common::{GapEditOp, Hit};
use crate::config::ScoringMatrix;
use crate::ncbi_str::truncate_spaces;
use crate::stats::{BlastGumbelBlk, KarlinParams};
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3563-3581
use crate::report::ostream::ReportWrite;
use crate::utils::matrix::protein_display_score;
use std::io::{self, Write};
use std::sync::Arc;

// psiblast's report of one round (line PS, PS5a): a child of this module, which writes through
// its parts.
#[path = "psiblast_round.rs"]
pub(crate) mod psiblast_round;

// =============================================================================
// NCBI Pairwise Output Format (outfmt 0)
// =============================================================================

/// Line length for alignment display (NCBI default: 60)
pub const DEFAULT_LINE_LENGTH: usize = 60;

/// Configuration for pairwise output
#[derive(Debug, Clone)]
pub struct PairwiseConfig {
    /// Line length for sequence display
    pub line_length: usize,
    /// Show frame information (for translated searches)
    pub show_frame: bool,
    /// Program type for proper formatting
    pub program: String,
    /// Protein scoring matrix for positives/midline rendering.
    pub protein_matrix: ScoringMatrix,
    /// `-sorthits` and `-sorthsps` (NCBI's `m_HitsSortOption`, `m_HspsSortOption`): the order of
    /// each query's hits and HSPs (`hit_sort::sort_query_hits`) and, with `-sorthits`, the
    /// columns of the description table (`SortColumns`). None for psiblast and rpsblast, whose
    /// applications do not hand them to the formatter.
    pub sort: super::hit_sort::ReportSort,
    /// The alignment section: outfmt 0's alignments, or the query-anchored section of outfmt 1
    /// to 4. The rest of the report is the same.
    pub view: PairwiseView,
    /// `-html` (NCBI's `CBlastFormat::m_IsHTML`, blast_format.cpp:90,224): the HTML form of the
    /// report, with the run's links and the report's OIDs (`html_ids::HtmlReport`, made only with
    /// `-html`); `None`: the plain report.
    pub html: Option<super::html_ids::HtmlReport>,
    /// A query whose sequence NCBI's alignment display cannot read (psiblast's and tblastn's PSSM
    /// query of the molecule `other`): the alignments of outfmt 0 are one line each, the
    /// query-anchored section of outfmt 1 to 4 a line feed where NCBI's merge of its HSPs
    /// returns, fails with NCBI's `CAlnException` where that merge throws, and is refused where
    /// that merge reads past a row's end, crashes, reads freed memory, or wraps a segment length
    /// (`query_anchored::write_unreadable_section`).
    /// `None` for every other report.
    pub unreadable_query: Option<UnreadableQuery>,
}

impl Default for PairwiseConfig {
    fn default() -> Self {
        Self {
            line_length: DEFAULT_LINE_LENGTH,
            show_frame: true,
            program: "tblastx".to_string(),
            protein_matrix: ScoringMatrix::Blosum62,
            sort: Default::default(),
            view: PairwiseView::Pairwise,
            html: None,
            unreadable_query: None,
        }
    }
}

impl PairwiseConfig {
    /// The Seq-id text of NCBI's line in place of each alignment of the subject `s_idx` where the
    /// query's sequence cannot be read (`UnreadableQuery`); `None` for a readable query.
    fn unreadable_query_id(
        &self,
        parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
        s_idx: u32,
    ) -> Option<Vec<u8>> {
        self.unreadable_query
            .as_ref()
            .map(|query| query.skipped_id(parsed, s_idx))
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:108-140
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3956-3966
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1950-1981
/// A query whose sequence NCBI's alignment display cannot read: a PSSM query of the molecule
/// `other` has the sequence data of a protein under a molecule that is none, so the sequence
/// vector of its row, which the display reads before it writes anything of the alignment
/// (`x_PrepareRowData`), cannot convert its letters and throws. The catch of the display writes
/// `Sequence with id <id> no longer exists in database...alignment skipped` for each HSP, no
/// heading and no score lines: the same line as where a heading's title cannot be made
/// (`write_skipped_alignment`), for every subject, so the Seq-id text of each subject's
/// Seq-align id is the one thing to know of it (`skipped_id`).
#[derive(Clone, Debug, Default)]
pub struct UnreadableQuery {
    /// The OIDs of a database's listed subjects by subject index (`ListedSubjects::ordinals` in
    /// the opened database), for the ones without parsed Seq-ids, whose id is
    /// `gnl|BL_ORD_ID|<oid>`. `None` for `-subject` records.
    pub oids: Option<Arc<[u32]>>,
}

impl UnreadableQuery {
    /// The Seq-id text of the Seq-align id of the subject `s_idx` (`CSeq_id::GetSeqIdString`):
    /// its parsed id; else `BL_ORD_ID:<oid>` of a database sequence, whose only id is the
    /// general one of its OID; else the local id `Subject_N` of a `-subject` record
    /// (`subject_skipped_id`).
    fn skipped_id(
        &self,
        parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
        s_idx: u32,
    ) -> Vec<u8> {
        subject_skipped_id(parsed, self.oids.as_deref(), s_idx)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1410-1590
// (`PrintOneResultSet`: formats 0 to 4 write the same query header, description table and
// footer; `s_SetFlags`, :633-649, gives formats 1 to 4 `eMergeAlign` in place of outfmt 0's
// alignments, cited at `report::query_anchored`.)
/// What a pairwise report (outfmt 0 to 4) writes after each query's description table: outfmt
/// 0's alignment of each HSP, or one query-anchored section of all of them
/// (`report::query_anchored`): 1 with identities (`.` for a letter equal to the query's) and only
/// the query's columns, 2 without identities, 3 and 4 likewise with every column ("flat").
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PairwiseView {
    #[default]
    Pairwise,
    QueryAnchored {
        /// `eShowIdentity` (formats 1 and 3).
        identity: bool,
        /// `eMasterAnchored` (formats 1 and 2).
        anchored: bool,
    },
}

impl PairwiseView {
    /// The view of a report format: outfmt 0 to 4; `None` for the others.
    pub fn of(format: crate::blastinput::app::ReportFormat) -> Option<Self> {
        use crate::blastinput::app::ReportFormat;
        Some(match format {
            ReportFormat::Pairwise => Self::Pairwise,
            ReportFormat::QueryAnchoredIdentities => Self::QueryAnchored {
                identity: true,
                anchored: true,
            },
            ReportFormat::QueryAnchoredNoIdentities => Self::QueryAnchored {
                identity: false,
                anchored: true,
            },
            ReportFormat::FlatQueryAnchoredIdentities => Self::QueryAnchored {
                identity: true,
                anchored: false,
            },
            ReportFormat::FlatQueryAnchoredNoIdentities => Self::QueryAnchored {
                identity: false,
                anchored: false,
            },
            _ => return None,
        })
    }

    /// The format's number (0 to 4).
    pub fn outfmt(self) -> i32 {
        match self {
            Self::Pairwise => 0,
            Self::QueryAnchored {
                identity: true,
                anchored: true,
            } => 1,
            Self::QueryAnchored {
                identity: false,
                anchored: true,
            } => 2,
            Self::QueryAnchored {
                identity: true,
                anchored: false,
            } => 3,
            Self::QueryAnchored {
                identity: false,
                anchored: false,
            } => 4,
        }
    }

    /// What the format's report does after a subject whose lookup throws where its batch's
    /// prefetch is off (`report::db_subject::AfterThrow`), with `-html` where `html`: outfmt 0
    /// skips the alignments of such a subject (or ends the run at its first alignment), 1 to 4 end
    /// the run at the first; with `-html` the run ends at its description line
    /// (`AfterThrow::description_end`). The one answer of the model of the loads
    /// (`db_subject::FormatThrow::of`) and of the writers (`PairwiseConfig::after_throw`).
    pub fn after_throw(self, html: bool) -> crate::report::db_subject::AfterThrow {
        use crate::report::db_subject::AfterThrow;
        match self {
            Self::Pairwise => AfterThrow::SkipsAlignment { html },
            Self::QueryAnchored { .. } => AfterThrow::EndsAtAlignment { html },
        }
    }

    /// How a report of this view reads its subjects (`PairwiseReading`): `sort` as it reads it
    /// (`report_sort`) and what it shows of each query (`shown`, `-num_descriptions` and
    /// `-num_alignments`).
    pub fn reading(
        self,
        sort: super::hit_sort::ReportSort,
        shown: crate::report::db_subject::Shown,
    ) -> PairwiseReading {
        PairwiseReading {
            sort: self.report_sort(sort),
            shown,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2980-2987
    /// `sort` as this format reads it: `-sorthits` for every view (formats 0 to 4), `-sorthsps`
    /// for outfmt 0 only.
    pub fn report_sort(self, sort: super::hit_sort::ReportSort) -> super::hit_sort::ReportSort {
        match self {
            Self::Pairwise => sort,
            Self::QueryAnchored { .. } => super::hit_sort::ReportSort {
                hits: sort.hits,
                hsps: None,
            },
        }
    }
}

/// How a pairwise report (outfmt 0 to 4) reads the subjects of each query, for the models of its
/// reads made before it is written (`x_pairwise_par::QueryGroup::read_order`): the order of its
/// sorted hits (`sort`, `PairwiseView::report_sort` of `-sorthits` and `-sorthsps`) and what it
/// shows of them (`shown`, `-num_descriptions` and `-num_alignments`), as its view's writer reads
/// them (`PairwiseView::reading`). Where `-html`'s description lines end the run is what the
/// report does after a throw (`PairwiseView::after_throw`, `AfterThrow::description_end`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PairwiseReading {
    pub sort: super::hit_sort::ReportSort,
    pub shown: crate::report::db_subject::Shown,
}

impl PairwiseConfig {
    /// The run's `-html` data of the report's prolog (`write_prolog`), `None` for the plain report.
    pub(crate) fn html_links(&self) -> Option<&super::html_ids::HtmlLinks> {
        self.html.as_ref().map(|report| &*report.links)
    }

    /// What the report does after a subject whose lookup throws (`PairwiseView::after_throw` of
    /// its view, with `-html` where it has its links).
    pub(crate) fn after_throw(&self) -> crate::report::db_subject::AfterThrow {
        self.view.after_throw(self.html.is_some())
    }

    /// The `-html` report with what it does after a throw (`after_throw`), for its description
    /// table (`HtmlRows`); `None` for the plain report.
    pub(crate) fn html_rows(
        &self,
    ) -> Option<(
        &super::html_ids::HtmlReport,
        crate::report::db_subject::AfterThrow,
    )> {
        self.html
            .as_ref()
            .map(|report| (report, self.after_throw()))
    }
}

/// `write_prolog`'s `html` of a report of `-html` links `links` (`None`: the plain report): the
/// links and the sink of the report's posts (`QueryWarnings::html_posts_sink`: the run's
/// diagnostics), which a `-html` report of outfmt 0 to 4 gets at any place among the run's
/// formats (`query_warnings::report_posts`, `ReportFormat::writes_html_links`). One without them
/// is LOSAT's internal error (its posts would be lost).
pub(crate) fn html_with_diagnostics<'w>(
    links: Option<&'w super::html_ids::HtmlLinks>,
    warnings: Option<&'w mut super::query_warnings::QueryWarnings<'_>>,
) -> io::Result<Option<(&'w super::html_ids::HtmlLinks, &'w mut dyn Write)>> {
    match (links, warnings) {
        (None, _) => Ok(None),
        (Some(links), Some(warnings)) => Ok(Some((links, warnings.html_posts_sink()))),
        (Some(_), None) => Err(io::Error::other(
            super::query_warnings::HTML_POSTS_WITHOUT_WARNINGS,
        )),
    }
}

impl PairwiseConfig {
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1525-1528
    /// Whether the report is TBLASTX's (`m_Program` is "tblastx"): the sorts and the description
    /// table's identities read its alignment lengths in nucleotides (`HspScores::alignment_length`).
    pub(crate) fn tblastx(&self) -> bool {
        self.program == "tblastx"
    }

    /// The columns that `-sorthits` adds to the description table of a query of `query_length`
    /// letters (`ReportSort::shows_columns`), or `None`.
    pub(crate) fn sort_columns(&self, query_length: usize) -> Option<SortColumns> {
        self.sort.shows_columns().then(|| SortColumns {
            query_length,
            tblastx: self.tblastx(),
        })
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:518-522
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1527-1530
/// The columns of a description table with `-sorthits` (NCBI's `eShowTotalScore`,
/// `eShowQueryCoverage` and `eShowPercentIdent`): the length of the query (`m_QueryLength`, the
/// query Bioseq's) and whether the report is TBLASTX's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct SortColumns {
    pub(crate) query_length: usize,
    pub(crate) tblastx: bool,
}

/// Extended hit information for pairwise display
///
/// For full pairwise output, we need the actual aligned sequences.
/// This struct extends Hit with optional sequence data.
#[derive(Debug, Clone)]
pub struct PairwiseHit {
    /// Base hit information
    pub hit: Hit,
    /// Aligned query sequence (if available)
    pub query_seq: Option<String>,
    /// Aligned subject sequence (if available)
    pub subject_seq: Option<String>,
    /// Query frame (for translated searches)
    pub query_frame: Option<i8>,
    /// Subject frame (for translated searches)
    pub subject_frame: Option<i8>,
    /// Number of positive matches (for protein/translated)
    pub positives: Option<usize>,
    /// Number of gaps
    pub gaps: Option<usize>,
    /// Subject sequence length
    pub subject_length: Option<usize>,
    /// Subject description/title
    pub subject_title: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:331-342;
    // c++/src/objtools/align_format/showalign.cpp:3599-3604:
    // eCompoScaleOldMatrix => Method 1; other adjusted matrices => Method 2.
    pub comp_adjust_method: Option<u8>,
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3595-3598:
    pub sum_n: Option<i32>,
    /// The columns (0-based) of `subject_seq` whose letter is a stored gap character (`-`,
    /// NCBIstdaa 0, which makeblastdb stores for a `-` of a protein FASTA): residues of the
    /// subject that the row shows as the letter of a gap. Empty when the subject holds none.
    /// NCBI counts the coordinates of the rows by the Dense-seg, where these columns are
    /// residues (`write_hsp_rows`).
    pub subject_stored_gaps: Vec<u32>,
    /// The columns (0-based) of `query_seq` whose letter is the gap residue of the query (`-` of
    /// ncbieaa, 0 of ncbistdaa, which a PSSM query of `-in_pssm` can hold): residues of the query
    /// that the row shows as the letter of a gap, as `subject_stored_gaps` are the subject's.
    /// Empty when the query holds none.
    pub query_stored_gaps: Vec<u32>,
    /// How the reports compare the letters of `query_seq` (`row_letters::QueryCase`): a record's
    /// row, whose lower case is a mask, or a PSSM query's Bioseq letters as they are.
    pub(crate) query_case: super::row_letters::QueryCase,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1199-1203
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:828-830
    /// The Seq-align's `num_ident` score of the HSP: the engine's identities (`BlastHSP::num_ident`),
    /// the value that the program's Seq-align writer (outfmt 8, 11, 12) writes, which the report
    /// reads back with `GetAlnScores` for the sorts of `-sorthits` and `-sorthsps` and the
    /// description table's "Max ident" (`hit_sort::HspScores`). It can differ from `hit.num_ident`,
    /// the identities of the displayed rows: TBLASTN's engine counts a SEG-masked query letter
    /// as a mismatch; BLASTP's display recount (pyrrolysine, stored gap characters), BLASTX's and
    /// TBLASTX's rows translated by the display's tables count others. Every count is at least
    /// 0, so every Seq-align carries it.
    pub(crate) seqalign_num_ident: usize,
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1950-1982
    // NCBI reference (598d8ae6): c++/include/objects/seqalign/Dense_diag.hpp:98-121
    /// Program-specific report edge marker packed into one byte so `PairwiseHit` stays at 376
    /// bytes: [`PairwiseEdge::QueryLeftOfRecord`] for BLASTN `-db` class N HSPs (carrying record
    /// query start/len in `hit.sort_query_offset`/`hit.sort_query_end`), or the two TBLASTN
    /// `-in_pssm -window_size 0` NULLB-start stop/skip markers. [`PairwiseEdge::Normal`] for all
    /// ordinary HSPs.
    pub(crate) edge: PairwiseEdge,
}

/// One-byte mutually exclusive report edge marker for [`PairwiseHit::edge`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum PairwiseEdge {
    /// Ordinary HSP.
    #[default]
    Normal,
    /// BLASTN `-db` HSP whose Seq-align query start lies below 0 in record coordinates (class N):
    /// outfmt 0 skips it at `x_PrepareRowData` without starting a heading or changing `previousId`,
    /// outfmt 1-4 feeds its negative query start into `CAlnMix` (`alnmerger.cpp:372`), and
    /// `HspScores` uses its unsigned `TSeqPos` query range for coverage and signed minimum
    /// coordinate for `-sorthsps 2`.
    QueryLeftOfRecord,
    /// TBLASTN `-in_pssm -window_size 0` `-subject`, positive frame, `rem in {1, 2}`, `k >= 2`:
    /// `CreateDensegFromStdseg` throws `Std-seg segment lengths not accurate!`
    /// (`EmptyRowsStop::StdSegLengths`).
    NullbStdSegLengths,
    /// TBLASTN `-in_pssm -window_size 0`, positive frame, `rem in {1, 2}`, `-db` or (`-subject`
    /// with `k == 1`): `CAlnVec::TranslateNAToAA` throws inside outfmt 0's per-HSP `try` block,
    /// writing `alignment skipped` with `subject_skipped_id` without emitting the heading yet.
    NullbSkipTranslate,
    /// TBLASTN `-in_pssm -window_size 0` `-db`, negative frame, every subject residue in the
    /// trailing `Sub` segment (`NullbReport::subject_row_empty`): `GetSeqRange` on row 1 throws
    /// `Row is empty` (`EmptyRowsStop::RowRange(RowRangeError::Empty)`) in `-sorthsps 4` and the
    /// description table; with `-num_descriptions 0`, outfmt 0 skips the HSP inside its per-HSP
    /// `try` block and outfmt 1-4 merges the query segment without a subject row.
    NullbEmptySubjectRow,
}

impl PairwiseHit {
    /// Whether this is a BLASTN `-db` class N HSP (`PairwiseEdge::QueryLeftOfRecord`).
    #[inline]
    pub(crate) fn query_left_of_record(&self) -> bool {
        self.edge == PairwiseEdge::QueryLeftOfRecord
    }

    /// Whether outfmt 0 skips this HSP inside NCBI's per-HSP `try` block before `x_PrintDefLine`
    /// (`PairwiseEdge::QueryLeftOfRecord` for BLASTN class N, or `PairwiseEdge::NullbSkipTranslate`
    /// / `PairwiseEdge::NullbEmptySubjectRow` for TBLASTN `-in_pssm -window_size 0`).
    #[inline]
    fn skips_pairwise_hsp(&self) -> bool {
        matches!(
            self.edge,
            PairwiseEdge::QueryLeftOfRecord
                | PairwiseEdge::NullbSkipTranslate
                | PairwiseEdge::NullbEmptySubjectRow
        )
    }

    /// Record-coordinate `(query_start, query_len)` of a BLASTN `-db` class N HSP
    /// (`query_left_of_record`), or `None` for normal in-context HSPs.
    #[inline]
    pub(crate) fn left_of_record(&self) -> Option<(i32, u32)> {
        self.query_left_of_record().then_some((
            self.hit.sort_query_offset as u32 as i32,
            self.hit.sort_query_end as u32,
        ))
    }
    /// Whether the rows of the HSP are empty: BLASTX's ungapped HSPs of length 0
    /// (`algorithm::blastx::empty_rows::rows_empty`), the only ones without a column, whose
    /// `CSeq_align::GetSeqRange` throws `Row is empty` (Seq_align.cpp:227). It reads the hit's
    /// offsets, which a BLASTX hit takes from its Seq-align (`Hsp::seqalign_ends`).
    pub(crate) fn rows_empty(&self) -> bool {
        self.hit.has_sort_offsets
            && super::serial::seqalign::ungapped_rows_empty(
                self.hit.sort_query_offset as i64,
                self.hit.sort_query_end as i64,
            )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1100-1122
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2576-2613
    /// The HSP's columns, gap letters and gap openings inside the query row's aligned range
    /// (`serial::seqalign::AlignLengths::within_query_range`: the query row's end gaps of a
    /// Smith-Waterman alignment left out), as NCBI's tables count them (`GetAlignLengths`) and
    /// outfmt 0's HSP header its `Gaps` (`x_GetNumGaps`, in `write_hsp_identities` and
    /// `write_tblastn_hsp_info`; the line's denominators stay the alignment's columns). The
    /// library's `HspRecord` holds them (`api::run::records`).
    pub(crate) fn query_range_lengths(&self) -> super::serial::seqalign::AlignLengths {
        super::serial::seqalign::AlignLengths {
            length: self.hit.length,
            gaps: self.gaps.unwrap_or_else(|| self.hit.gap_letters()),
            gap_opens: self.hit.gapopen,
        }
        .within_query_range(self.hit.gap_info.as_deref().unwrap_or_default())
    }
}

/// What `CSeq_align::GetSeqRange` of the row finds on the Seq-align of `hit`, of the subject whose hits are
/// `hits`, as NCBI's formatters read it (`serial::seqalign::reported_std_seg_row_range`):
/// `-subject`'s Seq-loc mapper (`remaps`) rewrites the empty query row of a subject with rows of
/// letters into one nucleotide (`mapper_rewrites`; the empty subject row stays). Its leaving out of
/// an end segment with a gap row is in the hit already (`serial::seqalign::mapped_script`). The
/// hits that keep their Seq-align's offsets (`has_sort_offsets`: BLASTX, TBLASTN and TBLASTX) have
/// Std-segs; the others' Dense-seg has a range in every row.
fn seen_row_range(
    remaps: bool,
    hits: &[&PairwiseHit],
    hit: &PairwiseHit,
    row: AlignRow,
) -> Option<RowRangeError> {
    if !hit.hit.has_sort_offsets {
        return None;
    }
    // NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:176-230
    // A TBLASTN `-db` negative-frame NULLB-start HSP with no backed subject residues has an empty
    // subject range (`Row is empty`), even when `rem in {0, 1}` keeps its trailing `Sub` op.
    if row == AlignRow::Subject && hit.edge == PairwiseEdge::NullbEmptySubjectRow {
        return Some(RowRangeError::Empty);
    }
    super::serial::seqalign::reported_std_seg_row_range(
        hit.hit.sort_query_offset as i64,
        hit.hit.sort_query_end as i64,
        hit.hit.gap_info.as_deref().unwrap_or_default(),
        row,
        mapper_rewrites(remaps, hits),
    )
}

/// The first exception of `GetSeqRange` where a formatter reads the ranges of `row` of every HSP of
/// the subjects `oids` (their hits in `grouped`), in their order (`seen_row_range`).
fn read_row_range<'a>(
    remaps: bool,
    mut oids: impl Iterator<Item = &'a u32>,
    grouped: &std::collections::HashMap<u32, Vec<&PairwiseHit>>,
    row: AlignRow,
) -> Option<RowRangeError> {
    oids.find_map(|oid| {
        let hits = &grouped[oid];
        hits.iter()
            .find_map(|hit| seen_row_range(remaps, hits, hit, row))
    })
}

/// Where `write_blastx_pairwise_report` stopped a query's report at an HSP with empty rows
/// (`PairwiseHit::rows_empty`) or a row that starts with a gap (`PairwiseHit::row_range`), as the
/// error of its writes: which of NCBI's exceptions ends the report there
/// (`algorithm::blastx::empty_rows`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum EmptyRowsStop {
    /// `CSeq_align::GetSeqRange`'s exception (`Row is empty`, Seq_align.cpp:227; `Row seqids not
    /// consistent`, Seq_align.cpp:186): the description table, or a sort's comparator, reads a
    /// row whose range it cannot take (after the query's `Length=`).
    RowRange(RowRangeError),
    /// `CreateDensegFromStdseg(): Std-seg segment lengths not accurate!` (Seq_align.cpp:801): a
    /// sort's comparator or the alignment section makes a Dense-seg of `-subject`'s mapped
    /// Std-seg (before anything of the HSP), also the query-anchored section of outfmt 1 to 4
    /// (before its first byte).
    StdSegLengths,
    /// Not an exception: the merge of a query-anchored section (outfmt 1 to 4) frees a segment
    /// that it still reads (`blastx_anchored_stop`), which LOSAT refuses before the call's first
    /// byte (`UNPORTED_ANCHORED_FREE`).
    MergeFreesSegment,
}

/// The error of outfmt 0's writer (and of outfmt 1 to 4's): its stop at an HSP of length 0 or at a
/// row that starts with a gap (`EmptyRowsStop`) is NCBI's exception, `CSeq_align::GetSeqRange`'s
/// where the sorts or the description table read such a row (after the query's `Length=`;
/// `serial::seqalign::row_range_exception`), `CreateDensegFromStdseg` where they or
/// `-subject`'s alignment section make a Dense-seg of the Seq-loc mapper's rows, or LOSAT's
/// refusal of a query-anchored section whose merge NCBI frees a segment in (before the call's
/// first byte); the exception of a subject's load that ends the run (`cli::NativeError::of_writer`,
/// `x_pairwise_par::alignment_end`); a write error stays. The one map of the writers' stops: BLASTX's
/// and TBLASTN's outfmt 0 to 4 map their writer's errors with it.
pub(crate) fn empty_rows_stop(error: std::io::Error) -> anyhow::Error {
    match error
        .get_ref()
        .and_then(|e| e.downcast_ref::<EmptyRowsStop>())
    {
        Some(EmptyRowsStop::RowRange(error)) => crate::cli::NativeError::exception(
            &super::serial::seqalign::row_range_exception(*error),
        )
        .into(),
        Some(EmptyRowsStop::StdSegLengths) => crate::cli::NativeError::exception(
            &super::serial::seqalign::std_seg_lengths_not_accurate(),
        )
        .into(),
        Some(EmptyRowsStop::MergeFreesSegment) => {
            EmptyRowsRefusal(UNPORTED_ANCHORED_FREE.to_string()).into()
        }
        None => crate::cli::NativeError::of_writer(error),
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

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmerger.cpp:598-606 (cited at
// `report::query_anchored::aln_mix::Mix::place`), c++/src/objtools/alnmgr/alnsegments.cpp:
// 310-327 (cited at `report::query_anchored::aln_mix::Mix::fill_unaligned_regions`)
/// The refusal of a query-anchored report (outfmt 1 to 4) whose section NCBI's merge frees a
/// segment in (`EmptyRowsStop::MergeFreesSegment`,
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

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1522-1530
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3074-3094
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4334-4335
/// The exception at which NCBI's sorts of outfmt 0 (`-sorthits`, `-sorthsps` above 0, before the
/// description table) stop a query's report whose hit list (`order`, the subjects' hits in
/// `grouped`) holds an HSP with a row whose range they cannot take (`seen_row_range`: empty rows,
/// a row that starts with a gap); `remaps`: `-subject`, whose Seq-loc mapper gives an HSP of length
/// 0 a query row of one nucleotide in a subject with rows of letters (`mapper_rewrites`; the rows
/// of a subject of empty rows only stay as `-db`'s). `std::list::sort` compares every element of a
/// list of two or more, so a comparator that reads an HSP's ranges reaches it:
/// - the subjects by percent identity (`-sorthits 3`: `GetSeqAlignSetCalcPercentIdent` makes a
///   Dense-seg of every HSP, `GetAlignmentLength`): `-subject`'s mapped Std-seg throws (801), the
///   raw empty rows of `-db` give a length of 0;
/// - the subjects by query coverage (`-sorthits 4`: `GetMasterCoverage`): the query row (`-db`'s
///   empty one, 227; one that starts with a gap, 186), `-subject`'s nucleotide is not empty;
/// - then each subject's HSPs, a list of two or more: by query start (`-sorthsps 2`,
///   `GetSeqStart` of row 0) the query row as `-sorthits 4`; by percent identity (`-sorthsps 3`)
///   `-subject`'s Dense-seg (801); by subject start (`-sorthsps 4`, `GetSeqStart` of row 1) the
///   subject row (the empty one of both, 227; one that starts with a gap, 186).
///
/// The scores' sorts (`-sorthits 1` and `2`, `-sorthsps 1`) read no range. Probes: survey set K
/// (`X/SURVEY_X2g_formats.md` 6), and line PS item PS12d's (`PR/ps12d/scripts/gcmp.sh`).
fn sort_stop(
    sort: super::hit_sort::ReportSort,
    remaps: bool,
    order: &[u32],
    grouped: &std::collections::HashMap<u32, Vec<&PairwiseHit>>,
) -> Option<EmptyRowsStop> {
    // Whether the subject holds an HSP with empty rows that the mapper rewrote: its Dense-seg
    // throws.
    let mapped_empty = |oid: &u32| {
        let hits = &grouped[oid];
        hits.iter().any(|hit| hit.rows_empty()) && mapper_rewrites(remaps, hits)
    };
    if order.len() >= 2 {
        match sort.hits {
            Some(3) if order.iter().any(mapped_empty) => return Some(EmptyRowsStop::StdSegLengths),
            Some(4) => {
                if let Some(error) = read_row_range(remaps, order.iter(), grouped, AlignRow::Query)
                {
                    return Some(EmptyRowsStop::RowRange(error));
                }
            }
            _ => {}
        }
    }
    let listed = || order.iter().filter(|oid| grouped[*oid].len() >= 2);
    let row = match sort.hsps {
        Some(2) => AlignRow::Query,
        Some(3) if listed().any(mapped_empty) => return Some(EmptyRowsStop::StdSegLengths),
        Some(4) => AlignRow::Subject,
        _ => return None,
    };
    read_row_range(remaps, listed(), grouped, row).map(EmptyRowsStop::RowRange)
}

/// Whether `-subject`'s Seq-loc mapper (`remaps`) rewrote the empty rows of the subject whose
/// hits are `hits` (`serial::seqalign::mapper_rewrites_empty_rows`, whose frames make an empty
/// query row one nucleotide, `x_GetDstStd`, seq_align_mapper_base.cpp:1419-1446); otherwise they
/// stay as `-db`'s.
fn mapper_rewrites(remaps: bool, hits: &[&PairwiseHit]) -> bool {
    super::serial::seqalign::mapper_rewrites_empty_rows(
        remaps,
        hits.iter().map(|hit| hit.rows_empty()),
    )
}

impl std::fmt::Display for EmptyRowsStop {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for EmptyRowsStop {}

impl From<Hit> for PairwiseHit {
    fn from(hit: Hit) -> Self {
        let positives = hit.num_positives;
        let gaps = hit.gap_letters();
        let seqalign_num_ident = hit.num_ident;
        Self {
            hit,
            query_seq: None,
            subject_seq: None,
            query_frame: None,
            subject_frame: None,
            positives: Some(positives),
            gaps: Some(gaps),
            subject_length: None,
            subject_title: None,
            comp_adjust_method: None,
            sum_n: None,
            subject_stored_gaps: Vec::new(),
            query_stored_gaps: Vec::new(),
            query_case: Default::default(),
            seqalign_num_ident,
            edge: PairwiseEdge::Normal,
        }
    }
}

/// NCBI BLASTP pairwise per-query ancillary data.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:72-115
#[derive(Debug, Clone)]
pub struct BlastpPairwiseQuery {
    /// BLASTX's `Query=` text (the others read their queries' title bytes by query index).
    pub query_name: Vec<u8>,
    pub query_length: usize,
    pub ungapped_karlin: KarlinParams,
    pub effective_search_space: i64,
    /// Whether the query's Karlin-Altschul parameters could be computed (`is_valid`).
    pub valid: bool,
    /// Whether every query of the query's batch is invalid, so NCBI did not search the
    /// batch (local_blast.cpp:177-207: -1 parameters and no Seq-align set).
    pub batch_skipped: bool,
}

/// NCBI BLASTP pairwise run-level footer/header data.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2248-2285
#[derive(Debug, Clone)]
pub struct BlastpPairwiseReport {
    pub version: String,
    pub database_name: Vec<u8>,
    /// The `Database:` lines' counts of a FASTA subject set (its records and letters), NCBI's
    /// `Int8` (`tot_num_seqs`, `tot_length`); a `-db` search prints `PairwiseDatabase` and keeps
    /// its sums here (an alias file's `LENGTH` can pass 2^32 on a 32-bit target).
    pub database_num_sequences: i64,
    pub database_total_letters: i64,
    /// `options.GetMatrixName()`; `None` where NCBI's options name no matrix (a null `const
    /// char*`: an RPS search that set up no query batch, `write_matrix_line`).
    pub matrix_name: Option<String>,
    pub gap_open: i32,
    pub gap_extend: i32,
    /// `GetWordThreshold()`, a double.
    pub word_threshold: f64,
    pub window_size: i32,
    pub gapped_karlin: KarlinParams,
    /// None with `-searchsp` (`PrintKAParameters` without `gbp`).
    pub gumbel: Option<BlastGumbelBlk>,
    /// The options' gapped mode: an ungapped search's footers print the ungapped Karlin block
    /// only and its epilog no `Gap Penalties` line (`gapped_karlin` is then unused).
    pub gapped: bool,
    /// NCBI's `m_ShowLinkedSetSize` (sum statistics in an ungapped search): the description
    /// table's `N` column.
    pub show_linked_set_size: bool,
    /// Subjects in the description table and with alignments, per query (NCBI's
    /// `m_NumDescriptions` and `m_NumAlignments`).
    pub num_descriptions: usize,
    pub num_alignments: usize,
}

/// One database of a `-db` list as the reports show it (NCBI's `SDbInfo` of a BLAST
/// database, `seqdb::NameRecord`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:385-391
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DatabaseInfo {
    /// `definition`: the title of the database, or its name for an empty title.
    pub definition: Vec<u8>,
    /// `date`.
    pub date: String,
    /// `number_seqs`: NCBI's `int` of `GetNumSeqs` (negative for an alias file's negative
    /// `NSEQ`, or one of 2^31 or more).
    pub num_sequences: i32,
    /// `total_length`: NCBI's `Int8` of `GetTotalLength` (negative for an alias file's `LENGTH`
    /// of 2^63 or more).
    pub total_letters: i64,
    /// `filt_algorithm_name`: the database masking's algorithm (`Masked using:`), empty for a
    /// search without it.
    pub filt_algorithm_name: String,
    /// `filt_algorithm_options`: its options, empty for its defaults.
    pub filt_algorithm_options: String,
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:511-523
/// The database of a `-db` search as the pairwise report shows it: one `DatabaseInfo` per
/// name of the list (NCBI's `m_DbInfo`), and the values that the `Database:` line of the
/// prolog, the database counts of XML and JSON and the title of the Seq-annot formats take
/// of them all. A search of a FASTA subject set has none: its `Database:` lines are
/// `BlastpPairwiseReport::database_name` with a period, and its date is `Unknown`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairwiseDatabase {
    /// The definitions of the databases joined with `"; "`.
    pub definition: Vec<u8>,
    /// The sum of their numbers of sequences, NCBI's `Int8` (`tot_num_seqs`).
    pub num_sequences: i64,
    /// The sum of their total lengths, NCBI's `Int8` (`tot_length`). A sum that overflows the
    /// `Int8` (`total_letters_overflow`) is refused where NCBI's formatter would sum it
    /// (`common::database::report_sum_refusal`).
    pub total_letters: i64,
    /// Whether the sum of the total lengths passes NCBI's `Int8` (a signed overflow).
    pub total_letters_overflow: bool,
    /// The databases in the order of the list, a name given twice twice: the blocks of the
    /// epilog.
    pub databases: Vec<DatabaseInfo>,
}

impl PairwiseDatabase {
    /// The report's database for the databases `databases` of the list.
    pub fn new(databases: Vec<DatabaseInfo>) -> Self {
        let mut definition = Vec::new();
        for (at, database) in databases.iter().enumerate() {
            if at > 0 {
                definition.extend_from_slice(b"; ");
            }
            definition.extend_from_slice(&database.definition);
        }
        Self {
            definition,
            num_sequences: databases
                .iter()
                .map(|database| i64::from(database.num_sequences))
                .fold(0, i64::wrapping_add),
            total_letters: databases
                .iter()
                .map(|database| database.total_letters)
                .fold(0, i64::wrapping_add),
            total_letters_overflow: databases
                .iter()
                .try_fold(0i64, |sum, database| {
                    sum.checked_add(database.total_letters)
                })
                .is_none(),
            databases,
        }
    }
}

// =============================================================================
// Pairwise Output Writers
// =============================================================================

/// Write database/subject information header
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp
///
/// Format:
/// ```text
/// >subject_id subject_description
/// Length=XXXX
/// ```
pub fn write_subject_header<W: Write>(
    writer: &mut W,
    subject_id: &[u8],
    subject_title: Option<&[u8]>,
    subject_length: Option<usize>,
) -> io::Result<()> {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2385-2388,2461-2471,340-359
    // `s_WrapOutputLine` counts and writes bytes (`str[i]`, `isspace` of `str[i]` cast to
    // `unsigned char`), so a line can end inside a UTF-8 sequence.
    let mut title = subject_id.to_vec();
    if let Some(subject_title) = subject_title.filter(|title| !title.is_empty()) {
        title.push(b' ');
        title.extend_from_slice(subject_title);
    }
    writer.write_all(b"> ")?;
    write_wrapped(writer, &title)?;
    writeln!(writer)?;

    // Length line
    if let Some(len) = subject_length {
        writeln!(writer, "Length={}", len)?;
    }

    writeln!(writer)?;
    Ok(())
}

/// The alignments of one subject in outfmt 0: NCBI's `x_DisplayAlnvecInfo` for each HSP,
/// the first one with the subject's heading (`GenerateDefline` of its title with flags 0,
/// `defline.rs`). `write_hsp` writes the rest of one call (the score lines, the alignment
/// and what follows it) of an HSP and its index in the report's hits (`HspIndex`). `probe`
/// marks the heading and each HSP without changing the written bytes.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3613-3632
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2273
/// Where `GenerateDefline` throws (a title whose encoding NCBI cannot guess,
/// `UnknownEncoding`), NCBI catches the exception for that HSP and writes one line in its
/// place. `previousId` is not set then, so the next HSP of the subject tries the heading
/// again: every HSP of the subject gets the line. `GetSeqIdString` of the local ID
/// `lcl|Subject_N` is `Subject_N`, N being the subject's position in the input, which is
/// `s_idx + 1` (every subject record keeps its index). The probe gets an empty heading.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1950-1981
/// The alignments of the HSPs `part` of the subject's (the heading with the first); the parts of
/// `0..hits.len()` in order write its bytes (`x_pairwise_par`'s pieces).
///
/// With `-html` (`html`) the heading is `>`, the link of `lcl|Subject_N` (the registry's `[LOCAL_ID]
/// TOOL_URL`, if any) and `</a>` after one, `<a name=Subject_N></a>`, a space and the
/// HTML-encoded title wrapped (`HtmlHeading::unparsed_start`; the link precedes the title,
/// `x_FillAlnDispParams` with the Bioseq's handle), and where the registry throws NCBI's line is
/// written in place of each alignment too (`HtmlHeading::skips`).
#[allow(clippy::too_many_arguments)]
fn write_subject_alignments_part<W: Write>(
    writer: &mut W,
    mut probe: Option<&mut FormatProbe<'_>>,
    s_idx: u32,
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    title: &[u8],
    protein: bool,
    hits: &[&PairwiseHit],
    hsp_indices: &[HspIndex],
    part: std::ops::Range<usize>,
    html: Option<HtmlHeading<'_>>,
    place: super::query_warnings::ReportHook<'_, W>,
    write_hsp: impl FnMut(&mut W, &PairwiseHit, HspIndex) -> io::Result<()>,
) -> io::Result<()> {
    if hits.iter().any(|hit| hit.skips_pairwise_hsp()) {
        return write_subject_with_skips(
            writer,
            probe,
            SubjectHeading::Record(title),
            protein,
            parsed,
            None,
            s_idx,
            hits,
            hsp_indices,
            part,
            html,
            place,
            write_hsp,
        );
    }
    let first_index = hsp_indices[0];
    let heading = super::defline::generate_defline(title, protein, false);
    let attempt = |html: HtmlHeading<'_>| record_html_lines(&heading, html, parsed, s_idx);
    let skipped = match html {
        Some(html) if part.start == 0 => {
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.subject_begin(first_index);
            }
            let attempted = attempt(html);
            place_heading_posts(html, true, writer, &mut *place)?;
            let skipped = match attempted? {
                Ok(lines) => {
                    write_heading_lines(writer, &lines, hits[0].subject_length)?;
                    None
                }
                Err(id) => Some(id),
            };
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.subject_end(first_index);
            }
            skipped
        }
        Some(html) => {
            let skipped = attempt(html).map(Result::err);
            later_part_posts(html, &skipped, writer, &mut *place)?;
            skipped?
        }
        None => {
            if part.start == 0 {
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.subject_begin(first_index);
                }
                if let Ok(heading) = &heading {
                    write_subject_header(writer, &heading.text, None, hits[0].subject_length)?;
                }
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.subject_end(first_index);
                }
            }
            heading
                .as_ref()
                .err()
                .map(|_| subject_skipped_id(parsed, None, s_idx))
        }
    };
    let retry = html.map(|html| {
        let retry = move || -> io::Result<Option<Vec<u8>>> { Ok(attempt(html)?.err()) };
        (retry, html)
    });
    write_hsps_or_skips(
        writer,
        probe,
        hits,
        hsp_indices,
        part,
        skipped,
        retry,
        place,
        write_hsp,
    )
}

/// Writes at `place` the reads of NCBI's formatter registries that the heading links of `html`
/// made (`HtmlHeading::take_posts`, `ReportPlace::Html`), before the heading's bytes or NCBI's
/// line in place of its alignment; `made` false drops them (a heading made by a later piece of
/// the subject's alignments only to tell whether they are skipped, `later_part_posts`).
fn place_heading_posts<W: Write>(
    html: HtmlHeading<'_>,
    made: bool,
    writer: &mut W,
    place: super::query_warnings::ReportHook<'_, W>,
) -> io::Result<()> {
    for post in html.take_posts() {
        if made {
            place(super::query_warnings::ReportPlace::Html(post), writer)?;
        }
    }
    Ok(())
}

/// The reads of the registries of the heading that a later part of a subject's alignments
/// (EXPERIMENT LOSAT_X_PAIRWISEPAR's pieces) makes to tell whether they are skipped (`skipped`,
/// or its error): the serial writer makes that heading again only for a skipped subject
/// (`write_hsps_or_skips`' retry of the part's first HSP), so its reads are written only then, or
/// before an error that ends the report.
fn later_part_posts<W: Write>(
    html: HtmlHeading<'_>,
    skipped: &io::Result<Option<Vec<u8>>>,
    writer: &mut W,
    place: super::query_warnings::ReportHook<'_, W>,
) -> io::Result<()> {
    let made = !matches!(skipped, Ok(None));
    place_heading_posts(html, made, writer, place)
}

/// HTML heading lines of a `-subject` record without parsed Seq-ids, or the Seq-id text of NCBI's
/// `alignment skipped` line when the heading's title or HTML link throws.
fn record_html_lines(
    heading: &Result<super::defline::Defline, super::defline::UnknownEncoding>,
    html: HtmlHeading<'_>,
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    s_idx: u32,
) -> io::Result<Result<Vec<u8>, Vec<u8>>> {
    let mut lines = Vec::new();
    if let Err(error) = html.unparsed_start(&mut lines) {
        return HtmlHeading::skips(error, || html.unparsed_skipped_id()).map(Err);
    }
    Ok(match heading {
        Ok(heading) => {
            super::db_subject::wrap_title(&super::html::html_encode(&heading.text), &mut lines);
            lines.push(b'\n');
            Ok(lines)
        }
        Err(super::defline::UnknownEncoding) => Err(subject_skipped_id(parsed, None, s_idx)),
    })
}

type ParsedSubjects<'a> = &'a [Option<Arc<crate::report::db_subject::ParsedSubject>>];

/// The alignment heading of one subject of a report: a database sequence's or a `-subject`
/// record's with parsed Seq-ids (`DatabaseHeading`), or a `-subject` record's without them, with
/// its title (`write_subject_alignments_part`).
#[derive(Clone, Copy)]
enum SubjectHeading<'a> {
    Database(DatabaseHeading<'a>),
    Record(&'a [u8]),
}

impl SubjectHeading<'_> {
    /// How many headings NCBI starts for the subject (`QueryGroup::heading_rank`): none where its
    /// load throws (`DatabaseHeading::Skipped`), one for each alignment where the heading's title
    /// cannot be made (NCBI's line in place of each alignment), else one.
    fn starts(self, protein: bool) -> super::x_pairwise_par::HeadingStarts {
        use super::x_pairwise_par::HeadingStarts;
        let throws = match self {
            SubjectHeading::Database(DatabaseHeading::Skipped(_)) => return HeadingStarts::Never,
            SubjectHeading::Database(heading) => {
                matches!(heading.skipped_id(protein, None), Ok(Some(_)))
            }
            SubjectHeading::Record(title) => {
                super::defline::generate_defline(title, protein, false).is_err()
            }
        };
        if throws {
            HeadingStarts::EachAlignment
        } else {
            HeadingStarts::Once
        }
    }

    /// Seq-id text written by `write_skipped_alignment` for a per-HSP skip (`skips_pairwise_hsp`).
    fn skip_id(self, parsed: ParsedSubjects<'_>, oids: Option<&[u32]>, s_idx: u32) -> Vec<u8> {
        match self {
            SubjectHeading::Database(DatabaseHeading::Skipped(id)) => id.to_vec(),
            SubjectHeading::Database(DatabaseHeading::Parsed(subject, _)) => {
                subject.skipped_id().to_vec()
            }
            SubjectHeading::Database(DatabaseHeading::Title(..)) => {
                subject_skipped_id(parsed, oids, s_idx)
            }
            SubjectHeading::Record(_) => subject_skipped_id(parsed, None, s_idx),
        }
    }

    /// Writes the subject heading (`Database` or `Record`), returning `Some(skip_id)` when the
    /// heading itself throws and each alignment is skipped instead.
    #[allow(clippy::too_many_arguments)]
    fn write_heading<W: Write>(
        self,
        writer: &mut W,
        protein: bool,
        subject_length: Option<usize>,
        html: Option<HtmlHeading<'_>>,
        parsed: ParsedSubjects<'_>,
        s_idx: u32,
        place: super::query_warnings::ReportHook<'_, W>,
    ) -> io::Result<Option<Vec<u8>>> {
        match self {
            SubjectHeading::Database(db_heading) => write_database_subject_heading(
                writer,
                db_heading,
                protein,
                subject_length,
                html,
                place,
            ),
            SubjectHeading::Record(title) => {
                let defline = super::defline::generate_defline(title, protein, false);
                match html {
                    Some(html) => {
                        let attempted = record_html_lines(&defline, html, parsed, s_idx);
                        place_heading_posts(html, true, writer, place)?;
                        match attempted? {
                            Ok(lines) => {
                                write_heading_lines(writer, &lines, subject_length)?;
                                Ok(None)
                            }
                            Err(id) => Ok(Some(id)),
                        }
                    }
                    None => match defline {
                        Ok(defline) => {
                            write_subject_header(writer, &defline.text, None, subject_length)?;
                            Ok(None)
                        }
                        Err(super::defline::UnknownEncoding) => {
                            Ok(Some(subject_skipped_id(parsed, None, s_idx)))
                        }
                    },
                }
            }
        }
    }

    /// Checks whether the heading throws (for a later split part or an `-html` retry).
    fn skipped_id(
        self,
        protein: bool,
        html: Option<HtmlHeading<'_>>,
        parsed: ParsedSubjects<'_>,
        s_idx: u32,
    ) -> io::Result<Option<Vec<u8>>> {
        match self {
            SubjectHeading::Database(db_heading) => db_heading.skipped_id(protein, html),
            SubjectHeading::Record(title) => {
                let generated = super::defline::generate_defline(title, protein, false);
                match html {
                    Some(html) => {
                        record_html_lines(&generated, html, parsed, s_idx).map(Result::err)
                    }
                    None => Ok(generated
                        .err()
                        .map(|_| subject_skipped_id(parsed, None, s_idx))),
                }
            }
        }
    }
}

/// The heading of subject `s_idx` of `group` in a report whose subjects' parsed Seq-ids are
/// `parsed` (by subject index): the definition lines of one with them (`DatabaseHeading::loaded`
/// with `gapped`, `Skipped` where its load throws: `thrown`, `QuerySections::thrown`), else a
/// database sequence's stored title (`database`, the sources `sources`) or a `-subject` record's
/// title.
fn group_heading<'g>(
    group: &'g super::x_pairwise_par::QueryGroup<'_>,
    parsed: &'g [Option<Arc<crate::report::db_subject::ParsedSubject>>],
    sources: &'g [super::defline::SourceTaxname],
    database: bool,
    gapped: bool,
    thrown: bool,
    s_idx: u32,
) -> SubjectHeading<'g> {
    let title = &group.subject_titles[&s_idx];
    match parsed_subject(parsed, s_idx) {
        Some(subject) => SubjectHeading::Database(DatabaseHeading::loaded(
            subject,
            gapped,
            group.own_hsps(s_idx),
            thrown,
        )),
        None if database => {
            SubjectHeading::Database(DatabaseHeading::Title(title, stored_source(sources, s_idx)))
        }
        None => SubjectHeading::Record(title),
    }
}

/// The alignments of the HSPs `part` of subject `s_idx` of `group` (its heading with the first):
/// the subject's heading is `heading_of(s_idx)` (a database heading or a `-subject` record's,
/// `write_database_subject_alignments_part`, `write_subject_alignments_part`; `records`: the
/// `-subject` records' parsed Seq-ids), each HSP by `write_hsp`. With `-html` (`html`, the
/// report's links) the heading's rank counts the headings of the subjects before it
/// (`QueryGroup::heading_rank` with `heading_of`), made only then.
#[allow(clippy::too_many_arguments)]
fn write_subject_part<'h, W: Write>(
    writer: &mut W,
    probe: Option<&mut FormatProbe<'_>>,
    group: &super::x_pairwise_par::QueryGroup<'_>,
    s_idx: u32,
    part: std::ops::Range<usize>,
    protein: bool,
    html: Option<&super::html_ids::HtmlReport>,
    records: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    oids: Option<&[u32]>,
    heading_of: impl Fn(u32) -> SubjectHeading<'h>,
    place: super::query_warnings::ReportHook<'_, W>,
    write_hsp: impl FnMut(&mut W, &PairwiseHit, HspIndex) -> io::Result<()>,
) -> io::Result<()> {
    let posts = std::cell::RefCell::new(Vec::new());
    let html = html.map(|report| HtmlHeading {
        report,
        s_idx,
        rank: group.heading_rank(s_idx, |subject| heading_of(subject).starts(protein)),
        posts: &posts,
    });
    let (hits, indices) = (
        &group.subject_hits[&s_idx],
        &group.subject_hit_indices[&s_idx],
    );
    match heading_of(s_idx) {
        SubjectHeading::Database(heading) => write_database_subject_alignments_part(
            writer, probe, heading, protein, records, oids, s_idx, hits, indices, part, html,
            place, write_hsp,
        ),
        SubjectHeading::Record(title) => write_subject_alignments_part(
            writer, probe, s_idx, records, title, protein, hits, indices, part, html, place,
            write_hsp,
        ),
    }
}

/// One subject's alignments in a BLASTP or PSIBLAST pairwise report (`BlastpSections::subject`
/// and `psiblast_round::RoundSections::subject`): an unreadable query's skipped lines, or the
/// subject's heading (`group_heading`) and HSPs (`write_hsp_info` + `write_alignment`).
#[allow(clippy::too_many_arguments)]
fn write_blastp_subject_part<W: crate::report::ostream::ReportWrite>(
    writer: &mut W,
    probe: Option<&mut FormatProbe<'_>>,
    config: &PairwiseConfig,
    database: bool,
    parsed_subjects: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    subject_sources: &[super::defline::SourceTaxname],
    loads: &crate::report::db_subject::HitLoads,
    group: &super::x_pairwise_par::QueryGroup<'_>,
    s_idx: u32,
    part: std::ops::Range<usize>,
    row_masks: HtmlRowMasks<'_>,
    place: super::query_warnings::ReportHook<'_, W>,
) -> io::Result<()> {
    if let Some(id) = config.unreadable_query_id(parsed_subjects, s_idx) {
        return write_database_subject_alignments_part(
            writer,
            probe,
            DatabaseHeading::Skipped(&id),
            true,
            parsed_subjects,
            None,
            s_idx,
            &group.subject_hits[&s_idx],
            &group.subject_hit_indices[&s_idx],
            part,
            None,
            place,
            |_, _, _| unreachable!("a skipped alignment is written by its line"),
        );
    }
    write_subject_part(
        writer,
        probe,
        group,
        s_idx,
        part,
        true,
        config.html.as_ref(),
        parsed_subjects,
        None,
        |subject| {
            group_heading(
                group,
                parsed_subjects,
                subject_sources,
                database,
                true,
                thrown_load(loads, parsed_subjects, subject).is_some(),
                subject,
            )
        },
        place,
        |writer, hit, hsp_index| {
            write_hsp_info(writer, hit, config)?;
            write_alignment(writer, hit, row_masks.of(hsp_index), config)
        },
    )
}

/// The heading `lines` of a subject (each line ending with a newline), `Length=` when the length
/// is known, and one blank line.
fn write_heading_lines<W: Write>(
    writer: &mut W,
    lines: &[u8],
    length: Option<usize>,
) -> io::Result<()> {
    writer.write_all(lines)?;
    if let Some(length) = length {
        writeln!(writer, "Length={length}")?;
    }
    writeln!(writer)
}

/// The HSPs `part` of a subject (`hits`, `hsp_indices`) after its heading: each by `write_hsp`,
/// or NCBI's line in place of each where the heading was not made (`skipped`, the Seq-id text of
/// the line; showalign.cpp:1976-1981, cited at `write_subject_alignments_part`). NCBI sets
/// `previousId` only after a heading, so each alignment of a skipped subject makes its heading
/// again: with `-html` the next one's attempt is `retry` with the heading's `-html` data (its id
/// links, whose reads of the registries `place` writes before the line, `place_heading_posts`).
#[allow(clippy::too_many_arguments)]
fn write_hsps_or_skips<W: Write>(
    writer: &mut W,
    mut probe: Option<&mut FormatProbe<'_>>,
    hits: &[&PairwiseHit],
    hsp_indices: &[HspIndex],
    part: std::ops::Range<usize>,
    mut skipped: Option<Vec<u8>>,
    mut retry: Option<(impl FnMut() -> io::Result<Option<Vec<u8>>>, HtmlHeading<'_>)>,
    place: super::query_warnings::ReportHook<'_, W>,
    mut write_hsp: impl FnMut(&mut W, &PairwiseHit, HspIndex) -> io::Result<()>,
) -> io::Result<()> {
    for (at, (hit, &hsp_index)) in hits[part.clone()]
        .iter()
        .zip(&hsp_indices[part])
        .enumerate()
    {
        if let Some(probe) = probe.as_mut() {
            writer.flush()?;
            probe.begin(hsp_index);
        }
        if at > 0 && skipped.is_some() {
            if let Some((retry, html)) = retry.as_mut() {
                let retried = retry();
                place_heading_posts(*html, true, writer, &mut *place)?;
                skipped = retried?;
            }
        }
        match &skipped {
            None => write_hsp(writer, hit, hsp_index)?,
            Some(id) => write_skipped_alignment(writer, id)?,
        }
        if let Some(probe) = probe.as_mut() {
            writer.flush()?;
            probe.end(hsp_index);
        }
    }
    Ok(())
}

/// The alignments of one database subject (an `-db` search) in outfmt 0: each HSP by
/// `write_hsp`, the first one after the subject's heading. The heading of a sequence that
/// has a BLAST definition line with the id `gnl|BL_ORD_ID|n` is `>` and the stored title
/// (no id and no space), or the title that `GenerateDefline` makes when the stored one is
/// empty, wrapped as every heading, then `Length=`. `probe` marks the heading and each
/// HSP without changing the written bytes.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2249-2254
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2383-2386
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2421-2422
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2461-2479
/// Nothing was written after `>` for a `BL_ORD_ID` id, so `tellp()` is 1 and no space
/// follows it. A stored title is a VisibleString (bytes from ' ' to '~'); the title that
/// `GenerateDefline` makes of an empty one from a name of the taxonomy database that NCBI
/// cannot decode throws, and NCBI's line in place of each alignment is written instead
/// (`DatabaseHeading::skipped_id`). `protein` is the database's molecule type (BLASTP and
/// BLASTN share this heading).
///
/// The alignments of the HSPs `part` of the subject's (the heading with the first), as
/// `write_subject_alignments_part`. The part with the heading makes its title
/// once (`write_database_subject_heading`), which also tells whether the alignments are skipped;
/// a later part (EXPERIMENT LOSAT_X_PAIRWISEPAR's pieces) asks only that
/// (`DatabaseHeading::skipped_id`).
#[allow(clippy::too_many_arguments)]
fn write_database_subject_alignments_part<W: Write>(
    writer: &mut W,
    mut probe: Option<&mut FormatProbe<'_>>,
    title: DatabaseHeading<'_>,
    protein: bool,
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    oids: Option<&[u32]>,
    s_idx: u32,
    hits: &[&PairwiseHit],
    hsp_indices: &[HspIndex],
    part: std::ops::Range<usize>,
    html: Option<HtmlHeading<'_>>,
    place: super::query_warnings::ReportHook<'_, W>,
    write_hsp: impl FnMut(&mut W, &PairwiseHit, HspIndex) -> io::Result<()>,
) -> io::Result<()> {
    if hits.iter().any(|hit| hit.skips_pairwise_hsp()) {
        return write_subject_with_skips(
            writer,
            probe,
            SubjectHeading::Database(title),
            protein,
            parsed,
            oids,
            s_idx,
            hits,
            hsp_indices,
            part,
            html,
            place,
            write_hsp,
        );
    }
    let first_index = hsp_indices[0];
    let skipped = if part.start == 0 {
        if let Some(probe) = probe.as_mut() {
            writer.flush()?;
            probe.subject_begin(first_index);
        }
        let skipped = write_database_subject_heading(
            writer,
            title,
            protein,
            hits[0].subject_length,
            html,
            &mut *place,
        )?;
        if let Some(probe) = probe.as_mut() {
            writer.flush()?;
            probe.subject_end(first_index);
        }
        skipped
    } else {
        let skipped = title.skipped_id(protein, html);
        if let Some(html) = html {
            later_part_posts(html, &skipped, writer, &mut *place)?;
        }
        skipped?
    };
    let retry = html.map(|html| (move || title.skipped_id(protein, Some(html)), html));
    write_hsps_or_skips(
        writer,
        probe,
        hits,
        hsp_indices,
        part,
        skipped,
        retry,
        place,
        write_hsp,
    )
}

/// Per-HSP skip handling for a subject (`Database` or `Record` heading) that contains at least
/// one HSP skipped inside NCBI's per-HSP `try` block before `x_PrintDefLine`
/// (`PairwiseHit::skips_pairwise_hsp`: BLASTN class N `QueryLeftOfRecord` or TBLASTN
/// `NullbSkipTranslate`). Each skipped HSP writes NCBI's `...alignment skipped` line without
/// starting a heading or updating `previousId`; the first non-skipped HSP writes the heading
/// (and sets `previousId`), and a later skipped HSP leaves `previousId` unchanged so subsequent
/// non-skipped HSPs do not repeat the heading. Handles split parts (`part.start > first_shown`)
/// and `-html` heading posts (`later_part_posts` and retries) uniformly for both heading kinds.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1941-1984
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3956-3981
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn write_subject_with_skips<W: Write>(
    writer: &mut W,
    mut probe: Option<&mut FormatProbe<'_>>,
    heading: SubjectHeading<'_>,
    protein: bool,
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    oids: Option<&[u32]>,
    s_idx: u32,
    hits: &[&PairwiseHit],
    hsp_indices: &[HspIndex],
    part: std::ops::Range<usize>,
    html: Option<HtmlHeading<'_>>,
    place: super::query_warnings::ReportHook<'_, W>,
    mut write_hsp: impl FnMut(&mut W, &PairwiseHit, HspIndex) -> io::Result<()>,
) -> io::Result<()> {
    let skip_id = heading.skip_id(parsed, oids, s_idx);
    let Some(first_shown) = hits.iter().position(|hit| !hit.skips_pairwise_hsp()) else {
        for &hsp_index in &hsp_indices[part] {
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.begin(hsp_index);
            }
            write_skipped_alignment(writer, &skip_id)?;
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.end(hsp_index);
            }
        }
        return Ok(());
    };
    let mut heading_done = part.start > first_shown;
    let mut heading_skipped: Option<Vec<u8>> = None;
    let mut seen_shown_in_part = false;
    for idx in part {
        let (hit, hsp_index) = (hits[idx], hsp_indices[idx]);
        if hit.skips_pairwise_hsp() {
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.begin(hsp_index);
            }
            write_skipped_alignment(writer, &skip_id)?;
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.end(hsp_index);
            }
            continue;
        }
        if !heading_done {
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.subject_begin(hsp_index);
            }
            heading_skipped = heading.write_heading(
                writer,
                protein,
                hit.subject_length,
                html,
                parsed,
                s_idx,
                &mut *place,
            )?;
            if let Some(probe) = probe.as_mut() {
                writer.flush()?;
                probe.subject_end(hsp_index);
            }
            heading_done = true;
            seen_shown_in_part = true;
        } else if !seen_shown_in_part {
            let skipped = heading.skipped_id(protein, html, parsed, s_idx);
            if let Some(html) = html {
                later_part_posts(html, &skipped, writer, &mut *place)?;
            }
            heading_skipped = skipped?;
            seen_shown_in_part = true;
        } else if heading_skipped.is_some() {
            if let Some(html) = html {
                let retried = heading.skipped_id(protein, Some(html), parsed, s_idx);
                place_heading_posts(html, true, writer, &mut *place)?;
                heading_skipped = retried?;
            }
        }
        if let Some(probe) = probe.as_mut() {
            writer.flush()?;
            probe.begin(hsp_index);
        }
        match &heading_skipped {
            None => write_hsp(writer, hit, hsp_index)?,
            Some(id) => write_skipped_alignment(writer, id)?,
        }
        if let Some(probe) = probe.as_mut() {
            writer.flush()?;
            probe.end(hsp_index);
        }
    }
    Ok(())
}

/// NCBI's line in place of an alignment whose heading it cannot make (c++/src/objtools/align_format/showalign.cpp:1976-1981,
/// cited at `write_subject_alignments_part`): of a database subject (`DatabaseHeading::skipped_id`)
/// and of a `-subject` record (`Subject_N`), `id` the Seq-id text of the Seq-align's subject.
fn write_skipped_alignment<W: Write>(writer: &mut W, id: &[u8]) -> io::Result<()> {
    writer.write_all(b"Sequence with id ")?;
    writer.write_all(id)?;
    writer.write_all(b" no longer exists in database...alignment skipped\n")
}

/// The definition lines of a database sequence's alignment heading: the stored title of a
/// sequence without parsed Seq-ids, with the source from which NCBI makes an empty one
/// (`report::defline::SourceTaxname`), or the definition lines of one with them
/// (`ParsedSubject::heading`).
#[derive(Clone, Copy)]
pub enum DatabaseHeading<'a> {
    Title(&'a [u8], &'a super::defline::SourceTaxname),
    /// With whether every definition line is shown (`ParsedSubject::heading`).
    Parsed(&'a crate::report::db_subject::ParsedSubject, bool),
    /// No heading and one line in place of each alignment, whose Seq-id text this is: where the
    /// query cannot be read (`UnreadableQuery`), or the subject's load throws where the alignment
    /// loads it (its batch's prefetch is off, `report::db_subject::HitLoad::Throws`;
    /// `ParsedSubject::skipped_id`).
    Skipped(&'a [u8]),
}

impl<'a> DatabaseHeading<'a> {
    /// The Seq-id text that outfmt 0 writes in place of the subject's alignments (`Sequence with
    /// id <id> no longer exists in database...alignment skipped`, showalign.cpp:1976-1981, quoted
    /// at `ParsedSubject::heading_skipped_id`) where the heading's generated title throws: an empty
    /// stored title made from a name that NCBI cannot decode (`SourceTaxname::undecodable_id`), or
    /// the generated title of a sequence with parsed Seq-ids (`heading_skipped_id`). With `-html`
    /// (`html`) the heading is made (`html_lines`: its id links come first and the registry can
    /// throw), as NCBI makes it for each alignment until it succeeds.
    fn skipped_id(
        self,
        protein: bool,
        html: Option<HtmlHeading<'_>>,
    ) -> io::Result<Option<Vec<u8>>> {
        if let Some(html) = html {
            return Ok(self.html_lines(protein, html)?.err());
        }
        Ok(match self {
            DatabaseHeading::Title([], source) => source.undecodable_id().map(<[u8]>::to_vec),
            DatabaseHeading::Title(..) => None,
            DatabaseHeading::Parsed(subject, every_line) => {
                subject.heading_skipped_id(protein, every_line)
            }
            DatabaseHeading::Skipped(id) => Some(id.to_vec()),
        })
    }

    /// The `-html` heading's lines (`write_lines`' with links, anchors and HTML-encoded titles),
    /// or the Seq-id text of NCBI's line in place of each alignment: the definition lines of a
    /// sequence with parsed Seq-ids (`ParsedSubject::heading`), else `>`, `<a name=BL_ORD_ID:N></a>`,
    /// a space and the stored title (`HtmlHeading::unparsed_start`; `bp.dp.f0.h`).
    fn html_lines(
        self,
        protein: bool,
        html: HtmlHeading<'_>,
    ) -> io::Result<Result<Vec<u8>, Vec<u8>>> {
        let (title, source) = match self {
            DatabaseHeading::Parsed(subject, every_line) => {
                return subject.heading(protein, every_line, Some(html));
            }
            // The heading is not made (`UnreadableQuery`), so neither are its links.
            DatabaseHeading::Skipped(id) => return Ok(Err(id.to_vec())),
            DatabaseHeading::Title(title, source) => (title, source),
        };
        let mut lines = Vec::new();
        if let Err(error) = html.unparsed_start(&mut lines) {
            return HtmlHeading::skips(error, || html.unparsed_skipped_id()).map(Err);
        }
        Ok(match stored_heading_title(title, source, protein) {
            Ok(title) => {
                super::db_subject::wrap_title(&super::html::html_encode(&title), &mut lines);
                lines.push(b'\n');
                Ok(lines)
            }
            Err(id) => Err(id.to_vec()),
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:340-359
    /// Writes the definition-line lines of the heading (each ending with a newline; `Length=`
    /// follows): `>` and the stored title wrapped as `s_WrapOutputLine` wraps it, an empty one
    /// made by `GenerateDefline` from the source (`showalign.cpp:2250-2255`), or of a sequence with
    /// parsed Seq-ids the line of each definition line (`ParsedSubject::heading`). Where the title
    /// made for the heading throws it writes nothing and gives the Seq-id text of NCBI's line in
    /// place of each alignment (`skipped_id`). The title is made once.
    /// With `-html` (`html`) the lines are `html_lines`.
    fn write_lines<W: Write>(
        self,
        writer: &mut W,
        protein: bool,
        html: Option<HtmlHeading<'_>>,
        place: super::query_warnings::ReportHook<'_, W>,
    ) -> io::Result<Option<Vec<u8>>> {
        let lines = match (self, html) {
            (_, Some(html)) => {
                let lines = self.html_lines(protein, html);
                place_heading_posts(html, true, writer, place)?;
                lines?
            }
            (DatabaseHeading::Parsed(subject, every_line), None) => {
                subject.heading(protein, every_line, None)?
            }
            (DatabaseHeading::Title(title, source), None) => {
                return match stored_heading_title(title, source, protein) {
                    Ok(heading) => {
                        writer.write_all(b">")?;
                        write_wrapped(writer, &heading)?;
                        writer.write_all(b"\n").map(|()| None)
                    }
                    Err(id) => Ok(Some(id.to_vec())),
                };
            }
            (DatabaseHeading::Skipped(id), None) => return Ok(Some(id.to_vec())),
        };
        match lines {
            Ok(lines) => writer.write_all(&lines).map(|()| None),
            Err(id) => Ok(Some(id)),
        }
    }

    /// The heading of a subject with parsed Seq-ids and `hsps` HSPs in a search that is
    /// `gapped` or not (`ParsedSubject::heading`): the subject's own HSPs, not those of the
    /// copies shown in its hit (`QueryGroup::own_hsps`), as the Seq-align of each subject's
    /// HSPs is split.
    fn parsed(
        subject: &'a crate::report::db_subject::ParsedSubject,
        gapped: bool,
        hsps: usize,
    ) -> Self {
        DatabaseHeading::Parsed(subject, !gapped && hsps > 1)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1952-1981
    /// The heading of a subject with parsed Seq-ids as `parsed` makes it, or `Skipped` with its
    /// Seq-align id where its load throws (`thrown`, `x_pairwise_par::QuerySections::thrown`) or
    /// its scope gets no Bioseq (`ParsedSubject::scope_finds_none`: the alignment's handle throws
    /// before any heading starts, so that none counts in the ranks of the headings after it).
    fn loaded(
        subject: &'a crate::report::db_subject::ParsedSubject,
        gapped: bool,
        hsps: usize,
        thrown: bool,
    ) -> Self {
        if thrown || subject.scope_finds_none() {
            return DatabaseHeading::Skipped(subject.skipped_id());
        }
        Self::parsed(subject, gapped, hsps)
    }
}

/// `x_pairwise_par::QuerySections::thrown` of a report whose subjects' loads are `loads` and whose
/// subjects with parsed Seq-ids are `parsed`, by subject index: the exception of subject `s_idx`'s
/// load where it throws (`HitLoads::thrown`), and whether its Seq-align id is local
/// (`ParsedSubject::alignment_id_is_local`).
fn thrown_load<'l>(
    loads: &'l crate::report::db_subject::HitLoads,
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    s_idx: u32,
) -> Option<(&'l crate::ncbi_exception::NcbiException, bool)> {
    let exception = loads.thrown(s_idx as usize)?;
    Some((exception, alignment_id_local(parsed, s_idx)))
}

/// Whether subject `s_idx` of a report whose subjects with parsed Seq-ids are `parsed` (by subject
/// index) has a local Seq-align id (`ParsedSubject::alignment_id_is_local`; one without parsed
/// Seq-ids has `gnl|BL_ORD_ID|n` or, a `-subject` record, no load to order).
fn alignment_id_local(
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    s_idx: u32,
) -> bool {
    parsed
        .get(s_idx as usize)
        .and_then(Option::as_deref)
        .is_some_and(crate::report::db_subject::ParsedSubject::alignment_id_is_local)
}

/// The title of the heading of a database sequence without parsed Seq-ids: its stored `title`, or
/// of an empty one the title that NCBI makes from the Bioseq's `source`
/// (`c++/src/objtools/align_format/showalign.cpp:2250-2255`, cited at
/// `ParsedSubject::generated_title_fetched`; LOSAT refuses a source that it refuses
/// before the report, `StoredTitle::check_shown_source`); `Err` with the Seq-id text of NCBI's
/// line in place of each alignment where that throws (`SourceTaxname::undecodable_id`).
fn stored_heading_title<'t>(
    title: &'t [u8],
    source: &'t super::defline::SourceTaxname,
    protein: bool,
) -> Result<std::borrow::Cow<'t, [u8]>, &'t [u8]> {
    if !title.is_empty() {
        return Ok(title.into());
    }
    if let Some(id) = source.undecodable_id() {
        return Err(id);
    }
    let facts = super::defline::BioseqFacts::of_source(source.name().ok().flatten());
    Ok(
        super::defline::generate_bioseq_defline(title, protein, false, &facts)
            .map(|defline| defline.text)
            .unwrap_or_default()
            .into(),
    )
}

/// The source of the subject `s_idx` without parsed Seq-ids among `sources`, the sources of a
/// report's database sequences by subject index (`algorithm::common::database::keep_source`;
/// empty where none has one).
pub fn stored_source(
    sources: &[super::defline::SourceTaxname],
    s_idx: u32,
) -> &super::defline::SourceTaxname {
    sources
        .get(s_idx as usize)
        .unwrap_or(&super::defline::UNNAMED_SOURCE)
}

/// The heading `heading` of a database sequence's alignments
/// (`write_database_subject_alignments_part`): its definition-line lines (`DatabaseHeading::write_lines`), then `Length=` when the length is
/// known, and one blank line; or nothing, and the Seq-id text of NCBI's line in place of each
/// alignment, where NCBI cannot make the heading's title.
fn write_database_subject_heading<W: Write>(
    writer: &mut W,
    heading: DatabaseHeading<'_>,
    protein: bool,
    length: Option<usize>,
    html: Option<HtmlHeading<'_>>,
    place: super::query_warnings::ReportHook<'_, W>,
) -> io::Result<Option<Vec<u8>>> {
    if let Some(id) = heading.write_lines(writer, protein, html, place)? {
        return Ok(Some(id));
    }
    if let Some(length) = length {
        writeln!(writer, "Length={length}")?;
    }
    writeln!(writer).map(|()| None)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
fn ncbi_percent_match(numerator: usize, denominator: usize) -> usize {
    if numerator == denominator {
        100
    } else {
        ((0.5 + 100.0 * numerator as f64 / denominator as f64) as usize).min(99)
    }
}

/// The `Expect` label of an HSP's header: `Expect` followed by n in parentheses for an HSP of a
/// linked set of n > 1 HSPs
/// (`sum_n`, `BlastHSP::num`), else `Expect` (showalign.cpp:3593-3598, cited at
/// `write_hsp_info`).
#[inline]
fn expect_label(sum_n: Option<i32>) -> String {
    sum_n
        .filter(|&count| count > 1)
        .map_or_else(|| "Expect".to_string(), |count| format!("Expect({count})"))
}

/// Write HSP score information
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp x_DisplayAlignInfo()
///
/// Format:
/// ```text
///  Score = XXX bits (YYY),  Expect = Z.Ze-NN
///  Identities = AA/BB (CC%), Positives = DD/BB (EE%), Gaps = FF/BB (GG%)
///  Frame = +X/+Y
/// ```
pub fn write_hsp_info<W: Write>(
    writer: &mut W,
    hit: &PairwiseHit,
    config: &PairwiseConfig,
) -> io::Result<()> {
    let h = &hit.hit;

    // Score line
    // NCBI: " Score = XXX bits (YYY),  Expect = Z.Ze-NN"
    let bit_score_str = format_bitscore_ncbi(h.bit_score);
    let evalue_str = format_evalue_ncbi(h.e_value);
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3593-3598
    // An HSP of a linked set (psiblast's `-sum_stats true`) has its "sum_n" score, the set's size.
    let expect = expect_label(hit.sum_n);
    // psiblast's first round has blastp's methods; rpsblast's redo marks its HSPs likewise.
    if matches!(config.program.as_str(), "blastp" | "psiblast" | "rpsblast") {
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3600-3603
        let method = match hit.comp_adjust_method {
            Some(1) => ", Method: Composition-based stats.",
            Some(2) => ", Method: Compositional matrix adjust.",
            _ => "",
        };
        writeln!(
            writer,
            " Score = {} bits ({}),  {} = {}{}",
            bit_score_str, h.raw_score, expect, evalue_str, method
        )?;
    } else {
        writeln!(
            writer,
            " Score = {} bits ({}),  {} = {}",
            bit_score_str, h.raw_score, expect, evalue_str
        )?;
    }
    write_hsp_identities(writer, hit, config)
}

/// The header of a BLASTN HSP under matrix-only scoring (reward 0, penalty 0: NCBI's
/// `m_DisableKAStats`, which sets `eShowRawScoreOnly`): the raw score alone and an empty line,
/// then the lines of `write_hsp_identities`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3586-3591,3605
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:654-655
pub fn write_hsp_info_raw_score<W: Write>(
    writer: &mut W,
    hit: &PairwiseHit,
    config: &PairwiseConfig,
) -> io::Result<()> {
    writeln!(writer, " Score = {}", hit.hit.raw_score)?;
    writeln!(writer)?;
    write_hsp_identities(writer, hit, config)
}

/// The lines of an HSP's header after its score line: identities (and positives), gaps, the
/// strand or the frames, and the empty line before the alignment (`write_hsp_info`).
fn write_hsp_identities<W: Write>(
    writer: &mut W,
    hit: &PairwiseHit,
    config: &PairwiseConfig,
) -> io::Result<()> {
    let h = &hit.hit;
    // Identity/Positives/Gaps line
    let align_len = h.length;
    let num_ident = h.num_ident;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
    let ident_pct = ncbi_percent_match(num_ident, align_len);

    let positives = hit.positives.unwrap_or(h.num_positives.max(num_ident));
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
    let pos_pct = ncbi_percent_match(positives, align_len);

    let gaps = hit.query_range_lengths().gaps;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
    let gap_pct = ncbi_percent_match(gaps, align_len);

    // For protein/translated: show Identities, Positives, Gaps
    // For nucleotide: show Identities, Gaps (no Positives)
    if config.program.contains("blast") && !config.program.contains("blastn") {
        writeln!(
            writer,
            " Identities = {}/{} ({}%), Positives = {}/{} ({}%), Gaps = {}/{} ({}%)",
            num_ident,
            align_len,
            ident_pct,
            positives,
            align_len,
            pos_pct,
            gaps,
            align_len,
            gap_pct
        )?;
    } else if config.program == "blastn" {
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1794-1805,2131-2141
        // The displayed counts come from the two rows (before the lowercase masking of
        // the display), and the gaps are the gap columns of the alignment (`x_GetNumGaps`'s
        // chunks): `gaps` above, the hit's count of its edit script's gaps, which BLASTN's rows
        // always give (`algorithm::blastn::pairwise::Rows`; a query letter past the query's
        // record shows the gap character too, remaining-gaps X1b-p).
        let (query_row, subject_row) = (
            hit.query_seq.as_deref().unwrap_or_default().as_bytes(),
            hit.subject_seq.as_deref().unwrap_or_default().as_bytes(),
        );
        let columns = query_row.len();
        let matches = query_row
            .iter()
            .zip(subject_row)
            .filter(|(q, s)| q.eq_ignore_ascii_case(s))
            .count();
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:310-320
        // The query is always shown on its plus strand (showalign.cpp:1852-1858).
        writeln!(
            writer,
            " Identities = {}/{} ({}%), Gaps = {}/{} ({}%)",
            matches,
            columns,
            ncbi_percent_match(matches, columns),
            gaps,
            columns,
            ncbi_percent_match(gaps, columns)
        )?;
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4011-4018
        // The strand is the HSP's (its query frame; the shown query is on its plus
        // strand), not the order of the coordinates, which is the same for an HSP of one
        // letter.
        let subject_strand = if h.query_frame < 0 { "Minus" } else { "Plus" };
        writeln!(writer, " Strand=Plus/{subject_strand}")?;
    } else {
        writeln!(
            writer,
            " Identities = {}/{} ({:.0}%), Gaps = {}/{} ({:.0}%)",
            num_ident, align_len, ident_pct, gaps, align_len, gap_pct
        )?;
    }

    // Frame line (for translated searches)
    if config.show_frame {
        if let (Some(qf), Some(sf)) = (hit.query_frame, hit.subject_frame) {
            let qf_str = if qf > 0 {
                format!("+{}", qf)
            } else {
                format!("{}", qf)
            };
            let sf_str = if sf > 0 {
                format!("+{}", sf)
            } else {
                format!("{}", sf)
            };
            writeln!(writer, " Frame = {}/{}", qf_str, sf_str)?;
        }
    }

    writeln!(writer)?;
    Ok(())
}

/// Write alignment rows
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp x_DisplayRowData()
///
/// Format:
/// ```text
/// Query  1    MVLSPADKTN VKAAWGKVGA HAGEYGAEAL ERMFLSFPTT KTYFPHFDLS H  60
///             MVLSPADKTN VKAAWGKVGA HAGEYGAEAL ERMFLSFPTT KTYFPHFDLS H
/// Sbjct  1    MVLSPADKTN VKAAWGKVGA HAGEYGAEAL ERMFLSFPTT KTYFPHFDLS H  60
/// ```
///
/// `masks`: the rows' masks under `-html` (`HtmlRowMasks`), `None` without.
pub fn write_alignment<W: Write>(
    writer: &mut W,
    hit: &PairwiseHit,
    masks: Option<&super::query_anchored::HspMasks>,
    config: &PairwiseConfig,
) -> io::Result<()> {
    let h = &hit.hit;

    // If we have actual sequences, display them
    if let (Some(ref qseq), Some(ref sseq)) = (&hit.query_seq, &hit.subject_seq) {
        write_alignment_with_sequences(writer, h, qseq, sseq, hit.query_case, masks, config)?;
    } else {
        // No sequences available - show placeholder
        writeln!(
            writer,
            "Query  {}  [... {} aa ...]  {}",
            h.q_start, h.length, h.q_end
        )?;
        writeln!(writer)?;
        writeln!(
            writer,
            "Sbjct  {}  [... {} aa ...]  {}",
            h.s_start, h.length, h.s_end
        )?;
    }

    writeln!(writer)?;
    Ok(())
}

/// Write alignment with actual sequences: the rows of BLASTP, RPSBLAST, PSIBLAST and DELTABLAST
/// (`seqalign::protein_row_segments`) and of BLASTN (`seqalign::nucleotide_row_segments`; a
/// minus-strand HSP's Dense-seg is reversed, `query_anchored::HspAlnVec`, so its subject row runs
/// from high to low), rows of letters (`write_hsp_rows`), the columns of `hit`'s edit script.
/// `query_case` is how the query row's letters compare (`PairwiseHit::query_case`), `masks` the
/// rows' masks under `-html` (`write_hsp_rows`).
fn write_alignment_with_sequences<W: Write>(
    writer: &mut W,
    hit: &Hit,
    query_seq: &str,
    subject_seq: &str,
    query_case: super::row_letters::QueryCase,
    masks: Option<&super::query_anchored::HspMasks>,
    config: &PairwiseConfig,
) -> io::Result<()> {
    let line_len = config.line_length;
    let q_chars: Vec<char> = query_seq.chars().collect();
    let s_chars: Vec<char> = subject_seq.chars().collect();

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:96-98
    //
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1600-1625
    let max_start_len = coordinate_width(hit);
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1570-1573
    let nucleotide = config.program == "blastn";
    // The rows' end gaps (`row_spans`). A BLASTN alignment (X-drop or greedy) starts and ends
    // with a substitution and has none: a `-` at the end of its query row fills a segment past the
    // query's end, the gap character (alnvec.cpp:288-297, cited at `row_spans`;
    // `algorithm::blastn::pairwise::context_letters`).
    let ops = match nucleotide {
        true => &[][..],
        false => hit.gap_info.as_deref().unwrap_or_default(),
    };
    let [query_span, subject_span] = row_spans(ops, [q_chars.len(), s_chars.len()]);

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2131-2155
    // The rows may carry lowercase masking, which NCBI applies only when it prints a
    // row (showalign.cpp:2495-2521), so they are compared without case, and the protein
    // middle line shows the uppercase residue; a PSSM query's Bioseq letters compare as they are
    // (`row_letters::QueryCase`).
    let compared = query_case.compared_letters();
    let write_middle = |writer: &mut W, columns: std::ops::Range<usize>| {
        let middle = q_chars[columns.clone()]
            .iter()
            .zip(s_chars[columns].iter())
            .map(|(&q, s)| {
                let q = char::from(compared[usize::from(q as u8)]);
                let s = s.to_ascii_uppercase();
                if nucleotide {
                    if q == s {
                        '|'
                    } else {
                        ' '
                    }
                } else if q == s {
                    q // Identity: show the character
                } else if is_positive_match(q, s, config.protein_matrix) {
                    '+' // Positive: show +
                } else {
                    ' ' // Mismatch: show space
                }
            });
        write_buffered(writer, middle)
    };

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1497-1579
    // A row is written on a line that meets its columns between its end gaps
    // (`rowRng` is set from `m_AV->GetSeqAlnRange` of the row, showalign.cpp:1272; `write_hsp_rows`), the middle
    // line with the query's; its end gaps are NCBI's end character (`row_string`).
    write_hsp_rows(
        writer,
        |segments| {
            if nucleotide {
                super::serial::seqalign::nucleotide_row_segments(hit, segments)
            } else {
                blastp_row_segments(hit, segments)
            }
        },
        [1, 1],
        line_len,
        [
            &row_string(query_seq, &q_chars, query_span),
            &row_string(subject_seq, &s_chars, subject_span),
        ],
        masks,
        max_start_len,
        write_middle,
    )
}

/// The rows' masks of a report's HSPs under `-html` (`HspMasks` by HSP, the masked regions of
/// NCBI's row information of each HSP): `query_anchored::output_seq` puts each mask
/// of a row's line in NCBI's span. The rows' builders keep them only then (the rows hold the
/// masks' letters in lower case in every report, and `PairwiseHit` holds no list: the plain
/// report reads none). `masks[i]` is that of HSP `first + i` of the run's hits (`HspIndex`);
/// without `-html` the table is empty.
#[derive(Clone, Copy, Default)]
pub(crate) struct HtmlRowMasks<'a> {
    masks: &'a [super::query_anchored::HspMasks],
    first: HspIndex,
}

impl<'a> HtmlRowMasks<'a> {
    /// The table of the masks `masks` (`None` without `-html`) of the hits from HSP `first` on.
    pub(crate) fn new(
        masks: Option<&'a [super::query_anchored::HspMasks]>,
        first: HspIndex,
    ) -> Self {
        Self {
            masks: masks.unwrap_or_default(),
            first,
        }
    }

    /// The masks of the rows of HSP `hsp_index`; `None` without `-html`.
    #[inline]
    pub(crate) fn of(self, hsp_index: HspIndex) -> Option<&'a super::query_anchored::HspMasks> {
        self.masks.get(hsp_index.checked_sub(self.first)?)
    }
}

thread_local! {
    /// The vectors of outfmt 0's rows (`write_hsp_rows`): an HSP's segments and its `CAlnVec`,
    /// kept from HSP to HSP of the thread's reports. No NCBI counterpart: NCBI makes a new
    /// Seq-align and `CAlnVec` per HSP (showalign.cpp:1823-1861).
    static HSP_ROWS: std::cell::RefCell<(
        super::serial::seqalign::Segments,
        super::query_anchored::HspAlnVec,
        Vec<u8>,
    )> = std::cell::RefCell::default();
}

/// The labels of an HSP's rows in outfmt 0 (`eShowBlastStyleId`, `x_FillSeqid`).
const ROW_LABELS: [&str; 2] = ["Query", "Sbjct"];

/// Writes an HSP's rows in outfmt 0, `line_length` columns a line: each line's query row, the
/// middle line, the subject row and a blank line. The rows' coordinates come from the HSP's
/// `CAlnVec` (`query_anchored::HspAlnVec`) of its Seq-align's segments, which `segments_of` makes
/// (the producer of the program, `report::serial::seqalign`, as the query-anchored reports'
/// `section_hsps`), of row widths `widths` (`HspStrands::no_widths`: none); `letters` are the rows'
/// letters (NCBI's row strings: an end gap is the end character, a space; a mask's letters in
/// lower case), a line of a row written through `x_OutputSeq` (`query_anchored::write_output_seq`):
/// under `-html` with the row's masks of `masks` (`HtmlRowMasks`), each in NCBI's span, else as
/// the letters are (`masks` is `None`), and `write_middle(writer, columns)` writes the middle
/// line's. `width` is the width
/// of the coordinates (`maxStartLen`, `coordinate_width`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1435-1438,1452-1453
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1486-1497
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1538-1543,1550-1551
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1781-1783
/// The lines are those of the HSP's Dense-seg (`GetAlnStop`), whose columns the rows' letters
/// are: a difference is an internal error. An HSP of length 0 (BLASTX's empty rows) has none
/// (`aln_stop` cast to `int` is -1). `prev_stop` starts at 0 (`vector<int>`).
#[allow(clippy::too_many_arguments)]
fn write_hsp_rows<W: Write>(
    writer: &mut W,
    segments_of: impl FnOnce(
        &mut super::serial::seqalign::Segments,
    ) -> super::query_anchored::HspStrands,
    widths: [u32; 2],
    line_length: usize,
    letters: [&[u8]; 2],
    masks: Option<&super::query_anchored::HspMasks>,
    width: usize,
    mut write_middle: impl FnMut(&mut W, std::ops::Range<usize>) -> io::Result<()>,
) -> io::Result<()> {
    let columns = letters[0].len();
    HSP_ROWS.with_borrow_mut(|(segments, rows, buffer)| {
        let strands = segments_of(segments);
        let widths = if strands.no_widths { [1, 1] } else { widths };
        rows.set(
            segments.starts(),
            segments.lens(),
            segments.strands(),
            widths,
            line_length,
        )?;
        let aln_len = usize::try_from(rows.aln_stop() as i32).map_or(0, |stop| stop + 1);
        if aln_len != columns {
            return Err(io::Error::other(
                "LOSAT internal error in the pairwise report: an HSP's rows and its Dense-seg \
                 have different lengths",
            ));
        }
        let mut prev_stop = [0i32; 2];
        for (index, aln_start) in (0..aln_len).step_by(line_length).enumerate() {
            let line = aln_start..(aln_start + line_length).min(aln_len);
            let cur_range = (aln_start as i64, line.end as i64 - 1);
            for (row, label) in ROW_LABELS.into_iter().enumerate() {
                let Some((left, right)) = rows.line(row, index, cur_range.0, cur_range.1)? else {
                    continue;
                };
                let (start, end) = (left.wrapping_add(1) as i32, right.wrapping_add(1) as i32);
                write_sequence_row(
                    writer,
                    label,
                    start,
                    end,
                    aln_start,
                    prev_stop[row],
                    width,
                    |writer| {
                        super::query_anchored::write_output_seq(
                            writer,
                            buffer,
                            letters[row],
                            aln_start,
                            line.len(),
                            masks.and_then(|masks| masks.row(row)),
                            masks.is_some(),
                        )
                    },
                )?;
                if row == 0 {
                    write_spaces(writer, label.len() + 2 + width + 2)?;
                    write_middle(writer, line.clone())?;
                    writer.write_all(b"\n")?;
                }
                prev_stop[row] = end;
            }
            writer.write_all(b"\n")?;
        }
        Ok(())
    })
}

/// Writes one row of one line of an HSP (`write_hsp_rows`): its label, its first and last
/// position on the line (`start` and `end`, 1-based) and its letters (`write_letters`), on the
/// line of first column `j`; `prev_stop` is the row's `end` on its last line shown.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1597-1626
/// A later line without residues ends where the row's previous line ended (`end == prev_stop`).
/// `AddSpace` takes a `size_t`, so an omitted start gives two more than `maxStartLen` spaces. The
/// labels of the pairwise rows are equally long (`maxIdLen`): two spaces follow. A first line of
/// one residue at 1 of a row of letters shows no coordinates (the `-line_length 1` rows of an HSP
/// from a sequence's first letter); a translated row's first residue spans a codon, so its start
/// and end differ.
#[allow(clippy::too_many_arguments)]
fn write_sequence_row<W: Write>(
    writer: &mut W,
    label: &str,
    start: i32,
    end: i32,
    j: usize,
    prev_stop: i32,
    width: usize,
    write_letters: impl FnOnce(&mut W) -> io::Result<()>,
) -> io::Result<()> {
    let empty = (j > 0 && end == prev_stop) || (j == 0 && start == 1 && end == 1);
    let (mut start_text, mut end_text) = (itoa_i32(start), itoa_i32(end));
    let (start_text, end_text): (&[u8], &[u8]) = if empty {
        (b"", b"")
    } else {
        (start_text.text(), end_text.text())
    };
    // `startLen`: the length of the start's text (the size of `NStr::IntToString` of `start`), 0 without.
    write_joined(
        writer,
        &[label.as_bytes(), b"  ", start_text],
        width + 2 - start_text.len(),
    )?;
    write_letters(writer)?;
    write_joined(writer, &[b"  ", end_text, b"\n"], 0)
}

/// Writes `parts` and then `spaces` spaces (a row's label, start and the spaces before its
/// letters, or the spaces, end and line end after them, `write_sequence_row`), in one write
/// where they fit in a buffer of 128 bytes.
fn write_joined<W: Write>(writer: &mut W, parts: &[&[u8]], spaces: usize) -> io::Result<()> {
    let mut buffer = [b' '; 128];
    let total = parts.iter().map(|part| part.len()).sum::<usize>() + spaces;
    if total > buffer.len() {
        for part in parts {
            writer.write_all(part)?;
        }
        return write_spaces(writer, spaces);
    }
    let mut used = 0;
    for part in parts {
        buffer[used..used + part.len()].copy_from_slice(part);
        used += part.len();
    }
    writer.write_all(&buffer[..total])
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2491-2493
/// A letter that `write_buffered` writes: a `char` of a row kept as `char`s (in UTF-8) or a byte
/// of a middle line of rows kept as bytes.
trait RowLetter: Copy {
    /// The most bytes of one letter.
    const MAX_LEN: usize;
    /// Writes its bytes into `buffer` from `at` (`at + MAX_LEN` is at most its length) and gives
    /// their number.
    fn encode_at(self, buffer: &mut [u8], at: usize) -> usize;
}

impl RowLetter for char {
    const MAX_LEN: usize = 4;

    #[inline(always)]
    fn encode_at(self, buffer: &mut [u8], at: usize) -> usize {
        self.encode_utf8(&mut buffer[at..]).len()
    }
}

impl RowLetter for u8 {
    const MAX_LEN: usize = 1;

    #[inline(always)]
    fn encode_at(self, buffer: &mut [u8], at: usize) -> usize {
        buffer[at] = self;
        1
    }
}

/// Writes `letters` (a middle line,
/// `write_alignment_with_sequences` and `write_protein_middle`) a buffer of 256 bytes at a time
/// (written when a letter of the most bytes might not fit), in its caller's loop.
#[inline(always)]
fn write_buffered<W: Write, L: RowLetter>(
    writer: &mut W,
    letters: impl IntoIterator<Item = L>,
) -> io::Result<()> {
    let mut buffer = [0u8; 256];
    let mut used = 0;
    for letter in letters {
        if used > buffer.len() - L::MAX_LEN {
            writer.write_all(&buffer[..used])?;
            used = 0;
        }
        used += letter.encode_at(&mut buffer, used);
    }
    writer.write_all(&buffer[..used])
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2132-2155
/// Writes the protein middle line of the columns `columns` of the rows of bytes `query` and
/// `subject` (BLASTX's, TBLASTN's and TBLASTX's rows, `write_hsp_rows`): a column's residue where
/// the rows' residues are equal, `+` where `matrix` scores them above 0 (`is_positive_match`),
/// else a space. The rows may carry lower case, which NCBI applies only when it prints a row
/// (showalign.cpp:2495-2521), so the residues are compared in upper case: the query row's (its
/// masks), and the subject row's where `upper_subject` (BLASTX's row of letters: a `-subject`'s
/// lower case with `-lcase_masking`, `blastx/x_soft_subject.rs`). TBLASTN's and TBLASTX's
/// subject rows are translations (`tblastx::report::DisplayCode`), in upper case.
#[inline(always)]
fn write_protein_middle<W: Write>(
    writer: &mut W,
    query: &[u8],
    subject: &[u8],
    columns: std::ops::Range<usize>,
    matrix: ScoringMatrix,
    upper_subject: bool,
    query_case: super::row_letters::QueryCase,
) -> io::Result<()> {
    // A record's mask compares as its capital, a PSSM query's Bioseq letter as it is
    // (`row_letters::QueryCase`).
    let compared = query_case.compared_letters();
    let middle = query[columns.clone()]
        .iter()
        .zip(&subject[columns])
        .map(|(&q, &s)| {
            let q = compared[usize::from(q)];
            let s = if upper_subject {
                s.to_ascii_uppercase()
            } else {
                s
            };
            if q == s {
                q
            } else if is_positive_match(q as char, s as char, matrix) {
                b'+'
            } else {
                b' '
            }
        });
    write_buffered(writer, middle)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:340-359
/// Writes a subject's heading `text` as `s_WrapOutputLine` does: after the 60th byte, a line ends
/// at the next space (`isspace`: also a vertical tab, 11), a run of bytes at a time.
fn write_wrapped<W: Write>(writer: &mut W, text: &[u8]) -> io::Result<()> {
    let mut do_wrap = false;
    let mut from = 0;
    for (i, &c) in text.iter().enumerate() {
        if i > 0 && i % 60 == 0 {
            do_wrap = true;
        }
        if do_wrap && (c.is_ascii_whitespace() || c == 11) {
            writer.write_all(&text[from..=i])?;
            writer.write_all(b"\n")?;
            from = i + 1;
            do_wrap = false;
        }
    }
    writer.write_all(&text[from..])
}

/// The row `row` (`text`, its `char`s `chars`) as NCBI's row string (`GetWholeAlnSeqString`):
/// its end gaps (the columns outside `span`, `row_spans`) the end character, a space. A row
/// of ASCII letters without end gaps is its text as it is.
fn row_string<'a>(
    text: &'a str,
    chars: &[char],
    span: (usize, usize),
) -> std::borrow::Cow<'a, [u8]> {
    if text.is_ascii() {
        return ascii_row_string(text.as_bytes(), span);
    }
    let mut row = String::with_capacity(text.len());
    for (column, &letter) in chars.iter().enumerate() {
        row.push(if column < span.0 || column >= span.1 {
            ' '
        } else {
            letter
        });
    }
    std::borrow::Cow::Owned(row.into_bytes())
}

/// The row `row` of ASCII letters as NCBI's row string (`row_string`): its columns outside `span`
/// (`row_spans`) a space.
fn ascii_row_string(row: &[u8], span: (usize, usize)) -> std::borrow::Cow<'_, [u8]> {
    if span == (0, row.len()) {
        return std::borrow::Cow::Borrowed(row);
    }
    let mut owned = row.to_vec();
    owned[..span.0].fill(b' ');
    owned[span.1..].fill(b' ');
    std::borrow::Cow::Owned(owned)
}

/// The rows `[query, subject]` of ASCII letters of a translated program's HSP (BLASTX's,
/// TBLASTN's), the columns of its edit script `ops`, as NCBI's row strings (`ascii_row_string`):
/// a Smith-Waterman alignment that starts with a gap (`-use_sw_tback`, line PS item PS12d) shows
/// it as spaces, as the end gaps of BLASTP's rows (`write_alignment_with_sequences`).
fn ascii_row_strings<'a>(
    rows: [&'a [u8]; 2],
    ops: &[GapEditOp],
) -> [std::borrow::Cow<'a, [u8]>; 2] {
    let spans = row_spans(ops, rows.map(<[u8]>::len));
    [0, 1].map(|row| ascii_row_string(rows[row], spans[row]))
}

// NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnvec.hpp:386-393
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:362-370 (`GetWholeAlnSeqString`)
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:279-297
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1209-1216
// (`eShowEndGaps` is set by the query-anchored views only, blast_format.cpp:1754-1758; the
// tabular writer sets the gap character alone, tabular.cpp:982-986.)
/// The columns `lo..hi` of the query's and the subject's rows (of `lens` columns) of an HSP whose
/// edit script is `ops` between their end gaps (`serial::seqalign::row_end_gaps`): a row's gap
/// operations before its first letter and after its last, which NCBI's alignment strings fill
/// with the end character (a space in the pairwise view and in the tabular `qseq`, `sseq` and
/// `btop`, `report::tabular::push_alignment_row`) where an inner gap is `-`. NCBI decides by where
/// a `-` comes from: a `-` that is a residue (a database's stored `-`, a PSSM query's gap
/// residue) is a letter of the script, as is the gap character that fills a segment past a
/// sequence's end (BLASTN's query row of an HSP that NCBI keeps outside its query context,
/// `algorithm::blastn::pairwise::context_letters`, whose rows are shown without a script,
/// remaining-gaps X1b-p). Only a Smith-Waterman alignment has end gaps (`-use_sw_tback`:
/// `ALIGN_EX` runs from the Smith-Waterman start over a rectangle whose best path may start with
/// a gap); an X-drop alignment starts and ends with a substitution.
fn row_spans(ops: &[GapEditOp], lens: [usize; 2]) -> [(usize, usize); 2] {
    use super::serial::seqalign::row_end_gaps;
    [
        row_end_gaps(ops, AlignRow::Query).span(lens[0]),
        row_end_gaps(ops, AlignRow::Subject).span(lens[1]),
    ]
}

/// Check if two amino acids are a positive match (similar)
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2122-2149
fn is_positive_match(a: char, b: char, matrix: ScoringMatrix) -> bool {
    if a == b {
        return true;
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3563-3581
    protein_display_score(matrix, a as u8, b as u8) > 0
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:75-83
// `kBits` is a static initialized when the program starts, from whether the
// environment has CTOOLKIT_COMPATIBLE (any value, also an empty one).
fn ncbi_k_bits() -> &'static str {
    static K_BITS: std::sync::OnceLock<&'static str> = std::sync::OnceLock::new();
    K_BITS.get_or_init(|| ncbi_k_bits_for(std::env::var_os("CTOOLKIT_COMPATIBLE").is_some()))
}

fn ncbi_k_bits_for(ctoolkit_compatible: bool) -> &'static str {
    if ctoolkit_compatible {
        "(bits)"
    } else {
        "(Bits)"
    }
}

/// BLASTX's description table: `subject_labels` holds each subject's line (`database_description_label`
/// of a database sequence or a `-subject` record).
#[allow(clippy::too_many_arguments)]
fn write_subject_summary_table_with_sum_n<W: Write>(
    writer: &mut W,
    subject_order: &[u32],
    last_row_counted: bool,
    subject_hits: &std::collections::HashMap<u32, Vec<&PairwiseHit>>,
    subject_labels: &[DescriptionLabel],
    show_sum_n: bool,
    columns: Option<SortColumns>,
    place: super::query_warnings::ReportHook<'_, W>,
    // `-html`: the rows' links with what the report does after a throw, the subjects with parsed
    // Seq-ids by subject index and what the report's loads do at each subject (`HtmlRows`).
    html: Option<(
        (
            &super::html_ids::HtmlReport,
            crate::report::db_subject::AfterThrow,
        ),
        &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
        &crate::report::db_subject::HitLoads,
    )>,
    // Where the rows' ranges throw (`write_description_table`).
    rows_throw: Option<EmptyRowsStop>,
) -> io::Result<()> {
    write_description_table(
        writer,
        subject_order,
        last_row_counted,
        subject_hits,
        show_sum_n,
        columns,
        |s_idx| {
            subject_labels.get(s_idx as usize).map_or_else(
                || DescriptionLabel::title_only(b"unknown".to_vec()),
                DescriptionLabel::clone,
            )
        },
        place,
        html.map(|((report, after), parsed, loads)| HtmlRows {
            report,
            after,
            parsed,
            loads,
        }),
        rows_throw,
    )
}

/// The number of decimal digits of `value` (the length of `value.to_string()`).
#[inline]
fn digit_count(value: usize) -> usize {
    value.checked_ilog10().map_or(1, |log| log as usize + 1)
}

/// Width of the start-coordinate column of an alignment block.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1366-1367
/// `GetSeqStart` and `GetSeqStop` are 0-based, so the width is that of the largest
/// 1-based coordinate minus one, and a printed 1-based start can have one digit
/// more. NCBI then pads with the `size_t` value of `maxStartLen` less `startLen`, plus
/// `k_StartSequenceMargin`, which is one space in that case; callers write
/// `width + 2 - digits` (showalign.cpp:1614; align_format_util.cpp:932-936).
fn coordinate_width(hit: &Hit) -> usize {
    digit_count(hit.q_start.max(hit.q_end).max(hit.s_start).max(hit.s_end) - 1)
}

#[inline]
fn write_spaces<W: Write>(writer: &mut W, count: usize) -> io::Result<()> {
    const SPACES: [u8; 64] = [b' '; 64];
    let mut left = count;
    while left > 0 {
        let chunk = left.min(SPACES.len());
        writer.write_all(&SPACES[..chunk])?;
        left -= chunk;
    }
    Ok(())
}

#[inline]
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:797-801
/// The database name of a search of a `-subject` file (`m_IsDbScan`): `User specified sequence
/// set (Input: <name>)` with the bytes of the file's name (`utils::os_bytes::path_bytes`), as
/// NCBI keeps the `std::string` of the command line (a name that is not UTF-8, BLASTX's
/// `-subject`: probe R2 of line LIB-B).
pub(crate) fn subject_set_name(subject: &std::path::Path) -> Vec<u8> {
    [
        b"User specified sequence set (Input: ".as_slice(),
        &crate::utils::os_bytes::path_bytes(subject),
        b")",
    ]
    .concat()
}

/// The definition of a subject set in the pairwise report: [`subject_set_name`], or, for the
/// subjects of an imported strategy (an empty path: they name no input), [`scan_mode_definition`]'s
/// definition without an input.
pub(crate) fn subject_set_definition(subject: &std::path::Path) -> Vec<u8> {
    if subject.as_os_str().is_empty() {
        return scan_mode_definition("").into_bytes();
    }
    subject_set_name(subject)
}

fn ensure_trailing_period(text: &[u8]) -> std::borrow::Cow<'_, [u8]> {
    if text.ends_with(b".") {
        std::borrow::Cow::Borrowed(text)
    } else {
        std::borrow::Cow::Owned([text, b"."].concat())
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:584-603
// C's `%#.3g` (C99 7.21.6.1): with precision P = 3 and X the exponent of the
// `%e` conversion (after rounding to P digits), `%f` with precision P - 1 - X
// when P > X >= -4, otherwise `%e` with precision P - 1 and an exponent of a
// sign and at least two digits; `#` keeps the decimal point and the trailing
// zeros. Rust's precision formatting rounds the exact binary value half to
// even, as glibc does.
fn format_ncbi_ka_value(value: f64) -> String {
    const PRECISION: i32 = 3;
    if !value.is_finite() {
        let text = if value.is_nan() { "nan" } else { "inf" };
        return if value.is_sign_negative() {
            format!("-{text}")
        } else {
            text.to_string()
        };
    }
    let e_style = format!("{:.*e}", (PRECISION - 1) as usize, value);
    let (mantissa, exponent) = e_style
        .split_once('e')
        .expect("Rust exponent formatting has an 'e'");
    let exponent: i32 = exponent
        .parse()
        .expect("Rust exponent formatting has an integer exponent");
    if (-4..PRECISION).contains(&exponent) {
        let mut formatted = format!("{:.*}", (PRECISION - 1 - exponent) as usize, value);
        if !formatted.contains('.') {
            formatted.push('.');
        }
        formatted
    } else {
        let sign = if exponent < 0 { '-' } else { '+' };
        format!("{mantissa}e{sign}{:02}", exponent.unsigned_abs())
    }
}

#[inline]
fn write_ncbi_ka_field<W: Write>(writer: &mut W, value: f64) -> io::Result<()> {
    write!(writer, "{:>8} ", format_ncbi_ka_value(value))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:391-395
/// The two `NcbiEndl` after the version line of `PrintProlog`'s references: two line feeds,
/// each a flush of NCBI's (`ReportWrite::ncbi_flush`).
fn write_prolog_endl(writer: &mut (impl ReportWrite + ?Sized)) -> io::Result<()> {
    writeln!(writer)?;
    writer.ncbi_flush()?;
    writeln!(writer)?;
    writer.ncbi_flush()
}

/// `CBlastFormat::PrintProlog` of the reports of formats 0-4 of every program (the one place of
/// the start of their output): the version line (`reference::write_version_line`), the `NcbiEndl`
/// pair (`write_prolog_endl`) and the references (`reference::write_prolog_references`: the
/// Gapped BLAST or, for `megablast`, the megablast reference and, with `indexed_megablast`, the
/// database indexing one after it, psiblast's compositionally adjusted matrices and the
/// composition-based statistics of psiblast and blastp). `program` is NCBI's `m_Program` (any
/// case; `blastp_family_program` for blastp's driver).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:360-395
/// NCBI writes no prolog for the XML formats and those after them (they write their own
/// headers); an HTML report starts with the HTML prefix line; then comes the version line, and,
/// unless a bl2seq search without a database scan ends there, the two `NcbiEndl`. (The callers
/// are those of formats 0-4 that write the references.) `disable_ka_stats` is
/// NCBI's `m_DisableKAStats` (blastn's matrix-only scoring, `BlastnIntro::matrix_only`): a
/// blastn search then names `rmblastn` in the version line and writes RMBlast's reference, three
/// lines with no line end after the last, in NCBI's four insertions, before the `NcbiEndl` pair;
/// the references after it are blastn's. `NcbiEndl` is `std::endl`, which flushes in every NCBI
/// program: each is a flush of NCBI's here, in every program and also in an `-html` prolog (line
/// RG item T6a: the flush is where these bytes reach the device and where a failed write shows;
/// line RG T4a2 review B4: blastp's `-html -out /dev/full` with `GETURL_DEBUG` wrote
/// `REGISTRY: .ncbirc` before NCBI's write error). The references and the database lines after
/// the pair stay in NCBI's buffer (`ReportStream`) until a post on `cerr` (tied to `cout`) or its
/// buffer writes them, so `GETURL_DEBUG`'s `REGISTRY:` line, which is not a post
/// (`QueryWarnings::html_stderr_line`), comes before them where both streams go to one pipe
/// (line RG T4a2 fix round 2, review c C4: NCBI BLAST+ 2.17.0 after the version line's 159 bytes
/// in every program, probes `T/t4/fix2/out/c4.*`).
/// `html`: the run's `-html` data (`html_ids::HtmlLinks`; `None`: the plain report) and the run's
/// diagnostics: the HTML head and `<PRE>` (`HTML_PREFIX`) first, the version line in bold, and the
/// references with their links, each with its own read of `./.ncbirc`, whose posts go to the
/// diagnostics (`HtmlLinks::ncbirc`); a malformed file ends the prolog after the version line and
/// its two newlines with NCBI's exception (line RG T4 probes `bp.di.f0.rc_bad.h`,
/// `rc_bad_dontuse.h`: exit status 255).
pub(crate) fn write_prolog<W: ReportWrite + ?Sized>(
    writer: &mut W,
    program: &str,
    version: &str,
    megablast: bool,
    indexed_megablast: bool,
    disable_ka_stats: bool,
    html: Option<(&super::html_ids::HtmlLinks, &mut dyn Write)>,
) -> io::Result<()> {
    if html.is_some() {
        writer.write_all(HTML_PREFIX)?;
    }
    if disable_ka_stats && program.eq_ignore_ascii_case("blastn") {
        super::reference::write_version_line(writer, html.is_some(), "rmblastn", version)?;
        for insertion in [
            &b"\n\n"[..],
            b"Reference: Robert M. Hubley, Arian Smit\n",
            b"RMBlast - RepeatMasker Search Engine\n",
            b"2010 <http://www.repeatmasker.org>",
        ] {
            writer.write_all(insertion)?;
        }
    } else {
        super::reference::write_version_line(writer, html.is_some(), program, version)?;
    }
    write_prolog_endl(writer)?;
    super::reference::write_prolog_references(
        writer,
        html.map(|(links, diagnostics)| (&links.ncbirc, diagnostics)),
        program,
        megablast,
        indexed_megablast,
    )
}

/// The start of an HTML report of formats 0-4 (`kHTML_Prefix`, its four lines, and the newline
/// that `PrintProlog` writes after it), before the version line.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:326-330
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:371-373
const HTML_PREFIX: &[u8] = b"<HTML>\n<HEAD><TITLE>BLAST Search Results</TITLE></HEAD>\n<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" VLINK=\"#660099\" ALINK=\"#660099\">\n<PRE>\n\n";

/// The end of an HTML report of formats 0-4 (`kHTML_Suffix`, its three lines, and the newline that
/// `PrintEpilog` writes after it), after the epilog's last line.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:332-335
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2290-2292
const HTML_SUFFIX: &[u8] = b"</PRE>\n</BODY>\n</HTML>\n";

/// NCBI's `m_Program` of a report of blastp's driver: `psiblast`, `rpsblast`, or `blastp` for
/// any other name (its tasks).
pub(crate) fn blastp_family_program(program: &str) -> &str {
    match program {
        "psiblast" | "rpsblast" => program,
        _ => "blastp",
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:125-138
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:5088-5340 (WrapIt, fWrap_FlatFile)
// Breakpoint priority (`EScore`): forced < punctuation < comma < space < newline;
// keep the rightmost candidate position past the line start with score at least the best seen.
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:288-293
// `NStr::Wrap` measures the `std::string` in bytes, so a line can end inside a UTF-8
// sequence (`\xe6\xbc` and `\xa2` on two lines). Each line and its newline are NCBI's
// insertions (`cli::insert`).
fn write_flatfile_wrapped<W: Write + ?Sized>(
    writer: &mut W,
    text: &[u8],
    width: usize,
) -> io::Result<()> {
    write_flatfile_wrapped_with(writer, text, width, TextFacts::of(text))
}

/// Whether a text holds a newline or a backspace, which `write_flatfile_wrapped_with` needs to know (a text
/// without either has a column for each byte and no line of its own).
#[derive(Clone, Copy)]
pub(crate) struct TextFacts {
    /// The text holds a `\n`.
    pub(crate) newlines: bool,
    /// The text holds a backspace (8).
    pub(crate) backspaces: bool,
}

impl TextFacts {
    /// The facts of `text`.
    pub(crate) fn of(text: &[u8]) -> Self {
        TextFacts {
            newlines: text.contains(&b'\n'),
            backspaces: text.contains(&8),
        }
    }
}

/// `write_flatfile_wrapped` of a text whose facts (`TextFacts`) are `facts`: the caller of a
/// constant text knows them (`report::reference`).
pub(crate) fn write_flatfile_wrapped_with<W: Write + ?Sized>(
    writer: &mut W,
    text: &[u8],
    width: usize,
    facts: TextFacts,
) -> io::Result<()> {
    let bytes = text;
    let mut pos = 0;
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:5100,5153-5157,5137
    let mut newline = 0;
    while pos < bytes.len() {
        let mut has_backspace = false;
        if newline <= pos {
            newline = if facts.newlines {
                bytes[pos..]
                    .iter()
                    .position(|&c| c == b'\n')
                    .map_or(bytes.len(), |n| pos + n)
            } else {
                bytes.len()
            };
        }
        let (best, best_score) = if !facts.backspaces && width > 0 && newline - pos > width {
            // The scan below on a text without backspaces, from a line start whose next newline
            // is past `width` columns: it reads the columns `pos..=pos + width`, which hold no
            // newline, and ends past them (`column > width`); its break is the last of the
            // highest score, which `last_break` finds from the end.
            last_break(bytes, pos, pos + width)
        } else {
            let pos0 = if newline - pos <= width { newline } else { pos };
            let mut scan = pos0;
            let mut column = 0;
            let mut best = pos;
            let mut best_score = 0;
            while scan < bytes.len() && column <= width {
                let c = bytes[scan];
                let mut score = 0;
                let mut score_pos = scan;
                if c == b'\n' {
                    best = scan;
                    best_score = 4;
                    break;
                }
                if c.is_ascii_whitespace() || c == 11 {
                    score = 3;
                } else if c == b',' && column < width && scan + 1 < bytes.len() {
                    score = 2;
                    score_pos += 1;
                } else if c == b'-' && column < width && scan + 1 < bytes.len() {
                    score = 1;
                    score_pos += 1;
                }
                if score >= best_score && score_pos > pos0 {
                    best = score_pos;
                    best_score = score;
                }
                while scan + 1 < bytes.len() && bytes[scan + 1] == 8 {
                    scan += 1;
                    column = column.saturating_sub(1);
                    has_backspace = true;
                }
                scan += 1;
                column += 1;
            }
            // NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:5260-5266,5277-5300
            // When the remaining bytes fit on the line, extend the break to the end and mark
            // backspace processing so backspaces erase their preceding characters.
            if best_score != 4 && column <= width && best != bytes.len() {
                best = bytes.len();
                has_backspace = true;
            }
            (best, best_score)
        };
        // Eating backspaces leaves a line without them as it is.
        if has_backspace && facts.backspaces {
            let mut line = Vec::with_capacity(best - pos);
            for &c in &bytes[pos..best] {
                if c == 8 {
                    line.pop();
                } else {
                    line.push(c);
                }
            }
            crate::cli::insert(writer, &line)?;
        } else {
            crate::cli::insert(writer, &bytes[pos..best])?;
        }
        crate::cli::insert(writer, b"\n")?;
        pos = best;
        if best_score == 3 {
            while pos < bytes.len() && bytes[pos] == b' ' {
                pos += 1;
            }
            if pos < bytes.len() && bytes[pos] == b'\n' {
                pos += 1;
            }
        }
        if best_score == 4 {
            pos += 1;
        }
        while pos < bytes.len() && bytes[pos] == 8 {
            pos += 1;
        }
    }
    Ok(())
}

/// The break of the line from `pos` of `write_flatfile_wrapped_with`'s scan where its columns
/// `pos..=end` hold no newline and no backspace and the scan ends past `end`: the last of the
/// highest score, the running maximum of `score >= best_score && score_pos > pos0` (`pos0` is
/// `pos`): a blank after `pos` (score 3, the line ends before it), else a `,` or else a `-` before
/// `end` (scores 2 and 1, the line ends after it), else `end` itself (score 0: every column after
/// `pos` scores 0, `end` last).
fn last_break(bytes: &[u8], pos: usize, end: usize) -> (usize, i32) {
    let blank = |c: &u8| c.is_ascii_whitespace() || *c == 11;
    if let Some(at) = bytes[pos + 1..=end].iter().rposition(blank) {
        return (pos + 1 + at, 3);
    }
    if let Some(at) = bytes[pos..end].iter().rposition(|&c| c == b',') {
        return (pos + at + 1, 2);
    }
    if let Some(at) = bytes[pos..end].iter().rposition(|&c| c == b'-') {
        return (pos + at + 1, 1);
    }
    (end, 0)
}

/// The outfmt 0 prolog of blastp (the version, the references and the database), which
/// NCBI writes before it reads the first query batch (c++/src/algo/blast/format/blast_format.cpp, `PrintProlog`).
pub fn write_blastp_pairwise_prolog<W: ReportWrite>(
    writer: &mut W,
    version: &str,
    program: &str,
    database_name: &[u8],
    database_num_sequences: i64,
    database_total_letters: i64,
    database: Option<&PairwiseDatabase>,
    html: Option<(&super::html_ids::HtmlLinks, &mut dyn Write)>,
) -> io::Result<()> {
    write_prolog(
        writer,
        blastp_family_program(program),
        version,
        false,
        false,
        false,
        html,
    )?;
    if let Some(database) = database {
        return write_blast_database_header(writer, database);
    }
    // The blank lines before the first query are the query's (as TBLASTX's prolog).
    write_blastp_database_header_spacing(
        writer,
        database_name,
        database_num_sequences,
        database_total_letters,
        1,
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:511-538
/// The `Database:` lines of the prolog for the database of a `-db` search: the titles of the
/// databases of its list joined with `; ` (`PairwiseDatabase::new`), wrapped at 68, then the
/// summed counts and one blank line.
fn write_blast_database_header<W: ReportWrite>(
    writer: &mut W,
    database: &PairwiseDatabase,
) -> io::Result<()> {
    write!(writer, "Database: ")?;
    write_flatfile_wrapped(writer, &database.definition, 68)?;
    if let Some(first) = database.databases.first() {
        write_masked_using(writer, "", first)?;
    }
    writeln!(
        writer,
        "           {} sequences; {} total letters",
        crate::ncbi_str::int8_with_commas(database.num_sequences),
        crate::ncbi_str::int8_with_commas(database.total_letters)
    )?;
    writeln!(writer)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:525-553
/// The `Masked using:` line of a database searched with database masking (`indent` before it:
/// none in the prolog, two spaces in the epilog), not wrapped; nothing without it.
fn write_masked_using<W: ReportWrite>(
    writer: &mut W,
    indent: &str,
    database: &DatabaseInfo,
) -> io::Result<()> {
    if database.filt_algorithm_name.is_empty() {
        return Ok(());
    }
    write!(
        writer,
        "{indent}Masked using: '{}'",
        database.filt_algorithm_name
    )?;
    if !database.filt_algorithm_options.is_empty() {
        write!(writer, ", options: '{}'", database.filt_algorithm_options)?;
    }
    // An insertion of `endl` (align_format_util.cpp:531, :551): a flush.
    crate::report::ostream::ncbi_endl(writer)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:541-574
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2255-2267
/// The database blocks of the epilog, one per database of the list, and the blank lines
/// before the matrix line (as `write_final_database_report`).
fn write_final_database_block<W: ReportWrite>(
    writer: &mut W,
    database: &PairwiseDatabase,
) -> io::Result<()> {
    for database in &database.databases {
        write!(writer, "  Database: ")?;
        write_flatfile_wrapped(writer, &database.definition, 68)?;
        write_masked_using(writer, "  ", database)?;
        writeln!(writer, "    Posted date:  {}", database.date)?;
        writeln!(
            writer,
            "  Number of letters in database: {}",
            crate::ncbi_str::int8_with_commas(database.total_letters)
        )?;
        writeln!(
            writer,
            "  Number of sequences in database:  {}",
            crate::ncbi_str::int_with_commas(database.num_sequences)
        )?;
        writeln!(writer)?;
    }
    writeln!(writer)?;
    writeln!(writer)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:533-538
fn write_blastp_database_header_spacing<W: Write>(
    writer: &mut W,
    database_name: &[u8],
    database_num_sequences: i64,
    database_total_letters: i64,
    trailing: usize,
) -> io::Result<()> {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:513-525
    write!(writer, "Database: ")?;
    write_flatfile_wrapped(writer, &ensure_trailing_period(database_name), 68)?;
    writeln!(
        writer,
        "           {} sequences; {} total letters",
        crate::ncbi_str::int8_with_commas(database_num_sequences),
        crate::ncbi_str::int8_with_commas(database_total_letters)
    )?;
    for _ in 0..trailing {
        writeln!(writer)?;
    }
    Ok(())
}

/// The `Query=` block of a query of title `query_title` and length `query_length`, as NCBI's
/// insertions (`cli::insert`; the organism report's stream counts them, outfmt 0's buffers them).
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:729
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:741-745
pub(crate) fn write_blastp_query_header<W: Write + ?Sized>(
    writer: &mut W,
    query_title: &[u8],
    query_length: usize,
    html: bool,
) -> io::Result<()> {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:732-734
    // The ID of a local query is not shown (`believe_query` false), so the line is the
    // title, raw (no decoding), without the white space at its ends.
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:724-731
    // `-html` (D4 of the T4 survey): the label in bold and the line HTML-encoded and wrapped
    // with `fWrap_HTMLPre` (`x_WrapOutputLine`, `report::html::write_html_wrapped`).
    if html {
        crate::cli::insert(writer, b"<b>Query=</b> ")?;
        super::html::write_html_wrapped(writer, truncate_spaces(query_title), 68)?;
    } else {
        crate::cli::insert(writer, b"Query")?;
        crate::cli::insert(writer, b"= ")?;
        write_flatfile_wrapped(writer, truncate_spaces(query_title), 68)?;
    }
    crate::cli::insert(writer, b"\nLength=")?;
    crate::cli::insert(writer, query_length.to_string().as_bytes())?;
    crate::cli::insert(writer, b"\n")
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1510-1514
fn write_no_hits_found<W: Write>(writer: &mut W) -> io::Result<()> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1510-1514
    // Write `kNoHitsFound` between `"***** "` and `" *****"`, after two newlines and followed by three newlines.
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer, "***** No hits found *****")?;
    writeln!(writer)?;
    writeln!(writer)?;
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:445-477
// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180,204-208:
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:451-478 and
// c++/src/objtools/align_format/align_format_util.cpp:581-613:
// PrintKAParameters emits both -1 Karlin blocks with no Gumbel columns.
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:477-478
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1491-1491
fn write_tblastn_unsearched_query_footer_spacing<W: crate::report::ostream::ReportWrite>(
    writer: &mut W,
    trailing: bool,
) -> io::Result<()> {
    // `x_PrintOneQueryFooter`'s first `NcbiEndl` flushes the stream.
    crate::report::ostream::ncbi_endl(writer)?;
    writeln!(writer, "Lambda      K        H")?;
    for _ in 0..3 {
        write_ncbi_ka_field(writer, -1.0)?;
    }
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer, "Gapped")?;
    writeln!(writer, "Lambda      K        H")?;
    for _ in 0..3 {
        write_ncbi_ka_field(writer, -1.0)?;
    }
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer, "Effective search space used: 0")?;
    if trailing {
        writeln!(writer)?;
        writeln!(writer)?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:445-478:
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:477-478
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1491-1491
fn write_blastp_query_footer_spacing<W: crate::report::ostream::ReportWrite>(
    writer: &mut W,
    ungapped_karlin: KarlinParams,
    gapped_karlin: KarlinParams,
    gumbel: Option<BlastGumbelBlk>,
    gapped: bool,
    effective_search_space: i64,
    trailing: bool,
) -> io::Result<()> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:445-478 (cited
    // above): `kbp_gap` (the `Gapped` block) is printed only where the search made it; an
    // ungapped search has neither it nor a Gumbel block.
    let Some(gumbel) = gumbel.filter(|_| gapped) else {
        return write_nucleotide_query_footer_spacing(
            writer,
            Some((ungapped_karlin, gapped_karlin)),
            gapped,
            effective_search_space,
            trailing,
        );
    };
    // `x_PrintOneQueryFooter`'s first `NcbiEndl` flushes the stream.
    crate::report::ostream::ncbi_endl(writer)?;
    // A block that the query's results do not keep (`stats::tables::ancillary_kept`) has no
    // part in the footer.
    if crate::stats::tables::ancillary_kept(ungapped_karlin.lambda) {
        writeln!(writer, "Lambda      K        H        a         alpha")?;
        write_ncbi_ka_field(writer, ungapped_karlin.lambda)?;
        write_ncbi_ka_field(writer, ungapped_karlin.k)?;
        write_ncbi_ka_field(writer, ungapped_karlin.h)?;
        write_ncbi_ka_field(writer, gumbel.a_un)?;
        write_ncbi_ka_field(writer, gumbel.alpha_un)?;
        writeln!(writer)?;
    }
    writeln!(writer)?;
    if crate::stats::tables::ancillary_kept(gapped_karlin.lambda) {
        writeln!(writer, "Gapped")?;
        writeln!(
            writer,
            "Lambda      K        H        a         alpha    sigma"
        )?;
        write_ncbi_ka_field(writer, gapped_karlin.lambda)?;
        write_ncbi_ka_field(writer, gapped_karlin.k)?;
        write_ncbi_ka_field(writer, gapped_karlin.h)?;
        write_ncbi_ka_field(writer, gumbel.a)?;
        write_ncbi_ka_field(writer, gumbel.alpha)?;
        write_ncbi_ka_field(writer, gumbel.sigma)?;
        writeln!(writer)?;
    }
    writeln!(writer)?;
    writeln!(
        writer,
        "Effective search space used: {}",
        effective_search_space
    )?;
    if trailing {
        writeln!(writer)?;
        writeln!(writer)?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2266-2268
// Observed, not taken from a library's source: inserting a null `const char*` into NCBI's
// stream sets `badbit` (see `write_matrix_line`).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:96-103
// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:73-83 (cited at
// `rps_search::volumes::search_sets`: `CreateRpsStructures` names the database's matrix in the
// options of each RPS query batch it sets up)
/// The epilog's matrix line after its blank lines (`options.GetMatrixName()`). `None`: the
/// options name no matrix, as RPS-BLAST's (`CBlastRPSOptionsHandle` sets none) before a query
/// batch's set-up names the database's; the null `const char*` sets `badbit` on the formatter's
/// stream, which throws (`ReportWrite::badbit`): the report ends after `Matrix: ` (NCBI BLAST+
/// 2.17.0: rpsblast and rpstblastn with an empty pipe, outfmt 0 to 4, exit status 6).
fn write_matrix_line(writer: &mut impl ReportWrite, matrix_name: Option<&str>) -> io::Result<()> {
    writer.write_all(b"Matrix: ")?;
    let Some(name) = matrix_name else {
        return writer.badbit();
    };
    writeln!(writer, "{name}")
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:503-579
//
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2258-2285
fn write_blastp_final_footer<W: ReportWrite>(
    writer: &mut W,
    report: &BlastpPairwiseReport,
    database: Option<&PairwiseDatabase>,
    html: bool,
) -> io::Result<()> {
    match database {
        Some(database) => write_final_database_block(writer, database)?,
        None => write_final_database_report(
            writer,
            &report.database_name,
            report.database_num_sequences,
            report.database_total_letters,
        )?,
    }
    write_matrix_line(writer, report.matrix_name.as_deref())?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2270-2288 (cited at
    // `write_blastn_final_footer`): with `options.GetGappedMode()` true it writes "Gap Penalties: ".
    if report.gapped {
        writeln!(
            writer,
            "Gap Penalties: Existence: {}, Extension: {}",
            report.gap_open, report.gap_extend
        )?;
    }
    write_epilog_end(
        writer,
        report.word_threshold,
        report.window_size as i64,
        html,
    )
}

/// The end of `CBlastFormat::PrintEpilog` of formats 0-4, the one place of the end of their output
/// for every program: the word threshold's and the window's lines, each when it is not 0
/// (BLASTN's threshold is 0 for every task, BLAST_WORD_THRESHOLD_BLASTN and
/// BLAST_WORD_THRESHOLD_MEGABLAST, blast_options.h). `GetWordThreshold()` is a double, written
/// with the stream's default format (`cpp_default_double`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2281-2292
/// The two lines, then, for an HTML report (`html`, `-html`), the end of its `<PRE>`, body and
/// document (`HTML_SUFFIX`).
fn write_epilog_end<W: Write + ?Sized>(
    writer: &mut W,
    word_threshold: f64,
    window_size: i64,
    html: bool,
) -> io::Result<()> {
    if word_threshold != 0.0 {
        writeln!(
            writer,
            "Neighboring words threshold: {}",
            cpp_default_double(word_threshold)
        )?;
    }
    if window_size != 0 {
        writeln!(writer, "Window for multiple hits: {window_size}")?;
    }
    if html {
        writer.write_all(HTML_SUFFIX)?;
    }
    Ok(())
}

/// The database block of the epilog and the blank lines before the matrix line
/// (align_format_util.cpp:503-579 and blast_format.cpp:2258-2262, cited at
/// `write_blastp_final_footer`).
fn write_final_database_report<W: Write>(
    writer: &mut W,
    database_name: &[u8],
    database_num_sequences: i64,
    database_total_letters: i64,
) -> io::Result<()> {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:543-544
    write!(writer, "  Database: ")?;
    write_flatfile_wrapped(writer, &ensure_trailing_period(database_name), 68)?;
    writeln!(writer, "    Posted date:  Unknown")?;
    writeln!(
        writer,
        "  Number of letters in database: {}",
        crate::ncbi_str::int8_with_commas(database_total_letters)
    )?;
    writeln!(
        writer,
        "  Number of sequences in database:  {}",
        crate::ncbi_str::int_with_commas(i32::try_from(database_num_sequences).unwrap_or(0))
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:387-441
//
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1520-1589
/// The report of `queries`: the prolog when `prolog` is set (EXPERIMENT (LOSAT_X_BPQBLOCK):
/// false for the query blocks after the first of blastp's block search, whose first block wrote
/// it), each query's section, and the epilog when `epilog` is set. `query_gapped_k` is
/// rpsblast's gapped K of each query's footer (`rps_search::engine::batch_gapped_k`), `None` for
/// blastp and psiblast (the report's `gapped_karlin`).
pub fn write_blastp_pairwise_report<W: crate::report::ostream::ReportWrite>(
    hits: &[PairwiseHit],
    writer: &mut W,
    config: &PairwiseConfig,
    queries: &[BlastpPairwiseQuery],
    query_titles: &[Arc<[u8]>],
    subject_titles: &[Arc<[u8]>],
    parsed_subjects: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    subject_sources: &[super::defline::SourceTaxname],
    report: &BlastpPairwiseReport,
    query_gapped_k: Option<&[f64]>,
    database: Option<&PairwiseDatabase>,
    probe: Option<&mut FormatProbe<'_>>,
    warnings: Option<&mut super::query_warnings::QueryWarnings<'_>>,
    prolog: bool,
    epilog: bool,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
    loads: &crate::report::db_subject::HitLoads,
) -> io::Result<()> {
    write_blastp_pairwise_window(
        hits,
        0..queries.len(),
        writer,
        config,
        queries,
        query_titles,
        subject_titles,
        parsed_subjects,
        subject_sources,
        report,
        query_gapped_k,
        database,
        probe,
        warnings,
        prolog,
        epilog,
        x_pool,
        None,
        None,
        loads,
    )
}

/// The queries `query_range` of a BLASTP pairwise report, whose `hits` are those queries' (the
/// rows of a caller that makes them a window of queries at a time, as BLASTN's and TBLASTX's
/// do, `x_pairwise_par::query_windows`): the prolog with the report's first query, the epilog
/// with its last. `anchored`: the sequences of the query-anchored sections (`config.view` of
/// outfmt 1 to 4; `None` for outfmt 0).
#[allow(clippy::too_many_arguments)]
pub fn write_blastp_pairwise_window<W: crate::report::ostream::ReportWrite>(
    hits: &[PairwiseHit],
    query_range: std::ops::Range<usize>,
    writer: &mut W,
    config: &PairwiseConfig,
    queries: &[BlastpPairwiseQuery],
    query_titles: &[Arc<[u8]>],
    subject_titles: &[Arc<[u8]>],
    parsed_subjects: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    subject_sources: &[super::defline::SourceTaxname],
    report: &BlastpPairwiseReport,
    query_gapped_k: Option<&[f64]>,
    database: Option<&PairwiseDatabase>,
    probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut super::query_warnings::QueryWarnings<'_>>,
    prolog: bool,
    epilog: bool,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
    anchored: Option<&dyn super::query_anchored::ReportSequences>,
    // The rows' masks of `hits` under `-html` (`HtmlRowMasks`; `None` without).
    row_masks: Option<&[super::query_anchored::HspMasks]>,
    loads: &crate::report::db_subject::HitLoads,
) -> io::Result<()> {
    let mut buffered = io::BufWriter::new(writer);
    let writer = &mut buffered;

    if prolog && query_range.start == 0 {
        write_prolog(
            writer,
            blastp_family_program(&config.program),
            &report.version,
            false,
            false,
            false,
            html_with_diagnostics(config.html_links(), warnings.as_deref_mut())?,
        )?;
        // The prolog ends with one blank line; each query's preamble starts with two
        // (blast_format.cpp:1491), after the query's warnings.
        match database {
            Some(database) => write_blast_database_header(writer, database)?,
            None => write_blastp_database_header_spacing(
                writer,
                &report.database_name,
                report.database_num_sequences,
                report.database_total_letters,
                1,
            )?,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411
    // Each hit keeps its index in `hits` (the final result order), which identifies it
    // to the probe.
    let mut hits_by_query: Vec<Vec<(HspIndex, &PairwiseHit)>> = vec![Vec::new(); query_range.len()];
    for (hsp_index, hit) in hits.iter().enumerate() {
        if let Some(bucket) = (hit.hit.q_idx as usize)
            .checked_sub(query_range.start)
            .and_then(|offset| hits_by_query.get_mut(offset))
        {
            bucket.push((hsp_index, hit));
        }
    }

    // The queries' reports, on `x_pool` with LOSAT_X_PAIRWISEPAR (`x_pairwise_par`).
    super::x_pairwise_par::write(
        writer,
        &BlastpSections {
            config,
            queries,
            query_titles,
            report,
            query_gapped_k,
            database,
            parsed_subjects,
            subject_sources,
            anchored,
            row_masks: HtmlRowMasks::new(row_masks, 0),
            loads,
        },
        query_range.start,
        &hits_by_query,
        subject_titles,
        probe,
        warnings.as_deref_mut(),
        x_pool,
    )?;
    // The warnings after the last query's report (`QueryWarnings::after_last_query`).
    if let Some(warnings) = warnings.filter(|_| query_range.end == queries.len()) {
        warnings.after_last_query(queries.len(), writer)?;
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:294-295
    // A search that stops at a query batch (`Empty CBlastQueryVector` with `-query_loc`)
    // writes no epilog.
    if !(epilog && query_range.end == queries.len()) {
        return writer.flush();
    }

    // The blank lines before the epilog (blast_format.cpp:2249): two insertions of `NcbiEndl`,
    // each a flush.
    crate::report::ostream::ncbi_endl(writer)?;
    crate::report::ostream::ncbi_endl(writer)?;
    write_blastp_final_footer(writer, report, database, config.html.is_some())?;
    writer.flush()
}

/// The parts of one query's BLASTP pairwise report (`x_pairwise_par::QuerySections`).
struct BlastpSections<'a> {
    config: &'a PairwiseConfig,
    queries: &'a [BlastpPairwiseQuery],
    query_titles: &'a [Arc<[u8]>],
    report: &'a BlastpPairwiseReport,
    /// rpsblast's gapped K by query (`write_blastp_pairwise_report`).
    query_gapped_k: Option<&'a [f64]>,
    database: Option<&'a PairwiseDatabase>,
    /// The subjects with parsed Seq-ids by subject (empty without them).
    parsed_subjects: &'a [Option<Arc<crate::report::db_subject::ParsedSubject>>],
    /// The sources of the others' Bioseqs by subject (`stored_source`; empty without them).
    subject_sources: &'a [super::defline::SourceTaxname],
    /// The sequences of the query-anchored sections (outfmt 1 to 4).
    anchored: Option<&'a dyn super::query_anchored::ReportSequences>,
    /// The rows' masks of the HSPs under `-html` (outfmt 0).
    row_masks: HtmlRowMasks<'a>,
    /// What the reports' loads do at each subject, by subject (`thrown_load`).
    loads: &'a crate::report::db_subject::HitLoads,
}

/// The HSPs of a query-anchored section (outfmt 1 to 4), added to `hsps`: the first
/// `num_alignments` subjects of the query's sorted hits (`QueryGroup::aligned`), each with all its
/// HSPs in their order. An HSP's Dense-seg is its Seq-align's, which `segments_of` makes and gives
/// its rows' strands and query frame of (`seqalign::protein_row_segments`,
/// `nucleotide_row_segments`, `translated_row_segments`); the subjects shown in one hit (copies of
/// a database record, `joins`) are the hit's subject, one Bioseq for NCBI's scope; `is_query`
/// tells the subject that is the query's own Bioseq. A subject whose Bioseq NCBI's scope does not
/// get (`parsed`, by subject index: `ParsedSubject::scope_finds_none`) is left out, and counted
/// among the first `num_alignments`.
fn section_hsps(
    mut hsps: super::query_anchored::Hsps,
    num_alignments: usize,
    group: &super::x_pairwise_par::QueryGroup<'_>,
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    is_query: impl Fn(u32) -> bool,
    mut segments_of: impl FnMut(
        &PairwiseHit,
        &mut super::serial::seqalign::Segments,
    ) -> super::query_anchored::HspStrands,
) -> super::query_anchored::Hsps {
    let mut segments = super::serial::seqalign::Segments::default();
    for &s_idx in group.aligned(num_alignments) {
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1997-2003
        if finds_no_bioseq(parsed, s_idx) {
            continue;
        }
        if is_query(s_idx) {
            hsps.set_query_bioseq(s_idx);
        }
        for hit in &group.subject_hits[&s_idx] {
            let strands = segments_of(hit, &mut segments);
            hsps.push_dense_seg(s_idx, segments.starts(), segments.lens(), 0, 0, strands);
        }
    }
    hsps
}

/// `section_hsps` of a BLASTP section (`seqalign::protein_row_segments`: no strands).
fn blastp_section_hsps(
    num_alignments: usize,
    group: &super::x_pairwise_par::QueryGroup<'_>,
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    is_query: impl Fn(u32) -> bool,
) -> super::query_anchored::Hsps {
    section_hsps(
        super::query_anchored::Hsps::default(),
        num_alignments,
        group,
        parsed,
        is_query,
        |hit, segments| blastp_row_segments(&hit.hit, segments),
    )
}

/// The Dense-seg segments of a BLASTP, RPSBLAST, PSIBLAST or DELTABLAST HSP, for its rows in
/// outfmt 0 (`write_alignment_with_sequences`) and 1 to 4 (`blastp_section_hsps`):
/// `seqalign::protein_row_segments`, whose rows have no strands.
fn blastp_row_segments(
    hit: &Hit,
    segments: &mut super::serial::seqalign::Segments,
) -> super::query_anchored::HspStrands {
    super::serial::seqalign::protein_row_segments(hit, segments);
    super::query_anchored::HspStrands::default()
}

/// `section_hsps` of a BLASTN section (`seqalign::nucleotide_row_segments`).
pub(crate) fn blastn_section_hsps(
    num_alignments: usize,
    group: &super::x_pairwise_par::QueryGroup<'_>,
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    is_query: impl Fn(u32) -> bool,
) -> super::query_anchored::Hsps {
    section_hsps(
        super::query_anchored::Hsps::default(),
        num_alignments,
        group,
        parsed,
        is_query,
        |hit, segments| match hit.left_of_record() {
            Some(left) if hit.hit.gap_info.is_none() => {
                super::serial::seqalign::nucleotide_left_of_record_segments(
                    &hit.hit, left, segments,
                )
            }
            _ => super::serial::seqalign::nucleotide_row_segments(&hit.hit, segments),
        },
    )
}

/// `section_hsps` of a BLASTX section (`seqalign::translated_row_segments`: a nucleotide query row
/// of width 3 in the hit's frame of the query's record, a protein subject row; an HSP of length 0,
/// whose rows the writer reaches as the engine made them, `blastx_mapped_empty_rows`: one segment
/// of length 0 without widths, `seqalign::empty_row_segments`).
fn blastx_section_hsps(
    num_alignments: usize,
    group: &super::x_pairwise_par::QueryGroup<'_>,
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
) -> super::query_anchored::Hsps {
    section_hsps(
        super::query_anchored::Hsps::of_translated_query(),
        num_alignments,
        group,
        parsed,
        |_| false,
        blastx_row_segments,
    )
}

/// The Dense-seg segments of a BLASTX HSP, for its rows in outfmt 0 (`write_blastx_alignment`)
/// and 1 to 4 (`blastx_section_hsps`): `seqalign::translated_row_segments` of a nucleotide query
/// row of width 3 in the hit's frame of the query's record and a protein subject row; an HSP of
/// length 0, whose rows the writer reaches as the engine made them (`blastx_mapped_empty_rows`):
/// one segment of length 0 without widths (`seqalign::empty_row_segments`).
fn blastx_row_segments(
    hit: &PairwiseHit,
    segments: &mut super::serial::seqalign::Segments,
) -> super::query_anchored::HspStrands {
    let frame = hit.query_frame.unwrap_or(1);
    if hit.rows_empty() {
        super::serial::seqalign::empty_row_segments(&hit.hit, frame, segments)
    } else {
        super::serial::seqalign::translated_row_segments(
            &hit.hit,
            Some((frame, hit.hit.query_length)),
            None,
            segments,
        )
    }
}

/// The Dense-seg segments of a TBLASTN HSP, for its rows in outfmt 0 (`write_tblastn_alignment`)
/// and 1 to 4 (`TblastnSections::alignments`): `seqalign::translated_row_segments` of a protein
/// query row and a nucleotide subject row of width 3 in the hit's frame of the subject's record.
fn tblastn_row_segments(
    hit: &PairwiseHit,
    segments: &mut super::serial::seqalign::Segments,
) -> super::query_anchored::HspStrands {
    if hit.edge == PairwiseEdge::NullbEmptySubjectRow {
        return super::serial::seqalign::empty_subject_row_segments(
            &hit.hit,
            hit.subject_frame.unwrap_or(1),
            segments,
        );
    }
    super::serial::seqalign::translated_row_segments(
        &hit.hit,
        None,
        Some((
            hit.subject_frame.unwrap_or(1),
            hit.subject_length.unwrap_or(0),
        )),
        segments,
    )
}

/// The Dense-seg segments of a TBLASTX HSP of a query whose record has `query_length`
/// nucleotides, for its rows in outfmt 0 (`write_tblastx_alignment`) and 1 to 4
/// (`TblastxSections::alignments`): `seqalign::translated_row_segments` of two nucleotide rows of
/// width 3 in the hit's frames of their records.
fn tblastx_row_segments(
    hit: &PairwiseHit,
    query_length: usize,
    segments: &mut super::serial::seqalign::Segments,
) -> super::query_anchored::HspStrands {
    super::serial::seqalign::translated_row_segments(
        &hit.hit,
        Some((hit.query_frame.unwrap_or(1), query_length)),
        Some((
            hit.subject_frame.unwrap_or(1),
            hit.subject_length.unwrap_or(0),
        )),
        segments,
    )
}

/// The exception at which the description table of the subjects `descriptions` (their hits in
/// `grouped`; `remaps`: `-subject`) stops: it reads the query and then the subject range of every
/// HSP of each described subject before its first line (`x_InitDeflineTable` ->
/// `GetSeqAlignSetCalcParams` -> `s_ProcessAlignSet`), which throw on a row whose range they cannot
/// take (`seen_row_range`): `Row is empty` on an HSP with empty rows, `Row seqids not consistent` on
/// one that starts with a gap (RPSTBLASTN's `-use_sw_tback`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:2913-2921
fn described_row_range(
    remaps: bool,
    descriptions: &[u32],
    grouped: &std::collections::HashMap<u32, Vec<&PairwiseHit>>,
) -> Option<RowRangeError> {
    descriptions.iter().find_map(|oid| {
        let hits = &grouped[oid];
        hits.iter().find_map(|hit| {
            seen_row_range(remaps, hits, hit, AlignRow::Query)
                .or_else(|| seen_row_range(remaps, hits, hit, AlignRow::Subject))
        })
    })
}

/// Whether a BLASTX subject's HSPs `hits` hold one with empty rows whose Std-seg `-subject`'s
/// Seq-loc mapper (`remaps`) rewrote (`mapper_rewrites`): its query row of one nucleotide against
/// the empty subject row, on which `CreateDensegFromStdseg` throws (`EmptyRowsStop::StdSegLengths`).
fn blastx_mapped_empty_rows(remaps: bool, hits: &[&PairwiseHit]) -> bool {
    hits.iter().any(|hit| hit.rows_empty()) && mapper_rewrites(remaps, hits)
}

/// Where NCBI BLAST+'s report of a BLASTX query whose hits are `group` (`BlastxSections`' group) in
/// the query-anchored view of `config` stops at an HSP with empty rows, in the order it meets them:
/// its sorts (`sort_stop`), its description table (`described_row_range`), the Dense-segs of its
/// aligned subjects (`blastx_mapped_empty_rows`, with `-subject`, `remaps`), all of which throw
/// (`write_blastx_pairwise_report` writes NCBI's report up to there and its exception), and the
/// merge of its section, which frees a segment that it still reads
/// (`EmptyRowsStop::MergeFreesSegment`, `query_anchored::section_frees`: the split of a segment at
/// an HSP of length 0 at a codon of an alignment in the query's reading frame +1, or
/// `FillUnalignedRegions` after a segment of length 0). `None` for a report without such an
/// HSP among its aligned subjects. What NCBI writes before the crash is what its stream had flushed
/// (no prefix of it is NCBI's, as for `query_anchored::ncbi_aborts`): the program refuses the
/// report before the call's first byte (`UNPORTED_ANCHORED_FREE`), unless a query
/// before it ends the run with an exception.
fn blastx_anchored_stop(
    config: &PairwiseConfig,
    options: &BlastxPairwiseOptions,
    group: &super::x_pairwise_par::QueryGroup<'_>,
    remaps: bool,
) -> io::Result<Option<EmptyRowsStop>> {
    let (grouped, order) = (&group.subject_hits, &group.subject_order);
    let aligned = group.aligned(options.num_alignments);
    if let Some(stop) = sort_stop(config.sort, remaps, order, grouped) {
        return Ok(Some(stop));
    }
    if let Some(error) = described_row_range(
        remaps,
        &order[..order.len().min(options.num_descriptions)],
        grouped,
    ) {
        return Ok(Some(EmptyRowsStop::RowRange(error)));
    }
    if aligned
        .iter()
        .any(|s_idx| blastx_mapped_empty_rows(remaps, &grouped[s_idx]))
    {
        return Ok(Some(EmptyRowsStop::StdSegLengths));
    }
    let frees = aligned
        .iter()
        .any(|s_idx| grouped[s_idx].iter().any(|hit| hit.rows_empty()))
        && super::query_anchored::section_frees(&blastx_section_hsps(
            options.num_alignments,
            group,
            &[],
        ))?;
    Ok(frees.then_some(EmptyRowsStop::MergeFreesSegment))
}

/// Whether NCBI BLAST+ aborts on the query-anchored section (`view`) of query `q_idx` of a
/// BLASTP `-subject` search (`section_aborts` of `blastp_section_aborts`).
pub(crate) fn blastp_section_aborts(
    view: PairwiseView,
    sort: super::hit_sort::ReportSort,
    num_alignments: usize,
    q_idx: usize,
    rows: &[PairwiseHit],
    subject_titles: &[Arc<[u8]>],
    is_query: impl Fn(u32) -> bool,
) -> io::Result<bool> {
    section_aborts(
        view,
        sort,
        num_alignments,
        q_idx,
        rows,
        subject_titles,
        &is_query,
        |group| blastp_section_hsps(num_alignments, group, &[], &is_query),
    )
}

/// Whether the merge of the query-anchored section (`view`) of query `q_idx` of a BLASTP
/// `-ungapped` search frees a segment that it still reads (`query_anchored::section_frees`, which
/// NCBI BLAST+ crashes on): an HSP of length 0 among the section's aligned subjects, `rows` the
/// rows of the call's hits in the final order and `joins` the copies of one database record
/// (`blastp_report_group`).
pub(crate) fn blastp_section_frees(
    view: PairwiseView,
    sort: super::hit_sort::ReportSort,
    num_alignments: usize,
    q_idx: usize,
    rows: &[PairwiseHit],
    subject_titles: &[Arc<[u8]>],
    joins: impl Fn(u32, u32) -> bool,
) -> io::Result<bool> {
    if view == PairwiseView::Pairwise {
        return Ok(false);
    }
    let Some(group) =
        blastp_report_group(view.report_sort(sort), q_idx, rows, subject_titles, joins)
    else {
        return Ok(false);
    };
    super::query_anchored::section_frees(&blastp_section_hsps(num_alignments, &group, &[], |_| {
        false
    }))
}

/// Whether NCBI BLAST+ aborts on the query-anchored section (`view`) of query `q_idx` of a
/// `-subject` search (`query_anchored::ncbi_aborts`), `rows` the rows of its hits in the final
/// order: the section of the report (`QuerySections::alignments`), whose group is sorted by
/// `sort` (`PairwiseView::report_sort`), whose aligned subjects are the first `num_alignments`
/// and whose HSPs `hsps_of` makes from it (the program's `section_hsps`). Asked before any format of the
/// search's call is written, for a query whose own record is one of the subjects (`is_query`,
/// `-parse_deflines`), so that the refusal comes before the first byte of the call's reports (in
/// a block search, of the query's block), as the other refusals of such a record
/// (`believed::refuse_shared_masks`). Only a `-subject` search has such a subject, whose hits
/// `joins` never shows in another subject's.
#[allow(clippy::too_many_arguments)]
pub(crate) fn section_aborts(
    view: PairwiseView,
    sort: super::hit_sort::ReportSort,
    num_alignments: usize,
    q_idx: usize,
    rows: &[PairwiseHit],
    subject_titles: &[Arc<[u8]>],
    is_query: impl Fn(u32) -> bool,
    hsps_of: impl FnOnce(&super::x_pairwise_par::QueryGroup<'_>) -> super::query_anchored::Hsps,
) -> io::Result<bool> {
    if view == PairwiseView::Pairwise {
        return Ok(false);
    }
    let Some(group) = blastp_report_group(
        view.report_sort(sort),
        q_idx,
        rows,
        subject_titles,
        |_, _| false,
    ) else {
        return Ok(false);
    };
    if !group
        .aligned(num_alignments)
        .iter()
        .any(|&s_idx| is_query(s_idx))
    {
        return Ok(false);
    }
    super::query_anchored::ncbi_aborts(&hsps_of(&group))
}

/// The group of query `q_idx` of a report (`QueryGroup::sorted`, which `for_sections` also
/// calls; BLASTP's and BLASTN's), `rows` the rows
/// of the call's hits in the final order: its subjects grouped by `joins` (copies of a database
/// record, `QuerySections::joins`) and sorted by `sort` (`PairwiseConfig::sort`). `None` for a
/// query without hits.
fn blastp_report_group<'h>(
    sort: super::hit_sort::ReportSort,
    q_idx: usize,
    rows: &'h [PairwiseHit],
    subject_titles: &[Arc<[u8]>],
    joins: impl Fn(u32, u32) -> bool,
) -> Option<super::x_pairwise_par::QueryGroup<'h>> {
    let query_hits: Vec<(HspIndex, &PairwiseHit)> = rows
        .iter()
        .enumerate()
        .filter(|(_, hit)| hit.hit.q_idx as usize == q_idx)
        .collect();
    super::x_pairwise_par::QueryGroup::sorted(&query_hits, subject_titles, joins, sort, false)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:686-700
/// The subjects (by subject index) whose description lines the report of query `q_idx` of a
/// BLASTP search shows (outfmt 0 to 4, psiblast's rounds), `rows` the rows of the call's hits in
/// the final order, in the order of the lines: `QueryGroup::described_subjects` of the report's
/// group, grouped by `joins` and sorted by `sort` as the writer's (the first `num_descriptions`
/// subjects, or with `found_before`, a psiblast round with the previous round's set,
/// `psiblast_round`, those of the lines of its two-part table). The description line of each
/// looks its title up again by its display id (`common::database::refuse_display_id_of_a_query`);
/// asked before any format of the call is written, only for a query that reports a sequence whose
/// display id a query of its batch has.
pub(crate) fn blastp_described_subjects(
    sort: super::hit_sort::ReportSort,
    found_before: Option<&[bool]>,
    num_descriptions: usize,
    q_idx: usize,
    rows: &[PairwiseHit],
    joins: impl Fn(u32, u32) -> bool,
) -> Vec<u32> {
    blastp_report_group(sort, q_idx, rows, &[], joins).map_or_else(Vec::new, |group| {
        group.described_subjects(num_descriptions, found_before)
    })
}

impl BlastpSections<'_> {
    /// The gapped Karlin block of query `q_idx`'s footer: the report's, with rpsblast's K of
    /// the query (the first query of a batch has the K that the RPS traceback wrote).
    fn gapped_karlin(&self, q_idx: usize) -> KarlinParams {
        query_gapped_karlin(self.report, self.query_gapped_k, q_idx)
    }
}

/// The gapped Karlin block of query `q_idx`'s footer: the report's, with the query's K where
/// the search gives each query its own (`query_gapped_k`: rpsblast's and RPSTBLASTN's).
fn query_gapped_karlin(
    report: &BlastpPairwiseReport,
    query_gapped_k: Option<&[f64]>,
    q_idx: usize,
) -> KarlinParams {
    match query_gapped_k {
        Some(gapped_k) => KarlinParams {
            k: gapped_k[q_idx],
            ..report.gapped_karlin
        },
        None => report.gapped_karlin,
    }
}

impl super::x_pairwise_par::QuerySections for BlastpSections<'_> {
    fn num_alignments(&self) -> usize {
        self.report.num_alignments
    }

    fn line_length(&self) -> usize {
        self.config.line_length
    }

    fn sort(&self) -> (super::hit_sort::ReportSort, bool) {
        (self.config.sort, self.config.tblastx())
    }

    fn joins(&self, previous: u32, s_idx: u32) -> bool {
        // `-subject` records with one id are refused before any report.
        self.database.is_some()
            && crate::report::db_subject::rank_joins_as_copy(
                self.parsed_subjects,
                previous as usize,
                s_idx as usize,
            )
    }

    fn head<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: Option<&super::x_pairwise_par::QueryGroup<'_>>,
        place: super::query_warnings::ReportHook<'_, W>,
    ) -> io::Result<()> {
        let (query, report) = (&self.queries[q_idx], self.report);
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1474-1476
        // psiblast's results carry their round (`PrintOneResultSet` with
        // `itr.GetIterationNumber()`, psiblast_app.cpp:283-285); LOSAT's psiblast runs its
        // first round only.
        if self.config.program == "psiblast" {
            writeln!(writer, "Results from round 1")?;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1491
        // Emit two blank lines before the query header.
        writeln!(writer)?;
        writeln!(writer)?;
        // The title's bytes by query index (`query_name` is BLASTX's).
        let query_title = self
            .query_titles
            .get(q_idx)
            .map_or(&b""[..], |title| title.as_ref());
        write_blastp_query_header(
            writer,
            query_title,
            query.query_length,
            self.config.html.is_some(),
        )?;
        let Some(group) = group else {
            write_no_hits_found(writer)?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180,204-208:
            // The all-invalid batch uses -1 sentinel blocks. An invalid query
            // within a searched batch has null blocks (blast_results.cpp:82-103),
            // so only the blank lines and the zero search space are written
            // (blast_format.cpp:445-477), as TBLASTN does.
            if query.batch_skipped {
                write_tblastn_unsearched_query_footer_spacing(writer, false)?;
            } else if !query.valid {
                // `x_PrintOneQueryFooter`'s first `NcbiEndl` flushes the stream.
                crate::report::ostream::ncbi_endl(writer)?;
                writeln!(writer)?;
                writeln!(writer)?;
                writeln!(writer, "Effective search space used: 0")?;
            } else {
                write_blastp_query_footer_spacing(
                    writer,
                    query.ungapped_karlin,
                    self.gapped_karlin(q_idx),
                    report.gumbel,
                    report.gapped,
                    query.effective_search_space,
                    false,
                )?;
            }
            return Ok(());
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:606-611
        // The table follows `x_InitDeflineTable` (`write_blastn_description_table`): the
        // subject's highest bit score and that HSP's E-value, and the protein title of the
        // record read (by subject index), as TBLASTN's.
        write_description_table(
            writer,
            group.described(report.num_descriptions),
            group.subject_order.len() <= report.num_descriptions,
            &group.subject_hits,
            false,
            self.config.sort_columns(query.query_length),
            |s_idx| {
                let title = group
                    .subject_titles
                    .get(&s_idx)
                    .map_or(&[][..], Vec::as_slice);
                database_description_label(
                    title,
                    stored_source(self.subject_sources, s_idx),
                    parsed_subject(self.parsed_subjects, s_idx),
                    true,
                    self.thrown(s_idx).is_some(),
                )
            },
            place,
            self.config.html_rows().map(|(report, after)| HtmlRows {
                report,
                after,
                parsed: self.parsed_subjects,
                loads: self.loads,
            }),
            None,
        )?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1520-1589
        writeln!(writer)
    }

    fn subject<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        probe: Option<&mut FormatProbe<'_>>,
        group: &super::x_pairwise_par::QueryGroup<'_>,
        s_idx: u32,
        part: std::ops::Range<usize>,
        place: super::query_warnings::ReportHook<'_, W>,
    ) -> io::Result<()> {
        write_blastp_subject_part(
            writer,
            probe,
            self.config,
            self.database.is_some(),
            self.parsed_subjects,
            self.subject_sources,
            self.loads,
            group,
            s_idx,
            part,
            self.row_masks,
            place,
        )
    }

    fn thrown(&self, s_idx: u32) -> Option<(&crate::ncbi_exception::NcbiException, bool)> {
        thrown_load(self.loads, self.parsed_subjects, s_idx)
    }

    fn view(&self) -> PairwiseView {
        self.config.view
    }

    fn html_links(&self) -> Option<&super::html_ids::HtmlLinks> {
        self.config.html_links()
    }

    fn alignment_id_is_local(&self, s_idx: u32) -> bool {
        alignment_id_local(self.parsed_subjects, s_idx)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1550-1552
    /// The query-anchored section of the query's aligned subjects (`blastp_section_hsps`).
    fn alignments<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: &super::x_pairwise_par::QueryGroup<'_>,
    ) -> io::Result<()> {
        write_anchored_section(writer, self.config, self.anchored, q_idx, |sequences| {
            blastp_section_hsps(
                self.report.num_alignments,
                group,
                self.parsed_subjects,
                |s_idx| sequences.subject_is_query(q_idx, s_idx),
            )
        })
    }

    fn foot<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
    ) -> io::Result<()> {
        let query = &self.queries[q_idx];
        write_blastp_query_footer_spacing(
            writer,
            query.ungapped_karlin,
            self.gapped_karlin(q_idx),
            self.report.gumbel,
            self.report.gapped,
            query.effective_search_space,
            false,
        )
    }
}

// =============================================================================
// BLASTN pairwise report
// =============================================================================

/// What the BLASTN prolog's program line and references follow: `-task megablast` (NCBI's
/// `m_Megablast`; the program name stays `blastn`), a search with a MegaBLAST database index
/// (`m_IndexedMegablast`, the `-use_index true` that NCBI's loader kept), and matrix-only scoring
/// (the option reward and penalty both 0, NCBI's `m_DisableKAStats`), which also leaves out the
/// description table, the query footers and the bit scores and e-values of the HSP headers.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:168-177
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct BlastnIntro {
    pub megablast: bool,
    pub indexed_megablast: bool,
    pub matrix_only: bool,
}

/// One BLASTN query of the pairwise report (NCBI `CBlastAncillaryData` and the query's
/// validity).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-104
#[derive(Debug, Clone)]
pub struct BlastnPairwiseQuery {
    /// The query's title bytes (`Query=`; empty without a title).
    pub query_name: Vec<u8>,
    pub query_length: usize,
    /// The ungapped and the gapped block of the query's first valid context; `None` for
    /// an invalid query.
    pub karlin: Option<(KarlinParams, KarlinParams)>,
    pub effective_search_space: i64,
}

/// Run-level data of the BLASTN pairwise report.
#[derive(Debug, Clone)]
pub struct BlastnPairwiseReport {
    pub version: String,
    /// The program line and references of the prolog (`-task megablast`, matrix-only scoring).
    pub intro: BlastnIntro,
    pub database_name: Vec<u8>,
    /// NCBI's `Int8` (`BlastpPairwiseReport::database_num_sequences`).
    pub database_num_sequences: i64,
    pub database_total_letters: i64,
    /// The database of a `-db` search (its `Database:` lines); `None` for a subject set,
    /// whose lines are `database_name` with a period and the counts.
    pub database: Option<PairwiseDatabase>,
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    /// The task is megablast or blastn, for which a gap extension cost of 0 is printed as
    /// the PMID 10890397 value (NCBI compares the task name, `m_Program`).
    pub zero_gap_extension_formula: bool,
    /// The two-hit window (`options.GetWindowSize()`): 40 for dc-megablast.
    pub window_size: i32,
    /// A gapped search (the options' gapped mode); an ungapped (`-ungapped`) search has no
    /// Gapped Karlin block in a query's footer and no `Gap Penalties` line in the epilog.
    pub gapped: bool,
    /// Subjects in the description table and with alignments, per query.
    pub num_descriptions: usize,
    pub num_alignments: usize,
    /// Every query of the batch is invalid, so NCBI did not search it (all queries of such
    /// a batch get the `-1` footer, local_blast.cpp:177-207). Indexed like the queries.
    pub unsearched: Vec<bool>,
    /// Whether NCBI's epilog (the database and statistics footer, `PrintEpilog`) ends the
    /// report.
    pub epilog: bool,
    /// Whether the report starts with NCBI's prolog (`PrintProlog`); false when the caller
    /// wrote it before the search, as NCBI does.
    pub prolog: bool,
    /// The database sequences with parsed Seq-ids by subject index (empty without them): their
    /// descriptions and headings come from their definition lines (`db_subject::ParsedSubject`).
    pub parsed_subjects: Vec<Option<Arc<crate::report::db_subject::ParsedSubject>>>,
    /// The sources of the others' Bioseqs by subject index, from which NCBI makes an empty stored
    /// title (`stored_source`; empty where none has one).
    pub subject_sources: Vec<super::defline::SourceTaxname>,
    /// The OIDs of a database search's listed subjects by subject index (`BL_ORD_ID:<oid>` for
    /// outfmt 0 `...alignment skipped` lines on unparsed database subjects); `None` for `-subject`.
    pub subject_oids: Option<Arc<[u32]>>,
}

/// The HSP of a subject's HSPs `hits` whose values the description table shows
/// (`GetSeqAlignSetCalcParams`): NCBI keeps the highest bit score from 0 and takes the bit score,
/// E-value, identities and length of each HSP that raises it, so the HSP is the first of the
/// highest bit score above 0, and a subject whose HSPs all have a bit score of 0 or less or NaN (a
/// PSSM whose gapped Karlin-Altschul block is not valid gives -1 or NaN) has none. The one
/// predicate of the table's scores (`table_scores`) and of `-sorthits`' identity
/// (`SortColumnValues::of`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4261-4309
pub(super) fn shown_hsp<'a>(hits: &[&'a PairwiseHit]) -> Option<&'a PairwiseHit> {
    let mut shown: Option<&PairwiseHit> = None;
    for &hit in hits {
        if hit.hit.bit_score > shown.map_or(0.0, |shown| shown.hit.bit_score) {
            shown = Some(hit);
        }
    }
    shown
}

/// The bit score and E-value that the description table shows for a subject whose shown HSP is
/// `shown` (`shown_hsp`): 0.0 and 0.0 without one.
pub(super) fn table_scores(shown: Option<&PairwiseHit>) -> (f64, f64) {
    shown.map_or((0.0, 0.0), |hit| (hit.hit.bit_score, hit.hit.e_value))
}

// The description table of a report: BLASTN's, BLASTP's, TBLASTN's and TBLASTX's
// (`write_blastn_description_table`), BLASTX's
// (`write_subject_summary_table_with_sum_n`), which differ only in the line that names a
// subject (`label_of`).
//
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1004-1011
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1100-1142
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4285-4310
// A row shows the highest bit score of the subject's HSPs and that HSP's E-value. The
// widths start from "(Bits)" and "Value" (showdefline.cpp:1059-1062). The last row's
// total score sets the score width when it is wider than every earlier total (the
// width that NCBI assigns there is the score width). The last row exists only when the
// subjects fit in the table: otherwise the loop stops with an empty last group. A total
// is wider than "Total" only above 99999 bits (`%5.3le`), which `format_bitscore_ncbi`
// formats alike. The highest bit score is the first of equal maxima (`>`), and the
// subject's first HSP (the first by E-value) can have a lower one: sum statistics give a
// linked set of lower scores a lower E-value than a single HSP of a higher score.
//
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4229-4230,4259
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1061,1118-1120
// With `show_sum_n` (NCBI's `eShowSumN`: ungapped sum statistics, blast_format.cpp:165-167,
// 503-504) a last column shows the linked-set size of the subject's first HSP; the HSP
// without a "sum_n" score (`num` of 1) shows 1.
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:771-794
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:809-840
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:66-104,1051-1065
// With `-sorthits` (`columns`, NCBI's `eShowTotalScore`, `eShowQueryCoverage` and
// `eShowPercentIdent`) the header and the rows add the total score, the query cover and the
// identity of the subject's best HSP. `m_MaxPercentIdentityLen` starts at the 9 letters of
// "Max ident" and grows only for a percentage of more digits, which none has.
#[allow(clippy::too_many_arguments)]
fn write_description_table<W: Write>(
    writer: &mut W,
    described: &[u32],
    last_row_counted: bool,
    subject_hits: &std::collections::HashMap<u32, Vec<&PairwiseHit>>,
    show_sum_n: bool,
    columns: Option<SortColumns>,
    mut label_of: impl FnMut(u32) -> DescriptionLabel,
    place: super::query_warnings::ReportHook<'_, W>,
    // `-html`: the rows' links (`HtmlRows`); `None`: the plain table.
    html: Option<HtmlRows<'_>>,
    // Where the table's read of its rows' ranges throws (`described_row_range`: BLASTX's HSPs with
    // empty rows, BLASTX's and TBLASTN's rows that start with a gap): that exception; `None` for
    // every other table.
    rows_throw: Option<EmptyRowsStop>,
) -> io::Result<()> {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1068-1070
    // An `-html` table reads `./.ncbirc` first (`x_InitDeflineTable`), also one without rows
    // (`-num_descriptions 0`): its posts come before the table's bytes (line RG T4a2 probes
    // `bp_nd0`, `bx`: after the query's `Length=` line).
    if html.is_some() {
        place(
            super::query_warnings::ReportPlace::Html(super::html_ids::HtmlPost::Ncbirc),
            writer,
        )?;
    }
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1080-1100
    // NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:225-228
    // Then its rows' scores read the ranges of every HSP of each described subject
    // (`x_GetScoreInfoForTable`), where an empty row throws (`Row is empty`), before the table's
    // bytes and after that read (line RG T4b review a B2).
    if let Some(stop) = rows_throw {
        writer.flush()?;
        return Err(io::Error::other(stop));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:513-514
    // c++/src/objtools/align_format/showdefline.cpp:757,848: with `-num_descriptions 0`
    // the header block (the two header lines and the blank line) and the rows are not written;
    // the blank line that closes the table is (a line feed written to `out` at blast_format.cpp:613).
    if described.is_empty() {
        return writeln!(writer);
    }
    let rows: Vec<(
        u32,
        Option<&PairwiseHit>,
        f64,
        i32,
        Option<SortColumnValues>,
    )> = described
        .iter()
        .filter_map(|s_idx| {
            let hits = subject_hits.get(s_idx)?;
            let first = *hits.first()?;
            let best = shown_hsp(hits);
            let sum_n = first.sum_n.filter(|&num| num > 1).unwrap_or(1);
            Some((
                *s_idx,
                best,
                hits.iter()
                    .fold(0.0, |total, hit| total + hit.hit.bit_score),
                sum_n,
                columns.map(|columns| SortColumnValues::of(hits, best, columns)),
            ))
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1059
    let mut max_score = ncbi_k_bits().len();
    let mut max_evalue = "Value".len();
    let mut max_total = "Total".len();
    let mut max_sum_n = 1;
    for (index, (_, best, total, sum_n, _)) in rows.iter().enumerate() {
        max_score = max_score.max(format_bitscore_ncbi(table_scores(*best).0).len());
        max_evalue = max_evalue.max(format_evalue_ncbi(table_scores(*best).1).len());
        max_sum_n = max_sum_n.max(sum_n.to_string().len());
        let total_len = format_total_bitscore_ncbi(*total).len();
        if index + 1 == rows.len() && last_row_counted {
            if max_total < total_len {
                max_score = total_len;
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:940-944
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:932-938
                // The last row's total, wider than the column (above 99999 bits, `%5.3le`),
                // pads with as many spaces as `m_MaxTotalScoreLen` less `size`, a `size_t` below 0: NCBI
                // writes spaces without end. No valid report exists to match, so the run stops
                // with an error (an NCBI failure is an explicit rejection, AGENTS.md rule 3).
                if columns.is_some() {
                    return Err(TotalScoreRefusal::error());
                }
            }
        } else {
            max_total = max_total.max(total_len);
        }
    }
    write_spaces(writer, 70)?;
    write!(writer, "{:<max_score$}  ", "Score")?;
    if columns.is_some() {
        write!(writer, "{:<max_total$}  Query  ", "Total")?;
    }
    write!(writer, "  E")?;
    if columns.is_some() {
        write_spaces(writer, max_evalue - "Value".len() + 5)?;
        write!(writer, "Max")?;
    }
    writeln!(writer)?;
    write!(
        writer,
        "{:<69}",
        "Sequences producing significant alignments:"
    )?;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:811-813
    write!(writer, "{:<max_score$}  ", ncbi_k_bits())?;
    if columns.is_some() {
        write!(writer, " {:<max_total$}  cover   ", "Score")?;
    }
    write!(writer, "Value")?;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:829-836
    if show_sum_n || columns.is_some() {
        write_spaces(writer, max_evalue - "Value".len() + 2)?;
    }
    if show_sum_n {
        write!(writer, "N")?;
    }
    if columns.is_some() {
        write!(writer, "Ident")?;
    }
    writeln!(writer)?;
    writeln!(writer)?;
    // Where the `-html` rows end the run (`HtmlRows::end`).
    let end = html.and_then(|html| html.end(&rows.iter().map(|row| row.0).collect::<Vec<_>>()));
    for (index, (s_idx, best, total, sum_n, values)) in rows.iter().enumerate() {
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1587-1589
        // Each line loads its subject (`x_GetDeflineInfo` of `x_DisplayDeflineTableBody`) after
        // the lines before it are written (`QueryWarnings::subject_read`).
        place(super::query_warnings::ReportPlace::Load(*s_idx), writer)?;
        super::x_pairwise_par::end_here(end, index)?;
        let bits = format_bitscore_ncbi(table_scores(*best).0);
        write_description_row_start(
            writer,
            &label_of(*s_idx),
            (*s_idx, index + 1),
            &bits,
            max_score,
            html,
            &mut *place,
        )?;
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:940-956
        if let Some(values) = values {
            write!(
                writer,
                "   {:<max_total$}",
                format_total_bitscore_ncbi(*total)
            )?;
            let coverage = format!("{}%", values.percent_coverage);
            write!(writer, "  {coverage:<5}")?;
        }
        write!(
            writer,
            "  {:<max_evalue$}",
            format_evalue_ncbi(table_scores(*best).1)
        )?;
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:961-965
        if show_sum_n {
            write!(writer, "  {sum_n:<max_sum_n$}")?;
        }
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:966-974
        if let Some(values) = values {
            let identity = values.percent_identity.to_string();
            write!(writer, "  {identity}%")?;
            write_spaces(writer, "Max ident".len().saturating_sub(identity.len()))?;
        }
        writeln!(writer)?;
    }
    writeln!(writer)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:940-944
/// LOSAT's refusal of a `-sorthits` description table whose last row's total score is wider than
/// its column (`write_description_table`), where NCBI BLAST+ writes spaces without end: no
/// error of NCBI's, so no program reports it as NCBI's "BLAST failed to write output". It
/// travels as the `io::Error` of the report's writes (`TotalScoreRefusal::error`), and a
/// program that turns its report's write errors into NCBI's passes it on unchanged
/// (`is_total_score_refusal`, BLASTX's `report_error`), so every program stops with this text
/// and exit status 1.
#[derive(Debug)]
pub(crate) struct TotalScoreRefusal;

impl TotalScoreRefusal {
    /// The refusal as the error of the report's writes.
    pub(crate) fn error() -> io::Error {
        io::Error::new(io::ErrorKind::Unsupported, TotalScoreRefusal)
    }
}

impl std::fmt::Display for TotalScoreRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(
            "-sorthits: the total score of the description table's last subject is wider than its column (more than 99999 bits), where NCBI BLAST+ 2.17.0 writes spaces without end (CAlignFormatUtil::AddSpace of a negative size); run without -sorthits",
        )
    }
}

impl std::error::Error for TotalScoreRefusal {}

/// Whether `error` is the report's `TotalScoreRefusal` (an `io::Error` that holds it).
pub(crate) fn is_total_score_refusal(error: &anyhow::Error) -> bool {
    error
        .downcast_ref::<io::Error>()
        .and_then(io::Error::get_ref)
        .is_some_and(|inner| inner.is::<TotalScoreRefusal>())
}

/// The values of one row's columns of `-sorthits` (`SortColumns`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct SortColumnValues {
    /// `percent_coverage`.
    percent_coverage: i64,
    /// The shown identity: `percent_identity` plus 0.5, cast to `int`, 99 above 100.
    percent_identity: i64,
}

impl SortColumnValues {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4271-4306
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1525-1530
    /// The columns of the subject whose HSPs are `hits` and whose shown HSP is `best`
    /// (`shown_hsp`; without one NCBI's identities 0 over a length of 1): `GetSeqAlignSetCalcParams`
    /// with `m_TranslatedNucAlignment` false (TBLASTX's identities over its nucleotide lengths). The
    /// cover is computed without NCBI's `int` product `100*master_covered_length`, which
    /// overflows above about 21 million covered letters.
    fn of(hits: &[&PairwiseHit], best: Option<&PairwiseHit>, columns: SortColumns) -> Self {
        use super::hit_sort::{percent_identity, query_coverage, HspScores};
        let scores: Vec<HspScores> = hits.iter().map(|hit| HspScores::of(hit)).collect();
        let covered = query_coverage(&scores);
        let mut percent_coverage =
            (0.5 + (100 * covered) as f64 / columns.query_length as f64) as i64;
        if percent_coverage > 100 {
            percent_coverage = percent_coverage.min(99);
        }
        let (num_ident, length) = best.map_or((0, 1), |best| {
            let best = HspScores::of(best);
            (
                best.num_ident,
                best.alignment_length(columns.tblastx, false),
            )
        });
        let identity = percent_identity(num_ident, length);
        let mut percent_identity = (0.5 + identity) as i64;
        if percent_identity > 100 {
            percent_identity = percent_identity.min(99);
        }
        Self {
            percent_coverage,
            percent_identity,
        }
    }
}

/// The subject of `s_idx` among `parsed` when it has parsed Seq-ids (not a `-subject` record shown
/// as one without `-parse_deflines`, `ParsedSubject::shown_as_record`).
/// Whether NCBI's scope gets no Bioseq for subject `s_idx` of a report whose subjects with parsed
/// Seq-ids are `parsed` (by subject index; `ParsedSubject::scope_finds_none`).
fn finds_no_bioseq(
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    s_idx: u32,
) -> bool {
    parsed_subject(parsed, s_idx)
        .is_some_and(crate::report::db_subject::ParsedSubject::scope_finds_none)
}

fn parsed_subject(
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    s_idx: u32,
) -> Option<&crate::report::db_subject::ParsedSubject> {
    parsed
        .get(s_idx as usize)
        .and_then(Option::as_deref)
        .filter(|subject| !subject.shown_as_record())
}

/// The id of NCBI's line in place of each alignment of the subject `s_idx` whose alignment NCBI
/// cannot write (`write_subject_alignments_part`, `UnreadableQuery`): `GetSeqIdString` of the
/// Seq-align id. A subject with parsed ids has its parsed id (`ParsedSubject::skipped_id`,
/// `Subject_9` of `lcl|Subject_9`); a database sequence without them (`oids`, the OIDs of the
/// listed subjects) the general id `BL_ORD_ID:<oid>`; a `-subject` record the local id
/// `Subject_N` (N = `s_idx` + 1).
fn subject_skipped_id(
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    oids: Option<&[u32]>,
    s_idx: u32,
) -> Vec<u8> {
    match (parsed.get(s_idx as usize).and_then(Option::as_deref), oids) {
        (Some(subject), _) => subject.skipped_id().to_vec(),
        (None, Some(oids)) => format!("BL_ORD_ID:{}", oids[s_idx as usize]).into_bytes(),
        (None, None) => format!("Subject_{}", u64::from(s_idx) + 1).into_bytes(),
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2145-2154
/// The Seq-id text of the Seq-align id of the `-subject` record `s_idx` without and with its
/// version (`CSeq_id::GetSeqIdString` with false and true, `x_pairwise_par::QuerySections::
/// alignment_id_texts`): `subject_skipped_id`, and with the version the label of a parsed id
/// (`seqid::label::content`); the local id `Subject_N` of a record without one has no version.
fn record_alignment_id_texts(
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    s_idx: u32,
) -> [Vec<u8>; 2] {
    let id = subject_skipped_id(parsed, None, s_idx);
    let versioned = match parsed.get(s_idx as usize).and_then(Option::as_deref) {
        Some(subject) => {
            let mut text = Vec::new();
            crate::seqid::label::content(subject.alignment_id(), true, &mut text);
            text
        }
        None => id.clone(),
    };
    [id, versioned]
}

/// A description line's label in NCBI's two parts (`x_DisplayDeflineTableBody`): the id text
/// (`gi|N|` and the id, the `line_length` that NCBI counts before the title; under `-html` the
/// text of the sequence's link) and the title (`line_component`: the defline, after two spaces
/// where a gi was written; under `-html` HTML-encoded). The plain line writes them as one string
/// (`write_cut`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:883-913
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DescriptionLabel {
    /// The id text (`ParsedSubject::description`'s start; empty for a title alone).
    pub id: Vec<u8>,
    /// The title.
    pub title: Vec<u8>,
    /// Whether the line is that of NCBI's `catch` (`x_GetDeflineInfo`: the subject's title or its
    /// load throws; the id text is the Seq-align id's and the title `Unknown`): its `-html` row
    /// has no score link and links the Seq-align id (`HtmlReport::description_links`).
    pub caught: bool,
}

impl DescriptionLabel {
    /// A label of a title alone (a generated defline).
    pub fn title_only(title: Vec<u8>) -> Self {
        DescriptionLabel {
            id: Vec::new(),
            title,
            caught: false,
        }
    }

    /// The plain label on the description line: the id text and the title as one string cut to
    /// 68 columns (a longer one is its first 65 bytes and `...`), or padded to them.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:915-918,930
    /// String widths are byte counts, including non-ASCII FASTA titles. (NCBI cuts the title
    /// alone, which is the same cut for an id text of up to 65 bytes; a longer id is refused
    /// before any output (`common/database.rs`, `report/db_subject.rs`), so the cut of the
    /// whole label below never cuts the id.)
    fn write_cut<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        // A label of a title alone is written as the joined label was: one write of its letters.
        let (id, title, cut) = self.cut();
        if id > 0 {
            writer.write_all(&self.id[..id])?;
        }
        writer.write_all(&self.title[..title])?;
        if cut {
            writer.write_all(b"...")
        } else {
            write_spaces(writer, 68 - id - title)
        }
    }

    /// The bytes of the id text and of the title that the line keeps and whether `...` follows
    /// (`write_cut`): a label longer than 68 columns keeps 65 bytes.
    fn cut(&self) -> (usize, usize, bool) {
        if self.id.len() + self.title.len() > 68 {
            let id = self.id.len().min(65);
            (id, 65 - id, true)
        } else {
            (self.id.len(), self.title.len(), false)
        }
    }

    /// The `-html` label (`write_cut` with `eHtml`): the link `link` (the line's id URL), the id
    /// text, `</a>` after a link, then the cut title `HTMLEncode`d; the cut and the padding count
    /// the raw bytes (`line_length`).
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:883-885, 909-931
    fn write_html<W: Write>(&self, writer: &mut W, link: &[u8]) -> io::Result<()> {
        let (id, title, cut) = self.cut();
        writer.write_all(link)?;
        writer.write_all(&self.id[..id])?;
        if !link.is_empty() {
            writer.write_all(b"</a>")?;
        }
        let mut encoded = Vec::new();
        super::html::html_encode_into(&self.title[..title], &mut encoded);
        writer.write_all(&encoded)?;
        if cut {
            writer.write_all(b"...")
        } else {
            write_spaces(writer, 68 - id - title)
        }
    }
}

/// The label and bit-score columns of one description row (`write_description_table` and
/// `psiblast_round::RoundSections::write_part`): plain (`DescriptionLabel::write_cut` and the
/// score), or `-html` (`HtmlRows::write_start`).
fn write_description_row_start<W: Write>(
    writer: &mut W,
    label: &DescriptionLabel,
    rank: (u32, usize),
    bits: &str,
    max_score: usize,
    html: Option<HtmlRows<'_>>,
    place: super::query_warnings::ReportHook<'_, W>,
) -> io::Result<()> {
    match html {
        None => {
            label.write_cut(writer)?;
            write!(writer, "  {bits:<max_score$}")
        }
        Some(html) => html.write_start(writer, label, rank, bits, max_score, place),
    }
}

/// The `-html` inputs of a description table (`write_description_table`): the report's links and
/// OIDs (`PairwiseConfig::html`) and the subjects with parsed Seq-ids by subject index (those of
/// the labels).
#[derive(Clone, Copy)]
struct HtmlRows<'a> {
    report: &'a super::html_ids::HtmlReport,
    /// What the report does after a throw (`PairwiseConfig::after_throw`), where its description
    /// lines end the run (`AfterThrow::description_end`).
    after: crate::report::db_subject::AfterThrow,
    parsed: &'a [Option<Arc<crate::report::db_subject::ParsedSubject>>],
    /// What the report's loads do at each subject, by subject (`HitLoads::thrown`).
    loads: &'a crate::report::db_subject::HitLoads,
}

impl<'a> HtmlRows<'a> {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1588-1611
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3904-3907
    /// Where the `-html` rows of the subjects `described` (in their order) end the run: the
    /// position of the first whose load throws (`HitLoads::thrown`) and its exception, as the
    /// report's answer says (`AfterThrow::description_end`, which the model of the loads asks
    /// too): the `catch` of the row's load makes the row's link, whose `GetIDUrl` loads the
    /// subject again with the scope, and that throw leaves the report. The writer returns the
    /// exception there as its error (`x_pairwise_par::end_here`), after the rows before it and the
    /// row's own read; the plain row is the `catch`'s `Unknown`.
    fn end(self, described: &[u32]) -> Option<(usize, &'a crate::ncbi_exception::NcbiException)> {
        let at = self.after.description_end(described.iter(), |&s_idx| {
            self.loads.thrown(s_idx as usize).is_some()
        })?;
        Some((at, self.loads.thrown(described[at] as usize)?))
    }

    /// The start of the `-html` description row of rank `rank` (from 1) for the subject `s_idx`
    /// with the label `label` and the bit score `bits` (in a column of `max_score`): the label
    /// (`DescriptionLabel::write_html`), two spaces, then the bit score inside the score link (none
    /// on the row of NCBI's `catch`, `DescriptionLabel::caught`). The row's links read NCBI's
    /// formatter registry (`HtmlPost`), which `place` writes before the row's bytes, also before
    /// the registry's exception.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:932-939
    fn write_start<W: Write>(
        self,
        writer: &mut W,
        label: &DescriptionLabel,
        (s_idx, rank): (u32, usize),
        bits: &str,
        max_score: usize,
        place: super::query_warnings::ReportHook<'_, W>,
    ) -> io::Result<()> {
        let parsed = self.parsed.get(s_idx as usize).and_then(Option::as_deref);
        let mut posts = Vec::new();
        let links = self
            .report
            .description_links(parsed, s_idx, rank, label.caught, &mut posts);
        for post in posts {
            place(super::query_warnings::ReportPlace::Html(post), writer)?;
        }
        let (link, score) = links?;
        label.write_html(writer, &link)?;
        writer.write_all(b"  ")?;
        writer.write_all(&score)?;
        writer.write_all(bits.as_bytes())?;
        if !score.is_empty() {
            writer.write_all(b"</a>")?;
        }
        write_spaces(writer, max_score.saturating_sub(bits.len()))
    }
}

/// The description-table label of a database sequence or a `-subject` record (before the cut to
/// the line): the title that `GenerateDefline` makes of the stored title `title` and the source
/// `source` of the Bioseq, which an empty title is made from, or for a sequence with parsed
/// Seq-ids its id and title (`ParsedSubject::description_line`); where the title throws
/// (`Unknown`: NCBI's `HtmlDecode` cannot decode it) or the load does (`thrown`), the line of
/// NCBI's `catch` (`DescriptionLabel::caught`; the Seq-align id of a sequence with parsed Seq-ids,
/// `ParsedSubject::thrown_description`; none is shown of one without).
pub fn database_description_label(
    title: &[u8],
    source: &super::defline::SourceTaxname,
    parsed: Option<&crate::report::db_subject::ParsedSubject>,
    protein: bool,
    thrown: bool,
) -> DescriptionLabel {
    if let Some(subject) = parsed {
        let (id, title, caught) = if thrown {
            let (id, title) = subject.thrown_description();
            (id, title, true)
        } else {
            subject.description_line(protein)
        };
        return DescriptionLabel {
            id,
            title: title.text,
            caught,
        };
    }
    let facts = super::defline::BioseqFacts::of_source(source.name().ok().flatten());
    match super::defline::generate_bioseq_defline(title, protein, true, &facts) {
        Ok(defline) => DescriptionLabel::title_only(defline.text),
        Err(super::defline::UnknownEncoding) => DescriptionLabel {
            id: Vec::new(),
            title: b"Unknown".to_vec(),
            caught: true,
        },
    }
}

// The description table of BLASTN, BLASTP, TBLASTN and TBLASTX: a subject's line is the
// title of its record (`write_description_table`).
fn write_blastn_description_table<W: Write>(
    writer: &mut W,
    described: &[u32],
    last_row_counted: bool,
    subject_hits: &std::collections::HashMap<u32, Vec<&PairwiseHit>>,
    subject_titles: &std::collections::HashMap<u32, Vec<u8>>,
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    sources: &[super::defline::SourceTaxname],
    show_sum_n: bool,
    protein: bool,
    columns: Option<SortColumns>,
    place: super::query_warnings::ReportHook<'_, W>,
    // `-html`: the report's links and OIDs (`PairwiseConfig::html`), with what it does after a
    // throw (`PairwiseConfig::html_rows`).
    html: Option<(
        &super::html_ids::HtmlReport,
        crate::report::db_subject::AfterThrow,
    )>,
    loads: &crate::report::db_subject::HitLoads,
    // Where the rows' ranges throw (`write_description_table`; TBLASTN's, `described_row_range`).
    rows_throw: Option<EmptyRowsStop>,
) -> io::Result<()> {
    write_description_table(
        writer,
        described,
        last_row_counted,
        subject_hits,
        show_sum_n,
        columns,
        |s_idx| {
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:498
            // The description keeps the prefixes (`fLeavePrefixSuffix`, `defline.rs`).
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1579-1592
            // A subject whose load throws where the report reads it (`loads`,
            // `thrown_description`) takes the `catch` of that line too: its description is the
            // Seq-align id and `Unknown`.
            let title = subject_titles.get(&s_idx).map_or(&[][..], Vec::as_slice);
            database_description_label(
                title,
                stored_source(sources, s_idx),
                parsed_subject(parsed, s_idx),
                protein,
                loads.thrown(s_idx as usize).is_some(),
            )
        },
        place,
        html.map(|(report, after)| HtmlRows {
            report,
            after,
            parsed,
            loads,
        }),
        rows_throw,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:451-478
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:472-480
// Without a Gumbel block, PrintKAParameters (align_format_util.cpp:585-619) prints only
// the Lambda, K and H columns. An invalid query has neither Karlin block (`karlin` is
// `None`) and a search space of 0.
/// The footer of a nucleotide query, with the two blank lines that follow it (the next
/// query's preamble or the epilog) when `trailing`.
fn write_nucleotide_query_footer_spacing<W: crate::report::ostream::ReportWrite>(
    writer: &mut W,
    karlin: Option<(KarlinParams, KarlinParams)>,
    gapped: bool,
    effective_search_space: i64,
    trailing: bool,
) -> io::Result<()> {
    // `x_PrintOneQueryFooter`'s first `NcbiEndl` flushes the stream.
    crate::report::ostream::ncbi_endl(writer)?;
    // (A block that the results do not keep is not written, `stats::tables::ancillary_kept`.)
    if let Some((ungapped, _)) =
        karlin.filter(|(ungapped, _)| crate::stats::tables::ancillary_kept(ungapped.lambda))
    {
        writeln!(writer, "Lambda      K        H")?;
        write_ncbi_ka_field(writer, ungapped.lambda)?;
        write_ncbi_ka_field(writer, ungapped.k)?;
        write_ncbi_ka_field(writer, ungapped.h)?;
        writeln!(writer)?;
    }
    writeln!(writer)?;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1130-1159
    // An ungapped search has no gapped Karlin block (`kbp_gap` is not made).
    if let Some((_, gapped_block)) = karlin
        .filter(|_| gapped)
        .filter(|(_, gapped_block)| crate::stats::tables::ancillary_kept(gapped_block.lambda))
    {
        writeln!(writer, "Gapped")?;
        writeln!(writer, "Lambda      K        H")?;
        write_ncbi_ka_field(writer, gapped_block.lambda)?;
        write_ncbi_ka_field(writer, gapped_block.k)?;
        write_ncbi_ka_field(writer, gapped_block.h)?;
        writeln!(writer)?;
    }
    writeln!(writer)?;
    writeln!(
        writer,
        "Effective search space used: {}",
        effective_search_space
    )?;
    // The two blank lines of the next query's preamble or of the epilog
    // (blast_format.cpp:1491, 2249).
    if trailing {
        writeln!(writer)?;
        writeln!(writer)?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2261-2279
// The word threshold and the window size are 0 for blastn (blast_format.cpp:2281-2288),
// so those lines are not printed. The extension is a whole or half number, which the
// C++ stream and Rust print alike (`2`, `2.5`).
fn write_blastn_final_footer<W: ReportWrite>(
    writer: &mut W,
    report: &BlastnPairwiseReport,
    html: bool,
) -> io::Result<()> {
    match report.database.as_ref() {
        Some(database) => write_final_database_block(writer, database)?,
        None => write_final_database_report(
            writer,
            &report.database_name,
            report.database_num_sequences,
            report.database_total_letters,
        )?,
    }
    writeln!(
        writer,
        "Matrix: blastn matrix {} {}",
        report.reward, report.penalty
    )?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2270-2288
    let gap_extension = if report.zero_gap_extension_formula && report.gap_extend == 0 {
        f64::from(-2 * report.penalty + report.reward) / 2.0
    } else {
        f64::from(report.gap_extend)
    };
    if report.gapped {
        writeln!(
            writer,
            "Gap Penalties: Existence: {}, Extension: {}",
            report.gap_open, gap_extension
        )?;
    }
    write_epilog_end(writer, 0.0, report.window_size as i64, html)
}

/// The start of the BLASTN report, which NCBI writes before it searches: the program
/// and its reference, then the subjects (the subject set's name and counts, or the
/// `Database:` lines of a database, `write_blast_database_header`).
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:258-261
pub fn write_blastn_pairwise_prolog<W: ReportWrite>(
    writer: &mut W,
    version: &str,
    intro: BlastnIntro,
    database_name: &[u8],
    database_num_sequences: i64,
    database_total_letters: i64,
    database: Option<&PairwiseDatabase>,
    html: Option<(&super::html_ids::HtmlLinks, &mut dyn Write)>,
) -> io::Result<()> {
    if intro.matrix_only {
        write_prolog(
            writer,
            "blastn",
            version,
            intro.megablast,
            intro.indexed_megablast,
            true,
            html,
        )?;
    } else if intro.megablast {
        write_prolog(
            writer,
            "blastn",
            version,
            true,
            intro.indexed_megablast,
            false,
            html,
        )?;
    } else {
        write_prolog(writer, "blastn", version, false, false, false, html)?;
    }
    if let Some(database) = database {
        return write_blast_database_header(writer, database);
    }
    write_blastp_database_header_spacing(
        writer,
        database_name,
        database_num_sequences,
        database_total_letters,
        1,
    )
}

/// Writes the BLASTN pairwise report (outfmt 0) of the queries `query_range` of `queries`.
///
/// `hits` is the final hit list of those queries in its order (queries in input order, then
/// subjects and HSPs in result order); each hit's index in it identifies it to the probe. The
/// prolog comes with the first query of `queries`, the epilog with the last: the reports of
/// consecutive ranges written in turn are the report of all of them.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1491-1589
// Kept out of line: this runs only for outfmt 0 or hit records, and inlining it into the
// shared post-processing makes every BLASTN run compile it in Wasm hosts.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub fn write_blastn_pairwise_report<W: crate::report::ostream::ReportWrite>(
    hits: &[PairwiseHit],
    writer: &mut W,
    config: &PairwiseConfig,
    queries: &[BlastnPairwiseQuery],
    query_range: std::ops::Range<usize>,
    subject_titles: &[Arc<[u8]>],
    report: &BlastnPairwiseReport,
    probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut super::query_warnings::QueryWarnings<'_>>,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
    anchored: Option<&dyn super::query_anchored::ReportSequences>,
    // The rows' masks of `hits` under `-html` (`HtmlRowMasks`; `None` without).
    row_masks: Option<&[super::query_anchored::HspMasks]>,
    loads: &crate::report::db_subject::HitLoads,
) -> io::Result<()> {
    let mut buffered = io::BufWriter::new(writer);
    let writer = &mut buffered;

    // The prolog with the report's first query, the epilog with its last (`query_range` is all of
    // `queries`, or one window of a caller that makes the rows a window at a time,
    // `x_pairwise_par::query_windows`).
    if report.prolog && query_range.start == 0 {
        write_blastn_pairwise_prolog(
            writer,
            &report.version,
            report.intro,
            &report.database_name,
            report.database_num_sequences,
            report.database_total_letters,
            report.database.as_ref(),
            html_with_diagnostics(config.html_links(), warnings.as_deref_mut())?,
        )?;
    }
    let mut hits_by_query: Vec<Vec<(HspIndex, &PairwiseHit)>> = vec![Vec::new(); query_range.len()];
    for (hsp_index, hit) in hits.iter().enumerate() {
        if let Some(bucket) = (hit.hit.q_idx as usize)
            .checked_sub(query_range.start)
            .and_then(|offset| hits_by_query.get_mut(offset))
        {
            bucket.push((hsp_index, hit));
        }
    }
    // The queries' reports, on `x_pool` with LOSAT_X_PAIRWISEPAR (`x_pairwise_par`).
    super::x_pairwise_par::write(
        writer,
        &BlastnSections {
            config,
            queries,
            report,
            anchored,
            row_masks: HtmlRowMasks::new(row_masks, 0),
            loads,
        },
        query_range.start,
        &hits_by_query,
        subject_titles,
        probe,
        warnings,
        x_pool,
    )?;

    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:317-318
    // An error in a later query batch skips the epilog, and nothing follows the last
    // query's footer. The epilog starts with two blank lines (blast_format.cpp:2249).
    if report.epilog && query_range.end == queries.len() {
        crate::report::ostream::ncbi_endl(writer)?;
        crate::report::ostream::ncbi_endl(writer)?;
        write_blastn_final_footer(writer, report, config.html.is_some())?;
    }
    writer.flush()
}

/// The parts of one query's BLASTN pairwise report (`x_pairwise_par::QuerySections`).
struct BlastnSections<'a> {
    config: &'a PairwiseConfig,
    queries: &'a [BlastnPairwiseQuery],
    report: &'a BlastnPairwiseReport,
    /// The sequences of the query-anchored sections (outfmt 1 to 4; `None` for outfmt 0).
    anchored: Option<&'a dyn super::query_anchored::ReportSequences>,
    /// The rows' masks of the HSPs under `-html` (outfmt 0).
    row_masks: HtmlRowMasks<'a>,
    /// What the reports' loads do at each subject, by subject (`thrown_load`).
    loads: &'a crate::report::db_subject::HitLoads,
}

impl super::x_pairwise_par::QuerySections for BlastnSections<'_> {
    fn num_alignments(&self) -> usize {
        self.report.num_alignments
    }

    fn line_length(&self) -> usize {
        self.config.line_length
    }

    fn sort(&self) -> (super::hit_sort::ReportSort, bool) {
        (self.config.sort, self.config.tblastx())
    }

    fn joins(&self, previous: u32, s_idx: u32) -> bool {
        // `-subject` records with one id are refused before any report.
        self.report.database.is_some()
            && crate::report::db_subject::rank_joins_as_copy(
                &self.report.parsed_subjects,
                previous as usize,
                s_idx as usize,
            )
    }

    fn head<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: Option<&super::x_pairwise_par::QueryGroup<'_>>,
        place: super::query_warnings::ReportHook<'_, W>,
    ) -> io::Result<()> {
        let (query, report) = (&self.queries[q_idx], self.report);
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1489-1491
        writer.write_all(b"\n\n")?;
        write_blastp_query_header(
            writer,
            &query.query_name,
            query.query_length,
            self.config.html.is_some(),
        )?;
        let Some(group) = group else {
            write_no_hits_found(writer)?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:445-449
            if report.intro.matrix_only {
                return Ok(());
            }
            if report.unsearched.get(q_idx).copied().unwrap_or(false) {
                write_tblastn_unsearched_query_footer_spacing(writer, false)?;
            } else {
                write_nucleotide_query_footer_spacing(
                    writer,
                    query.karlin,
                    report.gapped,
                    query.effective_search_space,
                    false,
                )?;
            }
            return Ok(());
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1523-1530
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3074-3096
        // NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:235-237
        // For gapped alignments (`Dense-seg`), `-sorthits 4` (`SortHitByMasterCoverageDescending`
        // over >= 2 hits) and `-sorthsps 2` (`SortHspByMasterStartAscending` over a hit with >= 2
        // HSPs) read row 0 via `GetSeqRange` or `GetSeqStop` before the description table and
        // throw `Row is empty` on a class N HSP, even with `-num_descriptions 0` or matrix-only.
        if report.gapped {
            let sort = self.config.sort;
            let hits_throw = sort.hits == Some(4)
                && group.subject_order.len() >= 2
                && group.subject_order.iter().any(|s| {
                    group.subject_hits[s]
                        .iter()
                        .any(|h| h.query_left_of_record())
                });
            let hsps_throw = sort.hsps == Some(2)
                && group.subject_order.iter().any(|s| {
                    let hsps = &group.subject_hits[s];
                    hsps.len() >= 2 && hsps.iter().any(|h| h.query_left_of_record())
                });
            if hits_throw || hsps_throw {
                writer.flush()?;
                return Err(
                    crate::algorithm::blastn::left_of_query::dense_seg_empty_row_stop()
                        .into_io_error(),
                );
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:606-611
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1533-1543
        if report.intro.matrix_only {
            return writeln!(writer);
        }
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1068-1155
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:2914,4271
        // `CShowBlastDefline::x_InitDeflineTable` reads row 0 via `GetSeqRange` for every HSP of
        // the first `num_descriptions` subjects before writing the table; a gapped (`Dense-seg`)
        // class N HSP throws `CDense_seg::GetSeqStop(): Row is empty` after `-html`'s `.ncbirc`
        // lookup.
        let described = group.described(report.num_descriptions);
        if report.gapped
            && described.iter().any(|s| {
                group.subject_hits[s]
                    .iter()
                    .any(|h| h.query_left_of_record())
            })
        {
            if self.config.html.is_some() {
                place(
                    super::query_warnings::ReportPlace::Html(super::html_ids::HtmlPost::Ncbirc),
                    writer,
                )?;
            }
            writer.flush()?;
            return Err(
                crate::algorithm::blastn::left_of_query::dense_seg_empty_row_stop().into_io_error(),
            );
        }
        // The titles of the subjects of the records read (by subject index).
        write_blastn_description_table(
            writer,
            described,
            group.subject_order.len() <= report.num_descriptions,
            &group.subject_hits,
            &group.subject_titles,
            &report.parsed_subjects,
            &report.subject_sources,
            false,
            false,
            self.config.sort_columns(query.query_length),
            place,
            self.config.html_rows(),
            self.loads,
            None,
        )?;
        writeln!(writer)
    }

    fn subject<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        probe: Option<&mut FormatProbe<'_>>,
        group: &super::x_pairwise_par::QueryGroup<'_>,
        s_idx: u32,
        part: std::ops::Range<usize>,
        place: super::query_warnings::ReportHook<'_, W>,
    ) -> io::Result<()> {
        let config = self.config;
        // The HSP header: the raw score alone under matrix-only scoring
        // (`write_hsp_info_raw_score`), chosen once for the subject.
        let hsp_info = if self.report.intro.matrix_only {
            write_hsp_info_raw_score::<W>
        } else {
            write_hsp_info::<W>
        };
        // A database sequence's heading is its stored title after `>` (a `BL_ORD_ID`
        // id is not shown), as BLASTP's (`write_database_subject_alignments_part`), or the lines
        // of the definition lines of one with parsed Seq-ids.
        let report = self.report;
        let (parsed, sources) = (&report.parsed_subjects[..], &report.subject_sources[..]);
        let (database, gapped) = (report.database.is_some(), report.gapped);
        write_subject_part(
            writer,
            probe,
            group,
            s_idx,
            part,
            false,
            config.html.as_ref(),
            parsed,
            report.subject_oids.as_deref(),
            |subject| {
                group_heading(
                    group,
                    parsed,
                    sources,
                    database,
                    gapped,
                    self.thrown(subject).is_some(),
                    subject,
                )
            },
            place,
            |writer, hit, hsp_index| {
                hsp_info(writer, hit, config)?;
                write_alignment(writer, hit, self.row_masks.of(hsp_index), config)
            },
        )
    }

    fn thrown(&self, s_idx: u32) -> Option<(&crate::ncbi_exception::NcbiException, bool)> {
        thrown_load(self.loads, &self.report.parsed_subjects, s_idx)
    }

    fn view(&self) -> PairwiseView {
        self.config.view
    }

    fn html_links(&self) -> Option<&super::html_ids::HtmlLinks> {
        self.config.html_links()
    }

    fn alignment_id_is_local(&self, s_idx: u32) -> bool {
        alignment_id_local(&self.report.parsed_subjects, s_idx)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1550-1552
    /// The query-anchored section of the query's aligned subjects (`blastn_section_hsps`).
    fn alignments<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: &super::x_pairwise_par::QueryGroup<'_>,
    ) -> io::Result<()> {
        write_anchored_section(writer, self.config, self.anchored, q_idx, |sequences| {
            blastn_section_hsps(
                self.report.num_alignments,
                group,
                &self.report.parsed_subjects,
                |s_idx| sequences.subject_is_query(q_idx, s_idx),
            )
        })
    }

    fn foot<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
    ) -> io::Result<()> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:445-449
        if self.report.intro.matrix_only {
            return Ok(());
        }
        let query = &self.queries[q_idx];
        write_nucleotide_query_footer_spacing(
            writer,
            query.karlin,
            self.report.gapped,
            query.effective_search_space,
            false,
        )
    }
}

/// The query-anchored section (`query_anchored::write_section`) of query `q_idx` in the view of
/// `config` (outfmt 1 to 4; with `-html`, `PairwiseConfig::html`, the rows' mask spans), of the
/// HSPs that `hsps_of` makes with the report's sequences `anchored` (`QuerySections::alignments`
/// of each program).
fn write_anchored_section<W: Write>(
    writer: &mut W,
    config: &PairwiseConfig,
    anchored: Option<&dyn super::query_anchored::ReportSequences>,
    q_idx: usize,
    hsps_of: impl FnOnce(&dyn super::query_anchored::ReportSequences) -> super::query_anchored::Hsps,
) -> io::Result<()> {
    let (PairwiseView::QueryAnchored { identity, anchored }, Some(sequences)) =
        (config.view, anchored)
    else {
        return Err(io::Error::other(
            "LOSAT internal error: a query-anchored report without its sequences",
        ));
    };
    // A query whose sequence the display cannot read (`UnreadableQuery`): NCBI merges its HSPs
    // as a nucleotide query's (its Bioseq is not a protein), then each merged frame's display
    // throws and is caught; where the merge throws `CAlnException`, LOSAT reproduces it, and where
    // it hits undefined behavior, crashes, or wraps a segment length, LOSAT refuses
    // (`query_anchored::write_unreadable_section`).
    if config.unreadable_query.is_some() {
        let mut hsps = hsps_of(sequences);
        hsps.set_query_not_protein();
        return super::query_anchored::write_unreadable_section(
            writer,
            &hsps,
            &config.program.to_ascii_uppercase(),
        );
    }
    super::query_anchored::write_section(
        writer,
        super::query_anchored::View {
            identity,
            anchored,
            html: config.html.is_some(),
        },
        config.line_length,
        &hsps_of(sequences),
        &super::query_anchored::OfQuery { sequences, q_idx },
    )
}

// =============================================================================
// TBLASTX pairwise report (outfmt 0)
// =============================================================================

/// Per-query data of the TBLASTX pairwise report.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-100
#[derive(Debug, Clone)]
pub struct TblastxPairwiseQuery {
    /// The query's title bytes (`Query=`; empty without a title).
    pub query_name: Vec<u8>,
    pub query_length: usize,
    /// The ungapped block of the query's first valid context; `None` for an invalid query.
    pub karlin: Option<KarlinParams>,
    pub effective_search_space: i64,
}

/// Run-level data of the TBLASTX pairwise report.
#[derive(Debug, Clone)]
pub struct TblastxPairwiseReport {
    pub version: String,
    pub database_name: Vec<u8>,
    /// NCBI's `Int8` (`BlastpPairwiseReport::database_num_sequences`).
    pub database_num_sequences: i64,
    pub database_total_letters: i64,
    /// The `-matrix` name as typed (`options.GetMatrixName()`), BLOSUM62 when omitted.
    pub matrix_name: String,
    /// `-threshold` and `-window_size` (the epilog prints them when they are not 0).
    pub word_threshold: f64,
    pub window_size: i32,
    /// Subjects in the description table and with alignments, per query.
    pub num_descriptions: usize,
    pub num_alignments: usize,
    /// The queries of the batches without a valid context, which NCBI does not search
    /// (local_blast.cpp:177-208). Indexed like the queries.
    pub unsearched: Vec<bool>,
    /// Whether NCBI's epilog (`PrintEpilog`) ends the report.
    pub epilog: bool,
    /// Whether the report starts with NCBI's prolog (`PrintProlog`); false when the caller
    /// wrote it before the search, as NCBI does.
    pub prolog: bool,
    /// The database sequences with parsed Seq-ids by subject index (empty without them), as
    /// BLASTN's (`BlastnPairwiseReport::parsed_subjects`).
    pub parsed_subjects: Vec<Option<Arc<crate::report::db_subject::ParsedSubject>>>,
    /// The sources of the others' Bioseqs by subject index, as BLASTN's
    /// (`BlastnPairwiseReport::subject_sources`).
    pub subject_sources: Vec<super::defline::SourceTaxname>,
    /// Whether the search linked its HSPs with sum statistics (`-sum_stats`, true by default):
    /// the description table then has its `N` column (`TblastxSections::head`).
    pub sum_stats: bool,
}

/// The start of the TBLASTX report, which NCBI writes before it searches: the program and
/// its reference, then the subjects.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:173-178
pub fn write_tblastx_pairwise_prolog<W: ReportWrite>(
    writer: &mut W,
    version: &str,
    database_name: &[u8],
    database_num_sequences: i64,
    database_total_letters: i64,
    database: Option<&PairwiseDatabase>,
    html: Option<(&super::html_ids::HtmlLinks, &mut dyn Write)>,
) -> io::Result<()> {
    write_prolog(writer, "tblastx", version, false, false, false, html)?;
    if let Some(database) = database {
        return write_blast_database_header(writer, database);
    }
    write_blastp_database_header_spacing(
        writer,
        database_name,
        database_num_sequences,
        database_total_letters,
        1,
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3593-3606
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:310-332
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4229-4230
// The Seq-align has a "sum_n" score only for an HSP of a linked set of two or more
// (blast_seqalign.cpp:1180-1183); without it the display's `sum_n` is -1, which prints no
// count. The counts are those of the displayed rows (`x_FillIdentityInfo`), which the
// caller puts in the record; a TBLASTX HSP has no gaps.
fn write_tblastx_hsp_info<W: Write>(writer: &mut W, hit: &PairwiseHit) -> io::Result<()> {
    let h = &hit.hit;
    let bits = format_bitscore_ncbi(h.bit_score);
    let evalue = format_evalue_ncbi(h.e_value);
    write!(writer, " Score = {bits} bits ({}),  Expect", h.raw_score)?;
    if let Some(sum_n) = hit.sum_n.filter(|&num| num > 1) {
        write!(writer, "({sum_n})")?;
    }
    writeln!(writer, " = {evalue}")?;
    let length = h.length;
    let positives = hit.positives.unwrap_or(h.num_positives);
    let gaps = hit.gaps.unwrap_or(0);
    writeln!(
        writer,
        " Identities = {}/{length} ({}%), Positives = {positives}/{length} ({}%), Gaps = {gaps}/{length} ({}%)",
        h.num_ident,
        ncbi_percent_match(h.num_ident, length),
        ncbi_percent_match(positives, length),
        ncbi_percent_match(gaps, length),
    )?;
    let frame = |frame: i8| {
        if frame > 0 {
            format!("+{frame}")
        } else {
            frame.to_string()
        }
    };
    writeln!(
        writer,
        " Frame = {}/{}",
        frame(hit.query_frame.unwrap_or(1)),
        frame(hit.subject_frame.unwrap_or(1))
    )?;
    writeln!(writer)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1598-1626
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:230,305,358
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1776-1784
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2131-2150
// Both rows are translated (widths 3, `CreateTranslatedDensegFromNADenseg`): a 60-residue line
// spans 180 nucleotides, each row runs in the direction of its frame, and a minus-frame row is
// not reversed (showalign.cpp:1852-1858 applies only without widths; `write_hsp_rows`). The
// query's record has `query_length` nucleotides. The coordinate width is that of the largest
// 0-based nucleotide coordinate (`coordinate_width`). The query row may carry lowercase masks,
// which NCBI applies only when it prints the row (showalign.cpp:2495-2521), so the middle line
// compares the uppercase residues.
fn write_tblastx_alignment<W: Write>(
    writer: &mut W,
    hit: &PairwiseHit,
    masks: Option<&super::query_anchored::HspMasks>,
    query_length: usize,
    config: &PairwiseConfig,
) -> io::Result<()> {
    let (Some(qseq), Some(sseq)) = (&hit.query_seq, &hit.subject_seq) else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing TBLASTX alignment strings",
        ));
    };
    let q = qseq.as_bytes();
    let s = sseq.as_bytes();
    if q.len() != s.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "TBLASTX alignment string lengths differ",
        ));
    }
    let width = coordinate_width(&hit.hit);
    write_hsp_rows(
        writer,
        |segments| tblastx_row_segments(hit, query_length, segments),
        [3, 3],
        config.line_length,
        [q, s],
        masks,
        width,
        |writer, columns| {
            write_protein_middle(
                writer,
                q,
                s,
                columns,
                config.protein_matrix,
                false,
                hit.query_case,
            )
        },
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:451-478
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:583-603
// TBLASTX is ungapped and has no Gumbel block: only the ungapped Karlin block of the
// query's first valid context is printed. A query without a valid context in a searched
// batch prints neither block and a search space of 0 (blast_results.cpp:92-94).
fn write_tblastx_query_footer<W: crate::report::ostream::ReportWrite>(
    writer: &mut W,
    query: &TblastxPairwiseQuery,
) -> io::Result<()> {
    // `x_PrintOneQueryFooter`'s first `NcbiEndl` flushes the stream.
    crate::report::ostream::ncbi_endl(writer)?;
    if let Some(karlin) = query.karlin {
        writeln!(writer, "Lambda      K        H")?;
        write_ncbi_ka_field(writer, karlin.lambda)?;
        write_ncbi_ka_field(writer, karlin.k)?;
        write_ncbi_ka_field(writer, karlin.h)?;
        writeln!(writer)?;
    }
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Effective search space used: {}",
        if query.karlin.is_some() {
            query.effective_search_space
        } else {
            0
        }
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2249,2256-2259
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2266-2288
// TBLASTX is ungapped, so the gap penalties are not printed.
fn write_tblastx_epilog<W: crate::report::ostream::ReportWrite>(
    writer: &mut W,
    report: &TblastxPairwiseReport,
    database: Option<&PairwiseDatabase>,
    html: bool,
) -> io::Result<()> {
    // Two insertions of `NcbiEndl`, each a flush.
    crate::report::ostream::ncbi_endl(writer)?;
    crate::report::ostream::ncbi_endl(writer)?;
    match database {
        Some(database) => write_final_database_block(writer, database)?,
        None => write_final_database_report(
            writer,
            &report.database_name,
            report.database_num_sequences,
            report.database_total_letters,
        )?,
    }
    write_matrix_line(writer, Some(&report.matrix_name))?;
    write_epilog_end(
        writer,
        report.word_threshold,
        report.window_size as i64,
        html,
    )
}

/// A `double` written to a C++ stream with the default flags: `%g` with precision 6
/// (`1e+06` for 1000000, `1.23457e+06` for 1234567).
pub(crate) fn cpp_default_double(value: f64) -> String {
    fn trim(text: &str) -> &str {
        if text.contains('.') {
            text.trim_end_matches('0').trim_end_matches('.')
        } else {
            text
        }
    }
    // glibc's `%g` keeps the sign of a zero (segmasker's `-locut -0` in `algo-options`).
    if value == 0.0 {
        return if value.is_sign_negative() { "-0" } else { "0" }.to_string();
    }
    // glibc's `%g` of an infinity or a NaN.
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.to_string();
    }
    if value.is_nan() {
        return if value.is_sign_negative() {
            "-nan"
        } else {
            "nan"
        }
        .to_string();
    }
    let scientific = format!("{value:.5e}");
    let (mantissa, exponent) = scientific.split_once('e').expect("exponent");
    let exponent: i32 = exponent.parse().expect("exponent digits");
    if !(-4..6).contains(&exponent) {
        let sign = if exponent < 0 { '-' } else { '+' };
        format!("{}e{sign}{:02}", trim(mantissa), exponent.abs())
    } else {
        trim(&format!("{value:.*}", (5 - exponent) as usize)).to_string()
    }
}

/// Writes the TBLASTX pairwise report (outfmt 0) of the queries `query_range` of `queries`.
///
/// `hits` is the final hit list of those queries in its order (queries in input order, then
/// subjects and HSPs in result order); each hit's index in it identifies it to the probe. The
/// prolog comes with the first query of `queries`, the epilog with the last: the reports of
/// consecutive ranges written in turn are the report of all of them. The probe numbers the HSPs
/// over the whole report, so a caller with a probe passes every query (`0..queries.len()`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1491-1589
// Kept out of line: this runs only for outfmt 0.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub fn write_tblastx_pairwise_report<W: crate::report::ostream::ReportWrite>(
    hits: &[PairwiseHit],
    writer: &mut W,
    config: &PairwiseConfig,
    queries: &[TblastxPairwiseQuery],
    query_range: std::ops::Range<usize>,
    subject_titles: &[Arc<[u8]>],
    report: &TblastxPairwiseReport,
    database: Option<&PairwiseDatabase>,
    probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut super::query_warnings::QueryWarnings<'_>>,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
    // The sequences of the query-anchored sections (outfmt 1 to 4, `config.view`; `None` for
    // outfmt 0).
    anchored: Option<&dyn super::query_anchored::ReportSequences>,
    // The rows' masks of `hits` under `-html` (`HtmlRowMasks`; `None` without).
    row_masks: Option<&[super::query_anchored::HspMasks]>,
    loads: &crate::report::db_subject::HitLoads,
) -> io::Result<()> {
    let mut buffered = io::BufWriter::new(writer);
    let writer = &mut buffered;

    // The prolog with the report's first query, the epilog with its last (`query_range` is all of
    // `queries`, or one window of a caller that makes the rows a window at a time).
    if report.prolog && query_range.start == 0 {
        write_tblastx_pairwise_prolog(
            writer,
            &report.version,
            &report.database_name,
            report.database_num_sequences,
            report.database_total_letters,
            database,
            html_with_diagnostics(config.html_links(), warnings.as_deref_mut())?,
        )?;
    }
    let mut hits_by_query: Vec<Vec<(HspIndex, &PairwiseHit)>> = vec![Vec::new(); query_range.len()];
    for (hsp_index, hit) in hits.iter().enumerate() {
        if let Some(bucket) = (hit.hit.q_idx as usize)
            .checked_sub(query_range.start)
            .and_then(|offset| hits_by_query.get_mut(offset))
        {
            bucket.push((hsp_index, hit));
        }
    }
    // The queries' reports, on `x_pool` with LOSAT_X_PAIRWISEPAR (`x_pairwise_par`).
    super::x_pairwise_par::write(
        writer,
        &TblastxSections {
            config,
            queries,
            report,
            database,
            anchored,
            row_masks: HtmlRowMasks::new(row_masks, 0),
            loads,
        },
        query_range.start,
        &hits_by_query,
        subject_titles,
        probe,
        warnings,
        x_pool,
    )?;

    // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:209
    // An error in a later query batch skips the epilog, and nothing follows the last
    // query's footer.
    if report.epilog && query_range.end == queries.len() {
        write_tblastx_epilog(writer, report, database, config.html.is_some())?;
    }
    writer.flush()
}

/// The parts of one query's TBLASTX pairwise report (`x_pairwise_par::QuerySections`).
struct TblastxSections<'a> {
    config: &'a PairwiseConfig,
    queries: &'a [TblastxPairwiseQuery],
    report: &'a TblastxPairwiseReport,
    database: Option<&'a PairwiseDatabase>,
    /// The sequences of the query-anchored sections (outfmt 1 to 4; `None` for outfmt 0).
    anchored: Option<&'a dyn super::query_anchored::ReportSequences>,
    /// The rows' masks of the HSPs under `-html` (outfmt 0).
    row_masks: HtmlRowMasks<'a>,
    /// What the reports' loads do at each subject, by subject (`thrown_load`).
    loads: &'a crate::report::db_subject::HitLoads,
}

impl super::x_pairwise_par::QuerySections for TblastxSections<'_> {
    fn num_alignments(&self) -> usize {
        self.report.num_alignments
    }

    fn line_length(&self) -> usize {
        self.config.line_length
    }

    fn sort(&self) -> (super::hit_sort::ReportSort, bool) {
        (self.config.sort, self.config.tblastx())
    }

    fn joins(&self, previous: u32, s_idx: u32) -> bool {
        // `-subject` records with one id are refused before any report.
        self.database.is_some()
            && crate::report::db_subject::rank_joins_as_copy(
                &self.report.parsed_subjects,
                previous as usize,
                s_idx as usize,
            )
    }

    fn head<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: Option<&super::x_pairwise_par::QueryGroup<'_>>,
        place: super::query_warnings::ReportHook<'_, W>,
    ) -> io::Result<()> {
        let (query, report) = (&self.queries[q_idx], self.report);
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1489-1491
        writer.write_all(b"\n\n")?;
        write_blastp_query_header(
            writer,
            &query.query_name,
            query.query_length,
            self.config.html.is_some(),
        )?;
        let Some(group) = group else {
            write_no_hits_found(writer)?;
            if report.unsearched.get(q_idx).copied().unwrap_or(false) {
                write_tblastn_unsearched_query_footer_spacing(writer, false)?;
            } else {
                write_tblastx_query_footer(writer, query)?;
            }
            return Ok(());
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:606-613
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:165-167
        // NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:157 (the formatter's
        // `use_sum_statistics`) and c++/src/algo/blast/format/blast_format.cpp:503-504
        // The `N` column follows `-sum_stats` (`report.sum_stats`).
        // The titles of the subjects of the records read (by subject index), as BLASTN's.
        write_blastn_description_table(
            writer,
            group.described(report.num_descriptions),
            group.subject_order.len() <= report.num_descriptions,
            &group.subject_hits,
            &group.subject_titles,
            &report.parsed_subjects,
            &report.subject_sources,
            report.sum_stats,
            false,
            self.config.sort_columns(query.query_length),
            place,
            self.config.html_rows(),
            self.loads,
            None,
        )?;
        writeln!(writer)
    }

    fn subject<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        probe: Option<&mut FormatProbe<'_>>,
        group: &super::x_pairwise_par::QueryGroup<'_>,
        s_idx: u32,
        part: std::ops::Range<usize>,
        place: super::query_warnings::ReportHook<'_, W>,
    ) -> io::Result<()> {
        let config = self.config;
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3971-3980
        // The query's record length (`Length=`), also with `-query_loc`: the hits of the query's
        // group are its hits (`x_pairwise_par::write` groups them by `Hit::q_idx`).
        let write_hsp = |writer: &mut W, hit: &PairwiseHit, hsp_index: HspIndex| {
            let query_length = self.queries[hit.hit.q_idx as usize].query_length;
            write_tblastx_hsp_info(writer, hit)?;
            let masks = self.row_masks.of(hsp_index);
            write_tblastx_alignment(writer, hit, masks, query_length, config)?;
            writeln!(writer)
        };
        // A database sequence's heading is its stored title after `>` (as BLASTP's,
        // BLASTN's and TBLASTN's `-db` reports), or the lines of the definition lines of one
        // with parsed Seq-ids.
        let (parsed, sources) = (
            &self.report.parsed_subjects[..],
            &self.report.subject_sources[..],
        );
        let database = self.database.is_some();
        write_subject_part(
            writer,
            probe,
            group,
            s_idx,
            part,
            false,
            config.html.as_ref(),
            parsed,
            None,
            |subject| {
                group_heading(
                    group,
                    parsed,
                    sources,
                    database,
                    false,
                    self.thrown(subject).is_some(),
                    subject,
                )
            },
            place,
            write_hsp,
        )
    }

    fn thrown(&self, s_idx: u32) -> Option<(&crate::ncbi_exception::NcbiException, bool)> {
        thrown_load(self.loads, &self.report.parsed_subjects, s_idx)
    }

    fn view(&self) -> PairwiseView {
        self.config.view
    }

    fn html_links(&self) -> Option<&super::html_ids::HtmlLinks> {
        self.config.html_links()
    }

    fn alignment_id_is_local(&self, s_idx: u32) -> bool {
        alignment_id_local(&self.report.parsed_subjects, s_idx)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1550-1552
    /// The query-anchored section of the query's aligned subjects: Dense-segs of two translated
    /// nucleotide rows in their frames of their records (`seqalign::translated_row_segments`),
    /// grouped by the query's frame and added with `fForceTranslation`
    /// (`eTranslateNucToNucAlignment`, `Hsps::of_translated_nucleotides`).
    fn alignments<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: &super::x_pairwise_par::QueryGroup<'_>,
    ) -> io::Result<()> {
        // The query record's length (`Length=`), also with `-query_loc`.
        let query_length = self.queries[q_idx].query_length;
        write_anchored_section(writer, self.config, self.anchored, q_idx, |sequences| {
            section_hsps(
                super::query_anchored::Hsps::of_translated_nucleotides(),
                self.report.num_alignments,
                group,
                &self.report.parsed_subjects,
                |s_idx| sequences.subject_is_query(q_idx, s_idx),
                |hit, segments| tblastx_row_segments(hit, query_length, segments),
            )
        })
    }

    fn foot<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
    ) -> io::Result<()> {
        write_tblastx_query_footer(writer, &self.queries[q_idx])
    }
}

// =============================================================================
// Main pairwise output function
// =============================================================================

#[cfg(test)]
#[path = "pairwise_tests.rs"]
mod tests;

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:348-445,1490-1589:
// This uses the local-subject database header and TBLASTN translated alignment.
/// The outfmt 0 prolog of tblastn (the version, the references and the database), which
/// NCBI writes before it reads the first query batch (as `write_tblastx_pairwise_prolog`).
/// `program` is the program line's name: `TBLASTN`, or `PSITBLASTN` with `-in_pssm`
/// (`Blast_ProgramNameFromType`, blast_format.cpp:94).
pub fn write_tblastn_pairwise_prolog<W: ReportWrite>(
    writer: &mut W,
    program: &str,
    version: &str,
    database_name: &[u8],
    database_num_sequences: i64,
    database_total_letters: i64,
    database: Option<&PairwiseDatabase>,
    html: Option<(&super::html_ids::HtmlLinks, &mut dyn Write)>,
) -> io::Result<()> {
    write_prolog(writer, program, version, false, false, false, html)?;
    if let Some(database) = database {
        return write_blast_database_header(writer, database);
    }
    write_blastp_database_header_spacing(
        writer,
        database_name,
        database_num_sequences,
        database_total_letters,
        1,
    )
}

/// `query_titles` and `subject_titles` are the title bytes of the records by query and
/// subject index (`Query=`, the description table and the alignment headings; `subject_sources`
/// the sources from which NCBI makes an empty stored title, `stored_source`); the
/// `query_name` of `queries` is not read. `program` is the program line's name (`TBLASTN`, or
/// `PSITBLASTN` with `-in_pssm`).
pub fn write_tblastn_pairwise_report<W: crate::report::ostream::ReportWrite>(
    program: &str,
    hits: &[PairwiseHit],
    writer: &mut W,
    config: &PairwiseConfig,
    queries: &[BlastpPairwiseQuery],
    query_titles: &[Arc<[u8]>],
    query_validity: &[bool],
    query_batch_skipped: &[bool],
    subject_titles: &[Arc<[u8]>],
    parsed_subjects: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    subject_sources: &[super::defline::SourceTaxname],
    report: &BlastpPairwiseReport,
    database: Option<&PairwiseDatabase>,
    probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut super::query_warnings::QueryWarnings<'_>>,
    epilog: bool,
    anchored: Option<&dyn super::query_anchored::ReportSequences>,
    // The rows' masks of `hits` under `-html` (`HtmlRowMasks`; `None` without).
    row_masks: Option<&[super::query_anchored::HspMasks]>,
    loads: &crate::report::db_subject::HitLoads,
    subject_oids: Option<&[u32]>,
) -> io::Result<()> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-224:
    // an all-invalid Run() batch carries -1 Karlin sentinel blocks only for
    // its own queries, even if another batch of the input has valid contexts.
    if query_batch_skipped.len() != queries.len()
        || query_validity.len() != queries.len()
        || query_titles.len() != queries.len()
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "TBLASTN batch validity count mismatch",
        ));
    }
    let mut writer = io::BufWriter::new(writer);
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:372-424:
    write_prolog(
        &mut writer,
        program,
        &report.version,
        false,
        false,
        false,
        html_with_diagnostics(config.html_links(), warnings.as_deref_mut())?,
    )?;
    // The prolog ends with one blank line; each query's preamble starts with two
    // (blast_format.cpp:1491), after the query's warnings.
    match database {
        Some(database) => write_blast_database_header(&mut writer, database)?,
        None => write_blastp_database_header_spacing(
            &mut writer,
            &report.database_name,
            report.database_num_sequences,
            report.database_total_letters,
            1,
        )?,
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411
    // Each hit keeps its index in `hits` (the final result order), which identifies it
    // to the probe.
    let mut hits_by_query: Vec<Vec<(HspIndex, &PairwiseHit)>> = vec![Vec::new(); queries.len()];
    for (hsp_index, hit) in hits.iter().enumerate() {
        if let Some(bucket) = hits_by_query.get_mut(hit.hit.q_idx as usize) {
            bucket.push((hsp_index, hit));
        }
    }
    // The queries' reports, by the shared writer of the pairwise reports (`x_pairwise_par`; its
    // serial writer: TBLASTN gives it no pool).
    super::x_pairwise_par::write(
        &mut writer,
        &TblastnSections {
            config,
            queries,
            query_titles,
            query_validity,
            query_batch_skipped,
            parsed_subjects,
            subject_sources,
            report,
            database,
            anchored,
            row_masks: HtmlRowMasks::new(row_masks, 0),
            loads,
            subject_oids,
        },
        0,
        &hits_by_query,
        subject_titles,
        probe,
        warnings,
        None,
    )?;
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:342
    // A search that stops at a query batch (`Empty CBlastQueryVector` with `-query_loc`)
    // writes no epilog.
    if !epilog {
        return writer.flush();
    }
    // The blank lines before the epilog (blast_format.cpp:2249), each a flush.
    crate::report::ostream::ncbi_endl(&mut writer)?;
    crate::report::ostream::ncbi_endl(&mut writer)?;
    write_blastp_final_footer(&mut writer, report, database, config.html.is_some())?;
    writer.flush()
}

/// The parts of one query's TBLASTN pairwise report (`x_pairwise_par::QuerySections`).
struct TblastnSections<'a> {
    config: &'a PairwiseConfig,
    queries: &'a [BlastpPairwiseQuery],
    /// The title bytes of the records by query index (`Query=`).
    query_titles: &'a [Arc<[u8]>],
    query_validity: &'a [bool],
    query_batch_skipped: &'a [bool],
    parsed_subjects: &'a [Option<Arc<crate::report::db_subject::ParsedSubject>>],
    /// What NCBI makes an empty stored title of, by subject (`stored_source`).
    subject_sources: &'a [super::defline::SourceTaxname],
    report: &'a BlastpPairwiseReport,
    database: Option<&'a PairwiseDatabase>,
    /// The sequences of the query-anchored sections (outfmt 1 to 4; `None` for outfmt 0).
    anchored: Option<&'a dyn super::query_anchored::ReportSequences>,
    /// The rows' masks of the HSPs under `-html` (outfmt 0).
    row_masks: HtmlRowMasks<'a>,
    /// What the reports' loads do at each subject, by subject (`thrown_load`).
    loads: &'a crate::report::db_subject::HitLoads,
    /// Database OIDs by subject rank (`None` for `-subject`).
    subject_oids: Option<&'a [u32]>,
}

impl super::x_pairwise_par::QuerySections for TblastnSections<'_> {
    fn num_alignments(&self) -> usize {
        self.report.num_alignments
    }

    fn line_length(&self) -> usize {
        self.config.line_length
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1522-1530
    // `-sorthits` and `-sorthsps` before the table and the alignments (`QueryGroup::sort`).
    fn sort(&self) -> (super::hit_sort::ReportSort, bool) {
        (self.config.sort, self.config.tblastx())
    }

    // The subjects' hits; a copy of one record that follows another is shown in its hit
    // (`QueryGroup::new`, `db_subject::joins_as_copy`).
    fn joins(&self, previous: u32, s_idx: u32) -> bool {
        // `-subject` records with one id are refused before any report.
        self.database.is_some()
            && crate::report::db_subject::rank_joins_as_copy(
                self.parsed_subjects,
                previous as usize,
                s_idx as usize,
            )
    }

    fn head<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: Option<&super::x_pairwise_par::QueryGroup<'_>>,
        place: super::query_warnings::ReportHook<'_, W>,
    ) -> io::Result<()> {
        let query = &self.queries[q_idx];
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1491
        // Emit two blank lines before the query header.
        writeln!(writer)?;
        writeln!(writer)?;
        write_blastp_query_header(
            writer,
            &self.query_titles[q_idx],
            query.query_length,
            self.config.html.is_some(),
        )?;
        let Some(group) = group else {
            write_no_hits_found(writer)?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-97 leaves both Karlin blocks null
            // when there is no valid context. blast_format.cpp:445-477
            // writes only the blank lines and zero search space in that case.
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180,204-208:
            // The all-invalid batch uses -1 sentinel blocks. An invalid query
            // within a searched batch has null blocks (blast_results.cpp:82-103).
            if self.query_batch_skipped[q_idx] {
                write_tblastn_unsearched_query_footer_spacing(writer, false)?;
            } else if !self.query_validity[q_idx] {
                // `x_PrintOneQueryFooter`'s first `NcbiEndl` flushes the stream.
                crate::report::ostream::ncbi_endl(writer)?;
                writeln!(writer)?;
                writeln!(writer)?;
                writeln!(writer, "Effective search space used: 0")?;
            } else {
                write_translated_query_footer(
                    writer,
                    query,
                    self.report,
                    self.report.gapped_karlin,
                    self.report.gapped,
                )?;
            }
            return Ok(());
        };
        // NCBI's sorts read the HSPs' ranges before the table (`sort_stop`): a query row that
        // starts with a gap (a Smith-Waterman alignment of `-use_sw_tback`, `-db`: `-subject`'s
        // records come mapped) throws there.
        let (grouped, order) = (&group.subject_hits, &group.subject_order);
        let remaps = self.database.is_none();
        if let Some(stop) = sort_stop(self.config.sort, remaps, order, grouped) {
            writer.flush()?;
            return Err(io::Error::other(stop));
        }
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1004-1011
        // The table of every program follows `x_InitDeflineTable`
        // (`write_blastn_description_table`): the subject's highest bit score and that HSP's
        // E-value, and the title of the subject (`CDeflineGenerator`), for the first
        // `m_NumDescriptions` subjects. The titles are those of the records read (by
        // subject index), as BLASTN's and TBLASTX's. It reads the ranges of every HSP of each
        // described subject first, as BLASTX's (`described_row_range`).
        let descriptions = group.described(self.report.num_descriptions);
        let rows_throw = described_row_range(remaps, descriptions, grouped)
            .map(EmptyRowsStop::RowRange)
            .or_else(|| {
                descriptions
                    .iter()
                    .any(|oid| {
                        grouped[oid]
                            .iter()
                            .any(|hit| hit.edge == PairwiseEdge::NullbStdSegLengths)
                    })
                    .then_some(EmptyRowsStop::StdSegLengths)
            });
        write_blastn_description_table(
            writer,
            descriptions,
            order.len() <= self.report.num_descriptions,
            grouped,
            &group.subject_titles,
            self.parsed_subjects,
            self.subject_sources,
            self.report.show_linked_set_size,
            false,
            self.config.sort_columns(query.query_length),
            place,
            self.config.html_rows(),
            self.loads,
            rows_throw,
        )?;
        writeln!(writer)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1014-1040
    // (PruneSeqalign keeps the alignments of the first `m_NumAlignments` subjects).
    fn subject<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        probe: Option<&mut FormatProbe<'_>>,
        group: &super::x_pairwise_par::QueryGroup<'_>,
        s_idx: u32,
        part: std::ops::Range<usize>,
        place: super::query_warnings::ReportHook<'_, W>,
    ) -> io::Result<()> {
        let config = self.config;
        // A PSSM query whose letters the display cannot read (`UnreadableQuery`): NCBI's line in
        // place of each alignment.
        if let Some(id) = config.unreadable_query_id(self.parsed_subjects, s_idx) {
            return write_database_subject_alignments_part(
                writer,
                probe,
                DatabaseHeading::Skipped(&id),
                false,
                self.parsed_subjects,
                None,
                s_idx,
                &group.subject_hits[&s_idx],
                &group.subject_hit_indices[&s_idx],
                part,
                None,
                place,
                |_, _, _| unreachable!("a skipped alignment is written by its line"),
            );
        }
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3650-3668:
        // display leaves an extra blank line after each translated HSP.
        let write_hsp = |writer: &mut W, hit: &PairwiseHit, hsp_index: HspIndex| {
            write_tblastn_hsp_info(writer, hit)?;
            write_tblastn_alignment(writer, hit, self.row_masks.of(hsp_index), config)?;
            writeln!(writer)
        };
        // A database sequence's heading is its stored title after `>` (as BLASTP's
        // and BLASTN's `-db` reports), or the lines of the definition lines of one with
        // parsed Seq-ids.
        let (parsed, sources) = (self.parsed_subjects, self.subject_sources);
        let database = self.database.is_some();
        let hits = &group.subject_hits[&s_idx];
        let stop = part
            .clone()
            .find(|&hsp| hits[hsp].edge == PairwiseEdge::NullbStdSegLengths);
        let written = part.start..stop.unwrap_or(part.end);
        let heading = |subject| {
            group_heading(
                group,
                parsed,
                sources,
                database,
                self.report.gapped,
                self.thrown(subject).is_some(),
                subject,
            )
        };
        if !written.is_empty() {
            write_subject_part(
                writer,
                probe,
                group,
                s_idx,
                written,
                false,
                config.html.as_ref(),
                parsed,
                self.subject_oids,
                heading,
                place,
                write_hsp,
            )?;
        }
        if stop.is_some() {
            writer.flush()?;
            return Err(io::Error::other(EmptyRowsStop::StdSegLengths));
        }
        Ok(())
    }

    fn thrown(&self, s_idx: u32) -> Option<(&crate::ncbi_exception::NcbiException, bool)> {
        thrown_load(self.loads, self.parsed_subjects, s_idx)
    }

    fn view(&self) -> PairwiseView {
        self.config.view
    }

    fn html_links(&self) -> Option<&super::html_ids::HtmlLinks> {
        self.config.html_links()
    }

    fn alignment_id_is_local(&self, s_idx: u32) -> bool {
        alignment_id_local(self.parsed_subjects, s_idx)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1550-1552
    /// The query-anchored section of the query's aligned subjects: Dense-segs of a protein query
    /// row and a translated nucleotide subject row (`seqalign::translated_row_segments`).
    fn alignments<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: &super::x_pairwise_par::QueryGroup<'_>,
    ) -> io::Result<()> {
        let stops = group
            .subject_order
            .iter()
            .take(self.report.num_alignments)
            .any(|s_idx| {
                group.subject_hits[s_idx]
                    .iter()
                    .any(|hit| hit.edge == PairwiseEdge::NullbStdSegLengths)
            });
        if stops {
            writer.flush()?;
            return Err(io::Error::other(EmptyRowsStop::StdSegLengths));
        }
        write_anchored_section(writer, self.config, self.anchored, q_idx, |_| {
            section_hsps(
                super::query_anchored::Hsps::of_translated_subjects(),
                self.report.num_alignments,
                group,
                self.parsed_subjects,
                |_| false,
                tblastn_row_segments,
            )
        })
    }

    fn foot<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
    ) -> io::Result<()> {
        // An ungapped search's footer has the ungapped Karlin block only (`x_PrintOneQueryFooter`,
        // cited at `write_translated_query_footer`).
        write_translated_query_footer(
            writer,
            &self.queries[q_idx],
            self.report,
            self.report.gapped_karlin,
            self.report.gapped,
        )
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3578-3604,2122-2149:
// Format the identities, positives, and gaps summary line.
// TBLASTN reports one translated subject frame, without a query frame.
fn write_tblastn_hsp_info<W: Write>(writer: &mut W, hit: &PairwiseHit) -> io::Result<()> {
    let h = &hit.hit;
    let bits = format_bitscore_ncbi(h.bit_score);
    let evalue = format_evalue_ncbi(h.e_value);
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3599-3604:
    // Append the composition-adjustment method suffix for methods 1 and 2.
    let method = match hit.comp_adjust_method {
        Some(1) => ", Method: Composition-based stats.",
        Some(2) => ", Method: Compositional matrix adjust.",
        _ => "",
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1180-1182:
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3595-3598:
    let expect = expect_label(hit.sum_n);
    writeln!(
        writer,
        " Score = {} bits ({}),  {} = {}{}",
        bits, h.raw_score, expect, evalue, method
    )?;
    let length = h.length;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161:
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
    let percent = |n: usize| ncbi_percent_match(n, length);
    let positives = hit.positives.unwrap_or(h.num_positives);
    let gaps = hit.query_range_lengths().gaps;
    writeln!(
        writer,
        " Identities = {}/{} ({:.0}%), Positives = {}/{} ({:.0}%), Gaps = {}/{} ({:.0}%)",
        h.num_ident,
        length,
        percent(h.num_ident),
        positives,
        length,
        percent(positives),
        gaps,
        length,
        percent(gaps)
    )?;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:321-332
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:321-331
    // The frames come from the alignment's rows: an HSP with empty rows (BLASTX's ungapped HSPs
    // of length 0 with `-db`, `PairwiseHit::rows_empty`) has none, and no `Frame` line.
    if let Some(frame) = hit
        .query_frame
        .or(hit.subject_frame)
        .filter(|_| !hit.rows_empty())
    {
        if frame > 0 {
            writeln!(writer, " Frame = +{frame}")?;
        } else {
            writeln!(writer, " Frame = {frame}")?;
        }
    }
    writeln!(writer)?;
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1600-1625,2122-2149:
// The query's row of letters and the subject's translated row (width 3, three nucleotides per
// residue, in the direction of its frame; `write_hsp_rows`).
fn write_tblastn_alignment<W: Write>(
    writer: &mut W,
    hit: &PairwiseHit,
    masks: Option<&super::query_anchored::HspMasks>,
    config: &PairwiseConfig,
) -> io::Result<()> {
    let (Some(qseq), Some(sseq)) = (&hit.query_seq, &hit.subject_seq) else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing TBLASTN alignment strings",
        ));
    };
    let q = qseq.as_bytes();
    let s = sseq.as_bytes();
    if q.len() != s.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "TBLASTN alignment string lengths differ",
        ));
    }
    let width = coordinate_width(&hit.hit);
    let rows = ascii_row_strings([q, s], hit.hit.gap_info.as_deref().unwrap_or_default());
    write_hsp_rows(
        writer,
        |segments| tblastn_row_segments(hit, segments),
        [1, 3],
        config.line_length,
        [&rows[0], &rows[1]],
        masks,
        width,
        |writer, columns| {
            write_protein_middle(
                writer,
                q,
                s,
                columns,
                config.protein_matrix,
                false,
                hit.query_case,
            )
        },
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:275-282
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1548-1565
pub struct BlastxPairwiseOptions {
    pub word_threshold: f64,
    pub gapped: bool,
    pub sum_stats: bool,
    pub num_descriptions: usize,
    pub num_alignments: usize,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1490-1519
/// The subjects of a BLASTX `-db` report (`write_blastx_pairwise_report`): the database, and by
/// subject index the stored titles, the parsed subjects and the sources of the empty titles
/// (`stored_source`).
pub type BlastxDatabaseSubjects<'a> = (
    &'a PairwiseDatabase,
    &'a [Arc<[u8]>],
    &'a [Option<Arc<crate::report::db_subject::ParsedSubject>>],
    &'a [super::defline::SourceTaxname],
);

/// Writes the queries' part of the BLASTX pairwise report (outfmt 0, and the query-anchored 1 to 4
/// of `config.view`) of `queries`, whose hits are `hits` in the order of the reports (`q_idx` is
/// the index in `queries`), by the shared writer of the pairwise reports (`x_pairwise_par::write`;
/// on `pool` without a `probe`): each query's warnings and the texts of its subjects' reads
/// (`warnings`, `QueryWarnings::before_query` and `subject_read`), its head, its alignments and its
/// footer (`BlastxSections`). The prolog and the epilog are the run's (`write_blastx_prolog`,
/// `write_blastx_epilog`). `probe` marks each subject's heading and each HSP, `hits[i]` being HSP
/// `first_hsp + i`, without changing the written bytes. `database` holds a `-db` report's subjects;
/// a `-subject` report's are `records_parsed` and `record_titles` (by subject index).
#[allow(clippy::too_many_arguments)]
pub fn write_blastx_pairwise_report<W: crate::report::ostream::ReportWrite>(
    hits: &[PairwiseHit],
    writer: &mut W,
    mut warnings: Option<&mut super::query_warnings::QueryWarnings<'_>>,
    config: &PairwiseConfig,
    queries: &[BlastpPairwiseQuery],
    subject_labels: &[DescriptionLabel],
    report: &BlastpPairwiseReport,
    options: &BlastxPairwiseOptions,
    database: Option<BlastxDatabaseSubjects<'_>>,
    records_parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    record_titles: &[&[u8]],
    // Each query's gapped K where it is not the report's (`BatchResults::gapped_k`).
    query_gapped_k: Option<&[f64]>,
    probe: Option<&mut FormatProbe<'_>>,
    first_hsp: HspIndex,
    // The sequences of the query-anchored sections (outfmt 1 to 4, `config.view`; `None` for
    // outfmt 0).
    anchored: Option<&dyn super::query_anchored::ReportSequences>,
    pool: Option<&crate::utils::threading::SearchPool<'_>>,
    // The rows' masks of `hits` under `-html` (`HtmlRowMasks`; `None` without).
    row_masks: Option<&[super::query_anchored::HspMasks]>,
    loads: &crate::report::db_subject::HitLoads,
) -> io::Result<()> {
    // The report's insertions reach the stream in chunks, as in the other pairwise writers: each
    // region change of a pairwise report is at a flush of NCBI's, which writes the chunk through.
    let mut buffered = io::BufWriter::new(writer);
    let writer = &mut buffered;
    let sections = BlastxSections {
        config,
        options,
        queries,
        subject_labels,
        report,
        query_gapped_k,
        database,
        records_parsed,
        record_titles,
        anchored,
        row_masks: HtmlRowMasks::new(row_masks, first_hsp),
        loads,
    };
    // `hits[i]` is HSP `first_hsp + i` of the run's hit list (LOSAT Web's `run_local` writes a run
    // batch by batch).
    let mut by_query: Vec<Vec<(HspIndex, &PairwiseHit)>> = vec![Vec::new(); queries.len()];
    for (index, hit) in hits.iter().enumerate() {
        by_query[hit.hit.q_idx as usize].push((first_hsp + index, hit));
    }
    // A query-anchored section (outfmt 1 to 4) whose merge NCBI BLAST+ frees a segment in is
    // refused before the call's first byte, unless the report of a query before it ends with
    // NCBI's exception (`blastx_anchored_stop`); only HSPs of length 0 get there.
    if config.view != PairwiseView::Pairwise && hits.iter().any(PairwiseHit::rows_empty) {
        for query_hits in &by_query {
            let Some(group) =
                super::x_pairwise_par::QueryGroup::for_sections(query_hits, &[], &sections)
            else {
                continue;
            };
            match blastx_anchored_stop(config, options, &group, database.is_none())? {
                Some(EmptyRowsStop::MergeFreesSegment) => {
                    return Err(io::Error::other(EmptyRowsStop::MergeFreesSegment))
                }
                Some(_) => break,
                None => {}
            }
        }
    }
    // The queries' reports, by the shared writer of the pairwise reports (`x_pairwise_par`). The
    // titles are the sections' own (`BlastxSections`), so the groups copy none. A report that
    // stops (NCBI's exception at an HSP of length 0, `EmptyRowsStop`; a failed write) ends the
    // run there: the reads after it make no text (`x_pairwise_par::write`).
    super::x_pairwise_par::write(
        writer,
        &sections,
        0,
        &by_query,
        &[],
        probe,
        warnings.as_deref_mut(),
        pool,
    )?;
    // The texts of the reads that the reports did not place come after them
    // (`QueryWarnings::after_last_query`).
    if let Some(warnings) = warnings {
        warnings.after_last_query(queries.len(), writer)?;
    }
    writer.flush()
}

/// The parts of one query's BLASTX pairwise report (`x_pairwise_par::QuerySections`), RPSTBLASTN's
/// too (`write_blastx_pairwise_report`).
struct BlastxSections<'a> {
    config: &'a PairwiseConfig,
    options: &'a BlastxPairwiseOptions,
    queries: &'a [BlastpPairwiseQuery],
    /// The description line of each subject, by subject index (`database_description_label`).
    subject_labels: &'a [DescriptionLabel],
    report: &'a BlastpPairwiseReport,
    /// Each query's gapped K where it is not the report's (an RPSTBLASTN search with
    /// composition-based statistics, `BatchResults::gapped_k`), as rpsblast's.
    query_gapped_k: Option<&'a [f64]>,
    /// A `-db` report's subjects (`None` for `-subject`).
    database: Option<BlastxDatabaseSubjects<'a>>,
    /// A `-subject` report's records with parsed Seq-ids, by subject index.
    records_parsed: &'a [Option<Arc<crate::report::db_subject::ParsedSubject>>],
    /// A `-subject` report's record titles (their headings' `GenerateDefline`), by subject index.
    record_titles: &'a [&'a [u8]],
    /// The sequences of the query-anchored sections (outfmt 1 to 4; `None` for outfmt 0).
    anchored: Option<&'a dyn super::query_anchored::ReportSequences>,
    /// The rows' masks of the HSPs under `-html` (outfmt 0).
    row_masks: HtmlRowMasks<'a>,
    /// What the reports' loads do at each database subject, by subject (`thrown_load`).
    loads: &'a crate::report::db_subject::HitLoads,
}

impl BlastxSections<'_> {
    /// Whether `-subject`'s Seq-loc mapper rewrites the Std-segs of the HSPs of length 0
    /// (`mapper_rewrites`; a `-db` report keeps them as the engine made them).
    fn remaps(&self) -> bool {
        self.database.is_none()
    }

    /// The alignments of the HSPs `part` of subject `s_idx` of `group` (the heading with the
    /// first), each HSP's score lines, rows and blank line (`write_tblastn_hsp_info`,
    /// `write_blastx_alignment`): a database sequence's heading is its stored title after `>` (as
    /// BLASTP's, BLASTN's and TBLASTN's `-db` reports), or the lines of the definition lines of a
    /// sequence with parsed Seq-ids (also a `-subject` record's), or NCBI's line in place of each
    /// alignment where it cannot make the heading's title
    /// (`write_database_subject_alignments_part`); a `-subject` record's is the title that
    /// `GenerateDefline` makes (flags 0, `write_subject_alignments_part`).
    fn write_part<W: Write>(
        &self,
        writer: &mut W,
        probe: Option<&mut FormatProbe<'_>>,
        group: &super::x_pairwise_par::QueryGroup<'_>,
        s_idx: u32,
        part: std::ops::Range<usize>,
        place: super::query_warnings::ReportHook<'_, W>,
    ) -> io::Result<()> {
        let config = self.config;
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3650-3668:
        // display leaves an extra blank line after each translated HSP.
        let write_hsp = |writer: &mut W, hit: &PairwiseHit, hsp_index: HspIndex| {
            write_tblastn_hsp_info(writer, hit)?;
            write_blastx_alignment(writer, hit, self.row_masks.of(hsp_index), config)?;
            writeln!(writer)
        };
        let heading_of = |s_idx: u32| {
            let parsed = match self.database {
                Some((_, _, parsed, _)) => parsed_subject(parsed, s_idx),
                None => parsed_subject(self.records_parsed, s_idx),
            };
            match (parsed, self.database) {
                (Some(subject), _) => SubjectHeading::Database(DatabaseHeading::loaded(
                    subject,
                    self.options.gapped,
                    group.own_hsps(s_idx),
                    super::x_pairwise_par::QuerySections::thrown(self, s_idx).is_some(),
                )),
                (None, Some((_, titles, _, sources))) => SubjectHeading::Database(
                    DatabaseHeading::Title(&titles[s_idx as usize], stored_source(sources, s_idx)),
                ),
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2273
                (None, None) => SubjectHeading::Record(
                    self.record_titles
                        .get(s_idx as usize)
                        .copied()
                        .unwrap_or_default(),
                ),
            }
        };
        write_subject_part(
            writer,
            probe,
            group,
            s_idx,
            part,
            true,
            config.html.as_ref(),
            self.records_parsed,
            None,
            heading_of,
            place,
            write_hsp,
        )
    }
}

impl BlastxSections<'_> {
    /// The subjects with parsed Seq-ids by subject index: the database's, or the `-subject`
    /// records'.
    fn database_parsed(&self) -> &[Option<Arc<crate::report::db_subject::ParsedSubject>>] {
        match self.database {
            Some((_, _, parsed, _)) => parsed,
            None => self.records_parsed,
        }
    }
}

impl super::x_pairwise_par::QuerySections for BlastxSections<'_> {
    fn num_alignments(&self) -> usize {
        self.options.num_alignments
    }

    fn line_length(&self) -> usize {
        self.config.line_length
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1522-1530
    // `-sorthits` and `-sorthsps` before the table and the alignments (`QueryGroup::sort`).
    fn sort(&self) -> (super::hit_sort::ReportSort, bool) {
        (self.config.sort, self.config.tblastx())
    }

    // A copy of one record of the database that follows another is shown in its hit
    // (`QueryGroup::new`, `db_subject::joins_as_copy`); `-subject` records with one id are
    // refused before any report.
    fn joins(&self, previous: u32, s_idx: u32) -> bool {
        self.database.is_some_and(|(_, _, parsed, _)| {
            crate::report::db_subject::rank_joins_as_copy(parsed, previous as usize, s_idx as usize)
        })
    }

    fn head<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: Option<&super::x_pairwise_par::QueryGroup<'_>>,
        place: super::query_warnings::ReportHook<'_, W>,
    ) -> io::Result<()> {
        let query = &self.queries[q_idx];
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1491-1492
        writeln!(writer)?;
        writeln!(writer)?;
        write_blastp_query_header(
            writer,
            &query.query_name,
            query.query_length,
            self.config.html.is_some(),
        )?;
        let Some(group) = group else {
            write_no_hits_found(writer)?;
            return self.foot(writer, q_idx);
        };
        let (grouped, order) = (&group.subject_hits, &group.subject_order);
        // NCBI's sorts read the HSPs' ranges before the table (`sort_stop`).
        if let Some(stop) = sort_stop(self.config.sort, self.remaps(), order, grouped) {
            writer.flush()?;
            return Err(io::Error::other(stop));
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:606-611
        // The table follows `x_InitDeflineTable` (`write_description_table`): the subject's
        // highest bit score and that HSP's E-value, for the first `m_NumDescriptions` subjects;
        // the last one is NCBI's "last hit" when every subject fits.
        let descriptions = group.described(self.options.num_descriptions);
        // The table reads the ranges of every HSP of each described subject (after its `-html`
        // read of `./.ncbirc`): a row whose range it cannot take stops the report after the
        // query's `Length=` (`described_row_range`).
        let rows_throw =
            described_row_range(self.remaps(), descriptions, grouped).map(EmptyRowsStop::RowRange);
        write_subject_summary_table_with_sum_n(
            writer,
            descriptions,
            order.len() <= self.options.num_descriptions,
            grouped,
            self.subject_labels,
            !self.options.gapped && self.options.sum_stats,
            self.config.sort_columns(query.query_length),
            place,
            self.config.html_rows().map(|report| {
                let parsed = match self.database {
                    Some((_, _, parsed, _)) => parsed,
                    None => self.records_parsed,
                };
                (report, parsed, self.loads)
            }),
            rows_throw,
        )?;
        writeln!(writer)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1948-1973
    // Each HSP's alignment is made before its heading and block: the Std-seg of an HSP with empty
    // rows that `-subject`'s Seq-loc mapper rewrote throws there (`CreateDensegFromStdseg`), before
    // anything of the HSP (`EmptyRowsStop::StdSegLengths`), the subject's heading for its first.
    fn subject<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        probe: Option<&mut FormatProbe<'_>>,
        group: &super::x_pairwise_par::QueryGroup<'_>,
        s_idx: u32,
        part: std::ops::Range<usize>,
        place: super::query_warnings::ReportHook<'_, W>,
    ) -> io::Result<()> {
        let hits = &group.subject_hits[&s_idx];
        let stop = mapper_rewrites(self.remaps(), hits)
            .then(|| part.clone().find(|&hsp| hits[hsp].rows_empty()))
            .flatten();
        let written = part.start..stop.unwrap_or(part.end);
        if !written.is_empty() {
            self.write_part(writer, probe, group, s_idx, written, place)?;
        }
        if stop.is_some() {
            writer.flush()?;
            return Err(io::Error::other(EmptyRowsStop::StdSegLengths));
        }
        Ok(())
    }

    fn thrown(&self, s_idx: u32) -> Option<(&crate::ncbi_exception::NcbiException, bool)> {
        let parsed = self.database.map_or(&[][..], |(_, _, parsed, _)| parsed);
        thrown_load(self.loads, parsed, s_idx)
    }

    fn view(&self) -> PairwiseView {
        self.config.view
    }

    fn html_links(&self) -> Option<&super::html_ids::HtmlLinks> {
        self.config.html_links()
    }

    fn alignment_id_is_local(&self, s_idx: u32) -> bool {
        alignment_id_local(self.database_parsed(), s_idx)
    }

    // Each HSP's Dense-seg is made after its subject's read, before the section's first byte and
    // outside its `try` (with `-html` in outfmt 0, before the first alignment): the Std-seg of an
    // HSP with empty rows that `-subject`'s Seq-loc mapper rewrote throws there
    // (`EmptyRowsStop::StdSegLengths`, as outfmt 0's plain alignments at the HSP).
    fn subject_densegs<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        group: &super::x_pairwise_par::QueryGroup<'_>,
        s_idx: u32,
    ) -> io::Result<()> {
        if blastx_mapped_empty_rows(self.remaps(), &group.subject_hits[&s_idx]) {
            writer.flush()?;
            return Err(io::Error::other(EmptyRowsStop::StdSegLengths));
        }
        Ok(())
    }

    // Only `-subject`'s Dense-segs throw (`remaps`, the test of `subject_densegs`), so a `-db`
    // report gives no texts and makes none: where `subject_densegs` throws, this gives texts.
    fn alignment_id_texts(&self, s_idx: u32) -> Option<[Vec<u8>; 2]> {
        self.remaps()
            .then(|| record_alignment_id_texts(self.records_parsed, s_idx))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1550-1552
    /// The query-anchored section of the query's aligned subjects (`blastx_section_hsps`).
    fn alignments<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: &super::x_pairwise_par::QueryGroup<'_>,
    ) -> io::Result<()> {
        write_anchored_section(writer, self.config, self.anchored, q_idx, |_| {
            blastx_section_hsps(self.options.num_alignments, group, self.database_parsed())
        })
    }

    fn foot<W: crate::report::ostream::ReportWrite>(
        &self,
        writer: &mut W,
        q_idx: usize,
    ) -> io::Result<()> {
        write_translated_query_footer(
            writer,
            &self.queries[q_idx],
            self.report,
            query_gapped_karlin(self.report, self.query_gapped_k, q_idx),
            self.options.gapped,
        )
    }
}
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254-254
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:297-297
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:387-388
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:406-407
/// `program` is the program line's name: `BLASTX`, or `RPSTBLASTN`, whose prolog is BLASTX's (the
/// Gapped BLAST reference only, blast_format.cpp:416-426, `write_prolog`).
pub fn write_blastx_prolog(
    writer: &mut impl ReportWrite,
    program: &str,
    report: &BlastpPairwiseReport,
    database: Option<&PairwiseDatabase>,
    html: Option<(&super::html_ids::HtmlLinks, &mut dyn Write)>,
) -> io::Result<()> {
    write_prolog(writer, program, &report.version, false, false, false, html)?;
    if let Some(database) = database {
        return write_blast_database_header(writer, database);
    }
    write_blastp_database_header_spacing(
        writer,
        &report.database_name,
        report.database_num_sequences,
        report.database_total_letters,
        1,
    )
}
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2264-2288
pub fn write_blastx_epilog(
    writer: &mut impl crate::report::ostream::ReportWrite,
    report: &BlastpPairwiseReport,
    options: &BlastxPairwiseOptions,
    database: Option<&PairwiseDatabase>,
    html: bool,
) -> io::Result<()> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2264-2288
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2249-2249
    crate::report::ostream::ncbi_endl(writer)?;
    crate::report::ostream::ncbi_endl(writer)?;
    if let Some(database) = database {
        write_final_database_block(writer, database)?;
    } else {
        write!(writer, "  Database: ")?;
        write_flatfile_wrapped(writer, &ensure_trailing_period(&report.database_name), 68)?;
        writeln!(writer, "    Posted date:  Unknown")?;
        writeln!(
            writer,
            "  Number of letters in database: {}",
            crate::ncbi_str::int8_with_commas(report.database_total_letters)
        )?;
        writeln!(
            writer,
            "  Number of sequences in database:  {}",
            crate::ncbi_str::int_with_commas(report.database_num_sequences as i32)
        )?;
        writeln!(writer)?;
        writeln!(writer)?;
        writeln!(writer)?;
    }
    write_matrix_line(writer, report.matrix_name.as_deref())?;
    if options.gapped {
        writeln!(
            writer,
            "Gap Penalties: Existence: {}, Extension: {}",
            report.gap_open, report.gap_extend
        )?;
    }
    write_epilog_end(
        writer,
        options.word_threshold,
        report.window_size as i64,
        html,
    )
}
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:445-478
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:583-613
/// The footer of a BLASTX or TBLASTN query (`BlastpPairwiseQuery::valid` and `batch_skipped`: its
/// Karlin blocks, `gapped`: the gapped one, `gapped_karlin`: the query's, `query_gapped_karlin`).
fn write_translated_query_footer(
    writer: &mut impl crate::report::ostream::ReportWrite,
    query: &BlastpPairwiseQuery,
    report: &BlastpPairwiseReport,
    gapped_karlin: KarlinParams,
    gapped: bool,
) -> io::Result<()> {
    let (valid, skipped) = (query.valid, query.batch_skipped);
    if skipped {
        return write_tblastn_unsearched_query_footer_spacing(writer, false);
    }
    if valid && gapped {
        return write_blastp_query_footer_spacing(
            writer,
            query.ungapped_karlin,
            gapped_karlin,
            report.gumbel,
            true,
            query.effective_search_space,
            false,
        );
    }
    // `x_PrintOneQueryFooter`'s first `NcbiEndl` flushes the stream.
    crate::report::ostream::ncbi_endl(writer)?;
    if valid {
        writeln!(writer, "Lambda      K        H")?;
        for v in [
            query.ungapped_karlin.lambda,
            query.ungapped_karlin.k,
            query.ungapped_karlin.h,
        ] {
            write_ncbi_ka_field(writer, v)?;
        }
        writeln!(writer)?;
    }
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Effective search space used: {}",
        if valid {
            query.effective_search_space
        } else {
            0
        }
    )?;
    Ok(())
}
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1596-1627
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:120-134
fn write_blastx_alignment(
    writer: &mut impl Write,
    hit: &PairwiseHit,
    masks: Option<&super::query_anchored::HspMasks>,
    config: &PairwiseConfig,
) -> io::Result<()> {
    let (Some(qseq), Some(sseq)) = (&hit.query_seq, &hit.subject_seq) else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing BLASTX alignment strings",
        ));
    };
    let q = qseq.as_bytes();
    let s = sseq.as_bytes();
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1364-1367
    let width = digit_count(
        hit.hit
            .q_start
            .max(hit.hit.q_end)
            .max(hit.hit.s_start)
            .max(hit.hit.s_end)
            - 1,
    );
    // The query's translated row (width 3: three nucleotides per residue, in the direction of its
    // frame) and the subject's row of letters, whose stored gap characters are residues of the
    // Dense-seg (`write_hsp_rows`); their end gaps are spaces (`ascii_row_strings`).
    if q.len() != s.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "BLASTX alignment string lengths differ",
        ));
    }
    let rows = ascii_row_strings([q, s], hit.hit.gap_info.as_deref().unwrap_or_default());
    write_hsp_rows(
        writer,
        |segments| blastx_row_segments(hit, segments),
        [3, 1],
        config.line_length,
        [&rows[0], &rows[1]],
        masks,
        width,
        |writer, columns| {
            write_protein_middle(
                writer,
                q,
                s,
                columns,
                config.protein_matrix,
                true,
                hit.query_case,
            )
        },
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
#[cfg(test)]
#[test]
fn percent_match_matches_independent_pinned_cpp_edges() {
    let expected = include_str!("../../tests/unit/blastx_stage_e_percent_expected.tsv");
    for line in expected.lines() {
        let fields: Vec<usize> = line
            .split('\t')
            .map(|value| value.parse().unwrap())
            .collect();
        assert_eq!(fields.len(), 3);
        assert_eq!(
            ncbi_percent_match(fields[0], fields[1]),
            fields[2],
            "{line}"
        );
    }
}

#[cfg(test)]
#[path = "pairwise_tblastx_tests.rs"]
mod tblastx_tests;

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:418-424
/// The definition of a subject set in the pairwise report (`FillScanModeBlastDbInfo`), without
/// its period, which the writers add: an empty tag (an imported strategy's subjects) names no
/// input.
pub(crate) fn scan_mode_definition(tag: &str) -> String {
    if tag.is_empty() {
        "User specified sequence set".to_string()
    } else {
        format!("User specified sequence set (Input: {tag})")
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:418-424
/// Byte-preserving form of [`scan_mode_definition`] for subject tags that may not be UTF-8.
pub(crate) fn scan_mode_definition_bytes(tag: &[u8]) -> Vec<u8> {
    if tag.is_empty() {
        b"User specified sequence set".to_vec()
    } else {
        [
            b"User specified sequence set (Input: ".as_slice(),
            tag,
            b")",
        ]
        .concat()
    }
}
