//! Pairwise alignment output (outfmt 0)
//!
//! This module implements the traditional BLAST pairwise alignment output format.
//!
//! Reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp

use super::outfmt6::{format_bitscore_ncbi, format_evalue_ncbi, format_total_bitscore_ncbi};
use super::serial::asn_text::itoa_i32;
use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::common::Hit;
use crate::config::ScoringMatrix;
use crate::ncbi_str::truncate_spaces;
use crate::stats::{BlastGumbelBlk, KarlinParams};
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3563-3581
// ```c++
//     retval.Resize(k_NumAsciiChar, k_NumAsciiChar, -1000);
//
//     SNCBIFullScoreMatrix mtx;
//     NCBISM_Unpack(packed_mtx, &mtx);
//
//     for(int i = 0; i < ePMatrixSize; ++i){
//         for(int j = 0; j < ePMatrixSize; ++j){
//             retval((size_t)k_PSymbol[i], (size_t)k_PSymbol[j]) =
//                 mtx.s[(size_t)k_PSymbol[i]][(size_t)k_PSymbol[j]];
//         }
//     }
//     for(int i = 0; i < ePMatrixSize; ++i) {
//         retval((size_t)k_PSymbol[i], '*') = retval('*',(size_t)k_PSymbol[i]) = -4;
//     }
//     retval('*', '*') = 1;
//     // this is to count Selenocysteine to Cysteine matches as positive
//     retval('U', 'U') = retval('C', 'C');
//     retval('U', 'C') = retval('C', 'C');
//     retval('C', 'U') = retval('C', 'C');
// ```
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
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1410-1590
// (`PrintOneResultSet`: formats 0 to 4 write the same query header, description table and
// footer; `s_SetFlags`, :633-649, gives formats 1 to 4 `eMergeAlign` in place of outfmt 0's
// alignments, quoted at `report::query_anchored`.)
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
    /// prefetch is off (`report::db_subject::AfterThrow`): outfmt 0 skips the alignments of such a
    /// subject (or ends the run at its first alignment), 1 to 4 end the run at the first.
    pub fn after_throw(self) -> crate::report::db_subject::AfterThrow {
        use crate::report::db_subject::AfterThrow;
        match self {
            Self::Pairwise => AfterThrow::SkipsAlignment,
            Self::QueryAnchored { .. } => AfterThrow::EndsAtAlignment,
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
    // ```c++
    //     if(args.Exist(kArgSortHSPs) && args[kArgSortHSPs])
    //     {
    //         int hspsSortOption = args[kArgSortHSPs].AsInteger();
    //         if(m_OutputFormat == ePairwise) {
    //        	    m_HspsSortOption = hspsSortOption;
    //         }
    //         else {
    //             ERR_POST(Warning << "The parameter -sorthsps is ignored for output formats != 0.");
    // ```
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
/// them (`PairwiseView::reading`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PairwiseReading {
    pub sort: super::hit_sort::ReportSort,
    pub shown: crate::report::db_subject::Shown,
}

impl PairwiseConfig {
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1525-1528
    // ```c++
    //         aln_set = CBlastFormatUtil::SortSeqalignForSortableFormat(
    //                                       *(const_cast<CSeq_align_set*>(aln_set.GetPointer())),
    //                                       (m_Program == "tblastx") ? true : false,
    // ```
    /// Whether the report is TBLASTX's (`m_Program == "tblastx"`): the sorts and the description
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
// ```c++
//     if(m_HitsSortOption >= 0) {
//         flags |= CShowBlastDefline::eShowPercentIdent;
//         flags |= CShowBlastDefline::eShowTotalScore;
//         flags |= CShowBlastDefline::eShowQueryCoverage;
//     }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1527-1530
// ```c++
//     unique_ptr<CAlignFormatUtil::SSeqAlignSetCalcParams> seqSetInfo( CAlignFormatUtil::GetSeqAlignSetCalcParamsFromASN(aln));
//     if(seqSetInfo->hspNum == 0) {//calulated params are not in ASN - calculate now
//         seqSetInfo.reset( CAlignFormatUtil::GetSeqAlignSetCalcParams(aln,m_QueryLength,m_TranslatedNucAlignment));
//     }
// ```
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
    // NCBI c++/src/algo/blast/core/blast_kappa.c:331-342;
    // c++/src/objtools/align_format/showalign.cpp:3599-3604:
    // eCompoScaleOldMatrix => Method 1; other adjusted matrices => Method 2.
    pub comp_adjust_method: Option<u8>,
    // NCBI c++/src/objtools/align_format/showalign.cpp:3595-3598:
    // if (aln_vec_info->sum_n > 0) out << "(" << sum_n << ")";
    pub sum_n: Option<i32>,
    /// The columns (0-based) of `subject_seq` whose letter is a stored gap character (`-`,
    /// NCBIstdaa 0, which makeblastdb stores for a `-` of a protein FASTA): residues of the
    /// subject that the row shows as the letter of a gap. Empty when the subject holds none.
    /// NCBI counts the coordinates of the rows by the Dense-seg, where these columns are
    /// residues (`write_hsp_rows`).
    pub subject_stored_gaps: Vec<u32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1199-1203
    // ```c++
    //     // Set the identity score
    //     if (hsp->num_ident >= 0) {
    //         static const string kNumIdent("num_ident");
    //         scores.push_back(s_MakeScore(kNumIdent, 0.0, hsp->num_ident, true));
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:828-830
    // ```c++
    //     //look for scores at seqalign level first
    //     hasScore = s_GetBlastScore(aln.GetScore(), score, bits, evalue,
    //                                sum_n, num_ident, use_this_gi, comp_adj_method);
    // ```
    /// The Seq-align's `num_ident` score of the HSP: the engine's identities (`BlastHSP::num_ident`),
    /// the value that the program's Seq-align writer (outfmt 8, 11, 12) writes, which the report
    /// reads back with `GetAlnScores` for the sorts of `-sorthits` and `-sorthsps` and the
    /// description table's "Max ident" (`hit_sort::HspScores`). It can differ from `hit.num_ident`,
    /// the identities of the displayed rows: TBLASTN's engine counts a SEG-masked query letter
    /// as a mismatch; BLASTP's display recount (pyrrolysine, stored gap characters), BLASTX's and
    /// TBLASTX's rows translated by the display's tables count others. Every count is at least
    /// 0, so every Seq-align carries it.
    pub(crate) seqalign_num_ident: usize,
}

impl PairwiseHit {
    /// Whether the rows of the HSP are empty: BLASTX's ungapped HSPs of length 0
    /// (`algorithm::blastx::empty_rows::rows_empty`), the only ones without a column, whose
    /// `CSeq_align::GetSeqRange` throws `Row is empty` (Seq_align.cpp:227).
    pub(crate) fn rows_empty(&self) -> bool {
        self.hit.has_sort_offsets
            && super::serial::seqalign::ungapped_rows_empty(
                self.hit.sort_query_offset as i64,
                self.hit.sort_query_end as i64,
            )
    }
}

/// Where `write_blastx_pairwise_report` stopped a query's report at an HSP with empty rows
/// (`PairwiseHit::rows_empty`), as the error of its writes: which of NCBI's exceptions ends the
/// report there (`algorithm::blastx::empty_rows`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum EmptyRowsStop {
    /// `CSeq_align::GetSeqRange(): Row is empty` (Seq_align.cpp:227): the description table, or a
    /// sort's comparator, reads an empty row (after the query's `Length=`).
    RowIsEmpty,
    /// `CreateDensegFromStdseg(): Std-seg segment lengths not accurate!` (Seq_align.cpp:801): a
    /// sort's comparator or the alignment section makes a Dense-seg of `-subject`'s mapped
    /// Std-seg (before anything of the HSP), also the query-anchored section of outfmt 1 to 4
    /// (before its first byte).
    StdSegLengths,
    /// Not an exception: the merge of a query-anchored section (outfmt 1 to 4) frees a segment
    /// that it still reads (`blastx_anchored_stop`), which LOSAT refuses before the call's first
    /// byte (`algorithm::blastx::empty_rows::UNPORTED_ANCHORED_FREE`).
    MergeFreesSegment,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1522-1530
// ```c++
//     //invoke sorting only for m_HitsSortOption > CAlignFormatUtil::eEvalue or m_HspsSortOption > CAlignFormatUtil::eHspEvalue
//     if(m_HitsSortOption > 0 || m_HspsSortOption > 0) {
//         aln_set = CBlastFormatUtil::SortSeqalignForSortableFormat(
//                                       *(const_cast<CSeq_align_set*>(aln_set.GetPointer())),
//                                       (m_Program == "tblastx") ? true : false,
//                                       m_HitsSortOption,
//                                       m_HspsSortOption);
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3074-3094
// ```c++
//     if (hit_sort == eTotalScore) {
//         seqalign_hit_list.sort(SortHitByTotalScoreDescending);
//     } else if (hit_sort == eHighestScore) {
//         seqalign_hit_list.sort(CAlignFormatUtil::SortHitByScoreDescending);
//     } else if (hit_sort == ePercentIdentity) {
//         SortHitByPercentIdentityDescending(seqalign_hit_list,
//                                                nuc_to_nuc_translation);
//     } else if (hit_sort == eQueryCoverage) {
//         seqalign_hit_list.sort(SortHitByMasterCoverageDescending);
//     }
//
//     ITERATE(list< CRef<CSeq_align_set> >, iter2, seqalign_hit_list) {
//         CRef<CSeq_align_set> temp(*iter2);
//         if (hsp_sort == eQueryStart) {
//             temp->Set().sort(SortHspByMasterStartAscending);
//         } else if (hsp_sort == eHspPercentIdentity) {
//             temp->Set().sort(SortHspByPercentIdentityDescending);
//         } else if (hsp_sort == eScore) {
//             temp->Set().sort(SortHspByScoreDescending);
//         } else if (hsp_sort == eSubjectStart) {
//             temp->Set().sort(SortHspBySubjectStartAscending);
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4334-4335
// ```c++
//     ITERATE(CSeq_align_set::Tdata, iter, aln.Get()) {
//         int align_length = CAlignFormatUtil::GetAlignmentLength(**iter, do_translation);
// ```
/// The exception at which NCBI's sorts of outfmt 0 (`-sorthits`, `-sorthsps` above 0, before the
/// description table) stop a query's report whose hit list (`order`, the subjects' hits in
/// `grouped`) holds an HSP with empty rows; `remaps`: `-subject`, whose Seq-loc mapper gives that
/// HSP a query row of one nucleotide in a subject with rows of letters (`mapper_rewrites`; the
/// rows of a subject of empty rows only stay as `-db`'s). `std::list::sort` compares every
/// element of a list of two or more, so a comparator that reads an HSP's ranges reaches it:
/// - the subjects by percent identity (`-sorthits 3`: `GetSeqAlignSetCalcPercentIdent` makes a
///   Dense-seg of every HSP, `GetAlignmentLength`): `-subject`'s mapped Std-seg throws (801), the
///   raw empty rows of `-db` give a length of 0;
/// - the subjects by query coverage (`-sorthits 4`: `GetMasterCoverage`): `-db`'s empty query row
///   (227), `-subject`'s nucleotide is not empty;
/// - then each subject's HSPs, a list of two or more: by query start (`-sorthsps 2`,
///   `GetSeqStart(0)`) the empty query row of `-db` (227); by percent identity (`-sorthsps 3`)
///   `-subject`'s Dense-seg (801); by subject start (`-sorthsps 4`, `GetSeqStart(1)`) the empty
///   subject row of both (227).
///
/// The scores' sorts (`-sorthits 1` and `2`, `-sorthsps 1`) read no range. Probes: survey set K
/// (`X/SURVEY_X2g_formats.md` 6).
fn sort_stop(
    sort: super::hit_sort::ReportSort,
    remaps: bool,
    order: &[u32],
    grouped: &std::collections::HashMap<u32, Vec<&PairwiseHit>>,
) -> Option<EmptyRowsStop> {
    // Whether the subject holds an HSP with empty rows that the mapper rewrote (`Some(true)`),
    // left as made (`Some(false)`), or either (`None`).
    let holds = |oid: &u32, mapped: Option<bool>| {
        let hits = &grouped[oid];
        hits.iter().any(|hit| hit.rows_empty())
            && mapped.is_none_or(|mapped| mapper_rewrites(remaps, hits) == mapped)
    };
    let hit_stop = match sort.hits {
        Some(3) => Some((EmptyRowsStop::StdSegLengths, Some(true))),
        Some(4) => Some((EmptyRowsStop::RowIsEmpty, Some(false))),
        _ => None,
    };
    if let Some((stop, mapped)) = hit_stop {
        if order.len() >= 2 && order.iter().any(|oid| holds(oid, mapped)) {
            return Some(stop);
        }
    }
    let hsp_stop = match sort.hsps {
        Some(2) => Some((EmptyRowsStop::RowIsEmpty, Some(false))),
        Some(3) => Some((EmptyRowsStop::StdSegLengths, Some(true))),
        Some(4) => Some((EmptyRowsStop::RowIsEmpty, None)),
        _ => None,
    };
    if let Some((stop, mapped)) = hsp_stop {
        if order
            .iter()
            .any(|oid| grouped[oid].len() >= 2 && holds(oid, mapped))
        {
            return Some(stop);
        }
    }
    None
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
            seqalign_num_ident,
        }
    }
}

/// NCBI BLASTP pairwise per-query ancillary data.
///
/// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/blast_results.cpp:72-115
/// ```c
/// CBlastAncillaryData::CBlastAncillaryData(EBlastProgramType program_type,
///                     int query_number,
///                     const BlastScoreBlk *sbp,
///                     const BlastQueryInfo *query_info)
/// {
///     ...
///     m_SearchSpace = ctx->eff_searchsp;
///     ...
///     s_InitializeKarlinBlk(sbp->kbp_std[ctx_index], &m_UngappedKarlinBlk);
/// }
/// ```
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
/// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:2248-2285
/// ```c
/// if ( !m_IsBl2Seq || m_IsDbScan) {
///     CBlastFormatUtil::PrintDbReport(m_DbInfo, kFormatLineLength,
///                                     m_Outfile, false);
/// }
/// ...
/// m_Outfile << "\n\nMatrix: " << options.GetMatrixName() << "\n";
/// ...
/// m_Outfile << "Neighboring words threshold: "
///           << options.GetWordThreshold() << "\n";
/// m_Outfile << "Window for multiple hits: "
///           << options.GetWindowSize() << "\n";
/// ```
#[derive(Debug, Clone)]
pub struct BlastpPairwiseReport {
    pub version: String,
    pub database_name: String,
    pub database_num_sequences: usize,
    pub database_total_letters: usize,
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
    /// Subjects in the description table and with alignments, per query (NCBI's
    /// `m_NumDescriptions` and `m_NumAlignments`).
    pub num_descriptions: usize,
    pub num_alignments: usize,
}

/// One database of a `-db` list as the reports show it (NCBI's `SDbInfo` of a BLAST
/// database, `seqdb::NameRecord`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:385-391
/// ```c++
///     info.name = seqdb->GetDBNameList();
///     info.definition = seqdb->GetTitle();
///     if (info.definition.empty())
///         info.definition = info.name;
///     info.date = seqdb->GetDate();
///     info.total_length = seqdb->GetTotalLength();
///     info.number_seqs = seqdb->GetNumSeqs();
/// ```
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
// ```c++
//     if (top) {
//         const CAlignFormatUtil::SDbInfo* dbinfo = &(dbinfo_list.front());
//         out << "Database: ";
//
//         string db_titles = dbinfo->definition;
//         Int8 tot_num_seqs = static_cast<Int8>(dbinfo->number_seqs);
//         Int8 tot_length = dbinfo->total_length;
//
//         for (size_t i = 1; i < dbinfo_list.size(); i++) {
//             db_titles += "; " + dbinfo_list[i].definition;
//             tot_num_seqs += static_cast<Int8>(dbinfo_list[i].number_seqs);
//             tot_length += dbinfo_list[i].total_length;
//         }
// ```
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
/// Reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp
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
    // NCBI reference: c++/src/objtools/align_format/showalign.cpp:2385-2388,2461-2471,340-359
    // out << ">"; if (out.tellp() > 1L) out << " "; s_WrapOutputLine(out, alnDispParams->title);
    // `s_WrapOutputLine` counts and writes bytes (`str[i]`, `isspace((unsigned char)
    // str[i])`), so a line can end inside a UTF-8 sequence.
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
/// and what follows it). `probe` marks the heading and each HSP without changing the
/// written bytes.
///
/// NCBI reference: c++/src/objtools/align_format/showalign.cpp:3613-3632
/// ```c++
///     if(show_defline) {
/// 		...
/// 				string deflines = x_PrintDefLine(bsp_handle, aln_vec_info);
/// 				out<< deflines;
/// 		...
/// 			out << "\n";
/// ```
/// NCBI reference: c++/src/objtools/align_format/showalign.cpp:2273
/// ```c++
///             alnDispParams->title = CDeflineGenerator().GenerateDefline(bsp_handle);
/// ```
/// Where `GenerateDefline` throws (a title whose encoding NCBI cannot guess,
/// `UnknownEncoding`), NCBI catches the exception for that HSP and writes one line in its
/// place. `previousId` is not set then, so the next HSP of the subject tries the heading
/// again: every HSP of the subject gets the line. `GetSeqIdString` of the local ID
/// `lcl|Subject_N` is `Subject_N`, N being the subject's position in the input, which is
/// `s_idx + 1` (every subject record keeps its index). The probe gets an empty heading.
///
/// NCBI reference: c++/src/objtools/align_format/showalign.cpp:1950-1981
/// ```c++
///                 try{
///                     const CBioseq_Handle& handle = avRef->GetBioseqHandle(1);
///                     if(handle){
///                         ...
///                         subid=&(avRef->GetSeqId(1));
///                         bool showDefLine = previousId.Empty() || !subid->Match(*previousId);
///                         x_DisplayAlnvecInfo(out, alnvecInfo,showDefLine);
///
///                         previousId = subid;
///                     }
///                 } catch (const CException&){
///                     out << "Sequence with id "
///                         << (avRef->GetSeqId(1)).GetSeqIdString().c_str()
///                         <<" no longer exists in database...alignment skipped\n";
///                     continue;
///                 }
/// ```
/// The alignments of the HSPs `part` of the subject's (the heading with the first); the parts of
/// `0..hits.len()` in order write its bytes (`x_pairwise_par`'s pieces).
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
    mut write_hsp: impl FnMut(&mut W, &PairwiseHit) -> io::Result<()>,
) -> io::Result<()> {
    let first_index = hsp_indices[0];
    let heading = super::defline::generate_defline(title, protein, false);
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
    for (hit, &hsp_index) in hits[part.clone()].iter().zip(&hsp_indices[part]) {
        if let Some(probe) = probe.as_mut() {
            writer.flush()?;
            probe.begin(hsp_index);
        }
        match heading {
            Ok(_) => write_hsp(writer, hit)?,
            Err(super::defline::UnknownEncoding) => {
                write_skipped_alignment(writer, &record_skipped_id(parsed, s_idx))?
            }
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
/// ```c++
/// 		if(bdl->IsSetTitle()){
/// 			alnDispParams->title = bdl->GetTitle();
/// 		}
///         if(alnDispParams->title.empty()) {
///             alnDispParams->title = CDeflineGenerator().GenerateDefline(bsp_handle);
///         }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2383-2386
/// ```c++
///                 if(alnDispParams) {
///                     numBdl++;
///                     if(isFirst){
///                         out << ">";
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2421-2422
/// ```c++
///                     if(!(alnDispParams->seqID->AsFastaString().find("gnl|BL_ORD_ID") != string::npos) ||
/// 			alnDispParams->seqID->AsFastaString().find("lcl|Subject_") != string::npos){
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2461-2479
/// ```c++
///                     if (out.tellp() > 1L) {
///                         out << " ";
///                     }
/// 					if(!alnDispParams->title.empty()) {
///                         s_WrapOutputLine(out, (m_AlignOption&eHtml) ?
///                                          CHTMLHelper::
///                                          HTMLEncode(alnDispParams->title) :
///                                          alnDispParams->title);
///                     }
///                     out<<"\n";
///                     isFirst = false;
///                 }
///             }
///             ...
///     out<<kLengthString<<bsp_handle.GetBioseqLength()<<"\n";
/// ```
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
    hits: &[&PairwiseHit],
    hsp_indices: &[HspIndex],
    part: std::ops::Range<usize>,
    mut write_hsp: impl FnMut(&mut W, &PairwiseHit) -> io::Result<()>,
) -> io::Result<()> {
    let first_index = hsp_indices[0];
    let skipped = if part.start == 0 {
        if let Some(probe) = probe.as_mut() {
            writer.flush()?;
            probe.subject_begin(first_index);
        }
        let skipped =
            write_database_subject_heading(writer, title, protein, hits[0].subject_length)?;
        if let Some(probe) = probe.as_mut() {
            writer.flush()?;
            probe.subject_end(first_index);
        }
        skipped
    } else {
        title.skipped_id(protein)
    };
    for (hit, &hsp_index) in hits[part.clone()].iter().zip(&hsp_indices[part]) {
        if let Some(probe) = probe.as_mut() {
            writer.flush()?;
            probe.begin(hsp_index);
        }
        match &skipped {
            None => write_hsp(writer, hit)?,
            Some(id) => write_skipped_alignment(writer, id)?,
        }
        if let Some(probe) = probe.as_mut() {
            writer.flush()?;
            probe.end(hsp_index);
        }
    }
    Ok(())
}

/// NCBI's line in place of an alignment whose heading it cannot make (showalign.cpp:1976-1981,
/// quoted at `write_subject_alignments_part`): of a database subject (`DatabaseHeading::skipped_id`)
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
    /// A sequence whose lookup throws where the alignment loads it (its batch's prefetch is off,
    /// `report::db_subject::HitLoad::Throws`): the `catch` of each alignment writes NCBI's line in
    /// its place (`skipped_id`, showalign.cpp:1976-1981), and no heading.
    Thrown(&'a crate::report::db_subject::ParsedSubject),
}

impl<'a> DatabaseHeading<'a> {
    /// The Seq-id text that outfmt 0 writes in place of the subject's alignments (`Sequence with
    /// id <id> no longer exists in database...alignment skipped`, showalign.cpp:1976-1981, quoted
    /// at `ParsedSubject::heading_skipped_id`) where the heading's generated title throws: an empty
    /// stored title made from a name that NCBI cannot decode (`SourceTaxname::undecodable_id`), or
    /// the generated title of a sequence with parsed Seq-ids (`heading_skipped_id`).
    fn skipped_id(self, protein: bool) -> Option<Vec<u8>> {
        match self {
            DatabaseHeading::Title([], source) => source.undecodable_id().map(<[u8]>::to_vec),
            DatabaseHeading::Title(..) => None,
            DatabaseHeading::Parsed(subject, every_line) => {
                subject.heading_skipped_id(protein, every_line)
            }
            DatabaseHeading::Thrown(subject) => Some(subject.skipped_id()),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:340-359
    // ```c++
    // static void s_WrapOutputLine(CNcbiOstream& out, const string& str)
    // {
    //     const int line_len = 60;
    //     bool do_wrap = false;
    //     int length = (int) str.size();
    //     if (length > line_len) {
    //         for (int i = 0; i < length; i ++){
    //             if(i > 0 && i % line_len == 0){
    //                 do_wrap = true;
    //             }
    //             out << str[i];
    //             if(do_wrap && isspace((unsigned char) str[i])){
    //                 out << "\n";
    //                 do_wrap = false;
    //             }
    //         }
    //     } else {
    //         out << str;
    //     }
    // }
    // ```
    /// Writes the definition-line lines of the heading (each ending with a newline; `Length=`
    /// follows): `>` and the stored title wrapped as `s_WrapOutputLine` wraps it, an empty one
    /// made by `GenerateDefline` from the source (`showalign.cpp:2250-2255`), or of a sequence with
    /// parsed Seq-ids the line of each definition line (`ParsedSubject::heading`). Where the title
    /// made for the heading throws it writes nothing and gives the Seq-id text of NCBI's line in
    /// place of each alignment (`skipped_id`). The title is made once.
    fn write_lines<W: Write>(self, writer: &mut W, protein: bool) -> io::Result<Option<Vec<u8>>> {
        let (title, source) = match self {
            DatabaseHeading::Title(title, source) => (title, source),
            DatabaseHeading::Parsed(subject, every_line) => {
                return match subject.heading(protein, every_line) {
                    Ok(lines) => writer.write_all(&lines).map(|()| None),
                    Err(id) => Ok(Some(id)),
                };
            }
            DatabaseHeading::Thrown(subject) => return Ok(Some(subject.skipped_id())),
        };
        let generated;
        let heading: &[u8] = if title.is_empty() {
            if let Some(id) = source.undecodable_id() {
                return Ok(Some(id.to_vec()));
            }
            // NCBI makes an empty title from the Bioseq's source (`showalign.cpp:2250-2255`,
            // quoted at `ParsedSubject::generated_title_fetched`; LOSAT refuses a source that it
            // refuses before the report, `StoredTitle::check_shown_source`).
            let facts = super::defline::BioseqFacts::of_source(source.name().ok().flatten());
            generated = super::defline::generate_bioseq_defline(title, protein, false, &facts)
                .map(|defline| defline.text)
                .unwrap_or_default();
            &generated
        } else {
            title
        };
        writer.write_all(b">")?;
        write_wrapped(writer, heading)?;
        writer.write_all(b"\n")?;
        Ok(None)
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

    /// The heading of a subject with parsed Seq-ids as `parsed` makes it, or `Thrown` where its
    /// load throws (`thrown`, `x_pairwise_par::QuerySections::thrown`).
    fn loaded(
        subject: &'a crate::report::db_subject::ParsedSubject,
        gapped: bool,
        hsps: usize,
        thrown: bool,
    ) -> Self {
        if thrown {
            return DatabaseHeading::Thrown(subject);
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
    let local = parsed
        .get(s_idx as usize)
        .and_then(Option::as_deref)
        .is_some_and(crate::report::db_subject::ParsedSubject::alignment_id_is_local);
    Some((exception, local))
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
) -> io::Result<Option<Vec<u8>>> {
    if let Some(id) = heading.write_lines(writer, protein)? {
        return Ok(Some(id));
    }
    if let Some(length) = length {
        writeln!(writer, "Length={length}")?;
    }
    writeln!(writer).map(|()| None)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
// ```c++
// int CAlignFormatUtil::GetPercentMatch(int numerator, int denominator)
// {
//      if (numerator == denominator)
//         return 100;
//      else {
//        int retval =(int) (0.5 + 100.0*((double)numerator)/((double)denominator));
//        retval = min(99, retval);
//        return retval;
//      }
// }
// ```
fn ncbi_percent_match(numerator: usize, denominator: usize) -> usize {
    if numerator == denominator {
        100
    } else {
        ((0.5 + 100.0 * numerator as f64 / denominator as f64) as usize).min(99)
    }
}

/// The `Expect` label of an HSP's header: `Expect(n)` for an HSP of a linked set of n > 1 HSPs
/// (`sum_n`, `BlastHSP::num`), else `Expect` (showalign.cpp:3593-3598, quoted at
/// `write_hsp_info`).
#[inline]
fn expect_label(sum_n: Option<i32>) -> String {
    sum_n
        .filter(|&count| count > 1)
        .map_or_else(|| "Expect".to_string(), |count| format!("Expect({count})"))
}

/// Write HSP score information
///
/// Reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp x_DisplayAlignInfo()
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
    // ```c++
    //             out<<" Score = "<<bit_score_buf<<" ";
    //             out<<"bits ("<<aln_vec_info->score<<"),"<<"  ";
    //             out<<"Expect";
    //             if (aln_vec_info->sum_n > 0) {
    //             out << "(" << aln_vec_info->sum_n << ")";
    //             }
    // ```
    // An HSP of a linked set (psiblast's `-sum_stats true`) has its "sum_n" score, the set's size.
    let expect = expect_label(hit.sum_n);
    // psiblast's first round has blastp's methods; rpsblast's redo marks its HSPs likewise.
    if matches!(config.program.as_str(), "blastp" | "psiblast" | "rpsblast") {
        // NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:3600-3603
        // ```c
        //             if (aln_vec_info->comp_adj_method == 1)
        //             out << ", Method: Composition-based stats.";
        //             else if (aln_vec_info->comp_adj_method == 2)
        //             out << ", Method: Compositional matrix adjust.";
        // ```
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

    // Identity/Positives/Gaps line
    let align_len = h.length;
    let num_ident = h.num_ident;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
    // ```c++
    // int CAlignFormatUtil::GetPercentMatch(int numerator, int denominator)
    // {
    //      if (numerator == denominator)
    //         return 100;
    //      else {
    //        int retval =(int) (0.5 + 100.0*((double)numerator)/((double)denominator));
    //        retval = min(99, retval);
    //        return retval;
    //      }
    // }
    // ```
    let ident_pct = ncbi_percent_match(num_ident, align_len);

    let positives = hit.positives.unwrap_or(h.num_positives.max(num_ident));
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
    // ```c++
    // int CAlignFormatUtil::GetPercentMatch(int numerator, int denominator)
    // {
    //      if (numerator == denominator)
    //         return 100;
    //      else {
    //        int retval =(int) (0.5 + 100.0*((double)numerator)/((double)denominator));
    //        retval = min(99, retval);
    //        return retval;
    //      }
    // }
    // ```
    let pos_pct = ncbi_percent_match(positives, align_len);

    let gaps = hit.gaps.unwrap_or_else(|| h.gap_letters());
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
    // ```c++
    // int CAlignFormatUtil::GetPercentMatch(int numerator, int denominator)
    // {
    //      if (numerator == denominator)
    //         return 100;
    //      else {
    //        int retval =(int) (0.5 + 100.0*((double)numerator)/((double)denominator));
    //        retval = min(99, retval);
    //        return retval;
    //      }
    // }
    // ```
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
        // NCBI reference: c++/src/objtools/align_format/showalign.cpp:1794-1805,2131-2141
        // ```c++
        //     x_FillIdentityInfo(aln_vec_info->alnRowInfo->sequence[0],
        //                        aln_vec_info->alnRowInfo->sequence[1],
        //                        aln_vec_info->match,
        //     ...
        //         aln_vec_info->gap = x_GetNumGaps();
        // ...
        //         if(sequence_standard[i]==sequence[i]){
        //             ...
        //             match ++;
        // ```
        // The displayed counts come from the two rows (before the lowercase masking of
        // the display), and the gaps are the gap columns of both rows.
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
        let gaps = query_row
            .iter()
            .chain(subject_row)
            .filter(|&&residue| residue == b'-')
            .count();
        // NCBI reference: c++/src/objtools/align_format/showalign.cpp:310-320
        // ```c++
        //     out<<" Identities = "<<match<<"/"<<(aln_stop+1)<<" ("<<identity<<"%"<<")";
        //     ...
        //     out<<", Gaps = "<<gap<<"/"<<(aln_stop+1)
        //        <<" ("<<CAlignFormatUtil::GetPercentMatch(gap, aln_stop+1)<<"%"<<")"<<"\n";
        //     if (!aln_is_prot){
        //         out<<" Strand="<<(master_strand==1 ? "Plus" : "Minus")
        //            <<"/"<<(slave_strand==1? "Plus" : "Minus")<<"\n";
        // ```
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
        // NCBI reference: c++/src/objtools/align_format/showalign.cpp:4011-4018
        // ```c++
        //         s_DisplayIdentityInfo(out,
        //                               ...
        //                               m_AV->StrandSign(0),
        //                               m_AV->StrandSign(1),
        // ```
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
/// Reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp x_DisplayRowData()
///
/// Format:
/// ```text
/// Query  1    MVLSPADKTN VKAAWGKVGA HAGEYGAEAL ERMFLSFPTT KTYFPHFDLS H  60
///             MVLSPADKTN VKAAWGKVGA HAGEYGAEAL ERMFLSFPTT KTYFPHFDLS H
/// Sbjct  1    MVLSPADKTN VKAAWGKVGA HAGEYGAEAL ERMFLSFPTT KTYFPHFDLS H  60
/// ```
pub fn write_alignment<W: Write>(
    writer: &mut W,
    hit: &PairwiseHit,
    config: &PairwiseConfig,
) -> io::Result<()> {
    let h = &hit.hit;

    // If we have actual sequences, display them
    if let (Some(ref qseq), Some(ref sseq)) = (&hit.query_seq, &hit.subject_seq) {
        write_alignment_with_sequences(writer, h, qseq, sseq, &hit.subject_stored_gaps, config)?;
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
/// from high to low), rows of letters (`write_hsp_rows`).
fn write_alignment_with_sequences<W: Write>(
    writer: &mut W,
    hit: &Hit,
    query_seq: &str,
    subject_seq: &str,
    subject_stored_gaps: &[u32],
    config: &PairwiseConfig,
) -> io::Result<()> {
    let line_len = config.line_length;
    let q_chars: Vec<char> = query_seq.chars().collect();
    let s_chars: Vec<char> = subject_seq.chars().collect();

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/showalign.cpp:96-98
    // ```c
    // static const int k_IdStartMargin = 2;
    // static const int k_SeqStopMargin = 2;
    // static const int k_StartSequenceMargin = 2;
    // ```
    //
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/showalign.cpp:1600-1625
    // ```c
    // CAlignFormatUtil::AddSpace(out, alnRoInfo->maxIdLen-alnRoInfo->seqidArray[row].size()
    //                            + k_IdStartMargin);
    // out << start;
    // CAlignFormatUtil::AddSpace(out, alnRoInfo->maxStartLen-startLen + k_StartSequenceMargin);
    // ...
    // CAlignFormatUtil::AddSpace(out, k_SeqStopMargin);
    // out << end;
    // ```
    let max_start_len = coordinate_width(hit);
    // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1570-1573
    // ```c++
    //     if (m_Program == "blastn" || m_Program == "megablast") {
    //         display.SetMiddleLineStyle(CDisplaySeqalign::eBar);
    //         display.SetAlignType(CDisplaySeqalign::eNuc);
    // ```
    let nucleotide = config.program == "blastn";
    // The rows' end gaps (`end_gap_span`; the subject's gaps stored as residues are residues).
    let query_span = end_gap_span(&q_chars, |_| false);
    let subject_span = end_gap_span(&s_chars, |column| {
        subject_stored_gaps.binary_search(&(column as u32)).is_ok()
    });

    // NCBI reference: c++/src/objtools/align_format/showalign.cpp:2131-2155
    // ```c++
    //         if(sequence_standard[i]==sequence[i]){
    //             if(m_AlignOption & eShowMiddleLine){
    //                 if(m_MidLineStyle == eBar ) {
    //                     middle_line[i] = '|';
    //                 } else if (m_MidLineStyle == eChar){
    //                     middle_line[i] = sequence[i];
    //         ...
    //                 if (m_AlignOption & eShowMiddleLine){
    //                     middle_line[i] = ' ';
    // ```
    // The rows may carry lowercase masking, which NCBI applies only when it prints a
    // row (showalign.cpp:2495-2521), so they are compared without case, and the protein
    // middle line shows the uppercase residue.
    let write_middle = |writer: &mut W, columns: std::ops::Range<usize>| {
        let middle = q_chars[columns.clone()]
            .iter()
            .zip(s_chars[columns].iter())
            .map(|(q, s)| {
                let (q, s) = (q.to_ascii_uppercase(), s.to_ascii_uppercase());
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

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1497-1502,1535-1538
    // ```c++
    //     for (int row=0; row<rowNum; row++) {
    //         bool hasSequence = true;
    //         if (!(m_AlignOption & eShowGapOnlyLines)) {
    //             hasSequence = curRange.IntersectingWith(alnRoInfo->rowRng[row]);
    //         }
    //         //only output rows that have sequence
    //         if (hasSequence){
    //     ...
    //             //display middle line for pairwise
    //             if (row == 0 && ((m_AlignOption & eShowMiddleLine)) && !(m_AlignOption&eMergeAlign)) {
    //                 x_DisplayMiddLine(alnRoInfo, row,out);
    //             }
    // ```
    // A row is written on a line that meets its columns between its end gaps
    // (`rowRng = m_AV->GetSeqAlnRange(row)`, showalign.cpp:1272; `write_hsp_rows`), the middle
    // line with the query's; its end gaps are NCBI's end character (`write_row_chars`).
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
        q_chars.len(),
        max_start_len,
        |writer, row, columns| {
            if row == 0 {
                write_row_chars(writer, &q_chars, columns, query_span)
            } else {
                write_row_chars(writer, &s_chars, columns, subject_span)
            }
        },
        write_middle,
    )
}

thread_local! {
    /// The vectors of outfmt 0's rows (`write_hsp_rows`): an HSP's segments and its `CAlnVec`,
    /// kept from HSP to HSP of the thread's reports. No NCBI counterpart: NCBI makes a new
    /// Seq-align and `CAlnVec` per HSP (showalign.cpp:1823-1861).
    static HSP_ROWS: std::cell::RefCell<(
        super::serial::seqalign::Segments,
        super::query_anchored::HspAlnVec,
    )> = std::cell::RefCell::default();
}

/// The labels of an HSP's rows in outfmt 0 (`eShowBlastStyleId`, `x_FillSeqid`).
const ROW_LABELS: [&str; 2] = ["Query", "Sbjct"];

/// Writes an HSP's rows in outfmt 0, `line_length` columns a line: each line's query row, the
/// middle line, the subject row and a blank line. The rows' coordinates come from the HSP's
/// `CAlnVec` (`query_anchored::HspAlnVec`) of its Seq-align's segments, which `segments_of` makes
/// (the producer of the program, `report::serial::seqalign`, as the query-anchored reports'
/// `section_hsps`), of row widths `widths` (`HspStrands::no_widths`: none); `columns` is the
/// number of columns of the rows' letters, `write_letters(writer, row, columns)` writes a row's
/// letters of a line and `write_middle(writer, columns)` the middle line's. `width` is the width
/// of the coordinates (`maxStartLen`, `coordinate_width`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1435-1438,1452-1453
/// ```c++
/// void CDisplaySeqalign::x_DisplayRowData(SAlnRowInfo *alnRoInfo,CNcbiOstream& out)
/// {
///     size_t aln_stop=m_AV->GetAlnStop();
///     int rowNum = alnRoInfo->rowNum;
///     ...
///     for(int j=0; j<=(int)aln_stop; j+=(int)m_LineLen){
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1486-1497
/// ```c++
///     //output according to aln coordinates
///     if(aln_stop-aln_start+1<m_LineLen) {
///         actualLineLen=aln_stop-aln_start+1;
///     } else {
///         actualLineLen=m_LineLen;
///     }
///     CAlnMap::TSignedRange curRange(aln_start, aln_start+(int)actualLineLen-1);
///     alnRoInfo->currPrintSegment = aln_start;
///     alnRoInfo->currActualLineLen = actualLineLen;
///     alnRoInfo->currRange = curRange;
///     //here is each row
///     for (int row=0; row<rowNum; row++) {
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1538-1543,1550-1551
/// ```c++
///             //display middle line for pairwise
///             if (row == 0 && ((m_AlignOption & eShowMiddleLine)) && !(m_AlignOption&eMergeAlign)) {
///                 x_DisplayMiddLine(alnRoInfo, row,out);
///             }
///             prev_stop[row] = end;
///         }
///     ...
///     }//end of displaying rows
///     out<<"\n";
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1781-1783
/// ```c++
///     CAlignFormatUtil:: AddSpace(out, alnRoInfo->maxIdLen + k_IdStartMargin + alnRoInfo->maxStartLen + k_StartSequenceMargin);
///     x_OutputSeq(alnRoInfo->middleLine, no_id, j, (int)actualLineLen, 0, row, false, alnRoInfo->masked_regions[row], out);
///     out<<"\n";
/// ```
/// The lines are those of the HSP's Dense-seg (`GetAlnStop`), whose columns the rows' letters
/// are: a difference is an internal error. An HSP of length 0 (BLASTX's empty rows) has none
/// (`(int)aln_stop` is -1). `prev_stop` starts at 0 (`vector<int>`).
#[allow(clippy::too_many_arguments)]
fn write_hsp_rows<W: Write>(
    writer: &mut W,
    segments_of: impl FnOnce(
        &mut super::serial::seqalign::Segments,
    ) -> super::query_anchored::HspStrands,
    widths: [u32; 2],
    line_length: usize,
    columns: usize,
    width: usize,
    mut write_letters: impl FnMut(&mut W, usize, std::ops::Range<usize>) -> io::Result<()>,
    mut write_middle: impl FnMut(&mut W, std::ops::Range<usize>) -> io::Result<()>,
) -> io::Result<()> {
    HSP_ROWS.with_borrow_mut(|(segments, rows)| {
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
                    |writer| write_letters(writer, row, line.clone()),
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
/// ```c++
///     size_t startLen=0;
///     int start = alnRoInfo->seqStarts[row].front() + 1;  //+1 for 1 based
///     int end = alnRoInfo->seqStops[row].front() + 1;
///     int j = alnRoInfo->currPrintSegment;
///     int actualLineLen = alnRoInfo->currActualLineLen;
///     //print out sequence line
///     //adjust space between id and start
///     CAlignFormatUtil::AddSpace(out, alnRoInfo->maxIdLen-alnRoInfo->seqidArray[row].size() + k_IdStartMargin);
///     //not to display start and stop number for empty row
///     if ((j > 0 && end == prev_stop)
///         || (j == 0 && start == 1 && end == 1)) {
///         startLen = 0;
///     } else {
///         out << start;
///         startLen=NStr::IntToString(start).size();
///     }
///
///     CAlignFormatUtil::AddSpace(out, alnRoInfo->maxStartLen-startLen + k_StartSequenceMargin);
///     x_OutputSeq(alnRoInfo->sequence[row], m_AV->GetSeqId(row), j,
///                 (int)actualLineLen, alnRoInfo->frame[row], row,
///                 (row > 0 && alnRoInfo->colorMismatch)?true:false,
///                 alnRoInfo->masked_regions[row], out);
///     CAlignFormatUtil::AddSpace(out, k_SeqStopMargin);
///
///      //not to display stop number for empty row in the middle
///     if (!(j > 0 && end == prev_stop)
///         && !(j == 0 && start == 1 && end == 1)) {
///         out << end;
///     }
///     out<<"\n";
/// ```
/// A later line without residues ends where the row's previous line ended (`end ==
/// prev_stop`). `AddSpace` takes a `size_t`, so an omitted start gives `maxStartLen + 2`
/// spaces. The labels of the pairwise rows are equally long (`maxIdLen`): two spaces follow. A
/// first line of one residue at 1 of a row of letters shows no coordinates (the `-line_length 1`
/// rows of an HSP from a sequence's first letter); a translated row's first residue spans a
/// codon, so its start and end differ.
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
    // `startLen`: the length of the start's text (`NStr::IntToString(start).size()`), 0 without.
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
// ```c++
//     _ASSERT((int)sequence.size() > start);
//     list<CRange<int> > actualSeqloc;
//     string actualSeq = sequence.substr(start, len);
// ```
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

/// Writes `letters` (a chunk of a row kept as `char`s, `write_row_chars`, or a middle line,
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
// ```c++
//         if(sequence_standard[i]==sequence[i]){
//             if(m_AlignOption & eShowMiddleLine){
//                 if(m_MidLineStyle == eBar ) {
//                     middle_line[i] = '|';
//                 } else if (m_MidLineStyle == eChar){
//                     middle_line[i] = sequence[i];
//                 }
//             }
//             match ++;
//         } else {
//             if ((m_AlignType&eProt)
//                 && m_Matrix[(int)sequence_standard[i]][(int)sequence[i]] > 0){
//                 positive ++;
//                 if(m_AlignOption & eShowMiddleLine){
//                     if (m_MidLineStyle == eChar){
//                         middle_line[i] = '+';
//                     }
//                 }
//             } else {
//                 if (m_AlignOption & eShowMiddleLine){
//                     middle_line[i] = ' ';
//                 }
//             }
//         }
// ```
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
) -> io::Result<()> {
    let middle = query[columns.clone()]
        .iter()
        .zip(&subject[columns])
        .map(|(&q, &s)| {
            let q = q.to_ascii_uppercase();
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
// ```c++
// static void s_WrapOutputLine(CNcbiOstream& out, const string& str)
// {
//     const int line_len = 60;
//     bool do_wrap = false;
//     int length = (int) str.size();
//     if (length > line_len) {
//         for (int i = 0; i < length; i ++){
//             if(i > 0 && i % line_len == 0){
//                 do_wrap = true;
//             }
//             out << str[i];
//             if(do_wrap && isspace((unsigned char) str[i])){
//                 out << "\n";
//                 do_wrap = false;
//             }
//         }
//     } else {
//         out << str;
//     }
// }
// ```
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

/// Writes the columns `columns` of the row `row` kept as `char`s, its end gaps (the columns
/// outside `span`, `end_gap_span`) as NCBI's end character, a space.
fn write_row_chars<W: Write>(
    writer: &mut W,
    row: &[char],
    columns: std::ops::Range<usize>,
    span: (usize, usize),
) -> io::Result<()> {
    if span.0 <= columns.start && columns.end <= span.1 {
        return write_buffered(writer, row[columns].iter().copied());
    }
    write_buffered(
        writer,
        columns.map(|column| {
            if column < span.0 || column >= span.1 {
                ' '
            } else {
                row[column]
            }
        }),
    )
}

// NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnvec.hpp:386-393
// ```c++
// inline
// CSeqVector::TResidue CAlnVec::GetEndChar() const
// {
//     if (IsSetEndChar()) {
//         return m_EndChar;
//     } else {
//         return ' ';
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:362-370 (`GetWholeAlnSeqString`)
// ```c++
//                 if (seg < left_seg  ||  seg > right_seg) {
//                     fill_ch = GetEndChar();
//                 } else {
//                     fill_ch = GetGapChar(row);
//                 }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1209-1216
// ```c++
//     m_AV->SetGapChar('-');
//
//     if (m_AlignOption & eShowEndGaps) {
//         m_AV->SetEndChar('-');
//     }
//     else {
//         m_AV->SetEndChar(' ');
//     }
// ```
// (`eShowEndGaps` is set by the query-anchored views only, blast_format.cpp:1754-1758; the
// tabular writer sets the gap character alone, tabular.cpp:982-986.)
/// The columns `lo..hi` of a row of a pairwise alignment between its end gaps: its columns before
/// its first residue and after its last, which NCBI's alignment strings fill with the end
/// character (a space in the pairwise view and in the tabular `qseq`, `sseq` and `btop`) where an
/// inner gap is `-`. `residue_dash(column)` says whether a `-` of the row is a residue (a gap
/// stored in the sequence). Only a Smith-Waterman alignment has end gaps (rpsblast's
/// `-use_sw_tback`: `ALIGN_EX` runs from the Smith-Waterman start over a rectangle whose best path
/// may start with a gap); an X-drop alignment starts and ends with a substitution. The `-` that
/// fills a segment past a sequence's end is the gap character (alnvec.cpp:279-297): `-subject`'s
/// mapped row of BLASTX, which the tables write as it is (`tabular::TabularHsp::end_gaps`) and
/// outfmt 0 refuses (decision X2g-5 of line RG).
pub(crate) fn end_gap_span<T: Copy + PartialEq + From<u8>>(
    row: &[T],
    residue_dash: impl Fn(usize) -> bool,
) -> (usize, usize) {
    let gap = T::from(b'-');
    let end_gap = |column: usize| row[column] == gap && !residue_dash(column);
    let lo = (0..row.len())
        .find(|&column| !end_gap(column))
        .unwrap_or(row.len());
    let hi = (lo..row.len())
        .rev()
        .find(|&column| !end_gap(column))
        .map_or(lo, |column| column + 1);
    (lo, hi)
}

/// Check if two amino acids are a positive match (similar)
///
/// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:2122-2149
/// ```c
/// if(sequence_standard[i]==sequence[i]){
///     match ++;
/// } else {
///     if ((m_AlignType&eProt)
///         && m_Matrix[(int)sequence_standard[i]][(int)sequence[i]] > 0){
///         positive ++;
///         if(m_AlignOption & eShowMiddleLine){
///             middle_line[i] = '+';
///         }
///     }
/// }
/// ```
fn is_positive_match(a: char, b: char, matrix: ScoringMatrix) -> bool {
    if a == b {
        return true;
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3563-3581
    // ```c++
    //     retval.Resize(k_NumAsciiChar, k_NumAsciiChar, -1000);
    //
    //     SNCBIFullScoreMatrix mtx;
    //     NCBISM_Unpack(packed_mtx, &mtx);
    //
    //     for(int i = 0; i < ePMatrixSize; ++i){
    //         for(int j = 0; j < ePMatrixSize; ++j){
    //             retval((size_t)k_PSymbol[i], (size_t)k_PSymbol[j]) =
    //                 mtx.s[(size_t)k_PSymbol[i]][(size_t)k_PSymbol[j]];
    //         }
    //     }
    //     for(int i = 0; i < ePMatrixSize; ++i) {
    //         retval((size_t)k_PSymbol[i], '*') = retval('*',(size_t)k_PSymbol[i]) = -4;
    //     }
    //     retval('*', '*') = 1;
    //     // this is to count Selenocysteine to Cysteine matches as positive
    //     retval('U', 'U') = retval('C', 'C');
    //     retval('U', 'C') = retval('C', 'C');
    //     retval('C', 'U') = retval('C', 'C');
    // ```
    protein_display_score(matrix, a as u8, b as u8) > 0
}

// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showdefline.cpp:75-83
// ```c
// static const char*  kHeader = "Sequences producing significant alignments:";
// ...
// static const char*  kBits = (getenv("CTOOLKIT_COMPATIBLE") ? "(bits)" : "(Bits)");
// static const size_t kBits_size = strlen(kBits);
// ...
// static const char*  kValue = "Value";
// ```
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
    read: super::query_warnings::SubjectReadHook<'_, W>,
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
        read,
    )
}

/// The number of decimal digits of `value` (the length of `value.to_string()`).
#[inline]
fn digit_count(value: usize) -> usize {
    value.checked_ilog10().map_or(1, |log| log as usize + 1)
}

/// Width of the start-coordinate column of an alignment block.
///
/// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:1366-1367
/// ```c++
///         size_t maxCood=max<size_t>(m_AV->GetSeqStart(row), m_AV->GetSeqStop(row));
///         maxStartLen = max<size_t>(NStr::SizetToString(maxCood).size(), maxStartLen);
/// ```
/// `GetSeqStart` and `GetSeqStop` are 0-based, so the width is that of the largest
/// 1-based coordinate minus one, and a printed 1-based start can have one digit
/// more. NCBI then pads with the `size_t` value `maxStartLen - startLen +
/// k_StartSequenceMargin`, which is one space in that case; callers write
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
fn ensure_trailing_period(text: &str) -> String {
    if text.ends_with('.') {
        text.to_string()
    } else {
        format!("{text}.")
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:584-603
// ```c
// if (gapped) {
//     out << "Gapped" << "\n";
// }
// out << "Lambda      K        H";
// if (gbp) {
//     if (gapped) {
//         out << "        a         alpha    sigma";
//     } else {
//         out << "        a         alpha";
//     }
// }
// out << "\n";
// sprintf(buffer, "%#8.3g ", lambda);
// ```
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

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:387-395
// ```c++
//       CBlastFormatUtil::BlastPrintVersionInfo(m_Program, m_IsHTML,
//                                               m_Outfile);
//     }
//
//     if (m_IsBl2Seq && !m_IsDbScan) {
//         return;
//     }
//
//     m_Outfile << NcbiEndl << NcbiEndl;
// ```
/// The prolog of blastp's, psiblast's and rpsblast's pairwise reports before the database: the
/// version line (`reference::write_version_line`), two newlines and the references of
/// `PrintProlog` (`reference::write_prolog_references`: the Gapped BLAST reference, psiblast's
/// compositionally adjusted matrices and the composition-based statistics of psiblast and
/// blastp). `program` is `psiblast`, `rpsblast`, or blastp's for any other.
fn write_blastp_pairwise_intro<W: Write>(
    writer: &mut W,
    version: &str,
    program: &str,
) -> io::Result<()> {
    let program = match program {
        "psiblast" | "rpsblast" => program,
        _ => "blastp",
    };
    super::reference::write_version_line(writer, false, program, version)?;
    writeln!(writer)?;
    writeln!(writer)?;
    super::reference::write_prolog_references(
        writer,
        false,
        super::reference::DEFAULT_PROTOCOL,
        program,
        false,
    )
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blastfmtutil.cpp:125-138
// ```c
// string dbString = (html) ? "<b>Database:</b> " : "Database: ";
// str << dbString << definition_line << endl;
// ...
// out << "           " << NStr::IntToString(nNumSeqs,NStr::fWithCommas)
//     << " sequences; " << NStr::UInt8ToString(nTotalLength,NStr::fWithCommas)
//     << " total letters" << endl;
// ```
// NCBI reference: c++/src/corelib/ncbistr.cpp:5088-5340 (WrapIt, fWrap_FlatFile)
// enum EScore { eForced, ePunct, eComma, eSpace, eNewline };
// if (score >= best_score && score_pos > pos0) { best_pos = score_pos; best_score = score; }
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:288-293
// ```c++
//     NStr::Wrap(str, line_len, string_l, flags);
//     list<string>::iterator iter = string_l.begin();
//     while(iter != string_l.end())
//     {
//        out << *iter;
//        out << "\n";
// ```
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
    // NCBI reference: c++/src/corelib/ncbistr.cpp:5100,5153-5157,5137
    // SIZE_TYPE pos=0,len=str.size(),nl_pos=0; if(nl_pos<=pos) nl_pos=str.find('\n',pos);
    // bool thisPartHasBackspace = false;
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
            // NCBI reference: c++/src/corelib/ncbistr.cpp:5260-5266,5277-5300
            // if (best_pos != len) { best_pos=len; thisPartHasBackspace=true; }
            // if (thisPartHasBackspace) { /* eat backspaces and preceding characters */ }
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
/// NCBI writes before it reads the first query batch (blast_format.cpp, `PrintProlog`).
pub fn write_blastp_pairwise_prolog<W: Write>(
    writer: &mut W,
    version: &str,
    program: &str,
    database_name: &str,
    database_num_sequences: usize,
    database_total_letters: usize,
    database: Option<&PairwiseDatabase>,
) -> io::Result<()> {
    write_blastp_pairwise_intro(writer, version, program)?;
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
// ```c++
//     if (top) {
//         const CAlignFormatUtil::SDbInfo* dbinfo = &(dbinfo_list.front());
//         out << "Database: ";
//
//         string db_titles = dbinfo->definition;
//         Int8 tot_num_seqs = static_cast<Int8>(dbinfo->number_seqs);
//         Int8 tot_length = dbinfo->total_length;
//
//         for (size_t i = 1; i < dbinfo_list.size(); i++) {
//             db_titles += "; " + dbinfo_list[i].definition;
//             tot_num_seqs += static_cast<Int8>(dbinfo_list[i].number_seqs);
//             tot_length += dbinfo_list[i].total_length;
//         }
//
//         x_WrapOutputLine(db_titles, line_length, out);
//         ...
//         CAlignFormatUtil::AddSpace(out, 11);
//         out << NStr::Int8ToString(tot_num_seqs, NStr::fWithCommas) <<
//             " sequences; " <<
//             NStr::Int8ToString(tot_length, NStr::fWithCommas) <<
//             " total letters\n\n";
//         return;
//     }
// ```
/// The `Database:` lines of the prolog for the database of a `-db` search: the titles of the
/// databases of its list joined with `; ` (`PairwiseDatabase::new`), wrapped at 68, then the
/// summed counts and one blank line.
fn write_blast_database_header<W: Write>(
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

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:526-532,546-553
// ```c++
//         x_WrapOutputLine(db_titles, line_length, out);
//         if ( !dbinfo->filt_algorithm_name.empty() ) {
//             out << "Masked using: '" << dbinfo->filt_algorithm_name << "'";
//             if ( !dbinfo->filt_algorithm_options.empty() ) {
//                 out << ", options: '" << dbinfo->filt_algorithm_options << "'";
//             }
//             out << endl;
//         }
//     ...
//             x_WrapOutputLine(dbinfo->definition, line_length, out);
//
//             if ( !dbinfo->filt_algorithm_name.empty() ) {
//                 out << "  Masked using: '" << dbinfo->filt_algorithm_name << "'";
//                 if ( !dbinfo->filt_algorithm_options.empty() ) {
//                     out << ", options: '" << dbinfo->filt_algorithm_options << "'";
//                 }
//                 out << endl;
//             }
// ```
/// The `Masked using:` line of a database searched with database masking (`indent` before it:
/// none in the prolog, two spaces in the epilog), not wrapped; nothing without it.
fn write_masked_using<W: Write>(
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
    writeln!(writer)
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:541-574
// ```c++
//     ITERATE(vector<SDbInfo>, dbinfo, dbinfo_list) {
//         if (dbinfo->subset == false) {
//             out << "  Database: ";
//             x_WrapOutputLine(dbinfo->definition, line_length, out);
//             ...
//             out << "    Posted date:  ";
//             out << dbinfo->date << "\n";
//
//             out << "  Number of letters in database: ";
//             out << NStr::Int8ToString(dbinfo->total_length,
//                                       NStr::fWithCommas) << "\n";
//             out << "  Number of sequences in database:  ";
//             out << NStr::IntToString(dbinfo->number_seqs,
//                                      NStr::fWithCommas) << "\n";
//             ...
//         }
//         out << "\n";
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2255-2266
// ```c++
//     if ( !m_IsBl2Seq || m_IsDbScan) {
//         CBlastFormatUtil::PrintDbReport(m_DbInfo, kFormatLineLength,
//                                         m_Outfile, false);
//     }
//
//     if (m_Program == "blastn" || m_Program == "megablast") {
//         m_Outfile << "\n\nMatrix: " << "blastn matrix " <<
//     ...
//     else {
//         m_Outfile << "\n\nMatrix: " << options.GetMatrixName() << "\n";
//     }
// ```
/// The database blocks of the epilog, one per database of the list, and the blank lines
/// before the matrix line (as `write_final_database_report`).
fn write_final_database_block<W: Write>(
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
// ```c++
//         CAlignFormatUtil::AddSpace(out, 11);
//         out << NStr::Int8ToString(tot_num_seqs, NStr::fWithCommas) <<
//             " sequences; " <<
//             NStr::Int8ToString(tot_length, NStr::fWithCommas) <<
//             " total letters\n\n";
//         return;
// ```
fn write_blastp_database_header_spacing<W: Write>(
    writer: &mut W,
    database_name: &str,
    database_num_sequences: usize,
    database_total_letters: usize,
    trailing: usize,
) -> io::Result<()> {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:513-525
    // ```c++
    //         out << "Database: ";
    //
    //         string db_titles = dbinfo->definition;
    //         Int8 tot_num_seqs = static_cast<Int8>(dbinfo->number_seqs);
    //         Int8 tot_length = dbinfo->total_length;
    //
    //         for (size_t i = 1; i < dbinfo_list.size(); i++) {
    //             db_titles += "; " + dbinfo_list[i].definition;
    //             tot_num_seqs += static_cast<Int8>(dbinfo_list[i].number_seqs);
    //             tot_length += dbinfo_list[i].total_length;
    //         }
    //
    //         x_WrapOutputLine(db_titles, line_length, out);
    // ```
    write!(writer, "Database: ")?;
    write_flatfile_wrapped(writer, ensure_trailing_period(database_name).as_bytes(), 68)?;
    writeln!(
        writer,
        "           {} sequences; {} total letters",
        crate::ncbi_str::int8_with_commas(i64::try_from(database_num_sequences).unwrap_or(0)),
        crate::ncbi_str::int8_with_commas(i64::try_from(database_total_letters).unwrap_or(0))
    )?;
    for _ in 0..trailing {
        writeln!(writer)?;
    }
    Ok(())
}

/// The `Query=` block of a query of title `query_title` and length `query_length`, as NCBI's
/// insertions (`cli::insert`; the organism report's stream counts them, outfmt 0's buffers them).
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:729
// ```c++
//         out << label << "= ";
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:741-745
// ```c++
//         x_WrapOutputLine(all_id_str, line_len, out, html);
//         if(cbs.IsSetInst() && cbs.GetInst().CanGetLength()){
//             out << "\nLength=";
//             out << cbs.GetInst().GetLength() <<"\n";
//         }
// ```
pub(crate) fn write_blastp_query_header<W: Write + ?Sized>(
    writer: &mut W,
    query_title: &[u8],
    query_length: usize,
) -> io::Result<()> {
    // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:732-734
    // ```c++
    //     string all_id_str = GetSeqIdString(cbs, believe_query);
    //     all_id_str += " ";
    //     all_id_str = NStr::TruncateSpaces(all_id_str + GetSeqDescrString(cbs));
    // ```
    // The ID of a local query is not shown (`believe_query` false), so the line is the
    // title, raw (no decoding), without the white space at its ends.
    crate::cli::insert(writer, b"Query")?;
    crate::cli::insert(writer, b"= ")?;
    write_flatfile_wrapped(writer, truncate_spaces(query_title), 68)?;
    crate::cli::insert(writer, b"\nLength=")?;
    crate::cli::insert(writer, query_length.to_string().as_bytes())?;
    crate::cli::insert(writer, b"\n")
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1510-1514
// ```c
// m_Outfile << "\n\n"
//           << "***** " << CBlastFormatUtil::kNoHitsFound << " *****" << "\n"
//           << "\n\n";
// ```
fn write_no_hits_found<W: Write>(writer: &mut W) -> io::Result<()> {
    // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1510-1514
    // m_Outfile << "\n\n" << "***** " << kNoHitsFound << " *****\n\n\n";
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer, "***** No hits found *****")?;
    writeln!(writer)?;
    writeln!(writer)?;
    Ok(())
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:445-477
// ```c
// m_Outfile << NcbiEndl;
// if (kbp_ungap) {
//     CBlastFormatUtil::PrintKAParameters(..., false, gbp);
// }
// ...
// if (kbp_gap) {
//     CBlastFormatUtil::PrintKAParameters(..., true, gbp);
// }
// m_Outfile << "\n";
// m_Outfile << "Effective search space used: "
//           << summary.GetSearchSpace() << "\n";
// ```
// NCBI c++/src/algo/blast/api/local_blast.cpp:177-180,204-208:
// if (status != 0) {
//     pair<double, double> tmp_pair(-1.0, -1.0);
//     CRef<CBlastAncillaryData> tmp_ancillary_data(
//         new CBlastAncillaryData(tmp_pair, tmp_pair, tmp_pair, 0));
// }
// NCBI c++/src/algo/blast/format/blast_format.cpp:451-478 and
// c++/src/objtools/align_format/align_format_util.cpp:581-613:
// PrintKAParameters emits both -1 Karlin blocks with no Gumbel columns.
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:477-478
// ```c++
//     m_Outfile << "Effective search space used: " <<
//                         summary.GetSearchSpace() << "\n";
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1491-1491
// ```c++
//     m_Outfile << "\n\n";
// ```
fn write_tblastn_unsearched_query_footer_spacing<W: Write>(
    writer: &mut W,
    trailing: bool,
) -> io::Result<()> {
    writeln!(writer)?;
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

// NCBI c++/src/algo/blast/format/blast_format.cpp:445-478:
// if (kbp_ungap) CBlastFormatUtil::PrintKAParameters(..., false, gbp);
// if (kbp_gap) CBlastFormatUtil::PrintKAParameters(..., true, gbp);
// m_Outfile << "Effective search space used: " << summary.GetSearchSpace() << "\n";
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:477-478
// ```c++
//     m_Outfile << "Effective search space used: " <<
//                         summary.GetSearchSpace() << "\n";
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1491-1491
// ```c++
//     m_Outfile << "\n\n";
// ```
fn write_blastp_query_footer_spacing<W: Write>(
    writer: &mut W,
    ungapped_karlin: KarlinParams,
    gapped_karlin: KarlinParams,
    gumbel: Option<BlastGumbelBlk>,
    effective_search_space: i64,
    trailing: bool,
) -> io::Result<()> {
    let Some(gumbel) = gumbel else {
        return write_nucleotide_query_footer_spacing(
            writer,
            Some((ungapped_karlin, gapped_karlin)),
            true,
            effective_search_space,
            trailing,
        );
    };
    writeln!(writer)?;
    writeln!(writer, "Lambda      K        H        a         alpha")?;
    write_ncbi_ka_field(writer, ungapped_karlin.lambda)?;
    write_ncbi_ka_field(writer, ungapped_karlin.k)?;
    write_ncbi_ka_field(writer, ungapped_karlin.h)?;
    write_ncbi_ka_field(writer, gumbel.a_un)?;
    write_ncbi_ka_field(writer, gumbel.alpha_un)?;
    writeln!(writer)?;
    writeln!(writer)?;
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
// ```c++
//     else {
//         m_Outfile << "\n\nMatrix: " << options.GetMatrixName() << "\n";
//     }
// ```
// NCBI reference (GCC 13 libstdc++): include/c++/13/ostream:661-670
// ```c++
//     inline basic_ostream<char, _Traits>&
//     operator<<(basic_ostream<char, _Traits>& __out, const char* __s)
//     {
//       if (!__s)
// 	__out.setstate(ios_base::badbit);
//       else
// 	__ostream_insert(__out, __s,
// 			 static_cast<streamsize>(_Traits::length(__s)));
//       return __out;
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:96-103
// ```c++
// void
// CBlastRPSOptionsHandle::SetScoringOptionsDefaults()
// {
//     SetGappedMode();
//     // set invalid values for options that are not applicable
//     m_Opts->SetOutOfFrameMode(false);
//     m_Opts->SetFrameShiftPenalty(INT2_MAX);
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:73-83 (quoted at
// `rps_search::volumes::search_sets`: `CreateRpsStructures` names the database's matrix in the
// options of each RPS query batch it sets up)
/// The epilog's matrix line after its blank lines (`options.GetMatrixName()`). `None`: the
/// options name no matrix, as RPS-BLAST's (`CBlastRPSOptionsHandle` sets none) before a query
/// batch's set-up names the database's; the null `const char*` sets `badbit` on the formatter's
/// stream, which throws (`cli::BadbitInsertion`): the report ends after `Matrix: ` (NCBI BLAST+
/// 2.17.0: rpsblast and rpstblastn with an empty pipe, outfmt 0 to 4, exit status 6).
fn write_matrix_line(writer: &mut impl Write, matrix_name: Option<&str>) -> io::Result<()> {
    writer.write_all(b"Matrix: ")?;
    let Some(name) = matrix_name else {
        return Err(crate::cli::BadbitInsertion::error());
    };
    writeln!(writer, "{name}")
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/align_format_util.cpp:503-579
// ```c
// out << "  Database: ";
// ...
// out << "  Number of letters in database: ";
// out << NStr::Int8ToString(dbinfo->total_length, NStr::fWithCommas) << "\n";
// out << "  Number of sequences in database:  ";
// out << NStr::IntToString(dbinfo->number_seqs, NStr::fWithCommas) << "\n";
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:2258-2285
// ```c
// m_Outfile << "\n\nMatrix: " << options.GetMatrixName() << "\n";
// ...
// m_Outfile << "Gap Penalties: Existence: "
//         << options.GetGapOpeningCost() << ", Extension: "
//         << gap_extension << "\n";
// ...
// m_Outfile << "Neighboring words threshold: "
//         << options.GetWordThreshold() << "\n";
// m_Outfile << "Window for multiple hits: "
//         << options.GetWindowSize() << "\n";
// ```
fn write_blastp_final_footer<W: Write>(
    writer: &mut W,
    report: &BlastpPairwiseReport,
    database: Option<&PairwiseDatabase>,
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
    writeln!(
        writer,
        "Gap Penalties: Existence: {}, Extension: {}",
        report.gap_open, report.gap_extend
    )?;
    if report.word_threshold != 0.0 {
        // `GetWordThreshold()` is a double, written with the stream's default format.
        writeln!(
            writer,
            "Neighboring words threshold: {}",
            cpp_default_double(report.word_threshold)
        )?;
    }
    if report.window_size != 0 {
        writeln!(writer, "Window for multiple hits: {}", report.window_size)?;
    }
    Ok(())
}

/// The database block of the epilog and the blank lines before the matrix line
/// (align_format_util.cpp:503-579 and blast_format.cpp:2258-2262, quoted at
/// `write_blastp_final_footer`).
fn write_final_database_report<W: Write>(
    writer: &mut W,
    database_name: &str,
    database_num_sequences: usize,
    database_total_letters: usize,
) -> io::Result<()> {
    // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:543-544
    // out << "  Database: "; x_WrapOutputLine(dbinfo->definition, line_length, out);
    write!(writer, "  Database: ")?;
    write_flatfile_wrapped(writer, ensure_trailing_period(database_name).as_bytes(), 68)?;
    writeln!(writer, "    Posted date:  Unknown")?;
    writeln!(
        writer,
        "  Number of letters in database: {}",
        crate::ncbi_str::int8_with_commas(i64::try_from(database_total_letters).unwrap_or(0))
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

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:372-424
// ```c
// CBlastFormatUtil::BlastPrintVersionInfo(m_Program, m_IsHTML, m_Outfile);
// ...
// CBlastFormatUtil::BlastPrintReference(...);
// ...
// CBlastFormatUtil::BlastPrintReference(..., CReference::eCompBasedStats, ...);
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1520-1589
// ```c
// if ( (!m_IsBl2Seq || m_IsDbScan) && !(m_DisableKAStats || kIsGlobal) ) {
//     x_DisplayDeflines(aln_set, itr_num, prev_seqids);
// }
// ...
// display.DisplaySeqalign(m_Outfile);
// x_PrintOneQueryFooter(*results.GetAncillaryData());
// ```
/// The report of `queries`: the prolog when `prolog` is set (EXPERIMENT (LOSAT_X_BPQBLOCK):
/// false for the query blocks after the first of blastp's block search, whose first block wrote
/// it), each query's section, and the epilog when `epilog` is set. `query_gapped_k` is
/// rpsblast's gapped K of each query's footer (`rps_search::engine::batch_gapped_k`), `None` for
/// blastp and psiblast (the report's `gapped_karlin`).
pub fn write_blastp_pairwise_report<W: Write>(
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
        loads,
    )
}

/// The queries `query_range` of a BLASTP pairwise report, whose `hits` are those queries' (the
/// rows of a caller that makes them a window of queries at a time, as BLASTN's and TBLASTX's
/// do, `x_pairwise_par::query_windows`): the prolog with the report's first query, the epilog
/// with its last. `anchored`: the sequences of the query-anchored sections (`config.view` of
/// outfmt 1 to 4; `None` for outfmt 0).
#[allow(clippy::too_many_arguments)]
pub fn write_blastp_pairwise_window<W: Write>(
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
    loads: &crate::report::db_subject::HitLoads,
) -> io::Result<()> {
    let mut buffered = io::BufWriter::new(writer);
    let writer = &mut buffered;

    if prolog && query_range.start == 0 {
        write_blastp_pairwise_intro(writer, &report.version, &config.program)?;
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

    // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1411
    // ```c++
    // CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
    // ```
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

    // NCBI reference: c++/src/app/blast/blastp_app.cpp:294-295
    // ```c++
    //         BLAST_PROF_START( APP.POST );
    //         formatter.PrintEpilog(opt);
    // ```
    // A search that stops at a query batch (`Empty CBlastQueryVector` with `-query_loc`)
    // writes no epilog.
    if !(epilog && query_range.end == queries.len()) {
        return writer.flush();
    }

    // The blank lines before the epilog (blast_format.cpp:2249).

    writeln!(writer)?;

    writeln!(writer)?;

    write_blastp_final_footer(writer, report, database)?;
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
    /// What the reports' loads do at each subject, by subject (`thrown_load`).
    loads: &'a crate::report::db_subject::HitLoads,
}

/// The HSPs of a query-anchored section (outfmt 1 to 4), added to `hsps`: the first
/// `num_alignments` subjects of the query's sorted hits (`QueryGroup::aligned`), each with all its
/// HSPs in their order. An HSP's Dense-seg is its Seq-align's, which `segments_of` makes and gives
/// its rows' strands and query frame of (`seqalign::protein_row_segments`,
/// `nucleotide_row_segments`, `translated_row_segments`); the subjects shown in one hit (copies of
/// a database record, `joins`) are the hit's subject, one Bioseq for NCBI's scope; `is_query`
/// tells the subject that is the query's own Bioseq.
fn section_hsps(
    mut hsps: super::query_anchored::Hsps,
    num_alignments: usize,
    group: &super::x_pairwise_par::QueryGroup<'_>,
    is_query: impl Fn(u32) -> bool,
    mut segments_of: impl FnMut(
        &PairwiseHit,
        &mut super::serial::seqalign::Segments,
    ) -> super::query_anchored::HspStrands,
) -> super::query_anchored::Hsps {
    let mut segments = super::serial::seqalign::Segments::default();
    for &s_idx in group.aligned(num_alignments) {
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
    is_query: impl Fn(u32) -> bool,
) -> super::query_anchored::Hsps {
    section_hsps(
        super::query_anchored::Hsps::default(),
        num_alignments,
        group,
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
    is_query: impl Fn(u32) -> bool,
) -> super::query_anchored::Hsps {
    section_hsps(
        super::query_anchored::Hsps::default(),
        num_alignments,
        group,
        is_query,
        |hit, segments| super::serial::seqalign::nucleotide_row_segments(&hit.hit, segments),
    )
}

/// `section_hsps` of a BLASTX section (`seqalign::translated_row_segments`: a nucleotide query row
/// of width 3 in the hit's frame of the query's record, a protein subject row; an HSP of length 0,
/// whose rows the writer reaches as the engine made them, `blastx_mapped_empty_rows`: one segment
/// of length 0 without widths, `seqalign::empty_row_segments`).
fn blastx_section_hsps(
    num_alignments: usize,
    group: &super::x_pairwise_par::QueryGroup<'_>,
) -> super::query_anchored::Hsps {
    section_hsps(
        super::query_anchored::Hsps::of_translated_query(),
        num_alignments,
        group,
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

/// Whether the description table of the subjects `descriptions` (their hits in `grouped`) reads
/// an HSP with empty rows (`PairwiseHit::rows_empty`): it reads the query and subject ranges of
/// every HSP of each described subject before its first line (`x_InitDeflineTable` ->
/// `GetSeqAlignSetCalcParams` -> `s_ProcessAlignSet`), which throw `Row is empty` on such an HSP.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:2913-2921
/// ```c++
///     ITERATE(CSeq_align_set::Tdata, iter, alnset.Get()) {
///         CRange<TSeqPos> query_range = (*iter)->GetSeqRange(0);
///         ...
///         CRange<TSeqPos> subject_range = (*iter)->GetSeqRange(1);
/// ```
fn describes_empty_rows(
    descriptions: &[u32],
    grouped: &std::collections::HashMap<u32, Vec<&PairwiseHit>>,
) -> bool {
    descriptions
        .iter()
        .any(|oid| grouped[oid].iter().any(|hit| hit.rows_empty()))
}

/// Whether a BLASTX subject's HSPs `hits` hold one with empty rows whose Std-seg `-subject`'s
/// Seq-loc mapper (`remaps`) rewrote (`mapper_rewrites`): its query row of one nucleotide against
/// the empty subject row, on which `CreateDensegFromStdseg` throws (`EmptyRowsStop::StdSegLengths`).
fn blastx_mapped_empty_rows(remaps: bool, hits: &[&PairwiseHit]) -> bool {
    hits.iter().any(|hit| hit.rows_empty()) && mapper_rewrites(remaps, hits)
}

/// Where NCBI BLAST+'s report of a BLASTX query whose hits are `group` (`BlastxSections`' group) in
/// the query-anchored view of `config` stops at an HSP with empty rows, in the order it meets them:
/// its sorts (`sort_stop`), its description table (`describes_empty_rows`), the Dense-segs of its
/// aligned subjects (`blastx_mapped_empty_rows`, with `-subject`, `remaps`), all of which throw
/// (`write_blastx_pairwise_report` writes NCBI's report up to there and its exception), and the
/// merge of its section, which frees a segment that it still reads
/// (`EmptyRowsStop::MergeFreesSegment`, `query_anchored::section_frees`: the split of a segment at
/// an HSP of length 0 at a codon of an alignment in the query's reading frame +1, or
/// `FillUnalignedRegions` after a segment of length 0). `None` for a report without such an
/// HSP among its aligned subjects. What NCBI writes before the crash is what its stream had flushed
/// (no prefix of it is NCBI's, as for `query_anchored::ncbi_aborts`): the program refuses the
/// report before the call's first byte (`empty_rows::UNPORTED_ANCHORED_FREE`), unless a query
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
    if describes_empty_rows(&order[..order.len().min(options.num_descriptions)], grouped) {
        return Ok(Some(EmptyRowsStop::RowIsEmpty));
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
        |group| blastp_section_hsps(num_alignments, group, &is_query),
    )
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
// ```c++
//     for (CSeq_align_set::Tdata::const_iterator
//              iter = m_AlnSetRef->Get().begin();
//          iter != m_AlnSetRef->Get().end() && num_align < m_NumToShow;
//          iter++, ialn++){
//         ...
//         subid = &((*iter)->GetSeq_id(1));
//         if(is_first_aln || (!is_first_aln && !subid->Match(*previous_id))) {
//             SScoreInfo* sci = x_GetScoreInfo(**iter, static_cast<int>(num_align));
// ```
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
        match self.query_gapped_k {
            Some(gapped_k) => KarlinParams {
                k: gapped_k[q_idx],
                ..self.report.gapped_karlin
            },
            None => self.report.gapped_karlin,
        }
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

    fn head<W: Write>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: Option<&super::x_pairwise_par::QueryGroup<'_>>,
        read: super::query_warnings::SubjectReadHook<'_, W>,
    ) -> io::Result<()> {
        let (query, report) = (&self.queries[q_idx], self.report);
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1474-1476
        // ```c++
        //     if (itr_num != numeric_limits<unsigned int>::max()) {
        //         m_Outfile << "Results from round " << itr_num << "\n";
        //     }
        // ```
        // psiblast's results carry their round (`PrintOneResultSet` with
        // `itr.GetIterationNumber()`, psiblast_app.cpp:283-285); LOSAT's psiblast runs its
        // first round only.
        if self.config.program == "psiblast" {
            writeln!(writer, "Results from round 1")?;
        }
        // NCBI blast_format.cpp:1491: m_Outfile << "\n\n";
        writeln!(writer)?;
        writeln!(writer)?;
        // The title's bytes by query index (`query_name` is BLASTX's).
        let query_title = self
            .query_titles
            .get(q_idx)
            .map_or(&b""[..], |title| title.as_ref());
        write_blastp_query_header(writer, query_title, query.query_length)?;
        let Some(group) = group else {
            write_no_hits_found(writer)?;
            // NCBI c++/src/algo/blast/api/local_blast.cpp:177-180,204-208:
            // if (m_PrelimSearch->CheckInternalData() != 0)
            //     new CBlastAncillaryData(tmp_pair, tmp_pair, tmp_pair, 0);
            // The all-invalid batch uses -1 sentinel blocks. An invalid query
            // within a searched batch has null blocks (blast_results.cpp:82-103),
            // so only the blank lines and the zero search space are written
            // (blast_format.cpp:445-477), as TBLASTN does.
            if query.batch_skipped {
                write_tblastn_unsearched_query_footer_spacing(writer, false)?;
            } else if !query.valid {
                writeln!(writer)?;
                writeln!(writer)?;
                writeln!(writer)?;
                writeln!(writer, "Effective search space used: 0")?;
            } else {
                write_blastp_query_footer_spacing(
                    writer,
                    query.ungapped_karlin,
                    self.gapped_karlin(q_idx),
                    report.gumbel,
                    query.effective_search_space,
                    false,
                )?;
            }
            return Ok(());
        };
        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:606-611
        // ```c++
        //     CShowBlastDefline showdef(*aln_set, *m_Scope,
        //                               defline_length == -1 ? kFormatLineLength:defline_length,
        //                               m_NumSummary + additional);
        // ```
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
            read,
        )?;
        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1520-1589
        // x_DisplayDeflines(aln_set, ...); ... display.DisplaySeqalign(m_Outfile);
        writeln!(writer)
    }

    fn subject<W: Write>(
        &self,
        writer: &mut W,
        probe: Option<&mut FormatProbe<'_>>,
        group: &super::x_pairwise_par::QueryGroup<'_>,
        s_idx: u32,
        part: std::ops::Range<usize>,
    ) -> io::Result<()> {
        let config = self.config;
        let parsed = parsed_subject(self.parsed_subjects, s_idx);
        if self.database.is_some() || parsed.is_some() {
            let heading = match parsed {
                Some(subject) => DatabaseHeading::loaded(
                    subject,
                    true,
                    group.own_hsps(s_idx),
                    self.thrown(s_idx).is_some(),
                ),
                None => DatabaseHeading::Title(
                    &group.subject_titles[&s_idx],
                    stored_source(self.subject_sources, s_idx),
                ),
            };
            return write_database_subject_alignments_part(
                writer,
                probe,
                heading,
                true,
                &group.subject_hits[&s_idx],
                &group.subject_hit_indices[&s_idx],
                part,
                |writer, hit| {
                    write_hsp_info(writer, hit, config)?;
                    write_alignment(writer, hit, config)
                },
            );
        }
        write_subject_alignments_part(
            writer,
            probe,
            s_idx,
            self.parsed_subjects,
            &group.subject_titles[&s_idx],
            true,
            &group.subject_hits[&s_idx],
            &group.subject_hit_indices[&s_idx],
            part,
            |writer, hit| {
                write_hsp_info(writer, hit, config)?;
                write_alignment(writer, hit, config)
            },
        )
    }

    fn thrown(&self, s_idx: u32) -> Option<(&crate::ncbi_exception::NcbiException, bool)> {
        thrown_load(self.loads, self.parsed_subjects, s_idx)
    }

    fn view(&self) -> PairwiseView {
        self.config.view
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1550-1552
    // ```c++
    //     CSeq_align_set copy_aln_set;
    //     CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_NumAlignments);
    // ```
    /// The query-anchored section of the query's aligned subjects (`blastp_section_hsps`).
    fn alignments<W: Write>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: &super::x_pairwise_par::QueryGroup<'_>,
    ) -> io::Result<()> {
        write_anchored_section(writer, self.config, self.anchored, q_idx, |sequences| {
            blastp_section_hsps(self.report.num_alignments, group, |s_idx| {
                sequences.subject_is_query(q_idx, s_idx)
            })
        })
    }

    fn foot<W: Write>(&self, writer: &mut W, q_idx: usize) -> io::Result<()> {
        let query = &self.queries[q_idx];
        write_blastp_query_footer_spacing(
            writer,
            query.ungapped_karlin,
            self.gapped_karlin(q_idx),
            self.report.gumbel,
            query.effective_search_space,
            false,
        )
    }
}

// =============================================================================
// BLASTN pairwise report
// =============================================================================

/// One BLASTN query of the pairwise report (NCBI `CBlastAncillaryData` and the query's
/// validity).
///
/// NCBI reference: c++/src/algo/blast/api/blast_results.cpp:82-104
/// ```c++
///     // find the first valid context corresponding to this query
///     ...
///     m_SearchSpace = ctx->eff_searchsp;
///     ...
///     s_InitializeKarlinBlk(sbp->kbp_std[ctx_index], &m_UngappedKarlinBlk);
/// ```
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
    /// `-task megablast` (NCBI's `m_Megablast`; the program name stays `blastn`).
    pub megablast: bool,
    pub database_name: String,
    pub database_num_sequences: usize,
    pub database_total_letters: usize,
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
    pub window_size: usize,
    /// A gapped search (`options.GetGappedMode()`); an ungapped (`-ungapped`) search has no
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
}

// The description table of a report: BLASTN's, BLASTP's, TBLASTN's and TBLASTX's
// (`write_blastn_description_table`), BLASTX's
// (`write_subject_summary_table_with_sum_n`), which differ only in the line that names a
// subject (`label_of`).
//
// NCBI reference: c++/src/objtools/align_format/showdefline.cpp:1004-1011
// ```c++
// void CShowBlastDefline::DisplayBlastDefline(CNcbiOstream & out)
// {
//     x_InitDeflineTable();
//     ...
//     x_DisplayDefline(out);
// ```
// NCBI reference: c++/src/objtools/align_format/showdefline.cpp:1100-1142
// ```c++
//         if(!is_first_aln && !(subid->Match(*previous_id))) {
//             SScoreInfo* sci = x_GetScoreInfoForTable(hit, num_align);
//             if(sci){
//                 m_ScoreList.push_back(sci);
//                 if(m_MaxScoreLen < sci->bit_string.size()){
//                     m_MaxScoreLen = sci->bit_string.size();
//                 }
//                 if(m_MaxTotalScoreLen < sci->total_bit_string.size()){
//                     m_MaxTotalScoreLen = sci->total_bit_string.size();
//                 }
//     ...
//         if (num_align < m_NumToShow) { //no adding if number to show already reached
//             hit.Set().push_back(*iter);
//         }
//     ...
//     //the last hit
//     SScoreInfo* sci = x_GetScoreInfoForTable(hit, num_align);
//     if(sci){
//          m_ScoreList.push_back(sci);
//         if(m_MaxScoreLen < sci->bit_string.size()){
//             m_MaxScoreLen = sci->bit_string.size();
//         }
//         if(m_MaxTotalScoreLen < sci->total_bit_string.size()){
//             m_MaxScoreLen = sci->total_bit_string.size();
//         }
// ```
// NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:4285-4310
// ```c++
//         total_bits += bits;
//     ...
//         if (bits > highest_bits) {
//             highest_bits = bits;
//             lowest_evalue = evalue;
//         }
//     ...
//     seqSetInfo->total_bit_score = total_bits;
//     seqSetInfo->bit_score = highest_bits;
//     seqSetInfo->evalue = lowest_evalue;
// ```
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
// NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:4229-4230,4259
// ```c++
//     seqSetInfo->sum_n = sum_n == -1 ? 1:sum_n ;
//     ...
//     seqSetInfo = GetSeqAlignCalcParams(*(aln.Get().front()));
// ```
// NCBI reference: c++/src/objtools/align_format/showdefline.cpp:1061,1118-1120
// ```c++
//     m_MaxSumNLen =1;
//     ...
//                 if( m_MaxSumNLen < NStr::IntToString(sci->sum_n).size()){
//                     m_MaxSumNLen = NStr::IntToString(sci->sum_n).size();
//                 }
// ```
// With `show_sum_n` (NCBI's `eShowSumN`: ungapped sum statistics, blast_format.cpp:165-167,
// 503-504) a last column shows the linked-set size of the subject's first HSP; the HSP
// without a "sum_n" score (`num` of 1) shows 1.
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:771-791
// ```c++
//             out << kScore;
//             CAlignFormatUtil::AddSpace(out, m_MaxScoreLen - kScore_size);
//             CAlignFormatUtil::AddSpace(out, kTwoSpaceMargin_size);
//
//
//             if (m_Option & eShowTotalScore) {
//                 out << kTotal;
//                 CAlignFormatUtil::AddSpace(out, m_MaxTotalScoreLen - kTotal_size);
//                 CAlignFormatUtil::AddSpace(out, kTwoSpaceMargin_size);
//             }
//             if (m_Option & eShowQueryCoverage) {
//                 out << kQueryCov;
//                 CAlignFormatUtil::AddSpace(out, kTwoSpaceMargin_size);
//             }
//
//             CAlignFormatUtil::AddSpace(out, 2); //E align to l of value
//             out << kE;
//
//             if (m_Option & eShowPercentIdent) {
//                 CAlignFormatUtil::AddSpace(out, m_MaxEvalueLen - kValue_size);
//                 CAlignFormatUtil::AddSpace(out, kTwoSpaceMargin_size);
//                 CAlignFormatUtil::AddSpace(out, kTwoSpaceMargin_size);
//                 CAlignFormatUtil::AddSpace(out, kOneSpaceMargin_size);
//                 out << kMax;//"Max" - "ident" -second line
//             }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:809-840
// ```c++
//             CAlignFormatUtil::AddSpace(out, m_LineLen - kHeader_size);
//             CAlignFormatUtil::AddSpace(out, kOneSpaceMargin_size);
//             out << kBits;
//             //in case m_MaxScoreLen > kBits.size()
//             CAlignFormatUtil::AddSpace(out, m_MaxScoreLen - kBits_size);
//             CAlignFormatUtil::AddSpace(out, kTwoSpaceMargin_size);
//
//
//             if (m_Option & eShowTotalScore) {
//                 CAlignFormatUtil::AddSpace(out, kOneSpaceMargin_size);
//                 out << kScoreLine2;//"score"
//                 CAlignFormatUtil::AddSpace(out, m_MaxTotalScoreLen - kTotal_size);
//                 CAlignFormatUtil::AddSpace(out, kTwoSpaceMargin_size);
//             }
//             if (m_Option & eShowQueryCoverage) {
//                 out << kQueryCovLine2;//"cov"
//                 CAlignFormatUtil::AddSpace(out, kTwoSpaceMargin_size);
//                 CAlignFormatUtil::AddSpace(out, kOneSpaceMargin_size);
//             }
//
//             out << kValue;
//             if((m_Option & eShowSumN) || (m_Option & eShowPercentIdent)){
//                 CAlignFormatUtil::AddSpace(out, m_MaxEvalueLen - kValue_size);
//                 CAlignFormatUtil::AddSpace(out, kTwoSpaceMargin_size);
//             }
//             if(m_Option & eShowSumN){
//                 out << kN;
//             }
//             if (m_Option & eShowPercentIdent) {
//                 out << kIdentLine2;//"ident"
//             }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:66-104,1051-1065
// ```c++
// static const char*  kTotal = "Total";
// ...
// static const char*  kIdentity = "Max ident";
// static const size_t kIdentity_size = strlen(kIdentity);
// ...
// static const char*  kMax = "Max";
// static const char*  kIdentLine2 = "Ident";
// static const char*  kScoreLine2 = "Score";
// static const char*  kQueryCov = "Query";
// static const size_t kQueryCov_size = strlen(kQueryCov);
// static const char*  kQueryCovLine2 = "cover";
// ...
//     m_MaxTotalScoreLen = kTotal_size;
//     m_MaxPercentIdentityLen = kIdentity_size;
// ```
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
    read: super::query_warnings::SubjectReadHook<'_, W>,
) -> io::Result<()> {
    // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:513-514
    // ```c++
    //     if (num_descriptions_to_show == 0)
    //         flags |= CShowBlastDefline::eNoShowHeader;
    // ```
    // c++/src/objtools/align_format/showdefline.cpp:757,848: with `-num_descriptions 0`
    // the header block (the two header lines and the blank line) and the rows are not written;
    // the blank line that closes the table is (`out << "\n"` at blast_format.cpp:613).
    if described.is_empty() {
        return writeln!(writer);
    }
    let rows: Vec<(u32, &PairwiseHit, f64, i32, Option<SortColumnValues>)> = described
        .iter()
        .filter_map(|s_idx| {
            let hits = subject_hits.get(s_idx)?;
            let first = *hits.first()?;
            let mut best = first;
            for hit in hits {
                if hit.hit.bit_score > best.hit.bit_score {
                    best = hit;
                }
            }
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
    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/showdefline.cpp:1059
    // ```c
    // m_MaxScoreLen = kBits_size;
    // ```
    let mut max_score = ncbi_k_bits().len();
    let mut max_evalue = "Value".len();
    let mut max_total = "Total".len();
    let mut max_sum_n = 1;
    for (index, (_, best, total, sum_n, _)) in rows.iter().enumerate() {
        max_score = max_score.max(format_bitscore_ncbi(best.hit.bit_score).len());
        max_evalue = max_evalue.max(format_evalue_ncbi(best.hit.e_value).len());
        max_sum_n = max_sum_n.max(sum_n.to_string().len());
        let total_len = format_total_bitscore_ncbi(*total).len();
        if index + 1 == rows.len() && last_row_counted {
            if max_total < total_len {
                max_score = total_len;
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:940-944
                // ```c++
                //         if (m_Option & eShowTotalScore) {
                //             out << kTwoSpaceMargin << kOneSpaceMargin << (*iter)->total_bit_string;
                //             CAlignFormatUtil::AddSpace(out, m_MaxTotalScoreLen -
                //                                    (*iter)->total_bit_string.size());
                //         }
                // ```
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:932-938
                // ```c++
                // void CAlignFormatUtil::AddSpace(CNcbiOstream& out, size_t number)
                //
                // {
                //     for(auto i=0; i<number; i++){
                //         out<<" ";
                //     }
                // }
                // ```
                // The last row's total, wider than the column (above 99999 bits, `%5.3le`),
                // pads with `m_MaxTotalScoreLen - size` spaces as a `size_t` below 0: NCBI
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
    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/showdefline.cpp:811-813
    // ```c
    // out << kBits;
    // //in case m_MaxScoreLen > kBits.size()
    // CAlignFormatUtil::AddSpace(out, m_MaxScoreLen - kBits_size);
    // ```
    write!(writer, "{:<max_score$}  ", ncbi_k_bits())?;
    if columns.is_some() {
        write!(writer, " {:<max_total$}  cover   ", "Score")?;
    }
    write!(writer, "Value")?;
    // NCBI reference: c++/src/objtools/align_format/showdefline.cpp:829-836
    // ```c++
    //             out << kValue;
    //             if((m_Option & eShowSumN) || (m_Option & eShowPercentIdent)){
    //                 CAlignFormatUtil::AddSpace(out, m_MaxEvalueLen - kValue_size);
    //                 CAlignFormatUtil::AddSpace(out, kTwoSpaceMargin_size);
    //             }
    //             if(m_Option & eShowSumN){
    //                 out << kN;
    //             }
    // ```
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
    for (s_idx, best, total, sum_n, values) in &rows {
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1587-1589
        // ```c++
        //     try{
        //         const CBioseq_Handle& handle = m_ScopeRef->GetBioseqHandle(*id);
        //         x_FillDeflineAndId(handle, *id, use_this_seqid, sdl, blast_rank);
        // ```
        // Each line loads its subject (`x_GetDeflineInfo` of `x_DisplayDeflineTableBody`) after
        // the lines before it are written (`QueryWarnings::subject_read`).
        read(*s_idx, writer)?;
        label_of(*s_idx).write_cut(writer)?;
        write!(
            writer,
            "  {:<max_score$}",
            format_bitscore_ncbi(best.hit.bit_score)
        )?;
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:940-956
        // ```c++
        //         if (m_Option & eShowTotalScore) {
        //             out << kTwoSpaceMargin << kOneSpaceMargin << (*iter)->total_bit_string;
        //             CAlignFormatUtil::AddSpace(out, m_MaxTotalScoreLen -
        //                                    (*iter)->total_bit_string.size());
        //         }
        //
        //         if (m_Option & eShowQueryCoverage) {
        //             //int percent_coverage = 100*(*iter)->master_covered_length/m_QueryLength;
        //             int percent_coverage = (*iter)->percent_coverage;
        //
        //             out << kTwoSpaceMargin << percent_coverage << "%";
        //             //minus one due to % sign
        //             CAlignFormatUtil::AddSpace(out, kQueryCov_size -
        //                                        NStr::IntToString(percent_coverage).size() - 1);
        //         }
        // ```
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
            format_evalue_ncbi(best.hit.e_value)
        )?;
        // NCBI reference: c++/src/objtools/align_format/showdefline.cpp:961-965
        // ```c++
        //         if(m_Option & eShowSumN){
        //             out << kTwoSpaceMargin << (*iter)->sum_n;
        //             CAlignFormatUtil::AddSpace(out, m_MaxSumNLen -
        //                      NStr::IntToString((*iter)->sum_n).size());
        //         }
        // ```
        if show_sum_n {
            write!(writer, "  {sum_n:<max_sum_n$}")?;
        }
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:966-974
        // ```c++
        //         if(m_Option & eShowPercentIdent){
        //             int percent_identity =(int) (0.5 + (*iter)->percent_identity);
        //             if(percent_identity > 100) {
        //                 percent_identity = min(99, percent_identity);
        //             }
        //             out << kTwoSpaceMargin << percent_identity <<"%";
        //             CAlignFormatUtil::AddSpace(out, m_MaxPercentIdentityLen -
        //                                            NStr::IntToString(percent_identity).size());
        //         }
        // ```
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
// ```c++
//         if (m_Option & eShowTotalScore) {
//             out << kTwoSpaceMargin << kOneSpaceMargin << (*iter)->total_bit_string;
//             CAlignFormatUtil::AddSpace(out, m_MaxTotalScoreLen -
//                                    (*iter)->total_bit_string.size());
//         }
// ```
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
    /// The shown identity: `(int)(0.5 + percent_identity)`, 99 above 100.
    percent_identity: i64,
}

impl SortColumnValues {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4271-4306
    // ```c++
    //     seqSetInfo->subjRange = CAlignFormatUtil::GetSeqAlignCoverageParams(aln,&seqSetInfo->master_covered_length,&seqSetInfo->flip);
    //     seqSetInfo->percent_coverage = (int)(0.5 + 100*seqSetInfo->master_covered_length/(double)queryLength);
    //     if(seqSetInfo->percent_coverage > 100) {
    //         seqSetInfo->percent_coverage = min(99, seqSetInfo->percent_coverage);
    //     }
    //
    //     ITERATE(CSeq_align_set::Tdata, iter, aln.Get()) {
    //         int align_length = CAlignFormatUtil::GetAlignmentLength(**iter, do_translation);
    //         ...
    //         if (bits > highest_bits) {  // this is the replacement condition (WB-1175)
    //             highest_length = align_length;
    //             highest_ident = num_ident;
    // ///         highest_identity = 100*num_ident/align_length;
    //         }
    //     ...
    //     seqSetInfo->match = highest_ident;
    //     seqSetInfo->align_length = highest_length;
    //     seqSetInfo->percent_identity = CAlignFormatUtil::GetPercentIdentity(seqSetInfo->match, seqSetInfo->align_length);
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1525-1530
    // ```c++
    //     unique_ptr<CAlignFormatUtil::SSeqAlignSetCalcParams> seqSetInfo( CAlignFormatUtil::GetSeqAlignSetCalcParamsFromASN(aln));
    //     if(seqSetInfo->hspNum == 0) {//calulated params are not in ASN - calculate now
    //         seqSetInfo.reset( CAlignFormatUtil::GetSeqAlignSetCalcParams(aln,m_QueryLength,m_TranslatedNucAlignment));
    //     }
    // ```
    /// The columns of the subject whose HSPs are `hits` and whose best HSP (the first of the
    /// highest bit score, `write_description_table`) is `best`: `GetSeqAlignSetCalcParams` with
    /// `m_TranslatedNucAlignment` false (TBLASTX's identities over its nucleotide lengths). The
    /// cover is computed without NCBI's `int` product `100*master_covered_length`, which
    /// overflows above about 21 million covered letters.
    fn of(hits: &[&PairwiseHit], best: &PairwiseHit, columns: SortColumns) -> Self {
        use super::hit_sort::{percent_identity, query_coverage, HspScores};
        let scores: Vec<HspScores> = hits.iter().map(|hit| HspScores::of(hit)).collect();
        let covered = query_coverage(&scores);
        let mut percent_coverage =
            (0.5 + (100 * covered) as f64 / columns.query_length as f64) as i64;
        if percent_coverage > 100 {
            percent_coverage = percent_coverage.min(99);
        }
        let best = HspScores::of(best);
        let identity = percent_identity(
            best.num_ident,
            best.alignment_length(columns.tblastx, false),
        );
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
fn parsed_subject(
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    s_idx: u32,
) -> Option<&crate::report::db_subject::ParsedSubject> {
    parsed
        .get(s_idx as usize)
        .and_then(Option::as_deref)
        .filter(|subject| !subject.shown_as_record())
}

/// The id of NCBI's line in place of each alignment of `-subject` record `s_idx` whose heading
/// it cannot make (`write_subject_alignments_part`): `GetSeqIdString` of the Seq-align id, the local id
/// `Subject_N` (N = `s_idx` + 1) of a record without parsed ids, else its parsed id (`parsed`:
/// `ParsedSubject::skipped_id`, `Subject_9` of `lcl|Subject_9`).
fn record_skipped_id(
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    s_idx: u32,
) -> Vec<u8> {
    match parsed.get(s_idx as usize).and_then(Option::as_deref) {
        Some(subject) => subject.skipped_id(),
        None => format!("Subject_{}", u64::from(s_idx) + 1).into_bytes(),
    }
}

/// A description line's label in NCBI's two parts (`x_DisplayDeflineTableBody`): the id text
/// (`gi|N|` and the id, the `line_length` that NCBI counts before the title; under `-html` the
/// text of the sequence's link) and the title (`line_component`: the defline, after two spaces
/// where a gi was written; under `-html` HTML-encoded). The plain line writes them as one string
/// (`write_cut`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:883-913
/// ```c++
///         if((m_Option & eHtml) && (sdl->id_url != NcbiEmptyString)) {
///             out << sdl->id_url;
///         }
///         if(m_Option & eShowGi){
///             if(sdl->gi > ZERO_GI){
///                 line_component = "gi|" + NStr::NumericToString(sdl->gi) + "|";
///                 out << line_component;
///                 line_length += line_component.size();
///             }
///         }
///         if(!sdl->id.Empty()){
///             ...
///                 out << idStr;
///                 line_length += idStr.size();
///             }
///         }
///         if((m_Option & eHtml) && (sdl->id_url != NcbiEmptyString)) {
///             out << "</a>";
///         }
///         line_component = (line_component.empty() ? "" : "  ") + sdl->defline;
/// ```
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DescriptionLabel {
    /// The id text (`ParsedSubject::description`'s start; empty for a title alone).
    pub id: Vec<u8>,
    /// The title.
    pub title: Vec<u8>,
}

impl DescriptionLabel {
    /// A label of a title alone (a generated defline, `Unknown`).
    pub fn title_only(title: Vec<u8>) -> Self {
        DescriptionLabel {
            id: Vec::new(),
            title,
        }
    }

    /// The plain label on the description line: the id text and the title as one string cut to
    /// 68 columns (a longer one is its first 65 bytes and `...`), or padded to them.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:915-918,930
    /// ```c++
    ///         if(line_component.size()+line_length > m_LineLen){
    ///             actual_line_component = line_component.substr(0, m_LineLen -
    ///                                                           line_length - 3);
    ///             actual_line_component += kEllipsis;
    ///     ...
    ///         CAlignFormatUtil::AddSpace(out, m_LineLen - line_length);
    /// ```
    /// String widths are byte counts, including non-ASCII FASTA titles. (NCBI cuts the title
    /// alone, which is the same cut for an id text of up to 65 bytes; a longer id is refused
    /// before any output (`common/database.rs`, `report/db_subject.rs`), so the cut of the
    /// whole label below never cuts the id.)
    fn write_cut<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        // A label of a title alone is written as the joined label was: one write of its letters.
        let len = self.id.len() + self.title.len();
        if len > 68 {
            let id = self.id.len().min(65);
            if id > 0 {
                writer.write_all(&self.id[..id])?;
            }
            writer.write_all(&self.title[..65 - id])?;
            writer.write_all(b"...")
        } else {
            if !self.id.is_empty() {
                writer.write_all(&self.id)?;
            }
            writer.write_all(&self.title)?;
            write_spaces(writer, 68 - len)
        }
    }
}

/// The description-table label of a database sequence (before the cut to the line): the title
/// that `GenerateDefline` makes of the stored title `title` and the source `source` of the
/// Bioseq, which an empty title is made from (`Unknown` where NCBI's `HtmlDecode` throws), or for
/// a sequence with parsed Seq-ids its id and title (`ParsedSubject::description`; where its load
/// throws, `thrown`, the `catch`'s Seq-align id and `Unknown`, `ParsedSubject::thrown_description`).
pub fn database_description_label(
    title: &[u8],
    source: &super::defline::SourceTaxname,
    parsed: Option<&crate::report::db_subject::ParsedSubject>,
    protein: bool,
    thrown: bool,
) -> DescriptionLabel {
    if let Some(subject) = parsed {
        let (id, title) = if thrown {
            subject.thrown_description()
        } else {
            subject.description(protein)
        };
        return DescriptionLabel {
            id,
            title: title.text,
        };
    }
    let facts = super::defline::BioseqFacts::of_source(source.name().ok().flatten());
    DescriptionLabel::title_only(
        match super::defline::generate_bioseq_defline(title, protein, true, &facts) {
            Ok(defline) => defline.text,
            Err(super::defline::UnknownEncoding) => b"Unknown".to_vec(),
        },
    )
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
    read: super::query_warnings::SubjectReadHook<'_, W>,
    loads: &crate::report::db_subject::HitLoads,
) -> io::Result<()> {
    write_description_table(
        writer,
        described,
        last_row_counted,
        subject_hits,
        show_sum_n,
        columns,
        |s_idx| {
            // NCBI reference: c++/src/objtools/align_format/showdefline.cpp:498
            // The description keeps the prefixes (`fLeavePrefixSuffix`, `defline.rs`).
            // NCBI reference: c++/src/objtools/align_format/showdefline.cpp:1579-1592
            // ```c++
            //     sdl->defline = "Unknown";
            //     ...
            //     try{
            //         const CBioseq_Handle& handle = m_ScopeRef->GetBioseqHandle(*id);
            //         x_FillDeflineAndId(handle, *id, use_this_seqid, sdl, blast_rank);
            //     } catch (const CException&){
            //         sdl->defline = "Unknown";
            // ```
            // (also the `catch` of a load that throws, `loads`, `thrown_description`)
            let title = subject_titles.get(&s_idx).map_or(&[][..], Vec::as_slice);
            database_description_label(
                title,
                stored_source(sources, s_idx),
                parsed_subject(parsed, s_idx),
                protein,
                loads.thrown(s_idx as usize).is_some(),
            )
        },
        read,
    )
}

// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:387-407
// ```c++
//       CBlastFormatUtil::BlastPrintVersionInfo(m_Program, m_IsHTML,
//                                               m_Outfile);
//     ...
//     m_Outfile << NcbiEndl << NcbiEndl;
//     ...
//     if (m_Megablast)
//         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
//                                           m_Outfile, CReference::eMegaBlast);
// ```
// The prolog of a megablast search: the reference is the megablast one
// (`reference::write_prolog_references`); the layout is that of
// `write_translated_pairwise_intro`.
fn write_megablast_pairwise_intro(writer: &mut impl Write, version: &str) -> io::Result<()> {
    super::reference::write_version_line(writer, false, "blastn", version)?;
    writeln!(writer)?;
    writeln!(writer)?;
    super::reference::write_prolog_references(
        writer,
        false,
        super::reference::DEFAULT_PROTOCOL,
        "blastn",
        true,
    )
}

// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:451-478
// ```c++
//     m_Outfile << NcbiEndl;
//     if (kbp_ungap) {
//         CBlastFormatUtil::PrintKAParameters(kbp_ungap->Lambda, kbp_ungap->K,
//                                             kbp_ungap->H, kFormatLineLength,
//                                             m_Outfile, false, gbp);
//     }
//
//     m_Outfile << "\n";
//     if (kbp_gap) {
//         ...
//     }
//
//     m_Outfile << "Effective search space used: " <<
//         summary.GetSearchSpace() << "\n";
// ```
// NCBI reference: c++/src/algo/blast/core/blast_setup.c:472-480
// ```c
//     if (program_number == eBlastTypeBlastn ||
//     ...
//         /* disable new FSC rules for nucleotide case for now */
//         if (sbp && sbp->gbp) {
//             sfree(sbp->gbp);
// ```
// Without a Gumbel block, PrintKAParameters (align_format_util.cpp:585-619) prints only
// the Lambda, K and H columns. An invalid query has neither Karlin block (`karlin` is
// `None`) and a search space of 0.
/// The footer of a nucleotide query, with the two blank lines that follow it (the next
/// query's preamble or the epilog) when `trailing`.
fn write_nucleotide_query_footer_spacing<W: Write>(
    writer: &mut W,
    karlin: Option<(KarlinParams, KarlinParams)>,
    gapped: bool,
    effective_search_space: i64,
    trailing: bool,
) -> io::Result<()> {
    writeln!(writer)?;
    if let Some((ungapped, _)) = karlin {
        writeln!(writer, "Lambda      K        H")?;
        write_ncbi_ka_field(writer, ungapped.lambda)?;
        write_ncbi_ka_field(writer, ungapped.k)?;
        write_ncbi_ka_field(writer, ungapped.h)?;
        writeln!(writer)?;
    }
    writeln!(writer)?;
    // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:1130-1159
    // ```c++
    //     if (kbp_ungap) ... PrintKAParameters(kbp_ungap->Lambda, kbp_ungap->K, kbp_ungap->H, ...,  false)
    //     if (kbp_gap) ... PrintKAParameters(..., true)   ("Gapped")
    // ```
    // An ungapped search has no gapped Karlin block (`kbp_gap` is not made).
    if let Some((_, gapped_block)) = karlin.filter(|_| gapped) {
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

// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:2261-2279
// ```c++
//     if (m_Program == "blastn" || m_Program == "megablast") {
//         m_Outfile << "\n\nMatrix: " << "blastn matrix " <<
//                         options.GetMatchReward() << " " <<
//                         options.GetMismatchPenalty() << "\n";
//     }
//     ...
//     double gap_extension = (double) options.GetGapExtensionCost();
//     if ((m_Program == "megablast" || m_Program == "blastn") && options.GetGapExtensionCost() == 0)
//     {
//         // Formula from PMID 10890397 applies if both gap values are zero.
//         gap_extension = -2*options.GetMismatchPenalty() + options.GetMatchReward();
//         gap_extension /= 2.0;
//     }
//     m_Outfile << "Gap Penalties: Existence: "
//             << options.GetGapOpeningCost() << ", Extension: "
//             << gap_extension << "\n";
// ```
// The word threshold and the window size are 0 for blastn (blast_format.cpp:2281-2288),
// so those lines are not printed. The extension is a whole or half number, which the
// C++ stream and Rust print alike (`2`, `2.5`).
fn write_blastn_final_footer<W: Write>(
    writer: &mut W,
    report: &BlastnPairwiseReport,
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
    // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:2270-2288
    // ```c
    //     if (options.GetGappedMode() == true) {
    //         double gap_extension = (double) options.GetGapExtensionCost();
    //         if ((m_Program == "megablast" || m_Program == "blastn") && options.GetGapExtensionCost() == 0)
    //         { // Formula from PMID 10890397 applies if both gap values are zero.
    //                gap_extension = -2*options.GetMismatchPenalty() + options.GetMatchReward();
    //                gap_extension /= 2.0;
    //         }
    //         m_Outfile << "Gap Penalties: Existence: "
    //                 << options.GetGapOpeningCost() << ", Extension: "
    //                 << gap_extension << "\n";
    //     }
    //     if (options.GetWordThreshold()) {
    //         m_Outfile << "Neighboring words threshold: " <<
    //                         options.GetWordThreshold() << "\n";
    //     }
    //     if (options.GetWindowSize()) {
    //         m_Outfile << "Window for multiple hits: " <<
    //                         options.GetWindowSize() << "\n";
    //     }
    // ```
    // BLASTN's word threshold is 0 for every task (BLAST_WORD_THRESHOLD_BLASTN and
    // BLAST_WORD_THRESHOLD_MEGABLAST, blast_options.h).
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
    if report.window_size != 0 {
        writeln!(writer, "Window for multiple hits: {}", report.window_size)?;
    }
    Ok(())
}

/// The start of the BLASTN report, which NCBI writes before it searches: the program
/// and its reference, then the subjects (the subject set's name and counts, or the
/// `Database:` lines of a database, `write_blast_database_header`).
///
/// NCBI reference: c++/src/app/blast/blastn_app.cpp:258-261
/// ```c
///         formatter.PrintProlog();
///
///         /*** Process the input ***/
///         CBatchSizeMixer mixer(SplitQuery_GetChunkSize(opt.GetProgram())-1000);
/// ```
pub fn write_blastn_pairwise_prolog<W: Write>(
    writer: &mut W,
    version: &str,
    megablast: bool,
    database_name: &str,
    database_num_sequences: usize,
    database_total_letters: usize,
    database: Option<&PairwiseDatabase>,
) -> io::Result<()> {
    if megablast {
        write_megablast_pairwise_intro(writer, version)?;
    } else {
        write_translated_pairwise_intro(writer, "BLASTN", version)?;
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
/// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1491-1589
/// ```c++
///     m_Outfile << "\n\n";
///     CBlastFormatUtil::AcknowledgeBlastQuery(*bioseq, kFormatLineLength,
///     ...
///         m_Outfile << "\n\n"
///                   << "***** " << CBlastFormatUtil::kNoHitsFound << " *****" << "\n"
///                   << "\n\n";
///         x_PrintOneQueryFooter(*results.GetAncillaryData());
///     ...
///         x_DisplayDeflines(aln_set, itr_num, prev_seqids);
///     ...
///     m_Outfile << "\n";
///     ...
///     CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_NumAlignments);
///     ...
///     display.DisplaySeqalign(m_Outfile);
///     ...
///     x_PrintOneQueryFooter(*results.GetAncillaryData());
/// ```
// Kept out of line: this runs only for outfmt 0 or hit records, and inlining it into the
// shared post-processing makes every BLASTN run compile it in Wasm hosts.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub fn write_blastn_pairwise_report<W: Write>(
    hits: &[PairwiseHit],
    writer: &mut W,
    config: &PairwiseConfig,
    queries: &[BlastnPairwiseQuery],
    query_range: std::ops::Range<usize>,
    subject_titles: &[Arc<[u8]>],
    report: &BlastnPairwiseReport,
    probe: Option<&mut FormatProbe<'_>>,
    warnings: Option<&mut super::query_warnings::QueryWarnings<'_>>,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
    anchored: Option<&dyn super::query_anchored::ReportSequences>,
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
            report.megablast,
            &report.database_name,
            report.database_num_sequences,
            report.database_total_letters,
            report.database.as_ref(),
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
            loads,
        },
        query_range.start,
        &hits_by_query,
        subject_titles,
        probe,
        warnings,
        x_pool,
    )?;

    // NCBI reference: ncbi-blast/c++/src/app/blast/blastn_app.cpp:317-318
    // ```c
    //         BLAST_PROF_START( APP.POST );
    //         formatter.PrintEpilog(opt);
    // ```
    // An error in a later query batch skips the epilog, and nothing follows the last
    // query's footer. The epilog starts with two blank lines (blast_format.cpp:2249).
    if report.epilog && query_range.end == queries.len() {
        writer.write_all(b"\n\n")?;
        write_blastn_final_footer(writer, report)?;
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

    fn head<W: Write>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: Option<&super::x_pairwise_par::QueryGroup<'_>>,
        read: super::query_warnings::SubjectReadHook<'_, W>,
    ) -> io::Result<()> {
        let (query, report) = (&self.queries[q_idx], self.report);
        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1489-1491
        // ```c++
        //     // print the preamble for this query
        //
        //     m_Outfile << "\n\n";
        // ```
        writer.write_all(b"\n\n")?;
        write_blastp_query_header(writer, &query.query_name, query.query_length)?;
        let Some(group) = group else {
            write_no_hits_found(writer)?;
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
        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:606-611
        // ```c++
        //     CShowBlastDefline showdef(*aln_set, *m_Scope,
        //                               defline_length == -1 ? kFormatLineLength:defline_length,
        //                               m_NumSummary + additional);
        // ```
        // The titles of the subjects of the records read (by subject index).
        write_blastn_description_table(
            writer,
            group.described(report.num_descriptions),
            group.subject_order.len() <= report.num_descriptions,
            &group.subject_hits,
            &group.subject_titles,
            &report.parsed_subjects,
            &report.subject_sources,
            false,
            false,
            self.config.sort_columns(query.query_length),
            read,
            self.loads,
        )?;
        writeln!(writer)
    }

    fn subject<W: Write>(
        &self,
        writer: &mut W,
        probe: Option<&mut FormatProbe<'_>>,
        group: &super::x_pairwise_par::QueryGroup<'_>,
        s_idx: u32,
        part: std::ops::Range<usize>,
    ) -> io::Result<()> {
        let config = self.config;
        // A database sequence's heading is its stored title after `>` (a `BL_ORD_ID`
        // id is not shown), as BLASTP's (`write_database_subject_alignments_part`), or the lines
        // of the definition lines of one with parsed Seq-ids.
        let parsed = parsed_subject(&self.report.parsed_subjects, s_idx);
        if self.report.database.is_some() || parsed.is_some() {
            let heading = match parsed {
                Some(subject) => DatabaseHeading::loaded(
                    subject,
                    self.report.gapped,
                    group.own_hsps(s_idx),
                    self.thrown(s_idx).is_some(),
                ),
                None => DatabaseHeading::Title(
                    &group.subject_titles[&s_idx],
                    stored_source(&self.report.subject_sources, s_idx),
                ),
            };
            return write_database_subject_alignments_part(
                writer,
                probe,
                heading,
                false,
                &group.subject_hits[&s_idx],
                &group.subject_hit_indices[&s_idx],
                part,
                |writer, hit| {
                    write_hsp_info(writer, hit, config)?;
                    write_alignment(writer, hit, config)
                },
            );
        }
        write_subject_alignments_part(
            writer,
            probe,
            s_idx,
            &self.report.parsed_subjects,
            &group.subject_titles[&s_idx],
            false,
            &group.subject_hits[&s_idx],
            &group.subject_hit_indices[&s_idx],
            part,
            |writer, hit| {
                write_hsp_info(writer, hit, config)?;
                write_alignment(writer, hit, config)
            },
        )
    }

    fn thrown(&self, s_idx: u32) -> Option<(&crate::ncbi_exception::NcbiException, bool)> {
        thrown_load(self.loads, &self.report.parsed_subjects, s_idx)
    }

    fn view(&self) -> PairwiseView {
        self.config.view
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1550-1552
    // ```c++
    //     CSeq_align_set copy_aln_set;
    //     CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_NumAlignments);
    // ```
    /// The query-anchored section of the query's aligned subjects (`blastn_section_hsps`).
    fn alignments<W: Write>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: &super::x_pairwise_par::QueryGroup<'_>,
    ) -> io::Result<()> {
        write_anchored_section(writer, self.config, self.anchored, q_idx, |sequences| {
            blastn_section_hsps(self.report.num_alignments, group, |s_idx| {
                sequences.subject_is_query(q_idx, s_idx)
            })
        })
    }

    fn foot<W: Write>(&self, writer: &mut W, q_idx: usize) -> io::Result<()> {
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
/// `config` (outfmt 1 to 4), of the HSPs that `hsps_of` makes with the report's sequences
/// `anchored` (`QuerySections::alignments` of each program).
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
    super::query_anchored::write_section(
        writer,
        super::query_anchored::View { identity, anchored },
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
/// NCBI reference: c++/src/algo/blast/api/blast_results.cpp:82-100
/// ```c++
///     // find the first valid context corresponding to this query
///     for (i = 0; i < context_per_query; i++) {
///         BlastContextInfo *ctx = query_info->contexts +
///                                 query_number * context_per_query + i;
///         if (ctx->is_valid) {
///             m_SearchSpace = ctx->eff_searchsp;
///     ...
///     if (sbp->kbp_std) {
///         s_InitializeKarlinBlk(sbp->kbp_std[ctx_index], &m_UngappedKarlinBlk);
///     }
/// ```
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
    pub database_name: String,
    pub database_num_sequences: usize,
    pub database_total_letters: usize,
    /// The `-matrix` name as typed (`options.GetMatrixName()`), BLOSUM62 when omitted.
    pub matrix_name: String,
    /// `-threshold` and `-window_size` (the epilog prints them when they are not 0).
    pub word_threshold: f64,
    pub window_size: usize,
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
}

/// The start of the TBLASTX report, which NCBI writes before it searches: the program and
/// its reference, then the subjects.
///
/// NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:173-178
/// ```c
///         formatter.PrintProlog();
///
///         /*** Process the input ***/
///         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
///
///             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
/// ```
pub fn write_tblastx_pairwise_prolog<W: Write>(
    writer: &mut W,
    version: &str,
    database_name: &str,
    database_num_sequences: usize,
    database_total_letters: usize,
    database: Option<&PairwiseDatabase>,
) -> io::Result<()> {
    write_translated_pairwise_intro(writer, "TBLASTX", version)?;
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

// NCBI reference: c++/src/objtools/align_format/showalign.cpp:3593-3606
// ```c++
//             out<<" Score = "<<bit_score_buf<<" ";
//             out<<"bits ("<<aln_vec_info->score<<"),"<<"  ";
//             out<<"Expect";
//             if (aln_vec_info->sum_n > 0) {
//             out << "(" << aln_vec_info->sum_n << ")";
//             }
//             out << " = " << evalue_buf;
//     ...
//     out << "\n";
// ```
// NCBI reference: c++/src/objtools/align_format/showalign.cpp:310-332
// ```c++
//     out<<" Identities = "<<match<<"/"<<(aln_stop+1)<<" ("<<identity<<"%"<<")";
//     if(aln_is_prot){
//         out<<", Positives = "<<(positive + match)<<"/"<<(aln_stop+1)
// 			<<" ("<<CAlignFormatUtil::GetPercentMatch(positive + match, aln_stop+1)<<"%"<<")";
//     }
//     out<<", Gaps = "<<gap<<"/"<<(aln_stop+1)
//        <<" ("<<CAlignFormatUtil::GetPercentMatch(gap, aln_stop+1)<<"%"<<")"<<"\n";
//     ...
//     if(master_frame != 0 && slave_frame != 0) {
//         out <<" Frame = " << ((master_frame > 0) ? "+" : "")
//             << master_frame <<"/"<<((slave_frame > 0) ? "+" : "")
//             << slave_frame<<"\n";
//     ...
//     out<<"\n";
// ```
// NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:4229-4230
// ```c++
//     seqSetInfo->sum_n = sum_n == -1 ? 1:sum_n ;
// ```
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

// NCBI reference: c++/src/objtools/align_format/showalign.cpp:1598-1626
// ```c++
//     int start = alnRoInfo->seqStarts[row].front() + 1;  //+1 for 1 based
//     int end = alnRoInfo->seqStops[row].front() + 1;
//     ...
//         out << start;
//         startLen=NStr::IntToString(start).size();
//     ...
//     CAlignFormatUtil::AddSpace(out, alnRoInfo->maxStartLen-startLen + k_StartSequenceMargin);
//     x_OutputSeq(alnRoInfo->sequence[row], m_AV->GetSeqId(row), j,
//     ...
//     CAlignFormatUtil::AddSpace(out, k_SeqStopMargin);
//     ...
//         out << end;
//     ...
//     out<<"\n";
// ```
// NCBI reference: c++/src/objtools/alnmgr/alnvec.cpp:230,305,358
// ```c++
//     scrn_width *= width;
//     ...
//                         scrn_lft_seq_pos = plus ? start : stop;
//     ...
//                     scrn_rgt_seq_pos = plus ? stop : start;
// ```
// NCBI reference: c++/src/objtools/align_format/showalign.cpp:1776-1784
// ```c++
//     CAlignFormatUtil:: AddSpace(out, alnRoInfo->maxIdLen + k_IdStartMargin + alnRoInfo->maxStartLen + k_StartSequenceMargin);
//     x_OutputSeq(alnRoInfo->middleLine, no_id, j, (int)actualLineLen, 0, row, false, alnRoInfo->masked_regions[row], out);
// ```
// NCBI reference: c++/src/objtools/align_format/showalign.cpp:2131-2150
// ```c++
//         if(sequence_standard[i]==sequence[i]){
//             if(m_AlignOption & eShowMiddleLine){
//                 if(m_MidLineStyle == eBar ) {
//                     middle_line[i] = '|';
//                 } else if (m_MidLineStyle == eChar){
//                     middle_line[i] = sequence[i];
//                 }
//             }
//             match ++;
//         } else {
//             if ((m_AlignType&eProt)
//                 && m_Matrix[(int)sequence_standard[i]][(int)sequence[i]] > 0){
//                 positive ++;
//                 if(m_AlignOption & eShowMiddleLine){
//                     if (m_MidLineStyle == eChar)
//                         middle_line[i] = '+';
//                 }
// ```
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
        q.len(),
        width,
        |writer, row, columns| writer.write_all(&[q, s][row][columns]),
        |writer, columns| write_protein_middle(writer, q, s, columns, config.protein_matrix, false),
    )
}

// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:451-478
// ```c++
//     m_Outfile << NcbiEndl;
//     if (kbp_ungap) {
//         CBlastFormatUtil::PrintKAParameters(kbp_ungap->Lambda,
//                                             kbp_ungap->K, kbp_ungap->H,
//                                             kFormatLineLength, m_Outfile,
//                                             false, gbp);
//     }
//     ...
//     m_Outfile << "\n";
//     if (kbp_gap) {
//     ...
//     m_Outfile << "\n";
//     m_Outfile << "Effective search space used: " <<
//                         summary.GetSearchSpace() << "\n";
// ```
// NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:583-603
// ```c++
//     char buffer[256];
//     if (gapped) {
//         out << "Gapped" << "\n";
//     }
//     out << "Lambda      K        H";
//     if (gbp) {
//     ...
//     out << "\n";
//     sprintf(buffer, "%#8.3g ", lambda);
//     out << buffer;
// ```
// TBLASTX is ungapped and has no Gumbel block: only the ungapped Karlin block of the
// query's first valid context is printed. A query without a valid context in a searched
// batch prints neither block and a search space of 0 (blast_results.cpp:92-94).
fn write_tblastx_query_footer<W: Write>(
    writer: &mut W,
    query: &TblastxPairwiseQuery,
) -> io::Result<()> {
    writeln!(writer)?;
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

// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:2249,2256-2259
// ```c++
//     m_Outfile << NcbiEndl << NcbiEndl;
//     ...
//     if ( !m_IsBl2Seq || m_IsDbScan) {
//         CBlastFormatUtil::PrintDbReport(m_DbInfo, kFormatLineLength,
//                                         m_Outfile, false);
//     }
// ```
// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:2266-2288
// ```c++
//         m_Outfile << "\n\nMatrix: " << options.GetMatrixName() << "\n";
//     }
//
//     if (options.GetGappedMode() == true) {
//     ...
//     }
//     if (options.GetWordThreshold()) {
//         m_Outfile << "Neighboring words threshold: " <<
//                         options.GetWordThreshold() << "\n";
//     }
//     if (options.GetWindowSize()) {
//         m_Outfile << "Window for multiple hits: " <<
//                         options.GetWindowSize() << "\n";
//     }
// ```
// TBLASTX is ungapped, so the gap penalties are not printed.
fn write_tblastx_epilog<W: Write>(
    writer: &mut W,
    report: &TblastxPairwiseReport,
    database: Option<&PairwiseDatabase>,
) -> io::Result<()> {
    writer.write_all(b"\n\n")?;
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
    if report.word_threshold != 0.0 {
        // `GetWordThreshold()` is a double, written with the stream's default format.
        writeln!(
            writer,
            "Neighboring words threshold: {}",
            cpp_default_double(report.word_threshold)
        )?;
    }
    if report.window_size != 0 {
        writeln!(writer, "Window for multiple hits: {}", report.window_size)?;
    }
    Ok(())
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
/// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1491-1589
/// ```c++
///     m_Outfile << "\n\n";
///     CBlastFormatUtil::AcknowledgeBlastQuery(*bioseq, kFormatLineLength,
///     ...
///         m_Outfile << "\n\n"
///                   << "***** " << CBlastFormatUtil::kNoHitsFound << " *****" << "\n"
///                   << "\n\n";
///         x_PrintOneQueryFooter(*results.GetAncillaryData());
///     ...
///         x_DisplayDeflines(aln_set, itr_num, prev_seqids);
///     ...
///     m_Outfile << "\n";
///     ...
///     CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_NumAlignments);
///     ...
///     display.DisplaySeqalign(m_Outfile);
///     ...
///     x_PrintOneQueryFooter(*results.GetAncillaryData());
/// ```
// Kept out of line: this runs only for outfmt 0.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub fn write_tblastx_pairwise_report<W: Write>(
    hits: &[PairwiseHit],
    writer: &mut W,
    config: &PairwiseConfig,
    queries: &[TblastxPairwiseQuery],
    query_range: std::ops::Range<usize>,
    subject_titles: &[Arc<[u8]>],
    report: &TblastxPairwiseReport,
    database: Option<&PairwiseDatabase>,
    probe: Option<&mut FormatProbe<'_>>,
    warnings: Option<&mut super::query_warnings::QueryWarnings<'_>>,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
    // The sequences of the query-anchored sections (outfmt 1 to 4, `config.view`; `None` for
    // outfmt 0).
    anchored: Option<&dyn super::query_anchored::ReportSequences>,
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
            loads,
        },
        query_range.start,
        &hits_by_query,
        subject_titles,
        probe,
        warnings,
        x_pool,
    )?;

    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:209
    // ```c
    //         formatter.PrintEpilog(opt);
    // ```
    // An error in a later query batch skips the epilog, and nothing follows the last
    // query's footer.
    if report.epilog && query_range.end == queries.len() {
        write_tblastx_epilog(writer, report, database)?;
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

    fn head<W: Write>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: Option<&super::x_pairwise_par::QueryGroup<'_>>,
        read: super::query_warnings::SubjectReadHook<'_, W>,
    ) -> io::Result<()> {
        let (query, report) = (&self.queries[q_idx], self.report);
        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1489-1491
        // ```c++
        //     // print the preamble for this query
        //
        //     m_Outfile << "\n\n";
        // ```
        writer.write_all(b"\n\n")?;
        write_blastp_query_header(writer, &query.query_name, query.query_length)?;
        let Some(group) = group else {
            write_no_hits_found(writer)?;
            if report.unsearched.get(q_idx).copied().unwrap_or(false) {
                write_tblastn_unsearched_query_footer_spacing(writer, false)?;
            } else {
                write_tblastx_query_footer(writer, query)?;
            }
            return Ok(());
        };
        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:606-613
        // ```c++
        //     CShowBlastDefline showdef(*aln_set, *m_Scope,
        //                               defline_length == -1 ? kFormatLineLength:defline_length,
        //                               m_NumSummary + additional);
        //     ...
        //     m_Outfile << "\n";
        // ```
        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:165-167
        // ```c++
        //     if (use_sum_statistics && m_IsUngappedSearch) {
        //         m_ShowLinkedSetSize = true;
        //     }
        // ```
        // The titles of the subjects of the records read (by subject index), as BLASTN's.
        write_blastn_description_table(
            writer,
            group.described(report.num_descriptions),
            group.subject_order.len() <= report.num_descriptions,
            &group.subject_hits,
            &group.subject_titles,
            &report.parsed_subjects,
            &report.subject_sources,
            true,
            false,
            self.config.sort_columns(query.query_length),
            read,
            self.loads,
        )?;
        writeln!(writer)
    }

    fn subject<W: Write>(
        &self,
        writer: &mut W,
        probe: Option<&mut FormatProbe<'_>>,
        group: &super::x_pairwise_par::QueryGroup<'_>,
        s_idx: u32,
        part: std::ops::Range<usize>,
    ) -> io::Result<()> {
        let config = self.config;
        // NCBI reference: c++/src/objtools/align_format/showalign.cpp:3971-3980
        // ```c++
        // 		x_ShowAlnvecInfo(out,aln_vec_info,show_defline);
        // 	...
        //     out<<"\n";
        // ```
        // The query's record length (`Length=`), also with `-query_loc`: the hits of the query's
        // group are its hits (`x_pairwise_par::write` groups them by `Hit::q_idx`).
        let write_hsp = |writer: &mut W, hit: &PairwiseHit| {
            let query_length = self.queries[hit.hit.q_idx as usize].query_length;
            write_tblastx_hsp_info(writer, hit)?;
            write_tblastx_alignment(writer, hit, query_length, config)?;
            writeln!(writer)
        };
        // A database sequence's heading is its stored title after `>` (as BLASTP's,
        // BLASTN's and TBLASTN's `-db` reports), or the lines of the definition lines of one
        // with parsed Seq-ids.
        let parsed = parsed_subject(&self.report.parsed_subjects, s_idx);
        if self.database.is_some() || parsed.is_some() {
            let heading = match parsed {
                Some(subject) => DatabaseHeading::loaded(
                    subject,
                    false,
                    group.own_hsps(s_idx),
                    self.thrown(s_idx).is_some(),
                ),
                None => DatabaseHeading::Title(
                    &group.subject_titles[&s_idx],
                    stored_source(&self.report.subject_sources, s_idx),
                ),
            };
            return write_database_subject_alignments_part(
                writer,
                probe,
                heading,
                false,
                &group.subject_hits[&s_idx],
                &group.subject_hit_indices[&s_idx],
                part,
                write_hsp,
            );
        }
        write_subject_alignments_part(
            writer,
            probe,
            s_idx,
            &self.report.parsed_subjects,
            &group.subject_titles[&s_idx],
            false,
            &group.subject_hits[&s_idx],
            &group.subject_hit_indices[&s_idx],
            part,
            write_hsp,
        )
    }

    fn thrown(&self, s_idx: u32) -> Option<(&crate::ncbi_exception::NcbiException, bool)> {
        thrown_load(self.loads, &self.report.parsed_subjects, s_idx)
    }

    fn view(&self) -> PairwiseView {
        self.config.view
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1550-1552
    // ```c++
    //     CSeq_align_set copy_aln_set;
    //     CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_NumAlignments);
    // ```
    /// The query-anchored section of the query's aligned subjects: Dense-segs of two translated
    /// nucleotide rows in their frames of their records (`seqalign::translated_row_segments`),
    /// grouped by the query's frame and added with `fForceTranslation`
    /// (`eTranslateNucToNucAlignment`, `Hsps::of_translated_nucleotides`).
    fn alignments<W: Write>(
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
                |s_idx| sequences.subject_is_query(q_idx, s_idx),
                |hit, segments| tblastx_row_segments(hit, query_length, segments),
            )
        })
    }

    fn foot<W: Write>(&self, writer: &mut W, q_idx: usize) -> io::Result<()> {
        write_tblastx_query_footer(writer, &self.queries[q_idx])
    }
}

// =============================================================================
// Main pairwise output function
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// `write_flatfile_wrapped_with` finds the break of a line without newlines or backspaces from the window's
    /// end; it equals the scan of every column (the form it takes for a text that may hold a
    /// backspace) on texts of every scored class, at every width.
    #[test]
    fn the_break_found_from_the_end_is_the_scans() {
        let alphabet = b"ab ,-\n\t\x0bxyz  \x08";
        let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for _ in 0..20_000 {
            let len = (next() % 160) as usize;
            let backspace = next() % 4 == 0;
            let text: Vec<u8> = (0..len)
                .map(|_| alphabet[(next() % alphabet.len() as u64) as usize])
                .map(|c| if c == 8 && !backspace { b'q' } else { c })
                .collect();
            let width = 1 + (next() % 70) as usize;
            let mut found = Vec::new();
            write_flatfile_wrapped(&mut found, &text, width).unwrap();
            let mut scanned = Vec::new();
            let every_column = TextFacts {
                newlines: true,
                backspaces: true,
            };
            write_flatfile_wrapped_with(&mut scanned, &text, width, every_column).unwrap();
            assert_eq!(
                found,
                scanned,
                "{:?} at width {width}",
                String::from_utf8_lossy(&text)
            );
        }
    }

    /// The plain description label is the id text and the title as one string, cut to 68 columns
    /// or padded to them, wherever the parts split it (line RG, item T4a0).
    #[test]
    fn a_description_label_in_parts_writes_the_joined_line() {
        let cut = |id: &[u8], title: &[u8]| {
            let mut out = Vec::new();
            DescriptionLabel {
                id: id.to_vec(),
                title: title.to_vec(),
            }
            .write_cut(&mut out)
            .unwrap();
            String::from_utf8(out).unwrap()
        };
        let joined = |label: &[u8]| cut(b"", label);
        assert_eq!(cut(b"P1 ", b"title"), format!("P1 title{}", " ".repeat(60)));
        assert_eq!(joined(b"P1 title"), cut(b"P1 ", b"title"));
        let long = [b'x'; 70];
        for split in [0, 3, 64, 65, 66, 70] {
            let line = cut(&long[..split], &long[split..]);
            assert_eq!(line, format!("{}...", "x".repeat(65)), "split {split}");
        }
        let exact = [b'y'; 68];
        assert_eq!(cut(&exact[..10], &exact[10..]), "y".repeat(68));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2266-2268 (quoted at
    // `write_matrix_line`)
    // The matrix line of a named matrix; without one (RPS-BLAST before a query batch's set-up) the
    // line stops after `Matrix: ` with the stream's failure, which a program reports as NCBI's
    // "BLAST failed to write output" (NCBI BLAST+ 2.17.0: rpsblast and rpstblastn, an empty pipe,
    // outfmt 0: stdout ends in `Matrix: `, exit status 6; rows `db.rps.stdin.empty.0`,
    // `db.rpstbn.stdin.empty.*`).
    #[test]
    fn a_null_matrix_name_stops_the_epilog_after_its_label() {
        let mut named = Vec::new();
        write_matrix_line(&mut named, Some("BLOSUM62")).unwrap();
        assert_eq!(named, b"Matrix: BLOSUM62\n");
        let mut unnamed = Vec::new();
        let error = write_matrix_line(&mut unnamed, None).unwrap_err();
        assert_eq!(unnamed, b"Matrix: ");
        assert!(crate::cli::BadbitInsertion::is(&anyhow::Error::from(error)));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:642-649 (`s_SetFlags`,
    // quoted at `report::query_anchored`) and c++/src/algo/blast/blastinput/blast_args.cpp:2980-2987
    // (quoted at `PairwiseView::report_sort`).
    #[test]
    fn the_views_of_formats_0_to_4_and_their_sorts() {
        use crate::blastinput::app::ReportFormat;
        let sort = crate::report::hit_sort::ReportSort {
            hits: Some(3),
            hsps: Some(1),
        };
        for number in 0..=4 {
            let view = PairwiseView::of(ReportFormat::from_number(number).unwrap()).unwrap();
            assert_eq!(view.outfmt(), number);
            let expected = match number {
                0 => sort,
                _ => crate::report::hit_sort::ReportSort {
                    hits: Some(3),
                    hsps: None,
                },
            };
            assert_eq!(view.report_sort(sort), expected);
        }
        assert_eq!(
            PairwiseView::of(ReportFormat::QueryAnchoredIdentities),
            Some(PairwiseView::QueryAnchored {
                identity: true,
                anchored: true
            })
        );
        assert_eq!(
            PairwiseView::of(ReportFormat::FlatQueryAnchoredNoIdentities),
            Some(PairwiseView::QueryAnchored {
                identity: false,
                anchored: false
            })
        );
        assert_eq!(PairwiseView::of(ReportFormat::Tabular), None);
    }

    fn long_titled_database() -> PairwiseDatabase {
        PairwiseDatabase::new(vec![DatabaseInfo {
            definition: [
                &b"a database title that is much longer than sixty-eight letters, with commas, "[..],
                b"hyphen-ated words and a-very-long-unbroken-token-",
                &[b'z'; 80],
                b" at its end",
            ]
            .concat(),
            date: "Oct 13, 2026  1:18 PM".to_string(),
            num_sequences: 17,
            total_letters: 2287,
            filt_algorithm_name: String::new(),
            filt_algorithm_options: String::new(),
        }])
    }

    // NCBI BLAST+ 2.17.0 `blastp -db longtitle -query hba.faa -num_threads 1` (campaign
    // dbsearch-20261012, scripts/db3_compare.py case shapes.longtitle.0): the prolog's
    // `Database:` lines. The title is wrapped at 68 without the period that the text of
    // a FASTA subject set gets.
    #[test]
    fn database_header_is_the_stored_title_wrapped_and_the_real_counts() {
        let mut written = Vec::new();
        write_blast_database_header(&mut written, &long_titled_database()).unwrap();
        assert_eq!(
            String::from_utf8(written).unwrap(),
            "Database: a database title that is much longer than sixty-eight letters, with\n\
             commas, hyphen-ated words and\n\
             a-very-long-unbroken-token-\n\
             zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz\n\
             zzzzzzzzzzzz at its end\n           17 sequences; 2,287 total letters\n\n"
        );
        let mut written = Vec::new();
        write_blast_database_header(
            &mut written,
            &PairwiseDatabase::new(vec![DatabaseInfo {
                definition: b"Swiss-Prot".to_vec(),
                date: "Oct 12, 2026  1:01 AM".to_string(),
                num_sequences: 575_748,
                total_letters: 209_017_843,
                filt_algorithm_name: String::new(),
                filt_algorithm_options: String::new(),
            }]),
        )
        .unwrap();
        assert_eq!(
            String::from_utf8(written).unwrap(),
            "Database: Swiss-Prot\n           575,748 sequences; 209,017,843 total letters\n\n"
        );
    }

    // NCBI BLAST+ 2.17.0 `blastn -db db/b5\ db/a5 -db_soft_mask 11` (probe m04 of the M4a-2 survey):
    // the prolog shows the first name's algorithm after all titles, the epilog each name's own;
    // `, options:` only for an algorithm with options (probe m06c, `'windowmasker'`).
    #[test]
    fn masked_databases_show_their_algorithms() {
        let database = |definition: &[u8], name: &str, options: &str| DatabaseInfo {
            definition: definition.to_vec(),
            date: "Nov 14, 2026  4:29 AM".to_string(),
            num_sequences: 5,
            total_letters: 7360,
            filt_algorithm_name: name.to_string(),
            filt_algorithm_options: options.to_string(),
        };
        let shown = PairwiseDatabase::new(vec![
            database(b"dust40 only", "dust", "window=32; level=40; linker=1"),
            database(b"mask probe v5", "windowmasker", ""),
        ]);
        let mut written = Vec::new();
        write_blast_database_header(&mut written, &shown).unwrap();
        assert_eq!(
            String::from_utf8(written).unwrap(),
            "Database: dust40 only; mask probe v5\n\
             Masked using: 'dust', options: 'window=32; level=40; linker=1'\n\
             \x20          10 sequences; 14,720 total letters\n\n"
        );
        let mut written = Vec::new();
        write_final_database_block(&mut written, &shown).unwrap();
        assert_eq!(
            String::from_utf8(written).unwrap(),
            "  Database: dust40 only\n\
             \x20 Masked using: 'dust', options: 'window=32; level=40; linker=1'\n\
             \x20   Posted date:  Nov 14, 2026  4:29 AM\n\
             \x20 Number of letters in database: 7,360\n\
             \x20 Number of sequences in database:  5\n\n\
             \x20 Database: mask probe v5\n\
             \x20 Masked using: 'windowmasker'\n\
             \x20   Posted date:  Nov 14, 2026  4:29 AM\n\
             \x20 Number of letters in database: 7,360\n\
             \x20 Number of sequences in database:  5\n\n\n\n"
        );
    }

    // The same report's epilog: the posted date of the database after its title.
    #[test]
    fn database_footer_has_the_posted_date() {
        let mut written = Vec::new();
        write_final_database_block(&mut written, &long_titled_database()).unwrap();
        assert_eq!(
            String::from_utf8(written).unwrap(),
            "  Database: a database title that is much longer than sixty-eight letters, with\n\
             commas, hyphen-ated words and\n\
             a-very-long-unbroken-token-\n\
             zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz\n\
             zzzzzzzzzzzz at its end\n    Posted date:  Oct 13, 2026  1:18 PM\n  \
             Number of letters in database: 2,287\n  \
             Number of sequences in database:  17\n\n\n\n"
        );
    }

    fn heading_of(title: &[u8]) -> String {
        let mut hit = PairwiseHit::from(Hit {
            identity: 100.0,
            length: 10,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end: 10,
            s_start: 1,
            s_end: 10,
            e_value: 1e-5,
            bit_score: 30.0,
            num_ident: 10,
            query_frame: 0,
            query_length: 10,
            q_idx: 0,
            s_idx: 0,
            raw_score: 60,
            sort_query_offset: 0,
            sort_query_end: 10,
            sort_subject_offset: 0,
            sort_subject_end: 10,
            has_sort_offsets: false,
            gap_info: None,
            num_positives: 10,
        });
        hit.subject_length = Some(142);
        let mut written = Vec::new();
        write_database_subject_alignments_part(
            &mut written,
            None,
            DatabaseHeading::Title(title, &crate::report::defline::UNNAMED_SOURCE),
            true,
            &[&hit],
            &[0],
            0..1,
            |writer: &mut Vec<u8>, _| writer.write_all(b"<hsp>\n"),
        )
        .unwrap();
        String::from_utf8(written).unwrap()
    }

    // NCBI BLAST+ 2.17.0 `blastp -db titles -query hba.faa -num_threads 1` (case
    // shapes.titles.0 of the same script): the heading of a database subject is `>` and
    // the stored title with no space between them, wrapped after the first white space
    // from every sixtieth letter; an empty title is the generated one.
    #[test]
    fn database_subject_heading_is_the_stored_title_after_the_mark() {
        assert_eq!(
            heading_of(b"s_plain plain copy"),
            ">s_plain plain copy\nLength=142\n\n<hsp>\n"
        );
        assert_eq!(
            heading_of(
                b"a title that is longer than sixty letters and has spaces so that the heading wraps at a space after the sixtieth"
            ),
            ">a title that is longer than sixty letters and has spaces so that \n\
             the heading wraps at a space after the sixtieth\nLength=142\n\n<hsp>\n"
        );
        assert_eq!(
            heading_of(
                b"a title of exactly sixty letters long 123456789 123456789 1234"
            ),
            ">a title of exactly sixty letters long 123456789 123456789 1234\nLength=142\n\n<hsp>\n"
        );
        let unbroken = b"a_title_without_any_space_that_is_longer_than_sixty_letters_so_that_the_heading_cannot_wrap_at_a_space_anywhere";
        assert_eq!(
            heading_of(unbroken),
            format!(
                ">{}\nLength=142\n\n<hsp>\n",
                String::from_utf8_lossy(unbroken)
            )
        );
        let seventy = [
            &[b'x'; 70][..],
            b" then words after seventy letters without a space before them at all in the title",
        ]
        .concat();
        assert_eq!(
            heading_of(&seventy),
            format!(
                ">{} \nthen words after seventy letters without a space before \n\
                 them at all in the title\nLength=142\n\n<hsp>\n",
                "x".repeat(70)
            )
        );
        // The prefixes that the generated title loses stay in the stored one.
        assert_eq!(
            heading_of(b"TPA: a third party annotation prefix"),
            ">TPA: a third party annotation prefix\nLength=142\n\n<hsp>\n"
        );
        // case shapes.punct00.0: an empty stored title.
        assert_eq!(
            heading_of(b""),
            ">unnamed protein product\nLength=142\n\n<hsp>\n"
        );
    }

    #[test]
    fn k_bits_follows_ctoolkit_compatible() {
        assert_eq!(ncbi_k_bits_for(true), "(bits)");
        assert_eq!(ncbi_k_bits_for(false), "(Bits)");
        assert_eq!(ncbi_k_bits_for(true).len(), ncbi_k_bits_for(false).len());
    }

    #[test]
    fn ka_value_matches_c_percent_hash_8_3g() {
        // Expected strings from C's sprintf(buffer, "%#8.3g ", value) semantics
        // (Python's '%#.3g' % value gives the same strings).
        let cases: [(f64, &str); 32] = [
            (0.0, "0.00"),
            (-1.0, "-1.00"),
            (0.625, "0.625"),
            (0.41, "0.410"),
            (0.78, "0.780"),
            (1.37, "1.37"),
            (1.28, "1.28"),
            (0.46, "0.460"),
            (0.85, "0.850"),
            (0.99996, "1.00"),
            (9.9996, "10.0"),
            (999.5, "1.00e+03"),
            (1.23e-05, "1.23e-05"),
            (0.0001, "0.000100"),
            (9.999e-05, "0.000100"),
            (9.9994e-05, "0.000100"),
            (1000.0, "1.00e+03"),
            (99.95, "100."),
            (99.94, "99.9"),
            (0.125, "0.125"),
            (0.375, "0.375"),
            (12.5, "12.5"),
            (125.0, "125."),
            (1.5, "1.50"),
            (2.25e-07, "2.25e-07"),
            (123456.0, "1.23e+05"),
            (-0.00042, "-0.000420"),
            (5e-324, "4.94e-324"),
            (f64::MAX, "1.80e+308"),
            (0.0009995, "0.000999"),
            (0.00099949, "0.000999"),
            (-999.49, "-999."),
        ];
        for (value, expected) in cases {
            assert_eq!(format_ncbi_ka_value(value), expected, "value {value:e}");
            let mut field = Vec::new();
            write_ncbi_ka_field(&mut field, value).unwrap();
            assert_eq!(String::from_utf8(field).unwrap(), format!("{expected:>8} "));
        }
        assert_eq!(format_ncbi_ka_value(f64::NAN), "nan");
        assert_eq!(format_ncbi_ka_value(f64::INFINITY), "inf");
        assert_eq!(format_ncbi_ka_value(f64::NEG_INFINITY), "-inf");
    }

    fn make_hit() -> Hit {
        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
        // ```c
        // typedef struct BlastHSPList {
        //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
        //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
        //                       Set to 0 if not applicable */
        //    BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
        //    Int4 hspcnt; /**< Number of HSPs saved */
        //    ...
        // } BlastHSPList;
        // ```
        Hit {
            identity: 95.0,
            length: 100,
            mismatch: 5,
            gapopen: 0,
            q_start: 1,
            q_end: 100,
            s_start: 1,
            s_end: 100,
            e_value: 1e-50,
            bit_score: 185.5,
            num_ident: 95,
            // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
            // ```c
            // if (hsp->query.frame != hsp->subject.frame) {
            //    *q_end = query_length - hsp->query.offset;
            //    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
            // }
            // ```
            query_frame: 1,
            query_length: 0,
            q_idx: 0,
            s_idx: 0,
            raw_score: 200,
            // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-143
            // BlastHSP stores query/subject BlastSeg offsets used by HSP comparators.
            sort_query_offset: 0,
            sort_query_end: 0,
            sort_subject_offset: 0,
            sort_subject_end: 0,
            has_sort_offsets: false,
            gap_info: None,
            num_positives: 95,
        }
    }

    // NCBI reference: c++/src/corelib/ncbistr.cpp:5088-5340
    // enum EScore { eForced, ePunct, eComma, eSpace, eNewline };
    #[test]
    fn flatfile_wrapping_preserves_ncbi_break_boundaries() {
        for (input, width, expected) in [
            ("abcd", 4, "abcd\n"),
            ("abc-def", 4, "abc-\ndef\n"),
            ("a b c", 4, "a b\nc\n"),
            ("abcd ef", 4, "abcd\nef\n"),
            ("a\n\nb", 4, "a\n\nb\n"),
            ("abc\u{8}d", 4, "abd\n"),
            ("ab\u{8}cd\nz", 68, "ab\u{8}cd\nz\n"),
            ("abc,def", 4, "abc,\ndef\n"),
        ] {
            let mut output = Vec::new();
            write_flatfile_wrapped(&mut output, input.as_bytes(), width).unwrap();
            assert_eq!(output, expected.as_bytes(), "{input:?}");
        }
    }

    #[test]
    fn test_write_subject_header() {
        let mut output = Vec::new();
        write_subject_header(&mut output, b"seq1", Some(&b"Test sequence"[..]), Some(500)).unwrap();
        let output_str = String::from_utf8(output).unwrap();

        assert!(output_str.contains("> seq1 Test sequence"));
        assert!(output_str.contains("Length=500"));
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1598-1626 (quoted at
    // `write_sequence_row`): each row of `m_LineLen` columns; rpsblast's `-line_length` (1 or
    // more) sets the width. A width of 1 prints one column per row with its own coordinates,
    // the coordinates' column as wide as the HSP's widest coordinate.
    #[test]
    fn alignment_rows_of_one_column_step_one_letter() {
        let hit = Hit {
            length: 5,
            q_start: 8,
            q_end: 12,
            s_start: 98,
            s_end: 102,
            ..make_hit()
        };
        let row_of = |line_length: usize| {
            let config = PairwiseConfig {
                line_length,
                program: "rpsblast".to_string(),
                show_frame: false,
                ..PairwiseConfig::default()
            };
            let mut output = Vec::new();
            write_alignment_with_sequences(&mut output, &hit, "ACDeF", "ACDEW", &[], &config)
                .unwrap();
            String::from_utf8(output).unwrap()
        };
        let one = row_of(1);
        let queries: Vec<&str> = one.lines().filter(|l| l.starts_with("Query")).collect();
        assert_eq!(
            queries,
            [
                "Query  8    A  8",
                "Query  9    C  9",
                "Query  10   D  10",
                "Query  11   e  11",
                "Query  12   F  12"
            ]
        );
        let subjects: Vec<&str> = one.lines().filter(|l| l.starts_with("Sbjct")).collect();
        assert_eq!(subjects[0], "Sbjct  98   A  98");
        assert_eq!(subjects[4], "Sbjct  102  W  102");
        // The middle line reads the letters without case (the lower-case mask is a display).
        assert_eq!(one.lines().nth(13).unwrap(), "            E");
        // Widths at and above the HSP's length give the one row of the default.
        assert_eq!(row_of(5), row_of(60));
        assert_eq!(row_of(1000), row_of(60));
        // A first row of one residue at 1 shows no coordinates (`j == 0 && start == 1 && end ==
        // 1`), as NCBI's `-line_length 1` rows of a query that aligns from its first letter.
        let hit = Hit {
            length: 2,
            q_start: 1,
            q_end: 2,
            s_start: 1,
            s_end: 2,
            ..make_hit()
        };
        let config = PairwiseConfig {
            line_length: 1,
            program: "rpsblast".to_string(),
            show_frame: false,
            ..PairwiseConfig::default()
        };
        let mut output = Vec::new();
        write_alignment_with_sequences(&mut output, &hit, "MQ", "MQ", &[], &config).unwrap();
        assert_eq!(
            String::from_utf8(output).unwrap(),
            "Query     M  \n          M\nSbjct     M  \n\nQuery  2  Q  2\n          Q\nSbjct  2  Q  2\n\n"
        );
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1598-1626 (quoted at
    // `write_sequence_row`): the rule `j == 0 && start == 1 && end == 1` hides the coordinates of
    // a first row of one residue at 1 in BLASTX's subject and TBLASTN's query, the rows of letters
    // (NCBI 2.17.0 `-line_length 1`: "Sbjct       M  ", "Query        M  "); a translated row's
    // first residue spans the nucleotides 1 to 3, so it keeps them, in TBLASTX's rows too
    // ("Query  1    F  3").
    #[test]
    fn translated_alignment_rows_of_one_column_hide_only_a_letter_row_at_1() {
        let config = PairwiseConfig {
            line_length: 1,
            ..PairwiseConfig::default()
        };
        let hit = |q_end: usize, s_end: usize, query_frame, subject_frame| {
            let mut hit = PairwiseHit::from(Hit {
                length: 2,
                q_start: 1,
                q_end,
                s_start: 1,
                s_end,
                gap_info: None,
                ..make_hit()
            });
            hit.query_seq = Some("MR".to_string());
            hit.subject_seq = Some("MR".to_string());
            hit.query_frame = query_frame;
            hit.subject_frame = subject_frame;
            hit
        };
        let text = |write: &dyn Fn(&mut Vec<u8>) -> io::Result<()>| {
            let mut output = Vec::new();
            write(&mut output).unwrap();
            String::from_utf8(output).unwrap()
        };
        let blastx = hit(6, 2, Some(1), None);
        assert_eq!(
            text(&|out| write_blastx_alignment(out, &blastx, &config)),
            "Query  1  M  3\n          M\nSbjct     M  \n\n\
             Query  4  R  6\n          R\nSbjct  2  R  2\n\n"
        );
        let tblastn = hit(2, 6, None, Some(1));
        assert_eq!(
            text(&|out| write_tblastn_alignment(out, &tblastn, &config)),
            "Query     M  \n          M\nSbjct  1  M  3\n\n\
             Query  2  R  2\n          R\nSbjct  4  R  6\n\n"
        );
        let tblastx = hit(6, 6, Some(1), Some(1));
        assert_eq!(
            text(&|out| write_tblastx_alignment(out, &tblastx, 6, &config)),
            "Query  1  M  3\n          M\nSbjct  1  M  3\n\n\
             Query  4  R  6\n          R\nSbjct  4  R  6\n\n"
        );
        // On a minus frame the first residue of a translated row at 1 to 3 runs from 3 to 1.
        // The query record's length is the HSP's (frame -1 ends at its last nucleotide).
        let mut minus = hit(1, 2, Some(-1), None);
        (minus.hit.q_start, minus.hit.query_length) = (6, 6);
        assert_eq!(
            text(&|out| write_blastx_alignment(out, &minus, &config)),
            "Query  6  M  4\n          M\nSbjct     M  \n\n\
             Query  3  R  1\n          R\nSbjct  2  R  2\n\n"
        );
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1598-1599
    // ```c++
    //     int start = alnRoInfo->seqStarts[row].front() + 1;  //+1 for 1 based
    //     int end = alnRoInfo->seqStops[row].front() + 1;
    // ```
    // NCBI BLAST+ 2.17.0 (line RG T1, cases d3.*): with a wrapped screen only the end of a
    // translated row moves, to the end of its first screen; the letters are the whole row's
    // and the rows of letters do not change. blastx (minus frame) `Query  2707 ... 2706` at
    // 1431655766, `... 2703` at 1431655767, the row's `... 1325` at 1431655765; tblastn's
    // subject (plus frame) 574..575, 574..578 and 574..1431; tblastx's query 2..3, 2..6 and
    // 2..598 and minus-frame rows 598..597 and 598..594.
    #[test]
    fn a_translated_row_ends_where_ncbi_s_wrapped_screen_ends() {
        let config = |line_length| PairwiseConfig {
            line_length,
            ..PairwiseConfig::default()
        };
        let hit = |q: (usize, usize), s: (usize, usize), query_frame, subject_frame| {
            let mut hit = PairwiseHit::from(Hit {
                length: 5,
                q_start: q.0,
                q_end: q.1,
                s_start: s.0,
                s_end: s.1,
                gap_info: None,
                ..make_hit()
            });
            hit.query_seq = Some("MRKVL".to_string());
            hit.subject_seq = Some("MRKVL".to_string());
            hit.query_frame = query_frame;
            hit.subject_frame = subject_frame;
            hit
        };
        let text = |write: &dyn Fn(&mut Vec<u8>) -> io::Result<()>| {
            let mut output = Vec::new();
            write(&mut output).unwrap();
            String::from_utf8(output).unwrap()
        };
        // The records' lengths match the frames: the BLASTX query's of frame -1 at 2707, the
        // TBLASTX subject's of frame -1 at 598 (and the TBLASTX query's frame at 2 is +2).
        let mut blastx = hit((2707, 2693), (4, 8), Some(-1), None);
        blastx.hit.query_length = 2707;
        for (line_length, end) in [(1431655765, 2693), (1431655766, 2706), (1431655767, 2703)] {
            assert_eq!(
                text(&|out| write_blastx_alignment(out, &blastx, &config(line_length))),
                format!("Query  2707  MRKVL  {end}\n             MRKVL\nSbjct  4     MRKVL  8\n\n"),
                "-line_length {line_length}"
            );
        }
        let tblastn = hit((1, 5), (574, 588), None, Some(1));
        for (line_length, end) in [(1431655765, 588), (1431655766, 575), (1431655767, 578)] {
            assert_eq!(
                text(&|out| write_tblastn_alignment(out, &tblastn, &config(line_length))),
                format!("Query  1    MRKVL  5\n            MRKVL\nSbjct  574  MRKVL  {end}\n\n"),
                "-line_length {line_length}"
            );
        }
        let mut tblastx = hit((2, 16), (598, 584), Some(2), Some(-1));
        tblastx.subject_length = Some(598);
        for (line_length, ends) in [
            (1431655765, (16, 584)),
            (1431655766, (3, 597)),
            (1431655767, (6, 594)),
            (2147483647, (16, 584)),
        ] {
            assert_eq!(
                text(&|out| write_tblastx_alignment(out, &tblastx, 16, &config(line_length))),
                format!(
                    "Query  2    MRKVL  {}\n            MRKVL\nSbjct  598  MRKVL  {}\n\n",
                    ends.0, ends.1
                ),
                "-line_length {line_length}"
            );
        }
        // Rows of letters keep their ends (BLASTP's rows).
        let blastp = hit((1, 5), (11, 15), None, None);
        let blastp_text = |line_length| {
            text(&|out| {
                write_alignment_with_sequences(
                    out,
                    &blastp.hit,
                    "MRKVL",
                    "MRKVL",
                    &[],
                    &PairwiseConfig {
                        program: "blastp".into(),
                        ..config(line_length)
                    },
                )
            })
        };
        assert_eq!(blastp_text(1431655766), blastp_text(60));
        assert!(blastp_text(1431655766).starts_with("Query  1   MRKVL  5\n"));
    }

    // `write_buffered` writes every letter once and in order across its buffer's ends, `char`s of
    // two to four bytes included.
    #[test]
    fn buffered_letters_are_written_in_order_across_buffers() {
        let chars: String = (0..1000).map(|i| ['A', 'é', '€', '𝄞'][i % 7 % 4]).collect();
        let mut out = Vec::new();
        write_buffered(&mut out, chars.chars()).unwrap();
        assert_eq!(out, chars.as_bytes());
        let bytes: Vec<u8> = (0..1000u32).map(|i| (i % 251) as u8).collect();
        let mut out = Vec::new();
        write_buffered(&mut out, bytes.iter().copied()).unwrap();
        assert_eq!(out, bytes);
    }

    // The parallel writer's estimate of an HSP's text (`x_pairwise_par::hit_estimate`) holds at
    // every line length: the rows of a 300-column HSP whose coordinates have 4 digits (55 bytes
    // besides the letters per line) and about 300 bytes of its header lines fit in it, and it is
    // at most a tenth more than the rows and its 400 bytes of header lines.
    #[test]
    fn the_writer_s_estimate_of_an_hsp_holds_at_every_line_length() {
        let hit = Hit {
            length: 300,
            q_start: 1001,
            q_end: 1300,
            s_start: 2001,
            s_end: 2300,
            gap_info: None,
            ..make_hit()
        };
        let letters = "ACDEFGHIKLMNPQRSTVWY".repeat(15);
        for line_length in [1, 2, 3, 7, 59, 60, 61, 299, 300, 1000, 2147483647] {
            let mut rows = Vec::new();
            write_alignment_with_sequences(
                &mut rows,
                &hit,
                &letters,
                &letters,
                &[],
                &PairwiseConfig {
                    line_length,
                    program: "blastp".into(),
                    ..PairwiseConfig::default()
                },
            )
            .unwrap();
            let estimate = crate::report::x_pairwise_par::hit_estimate(300, line_length);
            assert!(
                rows.len() + 300 <= estimate && estimate <= rows.len() * 11 / 10 + 400,
                "-line_length {line_length}: rows {} B, estimate {estimate} B",
                rows.len()
            );
        }
    }

    /// NCBI's `GetWholeAlnSeqString` (alnvec.cpp:199-397, quoted at
    /// `query_anchored::aln_vec::AlnMap::aln_seq_string`), a second transcription, for one row of
    /// `segments` (a 0-based Dense-seg start, -1 for a gap, and a length in columns), as
    /// written, unsigned arithmetic included: the first entry of `scrn_rights`.
    fn ncbi_first_screen_right(
        segments: &[(i64, u32)],
        plus: bool,
        width: u32,
        scrn_width: u32,
    ) -> i64 {
        let (mut aln_pos, mut scrn_pos) = (0u32, 0u32);
        let (mut lft, mut rgt) = (-1i64, -1i64);
        let mut rights = Vec::new();
        let left_seg = segments.iter().position(|seg| seg.0 >= 0).unwrap();
        let toward = if plus { -1 } else { 1 };
        for (seg, &(start, seg_len)) in segments.iter().enumerate() {
            let len = seg_len.wrapping_mul(width);
            if start >= 0 {
                let stop = start + i64::from(len) - 1;
                let first = if plus { start } else { stop };
                if lft < 0 {
                    lft = first;
                    if rgt < 0 {
                        rgt = lft;
                    }
                }
                let nscrns = (aln_pos.wrapping_sub(scrn_pos) / scrn_width) as i32;
                for _ in 0..nscrns {
                    rights.push(rgt);
                    scrn_pos = scrn_pos.wrapping_add(scrn_width);
                }
                if nscrns > 0 {
                    lft = first;
                }
                let nscrns = (aln_pos.wrapping_add(len).wrapping_sub(scrn_pos) / scrn_width) as i32;
                let mut curr_pos = aln_pos;
                for _ in 0..nscrns {
                    let delta = if plus {
                        scrn_width.wrapping_sub(curr_pos.wrapping_sub(scrn_pos))
                    } else {
                        curr_pos.wrapping_sub(scrn_pos).wrapping_sub(scrn_width)
                    } as i32 as i64;
                    if if plus { lft < start } else { lft > stop } {
                        lft = first + delta;
                        rgt = lft + toward;
                    } else {
                        rgt = lft + toward + delta;
                        lft += delta;
                    }
                    if seg == left_seg && lft == rgt {
                        rgt += toward;
                    }
                    rights.push(rgt);
                    scrn_pos = scrn_pos.wrapping_add(scrn_width);
                    curr_pos = scrn_pos;
                }
                if aln_pos.wrapping_add(len) <= scrn_pos {
                    lft = -1;
                }
                rgt = if plus { stop } else { start };
            }
            aln_pos = aln_pos.wrapping_add(len);
        }
        let pos_diff = aln_pos.wrapping_sub(scrn_pos);
        if pos_diff > 0 {
            let nscrns = pos_diff / scrn_width + u32::from(pos_diff % scrn_width != 0);
            for _ in 0..nscrns {
                rights.push(rgt);
            }
        }
        rights[0]
    }

    // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:310-358 (quoted at
    // `query_anchored::aln_vec::AlnMap::aln_seq_string`): subject rows of residues and gaps, on
    // both strands, of width 1 or 3, with screens narrower than the row (wrapped ones: 3 * L -
    // 2^32 nucleotides, 2, 5, 8, ...) and wider, and residue runs split into segments at random
    // (the other row's gaps split them): the end of the first line that `write_hsp_rows` prints
    // (`query_anchored::HspAlnVec`, the owner of the rows' coordinates) is the first screen's
    // right of a second transcription of NCBI's computation (`ncbi_first_screen_right`).
    #[test]
    fn the_first_line_ends_at_ncbi_s_first_screen() {
        use crate::report::serial::seqalign::Strand;
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut random = |bound: u64| {
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            state.wrapping_mul(0x2545_f491_4f6c_dd1d) % bound
        };
        let mut rows = crate::report::query_anchored::HspAlnVec::default();
        for case in 0..3000 {
            let letters: u32 = if case % 5 == 0 { 1 } else { 3 };
            let columns = 1 + random(60) as usize;
            // Residues, with runs of gaps after the first one (an HSP's row begins with one).
            let mut row: Vec<u8> = Vec::with_capacity(columns);
            while row.len() < columns {
                if row.is_empty() || random(4) != 0 {
                    row.push(b'A');
                } else {
                    let gaps = 1 + random(4) as usize;
                    row.resize(row.len() + gaps, b'-');
                }
            }
            row.truncate(columns);
            let plus = random(2) == 0;
            let start = 1 + random(5000) as i64;
            let start = if plus {
                start
            } else {
                start + i64::from(letters) * columns as i64
            };
            let span = u64::from(letters) * columns as u64;
            let line_length = match (letters, random(3)) {
                // A wrapped screen narrower than the row: (W + 2^32) / 3 for a W of 2 mod 3.
                (3, 0 | 1) if span > 2 => {
                    let screen = 2 + 3 * random((span - 2).div_ceil(3));
                    ((screen + (1u64 << 32)) / 3) as usize
                }
                (3, 2) => [1431655765, 2147483647][random(2) as usize],
                _ => columns + random(100) as usize,
            };
            // The row's segments, runs of residues split at random, under a query row of
            // residues (width 1, plus strand).
            let mut segments = Vec::new();
            let mut next = start - 1;
            let mut column = 0;
            while column < row.len() {
                let gap = row[column] == b'-';
                let mut run = 1;
                while column + run < row.len()
                    && (row[column + run] == b'-') == gap
                    && (gap || random(3) != 0)
                {
                    run += 1;
                }
                let nucleotides = i64::from(letters) * run as i64;
                if gap {
                    segments.push((-1, run as u32));
                } else if plus {
                    segments.push((next, run as u32));
                    next += nucleotides;
                } else {
                    segments.push((next - nucleotides + 1, run as u32));
                    next -= nucleotides;
                }
                column += run;
            }
            let (mut starts, mut lens, mut strands) = (Vec::new(), Vec::new(), Vec::new());
            let mut query = 0;
            for &(subject, len) in &segments {
                starts.extend([query, subject as i32]);
                lens.push(len as i32);
                strands.extend([
                    Strand::Plus,
                    if plus { Strand::Plus } else { Strand::Minus },
                ]);
                query += len as i32;
            }
            rows.set(&starts, &lens, &strands, [1, letters], line_length)
                .unwrap();
            let last = line_length.min(columns) as i64 - 1;
            let (_, right) = rows.line(1, 0, 0, last).unwrap().unwrap();
            let scrn_width = (line_length as u32).wrapping_mul(letters);
            assert_eq!(
                i64::from(right) + 1,
                ncbi_first_screen_right(&segments, plus, letters, scrn_width) + 1,
                "plus={plus} letters={letters} L={line_length} segments={segments:?}"
            );
        }
    }

    #[test]
    fn test_write_hsp_info() {
        let hit = PairwiseHit::from(make_hit());
        let config = PairwiseConfig::default();

        let mut output = Vec::new();
        write_hsp_info(&mut output, &hit, &config).unwrap();
        let output_str = String::from_utf8(output).unwrap();

        assert!(output_str.contains("Score ="));
        assert!(output_str.contains("Expect ="));
        assert!(output_str.contains("Identities ="));
    }

    // NCBI reference: c++/src/objtools/align_format/showalign.cpp:4017-4018
    // ```c++
    //                               m_AV->StrandSign(0),
    //                               m_AV->StrandSign(1),
    // ```
    // A BLASTN HSP of one letter has the same start and end; its strand is its frame's.
    #[test]
    fn blastn_strand_comes_from_the_hsp_frame() {
        let config = PairwiseConfig {
            program: "blastn".to_string(),
            show_frame: false,
            ..PairwiseConfig::default()
        };
        for (query_frame, strand) in [(1, " Strand=Plus/Plus\n"), (-1, " Strand=Plus/Minus\n")] {
            let mut hit = PairwiseHit::from(Hit {
                length: 1,
                s_start: 10,
                s_end: 10,
                query_frame,
                ..make_hit()
            });
            hit.query_seq = Some("G".to_string());
            hit.subject_seq = Some("G".to_string());
            let mut output = Vec::new();
            write_hsp_info(&mut output, &hit, &config).unwrap();
            let output = String::from_utf8(output).unwrap();
            assert!(output.contains(strand), "{output}");
        }
    }

    // Line DBID item ID4f: BLASTN's outfmt 0 makes its rows a window of queries at a time and
    // writes each window's report (`algorithm::blastn::blast_engine::run`); the reports of
    // consecutive query ranges, each given only its own hits, are the report of all queries,
    // on one thread and on a pool.
    #[test]
    #[cfg_attr(
        not(feature = "parallel"),
        ignore = "needs threads: a build without the parallel feature has none"
    )]
    fn blastn_reports_of_consecutive_query_ranges_are_the_whole_report() {
        let config = PairwiseConfig {
            program: "blastn".to_string(),
            show_frame: false,
            ..PairwiseConfig::default()
        };
        let karlin = KarlinParams {
            lambda: 1.28,
            k: 0.46,
            h: 0.85,
            alpha: 0.0,
            beta: 0.0,
        };
        let queries: Vec<BlastnPairwiseQuery> = (0..4)
            .map(|q_idx| BlastnPairwiseQuery {
                query_name: format!("q{q_idx} query {q_idx}").into_bytes(),
                query_length: 100,
                karlin: Some((karlin, karlin)),
                effective_search_space: 10_000,
            })
            .collect();
        let report = BlastnPairwiseReport {
            version: "2.17.0+".to_string(),
            megablast: true,
            database_name: "subjects".to_string(),
            database_num_sequences: 3,
            database_total_letters: 300,
            database: None,
            reward: 1,
            penalty: -2,
            gap_open: 0,
            gap_extend: 0,
            zero_gap_extension_formula: true,
            window_size: 0,
            gapped: true,
            num_descriptions: 500,
            num_alignments: 250,
            unsearched: vec![false; queries.len()],
            epilog: true,
            prolog: true,
            parsed_subjects: Vec::new(),
            subject_sources: Vec::new(),
        };
        // Query 0: subjects 1 (2 HSPs) and 2; query 1: no hits; queries 2 and 3: subject 1.
        let hits: Vec<PairwiseHit> = [(0, 1), (0, 1), (0, 2), (2, 1), (3, 1)]
            .into_iter()
            .map(|(q_idx, s_idx)| {
                let mut hit = PairwiseHit::from(Hit {
                    q_idx,
                    s_idx,
                    ..make_hit()
                });
                hit.query_seq = Some("A".repeat(100));
                hit.subject_seq = Some("A".repeat(100));
                hit.gaps = Some(0);
                hit.subject_length = Some(100);
                hit
            })
            .collect();
        let titles: Vec<Arc<[u8]>> = (0..3)
            .map(|s_idx| Arc::from(format!("s{s_idx} subject {s_idx}").as_bytes()))
            .collect();
        crate::utils::threading::with_search_pool(4, "test", |pool| {
            for x_pool in [None, Some(pool)] {
                let write = |range: std::ops::Range<usize>| {
                    let range_hits: Vec<PairwiseHit> = hits
                        .iter()
                        .filter(|hit| range.contains(&(hit.hit.q_idx as usize)))
                        .cloned()
                        .collect();
                    let mut written = Vec::new();
                    write_blastn_pairwise_report(
                        &range_hits,
                        &mut written,
                        &config,
                        &queries,
                        range,
                        &titles,
                        &report,
                        None,
                        None,
                        x_pool,
                        None,
                        &crate::report::db_subject::HitLoads::default(),
                    )
                    .unwrap();
                    written
                };
                let whole = write(0..queries.len());
                let text = String::from_utf8(whole.clone()).unwrap();
                assert!(text.starts_with("BLASTN 2.17.0+"), "{text}");
                assert_eq!(text.matches("Query= ").count(), queries.len(), "{text}");
                assert!(text.contains("No hits found"), "{text}");
                for cuts in [vec![1], vec![2], vec![3], vec![1, 2, 3], vec![1, 3]] {
                    let mut joined = Vec::new();
                    let mut start = 0;
                    for end in cuts.into_iter().chain([queries.len()]) {
                        joined.extend(write(start..end));
                        start = end;
                    }
                    assert_eq!(String::from_utf8(joined).unwrap(), text);
                }
            }
            Ok(())
        })
        .unwrap();
    }
}

// NCBI c++/src/algo/blast/format/blast_format.cpp:348-445,1490-1589:
// PrintProlog(); AcknowledgeBlastQuery(...); x_DisplayDeflines(...);
// display.DisplaySeqalign(...); x_PrintOneQueryFooter(...);
// This uses the local-subject database header and TBLASTN translated alignment.
/// The outfmt 0 prolog of tblastn (the version, the references and the database), which
/// NCBI writes before it reads the first query batch (as `write_tblastx_pairwise_prolog`).
/// `program` is the program line's name: `TBLASTN`, or `PSITBLASTN` with `-in_pssm`
/// (`Blast_ProgramNameFromType`, blast_format.cpp:94).
pub fn write_tblastn_pairwise_prolog<W: Write>(
    writer: &mut W,
    program: &str,
    version: &str,
    database_name: &str,
    database_num_sequences: usize,
    database_total_letters: usize,
    database: Option<&PairwiseDatabase>,
) -> io::Result<()> {
    write_translated_pairwise_intro(writer, program, version)?;
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
pub fn write_tblastn_pairwise_report<W: Write>(
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
    warnings: Option<&mut super::query_warnings::QueryWarnings<'_>>,
    epilog: bool,
    anchored: Option<&dyn super::query_anchored::ReportSequences>,
    loads: &crate::report::db_subject::HitLoads,
) -> io::Result<()> {
    // NCBI c++/src/algo/blast/api/local_blast.cpp:177-224:
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
    // NCBI c++/src/algo/blast/format/blast_format.cpp:372-424:
    // BlastPrintVersionInfo(m_Program,...); BlastPrintReference(...);
    write_translated_pairwise_intro(&mut writer, program, &report.version)?;
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

    // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1411
    // ```c++
    // CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
    // ```
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
            loads,
        },
        0,
        &hits_by_query,
        subject_titles,
        probe,
        warnings,
        None,
    )?;
    // NCBI reference: c++/src/app/blast/tblastn_app.cpp:342
    // ```c++
    //         formatter.PrintEpilog(opt);
    // ```
    // A search that stops at a query batch (`Empty CBlastQueryVector` with `-query_loc`)
    // writes no epilog.
    if !epilog {
        return writer.flush();
    }
    // The blank lines before the epilog (blast_format.cpp:2249).
    writeln!(&mut writer)?;
    writeln!(&mut writer)?;
    write_blastp_final_footer(&mut writer, report, database)?;
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
    /// What the reports' loads do at each subject, by subject (`thrown_load`).
    loads: &'a crate::report::db_subject::HitLoads,
}

impl super::x_pairwise_par::QuerySections for TblastnSections<'_> {
    fn num_alignments(&self) -> usize {
        self.report.num_alignments
    }

    fn line_length(&self) -> usize {
        self.config.line_length
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1522-1530
    // ```c++
    //     if(m_HitsSortOption > 0 || m_HspsSortOption > 0) {
    //         aln_set = CBlastFormatUtil::SortSeqalignForSortableFormat(
    // ```
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

    fn head<W: Write>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: Option<&super::x_pairwise_par::QueryGroup<'_>>,
        read: super::query_warnings::SubjectReadHook<'_, W>,
    ) -> io::Result<()> {
        let query = &self.queries[q_idx];
        // NCBI blast_format.cpp:1491: m_Outfile << "\n\n";
        writeln!(writer)?;
        writeln!(writer)?;
        write_blastp_query_header(writer, &self.query_titles[q_idx], query.query_length)?;
        let Some(group) = group else {
            write_no_hits_found(writer)?;
            // NCBI blast_results.cpp:82-97 leaves both Karlin blocks null
            // when there is no valid context. blast_format.cpp:445-477
            // writes only the blank lines and zero search space in that case.
            // NCBI c++/src/algo/blast/api/local_blast.cpp:177-180,204-208:
            // if (m_PrelimSearch->CheckInternalData() != 0)
            //     new CBlastAncillaryData(tmp_pair, tmp_pair, tmp_pair, 0);
            // The all-invalid batch uses -1 sentinel blocks. An invalid query
            // within a searched batch has null blocks (blast_results.cpp:82-103).
            if self.query_batch_skipped[q_idx] {
                write_tblastn_unsearched_query_footer_spacing(writer, false)?;
            } else if !self.query_validity[q_idx] {
                writeln!(writer)?;
                writeln!(writer)?;
                writeln!(writer)?;
                writeln!(writer, "Effective search space used: 0")?;
            } else {
                write_blastp_query_footer_spacing(
                    writer,
                    query.ungapped_karlin,
                    self.report.gapped_karlin,
                    self.report.gumbel,
                    query.effective_search_space,
                    false,
                )?;
            }
            return Ok(());
        };
        // NCBI reference: c++/src/objtools/align_format/showdefline.cpp:1004-1011
        // ```c++
        // void CShowBlastDefline::DisplayBlastDefline(CNcbiOstream & out)
        // {
        //     x_InitDeflineTable();
        //     ...
        //     x_DisplayDefline(out);
        // ```
        // The table of every program follows `x_InitDeflineTable`
        // (`write_blastn_description_table`): the subject's highest bit score and that HSP's
        // E-value, and the title of the subject (`CDeflineGenerator`), for the first
        // `m_NumDescriptions` subjects. The titles are those of the records read (by
        // subject index), as BLASTN's and TBLASTX's.
        write_blastn_description_table(
            writer,
            group.described(self.report.num_descriptions),
            group.subject_order.len() <= self.report.num_descriptions,
            &group.subject_hits,
            &group.subject_titles,
            self.parsed_subjects,
            self.subject_sources,
            false,
            false,
            self.config.sort_columns(query.query_length),
            read,
            self.loads,
        )?;
        writeln!(writer)
    }

    // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:1014-1040
    // (PruneSeqalign keeps the alignments of the first `m_NumAlignments` subjects).
    fn subject<W: Write>(
        &self,
        writer: &mut W,
        probe: Option<&mut FormatProbe<'_>>,
        group: &super::x_pairwise_par::QueryGroup<'_>,
        s_idx: u32,
        part: std::ops::Range<usize>,
    ) -> io::Result<()> {
        let config = self.config;
        // NCBI c++/src/objtools/align_format/showalign.cpp:3650-3668:
        // display leaves an extra blank line after each translated HSP.
        let write_hsp = |writer: &mut W, hit: &PairwiseHit| {
            write_tblastn_hsp_info(writer, hit)?;
            write_tblastn_alignment(writer, hit, config)?;
            writeln!(writer)
        };
        // A database sequence's heading is its stored title after `>` (as BLASTP's
        // and BLASTN's `-db` reports), or the lines of the definition lines of one with
        // parsed Seq-ids.
        let parsed = parsed_subject(self.parsed_subjects, s_idx);
        if self.database.is_some() || parsed.is_some() {
            let heading = match parsed {
                Some(subject) => DatabaseHeading::loaded(
                    subject,
                    true,
                    group.own_hsps(s_idx),
                    self.thrown(s_idx).is_some(),
                ),
                None => DatabaseHeading::Title(
                    &group.subject_titles[&s_idx],
                    stored_source(self.subject_sources, s_idx),
                ),
            };
            return write_database_subject_alignments_part(
                writer,
                probe,
                heading,
                false,
                &group.subject_hits[&s_idx],
                &group.subject_hit_indices[&s_idx],
                part,
                write_hsp,
            );
        }
        write_subject_alignments_part(
            writer,
            probe,
            s_idx,
            self.parsed_subjects,
            &group.subject_titles[&s_idx],
            false,
            &group.subject_hits[&s_idx],
            &group.subject_hit_indices[&s_idx],
            part,
            write_hsp,
        )
    }

    fn thrown(&self, s_idx: u32) -> Option<(&crate::ncbi_exception::NcbiException, bool)> {
        thrown_load(self.loads, self.parsed_subjects, s_idx)
    }

    fn view(&self) -> PairwiseView {
        self.config.view
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1550-1552
    // ```c++
    //     CSeq_align_set copy_aln_set;
    //     CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_NumAlignments);
    // ```
    /// The query-anchored section of the query's aligned subjects: Dense-segs of a protein query
    /// row and a translated nucleotide subject row (`seqalign::translated_row_segments`).
    fn alignments<W: Write>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: &super::x_pairwise_par::QueryGroup<'_>,
    ) -> io::Result<()> {
        write_anchored_section(writer, self.config, self.anchored, q_idx, |_| {
            section_hsps(
                super::query_anchored::Hsps::of_translated_subjects(),
                self.report.num_alignments,
                group,
                |_| false,
                tblastn_row_segments,
            )
        })
    }

    fn foot<W: Write>(&self, writer: &mut W, q_idx: usize) -> io::Result<()> {
        let query = &self.queries[q_idx];
        write_blastp_query_footer_spacing(
            writer,
            query.ungapped_karlin,
            self.report.gapped_karlin,
            self.report.gumbel,
            query.effective_search_space,
            false,
        )
    }
}

// NCBI c++/src/objtools/align_format/showalign.cpp:3578-3604,2122-2149:
// out << " Score = " << bit_score << " bits (" << score << ")";
// if (comp_adjust) out << ", Method: Compositional matrix adjust.";
// out << " Identities = ..." << " Positives = ..." << " Gaps = ...";
// TBLASTN reports one translated subject frame, without a query frame.
fn write_tblastn_hsp_info<W: Write>(writer: &mut W, hit: &PairwiseHit) -> io::Result<()> {
    let h = &hit.hit;
    let bits = format_bitscore_ncbi(h.bit_score);
    let evalue = format_evalue_ncbi(h.e_value);
    // NCBI c++/src/objtools/align_format/showalign.cpp:3599-3604:
    // if (comp_adj_method == 1) out << ", Method: Composition-based stats.";
    // else if (comp_adj_method == 2) out << ", Method: Compositional matrix adjust.";
    let method = match hit.comp_adjust_method {
        Some(1) => ", Method: Composition-based stats.",
        Some(2) => ", Method: Compositional matrix adjust.",
        _ => "",
    };
    // NCBI c++/src/algo/blast/api/blast_seqalign.cpp:1180-1182:
    // if (hsp->num > 1) scores.push_back(s_MakeScore("sum_n", ..., hsp->num));
    // NCBI c++/src/objtools/align_format/showalign.cpp:3595-3598:
    // if (aln_vec_info->sum_n > 0) out << "(" << aln_vec_info->sum_n << ")";
    let expect = expect_label(hit.sum_n);
    writeln!(
        writer,
        " Score = {} bits ({}),  {} = {}{}",
        bits, h.raw_score, expect, evalue, method
    )?;
    let length = h.length;
    // NCBI align_format_util.cpp:3152-3161:
    // if (numerator == denominator) return 100;
    // int retval = (int)(0.5 + 100.0*numerator/denominator);
    // return min(99, retval);
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
    // ```c++
    // int CAlignFormatUtil::GetPercentMatch(int numerator, int denominator)
    // {
    //      if (numerator == denominator)
    //         return 100;
    //      else {
    //        int retval =(int) (0.5 + 100.0*((double)numerator)/((double)denominator));
    //        retval = min(99, retval);
    //        return retval;
    //      }
    // }
    // ```
    let percent = |n: usize| ncbi_percent_match(n, length);
    let positives = hit.positives.unwrap_or(h.num_positives);
    let gaps = hit.gaps.unwrap_or_else(|| h.gap_letters());
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
    // ```c++
    //     if(master_frame != 0 && slave_frame != 0) {
    //         out <<" Frame = " << ((master_frame > 0) ? "+" : "")
    //             << master_frame <<"/"<<((slave_frame > 0) ? "+" : "")
    //             << slave_frame<<"\n";
    //     } else if (master_frame != 0){
    //         out <<" Frame = " << ((master_frame > 0) ? "+" : "")
    //             << master_frame << "\n";
    //     }  else if (slave_frame != 0){
    //         out <<" Frame = " << ((slave_frame > 0) ? "+" : "")
    //             << slave_frame <<"\n";
    //     }
    //     out<<"\n";
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:321-331
    // ```c++
    //     if(master_frame != 0 && slave_frame != 0) {
    //         ...
    //     } else if (master_frame != 0){
    //         out <<" Frame = " << ((master_frame > 0) ? "+" : "")
    //             << master_frame << "\n";
    // ```
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

// NCBI c++/src/objtools/align_format/showalign.cpp:1600-1625,2122-2149:
// AddSpace(out,maxStartLen-startLen+k_StartSequenceMargin); out << sequence;
// if (sequence_standard[i]==sequence[i]) middle_line[i]=sequence[i];
// else if (m_Matrix[...]>0) middle_line[i]='+';
// The query's row of letters and the subject's translated row (width 3, three nucleotides per
// residue, in the direction of its frame; `write_hsp_rows`).
fn write_tblastn_alignment<W: Write>(
    writer: &mut W,
    hit: &PairwiseHit,
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
    write_hsp_rows(
        writer,
        |segments| tblastn_row_segments(hit, segments),
        [1, 3],
        config.line_length,
        q.len(),
        width,
        |writer, row, columns| writer.write_all(&[q, s][row][columns]),
        |writer, columns| write_protein_middle(writer, q, s, columns, config.protein_matrix, false),
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:387-407
// ```c++
//       CBlastFormatUtil::BlastPrintVersionInfo(m_Program, m_IsHTML,
//                                               m_Outfile);
//     }
//
//     if (m_IsBl2Seq && !m_IsDbScan) {
//         return;
//     }
//
//     m_Outfile << NcbiEndl << NcbiEndl;
//     if (m_Program == "deltablast") {
//         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
//                               m_Outfile, CReference::eDeltaBlast);
//         m_Outfile << "\n";
//     }
//
//     if (m_Megablast)
//         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
//                                           m_Outfile, CReference::eMegaBlast);
//     else
//         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
//                                           m_Outfile);
// ```
fn write_translated_pairwise_intro(
    writer: &mut impl Write,
    program: &str,
    version: &str,
) -> io::Result<()> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:66-73
    // ```c++
    // void CBlastFormatUtil::BlastPrintVersionInfo(const string program, bool html,
    //                                              CNcbiOstream& out)
    // {
    //     if (html)
    //         out << "<b>" << BlastGetVersion(program) << "</b>" << "\n";
    //     else
    //         out << BlastGetVersion(program) << "\n";
    // }
    // ```
    super::reference::write_version_line(writer, false, program, version)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:395-395
    // ```c++
    //     m_Outfile << NcbiEndl << NcbiEndl;
    // ```
    writeln!(writer)?;
    writer.flush()?;
    writeln!(writer)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:395-395
    // ```c++
    //     m_Outfile << NcbiEndl << NcbiEndl;
    // ```
    writer.flush()?;
    super::reference::write_prolog_references(
        writer,
        false,
        super::reference::DEFAULT_PROTOCOL,
        program,
        false,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:275-282
// ```c++
//     if ((m_FormatType == CFormattingArgs::eXml2) || (m_FormatType == CFormattingArgs::eJson) ||
//         (m_FormatType == CFormattingArgs::eXml2_S) || (m_FormatType == CFormattingArgs::eJson_S)) {
//            m_AccumulatedQueries.Reset(new CBlastQueryVector());
//     }
//
//     if (opts.GetSumStatisticsMode() && m_IsUngappedSearch) {
//         m_ShowLinkedSetSize = true;
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1548-1565
// ```c++
//     results.GetMaskedQueryRegions(masklocs);
//
//     CSeq_align_set copy_aln_set;
//     CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_NumAlignments);
//
//     int flags = s_SetFlags(m_Program, m_FormatType, m_IsHTML, m_ShowGi,
//                            (m_IsBl2Seq && !m_IsDbScan), (m_DisableKAStats || kIsGlobal));
//
//     CDisplaySeqalign display(copy_aln_set, *m_Scope, &masklocs, NULL, m_MatrixName);
//     display.SetDbName(m_DbName);
//     display.SetDbType(!m_DbIsAA);
//     display.SetLineLen(m_LineLength);
//     int kAlignToShow=2000000000;  // Nice large number per SB-1817
//     display.SetNumAlignToShow(kAlignToShow);
//
//     // set the alignment flags
//     display.SetAlignOption(flags);
//
// ```
pub struct BlastxPairwiseOptions {
    pub word_threshold: f64,
    pub gapped: bool,
    pub sum_stats: bool,
    pub num_descriptions: usize,
    pub num_alignments: usize,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1490-1519
// ```c++
//
//     m_Outfile << "\n\n";
//     CBlastFormatUtil::AcknowledgeBlastQuery(*bioseq, kFormatLineLength,
//                                             m_Outfile, m_BelieveQuery,
//                                             m_IsHTML, kIsTabularOutput,
//                                             results.GetRID());
//
//     if (m_IsBl2Seq && !m_IsDbScan) {
//         m_Outfile << "\n";
//         // FIXME: this might be configurable in the future
//         const bool kBelieveSubject = false;
//         CConstRef<CBioseq> subject_bioseq = x_CreateSubjectBioseq();
//         CBlastFormatUtil::AcknowledgeBlastSubject(*subject_bioseq,
//                                                   kFormatLineLength,
//                                                   m_Outfile, kBelieveSubject,
//                                                   m_IsHTML, kIsTabularOutput);
//     }
//
//     // quit early if there are no hits
//     if ( !results.HasAlignments() ) {
//         m_Outfile << "\n\n"
//               << "***** " << CBlastFormatUtil::kNoHitsFound << " *****" << "\n"
//               << "\n\n";
//         x_PrintOneQueryFooter(*results.GetAncillaryData());
//         return;
//     }
//
//     CConstRef<CSeq_align_set> aln_set = results.GetSeqAlign();
//     _ASSERT(results.HasAlignments());
//     if (m_IsUngappedSearch) {
// ```
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
pub fn write_blastx_pairwise_report<W: Write>(
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
    probe: Option<&mut FormatProbe<'_>>,
    first_hsp: HspIndex,
    // The sequences of the query-anchored sections (outfmt 1 to 4, `config.view`; `None` for
    // outfmt 0).
    anchored: Option<&dyn super::query_anchored::ReportSequences>,
    pool: Option<&crate::utils::threading::SearchPool<'_>>,
    loads: &crate::report::db_subject::HitLoads,
) -> io::Result<()> {
    let sections = BlastxSections {
        config,
        options,
        queries,
        subject_labels,
        report,
        database,
        records_parsed,
        record_titles,
        anchored,
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
    /// A `-db` report's subjects (`None` for `-subject`).
    database: Option<BlastxDatabaseSubjects<'a>>,
    /// A `-subject` report's records with parsed Seq-ids, by subject index.
    records_parsed: &'a [Option<Arc<crate::report::db_subject::ParsedSubject>>],
    /// A `-subject` report's record titles (their headings' `GenerateDefline`), by subject index.
    record_titles: &'a [&'a [u8]],
    /// The sequences of the query-anchored sections (outfmt 1 to 4; `None` for outfmt 0).
    anchored: Option<&'a dyn super::query_anchored::ReportSequences>,
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
    ) -> io::Result<()> {
        let config = self.config;
        // NCBI c++/src/objtools/align_format/showalign.cpp:3650-3668:
        // display leaves an extra blank line after each translated HSP.
        let write_hsp = |writer: &mut W, hit: &PairwiseHit| {
            write_tblastn_hsp_info(writer, hit)?;
            write_blastx_alignment(writer, hit, config)?;
            writeln!(writer)
        };
        let (hits, indices) = (
            &group.subject_hits[&s_idx],
            &group.subject_hit_indices[&s_idx],
        );
        let parsed = match self.database {
            Some((_, _, parsed, _)) => parsed_subject(parsed, s_idx),
            None => parsed_subject(self.records_parsed, s_idx),
        };
        let heading = match (parsed, self.database) {
            (Some(subject), _) => Some(DatabaseHeading::loaded(
                subject,
                self.options.gapped,
                group.own_hsps(s_idx),
                super::x_pairwise_par::QuerySections::thrown(self, s_idx).is_some(),
            )),
            (None, Some((_, titles, _, sources))) => Some(DatabaseHeading::Title(
                &titles[s_idx as usize],
                stored_source(sources, s_idx),
            )),
            (None, None) => None,
        };
        if let Some(heading) = heading {
            return write_database_subject_alignments_part(
                writer, probe, heading, true, hits, indices, part, write_hsp,
            );
        }
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2273
        // ```c++
        // 	alnDispParams->title = CDeflineGenerator().GenerateDefline(bsp_handle);
        // ```
        let title = self
            .record_titles
            .get(s_idx as usize)
            .copied()
            .unwrap_or_default();
        write_subject_alignments_part(
            writer,
            probe,
            s_idx,
            self.records_parsed,
            title,
            true,
            hits,
            indices,
            part,
            write_hsp,
        )
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
    // ```c++
    //     if(m_HitsSortOption > 0 || m_HspsSortOption > 0) {
    //         aln_set = CBlastFormatUtil::SortSeqalignForSortableFormat(
    // ```
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

    fn head<W: Write>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: Option<&super::x_pairwise_par::QueryGroup<'_>>,
        read: super::query_warnings::SubjectReadHook<'_, W>,
    ) -> io::Result<()> {
        let query = &self.queries[q_idx];
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1491-1492
        // ```c++
        //     m_Outfile << "\n\n";
        //     CBlastFormatUtil::AcknowledgeBlastQuery(*bioseq, kFormatLineLength,
        // ```
        writeln!(writer)?;
        writeln!(writer)?;
        write_blastp_query_header(writer, &query.query_name, query.query_length)?;
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
        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:606-611
        // ```c++
        //     CShowBlastDefline showdef(*aln_set, *m_Scope,
        //                               defline_length == -1 ? kFormatLineLength:defline_length,
        //                               m_NumSummary + additional);
        // ```
        // The table follows `x_InitDeflineTable` (`write_description_table`): the subject's
        // highest bit score and that HSP's E-value, for the first `m_NumDescriptions` subjects;
        // the last one is NCBI's "last hit" when every subject fits.
        let descriptions = group.described(self.options.num_descriptions);
        // The table reads the ranges of every HSP of each described subject first: an HSP with
        // empty rows stops the report after the query's `Length=` (`describes_empty_rows`).
        if describes_empty_rows(descriptions, grouped) {
            writer.flush()?;
            return Err(io::Error::other(EmptyRowsStop::RowIsEmpty));
        }
        write_subject_summary_table_with_sum_n(
            writer,
            descriptions,
            order.len() <= self.options.num_descriptions,
            grouped,
            self.subject_labels,
            !self.options.gapped && self.options.sum_stats,
            self.config.sort_columns(query.query_length),
            read,
        )?;
        writeln!(writer)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1948-1973
    // ```c++
    //             CRef<CAlnVec> avRef = x_GetAlnVecForSeqalign(**iter);
    //             ...
    //                         bool showDefLine = previousId.Empty() || !subid->Match(*previousId);
    //                         x_DisplayAlnvecInfo(out, alnvecInfo,showDefLine);
    // ```
    // Each HSP's alignment is made before its heading and block: the Std-seg of an HSP with empty
    // rows that `-subject`'s Seq-loc mapper rewrote throws there (`CreateDensegFromStdseg`), before
    // anything of the HSP (`EmptyRowsStop::StdSegLengths`), the subject's heading for its first.
    fn subject<W: Write>(
        &self,
        writer: &mut W,
        probe: Option<&mut FormatProbe<'_>>,
        group: &super::x_pairwise_par::QueryGroup<'_>,
        s_idx: u32,
        part: std::ops::Range<usize>,
    ) -> io::Result<()> {
        let hits = &group.subject_hits[&s_idx];
        let stop = mapper_rewrites(self.remaps(), hits)
            .then(|| part.clone().find(|&hsp| hits[hsp].rows_empty()))
            .flatten();
        let written = part.start..stop.unwrap_or(part.end);
        if !written.is_empty() {
            self.write_part(writer, probe, group, s_idx, written)?;
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

    // Each HSP's Dense-seg is made after its subject's read, before the section's first byte and
    // outside its `try`: the Std-seg of an HSP with empty rows that `-subject`'s Seq-loc mapper
    // rewrote throws there (`EmptyRowsStop::StdSegLengths`, as outfmt 0's alignments).
    fn merged_subject<W: Write>(
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1550-1552
    // ```c++
    //     CSeq_align_set copy_aln_set;
    //     CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_NumAlignments);
    // ```
    /// The query-anchored section of the query's aligned subjects (`blastx_section_hsps`).
    fn alignments<W: Write>(
        &self,
        writer: &mut W,
        q_idx: usize,
        group: &super::x_pairwise_par::QueryGroup<'_>,
    ) -> io::Result<()> {
        write_anchored_section(writer, self.config, self.anchored, q_idx, |_| {
            blastx_section_hsps(self.options.num_alignments, group)
        })
    }

    fn foot<W: Write>(&self, writer: &mut W, q_idx: usize) -> io::Result<()> {
        write_blastx_query_footer(
            writer,
            &self.queries[q_idx],
            self.report,
            self.options.gapped,
        )
    }
}
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:254-254
// ```c++
//         formatter.PrintProlog();
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:297-297
// ```c++
//         formatter.PrintEpilog(opt);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:387-388
// ```c++
//       CBlastFormatUtil::BlastPrintVersionInfo(m_Program, m_IsHTML,
//                                               m_Outfile);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:406-407
// ```c++
//         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
//                                           m_Outfile);
// ```
/// `program` is the program line's name: `BLASTX`, or `RPSTBLASTN`, whose prolog is BLASTX's (the
/// Gapped BLAST reference only, blast_format.cpp:416-426 quoted at `write_translated_pairwise_intro`).
pub fn write_blastx_prolog(
    writer: &mut impl Write,
    program: &str,
    report: &BlastpPairwiseReport,
    database: Option<&PairwiseDatabase>,
) -> io::Result<()> {
    write_translated_pairwise_intro(writer, program, &report.version)?;
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
// ```c++
//                         options.GetMismatchPenalty() << "\n";
//     }
//     else {
//         m_Outfile << "\n\nMatrix: " << options.GetMatrixName() << "\n";
//     }
//
//     if (options.GetGappedMode() == true) {
//         double gap_extension = (double) options.GetGapExtensionCost();
//         if ((m_Program == "megablast" || m_Program == "blastn") && options.GetGapExtensionCost() == 0)
//         { // Formula from PMID 10890397 applies if both gap values are zero.
//                gap_extension = -2*options.GetMismatchPenalty() + options.GetMatchReward();
//                gap_extension /= 2.0;
//         }
//         m_Outfile << "Gap Penalties: Existence: "
//                 << options.GetGapOpeningCost() << ", Extension: "
//                 << gap_extension << "\n";
//     }
//     if (options.GetWordThreshold()) {
//         m_Outfile << "Neighboring words threshold: " <<
//                         options.GetWordThreshold() << "\n";
//     }
//     if (options.GetWindowSize()) {
//         m_Outfile << "Window for multiple hits: " <<
//                         options.GetWindowSize() << "\n";
//     }
// ```
pub fn write_blastx_epilog(
    writer: &mut impl Write,
    report: &BlastpPairwiseReport,
    options: &BlastxPairwiseOptions,
    database: Option<&PairwiseDatabase>,
) -> io::Result<()> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2264-2288
    // ```c++
    //                         options.GetMismatchPenalty() << "\n";
    //     }
    //     else {
    //         m_Outfile << "\n\nMatrix: " << options.GetMatrixName() << "\n";
    //     }
    //
    //     if (options.GetGappedMode() == true) {
    //         double gap_extension = (double) options.GetGapExtensionCost();
    //         if ((m_Program == "megablast" || m_Program == "blastn") && options.GetGapExtensionCost() == 0)
    //         { // Formula from PMID 10890397 applies if both gap values are zero.
    //                gap_extension = -2*options.GetMismatchPenalty() + options.GetMatchReward();
    //                gap_extension /= 2.0;
    //         }
    //         m_Outfile << "Gap Penalties: Existence: "
    //                 << options.GetGapOpeningCost() << ", Extension: "
    //                 << gap_extension << "\n";
    //     }
    //     if (options.GetWordThreshold()) {
    //         m_Outfile << "Neighboring words threshold: " <<
    //                         options.GetWordThreshold() << "\n";
    //     }
    //     if (options.GetWindowSize()) {
    //         m_Outfile << "Window for multiple hits: " <<
    //                         options.GetWindowSize() << "\n";
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2249-2249
    // ```c++
    //     m_Outfile << NcbiEndl << NcbiEndl;
    // ```
    writeln!(writer)?;
    writeln!(writer)?;
    if let Some(database) = database {
        write_final_database_block(writer, database)?;
    } else {
        write!(writer, "  Database: ")?;
        write_flatfile_wrapped(
            writer,
            ensure_trailing_period(&report.database_name).as_bytes(),
            68,
        )?;
        writeln!(writer, "    Posted date:  Unknown")?;
        writeln!(
            writer,
            "  Number of letters in database: {}",
            crate::ncbi_str::int8_with_commas(report.database_total_letters as i64)
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
    if options.word_threshold != 0.0 {
        // `GetWordThreshold()` is a double, written with the stream's default format.
        writeln!(
            writer,
            "Neighboring words threshold: {}",
            cpp_default_double(options.word_threshold)
        )?;
    }
    if report.window_size != 0 {
        writeln!(writer, "Window for multiple hits: {}", report.window_size)?;
    }
    Ok(())
}
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:445-478
// ```c++
// CBlastFormat::x_PrintOneQueryFooter(const blast::CBlastAncillaryData& summary)
// {
//     /* Skip printing KA parameters if the program is rmblastn -RMH- */
//     if ( m_DisableKAStats )
//       return;
//
//     const Blast_KarlinBlk *kbp_ungap =
//         (m_Program == "psiblast" || m_Program == "deltablast")
//         ? summary.GetPsiUngappedKarlinBlk()
//         : summary.GetUngappedKarlinBlk();
//     const Blast_GumbelBlk *gbp = summary.GetGumbelBlk();
//     m_Outfile << NcbiEndl;
//     if (kbp_ungap) {
//         CBlastFormatUtil::PrintKAParameters(kbp_ungap->Lambda,
//                                             kbp_ungap->K, kbp_ungap->H,
//                                             kFormatLineLength, m_Outfile,
//                                             false, gbp);
//     }
//
//     const Blast_KarlinBlk *kbp_gap =
//         (m_Program == "psiblast" || m_Program == "deltablast")
//         ? summary.GetPsiGappedKarlinBlk()
//         : summary.GetGappedKarlinBlk();
//     m_Outfile << "\n";
//     if (kbp_gap) {
//         CBlastFormatUtil::PrintKAParameters(kbp_gap->Lambda,
//                                             kbp_gap->K, kbp_gap->H,
//                                             kFormatLineLength, m_Outfile,
//                                             true, gbp);
//     }
//
//     m_Outfile << "\n";
//     m_Outfile << "Effective search space used: " <<
//                         summary.GetSearchSpace() << "\n";
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:583-613
// ```c++
//
//     char buffer[256];
//     if (gapped) {
//         out << "Gapped" << "\n";
//     }
//     out << "Lambda      K        H";
//     if (gbp) {
//         if (gapped) {
//             out << "        a         alpha    sigma";
//         } else {
//             out << "        a         alpha";
//         }
//     }
//     out << "\n";
//     sprintf(buffer, "%#8.3g ", lambda);
//     out << buffer;
//     sprintf(buffer, "%#8.3g ", k);
//     out << buffer;
//     sprintf(buffer, "%#8.3g ", h);
//     out << buffer;
//     if (gbp) {
//         if (gapped) {
//             sprintf(buffer, "%#8.3g ", gbp->a);
//             out << buffer;
//             sprintf(buffer, "%#8.3g ", gbp->Alpha);
//             out << buffer;
//             sprintf(buffer, "%#8.3g ", gbp->Sigma);
//             out << buffer;
//         } else {
//             sprintf(buffer, "%#8.3g ", gbp->a_un);
//             out << buffer;
// ```
/// The footer of a BLASTX query (`BlastpPairwiseQuery::valid` and `batch_skipped`: its Karlin
/// blocks, `gapped`: the gapped one).
fn write_blastx_query_footer(
    writer: &mut impl Write,
    query: &BlastpPairwiseQuery,
    report: &BlastpPairwiseReport,
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
            report.gapped_karlin,
            report.gumbel,
            query.effective_search_space,
            false,
        );
    }
    writeln!(writer)?;
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
// ```c++
// {
//     size_t startLen=0;
//     int start = alnRoInfo->seqStarts[row].front() + 1;  //+1 for 1 based
//     int end = alnRoInfo->seqStops[row].front() + 1;
//     int j = alnRoInfo->currPrintSegment;
//     int actualLineLen = alnRoInfo->currActualLineLen;
//     //print out sequence line
//     //adjust space between id and start
//     CAlignFormatUtil::AddSpace(out, alnRoInfo->maxIdLen-alnRoInfo->seqidArray[row].size() + k_IdStartMargin);
//     //not to display start and stop number for empty row
//     if ((j > 0 && end == prev_stop)
//         || (j == 0 && start == 1 && end == 1)) {
//         startLen = 0;
//     } else {
//         out << start;
//         startLen=NStr::IntToString(start).size();
//     }
//
//     CAlignFormatUtil::AddSpace(out, alnRoInfo->maxStartLen-startLen + k_StartSequenceMargin);
//     x_OutputSeq(alnRoInfo->sequence[row], m_AV->GetSeqId(row), j,
//                 (int)actualLineLen, alnRoInfo->frame[row], row,
//                 (row > 0 && alnRoInfo->colorMismatch)?true:false,
//                 alnRoInfo->masked_regions[row], out);
//     CAlignFormatUtil::AddSpace(out, k_SeqStopMargin);
//
//      //not to display stop number for empty row in the middle
//     if (!(j > 0 && end == prev_stop)
//         && !(j == 0 && start == 1 && end == 1)) {
//         out << end;
//     }
//     out<<"\n";
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:120-134
// ```c++
//
//     if (strand == eNa_strand_minus) {
//
//         if (translate)
//             retval = original_length -
//                 CODON_LENGTH*(s_GetCurrPos(curr_pos, num) + num)
//                 + frame + 1;
//         else
//             retval = length - s_GetCurrPos(curr_pos, num) - num;
//
//     } else {
//
//         if (translate)
//             retval = frame - 1 + CODON_LENGTH*s_GetCurrPos(curr_pos, num);
//         else
// ```
fn write_blastx_alignment(
    writer: &mut impl Write,
    hit: &PairwiseHit,
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
    // ```c++
    //         x_FillSeqid(seqidArray[row], row);
    //         maxIdLen=max<size_t>(seqidArray[row].size(), maxIdLen);
    //         size_t maxCood=max<size_t>(m_AV->GetSeqStart(row), m_AV->GetSeqStop(row));
    //         maxStartLen = max<size_t>(NStr::SizetToString(maxCood).size(), maxStartLen);
    // ```
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
    // Dense-seg (`write_hsp_rows`).
    if q.len() != s.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "BLASTX alignment string lengths differ",
        ));
    }
    write_hsp_rows(
        writer,
        |segments| blastx_row_segments(hit, segments),
        [3, 1],
        config.line_length,
        q.len(),
        width,
        |writer, row, columns| writer.write_all(&[q, s][row][columns]),
        |writer, columns| write_protein_middle(writer, q, s, columns, config.protein_matrix, true),
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3152-3161
// ```c++
// int CAlignFormatUtil::GetPercentMatch(int numerator, int denominator)
// {
//      if (numerator == denominator)
//         return 100;
//      else {
//        int retval =(int) (0.5 + 100.0*((double)numerator)/((double)denominator));
//        retval = min(99, retval);
//        return retval;
//      }
// }
// ```
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
mod tblastx_tests {
    use super::*;

    #[test]
    fn thresholds_are_written_as_a_cpp_stream_writes_a_double() {
        for (value, text) in [
            (11.0, "11"),
            (999999.0, "999999"),
            (1000000.0, "1e+06"),
            (1234567.0, "1.23457e+06"),
            (2147483647.0, "2.14748e+09"),
            (12.5, "12.5"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
        ] {
            assert_eq!(cpp_default_double(value), text, "{value}");
        }
    }

    /// The rows of a protein alignment `query` x `subject` (`-` gaps) from (`q_start`, `s_start`),
    /// with the edit script of their columns (`gap_info`: a gap in the subject is an insertion,
    /// one in the query a deletion), from which the report makes the HSP's Dense-seg.
    fn rows_hit(query: &str, subject: &str, q_start: usize, s_start: usize) -> Hit {
        use crate::common::GapEditOp;
        let residues = |row: &str| row.bytes().filter(|&b| b != b'-').count();
        let mut ops: Vec<GapEditOp> = Vec::new();
        for (q, s) in query.bytes().zip(subject.bytes()) {
            let op = match (q == b'-', s == b'-') {
                (false, true) => GapEditOp::Ins(1),
                (true, false) => GapEditOp::Del(1),
                _ => GapEditOp::Sub(1),
            };
            match (ops.last_mut(), op) {
                (Some(GapEditOp::Sub(num)), GapEditOp::Sub(_))
                | (Some(GapEditOp::Ins(num)), GapEditOp::Ins(_))
                | (Some(GapEditOp::Del(num)), GapEditOp::Del(_)) => *num += 1,
                _ => ops.push(op),
            }
        }
        Hit {
            identity: 0.0,
            length: query.len(),
            mismatch: 0,
            gapopen: 1,
            q_start,
            q_end: q_start + residues(query) - 1,
            s_start,
            s_end: s_start + residues(subject) - 1,
            e_value: 1e-10,
            bit_score: 50.0,
            num_ident: 0,
            query_frame: 0,
            query_length: 0,
            q_idx: 0,
            s_idx: 0,
            raw_score: 100,
            sort_query_offset: 0,
            sort_query_end: 0,
            sort_subject_offset: 0,
            sort_subject_end: 0,
            has_sort_offsets: false,
            gap_info: Some(ops),
            num_positives: 0,
        }
    }

    fn rows_text(
        query: &str,
        subject: &str,
        q_start: usize,
        s_start: usize,
        line: usize,
    ) -> String {
        let hit = rows_hit(query, subject, q_start, s_start);
        let config = PairwiseConfig {
            line_length: line,
            program: "rpsblast".to_string(),
            ..PairwiseConfig::default()
        };
        let mut out = Vec::new();
        write_alignment_with_sequences(&mut out, &hit, query, subject, &[], &config).unwrap();
        String::from_utf8(out).unwrap()
    }

    // showalign.cpp:1209-1216 and 1497-1502 (quoted at `end_gap_span` and
    // `write_alignment_with_sequences`): a row's end gaps are spaces, a chunk without the row's
    // letters leaves the row out (the Query row with its middle line), and the next chunk's row
    // shows its first residue's number.
    #[test]
    fn end_gaps_are_spaces_and_rows_without_letters_are_left_out() {
        let middle = |letters: &str| format!("{}{letters}\n", " ".repeat(10));
        // The subject starts with six end gaps: the first chunk has no Sbjct line, the second
        // starts it at its residue 1 after two spaces.
        let text = rows_text("MKVLWAALLV", "------AALL", 1, 1, 4);
        let expected = [
            "Query  1  MKVL  4\n".to_string(),
            middle("    "),
            "\n".to_string(),
            "Query  5  WAAL  8\n".to_string(),
            middle("  A "),
            "Sbjct  1    AA  2\n\n".to_string(),
            "Query  9  LV  10\n".to_string(),
            middle("L+"),
            "Sbjct  3  LL  4\n\n".to_string(),
        ]
        .concat();
        assert_eq!(text, expected);
        // The query starts with three end gaps.
        let text = rows_text("---KVLW", "MQLKVLW", 7, 1, 60);
        let expected = [
            "Query  7     KVLW  10\n".to_string(),
            middle("   KVLW"),
            "Sbjct  1  MQLKVLW  7\n\n".to_string(),
        ]
        .concat();
        assert_eq!(text, expected);
        // A query chunk of end gaps only: no Query line and no middle line.
        let text = rows_text("----KVLW", "MQLWKVLW", 1, 1, 4);
        let expected = [
            "Sbjct  1  MQLW  4\n\n".to_string(),
            "Query  1  KVLW  4\n".to_string(),
            middle("KVLW"),
            "Sbjct  5  KVLW  8\n\n".to_string(),
        ]
        .concat();
        assert_eq!(text, expected);
    }

    fn tblastx_hit(sum_n: i32) -> PairwiseHit {
        let query: String = "ACDEFGHIKL".repeat(6) + "M";
        let subject: String = "ACDEFGHIKL".repeat(6) + "N";
        PairwiseHit {
            hit: Hit {
                identity: 0.0,
                length: 61,
                mismatch: 1,
                gapopen: 0,
                q_start: 200,
                q_end: 18,
                s_start: 3,
                s_end: 185,
                e_value: 2e-30,
                bit_score: 120.5,
                num_ident: 60,
                query_frame: -2,
                query_length: 0,
                q_idx: 0,
                s_idx: 0,
                raw_score: 255,
                sort_query_offset: 0,
                sort_query_end: 61,
                sort_subject_offset: 0,
                sort_subject_end: 61,
                has_sort_offsets: true,
                gap_info: None,
                num_positives: 60,
            },
            query_seq: Some(query),
            subject_seq: Some(subject),
            query_frame: Some(-2),
            subject_frame: Some(3),
            positives: Some(60),
            gaps: Some(0),
            subject_length: Some(400),
            subject_title: None,
            comp_adjust_method: None,
            sum_n: Some(sum_n),
            subject_stored_gaps: Vec::new(),
            seqalign_num_ident: 60,
        }
    }

    // Line DBID item ID4t: TBLASTX's outfmt 0 makes its rows a window of queries at a time and
    // writes each window's report (`algorithm::tblastx::blast_engine::run_impl`), as BLASTN's
    // (ID4f, above): the reports of consecutive query ranges, each given only its own hits, are
    // the report of all queries, on one thread and on a pool.
    #[test]
    #[cfg_attr(
        not(feature = "parallel"),
        ignore = "needs threads: a build without the parallel feature has none"
    )]
    fn tblastx_reports_of_consecutive_query_ranges_are_the_whole_report() {
        let config = PairwiseConfig {
            program: "tblastx".to_string(),
            ..PairwiseConfig::default()
        };
        let karlin = KarlinParams {
            lambda: 0.3176,
            k: 0.134,
            h: 0.4012,
            alpha: 0.0,
            beta: 0.0,
        };
        let queries: Vec<TblastxPairwiseQuery> = (0..4)
            .map(|q_idx| TblastxPairwiseQuery {
                query_name: format!("q{q_idx} query {q_idx}").into_bytes(),
                query_length: 300,
                karlin: Some(karlin),
                effective_search_space: 10_000,
            })
            .collect();
        let report = TblastxPairwiseReport {
            version: "2.17.0+".to_string(),
            database_name: "subjects".to_string(),
            database_num_sequences: 3,
            database_total_letters: 900,
            matrix_name: "BLOSUM62".to_string(),
            word_threshold: 13.0,
            window_size: 40,
            num_descriptions: 500,
            num_alignments: 250,
            // Query 1 is not searched (a batch without a valid context).
            unsearched: vec![false, true, false, false],
            epilog: true,
            prolog: true,
            parsed_subjects: Vec::new(),
            subject_sources: Vec::new(),
        };
        // Query 0: subjects 1 (2 HSPs) and 2; query 1: no hits; queries 2 and 3: subject 1.
        let hits: Vec<PairwiseHit> = [(0, 1, 4), (0, 1, 4), (0, 2, 1), (2, 1, 1), (3, 1, 1)]
            .into_iter()
            .map(|(q_idx, s_idx, sum_n)| {
                let mut hit = tblastx_hit(sum_n);
                hit.hit.q_idx = q_idx;
                hit.hit.s_idx = s_idx;
                hit
            })
            .collect();
        let titles: Vec<Arc<[u8]>> = (0..3)
            .map(|s_idx| Arc::from(format!("s{s_idx} subject {s_idx}").as_bytes()))
            .collect();
        // The warnings written before queries 0 and 2 (`QueryWarnings`, by query index).
        let before: Vec<crate::report::diag_log::Posts> = (0..queries.len())
            .map(|q_idx| match q_idx {
                0 | 2 => crate::report::diag_log::Posts::of_lines(
                    format!(
                        "Warning: [tblastx] Query_{} q{q_idx}: a warning\n",
                        q_idx + 1
                    )
                    .into_bytes(),
                ),
                _ => crate::report::diag_log::Posts::new(),
            })
            .collect();
        crate::utils::threading::with_search_pool(4, "test", |pool| {
            for x_pool in [None, Some(pool)] {
                // The report and the warnings of the queries `range`.
                let write = |range: std::ops::Range<usize>| {
                    let range_hits: Vec<PairwiseHit> = hits
                        .iter()
                        .filter(|hit| range.contains(&(hit.hit.q_idx as usize)))
                        .cloned()
                        .collect();
                    let mut written = Vec::new();
                    let mut diagnostics = Vec::new();
                    write_tblastx_pairwise_report(
                        &range_hits,
                        &mut written,
                        &config,
                        &queries,
                        range,
                        &titles,
                        &report,
                        None,
                        None,
                        Some(&mut crate::report::query_warnings::QueryWarnings::of_posts(
                            &before,
                            None,
                            &mut diagnostics,
                        )),
                        x_pool,
                        None,
                        &crate::report::db_subject::HitLoads::default(),
                    )
                    .unwrap();
                    (written, diagnostics)
                };
                let (whole, whole_warnings) = write(0..queries.len());
                let text = String::from_utf8(whole.clone()).unwrap();
                assert!(text.starts_with("TBLASTX 2.17.0+"), "{text}");
                assert_eq!(text.matches("Query= ").count(), queries.len(), "{text}");
                assert!(text.contains("No hits found"), "{text}");
                assert!(text.contains("Expect(4)"), "{text}");
                assert!(text.ends_with("Window for multiple hits: 40\n"), "{text}");
                assert_eq!(text.matches("Window for multiple hits").count(), 1, "{text}");
                assert_eq!(
                    String::from_utf8(whole_warnings.clone()).unwrap(),
                    "Warning: [tblastx] Query_1 q0: a warning\nWarning: [tblastx] Query_3 q2: a warning\n"
                );
                for cuts in [vec![1], vec![2], vec![3], vec![1, 2, 3], vec![1, 3]] {
                    let mut joined = Vec::new();
                    let mut joined_warnings = Vec::new();
                    let mut start = 0;
                    for end in cuts.into_iter().chain([queries.len()]) {
                        let (written, warnings) = write(start..end);
                        joined.extend(written);
                        joined_warnings.extend(warnings);
                        start = end;
                    }
                    assert_eq!(String::from_utf8(joined).unwrap(), text);
                    assert_eq!(joined_warnings, whole_warnings);
                }
            }
            Ok(())
        })
        .unwrap();
    }

    // NCBI c++/src/objtools/align_format/showalign.cpp:3593-3606 and 310-332: Expect(n)
    // only for a linked set of two or more, both frames, no Strand line.
    #[test]
    fn tblastx_hsp_info_prints_the_linked_set_and_both_frames() {
        let mut out = Vec::new();
        write_tblastx_hsp_info(&mut out, &tblastx_hit(4)).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            " Score = 120 bits (255),  Expect(4) = 2e-30\n Identities = 60/61 (98%), Positives = 60/61 (98%), Gaps = 0/61 (0%)\n Frame = -2/+3\n\n"
        );
        let mut out = Vec::new();
        write_tblastx_hsp_info(&mut out, &tblastx_hit(1)).unwrap();
        assert!(String::from_utf8(out)
            .unwrap()
            .starts_with(" Score = 120 bits (255),  Expect = 2e-30\n"));
    }

    // NCBI c++/src/objtools/alnmgr/alnvec.cpp:230,305,358 and showalign.cpp:1598-1626: 60
    // residues span 180 nucleotides, each row in the direction of its frame, and the width
    // is that of the largest 0-based coordinate (199: 3 digits).
    #[test]
    fn tblastx_rows_step_three_nucleotides_in_each_frame_direction() {
        let mut out = Vec::new();
        // The query's record has 300 nucleotides (frame -2 at 200).
        write_tblastx_alignment(&mut out, &tblastx_hit(1), 300, &PairwiseConfig::default())
            .unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        let row = "ACDEFGHIKL".repeat(6);
        assert_eq!(lines[0], format!("Query  200  {row}  21"));
        assert_eq!(lines[1], format!("            {row}"));
        assert_eq!(lines[2], format!("Sbjct  3    {row}  182"));
        assert_eq!(lines[3], "");
        assert_eq!(lines[4], "Query  20   M  18");
        assert_eq!(lines[5], "             ");
        assert_eq!(lines[6], "Sbjct  183  N  185");
    }

    // NCBI c++/src/objtools/align_format/showdefline.cpp:829-836,961-965: the N column of
    // an ungapped sum-statistics search, from the first HSP of each subject.
    #[test]
    fn description_table_shows_the_linked_set_column() {
        let first = tblastx_hit(12);
        let mut second = tblastx_hit(1);
        second.hit.s_idx = 1;
        second.hit.bit_score = 30.2;
        second.hit.e_value = 0.5;
        let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
            [(0, vec![&first]), (1, vec![&second])]
                .into_iter()
                .collect();
        let titles: std::collections::HashMap<u32, Vec<u8>> =
            [(0, b"s0".to_vec()), (1, b"s1".to_vec())]
                .into_iter()
                .collect();
        let mut out = Vec::new();
        write_blastn_description_table(
            &mut out,
            &[0, 1],
            true,
            &hits,
            &titles,
            &[],
            &[],
            true,
            false,
            None,
            &mut crate::report::query_warnings::no_subject_reads,
            &crate::report::db_subject::HitLoads::default(),
        )
        .unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert!(lines[0].ends_with("Score     E"), "{text}");
        assert!(lines[1].ends_with("(Bits)  Value  N"), "{text}");
        assert!(lines[3].ends_with("  120     2e-30  12"), "{text}");
        assert!(lines[4].ends_with("  30.2    0.50   1 "), "{text}");
    }

    // NCBI c++/src/objtools/align_format/showdefline.cpp:771-840,930-975 (`-sorthits`): the
    // Total score, Query cover and Max ident columns, as NCBI BLAST+ 2.17.0 `blastp -sorthits 0`
    // writes them (probe t1 `blastp.s.h0kx`): the total in the HSPs' order, the cover of the
    // union of the query ranges, the identity of the best HSP (its Seq-align's `num_ident`,
    // `PairwiseHit::seqalign_num_ident`), and the trailing spaces of the nine letters of
    // "Max ident".
    #[test]
    fn sorthits_adds_the_total_cover_and_identity_columns() {
        let mut a = tblastx_hit(1);
        (
            a.hit.bit_score,
            a.hit.e_value,
            a.seqalign_num_ident,
            a.hit.length,
        ) = (28.1, 0.18, 27, 99);
        (a.hit.q_start, a.hit.q_end) = (1, 99);
        let mut b1 = a.clone();
        b1.hit.s_idx = 1;
        (
            b1.hit.bit_score,
            b1.hit.e_value,
            b1.seqalign_num_ident,
            b1.hit.length,
        ) = (24.6, 2.3, 16, 25);
        (b1.hit.q_start, b1.hit.q_end) = (10, 34);
        let mut b2 = b1.clone();
        (
            b2.hit.bit_score,
            b2.hit.e_value,
            b2.seqalign_num_ident,
            b2.hit.length,
        ) = (24.6, 2.4, 10, 20);
        (b2.hit.q_start, b2.hit.q_end) = (49, 30);
        let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
            [(0, vec![&a]), (1, vec![&b1, &b2])].into_iter().collect();
        let titles: std::collections::HashMap<u32, Vec<u8>> =
            [(0, b"s0".to_vec()), (1, b"s1".to_vec())]
                .into_iter()
                .collect();
        let columns = Some(SortColumns {
            query_length: 897,
            tblastx: false,
        });
        let mut out = Vec::new();
        write_blastn_description_table(
            &mut out,
            &[0, 1],
            true,
            &hits,
            &titles,
            &[],
            &[],
            false,
            true,
            columns,
            &mut crate::report::query_warnings::no_subject_reads,
            &crate::report::db_subject::HitLoads::default(),
        )
        .unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(
            lines[0],
            format!("{}Score   Total  Query    E     Max", " ".repeat(70))
        );
        assert_eq!(
            lines[1],
            "Sequences producing significant alignments:                          (Bits)   Score  cover   Value  Ident"
        );
        assert_eq!(lines[2], "");
        let label = |name: &str| format!("{name:<68}");
        // 99 of 897 letters (11 %), 27 of 99 identities (27 %).
        assert_eq!(
            lines[3],
            format!("{}  28.1     28.1   11%    0.18   27%       ", label("s0"))
        );
        // 24.6 + 24.6; positions 10-49 (40 of 897: 4 %); the first of the equal best HSPs.
        assert_eq!(
            lines[4],
            format!("{}  24.6     49.2   4%     2.3    64%       ", label("s1"))
        );
        // TBLASTX: the identity over three nucleotides per amino acid (60 of 183: 33 %), after
        // the N column ("NIdent").
        let t = tblastx_hit(2);
        let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
            [(0, vec![&t])].into_iter().collect();
        let titles: std::collections::HashMap<u32, Vec<u8>> =
            [(0, b"t".to_vec())].into_iter().collect();
        let columns = Some(SortColumns {
            query_length: 300,
            tblastx: true,
        });
        let mut out = Vec::new();
        write_blastn_description_table(
            &mut out,
            &[0],
            true,
            &hits,
            &titles,
            &[],
            &[],
            true,
            false,
            columns,
            &mut crate::report::query_warnings::no_subject_reads,
            &crate::report::db_subject::HitLoads::default(),
        )
        .unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert!(
            lines[1].ends_with("(Bits)   Score  cover   Value  NIdent"),
            "{text}"
        );
        // Positions 18-200: 183 of 300 (61 %).
        assert!(
            lines[3].ends_with("  120      120    61%    2e-30  2  33%       "),
            "{text}"
        );
        // The identity is the Seq-align's, not that of the rows: NCBI BLAST+ 2.17.0 `tblastn
        // -query` (the first record of t1 `tn_q.faa`) `-subject tn_s.fna -sorthits 0` shows
        // "94%" for an HSP of "Identities = 284/286" whose engine counted the query's SEG-masked
        // letters as mismatches, and "99%" with `-seg no`; here 270 engine identities of 286.
        let mut masked = tblastx_hit(1);
        (
            masked.hit.num_ident,
            masked.seqalign_num_ident,
            masked.hit.length,
        ) = (284, 270, 286);
        (masked.hit.q_start, masked.hit.q_end) = (1, 286);
        let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
            [(0, vec![&masked])].into_iter().collect();
        let columns = Some(SortColumns {
            query_length: 286,
            tblastx: false,
        });
        let mut out = Vec::new();
        write_blastn_description_table(
            &mut out,
            &[0],
            true,
            &hits,
            &titles,
            &[],
            &[],
            false,
            true,
            columns,
            &mut crate::report::query_warnings::no_subject_reads,
            &crate::report::db_subject::HitLoads::default(),
        )
        .unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(
            text.lines()
                .nth(3)
                .unwrap()
                .ends_with("  100%   2e-30  94%       "),
            "{text}"
        );
    }

    // NCBI c++/src/objtools/align_format/showdefline.cpp:1133-1138,940-944 and
    // align_format_util.cpp:932-938: a last row whose total is wider than every earlier total
    // pads the total column with a negative `size_t`, where NCBI writes spaces without end; the
    // table without `-sorthits` only widens the score column.
    #[test]
    fn sorthits_refuses_a_last_total_wider_than_its_column() {
        let mut a = tblastx_hit(1);
        a.hit.bit_score = 60000.0;
        let mut b = a.clone();
        b.hit.s_idx = 1;
        let mut c = a.clone();
        c.hit.s_idx = 1;
        let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
            [(0, vec![&a]), (1, vec![&b, &c])].into_iter().collect();
        let titles: std::collections::HashMap<u32, Vec<u8>> =
            [(0, b"s0".to_vec()), (1, b"s1".to_vec())]
                .into_iter()
                .collect();
        let columns = Some(SortColumns {
            query_length: 300,
            tblastx: false,
        });
        let table = |order: &[u32], last_row_counted: bool, columns: Option<SortColumns>| {
            let mut out = Vec::new();
            write_blastn_description_table(
                &mut out,
                order,
                last_row_counted,
                &hits,
                &titles,
                &[],
                &[],
                false,
                true,
                columns,
                &mut crate::report::query_warnings::no_subject_reads,
                &crate::report::db_subject::HitLoads::default(),
            )
            .map(|()| out)
        };
        let error = table(&[0, 1], true, columns).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::Unsupported);
        assert_eq!(error.to_string(), TotalScoreRefusal.to_string());
        assert!(is_total_score_refusal(&error.into()));
        assert!(!is_total_score_refusal(
            &io::Error::new(io::ErrorKind::Unsupported, "other").into()
        ));
        assert!(!is_total_score_refusal(&anyhow::anyhow!("other")));
        assert!(table(&[0, 1], true, None).is_ok());
        // Not the last subject of the list, or an earlier total as wide: written.
        assert!(table(&[0, 1], false, columns).is_ok());
        assert!(table(&[1, 1], true, columns).is_ok());
    }

    // NCBI c++/src/objtools/align_format/align_format_util.cpp:4259,4299-4302: BLASTX's row
    // shows the subject's highest bit score and that HSP's E-value (the first of equal
    // maxima), not its first HSP by E-value, and the N column the first HSP's linked-set
    // size. NCBI BLAST+ 2.17.0 `blastx -query tests/fasta/outfmt0/bx_best_hsp_query.fna
    // -subject tests/fasta/outfmt0/bx_best_hsp_subject.faa` (fixture `best_hsp.blastx`):
    // HSPs of 119 bits (Expect(2) = 8e-38), 45.4 and 131 (5e-35); the row is `131     5e-35`.
    #[test]
    fn blastx_description_row_is_the_highest_bit_score() {
        let mut first = tblastx_hit(2);
        first.hit.bit_score = 119.2;
        first.hit.e_value = 8e-38;
        let mut best = tblastx_hit(1);
        best.hit.bit_score = 131.3;
        best.hit.e_value = 5e-35;
        let mut tie = best.clone();
        tie.hit.e_value = 0.5;
        for hit in [&mut first, &mut best, &mut tie] {
            hit.subject_title = Some("pep YOL075C".to_string());
        }
        let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
            [(0, vec![&first, &best, &tie])].into_iter().collect();
        let ids = [DescriptionLabel::title_only(
            b"YOL075C pep YOL075C".to_vec(),
        )];
        for (show_sum_n, end) in [(false, ""), (true, "  2")] {
            let mut out = Vec::new();
            write_subject_summary_table_with_sum_n(
                &mut out,
                &[0],
                true,
                &hits,
                &ids,
                show_sum_n,
                None,
                &mut crate::report::query_warnings::no_subject_reads,
            )
            .unwrap();
            let text = String::from_utf8(out).unwrap();
            let lines: Vec<&str> = text.lines().collect();
            assert_eq!(
                lines[3],
                format!("{:<68}  131     5e-35{end}", "YOL075C pep YOL075C"),
                "{text}"
            );
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452 (quoted at
    // `QueryWarnings`): a query's warning goes to `cerr`, whose tie flushes the report before it.
    // Line RG T3d: BLASTX's report goes through the shared per-query writer
    // (`x_pairwise_par::write`): a query's warning comes after the reports written before it,
    // also from a writer that holds its bytes until it is flushed (LOSAT's buffer of the command
    // line, `cli::NcbiOstream::buffered`); a query whose description table reads an HSP with empty
    // rows stops after its `Length=` (`EmptyRowsStop::RowIsEmpty`): its error comes after the
    // report's last bytes, and nothing is written after it, no text of a later read either
    // (`QueryWarnings::end_run`); the parallel writer (`x_pairwise_par::write` given a pool, which
    // `LOSAT_X_PAIRWISEPAR` does not turn off here) writes the same bytes.
    #[test]
    #[cfg_attr(
        not(feature = "parallel"),
        ignore = "needs threads: a build without the parallel feature has none"
    )]
    fn blastx_query_warnings_follow_the_reports_before_them() {
        use std::sync::{Arc as SyncArc, Mutex};
        #[derive(Clone, Default)]
        struct Log(SyncArc<Mutex<Vec<u8>>>);
        impl Write for Log {
            fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
                self.0.lock().unwrap().extend_from_slice(bytes);
                Ok(bytes.len())
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        /// Holds what is written until it is flushed into the log.
        struct Held(Vec<u8>, Log);
        impl Write for Held {
            fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
                self.0.extend_from_slice(bytes);
                Ok(bytes.len())
            }
            fn flush(&mut self) -> io::Result<()> {
                self.1.write_all(&std::mem::take(&mut self.0))
            }
        }
        let first = tblastx_hit(1);
        let mut empty = tblastx_hit(1);
        empty.hit.q_idx = 2;
        empty.hit.sort_query_end = empty.hit.sort_query_offset;
        let query = |name: &[u8]| BlastpPairwiseQuery {
            query_name: name.to_vec(),
            query_length: 300,
            ungapped_karlin: Default::default(),
            effective_search_space: 1,
            valid: true,
            batch_skipped: false,
        };
        let queries = [query(b"q0"), query(b"q1"), query(b"q2")];
        let report = BlastpPairwiseReport {
            version: String::new(),
            database_name: String::new(),
            database_num_sequences: 1,
            database_total_letters: 400,
            matrix_name: Some("BLOSUM62".into()),
            gap_open: 11,
            gap_extend: 1,
            word_threshold: 12.0,
            window_size: 40,
            gapped_karlin: Default::default(),
            gumbel: None,
            num_descriptions: usize::MAX,
            num_alignments: usize::MAX,
        };
        let options = BlastxPairwiseOptions {
            word_threshold: 12.0,
            gapped: true,
            sum_stats: false,
            num_descriptions: 500,
            num_alignments: 0,
        };
        let posts: Vec<crate::report::diag_log::Posts> = [&b""[..], b"W1\n", b"W2\n"]
            .map(|lines| crate::report::diag_log::Posts::of_lines(lines.to_vec()))
            .into();
        // Query 2's table would read subject 0 (its text `T0`): the stop comes before it.
        let reads = vec![
            Vec::new(),
            Vec::new(),
            vec![crate::report::query_warnings::SubjectRead::of_text(
                0, b"T0\n",
            )],
        ];
        let config = PairwiseConfig {
            program: "blastx".into(),
            ..PairwiseConfig::default()
        };
        let run = |pool: Option<&crate::utils::threading::SearchPool<'_>>| -> (String, String) {
            let log = Log::default();
            let mut sink = log.clone();
            let mut held = Held(Vec::new(), log.clone());
            let mut warnings =
                crate::report::query_warnings::QueryWarnings::of_posts(&posts, None, &mut sink)
                    .with_subject_reads(&reads);
            let error = write_blastx_pairwise_report(
                &[first.clone(), empty.clone()],
                &mut held,
                Some(&mut warnings),
                &config,
                &queries,
                &[DescriptionLabel::title_only(b"s0 title".to_vec())],
                &report,
                &options,
                None,
                &[],
                &[],
                None,
                0,
                None,
                pool,
                &crate::report::db_subject::HitLoads::default(),
            )
            .unwrap_err();
            // The merged stream when the error comes back: the report's last bytes are there
            // (written to the stream, or held by its writer until the flush that comes before the
            // error's post, as `cerr`'s tie flushes `cout`).
            let at_error = [log.0.lock().unwrap().as_slice(), held.0.as_slice()].concat();
            drop(warnings);
            held.flush().unwrap();
            let text = String::from_utf8(log.0.lock().unwrap().clone()).unwrap();
            assert_eq!(
                text.as_bytes(),
                at_error.as_slice(),
                "bytes written after the error (pool {})",
                pool.is_some()
            );
            (text, error.to_string())
        };
        let (serial, error) = run(None);
        assert_eq!(error, "RowIsEmpty");
        let warning = serial.find("W1\n").expect("query 1's warning");
        assert!(
            serial[..warning].contains("Query= q0")
                && serial[..warning].ends_with("Effective search space used: 1\n"),
            "{serial}"
        );
        assert!(
            serial[warning..].starts_with("W1\n\n\nQuery= q1"),
            "{serial}"
        );
        assert!(
            serial.ends_with("W2\n\n\nQuery= q2\n\nLength=300\n"),
            "{serial}"
        );
        crate::utils::threading::x_with_search_pool_local(4, "test", |pool| {
            assert_eq!(run(Some(pool)), (serial.clone(), error.clone()));
            Ok(())
        })
        .unwrap();
    }

    /// The title that the FASTA reader gives a record with the defline `defline`.
    fn reader_title(defline: &[u8]) -> Vec<u8> {
        use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};
        let bytes = [&b">"[..], defline, b"\nACGT\n"].concat();
        let mut source =
            FastaInputSource::from_bytes(&bytes, ReaderConfig::query("BLASTN", false, false));
        read_all(&mut source, &mut |_| Ok(())).unwrap()[0]
            .title
            .clone()
    }

    // NCBI BLAST+ 2.17.0, `blastn -task blastn` outfmt 0 (inventory range RP,
    // `scratch_RP/out/<case>/blastn_f0.out`): the `Query= ` line is the raw title wrapped
    // by `NStr::Wrap` at 68 bytes (a line can end inside a UTF-8 sequence), and the
    // heading is the title wrapped by `s_WrapOutputLine` (60 bytes, then the next white
    // space). `w_<kind>_<k>`: `id `, k letters `a`, one character, `zz tail words`.
    #[test]
    fn query_lines_and_headings_wrap_bytes_as_ncbi() {
        let cjk = &b"\xe6\xbc\xa2"[..];
        let eacute = &b"\xc3\xa9"[..];
        let a = |count: usize| vec![b'a'; count];
        for (defline, query_block) in [
            (b"id desc".to_vec(), b"id desc\n".to_vec()),
            (b"id\tdesc".to_vec(), b"id\n".to_vec()),
            (b"".to_vec(), b"".to_vec()),
            (b"id\xe9x desc".to_vec(), b"id\xe9x desc\n".to_vec()),
            (
                b"id &amp; &lt;b&gt; \"q\" desc".to_vec(),
                b"id &amp; &lt;b&gt; \"q\" desc\n".to_vec(),
            ),
            (b"\x01id desc".to_vec(), b"\x01id desc\n".to_vec()),
            (b"id desc   ".to_vec(), b"id desc\n".to_vec()),
            (b"TPA: id desc".to_vec(), b"TPA: id desc\n".to_vec()),
            (
                [&b"id "[..], &a(59), cjk, b"zz tail words"].concat(),
                [&b"id "[..], &a(59), cjk, b"zz\ntail words\n"].concat(),
            ),
            (
                [&b"id "[..], &a(60), cjk, b"zz tail words"].concat(),
                [&b"id "[..], &a(60), cjk, b"zz\ntail words\n"].concat(),
            ),
            (
                [&b"id "[..], &a(63), cjk, b"zz tail words"].concat(),
                [&b"id\n"[..], &a(63), cjk, b"zz\ntail words\n"].concat(),
            ),
            (
                [&b"id "[..], &a(64), cjk, b"zz tail words"].concat(),
                [&b"id\n"[..], &a(64), cjk, b"z\nz tail words\n"].concat(),
            ),
            (
                [&b"id "[..], &a(65), cjk, b"zz tail words"].concat(),
                [&b"id\n"[..], &a(65), cjk, b"\nzz tail words\n"].concat(),
            ),
            (
                [&b"id "[..], &a(66), cjk, b"zz tail words"].concat(),
                [&b"id\n"[..], &a(66), b"\xe6\xbc\n\xa2zz tail words\n"].concat(),
            ),
            (
                [&b"id "[..], &a(67), cjk, b"zz tail words"].concat(),
                [&b"id\n"[..], &a(67), b"\xe6\n\xbc\xa2zz tail words\n"].concat(),
            ),
            (
                [&b"id "[..], &a(66), eacute, b"zz tail words"].concat(),
                [&b"id\n"[..], &a(66), eacute, b"\nzz tail words\n"].concat(),
            ),
            (
                [&b"id "[..], &a(67), eacute, b"zz tail words"].concat(),
                [&b"id\n"[..], &a(67), b"\xc3\n\xa9zz tail words\n"].concat(),
            ),
            (
                [&b"id "[..], &a(66), b"Zzz tail words"].concat(),
                [&b"id\n"[..], &a(66), b"Zz\nz tail words\n"].concat(),
            ),
        ] {
            let title = reader_title(&defline);
            let mut out = Vec::new();
            write_blastp_query_header(&mut out, &title, 330).unwrap();
            assert_eq!(
                out,
                [&b"Query= "[..], &query_block, b"\nLength=330\n"].concat(),
                "{defline:?}"
            );
        }
        for (middle, heading) in [
            (
                [&a(59), cjk].concat(),
                [&a(59), cjk, b"zz \ntail words"].concat(),
            ),
            (
                [&a(66), cjk].concat(),
                [&a(66), cjk, b"zz \ntail words"].concat(),
            ),
            (
                [&a(67), eacute].concat(),
                [&a(67), eacute, b"zz \ntail words"].concat(),
            ),
            (
                [&a(66)[..], b"Z"].concat(),
                [&a(66)[..], b"Zzz \ntail words"].concat(),
            ),
        ] {
            let title = [&b"id "[..], &middle, b"zz tail words"].concat();
            let mut out = Vec::new();
            write_subject_header(&mut out, &title, None, Some(450)).unwrap();
            assert_eq!(
                out,
                [&b"> id "[..], &heading, b"\nLength=450\n\n"].concat(),
                "{title:?}"
            );
        }
    }

    // NCBI BLAST+ 2.17.0 outfmt 0 for a subject whose title NCBI cannot decode (inventory
    // range RP, `scratch_RP/out/u06_undef81/tblastn_f0.out`: `id \x81 x`, three HSPs):
    // the description is `Unknown`, and each HSP is replaced by one line, without a heading.
    #[test]
    fn a_title_of_unknown_encoding_is_unknown_and_its_alignments_are_skipped() {
        let first = tblastx_hit(1);
        let mut second = tblastx_hit(1);
        second.hit.bit_score = 13.9;
        let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
            [(0, vec![&first, &second])].into_iter().collect();
        let titles: std::collections::HashMap<u32, Vec<u8>> =
            [(0, b"id \x81 x".to_vec())].into_iter().collect();
        let mut out = Vec::new();
        write_blastn_description_table(
            &mut out,
            &[0],
            true,
            &hits,
            &titles,
            &[],
            &[],
            false,
            false,
            None,
            &mut crate::report::query_warnings::no_subject_reads,
            &crate::report::db_subject::HitLoads::default(),
        )
        .unwrap();
        let row = out.split(|&byte| byte == b'\n').nth(3).unwrap().to_vec();
        assert!(row.starts_with(&[&b"Unknown"[..], &[b' '; 61], b"  120"].concat()));
        let mut out = Vec::new();
        write_subject_alignments_part(
            &mut out,
            None,
            0,
            &[],
            &titles[&0],
            false,
            &hits[&0],
            &[0, 1],
            0..hits[&0].len(),
            |_, _| panic!("no HSP is written"),
        )
        .unwrap();
        assert_eq!(
            out,
            b"Sequence with id Subject_1 no longer exists in database...alignment skipped\n\
              Sequence with id Subject_1 no longer exists in database...alignment skipped\n"
        );
    }

    // DBID pfoff-3 (DECISIONS 62 (13); NCBI BLAST+ 2.17.0, probe `R/pfoff/p3`, `tx.U2.b1.dflt.0`):
    // a subject whose load throws where its batch's prefetch is off has the description line of
    // its Seq-align id and `Unknown`, and in outfmt 0 each of its HSPs one line that names the
    // id's accession without its version, with no heading (`DatabaseHeading::Thrown`).
    #[test]
    fn a_subject_whose_load_throws_is_unknown_and_its_alignments_are_skipped() {
        use crate::report::db_subject::{
            AfterThrow, FoundInstead, HitLoads, LookupFailure, ParsedSubject, RunThrow,
        };
        let subject = |accession: &[u8], throws: bool| {
            let subject = ParsedSubject::new(crate::seqid::DeflineSet {
                deflines: vec![crate::seqid::BlastDefline {
                    title: Some(b"ntitle".to_vec()),
                    ids: vec![crate::seqid::ParsedSeqId::Other(
                        crate::seqid::model::TextseqId {
                            accession: Some(accession.to_vec()),
                            version: Some(1),
                            ..Default::default()
                        },
                    )],
                    ..Default::default()
                }],
            })
            .unwrap();
            Some(Arc::new(if throws {
                subject.found_instead(Some(FoundInstead::Failed(LookupFailure::ncbi(
                    crate::seqdb::tax_files::NO_ACCESSION_INFO,
                ))))
            } else {
                subject
            }))
        };
        // 0 loads, 1 throws.
        let parsed = vec![subject(b"NM_100000", false), subject(b"NM_100001", true)];
        let loads = HitLoads::of(
            Some(RunThrow::every_row(AfterThrow::SkipsAlignment)),
            [(true, vec![Some(0), Some(1)])],
            |key| parsed[key].as_deref(),
            |_, _| 0,
        );
        assert!(loads.loaded(0));
        assert!(loads.thrown(1).is_some());
        let mut second = tblastx_hit(1);
        second.hit.s_idx = 1;
        let mut third = second.clone();
        third.hit.bit_score = 13.9;
        let first = tblastx_hit(1);
        let hits: std::collections::HashMap<u32, Vec<&PairwiseHit>> =
            [(0, vec![&first]), (1, vec![&second, &third])]
                .into_iter()
                .collect();
        let titles: std::collections::HashMap<u32, Vec<u8>> =
            [(0, b"ntitle".to_vec()), (1, b"ntitle".to_vec())]
                .into_iter()
                .collect();
        let mut out = Vec::new();
        write_blastn_description_table(
            &mut out,
            &[0, 1],
            true,
            &hits,
            &titles,
            &parsed,
            &[],
            false,
            false,
            None,
            &mut crate::report::query_warnings::no_subject_reads,
            &loads,
        )
        .unwrap();
        let rows: Vec<&[u8]> = out.split(|&byte| byte == b'\n').skip(3).take(2).collect();
        assert!(
            rows[0].starts_with(b"NM_100000.1 ntitle "),
            "{:?}",
            String::from_utf8_lossy(rows[0])
        );
        assert!(
            rows[1].starts_with(b"NM_100001.1 Unknown "),
            "{:?}",
            String::from_utf8_lossy(rows[1])
        );
        let mut out = Vec::new();
        write_database_subject_alignments_part(
            &mut out,
            None,
            DatabaseHeading::loaded(parsed[1].as_deref().unwrap(), false, 2, true),
            false,
            &hits[&1],
            &[1, 2],
            0..2,
            |_, _| panic!("no HSP is written"),
        )
        .unwrap();
        assert_eq!(
            out,
            b"Sequence with id NM_100001 no longer exists in database...alignment skipped\n\
              Sequence with id NM_100001 no longer exists in database...alignment skipped\n"
        );
    }

    /// A BLASTP pairwise report written a window of queries at a time (the rows of a caller that
    /// makes them window by window) is the report written at once: the prolog with the first
    /// window, the epilog with the last, the query indices global. So is a psiblast round's report
    /// (`psiblast_round`, the prolog with the first window, no epilog), also in two parts.
    #[test]
    fn a_blastp_report_written_in_windows_is_the_report_written_at_once() {
        let karlin = KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            ..KarlinParams::default()
        };
        let query = |length: usize| BlastpPairwiseQuery {
            query_name: Vec::new(),
            query_length: length,
            ungapped_karlin: karlin,
            effective_search_space: 1_000_000,
            valid: true,
            batch_skipped: false,
        };
        let queries = [query(30), query(40), query(50), query(60)];
        let query_titles: Vec<Arc<[u8]>> = (0..4)
            .map(|i| Arc::from(format!("Query_{i}").into_bytes().into_boxed_slice()))
            .collect();
        let subject_titles: Vec<Arc<[u8]>> = vec![Arc::from(&b"gnl|CDD|1 profile"[..])];
        let parsed: Vec<Option<Arc<crate::report::db_subject::ParsedSubject>>> = vec![None];
        let hit = |q_idx: u32| {
            let mut pairwise = PairwiseHit::from(Hit {
                identity: 100.0,
                length: 10,
                mismatch: 0,
                gapopen: 0,
                q_start: 1,
                q_end: 10,
                s_start: 1,
                s_end: 10,
                e_value: 1e-5,
                bit_score: 25.0,
                num_ident: 10,
                query_frame: 0,
                query_length: 30,
                q_idx,
                s_idx: 0,
                raw_score: 50,
                sort_query_offset: 0,
                sort_query_end: 0,
                sort_subject_offset: 0,
                sort_subject_end: 0,
                has_sort_offsets: false,
                gap_info: None,
                num_positives: 10,
            });
            pairwise.query_seq = Some("MKVLAAGIVG".to_string());
            pairwise.subject_seq = Some("MKVLAAGIVG".to_string());
            pairwise.subject_length = Some(180);
            pairwise.subject_title = Some("profile".to_string());
            pairwise
        };
        // Queries 0 and 2 have no hit.
        let hits = [hit(1), hit(3)];
        let report = BlastpPairwiseReport {
            version: "BLASTP 2.17.0+".to_string(),
            database_name: "db".to_string(),
            database_num_sequences: 1,
            database_total_letters: 180,
            matrix_name: Some("BLOSUM62".to_string()),
            gap_open: 11,
            gap_extend: 1,
            word_threshold: 0.0,
            window_size: 40,
            gapped_karlin: karlin,
            gumbel: None,
            num_descriptions: 500,
            num_alignments: 250,
        };
        let config = PairwiseConfig {
            line_length: 60,
            show_frame: false,
            program: "rpsblast".to_string(),
            protein_matrix: ScoringMatrix::Blosum62,
            sort: Default::default(),
            view: PairwiseView::Pairwise,
        };
        let gapped_k = [0.041, 0.047, 0.041, 0.041];
        let write = |range: std::ops::Range<usize>, rows: &[PairwiseHit], out: &mut Vec<u8>| {
            write_blastp_pairwise_window(
                rows,
                range,
                out,
                &config,
                &queries,
                &query_titles,
                &subject_titles,
                &parsed,
                &[],
                &report,
                Some(&gapped_k),
                None,
                None,
                None,
                true,
                true,
                None,
                None,
                &crate::report::db_subject::HitLoads::default(),
            )
            .unwrap();
        };
        let mut at_once = Vec::new();
        write(0..4, &hits, &mut at_once);
        assert!(String::from_utf8_lossy(&at_once).contains("Query_3"));
        for windows in [
            vec![0..1, 1..4],
            vec![0..2, 2..3, 3..4],
            vec![0..1, 1..2, 2..3, 3..4],
        ] {
            let mut windowed = Vec::new();
            for window in windows {
                let rows: Vec<PairwiseHit> = hits
                    .iter()
                    .filter(|hit| window.contains(&(hit.hit.q_idx as usize)))
                    .cloned()
                    .collect();
                write(window, &rows, &mut windowed);
            }
            assert_eq!(
                String::from_utf8_lossy(&windowed),
                String::from_utf8_lossy(&at_once)
            );
        }
        let round_config = PairwiseConfig {
            program: "psiblast".to_string(),
            ..config.clone()
        };
        for found_before in [None, Some(&[true][..])] {
            let round = |range: std::ops::Range<usize>, rows: &[PairwiseHit], out: &mut Vec<u8>| {
                psiblast_round::write_psiblast_round_report(
                    rows,
                    range,
                    out,
                    &round_config,
                    &queries,
                    &query_titles,
                    &subject_titles,
                    &parsed,
                    &[],
                    &crate::report::db_subject::HitLoads::default(),
                    &report,
                    None,
                    None,
                    None,
                    true,
                    psiblast_round::ReportRound {
                        number: 2,
                        found_before,
                    },
                )
                .unwrap();
            };
            let mut at_once = Vec::new();
            round(0..4, &hits, &mut at_once);
            let text = String::from_utf8_lossy(&at_once).into_owned();
            assert!(text.contains("Results from round 2") && text.contains("Query_3"));
            for windows in [vec![0..1, 1..4], vec![0..2, 2..3, 3..4]] {
                let mut windowed = Vec::new();
                for window in windows {
                    let rows: Vec<PairwiseHit> = hits
                        .iter()
                        .filter(|hit| window.contains(&(hit.hit.q_idx as usize)))
                        .cloned()
                        .collect();
                    round(window, &rows, &mut windowed);
                }
                assert_eq!(String::from_utf8_lossy(&windowed), text);
            }
        }
    }
}
