//! Pairwise alignment output (outfmt 0)
//!
//! This module implements the traditional BLAST pairwise alignment output format.
//!
//! Reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp

use super::outfmt6::{
    format_bitscore_ncbi, format_evalue_ncbi, format_total_bitscore_ncbi, truncate_spaces,
};
use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::common::Hit;
use crate::config::ScoringMatrix;
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
}

impl Default for PairwiseConfig {
    fn default() -> Self {
        Self {
            line_length: DEFAULT_LINE_LENGTH,
            show_frame: true,
            program: "tblastx".to_string(),
            protein_matrix: ScoringMatrix::Blosum62,
            sort: Default::default(),
        }
    }
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
    /// residues (`write_sequence_row`).
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
    pub query_name: String,
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
    pub matrix_name: String,
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
    /// `number_seqs`.
    pub num_sequences: u64,
    /// `total_length`.
    pub total_letters: u64,
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
    /// The sum of their numbers of sequences.
    pub num_sequences: u64,
    /// The sum of their total lengths.
    pub total_letters: u64,
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
                .map(|database| database.num_sequences)
                .sum(),
            total_letters: databases
                .iter()
                .map(|database| database.total_letters)
                .sum(),
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
    let mut do_wrap = false;
    for (i, &c) in title.iter().enumerate() {
        if i > 0 && i % 60 == 0 {
            do_wrap = true;
        }
        writer.write_all(&[c])?;
        if do_wrap && (c.is_ascii_whitespace() || c == 11) {
            writer.write_all(b"\n")?;
            do_wrap = false;
        }
    }
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
#[allow(clippy::too_many_arguments)]
fn write_subject_alignments<W: Write>(
    writer: &mut W,
    probe: Option<&mut FormatProbe<'_>>,
    s_idx: u32,
    title: &[u8],
    protein: bool,
    hits: &[&PairwiseHit],
    hsp_indices: &[HspIndex],
    write_hsp: impl FnMut(&mut W, &PairwiseHit) -> io::Result<()>,
) -> io::Result<()> {
    write_subject_alignments_part(
        writer,
        probe,
        s_idx,
        title,
        protein,
        hits,
        hsp_indices,
        0..hits.len(),
        write_hsp,
    )
}

/// `write_subject_alignments` for the HSPs `part` of the subject's (the heading with the first);
/// the parts of `0..hits.len()` in order write its bytes (EXPERIMENT LOSAT_X_PAIRWISEPAR's pieces,
/// `x_pairwise_par`).
#[allow(clippy::too_many_arguments)]
fn write_subject_alignments_part<W: Write>(
    writer: &mut W,
    mut probe: Option<&mut FormatProbe<'_>>,
    s_idx: u32,
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
            Err(super::defline::UnknownEncoding) => writeln!(
                writer,
                "Sequence with id Subject_{} no longer exists in database...alignment skipped",
                u64::from(s_idx) + 1
            )?,
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
/// follows it. A stored title is a VisibleString (bytes from ' ' to '~'), so the title
/// that `GenerateDefline` makes of an empty one (for a protein, `unnamed protein product`)
/// has an encoding that is never unknown. `protein` is the database's molecule type (BLASTP
/// and BLASTN share this heading).
fn write_database_subject_alignments<'h, W: Write>(
    writer: &mut W,
    probe: Option<&mut FormatProbe<'_>>,
    title: impl Into<DatabaseHeading<'h>>,
    protein: bool,
    hits: &[&PairwiseHit],
    hsp_indices: &[HspIndex],
    write_hsp: impl FnMut(&mut W, &PairwiseHit) -> io::Result<()>,
) -> io::Result<()> {
    write_database_subject_alignments_part(
        writer,
        probe,
        title,
        protein,
        hits,
        hsp_indices,
        0..hits.len(),
        write_hsp,
    )
}

/// `write_database_subject_alignments` for the HSPs `part` of the subject's (the heading with
/// the first), as `write_subject_alignments_part`.
#[allow(clippy::too_many_arguments)]
fn write_database_subject_alignments_part<'h, W: Write>(
    writer: &mut W,
    mut probe: Option<&mut FormatProbe<'_>>,
    title: impl Into<DatabaseHeading<'h>>,
    protein: bool,
    hits: &[&PairwiseHit],
    hsp_indices: &[HspIndex],
    part: std::ops::Range<usize>,
    mut write_hsp: impl FnMut(&mut W, &PairwiseHit) -> io::Result<()>,
) -> io::Result<()> {
    let first_index = hsp_indices[0];
    let title = title.into();
    let skipped = title.skipped_id(protein);
    if part.start == 0 {
        if let Some(probe) = probe.as_mut() {
            writer.flush()?;
            probe.subject_begin(first_index);
        }
        if skipped.is_none() {
            write_database_subject_heading(writer, title, protein, hits[0].subject_length)?;
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

/// NCBI's line in place of an alignment whose heading it cannot make (`DatabaseHeading::skipped_id`).
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
}

impl<'a> DatabaseHeading<'a> {
    /// The Seq-id text that outfmt 0 writes in place of the subject's alignments ("Sequence with
    /// id <id> no longer exists in database...alignment skipped", showalign.cpp:1976-1981, quoted
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
}

impl<'a> From<&'a [u8]> for DatabaseHeading<'a> {
    fn from(title: &'a [u8]) -> Self {
        DatabaseHeading::Title(title, &super::defline::UNNAMED_SOURCE)
    }
}

impl<'a> From<&'a Vec<u8>> for DatabaseHeading<'a> {
    fn from(title: &'a Vec<u8>) -> Self {
        DatabaseHeading::Title(title, &super::defline::UNNAMED_SOURCE)
    }
}

impl<'a> From<&'a Arc<[u8]>> for DatabaseHeading<'a> {
    fn from(title: &'a Arc<[u8]>) -> Self {
        DatabaseHeading::Title(title, &super::defline::UNNAMED_SOURCE)
    }
}

/// The source of the subject `s_idx` without parsed Seq-ids among `sources`, the sources of a
/// report's database sequences by subject index (`report::defline::source_facts`; empty where
/// none has one).
pub fn stored_source(
    sources: &[super::defline::SourceTaxname],
    s_idx: u32,
) -> &super::defline::SourceTaxname {
    sources
        .get(s_idx as usize)
        .unwrap_or(&super::defline::UNNAMED_SOURCE)
}

/// The heading of a database sequence's alignments (`write_database_subject_alignments`):
/// `>` and the stored title wrapped as `s_WrapOutputLine` wraps it (a sequence with parsed
/// Seq-ids: the line of each definition line, `ParsedSubject::heading`), then `Length=` when the
/// length is known, and one blank line. BLASTX writes it with its own HSP blocks.
fn write_database_subject_heading<'h, W: Write>(
    writer: &mut W,
    title: impl Into<DatabaseHeading<'h>>,
    protein: bool,
    length: Option<usize>,
) -> io::Result<()> {
    let (title, source) = match title.into() {
        DatabaseHeading::Title(title, source) => (title, source),
        DatabaseHeading::Parsed(subject, every_line) => {
            writer.write_all(&subject.heading(protein, every_line))?;
            if let Some(length) = length {
                writeln!(writer, "Length={length}")?;
            }
            return writeln!(writer);
        }
    };
    let generated;
    let heading: &[u8] = if title.is_empty() {
        // NCBI makes an empty title from the Bioseq's source (`showalign.cpp:2250-2255`, quoted at
        // `ParsedSubject::generated_title_fetched`; LOSAT refuses a source that it refuses before
        // the report, `StoredTitle::check_shown_source`).
        let facts = super::defline::BioseqFacts::of_source(source.name().ok().flatten());
        generated = super::defline::generate_bioseq_defline(title, protein, false, &facts)
            .map(|defline| defline.text)
            .unwrap_or_default();
        &generated
    } else {
        title
    };
    writer.write_all(b">")?;
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
    let mut do_wrap = false;
    for (i, &c) in heading.iter().enumerate() {
        if i > 0 && i % 60 == 0 {
            do_wrap = true;
        }
        writer.write_all(&[c])?;
        if do_wrap && (c.is_ascii_whitespace() || c == 11) {
            writer.write_all(b"\n")?;
            do_wrap = false;
        }
    }
    writeln!(writer)?;
    if let Some(length) = length {
        writeln!(writer, "Length={length}")?;
    }
    writeln!(writer)
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

/// Write alignment with actual sequences
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
    let max_id_len = "Query".len().max("Sbjct".len());
    // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1570-1573
    // ```c++
    //     if (m_Program == "blastn" || m_Program == "megablast") {
    //         display.SetMiddleLineStyle(CDisplaySeqalign::eBar);
    //         display.SetAlignType(CDisplaySeqalign::eNuc);
    // ```
    let nucleotide = config.program == "blastn";
    // A minus-strand subject is shown with decreasing coordinates (the query is always
    // on its plus strand, showalign.cpp:1852-1858).
    let subject_step: isize = if hit.s_start > hit.s_end { -1 } else { 1 };

    let mut query_row = RowCursor::new(hit.q_start as isize, 1, 1, line_len);
    let mut subject_row = RowCursor::new(hit.s_start as isize, subject_step, 1, line_len);
    let mut offset = 0;

    while offset < q_chars.len() {
        let end = (offset + line_len).min(q_chars.len());

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
        let middle: String = q_chars[offset..end]
            .iter()
            .zip(s_chars[offset..end].iter())
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
            })
            .collect();

        write_sequence_row(
            writer,
            "Query",
            &mut query_row,
            &q_chars[offset..end],
            0,
            max_start_len,
            |writer| write_chars(writer, &q_chars[offset..end]),
        )?;
        write_spaces(writer, max_id_len + 2 + max_start_len + 2)?;
        writeln!(writer, "{}", middle)?;
        write_sequence_row(
            writer,
            "Sbjct",
            &mut subject_row,
            &s_chars[offset..end],
            stored_gaps_between(subject_stored_gaps, offset, end),
            max_start_len,
            |writer| write_chars(writer, &s_chars[offset..end]),
        )?;
        writeln!(writer)?;

        offset = end;
    }

    Ok(())
}

/// Writes one row of one alignment chunk and advances `next`, the position of the row's
/// next residue, by `step` per residue.
///
/// NCBI reference: c++/src/objtools/align_format/showalign.cpp:1598-1626
/// ```c++
///     int start = alnRoInfo->seqStarts[row].front() + 1;  //+1 for 1 based
///     int end = alnRoInfo->seqStops[row].front() + 1;
///     ...
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
///     ...
///     CAlignFormatUtil::AddSpace(out, k_SeqStopMargin);
///
///      //not to display stop number for empty row in the middle
///     if (!(j > 0 && end == prev_stop)
///         && !(j == 0 && start == 1 && end == 1)) {
///         out << end;
///     }
///     out<<"\n";
/// ```
/// A later chunk without residues ends where the previous chunk ended
/// (`end == prev_stop`). `AddSpace` takes a `size_t`, so an omitted start gives
/// `maxStartLen + 2` spaces. Every program's rows are written here: `columns` are the chunk's
/// letters (`-` a gap, or of `stored_gaps` of them a stored gap character, a residue) and
/// `write_letters` writes them. They are the rows of letters (BLASTN, BLASTP, RPSBLAST, PSIBLAST,
/// and the protein rows of BLASTX and TBLASTN) and the rows translated from nucleotides
/// (BLASTX's query, TBLASTN's subject, both rows of TBLASTX), whose chunk of residues runs over
/// three nucleotides per residue (`RowCursor`). A translated row's first chunk spans a codon, so
/// its start and end differ and it always shows them; a first chunk of one residue at 1 of a
/// row of letters shows none (the `-line_length 1` rows of an HSP from a sequence's first
/// letter). The end is that of NCBI's first screen (`RowCursor::first_screen`): the chunk's
/// last letter, but for a translated row whose screen of `-line_length` 1431655766 or more
/// wraps, where only the end moves (the letters are the whole row's).
fn write_sequence_row<W: Write, T: Copy + PartialEq + From<u8>>(
    writer: &mut W,
    label: &str,
    row: &mut RowCursor,
    columns: &[T],
    stored_gaps: usize,
    width: usize,
    write_letters: impl FnOnce(&mut W) -> io::Result<()>,
) -> io::Result<()> {
    let count = residue_count(columns, stored_gaps) as isize;
    let first_chunk = row.first;
    let start = row.next;
    let in_screen = match first_chunk.then(|| row.first_screen(columns)).flatten() {
        Some(in_screen) => in_screen,
        None => row.letters * count,
    };
    let end = start + row.step * (in_screen - 1);
    let shown = !(!first_chunk && count == 0 || first_chunk && start == 1 && end == 1);
    write!(writer, "{label}")?;
    write_spaces(writer, 2)?;
    if shown {
        write!(writer, "{start}")?;
        write_spaces(writer, width + 2 - digit_count(start.unsigned_abs()))?;
    } else {
        write_spaces(writer, width + 2)?;
    }
    write_letters(writer)?;
    write_spaces(writer, 2)?;
    if shown {
        write!(writer, "{end}")?;
    }
    writeln!(writer)?;
    row.next += row.step * row.letters * count;
    row.first = false;
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1281-1283
// ```c++
//         m_AV->GetWholeAlnSeqString(row, sequence[row], &insertAlnStart[row],
//                                    &insertStart[row], &insertLength[row],
//                                    (int)m_LineLen, &seqStarts[row], &seqStops[row]);
// ```
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:199-206,227-230,305,358
// ```c++
// string& CAlnVec::GetWholeAlnSeqString(TNumrow       row,
//                                       string&       buffer,
//                                       TSeqPosList * insert_aln_starts,
//                                       TSeqPosList * insert_starts,
//                                       TSeqPosList * insert_lens,
//                                       unsigned int  scrn_width,
//                                       TSeqPosList * scrn_lefts,
//                                       TSeqPosList * scrn_rights) const
//     ...
//     bool plus     = IsPositiveStrand(row);
//     int  width    = GetWidth(row);
//
//     scrn_width *= width;
//     ...
//                         scrn_lft_seq_pos = plus ? start : stop;
//     ...
//                     scrn_rgt_seq_pos = plus ? stop : start;
// ```
/// Where one row of an alignment is (`write_sequence_row`): the coordinate of its next residue,
/// the direction of its coordinates (1, or -1 on a minus strand or frame), the letters of its
/// sequence per residue (`GetWidth`: 1, or 3 for a row translated from nucleotides, whose
/// screen of `m_LineLen` residues spans three times as many nucleotides), NCBI's screen in
/// letters of its sequence (`scrn_width` after `*= width`) and whether its first chunk is still
/// to be written (NCBI's `j == 0`).
struct RowCursor {
    next: isize,
    step: isize,
    letters: isize,
    screen: u32,
    first: bool,
}

impl RowCursor {
    /// The row from `next` in chunks of `line_length` columns. `(int)m_LineLen` comes in as an
    /// `unsigned int` (`-line_length` is at most 2^31 - 1) and `*= width` wraps it modulo 2^32,
    /// an unsigned product that C++ defines: a translated row's screen is `3 * line_length`
    /// nucleotides up to 1431655765 (2^32 - 1), and `3 * line_length - 2^32` (2, 5, ... up to
    /// 2147483645) from 1431655766 on. A row of letters' screen is `line_length`.
    fn new(next: isize, step: isize, letters: isize, line_length: usize) -> Self {
        Self {
            next,
            step,
            letters,
            screen: (line_length as u32).wrapping_mul(letters as u32),
            first: true,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:310-343
    // ```c++
    //                     // previous scrns
    //                     nscrns = (aln_pos - scrn_pos) / scrn_width;
    //                     for (int i = 0; i < nscrns; i++) {
    //                         scrn_lefts->push_back(scrn_lft_seq_pos);
    //                         scrn_rights->push_back(scrn_rgt_seq_pos);
    //     ...
    //                     // current scrns
    //                     nscrns = (aln_pos + len - scrn_pos) / scrn_width;
    //                     curr_pos = aln_pos;
    //                     for (int i = 0; i < nscrns; i++) {
    //                         delta = (plus ?
    //                                  scrn_width - (curr_pos - scrn_pos) :
    //                                  curr_pos - scrn_pos - scrn_width);
    //
    //                         scrn_lefts->push_back(scrn_lft_seq_pos);
    //                         if (plus ?
    //                             scrn_lft_seq_pos < start :
    //                             scrn_lft_seq_pos > stop) {
    //                             scrn_lft_seq_pos = (plus ? start : stop) +
    //                                 delta;
    //                             scrn_rgt_seq_pos = scrn_lft_seq_pos +
    //                                 (plus ? -1 : 1);
    //                         } else {
    //                             scrn_rgt_seq_pos = scrn_lft_seq_pos + (plus ? -1 : 1)
    //                                 + delta;
    //                             scrn_lft_seq_pos += delta;
    //                         }
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1598-1599
    // ```c++
    //     int start = alnRoInfo->seqStarts[row].front() + 1;  //+1 for 1 based
    //     int end = alnRoInfo->seqStops[row].front() + 1;
    // ```
    /// The letters of the row's sequence in NCBI's first screen when the screen is narrower
    /// than the first chunk of `columns` (`None` otherwise: the chunk ends at its last letter).
    /// A chunk prints the right end of its screen (`seqStops.front()`); a wrapped screen makes
    /// the screens narrower than the chunks, of which the row then has one (a chunk of
    /// 1431655766 or more residues would span more than 2^32 nucleotides). The screen covers
    /// the alignment's first `scrn_width` letters of the row, a column counting `width` of them
    /// (`aln_pos += len`, `len = seg_len * width`, residues and gaps alike), and its right end
    /// is the row's last letter there: within a segment the end of `delta` letters from its
    /// start, in a gap the previous segment's stop (`scrn_rgt_seq_pos = plus ? stop : start`).
    /// A row's letters are contiguous across its gaps, so that end is the start moved by the
    /// letters in the screen less one: whole residues of the first `scrn_width / width`
    /// columns and `scrn_width % width` letters of the next when it is a residue. Without any
    /// (a row that begins with gaps, which no HSP does), the screen's right is its left (lines
    /// 304-308 and 312-314: one letter). Only translated rows' screens wrap, and their `-` is
    /// a gap: a stored gap character comes only in protein sequences.
    fn first_screen<T: Copy + PartialEq + From<u8>>(&self, columns: &[T]) -> Option<isize> {
        let screen = u64::from(self.screen);
        let width = self.letters as u64;
        if screen >= width * columns.len() as u64 {
            return None;
        }
        let whole = (screen / width) as usize;
        let part = (screen % width) as usize;
        let partial = if part > 0 && columns[whole] != T::from(b'-') {
            part
        } else {
            0
        };
        let letters = self.letters as usize * residue_count(&columns[..whole], 0) + partial;
        Some(letters.max(1) as isize)
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:302-305
// ```c++
//                 // take care of coords if necessary
//                 if (record_coords) {
//                     if (scrn_lft_seq_pos < 0) {
//                         scrn_lft_seq_pos = plus ? start : stop;
// ```
/// The residues of a chunk of a row, the segments of the Dense-seg that its coordinates count
/// (a gap segment has none): its letters other than `-` and its `stored_gaps` columns whose `-`
/// is a stored gap character (a residue of the Dense-seg).
fn residue_count<T: Copy + PartialEq + From<u8>>(letters: &[T], stored_gaps: usize) -> usize {
    let gap = T::from(b'-');
    letters.iter().filter(|&&letter| letter != gap).count() + stored_gaps
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2491-2493
// ```c++
//     _ASSERT((int)sequence.size() > start);
//     list<CRange<int> > actualSeqloc;
//     string actualSeq = sequence.substr(start, len);
// ```
/// Writes the letters of a chunk of a row kept as `char`s (`write_alignment_with_sequences`).
fn write_chars<W: Write>(writer: &mut W, letters: &[char]) -> io::Result<()> {
    for &letter in letters {
        write!(writer, "{letter}")?;
    }
    Ok(())
}

/// The number of the sorted columns `stored_gaps` in `start..end`.
fn stored_gaps_between(stored_gaps: &[u32], start: usize, end: usize) -> usize {
    let below = |column: usize| stored_gaps.partition_point(|&c| (c as usize) < column);
    below(end) - below(start)
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

#[allow(clippy::too_many_arguments)]
fn write_subject_summary_table_with_sum_n<W: Write>(
    writer: &mut W,
    subject_order: &[u32],
    last_row_counted: bool,
    subject_hits: &std::collections::HashMap<u32, Vec<&PairwiseHit>>,
    subject_ids: &[Arc<str>],
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
            let subject_id = subject_ids
                .get(s_idx as usize)
                .map(|id| id.as_ref())
                .unwrap_or("unknown");
            let mut label = subject_id.as_bytes().to_vec();
            let title = subject_hits
                .get(&s_idx)
                .and_then(|hits| hits.first())
                .and_then(|hit| hit.subject_title.as_deref());
            if let Some(title) = title {
                label.push(b' ');
                label.extend_from_slice(title.as_bytes());
            }
            label
        },
        read,
    )
}

#[inline]
fn digit_count(value: usize) -> usize {
    value.to_string().len()
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
    for _ in 0..count {
        writer.write_all(b" ")?;
    }
    Ok(())
}

#[inline]
fn format_count_with_commas(value: i64) -> String {
    let digits = value.to_string();
    let mut formatted = String::with_capacity(digits.len() + digits.len() / 3);
    let len = digits.len();
    for (index, ch) in digits.chars().enumerate() {
        if index > 0 && (len - index).is_multiple_of(3) {
            formatted.push(',');
        }
        formatted.push(ch);
    }
    formatted
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

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/format/blastfmtutil.cpp:73-111
// ```c
// if (m_Program == "psiblast" || m_Program == "blastp") {
//     CBlastFormatUtil::BlastPrintReference(..., CReference::eCompBasedStats, ...);
// }
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/api/version.cpp:46-69
// ```c
// "Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schaffer, ..."
// "Alejandro A. Schaffer, L. Aravind, Thomas L. Madden, ..."
// ```
fn write_blastp_pairwise_intro<W: Write>(
    writer: &mut W,
    version: &str,
    program: &str,
) -> io::Result<()> {
    if program == "psiblast" {
        return write_psiblast_pairwise_intro(writer, version);
    }
    if program == "rpsblast" {
        return write_rpsblast_pairwise_intro(writer, version);
    }
    writeln!(writer, "BLASTP {}", version)?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A."
    )?;
    writeln!(
        writer,
        "Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J."
    )?;
    writeln!(
        writer,
        "Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of"
    )?;
    writeln!(
        writer,
        "protein database search programs\", Nucleic Acids Res. 25:3389-3402."
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Reference for composition-based statistics: Alejandro A. Schaffer,"
    )?;
    writeln!(
        writer,
        "L. Aravind, Thomas L. Madden, Sergei Shavirin, John L. Spouge, Yuri"
    )?;
    writeln!(
        writer,
        "I. Wolf, Eugene V. Koonin, and Stephen F. Altschul (2001),"
    )?;
    writeln!(
        writer,
        "\"Improving the accuracy of PSI-BLAST protein database searches with"
    )?;
    writeln!(
        writer,
        "composition-based statistics and other refinements\", Nucleic Acids"
    )?;
    writeln!(writer, "Res. 29:2994-3005.")?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)?;
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:386-405
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
//     ...
//     else
//         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
//                                           m_Outfile);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:416-426 (quoted below:
// the references of compositional matrix adjustment and of composition-based statistics are
// psiblast's, deltablast's and blastp's only)
/// The prolog of rpsblast's pairwise report: the program line and the Gapped BLAST reference
/// only.
fn write_rpsblast_pairwise_intro<W: Write>(writer: &mut W, version: &str) -> io::Result<()> {
    writeln!(writer, "RPSBLAST {}", version)?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A."
    )?;
    writeln!(
        writer,
        "Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J."
    )?;
    writeln!(
        writer,
        "Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of"
    )?;
    writeln!(
        writer,
        "protein database search programs\", Nucleic Acids Res. 25:3389-3402."
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)?;
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:416-426
// ```c++
//     if (m_Program == "psiblast" || m_Program == "deltablast") {
//         m_Outfile << "\n";
//         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
//                               m_Outfile, CReference::eCompAdjustedMatrices);
//     }
//     if (m_Program == "psiblast" || m_Program == "blastp") {
//         m_Outfile << "\n";
//         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
//                               m_Outfile, CReference::eCompBasedStats,
//                               (bool)(m_Program == "psiblast"));
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:81-88
// ```c++
//     string reference("Reference");
//     if (pub == blast::CReference::eCompAdjustedMatrices) {
//         reference += " for compositional score matrix adjustment";
//     } else if (pub == blast::CReference::eCompBasedStats) {
//         reference += " for composition-based statistics";
//         if (is_psiblast) {
//             reference += " starting in round 2";
//         }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/version.cpp:62-72
// ```c
//     // eCompBasedStats
//     "Alejandro A. Sch&auml;ffer, L. Aravind, Thomas L. Madden, Sergei Shavirin, \
// John L. Spouge, Yuri I. Wolf, Eugene V. Koonin, and Stephen F. Altschul \
// (2001), \"Improving the accuracy of PSI-BLAST protein database searches \
// with composition-based statistics and other refinements\", Nucleic Acids \
// Res. 29:2994-3005.",
//     // eCompAdjustedMatrices
//     "Stephen F. Altschul, John C. Wootton, E. Michael Gertz, Richa Agarwala, \
// Aleksandr Morgulis, Alejandro A. Sch&auml;ffer, and Yi-Kuo Yu (2005) \"Protein \
// database searches using compositionally adjusted substitution matrices\", \
// FEBS J. 272:5101-5109.",
// ```
/// psiblast's outfmt 0 intro: the program, the gapped BLAST reference, the reference of the
/// compositionally adjusted matrices and the reference of composition-based statistics "starting
/// in round 2", each wrapped at 68 columns (`x_WrapOutputLine`) as NCBI writes them.
fn write_psiblast_pairwise_intro<W: Write>(writer: &mut W, version: &str) -> io::Result<()> {
    writeln!(writer, "PSIBLAST {}", version)?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A."
    )?;
    writeln!(
        writer,
        "Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J."
    )?;
    writeln!(
        writer,
        "Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of"
    )?;
    writeln!(
        writer,
        "protein database search programs\", Nucleic Acids Res. 25:3389-3402."
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Reference for compositional score matrix adjustment: Stephen F."
    )?;
    writeln!(
        writer,
        "Altschul, John C. Wootton, E. Michael Gertz, Richa Agarwala,"
    )?;
    writeln!(
        writer,
        "Aleksandr Morgulis, Alejandro A. Schaffer, and Yi-Kuo Yu (2005)"
    )?;
    writeln!(
        writer,
        "\"Protein database searches using compositionally adjusted"
    )?;
    writeln!(writer, "substitution matrices\", FEBS J. 272:5101-5109.")?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Reference for composition-based statistics starting in round 2:"
    )?;
    writeln!(
        writer,
        "Alejandro A. Schaffer, L. Aravind, Thomas L. Madden, Sergei"
    )?;
    writeln!(
        writer,
        "Shavirin, John L. Spouge, Yuri I. Wolf, Eugene V. Koonin, and"
    )?;
    writeln!(
        writer,
        "Stephen F. Altschul (2001), \"Improving the accuracy of PSI-BLAST"
    )?;
    writeln!(
        writer,
        "protein database searches with composition-based statistics and"
    )?;
    writeln!(
        writer,
        "other refinements\", Nucleic Acids Res. 29:2994-3005."
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)?;
    Ok(())
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
    let bytes = text;
    let mut pos = 0;
    // NCBI reference: c++/src/corelib/ncbistr.cpp:5100,5153-5157,5137
    // SIZE_TYPE pos=0,len=str.size(),nl_pos=0; if(nl_pos<=pos) nl_pos=str.find('\n',pos);
    // bool thisPartHasBackspace = false;
    let mut newline = 0;
    while pos < bytes.len() {
        let mut has_backspace = false;
        if newline <= pos {
            newline = bytes[pos..]
                .iter()
                .position(|&c| c == b'\n')
                .map_or(bytes.len(), |n| pos + n);
        }
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
        let mut line = Vec::with_capacity(best - pos);
        for &c in &bytes[pos..best] {
            if has_backspace && c == 8 {
                line.pop();
            } else {
                line.push(c);
            }
        }
        crate::cli::insert(writer, &line)?;
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
    writeln!(
        writer,
        "           {} sequences; {} total letters",
        format_count_with_commas(i64::try_from(database.num_sequences).unwrap_or(i64::MAX)),
        format_count_with_commas(i64::try_from(database.total_letters).unwrap_or(i64::MAX))
    )?;
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
        writeln!(writer, "    Posted date:  {}", database.date)?;
        writeln!(
            writer,
            "  Number of letters in database: {}",
            format_count_with_commas(i64::try_from(database.total_letters).unwrap_or(i64::MAX))
        )?;
        writeln!(
            writer,
            "  Number of sequences in database:  {}",
            format_count_with_commas(i64::try_from(database.num_sequences).unwrap_or(i64::MAX))
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
        format_count_with_commas(i64::try_from(database_num_sequences).unwrap_or(0)),
        format_count_with_commas(i64::try_from(database_total_letters).unwrap_or(0))
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
    writeln!(writer, "Matrix: {}", report.matrix_name)?;
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
        format_count_with_commas(i64::try_from(database_total_letters).unwrap_or(0))
    )?;
    writeln!(
        writer,
        "  Number of sequences in database:  {}",
        format_count_with_commas(i64::try_from(database_num_sequences).unwrap_or(0))
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
    )
}

/// The queries `query_range` of a BLASTP pairwise report, whose `hits` are those queries' (the
/// rows of a caller that makes them a window of queries at a time, as BLASTN's and TBLASTX's
/// do, `x_pairwise_par::query_windows`): the prolog with the report's first query, the epilog
/// with its last.
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
                Some(subject) => DatabaseHeading::parsed(subject, true, group.own_hsps(s_idx)),
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
    mut label_of: impl FnMut(u32) -> Vec<u8>,
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
        let label = label_of(*s_idx);
        // NCBI reference: c++/src/objtools/align_format/showdefline.cpp:915-918,930
        // ```c++
        //         if(line_component.size()+line_length > m_LineLen){
        //             actual_line_component = line_component.substr(0, m_LineLen -
        //                                                           line_length - 3);
        //             actual_line_component += kEllipsis;
        //     ...
        //         CAlignFormatUtil::AddSpace(out, m_LineLen - line_length);
        // ```
        // String widths are byte counts, including non-ASCII FASTA titles.
        if label.len() > 68 {
            writer.write_all(&label[..65])?;
            writer.write_all(b"...")?;
        } else {
            writer.write_all(&label)?;
            write_spaces(writer, 68 - label.len())?;
        }
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

/// The description-table label of a database sequence (before the cut to the line): the title
/// that `GenerateDefline` makes of the stored title `title` and the source `source` of the
/// Bioseq, which an empty title is made from (`Unknown` where NCBI's `HtmlDecode` throws), or for
/// a sequence with parsed Seq-ids its id and title (`ParsedSubject::description`).
pub fn database_description_label(
    title: &[u8],
    source: &super::defline::SourceTaxname,
    parsed: Option<&crate::report::db_subject::ParsedSubject>,
    protein: bool,
) -> Vec<u8> {
    if let Some(subject) = parsed {
        let (mut label, title) = subject.description(protein);
        label.extend_from_slice(&title.text);
        return label;
    }
    let facts = super::defline::BioseqFacts::of_source(source.name().ok().flatten());
    match super::defline::generate_bioseq_defline(title, protein, true, &facts) {
        Ok(defline) => defline.text,
        Err(super::defline::UnknownEncoding) => b"Unknown".to_vec(),
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
    read: super::query_warnings::SubjectReadHook<'_, W>,
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
            let title = subject_titles.get(&s_idx).map_or(&[][..], Vec::as_slice);
            database_description_label(
                title,
                stored_source(sources, s_idx),
                parsed_subject(parsed, s_idx),
                protein,
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
// NCBI reference: c++/src/algo/blast/api/version.cpp:58-61
// ```c++
//     // eMegaBlast
//     "Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), \
// \"A greedy algorithm for aligning DNA sequences\", \
// J Comput Biol 2000; 7(1-2):203-14.",
// ```
// The reference is printed after "Reference: " and wrapped at 68 columns
// (blastfmtutil.cpp:116-121); the layout is that of `write_translated_pairwise_intro`.
fn write_megablast_pairwise_intro(writer: &mut impl Write, version: &str) -> io::Result<()> {
    writeln!(writer, "BLASTN {version}")?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Reference: Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb"
    )?;
    writeln!(
        writer,
        "Miller (2000), \"A greedy algorithm for aligning DNA sequences\", J"
    )?;
    writeln!(writer, "Comput Biol 2000; 7(1-2):203-14.")?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)
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
}

impl super::x_pairwise_par::QuerySections for BlastnSections<'_> {
    fn num_alignments(&self) -> usize {
        self.report.num_alignments
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
        // id is not shown), as BLASTP's (`write_database_subject_alignments`), or the lines
        // of the definition lines of one with parsed Seq-ids.
        let parsed = parsed_subject(&self.report.parsed_subjects, s_idx);
        if self.report.database.is_some() || parsed.is_some() {
            let heading = match parsed {
                Some(subject) => {
                    DatabaseHeading::parsed(subject, self.report.gapped, group.own_hsps(s_idx))
                }
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
// Both rows are translated (widths 3): a 60-residue chunk spans 180 nucleotides, each row
// runs in the direction of its frame, and a minus-frame row is not reversed
// (showalign.cpp:1852-1858 applies only without widths). The coordinate width is that of
// the largest 0-based nucleotide coordinate (`coordinate_width`). The query row may carry
// lowercase masks, which NCBI applies only when it prints the row (showalign.cpp:2495-2521),
// so the middle line compares the uppercase residues.
fn write_tblastx_alignment<W: Write>(
    writer: &mut W,
    hit: &PairwiseHit,
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
    let direction = |frame: Option<i8>| -> isize {
        if frame.unwrap_or(1) < 0 {
            -1
        } else {
            1
        }
    };
    // Both rows are translated: three nucleotides per residue; no row has a gap.
    let line_length = config.line_length;
    let mut query_row = RowCursor::new(
        hit.hit.q_start as isize,
        direction(hit.query_frame),
        3,
        line_length,
    );
    let mut subject_row = RowCursor::new(
        hit.hit.s_start as isize,
        direction(hit.subject_frame),
        3,
        line_length,
    );
    for (qpart, spart) in q.chunks(line_length).zip(s.chunks(line_length)) {
        write_sequence_row(writer, "Query", &mut query_row, qpart, 0, width, |writer| {
            writer.write_all(qpart)
        })?;
        write_spaces(writer, 7 + width + 2)?;
        for (&qc, &sc) in qpart.iter().zip(spart) {
            let query_residue = qc.to_ascii_uppercase();
            let mid = if query_residue == sc {
                query_residue
            } else if is_positive_match(query_residue as char, sc as char, config.protein_matrix) {
                b'+'
            } else {
                b' '
            };
            writer.write_all(&[mid])?;
        }
        writeln!(writer)?;
        write_sequence_row(
            writer,
            "Sbjct",
            &mut subject_row,
            spart,
            0,
            width,
            |writer| writer.write_all(spart),
        )?;
        writeln!(writer)?;
    }
    Ok(())
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
    writeln!(writer, "Matrix: {}", report.matrix_name)?;
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
}

impl super::x_pairwise_par::QuerySections for TblastxSections<'_> {
    fn num_alignments(&self) -> usize {
        self.report.num_alignments
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
        let write_hsp = |writer: &mut W, hit: &PairwiseHit| {
            write_tblastx_hsp_info(writer, hit)?;
            write_tblastx_alignment(writer, hit, config)?;
            writeln!(writer)
        };
        // A database sequence's heading is its stored title after `>` (as BLASTP's,
        // BLASTN's and TBLASTN's `-db` reports), or the lines of the definition lines of one
        // with parsed Seq-ids.
        let parsed = parsed_subject(&self.report.parsed_subjects, s_idx);
        if self.database.is_some() || parsed.is_some() {
            let heading = match parsed {
                Some(subject) => DatabaseHeading::parsed(subject, false, group.own_hsps(s_idx)),
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
            &group.subject_titles[&s_idx],
            false,
            &group.subject_hits[&s_idx],
            &group.subject_hit_indices[&s_idx],
            part,
            write_hsp,
        )
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
            }]),
        )
        .unwrap();
        assert_eq!(
            String::from_utf8(written).unwrap(),
            "Database: Swiss-Prot\n           575,748 sequences; 209,017,843 total letters\n\n"
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
        write_database_subject_alignments(
            &mut written,
            None,
            title,
            true,
            &[&hit],
            &[0],
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
            text(&|out| write_tblastx_alignment(out, &tblastx, &config)),
            "Query  1  M  3\n          M\nSbjct  1  M  3\n\n\
             Query  4  R  6\n          R\nSbjct  4  R  6\n\n"
        );
        // On a minus frame the first residue of a translated row at 1 to 3 runs from 3 to 1.
        let mut minus = hit(1, 2, Some(-1), None);
        minus.hit.q_start = 6;
        assert_eq!(
            text(&|out| write_blastx_alignment(out, &minus, &config)),
            "Query  6  M  4\n          M\nSbjct     M  \n\n\
             Query  3  R  1\n          R\nSbjct  2  R  2\n\n"
        );
    }

    // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:204,230
    // ```c++
    //                                       unsigned int  scrn_width,
    //     ...
    //     scrn_width *= width;
    // ```
    // The screen of a translated row wraps modulo 2^32 from `-line_length` 1431655766 on: 2
    // nucleotides there, 5 at 1431655767, 2^32 - 1 at 1431655765; a row of letters keeps
    // `-line_length`.
    #[test]
    fn a_translated_row_s_screen_wraps_at_2_to_the_32() {
        let screen = |letters, line_length| RowCursor::new(1, 1, letters, line_length).screen;
        assert_eq!(screen(3, 1431655765), u32::MAX);
        assert_eq!(screen(3, 1431655766), 2);
        assert_eq!(screen(3, 1431655767), 5);
        assert_eq!(screen(3, 1431655800), 104);
        assert_eq!(screen(3, 2147483647), 2147483645);
        assert_eq!(screen(3, 60), 180);
        assert_eq!(screen(1, 2147483647), 2147483647);
        assert_eq!(screen(1, 1431655766), 1431655766);
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
        let blastx = hit((2707, 2693), (4, 8), Some(-1), None);
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
        let tblastx = hit((2, 16), (598, 584), Some(1), Some(-1));
        for (line_length, ends) in [
            (1431655765, (16, 584)),
            (1431655766, (3, 597)),
            (1431655767, (6, 594)),
            (2147483647, (16, 584)),
        ] {
            assert_eq!(
                text(&|out| write_tblastx_alignment(out, &tblastx, &config(line_length))),
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

    /// NCBI's `GetWholeAlnSeqString` (alnvec.cpp:199-397, quoted at `RowCursor`) for one row of
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
    // `RowCursor::first_screen`): rows of residues and gaps, on both strands, with screens
    // narrower than the row (wrapped ones: 3 * L - 2^32 nucleotides, 2, 5, 8, ...) and wider,
    // and residue runs split into segments at random (the other row's gaps split them): the
    // end that `write_sequence_row` writes is the first screen's right of NCBI's computation.
    #[test]
    fn the_first_chunk_ends_at_ncbi_s_first_screen() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut random = |bound: u64| {
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            state.wrapping_mul(0x2545_f491_4f6c_dd1d) % bound
        };
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
            let mut cursor = RowCursor::new(
                start as isize,
                if plus { 1 } else { -1 },
                letters as isize,
                line_length,
            );
            let mut out = Vec::new();
            write_sequence_row(&mut out, "Query", &mut cursor, &row, 0, 5, |out| {
                out.write_all(&row)
            })
            .unwrap();
            let line = String::from_utf8(out).unwrap();
            let words: Vec<&str> = line.split_whitespace().collect();
            // A first chunk of one letter at 1 shows no coordinates (its end is 1 too).
            let end: i64 = match words[..] {
                [_, _] => 1,
                _ => words.last().unwrap().parse().unwrap(),
            };
            // The row's segments, runs of residues split at random.
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
            let scrn_width = (line_length as u32).wrapping_mul(letters);
            assert_eq!(
                end,
                ncbi_first_screen_right(&segments, plus, letters, scrn_width) + 1,
                "{line:?} plus={plus} letters={letters} L={line_length} segments={segments:?}"
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
    mut probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut super::query_warnings::QueryWarnings<'_>>,
    epilog: bool,
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
    for (q_idx, query) in queries.iter().enumerate() {
        // The query's warnings come before its preamble (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut writer)?;
        }
        // NCBI blast_format.cpp:1491: m_Outfile << "\n\n";
        writeln!(&mut writer)?;
        writeln!(&mut writer)?;
        write_blastp_query_header(&mut writer, &query_titles[q_idx], query.query_length)?;
        let query_hits = &hits_by_query[q_idx];
        if query_hits.is_empty() {
            write_no_hits_found(&mut writer)?;
            // NCBI blast_results.cpp:82-97 leaves both Karlin blocks null
            // when there is no valid context. blast_format.cpp:445-477
            // writes only the blank lines and zero search space in that case.
            // NCBI c++/src/algo/blast/api/local_blast.cpp:177-180,204-208:
            // if (m_PrelimSearch->CheckInternalData() != 0)
            //     new CBlastAncillaryData(tmp_pair, tmp_pair, tmp_pair, 0);
            // The all-invalid batch uses -1 sentinel blocks. An invalid query
            // within a searched batch has null blocks (blast_results.cpp:82-103).
            if query_batch_skipped[q_idx] {
                write_tblastn_unsearched_query_footer_spacing(&mut writer, false)?;
            } else if !query_validity[q_idx] {
                writeln!(writer)?;
                writeln!(writer)?;
                writeln!(writer)?;
                writeln!(writer, "Effective search space used: 0")?;
            } else {
                write_blastp_query_footer_spacing(
                    &mut writer,
                    query.ungapped_karlin,
                    report.gapped_karlin,
                    report.gumbel,
                    query.effective_search_space,
                    false,
                )?;
            }
            continue;
        }
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
        // The subjects' hits; a copy of one record that follows another is shown in its hit
        // (`QueryGroup::new`, `db_subject::joins_as_copy`).
        let mut group = super::x_pairwise_par::QueryGroup::new(
            query_hits,
            subject_titles,
            |previous, s_idx| {
                // `-subject` records with one id are refused before any report.
                database.is_some()
                    && crate::report::db_subject::rank_joins_as_copy(
                        parsed_subjects,
                        previous as usize,
                        s_idx as usize,
                    )
            },
        )
        .expect("a query with hits has a group");
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1522-1530
        // ```c++
        //     if(m_HitsSortOption > 0 || m_HspsSortOption > 0) {
        //         aln_set = CBlastFormatUtil::SortSeqalignForSortableFormat(
        // ```
        // `-sorthits` and `-sorthsps` before the table and the alignments (`QueryGroup::sort`).
        group.sort(config.sort, config.tblastx());
        let super::x_pairwise_par::QueryGroup {
            subject_order,
            subject_hits,
            subject_hit_indices,
            subject_titles,
        } = &group;
        let described = &subject_order[..subject_order.len().min(report.num_descriptions)];
        // NCBI loads each subject that the query's batch did not prefetch at its description
        // line, or one the table does not show at its alignment's heading
        // (`QueryWarnings::subject_read`).
        let mut read =
            |s_idx: u32, writer: &mut io::BufWriter<&mut W>| match warnings.as_deref_mut() {
                Some(warnings) => warnings.subject_read(q_idx, s_idx, writer),
                None => Ok(()),
            };
        write_blastn_description_table(
            &mut writer,
            described,
            subject_order.len() <= report.num_descriptions,
            subject_hits,
            subject_titles,
            parsed_subjects,
            subject_sources,
            false,
            false,
            config.sort_columns(query.query_length),
            &mut read,
        )?;
        writeln!(writer)?;
        // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:1014-1040
        // (PruneSeqalign keeps the alignments of the first `m_NumAlignments` subjects).
        let aligned = &subject_order[..subject_order.len().min(report.num_alignments)];
        for &s_idx in aligned {
            read(s_idx, &mut writer)?;
            // NCBI c++/src/objtools/align_format/showalign.cpp:3650-3668:
            // display leaves an extra blank line after each translated HSP.
            let write_hsp = |writer: &mut io::BufWriter<&mut W>, hit: &PairwiseHit| {
                write_tblastn_hsp_info(writer, hit)?;
                write_tblastn_alignment(writer, hit, config)?;
                writeln!(writer)
            };
            // A database sequence's heading is its stored title after `>` (as BLASTP's
            // and BLASTN's `-db` reports), or the lines of the definition lines of one with
            // parsed Seq-ids.
            let parsed = parsed_subject(parsed_subjects, s_idx);
            if database.is_some() || parsed.is_some() {
                let heading = match parsed {
                    Some(subject) => DatabaseHeading::parsed(subject, true, group.own_hsps(s_idx)),
                    None => DatabaseHeading::Title(
                        &subject_titles[&s_idx],
                        stored_source(subject_sources, s_idx),
                    ),
                };
                write_database_subject_alignments(
                    &mut writer,
                    probe.as_deref_mut(),
                    heading,
                    false,
                    &subject_hits[&s_idx],
                    &subject_hit_indices[&s_idx],
                    write_hsp,
                )?;
                continue;
            }
            write_subject_alignments(
                &mut writer,
                probe.as_deref_mut(),
                s_idx,
                &subject_titles[&s_idx],
                false,
                &subject_hits[&s_idx],
                &subject_hit_indices[&s_idx],
                write_hsp,
            )?;
        }
        write_blastp_query_footer_spacing(
            &mut writer,
            query.ungapped_karlin,
            report.gapped_karlin,
            report.gumbel,
            query.effective_search_space,
            false,
        )?;
    }
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
    if let Some(frame) = hit.query_frame.or(hit.subject_frame) {
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
// NCBI c++/src/algo/blast/core/blast_hits.c:1087-1105:
// translated subject offsets advance by CODON_LENGTH per aligned residue.
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
    let direction: isize = if hit.subject_frame.unwrap_or(1) < 0 {
        -1
    } else {
        1
    };
    // The query's row of letters, and the subject's translated row: three nucleotides per
    // residue, in the direction of its frame.
    let line_length = config.line_length;
    let mut query_row = RowCursor::new(hit.hit.q_start as isize, 1, 1, line_length);
    let mut subject_row = RowCursor::new(hit.hit.s_start as isize, direction, 3, line_length);
    for (qpart, spart) in q.chunks(line_length).zip(s.chunks(line_length)) {
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1605-1627
        // ```c++
        //     //not to display start and stop number for empty row
        //     if ((j > 0 && end == prev_stop)
        //         || (j == 0 && start == 1 && end == 1)) {
        //         startLen = 0;
        //     } else {
        //         out << start;
        // ...
        //      //not to display stop number for empty row in the middle
        //     if (!(j > 0 && end == prev_stop)
        // ```
        // A row of gaps only shows no coordinates, nor a first row of one residue at 1
        // (`write_sequence_row`).
        write_sequence_row(writer, "Query", &mut query_row, qpart, 0, width, |writer| {
            writer.write_all(qpart)
        })?;
        write_spaces(writer, 7 + width + 2)?;
        for (&qc, &sc) in qpart.iter().zip(spart) {
            // NCBI c++/src/objtools/align_format/showalign.cpp:2122-2149:
            // the displayed query preserves lowercase masking, while the
            // protein matrix comparison uses its uppercase residue.
            let query_residue = qc.to_ascii_uppercase();
            let mid = if query_residue == sc {
                query_residue
            } else if is_positive_match(query_residue as char, sc as char, config.protein_matrix) {
                b'+'
            } else {
                b' '
            };
            writer.write_all(&[mid])?;
        }
        writeln!(writer)?;
        write_sequence_row(
            writer,
            "Sbjct",
            &mut subject_row,
            spart,
            0,
            width,
            |writer| writer.write_all(spart),
        )?;
        writeln!(writer)?;
    }
    Ok(())
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
    writeln!(writer, "{} {}", program, version)?;
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
    writeln!(
        writer,
        "Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A."
    )?;
    writeln!(
        writer,
        "Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J."
    )?;
    writeln!(
        writer,
        "Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of"
    )?;
    writeln!(
        writer,
        "protein database search programs\", Nucleic Acids Res. 25:3389-3402."
    )?;
    writeln!(writer)?;
    writeln!(writer)?;
    writeln!(writer)?;
    Ok(())
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
/// Writes the BLASTX pairwise report (outfmt 0) of `queries`, whose hits are `hits` in the
/// order of the reports (`q_idx` is the index in `queries`). `probe` marks each subject's
/// heading and each HSP, `hits[i]` being HSP `first_hsp + i`, without changing the written bytes.
#[allow(clippy::too_many_arguments)]
pub fn write_blastx_pairwise_report<W: Write>(
    hits: &[PairwiseHit],
    writer: &mut W,
    stderr: &mut dyn Write,
    warnings: &[Vec<u8>],
    reads: &[super::query_warnings::QueryReads],
    config: &PairwiseConfig,
    queries: &[BlastpPairwiseQuery],
    query_validity: &[bool],
    skipped: &[bool],
    subject_ids: &[Arc<str>],
    report: &BlastpPairwiseReport,
    options: &BlastxPairwiseOptions,
    database: Option<(
        &PairwiseDatabase,
        &[Arc<[u8]>],
        &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
        &[super::defline::SourceTaxname],
    )>,
    records_parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    mut probe: Option<&mut FormatProbe<'_>>,
    first_hsp: HspIndex,
) -> io::Result<()> {
    let mut writer = io::BufWriter::new(writer);
    // `hits[i]` is HSP `first_hsp + i` of the run's hit list (LOSAT Web's `run_local` writes a run
    // batch by batch).
    let mut by_query: Vec<Vec<(HspIndex, &PairwiseHit)>> = vec![Vec::new(); queries.len()];
    for (index, hit) in hits.iter().enumerate() {
        by_query[hit.hit.q_idx as usize].push((first_hsp + index, hit));
    }
    // The texts of the subjects that each query's report reads one by one, written where it
    // reads each (`ReadCursor`): description lines and alignment headings.
    let mut reads = super::query_warnings::ReadCursor::new(reads);
    for (qidx, query) in queries.iter().enumerate() {
        // A text whose read the reports before did not place comes after them.
        if let Some(reads) = reads.as_mut() {
            reads.before_query(qidx, &mut writer, &mut *stderr)?;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452
        // ```c++
        //     if (results.HasWarnings()) {
        //         ERR_POST(Warning << results.GetWarningStrings());
        //     }
        // ```
        stderr.write_all(&warnings[qidx])?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1491-1492
        // ```c++
        //     m_Outfile << "\n\n";
        //     CBlastFormatUtil::AcknowledgeBlastQuery(*bioseq, kFormatLineLength,
        // ```
        writeln!(writer)?;
        writeln!(writer)?;
        write_blastp_query_header(&mut writer, query.query_name.as_bytes(), query.query_length)?;
        let query_hits = &by_query[qidx];
        if query_hits.is_empty() {
            write_no_hits_found(&mut writer)?;
        } else {
            // The subjects' hits; a copy of one record of the database that follows another
            // is shown in its hit (`QueryGroup::new`, `db_subject::joins_as_copy`).
            let mut group =
                super::x_pairwise_par::QueryGroup::new(query_hits, &[], |previous, s_idx| {
                    database.is_some_and(|(_, _, parsed, _)| {
                        crate::report::db_subject::rank_joins_as_copy(
                            parsed,
                            previous as usize,
                            s_idx as usize,
                        )
                    })
                })
                .expect("a query with hits has a group");
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1522-1530
            // ```c++
            //     if(m_HitsSortOption > 0 || m_HspsSortOption > 0) {
            //         aln_set = CBlastFormatUtil::SortSeqalignForSortableFormat(
            // ```
            // `-sorthits` and `-sorthsps` before the table and the alignments
            // (`QueryGroup::sort`).
            group.sort(config.sort, config.tblastx());
            let (grouped, indices, order) = (
                &group.subject_hits,
                &group.subject_hit_indices,
                &group.subject_order,
            );
            // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:606-611
            // ```c++
            //     CShowBlastDefline showdef(*aln_set, *m_Scope,
            //                               defline_length == -1 ? kFormatLineLength:defline_length,
            //                               m_NumSummary + additional);
            // ```
            // The table follows `x_InitDeflineTable` (`write_description_table`): the
            // subject's highest bit score and that HSP's E-value, for the first
            // `m_NumDescriptions` subjects; the last one is NCBI's "last hit" when every
            // subject fits.
            let descriptions = &order[..order.len().min(options.num_descriptions)];
            write_subject_summary_table_with_sum_n(
                &mut writer,
                descriptions,
                order.len() <= options.num_descriptions,
                grouped,
                subject_ids,
                !options.gapped && options.sum_stats,
                config.sort_columns(query.query_length),
                &mut |s_idx, writer: &mut io::BufWriter<&mut W>| match reads.as_mut() {
                    Some(reads) => reads.subject_read(qidx, s_idx, writer, &mut *stderr),
                    None => Ok(()),
                },
            )?;
            writeln!(writer)?;
            for oid in order.iter().take(options.num_alignments) {
                if let Some(reads) = reads.as_mut() {
                    reads.subject_read(qidx, *oid, &mut writer, &mut *stderr)?;
                }
                let shits = &grouped[oid];
                let first = shits[0];
                let hsp_indices = indices.get(oid).map_or(&[][..], Vec::as_slice);
                if let Some(probe) = probe.as_deref_mut() {
                    writer.flush()?;
                    probe.subject_begin(hsp_indices[0]);
                }
                // A database sequence's heading is its stored title after `>` (as
                // BLASTP's, BLASTN's and TBLASTN's `-db` reports), or NCBI's line in place of
                // each alignment where it cannot make the heading's title
                // (`write_database_subject_alignments_part`).
                let heading = match database {
                    Some((_, titles, parsed, sources)) => {
                        Some(match parsed_subject(parsed, *oid) {
                            Some(subject) => DatabaseHeading::parsed(
                                subject,
                                options.gapped,
                                group.own_hsps(*oid),
                            ),
                            None => DatabaseHeading::Title(
                                &titles[*oid as usize],
                                stored_source(sources, *oid),
                            ),
                        })
                    }
                    None => parsed_subject(records_parsed, *oid).map(|subject| {
                        DatabaseHeading::parsed(subject, options.gapped, shits.len())
                    }),
                };
                let skipped = heading.and_then(|heading| heading.skipped_id(true));
                match heading {
                    Some(heading) if skipped.is_none() => write_database_subject_heading(
                        &mut writer,
                        heading,
                        true,
                        first.subject_length,
                    )?,
                    Some(_) => {}
                    None => write_subject_header(
                        &mut writer,
                        subject_ids[*oid as usize].as_bytes(),
                        first.subject_title.as_deref().map(str::as_bytes),
                        first.subject_length,
                    )?,
                }
                if let Some(probe) = probe.as_deref_mut() {
                    writer.flush()?;
                    probe.subject_end(hsp_indices[0]);
                }
                for (position, hit) in shits.iter().enumerate() {
                    if let Some(probe) = probe.as_deref_mut() {
                        writer.flush()?;
                        probe.begin(hsp_indices[position]);
                    }
                    match &skipped {
                        None => {
                            write_tblastn_hsp_info(&mut writer, hit)?;
                            write_blastx_alignment(&mut writer, hit, config)?;
                            writeln!(writer)?;
                        }
                        Some(id) => write_skipped_alignment(&mut writer, id)?,
                    }
                    if let Some(probe) = probe.as_deref_mut() {
                        writer.flush()?;
                        probe.end(hsp_indices[position]);
                    }
                }
            }
        }
        write_blastx_query_footer(
            &mut writer,
            query,
            report,
            query_validity[qidx],
            skipped[qidx],
            options.gapped,
        )?;
    }
    if let Some(reads) = reads.as_mut() {
        reads.before_query(queries.len(), &mut writer, &mut *stderr)?;
    }
    writer.flush()
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
            format_count_with_commas(report.database_total_letters as i64)
        )?;
        writeln!(
            writer,
            "  Number of sequences in database:  {}",
            format_count_with_commas(report.database_num_sequences as i64)
        )?;
        writeln!(writer)?;
        writeln!(writer)?;
        writeln!(writer)?;
    }
    writeln!(writer, "Matrix: {}", report.matrix_name)?;
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
fn write_blastx_query_footer(
    writer: &mut impl Write,
    query: &BlastpPairwiseQuery,
    report: &BlastpPairwiseReport,
    valid: bool,
    skipped: bool,
    gapped: bool,
) -> io::Result<()> {
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
    let direction: isize = if hit.query_frame.unwrap_or(1) < 0 {
        -1
    } else {
        1
    };
    // The query's translated row (three nucleotides per residue, in the direction of its frame)
    // and the subject's row of letters.
    let line_length = config.line_length;
    let mut query_row = RowCursor::new(hit.hit.q_start as isize, direction, 3, line_length);
    let mut subject_row = RowCursor::new(hit.hit.s_start as isize, 1, 1, line_length);
    for (chunk, (qpart, spart)) in q.chunks(line_length).zip(s.chunks(line_length)).enumerate() {
        let first = chunk * line_length;
        write_sequence_row(writer, "Query", &mut query_row, qpart, 0, width, |writer| {
            writer.write_all(qpart)
        })?;
        write_spaces(writer, 7 + width + 2)?;
        // Both rows may carry lower case (the query's and, with `-lcase_masking`, a `-subject`'s,
        // `blastx/x_soft_subject.rs`), which NCBI applies only when it prints a row
        // (showalign.cpp:2495-2521), so the residues are compared in upper case.
        for (&qc, &sc) in qpart.iter().zip(spart) {
            let (qu, su) = (qc.to_ascii_uppercase(), sc.to_ascii_uppercase());
            let middle = if qu == su {
                qu
            } else if is_positive_match(qu as char, su as char, config.protein_matrix) {
                b'+'
            } else {
                b' '
            };
            writer.write_all(&[middle])?;
        }
        writeln!(writer)?;
        // The stored gap characters of the subject are residues (`residue_count`).
        write_sequence_row(
            writer,
            "Sbjct",
            &mut subject_row,
            spart,
            stored_gaps_between(&hit.subject_stored_gaps, first, first + spart.len()),
            width,
            |writer| writer.write_all(spart),
        )?;
        writeln!(writer)?;
    }
    Ok(())
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
        let before: Vec<Vec<u8>> = (0..queries.len())
            .map(|q_idx| match q_idx {
                0 | 2 => format!(
                    "Warning: [tblastx] Query_{} q{q_idx}: a warning\n",
                    q_idx + 1
                )
                .into_bytes(),
                _ => Vec::new(),
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
                        Some(&mut crate::report::query_warnings::QueryWarnings::new(
                            &before,
                            &mut diagnostics,
                        )),
                        x_pool,
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
        write_tblastx_alignment(&mut out, &tblastx_hit(1), &PairwiseConfig::default()).unwrap();
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
        let ids = [Arc::<str>::from("YOL075C")];
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
        )
        .unwrap();
        let row = out.split(|&byte| byte == b'\n').nth(3).unwrap().to_vec();
        assert!(row.starts_with(&[&b"Unknown"[..], &[b' '; 61], b"  120"].concat()));
        let mut out = Vec::new();
        write_subject_alignments(
            &mut out,
            None,
            0,
            &titles[&0],
            false,
            &hits[&0],
            &[0, 1],
            |_, _| panic!("no HSP is written"),
        )
        .unwrap();
        assert_eq!(
            out,
            b"Sequence with id Subject_1 no longer exists in database...alignment skipped\n\
              Sequence with id Subject_1 no longer exists in database...alignment skipped\n"
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
            query_name: String::new(),
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
            matrix_name: "BLOSUM62".to_string(),
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
