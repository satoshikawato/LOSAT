use std::cmp::Ordering;
use std::io::{self, Write};
use std::sync::Arc;

use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::cli::NativeError;
use crate::common::{GapEditOp, Hit};
use crate::report::{write_hit_fields, OutputConfig};

use super::tracing as blastn_trace;

// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:770-782
// ```c
// if (m_FormatType == CFormattingArgs::eTabular ||
//     m_FormatType == CFormattingArgs::eTabularWithComments) {
//     CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec, kDelim);
// }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1410-1414
// ```c
// void
// CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
//                         CConstRef<blast::CBlastQueryVector> queries,
//                         unsigned int itr_num
// ```
// The pairwise report (outfmt 0) is the non-tabular branch of the same formatter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlastnOutputFormat {
    Pairwise,
    Tabular,
    TabularWithComments,
}

/// The output format of a `-outfmt` value, as NCBI parses it: white space around the
/// value is removed, the format number ends at the first space, and the rest (a custom
/// specification) counts only for the tabular formats, where LOSAT does not support it.
///
/// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2801-2809
/// ```c
///     if (args[kArgOutputFormat]) {
///         string fmt_choice =
///             NStr::TruncateSpaces(args[kArgOutputFormat].AsString());
///         string::size_type pos;
///         if ( (pos = fmt_choice.find_first_of(' ')) != string::npos) {
///             custom_fmt_spec.assign(fmt_choice, pos+1,
///                                    fmt_choice.size()-(pos+1));
///             fmt_choice.erase(pos);
///         }
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2828-2851
/// ```c
///         int val = 0;
///         try { val = NStr::StringToInt(fmt_choice); }
///         catch (const CStringException&) {   // probably a conversion error
///             CNcbiOstrstream os;
///             os << "'" << fmt_choice << "' is not a valid output format";
/// ...
///         fmt_type = static_cast<EOutputFormat>(val);
///         if ( !(fmt_type == eTabular ||
///                fmt_type == eTabularWithComments ||
///                fmt_type == eCommaSeparatedValues ||
///                fmt_type == eCommaSeparatedValuesWithHeader ||
///                fmt_type == eSAM) ) {
///                custom_fmt_spec.clear();
///         }
/// ```
/// NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.hpp:260-263
/// ```c
///     catch (const std::exception& e) {                                       \
///         LOG_POST(Error << "Error: " << e.what());                           \
///         exit_code = BLAST_UNKNOWN_ERROR;                                    \
///     }                                                                       \
/// ```
/// `NStr::TruncateSpaces` removes the characters of C's `isspace`, and `NStr::StringToInt`
/// reads what `i32::from_str` reads. NCBI's errors are its messages and exit statuses; the
/// formats and specifications that NCBI supports and LOSAT does not are rejected with
/// LOSAT's message.
pub fn parse_blastn_output_format(spec: &str) -> anyhow::Result<BlastnOutputFormat> {
    let is_space = |c: char| matches!(c, ' ' | '\t' | '\n' | '\x0b' | '\x0c' | '\r');
    let choice = spec.trim_matches(is_space);
    let (choice, custom) = choice.split_once(' ').unwrap_or((choice, ""));
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2810-2825
    // ```c
    //         if(!custom_fmt_spec.empty()) {
    //             if(NStr::StartsWith(custom_fmt_spec, "delim")) {
    //                 vector <string> tokens;
    //                 NStr::Split(custom_fmt_spec," ",tokens);
    //                 if(tokens.size() > 0) {
    //                     string tag;
    //                     bool isValid = NStr::SplitInTwo(tokens[0],"=",tag,custom_delim);
    //                     if(!isValid) {
    //                         string msg("Delimiter format is invalid. Valid format is delim=<delimiter value>");
    //                         NCBI_THROW(CInputException, eInvalidInput, msg);
    //                     }
    //                     else {
    //                         custom_fmt_spec = NStr::Replace(custom_fmt_spec,tokens[0],"");
    //                         custom_fmt_spec = NStr::TruncateSpaces(custom_fmt_spec);
    //                     }
    // ```
    // The delimiter is checked before the format number; an empty one is the default.
    let mut custom = custom.to_string();
    let mut custom_delimiter = false;
    if custom.starts_with("delim") {
        let token = custom.split(' ').next().unwrap_or_default().to_string();
        let Some((_, value)) = token.split_once('=') else {
            return Err(NativeError {
                exit: 1,
                message: "BLAST query/options error: Delimiter format is invalid. Valid format is delim=<delimiter value>\nPlease refer to the BLAST+ user manual.\n".to_string(),
            }
            .into());
        };
        custom_delimiter = !value.is_empty();
        custom = custom
            .replace(&token, "")
            .trim_matches(is_space)
            .to_string();
    }
    let Ok(format) = choice.parse::<i32>() else {
        return Err(NativeError {
            exit: 1,
            message: format!(
                "BLAST query/options error: '{choice}' is not a valid output format\nPlease refer to the BLAST+ user manual.\n"
            ),
        }
        .into());
    };
    if !(0..NCBI_OUTPUT_FORMAT_END).contains(&format) {
        return Err(NativeError {
            exit: 255,
            message: "Error: Formatting choice is out of range\n".to_string(),
        }
        .into());
    }
    let format = match format {
        0 => BlastnOutputFormat::Pairwise,
        6 => BlastnOutputFormat::Tabular,
        7 => BlastnOutputFormat::TabularWithComments,
        _ => anyhow::bail!("output format {format} is not supported by LOSAT's BLASTN"),
    };
    // The specification and the delimiter count only for the tabular formats.
    if format != BlastnOutputFormat::Pairwise && (!custom.is_empty() || custom_delimiter) {
        anyhow::bail!(
            "the custom output format specification {spec:?} (fields or a delimiter) is not supported by LOSAT's BLASTN"
        );
    }
    Ok(format)
}

/// NCBI reference: ncbi-blast/c++/include/algo/blast/blastinput/blast_args.hpp:1067-1072
/// ```c
///         eCommaSeparatedValuesWithHeader,
///
///         /// unaligned reads in magicblast
///         eFasta,
///         /// Sentinel value for error checking
///         eEndValue
/// ```
/// `eEndValue` is 22 (`ePairwise` is 0, `eFasta` 21).
const NCBI_OUTPUT_FORMAT_END: i32 = 22;

pub(crate) const NCBI_BLASTN_VERSION: &str = "2.17.0+";

// NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1264-1284
// ```c
// x_PrintQueryAndDbNames(program_version, bioseq, dbname, rid, iteration,
//                        subj_bioseq);
// if (align_set) {
//     int num_hits = align_set->Get().size();
//     if (num_hits != 0) PrintFieldNames(is_csv);
//     m_Ostream << "# " << num_hits << " hits found" << "\n";
// }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:762
// ```c
//     CConstRef<CSeq_align_set> aln_set = results.GetSeqAlign();
// ```
// A query of a batch that NCBI did not search has no alignment set (`num_hits` is
// `None`), so its header has no count (local_blast.cpp:177-207).
fn write_blastn_outfmt7_header<W: Write>(
    writer: &mut W,
    query_title: &str,
    subject_title: &str,
    num_hits: Option<usize>,
) -> io::Result<()> {
    writeln!(writer, "# BLASTN {NCBI_BLASTN_VERSION}")?;
    // NCBI reference: c++/src/objtools/align_format/tabular.cpp:1305-1308
    // ```c
    //     CAlignFormatUtil::AcknowledgeBlastQuery(bioseq, kLineLength, m_Ostream,
    //                                             m_ParseLocalIds, kHtmlFormat,
    //                                             kTabularFormat, rid);
    // ```
    // The title's bytes (`write_outfmt7_query_line`).
    crate::report::outfmt6::write_outfmt7_query_line(writer, query_title.as_bytes())?;
    writeln!(writer, "# Database: {subject_title}")?;
    let Some(num_hits) = num_hits else {
        return Ok(());
    };
    if num_hits > 0 {
        writeln!(
            writer,
            "# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score"
        )?;
    }
    writeln!(writer, "# {num_hits} hits found")
}

fn blastn_hsp_count(hit_list: Option<&BlastnHitList>) -> usize {
    hit_list
        .map(|list| {
            list.hsplist_array
                .iter()
                .map(|hsp_list| hsp_list.hsps.len())
                .sum()
        })
        .unwrap_or(0)
}

#[derive(Debug, Clone)]
// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-148
// ```c
// typedef struct BlastHSP {
//    Int4 score;           /**< This HSP's raw score */
//    Int4 num_ident;       /**< Number of identical base pairs in this HSP */
//    double bit_score;     /**< Bit score, calculated from score */
//    double evalue;        /**< This HSP's e-value */
//    BlastSeg query;       /**< Query sequence info. */
//    BlastSeg subject;     /**< Subject sequence info. */
//    Int4     context;     /**< Context number of query */
//    GapEditScript* gap_info;/**< ALL gapped alignment is here */
// } BlastHSP;
// ```
pub struct BlastnHsp {
    pub identity: f64,
    pub length: usize,
    pub mismatch: usize,
    pub gapopen: usize,
    pub q_start: usize,
    pub q_end: usize,
    pub s_start: usize,
    pub s_end: usize,
    pub e_value: f64,
    pub bit_score: f64,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-143
    // ```c
    // typedef struct BlastHSP {
    //    Int4 score;
    //    Int4 num_ident;
    //    double bit_score;
    //    double evalue;
    //    ...
    //    Int4 num_positives;
    // } BlastHSP;
    // ```
    pub num_ident: usize,
    pub query_frame: i32,
    pub query_length: usize,
    pub q_idx: u32,
    pub s_idx: u32,
    pub raw_score: i32,
    pub internal_q_offset_0: usize,
    pub internal_q_end_0: usize,
    pub internal_s_offset_0: usize,
    pub internal_s_end_0: usize,
    pub internal_query_context_offset: i32,
    pub gap_info: Option<Vec<GapEditOp>>,
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-143
    // ```c
    // typedef struct BlastHSP {
    //    ...
    //    Int4 num_positives;
    // } BlastHSP;
    // ```
    pub num_positives: usize,
}

#[derive(Debug)]
// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
// ```c
// typedef struct BlastHSPList {
//    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
//    Int4 query_index; /**< Index of the query which this HSPList corresponds to. */
//    BlastHSP** hsp_array;
//    Int4 hspcnt;
//    double best_evalue;
// } BlastHSPList;
// ```
pub struct BlastnHspList {
    pub oid: u32,
    pub query_index: u32,
    pub hsps: Vec<BlastnHsp>,
    pub best_evalue: f64,
}

/// The fields of NCBI's `BlastHSPList` that a hit list (`Blast_HitListUpdate`) reads: the
/// final HSP lists of the traceback (`BlastnHspList`) and the preliminary HSP lists of the
/// preliminary stage (`run.rs`).
pub trait HitListEntry {
    fn oid(&self) -> u32;
    fn hsp_count(&self) -> usize;
    fn best_evalue(&self) -> f64;
    /// `hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list)`.
    fn update_best_evalue(&mut self);
    /// `hsp_list->hsp_array[0]->score`.
    fn first_score(&self) -> Option<i32>;
    /// `Blast_HSPListSortByEvalue`.
    fn sort_by_evalue(&mut self);
}

impl HitListEntry for BlastnHspList {
    fn oid(&self) -> u32 {
        self.oid
    }

    fn hsp_count(&self) -> usize {
        self.hsps.len()
    }

    fn best_evalue(&self) -> f64 {
        self.best_evalue
    }

    fn update_best_evalue(&mut self) {
        update_best_evalue(self);
    }

    fn first_score(&self) -> Option<i32> {
        self.hsps.first().map(|h| h.raw_score)
    }

    fn sort_by_evalue(&mut self) {
        sort_hsplist_by_evalue(self);
    }
}

#[derive(Debug)]
// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:168-180
// ```c
// typedef struct BlastHitList {
//    Int4 hsplist_count; /**< Filled size of the HSP lists array */
//    Int4 hsplist_max; /**< Maximal allowed size of the HSP lists array */
//    double worst_evalue; /**< Highest of the best e-values among the HSP lists */
//    Int4 low_score; /**< The lowest of the best scores among the HSP lists */
//    Boolean heapified; /**< Is this hit list already heapified? */
//    BlastHSPList** hsplist_array; /**< Array of HSP lists for individual database hits */
//    Int4 hsplist_current; /**< Number of allocated HSP list arrays. */
//    Int4 num_hits; /**< Number of similar hits for the query (for mapping) */
// } BlastHitList;
// ```
pub struct HitList<L> {
    pub hsplist_count: usize,
    pub hsplist_max: usize,
    pub worst_evalue: f64,
    pub low_score: i32,
    pub heapified: bool,
    pub hsplist_array: Vec<L>,
    pub hsplist_current: usize,
    pub num_hits: usize,
}

/// The hit list of a query's final HSP lists.
pub type BlastnHitList = HitList<BlastnHspList>;

pub type BlastnHspCompare = fn(&BlastnHsp, &BlastnHsp) -> Ordering;

impl BlastnHsp {
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-148
    // ```c
    // typedef struct BlastHSP {
    //    Int4 score;
    //    Int4 num_ident;
    //    double bit_score;
    //    double evalue;
    //    BlastSeg query;
    //    BlastSeg subject;
    //    Int4 context;
    //    GapEditScript* gap_info;
    // } BlastHSP;
    // ```
    pub fn from_hit(hit: Hit) -> Self {
        let Hit {
            identity,
            length,
            mismatch,
            gapopen,
            q_start,
            q_end,
            s_start,
            s_end,
            e_value,
            bit_score,
            num_ident,
            query_frame,
            query_length,
            q_idx,
            s_idx,
            raw_score,
            gap_info,
            num_positives,
            ..
        } = hit;
        Self {
            identity,
            length,
            mismatch,
            gapopen,
            q_start,
            q_end,
            s_start,
            s_end,
            e_value,
            bit_score,
            num_ident,
            query_frame,
            query_length,
            q_idx,
            s_idx,
            raw_score,
            internal_q_offset_0: if query_length > 0 && query_frame < 0 {
                query_length.saturating_sub(q_end)
            } else {
                q_start.saturating_sub(1)
            },
            internal_q_end_0: if query_length > 0 && query_frame < 0 {
                query_length.saturating_sub(q_start).saturating_add(1)
            } else {
                q_end
            },
            internal_s_offset_0: s_start.min(s_end).saturating_sub(1),
            internal_s_end_0: s_start.max(s_end),
            internal_query_context_offset: 0,
            gap_info,
            num_positives,
        }
    }

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to. */
    // } BlastHSPList;
    // ```
    pub fn into_hit(self) -> Hit {
        Hit {
            identity: self.identity,
            length: self.length,
            mismatch: self.mismatch,
            gapopen: self.gapopen,
            q_start: self.q_start,
            q_end: self.q_end,
            s_start: self.s_start,
            s_end: self.s_end,
            e_value: self.e_value,
            bit_score: self.bit_score,
            num_ident: self.num_ident,
            query_frame: self.query_frame,
            query_length: self.query_length,
            q_idx: self.q_idx,
            s_idx: self.s_idx,
            raw_score: self.raw_score,
            // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-143
            // ```c
            // typedef struct BlastHSP {
            //    BlastSeg query;
            //    BlastSeg subject;
            // } BlastHSP;
            // ```
            sort_query_offset: 0,
            sort_query_end: 0,
            sort_subject_offset: 0,
            sort_subject_end: 0,
            has_sort_offsets: false,
            gap_info: self.gap_info,
            num_positives: self.num_positives,
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:43-70 (GetPrelimHitlistSize)
// ```c
// Int4
// GetPrelimHitlistSize(Int4 hitlist_size, Int4 compositionBasedStats, Boolean gapped_calculation)
// {
//     Int4 prelim_hitlist_size = hitlist_size;
//     char * ADAPTIVE_CBS_ENV = getenv("ADAPTIVE_CBS");
//     if (compositionBasedStats) {
//         if(ADAPTIVE_CBS_ENV != NULL) {
//             if(hitlist_size < 1000) {
//                 prelim_hitlist_size = MAX(prelim_hitlist_size + 1000, 1500);
//             }
//             else {
//                 prelim_hitlist_size = prelim_hitlist_size*2 + 50;
//             }
//         }
//         else {
//             if(hitlist_size <= 500) {
//                 prelim_hitlist_size = 1050;
//             }
//             else {
//                 prelim_hitlist_size = prelim_hitlist_size*2 + 50;
//             }
//         }
//     }
//     else if (gapped_calculation) {
//          prelim_hitlist_size = MIN(MAX(2 * prelim_hitlist_size, 10),
//                                   prelim_hitlist_size + 50);
//     }
//     return prelim_hitlist_size;
// }
// ```
/// NCBI computes the size in `Int4` (blast_hits.c:44-46), whose sums the compiled NCBI
/// wraps: a hit list size from 2^30 to 2^31 - 51 gives 10 with gapped search, and a larger
/// one a negative size, with which NCBI crashes (LOSAT rejects it before the search,
/// `scoring.rs` `check_losat_limits`).
pub fn get_prelim_hitlist_size(
    hitlist_size: usize,
    composition_based_stats: bool,
    gapped_calculation: bool,
) -> i32 {
    // The argument is an `int` (`CArg_Integer`); the default is 500.
    let hitlist_size = i32::try_from(hitlist_size).unwrap_or(i32::MAX);
    let mut prelim_hitlist_size = hitlist_size;
    let adaptive_cbs = std::env::var_os("ADAPTIVE_CBS").is_some();
    if composition_based_stats {
        if adaptive_cbs {
            if hitlist_size < 1000 {
                prelim_hitlist_size = std::cmp::max(prelim_hitlist_size.wrapping_add(1000), 1500);
            } else {
                prelim_hitlist_size = prelim_hitlist_size.wrapping_mul(2).wrapping_add(50);
            }
        } else if hitlist_size <= 500 {
            prelim_hitlist_size = 1050;
        } else {
            prelim_hitlist_size = prelim_hitlist_size.wrapping_mul(2).wrapping_add(50);
        }
    } else if gapped_calculation {
        prelim_hitlist_size = std::cmp::min(
            std::cmp::max(prelim_hitlist_size.wrapping_mul(2), 10),
            prelim_hitlist_size.wrapping_add(50),
        );
    }
    prelim_hitlist_size
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1389-1403
// ```c
// static int s_EvalueComp(double evalue1, double evalue2)
// {
//     const double epsilon = 1.0e-180;
//     if (evalue1 < epsilon && evalue2 < epsilon) { return 0; }
//     if (evalue1 < evalue2) return -1;
//     else if (evalue1 > evalue2) return 1;
//     else return 0;
// }
// ```
pub fn evalue_comp(evalue1: f64, evalue2: f64) -> Ordering {
    const EPSILON: f64 = 1.0e-180;
    if evalue1 < EPSILON && evalue2 < EPSILON {
        Ordering::Equal
    } else if evalue1 < evalue2 {
        Ordering::Less
    } else if evalue1 > evalue2 {
        Ordering::Greater
    } else {
        Ordering::Equal
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1353
// ```c
// int ScoreCompareHSPs(const void* h1, const void* h2) {
//    if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
//        0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
//        0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
//        0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
//        result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
//    }
// }
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1122-1132
// ```c
// if (hsp->query.frame != hsp->subject.frame) {
//    *q_end = query_length - hsp->query.offset;
//    *q_start = *q_end - hsp->query.end + hsp->query.offset + 1;
// }
// ```
pub fn score_compare_hsps(a: &BlastnHsp, b: &BlastnHsp) -> Ordering {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1353
    // ```c
    // if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
    //     0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
    //     0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
    //     0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
    //     result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
    // }
    // ```
    // ScoreCompareHSPs compares canonical internal offsets/endpoints directly.
    let query_offsets = |h: &BlastnHsp| (h.internal_q_offset_0, h.internal_q_end_0);
    let (a_q_offset, a_q_end) = query_offsets(a);
    let (b_q_offset, b_q_end) = query_offsets(b);
    let (a_s_offset, a_s_end) = (a.internal_s_offset_0, a.internal_s_end_0);
    let (b_s_offset, b_s_end) = (b.internal_s_offset_0, b.internal_s_end_0);

    match b.raw_score.cmp(&a.raw_score) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match a_s_offset.cmp(&b_s_offset) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match b_s_end.cmp(&a_s_end) {
        Ordering::Equal => {}
        ord => return ord,
    }
    match a_q_offset.cmp(&b_q_offset) {
        Ordering::Equal => {}
        ord => return ord,
    }
    b_q_end.cmp(&a_q_end)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1414-1435
// ```c
// static int
// s_EvalueCompareHSPs(const void* v1, const void* v2)
// {
//    BlastHSP* h1,* h2;
//    int retval = 0;
//    ...
//    if ((retval = s_EvalueComp(h1->evalue, h2->evalue)) != 0)
//       return retval;
//    return ScoreCompareHSPs(v1, v2);
// }
// ```
pub fn evalue_compare_hsps(a: &BlastnHsp, b: &BlastnHsp) -> Ordering {
    let cmp = evalue_comp(a.e_value, b.e_value);
    if cmp == Ordering::Equal {
        score_compare_hsps(a, b)
    } else {
        cmp
    }
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1379-1381
// ```c
// qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//       ScoreCompareHSPs);
// ```
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1453-1454
// ```c
// qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//       s_EvalueCompareHSPs);
// ```
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2478
// ```c
// qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryOffsetCompareHSPs);
// ```
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2504
// ```c
// qsort(hsp_array, hsp_count, sizeof(BlastHSP*), s_QueryEndCompareHSPs);
// ```
// Preserve NCBI's comparator and invocation timing while using one stable,
// target-neutral Rust sort. Comparator-equal HSPs retain their incoming order.
pub fn qsort_blastn_hsps_by(hsps: &mut [BlastnHsp], compare: BlastnHspCompare) {
    if hsps.len() <= 1 {
        return;
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1379-1381
    // ```c
    // qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
    //       ScoreCompareHSPs);
    // ```
    hsps.sort_by(compare);
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1374-1382
// ```c
// void Blast_HSPListSortByScore(BlastHSPList* hsp_list)
// {
//     if (!hsp_list || hsp_list->hspcnt <= 1)
//         return;
//     if (!Blast_HSPListIsSortedByScore(hsp_list)) {
//         qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//               ScoreCompareHSPs);
//     }
// }
// ```
pub fn sort_hsps_by_score(hsps: &mut [BlastnHsp]) {
    if hsps.len() <= 1 {
        return;
    }
    let sorted = hsps
        .windows(2)
        .all(|pair| score_compare_hsps(&pair[0], &pair[1]) != Ordering::Greater);
    if !sorted {
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1379-1381
        // ```c
        // qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
        //       ScoreCompareHSPs);
        // ```
        qsort_blastn_hsps_by(hsps, score_compare_hsps);
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1437-1455
// ```c
// void Blast_HSPListSortByEvalue(BlastHSPList* hsp_list)
// {
//     if (hsp_list->hspcnt > 1) {
//         Int4 index;
//         BlastHSP** hsp_array = hsp_list->hsp_array;
//         for (index = 0; index < hsp_list->hspcnt - 1; ++index) {
//             if (s_EvalueCompareHSPs(&hsp_array[index],
//                                     &hsp_array[index+1]) > 0) {
//                 break;
//             }
//         }
//         if (index < hsp_list->hspcnt - 1) {
//             qsort(hsp_list->hsp_array, hsp_list->hspcnt,
//                   sizeof(BlastHSP*), s_EvalueCompareHSPs);
//         }
//     }
// }
// ```
pub fn sort_hsplist_by_evalue(list: &mut BlastnHspList) {
    if list.hsps.len() > 1 {
        let mut index = 0usize;
        while index < list.hsps.len() - 1 {
            if evalue_compare_hsps(&list.hsps[index], &list.hsps[index + 1]) == Ordering::Greater {
                break;
            }
            index += 1;
        }
        if index < list.hsps.len() - 1 {
            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1453-1454
            // ```c
            // qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
            //       s_EvalueCompareHSPs);
            // ```
            qsort_blastn_hsps_by(&mut list.hsps, evalue_compare_hsps);
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1739-1748
// ```c
// static double s_BlastGetBestEvalue(const BlastHSPList* hsp_list)
// {
//     double best_evalue = (double) INT4_MAX;
//     for (index=0; index<hsp_list->hspcnt; index++)
//         best_evalue = MIN(hsp_list->hsp_array[index]->evalue, best_evalue);
//     return best_evalue;
// }
// ```
pub fn update_best_evalue(list: &mut BlastnHspList) {
    let mut best = i32::MAX as f64;
    for hsp in &list.hsps {
        if hsp.e_value < best {
            best = hsp.e_value;
        }
    }
    list.best_evalue = best;
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2049-2067
// ```c
// Int2 Blast_TrimHSPListByMaxHsps(BlastHSPList* hsp_list,
//                                const BlastHitSavingOptions* hit_options)
// {
//    if ((hsp_list == NULL) ||
//        (hit_options->max_hsps_per_subject == 0) ||
//        (hsp_list->hspcnt <= hit_options->max_hsps_per_subject))
//       return 0;
//    hsp_list->hspcnt = hsp_max;
// }
// ```
pub fn trim_by_max_hsps(list: &mut BlastnHspList, max_hsps_per_subject: usize) {
    if max_hsps_per_subject == 0 || list.hsps.len() <= max_hsps_per_subject {
        return;
    }
    list.hsps.truncate(max_hsps_per_subject);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3077-3106
// ```c
// static int s_EvalueCompareHSPLists(const void* v1, const void* v2)
// {
//    if (h1->hspcnt == 0 && h2->hspcnt == 0) return 0;
//    else if (h1->hspcnt == 0) return 1;
//    else if (h2->hspcnt == 0) return -1;
//    if ((retval = s_EvalueComp(h1->best_evalue, h2->best_evalue)) != 0)
//       return retval;
//    if (h1->hsp_array[0]->score > h2->hsp_array[0]->score) return -1;
//    if (h1->hsp_array[0]->score < h2->hsp_array[0]->score) return 1;
//    return BLAST_CMP(h2->oid, h1->oid);
// }
// ```
pub fn compare_hsp_lists<L: HitListEntry>(a: &L, b: &L) -> Ordering {
    if a.hsp_count() == 0 && b.hsp_count() == 0 {
        return Ordering::Equal;
    } else if a.hsp_count() == 0 {
        return Ordering::Greater;
    } else if b.hsp_count() == 0 {
        return Ordering::Less;
    }

    let cmp = evalue_comp(a.best_evalue(), b.best_evalue());
    if cmp != Ordering::Equal {
        return cmp;
    }
    let a_score = a.first_score().unwrap_or(0);
    let b_score = b.first_score().unwrap_or(0);
    match b_score.cmp(&a_score) {
        Ordering::Equal => {}
        ord => return ord,
    }
    b.oid().cmp(&a.oid())
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1627-1650
// ```c
// static void
// s_Heapify (char* base0, char* base, char* lim, char* last,
//            size_t width, int (*compar )(const void*, const void* ))
// {
//    ...
//    left_son = base0 + 2*(base-base0) + width;
//    while (base <= lim) {
//       ...
//       large_son = (*compar)(left_son, left_son+width) >= 0 ?
//          left_son : left_son+width;
//       if ((*compar)(base, large_son) < 0) {
//          ...
//       } else
//          break;
//    }
// }
// ```
fn heapify_hsplist_array<L: HitListEntry>(lists: &mut [L], start: usize, end: usize) {
    let mut root = start;
    loop {
        let left = root.saturating_mul(2).saturating_add(1);
        if left > end {
            break;
        }
        let mut large = left;
        let right = left + 1;
        if right <= end && compare_hsp_lists(&lists[left], &lists[right]) == Ordering::Less {
            large = right;
        }
        if compare_hsp_lists(&lists[root], &lists[large]) == Ordering::Less {
            lists.swap(root, large);
            root = large;
        } else {
            break;
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1659-1676
// ```c
// static void
// s_CreateHeap (void* b, size_t nel, size_t width,
//    int (*compar )(const void*, const void* ))
// {
//    if (nel < 2)
//       return;
//    ...
//    i = nel/2;
//    for (base = &base0[(i - 1)*width]; i > 0; base = base - width) {
//       s_Heapify(base0, base, lim, basef, width, compar);
//       i--;
//    }
// }
// ```
fn create_hsplist_heap<L: HitListEntry>(lists: &mut [L]) {
    let nel = lists.len();
    if nel < 2 {
        return;
    }
    let mut i = nel / 2;
    while i > 0 {
        i -= 1;
        heapify_hsplist_array(lists, i, nel - 1);
    }
}

impl<L: HitListEntry> HitList<L> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3125-3133
    // ```c
    // BlastHitList* Blast_HitListNew(Int4 hitlist_size)
    // {
    //    BlastHitList* new_hitlist = (BlastHitList*) calloc(1, sizeof(BlastHitList));
    //    new_hitlist->hsplist_max = hitlist_size;
    //    new_hitlist->low_score = INT4_MAX;
    //    new_hitlist->hsplist_count = 0;
    //    new_hitlist->hsplist_current = 0;
    //    return new_hitlist;
    // }
    // ```
    pub fn new(hitlist_size: usize) -> Self {
        Self {
            hsplist_count: 0,
            hsplist_max: hitlist_size,
            worst_evalue: 0.0,
            low_score: i32::MAX,
            heapified: false,
            hsplist_array: Vec::new(),
            hsplist_current: 0,
            num_hits: 0,
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3219-3240
    // ```c
    // static Int2 s_Blast_HitListGrowHSPListArray(BlastHitList* hit_list)
    // {
    //     const int kStartValue = 100;
    //     if (hit_list->hsplist_current >= hit_list->hsplist_max)
    //        return 1;
    //     if (hit_list->hsplist_current <= 0)
    //        hit_list->hsplist_current = kStartValue;
    //     else
    //        hit_list->hsplist_current =
    //           MIN(2*hit_list->hsplist_current, hit_list->hsplist_max);
    //     hit_list->hsplist_array =
    //        (BlastHSPList**) realloc(hit_list->hsplist_array,
    //                                 hit_list->hsplist_current*sizeof(BlastHSPList*));
    //     if (hit_list->hsplist_array == NULL)
    //        return BLASTERR_MEMORY;
    //     return 0;
    // }
    // ```
    fn grow_hsplist_array(&mut self) -> bool {
        const K_START_VALUE: usize = 100;
        if self.hsplist_current >= self.hsplist_max {
            return false;
        }
        if self.hsplist_current == 0 {
            self.hsplist_current = K_START_VALUE.min(self.hsplist_max);
        } else {
            let next = self.hsplist_current.saturating_mul(2);
            self.hsplist_current = next.min(self.hsplist_max);
        }
        if self.hsplist_array.capacity() < self.hsplist_current {
            self.hsplist_array
                .reserve(self.hsplist_current - self.hsplist_array.capacity());
        }
        true
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3196-3209
    // ```c
    // static void
    // s_BlastHitListInsertHSPListInHeap(BlastHitList* hit_list,
    //                                  BlastHSPList* hsp_list)
    // {
    //       Blast_HSPListFree(hit_list->hsplist_array[0]);
    //       hit_list->hsplist_array[0] = hsp_list;
    //       if (hit_list->hsplist_count >= 2) {
    //          s_Heapify((char*)hit_list->hsplist_array, (char*)hit_list->hsplist_array,
    //                  (char*)&hit_list->hsplist_array[hit_list->hsplist_count/2 - 1],
    //                  (char*)&hit_list->hsplist_array[hit_list->hsplist_count-1],
    //                  sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
    //       }
    //       hit_list->worst_evalue = hit_list->hsplist_array[0]->best_evalue;
    //       hit_list->low_score = hit_list->hsplist_array[0]->hsp_array[0]->score;
    // }
    // ```
    fn insert_hsplist_in_heap(&mut self, hsp_list: L) {
        if self.hsplist_array.is_empty() {
            self.hsplist_array.push(hsp_list);
            self.hsplist_count = 1;
        } else {
            self.hsplist_array[0] = hsp_list;
        }
        if self.hsplist_count >= 2 {
            heapify_hsplist_array(&mut self.hsplist_array, 0, self.hsplist_count - 1);
        }
        if let Some(root) = self.hsplist_array.first() {
            self.worst_evalue = root.best_evalue();
            if let Some(score) = root.first_score() {
                self.low_score = score;
            }
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3243-3297
    // ```c
    // Int2 Blast_HitListUpdate(BlastHitList* hit_list, BlastHSPList* hsp_list)
    // {
    //    hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
    //    if (hit_list->hsplist_count < hit_list->hsplist_max) {
    //       if (hit_list->hsplist_current == hit_list->hsplist_count)
    //       {
    //          Int2 status = s_Blast_HitListGrowHSPListArray(hit_list);
    //          if (status)
    //            return status;
    //       }
    //       hit_list->hsplist_array[hit_list->hsplist_count++] = hsp_list;
    //       hit_list->worst_evalue =
    //          MAX(hsp_list->best_evalue, hit_list->worst_evalue);
    //       hit_list->low_score =
    //          MIN(hsp_list->hsp_array[0]->score, hit_list->low_score);
    //    } else {
    //       if (!hit_list->heapified) {
    //           for (index =0; index < hit_list->hsplist_count; index++) {
    //               Blast_HSPListSortByEvalue(hit_list->hsplist_array[index]);
    //               hit_list->hsplist_array[index]->best_evalue =
    //                   s_BlastGetBestEvalue(hit_list->hsplist_array[index]);
    //           }
    //           s_CreateHeap(hit_list->hsplist_array, hit_list->hsplist_count,
    //                        sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
    //           hit_list->heapified = TRUE;
    //       }
    //       Blast_HSPListSortByEvalue(hsp_list);
    //       hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list);
    //       evalue_order = s_EvalueCompareHSPLists(&(hit_list->hsplist_array[0]), &hsp_list);
    //       if (evalue_order < 0) {
    //          Blast_HSPListFree(hsp_list);
    //       } else {
    //          s_BlastHitListInsertHSPListInHeap(hit_list, hsp_list);
    //       }
    //    }
    //    return 0;
    // }
    // ```
    pub fn update(&mut self, mut hsp_list: L) {
        hsp_list.update_best_evalue();

        if self.hsplist_count < self.hsplist_max {
            if self.hsplist_current == self.hsplist_count && !self.grow_hsplist_array() {
                return;
            }
            self.hsplist_array.push(hsp_list);
            self.hsplist_count += 1;
            self.worst_evalue = self
                .worst_evalue
                .max(self.hsplist_array.last().unwrap().best_evalue());
            if let Some(score) = self.hsplist_array.last().unwrap().first_score() {
                self.low_score = self.low_score.min(score);
            }
        } else {
            if !self.heapified {
                for list in &mut self.hsplist_array {
                    list.sort_by_evalue();
                    list.update_best_evalue();
                }
                create_hsplist_heap(&mut self.hsplist_array);
                self.heapified = true;
            }

            hsp_list.sort_by_evalue();
            hsp_list.update_best_evalue();
            let evalue_order = compare_hsp_lists(&self.hsplist_array[0], &hsp_list);
            if evalue_order == Ordering::Less {
                return;
            }
            self.insert_hsplist_in_heap(hsp_list);
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3172-3187
    // ```c
    // static void s_BlastHitListPurge(BlastHitList* hit_list)
    // {
    //    if (!hit_list) return;
    //    hsplist_count = hit_list->hsplist_count;
    //    for (index = 0; index < hsplist_count &&
    //            hit_list->hsplist_array[index]->hspcnt > 0; ++index);
    //    hit_list->hsplist_count = index;
    //    for ( ; index < hsplist_count; ++index) {
    //       Blast_HSPListFree(hit_list->hsplist_array[index]);
    //    }
    // }
    // ```
    fn purge(&mut self) {
        let mut index = 0usize;
        while index < self.hsplist_count {
            if self.hsplist_array[index].hsp_count() == 0 {
                break;
            }
            index += 1;
        }
        self.hsplist_count = index;
        self.hsplist_array.truncate(index);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3331-3337
    // ```c
    // Int2 Blast_HitListSortByEvalue(BlastHitList* hit_list)
    // {
    //    if (hit_list && hit_list->hsplist_count > 1) {
    //       qsort(hit_list->hsplist_array, hit_list->hsplist_count,
    //             sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
    //    }
    //    s_BlastHitListPurge(hit_list);
    //    return 0;
    // }
    // ```
    pub fn sort_by_evalue(&mut self) {
        if self.hsplist_count > 1 {
            self.hsplist_array.sort_by(compare_hsp_lists);
        }
        self.purge();
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:877-892
    // ```c
    // static void s_BlastPruneExtraHits(BlastHSPResults* results, Int4 hitlist_size)
    // {
    //    for (subject_index = hitlist_size;
    //         subject_index < hit_list->hsplist_count; ++subject_index) {
    //       hit_list->hsplist_array[subject_index] =
    //       Blast_HSPListFree(hit_list->hsplist_array[subject_index]);
    //    }
    //    hit_list->hsplist_count = MIN(hit_list->hsplist_count, hitlist_size);
    // }
    // ```
    pub fn prune_by_size(&mut self, hitlist_size: usize) {
        if hitlist_size == 0 {
            self.hsplist_array.clear();
            self.hsplist_count = 0;
            return;
        }
        if self.hsplist_count > hitlist_size {
            self.hsplist_array.truncate(hitlist_size);
            self.hsplist_count = hitlist_size;
        }
    }
}

/// Write BLASTN hit lists to an existing writer.
///
/// NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
/// ```c
/// void CBlastTabularInfo::Print()
/// {
///     ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
///         if (iter != m_FieldsToShow.begin())
///             m_Ostream << m_FieldDelimiter;
///         x_PrintField(*iter);
///     }
///     m_Ostream << "\n";
/// }
/// ```
#[allow(clippy::too_many_arguments)]
pub fn write_output_blastn_hitlists_to_writer<W: Write>(
    hit_lists: &[Option<BlastnHitList>],
    writer: &mut W,
    query_ids: &[Arc<str>],
    subject_ids: &[Arc<str>],
    output_format: BlastnOutputFormat,
    query_titles: &[Arc<str>],
    subject_title: &str,
    unsearched: &[bool],
    epilog: bool,
    mut probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
) -> io::Result<()> {
    let config = OutputConfig::ncbi_compat();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1411
    // ```c
    // CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
    // ```
    // The rows are printed in the order of the final HSP list, so a running count is
    // each row's HSP index.
    let mut hsp_index: HspIndex = 0;

    for (q_idx, hit_list_opt) in hit_lists.iter().enumerate() {
        // The query's warnings come before its lines (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        if output_format == BlastnOutputFormat::TabularWithComments {
            let query_title = query_titles
                .get(q_idx)
                .map(|title| title.as_ref())
                .unwrap_or("unknown");
            write_blastn_outfmt7_header(
                writer,
                query_title,
                subject_title,
                (!unsearched.get(q_idx).copied().unwrap_or(false))
                    .then(|| blastn_hsp_count(hit_list_opt.as_ref())),
            )?;
        }
        let hit_list = match hit_list_opt {
            Some(value) => value,
            None => continue,
        };
        for hsp_list in &hit_list.hsplist_array {
            for hsp in &hsp_list.hsps {
                let query_id = query_ids
                    .get(hsp.q_idx as usize)
                    .map(|id| id.as_ref())
                    .unwrap_or("unknown");
                let subject_id = subject_ids
                    .get(hsp.s_idx as usize)
                    .map(|id| id.as_ref())
                    .unwrap_or("unknown");
                // NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
                // ```c
                // ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
                //     // Add tab in front of field, except for the first field.
                //     if (iter != m_FieldsToShow.begin())
                //         m_Ostream << m_FieldDelimiter;
                //     x_PrintField(*iter);
                // }
                // m_Ostream << "\n";
                // ```
                if blastn_trace::should_trace_range(
                    "hitlist",
                    hsp.q_idx * 2 + if hsp.query_frame < 0 { 1 } else { 0 },
                    hsp.s_idx as usize,
                    subject_id,
                    hsp.internal_q_offset_0,
                    hsp.internal_q_end_0,
                    hsp.internal_s_offset_0,
                    hsp.internal_s_end_0,
                    hsp.query_length,
                    hsp.query_frame,
                ) {
                    blastn_trace::log(
                        "hitlist",
                        format!(
                            "query={} subject={} q={}..{} s={}..{} raw_score={} bit_score={:.12} evalue={:.12e} length={} identities={} mismatches={} gapopen={}",
                            query_id,
                            subject_id,
                            hsp.q_start,
                            hsp.q_end,
                            hsp.s_start,
                            hsp.s_end,
                            hsp.raw_score,
                            hsp.bit_score,
                            hsp.e_value,
                            hsp.length,
                            hsp.num_ident,
                            hsp.mismatch,
                            hsp.gapopen
                        ),
                    );
                }
                // NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
                // ```c
                // void CBlastTabularInfo::Print()
                // {
                //     ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
                //         if (iter != m_FieldsToShow.begin())
                //             m_Ostream << m_FieldDelimiter;
                //         x_PrintField(*iter);
                //     }
                //     m_Ostream << "\n";
                // }
                // ```
                // One printed row is one HSP; the probe marks it without changing it.
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.begin(hsp_index);
                }
                write_hit_fields(
                    writer,
                    query_id,
                    subject_id,
                    hsp.identity,
                    hsp.num_ident,
                    hsp.length,
                    hsp.mismatch,
                    hsp.gapopen,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                    hsp.e_value,
                    hsp.bit_score,
                    &config,
                )?;
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.end(hsp_index);
                }
                hsp_index += 1;
            }
        }
    }

    // `PrintEpilog` writes the last line; an error in a later query batch skips it.
    if output_format == BlastnOutputFormat::TabularWithComments && epilog {
        // NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1324
        // ```c
        // m_Ostream << "# BLAST processed " << num_queries << " queries\n";
        // ```
        writeln!(writer, "# BLAST processed {} queries", hit_lists.len())?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-148
    // ```c
    // typedef struct BlastHSP {
    //    Int4 score;
    //    double evalue;
    //    BlastSeg query;
    //    BlastSeg subject;
    //    Int4 context;
    // } BlastHSP;
    // ```
    fn make_hsp(e_value: f64, raw_score: i32, s_idx: u32) -> BlastnHsp {
        BlastnHsp {
            identity: 0.0,
            length: 0,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end: 1,
            s_start: 1,
            s_end: 1,
            e_value,
            bit_score: 0.0,
            num_ident: 0,
            query_frame: 1,
            query_length: 1,
            q_idx: 0,
            s_idx,
            raw_score,
            internal_q_offset_0: 0,
            internal_q_end_0: 1,
            internal_s_offset_0: 0,
            internal_s_end_0: 1,
            internal_query_context_offset: 0,
            gap_info: None,
            num_positives: 0,
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1353
    // ```c
    // if (0 == (result = BLAST_CMP(hsp2->score, hsp1->score)) && ...) {
    //     result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
    // }
    // return result;
    // ```
    #[test]
    fn stable_hsp_sort_preserves_comparator_equal_input_order() {
        let mut hsps = vec![
            make_hsp(1.0, 100, 7),
            make_hsp(1.0, 100, 3),
            make_hsp(1.0, 100, 9),
        ];

        qsort_blastn_hsps_by(&mut hsps, |_, _| Ordering::Equal);

        assert_eq!(
            hsps.iter().map(|hsp| hsp.s_idx).collect::<Vec<_>>(),
            vec![7, 3, 9]
        );
    }

    #[test]
    fn test_prelim_hitlist_size_gapped() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:43-70
        // ```c
        // else if (gapped_calculation) {
        //      prelim_hitlist_size = MIN(MAX(2 * prelim_hitlist_size, 10),
        //                               prelim_hitlist_size + 50);
        // }
        // ```
        assert_eq!(get_prelim_hitlist_size(1, false, true), 10);
        // NCBI's Int4 arithmetic wraps (blast_hits.c:68).
        assert_eq!(get_prelim_hitlist_size(1 << 30, false, true), 10);
        assert_eq!(
            get_prelim_hitlist_size((1 << 30) - 1, false, true),
            (1 << 30) + 49
        );
        assert!(get_prelim_hitlist_size(i32::MAX as usize - 49, false, true) < 0);
        assert_eq!(get_prelim_hitlist_size(30, false, true), 60);
        assert_eq!(get_prelim_hitlist_size(1000, false, true), 1050);
    }

    #[test]
    fn test_hitlist_update_keeps_best() {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3243-3297
        // ```c
        // if (hit_list->hsplist_count < hit_list->hsplist_max) { ... }
        // else {
        //    if (!hit_list->heapified) { ... }
        //    evalue_order = s_EvalueCompareHSPLists(&(hit_list->hsplist_array[0]), &hsp_list);
        //    if (evalue_order < 0) { Blast_HSPListFree(hsp_list); }
        //    else { s_BlastHitListInsertHSPListInHeap(hit_list, hsp_list); }
        // }
        // ```
        let mut hit_list = BlastnHitList::new(1);

        let hsp_list_a = BlastnHspList {
            oid: 10,
            query_index: 0,
            hsps: vec![make_hsp(5.0, 50, 10)],
            best_evalue: i32::MAX as f64,
        };
        hit_list.update(hsp_list_a);
        assert_eq!(hit_list.hsplist_count, 1);
        assert_eq!(hit_list.hsplist_array[0].oid, 10);

        let hsp_list_b = BlastnHspList {
            oid: 20,
            query_index: 0,
            hsps: vec![make_hsp(1.0, 80, 20)],
            best_evalue: i32::MAX as f64,
        };
        hit_list.update(hsp_list_b);
        assert_eq!(hit_list.hsplist_count, 1);
        assert_eq!(hit_list.hsplist_array[0].oid, 20);

        let hsp_list_c = BlastnHspList {
            oid: 30,
            query_index: 0,
            hsps: vec![make_hsp(10.0, 40, 30)],
            best_evalue: i32::MAX as f64,
        };
        hit_list.update(hsp_list_c);
        assert_eq!(hit_list.hsplist_count, 1);
        assert_eq!(hit_list.hsplist_array[0].oid, 20);
    }

    fn make_list(oid: u32, hsps: &[(f64, i32)]) -> BlastnHspList {
        BlastnHspList {
            oid,
            query_index: 0,
            hsps: hsps
                .iter()
                .map(|&(e_value, score)| make_hsp(e_value, score, oid))
                .collect(),
            best_evalue: i32::MAX as f64,
        }
    }

    fn list_oids(hit_list: &BlastnHitList) -> Vec<u32> {
        hit_list.hsplist_array.iter().map(|list| list.oid).collect()
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3095-3106
    // ```c
    //    if ((retval = s_EvalueComp(h1->best_evalue,
    //                                    h2->best_evalue)) != 0)
    //       return retval;
    //
    //    if (h1->hsp_array[0]->score > h2->hsp_array[0]->score)
    //       return -1;
    //    if (h1->hsp_array[0]->score < h2->hsp_array[0]->score)
    //       return 1;
    //
    //    /* In case of equal best E-values and scores, order will be determined
    //       by ordinal ids of the subject sequences */
    //    return BLAST_CMP(h2->oid, h1->oid);
    // ```
    // An equal list replaces the heap root (`evalue_order < 0` frees only a worse one,
    // blast_hits.c:3285-3291), so equal lists keep the higher ordinal ids, first.
    #[test]
    fn hit_list_equal_evalue_and_score_keeps_higher_oids() {
        let mut hit_list = BlastnHitList::new(3);
        for oid in 0..10 {
            hit_list.update(make_list(oid, &[(1e-10, 100)]));
        }
        assert!(hit_list.heapified);
        hit_list.sort_by_evalue();
        assert_eq!(list_oids(&hit_list), vec![9, 8, 7]);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1392-1395
    // ```c
    //     const double epsilon = 1.0e-180;
    //     if (evalue1 < epsilon && evalue2 < epsilon) {
    //         return 0;
    //     }
    // ```
    #[test]
    fn hit_list_evalues_below_1e180_compare_equal() {
        let mut hit_list = BlastnHitList::new(1);
        hit_list.update(make_list(1, &[(1e-200, 50)]));
        hit_list.update(make_list(2, &[(1e-190, 60)]));
        assert_eq!(list_oids(&hit_list), vec![2]);
        let mut hit_list = BlastnHitList::new(1);
        hit_list.update(make_list(1, &[(1e-100, 50)]));
        hit_list.update(make_list(2, &[(1e-90, 60)]));
        assert_eq!(list_oids(&hit_list), vec![1]);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3099-3102
    // ```c
    //    if (h1->hsp_array[0]->score > h2->hsp_array[0]->score)
    //       return -1;
    //    if (h1->hsp_array[0]->score < h2->hsp_array[0]->score)
    //       return 1;
    // ```
    #[test]
    fn hit_list_equal_evalue_higher_first_score_wins() {
        let mut hit_list = BlastnHitList::new(1);
        hit_list.update(make_list(5, &[(1e-5, 40)]));
        hit_list.update(make_list(1, &[(1e-5, 41)]));
        assert_eq!(list_oids(&hit_list), vec![1]);
        hit_list.update(make_list(9, &[(1e-5, 40)]));
        assert_eq!(list_oids(&hit_list), vec![1]);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3276-3287
    // ```c
    //           s_CreateHeap(hit_list->hsplist_array, hit_list->hsplist_count,
    //                        sizeof(BlastHSPList*), s_EvalueCompareHSPLists);
    //           hit_list->heapified = TRUE;
    //       }
    //       ...
    //       evalue_order = s_EvalueCompareHSPLists(&(hit_list->hsplist_array[0]), &hsp_list);
    //       if (evalue_order < 0) {
    // ```
    #[test]
    fn hit_list_heap_keeps_the_best_lists_in_comparator_order() {
        let keys: Vec<(u32, f64, i32)> = (0..40u32)
            .map(|oid| {
                let e_value = [1e-30, 1e-10, 1e-3, 0.5, 2.0][(oid * 7 % 5) as usize];
                (oid, e_value, 30 + (oid * 13 % 11) as i32)
            })
            .collect();
        let mut hit_list = BlastnHitList::new(7);
        for &(oid, e_value, score) in &keys {
            hit_list.update(make_list(oid, &[(e_value, score)]));
        }
        hit_list.sort_by_evalue();
        let mut all: Vec<BlastnHspList> = keys
            .iter()
            .map(|&(oid, e_value, score)| {
                let mut list = make_list(oid, &[(e_value, score)]);
                list.update_best_evalue();
                list
            })
            .collect();
        all.sort_by(compare_hsp_lists);
        let best: Vec<u32> = all.iter().take(7).map(|list| list.oid).collect();
        assert_eq!(list_oids(&hit_list), best);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3267-3283
    // ```c
    // if (!hit_list->heapified) {
    //     for (index =0; index < hit_list->hsplist_count; index++) {
    //         Blast_HSPListSortByEvalue(hit_list->hsplist_array[index]);
    //         hit_list->hsplist_array[index]->best_evalue =
    //             s_BlastGetBestEvalue(hit_list->hsplist_array[index]);
    //     }
    // ```
    #[test]
    fn hit_list_heapify_sorts_each_list_by_evalue() {
        let mut hit_list = BlastnHitList::new(1);
        hit_list.update(make_list(3, &[(1e-3, 10), (1e-20, 50)]));
        assert_eq!(hit_list.hsplist_array[0].hsps[0].raw_score, 10);
        hit_list.update(make_list(4, &[(1e-2, 70)]));
        assert!(hit_list.heapified);
        assert_eq!(list_oids(&hit_list), vec![3]);
        assert_eq!(hit_list.hsplist_array[0].hsps[0].raw_score, 50);
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3172-3187
    // ```c
    // for (index = 0; index < hsplist_count &&
    //         hit_list->hsplist_array[index]->hspcnt > 0; ++index);
    // hit_list->hsplist_count = index;
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:877-892
    // ```c
    // hit_list->hsplist_count = MIN(hit_list->hsplist_count, hitlist_size);
    // ```
    #[test]
    fn hit_list_sort_purges_empty_lists_and_prune_caps_the_size() {
        let mut hit_list = BlastnHitList::new(10);
        hit_list.update(make_list(1, &[]));
        hit_list.update(make_list(2, &[(1e-5, 50)]));
        hit_list.update(make_list(3, &[]));
        hit_list.update(make_list(4, &[(1e-6, 60)]));
        for oid in 5..9 {
            hit_list.update(make_list(oid, &[(1e-4 * f64::from(oid), 40)]));
        }
        hit_list.sort_by_evalue();
        assert_eq!(list_oids(&hit_list), vec![4, 2, 5, 6, 7, 8]);
        hit_list.prune_by_size(4);
        assert_eq!(list_oids(&hit_list), vec![4, 2, 5, 6]);
        assert_eq!(hit_list.hsplist_count, 4);
    }

    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1264-1284
    // ```c
    // x_PrintQueryAndDbNames(program_version, bioseq, dbname, rid, iteration,
    //                        subj_bioseq);
    // if (align_set) {
    //     int num_hits = align_set->Get().size();
    //     if (num_hits != 0) PrintFieldNames(is_csv);
    //     m_Ostream << "# " << num_hits << " hits found" << "\n";
    // }
    // ```
    #[test]
    fn test_blastn_outfmt7_header_matches_ncbi_shape() {
        let hit_lists = vec![None];
        let query_ids = vec![Arc::<str>::from("query")];
        let subject_ids = vec![Arc::<str>::from("subject")];
        let query_titles = vec![Arc::<str>::from("query full description")];
        let mut output = Vec::new();

        write_output_blastn_hitlists_to_writer(
            &hit_lists,
            &mut output,
            &query_ids,
            &subject_ids,
            BlastnOutputFormat::TabularWithComments,
            &query_titles,
            "User specified sequence set (Input: subject.fasta)",
            &[],
            true,
            None,
            None,
        )
        .unwrap();

        assert_eq!(
            String::from_utf8(output).unwrap(),
            "# BLASTN 2.17.0+\n# Query: query full description\n# Database: User specified sequence set (Input: subject.fasta)\n# 0 hits found\n# BLAST processed 1 queries\n"
        );
        // A query of a batch that NCBI did not search has no count.
        let mut output = Vec::new();
        write_output_blastn_hitlists_to_writer(
            &hit_lists,
            &mut output,
            &query_ids,
            &subject_ids,
            BlastnOutputFormat::TabularWithComments,
            &query_titles,
            "User specified sequence set (Input: subject.fasta)",
            &[true],
            false,
            None,
            None,
        )
        .unwrap();
        let text = String::from_utf8(output).unwrap();
        assert!(!text.contains("hits found"));
        // Without NCBI's epilog (an error in a later query batch), no last line.
        assert!(!text.contains("BLAST processed"));
    }

    #[test]
    fn test_blastn_rejects_custom_outfmt_until_fields_are_ported() {
        assert!(parse_blastn_output_format("6").is_ok());
        assert!(parse_blastn_output_format("7").is_ok());
        assert!(parse_blastn_output_format("6 qaccver saccver").is_err());
        // NCBI's parse (blast_args.cpp:2801-2851): white space is trimmed, the number is a
        // signed decimal, and a custom specification counts only for the tabular formats.
        for (spec, format) in [
            (" 6\t", BlastnOutputFormat::Tabular),
            ("+6", BlastnOutputFormat::Tabular),
            ("07", BlastnOutputFormat::TabularWithComments),
            ("0 qaccver", BlastnOutputFormat::Pairwise),
        ] {
            assert_eq!(
                parse_blastn_output_format(spec).unwrap(),
                format,
                "{spec:?}"
            );
        }
        let native = |spec: &str| {
            let error = parse_blastn_output_format(spec)
                .unwrap_err()
                .downcast::<NativeError>()
                .unwrap();
            (error.exit, error.message)
        };
        assert_eq!(
            native("6\u{a0}"),
            (1, "BLAST query/options error: '6\u{a0}' is not a valid output format\nPlease refer to the BLAST+ user manual.\n".to_string())
        );
        assert_eq!(
            native(" ").1.lines().next(),
            Some("BLAST query/options error: '' is not a valid output format")
        );
        for spec in ["22", "-1", "99"] {
            assert_eq!(
                native(spec),
                (
                    255,
                    "Error: Formatting choice is out of range\n".to_string()
                )
            );
        }
        // NCBI checks a delimiter before the format number (blast_args.cpp:2810-2825).
        for spec in ["0 delim", "abc delim", "99 delimiter qaccver"] {
            assert_eq!(native(spec).0, 1, "{spec:?}");
            assert!(
                native(spec).1.contains("Delimiter format is invalid"),
                "{spec:?}"
            );
        }
        for spec in ["0 delim=,", "6 delim=", "7 delim= "] {
            assert!(parse_blastn_output_format(spec).is_ok(), "{spec:?}");
        }
        for spec in ["5", "21", "6 delim=,", "6 delimiter=;", "7  qaccver"] {
            let error = parse_blastn_output_format(spec).unwrap_err().to_string();
            assert!(
                error.contains("not supported by LOSAT's BLASTN"),
                "{spec:?}: {error}"
            );
        }
    }
}
