use std::cmp::Ordering;
use std::io::{self, Write};
use std::sync::Arc;

use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::blastinput::app::ReportFormat;
use crate::common::{GapEditOp, Hit};
use crate::report::outfmt6::write_hit_fields_bytes;
use crate::report::OutputConfig;

use super::tracing as blastn_trace;

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:770-782
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1410-1414
// The pairwise report (outfmt 0) is the non-tabular branch of the same formatter, as are the
// query-anchored formats 1 to 4, which differ from it only in their alignment section
// (`report::pairwise::PairwiseView`, `report::query_anchored`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlastnOutputFormat {
    /// Outfmt 0 to 4: the pairwise report with outfmt 0's alignments or a query-anchored section.
    Pairwise(crate::report::pairwise::PairwiseView),
    Tabular,
    TabularWithComments,
    /// The tabular reports with a field list or a delimiter (outfmt 6 and 7 other than the
    /// default `std` fields with a tab) and outfmt 10 and 20 (`report::tabular`, through
    /// `algorithm::blastn::tabular`). The default outfmt 6 and 7 keep `Tabular` and
    /// `TabularWithComments`, whose writers this report never enters.
    CustomTabular(&'static crate::report::tabular::TabularSpec),
    /// The Seq-annot of Seq-aligns of each query: outfmt 8 (ASN.1 text) or 12 (JSON)
    /// (`report::serial::write_query`, through `algorithm::blastn::seqalign`).
    Seqalign(crate::report::serial::SerialFormat),
    /// BLAST XML (outfmt 5; `report::blast_xml`, through `algorithm::blastn::blast_xml`).
    Xml,
    /// The XML2 / JSON documents of outfmt 16, 15, 14 and 13 (`report::blast_xml2`, through the same
    /// `algorithm::blastn::blast_xml`).
    Xml2(crate::report::blast_xml2::Xml2Kind),
    /// The BLAST archive (outfmt 11; `report::archive`, through `algorithm::blastn::archive`).
    Archive,
    /// SAM (outfmt 17; `report::sam`, through `algorithm::blastn::sam`).
    Sam(crate::report::sam::SamSpec),
}

impl BlastnOutputFormat {
    /// The BLASTN report of a format that `report_format` accepts for BLASTN
    /// (`Program::writes`; NCBI's `-outfmt` reading is `blastinput/app.rs`
    /// `parse_formatting_string`).
    ///
    /// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_args.hpp:1025-1040
    /// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_args.hpp:1041-1050
    ///
    /// `choice` is the format choice: SAM reads its custom specification (`SQ`, `SR`), and the
    /// tabular formats with a field list or a delimiter, and outfmt 10 and 20, are the custom
    /// tabular report (`report::tabular::TabularSpec::of_app`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:770-786
    pub(crate) fn of(
        format: ReportFormat,
        choice: &crate::blastinput::app::FormatChoice,
    ) -> anyhow::Result<Self> {
        if let Some(spec) = crate::report::tabular::TabularSpec::of_app(
            format,
            choice,
            crate::report::archive::options::RequestKind::Blastn,
        ) {
            return Ok(BlastnOutputFormat::CustomTabular(spec.intern()));
        }
        if let Some(view) = crate::report::pairwise::PairwiseView::of(format) {
            return Ok(BlastnOutputFormat::Pairwise(view));
        }
        match format {
            ReportFormat::Sam => Ok(BlastnOutputFormat::Sam(crate::report::sam::SamSpec::parse(
                &choice.spec,
            ))),
            ReportFormat::Tabular => Ok(BlastnOutputFormat::Tabular),
            ReportFormat::TabularWithComments => Ok(BlastnOutputFormat::TabularWithComments),
            ReportFormat::Xml => Ok(BlastnOutputFormat::Xml),
            ReportFormat::ArchiveFormat => Ok(BlastnOutputFormat::Archive),
            // 16, 15, 14 and 13 (the last two write files next to `-out`: only the command line
            // reaches them, `run_local` has no output file).
            ReportFormat::Xml2S | ReportFormat::JsonS | ReportFormat::Xml2 | ReportFormat::Json => {
                Ok(BlastnOutputFormat::Xml2(
                    crate::report::blast_xml2::Xml2Kind::of(format)
                        .expect("an XML2 or JSON document format"),
                ))
            }
            other => match crate::report::serial::SerialFormat::of(other) {
                Some(serial) => Ok(BlastnOutputFormat::Seqalign(serial)),
                None => anyhow::bail!(
                    "output format {} is not supported by LOSAT's BLASTN",
                    other.number()
                ),
            },
        }
    }

    /// The Seq-annot format of outfmt 8 and 12.
    pub(crate) fn serial(self) -> Option<crate::report::serial::SerialFormat> {
        match self {
            BlastnOutputFormat::Seqalign(serial) => Some(serial),
            _ => None,
        }
    }

    /// The SAM document of outfmt 17.
    pub(crate) fn sam_spec(self) -> Option<crate::report::sam::SamSpec> {
        match self {
            BlastnOutputFormat::Sam(spec) => Some(spec),
            _ => None,
        }
    }

    /// Whether the report is a structured XML / JSON document of the queries' final hit lists:
    /// BLAST XML (outfmt 5) or the XML2 / JSON documents (outfmt 15 and 16), which read the same
    /// per-query records and the subjects' letters.
    pub(crate) fn is_structured_xml(self) -> bool {
        matches!(self, BlastnOutputFormat::Xml | BlastnOutputFormat::Xml2(_))
    }

    /// The XML2 / JSON document of outfmt 13 to 16.
    pub(crate) fn xml2_kind(self) -> Option<crate::report::blast_xml2::Xml2Kind> {
        match self {
            BlastnOutputFormat::Xml2(kind) => Some(kind),
            _ => None,
        }
    }

    /// What the format does with a subject whose lookup throws where its batch's prefetch is off
    /// (`report::db_subject::FormatThrow::of`): the default tables write the row (their `saccver`
    /// is a field of the subject block), a table with a field list as its fields say
    /// (`TabularSpec::format_throw`), the XML2 and JSON documents as their form says; the
    /// Seq-annots and the archive read nothing of the subject; LOSAT refuses it in the others.
    pub(crate) fn format_throw(&self) -> crate::report::db_subject::FormatThrow {
        use crate::report::db_subject::{AfterThrow, FormatThrow, RunThrow};
        match *self {
            BlastnOutputFormat::Tabular | BlastnOutputFormat::TabularWithComments => {
                FormatThrow::Writes(RunThrow::every_row(AfterThrow::WritesRow))
            }
            BlastnOutputFormat::CustomTabular(spec) => {
                spec.format_throw(crate::blastinput::app::Program::Blastn)
            }
            BlastnOutputFormat::Xml2(kind) => {
                FormatThrow::Writes(RunThrow::every_row(kind.after_throw()))
            }
            BlastnOutputFormat::Seqalign(_) | BlastnOutputFormat::Archive => FormatThrow::Unread,
            BlastnOutputFormat::Pairwise(_)
            | BlastnOutputFormat::Xml
            | BlastnOutputFormat::Sam(_) => FormatThrow::Refused,
        }
    }

    /// Whether NCBI prefetches the subjects that each query batch reports before the batch's
    /// reports and the report shows no taxonomy of them (outfmt 0, 5, 6, 7 and 17; with
    /// `PRE_FETCH_SEQS_LIMIT=2147483647` also 8, 9 and 12, `database::prefetches_subjects`): as
    /// the call's first report, the taxonomy database's messages of that prefetch go with the
    /// query warnings.
    pub(crate) fn prefetches_subjects(self) -> bool {
        crate::algorithm::common::database::prefetches_subjects(match self {
            BlastnOutputFormat::Pairwise(view) => Some(view.outfmt()),
            BlastnOutputFormat::Xml => Some(5),
            BlastnOutputFormat::Tabular => Some(6),
            BlastnOutputFormat::TabularWithComments => Some(7),
            // NCBI's `s_IsPrefetchFormat` prefetches for every tabular format (6, 7, 10 and 20):
            // they read the subjects as 6 and 7 do; the rows' own lookups of the taxonomy names
            // come with the rows (`report::tabular::TabularNames`).
            BlastnOutputFormat::CustomTabular(spec) => Some(match spec.kind {
                crate::report::tabular::TabularKind::Commented => 7,
                _ => 6,
            }),
            BlastnOutputFormat::Sam(_) => Some(17),
            BlastnOutputFormat::Seqalign(serial) => Some(serial.outfmt()),
            // The archive fetches nothing; XML2 and JSON write the messages at each `Hit`
            // (`TaxdbPrinter`).
            BlastnOutputFormat::Archive | BlastnOutputFormat::Xml2(_) => None,
        })
    }

    /// Whether the format is the pairwise report of outfmt 0 to 4 (`Pairwise`): the same prolog,
    /// query headers, description tables, footers and epilog, and the same reads of the subjects.
    pub(crate) fn is_pairwise(self) -> bool {
        matches!(self, BlastnOutputFormat::Pairwise(_))
    }

    /// Whether the format is one of the tabular formats (outfmt 6, 7, 10 and 20), which read the
    /// tabular subject IDs (`report::tabular::writes_rows`: not the organism report).
    pub(crate) fn is_tabular(self) -> bool {
        match self {
            BlastnOutputFormat::Tabular | BlastnOutputFormat::TabularWithComments => true,
            BlastnOutputFormat::CustomTabular(spec) => spec.writes_rows(),
            _ => false,
        }
    }

    /// Whether the format is a table whose fields read the sequence that NCBI's formatter finds
    /// by a subject's Seq-align id (`report::tabular::reads_subject`: the default fields do).
    pub(crate) fn tabular_reads_subject(self) -> bool {
        self.is_tabular()
            && crate::report::tabular::reads_subject(
                self.tabular_spec(),
                crate::blastinput::app::Program::Blastn,
            )
    }

    /// The field list of the custom tabular report (`CustomTabular`).
    pub(crate) fn tabular_spec(self) -> Option<&'static crate::report::tabular::TabularSpec> {
        match self {
            BlastnOutputFormat::CustomTabular(spec) => Some(spec),
            _ => None,
        }
    }

    /// Whether the report reads the queries' letters (EXPERIMENT, LOSAT_X_BNQPACK's rule): the
    /// HSP rows of outfmt 0 and of the XML / JSON documents, the archive's Bioseqs and DUST, and
    /// SAM's SEQ where the query is the target (`SQ SR`). The default tabular formats, the
    /// Seq-aligns and the other SAM documents show none; a custom tabular report shows them
    /// where its fields need the alignment rows (`TabularNeeds::rows`).
    pub(crate) fn reads_query_letters(self) -> bool {
        match self {
            BlastnOutputFormat::Tabular
            | BlastnOutputFormat::TabularWithComments
            | BlastnOutputFormat::Seqalign(_) => false,
            BlastnOutputFormat::CustomTabular(spec) => spec.needs().rows,
            BlastnOutputFormat::Sam(spec) => spec.sequence && spec.subject_as_reference,
            BlastnOutputFormat::Pairwise(_)
            | BlastnOutputFormat::Xml
            | BlastnOutputFormat::Xml2(_)
            | BlastnOutputFormat::Archive => true,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1264-1284
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:762
// A query of a batch that NCBI did not search has no alignment set (`num_hits` is
// `None`), so its header has no count (local_blast.cpp:177-207).
fn write_blastn_outfmt7_header<W: Write>(
    writer: &mut W,
    query_title: &[u8],
    subject_title: &str,
    num_hits: Option<usize>,
) -> io::Result<()> {
    writeln!(writer, "# BLASTN {}", crate::NCBI_BLAST_VERSION)?;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1305-1308
    // The title's bytes (`write_outfmt7_query_line`).
    crate::report::outfmt6::write_outfmt7_query_line(writer, query_title)?;
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

pub(crate) fn blastn_hsp_count(hit_list: Option<&BlastnHitList>) -> usize {
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
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
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
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
    pub num_ident: usize,
    pub query_frame: i32,
    pub query_length: usize,
    pub q_idx: u32,
    pub s_idx: u32,
    pub raw_score: i32,
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:96-101
    /// NCBI's `Int4` `query.offset` and `query.end` in the HSP's context. They lie outside the
    /// context (`internal_q_offset_0 < 0`, or `internal_q_end_0` past the context's length,
    /// `BlastnHsp::outside_query_context`) for an HSP from a hit that passes an end of its
    /// context, which LOSAT refuses where its report starts
    /// (`blast_engine::run::post_process_hits_and_write`).
    pub internal_q_offset_0: i32,
    pub internal_q_end_0: i32,
    pub internal_s_offset_0: usize,
    pub internal_s_end_0: usize,
    pub internal_query_context_offset: i32,
    pub gap_info: Option<Vec<GapEditOp>>,
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
    pub num_positives: usize,
}

// No NCBI counterpart: the size of the record of a final HSP (64-bit targets; `usize` fields).
#[cfg(target_pointer_width = "64")]
const _: () = assert!(std::mem::size_of::<BlastnHsp>() == 176);

impl BlastnHsp {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2390-2391
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1126-1129
    /// Whether the HSP lies outside its query context (`query.offset < 0`, which NCBI's debug
    /// build asserts against, or `query.end` past the context's length `query_length`), where
    /// NCBI's report coordinates (`Blast_HSPGetAdjustedOffsets`) leave the query: an HSP from a
    /// hit whose ungapped extension passed an end of its context (an ungapped X-drop of 2^30
    /// raw or more, remaining-gaps X1b-o).
    #[inline]
    pub fn outside_query_context(&self) -> bool {
        self.internal_q_offset_0 < 0 || self.internal_q_end_0 as i64 > self.query_length as i64
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1122-1136
/// NCBI's 1-based report offsets of a blastn HSP (`Blast_HSPGetAdjustedOffsets`) from its offsets
/// in its context: `query_offset` and `query_end` are NCBI's `Int4` ones (`BlastnHsp::
/// internal_q_offset_0`), `query_length` the context's length and `query_frame` its strand. The
/// query offsets are computed in `i64`; for an HSP outside its query context
/// (`BlastnHsp::outside_query_context`) one below 0 is kept as 0: such an HSP never reaches a
/// report (`run::post_process_hits_and_write` refuses the search).
pub(crate) fn adjust_blastn_offsets(
    query_offset: i32,
    query_end: i32,
    subject_offset: usize,
    subject_end: usize,
    query_length: usize,
    query_frame: i32,
) -> (usize, usize, usize, usize) {
    let (offset, end, length) = (
        i64::from(query_offset),
        i64::from(query_end),
        query_length as i64,
    );
    let reported = |value: i64| value.max(0) as usize;
    if query_frame < 0 {
        let q_end = length - offset;
        let q_start = q_end - end + offset + 1;
        (
            reported(q_start),
            reported(q_end),
            subject_end,
            subject_offset + 1,
        )
    } else {
        (
            reported(offset + 1),
            reported(end),
            subject_offset + 1,
            subject_end,
        )
    }
}

#[derive(Debug, Clone)]
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:153-166
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

#[derive(Debug, Clone)]
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:168-180
// `num_hits` (for mapping) is not stored: nothing reads it.
pub struct HitList<L> {
    pub hsplist_count: usize,
    pub hsplist_max: usize,
    pub worst_evalue: f64,
    pub low_score: i32,
    pub heapified: bool,
    pub hsplist_array: Vec<L>,
    pub hsplist_current: usize,
}

/// The hit list of a query's final HSP lists.
pub type BlastnHitList = HitList<BlastnHspList>;

/// No NCBI counterpart: reads LOSAT_X_BNHLGROW once. It is a strict-set switch: on unless it is
/// `0`. With it a hit list's array is not reserved to NCBI's array size (`kStartValue` = 100
/// lists, then doubled) when a list comes; it grows as lists come (`HitList::x_room_for_one`),
/// up to that size. NCBI's array holds pointers (800 bytes for 100), LOSAT's holds the lists
/// (40 bytes each, preliminary and final), so a search of many short queries with few subjects
/// each kept about 4 KB per query. It changes no value NCBI computes: the array's capacity is
/// read by nothing but `push`.
#[inline(always)]
fn x_bnhlgrow() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNHLGROW"))
}

pub type BlastnHspCompare = fn(&BlastnHsp, &BlastnHsp) -> Ordering;

impl BlastnHsp {
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:153-166
    pub fn into_hit(self) -> Hit {
        Hit {
            identity: self.identity,
            length: self.length,
            mismatch: self.mismatch as i64,
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
            // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:43-70 (GetPrelimHitlistSize)
/// The preliminary hit list size, computed without overflow. NCBI computes it in `Int4`
/// (c++/src/algo/blast/core/blast_hits.c:46-68), a signed overflow (undefined in C) from a hit
/// list size of 2^30 on with gapped search and of 2^30 - 25 on with composition-based
/// statistics, where NCBI BLAST+ 2.17.0's build gives 10, a size of 0 or less (NCBI crashes) or
/// 2 to 48. LOSAT does not reproduce a signed `Int4` overflow (PD-LOSAT-NCBI-DEFECTS, the
/// approved exception for hit list sizes past `Int4`): such a size keeps every hit, as NCBI's
/// 2^30 - 26 does. On a 32-bit target a size past `usize` is `usize::MAX`, which no hit list
/// reaches (NCBI's OIDs are `Int4`).
pub fn get_prelim_hitlist_size(
    hitlist_size: usize,
    composition_based_stats: bool,
    gapped_calculation: bool,
) -> usize {
    // The argument is an `int` (`CArg_Integer`, at most 2^31 - 1); the default is 500.
    let hitlist_size = hitlist_size as u64;
    let prelim_hitlist_size = if composition_based_stats {
        if crate::core::blast_hspstream::adaptive_cbs_environment() {
            if hitlist_size < 1000 {
                std::cmp::max(hitlist_size + 1000, 1500)
            } else {
                hitlist_size.saturating_mul(2).saturating_add(50)
            }
        } else if hitlist_size <= 500 {
            1050
        } else {
            hitlist_size.saturating_mul(2).saturating_add(50)
        }
    } else if gapped_calculation {
        std::cmp::min(
            std::cmp::max(hitlist_size.saturating_mul(2), 10),
            hitlist_size.saturating_add(50),
        )
    } else {
        hitlist_size
    };
    usize::try_from(prelim_hitlist_size).unwrap_or(usize::MAX)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1389-1403
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1353
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1122-1132
pub fn score_compare_hsps(a: &BlastnHsp, b: &BlastnHsp) -> Ordering {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1353
    // ScoreCompareHSPs compares canonical internal offsets/endpoints directly (NCBI's `Int4`
    // offsets, which may lie outside the context: `BlastnHsp::internal_q_offset_0`).
    let (a_q_offset, a_q_end) = (a.internal_q_offset_0, a.internal_q_end_0);
    let (b_q_offset, b_q_end) = (b.internal_q_offset_0, b.internal_q_end_0);
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1414-1435
pub fn evalue_compare_hsps(a: &BlastnHsp, b: &BlastnHsp) -> Ordering {
    let cmp = evalue_comp(a.e_value, b.e_value);
    if cmp == Ordering::Equal {
        score_compare_hsps(a, b)
    } else {
        cmp
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1381
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1453-1454
// Preserve NCBI's comparator and invocation timing while using one stable,
// target-neutral Rust sort. Comparator-equal HSPs retain their incoming order.
// NCBI sorts its array of HSP pointers; the HSPs are sorted through an order of
// their positions (`x_sort_stable_by_order`), so no scratch copy of the HSPs is made.
pub fn qsort_blastn_hsps_by(hsps: &mut [BlastnHsp], compare: BlastnHspCompare) {
    if hsps.len() <= 1 {
        return;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1381
    crate::utils::x_order::x_sort_stable_by_order(hsps, compare);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1382
pub fn sort_hsps_by_score(hsps: &mut [BlastnHsp]) {
    if hsps.len() <= 1 {
        return;
    }
    let sorted = hsps
        .windows(2)
        .all(|pair| score_compare_hsps(&pair[0], &pair[1]) != Ordering::Greater);
    if !sorted {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1381
        qsort_blastn_hsps_by(hsps, score_compare_hsps);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1437-1455
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
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1453-1454
            qsort_blastn_hsps_by(&mut list.hsps, evalue_compare_hsps);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1739-1748
pub fn update_best_evalue(list: &mut BlastnHspList) {
    let mut best = i32::MAX as f64;
    for hsp in &list.hsps {
        if hsp.e_value < best {
            best = hsp.e_value;
        }
    }
    list.best_evalue = best;
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2049-2067
// The filter of every program (`core::blast_hits::trim_by_max_hsps`).
pub fn trim_by_max_hsps(list: &mut BlastnHspList, max_hsps_per_subject: usize) {
    crate::core::blast_hits::trim_by_max_hsps(&mut list.hsps, max_hsps_per_subject);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3077-3106
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1626-1650
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1659-1676
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3125-3133
    pub fn new(hitlist_size: usize) -> Self {
        Self {
            hsplist_count: 0,
            hsplist_max: hitlist_size,
            worst_evalue: 0.0,
            low_score: i32::MAX,
            heapified: false,
            hsplist_array: Vec::new(),
            hsplist_current: 0,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3219-3240
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
        // EXPERIMENT (LOSAT_X_BNHLGROW): the array is not reserved to `hsplist_current` here;
        // `x_room_for_one` grows it as the lists come.
        if self.hsplist_array.capacity() < self.hsplist_current && !x_bnhlgrow() {
            self.hsplist_array
                .reserve(self.hsplist_current - self.hsplist_array.capacity());
        }
        true
    }

    /// No NCBI counterpart (LOSAT_X_BNHLGROW): room for one more list in a full array below
    /// `hsplist_current`, NCBI's array size: doubled from 4, never past that size. Without the
    /// switch `grow_hsplist_array` has reserved `hsplist_current` lists, so no such array is full.
    /// An array filled by its owner without `update` (its `hsplist_current` stays 0) is left to
    /// `push`, as before.
    #[inline]
    fn x_room_for_one(&mut self) {
        let len = self.hsplist_array.len();
        if len == self.hsplist_array.capacity() && len < self.hsplist_current {
            let room = (2 * len).max(4).min(self.hsplist_current);
            self.hsplist_array.reserve_exact(room - len);
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3195-3209
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3243-3298
    pub fn update(&mut self, mut hsp_list: L) {
        hsp_list.update_best_evalue();

        if self.hsplist_count < self.hsplist_max {
            if self.hsplist_current == self.hsplist_count && !self.grow_hsplist_array() {
                return;
            }
            self.x_room_for_one();
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3172-3187
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3328-3337
    pub fn sort_by_evalue(&mut self) {
        if self.hsplist_count > 1 {
            self.hsplist_array.sort_by(compare_hsp_lists);
        }
        self.purge();
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:877-892
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1700-1707, 1768-1772
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:591-608 (the traceback
    // reads the preliminary lists by OID, `BlastHSPStreamBatchRead`) and
    // c++/src/algo/blast/core/blast_traceback_mt_priv.c:173-188 (one thread's hit list is moved
    // to the results as it is).
    /// The hit list as NCBI's one-thread traceback leaves it where it does not sort the results
    /// at the end (`DbStatistics::sorts_results`: a database length of 0, or of 2^63 or more):
    /// the lists that kept HSPs, inserted in the order of their OIDs into a hit list of
    /// `hitlist_size` (`Blast_HitListUpdate`: the order of insertion, then the heap's order once
    /// a list comes to a full hit list). The output of every thread count is this order (Owner
    /// 2026-10-13: NCBI's `-num_threads 1` bytes; NCBI's own order differs with threads, probe
    /// p5_L0#41 of line RG D1a-3).
    pub fn put_in_traceback_order(&mut self, hitlist_size: usize) {
        let mut lists = std::mem::take(&mut self.hsplist_array);
        lists.retain(|list| list.hsp_count() > 0);
        // One list per subject of a query: the OIDs are distinct.
        lists.sort_by_key(|list| list.oid());
        let mut rebuilt = Self::new(hitlist_size);
        for list in lists {
            rebuilt.update(list);
        }
        *self = rebuilt;
    }

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

/// Write BLASTN hit lists to an existing writer. `epilog` is `Some(n)` when NCBI's epilog
/// (`PrintEpilog`) follows these reports, `n` being the number of queries that the search's
/// reports formatted (outfmt 7's last line); `None` without it (an error in a later query batch).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1100-1108
#[allow(clippy::too_many_arguments)]
pub fn write_output_blastn_hitlists_to_writer<W: crate::report::ostream::ReportWrite>(
    hit_lists: &[Option<BlastnHitList>],
    writer: &mut W,
    query_ids: &[Arc<[u8]>],
    subject_ids: &[crate::algorithm::common::tabular::SubjectId],
    loads: &crate::report::db_subject::HitLoads,
    output_format: BlastnOutputFormat,
    query_titles: &[Arc<[u8]>],
    subject_title: &str,
    unsearched: &[bool],
    epilog: Option<usize>,
    mut probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
) -> io::Result<()> {
    let config = OutputConfig::ncbi_compat();
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411
    // The rows are printed in the order of the final HSP list, so a running count is
    // each row's HSP index.
    let mut hsp_index: HspIndex = 0;

    for (q_idx, hit_list_opt) in hit_lists.iter().enumerate() {
        // The query's warnings come before its lines (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        // NCBI loads each subject that the query's batch did not prefetch at its first row
        // (`QueryWarnings::subject_read`).
        let mut read = |s_idx: u32, writer: &mut W| match warnings.as_deref_mut() {
            Some(warnings) => warnings.subject_read(q_idx, s_idx, writer),
            None => Ok(()),
        };
        write_query_rows(
            writer,
            q_idx,
            hit_list_opt.as_ref(),
            query_ids,
            subject_ids,
            loads,
            output_format,
            query_titles,
            subject_title,
            unsearched,
            &config,
            &mut hsp_index,
            probe.as_deref_mut(),
            &mut read,
        )?;
        // `~CBlastTabularInfo` of the query's report flushes the stream.
        writer.ncbi_flush()?;
    }
    write_tabular_epilog(writer, output_format, epilog)
}

/// One query's lines of `write_output_blastn_hitlists_to_writer` (outfmt 7's header, then one row
/// per HSP of its final hit list), without the warnings that come before them. `hsp_index` counts
/// the rows written so far (the format observer's HSP index). The EXPERIMENT LOSAT_X_BNTABPAR
/// (`x_tab_par.rs`) formats pieces of queries into buffers with `write_query_header` and
/// `write_hsp_rows`, the two parts of this function, on the search's pool.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:828-832
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_query_rows<W: Write>(
    writer: &mut W,
    q_idx: usize,
    hit_list_opt: Option<&BlastnHitList>,
    query_ids: &[Arc<[u8]>],
    subject_ids: &[crate::algorithm::common::tabular::SubjectId],
    loads: &crate::report::db_subject::HitLoads,
    output_format: BlastnOutputFormat,
    query_titles: &[Arc<[u8]>],
    subject_title: &str,
    unsearched: &[bool],
    config: &OutputConfig,
    hsp_index: &mut HspIndex,
    mut probe: Option<&mut FormatProbe<'_>>,
    read: crate::report::query_warnings::SubjectReadHook<'_, W>,
) -> io::Result<()> {
    write_query_header(
        writer,
        q_idx,
        hit_list_opt,
        output_format,
        query_titles,
        subject_title,
        unsearched,
    )?;
    let hit_list = match hit_list_opt {
        Some(value) => value,
        None => return Ok(()),
    };
    // The query's formatter is new: its rows keep no subject before the first one that loads
    // it (`TabularRows`).
    let mut rows = crate::report::tabular::TabularRows::default();
    for hsp_list in &hit_list.hsplist_array {
        let Some(first) = hsp_list.hsps.first() else {
            continue;
        };
        // The subject's load throws where the query's batch has its prefetch off
        // (`report::db_subject::HitLoads`): each row is written as the catch leaves it
        // (`TabularRows::write_thrown_default`).
        if loads.thrown(first.s_idx as usize).is_some() {
            for hsp in &hsp_list.hsps {
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.begin(*hsp_index);
                }
                let query_id = query_ids
                    .get(hsp.q_idx as usize)
                    .map_or(&b"unknown"[..], |id| id.as_ref());
                rows.write_thrown_default(
                    writer,
                    query_id,
                    |kept| crate::algorithm::common::tabular::printed_id(subject_ids, kept),
                    hsp.e_value,
                    hsp.bit_score,
                    config,
                )?;
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.end(*hsp_index);
                }
                *hsp_index += 1;
            }
            continue;
        }
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-838
        // Every row of BLASTN's tables reads its subject: the first loads it.
        read(hsp_list.oid, writer)?;
        rows.found(first.s_idx as usize);
        write_hsp_rows(
            writer,
            &hsp_list.hsps,
            query_ids,
            subject_ids,
            config,
            hsp_index,
            probe.as_deref_mut(),
        )?;
    }
    Ok(())
}

/// outfmt 7's header lines of query `q_idx` (nothing for outfmt 6), the first lines of
/// `write_query_rows`.
pub(crate) fn write_query_header<W: Write>(
    writer: &mut W,
    q_idx: usize,
    hit_list_opt: Option<&BlastnHitList>,
    output_format: BlastnOutputFormat,
    query_titles: &[Arc<[u8]>],
    subject_title: &str,
    unsearched: &[bool],
) -> io::Result<()> {
    if output_format == BlastnOutputFormat::TabularWithComments {
        let query_title = query_titles
            .get(q_idx)
            .map(|title| title.as_ref())
            .unwrap_or(b"unknown");
        write_blastn_outfmt7_header(
            writer,
            query_title,
            subject_title,
            (!unsearched.get(q_idx).copied().unwrap_or(false))
                .then(|| blastn_hsp_count(hit_list_opt)),
        )?;
    }
    Ok(())
}

/// One row per HSP of `hsps` (`write_query_rows`' rows of one HSP list, or a run of them).
/// `hsp_index` counts the rows written so far (the format observer's HSP index).
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_hsp_rows<W: Write>(
    writer: &mut W,
    hsps: &[BlastnHsp],
    query_ids: &[Arc<[u8]>],
    subject_ids: &[crate::algorithm::common::tabular::SubjectId],
    config: &OutputConfig,
    hsp_index: &mut HspIndex,
    mut probe: Option<&mut FormatProbe<'_>>,
) -> io::Result<()> {
    for hsp in hsps {
        let query_id = query_ids
            .get(hsp.q_idx as usize)
            .map(|id| id.as_ref())
            .unwrap_or(b"unknown");
        // A failed lookup's rows are refused before (`check_failed_lookups`).
        let subject_id = subject_ids
            .get(hsp.s_idx as usize)
            .map(|id| id.as_deref().unwrap_or_default())
            .unwrap_or(b"unknown");
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1100-1108
        if blastn_trace::should_trace_range(
            "hitlist",
            hsp.q_idx * 2 + if hsp.query_frame < 0 { 1 } else { 0 },
            hsp.s_idx as usize,
            &String::from_utf8_lossy(subject_id),
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
                    String::from_utf8_lossy(query_id),
                    String::from_utf8_lossy(subject_id),
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
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1100-1108
        // One printed row is one HSP; the probe marks it without changing it.
        if let Some(probe) = probe.as_mut() {
            writer.flush()?;
            probe.begin(*hsp_index);
        }
        write_hit_fields_bytes(
            writer,
            query_id,
            subject_id,
            hsp.identity,
            hsp.num_ident,
            hsp.length,
            hsp.mismatch as i64,
            hsp.gapopen,
            hsp.q_start,
            hsp.q_end,
            hsp.s_start,
            hsp.s_end,
            hsp.e_value,
            hsp.bit_score,
            config,
        )?;
        if let Some(probe) = probe.as_mut() {
            writer.flush()?;
            probe.end(*hsp_index);
        }
        *hsp_index += 1;
    }
    Ok(())
}

/// The line after the last query's lines (outfmt 7's count of queries, when `epilog`).
pub(crate) fn write_tabular_epilog<W: crate::report::ostream::ReportWrite + ?Sized>(
    writer: &mut W,
    output_format: BlastnOutputFormat,
    epilog: Option<usize>,
) -> io::Result<()> {
    // `PrintEpilog` writes the last line; an error in a later query batch skips it.
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1429-1430
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2233-2235
    // The count is that of every query the search's reports formatted (`epilog`).
    if let (BlastnOutputFormat::TabularWithComments, Some(num_queries)) = (output_format, epilog) {
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1324
        writeln!(writer, "# BLAST processed {} queries", num_queries)?;
        // The epilog's `CBlastTabularInfo` flushes the stream in its destructor.
        writer.ncbi_flush()?;
    }

    Ok(())
}

#[cfg(test)]
#[path = "hsp_tests.rs"]
mod tests;
