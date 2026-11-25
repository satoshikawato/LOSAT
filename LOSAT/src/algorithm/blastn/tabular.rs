//! BLASTN's tabular reports with a field list or a delimiter (outfmt 6 and 7 other than the
//! default `std` fields with a tab) and outfmt 10 and 20: the adapter of the shared writer
//! (`report::tabular`, NCBI's `CBlastTabularInfo`). For each query in order it writes the query's
//! warnings, outfmt 7's comment lines and one row per HSP of the query's final hit list, in the
//! list's order; outfmt 20's header line comes before the first query of the run and outfmt 7's
//! count of queries after the last.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:811-834
//! What `SetFields` reads of an HSP is its Seq-align (`algorithm::blastn::blast_xml` has the
//! coordinates): the query on its plus strand, the subject descending for a minus-strand HSP
//! (frames `1/-1`), the score, e-value and bit score, the alignment's lengths and gaps, and for
//! `qseq`, `sseq`, `btop`, `positive` and `ppos` the two rows, from which it counts the
//! identities again (`report::tabular::TabularWriter::row`). BLASTN's formatter reads no other
//! sequence data (`SetNoFetch` of true) and has no scoring matrix: positives are identities.

use std::io::Write;
use std::sync::Arc;

use anyhow::{Context, Result};

use super::hsp::{blastn_hsp_count, BlastnHitList, BlastnHsp};
use super::pairwise::xml_rows_into;
use super::x_subject_store::SubjectLetters;
use crate::algorithm::common::tabular::TabularSubjects;
use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::common::GapEditOp;
use crate::report::query_warnings::QueryWarnings;
use crate::report::tabular::{
    coverage_percent, union_length, QueryCoveragePerHsp, TabularHsp, TabularKind, TabularNames,
    TabularQuery, TabularRows, TabularSpec, TabularWriter,
};
use crate::utils::lower_bound::lower_bound_by;

const PROGRAM: &str = "BLASTN";

/// What the report needs besides the final hit lists.
pub(crate) struct BlastnTabular<'a> {
    pub spec: &'a TabularSpec,
    /// The queries of this call.
    pub queries: &'a [FastaRecord],
    /// The records whose letters the rows read (`pairwise::QueryRecords`).
    pub query_records: super::pairwise::QueryRecords<'a>,
    /// By query: the shown id (the default outfmt 6 query id), the title (outfmt 7's
    /// `# Query:` line) and the searched length.
    pub query_ids: &'a [Arc<[u8]>],
    pub query_titles: &'a [Arc<[u8]>],
    pub searched_lengths: &'a [usize],
    /// The length of the `-query_loc` range as given (NCBI's `m_QueryRange`), `None` without it.
    pub range_length: Option<usize>,
    /// An ungapped search (`-ungapped`), whose Seq-aligns are Dense-diags (`qcovhsp`,
    /// `report::tabular::QueryCoveragePerHsp`).
    pub ungapped: bool,
    /// Outfmt 7's `# Database:` text.
    pub database: &'a str,
    pub subjects: TabularSubjects<'a>,
    /// The letters of the subjects, which the rows read.
    pub letters: SubjectLetters<'a>,
    /// By query: whether NCBI did not search its batch (outfmt 7's lines without a count).
    pub unsearched: &'a [bool],
    /// Whether outfmt 20's header line comes before the first query (the run's first report).
    pub prolog: bool,
    /// The number of queries formatted (outfmt 7's last line), `None` without the epilog.
    pub epilog: Option<usize>,
    /// Whether `taxdb.bti` is on NCBI's search path (`SeqDB_ResolveDbPath`).
    pub taxdb_found: bool,
    /// The names of the taxonomy database that the rows look up (`None`: every lookup fails).
    pub names: Option<&'a TabularNames>,
}

/// Writes outfmt 20's line of the field names (`PrintProlog`; nothing for the other kinds).
/// The taxonomy database's warning of the formatter that NCBI makes for it is written before the
/// search (`write_tabular_prolog_warning`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:348-357
pub(crate) fn write_tabular_prolog(spec: &TabularSpec, writer: &mut dyn Write) -> Result<()> {
    let version = format!("BLASTN {}", crate::wording::REPORT_VERSION);
    TabularWriter::new(spec, &version, None).prolog(writer)?;
    Ok(())
}

/// Writes the warning of the formatter that NCBI makes for outfmt 20's header line, which it
/// writes before it reads the first query batch (`PrintProlog`, blastn_app.cpp:258): the
/// taxonomy database's warning when the field list checks it and `taxdb.bti` is not found
/// (`x_CheckTaxDB`, `TabularSpec::checks_taxdb`). A run that finds it writes nothing here, so
/// the line itself can come with the first report, after LOSAT's refusals.
pub(crate) fn write_tabular_prolog_warning(
    spec: &TabularSpec,
    taxdb_found: bool,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    if let Some(line) = spec
        .formatter_warning(taxdb_found, "blastn")
        .filter(|_| spec.kind == TabularKind::CsvWithHeader)
    {
        diagnostics.write_all(line.as_bytes())?;
        diagnostics.flush()?;
    }
    Ok(())
}

/// Writes the report of every query of `hit_lists` (indexed by query) to `writer`, each query's
/// warnings (and the taxonomy database's warning of its formatter) before its lines, outfmt 20's
/// header line first when `tab.prolog` and outfmt 7's last line after the last query when
/// `tab.epilog`. `probe` marks each row (the format observer).
pub(crate) fn write_blastn_tabular(
    tab: &BlastnTabular<'_>,
    hit_lists: &[Option<BlastnHitList>],
    writer: &mut dyn Write,
    mut warnings: Option<&mut QueryWarnings<'_>>,
    mut probe: Option<&mut FormatProbe<'_>>,
) -> Result<()> {
    let spec = tab.spec;
    let needs = spec.needs();
    // Whether the rows read their subject (`TabularSpec::subject_fetch`).
    let fetches = spec.subject_fetch(crate::blastinput::app::Program::Blastn)
        != crate::report::tabular::SubjectFetch::Never;
    // The refusals of `sall*` (`ParsedSubject::all_tabular_ids`) come before the first byte;
    // those of the other ids are the reader's (`ShownSubjects`).
    for hit_list in hit_lists.iter().flatten().filter(|_| needs.all_deflines) {
        for hsp_list in &hit_list.hsplist_array {
            if let Some(first) = hsp_list.hsps.first() {
                tab.subjects.values(first.s_idx as usize, spec, PROGRAM)?;
            }
        }
    }
    let version = format!("BLASTN {}", crate::wording::REPORT_VERSION);
    let table = TabularWriter::new(spec, &version, None);
    // NCBI makes a formatter for each query's report and for outfmt 7's epilog, and each warns
    // when the field list asks for a name of the taxonomy database that it cannot find.
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:139-158
    let taxdb_warning = spec.formatter_warning(tab.taxdb_found, "blastn");
    let warn = |writer: &mut dyn Write, warnings: Option<&mut QueryWarnings<'_>>| -> Result<()> {
        if let (Some(line), Some(warnings)) = (taxdb_warning.as_deref(), warnings) {
            writer.flush()?;
            warnings.sink.write_all(line.as_bytes())?;
            warnings.sink.flush()?;
        }
        Ok(())
    };
    if tab.prolog {
        table.prolog(writer)?;
    }
    let mut buf: Vec<u8> = Vec::new();
    let (mut query_row, mut subject_row): (Vec<u8>, Vec<u8>) = (Vec::new(), Vec::new());
    let mut query_ranges: Vec<(u64, u64)> = Vec::new();
    let mut hsp_index: HspIndex = 0;
    for (q_idx, hit_list) in hit_lists.iter().enumerate() {
        // The query's warnings come before its lines (`QueryWarnings`), then the warning of its
        // formatter.
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        warn(&mut *writer, warnings.as_deref_mut())?;
        let hit_list = hit_list.as_ref();
        let searched_length = tab.searched_lengths.get(q_idx).copied().unwrap_or_default();
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:787-790
        let query_length = tab.range_length.unwrap_or(searched_length);
        let query = TabularQuery::of_ids(
            tab.queries.get(q_idx).and_then(FastaRecord::believed),
            tab.query_ids
                .get(q_idx)
                .map_or(&b"unknown"[..], |id| id.as_ref()),
            query_length,
        );
        buf.clear();
        // A query of a batch that NCBI did not search has no count (`write_query_header`).
        table.query_header(
            &mut buf,
            tab.query_titles
                .get(q_idx)
                .map_or(&b"unknown"[..], |title| title.as_ref()),
            tab.database.as_bytes(),
            (!tab.unsearched.get(q_idx).copied().unwrap_or(false))
                .then(|| blastn_hsp_count(hit_list)),
        )?;
        if probe.is_some() {
            writer.write_all(&buf)?;
            buf.clear();
        }
        let Some(hit_list) = hit_list else {
            writer.write_all(&buf)?;
            continue;
        };
        // The query's formatter (`TabularRows`).
        let mut rows = TabularRows::default();
        for hsp_list in &hit_list.hsplist_array {
            let Some(first) = hsp_list.hsps.first() else {
                continue;
            };
            let s_idx = first.s_idx as usize;
            let values = tab.subjects.values(s_idx, spec, PROGRAM)?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:595-605
            // The HSPs of one subject are one group of `InsertSubjectScores` (a hit list holds a
            // subject once).
            let qcovs = if needs.qcovs {
                query_ranges.clear();
                query_ranges.extend(
                    hsp_list
                        .hsps
                        .iter()
                        .map(|hsp| (hsp.q_start as u64 - 1, hsp.q_end as u64 - 1)),
                );
                coverage_percent(union_length(&mut query_ranges), query_length)
            } else {
                -1
            };
            let qcovus = if needs.qcovus {
                uniq_subject_coverage(&hsp_list.hsps, query_length)
            } else {
                -1
            };
            let subject = values.nucleotide_subject(tab.names, qcovs, qcovus);
            // A subject whose lookups fail shows what the last found subject set
            // (`TabularRows::last_found`), made again from its index.
            let kept_values = match rows.last_found() {
                Some(kept) if table.fails(&subject) => {
                    Some(tab.subjects.values(kept, spec, PROGRAM)?)
                }
                _ => None,
            };
            let kept = kept_values
                .as_ref()
                .map(|values| values.nucleotide_subject(tab.names, -1, -1));
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-838
            // A table whose rows read the subject (every row or none in BLASTN,
            // `TabularSpec::subject_fetch`) loads it at its first row where the query's batch
            // did not prefetch it: the texts of its taxonomy lookups follow the rows before it
            // (`QueryWarnings::subject_read`), before those of the row's own lookups. A table
            // that reads a subject whose load throws outside the subject block ends the run at
            // its first row (`TabularWriter::ends_run`).
            if fetches {
                if let Some(warnings) = warnings.as_deref_mut() {
                    if warnings.reads_with_texts(q_idx) {
                        writer.write_all(&buf)?;
                        buf.clear();
                        warnings.subject_read(q_idx, hsp_list.oid, &mut *writer)?;
                    }
                }
            }
            if let Some(end) = table.ends_run(
                &subject,
                crate::blastinput::app::Program::Blastn,
                first.num_ident,
            ) {
                return Err(crate::report::tabular::run_end(writer, &buf, end, || {
                    if let Some(warnings) = warnings.as_deref_mut() {
                        warnings.end_run();
                    }
                }));
            }
            // The subject's `qcovhsp` rule: its rows only divide.
            let qcovhsp = QueryCoveragePerHsp::of_subject(
                tab.ungapped,
                hsp_list.hsps.len(),
                searched_length as i32,
                query_length as i32,
            );
            for hsp in &hsp_list.hsps {
                // The texts of the row's lookups of the taxonomy database (`x_SetTaxInfoAll`,
                // `x_SetTaxInfo`) follow the rows before it.
                if let (Some(names), Some(warnings), true) = (
                    tab.names,
                    warnings.as_deref_mut(),
                    table.sets_subject(&subject),
                ) {
                    let text = names.row_messages(spec, subject.taxonomy);
                    if !text.is_empty() {
                        writer.write_all(&buf)?;
                        buf.clear();
                        writer.flush()?;
                        warnings.sink.write_all(&text)?;
                        warnings.sink.flush()?;
                    }
                }
                // A row whose subject's load throws prints none of its letters.
                if needs.rows && !subject.throws() {
                    xml_rows_into(
                        hsp,
                        tab.query_records,
                        tab.letters,
                        &mut query_row,
                        &mut subject_row,
                    )
                    .context("BLASTN tabular rows of a final HSP")?;
                }
                let span = (hsp.q_end - hsp.q_start + 1) as i32;
                let row = tabular_hsp(
                    hsp,
                    qcovhsp.of_hsp(span, span),
                    needs
                        .rows
                        .then_some((query_row.as_slice(), subject_row.as_slice())),
                );
                match probe.as_deref_mut() {
                    // One printed row is one HSP; the probe marks it without changing it.
                    Some(probe) => {
                        buf.clear();
                        table.row(&mut buf, &mut rows, &query, &subject, kept.as_ref(), &row);
                        writer.flush()?;
                        probe.begin(hsp_index);
                        writer.write_all(&buf)?;
                        writer.flush()?;
                        probe.end(hsp_index);
                        buf.clear();
                    }
                    None => table.row(&mut buf, &mut rows, &query, &subject, kept.as_ref(), &row),
                }
                hsp_index += 1;
            }
        }
        writer.write_all(&buf)?;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2233-2236
    if let (TabularKind::Commented, Some(num_queries)) = (spec.kind, tab.epilog) {
        warn(&mut *writer, warnings)?;
        table.epilog(writer, num_queries)?;
    }
    Ok(())
}

/// The values of one HSP's row: its Seq-align's coordinates, scores and lengths (`SetFields`),
/// its `qcovhsp` (`report::tabular::QueryCoveragePerHsp`), `rows` when a field shows them.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1061-1091
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1100-1121
/// The gaps are the letters of the edit script's gap operations; the gap openings and the
/// identities are the default row's (`write_hsp_rows`).
fn tabular_hsp<'r>(
    hsp: &BlastnHsp,
    qcovhsp: i32,
    rows: Option<(&'r [u8], &'r [u8])>,
) -> TabularHsp<'r> {
    let gaps = hsp.gap_info.as_deref().map_or(0, |ops| {
        ops.iter()
            .filter(|op| !matches!(op, GapEditOp::Sub(_)))
            .map(|op| op.num() as usize)
            .sum()
    });
    TabularHsp {
        q_start: hsp.q_start as i64,
        q_end: hsp.q_end as i64,
        s_start: hsp.s_start as i64,
        s_end: hsp.s_end as i64,
        evalue: hsp.e_value,
        bit_score: hsp.bit_score,
        score: hsp.raw_score,
        length: hsp.length,
        gap_opens: hsp.gapopen,
        gaps,
        num_ident: hsp.num_ident,
        q_frame: 1,
        s_frame: if hsp.s_start > hsp.s_end { -1 } else { 1 },
        rows,
        // NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:288-297
        // A BLASTN alignment (X-drop or greedy) starts and ends with a substitution, so its rows
        // have no end gaps: a `-` at a row's end fills a segment past the query's end (an HSP that
        // NCBI keeps outside its query context, `pairwise::context_letters`), the gap character.
        end_gaps: false,
        stored_gaps: [&[], &[]],
        qcovhsp,
    }
}

/// `qcovus` of one subject's HSPs (in the Seq-align order) for a query of `query_length`:
/// `GetUniqSeqCoverage` in percent as `InsertSubjectScores` keeps it, -1 (`N/A`) when nothing is
/// covered.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:646-657
fn uniq_subject_coverage(hsps: &[BlastnHsp], query_length: usize) -> i32 {
    if query_length == 0 {
        return -1;
    }
    // `uniq_coverage`, an `int`, of the `TSeqPos` length.
    let covered = uniq_seq_coverage(hsps) as i32;
    if covered == 0 {
        return -1;
    }
    let mut pct = 100.0 * f64::from(covered) / query_length as f64;
    if pct < 99.0 {
        pct += 0.5;
    }
    pct as i32
}

/// NCBI's `CRange<TSeqPos>`: `from` and `to_open` as `TSeqPos` (unsigned 32-bit, wrapping).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct SeqRange {
    from: u32,
    to_open: u32,
}

impl SeqRange {
    /// `CRange` from `from` to `to` (`Set` calls `SetFrom` with `from` and `SetTo` with `to`,
    /// which is `SetToOpen` with `to` plus 1).
    fn new(from: u32, to: u32) -> Self {
        SeqRange {
            from,
            to_open: to.wrapping_add(1),
        }
    }

    fn to(self) -> u32 {
        self.to_open.wrapping_sub(1)
    }

    fn is_empty(self) -> bool {
        self.to_open <= self.from
    }

    fn length(self) -> u32 {
        self.to_open.saturating_sub(self.from)
    }

    /// `CombineWith`.
    fn combine_with(&mut self, r: SeqRange) {
        if !r.is_empty() {
            if !self.is_empty() {
                self.from = self.from.min(r.from);
                self.to_open = self.to_open.max(r.to_open);
            } else {
                *self = r;
            }
        }
    }

    /// `s_FixMinusStrandRange`.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4829-4836
    fn fix_minus_strand(self) -> Self {
        if self.from > self.to() {
            SeqRange::new(self.to(), self.from)
        } else {
            self
        }
    }
}

/// NCBI's `CRangeCollection<TSeqPos>`: sorted ranges that do not overlap or touch (the parts of
/// `x_CombineWith`, `x_Subtract` and `GetCoveredLength` that `GetUniqSeqCoverage` calls, with
/// their `TSeqPos` arithmetic).
#[derive(Clone, Debug, Default)]
struct RangeCollection(Vec<SeqRange>);

impl RangeCollection {
    /// `lower_bound(begin, end, pos, p)` of `x_CombineWith` and `x_Subtract` (cited below) from
    /// position `start`: the first position of `start..len` whose range is not less than `pos` by
    /// `p` (less: `to_open <= pos`), or `len` when there is none.
    ///
    /// NCBI reference (598d8ae6): c++/include/util/range_coll.hpp:44-51
    /// The ranges are not always in the order of `to_open`: a query gap at the first or last subject
    /// letter of a part gives a query range whose `TSeqPos` arithmetic wrapped (`uniq_seq_coverage`).
    /// The position found then depends on which ranges are compared, so they are compared in the
    /// order NCBI's `lower_bound` compares them (`utils::lower_bound`).
    fn lower_bound(&self, start: usize, pos: u32) -> usize {
        lower_bound_by(start, self.0.len(), |i| self.0[i].to_open <= pos)
    }

    /// `x_CombineWith` of a `TRange`.
    ///
    /// NCBI reference (598d8ae6): c++/include/util/range_coll.hpp:321-346
    fn combine_with(&mut self, r: SeqRange) {
        let mut pos_from = r.from;
        let pos_to_open = r.to_open;
        if pos_from == 0 {
            pos_from += 1;
        }
        let begin = self.lower_bound(0, pos_from - 1);
        if begin < self.0.len() && self.0[begin].from <= pos_to_open {
            let mut end = self.lower_bound(begin, pos_to_open);
            self.0[begin].combine_with(r);
            if end < self.0.len() && self.0[end].from <= pos_to_open {
                self.0[begin].to_open = self.0[end].to_open;
                end += 1;
            }
            self.0.drain(begin + 1..end);
        } else {
            self.0.insert(begin, r);
        }
    }

    /// `x_Subtract` of a `TRange`.
    ///
    /// NCBI reference (598d8ae6): c++/include/util/range_coll.hpp:348-376
    fn subtract(&mut self, r: SeqRange) {
        let pos_from = r.from;
        let pos_to_open = r.to_open;
        let mut begin = self.lower_bound(0, pos_from);
        if begin == self.0.len() {
            return;
        }
        if self.0[begin].from < pos_from && self.0[begin].to_open > pos_to_open {
            let range = self.0[begin];
            self.0.insert(begin, range);
            self.0[begin].to_open = pos_from;
            self.0[begin + 1].from = pos_to_open;
        } else {
            if self.0[begin].from < pos_from {
                self.0[begin].to_open = pos_from;
                begin += 1;
            }
            let end = self.lower_bound(begin, pos_to_open);
            if end < self.0.len() && self.0[end].from < pos_to_open {
                self.0[end].from = pos_to_open;
            }
            self.0.drain(begin..end);
        }
    }

    /// `GetCoveredLength`: the sum of the lengths, in `TSeqPos`.
    fn covered_length(&self) -> u32 {
        self.0
            .iter()
            .fold(0u32, |total, range| total.wrapping_add(range.length()))
    }
}

/// The subject range of an HSP's Seq-align (`GetSeqRange` of row 1, record coordinates, 0-based).
fn subject_range(hsp: &BlastnHsp) -> SeqRange {
    let (low, high) = (hsp.s_start.min(hsp.s_end), hsp.s_start.max(hsp.s_end));
    SeqRange::new(low as u32 - 1, high as u32 - 1)
}

/// The query range of an HSP's Seq-align (`GetSeqRange` of row 0: the query always ascending).
fn query_range(hsp: &BlastnHsp) -> SeqRange {
    SeqRange::new(hsp.q_start as u32 - 1, hsp.q_end as u32 - 1)
}

/// The query position (record coordinates, 0-based) in the column of the HSP's alignment that
/// holds subject position `subject_pos`, -1 where the query has a gap there: `CAlnMap`'s
/// `GetSeqPosFromAlnPos` in row 0 of the column that `GetAlnPosFromSeqPos` gives for
/// `subject_pos` in row 1, on the HSP's Dense-seg, whose subject row ascends and whose query row
/// descends for a minus-strand HSP
/// (`s_CollectSeqAlignData`; the columns are the edit script's).
///
/// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnmap.cpp:663-686
fn query_position(hsp: &BlastnHsp, subject_pos: u32) -> i64 {
    let subject_first = hsp.s_start.min(hsp.s_end) as i64 - 1;
    let target = i64::from(subject_pos) - subject_first;
    let ungapped = [GapEditOp::Sub(
        (hsp.internal_q_end_0 - hsp.internal_q_offset_0) as u32,
    )];
    let ops = hsp.gap_info.as_deref().unwrap_or(&ungapped);
    // The query letters and subject letters before the current operation.
    let (mut query, mut subject) = (0i64, 0i64);
    for op in ops {
        let n = i64::from(op.num());
        match op {
            GapEditOp::Sub(_) => {
                if target < subject + n {
                    let offset = query + (target - subject);
                    return if hsp.query_frame < 0 {
                        hsp.q_end as i64 - 1 - offset
                    } else {
                        hsp.q_start as i64 - 1 + offset
                    };
                }
                query += n;
                subject += n;
            }
            GapEditOp::Del(_) => {
                if target < subject + n {
                    return -1;
                }
                subject += n;
            }
            GapEditOp::Ins(_) => query += n,
        }
    }
    -1
}

/// `GetUniqSeqCoverage` of one subject's HSPs in the Seq-align order: the query letters that
/// the HSPs cover, where a part of an HSP whose subject letters an HSP before it covered does
/// not count (its query positions mapped through the alignment's columns).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4838-4896
/// A query gap in the column of a part's first or last subject letter gives NCBI the position
/// -1, which the `TSeqPos` range keeps as 2^32 - 1 (the same arithmetic here). LOSAT's BLASTN
/// alignments are always gapped Dense-segs (no Dense-diag).
fn uniq_seq_coverage(hsps: &[BlastnHsp]) -> u32 {
    let Some((first, rest)) = hsps.split_first() else {
        return 0;
    };
    let mut subject_ranges = RangeCollection(vec![subject_range(first)]);
    let mut query_ranges = RangeCollection(vec![query_range(first).fix_minus_strand()]);
    for hsp in rest {
        let align_subject = subject_range(hsp);
        let mut parts = RangeCollection(vec![align_subject]);
        for &covered in &subject_ranges.0 {
            parts.subtract(covered);
        }
        if parts.0.is_empty() {
            continue;
        }
        if parts.0[0] == align_subject {
            query_ranges.combine_with(query_range(hsp).fix_minus_strand());
            subject_ranges.combine_with(align_subject);
        } else {
            for part in parts.0 {
                // `SetFrom` and `SetTo` of a `TSignedSeqPos` into the `TSeqPos` range.
                let from = query_position(hsp, part.from) as u32;
                let to = query_position(hsp, part.to()) as u32;
                let query_part = SeqRange {
                    from,
                    to_open: to.wrapping_add(1),
                };
                query_ranges.combine_with(query_part.fix_minus_strand());
                subject_ranges.combine_with(part);
            }
        }
    }
    query_ranges.covered_length()
}

#[cfg(test)]
#[path = "tabular_tests.rs"]
mod tests;
