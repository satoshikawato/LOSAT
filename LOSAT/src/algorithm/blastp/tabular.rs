//! BLASTP's (and psiblast's) tabular reports with a field list or a delimiter (outfmt 6 and 7
//! other than the default `std` fields with a tab) and outfmt 10 and 20: the adapter of the shared
//! writer (`report::tabular`, NCBI's `CBlastTabularInfo`). For each query in order it writes the
//! query's warnings, outfmt 7's comment lines (psiblast's with its round) and one row per HSP of
//! the query's final hit list, in the list's order; outfmt 20's header line comes before the
//! first query of the run and outfmt 7's count of queries after the last.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:811-834
//! ```c++
//!         if (results.HasAlignments()) {
//!     	    CSeq_align_set copy_aln_set;
//!             CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_HitlistSize);
//!
//!             {
//!             	unsigned int scores = CBlastFormatUtil::eNoQuerySubjCov;
//!             	if(string::npos != m_CustomOutputFormatSpec.find("qcovs"))
//!             		scores |= CBlastFormatUtil::eQueryCovPerSubj ;
//!             	if(string::npos != m_CustomOutputFormatSpec.find("qcovus") &&
//!             			ncbi::NStr::ToLower(m_Program) == string("blastn"))
//!             		scores |= CBlastFormatUtil::eQueryCovPerUniqSubj ;
//!             ...
//!             ITERATE(CSeq_align_set::Tdata, itr, copy_aln_set.Get()) {
//!                     const CSeq_align& s = **itr;
//!                     tabinfo.SetFields(s, *m_Scope, &m_ScoringMatrix);
//!                     tabinfo.Print();
//!             }
//! ```
//! What `SetFields` reads of a BLASTP HSP is its Seq-align (a Dense-seg of two proteins): the
//! coordinates, the score, e-value and bit score, the alignment's lengths and gaps, and for
//! `qseq`, `sseq`, `btop`, `positive` and `ppos` the two rows (the rendered alignment that outfmt
//! 0 shows, `PairwiseHit`), from which it counts the identities and, with the formatter's matrix
//! (`GetAsciiProteinMatrix`), the positives again (`report::tabular::TabularWriter::row`). The
//! frames are 1, `sstrand` is `N/A` (a protein subject) and `qcovus` is `N/A` (NCBI computes it
//! for blastn only).

use std::borrow::Cow;
use std::io::Write;
use std::sync::Arc;

use anyhow::Result;

use super::args::BlastpApp;
use super::hsp::{BlastpHitList, BlastpHsp};
use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::config::ScoringMatrix;
use crate::report::db_subject::ParsedSubject;
use crate::report::pairwise::PairwiseHit;
use crate::report::query_warnings::QueryWarnings;
use crate::report::tabular::{
    coverage_percent, hsp_coverage, taxdb_warning_line, union_length, DeflineIds, TabularField,
    TabularHsp, TabularKind, TabularNames, TabularQuery, TabularSpec, TabularSubject,
    TabularTaxonomy, TabularWriter,
};

/// What the report needs besides the final hit lists.
pub(super) struct BlastpTabular<'a> {
    pub spec: &'a TabularSpec,
    /// blastp or psiblast: NCBI's messages name it (`ncbi_name`), LOSAT's refusals too (`name`).
    pub app: BlastpApp,
    /// outfmt 7's `# BLASTP 2.17.0+` line (`# PSIBLAST 2.17.0+`).
    pub program_version: &'a str,
    /// psiblast's round (`# Iteration: 1`), `None` for blastp.
    pub iteration: Option<u32>,
    /// The search's matrix, whose formatting matrix counts the positives of the rows.
    pub matrix: ScoringMatrix,
    /// By query: the shown id (the default outfmt 6 query id), the title (outfmt 7's
    /// `# Query:` line) and the searched length.
    pub query_ids: &'a [Arc<[u8]>],
    pub query_titles: &'a [Arc<[u8]>],
    pub searched_lengths: &'a [usize],
    /// The length of the `-query_loc` range as given (NCBI's `m_QueryRange`), `None` without it.
    pub range_length: Option<usize>,
    /// Outfmt 7's `# Database:` text.
    pub database: &'a str,
    pub subjects: TabularSubjects<'a>,
    /// By query: whether NCBI did not search its batch (outfmt 7's lines without a count).
    pub unsearched: &'a [bool],
    /// psiblast's count of outfmt 7's `# N hits found` line by query (the HSPs of its larger
    /// traceback hit list), `None` for blastp (the HSPs of the final hit list).
    pub hits_found: Option<&'a [usize]>,
    /// The rendered alignments of the final HSPs in the order of the hit lists, when a field
    /// shows the rows (`TabularNeeds::rows`).
    pub rows: Option<&'a [PairwiseHit]>,
    /// Whether outfmt 20's header line comes before the first query (the run's first report).
    pub prolog: bool,
    /// The number of queries formatted (outfmt 7's last line), `None` without the epilog.
    pub epilog: Option<usize>,
    /// Whether `taxdb.bti` is on NCBI's search path (`SeqDB_ResolveDbPath`).
    pub taxdb_found: bool,
    /// The names of the taxonomy database that the rows look up (`None`: every lookup fails).
    pub names: Option<&'a TabularNames>,
}

/// The subjects as the report shows them, by subject index (a subject set) or rank (the
/// reported sequences of a database).
pub(super) enum TabularSubjects<'a> {
    /// The records of `-subject`: the tabular id of each (`sseqid`, `sacc`, `saccver`; the
    /// default report's) and the length of its record.
    Records {
        ids: &'a [Arc<[u8]>],
        length: &'a dyn Fn(usize) -> usize,
    },
    /// The reported sequences of a database: their tabular ids (`sseqid`, and `sacc` and
    /// `saccver` of a sequence with parsed Seq-ids; the default report's), ordinal ids, stored
    /// titles (the title of the first definition line of a sequence with parsed Seq-ids), parsed
    /// Seq-ids (empty when none has them), lengths, and taxonomy (empty when no field shows it).
    Database {
        ids: &'a [Arc<[u8]>],
        accessions: &'a [Option<[Arc<[u8]>; 2]>],
        ordinals: &'a [u32],
        titles: &'a [Arc<[u8]>],
        parsed: &'a [Option<Arc<ParsedSubject>>],
        length: &'a dyn Fn(usize) -> usize,
        taxonomy: &'a [TabularTaxonomy],
    },
}

/// The values of one subject that its rows show (`TabularSubject` without the coverages).
struct SubjectValues<'a> {
    id: &'a [u8],
    accession: &'a [u8],
    accession_version: &'a [u8],
    gi: i64,
    all_ids: Vec<DeflineIds>,
    titles: Option<Vec<&'a [u8]>>,
    length: usize,
    taxonomy: Cow<'a, TabularTaxonomy>,
}

impl<'a> TabularSubjects<'a> {
    /// The values of subject `s_idx` that the fields of `spec` show.
    fn values(
        &self,
        s_idx: usize,
        spec: &TabularSpec,
        app: BlastpApp,
    ) -> Result<SubjectValues<'a>> {
        let needs = spec.needs();
        let shows_length = spec.has(TabularField::SubjectLength);
        match *self {
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:577-582
            // ```c++
            //         list<CRef<objects::CSeq_id> > subject_id_list;
            //         ITERATE(CBioseq_Handle::TId, itr, bh.GetId()) {
            //             CRef<CSeq_id> next_id = s_ReplaceLocalId(bh, itr->GetSeqId(), m_ParseLocalIds);
            //             subject_id_list.push_back(next_id);
            //         }
            //         m_SubjectIds.push_back(subject_id_list);
            // ```
            // A FASTA subject has no definition lines: its `sall*` are its own ids, its title
            // and taxonomy are empty (`N/A`).
            TabularSubjects::Records { ids, length } => {
                let id: &[u8] = ids.get(s_idx).map_or(&b""[..], |id| id.as_ref());
                Ok(SubjectValues {
                    id,
                    accession: id,
                    accession_version: id,
                    gi: 0,
                    all_ids: Vec::new(),
                    titles: None,
                    length: if shows_length { length(s_idx) } else { 0 },
                    taxonomy: Cow::Owned(TabularTaxonomy::default()),
                })
            }
            TabularSubjects::Database {
                ids,
                accessions,
                ordinals,
                titles,
                parsed,
                length,
                taxonomy,
            } => {
                let id: &[u8] = ids.get(s_idx).map_or(&b""[..], |id| id.as_ref());
                let (accession, accession_version): (&[u8], &[u8]) =
                    match accessions.get(s_idx).and_then(Option::as_ref) {
                        Some([accession, version]) => (accession, version),
                        None => (id, id),
                    };
                let taxonomy = match taxonomy.get(s_idx) {
                    Some(taxonomy) => Cow::Borrowed(taxonomy),
                    None => Cow::Owned(TabularTaxonomy::default()),
                };
                let length = if shows_length { length(s_idx) } else { 0 };
                match parsed.get(s_idx).and_then(Option::as_deref) {
                    // A sequence with parsed Seq-ids: the ids of each definition line and their
                    // titles.
                    Some(subject) => {
                        let all_ids = if needs.all_deflines {
                            subject.all_tabular_ids(true).map_err(|problem| {
                                crate::algorithm::common::database::parsed_refusal(
                                    ordinals.get(s_idx).copied().unwrap_or_default(),
                                    &format!("has {problem}"),
                                    app.name(),
                                )
                            })?
                        } else {
                            Vec::new()
                        };
                        Ok(SubjectValues {
                            id,
                            accession,
                            accession_version,
                            gi: subject.tabular_gi(),
                            all_ids,
                            titles: Some(
                                subject
                                    .deflines()
                                    .iter()
                                    .map(|defline| defline.title.as_deref().unwrap_or_default())
                                    .collect(),
                            ),
                            length,
                            taxonomy,
                        })
                    }
                    // A sequence without parsed Seq-ids: its one definition line's title.
                    None => Ok(SubjectValues {
                        id,
                        accession,
                        accession_version,
                        gi: 0,
                        all_ids: Vec::new(),
                        titles: Some(vec![titles.get(s_idx).map_or(&[][..], |title| title)]),
                        length,
                        taxonomy,
                    }),
                }
            }
        }
    }
}

/// The spec of a tabular report with a field list or a delimiter, or of outfmt 10 or 20, of a
/// `-outfmt` value read as NCBI reads it (`app::parse_formatting_string`); `None` for the other
/// formats and for the default outfmt 6 and 7 (`TabularSpec::of`).
pub(super) fn custom_tabular_spec(outfmt: &str) -> Option<&'static TabularSpec> {
    let choice = crate::blastinput::app::parse_formatting_string(outfmt).ok()?;
    let format = crate::blastinput::app::ReportFormat::from_number(choice.number)?;
    TabularSpec::of(format, &choice).map(TabularSpec::intern)
}

/// Writes the warning of the formatter that NCBI makes for outfmt 20's header line, which it
/// writes before it reads the first query batch (`PrintProlog`, blastp_app.cpp:253): the
/// taxonomy database's warning when the field list checks it and `taxdb.bti` is not found
/// (`x_CheckTaxDB`, `TabularSpec::checks_taxdb`). The line itself comes with the first report
/// (`write_tabular_prolog`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:348-357
/// ```c++
///     if (m_FormatType == CFormattingArgs::eCommaSeparatedValuesWithHeader) {
///         const CBlastTabularInfo::EFieldDelimiter kDelim = CBlastTabularInfo::eComma;
///         CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec, kDelim);
///         if(!m_CustomDelim.empty()) {
///             tabinfo.SetCustomDelim(m_CustomDelim);
///         }
///         tabinfo.PrintFieldSpecs();
///         return;
///     }
/// ```
pub(super) fn write_tabular_prolog_warning(
    spec: &TabularSpec,
    taxdb_found: bool,
    program: &str,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    if spec.kind == TabularKind::CsvWithHeader && spec.checks_taxdb() && !taxdb_found {
        diagnostics.write_all(taxdb_warning_line(program).as_bytes())?;
        diagnostics.flush()?;
    }
    Ok(())
}

/// Writes outfmt 20's line of the field names (`PrintProlog`; nothing for the other kinds).
pub(super) fn write_tabular_prolog(spec: &TabularSpec, writer: &mut dyn Write) -> Result<()> {
    TabularWriter {
        spec,
        program_version: "",
        matrix: None,
    }
    .prolog(writer)?;
    Ok(())
}

/// The reports of a run without a query batch: outfmt 20's line of the field names
/// (`PrintProlog`), outfmt 7's count of 0 queries with the warning of the epilog's formatter
/// (`PrintEpilog`), nothing for outfmt 6 and 10.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2233-2238
/// ```c++
///     if (m_FormatType == CFormattingArgs::eTabularWithComments) {
///         CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec);
///         tabinfo.PrintNumProcessed(m_QueriesFormatted);
///         return;
///     } else if (m_FormatType >= CFormattingArgs::eTabular)
///         return;  // No footer for these.
/// ```
pub(super) fn write_tabular_without_queries(
    spec: &TabularSpec,
    taxdb_found: bool,
    program: &str,
    writer: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    write_tabular_prolog(spec, writer)?;
    if spec.kind == TabularKind::Commented {
        if spec.checks_taxdb() && !taxdb_found {
            writer.flush()?;
            diagnostics.write_all(taxdb_warning_line(program).as_bytes())?;
            diagnostics.flush()?;
        }
        TabularWriter {
            spec,
            program_version: "",
            matrix: None,
        }
        .epilog(writer, 0)?;
    }
    Ok(())
}

/// Writes the report of every query of `hit_lists` (indexed by query) to `writer`, each query's
/// warnings (and the taxonomy database's warning of its formatter) before its lines, outfmt 20's
/// header line first when `tab.prolog` and outfmt 7's last line after the last query when
/// `tab.epilog`. `probe` marks each row (the format observer).
pub(super) fn write_blastp_tabular(
    tab: &BlastpTabular<'_>,
    hit_lists: &[Option<BlastpHitList>],
    writer: &mut dyn Write,
    mut warnings: Option<&mut QueryWarnings<'_>>,
    mut probe: Option<&mut FormatProbe<'_>>,
) -> Result<()> {
    let spec = tab.spec;
    let needs = spec.needs();
    // The refusals of `sall*` (`ParsedSubject::all_tabular_ids`) come before the first byte;
    // those of the other ids are made with them (`tabular_subject_ids`).
    if needs.all_deflines {
        for hit_list in hit_lists.iter().flatten() {
            for hsp_list in reported(hit_list) {
                if let Some(first) = hsp_list.hsps.first() {
                    tab.subjects.values(first.s_idx as usize, spec, tab.app)?;
                }
            }
        }
    }
    let table = TabularWriter {
        spec,
        program_version: tab.program_version,
        matrix: Some(tab.matrix),
    };
    // NCBI makes a formatter for each query's report and for outfmt 7's epilog, and each warns
    // when the field list asks for a name of the taxonomy database that it cannot find.
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:139-158
    // ```c++
    // CBlastTabularInfo::CBlastTabularInfo(CNcbiOstream& ostr, const string& format,
    // ...
    //     x_CheckTaxDB();
    // }
    // ```
    let taxdb_warning =
        (spec.checks_taxdb() && !tab.taxdb_found).then(|| taxdb_warning_line(tab.app.ncbi_name()));
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
    let mut query_ranges: Vec<(u64, u64)> = Vec::new();
    // The index of the next HSP in the final hit lists (the order of `tab.rows`).
    let mut hsp_index: HspIndex = 0;
    for (q_idx, hit_list) in hit_lists.iter().enumerate() {
        // The query's warnings come before its lines (`QueryWarnings`), then the warning of its
        // formatter.
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        warn(&mut *writer, warnings.as_deref_mut())?;
        let searched_length = tab.searched_lengths.get(q_idx).copied().unwrap_or_default();
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:787-790
        // ```c++
        //             if(m_QueryRange.NotEmpty())
        //             	m_QueryLength = m_QueryRange.GetLength();
        //             else
        //             	m_QueryLength = query_bh.GetBioseqLength();
        // ```
        let query_length = tab.range_length.unwrap_or(searched_length);
        let query = TabularQuery {
            id: tab.query_ids.get(q_idx).map_or(&b""[..], |id| id.as_ref()),
            length: query_length,
        };
        buf.clear();
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1277-1278
        // ```c++
        //     if (align_set) {
        //        int num_hits = align_set->Get().size();
        // ```
        // A query of a batch that NCBI did not search has no Seq-align set, so no count.
        let num_hits = (!tab.unsearched.get(q_idx).copied().unwrap_or(false)).then(|| {
            tab.hits_found.map_or_else(
                || {
                    hit_list.as_ref().map_or(0, |hit_list| {
                        reported(hit_list).map(|hsp_list| hsp_list.hsps.len()).sum()
                    })
                },
                |found| found.get(q_idx).copied().unwrap_or_default(),
            )
        });
        table.round_query_header(
            &mut buf,
            tab.iteration,
            tab.query_titles
                .get(q_idx)
                .map_or(&b""[..], |title| title.as_ref()),
            tab.database.as_bytes(),
            num_hits,
        )?;
        if probe.is_some() {
            writer.write_all(&buf)?;
            buf.clear();
        }
        let Some(hit_list) = hit_list else {
            writer.write_all(&buf)?;
            continue;
        };
        for hsp_list in reported(hit_list) {
            let Some(first) = hsp_list.hsps.first() else {
                continue;
            };
            let values = tab.subjects.values(first.s_idx as usize, spec, tab.app)?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:595-605
            // ```c++
            //     CConstRef<CBioseq> query_bioseq = query_handle.GetCompleteBioseq();
            //     int query_len = 0;
            //     if(query_range.NotEmpty())
            //     {
            //     	query_len = query_range.GetLength();
            //     }
            //     else if (!query_bioseq.Empty() && query_bioseq->IsSetLength())
            //     {
            //     	query_len = query_bioseq->GetLength();
            //     }
            // ```
            // The HSPs of one subject are one group of `InsertSubjectScores` (a hit list holds a
            // subject once); `qcovus` is blastn's only (blast_format.cpp:818-820).
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
            let subject = TabularSubject {
                id: values.id,
                accession: values.accession,
                accession_version: values.accession_version,
                gi: values.gi,
                all_ids: &values.all_ids,
                length: values.length,
                titles: values.titles.as_deref(),
                taxonomy: &values.taxonomy,
                names: tab.names,
                nucleotide: false,
                qcovs,
                qcovus: -1,
            };
            for hsp in &hsp_list.hsps {
                // The texts of the row's lookups of the taxonomy database (`x_SetTaxInfoAll`,
                // `x_SetTaxInfo`) follow the rows before it.
                if let (Some(names), Some(warnings)) = (tab.names, warnings.as_deref_mut()) {
                    let text = names.row_messages(spec, &values.taxonomy);
                    if !text.is_empty() {
                        writer.write_all(&buf)?;
                        buf.clear();
                        writer.flush()?;
                        warnings.sink.write_all(&text)?;
                        warnings.sink.flush()?;
                    }
                }
                let rows = match tab.rows.filter(|_| needs.rows) {
                    Some(rows) => Some(hsp_rows(rows, hsp_index, hsp)?),
                    None => None,
                };
                let row = tabular_hsp(hsp, searched_length, rows);
                match probe.as_deref_mut() {
                    // One printed row is one HSP; the probe marks it without changing it.
                    Some(probe) => {
                        buf.clear();
                        table.row(&mut buf, &query, &subject, &row);
                        writer.flush()?;
                        probe.begin(hsp_index);
                        writer.write_all(&buf)?;
                        writer.flush()?;
                        probe.end(hsp_index);
                        buf.clear();
                    }
                    None => table.row(&mut buf, &query, &subject, &row),
                }
                hsp_index += 1;
            }
        }
        writer.write_all(&buf)?;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2233-2236
    // ```c++
    //     if (m_FormatType == CFormattingArgs::eTabularWithComments) {
    //         CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec);
    //         tabinfo.PrintNumProcessed(m_QueriesFormatted);
    //         return;
    // ```
    if let (TabularKind::Commented, Some(num_queries)) = (spec.kind, tab.epilog) {
        warn(&mut *writer, warnings.as_deref_mut())?;
        table.epilog(writer, num_queries)?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:271-327
// ```c
// *hsp_list_out = hit_list->hsplist_array[last_hsplist_index];
// (*hsp_list_out)->query_index = index;
// --hit_list->hsplist_count;
// ```
/// The HSP lists of a final hit list that the reports show (its first `hsplist_count`), as the
/// default report walks them (`write_blastp_tabular_hit_lists`).
fn reported(hit_list: &BlastpHitList) -> impl Iterator<Item = &super::hsp::BlastpHspList> {
    hit_list.hsplist_array.iter().take(hit_list.hsplist_count)
}

/// The rows of the HSP `hsp` (the `index`-th of the final hit lists): its rendered alignment.
fn hsp_rows<'r>(
    rows: &'r [PairwiseHit],
    index: HspIndex,
    hsp: &BlastpHsp,
) -> Result<(&'r [u8], &'r [u8])> {
    let rendered = rows
        .get(index)
        .filter(|row| row.hit.q_idx == hsp.q_idx && row.hit.s_idx == hsp.s_idx);
    match rendered.and_then(|row| Some((row.query_seq.as_deref()?, row.subject_seq.as_deref()?))) {
        Some((query, subject)) => Ok((query.as_bytes(), subject.as_bytes())),
        None => anyhow::bail!("unimplemented BLASTP tabular rows of a final HSP"),
    }
}

/// The values of one HSP's row: its Seq-align's coordinates, scores and lengths (`SetFields`),
/// `rows` when a field shows them.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1071-1094
/// ```c++
///     int query_frame = 1, subject_frame = 1;
///     if (kTranslated) {
///     ...
///     else {
///     	if (x_IsFieldRequested(eSubjFrame) || x_IsFieldRequested(eFrames)) {
///     		if ((s_start - s_end) > 0 ) {
///     			subject_frame = -1;
///     		}
///     	}
///     }
///     SetCounts(num_ident, align_length, num_gaps, num_gap_opens, num_positives,
///               query_frame, subject_frame);
/// ```
/// NCBI reference (598d8ae6): c++/include/objtools/align_format/tabular.hpp:543-547
/// ```c++
///     int num_mismatches = m_AlignLength - m_NumIdent - m_NumGaps;
/// ```
/// The gaps are the alignment's columns that are neither identities nor mismatches, as the
/// default report's mismatches are the columns of neither (`Hit::gap_letters`); the gap
/// openings and the identities are the default row's. A protein HSP runs forward on both
/// sequences: the frames are 1.
fn tabular_hsp<'r>(
    hsp: &BlastpHsp,
    searched_length: usize,
    rows: Option<(&'r [u8], &'r [u8])>,
) -> TabularHsp<'r> {
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
        gaps: hsp
            .length
            .saturating_sub(hsp.num_ident.saturating_add(hsp.mismatch)),
        num_ident: hsp.num_ident,
        q_frame: 1,
        s_frame: 1,
        rows,
        // `Blast_HSPGetQueryCoverage` of the HSP's context, which has the searched length.
        qcovhsp: hsp_coverage((hsp.q_end - hsp.q_start + 1) as i32, searched_length as i32),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blastinput::app::{FormatChoice, ReportFormat};
    use crate::common::Hit;

    fn hit(q: (usize, usize), s: (usize, usize), s_idx: u32) -> Hit {
        Hit {
            identity: 50.0,
            length: q.1 - q.0 + 1,
            mismatch: 2,
            gapopen: 0,
            q_start: q.0,
            q_end: q.1,
            s_start: s.0,
            s_end: s.1,
            e_value: 1.0e-10,
            bit_score: 40.5,
            num_ident: q.1 - q.0 - 1,
            query_frame: 0,
            query_length: 100,
            q_idx: 0,
            s_idx,
            raw_score: 90,
            sort_query_offset: 0,
            sort_query_end: 0,
            sort_subject_offset: 0,
            sort_subject_end: 0,
            has_sort_offsets: false,
            gap_info: None,
            num_positives: 0,
        }
    }

    fn hit_list(hsp_lists: Vec<Vec<Hit>>, count: usize) -> BlastpHitList {
        BlastpHitList {
            hsplist_count: count,
            hsplist_max: hsp_lists.len(),
            worst_evalue: 0.0,
            low_score: 0,
            heapified: false,
            hsplist_array: hsp_lists
                .into_iter()
                .map(|hits| super::super::hsp::BlastpHspList {
                    oid: hits[0].s_idx,
                    query_index: 0,
                    hsps: hits.into_iter().map(BlastpHsp::from_hit).collect(),
                    best_evalue: 0.0,
                })
                .collect(),
            hsplist_current: 0,
            num_hits: 0,
        }
    }

    fn report(
        number: i32,
        fields: &str,
        app: BlastpApp,
        hit_lists: &[Option<BlastpHitList>],
    ) -> (String, String) {
        let format = ReportFormat::from_number(number).unwrap();
        let spec = TabularSpec::of(
            format,
            &FormatChoice {
                number,
                spec: fields.to_string(),
                delimiter: String::new(),
            },
        )
        .unwrap();
        let ids: Vec<Arc<[u8]>> = vec![Arc::from(&b"s1"[..]), Arc::from(&b"s2"[..])];
        let length = |s_idx: usize| 200 + s_idx;
        let tab = BlastpTabular {
            spec: &spec,
            app,
            program_version: "PSIBLAST 2.17.0+",
            iteration: (app == BlastpApp::Psiblast).then_some(1),
            matrix: ScoringMatrix::Blosum62,
            query_ids: &[Arc::from(&b"q1"[..]), Arc::from(&b"q2"[..])],
            query_titles: &[Arc::from(&b"q1 first"[..]), Arc::from(&b"q2"[..])],
            searched_lengths: &[100, 50],
            range_length: None,
            database: "User specified sequence set (Input: s.faa)",
            subjects: TabularSubjects::Records {
                ids: &ids,
                length: &length,
            },
            unsearched: &[false, true],
            hits_found: None,
            rows: None,
            prolog: true,
            epilog: Some(2),
            taxdb_found: false,
            names: None,
        };
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let mut warnings = QueryWarnings {
            before: &[b"Warning: q1\n".to_vec()],
            sink: &mut err,
        };
        write_blastp_tabular(&tab, hit_lists, &mut out, Some(&mut warnings), None).unwrap();
        (
            String::from_utf8(out).unwrap(),
            String::from_utf8(err).unwrap(),
        )
    }

    // NCBI tabular.cpp:1061-1091 (frames 1, `sstrand` N/A of a protein subject), blastfmtutil.cpp
    // :629-640 (`qcovs`, the union of the subject's HSPs; `qcovus` is blastn's only), tabular.cpp
    // :1264-1298 (outfmt 7's lines: psiblast's round, no count for an unsearched batch) and the
    // taxonomy database's warning of each formatter (tabular.cpp:1093-1104).
    #[test]
    fn rows_and_comment_lines_are_ncbis() {
        let lists = [
            Some(hit_list(
                vec![
                    vec![hit((1, 40), (1, 40), 0), hit((31, 60), (101, 130), 0)],
                    vec![hit((5, 14), (5, 14), 1)],
                ],
                1,
            )),
            None,
        ];
        let (out, err) = report(
            6,
            "qseqid sseqid slen qcovs qcovhsp qcovus sstrand frames stitle",
            BlastpApp::Blastp,
            &lists,
        );
        assert_eq!(
            out,
            "q1\ts1\t200\t60\t40\tN/A\tN/A\t1/1\tN/A\nq1\ts1\t200\t60\t30\tN/A\tN/A\t1/1\tN/A\n"
        );
        assert_eq!(err, "Warning: q1\n");
        let (out, err) = report(7, "qseqid sseqid sscinames", BlastpApp::Psiblast, &lists);
        assert_eq!(
            out,
            "# PSIBLAST 2.17.0+\n# Iteration: 1\n# Query: q1 first\n# Database: User specified sequence set (Input: s.faa)\n# Fields: query id, subject id, subject sci names\n# 2 hits found\nq1\ts1\tN/A\nq1\ts1\tN/A\n\
             # PSIBLAST 2.17.0+\n# Iteration: 1\n# Query: q2\n# Database: User specified sequence set (Input: s.faa)\n# BLAST processed 2 queries\n"
        );
        let warning = taxdb_warning_line("psiblast");
        assert_eq!(err, format!("Warning: q1\n{warning}{warning}{warning}"));
        let (out, _) = report(20, "qseqid evalue", BlastpApp::Blastp, &lists);
        assert_eq!(out, "qseqid,evalue\nq1,1.00e-10\nq1,1.00e-10\n");
    }
}
