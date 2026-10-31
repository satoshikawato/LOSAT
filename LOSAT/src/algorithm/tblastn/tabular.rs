//! TBLASTN's tabular reports with a field list or a delimiter (outfmt 6 and 7 other than the
//! default `std` fields with a tab) and outfmt 10 and 20: the adapter of the shared writer
//! (`report::tabular`, NCBI's `CBlastTabularInfo`). For each query in order it writes the query's
//! warnings, outfmt 7's comment lines and one row per HSP of the query's final hit list, in the
//! order of the default outfmt 6 rows (`stage_e_report::write_tabular`); outfmt 20's header line
//! comes before the first query and outfmt 7's count of queries after the last. `render` calls it
//! once per run, after every query batch was searched, so LOSAT's refusals (the taxonomy database's
//! names, `TabularInputs::names`, the ids of a database sequence) come before the first byte of the
//! report.
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
//!             tabinfo.SetQueryGeneticCode(m_QueryGenCode);
//!             tabinfo.SetDbGeneticCode(m_DbGenCode);
//!             ITERATE(CSeq_align_set::Tdata, itr, copy_aln_set.Get()) {
//!                     const CSeq_align& s = **itr;
//!                     tabinfo.SetFields(s, *m_Scope, &m_ScoringMatrix);
//!                     tabinfo.Print();
//!             }
//! ```
//! What `SetFields` reads of an HSP is its Seq-align (`blast_xml` has the coordinates): the protein
//! query ascending, the subject in nucleotides (descending on a minus frame), the score, the
//! e-value (the sum e-value of a linked set) and the bit score, the alignment's lengths and gaps,
//! and the two rows of the alignment, which the TBLASTN formatter reads (no `SetNoFetch`) for
//! `qseq`, `sseq`, `btop`, `positive`, `ppos` and the identities: the protein query, unmasked, and
//! the subject translated with `-db_gencode` (the XML's rows without its masks,
//! `aligned_sequences_into`). Without those fields the identities are the default row's
//! (`report_num_ident`, counted from the same rows). `qcovus` is BLASTN's alone (`N/A`).

use std::io::Write;
use std::sync::Arc;

use anyhow::{ensure, Context, Result};

use super::search_gapped::TargetTranslation;
use super::stage_d_results::KappaResultHitList;
use super::stage_e_report::{
    aligned_sequences_into, record_coordinates, ReportRanges, ShownSubjects,
};
use super::x_subject_store::SubjectLetters;
use crate::algorithm::common::tabular::{TabularSubjects, COUNTS_BLOCK};
use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::config::ScoringMatrix;
use crate::report::blast_xml::hsp::translation_frame;
use crate::report::query_warnings::QueryWarnings;
use crate::report::tabular::{
    coverage_percent, hsp_coverage, taxdb_warning_line, union_length, TabularField, TabularHsp,
    TabularKind, TabularNames, TabularQuery, TabularSpec, TabularTaxonomy, TabularWriter,
};
use crate::utils::genetic_code::GeneticCode;

const PROGRAM: &str = "TBLASTN";
/// The `# TBLASTN 2.17.0+` line of outfmt 7 (`NStr::ToUpper(m_Program) + " " + version`).
const PROGRAM_VERSION: &str = "TBLASTN 2.17.0+";

/// Whether a report of `spec` can read the rows of an HSP (`write_tblastn_tabular`: the letters of
/// its subject): its fields read them, or `nident` is the only field of the counts' block, which
/// reads them for some HSPs.
pub(super) fn reads_rows(spec: &TabularSpec) -> bool {
    spec.needs().rows
        || (spec.has(TabularField::NumIdentical)
            && !spec.fields.iter().any(|field| COUNTS_BLOCK.contains(field)))
}

/// What the custom tabular reports of a run need beyond `render`'s other inputs, made by the
/// search before the reports.
#[derive(Default)]
pub(super) struct TabularInputs {
    /// Whether `taxdb.bti` is on NCBI's search path (`TblastnArgs::tabular_taxdb`; not found
    /// where it was not asked).
    pub taxdb_found: bool,
    /// `-db`: the taxonomy of the reported sequences by rank, read where a field shows it
    /// (`TblastnDatabase::tabular_taxonomy`); empty otherwise.
    pub taxonomy: Vec<TabularTaxonomy>,
    /// `-db`: the names of the taxonomy database that the rows of the run's custom report (the
    /// command line's one `-outfmt`) look up (`TblastnDatabase::tabular_names`); `None`: every
    /// lookup fails (no name field, no `taxdb.bti`, `-subject`).
    pub names: Option<TabularNames>,
    /// The length of the `-query_loc` range as given (NCBI's `m_QueryRange`), `None` without it.
    pub range_length: Option<usize>,
}

/// One custom tabular report of a run: `render`'s inputs that it reads.
pub(super) struct TblastnTabular<'a> {
    pub spec: &'a TabularSpec,
    pub inputs: &'a TabularInputs,
    /// The queries (their searched letters, which the rows read).
    pub queries: &'a [FastaRecord],
    /// The letters of the subjects by subject index (`-subject`) or rank (`-db`).
    pub letters: SubjectLetters<'a>,
    /// The default report's tabular subject ids (`sseqid`, `sacc`, `saccver`), by subject index
    /// or rank (`stage_e_report::tabular_subject_ids`, `ReportedSubjects::tabular_subject_ids`).
    pub subject_ids: &'a [Arc<[u8]>],
    pub shown: &'a ShownSubjects<'a>,
    /// By query: whether NCBI did not search its batch (outfmt 7's lines without a count).
    pub query_batch_skipped: &'a [bool],
    pub ranges: &'a ReportRanges<'a>,
    /// `-db_gencode`: the code that translates the subject rows (`SetDbGeneticCode`).
    pub genetic_code: u8,
    /// The run's matrix (`m_ScoringMatrix`, the positives).
    pub matrix: ScoringMatrix,
}

/// The custom tabular report of a checked `-outfmt` value (`TabularSpec::of`): a field list or a
/// delimiter of outfmt 6 and 7, and every outfmt 10 and 20; `None` for the default 6 and 7 (the
/// `std` fields with a tab, `stage_e_report::write_tabular`) and the other formats.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:770-786
/// ```c++
///     if (m_FormatType == CFormattingArgs::eTabular ||
///         m_FormatType == CFormattingArgs::eTabularWithComments ||
///         m_FormatType == CFormattingArgs::eCommaSeparatedValues ||
///         m_FormatType == CFormattingArgs::eCommaSeparatedValuesWithHeader) {
///     ...
///         CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec, kDelim);
///         if(!m_CustomDelim.empty()) {
///             tabinfo.SetCustomDelim(m_CustomDelim);
///         }
/// ```
pub(super) fn tabular_spec(outfmt: &str) -> Option<TabularSpec> {
    let choice = crate::blastinput::app::parse_formatting_string(outfmt).ok()?;
    TabularSpec::of(
        crate::blastinput::app::ReportFormat::from_number(choice.number)?,
        &choice,
    )
}

/// Writes outfmt 20's line of the field names and its formatter's warning
/// (`common::tabular::write_tabular_prolog`); nothing for the other kinds.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:348-357
/// ```c++
///     if (m_FormatType == CFormattingArgs::eCommaSeparatedValuesWithHeader) {
///     ...
///         tabinfo.PrintFieldSpecs();
/// ```
pub(super) fn write_tabular_prolog(
    spec: &TabularSpec,
    taxdb_found: bool,
    writer: &mut dyn Write,
    diagnostics: &mut dyn Write,
) -> Result<()> {
    crate::algorithm::common::tabular::write_tabular_prolog(
        spec,
        taxdb_found,
        "tblastn",
        PROGRAM_VERSION,
        writer,
        diagnostics,
    )
}

/// The subjects of the report's rows in report order (subject index or rank of each hit list
/// with HSPs).
fn reported_ranks(hitlists: &[KappaResultHitList]) -> impl Iterator<Item = usize> + '_ {
    hitlists
        .iter()
        .flat_map(KappaResultHitList::lists)
        .filter(|list| !list.hsps.hsps.is_empty())
        .map(|list| list.oid as usize)
}

/// Writes the report of every query of `hitlists` (indexed by query, sorted for the report by
/// `render`) to `writer`: outfmt 20's header line first, then for each query its warnings (and the
/// taxonomy database's warning of its formatter) before its lines, and outfmt 7's last line after
/// the last query when `ranges.epilog`. `probe` marks each row (the format observer).
pub(super) fn write_tblastn_tabular(
    tab: &TblastnTabular<'_>,
    hitlists: &[KappaResultHitList],
    writer: &mut dyn Write,
    mut warnings: Option<&mut QueryWarnings<'_>>,
    mut probe: Option<&mut FormatProbe<'_>>,
) -> Result<()> {
    ensure!(
        tab.queries.len() == hitlists.len() && tab.query_batch_skipped.len() == hitlists.len(),
        "TBLASTN tabular report query count mismatch"
    );
    let spec = tab.spec;
    let needs = spec.needs();
    let counts_block = spec.fields.iter().any(|field| COUNTS_BLOCK.contains(field));
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-842
    // ```c++
    //        		const CBioseq_Handle& subject_bh =
    //                 scope.GetBioseqHandle(align.GetSeq_id(1));
    //             if(setSubjectId) {
    //             	SetSubjectId(subject_bh);
    //             }
    //             m_SubjectLength = subject_bh.GetBioseqLength();
    // ```
    // `slen` is the length of the subject's Bioseq: the whole record with `-subject_loc` (the
    // XML's `Hit_len`), a database sequence's stored length.
    let record_length = |s_idx: usize| tab.ranges.subjects.length(s_idx, tab.letters.length(s_idx));
    let lengths: Vec<usize> = match tab.shown {
        ShownSubjects::Database { reported, .. } if spec.has(TabularField::SubjectLength) => {
            (0..reported.records.len()).map(record_length).collect()
        }
        _ => Vec::new(),
    };
    let subjects = match tab.shown {
        ShownSubjects::Fasta { .. } => TabularSubjects::Records {
            ids: tab.subject_ids,
            length: &record_length,
        },
        ShownSubjects::Database { reported, .. } => TabularSubjects::Database {
            ids: tab.subject_ids,
            ordinals: &reported.ordinals,
            records: &reported.records,
            parsed: &reported.parsed,
            lengths: &lengths,
            taxonomy: &tab.inputs.taxonomy,
        },
    };
    let names = tab.inputs.names.as_ref();
    // The refusals of `sall*` (`ParsedSubject::all_tabular_ids`) come before the first byte;
    // those of the other ids are `ReportedSubjects::tabular_subject_ids`'.
    if needs.all_deflines {
        for rank in reported_ranks(hitlists) {
            subjects.values(rank, spec, PROGRAM)?;
        }
    }
    let code = GeneticCode::try_from_id(tab.genetic_code).map_err(anyhow::Error::msg)?;
    // The positives of the rows: the formatter's `m_ScoringMatrix` (`GetAsciiProteinMatrix`).
    let table = TabularWriter {
        spec,
        program_version: PROGRAM_VERSION,
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
        (spec.checks_taxdb() && !tab.inputs.taxdb_found).then(|| taxdb_warning_line("tblastn"));
    let warn = |writer: &mut dyn Write, warnings: Option<&mut QueryWarnings<'_>>| -> Result<()> {
        if let (Some(line), Some(warnings)) = (taxdb_warning.as_deref(), warnings) {
            writer.flush()?;
            warnings.sink.write_all(line.as_bytes())?;
            warnings.sink.flush()?;
        }
        Ok(())
    };
    // Outfmt 20's line and its formatter's warning, before the first query's warnings (NCBI's
    // `PrintProlog` comes before it reads the first query batch).
    // (`render` gives the diagnostics wherever the warning is written.)
    let mut no_diagnostics = std::io::sink();
    let diagnostics: &mut dyn Write = match warnings.as_deref_mut() {
        Some(warnings) => &mut *warnings.sink,
        None => &mut no_diagnostics,
    };
    write_tabular_prolog(spec, tab.inputs.taxdb_found, writer, diagnostics)?;
    let mut buf: Vec<u8> = Vec::new();
    let mut query_ranges: Vec<(u64, u64)> = Vec::new();
    let mut x_window: Vec<u8> = Vec::new();
    let (mut qrow, mut srow) = (String::new(), String::new());
    // The rows are printed in the order of the final HSP lists, so a running count is each
    // row's HSP index.
    let mut hsp_index: HspIndex = 0;
    for (q_idx, (query, hitlist)) in tab.queries.iter().zip(hitlists).enumerate() {
        // The query's warnings come before its lines (`QueryWarnings`), then the warning of its
        // formatter.
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        warn(&mut *writer, warnings.as_deref_mut())?;
        let query_letters = query.seq();
        let searched_length = query_letters.len();
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:787-790
        // ```c++
        //             if(m_QueryRange.NotEmpty())
        //             	m_QueryLength = m_QueryRange.GetLength();
        //             else
        //             	m_QueryLength = query_bh.GetBioseqLength();
        // ```
        let query_length = tab.inputs.range_length.unwrap_or(searched_length);
        let tabular_query = TabularQuery::of_ids(query.believed(), query.shown_id(), query_length);
        let hit_count: usize = hitlist
            .lists()
            .iter()
            .map(|list| list.hsps.hsps.len())
            .sum();
        buf.clear();
        // NCBI c++/src/algo/blast/api/local_blast.cpp:177-180,204-208: a skipped search has a
        // null alignment set, so outfmt 7 has no count line (`write_tabular`).
        table.query_header(
            &mut buf,
            &query.query_line(),
            tab.shown.name().as_bytes(),
            (!tab.query_batch_skipped[q_idx]).then_some(hit_count),
        )?;
        if probe.is_some() {
            writer.write_all(&buf)?;
            buf.clear();
        }
        for list in hitlist.lists() {
            if list.hsps.hsps.is_empty() {
                continue;
            }
            let oid = usize::try_from(list.oid).context("negative TBLASTN subject OID")?;
            ensure!(
                list.payloads.len() == list.hsps.hsps.len(),
                "TBLASTN report payload count mismatch"
            );
            let values = subjects.values(oid, spec, PROGRAM)?;
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
            // subject once); their query intervals in residues, which a `-query_loc` shift moves
            // alike.
            let qcovs = if needs.qcovs {
                query_ranges.clear();
                query_ranges.extend(list.hsps.hsps.iter().map(|linked| {
                    (
                        linked.hsp.q_start as u64,
                        (linked.hsp.q_end as u64).saturating_sub(1),
                    )
                }));
                coverage_percent(union_length(&mut query_ranges), query_length)
            } else {
                -1
            };
            // `qcovus` is inserted for BLASTN alone (`eQueryCovPerUniqSubj`): `N/A`.
            let subject = values.nucleotide_subject(names, qcovs, -1);
            let subject_len = i32::try_from(tab.letters.length(oid))?;
            let whole_length = i32::try_from(record_length(oid))?;
            let mut translation = None;
            for (linked, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                // The texts of the row's lookups of the taxonomy database (`x_SetTaxInfoAll`,
                // `x_SetTaxInfo`) follow the rows before it.
                if let (Some(names), Some(warnings)) = (names, warnings.as_deref_mut()) {
                    let text = names.row_messages(spec, subject.taxonomy);
                    if !text.is_empty() {
                        writer.write_all(&buf)?;
                        buf.clear();
                        writer.flush()?;
                        warnings.sink.write_all(&text)?;
                        warnings.sink.flush()?;
                    }
                }
                let hsp = &linked.hsp;
                let frame = i32::from(hsp.frame);
                // NCBI c++/src/algo/blast/core/blast_hits.c:1087-1105: the nucleotide
                // coordinates of a translated frame, descending on a minus frame (the default
                // row's, `write_tabular`).
                let (s_start, s_end) = if frame > 0 {
                    (3 * hsp.s_start + frame, 3 * hsp.s_end + frame - 1)
                } else {
                    (
                        subject_len - 3 * hsp.s_start + frame + 1,
                        subject_len - 3 * hsp.s_end + frame + 2,
                    )
                };
                let (q_start, q_end, s_start, s_end) = record_coordinates(
                    tab.ranges,
                    q_idx,
                    oid,
                    (hsp.q_start + 1, hsp.q_end, s_start, s_end),
                )?;
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-908
                // ```c++
                //         x_IsFieldRequested(eSubjectSeq) || (x_IsFieldRequested(eNumIdentical) && num_ident > 0 ) ||
                // ```
                // `num_ident` there is the Seq-align's `num_ident` score (the engine's identities,
                // `payload.num_ident`), which can differ from the report's identities. With
                // `nident` and no other field of the block, NCBI counts the identities only
                // when the engine's are positive and prints 0 otherwise. The shared writer
                // decides with `TabularHsp::num_ident` and counts again from the rows when
                // they are given: where the engine's identities are positive and the report's
                // are 0, the rows are read and `num_ident` is the engine's.
                let engine_ident = usize::try_from(payload.num_ident).unwrap_or(0);
                let nident_alone = spec.has(TabularField::NumIdentical) && !counts_block;
                let rows_read = needs.rows
                    || (nident_alone && engine_ident > 0 && payload.report_num_ident == 0);
                if rows_read {
                    let translation = translation
                        .get_or_insert_with(|| TargetTranslation::of(tab.letters.nt(oid), &code));
                    // The rows without the lower case of the masks (`CAlnVec` of the query's
                    // Bioseq, whose letters the FASTA reader made upper case).
                    aligned_sequences_into(
                        query_letters,
                        tab.letters.nt(oid),
                        &mut x_window,
                        &code,
                        translation,
                        hsp.frame,
                        hsp.q_start,
                        hsp.s_start,
                        hsp.s_end,
                        &payload.edit_script,
                        &[],
                        &mut qrow,
                        &mut srow,
                    )?;
                }
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1072-1084
                // ```c++
                //     if (kTranslated) {
                //         if (x_IsFieldRequested(eQueryFrame) || x_IsFieldRequested(eFrames)) {
                //             query_frame = CAlignFormatUtil::
                //                 GetFrame (q_start - 1, ds.GetSeqStrand(kQueryRow),
                //                           scope.GetBioseqHandle(align.GetSeq_id(0)));
                //         }
                //
                //         if (x_IsFieldRequested(eSubjFrame) || x_IsFieldRequested(eFrames)) {
                //             subject_frame = CAlignFormatUtil::
                //                 GetFrame (s_start - 1, ds.GetSeqStrand(kSubjectRow),
                //                           scope.GetBioseqHandle(align.GetSeq_id(1)));
                //         }
                // ```
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1233-1244
                // ```c++
                //     int frame = 0;
                //     if (strand == eNa_strand_plus) {
                //         frame = (start % 3) + 1;
                //     } else if (strand == eNa_strand_minus) {
                //         frame = -(((int)handle.GetBioseqLength() - start - 1)
                //                   % 3 + 1);
                // ```
                // The protein query has no strand (frame 0); the subject's frame is read in its
                // whole record, as the XML's `Hit_frame` (`translation_frame`).
                let s_frame = translation_frame(
                    frame > 0,
                    s_start.min(s_end),
                    s_start.max(s_end),
                    whole_length,
                );
                let row = TabularHsp {
                    q_start: i64::from(q_start),
                    q_end: i64::from(q_end),
                    s_start: i64::from(s_start),
                    s_end: i64::from(s_end),
                    // The e-value of the default row: a linked set's sum e-value (`sum_e`).
                    evalue: linked.evalue,
                    bit_score: payload.bit_score,
                    score: hsp.score,
                    length: payload.align_length,
                    gap_opens: payload.gap_opens,
                    gaps: payload.gap_letters,
                    num_ident: if rows_read || (nident_alone && engine_ident == 0) {
                        engine_ident
                    } else {
                        payload.report_num_ident
                    },
                    q_frame: 0,
                    s_frame,
                    rows: rows_read.then_some((qrow.as_bytes(), srow.as_bytes())),
                    // `Blast_HSPGetQueryCoverage` of the HSP's context, which has the searched
                    // length (blast_seqalign.cpp:1222-1225, `hsp_percent_coverage`).
                    qcovhsp: hsp_coverage(hsp.q_end - hsp.q_start, i32::try_from(searched_length)?),
                };
                match probe.as_deref_mut() {
                    // One printed row is one HSP; the probe marks it without changing it.
                    Some(probe) => {
                        buf.clear();
                        table.row(&mut buf, &tabular_query, &subject, &row);
                        writer.flush()?;
                        probe.begin(hsp_index);
                        writer.write_all(&buf)?;
                        writer.flush()?;
                        probe.end(hsp_index);
                        buf.clear();
                    }
                    None => table.row(&mut buf, &tabular_query, &subject, &row),
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
    // A search that stops at a query batch writes no epilog (`write_tabular`).
    if spec.kind == TabularKind::Commented && tab.ranges.epilog {
        warn(&mut *writer, warnings.as_deref_mut())?;
        table.epilog(writer, tab.queries.len())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI blast_format.cpp:770-786: every tabular kind goes to `CBlastTabularInfo`; LOSAT's
    // default 6 and 7 (the `std` fields with a tab) keep `stage_e_report::write_tabular`.
    #[test]
    fn only_a_field_list_a_delimiter_and_10_and_20_take_the_shared_writer() {
        for outfmt in [
            "6",
            "7",
            "6 std",
            " 7 std ",
            "6 -qseqid",
            "6 delim=",
            "0 qseqid",
            "5",
        ] {
            assert!(tabular_spec(outfmt).is_none(), "{outfmt}");
        }
        for (outfmt, kind) in [
            ("6 qseqid", TabularKind::Tabular),
            ("6 delim=, std", TabularKind::Tabular),
            ("7 qaccver foo", TabularKind::Commented),
            ("10", TabularKind::Csv),
            ("20 qseqid", TabularKind::CsvWithHeader),
        ] {
            assert_eq!(
                tabular_spec(outfmt).map(|spec| spec.kind),
                Some(kind),
                "{outfmt}"
            );
        }
    }

    #[test]
    fn outfmt_20_writes_its_line_and_its_formatter_s_warning() {
        let spec = tabular_spec("20 qseqid sscinames").unwrap();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        write_tabular_prolog(&spec, false, &mut out, &mut err).unwrap();
        assert_eq!(out, b"qseqid,sscinames\n");
        assert_eq!(err, taxdb_warning_line("tblastn").as_bytes());
        let (mut out, mut err) = (Vec::new(), Vec::new());
        write_tabular_prolog(&spec, true, &mut out, &mut err).unwrap();
        assert_eq!(
            (out.as_slice(), err.as_slice()),
            (&b"qseqid,sscinames\n"[..], &b""[..])
        );
        let (mut out, mut err) = (Vec::new(), Vec::new());
        write_tabular_prolog(
            &tabular_spec("10 sscinames").unwrap(),
            false,
            &mut out,
            &mut err,
        )
        .unwrap();
        assert!(out.is_empty() && err.is_empty());
    }
}
