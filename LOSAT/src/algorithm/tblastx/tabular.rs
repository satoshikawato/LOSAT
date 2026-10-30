//! TBLASTX's tabular reports with a field list or a delimiter (outfmt 6 and 7 other than the
//! default `std` fields with a tab) and outfmt 10 and 20: the adapter of the shared writer
//! (`report::tabular`, NCBI's `CBlastTabularInfo`). For each query in order it writes the query's
//! warnings, outfmt 7's comment lines and one row per HSP of the final hit list, in the order of
//! the default outfmt 6 rows (`report::write_tabular`); outfmt 20's header line comes before the
//! first query and outfmt 7's count of queries after the last. The report is written once per run
//! after every query batch was searched (`report_stream`'s `Off::CustomTabular`), so LOSAT's
//! refusals (a name of the taxonomy database, the ids of a database sequence) come before the
//! first byte of the report.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:763-834
//! ```c++
//!     CConstRef<CSeq_align_set> aln_set = results.GetSeqAlign();
//!     if (m_IsUngappedSearch && results.HasAlignments()) {
//!         aln_set.Reset(CDisplaySeqalign::PrepareBlastUngappedSeqalign(*aln_set));
//!     }
//!     ...
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
//! TBLASTX is ungapped: each HSP is a Seq-align of one Std-seg (`PrepareBlastUngappedSeqalign`
//! splits the subject's `diags` Seq-align), with the scores of `s_BuildScoreList` (the score, the
//! e-value or the sum e-value of a linked set, the bit score, the engine's identities and
//! `hsp_percent_coverage`). What `SetFields` reads of it is what the default row shows (the
//! nucleotide coordinates in the records, descending on a minus frame, the length, no gaps) and
//! the two rows of the alignment, both translated (`CreateTranslatedDensegFromNADenseg`, the query
//! with `-query_gencode`, the subject with `-db_gencode`): the rows of outfmt 0 and of the XML
//! without the masks' lower case (`SubjectLetters::rows_into`). TBLASTX does not set the
//! formatter's no-fetch mode, so `nident`, `mismatch` and `pident` read the rows too. `qcovus`
//! is BLASTN's alone (`N/A`).

use std::io::Write;
use std::sync::Arc;

use anyhow::Result;

use super::args::TblastxMatrix;
use super::blast_engine::{TblastxHsp, TblastxReportRun};
use crate::algorithm::common::tabular::{check_taxonomy_names, TabularSubjects, COUNTS_BLOCK};
use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::config::ScoringMatrix;
use crate::report::blast_xml::hsp::translation_frame;
use crate::report::blast_xml::matrix::FormattingMatrix;
use crate::report::query_warnings::QueryWarnings;
use crate::report::tabular::{
    coverage_percent, hsp_coverage, taxdb_warning_line, union_length, TabularField, TabularHsp,
    TabularKind, TabularQuery, TabularSpec, TabularWriter,
};
use crate::utils::genetic_code::GeneticCode;

const PROGRAM: &str = "TBLASTX";
/// The `# TBLASTX 2.17.0+` line of outfmt 7 (`NStr::ToUpper(m_Program) + " " + version`).
const PROGRAM_VERSION: &str = "TBLASTX 2.17.0+";

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:970-986
// ```c++
//     if (x_IsFieldRequested(eQuerySeq) ||
//         x_IsFieldRequested(eSubjectSeq) ||
//         x_IsFieldRequested(ePositives) ||
//         x_IsFieldRequested(ePercentPositives) ||
//         x_IsFieldRequested(eBTOP) ||
//         (x_IsFieldRequested(eNumIdentical) && !kNoFetchSequence) ||
//         (x_IsFieldRequested(eMismatches) && !kNoFetchSequence) ||
//         (x_IsFieldRequested(ePercentIdentical) && !kNoFetchSequence)) {
//
//         alnVec->SetGapChar('-');
//         alnVec->SetGenCode(m_QueryGeneticCode, 0);
//         alnVec->SetGenCode(m_DbGeneticCode, 1);
//         alnVec->GetWholeAlnSeqString(0, m_QuerySeq);
//         alnVec->GetWholeAlnSeqString(1, m_SubjectSeq);
// ```
/// The fields that make that block read the two rows, from which it counts the identities and
/// positives again (TBLASTX fetches the sequences: blast_format.cpp:791-792 sets no-fetch for
/// BLASTN alone).
const ROW_FIELDS: [TabularField; 8] = [
    TabularField::QuerySeq,
    TabularField::SubjectSeq,
    TabularField::Positives,
    TabularField::PercentPositives,
    TabularField::Btop,
    TabularField::NumIdentical,
    TabularField::Mismatches,
    TabularField::PercentIdentical,
];

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
        "tblastx",
        PROGRAM_VERSION,
        writer,
        diagnostics,
    )
}

/// One custom tabular report of a run: what `write_tblastx_outputs` gives it besides the run.
pub(super) struct TblastxTabular<'a, 'r> {
    pub spec: &'a TabularSpec,
    /// The final list after `shift_to_records` (record coordinates), in query order, the HSPs of
    /// a subject following each other.
    pub hits: &'a [TblastxHsp],
    pub run: &'a TblastxReportRun<'r>,
    /// The `# Query:` title of each query.
    pub query_titles: &'a [&'a [u8]],
    /// The `# Database:` name.
    pub database: &'a str,
    /// The default report's tabular subject ids (`sseqid`, `sacc`, `saccver`), by subject index
    /// or rank (`tabular_subject_ids`, `ReportedSubjects::tabular_subject_ids`).
    pub subject_ids: &'a [Arc<[u8]>],
    /// By query: whether NCBI did not search its batch (outfmt 7's lines without a count).
    pub unsearched: &'a [bool],
}

/// Writes the report of every query of `tab.run` to `writer`: outfmt 20's header line first
/// (for the run's first queries), then for each query its warnings (and the taxonomy database's
/// warning of its formatter) before its lines, and outfmt 7's last line after the last query when
/// `run.epilog`. `probe` marks each row (the format observer).
pub(super) fn write_tblastx_tabular(
    tab: &TblastxTabular<'_, '_>,
    writer: &mut dyn Write,
    mut warnings: Option<&mut QueryWarnings<'_>>,
    mut probe: Option<&mut FormatProbe<'_>>,
) -> Result<()> {
    let (spec, run, hits) = (tab.spec, tab.run, tab.hits);
    let args = run.args;
    let needs = spec.needs();
    let taxdb_found = args.tabular_taxdb.unwrap_or(false);
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
    let record_length = |s_idx: usize| {
        run.ranges
            .subject_placements
            .length(s_idx, run.subject_letters.length(s_idx))
    };
    let lengths: Vec<usize> = match &run.database {
        Some(database) if spec.has(TabularField::SubjectLength) => {
            (0..database.reported.records.len())
                .map(record_length)
                .collect()
        }
        _ => Vec::new(),
    };
    let subjects = match &run.database {
        None => TabularSubjects::Records {
            ids: tab.subject_ids,
            length: &record_length,
        },
        Some(database) => TabularSubjects::Database {
            ids: tab.subject_ids,
            ordinals: &database.reported.ordinals,
            records: &database.reported.records,
            parsed: &database.reported.parsed,
            lengths: &lengths,
            taxonomy: run.tabular_taxonomy,
        },
    };
    // The subjects of the report's rows in report order, each once per run of its HSPs.
    let reported_ranks = || {
        hits.iter()
            .enumerate()
            .filter(|&(index, hsp)| index == 0 || hits[index - 1].hit.s_idx != hsp.hit.s_idx)
            .map(|(_, hsp)| hsp.hit.s_idx as usize)
    };
    check_taxonomy_names(spec, taxdb_found, &subjects, reported_ranks(), PROGRAM)?;
    // The refusals of `sall*` (`ParsedSubject::all_tabular_ids`) come before the first byte;
    // those of the other ids are `ReportedSubjects::tabular_subject_ids`'.
    if needs.all_deflines {
        for rank in reported_ranks() {
            subjects.values(rank, spec, PROGRAM)?;
        }
    }
    let query_code = GeneticCode::from_id(args.query_gencode);
    let db_code = GeneticCode::from_id(args.db_gencode);
    // The positives of the rows: the formatter's `m_ScoringMatrix` (`GetAsciiProteinMatrix` of
    // `-matrix`). A name without a table stops the search before any report has a hit
    // (`run_batch_loop`).
    let matrix = match args.scoring_matrix() {
        TblastxMatrix::Table(matrix) => matrix,
        TblastxMatrix::Unknown | TblastxMatrix::Aborts => ScoringMatrix::Blosum62,
    };
    let positives = FormattingMatrix::ascii_protein(matrix);
    let table = TabularWriter {
        spec,
        program_version: PROGRAM_VERSION,
        matrix: Some(&positives),
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
        (spec.checks_taxdb() && !taxdb_found).then(|| taxdb_warning_line("tblastx"));
    let warn = |writer: &mut dyn Write, warnings: Option<&mut QueryWarnings<'_>>| -> Result<()> {
        if let (Some(line), Some(warnings)) = (taxdb_warning.as_deref(), warnings) {
            writer.flush()?;
            warnings.sink.write_all(line.as_bytes())?;
            warnings.sink.flush()?;
        }
        Ok(())
    };
    // Outfmt 20's line and its formatter's warning, before the first query's warnings (NCBI's
    // `PrintProlog` comes before it reads the first query batch); the custom reports are written
    // in one call (`report_stream`), so this is the run's first.
    // (`write_tblastx_outputs` gives the diagnostics wherever the warning is written.)
    if run.queries_before == 0 {
        let mut no_diagnostics = std::io::sink();
        let diagnostics: &mut dyn Write = match warnings.as_deref_mut() {
            Some(warnings) => &mut *warnings.sink,
            None => &mut no_diagnostics,
        };
        write_tabular_prolog(spec, taxdb_found, writer, diagnostics)?;
    }
    // The HSPs of each query in report order, numbered by their place in the final list (the
    // HSP numbers of the format observer, `report::write_tabular`).
    let mut by_query: Vec<Vec<HspIndex>> = vec![Vec::new(); run.queries.len()];
    for (hsp_index, hsp) in hits.iter().enumerate() {
        by_query[hsp.hit.q_idx as usize].push(hsp_index);
    }
    let reads_rows = spec.fields.iter().any(|field| ROW_FIELDS.contains(field));
    let counts_block = spec.fields.iter().any(|field| COUNTS_BLOCK.contains(field));
    let mut buf: Vec<u8> = Vec::new();
    let mut query_ranges: Vec<(u64, u64)> = Vec::new();
    let (mut window, mut qrow, mut srow): (Vec<u8>, Vec<u8>, Vec<u8>) =
        (Vec::new(), Vec::new(), Vec::new());
    for (q_idx, query) in run.queries.iter().enumerate() {
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
        // `qlen` and the `qcovs` of the subjects: the length of the `-query_loc` range as given,
        // else the nucleotide record's.
        let query_length = run.ranges.input.range_length.unwrap_or(searched_length);
        // The record's length, which the frames are read in (`GetFrame`).
        let query_record_length = run.ranges.input.placements.length(q_idx, searched_length);
        let tabular_query = TabularQuery {
            id: run.query_ids[q_idx].as_ref(),
            length: query_length,
        };
        let query_hits = &by_query[q_idx];
        buf.clear();
        // NCBI c++/src/algo/blast/api/local_blast.cpp:177-180,204-208: a skipped search has a
        // null alignment set, so outfmt 7 has no count line (`report::write_tabular`).
        table.query_header(
            &mut buf,
            tab.query_titles[q_idx],
            tab.database.as_bytes(),
            (!tab.unsearched[q_idx]).then_some(query_hits.len()),
        )?;
        if probe.is_some() {
            writer.write_all(&buf)?;
            buf.clear();
        }
        let mut first = 0usize;
        while first < query_hits.len() {
            let s_idx = hits[query_hits[first]].hit.s_idx;
            let end = first
                + query_hits[first..]
                    .iter()
                    .take_while(|&&index| hits[index].hit.s_idx == s_idx)
                    .count();
            let group = &query_hits[first..end];
            first = end;
            let s_idx = s_idx as usize;
            let values = subjects.values(s_idx, spec, PROGRAM)?;
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:2969-2985
            // ```c++
            //     ITERATE(CSeq_align_set::Tdata, iter, alnset.Get()) {
            //         CRange<TSeqPos> seq_range = (*iter)->GetSeqRange(0);
            //         //for minus strand
            //         if(seq_range.GetFrom() > seq_range.GetTo()){
            //             seq_range.Set(seq_range.GetTo(), seq_range.GetFrom());
            //         }
            //         temp.push_back(seq_range);
            //     }
            // ```
            // The HSPs of one subject are one group of `InsertSubjectScores` (consecutive
            // Seq-aligns of one subject); their query intervals in nucleotides (the Std-seg's
            // query location), which a `-query_loc` shift moves alike.
            let qcovs = if needs.qcovs {
                query_ranges.clear();
                query_ranges.extend(group.iter().map(|&index| {
                    let hit = &hits[index].hit;
                    (
                        hit.q_start.min(hit.q_end) as u64 - 1,
                        hit.q_start.max(hit.q_end) as u64 - 1,
                    )
                }));
                coverage_percent(union_length(&mut query_ranges), query_length)
            } else {
                -1
            };
            // `qcovus` is inserted for BLASTN alone (`eQueryCovPerUniqSubj`): `N/A`.
            let subject = values.nucleotide_subject(qcovs, -1);
            let subject_record_length = i32::try_from(record_length(s_idx))?;
            for &hsp_index in group {
                let hsp = &hits[hsp_index];
                let hit = &hsp.hit;
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:762-768
                // ```c++
                //         CAlignFormatUtil::GetAlnScores(align, score, bit_score, evalue, sum_n,
                //                                        num_ident, use_this_gi);
                // ```
                // `num_ident` decides whether `nident` alone enters the block above: the
                // Seq-align's (the engine's identities, `hsp_engine_counts`), counted only then.
                let engine_ident = if spec.has(TabularField::NumIdentical) && !counts_block {
                    run.subject_letters.engine_identities(
                        hsp,
                        query_letters,
                        &query_code,
                        &db_code,
                        &mut window,
                    )
                } else {
                    0
                };
                let rows_read = reads_rows && (counts_block || engine_ident > 0);
                if rows_read {
                    run.subject_letters.rows_into(
                        hsp,
                        query_letters,
                        &query_code,
                        &db_code,
                        &mut window,
                        &mut qrow,
                        &mut srow,
                    );
                }
                let (q_low, q_high) = (hit.q_start.min(hit.q_end), hit.q_start.max(hit.q_end));
                let (s_low, s_high) = (hit.s_start.min(hit.s_end), hit.s_start.max(hit.s_end));
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
                // The start of a minus strand is its higher end (the row's start), so both frames
                // are read in their whole records as the XML's `Hsp_query-frame` and
                // `Hsp_hit-frame` (`translation_frame`).
                let q_frame = translation_frame(
                    hit.query_frame > 0,
                    i32::try_from(q_low)?,
                    i32::try_from(q_high)?,
                    i32::try_from(query_record_length)?,
                );
                let s_frame = translation_frame(
                    hsp.subject_frame > 0,
                    i32::try_from(s_low)?,
                    i32::try_from(s_high)?,
                    subject_record_length,
                );
                let row = TabularHsp {
                    // The record coordinates of the default row (`shift_to_records`).
                    q_start: i64::try_from(hit.q_start)?,
                    q_end: i64::try_from(hit.q_end)?,
                    s_start: i64::try_from(hit.s_start)?,
                    s_end: i64::try_from(hit.s_end)?,
                    // The e-value of the default row: a linked set's sum e-value (`sum_e`).
                    evalue: hit.e_value,
                    bit_score: hit.bit_score,
                    score: hit.raw_score,
                    // An ungapped HSP: its residues, no gap.
                    length: hit.length,
                    gap_opens: 0,
                    gaps: 0,
                    // Counted again from the rows when they are read (the writer).
                    num_ident: engine_ident,
                    q_frame,
                    s_frame,
                    rows: rows_read.then_some((qrow.as_slice(), srow.as_slice())),
                    qcovhsp: if group.len() > 1 {
                        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1222-1225
                        // ```c++
                        //     if(query_length > 0) {
                        //             static const string kQueryCovHsp("hsp_percent_coverage");
                        //             double hsp_coverage = Blast_HSPGetQueryCoverage( hsp, query_length);
                        // ```
                        // The Seq-align of an HSP split from its subject's (the scores of its
                        // Std-seg): `Blast_HSPGetQueryCoverage` divides the HSP's residues in its
                        // query frame by the length that `BlastHitList2SeqAlign_OMF` gives, the
                        // searched query's in nucleotides (`query.GetSeqLength`), about a third
                        // of the HSP's nucleotide coverage, as NCBI writes it.
                        hsp_coverage(
                            i32::try_from(hit.sort_query_end - hit.sort_query_offset)?,
                            i32::try_from(searched_length)?,
                        )
                    } else {
                        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:737-746
                        // ```c++
                        // void CBlastTabularInfo::x_SetQueryCovSeqalign(const CSeq_align & align, int query_len)
                        // {
                        // 	double tmp = 0;
                        // 	if(!align.GetNamedScore("hsp_percent_coverage", tmp)) {
                        // 		int len = abs((int) (align.GetSeqStop(0) - align.GetSeqStart(0))) + 1;
                        // 		tmp = 100.0 * len/(double) query_len;
                        // 		if(tmp  < 99)
                        // 			tmp +=0.5;
                        // 	}
                        // 	m_QueryCovSeqalign = (int)tmp;
                        // ```
                        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3170-3183
                        // ```c++
                        //         if(seg.Which() == CSeq_align::C_Segs::e_Std){
                        //             if(seg.GetStd().size() > 1){
                        //     ...
                        //                     if((*iterStdseg)->IsSetScores()){
                        //                         aln->SetScore() = (*iterStdseg)->GetScores();
                        //                     }
                        //     ...
                        //             } else {
                        //                 alnSetRef->Set().push_back(*iter);
                        //             }
                        // ```
                        // A subject with one HSP keeps its `diags` Seq-align, whose scores are on
                        // its Std-seg, and `GetNamedScore` reads the Seq-align's own: the HSP's
                        // nucleotide span over the query length of `qlen`.
                        hsp_coverage(
                            i32::try_from(q_high - q_low + 1)?,
                            i32::try_from(query_length)?,
                        )
                    },
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
    // A search that stops at a query batch writes no epilog (`report::write_tabular`).
    if spec.kind == TabularKind::Commented && run.epilog {
        warn(&mut *writer, warnings.as_deref_mut())?;
        table.epilog(writer, run.queries_before + run.queries.len())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI blast_format.cpp:348-357: outfmt 20's line of the field names, after the warning of its
    // formatter when `taxdb.bti` is not found and a plural name is shown.
    #[test]
    fn outfmt_20_writes_its_line_and_its_formatter_s_warning() {
        let spec = |outfmt: &str| match crate::algorithm::tblastx::report::output_format(outfmt) {
            crate::algorithm::tblastx::report::TblastxOutputFormat::CustomTabular(spec) => spec,
            other => panic!("{outfmt}: {other:?}"),
        };
        let (mut out, mut err) = (Vec::new(), Vec::new());
        write_tabular_prolog(spec("20 qseqid sscinames"), false, &mut out, &mut err).unwrap();
        assert_eq!(out, b"qseqid,sscinames\n");
        assert_eq!(err, taxdb_warning_line("tblastx").as_bytes());
        let (mut out, mut err) = (Vec::new(), Vec::new());
        write_tabular_prolog(spec("20 delim=| qseqid"), true, &mut out, &mut err).unwrap();
        assert_eq!(
            (out.as_slice(), err.as_slice()),
            (&b"qseqid\n"[..], &b""[..])
        );
        let (mut out, mut err) = (Vec::new(), Vec::new());
        write_tabular_prolog(spec("10 sscinames"), false, &mut out, &mut err).unwrap();
        assert!(out.is_empty() && err.is_empty());
    }

    // NCBI tabular.cpp:900-908 and 970-986: `nident` alone enters the block (and reads the rows)
    // with a positive engine count only; the other count fields always, the row fields read the
    // rows.
    #[test]
    fn the_counts_block_and_the_row_fields_are_ncbis() {
        for field in ROW_FIELDS {
            assert!(
                field == TabularField::NumIdentical || COUNTS_BLOCK.contains(&field),
                "{field:?}"
            );
        }
        assert!(!COUNTS_BLOCK.contains(&TabularField::NumIdentical));
        assert!(!COUNTS_BLOCK.contains(&TabularField::Frames));
        assert!(!COUNTS_BLOCK.contains(&TabularField::SubjectFrame));
    }
}
