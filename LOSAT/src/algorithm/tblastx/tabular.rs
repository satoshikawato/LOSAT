//! TBLASTX's tabular reports with a field list or a delimiter (outfmt 6 and 7 other than the
//! default `std` fields with a tab) and outfmt 10 and 20: the adapter of the shared writer
//! (`report::tabular`, NCBI's `CBlastTabularInfo`). For each query in order it writes the query's
//! warnings, outfmt 7's comment lines and one row per HSP of the final hit list, in the order of
//! the default outfmt 6 rows (`report::write_tabular`); outfmt 20's header line comes before the
//! first query and outfmt 7's count of queries after the last. The report is written once per run
//! after every query batch was searched (`report_stream`'s `Off::CustomTabular`), so LOSAT's
//! refusals (the taxonomy database's names, `TblastxReportRun::tabular_names`, the ids of a
//! database sequence) come before the first byte of the report.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:762-834
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

use anyhow::Result;

use super::args::TblastxMatrix;
use super::blast_engine::{TblastxHsp, TblastxReportRun};
use crate::algorithm::common::tabular::{TabularSubjects, COUNTS_BLOCK};
use crate::api::local_blast::{FormatProbe, HspIndex};
use crate::config::ScoringMatrix;
use crate::report::blast_xml::hsp::translation_frame;
use crate::report::query_warnings::QueryWarnings;
use crate::report::tabular::{
    coverage_percent, union_length, QueryCoveragePerHsp, TabularField, TabularHsp, TabularKind,
    TabularQuery, TabularRows, TabularSpec, TabularWriter,
};
use crate::utils::genetic_code::GeneticCode;

const PROGRAM: &str = "TBLASTX";
/// The `# TBLASTX 2.17.0+` line of outfmt 7 (`NStr::ToUpper` of `m_Program`, a space and `version`).
const PROGRAM_VERSION: &str = concat!("TBLASTX ", crate::wording::report_version!());

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:970-986
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

/// Whether the custom tabular report of `fields` reads the letters of the subjects: a field of
/// `ROW_FIELDS` or of `COUNTS_BLOCK` (`write_tblastx_tabular`). Its other fields read the HSP,
/// the titles and the lengths only.
pub(crate) fn reads_letters(fields: &[TabularField]) -> bool {
    fields
        .iter()
        .any(|field| ROW_FIELDS.contains(field) || COUNTS_BLOCK.contains(field))
}

/// Writes outfmt 20's line of the field names and its formatter's warning
/// (`common::tabular::write_tabular_prolog`); nothing for the other kinds.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:348-357
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
    pub subject_ids: &'a [crate::algorithm::common::tabular::SubjectId],
    /// By query: whether NCBI did not search its batch (outfmt 7's lines without a count).
    pub unsearched: &'a [bool],
}

/// Writes the report of every query of `tab.run` to `writer`: outfmt 20's header line first
/// (for the run's first queries), then for each query its warnings (and the taxonomy database's
/// warning of its formatter) before its lines, and outfmt 7's last line after the last query when
/// `run.epilog`. `probe` marks each row (the format observer).
pub(super) fn write_tblastx_tabular(
    tab: &TblastxTabular<'_, '_>,
    writer: &mut crate::cli::NcbiOstream<'_>,
    mut warnings: Option<&mut QueryWarnings<'_>>,
    mut probe: Option<&mut FormatProbe<'_>>,
) -> Result<()> {
    let (spec, run, hits) = (tab.spec, tab.run, tab.hits);
    let args = run.args;
    let needs = spec.needs();
    let taxdb_found = args.tabular_taxdb.unwrap_or(false);
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-842
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
    // The records of a subject set whose ids `-parse_deflines` read, by subject index.
    let records_parsed = if run.args.parse_deflines && run.database.is_none() {
        crate::report::db_subject::parsed_records(
            run.subjects
                .iter()
                .map(|record| (record.believed(), record.title.as_slice())),
            run.args.show_gis,
        )
    } else {
        Vec::new()
    };
    let subjects = match &run.database {
        None => TabularSubjects::Records {
            ids: tab.subject_ids,
            length: &record_length,
            parsed: &records_parsed,
            records: run.subjects,
        },
        Some(database) => TabularSubjects::Database {
            ids: tab.subject_ids,
            ordinals: &database.reported.ordinals,
            oids: database.oids,
            records: &database.reported.records,
            parsed: &database.reported.parsed,
            sources: &database.reported.sources,
            lengths: &lengths,
            taxonomy: run.tabular_taxonomy,
            loads: &database.reported.loads,
        },
    };
    // The subjects of the report's rows in report order, each once per run of its HSPs.
    let reported_ranks = || {
        hits.iter()
            .enumerate()
            .filter(|&(index, hsp)| index == 0 || hits[index - 1].hit.s_idx != hsp.hit.s_idx)
            .map(|(_, hsp)| hsp.hit.s_idx as usize)
    };
    let names = run.tabular_names;
    // Outfmt 18, the organism report: the `Query=` block shows the query record's length, also
    // with `-query_loc` (as outfmt 0's).
    if spec.kind == crate::report::tabular::TabularKind::TaxReport {
        let mut queries: Vec<crate::report::tax_report::TaxQuery<'_>> = run
            .queries
            .iter()
            .enumerate()
            .map(|(q_idx, query)| crate::report::tax_report::TaxQuery {
                title: tab.query_titles[q_idx],
                length: run.ranges.input.placements.length(q_idx, query.seq().len()),
                // A query of a batch that NCBI did not search has no Seq-align-set.
                hsps: (!tab.unsearched[q_idx]).then(Vec::new),
            })
            .collect();
        // The HSPs of each query (a query that NCBI did not search has none).
        for hsp in hits {
            if let Some(hsps) = &mut queries[hsp.hit.q_idx as usize].hsps {
                hsps.push(crate::report::tax_report::TaxHsp {
                    subject: hsp.hit.s_idx as usize,
                    evalue: hsp.hit.e_value,
                    bit_score: hsp.hit.bit_score,
                });
            }
        }
        return crate::algorithm::common::tabular::write_tax_report(
            &subjects,
            run.subjects,
            false,
            names,
            PROGRAM,
            &queries,
            writer,
            warnings,
            args.html,
        );
    }
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
    let table = TabularWriter::new(spec, PROGRAM_VERSION, Some(matrix));
    // NCBI makes a formatter for each query's report and for outfmt 7's epilog, and each warns
    // when the field list asks for a name of the taxonomy database that it cannot find.
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:139-158
    let taxdb_warning = spec.formatter_warning(taxdb_found, "tblastx");
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
    // Which rows read their subject (`TabularSpec::subject_fetch`).
    let fetch = spec.subject_fetch(crate::blastinput::app::Program::Tblastx);
    // The bytes of one comment block or one row, written to the report stream as soon as they
    // are made (as NCBI's insertions reach its stream): nothing of a query's report is held.
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:595-605
        // `qlen` and the `qcovs` of the subjects: the length of the `-query_loc` range as given,
        // else the nucleotide record's.
        let query_length = run.ranges.input.range_length.unwrap_or(searched_length);
        // The record's length, which the frames are read in (`GetFrame`).
        let query_record_length = run.ranges.input.placements.length(q_idx, searched_length);
        let tabular_query = TabularQuery::of_ids(
            run.queries[q_idx].believed(),
            run.query_ids[q_idx].as_ref(),
            query_length,
        );
        let query_hits = &by_query[q_idx];
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180,204-208: a skipped search has a
        // null alignment set, so outfmt 7 has no count line (`report::write_tabular`).
        table.query_header(
            &mut buf,
            tab.query_titles[q_idx],
            tab.database.as_bytes(),
            (!tab.unsearched[q_idx]).then_some(query_hits.len()),
        )?;
        writer.write_all(&buf)?;
        buf.clear();
        // The query's formatter (`TabularRows`).
        let mut rows = TabularRows::default();
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
            let subject = values.nucleotide_subject(names, qcovs, -1);
            // A subject whose lookups fail shows what the last found subject set
            // (`TabularRows::last_found`), made again from its index.
            let kept_values = match rows.last_found() {
                Some(kept) if table.fails(&subject) => Some(subjects.values(kept, spec, PROGRAM)?),
                _ => None,
            };
            let kept = kept_values
                .as_ref()
                .map(|values| values.nucleotide_subject(names, -1, -1));
            let subject_record_length = i32::try_from(record_length(s_idx))?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1222-1225
            // The subject's `qcovhsp` rule (an ungapped search; its rows only divide): the HSP's
            // residues in its query frame against the length that `BlastHitList2SeqAlign_OMF`
            // gives, the searched query's in nucleotides (`query.GetSeqLength`), about a third of
            // the HSP's nucleotide coverage, as NCBI writes it; a subject's single HSP has the
            // Seq-align's nucleotide span over `qlen`.
            let qcovhsp = QueryCoveragePerHsp::of_subject(
                true,
                group.len(),
                i32::try_from(searched_length)?,
                i32::try_from(query_length)?,
            );
            for &hsp_index in group {
                let hsp = &hits[hsp_index];
                let hit = &hsp.hit;
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:762-768,771-772
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
                // A row whose subject's load throws prints none of its letters.
                let rows_read =
                    reads_rows && !subject.throws() && (counts_block || engine_ident > 0);
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-838
                // A row that reads the subject (`TabularSpec::subject_fetch`, by the engine's
                // identities for `nident` alone) loads it first where the query's batch did not
                // prefetch it: the texts of its taxonomy lookups follow the rows before it
                // (`QueryWarnings::subject_read`), before those of the row's own lookups. A row
                // that reads a subject whose load throws outside the subject block ends the run
                // (`TabularWriter::ends_run`).
                if let Some(warnings) = warnings.as_deref_mut() {
                    if warnings.reads_with_texts(q_idx) && fetch.reads(engine_ident) {
                        warnings.subject_read(q_idx, hit.s_idx, &mut *writer)?;
                    }
                }
                if let Some(end) = table.ends_run(
                    &subject,
                    crate::blastinput::app::Program::Tblastx,
                    engine_ident,
                ) {
                    return Err(crate::report::tabular::run_end(writer, &buf, end, || {
                        if let Some(warnings) = warnings.as_deref_mut() {
                            warnings.end_run();
                        }
                    }));
                }
                // The texts of the row's lookups of the taxonomy database (`x_SetTaxInfoAll`,
                // `x_SetTaxInfo`) follow the rows before it.
                if let (Some(names), Some(warnings), true) =
                    (names, warnings.as_deref_mut(), table.sets_subject(&subject))
                {
                    let text = names.row_messages(spec, subject.taxonomy);
                    if !text.is_empty() {
                        writer.flush()?;
                        warnings.sink.write_all(&text)?;
                        warnings.sink.flush()?;
                    }
                }
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
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1233-1244
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
                    end_gaps: true,
                    stored_gaps: [&[], &[]],
                    qcovhsp: qcovhsp.of_hsp(
                        i32::try_from(hit.sort_query_end - hit.sort_query_offset)?,
                        i32::try_from(q_high - q_low + 1)?,
                    ),
                };
                // One printed row is one HSP; the probe marks it without changing it.
                if let Some(probe) = probe.as_deref_mut() {
                    writer.flush()?;
                    probe.begin(hsp_index);
                }
                table.row(
                    &mut buf,
                    &mut rows,
                    &tabular_query,
                    &subject,
                    kept.as_ref(),
                    &row,
                );
                writer.write_all(&buf)?;
                buf.clear();
                if let Some(probe) = probe.as_deref_mut() {
                    writer.flush()?;
                    probe.end(hsp_index);
                }
            }
        }
        // `~CBlastTabularInfo` of the query's report flushes the stream.
        crate::report::ostream::ReportWrite::ncbi_flush(writer)?;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2233-2236
    // A search that stops at a query batch writes no epilog (`report::write_tabular`).
    if spec.kind == TabularKind::Commented && run.epilog {
        warn(&mut *writer, warnings)?;
        table.epilog(writer, run.queries_before + run.queries.len())?;
        crate::report::ostream::ReportWrite::ncbi_flush(writer)?;
    }
    Ok(())
}

#[cfg(test)]
#[path = "tabular_tests.rs"]
mod tests;
